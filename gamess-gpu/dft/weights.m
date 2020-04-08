c***********************************************************************
c
c  Weights.m 
c
c  This file contains all the subroutines that apply so called weighting
c  schemes. See Becke [1] for a discussion of how weighting schemes
c  are used to construct a molecular quadrature grid from atomic grids
c  in polar coordinates.
c
c  [1] "A multicenter numerical integration scheme for polyatomic 
c       molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547
c
c***********************************************************************
c
      subroutine ver_dft_weights(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/weights.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
c
c---- memory counting routines -----------------------------------------
c
      subroutine memreq_calc_weights(igrid,mxp,npts,natoms,ngridc,
     &                               latm,idwght,num_near_atoms,
     &                               memory_fp,memory_int)
      implicit none
c
c     Description:
c
c     This routine estimates the memory requirements for the weights
c     calculation.
c
c     This subroutine is responsible for selecting the appropriate
c     weighting scheme, allocating the required workspaces, invoking
c     the weighting scheme, and releasing the workspaces.
c
c     This routine assumes that under all circumstances there will be
c     a near_atom_list (even when screening is not used).
c
c     For memory allocation the maximum number of points is always used.
c     This guarantees that if the first batch of grid point can be
c     processed then we have enough memory to process all batches. 
c     This way we do not risk a crash half way through the integration.
c
c     Nasty: We cannot afford to store the 2nd derivatives of the 
c            weights, the memory requirements are simply excessive. 
c            As these derivatives are used only in the final Hessian
c            term we can compute that contribution in the weights
c            routine. This is not neat but reality leaves us no other
c            option.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
c
c     Commons:
c
INCLUDE(common/dft_module_comm)
c
c     Input:
c
      integer igrid          ! the grid identifier
      integer mxp            ! the maximum number of points
      integer npts           ! the number of points in the current batch
      integer natoms         ! the number of atoms in the molecule
      integer ngridc         ! the number of grid centres in the system
      integer latm           ! the number of the atom on which the 
                             ! current grid points are centered.
      integer idwght         ! the number of derivative required
      integer num_near_atoms ! the number of atoms near to the grid 
                             ! grid points
c     integer near_atom_list(num_near_atoms)
      REAL rshell            ! the radius of the current batch
      REAL rnear             ! distance of the closest atom to the
                             ! current atom
c
c     Output:
c
      logical dwt_avail_sw   ! true if gradients of the weights have
                             ! been calculated
c
c     Workspace:
c
      REAL memory_fp(*)
      integer memory_int(*)
c
c     Local variables:
c
      integer sk_pt, awt_pt, twt_pt
      integer gZ_pt, gpa_pt, ti_pt, tj_pt
      integer ra2_val_pt, ra2_com_pt, near_pt, indx_pt
      integer aij_g_pt, xij_g_pt, rij_g_pt
      integer aij_pt,   uij_pt,   xij_pt,   rij_pt, pk_pt
      integer duijdi_pt, duijdj_pt, duijdk_pt
      integer dsijdi_pt, dsijdj_pt, dsijdk_pt
      integer d2uijdi2_pt,  d2uijdj2_pt,  d2uijdk2_pt
      integer d2uijdidj_pt, d2uijdidk_pt, d2uijdjdk_pt
      integer d2sijdi2_pt,  d2sijdj2_pt,  d2sijdk2_pt
      integer d2sijdidj_pt, d2sijdidk_pt, d2sijdjdk_pt
      integer g2Z_pt
      integer gwt_pt, g2wt_pt
      integer nnatm
      integer iwrk3,mmax0
      character*9 fnm
      character*12 snm
c
c     Functions:
c
      integer incr_memory
      external incr_memory
c
c     Code:
c
      fnm = 'weights.m'
      snm = 'calc_weights'
      dwt_avail_sw = (idwght.gt.0)
c
c     Allocate workspaces:
c
      nnatm = num_near_atoms
      if (idwght.le.0) then
        sk_pt  = incr_memory(mxp,'d')
        awt_pt = incr_memory(mxp,'d')
        twt_pt = incr_memory(mxp,'d')
      else if (idwght.le.1) then
        sk_pt  = incr_memory(mxp,'d')
        awt_pt = incr_memory(mxp,'d')
        twt_pt = incr_memory(mxp,'d')
        gZ_pt  = incr_memory(3*mxp*nnatm,'d')
        gpa_pt = incr_memory(3*mxp*nnatm,'d')
        ti_pt  = incr_memory(3*mxp*nnatm,'d')
        tj_pt  = incr_memory(3*mxp*nnatm,'d')
      else if (idwght.le.2) then
        aij_g_pt     = incr_memory(nnatm*nnatm,'d')
        xij_g_pt     = incr_memory(3*nnatm*nnatm,'d')
        rij_g_pt     = incr_memory(2*nnatm*nnatm,'d')
        aij_pt       = incr_memory(nnatm*nnatm,'d')
        uij_pt       = incr_memory(nnatm*nnatm,'d')
        xij_pt       = incr_memory(3*nnatm*nnatm,'d')
        rij_pt       = incr_memory(2*nnatm*nnatm,'d')
        ra2_val_pt   = incr_memory(nnatm,'d')
        ra2_com_pt   = incr_memory(3*nnatm,'d')
        near_pt      = incr_memory(nnatm,'i')
        indx_pt      = incr_memory(nnatm,'i')
        mmax0        = max(mxp,nnatm)
        pk_pt        = incr_memory(mmax0,'d')
        duijdi_pt    = incr_memory(3*nnatm*nnatm,'d')
        duijdj_pt    = incr_memory(3*nnatm*nnatm,'d')
        duijdk_pt    = incr_memory(3*nnatm*nnatm,'d')
        dsijdi_pt    = incr_memory(3*nnatm*nnatm,'d')
        dsijdj_pt    = incr_memory(3*nnatm*nnatm,'d')
        dsijdk_pt    = incr_memory(3*nnatm*nnatm,'d')
        d2uijdi2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2uijdj2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2uijdk2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2uijdidj_pt = incr_memory(9*nnatm*nnatm,'d')
        d2uijdidk_pt = incr_memory(9*nnatm*nnatm,'d')
        d2uijdjdk_pt = incr_memory(9*nnatm*nnatm,'d')
        d2sijdi2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2sijdj2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2sijdk2_pt  = incr_memory(9*nnatm*nnatm,'d')
        d2sijdidj_pt = incr_memory(9*nnatm*nnatm,'d')
        d2sijdidk_pt = incr_memory(9*nnatm*nnatm,'d')
        d2sijdjdk_pt = incr_memory(9*nnatm*nnatm,'d')
        gZ_pt        = incr_memory(3*nnatm,'d')
        g2Z_pt       = incr_memory(3*nnatm*3*nnatm,'d')
        gwt_pt       = incr_memory(3*nnatm,'d')
        g2wt_pt      = incr_memory(3*nnatm*3*nnatm,'d')
      else
        call caserr('no 3rd or higher derivatives implemented')
      endif
c
c     Call the appropriate weighting scheme:
c
      if(weight_scheme(igrid) .eq. WT_BECKE)then
c     ...
      else if(weight_scheme(igrid) .eq. WT_SSFSCR)then
c
        iwrk3 = incr_memory(mxp,'i')
        if (idwght.le.0) then
c         call ssfwt_scr(ra2_val,latm,arad,wt,
        else if (idwght.le.1) then
c         call ssfgwt_scr(ra2_val,ra2_comp,latm,
        else if (idwght.le.2) then
c         call ssfg2wt_scr(ra2_val,ra2_comp,latm,
        endif
        call decr_memory(iwrk3,'i')
c
      else if(weight_scheme(igrid) .eq. WT_MHL4SSFSCR)then
c
        iwrk3 = incr_memory(mxp,'i')
        if (idwght.le.0) then
c         call mhl4ssfwt_scr(ra2_val,latm,arad,wt,memory_fp(sk_pt),
        else if (idwght.le.1) then
c         call mhl4ssfgwt_scr(ra2_val,ra2_comp,
        else if (idwght.le.2) then
c         call mhl4ssfg2wt_scr(ra2_val,ra2_comp,latm,num_near_atoms,
        endif
        call decr_memory(iwrk3,'i')
c
      else if(weight_scheme(igrid) .eq. WT_MHL8SSFSCR)then
c
        iwrk3 = incr_memory(mxp,'i')
        if (idwght.le.0) then
c         call mhl8ssfwt_scr(ra2_val,latm,arad,wt,memory_fp(sk_pt),
        else if (idwght.le.1) then
c         call mhl8ssfgwt_scr(ra2_val,ra2_comp,latm,arad,wt,dwt,
        else if (idwght.le.2) then
c         call mhl8ssfg2wt_scr(ra2_val,ra2_comp,latm,num_near_atoms,
        endif
        call decr_memory(iwrk3,'i')
c
      else
c        call caserr('calc_weights: unknown weighting scheme!?')
      endif
c
c     Release workspace memory:
c
      if (idwght.le.0) then
        call decr_memory(twt_pt,'d')
        call decr_memory(awt_pt,'d')
        call decr_memory(sk_pt,'d')
      else if (idwght.le.1) then
        call decr_memory(tj_pt,'d')
        call decr_memory(ti_pt,'d')
        call decr_memory(gpa_pt,'d')
        call decr_memory(gZ_pt,'d')
        call decr_memory(twt_pt,'d')
        call decr_memory(awt_pt,'d')
        call decr_memory(sk_pt,'d')
      else if (idwght.le.2) then
        call decr_memory(g2wt_pt,'d')
        call decr_memory(gwt_pt,'d')
        call decr_memory(g2Z_pt,'d')
        call decr_memory(gZ_pt,'d')
        call decr_memory(d2sijdjdk_pt,'d')
        call decr_memory(d2sijdidk_pt,'d')
        call decr_memory(d2sijdidj_pt,'d')
        call decr_memory(d2sijdk2_pt,'d')
        call decr_memory(d2sijdj2_pt,'d')
        call decr_memory(d2sijdi2_pt,'d')
        call decr_memory(d2uijdjdk_pt,'d')
        call decr_memory(d2uijdidk_pt,'d')
        call decr_memory(d2uijdidj_pt,'d')
        call decr_memory(d2uijdk2_pt,'d')
        call decr_memory(d2uijdj2_pt,'d')
        call decr_memory(d2uijdi2_pt,'d')
        call decr_memory(dsijdk_pt,'d')
        call decr_memory(dsijdj_pt,'d')
        call decr_memory(dsijdi_pt,'d')
        call decr_memory(duijdk_pt,'d')
        call decr_memory(duijdj_pt,'d')
        call decr_memory(duijdi_pt,'d')
        call decr_memory(pk_pt,'d')
        call decr_memory(indx_pt,'i')
        call decr_memory(near_pt,'i')
        call decr_memory(ra2_com_pt,'d')
        call decr_memory(ra2_val_pt,'d')
        call decr_memory(rij_pt,'d')
        call decr_memory(xij_pt,'d')
        call decr_memory(uij_pt,'d')
        call decr_memory(aij_pt,'d')
        call decr_memory(rij_g_pt,'d')
        call decr_memory(xij_g_pt,'d')
        call decr_memory(aij_g_pt,'d')
      endif
c
      end
c
c---- the routines that do the real work -------------------------------
c
      subroutine set_weight_derivative_level(idwght,gradwght_sw,grad_sw,
     &                                       dksm_exp_sw,hess_sw)
      implicit none
c
c     Description:
c
c     Determines how many derivatives of the weights are needed to
c     evaluate the terms currently requested. 
c
c     Input:
c
      logical gradwght_sw  ! are gradients of the weight requested.
      logical grad_sw      ! is the gradient of the energy requested.
      logical dksm_exp_sw  ! are gradients of the Kohn-Sham matrix
                           ! requested.
      logical hess_sw      ! are the second derivatives of the energy
                           ! requested.
c
c     Output:
c
      integer idwght       ! the number of derivatives of the weights
                           ! needed. I.e.
                           ! 0 no derivatives needed
                           ! 1 gradients of the weights needed
                           ! 2 seconds derivatives of the weights needed
c
c     Code:
c
      if (gradwght_sw) then
        if (hess_sw) then
          idwght = 2
        else if (grad_sw.or.dksm_exp_sw) then
          idwght = 1
        else
          idwght = 0
        endif
      else
        idwght = 0
      endif
c
      end
c
c-----------------------------------------------------------------------
c
c     CALC_WEIGHTS stuff
c
c-----------------------------------------------------------------------
c
      subroutine calc_weights(igrid,mxp,npts,natoms,ngridc,latm,idwght,
     &                        num_near_atoms,near_atom_list,rshell,
     &                        rnear,ra2_val,ra2_comp,
     &                        arad,wt,dwt,dwt_avail_sw,
     &                        xc_ept,hess,
     &                        memory_fp,memory_int)
      implicit none
c
c     Description:
c
c     This subroutine is responsible for selecting the appropriate
c     weighting scheme, allocating the required workspaces, invoking
c     the weighting scheme, and releasing the workspaces.
c
c     This routine assumes that under all circumstances there will be
c     a near_atom_list (even when screening is not used).
c
c     For memory allocation the maximum number of points is always used.
c     This guarantees that if the first batch of grid point can be
c     processed then we have enough memory to process all batches. 
c     This way we do not risk a crash half way through the integration.
c
c     Nasty: We cannot afford to store the 2nd derivatives of the 
c            weights, the memory requirements are simply excessive. 
c            As these derivatives are used only in the final Hessian
c            term we can compute that contribution in the weights
c            routine. This is not neat but reality leaves us no other
c            option.
c
c     Parameters:
c
INCLUDE(common/dft_parameters)
c
c     Commons:
c
INCLUDE(common/dft_module_comm)
c
c     Input:
c
      integer igrid          ! the grid identifier
      integer mxp            ! the maximum number of points
      integer npts           ! the number of points in the current batch
      integer natoms         ! the number of atoms in the molecule
      integer ngridc         ! the number of grid centres in the system
      integer latm           ! the number of the atom on which the 
                             ! current grid points are centred.
      integer idwght         ! the number of derivative required
      integer num_near_atoms ! the number of atoms near to the grid 
                             ! grid points
      integer near_atom_list(num_near_atoms)
      REAL rshell            ! the radius of the current batch
      REAL rnear             ! distance of the closest atom to the
                             ! current atom
      REAL ra2_val(mxp,natoms,2)    ! distance(point,atom)
      REAL ra2_comp(mxp,natoms,3) ! x_point-x_atom
      REAL arad(max_atom)    ! radii of the atoms
      REAL xc_ept(mxp)       ! the value of the functional
c
c     Input/Output:
c
      REAL wt(mxp)           ! on entry the atomic weights
                             ! on return the molecular weights
      REAL hess(3*natoms,3*natoms) ! the hessian matrix
c
c     Output:
c
      REAL dwt(3,mxp,num_near_atoms)
c     REAL d2wt(mxp,3*num_near_atoms,3*num_near_atoms)
      logical dwt_avail_sw   ! true if gradients of the weights have
                             ! been calculated
c
c     Workspace:
c
      REAL memory_fp(*)
      integer memory_int(*)
c
c     Local variables:
c
      integer sk_pt, awt_pt, twt_pt
      integer gZ_pt, gpa_pt, ti_pt, tj_pt
      integer ra2_val_pt, ra2_com_pt, near_pt, indx_pt
      integer aij_g_pt, xij_g_pt, rij_g_pt
      integer aij_pt,   uij_pt,   xij_pt,   rij_pt, pk_pt
      integer duijdi_pt, duijdj_pt, duijdk_pt
      integer dsijdi_pt, dsijdj_pt, dsijdk_pt
      integer d2uijdi2_pt,  d2uijdj2_pt,  d2uijdk2_pt
      integer d2uijdidj_pt, d2uijdidk_pt, d2uijdjdk_pt
      integer d2sijdi2_pt,  d2sijdj2_pt,  d2sijdk2_pt
      integer d2sijdidj_pt, d2sijdidk_pt, d2sijdjdk_pt
      integer g2Z_pt
      integer gwt_pt, g2wt_pt
      integer nnatm, nprt
      integer iwrk3, mmax0
      character*9 fnm
      character*12 snm
c
c     Functions:
c
      integer allocate_memory2
      external allocate_memory2
c
c     Code:
c
      fnm = 'weights.m'
      snm = 'calc_weights'
      dwt_avail_sw = (idwght.gt.0)
c
c     Allocate workspaces:
c
      nprt = 3*natoms
      nnatm = num_near_atoms
      if (idwght.le.0) then
        sk_pt  = allocate_memory2(mxp,'d',fnm,snm,'sk')
        awt_pt = allocate_memory2(mxp,'d',fnm,snm,'awt')
        twt_pt = allocate_memory2(mxp,'d',fnm,snm,'twt')
      else if (idwght.le.1) then
        sk_pt  = allocate_memory2(mxp,'d',fnm,snm,'sk')
        awt_pt = allocate_memory2(mxp,'d',fnm,snm,'awt')
        twt_pt = allocate_memory2(mxp,'d',fnm,snm,'twt')
        gZ_pt  = allocate_memory2(3*mxp*nnatm,'d',fnm,snm,'gZ')
        gpa_pt = allocate_memory2(3*mxp*nnatm,'d',fnm,snm,'gpa')
        ti_pt  = allocate_memory2(3*mxp*nnatm,'d',fnm,snm,'ti')
        tj_pt  = allocate_memory2(3*mxp*nnatm,'d',fnm,snm,'tj')
      else if (idwght.le.2) then
        aij_g_pt   = allocate_memory2(nnatm*nnatm,'d',fnm,snm,'aij_g')
        xij_g_pt   = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,'xij_g')
        rij_g_pt   = allocate_memory2(2*nnatm*nnatm,'d',fnm,snm,'rij_g')
        aij_pt     = allocate_memory2(nnatm*nnatm,'d',fnm,snm,'aij')
        uij_pt     = allocate_memory2(nnatm*nnatm,'d',fnm,snm,'uij')
        xij_pt     = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,'xij')
        rij_pt     = allocate_memory2(2*nnatm*nnatm,'d',fnm,snm,'rij')
        ra2_val_pt = allocate_memory2(nnatm,'d',fnm,snm,'ra2_val_pt')
        ra2_com_pt = allocate_memory2(3*nnatm,'d',fnm,snm,'ra2_com_pt')
        near_pt    = allocate_memory2(nnatm,'i',fnm,snm,'near_pt')
        indx_pt    = allocate_memory2(nnatm,'i',fnm,snm,'indx_pt')
        mmax0      =   max(mxp,nnatm)
        pk_pt        = allocate_memory2(mmax0,'d',fnm,snm,
     &                                  'pk')
        duijdi_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'duijdi')
        duijdj_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'duijdj')
        duijdk_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'duijdk')
        dsijdi_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'dsijdi')
        dsijdj_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'dsijdj')
        dsijdk_pt    = allocate_memory2(3*nnatm*nnatm,'d',fnm,snm,
     &                                  'dsijdk')
        d2uijdi2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdi2')
        d2uijdj2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdj2')
        d2uijdk2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdk2')
        d2uijdidj_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdidj')
        d2uijdidk_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdidk')
        d2uijdjdk_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2uijdjdk')
        d2sijdi2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdi2')
        d2sijdj2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdj2')
        d2sijdk2_pt  = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdk2')
        d2sijdidj_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdidj')
        d2sijdidk_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdidk')
        d2sijdjdk_pt = allocate_memory2(9*nnatm*nnatm,'d',fnm,snm,
     &                                  'd2sijdjdk')
        gZ_pt        = allocate_memory2(3*nnatm,'d',fnm,snm,'gZ')
        g2Z_pt       = allocate_memory2(3*nnatm*3*nnatm,'d',fnm,snm,
     &                                  'g2Z')
        gwt_pt       = allocate_memory2(3*nnatm,'d',fnm,snm,'gwt')
        g2wt_pt      = allocate_memory2(3*nnatm*3*nnatm,'d',fnm,snm,
     &                                  'g2wt')
      else
        call caserr('no 3rd or higher derivatives implemented')
      endif
c
c     Call the appropriate weighting scheme:
c
      if(weight_scheme(igrid) .eq. WT_BECKE)then
c     
        if (idwght.le.0) then
          call beckewt(ra2_val(1,1,2),latm,wt,memory_fp(sk_pt),
     &                 memory_fp(awt_pt),memory_fp(twt_pt),npts,mxp,
     &                 ngridc,igrid)
        else if (idwght.le.1) then
          call beckegwt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                  memory_fp(sk_pt),
     &                  memory_fp(awt_pt),memory_fp(twt_pt),
     &                  memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                  memory_fp(ti_pt),memory_fp(tj_pt),
     &                  npts,mxp,ngridc,igrid)
        else if (idwght.le.2) then
          call beckeg2wt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                   memory_fp(aij_pt),memory_fp(uij_pt),
     &                   memory_fp(xij_pt),memory_fp(rij_pt),
     &                   memory_fp(pk_pt),
     &                   memory_fp(duijdi_pt),
     &                   memory_fp(duijdj_pt),
     &                   memory_fp(duijdk_pt),
     &                   memory_fp(dsijdi_pt),
     &                   memory_fp(dsijdj_pt),
     &                   memory_fp(dsijdk_pt),
     &                   memory_fp(gZ_pt),
     &                   memory_fp(d2uijdi2_pt),
     &                   memory_fp(d2uijdj2_pt),
     &                   memory_fp(d2uijdk2_pt),
     &                   memory_fp(d2uijdidj_pt),
     &                   memory_fp(d2uijdidk_pt),
     &                   memory_fp(d2uijdjdk_pt),
     &                   memory_fp(d2sijdi2_pt),
     &                   memory_fp(d2sijdj2_pt),
     &                   memory_fp(d2sijdk2_pt),
     &                   memory_fp(d2sijdidj_pt),
     &                   memory_fp(d2sijdidk_pt),
     &                   memory_fp(d2sijdjdk_pt),
     &                   memory_fp(g2Z_pt),
     &                   memory_fp(g2wt_pt),
     &                   xc_ept,hess,nprt,
     &                   npts,mxp,ngridc,igrid)
        endif
c
      else if(weight_scheme(igrid) .eq. WT_BECKESCR)then
c
        if (idwght.le.0) then
          call beckewt_scr(ra2_val(1,1,2),latm,num_near_atoms,
     &                     near_atom_list,
     &                     arad,wt,memory_fp(sk_pt),memory_fp(awt_pt),
     &                     memory_fp(twt_pt),npts,mxp,ngridc,igrid)
        else if (idwght.le.1) then
          call beckegwt_scr(ra2_val(1,1,2),ra2_comp,latm,num_near_atoms,
     &                      near_atom_list,arad,wt,dwt,dwt_avail_sw,
     &                      memory_fp(sk_pt),memory_fp(awt_pt),
     &                      memory_fp(twt_pt),memory_fp(gpa_pt),
     &                      memory_fp(gZ_pt),memory_fp(ti_pt),
     &                      memory_fp(tj_pt),npts,mxp,ngridc,igrid)
        else if (idwght.le.2) then
          call beckeg2wt_scr(ra2_val(1,1,2),ra2_comp,
     &                       memory_fp(ra2_val_pt),
     &                       memory_fp(ra2_com_pt),latm,
     &                       num_near_atoms,near_atom_list,
     &                       memory_int(near_pt),
     &                       memory_int(indx_pt),
     &                       arad,
     &                       wt,dwt,dwt_avail_sw,
     &                       memory_fp(aij_g_pt),
     &                       memory_fp(xij_g_pt),memory_fp(rij_g_pt),
     &                       memory_fp(aij_pt),memory_fp(uij_pt),
     &                       memory_fp(xij_pt),memory_fp(rij_pt),
     &                       memory_fp(pk_pt),
     &                       memory_fp(duijdi_pt),
     &                       memory_fp(duijdj_pt),
     &                       memory_fp(duijdk_pt),
     &                       memory_fp(dsijdi_pt),
     &                       memory_fp(dsijdj_pt),
     &                       memory_fp(dsijdk_pt),
     &                       memory_fp(gZ_pt),
     &                       memory_fp(gwt_pt),
     &                       memory_fp(d2uijdi2_pt),
     &                       memory_fp(d2uijdj2_pt),
     &                       memory_fp(d2uijdk2_pt),
     &                       memory_fp(d2uijdidj_pt),
     &                       memory_fp(d2uijdidk_pt),
     &                       memory_fp(d2uijdjdk_pt),
     &                       memory_fp(d2sijdi2_pt),
     &                       memory_fp(d2sijdj2_pt),
     &                       memory_fp(d2sijdk2_pt),
     &                       memory_fp(d2sijdidj_pt),
     &                       memory_fp(d2sijdidk_pt),
     &                       memory_fp(d2sijdjdk_pt),
     &                       memory_fp(g2Z_pt),
     &                       memory_fp(g2wt_pt),
     &                       xc_ept,hess,nprt,
     &                       npts,mxp,ngridc,igrid)
        endif
c
      else if(weight_scheme(igrid) .eq. WT_MHL)then
c
        if (idwght.le.0) then
          call mhlwt(ra2_val(1,1,2),latm,wt,memory_fp(sk_pt),
     &               memory_fp(awt_pt),memory_fp(twt_pt),
     &               npts,mxp,ngridc,igrid)
        else if (idwght.le.1) then
          call mhlgwt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                memory_fp(sk_pt),
     &                memory_fp(awt_pt),memory_fp(twt_pt),
     &                memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                memory_fp(ti_pt),memory_fp(tj_pt),
     &                npts,mxp,ngridc,igrid)
        else if (idwght.le.2) then
          call mhlg2wt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                 memory_fp(aij_pt),memory_fp(uij_pt),
     &                 memory_fp(xij_pt),memory_fp(rij_pt),
     &                 memory_fp(pk_pt),
     &                 memory_fp(duijdi_pt),
     &                 memory_fp(duijdj_pt),
     &                 memory_fp(duijdk_pt),
     &                 memory_fp(dsijdi_pt),
     &                 memory_fp(dsijdj_pt),
     &                 memory_fp(dsijdk_pt),
     &                 memory_fp(gZ_pt),
     &                 memory_fp(d2uijdi2_pt),
     &                 memory_fp(d2uijdj2_pt),
     &                 memory_fp(d2uijdk2_pt),
     &                 memory_fp(d2uijdidj_pt),
     &                 memory_fp(d2uijdidk_pt),
     &                 memory_fp(d2uijdjdk_pt),
     &                 memory_fp(d2sijdi2_pt),
     &                 memory_fp(d2sijdj2_pt),
     &                 memory_fp(d2sijdk2_pt),
     &                 memory_fp(d2sijdidj_pt),
     &                 memory_fp(d2sijdidk_pt),
     &                 memory_fp(d2sijdjdk_pt),
     &                 memory_fp(g2Z_pt),
     &                 memory_fp(g2wt_pt),
     &                 xc_ept,hess,nprt,
     &                 npts,mxp,ngridc,igrid)
        endif
c
      else if(weight_scheme(igrid) .eq. WT_MHLSCR)then
c
        if (idwght.le.0) then
          call mhlwt_scr(ra2_val(1,1,2),latm,num_near_atoms,
     &                   near_atom_list,
     &                   arad,wt,memory_fp(sk_pt),memory_fp(awt_pt),
     &                   memory_fp(twt_pt),npts,mxp,ngridc,igrid)
        else if (idwght.le.1) then
          call mhlgwt_scr(ra2_val(1,1,2),ra2_comp,latm,num_near_atoms,
     &                    near_atom_list,arad,wt,dwt,dwt_avail_sw,
     &                    memory_fp(sk_pt),memory_fp(awt_pt),
     &                    memory_fp(twt_pt),
     &                    memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                    memory_fp(ti_pt),memory_fp(tj_pt),
     &                    npts,mxp,ngridc,igrid)
        else if (idwght.le.2) then
          call mhlg2wt_scr(ra2_val(1,1,2),ra2_comp,
     &                     memory_fp(ra2_val_pt),
     &                     memory_fp(ra2_com_pt),
     &                     latm,num_near_atoms,near_atom_list,
     &                     memory_int(near_pt),
     &                     memory_int(indx_pt),
     &                     arad,wt,dwt,dwt_avail_sw,
     &                     memory_fp(aij_g_pt),
     &                     memory_fp(xij_g_pt),memory_fp(rij_g_pt),
     &                     memory_fp(aij_pt),memory_fp(uij_pt),
     &                     memory_fp(xij_pt),memory_fp(rij_pt),
     &                     memory_fp(pk_pt),
     &                     memory_fp(duijdi_pt),
     &                     memory_fp(duijdj_pt),
     &                     memory_fp(duijdk_pt),
     &                     memory_fp(dsijdi_pt),
     &                     memory_fp(dsijdj_pt),
     &                     memory_fp(dsijdk_pt),
     &                     memory_fp(gZ_pt),
     &                     memory_fp(gwt_pt),
     &                     memory_fp(d2uijdi2_pt),
     &                     memory_fp(d2uijdj2_pt),
     &                     memory_fp(d2uijdk2_pt),
     &                     memory_fp(d2uijdidj_pt),
     &                     memory_fp(d2uijdidk_pt),
     &                     memory_fp(d2uijdjdk_pt),
     &                     memory_fp(d2sijdi2_pt),
     &                     memory_fp(d2sijdj2_pt),
     &                     memory_fp(d2sijdk2_pt),
     &                     memory_fp(d2sijdidj_pt),
     &                     memory_fp(d2sijdidk_pt),
     &                     memory_fp(d2sijdjdk_pt),
     &                     memory_fp(g2Z_pt),
     &                     memory_fp(g2wt_pt),
     &                     xc_ept,hess,nprt,
     &                     npts,mxp,ngridc,igrid)
        endif
c
      else if(weight_scheme(igrid) .eq. WT_SSF)then
c
        if (idwght.le.0) then
          call ssfwt(ra2_val(1,1,2),latm,wt,memory_fp(sk_pt),
     &               memory_fp(awt_pt),memory_fp(twt_pt),
     &               npts,mxp,ngridc,rshell,rnear)
        else if (idwght.le.1) then
          call ssfgwt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                memory_fp(sk_pt),memory_fp(awt_pt),
     &                memory_fp(twt_pt),
     &                memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                memory_fp(ti_pt),memory_fp(tj_pt),
     &                npts,mxp,ngridc,rshell,rnear)
        else if (idwght.le.2) then
          call ssfg2wt(ra2_val(1,1,2),ra2_comp,latm,wt,dwt,
     &                 memory_fp(uij_pt),memory_fp(xij_pt),
     &                 memory_fp(rij_pt),memory_fp(pk_pt),
     &                 memory_fp(duijdi_pt),
     &                 memory_fp(duijdj_pt),
     &                 memory_fp(duijdk_pt),
     &                 memory_fp(dsijdi_pt),
     &                 memory_fp(dsijdj_pt),
     &                 memory_fp(dsijdk_pt),
     &                 memory_fp(gZ_pt),
     &                 memory_fp(d2uijdi2_pt),
     &                 memory_fp(d2uijdj2_pt),
     &                 memory_fp(d2uijdk2_pt),
     &                 memory_fp(d2uijdidj_pt),
     &                 memory_fp(d2uijdidk_pt),
     &                 memory_fp(d2uijdjdk_pt),
     &                 memory_fp(d2sijdi2_pt),
     &                 memory_fp(d2sijdj2_pt),
     &                 memory_fp(d2sijdk2_pt),
     &                 memory_fp(d2sijdidj_pt),
     &                 memory_fp(d2sijdidk_pt),
     &                 memory_fp(d2sijdjdk_pt),
     &                 memory_fp(g2Z_pt),
     &                 memory_fp(g2wt_pt),
     &                 xc_ept,hess,nprt,
     &                 npts,mxp,ngridc,rshell,rnear)
        endif
c
      else if(weight_scheme(igrid) .eq. WT_SSFSCR)then
c
        iwrk3 = allocate_memory2(npts,'i',fnm,snm,'iwrk3')
        if (idwght.le.0) then
          call ssfwt_scr(ra2_val(1,1,2),latm,arad,wt,
     &                   memory_fp(sk_pt),memory_fp(awt_pt),
     &                   memory_fp(twt_pt),npts,mxp,ngridc,
     &                   near_atom_list, num_near_atoms,
     &                   rshell, rnear, memory_int(iwrk3))
        else if (idwght.le.1) then
          call ssfgwt_scr(ra2_val(1,1,2),ra2_comp,latm,
     &                    arad,wt,dwt,dwt_avail_sw,
     &                    memory_fp(sk_pt),memory_fp(awt_pt),
     &                    memory_fp(twt_pt),
     &                    memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                    memory_fp(ti_pt),memory_fp(tj_pt),
     &                    npts,mxp,ngridc,near_atom_list,num_near_atoms,
     &                    rshell, rnear, memory_int(iwrk3))
        else if (idwght.le.2) then
          call ssfg2wt_scr(ra2_val(1,1,2),ra2_comp,
     &                     memory_fp(ra2_val_pt),
     &                     memory_fp(ra2_com_pt),
     &                     latm,num_near_atoms,near_atom_list,
     &                     memory_int(near_pt),
     &                     memory_int(indx_pt),
     &                     arad,wt,dwt,dwt_avail_sw,
     &                     memory_fp(xij_g_pt),memory_fp(rij_g_pt),
     &                     memory_fp(uij_pt),
     &                     memory_fp(xij_pt),memory_fp(rij_pt),
     &                     memory_fp(pk_pt),
     &                     memory_fp(duijdi_pt),
     &                     memory_fp(duijdj_pt),
     &                     memory_fp(duijdk_pt),
     &                     memory_fp(dsijdi_pt),
     &                     memory_fp(dsijdj_pt),
     &                     memory_fp(dsijdk_pt),
     &                     memory_fp(gZ_pt),
     &                     memory_fp(gwt_pt),
     &                     memory_fp(d2uijdi2_pt),
     &                     memory_fp(d2uijdj2_pt),
     &                     memory_fp(d2uijdk2_pt),
     &                     memory_fp(d2uijdidj_pt),
     &                     memory_fp(d2uijdidk_pt),
     &                     memory_fp(d2uijdjdk_pt),
     &                     memory_fp(d2sijdi2_pt),
     &                     memory_fp(d2sijdj2_pt),
     &                     memory_fp(d2sijdk2_pt),
     &                     memory_fp(d2sijdidj_pt),
     &                     memory_fp(d2sijdidk_pt),
     &                     memory_fp(d2sijdjdk_pt),
     &                     memory_fp(g2Z_pt),
     &                     memory_fp(g2wt_pt),
     &                     xc_ept,hess,nprt,
     &                     npts,mxp,ngridc,rshell,rnear)
        endif
        call free_memory2(iwrk3,'i',fnm,snm,'iwrk3')
c
      else if(weight_scheme(igrid) .eq. WT_MHL4SSFSCR)then
c
        iwrk3 = allocate_memory2(npts,'i',fnm,snm,'iwrk3')
        if (idwght.le.0) then
          call mhl4ssfwt_scr(ra2_val(1,1,2),latm,arad,wt,
     &                       memory_fp(sk_pt),
     &                       memory_fp(awt_pt),memory_fp(twt_pt),
     &                       npts,mxp,ngridc,near_atom_list, 
     &                       num_near_atoms,
     &                       rshell, rnear, memory_int(iwrk3),igrid)
        else if (idwght.le.1) then
          call mhl4ssfgwt_scr(ra2_val(1,1,2),ra2_comp,
     &                        latm,arad,wt,dwt,dwt_avail_sw,
     &                        memory_fp(sk_pt),memory_fp(awt_pt),
     &                        memory_fp(twt_pt),
     &                        memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                        memory_fp(ti_pt),memory_fp(tj_pt),
     &                        npts,mxp,ngridc,near_atom_list,
     &                        num_near_atoms,
     &                        rshell, rnear, memory_int(iwrk3),igrid)
        else if (idwght.le.2) then
          call mhl4ssfg2wt_scr(ra2_val(1,1,2),ra2_comp,
     &                         memory_fp(ra2_val_pt),
     &                         memory_fp(ra2_com_pt),
     &                         latm,num_near_atoms,near_atom_list,
     &                         memory_int(near_pt),
     &                         memory_int(indx_pt),
     &                         arad,wt,dwt,
     &                         dwt_avail_sw,
     &                         memory_fp(aij_g_pt),
     &                         memory_fp(xij_g_pt),memory_fp(rij_g_pt),
     &                         memory_fp(aij_pt),memory_fp(uij_pt),
     &                         memory_fp(xij_pt),memory_fp(rij_pt),
     &                         memory_fp(pk_pt),
     &                         memory_fp(duijdi_pt),
     &                         memory_fp(duijdj_pt),
     &                         memory_fp(duijdk_pt),
     &                         memory_fp(dsijdi_pt),
     &                         memory_fp(dsijdj_pt),
     &                         memory_fp(dsijdk_pt),
     &                         memory_fp(gZ_pt),
     &                         memory_fp(gwt_pt),
     &                         memory_fp(d2uijdi2_pt),
     &                         memory_fp(d2uijdj2_pt),
     &                         memory_fp(d2uijdk2_pt),
     &                         memory_fp(d2uijdidj_pt),
     &                         memory_fp(d2uijdidk_pt),
     &                         memory_fp(d2uijdjdk_pt),
     &                         memory_fp(d2sijdi2_pt),
     &                         memory_fp(d2sijdj2_pt),
     &                         memory_fp(d2sijdk2_pt),
     &                         memory_fp(d2sijdidj_pt),
     &                         memory_fp(d2sijdidk_pt),
     &                         memory_fp(d2sijdjdk_pt),
     &                         memory_fp(g2Z_pt),
     &                         memory_fp(g2wt_pt),
     &                         xc_ept,hess,nprt,
     &                         npts,mxp,ngridc,igrid)
        endif
        call free_memory2(iwrk3,'i',fnm,snm,'iwrk3')
c
      else if(weight_scheme(igrid) .eq. WT_MHL8SSFSCR)then
c
        iwrk3 = allocate_memory2(npts,'i',fnm,snm,'iwrk3')
        if (idwght.le.0) then
          call mhl8ssfwt_scr(ra2_val(1,1,2),latm,arad,wt,
     &                       memory_fp(sk_pt),
     &                       memory_fp(awt_pt),memory_fp(twt_pt),
     &                       npts,mxp,ngridc,near_atom_list,
     &                       num_near_atoms,
     &                       rshell,rnear,memory_int(iwrk3),igrid)
        else if (idwght.le.1) then
          call mhl8ssfgwt_scr(ra2_val(1,1,2),ra2_comp,latm,arad,wt,dwt,
     &                        dwt_avail_sw,memory_fp(sk_pt),
     &                        memory_fp(awt_pt),memory_fp(twt_pt),
     &                        memory_fp(gpa_pt),memory_fp(gZ_pt),
     &                        memory_fp(ti_pt),memory_fp(tj_pt),
     &                        npts,mxp,ngridc,near_atom_list,
     &                        num_near_atoms,
     &                        rshell, rnear, memory_int(iwrk3),igrid)
        else if (idwght.le.2) then
          call mhl8ssfg2wt_scr(ra2_val(1,1,2),ra2_comp,
     &                         memory_fp(ra2_val_pt),
     &                         memory_fp(ra2_com_pt),
     &                         latm,num_near_atoms,near_atom_list,
     &                         memory_int(near_pt),
     &                         memory_int(indx_pt),
     &                         arad,wt,dwt,dwt_avail_sw,
     &                         memory_fp(aij_g_pt),
     &                         memory_fp(xij_g_pt),memory_fp(rij_g_pt),
     &                         memory_fp(aij_pt),memory_fp(uij_pt),
     &                         memory_fp(xij_pt),memory_fp(rij_pt),
     &                         memory_fp(pk_pt),
     &                         memory_fp(duijdi_pt),
     &                         memory_fp(duijdj_pt),
     &                         memory_fp(duijdk_pt),
     &                         memory_fp(dsijdi_pt),
     &                         memory_fp(dsijdj_pt),
     &                         memory_fp(dsijdk_pt),
     &                         memory_fp(gZ_pt),
     &                         memory_fp(gwt_pt),
     &                         memory_fp(d2uijdi2_pt),
     &                         memory_fp(d2uijdj2_pt),
     &                         memory_fp(d2uijdk2_pt),
     &                         memory_fp(d2uijdidj_pt),
     &                         memory_fp(d2uijdidk_pt),
     &                         memory_fp(d2uijdjdk_pt),
     &                         memory_fp(d2sijdi2_pt),
     &                         memory_fp(d2sijdj2_pt),
     &                         memory_fp(d2sijdk2_pt),
     &                         memory_fp(d2sijdidj_pt),
     &                         memory_fp(d2sijdidk_pt),
     &                         memory_fp(d2sijdjdk_pt),
     &                         memory_fp(g2Z_pt),
     &                         memory_fp(g2wt_pt),
     &                         xc_ept,hess,nprt,
     &                         npts,mxp,ngridc,igrid)
        endif
        call free_memory2(iwrk3,'i',fnm,snm,'iwrk3')
c
      else
         call caserr('calc_weights: unknown weighting scheme!?')
      endif
c
c     Release workspace memory:
c
      if (idwght.le.0) then
        call free_memory2(twt_pt,'d',fnm,snm,'twt')
        call free_memory2(awt_pt,'d',fnm,snm,'awt')
        call free_memory2(sk_pt,'d',fnm,snm,'sk')
      else if (idwght.le.1) then
        call free_memory2(tj_pt,'d',fnm,snm,'tj')
        call free_memory2(ti_pt,'d',fnm,snm,'ti')
        call free_memory2(gpa_pt,'d',fnm,snm,'gpa')
        call free_memory2(gZ_pt,'d',fnm,snm,'gZ')
        call free_memory2(twt_pt,'d',fnm,snm,'twt')
        call free_memory2(awt_pt,'d',fnm,snm,'awt')
        call free_memory2(sk_pt,'d',fnm,snm,'sk')
      else if (idwght.le.2) then
        call free_memory2(g2wt_pt,'d',fnm,snm,'g2wt')
        call free_memory2(gwt_pt,'d',fnm,snm,'gwt')
        call free_memory2(g2Z_pt,'d',fnm,snm,'g2Z')
        call free_memory2(gZ_pt,'d',fnm,snm,'gZ')
        call free_memory2(d2sijdjdk_pt,'d',fnm,snm,'d2sijdjdk')
        call free_memory2(d2sijdidk_pt,'d',fnm,snm,'d2sijdidk')
        call free_memory2(d2sijdidj_pt,'d',fnm,snm,'d2sijdidj')
        call free_memory2(d2sijdk2_pt,'d',fnm,snm,'d2sijdk2')
        call free_memory2(d2sijdj2_pt,'d',fnm,snm,'d2sijdj2')
        call free_memory2(d2sijdi2_pt,'d',fnm,snm,'d2sijdi2')
        call free_memory2(d2uijdjdk_pt,'d',fnm,snm,'d2uijdjdk')
        call free_memory2(d2uijdidk_pt,'d',fnm,snm,'d2uijdidk')
        call free_memory2(d2uijdidj_pt,'d',fnm,snm,'d2uijdidj')
        call free_memory2(d2uijdk2_pt,'d',fnm,snm,'d2uijdk2')
        call free_memory2(d2uijdj2_pt,'d',fnm,snm,'d2uijdj2')
        call free_memory2(d2uijdi2_pt,'d',fnm,snm,'d2uijdi2')
        call free_memory2(dsijdk_pt,'d',fnm,snm,'dsijdk')
        call free_memory2(dsijdj_pt,'d',fnm,snm,'dsijdj')
        call free_memory2(dsijdi_pt,'d',fnm,snm,'dsijdi')
        call free_memory2(duijdk_pt,'d',fnm,snm,'duijdk')
        call free_memory2(duijdj_pt,'d',fnm,snm,'duijdj')
        call free_memory2(duijdi_pt,'d',fnm,snm,'duijdi')
        call free_memory2(pk_pt,'d',fnm,snm,'pk')
        call free_memory2(indx_pt,'i',fnm,snm,'indx_pt')
        call free_memory2(near_pt,'i',fnm,snm,'near_pt')
        call free_memory2(ra2_com_pt,'d',fnm,snm,'ra2_com_pt')
        call free_memory2(ra2_val_pt,'d',fnm,snm,'ra2_val_pt')
        call free_memory2(rij_pt,'d',fnm,snm,'rij')
        call free_memory2(xij_pt,'d',fnm,snm,'xij')
        call free_memory2(uij_pt,'d',fnm,snm,'uij')
        call free_memory2(aij_pt,'d',fnm,snm,'aij')
        call free_memory2(rij_g_pt,'d',fnm,snm,'rij_g')
        call free_memory2(xij_g_pt,'d',fnm,snm,'xij_g')
        call free_memory2(aij_g_pt,'d',fnm,snm,'aij_g')
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine beckewt(ra2_val,latm,wt,  
     &     sk, awt, totwt, npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:		                                       *
C  *Calculates Becke weight at point                                   *
c  * simplistic version with 3x loop unroll                            *
C  *********************************************************************

      implicit none
C  *********************************************************************
C  *Declarations                                                       *
C  *                                                                   *
C  *Parameters                                                         *
INCLUDE(common/dft_parameters)
C  *In variables		                                       *
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      integer latm
      integer igrd
c
C  *In/out variables                                                   *
      REAL wt(npts)
C  *Work space                                                         *
      REAL awt(npts), totwt(npts), sk(npts)
C  *Local variables                                                    *
      integer i,j,ipt,igrid,jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL uij1,fk1
      integer ipt1, ipt2
      REAL uij2,fk2
C  *End declarations                                                   *
C  *********************************************************************

      call aclear_dp(totwt,npts,0.0d0)
      do 10 i=1,ngridcentres

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10

         call aclear_dp(sk,npts,1.0d0)
         do 20 j=1,ngridcentres
            if(j.eq.i) goto 20

            jgrid = gtype_num(j)
            if (jgrid.eq.0) goto 20

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(i,j)

            do ipt=1,npts
               uij=diji*(ra2_val(ipt,i)-ra2_val(ipt,j))
               fk=uij+aij*(1.0d0-(uij*uij))
               fk=1.5d0*fk-0.5d0*fk*fk*fk
               fk=1.5d0*fk-0.5d0*fk*fk*fk
               fk=1.5d0*fk-0.5d0*fk*fk*fk
               sk(ipt)=sk(ipt)*0.5d0*(1.0d0-fk)
            enddo
 20      continue
         call daxpy(npts,1.0d0,sk,1,totwt,1)
         if(i.eq.latm)call dcopy(npts,sk,1,awt,1)
 10   continue
      
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine beckegwt(ra2_val,ra2_comp,latm,wt,gwt,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:		                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights.                                          *
C  *********************************************************************

      implicit none
C  *********************************************************************
C  *Declarations                                                       *
C  *                                                                   *
C  *Parameters                                                         *
INCLUDE(common/dft_parameters)
C  *In variables		                                       *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      integer latm
      integer igrd
c
C  *In/out variables                                                   *
      REAL wt(npts), gwt(3,mxp,natm)
C  *Work space                                                         *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,natm),gZ(3,npts,natm)
      REAL ti_ij(3,npts,natm), tj_ij(3,npts,natm)
C  *Local variables                                                    *
      integer i,j,ipt,loopct,igrid,jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL uij1,fk1
      integer ipt1, ipt2
      REAL uij2,fk2
      REAL dui_ijx,dui_ijy,dui_ijz
      REAL duj_ijx,duj_ijy,duj_ijz
      REAL vij,dvij,dfk,sk,dsk
C  *Functions                                                          *
c     REAL srad
C  *End declarations                                                   *
C  *********************************************************************
      
      call aclear_dp(Z,npts,0.0d0)
      call aclear_dp(gZ,3*npts*natm,0.0d0)
      call aclear_dp(gpa,3*npts*natm,0.0d0)

      do 10 i=1,natm

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10

         call aclear_dp(pk,npts,1.0d0)
         do 20 j=1,natm
            if(j.eq.i) goto 20

            jgrid = gtype_num(j)
            if (jgrid.eq.0) then
               do ipt=1,npts
                  ti_ij(1,ipt,j)=0.0d0
                  ti_ij(2,ipt,j)=0.0d0
                  ti_ij(3,ipt,j)=0.0d0
                  tj_ij(1,ipt,j)=0.0d0
                  tj_ij(2,ipt,j)=0.0d0
                  tj_ij(3,ipt,j)=0.0d0
               enddo
               goto 20
            endif

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(i,j)

            do ipt=1,npts
               uij=diji*(ra2_val(ipt,i)-ra2_val(ipt,j))
               dui_ijx=-ra2_comp(ipt,i,1)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
               dui_ijy=-ra2_comp(ipt,i,2)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
               dui_ijz=-ra2_comp(ipt,i,3)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
               duj_ijx= ra2_comp(ipt,j,1)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
               duj_ijy= ra2_comp(ipt,j,2)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
               duj_ijz= ra2_comp(ipt,j,3)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
c
               vij=uij+aij*(1.0d0-(uij*uij))
               dvij=1.0d0-2.0d0*aij*uij
c
               dfk=(1.5d0-1.5d0*vij*vij)*dvij
               fk=1.5d0*vij-0.5d0*vij*vij*vij
               dfk=(1.5d0-1.5d0*fk*fk)*dfk
               fk=1.5d0*fk-0.5d0*fk*fk*fk
               dfk=(1.5d0-1.5d0*fk*fk)*dfk
               fk=1.5d0*fk-0.5d0*fk*fk*fk
c
               sk=0.5d0*(1.0d0-fk)
               if (dabs(sk).lt.1.0d-15) then
                  dsk= 0.0d0
               else
                  dsk=-0.5d0*dfk/sk
               endif
c
               pk(ipt)=pk(ipt)*sk
               ti_ij(1,ipt,j)=dsk*dui_ijx
               ti_ij(2,ipt,j)=dsk*dui_ijy
               ti_ij(3,ipt,j)=dsk*dui_ijz
               tj_ij(1,ipt,j)=dsk*duj_ijx
               tj_ij(2,ipt,j)=dsk*duj_ijy
               tj_ij(3,ipt,j)=dsk*duj_ijz
            enddo
 20      continue
c
         call daxpy(npts,1.0d0,pk,1,Z,1)
c
         if(i.eq.latm) then
c
c...        Save Pa
c
            call dcopy(npts,pk,1,pa,1)
c
c...        Construct grad Pa
c
            do j=1,natm
               if (j.ne.i) then
                  do ipt=1,npts
                     gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
                     gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
                     gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
                  enddo
               endif
            enddo
         endif
c
c...     Construct grad Z
c
         do j=1,natm
            if (j.ne.i) then
               do ipt=1,npts
                  gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                  gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                  gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                  gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                  gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                  gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
               enddo
            endif
         enddo
 10   continue
      
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      call aclear_dp(gwt,3*mxp*natm,0.0d0)
      do j=1,natm
         if (j.ne.latm) then
            do ipt=1,npts
               gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(1,ipt,j)
     +                         -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
               gwt(1,ipt,latm) = gwt(1,ipt,latm) - gwt(1,ipt,j)
               gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(2,ipt,j)
     +                         -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
               gwt(2,ipt,latm) = gwt(2,ipt,latm) - gwt(2,ipt,j)
               gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(3,ipt,j)
     +                         -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
               gwt(3,ipt,latm) = gwt(3,ipt,latm) - gwt(3,ipt,j)
            enddo
         endif
      enddo
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine beckeg2wt(ra2_val,ra2_comp,latm,wt,gwt,
     &     aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, 
     &     d2uijdi2,  d2uijdj2,  d2uijdk2, 
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2, 
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, nprt
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      REAL xc_ept(mxp)
      integer latm
      integer igrd
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt(3,mxp,natm)
c
c...  Workspace
c
      REAL aij(natm,natm)
      REAL uij(natm,natm)
      REAL xij(3,natm,natm)
      REAL rij(2,natm,natm) ! rij(i,*) stores dij**(-i)
      REAL pk(natm)
      REAL duijdi(3,natm,natm)
      REAL duijdj(3,natm,natm)
      REAL duijdk(3,natm,natm)
      REAL d2uijdi2(3,3,natm,natm)
      REAL d2uijdj2(3,3,natm,natm)
      REAL d2uijdk2(3,3,natm,natm)
      REAL d2uijdidj(3,3,natm,natm)
      REAL d2uijdidk(3,3,natm,natm)
      REAL d2uijdjdk(3,3,natm,natm)
      REAL gZ(3,natm)
      REAL g2Z(3*natm,3*natm)
      REAL g2wt(3*natm,3*natm)
      REAL dsijdi(3,natm,natm)
      REAL dsijdj(3,natm,natm)
      REAL dsijdk(3,natm,natm)
      REAL d2sijdi2(3,3,natm,natm)
      REAL d2sijdj2(3,3,natm,natm)
      REAL d2sijdk2(3,3,natm,natm)
      REAL d2sijdidj(3,3,natm,natm)
      REAL d2sijdidk(3,3,natm,natm)
      REAL d2sijdjdk(3,3,natm,natm)
c
c...  Local variables
c
      integer i, j, k, l, ipt
      integer ix, jx, kx, lx
      integer ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, f1ij, f2ij, f3ij
      REAL df2ijdf1, df3ijdf2
      REAL dvijdu, df1ijdv, df2ijdv, df3ijdv
      REAL d2vijdu2, d2f1ijdv2, d2f2ijdv2, d2f3ijdv2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
c
c...  Functions
c
c     REAL srad
c
c...  Statement functions
c
c...  iks: small triangle i.e. without diagonal
c...  ikb: big   triangle i.e. with    diagonal
c...  ikm: ikb with min's and max's
c     integer iks, ikb, ikm
c     iks(i,j) = iky(i-1)+j
c     ikb(i,j) = iky(i)+j
c     ikm(i,j) = ikb(max(i,j),min(i,j))
c
c...  Data statements
c
c     data ikx/1,2,3,4,5,6,7,8,9/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      do 5 i=1,natm
         igrid = gtype_num(i)
         if (igrid.eq.0) goto 5
         do 10 j=1,natm
            if (i.eq.j) goto 10
            jgrid = gtype_num(j)
            if (jgrid.eq.0) goto 10
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij(i,j)=tij
c
            xij(1,i,j)=atom_c(i,1)-atom_c(j,1)
            xij(2,i,j)=atom_c(i,2)-atom_c(j,2)
            xij(3,i,j)=atom_c(i,3)-atom_c(j,3)
c
            tij=1.0d0/dij(i,j)
            rij(1,i,j)=tij
            rij(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt,3*natm*mxp,0.0d0)
      do ipt=1,npts
c
         call aclear_dp(duijdi,3*natm*natm,0.0d0)
         call aclear_dp(duijdk,3*natm*natm,0.0d0)
         call aclear_dp(d2uijdi2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdk2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdidk,3*3*natm*natm,0.0d0)
c
c...     Construct uij and its derivatives
c
         do 95 i=1,natm
            if (gtype_num(i).eq.0) goto 95
            do k=1,3
               rix(k)=ra2_comp(ipt,i,k)
            enddo
            do 100 j=1,natm
               if (i.eq.j) goto 100
               if (gtype_num(j).eq.0) goto 100
               ri=ra2_val(ipt,i)
               rj=ra2_val(ipt,j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(ipt,j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
 100        continue ! j=1,i-1
  95     continue ! i=1,natm
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,natm,1.0d0)
         do 195 i=1,natm
            if (gtype_num(i).eq.0) goto 195
            do 200 j=1,natm
               if (i.eq.j) goto 200
               if (gtype_num(j).eq.0) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
               f1ij=(1.5d0-0.5d0*vij*vij)*vij
               f2ij=(1.5d0-0.5d0*f1ij*f1ij)*f1ij
               f3ij=(1.5d0-0.5d0*f2ij*f2ij)*f2ij
               sij=0.5d0*(1.0d0-f3ij)
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0-2.0d0*aij(i,j)*tij
               df1ijdv=1.5d0-1.5d0*vij*vij
               df2ijdf1=1.5d0-1.5d0*f1ij*f1ij
               df3ijdf2=1.5d0-1.5d0*f2ij*f2ij
               df2ijdv=df2ijdf1*df1ijdv
               df3ijdv=df3ijdf2*df2ijdv
               dsijdu=-0.5d0*df3ijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)
               d2f1ijdv2=-3.0d0*vij
               d2f2ijdv2=-3.0d0*f1ij*df1ijdv*df1ijdv
     &                   +df2ijdf1*d2f1ijdv2
               d2f3ijdv2=-3.0d0*f2ij*df2ijdv*df2ijdv
     &                   +df3ijdf2*d2f2ijdv2
               d2sijdu2=-0.5d0*d2f3ijdv2*sij*dvijdu*dvijdu
     &                  -0.5d0*df3ijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
c
 200        continue ! j=1,i-1
 195     continue ! i=1,natm
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,natm
            if(gtype_num(i).ne.0) Z=Z+pk(i)
         enddo ! i=1,natm
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*natm,0.0d0)
         do 295 i=1,natm
            if (gtype_num(i).eq.0) goto 295
            do 300 j=1,natm
               if (i.eq.j) goto 300
               if (gtype_num(j).eq.0) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,latm)=gZ(k,latm)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,i-1
 295     enddo !i=1,natm
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*natm*3*natm,0.0d0)
         lx=3*(latm-1)
         do 520 i=1,natm
            if (gtype_num(i).eq.0) goto 520
            ix=3*(i-1)
            do 500 j=1,natm
               if (j.eq.i) goto 500
               if (gtype_num(i).eq.0) goto 500
               jx=3*(j-1)
               do 510 k=1,natm
                  if (k.eq.i.or.k.eq.j) goto 510
                  if (gtype_num(k).eq.0) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,natm
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(latm)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,natm
            if (i.eq.latm) goto 400
            do k=1,3
               gwt(k,ipt,i)=dsijdj(k,latm,i)-gZ(k,i)
               gwt(k,ipt,latm)=gwt(k,ipt,latm)
     &                        +dsijdi(k,latm,i)+dsijdk(k,latm,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,ipt,latm)=gwt(k,ipt,latm)-gZ(k,latm)
         enddo
         do i=1,natm
            do k=1,3
               gwt(k,ipt,i)=wt(ipt)*gwt(k,ipt,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         call aclear_dp(g2wt,3*natm*3*natm,0.0d0)
         lx=3*(latm-1) 
         do 600 i=1,natm 
            if (i.eq.latm) goto 600
            if (gtype_num(i).eq.0) goto 600
            ix=3*(i-1)
            do 610 j=1,natm
               if (j.eq.i.or.j.eq.latm) goto 610
               if (gtype_num(j).eq.0) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &               *(dsijdi(jc,latm,j)+dsijdk(jc,latm,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,latm,i)*dsijdj(jc,latm,j)
     &               -dsijdj(ic,latm,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,latm,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,latm)
     &            -gZ(ic,latm)
     &            *(dsijdi(jc,latm,i)+dsijdk(jc,latm,i))
     &            +(d2sijdi2(ic,jc,latm,i)
     &             +d2sijdk2(ic,jc,latm,i)
     &             +d2sijdidk(ic,jc,latm,i)
     &             +d2sijdidk(jc,ic,latm,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,latm,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,latm,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,latm,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,latm)*gZ(jc,latm)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               hess(j,i)=hess(j,i)+g2wt(j,i)*wt(ipt)*xc_ept(ipt)
            enddo
         enddo
      enddo ! ipt

      end
c
c-----------------------------------------------------------------------
c
      subroutine beckewt_scr(ra2_val,latm,num_near,near_atom, arad, wt, 
     &     sk, awt, totwt, npts, mxp, natm, igrd)
C     ******************************************************************
C     *Description:						       *
C     *Calculates Becke weight at point                                *
c
c    screened version not currently in use
c
C     ******************************************************************
      implicit none
C     ******************************************************************
C     *Declarations                                                    *
C     *								       *
C     *Parameters						       *
INCLUDE(common/dft_parameters)
C     *In variables						       *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, num_near, natm
      REAL ra2_val(mxp,natm)
      REAL arad(natm)
      integer near_atom(num_near)
      integer latm
      integer igrd

C     *In/out variables						       *
      REAL wt(npts)
C     *Work space
      REAL awt(npts), totwt(npts), sk(npts)
C     *Local variables						       *
      integer i,j,ipt,imin,loopct,igrid,jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL uij1,fk1
      integer ipt1, ipt2, iatm, jatm
      REAL uij2,fk2
C     *Functions						       *
c     REAL srad
_IF(single)
      integer isamin
_ELSE
      integer idamin
_ENDIF

C     *End declarations                                                *
C     ******************************************************************
C
      if (num_near.le.1) return
c
      call aclear_dp(totwt,npts,0.0d0)
      do 10 i=1,num_near
         iatm = near_atom(i)

         igrid = gtype_num(iatm)

         call aclear_dp(sk,npts,1.0d0)
         do 20 j=1,num_near
            jatm = near_atom(j)
            if(j.eq.i) goto 20

            jgrid = gtype_num(jatm)

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(iatm,jatm)

            do ipt=1,npts
               if (ra2_val(ipt,iatm).le.arad(iatm).and.
     &             ra2_val(ipt,jatm).le.arad(jatm)) then
                  uij=diji*(ra2_val(ipt,iatm)-ra2_val(ipt,jatm))
                  fk=uij+aij*(1.0d0-(uij*uij))
                  fk=1.5d0*fk-0.5d0*fk*fk*fk
                  fk=1.5d0*fk-0.5d0*fk*fk*fk
                  fk=1.5d0*fk-0.5d0*fk*fk*fk
                  sk(ipt)=sk(ipt)*0.5d0*(1.0d0-fk)
               endif
            enddo
 20      continue
         do ipt=1,npts
            if (ra2_val(ipt,iatm).le.arad(iatm)) then
               totwt(ipt) = totwt(ipt) + sk(ipt)
            endif
         enddo
c        call daxpy(npts,1.0d0,sk,1,totwt,1)
         if(iatm.eq.latm)call dcopy(npts,sk,1,awt,1)
 10   continue
      
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine beckegwt_scr(ra2_val,ra2_comp,latm,
     &     num_near,near_atom,arad,wt,gwt,gwt_avail_sw,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:		                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights.                                          *
C  *********************************************************************

      implicit none
C  *********************************************************************
C  *Declarations                                                       *
C  *                                                                   *
C  *Parameters                                                         *
INCLUDE(common/dft_parameters)
C  *In variables		                                       *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm, num_near
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      REAL arad(natm)
      integer near_atom(num_near)
      integer latm
      integer igrd
      logical gwt_avail_sw
c
C  *In/out variables                                                   *
      REAL wt(npts), gwt(3,mxp,num_near)
C  *Work space                                                         *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,num_near),gZ(3,npts,num_near)
      REAL ti_ij(3,npts,num_near), tj_ij(3,npts,num_near)
C  *Local variables                                                    *
      integer i,j,ipt,loopct, iat, jat, igrid, jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL uij1,fk1
      integer ipt1, ipt2
      REAL uij2,fk2
      REAL dui_ijx,dui_ijy,dui_ijz
      REAL duj_ijx,duj_ijy,duj_ijz
      REAL vij,dvij,dfk,sk,dsk
C  *Functions                                                          *
c     REAL srad
C  *End declarations                                                   *
C  *********************************************************************
      
      if (num_near.le.1) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      call aclear_dp(Z,npts,0.0d0)
      call aclear_dp(gZ,3*npts*num_near,0.0d0)
      call aclear_dp(gpa,3*npts*num_near,0.0d0)

      do 10 i=1,num_near
         iat = near_atom(i)

         igrid = gtype_num(iat)

         call aclear_dp(pk,npts,1.0d0)
         do 20 j=1,num_near
            jat = near_atom(j)
            if(j.eq.i) goto 20

            jgrid = gtype_num(jat)

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(iat,jat)

            do ipt=1,npts
               if (ra2_val(ipt,iat).le.arad(iat).and.
     +             ra2_val(ipt,jat).le.arad(jat)) then
                  uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  vij=uij+aij*(1.0d0-(uij*uij))
                  dvij=1.0d0-2.0d0*aij*uij
c
                  dfk=(1.5d0-1.5d0*vij*vij)*dvij
                  fk=1.5d0*vij-0.5d0*vij*vij*vij
                  dfk=(1.5d0-1.5d0*fk*fk)*dfk
                  fk=1.5d0*fk-0.5d0*fk*fk*fk
                  dfk=(1.5d0-1.5d0*fk*fk)*dfk
                  fk=1.5d0*fk-0.5d0*fk*fk*fk
c
                  sk=0.5d0*(1.0d0-fk)
                  if (dabs(sk).lt.1.0d-15) then
                     dsk= 0.0d0
                  else
                     dsk=-0.5d0*dfk/sk
                  endif
c
                  pk(ipt)=pk(ipt)*sk
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               else
                  ti_ij(1,ipt,j)=0.0d0
                  ti_ij(2,ipt,j)=0.0d0
                  ti_ij(3,ipt,j)=0.0d0
                  tj_ij(1,ipt,j)=0.0d0
                  tj_ij(2,ipt,j)=0.0d0
                  tj_ij(3,ipt,j)=0.0d0
               endif
            enddo
 20      continue
c
         do ipt = 1, npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               Z(ipt) = Z(ipt) + pk(ipt)
            endif
         enddo
c
         if(iat.eq.latm) then
c
c...        Save Pa
c
            call dcopy(npts,pk,1,pa,1)
c
c...        Construct grad Pa
c
            do j=1,num_near
               if (j.ne.i) then
                  do ipt=1,npts
                     gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
                     gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
                     gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
                  enddo
               endif
            enddo
         endif
c
c...     Construct grad Z
c
         do j=1,num_near
            if (j.ne.i) then
               do ipt=1,npts
                  if (ra2_val(ipt,iat).le.arad(iat)) then
                     gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                     gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                     gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                     gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                     gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                     gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
                  endif
               enddo
            endif
         enddo
 10   continue
      
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      call aclear_dp(gwt,3*mxp*num_near,0.0d0)
      i = num_near
      do j=1,num_near-1
         if (latm.eq.near_atom(j)) i=j
      enddo
      do j=1,num_near
         if (j.ne.i) then
            do ipt=1,npts
               gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(1,ipt,j)
     +                         -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
               gwt(1,ipt,i) = gwt(1,ipt,i) - gwt(1,ipt,j)
               gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(2,ipt,j)
     +                         -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
               gwt(2,ipt,i) = gwt(2,ipt,i) - gwt(2,ipt,j)
               gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(3,ipt,j)
     +                         -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
               gwt(3,ipt,i) = gwt(3,ipt,i) - gwt(3,ipt,j)
            enddo
         endif
      enddo
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine beckeg2wt_scr(ra2_val_g,ra2_comp_g,ra2_val,ra2_comp,
     &     latm,num_near_g,near_atom_g,near_atoms,indx_atoms,arad,
     &     wt,gwt_g,gwt_avail_sw,
     &     aij_g, xij_g, rij_g, aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, gwt,
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, latm, nprt
      REAL ra2_val_g(mxp,natoms)
      REAL ra2_comp_g(mxp,natoms,3)
      REAL arad(natm)
      REAL xc_ept(mxp)
      integer num_near_g
      integer near_atom_g(num_near_g)
      integer igrd
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt_g(3,mxp,num_near_g)
      logical gwt_avail_sw
c
c...  Workspace
c
      REAL aij_g(num_near_g,num_near_g)
      REAL xij_g(3,num_near_g,num_near_g)
      REAL rij_g(2,num_near_g,num_near_g) ! rij(i,*) = dij**(-i)
      REAL aij(num_near_g,num_near_g)
      REAL xij(3,num_near_g,num_near_g)
      REAL rij(2,num_near_g,num_near_g) ! rij(i,*) stores dij**(-i)
      REAL uij(num_near_g,num_near_g)
      REAL pk(num_near_g)
      REAL ra2_val(num_near_g)
      REAL ra2_comp(num_near_g,3)
      REAL duijdi(3,num_near_g,num_near_g)
      REAL duijdj(3,num_near_g,num_near_g)
      REAL duijdk(3,num_near_g,num_near_g)
      REAL d2uijdi2(3,3,num_near_g,num_near_g)
      REAL d2uijdj2(3,3,num_near_g,num_near_g)
      REAL d2uijdk2(3,3,num_near_g,num_near_g)
      REAL d2uijdidj(3,3,num_near_g,num_near_g)
      REAL d2uijdidk(3,3,num_near_g,num_near_g)
      REAL d2uijdjdk(3,3,num_near_g,num_near_g)
      REAL gZ(3,num_near_g)
      REAL g2Z(3*num_near_g,3*num_near_g)
      REAL gwt(3,num_near_g)
      REAL g2wt(3*num_near_g,3*num_near_g)
      REAL dsijdi(3,num_near_g,num_near_g)
      REAL dsijdj(3,num_near_g,num_near_g)
      REAL dsijdk(3,num_near_g,num_near_g)
      REAL d2sijdi2(3,3,num_near_g,num_near_g)
      REAL d2sijdj2(3,3,num_near_g,num_near_g)
      REAL d2sijdk2(3,3,num_near_g,num_near_g)
      REAL d2sijdidj(3,3,num_near_g,num_near_g)
      REAL d2sijdidk(3,3,num_near_g,num_near_g)
      REAL d2sijdjdk(3,3,num_near_g,num_near_g)
      integer near_atoms(num_near_g) ! per pt the real near atoms
      integer indx_atoms(num_near_g) ! per pt the index in near_atom_g
c
c...  Local variables
c
      integer indx_latm
      integer num_near
      integer i, j, k, l, n, ipt
      integer ix, jx, kx, lx
      integer iy, jy
      integer iatm, jatm, ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, f1ij, f2ij, f3ij
      REAL df2ijdf1, df3ijdf2
      REAL dvijdu, df1ijdv, df2ijdv, df3ijdv
      REAL d2vijdu2, d2f1ijdv2, d2f2ijdv2, d2f3ijdv2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
c
c...  Functions
c
c     REAL srad
c
c...  Statement functions
c
c...  iks: small triangle i.e. without diagonal
c...  ikb: big   triangle i.e. with    diagonal
c...  ikm: ikb with min's and max's
c     integer iks, ikb, ikm
c     iks(i,j) = iky(i-1)+j
c     ikb(i,j) = iky(i)+j
c     ikm(i,j) = ikb(max(i,j),min(i,j))
c
c...  Data statements
c
c     data ikx/1,2,3,4,5,6,7,8,9/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      if (num_near_g.eq.0) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      do 5 i=1,num_near_g
         iatm = near_atom_g(i)
         igrid = gtype_num(iatm)
         if (igrid.eq.0) goto 5
         do 10 j=1,num_near_g
            if (i.eq.j) goto 10
            jatm = near_atom_g(j)
            jgrid = gtype_num(jatm)
            if (jgrid.eq.0) goto 10
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij_g(i,j)=tij
c
            xij_g(1,i,j)=atom_c(iatm,1)-atom_c(jatm,1)
            xij_g(2,i,j)=atom_c(iatm,2)-atom_c(jatm,2)
            xij_g(3,i,j)=atom_c(iatm,3)-atom_c(jatm,3)
c
            tij=1.0d0/dij(iatm,jatm)
            rij_g(1,i,j)=tij
            rij_g(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt_g,3*num_near_g*mxp,0.0d0)
      do 20 ipt=1,npts
c
c...     Construct the near atom list specific to the current point
c
         num_near = 0
         indx_latm = 0
         do i=1,num_near_g
            iatm = near_atom_g(i)
            if (ra2_val_g(ipt,iatm).le.arad(iatm).and.
     &          gtype_num(i).ne.0) then
               if (iatm.eq.latm) then
                  indx_latm = i
               else
                  num_near = num_near+1
                  near_atoms(num_near) = iatm
                  indx_atoms(num_near) = i
               endif
            endif
         enddo
         if (indx_latm.ne.0) then
            num_near = num_near+1
            near_atoms(num_near) = latm
            indx_atoms(num_near) = indx_latm
         endif
         if (num_near.le.1) goto 20
c
         n = num_near_g
         call aclear_dp(gwt,3*n,0.0d0)
         call aclear_dp(g2wt,3*n*3*n,0.0d0)
         call aclear_dp(duijdi,3*n*n,0.0d0)
         call aclear_dp(duijdk,3*n*n,0.0d0)
         call aclear_dp(d2uijdi2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdk2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdidk,3*3*n*n,0.0d0)
c
c...     Gather the data required for this point
c
         do i=1,num_near
            iatm=near_atoms(i)
            do k=1,3
               ra2_comp(i,k)=ra2_comp_g(ipt,iatm,k)
            enddo
            ra2_val(i)=ra2_val_g(ipt,iatm)
         enddo
         do i=1,num_near
            iatm=indx_atoms(i)
            do j=1,num_near
               jatm=indx_atoms(j)
               aij(j,i)=aij_g(jatm,iatm)
               do k=1,3
                  xij(k,j,i)=xij_g(k,jatm,iatm)
               enddo
               do k=1,2
                  rij(k,j,i)=rij_g(k,jatm,iatm)
               enddo
            enddo
         enddo
c
c...     Construct uij and its derivatives
c
         do i=1,num_near
            do k=1,3
               rix(k)=ra2_comp(i,k)
            enddo
            do 100 j=1,num_near
               if (i.eq.j) goto 100
               ri=ra2_val(i)
               rj=ra2_val(j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
 100        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,num_near,1.0d0)
         do i=1,num_near
            do 200 j=1,num_near
               if (i.eq.j) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
               f1ij=(1.5d0-0.5d0*vij*vij)*vij
               f2ij=(1.5d0-0.5d0*f1ij*f1ij)*f1ij
               f3ij=(1.5d0-0.5d0*f2ij*f2ij)*f2ij
               sij=0.5d0*(1.0d0-f3ij)
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0-2.0d0*aij(i,j)*tij
               df1ijdv=1.5d0-1.5d0*vij*vij
               df2ijdf1=1.5d0-1.5d0*f1ij*f1ij
               df3ijdf2=1.5d0-1.5d0*f2ij*f2ij
               df2ijdv=df2ijdf1*df1ijdv
               df3ijdv=df3ijdf2*df2ijdv
               dsijdu=-0.5d0*df3ijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)
               d2f1ijdv2=-3.0d0*vij
               d2f2ijdv2=-3.0d0*f1ij*df1ijdv*df1ijdv
     &                   +df2ijdf1*d2f1ijdv2
               d2f3ijdv2=-3.0d0*f2ij*df2ijdv*df2ijdv
     &                   +df3ijdf2*d2f2ijdv2
               d2sijdu2=-0.5d0*d2f3ijdv2*sij*dvijdu*dvijdu
     &                  -0.5d0*df3ijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
c
 200        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,num_near
            Z=Z+pk(i)
         enddo ! i=1,num_near
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*num_near_g,0.0d0)
         do i=1,num_near
            do 300 j=1,num_near
               if (i.eq.j) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,num_near)=gZ(k,num_near)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*num_near_g*3*num_near_g,0.0d0)
         lx=3*(num_near-1)
         do 520 i=1,num_near
            ix=3*(i-1)
            do 500 j=1,num_near
               if (j.eq.i) goto 500
               jx=3*(j-1)
               do 510 k=1,num_near
                  if (k.eq.i.or.k.eq.j) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,num_near
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(num_near)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,num_near-1
            do k=1,3
               gwt(k,i)=dsijdj(k,num_near,i)-gZ(k,i)
               gwt(k,num_near)=gwt(k,num_near)
     &                        +dsijdi(k,num_near,i)+dsijdk(k,num_near,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,num_near)=gwt(k,num_near)-gZ(k,num_near)
         enddo
         do i=1,num_near
            do k=1,3
               gwt(k,i)=wt(ipt)*gwt(k,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         lx=3*(num_near-1) 
         do 600 i=1,num_near-1
            ix=3*(i-1)
            do 610 j=1,num_near-1
               if (j.eq.i) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &               *(dsijdi(jc,num_near,j)+dsijdk(jc,num_near,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,num_near,i)*dsijdj(jc,num_near,j)
     &               -dsijdj(ic,num_near,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,num_near,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,num_near)
     &            -gZ(ic,num_near)
     &            *(dsijdi(jc,num_near,i)+dsijdk(jc,num_near,i))
     &            +(d2sijdi2(ic,jc,num_near,i)
     &             +d2sijdk2(ic,jc,num_near,i)
     &             +d2sijdidk(ic,jc,num_near,i)
     &             +d2sijdidk(jc,ic,num_near,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,num_near,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,num_near,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,num_near,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,num_near)*gZ(jc,num_near)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2wt(j,i)=g2wt(j,i)*wt(ipt)
            enddo
         enddo
c
c...     Scatter the results into their final location
c
         do i=1,num_near
            iatm=indx_atoms(i)
            do k=1,3
               gwt_g(k,ipt,iatm)=gwt(k,i)
            enddo
         enddo
         do i=1,num_near
            ix=3*(near_atom_g(indx_atoms(i))-1)
            iy=3*(i-1)
            do ic=1,3
               do j=1,num_near
                  jx=3*(near_atom_g(indx_atoms(j))-1)
                  jy=3*(j-1)
                  do jc=1,3
                     hess(jx+jc,ix+ic)=hess(jx+jc,ix+ic)
     &               +g2wt(jy+jc,iy+ic)*xc_ept(ipt)
                  enddo ! jc
               enddo ! j
            enddo ! ic
         enddo ! i
c        
 20   continue ! ipt

      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlwt(ra2_val,latm,wt,  
     &     sk, awt, totwt, npts, mxp, natm, igrd)

c **********************************************************************
c *   Huub van Dam, 1998                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates Becke weight [1] at a point but using Murray-s        *
c *   function for s(u) [2]. This function is based on the equation    *
c *   ds(u)/du = a (1-u**2)**m with m = 10 leading to a polynomial:    *
c *                                                                    *
c *       s(u) = sum(i=0:m) b_i u**(1+2i)                              *
c *                                                                    *
c *   where                                                            *
c *                                                                    *
c *                     ( m )  a                                       *
c *       b_i = (-1)**i (   ) ----,  i = 0,1,...,m                     *
c *                     ( i ) 1+2i                                     *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming,              *
c *        Mol.Phys. Vol 78 (1993) 997 (in particular pages 1003,1007) *
c *                                                                    *
c **********************************************************************

      implicit none
C **********************************************************************
C *   Declarations                                                     *
C *                                                                    *
C *   Parameters                                                       *
C *                                                                    *
INCLUDE(common/dft_parameters)
C *                                                                    *
C *   In variables		                                       *
C *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      integer latm
      integer igrd
C *                                                                    *
C *   In/out variables                                                 *
C *                                                                    *
      REAL wt(npts)
C *                                                                    *
C *   Work space                                                       *
C *                                                                    *
      REAL awt(npts), totwt(npts), sk(npts)
C *                                                                    *
C *   Local variables                                                  *
C *                                                                    *
      integer i,j,ipt,m, igrid, jgrid
      REAL radi,radj,ratij,uij,uij2,uijm,aij,fk,diji
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
C *                                                                    *
C *   Data statements
C *                                                                    *
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/

C *                                                                    *
C *   End declarations                                                 *
C *                                                                    *
C **********************************************************************
      call aclear_dp(totwt,npts,0.0d0)
      do 10 i=1,ngridcentres

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10
         call aclear_dp(sk,npts,1.0d0)

         do 20 j=1,ngridcentres

            if(j.eq.i) goto 20
            jgrid = gtype_num(j)
            if (jgrid.eq.0) goto 20

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(i,j)

            do ipt = 1, npts
               uij  = diji*(ra2_val(ipt,i)-ra2_val(ipt,j))
               uij  = uij + aij*(1.0d0-uij*uij)
               uij2 = uij*uij
               uijm = uij
               fk = 0.5d0 + b(0)*uijm
               do m = 1, mmax
                  uijm = uijm * uij2
                  fk = fk + b(m)*uijm
               enddo
               sk(ipt) = sk(ipt)*fk
            enddo

 20      continue
         call daxpy(npts,1.0d0,sk,1,totwt,1)
         if(i.eq.latm)call dcopy(npts,sk,1,awt,1)
 10   continue
      
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlgwt(ra2_val,ra2_comp,latm,wt,gwt,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm, igrd)

c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates Becke weight [1] and the gradient of the weight       *
c *   according to Johnson [2] at a point but using Murray-s           *
c *   function for s(u) [3]. This function is based on the equation    *
c *   ds(u)/du = a (1-u**2)**m with m = 10 leading to a polynomial:    *
c *                                                                    *
c *       s(u) = sum(i=0:m) b_i u**(1+2i)                              *
c *                                                                    *
c *   where                                                            *
c *                                                                    *
c *                     ( m )  a                                       *
c *       b_i = (-1)**i (   ) ----,  i = 0,1,...,m                     *
c *                     ( i ) 1+2i                                     *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. 98      *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming,              *
c *        Mol.Phys. Vol 78 (1993) 997 (in particular pages 1003,1007) *
c *                                                                    *
c **********************************************************************
      implicit none
C **********************************************************************
C *                                                                    *
C *   Declarations                                                     *
C *                                                                    *
C *   Parameters                                                       *
C *                                                                    *
INCLUDE(common/dft_parameters)
C *                                                                    *
C *   In variables                                                     *
C *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      integer latm
      integer igrd
c
C *   In/out variables                                                 *
C *                                                                    *
      REAL wt(npts), gwt(3,mxp,natm)
C *                                                                    *
C *   Work space                                                       *
C *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,natm),gZ(3,npts,natm)
      REAL ti_ij(3,npts,natm), tj_ij(3,npts,natm)
C *                                                                    *
C *   Local variables                                                  *
C *                                                                    *
      integer i,j,ipt,loopct, igrid, jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL dui_ijx,dui_ijy,dui_ijz
      REAL duj_ijx,duj_ijy,duj_ijz
      REAL vij,dvij,dfk,sk,dsk
      REAL vij2,vijm,vijm1
      integer m
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
C *                                                                    *
C *   Data statements
C *                                                                    *
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/
C *                                                                    *
C *   End declarations                                                 *
C **********************************************************************
      
      call aclear_dp(Z,npts,0.0d0)
      call aclear_dp(gZ,3*npts*natm,0.0d0)
      call aclear_dp(gpa,3*npts*natm,0.0d0)

      do 10 i=1,natm

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10

         call aclear_dp(pk,npts,1.0d0)
         do 20 j=1,natm
            if(j.eq.i) goto 20

            jgrid = gtype_num(j)
            if (jgrid.eq.0) then
               do ipt=1,npts
                  ti_ij(1,ipt,j)=0.0d0
                  ti_ij(2,ipt,j)=0.0d0
                  ti_ij(3,ipt,j)=0.0d0
                  tj_ij(1,ipt,j)=0.0d0
                  tj_ij(2,ipt,j)=0.0d0
                  tj_ij(3,ipt,j)=0.0d0
               enddo
               goto 20
            endif

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(i,j)

            do ipt=1,npts
               uij=diji*(ra2_val(ipt,i)-ra2_val(ipt,j))
               dui_ijx=-ra2_comp(ipt,i,1)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
               dui_ijy=-ra2_comp(ipt,i,2)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
               dui_ijz=-ra2_comp(ipt,i,3)*diji/ra2_val(ipt,i)
     +                 -(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
               duj_ijx= ra2_comp(ipt,j,1)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
               duj_ijy= ra2_comp(ipt,j,2)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
               duj_ijz= ra2_comp(ipt,j,3)*diji/ra2_val(ipt,j)
     +                 +(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
c
               vij=uij+aij*(1.0d0-(uij*uij))
               dvij=1.0d0-2.0d0*aij*uij
c
               vij2=vij*vij
               vijm=vij
c
               dfk = dvij*b(0)*(1.0d0-vij2)**mmax
               fk  = 0.5d0 + b(0)*vijm
               do m = 1, mmax
                  vijm  = vijm  * vij2
c                 vijm1 = vijm1 * vij2
                  fk  = fk  + b(m)*vijm
c                 dfk = dfk + (2*m+1)*b(m)*vijm1
               enddo
c              dfk = dfk*dvij
c
               sk=fk
               if (dabs(sk).lt.1.0d-15) then
                  dsk= 0.0d0
               else
                  dsk=dfk/sk
               endif
c
               pk(ipt)=pk(ipt)*sk
               ti_ij(1,ipt,j)=dsk*dui_ijx
               ti_ij(2,ipt,j)=dsk*dui_ijy
               ti_ij(3,ipt,j)=dsk*dui_ijz
               tj_ij(1,ipt,j)=dsk*duj_ijx
               tj_ij(2,ipt,j)=dsk*duj_ijy
               tj_ij(3,ipt,j)=dsk*duj_ijz
            enddo
 20      continue
c
         call daxpy(npts,1.0d0,pk,1,Z,1)
c
         if(i.eq.latm) then
c
c...        Save Pa
c
            call dcopy(npts,pk,1,pa,1)
c
c...        Construct grad Pa
c
            do j=1,natm
               if (j.ne.i) then
                  do ipt=1,npts
                     gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
                     gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
                     gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
                  enddo
               endif
            enddo
         endif
c
c...     Construct grad Z
c
         do j=1,natm
            if (j.ne.i) then
               do ipt=1,npts
                  gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                  gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                  gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                  gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                  gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                  gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
               enddo
            endif
         enddo
 10   continue
      
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      call aclear_dp(gwt,3*mxp*natm,0.0d0)
      do j=1,natm
         if (j.ne.latm) then
            do ipt=1,npts
               gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(1,ipt,j)
     +                         -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
               gwt(1,ipt,latm) = gwt(1,ipt,latm) - gwt(1,ipt,j)
               gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(2,ipt,j)
     +                         -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
               gwt(2,ipt,latm) = gwt(2,ipt,latm) - gwt(2,ipt,j)
               gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(3,ipt,j)
     +                         -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
               gwt(3,ipt,latm) = gwt(3,ipt,latm) - gwt(3,ipt,j)
            enddo
         endif
      enddo
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlg2wt(ra2_val,ra2_comp,latm,wt,gwt,
     &     aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, 
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates MHL weights at the grid points and also calculates      *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, nprt
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      REAL xc_ept(mxp)
      integer latm
      integer igrd
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt(3,mxp,natm)
c     REAL g2wt(mxp,3*natm,3*natm)
c
c...  Workspace
c
      REAL aij(natm,natm)
      REAL uij(natm,natm)
      REAL xij(3,natm,natm)
      REAL rij(2,natm,natm) ! rij(i,*) stores dij**(-i)
      REAL pk(natm)
      REAL duijdi(3,natm,natm)
      REAL duijdj(3,natm,natm)
      REAL duijdk(3,natm,natm)
      REAL d2uijdi2(3,3,natm,natm)
      REAL d2uijdj2(3,3,natm,natm)
      REAL d2uijdk2(3,3,natm,natm)
      REAL d2uijdidj(3,3,natm,natm)
      REAL d2uijdidk(3,3,natm,natm)
      REAL d2uijdjdk(3,3,natm,natm)
      REAL gZ(3,natm)
      REAL g2Z(3*natm,3*natm)
      REAL g2wt(3*natm,3*natm)
      REAL dsijdi(3,natm,natm)
      REAL dsijdj(3,natm,natm)
      REAL dsijdk(3,natm,natm)
      REAL d2sijdi2(3,3,natm,natm)
      REAL d2sijdj2(3,3,natm,natm)
      REAL d2sijdk2(3,3,natm,natm)
      REAL d2sijdidj(3,3,natm,natm)
      REAL d2sijdidk(3,3,natm,natm)
      REAL d2sijdjdk(3,3,natm,natm)
c
c...  Local variables
c
      integer i, j, k, l, ipt
      integer ix, jx, kx, lx
      integer ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, vij2, vijm
      REAL dvijdu, dfijdv, d2fijdv2
      REAL d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
      integer m
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
c
c...  Functions
c
c     REAL srad
c
c...  Data statements
c
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      do 5 i=1,natm
         igrid = gtype_num(i)
         if (igrid.eq.0) goto 5
         do 10 j=1,natm
            if (i.eq.j) goto 10
            jgrid = gtype_num(j)
            if (jgrid.eq.0) goto 10
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij(i,j)=tij
c
            xij(1,i,j)=atom_c(i,1)-atom_c(j,1)
            xij(2,i,j)=atom_c(i,2)-atom_c(j,2)
            xij(3,i,j)=atom_c(i,3)-atom_c(j,3)
c
            tij=1.0d0/dij(i,j)
            rij(1,i,j)=tij
            rij(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt,3*natm*mxp,0.0d0)
c     call aclear_dp(g2wt,3*natm*3*natm*mxp,0.0d0)
      do ipt=1,npts
c
         call aclear_dp(duijdi,3*natm*natm,0.0d0)
         call aclear_dp(duijdk,3*natm*natm,0.0d0)
         call aclear_dp(d2uijdi2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdk2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdidk,3*3*natm*natm,0.0d0)
c
c...     Construct uij and its derivatives
c
         do 95 i=1,natm
            if (gtype_num(i).eq.0) goto 95
            do k=1,3
               rix(k)=ra2_comp(ipt,i,k)
            enddo
            do 100 j=1,natm
               if (i.eq.j) goto 100
               if (gtype_num(j).eq.0) goto 100
               ri=ra2_val(ipt,i)
               rj=ra2_val(ipt,j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(ipt,j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
 100        continue ! j=1,i-1
  95     continue ! i=1,natm
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,natm,1.0d0)
         do 195 i=1,natm
            if (gtype_num(i).eq.0) goto 195
            do 200 j=1,natm
               if (i.eq.j) goto 200
               if (gtype_num(j).eq.0) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
               vij2=vij*vij
               vijm=vij
               sij=0.5d0 + b(0)*vijm
               do m = 1, mmax
                  vijm = vijm * vij2
                  sij  = sij  + b(m)*vijm
               enddo

               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0-2.0d0*aij(i,j)*tij
               dfijdv=b(0)*(1.0d0-vij2)**mmax
               dsijdu=dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)
               d2fijdv2=-2.0d0*mmax*b(0)*vij*(1.0d0-vij2)**(mmax-1)
               d2sijdu2=d2fijdv2*sij*dvijdu*dvijdu
     &                 +dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
c
 200        continue ! j=1,i-1
 195     continue ! i=1,natm
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,natm
            if (gtype_num(i).ne.0) Z=Z+pk(i)
         enddo ! i=1,natm
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*natm,0.0d0)
         do 295 i=1,natm
            if (gtype_num(i).eq.0) goto 295
            do 300 j=1,natm
               if (i.eq.j) goto 300
               if (gtype_num(j).eq.0) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,latm)=gZ(k,latm)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,i-1
 295     continue !i=1,natm
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*natm*3*natm,0.0d0)
         lx=3*(latm-1)
         do 520 i=1,natm
            if (gtype_num(i).eq.0) goto 520
            ix=3*(i-1)
            do 500 j=1,natm
               if (j.eq.i) goto 500
               if (gtype_num(j).eq.0) goto 500
               jx=3*(j-1)
               do 510 k=1,natm
                  if (k.eq.i.or.k.eq.j) goto 510
                  if (gtype_num(k).eq.0) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,natm
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(latm)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,natm
            if (i.eq.latm.or.gtype_num(i).eq.0) goto 400
            do k=1,3
               gwt(k,ipt,i)=dsijdj(k,latm,i)-gZ(k,i)
               gwt(k,ipt,latm)=gwt(k,ipt,latm)
     &                        +dsijdi(k,latm,i)+dsijdk(k,latm,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,ipt,latm)=gwt(k,ipt,latm)-gZ(k,latm)
         enddo
         do 410 i=1,natm
            if (gtype_num(i).eq.0) goto 410
            do k=1,3
               gwt(k,ipt,i)=wt(ipt)*gwt(k,ipt,i)
            enddo
 410     continue
c
c...     Calculate the hessians of the weights
c 
         call aclear_dp(g2wt,3*natm*3*natm,0.0d0)
         lx=3*(latm-1) 
         do 600 i=1,natm 
            if (i.eq.latm) goto 600
            if (gtype_num(i).eq.0) goto 600
            ix=3*(i-1)
            do 610 j=1,natm
               if (j.eq.i.or.j.eq.latm) goto 610
               if (gtype_num(j).eq.0) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &               *(dsijdi(jc,latm,j)+dsijdk(jc,latm,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,latm,i)*dsijdj(jc,latm,j)
     &               -dsijdj(ic,latm,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,latm,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,latm)
     &            -gZ(ic,latm)
     &            *(dsijdi(jc,latm,i)+dsijdk(jc,latm,i))
     &            +(d2sijdi2(ic,jc,latm,i)
     &             +d2sijdk2(ic,jc,latm,i)
     &             +d2sijdidk(ic,jc,latm,i)
     &             +d2sijdidk(jc,ic,latm,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,latm,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,latm,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,latm,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,latm)*gZ(jc,latm)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               hess(j,i)=hess(j,i)+g2wt(j,i)*wt(ipt)*xc_ept(ipt)
            enddo
         enddo
      enddo ! ipt
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlwt_scr(ra2_val,latm,num_near,near_atom, arad,wt,  
     &     sk, awt, totwt, npts, mxp, natm, igrd)

c **********************************************************************
c *   Huub van Dam, 1998                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates Becke weight [1] at a point but using Murray-s        *
c *   function for s(u) [2]. This function is based on the equation    *
c *   ds(u)/du = a (1-u**2)**m with m = 10 leading to a polynomial:    *
c *                                                                    *
c *       s(u) = sum(i=0:m) b_i u**(1+2i)                              *
c *                                                                    *
c *   where                                                            *
c *                                                                    *
c *                     ( m )  a                                       *
c *       b_i = (-1)**i (   ) ----,  i = 0,1,...,m                     *
c *                     ( i ) 1+2i                                     *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming,              *
c *        Mol.Phys. Vol 78 (1993) 997 (in particular pages 1003,1007) *
c *                                                                    *
c **********************************************************************

      implicit none
C **********************************************************************
C *   Declarations                                                     *
C *                                                                    *
C *   Parameters                                                       *
C *                                                                    *
INCLUDE(common/dft_parameters)
C *                                                                    *
C *   In variables		                                       *
C *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, num_near, natm
      REAL ra2_val(mxp,natoms)
      REAL arad(natm)
      integer near_atom(num_near)
      integer latm
      integer igrd
C *                                                                    *
C *   In/out variables                                                 *
C *                                                                    *
      REAL wt(npts)
C *                                                                    *
C *   Work space                                                       *
C *                                                                    *
      REAL awt(npts), totwt(npts), sk(npts)
C *                                                                    *
C *   Local variables                                                  *
C *                                                                    *
      integer i,j,ipt,m, iatm, jatm, igrid, jgrid
      REAL radi,radj,ratij,uij,uij2,uijm,aij,fk,diji
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
C *                                                                    *
C *   Data statements
C *                                                                    *
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/

C *                                                                    *
C *   End declarations                                                 *
C *                                                                    *
C **********************************************************************
      if (num_near.eq.0) then
c        nothing to be done here
         return
      endif
      call aclear_dp(totwt,npts,0.0d0)
      do 10 i=1,num_near
         iatm = near_atom(i)
         igrid = gtype_num(iatm)

         call aclear_dp(sk,npts,1.0d0)
         do 20 j=1,num_near
            jatm = near_atom(j)
            if(j.eq.i) goto 20

            jgrid = gtype_num(jatm)

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(iatm,jatm)

            do ipt = 1, npts
               if (ra2_val(ipt,iatm).le.arad(iatm).and.
     &             ra2_val(ipt,jatm).le.arad(jatm)) then
                  uij  = diji*(ra2_val(ipt,iatm)-ra2_val(ipt,jatm))
                  uij  = uij + aij*(1.0d0-uij*uij)
                  uij2 = uij*uij
                  uijm = uij
                  fk = 0.5d0 + b(0)*uijm
                  do m = 1, mmax
                     uijm = uijm * uij2
                     fk = fk + b(m)*uijm
                  enddo
                  sk(ipt) = sk(ipt)*fk
               endif
            enddo

 20      continue
         do ipt=1,npts
            if (ra2_val(ipt,iatm).le.arad(iatm)) then
               totwt(ipt) = totwt(ipt) + sk(ipt)
            endif
         enddo
         if(iatm.eq.latm)call dcopy(npts,sk,1,awt,1)
 10   continue
c     
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlgwt_scr(ra2_val,ra2_comp,latm,
     &     num_near,near_atom,arad,wt,gwt,gwt_avail_sw,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm, igrd)

c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates Becke weight [1] and the gradient of the weight       *
c *   according to Johnson [2] at a point but using Murray-s           *
c *   function for s(u) [3]. This function is based on the equation    *
c *   ds(u)/du = a (1-u**2)**m with m = 10 leading to a polynomial:    *
c *                                                                    *
c *       s(u) = sum(i=0:m) b_i u**(1+2i)                              *
c *                                                                    *
c *   where                                                            *
c *                                                                    *
c *                     ( m )  a                                       *
c *       b_i = (-1)**i (   ) ----,  i = 0,1,...,m                     *
c *                     ( i ) 1+2i                                     *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. 98      *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming,              *
c *        Mol.Phys. Vol 78 (1993) 997 (in particular pages 1003,1007) *
c *                                                                    *
c **********************************************************************
      implicit none
C **********************************************************************
C *                                                                    *
C *   Declarations                                                     *
C *                                                                    *
C *   Parameters                                                       *
C *                                                                    *
INCLUDE(common/dft_parameters)
C *                                                                    *
C *   In variables                                                     *
C *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, num_near, natm
      REAL ra2_val(mxp,natm)
      REAL ra2_comp(mxp,natoms,3)
      REAL arad(natm)
      integer near_atom(num_near)
      integer latm
      integer igrd
      logical gwt_avail_sw
c
C *   In/out variables                                                 *
C *                                                                    *
      REAL wt(npts), gwt(3,mxp,num_near)
C *                                                                    *
C *   Work space                                                       *
C *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,num_near),gZ(3,npts,num_near)
      REAL ti_ij(3,npts,num_near), tj_ij(3,npts,num_near)
C *                                                                    *
C *   Local variables                                                  *
C *                                                                    *
      integer i,j,ipt,loopct, iat, jat, igrid, jgrid
      REAL radi,radj,ratij,uij,aij,fk,diji
      REAL dui_ijx,dui_ijy,dui_ijz
      REAL duj_ijx,duj_ijy,duj_ijz
      REAL vij,dvij,dfk,sk,dsk
      REAL vij2,vijm,vijm1
      integer m
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
C *                                                                    *
C *   Data statements
C *                                                                    *
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/
C *                                                                    *
C *   End declarations                                                 *
C **********************************************************************
      if (num_near.eq.0) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c 
      call aclear_dp(Z,npts,0.0d0)
      call aclear_dp(gZ,3*npts*num_near,0.0d0)
      call aclear_dp(gpa,3*npts*num_near,0.0d0)

      do 10 i=1,num_near
         iat = near_atom(i)

         igrid = gtype_num(iat)

         call aclear_dp(pk,npts,1.0d0)
         do 20 j=1,num_near
            jat = near_atom(j)
            if(j.eq.i) goto 20

            jgrid = gtype_num(jat)

            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            ratij=radi/radj
            uij=(ratij-1.0d0)/(ratij+1.0d0)
            aij=uij/((uij*uij)-1.0d0)
            if(abs(aij).gt.0.5d0)aij=0.5d0*abs(aij)/aij

            diji = 1.0d0/dij(iat,jat)

            do ipt=1,npts
               if (ra2_val(ipt,iat).le.arad(iat).and.
     +             ra2_val(ipt,jat).le.arad(jat)) then
                  uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  vij=uij+aij*(1.0d0-(uij*uij))
                  dvij=1.0d0-2.0d0*aij*uij
c
                  vij2=vij*vij
                  vijm=vij
c
                  dfk = dvij*b(0)*(1.0d0-vij2)**mmax
                  fk  = 0.5d0 + b(0)*vijm
                  do m = 1, mmax
                     vijm  = vijm  * vij2
                     fk  = fk  + b(m)*vijm
                  enddo
c
                  sk=fk
                  if (dabs(sk).lt.1.0d-15) then
                     dsk= 0.0d0
                  else
                     dsk=dfk/sk
                  endif
c
                  pk(ipt)=pk(ipt)*sk
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               else
                  ti_ij(1,ipt,j)=0.0d0
                  ti_ij(2,ipt,j)=0.0d0
                  ti_ij(3,ipt,j)=0.0d0
                  tj_ij(1,ipt,j)=0.0d0
                  tj_ij(2,ipt,j)=0.0d0
                  tj_ij(3,ipt,j)=0.0d0
               endif
            enddo
 20      continue
c
         do ipt = 1, npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               Z(ipt) = Z(ipt) + pk(ipt)
            endif
         enddo
c
         if(iat.eq.latm) then
c
c...        Save Pa
c
            call dcopy(npts,pk,1,pa,1)
c
c...        Construct grad Pa
c
            do j=1,num_near
               if (j.ne.i) then
                  do ipt=1,npts
                     gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
                     gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
                     gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
                  enddo
               endif
            enddo
         endif
c
c...     Construct grad Z
c
         do j=1,num_near
            if (j.ne.i) then
               do ipt=1,npts
                  if (ra2_val(ipt,iat).le.arad(iat)) then
                     gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                     gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                     gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                     gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                     gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                     gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
                  endif
               enddo
            endif
         enddo
 10   continue
      
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      call aclear_dp(gwt,3*mxp*num_near,0.0d0)
      i = num_near
      do j=1,num_near-1
         if (latm.eq.near_atom(j)) i=j
      enddo
      do j=1,num_near
         if (j.ne.i) then
            do ipt=1,npts
               gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(1,ipt,j)
     +                         -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
               gwt(1,ipt,i) = gwt(1,ipt,i) - gwt(1,ipt,j)
               gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(2,ipt,j)
     +                         -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
               gwt(2,ipt,i) = gwt(2,ipt,i) - gwt(2,ipt,j)
               gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(3,ipt,j)
     +                         -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
               gwt(3,ipt,i) = gwt(3,ipt,i) - gwt(3,ipt,j)
            enddo
         endif
      enddo
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhlg2wt_scr(ra2_val_g,ra2_comp_g,ra2_val,ra2_comp,
     &     latm,num_near_g,near_atom_g,near_atoms,indx_atoms,arad,
     &     wt,gwt_g,gwt_avail_sw,
     &     aij_g, xij_g, rij_g, aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, gwt,
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, latm, igrd, nprt
      REAL ra2_val_g(mxp,natoms)
      REAL ra2_comp_g(mxp,natoms,3)
      REAL arad(natm)
      REAL xc_ept(mxp)
      integer num_near_g
      integer near_atom_g(num_near_g)
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt_g(3,mxp,num_near_g)
      logical gwt_avail_sw
c
c...  Workspace
c
      REAL aij_g(num_near_g,num_near_g)
      REAL xij_g(3,num_near_g,num_near_g)
      REAL rij_g(2,num_near_g,num_near_g) ! rij(i,*) = dij**(-i)
      REAL aij(num_near_g,num_near_g)
      REAL xij(3,num_near_g,num_near_g)
      REAL rij(2,num_near_g,num_near_g) ! rij(i,*) stores dij**(-i)
      REAL uij(num_near_g,num_near_g)
      REAL pk(num_near_g)
      REAL ra2_val(num_near_g)
      REAL ra2_comp(num_near_g,3)
      REAL duijdi(3,num_near_g,num_near_g)
      REAL duijdj(3,num_near_g,num_near_g)
      REAL duijdk(3,num_near_g,num_near_g)
      REAL d2uijdi2(3,3,num_near_g,num_near_g)
      REAL d2uijdj2(3,3,num_near_g,num_near_g)
      REAL d2uijdk2(3,3,num_near_g,num_near_g)
      REAL d2uijdidj(3,3,num_near_g,num_near_g)
      REAL d2uijdidk(3,3,num_near_g,num_near_g)
      REAL d2uijdjdk(3,3,num_near_g,num_near_g)
      REAL gZ(3,num_near_g)
      REAL g2Z(3*num_near_g,3*num_near_g)
      REAL gwt(3,num_near_g)
      REAL g2wt(3*num_near_g,3*num_near_g)
      REAL dsijdi(3,num_near_g,num_near_g)
      REAL dsijdj(3,num_near_g,num_near_g)
      REAL dsijdk(3,num_near_g,num_near_g)
      REAL d2sijdi2(3,3,num_near_g,num_near_g)
      REAL d2sijdj2(3,3,num_near_g,num_near_g)
      REAL d2sijdk2(3,3,num_near_g,num_near_g)
      REAL d2sijdidj(3,3,num_near_g,num_near_g)
      REAL d2sijdidk(3,3,num_near_g,num_near_g)
      REAL d2sijdjdk(3,3,num_near_g,num_near_g)
      integer near_atoms(num_near_g) ! per pt the real near atoms
      integer indx_atoms(num_near_g) ! per pt the index in near_atom_g
c
c...  Local variables
c
      integer indx_latm
      integer num_near
      integer i, j, k, l, n, ipt
      integer ix, jx, kx, lx
      integer iy, jy
      integer iatm, jatm, ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, vij2, vijm
      REAL dvijdu, dfijdv, d2fijdv2
      REAL d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
      integer m
      integer mmax
      parameter (mmax=10)
      REAL b(0:mmax)
c
c...  Functions
c
c     REAL srad
c
c...  Data statements
c
      data b/
     +      -0.18500690460205271d+01,
     +       0.61668968200684242d+01,
     +      -0.16650621414184744d+02,
     +       0.31715469360351893d+02,
     +      -0.43168277740478963d+02,
     +       0.42383399963379354d+02,
     +      -0.29885730743408516d+02,
     +       0.14800552368164217d+02,
     +      -0.48972415924072781d+01,
     +       0.97372055053711959d+00,
     +      -0.88098526000977478d-01/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      if (num_near_g.eq.0) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      do 5 i=1,num_near_g
         iatm = near_atom_g(i)
         igrid = gtype_num(iatm)
         if (igrid.eq.0) goto 5
         do 10 j=1,num_near_g
            if (i.eq.j) goto 10
            jatm = near_atom_g(j)
            jgrid = gtype_num(jatm)
            if (igrid.eq.0) goto 10
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij_g(i,j)=tij
c
            xij_g(1,i,j)=atom_c(iatm,1)-atom_c(jatm,1)
            xij_g(2,i,j)=atom_c(iatm,2)-atom_c(jatm,2)
            xij_g(3,i,j)=atom_c(iatm,3)-atom_c(jatm,3)
c
            tij=1.0d0/dij(iatm,jatm)
            rij_g(1,i,j)=tij
            rij_g(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt_g,3*num_near_g*mxp,0.0d0)
      do 20 ipt=1,npts
c
c...     Construct the near atom list specific to the current point
c
         indx_latm = 0
         num_near = 0
         do i=1,num_near_g
            iatm = near_atom_g(i)
            if (ra2_val_g(ipt,iatm).le.arad(iatm).and.
     &          gtype_num(i).ne.0) then
               if (iatm.eq.latm) then
                  indx_latm = i
               else
                  num_near = num_near+1
                  near_atoms(num_near) = iatm
                  indx_atoms(num_near) = i
               endif
            endif
         enddo
         if (indx_latm.ne.0) then
            num_near = num_near+1
            near_atoms(num_near) = latm
            indx_atoms(num_near) = indx_latm
         endif
         if (num_near.le.1) goto 20
c
         n = num_near_g
         call aclear_dp(gwt,3*n,0.0d0)
         call aclear_dp(g2wt,3*n*3*n,0.0d0)
         call aclear_dp(duijdi,3*n*n,0.0d0)
         call aclear_dp(duijdk,3*n*n,0.0d0)
         call aclear_dp(d2uijdi2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdk2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdidk,3*3*n*n,0.0d0)
c
c...     Gather the data required for this point
c
         do i=1,num_near
            iatm=near_atoms(i)
            do k=1,3
               ra2_comp(i,k)=ra2_comp_g(ipt,iatm,k)
            enddo
            ra2_val(i)=ra2_val_g(ipt,iatm)
         enddo
         do i=1,num_near
            iatm=indx_atoms(i)
            do j=1,num_near
               jatm=indx_atoms(j)
               aij(j,i)=aij_g(jatm,iatm)
               do k=1,3
                  xij(k,j,i)=xij_g(k,jatm,iatm)
               enddo
               do k=1,2
                  rij(k,j,i)=rij_g(k,jatm,iatm)
               enddo
            enddo
         enddo
c
c...     Construct uij and its derivatives
c
         do i=1,num_near
            do k=1,3
               rix(k)=ra2_comp(i,k)
            enddo
            do 100 j=1,num_near
               if (i.eq.j) goto 100
               ri=ra2_val(i)
               rj=ra2_val(j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
 100        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,num_near,1.0d0)
         do i=1,num_near
            do 200 j=1,num_near
               if (i.eq.j) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
               vij2=vij*vij
               vijm=vij
               sij=0.5d0 + b(0)*vijm
               do m = 1, mmax
                  vijm = vijm * vij2
                  sij  = sij  + b(m)*vijm
               enddo
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0-2.0d0*aij(i,j)*tij
               dfijdv=b(0)*(1.0d0-vij2)**mmax
               dsijdu=dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)
               d2fijdv2=-2.0d0*b(0)*mmax*vij*(1.0d0-vij2)**(mmax-1)
               d2sijdu2=d2fijdv2*sij*dvijdu*dvijdu
     &                 +dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
c
 200        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,num_near
            Z=Z+pk(i)
         enddo ! i=1,num_near
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*num_near_g,0.0d0)
         do i=1,num_near
            do 300 j=1,num_near
               if (i.eq.j) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,num_near)=gZ(k,num_near)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*num_near_g*3*num_near_g,0.0d0)
         lx=3*(num_near-1)
         do 520 i=1,num_near
            ix=3*(i-1)
            do 500 j=1,num_near
               if (j.eq.i) goto 500
               jx=3*(j-1)
               do 510 k=1,num_near
                  if (k.eq.i.or.k.eq.j) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,num_near
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(num_near)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,num_near-1
            do k=1,3
               gwt(k,i)=dsijdj(k,num_near,i)-gZ(k,i)
               gwt(k,num_near)=gwt(k,num_near)
     &                        +dsijdi(k,num_near,i)+dsijdk(k,num_near,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,num_near)=gwt(k,num_near)-gZ(k,num_near)
         enddo
         do i=1,num_near
            do k=1,3
               gwt(k,i)=wt(ipt)*gwt(k,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         lx=3*(num_near-1) 
         do 600 i=1,num_near-1
            ix=3*(i-1)
            do 610 j=1,num_near-1
               if (j.eq.i) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &               *(dsijdi(jc,num_near,j)+dsijdk(jc,num_near,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,num_near,i)*dsijdj(jc,num_near,j)
     &               -dsijdj(ic,num_near,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,num_near,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,num_near)
     &            -gZ(ic,num_near)
     &            *(dsijdi(jc,num_near,i)+dsijdk(jc,num_near,i))
     &            +(d2sijdi2(ic,jc,num_near,i)
     &             +d2sijdk2(ic,jc,num_near,i)
     &             +d2sijdidk(ic,jc,num_near,i)
     &             +d2sijdidk(jc,ic,num_near,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,num_near,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,num_near,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,num_near,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,num_near)*gZ(jc,num_near)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2wt(j,i)=g2wt(j,i)*wt(ipt)
            enddo
         enddo
c
c...     Scatter the results into their final location
c
         do i=1,num_near
            iatm=indx_atoms(i)
            do k=1,3
               gwt_g(k,ipt,iatm)=gwt(k,i)
            enddo
         enddo
         do i=1,num_near
            ix=3*(near_atom_g(indx_atoms(i))-1)
            iy=3*(i-1)
            do ic=1,3
               do j=1,num_near
                  jx=3*(near_atom_g(indx_atoms(j))-1)
                  jy=3*(j-1)
                  do jc=1,3
                     hess(jx+jc,ix+ic)=hess(jx+jc,ix+ic)
     &               +g2wt(jy+jc,iy+ic)*xc_ept(ipt)
                  enddo ! jc
               enddo ! j
            enddo ! ic
         enddo ! i
c        
 20   continue ! ipt
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfwt(ra2_val,latm,wt,  
     &     sk, awt, totwt, npts, mxp, natm,
     &     rshell, rnear)
c
C***********************************************************************
C Description:	
C   Calculates weight at point                                          
c   Following Stratmann, Scuseria, Frisch, CPL, 257 (1996) p213
c
c   version 1 : old becke loop structure, with no screening
c
      implicit none

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      integer latm
c
      REAL rshell, rnear
c
C     In/out variables
c
      REAL wt(npts)
c
C     Work space
c
      REAL awt(npts), totwt(npts), sk(npts)
c
C     Local variables	
c
      integer i,j,ipt,igrid,jgrid
      REAL radi,radj,ratij,uij,fk,diji, aij
      REAL a, uoa, uoa2, uoa3, uoa5, uoa7
c
      a = 0.60d0
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
      if(rshell .lt. 0.5d0*(1.0d0-a)*rnear)return
c
      call aclear_dp(totwt,npts,0.0d0)

      do 10 i=1,ngridcentres

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10

         call aclear_dp(sk,npts,1.0d0)
         do 20 j=1,ngridcentres
            jgrid = gtype_num(j)
            if (jgrid.eq.0) goto 20

            if(j.ne.i) then

               diji = 1.0d0/dij(i,j)

               do ipt=1,npts

                  uij=diji*(ra2_val(ipt,i)-ra2_val(ipt,j))

                  if(uij .lt. -a)then
c                    sk values unchanged
                  else if(uij .gt. a)then
                     sk(ipt) = 0.0d0
                  else
                     uoa = uij / a
                     uoa2 = uoa*uoa
                     uoa3 = uoa*uoa2
                     uoa5 = uoa3*uoa2
                     uoa7 = uoa5*uoa2
                     fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                     fk = fk * 0.0625d0
                     sk(ipt)=sk(ipt)*0.5d0*(1.0d0-fk)
                  endif
               enddo
            endif
 20      continue
         call daxpy(npts,1.0d0,sk,1,totwt,1)
         if(i.eq.latm)call dcopy(npts,sk,1,awt,1)
 10   continue
c     
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfgwt(ra2_val,ra2_comp,latm,wt,gwt,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm,
     &     rshell, rnear )
c
c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates the Becke weight [1] and the gradient of the weight   *
c *   according to Johnson et al. [2] but using the Stratmann-s et al. *
c *   [3] function for s(u). In this subroutine no screening is used   *
c *   to save computation (see ssfgwt_scr for that).                   *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. Vol. 98 *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Achieving linear scaling in exchange-correlation density    *
c *        functional quadratures", R.E. Stratmann, G.E. Scuseria,     *
c *        M.J. Frisch, Chem.Phys.Lett. Vol. 257 (1996) 213            *
c *                                                                    *
c **********************************************************************
      implicit none
c **********************************************************************
c *                                                                    *
c *   Declarations                                                     *
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
INCLUDE(common/dft_parameters)
c *                                                                    *
c *   In variables                                                     *
c *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      integer latm, natm
      REAL rshell, rnear
c *                                                                    *
c *   In/out variables                                                 *
c *                                                                    *
      REAL wt(npts), gwt(3,mxp,natm)
c *                                                                    *
c *   Work space                                                       *
c *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,natm),gZ(3,npts,natm)
      REAL ti_ij(3,npts,natm), tj_ij(3,npts,natm)
c *                                                                    *
c *   Local variables                                                  *
c *                                                                    *
      integer i,j,ipt, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji, aij
      REAL a, uoa, uoa2, uoa3, uoa4, uoa5, uoa6, uoa7
      REAL vij, dvij, dfk, dsk, sk
      REAL dui_ijx, dui_ijy, dui_ijz, duj_ijx, duj_ijy, duj_ijz
c
      a = 0.60d0
c
c     If the shell is within sphere defined by eq 15. in SSF [3]
c     there is no change to the weight for any of the points
c     in this batch, this means that the gradient of these weights
c     will be zero.
c
      if (rshell .lt. 0.5d0*(1.0d0-a)*rnear) then
         call aclear_dp(gwt,3*mxp*natm,0.0d0)
         return
      endif
c
      call aclear_dp(Z,npts,0.0d0)
      call aclear_dp(gZ,3*npts*natm,0.0d0)
      call aclear_dp(gpa,3*npts*natm,0.0d0)

      do 10 i=1,natm

         igrid = gtype_num(i)
         if (igrid.eq.0) goto 10

         call aclear_dp(pk,npts,1.0d0)
         do 20 j=1,natm
            jgrid = gtype_num(j)
            if (jgrid.eq.0) then
               do ipt=1,npts
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
               enddo
               goto 20
            endif

            if(j.ne.i) then

               diji = 1.0d0/dij(i,j)

               do ipt=1,npts

                  uij=diji*(ra2_val(ipt,i)-ra2_val(ipt,j))
                  vij=uij

                  if(vij .lt. -a)then
c
c...                 sk =1.0 therefore pk remains unchanged
c...                 dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                     ti_ij(1,ipt,j) = 0.0d0
                     ti_ij(2,ipt,j) = 0.0d0
                     ti_ij(3,ipt,j) = 0.0d0
                     tj_ij(1,ipt,j) = 0.0d0
                     tj_ij(2,ipt,j) = 0.0d0
                     tj_ij(3,ipt,j) = 0.0d0
c
                  else if(vij .gt. a)then
c
c...                 sk =0.0 therefore pk drops to zero
c...                 dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                     pk(ipt) = 0.0d0
                     ti_ij(1,ipt,j) = 0.0d0
                     ti_ij(2,ipt,j) = 0.0d0
                     ti_ij(3,ipt,j) = 0.0d0
                     tj_ij(1,ipt,j) = 0.0d0
                     tj_ij(2,ipt,j) = 0.0d0
                     tj_ij(3,ipt,j) = 0.0d0
c
                  else
c
                     dui_ijx=-ra2_comp(ipt,i,1)*diji/ra2_val(ipt,i)
     +                       -(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
                     dui_ijy=-ra2_comp(ipt,i,2)*diji/ra2_val(ipt,i)
     +                       -(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
                     dui_ijz=-ra2_comp(ipt,i,3)*diji/ra2_val(ipt,i)
     +                       -(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
                     duj_ijx= ra2_comp(ipt,j,1)*diji/ra2_val(ipt,j)
     +                       +(atom_c(i,1)-atom_c(j,1))*uij*diji*diji
                     duj_ijy= ra2_comp(ipt,j,2)*diji/ra2_val(ipt,j)
     +                       +(atom_c(i,2)-atom_c(j,2))*uij*diji*diji
                     duj_ijz= ra2_comp(ipt,j,3)*diji/ra2_val(ipt,j)
     +                       +(atom_c(i,3)-atom_c(j,3))*uij*diji*diji
c
                     dvij=1.0d0
c
                     uoa = vij / a
                     uoa2 = uoa*uoa
                     uoa3 = uoa*uoa2
c                    uoa4 = uoa2*uoa2
                     uoa5 = uoa3*uoa2
c                    uoa6 = uoa4*uoa2
                     uoa7 = uoa5*uoa2
c
                     fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                     fk = fk * 0.0625d0
                     sk = 0.5d0*(1.0d0-fk)
                     pk(ipt)=pk(ipt)*sk
c
                     dfk = dvij*35.0d0*(1.0d0-uoa2)**3
                     dfk = dfk * 0.0625d0 / a
                     if (sk.lt.1.0d-15) then
                        dsk = 0.0d0
                     else
                        dsk =-0.5d0*dfk/sk
                     endif
                     ti_ij(1,ipt,j)=dsk*dui_ijx
                     ti_ij(2,ipt,j)=dsk*dui_ijy
                     ti_ij(3,ipt,j)=dsk*dui_ijz
                     tj_ij(1,ipt,j)=dsk*duj_ijx
                     tj_ij(2,ipt,j)=dsk*duj_ijy
                     tj_ij(3,ipt,j)=dsk*duj_ijz
                  endif
               enddo
            endif
 20      continue
c
         call daxpy(npts,1.0d0,pk,1,Z,1)
c
         if(i.eq.latm) then
c
c...        Save Pa
c
            call dcopy(npts,pk,1,pa,1)
c
c...        Construct grad Pa
c
            do j=1,natm
               if (j.ne.i) then
                  do ipt=1,npts
                     gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
                     gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
                     gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
                  enddo
               endif
            enddo
         endif
c
c...     Construct grad Z
c
         do j=1,natm
            if (j.ne.i) then
               do ipt=1,npts
                  gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                  gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                  gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                  gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                  gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                  gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
               enddo
            endif
         enddo
 10   continue
c
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      call aclear_dp(gwt(1,1,latm),3*mxp,0.0d0)
      do j=1,natm
         if (j.ne.latm) then
            do ipt=1,npts
               gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(1,ipt,j)
     +                         -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
               gwt(1,ipt,latm) = gwt(1,ipt,latm) - gwt(1,ipt,j)
               gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(2,ipt,j)
     +                         -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
               gwt(2,ipt,latm) = gwt(2,ipt,latm) - gwt(2,ipt,j)
               gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(
     +                          gpa(3,ipt,j)
     +                         -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
               gwt(3,ipt,latm) = gwt(3,ipt,latm) - gwt(3,ipt,j)
            enddo
         endif
      enddo
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfg2wt(ra2_val,ra2_comp,latm,wt,gwt,
     &     uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, 
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, rshell, rnear)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates MHL weights at the grid points and also calculates      *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, nprt
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      REAL xc_ept(mxp)
      REAL rshell, rnear
      integer latm
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt(3,mxp,natm)
c     REAL g2wt(mxp,3*natm,3*natm)
c
c...  Workspace
c
      REAL uij(natm,natm)
      REAL xij(3,natm,natm)
      REAL rij(2,natm,natm) ! rij(i,*) stores dij**(-i)
      REAL pk(natm)
      REAL duijdi(3,natm,natm)
      REAL duijdj(3,natm,natm)
      REAL duijdk(3,natm,natm)
      REAL d2uijdi2(3,3,natm,natm)
      REAL d2uijdj2(3,3,natm,natm)
      REAL d2uijdk2(3,3,natm,natm)
      REAL d2uijdidj(3,3,natm,natm)
      REAL d2uijdidk(3,3,natm,natm)
      REAL d2uijdjdk(3,3,natm,natm)
      REAL gZ(3,natm)
      REAL g2Z(3*natm,3*natm)
      REAL g2wt(3*natm,3*natm)
      REAL dsijdi(3,natm,natm)
      REAL dsijdj(3,natm,natm)
      REAL dsijdk(3,natm,natm)
      REAL d2sijdi2(3,3,natm,natm)
      REAL d2sijdj2(3,3,natm,natm)
      REAL d2sijdk2(3,3,natm,natm)
      REAL d2sijdidj(3,3,natm,natm)
      REAL d2sijdidk(3,3,natm,natm)
      REAL d2sijdjdk(3,3,natm,natm)
c
c...  Local variables
c
      integer i, j, k, l, ipt
      integer ix, jx, kx, lx
      integer ic, jc, kc
      integer igrid, jgrid
      REAL vij, vij2, vij3, vij5, vij7
      REAL fij, dfijdv, d2fijdv2
      REAL dvijdu, d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
c
c...  Parameters
c
      REAL a
      parameter(a=0.60d0)
c
c     ==================================================================
c
      if (rshell.lt.0.5d0*(1.0d0-a)*rnear) then
         call aclear_dp(gwt,3*natm*mxp,0.0d0)
c        call aclear_dp(g2wt,3*natm*3*natm*mxp,0.0d0)
         return
      endif
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      do 5 i=1,natm
         if (gtype_num(i).eq.0) goto 5
         do 10 j=1,natm
            if (i.eq.j) goto 10
            if (gtype_num(j).eq.0) goto 10
            xij(1,i,j)=atom_c(i,1)-atom_c(j,1)
            xij(2,i,j)=atom_c(i,2)-atom_c(j,2)
            xij(3,i,j)=atom_c(i,3)-atom_c(j,3)
c
            tij=1.0d0/dij(i,j)
            rij(1,i,j)=tij
            rij(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt,3*natm*mxp,0.0d0)
c     call aclear_dp(g2wt,3*natm*3*natm*mxp,0.0d0)
      do ipt=1,npts
c
         call aclear_dp(duijdi,3*natm*natm,0.0d0)
         call aclear_dp(duijdk,3*natm*natm,0.0d0)
         call aclear_dp(d2uijdi2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdk2,3*3*natm*natm,0.0d0)
         call aclear_dp(d2uijdidk,3*3*natm*natm,0.0d0)
c
c...     Construct uij and its derivatives
c
         do 95 i=1,natm
            if (gtype_num(i).eq.0) goto 95
            do k=1,3
               rix(k)=ra2_comp(ipt,i,k)
            enddo
            do 100 j=1,natm
               if (i.eq.j) goto 100
               if (gtype_num(i).eq.0) goto 100
               ri=ra2_val(ipt,i)
               rj=ra2_val(ipt,j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(ipt,j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               if (tij.ge.-a.and.tij.le.a) then
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
               endif
 100        continue ! j=1,i-1
  95     continue ! i=1,natm
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,natm,1.0d0)
         do 195 i=1,natm
            if (gtype_num(i).eq.0) goto 195
            do 200 j=1,natm
               if (i.eq.j) goto 200
               if (gtype_num(j).eq.0) goto 200
               tij=uij(i,j)
c
               if (tij.ge.-a.and.tij.le.a) then
               vij=tij/a
               vij2=vij*vij
               vij3=vij2*vij
               vij5=vij3*vij2
               vij7=vij5*vij2

               fij=0.0625d0*(35.0d0*(vij-vij3)+21.0d0*vij5-5.0d0*vij7)
               sij=0.5d0*(1.0d0-fij)
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0/a
               dfijdv=0.0625d0*35.0d0*(1.0d0-vij2)**3
               dsijdu=-0.5d0*dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=0.0d0
               d2fijdv2=0.0625d0*35.0d0*3*(-2)*vij*(1.0d0-vij2)**2
               d2sijdu2=-0.5d0*d2fijdv2*sij*dvijdu*dvijdu
     &                  -0.5d0*dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
               else
                  if (tij.gt.a) pk(i)=0.0d0
                  do k=1,3
                     dsijdi(k,i,j)=0.0d0
                     dsijdj(k,i,j)=0.0d0
                     dsijdk(k,i,j)=0.0d0
                  enddo
                  do k=1,3
                     do l=1,3
                        d2sijdi2(k,l,i,j)=0.0d0
                        d2sijdj2(k,l,i,j)=0.0d0
                        d2sijdk2(k,l,i,j)=0.0d0
                        d2sijdidj(k,l,i,j)=0.0d0
                        d2sijdidk(k,l,i,j)=0.0d0
                        d2sijdjdk(k,l,i,j)=0.0d0
                     enddo
                  enddo
               endif
c
 200        continue ! j=1,i-1
 195     continue ! i=1,natm
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,natm
            if (gtype_num(i).ne.0) Z=Z+pk(i)
         enddo ! i=1,natm
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*natm,0.0d0)
         do 295 i=1,natm
            if (gtype_num(i).eq.0) goto 295
            do 300 j=1,natm
               if (i.eq.j) goto 300
               if (gtype_num(j).eq.0) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,latm)=gZ(k,latm)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,i-1
 295     continue !i=1,natm
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*natm*3*natm,0.0d0)
         lx=3*(latm-1)
         do 520 i=1,natm
            if (gtype_num(i).eq.0) goto 520
            ix=3*(i-1)
            do 500 j=1,natm
               if (j.eq.i) goto 500
               if (gtype_num(j).eq.0) goto 500
               jx=3*(j-1)
               do 510 k=1,natm
                  if (k.eq.i.or.k.eq.j) goto 510
                  if (gtype_num(k).eq.0) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,natm
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(latm)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,natm
            if (i.eq.latm.or.gtype_num(i).eq.0) goto 400
            do k=1,3
               gwt(k,ipt,i)=dsijdj(k,latm,i)-gZ(k,i)
               gwt(k,ipt,latm)=gwt(k,ipt,latm)
     &                        +dsijdi(k,latm,i)+dsijdk(k,latm,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,ipt,latm)=gwt(k,ipt,latm)-gZ(k,latm)
         enddo
         do i=1,natm
            do k=1,3
               gwt(k,ipt,i)=wt(ipt)*gwt(k,ipt,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         call aclear_dp(g2wt,3*natm*3*natm,0.0d0)
         lx=3*(latm-1) 
         do 600 i=1,natm 
            if (i.eq.latm.or.gtype_num(i).eq.0) goto 600
            ix=3*(i-1)
            do 610 j=1,natm
               if (j.eq.i.or.j.eq.latm.or.gtype_num(j).eq.0) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &               *(dsijdi(jc,latm,j)+dsijdk(jc,latm,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,latm,i)*dsijdj(jc,latm,j)
     &               -dsijdj(ic,latm,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,latm,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *dsijdj(jc,latm,i)
     &               -(dsijdi(ic,latm,j)+dsijdk(ic,latm,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,latm)
     &            -gZ(ic,latm)
     &            *(dsijdi(jc,latm,i)+dsijdk(jc,latm,i))
     &            +(d2sijdi2(ic,jc,latm,i)
     &             +d2sijdk2(ic,jc,latm,i)
     &             +d2sijdidk(ic,jc,latm,i)
     &             +d2sijdidk(jc,ic,latm,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,latm,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,latm,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,latm,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,latm,i)+dsijdk(ic,latm,i))
     &            *gZ(jc,i)
     &            -gZ(ic,latm)*dsijdj(jc,latm,i)
     &            +(d2sijdidj(ic,jc,latm,i)+d2sijdjdk(jc,ic,latm,i))
     &            +2*gZ(ic,latm)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,latm)*gZ(jc,latm)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*natm
            do j=1,3*natm
               hess(j,i)=hess(j,i)+g2wt(j,i)*wt(ipt)*xc_ept(ipt)
            enddo
         enddo
      enddo ! ipt
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfwt_scr(ra2_val,latm,arad,wt,  
     &     sk, awt, totwt, npts, mxp, natm,
     &     near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff )

      implicit none

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
c
c     in variables
c
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      REAL arad(natm)

      integer latm
      integer num_near_atoms, near_atom_list(*)
      REAL rshell, rnear
c
c     in/out variables
c
      REAL wt(npts)
c
c     Work space
c
      REAL awt(npts), totwt(npts), sk(npts)
      integer ibuff(npts)
c
c     Local variables
c
      integer i,j,ipt, iat, jat, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, uoa3, uoa5, uoa7, aij
c
c     End declarations
c     ***********************************************************
c
      a = 0.60d0
      ainv = 1.0d0/a
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
      if(rshell .lt. 0.5d0*(1.0d0-a)*rnear)return
c
c     first compute cell factor for local atom
c
      call aclear_dp(sk,npts,1.0d0)
      iat = latm
      igrid = gtype_num(iat)
      do 10 j=1,num_near_atoms

         jat = near_atom_list(j)
         if (jat.eq.iat) goto 10
         jgrid = gtype_num(jat)

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               if(uij .lt. -a)then
c                 sk values unchanged
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
               else
                  uoa = uij * ainv
                  uoa2 = uoa*uoa
                  uoa3 = uoa*uoa2
                  uoa5 = uoa3*uoa2
                  uoa7 = uoa5*uoa2
                  fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                  fk = fk * 0.0625d0
                  sk(ipt)=sk(ipt)*0.5d0*(1.0d0-fk)
               endif
            endif
         enddo
 10   continue

      call dcopy(npts,sk,1,awt,1)
      call dcopy(npts,sk,1,totwt,1)
c
c     loop over neighbours
c
      do 20 i=1,num_near_atoms
         iat = near_atom_list(i)
         jat = latm
         if (iat.eq.latm) goto 20
         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)
         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               if(uij .lt. -a)then
                  sk(ipt) = 1.0d0
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
                  ibuff(ipt)=0
               else
                  uoa = uij * ainv
                  uoa2 = uoa*uoa
                  uoa3 = uoa*uoa2
                  uoa5 = uoa3*uoa2
                  uoa7 = uoa5*uoa2
                  fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                  fk = fk * 0.0625d0
                  sk(ipt)=0.5d0*(1.0d0-fk)
               endif
            endif
         enddo

         do 30 j=1,num_near_atoms
            jat = near_atom_list(j)
            if(iat.ne.jat.and.jat.ne.latm)then
               jgrid = gtype_num(jat)
            
               diji = 1.0d0/dij(iat,jat)

               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        if(uij .lt. -a)then
c                          sk values unchanged
                        else if(uij .gt. a)then
                           sk(ipt) = 0.0d0
                           ibuff(ipt) = 0
                        else
                           uoa = uij * ainv
                           uoa2 = uoa*uoa
                           uoa3 = uoa*uoa2
                           uoa5 = uoa3*uoa2
                           uoa7 = uoa5*uoa2
                           fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -
     &                          5.0d0*uoa7
                           fk = fk * 0.0625d0
                           sk(ipt)=sk(ipt)*0.5d0*(1.0d0-fk)
                        endif
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               totwt(ipt) = totwt(ipt) + sk(ipt)
            endif
         enddo
c
 20   continue
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfgwt_scr(ra2_val,ra2_comp,latm,arad,wt,gwt,
     &     gwt_avail_sw, pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp,
     &     natm, near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff )
c
c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates the Becke weight [1] and the gradient of the weight   *
c *   according to Johnson et al. [2] but using the Stratmann-s et al. *
c *   [3] function for s(u).                                           *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. Vol. 98 *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Achieving linear scaling in exchange-correlation density    *
c *        functional quadratures", R.E. Stratmann, G.E. Scuseria,     *
c *        M.J. Frisch, Chem.Phys.Lett. Vol. 257 (1996) 213            *
c *                                                                    *
c **********************************************************************
      implicit none
c **********************************************************************
c *                                                                    *
c *   Declarations                                                     *
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
INCLUDE(common/dft_parameters)
c *                                                                    *
c *   In variables                                                     *
c *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      REAL ra2_comp(mxp,natoms,3)
      REAL arad(natm)
      integer latm
      integer num_near_atoms, near_atom_list(*)
      REAL rshell, rnear
c *                                                                    *
c *   In/out variables                                                 *
c *                                                                    *
      REAL wt(npts)
      REAL gwt(3,mxp,num_near_atoms)
      logical gwt_avail_sw
c *                                                                    *
c *   Work space                                                       *
c *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,num_near_atoms),gZ(3,npts,num_near_atoms)
      REAL ti_ij(3,npts,num_near_atoms)
      REAL tj_ij(3,npts,num_near_atoms)
      integer ibuff(npts)
c *                                                                    *
c *   Local variables                                                  *
c *                                                                    *
      integer i,j,ipt, iat, jat, igrid, jgrid, ii
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, uoa3, uoa5, uoa7, aij
      REAL vij, dvij, dfk, dsk, sk
      REAL dui_ijx, dui_ijy, dui_ijz, duj_ijx, duj_ijy, duj_ijz
c *                                                                    *
c *   End declarations                                                 *
c **********************************************************************
c
      a = 0.60d0
      ainv = 1.0d0/a
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
      if(rshell .lt. 0.5d0*(1.0d0-a)*rnear) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
c     first compute cell factor for local atom
c
      do j=1,num_near_atoms
         if (near_atom_list(j).eq.latm) ii=j
      enddo
      call aclear_dp(gZ(1,1,ii),3*npts,0.0d0)
      call aclear_dp(pa,npts,1.0d0)
c
      iat = near_atom_list(ii)
      igrid = gtype_num(iat)
      do 10 j=1,num_near_atoms

         jat = near_atom_list(j)
         if (jat.eq.iat) goto 10
         jgrid = gtype_num(jat)

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij

               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pa(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  uoa3 = uoa*uoa2
                  uoa5 = uoa3*uoa2
                  uoa7 = uoa5*uoa2
                  fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                  fk = fk * 0.0625d0
                  sk = 0.5d0*(1.0d0-fk)
                  pa(ipt)=pa(ipt)*sk
c
                  dfk = dvij*35.0d0*(1.0d0-uoa2)**3
                  dfk = dfk * 0.0625d0 * ainv
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk =-0.5d0*dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo
 10   continue

      call dcopy(npts,pa,1,Z,1)
      i = ii
      do 15 j=1,num_near_atoms
         if (j.eq.i) goto 15
         do ipt=1,npts
            gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
            gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
            gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
            gZ(1,ipt,j)  = gpa(1,ipt,j)
            gZ(2,ipt,j)  = gpa(2,ipt,j)
            gZ(3,ipt,j)  = gpa(3,ipt,j)
            gZ(1,ipt,i)  = gZ(1,ipt,i) + pa(ipt)*ti_ij(1,ipt,j)
            gZ(2,ipt,i)  = gZ(2,ipt,i) + pa(ipt)*ti_ij(2,ipt,j)
            gZ(3,ipt,i)  = gZ(3,ipt,i) + pa(ipt)*ti_ij(3,ipt,j)
         enddo
 15   continue
c
c     loop over neighbours
c
      do 20 i=1,num_near_atoms
         j   = ii
         iat = near_atom_list(i)
         if (iat.eq.latm) goto 20
         jat = near_atom_list(j)
         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij
c
               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged (i.e. 1)
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pk(ipt) = 1.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ibuff(ipt)=0
                  pk(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  uoa3 = uoa*uoa2
                  uoa5 = uoa3*uoa2
                  uoa7 = uoa5*uoa2
                  fk = 35.0d0*(uoa - uoa3) + 21.0d0*uoa5 -5.0d0*uoa7
                  fk = fk * 0.0625d0
                  sk = 0.5d0*(1.0d0-fk)
                  pk(ipt)=sk
c
                  dfk = dvij*35.0d0*(1.0d0-uoa2)**3
                  dfk = dfk * 0.0625d0 * ainv
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk =-0.5d0*dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo

         do 30 j=1,num_near_atoms
            jat = near_atom_list(j)
            if(iat.ne.jat.and.jat.ne.latm)then
               jgrid = gtype_num(jat)
            
               diji = 1.0d0/dij(iat,jat)

               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        vij=uij
c
                        if(vij .lt. -a)then
c
c...                       sk =1.0 therefore pk remains unchanged
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else if(vij .gt. a)then
c
c...                       sk =0.0 therefore pk drops to zero
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ibuff(ipt) = 0
                           pk(ipt) = 0.0d0
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else
c
                           dui_ijx=-ra2_comp(ipt,iat,1)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           dui_ijy=-ra2_comp(ipt,iat,2)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           dui_ijz=-ra2_comp(ipt,iat,3)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
                           duj_ijx= ra2_comp(ipt,jat,1)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           duj_ijy= ra2_comp(ipt,jat,2)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           duj_ijz= ra2_comp(ipt,jat,3)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
c
                           dvij=1.0d0
c
                           uoa = vij * ainv
                           uoa2 = uoa*uoa
                           uoa3 = uoa*uoa2
                           uoa5 = uoa3*uoa2
                           uoa7 = uoa5*uoa2
                           fk = 1.09375d0*(uoa - uoa3) + 0.65625d0*uoa5
     +                        - 0.15625d0*uoa7
                           sk=(0.5d0-fk)
                           pk(ipt)=pk(ipt)*sk
c
                           dfk = dvij*35.0d0*(1.0d0-uoa2)**3
                           dfk = dfk * 0.0625d0 * ainv
                           if (sk.lt.1.0d-15) then
                              dsk = 0.0d0
                           else
                              dsk =-0.5d0*dfk/sk
                           endif
                           ti_ij(1,ipt,j)=dsk*dui_ijx
                           ti_ij(2,ipt,j)=dsk*dui_ijy
                           ti_ij(3,ipt,j)=dsk*dui_ijz
                           tj_ij(1,ipt,j)=dsk*duj_ijx
                           tj_ij(2,ipt,j)=dsk*duj_ijy
                           tj_ij(3,ipt,j)=dsk*duj_ijz
                        endif
                     else
                        ti_ij(1,ipt,j) = 0.0d0
                        ti_ij(2,ipt,j) = 0.0d0
                        ti_ij(3,ipt,j) = 0.0d0
                        tj_ij(1,ipt,j) = 0.0d0
                        tj_ij(2,ipt,j) = 0.0d0
                        tj_ij(3,ipt,j) = 0.0d0
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               Z(ipt) = Z(ipt) + pk(ipt)
            endif
         enddo
c
         do j=1,num_near_atoms
            if (j.ne.i) then
               do ipt=1,npts
                  if (ra2_val(ipt,iat).le.arad(iat)) then
                     gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                     gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                     gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                     gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                     gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                     gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
                  endif
               enddo
            endif
         enddo
c
 20   continue
c
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      i = ii
      call aclear_dp(gwt(1,1,i),3*mxp,0.0d0)
      do 25 j=1,num_near_atoms
         if (i.eq.j) goto 25
         do ipt=1,npts
            gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(gpa(1,ipt,j)
     +                      -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
            gwt(1,ipt,i) = gwt(1,ipt,i) - gwt(1,ipt,j)
            gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(gpa(2,ipt,j)
     +                      -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
            gwt(2,ipt,i) = gwt(2,ipt,i) - gwt(2,ipt,j)
            gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(gpa(3,ipt,j)
     +                      -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
            gwt(3,ipt,i) = gwt(3,ipt,i) - gwt(3,ipt,j)
         enddo
 25   continue
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ssfg2wt_scr(ra2_val_g,ra2_comp_g,ra2_val,ra2_comp,
     &     latm,num_near_g,near_atom_g,near_atoms,indx_atoms,arad,
     &     wt,gwt_g,gwt_avail_sw, 
     &     xij_g, rij_g, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, gwt,
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, rshell, rnear)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, latm, nprt
      REAL ra2_val_g(mxp,natoms)
      REAL ra2_comp_g(mxp,natoms,3)
      REAL xc_ept(mxp)
      REAL arad(natm)
      REAL rshell, rnear
      integer num_near_g
      integer near_atom_g(num_near_g)
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt_g(3,mxp,num_near_g)
      logical gwt_avail_sw
c
c...  Workspace
c
      REAL xij_g(3,num_near_g,num_near_g)
      REAL rij_g(2,num_near_g,num_near_g) ! rij(i,*) = dij**(-i)
      REAL xij(3,num_near_g,num_near_g)
      REAL rij(2,num_near_g,num_near_g) ! rij(i,*) stores dij**(-i)
      REAL uij(num_near_g,num_near_g)
      REAL pk(num_near_g)
      REAL ra2_val(num_near_g)
      REAL ra2_comp(num_near_g,3)
      REAL duijdi(3,num_near_g,num_near_g)
      REAL duijdj(3,num_near_g,num_near_g)
      REAL duijdk(3,num_near_g,num_near_g)
      REAL d2uijdi2(3,3,num_near_g,num_near_g)
      REAL d2uijdj2(3,3,num_near_g,num_near_g)
      REAL d2uijdk2(3,3,num_near_g,num_near_g)
      REAL d2uijdidj(3,3,num_near_g,num_near_g)
      REAL d2uijdidk(3,3,num_near_g,num_near_g)
      REAL d2uijdjdk(3,3,num_near_g,num_near_g)
      REAL gZ(3,num_near_g)
      REAL g2Z(3*num_near_g,3*num_near_g)
      REAL gwt(3,num_near_g)
      REAL g2wt(3*num_near_g,3*num_near_g)
      REAL dsijdi(3,num_near_g,num_near_g)
      REAL dsijdj(3,num_near_g,num_near_g)
      REAL dsijdk(3,num_near_g,num_near_g)
      REAL d2sijdi2(3,3,num_near_g,num_near_g)
      REAL d2sijdj2(3,3,num_near_g,num_near_g)
      REAL d2sijdk2(3,3,num_near_g,num_near_g)
      REAL d2sijdidj(3,3,num_near_g,num_near_g)
      REAL d2sijdidk(3,3,num_near_g,num_near_g)
      REAL d2sijdjdk(3,3,num_near_g,num_near_g)
      integer near_atoms(num_near_g) ! per pt the real near atoms
      integer indx_atoms(num_near_g) ! per pt the index in near_atom_g
c
c...  Local variables
c
      integer indx_latm
      integer num_near
      integer i, j, k, l, n, ipt
      integer ix, jx, kx, lx
      integer iy, jy
      integer iatm, jatm, ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, vij2, vij3, vij5, vij7
      REAL fij, dfijdv, d2fijdv2
      REAL dvijdu, d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
c
c...  Parameters
c
      REAL a
      parameter(a=0.60d0)
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      if (num_near_g.eq.0) then
         gwt_avail_sw = .false.
         return
      endif
      if (rshell.lt.0.5d0*(1.0d0-a)*rnear) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      do 5 i=1,num_near_g
         iatm = near_atom_g(i)
         if (gtype_num(iatm).eq.0) goto 5
         do 10 j=1,num_near_g
            if (i.eq.j) goto 10
            jatm = near_atom_g(j)
            if (gtype_num(jatm).eq.0) goto 10
c
            xij_g(1,i,j)=atom_c(iatm,1)-atom_c(jatm,1)
            xij_g(2,i,j)=atom_c(iatm,2)-atom_c(jatm,2)
            xij_g(3,i,j)=atom_c(iatm,3)-atom_c(jatm,3)
c
            tij=1.0d0/dij(iatm,jatm)
            rij_g(1,i,j)=tij
            rij_g(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt_g,3*num_near_g*mxp,0.0d0)
      do 20 ipt=1,npts
c
c...     Construct the near atom list specific to the current point
c
         indx_latm = 0
         num_near = 0
         do i=1,num_near_g
            iatm = near_atom_g(i)
            if (ra2_val_g(ipt,iatm).le.arad(iatm).and.
     &          gtype_num(iatm).ne.0) then
               if (iatm.eq.latm) then
                  indx_latm = i
               else
                  num_near = num_near+1
                  near_atoms(num_near) = iatm
                  indx_atoms(num_near) = i
               endif
            endif
         enddo
         if (indx_latm.ne.0) then
            num_near = num_near+1
            near_atoms(num_near) = latm
            indx_atoms(num_near) = indx_latm
         endif
         if (num_near.le.1) goto 20
c
         n = num_near_g
         call aclear_dp(gwt,3*n,0.0d0)
         call aclear_dp(g2wt,3*n*3*n,0.0d0)
         call aclear_dp(duijdi,3*n*n,0.0d0)
         call aclear_dp(duijdk,3*n*n,0.0d0)
         call aclear_dp(d2uijdi2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdk2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdidk,3*3*n*n,0.0d0)
c
c...     Gather the data required for this point
c
         do i=1,num_near
            iatm=near_atoms(i)
            do k=1,3
               ra2_comp(i,k)=ra2_comp_g(ipt,iatm,k)
            enddo
            ra2_val(i)=ra2_val_g(ipt,iatm)
         enddo
         do i=1,num_near
            iatm=indx_atoms(i)
            do j=1,num_near
               jatm=indx_atoms(j)
               do k=1,3
                  xij(k,j,i)=xij_g(k,jatm,iatm)
               enddo
               do k=1,2
                  rij(k,j,i)=rij_g(k,jatm,iatm)
               enddo
            enddo
         enddo
c
c...     Construct uij and its derivatives
c
         do i=1,num_near
            do k=1,3
               rix(k)=ra2_comp(i,k)
            enddo
            do 100 j=1,num_near
               if (i.eq.j) goto 100
               ri=ra2_val(i)
               rj=ra2_val(j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
c
               if (tij.ge.-a.and.tij.le.a) then
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
               endif
 100        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,num_near,1.0d0)
         do i=1,num_near
            do 200 j=1,num_near
               if (i.eq.j) goto 200
               tij=uij(i,j)
c
               if (tij.ge.-a.and.tij.le.a) then
               vij=tij/a
               vij2=vij*vij
               vij3=vij2*vij
               vij5=vij3*vij2
               vij7=vij5*vij2

               fij=0.0625d0*(35.0d0*(vij-vij3)+21.0d0*vij5-5.0d0*vij7)
               sij=0.5d0*(1.0d0-fij)
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=1.0d0/a
               dfijdv=0.0625d0*35.0d0*(1.0d0-vij2)**3
               dsijdu=-0.5d0*dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=0.0d0
               d2fijdv2=0.0625d0*35.0d0*3*(-2)*vij*(1.0d0-vij2)**2
               d2sijdu2=-0.5d0*d2fijdv2*sij*dvijdu*dvijdu
     &                  -0.5d0*dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
               else
                  if (tij.gt.a) pk(i)=0.0d0
                  do k=1,3
                     dsijdi(k,i,j)=0.0d0
                     dsijdj(k,i,j)=0.0d0
                     dsijdk(k,i,j)=0.0d0
                  enddo
                  do k=1,3
                     do l=1,3
                        d2sijdi2(k,l,i,j)=0.0d0
                        d2sijdj2(k,l,i,j)=0.0d0
                        d2sijdk2(k,l,i,j)=0.0d0
                        d2sijdidj(k,l,i,j)=0.0d0
                        d2sijdidk(k,l,i,j)=0.0d0
                        d2sijdjdk(k,l,i,j)=0.0d0
                     enddo
                  enddo
               endif
c
 200        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,num_near
            Z=Z+pk(i)
         enddo ! i=1,num_near
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*num_near_g,0.0d0)
         do i=1,num_near
            do 300 j=1,num_near
               if (i.eq.j) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,num_near)=gZ(k,num_near)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*num_near_g*3*num_near_g,0.0d0)
         lx=3*(num_near-1)
         do 520 i=1,num_near
            ix=3*(i-1)
            do 500 j=1,num_near
               if (j.eq.i) goto 500
               jx=3*(j-1)
               do 510 k=1,num_near
                  if (k.eq.i.or.k.eq.j) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,num_near
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(num_near)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,num_near-1
            do k=1,3
               gwt(k,i)=dsijdj(k,num_near,i)-gZ(k,i)
               gwt(k,num_near)=gwt(k,num_near)
     &                        +dsijdi(k,num_near,i)+dsijdk(k,num_near,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,num_near)=gwt(k,num_near)-gZ(k,num_near)
         enddo
         do i=1,num_near
            do k=1,3
               gwt(k,i)=wt(ipt)*gwt(k,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         lx=3*(num_near-1) 
         do 600 i=1,num_near-1
            ix=3*(i-1)
            do 610 j=1,num_near-1
               if (j.eq.i) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &               *(dsijdi(jc,num_near,j)+dsijdk(jc,num_near,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,num_near,i)*dsijdj(jc,num_near,j)
     &               -dsijdj(ic,num_near,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,num_near,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,num_near)
     &            -gZ(ic,num_near)
     &            *(dsijdi(jc,num_near,i)+dsijdk(jc,num_near,i))
     &            +(d2sijdi2(ic,jc,num_near,i)
     &             +d2sijdk2(ic,jc,num_near,i)
     &             +d2sijdidk(ic,jc,num_near,i)
     &             +d2sijdidk(jc,ic,num_near,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,num_near,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,num_near,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,num_near,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,num_near)*gZ(jc,num_near)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2wt(j,i)=g2wt(j,i)*wt(ipt)
            enddo
         enddo
c
c...     Scatter the results into their final location
c
         do i=1,num_near
            iatm=indx_atoms(i)
            do k=1,3
               gwt_g(k,ipt,iatm)=gwt(k,i)
            enddo
         enddo
         do i=1,num_near
            ix=3*(near_atom_g(indx_atoms(i))-1)
            iy=3*(i-1)
            do ic=1,3
               do j=1,num_near
                  jx=3*(near_atom_g(indx_atoms(j))-1)
                  jy=3*(j-1)
                  do jc=1,3
                     hess(jx+jc,ix+ic)=hess(jx+jc,ix+ic)
     &               +g2wt(jy+jc,iy+ic)*xc_ept(ipt)
                  enddo ! jc
               enddo ! j
            enddo ! ic
         enddo ! i
c        
 20   continue ! ipt
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl4ssfwt_scr(ra2_val,latm,arad,wt,  
     &     sk, awt, totwt, npts, mxp, natm, 
     &     near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff, igrd )

      implicit none
c
c     This subroutine uses the screening technique developed by 
c
c         Stratmann, Scuseria, Frisch, CPL, 257 (1996) p213
c
c     but instead of their cell function it uses the cell function
c     suggested by
c
c         Murray, Handy, Laming, Mol.Phys. 78 (1993) p997
c
c     with m = 4 and the scaling factor a chosen such that
c     ds(q;a,m=4)/dq = ds(q;m=10)/dq.
c
c     So the cell function should be approximately as steep as the one
c     used with the MHL weighting scheme, have 4 derivatives equal 0
c     at q = -a or q = a, and has a value for a < 1 so that the 
c     screening technique by SSF can be applied.
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
c
c     in
c
      integer npts, mxp, natm
      REAL ra2_val(mxp,natoms)
      REAL arad(natm)
      integer latm, igrd
      integer num_near_atoms, near_atom_list(*)
      REAL rshell, rnear
c
c     in/out variables
c
      REAL wt(npts)
c
c     work space
c
      REAL awt(npts), totwt(npts), sk(npts)
      integer ibuff(npts)
c
c     local variables
c
      integer i,j,ipt, iat, jat, m, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, aij, b(0:4)
c
c     End declarations
c     ***********************************************************
c
      data b/
     +    -0.123046875d+01,
     +     0.164062500d+01,
     +    -0.147656250d+01,
     +     0.703125000d+00,
     +    -0.136718750d+00/
      a = 0.6651d0
      ainv = 1.0d0/a
c
      if (num_near_atoms.eq.0) return
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
c     if(rshell .lt. 0.5d0*(1.0d0-a)*rnear)return
c
c     first compute cell factor for local atom
c
      call aclear_dp(sk,npts,1.0d0)
      iat = latm
      igrid = gtype_num(iat)

      do 10 j=1,num_near_atoms

         jat = near_atom_list(j)
         if (jat.eq.iat) goto 10
         jgrid = gtype_num(jat)
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               uij=uij+aij*(1.0d0-uij*uij)
               if(uij .lt. -a)then
c                 sk values unchanged
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
               else
                  uoa  = uij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m = 1, 4
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk(ipt)=sk(ipt)*fk
               endif
            endif
         enddo
 10   continue

      call dcopy(npts,sk,1,awt,1)
      call dcopy(npts,sk,1,totwt,1)
c
c     loop over neighbours
c
      do 20 i=1,num_near_atoms
         iat = near_atom_list(i)
         if (iat.eq.latm) goto 20
         jat = latm

         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)

         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               uij=uij+aij*(1.0d0-uij*uij)
               if(uij .lt. -a)then
                  sk(ipt) = 1.0d0
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
                  ibuff(ipt)=0
               else
                  uoa  = uij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m = 1, 4
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk(ipt)=fk
               endif
            endif
         enddo

         do 30 j=1,num_near_atoms
            jat = near_atom_list(j)
            if(i.ne.j.and.jat.ne.latm)then
               jgrid = gtype_num(jat)
               radi=weight_atom_radius(igrid,igrd)
               radj=weight_atom_radius(jgrid,igrd)
               ratij=radi/radj
               uij=(ratij-1.0d0)/(ratij+1.0d0)
               aij=uij/(uij*uij-1.0d0)
               if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij
               
               diji = 1.0d0/dij(iat,jat)

               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        uij=uij+aij*(1.0d0-uij*uij)
                        if(uij .lt. -a)then
c                          sk values unchanged
                        else if(uij .gt. a)then
                           sk(ipt) = 0.0d0
                           ibuff(ipt) = 0
                        else
                           uoa  = uij * ainv
                           uoa2 = uoa*uoa
                           fk = 0.5d0 + b(0)*uoa
                           do m = 1, 4
                              uoa = uoa*uoa2
                              fk = fk + b(m)*uoa
                           enddo
                           sk(ipt)=sk(ipt)*fk
                        endif
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               totwt(ipt) = totwt(ipt) + sk(ipt)
            endif
         enddo
 20   continue
c     
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl4ssfgwt_scr(ra2_val,ra2_comp,latm,arad,wt,gwt,
     &     gwt_avail_sw,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm,
     &     near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff, igrd )
c
c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates the Becke weight [1] and the gradient of the weight   *
c *   according to Johnson et al. [2] but using the Stratmann-s et al. *
c *   [3] approach to save computation and using Murray-s et al. [4]   *
c *   function for s(u), with m=4 and the scaling factor a chosen such *
c *   that ds(u;a,m=4)/du = ds(u;m=10)/du.                             *
c *                                                                    *
c *   This scheme leads to a cut-off profile that is approximately as  *
c *   steep as the one used with the MHL weighting scheme [4], has 4   *
c *   derivatives equal 0 at u = -a and u = a, and has values 1 if     *
c *   u < -a and 0 if u > a so that the SSF screening technique can be *
c *   applied.                                                         *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. Vol. 98 *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Achieving linear scaling in exchange-correlation density    *
c *        functional quadratures", R.E. Stratmann, G.E. Scuseria,     *
c *        M.J. Frisch, Chem.Phys.Lett. Vol. 257 (1996) 213            *
c *                                                                    *
c *   [4] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming, Mol.Phys.    *
c *        Vol. 78 (1993) 997                                          *
c *                                                                    *
c **********************************************************************
      implicit none
c **********************************************************************
c *                                                                    *
c *   Declarations                                                     *
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
INCLUDE(common/dft_parameters)
      integer mm
      parameter (mm=4)
c *                                                                    *
c *   In variables                                                     *
c *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natoms)
      REAL ra2_comp(mxp,natoms,3)
      REAL arad(natm)
      integer latm, igrd
      integer num_near_atoms, near_atom_list(*)
      REAL rshell, rnear
c *                                                                    *
c *   In/out variables                                                 *
c *                                                                    *
      REAL wt(npts)
      REAL gwt(3,mxp,num_near_atoms)
      logical gwt_avail_sw
c *                                                                    *
c *   Work space                                                       *
c *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,num_near_atoms),gZ(3,npts,num_near_atoms)
      REAL ti_ij(3,npts,num_near_atoms)
      REAL tj_ij(3,npts,num_near_atoms)
      integer ibuff(npts)
c *                                                                    *
c *   Local variables                                                  *
c *                                                                    *
      integer i, j, ii, ipt, iat, jat, m, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, uoa3, uoa5, uoa7, aij
      REAL vij, dvij, dfk, dsk, sk
      REAL dui_ijx, dui_ijy, dui_ijz, duj_ijx, duj_ijy, duj_ijz
      REAL b(0:mm)
c *                                                                    *
c *   End declarations                                                 *
c **********************************************************************
c
      data b/
     +    -0.123046875d+01,
     +     0.164062500d+01,
     +    -0.147656250d+01,
     +     0.703125000d+00,
     +    -0.136718750d+00/
      a = 0.6651d0
      ainv = 1.0d0/a
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
c     if(rshell .lt. 0.5d0*(1.0d0-a)*rnear) then
c        gwt_avail_sw = .false.
c        return
c     endif
      if(num_near_atoms.eq.0) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
c     first compute cell factor for local atom
c
      do j=1,num_near_atoms
         if (near_atom_list(j).eq.latm) ii=j
      enddo
      i=num_near_atoms
      call aclear_dp(gZ(1,1,ii),3*npts,0.0d0)
      call aclear_dp(pa,npts,1.0d0)
c
      iat = near_atom_list(ii)
      igrid = gtype_num(iat)
      do 10 j=1,num_near_atoms

         jat = near_atom_list(j)
         if (jat.eq.iat) goto 10
         jgrid = gtype_num(jat)

         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij+aij*(1.0d0-uij*uij)

               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pa(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0-2.0d0*aij*uij
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m=1,mm
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk = fk
                  pa(ipt)=pa(ipt)*sk
c
                  dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk = dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo
 10   continue

      i = ii
      call dcopy(npts,pa,1,Z,1)
      do 15 j=1,num_near_atoms
         if (j.eq.i) goto 15
         do ipt=1,npts
            gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
            gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
            gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
            gZ(1,ipt,j)  = gpa(1,ipt,j)
            gZ(2,ipt,j)  = gpa(2,ipt,j)
            gZ(3,ipt,j)  = gpa(3,ipt,j)
            gZ(1,ipt,i)  = gZ(1,ipt,i) + pa(ipt)*ti_ij(1,ipt,j)
            gZ(2,ipt,i)  = gZ(2,ipt,i) + pa(ipt)*ti_ij(2,ipt,j)
            gZ(3,ipt,i)  = gZ(3,ipt,i) + pa(ipt)*ti_ij(3,ipt,j)
         enddo
 15   continue
c
c loop over neighbours
c
      do 20 i=1,num_near_atoms
         j   = ii
         if (i.eq.j) goto 20
         iat = near_atom_list(i)
         jat = near_atom_list(j)
         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij+aij*(1.0d0-uij*uij)
c
               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged (i.e. 1)
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pk(ipt) = 1.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ibuff(ipt)=0
                  pk(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0-2.0d0*aij*uij
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m=1,mm
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk = fk
                  pk(ipt)=sk
c
                  dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk = dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo

         do 30 j=1,num_near_atoms
            if(i.ne.j.and.j.ne.ii)then
               jat = near_atom_list(j)
               jgrid = gtype_num(jat)
               radi=weight_atom_radius(igrid,igrd)
               radj=weight_atom_radius(jgrid,igrd)
               ratij=radi/radj
               uij=(ratij-1.0d0)/(ratij+1.0d0)
               aij=uij/(uij*uij-1.0d0)
               if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij
               
               diji = 1.0d0/dij(iat,jat)
   
               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
c
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        vij=uij+aij*(1.0d0-uij*uij)
c
                        if(vij .lt. -a)then
c
c...                       sk =1.0 therefore pk remains unchanged
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else if(vij .gt. a)then
c
c...                       sk =0.0 therefore pk drops to zero
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ibuff(ipt) = 0
                           pk(ipt) = 0.0d0
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else
c
                           dui_ijx=-ra2_comp(ipt,iat,1)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           dui_ijy=-ra2_comp(ipt,iat,2)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           dui_ijz=-ra2_comp(ipt,iat,3)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
                           duj_ijx= ra2_comp(ipt,jat,1)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           duj_ijy= ra2_comp(ipt,jat,2)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           duj_ijz= ra2_comp(ipt,jat,3)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
c
                           dvij=1.0d0-2.0d0*aij*uij
c
                           uoa = vij * ainv
                           uoa2 = uoa*uoa
                           fk = 0.5d0 + b(0)*uoa
                           do m=1,mm
                              uoa = uoa*uoa2
                              fk = fk + b(m)*uoa
                           enddo
                           sk = fk
                           pk(ipt)=pk(ipt)*sk
c
                           dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                           if (sk.lt.1.0d-15) then
                              dsk = 0.0d0
                           else
                              dsk = dfk/sk
                           endif
                           ti_ij(1,ipt,j)=dsk*dui_ijx
                           ti_ij(2,ipt,j)=dsk*dui_ijy
                           ti_ij(3,ipt,j)=dsk*dui_ijz
                           tj_ij(1,ipt,j)=dsk*duj_ijx
                           tj_ij(2,ipt,j)=dsk*duj_ijy
                           tj_ij(3,ipt,j)=dsk*duj_ijz
                        endif
                     else
                        ti_ij(1,ipt,j) = 0.0d0
                        ti_ij(2,ipt,j) = 0.0d0
                        ti_ij(3,ipt,j) = 0.0d0
                        tj_ij(1,ipt,j) = 0.0d0
                        tj_ij(2,ipt,j) = 0.0d0
                        tj_ij(3,ipt,j) = 0.0d0
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               Z(ipt) = Z(ipt) + pk(ipt)
            endif
         enddo
c
         do j=1,num_near_atoms
            if (j.ne.i) then
               do ipt=1,npts
                  if (ra2_val(ipt,iat).le.arad(iat)) then
                     gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                     gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                     gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                     gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                     gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                     gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
                  endif
               enddo
            endif
         enddo
 20   continue
c
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      i = ii
      call aclear_dp(gwt(1,1,i),3*mxp,0.0d0)
      do 25 j=1,num_near_atoms
         if (j.eq.ii) goto 25
         do ipt=1,npts
            gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(gpa(1,ipt,j)
     +                      -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
            gwt(1,ipt,i) = gwt(1,ipt,i) - gwt(1,ipt,j)
            gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(gpa(2,ipt,j)
     +                      -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
            gwt(2,ipt,i) = gwt(2,ipt,i) - gwt(2,ipt,j)
            gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(gpa(3,ipt,j)
     +                      -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
            gwt(3,ipt,i) = gwt(3,ipt,i) - gwt(3,ipt,j)
         enddo
 25   continue
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl4ssfg2wt_scr(ra2_val_g,ra2_comp_g,ra2_val,ra2_comp,
     &     latm,num_near_g,near_atom_g,near_atoms,indx_atoms,arad,
     &     wt,gwt_g,gwt_avail_sw, 
     &     aij_g, xij_g, rij_g, aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, gwt,
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, latm, igrd, nprt
      REAL ra2_val_g(mxp,natoms)
      REAL ra2_comp_g(mxp,natoms,3)
      REAL xc_ept(mxp)
      REAL arad(natm)
      integer num_near_g
      integer near_atom_g(num_near_g)
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt_g(3,mxp,num_near_g)
      logical gwt_avail_sw
c
c...  Workspace
c
      REAL aij_g(num_near_g,num_near_g)
      REAL xij_g(3,num_near_g,num_near_g)
      REAL rij_g(2,num_near_g,num_near_g) ! rij(i,*) = dij**(-i)
      REAL aij(num_near_g,num_near_g)
      REAL xij(3,num_near_g,num_near_g)
      REAL rij(2,num_near_g,num_near_g) ! rij(i,*) stores dij**(-i)
      REAL uij(num_near_g,num_near_g)
      REAL pk(num_near_g)
      REAL ra2_val(num_near_g)
      REAL ra2_comp(num_near_g,3)
      REAL duijdi(3,num_near_g,num_near_g)
      REAL duijdj(3,num_near_g,num_near_g)
      REAL duijdk(3,num_near_g,num_near_g)
      REAL d2uijdi2(3,3,num_near_g,num_near_g)
      REAL d2uijdj2(3,3,num_near_g,num_near_g)
      REAL d2uijdk2(3,3,num_near_g,num_near_g)
      REAL d2uijdidj(3,3,num_near_g,num_near_g)
      REAL d2uijdidk(3,3,num_near_g,num_near_g)
      REAL d2uijdjdk(3,3,num_near_g,num_near_g)
      REAL gZ(3,num_near_g)
      REAL g2Z(3*num_near_g,3*num_near_g)
      REAL gwt(3,num_near_g)
      REAL g2wt(3*num_near_g,3*num_near_g)
      REAL dsijdi(3,num_near_g,num_near_g)
      REAL dsijdj(3,num_near_g,num_near_g)
      REAL dsijdk(3,num_near_g,num_near_g)
      REAL d2sijdi2(3,3,num_near_g,num_near_g)
      REAL d2sijdj2(3,3,num_near_g,num_near_g)
      REAL d2sijdk2(3,3,num_near_g,num_near_g)
      REAL d2sijdidj(3,3,num_near_g,num_near_g)
      REAL d2sijdidk(3,3,num_near_g,num_near_g)
      REAL d2sijdjdk(3,3,num_near_g,num_near_g)
      integer near_atoms(num_near_g) ! per pt the real near atoms
      integer indx_atoms(num_near_g) ! per pt the index in near_atom_g
c
c...  Local variables
c
      integer indx_latm
      integer num_near
      integer i, j, k, l, n, ipt
      integer ix, jx, kx, lx
      integer iy, jy
      integer iatm, jatm, ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, vij2, vijm
      REAL fij, dfijdv, d2fijdv2
      REAL dvijdu, d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
      integer m
c
c...  Functions
c
c     REAL srad
c
c...  Parameters
c
      REAL a
      parameter(a=0.6651d0)
      integer mmax
      parameter (mmax=4)
      REAL b(0:mmax)
      data b/
     +    -0.123046875d+01,
     +     0.164062500d+01,
     +    -0.147656250d+01,
     +     0.703125000d+00,
     +    -0.136718750d+00/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      if (num_near_g.eq.1) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      do 5 i=1,num_near_g
         iatm = near_atom_g(i)
         igrid = gtype_num(iatm)
         if (igrid.eq.0) goto 5
         do 10 j=1,num_near_g
            if (i.eq.j) goto 10
            jatm = near_atom_g(j)
            jgrid = gtype_num(jatm)
            if (jgrid.eq.0) goto 10
c
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij_g(i,j)=tij
c
            xij_g(1,i,j)=atom_c(iatm,1)-atom_c(jatm,1)
            xij_g(2,i,j)=atom_c(iatm,2)-atom_c(jatm,2)
            xij_g(3,i,j)=atom_c(iatm,3)-atom_c(jatm,3)
c
            tij=1.0d0/dij(iatm,jatm)
            rij_g(1,i,j)=tij
            rij_g(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt_g,3*num_near_g*mxp,0.0d0)
      do 20 ipt=1,npts
c
c...     Construct the near atom list specific to the current point
c
         num_near = 0
         indx_latm = 0
         do i=1,num_near_g
            iatm = near_atom_g(i)
            if (ra2_val_g(ipt,iatm).le.arad(iatm).and.
     &          gtype_num(iatm).ne.0) then
               if (iatm.eq.latm) then
                  indx_latm = i
               else
                  num_near = num_near+1
                  near_atoms(num_near) = iatm
                  indx_atoms(num_near) = i
               endif
            endif
         enddo
         if (indx_latm.ne.0) then
            num_near = num_near+1
            near_atoms(num_near) = latm
            indx_atoms(num_near) = indx_latm
         endif
         if (num_near.le.1) goto 20
c
         n = num_near_g
         call aclear_dp(gwt,3*n,0.0d0)
         call aclear_dp(g2wt,3*n*3*n,0.0d0)
         call aclear_dp(duijdi,3*n*n,0.0d0)
         call aclear_dp(duijdk,3*n*n,0.0d0)
         call aclear_dp(d2uijdi2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdk2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdidk,3*3*n*n,0.0d0)
c
c...     Gather the data required for this point
c
         do i=1,num_near
            iatm=near_atoms(i)
            do k=1,3
               ra2_comp(i,k)=ra2_comp_g(ipt,iatm,k)
            enddo
            ra2_val(i)=ra2_val_g(ipt,iatm)
         enddo
         do i=1,num_near
            iatm=indx_atoms(i)
            do j=1,num_near
               jatm=indx_atoms(j)
               aij(j,i)=aij_g(jatm,iatm)
               do k=1,3
                  xij(k,j,i)=xij_g(k,jatm,iatm)
               enddo
               do k=1,2
                  rij(k,j,i)=rij_g(k,jatm,iatm)
               enddo
            enddo
         enddo
c
c...     Construct uij and its derivatives
c
         do i=1,num_near
            do k=1,3
               rix(k)=ra2_comp(i,k)
            enddo
            do 100 j=1,num_near
               if (i.eq.j) goto 100
               ri=ra2_val(i)
               rj=ra2_val(j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
c
               if (vij.ge.-a.and.vij.le.a) then
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
               endif
 100        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,num_near,1.0d0)
         do i=1,num_near
            do 200 j=1,num_near
               if (i.eq.j) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
c
               if (vij.ge.-a.and.vij.le.a) then
               vij=vij/a
               vij2=vij*vij
               vijm=vij
               sij=0.5d0 + b(0)*vijm
               do m = 1, mmax
                  vijm = vijm * vij2
                  sij  = sij  + b(m)*vijm
               enddo
c
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=(1.0d0-2.0d0*aij(i,j)*tij)/a
               dfijdv=b(0)*(1.0d0-vij2)**mmax
               dsijdu=dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)/a
               d2fijdv2=-2.0d0*mmax*b(0)*vij*(1.0d0-vij2)**(mmax-1)
               d2sijdu2=d2fijdv2*sij*dvijdu*dvijdu
     &                 +dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
               else
                  if (vij.gt.a) pk(i)=0.0d0
                  do k=1,3
                     dsijdi(k,i,j)=0.0d0
                     dsijdj(k,i,j)=0.0d0
                     dsijdk(k,i,j)=0.0d0
                  enddo
                  do k=1,3
                     do l=1,3
                        d2sijdi2(k,l,i,j)=0.0d0
                        d2sijdj2(k,l,i,j)=0.0d0
                        d2sijdk2(k,l,i,j)=0.0d0
                        d2sijdidj(k,l,i,j)=0.0d0
                        d2sijdidk(k,l,i,j)=0.0d0
                        d2sijdjdk(k,l,i,j)=0.0d0
                     enddo
                  enddo
               endif
c
 200        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,num_near
            Z=Z+pk(i)
         enddo ! i=1,num_near
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*num_near_g,0.0d0)
         do i=1,num_near
            do 300 j=1,num_near
               if (i.eq.j) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,num_near)=gZ(k,num_near)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*num_near_g*3*num_near_g,0.0d0)
         lx=3*(num_near-1)
         do 520 i=1,num_near
            ix=3*(i-1)
            do 500 j=1,num_near
               if (j.eq.i) goto 500
               jx=3*(j-1)
               do 510 k=1,num_near
                  if (k.eq.i.or.k.eq.j) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,num_near
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(num_near)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,num_near-1
            do k=1,3
               gwt(k,i)=dsijdj(k,num_near,i)-gZ(k,i)
               gwt(k,num_near)=gwt(k,num_near)
     &                        +dsijdi(k,num_near,i)+dsijdk(k,num_near,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,num_near)=gwt(k,num_near)-gZ(k,num_near)
         enddo
         do i=1,num_near
            do k=1,3
               gwt(k,i)=wt(ipt)*gwt(k,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         lx=3*(num_near-1) 
         do 600 i=1,num_near-1
            ix=3*(i-1)
            do 610 j=1,num_near-1
               if (j.eq.i) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &               *(dsijdi(jc,num_near,j)+dsijdk(jc,num_near,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,num_near,i)*dsijdj(jc,num_near,j)
     &               -dsijdj(ic,num_near,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,num_near,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,num_near)
     &            -gZ(ic,num_near)
     &            *(dsijdi(jc,num_near,i)+dsijdk(jc,num_near,i))
     &            +(d2sijdi2(ic,jc,num_near,i)
     &             +d2sijdk2(ic,jc,num_near,i)
     &             +d2sijdidk(ic,jc,num_near,i)
     &             +d2sijdidk(jc,ic,num_near,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,num_near,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,num_near,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,num_near,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,num_near)*gZ(jc,num_near)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2wt(j,i)=g2wt(j,i)*wt(ipt)
            enddo
         enddo
c
c...     Scatter the results into their final location
c
         do i=1,num_near
            iatm=indx_atoms(i)
            do k=1,3
               gwt_g(k,ipt,iatm)=gwt(k,i)
            enddo
         enddo
         do i=1,num_near
            ix=3*(near_atom_g(indx_atoms(i))-1)
            iy=3*(i-1)
            do ic=1,3
               do j=1,num_near
                  jx=3*(near_atom_g(indx_atoms(j))-1)
                  jy=3*(j-1)
                  do jc=1,3
                     hess(jx+jc,ix+ic)=hess(jx+jc,ix+ic)
     &               +g2wt(jy+jc,iy+ic)*xc_ept(ipt)
                  enddo ! jc
               enddo ! j
            enddo ! ic
         enddo ! i
c        
 20   continue ! ipt
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl8ssfwt_scr(ra2_val,latm,arad,wt,  
     &     sk, awt, totwt, npts, mxp, natm,
     &     near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff, igrd )

      implicit none
c
c     This subroutine uses the screening technique developed by 
c
c         Stratmann, Scuseria, Frisch, CPL, 257 (1996) p213
c
c     but instead of their cell function it uses the cell function
c     suggested by
c
c         Murray, Handy, Laming, Mol.Phys. 78 (1993) p997
c
c     with m = 8 and the scaling factor a chosen such that
c     ds(q;a,m=8)/dq = ds(q;m=10)/dq.
c
c     So the cell function should be approximately as steep as the one
c     used with the MHL weighting scheme, have 8 derivatives equal 0
c     at q = -a or q = a, and has a value for a < 1 so that the 
c     screening technique by SSF can be applied.
c
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
c
c     in
c
      integer npts, mxp, natm
      REAL ra2_val(mxp,natoms)
      REAL arad(natm)
      integer latm
      integer igrd
      integer num_near_atoms, near_atom_list(num_near_atoms)
      REAL rshell, rnear
c
c     in/out variables
c
      REAL wt(npts)
c
c     work space
c
      REAL awt(npts), totwt(npts), sk(npts)
      integer ibuff(npts)
c
c     local variables
c
      integer i,j,ipt, iat, jat, m, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, aij, b(0:8)
c
c     end declarations
c     ***********************************************************
c
      data b/
     +    -0.16692352294921837d+01,
     +     0.44512939453124893d+01,
     +    -0.93477172851562287d+01,
     +     0.13353881835937466d+02,
     +    -0.12982940673828091d+02,
     +     0.84979248046874787d+01,
     +    -0.35952758789062407d+01,
     +     0.89025878906249778d+00,
     +    -0.98190307617187250d-01/
      a = 0.902256d0
      ainv = 1.0d0/a
c
      if (num_near_atoms.eq.0) return
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
c     if(rshell .lt. 0.5d0*(1.0d0-a)*rnear)return
c
c     first compute cell factor for local atom
c
      call aclear_dp(sk,npts,1.0d0)
      iat = latm
      igrid = gtype_num(iat)
      do 10 j=1,num_near_atoms

         jat = near_atom_list(j)
         if (jat.eq.iat) goto 10
         jgrid = gtype_num(jat)
         if (jgrid.eq.0) goto 10
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               uij=uij+aij*(1.0d0-uij*uij)
               if(uij .lt. -a)then
c                 sk values unchanged
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
               else
                  uoa  = uij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m = 1, 8
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk(ipt)=sk(ipt)*fk
               endif
            endif
         enddo
 10   continue

      call dcopy(npts,sk,1,awt,1)
      call dcopy(npts,sk,1,totwt,1)
c
c loop over neighbours
c
      do 20 i=1,num_near_atoms
         iat = near_atom_list(i)
         jat = latm
         if (iat.eq.jat) goto 20
         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               uij=uij+aij*(1.0d0-uij*uij)
               if(uij .lt. -a)then
                  sk(ipt) = 1.0d0
               else if(uij .gt. a)then
                  sk(ipt) = 0.0d0
                  ibuff(ipt)=0
               else
                  uoa  = uij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m = 1, 8
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk(ipt)=fk
               endif
            endif
         enddo

         do 30 j=1,num_near_atoms
            jat = near_atom_list(j)
            if(i.ne.j.and.jat.ne.latm)then
               jgrid = gtype_num(jat)
               radi=weight_atom_radius(igrid,igrd)
               radj=weight_atom_radius(jgrid,igrd)
               ratij=radi/radj
               uij=(ratij-1.0d0)/(ratij+1.0d0)
               aij=uij/(uij*uij-1.0d0)
               if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij
            
               diji = 1.0d0/dij(iat,jat)

               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        uij=uij+aij*(1.0d0-uij*uij)
                        if(uij .lt. -a)then
c                          sk values unchanged
                        else if(uij .gt. a)then
                           sk(ipt) = 0.0d0
                           ibuff(ipt) = 0
                        else
                           uoa  = uij * ainv
                           uoa2 = uoa*uoa
                           fk = 0.5d0 + b(0)*uoa
                           do m = 1, 8
                              uoa = uoa*uoa2
                              fk = fk + b(m)*uoa
                           enddo
                           sk(ipt)=sk(ipt)*fk
                        endif
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               totwt(ipt) = totwt(ipt) + sk(ipt)
            endif
         enddo
 20   continue
c     
      do ipt=1,npts
         wt(ipt)=wt(ipt)*awt(ipt)/totwt(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl8ssfgwt_scr(ra2_val,ra2_comp,latm,arad,wt,gwt,
     &     gwt_avail_sw,
     &     pk, pa, Z, gpa, gZ, ti_ij, tj_ij, npts, mxp, natm,
     &     near_atom_list, num_near_atoms,
     &     rshell, rnear, ibuff, igrd )
c
c **********************************************************************
c *   Huub van Dam, 1999                                               *
c *                                                                    *
c *   Description:                                                     *
c *                                                                    *
c *   Calculates the Becke weight [1] and the gradient of the weight   *
c *   according to Johnson et al. [2] but using the Stratmann-s et al. *
c *   [3] approach to save computation and using Murray-s et al. [4]   *
c *   function for s(u), with m=8 and the scaling factor a chosen such *
c *   that ds(u;a,m=8)/du = ds(u;m=10)/du.                             *
c *                                                                    *
c *   This scheme leads to a cut-off profile that is approximately as  *
c *   steep as the one used with the MHL weighting scheme [4], has 4   *
c *   derivatives equal 0 at u = -a and u = a, and has values 1 if     *
c *   u < -a and 0 if u > a so that the SSF screening technique can be *
c *   applied.                                                         *
c *                                                                    *
c *   [1] "A multicenter numerical integration scheme for polyatomic   *
c *        molecules", A.D. Becke, J.Chem.Phys. Vol. 88 (1988) 2547    *
c *                                                                    *
c *   [2] "The performance of a family of density functional methods"  *
c *        B.G. Johnson, P.M.W. Gill, J.A. Pople, J.Chem.Phys. Vol. 98 *
c *        (1993) 5612                                                 *
c *                                                                    *
c *   [3] "Achieving linear scaling in exchange-correlation density    *
c *        functional quadratures", R.E. Stratmann, G.E. Scuseria,     *
c *        M.J. Frisch, Chem.Phys.Lett. Vol. 257 (1996) 213            *
c *                                                                    *
c *   [4] "Quadrature schemes for integrals of density functional      *
c *        theory", C.W. Murray, N.C. Handy, G.J. Laming, Mol.Phys.    *
c *        Vol. 78 (1993) 997                                          *
c *                                                                    *
c **********************************************************************
      implicit none
c **********************************************************************
c *                                                                    *
c *   Declarations                                                     *
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
INCLUDE(common/dft_parameters)
      integer mm
      parameter (mm=8)
c *                                                                    *
c *   In variables                                                     *
c *                                                                    *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_dij)
      integer npts, mxp, natm
      REAL ra2_val(mxp,natm)
      REAL ra2_comp(mxp,natoms,3)
      REAL arad(natm)
      integer latm
      integer igrd
      integer num_near_atoms
      integer near_atom_list(num_near_atoms)
      REAL rshell, rnear
c *                                                                    *
c *   In/out variables                                                 *
c *                                                                    *
      REAL wt(npts)
      REAL gwt(3,mxp,num_near_atoms)
      logical gwt_avail_sw
c *                                                                    *
c *   Work space                                                       *
c *                                                                    *
      REAL pa(npts), Z(npts), pk(npts)
      REAL gpa(3,npts,num_near_atoms),gZ(3,npts,num_near_atoms)
      REAL ti_ij(3,npts,num_near_atoms)
      REAL tj_ij(3,npts,num_near_atoms)
      integer ibuff(npts)
c *                                                                    *
c *   Local variables                                                  *
c *                                                                    *
      integer i, j, ii, ipt, iat, jat, m, igrid, jgrid
      REAL radi,radj,ratij,uij,fk,diji
      REAL a, ainv, uoa, uoa2, uoa3, uoa5, uoa7, aij
      REAL vij, dvij, dfk, dsk, sk
      REAL dui_ijx, dui_ijy, dui_ijz, duj_ijx, duj_ijy, duj_ijz
      REAL b(0:mm)
c *                                                                    *
c *   End declarations                                                 *
c **********************************************************************
c
      data b/
     +    -0.16692352294921837d+01,
     +     0.44512939453124893d+01,
     +    -0.93477172851562287d+01,
     +     0.13353881835937466d+02,
     +    -0.12982940673828091d+02,
     +     0.84979248046874787d+01,
     +    -0.35952758789062407d+01,
     +     0.89025878906249778d+00,
     +    -0.98190307617187250d-01/
      a = 0.902256d0
      ainv = 1.0d0/a
c
c     if the shell is within sphere defined by eq 15. in SSF
c     there is no change to the weight for any of the points
c     in this batch
c
c     if(rshell .lt. 0.5d0*(1.0d0-a)*rnear) then
c        gwt_avail_sw = .false.
c        return
c     endif
      if(num_near_atoms.eq.1) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
c     first compute cell factor for local atom
c
      do j=1,num_near_atoms
         if (near_atom_list(j).eq.latm) ii=j
      enddo
      call aclear_dp(gZ(1,1,ii),3*npts,0.0d0)
      call aclear_dp(pa,npts,1.0d0)
c
      iat = near_atom_list(ii)
      igrid = gtype_num(iat)
      do 10 j=1,num_near_atoms

         if (j.eq.ii) goto 10
         jat = near_atom_list(j)
         jgrid = gtype_num(jat)
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij

         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts

            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij+aij*(1.0d0-uij*uij)

               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pa(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0-2.0d0*aij*uij
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m=1,mm
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk = fk
                  pa(ipt)=pa(ipt)*sk
c
                  dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk = dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo
 10   continue
c
      i = ii
      call dcopy(npts,pa,1,Z,1)
      do 15 j=1,num_near_atoms
         if (j.eq.i) goto 15
         do ipt=1,npts
            gpa(1,ipt,j) = pa(ipt)*tj_ij(1,ipt,j)
            gpa(2,ipt,j) = pa(ipt)*tj_ij(2,ipt,j)
            gpa(3,ipt,j) = pa(ipt)*tj_ij(3,ipt,j)
            gZ(1,ipt,j)  = gpa(1,ipt,j)
            gZ(2,ipt,j)  = gpa(2,ipt,j)
            gZ(3,ipt,j)  = gpa(3,ipt,j)
            gZ(1,ipt,i)  = gZ(1,ipt,i) + pa(ipt)*ti_ij(1,ipt,j)
            gZ(2,ipt,i)  = gZ(2,ipt,i) + pa(ipt)*ti_ij(2,ipt,j)
            gZ(3,ipt,i)  = gZ(3,ipt,i) + pa(ipt)*ti_ij(3,ipt,j)
         enddo
 15   continue
c
c loop over neighbours
c
      do 20 i=1,num_near_atoms
         j   = ii
         if (i.eq.j) goto 20
         iat = near_atom_list(i)
         jat = near_atom_list(j)
         igrid = gtype_num(iat)
         jgrid = gtype_num(jat)
         radi=weight_atom_radius(igrid,igrd)
         radj=weight_atom_radius(jgrid,igrd)
         ratij=radi/radj
         uij=(ratij-1.0d0)/(ratij+1.0d0)
         aij=uij/(uij*uij-1.0d0)
         if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij
         diji = 1.0d0/dij(iat,jat)

         do ipt=1,npts
            ibuff(ipt)=1
            if (ra2_val(ipt,iat).le.arad(iat).and.
     &          ra2_val(ipt,jat).le.arad(jat)) then
               uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
               vij=uij+aij*(1.0d0-uij*uij)
c
               if(vij .lt. -a)then
c
c...              sk =1.0 therefore pk remains unchanged (i.e. 1)
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  pk(ipt) = 1.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else if(vij .gt. a)then
c
c...              sk =0.0 therefore pk drops to zero
c...              dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                  ibuff(ipt)=0
                  pk(ipt) = 0.0d0
                  ti_ij(1,ipt,j) = 0.0d0
                  ti_ij(2,ipt,j) = 0.0d0
                  ti_ij(3,ipt,j) = 0.0d0
                  tj_ij(1,ipt,j) = 0.0d0
                  tj_ij(2,ipt,j) = 0.0d0
                  tj_ij(3,ipt,j) = 0.0d0
c
               else
c
                  dui_ijx=-ra2_comp(ipt,iat,1)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  dui_ijy=-ra2_comp(ipt,iat,2)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  dui_ijz=-ra2_comp(ipt,iat,3)*diji/ra2_val(ipt,iat)
     +                    -(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
                  duj_ijx= ra2_comp(ipt,jat,1)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,1)-atom_c(jat,1))*uij*diji*diji
                  duj_ijy= ra2_comp(ipt,jat,2)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,2)-atom_c(jat,2))*uij*diji*diji
                  duj_ijz= ra2_comp(ipt,jat,3)*diji/ra2_val(ipt,jat)
     +                    +(atom_c(iat,3)-atom_c(jat,3))*uij*diji*diji
c
                  dvij=1.0d0-2.0d0*aij*uij
c
                  uoa = vij * ainv
                  uoa2 = uoa*uoa
                  fk = 0.5d0 + b(0)*uoa
                  do m=1,mm
                     uoa = uoa*uoa2
                     fk = fk + b(m)*uoa
                  enddo
                  sk = fk
                  pk(ipt)=sk
c
                  dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                  if (sk.lt.1.0d-15) then
                     dsk = 0.0d0
                  else
                     dsk = dfk/sk
                  endif
                  ti_ij(1,ipt,j)=dsk*dui_ijx
                  ti_ij(2,ipt,j)=dsk*dui_ijy
                  ti_ij(3,ipt,j)=dsk*dui_ijz
                  tj_ij(1,ipt,j)=dsk*duj_ijx
                  tj_ij(2,ipt,j)=dsk*duj_ijy
                  tj_ij(3,ipt,j)=dsk*duj_ijz
               endif
            else
               ti_ij(1,ipt,j) = 0.0d0
               ti_ij(2,ipt,j) = 0.0d0
               ti_ij(3,ipt,j) = 0.0d0
               tj_ij(1,ipt,j) = 0.0d0
               tj_ij(2,ipt,j) = 0.0d0
               tj_ij(3,ipt,j) = 0.0d0
            endif
         enddo

         do 30 j=1,num_near_atoms
            if(i.ne.j.and.j.ne.ii)then
               jat = near_atom_list(j)
               jgrid = gtype_num(jat)
               radi=weight_atom_radius(igrid,igrd)
               radj=weight_atom_radius(jgrid,igrd)
               ratij=radi/radj
               uij=(ratij-1.0d0)/(ratij+1.0d0)
               aij=uij/(uij*uij-1.0d0)
               if(abs(aij).gt.0.5d0) aij=0.5d0*abs(aij)/aij
            
               diji = 1.0d0/dij(iat,jat)

               do ipt=1,npts
                  if(ibuff(ipt) .ne. 0)then
c
                     if (ra2_val(ipt,iat).le.arad(iat).and.
     &                   ra2_val(ipt,jat).le.arad(jat)) then
                        uij=diji*(ra2_val(ipt,iat)-ra2_val(ipt,jat))
                        vij=uij+aij*(1.0d0-uij*uij)
c
                        if(vij .lt. -a)then
c
c...                       sk =1.0 therefore pk remains unchanged
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else if(vij .gt. a)then
c
c...                       sk =0.0 therefore pk drops to zero
c...                       dsk=0.0 therefore ti_ij=0.0 and tj_ij=0.0
c
                           ibuff(ipt) = 0
                           pk(ipt) = 0.0d0
                           ti_ij(1,ipt,j) = 0.0d0
                           ti_ij(2,ipt,j) = 0.0d0
                           ti_ij(3,ipt,j) = 0.0d0
                           tj_ij(1,ipt,j) = 0.0d0
                           tj_ij(2,ipt,j) = 0.0d0
                           tj_ij(3,ipt,j) = 0.0d0
c
                        else
c
                           dui_ijx=-ra2_comp(ipt,iat,1)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           dui_ijy=-ra2_comp(ipt,iat,2)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           dui_ijz=-ra2_comp(ipt,iat,3)*diji
     +                            /ra2_val(ipt,iat)
     +                            -(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
                           duj_ijx= ra2_comp(ipt,jat,1)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,1)-atom_c(jat,1))
     +                            *uij*diji*diji
                           duj_ijy= ra2_comp(ipt,jat,2)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,2)-atom_c(jat,2))
     +                            *uij*diji*diji
                           duj_ijz= ra2_comp(ipt,jat,3)*diji
     +                            /ra2_val(ipt,jat)
     +                            +(atom_c(iat,3)-atom_c(jat,3))
     +                            *uij*diji*diji
c
                           dvij=1.0d0-2.0d0*aij*uij
c
                           uoa = vij * ainv
                           uoa2 = uoa*uoa
                           fk = 0.5d0 + b(0)*uoa
                           do m=1,mm
                              uoa = uoa*uoa2
                              fk = fk + b(m)*uoa
                           enddo
                           sk = fk
                           pk(ipt)=pk(ipt)*sk
c
                           dfk = ainv*dvij*b(0)*(1.0d0-uoa2)**mm
                           if (sk.lt.1.0d-15) then
                              dsk = 0.0d0
                           else
                              dsk = dfk/sk
                           endif
                           ti_ij(1,ipt,j)=dsk*dui_ijx
                           ti_ij(2,ipt,j)=dsk*dui_ijy
                           ti_ij(3,ipt,j)=dsk*dui_ijz
                           tj_ij(1,ipt,j)=dsk*duj_ijx
                           tj_ij(2,ipt,j)=dsk*duj_ijy
                           tj_ij(3,ipt,j)=dsk*duj_ijz
                        endif
                     else
                        ti_ij(1,ipt,j) = 0.0d0
                        ti_ij(2,ipt,j) = 0.0d0
                        ti_ij(3,ipt,j) = 0.0d0
                        tj_ij(1,ipt,j) = 0.0d0
                        tj_ij(2,ipt,j) = 0.0d0
                        tj_ij(3,ipt,j) = 0.0d0
                     endif
                  endif
               enddo
            endif
 30      continue
c
         do ipt=1,npts
            if (ra2_val(ipt,iat).le.arad(iat)) then
               Z(ipt) = Z(ipt) + pk(ipt)
            endif
         enddo
c
         do j=1,num_near_atoms
            if (j.ne.i) then
               do ipt=1,npts
                  if (ra2_val(ipt,iat).le.arad(iat)) then
                     gZ(1,ipt,j) = gZ(1,ipt,j) + pk(ipt)*tj_ij(1,ipt,j)
                     gZ(2,ipt,j) = gZ(2,ipt,j) + pk(ipt)*tj_ij(2,ipt,j)
                     gZ(3,ipt,j) = gZ(3,ipt,j) + pk(ipt)*tj_ij(3,ipt,j)
                     gZ(1,ipt,i) = gZ(1,ipt,i) + pk(ipt)*ti_ij(1,ipt,j)
                     gZ(2,ipt,i) = gZ(2,ipt,i) + pk(ipt)*ti_ij(2,ipt,j)
                     gZ(3,ipt,i) = gZ(3,ipt,i) + pk(ipt)*ti_ij(3,ipt,j)
                  endif
               enddo
            endif
         enddo
c
 20   continue
c
      do ipt=1,npts
         if (Z(ipt).lt.1.0d-15) then
            Z(ipt)=0.0d0
         else
            Z(ipt)=1.0d0/Z(ipt)
         endif
      enddo
c
      i = ii
      call aclear_dp(gwt(1,1,i),3*mxp,0.0d0)
      do 25 j=1,num_near_atoms
         if (j.eq.i) goto 25
         do ipt=1,npts
            gwt(1,ipt,j) = wt(ipt)*Z(ipt)*(gpa(1,ipt,j)
     +                      -pa(ipt)*gZ(1,ipt,j)*Z(ipt))
            gwt(1,ipt,i) = gwt(1,ipt,i) - gwt(1,ipt,j)
            gwt(2,ipt,j) = wt(ipt)*Z(ipt)*(gpa(2,ipt,j)
     +                      -pa(ipt)*gZ(2,ipt,j)*Z(ipt))
            gwt(2,ipt,i) = gwt(2,ipt,i) - gwt(2,ipt,j)
            gwt(3,ipt,j) = wt(ipt)*Z(ipt)*(gpa(3,ipt,j)
     +                      -pa(ipt)*gZ(3,ipt,j)*Z(ipt))
            gwt(3,ipt,i) = gwt(3,ipt,i) - gwt(3,ipt,j)
         enddo
 25   continue
c
      do ipt=1,npts
         wt(ipt)=wt(ipt)*pa(ipt)*Z(ipt)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine mhl8ssfg2wt_scr(ra2_val_g,ra2_comp_g,ra2_val,ra2_comp,
     &     latm,num_near_g,near_atom_g,near_atoms,indx_atoms,arad,
     &     wt,gwt_g,gwt_avail_sw, 
     &     aij_g, xij_g, rij_g, aij, uij, xij, rij, pk,
     &     duijdi, duijdj, duijdk, dsijdi, dsijdj, dsijdk, gZ, gwt,
     &     d2uijdi2,  d2uijdj2,  d2uijdk2,
     &     d2uijdidj, d2uijdidk, d2uijdjdk,
     &     d2sijdi2,  d2sijdj2,  d2sijdk2,
     &     d2sijdidj, d2sijdidk, d2sijdjdk,
     &     g2Z, g2wt,
     &     xc_ept, hess, nprt,
     &     npts, mxp, natm, igrd)

C  *********************************************************************
C  *Description:                                                       *
C  *Calculates Becke weights at the grid points and also calculates    *
C  *gradients of the weights and the hessian of the weights            *
C  *********************************************************************
c
c     Due to total panic trying Tozer approach now
c
      implicit none
c
c...  Parameters
c
INCLUDE(common/dft_parameters)
INCLUDE(../m4/common/sizes)
c
c...  Inputs
c
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_dij)
INCLUDE(common/dft_basder)
INCLUDE(../m4/common/mapper)
c
      integer npts, mxp, natm, latm, igrd, nprt
      REAL ra2_val_g(mxp,natoms)
      REAL ra2_comp_g(mxp,natoms,3)
      REAL xc_ept(mxp)
      REAL arad(natm)
      integer num_near_g
      integer near_atom_g(num_near_g)
c
c...  Input/Outputs
c
      REAL wt(npts)
      REAL hess(nprt,nprt)
c
c...  Outputs
c
      REAL gwt_g(3,mxp,num_near_g)
      logical gwt_avail_sw
c
c...  Workspace
c
      REAL aij_g(num_near_g,num_near_g)
      REAL xij_g(3,num_near_g,num_near_g)
      REAL rij_g(2,num_near_g,num_near_g) ! rij(i,*) = dij**(-i)
      REAL aij(num_near_g,num_near_g)
      REAL xij(3,num_near_g,num_near_g)
      REAL rij(2,num_near_g,num_near_g) ! rij(i,*) stores dij**(-i)
      REAL uij(num_near_g,num_near_g)
      REAL pk(num_near_g)
      REAL ra2_val(num_near_g)
      REAL ra2_comp(num_near_g,3)
      REAL duijdi(3,num_near_g,num_near_g)
      REAL duijdj(3,num_near_g,num_near_g)
      REAL duijdk(3,num_near_g,num_near_g)
      REAL d2uijdi2(3,3,num_near_g,num_near_g)
      REAL d2uijdj2(3,3,num_near_g,num_near_g)
      REAL d2uijdk2(3,3,num_near_g,num_near_g)
      REAL d2uijdidj(3,3,num_near_g,num_near_g)
      REAL d2uijdidk(3,3,num_near_g,num_near_g)
      REAL d2uijdjdk(3,3,num_near_g,num_near_g)
      REAL gZ(3,num_near_g)
      REAL g2Z(3*num_near_g,3*num_near_g)
      REAL gwt(3,num_near_g)
      REAL g2wt(3*num_near_g,3*num_near_g)
      REAL dsijdi(3,num_near_g,num_near_g)
      REAL dsijdj(3,num_near_g,num_near_g)
      REAL dsijdk(3,num_near_g,num_near_g)
      REAL d2sijdi2(3,3,num_near_g,num_near_g)
      REAL d2sijdj2(3,3,num_near_g,num_near_g)
      REAL d2sijdk2(3,3,num_near_g,num_near_g)
      REAL d2sijdidj(3,3,num_near_g,num_near_g)
      REAL d2sijdidk(3,3,num_near_g,num_near_g)
      REAL d2sijdjdk(3,3,num_near_g,num_near_g)
      integer near_atoms(num_near_g) ! per pt the real near atoms
      integer indx_atoms(num_near_g) ! per pt the index in near_atom_g
c
c...  Local variables
c
      integer indx_latm
      integer num_near
      integer i, j, k, l, n, ipt
      integer ix, jx, kx, lx
      integer iy, jy
      integer iatm, jatm, ic, jc, kc
      integer igrid, jgrid
      REAL radi, radj
      REAL vij, vij2, vijm
      REAL fij, dfijdv, d2fijdv2
      REAL dvijdu, d2vijdu2
      REAL sij
      REAL dsijdu
      REAL d2sijdu2
      REAL tij
      REAL tijk, tik, tjk
      REAL tijl, til, tjl
      REAL rix(3), rjx(3), ri, rj, ri1, rj1
      REAL rijx(3), rij1, rij2
      REAL Z
      integer m
c
c...  Functions
c
c     REAL srad
c
c...  Parameters
c
      REAL a
      parameter(a=0.902256d0)
      integer mmax
      parameter (mmax=8)
      REAL b(0:mmax)
      data b/
     +    -0.16692352294921837d+01,
     +     0.44512939453124893d+01,
     +    -0.93477172851562287d+01,
     +     0.13353881835937466d+02,
     +    -0.12982940673828091d+02,
     +     0.84979248046874787d+01,
     +    -0.35952758789062407d+01,
     +     0.89025878906249778d+00,
     +    -0.98190307617187250d-01/
c
c     ==================================================================
c
c     First compute tables per atoms. This includes aij and xij and dij.
c     Note that aij = -aji and xij = -xji and dij = dji
c
      if (num_near_g.eq.1) then
         gwt_avail_sw = .false.
         return
      endif
      gwt_avail_sw = .true.
c
      do 5 i=1,num_near_g
         iatm = near_atom_g(i)
         igrid = gtype_num(iatm)
         if (igrid.eq.0) goto 5
         do 10 j=1,num_near_g
            if (i.eq.j) goto 10
            jatm = near_atom_g(j)
            jgrid = gtype_num(jatm)
            if (jgrid.eq.0) goto 10
c
            radi=weight_atom_radius(igrid,igrd)
            radj=weight_atom_radius(jgrid,igrd)
            tij=radi/radj
            tij=(tij-1.0d0)/(tij+1.0d0)
            tij=tij/((tij*tij)-1.0d0)
            if(abs(tij).gt.0.5d0)tij=0.5d0*abs(tij)/tij
            aij_g(i,j)=tij
c
            xij_g(1,i,j)=atom_c(iatm,1)-atom_c(jatm,1)
            xij_g(2,i,j)=atom_c(iatm,2)-atom_c(jatm,2)
            xij_g(3,i,j)=atom_c(iatm,3)-atom_c(jatm,3)
c
            tij=1.0d0/dij(iatm,jatm)
            rij_g(1,i,j)=tij
            rij_g(2,i,j)=tij*tij
 10      continue
  5   continue
c
c...  Compute the weights and the 1st and 2nd derivatives of the
c...  weights per grid point.
c
      call aclear_dp(gwt_g,3*num_near_g*mxp,0.0d0)
      do 20 ipt=1,npts
c
c...     Construct the near atom list specific to the current point
c
         num_near = 0
         indx_latm = 0
         do i=1,num_near_g
            iatm = near_atom_g(i)
            if (ra2_val_g(ipt,iatm).le.arad(iatm).and.
     &          gtype_num(iatm).ne.0) then
               if (iatm.eq.latm) then
                  indx_latm = i
               else
                  num_near = num_near+1
                  near_atoms(num_near) = iatm
                  indx_atoms(num_near) = i
               endif
            endif
         enddo
         if (indx_latm.ne.0) then
            num_near = num_near+1
            near_atoms(num_near) = latm
            indx_atoms(num_near) = indx_latm
         endif
         if (num_near.le.1) goto 20
c
         n = num_near_g
         call aclear_dp(gwt,3*n,0.0d0)
         call aclear_dp(g2wt,3*n*3*n,0.0d0)
         call aclear_dp(duijdi,3*n*n,0.0d0)
         call aclear_dp(duijdk,3*n*n,0.0d0)
         call aclear_dp(d2uijdi2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdk2,3*3*n*n,0.0d0)
         call aclear_dp(d2uijdidk,3*3*n*n,0.0d0)
c
c...     Gather the data required for this point
c
         do i=1,num_near
            iatm=near_atoms(i)
            do k=1,3
               ra2_comp(i,k)=ra2_comp_g(ipt,iatm,k)
            enddo
            ra2_val(i)=ra2_val_g(ipt,iatm)
         enddo
         do i=1,num_near
            iatm=indx_atoms(i)
            do j=1,num_near
               jatm=indx_atoms(j)
               aij(j,i)=aij_g(jatm,iatm)
               do k=1,3
                  xij(k,j,i)=xij_g(k,jatm,iatm)
               enddo
               do k=1,2
                  rij(k,j,i)=rij_g(k,jatm,iatm)
               enddo
            enddo
         enddo
c
c...     Construct uij and its derivatives
c
         do i=1,num_near
            do k=1,3
               rix(k)=ra2_comp(i,k)
            enddo
            do 100 j=1,num_near
               if (i.eq.j) goto 100
               ri=ra2_val(i)
               rj=ra2_val(j)
               ri1=1.0d0/ri
               rj1=1.0d0/rj
               rij1=rij(1,i,j)
               rij2=rij(2,i,j)
c
               do k=1,3
                  rijx(k)=xij(k,i,j)
                  rjx(k)=ra2_comp(j,k)
               enddo
c
               tij=(ri-rj)*rij1
               uij(i,j)=tij
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
c
               if (vij.ge.-a.and.vij.le.a) then
c
               do k=1,3
                  duijdi(k,i,j)=-rix(k)*ri1*rij1
     &                          -rijx(k)*tij*rij2
                  duijdj(k,i,j)= rjx(k)*rj1*rij1
     &                          +rijx(k)*tij*rij2
                  duijdk(k,i,j)= rix(k)*ri1*rij1
     &                          -rjx(k)*rj1*rij1
               enddo
c
               do k=1,3
                  tijk=rijx(k)*rij2
                  tik=rix(k)*ri1
                  tjk=rjx(k)*rj1
                  do 110 l=1,3
                     if (l.eq.k) goto 110
                     tijl=rijx(l)*rij2
                     til=rix(l)*ri1
                     tjl=rjx(l)*rj1
                     d2uijdi2(k,l,i,j)=
     &                   rij1*(-tik*til*ri1+til*tijk+tik*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdj2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1+tjl*tijk+tjk*tijl)
     &                  +3*tijk*tijl*tij
                     d2uijdk2(k,l,i,j)=
     &                   rij1*(tjk*tjl*rj1-tik*til*ri1)
                     d2uijdidj(k,l,i,j)=
     &                   rij1*(-tik*tijl-tjl*tijk)
     &                  -3*tijk*tijl*tij
                     d2uijdidk(k,l,i,j)=
     &                   rij1*(ri1*tik*til-til*tijk+tjl*tijk)
                     d2uijdjdk(k,l,i,j)=
     &                   rij1*(-rj1*tjk*tjl+til*tijk-tjl*tijk)
 110              continue
                  d2uijdi2(k,k,i,j)=
     &                rij1*(ri1*(1.0d0-tik*tik)
     &                      +2*tik*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdj2(k,k,i,j)=
     &                rij1*(rj1*(-1.0d0+tjk*tjk)
     &                      +2*tjk*tijk)
     &               +tij*(3*tijk*tijk-rij2)
                  d2uijdk2(k,k,i,j)=
     &                ri1*rij1*(1.0d0-tik*tik)
     &               -rj1*rij1*(1.0d0-tjk*tjk)
                  d2uijdidj(k,k,i,j)=
     &                rij1*tijk*(-tik-tjk)
     &               +tij*(rij2-3*tijk*tijk)
                  d2uijdidk(k,k,i,j)=
     &                rij1*(ri1*(tik*tik-1.0d0)
     &                     +tijk*(tjk-tik))
                  d2uijdjdk(k,k,i,j)=
     &                rij1*(rj1*(1.0d0-tjk*tjk)
     &                     +tijk*(tik-tjk))
               enddo
               endif
 100        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct Pc, (ds(C,B))/s(C,B), (dds(C,B))/s(C,B)
c
         call aclear_dp(pk,num_near,1.0d0)
         do i=1,num_near
            do 200 j=1,num_near
               if (i.eq.j) goto 200
               tij=uij(i,j)
               vij=tij+aij(i,j)*(1.0d0-tij*tij)
c
               if (vij.ge.-a.and.vij.le.a) then
               vij=vij/a
               vij2=vij*vij
               vijm=vij
               sij=0.5d0 + b(0)*vijm
               do m = 1, mmax
                  vijm = vijm * vij2
                  sij  = sij  + b(m)*vijm
               enddo
c
               pk(i)=pk(i)*sij
               if (sij.lt.1.0d-15) then
                  sij=0.0d0
               else
                  sij=1.0d0/sij
               endif
c
               dvijdu=(1.0d0-2.0d0*aij(i,j)*tij)/a
               dfijdv=b(0)*(1.0d0-vij2)**mmax
               dsijdu=dfijdv*sij*dvijdu
               do k=1,3
                  dsijdi(k,i,j)=dsijdu*duijdi(k,i,j)
                  dsijdj(k,i,j)=dsijdu*duijdj(k,i,j)
                  dsijdk(k,i,j)=dsijdu*duijdk(k,i,j)
               enddo
c
               d2vijdu2=-2.0d0*aij(i,j)/a
               d2fijdv2=-2.0d0*mmax*b(0)*vij*(1.0d0-vij2)**(mmax-1)
               d2sijdu2=d2fijdv2*sij*dvijdu*dvijdu
     &                 +dfijdv*sij*d2vijdu2
c
               do k=1,3
                  do l=1,3
                     d2sijdj2(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdj2(k,l,i,j)
                     d2sijdi2(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdi(l,i,j)
     &                  +dsijdu*d2uijdi2(k,l,i,j)
                     d2sijdk2(k,l,i,j)
     &                  =d2sijdu2*duijdk(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdk2(k,l,i,j)
                     d2sijdidj(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdj(l,i,j)
     &                  +dsijdu*d2uijdidj(k,l,i,j)
                     d2sijdidk(k,l,i,j)
     &                  =d2sijdu2*duijdi(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdidk(k,l,i,j)
                     d2sijdjdk(k,l,i,j)
     &                  =d2sijdu2*duijdj(k,i,j)*duijdk(l,i,j)
     &                  +dsijdu*d2uijdjdk(k,l,i,j)
                  enddo
               enddo
               else
                  if (vij.gt.a) pk(i)=0.0d0
                  do k=1,3
                     dsijdi(k,i,j)=0.0d0
                     dsijdj(k,i,j)=0.0d0
                     dsijdk(k,i,j)=0.0d0
                  enddo
                  do k=1,3
                     do l=1,3
                        d2sijdi2(k,l,i,j)=0.0d0
                        d2sijdj2(k,l,i,j)=0.0d0
                        d2sijdk2(k,l,i,j)=0.0d0
                        d2sijdidj(k,l,i,j)=0.0d0
                        d2sijdidk(k,l,i,j)=0.0d0
                        d2sijdjdk(k,l,i,j)=0.0d0
                     enddo
                  enddo
               endif
c
 200        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Calculate Z
c
         Z=0.0d0
         do i=1,num_near
            Z=Z+pk(i)
         enddo ! i=1,num_near
c
c...     Construct gradients of Z
c
         call aclear_dp(gZ,3*num_near_g,0.0d0)
         do i=1,num_near
            do 300 j=1,num_near
               if (i.eq.j) goto 300
               do k=1,3
                  gZ(k,j)=gZ(k,j)+pk(i)*dsijdj(k,i,j)
                  gZ(k,i)=gZ(k,i)+pk(i)*dsijdi(k,i,j)
                  gZ(k,num_near)=gZ(k,num_near)+pk(i)*dsijdk(k,i,j)
               enddo
 300        continue ! j=1,num_near
         enddo ! i=1,num_near
c
c...     Construct hessian of Z
c
         call aclear_dp(g2Z,3*num_near_g*3*num_near_g,0.0d0)
         lx=3*(num_near-1)
         do 520 i=1,num_near
            ix=3*(i-1)
            do 500 j=1,num_near
               if (j.eq.i) goto 500
               jx=3*(j-1)
               do 510 k=1,num_near
                  if (k.eq.i.or.k.eq.j) goto 510
                  kx=3*(k-1)
                  do jc=1,3
                     do kc=1,3
                        g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(ix+jc,kx+kc)=g2Z(ix+jc,kx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                                  +pk(i)*dsijdi(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(jx+jc,ix+kc)=g2Z(jx+jc,ix+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(jx+jc,kx+kc)=g2Z(jx+jc,kx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                                  +pk(i)*dsijdj(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                        g2Z(lx+jc,ix+kc)=g2Z(lx+jc,ix+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdi(kc,i,k)
                        g2Z(lx+jc,kx+kc)=g2Z(lx+jc,kx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdj(kc,i,k)
                        g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                                  +pk(i)*dsijdk(jc,i,j)
     &                                        *dsijdk(kc,i,k)
                     enddo
                  enddo
 510           continue ! k
               kx=jx
               do jc=1,3
                  do kc=1,3
                    g2Z(ix+jc,ix+kc)=g2Z(ix+jc,ix+kc)
     &                              +pk(i)*d2sijdi2(jc,kc,i,j)
                    g2Z(jx+jc,jx+kc)=g2Z(jx+jc,jx+kc)
     &                              +pk(i)*d2sijdj2(jc,kc,i,j)
                    g2Z(lx+jc,lx+kc)=g2Z(lx+jc,lx+kc)
     &                              +pk(i)*d2sijdk2(jc,kc,i,j)
                    g2Z(ix+jc,lx+kc)=g2Z(ix+jc,lx+kc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(lx+kc,ix+jc)=g2Z(lx+kc,ix+jc)
     &                              +pk(i)*d2sijdidk(jc,kc,i,j)
                    g2Z(ix+jc,jx+kc)=g2Z(ix+jc,jx+kc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+kc,ix+jc)=g2Z(jx+kc,ix+jc)
     &                              +pk(i)*d2sijdidj(jc,kc,i,j)
                    g2Z(jx+jc,lx+kc)=g2Z(jx+jc,lx+kc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                    g2Z(lx+kc,jx+jc)=g2Z(lx+kc,jx+jc)
     &                              +pk(i)*d2sijdjdk(jc,kc,i,j)
                  enddo
               enddo ! jc
 500         continue ! j
 520     continue
c
         if (Z.lt.1.0d-15) then
            Z=0.0d0
         else
            Z=1.0d0/Z
         endif
         do i=1,num_near
            do k=1,3
               gZ(k,i)=gZ(k,i)*Z
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2Z(j,i)=g2Z(j,i)*Z
            enddo
         enddo
c
c...     Calculate the weights
c
         wt(ipt)=wt(ipt)*pk(num_near)*Z
c
c...     Calculate the gradients of the weights
c
         do 400 i=1,num_near-1
            do k=1,3
               gwt(k,i)=dsijdj(k,num_near,i)-gZ(k,i)
               gwt(k,num_near)=gwt(k,num_near)
     &                        +dsijdi(k,num_near,i)+dsijdk(k,num_near,i)
            enddo
 400     continue
         do k=1,3
            gwt(k,num_near)=gwt(k,num_near)-gZ(k,num_near)
         enddo
         do i=1,num_near
            do k=1,3
               gwt(k,i)=wt(ipt)*gwt(k,i)
            enddo
         enddo
c
c...     Calculate the hessians of the weights
c 
         lx=3*(num_near-1) 
         do 600 i=1,num_near-1
            ix=3*(i-1)
            do 610 j=1,num_near-1
               if (j.eq.i) goto 610
               jx=3*(j-1)
               do ic=1,3
                  do jc=1,3
                     g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &               +(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &               *(dsijdi(jc,num_near,j)+dsijdk(jc,num_near,j))
                     g2wt(ix+ic,jx+jc)=g2wt(ix+ic,jx+jc)
     &               +dsijdj(ic,num_near,i)*dsijdj(jc,num_near,j)
     &               -dsijdj(ic,num_near,i)*gZ(jc,j)
     &               -gZ(ic,i)*dsijdj(jc,num_near,j)
     &               -g2Z(ix+ic,jx+jc)
     &               +2*gZ(ic,i)*gZ(jc,j)
                     g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
c
                     g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &               +(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *dsijdj(jc,num_near,i)
     &               -(dsijdi(ic,num_near,j)+dsijdk(ic,num_near,j))
     &               *gZ(jc,i)
                  enddo
               enddo
 610        continue
            do ic=1,3
               do jc=1,3
                  g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,num_near)
     &            -gZ(ic,num_near)
     &            *(dsijdi(jc,num_near,i)+dsijdk(jc,num_near,i))
     &            +(d2sijdi2(ic,jc,num_near,i)
     &             +d2sijdk2(ic,jc,num_near,i)
     &             +d2sijdidk(ic,jc,num_near,i)
     &             +d2sijdidk(jc,ic,num_near,i))
                  g2wt(ix+ic,ix+jc)=g2wt(ix+ic,ix+jc)
     &            -dsijdj(ic,num_near,i)*gZ(jc,i)
     &            -gZ(ic,i)*dsijdj(jc,num_near,i)
     &            -g2Z(ix+ic,ix+jc)
     &            +d2sijdj2(ic,jc,num_near,i)
     &            +2*gZ(ic,i)*gZ(jc,i)
                  g2wt(lx+ic,ix+jc)=g2wt(lx+ic,ix+jc)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
                  g2wt(ix+jc,lx+ic)=g2wt(ix+jc,lx+ic)
     &            -(dsijdi(ic,num_near,i)+dsijdk(ic,num_near,i))
     &            *gZ(jc,i)
     &            -gZ(ic,num_near)*dsijdj(jc,num_near,i)
     &            +(d2sijdidj(ic,jc,num_near,i)+
     &              d2sijdjdk(jc,ic,num_near,i))
     &            +2*gZ(ic,num_near)*gZ(jc,i)
     &            -g2Z(lx+ic,ix+jc)
               enddo
            enddo
 600     continue
         do ic=1,3
            do jc=1,3
               g2wt(lx+ic,lx+jc)=g2wt(lx+ic,lx+jc)
     &         +2*gZ(ic,num_near)*gZ(jc,num_near)
     &         -g2Z(lx+ic,lx+jc)
            enddo
         enddo
         do i=1,3*num_near
            do j=1,3*num_near
               g2wt(j,i)=g2wt(j,i)*wt(ipt)
            enddo
         enddo
c
c...     Scatter the results into their final location
c
         do i=1,num_near
            iatm=indx_atoms(i)
            do k=1,3
               gwt_g(k,ipt,iatm)=gwt(k,i)
            enddo
         enddo
         do i=1,num_near
            ix=3*(near_atom_g(indx_atoms(i))-1)
            iy=3*(i-1)
            do ic=1,3
               do j=1,num_near
                  jx=3*(near_atom_g(indx_atoms(j))-1)
                  jy=3*(j-1)
                  do jc=1,3
                     hess(jx+jc,ix+ic)=hess(jx+jc,ix+ic)
     &               +g2wt(jy+jc,iy+ic)*xc_ept(ipt)
                  enddo ! jc
               enddo ! j
            enddo ! ic
         enddo ! i
c        
 20   continue ! ipt
c
      end
