c
c Code to implement the setup of a frozen density
c expressed as a linear combination of gaussian functions
c
c It is assumed that GAMESS-UK is running and that the
c energy and force contributions will be computed there
c via the normal bq terms (for point charges) and the 
c additional calls to the 3-centre integral and gradient
c integral codes (implemented here).
c
c  - issues
c      - what happens when the coulomb fit code is also
c        in use ???
c            a) hold both sets of gaussians in the same list..
c            b) set up two basis sets and compute 2 independent
c               contributins to the energy/potential
c
c The following routines are provided
c
c Used gamini(CH) -> gamess(GAM) -> chmdat2(CH) 
c
c  gden_init
c
c  gden_add_gauss(ncentr, iout, debug)
c
c    Called once, to define the exponents. 
c    It is assumed that the coordinates have already been added
c    to the atom lists
c
c  Called from within GAMESS-UK to get terms in the 
c  energy/forces
c
c  gden_energy()
c    Compute the energy contributions
c
c  gden_forces()
c    Compute the contributions, to both AO and expansion
c    centres
c
c**********************************************************************
c
c  gden_init
c
c  by this point any control settings should have been made
c
c**********************************************************************

      function gden_init(memory_fp,iout,debug)

      implicit none

C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
C *Inter-module communication variables.
C *Suffixes and what they mean		
C *---------------------------	
C *_sw			-		switch		
C *_num 		-		numbers	
C *_ch			-		Input/output channels
      integer out_ch,in_ch
      common/io_channels/out_ch,in_ch
C *
C *Global switches
C *
      logical debug_sw
      common/global_switches/debug_sw
C *
C *Switches and numbers used in dft routines
C *
      logical optim_sw,triangle_sw
      common/scf_control_switch/optim_sw,triangle_sw

      logical jfit_sw,jfitg_sw,cmm_sw,dunlap_sw,potential_sw
      logical kqua_sw,kfit_sw
      logical rks_sw
      logical ludm_sw,svdm_sw
      logical jown_sw,dega_sw,kown_sw
      logical mult_sw, dft2e_sw
      common/scftype/rks_sw
      common/j_switch/jfit_sw,jfitg_sw,cmm_sw,mult_sw,
     &                dunlap_sw,potential_sw,dft2e_sw

      common/xc_switch/kqua_sw,kfit_sw,
     &     ludm_sw,svdm_sw,jown_sw,dega_sw,kown_sw
c
c The grid parameters
c
c     1) SG1 fully specifies everything
c     2) rad_grid_scheme specifies for each type which radial grid
c        is used
c       -1) if RG_MK then
c              radm_num specifies m
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -2) if RG_EML then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -3) if RG_B then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c     3) ang_grid_scheme specifies which angular grid to use
c       -1) if AG_LEB then
c              angupt_num specifies the maximum number of angular grid
c              points
c       -2) if AG_LEG then
c              thetpt_num specifies the maximum number of theta points
c              phipt_num specifies the maximum number of phi points
c     4) ang_prune_scheme specifies which scheme to use for pruning 
c        the angular grid as a function of the radius.
c       -1) if AP_MHL (no other info needed)
c       -2) if AP_RADZONE then
c              radzones_num specifies the number of radial zones
c              bnd_radzn specifies the location of zone boundaries
c              angpt_radzn_num specifies the number of angular grid 
c              points per zone.
c              
c     integer angupt_num,thetpt_num,phipt_num,radpt_num
      integer radpt_num
      integer weight_scheme, radzones_num, angpt_radzn_num
      integer thetpt_radzn_num, phipt_radzn_num
      integer ang_prune_scheme
      integer rad_grid_scheme, ang_grid_scheme
      integer gtype_num, ngtypes, gaccu_num
      integer iauto_cnt
      integer rad_scale_scheme
      integer radnpt_row
      integer angnpt_row
      integer grid_generation
      real*8 grid_scale, radm_num, bnd_radzn
      real*8 grid_atom_radius
      real*8 weight_atom_radius
      real*8 prune_atom_radius
      real*8 screen_atom_radius
c
c     subtle difference here: 
c
c     - grid_atom_radius:   used as scale factor of the radial grids.
c     - weight_atom_radius: used for atom-size-adjustments in the
c                           weighting scheme.
c     - prune_atom_radius:  used for pruning the angular grid in the
c                           MHL pruning scheme
c     - screen_atom_radius: used for screening of the radial grids.
c
      real*8 psitol, warntol
      logical conv_prune_sw, gradwght_sw, sort_points_sw
      logical ignore_accuracy_sw
      common/dft_grid_parameters/
     &     psitol(0:max_gtype,max_grids),
     &     warntol(0:max_gtype,max_grids),
     &     radm_num(0:max_gtype,max_grids),
     &     grid_scale(0:max_gtype,max_grids),
     &     grid_atom_radius(0:max_gtype,max_grids),
     &     weight_atom_radius(0:max_gtype,max_grids),
     &     prune_atom_radius(0:max_gtype,max_grids),
     &     screen_atom_radius(0:max_gtype,max_grids),
     &     bnd_radzn(maxradzn-1,0:max_gtype,max_grids),
     &     radnpt_row(7),angnpt_row(7),
     &     angpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     thetpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     phipt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     radzones_num(0:max_gtype,max_grids),
     &     ang_prune_scheme(0:max_gtype,max_grids),
     &     rad_grid_scheme(0:max_gtype,max_grids),
     &     ang_grid_scheme(0:max_gtype,max_grids),
     &     radpt_num(0:max_gtype,max_grids),
     &     gaccu_num(0:max_gtype,max_grids),
     &     gtype_num(max_atom),
     &     grid_generation,
     &     ngtypes,iauto_cnt,
     &     rad_scale_scheme,
     &     weight_scheme(max_grids),
     &     conv_prune_sw,
     &     gradwght_sw,
     &     sort_points_sw,
     &     ignore_accuracy_sw

      integer poleexp_num,over_tol,pener_tol,schwarz_tol
      real*8  tttt2
      common/pole_options/tttt2,poleexp_num,over_tol,pener_tol,
     &                    schwarz_tol

      integer    MAX_DEBUG
      parameter (MAX_DEBUG=25)

      logical print_sw(MAX_DEBUG)
      common/debugpr/print_sw
c
c debug array indices
c
      integer    DEBUG_KSMATRIX
      parameter (DEBUG_KSMATRIX = 1)
      integer    DEBUG_TR
      parameter (DEBUG_TR       = 2)
      integer    DEBUG_NORM
      parameter (DEBUG_NORM     = 3)
      integer    DEBUG_DENSITY
      parameter (DEBUG_DENSITY  = 4)
      integer    DEBUG_JFIT
      parameter (DEBUG_JFIT     = 5)
      integer    DEBUG_NR
      parameter (DEBUG_NR       = 6)

      integer    DEBUG_JBAS
      parameter (DEBUG_JBAS     = 7)
      integer    DEBUG_KBAS
      parameter (DEBUG_KBAS     = 8)
      integer    DEBUG_AOBAS
      parameter (DEBUG_AOBAS    = 9)

      integer    DEBUG_FORCES
      parameter (DEBUG_FORCES   = 10)

      integer    DEBUG_TIMING
      parameter (DEBUG_TIMING   = 11)

      integer    DEBUG_CONTROL
      parameter (DEBUG_CONTROL  = 12)

      integer    DEBUG_MEMORY
      parameter (DEBUG_MEMORY   = 13)

      integer    DEBUG_QUAD
      parameter (DEBUG_QUAD     = 14)

      integer    DEBUG_PARALLEL
      parameter (DEBUG_PARALLEL = 15)

      integer    DEBUG_CHF_RHS
      parameter (DEBUG_CHF_RHS  = 16)
      integer    DEBUG_CHF_LHS
      parameter (DEBUG_CHF_LHS  = 17)
      integer    DEBUG_CHF_DKSM
      parameter (DEBUG_CHF_DKSM = 18)
      integer    DEBUG_DKSM_EXP
      parameter (DEBUG_DKSM_EXP = 19)
      integer    DEBUG_HESS
      parameter (DEBUG_HESS     = 20)

      logical active_sw
      logical ccpdft_sw
      logical abort_sw
      integer print_stack_depth
      integer MAX_PRINT_STACK
      parameter (MAX_PRINT_STACK=10)
      integer current_print_level(MAX_PRINT_STACK)
      common/pauls/
     &     current_print_level,
     &	   print_stack_depth,
     &     active_sw, ccpdft_sw, abort_sw
c
c we need a parameter for stating that something is undefined
c
      integer DFT_UNDEF
      parameter (DFT_UNDEF=-1)
c
c legitimate choices for weight_scheme
c
      integer WT_BECKE, WT_BECKESCR, WT_SSF, WT_SSFSCR, WT_MHL,
     +        WT_MHL4SSFSCR, WT_MHL8SSFSCR, WT_MHLSCR
      parameter (WT_BECKE=1)
      parameter (WT_BECKESCR=2)
      parameter (WT_SSF=3)
      parameter (WT_SSFSCR=4)
      parameter (WT_MHL=5)
      parameter (WT_MHLSCR=6)
      parameter (WT_MHL4SSFSCR = 7)
      parameter (WT_MHL8SSFSCR = 8)
c
c legitimate choices for grids (ie. based on the terms for which they
c                               are used).
c     G_KS:   the "normal" Kohn-Sham grid
c     G_CPKS: the grid to used for the Coupled Perturbed Kohn-Sham 
c             equations
c
      integer G_KS, G_CPKS
      parameter (G_KS=1)
      parameter (G_CPKS=2)
c
c legitimate choices for the angular grid pruning schemes
c
c     DFT_UNDEF: Undefined pruning scheme
c     AP_NONE: No pruning of the angular grid (has been replaced by
c              AP_RADZONE with 1 radial zone).
c     AP_MHL:  Pruning of angular grid as suggested by Murray, Handy 
c              and Laming
c     AP_AUTO: Pruning of angular grid according to obtained energies
c              (automatic)
c     AP_RADZONE: Pruning of angular grid using user specified numbers
c              of angular grid points for each radial domain.
c     AP_SG1:  Pruning of angular grid according to SG1 specification
c     AP_SG1a: Pruning of angular grid according to modified SG1 
c              specification
c     
      integer AP_MHL, AP_RADZONE, AP_SG1, AP_SG1a, AP_AUTO
      parameter (AP_MHL=11)
      parameter (AP_RADZONE=12)
      parameter (AP_SG1=13)
      parameter (AP_SG1a=14)
      parameter (AP_AUTO=15)
c
c legitimate choices for the radial grid schemes
c
c     DFT_UNDEF: Undefined radial grid
c     RG_MK: Mura & Knowles logarithmic grid
c     RG_EML: Murray, Handy and Lamings Euler-MacLaurin grid
c     RG_B: The Becke radial grid
c     RG_SG1: The SG1 radial grid (which is EML with special scale 
c             factors)
c
      integer RG_MK, RG_EML, RG_B, RG_SG1
      parameter (RG_MK=21)
      parameter (RG_EML=22)
      parameter (RG_B=23)
      parameter (RG_SG1=24)
c
c legitimate choices for the radial grid scale factor (grid_atom_radius)
c
      integer SC_MK, SC_GAM1, SC_GAM2
      parameter (SC_MK=31)
      parameter (SC_GAM1=32)
      parameter (SC_GAM2=33)
c
c legitimate choices for the angular grid schemes
c
c     DFT_UNDEF: Undefined angular grid
c     AG_LEB: Lebedev-Laikov angular grids
c     AG_LEG: Gauss-Legendre angular grids.
c
      integer AG_LEB, AG_LEG
      parameter (AG_LEB=41)
      parameter (AG_LEG=42)
c
c legitimate choices for the grid accuracy schemes
c
c     DFT_UNDEF:       Undefined grid accuracy
c     GACC_LOW:        Low accuracy predefined grid
c     GACC_LOWMEDIUM:  Low-medium accuracy predefined grid
c     GACC_MEDIUM:     Medium accuracy predefined grid
c     GACC_MEDIUMHIGH: Medium-high accuracy predefined grid
c     GACC_HIGH:       High accuracy predefined grid
c     GACC_VERYHIGH:   Very high accuracy predefined grid
c     GACC_REF:        Reference grid
c     GACC_SG1:        SG1 grid
c
      integer GACC_LOW, GACC_LOWMEDIUM, GACC_MEDIUM
      integer GACC_MEDIUMHIGH, GACC_HIGH, GACC_VERYHIGH, GACC_REF
      integer GACC_SG1
      parameter (GACC_LOW       = 51)
      parameter (GACC_LOWMEDIUM = 52)
      parameter (GACC_MEDIUM    = 53)
      parameter (GACC_MEDIUMHIGH= 54)
      parameter (GACC_HIGH      = 55)
      parameter (GACC_VERYHIGH  = 56)
      parameter (GACC_REF       = 57)
      parameter (GACC_SG1       = 58)
      integer iwr
      common/dft_iofile/iwr
C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
c
c integral accuracy parameters
c
c for energy integrals
c (if -1 use generic parameters)
c
c  icut_2c
c  icut_3c
c  icut_4c [ not yet used ]
c
c  itol_2c
c  itol_3c
c  itol_4c [ not yet used ]
c
c For derivative integrals code
c (if -1 use generic parameters)
c
c  icutd_2c
c  icutd_3c
c  icutd_4c
c
c  itold_2c
c  itold_3c
c  itold_4c
c
c  Generic parameters
c
c  icut_dft
c  itol_dft
c  icutd_dft
c  itold_dft

      integer icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

      common/dft_intctl/ icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

c
c for documentation of these functions see the start
c of dft/basis.m
c
      integer BL_create_atomtag
      external BL_create_atomtag

      integer BL_find_atomtag
      external BL_find_atomtag

      integer BL_import_shell
      external BL_import_shell

      integer BL_assign_types_by_z
      external BL_assign_types_by_z

      integer BL_assign_type
      external BL_assign_type

      integer BL_write_basis
      external BL_write_basis

      integer BL_clear_basis_set
      external BL_clear_basis_set

      integer BL_maxang_on_atom
      external BL_maxang_on_atom

      integer BL_basis_size
      external BL_basis_size

      integer BL_max_shell_count
      external BL_max_shell_count

      integer BL_num_sets
      external BL_num_sets

      integer BL_num_types
      external BL_num_types

      integer BL_get_atom_type
      external BL_get_atom_type

      logical BL_atomtyp_exist
      external BL_atomtyp_exist

      integer BL_summarise
      external BL_summarise
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_ao
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)

c
      integer ao_bas, cd_bas, nu_bas
      common/mdendata/ao_bas, cd_bas, nu_bas

      real*8 memory_fp(*)
      integer iout
      logical debug
      integer allocate_memory

      integer ierror, i

      logical opg_root
      logical rhf_sw, opshell_sw

      rhf_sw = .true.
      opshell_sw = .false.
      rks_sw  = .true.

c      if (.not. active_sw) then
c        gden_init = 0
c        return
c      endif

      icut_dft  = 9
      itol_dft  = 20
      icutd_dft = 10
      itold_dft = 20

      icut_2c = 0  
      itol_2c = 0
      icutd_2c = 0
      itold_2c = 0
      icut_3c = 0 
      itol_3c = 0
      icutd_3c = 0
      itold_3c = 0
      icut_4c = 0
      itol_4c = 0
      icutd_4c = 0
      itold_4c = 0
c
c  set default tolerances for integral codes
c  if the user hasnt explicit set them
c
      if(icut_2c .eq. 0)icut_2c = icut_dft
      if(itol_2c .eq. 0)itol_2c = itol_dft
      if(icut_3c .eq. 0)icut_3c = icut_dft
      if(itol_3c .eq. 0)itol_3c = itol_dft
      if(icut_4c .eq. 0)icut_4c = icut_dft
      if(itol_4c .eq. 0)itol_4c = itol_dft
      if(icutd_2c .eq. 0)icutd_2c = icutd_dft
      if(itold_2c .eq. 0)itold_2c = itold_dft
      if(icutd_3c .eq. 0)icutd_3c = icutd_dft
      if(itold_3c .eq. 0)itold_3c = itold_dft
      if(icutd_4c .eq. 0)icutd_4c = icutd_dft
      if(itold_4c .eq. 0)itold_4c = itold_dft

      if(debug)write(iout,*)'3-centre tolerances',icut_3c, itol_3c
c
c  print control
c
      do i=1,MAX_DEBUG
        print_sw(i) = .false.
      enddo
      debug_sw = .false.
c
c  print stack
c
      print_stack_depth = 1
      current_print_level(1) = PRINT_DEFAULT
c
c output stream
c
      iwr = iout
C 
C Write out basis sets, if requested.
c
c !!! Need to set these tags more carefully 
c
      ao_bas = 1
c model density will be basis 2
      cd_bas = 2
c nuclear density will be basis 3
      nu_bas = 3

      call order_fill
      call BL_init
      call interface_gamess(rhf_sw,opshell_sw,ao_bas,iout)
C 
C Write out AO basis sets, if requested.
c
ccc   if(opg_root() .and. print_sw(DEBUG_AOBAS))then

      if(debug)then
         write(iout,*)'AO  Basis'
         ierror=BL_write_basis(ao_bas,iout)
      endif

      if(debug)write(iout,*)'3-centre tolerances end',icut_3c, itol_3c

      gden_init = 0
      return
      end

c**********************************************************************
c
c  gden_add_gauss(ncent, iout, debug)
c
c   ncent   (in)   number of sites to consider
c   iout    (in)   unit number for listing 
c   debug   (in)   debug output flag 
c
c  in common/blur
c   alpha   (in)   exponents (negative values denote no model fn
c                  on the site)
c   wght    (in)   weight - i.e. charge being blurred
c  
c  by this point any control settings should have been made
c
c**********************************************************************

      subroutine gden_add_gauss(ncent,iout,debug)

      implicit none

      integer ncent,iout
      logical debug

C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
c for documentation of these functions see the start
c of dft/basis.m
c
      integer BL_create_atomtag
      external BL_create_atomtag

      integer BL_find_atomtag
      external BL_find_atomtag

      integer BL_import_shell
      external BL_import_shell

      integer BL_assign_types_by_z
      external BL_assign_types_by_z

      integer BL_assign_type
      external BL_assign_type

      integer BL_write_basis
      external BL_write_basis

      integer BL_clear_basis_set
      external BL_clear_basis_set

      integer BL_maxang_on_atom
      external BL_maxang_on_atom

      integer BL_basis_size
      external BL_basis_size

      integer BL_max_shell_count
      external BL_max_shell_count

      integer BL_num_sets
      external BL_num_sets

      integer BL_num_types
      external BL_num_types

      integer BL_get_atom_type
      external BL_get_atom_type

      logical BL_atomtyp_exist
      external BL_atomtyp_exist

      integer BL_summarise
      external BL_summarise
C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres

      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer mxbltyp
      parameter(mxbltyp = 50)

      logical oblur
      common/blurrr/oblur

      integer nbltyp, bltab(maxat), nutab(maxat)

      real*8 blexpo(maxat), blwght(maxat)
      real*8 blexpo2(mxbltyp), blwght2(mxbltyp)

      common/blurr2/blexpo, blwght,
     &     blexpo2, blwght2,
     &     bltab,nbltyp,nutab

      character*8 ztagbl
      common/blurrc/ztagbl(mxbltyp)

c
c integral accuracy parameters
c
c for energy integrals
c (if -1 use generic parameters)
c
c  icut_2c
c  icut_3c
c  icut_4c [ not yet used ]
c
c  itol_2c
c  itol_3c
c  itol_4c [ not yet used ]
c
c For derivative integrals code
c (if -1 use generic parameters)
c
c  icutd_2c
c  icutd_3c
c  icutd_4c
c
c  itold_2c
c  itold_3c
c  itold_4c
c
c  Generic parameters
c
c  icut_dft
c  itol_dft
c  icutd_dft
c  itold_dft

      integer icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

      common/dft_intctl/ icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

c
      integer ao_bas, cd_bas, nu_bas
      common/mdendata/ao_bas, cd_bas, nu_bas

      integer ibas, i, type, atomno, angm, hybr
      integer nprm_count, itype, ierror
      real*8 ex(10),cs(10),cp(10),cd(10),cf(10),cg(10)

      integer itest, idum, icount, icount2

      integer inuc, inullnuc, inullcd

c this is kept internally just for houskeeping
      real*8 expon(max_atype)
c
c External functions
c
      logical opg_root

      if(debug) write(iout,*) 'Entered gden_add_gauss'
c
      ierror = BL_clear_basis_set(cd_bas)
      ierror = BL_clear_basis_set(nu_bas)

c
c Create three atom tags 
c
c inullnuc - nuclear basis tag for atoms which to 
c            not have nuclear charges (blurred centres)
c
c inuc     - a single S function for nuclear charges
c            represented as a nuclear basis
c 
c inullcd  - a CD basis for atoms with no CD basis
c            functions - all un-blurred
c            point charges and real atoms
c

      atomno = -1

      inullnuc=BL_create_atomtag(nu_bas,atomno)
      inuc=BL_create_atomtag(nu_bas,atomno)

      inullcd=BL_create_atomtag(cd_bas,atomno)

c expon is used to remove the need for duplicate
c CD tags if atoms share a blur exponent
c the null tag won't be used for this

      expon(inullcd) = -1.0

      if(debug)then
         write(iout,*)'atom tags'
         write(iout,*)'null cd tag  = ',inullcd
         write(iout,*)'nuclear tag  = ',inullcd
         write(iout,*)'null nuc tag = ',inullnuc
      endif
c
c Check in a single S basis function, we will use
c this for all nuclei
c
      nprm_count = 1
      angm = 1
      hybr = 1
      cs(nprm_count)=1.0d0
      cp(nprm_count)=1.0d0
      cd(nprm_count)=1.0d0
      cf(nprm_count)=1.0d0
      cg(nprm_count)=1.0d0
      ex(nprm_count)=100.0d0
      ierror=BL_import_shell(nu_bas,inuc,nprm_count,angm,hybr,
     &     ex,cs,cp,cd,cf,cg)

      icount = 0
      icount2 = 0

      do i = 1,ncent
c
c Is there a blurred charge on the centre
c
         if(debug)write(iout,*)'atom',i,' blur ',blexpo(i)

         if(blexpo(i).gt.-1.0d-10)then
c
c Check if this exponent has been used before
c
            type = -1
            do itype = 1,BL_num_types(cd_bas)
               if (abs(expon(itype)-blexpo(i)) .lt. 1.0d-10)then
                  type=itype
               endif
            enddo

            if(type .eq. -1)then
c
c Create new atom type (use -1 as this is not a real atom)
c
               atomno = -1
               type=BL_create_atomtag(cd_bas,atomno)
               if(debug)
     &         write(iout,*) 'creating CD function exp=',blexpo(i),type
c
c Check in a single S basis function
c Use weight as the coefficient
c
               nprm_count = 1
               angm = 1
               hybr = 1
               cs(nprm_count)=1.0d0
               cp(nprm_count)=1.0d0
               cd(nprm_count)=1.0d0
               cf(nprm_count)=1.0d0
               cg(nprm_count)=1.0d0
               ex(nprm_count)=blexpo(i)
               ierror=BL_import_shell(cd_bas,type,nprm_count,angm,hybr,
     &              ex,cs,cp,cd,cf,cg)
               expon(type) = blexpo(i)
            endif
            if(debug)write(iout,*) 'assigning cd type', type
            idum = BL_assign_type(cd_bas,i,type)
            icount = icount + 1
            bltab(icount) = i
         else
c
c assign null type 
c     
            if(debug)write(iout,*) 'assigning null cd'
            idum = BL_assign_type(cd_bas,i,inullcd)
         endif

c
c assign nuclear density basis
c
         if(ian(i) .gt. 0)then
            if(debug)write(iout,*) 'assigning nuc type',i
            idum = BL_assign_type(nu_bas,i,inuc)
            icount2 = icount2 + 1
            nutab(icount2) = i
         else
            if(debug)write(iout,*) 'assigning null nuc type',i
            idum = BL_assign_type(nu_bas,i,inullnuc)
         endif
      enddo

      write(iout,*)'blurred functions assigned=',icount
      write(iout,*)'nuclear functions assigned=',icount2

      call checkin_basis(cd_bas, debug)
      call checkin_basis(nu_bas, debug)

      if(debug)then
         write(iout,*)'Blurred charge  Basis'
         if(opg_root())ierror=BL_write_basis(cd_bas,iout)

         write(iout,*)'Nuclear charge  Basis'
         if(opg_root())ierror=BL_write_basis(nu_bas,iout)

         ierror = BL_summarise()
      endif

      return
      end

***********************************************************************
*
*  gden_energy
*
* Description:                                                        *
*                                                                     *
* Form the Kohn-Sham matrix contributions from the blurred charges    *
* Also construct contributions to the energy                          *
*
*  e_nuc ... blurred-charge... nuclei
*  e_elec .. blurred-charge... electrons
*
*
* NB - e_elec is diagnostic only as ehf1 will include this            *
* 
***********************************************************************

      integer function gden_energy(
     &     coords,
     &     kma,kmb,adenm,bdenm,e_elec,e_nuc,
     &     memory_int,memory_fp,print,debug,iout)

      implicit none
c
C In variables
c
      real*8 coords(3,*)
      real*8 kma(*),kmb(*),adenm(*),bdenm(*)
      logical print, debug
      integer iout
c
c Out variables
c
      real*8 e_elec, e_nuc

C *Scratch space and pointers
      integer memory_int(*)
      real*8 memory_fp(*)

C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
C *API interface data storage
      real*8 alpha_Den 
      real*8 beta_Den
      real*8 J_energy
      real*8 lmult
      real*8 totDen
      integer totPts
      real*8 XC_energy
      common/results_api/alpha_Den,beta_Den,J_energy,lmult,totDen,
     &                   XC_energy,totpts

C *Inter-module communication variables.
C *Suffixes and what they mean		
C *---------------------------	
C *_sw			-		switch		
C *_num 		-		numbers	
C *_ch			-		Input/output channels
      integer out_ch,in_ch
      common/io_channels/out_ch,in_ch
C *
C *Global switches
C *
      logical debug_sw
      common/global_switches/debug_sw
C *
C *Switches and numbers used in dft routines
C *
      logical optim_sw,triangle_sw
      common/scf_control_switch/optim_sw,triangle_sw

      logical jfit_sw,jfitg_sw,cmm_sw,dunlap_sw,potential_sw
      logical kqua_sw,kfit_sw
      logical rks_sw
      logical ludm_sw,svdm_sw
      logical jown_sw,dega_sw,kown_sw
      logical mult_sw, dft2e_sw
      common/scftype/rks_sw
      common/j_switch/jfit_sw,jfitg_sw,cmm_sw,mult_sw,
     &                dunlap_sw,potential_sw,dft2e_sw

      common/xc_switch/kqua_sw,kfit_sw,
     &     ludm_sw,svdm_sw,jown_sw,dega_sw,kown_sw
c
c The grid parameters
c
c     1) SG1 fully specifies everything
c     2) rad_grid_scheme specifies for each type which radial grid
c        is used
c       -1) if RG_MK then
c              radm_num specifies m
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -2) if RG_EML then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -3) if RG_B then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c     3) ang_grid_scheme specifies which angular grid to use
c       -1) if AG_LEB then
c              angupt_num specifies the maximum number of angular grid
c              points
c       -2) if AG_LEG then
c              thetpt_num specifies the maximum number of theta points
c              phipt_num specifies the maximum number of phi points
c     4) ang_prune_scheme specifies which scheme to use for pruning 
c        the angular grid as a function of the radius.
c       -1) if AP_MHL (no other info needed)
c       -2) if AP_RADZONE then
c              radzones_num specifies the number of radial zones
c              bnd_radzn specifies the location of zone boundaries
c              angpt_radzn_num specifies the number of angular grid 
c              points per zone.
c              
c     integer angupt_num,thetpt_num,phipt_num,radpt_num
      integer radpt_num
      integer weight_scheme, radzones_num, angpt_radzn_num
      integer thetpt_radzn_num, phipt_radzn_num
      integer ang_prune_scheme
      integer rad_grid_scheme, ang_grid_scheme
      integer gtype_num, ngtypes, gaccu_num
      integer iauto_cnt
      integer rad_scale_scheme
      integer radnpt_row
      integer angnpt_row
      integer grid_generation
      real*8 grid_scale, radm_num, bnd_radzn
      real*8 grid_atom_radius
      real*8 weight_atom_radius
      real*8 prune_atom_radius
      real*8 screen_atom_radius
c
c     subtle difference here: 
c
c     - grid_atom_radius:   used as scale factor of the radial grids.
c     - weight_atom_radius: used for atom-size-adjustments in the
c                           weighting scheme.
c     - prune_atom_radius:  used for pruning the angular grid in the
c                           MHL pruning scheme
c     - screen_atom_radius: used for screening of the radial grids.
c
      real*8 psitol, warntol
      logical conv_prune_sw, gradwght_sw, sort_points_sw
      logical ignore_accuracy_sw
      common/dft_grid_parameters/
     &     psitol(0:max_gtype,max_grids),
     &     warntol(0:max_gtype,max_grids),
     &     radm_num(0:max_gtype,max_grids),
     &     grid_scale(0:max_gtype,max_grids),
     &     grid_atom_radius(0:max_gtype,max_grids),
     &     weight_atom_radius(0:max_gtype,max_grids),
     &     prune_atom_radius(0:max_gtype,max_grids),
     &     screen_atom_radius(0:max_gtype,max_grids),
     &     bnd_radzn(maxradzn-1,0:max_gtype,max_grids),
     &     radnpt_row(7),angnpt_row(7),
     &     angpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     thetpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     phipt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     radzones_num(0:max_gtype,max_grids),
     &     ang_prune_scheme(0:max_gtype,max_grids),
     &     rad_grid_scheme(0:max_gtype,max_grids),
     &     ang_grid_scheme(0:max_gtype,max_grids),
     &     radpt_num(0:max_gtype,max_grids),
     &     gaccu_num(0:max_gtype,max_grids),
     &     gtype_num(max_atom),
     &     grid_generation,
     &     ngtypes,iauto_cnt,
     &     rad_scale_scheme,
     &     weight_scheme(max_grids),
     &     conv_prune_sw,
     &     gradwght_sw,
     &     sort_points_sw,
     &     ignore_accuracy_sw

      integer poleexp_num,over_tol,pener_tol,schwarz_tol
      real*8  tttt2
      common/pole_options/tttt2,poleexp_num,over_tol,pener_tol,
     &                    schwarz_tol

      integer    MAX_DEBUG
      parameter (MAX_DEBUG=25)

      logical print_sw(MAX_DEBUG)
      common/debugpr/print_sw
c
c debug array indices
c
      integer    DEBUG_KSMATRIX
      parameter (DEBUG_KSMATRIX = 1)
      integer    DEBUG_TR
      parameter (DEBUG_TR       = 2)
      integer    DEBUG_NORM
      parameter (DEBUG_NORM     = 3)
      integer    DEBUG_DENSITY
      parameter (DEBUG_DENSITY  = 4)
      integer    DEBUG_JFIT
      parameter (DEBUG_JFIT     = 5)
      integer    DEBUG_NR
      parameter (DEBUG_NR       = 6)

      integer    DEBUG_JBAS
      parameter (DEBUG_JBAS     = 7)
      integer    DEBUG_KBAS
      parameter (DEBUG_KBAS     = 8)
      integer    DEBUG_AOBAS
      parameter (DEBUG_AOBAS    = 9)

      integer    DEBUG_FORCES
      parameter (DEBUG_FORCES   = 10)

      integer    DEBUG_TIMING
      parameter (DEBUG_TIMING   = 11)

      integer    DEBUG_CONTROL
      parameter (DEBUG_CONTROL  = 12)

      integer    DEBUG_MEMORY
      parameter (DEBUG_MEMORY   = 13)

      integer    DEBUG_QUAD
      parameter (DEBUG_QUAD     = 14)

      integer    DEBUG_PARALLEL
      parameter (DEBUG_PARALLEL = 15)

      integer    DEBUG_CHF_RHS
      parameter (DEBUG_CHF_RHS  = 16)
      integer    DEBUG_CHF_LHS
      parameter (DEBUG_CHF_LHS  = 17)
      integer    DEBUG_CHF_DKSM
      parameter (DEBUG_CHF_DKSM = 18)
      integer    DEBUG_DKSM_EXP
      parameter (DEBUG_DKSM_EXP = 19)
      integer    DEBUG_HESS
      parameter (DEBUG_HESS     = 20)

      logical active_sw
      logical ccpdft_sw
      logical abort_sw
      integer print_stack_depth
      integer MAX_PRINT_STACK
      parameter (MAX_PRINT_STACK=10)
      integer current_print_level(MAX_PRINT_STACK)
      common/pauls/
     &     current_print_level,
     &	   print_stack_depth,
     &     active_sw, ccpdft_sw, abort_sw
c
c we need a parameter for stating that something is undefined
c
      integer DFT_UNDEF
      parameter (DFT_UNDEF=-1)
c
c legitimate choices for weight_scheme
c
      integer WT_BECKE, WT_BECKESCR, WT_SSF, WT_SSFSCR, WT_MHL,
     +        WT_MHL4SSFSCR, WT_MHL8SSFSCR, WT_MHLSCR
      parameter (WT_BECKE=1)
      parameter (WT_BECKESCR=2)
      parameter (WT_SSF=3)
      parameter (WT_SSFSCR=4)
      parameter (WT_MHL=5)
      parameter (WT_MHLSCR=6)
      parameter (WT_MHL4SSFSCR = 7)
      parameter (WT_MHL8SSFSCR = 8)
c
c legitimate choices for grids (ie. based on the terms for which they
c                               are used).
c     G_KS:   the "normal" Kohn-Sham grid
c     G_CPKS: the grid to used for the Coupled Perturbed Kohn-Sham 
c             equations
c
      integer G_KS, G_CPKS
      parameter (G_KS=1)
      parameter (G_CPKS=2)
c
c legitimate choices for the angular grid pruning schemes
c
c     DFT_UNDEF: Undefined pruning scheme
c     AP_NONE: No pruning of the angular grid (has been replaced by
c              AP_RADZONE with 1 radial zone).
c     AP_MHL:  Pruning of angular grid as suggested by Murray, Handy 
c              and Laming
c     AP_AUTO: Pruning of angular grid according to obtained energies
c              (automatic)
c     AP_RADZONE: Pruning of angular grid using user specified numbers
c              of angular grid points for each radial domain.
c     AP_SG1:  Pruning of angular grid according to SG1 specification
c     AP_SG1a: Pruning of angular grid according to modified SG1 
c              specification
c     
      integer AP_MHL, AP_RADZONE, AP_SG1, AP_SG1a, AP_AUTO
      parameter (AP_MHL=11)
      parameter (AP_RADZONE=12)
      parameter (AP_SG1=13)
      parameter (AP_SG1a=14)
      parameter (AP_AUTO=15)
c
c legitimate choices for the radial grid schemes
c
c     DFT_UNDEF: Undefined radial grid
c     RG_MK: Mura & Knowles logarithmic grid
c     RG_EML: Murray, Handy and Lamings Euler-MacLaurin grid
c     RG_B: The Becke radial grid
c     RG_SG1: The SG1 radial grid (which is EML with special scale 
c             factors)
c
      integer RG_MK, RG_EML, RG_B, RG_SG1
      parameter (RG_MK=21)
      parameter (RG_EML=22)
      parameter (RG_B=23)
      parameter (RG_SG1=24)
c
c legitimate choices for the radial grid scale factor (grid_atom_radius)
c
      integer SC_MK, SC_GAM1, SC_GAM2
      parameter (SC_MK=31)
      parameter (SC_GAM1=32)
      parameter (SC_GAM2=33)
c
c legitimate choices for the angular grid schemes
c
c     DFT_UNDEF: Undefined angular grid
c     AG_LEB: Lebedev-Laikov angular grids
c     AG_LEG: Gauss-Legendre angular grids.
c
      integer AG_LEB, AG_LEG
      parameter (AG_LEB=41)
      parameter (AG_LEG=42)
c
c legitimate choices for the grid accuracy schemes
c
c     DFT_UNDEF:       Undefined grid accuracy
c     GACC_LOW:        Low accuracy predefined grid
c     GACC_LOWMEDIUM:  Low-medium accuracy predefined grid
c     GACC_MEDIUM:     Medium accuracy predefined grid
c     GACC_MEDIUMHIGH: Medium-high accuracy predefined grid
c     GACC_HIGH:       High accuracy predefined grid
c     GACC_VERYHIGH:   Very high accuracy predefined grid
c     GACC_REF:        Reference grid
c     GACC_SG1:        SG1 grid
c
      integer GACC_LOW, GACC_LOWMEDIUM, GACC_MEDIUM
      integer GACC_MEDIUMHIGH, GACC_HIGH, GACC_VERYHIGH, GACC_REF
      integer GACC_SG1
      parameter (GACC_LOW       = 51)
      parameter (GACC_LOWMEDIUM = 52)
      parameter (GACC_MEDIUM    = 53)
      parameter (GACC_MEDIUMHIGH= 54)
      parameter (GACC_HIGH      = 55)
      parameter (GACC_VERYHIGH  = 56)
      parameter (GACC_REF       = 57)
      parameter (GACC_SG1       = 58)
c
c for documentation of these functions see the start
c of dft/basis.m
c
      integer BL_create_atomtag
      external BL_create_atomtag

      integer BL_find_atomtag
      external BL_find_atomtag

      integer BL_import_shell
      external BL_import_shell

      integer BL_assign_types_by_z
      external BL_assign_types_by_z

      integer BL_assign_type
      external BL_assign_type

      integer BL_write_basis
      external BL_write_basis

      integer BL_clear_basis_set
      external BL_clear_basis_set

      integer BL_maxang_on_atom
      external BL_maxang_on_atom

      integer BL_basis_size
      external BL_basis_size

      integer BL_max_shell_count
      external BL_max_shell_count

      integer BL_num_sets
      external BL_num_sets

      integer BL_num_types
      external BL_num_types

      integer BL_get_atom_type
      external BL_get_atom_type

      logical BL_atomtyp_exist
      external BL_atomtyp_exist

      integer BL_summarise
      external BL_summarise
C *Contains pointers and information about issued memory.
c  Used for general book keeping of memory.
      integer iblock_info,start_ifreemem,inum_issued,iblock_amt
      integer dblock_info,start_dfreemem,dnum_issued,dblock_amt
      integer last_block
      common/memory_info/iblock_info(max_block),start_ifreemem,
     &                   inum_issued,iblock_amt(max_block),
     &                   dblock_info(max_block),start_dfreemem,
     &                   dnum_issued,dblock_amt(max_block),
     &                   last_block
c
c integral accuracy parameters
c
c for energy integrals
c (if -1 use generic parameters)
c
c  icut_2c
c  icut_3c
c  icut_4c [ not yet used ]
c
c  itol_2c
c  itol_3c
c  itol_4c [ not yet used ]
c
c For derivative integrals code
c (if -1 use generic parameters)
c
c  icutd_2c
c  icutd_3c
c  icutd_4c
c
c  itold_2c
c  itold_3c
c  itold_4c
c
c  Generic parameters
c
c  icut_dft
c  itol_dft
c  icutd_dft
c  itold_dft

      integer icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

      common/dft_intctl/ icut_dft,  itol_dft, icutd_dft, itold_dft,
     &     icut_2c,  itol_2c, icutd_2c, itold_2c,
     &     icut_3c,  itol_3c, icutd_3c, itold_3c,
     &     icut_4c,  itol_4c, icutd_4c, itold_4c

C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer mxbltyp
      parameter(mxbltyp = 50)

      logical oblur
      common/blurrr/oblur

      integer nbltyp, bltab(maxat), nutab(maxat)

      real*8 blexpo(maxat), blwght(maxat)
      real*8 blexpo2(mxbltyp), blwght2(mxbltyp)

      common/blurr2/blexpo, blwght,
     &     blexpo2, blwght2,
     &     bltab,nbltyp,nutab

      character*8 ztagbl
      common/blurrc/ztagbl(mxbltyp)
c
      integer ao_bas, cd_bas, nu_bas
      common/mdendata/ao_bas, cd_bas, nu_bas

c to suppres output from integral routines
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic

C *Local variables
      integer i,j,iscr,iscrb,ltri
      real*8 grad_dum
      integer xbfn_num
      integer idum, nsh, n_ao, n_cd, n_nu
      integer gout_pt,eri_pt,iso_pt,coef_pt, coef2_pt
      integer imode
      real*8 dum
      real*8 gamma, fac, charge
      integer i2c2e_pt
      integer nprint_keep
C *Functions
      integer allocate_memory
      logical opg_root
      real*8 tracep      
      integer cd_update_geom

      nprint_keep = nprint

      if(debug)then
         write(iout,*) 'Entering gden_energy'
         write(iout,*)'3-centre tolerances',icut_3c, itol_3c
         write(iout,*)'2-centre tolerances',icut_2c, itol_2c
      else
         nprint = -5
      endif
c
c  Refresh coordinates
c
      idum = cd_update_geom(coords)

      nsh = BL_max_shell_count()
      n_ao = BL_basis_size(ao_bas)
      n_cd = BL_basis_size(cd_bas)
      n_nu = BL_basis_size(nu_bas)

      ltri=((n_ao+1)*n_ao)/2

      if(opg_root() .and. debug)then
         write(iout,*)'** Input Density **'
         call prtri(adenm,n_ao)
      endif

      iscr = allocate_memory(ltri,'d')
      call aclear_dp(memory_fp(iscr),ltri,0.0d0)

      iso_pt  = allocate_memory(nsh,'d')
      gout_pt = allocate_memory(50625,'d')
      eri_pt  = allocate_memory(225*n_cd,'d')
      coef_pt = allocate_memory(n_cd,'d')
      coef2_pt = allocate_memory(n_nu,'d')
c 
c Note -1 so that model density is stored in units of
c the electron charge
c
      do i = 1, n_cd
         gamma = blexpo(bltab(i))
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef_pt + i -1) = - fac*blwght(bltab(i))
         if(debug)write(iout,*)'Coef ',i, memory_fp(coef_pt + i -1)
      enddo

c
      do i = 1, n_nu
c @@ pseudo alert
	 if(debug)write(iout,*)'Nutab',i,nutab(i)
         charge = ian(nutab(i))
         if(debug)write(iout,*)'ian',i,ian(nutab(i)),charge
         gamma = 100.0d0
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef2_pt + i -1) = -fac*charge
         if(debug)write(iout,*)'Coef2 ',i, memory_fp(coef2_pt + i -1)
      enddo

      i2c2e_pt = allocate_memory(n_cd*n_nu,'d')
      imode = 6
      call jkint_dft(memory_fp(iso_pt),
     &     memory_fp(gout_pt),nsh,
     &     cd_bas, -1, nu_bas, -1, 
     &     imode,memory_fp(i2c2e_pt),dum,
     &     dum,dum,
     &     dum,dum,dum)
c

      if(debug)then
      write(iout,*)'2-centre integrals'
      do i = 1,n_cd
         write(iout,100) (memory_fp(i2c2e_pt + 
     &        n_cd*(j-1) + (i-1) ) *
     &        memory_fp(coef_pt + i -1) * 
     &        memory_fp(coef2_pt + j -1),
     &        j=1,n_nu)
      enddo
      endif

      e_nuc = 0.0d0
      do i = 1,n_cd
         do j = 1,n_nu
            e_nuc= e_nuc + 
     &           memory_fp(i2c2e_pt + n_cd*(j-1) + (i-1) ) *
     &           memory_fp(coef_pt + i-1) * 
     &           memory_fp(coef2_pt + j-1)
         enddo
      enddo

      if(debug)write(iout,*)'e_nuc',e_nuc

 100  format(1x,10f15.8)

      call free_memory(i2c2e_pt,'d')

      imode = 3
      call jkint_dft(memory_fp(iso_pt),
     &     memory_fp(gout_pt),nsh,
     &     ao_bas, ao_bas, cd_bas, -1, 
     &     imode,memory_fp(coef_pt),dum,
     &     memory_fp(iscr),dum,
     &     memory_fp(eri_pt),dum,dum)

      call free_memory(coef2_pt,'d')
      call free_memory(coef_pt,'d')
      call free_memory(eri_pt,'d')
      call free_memory(gout_pt,'d')
      call free_memory(iso_pt,'d')
c
c Compute 3C energy contribution
c
      e_elec=tracep(adenm,memory_fp(iscr),n_ao)
c
c print Vxc
c
      if(opg_root() .and. debug)then
         if(rks_sw)then
            write(iout,*)'** V coul **'
            call prtriprec(memory_fp(iscr),n_ao)
         else
            write(iout,*)'** Vxc alpha **'
            call prtri(memory_fp(iscr),n_ao)
            write(iout,*)'** Vxc beta **'
            call prtri(memory_fp(iscrb),n_ao)
         endif
      endif

c
c sum in Vxc
c
      call daxpy(ltri,1.0d0,memory_fp(iscr),1,kma,1)
      call free_memory(iscr,'d')

      if(.not. rks_sw)then
         call daxpy(ltri,1.0d0,memory_fp(iscrb),1,kmb,1)
         call free_memory(iscrb,'d')
      endif
c
c exit code if requested
c

      if(abort_sw)then
         if(opg_root())then
            write(iout,*)'EXIT REQUESTED'
            call list_time_periods(.false.,.false.)
         endif
         call pg_end
         call exitc(0)
      endif

      GDEN_energy = 0

      if(debug)call gmem_check_guards('end gden_energy')

      if(debug)then
         nprint = nprint_keep
      endif

      end


      subroutine prtriprec(d,n)
c
c     ----- print out a triangular matrix -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
      dimension d(*),dd(12)
      mmax=7
      imax = 0
  100 imin = imax+1
      imax = imax+mmax
      if (imax .gt. n) imax = n
      write (iwr,9008)
      write (iwr,9028) (i,i = imin,imax)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m = max(i,j)*(max(i,j)-1)/2 + min(i,j)
  140 dd(k) = d(m)
      write (iwr,9048) j,(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.4)
      end


C *****************************************************************************
c
c   integer function gden_forces - Driver for Force evaluation
c
c    Arguments:
c
c     coords,        I   current coordinates
c     adenm,         I   alpha or total density 
c     bdenm,         I   beta density, not used.
c     memory_int,    S   base for dynamic memory allocation (integer)
c     memory_fp,     S   base for dynamic memory allocation (real*8)
c     grad,         I/O  gradient 
c
C *****************************************************************************

      integer function gden_forces(coords,adens,bdens,
     &     memory_int,memory_fp,grad,print,debug,iout)
      implicit none

      real*8 coords(3,*)
      real*8 adens(*),bdens(*)
      integer memory_int(*)
      real*8 memory_fp(*)
C *Out variables
      real*8 grad(3,*)

      logical print, debug
      integer iout

C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
C *Inter-module communication variables.
C *Suffixes and what they mean		
C *---------------------------	
C *_sw			-		switch		
C *_num 		-		numbers	
C *_ch			-		Input/output channels
      integer out_ch,in_ch
      common/io_channels/out_ch,in_ch
C *
C *Global switches
C *
      logical debug_sw
      common/global_switches/debug_sw
C *
C *Switches and numbers used in dft routines
C *
      logical optim_sw,triangle_sw
      common/scf_control_switch/optim_sw,triangle_sw

      logical jfit_sw,jfitg_sw,cmm_sw,dunlap_sw,potential_sw
      logical kqua_sw,kfit_sw
      logical rks_sw
      logical ludm_sw,svdm_sw
      logical jown_sw,dega_sw,kown_sw
      logical mult_sw, dft2e_sw
      common/scftype/rks_sw
      common/j_switch/jfit_sw,jfitg_sw,cmm_sw,mult_sw,
     &                dunlap_sw,potential_sw,dft2e_sw

      common/xc_switch/kqua_sw,kfit_sw,
     &     ludm_sw,svdm_sw,jown_sw,dega_sw,kown_sw
c
c The grid parameters
c
c     1) SG1 fully specifies everything
c     2) rad_grid_scheme specifies for each type which radial grid
c        is used
c       -1) if RG_MK then
c              radm_num specifies m
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -2) if RG_EML then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c       -3) if RG_B then
c              radpt_num specifies the number of grid points
c              grid_scale specifies a scale factor
c     3) ang_grid_scheme specifies which angular grid to use
c       -1) if AG_LEB then
c              angupt_num specifies the maximum number of angular grid
c              points
c       -2) if AG_LEG then
c              thetpt_num specifies the maximum number of theta points
c              phipt_num specifies the maximum number of phi points
c     4) ang_prune_scheme specifies which scheme to use for pruning 
c        the angular grid as a function of the radius.
c       -1) if AP_MHL (no other info needed)
c       -2) if AP_RADZONE then
c              radzones_num specifies the number of radial zones
c              bnd_radzn specifies the location of zone boundaries
c              angpt_radzn_num specifies the number of angular grid 
c              points per zone.
c              
c     integer angupt_num,thetpt_num,phipt_num,radpt_num
      integer radpt_num
      integer weight_scheme, radzones_num, angpt_radzn_num
      integer thetpt_radzn_num, phipt_radzn_num
      integer ang_prune_scheme
      integer rad_grid_scheme, ang_grid_scheme
      integer gtype_num, ngtypes, gaccu_num
      integer iauto_cnt
      integer rad_scale_scheme
      integer radnpt_row
      integer angnpt_row
      integer grid_generation
      real*8 grid_scale, radm_num, bnd_radzn
      real*8 grid_atom_radius
      real*8 weight_atom_radius
      real*8 prune_atom_radius
      real*8 screen_atom_radius
c
c     subtle difference here: 
c
c     - grid_atom_radius:   used as scale factor of the radial grids.
c     - weight_atom_radius: used for atom-size-adjustments in the
c                           weighting scheme.
c     - prune_atom_radius:  used for pruning the angular grid in the
c                           MHL pruning scheme
c     - screen_atom_radius: used for screening of the radial grids.
c
      real*8 psitol, warntol
      logical conv_prune_sw, gradwght_sw, sort_points_sw
      logical ignore_accuracy_sw
      common/dft_grid_parameters/
     &     psitol(0:max_gtype,max_grids),
     &     warntol(0:max_gtype,max_grids),
     &     radm_num(0:max_gtype,max_grids),
     &     grid_scale(0:max_gtype,max_grids),
     &     grid_atom_radius(0:max_gtype,max_grids),
     &     weight_atom_radius(0:max_gtype,max_grids),
     &     prune_atom_radius(0:max_gtype,max_grids),
     &     screen_atom_radius(0:max_gtype,max_grids),
     &     bnd_radzn(maxradzn-1,0:max_gtype,max_grids),
     &     radnpt_row(7),angnpt_row(7),
     &     angpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     thetpt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     phipt_radzn_num(maxradzn,0:max_gtype,max_grids),
     &     radzones_num(0:max_gtype,max_grids),
     &     ang_prune_scheme(0:max_gtype,max_grids),
     &     rad_grid_scheme(0:max_gtype,max_grids),
     &     ang_grid_scheme(0:max_gtype,max_grids),
     &     radpt_num(0:max_gtype,max_grids),
     &     gaccu_num(0:max_gtype,max_grids),
     &     gtype_num(max_atom),
     &     grid_generation,
     &     ngtypes,iauto_cnt,
     &     rad_scale_scheme,
     &     weight_scheme(max_grids),
     &     conv_prune_sw,
     &     gradwght_sw,
     &     sort_points_sw,
     &     ignore_accuracy_sw

      integer poleexp_num,over_tol,pener_tol,schwarz_tol
      real*8  tttt2
      common/pole_options/tttt2,poleexp_num,over_tol,pener_tol,
     &                    schwarz_tol

      integer    MAX_DEBUG
      parameter (MAX_DEBUG=25)

      logical print_sw(MAX_DEBUG)
      common/debugpr/print_sw
c
c debug array indices
c
      integer    DEBUG_KSMATRIX
      parameter (DEBUG_KSMATRIX = 1)
      integer    DEBUG_TR
      parameter (DEBUG_TR       = 2)
      integer    DEBUG_NORM
      parameter (DEBUG_NORM     = 3)
      integer    DEBUG_DENSITY
      parameter (DEBUG_DENSITY  = 4)
      integer    DEBUG_JFIT
      parameter (DEBUG_JFIT     = 5)
      integer    DEBUG_NR
      parameter (DEBUG_NR       = 6)

      integer    DEBUG_JBAS
      parameter (DEBUG_JBAS     = 7)
      integer    DEBUG_KBAS
      parameter (DEBUG_KBAS     = 8)
      integer    DEBUG_AOBAS
      parameter (DEBUG_AOBAS    = 9)

      integer    DEBUG_FORCES
      parameter (DEBUG_FORCES   = 10)

      integer    DEBUG_TIMING
      parameter (DEBUG_TIMING   = 11)

      integer    DEBUG_CONTROL
      parameter (DEBUG_CONTROL  = 12)

      integer    DEBUG_MEMORY
      parameter (DEBUG_MEMORY   = 13)

      integer    DEBUG_QUAD
      parameter (DEBUG_QUAD     = 14)

      integer    DEBUG_PARALLEL
      parameter (DEBUG_PARALLEL = 15)

      integer    DEBUG_CHF_RHS
      parameter (DEBUG_CHF_RHS  = 16)
      integer    DEBUG_CHF_LHS
      parameter (DEBUG_CHF_LHS  = 17)
      integer    DEBUG_CHF_DKSM
      parameter (DEBUG_CHF_DKSM = 18)
      integer    DEBUG_DKSM_EXP
      parameter (DEBUG_DKSM_EXP = 19)
      integer    DEBUG_HESS
      parameter (DEBUG_HESS     = 20)

      logical active_sw
      logical ccpdft_sw
      logical abort_sw
      integer print_stack_depth
      integer MAX_PRINT_STACK
      parameter (MAX_PRINT_STACK=10)
      integer current_print_level(MAX_PRINT_STACK)
      common/pauls/
     &     current_print_level,
     &	   print_stack_depth,
     &     active_sw, ccpdft_sw, abort_sw
c
c we need a parameter for stating that something is undefined
c
      integer DFT_UNDEF
      parameter (DFT_UNDEF=-1)
c
c legitimate choices for weight_scheme
c
      integer WT_BECKE, WT_BECKESCR, WT_SSF, WT_SSFSCR, WT_MHL,
     +        WT_MHL4SSFSCR, WT_MHL8SSFSCR, WT_MHLSCR
      parameter (WT_BECKE=1)
      parameter (WT_BECKESCR=2)
      parameter (WT_SSF=3)
      parameter (WT_SSFSCR=4)
      parameter (WT_MHL=5)
      parameter (WT_MHLSCR=6)
      parameter (WT_MHL4SSFSCR = 7)
      parameter (WT_MHL8SSFSCR = 8)
c
c legitimate choices for grids (ie. based on the terms for which they
c                               are used).
c     G_KS:   the "normal" Kohn-Sham grid
c     G_CPKS: the grid to used for the Coupled Perturbed Kohn-Sham 
c             equations
c
      integer G_KS, G_CPKS
      parameter (G_KS=1)
      parameter (G_CPKS=2)
c
c legitimate choices for the angular grid pruning schemes
c
c     DFT_UNDEF: Undefined pruning scheme
c     AP_NONE: No pruning of the angular grid (has been replaced by
c              AP_RADZONE with 1 radial zone).
c     AP_MHL:  Pruning of angular grid as suggested by Murray, Handy 
c              and Laming
c     AP_AUTO: Pruning of angular grid according to obtained energies
c              (automatic)
c     AP_RADZONE: Pruning of angular grid using user specified numbers
c              of angular grid points for each radial domain.
c     AP_SG1:  Pruning of angular grid according to SG1 specification
c     AP_SG1a: Pruning of angular grid according to modified SG1 
c              specification
c     
      integer AP_MHL, AP_RADZONE, AP_SG1, AP_SG1a, AP_AUTO
      parameter (AP_MHL=11)
      parameter (AP_RADZONE=12)
      parameter (AP_SG1=13)
      parameter (AP_SG1a=14)
      parameter (AP_AUTO=15)
c
c legitimate choices for the radial grid schemes
c
c     DFT_UNDEF: Undefined radial grid
c     RG_MK: Mura & Knowles logarithmic grid
c     RG_EML: Murray, Handy and Lamings Euler-MacLaurin grid
c     RG_B: The Becke radial grid
c     RG_SG1: The SG1 radial grid (which is EML with special scale 
c             factors)
c
      integer RG_MK, RG_EML, RG_B, RG_SG1
      parameter (RG_MK=21)
      parameter (RG_EML=22)
      parameter (RG_B=23)
      parameter (RG_SG1=24)
c
c legitimate choices for the radial grid scale factor (grid_atom_radius)
c
      integer SC_MK, SC_GAM1, SC_GAM2
      parameter (SC_MK=31)
      parameter (SC_GAM1=32)
      parameter (SC_GAM2=33)
c
c legitimate choices for the angular grid schemes
c
c     DFT_UNDEF: Undefined angular grid
c     AG_LEB: Lebedev-Laikov angular grids
c     AG_LEG: Gauss-Legendre angular grids.
c
      integer AG_LEB, AG_LEG
      parameter (AG_LEB=41)
      parameter (AG_LEG=42)
c
c legitimate choices for the grid accuracy schemes
c
c     DFT_UNDEF:       Undefined grid accuracy
c     GACC_LOW:        Low accuracy predefined grid
c     GACC_LOWMEDIUM:  Low-medium accuracy predefined grid
c     GACC_MEDIUM:     Medium accuracy predefined grid
c     GACC_MEDIUMHIGH: Medium-high accuracy predefined grid
c     GACC_HIGH:       High accuracy predefined grid
c     GACC_VERYHIGH:   Very high accuracy predefined grid
c     GACC_REF:        Reference grid
c     GACC_SG1:        SG1 grid
c
      integer GACC_LOW, GACC_LOWMEDIUM, GACC_MEDIUM
      integer GACC_MEDIUMHIGH, GACC_HIGH, GACC_VERYHIGH, GACC_REF
      integer GACC_SG1
      parameter (GACC_LOW       = 51)
      parameter (GACC_LOWMEDIUM = 52)
      parameter (GACC_MEDIUM    = 53)
      parameter (GACC_MEDIUMHIGH= 54)
      parameter (GACC_HIGH      = 55)
      parameter (GACC_VERYHIGH  = 56)
      parameter (GACC_REF       = 57)
      parameter (GACC_SG1       = 58)
C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
c
c for documentation of these functions see the start
c of dft/basis.m
c
      integer BL_create_atomtag
      external BL_create_atomtag

      integer BL_find_atomtag
      external BL_find_atomtag

      integer BL_import_shell
      external BL_import_shell

      integer BL_assign_types_by_z
      external BL_assign_types_by_z

      integer BL_assign_type
      external BL_assign_type

      integer BL_write_basis
      external BL_write_basis

      integer BL_clear_basis_set
      external BL_clear_basis_set

      integer BL_maxang_on_atom
      external BL_maxang_on_atom

      integer BL_basis_size
      external BL_basis_size

      integer BL_max_shell_count
      external BL_max_shell_count

      integer BL_num_sets
      external BL_num_sets

      integer BL_num_types
      external BL_num_types

      integer BL_get_atom_type
      external BL_get_atom_type

      logical BL_atomtyp_exist
      external BL_atomtyp_exist

      integer BL_summarise
      external BL_summarise
C *API interface data storage
      real*8 alpha_Den 
      real*8 beta_Den
      real*8 J_energy
      real*8 lmult
      real*8 totDen
      integer totPts
      real*8 XC_energy
      common/results_api/alpha_Den,beta_Den,J_energy,lmult,totDen,
     &                   XC_energy,totpts


cccINCLUDE(common/dft_mbasis)

c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c

      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      integer mxbltyp
      parameter(mxbltyp = 50)

      logical oblur
      common/blurrr/oblur

      integer nbltyp, bltab(maxat), nutab(maxat)

      real*8 blexpo(maxat), blwght(maxat)
      real*8 blexpo2(mxbltyp), blwght2(mxbltyp)

      common/blurr2/blexpo, blwght,
     &     blexpo2, blwght2,
     &     bltab,nbltyp,nutab

      character*8 ztagbl
      common/blurrc/ztagbl(mxbltyp)
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer ao_bas, cd_bas, nu_bas
      common/mdendata/ao_bas, cd_bas, nu_bas
     
c
c Local variables
c
      integer lbasf,i,j,col,ico
      real*8 gamma, fac, charge, dum
      integer coef_pt, coef2_pt
      integer idum
      integer iiso, nsh, n_ao, n_cd, n_nu
      character*8 zrhf
      real*8 grad2(3,maxat)
      real*8 grad3(3,maxat)
      integer nprint_keep

c
c Functions
c
      integer allocate_memory
      integer cd_update_geom
      logical opg_root
      integer ipg_nodeid

      nprint_keep = nprint
      if(debug)then
         write(iout,*) 'Entering model density grad'
      else
         nprint = -5
      endif
c
C *End declarations
c
C    ===== Store current coordinates ====
c
      idum = cd_update_geom(coords)
c
C    ===== Form the gradients due to the coulomb part   =====

      n_cd=BL_basis_size(cd_bas)
      n_ao=BL_basis_size(ao_bas)
      n_nu=BL_basis_size(nu_bas)
C
C Assign atomic charges as fitting coefficients
C (probably will need a normalisation term also)
c
      coef_pt = allocate_memory(n_cd,'d')
      coef2_pt = allocate_memory(n_nu,'d')

      do i = 1, n_cd
         gamma = blexpo(bltab(i))
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef_pt + i -1) = fac*blwght(bltab(i))
         if(debug)write(iout,*)'Coef ',i, memory_fp(coef_pt + i -1)
      enddo

      do i = 1, n_nu
c @@ pseudo alert
         charge = ian(nutab(i))
         gamma = 100.0d0
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef2_pt + i -1) = -fac*charge
         if(debug)write(iout,*)'Coef2 ',i, memory_fp(coef2_pt + i -1)
      enddo

c
c NB iso should be redundant as symmetry is not used, but array is
c still referenced
c
      nsh = BL_max_shell_count()
      iiso = allocate_memory(nsh,'d')

      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
         write(iout,*) 'Forces before blurred charge terms:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad(j,i),j=1,3)
         enddo
 1000    format(i3,3f16.8)
      endif

      do i=1,natoms
         do j = 1,3
            grad2(j,i) = 0.0d0
            grad3(j,i) = 0.0d0
         enddo
      enddo
c
c 3 centre term
c
      zrhf = 'rhf'
      call jkder_dft(zrhf,
     &     memory_fp,
     &     memory_fp(iiso),
     &     nsh,
     &     ao_bas, ao_bas, cd_bas, -1,
     &     adens,bdens,
     &     memory_fp(coef_pt),dum,
     &     grad3)

cc      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
        if(debug)then
         write(iout,*) '3-centre blurred charge forces:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad3(j,i),j=1,3)
         enddo
        endif
c
c 2 centre term (blur, nuclear)
c
      zrhf = 'rhf'
      call jkder_dft(zrhf,
     &     memory_fp,
     &     memory_fp(iiso),
     &     nsh,
     &     cd_bas, -1, nu_bas, -1,
     &     adens,bdens,
     &     memory_fp(coef_pt),
     &     memory_fp(coef2_pt),
     &     grad2)

      if(debug)then
        write(iout,*) '2-centre blurred charge terms:'
        write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
        do i=1,natoms
           write(iout,1000)i,(grad2(j,i),j=1,3)
        enddo
      endif

      do i=1,natoms
         do j = 1,3
            grad(j,i) = grad(j,i) + grad2(j,i) - grad3(j,i)
c            grad(j,i) = grad(j,i) + grad2(j,i)
         enddo
      enddo

      if(debug)then
        write(iout,*) 'cumulated forces:'
        write(iout,*) 
     &     'atom         de/dx             de/y          de/dz'
        do i=1,natoms
           write(iout,1000)i,(grad(j,i),j=1,3)
       enddo
      endif

      call free_memory(iiso,'d')
      call free_memory(coef2_pt,'d')
      call free_memory(coef_pt,'d')

      gden_forces = 0

      if(debug)call gmem_check_guards('end gmem_forces')

      nprint = nprint_keep

      return
      end

      subroutine ver_dft_gden(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gden.m,v $
     +     "/
      data revision /
     +     "$Revision: 6286 $"
     +      /
      data date /
     +     "$Date: 2013-04-05 23:36:00 +0200 (Fri, 05 Apr 2013) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
