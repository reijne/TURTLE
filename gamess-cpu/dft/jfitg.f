c---- memory counting routines -----------------------------------------
c
c  Coulomb fit gradient
c
      subroutine memreq_jfitg(memory_fp, memory_int)
      implicit none
c
c arguments
c      
      real*8 memory_fp(*)         ! Core memory
      integer memory_int(*)     ! Core memory

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

C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
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
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

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
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
c     The elements in this common block are used with the incore
c     coulomb fit approach.
c
c     ocfit_sw     : is .true. if the incore approach is to be used.
c
c     icfit_pt     : is the "pointer" to the fitting coefficient store.
c     ncfit        : is the size of the fitting coefficient store.
c
c     iite2c_stored: is the "pointer" to the 2-center Schwarz table
c     nte2c_shl    : is the size of the 2-center Schwarz table on each 
c                    node in integers.
c     iite3c_stored: is the "pointer" to the 3-center Schwarz table
c     nte3c_shl    : is the size of the 3-center Schwarz table on each
c                    node in integers.
c
c     ite2c_store  : is the "pointer" to the 2-center integral store
c     nte2c_int    : is the size of the 2-center integral store
c     ite3c_store  : is the "pointer" to the 3-center integral store
c     nte3c_int    : is the size of the 3-center integral store
c
c     idunlap_called: counts the number of times the dunlap fitting
c                     routines have been called since the last geometry
c                     change.
c
      integer icfit_pt, ncfit
      integer iite2c_stored, iite3c_stored, nte2c_shl, nte3c_shl
      integer ite2c_store, ite3c_store, nte2c_int, nte3c_int
      integer idunlap_called
      logical ocfit_sw
      common/dft_dunlap/icfit_pt, ncfit,
     +          iite2c_stored, iite3c_stored, nte2c_shl, nte3c_shl,
     +          ite2c_store, ite3c_store, nte2c_int, nte3c_int,
     +          idunlap_called,
     +          ocfit_sw

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

     
C Scratch space and pointers
      integer iso_pt,tr_pt,nr_pt, schwarz_ao, schwarz_cd, nsh
      integer gout_pt

      integer ntmp,n2c,n3c
      integer nnodes
c
c Local variables
c
      integer lbasf,i,j,col,ico
      integer ao_tag,cd_tag,cdbf_num
      integer ncentres,ao_basfn,tot_basfn
      integer max_shells,tot_prm
      integer size_required,size_freed
      real*8 iVn_sum,sum1,sum2,lambda,t_lambdan
      real*8 e_coul,e_self,rho
      integer basi,basj,bask,basl,imode
      real*8 dum(2)
      real*8 adens,bdens,grad
      integer iiso
      character*8 zrhf
      integer vqr_pt
c
c Functions
c
      integer incr_memory, lenwrd, null_memory
      integer ipg_nnodes
      real*8 ddot

      integer matrix_pt
      logical opg_root

      integer idum, cd_4c2eon, cd_jfitgon

c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)

      if(debug_sw .and. opg_root()) write(6,*) 'Entering Jfit_dunlap.f',
     &     num_bset

      potential_sw = .false.
      dunlap_sw = .true.
      e_coul   = 0.0d0
      rho      = 0.0d0

      ao_tag=bset_tags(1)
      cd_tag=bset_tags(2)

      cdbf_num=totbfn(cd_tag)
      ao_basfn=totbfn(ao_tag)
c
c     Storage for Schwarz integrals
c
      if (schwarz_tol.ge.0.and.
     &    ((idunlap_called.eq.0.and.ocfit_sw).or.
     &     (.not.ocfit_sw))) then
         schwarz_ao = incr_memory((nshell(ao_tag)+1)*
     &                             nshell(ao_tag)/2,'d')
         schwarz_cd = incr_memory(nshell(cd_tag),'d')
      else
         schwarz_ao = null_memory()
         schwarz_cd = null_memory()
      endif
c
c     Fit coefs
c
      if (ocfit_sw) then
c        The fitting coefficients are stored in the memory arranged by
c        CD_jfit_init1.
      else
         icfit_pt = incr_memory(cdbf_num,'d')
      endif
      tot_basfn=cdbf_num*ao_basfn
c
c     If we did not keep the fitting coefficients from the SCF
c     or if we restarted the calculation at the gradient evaluation
c     then we need to reconstruct the fitting coefficients. 
c     Otherwise we skip to the gradient evaluation straight away.
c
      if (.not.ocfit_sw.or.idunlap_called.eq.0) then
c
c        Tr
c
         tr_pt   = incr_memory(cdbf_num,'d')
c
c        Nr
c
         nr_pt   = incr_memory(cdbf_num,'d')
c
c        Vqr (2c2e integrals)
c
         size_required=cdbf_num*cdbf_num
         vqr_pt  = incr_memory(size_required,'d')

C *End memory allocation
C **********************************************************************
C *
C *Obtain fitting coefficients for the density
C *
C *
C *Calculate 2 centre 2 electron repulsion integrals and invert matrix
C *

         call memreq_Jfit_iVform(memory_int,
     &                    memory_fp,
     &                    cdbf_num,
     &                    .true.,
     &                    memory_fp(vqr_pt),
     &                    memory_fp(schwarz_cd))

C *
C *Form the tr matrix
C *
c        call aclear_dp(memory_fp(tr_pt),cdbf_num,0.0d0)

         nsh    = BL_max_shell_count()
         iso_pt = incr_memory(nsh,'d')

         gout_pt=incr_memory(50625,'d')
c        call aclear_dp(memory_fp(gout_pt),50625,0.0d0)
         
         if (schwarz_tol.ge.0.and.
     &      ((idunlap_called.eq.0.and.ocfit_sw).or.
     &       (.not.ocfit_sw))) then
c           call coulmb(memory_fp(schwarz_ao),memory_fp(gout_pt))
c           call te2c_rep_schwarz(cd_tag,memory_fp(gout_pt),
c    &                         memory_fp(schwarz_cd))
         endif

         basi = ao_tag
         basj = ao_tag
         bask = cd_tag
         basl = -1

         imode = 4

c        if (.not.ocfit_sw) then
c           call jkint_dft(memory_fp(iso_pt),
c    &           memory_fp(gout_pt),nsh,
c    &           basi, basj, bask, basl,
c    &           imode,adens,bdens,memory_fp(tr_pt),dum,
c    &           memory_fp(schwarz_ao),
c    &           memory_fp(schwarz_cd))
c        else if (idunlap_called.eq.0) then
c           call jkint_dft_cntgen(memory_fp(iso_pt),
c    &           memory_fp(gout_pt),nsh,
c    &           basi, basj, bask, basl,
c    &           imode,adens,bdens,memory_fp(tr_pt),dum,
c    &           memory_fp(schwarz_ao),
c    &           memory_fp(schwarz_cd),
c    &           memory_fp(ite3c_store),nte3c_int,
c    &           memory_int(iite3c_stored),nte3c_shl,
c    &           memory_fp(ite2c_store),nte2c_int,
c    &           memory_int(iite2c_stored),nte2c_shl)
c        else
c           call caserr('jfitg: logically impossible to get here!')
c        endif

         call decr_memory(gout_pt,'d')
         call decr_memory(iso_pt,'d')

c        idunlap_called = idunlap_called + 1


C *
C *Form 1 centre integral array from fitting functions, ie nr
C *
c        call aclear_dp(memory_fp(nr_pt),cdbf_num-1,0.0d0)

         call memreq_intPack_drv(memory_int,
     &                    memory_fp, 
     &                    memory_fp(icfit_pt),
     &                    .true.,.false.,
     &                    .true.,
     &                    .false.,.false.,
     &                    .false.,.false.,
     &                    .false.,
     &                    memory_fp(nr_pt)) 
C *
C *Calculate Langrange Multiplier lambda
C *
         col  = 0
         sum1 = 0.0d0
         sum2 = 0.0d0

c        do lbasf=0,cdbf_num-1
c           iVn_sum=ddot(cdbf_num,memory_fp(vqr_pt+col)
c    &      ,1,memory_fp(nr_pt),1)
            col=col+cdbf_num
c           sum1=sum1+iVn_sum*memory_fp(tr_pt+lbasf)
c           sum2=sum2+iVn_sum*memory_fp(nr_pt+lbasf)
c        enddo
c        lambda=(nelectrons-sum1)/sum2

C *
C *api entry point
C
c        lmult=lambda

C *
C *Calculate fitting coefficients themselves
C *
c        call aclear_dp(memory_fp(icfit_pt),cdbf_num,0.0d0)
         col = 0
c        do lbasf=0,cdbf_num-1
c           t_lambdan=memory_fp(tr_pt+lbasf)+
c    &                lambda*memory_fp(nr_pt+lbasf)
c           call daxpy(cdbf_num,t_lambdan
c    &           ,memory_fp(vqr_pt+col),1
c    &           ,memory_fp(icfit_pt),1)
            col=col+cdbf_num
c        enddo

C 
C Calculate matrix elements and energy
c 
         call decr_memory(vqr_pt,'d')
         call decr_memory(nr_pt,'d')
         call decr_memory(tr_pt,'d')
      endif
c
c NB should be redundant as symmetry is not used, but array is
c still referenced
c
      nsh  = max(nw196(5),BL_max_shell_count())
      iiso = incr_memory(nsh,'d')
c     call rdedx(memory_fp(iiso),nw196(5),ibl196(5),idaf)
      zrhf = 'rhf'
c
c This is the 4-centre force evaluation used for
c checking
c
c      call jkder_dft(zrhf,memory_fp,memory_fp(iiso),nshell(1),
c     &     ao_tag, ao_tag, ao_tag, ao_tag,
c     &     adens,bdens,memory_fp(icfit_pt),dum,grad)
c

c
c 2 centre term
c
      if (.not.ocfit_sw) then
         call memreq_jkder_dft(zrhf,memory_fp,memory_fp(iiso),nsh,
     &        cd_tag, -1, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_fp(schwarz_ao),memory_fp(schwarz_cd))
      else if (idunlap_called.eq.0) then
         call caserr('Jfitg needs jfit to be run first')
      else 
         call memreq_jkder_dft_genuse(zrhf,memory_fp,memory_fp(iiso),
     &        nsh, cd_tag, -1, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_int(iite3c_stored),nte3c_shl,
     &        memory_int(iite2c_stored),nte2c_shl)
      endif
c
c 3 centre term
c
      if (.not.ocfit_sw) then
         call memreq_jkder_dft(zrhf,memory_fp,memory_fp(iiso),nsh,
     &        ao_tag, ao_tag, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_fp(schwarz_ao),memory_fp(schwarz_cd))
      else if (idunlap_called.eq.0) then
         call caserr('Jfitg needs jfit to be run first')
      else 
         call memreq_jkder_dft_genuse(zrhf,memory_fp,memory_fp(iiso),
     &        nsh, ao_tag, ao_tag, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_int(iite3c_stored),nte3c_shl,
     &        memory_int(iite2c_stored),nte2c_shl)
      endif

      call decr_memory(iiso,'d')


      if (ocfit_sw) then
c        Fitting coefficient memory taken care of by CD_jfit_clean1
      else
         call decr_memory(icfit_pt,'d')
      endif
      if (schwarz_tol.ge.0.and.
     &    ((idunlap_called.eq.1.and.ocfit_sw).or.
     &     (.not.ocfit_sw))) then
         call decr_memory(schwarz_cd,'d')
         call decr_memory(schwarz_ao,'d')
      endif

      return

      end
c
c---- the routines that do the real work -------------------------------
c
c  Coulomb fit gradient
c
      subroutine jfitg(memory_fp, memory_int, adens, bdens, grad, iout)
      implicit none
c
c arguments
c      
      integer iout              ! Unit for standard out.
      real*8 memory_fp(*)         ! Core memory
      integer memory_int(*)     ! Core memory
      real*8 adens(*), bdens(*)   ! Density matrix     (in)
      real*8 grad(3,*)            ! Accumulated forces (in/out)

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

C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
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
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

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
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
c     The elements in this common block are used with the incore
c     coulomb fit approach.
c
c     ocfit_sw     : is .true. if the incore approach is to be used.
c
c     icfit_pt     : is the "pointer" to the fitting coefficient store.
c     ncfit        : is the size of the fitting coefficient store.
c
c     iite2c_stored: is the "pointer" to the 2-center Schwarz table
c     nte2c_shl    : is the size of the 2-center Schwarz table on each 
c                    node in integers.
c     iite3c_stored: is the "pointer" to the 3-center Schwarz table
c     nte3c_shl    : is the size of the 3-center Schwarz table on each
c                    node in integers.
c
c     ite2c_store  : is the "pointer" to the 2-center integral store
c     nte2c_int    : is the size of the 2-center integral store
c     ite3c_store  : is the "pointer" to the 3-center integral store
c     nte3c_int    : is the size of the 3-center integral store
c
c     idunlap_called: counts the number of times the dunlap fitting
c                     routines have been called since the last geometry
c                     change.
c
      integer icfit_pt, ncfit
      integer iite2c_stored, iite3c_stored, nte2c_shl, nte3c_shl
      integer ite2c_store, ite3c_store, nte2c_int, nte3c_int
      integer idunlap_called
      logical ocfit_sw
      common/dft_dunlap/icfit_pt, ncfit,
     +          iite2c_stored, iite3c_stored, nte2c_shl, nte3c_shl,
     +          ite2c_store, ite3c_store, nte2c_int, nte3c_int,
     +          idunlap_called,
     +          ocfit_sw

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

     
C Scratch space and pointers
      integer iso_pt,tr_pt,nr_pt, schwarz_ao, schwarz_cd, nsh
      integer gout_pt

      integer ntmp,n2c,n3c
      integer nnodes
c
c Local variables
c
      integer lbasf,i,j,col,ico
      integer ao_tag,cd_tag,cdbf_num
      integer ncentres,ao_basfn,tot_basfn
      integer max_shells,tot_prm
      integer size_required,size_freed
      real*8 iVn_sum,sum1,sum2,lambda,t_lambdan
      real*8 e_coul,e_self,rho
      integer basi,basj,bask,basl,imode
      real*8 dum(2)
      integer iiso
      character*8 zrhf
c
      integer push_memory_count,    pop_memory_count
      integer push_memory_estimate, pop_memory_estimate
      integer imemcount, imemusage, imemestimate
c
      integer vqr_pt
c
c Functions
c
      integer allocate_memory, lenwrd, null_memory
      integer ipg_nnodes
      real*8 ddot

      integer matrix_pt
      logical opg_root

      integer idum, cd_4c2eon, cd_jfitgon

      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      imemcount = push_memory_count()

      if(debug_sw .and. opg_root()) write(iout,*) 
     &     'Entering Jfit_dunlap.f',num_bset

      call start_time_period(TP_DFT_JFITG)

      potential_sw = .false.
      dunlap_sw = .true.
      e_coul   = 0.0d0
      rho      = 0.0d0

      ao_tag=bset_tags(1)
      cd_tag=bset_tags(2)

      cdbf_num=totbfn(cd_tag)
      ao_basfn=totbfn(ao_tag)
c
c     Storage for Schwarz integrals
c
      if (schwarz_tol.ge.0.and.
     &    ((idunlap_called.eq.0.and.ocfit_sw).or.
     &     (.not.ocfit_sw))) then
         schwarz_ao = allocate_memory((nshell(ao_tag)+1)*
     &                                 nshell(ao_tag)/2,'d')
         schwarz_cd = allocate_memory(nshell(cd_tag),'d')
      else
         schwarz_ao = null_memory()
         schwarz_cd = null_memory()
      endif
c
c     Fit coefs
c
      if (ocfit_sw) then
c        The fitting coefficients are stored in the memory arranged by
c        CD_jfit_init1.
      else
         icfit_pt = allocate_memory(cdbf_num,'d')
      endif
      tot_basfn=cdbf_num*ao_basfn
c
c     If we did not keep the fitting coefficients from the SCF
c     or if we restarted the calculation at the gradient evaluation
c     then we need to reconstruct the fitting coefficients. 
c     Otherwise we skip to the gradient evaluation straight away.
c
      if (.not.ocfit_sw.or.idunlap_called.eq.0) then
c
c        Tr
c
         tr_pt   = allocate_memory(cdbf_num,'d')
c
c        Nr
c
         nr_pt   = allocate_memory(cdbf_num,'d')
c
c        Vqr (2c2e integrals)
c
         size_required=cdbf_num*cdbf_num
         vqr_pt  = allocate_memory(size_required,'d')

C *End memory allocation
C **********************************************************************
C *
C *Obtain fitting coefficients for the density
C *
         if(debug_sw .and. opg_root()) then
           write(iout,*)
           if(triangle_sw) write(iout,*) 'Using triangle matrices'
           write(iout,*) 'Number of fitting basis functions :',cdbf_num
           write(iout,*) 'AO matrix size                    :',ao_basfn
           write(iout,*) 'Matrix size                       :',tot_basfn
         endif
C *
C *Calculate 2 centre 2 electron repulsion integrals and invert matrix
C *

         call start_time_period(TP_DFT_JFITG_VFORM)


         call Jfit_iVform(memory_int,
     &                    memory_fp,
     &                    cdbf_num,
     &                    .true.,
     &                    memory_fp(vqr_pt),
     &                    memory_fp(schwarz_cd))
         call end_time_period(TP_DFT_JFITG_VFORM)

         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(6,*)'check guards after ivform'
            call gmem_check_guards
         endif

C *
C *Form the tr matrix
C *
         call aclear_dp(memory_fp(tr_pt),cdbf_num,0.0d0)

         call start_time_period(TP_DFT_JFITG_TR)

         nsh    = BL_max_shell_count()
         iso_pt = allocate_memory(nsh,'d')

         gout_pt=allocate_memory(50625,'d')
         call aclear_dp(memory_fp(gout_pt),50625,0.0d0)
         
         if (schwarz_tol.ge.0.and.
     &      ((idunlap_called.eq.0.and.ocfit_sw).or.
     &       (.not.ocfit_sw))) then
            call coulmb_dft(memory_fp(schwarz_ao),memory_fp(gout_pt),
     &                      ao_tag)
            call te2c_rep_schwarz(cd_tag,memory_fp(gout_pt),
     &                         memory_fp(schwarz_cd))
            if(print_sw(DEBUG_MEMORY))then
               if(opg_root())write(6,*)'check guards after coulmb'
               call gmem_check_guards
            endif
         endif

         basi = ao_tag
         basj = ao_tag
         bask = cd_tag
         basl = -1

         imode = 4

         if (.not.ocfit_sw) then
            call jkint_dft(memory_fp(iso_pt),
     &           memory_fp(gout_pt),nsh,
     &           basi, basj, bask, basl,
     &           imode,adens,bdens,memory_fp(tr_pt),dum,
     &           memory_fp(schwarz_ao),
     &           memory_fp(schwarz_cd))
         else if (idunlap_called.eq.0) then
            call jkint_dft_cntgen(memory_fp(iso_pt),
     &           memory_fp(gout_pt),nsh,
     &           basi, basj, bask, basl,
     &           imode,adens,bdens,memory_fp(tr_pt),dum,
     &           memory_fp(schwarz_ao),
     &           memory_fp(schwarz_cd),
     &           memory_fp(ite3c_store),nte3c_int,
     &           memory_int(iite3c_stored),nte3c_shl,
     &           memory_fp(ite2c_store),nte2c_int,
     &           memory_int(iite2c_stored),nte2c_shl)
         else
            call caserr('jfitg: logically impossible to get here!')
         endif

         call free_memory(gout_pt,'d')
         call free_memory(iso_pt,'d')

         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(iout,*)'check guards after tr'
            call gmem_check_guards
         endif

         if(opg_root() .and. print_sw(DEBUG_TR))then
            do ico=0,cdbf_num-1
               write(iout,*) 'Tr Vector:',ico+1,memory_fp(tr_pt+ico)
            enddo
         endif

         idunlap_called = idunlap_called + 1

         call end_time_period(TP_DFT_JFITG_TR)

C *
C *Form 1 centre integral array from fitting functions, ie nr
C *
         call aclear_dp(memory_fp(nr_pt),cdbf_num-1,0.0d0)

         call start_time_period(TP_DFT_JFITG_NR)

         call intPack_drv(memory_int,
     &                    memory_fp, 
     &                    memory_fp(icfit_pt),
     &                    .true.,.false.,
     &                    .true.,
     &                    .false.,.false.,
     &                    .false.,.false.,
     &                    .false.,
     &                    memory_fp(nr_pt)) 

         call end_time_period(TP_DFT_JFITG_NR)


         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(iout,*)'check guards after nr'
            call gmem_check_guards
         endif

         if(opg_root() .and. print_sw(DEBUG_NR))then
            do ico=0,cdbf_num-1
               write(iout,*) 'nr:',ico+1,memory_fp(nr_pt+ico)
            enddo
         endif
C *
C *Calculate Langrange Multiplier lambda
C *
         col  = 0
         sum1 = 0.0d0
         sum2 = 0.0d0

         call start_time_period(TP_DFT_JFITG_COEF)
         do lbasf=0,cdbf_num-1
            iVn_sum=ddot(cdbf_num,memory_fp(vqr_pt+col)
     &      ,1,memory_fp(nr_pt),1)
            col=col+cdbf_num
            sum1=sum1+iVn_sum*memory_fp(tr_pt+lbasf)
            sum2=sum2+iVn_sum*memory_fp(nr_pt+lbasf)
c           write(iout,*) 'iVn_sum:',iVn_sum,'SUMS',sum1,sum2
         enddo
         lambda=(nelectrons-sum1)/sum2

C *
C *api entry point
C
         lmult=lambda

         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(iout,*)'check guards after lambda'
            call gmem_check_guards
         endif
C *
C *Calculate fitting coefficients themselves
C *
         call aclear_dp(memory_fp(icfit_pt),cdbf_num,0.0d0)
         col = 0
         do lbasf=0,cdbf_num-1
            t_lambdan=memory_fp(tr_pt+lbasf)+
     &                lambda*memory_fp(nr_pt+lbasf)
            call daxpy(cdbf_num,t_lambdan
     &           ,memory_fp(vqr_pt+col),1
     &           ,memory_fp(icfit_pt),1)
            col=col+cdbf_num
         enddo

         call end_time_period(TP_DFT_JFITG_COEF)

         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(iout,*)'check guards after fit'
            call gmem_check_guards
         endif
C 
C Calculate matrix elements and energy
c 
         if(opg_root() .and. print_sw(DEBUG_JFIT))then
            write(iout,*)'Cfit:'
            do col=0,cdbf_num-1
               write(iout,101)memory_fp(icfit_pt+col)
 101           format(1x,f16.5)
            enddo
         endif

         if(print_sw(DEBUG_MEMORY))then
            if(opg_root())write(iout,*)'check guards after energy'
            call gmem_check_guards
         endif

         call free_memory(vqr_pt,'d')
         call free_memory(nr_pt,'d')
         call free_memory(tr_pt,'d')
      endif
c
c NB should be redundant as symmetry is not used, but array is
c still referenced
c
      nsh  = max(nw196(5),BL_max_shell_count())
      iiso = allocate_memory(nsh,'d')
      call rdedx(memory_fp(iiso),nw196(5),ibl196(5),idaf)
      zrhf = 'rhf'
c
c This is the 4-centre force evaluation used for
c checking
c
c      call jkder_dft(zrhf,memory_fp,memory_fp(iiso),nshell(1),
c     &     ao_tag, ao_tag, ao_tag, ao_tag,
c     &     adens,bdens,memory_fp(icfit_pt),dum,grad)
c

      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
         write(iout,*) 'Forces before coulomb:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad(j,i),j=1,3)
         enddo
 1000    format(i3,3f16.8)
      endif
c
c 2 centre term
c
      call start_time_period(TP_DFT_JFITG_2C)
      if (.not.ocfit_sw) then
         call jkder_dft(zrhf,memory_fp,memory_fp(iiso),nsh,
     &        cd_tag, -1, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_fp(schwarz_ao),memory_fp(schwarz_cd))
      else if (idunlap_called.eq.0) then
         call caserr('Jfitg needs jfit to be run first')
      else 
         call jkder_dft_genuse(zrhf,memory_fp,memory_fp(iiso),
     &        nsh, cd_tag, -1, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_int(iite3c_stored),nte3c_shl,
     &        memory_int(iite2c_stored),nte2c_shl)
      endif

      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
         write(iout,*) 'Forces after 2c:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad(j,i),j=1,3)
         enddo
      endif
      call end_time_period(TP_DFT_JFITG_2C)
c
c 3 centre term
c
      call start_time_period(TP_DFT_JFITG_3C)
      if (.not.ocfit_sw) then
         call jkder_dft(zrhf,memory_fp,memory_fp(iiso),nsh,
     &        ao_tag, ao_tag, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_fp(schwarz_ao),memory_fp(schwarz_cd))
      else if (idunlap_called.eq.0) then
         call caserr('Jfitg needs jfit to be run first')
      else 
         call jkder_dft_genuse(zrhf,memory_fp,memory_fp(iiso),
     &        nsh, ao_tag, ao_tag, cd_tag, -1,
     &        adens,bdens,memory_fp(icfit_pt),dum,grad,
     &        memory_int(iite3c_stored),nte3c_shl,
     &        memory_int(iite2c_stored),nte2c_shl)
      endif

      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
         write(iout,*) 'Forces after 3c:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad(j,i),j=1,3)
         enddo
      endif

      call free_memory(iiso,'d')
      call end_time_period(TP_DFT_JFITG_3C)


      call end_time_period(TP_DFT_JFITG)

      if (ocfit_sw) then
c        Fitting coefficient memory taken care of by CD_jfit_clean1
      else
         call free_memory(icfit_pt,'d')
      endif
      if (schwarz_tol.ge.0.and.
     &    ((idunlap_called.eq.1.and.ocfit_sw).or.
     &     (.not.ocfit_sw))) then
         call free_memory(schwarz_cd,'d')
         call free_memory(schwarz_ao,'d')
      endif

      imemusage = pop_memory_count(imemcount)
      imemcount = push_memory_estimate()
      call memreq_jfitg(memory_fp,memory_int)
      imemestimate = pop_memory_estimate(imemcount)
      if (opg_root().and.imemusage.ne.imemestimate) then
         write(iout,*)'*** estimated memory usage = ',imemestimate,
     &                ' words'
         write(iout,*)'*** actual    memory usage = ',imemusage   ,
     &                ' words'
         write(iout,*)'*** WARNING: the memory usage estimates for ',
     &                'jfitg seem to be incorrect'
      endif

      return

      end

      subroutine ver_dft_jfitg(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/jfitg.m,v $
     +     "/
      data revision /
     +     "$Revision: 5919 $"
     +      /
      data date /
     +     "$Date: 2009-03-31 15:11:24 +0200 (Tue, 31 Mar 2009) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
