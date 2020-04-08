c
c  implementation of DFT memory routines using GAMESS-UK
c  memory manager
c
      integer function push_memory_estimate()
      implicit none
      integer igmem_push_estimate
      push_memory_estimate = igmem_push_estimate()
      return
      end

      integer function pop_memory_estimate(icnt)
      implicit none
      integer icnt
      integer igmem_pop_estimate
      pop_memory_estimate = igmem_pop_estimate(icnt)
      return
      end

      integer function push_memory_count()
      implicit none
      integer igmem_push_usage
      push_memory_count = igmem_push_usage()
      return
      end

      integer function pop_memory_count(icnt)
      implicit none
      integer icnt
      integer igmem_pop_usage
      pop_memory_count = igmem_pop_usage(icnt)
      return
      end

c     subroutine init_memory_count
c     implicit none
c     call gmem_count_init
c     return
c     end

c     subroutine reset_memory_count
c     implicit none
c     call gmem_count_reset
c     return
c     end

      integer function max_memory_count()
      implicit none
      integer igmem_count_max
      max_memory_count = igmem_count_max()
      return
      end

      integer function current_memory_count()
      implicit none
      integer igmem_count_current
      current_memory_count = igmem_count_current()
      return
      end

      integer function incr_memory(length,type)
      implicit none
      integer iwr
      common/dft_iofile/iwr
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_incr
      integer null_memory
      character*(*) type
      if(type .eq. 'd') then
c
         itemp = igmem_incr(length)
         incr_memory = itemp
         return 
c
      else if (type .eq. 'i') then
c 
         ilen = lenwrd()
c        itemp = igmem_incr((length-1)/ilen + 1)
c        the line above and line below are not the same in fortran77.
         incr_memory = igmem_incr((length+ilen-1)/ilen)
         return 
c
      else
         write(iwr,*)type
         incr_memory = -1
         call caserr('incr_memory: unknown data type')
         incr_memory = null_memory()
      endif
c
      end
      

      subroutine decr_memory(address,type)
      integer iwr
      common/dft_iofile/iwr
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory
c
      if(type .eq. 'd') then
c
         call gmem_decr(address)
c
      else if (type .eq. 'i') then
c
         ilen = lenwrd()
c        call gmem_decr( (address + (ilen-1)) /ilen)
         call gmem_decr(address)

      else
         write(6,*)type
         call caserr('decr_memory: unknown type')
      endif
c
      address = null_memory()
c
      end


      integer function incr_memory2(length,type,filename,subrname,
     &                              varid)
      implicit none
      integer length
      integer incr_memory
      character*(*) filename,subrname,varid
      character*(*) type
      incr_memory2=incr_memory(length,type)
      return
      end


      subroutine decr_memory2(address,type,filename,subrname,varid)
      implicit none
      integer address
      character*(*) type
      character*(*) filename,subrname,varid
      call decr_memory(address,type)
      end


      integer function memory_overhead()
      implicit none
      integer igmem_overhead
      memory_overhead = igmem_overhead()
      return
      end

      integer function null_memory()
      implicit none
      integer igmem_null
      null_memory = igmem_null()
      return
      end

      integer function allocate_memory(length,type)
      implicit none
      integer iwr
      common/dft_iofile/iwr
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_alloc
      integer null_memory

      character*(*) type

      if(type .eq. 'd') then

         itemp = igmem_alloc(length)
         allocate_memory = itemp
         return 

      else if (type .eq. 'i') then

         ilen = lenwrd()
c        itemp = igmem_alloc((length-1)/ilen + 1)
c        the line above and line below are not the same in fortran77.
         itemp = igmem_alloc((length+ilen-1)/ilen)
         allocate_memory = itemp * ilen - (ilen -1)
         return 

      else
         write(6,*)type
         allocate_memory = -1
         call caserr('allocate_memory: unknown data type')
         allocate_memory = null_memory()
      endif

      end
      
      subroutine free_memory(address,type)
      integer iwr
      common/dft_iofile/iwr
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory

      if(type .eq. 'd') then

         call gmem_free(address)

      else if (type .eq. 'i') then

         ilen = lenwrd()
         call gmem_free( (address + (ilen-1)) /ilen)

      else
         write(6,*)type
         call caserr('free_memory: unknown type')
      endif

      address = null_memory()

      end

      integer function allocate_memory2(length,type,filename,subrname,
     &                                  varid)
      implicit none
      integer iwr
      common/dft_iofile/iwr
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer length
      integer lenwrd
      integer itemp
      integer ilen
      integer igmem_alloc_inf
      integer null_memory
      character*(*) filename,subrname,varid

      character*(*) type

      if(type .eq. 'd') then

         itemp = igmem_alloc_inf(length,filename,subrname,
     &                           varid,IGMEM_DEBUG)
         allocate_memory2 = itemp
         return 

      else if (type .eq. 'i') then

         ilen = lenwrd()
c        itemp = igmem_alloc_inf((length-1)/ilen + 1,filename,subrname,
c    &                           varid,IGMEM_DEBUG)
c        the line above and line below are not the same in fortran77.
         itemp = igmem_alloc_inf((length+ilen-1)/ilen,filename,subrname,
     &                           varid,IGMEM_DEBUG)
         allocate_memory2 = itemp * ilen - (ilen -1)
         return 

      else
         write(6,*)type
         allocate_memory2 = -1
         call caserr('allocate_memory: unknown data type')
         allocate_memory2 = null_memory()
      endif

      end
      
      subroutine free_memory2(address,type,filename,subrname,varid)
      implicit none
      integer iwr
      common/dft_iofile/iwr
      integer ilen
      integer address
      character*(*) type
      integer lenwrd
      integer null_memory
      character*(*) filename,subrname,varid

      if(type .eq. 'd') then

         call gmem_free_inf(address,filename,subrname,varid)

      else if (type .eq. 'i') then

         ilen = lenwrd()
         call gmem_free_inf( (address + (ilen-1)) /ilen,filename,
     &                      subrname,varid)

      else
         write(6,*)type
         call caserr('free_memory: unknown type')
      endif

      address = null_memory()

      end

      subroutine interface_gamess(rhf_sw,opshell_sw,ao_tag,iout)
C **********************************************************************
C *Description:							       *
C *Interface into gamess for ccp1 dft modules                          *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations
C *
C *Parameters
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
C *In variables
      logical rhf_sw
      integer ao_tag
      integer iout
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
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
C *Out variables
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
C *System information
      integer natoms,nelectrons,ian,ngridcentres
      real*8  atom_c
      common/sysinf/atom_c(max_atom,3),ian(max_atom),natoms,nelectrons,
     +              ngridcentres
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
      logical opshell_sw
C * Functions
      logical     opg_root
      integer     isubst
C *Local variables
      integer     necheck
      integer     atmnum(maxat)
      integer     latm,lgshl,lshl,kgshl,kshl
      logical     atomexist_sw, osame
      integer     atom_id(maxat),minshl,maxshl,ploc
c     character*8 z_at_id(max_atype)
c     integer     n_at_id,i_id, i, j
      integer     i_id, i, j
      integer     hybri
      integer     latmno, katmno
      integer     nlshl, nkshl, plloc, pkloc, nfunc
      integer     ierror
C *End declarations                                                    *
C **********************************************************************
c     if(opg_root())call ccpdft_banner(iout)
      necheck=mod(ne,2)
      if(necheck.ne.0.and.rhf_sw) opshell_sw = .true.

      ierror = BL_clear_basis_set(ao_tag)

      do i = 1,natoms
        atom_id(i) = 0
      enddo
C 
C Export AO basis set
c     n_at_id = 0
      do lgshl=1,nshell
        atomexist_sw=.false.
        latmno = isubst(zaname(katom(lgshl)))
        do kgshl=1,lgshl-1
           if(zaname(katom(kgshl)).eq.zaname(katom(lgshl))) then
              atomexist_sw = .true.
              atom_id(katom(lgshl)) = atom_id(katom(kgshl))
              goto 5
           endif
        enddo
 5      continue
c
        if(.not.atomexist_sw) then
c
c...       find all shells on centre l
c
           nlshl=lgshl
           do lshl=lgshl+1,nshell
             if(katom(lshl).eq.katom(lgshl)) nlshl=nlshl+1
           enddo
c
           if (lgshl.eq.1) then
              atomexist_sw = .false.
           else
              kgshl=1
c
c...          find all shells on centre k
c
 10           nkshl=kgshl
              do kshl=kgshl+1,nshell
                 if(katom(kshl).eq.katom(kgshl)) nkshl=nkshl+1
              enddo
c
c...          check if all exponents and coefficients are the same on 
c...          both centres (if the nuclear charges are different we
c...          will assume that the basis is different too).
c
              katmno = isubst(zaname(katom(kgshl)))
              if (nkshl-kgshl.eq.nlshl-lgshl.and.latmno.eq.katmno) then
                 osame = .true.
                 do i = 0,nkshl-kgshl
c
c...                check if the number of primitives is the same
c
                    osame = osame.and.(kng(lgshl+i).eq.kng(kgshl+i))
     +                           .and.(kmin(lgshl+i).eq.kmin(kgshl+i))
     +                           .and.(kmax(lgshl+i).eq.kmax(kgshl+i))
                    if (osame) then
                       nfunc = kng(lgshl+i)
                       plloc = kstart(lgshl+i)
                       pkloc = kstart(kgshl+i)
                       do j = 0, nfunc-1
                          osame = osame.and.(ex(plloc+j).eq.ex(pkloc+j))
                          osame = osame.and.(
     +                               cs(plloc+j).eq.cs(pkloc+j).or.
     +                               1.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cp(plloc+j).eq.cp(pkloc+j).or.
     +                               2.gt.kmax(lgshl+i).or.
     +                               4.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cd(plloc+j).eq.cd(pkloc+j).or.
     +                               5.gt.kmax(lgshl+i).or.
     +                               10.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cf(plloc+j).eq.cf(pkloc+j).or.
     +                               11.gt.kmax(lgshl+i).or.
     +                               20.lt.kmin(lgshl+i) )
                          osame = osame.and.(
     +                               cg(plloc+j).eq.cg(pkloc+j).or.
     +                               21.gt.kmax(lgshl+i))
                       enddo
                    endif
                 enddo
              else
                 osame = .false.
              endif
              atomexist_sw = atomexist_sw.or.osame
              if (.not.atomexist_sw) then
                 kgshl=nkshl+1
                 if (kgshl.lt.lgshl) goto 10
              else
                 atom_id(katom(lgshl)) = atom_id(katom(kgshl))
              endif
           endif
c
c...       
           if(.not.atomexist_sw) then
c
c step 1 - create new atom tag
c            if (n_at_id.ge.max_atype) then
c               call caserr('Maximum number of atom types exceeded')
c            endif
c            n_at_id=n_at_id+1
             atom_id(katom(lgshl))=BL_create_atomtag(ao_tag,
     +                                               ian(katom(lgshl)))
c            z_at_id(n_at_id)=zaname(katom(kgshl))
             if (atom_id(katom(lgshl)).lt.0) then
                call caserr('Failed to create a new atom tag')
             endif
c
c step 2 - find out how many shells exist on this centre
c            nshl=lgshl
c            do lshl=lgshl+1,nshell
c               if(katom(lshl).eq.katom(lgshl)) nshl=nshl+1
c            enddo
c
c step3 - loop over shells on this centre and import the information
             minshl = lgshl
             maxshl = nlshl
             do lshl=minshl,maxshl
               ploc=kstart(lshl)
               hybri=ktype(lshl)
               if((ktype(lshl).eq.2).and.(kmax(lshl)-kmin(lshl).eq.3))
     &           hybri=1
               ierror=BL_import_shell(ao_tag,atom_id(katom(lshl)),
     &                              kng(lshl),ktype(lshl),hybri,
     &                              ex(ploc),cs(ploc),cp(ploc),cd(ploc),
     &                              cf(ploc),cg(ploc))
             enddo
           endif
        endif
      enddo
c
c     Assign basis functions to centres
c
c     do lgshl=1,nshell
c       do i_id=1,n_at_id
c          if (z_at_id(i_id).eq.zaname(katom(lgshl))) then
c             ierror = BL_assign_type(ao_tag,katom(lgshl),atom_id(i_id))
c             if(ierror .ne. 0) call caserr(
c    &             'basis set assignment failed for a.o. basis')
c          endif
c       enddo

      do i = 1, natoms
        ierror = BL_assign_type(ao_tag,i,atom_id(i))
        if(ierror .ne. 0) call caserr(
     &             'basis set assignment failed for a.o. basis')
      enddo
C 
C Assign basis functions to centres
c Currently only the atomic number is used
c
c     ierror=BL_assign_types_by_z(ao_tag)
c     if(ierror .ne. 0) call caserr(
c    &     'basis set assignment failed for a.o. basis')
C
C Check in basis
      call checkin_basis(ao_tag,.false.)
      bset_tags(ao_tag) = 1

      return
      end
c
c  $Author: hvd $ $Revision: 5774 $ $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gamess.m,v $
c

      subroutine ver_dft_gamess(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gamess.m,v $
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
