      subroutine jkint_dft(iso,gout,nshels,
     &     basi, basj, bask, basl,
     &     imode,aaa,bbb,ccc,ddd, schwarz_ao, schwarz_cd)

      implicit none

      real*8 gout, aaa, bbb, ccc, ddd, schwarz_ao, schwarz_cd
      dimension gout(*), aaa(*), bbb(*), ccc(*), ddd(*)
      dimension schwarz_ao(*), schwarz_cd(*)
      integer basi, basj, bask, basl
      integer imode, nshels
      integer iso(nshels,*)

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
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)


cINCLUDE(../m4/common/symtry)
       integer nt
c
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
      integer iwr
      common/dft_iofile/iwr
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
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
cINCLUDE(common/dft_auxvar)
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
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
c inclusion of dft_indez needed for new inlined version of mat_load_3c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      logical odbg
      common/dbgdbg/odbg
c
      logical opg_root
      real*8 cpulft
c
c Local variables
c 
      integer i, j, k
      integer ld, ii, jj, kk, ll, j0,l0, k0
      integer id, jd, kd, nd
      integer ish, jsh, ksh, lsh
      integer ist0, jst0, kst0, lst0
      integer it, itrij, itrkl, ikykk, n4, mij, mjk
      integer ikyii, loop
      integer maxjj, maxkk, maxll
      integer icut0, nschwz
      integer ann, nn, knn, ltyi, ltyj, ltyk, imm, ijn

      integer imc, nijx, ibas_num, jbas_num
      integer locij_t, locij_km, locij_dm, idum
      real*8 fac

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      logical oskipp
      real*8 q4, tijl, pi, schwarz_lim, test, schw_max

      real*8 dt0, dt1, tim0, time, tim1

      integer mi, mj, mk, m0
      dimension mi(48),mj(48),mk(48),m0(48)

      integer ib(4,4)

      real*8 scr(225)  

      integer m25, itmp
      real*8 done, two, twopt5, four, e
      character*1 xn,xt
      integer nlen

      data  nt/1/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data m25/25/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data xn,xt/'n','t'/
      data e /2.30258d0/
c
c ***
c *** establish arrays normally set in main line code
c ***

c      dummy iso array
c      note that nt is hardwired to 1
      idum = 0
      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
c
c ... /mapperx/
c
      do i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
      enddo
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue

c
c ... /auxvarx/
c
c      ofast = .false.
c
cc      call setsto(10,0,intcut)
cc      call setsto(1060,0,intmag)
cc      nopk = 1

      icut0 = iabs(icut_3c)
      if (schwarz_tol.ge.0) then
         schwarz_lim = done/(10.0d0**schwarz_tol)
      else
         schwarz_lim = -1.0d0
      endif
      cutoff = done/(10.0d0**icut0)
      tol = e*itol_3c
c
c***
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_dft
      nopkr = nopk
      iofrst = iofsym
c
      nindmx = 1
      call aux_find(basi)
c
c     ----- two-electron integrals -----
c

c
c Determine which term to compute
c
      ncentr = 0

      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         write(iwr,*)'fitting coefficients',
     &        (ccc(i),i=1,BL_basis_size(bask))
      endif

      if(opg_root().and.(nprint.ne.-5))then
         write(iwr,12345)ncentr,imode,cpulft(1)
         write(iwr,12346)basi,basj,bask,basl
         write(iwr,12347)icut_3c,itol_3c
      endif

c-schw      ltri = ikyp(nshels)
c      l2 = ikyp(numorb)
      oskipp = .false.

c fail

      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
         igt(loop) = ib(1,ii)
         jgt(loop) = ib(2,ii)
         kgt(loop) = ib(3,ii)
         lgt(loop) = ib(4,ii)
 30   continue


c imode 6,61
      imc=0
c 3centr
      locij_dm=1
      locij_km=0

c
c     ----- set some parameters -----
c
c     ----- allocate core memory
c
      nschwz = 0
c
      ist0 = 1
      jst0 = 1
      kst0 = 1
      lst0 = 1

c
c     ----- ishell -----
c

      if (ist0.le.nshell(basi)) then

c@@ invert for parallel

         do 150 ii = ist0 , nshell(basi)
c
c     ----- eliminate ishell -----
c
            do 50 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) then
                  write(iwr,*)'symm skip on ii',id,ii,nt
                  go to 150
               endif
               m0(it) = id
 50         continue
            ikyii = iky(ii)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst0
               maxjj = ii
            endif

c@@ parallel invert

            do 140 jj = j0 , maxjj

               jst0 = 1
               itrij = ikyii + jj
c
c apply i/j tests only when i and j are AO basis fns
c @@ Note this still dangerous - see comment on triangulation
c in ijprim
c
               if (ncentr .eq. 3  .or. ncentr .eq. 4)then
               do 70 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) then
                     write(iwr,*)'symm skip on jj',jj,jd,nt
                     go to 140
                  endif
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) then
                     write(iwr,*)'id.eq.ii .and. jd.gt.jj',jj,jd,nt
                     go to 140
                  endif
                  mi(it) = id
                  mj(it) = jd
 70            continue
               endif
c
c     ----- get information about i-shell and j-shell -----
c
               ish = ii
               jsh = jj
               call shells_dft(gout,1,
     &              ish,jsh,ksh,lsh,
     &              basi, basj, bask, basl, ncentr, 
     &              1,idum)

               call ijprim_dft(ncentr)

               if (nij.ne.0) then
                  if(imode.eq.3 .or. imode .eq. 4) then
                     imc=0
                  endif
c
c     ----- kshell -----
c
                  if(ncentr .eq. 4 .or. 
     &              (ncentr .eq. 2) .and. (basi .eq. bask))then
                     k0 = kst0
                     maxkk = ii
                  else
                     k0 = kst0
                     maxkk = nshell(bask)
                  endif

c - used in 3c processing, size of i*j rectangle
                  ibas_num=(maxi-mini)+1
                  jbas_num=(maxj-minj)+1
                  nijx=ibas_num*jbas_num
c
                  if (imode.eq.3) then
c clear scr -- just outside the kk/ll loops: preparation
c for integral contraction
                    call aclear_dp(scr,nijx,0.0d0)
c
                  elseif (imode.eq.4) then
c dvec_fill call moved up to this point from after the kk/ll loops:
c
c generate a rectangular block of density in scr: preparation
c for integral contraction
                     if(rks_sw) then
                        call dvec_fill(locij_dm,
     &                       aaa,ii,jj,mini,maxi,minj,maxj,scr)
                     else
                        call dvec_fill2(locij_dm,aaa,bbb,
     &                       ii,jj,mini,maxi,minj,maxj,scr)
                     endif
c
                     fac = 2.0d0
                     if(ii.eq.jj) fac = 1.0d0
c
                  endif
c
                  do 130 kk = k0 , maxkk
                     kst0 = 1

                     if( ncentr .eq. 4) then
                     do 90 it = 1 , nt
                        kd = iso(kk,it)
                        if (kd.gt.ii) go to 130
                        mk(it) = kd
 90                  continue
                     endif
c
c     ----- lshell ----
c
                     if(basl .lt. 0)then
                        l0 = 1
                        maxll = 1
                     else
                        l0 = lst0
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                     endif

                     do 120 ll = l0 , maxll
                        lst0 = 1
                        if (schwarz_tol.ge.0.and.
     &                      (imode.eq.3.or.imode.eq.4)) then
                           test = schwarz_ao(itrij) * schwarz_cd(kk)
                           oskipp = test.lt.schwarz_lim
                           if(oskipp) then
                              mink = kmin(bask,kk)
                              maxk = kmax(bask,kk)
                              imc=imc+maxk-mink+1
                              nschwz = nschwz + 1
                              go to 120
                           endif
                        endif
                           if( ncentr .eq. 4 ) then
                           n4 = 0
                           do 110 it = 1 , nt
                              ld = iso(ll,it)
                              if (ld.gt.ii) go to 120
                              kd = mk(it)
                              if (kd.lt.ld) then
                                 nd = kd
                                 kd = ld
                                 ld = nd
                              end if
                              id = mi(it)
                              jd = mj(it)
                              if (id.eq.ii .or. kd.eq.ii) then
                                 if (kd.ge.id) then
                                    if (kd.ne.id .or. ld.gt.jd) then
                                       nd = id
                                       id = kd
                                       kd = nd
                                       nd = jd
                                       jd = ld
                                       ld = nd
                                    end if
                                 end if
                                 if (jd.ge.jj) then
                                    if (jd.gt.jj) go to 120
                                    if (kd.ge.kk) then
                                       if (kd.gt.kk) go to 120
                                       if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 120
                                         n4 = n4 + 1
                                       end if
                                    end if
                                 end if
                              end if
 110                       continue
                           q4 = dble(nt)/dble(n4)
                        elseif (ncentr .eq. 3) then

c     rather empirical -it seems triangulation effects
c     are already corrected for - following derivatives
c
                           q4 = 1.0d0

                        elseif (ncentr .eq. 2) then
                           q4 = 1.0d0
                        endif

c
c     ----- (ii,jj//kk,ll) -----
c
                           ish = ii
                           jsh = jj
                           ksh = kk
                           lsh = ll
                           qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- get information about ksh and lsh -----
c
                           call shells_dft(gout,2,
     &                          ish,jsh,ksh,lsh,
     &                          basi, basj, bask, basl,ncentr,
     &                          1,idum)
c
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
                           call genral_dft(gout,ncentr)
c
c ................... use integrals here
c
                           if(imode.eq.1)then


                           else if (imode.eq.3) then
c
c code for imode = 3 / 4 cases rewritten such that the integrals
c are not only used but are also right away contracted as required 
c (the contraction part happened outside the kk/ll loops before)
c
c following code taken from mat_form_3c, but integral contraction inserted
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                scr(ann) = scr(ann) +
     &                                     gout(nn) * aaa(imm)
                              enddo
                            enddo
                          enddo
c update imc:
                              imc=imc+(maxk-mink)+1
c
                        else if (imode.eq.4) then
c
c following code taken from mat_form_3c, but integral contraction inserted
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                ccc(imm) = ccc(imm) +
     &                                     gout(nn)*fac*scr(ann)
                              enddo
                            enddo
                          enddo
c update imc:
                              imc=imc+(maxk-mink)+1

                           else if(imode.eq.6)then
c
c 2 centre cases, just copy integrals into aaa indexed (i,k)
c
                              call mat_form_2c(gout,ii,kk,
     &                             basi, bask,
     &                             imc,aaa,nbasfn(basi))
                           else
                              call caserr('invalid imode')
                           endif

c
c     ----- check cpu time/ maxblock condition -----
c
c-tim                           call chkout(ii,jj,kk,ll,q)
c-tim                           if (omaxb .or. tim.gt.timlim) go to 160

c
 120                 continue
 130              continue


                  if(imode .eq. 3)then 
c 3 centre fock build, cfit in aaa, fock in ccc
c --> contraction moved to before kk/ll loops as well as into them
c
c copy to triangular array (NB, specific to AO basis)
                     locij_t=locij_km
                     call km_fill(locij_t,ii,jj,ccc,scr)
                     locij_km=locij_km+locij_t
                  endif
                  
                  if(imode .eq. 4)then 
c 3 centre density expansion
c --> contraction moved to before kk/ll loops as well as into them
c
                  endif
               end if
               if(imode .eq. 6)then
                  imc=imc+(maxi-mini)+1
               endif
 140        continue
 150     continue

      if (imode .eq. 3 .and. .not.rks_sw) then
         nlen =  ((nbasfn(basi)+1)*nbasfn(basi))/2
         call dcopy(nlen,ccc,1,ddd,1)
      endif
      call final_dft

      if(nprint.ne.-5.and.nschwz.ne.0) write(iwr,6030) nschwz
      end if
c
c     ----- reset core memory
c
c-schw      call gmem_free(ischw)

c ***
c      if (outv) then
c         write (iwr,7010) (intcut(i),i=1,3)
c         if (oimag) then
c            do k = 1 , 11
c               isum = 0
c               do kk = (k-1)*4 + 1 , (k-1)*4 + 4
c                  isum = isum + intmag(kk)
c               enddo
c               intmag(k) = isum
c            enddo
c            write (iwr,7020) (k-34,k=1,42,4)
c            write (iwr,7030) (intmag(k),k=1,11)
c         end if
c      end if
c      if (nprint.ne.-5) write (iwr,7040) cpu
      

c6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
c6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
 6030 format (1x,'schwarz inequality test skipped ',i10,
     +           ' integral blocks')
c7010 format (/,' integral test counts'/1x,30('=')
c    +        /' on ij shell         ',i10/' on ijkl shells      ',
c    +        i10/' on ijkl shells & den',i10/1x,30('='))
c7020 format (' magnitudes of computed integrals '//' 2**',4x,11i8)
c7030 format (' ',8x,11i8)
c7040 format (/' end of 2-electron integral evaluation at ',f8.2)
12345 format (/1x,'Jkint: Computing',i2,'-centre integrals mode ',i2,
     +         ' at ',f8.2,' seconds')
12346 format (1x,'basis sets: ',4i4)
12347 format (1x,'Using cut, tol = ',2i4)

      end

      subroutine jkint_dft_cntgen(iso,gout,nshels,
     &     basi, basj, bask, basl,
     &     imode,aaa,bbb,ccc,ddd, schwarz_ao, schwarz_cd,
     &     te3c_int, nte3c_int, ite3c_stored, nte3c_shl,
     &     te2c_int, nte2c_int, ite2c_stored, nte2c_shl)

      implicit none

      real*8 gout, aaa, bbb, ccc, ddd, schwarz_ao, schwarz_cd
      dimension gout(*), aaa(*), bbb(*), ccc(*), ddd(*)
      dimension schwarz_ao(*), schwarz_cd(*)
      integer basi, basj, bask, basl
      integer imode, nshels
      integer iso(nshels,*)
      integer nte3c_int, nte3c_shl, nte2c_int, nte2c_shl
      real*8 te3c_int, te2c_int
      integer ite3c_stored, ite2c_stored
      dimension te3c_int(nte3c_int), te2c_int(nte2c_int)
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)
c
c...  This subroutine is handed fixed size arrays te3c_store and
c...  te2c_store. Per set of shell labels it determines if it has 
c...  enough space left to store the integrals, and whether the
c...  integrals differ significantly from 0.
c...
c...  If the integrals are 0 they are not stored and the set of shell
c...  labels is flagged to rule out recalculating the integrals.
c...
c...  If the integrals are non-0 and can be stored then they will be
c...  stored on the appropriate array so that they can be loaded from
c...  memory in a later stage.
c...
c...  If the integrals are non-0 and can not be stored the set of shell
c...  labels will be flagged to issue recalculation where needed.
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
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)


cINCLUDE(../m4/common/symtry)
       integer nt
c
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
      integer iwr
      common/dft_iofile/iwr
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
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
cINCLUDE(common/dft_auxvar)
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
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
c inclusion of dft_indez needed for new inlined version of mat_load_3c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      logical odbg
      common/dbgdbg/odbg
c
      real*8 cpulft
      logical opg_root
      integer ipg_nodeid
c
c Local variables
c 
      integer i, j, k
      integer ld, ii, jj, kk, ll, j0,l0, k0
      integer id, jd, kd, nd
      integer ish, jsh, ksh, lsh
      integer ist0, jst0, kst0, lst0
      integer it, itrij, itrkl, ikykk, n4, mij, mjk
      integer ikyii, loop
      integer maxjj, maxkk, maxll
      integer icut0, nschwz, nlen
      integer ann, nn, knn, ltyi, ltyj, ltyk, imm, ijn

      integer ite3c_int, ite2c_int, ntexc_int, ite3c_shl, ite2c_shl

      integer imc, nijx, ibas_num, jbas_num
      integer locij_t, locij_km, locij_dm
      real*8 fac

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      logical oskipp
      real*8 q4, tijl, pi, schwarz_lim, test, schw_max

      real*8 dt0, dt1, tim0, time, tim1

      integer mi, mj, mk, m0
      dimension mi(48),mj(48),mk(48),m0(48)

      integer ib(4,4)

      real*8 scr(225)  

      integer nshl_drct, nshl_incore, nint_nonzero

      integer m25, itmp
      real*8 done, two, twopt5, four, e
      character*1 xn,xt
c
c...  Used in the parallel version
c
      integer liw, ibuff(3), ilen, ifrom
      integer nshl_drct_t, nshl_incore_t, nint_nonzero_t

      data  nt/1/
      data m25/25/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data xn,xt/'n','t'/
      data e /2.30258d0/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /

c
      nshl_incore = 0
      nshl_drct = 0
      nint_nonzero = 0
c
c ***
c *** establish arrays normally set in main line code
c ***

c      dummy iso array
c      note that nt is hardwired to 1
      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
c
c ... /mapperx/
c
      do i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
      enddo
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue

      icut0 = iabs(icut_3c)
      if (schwarz_tol.ge.0) then
         schwarz_lim = done/(10.0d0**schwarz_tol)
      else
         schwarz_lim = -1.0d0
      endif
      cutoff = done/(10.0d0**icut0)
      tol = e*itol_3c
c
c
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_dft
      nopkr = nopk
      iofrst = iofsym
c
      nindmx = 1

      call aux_find(basi)
c
c     ----- two-electron integrals -----
c

c
c Determine which term to compute
c
      ncentr = 0

      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         write(iwr,*)'fitting coefficients',
     &        (ccc(i),i=1,BL_basis_size(bask))
      endif

      if(opg_root().and.(nprint.ne.-5))then
         write(iwr,12345)ncentr,imode,cpulft(1)
         write(iwr,12346)basi,basj,bask,basl
         write(iwr,12347)icut_3c,itol_3c
      endif

      oskipp = .false.

c fail

      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
         igt(loop) = ib(1,ii)
         jgt(loop) = ib(2,ii)
         kgt(loop) = ib(3,ii)
         lgt(loop) = ib(4,ii)
 30   continue


c imode 6
      imc=0
c 3centr
      locij_dm=1
      locij_km=0
      ite3c_int = 0
      ite2c_int = 0
      ntexc_int = 0
      ite3c_shl = 0
      ite2c_shl = 0
c
c     ----- set some parameters -----
c
c     ----- allocate core memory
c
      nschwz = 0
c
      ist0 = 1
      jst0 = 1
      kst0 = 1
      lst0 = 1
c
c     ----- ishell -----
c

c     write(6,*)'Loop i 2'

      if (ist0.le.nshell(basi)) then

c@@ invert for parallel

         do 150 ii = ist0 , nshell(basi)
c
c     ----- eliminate ishell -----
c
            do 50 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) then
                  write(iwr,*)'symm skip on ii',id,ii,nt
                  go to 150
               endif
               m0(it) = id
 50         continue
            ikyii = iky(ii)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst0
               maxjj = ii
            endif

c@@ parallel invert

            do 140 jj = j0 , maxjj

               jst0 = 1
               itrij = ikyii + jj
c
c apply i/j tests only when i and j are AO basis fns
c @@ Note this still dangerous - see comment on triangulation
c in ijprim
c
               if (ncentr .eq. 3  .or. ncentr .eq. 4)then
               do 70 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) then
                     write(iwr,*)'symm skip on jj',jj,jd,nt
                     go to 140
                  endif
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) then
                     write(iwr,*)'id.eq.ii .and. jd.gt.jj',jj,jd,nt
                     go to 140
                  endif
                  mi(it) = id
                  mj(it) = jd
 70            continue
               endif
c
c     ----- get information about i-shell and j-shell -----
c
               ish = ii
               jsh = jj
               call shells_dft(gout,1,
     &              ish,jsh,ksh,lsh,
     &              basi, basj, bask, basl, ncentr, 
     &              1,ntexc_int)

               call ijprim_dft(ncentr)

               if (nij.ne.0) then
                  if(imode.eq.3 .or. imode .eq. 4) then
                     imc=0
                  endif
c
c     ----- kshell -----
c
                  if(ncentr .eq. 4 .or. 
     &              (ncentr .eq. 2) .and. (basi .eq. bask))then
                     k0 = kst0
                     maxkk = ii
                  else
                     k0 = kst0
                     maxkk = nshell(bask)
                  endif

c - used in 3c processing, size of i*j rectangle
                  ibas_num=(maxi-mini)+1
                  jbas_num=(maxj-minj)+1
                  nijx=ibas_num*jbas_num
c
                  if (imode.eq.3) then
c clear scr -- just outside the kk/ll loops: preparation
c for integral contraction
                    call aclear_dp(scr,nijx,0.0d0)
c
                  elseif (imode.eq.4) then
c dvec_fill call moved up to this point from after the kk/ll loops:
c
c generate a rectangular block of density in scr: preparation
c for integral contraction
                     if (rks_sw) then
                        call dvec_fill(locij_dm,
     &                       aaa,ii,jj,mini,maxi,minj,maxj,scr)
                     else
                        call dvec_fill2(locij_dm,aaa,bbb,
     &                       ii,jj,mini,maxi,minj,maxj,scr)
                     endif
c
                     fac = 2.0d0
                     if(ii.eq.jj) fac = 1.0d0
c
                  endif
c
                  do 130 kk = k0 , maxkk
                     kst0 = 1

                     if( ncentr .eq. 4) then
                     do 90 it = 1 , nt
                        kd = iso(kk,it)
                        if (kd.gt.ii) go to 130
                        mk(it) = kd
 90                  continue
                     endif
c
c     ----- lshell ----
c
                     if(basl .lt. 0)then
                        l0 = 1
                        maxll = 1
                     else
                        l0 = lst0
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                     endif

                     do 120 ll = l0 , maxll
                        lst0 = 1
                        ite3c_shl = ite3c_shl + 1
                        if (ite3c_shl.gt.nte3c_shl) then
                           write(*,*)'*** JKINT_DFT_CNTGEN ite3c_shl,',
     &                        'nte3c_shl:',ite3c_shl,nte3c_shl
                           call caserr('ite3c_shl exceeds nte3c_shl')
                        endif
                        if (schwarz_tol.ge.0.and.
     &                      (imode.eq.3.or.imode.eq.4)) then
                           test = schwarz_ao(itrij) * schwarz_cd(kk)
                           oskipp = test.lt.schwarz_lim
                           if(oskipp) then
                              ite3c_stored(ite3c_shl) = 0
                              mink = kmin(bask,kk)
                              maxk = kmax(bask,kk)
                              imc=imc+maxk-mink+1
                              nschwz = nschwz + 1
                              go to 120
                           endif
                        endif
                           if( ncentr .eq. 4 ) then
                           n4 = 0
                           do 110 it = 1 , nt
                              ld = iso(ll,it)
                              if (ld.gt.ii) go to 120
                              kd = mk(it)
                              if (kd.lt.ld) then
                                 nd = kd
                                 kd = ld
                                 ld = nd
                              end if
                              id = mi(it)
                              jd = mj(it)
                              if (id.eq.ii .or. kd.eq.ii) then
                                 if (kd.ge.id) then
                                    if (kd.ne.id .or. ld.gt.jd) then
                                       nd = id
                                       id = kd
                                       kd = nd
                                       nd = jd
                                       jd = ld
                                       ld = nd
                                    end if
                                 end if
                                 if (jd.ge.jj) then
                                    if (jd.gt.jj) go to 120
                                    if (kd.ge.kk) then
                                       if (kd.gt.kk) go to 120
                                       if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 120
                                         n4 = n4 + 1
                                       end if
                                    end if
                                 end if
                              end if
 110                       continue
                           q4 = dble(nt)/dble(n4)
                        elseif (ncentr .eq. 3) then

c     rather empirical -it seems triangulation effects
c     are already corrected for - following derivatives
c
                           q4 = 1.0d0

                        elseif (ncentr .eq. 2) then
                           q4 = 1.0d0
                        endif

c
c     ----- (ii,jj//kk,ll) -----
c
                           ish = ii
                           jsh = jj
                           ksh = kk
                           lsh = ll
                           qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- get information about ksh and lsh -----
c
                           call shells_dft(gout,2,
     &                          ish,jsh,ksh,lsh,
     &                          basi, basj, bask, basl,ncentr,
     &                          1,ntexc_int)
c
                           if(imode.eq.1)then


                           else if (imode.eq.3 .or. imode.eq.4) then
c
                              call genral_dft(gout,ncentr)
                              if (ntexc_int.gt.nte3c_int) then
                                 ite3c_stored(ite3c_shl) = 2
                                 nshl_drct = nshl_drct + 1
                              else
                                 ite3c_stored(ite3c_shl) = 1
                                 call mat_save_3c(gout,mini,maxi,
     &                             minj,maxj,mink,maxk,
     &                             te3c_int,ite3c_int)
                                 nshl_incore = nshl_incore + 1
                              endif
                              nint_nonzero = nint_nonzero +
     &                                       (maxi-mini+1)*
     &                                       (maxj-minj+1)*
     &                                       (maxk-mink+1)
c
c code for imode = 3 / 4 cases rewritten such that the integrals
c are not only used but are also right away contracted as required 
c (the contraction part happened outside the kk/ll loops before)
c
                           if (imode.eq.3) then
c
c following code taken from mat_form_3c, but integral contraction inserted
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                scr(ann) = scr(ann) +
     &                                     gout(nn) * aaa(imm)
                              enddo
                            enddo
                          enddo
c update imc:
                              imc=imc+(maxk-mink)+1
c
                        else if (imode.eq.4) then
c
c following code taken from mat_form_3c, but integral contraction inserted
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                ccc(imm) = ccc(imm) +
     &                                     gout(nn)*fac*scr(ann)
                              enddo
                            enddo
                          enddo
c
                         endif
c update imc:
                              imc=imc+(maxk-mink)+1
c
                           else if(imode.eq.6)then
c
c 2 centre cases, just copy integrals into aaa indexed (i,k)
c
                              call genral_dft(gout,ncentr)
                              ite2c_shl = ite2c_shl + 1
                              if (ite2c_shl.gt.nte2c_shl) then
                                 write(*,*)'*** JKINT_DFT_CNTGEN ',
     &                              'ite2c_shl,nte2c_shl:',
     &                               ite2c_shl,nte2c_shl
                                 call caserr('ite2c_shl exceeds nte2c_sh
     &l')
                              endif
                              if (ntexc_int.gt.nte2c_int) then
                                 ite2c_stored(ite2c_shl) = 2
                              else
                                 ite2c_stored(ite2c_shl) = 1
                                 call mat_save_2c(gout,ii,kk,
     &                                te2c_int,ite2c_int)
                              endif
                              call mat_form_2c(gout,ii,kk,
     &                             basi, bask,
     &                             imc,aaa,nbasfn(basi))
                           else
                              call caserr('invalid imode')
                           endif
c
 120                 continue
 130              continue

                  if(imode .eq. 3)then 
c 3 centre fock build, cfit in aaa, fock in ccc
c --> contraction moved to before kk/ll loops as well as into them
c
c copy to triangular array (NB, specific to AO basis)
                     locij_t=locij_km
                     call km_fill(locij_t,ii,jj,ccc,scr)
                     locij_km=locij_km+locij_t
                  endif
                  
                  if(imode .eq. 4)then 
c 3 centre density expansion
c --> contraction moved to before kk/ll loops as well as into them
c
                  endif
               end if
               if(imode .eq. 6)then
                  imc=imc+(maxi-mini)+1
               endif
 140        continue
 150     continue

      if (imode .eq. 3 .and. .not.rks_sw) then
         nlen =  ((nbasfn(basi)+1)*nbasfn(basi))/2
         call dcopy(nlen,ccc,1,ddd,1)
      endif
      call final_dft

      if(opg_root())then
        if(nprint.ne.-5.and.nschwz.ne.0) write(iwr,6030) nschwz
      endif

      if (imode.eq.3 .or. imode.eq.4) then
       write(iwr,12348) 
     +    nshl_drct, 
     +    100.0d0*dble(nshl_drct)/(nshl_drct+nshl_incore),
     +    nshl_incore,
     +    100.0d0*dble(nshl_incore)/(nshl_drct+nshl_incore),
     +    nint_nonzero
12348  format(1x,'Direct =',i10,' (',f5.1,'%) : InCore = ',i10,
     +        ' (',f5.1,'%) : NonZeroInts =',i15)
      endif

      end if
      

c6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
c6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
 6030 format (1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')

c7010 format (/,' integral test counts'/1x,30('=')
c    +        /' on ij shell         ',i10/' on ijkl shells      ',
c    +        i10/' on ijkl shells & den',i10/1x,30('='))
c7020 format (' magnitudes of computed integrals '//' 2**',4x,11i8)
c7030 format (' ',8x,11i8)
c7040 format (/' end of 2-electron integral evaluation at ',f8.2)
12347 format (1x,'Using cut, tol = ',2i4)
12346 format (1x,'basis sets: ',4i4)
12345 format (/1x,'Jkint: Computing',i2,'-centre integrals mode ',i2,
     +         ' at ',f8.2,' seconds')

      end

      subroutine jkint_dft_genuse(iso,gout,nshels,
     &     basi, basj, bask, basl,
     &     imode,aaa,bbb,ccc,ddd, 
     &     te3c_int, nte3c_int, ite3c_stored, nte3c_shl,
     &     te2c_int, nte2c_int, ite2c_stored, nte2c_shl)

      implicit none

      real*8 gout, aaa, bbb, ccc, ddd, schwarz_ao, schwarz_cd
      dimension gout(*), aaa(*), bbb(*), ccc(*), ddd(*)
      integer basi, basj, bask, basl
      integer imode, nshels
      integer iso(nshels,*)
      integer nte3c_int, nte2c_int, nte3c_shl, nte2c_shl
      integer ite3c_stored, ite2c_stored
      real*8 te3c_int, te2c_int
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)
      dimension te3c_int(nte3c_int), te2c_int(nte2c_int)
c
c...  This subroutine generates/uses the 3 centre and 2 centre
c...  2 electron integrals. 
c
c...  For a given set of shell labels ite3c_stored holds
c...     0  if the 3 centre 2 electron integrals were discarded 
c...        altogether.
c...     1  if the 3 centre 2 electron integrals are needed and stored
c...        on te3c_store.
c...     2  if the 3 centre 2 electron integrals are needed but not 
c...        stored and therefore have to be recomputed.
c...
c...  For a given set of shell labels ite2c_stored holds
c...     0  if the 2 centre 2 electron integrals were discarded 
c...        altogether.
c...     1  if the 2 centre 2 electron integrals are needed and stored
c...        on te3c_store.
c...     2  if the 2 centre 2 electron integrals are needed but not 
c...        stored and therefore have to be recomputed.
c
c...  The tables ite3c_stored, ite2c_stored, ite3c_store and te2c_store
c...  are generated by jkint_dft_cntgen.
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
c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      integer idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif
c
c     ipdiagmode - which diag (IDIAG_PEIGS/_PDSYEV/_PDSYEVD....)
c
c     ipdiagif   - which interface (IDIAG_NO_GAWRAP/IDIAG_GAWRAP)
c
c     idpdiag = dimension for parallel diag
c
c     ipinvmode : what matrix inversion algorithm to use
c                 either the old diag based approach or the more
c                 appropriate ScaLAPACK Cholesky factorisation can 
c                 be used (INV_CHOLESKY or INV_DIAG )
c
      logical odebugp

      common /parcntl/ idpdiag, ipdiagmode, ntchnk, limchnk, iparapr,
     & idpdiis, idpmult2, idporth, ipiomode, iptest, ipinvmode,
     & ipdiagif, odebugp
c
      logical ga_initted
      common/ gainit/ ga_initted

      integer IO_NZ,IO_NZ_S,IO_A
      parameter(IO_NZ=1,IO_NZ_S=2,IO_A=3)

      integer INV_CHOLESKY, INV_DIAG
      parameter(INV_CHOLESKY = 100)
      parameter(INV_DIAG     = INV_CHOLESKY + 1)
c
      integer IDIAG_PEIGS,   IDIAG_PDSYEV, IDIAG_PDSYEVX
      integer IDIAG_PDSYEVD, IDIAG_PDSYEVR
      parameter(IDIAG_PEIGS   = 10)
      parameter(IDIAG_PDSYEV  = IDIAG_PEIGS   + 1)
      parameter(IDIAG_PDSYEVX = IDIAG_PDSYEV  + 1)
      parameter(IDIAG_PDSYEVD = IDIAG_PDSYEVX + 1)
      parameter(IDIAG_PDSYEVR = IDIAG_PDSYEVD + 1)

      integer IDIAG_NO_GAWRAP, IDIAG_GAWRAP
      parameter(IDIAG_NO_GAWRAP=200)
      parameter(IDIAG_GAWRAP=201)


cINCLUDE(../m4/common/symtry)
       integer nt
c
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
      integer iwr
      common/dft_iofile/iwr
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
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
cINCLUDE(common/dft_auxvar)
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
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
c inclusion of dft_indez needed for new inlined version of mat_load_3c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      logical odbg
      common/dbgdbg/odbg
c
      real*8 cpulft
      logical opg_root
c
c Local variables
c 
      integer i, j, k
      integer ld, ii, jj, kk, ll, j0,l0, k0
      integer id, jd, kd, nd
      integer ish, jsh, ksh, lsh
      integer ist0, jst0, kst0, lst0
      integer it, itrij, itrkl, ikykk, n4, mij, mjk
      integer ikyii, loop
      integer maxjj, maxkk, maxll
      integer icut0, nschwz, nlen
      integer ann, nn, knn, ltyi, ltyj, ltyk, imm, ijn

      integer ndum
      integer ite3c_int, ite2c_int, ite3c_shl, ite2c_shl

      integer imc, nijx, ibas_num, jbas_num
      integer locij_t, locij_km, locij_dm
      real*8 fac

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      logical oskipp
      real*8 q4, tijl, pi, schwarz_lim, test, schw_max

      real*8 dt0, dt1, tim0, time, tim1

      integer mi, mj, mk, m0
      dimension mi(48),mj(48),mk(48),m0(48)

      integer ib(4,4)

      real*8 scr(225)  

      integer m25, itmp
      real*8 done, two, twopt5, four, e
      character*1 xn,xt

      data  nt/1/
      data ib/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data m25/25/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data xn,xt/'n','t'/
      data e /2.30258d0/
c
c
c ***
c *** establish arrays normally set in main line code
c ***

c      dummy iso array
c      note that nt is hardwired to 1
      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
c
c ... /mapperx/
c
      do i = 1 , maxorb
         k = i*(i-1)/2
         iky(i) = k
      enddo
c
c ... /piconx/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
c
c ... /rootx/
c
      do 20 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 20   continue

      icut0 = iabs(icut_3c)
      if (schwarz_tol.ge.0) then
         schwarz_lim = done/(10.0d0**schwarz_tol)
      else
         schwarz_lim = -1.0d0
      endif
      cutoff = done/(10.0d0**icut0)
      tol = e*itol_3c
c
c
c     check for pure sp basis. if so, do gaussian integrals.
      call spchck_dft
      nopkr = nopk
      iofrst = iofsym
c
      nindmx = 1

      call aux_find(basi)
c
c     ----- two-electron integrals -----
c

c
c Determine which term to compute
c
      ncentr = 0

      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         write(iwr,*)'fitting coefficients',
     &        (ccc(i),i=1,BL_basis_size(bask))
      endif

      if(opg_root().and.(nprint.ne.-5))then
         write(iwr,12345)ncentr,imode,cpulft(1)
         write(iwr,12346)basi,basj,bask,basl
         write(iwr,12347)icut_3c,itol_3c
      endif
      oskipp = .false.

c fail

      ii = 1
      if (opdbas) then
         ii = 2
      end if
      if (opfbas) then
         ii = 3
      end if
      if (opgbas) then
         ii = 4
      end if
      do 30 loop = 1 , 3
         igt(loop) = ib(1,ii)
         jgt(loop) = ib(2,ii)
         kgt(loop) = ib(3,ii)
         lgt(loop) = ib(4,ii)
 30   continue


c imode 6
      imc=0
c 3centr
      locij_dm=1
      locij_km=0

c
c     ----- set some parameters -----
c
c     ----- allocate core memory
c
      nschwz = 0
      ite3c_int = 0
      ite2c_int = 0
      ite3c_shl = 0
      ite2c_shl = 0
c
      ist0 = 1
      jst0 = 1
      kst0 = 1
      lst0 = 1
c
c     ----- ishell -----
c
c     write(6,*)'Loop i 3'


      if (ist0.le.nshell(basi)) then

c@@ invert for parallel

         do 150 ii = ist0 , nshell(basi)
c
c     ----- eliminate ishell -----
c
            do 50 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) then
                  write(iwr,*)'symm skip on ii',id,ii,nt
                  go to 150
               endif
               m0(it) = id
 50         continue
            ikyii = iky(ii)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst0
               maxjj = ii
            endif

c@@ parallel invert

            do 140 jj = j0 , maxjj

               jst0 = 1
               itrij = ikyii + jj
c
c apply i/j tests only when i and j are AO basis fns
c @@ Note this still dangerous - see comment on triangulation
c in ijprim
c
               if (ncentr .eq. 3  .or. ncentr .eq. 4)then
               do 70 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) then
                     write(iwr,*)'symm skip on jj',jj,jd,nt
                     go to 140
                  endif
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) then
                     write(iwr,*)'id.eq.ii .and. jd.gt.jj',jj,jd,nt
                     go to 140
                  endif
                  mi(it) = id
                  mj(it) = jd
 70            continue
               endif
c
c     ----- get information about i-shell and j-shell -----
c
               ish = ii
               jsh = jj
               call shells_dft(gout,1,
     &              ish,jsh,ksh,lsh,
     &              basi, basj, bask, basl, ncentr, 
     &              1,ndum)

               call ijprim_dft(ncentr)

               if (nij.ne.0) then
                  if(imode.eq.3 .or. imode .eq. 4) then
                     imc=0
                  endif
c
c     ----- kshell -----
c
                  if(ncentr .eq. 4 .or. 
     &              (ncentr .eq. 2) .and. (basi .eq. bask))then
                     k0 = kst0
                     maxkk = ii
                  else
                     k0 = kst0
                     maxkk = nshell(bask)
                  endif

c - used in 3c processing, size of i*j rectangle
                  ibas_num=(maxi-mini)+1
                  jbas_num=(maxj-minj)+1
                  nijx=ibas_num*jbas_num
c
                  if (imode.eq.3) then
c clear scr -- just outside the kk/ll loops: preparation
c for integral contraction
                    call aclear_dp(scr,nijx,0.0d0)
c
                  elseif (imode.eq.4) then
c dvec_fill call moved up to this point from after the kk/ll loops:
c
c generate a rectangular block of density in scr: preparation
c for integral contraction
                     if (rks_sw) then
                        call dvec_fill(locij_dm,
     &                       aaa,ii,jj,mini,maxi,minj,maxj,scr)
                     else
                        call dvec_fill2(locij_dm,aaa,bbb,
     &                       ii,jj,mini,maxi,minj,maxj,scr)
                     endif
c
                  endif
c
                  do 130 kk = k0 , maxkk
                     kst0 = 1

                     if( ncentr .eq. 4) then
                     do 90 it = 1 , nt
                        kd = iso(kk,it)
                        if (kd.gt.ii) go to 130
                        mk(it) = kd
 90                  continue
                     endif
c
c     ----- lshell ----
c
                     if(basl .lt. 0)then
                        l0 = 1
                        maxll = 1
                     else
                        l0 = lst0
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                     endif

                     do 120 ll = l0 , maxll
                        ite3c_shl = ite3c_shl + 1
                        lst0 = 1
                        if (imode.eq.3.or.imode.eq.4) then
                           if(ite3c_stored(ite3c_shl).eq.0) then
                              mink = kmin(bask,kk)
                              maxk = kmax(bask,kk)
                              imc=imc+maxk-mink+1
                              nschwz = nschwz + 1
                              go to 120
                           endif
                        endif
                        if( ncentr .eq. 4 ) then
                           n4 = 0
                           do 110 it = 1 , nt
                              ld = iso(ll,it)
                              if (ld.gt.ii) go to 120
                              kd = mk(it)
                              if (kd.lt.ld) then
                                 nd = kd
                                 kd = ld
                                 ld = nd
                              end if
                              id = mi(it)
                              jd = mj(it)
                              if (id.eq.ii .or. kd.eq.ii) then
                                 if (kd.ge.id) then
                                    if (kd.ne.id .or. ld.gt.jd) then
                                       nd = id
                                       id = kd
                                       kd = nd
                                       nd = jd
                                       jd = ld
                                       ld = nd
                                    end if
                                 end if
                                 if (jd.ge.jj) then
                                    if (jd.gt.jj) go to 120
                                    if (kd.ge.kk) then
                                       if (kd.gt.kk) go to 120
                                       if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 120
                                         n4 = n4 + 1
                                       end if
                                    end if
                                 end if
                              end if
 110                       continue
                           q4 = dble(nt)/dble(n4)
                        elseif (ncentr .eq. 3) then

c     rather empirical -it seems triangulation effects
c     are already corrected for - following derivatives
c
                           q4 = 1.0d0

                        elseif (ncentr .eq. 2) then
                           q4 = 1.0d0
                        endif

c
c     ----- (ii,jj//kk,ll) -----
c
                        ksh = kk
c
c     ----- initialize gout to zero -----
c     ----- get information about ksh and lsh -----
c
ccc  --> shells_dft call moved from here into the following alternative
c
c ................... generate/load and use integrals here
c
c code for imode = 3 / 4 cases rewritten such that the integrals
c are not only used but are also right away contracted as required 
c (the contraction part happened outside the kk/ll loops before)
c
                        if (imode.eq.3 .or. imode.eq.4) then
c
                          if (ite3c_stored(ite3c_shl).eq.1) then
c
c generate pointers that came from shells_dft before (only maxk, mink)
                            mink = kmin(bask,ksh)
                            maxk = kmax(bask,ksh)
c
                          else
c
c shells_dft call moved into this case from outside this alternative
                           ish = ii
                           jsh = jj
                           lsh = ll
                           qq4 = q4
                            call shells_dft(gout,2,
     &                           ish,jsh,ksh,lsh,
     &                           basi, basj, bask, basl,ncentr,
     &                           1,ndum)
c
                            call genral_dft(gout,ncentr)
                          endif
c
                        endif
c
                        if (imode.eq.3) then
c
c following code taken from mat_load_3c / mat_form_3c, respectively
c but integral contraction inserted
c
                         if (ite3c_stored(ite3c_shl).eq.1) then
c
c in this case, code from mat_load_3c
c
                          ann=0
                          nn=ite3c_int
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                imm=imm+1
                                nn=nn+1
c here the contraction
                                scr(ann) = scr(ann) +
     &                                     te3c_int(nn) * aaa(imm)
                              enddo
                            enddo
                          enddo
                          ite3c_int=nn
c
                         else
c
c in this case, code from mat_form_3c
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                scr(ann) = scr(ann) +
     &                                     gout(nn) * aaa(imm)
                              enddo
                            enddo
                          enddo
c
                         endif
c update imc:
                         imc = imc+(maxk-mink)+1
c
                        else if (imode.eq.4) then
c
                         fac = 2.0d0
                         if(ii.eq.jj) fac=1.0d0
c
                         if (ite3c_stored(ite3c_shl).eq.1) then
c
c in this case, code from mat_load_3c
c
                          ann=0
                          nn=ite3c_int
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                imm=imm+1
                                nn=nn+1
c here the contraction
                                ccc(imm) = ccc(imm) +
     &                                     te3c_int(nn)*fac*scr(ann)
                              enddo
                            enddo
                          enddo
                          ite3c_int=nn
c
                         else
c
c in this case, code from mat_form_3c
c
                          ijn=0
                          ann=0
                          do ltyi=mini,maxi
                            do ltyj=minj,maxj
                              ijn=ijn+1
                              knn=0
                              ann=ann+1
                              imm=imc
                              do ltyk=mink,maxk
                                knn=knn+1
                                imm=imm+1
                                nn=ijgt(ijn)+klgt(knn)
c here the contraction
                                ccc(imm) = ccc(imm) +
     &                                     gout(nn)*fac*scr(ann)
                              enddo
                            enddo
                          enddo
c
                         endif
c update imc:
                         imc = imc+(maxk-mink)+1
c
c
                        elseif(imode.eq.1)then
c
c
                        else if(imode.eq.6)then

c
c shells_dft call moved into this case from outside the 'imode' 
c alternative (I am not sure whether the whole routine is needed here)
c
                           ish = ii
                           jsh = jj
                           lsh = ll
                           qq4 = q4
                            call shells_dft(gout,2,
     &                           ish,jsh,ksh,lsh,
     &                           basi, basj, bask, basl,ncentr,
     &                           1,ndum)
c
c 2 centre cases, just copy integrals into aaa indexed (i,k)
c
                           ite2c_shl = ite2c_shl+1
                           if (ite2c_stored(ite2c_shl).eq.0) then
                           elseif (ite2c_stored(ite2c_shl).eq.1) then
                              call mat_load_2c(ii,kk,basi,bask,
     &                             imc,aaa,nbasfn(basi),
     &                             te2c_int,ite2c_int)
                           else
                              call genral_dft(gout,ncentr)
                              call mat_form_2c(gout,ii,kk,
     &                             basi, bask,
     &                             imc,aaa,nbasfn(basi))
                           endif
                        else
                           call caserr('invalid imode')
                        endif
c
c     ----- end of kk/ll loops -----
 120                 continue
 130              continue

                  if(imode .eq. 3)then 
c 3 centre fock build, cfit in aaa, fock in ccc
c --> contraction moved to before kk/ll loops as well as into them
c
c copy to triangular array (NB, specific to AO basis)
                     locij_t=locij_km
                     call km_fill(locij_t,ii,jj,ccc,scr)
                     locij_km=locij_km+locij_t
                  endif
                  
                  if(imode .eq. 4)then 
c 3 centre density expansion
c --> contraction moved to before kk/ll loops as well as into them
c
                  endif
               end if
               if(imode .eq. 6)then
                  imc=imc+(maxi-mini)+1
               endif
c
c     ----- end of ii/jj loops -----
 140        continue
 150     continue

      if (imode .eq. 3 .and. .not.rks_sw) then
         nlen =  ((nbasfn(basi)+1)*nbasfn(basi))/2
         call dcopy(nlen,ccc,1,ddd,1)
      endif
      call final_dft

      if(opg_root())then
         if(nprint.ne.-5.and.nschwz.ne.0) write(iwr,6030) nschwz
      endif
      end if
      

c6010 format (i4,3i5,1x,i10,i9,f11.2,f9.2)
c6020 format (i4,3i5,1x,i10,9x,f11.2,f9.2)
 6030 format (1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')

c7010 format (/,' integral test counts'/1x,30('=')
c    +        /' on ij shell         ',i10/' on ijkl shells      ',
c    +        i10/' on ijkl shells & den',i10/1x,30('='))
c7020 format (' magnitudes of computed integrals '//' 2**',4x,11i8)
c7030 format (' ',8x,11i8)
c7040 format (/' end of 2-electron integral evaluation at ',f8.2)
12347 format(1x,'Using cut, tol = ',2i4)
12346 format(1x,'basis sets: ',4i4)
12345 format(/1x,'Jkint: Computing',i2,'-centre integrals mode ',i2,
     +         ' at ',f8.2,' seconds')

      end
