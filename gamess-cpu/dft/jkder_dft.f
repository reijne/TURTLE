      subroutine jkder_dft(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,schwarz_ao,schwarz_cd)

      implicit none
c
c arguments
c
      character*8 zscftp
      integer iso
      real*8 q
      integer nshels
      dimension iso(nshels,*),q(*)
      integer basi, basj, bask, basl
      real*8 adens(*), bdens(*)
      real*8 cfit(*), cfit2(*)
      real*8 grad(3,*)
      real*8 schwarz_ao, schwarz_cd
      dimension schwarz_ao(*), schwarz_cd(*)

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

cc      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w

c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
cINCLUDE(../m4/common/restar)
cINCLUDE(../m4/common/restri)
c
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb

        real*8 dummy

c
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c

cINCLUDE(../m4/common/symtry)
      integer nt
      data nt/1/

c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c

      real*8 dgout
      common/tgrad/dgout(9)

      integer mxp2
      parameter (mxp2 = mxprms*mxprms)

      integer ncmax
      parameter(ncmax=65)

      real*8 ddij,ddkl,aei,aej,aek,ael,aa,r,x1,y1,z1,dd
      integer ijden,ik,ijx,ijy,ijz,klx,kly,klz
      real*8 dij,dkl
      integer ijgt,klgt

      common/dft_d2escr/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      integer iwr
      common/dft_iofile/iwr
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
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c


c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c

c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb

      logical odbg
      common/dbgdbg/odbg
c
c local variables
c      
      integer m1, m2, m3, m0
      dimension m0(48),m1(48),m2(48),m3(48)

      integer najkl, nabcl, nd2str, nactp, nshdm
      integer numcmo, n4
      integer id3, id4, id5, id6, id7, id8, id9
      integer i10, i20, i30, i40, i00
      integer iabd1, iabd2, iabd3
      integer ijshel
      integer i, j, issi, lnddm, lensec, m
      integer kt_max, inc1_max, ncmmm_max
      integer maxll, maxjj, maxkk, it
      logical ofpres, odpres, ogpres

      logical olab, olabc, olabcd
      integer ii, jj, kk, ll
      integer i0, j0, k0, l0
      integer l2,  m00
      integer klshel, kadi, id, jd, kd, ld, nd
      integer iblok
      integer nav, ifmp1, ifmp2, ifmp3
      integer iwor1, iwor2, iwor3
      integer ndgout, next
      integer kadij, kadijk, iceni
      logical omp2, ocifor, omp2w, ohf, ociopt, oskipp
      real*8 tolij, tolijk, abmax
      real*8 dtim, tim0
      real*8 schwarz_lim, test
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf

      character*9 fnm,snm

      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      real*8 gmax

      integer itmp1, itmp2
      integer itmp(4), bas
c     
c functions/stmnt fns
c
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, allocate_memory2
      logical opg_root
c-      integer indq

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft'/
c
c-      indq(m,n) = (m-1)*n
c
c-      ompir = runtyp.eq.zdipd .or. runtyp.eq.zinfra
c-      ompir = ompir .and. omp2
c
c     ----- check for grhf or gvb cases -----
c
c-      ouhf = zscftp.eq.zuhf
c-      orgvb = zscftp.eq.zgvb
c-      ogrhf = zscftp.eq.zgrhf
c-      oclos = zscftp.eq.zrhf
c-      omcscf = zscftp.eq.zmcscf
c
c     ----- check for casscf
c
c-      ocas = zscftp.eq.zcas
c-      if (ocas .or. mp3 .or. (omp2 .and. .not.ompir)) then
c-         call setsto(1360,0,ijkl)
c-      end if

c
c disable prefactor testing
c
c      if (nprint.ne.-5 .and. oprint(57)) then
c         write (iwr,6010)
c         call writel(q(iprefa),nshels)
c      end if

c      dummy iso
c@@ need to decide how big this array should
c   be, and only use it when centres are interchangeable
c

      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
      if (schwarz_tol.ge.0) then
         schwarz_lim = 0.1d0**schwarz_tol
      else
         schwarz_lim = -1.0d0
      endif
      nschwz = 0
c
c this routine forces the load of the block data
c
       call dummy_intgx
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
         if(basi .eq. bask)then
         write(6,*)'fitting coefficients',(cfit(i),i=1,nbasfn(basi))
         else
         write(6,*)'fitting coefficients 1',(cfit(i),i=1,nbasfn(basi))
         write(6,*)'fitting coefficients 2',(cfit2(i),i=1,nbasfn(basi))
         endif
      endif

      if(print_sw(DEBUG_FORCES) .and. opg_root())then
        write(6,*)'Computing',ncentr,' centre deriv integrals'
      endif
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

c-      oeof = .false.

cc      if( icontij .eq. DENSITY ) then

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

cc      endif


      nat3 = natoms*3
c      nbsq = num*num
c      lenb = lensec(nx)

      odpres = .false.
      ofpres = .false.
      ogpres = .false.

      itmp(1) = basi
      itmp(2) = basj
      itmp(3) = bask
      itmp(4) = basl
      do j=1,4
      if(itmp(j) .lt. 0) itmp(j) = itmp(1)
      do 20 i = 1 , nshell(itmp(j))
         bas = itmp(j)
         if (ktype(bas,i).eq.3) odpres = .true.
         if (ktype(bas,i).eq.4) ofpres = .true.
         if (ktype(bas,i).eq.5) ogpres = .true.
 20   continue
      enddo
c

      m = ntpdm + 9
c     if (omp2w .or. ompir) m = m + 3
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c     code for dgenrl preallocation 
c     buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do j=1,4
      bas = itmp(j)
      do 21 i = 1 , nshell(bas)
         kt_max=max(kt_max,ktype(bas,i))
 21   continue
      enddo
      inc1_max=(kt_max+1)**4
c
c     vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6

c      write(6,*)'dfg',odpres, ofpres,ogpres
c
c     ----- set pointers for partitioning of core -----
c

c i10 used for nconf

      i10 = null_memory()
cc      i20 = igmem_alloc(l2)
cc      i30 = igmem_alloc(l2)

      ndens = 1
c     if (ouhf) then
c        ndens = 2
c     end if
c     if (orgvb) then
c        ndens = 3
c     end if
c     if (ogrhf) then
c        ndens = njk
c     end if
c     if (ocas) then
c        ndens = 1
c     end if
c     if (ofock .or. ompir) then
c        ndens = ndens + ndenin
c     end if

cc      if( icontij .eq. DENSITY ) then

         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

cc      else
cc         ida = 0
cc         i20 = 0
cc         i30 = 0
cc      endif

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
c-         if (ogrhf) nfok = njk*natoms*6
      end if
      if (ofokab) nfok = nat3*ndenin
c

c-      ifok = igmem_alloc(nx*nfok)
c-      ifok = ifok - 1

c  core at i00 is for indexing
      i00 = allocate_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = allocate_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

c      write(6,*)'core i20',i20
c      call chkadr(q(i20))

      ic7 = allocate_memory2(lnddm*9+inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

      call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
     +     ist,jst,kst,lst,ncentr,q)


      next = ipg_dlbtask()

      if(print_sw(DEBUG_PARALLEL))then
         write(6,*)'Task info on entry: ',ipg_nodeid(),next,icount_dlb
      endif

      if (ist.le.nshels) then
c
c     ----- ishell -----
c
         do 140 ii = ist , nshell(basi)

c-            kadi = kad(ii)
            ijshel = ii*(ii-1)/2

c
            if(ncentr .eq. 3 .or. ncentr .eq. 4)then
            do 40 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 140
               m0(it) = id
 40         continue
            endif

            iceni = katom(basi,ii)

c-            if (omcscf) call mcajkl(q(id3),q(id5),q(id7),q(id8),q(id9),
c-     +                              ii,nactp)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst
               maxjj = ii
            endif

            do 130 jj = j0 , maxjj
c-               kadij = kadi + kad(jj)
               jst = 1
               itrij = ijshel+jj

c-               tolij = dlntol + q(iprefa+ijshel+jj-1)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
c                  - DISABLED -
c-               if (tolij.gt.-3.401d0) then
                 if(.true.)then
c
c apply i/j tests only when i and j are AO basis fns
c
                  if (ncentr .eq. 3  .or. ncentr .eq. 4)then
                  do 60 it = 1 , nt
                     id = m0(it)
                     jd = iso(jj,it)
                     if (jd.gt.ii) go to 130
                     if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (id.eq.ii .and. jd.gt.jj) go to 130
                     m1(it) = id
                     m2(it) = jd
 60               continue
                  endif

                  if(basj.lt.0)then
                     olab = .true.
                  else
                     olab = katom(basj,jj).eq.iceni
                  endif

c
c     store information about the pair (ij)
c
                  call dshell_dft(1,ii,jj,kk,ll,
     +                 basi, basj, bask, basl,ncentr)

                  call dprim_dft

                  if (nij.ne.0) then

c-                     if (omcscf) call mcabkl(q(id3),q(id7),q(id4),ii,jj,
c-     +                   nactp)

                     icount_dlb = icount_dlb + 1
                     if(icount_dlb . eq. next) then

                        if(print_sw(DEBUG_PARALLEL))then
                         write(6,*)'Task ',next,' on node',ipg_nodeid()
                        endif

c
c     ----- kshell -----
c
                     if(ncentr .eq. 4 .or. 
     &                 (ncentr .eq. 2 .and. basi .eq. bask))then
                        k0 = kst
                        maxkk = ii
                     else
                        k0 = kst
                        maxkk = nshell(bask)
                     endif

                     if(odbg)write(6,*)'kk loop',k0, maxkk

                     do 120 kk = k0 , maxkk

c-                        kadijk = kadij + kad(kk)
                        kst = 1
                        klshel = kk*(kk-1)/2

                        if( ncentr .eq. 4) then
                        do 80 it = 1 , nt
                           kd = iso(kk,it)
                           if (kd.gt.ii) go to 120
                           m3(it) = kd
 80                     continue
                        endif

                        olabc = olab .and. katom(bask,kk).eq.iceni

c-                        if (omcscf)
c-     +                      call mcabcl(q(id3),q(id4),q(id6),q(id8),
c-     +                      q(id9),ii,jj,kk,nactp)

                        if(basl .lt. 0)then
                           l0 = 1
                           maxll = 1
                        else
                           l0 = lst
                           maxll = kk
                           if (kk.eq.ii) maxll = jj
                        endif

                        do 110 ll = l0 , maxll

                           lst = 1
                           if (ncentr.eq.3.and.schwarz_tol.ge.0) then
                              test = schwarz_ao(itrij) * schwarz_cd(kk)
                              oskipp = test.lt.schwarz_lim
                              if(oskipp) then
c                                mink = kmin(bask,kk)
c                                maxk = kmax(bask,kk)
c                                imc=imc+maxk-mink+1
                                 nschwz = nschwz + 1
                                 go to 110
                              endif
                              
                           endif

c-                           if (kadijk+kad(ll).lt.0) then
c-                              tolijk = tolij + q(iprefa+klshel+ll-1)
c-                              if (tolijk.gt.0.0d0) then

                           if(basl.lt.0)then
                              olabcd = olabc 
                           else
                              olabcd = olabc .and. 
     &                             katom(basl,ll).eq.iceni
                           endif

                           if(odbg)write(6,*)'olab',ii,jj,kk,ll,olabcd

                                 if (.not.(olabcd)) then

                                 if( ncentr .eq. 4 ) then

                                    n4 = 0
                                    do 100 it = 1 , nt
                                       ld = iso(ll,it)
                                       if (ld.gt.ii) go to 110
                                       kd = m3(it)
                                       if (kd.lt.ld) then
                                         nd = kd ! swap k,l
                                         kd = ld
                                         ld = nd
                                       end if
                                       id = m1(it)
                                       jd = m2(it)
                                       if (id.eq.ii .or. kd.eq.ii) then
                                         if (kd.ge.id) then
                                         if (kd.ne.id .or. ld.gt.jd)
     +                                      then
                                         nd = id  ! swap i,k
                                         id = kd
                                         kd = nd
                                         nd = jd  ! swap j,l
                                         jd = ld
                                         ld = nd
                                         end if
                                         end if
                                         if (jd.ge.jj) then
                                         if (jd.gt.jj) go to 110
                                         if (kd.ge.kk) then
                                         if (kd.gt.kk) go to 110
                                         if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 110
                                         n4 = n4 + 1
                                         end if
                                         end if
                                         end if
                                       end if
 100                                continue
c
c     ----- calculate q4 factor for this group of shells -----
c
                                    q4 = dble(nt)/dble(n4)

                                 elseif (ncentr .eq. 3) then

c rather empirical -it seems triangulation effects
c are already corrected for
c a factor of 2 is applied in dabab
                                    q4 = 1.0d0

                                 elseif (ncentr .eq. 2) then
                                    q4 = 1.0d0
                                 endif
c
c     ----- check for redundant combinations -----
c
                              call redund_dft(ii,jj,kk,ll,
     +                             basi, basj, bask, basl, 
     +                             iwr)


                              if (npass.eq.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk
                                 else if(ncentr .eq. 2)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,kk
                                 endif

                              endif
                              
                              if (npass.ne.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,q4
                                 else if(ncentr .eq. 2)then
                                    if(ii.eq.14.and.kk.eq.8)then
c                                       odbg=.true.
                                    else
c                                       odbg = .false.
                                    endif
                                    if(odbg)
     &               write(6,*)'Shell block',ii,kk, q4
                                 endif
c
c     ----- initialize dgout to zero -----
c
                          call vclr(dgout,1,ndgout)

                          call dshell_dft(2,ii,jj,kk,ll,
     +                         basi, basj, bask, basl,ncentr)
c
c     ----- form products of density matrix elements -----
c
                          call vclr(q(iabd+1),1,lendd*ntpdm)


                          call dabab_dft(ii,jj,kk,ll,
     &                         basi, basj, bask, basl,
     +                         ncentr,q4,
     +                         zscftp,adens,bdens,cfit,cfit2,
     +                         q(iabd+1))

c                          call chkadr2(q(i20),itmp1)
c                          call chkadr2(q(iabd+1),itmp2)
c                          write(6,*)'after dabab',itmp1,itmp2,
c     &                         q(i20),q(iabd+1)



c
c     ----mess about with the density matix by eliminating
c     zero elements
c

                                 call delim_dft(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)

                                 
                                 if (ijkld.ne.0) then
                                    call dgenrl_dft(q(1),q(i00),q(i00),
     +                                   abmax,gmax)

c-                                   if (ofock)
c-     +                                call fockd2(q,q(i00))
c-                                   if (ofokab) then
c-                                     call fokabd(q,q(i00))
c-                                   end if
c
c     ----- generate all 4 partial contributions to the gradient ----
c
                                   call formeg_dft
                                 end if              ! ijkld.ne.o
                              end if                 ! olabcd
                                 end if              ! npass .ne. 0

c-                              end if                 ! (tolijk.gt.0.0d0)
c-                           end if                    ! (kadijk+kad(ll).lt.0)

 110                    continue                     ! ll   loop
 120                 continue
                     next = ipg_dlbtask()
                     endif
                  end if
               end if
 130        continue

c
c     ----- save gradient and restart data -----
c                ==== disabled =====
c            call dfinal_dft(q,0,ii,basi,grad,nat)
c            if (tim.ge.timlim) go to 150

 140     continue
         call pg_dlbpush

      end if
c
c     ----- end of *shell* loops -----
c

      if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
       call caserr('no symmetry in jkder_dft')
ccc       isymd = igmem_alloc(nw196(6))
c
ccc       call symde(q(isymd),natoms)

c
c     ----- reset core memory from symde
c
ccc      call gmem_free(isymd)
c
      endif
      call dfinal_dft(q,1,ii,basi,grad, natoms)

c150  continue
      call timit(0)
c     dtim = tim - tim0
c
c     ----- reset core memory from jkder -----
c

c-      if(omcscf)then
c-         call gmem_free(id9)
c-         call gmem_free(id8)
c-         call gmem_free(id7)
c-         call gmem_free(id6)
c-         call gmem_free(id5)
c-         call gmem_free(id4)
c-         call gmem_free(id3)
c-         call gmem_free(i40)
c-      endif

cc
cc revert
cc      call gmem_free(ic1)

      ic7 = ic7 + 1
      call free_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call free_memory2(iabd,'d',fnm,snm,'iabd')
      call free_memory2(i00,'d',fnm,snm,'i00')
c-      ifok = ifok + 1
c-      call gmem_free(ifok)

cc      if( icontij .eq. DENSITY ) then
        ida = ida + 1
cc      endif

cc      call gmem_free(i40)
cc      call gmem_free(i30)
cc      call gmem_free(i20)

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c
