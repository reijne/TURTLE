c---- memory counting routines -----------------------------------------
      subroutine memreq_jkder_dft(zscftp,q,iso,nshels,
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
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
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

      character*16 fnm,snm

      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      real*8 gmax

      integer itmp(4), itmp1, itmp2, bas, idum
c     
c functions/stmnt fns
c
      integer memreq_pg_dgop
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, incr_memory2
      logical opg_root

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e_memory.m'/
      data snm/'memreq_jkder_dft'/
c
c     do i = 1,nshels
c        iso(i,1) = 1
c     enddo

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
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

      nat3 = natoms*3

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
      i10 = null_memory()

      ndens = 1
         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
      end if
      if (ofokab) nfok = nat3*ndenin
c
c  core at i00 is for indexing
      i00 = incr_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = incr_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

c      write(6,*)'core i20',i20
c      call chkadr(q(i20))

      ic7 = incr_memory2(lnddm*9+inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

c     call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
c    +     ist,jst,kst,lst,ncentr,q)

c     next = ipg_dlbtask()
c
c     ----- ishell -----
c
c     ----- jshell -----
c
c     ----- kshell -----
c
c     ----- calculate q4 factor for this group of shells -----
c
c     ----- check for redundant combinations -----
c
c     ----- initialize dgout to zero -----
c
c     ----- form products of density matrix elements -----
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
c     ----- generate all 4 partial contributions to the gradient ----
c
c     ----- save gradient and restart data -----
c
c     call dfinal_dft(q,1,ii,basi,grad, natoms)
         idum = memreq_pg_dgop(3*natoms,'+')
c
      ic7 = ic7 + 1
      call decr_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call decr_memory2(iabd,'d',fnm,snm,'iabd')
      call decr_memory2(i00,'d',fnm,snm,'i00')

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c
      subroutine memreq_jkder_dft_genuse(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,
     &  ite3c_stored, nte3c_shl, ite2c_stored, nte2c_shl)

      implicit none
c
c     This routine is the derivative equivalent of jkint_dft_genuse.
c     The results of the Schwarz inequality are stored in tables
c     so that we can lookup whether a block of integrals was used in
c     the energy evaluation. 
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
      integer nte3c_shl, nte2c_shl, ite3c_stored, ite2c_stored
      integer ite3c_shl, ite2c_shl
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)

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
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
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
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf
      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      real*8 gmax

      integer itmp(4), itmp1, itmp2, bas, idum

      character *9  fnm
      character *16 snm
c     
c functions/stmnt fns
c
      integer memreq_pg_dgop
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, incr_memory2
      integer iipsci
      logical opg_root, oipsci

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft_genuse'/
c
c     do i = 1,nshels
c        iso(i,1) = 1
c     enddo

      odbg = .false.
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
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

      nat3 = natoms*3

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
      i10 = null_memory()

      ndens = 1

         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
      end if
      if (ofokab) nfok = nat3*ndenin
c
c  core at i00 is for indexing
      i00 = incr_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = incr_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

      ic7 = incr_memory2(lnddm*9 + inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                   'ic7')
      ic7 = ic7 - 1

c     call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
c    +     ist,jst,kst,lst,ncentr,q)
c
c     ----- ishell -----
c
c     ----- jshell -----
c
c     ----- kshell -----
c
c     ----- calculate q4 factor for this group of shells -----
c
c     ----- check for redundant combinations -----
c
c     ----- initialize dgout to zero -----
c
c     ----- form products of density matrix elements -----
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
c     ----- generate all 4 partial contributions to the gradient ----
c
c     ----- save gradient and restart data -----
c
c     ----- end of *shell* loops -----
c
c     ----- allocate core memory for symde
c
c     ----- reset core memory from symde
c
c     call dfinal_dft(q,1,ii,basi,grad, natoms)
         idum = memreq_pg_dgop(3*natoms,'+')
c
c     ----- reset core memory from jkder -----
c
      ic7 = ic7 + 1
      call decr_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call decr_memory2(iabd,'d',fnm,snm,'iabd')
      call decr_memory2(i00,'d',fnm,snm,'i00')

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c---- routines that do the real work -----------------------------------
      subroutine jkder_dft_genuse(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,
     &  ite3c_stored, nte3c_shl, ite2c_stored, nte2c_shl)

      implicit none
c
c     This routine is the derivative equivalent of jkint_dft_genuse.
c     The results of the Schwarz inequality are stored in tables
c     so that we can lookup whether a block of integrals was used in
c     the energy evaluation. 
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
      integer nte3c_shl, nte2c_shl, ite3c_stored, ite2c_stored
      integer ite3c_shl, ite2c_shl
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)

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
c-      real*8 gijkl
c-      integer mword,nlenx,kworx,kworxx
c-      common/blkin/gijkl(510),mword,nlenx,kworx,kworxx
c-      integer ijkl
c-      common/craypk/ijkl(1360)

c-      real*8 gijkl1
c-      integer mwor1,nlen1,kwor1,kwor11
c-      common/blk1/gijkl1(510),mwor1,nlen1,kwor1,kwor11

c-      integer ijkl1
c-      common/sortpk/ijkl1(1360)

c-      real*8 gijkl2
c-      integer mwor2,nlen2,kwor2,kwor22
c-      common/bufc/gijkl2(510),mwor2,nlen2,kwor2,kwor22

c-      integer ijkl2
c-      common/three/ijkl2(1360)

c-      real*8 gijkl3
c-      integer mwor3,nlen3,kwor3,kwor33
c-      common/bufd/gijkl3(510),mwor3,nlen3,kwor3,kwor33

c-      integer ijkl3
c-      common/lsort/ijkl3(1360)

c-      real*8 dipd, dipn, dipi
c-      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)

c-      real*8 val, vall
c-      integer icnt, mxtr
c-      common/dbuf/val(312),vall(195),icnt,mxtr(4)

c-      integer iao, jao, kao, lao, iprt
c-      common/dlabs/iao(312),jao(312),kao(312),lao(312),iprt(312)

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
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf
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

      character *9  fnm
      character *16 snm
c     
c functions/stmnt fns
c
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, allocate_memory2
      integer iipsci
      logical opg_root, oipsci
c-      integer indq

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft_genuse'/
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

c-      if (ompir) then
c-         ndenin = 3
c-         ntpdm = 4
c-         iflden = idaf
c-         call secget(isect(31),31,ibden)
c-         iwor1 = 0
c-         iwor2 = 0
c-         iwor3 = 0
c-         mwor1 = 0
c-         mwor2 = 0
c-         mwor3 = 0
c-         ifmp1 = 20
c-         ifmp2 = 21
c-         ifmp3 = 22
c-         ib1 = 1
c-         call search(ib1,ifmp1)
c-         call search(ib1,ifmp2)
c-         call search(ib1,ifmp3)
c-         call setsto(1360,0,ijkl1)
c-         call setsto(1360,0,ijkl2)
c-         call setsto(1360,0,ijkl3)
c-         call secget(isect(57),57,iblok)
c-         call rdedx(dipd,lds(isect(57)),iblok,ifild)
c-      end if
c-      if (ogrhf) then
c-         m = 0
c-         call secget(isect(53),m,iblok)
c-         call rdedx(nact,lds(isect(53)),iblok,idaf)
c-      end if

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

      ic7 = allocate_memory2(lnddm*9 + inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

c-      if(omcscf)then
c-
c-         i40 = igmem_alloc(l2)
c-
c-ccc         id3 = i40 + nx
c-         numcmo = num*ncoorb
c-         nactp = ncact*(ncact+1)/2
c-         nd2str = indq(nactp+1,16)
c-         if (odpres) nd2str = indq(nactp+1,36)
c-         if (ofpres) nd2str = indq(nactp+1,100)
c-         if (ogpres) nd2str = indq(nactp+1,225)
c-         id3 = igmem_alloc(numcmo)
c-         i30 = i20
c-ccc         id4 = id3 + numcmo
c-         id4 = id3 + numcmo
c-c     nd2mo=ind(nactp,nactp)
c-         nd2mo = nactp*(nactp+1)/2
c-c     write(6,*)' length of tpdm',nd2mo
c-
c-         id4 = igmem_alloc(nd2str)
c-ccc         id5 = id4 + nd2str
c-         id5 = igmem_alloc(nd2mo)
c-ccc         id6 = id5 + nd2mo
c-c
c-         najkl = indq(ncact*4+1,nactp)
c-         if (odpres) najkl = indq(ncact*6+1,nactp)
c-         if (ofpres) najkl = indq(ncact*10+1,nactp)
c-         if (ogpres) najkl = indq(ncact*15+1,nactp)
c-c
c-         nabcl = indq(16+1,4*ncact)
c-         if (odpres) nabcl = indq(36+1,6*ncact)
c-         if (ofpres) nabcl = indq(100+1,10*ncact)
c-         if (ogpres) nabcl = indq(225+1,10*ncact)
c-c
c-         id6 = igmem_alloc(nabcl)
c-cc         id7 = id6 + nabcl
c-         id7 = igmem_alloc(najkl)
c-cc         id8 = id7 + najkl
c-         id8 = igmem_alloc(nactp)
c-cc         id9 = id8 + nactp
c-c
c-         nshdm = max(ncact,4)
c-         if (odpres) nshdm = max(ncact,6)
c-         if (ofpres) nshdm = max(ncact,10)
c-         if (ogpres) nshdm = max(ncact,15)
c-c
c-         id9 = igmem_alloc(nshdm*nshdm)
c-cc
c-cc         iabd = id9 + nshdm*nshdm
c-c
c-c   iabd should not be used for mcscf calculations - tpdm in id5
c-c
c-cc         i00 = iabd + lnddm*ntpdm
c-cc         ic7 = i00 + lnddm*(3/nav+1) - 1
c-      endif

c     tim0 = cpulft(1)
c
c-      if (ocas) then
c-         call dbutci(ist,jst,kst,lst)
c-      else if (omcscf) then
c-         call ddebut(zscftp,q(id3),q(i30),q(i10),q(i40),q(id5),
c-     +               ist,jst,kst,lst,q)
c-      else

      call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
     +     ist,jst,kst,lst,ncentr,q)

c-      end if
c-      if (mp3 .or. (omp2 .and. .not.ompir)) call search(iblk2d,ifil2d)

c-      mword = 0
c-      iword = 0

c-      kworx = 999
c-      kwor1 = 999
c-      kwor2 = 999
c-      kwor3 = 999
c
c check 
c??      kloc(nshels+1) = num + 1

c
c-      icnt = 0
c-      ib1 = 1
c-      if (omp2w) then
c-         call search(ib1,mpstrm(1))
c-         if (odebug(30)) write (iwr,6020)
c-      end if

c
c     ----- initialise static load balancing counter ------
c
      icount_dlb = iipsci()

      if(print_sw(DEBUG_PARALLEL))then
         write(6,*)'Task info on entry: ',ipg_nodeid(),icount_dlb
      endif

      ite3c_shl = 0
      ite2c_shl = 0

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
                     if(.not.oipsci()) then

                        if(print_sw(DEBUG_PARALLEL))then
                         write(6,*)'Task ',icount_dlb,' on node',
     &                             ipg_nodeid()
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
                           ite3c_shl = ite3c_shl + 1
                           lst = 1
                           if (ncentr.eq.3) then
                              if(ite3c_stored(ite3c_shl).eq.0) then
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
c                                   ite2c_shl = ite2c_shl+1
c                                   if (ite2c_stored(ite2c_shl).eq.0)
c    +                                 goto 110
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

c-                                 if (omcscf)
c-     +                              call mcabcd(q(id3),q(id6),
c-     +                              q(iabd+1),q(id9),ii,jj,kk,ll,
c-     +                              q4)
c-                                 if (mp3 .or.
c-     +                              (omp2 .and. .not.ompir))
c-     +                              call mcdab(q(iabd+1),ii,jj,kk,
c-     +                              ll,q4)
c-                                 if (ompir) then
c-                                   iabd1 = iabd + lendd + 1
c-                                   iabd2 = iabd1 + lendd
c-                                   iabd3 = iabd2 + lendd
c-                                   call dpdab1(q(iabd1),ii,jj,kk,
c-     +                                ll,q4,ifmp1,iwor1)
c-                                   call dpdab2(q(iabd2),ii,jj,kk,
c-     +                                ll,q4,ifmp2,iwor2)
c-                                   call dpdab3(q(iabd3),ii,jj,kk,
c-     +                                ll,q4,ifmp3,iwor3)
c-                                   call tpdder(q(iabd1),lendd,
c-     +                                q(i20),q(i30),nx,ndenin,ii,
c-     +                                jj,kk,ll,q4)
c-                                 end if
c-                                 if (.not.orgvb .and.
c-     +                              .not.ocas .and. .not.omcscf)

                          call dabab_dft(ii,jj,kk,ll,
     &                         basi, basj, bask, basl,
     +                         ncentr,q4,
     +                         zscftp,adens,bdens,cfit,cfit2,
     +                         q(iabd+1))

c                          call chkadr2(q(i20),itmp1)
c                          call chkadr2(q(iabd+1),itmp2)
c                          write(6,*)'after dabab',itmp1,itmp2,
c     &                         q(i20),q(iabd+1)


c-                                 if (orgvb)
c-     +                              call dabg(ii,jj,kk,ll,l1,norb,
c-     +                              q4,q(i20),q(i30),q(i10)
c-     +                              ,onocor,onopen,q(iabd+1))
c-                                 if (ocas)
c-     +                              call dabci(ii,jj,kk,ll,q4,
c-     +                              oeof,q(iabd+1))
c-                                 if (omcscf)
c-     +                              call dabmc(ii,jj,kk,ll,q4,
c-     +                              q(i30),q(i40),q(iabd+1))

c
c     ----mess about with the density matix by eliminating
c     zero elements
c

                                 call delim_dft(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)

c-                                 if (ompir) then
c-                                   call delim2(q(iabd1),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd2),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd3),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                 end if
                                 
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
c-      if (omp2w) then
c-         if (icnt.ne.0) then
c-            call pack(vall,8,iao,1560)
c-            call wrt3s(val,511,mpstrm(1))
c-         end if
c-         icnt = 0
c-         m00 = 0
c-         call put(val,m00,mpstrm(1))
c-         call shut1(mpstrm(1))
c-         if (odebug(30)) write (iwr,6030)
c-      end if

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
cc      endif

cc      call gmem_free(i40)
cc      call gmem_free(i30)
cc      call gmem_free(i20)

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
      subroutine ddebut_dft(zscftp,da,db,nconf,dd,fock,
     + ista,jsta,ksta,lsta,ncentr,q)
c
c      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
c      implicit character *8 (z),character *1 (x)
c      implicit character *4 (y)

      implicit none
c
c arguments
c
      character*8 zscftp
      real*8 da, db, fock, dd, q
      integer nconf
      dimension da(*),db(*),nconf(*),fock(*),dd(*), q(*)
      integer  ista,jsta,ksta,lsta,ncentr

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

c-      logical ociopt, omp2, ohf
c-      common/restrl/ociopt(2),omp2

c retained for normf normp itol nprint icut 
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

      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      integer iwr
      common/dft_iofile/iwr
c
      real*8 de
      common /dft_grad2/ de(3,maxat)
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
c functions
c
      integer ipg_nnodes
c
c local variables
c
      integer ioffv, ioffd, iofft
      integer itold, loop, ioff
      integer nbase, ncol, nop
      integer i, j, ij, k, kk, is, nsi, ni, ic, nc, np2, m, ils, icutd
      integer iblok
      
      real*8 factor, dum, duma, dumb

      character*8 zrhf,zuhf,zgrhf,zgvb,zmcscf
      real*8 done, ten, e
c
c  Functions
c
      logical opg_root

c
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
      data zmcscf/'mcscf'/
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
c
c      if (onocnt) write (iwr,6010)
c
c-      if (zscftp.eq.zgvb) then
c-c
c-c     ----- set up core density matrix, read in  eigenvectors for
c-c            gvb -----
c-c
c-         call rdedx(db,l3,ibl3qa,idaf)
c-         call tdown(db,ilifq,db,ilifq,l1)
c-         call dencor(da,db,l1)
c-c
c-c     ----- set up mo to fock operator pointers in nconf -----
c-c
c-         if (nco.ne.0) call setsto(nco,1,nconf)
c-         ic = 0
c-         if (nseto.gt.0) then
c-            nbase = ncores
c-            do 30 i = 1 , nseto
c-               nop = no(i)
c-               do 20 j = 1 , nop
c-                  nconf(ic+nco+j) = nbase + 1
c- 20            continue
c-               ic = ic + nop
c-               nbase = nbase + 1
c- 30         continue
c-         end if
c-         if (npair.gt.0) then
c-            np2 = npair + npair
c-            do 40 i = 1 , np2
c-               nconf(i+nco+ic) = ncores + nseto + i
c- 40         continue
c-         end if
c-         norb = nco + npair + npair + ic
c-         onocor = .false.
c-         onopen = .false.
c-         onocor = nco.eq.0
c-         onopen = nseto.eq.0 .and. npair.eq.0
c-c
c-      else if (zscftp.eq.zuhf) then
c-c
c-c     ----- read in density matrices (alpha+beta) in uhf -----
c-c
c-         call rdedx(da,l2,ibl3pa,idaf)
c-         call rdedx(db,l2,ibl3pb,idaf)
c-         do 50 i = 1 , l2
c-            duma = da(i)
c-            dumb = db(i)
c-            da(i) = duma + dumb
c-            db(i) = duma - dumb
c- 50      continue
c-c
c-      else if (zscftp.eq.zrhf) then
c-c
c-c     ----- read in density matrix  in rhf  -----
c-c
c-c-         call rdedx(da,l2,ibl3pa,idaf)
c-c      write(6,*)'input density',(da(i),i=1,l2)
c-c-c
c-      else if (zscftp.eq.zgrhf) then
c-c
c-c     general scf
c-c
c-         call rdedx(fock,l3,ibl3qa,idaf)
c-         call tdown(fock,ilifq,fock,ilifq,l1)
c-         m = 0
c-         do 90 is = 1 , njk
c-            call vclr(da(m+1),1,l2)
c-            nsi = nbshel(is)
c-            ils = ilfshl(is)
c-            do 80 ni = 1 , nsi
c-               nc = iactiv(ils+ni)
c-               ncol = (nc-1)*l1
c-               ij = 0
c-               do 70 i = 1 , l1
c-                  dum = fock(i+ncol)
c-                  do 60 j = 1 , i
c-                     ij = ij + 1
c-                     da(ij+m) = da(ij+m) + dum*fock(j+ncol)
c- 60               continue
c- 70            continue
c- 80         continue
c-            m = m + l2
c- 90      continue
c-c
c-      else if (zscftp.eq.zmcscf) then
c-c
c-c     mcscf/multi
c-c
c-         m = 0
c-         call secget(isecmo,m,iblok)
c-         call rdedx(da,l1*ncoorb,iblok+mvadd,idaf)
c-         m = 0
c-         call secget(isecdd,m,iblok)
c-         if(odebug(30)) write (iwr,6040) iblok
c-         call rdedx(dd,l2,iblok,idaf)
c-         if(odebug(30)) write (iwr,6050) ifil2d,iblk2d
c-         call rdedx(fock,nd2mo,iblk2d,ifil2d)
c-         if (ncore.gt.0) then
c-            call vclr(db,1,l2)
c-            ij = 0
c-            do 120 i = 1 , l1
c-               do 110 j = 1 , i
c-                  ij = ij + 1
c-                  do 100 k = 1 , ncore
c-                     kk = (k-1)*l1
c-                     db(ij) = db(ij) + da(i+kk)*da(j+kk)
c- 100              continue
c- 110           continue
c- 120        continue
c-            write (iwr,6020) ncore
c-         end if
c-      else
c-         call caserr('invalid scftype detected in gradient code')
c-      end if


c
c     ----- set starting parameters -----
c
      outd= nprint.eq. - 4

c-      if (ofokab .or. ompir) then
c-c
c-c      read in ndenin density matrices from some external file
c-c      matrices are in mo basis
c-c
c-         ioff = 1
c-         call search(ibden,iflden)
c-         do 130 loop = 1 , ndenin
c-            call reads(db(ioff),l2,iflden)
c-            if (odebug(21)) call prtris(db(ioff),l1,iwr)
c-            ioff = ioff + l2
c- 130     continue
c-
c-         ioffd = igmem_alloc(l2)
c-         iofft = igmem_alloc(l1)
c-         ioffv = igmem_alloc(l3)
c-         m = 0
c-         call rdedx(q(ioffv),l3,ibl3qa,idaf)
c-         call tdown(q(ioffv),ilifq,q(ioffv),ilifq,l1)
c-         ioff = 1
c-         do 140 loop = 1 , ndenin
c-            call dcopy(l2,db(ioff),1,q(ioffd),1)
c-            call demoao(q(ioffd),db(ioff),q(ioffv),q(iofft),l1,
c-     +  ncoorb,l1)
c-            ioff = ioff + l2
c- 140     continue
c-
c-         call gmem_free(ioffv)
c-         call gmem_free(iofft)
c-         call gmem_free(ioffd)
c-
c-      end if

         icutd = icut
         icutd = max(icutd,8)
         itold = itol + 1
c
c        if(opg_root())then
c           write(6,*)'old toler parameters',itold, icutd
c        endif

      if(ncentr .eq. 2) then
         icutd = icutd_2c
         itold = itold_2c
      else if(ncentr .eq. 3) then
         icutd = icutd_3c
         itold = itold_3c
      else if(ncentr .eq. 4) then
         icutd = icutd_4c
         itold = itold_4c
      endif

c     if(opg_root())then
c       write(6,*)'new toler parameters',itold, icutd
c     endif

      itold = max(itold,16)
      dcut = done/(ten**icutd)
      tol1 = e*(itold+2)
      tol2 = e*itold
      tol3 = done/ten**(itold+2)
      tol4 = done/ten**itold
      onorm = normf.ne.1 .or. normp.ne.1

c
c =====  Starting points for shell loops  =====
c      ( replace if implementing restarts)

      ista = 1
      jsta = 1
      ksta = 1
      lsta = 1
c
c ====  zero 2e dft gradient accumulator =====
c
      call dscal(natoms*3,0.0d0,de,1)

c-      if (ofock .or. ofokab) then
c-         if(odebug(30)) write (iwr,6030)
c-         call vclr(fock,1,nfok*l2)
c-      end if

      return
c6010 format (/1x,22('=')/1x,'gradient of the energy'/1x,22('='))
c6020 format (//1x,' number of frozen and core orbitals',i5//)
c6030 format (/1x,'zeroing storage for fock matrices')
c6040 format (/1x,'reading density matrix from block',i5)
c6050 format (/1x,'tpdm from file, block', 2i5)
      end
      subroutine delim_dft(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd,
     * abmax)

      implicit none

      logical oform
      integer ijgt, klgt
      real*8 qa, qb
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      integer ij, kl, ijkl
      integer lendd

      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen

c-      logical ociopt,omptwo,ohf,omp2w,ojunk
c-      common/restrl/ociopt(2),omptwo,ohf(7),omp2w

ccINCLUDE(../m4/common/specal)
c
c local variables
c
      integer i, k, n, nn
      real*8 ab, abmax

      integer itmp1

      do 20 i = 1 , lendd
         oform(i) = .false.
 20   continue
      abmax = 1.0d0

c-      if (.not.(ofock .or. ofokab .or. ompir .or. omp2w)) then
         abmax = 0.0d0
         nn = 0
         do 40 i = 1 , ij
            do 30 k = 1 , kl
               nn = nn + 1
               n = ijgt(i) + klgt(k)
               ab = dabs(qa(n))
ccccc               if (ab.lt.dcut) oform(nn) = .true.
               if (ab.gt.abmax) abmax = ab
 30         continue
 40      continue
c-      end if

      call dcopy(lendd,qa,1,qb,1)
      nn = 0
      ijkl = 0
      do 60 i = 1 , ij
         do 50 k = 1 , kl
            nn = nn + 1
            if (.not.(oform(nn))) then
               n = ijgt(i) + klgt(k)
               ijkl = ijkl + 1
               qa(ijkl) = qb(n)
            end if
 50      continue
 60   continue

      return
      end
      subroutine denfac_dft(dkl,csk,cpk,cdk,cfk,cgk,
     +                      csl,cpl,cdl,cfl,cgl,
     +                     mink,maxk,minl,maxl,okandl,double)
      implicit real*8  (a-h,o-z)
       logical okandl,double
       dimension dkl(*)

      logical odbg
      common/dbgdbg/odbg

      if (.not.double) then
         n = 0
         max = maxl
         do 130 k = mink , maxk
            if (okandl) max = k
            go to (20,30,60,60,
     +             40,60,60,60,60,60,
     +             50,60,60,60,60,60,60,60,60,60,
     +             55,60,60,60,60,60,60,60,60,60,
     +             60,60,60,60,60) , k
 20         dum1 = csk
            go to 60
 30         dum1 = cpk
            go to 60
 40         dum1 = cdk
            go to 60
 50         dum1 = cfk
            go to 60
 55         dum1 = cgk
 60         do 120 l = minl , max
               go to (70,80,110,110,
     +                90,110,110,110,110,110,
     +               100,110,110,110,110,110,110,110,110,110,
     +               105,110,110,110,110,110,110,110,110,110,
     +               110,110,110,110,110) , l
 70            dum2 = dum1*csl
               go to 110
 80            dum2 = dum1*cpl
               go to 110
 90            dum2 = dum1*cdl
               go to 110
 100           dum2 = dum1*cfl
               go to 110
 105           dum2 = dum1*cgl
 110           n = n + 1
               dkl(n) = dum2

               if(odbg)write(6,*)'dkl',n,dum2

 120        continue
 130     continue
      else
         n = 0
         max = maxl
         do 250 k = mink , maxk
            if (okandl) max = k
            go to (140,150,180,180,
     +             160,180,180,180,180,180,
     +             170,180,180,180,180,180,180,180,180,180,
     +             175,180,180,180,180,180,180,180,180,180,
     +             180,180,180,180,180) , k
 140        dum1 = csk
            go to 180
 150        dum1 = cpk
            go to 180
 160        dum1 = cdk
            go to 180
 170        dum1 = cfk
            go to 180
 175        dum1 = cgk
 180        do 240 l = minl , max
               go to (190,200,230,230,
     +                210,230,230,230,230,230,
     +                220,230,230,230,230,230,230,230,230,230,
     +                225,230,230,230,230,230,230,230,230,230,
     +                230,230,230,230,230) , l
 190           dum2 = dum1*csl
               if (k.gt.1) then
                  dum2 = dum2 + csk*cpl
               else
                  dum2 = dum2 + dum2
               end if
               go to 230
 200           dum2 = dum1*(cpl+cpl)
               go to 230
 210           dum2 = dum1*(cdl+cdl)
               go to 230
 220           dum2 = dum1*(cfl+cfl)
               go to 230
 225           dum2 = dum1*(cgl+cgl)
 230           n = n + 1
               dkl(n) = dum2
               if(odbg)write(6,*)'dkl',n,dum2
 240        continue
 250     continue
      end if
      return
      end
      subroutine dfinal_dft(q,index,ii,basi,grad, nat)

      implicit none

      real*8 q(*)
      integer index, ii, basi
      real*8 grad(3,*)
      integer nat

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


c***   ***node-MPP***
c***    writing (and gopping) for index ne 0 is disabled
c***    (or if other integrals still follow)
c***    so no restarts (nb. and no d-gradients!!)
c***    this saves quite a bit of bother
c***    restarts may be reeabled by dividing the partial grads by
c***    the number of nodes after gopping before writing
c***    (the restarting job must have same number then (tricky)
c***   ***node-MPP***


      integer iwr
      common/dft_iofile/iwr
c
      real*8 de
      common /dft_grad2/ de(3,maxat)
c
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
c functions
c
      logical opg_root
c
c local variables
c
      integer n, i, j
      integer min, max, ncoord

      character*4 ydnam
      dimension ydnam(3)

c-      real*8 cpu, cpulft
c-      real*8 f
c-      integer lensec
c-      integer istnu, jstnu, kstnu, lstnu
c-      integer iblok
c-      integer k
c 
      data ydnam /'e/x','e/y','e/z'/
c
c
c Restart block - commented out in jkder_dft
c
      if (index.ne.1) then
c
c     ----- get restart data -----
c
c         irest = 6
c         istnu = 1 + ii
c         jstnu = 1
c         kstnu = 1
c         lstnu = 1
c
c     ----- save gradient + restart data -----
c
c       call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
c
c     ----- check cpu time -----
c
c        if (istnu.gt.nshell(basi)) return
c          call texit(0,irest)
c          if (tim.lt.timlim) return
c         write (iwr,6010) tim , istnu , jstnu , kstnu , lstnu
c
c-         if (ofokab) call dabout(q,odebug(21),iwr)
c-         if (ofock) call hsymd(q,iwr)
c         write (iwr,6030)
c         return

c-      else if (onocnt) then
      else
c
c     ----- transfer gradient into -grad- -----
c
         ncoord = 3*nat
c
c sum contributions
c
         call pg_dgop(902,de,ncoord,'+')

         if(opg_root() .and. print_sw(DEBUG_FORCES) )then

            max = 0
 20         min = max + 1
            max = max + 8
            if (max.gt.nat) max = nat
            write(iwr,*)'Coulomb gradient from DFT'
            write (iwr,6020)
            write (iwr,6050) (i,i=min,max)
            write (iwr,6020)
            do 30 n = 1 , 3
               write (iwr,6060) ydnam(n) , (de(n,i),i=min,max)
 30         continue

            if (max.lt.nat) go to 20
         endif
c
c Add into the total gradient
c
         do i = 1, nat
            do j = 1,3
               grad(j,i) = grad(j,i) + de(j,i)
            enddo
         enddo
c
c-         if (ofokab) then
c-            call dabout(q,odebug(21),iwr)
c-            iochf(1) = iochf(1) + nat*ndenin*3*lensec(nx)
c-         end if
c-         if (ofock) then
c-            call hsymd(q,iwr)
c-            call clredx
c-         end if
c-         if (ompir) then
c-            call secget(isect(57),57,iblok)
c-            do 60 i = 1 , 3
c-               do 50 j = 1 , 3
c-                  do 40 k = 1 , nat
c-                     dipd(i,j,k) = dipd(i,j,k) + dipi(i,j,k)
c- 40               continue
c- 50            continue
c- 60         continue
c-            call wrt3(dipd,lds(isect(57)),iblok,idaf)
c-         end if
c-         call revind
c-         call clredx
c-         cpu = cpulft(1)
c-         if (nprint.ne.-5) write (iwr,6040) cpu
         return
      end if

c6010 format (/' insufficient time to complete evaluation of 2-electron'
c    +        ,' contribution to gradient'//' job dumped at ',f10.2,
c    +        ' seconds'//' you can forget next batch ',4i5)
 6020 format (/)
c6030 format (//10x,27('*')/10x,'*** warning ***'/10x,
c    +        'this job must be restarted'/10x,27('*')/)
c6040 format (/' end of calculation of the energy gradient at ',f8.2,
c    +        ' seconds'/)
 6050 format (1x,'atom',8(6x,i3,6x))
 6060 format (3x,a3,8f15.7)
      end
      subroutine dform_dft(x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
c      implicit real*8  (a-h,o-z)
      implicit none

c
c arguments
c
      integer ix(*),iy(*),iz(*)
      integer ncdim
      real*8 x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     1 yd(ncdim,*),zd(ncdim,*),g(*)

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
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
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

      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c

      real*8 t1, t2, t3
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
c local variables
c
      logical unroll
      integer mx, my, mz
      integer n1, n2, n3, nr, n
      real*8 dsum
c
cvd$r assoc
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5 .and. ijkld.ge.16
      if (.not.unroll .or. ncontr.gt.5) then
         do 30 n = 1 , ijkld
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            do 20 nr = 1 , ncontr
               t1(nr) = xd(nr,mx)*y(nr,my)*z(nr,mz)
               t2(nr) = x(nr,mx)*yd(nr,my)*z(nr,mz)
               t3(nr) = x(nr,mx)*y(nr,my)*zd(nr,mz)
 20         continue
            g(n1+n) = g(n1+n) + dsum(ncontr,t1,1)
            g(n2+n) = g(n2+n) + dsum(ncontr,t2,1)
            g(n3+n) = g(n3+n) + dsum(ncontr,t3,1)
 30      continue
         return
      else
         go to (40,60,80,100,120) , ncontr
      end if
 40   do 50 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz)
 50   continue
      return
 60   do 70 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz)
 70   continue
      return
 80   do 90 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz)
 90   continue
      return
 100  do 110 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +            *y(4,my)*z(4,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +            *yd(4,my)*z(4,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +            *y(4,my)*zd(4,mz)
 110  continue
      return
 120  do 130 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +            *y(4,my)*z(4,mz) + xd(5,mx)*y(5,my)*z(5,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +            *yd(4,my)*z(4,mz) + x(5,mx)*yd(5,my)*z(5,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +            *y(4,my)*zd(4,mz) + x(5,mx)*y(5,my)*zd(5,mz)
 130  continue
      return
      end
      subroutine dforma_dft(spij,spkl,noform,
     *                  x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
c      implicit real*8  (a-h,o-z)
      implicit none
c
c arguments
c
      real*8 x,y,z,xd,yd,zd,g
      integer ix,iy,iz
      dimension ix(*),iy(*),iz(*)
      logical noform,spij,spkl
      dimension noform(*)
      integer ncdim
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     & yd(ncdim,*),zd(ncdim,*),g(*)

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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
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

      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c

      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
c local variables
c
      logical unroll
      integer n, n1, n2, n3, nr, nn, i, k
      integer mx, my, mz
      real*8 t1, t2, t3, s1, s2, s3
      real*8 dsum
c
cvd$r assoc
c
      n = 0
      nn = 0
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5
      if (.not.unroll) then
c
         if (.not.spij) then
            do 40 i = 1 , ijd
               do 30 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     do 20 nr = 1 , ncontr
                        t1(nr) = ddkl(nr,k)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        t2(nr) = ddkl(nr,k)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        t3(nr) = ddkl(nr,k)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 20                  continue
                     g(n+n1) = dsum(ncontr,t1,1) + g(n+n1)
                     g(n+n2) = dsum(ncontr,t2,1) + g(n+n2)
                     g(n+n3) = dsum(ncontr,t3,1) + g(n+n3)
                  end if
 30            continue
 40         continue
         else if (.not.spkl) then
            do 70 i = 1 , ijd
               do 60 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 50 nr = 1 , ncontr
                        s1 = s1 + ddij(nr,i)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 50                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 60            continue
 70         continue
         else
            do 100 i = 1 , ijd
               do 90 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 80 nr = 1 , ncontr
                        s1 = s1 + (ddij(nr,i)*ddkl(nr,k))*xd(nr,mx)
     +                       *y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *y(nr,my)*zd(nr,mz)
 80                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 90            continue
 100        continue
         end if
         return
      else
         go to (110,140,170,200,230) , ncontr
      end if
 110  do 130 i = 1 , ijd
         do 120 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz)
            end if
 120     continue
 130  continue
      return
 140  do 160 i = 1 , ijd
         do 150 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz)
            end if
 150     continue
 160  continue
      return
 170  do 190 i = 1 , ijd
         do 180 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz)
            end if
 180     continue
 190  continue
      return
 200  do 220 i = 1 , ijd
         do 210 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz)
            end if
 210     continue
 220  continue
      return
 230  do 250 i = 1 , ijd
         do 240 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz) + (ddij(5,i)*ddkl(5,k))
     +                   *xd(5,mx)*y(5,my)*z(5,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*yd(5,my)
     +                   *z(5,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*y(5,my)
     +                   *zd(5,mz)
            end if
 240     continue
 250  continue
      return
      end
      subroutine dgenrl_dft(qq,iqq,noform,abmax,gmax)

      implicit none
c
c arguments
c
      integer iqq(*)
      real*8 qq(*), abmax
      logical noform(*)

      real*8 gmax

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
      real*8 ag,csa,cpa,cda,cfa,cga,bg,csb,cpb,cdb,cfb,cgb
      real*8 cgg,csc,cpc,cdc,cfc,cgc,dg,csd,cpd,cdd,cfd,cgd
      real*8 xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk
      real*8 exij,rsmall
      integer nga,ngb,ngc,ngd
      common/dshlnf/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +              cfa(mxprms),cga(mxprms),
     +              bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +              cfb(mxprms),cgb(mxprms),
     +             cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +              cfc(mxprms),cgc(mxprms),
     +              dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +              cfd(mxprms),cgd(mxprms),
     +              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     +              nga,ngb,ngc,ngd,exij(mxprms*mxprms),rsmall
c
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c

      real*8 bp01,b00,b10,xcp00,xc00,ycp00,yc00,zcp00
      real*8 zc00,f00,dxij,dyij,dzij,dxkl,dykl,dzkl
      integer in,kn,ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     +   xc00(ncmax),ycp00(ncmax),yc00(ncmax),zcp00(ncmax),
     +   zc00(ncmax),f00(ncmax),
     +   dxij,dyij,dzij,dxkl,dykl,dzkl,
     +   in(12),kn(12),ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2

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

      real*8 axak,ayak,azak,axai,
     1  ayai,azai,abv,aandbv,rhov,
     2  xxv,c1xv,c2xv,c3xv,c4xv,
     3  c1yv,c2yv,c3yv,c4yv,
     4  c1zv,c2zv,c3zv,c4zv,expev
      common/bufb/axak(mxp2),ayak(mxp2),azak(mxp2),axai(mxp2),
     1  ayai(mxp2),azai(mxp2),abv(mxp2),aandbv(mxp2),rhov(mxp2),
     2  xxv(mxp2),c1xv(mxp2),c2xv(mxp2),c3xv(mxp2),c4xv(mxp2),
     3  c1yv(mxp2),c2yv(mxp2),c3yv(mxp2),c4yv(mxp2),
     4  c1zv(mxp2),c2zv(mxp2),c3zv(mxp2),c4zv(mxp2),expev(mxp2)

ccccINCLUDE(common/segm)

c
c local variables
c
      real*8 xb, yb, zb
      real*8 bxbi, bxbk
      real*8 bybi, bybk
      real*8 bzbi, bzbk
      real*8 brrk, bbrrk
      real*8 csl, cpl, cdl, cgl, cfl
      real*8 csk, cpk, cdk, cgk, cfk
      real*8 akxk, akyk, akzk
      real*8 xc, yc, zc
      real*8 xd, yd, zd
      real*8 exkl, dddd, dijd, dkld
      real*8 ai, aj, ak, al, u2, dum, dum2, expe
      real*8 b, binv
      integer ig, jg, lg, i, ii, iii, is, inc
      integer k, n, ncmmm, nn, nnn0, m
      integer max, maxlg, kg, jgmax

      integer itmp1, itmp2

      real*8 one
      real*8 pi252

      logical trduij,spij,spkl,trdukl
      logical double
      logical sptru

      logical odbg
      common/dbgdbg/odbg


      data one/1.0d0/
c     data zero,pt5/0.0d0,0.5d0/
      data pi252/34.986836655250d0/
c
      if (ijkld.eq.0) return
      if (ijd.eq.1 .and. kld.eq.1) then
         call ssdss_dft(qq)
      else
         ni = lit
         if (oskip(1)) ni = lit - 1
         nj = ljt
         if (oskip(2)) nj = ljt - 1
         nk = lkt
         if (oskip(3)) nk = lkt - 1
         nl = llt
         if (oskip(4)) nl = llt - 1
         kln2 = 1
         kln1 = nl + 1
         ijn2 = kln1*(nk+1)
         ijn1 = ijn2*(nj+1)
         inc1 = ijn1*(ni+1)
c     if(mod(inc1,4).eq.0)inc1=inc1+1
         if (ni.lt.nj) then
            is = ni
            ni = nj
            nj = is
            ij1 = ijn2
            ij2 = ijn1
            xc = xj
            yc = yj
            zc = zj
            dxij = xj - xi
            dyij = yj - yi
            dzij = zj - zi
         else
            ij1 = ijn1
            ij2 = ijn2
            xc = xi
            yc = yi
            zc = zi
            dxij = xi - xj
            dyij = yi - yj
            dzij = zi - zj
         end if
         if (nk.lt.nl) then
            is = nl
            nl = nk
            nk = is
            kl1 = kln2
            kl2 = kln1
            xd = xl
            yd = yl
            zd = zl
            dxkl = xl - xk
            dykl = yl - yk
            dzkl = zl - zk
         else
            xd = xk
            yd = yk
            zd = zk
            dxkl = xk - xl
            dykl = yk - yl
            dzkl = zk - zl
            kl1 = kln1
            kl2 = kln2
         end if
         nmax = ni + nj
         mmax = nk + nl
         max = nmax + 1
         do 20 i = 1 , max
            n = i - 1
            if (n.le.ni) in(i) = ij1*n + 1
            if (n.gt.ni) in(i) = ij1*ni + ij2*(n-ni) + 1
 20      continue
         max = mmax + 1
         do 30 k = 1 , max
            n = k - 1
            if (n.le.nk) kn(k) = kl1*n
            if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 30      continue
c
c     indexing
c
         call indxa_dft(ijx,ijy,ijz,ijd,mini,maxi,minj,maxj,oiandj,
     +        ijn1, ijn2,1)
         call indxa_dft(klx,kly,klz,kld,mink,maxk,minl,maxl,okandl,
     +        kln1, kln2,0)

         nn = 0
         ijkld = 0
         do 50 i = 1 , ijd
            do 40 k = 1 , kld
               nn = nn + 1
               if (.not.(noform(nn))) then
                  ijkld = ijkld + 1
                  iqq(ijkld+ixi-1) = ijx(i) + klx(k)
                  iqq(ijkld+iyi-1) = ijy(i) + kly(k)
                  iqq(ijkld+izi-1) = ijz(i) + klz(k)
               end if
 40         continue
 50      continue
c
         do 60 n = 1 , nij
            axak(n) = aa(n)*(x1(n)-xd)
            ayak(n) = aa(n)*(y1(n)-yd)
            azak(n) = aa(n)*(z1(n)-zd)
            axai(n) = aa(n)*(x1(n)-xc)
            ayai(n) = aa(n)*(y1(n)-yc)
            azai(n) = aa(n)*(z1(n)-zc)

            if(odbg)write(6,*)'axak etc',n,axak(n),ayak(n),azak(n)
            if(odbg)write(6,*)'axai etc',n,axai(n),ayai(n),azai(n)

 60      continue

c
c
         trduij = lit.ge.3 .or. ljt.ge.3
         trdukl = lkt.ge.3 .or. llt.ge.3
         spkl = (mink.eq.1 .and. maxk.eq.4) .or.
     +          (minl.eq.1 .and. maxl.eq.4)
         spij = (mini.eq.1 .and. maxi.eq.4) .or.
     +          (minj.eq.1 .and. maxj.eq.4)
         sptru = spkl .or. spij
c
c Now dgenrl uses preallocated memory, max size defined 
c in jkder only pointers ic2 upwards are allocated here
c vect factor of 32 is hardwired
c
c NB there is room above ic7 assuming npass=3, lendd,=lnddm
c
          ncmmm = 32
          ncmmm = (ncmmm/nroots)*nroots

cccccc old memory algorithm cccccccccccccc
c
c     integrals stored at (ic7+1)
c     subsidiary integrals from ic1 onwards
c
cc         ic1 = ic7 + 1 + npass*lendd*3
cc         ncmmm = (nmaxly-ic1-1)/(inc1*6)
cc         ncnnn = ncmax - 1
cc         if (ncmmm.gt.ncnnn) ncmmm = ncnnn
cc         ncmmm = (ncmmm/nroots)*nroots
cc         if (ncmmm.lt.nroots) call caserr('insufficient core in dgenrl')
cc
ccccccc end old memory algorithm cccccccccccc
c
         ic1 = ic7 + 1 + npass*lendd*3
         ic2 = inc1*ncmmm + ic1
         ic3 = inc1*ncmmm + ic2
         ic4 = inc1*ncmmm + ic3
         ic5 = inc1*ncmmm + ic4
         ic6 = inc1*ncmmm + ic5
c
         inc = npass*lendd*3
         call vclr(qq(ic7+1),1,inc)
         ncontr = 0
c
         maxlg = ngd
         do 170 kg = 1 , ngc
            ak = cgg(kg)
            brrk = ak*rrk
            akxk = ak*xk
            akyk = ak*yk
            akzk = ak*zk
            csk = csc(kg)*pi252
            cpk = cpc(kg)*pi252
            cdk = cdc(kg)*pi252
            cfk = cfc(kg)*pi252
            cgk = cgc(kg)*pi252
c
c     ----- l primitive
c
            if (okandl) maxlg = kg
            do 160 lg = 1 , maxlg
               al = dg(lg)
               b = ak + al
               binv = one/b
               bbrrk = al*brrk*binv
               if ((bbrrk+rsmall).le.tol1) then
                  exkl = dexp(-bbrrk)
                  csl = csd(lg)*binv
                  cpl = cpd(lg)*binv
                  cdl = cdd(lg)*binv
                  cfl = cfd(lg)*binv
                  cgl = cgd(lg)*binv
                  xb = (akxk+al*xl)*binv
                  yb = (akyk+al*yl)*binv
                  zb = (akzk+al*zl)*binv
                  bxbk = b*(xb-xd)
                  bybk = b*(yb-yd)
                  bzbk = b*(zb-zd)
                  bxbi = b*(xb-xc)
                  bybi = b*(yb-yc)
                  bzbi = b*(zb-zc)
c
c     ----- density factor
c
                  double = okandl .and. kg.ne.lg

                  call denfac_dft(dkl,csk,cpk,cdk,cfk,cgk,
     +                        csl,cpl,cdl,cfl,cgl,
     +                        mink,maxk,minl,maxl,okandl,double)

                  dkld = dkl(1)
                  if (sptru) then
                     do 70 k = 1 , kld
                        dkl(k) = dkl(k)/dkld
 70                  continue
                  end if
                  dkld = dkld*exkl
                  if (dabs(dkld*abmax).ge.tol3) then
c
c     ----- pair of i,j primitives
c
                     do 80 n = 1 , nij
                        abv(n) = aa(n)*b
                        aandbv(n) = aa(n) + b
                        expev(n) = exij(n)/dsqrt(aa(n)+b)
                        rhov(n) = abv(n)/aandbv(n)
                        xxv(n) = rhov(n)
     +                           *((x1(n)-xb)**2+(y1(n)-yb)**2+(z1(n)
     +                           -zb)**2)
                        c1xv(n) = bxbk + axak(n)
                        c2xv(n) = bxbk*aa(n)
                        c3xv(n) = bxbi + axai(n)
                        c4xv(n) = b*axai(n)
                        c1yv(n) = bybk + ayak(n)
                        c2yv(n) = bybk*aa(n)
                        c3yv(n) = bybi + ayai(n)
                        c4yv(n) = b*ayai(n)
                        c1zv(n) = bzbk + azak(n)
                        c2zv(n) = bzbk*aa(n)
                        c3zv(n) = bzbi + azai(n)
                        c4zv(n) = b*azai(n)
 80                  continue
c
                     n = 0
                     nn = 0
                     jgmax = ngb
                     do 150 ig = 1 , nga
                        ai = ag(ig)
                        if (oiandj) jgmax = ig
                        do 140 jg = 1 , jgmax
                           n = n + 1
                           if ((bbrrk+r(n)).lt.tol2) then
                              aj = bg(jg)
                              dijd = dd(nn+1)
                              if (sptru) then
                                 dddd = one/dijd
                                 do 90 i = 1 , ijd
                                    dij(i) = dd(ijden(i)+nn)*dddd
 90                              continue
                              end if
                              expe = dkld*dijd*expev(n)
                              if (dabs(expe*abmax).ge.tol4) then
                                 pp = xxv(n)
c
c     ----- roots and weights for quadrature
c
                                 if (nroots.le.3) call rt123_dft
                                 if (nroots.eq.4) call roots4_dft
                                 if (nroots.eq.5) call roots5_dft
                                 if (nroots.gt.5) call rootss_dft
c
c     compute two-electron  integrals for each root
c
                                 nnn0 = ncontr
                                 do 100 m = 1 , nroots
                                    ncontr = ncontr + 1
                                    u2 = u(m)*rhov(n)
                                    f00(ncontr) = expe*w(m)
                                    dum2 = 0.5d0/(abv(n)+u2*aandbv(n))
                                    dum = dum2 + dum2
                                    bp01(ncontr) = (aa(n)+u2)*dum2
                                    b00(ncontr) = u2*dum2
                                    b10(ncontr) = (b+u2)*dum2
                                    xcp00(ncontr) = (u2*c1xv(n)+c2xv(n))
     +                                 *dum
                                    xc00(ncontr) = (u2*c3xv(n)+c4xv(n))
     +                                 *dum
                                    ycp00(ncontr) = (u2*c1yv(n)+c2yv(n))
     +                                 *dum
                                    yc00(ncontr) = (u2*c3yv(n)+c4yv(n))
     +                                 *dum
                                    zcp00(ncontr) = (u2*c1zv(n)+c2zv(n))
     +                                 *dum
                                    zc00(ncontr) = (u2*c3zv(n)+c4zv(n))
     +                                 *dum

                                    aei(ncontr) = ai
                                    aej(ncontr) = aj
                                    aek(ncontr) = ak
                                    ael(ncontr) = al

                            if(odbg)write(6,*)'aei etc',m, aei(ncontr),
     &                    aej(ncontr), aek(ncontr), ael(ncontr)
                            if(odbg)write(6,*)'xcp00',m, xcp00(ncontr),
     &                    ycp00(ncontr), zcp00(ncontr)
                            if(odbg)write(6,*)'xc00',m, xc00(ncontr),
     &                    yc00(ncontr), zc00(ncontr)


 100                             continue

                                 if (sptru) then
                                    ncontr = nnn0
                                    do 130 m = 1 , nroots
                                       ncontr = ncontr + 1
                                       do 110 iii = 1 , ijd
                                         ddij(ncontr,iii) = dij(iii)
 110                                   continue
                                       do 120 iii = 1 , kld
                                         ddkl(ncontr,iii) = dkl(iii)
 120                                   continue
 130                                continue
                                 end if

c
c
c----------------------------------------------
c    defer assembly stage until loop lengths
c    are long enough to vectorise effectively
c----------------------------------------------
c
                                 if (ncontr.ge.ncmmm) then
c     ----- form (i,j//k,l) integrals
c
                                 call dxyz_dft(qq(ic1),qq(ic2),qq(ic3),
     +                                   ncmmm)
                                    ioff = 0
                                    if (.not.(oskip(1))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aei,ljt,lkt,llt,lit,
     +                                    ijn2,kln1,kln2,ijn1,ncmmm)
                                       if (sptru) then
                                     call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                   qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                   iqq(izi),ncmmm)
                                       else
                                     call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(2))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aej,lit,lkt,llt,ljt,
     +                                    ijn1,kln1,kln2,ijn2,ncmmm)
                                       if (sptru) then
                                      call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                     qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(3))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aek,lit,ljt,llt,lkt,
     +                                    ijn1,ijn2,kln2,kln1,ncmmm)
                                       if (sptru) then
                                     call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                    qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(4))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),ael,lit,ljt,lkt,llt,
     +                                    ijn1,ijn2,kln1,kln2,ncmmm)
                                       if (sptru) then
                                      call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                    qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                    end if
                                    ncontr = 0
                                 end if
                              end if
                           end if
c
c
c     end of loops over primitives
c
                           nn = nn + 4
 140                    continue
 150                 continue
                  end if
               end if
 160        continue
 170     continue
c
c    tidy up any bits left unassembled
c
c
         if (ncontr.ne.0) then
            call dxyz_dft(qq(ic1),qq(ic2),qq(ic3),ncmmm)
            ioff = 0
            if (.not.(oskip(1))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,
     +                    ijn1,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(2))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,
     +                    ijn2,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(3))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,
     +                    kln1,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(4))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,
     +                    kln2,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                 iqq(izi),ncmmm)
               end if
            end if
            ncontr = 0
         end if
c
c ---------------------------------------------
c    fiddle about if first two centres are same
c ---------------------------------------------
c
         if (natomd(1).eq.natomd(2)) then
            inc = lendd*3
            ii = ic7 + inc
            do 180 i = 1 , inc
               qq(ic7+i) = qq(ic7+i) + qq(ii+i)
 180        continue
            iii = ii + inc
            do 190 i = 1 , inc
               qq(ii+i) = qq(iii+i)
 190        continue
            npass = npass - 1
            natomd(2) = natomd(3)
            natomd(3) = natomd(4)
            natomd(4) = 0
         end if
c
c    insert proper normalisation factors in d-functions
c    present (f- and g-functions also)
c

         if (trduij .or. trdukl) then
            if (onorm) then
               call dnorm_dft(qq(ic7+1),noform)
            end if
         end if
      end if
c
c    multiply by density matrix elements
c
      call grdcon_dft(qq,noform,gmax)
      return
      end
      subroutine dnorm_dft(qq,noform)

      implicit none
c
c arguments
c
      real*8 qq(*)
      logical noform(*)

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
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c

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

      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c

c
c  local variables
c
      integer max, i, j, k, l, n, nn, nnn, npp
      real*8 dum1, dum2, d1

      real*8 one
      data one/1.0d0/

      n = 0
      max = maxj
      dum1 = one
      do 30 i = mini , maxi
         if (i.eq.8)  dum1 = root3
         if (i.eq.14) dum1 = root5
         if (i.eq.20) dum1 = dum1*root3
         if (i.eq.24) dum1 = root7
         if (i.eq.30) dum1 = dum1*root53
         if (i.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (oiandj) max = i
         do 20 j = minj , max
            if (j.eq.8)  dum2 = dum2*root3
            if (j.eq.14) dum2 = dum2*root5
            if (j.eq.20) dum2 = dum2*root3
            if (j.eq.24) dum2 = dum2*root7
            if (j.eq.30) dum2 = dum2*root53
            if (j.eq.33) dum2 = dum2*root3
            n = n + 1
            dij(n) = dum2
 20      continue
 30   continue
      n = 0
      dum1 = one
      max = maxl
      do 50 k = mink , maxk
         if (k.eq.8)  dum1 = root3
         if (k.eq.14) dum1 = root5
         if (k.eq.20) dum1 = dum1*root3
         if (k.eq.24) dum1 = root7
         if (k.eq.30) dum1 = dum1*root53
         if (k.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (okandl) max = k
         do 40 l = minl , max
            if (l.eq.8)  dum2 = dum2*root3
            if (l.eq.14) dum2 = dum2*root5
            if (l.eq.20) dum2 = dum2*root3
            if (l.eq.24) dum2 = dum2*root7
            if (l.eq.30) dum2 = dum2*root53
            if (l.eq.33) dum2 = dum2*root3
            n = n + 1
            dkl(n) = dum2
 40      continue
 50   continue
      nn = 0
      nnn = 0
      do 80 i = 1 , ijd
         d1 = dij(i)
         do 70 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               nnn = nnn + 1
               ioff = 0
               do 60 npp = 1 , npass*3
                  qq(nnn+ioff) = qq(nnn+ioff)*d1*dkl(k)
                  ioff = ioff + lendd
 60            continue
            end if
 70      continue
 80   continue
      return
      end
c
c
c  generate  r, x1, y1, z1, dd in d2escr
c
c
      subroutine dprim_dft

      implicit none

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

      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      real*8 ag,csa,cpa,cda,cfa,cga,bg,csb,cpb,cdb,cfb,cgb
      real*8 cgg,csc,cpc,cdc,cfc,cgc,dg,csd,cpd,cdd,cfd,cgd
      real*8 xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk
      real*8 exij,rsmall
      integer nga,ngb,ngc,ngd
      common/dshlnf/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +              cfa(mxprms),cga(mxprms),
     +              bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +              cfb(mxprms),cgb(mxprms),
     +             cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +              cfc(mxprms),cgc(mxprms),
     +              dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +              cfd(mxprms),cgd(mxprms),
     +              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     +              nga,ngb,ngc,ngd,exij(mxprms*mxprms),rsmall
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
c
      logical odbg
      common/dbgdbg/odbg


      real*8 atmp, aj, ai, ainv, dum, azi, dum1, dum2
      real*8 csi, cpi, cdi, cfi, cgi
      real*8 csj, cpj, cdj, cfj, cgj
      integer n, i, j, jb, max
      real*8 one
      real*8 axi, ayi, arri
      integer nm, nn, ia, jbmax

      data one/1.0d0/
      max = maxj
      n = 0
      nn = 0
      do 70 i = mini , maxi
         go to (20,20,30,30,
     +          20,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          30,30,30,30,30) , i
 20      nm = nn
 30      nn = nm
         if (oiandj) max = i
         do 60 j = minj , max
            go to (40,40,50,50,
     +             40,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             50,50,50,50,50) , j
 40         nn = nn + 1
 50         n = n + 1
            ijden(n) = nn
 60      continue
 70   continue
c     ----- i primitive
      nij = 0
      jbmax = ngb
      do 230 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*xi
         ayi = ai*yi
         azi = ai*zi
         csi = csa(ia)
         cpi = cpa(ia)
         cdi = cda(ia)
         cfi = cfa(ia) 
         cgi = cga(ia) 
c     ----- j primitive
         if (oiandj) jbmax = ia
         do 220 jb = 1 , jbmax
            aj = bg(jb)
            atmp = ai + aj
            ainv = one/atmp
            dum = aj*arri*ainv
            csj = csb(jb)*ainv
            cpj = cpb(jb)*ainv
            cdj = cdb(jb)*ainv
            cfj = cfb(jb)*ainv
            cgj = cgb(jb)*ainv
            nm = (nij+nij) + (nij+nij)
            nn = nm
            nij = nij + 1
            r(nij) = dum
            aa(nij) = atmp
            x1(nij) = (axi+aj*xj)*ainv
            y1(nij) = (ayi+aj*yj)*ainv
            z1(nij) = (azi+aj*zj)*ainv
c     ----- density factor
            do 190 i = mini , maxi
               if (oiandj) max = i
               go to (80,90,190,190,
     +               100,190,190,190,190,190,
     +               110,190,190,190,190,190,190,190,190,190,
     +               115,190,190,190,190,190,190,190,190,190,
     +               190,190,190,190,190) , i
 80            dum1 = csi
               go to 120
 90            dum1 = cpi
               go to 120
 100           dum1 = cdi
               go to 120
 110           dum1 = cfi
               go to 120
 115           dum1 = cgi
 120           do 180 j = minj , max
                  go to (130,140,180,180,
     +                   150,180,180,180,180,180,
     +                   160,180,180,180,180,180,180,180,180,180,
     +                   165,180,180,180,180,180,180,180,180,180,
     +                   180,180,180,180,180) , j
 130              dum2 = dum1*csj
                  go to 170
 140              dum2 = dum1*cpj
                  go to 170
 150              dum2 = dum1*cdj
                  go to 170
 160              dum2 = dum1*cfj
                  go to 170
 165              dum2 = dum1*cgj
 170              nn = nn + 1
                  dd(nn) = dum2

c            if(odbg)write(6,*)'dprim dd',nn,dd(nn)

 180           continue
 190        continue


            if (.not.oiandj) go to 220
            if (ia.eq.jb) go to 220
            go to (210,200,210,210,210) , lit
 200        if (mini.ne.2) then
               dd(nm+2) = dd(nm+2) + csi*cpj
               dd(nm+3) = dd(nm+3) + dd(nm+3)
            end if
 210        dd(nm+1) = dd(nm+1) + dd(nm+1)
 220     continue
 230  continue

c      if(odbg)write(6,*)'dprim r',nij,(r(n),n=1,nij)
c      if(odbg)write(6,*)'dprim x1',(x1(n),n=1,nij)
c      if(odbg)write(6,*)'dprim y1',(y1(n),n=1,nij)
c      if(odbg)write(6,*)'dprim z1',(z1(n),n=1,nij)


      if (nij.eq.0) return
      rsmall = r(1)
      do 240 n = 1 , nij
         exij(n) = dexp(-r(n))
 240  continue
      do 250 n = 1 , nij
         if (rsmall.gt.r(n)) rsmall = r(n)
 250  continue
      if (rsmall.ge.tol1) nij = 0
      return
      end
      subroutine dshell_dft(nelec,ish,jsh,ksh,lsh,
     +     basi, basj, bask, basl, ncentr)

c      implicit real*8  (a-h,o-z)
      implicit none
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

      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
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

      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 ag,csa,cpa,cda,cfa,cga,bg,csb,cpb,cdb,cfb,cgb
      real*8 cgg,csc,cpc,cdc,cfc,cgc,dg,csd,cpd,cdd,cfd,cgd
      real*8 xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk
      real*8 exij,rsmall
      integer nga,ngb,ngc,ngd
      common/dshlnf/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +              cfa(mxprms),cga(mxprms),
     +              bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +              cfb(mxprms),cgb(mxprms),
     +             cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +              cfc(mxprms),cgc(mxprms),
     +              dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +              cfd(mxprms),cgd(mxprms),
     +              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     +              nga,ngb,ngc,ngd,exij(mxprms*mxprms),rsmall
c

      logical odbg
      common/dbgdbg/odbg

c
c arguments
c
      integer nelec
      integer ish, jsh, ksh, lsh
      integer basi, basj, bask, basl, ncentr
c
c local
c
      integer i, j, k, l
      integer i1, j1, k1, l1
      integer i2, j2, k2, l2
      integer ittt, max

      if (nelec.eq.2) then

         okandl = ksh.eq.lsh .and. ncentr .eq. 4
         osame = ish.eq.ksh .and. jsh.eq.lsh .and. ncentr .eq. 4

         k = katom(bask,ksh)
         xk = c(1,k)
         yk = c(2,k)
         zk = c(3,k)
         k1 = kstart(bask,ksh)
         k2 = k1 + kng(bask,ksh) - 1
         lkt = ktype(bask,ksh)
         mink = kmin(bask,ksh)
         maxk = kmax(bask,ksh)
         lock = kloc(bask,ksh) - mink
         ngc = 0
         if(odbg)write(6,*)'dshell k loop',k,lock
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex_m(bask,k)
            csc(ngc) = cs(bask,k)
            cpc(ngc) = cp(bask,k)
            cdc(ngc) = cd(bask,k)
            cfc(ngc) = cf(bask,k)
            cgc(ngc) = cg(bask,k)
            if(odbg)write(6,*)ngc,cgg(ngc),csc(ngc),cpc(ngc),cdc(ngc)
 20      continue

         if(basl .lt. 0) then
c
c  Dummy l centr
c
            l = k
            xl = c(1,l)
            yl = c(2,l)
            zl = c(3,l)
            l1 = 1
            l2 = 1
            llt = 1
            minl = 1
            maxl = 1
            locl = 0
            ngd = 0
            if(odbg)write(6,*)'dshell l loop',l,locl
            do 30 l = l1 , l2
               ngd = ngd + 1
               dg(ngd)  = 0.0d0
               csd(ngd) = 1.0d0
               cpd(ngd) = 1.0d0
               cdd(ngd) = 1.0d0
               cfd(ngd) = 1.0d0
               cgd(ngd) = 1.0d0
               if(odbg)write(6,*)ngc,dg(ngc),csd(ngc),cpd(ngc),cdd(ngc)
 30         continue
            rrk = 0.0d0
         else
            l = katom(basl,lsh)
            xl = c(1,l)
            yl = c(2,l)
            zl = c(3,l)
            l1 = kstart(basl,lsh)
            l2 = l1 + kng(basl,lsh) - 1
            llt = ktype(basl,lsh)
            minl = kmin(basl,lsh)
            maxl = kmax(basl,lsh)
            locl = kloc(basl,lsh) - minl
            ngd = 0
            if(odbg)write(6,*)'dshell l loop',locl
            do l = l1 , l2
               ngd = ngd + 1
               dg(ngd) = ex_m(basl,l)
               csd(ngd) = cs(basl,l)
               cpd(ngd) = cp(basl,l)
               cdd(ngd) = cd(basl,l)
               cfd(ngd) = cf(basl,l)
               cgd(ngd) = cg(basl,l)
               if(odbg)write(6,*)ngd,dg(ngd),csd(ngd),cpd(ngd),cdd(ngd)
            enddo
            rrk = ((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
         endif
         nroots = (lit+ljt+lkt+llt-1)/2

c
c     determine various offsets and indexing arrays
c
         inc2 = 1
         inc3 = inc2*(maxl-minl+1)
         inc4 = inc3*(maxk-mink+1)
         inc5 = inc4*(maxj-minj+1)
         lendd = inc5*(maxi-mini+1)

c         write(6,*)'test1',nelec,inc1,inc2,inc3,inc4,inc5,lendd
c         write(6,*)mini,maxi,minj,maxj

         if (mod(lendd,4).eq.0) lendd = lendd + 1
         ijd = 0
         max = maxj
         do 50 i = mini , maxi
            if (oiandj) max = i
            ittt = inc5*(i-mini) + 1
            do 40 j = minj , max
               ijd = ijd + 1
               ijgt(ijd) = ittt
               ittt = ittt + inc4
 40         continue
 50      continue
         kld = 0
         max = maxl
         do 70 k = mink , maxk
            if (okandl) max = k
            ittt = inc3*(k-mink)
            do 60 l = minl , max
               kld = kld + 1
               klgt(kld) = ittt
               ittt = ittt + inc2
 60         continue
 70      continue
         ijkld = ijd*kld
         ixi = lendd + 1
         iyi = ixi + lendd
         izi = iyi + lendd
         ioff = 0

c         if(odbg)write(6,*)ijkld,ixi,iyi,izi
c         if(odbg)write(6,*)ijgt
c         if(odbg)write(6,*)klgt

         return
      else
         oiandj = ish.eq.jsh .and. ncentr .gt. 2
         i = katom(basi,ish)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(basi,ish)
         i2 = i1 + kng(basi,ish) - 1
         lit = ktype(basi,ish)
         mini = kmin(basi,ish)
         maxi = kmax(basi,ish)
         loci = kloc(basi,ish) - mini
         nga = 0

         if(odbg)write(6,*)'dshell i loop',i, loci

         do 80 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex_m(basi,i)
            csa(nga) = cs(basi,i)
            cpa(nga) = cp(basi,i)
            cda(nga) = cd(basi,i)
            cfa(nga) = cf(basi,i)
            cga(nga) = cg(basi,i)

            if(odbg)write(6,*)nga,ag(nga),csa(nga),cpa(nga),cda(nga)

 80      continue

         if(basj .lt. 0)then
c
c Dummy j centre
c
            j = i
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = 1
            j2 = 1
            ljt = 1
            minj = 1
            maxj = 1
            locj = 0
            ngb = 0

            if(odbg)write(6,*)'dshell j loop',j, locj

            do 90 j = j1 , j2
               ngb = ngb + 1
               bg(ngb)  = 0.0d0
               csb(ngb) = 1.0d0
               cpb(ngb) = 1.0d0
               cdb(ngb) = 1.0d0
               cfb(ngb) = 1.0d0
               cgb(ngb) = 1.0d0

               if(odbg)write(6,*)ngb,bg(ngb),csb(ngb),cpb(ngb),cdb(ngb)

 90         continue
            rri = 0.0d0
         else
            j = katom(basj,jsh)
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = kstart(basj,jsh)
            j2 = j1 + kng(basj,jsh) - 1
            ljt = ktype(basj,jsh)
            minj = kmin(basj,jsh)
            maxj = kmax(basj,jsh)
            locj = kloc(basj,jsh) - minj
            ngb = 0
            if(odbg)write(6,*)'dshell j loop',j, locj

            do j = j1 , j2
               ngb = ngb + 1
               bg(ngb)  = ex_m(basj,j)
               csb(ngb) = cs(basj,j)
               cpb(ngb) = cp(basj,j)
               cdb(ngb) = cd(basj,j)
               cfb(ngb) = cf(basj,j)
               cgb(ngb) = cg(basj,j)
               if(odbg)write(6,*)ngb,bg(ngb),csb(ngb),cpb(ngb),cdb(ngb)
            enddo
            rri = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
         endif
         return
      end if
      end
      subroutine dxyz_dft(x,y,z,ncdim)
c      implicit real*8  (a-h,o-z)
      implicit none

      integer ncmax
      parameter (ncmax=65)
      real*8 x, y, z
      dimension x(*),y(*),z(*)
      integer ncdim

      logical n0,n1,m0,m1

      real*8 bp01,b00,b10,xcp00, xc00,ycp00,yc00,zcp00,zc00,f00,
     &     dxij,dyij,dzij,dxkl,dykl,dzkl
      integer iorg,korg,nimax,njmax,nkmax,nlmax,nmax,mmax,
     &     ij1x,ij2x,kl1x,kl2x

      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     1   xc00(ncmax),
     2ycp00(ncmax),yc00(ncmax),zcp00(ncmax),zc00(ncmax),f00(ncmax),dxij,
     3dyij,dzij,dxkl,dykl,dzkl,iorg(12),korg(12),
     4nimax,njmax,nkmax,nlmax,nmax,mmax,ij1x,ij2x,kl1x,kl2x

      real*8 ca, cb
      common/small/ca(ncmax),cb(ncmax)

c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c

      integer ni, nj, nk, nl
      integer i1, i2, i3, i4, i5, k2, k4
      integer ia, ib, ic, ink, nc, ij2
      integer i, ij1, k, n, m, km, k3, min
      integer kl1, kl2
      real*8 zero, one
c
      dimension i(12),k(12)
c
      data zero,one /0.0d+00,1.0d+00/
c
      do 20 n = 1 , nmax + 1
         i(n) = (iorg(n)-1)*ncdim + 1
 20   continue
      do 30 n = 1 , mmax + 1
         k(n) = korg(n)*ncdim
 30   continue
      ij1 = ij1x*ncdim
      ij2 = ij2x*ncdim
      kl1 = kl1x*ncdim
      kl2 = kl2x*ncdim
      ink = 1
c
      n0 = nmax.eq.0
      n1 = nmax.le.1
      m0 = mmax.eq.0
      m1 = mmax.le.1
      if (n0) then
         i1 = i(1)
         ia = 0
         do 40 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 40      continue
         if (m0) go to 670
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
         do 50 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 50      continue
         if (m1) go to 670
c     ----- i(0,nk) -----
c
         do 60 nc = 1 , ncontr
            ca(nc) = zero
 60      continue
         i3 = i1
         i4 = i1 + k2
         do 90 nk = 2 , mmax
c
            do 70 nc = 1 , ncontr
               ca(nc) = ca(nc) + bp01(nc)
 70         continue
            i5 = i1 + k(nk+1)
            ia = 0
            do 80 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
               ia = ia + ink
 80         continue
            i3 = i4
            i4 = i5
 90      continue
         if (nlmax.eq.0) go to 670
c     ----- i(0,0,nk,nl) -----
         i5 = k(mmax+1)
         min = nkmax
         go to 610
      else if (m0) then
         i1 = i(1)
         ia = 0
         do 100 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 100     continue
         if (n0) go to 600
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
         do 110 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 110     continue
         if (n1) go to 600
c     ----- i(ni,0) -----
c
         do 120 nc = 1 , ncontr
            ca(nc) = zero
 120     continue
         i3 = i1
         i4 = i2
         do 150 ni = 2 , nmax
            do 130 nc = 1 , ncontr
               ca(nc) = ca(nc) + b10(nc)
 130        continue
            i5 = i(ni+1)
            ia = 0
            do 140 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
               ia = ia + ink
 140        continue
            i3 = i4
            i4 = i5
 150     continue
         if (njmax.eq.0) go to 600
c     ----- i(ni,nj,0,0) -----
         i5 = i(nmax+1)
         min = nimax
         go to 540
      else
c     ----- i(0,0) -----
         i1 = i(1)
         ia = 0
c
         do 160 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 160     continue
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
c
         do 170 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 170     continue
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
c
         do 180 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 180     continue
c     ----- i(1,1) -----
         i3 = i2 + k2
         ia = 0
         do 190 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)*x(i2+ia) + b00(nc)
            y(i3+ia) = ycp00(nc)*y(i2+ia) + b00(nc)
            z(i3+ia) = zcp00(nc)*z(i2+ia) + b00(nc)*f00(nc)
            ia = ia + ink
 190     continue
         if (.not.(n1)) then
            do 200 nc = 1 , ncontr
               cb(nc) = b00(nc)
               ca(nc) = zero
 200        continue
            i3 = i1
            i4 = i2
            do 250 n = 2 , nmax
c     ----- i(n,0) -----
               i5 = i(n+1)
               do 210 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b10(nc)
 210           continue
               ia = 0
               do 220 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
                  ia = ia + ink
 220           continue
               do 230 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 230           continue
c     ----- i(n,1) -----
               i3 = i5 + k2
               ia = 0
               do 240 nc = 1 , ncontr
                  x(i3+ia) = xcp00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = ycp00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zcp00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 240           continue
               i3 = i4
               i4 = i5
 250        continue
         end if
         if (.not.(m1)) then
            do 260 nc = 1 , ncontr
               ca(nc) = zero
               cb(nc) = b00(nc)
 260        continue
            i3 = i1
            i4 = i1 + k2
            do 310 m = 2 , mmax
               do 270 nc = 1 , ncontr
                  ca(nc) = ca(nc) + bp01(nc)
 270           continue
c     ----- i(0,m) -----
               i5 = i1 + k(m+1)
               ia = 0
               do 280 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
                  ia = ia + ink
 280           continue
               do 290 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 290           continue
c     ----- i(1,m) -----
               i3 = i2 + k(m+1)
               ia = 0
               do 300 nc = 1 , ncontr
                  x(i3+ia) = xc00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = yc00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zc00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 300           continue
               i3 = i4
               i4 = i5
 310        continue
         end if
         if (.not.(n1 .or. m1)) then
c     ----- i(n,m) -----
c
            do 320 nc = 1 , ncontr
               ca(nc) = b00(nc)
 320        continue
            k3 = k2
            do 370 m = 2 , mmax
               k4 = k(m+1)
               do 330 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b00(nc)
                  cb(nc) = b10(nc)
 330           continue
               i3 = i1
               i4 = i2
               do 360 n = 2 , nmax
                  i5 = i(n+1)
                  ia = 0
                  do 340 nc = 1 , ncontr
                     x(i5+k4+ia) = cb(nc)*x(i3+k4+ia) + xc00(nc)
     +                             *x(i4+k4+ia) + ca(nc)*x(i4+k3+ia)
                     y(i5+k4+ia) = cb(nc)*y(i3+k4+ia) + yc00(nc)
     +                             *y(i4+k4+ia) + ca(nc)*y(i4+k3+ia)
                     z(i5+k4+ia) = cb(nc)*z(i3+k4+ia) + zc00(nc)
     +                             *z(i4+k4+ia) + ca(nc)*z(i4+k3+ia)
                     ia = ia + ink
 340              continue
                  do 350 nc = 1 , ncontr
                     cb(nc) = cb(nc) + b10(nc)
 350              continue
                  i3 = i4
                  i4 = i5
 360           continue
               k3 = k4
 370        continue
         end if
         if (njmax.eq.0) go to 450
c     ----- i(ni,nj,m) -----
         m = 0
         i5 = i(nmax+1)
      end if
 380  min = nimax
      km = k(m+1)
 390  n = nmax
      i3 = i5 + km
 400  i4 = i(n) + km
      ia = 0
      do 410 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 410  continue
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 400
      min = min + 1
      if (min.lt.nmax) go to 390
      if (nimax.ne.0) then
         i3 = ij2 + km + i1
         do 440 nj = 1 , njmax
            i4 = i3
            do 430 ni = 1 , nimax
               ia = 0
               do 420 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 420           continue
               i4 = i4 + ij1
 430        continue
            i3 = i3 + ij2
 440     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 380
 450  if (nlmax.eq.0) go to 530
c
c     ----- i(ni,nj,nk,nl) -----
c
      i5 = k(mmax+1)
      ia = i1
      ni = 0
 460  nj = 0
      ib = ia
      min = nkmax
 470  m = mmax
      i3 = ib + i5
 480  i4 = ib + k(m)
      ic = 0
      do 490 nc = 1 , ncontr
         x(i3+ic) = x(i3+ic) + dxkl*x(i4+ic)
         y(i3+ic) = y(i3+ic) + dykl*y(i4+ic)
         z(i3+ic) = z(i3+ic) + dzkl*z(i4+ic)
         ic = ic + ink
 490  continue
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 480
      min = min + 1
      if (min.lt.mmax) go to 470
      if (nkmax.ne.0) then
         i3 = ib + kl2
         do 520 nl = 1 , nlmax
            i4 = i3
            do 510 nk = 1 , nkmax
               ic = 0
               do 500 nc = 1 , ncontr
                  x(i4+ic) = x(i4+ic+kl1-kl2) + dxkl*x(i4+ic-kl2)
                  y(i4+ic) = y(i4+ic+kl1-kl2) + dykl*y(i4+ic-kl2)
                  z(i4+ic) = z(i4+ic+kl1-kl2) + dzkl*z(i4+ic-kl2)
                  ic = ic + ink
 500           continue
               i4 = i4 + kl1
 510        continue
            i3 = i3 + kl2
 520     continue
      end if
      nj = nj + 1
      ib = ib + ij2
      if (nj.le.njmax) then
         min = nkmax
         go to 470
      else
         ni = ni + 1
         ia = ia + ij1
         if (ni.le.nimax) go to 460
      end if
 530  return
 540  ni = nmax
      i3 = i5
 550  i4 = i(ni)
      ia = 0
      do 560 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 560  continue
      i3 = i4
      ni = ni - 1
      if (ni.gt.min) go to 550
      min = min + 1
      if (min.lt.nmax) go to 540
      if (nimax.ne.0) then
         i3 = ij2 + i1
         do 590 nj = 1 , njmax
            i4 = i3
            do 580 ni = 1 , nimax
               ia = 0
               do 570 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 570           continue
               i4 = i4 + ij1
 580        continue
            i3 = i3 + ij2
 590     continue
      end if
 600  return
 610  nk = mmax
      i3 = i1 + i5
 620  i4 = i1 + k(nk)
      ia = 0
      do 630 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxkl*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dykl*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzkl*z(i4+ia)
         ia = ia + ink
 630  continue
      i3 = i4
      nk = nk - 1
      if (nk.gt.min) go to 620
      min = min + 1
      if (min.lt.mmax) go to 610
      if (nkmax.ne.0) then
         i3 = i1 + kl2
         do 660 nl = 1 , nlmax
            i4 = i3
            do 650 nk = 1 , nkmax
               ia = 0
               do 640 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+kl1-kl2) + dxkl*x(i4+ia-kl2)
                  y(i4+ia) = y(i4+ia+kl1-kl2) + dykl*y(i4+ia-kl2)
                  z(i4+ia) = z(i4+ia+kl1-kl2) + dzkl*z(i4+ia-kl2)
                  ia = ia + ink
 640           continue
               i4 = i4 + kl1
 650        continue
            i3 = i3 + kl2
 660     continue
      end if
 670  return
      end
      subroutine subsd_dft(x,y,z,xd,yd,zd,
     *a,m1,m2,m3,m4,i1,i2,k1,k2,ncdim)

c      implicit real*8  (a-h,o-z)
      implicit none

      integer ncmax
      parameter (ncmax=65)
      real*8 x,y,z,xd,yd,zd,a
      integer ncdim
      dimension a(ncmax)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),
     *         xd(ncdim,*),yd(ncdim,*),zd(ncdim,*)
      integer m1,m2,m3,m4,i1,i2,k1,k2

      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
      logical odbg
      common/dbgdbg/odbg

      integer i, j, k, l
      integer n1, n2, n3, n4, nr
      real*8 fac

      do 20 i = 1 , ncontr
         a(i) = a(i) + a(i)
 20   continue
      n1 = 1
c
      do 120 i = 1 , m1
         n2 = n1
         do 110 j = 1 , m2
            n3 = n2
            do 100 k = 1 , m3
               n4 = n3
               do 90 l = 1 , m4
                  go to (30,50,70,70,70,70,70) , l
c
 30               do 40 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2)

                     
       if(odbg)write(6,*)'subsd1',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)


 40               continue
                  n4 = n4 + k2
                  go to 90
 50               do 60 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) - x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) - y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) - z(nr,n4-k2)

       if(odbg)write(6,*)'subsd2',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)

 60               continue
                  n4 = n4 + k2
                  go to 90
 70               fac = -dble(l-1)
                  do 80 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) + fac*x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) + fac*y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) + fac*z(nr,n4-k2)

       if(odbg)write(6,*)'subsd3',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)

 80               continue
c
                  n4 = n4 + k2
 90            continue
               n3 = n3 + k1
 100        continue
            n2 = n2 + i2
 110     continue
         n1 = n1 + i1
 120  continue
      return
      end
      subroutine redund_dft(ii,jj,kk,ll,
     +     basi, basj, bask, basl, 
     +     iw)

c      implicit real*8  (a-h,o-z)
      implicit none

      integer ii,jj,kk,ll,iw
      integer basi, basj, bask, basl

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
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
      integer iat, jat, kat, lat
      integer i,lll,lit,ljt,llt,lkt,iper
      integer min, imin
      integer n1, n2
      logical inej,inek,inel,jnek,jnel,knel

      dimension lll(4)
      equivalence (lll(1),lit)
      oskip(1) = .true.
      oskip(2) = .true.
      oskip(3) = .true.
      oskip(4) = .true.
      npass = 0
      do 20 i = 1 , 4
         natomd(i) = 0
 20   continue

      lit = ktype(basi,ii)
      iat = katom(basi,ii)
      if(basj.lt.0)then
         ljt = 1
         jat = iat
      else
         ljt = ktype(basj,jj)
         jat = katom(basj,jj)
      endif

      lkt = ktype(bask,kk)
      kat = katom(bask,kk)
      if(basl.lt.0)then
         llt = 1
         lat = kat
      else
         llt = ktype(basl,ll)
         lat = katom(basl,ll)
      endif

      inej = iat.ne.jat
      inek = iat.ne.kat
      inel = iat.ne.lat
      jnek = jat.ne.kat
      jnel = jat.ne.lat
      knel = kat.ne.lat

      if (inej) then
         if (.not.(inek)) then
            if (.not.(inel)) go to 40
c      iat=kat    jat=lat
            if (jnel) go to 50
            if (ii.ne.kk .or. jj.ne.ll) then
               n1 = (lit+1)*(lkt+1)*ljt*llt
               n2 = lit*lkt*(ljt+1)*(llt+1)
               if (n1.ge.n2) go to 50
               go to 70
            else
               if (ljt.le.lit) go to 60
               go to 40
            end if
         else if (jnek) then
            if (.not.(jnel)) go to 70
            if (.not.(knel)) go to 80
c     iat # jat # kat # lat  -- omit one centre
            min = lit
            imin = 1
            do 30 iper = 2 , 4
               if (lll(iper).lt.min) then
                  min = lll(iper)
                  imin = iper
               end if
 30         continue
            go to (90,100,110,120) , imin
            go to 90
         else
            if (.not.(jnel)) go to 60
c     ----- jat = kat # iat # lat -----
            oskip(1) = .false.
            oskip(4) = .false.
            natomd(1) = iat
            natomd(2) = lat
            natomd(3) = jat
            npass = 2
            go to 130
         end if
      else if (inek) then
         if (.not.(knel)) then
c     iat=jat  ,  kat=lat   differentiate one pair
            n1 = lit*ljt*(lkt+1)*(llt+1)
            n2 = (lit+1)*(ljt+1)*lkt*llt
            if (n2.lt.n1) go to 80
         end if
      else
         if (inel) then
c     iat=jat=kat   derivative (ij/kl')
            oskip(4) = .false.
            natomd(1) = lat
            natomd(2) = iat
            npass = 1
         end if
         go to 130
      end if
c     iat=jat   derivatives (ij/k'l) and (ij/kl')
      oskip(3) = .false.
      oskip(4) = .false.
      natomd(1) = kat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     iat=kat=lat   derivative (ij'/kl)
 40   oskip(2) = .false.
      natomd(1) = jat
      natomd(2) = iat
      npass = 1
      go to 130
c     iat=kat   derivatives (ij'/kl) and (ij/kl')
 50   oskip(2) = .false.
      oskip(4) = .false.
      natomd(1) = jat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     jat=kat=lat    (i'j/kl)
 60   oskip(1) = .false.
      natomd(1) = iat
      natomd(2) = jat
      npass = 1
      go to 130
c      jat=lat    derivatives (i'j/kl) and (ij/k'l)
 70   oskip(1) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = kat
      natomd(3) = jat
      npass = 2
      go to 130
c     kat=lat   derivatives (i'j/kl) and (ij'/kl)
 80   oskip(1) = .false.
      oskip(2) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      npass = 2
      go to 130
 90   natomd(1) = jat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = iat
      npass = 3
      oskip(2) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 100  natomd(1) = iat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = jat
      npass = 3
      oskip(1) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 110  natomd(1) = iat
      natomd(2) = jat
      natomd(3) = lat
      natomd(4) = kat
      npass = 3
      oskip(1) = .false.
      oskip(2) = .false.
      oskip(4) = .false.
      go to 130
 120  oskip(1) = .false.
      oskip(2) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      natomd(4) = lat
      npass = 3
c     -----
 130  if (.not.outd) return
      write (iw,6010) ii , jj , kk , ll , 
     +    oskip(1) , oskip(2) , oskip(3) , oskip(4) , 
     +    npass , (natomd(i),i=1,4)
      return
 6010 format (/,' ********  ii,jj,kk,ll =',4i3,' skip1,2,3,4 =',4l3,
     +        ' npass =',i2,' centres =',4i5,/)
      end

      subroutine grdcon_dft(qq,noform,gmax)
c      implicit real*8  (a-h,o-z)
      implicit none
c
c arguments
c
      real*8 qq
      logical noform
      dimension noform(*),qq(*)

      real*8 gmax


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

ccINCLUDE(../m4/common/specal)
  
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

c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dft_dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      integer iwr
      common/dft_iofile/iwr
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen

c-      real*8 dipd,dipn,dipi
c-      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
c
      logical ciopt, ciforc, mp2, hfgr, bfgs, ump2, lmeth2
      logical ump3, rmp3, ordmo, mp2w, loptor, ladp, lcpf
      logical lopti, lmcscf, lforce, lci, lcart, lmcdat
      logical lfdtrn, unit7, lcontr, lvcd, lgten
      logical ldenom, ignore, ldens, lset, ladapt, lsym, latmol
      logical berny, llibry, limpt, fpres, oss, ldiag, lskip
      logical opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common /restrl/ ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,
     +rmp3,ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c

      logical odbg
      common/dbgdbg/odbg

c
c functions
c
      real*8 ddot
c
c local variables
c
      integer ii, jj, kk, ll
      integer i1, j1, k1, l1
      integer n0, n1, n2, n3, nn, nnn, np3, npp
      integer lmax, jmax
      real*8 dumx, dumy, dumz, densty

      integer iii

c-      integer i, j, k, l
c-      logical nofk
c-      integer ioa, iob
c-      real*8 t1

      real*8 half
      data half/0.5d0/

c-      nofk = .not.(ofock .or. ofokab)
c
c
c-      if ((mp2w .or. ompir)) then
c-         ioa = ic7 + lendd*npass*3
c-         call vclr(qq(ioa+1),1,3*lendd)
c-         do 30 i = 1 , npass
c-            iob = ic7 + (i-1)*lendd*3
c-c-c-c-            do 20 k = 1 , lendd*3
c-               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
c- 20         continue
c- 30      continue
c-      end if
c-      if (ompir) then
c-         n2 = iabd + lendd + 1
c-         do 60 i = 1 , 3
c-            n1 = ic7 + 1
c-            do 50 k = 1 , npass + 1
c-               do 40 j = 1 , 3
c-                  t1 = ddot(ijkld,qq(n2),1,qq(n1),1)
c-                  dipi(i,j,natomd(k)) = dipi(i,j,natomd(k)) - t1
c-                  n1 = n1 + lendd
c- 40            continue
c- 50         continue
c-            n2 = n2 + lendd
c- 60      continue
c-      end if
c
c-      if (outd .or. ofock .or. ofokab) then
        if (outd) then
c
c     used only if extra output required or if need
c     derivatives of fock matrices
c
         jmax = maxj
         nn = 0
         nnn = 0
         do 110 ii = mini , maxi
            i1 = loci + ii
            if (oiandj) jmax = ii
            do 100 jj = minj , jmax
               j1 = locj + jj
c
c
               lmax = maxl
               do 90 kk = mink , maxk
                  k1 = lock + kk
                  if (okandl) lmax = kk
                  do 80 ll = minl , lmax
                     nn = nn + 1
                     if (.not.(noform(nn))) then
                        nnn = nnn + 1
                        l1 = locl + ll
                        densty = qq(iabd+nnn)
                        n0 = 1
                        n1 = ic7
                        do 70 npp = 1 , npass
                           n2 = n1 + lendd
                           n3 = n2 + lendd
                           dumx = qq(n1+nnn)
                           dumy = qq(n2+nnn)
                           dumz = qq(n3+nnn)
                           dgout(n0) = dgout(n0) + densty*dumx
                           dgout(n0+1) = dgout(n0+1) + densty*dumy
                           dgout(n0+2) = dgout(n0+2) + densty*dumz
                           if (outd) write(iwr,6010) npp , i1 , j1 , 
     +                         k1 , l1 , n1 , n2 , n3 , ic7 , lendd , 
     +                         nn , nnn , dumx , dumy , dumz , densty
c-                           if (.not.(nofk)) then
c-                              if (i1.eq.j1) then
c-                                 dumx = dumx*half
c-                                 dumy = dumy*half
c-                                 dumz = dumz*half
c-                              end if
c-                              if (k1.eq.l1) then
c-                                 dumx = dumx*half
c-                                 dumy = dumy*half
c-                                 dumz = dumz*half
c-                              end if
c-                              qq(n1+nnn) = q4*dumx
c-                              qq(n2+nnn) = q4*dumy
c-                              qq(n3+nnn) = q4*dumz
c-                           end if
                           n0 = n0 + 3
                           n1 = n3 + lendd
 70                     continue
                     end if
 80               continue
 90            continue
 100        continue
 110     continue
      else
c
c gmax used for debug purposes
c
c         n1 = ic7 + 1
c         gmax = 0.0
c         do n0 = 1 , np3
c            gmax = max(gmax, abs(ddot(ijkld,qq(iabd+1),1,qq(n1),1)))
c            n1 = n1 + lendd
c         enddo
c         if(abs(gmax) .gt. 1.0d0)odbg =  .true.

         n1 = ic7 + 1
         np3 = npass*3

         if(odbg)write(6,*)'dens',(qq(iabd+iii),iii=1,ijkld)
         do 120 n0 = 1 , np3

            if(odbg)write(6,*)'integ',n0,(qq(n1+iii-1),iii=1,ijkld)

            dgout(n0) = dgout(n0) + ddot(ijkld,qq(iabd+1),1,qq(n1),1)
            if(odbg)write(6,*)'dgout',n0,
     &           ddot(ijkld,qq(iabd+1),1,qq(n1),1)

            n1 = n1 + lendd
 120     continue

      end if
 6010 format (1x,12i6/30x,3e16.8,5x,e16.8)
      end
      subroutine formeg_dft

      implicit none

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
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      real*8 de
      common /dft_grad2/ de(3,maxat)
c


      logical odbg
      common/dbgdbg/odbg

      real*8 dgout
      common/tgrad/dgout(3,3)

      real*8 dum, dumx, dumy, dumz
      integer ipass, iat

      real*8 zero
      data zero /0.0d0/

      dumx = zero
      dumy = zero
      dumz = zero

      if(odbg)write(6,*)'formeg',dgout

      do 20 ipass = 1 , npass
         iat = natomd(ipass)
         dum = dgout(1,ipass)
         dumx = dumx + dum
         de(1,iat) = de(1,iat) + dum
         dum = dgout(2,ipass)
         dumy = dumy + dum
         de(2,iat) = de(2,iat) + dum
         dum = dgout(3,ipass)
         dumz = dumz + dum
         de(3,iat) = de(3,iat) + dum
 20   continue
      iat = natomd(npass+1)
      de(1,iat) = de(1,iat) - dumx
      de(2,iat) = de(2,iat) - dumy
      de(3,iat) = de(3,iat) - dumz
      return
      end
c
c normalisation factors for contraction with fitting coefficients
c
      block data dft_dfgdat

      real*8 fac
      common/dft_dfgfac/fac(35)

      data fac/7*1.0d0,
     &     3*0.5773502691896258d0,
     &     3*1.0d0,
     &     6*0.4472135954999579d0,
     &     0.2581988897471611d0,
     &     3*1.0d0,
     &     6*0.3779644730092272d0,
     &     3*0.2927700218845599d0,
     &     3*0.1690308509457033d0/
      end
c
c
      subroutine dabab_dft(ii,jj,kk,ll,
     &     basi, basj, bask, basl,
     &     ncentr,q4,
     &     zscftp,da,db,cfit,cfit2,abdens)

      implicit none
c
c arguments
c
      real*8 da(*), db(*), cfit(*), cfit2(*), abdens(*)
      integer ii,jj,kk,ll,basi, basj, bask, basl	
      real*8 q4
      character*8 zscftp
      integer ncentr

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
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
c local variables
c
      integer mjk, mkl, mil, mjl, mij, mik
      integer i, loci, locj, lock, locl, nn
      integer l, l1, ll1, ii1, j, jj1, j1, k1, k, kk1
      integer ni, nj, nk, nl
      integer mink, maxk, minl, maxl, mini, maxi, minj, maxj
      integer i1, i2, i3, i4
      integer n, is, js, isi, joff
      logical okleq, oijeq
      real*8 dfac
      integer itmp

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
      real*8 fac
      common/dft_dfgfac/fac(35)

      real*8 wght

      real*8 pt5,four

      logical odbg
      common/dbgdbg/odbg

      data pt5,four /0.5d0,4.0d0/

c
c      ouhf = zscftp.eq.zuhf
c      ogrhf = zscftp.eq.zgrhf

c
c make conditional on ncentr
c

      ni = 1
c
c-      if (.not.ogrhf) then
c

c         call chkadr2(da(1),itmp)
c         write(6,*)'addr dens input',itmp, da(1),da(2),da(3),da(4)

      if(ncentr .eq. 4) then

         mini = kmin(basi,ii)
         minj = kmin(basj,jj)
         mink = kmin(bask,kk)
         minl = kmin(basl,ll)

         maxi = kmax(basi,ii)
         maxj = kmax(basj,jj)
         maxk = kmax(bask,kk)
         maxl = kmax(basl,ll)

         loci = kloc(basi,ii) - mini
         locj = kloc(basj,jj) - minj
         lock = kloc(bask,kk) - mink
         locl = kloc(basl,ll) - minl

         do 60 i = mini , maxi
            nj = ni
            do 50 j = minj , maxj
               nk = nj
               do 40 k = mink , maxk
                  nl = nk
                  do 30 l = minl , maxl
                     nn = nl
                     i1 = loci + i
                     i2 = locj + j
                     i3 = lock + k
                     i4 = locl + l
                     if (i1.lt.i2) then
                        n = i1
                        i1 = i2
                        i2 = n
                     end if
                     if (i3.lt.i4) then
                        n = i3
                        i3 = i4
                        i4 = n
                     end if
                     if (i1.lt.i3) then
                     else if (i1.eq.i3) then
                        if (i2.ge.i4) go to 20
                     else
                        go to 20
                     end if
                     n = i1
                     i1 = i3
                     i3 = n
                     n = i2
                     i2 = i4
                     i4 = n
 20                  mij = iky(i1) + i2
                     mik = iky(i1) + i3
                     mil = iky(i1) + i4
                     mkl = iky(i3) + i4
                     if (i2.lt.i3) then
                        mjk = iky(i3) + i2
                        if (i2.lt.i4) then
                           mjl = iky(i4) + i2
                        else
                           mjl = iky(i2) + i4
                        end if
                     else
                        mjk = iky(i2) + i3
                        mjl = iky(i2) + i4
                     end if

c-                     if(.not. CD_active())then
c-                        dfac = da(mij)*da(mkl)*four - da(mik)*da(mjl)
c-     +                       - da(mil)*da(mjk)
c-                        if (ouhf) dfac = dfac - db(mik)*db(mjl) 
c-     +                       - db(mil)*db(mjk)
c-                     else
c
c optionally include coulomb and fraction of exclude exchange contribution
c
c-                        if (ouhf) call caserr('no uhf dft yet')

                        dfac=0.0d0

                        if(CD_HF_coulomb_deriv())then
                           dfac = da(mij)*da(mkl)*four 
                           if (.not.rks_sw) then
                              dfac = dfac + db(mij)*db(mkl)*four
                           endif
                        endif

                        if(CD_HF_exchange())then
                           wght = CD_HF_exchange_weight()
                           dfac = dfac - 
     &                         wght*(da(mik)*da(mjl) + da(mil)*da(mjk))
                           if (.not.rks_sw) then
                              dfac = dfac - 
     &                         wght*(db(mik)*db(mjl) + db(mil)*db(mjk))
                           endif
                        endif
c-                     endif

                     if (i1.eq.i2) dfac = dfac*pt5
                     if (i3.eq.i4) dfac = dfac*pt5
                     dfac = dfac*q4

c-                     if (omp2 .or. mp3) then
c-                        abdens(nn) = abdens(nn) + dfac
c-                     else
                        abdens(nn) = dfac
c-                     end if

                     nl = nl + inc2
 30               continue
                  nk = nk + inc3
 40            continue
               nj = nj + inc4
 50         continue
            ni = ni + inc5
 60      continue

      else if (ncentr .eq. 3) then

         mini = kmin(basi,ii)
         minj = kmin(basj,jj)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxj = kmax(basj,jj)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         locj = kloc(basj,jj) - minj
         lock = kloc(bask,kk) - mink

         ni = 1
         do  i = mini , maxi
            nj = ni
            do j = minj , maxj
               nk = nj

               i1 = loci + i
               i2 = locj + j

               if (i1.lt.i2) then
                  n = i1
                  i1 = i2
                  i2 = n
               end if

               mij = iky(i1) + i2

               do k = mink , maxk
                  nn = nk

                  dfac = 2.0d0*da(mij)*cfit(lock + k)*fac(k)
                  if (.not.rks_sw) then
                     dfac = dfac + 2.0d0*db(mij)*cfit(lock + k)*fac(k)
                  endif
                  if (i1.eq.i2) dfac = dfac*pt5
                  dfac = dfac*q4
                  abdens(nn) = dfac

                  nk = nk + inc3
               enddo
               nj = nj + inc4
            enddo
            ni = ni + inc5
         enddo

      else if (ncentr .eq. 2) then
c
         if(basi .eq. bask)then

         mini = kmin(basi,ii)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         lock = kloc(bask,kk) - mink

         ni = 1
         do i = mini , maxi
            nk = ni
            do k = mink , maxk
               nn = nk

               if(odbg)write(6,*)'dabab',nn,loci,lock,i,k,
     &	 cfit(loci + i),cfit(lock + k),fac(i),fac(k)

               abdens(nn) = - q4 * cfit(loci + i) * cfit(lock + k) *
     +              fac(i) * fac(k)
c NB diagonal terms give no gradients
c               if(i .eq. k)abdens(nn) = abdens(nn) * 0.5d0
               nk = nk + inc3
            enddo
            ni = ni + inc5
         enddo

         else
c
c Two distinct interacting charge distributions
c
         mini = kmin(basi,ii)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         lock = kloc(bask,kk) - mink

         ni = 1
         do i = mini , maxi
            nk = ni
            do k = mink , maxk
               nn = nk

c               if(odbg)write(6,*)'dabab',nn,loci,lock,i,k,
c     &	 cfit(loci + i),cfit(lock + k),fac(i),fac(k)

               abdens(nn) = - q4 * cfit(loci + i) * cfit2(lock + k) *
     +              fac(i) * fac(k)

               nk = nk + inc3
            enddo
            ni = ni + inc5
         enddo

         endif

      endif

c
c-      else
c-c
c-c     general case
c-c
c-         do 120 i = mini , maxi
c-            nj = ni
c-            i1 = loci + i
c-            ii1 = iky(i1)
c-            do 110 j = minj , maxj
c-               nk = nj
c-               j1 = locj + j
c-               jj1 = iky(j1)
c-               mij = ii1 + j1
c-               if (j1.gt.i1) mij = jj1 + i1
c-               oijeq = i1.eq.j1
c-               do 100 k = mink , maxk
c-                  nl = nk
c-                  k1 = lock + k
c-                  kk1 = iky(k1)
c-                  mik = ii1 + k1
c-                  if (k1.gt.i1) mik = kk1 + i1
c-                  mjk = jj1 + k1
c-                  if (k1.gt.j1) mjk = kk1 + j1
c-                  do 90 l = minl , maxl
c-                     nn = nl
c-                     l1 = locl + l
c-                     ll1 = iky(l1)
c-                     mkl = kk1 + l1
c-                     if (l1.gt.k1) mkl = ll1 + k1
c-                     mjl = jj1 + l1
c-                     if (l1.gt.j1) mjl = ll1 + j1
c-                     mil = ii1 + l1
c-                     if (l1.gt.i1) mil = ll1 + i1
c-                     okleq = k1.eq.l1
c-                     dfac = 0.0d0
c-                     ioff = 0
c-                     do 80 is = 1 , njk
c-                        isi = (is-1)*11
c-                        joff = 0
c-                        do 70 js = 1 , njk
c-                           dfac = dfac + 4.0d0*erga(isi+js)
c-     +                            *(da(ioff+mij)*da(joff+mkl)
c-     +                            +da(ioff+mkl)*da(joff+mij))
c-     +                            + 2.0d0*ergb(isi+js)
c-     +                            *(da(ioff+mik)*da(joff+mjl)
c-     +                            +da(ioff+mjl)*da(joff+mik)
c-     +                            +da(ioff+mil)*da(joff+mjk)
c-     +                            +da(ioff+mjk)*da(joff+mil))
c-                           joff = joff + nx
c- 70                     continue
c-                        ioff = ioff + nx
c- 80                  continue
c-                     if (oijeq) dfac = dfac*pt5
c-                     if (okleq) dfac = dfac*pt5
c-                     abdens(nn) = dfac*q4
c-                     nl = nl + inc2
c- 90               continue
c-                  nk = nk + inc3
c- 100           continue
c-               nj = nj + inc4
c- 110        continue
c-            ni = ni + inc5
c- 120     continue
c-c
c-      end if
c
      return
      end
      subroutine ver_dft_deriv2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/deriv2e.m,v $
     +     "/
      data revision /
     +     "$Revision: 6317 $"
     +      /
      data date /
     +     "$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
