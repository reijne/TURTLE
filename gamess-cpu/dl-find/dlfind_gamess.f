c
c        Interface between GAMESS-UK and DL-FIND
c
c
c
c
c   Copyright 2007 Johannes Kaestner (kaestner@theochem.uni-stuttgart.de),
c   Tom Keal (thomas.keal@stfc.ac.uk)
c 
c   This file is part of DL-FIND.
c 
c   DL-FIND is free software: you can redistribute it and/or modify
c   it under the terms of the GNU Lesser General Public License as 
c   published by the Free Software Foundation, either version 3 of the 
c   License, or (at your option) any later version.
c 
c   DL-FIND is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU Lesser General Public License for more details.
c 
c   You should have received a copy of the GNU Lesser General Public 
c   License along with DL-FIND.  If not, see 
c   <http://www.gnu.org/licenses/>.
c
c     ..................................................................
      subroutine dlfind_gamess(core)
      implicit none
      real*8 :: core(*)
c     local vars
      integer  :: nvar2, nprint_save,opg_root_int
      logical,external::  opg_root
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
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
c
      real*8 func0
      integer ncoord, npts, nserch, iupdat, icode, iupcod
      common /seerch/ func0,ncoord,npts,nserch,iupdat,icode,iupcod
c
c     runlab for zruntp so we call the wavefunction analysis
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
c     restar for nprint (needed by optend)
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
c     prnprn for print level (oprn)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c     iofile for iwr
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c     integer global_nodeid
c     external global_nodeid

      ! ****************************************************************
c      write(6,*)'calling dl-find icoord=',icoord
c      write(6,*)'calling dl-find ncons=',ncons_co
c      write(6,*)'calling dl-find nat=',nat
c      write(6,*)'calling dl-find nat=',nat

c     nprint is fiddled about with in the code to change the output
c     of the scf, so that there is less output after the first geometry
c     optimisation cycle. We therefore need to save the original nprint
      nprint_save=nprint
      npts = -1
      ncoord=3*nat
      nvar2=nat
      if(opg_root()) then
        opg_root_int=1
      else
        opg_root_int=0
      end if
      if(geom2.ne."") nvar2= nvar2 + 3*nat
c      write(6,*)'calling dl-find nvar=',3*nat
c      write(6,*)'calling dl-find nvar2=',nvar2
c      write(6,*)'calling dl-find nspec=',2*nat+5*ncons_co
      call dl_find(3*nat,nvar2,3*nat+5*ncons_co,opg_root_int,core)

c     write(6,*)'DL-FIND COMPLETED',global_nodeid()

c     restore print level
      nprint=nprint_save

c     header
      write(iwr,'(//1x,104("=")/)')

c     analysis of moment of inertia etc
      if(oprn(28)) call anamom(core)
c
c     analysis of bond lengths etc
      call intr(core)
c
c     Write out the orbitals
      call optend(core,nprint)
c
c     Set runtype to calculate wavefunction properties
      zruntp='prop'

      write(6,*)'DL-FIND RETURNING'

      return
      end
c
c     ..................................................................
      subroutine dlf_get_params(nvar,nvar2,nspec,
     +     coords,coords2,spec,ierr,
     +     tolerance,printl,maxcycle,maxene,
     +     tatoms,icoord_,iopt,iline,maxstep,
     +     scalestep,lbfgs_mem,nimage,nebk,
     +     dump,restart,nz,ncons,nconn,
     +     update,maxupd,delta,soft,inithessian,
     +     carthessian,tsrel,maxrot,tolrot,nframe,nmass,nweight,
     +     timestep,fric0,fricfac,fricp,
     +     imultistate,state_i,state_j,
     +     pf_c1,pf_c2,gp_c3,gp_c4,ln_t1,ln_t2,
     +     printf,tolerance_e,distort,massweight,minstep,maxdump,
     +     task, temperature, po_pop_size, po_radius,po_contraction,
     +     po_tolerance_r, po_tolerance_g, po_distribution,
     +     po_maxcycle, po_init_pop_size, po_reset, po_mutation_rate,
     +     po_death_rate, po_scalefac, po_nsave, ntasks, tdlf_farm, 
     +     n_po_scaling,
     +     neb_climb_test, neb_freeze_test,nzero, coupled_states, 
     +     qtsflag, imicroiter, maxmicrocycle, micro_esp_fit)
      implicit none
      integer   ,intent(in)      :: nvar 
      integer   ,intent(in)      :: nvar2
      integer   ,intent(in)      :: nspec
      real*8      ,intent(inout)   :: coords(nvar) ! start coordinates
      real*8      ,intent(inout)   :: coords2(nvar2) ! a real array that can be used
                                ! depending on the calculation
                                ! e.g. a second set of coordinates
      integer   ,intent(inout)   :: spec(nspec) ! specifications like fragment or frozen
      integer   ,intent(out)     :: ierr
      real*8      ,intent(inout)   :: tolerance
      integer   ,intent(inout)   :: printl
      integer   ,intent(inout)   :: maxcycle
      integer   ,intent(inout)   :: maxene
      integer   ,intent(inout)   :: tatoms
      integer   ,intent(inout)   :: icoord_
      integer   ,intent(inout)   :: iopt
      integer   ,intent(inout)   :: iline
      real*8      ,intent(inout)   :: maxstep
      real*8      ,intent(inout)   :: scalestep
      integer   ,intent(inout)   :: lbfgs_mem
      integer   ,intent(inout)   :: nimage
      real*8      ,intent(inout)   :: nebk
      integer   ,intent(inout)   :: dump
      integer   ,intent(inout)   :: restart
      integer   ,intent(inout)   :: nz
      integer   ,intent(inout)   :: ncons
      integer   ,intent(inout)   :: nconn
      integer   ,intent(inout)   :: update
      integer   ,intent(inout)   :: maxupd
      real*8      ,intent(inout)   :: delta
      real*8      ,intent(inout)   :: soft
      integer   ,intent(inout)   :: inithessian
      integer   ,intent(inout)   :: carthessian
      integer   ,intent(inout)   :: tsrel
      integer   ,intent(inout)   :: maxrot
      real*8      ,intent(inout)   :: tolrot
      integer   ,intent(inout)   :: nframe
      integer   ,intent(inout)   :: nmass
      integer   ,intent(inout)   :: nweight
      real*8      ,intent(inout)   :: timestep
      real*8      ,intent(inout)   :: fric0
      real*8      ,intent(inout)   :: fricfac
      real*8      ,intent(inout)   :: fricp
      integer   ,intent(inout)   :: imultistate
      integer   ,intent(inout)   :: state_i
      integer   ,intent(inout)   :: state_j
      real*8      ,intent(inout)   :: pf_c1
      real*8      ,intent(inout)   :: pf_c2
      real*8      ,intent(inout)   :: gp_c3
      real*8      ,intent(inout)   :: gp_c4
      real*8      ,intent(inout)   :: ln_t1
      real*8      ,intent(inout)   :: ln_t2
      integer   ,intent(inout)   :: printf
      real*8      ,intent(inout)   :: tolerance_e
      real*8      ,intent(inout)   :: distort
      integer   ,intent(inout)   :: massweight
      real*8      ,intent(inout)   :: minstep
      integer   ,intent(inout)   :: maxdump
c JMC new arguments
      integer   ,intent(inout)   :: task
      real*8      ,intent(inout)   :: temperature
      integer   ,intent(inout)   :: po_pop_size
      real*8      ,intent(inout)   :: po_radius
      real*8      ,intent(inout)   :: po_contraction
      real*8      ,intent(inout)   :: po_tolerance_r
      real*8      ,intent(inout)   :: po_tolerance_g
      integer   ,intent(inout)   :: po_distribution
      integer   ,intent(inout)   :: po_maxcycle
      integer   ,intent(inout)   :: po_init_pop_size
      integer   ,intent(inout)   :: po_reset
      real*8      ,intent(inout)   :: po_mutation_rate
      real*8      ,intent(inout)   :: po_death_rate
      real*8      ,intent(inout)   :: po_scalefac
      integer   ,intent(inout)   :: po_nsave
      integer   ,intent(inout)   :: ntasks
      integer   ,intent(inout)   :: tdlf_farm
      integer   ,intent(inout)   :: n_po_scaling
      real*8      ,intent(inout)   :: neb_climb_test
      real*8      ,intent(inout)   :: neb_freeze_test
      integer   ,intent(inout)   :: nzero
      integer   ,intent(inout)   :: coupled_states
      integer   ,intent(inout)   :: qtsflag
      integer   ,intent(inout)   :: imicroiter
      integer   ,intent(inout)   :: maxmicrocycle
      integer   ,intent(inout)   :: micro_esp_fit
c     local vars
      integer                    :: iat
      integer, external          :: jsubst
      real*8 , external            :: amass_get
      character(8)               :: ztag_(nvar/3)
      logical,external::  opg_root
      integer ifz
c     GAMESS common blocks
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
c infoa contains c(3,maxat)
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
c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
      ! ****************************************************************
c      write(6,*)'in get_params: nvar,nvar2,nspec',nvar,nvar2,nspec
      ierr=0
      if(nvar.ne.3*nat) call caserr2('nvar is not 3*nat in '//
     &     'dlf_get_params')
      coords(:)=reshape(c(1:3,1:nat),(/nvar/))
      nmass=nat
c     spec should contain: 
c       nvar spec (1 or constraint)
c       nvar atomic number (Z)
c       constraint data
      do iat=1,nat
        spec(nat+iat)=jsubst(ztag(iat))
      end do
c      print*,"spec: ",spec
      delta=delta_co
      icoord_=icoord
      if(icoord>=100.and.icoord<200) then
c       NEB
        nimage=nimage_co
        nebk=nebk_co
c       New options, not currently implemented in GAMESS
        neb_climb_test=-1.0d0
        neb_freeze_test=-1.0d0
      end if
c     Include all atoms into the first residue in case of internals
c     At the moment, input of HDLC residues is not possible
      if(mod(icoord,10)>0.and.mod(icoord,10)<5) then
        spec(1:nat)=1
      end if
      if (opg_root())write(6,*)' Number of frozen atoms::',nfreeze
      if (opg_root())write(6,*)' Frozen atoms are:'
      do ifz=1,nfreeze
        spec(ifreeze(ifz))=-1
        if (opg_root())write(6,*)ifreeze(ifz)
      enddo
c     Optimiser
      iopt=iopt_co
      lbfgs_mem=mem_co
      timestep=time_co
      fric0=fri0_co
      fricfac=frif_co
      fricp=frip_co
      if((iopt.eq.3.or.iopt.lt.0).and.icoord.lt.100) then
        iline=1
      end if
c     Update
      update=upd_co
      maxupd=rec_co
CTWK  This is correct if fd_co=1 -> one point, =2 -> two point FD:
      inithessian = fd_co
      soft=soft_co
      tolerance=4.D0/9.D0*tol_co
      tolerance_e=1.D0
      if(rst_co) restart=1
      dump=dump_co
c
      maxcycle=maxc_co
      maxstep=maxs_co
      printl=4
      nz=nat
      ncons=ncons_co
      printf=6
c     coords2
      if(geom2.ne."") then
        call read_xyz(geom2,nat,ztag_,coords2(1:nvar),iat)
        if(iat.ne.nat) call caserr2('number of atoms in second set of'//
     &       ' coordinates not equal to first')
        nframe=1
c       the masses
        do iat=1,nat
          coords2(iat+3*nat)=amass_get(1,iat)
        end do
      else
        nframe=0
c       the masses
        do iat=1,nat
          coords2(iat)=amass_get(1,iat)
        end do
      endif
      nweight=0
c      print*,"coords2",coords2
c     miscellaneous
      task=task_co
      temperature=temperature_co
c     parallel optimization
      po_pop_size=po_pop_size_co
      po_radius=po_radius_co
      po_contraction=po_contraction_co
      po_tolerance_r=po_tolerance_r_co
      po_tolerance_g=po_tolerance_g_co
      po_distribution=po_distribution_co
      po_maxcycle=po_maxcycle_co
      po_init_pop_size=po_init_pop_size_co
      po_reset=po_reset_co
      po_mutation_rate=po_mutation_rate_co
      po_death_rate=po_death_rate_co
      po_scalefac=po_scalefac_co
      po_nsave=po_nsave_co

      n_po_scaling=0 ! ??? for testing

c     taskfarming
      ntasks=ntasks_co

c JMC hardwire the tdlf_farm=0 for now -- is there any need for 
c tdlf_farm to be settable from within gamess?
c this will cause dl-find to load the split commss stuff from gamess

      tdlf_farm=0

      end 
c
c     ..................................................................
      subroutine dlf_get_gradient(nvar,coords,energy,gradient,iimage,
     +     kiter,core,status)
      implicit none
      integer   ,intent(in)    :: nvar
      integer   ,intent(in)    :: kiter
      real*8      ,intent(in)    :: coords(nvar)
      real*8      ,intent(out)   :: energy
      real*8      ,intent(out)   :: gradient(nvar)
      integer   ,intent(in)    :: iimage
      real*8                     :: core(*)
      integer   ,intent(out)   :: status
c     local vars
      real*8     :: co,g,dx,func
      integer,external :: lensec
      integer          :: isize,m17=17,natdlf
c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
c     it seems we have to rely on this common block from valopt (m4/optim.m)
      common/miscop/co(maxat*3),g(maxat*3),dx(maxat*3),func
c
c ----- allocate gradient section on dumpfile
c
      isize = lensec(nvar*nvar) + lensec(mach(7))
      call secput(isect(495),m17,isize,ibl3g)
      ibl3hs = ibl3g + lensec(mach(7))
c
c ----- set coordinates and calculate energy and gradient
c
      natdlf=nvar/3
      c(1:3,1:natdlf)=reshape(coords(:),(/3,natdlf/))

c     calculate energy and gradient in gamess
      call valopt(core)

c     get back the resulting energy and gradient
      energy=func
      gradient(:)=egrad(1:nvar)

      status=irest
      end
c
c     ..................................................................
      subroutine dlf_get_hessian(nvar,coords,hessian,status)
c     subroutine dlf_get_hessian(nvar,coords,hessian,status,core)
c                               !  get the hessian at a given geometry
c     This routine does not work for the moment, as core is not passed to it
      implicit none
      integer   ,intent(in)    :: nvar
      real*8      ,intent(in)    :: coords(nvar)
      real*8      ,intent(out)   :: hessian(nvar,nvar)
      integer   ,intent(out)   :: status
c     real*8                     :: core(*)
      real*8                     :: core
c     local vars
      integer          :: isize,nc2,len3,lenc,l3,m17=17
      integer,external :: lensec
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
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
c     ******************************************************************

      status=1
      return

      call flushout()
c     
c ----- allocate gradient section on dumpfile
c
      l3 = num*num
      len3 = lensec(l3)
      lenc = lensec(mach(9))
      nc2 = nvar**2
      isize = lensec(nc2) + lensec(mach(7))
      call secput(isect(495),m17,isize,ibl3g)

      c(1:3,1:nvar)=reshape(coords(:),(/3,nvar/))
      call flushout()
      if (omp2) then
         call mp2dd(core)
      else if( mp3) then
         call caserr2
     +   ('correlated second derivatives - only mp2 allowed')
      else
        call flushout()
        call scfdd(core)
        call flushout()
        call flushout()
        call anairr(core)
        call flushout()
        call flushout()
      end if
c
c     read cartesian force constant matrix
c
c      call rdfcm(vec,'iranal')
      call rdfcm(hessian,'dl-find')

c      hessian(:,:)=0.D0
      status=0
      end

c
c     ..................................................................
      subroutine dlf_get_multistate_gradients(nvar,coords,energy,
     +     gradient,coupling,needcoupling,iimage,status)
c     Dummy routine for multistate gradients. This option is 
c     not supported in Gamess
      implicit none
      integer     :: nvar
      real*8        :: coords
      real*8        :: energy
      real*8        :: gradient
      real*8        :: coupling
      integer     :: needcoupling
      integer     :: iimage
      integer     :: status
      end
c
c     ..................................................................
      subroutine dlf_put_coords(nvar,mode,energy,coords,iam)
      implicit none
      integer   ,intent(in)    :: nvar
      integer   ,intent(in)    :: mode
      integer   ,intent(in)    :: iam
      real*8      ,intent(in)    :: energy
      real*8      ,intent(in)    :: coords(nvar)
      end

c
c     ..................................................................
      subroutine dlf_update()
      implicit none
c     dummy routine here
      end

      subroutine dlf_error()
      implicit none
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
      call gamerr(
     &     'DL-FIND error',
     &     ERR_NO_CODE, ERR_UNLUCKY, ERR_SYNC, ERR_NO_SYS)
      end
c
c     ..................................................................
      subroutine dlf_geom(maxat,ztag,cat,nat)
c     Read in the geometrie(s)
c     at the moment, the file name is hardcoded to geom.xyz !!
c     a second geometry may be read in from a file specified after the keyword geom
      implicit none
      integer      ,intent(in) :: maxat
      character(8) ,intent(out):: ztag(maxat)
      real*8         ,intent(out):: cat(3,maxat)
      integer      ,intent(out):: nat
c     ******************************************************************
      call read_xyz("geom.xyz",maxat,ztag,cat,nat)
      end
c

      logical function isalpha(c)
c     true for ascii [A-Z,a-z]
      implicit none
      character(1), intent(in)::c
      character(1)            ::a,z
      logical                 ::alpha
      a='a'
      z='z'
      alpha=((iachar(c).ge.iachar(a)).and.(iachar(c).le.iachar(z)))
      a='A'
      z='Z'
      alpha=alpha.or.
     &((iachar(c).ge.iachar(a)).and.(iachar(c).le.iachar(z)))
      isalpha=alpha
      end function isalpha

c     ..................................................................
*     read_xyz now accepts both
*     neb .xyz output files (cartesian angstrom)
*         header           : # of atoms
*         with atom records: tag x y z
*     and gamess 'nuclear coordinates' blocks (cartesian,au or angstrom)
*         header           : <au/an> # of atoms
*         with atom records: x y z qq(ivoff+) tag
*     fname no longer a dummy.

      subroutine read_xyz(fname,maxat,ztag,cat,nat)
      implicit none
      integer      ,intent(in) :: maxat
      character(8) ,intent(out):: ztag(maxat)
      real*8         ,intent(out):: cat(3,maxat)
      integer      ,intent(out):: nat
c     local vars
      logical       :: tchk,guessunit,isalpha
      logical       :: dlfformat,gamessformat
      integer       :: iunit=5111,iat,wantnat
      character(2)  :: ch2
      real*8          :: ang=0.529177249d0,scalecoor,fdum
      character(32) :: theunits
      character(256):: buffer
c     old string type
      character*(*) fname
c     ******************************************************************
      guessunit=.false.
      dlfformat=.false.
      gamessformat=.false.
      inquire(FILE=fname,EXIST=tchk)
      if(.not.tchk) then
c       print*,"Input geometry NOT read from xyz file"
        nat=0
        return
      end if
c     print*,'dl-find: reading second geometry from xyz file'
c     print*,'maxat = ',maxat
      open (unit=iunit,file=fname,err=201)
c     possible header: 'angstrom/au #atoms'
      read(iunit,*,err=100,end=200) theunits,wantnat
      goto 150
100   continue
c     print*,'error on read 1 theunits,nat'
      guessunit=.true.
      rewind(iunit)
c     possible header:  '#atoms'
      read(iunit,*,err=110,end=200) wantnat
      goto 150
110   continue
c     no header a all...
c     print*,'error on read 2 nat'
      rewind(iunit)

150   continue
      if (.not.guessunit) then
*        print*,'dl-find: read_xyz, file uses ',theunits
*        print*,'dl-find: read_xyz, file has ',wantnat,' atoms'
         if (theunits(1:2).eq.'an'.or.theunits(1:2).eq.'AN') then
            scalecoor=ang
         else
            scalecoor=1d0
         endif
      else
*        original default: angstrom
         scalecoor=ang
      endif

      nat=0
      iat=1
      do while (.true.)
151      continue
         if (gamessformat) then
            read(iunit,'(A256)',end=161,err=200) buffer
            do while (len(trim(buffer)).lt.6) 
*              print*,'dl-find: read_xyz skipping line ',buffer
               read(iunit,'(A256)',end=161,err=200) buffer
            end do
            read(buffer,*,err=200) cat(:,iat),fdum,ch2
*           print*,'ch2 is ',ch2,' cat= ',cat(:,iat)
         else if (dlfformat) then
            read(iunit,'(A256)',end=161,err=200) buffer
            do while (len(trim(buffer)).lt.6) 
*              print*,'dl-find: read_xyz skipping line ',buffer
               read(iunit,'(A256)',end=161,err=200) buffer
            end do
            read(buffer,*,err=200) ch2,cat(:,iat)
*           print*,'ch2 is ',ch2,' cat= ',cat(:,iat)
         else
            read(iunit,*,end=161,err=200) ch2
            backspace(iunit)
            if (isalpha(ch2(1:1))) then
               dlfformat=.true.
*           print*,'dl-find: read_xyz, file uses default format'
            else
               gamessformat=.true.
*           print*,'dl-find: read_xyz, file uses gamess format'
            endif
            goto 151
         end if
         ztag(iat)=ch2
         iat=iat+1
      end do
161   continue
*     normal exit
      nat=iat-1
*     print*,'found ',nat,' atoms'
      close(iunit)
      cat=cat/scalecoor
      return

*     (premature) EOF
 201  print*,"dl-find read_xyz Error: end of input file reached"
      go to 203
*     read error on input record
 200  print*,"dl-find read_xyz Error reading input record"
 203  close(iunit)
      go to 205
*     open failed
 204  print*,"dl-find read_xyz Error opening input file"
 205  nat=0
      return
      end
c
c     ..................................................................
      subroutine dlf_output(guk_stdout,guk_stderr)
      use dlf_global, only: glob,stderr,stdout,keep_alloutput
      implicit none
      integer :: guk_stdout
      integer :: guk_stderr
c     ******************************************************************

      if (guk_stdout >= 0) stdout = guk_stdout
      if (guk_stderr >= 0) stderr = guk_stderr

      if (glob%nprocs > 1) then
        ! write some info on the parallelization
        write(stdout,'(1x,a,i10,a)')"I have rank ",glob%iam,
     &    " in mpi_comm_world"
        write(stdout,'(1x,a,i10)')"Total number of processors = ",
     &    glob%nprocs
        if (keep_alloutput) then
           write(stdout,'(1x,a)')"Keeping output from all processors"
        else
           write(stdout,'(1x,a)')
     &    "Not keeping output from processors /= 0"
        end if
      end if

      return
      end subroutine dlf_output
