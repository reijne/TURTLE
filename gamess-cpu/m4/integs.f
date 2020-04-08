c 
c  $Author: hvd $
c  $Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
c  $Locker:  $
c  $Revision: 6317 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integs.m,v $
c  $State: Exp $
c  
      subroutine jkin7a(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

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
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
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
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
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
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5
c
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0
      if(.not.odscf) then
      else
        l2 = nx
      endif
c
      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     ----- ishell -----
c
      call filmax
      do 400 ii = ist0,nshels
c
c     ----- print intermediate restart data -----
c
      if(kad(ii))400,401,401
 401  dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
       endif
c
c     ----- eliminate ishell -----
c
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
c
c     ----- jshell -----
c
      j0 = jst0
      do 380 jj = j0,ii
      jst0 = 1
      if(kad(jj))380,141,141
c ***
c ***
141       itrij=ikyii+jj
          tolij=dlntol+prefac(itrij)
c
c **** rejecting on basis of ij arguably too extreme
c ***  suppress
c         if(tolij .le. -3.401e0) then
c           if(odscf)intcut(1)=intcut(1)+1
c           goto 380
c         endif
c ***
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
      ishell = ii
      jshell = jj
      mij = itrij
c
c     ----- kshell -----
c
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
 361  do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ikykk=iky(kk)
      if (odscf) then
         itrik=ikyii+kk
         mik = itrik
         itrjk=iky(max(jj,kk))+min(jj,kk)
         mjk = itrjk
         tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
c
c     ----- lshell ----
c
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
      if(kad(ll).lt.0) goto 340
      itrkl=ikykk+ll
      tijkl=tolij+prefac(itrkl)
      if(tijkl.le.0.0d0) then
        if(odscf)intcut(2)=intcut(2)+1
        goto 340
      endif
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=tijkl+
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            intcut(3)=intcut(3)+1
            goto 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dble(nt)/ dble(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
      if (odscf) then
        if (zscftp.eq.'uhf')then
          call genr70(gout,1,.false.)
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout,
     +         fac1,fac2, facex, ocoul, oexch)
        else if (zscftp.eq.'gvb') then
          call genr70(gout,1,.false.)
          if(nsheld.le.1) then
            call dir_build_open_70(fock,exch,dens,gout)
          else
            call dir_build_open2_70(l2,fock,exch,dens,gout)
          endif
        else
          if(omorok) then
            call genr70(gout,1,.false.)
            call dbuild70_morok(fock,dens,gout)
          else
            call genr70(gout,1,.false.)
            call dbuild70(fock,dens,gout,
     &           fac1, fac2, facex, ocoul, oexch)
          endif
        endif
      else
        call genr70(gout,1,.false.)
        call qout70(gout)
      endif
c
c
c     ----- check cpu time/ maxblock condition -----
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim)go to 420
c
  340 continue
  360 continue
  380 continue
      time = cpulft(1)
  400 continue
      ist = 1
      jst = 1
      kst = 1
      lst = 1
      call final(q,fock,dens)
c
c     ----- reset core memory
c
  420 continue
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
      end
      subroutine jkin7s(zscftp,q,fock,fockb,exch,dens,densb,
     +                  prefac,rdmat,iso,gout,nshels,outvv)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *59 text
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

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
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
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
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
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
      dimension fock(*),fockb(*),exch(*),dens(*),densb(*)
      dimension prefac(*),rdmat(*)
      dimension iso(nshels,*),q(*),gout(*)
      dimension mi(48),mj(48),mk(48),m0(48)
      character *8 fnm
      character *6 snm
      data fnm/'integs.m'/
      data snm/'jkin7s'/
      data m25/25/
c     data m22,mword1/22,6000/
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
c
c
c for simplicity, only implement full coulomb version
c here for the moment
c
      if (CD_2e())then
        odft = .true.

        if(CD_request_multstate())
     &    call caserr('use integ high noschwarz for multipole runs')

        if(CD_HF_exchange())then
           facex = CD_HF_exchange_weight()
           oexch = .true.
        else
           facex = 0.0d0
           oexch = .false.
        endif
      else
         odft = .false.
         facex = 1.0d0
         oexch = .true.
      endif
c
c coulomb
      ocoul = CD_HF_coulomb() .or. .not. odft
c
c coulomb weighs, modified later in multipole case
c
      fac1 = 1.0d0
      fac2 = 1.0d0

      pidiv4 = datan(done)
      pi = four*pidiv4
      pito52 = two*pi**twopt5

      ltri = ikyp(nshels)
      l2 = nx
      oskipp = .false.

      dlncutoff = dlog(cutoff)

      nschwz = 0
      ischw = igmem_alloc_inf(ltri,fnm,snm,'schwarz',IGMEM_DEBUG)
c
      call sinset
c ***
      if(.not.odscf) then
        call rdmake(prefac)
        dlntol=-dlog(cutoff*0.1d0)
      endif
      time = cpulft(1)
      tim0 = time
      tim1 = time
c
c *** read in ints for schwarz inequality test
      call secget(isect(421),m25,iblk25)
      call rdedx(q(ischw),ltri,iblk25,idaf)
      dlnmxs = -1.0d50
      do ii = 0, ltri-1
         dlnmxs = dmax1(dlnmxs,q(ischw+ii))
      enddo
      ist0 = ist
      jst0 = jst
      kst0 = kst
      lst0 = lst
c
c     ----- ishell -----
c
      call filmax
      do 400 ii = ist0,nshels
c
c     ----- print intermediate restart data -----
c
      if(kad(ii))400,401,401
 401  dt0 = time-tim0
      dt1 = time-tim1
      tim1 = time
      if(outv ) then
        if(odscf) then
             if(iter.le.0)
     *       write(iwr,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             if(ologf) then
             write(text,9007)ii,jst0,kst0,lst0,nrec,dt1,dt0
             call sptchk(text)
             endif
        else
               write(iwr,9008)ii,jst0,kst0,lst0,nrec,icount,dt1,dt0
        endif
       endif
c
c     ----- eliminate ishell -----
c
      do 120 it = 1,nt
      id = iso(ii,it)
      if (id .gt. ii) go to 400
      m0(it) = id
  120 continue
      ikyii=iky(ii)
c
c     ----- jshell -----
c
      j0 = jst0
      do 380 jj = j0,ii
      jst0 = 1
      if(kad(jj))380,141,141
c ***
141   itrij=ikyii+jj
      mij = itrij
      do 180 it = 1,nt
      jd = iso(jj,it)
      if (jd .gt. ii) go to 380
      id = m0(it)
      if (id .ge. jd) go to 160
      nd = id
      id = jd
      jd = nd
  160 if (id .eq. ii .and. jd .gt. jj) go to 380
      mi(it) = id
      mj(it) = jd
  180 continue
      ishell = ii
      jshell = jj
c *** reject on basis of ij 
      ijij = itrij + ischw -1
c     call scrkl(q(ischw),iso,tolij,nshels)
c     test = q(ijij) + tolij
c     if(test.lt.dlncutoff) then
c      if(odscf)intcut(1)=intcut(1)+1
      if (odscf) then
       tijkl = dlntol + q(ijij) + dlnmxs + dlnmxd
        if (tijkl .le. 0.0d0) then
         intcut(1)=intcut(1)+1
         nschwz = nschwz + itrij
        goto 380
        endif
      else
       if (q(ijij)+dlnmxs.lt.dlncutoff) then
        nschwz = nschwz + itrij
       goto 380
       endif
      endif
c
c     ----- kshell -----
c
      k0 = kst0
      do 360 kk = k0,ii
      kst0 = 1
      if(kad(kk))360,361,361
 361  do 220 it = 1,nt
      kd = iso(kk,it)
      if (kd .gt. ii) go to 360
  220 mk(it) = kd
      maxll = kk
      if (kk .eq. ii) maxll = jj
      ikykk=iky(kk)
      if (odscf) then
         itrik=ikyii+kk
         mik = itrik
         itrjk=iky(max(jj,kk))+min(jj,kk)
         mjk = itrjk
         tijk = dmax1(rdmat(mij),rdmat(mik),rdmat(mjk))
      endif
c
c     ----- lshell ----
c
      l0 = lst0
      do 340 ll = l0,maxll
      lst0 = 1
c ***
      if(kad(ll).lt.0) goto 340
c ***
      itrkl=ikykk+ll
      klkl = itrkl + ischw -1
      test = q(ijij) + q(klkl)
      if(odscf) then
         mil = ikyii+ll
         mjl = iky(max(jj,ll)) + min(jj,ll)
         mkl = itrkl
         tijkl=dlntol + test +
     +        dmax1(tijk,rdmat(mil),rdmat(mjl),rdmat(mkl))
         if(tijkl.le.0.0d0) then
            nschwz = nschwz + 1
            intcut(3)=intcut(3)+1
            goto 340
         endif
      else
         oskipp = test.lt.dlncutoff
         if(oskipp) then
           nschwz = nschwz + 1
           go to 340
         endif
      endif
      n4 = 0
      do 300 it = 1,nt
      ld = iso(ll,it)
      if (ld .gt. ii) go to 340
      kd = mk(it)
      if (kd .ge. ld) go to 260
      nd = kd
      kd = ld
      ld = nd
  260 id = mi(it)
      jd = mj(it)
      if (id .ne. ii .and. kd .ne. ii) go to 300
      if (kd .lt. id) go to 280
      if (kd .eq. id .and. ld .le. jd) go to 280
      nd = id
      id = kd
      kd = nd
      nd = jd
      jd = ld
      ld = nd
  280 if (jd .lt. jj) go to 300
      if (jd .gt. jj) go to 340
      if (kd .lt. kk) go to 300
      if (kd .gt. kk) go to 340
      if (ld .lt. ll) go to 300
      if (ld .gt. ll) go to 340
      n4 = n4+1
  300 continue
      q4 =  dble(nt)/ dble(n4)
c
c     ----- (ii,jj//kk,ll) -----
c
      ishell = ii
      jshell = jj
      kshell = kk
      lshell = ll
      qq4 = q4
c
c     ----- initialize gout to zero -----
c     ----- compute two-electron integrals ----
c     ----- write them out on mainfile -----
c
      if (odscf) then
        if (zscftp.eq.'uhf')then
          call genr70(gout,1,.false.)
          call dir_build_uhf70(fock,fockb,
     +         dens,densb,gout,
     +         fac1,fac2, facex, ocoul, oexch)
        else if (zscftp.eq.'gvb'.or.zscftp.eq.'grhf') then
          call genr70(gout,1,.false.)
          if(nsheld.le.1) then
            call dir_build_open_70(fock,exch,dens,gout)
          else
            call dir_build_open2_70(l2,fock,exch,dens,gout)
          endif
        else
          if(omorok) then
            call genr70(gout,1,.false.)
            call dbuild70_morok(fock,dens,gout)
          else
            call genr70(gout,1,.false.)
            call dbuild70(fock,dens,gout,
     &           fac1, fac2, facex, ocoul, oexch)
          endif
        endif
      else
        call genr70(gout,1,.false.)
        call qout70(gout)
      endif
c
c
c     ----- check cpu time/ maxblock condition -----
      call chkout(ii,jj,kk,ll,fock,q)
      if(omaxb.or.tim.gt.timlim)go to 420
c
  340 continue
  360 continue
  380 continue
      time = cpulft(1)
  400 continue
      ist = 1
      jst = 1
      kst = 1
      lst = 1
      call final(q,fock,dens)
c
c     ----- reset core memory
c
*IJB Sum nschwz
c$$$  420 if(outvv) write(iwr,6030) nschwz
 420  call pg_igop( 160467, nschwz, 1, '+' )
      if(outvv) write(iwr,6030) nschwz

      call gmem_free_inf(ischw,fnm,snm,'schwarz')
      return
 9008 format(i4,3i5,1x,i10,i9,f11.2,f9.2)
 9007 format(i4,3i5,1x,i10,9x,  f11.2,f9.2)
 6030 format(/1x,'schwarz inequality test skipped ',i10,
     + ' integral blocks')
      end
      subroutine scrkl(schw,iso,scrnij,nshels)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c 
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      dimension schw(*),iso(nshels,*)
c
      scrnij = -1.0d50
      do 360 kk = 1,ishell
      if(kad(kk).ge.0) then
       do 220 it = 1,nt
       kd = iso(kk,it)
       if (kd .gt. ishell) go to 360
  220  maxll = kk
       if (kk .eq. ishell) maxll = jshell
       ikykk=iky(kk)
       do 340 ll = 1,maxll
       if(kad(ll).ge.0) then
        itrkl=ikykk+ll
        scrnij = max(scrnij,schw(itrkl))
       endif
  340  continue
      endif
  360 continue
      return
      end
      subroutine sp0000(gout)
c        *****  special fast routine for -p- loop for 0000 ****
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
      common/astore/qq,theta,n
      common/inttab/
     +  a0(333),b0(333),c0(333),abc(5001)
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
      common/miscg/mab,mcd,ngangb
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
c
      dimension gout(*)
c
      data sixty/60.0d0/
c
      g0000 = 0.d0
      do 220 k = 1,ngc
      gc = cgg(k)
      csck = csc(k)
      gcrcds = gc*rcdsq
      do 220 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.d0/gcd
      gdecd = gd*ecd
      ppp = -gdecd*gcrcds
      if(ppp+sixty)220,500,500
  500 pp = ecd* dexp(ppp)
  520 qqtest = pp*cmaxc(k)*cmaxd(l)
      if (qqtest .le. error1) go to 100
      ismlq = 0
      go to 120
  100 if (qqtest .le. error2) go to 220
      ismlq = 1
  120 cq = gdecd*rcd
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      h0000 = 0.d0
      do 200 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 200
      auxvar = var(isml+1)
      p = ((aqz-app(i))**2+qperp2)/(ep(i)+ecd)
      if (p .le. auxvar) go to 180
      h0000 = h0000+dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      go to 200
  180 continue
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2)
      h0000 = h0000+ dp00p(i)*f0/dsqrt(gp(i)+gcd)
  200 continue
      g0000 = g0000+h0000*csck*csd(l)*pp
  220 continue
      gout(1) = g0000
      return
      end
      subroutine filmax
c
c     --------------------------
c     --------------------------
c
c     finds maximum value of s and p coefficients
c     also sets limits determining how accurately a set of integrals
c     need be evaluated in order to guarantee an overall integral
c     accuracy of 10**-6
c
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
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
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
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      data five/5.0d0/,fiften/15.0d0/
      data ten,twenty/10.0d0,20.0d0/
c
      do 100 i = 1,nshell
      l = kstart(i)
      n = l+kng(i)-1
      do 100 j = l,n
      a1 =  dabs(cs(j))
      a2 =  dabs(cp(j))
      cmax(j) = dmax1(a1,a2)
  100 continue
c
      if (ifasti.eq.1) then
c     settings in gamess-us
       var(1) = twenty
       var(2) = ten
       error1 = dmin1(cutoff,1.0d-10)
       error2 = dmin1(cutoff,1.0d-10)
      else if (ifasti.eq.2) then
c     tighten thresholds
       var(1) = fiften
       var(2) = five
       error1 = dsqrt(cutoff) * 0.01d0
       error2 = cutoff * 0.001d0
      else
       error1 = dsqrt(cutoff) * 0.1d0
       var(1) = fiften
       var(2) = five
       error2 = cutoff * 0.01d0
      endif
      return
      end
      subroutine qout70(gout)
c
c ... qout70 now handles conventional integrals only
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
      common/craypk/integ(1)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
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
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      common/blkin/goutx(510),nword
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
      dimension gout(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c ***
c     data pt5 /0.5/
c
c     ----- pack the 4 indices of integral into one word
c     ----- write label + integral on mainfile
c
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      ijn = 0
      jmax = maxj
        do 2600 i = 1,maxi
          if (oianj) jmax = i
          i1 = loci + i
          ipack = i4096(i1)
          do 2400 j = 1,jmax
            ijn = ijn+1
            n1 = ib(ib1,i)+ib(jb1,j)+1
            lmax = maxl
            i2 = locj + j
            if(i1-i2)10, 20, 20
  10        lab1=i4096(i2) + i1
            go to 30
  20        lab1=ipack + i2
  30        kln = 0
            do 2200 k = 1,maxk
              if (okanl) lmax = k
              i3 = lock+k
              kpack = i4096(i3)
              do 2000 l = 1,lmax
                kln = kln+1
                if (oident .and. kln .gt. ijn) go to 2400
                nn = n1+ib(kb1,k)+ib(lb1,l)
                val = gout(nn)
                if ( dabs(val) .lt. cutoff) go to 2000
                i4 = locl+l
                if(i3-i4) 40, 50, 50
 40             lab2=i4096(i4) + i3
                go to 60
 50             lab2=kpack + i4
 60             goutx(icount) = val
                integ(ic4  )= max(lab1,lab2)
                integ(ic4+1)= min(lab1,lab2)
                ic4 = ic4 + 2
                icount = icount+1
                if (icount .le. nintmx) go to 2000
                call blocki
                if(omaxb)go to 261
2000          continue
2200        continue
2400      continue
2600    continue
  261 return
      end
      subroutine dbuild70(fock,dmat,gout,
     &     fac1,fac2, facex, ocoul, oexch)

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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
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
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c

      dimension ib(4,4)
      dimension fock(*),dmat(*),gout(*)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
c
c Coulomb terms
c
      if(ocoul)then
         val2=val+val
         val4=val2+val2
         f12 = facij*val4*dmat(itr34) + fock(itr12)
         fock(itr12) = f12
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
         endif
      endif
c
      if(oexch)then
c     
c Full exchange term for HF or weighted exchange term for b3lyp etc
c     
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end

      subroutine dir_build_uhf70(fock,ak,p,q,gout,
     &   fac1, fac2, facex, ocoul, oexch)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
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
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),fock(*),ak(*),p(*),q(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell UHF Fock builder 
c
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
      gik=val
      val2=gik+gik
      val4=val2+val2
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
c
c Coulomb term
c
      if(ocoul)then
         gik=val
         val2=gik+gik
         val4=val2+val2
         fock(itr12) = facij*val4*p(itr34) + fock(itr12)
         if(itr12 .ne. itr34)then
            fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
         endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
         gik=val*facex
         val2=gik+gik
         val4=val2+val2
         gil=gik
         if(i1.eq.i3.or.i2.eq.i4)gik=val2
         if(i2.eq.i3)gil=val2
         if(i2.ge.i3)goto 1
         itr23=ikyk+i2
         if(i2.ge.i4)goto 1
         itr24=iky(i4)+i2
 1       ajk=fock(itr23)-gil*p(itr14)
         bjk=ak(itr23)+gil*q(itr14)
         ail=fock(itr14)-gil*p(itr23)
         bil=ak(itr14)+gil*q(itr23)
         aik=fock(itr13)-gik*p(itr24)
         bik=ak(itr13)+gik*q(itr24)
         fock(itr24)=fock(itr24)-gik*p(itr13)
         ak(itr24)=ak(itr24)+gik*q(itr13)
         fock(itr23)=ajk
         ak(itr23)=bjk
         fock(itr14)=ail
         ak(itr14)=bil
         fock(itr13)=aik
         ak(itr13)=bik
      endif
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open_70(coul,exch,dens,gout)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
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
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell SCF J & K builder (nshell=1)
c
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        gik=val
        val2=val+val
        ikyi=iky(i1)
        ikyj=iky(i2)
        ikyk=iky(i3)
        itr13=ikyi+i3
        itr14=ikyi+i4
        itr12=ikyi+i2
        itr23=ikyj+i3
        itr24=ikyj+i4
        itr34=ikyk+i4
        gil=val
        if(i1.eq.i3.or.i2.eq.i4)gik=val2
        if(i2.eq.i3)gil=val2
        if(i2.ge.i3)goto 280
        itr23=ikyk+i2
        if(i2.ge.i4)goto 280
        itr24=iky(i4)+i2
  280   bij=val2*dens(itr34)+coul(itr12)
        coul(itr34)=val2*dens(itr12)+coul(itr34)
        coul(itr12)=bij
        bjk=exch(itr23)+gil*dens(itr14)
        bil=exch(itr14)+gil*dens(itr23)
        bik=exch(itr13)+gik*dens(itr24)
        exch(itr24)=exch(itr24)+gik*dens(itr13)
        exch(itr23)=bjk
        exch(itr14)=bil
        exch(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open2_70(l2,coul,exch,dens,gout)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
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
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension gout(*),coul(*),exch(*),dens(*)
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
c     ----- Direct open-shell SCF J & K builder (nshell > 1)
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=val+val
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      gil=val
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)go to 280
      itr23=ikyk+i2
      if(i2.ge.i4)go to 280
      itr24=iky(i4)+i2
  280 continue
      do 300 iiii=1,nsheld
      bij=val2*dens(itr34)+coul(itr12)
      coul(itr34)=val2*dens(itr12)+coul(itr34)
      coul(itr12)=bij
      bjk=exch(itr23)+gil*dens(itr14)
      bil=exch(itr14)+gil*dens(itr23)
      bik=exch(itr13)+gik*dens(itr24)
      exch(itr24)=exch(itr24)+gik*dens(itr13)
      exch(itr23)=bjk
      exch(itr14)=bil
      exch(itr13)=bik
      itr12=itr12+l2
      itr34=itr34+l2
      itr13=itr13+l2
      itr14=itr14+l2
      itr23=itr23+l2
      itr24=itr24+l2
  300 continue
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine sp0001(gout)
c        *****  special fast routine for -p- loop for 0001 ****
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
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
      common/astore/qq,theta,n
      common/inttab/
     +  a0(333),b0(333),c0(334),a1(333),b1(333),c1(334),
     +  a2(4000)
      dimension gout(*)
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
      common/miscg/mab,mcd,ngangb
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
      data sixty,tenm12/60.0d0,1.0d-12/
c
      gout1 = 0.0d0
      gout2 = 0.0d0
      gout3 = 0.0d0
      gout4 = 0.0d0
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = 1.0d0/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = 1.0d0
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      g = 1.d0/(ep(i)+ecd)
      p = (pqab*pqab+qperp2)*g
      if (p .le. auxvar) go to 140
      q0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      q1 = 0.5d0*q0/p
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      q0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      q1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*q
  160 u = g*q1
      h0000 = h0000+q0
      h0001 = h0001+u
      h0003 = h0003-u*pqab
  180 continue
      h0001 = h0001*ecd*qperp
      h0003 = h0003*ecd
      p = dq*h0000
      g0001 = h0001*cosp+p*sing
      g0002 = h0001*sinp
      g0003 = h0003+p*cosg
      gout1 = gout1+dq00*h0000
      gout2 = gout2+dq01*g0001
      gout3 = gout3+dq01*g0002
      gout4 = gout4+dq01*g0003
 940  continue
      t1 = gout2
      t2 = gout3
      t3 = gout4
      gout(1) = gout1
      gout(2) = p11*t1+p21*t2+p31*t3
      gout(3) = p12*t1+p22*t2+p32*t3
      gout(4) = p13*t1+p23*t2+p33*t3
      return
      end
      subroutine sp0101(gout)
c        *****  special fast routine for -p- loop for 0101 *****
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
      common/astore/qq,theta,n
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      dimension gout(*)
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
      common/miscg/mab,mcd,ngangb
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
c
      real*8 conp
c
c     conp = prefactors from pairs of primitives
c
      common /const/ conp(mxprms*mxprms)
c
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0100 = 0.d0
      h0101 = 0.d0
      h0103 = 0.d0
      h0300 = 0.d0
      h0301 = 0.d0
      h0303 = 0.d0
c        *****  begin -p- loop                   *****
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = (pqab2+qperp2)*g
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 continue
      g03 = -pqab*f1
      h0000 = h0000+f0 *dp00
      h0001 = h0001+f1 *dp00
      h0003 = h0003+g03*dp00
      h0100 = h0100-f1
      h0101 = h0101-f2
      h0103 = h0103+pqab*f2
      h0300 = h0300-g03+bp*f0
      h0301 = h0301+bp*f1
      h0303 = h0303-pqab2*f2+bp*g03
  180 continue
      p = qperp*ecd
      h0001 = h0001*p
      h0003 = h0003*ecd
      h0202 = -0.5d0*ecd*h0100
      h0100 = h0100*qperp
      h0101 = h0101*qperp2*ecd
      h0103 = h0103*p
      h0301 = h0301*p
      h0303 = h0303*ecd
      h0301 = h0301+h0103
      h0101 = h0101+h0202
      h0303 = h0303+h0202
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      g0101 = cosp*h0101
      g0102 = sinp*h0101
      g0201 = u12*h0202
      g0202 = cosp*h0202
      g0301 = cosp*h0301
      g0302 = sinp*h0301
      g0303 = h0303
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0300 = h0300
      g0000 = h0000
      h0101 = g0101
      h0102 = g0102
      h0201 = g0201
      h0202 = g0202
      g0101 = cosp*h0101+u12*h0201
      g0102 = cosp*h0102+u12*h0202
      g0103 = cosp*h0103
      g0201 = sinp*h0101+cosp*h0201
      g0202 = sinp*h0102+cosp*h0202
      g0203 = sinp*h0103
      g0100 = cosp*h0100
      g0200 = sinp*h0100
      go to 2000
  920 g0100 = h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      go to 2000
 1000 g0100 = -h0100
      g0101 = h0101
      g0102 = dzero
      g0103 = -h0103
      g0200 = dzero
      g0201 = dzero
      g0202 = h0202
      g0203 = dzero
      g0300 = h0300
      g0301 = -h0301
      g0302 = dzero
      g0303 = h0303
      g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
2000  r14 = dq*sing
      r34 = dq*cosg
      if (rcdsq) 720,720,700
  700 g0001 = g0001+r14*g0000
      g0101 = g0101+r14*g0100
      g0201 = g0201+r14*g0200
      g0301 = g0301+r14*g0300
      g0003 = g0003+r34*g0000
      g0103 = g0103+r34*g0100
      g0203 = g0203+r34*g0200
      g0303 = g0303+r34*g0300
  720 gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 17) = gout( 17)+g0100*dq00
      gout( 18) = gout( 18)+g0101*dq01
      gout( 19) = gout( 19)+g0102*dq01
      gout( 20) = gout( 20)+g0103*dq01
      gout( 33) = gout( 33)+g0200*dq00
      gout( 34) = gout( 34)+g0201*dq01
      gout( 35) = gout( 35)+g0202*dq01
      gout( 36) = gout( 36)+g0203*dq01
      gout( 49) = gout( 49)+g0300*dq00
      gout( 50) = gout( 50)+g0301*dq01
      gout( 51) = gout( 51)+g0302*dq01
      gout( 52) = gout( 52)+g0303*dq01
 940  continue
c ***
c ***
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 0101
c
c
c
      ind = 0
      do 1 loopkl = 1,16
      ind = ind+1
      i1 = 16+ind
      i2 = 32+ind
      i3 = 48+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   1  continue
      ind = -15
      do 2 j = 1,4
      ind = ind+16
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   2  continue
      return
      end
      subroutine sp0111(gout)
c        *****  special fast routine for -p- loop for 0111 *****
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
      dimension g(64),h(64)
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
      common/miscg/mab,mcd,ngangb
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
      real*8 conp
c
c     conp = prefactors from pairs of primitives
c
      common /const/ conp(mxprms*mxprms)
c
      common/astore/qq,theta,n
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
c
      dimension gout(*)
c
      data dzero/0.0d0/,done/1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
  660 p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
c        *****  begin -p- loop                   *****
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      eab = ep(i)
      dp00 = dp00p(i)
      bp = bpp(i)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(pqab2+qperp2)
      if (p .le. auxvar) go to 140
      f0 = conp(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      go to 160
  140 q = conp(i)/dsqrt(gp(i)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
      f3 = (a3(n+1)+theta*b3(n+1)-theta3*c3(n+1)+theta4*c3(n+2))*gggy
  160 continue
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f2pqa2 = f2*pqab2
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      q1 = q1+f0 *bp
      q2 = q2+f1 *bp
      q3 = q3+f2 *bp
      q4 = q4+f1pqab*bp
      q5 = q5+f2pqab*bp
      q6 = q6+f2pqa2*bp
      r1 = r1+f1
      r2 = r2+f2
      r3 = r3+f3
      r4 = r4+f1pqab
      r5 = r5+f2pqab
      r6 = r6+f3pqab
      r7 = r7+f2pqa2
      r8 = r8+f3*pqab2
      r9 = r9+f3pqab*pqab2
  180 continue
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      qecd = qperp*ecd
      qecd2 = qperp*ecd2
      q2ecd = qperp2*ecd
      q2ecd2 = qperp2*ecd2
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = 0.5d0*qecd2*r2
      h( 44) = -0.5d0*ecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q2ecd2*qperp*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
  120 u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      g( 22) = v44*h( 22)+v47*h( 27)
      g( 23) = v54*h( 22)+v57*h( 27)
      g( 27) = v74*h( 22)+v77*h( 27)
      g( 24) = cosp*h( 24)
      g( 28) = sinp*h( 24)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g( 64) = h( 64)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 16) = h( 16)
      g( 18) = cosp*h( 18)
      g( 19) = sinp*h( 18)
      g( 20) = h( 20)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g( 52) = h( 52)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g(  4) = h(  4)
      g( 49) = h( 49)
      g(  1) = h(  1)
      h( 22) = g( 22)
      h( 23) = g( 23)
      h( 24) = g( 24)
      h( 27) = g( 27)
      h( 28) = g( 28)
      h( 38) = g( 38)
      h( 39) = g( 39)
      h( 40) = g( 40)
      h( 43) = g( 43)
      h( 44) = g( 44)
      h( 18) = g( 18)
      h( 19) = g( 19)
      h( 34) = g( 34)
      h( 35) = g( 35)
      g( 22) = cosp*h( 22)+u12*h( 38)
      g( 23) = cosp*h( 23)+u12*h( 39)
      g( 24) = cosp*h( 24)+u12*h( 40)
      g( 27) = cosp*h( 27)+u12*h( 43)
      g( 28) = cosp*h( 28)+u12*h( 44)
      g( 32) = cosp*h( 32)
      g( 38) = sinp*h( 22)+cosp*h( 38)
      g( 39) = sinp*h( 23)+cosp*h( 39)
      g( 40) = sinp*h( 24)+cosp*h( 40)
      g( 43) = sinp*h( 27)+cosp*h( 43)
      g( 44) = sinp*h( 28)+cosp*h( 44)
      g( 48) = sinp*h( 32)
      g( 18) = cosp*h( 18)+u12*h( 34)
      g( 19) = cosp*h( 19)+u12*h( 35)
      g( 20) = cosp*h( 20)
      g( 34) = sinp*h( 18)+cosp*h( 34)
      g( 35) = sinp*h( 19)+cosp*h( 35)
      g( 36) = sinp*h( 20)
      g( 17) = cosp*h( 17)
      g( 33) = sinp*h( 17)
      go to 2000
  920 g( 17) = h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = h( 20)
      g( 22) = h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = h( 27)
      g( 28) = dzero
      g( 32) = h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
      go to 2000
 1000 g( 17) = -h( 17)
      g( 18) = h( 18)
      g( 19) = dzero
      g( 20) = -h( 20)
      g( 22) = -h( 22)
      g( 23) = dzero
      g( 24) = h( 24)
      g( 27) = -h( 27)
      g( 28) = dzero
      g( 32) = -h( 32)
      g( 33) = dzero
      g( 34) = dzero
      g( 35) = h( 35)
      g( 36) = dzero
      g( 38) = dzero
      g( 39) = -h( 39)
      g( 40) = dzero
      g( 43) = dzero
      g( 44) = h( 44)
      g( 48) = dzero
      g( 49) = h( 49)
      g( 50) = -h( 50)
      g( 51) = dzero
      g( 52) = h( 52)
      g( 54) = h( 54)
      g( 55) = dzero
      g( 56) = -h( 56)
      g( 59) = h( 59)
      g( 60) = dzero
      g( 64) = h( 64)
      g(  1) = h(  1)
      g(  2) = -h(  2)
      g(  3) = dzero
      g(  4) = h(  4)
      g(  6) = h(  6)
      g(  7) = dzero
      g(  8) = -h(  8)
      g( 11) = h( 11)
      g( 12) = dzero
      g( 16) = h( 16)
2000  r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,50,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 720,720,700
700   do 701 kq1=1,49,16
          t1=g(kq1)
          t5=g(kq1+4)
          t13=g(kq1+12)
          t5=t5+r13*t1
          t13=t13+r33*t1
          t2=g(kq1+1)
          g(kq1+5)=g(kq1+5)+r13*t2+r14*t5
          g(kq1+13)=g(kq1+13)+r33*t2+r14*t13
          g(kq1+1)=t2+r14*t1
          t3=g(kq1+2)
          g(kq1+6)=g(kq1+6)+t3*r13
          g(kq1+14)=g(kq1+14)+t3*r33
          t9=g(kq1+8)
          g(kq1+9)=g(kq1+9)+r14*t9
          g(kq1+11)=g(kq1+11)+r34*t9
          t4=g(kq1+3)
          g(kq1+7)=g(kq1+7)+t4*r13+t5*r34
          g(kq1+4)=t5
          g(kq1+15)=g(kq1+15)+t4*r33+t13*r34
          g(kq1+12)=t13
701       g(kq1+3)=t4+t1*r34
720   do 721 kq1=2,62,4
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq11
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq11
721       gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq11
      dq01dd=dq01-dq11
      do 722 kq1=1,49,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq00
          gout(kq1+1)=gout(kq1+1)+g(kq1+1)*dq01dd
          gout(kq1+2)=gout(kq1+2)+g(kq1+2)*dq01dd
722       gout(kq1+3)=gout(kq1+3)+g(kq1+3)*dq01dd
      do 723 kq1=5,53,16
          gout(kq1  )=gout(kq1  )+g(kq1  )*dq10
          gout(kq1+4)=gout(kq1+4)+g(kq1+4)*dq10
723       gout(kq1+8)=gout(kq1+8)+g(kq1+8)*dq10
  940 continue
      ind = 0
      do 1 loopkl = 1,16
      ind = ind+1
      i1 = 16+ind
      i2 = 32+ind
      i3 = 48+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   1  continue
      ind = -12
      do 2 j = 1,4
      ind = ind+12
      do 2 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   2  continue
      ind = -3
      do 3 loopjk = 1,16
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
   3  continue
      return
      end
c ******************************************************
c ******************************************************
c             =   sp1111  =
c ******************************************************
c ******************************************************
      subroutine sp1111(gout)
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
      parameter (dzero=0.0d0)
      dimension g(256),h(256)
      dimension kq1off(10),kq2off(6),kq3off(6),kq4off(4),kq5off(6)
      common/miscg/mab,mcd,ngangb
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
c
      real*8 conp
c
c     conp = prefactors from pairs of primitives
c
      common /const/ conp(mxprms*mxprms)
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
      common/astore/qq,theta,n
      common/inttab/
     +  aa(333),ba(333),ca(333),abc1,
     +  ab(333),bb(333),cb(333),abc2,
     +  ac(333),bc(333),cc(333),abc3,
     +  ad(333),bd(333),cd(333),abc4,
     +  ae(333),be(333),ce(333),abc5,
     +  af(333),bf(333),cf(333),abc6
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
c
      dimension gout(*)
c      data dzero/0.0e0/,done/1.0e0/
      data done/1.0d0/
      data kq1off/1,17,49,65,81,113,161,193,209,241/
      data kq2off/33,97,129,145,177,225/
      data kq3off/1,49,81,161,193,241/
      data kq4off/17,65,113,209/
      data kq5off/33,97,129,145,177,225/
      data sixty,tenm12/60.0d0,1.0d-12/
c ***
c *** this is the fps version of sp1111.
c ***
c *** as much code as possible reduced to loops (>=4)
c *** to avoid ps cache misses and to enhance compiler
c *** optimisation. will probably run like a drain on
c *** the cray-1s]
c ***
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
660   p1 = 0.d0
      p2 = 0.d0
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      q1 = 0.d0
      q2 = 0.d0
      q3 = 0.d0
      q4 = 0.d0
      q5 = 0.d0
      q6 = 0.d0
      r1 = 0.d0
      r2 = 0.d0
      r3 = 0.d0
      r4 = 0.d0
      r5 = 0.d0
      r6 = 0.d0
      r7 = 0.d0
      r8 = 0.d0
      r9 = 0.d0
      v1 = 0.d0
      v2 = 0.d0
      v3 = 0.d0
      v4 = 0.d0
      v5 = 0.d0
      v6 = 0.d0
      w1 = 0.d0
      w2 = 0.d0
      w3 = 0.d0
      w4 = 0.d0
      w5 = 0.d0
      w6 = 0.d0
      w7 = 0.d0
      w8 = 0.d0
      w9 = 0.d0
      s1 = 0.d0
      s2 = 0.d0
      s3 = 0.d0
      s4 = 0.d0
      s6 = 0.d0
      s7 = 0.d0
      s8 = 0.d0
      s9 = 0.d0
      s10 = 0.d0
      s11 = 0.d0
      s12 = 0.d0
      s13 = 0.d0
      s14 = 0.d0
      t1 = 0.d0
      t2 = 0.d0
      t3 = 0.d0
      t4 = 0.d0
      t5 = 0.d0
      t6 = 0.d0
      t7 = 0.d0
      t8 = 0.d0
      t9 = 0.d0
      t10 = 0.d0
      t11 = 0.d0
      t12 = 0.d0
      t13 = 0.d0
      t14 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      c3 = 0.d0
      c4 = 0.d0
      c5 = 0.d0
      c6 = 0.d0
      do 180 ind = 1,ngangb
      isml = ismlq+ismlp(ind)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
       eab = ep(ind)
      dp00 = dp00p(ind)
      dp01 = dp01p(ind)
      dp10 = dp10p(ind)
      ap = app(ind)
      bp = bpp(ind)
      pqab = aqz-ap
      pqab2 = pqab*pqab
      gabcd = 1.d0/(eab+ecd)
      p = gabcd*(qperp2+pqab2)
      if (p .le. auxvar) go to 140
      f0 = dsqrt(0.7853981625d0/(p*(gp(ind)+gcd)))*conp(ind)
      gtx = gabcd/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      f3 = 2.5d0*f2*gtx
      f4 = 3.5d0*f3*gtx
      go to 160
  140 q = conp(ind)/dsqrt(gp(ind)+gcd)
      gy = gabcd*q
      ggy = gabcd*gy
      gggy = gabcd*ggy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (aa(n+1)+theta*ba(n+1)-theta3*ca(n+1)+theta4*ca(n+2))*q
      f1 = (ab(n+1)+theta*bb(n+1)-theta3*cb(n+1)+theta4*cb(n+2))*gy
      f2 = (ac(n+1)+theta*bc(n+1)-theta3*cc(n+1)+theta4*cc(n+2))*ggy
      f3 = (ad(n+1)+theta*bd(n+1)-theta3*cd(n+1)+theta4*cd(n+2))*gggy
      f4 = (ae(n+1)+theta*be(n+1)-theta3*ce(n+1)+theta4*ce(n+2))*gggy*
     &     gabcd
  160 apbp = ap*bp
      eab2 = eab*eab
      bpdp01 = bp*dp01
      apdp10 = ap*dp10
      edp01 = eab*dp01
      edp10 = eab*dp10
      f1pqab = f1*pqab
      f2pqab = f2*pqab
      f3pqab = f3*pqab
      f4pqab = f4*pqab
      f1pqa2 = f1*pqab2
      f2pqa2 = f2*pqab2
      f3pqa2 = f3*pqab2
      f4pqa2 = f4*pqab2
      f2pqa3 = f2pqa2*pqab
      f3pqa3 = f3pqa2*pqab
      f4pqa3 = f4pqa2*pqab
      p1 = p1+f0 *dp00
      p2 = p2+f1 *dp00
      p3 = p3+f2 *dp00
      p4 = p4+f1pqab*dp00
      p5 = p5+f2pqab*dp00
      p6 = p6+f2pqa2*dp00
      r1 = r1+f1 *edp01
      r2 = r2+f2 *edp01
      r3 = r3+f3 *edp01
      r4 = r4+f1pqab *edp01
      r5 = r5+f2pqab *edp01
      r6 = r6+f3pqab *edp01
      r7 = r7+f2pqa2 *edp01
      r8 = r8+f3pqa2 *edp01
      r9 = r9+f3pqa3 *edp01
      w1 = w1+f1 *edp10
      w2 = w2+f2 *edp10
      w3 = w3+f3 *edp10
      w4 = w4+f1pqab *edp10
      w5 = w5+f2pqab *edp10
      w6 = w6+f3pqab *edp10
      w7 = w7+f2pqa2 *edp10
      w8 = w8+f3pqa2 *edp10
      w9 = w9+f3pqa3 *edp10
      s1 = s1+f0 *eab
      s2 = s2+f1 *eab
      s3 = s3+f2 *eab
      s4 = s4+f3 *eab
      s6 = s6+f1pqab*eab
      s7 = s7+f2pqab*eab
      s8 = s8+f3pqab*eab
      s9 = s9+f1pqa2*eab
      s10 = s10+f2pqa2*eab
      s11 = s11+f3pqa2*eab
      s12 = s12+f2pqa3*eab
      s13 = s13+f3pqa3*eab
      s14 = s14+f3pqa3*pqab*eab
      t1 = t1+f0 *eab2
      t2 = t2+f1 *eab2
      t3 = t3+f2 *eab2
      t4 = t4+f3 *eab2
      t5 = t5+f4 *eab2
      t6 = t6+f2pqab*eab2
      t7 = t7+f3pqab*eab2
      t8 = t8+f4pqab*eab2
      t9 = t9+f2pqa2*eab2
      t10 = t10+f3pqa2*eab2
      t11 = t11+f4pqa2*eab2
      t12 = t12+f3pqa3*eab2
      t13 = t13+f4pqa3*eab2
      t14 = t14+f4pqa3*pqab*eab2
      if (rabsq .eq. 0.0d0) go to 180
      q1 = q1+f0 *bpdp01
      q2 = q2+f1 *bpdp01
      q3 = q3+f2 *bpdp01
      q4 = q4+f1pqab*bpdp01
      q5 = q5+f2pqab*bpdp01
      q6 = q6+f2pqa2*bpdp01
      v1 = v1+f0 *apdp10
      v2 = v2+f1 *apdp10
      v3 = v3+f2 *apdp10
      v4 = v4+f1pqab*apdp10
      v5 = v5+f2pqab*apdp10
      v6 = v6+f2pqa2*apdp10
      c1 = c1+f0 *apbp
      c2 = c2+f1 *apbp
      c3 = c3+f2 *apbp
      c4 = c4+f1pqab*apbp
      c5 = c5+f2pqab*apbp
      c6 = c6+f2pqa2*apbp
  180 continue
      a1 = aqz*s2-s6
      a2 = aqz*s3-s7
      a3 = aqz*s4-s8
      a4 = aqz*s6-s9
      a5 = aqz*s7-s10
      a6 = aqz*s8-s11
      a8 = aqz*s10-s12
      a9 = aqz*s11-s13
      a10 = aqz*s13-s14
      bqz = aqz-rab
      b1 = bqz*s2-s6
      b2 = bqz*s3-s7
      b3 = bqz*s4-s8
      b4 = bqz*s6-s9
      b5 = bqz*s7-s10
      b6 = bqz*s8-s11
      b8 = bqz*s10-s12
      b9 = bqz*s11-s13
      b10 = bqz*s13-s14
      hecd = 0.5d0*ecd
      ecd2 = ecd*ecd
      hecd2 = 0.5d0*ecd2
      qecd = qperp*ecd
      hqecd = 0.5d0*qecd
      qecd2 = qperp*ecd2
      hqecd2 = 0.5d0*qecd2
      q2ecd = qperp2*ecd
      q3ecd = qperp*q2ecd
      q2ecd2 = qperp2*ecd2
      q3ecd2 = q2ecd2*qperp
      h(  1) = p1
      h(  2) = qecd*p2
      h(  4) = -ecd*p4
      h( 11) = hecd*(p1-ecd*p2)
      h(  6) = h( 11)+q2ecd2*p3
      h(  8) = -qecd2*p5
      h( 16) = h( 11)+ecd2*p6
      h( 17) = -qperp*r1
      h( 49) = r4+q1
      h( 35) = hecd*r1
      h( 18) = h( 35)-q2ecd*r2
      h( 20) = qecd*r5
      h( 50) = h( 20)+qecd*q2
      h( 52) = h( 35)-ecd*r7-ecd*q4
      h( 39) = hqecd2*r2
      h( 44) = -hecd2*r5
      h( 27) = h( 39)-qperp*h( 35)
      h( 59) = h( 44)+hecd*(h( 49)-ecd*q2)
      h( 24) = h( 44)+q2ecd2*r6
      h( 56) = h( 39)-qecd2*(r8+q5)
      h( 22) = h( 27)+h( 39)+h( 39)-q3ecd2*r3
      h( 32) = h( 27)-qecd2*r8
      h( 54) = h( 59)+q2ecd2*(r6+q3)
      h( 64) = h( 59)+h( 44)+h( 44)+ecd2*(r9+q6)
      h( 65) = -qperp*w1
      h(193) = w4+v1
      h(131) = hecd*w1
      h( 66) = h(131)-q2ecd*w2
      h( 68) = qecd*w5
      h(194) = h( 68)+qecd*v2
      h(196) = h(131)-ecd*w7-ecd*v4
      h(135) = hqecd2*w2
      h(140) = -hecd2*w5
      h( 75) = h(135)-qperp*h(131)
      h(203) = h(140)+hecd*(h(193)-ecd*v2)
      h( 72) = h(140)+q2ecd2*w6
      h(200) = h(135)-qecd2*(w8+v5)
      h( 70) = h( 75)+h(135)+h(135)-q3ecd2*w3
      h( 80) = h( 75)-qecd2*w8
      h(198) = h(203)+q2ecd2*(w6+v3)
      h(208) = h(203)+h(140)+h(140)+ecd2*(w9+v6)
      h(161) = 0.5d0*(s1-t2)
      h( 81) = h(161)+qperp2*t3
      h(113) = -qperp*(t6+b1)
      h(209) = -qperp*(t6+a1)
      h(241) = h(161)+t9+a4+b4+c1
      h(162) = hqecd*(s2-t3)
      h( 82) = h(162)-qecd*t3+q3ecd*t4
      temp = hecd*t6-q2ecd*t7
      h(114) = temp+hecd*b1-q2ecd*b2
      h(210) = temp+hecd*a1-q2ecd*a2
      h(242) = h(162)+qecd*(t10+a5+b5+c2)
      h( 99) = -hqecd*t3
      h(147) = h( 99)
      h(179) = hecd*(t6+b1)
      h(227) = hecd*(t6+a1)
      h(164) = hecd*(t6-s6)
      h( 84) = h(164)-q2ecd*t7
      temp = -hqecd*t3+qecd*t10
      h(116) = temp+qecd*b5
      h(212) = temp+qecd*a5
      h(244) = h(164)+ecd*(t6-t12-a8-b8-c4)+hecd*(a1+b1)
      h(103) = 0.25d0*ecd2*t3-0.5d0*q2ecd2*t4
      h(151) = h(103)
      h(183) = hqecd2*(t7+b2)
      h(231) = hqecd2*(t7+a2)
      h(108) = hqecd2*t7
      h(156) = h(108)
      h(188) = hecd2*(0.5d0*t3-t10-b5)
      h(236) = hecd2*(0.5d0*t3-t10-a5)
      hxxyy = 0.25d0*(ecd*(s1-t2)-ecd2*(s2-t3))
      h(171) = hxxyy+hecd2*t3
      h( 91) = hxxyy+0.5d0*(q2ecd*t3-q2ecd2*t4)
      temp = hqecd*(ecd*t7-t6)
      h(123) = temp+hqecd*(ecd*b2-b1)
      h(219) = temp+hqecd*(ecd*a2-a1)
      h(251) = hxxyy+hecd*(t9+a4+b4+c1)-hecd2*(t10+a5+b5+c2)
      h(166) = hxxyy+0.5d0*q2ecd2*(s3-t4)
      h( 86) = hxxyy+(hecd2+0.5d0*q2ecd)*t3+q2ecd2*(-3.d0*t4+
     +    0.5d0*s3+qperp2*t5)
      h(118) = 1.5d0*qecd2*(t7+b2)-hqecd*(t6+b1)-q3ecd2*(b3+t8)
      h(214) = 1.5d0*qecd2*(t7+a2)-hqecd*(t6+a1)-q3ecd2*(a3+t8)
      h(246) = hxxyy-hecd2*(qperp2*t4+t10+a5+b5)+hecd*(t9+a4+b4+c1-ecd*
     +    c2)+q2ecd2*(t11+0.5d0*s3+a6+b6+c3)
      h(168) = hqecd2*(t7-s7)
      h( 88) = 1.5d0*qecd2*t7-hqecd2*s7-q3ecd2*t8
      temp = hecd2*(0.5d0*t3-t10)+q2ecd2*(t11-0.5d0*t4)
      h(120) = temp-hecd2*b5+q2ecd2*b6
      h(216) = temp-hecd2*a5+q2ecd2*a6
      h(248) = qecd2*(1.5d0*t7-t13-a9-b9-c5)-hqecd2*(s7-a2-b2)
      h(176) = hxxyy+hecd2*(s10-t10)
      h( 96) = hxxyy-hecd2*(qperp2*t4+t10-s10)+0.5d0*q2ecd*t3+q2ecd2*
     +     t11
      h(128) = qecd2*(1.5d0*t7-t13-b9)-hqecd*(t6+b1)+hqecd2*b2
      h(224) = qecd2*(1.5d0*t7-t13-a9)-hqecd*(t6+a1)+hqecd2*a2
      h(256) = hxxyy+hecd2*(-3.d0*(a5+b5)+t3+s10-c2)+ecd2*(-3.d0*t10+
     +     t14+a10+b10+c6)+hecd*(t9+a4+b4+c1)
      if (sinp) 120,100,120
  100 if (cosp) 1000,120,920
 120  u12 = -sinp
      v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      v45 = v57+v57
      v55 = v44-v47
      do 103 kq1=22,214,48
          g(kq1  ) = v44*h(kq1) + v47*h(kq1+5)
          g(kq1+1) = v54*h(kq1) + v57*h(kq1+5)
103       g(kq1+5) = v74*h(kq1) + v77*h(kq1+5)
      do 101 kq1=24,216,48
          g(kq1  ) = cosp*h(kq1)
101       g(kq1+4) = sinp*h(kq1)
      do 102 kq1=18,210,48
          g(kq1  ) = cosp*h(kq1)
102       g(kq1+1) = sinp*h(kq1)
      g( 80) = h( 80)
      g( 86) = v44*h( 86)+v47*h( 91)
      g( 87) = v54*h( 86)+v57*h( 91)
      g( 91) = v74*h( 86)+v77*h( 91)
      g( 88) = cosp*h( 88)
      g( 92) = sinp*h( 88)
      g( 96) = h( 96)
      g(102) = v45*h(103)
      g(103) = v55*h(103)
      g(107) = -g(102)
      g(104) = u12*h(108)
      g(108) = cosp*h(108)
      g(112) = dzero
      g(128) = h(128)
      g(134) = v45*h(135)
      g(135) = v55*h(135)
      g(139) = -g(134)
      g(136) = u12*h(140)
      g(140) = cosp*h(140)
      g(144) = dzero
      g(150) = v45*h(151)
      g(151) = v55*h(151)
      g(155) = -g(150)
      g(152) = u12*h(156)
      g(156) = cosp*h(156)
      g(160) = dzero
      g(176) = h(176)
      g(182) = v45*h(183)
      g(183) = v55*h(183)
      g(187) = -g(182)
      g(184) = u12*h(188)
      g(188) = cosp*h(188)
      g(192) = dzero
      g(198) = v44*h(198)+v47*h(203)
      g(199) = v54*h(198)+v57*h(203)
      g(203) = v74*h(198)+v77*h(203)
      g(200) = cosp*h(200)
      g(204) = sinp*h(200)
      g(230) = v45*h(231)
      g(231) = v55*h(231)
      g(235) = -g(230)
      g(232) = u12*h(236)
      g(236) = cosp*h(236)
      g(240) = dzero
      g(246) = v44*h(246)+v47*h(251)
      g(247) = v54*h(246)+v57*h(251)
      g(251) = v74*h(246)+v77*h(251)
      g(248) = cosp*h(248)
      g(252) = sinp*h(248)
      g( 38) = v45*h( 39)
      g( 39) = v55*h( 39)
      g( 43) = -g( 38)
      g( 40) = u12*h( 44)
      g( 44) = cosp*h( 44)
      g( 48) = dzero
      g( 54) = v44*h( 54)+v47*h( 59)
      g( 55) = v54*h( 54)+v57*h( 59)
      g( 59) = v74*h( 54)+v77*h( 59)
      g( 56) = cosp*h( 56)
      g( 60) = sinp*h( 56)
      g(  6) = v44*h(  6)+v47*h( 11)
      g(  7) = v54*h(  6)+v57*h( 11)
      g( 11) = v74*h(  6)+v77*h( 11)
      g(  8) = cosp*h(  8)
      g( 12) = sinp*h(  8)
      g( 68) = h( 68)
      g( 82) = cosp*h( 82)
      g( 83) = sinp*h( 82)
      g( 84) = h( 84)
      g( 98) = u12*h( 99)
      g( 99) = cosp*h( 99)
      g(100) = dzero
      g(116) = h(116)
      g(130) = u12*h(131)
      g(131) = cosp*h(131)
      g(132) = dzero
      g(146) = u12*h(147)
      g(147) = cosp*h(147)
      g(148) = dzero
      g(164) = h(164)
      g(178) = u12*h(179)
      g(179) = cosp*h(179)
      g(180) = dzero
      g(194) = cosp*h(194)
      g(195) = sinp*h(194)
      g(226) = u12*h(227)
      g(227) = cosp*h(227)
      g(228) = dzero
      g(242) = cosp*h(242)
      g(243) = sinp*h(242)
      g( 34) = u12*h( 35)
      g( 35) = cosp*h( 35)
      g( 36) = dzero
      g( 50) = cosp*h( 50)
      g( 51) = sinp*h( 50)
      g(  2) = cosp*h(  2)
      g(  3) = sinp*h(  2)
      g( 65) = h( 65)
      g( 81) = h( 81)
      g( 97) = dzero
      g(113) = h(113)
      g(129) = dzero
      g(145) = dzero
      g(161) = h(161)
      g(177) = dzero
      g(225) = dzero
      g( 33) = dzero
      h( 80) = cosp*g( 80)
      h( 96) = cosp*g( 96)
      h(112) =           u12*g(176)
      h(128) = cosp*g(128)
      h(144) = sinp*g( 80)
      h(160) = sinp*g( 96)
      h(176) =           cosp*g(176)
      h(192) = sinp*g(128)
      do 121 kq1=70,118,16
          h(kq1   ) = cosp*g(kq1   ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1   ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+ 1) + u12*g(kq1+65)
          h(kq1+65) = sinp*g(kq1+ 1) + cosp*g(kq1+65)
          h(kq1+ 2) = cosp*g(kq1+ 2) + u12*g(kq1+66)
          h(kq1+66) = sinp*g(kq1+ 2) + cosp*g(kq1+66)
          h(kq1+ 5) = cosp*g(kq1+ 5) + u12*g(kq1+69)
          h(kq1+69) = sinp*g(kq1+ 5) + cosp*g(kq1+69)
          h(kq1+ 6) = cosp*g(kq1+ 6) + u12*g(kq1+70)
          h(kq1+70) = sinp*g(kq1+ 6) + cosp*g(kq1+70)
121   continue
      h( 68) = cosp*g( 68)
      h( 84) = cosp*g( 84)
      h(100) =           u12*g(164)
      h(116) = cosp*g(116)
      h(132) = sinp*g( 68)
      h(148) = sinp*g( 84)
      h(164) =           cosp*g(164)
      h(180) = sinp*g(116)
      do 122 kq1=66,114,16
          h(kq1   ) = cosp*g(kq1  ) + u12*g(kq1+64)
          h(kq1+64) = sinp*g(kq1  ) + cosp*g(kq1+64)
          h(kq1+ 1) = cosp*g(kq1+1) + u12*g(kq1+65)
122       h(kq1+65) = sinp*g(kq1+1) + cosp*g(kq1+65)
      do 1221 kq1=2,50,16
          h(kq1   ) = g(kq1   )
          h(kq1+ 1) = g(kq1+ 1)
          h(kq1+ 4) = g(kq1+ 4)
          h(kq1+ 5) = g(kq1+ 5)
          h(kq1+ 6) = g(kq1+ 6)
1221      h(kq1+ 9) = g(kq1+ 9)
      do 1222 kq1=12,60,16
          h(kq1    ) = g(kq1    )
          h(kq1+182) = g(kq1+182)
          h(kq1+183) = g(kq1+183)
          h(kq1+186) = g(kq1+186)
          h(kq1+187) = g(kq1+187)
1222      h(kq1+188) = g(kq1+188)
      do 1223 kq1=203,251,16
          h(kq1  ) = g(kq1  )
1223      h(kq1+1) = g(kq1+1)
      h( 65) = cosp*g( 65)
      h( 81) = cosp*g( 81)
      h( 97) =           u12*g(161)
      h(113) = cosp*g(113)
      h(129) = sinp*g( 65)
      h(145) = sinp*g( 81)
      h(161) =           cosp*g(161)
      h(177) = sinp*g(113)
      h( 48) = g( 48)
      h( 36) = g( 36)
      h(228) = g(228)
      h(240) = g(240)
      h(225) = g(225)
      h( 33) = g( 33)
      do 123 kq1=22,214,64
          g(kq1   ) = cosp*h(kq1   ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1   ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+ 1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+ 1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+ 2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+ 2) + cosp*h(kq1+18)
          g(kq1+ 5) = cosp*h(kq1+ 5) + u12* h(kq1+21)
          g(kq1+21) = sinp*h(kq1+ 5) + cosp*h(kq1+21)
          g(kq1+ 6) = cosp*h(kq1+ 6) + u12* h(kq1+22)
          g(kq1+22) = sinp*h(kq1+ 6) + cosp*h(kq1+22)
          g(kq1+10) = cosp*h(kq1+10) + u12* h(kq1+26)
123       g(kq1+26) = sinp*h(kq1+10) + cosp* h(kq1+26)
      do 124 kq1=17,209,64
          g(kq1   ) = cosp*h(kq1  ) + u12* h(kq1+16)
          g(kq1+16) = sinp*h(kq1  ) + cosp*h(kq1+16)
          g(kq1+ 1) = cosp*h(kq1+1) + u12* h(kq1+17)
          g(kq1+17) = sinp*h(kq1+1) + cosp*h(kq1+17)
          g(kq1+ 2) = cosp*h(kq1+2) + u12* h(kq1+18)
          g(kq1+18) = sinp*h(kq1+2) + cosp*h(kq1+18)
          g(kq1+ 3) = cosp*h(kq1+3) + u12* h(kq1+19)
124       g(kq1+19) = sinp*h(kq1+3) + cosp* h(kq1+19)
      do 125 kq1=49,177,64
          kkq1=kq1
          do 126 kkkq1=1,32
              g(kkq1)=h(kkq1)
              kkq1=kkq1+1
126       continue
125   continue
      do 127 kq1=1,16
          g(kq1)=h(kq1)
127       g(kq1+240)=h(kq1+240)
      goto 2000
920   do 921 kkq1=1,10
          kq1=kq1off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) = h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) = h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
921       g(kq1+15) = h(kq1+15)
      do 922 kkq1=1,6
          kq1=kq2off(kkq1)
          g(kq1   ) = dzero
          g(kq1+ 1) = dzero
          g(kq1+ 2) = h(kq1+ 2)
          g(kq1+ 3) = dzero
          g(kq1+ 5) = dzero
          g(kq1+ 6) = h(kq1+ 6)
          g(kq1+ 7) = dzero
          g(kq1+10) = dzero
          g(kq1+11) = h(kq1+11)
922       g(kq1+15) = dzero
      go to 2000
1000  do 1001 kkq1=1,6
          kq1=kq3off(kkq1)
          g(kq1   ) = h(kq1   )
          g(kq1+ 1) =-h(kq1+ 1)
          g(kq1+ 2) = dzero
          g(kq1+ 3) = h(kq1+ 3)
          g(kq1+ 5) = h(kq1+ 5)
          g(kq1+ 6) = dzero
          g(kq1+ 7) =-h(kq1+ 7)
          g(kq1+10) = h(kq1+10)
          g(kq1+11) = dzero
1001      g(kq1+15) = h(kq1+15)
      do 1002 kkq1=1,4
          kq1=kq4off(kkq1)
          g(kq1   ) = -h(kq1   )
          g(kq1+ 1) =  h(kq1+ 1)
          g(kq1+ 2) =  dzero
          g(kq1+ 3) = -h(kq1+ 3)
          g(kq1+ 5) = -h(kq1+ 5)
          g(kq1+ 6) =  dzero
          g(kq1+ 7) =  h(kq1+ 7)
          g(kq1+10) = -h(kq1+10)
          g(kq1+11) =  dzero
1002      g(kq1+15) = -h(kq1+15)
      do 1003 kkq1=1,6
          kq1=kq5off(kkq1)
          g(kq1   ) =  dzero
          g(kq1+ 1) =  dzero
          g(kq1+ 2) =  h(kq1+ 2)
          g(kq1+ 3) =  dzero
          g(kq1+ 5) =  dzero
          g(kq1+ 6) = -h(kq1+ 6)
          g(kq1+ 7) =  dzero
          g(kq1+10) =  dzero
          g(kq1+11) =  h(kq1+11)
1003      g(kq1+15) =  dzero
          g(99)=-g(99)
          g(108)=-g(108)
          g(147)=-g(147)
          g(156)=-g(156)
          g(103)=-g(103)
          g(151)=-g(151)
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      do 2001 kq1=2,242,16
          g(kq1+ 3) = g(kq1   )
          g(kq1+ 7) = g(kq1+ 1)
          g(kq1+11) = g(kq1+ 2)
          g(kq1+ 8) = g(kq1+ 5)
          g(kq1+12) = g(kq1+ 6)
2001      g(kq1+13) = g(kq1+10)
      if (rcdsq) 1200,1200,1300
1300  do 1301 kq1=1,4
          kkq1=kq1
          do 1302 jq1=1,16
              g(kkq1+4) = r13*g(kkq1) + g(kkq1+4)
              g(kkq1+12)= r33*g(kkq1) + g(kkq1+12)
1302          kkq1=kkq1+16
1301  continue
c ***
      do 1303 kq1=1,253,4
          g(kq1+1) = r14*g(kq1) + g(kq1+1)
1303      g(kq1+3) = r34*g(kq1) + g(kq1+3)
c ***
c1200  do 1201 kq1=1,256
c1201      gout(kq1) = dq11*g(kq1) + gout(kq1)
 1200  continue
       call daxpy(256,dq11,g,1,gout,1)
      dq01x=dq01-dq11
      do 1202 kq1=2,242,16
          gout(kq1  ) = dq01x*g(kq1  ) + gout(kq1  )
          gout(kq1+1) = dq01x*g(kq1+1) + gout(kq1+1)
1202      gout(kq1+2) = dq01x*g(kq1+2) + gout(kq1+2)
      dq10x=dq10-dq11
      dq00x=dq00-dq11
      do 1203 kq1=1,241,16
          gout(kq1   ) = dq00x*g(kq1   ) + gout(kq1   )
          gout(kq1+ 4) = dq10x*g(kq1+ 4) + gout(kq1+ 4)
          gout(kq1+ 8) = dq10x*g(kq1+ 8) + gout(kq1+ 8)
1203      gout(kq1+12) = dq10x*g(kq1+12) + gout(kq1+12)
940   continue
c
c     --------------------------
c     --------------------------
c
c     rotates up to 256 integrals to space fixed axes
c     incoming and outgoing integrals in common gout
c     indices in order 0000,0001,0002,...0010,0011,...0100,0101,...etc.
c     p11,...are direction cosines of space fixed axes wrt axes at p
c     q11,...are direction cosines of space fixed axes wrt axes at q
c     applies to case 1111
c
c
      i1 = 64
      i2 = 128
      i3 = 192
      do 1 jkl = 1,64
      i1 = i1+1
      i2 = i2+1
      i3 = i3+1
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    1 continue
      ind = -48
      do 2 i = 1,4
      ind = ind+48
      do 2 loopkl = 1,16
      ind = ind+1
      i1 = 16+ind
      i2 = 32+ind
      i3 = 48+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    2 continue
      ind = -12
      do 3 loopij = 1,16
      ind = ind+12
      do 3 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    3 continue
      ind = -3
      do 4 ijk = 1,64
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
    4 continue
c
      return
      end
      subroutine ver_integs(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integs.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
