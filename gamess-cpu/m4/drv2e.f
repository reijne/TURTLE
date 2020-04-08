c     deck=drv2e
c ******************************************************
c ******************************************************
c             =   drv2e   =
c ******************************************************
c ******************************************************
      subroutine dabab(ii,jj,kk,ll,q4,zscftp,da,db,abdens)
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
      common/restrl/ociopt,ocifor,omp2
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
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
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
      dimension da(*),db(*),abdens(*)
      data zuhf,zgrhf /'uhf','grhf'/
      data pt5,four /0.5d0,4.0d0/
c
      ouhf = zscftp.eq.zuhf
      ogrhf = zscftp.eq.zgrhf
      ohf_coul = (.not.CD_active()).or.CD_HF_coulomb_deriv()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
      hf_wght  = CD_HF_exchange_weight()
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
c
      if (.not.ogrhf) then
c
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
                     dfac=0.0d0
                     if(ohf_coul)then
                        dfac = da(mij)*da(mkl)*four 
                     endif
                     if(ohf_exch)then
                        dfac = dfac - hf_wght*(da(mik)*da(mjl) 
     &                                       + da(mil)*da(mjk))
                        if(ouhf)dfac = dfac - hf_wght*(db(mik)*db(mjl)
     &                                               + db(mil)*db(mjk))
                     endif
                     if (i1.eq.i2) dfac = dfac*pt5
                     if (i3.eq.i4) dfac = dfac*pt5
                     dfac = dfac*q4
                     if (omp2 .or. mp3) then
                        abdens(nn) = abdens(nn) + dfac
                     else
                        abdens(nn) = dfac
                     end if
                     nl = nl + inc2
 30               continue
                  nk = nk + inc3
 40            continue
               nj = nj + inc4
 50         continue
            ni = ni + inc5
 60      continue
c
      else
c
c     general case
c
         do 120 i = mini , maxi
            nj = ni
            i1 = loci + i
            ii1 = iky(i1)
            do 110 j = minj , maxj
               nk = nj
               j1 = locj + j
               jj1 = iky(j1)
               mij = ii1 + j1
               if (j1.gt.i1) mij = jj1 + i1
               oijeq = i1.eq.j1
               do 100 k = mink , maxk
                  nl = nk
                  k1 = lock + k
                  kk1 = iky(k1)
                  mik = ii1 + k1
                  if (k1.gt.i1) mik = kk1 + i1
                  mjk = jj1 + k1
                  if (k1.gt.j1) mjk = kk1 + j1
                  do 90 l = minl , maxl
                     nn = nl
                     l1 = locl + l
                     ll1 = iky(l1)
                     mkl = kk1 + l1
                     if (l1.gt.k1) mkl = ll1 + k1
                     mjl = jj1 + l1
                     if (l1.gt.j1) mjl = ll1 + j1
                     mil = ii1 + l1
                     if (l1.gt.i1) mil = ll1 + i1
                     okleq = k1.eq.l1
                     dfac = 0.0d0
                     ioff = 0
                     do 80 is = 1 , njk
                        isi = (is-1)*11
                        joff = 0
                        do 70 js = 1 , njk
                           dfac = dfac + 4.0d0*erga(isi+js)
     +                            *(da(ioff+mij)*da(joff+mkl)
     +                            +da(ioff+mkl)*da(joff+mij))
     +                            + 2.0d0*ergb(isi+js)
     +                            *(da(ioff+mik)*da(joff+mjl)
     +                            +da(ioff+mjl)*da(joff+mik)
     +                            +da(ioff+mil)*da(joff+mjk)
     +                            +da(ioff+mjk)*da(joff+mil))
                           joff = joff + nx
 70                     continue
                        ioff = ioff + nx
 80                  continue
                     if (oijeq) dfac = dfac*pt5
                     if (okleq) dfac = dfac*pt5
                     abdens(nn) = dfac*q4
                     nl = nl + inc2
 90               continue
                  nk = nk + inc3
 100           continue
               nj = nj + inc4
 110        continue
            ni = ni + inc5
 120     continue
c
      end if
c
      return
      end
      subroutine dabg(ii,jj,kk,ll,l1,norb,q4,da,v,nconf,
     + onocor,onopen,abdens)
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      dimension abdens(*),da(*),v(l1,*),nconf(*)
      data dzero,pt5,eight/0.0d0,0.5d0,8.0d0/
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
      do 100 i = mini , maxi
         nj = ni
         do 90 j = minj , maxj
            nk = nj
            do 80 k = mink , maxk
               nl = nk
               do 70 l = minl , maxl
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
 20               mij = iky(i1) + i2
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
                  dfac = dzero
                  if (.not.(onocor)) then
                     dfac = dfac + alpha(1)*da(mij)*da(mkl)
     +                      + pt5*beta(1)
     +                      *(da(mik)*da(mjl)+da(mil)*da(mjk))
                     if (onopen) go to 60
                     nco1 = nco + 1
                     do 30 io = nco1 , norb
                        iojo = iky(nconf(io)) + 1
                        dfac = dfac + alpha(iojo)
     +                         *(da(mij)*v(i3,io)*v(i4,io)+da(mkl)
     +                         *v(i1,io)*v(i2,io)) + pt5*beta(iojo)
     +                         *(v(i2,io)
     +                         *(da(mik)*v(i4,io)+da(mil)*v(i3,io))
     +                         +v(i1,io)
     +                         *(da(mjl)*v(i3,io)+da(mjk)*v(i4,io)))
 30                  continue
                  end if
                  if (.not.(onopen)) then
                     nco1 = nco + 1
                     do 50 io = nco1 , norb
                        do 40 jo = nco1 , norb
                           iof = nconf(io)
                           jof = nconf(jo)
                           iojo = iky(iof) + jof
                           if (jof.gt.iof) iojo = iky(jof) + iof
                           dfac = dfac + alpha(iojo)*v(i1,io)*v(i2,io)
     +                            *v(i3,jo)*v(i4,jo) + pt5*beta(iojo)
     +                            *v(i1,io)*v(i2,jo)
     +                            *(v(i3,io)*v(i4,jo)+v(i4,io)*v(i3,jo))
 40                     continue
 50                  continue
                  end if
 60               if (i1.eq.i2) dfac = dfac*pt5
                  if (i3.eq.i4) dfac = dfac*pt5
                  dfac = dfac*eight*q4
                  abdens(nn) = dfac
                  nl = nl + inc2
 70            continue
               nk = nk + inc3
 80         continue
            nj = nj + inc4
 90      continue
         ni = ni + inc5
 100  continue
      return
      end
      subroutine dabmc(ii,jj,kk,ll,q4,db,dc,abdens)
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
c
      dimension dc(*),db(*),abdens(*)
      data half,four /0.5d0,4.0d0/
c
      if (ncore.eq.0) return
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
c
c    additions for ci core
c
      ni = 1
      do 50 i = mini , maxi
         nj = ni
         i1 = loci + i
         ii1 = iky(i1)
         do 40 j = minj , maxj
            nk = nj
            j1 = locj + j
            jj1 = iky(j1)
            mij = ii1 + j1
            if (j1.gt.i1) mij = jj1 + i1
            oijeq = i1.eq.j1
            do 30 k = mink , maxk
               nl = nk
               k1 = lock + k
               kk1 = iky(k1)
               mik = ii1 + k1
               if (k1.gt.i1) mik = kk1 + i1
               mjk = jj1 + k1
               if (k1.gt.j1) mjk = kk1 + j1
               do 20 l = minl , maxl
                  nn = nl
                  l1 = locl + l
                  ll1 = iky(l1)
                  mkl = kk1 + l1
                  if (l1.gt.k1) mkl = ll1 + k1
                  mjl = jj1 + l1
                  if (l1.gt.j1) mjl = ll1 + j1
                  mil = ii1 + l1
                  if (l1.gt.i1) mil = ll1 + i1
                  okleq = k1.eq.l1
c
c     db ----- core density
c     dc ----- ci density
c
                  dfaca = db(mij)*db(mkl)*four - db(mik)*db(mjl)
     +                    - db(mjk)*db(mil)
                  dfacb = (db(mij)*dc(mkl)+db(mkl)*dc(mij))
     +                    *four - db(mil)*dc(mjk) - db(mjk)*dc(mil)
     +                    - db(mik)*dc(mjl) - db(mjl)*dc(mik)
                  dfac = dfaca*four + dfacb + dfacb
                  if (oijeq) dfac = dfac*half
                  if (okleq) dfac = dfac*half
                  abdens(nn) = abdens(nn) + dfac*q4
                  nl = nl + inc2
 20            continue
               nk = nk + inc3
 30         continue
            nj = nj + inc4
 40      continue
         ni = ni + inc5
 50   continue
      return
      end
      subroutine dabout(q,odebug,iw)
c -----------------------------------
c    write out results from fokabd
c ------------------------------------
      implicit real*8  (a-h,p-z),integer   (i-n),logical    (o)
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
      parameter (mxp2 = mxprms*mxprms)
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      dimension q(*)
      data half/0.5d0/
c
c     writes out accumulated derivative fock
c     onto fockfile.
c
      if (odebug) write (iw,6010)
      ioff = ifok
      ifocbl = ibdout
      do 50 n = 1 , nfok
         call rdedx(q(ida+1),nx,ifocbl,iflout)
         ij = 0
         do 30 k = 1 , num
            do 20 l = 1 , k
               ij = ij + 1
               if (k.ne.l) then
                  q(ioff+ij) = half*q(ioff+ij)
               end if
 20         continue
 30      continue
         do 40 k = 1 , nx
            q(ida+k) = q(ida+k) + q(ioff+k)
 40      continue
         call wrt3(q(ida+1),nx,ifocbl,iflout)
         if (odebug) call prtris(q(ida+1),num,iw)
         ioff = ioff + nx
         ifocbl = ifocbl + lenb
 50   continue
c
      return
 6010 format (///1x,'output from dabout'//)
      end
      subroutine ddebut(zscftp,da,db,nconf,dd,fock,
     + l1,l2,l3,ista,jsta,ksta,lsta,q)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/restrl/ociopt(2),omp2
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
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
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension da(*),db(*),nconf(*),fock(*),dd(*), q(*)
c
      real*8 de
      common/grad2/de(3,maxat)
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
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
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
      data zvb,zmcscf/'vb','mcscf'/
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
c
      if (onocnt) write (iwr,6010)
c
      if (zscftp.eq.zgvb) then
c
c     ----- set up core density matrix, read in  eigenvectors for
c            gvb -----
c
         call rdedx(db,l3,ibl3qa,idaf)
         call tdown(db,ilifq,db,ilifq,l1)
         call dencor(da,db,l1)
c
c     ----- set up mo to fock operator pointers in nconf -----
c
         if (nco.ne.0) call setsto(nco,1,nconf)
         ic = 0
         if (nseto.gt.0) then
            nbase = ncores
            do 30 i = 1 , nseto
               nop = no(i)
               do 20 j = 1 , nop
                  nconf(ic+nco+j) = nbase + 1
 20            continue
               ic = ic + nop
               nbase = nbase + 1
 30         continue
         end if
         if (npair.gt.0) then
            np2 = npair + npair
            do 40 i = 1 , np2
               nconf(i+nco+ic) = ncores + nseto + i
 40         continue
         end if
         norb = nco + npair + npair + ic
         onocor = .false.
         onopen = .false.
         onocor = nco.eq.0
         onopen = nseto.eq.0 .and. npair.eq.0
c
      else if (zscftp.eq.zuhf) then
c
c     ----- read in density matrices (alpha+beta) in uhf -----
c
         call rdedx(da,l2,ibl3pa,idaf)
         call rdedx(db,l2,ibl3pb,idaf)
         do 50 i = 1 , l2
            duma = da(i)
            dumb = db(i)
            da(i) = duma + dumb
            db(i) = duma - dumb
 50      continue
c
      else if (zscftp.eq.zrhf) then
c
c     ----- read in density matrix  in rhf  -----
c
         call rdedx(da,l2,ibl3pa,idaf)
c
      else if (zscftp.eq.zgrhf) then
c
c     general scf
c
         call rdedx(fock,l3,ibl3qa,idaf)
         call tdown(fock,ilifq,fock,ilifq,l1)
         m = 0
         do 90 is = 1 , njk
            call vclr(da(m+1),1,l2)
            nsi = nbshel(is)
            ils = ilfshl(is)
            do 80 ni = 1 , nsi
               nc = iactiv(ils+ni)
               ncol = (nc-1)*l1
               ij = 0
               do 70 i = 1 , l1
                  dum = fock(i+ncol)
                  do 60 j = 1 , i
                     ij = ij + 1
                     da(ij+m) = da(ij+m) + dum*fock(j+ncol)
 60               continue
 70            continue
 80         continue
            m = m + l2
 90      continue
c
      else if ((zscftp.eq.zmcscf).or.(zscftp.eq.zvb)) then
c
c     mcscf/multi
c
         m = 0
         call secget(isecmo,m,iblok)
         call rdedx(da,l1*ncoorb,iblok+mvadd,idaf)
         m = 0
         call secget(isecdd,m,iblok)
         if(odebug(30)) write (iwr,6040) iblok
         call rdedx(dd,l2,iblok,idaf)
         if(odebug(30)) write (iwr,6050) ifil2d,iblk2d
         call rdedx(fock,nd2mo,iblk2d,ifil2d)
         if (ncore.gt.0) then
            call vclr(db,1,l2)
            ij = 0
            do 120 i = 1 , l1
               do 110 j = 1 , i
                  ij = ij + 1
                  do 100 k = 1 , ncore
                     kk = (k-1)*l1
                     db(ij) = db(ij) + da(i+kk)*da(j+kk)
 100              continue
 110           continue
 120        continue
            write (iwr,6020) ncore
         end if
      else
         call caserr('invalid scftype detected in gradient code')
      end if
c
c     ----- set starting parameters -----
c
      outd= nprint.eq. - 4
      if (ofokab .or. ompir) then
c
c      read in ndenin density matrices from some external file
c      matrices are in mo basis
c
         ioff = 1
         call search(ibden,iflden)
         do 130 loop = 1 , ndenin
            call reads(db(ioff),l2,iflden)
            if (odebug(21)) call prtris(db(ioff),l1,iwr)
            ioff = ioff + l2
 130     continue

         ioffd = igmem_alloc(l2)
         iofft = igmem_alloc(l1)
         ioffv = igmem_alloc(l3)
         m = 0
         call rdedx(q(ioffv),l3,ibl3qa,idaf)
         call tdown(q(ioffv),ilifq,q(ioffv),ilifq,l1)
         ioff = 1
         do 140 loop = 1 , ndenin
            call dcopy(l2,db(ioff),1,q(ioffd),1)
            call demoao(q(ioffd),db(ioff),q(ioffv),q(iofft),l1,
     +  ncoorb,l1)
            ioff = ioff + l2
 140     continue

         call gmem_free(ioffv)
         call gmem_free(iofft)
         call gmem_free(ioffd)

      end if
      if (omp2 .or. ofokab .or. ofock) then
         icutd = icut + 1
      else
         icutd = icut
      end if
      icutd = max(icutd,8)
      if (omp2 .or. ofokab .or. ofock) then
         itold = itol + 2
      else
         itold = itol + 1
      end if
      itold = max(itold,16)
      dcut = done/(ten**icutd)
      tol1 = e*(itold+2)
      tol2 = e*itold
      tol3 = done/ten**(itold+2)
      tol4 = done/ten**itold
      onorm = normf.ne.1 .or. normp.ne.1
c
c     ----- read in 1e-gradient -----
c
      call rdgrd(de,irest,ista,jsta,ksta,lsta)
c
      if (ofock .or. ofokab) then
         if(odebug(30)) write (iwr,6030)
         call vclr(fock,1,nfok*l2)
      end if
      return
 6010 format (/1x,22('=')/1x,'gradient of the energy'/1x,22('='))
 6020 format (/1x,'number of frozen and core orbitals',i5/)
 6030 format (/1x,'zeroing storage for fock matrices')
 6040 format (/1x,'reading density matrix from block',i5)
 6050 format (/1x,'tpdm from file, block', 2i5)
      end
      subroutine delim(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd,
     * abmax)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      common/restrl/ociopt(2),omptwo,ohf(7),omp2w
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      do 20 i = 1 , lendd
         oform(i) = .false.
 20   continue
      abmax = 1.0d0
      if (.not.(ofock .or. ofokab .or. ompir .or. omp2w)) then
         abmax = 0.0d0
         nn = 0
         do 40 i = 1 , ij
            do 30 k = 1 , kl
               nn = nn + 1
               n = ijgt(i) + klgt(k)
               ab = dabs(qa(n))
               if (ab.lt.dcut) oform(nn) = .true.
               if (ab.gt.abmax) abmax = ab
 30         continue
 40      continue
      end if
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
      subroutine delim2(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      call dcopy(lendd,qa,1,qb,1)
      nn = 0
      ijkl = 0
      do 30 i = 1 , ij
         do 20 k = 1 , kl
            nn = nn + 1
            if (.not.(oform(nn))) then
               n = ijgt(i) + klgt(k)
               ijkl = ijkl + 1
               qa(ijkl) = qb(n)
            end if
 20      continue
 30   continue
      return
      end
      subroutine dencor(da,v,l1)
c     ----- calculate the density matrix for the core orbitals -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      dimension da(*),v(l1,*)
c
      if (nco.ne.0) then
         l2 = l1*(l1+1)/2
         call vclr(da,1,l2)
         j = 1
         do 30 i = 1 , l1
            do 20 k = 1 , nco
               fact = f(1)*v(i,k)
               if (fact.ne.0.0d0) then
                  call daxpy(i,fact,v(1,k),1,da(j),1)
               end if
 20         continue
            j = j + i
 30      continue
      end if
      return
      end
      subroutine denfac(dkl,csk,cpk,cdk,cfk,cgk,
     +                      csl,cpl,cdl,cfl,cgl,
     +                     mink,maxk,minl,maxl,okandl,double)
      implicit real*8  (a-h,o-z)
       logical okandl,double
       dimension dkl(*)
      if (.not.double) then
         n = 0
         mmax = maxl
         do 130 k = mink , maxk
            if (okandl) mmax = k
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
 60         do 120 l = minl , mmax
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
 120        continue
 130     continue
      else
         n = 0
         mmax = maxl
         do 250 k = mink , maxk
            if (okandl) mmax = k
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
 180        do 240 l = minl , mmax
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
 240        continue
 250     continue
      end if
      return
      end
      subroutine derwrt(qq)
      implicit real*8  (a-h,o-z)
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
      dimension qq(*)
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
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      logical xskip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dpqd(27,maxat),
     +    nunpr,xskip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/ipp(5,3120)
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      dimension mpert(12)
c
      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1
c
      cut = 1.0d0 / (10.0d0**(icut-1))
      npp = 0
      do 30 i = 1 , npass + 1
         do 20 j = 1 , 3
            npp = npp + 1
            mpert(npp) = (natomd(i)-1)*3 + j
 20      continue
 30   continue
c
      np3 = (npass+1)*3
      nn = 0
      jmax = maxj
      do 80 i = mini , maxi
         i1 = loci + i
         if (oiandj) jmax = i
         do 70 j = minj , jmax
            i2 = locj + j
            lmax = maxl
            if (osame) maxk = i
            do 60 k = mink , maxk
               i3 = lock + k
               if (okandl) lmax = k
               if (osame .and. i.eq.k) lmax = j
               do 50 l = minl , lmax
                  i4 = locl + l
                  nn = nn + 1
cc
                  if (osame) then
                     ijkl1 = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                       *inc3 + l - minl + 1
                     ijkl2 = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                       *inc3 + j - minj + 1
                  end if
                  fac = 1.0d0
                  if (i1.eq.i2) fac = fac + fac
                  if (i3.eq.i4) fac = fac + fac
cc
                  ioa = 0
                  do 40 npp = 1 , np3
                     derval = qq(nn+ioa)*fac
                     if (osame) derval = qq(ijkl1+ioa) + qq(ijkl2+ioa)
                     if (dabs(derval).ge.cut) then
                        if (iperm(mpert(npp)).ne.0) then
                           icnt = icnt + 1
                           nap = mpert(npp)
                           j1 = i1
                           j2 = i2
                           j3 = i3
                           j4 = i4
                           if (xskip(nap)) then
                              jop = iperm(nap)
                              j1 = ibasis(jop,i1)
                              j2 = ibasis(jop,i2)
                              j3 = ibasis(jop,i3)
                              j4 = ibasis(jop,i4)
                              nupert = mptr(jop,nap)
                              sigh = 
     +                           dsign(1.0d0,dble(j1*j2)*
     +                                       dble(j3*j4*nupert))
                              derval = derval*sigh
                              jj1 = iabs(j1)
                              jj2 = iabs(j2)
                              jj3 = iabs(j3)
                              jj4 = iabs(j4)
                              j1 = max(jj1,jj2)
                              j2 = min(jj1,jj2)
                              j3 = max(jj3,jj4)
                              j4 = min(jj3,jj4)
                              nap = iabs(nupert)
                           end if
                           val(icnt) = derval
                           if (j1.ge.j3) then
                              ipp(1,icnt) = j1
                              ipp(2,icnt) = j2
                              ipp(3,icnt) = j3
                              ipp(4,icnt) = j4
                           else
                              ipp(1,icnt) = j3
                              ipp(2,icnt) = j4
                              ipp(3,icnt) = j1
                              ipp(4,icnt) = j2
                           end if
                           if (j1.eq.j3 .and. j4.gt.j2) then
                              ipp(4,icnt) = j2
                              ipp(2,icnt) = j4
                           end if
                           ipp(5,icnt) = nap
                           if (icnt.eq.idblk) then
                              call pack(val(idblk1),lab816,ipp,idblk*5)
                              mderb = idblk + iidblk + 2
                              call wrt3s(dbufs,mderb,mpstrm(1))
                              icnt = 0
                           end if
                        end if
                     end if
                     ioa = ioa + lendd
 40               continue
 50            continue
 60         continue
 70      continue
 80   continue
      return
      end

      subroutine dfinal(q,index,ii)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      real*8 de
      common/grad2/de(3,maxat)
c
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
c
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      character*8 zopti,zfrozen
      common/runopt/zopti(maxat),zfrozen
      dimension f(3,maxat),ydnam(3),q(*)
      dimension ibuffi(8)
      character*10 charwall
      equivalence (egrad(1),f(1,1))
      data ydnam /'e/x','e/y','e/z'/
c
      if (index.eq.1) then
c
c...   this might be elimitated if the gradients are properly bypassed
c
         do n=1,nat
            if (zopti(n).eq.'no') then
                de(1,n) = 0.0d0
                de(2,n) = 0.0d0
                de(3,n) = 0.0d0
            end if
         end do
      end if
c
      if (index.ne.1) then
c
c     ----- get restart data -----
c
         irest = 6
         istnu = 1 + ii
         jstnu = 1
         kstnu = 1
         lstnu = 1
c
c     ----- save gradient + restart data -----
c
      call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
c
c     ----- check cpu time -----
c
         if (istnu.gt.nshell) return
         call texit(0,irest)
         if (tim.lt.timlim) return
         write (iwr,6010) tim , charwall(), istnu , jstnu , kstnu , lstnu
c
         if (ofokab) call dabout(q,odebug(21),iwr)
         if (ofock) call hsymd(q,iwr)
         write (iwr,6030)
         return
      else if (onocnt) then
c
c     ----- transfer gradient into -egrad- -----
c
         ncoord = 3*nat
         call dcopy(ncoord,de,1,egrad,1)
c
c     ----- save 1e+2e-gradient -----
c
         irest = 0
         ist = 1
         jst = 1
         kst = 1
         lst = 1
         call wrtgrd(egrad,irest,ist,jst,kst,lst)
c
c --- punch out the gradients
c
         call blkgra(nat,egrad)
c
         if (oprint(31)) then
c
c  print out non-bq centres only
c
         maxpr = 0
         do while (maxpr .lt. nat)
            omore = .true.
            nbuff = 0
            do while (omore)
               if (nbuff .eq. 8 .or. maxpr .ge. nat) then
                  omore = .false.
               else
                  maxpr = maxpr + 1
                  if (zaname(maxpr)(1:2) .ne. 'bq')then
                     nbuff = nbuff + 1
                     ibuffi(nbuff) = maxpr
                  endif
               endif
            enddo
            if (nbuff .gt. 0) then
               write (iwr,6020)
               write (iwr,6050) (ibuffi(i),i=1,nbuff)
               write (iwr,6020)
               do n = 1 , 3
                  write (iwr,6060) ydnam(n),(f(n,ibuffi(i)),i=1,nbuff)
               enddo
            endif
         enddo
         else

         maxpr = 0
 20      minpr = maxpr + 1
         maxpr = maxpr + 8
         if (maxpr.gt.nat) maxpr = nat
         write (iwr,6020)
         write (iwr,6050) (i,i=minpr,maxpr)
         write (iwr,6020)
         do 30 n = 1 , 3
            write (iwr,6060) ydnam(n) , (f(n,i),i=minpr,maxpr)
 30      continue
         if (maxpr.lt.nat) go to 20
c
         endif

         if (ofokab) then
            call dabout(q,odebug(21),iwr)
            iochf(1) = iochf(1) + nat*ndenin*3*lensec(nx)
         end if
         if (ofock) then
            call hsymd(q,iwr)
            call clredx
         end if
         if (ompir) then
            call secget(isect(57),57,iblok)
            do 60 i = 1 , 3
               do 50 j = 1 , 3
                  do 40 k = 1 , nat
                     dipd(i,j,k) = dipd(i,j,k) + dipi(i,j,k)
 40               continue
 50            continue
 60         continue
            call wrt3(dipd,lds(isect(57)),iblok,idaf)
         end if
         call revind
         call clredx
         cpu = cpulft(1)
         if (nprint.ne.-5) write (iwr,6040) cpu ,charwall()
         return
      else
         irest = 6
         istnu = istd
         jstnu = 1
         kstnu = 1
         lstnu = 1
         call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
         return
      end if
 6010 format (/' insufficient time to complete evaluation of 2-electron'
     +        ,' contribution to gradient'//' job dumped at ',f10.2,
     +        ' seconds',a10,' wall'//' next batch ',4i5)
 6020 format (/)
 6030 format (//10x,27('*')/10x,'*** warning ***'/10x,
     +        'this job must be restarted'/10x,27('*')/)
 6040 format (/' end of calculation of the energy gradient at ',f8.2,
     +        ' seconds',a10,' wall'/)
 6050 format (1x,'atom',8(6x,i3,6x))
 6060 format (3x,a3,8f15.7)
      end
      subroutine dform(x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
      implicit real*8  (a-h,o-z)
      logical unroll
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
      parameter (mxp2 = mxprms*mxprms)
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      dimension ix(*),iy(*),iz(*)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     1 yd(ncdim,*),zd(ncdim,*),g(*)
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
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
     +             *z(2,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz)
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
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +             *y(4,my)*z(4,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +             *yd(4,my)*z(4,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +             *y(4,my)*zd(4,mz)
 110  continue
      return
 120  do 130 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +             *y(4,my)*z(4,mz) + xd(5,mx)*y(5,my)*z(5,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +             *yd(4,my)*z(4,mz) + x(5,mx)*yd(5,my)*z(5,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +             *y(4,my)*zd(4,mz) + x(5,mx)*y(5,my)*zd(5,mz)
 130  continue
      return
      end
      subroutine dforma(spij,spkl,noform,
     *                  x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
      implicit real*8  (a-h,o-z)
      logical noform,spij,spkl,unroll
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
      parameter (mxp2 = mxprms*mxprms)
      parameter (ncmax=65)
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      dimension ix(*),iy(*),iz(*)
      dimension noform(*)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     & yd(ncdim,*),zd(ncdim,*),g(*)
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
      subroutine dgenrl(qq,iqq,noform,abmax)
c==========================================
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
      logical trduij,spij,spkl,trdukl
      logical double
      logical sptru,noform
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
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +            aej(ncmax),aek(ncmax),ael(ncmax),
     +            aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +            dd(4*mxp2),ijden(225),ik(225),
     +            ijx(225),ijy(225),ijz(225),klx(225),kly(225),
     +            klz(225),
     +            dij(225),dkl(225),ijgt(225),klgt(225)
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     +   xc00(ncmax),ycp00(ncmax),yc00(ncmax),zcp00(ncmax),
     +   zc00(ncmax),f00(ncmax),
     +   dxij,dyij,dzij,dxkl,dykl,dzkl,
     +   in(12),kn(12),ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/bufb/axak(mxp2),ayak(mxp2),azak(mxp2),axai(mxp2),
     1  ayai(mxp2),azai(mxp2),abv(mxp2),aandbv(mxp2),rhov(mxp2),
     2  xxv(mxp2),c1xv(mxp2),c2xv(mxp2),c3xv(mxp2),c4xv(mxp2),
     3  c1yv(mxp2),c2yv(mxp2),c3yv(mxp2),c4yv(mxp2),
     4  c1zv(mxp2),c2zv(mxp2),c3zv(mxp2),c4zv(mxp2),expev(mxp2)
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      dimension iqq(*),qq(*)
      dimension noform(*)
      data one/1.0d0/
c     data zero,pt5/0.0d0,0.5d0/
      data pi252/34.986836655250d0/
c
      if (ijkld.eq.0) return
      if (ijd.eq.1 .and. kld.eq.1) then
         call ssdss(qq)
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
         call indexa(ijx,ijy,ijz,ijd,mini,maxi,minj,maxj,oiandj,ijn1,
     +              ijn2,1)
         call indexa(klx,kly,klz,kld,mink,maxk,minl,maxl,okandl,kln1,
     +              kln2,0)
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
                  call denfac(dkl,csk,cpk,cdk,cfk,cgk,
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
                                 if (nroots.le.3) call rt123
                                 if (nroots.eq.4) call roots4
                                 if (nroots.eq.5) call roots5
                                 if (nroots.gt.5) call rootss
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
                                    call dxyz(qq(ic1),qq(ic2),qq(ic3),
     +                                 ncmmm)
                                    ioff = 0
                                    if (.not.(oskip(1))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aei,ljt,lkt,llt,lit,
     +                                    ijn2,kln1,kln2,ijn1,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(2))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aej,lit,lkt,llt,ljt,
     +                                    ijn1,kln1,kln2,ijn2,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(3))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aek,lit,ljt,llt,lkt,
     +                                    ijn1,ijn2,kln2,kln1,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(4))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),ael,lit,ljt,lkt,llt,
     +                                    ijn1,ijn2,kln1,kln2,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
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
            call dxyz(qq(ic1),qq(ic2),qq(ic3),ncmmm)
            ioff = 0
            if (.not.(oskip(1))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,
     +                    ijn1,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(2))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,
     +                    ijn2,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(3))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,
     +                    kln1,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(4))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,
     +                    kln2,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
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
               call dnorm(qq(ic7+1),noform)
            end if
         end if
      end if
c
c    multiply by density matrix elements
c
      call grdcon(qq,noform)
      return
      end
      subroutine dnorm(qq,noform)
      implicit real*8  (a-h,o-z)
      dimension qq(*)
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
      parameter (mxp2 = mxprms*mxprms)
c
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      logical noform
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      dimension noform(*)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
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
      subroutine dpdab1(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit real*8  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/blk1/gin(510),mword,nlll1,kwor1,kwor11
      common/sortpk/labs(1360)
c
      call vclr(abdens,1,lendd)
      ijdiff = ii.ne.jj
      if (kwor1.eq.0) return
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor1)
         iword = 0
         call unpack(gin(num2e+1),lab816,labs,numlab)
         if (kwor1.le.0) go to 30
      end if
      iword = iword + 1
      m = (iword+iword) + (iword+iword)
      i = labs(m-3)
      j = labs(m-2)
      k = labs(m-1)
      l = labs(m)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dpdab2(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit real*8  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/three/labs(1360)
      common/bufc/gin(510),mword,nlll2,kwor2,kwor22
c
      call vclr(abdens,1,lendd)
      if (kwor2.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor2)
         iword = 0
         call unpack(gin(num2e+1),lab816,labs,numlab)
         if (kwor2.le.0) go to 30
      end if
      iword = iword + 1
      m = (iword+iword) + (iword+iword)
      i = labs(m-3)
      j = labs(m-2)
      k = labs(m-1)
      l = labs(m)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dpdab3(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit real*8  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/lsort/labs(1360)
      common/bufd/gin(510),mword,nlll3,kwor3,kwor33
c
      call vclr(abdens,1,lendd)
      if (kwor3.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor3)
         iword = 0
         call unpack(gin(num2e+1),lab816,labs,numlab)
         if (kwor3.le.0) go to 30
      end if
      iword = iword + 1
      m = (iword+iword) + (iword+iword)
      i = labs(m-3)
      j = labs(m-2)
      k = labs(m-1)
      l = labs(m)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dprim
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
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
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +            aej(ncmax),aek(ncmax),ael(ncmax),
     +            a(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),dd(4*mxp2),
     +            ijden(225)
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
            aa = ai + aj
            ainv = one/aa
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
            a(nij) = aa
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
      subroutine dshell(nelec,ish,jsh,ksh,lsh)
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     +,aek(ncmax),ael(ncmax),
     + aaa(9*mxp2),ijden(225),
     + ik(225),ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
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
c
      if (nelec.eq.2) then
         okandl = ksh.eq.lsh
         osame = ish.eq.ksh .and. jsh.eq.lsh
         k = katom(ksh)
         xk = c(1,k)
         yk = c(2,k)
         zk = c(3,k)
         k1 = kstart(ksh)
         k2 = k1 + kng(ksh) - 1
         lkt = ktype(ksh)
         mink = kmin(ksh)
         maxk = kmax(ksh)
         lock = kloc(ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex(k)
            csc(ngc) = cs(k)
            cpc(ngc) = cp(k)
            cdc(ngc) = cd(k)
            cfc(ngc) = cf(k)
            cgc(ngc) = cg(k)
 20      continue
         l = katom(lsh)
         xl = c(1,l)
         yl = c(2,l)
         zl = c(3,l)
         l1 = kstart(lsh)
         l2 = l1 + kng(lsh) - 1
         llt = ktype(lsh)
         minl = kmin(lsh)
         maxl = kmax(lsh)
         locl = kloc(lsh) - minl
         ngd = 0
         do 30 l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex(l)
            csd(ngd) = cs(l)
            cpd(ngd) = cp(l)
            cdd(ngd) = cd(l)
            cfd(ngd) = cf(l)
            cgd(ngd) = cg(l)
 30      continue
         nroots = (lit+ljt+lkt+llt-1)/2
         rrk = ((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
c
c     determine various offsets and indexing arrays
c
         inc2 = 1
         inc3 = inc2*(maxl-minl+1)
         inc4 = inc3*(maxk-mink+1)
         inc5 = inc4*(maxj-minj+1)
         lendd = inc5*(maxi-mini+1)
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
         return
      else
         oiandj = ish.eq.jsh
         i = katom(ish)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ish)
         i2 = i1 + kng(ish) - 1
         lit = ktype(ish)
         mini = kmin(ish)
         maxi = kmax(ish)
         loci = kloc(ish) - mini
         nga = 0
         do 80 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex(i)
            csa(nga) = cs(i)
            cpa(nga) = cp(i)
            cda(nga) = cd(i)
            cfa(nga) = cf(i)
            cga(nga) = cg(i)
 80      continue
         j = katom(jsh)
         xj = c(1,j)
         yj = c(2,j)
         zj = c(3,j)
         j1 = kstart(jsh)
         j2 = j1 + kng(jsh) - 1
         ljt = ktype(jsh)
         minj = kmin(jsh)
         maxj = kmax(jsh)
         locj = kloc(jsh) - minj
         ngb = 0
         do 90 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex(j)
            csb(ngb) = cs(j)
            cpb(ngb) = cp(j)
            cdb(ngb) = cd(j)
            cfb(ngb) = cf(j)
            cgb(ngb) = cg(j)
 90      continue
         rri = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
         return
      end if
      end
      subroutine dxyz(x,y,z,ncdim)
      implicit real*8  (a-h,o-z)
      parameter (ncmax=65)
      dimension x(*),y(*),z(*)
      logical n0,n1,m0,m1
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     1   xc00(ncmax),
     2ycp00(ncmax),yc00(ncmax),zcp00(ncmax),zc00(ncmax),f00(ncmax),dxij,
     3dyij,dzij,dxkl,dykl,dzkl,iorg(12),korg(12),
     4nimax,njmax,nkmax,nlmax,nmax,mmax,ij1x,ij2x,kl1x,kl2x
      common/small/ca(ncmax),cb(ncmax)
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
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
      subroutine fockd2(qq,noform)
c============================================
c     derivative fock operators, j and k matrices
c     etc. for use in later chf equations
c================================================
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
c
c     derivative fock operators
c     additions to write integral derivatives
c
      logical noform
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  xij(225),xkl(225),ijgt(225),klgt(225)
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
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
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
      dimension qq(*),noform(*)
      data four /4.0d0/
      data small/1.0d-14/
c
c
      if (mp2w) then
         ioa = ic7 + lendd*npass*3
         call vclr(qq(ioa+1),1,3*lendd)
         do 30 i = 1 , npass
            iob = ic7 + (i-1)*lendd*3
            do 20 k = 1 , lendd*3
               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
 20         continue
 30      continue
      end if
      idb = ida + nx
      np1 = npass + 1
      nw3 = nx*3
      lenbig = nat*nw3
      len3 = lendd*3
      ntr = ifok + (natomd(np1)-1)*nw3
c
      hf_wght = CD_HF_exchange_weight()
c
      if (ogrhf) then
c
c         general case
c
         nn = 0
         nnn = 0
         jmax = maxj
         do 100 i = mini , maxi
            i1 = loci + i
            ia1 = iky(i1)
            if (oiandj) jmax = i
            do 90 j = minj , jmax
               i2 = locj + j
               ia2 = iky(i2)
               mij = ia1 + i2
               lmax = maxl
               do 80 k = mink , maxk
                  i3 = lock + k
                  ia3 = iky(i3)
                  mik = ia1 + i3
                  if (i3.gt.i1) mik = ia3 + i1
                  mjk = ia2 + i3
                  if (i3.gt.i2) mjk = ia3 + i2
                  if (okandl) lmax = k
                  do 70 l = minl , lmax
                     i4 = locl + l
                     ia4 = iky(i4)
                     mkl = ia3 + i4
                     mil = ia1 + i4
                     if (i4.gt.i1) mil = ia4 + i1
                     mjl = ia2 + i4
                     if (i4.gt.i2) mjl = ia4 + i2
                     nnn = nnn + 1
                     if (.not.(noform(nnn))) then
                        nn = nn + 1
c
                        do 60 npp = 1 , npass
                           ifa = ifok + (natomd(npp)-1)*nw3
                           itrf = ntr
                           ioa = ic7 + (npp-1)*len3
c
                           do 50 icomp = 1 , 3
                              g1 = qq(nn+ioa)
c
                              if (dabs(g1).ge.small) then
                                 g1 = g1 + g1
                                 g14 = g1 + g1
                                 itj = itrf
                                 ifj = ifa
                                 ifd = ida
                                 do 40 ns = 1 , ndens
                                    ifk = ifj + lenbig
                                    itk = itj + lenbig
                                    gdkl = g14*qq(ifd+mkl)
                                    qq(ifj+mij) = qq(ifj+mij) + gdkl
                                    gdij = g14*qq(ifd+mij)
                                    qq(ifj+mkl) = qq(ifj+mkl) + gdij
                                    gdjl = g1*qq(ifd+mjl)
                                    qq(ifk+mik) = qq(ifk+mik) + gdjl
                                    gdjk = g1*qq(ifd+mjk)
                                    qq(ifk+mil) = qq(ifk+mil) + gdjk
                                    gdil = g1*qq(ifd+mil)
                                    qq(ifk+mjk) = qq(ifk+mjk) + gdil
                                    gdik = g1*qq(ifd+mik)
                                    qq(ifk+mjl) = qq(ifk+mjl) + gdik
c     use translational invariance to get missing term
                                    qq(itj+mij) = qq(itj+mij) - gdkl
                                    qq(itj+mkl) = qq(itj+mkl) - gdij
                                    qq(itk+mik) = qq(itk+mik) - gdjl
                                    qq(itk+mil) = qq(itk+mil) - gdjk
                                    qq(itk+mjk) = qq(itk+mjk) - gdil
                                    qq(itk+mjl) = qq(itk+mjl) - gdik
                                    ifj = ifj + lenbig + lenbig
                                    itj = itj + lenbig + lenbig
                                    ifd = ifd + nx
 40                              continue
                              end if
                              ifa = ifa + nx
                              itrf = itrf + nx
                              ioa = ioa + lendd
 50                        continue
 60                     continue
                     end if
 70               continue
 80            continue
 90         continue
 100     continue
      else
         nn = 0
         nnn = 0
         jmax = maxj
         do 160 i = mini , maxi
            i1 = loci + i
            ia1 = iky(i1)
            if (oiandj) jmax = i
            do 150 j = minj , jmax
               i2 = locj + j
               ia2 = iky(i2)
               mij = ia1 + i2
               dij = four*qq(ida+mij)
               lmax = maxl
               do 140 k = mink , maxk
                  i3 = lock + k
                  ia3 = iky(i3)
                  mik = ia1 + i3
                  if (i3.gt.i1) mik = ia3 + i1
                  mjk = ia2 + i3
                  if (i3.gt.i2) mjk = ia3 + i2
                  dik = qq(ida+mik)*hf_wght
                  djk = qq(ida+mjk)*hf_wght
                  if (okandl) lmax = k
                  do 130 l = minl , lmax
                     i4 = locl + l
                     ia4 = iky(i4)
                     mkl = ia3 + i4
                     mil = ia1 + i4
                     if (i4.gt.i1) mil = ia4 + i1
                     mjl = ia2 + i4
                     if (i4.gt.i2) mjl = ia4 + i2
                     nnn = nnn + 1
                     if (.not.(noform(nnn))) then
                        nn = nn + 1
                        dkl = four*qq(ida+mkl)
                        dil = qq(ida+mil)*hf_wght
                        djl = qq(ida+mjl)*hf_wght
c
c
                        do 120 npp = 1 , npass
                           ifa = ifok + (natomd(npp)-1)*nw3
                           itrf = ntr
                           ioa = ic7 + (npp-1)*len3
                           do 110 icomp = 1 , 3
                              g1 = qq(nn+ioa)
                              if (dabs(g1).ge.small) then
                                 gdkl = g1*dkl
                                 qq(ifa+mij) = qq(ifa+mij) + gdkl
                                 gdij = g1*dij
                                 qq(ifa+mkl) = qq(ifa+mkl) + gdij
                                 gdjl = g1*djl
                                 qq(ifa+mik) = qq(ifa+mik) - gdjl
                                 gdjk = g1*djk
                                 qq(ifa+mil) = qq(ifa+mil) - gdjk
                                 gdil = g1*dil
                                 qq(ifa+mjk) = qq(ifa+mjk) - gdil
                                 gdik = g1*dik
                                 qq(ifa+mjl) = qq(ifa+mjl) - gdik
c     use translational invariance to get missing term
                                 qq(itrf+mij) = qq(itrf+mij) - gdkl
                                 qq(itrf+mkl) = qq(itrf+mkl) - gdij
                                 qq(itrf+mik) = qq(itrf+mik) + gdjl
                                 qq(itrf+mil) = qq(itrf+mil) + gdjk
                                 qq(itrf+mjk) = qq(itrf+mjk) + gdil
                                 qq(itrf+mjl) = qq(itrf+mjl) + gdik
                                 if (.not.(oclos)) then
                                    ika = ifa + lenbig
                                    itrk = itrf + lenbig
                                    gdjl = g1*qq(idb+mjl)
                                    qq(ika+mik) = qq(ika+mik) + gdjl
                                    gdil = g1*qq(idb+mil)
                                    qq(ika+mjk) = qq(ika+mjk) + gdil
                                    gdjk = g1*qq(idb+mjk)
                                    qq(ika+mil) = qq(ika+mil) + gdjk
                                    gdik = g1*qq(idb+mik)
                                    qq(ika+mjl) = qq(ika+mjl) + gdik
                                    qq(itrk+mik) = qq(itrk+mik) - gdjl
                                    qq(itrk+mjl) = qq(itrk+mjl) - gdik
                                    qq(itrk+mjk) = qq(itrk+mjk) - gdil
                                    qq(itrk+mil) = qq(itrk+mil) - gdjk
                                 end if
                              end if
                              ifa = ifa + nx
                              itrf = itrf + nx
                              ioa = ioa + lendd
 110                       continue
 120                    continue
                     end if
 130              continue
 140           continue
 150        continue
 160     continue
         if (mp2w) call derwrt(qq(ic7+1))
         return
      end if
      end
      subroutine fokabd(qq,noform)
c
c----------------------------------------------------
c    construct fock matrices from derivative integrals
c    and derivative density matrices.
c    used in polarizability derivative calculations
c-----------------------------------------------------
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
      logical noform
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +   aej(ncmax),aek(ncmax),ael(ncmax),
     +   aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +   dd(4*mxp2),ijden(225),ik(225),
     +   ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +   dij(225),dkl(225),ijgt(225),klgt(225)
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
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
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
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      dimension qq(*),noform(*)
      data small/1.0d-14/
c
      np1 = npass + 1
      nw3 = nx*3
      lenbig = nat*nw3
      len3 = lendd*3
      ntr = ifok + (natomd(np1)-1)*nw3
c
c       ndens density matrices for each nuclear displacement
c
      nn = 0
      nnn = 0
      jmax = maxj
      do 80 i = mini , maxi
         i1 = loci + i
         ia1 = iky(i1)
         if (oiandj) jmax = i
         do 70 j = minj , jmax
            i2 = locj + j
            ia2 = iky(i2)
            mij = ia1 + i2
            lmax = maxl
            do 60 k = mink , maxk
               i3 = lock + k
               ia3 = iky(i3)
               mik = ia1 + i3
               if (i3.gt.i1) mik = ia3 + i1
               mjk = ia2 + i3
               if (i3.gt.i2) mjk = ia3 + i2
               if (okandl) lmax = k
               do 50 l = minl , lmax
                  i4 = locl + l
                  ia4 = iky(i4)
                  mkl = ia3 + i4
                  mil = ia1 + i4
                  if (i4.gt.i1) mil = ia4 + i1
                  mjl = ia2 + i4
                  if (i4.gt.i2) mjl = ia4 + i2
                  nnn = nnn + 1
                  if (.not.(noform(nnn))) then
                     nn = nn + 1
c
c      loops over nuclei....
c
                     do 40 npp = 1 , npass
                        ifa = ifok + (natomd(npp)-1)*nw3
                        itrf = ntr
                        ioa = ic7 + (npp-1)*len3
                        do 30 icomp = 1 , 3
                           g1 = qq(nn+ioa)
c
                           if (dabs(g1).ge.small) then
                              g14 = (g1+g1) + (g1+g1)
                              itf = itrf
                              iff = ifa
                              ifd = ida + nx
                              do 20 ns = 1 , ndenin
                                 gdkl = g14*qq(ifd+mkl)
                                 qq(iff+mij) = qq(iff+mij) + gdkl
                                 gdij = g14*qq(ifd+mij)
                                 qq(iff+mkl) = qq(iff+mkl) + gdij
                                 gdjl = g1*qq(ifd+mjl)
                                 qq(iff+mik) = qq(iff+mik) - gdjl
                                 gdjk = g1*qq(ifd+mjk)
                                 qq(iff+mil) = qq(iff+mil) - gdjk
                                 gdil = g1*qq(ifd+mil)
                                 qq(iff+mjk) = qq(iff+mjk) - gdil
                                 gdik = g1*qq(ifd+mik)
                                 qq(iff+mjl) = qq(iff+mjl) - gdik
c     use translational invariance to get missing term
                                 qq(itf+mij) = qq(itf+mij) - gdkl
                                 qq(itf+mkl) = qq(itf+mkl) - gdij
                                 qq(itf+mik) = qq(itf+mik) + gdjl
                                 qq(itf+mil) = qq(itf+mil) + gdjk
                                 qq(itf+mjk) = qq(itf+mjk) + gdil
                                 qq(itf+mjl) = qq(itf+mjl) + gdik
                                 iff = iff + lenbig
                                 itf = itf + lenbig
                                 ifd = ifd + nx
 20                           continue
                           end if
                           ifa = ifa + nx
                           itrf = itrf + nx
                           ioa = ioa + lendd
 30                     continue
 40                  continue
                  end if
 50            continue
 60         continue
 70      continue
 80   continue
      return
      end
      subroutine formeg
      implicit real*8  (a-h,o-z)
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
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
c
      real*8 de
      common/grad2/de(3,maxat)
c
      common/tgrad/dgout(3,3)
      data zero /0.0d0/
      dumx = zero
      dumy = zero
      dumz = zero
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
      subroutine grdcon(qq,noform)
      implicit real*8  (a-h,o-z)
      logical nofk,noform
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
      parameter (mxp2 = mxprms*mxprms)
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
      common/tgrad/dgout(9)
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension noform(*),qq(*)
      data half/0.5d0/
      nofk = .not.(ofock .or. ofokab)
c
c
      if ((mp2w .or. ompir)) then
         ioa = ic7 + lendd*npass*3
         call vclr(qq(ioa+1),1,3*lendd)
         do 30 i = 1 , npass
            iob = ic7 + (i-1)*lendd*3
            do 20 k = 1 , lendd*3
               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
 20         continue
 30      continue
      end if
      if (ompir) then
         n2 = iabd + lendd + 1
         do 60 i = 1 , 3
            n1 = ic7 + 1
            do 50 k = 1 , npass + 1
               do 40 j = 1 , 3
                  t1 = ddot(ijkld,qq(n2),1,qq(n1),1)
                  dipi(i,j,natomd(k)) = dipi(i,j,natomd(k)) - t1
                  n1 = n1 + lendd
 40            continue
 50         continue
            n2 = n2 + lendd
 60      continue
      end if
c
      if (outd .or. ofock .or. ofokab) then
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
                           if (.not.(nofk)) then
                              if (i1.eq.j1) then
                                 dumx = dumx*half
                                 dumy = dumy*half
                                 dumz = dumz*half
                              end if
                              if (k1.eq.l1) then
                                 dumx = dumx*half
                                 dumy = dumy*half
                                 dumz = dumz*half
                              end if
                              qq(n1+nnn) = q4*dumx
                              qq(n2+nnn) = q4*dumy
                              qq(n3+nnn) = q4*dumz
                           end if
                           n0 = n0 + 3
                           n1 = n3 + lendd
 70                     continue
                     end if
 80               continue
 90            continue
 100        continue
 110     continue
      else
         n1 = ic7 + 1
         np3 = npass*3
         do 120 n0 = 1 , np3
            dgout(n0) = dgout(n0) + ddot(ijkld,qq(iabd+1),1,qq(n1),1)
            n1 = n1 + lendd
 120     continue
      end if
 6010 format (1x,12i6/30x,3e16.8,5x,e16.8)
      end
      subroutine hsymd(q,iwr)
c
c----------------------------------------------
c   write out derivative fock operators
c---------------------------------------------
c
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
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
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     + aej(ncmax),aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
      dimension q(*)
      data half/0.5d0/
c
c     writes out accumulated derivative fock and k matrices
c     onto hamfile.
c
      if (odebug(20)) write (iwr,6010)
      ioff = ifok
      ifocbl = iochf(11)
      if (ogrhf) ifocbl = ifocbl + nat3*lenb

      do 50 n = 1 , nfok
         call rdedx(q(ida+1),nx,ifocbl,ifockf)
         ij = 0
         do 30 k = 1 , num
            do 20 l = 1 , k
               ij = ij + 1
               if (k.ne.l) then
                  q(ioff+ij) = half*q(ioff+ij)
               end if
 20         continue
 30      continue
         do 40 k = 1 , nx
            q(ida+k) = q(ida+k) + q(ioff+k)
 40      continue

         call wrt3(q(ida+1),nx,ifocbl,ifockf)
         if (odebug(20)) call prtris(q(ida+1),num,iwr)
         ioff = ioff + nx
         ifocbl = ifocbl + lenb
c
c
 50   continue
      lfdtrn = .false.
      if(odebug(30)) write (iwr,6020) iochf(11), iochf(12)
      return
 6010 format (/1x,'**** output from hsymd ****'/)
 6020 format (/1x,'hamfile summary'/
     +         1x,'section 11 at block ',i5/
     +         1x,'section 12 at block ',i5/)
      end
      subroutine jkder(zscftp,q,prefac,iso,nshels)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      parameter (mxp2 = mxprms*mxprms)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w
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
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
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
      logical ofokab,ompir,ofock
      integer ndenin,iflden,ibden,iflout,ibdout, ntpdm
      integer nat3,ndens,nfok,nbsq,lenb
      common/specal/ndenin,iflden,ibden,iflout,ibdout,ofokab,ntpdm,
     1  ompir,ofock,nat3,ndens,nfok,nbsq,lenb
      common/blkin/gijkl(510),mword,nlenx,kworx,kworxx
      common/craypk/ijkl(1360)
      common/blk1/gijkl1(510),mwor1,nlen1,kwor1,kwor11
      common/sortpk/ijkl1(1360)
      common/bufc/gijkl2(510),mwor2,nlen2,kwor2,kwor22
      common/three/ijkl2(1360)
      common/bufd/gijkl3(510),mwor3,nlen3,kwor3,kwor33
      common/lsort/ijkl3(1360)
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/iao(5,3120)
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
      parameter (ncmax=65)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/tgrad/dgout(9)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     + aej(ncmax),aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
c
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      integer icount_dlb, icount_slb
      common /par_dlb/ icount_dlb, icount_slb
c
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
      dimension iso(nshels,*),q(*)
      dimension m0(48),m1(48),m2(48),m3(48)
      dimension prefac(*)
      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zvb,zmcscf/'casscf','vb','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
c
      indq(m,n) = (m-1)*n
c
      ompir = runtyp.eq.zdipd .or. runtyp.eq.zinfra
      ompir = ompir .and. omp2
c
c     ----- check for grhf or gvb cases -----
c
      ouhf = zscftp.eq.zuhf
      orgvb = zscftp.eq.zgvb
      ogrhf = zscftp.eq.zgrhf
      oclos = zscftp.eq.zrhf
      omcscf = zscftp.eq.zmcscf
c
c     ----- check for casscf
c
      ocas = zscftp.eq.zcas.or.zscftp.eq.zvb
      if (ocas .or. mp3 .or. (omp2 .and. .not.ompir)) then
         call setsto(1360,0,ijkl)
      end if
      if (nprint.ne.-5 .and. oprint(57)) then
         write (iwr,6010)
         call writel(prefac,nshels)
      end if
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()
      ntpdm = 1
      oeof = .false.
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      if (ompir) then
         ndenin = 3
         ntpdm = 4
         iflden = idaf
         call secget(isect(31),31,ibden)
         iwor1 = 0
         iwor2 = 0
         iwor3 = 0
         mwor1 = 0
         mwor2 = 0
         mwor3 = 0
         ifmp1 = 20
         ifmp2 = 21
         ifmp3 = 22
         ib1 = 1
         call search(ib1,ifmp1)
         call search(ib1,ifmp2)
         call search(ib1,ifmp3)
         call setsto(1360,0,ijkl1)
         call setsto(1360,0,ijkl2)
         call setsto(1360,0,ijkl3)
         call secget(isect(57),57,iblok)
         call rdedx(dipd,lds(isect(57)),iblok,ifild)
      end if

      if (ogrhf) then
         m = 0
         call secget(isect(53),m,iblok)
         call readi(nact,lds(isect(53))*nav,iblok,idaf)
      end if
      nat3 = nat*3
      nbsq = num*num
      lenb = lensec(nx)
      odpres = .false.
      ofpres = .false.
      ogpres = .false.
      do 20 i = 1 , nshels
         if (ktype(i).eq.3) odpres = .true.
         if (ktype(i).eq.4) ofpres = .true.
         if (ktype(i).eq.5) ogpres = .true.
 20   continue
c
      m = ntpdm + 9
      if (omp2w .or. ompir) m = m + 3
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
c code for dgenrl preallocation 
c   NOW all version (previously GA version only)
c
c add buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do 21 i = 1 , nshels
         kt_max=max(kt_max,ktype(i))
 21   continue
      inc1_max=(kt_max+1)**4
c
c vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6
c
c     ----- set pointers for partitioning of core -----
c
c
c revised scheme - 
c

      i10 = igmem_alloc_inf(maxorb,'drv2e.m','jkder','mo2fock',
     &                      IGMEM_DEBUG)
cc      i20 = igmem_alloc(l2)
cc      i30 = igmem_alloc(l2)

      ndens = 1
      if (ouhf) then
         ndens = 2
      end if
      if (orgvb) then
         ndens = 3
      end if
      if (ogrhf) then
         ndens = njk
      end if
      if (ocas) then
         ndens = 1
      end if
      if (ofock .or. ompir) then
         ndens = ndens + ndenin
      end if

      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1

      ida = igmem_alloc_inf(( ndens)*l2,'drv2e.m','jkder','dens-mats',
     &                      IGMEM_NORMAL)
      i20 = ida
      i30 = i20 + l2
      ida = ida - 1

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
         if (ogrhf) nfok = njk*nat*6
      end if
      if (ofokab) nfok = nat3*ndenin
c

      ifok = igmem_alloc_inf(nx*nfok,'drv2e.m','jkder','der-fock',
     &                       IGMEM_NORMAL)
      ifok = ifok - 1

c  core at i00 is for indexing
      i00 = igmem_alloc_inf(lnddm*(3/nav+1),'drv2e.m','jkder',
     &                      'oform-int',IGMEM_DEBUG)

      iabd = igmem_alloc_inf(lnddm*ntpdm,'drv2e.m','jkder',
     &                       'dens-mat-prod',IGMEM_NORMAL)
      iabd = iabd - 1

      ic7 = igmem_alloc_inf(lnddm*9 + inc1_max*ncmmm_max*6,'drv2e.m',
     &                      'jkder','tmp-int',IGMEM_NORMAL)
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

      if(omcscf)then

         i40 = igmem_alloc(l2)

ccc         id3 = i40 + nx
         numcmo = num*ncoorb
         nactp = ncact*(ncact+1)/2
         nd2str = indq(nactp+1,16)
         if (odpres) nd2str = indq(nactp+1,36)
         if (ofpres) nd2str = indq(nactp+1,100)
         if (ogpres) nd2str = indq(nactp+1,225)

         id3 = igmem_alloc(numcmo)
         i30 = i20
ccc         id4 = id3 + numcmo

c     nd2mo=ind(nactp,nactp)
         nd2mo = nactp*(nactp+1)/2
c     write(6,*)' length of tpdm',nd2mo

         id4 = igmem_alloc(nd2str)
ccc         id5 = id4 + nd2str
         id5 = igmem_alloc(nd2mo)
ccc         id6 = id5 + nd2mo
c
         najkl = indq(ncact*4+1,nactp)
         if (odpres) najkl = indq(ncact*6+1,nactp)
         if (ofpres) najkl = indq(ncact*10+1,nactp)
         if (ogpres) najkl = indq(ncact*15+1,nactp)
c
         nabcl = indq(16+1,4*ncact)
         if (odpres) nabcl = indq(36+1,6*ncact)
         if (ofpres) nabcl = indq(100+1,10*ncact)
         if (ogpres) nabcl = indq(225+1,15*ncact)
c
         id6 = igmem_alloc(nabcl)
cc         id7 = id6 + nabcl
         id7 = igmem_alloc(najkl)
cc         id8 = id7 + najkl
         id8 = igmem_alloc(nactp)
cc         id9 = id8 + nactp
c
         nshdm = max(ncact,4)
         if (odpres) nshdm = max(ncact,6)
         if (ofpres) nshdm = max(ncact,10)
         if (ogpres) nshdm = max(ncact,15)
c
         id9 = igmem_alloc(nshdm*nshdm)
cc
cc         iabd = id9 + nshdm*nshdm
c
c   iabd should not be used for mcscf calculations - tpdm in id5
c
cc         i00 = iabd + lnddm*ntpdm
cc         ic7 = i00 + lnddm*(3/nav+1) - 1
      endif


c     tim0 = cpulft(1)
c
      if (ocas) then
         call dbutci(ist,jst,kst,lst)
      else if (omcscf) then
         call ddebut(zscftp,q(id3),q(i30),q(i10),q(i40),q(id5),l1,
     +               l2,l3,ist,jst,kst,lst,q)
      else
         call ddebut(zscftp,q(i20),q(i30),q(i10),q(i10),q(ifok+1),
     +               l1,l2,l3,ist,jst,kst,lst,q)
      end if
      if (mp3 .or. (omp2 .and. .not.ompir)) call search(iblk2d,ifil2d)
      mword = 0
      iword = 0
      kworx = 999
      kwor1 = 999
      kwor2 = 999
      kwor3 = 999
      kloc(nshels+1) = num + 1
c
      icnt = 0
      ib1 = 1
      if (omp2w) then
         call search(ib1,mpstrm(1))
         if (odebug(30)) write (iwr,6020) mpstrm(1)
      end if
      if (ist.le.nshels) then
c
c     ----- ishell -----
c
         do 140 ii = ist , nshels
            kadi = kad(ii)
            ijshel = ii*(ii-1)/2
            do 40 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 140
               m0(it) = id
 40         continue
            iceni = katom(ii)
            if (omcscf) call mcajkl(q(id3),q(id5),q(id7),q(id8),q(id9),
     +                              ii,nactp)
c
c     ----- jshell -----
c
            j0 = jst
            do 130 jj = j0 , ii
               kadij = kadi + kad(jj)
               jst = 1
               tolij = dlntol + prefac(ijshel+jj)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
               if (tolij.gt.-3.401d0) then
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
                  olab = katom(jj).eq.iceni
c
c     store information about the pair (ij)
c
                  call dshell(1,ii,jj,kk,ll)
                  call dprim
                  if (nij.ne.0) then
                     if (omcscf) call mcabkl(q(id3),q(id7),q(id4),ii,jj,
     +                   nactp)
c
c     ----- kshell -----
c
                     k0 = kst
                     do 120 kk = k0 , ii
                        kadijk = kadij + kad(kk)
                        kst = 1
                        klshel = kk*(kk-1)/2
                        do 80 it = 1 , nt
                           kd = iso(kk,it)
                           if (kd.gt.ii) go to 120
                           m3(it) = kd
 80                     continue
                        olabc = olab .and. katom(kk).eq.iceni
                        if (omcscf)
     +                      call mcabcl(q(id3),q(id4),q(id6),q(id8),
     +                      q(id9),ii,jj,kk,nactp)
                        l0 = lst
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                        do 110 ll = l0 , maxll
                           lst = 1
                           if (kadijk+kad(ll).lt.0) then
                              tolijk = tolij + prefac(klshel+ll)
                              if (tolijk.gt.0.0d0) then
                                 olabcd = olabc .and. katom(ll).eq.iceni
                                 if (.not.(olabcd)) then
                                    n4 = 0
                                    do 100 it = 1 , nt
                                       ld = iso(ll,it)
                                       if (ld.gt.ii) go to 110
                                       kd = m3(it)
                                       if (kd.lt.ld) then
                                         nd = kd
                                         kd = ld
                                         ld = nd
                                       end if
                                       id = m1(it)
                                       jd = m2(it)
                                       if (id.eq.ii .or. kd.eq.ii) then
                                         if (kd.ge.id) then
                                         if (kd.ne.id .or. ld.gt.jd)
     +                                      then
                                         nd = id
                                         id = kd
                                         kd = nd
                                         nd = jd
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
c
c     ----- check for redundant combinations -----
c
                              call redund(ii,jj,kk,ll,iwr)
                              if (npass.ne.0) then
c
c     ----- initialize dgout to zero -----
c
                          call vclr(dgout,1,ndgout)
                                 call dshell(2,ii,jj,kk,ll)
c
c     ----- form products of density matrix elements -----
c
                          call vclr(q(iabd+1),1,lendd*ntpdm)
                                 if (omcscf)
     +                              call mcabcd(q(id3),q(id6),
     +                              q(iabd+1),q(id9),ii,jj,kk,ll,
     +                              q4)
                                 if (mp3 .or.
     +                              (omp2 .and. .not.ompir))
     +                              call mcdab(q(iabd+1),ii,jj,kk,
     +                              ll,q4)
                                 if (ompir) then
                                   iabd1 = iabd + lendd + 1
                                   iabd2 = iabd1 + lendd
                                   iabd3 = iabd2 + lendd
                                   call dpdab1(q(iabd1),ii,jj,kk,
     +                                ll,q4,ifmp1,iwor1)
                                   call dpdab2(q(iabd2),ii,jj,kk,
     +                                ll,q4,ifmp2,iwor2)
                                   call dpdab3(q(iabd3),ii,jj,kk,
     +                                ll,q4,ifmp3,iwor3)
                                   call tpdder(q(iabd1),lendd,
     +                                q(i20),q(i30),nx,ndenin,ii,
     +                                jj,kk,ll,q4)
                                 end if
                                 if (.not.orgvb .and.
     +                              .not.ocas .and. .not.omcscf)
     +                              call dabab(ii,jj,kk,ll,q4,
     +                              zscftp,q(i20),q(i30),q(iabd+1)
     +                              )
                                 if (orgvb)
     +                              call dabg(ii,jj,kk,ll,l1,norb,
     +                              q4,q(i20),q(i30),q(i10)
     +                              ,onocor,onopen,q(iabd+1))
                                 if (ocas)
     +                              call dabci(ii,jj,kk,ll,q4,
     +                              oeof,q(iabd+1))
                                 if (omcscf)
     +                              call dabmc(ii,jj,kk,ll,q4, 
     +                              q(i30),q(i40),q(iabd+1))
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
                                 call delim(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)
                                 if (ompir) then
                                   call delim2(q(iabd1),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                   call delim2(q(iabd2),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                   call delim2(q(iabd3),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                 end if
                                 if (ijkld.ne.0) then
                                   call dgenrl(q(1),q(i00),q(i00),
     +                                abmax)
                                   if (ofock)
     +                                call fockd2(q,q(i00))
                                   if (ofokab) then
                                   call fokabd(q,q(i00))
                                   end if
c
c     ----- generate all 4 partial contributions to the gradient ----
c
                                   call formeg
                                 end if
                              end if
                                 end if
                              end if
                           end if
 110                    continue
 120                 continue
                  end if
               end if
 130        continue
c
c     ----- save gradient and restart data -----
c
            call dfinal(q,0,ii)
            if (tim.ge.timlim) go to 150
 140     continue
      end if
c
c     ----- end of *shell* loops -----
c
      if (omp2w) then
         mderb = idblk + iidblk + 2
         if (icnt.ne.0) then
            call pack(val(idblk1),lab816,iao,idblk*5)
            call wrt3s(dbufs,mderb,mpstrm(1))
         end if
         icnt = 0
         call wrt3s(dbufs,mderb,mpstrm(1))
         call shut1(mpstrm(1))
         if (odebug(30)) write (iwr,6030) mpstrm(1)
      end if
      if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
       isymd = igmem_alloc_inf(nw196(6),'drv2e.m','jkder',
     &                         'trans-atoms-symm',IGMEM_DEBUG)
c
       call symde(q(isymd),nat)
c
c     ----- reset core memory from symde
c
      call gmem_free_inf(isymd,'drv2e.m','jkder','trans-atoms-symm')
c
      endif
      call dfinal(q,1,ii)
      nindmx = 0
 150  call timit(0)
c     dtim = tim - tim0
c
c     ----- reset core memory from jkder -----
c

      if(omcscf)then
         call gmem_free(id9)
         call gmem_free(id8)
         call gmem_free(id7)
         call gmem_free(id6)
         call gmem_free(id5)
         call gmem_free(id4)
         call gmem_free(id3)
         call gmem_free(i40)


      endif
cc
cc revert
cc      call gmem_free(ic1)
      ic7 = ic7 + 1
      call gmem_free_inf(ic7,'drv2e.m','jkder','tmp-int')
      iabd = iabd + 1
      call gmem_free_inf(iabd,'drv2e.m','jkder','dens-mat-prod')
      call gmem_free_inf(i00,'drv2e.m','jkder','oform-int')
      ifok = ifok + 1
      call gmem_free_inf(ifok,'drv2e.m','jkder','der-fock')
      ida = ida + 1
      call gmem_free_inf(ida,'drv2e.m','jkder','dens-mats')
      call gmem_free_inf(i10,'drv2e.m','jkder','mo2fock')
      return
 6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
 6020 format (/1x,'derivative integrals to be output to unit',i3)
 6030 format (/1x,'derivative integrals written to unit',i3)
      end
      subroutine mcabcd(cc,gabcl,gabcd,gsq,ish,jsh,ksh,lsh,q4)
      implicit real*8  (a-h,o-z)
      dimension cc(*),gabcl(*),gabcd(*),gsq(*)
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      data pt5/0.5d0/
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      minc = kloc(ksh)
      mind = kloc(lsh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      maxc = kloc(ksh+1) - 1
      maxd = kloc(lsh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nc = maxc - minc + 1
      nd = maxd - mind + 1
c     nab = na*nb
      ncd = nc*nd
      ncl = nc*ncact
      lmin = ncore*num + mind
      qq4 = 8.0d0*q4
c
      iab = 0
      do 50 ib = minb , maxb
         ib1 = (ib-minb)*inc4
         do 40 ia = mina , maxa
            iaa = (ia-mina)*inc5
            ibb = ib1 + iaa
            iab = iab + 1
c
            dd4 = qq4
            if (ia.eq.ib) dd4 = qq4*pt5
            iabcl = indq(iab,ncl) + 1
            call vclr(gsq(1),1,ncd)
            call mxmb(gabcl(iabcl),1,nc,cc(lmin),num,1,gsq,1,nc,nc,
     +                ncact,nd)
c
            icd = 0
            do 30 id = mind , maxd
               id1 = id - mind + inc2
               do 20 ic = minc , maxc
                  icc = (ic-minc)*inc3 + ibb
                  idd = id1 + icc
c
                  ddd4 = dd4
                  if (ic.eq.id) ddd4 = ddd4*pt5
                  icd = icd + 1
                  gabcd(idd) = gsq(icd)*ddd4
 20            continue
 30         continue
 40      continue
 50   continue
      return
      end
      subroutine mcabcl(cc,gabkl,gabcl,gtr,gsq,ish,jsh,ksh,nkl)
      implicit real*8  (a-h,o-z)
      dimension cc(*),gabkl(*),gabcl(*),gtr(*),gsq(*)
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
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      minc = kloc(ksh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      maxc = kloc(ksh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nc = maxc - minc + 1
      nab = na*nb
      ncl = nc*ncact
      kkmin = ncore*num + minc
c
      iab = 0
      call vclr(gabcl(1),1,nab*ncl)
      do 40 ib = 1 , nb
         do 30 ia = 1 , na
            iab = iab + 1
            icol2 = indq(iab,ncl) + 1
            ikl = iab
c
            do 20 kl = 1 , nkl
               gtr(kl) = gabkl(ikl)
               ikl = ikl + nab
 20         continue
            call squr(gtr(1),gsq(1),ncact)
            call mxmb(cc(kkmin),1,num,gsq,1,ncact,gabcl(icol2),1,nc,nc,
     +                ncact,ncact)
 30      continue
 40   continue
      return
      end
      subroutine mcabkl(cc,gajkl,gabkl,ish,jsh,nkl)
      implicit real*8  (a-h,o-z)
      dimension cc(*),gajkl(*),gabkl(*)
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
c
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nab = na*nb
      naj = na*ncact
      jmin = ncore*num + minb
c
      call vclr(gabkl(1),1,nab*nkl)
      do 20 kl = 1 , nkl
         icol1 = indq(kl,naj) + 1
         icol2 = indq(kl,nab) + 1
         call mxmb(gajkl(icol1),1,na,cc(jmin),num,1,gabkl(icol2),1,na,
     +             na,ncact,nb)
 20   continue
      return
      end
      subroutine mcajkl(cc,gijkl,gajkl,gtr,gsq,ish,nkl)
      implicit real*8  (a-h,o-z)
      dimension cc(*),gijkl(*),gajkl(*),gtr(*),gsq(*)
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
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      maxa = kloc(ish+1) - 1
      na = maxa - mina + 1
      naj = ncact*na
      imin = ncore*num + mina
c
      call vclr(gajkl(1),1,naj*nkl)
      do 30 kl = 1 , nkl
         icol1 = ind(kl,1)
         icol2 = indq(kl,naj) + 1
         call dcopy(kl,gijkl(icol1),1,gtr(1),1)
c
         icol1 = kl
         do 20 ij = kl + 1 , nkl
            ijkl = ind(ij,kl)
            gtr(icol1+1) = gijkl(ijkl)
            icol1 = icol1 + 1
 20      continue
c
         call squr(gtr(1),gsq(1),ncact)
         call mxmb(cc(imin),1,num,gsq,1,ncact,gajkl(icol2),1,na,na,
     +        ncact,ncact)
 30   continue
      return
      end
      subroutine mcdab(abdens,ii,jj,kk,ll,q4)
c================================================
c    read two particle density matrix from file
c================================================
      implicit real*8  (a-h,o-z)
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
      logical ijdiff,kldiff,ikdiff
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labs(1360)
      common/blkin/gin(510),mword,nlenx,kworx,kworxx
c
      dimension abdens(*)
      call vclr(abdens,1,lendd)
      if (kworx.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kworx)
         iword = 0
         call unpack(gin(num2e+1),lab816,labs,numlab)
         if (kworx.le.0) go to 30
      end if
      iword = iword + 1
      m = (iword+iword) + (iword+iword)
      i = labs(m-3)
      j = labs(m-2)
      k = labs(m-1)
      l = labs(m)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine redund(ii,jj,kk,ll,iw)
      implicit real*8  (a-h,o-z)
      logical inej,inek,inel,jnek,jnel,knel
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
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
      lit = ktype(ii)
      ljt = ktype(jj)
      lkt = ktype(kk)
      llt = ktype(ll)
      iat = katom(ii)
      jat = katom(jj)
      kat = katom(kk)
      lat = katom(ll)
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
     +        ' npass =',i2,' centers =',4i5,/)
      end
      subroutine ssdss(qq)
      implicit real*8  (a-h,o-z)
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
      parameter (mxp2 = mxprms*mxprms)
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
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
      dimension qq(*)
      data pt5,pt2/0.5d0,0.2d0/
      data third/0.33333333333333333d0/
      data one/1.0d0/
      data pi252/34.986836655250d0/
      data pie4 /7.85398163397448d-01/
      do 20 i = 1 , 12
         qq(ic7+i) = 0.0d0
 20   continue
c
c     lgmax = ngd
      do 150 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*xk
         akyk = ak*yk
         akzk = ak*zk
         csk = csc(kg)*pi252
c        if (okandl) lgmax = kg
         do 140 lg = 1 , ngd
            al = dg(lg)
            b = ak + al
            b1 = one/b
            bbrrk = al*brrk*b1
            if ((bbrrk+rsmall).le.tol1) then
               exkl = dexp(-bbrrk)*csd(lg)*csk*b1
               xb = (akxk+al*xl)*b1
               yb = (akyk+al*yl)*b1
               zb = (akzk+al*zl)*b1
               nn = 0
               n = 0
               jgmax = ngb
               do 130 ig = 1 , nga
                  ai = ag(ig)
                  if (oiandj) jgmax = ig
                  do 120 jg = 1 , jgmax
                     n = n + 1
                     aj = bg(jg)
                     a = ai + aj
                     if ((bbrrk+r(n)).le.tol2) then
                        ab = a*b
                        aandb = a + b
                        expe = dd(nn+1)*exkl*exij(n)/dsqrt(aandb)
                        if (dabs(expe).ge.tol4) then
                           rho = ab/aandb
                           xa = x1(n)
                           ya = y1(n)
                           za = z1(n)
                           x = rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
                           if (x.gt.5.0d0) then
                              xinv = one/x
                              if (x.le.15.0d0) then
                                 g = dexp(-x)
                                 if (x.gt.10.0d0) then
                                    ww1 = (((-1.8784686463512d-01*xinv+
     +                                 2.2991849164985d-01)
     +                                 *xinv-4.9893752514047d-01)
     +                                 *xinv-2.1916512131607d-05)
     +                                 *g + dsqrt(pie4*xinv)
                                    f1 = (ww1-g)*xinv*pt5
                                    rt1 = f1/(ww1-f1)
                                 else
                                    ww1 =
     +                                 ((((((4.6897511375022d-01*xinv-
     +                                 6.9955602298985d-01)
     +                                 *xinv+5.3689283271887d-01)
     +                                 *xinv-3.2883030418398d-01)
     +                                 *xinv+2.4645596956002d-01)
     +                                 *xinv-4.9984072848436d-01)
     +                                 *xinv-3.1501078774085d-06)
     +                                 *g + dsqrt(pie4*xinv)
                                    f1 = (ww1-g)*xinv*pt5
                                    rt1 = f1/(ww1-f1)
                                 end if
                              else if (x.gt.33.0d0) then
                                 ww1 = dsqrt(pie4*xinv)
                                 rt1 = pt5/(x-pt5)
                              else
                                 g = dexp(-x)
                                 ww1 = ((1.9623264149430d-01*xinv-
     +                                 4.9695241464490d-01)
     +                                 *xinv-6.0156581186481d-05)
     +                                 *g + dsqrt(pie4*xinv)
                                 f1 = (ww1-g)*xinv*pt5
                                 rt1 = f1/(ww1-f1)
                              end if
                           else if (x.gt.1.0d0) then
                              if (x.gt.3.0d0) then
                                 y = x - 4.0d0
                                 f1 = ((((((((((-2.62453564772299d-11*y+
     +                                3.24031041623823d-10)
     +                                *y-3.614965656163d-09)
     +                                *y+3.760256799971d-08)
     +                                *y-3.553558319675d-07)
     +                                *y+3.022556449731d-06)
     +                                *y-2.290098979647d-05)
     +                                *y+1.526537461148d-04)
     +                                *y-8.81947375894379d-04)
     +                                *y+4.33207949514611d-03)
     +                                *y-1.75257821619926d-02)
     +                                *y + 5.28406320615584d-02
                                 ww1 = (x+x)*f1 + dexp(-x)
                                 rt1 = f1/(ww1-f1)
                              else
                                 y = x - 2.0d0
                                 f1 = ((((((((((-1.61702782425558d-10*y+
     +                                1.96215250865776d-09)
     +                                *y-2.14234468198419d-08)
     +                                *y+2.17216556336318d-07)
     +                                *y-1.98850171329371d-06)
     +                                *y+1.62429321438911d-05)
     +                                *y-1.16740298039895d-04)
     +                                *y+7.24888732052332d-04)
     +                                *y-3.79490003707156d-03)
     +                                *y+1.61723488664661d-02)
     +                                *y-5.29428148329736d-02)
     +                                *y + 1.15702180856167d-01
                                 ww1 = (x+x)*f1 + dexp(-x)
                                 rt1 = f1/(ww1-f1)
                              end if
                           else if (x.gt.3.0d-07) then
                              f1 = ((((((((-8.36313918003957d-08*x+
     +                             1.21222603512827d-06)
     +                             *x-1.15662609053481d-05)
     +                             *x+9.25197374512647d-05)
     +                             *x-6.40994113129432d-04)
     +                             *x+3.78787044215009d-03)
     +                             *x-1.85185172458485d-02)
     +                             *x+7.14285713298222d-02)
     +                             *x-1.99999999997023d-01)
     +                             *x + 3.33333333333318d-01
                              ww1 = (x+x)*f1 + dexp(-x)
                              rt1 = f1/(ww1-f1)
                           else
                              rt1 = pt5 - x*pt2
                              ww1 = one - x*third
                           end if
                           u2 = rt1*rho
                           f00 = expe*ww1/(ab+u2*aandb)
                           joff = ic7
                           do 110 iper = 1 , 4
                              if (oskip(iper)) go to 110
                              go to (30,40,50,60) , iper
 30                           alpha = ai + ai
                              xd = xi
                              yd = yi
                              zd = zi
                              go to 70
 40                           alpha = aj + aj
                              xd = xj
                              yd = yj
                              zd = zj
                              go to 70
 50                           alpha = ak + ak
                              xd = xk
                              yd = yk
                              zd = zk
                              go to 70
 60                           alpha = al + al
                              xd = xl
                              yd = yl
                              zd = zl
 70                           bxbd = b*(xb-xd)
                              bybd = b*(yb-yd)
                              bzbd = b*(zb-zd)
                              axad = a*(xa-xd)
                              ayad = a*(ya-yd)
                              azad = a*(za-zd)
                              c1x = bxbd + axad
                              c1y = bybd + ayad
                              c1z = bzbd + azad
                              go to (90,90,80,80) , iper
 80                           c2x = a*bxbd
                              c2y = a*bybd
                              c2z = a*bzbd
                              go to 100
 90                           c2x = b*axad
                              c2y = b*ayad
                              c2z = b*azad
 100                          qq(joff+1) = qq(joff+1) + (u2*c1x+c2x)
     +                           *f00*alpha
                              qq(joff+2) = qq(joff+2) + (u2*c1y+c2y)
     +                           *f00*alpha
                              qq(joff+3) = qq(joff+3) + (u2*c1z+c2z)
     +                           *f00*alpha
                              joff = joff + 3
 110                       continue
                        end if
                     end if
                     nn = nn + 4
 120              continue
 130           continue
            end if
 140     continue
 150  continue
      if (natomd(1).ne.natomd(2)) return
      qq(ic7+1) = qq(ic7+1) + qq(ic7+4)
      qq(ic7+2) = qq(ic7+2) + qq(ic7+5)
      qq(ic7+3) = qq(ic7+3) + qq(ic7+6)
      qq(ic7+4) = qq(ic7+7)
      qq(ic7+5) = qq(ic7+8)
      qq(ic7+6) = qq(ic7+9)
      natomd(2) = natomd(3)
      natomd(3) = natomd(4)
      natomd(4) = 0
      npass = npass - 1
      return
      end
      subroutine subsd(x,y,z,xd,yd,zd,
     *a,m1,m2,m3,m4,i1,i2,k1,k2,ncdim)
      implicit real*8  (a-h,o-z)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      real*8 q4
      integer lit,ljt,lkt,llt,loci,locj,lock,locl
      integer mini,minj,mink,minl,maxi,maxj,maxk,maxl
      integer nij,ijd,kld,ijkld,ncontr
      common/dshlno/q4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ijd,kld,ijkld,ncontr
c
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),
     *         xd(ncdim,*),yd(ncdim,*),zd(ncdim,*)
      parameter (ncmax=65)
      dimension a(ncmax)
c
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
 40               continue
                  n4 = n4 + k2
                  go to 90
 50               do 60 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) - x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) - y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) - z(nr,n4-k2)
 60               continue
                  n4 = n4 + k2
                  go to 90
 70               fac = -dble(l-1)
                  do 80 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) + fac*x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) + fac*y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) + fac*z(nr,n4-k2)
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
      subroutine symde(ict,natoms)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 de
      common/grad2/de(3,maxat)
c
      common/junk/ptr(3,144)
      dimension ict(natoms,*)
      data dzero,done /0.0d0,1.0d0/
c
c     ----- read in tranformation matrices of coordinates. -----
c
c     ----- and
c
c     ----- read in transformation table- atoms versus symmetry operatio
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      nav = lenwrd()
      call readi(ict,nw196(6)*nav,ibl196(6),idaf)
c
c     ----- symmetryze gradient vector -----
c
      do 120 ic = 1 , nat
         do 50 it = 1 , nt
            if (ict(ic,it).gt.ic) go to 120
 50      continue
c
c ps ----start-----
c
c the following code zeros the gradient on atoms with no shells 
c it is suppressed as these gradients are needed for hybrid QM/MM
c calculations
c
c         do 60 ii = 1 , nshell
c            if (katom(ii).eq.ic) go to 80
c 60      continue
c         do 70 it = 1 , nt
c            knu = ict(ic,it)
c            de(1,knu) = 0.0d0
c            de(2,knu) = 0.0d0
c            de(3,knu) = 0.0d0
c 70      continue
c
c  ps ---end----
c
 80      dedx = dzero
         dedy = dzero
         dedz = dzero
         do 90 it = 1 , nt
            icnu = ict(ic,it)
            dedxp = de(1,icnu)
            dedyp = de(2,icnu)
            dedzp = de(3,icnu)
            n = 3*(it-1)
            dedx = dedx + dedxp*ptr(1,n+1) + dedyp*ptr(2,n+1)
     +             + dedzp*ptr(3,n+1)
            dedy = dedy + dedxp*ptr(1,n+2) + dedyp*ptr(2,n+2)
     +             + dedzp*ptr(3,n+2)
            dedz = dedz + dedxp*ptr(1,n+3) + dedyp*ptr(2,n+3)
     +             + dedzp*ptr(3,n+3)
 90      continue
         de(1,ic) = dedx
         de(2,ic) = dedy
         de(3,ic) = dedz
         do 110 it = 1 , nt
            icnu = ict(ic,it)
            if (icnu.ne.ic) then
               if (it.ne.nt) then
                  it1 = it + 1
                  do 100 jt = it1 , nt
                     if (ict(ic,jt).eq.icnu) go to 110
 100              continue
               end if
               jt = invt(it)
               n = 3*(jt-1)
               de(1,icnu) = de(1,ic)*ptr(1,n+1) + de(2,ic)*ptr(2,n+1)
     +                      + de(3,ic)*ptr(3,n+1)
               de(2,icnu) = de(1,ic)*ptr(1,n+2) + de(2,ic)*ptr(2,n+2)
     +                      + de(3,ic)*ptr(3,n+2)
               de(3,icnu) = de(1,ic)*ptr(1,n+3) + de(2,ic)*ptr(2,n+3)
     +                      + de(3,ic)*ptr(3,n+3)
            end if
 110     continue
 120  continue
      dum = done/dble(nt)
      do 140 n = 1 , nat
         do 130 i = 1 , 3
            de(i,n) = de(i,n)*dum
 130     continue
 140  continue
c
      return
      end
      subroutine tpdder(tpdm,lentpd,da,db,ltri,ndens,ii,jj,kk,ll,q4)
      implicit real*8  (a-h,o-z)
c
c     constructs scf tpdm's from cross products of two
c     sets of one particle density matrices. first , da ,
c     may be thought of as the unperturbed dm, second , in
c     db , may be thought of as a set of derivative density
c     matrices.  intended for use in property derivative
c     calculations
c
c     shells (ii,jj,kk,ll) : q4 =degeneracy (symmetry) factor
c
      logical ijeq,kleq
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
      dimension tpdm(lentpd,ndens),da(ltri),db(ltri,ndens)
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
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
      do 60 i = mini , maxi
         nj = ni
         i1 = loci + i
         ii1 = iky(i1)
         do 50 j = minj , maxj
            nk = nj
            j1 = locj + j
            jj1 = iky(j1)
            if (i1.ge.j1) then
               mij = ii1 + j1
            else
               mij = jj1 + i1
            end if
            ijeq = i1.eq.j1
            do 40 k = mink , maxk
               nl = nk
               k1 = lock + k
               kk1 = iky(k1)
               if (i1.ge.k1) then
                  mik = ii1 + k1
               else
                  mik = kk1 + i1
               end if
               if (j1.ge.k1) then
                  mjk = jj1 + k1
               else
                  mjk = kk1 + j1
               end if
               do 30 l = minl , maxl
                  nn = nl
                  l1 = locl + l
                  ll1 = iky(l1)
                  if (k1.ge.l1) then
                     mkl = kk1 + l1
                  else
                     mkl = ll1 + k1
                  end if
                  if (j1.ge.l1) then
                     mjl = jj1 + l1
                  else
                     mjl = ll1 + j1
                  end if
                  if (i1.ge.l1) then
                     mil = ii1 + l1
                  else
                     mil = ll1 + i1
                  end if
                  kleq = k1.eq.l1
                  do 20 ndd = 1 , ndens
                     dfac = (da(mij)*db(mkl,ndd)+da(mkl)*db(mij,ndd))
     +                      *4.0d0 - da(mik)*db(mjl,ndd) - da(mjl)
     +                      *db(mik,ndd) - da(mjk)*db(mil,ndd) - da(mil)
     +                      *db(mjk,ndd)
                     if (ijeq) dfac = dfac*0.5d0
                     if (kleq) dfac = dfac*0.5d0
                     tpdm(nn,ndd) = tpdm(nn,ndd) + dfac*q4
 20               continue
                  nl = nl + inc2
 30            continue
               nk = nk + inc3
 40         continue
            nj = nj + inc4
 50      continue
         ni = ni + inc5
 60   continue
      return
      end
      subroutine ver_drv2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/drv2e.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
