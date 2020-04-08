      subroutine mcx(xscm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       extension to hondo, drf + bem method to include limited
c       monte carlo sampling of discrete (polarizable) surroundings
c
c       written by alex de vries and piet van duijnen
c       dept. of chemistry
c       state university of groningen
c       nijenborgh 4
c       9747 ag groningen
c       the netherlands
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  this routine is the driver routine for the monte carlo
c       sampling of the discrete (polarizable) surroundings
c       it consists of the following calls:
c
c       mcset: sets up monte carlo run
c       mcrun: performs the actual monte carlo run
c       mcout: writes monte carlo output
c
c-----  declarations
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
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
      character*8 scftyp
      common /scfopt2/ scftyp
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      character*8 wfntyp
      common /wfnopt/ wfntyp
c
      character*8 runtyp
      common /runopt/ runtyp
c
      character*8 hndtyp
      common /hndopt/ hndtyp
c
c
cmw      include 'drf/dimpar'
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 gmdrf
      integer iguess
      common /mcinp/ gmdrf(mxst),iguess(mxst)
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
      dimension xscm(*)
c
      data zero /0.0d00/
c
caleko
c
c  Monte Carlo input has been read in by rfin
c  together with the rest of the drf-input.
c
c  collection of energies in the MC loop should 
c  only be done if wave-function converged;
c  this requires the logical excess to be set 
c  .true. in all scf routines
c  this is NOT YET IMPLEMENTED!
c
      excess = .false.
c
c  Now we proceed to the real work.
c
      l2 = num*(num+1)/2
      l3 = num*num
c
c-----  set the input vectors, etc... according to iguess
c
      do 1200, ist = 1, imcst
c 1-----
        if (iguess(ist) .eq. 0) then
c   2-----
          i10 = igmem_alloc(4*l2+2*l3)
          call daread(idafh,ioda,xscm(i10),l2,12)
          call daread(idafh,ioda,xscm(i10+l2),l2,13)
          call daread(idafh,ioda,xscm(i10+2*l2),l3,15)
          call daread(idafh,ioda,xscm(i10+2*l2+l3),l2,16)
c
          if (scftyp .eq. 'uhf') then
            call daread(idafh,ioda,xscm(i10+3*l2+l3),l3,19)
            call daread(idafh,ioda,xscm(i10+3*l2+2*l3),l2,20)
          endif
c
          call dawrit(idafdrf,iodadrf,xscm(i10),
     *      l2,181+ist,navdrf)
          call dawrit(idafdrf,iodadrf,xscm(i10+l2),
     *      l2,191+ist,navdrf)
          call dawrit(idafdrf,iodadrf,xscm(i10+2*l2),
     *      l3,201+ist,navdrf)
          call dawrit(idafdrf,iodadrf,xscm(i10+2*l2+l3),
     *      l2,211+ist,navdrf)
c
          if (scftyp .eq. 'uhf') then
            call dawrit(idafdrf,iodadrf,xscm(i10+3*l2+l3),l3,
     2                  221+ist,navdrf)
            call dawrit(idafdrf,iodadrf,xscm(i10+3*l2+2*l3),l2,
     2                  231+ist,navdrf)
          endif
c
          call gmem_free(i10)
c   2-----
        endif
c
        if (gmdrf(ist) .lt. zero) then
          gammc(ist) = gamdrf
        else
          gammc(ist) = gmdrf(ist)
        endif
c 1-----
 1200 continue
c
      call mcset(xscm)
c
      call mcrunr(xscm)
c
      call mcout
c
      return
      end
      subroutine exclset(excld)
c------
c      sets exclusion distance from solvent data
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      character*32 solnams
      dimension solnams(50)
      dimension srad(50)
c
      character *8 errmsg(3)
c
      data solnams /'water','formamide','ethylene glycol','methanol',
     1 'n-methylformamide','2-methoxyethanol','n-methylacetamide',
     2 'ethanol','acetic acid','propanol','isopropanol','nitromethane',
     3 'acetonitril','dimethyl sulphoxide','dmso','t-butanol',
     4 'dimethyl formamide','dimethyl acetamide','acetone',
     5 'nitrobenzene','dichloromethane','pyridine','acetophenon',
     6 'chloroform','1,2-dimethoxyethane','ethyl acetate',
     7 'tetrahydrofuran','1,4-dioxan','diethyl ether','benzene',
     8 'carbon disulphide','carbon tetrachloride','cyclohexane',
     9 'n-hexane','hexane','x','vacuum',13*' '/
c
c-----  solvent radii in angstrom from density
c
c       srad = ((3*mass*10-6)/(4*pi*rho*la))**1/3
c
c       where mass = molecular mass in amu
c             rho  = solvent density in g/cm
c             la   = avogadro's number
c
      data srad /
     + 1.927d0,  0.0d0,  0.0d0,2.522d0,  0.0d0,  0.0d0,  0.0d0,
     + 2.850d0,  0.0d0,3.126d0,  0.0d0,2.773d0,2.746d0,3.041d0,
     + 3.041d0,  0.0d0,  0.0d0,  0.0d0,3.076d0,3.435d0,  0.0d0,
     +   0.0d0,3.591d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,
     +   0.0d0,3.279d0,  0.0d0,3.369d0,  0.0d0,3.728d0,3.728d0,
     +  15*0.0d0/
c
      data errmsg /'program','stop in','-exclset'/
c
      if (solnam .eq. 'x') then
        write(iwr,1001) solnam
 1001   format(/,' solvent name ',a32,' incompatible with excld < 0.0:',
     1  /,' give exclusion distance in $montec',
     2    ', or explicit solvent name in $bem')
        call hnderr(3,errmsg)
        return
      endif
c
c-----  loop over standard solvents
c
      do 100, i = 1, 35
c 1-----
        if (solnam .eq. solnams(i)) then
c   2-----
          if (srad(i) .eq. 0.0) then
c     3-----
            write(iwr,1002) solnam
 1002       format(/,' no exclusion distance implemented for',a32,/,
     1  ' feel free to get it and insert it ',
     2  'in srad in subroutine exclset')
            call hnderr(3,errmsg)
            return
c     3-----
          endif
c
          excld = 2*srad(i)
          return
c   2-----
        endif
c 1-----  next solvent name
  100 continue
c
      write(iwr,1003) solnam
 1003 format(/,' solvent name ',a32,' not found:',/,
     1 ' give explicit data in $bem and $montec for this solvent')
c
      return
      end
      subroutine mcset(xscm)
c------
c      subroutine prepares monte carlo calculation
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
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
      integer ia
      common /ijpair/ ia(3*mxpts)
c
      dimension xscm(*)
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      character*8 wfntyp
      common /wfnopt/ wfntyp
c
      character*8 runtyp
      common /runopt/ runtyp
c
      character*8 hndtyp
      common /hndopt/ hndtyp
c
c
c
      real*8 enci, ec, eci
      integer nstat
      common /enrgci/ enci,ec,eci(10),nstat
c
      integer istate
      common /enrgc2/ istate(10)
c
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      integer neq2eps, momsav
      common /neqpar1/ neq2eps,momsav
c
      real*8 kapneq1, kapneq2, epsneq1, epsneq2
      common /neqpar2/ kapneq1,kapneq2,epsneq1,epsneq2
c
c
       integer idafinf, navinf, iodainf
       common /infdaf/ idafinf,navinf,iodainf(2,mxneq+1)
c
       integer idafsta, navsta, iodasta
       common /stadaf/ idafsta,navsta,iodasta(2,mxneq+1)
c
       integer idafcst, navcst, iodacst
       common /cstdaf/ idafcst,navcst,iodacst(2,mxneq+1)
c
       integer idafind, navind, iodaind
       common /inddaf/ idafind,navind,iodaind(2,mxneq+1)
c
       integer idafino, navino, iodaino
       common /inodaf/ idafino,navino,iodaino(2,mxneq+1)
c
       integer idafpol, navpol, iodapol
       common /poldaf/ idafpol,navpol,iodapol(2,mxneq+1)
c
       integer idafdis, navdis, iodadis
       common /disdaf/ idafdis,navdis,iodadis(2,mxneq+1)
c
       integer idafrep, navrep, iodarep
       common /repdaf/ idafrep,navrep,iodarep(2,mxneq+1)
c
       integer idafrqm, navrqm, iodarqm
       common /rqmdaf/ idafrqm,navrqm,iodarqm(2,mxneq+1)
c
c
      character*80 textrf
      common /texneq/ textrf
c
      character*80 textneq
      common /neqtex/ textneq(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
      real*8 realneq
      common /bufb/ realneq(mxneq*10)
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 ulow, geoml
      common /mcgeom/ ulow,geoml(3,mxpts)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
c
      dimension udiff(mxst*(mxst+1)/2)
c
      character*40 header, header2
      character*16 estring1, estring2, estring3
      character*8 endwrd
c
      dimension dumm(mxgran), dumx(mxgrpar)
c
c-----  declarations
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
      data boltz /3.16682973770242d-06/
      data incf /51/
c
      data header /'start configuration data '/
      data header2 /'for the second state '/
      data estring1, estring2, estring3 /
     1 'total energy  : ', 'class energy  : ', 'qm energy     : '/
      data endwrd /' $end'/
c
c-----  begin
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c-----  calculate -1/kt (in a.u.)
c
      onekt = - one / (boltz*temp)
c
      mcupdt = .false.
      ixnsurf = 0
c
      ist = iactst
      ngrpair = ngran*(ngran+1)/2
      gamdrf = gammc(ist)
c
c-----  read wavefunction data for actual state
c
      i10 = igmem_alloc(l3)
      call daread(idafdrf,iodadrf,xscm(i10),l2,181+ist)
      call dawrit(idafh,ioda,xscm(i10),l2,12,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l2,191+ist)
      call dawrit(idafh,ioda,xscm(i10),l2,13,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l3,201+ist)
      call dawrit(idafh,ioda,xscm(i10),l3,15,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l2,211+ist)
      call dawrit(idafh,ioda,xscm(i10),l2,16,navh)
c
      if (scftyp .eq. 'uhf') then
c 1-----
        call daread(idafdrf,iodadrf,xscm(i10),l3,221+ist)
        call dawrit(idafh,ioda,xscm(i10),l3,19,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,231+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,20,navh)
c 1-----
      endif
      call gmem_free(i10)
c
c-----  perform calculation of reference energy (input geometry)
c
cnot      call energx
      call hfscf(xscm)
c
c-----  perform analysis of reference energy
c
      ieps = 0
      call arfanal(ieps,unucrep,
     1  eoneel(ist),ekin(ist),enua(ist),etwoel(ist),
     1 uqm(ist),uscf(ist),snucnuc(ist),selel(ist),snua(ist),stwoel(ist),
     2 smolnuc(ist),smolel(ist),snucmol(ist),selmol(ist),smolmol(ist),
     3 suqm(ist),upolqm(ist),uneqnuc(ist),uneqel(ist),uneqqm(ist),
     4 ustanuc(ist),ustael(ist),ustaqm(ist),uclase,uclasd,uclasr,uclas,
     5 suclas,upolcl,uneqcl,ustacl,extnuc(ist),extel(ist),sextnuc(ist),
     6 sextel(ist),sextmol(ist),selext(ist),snucext(ist),smolext(ist),
     7 stotnuc(ist),stotel(ist),stotmol(ist),stotext(ist),stabtot(ist),
     8 uelst(ist),suint(ist),uint(ist),udisp(ist),rdisp(ist),
     9 repmod(ist),upoleq(ist),ucstst(ist),ucstpl(ist),uneq(ist),
     1 usta(ist),upolneq(ist),ustaneq(ist),uens(ist),uclasg,uclaseg,
     2 uclasdg,uclasrg,suclasg,upolclg,uneqclg,ustaclg,extnucg(1,ist),
     3 extelg(1,ist),sxtnucg(1,ist),sextelg(1,ist),sxtmolg(1,ist),
     4 selextg(1,ist),snucxtg(1,ist),smolxtg(1,ist),stotxtg(1,ist),
     5 uelstg(1,ist),suintg(1,ist),uintg(1,ist),repmodg(1,ist),xscm)
c
      if (itwoeps .eq. 1) then
c 1-----
        ieps = 1
        call arfanal(ieps,unucrepo,
     1 eoneelo(ist),ekino(ist),enuao(ist),etwoelo(ist),
     1 uqm(ist),uscfo(ist),snucno(ist),selel(ist),snuao(ist),
     1 stwoelo(ist),
     2 smolno(ist),smolelo(ist),snucmo(ist),selmolo(ist),smolmo(ist),
     3 suqmo(ist),upolqmo(ist),uneqno(ist),uneqelo(ist),uneqqmo(ist),
     4 ustano(ist),ustaelo(ist),ustaqmo(ist),uclaseo,uclasdo,
     4 uclasro,uclaso,
     5 suclaso,upolclo,uneqclo,ustaclo,extnuco(ist),extelo(ist),
     5 sextno(ist),
     6 sextelo(ist),sextmo(ist),selexto(ist),snucexo(ist),smolexo(ist),
     7 stotno(ist),stotelo(ist),stotmo(ist),stotexo(ist),stabto(ist),
     8 uelsto(ist),suinto(ist),uinto(ist),udisp(ist),rdispo(ist),
     9 repmo(ist),upoleqo(ist),ucststo(ist),ucstplo(ist),uneqo(ist),
     1 ustao(ist),upolno(ist),ustano(ist),uens(ist),uclasgo,uclasego,
     2 uclasdgo,uclasrgo,suclsog,uplclog,uneqclgo,ustaclgo,
     2 extnucgo(1,ist),
     3 extelgo(1,ist),sextnog(1,ist),sxtelog(1,ist),sextmog(1,ist),
     4 selxtog(1,ist),sncexog(1,ist),smlexog(1,ist),sttexog(1,ist),
     5 uelstgo(1,ist),suintog(1,ist),uintgo(1,ist),repmodgo(1,ist),
     5 xscm)
c 1-----
      endif
c
      call enanal
c
c-----  store reference energies uref
c
      ureft(ist) = uens(ist)
      uold(ist) = ureft(ist)
c
      eqmol(ist) = eqm(ist)
      eclasol = eclas
      ecleol = eclase
      ecldol = eclasd
      eclrol = eclasr
      eintol(ist) = eint(ist)
      eelstol(ist) = eelst(ist)
      edisol(ist) = edisp(ist)
      extnol(ist) = extnuc(ist)
      erepol(ist) = erep(ist)
      epolol(ist) = epol(ist)
      eneqol(ist) = eneq(ist)
      eneqpol(ist) = eneqp(ist)
c
      do 100, igr = 1, ngran
c 1-----
        do 50, jgr = 1, igr
c   2-----
          indxcl = ia(igr) + jgr
          eclsolg(indxcl) = eclasg(indxcl)
          ecleolg(indxcl) = eclaseg(indxcl)
          ecldolg(indxcl) = eclasdg(indxcl)
          eclrolg(indxcl) = eclasrg(indxcl)
c   2-----
   50   continue
c
        eintolg(igr,ist) = eintg(igr,ist)
        elstolg(igr,ist) = eelstg(igr,ist)
        erepolg(igr,ist) = erepg(igr,ist)
        extnolg(igr,ist) = extnucg(igr,ist)
c
c       NOTE: edispg array set to zero
c       At present, the QM/Classical dispersion 
c       estimate cannot be analysed in terms of 
c       classical group contributions 
c       (in fact, it may be done only approximately so
c       and is NOT IMPLEMENTED)
c
        edispg(igr,ist) = zero
c
c 1-----
  100 continue
c
      call drfout
c
c-----  store reference geometry as lowest energy conformation
c
      ulow = ureft(ist)
      do 200, ip = 1, nxtpts
        do 190, k = 1, 3
          geoml(k,ip) = xpts(k,ip)
  190   continue
  200 continue
c
c-----  copy wave function data
c
      i10 = igmem_alloc(l3)
      call daread(idafh,ioda,xscm(i10),l2,12)
      call dawrit(idafdrf,iodadrf,xscm(i10),l2,181+ist,navdrf)
c
      call daread(idafh,ioda,xscm(i10),l2,13)
      call dawrit(idafdrf,iodadrf,xscm(i10),l2,191+ist,navdrf)
c
      call daread(idafh,ioda,xscm(i10),l3,15)
      call dawrit(idafdrf,iodadrf,xscm(i10),l3,201+ist,navdrf)
c
      call daread(idafh,ioda,xscm(i10),l2,16)
      call dawrit(idafdrf,iodadrf,xscm(i10),l2,211+ist,navdrf)
c
      if (scftyp .eq. 'uhf') then
c 1-----
        call daread(idafh,ioda,xscm(i10),l3,19)
        call dawrit(idafdrf,iodadrf,xscm(i10),l3,221+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,20)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,231+ist,navdrf)
c 1-----
      endif
      call gmem_free(i10)
c
c-----  initialise collected potential at expansion centra if
c       required
c
c     if (isolsav .eq. 1) then
c 1-----
c       call clear(xscm,3*nqdim)
c
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,151,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,152,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,153,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,154,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,155,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqdim,156,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqcls,161,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqcls,162,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqcls,163,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,nqcls,164,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,1,165,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,1,166,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,1,167,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,1,168,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,3,169,navdrf)
c       call dawrit(idafdrf,iodadrf,xscm,3,170,navdrf)
c
c  -----  copy appropriate data
c
c       call daread(idafsta,iodasta,xscm,nqdim,index)
c       call dawrit(idafsta,iodasta,xscm,nqdim,index+1,navsta)
c
c       call daread(idafdis,iodadis,xscm,nqcls,index)
c       call dawrit(idafdis,iodadis,xscm,nqcls,index+1,navdis)
c
c       call daread(idafrep,iodarep,xscm,nqcls,index)
c       call dawrit(idafrep,iodarep,xscm,nqcls,index+1,navrep)
c
c       call daread(idafrqm,iodarqm,xscm,1,index)
c       call dawrit(idafrqm,iodarqm,xscm,1,index+1,navrqm)
c
c       call daread(idafcst,iodacst,xscm,3,index)
c       call dawrit(idafcst,iodacst,xscm,3,index+1,navcst)
c
c       if (field(5:) .ne. ' ') then
c   2-----
c         call daread(idafind,iodaind,xscm,nqdim,index)
c         call dawrit(idafind,iodaind,xscm,nqdim,index+1,navind)
c
c         if (itwoeps .eq. 1) then
c     3-----
c           call daread(idafino,iodaino,xscm,nqdim,index)
c           call dawrit(idafino,iodaino,xscm,nqdim,index+1,navino)
c     3-----
c         endif
c
c         call daread(idafpol,iodapol,xscm,7,index)
c         call dawrit(idafpol,iodapol,xscm,7,index+1,navpol)
c   2-----
c       endif
c
c       call addpot(nqdim,nqcls,itwoeps,field,
c    2              xscm(1),xscm(1+nqdim),xscm(1+2*nqdim))
c 1-----
c     endif
c
c-----  open tape51 and write reference data (if required)
c
      if (imcout .ge. 3) then
c 1-----
        open(unit=incf, file='mcdata', access='sequential',
     1       form='formatted',status='unknown')
c
cnot        call rewfil(incf)
        write(incf,1001) outfor
 1001   format(a16)
c
        write(incf,1002) header
 1002   format(2x,a40)
c
        write(incf,1003) estring1, ureft(ist), estring2, eclas
 1003   format(2x,a16,2x,e20.14,2x,a16,2x,e20.14)
c
c  -----  reference energy data
c
        if (outfor(5:7) .eq. 'all') then
          write(incf,1004) eqm(ist), eclase, eclasd, eclasr
          write(incf,1004) eint(ist), eelst(ist), edisp(ist), erep(ist)
          write(incf,1004) epol(ist), eneq(ist), eneqp(ist)
 1004     format(4e20.14)
        endif
c
        if (outfor(9:11) .eq. 'all') then
c   2-----
          do 300, k = 1, ngrpair
            write(incf,1004)
     1      eclasg(k), eclaseg(k), eclasdg(k), eclasrg(k)
  300     continue
c
          do 320, igr = 1, ngran
            write(incf,1004)
     1 eintg(igr,ist), eelstg(igr,ist), edispg(igr,ist), erepg(igr,ist)
  320     continue
c   2-----
        endif
c
c  -----  counter for group polarisabilities
c
        npols = 0
c
        do 400, ipart = 1, nxtpts
c   2-----
          write(incf,1005) nxcent(ipart),
     1         (xpts(k,ipart), k = 1, 3), chrg(ipart)
 1005     format(a16,4f16.6)
c
          if (nxcent(ipart)(1:5) .eq. 'group') then
c     3-----
            npols = npols + 1
            write(incf,1006) (grpol(k,npols), k = 1, 6)
 1006       format(6(f12.6))
c     3-----
          endif
c   2-----
  400   continue
c
        write(incf,1007) endwrd
 1007   format(a8)
c
        igrp = 0
        do 420, i = 1, npols, 10
          write(incf,1008) (igrpst(igrp+k), ngrppt(igrp+k), k = 1, 10)
 1008     format(20i4)
          igrp = igrp + 10
  420   continue
c 1-----
      endif
c
c-----  calculate further state reference energies
c
      do 500, ist = 2, imcst
c 1-----
c  -----  read wavefunction data for actual state
c
        i10 = igmem_alloc(l3)
        call daread(idafdrf,iodadrf,xscm(i10),l2,181+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,12,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,191+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,13,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l3,201+ist)
        call dawrit(idafh,ioda,xscm(i10),l3,15,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,211+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,16,navh)
c
        if (scftyp .eq. 'uhf') then
c   2-----
          call daread(idafdrf,iodadrf,xscm(i10),l3,221+ist)
          call dawrit(idafh,ioda,xscm(i10),l3,19,navh)
c
          call daread(idafdrf,iodadrf,xscm(i10),l2,231+ist)
          call dawrit(idafh,ioda,xscm(i10),l2,20,navh)
c   2-----
        endif
c
        secstat = .true.
        iactst = ist
        gamdrf = gammc(ist)
c
        call gmem_free(i10)
c
c  -----  calculate energy of second state
c
cnot        call energx
        call hfscf(xscm)
c
c  -----  perform analysis of reference energy
c
        ieps = 0
        call arfanal(ieps,unucrep,
     1  eoneel(ist),ekin(ist),enua(ist),etwoel(ist),
     1 uqm(ist),uscf(ist),snucnuc(ist),selel(ist),snua(ist),stwoel(ist),
     2 smolnuc(ist),smolel(ist),snucmol(ist),selmol(ist),smolmol(ist),
     3 suqm(ist),upolqm(ist),uneqnuc(ist),uneqel(ist),uneqqm(ist),
     4 ustanuc(ist),ustael(ist),ustaqm(ist),uclase,uclasd,uclasr,uclas,
     5 suclas,upolcl,uneqcl,ustacl,extnuc(ist),extel(ist),sextnuc(ist),
     6 sextel(ist),sextmol(ist),selext(ist),snucext(ist),smolext(ist),
     7 stotnuc(ist),stotel(ist),stotmol(ist),stotext(ist),stabtot(ist),
     8 uelst(ist),suint(ist),uint(ist),udisp(ist),rdisp(ist),
     9 repmod(ist),upoleq(ist),ucstst(ist),ucstpl(ist),uneq(ist),
     1 usta(ist),upolneq(ist),ustaneq(ist),uens(ist),uclasg,uclaseg,
     2 uclasdg,uclasrg,suclasg,upolclg,uneqclg,ustaclg,extnucg(1,ist),
     3 extelg(1,ist),sxtnucg(1,ist),sextelg(1,ist),sxtmolg(1,ist),
     4 selextg(1,ist),snucxtg(1,ist),smolxtg(1,ist),stotxtg(1,ist),
     5 uelstg(1,ist),suintg(1,ist),uintg(1,ist),repmodg(1,ist),xscm)
c
        if (itwoeps .eq. 1) then
          ieps = 1
        call arfanal(ieps,unucrepo,
     1 eoneelo(ist),ekino(ist),enuao(ist),etwoelo(ist),
     1 uqm(ist),uscfo(ist),snucno(ist),selel(ist),snuao(ist),
     1 stwoelo(ist),
     2 smolno(ist),smolelo(ist),snucmo(ist),selmolo(ist),smolmo(ist),
     3 suqmo(ist),upolqmo(ist),uneqno(ist),uneqelo(ist),uneqqmo(ist),
     4 ustano(ist),ustaelo(ist),ustaqmo(ist),uclaseo,uclasdo,
     4 uclasro,uclaso,
     5 suclaso,upolclo,uneqclo,ustaclo,extnuco(ist),extelo(ist),
     5 sextno(ist),
     6 sextelo(ist),sextmo(ist),selexto(ist),snucexo(ist),smolexo(ist),
     7 stotno(ist),stotelo(ist),stotmo(ist),stotexo(ist),stabto(ist),
     8 uelsto(ist),suinto(ist),uinto(ist),udisp(ist),rdispo(ist),
     9 repmo(ist),upoleqo(ist),ucststo(ist),ucstplo(ist),uneqo(ist),
     1 ustao(ist),upolno(ist),ustano(ist),uens(ist),uclasgo,uclasego,
     2 uclasdgo,uclasrgo,suclsog,uplclog,uneqclgo,ustaclgo,
     2 extnucgo(1,ist),
     3 extelgo(1,ist),sextnog(1,ist),sxtelog(1,ist),sextmog(1,ist),
     4 selxtog(1,ist),sncexog(1,ist),smlexog(1,ist),sttexog(1,ist),
     5 uelstgo(1,ist),suintog(1,ist),uintgo(1,ist),repmodgo(1,ist),
     5 xscm)
c
        endif
c
        call enanal
c
c  -----  store reference energies uref
c
        ureft(ist) = uens(ist)
        uold(ist) = ureft(ist)
        do 490, kst = 1, ist
          indxdf = ia(ist) + kst
          udiff(indxdf) = uens(ist) - uens(kst)
          edifol(indxdf) = udiff(indxdf)
  490   continue
c
        eqmol(ist) = eqm(ist)
        eintol(ist) = eint(ist)
        eelstol(ist) = eelst(ist)
        edisol(ist) = edisp(ist)
        extnol(ist) = extnuc(ist)
        erepol(ist) = erep(ist)
        epolol(ist) = epol(ist)
        eneqol(ist) = eneq(ist)
        eneqpol(ist) = eneqp(ist)
c
        do 470, igr = 1, ngran
c   2-----
          eintolg(igr,ist) = eintg(igr,ist)
          elstolg(igr,ist) = eelstg(igr,ist)
          erepolg(igr,ist) = erepg(igr,ist)
          extnolg(igr,ist) = extnucg(igr,ist)
c   2-----
  470 continue
c
c  -----  write wavefunction data for actual state
c
        i10 = igmem_alloc(i10)
        call daread(idafh,ioda,xscm(i10),l2,12)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,181+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,13)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,191+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l3,15)
        call dawrit(idafdrf,iodadrf,xscm(i10),l3,201+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,16)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,211+ist,navdrf)
c
        if (scftyp .eq. 'uhf') then
c   2-----
          call daread(idafh,ioda,xscm(i10),l3,19)
          call dawrit(idafdrf,iodadrf,xscm(i10),l3,221+ist,navdrf)
c
          call daread(idafh,ioda,xscm(i10),l2,20)
          call dawrit(idafdrf,iodadrf,xscm(i10),l2,231+ist,navdrf)
c   2-----
        endif
c
        call gmem_free(i10)
c
        if (imcout .ge. 3) then
c   2-----
          write(incf,1002) header2
c
          write(incf,1003) estring1, ureft(ist), estring3, eqm(ist)
c
c    -----  reference energy data
c
          if (outfor(5:7) .eq. 'all') then
            write(incf,1004) eint(ist),eelst(ist),edisp(ist),erep(ist)
            write(incf,1004) epol(ist), eneq(ist), eneqp(ist)
          endif
c
          if (outfor(9:11) .eq. 'all') then
c     3-----
            do 450, igr = 1, ngran
              write(incf,1004)
     1        eintg(igr,ist), eelstg(igr,ist), erepg(igr,ist)
  450       continue
c     3-----
          endif
c   2-----
        endif
c
        call drfout
c 1-----
  500 continue
c
      secstat = .false.
      iactst = 1
c
c-----  set update flag
c
      mcupdt = .true.
c
c-----  initialize random generator
c
      if (iseed .le. 0) iseed = 1234567
      call ranset(iseed)
c
      return
      end
      subroutine enanal
c------
c      performs analysis of energy contributions, for the benefit
c      of the monte carlo calculation
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
c-----  begin
c
      eqm(iactst) = uqm(iactst)
c
      eclas = uclas + suclas + upolcl + uneqcl
      eclase = uclase
      eclasd = uclasd
      eclasr = uclasr
c
      eelst(iactst) = uelst(iactst) + suint(iactst)
      edisp(iactst) = udisp(iactst) + rdisp(iactst)
      erep(iactst) = repmod(iactst)
c
      eint(iactst) = uint(iactst) + suint(iactst)
     1             + suqm(iactst) + edisp(iactst)
c
      if (field(5:) .ne. ' ') epol(iactst) = upoleq(iactst)
c
      eneq(iactst) = uneq(iactst)
      eneqp(iactst) = upolneq(iactst)
c
      do 100, igr = 1, ngran
c 1-----
        if (iactst .eq. 1) then
c   2-----
          do 200, jgr = 1, igr
c     3-----
            indxcl = ia(igr) + jgr
            eclaseg(indxcl) = uclaseg(indxcl)
            eclasdg(indxcl) = uclasdg(indxcl)
            eclasrg(indxcl) = uclasrg(indxcl)
c
            eclasg(indxcl) = uclasg(indxcl)
c
            if (igr .eq. jgr) then
              eclasg(indxcl) = eclasg(indxcl)
     1        + suclasg(igr,jgr) + upolclg(igr) + uneqclg(igr)
            else
              eclasg(indxcl) = eclasg(indxcl)
     1        + suclasg(igr,jgr) + suclasg(jgr,igr)
            endif
c     3-----
  200    continue
c   2-----
        endif
c
        eelstg(igr,iactst) = uelstg(igr,iactst) + suintg(igr,iactst)
        erepg(igr,iactst) = repmodg(igr,iactst)
        eintg(igr,iactst) = uintg(igr,iactst) + suintg(igr,iactst)
     1                    + edispg(igr,iactst)
c 1-----
  100 continue
c
      return
      end
      subroutine mcrunr(xscm)
c------
c      performs the actual monte carlo run
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
cahvINCLUDE(comdrf/scm)
      dimension xscm(*)
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
       integer idafinf, navinf, iodainf
       common /infdaf/ idafinf,navinf,iodainf(2,mxneq+1)
c
       integer idafsta, navsta, iodasta
       common /stadaf/ idafsta,navsta,iodasta(2,mxneq+1)
c
       integer idafcst, navcst, iodacst
       common /cstdaf/ idafcst,navcst,iodacst(2,mxneq+1)
c
       integer idafind, navind, iodaind
       common /inddaf/ idafind,navind,iodaind(2,mxneq+1)
c
       integer idafino, navino, iodaino
       common /inodaf/ idafino,navino,iodaino(2,mxneq+1)
c
       integer idafpol, navpol, iodapol
       common /poldaf/ idafpol,navpol,iodapol(2,mxneq+1)
c
       integer idafdis, navdis, iodadis
       common /disdaf/ idafdis,navdis,iodadis(2,mxneq+1)
c
       integer idafrep, navrep, iodarep
       common /repdaf/ idafrep,navrep,iodarep(2,mxneq+1)
c
       integer idafrqm, navrqm, iodarqm
       common /rqmdaf/ idafrqm,navrqm,iodarqm(2,mxneq+1)
c
c
      character*80 textrf
      common /texneq/ textrf
c
      character*80 textneq
      common /neqtex/ textneq(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
      real*8 realneq
      common /bufb/ realneq(mxneq*10)
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
      real*8 ulow, geoml
      common /mcgeom/ ulow,geoml(3,mxpts)
c
caleko
      logical obeen, obeen2, obeen3,obeen4
      common/nottwi/obeen,obeen2,obeen3,obeen4
caleko
c
c-----  local variables
c
      dimension udiff(mxst*(mxst+1)/2)
c
      logical accept
c
      character*40 string1,string2,string3,string4
c
      character*8 astring
      character*16 header, estring1, estring2
      character*1 comma
c
      dimension urunt(mxst),uruntq(mxst),urunti(mxst),urunte(mxst),
     1 uruntd(mxst),uruntr(mxst),uruntp(mxst),uruntn(mxst),urtnp(mxst),
     2 urtdf(mxst*(mxst+1)/2)
      dimension urunt2(mxst),uruntq2(mxst),urunti2(mxst),urunte2(mxst),
     1 uruntd2(mxst),uruntr2(mxst),uruntp2(mxst),uruntn2(mxst),
     2 urtnp2(mxst),urtdf2(mxst*(mxst+1)/2)
c
      dimension urtcg(mxgrpar),urtceg(mxgrpar),urtcdg(mxgrpar),
     1 urtcrg(mxgrpar)
      dimension urtcg2(mxgrpar),urtceg2(mxgrpar),urtcdg2(mxgrpar),
     1 urtcrg2(mxgrpar)
c
      dimension urtig(mxgran,mxst),urteg(mxgran,mxst),
     1 urtdg(mxgran,mxst),urtrg(mxgran,mxst)
      dimension urtig2(mxgran,mxst),urteg2(mxgran,mxst),
     1 urtdg2(mxgran,mxst),urtrg2(mxgran,mxst)
c
      dimension usamt(mxst),usamtq(mxst),usamti(mxst),usamte(mxst),
     1 usamtd(mxst),usamtr(mxst),usamtp(mxst),usamtn(mxst),
     2 usmtnp(mxst)
      dimension usmtdf(mxst*(mxst+1)/2),usmtdf2(mxst*(mxst+1)/2)
      dimension usamt2(mxst),usamtq2(mxst),usamti2(mxst),usamte2(mxst),
     1 usamtd2(mxst),usamtr2(mxst),usamtp2(mxst),usamtn2(mxst),
     2 usmtnp2(mxst)
c
      dimension usmtcg(mxgrpar),usmceg(mxgrpar),
     1 usmcdg(mxgrpar), usmcrg(mxgrpar)
      dimension usmtcg2(mxgrpar),usmceg2(mxgrpar),
     1 usmcdg2(mxgrpar), usmcrg2(mxgrpar)
c
      dimension usmtig(mxgran,mxst),usmteg(mxgran,mxst),
     1 usmtdg(mxgran,mxst),usmtrg(mxgran,mxst)
      dimension usmtig2(mxgran,mxst),usmteg2(mxgran,mxst),
     1 usmtdg2(mxgran,mxst),usmtrg2(mxgran,mxst)
c
      dimension ublt(mxst),ubltq(mxst),ublti(mxst),ublte(mxst),
     1 ubltd(mxst),ubltr(mxst),ubltp(mxst),ubltn(mxst),ubltnp(mxst),
     2 ubltdf(mxst*(mxst+1)/2)
      dimension ublt2(mxst),ubltq2(mxst),ublti2(mxst),ublte2(mxst),
     1 ubltd2(mxst),ubltr2(mxst),ubltp2(mxst),ubltn2(mxst),
     2 ubltnp2(mxst),ubltdf2(mxst*(mxst+1)/2)
c
      dimension ubltcg(mxgrpar),ublceg(mxgrpar),
     1 ublcdg(mxgrpar),ublcrg(mxgrpar)
      dimension ubltcg2(mxgrpar),ublceg2(mxgrpar),
     1 ublcdg2(mxgrpar),ublcrg2(mxgrpar)
c
      dimension ubltig(mxgran,mxst),ublteg(mxgran,mxst),
     1 ubltdg(mxgran,mxst),ubltrg(mxgran,mxst)
      dimension ubltig2(mxgran,mxst),ublteg2(mxgran,mxst),
     1 ubltdg2(mxgran,mxst),ubltrg2(mxgran,mxst)
c
      dimension utbl(mxst),uqbl(mxst),uibl(mxst),uebl(mxst),
     1 udbl(mxst), urbl(mxst),upbl(mxst),unbl(mxst),unpbl(mxst),
     2 udfbl(mxst*(mxst+1)/2)
      dimension utbl2(mxst),uqbl2(mxst),uibl2(mxst),uebl2(mxst),
     1 udbl2(mxst),urbl2(mxst),upbl2(mxst),unbl2(mxst),unpbl2(mxst),
     2 udfbl2(mxst*(mxst+1)/2)
c
      dimension ucgbl(mxgrpar),ucegbl(mxgrpar),
     1 ucdgbl(mxgrpar),ucrgbl(mxgrpar)
      dimension ucgbl2(mxgrpar),ucegbl2(mxgrpar),
     1 ucdgbl2(mxgrpar),ucrgbl2(mxgrpar)
c
      dimension uigbl(mxgran,mxst),uegbl(mxgran,mxst),
     1 udgbl(mxgran,mxst),urgbl(mxgran,mxst)
      dimension uigbl2(mxgran,mxst),uegbl2(mxgran,mxst),
     1 udgbl2(mxgran,mxst),urgbl2(mxgran,mxst)
c
      dimension umove(mxst),umqm(mxst),umint(mxst),umelst(mxst),
     1 umdis(mxst),umrep(mxst),umpol(mxst),umneq(mxst),umneqp(mxst)
      dimension umintg(mxgran,mxst),umelstg(mxgran,mxst),
     1 umdisg(mxgran,mxst),umrepg(mxgran,mxst)
c
      dimension umclasg(mxgrpar),umcleg(mxgrpar),
     1 umcldg(mxgrpar),umclrg(mxgrpar)
c
      dimension ublat(mxst),ublaq(mxst),ublai(mxst),ublae(mxst),
     1 ublad(mxst),ublar(mxst),ublap(mxst),ublan(mxst),ublanp(mxst),
     2 ubladf(mxst*(mxst+1)/2)
      dimension ublat2(mxst),ublaq2(mxst),ublai2(mxst),ublae2(mxst),
     1 ublad2(mxst),ublar2(mxst),ublap2(mxst),ublan2(mxst),
     2 ublanp2(mxst),ubladf2(mxst*(mxst+1)/2)
c
      dimension ubacg(mxgrpar),ubaceg(mxgrpar),
     1 ubacdg(mxgrpar),ubacrg(mxgrpar)
      dimension ubacg2(mxgrpar),ubaceg2(mxgrpar),
     1 ubacdg2(mxgrpar),ubacrg2(mxgrpar)
c
      dimension ublaig(mxgran,mxst),ublaeg(mxgran,mxst),
     1 ubladg(mxgran,mxst),ublarg(mxgran,mxst)
      dimension ublaig2(mxgran,mxst),ublaeg2(mxgran,mxst),
     1 ubladg2(mxgran,mxst),ublarg2(mxgran,mxst)
c
      dimension usamat(mxst),usamaq(mxst),usamai(mxst),usamae(mxst),
     1 usamad(mxst),usamar(mxst),usamap(mxst),usaman(mxst),
     2 usmanp(mxst),usmadf(mxst*(mxst+1)/2)
      dimension usamat2(mxst),usamaq2(mxst),usamai2(mxst),usamae2(mxst),
     1 usamad2(mxst),usamar2(mxst),usamap2(mxst),usaman2(mxst),
     2 usmanp2(mxst),usmadf2(mxst*(mxst+1)/2)
c
      dimension usacg(mxgrpar),usaceg(mxgrpar),
     1 usacdg(mxgrpar), usacrg(mxgrpar)
      dimension usacg2(mxgrpar),usaceg2(mxgrpar),
     1 usacdg2(mxgrpar), usacrg2(mxgrpar)
c
      dimension usmaig(mxgran,mxst),usmaeg(mxgran,mxst),
     1 usmadg(mxgran,mxst),usmarg(mxgran,mxst)
      dimension usmaig2(mxgran,mxst),usmaeg2(mxgran,mxst),
     1 usmadg2(mxgran,mxst),usmarg2(mxgran,mxst)
c
      dimension uacat(mxst),uacaq(mxst),uacai(mxst),uacae(mxst),
     1 uacad(mxst),uacar(mxst),uacap(mxst),uacan(mxst),uacanp(mxst),
     2 uacadf(mxst*(mxst+1)/2)
      dimension uacat2(mxst),uacaq2(mxst),uacai2(mxst),uacae2(mxst),
     1 uacad2(mxst),uacar2(mxst),uacap2(mxst),uacan2(mxst),
     2 uacanp2(mxst),uacadf2(mxst*(mxst+1)/2)
c
      dimension uaacg(mxgrpar), uaaceg(mxgrpar),
     1 uaacdg(mxgrpar),uaacrg(mxgrpar)
      dimension uaacg2(mxgrpar), uaaceg2(mxgrpar),
     1 uaacdg2(mxgrpar),uaacrg2(mxgrpar)
c
      dimension uaaig(mxgran,mxst),uaaeg(mxgran,mxst),
     1 uaadg(mxgran,mxst),uaarg(mxgran,mxst)
      dimension uaaig2(mxgran,mxst),uaaeg2(mxgran,mxst),
     1 uaadg2(mxgran,mxst),uaarg2(mxgran,mxst)
c
      data zero,one,two,four,pt3 /0.0d00,1.0d00,2.0d00,4.0d00,
     1       .3333333333333333333d00/
c
      data exptol /1.0d+02/
c
      data string1,string2,string3,string4 /' block averages ',
     1' accumulated block averages', ' sample averages ',
     2' accumulated sample averages '/
c
      data header, estring1, estring2 /'configuration # ',
     1 'total energy  : ', 'class energy  : '/
      data comma /','/
c
      data incf /51/
      data nstate /1/
c
c-----  start monte carlo run
c
      obeen = .false.
      obeen4 = .false.
      if (imcout .ge. 2) then
        write(iwr,1)
    1   format(//,' monte carlo run information:',/,
     1  ' #sample #block #moves     #accepted   average tot ener ',
     2  ' rms deviation')
      endif
c
      nstate = imcst
      ngrpair = ngran*(ngran+1)/2
c
c-----  initialize total number of attempted moves
c
      nmovt = 0
      iaccpt = 0
c
c-----  loop over samples
c
c-----  initialize accumulative totals
c
      itotacc = 0
c
      uruntc = zero
      uruntc2 = zero
c
      urtce = zero
      urtce2 = zero
c
      urtcd = zero
      urtcd2 = zero
c
      urtcr = zero
      urtcr2 = zero
c
      do 100, k = 1, ngrpair
c 1-----
        urtcg(k) = zero
        urtcg2(k) = zero
c
        urtceg(k) = zero
        urtceg2(k) = zero
c
        urtcdg(k) = zero
        urtcdg2(k) = zero
c
        urtcrg(k) = zero
        urtcrg2(k) = zero
c 1-----
  100 continue
c
      call clear (urtdf,imcst*(imcst+1)/2)
      call clear (urtdf2,imcst*(imcst+1)/2)
      do 200, ist = 1, imcst
c 1-----
        urunt(ist) = zero
        urunt2(ist) = zero
c
        uruntq(ist) = zero
        uruntq2(ist) = zero
c
        urunti(ist) = zero
        urunti2(ist) = zero
c
        urunte(ist) = zero
        urunte2(ist) = zero
c
        uruntd(ist) = zero
        uruntd2(ist) = zero
c
        uruntr(ist) = zero
        uruntr2(ist) = zero
c
        uruntp(ist) = zero
        uruntp2(ist) = zero
c
        uruntn(ist) = zero
        uruntn2(ist) = zero
c
        urtnp(ist) = zero
        urtnp2(ist) = zero
c
        q1tot(ist) = zero
c
        do 170, igr = 1, ngran
c   2-----
          urtig(igr,ist) = zero
          urtig2(igr,ist) = zero
c
          urteg(igr,ist) = zero
          urteg2(igr,ist) = zero
c
          urtdg(igr,ist) = zero
          urtdg2(igr,ist) = zero
c
          urtrg(igr,ist) = zero
          urtrg2(igr,ist) = zero
c   2-----
  170   continue
c 1-----
  200 continue
c
      do 9000, isamp = 1, nsamp
c 1-----
c  -----  initialize sample totals
c
        isamacc = 0
c
        usamtc = zero
        usamtc2 = zero
c
        usmce = zero
        usmce2 = zero
c
        usmcd = zero
        usmcd2 = zero
c
        usmcr = zero
        usmcr2 = zero
c
        do 8100, k = 1, ngrpair
c   2-----
          usmtcg(k) = zero
          usmtcg2(k) = zero
c
          usmceg(k) = zero
          usmceg2(k) = zero
c
          usmcdg(k) = zero
          usmcdg2(k) = zero
c
          usmcrg(k) = zero
          usmcrg2(k) = zero
c   2-----
 8100   continue
c
        call clear (usmtdf,imcst*(imcst+1)/2)
        call clear (usmtdf2,imcst*(imcst+1)/2)
        do 8150, ist = 1, imcst
c   2-----
          usamt(ist) = zero
          usamt2(ist) = zero
c
          usamtq(ist) = zero
          usamtq2(ist) = zero
c
          usamti(ist) = zero
          usamti2(ist) = zero
c
          usamte(ist) = zero
          usamte2(ist) = zero
c
          usamtd(ist) = zero
          usamtd2(ist) = zero
c
          usamtr(ist) = zero
          usamtr2(ist) = zero
c
          usamtp(ist) = zero
          usamtp2(ist) = zero
c
          usamtn(ist) = zero
          usamtn2(ist) = zero
c
          usmtnp(ist) = zero
          usmtnp2(ist) = zero
c
          do 8120, igr = 1, ngran
c     3-----
            usmtig(igr,ist) = zero
            usmtig2(igr,ist) = zero
c
            usmteg(igr,ist) = zero
            usmteg2(igr,ist) = zero
c
            usmtdg(igr,ist) = zero
            usmtdg2(igr,ist) = zero
c
            usmtrg(igr,ist) = zero
            usmtrg2(igr,ist) = zero
c     3-----
 8120     continue
c   2-----
 8150   continue
c
c  -----  loop over blocks
c
c  -----  initialize accumulative block totals
c
        iblacct = 0
c
        ubltc = zero
        ubltc2 = zero
c
        ublce = zero
        ublce2 = zero
c
        ublcd = zero
        ublcd2 = zero
c
        ublcr = zero
        ublcr2 = zero
c
        do 8200, k = 1, ngrpair
c   2-----
          ubltcg(k) = zero
          ubltcg2(k) = zero
c
          ublceg(k) = zero
          ublceg2(k) = zero
c
          ublcdg(k) = zero
          ublcdg2(k) = zero
c
          ublcrg(k) = zero
          ublcrg2(k) = zero
c   2-----
 8200   continue
c
        call clear (ubltdf,imcst*(imcst+1)/2)
        call clear (ubltdf2,imcst*(imcst+1)/2)
        do 8250, ist = 1, imcst
c   2-----
          ublt(ist) = zero
          ublt2(ist) = zero
c
          ubltq(ist) = zero
          ubltq2(ist) = zero
c
          ublti(ist) = zero
          ublti2(ist) = zero
c
          ublte(ist) = zero
          ublte2(ist) = zero
c
          ubltd(ist) = zero
          ubltd2(ist) = zero
c
          ubltr(ist) = zero
          ubltr2(ist) = zero
c
          ubltp(ist) = zero
          ubltp2(ist) = zero
c
          ubltn(ist) = zero
          ubltn2(ist) = zero
c
          ubltnp(ist) = zero
          ubltnp2(ist) = zero
c
          do 8220, igr = 1, ngran
c     3-----
            ubltig(igr,ist) = zero
            ubltig2(igr,ist) = zero
c
            ublteg(igr,ist) = zero
            ublteg2(igr,ist) = zero
c
            ubltdg(igr,ist) = zero
            ubltdg2(igr,ist) = zero
c
            ubltrg(igr,ist) = zero
            ubltrg2(igr,ist) = zero
c     3-----
 8220     continue
c   2-----
 8250   continue
c
        do 8300, iblock = 1, nblock
c   2-----
c    -----  initialize block totals
c
          iblacc = 0
c
          ucbl = zero
          ucbl2 = zero
c
          ucebl = zero
          ucebl2 = zero
c
          ucdbl = zero
          ucdbl2 = zero
c
          ucrbl = zero
          ucrbl2 = zero
c
          do 7200, k = 1, ngrpair
c     3-----
            ucgbl(k) = zero
            ucgbl2(k) = zero
c
            ucegbl(k) = zero
            ucegbl2(k) = zero
c
            ucdgbl(k) = zero
            ucdgbl2(k) = zero
c
            ucrgbl(k) = zero
            ucrgbl2(k) = zero
c     3-----
 7200     continue
c
          call clear (udfbl,imcst*(imcst+1)/2)
          call clear (udfbl2,imcst*(imcst+1)/2)
          do 7400, ist = 1, imcst
c     3-----
            utbl(ist) = zero
            utbl2(ist) = zero
c
            uqbl(ist) = zero
            uqbl2(ist) = zero
c
            uibl(ist) = zero
            uibl2(ist) = zero
c
            uebl(ist) = zero
            uebl2(ist) = zero
c
            udbl(ist) = zero
            udbl2(ist) = zero
c
            urbl(ist) = zero
            urbl2(ist) = zero
c
            upbl(ist) = zero
            upbl2(ist) = zero
c
            unbl(ist) = zero
            unbl2(ist) = zero
c
            unpbl(ist) = zero
            unpbl2(ist) = zero
c
            do 7370, igr = 1, ngran
c       4-----
              uigbl(igr,ist) = zero
              uigbl2(igr,ist) = zero
c
              uegbl(igr,ist) = zero
              uegbl2(igr,ist) = zero
c
              udgbl(igr,ist) = zero
              udgbl2(igr,ist) = zero
c
              urgbl(igr,ist) = zero
              urgbl2(igr,ist) = zero
c       4-----
 7370       continue
c     3-----
 7400     continue
c
c    -----  loop over moves
c
          do 7500, imove = 1, nmoves
c     3-----
c      -----  calculate attempted move
c             and corresponding energy,
c             and update wt, vr and omega if move is accepted
c
            call mcstepp(nmovt,nstate,umove,umqm,umcl,umcle,umcld,umclr,
     1           umint,umelst,umdis,umrep,umpol,umneq,umneqp,udiff,
     2           umintg,umelstg,umdisg,umrepg,umclasg,umcleg,umcldg,
     3           umclrg,accept,xscm)
c
c      -----  write attempted configuration, energy and
c             acceptance on tape51, if required (imcout.ge.3)
c
            if (imcout .ge. 3) then
c       4-----
              if (accept) then
                astring = 'accepted'
              else
                astring = 'rejected'
              endif
c
              write (incf,2001) header, nmovt+1 , astring
 2001         format(2x,a16,1x,i6,2x,a8)
c
              if ((accept) .or. (outfor(5:7) .eq. 'all')) then
                write (incf,2002) estring1, uens(1), estring2, eclas
 2002           format(2x,a16,2x,e20.14,2x,a16,2x,e20.14)
              endif
c
              if (nopes .eq. 1) then
                mcgrp = 1 + mod(nmovt,ngrpmc)
                ngrp = ibitmc(mcgrp)
              else
                ngrp = npesgrp
              endif
c
              if (accept) then
c         5-----
                if (outfor(5:7) .eq. 'all') then
c           6-----
c            -----  write all energy components
c
                  write(incf,2003) umqm(1),umcle,umcld,umclr
 2003             format(4e20.14)
c
                  write(incf,2003) umint(1),umelst(1),umdis(1),umrep(1)
                  write(incf,2003) umpol(1),umneq(1),umneqp(1)
c           6-----
                endif
c
                if (outfor(9:11) .eq. 'all') then
c           6-----
                  do 6100, igr = 1, ngran
c             7-----
                    do 6060, jgr = 1, igr
c               8-----
                      indxcl = ia(igr) + jgr
                      write(incf,2003)
     1 umclasg(indxcl), umcleg(indxcl), umcldg(indxcl), umclrg(indxcl)
c               8-----
 6060               continue
c
                    write(incf,2003)
     1     umintg(igr,1),umelstg(igr,1),umdisg(igr,1),umrepg(igr,1)
c             7-----
 6100             continue
c           6-----
                endif
c
                if (imcst .gt. 1) then
c           6-----
                  do 6150, ist = 2, imcst
c             7-----
                    if (outfor(5:7) .eq. 'all') then
c               8-----
                      write(incf,2003) uens(ist),umqm(ist)
                      write(incf,2003)
     1              umint(ist),umelst(ist),umdis(ist),umrep(ist)
                      write(incf,2003) umpol(ist),umneq(ist),umneqp(ist)
c               8-----
                    endif
c
                    if (outfor(9:11) .eq. 'all') then
c               8-----
                      do 6120, igr = 1, ngran
                        write(incf,2003)
     1                  umintg(igr,ist),umelstg(igr,ist),
     2                  umdisg(igr,ist),umrepg(igr,ist)
 6120                 continue
c               8-----
                    endif
c             7-----
 6150             continue
c           6-----
                endif
c
cmas added the acc possibility :
c
                if (outfor(1:4) .eq. 'long' .or.
     .              outfor(1:3) .eq. 'acc' ) then
c           6-----
c            -----  write complete conformation
c
                  do 6200, igrp = 1, ngrpol
c             7-----
                    istr = igrpst(igrp)
                    write(incf,9940) nxcent(istr),comma,chrg(istr),
     1    comma,xpts(1,istr),comma,xpts(2,istr),comma,xpts(3,istr)
c
c                    write(incf,2004) igrp, ngrppt(igrp),
c     1                   (xpts(k,istr), k = 1, 3)
c 2004               format(2i4,3f20.6)
c
                    do 6160, ipart = 1, ngrppt(igrp)
                      write(incf,9940) nxcent(istr+ipart),comma,
     1    chrg(istr+ipart),comma,xpts(1,istr+ipart),comma,
     2    xpts(2,istr+ipart),comma,xpts(3,istr+ipart),
     3    comma,alfext(istr+ipart)
c
c                      write(incf,2004) igrp, ipart+1,
c     1                   (xpts(k,istr+ipart), k = 1, 3)
 6160               continue
c             7-----
 6200             continue
c           6-----
                else
c           6-----
c            -----  write only updated group
c
                  istr = igrpst(ngrp)
                  write(incf,9940) nxcent(istr),comma,
     1    chrg(istr),comma,xpts(1,istr),comma,
     2    xpts(2,istr),comma,xpts(3,istr)
c
c                  write(incf,2004) ngrp, ngrppt(ngrp),
c     1                   (xpts(k,igrpst(ngrp)), k = 1, 3)
                  do 6250, ipart = 1, ngrppt(ngrp)
                    write(incf,9940) nxcent(istr+ipart),comma,
     1    chrg(istr+ipart),comma,xpts(1,istr+ipart),comma,
     2    xpts(2,istr+ipart),comma,xpts(3,istr+ipart),
     3    comma,alfext(istr+ipart)
c
c                    write(incf,2004) ngrp, ipart+1,
c     1                   (xpts(k,igrpst(ngrp)+ipart), k = 1, 3)
 6250             continue
c           6-----
                endif
c         5-----
              else if (outfor(1:4) .eq. 'long') then
c         5-----
                if (outfor(5:7) .eq. 'all') then
c           6-----
c            -----  write all energy components
c
                  write(incf,2003) umqm(1),umcle,umcld,umclr
                  write(incf,2003) umint(1),umelst(1),umdis(1),umrep(1)
                  write(incf,2003) umpol(1),umneq(1),umneqp(1)
c           6-----
                endif
                if (outfor(9:11) .eq. 'all') then
c           6-----
                  do 6300, igr = 1, ngran
c             7-----
                    do 6260, jgr = 1, igr
c               8-----
                      indxcl = ia(igr) + jgr
                      write(incf,2003)
     1 umclasg(indxcl), umcleg(indxcl), umcldg(indxcl), umclrg(indxcl)
c               8-----
 6260               continue
c
                    write(incf,2003)
     1     umintg(igr,1),umelstg(igr,1),umdisg(igr,1),umrepg(igr,1)
c             7-----
 6300             continue
c           6-----
                endif
c
                if (imcst .gt. 1) then
c           6-----
                  do 6350, ist = 2, imcst
c             7-----
                    if (outfor(5:7) .eq. 'all') then
c               8-----
                      write(incf,2003) uens(ist),umqm(ist)
                      write(incf,2003)
     1              umint(ist),umelst(ist),umdis(ist),umrep(ist)
                      write(incf,2003) umpol(ist),umneq(ist),umneqp(ist)
c               8-----
                    endif
c
                    if (outfor(9:11) .eq. 'all') then
c               8-----
                      do 6310, igr = 1, ngran
                        write(incf,2003)
     1                  umintg(igr,ist),umelstg(igr,ist),
     2                  umdisg(igr,ist),umrepg(igr,ist)
 6310                 continue
c               8-----
                    endif
c             7-----
 6350             continue
c           6-----
                endif
c
                do 6400, igrp = 1, ngrpol
c           6-----
                  istr = igrpst(igrp)
                  if (igrp .eq. ngrp) then
c             7-----
                    write(incf,9940) nxcent(istr),comma,
     1    chrg(istr),comma,xnpts(1,1),comma,
     2    xnpts(2,1),comma,xnpts(3,1)
c
c                    write(incf,2004) igrp, ngrppt(igrp),
c     1                   (xnpts(k,1), k = 1, 3)
                    do 6360, ipart = 2, ngrppt(igrp) + 1
                      write(incf,9940) nxcent(istr+ipart-1),comma,
     1    chrg(istr+ipart-1),comma,xnpts(1,ipart),comma,
     2    xnpts(2,ipart),comma,xnpts(3,ipart),
     3    comma,alfext(istr+ipart-1)
c
c                      write(incf,2004) igrp, ipart,
c     1                     (xnpts(k,ipart), k = 1, 3)
 6360               continue
c             7-----
                  else
c             7-----
                    istr = igrpst(igrp)
                    write(incf,9940) nxcent(istr),comma,
     1    chrg(istr),comma,xpts(1,istr),comma,
     2    xpts(2,istr),comma,xpts(3,istr)
c
c                    write(incf,2004) igrp, ngrppt(igrp),
c     1                   (xpts(k,istr), k = 1, 3)
                    do 6370, ipart = 1, ngrppt(igrp)
                    write(incf,9940) nxcent(istr+ipart),comma,
     1    chrg(istr+ipart),comma,xpts(1,istr+ipart),comma,
     2    xpts(2,istr+ipart),comma,xpts(3,istr+ipart),
     3    comma,alfext(istr+ipart)
c
c                      write(incf,2004) igrp, ipart+1,
c     1                     (xpts(k,istr+ipart), k = 1, 3)
 6370               continue
c             7-----
                  endif
c           6-----
 6400           continue
c         5-----
              endif
c       4-----
            endif
c
c      -----  if move is accepted, count and store this move
c
            if (accept) then
c       4-----
              iblacc = iblacc + 1
              iaccpt = iaccpt + 1
c       4-----
            endif
c
            if ((ncheck .ne. 0) .and.
     2          (mod(nmovt+1,ncheck) .eq. 0)) then
c       4-----
c        -----  adjust parameters if required
c
              ratio = dble(iaccpt)/dble(ncheck)
c
              if (ratio .gt. ratmax) then
                if (darot .lt. amxrot) darot = darot*1.05d0
                if (dtrans .lt. amxtrn) dtrans = dtrans*1.05d0
              else if (ratio .lt. ratmin) then
                if (darot .gt. amnrot) darot = darot*0.95d0
                if (dtrans .gt. amntrn) dtrans = dtrans*0.95d0
              endif
c
              iaccpt = 0
c       4-----
            endif
c
            ucbl = ucbl + umcl
            ucbl2 = ucbl2 + umcl**2
c
            ucebl = ucebl + umcle
            ucebl2 = ucebl2 + umcle**2
c
            ucdbl = ucdbl + umcld
            ucdbl2 = ucdbl2 + umcld**2
c
            ucrbl = ucrbl + umclr
            ucrbl2 = ucrbl2 + umclr**2
c
            do 6500, k = 1, ngrpair
c       4-----
              ucgbl(k) = ucgbl(k) + umclasg(k)
              ucgbl2(k) = ucgbl2(k) + umclasg(k)**2
c
              ucegbl(k) = ucegbl(k) + umcleg(k)
              ucegbl2(k) = ucegbl2(k) + umcleg(k)**2
c
              ucdgbl(k) = ucdgbl(k) + umcldg(k)
              ucdgbl2(k) = ucdgbl2(k) + umcldg(k)**2
c
              ucrgbl(k) = ucrgbl(k) + umclrg(k)
              ucrgbl2(k) = ucrgbl2(k) + umclrg(k)**2
c       4-----
 6500       continue
c
            do 6550, ist = 1, imcst
c       4-----
              utbl(ist) = utbl(ist) + umove(ist)
              utbl2(ist) = utbl2(ist) + umove(ist)**2
c
              uqbl(ist) = uqbl(ist) + umqm(ist)
              uqbl2(ist) = uqbl2(ist) + umqm(ist)**2
c
              uibl(ist) = uibl(ist) + umint(ist)
              uibl2(ist) = uibl2(ist) + umint(ist)**2
c
              uebl(ist) = uebl(ist) + umelst(ist)
              uebl2(ist) = uebl2(ist) + umelst(ist)**2
c
              udbl(ist) = udbl(ist) + umdis(ist)
              udbl2(ist) = udbl2(ist) + umdis(ist)**2
c
              urbl(ist) = urbl(ist) + umrep(ist)
              urbl2(ist) = urbl2(ist) + umrep(ist)**2
c
              upbl(ist) = upbl(ist) + umpol(ist)
              upbl2(ist) = upbl2(ist) + umpol(ist)**2
c
              unbl(ist) = unbl(ist) + umneq(ist)
              unbl2(ist) = unbl2(ist) + umneq(ist)**2
c
              unpbl(ist) = unpbl(ist) + umneqp(ist)
              unpbl2(ist) = unpbl2(ist) + umneqp(ist)**2
c
              if (onekt*(umove(ist)-ureft(ist)) .lt. exptol)
     2          q1tot(ist) = q1tot(ist) +
     3                     exp(onekt*(umove(ist)-ureft(ist)))
c
              do 6510, kst = 1, ist
                indxdf = ia(ist) + kst
                udfbl(indxdf) = udfbl(indxdf) + udiff(indxdf)
                udfbl2(indxdf) = udfbl2(indxdf) + udiff(indxdf)**2
 6510         continue
c
              do 6520, igr = 1, ngran
c         5-----
                uigbl(igr,ist) = uigbl(igr,ist) + umintg(igr,ist)
                uigbl2(igr,ist) = uigbl2(igr,ist) + umintg(igr,ist)**2
c
                uegbl(igr,ist) = uegbl(igr,ist) + umelstg(igr,ist)
                uegbl2(igr,ist) = uegbl2(igr,ist) + umelstg(igr,ist)**2
c
                udgbl(igr,ist) = udgbl(igr,ist) + umdisg(igr,ist)
                udgbl2(igr,ist) = udgbl2(igr,ist) + umdisg(igr,ist)**2
c
                urgbl(igr,ist) = urgbl(igr,ist) + umrepg(igr,ist)
                urgbl2(igr,ist) = urgbl2(igr,ist) + umrepg(igr,ist)**2
c         5-----
 6520         continue
c       4-----
 6550       continue
c
            if (isolsav .eq. 1) then
c       4-----
c        -----  collect total external potential and field
c               at the expansion centra
c
              call addpot(nqdim,nqcls,itwoeps,field,
     2                    xscm(1),xscm(1+nqdim),xscm(1+2*nqdim))
c       4-----
            endif
c
c      -----  end of move
c
            nmovt = nmovt + 1
c     3-----  next move
 7500     continue
c
c    -----  end of block: collect block averages if wanted
c
          if (imcout .ge. 4) then
c     3-----
            ublac = ucbl / nmoves
            ublac2 = ucbl2 / nmoves
c
            ublace = ucebl / nmoves
            ublace2 = ucebl2 / nmoves
c
            ublacd = ucdbl / nmoves
            ublacd2 = ucdbl2 / nmoves
c
            ublacr = ucrbl / nmoves
            ublacr2 = ucrbl2 / nmoves
c
            do 7550, k = 1, ngrpair
c       4-----
              ubacg(k) = ucgbl(k) / nmoves
              ubacg2(k) = ucgbl2(k) / nmoves
c
              ubaceg(k) = ucegbl(k) / nmoves
              ubaceg2(k) = ucegbl2(k) / nmoves
c
              ubacdg(k) = ucdgbl(k) / nmoves
              ubacdg2(k) = ucdgbl2(k) / nmoves
c
              ubacrg(k) = ucrgbl(k) / nmoves
              ubacrg2(k) = ucrgbl2(k) / nmoves
c       4-----
 7550       continue
c
            do 7600, ist = 1, imcst
c       4-----
              ublat(ist) = utbl(ist) / nmoves
              ublat2(ist) = utbl2(ist) / nmoves
c
              ublaq(ist) = uqbl(ist) / nmoves
              ublaq2(ist) = uqbl2(ist) / nmoves
c
              ublai(ist) = uibl(ist) / nmoves
              ublai2(ist) = uibl2(ist) / nmoves
c
              ublae(ist) = uebl(ist) / nmoves
              ublae2(ist) = uebl2(ist) / nmoves
c
              ublad(ist) = udbl(ist) / nmoves
              ublad2(ist) = udbl2(ist) / nmoves
c
              ublar(ist) = urbl(ist) / nmoves
              ublar2(ist) = urbl2(ist) / nmoves
c
              ublap(ist) = upbl(ist) / nmoves
              ublap2(ist) = upbl2(ist) / nmoves
c
              ublan(ist) = unbl(ist) / nmoves
              ublan2(ist) = unbl2(ist) / nmoves
c
              ublanp(ist) = unpbl(ist) / nmoves
              ublanp2(ist) = unpbl2(ist) / nmoves
c
              do 7560, kst = 1, ist
                indxdf = ia(ist) + kst
                ubladf(indxdf) = udfbl(indxdf) / nmoves
                ubladf2(indxdf) = udfbl2(indxdf) / nmoves
 7560         continue
c
              do 7570, igr= 1, ngran
c         5-----
                ublaig(igr,ist) = uigbl(igr,ist) / nmoves
                ublaig2(igr,ist) = uigbl2(igr,ist) / nmoves
c
                ublaeg(igr,ist) = uegbl(igr,ist) / nmoves
                ublaeg2(igr,ist) = uegbl2(igr,ist) / nmoves
c
                ubladg(igr,ist) = udgbl(igr,ist) / nmoves
                ubladg2(igr,ist) = udgbl2(igr,ist) / nmoves
c
                ublarg(igr,ist) = urgbl(igr,ist) / nmoves
                ublarg2(igr,ist) = urgbl2(igr,ist) / nmoves
c         5-----
 7570         continue
c       4-----
 7600       continue
c
c    -----  print block averages
c
            call mcblout(nstate,isamp,iblock,
     1        nmoves, iblacc, ublat, ublat2, ublaq, ublaq2,
     2        ublac, ublac2, ublace, ublace2, ublacd, ublacd2,
     3        ublacr,ublacr2,ublai,ublai2,ublae,ublae2,ublad,ublad2,
     4        ublar,ublar2,ublap,ublap2,ublan,ublan2,ublanp,ublanp2,
     5        ubladf,ubladf2,ubacg,ubacg2,ubaceg,ubace2g,ubacdg,ubacdg2,
     6        ubacrg,ubacrg2,ublaig,ublaig2,ublaeg,ublaeg2,ubladg,
     7        ubladg2,ublarg,ublarg2,string1)
c     3-----
          endif
c
c    -----  collect accumulated totals
c
          iblacct = iblacct + iblacc
c
          ubltc = ubltc + ucbl
          ubltc2 = ubltc2 + ucbl2
c
          ublce = ublce + ucebl
          ublce2 = ublce2 + ucebl2
c
          ublcd = ublcd + ucdbl
          ublcd2 = ublcd2 + ucdbl2
c
          ublcr = ublcr + ucrbl
          ublcr2 = ublcr2 + ucrbl2
c
          do 7700, k = 1, ngrpair
c     3-----
            ubltcg(k) = ubltcg(k) + ucgbl(k)
            ubltcg2(k) = ubltcg2(k) + ucgbl2(k)
c
            ublceg(k) = ublceg(k) + ucegbl(k)
            ublceg2(k) = ublceg2(k) + ucegbl2(k)
c
            ublcdg(k) = ublcdg(k) + ucdgbl(k)
            ublcdg2(k) = ublcdg2(k) + ucdgbl2(k)
c
            ublcrg(k) = ublcrg(k) + ucrgbl(k)
            ublcrg2(k) = ublcrg2(k) + ucrgbl2(k)
c     3-----
 7700     continue
c
          do 7800, ist = 1, imcst
c     3-----
            ublt(ist) = ublt(ist) + utbl(ist)
            ublt2(ist) = ublt2(ist) + utbl2(ist)
c
            ubltq(ist) = ubltq(ist) + uqbl(ist)
            ubltq2(ist) = ubltq2(ist) + uqbl2(ist)
c
            ublti(ist) = ublti(ist) + uibl(ist)
            ublti2(ist) = ublti2(ist) + uibl2(ist)
c
            ublte(ist) = ublte(ist) + uebl(ist)
            ublte2(ist) = ublte2(ist) + uebl2(ist)
c
            ubltd(ist) = ubltd(ist) + udbl(ist)
            ubltd2(ist) = ubltd2(ist) + udbl2(ist)
c
            ubltr(ist) = ubltr(ist) + urbl(ist)
            ubltr2(ist) = ubltr2(ist) + urbl2(ist)
c
            ubltp(ist) = ubltp(ist) + upbl(ist)
            ubltp2(ist) = ubltp2(ist) + upbl2(ist)
c
            ubltn(ist) = ubltn(ist) + unbl(ist)
            ubltn2(ist) = ubltn2(ist) + unbl2(ist)
c
            ubltnp(ist) = ubltnp(ist) + unpbl(ist)
            ubltnp2(ist) = ubltnp2(ist) + unpbl2(ist)
c
            do 7760, kst = 1, ist
              indxdf = ia(ist) + kst
              ubltdf(indxdf) = ubltdf(indxdf) + udfbl(indxdf)
              ubltdf2(indxdf) = ubltdf2(indxdf) + udfbl2(indxdf)
 7760       continue
c
            do 7770, igr = 1, ngran
c       4-----
              ubltig(igr,ist) = ubltig(igr,ist) + uigbl(igr,ist)
              ubltig2(igr,ist) = ubltig2(igr,ist) + uigbl2(igr,ist)
c
              ublteg(igr,ist) = ublteg(igr,ist) + uegbl(igr,ist)
              ublteg2(igr,ist) = ublteg2(igr,ist) + uegbl2(igr,ist)
c
              ubltdg(igr,ist) = ubltdg(igr,ist) + udgbl(igr,ist)
              ubltdg2(igr,ist) = ubltdg2(igr,ist) + udgbl2(igr,ist)
c
              ubltrg(igr,ist) = ubltrg(igr,ist) + urgbl(igr,ist)
              ubltrg2(igr,ist) = ubltrg2(igr,ist) + urgbl2(igr,ist)
c       4-----
 7770       continue
c     3-----
 7800     continue
c
c    -----  calculate accumulated block averages if wanted
c
          if (imcout .ge. 3) then
c     3-----
            uacac = ubltc / (iblock*nmoves)
            uacac2 = ubltc2 / (iblock*nmoves)
c
            uacace = ublce / (iblock*nmoves)
            uacace2 = ublce2 / (iblock*nmoves)
c
            uacacd = ublcd / (iblock*nmoves)
            uacacd2 = ublcd2 / (iblock*nmoves)
c
            uacacr = ublcr / (iblock*nmoves)
            uacacr2 = ublcr2 / (iblock*nmoves)
c
            do 7900, k = 1, ngrpair
c       5-----
              uaacg(k) = ubltcg(k) / (iblock*nmoves)
              uaacg2(k) = ubltcg2(k) / (iblock*nmoves)
c
              uaaceg(k) = ublceg(k) / (iblock*nmoves)
              uaaceg2(k) = ublceg2(k) / (iblock*nmoves)
c
              uaacdg(k) = ublcdg(k) / (iblock*nmoves)
              uaacdg2(k) = ublcdg2(k) / (iblock*nmoves)
c
              uaacrg(k) = ublcrg(k) / (iblock*nmoves)
              uaacrg2(k) = ublcrg2(k) / (iblock*nmoves)
c       5-----
 7900       continue
c
            do 7920, ist = 1, imcst
c       4-----
              uacat(ist) = ublt(ist) / (iblock*nmoves)
              uacat2(ist) = ublt2(ist) / (iblock*nmoves)
c
              uacaq(ist) = ubltq(ist) / (iblock*nmoves)
              uacaq2(ist) = ubltq2(ist) / (iblock*nmoves)
c
              uacai(ist) = ublti(ist) / (iblock*nmoves)
              uacai2(ist) = ublti2(ist) / (iblock*nmoves)
c
              uacae(ist) = ublte(ist) / (iblock*nmoves)
              uacae2(ist) = ublte2(ist) / (iblock*nmoves)
c
              uacad(ist) = ubltd(ist) / (iblock*nmoves)
              uacad2(ist) = ubltd2(ist) / (iblock*nmoves)
c
              uacar(ist) = ubltr(ist) / (iblock*nmoves)
              uacar2(ist) = ubltr2(ist) / (iblock*nmoves)
c
              uacap(ist) = ubltp(ist) / (iblock*nmoves)
              uacap2(ist) = ubltp2(ist) / (iblock*nmoves)
c
              uacan(ist) = ubltn(ist) / (iblock*nmoves)
              uacan2(ist) = ubltn2(ist) / (iblock*nmoves)
c
              uacanp(ist) = ubltnp(ist) / (iblock*nmoves)
              uacanp2(ist) = ubltnp2(ist) / (iblock*nmoves)
c
              do 7912, kst = 1, ist
                indxdf = ia(ist) + kst
                uacadf(indxdf) = ubltdf(indxdf) / (iblock*nmoves)
                uacadf2(indxdf) = ubltdf2(indxdf) / (iblock*nmoves)
 7912         continue
c
              do 7914, igr = 1, ngran
c         5-----
                uaaig(igr,ist) = ubltig(igr,ist) / (iblock*nmoves)
                uaaig2(igr,ist) = ubltig2(igr,ist) / (iblock*nmoves)
c
                uaaeg(igr,ist) = ublteg(igr,ist) / (iblock*nmoves)
                uaaeg2(igr,ist) = ublteg2(igr,ist) / (iblock*nmoves)
c
                uaadg(igr,ist) = ubltdg(igr,ist) / (iblock*nmoves)
                uaadg2(igr,ist) = ubltdg2(igr,ist) / (iblock*nmoves)
c
                uaarg(igr,ist) = ubltrg(igr,ist) / (iblock*nmoves)
                uaarg2(igr,ist) = ubltrg2(igr,ist) / (iblock*nmoves)
c         5-----
 7914         continue
c       4-----
 7920       continue
c
c    -----  print accumulated block averages
c
            call mcblout(nstate,isamp,iblock,
     1        iblock*nmoves,iblacct,uacat,uacat2,uacaq,uacaq2,
     2        uacac,uacac2,uacace,uacace2,uacacd,uacacd2,
     3        uacacr,uacacr2,uacai,uacai2,uacae,uacae2,uacad,uacad2,
     4        uacar,uacar2,uacap,uacap2,uacan,uacan2,uacanp,uacanp2,
     5        uacadf,uacadf2,uaacg,uaacg2,uaaceg,uaaceg2,uaacdg,uaacdg2,
     6        uaacrg,uaacrg2,uaaig,uaaig2,uaaeg,uaaeg2,uaadg,
     7        uaadg2,uaarg,uaarg2,string2)
c
c      -----  write current conformation
c
            write(iwr,8834)
 8834       format(/,'  -- current conformation')
            write(iwr,9935)
            imp = 1
            img = 1
c
            do 7930, ii = 1, nxtpts
c       4-----
              if (ii .eq. mpol(imp)) then
c         5-----
c                write(iwr,9940) nxcent(ii),comma,chrg(ii),comma,
c     1    xpts(1,ii),comma,xpts(2,ii),comma,xpts(3,ii),comma,
c     2                   polar(imp)
                imp = imp + 1
              else
                write(iwr,9940) nxcent(ii),comma,chrg(ii),comma,
     1    xpts(1,ii),comma,xpts(2,ii),comma,xpts(3,ii),comma,
     2    alfext(ii)
c         5-----
              endif
              if (ii .eq. igrlast(img)) then
                write (iwr,9950) grlabel(img)
                img = img + 1
              endif
c       4-----
 7930       continue
c     3-----
          endif
c   2-----  next block
 8300   continue
c
c  -----  end of sample: collect sample averages if wanted
c
        isamacc = iblacct
c
        usamtc = ubltc
        usamtc2 = ubltc2
c
        usmce = ublce
        usmce2 = ublce2
c
        usmcd = ublcd
        usmcd2 = ublcd2
c
        usmcr = ublcr
        usmcr2 = ublcr2
c
        do 8350, k = 1, ngrpair
c   2-----
          usmtcg(k) = ubltcg(k)
          usmtcg2(k) = ubltcg2(k)
c
          usmceg(k) = ublceg(k)
          usmceg2(k) = ublceg2(k)
c
          usmcdg(k) = ublcdg(k)
          usmcdg2(k) = ublcdg2(k)
c
          usmcrg(k) = ublcrg(k)
          usmcrg2(k) = ublcrg2(k)
c   2-----
 8350   continue
c
        do 8400, ist = 1, imcst
c   2-----
          usamt(ist) = ublt(ist)
          usamt2(ist) = ublt2(ist)
c
          usamtq(ist) = ubltq(ist)
          usamtq2(ist) = ubltq2(ist)
c
          usamti(ist) = ublti(ist)
          usamti2(ist) = ublti2(ist)
c
          usamte(ist) = ublte(ist)
          usamte2(ist) = ublte2(ist)
c
          usamtd(ist) = ubltd(ist)
          usamtd2(ist) = ubltd2(ist)
c
          usamtr(ist) = ubltr(ist)
          usamtr2(ist) = ubltr2(ist)
c
          usamtp(ist) = ubltp(ist)
          usamtp2(ist) = ubltp2(ist)
c
          usamtn(ist) = ubltn(ist)
          usamtn2(ist) = ubltn2(ist)
c
          usmtnp(ist) = ubltnp(ist)
          usmtnp2(ist) = ubltnp2(ist)
c
          do 8360, kst = 1, ist
            indxdf = ia(ist) + kst
            usmtdf(indxdf) = ubltdf(indxdf)
            usmtdf2(indxdf) = ubltdf2(indxdf)
 8360     continue
c
          do  8370, igr = 1, ngran
c     3-----
            usmtig(igr,ist) = ubltig(igr,ist)
            usmtig2(igr,ist) = ubltig2(igr,ist)
c
            usmteg(igr,ist) = ublteg(igr,ist)
            usmteg2(igr,ist) = ublteg2(igr,ist)
c
            usmtdg(igr,ist) = ubltdg(igr,ist)
            usmtdg2(igr,ist) = ubltdg2(igr,ist)
c
            usmtrg(igr,ist) = ubltrg(igr,ist)
            usmtrg2(igr,ist) = ubltrg2(igr,ist)
c     3-----
 8370     continue
c   2-----
 8400   continue
c
        if (imcout .ge. 2) then
c   2-----
          usamac = usamtc / (nblock*nmoves)
          usamac2 = usamtc2 / (nblock*nmoves)
c
          usmace = usmce / (nblock*nmoves)
          usmace2 = usmce2 / (nblock*nmoves)
c
          usmacd = usmcd / (nblock*nmoves)
          usmacd2 = usmcd2 / (nblock*nmoves)
c
          usmacr = usmcr / (nblock*nmoves)
          usmacr2 = usmcr2 / (nblock*nmoves)
c
          do 8450, k = 1, ngrpair
c     3-----
            usacg(k) = usmtcg(k) / (nblock*nmoves)
            usacg2(k) = usmtcg2(k) / (nblock*nmoves)
c
            usaceg(k) = usmceg(k) / (nblock*nmoves)
            usaceg2(k) = usmceg2(k) / (nblock*nmoves)
c
            usacdg(k) = usmcdg(k) / (nblock*nmoves)
            usacdg2(k) = usmcdg2(k) / (nblock*nmoves)
c
            usacrg(k) = usmcrg(k) / (nblock*nmoves)
            usacrg2(k) = usmcrg2(k) / (nblock*nmoves)
c     3-----
 8450     continue
c
          do 8500, ist = 1, imcst
c     3-----
            usamat(ist) = usamt(ist) / (nblock*nmoves)
            usamat2(ist) = usamt2(ist) / (nblock*nmoves)
c
            usamaq(ist) = usamtq(ist) / (nblock*nmoves)
            usamaq2(ist) = usamtq2(ist) / (nblock*nmoves)
c
            usamai(ist) = usamti(ist) / (nblock*nmoves)
            usamai2(ist) = usamti2(ist) / (nblock*nmoves)
c
            usamae(ist) = usamte(ist) / (nblock*nmoves)
            usamae2(ist) = usamte2(ist) / (nblock*nmoves)
c
            usamad(ist) = usamtd(ist) / (nblock*nmoves)
            usamad2(ist) = usamtd2(ist) / (nblock*nmoves)
c
            usamar(ist) = usamtr(ist) / (nblock*nmoves)
            usamar2(ist) = usamtr2(ist) / (nblock*nmoves)
c
            usamap(ist) = usamtp(ist) / (nblock*nmoves)
            usamap2(ist) = usamtp2(ist) / (nblock*nmoves)
c
            usaman(ist) = usamtn(ist) / (nblock*nmoves)
            usaman2(ist) = usamtn2(ist) / (nblock*nmoves)
c
            usmanp(ist) = usmtnp(ist) / (nblock*nmoves)
            usmanp2(ist) = usmtnp2(ist) / (nblock*nmoves)
c
            do 8460, kst = 1, ist
              indxdf = ia(ist) + kst
              usmadf(indxdf) = usmtdf(indxdf) / (nblock*nmoves)
              usmadf2(indxdf) = usmtdf2(indxdf) / (nblock*nmoves)
 8460       continue
c
            do 8470, igr = 1, ngran
c       4-----
              usmaig(igr,ist) = usmtig(igr,ist) / (nblock*nmoves)
              usmaig2(igr,ist) = usmtig2(igr,ist) / (nblock*nmoves)
c
              usmaeg(igr,ist) = usmteg(igr,ist) / (nblock*nmoves)
              usmaeg2(igr,ist) = usmteg2(igr,ist) / (nblock*nmoves)
c
              usmadg(igr,ist) = usmtdg(igr,ist) / (nblock*nmoves)
              usmadg2(igr,ist) = usmtdg2(igr,ist) / (nblock*nmoves)
c
              usmarg(igr,ist) = usmtrg(igr,ist) / (nblock*nmoves)
              usmarg2(igr,ist) = usmtrg2(igr,ist) / (nblock*nmoves)
c       4-----
 8470       continue
c     3-----
 8500     continue
c
c  -----  print sample averages
c
          call mcblout(nstate,isamp,iblock,
     1        nblock*nmoves,iblacct,usamat,usamat2,usamaq,usamaq2,
     2        usamac,usamac2,usmace,usmace2,usmacd,usmacd2,
     3        usmacr,usmacr2,usamai,usamai2,usamae,usamae2,
     4        usamad,usamad2,usamar,usamar2,usamap,usamap2,
     5        usaman,usaman2,usmanp,usmanp2,usmadf,usmadf2,
     6        usmacg,usmacg2,usmaceg,usmaceg2,usmacdg,usmacdg2,
     7        usmacrg,usmacrg2,usmaig,usmaig2,usmaeg,usmaeg2,
     8        usmadg,usmadg2,usmarg,usmarg2,string2)
c   2-----
        endif
c
c  -----  collect accumulated totals
c
        itotacc = itotacc + isamacc
c
        uruntc = uruntc + usamtc
        uruntc2 = uruntc2 + usamtc2
c
        urtce = urtce + usmce
        urtce2 = urtce2 + usmce2
c
        urtcd = urtcd + usmcd
        urtcd2 = urtcd2 + usmcd2
c
        urtcr = urtcr + usmcr
        urtcr2 = urtcr2 + usmcr2
c
        do 8550, k = 1, ngrpair
c   2-----
          urtcg(k) = urtcg(k) + usmtcg(k)
          urtcg2(k) = urtcg2(k) + usmtcg2(k)
c
          urtceg(k) = urtceg(k) + usmceg(k)
          urtceg2(k) = urtceg2(k) + usmceg2(k)
c
          urtcdg(k) = urtcdg(k) + usmcdg(k)
          urtcdg2(k) = urtcdg2(k) + usmcdg2(k)
c
          urtcrg(k) = urtcrg(k) + usmcrg(k)
          urtcrg2(k) = urtcrg2(k) + usmcrg2(k)
c   2-----
 8550   continue
c
        do 8600, ist = 1, imcst
c   2-----
          urunt(ist) = urunt(ist) + usamt(ist)
          urunt2(ist) = urunt2(ist) + usamt2(ist)
c
          uruntq(ist) = uruntq(ist) + usamtq(ist)
          uruntq2(ist) = uruntq2(ist) + usamtq2(ist)
c
          urunti(ist) = urunti(ist) + usamti(ist)
          urunti2(ist) = urunti2(ist) + usamti2(ist)
c
          urunte(ist) = urunte(ist) + usamte(ist)
          urunte2(ist) = urunte2(ist) + usamte2(ist)
c
          uruntd(ist) = uruntd(ist) + usamtd(ist)
          uruntd2(ist) = uruntd2(ist) + usamtd2(ist)
c
          uruntr(ist) = uruntr(ist) + usamtr(ist)
          uruntr2(ist) = uruntr2(ist) + usamtr2(ist)
c
          uruntp(ist) = uruntp(ist) + usamtp(ist)
          uruntp2(ist) = uruntp2(ist) + usamtp2(ist)
c
          uruntn(ist) = uruntn(ist) + usamtn(ist)
          uruntn2(ist) = uruntn2(ist) + usamtn2(ist)
c
          urtnp(ist) = urtnp(ist) + usmtnp(ist)
          urtnp2(ist) = urtnp2(ist) + usmtnp2(ist)
c
          do 8560, kst = 1, ist
            indxdf = ia(ist) + kst
            urtdf(indxdf) = urtdf(indxdf) + usmtdf(indxdf)
            urtdf2(indxdf) = urtdf2(indxdf) + usmtdf2(indxdf)
 8560     continue
c
          do 8570, igr = 1, ngran
c     3-----
          urtig(igr,ist) = urtig(igr,ist) + usmtig(igr,ist)
          urtig2(igr,ist) = urtig2(igr,ist) + usmtig2(igr,ist)
c
          urteg(igr,ist) = urteg(igr,ist) + usmteg(igr,ist)
          urteg2(igr,ist) = urteg2(igr,ist) + usmteg2(igr,ist)
c
          urtdg(igr,ist) = urtdg(igr,ist) + usmtdg(igr,ist)
          urtdg2(igr,ist) = urtdg2(igr,ist) + usmtdg2(igr,ist)
c
          urtrg(igr,ist) = urtrg(igr,ist) + usmtrg(igr,ist)
          urtrg2(igr,ist) = urtrg2(igr,ist) + usmtrg2(igr,ist)
c     3-----
 8570     continue
c   2-----
 8600   continue
c
c    -----  calculate accumulated sample averages if wanted
c
        if (imcout .ge. 2) then
c   2-----
          uacac = uruntc / (isamp*nblock*nmoves)
          uacac2 = uruntc2 / (isamp*nblock*nmoves)
c
          uacace = urtce / (isamp*nblock*nmoves)
          uacace2 = urtce2 / (isamp*nblock*nmoves)
c
          uacacd = urtcd / (isamp*nblock*nmoves)
          uacacd2 = urtcd2 / (isamp*nblock*nmoves)
c
          uacacr = urtcr / (isamp*nblock*nmoves)
          uacacr2 = urtcr2 / (isamp*nblock*nmoves)
c
          do 8700, k = 1, ngrpair
c     3-----
            uaacg(k) = urtcg(k) / (isamp*nblock*nmoves)
            uaacg2(k) = urtcg2(k) / (isamp*nblock*nmoves)
c
            uaaceg(k) = urtceg(k) / (isamp*nblock*nmoves)
            uaaceg2(k) = urtceg2(k) / (isamp*nblock*nmoves)
c
            uaacdg(k) = urtcdg(k) / (isamp*nblock*nmoves)
            uaacdg2(k) = urtcdg2(k) / (isamp*nblock*nmoves)
c
            uaacrg(k) = urtcrg(k) / (isamp*nblock*nmoves)
            uaacrg2(k) = urtcrg2(k) / (isamp*nblock*nmoves)
c     3-----
 8700     continue
c
          do 8800, ist = 1, imcst
c     3-----
            uacat(ist) = urunt(ist) / (isamp*nblock*nmoves)
            uacat2(ist) = urunt2(ist) /(isamp*nblock*nmoves)
c
            uacaq(ist) = uruntq(ist) / (isamp*nblock*nmoves)
            uacaq2(ist) = uruntq2(ist) / (isamp*nblock*nmoves)
c
            uacai(ist) = urunti(ist) / (isamp*nblock*nmoves)
            uacai2(ist) = urunti2(ist) / (isamp*nblock*nmoves)
c
            uacae(ist) = urunte(ist) / (isamp*nblock*nmoves)
            uacae2(ist) = urunte2(ist) / (isamp*nblock*nmoves)
c
            uacad(ist) = uruntd(ist) / (isamp*nblock*nmoves)
            uacad2(ist) = uruntd2(ist) / (isamp*nblock*nmoves)
c
            uacar(ist) = uruntr(ist) / (isamp*nblock*nmoves)
            uacar2(ist) = uruntr2(ist) / (isamp*nblock*nmoves)
c
            uacap(ist) = uruntp(ist) / (isamp*nblock*nmoves)
            uacap2(ist) = uruntp2(ist) / (isamp*nblock*nmoves)
c
            uacan(ist) = uruntn(ist) / (isamp*nblock*nmoves)
            uacan2(ist) = uruntn2(ist) / (isamp*nblock*nmoves)
c
            uacanp(ist) = urtnp(ist) / (isamp*nblock*nmoves)
            uacanp2(ist) = urtnp2(ist) / (isamp*nblock*nmoves)
c
            do 8760, kst = 1, ist
              indxdf = ia(ist) + kst
              uacadf(indxdf) = urtdf(indxdf) / (isamp*nblock*nmoves)
              uacadf2(indxdf) = urtdf2(indxdf) / (isamp*nblock*nmoves)
 8760       continue
c
            do 8770, igr = 1, ngran
c       4-----
              uaaig(igr,ist) = urtig(igr,ist) / (isamp*nblock*nmoves)
              uaaig2(igr,ist) = urtig2(igr,ist) / (isamp*nblock*nmoves)
c
              uaaeg(igr,ist) = urteg(igr,ist) / (isamp*nblock*nmoves)
              uaaeg2(igr,ist) = urteg2(igr,ist) / (isamp*nblock*nmoves)
c
              uaadg(igr,ist) = urtdg(igr,ist) / (isamp*nblock*nmoves)
              uaadg2(igr,ist) = urtdg2(igr,ist) / (isamp*nblock*nmoves)
c
              uaarg(igr,ist) = urtrg(igr,ist) / (isamp*nblock*nmoves)
              uaarg2(igr,ist) = urtrg2(igr,ist) / (isamp*nblock*nmoves)
c       4-----
 8770       continue
c     3-----
 8800     continue
c
c    -----  print sample averages
c
          call mcblout(nstate,isamp,iblock,
     1        isamp*nblock*nmoves,iblacct,uacat,uacat2,uacaq,uacaq2,
     2        uacac,uacac2,uacace,uacace2,uacacd,uacacd2,
     3        uacacr,uacacr2,uacai,uacai2,uacae,uacae2,uacad,uacad2,
     4        uacar,uacar2,uacap,uacap2,uacan,uacan2,uacanp,uacanp2,
     5        uacadf,uacadf2,uaacg,uaacg2,uaaceg,uaaceg2,uaacdg,uaacdg2,
     3        uaacrg,uaacrg2,uaaig,uaaig2,uaaeg,uaaeg2,uaadg,uaadg2,
     4        uaarg,uaarg2,string2)
c   2-----
        endif
c 1-----  next sample
 9000 continue
c
c-----  end of run: calculate run averages
c
      utavc = uruntc / (nsamp*nblock*nmoves)
      utavc2 = uruntc2 / (nsamp*nblock*nmoves)
c
      utavce = urtce / (nsamp*nblock*nmoves)
      utavce2 = urtce2 / (nsamp*nblock*nmoves)
c
      utavcd = urtcd / (nsamp*nblock*nmoves)
      utavcd2 = urtcd2 / (nsamp*nblock*nmoves)
c
      utavcr = urtcr / (nsamp*nblock*nmoves)
      utavcr2 = urtcr2 / (nsamp*nblock*nmoves)
c
      do 9100, k = 1, ngrpair
c 1-----
        utavcg(k) = urtcg(k) / (nsamp*nblock*nmoves)
        utavcg2(k) = urtcg2(k) / (nsamp*nblock*nmoves)
c
        utavceg(k) = urtceg(k) / (nsamp*nblock*nmoves)
        utvceg2(k) = urtceg2(k) / (nsamp*nblock*nmoves)
c
        utavcdg(k) = urtcdg(k) / (nsamp*nblock*nmoves)
        utvcdg2(k) = urtcdg2(k) / (nsamp*nblock*nmoves)
c
        utavcrg(k) = urtcrg(k) / (nsamp*nblock*nmoves)
        utvcrg2(k) = urtcrg2(k) / (nsamp*nblock*nmoves)
c 1-----
 9100 continue
c
      do 9200, ist = 1, imcst
c 1-----
        utav(ist) = urunt(ist) / (nsamp*nblock*nmoves)
        utav2(ist) = urunt2(ist) / (nsamp*nblock*nmoves)
c
        utavq(ist) = uruntq(ist) / (nsamp*nblock*nmoves)
        utavq2(ist) = uruntq2(ist) / (nsamp*nblock*nmoves)
c
        utavi(ist) = urunti(ist) / (nsamp*nblock*nmoves)
        utavi2(ist) = urunti2(ist) / (nsamp*nblock*nmoves)
c
        utave(ist) = urunte(ist) / (nsamp*nblock*nmoves)
        utave2(ist) = urunte2(ist) / (nsamp*nblock*nmoves)
c
        utavd(ist) = uruntd(ist) / (nsamp*nblock*nmoves)
        utavd2(ist) = uruntd2(ist) / (nsamp*nblock*nmoves)
c
        utavr(ist) = uruntr(ist) / (nsamp*nblock*nmoves)
        utavr2(ist) = uruntr2(ist) / (nsamp*nblock*nmoves)
c
        utavp(ist) = uruntp(ist) / (nsamp*nblock*nmoves)
        utavp2(ist) = uruntp2(ist) / (nsamp*nblock*nmoves)
c
        utavn(ist) = uruntn(ist) / (nsamp*nblock*nmoves)
        utavn2(ist) = uruntn2(ist) / (nsamp*nblock*nmoves)
c
        utavnp(ist) = urtnp(ist) / (nsamp*nblock*nmoves)
        utavnp2(ist) = urtnp2(ist) / (nsamp*nblock*nmoves)
c
        do 9160, kst = 1, ist
          indxdf = ia(ist) + kst
          utavdf(indxdf) = urtdf(indxdf) / (nsamp*nblock*nmoves)
          utavdf2(indxdf) = urtdf2(indxdf) / (nsamp*nblock*nmoves)
 9160   continue
c
        do 9170, igr = 1, ngran
c   2-----
          utavig(igr,ist) = urtig(igr,ist) / (nsamp*nblock*nmoves)
          utavig2(igr,ist) = urtig2(igr,ist) / (nsamp*nblock*nmoves)
c
          utaveg(igr,ist) = urteg(igr,ist) / (nsamp*nblock*nmoves)
          utaveg2(igr,ist) = urteg2(igr,ist) / (nsamp*nblock*nmoves)
c
          utavdg(igr,ist) = urtdg(igr,ist) / (nsamp*nblock*nmoves)
          utavdg2(igr,ist) = urtdg2(igr,ist) / (nsamp*nblock*nmoves)
c
          utavrg(igr,ist) = urtrg(igr,ist) / (nsamp*nblock*nmoves)
          utavrg2(igr,ist) = urtrg2(igr,ist) / (nsamp*nblock*nmoves)
c   2-----
 9170   continue
c 1-----
 9200 continue
c
c-----  write lowest energy conformation
c
      write(iwr,9934) ulow
 9934 format(/,'  -- lowest energy conformation: utot = ',e20.10)
      write(iwr,9935)
 9935 format(/'  name',t21,'charge',t31,'x',t41,'y',t51,'z',
     1           t63,'alfa (b**3)','    radius     ispec'/)
c
      imp = 1
      img = 1
c
      do 9300, ii = 1, nxtpts
c 1-----
        if (ii .eq. mpol(imp)) then
c   2-----
c          write(iwr,9940) nxcent(ii),comma,chrg(ii),comma,
c     1    geoml(1,ii),comma,geoml(2,ii),comma,geoml(3,ii),comma,
c     2                   polar(imp)
 9940     format(a16,a1,4(f10.6,a1),f8.4)
          imp = imp + 1
        else
          write(iwr,9940) nxcent(ii),comma,chrg(ii),comma,
     1    geoml(1,ii),comma,geoml(2,ii),comma,geoml(3,ii),comma,
     2    alfext(ii)
c   2----
        endif
c
        if (ii .eq. igrlast(img)) then
          write (iwr,9950) grlabel(img)
 9950     format(a)
          img = img + 1
        endif
c 1-----
 9300 continue
c
      close(unit=incf,status='keep')
c
c     if (isolsav .eq. 1) then
c 1-----
c       call potsav(nwtc,itwoeps,nsamp*nblock*nmoves,maxneq,
c    2              field,xscm(1),xscm(1+nwtc),xscm(1+2*nwtc),
c    3              xscm(1+3*nwtc),xscm(1+4*nwtc))
c 1-----
c     endif
c
      return
      end
      subroutine mcstepp(nmovt,nstat,umove,umqm,umcl,umcle,umcld,umclr,
     1           umint,umelst,umdis,umrep,umpol,umneq,umneqp,udiff,
     2           umintg,umelstg,umdisg,umrepg,umclasg,umcleg,umcldg,
     3           umclrg,accept,xscm)
c
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      logical accept, tooclos
c
      dimension umove(nstat), umqm(nstat), umint(nstat), umelst(nstat),
     1 umdis(nstat), umrep(nstat), umpol(nstat), umneq(nstat),
     2 umneqp(nstat)
c
      dimension udiff(nstat*(nstat+1)/2)
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
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
      dimension xscm(*)
c
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      character*8 wfntyp
      common /wfnopt/ wfntyp
c
      character*8 runtyp
      common /runopt/ runtyp
c
      character*8 hndtyp
      common /hndopt/ hndtyp
c
c
cmw      include 'drf/dimpar'
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
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
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
       integer idafinf, navinf, iodainf
       common /infdaf/ idafinf,navinf,iodainf(2,mxneq+1)
c
       integer idafsta, navsta, iodasta
       common /stadaf/ idafsta,navsta,iodasta(2,mxneq+1)
c
       integer idafcst, navcst, iodacst
       common /cstdaf/ idafcst,navcst,iodacst(2,mxneq+1)
c
       integer idafind, navind, iodaind
       common /inddaf/ idafind,navind,iodaind(2,mxneq+1)
c
       integer idafino, navino, iodaino
       common /inodaf/ idafino,navino,iodaino(2,mxneq+1)
c
       integer idafpol, navpol, iodapol
       common /poldaf/ idafpol,navpol,iodapol(2,mxneq+1)
c
       integer idafdis, navdis, iodadis
       common /disdaf/ idafdis,navdis,iodadis(2,mxneq+1)
c
       integer idafrep, navrep, iodarep
       common /repdaf/ idafrep,navrep,iodarep(2,mxneq+1)
c
       integer idafrqm, navrqm, iodarqm
       common /rqmdaf/ idafrqm,navrqm,iodarqm(2,mxneq+1)
c
c
      character*80 textrf
      common /texneq/ textrf
c
      character*80 textneq
      common /neqtex/ textneq(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
      real*8 realneq
      common /bufb/ realneq(mxneq*10)
c
      real*8 unucrep, eoneel, ekin, enua, etwoel, extel, extnuc
      real*8 snucnuc, selel, snua, snucel, selnuc, stwoel, sextel
      real*8 sextnuc, selext, snucext, upolqm 
      real*8 uneqnuc, uneqel, uclas, suclas, uclase, uclasd, uclasr
      real*8 upolcl, udisp, rdisp, repmod, ustanuc, ustael
      common /rfene1/ unucrep,eoneel(mxst),ekin(mxst),enua(mxst),
     +  etwoel(mxst),extel(mxst),extnuc(mxst),snucnuc(mxst),
     +  selel(mxst),snua(mxst),snucel(mxst),selnuc(mxst),stwoel(mxst),
     +  sextel(mxst),
     +  sextnuc(mxst),selext(mxst),snucext(mxst),upolqm(mxst),
     +  uneqnuc(mxst),uneqel(mxst),uclas,suclas,uclase,uclasd,uclasr,
     +  upolcl,udisp(mxst),rdisp(mxst),repmod(mxst),
     +  ustanuc(mxst),ustael(mxst)
c
      real*8 smolnuc, smolel, smolmol, snucmol, selmol, sextmol
      real*8 smolext
      real*8 stotnuc, stotel, stotext, stotmol, upoleq, ucstst, ucstpl
      common/rfene2/smolnuc(mxst),smolel(mxst),smolmol(mxst),
     +  snucmol(mxst),selmol(mxst),sextmol(mxst),smolext(mxst),
     +  stotnuc(mxst),stotel(mxst),stotext(mxst),stotmol(mxst),
     +  upoleq(mxst),ucstst(mxst),ucstpl(mxst)
c
      real*8 uscf, uqm, uelst, uint, suqm, suint, uneq, uneqqm, uneqcl
      real*8 upolneq, usta, ustaqm, ustacl, ustaneq, uens
      common /rfene3/ uscf(mxst),uqm(mxst),uelst(mxst),uint(mxst),
     +  suqm(mxst),suint(mxst),uneq(mxst),uneqqm(mxst),uneqcl,
     +  upolneq(mxst),usta(mxst),ustaqm(mxst),ustacl,ustaneq(mxst),
     +  uens(mxst)
c
      real*8 snucno, selelo, snuao, snucelo, selnuco, stwoelo, sextelo
      real*8 sextno, selexto, snucexo, upolqmo, suclaso, upolclo
      common/rfene4/snucno(mxst),selelo(mxst),snuao(mxst),snucelo(mxst),
     +  selnuco(mxst),stwoelo(mxst),sextelo(mxst),sextno(mxst),
     +  selexto(mxst),snucexo(mxst),upolqmo(mxst),suclaso,upolclo
c
      real*8 smolno, smolelo, smolmo, snucmo, selmolo, sextmo, smolexo
      real*8 stotno, stotelo, stotexo, stotmo, upoleqo, ucstplo
      common/rfene5/smolno(mxst),smolelo(mxst),smolmo(mxst),
     +  snucmo(mxst),selmolo(mxst),sextmo(mxst),smolexo(mxst),
     +  stotno(mxst),stotelo(mxst),stotexo(mxst),stotmo(mxst),
     +  upoleqo(mxst),ucstplo(mxst)
c
      real*8 suqmo, suinto, stabtot, stabto 
      common/rfene6/suqmo(mxst),suinto(mxst),stabtot(mxst),stabto(mxst)
c
      real*8 extelg, extnucg, sextelg, sxtnucg, selextg, snucxtg, uclasg
      real*8 uclaseg, uclasdg, uclasrg, suclasg, upolclg, repmodg
      common /rfene7/ extelg(mxgran,mxst),extnucg(mxgran,mxst),
     +  sextelg(mxgran,mxst),sxtnucg(mxgran,mxst),selextg(mxgran,mxst),
     +  snucxtg(mxgran,mxst),uclasg(mxgrpar),uclaseg(mxgrpar),
     +  uclasdg(mxgrpar),uclasrg(mxgrpar),suclasg(mxgran,mxgran),
     +  upolclg(mxgran),repmodg(mxgran,mxst)
c
      real*8 sxtmolg, smolxtg, stotxtg
      common /rfene8/ sxtmolg(mxgran,mxst),smolxtg(mxgran,mxst),
     +  stotxtg(mxgran,mxst)
c
      real*8 uelstg, uintg, suintg, uneqclg, ustaclg
      common /rfene9/ uelstg(mxgran,mxst),uintg(mxgran,mxst),
     +  suintg(mxgran,mxst),uneqclg(mxgran),ustaclg(mxgran)
c
      real*8 sxtelog, sextnog, selxtog, sncexog, suclsog, uplclog
      common /rfene10/ sxtelog(mxgran,mxst),sextnog(mxgran,mxst),
     +  selxtog(mxgran,mxst),sncexog(mxgran,mxst),
     +  suclsog(mxgran,mxgran),uplclog(mxgran)
c
      real*8 sextmog, smlexog, sttexog
      common/rfene11/sextmog(mxgran,mxst),smlexog(mxgran,mxst),
     +  sttexog(mxgran,mxst)
c
      real*8 suintog
      common /rfene12/ suintog(mxgran,mxst)
c
c
      real*8 unucrepo, eoneelo, ekino, enuao, etwoelo
      real*8 extelo, extnuco
      real*8 uneqno, uneqelo, uclaso, uclaseo, uclasdo, uclasro
      real*8 rdispo, repmo, ustano, ustaelo
      common /rfene21/ unucrepo,eoneelo(mxst),ekino(mxst),enuao(mxst),
     +  etwoelo(mxst),extelo(mxst),extnuco(mxst),
     +  uneqno(mxst),uneqelo(mxst),uclaso,uclaseo,
     +  uclasdo,uclasro,
     +  rdispo(mxst),repmo(mxst),
     +  ustano(mxst),ustaelo(mxst)
c
      real*8 ucststo
      common/rfene22/
     +  ucststo(mxst)
c
      real*8 uscfo, uelsto, uinto, uneqo, uneqqmo, uneqclo
      real*8 upolno, ustao, ustaqmo, ustaclo
      common /rfene23/ uscfo(mxst),uelsto(mxst),uinto(mxst),
     +  uneqo(mxst),uneqqmo(mxst),uneqclo,
     +  upolno(mxst),ustao(mxst),ustaqmo(mxst),ustaclo
c
      real*8 extelgo, extnucgo, uclasgo
      real*8 uclasego, uclasdgo, uclasrgo, repmodgo
      common /rfene27/ extelgo(mxgran,mxst),extnucgo(mxgran,mxst),
     +  uclasgo(mxgrpar),uclasego(mxgrpar),
     +  uclasdgo(mxgrpar),uclasrgo(mxgrpar),
     +  repmodgo(mxgran,mxst)
c
      real*8 uelstgo, uintgo, uneqclgo, ustaclgo
      common /rfene29/ uelstgo(mxgran,mxst),uintgo(mxgran,mxst),
     +  uneqclgo(mxgran),ustaclgo(mxgran)
c
c
      integer neq2eps, momsav
      common /neqpar1/ neq2eps,momsav
c
      real*8 kapneq1, kapneq2, epsneq1, epsneq2
      common /neqpar2/ kapneq1,kapneq2,epsneq1,epsneq2
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
      real*8 ulow, geoml
      common /mcgeom/ ulow,geoml(3,mxpts)
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
caleko
      logical obeen, obeen2, obeen3, obeen4
      common/nottwi/obeen,obeen2,obeen3,obeen4
caleko
c
      dimension umclasg(mxgrpar), umcleg(mxgrpar), umcldg(mxgrpar),
     1 umclrg(mxgrpar)
      dimension umintg(mxgran,mxst), umelstg(mxgran,mxst),
     1 umdisg(mxgran,mxst), umrepg(mxgran,mxst)
c
c-----  local declarations
c
c
      dimension alfa(3,3), alfan(3,3), r(3,3)
      dimension thspoln(6)
c
c      dimension clas(ngran*(ngran+1)*5 + 1)
c      dimension polat(mxnpts), polatt(6,mxnpts)
      common/scrtch/clas(3000),polat(3000),polatt(6,3000)
c
      dimension kpol(mxnpts)
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
      data deltol /1.0d0/
c
      obeen = .false.
      obeen4 = .false.
      ixrelay = 1
      ixomga = 1
      ixwtvr = 1
      ixamat = 1
      if (ibem .ne. 0) ixbmat = 1
c
      if (ngran*(ngran+1)*5 + 1.gt.3000) call caserr
     + ('mcsub.f: mcstepp: increase dimension of clas >3000')
c     check below redundant as mxnpts is hard-wired to 20 by parameter
c     statement in comdrf/sizesrf
c     if (mxnpts.gt.3000) call caserr
c    + ('mcstepp: mxnpts>3000 increase dimension of polat, polatt')
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
      ist = iactst
c
      ngrpair = ngran*(ngran+1)/2
c
      accept = .false.
      tooclos = .false.
c
c-----  perform random rotation (and translation)
c       in monte carlo run,
c       or controlled movement in classical pes run
c
c-----  select group to be moved
c
      if (nopes .eq. 1) then
c 1-----
        mcgrp = 1 + mod(nmovt,ngrpmc)
        igrp = ibitmc(mcgrp)
      else
        igrp = npesgrp
c 1-----
      endif
c
c-----  collect coordinates of group to be moved
c
      istrt = igrpst(igrp)
      nnpts = ngrppt(igrp) + 1
      kgrpol = igrpol(igrp)
      inpt(1) = kgrpol
c
      xnpts(1,1) = xpts(1,istrt)
      xnpts(2,1) = xpts(2,istrt)
      xnpts(3,1) = xpts(3,istrt)
c
      do 100, ipt = 2, nnpts
c
c  -----  set pointer to original coordinates that are now to be modifie
c
        inpt(ipt) = istrt - 1 + ipt
        do 100, k = 1, 3
          xnpts(k,ipt) = xpts(k,istrt-1+ipt)
  100 continue
c
c-----  get polarisability tensor components of group
c
      do 150, kl = 1, 6
        thispol(kl) = grpol(kl,kgrpol)
  150 continue
c
      if (nopes .eq. 1) then
c 1-----
        if (norot .eq. 0) then
c   2-----
c    -----  perform rotation of group 'round centre of polarisability
c           along all 3 cartesian axes
c
c    -----  rotation around z-axis
c
          rot = darot * (ranget() - 0.5d0)
c
          cc = cos(rot)
          ss = sin(rot)
c
          c2 = cc**2
          s2 = ss**2
          ccss = cc*ss
c
          do 221, ipt = 2, nnpts
            gx = xnpts(1,ipt) - xpts(1,istrt)
            gy = xnpts(2,ipt) - xpts(2,istrt)
c
            xnpts(1,ipt) = cc*gx + ss*gy + xpts(1,istrt)
            xnpts(2,ipt) = cc*gy - ss*gx + xpts(2,istrt)
c
  221     continue
c
c    -----  perform rotation of polarisability tensor
c           'round centre of polarisability
c
          thspoln(1) = c2*thispol(1) + s2*thispol(3)
     2               + 2.0d0*ccss*thispol(2)
          thspoln(2) = (c2-s2)*thispol(2)
     2               + ccss*(thispol(3)-thispol(1))
          thspoln(3) = c2*thispol(3) + s2*thispol(1)
     2               - 2.0d0*ccss*thispol(2)
          thspoln(4) = cc*thispol(4) + ss*thispol(5)
          thspoln(5) = cc*thispol(5) - ss*thispol(4)
          thspoln(6) = thispol(6)
c
c         call copyv (thspoln,thispol,6)
          call dcopy(6,thspoln,1,thispol,1)
c
c    -----  rotation around x-axis
c
          rot = darot * (ranget() - 0.5d0)
c
          cc = cos(rot)
          ss = sin(rot)
c
          c2 = cc**2
          s2 = ss**2
          ccss = cc*ss
          do 222, ipt = 2, nnpts
            gy = xnpts(2,ipt) - xpts(2,istrt)
            gz = xnpts(3,ipt) - xpts(3,istrt)
c
            xnpts(2,ipt) = cc*gy + ss*gz + xpts(2,istrt)
            xnpts(3,ipt) = cc*gz - ss*gy + xpts(3,istrt)
c
  222     continue
c
c    -----  perform rotation of polarisability tensor
c           'round centre of polarisability
c
          thspoln(3) = c2*thispol(3) + s2*thispol(6)
     2               + 2.0d0*ccss*thispol(5)
          thspoln(5) = (c2-s2)*thispol(5)
     2               + ccss*(thispol(6)-thispol(3))
          thspoln(6) = c2*thispol(6) + s2*thispol(3)
     2               - 2.0d0*ccss*thispol(5)
          thspoln(2) = cc*thispol(2) + ss*thispol(4)
          thspoln(4) = cc*thispol(4) - ss*thispol(2)
          thspoln(1) = thispol(1)
c
c         call copyv (thspoln,thispol,6)
          call dcopy(6,thspoln,1,thispol,1)
c
c    -----  rotation around y-axis
c
c
          rot = darot * (ranget() - 0.5d0)
c
          cc = cos(rot)
          ss = sin(rot)
c
          c2 = cc**2
          s2 = ss**2
          ccss = cc*ss
          do 223, ipt = 2, nnpts
            gz = xnpts(3,ipt) - xpts(3,istrt)
            gx = xnpts(1,ipt) - xpts(1,istrt)
c
            xnpts(3,ipt) = cc*gz - ss*gx + xpts(3,istrt)
            xnpts(1,ipt) = cc*gx + ss*gz + xpts(1,istrt)
c
  223     continue
c
c    -----  perform rotation of polarisability tensor
c           'round centre of polarisability
c
          thspoln(1) = c2*thispol(1) + s2*thispol(6)
     2               + 2.0d0*ccss*thispol(4)
          thspoln(4) = (c2-s2)*thispol(4)
     2               + ccss*(thispol(6)-thispol(1))
          thspoln(6) = c2*thispol(6) + s2*thispol(1)
     2               - 2.0d0*ccss*thispol(4)
          thspoln(2) = cc*thispol(2) + ss*thispol(5)
          thspoln(5) = cc*thispol(5) - ss*thispol(2)
          thspoln(3) = thispol(3)
c
c    -----   store updated polarisability tensor
c
c         call copyv (thspoln,thispol,6)
          call dcopy(6,thspoln,1,thispol,1)
c   2-----
        endif
c
c  -----  perform translation as well if required
c
        if (notrans .eq. 0) then
c   2-----
c
c    -----  calculate the distance through which to move
c
          trans = dtrans * (ranget() - 0.5d0)
c
c    -----  translation along z-axis
c
          do 521, ipt = 1, nnpts
            xnpts(3,ipt) = xnpts(3,ipt) + trans
  521     continue
c
c    -----  translation along x-axis
c
          trans = dtrans * (ranget() - 0.5d0)
c
          do 522, ipt = 1, nnpts
            xnpts(1,ipt) = xnpts(1,ipt) + trans
  522     continue
c
c    -----  translation along y-axis
c
          trans = dtrans * (ranget() - 0.5d0)
c
          do 523, ipt = 1, nnpts
            xnpts(2,ipt) = xnpts(2,ipt) + trans
  523     continue
c
c    -----  if a boundary is present, check whether the group has not mo
c           too close to it
c
          if (ibem .ne. 0) then
c     3-----
            call checdis(nnpts,xnpts,tooclos,xscm)
            if (tooclos) then
c       4-----
c        -----  the group has been moved too close to the boundary
c               and is rejected, attempt next move
c
              accept = .false.
              umove(iactst) = uold(iactst)
              umqm(iactst) = eqmol(iactst)
              umcl = eclasol
              umcle = ecleol
              umcld = ecldol
              umclr = eclrol
c
              do 500, k = 1, ngrpair
c         5-----
                umclasg(k) = eclsolg(k)
                umcleg(k) = ecleolg(k)
                umcldg(k) = ecldolg(k)
                umclrg(k) = eclrolg(k)
c         5-----
  500         continue
c
              umint(iactst) = eintol(iactst)
              umelst(iactst) = eelstol(iactst)
              umdis(iactst) = edisol(iactst)
              umrep(iactst) = erepol(iactst)
              umpol(iactst) = epolol(iactst)
              umneq(iactst) = eneqol(iactst)
              umneqp(iactst) = eneqpol(iactst)
              extnuc(iactst) = extnol(iactst)
              repmod(iactst) = erepol(iactst)
c
              do 600, igr = 1, ngran
c         5-----
                umintg(igr,iactst) = eintolg(igr,iactst)
                umelstg(igr,iactst) = elstolg(igr,iactst)
                umdisg(igr,iactst) = edisolg(igr,iactst)
                umrepg(igr,iactst) = erepolg(igr,iactst)
                extnucg(igr,iactst) = extnolg(igr,iactst)
                repmodg(igr,iactst) = erepolg(igr,iactst)
c         5-----
  600         continue
c
              do 700, ist = 2, imcst
c         5-----
                umove(ist) = uold(ist)
                umqm(ist) = eqmol(ist)
                umint(ist) = eintol(ist)
                umelst(ist) = eelstol(ist)
                umdis(ist) = edisol(ist)
                umrep(ist) = erepol(ist)
                umpol(ist) = epolol(ist)
                umneq(ist) = eneqol(ist)
                umneqp(ist) = eneqpol(ist)
c
                do 620, kst = 1, ist
                  indxdf = ia(ist) + kst
                  udiff(indxdf) = edifol(indxdf)
  620           continue
c
                extnuc(ist) = extnol(ist)
                repmod(ist) = erepol(ist)
c
                do 640, igr = 1, ngran
c           6-----
                  umintg(igr,ist) = eintolg(igr,ist)
                  umelstg(igr,ist) = elstolg(igr,ist)
                  umdisg(igr,ist) = edisolg(igr,ist)
                  umrepg(igr,ist) = erepolg(igr,ist)
                  extnucg(igr,ist) = extnolg(igr,ist)
                  repmodg(igr,ist) = erepolg(igr,ist)
c           6-----
  640           continue
c         5-----
  700         continue
c
              return
c       4-----
            endif
c     3-----
          endif
c
        if ((ieffpol .eq. 1) .and. (isodis .eq. 0)) then
c   2-----
c    -----  make the group polarisability matrix -amat- block diagonal
c           to get effective atom polarisabilities for use in
c           calculation of dispersion energy
c
          np = nnpts - 1
          np3 = 3*np
          nnp3 = np3*(np3+1)/2
c
c    -----  initialize a correct array kpol, denoting the
c           original polarizabilities of the group
c
          do 800, i = 2, nnpts
            kpol(i-1) = i
  800     continue
c
          i10 = igmem_alloc(nnp3)
          call clear (xscm(i10),nnp3)
          call drfamat
     1   (xscm(i10),1,.false.,np,xnpts,kpol,alfext(istrt+1),
     2    idrfout,afact,ithole)
cc
c  -----  invert -a-
c
c     subroutine linv3p can be replaced by standard 'vectorized'
c     matrix inversion subroutine (nag,linpack etc)
c
          call linv3p(xscm(i10),auxx,1,np3,ier)
          if (idrfout .ge. 3)
     1      call hatout(xscm(i10),np3,np3,3,'group-amat')
c
          i20 = igmem_alloc(nnp3)
          call effpol(np,np3,ia,xscm(i10),xscm(i20),
     1            polat,polatt)
c
          do 810, i = 2, nnpts
c     3-----
            atpol(istrt+i-1) = polat(i-1)
            do 820, k = 1, 6
              atpolt(k,istrt+i-1) = polatt(k,i-1)
  820       continue
c     3-----
  810     continue
          call gmem_free(i20)
          call gmem_free(i10)
c   2-----
        endif
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  get coordinates etc... for pes scan of classical group
c
c  -----  xscm work array used for reading in of coordinates,
c         and polarisability coupling matrix (if required)
c
        npt3 = 3*(nnpts - 1)
        lcor = (npt3*(npt3+1))/2
        ixpol = max(lpes,lcor) + 1
c
        i10 = igmem_alloc(max(lpes,lcor))
        i20 = igmem_alloc(nnpts)
        call getpesc(ipesout,nmovt+1,lpes,lpesc,intpesc,
     2               ieffpol,isodis,ithole,afact,istrt,nnpts,lcor,
     3               xnpts,thispol,xscm(i10),xscm(i20))
        call gmem_free(i20)
        call gmem_free(i10)
c 1-----
      endif
c
c-----  calculate energy of attempted move
c
c-----  set wave function data in order
c
      gamdrf = gammc(ist)
c
c-----  read wavefunction data for actual state
c
      i10 = igmem_alloc(l3)
      call daread(idafdrf,iodadrf,xscm(i10),l2,182)
      call dawrit(idafh,ioda,xscm(i10),l2,12,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l2,192)
      call dawrit(idafh,ioda,xscm(i10),l2,13,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l3,202)
      call dawrit(idafh,ioda,xscm(i10),l3,15,navh)
c
      call daread(idafdrf,iodadrf,xscm(i10),l2,212)
      call dawrit(idafh,ioda,xscm(i10),l2,16,navh)
c
      if (scftyp .eq. 'uhf') then
c 1-----
        call daread(idafdrf,iodadrf,xscm(i10),l3,222)
        call dawrit(idafh,ioda,xscm(i10),l3,19,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,232)
        call dawrit(idafh,ioda,xscm(i10),l2,20,navh)
c 1-----
      endif
c
      call gmem_free(i10)
c
cnot      call energx
      call hfscf(xscm)
      nprint = -5
c
c-----  perform analysis of updated energy if converged
c
      if (.not. excess) then
c 1-----
        ieps = 0
        call arfanal(ieps,unucrep,
     1 eoneel(ist),ekin(ist),enua(ist),etwoel(ist),
     1 uqm(ist),uscf(ist),snucnuc(ist),selel(ist),snua(ist),stwoel(ist),
     2 smolnuc(ist),smolel(ist),snucmol(ist),selmol(ist),smolmol(ist),
     3 suqm(ist),upolqm(ist),uneqnuc(ist),uneqel(ist),uneqqm(ist),
     4 ustanuc(ist),ustael(ist),ustaqm(ist),uclase,uclasd,uclasr,uclas,
     5 suclas,upolcl,uneqcl,ustacl,extnuc(ist),extel(ist),sextnuc(ist),
     6 sextel(ist),sextmol(ist),selext(ist),snucext(ist),smolext(ist),
     7 stotnuc(ist),stotel(ist),stotmol(ist),stotext(ist),stabtot(ist),
     8 uelst(ist),suint(ist),uint(ist),udisp(ist),rdisp(ist),
     9 repmod(ist),upoleq(ist),ucstst(ist),ucstpl(ist),uneq(ist),
     1 usta(ist),upolneq(ist),ustaneq(ist),uens(ist),uclasg,uclaseg,
     2 uclasdg,uclasrg,suclasg,upolclg,uneqclg,ustaclg,extnucg(1,ist),
     3 extelg(1,ist),sxtnucg(1,ist),sextelg(1,ist),sxtmolg(1,ist),
     4 selextg(1,ist),snucxtg(1,ist),smolxtg(1,ist),stotxtg(1,ist),
     5 uelstg(1,ist),suintg(1,ist),uintg(1,ist),repmodg(1,ist),xscm)
c
        if (itwoeps .eq. 1) then
c   2-----
          ieps = 1
        call arfanal(ieps,unucrepo,
     1 eoneelo(ist),ekino(ist),enuao(ist),etwoelo(ist),
     1 uqm(ist),uscfo(ist),snucno(ist),selel(ist),snuao(ist),
     1 stwoelo(ist),
     2 smolno(ist),smolelo(ist),snucmo(ist),selmolo(ist),smolmo(ist),
     3 suqmo(ist),upolqmo(ist),uneqno(ist),uneqelo(ist),uneqqmo(ist),
     4 ustano(ist),ustaelo(ist),ustaqmo(ist),uclaseo,uclasdo,
     4 uclasro,uclaso,
     5 suclaso,upolclo,uneqclo,ustaclo,extnuco(ist),extelo(ist),
     5 sextno(ist),
     6 sextelo(ist),sextmo(ist),selexto(ist),snucexo(ist),smolexo(ist),
     7 stotno(ist),stotelo(ist),stotmo(ist),stotexo(ist),stabto(ist),
     8 uelsto(ist),suinto(ist),uinto(ist),udisp(ist),rdispo(ist),
     9 repmo(ist),upoleqo(ist),ucststo(ist),ucstplo(ist),uneqo(ist),
     1 ustao(ist),upolno(ist),ustano(ist),uens(ist),uclasgo,uclasego,
     2 uclasdgo,uclasrgo,suclsog,uplclog,uneqclgo,ustaclgo,
     2 extnucgo(1,ist),
     3 extelgo(1,ist),sextnog(1,ist),sxtelog(1,ist),sextmog(1,ist),
     4 selxtog(1,ist),sncexog(1,ist),smlexog(1,ist),sttexog(1,ist),
     5 uelstgo(1,ist),suintog(1,ist),uintgo(1,ist),repmodgo(1,ist),
     5 xscm)
c   2-----
        endif
c
        call enanal
c
c  -----  copy wave function data
c
        i10 = igmem_alloc(l3)
        call daread(idafh,ioda,xscm(i10),l2,12)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,182,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,13)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,192,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l3,15)
        call dawrit(idafdrf,iodadrf,xscm(i10),l3,202,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,16)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,212,navdrf)
c
        if (scftyp .eq. 'uhf') then
c   2-----
          call daread(idafh,ioda,xscm(i10),l3,19)
          call dawrit(idafdrf,iodadrf,xscm(i10),l3,222,navdrf)
c
          call daread(idafh,ioda,xscm(i10),l2,20)
          call dawrit(idafdrf,iodadrf,xscm(i10),l2,232,navdrf)
c   2-----
        endif
c
      call gmem_free(i10)
c
        if ((nopes .eq. 0) .and. (ipesout .ge. 1)) call drfout
c
c  -----  accept or reject attempted move
c
        umove(iactst) = uens(iactst)
c
        delta = umove(iactst) - uold(iactst)
c 1-----
      else
c 1-----
        delta = 1.1d0*deltol
c 1-----
      endif
c
c-----  ensure acceptance for pes scan
c
      if (nopes .eq. 0) delta = -abs(delta)
c
      if ((delta .gt. zero) .or. (delta .lt. delmax)) then
c 1-----
        if ((delta .gt. deltol) .or. (exp(onekt*delta) .lt. ranget())
     2   .or. (umove(iactst) .lt. enmin) .or. (delta .lt. delmax)) then
c   2-----
c    -----  move is rejected, attempt next move
c
          accept = .false.
          umove(iactst) = uold(iactst)
          umqm(iactst) = eqmol(iactst)
          umcl = eclasol
          umcle = ecleol
          umcld = ecldol
          umclr = eclrol
c
          do 1000, k = 1, ngrpair
c     3-----
            umclasg(k) = eclsolg(k)
            umcleg(k) = ecleolg(k)
            umcldg(k) = ecldolg(k)
            umclrg(k) = eclrolg(k)
c     3-----
 1000     continue
c
          umint(iactst) = eintol(iactst)
          umelst(iactst) = eelstol(iactst)
          umdis(iactst) = edisol(iactst)
          umrep(iactst) = erepol(iactst)
          umpol(iactst) = epolol(iactst)
          umneq(iactst) = eneqol(iactst)
          umneqp(iactst) = eneqpol(iactst)
          extnuc(iactst) = extnol(iactst)
          repmod(iactst) = erepol(iactst)
c
          do 1100, igr = 1, ngran
c     3-----
            umintg(igr,iactst) = eintolg(igr,iactst)
            umelstg(igr,iactst) = elstolg(igr,iactst)
            umdisg(igr,iactst) = edisolg(igr,iactst)
            umrepg(igr,iactst) = erepolg(igr,iactst)
            extnucg(igr,iactst) = extnolg(igr,iactst)
            repmodg(igr,iactst) = erepolg(igr,iactst)
c     3-----
 1100     continue
c
          do 1200, ist = 2, imcst
c     3-----
            umove(ist) = uold(ist)
            umqm(ist) = eqmol(ist)
            umint(ist) = eintol(ist)
            umelst(ist) = eelstol(ist)
            umdis(ist) = edisol(ist)
            umrep(ist) = erepol(ist)
            umpol(ist) = epolol(ist)
            umneq(ist) = eneqol(ist)
            umneqp(ist) = eneqpol(ist)
c
            do 1120, kst = 1, ist
              indxdf = ia(ist) + kst
              udiff(indxdf) = edifol(indxdf)
 1120       continue
c
            extnuc(ist) = extnol(ist)
            repmod(ist) = erepol(ist)
c
            do 1140, igr = 1, ngran
c       4-----
              umintg(igr,ist) = eintolg(igr,ist)
              umelstg(igr,ist) = elstolg(igr,ist)
              umdisg(igr,ist) = edisolg(igr,ist)
              umrepg(igr,ist) = erepolg(igr,ist)
              extnucg(igr,ist) = extnolg(igr,ist)
              repmodg(igr,ist) = erepolg(igr,ist)
c       4-----
 1140       continue
c     3-----
 1200     continue
c
          return
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  total energy is lower than reference energy
c         (or acceptance of the configuration is forced)
c
c        if (umove(iactst) .lt. ulow) then
c   2-----
c    -----  save new lowest energy and geometry
c
c          ulow = umove(iactst)
c          do 2100, k = 1, 3
c            geoml(k,istrt) = xnpts(k,1)
c 2100     continue
c          do 2200, in = 2, nnpts
c            do 2200, k = 1, 3
c              geoml(k,inpt(in)) = xnpts(k,in)
c 2200     continue
c   2-----
c        endif
c 1-----
      endif
c
c-----  move is accepted, store updated relay, wt, vr
c                         zfa, dzfa and dpseu
c
      accept = .true.
c
      umqm(iactst) = eqm(iactst)
      uold(iactst) = umove(iactst)
      eqmol(iactst) = umqm(iactst)
c
      umcl = eclas
      umcle = eclase
      umcld = eclasd
      umclr = eclasr
c
      eclasol = umcl
      ecleol = umcle
      ecldol = umcld
      eclrol = umclr
c
      do 3100, k = 1, ngrpair
c 1-----
        umclasg(k) = eclasg(k)
        umcleg(k) = eclaseg(k)
        umcldg(k) = eclasdg(k)
        umclrg(k) = eclasrg(k)
c
        eclsolg(k) = umclasg(k)
        ecleolg(k) = umcleg(k)
        ecldolg(k) = umcldg(k)
        eclrolg(k) = umclrg(k)
c 1-----
 3100 continue
c
      umint(iactst) = eint(iactst)
      umelst(iactst) = eelst(iactst)
      umdis(iactst) = edisp(iactst)
      umrep(iactst) = erep(iactst)
      umpol(iactst) = epol(iactst)
      umneq(iactst) = eneq(iactst)
      umneqp(iactst) = eneqp(iactst)
c
      eintol(iactst) = umint(iactst)
      eelstol(iactst) = umelst(iactst)
      edisol(iactst) = umdis(iactst)
      extnol(iactst) = extnuc(iactst)
      erepol(iactst) = umrep(iactst)
      epolol(iactst) = umpol(iactst)
      eneqol(iactst) = umneq(iactst)
      eneqpol(iactst) = umneqp(iactst)
c
      do 3200, igr = 1, ngran
c 1-----
        umintg(igr,iactst) = eintg(igr,iactst)
        umelstg(igr,iactst) = eelstg(igr,iactst)
        umdisg(igr,iactst) = edispg(igr,iactst)
        umrepg(igr,iactst) = erepg(igr,iactst)
c
        eintolg(igr,iactst) = umintg(igr,iactst)
        elstolg(igr,iactst) = umelstg(igr,iactst)
        edisolg(igr,iactst) = umdisg(igr,iactst)
        extnolg(igr,iactst) = extnucg(igr,iactst)
        erepolg(igr,iactst) = umrepg(igr,iactst)
c 1-----
 3200 continue
c
c-----  calculate second state energies
c
      do 4000, ist = 2, imcst
c 1-----
c  -----  read wavefunction data for actual state
c
        i10 = igmem_alloc(l3)
        call daread(idafdrf,iodadrf,xscm(i10),l2,181+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,12,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,191+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,13,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l3,201+ist)
        call dawrit(idafh,ioda,xscm(i10),l3,15,navh)
c
        call daread(idafdrf,iodadrf,xscm(i10),l2,211+ist)
        call dawrit(idafh,ioda,xscm(i10),l2,16,navh)
c
        if (scftyp .eq. 'uhf') then
c   2-----
          call daread(idafdrf,iodadrf,xscm(i10),l3,221+ist)
          call dawrit(idafh,ioda,xscm(i10),l3,19,navh)
c
          call daread(idafdrf,iodadrf,xscm(i10),l2,231+ist)
          call dawrit(idafh,ioda,xscm(i10),l2,20,navh)
c   2-----
        endif
c
      call gmem_free(i10)
c
        secstat = .true.
        iactst = ist
c
        gamdrf = gammc(ist)
c
c  -----  calculate energy of second state
c
cnot        call energx
        call hfscf(xscm)
c
c  -----  perform analysis of energy
c
        ieps = 0
        call arfanal(ieps,unucrep,
     1 eoneel(ist),ekin(ist),enua(ist),etwoel(ist),
     1 uqm(ist),uscf(ist),snucnuc(ist),selel(ist),snua(ist),stwoel(ist),
     2 smolnuc(ist),smolel(ist),snucmol(ist),selmol(ist),smolmol(ist),
     3 suqm(ist),upolqm(ist),uneqnuc(ist),uneqel(ist),uneqqm(ist),
     4 ustanuc(ist),ustael(ist),ustaqm(ist),uclase,uclasd,uclasr,uclas,
     5 suclas,upolcl,uneqcl,ustacl,extnuc(ist),extel(ist),sextnuc(ist),
     6 sextel(ist),sextmol(ist),selext(ist),snucext(ist),smolext(ist),
     7 stotnuc(ist),stotel(ist),stotmol(ist),stotext(ist),stabtot(ist),
     8 uelst(ist),suint(ist),uint(ist),udisp(ist),rdisp(ist),
     9 repmod(ist),upoleq(ist),ucstst(ist),ucstpl(ist),uneq(ist),
     1 usta(ist),upolneq(ist),ustaneq(ist),uens(ist),uclasg,uclaseg,
     2 uclasdg,uclasrg,suclasg,upolclg,uneqclg,ustaclg,extnucg(1,ist),
     3 extelg(1,ist),sxtnucg(1,ist),sextelg(1,ist),sxtmolg(1,ist),
     4 selextg(1,ist),snucxtg(1,ist),smolxtg(1,ist),stotxtg(1,ist),
     5 uelstg(1,ist),suintg(1,ist),uintg(1,ist),repmodg(1,ist),xscm)
c
        if (itwoeps .eq. 1) then
c   2-----
          ieps = 1
        call arfanal(ieps,unucrepo,
     1 eoneelo(ist),ekino(ist),enuao(ist),etwoelo(ist),
     1 uqm(ist),uscfo(ist),snucno(ist),selel(ist),snuao(ist),
     1 stwoelo(ist),
     2 smolno(ist),smolelo(ist),snucmo(ist),selmolo(ist),smolmo(ist),
     3 suqmo(ist),upolqmo(ist),uneqno(ist),uneqelo(ist),uneqqmo(ist),
     4 ustano(ist),ustaelo(ist),ustaqmo(ist),uclaseo,uclasdo,
     4 uclasro,uclaso,
     5 suclaso,upolclo,uneqclo,ustaclo,extnuco(ist),extelo(ist),
     5 sextno(ist),
     6 sextelo(ist),sextmo(ist),selexto(ist),snucexo(ist),smolexo(ist),
     7 stotno(ist),stotelo(ist),stotmo(ist),stotexo(ist),stabto(ist),
     8 uelsto(ist),suinto(ist),uinto(ist),udisp(ist),rdispo(ist),
     9 repmo(ist),upoleqo(ist),ucststo(ist),ucstplo(ist),uneqo(ist),
     1 ustao(ist),upolno(ist),ustano(ist),uens(ist),uclasgo,uclasego,
     2 uclasdgo,uclasrgo,suclsog,uplclog,uneqclgo,ustaclgo,
     2 extnucgo(1,ist),
     3 extelgo(1,ist),sextnog(1,ist),sxtelog(1,ist),sextmog(1,ist),
     4 selxtog(1,ist),sncexog(1,ist),smlexog(1,ist),sttexog(1,ist),
     5 uelstgo(1,ist),suintog(1,ist),uintgo(1,ist),repmodgo(1,ist),
     5 xscm)
c   2-----
        endif
c
        call enanal
c
        umove(ist) = uens(ist)
        umqm(ist) = eqm(ist)
        umint(ist) = eint(ist)
        umelst(ist) = eelst(ist)
        umdis(ist) = edisp(ist)
        umrep(ist) = erep(ist)
        umpol(ist) = epol(ist)
        umneq(ist) = eneq(ist)
        umneqp(ist) = eneqp(ist)
        do 4100, kst = 1, ist
          indxdf = ia(ist) + kst
          udiff(indxdf) = uens(kst) - uens(ist)
          edifol(indxdf) = udiff(indxdf)
 4100   continue
c
        uold(ist) = umove(ist)
        eqmol(ist) = umqm(ist)
        eintol(ist) = umint(ist)
        eelstol(ist) = umelst(ist)
        edisol(ist) = umdis(ist)
        extnol(ist) = extnuc(ist)
        erepol(ist) = umrep(ist)
        epolol(ist) = umpol(ist)
        eneqol(ist) = umneq(ist)
        eneqpol(ist) = umneqp(ist)
c
        do 4200, igr = 1, ngran
c   2-----
          umintg(igr,ist) = eintg(igr,ist)
          umelstg(igr,ist) = eelstg(igr,ist)
          umdisg(igr,ist) = edispg(igr,ist)
          umrepg(igr,ist) = erepg(igr,ist)
          eintolg(igr,ist) = umintg(igr,ist)
          elstolg(igr,ist) = umelstg(igr,ist)
          edisolg(igr,ist) = umdisg(igr,ist)
          extnolg(igr,ist) = extnucg(igr,ist)
          erepolg(igr,ist) = umrepg(igr,ist)
c   2-----
 4200   continue
c
c  -----  write wavefunction data for actual state
c
        i10 = igmem_alloc(l3)
        call daread(idafh,ioda,xscm(i10),l2,12)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,181+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,13)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,191+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l3,15)
        call dawrit(idafdrf,iodadrf,xscm(i10),l3,201+ist,navdrf)
c
        call daread(idafh,ioda,xscm(i10),l2,16)
        call dawrit(idafdrf,iodadrf,xscm(i10),l2,211+ist,navdrf)
c
        if (scftyp .eq. 'uhf') then
c   2-----
          call daread(idafh,ioda,xscm(i10),l3,19)
          call dawrit(idafdrf,iodadrf,xscm(i10),l3,221+ist,navdrf)
c
          call daread(idafh,ioda,xscm(i10),l2,20)
          call dawrit(idafdrf,iodadrf,xscm(i10),l2,231+ist,navdrf)
c   2-----
        endif
        call gmem_free(i10)
c
        if ((nopes .eq. 0) .and. (ipesout .ge. 1)) call drfout
c 1-----
 4000 continue
c
      iactst = 1
      secstat =.false.
c
c-----  store new geometry, to be modified in the next move
c
      do 5000, k = 1, 3
        do 5000, ipt = 1, nnpts
          xpts(k,istrt-1+ipt) = xnpts(k,ipt)
 5000 continue
c
c-----  store new polarizability tensor components
c
      do 5100, kl = 1, 6
        grpol(kl,kgrpol) = thispol(kl)
 5100 continue
c
        if (umove(iactst) .lt. ulow) then
c   2-----
c    -----  save new lowest energy and geometry
c
          ulow = umove(iactst)
          do 2200, k = 1, 3
          do 2200, ipt = 1, nxtpts
              geoml(k,ipt) = xpts(k,ipt)
 2200     continue
c   2-----
        endif
c
c-----  debug printing
c
      if (imcout .eq. 5) then
        do 5200, i = 1, nxtpts
          write(iwr,1) i, (xpts(k,i), k = 1, 3)
    1     format(i4,3f12.6)
 5200   continue
      endif
c
c-----  store updated matrices for next move
c
c-----  classical energies
c
      call daread(idafdrf,iodadrf,clas,ngran*(ngran+1)*5+1,91)
      call dawrit(idafdrf,iodadrf,clas,ngran*(ngran+1)*5+1,21,navdrf)
c
c-----  relay
c
      i10 = igmem_alloc(ndim*ndim)
      call daread(idafdrf,iodadrf,xscm(i10),ndim*ndim,42)
      call dawrit(idafdrf,iodadrf,xscm(i10),ndim*ndim,41,navdrf)
      call gmem_free(i10)
c
c-----  wt
c
      i10 = igmem_alloc(nwtr*nwtc)
      call daread(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,27)
      call dawrit(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,26,navdrf)
      call gmem_free(i10)
c
c-----  vr
c
      i10 = igmem_alloc(nwtr*nwtc)
      call daread(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,29)
      call dawrit(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,28,navdrf)
      call gmem_free(i10)
c
c-----  zfa
c
      i10 = igmem_alloc(4*nzfa)
      call daread(idafdrf,iodadrf,xscm(i10),4*nzfa,76)
      call dawrit(idafdrf,iodadrf,xscm(i10),4*nzfa,3,navdrf)
      call gmem_free(i10)
c
c-----  zpa
c
      i10 = igmem_alloc(nzfa)
      call daread(idafdrf,iodadrf,xscm(i10),nzfa,77)
      call dawrit(idafdrf,iodadrf,xscm(i10),nzfa,4,navdrf)
      call gmem_free(i10)
c
c-----  dzfa
c
      i10 = igmem_alloc(12*nzfa)
      call daread(idafdrf,iodadrf,xscm(i10),12*nzfa,78)
      call dawrit(idafdrf,iodadrf,xscm(i10),12*nzfa,5,navdrf)
      call gmem_free(i10)
c
c-----  dpseu
c
      i10 = igmem_alloc(3*nzfa)
      call daread(idafdrf,iodadrf,xscm(i10),3*nzfa,79)
      call dawrit(idafdrf,iodadrf,xscm(i10),3*nzfa,6,navdrf)
      call gmem_free(i10)
c
      if (itwoeps .eq. 1) then
c
c  -----  relay
c
        i10 = igmem_alloc(ndim*ndim)
        call daread(idafdrf,iodadrf,xscm(i10),ndim*ndim,46)
        call dawrit(idafdrf,iodadrf,xscm(i10),ndim*ndim,45,navdrf)
        call gmem_free(i10)
c
c  -----  wt
c
        i10 = igmem_alloc(nwtr*nwtc)
        call daread(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,31)
        call dawrit(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,30,navdrf)
        call gmem_free(i10)
c
c  -----  vr
c
        i10 = igmem_alloc(nwtr*nwtc)
        call daread(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,33)
        call dawrit(idafdrf,iodadrf,xscm(i10),nwtr*nwtc,32,navdrf)
        call gmem_free(i10)
c
      endif
c
      if (neqrf .eq. 1) then
c
c  -----  vrneq
c
        i10 = igmem_alloc(ndim*ndim)
        call daread(idafdrf,iodadrf,xscm(i10),ndim*ndim,35)
        call dawrit(idafdrf,iodadrf,xscm(i10),ndim*ndim,34,navdrf)
        call gmem_free(i10)
c
        if (neq2eps .eq. 1) then
          i10 = igmem_alloc(ndim*ndim)
          call daread(idafdrf,iodadrf,xscm(i10),ndim*ndim,37)
          call dawrit(idafdrf,iodadrf,xscm(i10),ndim*ndim,36,navdrf)
          call gmem_free(i10)
        endif
c
      endif
c
c     if (isolsav .eq. 1) then
c 1-----
c  -----  static potential
c
c       call daread(idafsta,iodasta,xscm,nqdim,index)
c       call dawrit(idafsta,iodasta,xscm,nqdim,index+1,navsta)
c
c  -----  dispersion interactions
c
c       call daread(idafdis,iodadis,xscm,nqcls,index)
c       call dawrit(idafdis,iodadis,xscm,nqcls,index+1,navdis)
c
c  -----  repulsion interactions
c
c       call daread(idafrep,iodarep,xscm,nqcls,index)
c       call dawrit(idafrep,iodarep,xscm,nqcls,index+1,navrep)
c
c  -----  qm repulsion interaction
c
c       call daread(idafrqm,iodarqm,xscm,1,index)
c       call dawrit(idafrqm,iodarqm,xscm,1,index+1,navcls)
c
c  -----  classical energy costs
c
c       call daread(idafcst,iodacst,xscm,3,index)
c       call dawrit(idafcst,iodacst,xscm,3,index+1,navcst)
c
c       if (field(5:) .ne. ' ') then
c   2-----
c    -----  vrpot
c
c         call daread(idafind,iodaind,xscm,nqdim,index)
c         call dawrit(idafind,iodaind,xscm,nqdim,index+1,navind)
c
c         if (itwoeps .eq. 1) then
c     3-----
c           call daread(idafino,iodaino,xscm,nqdim,index)
c           call dawrit(idafino,iodaino,xscm,nqdim,index+1,navino)
c     3-----
c         endif
c   2-----
c       endif
c
c       call daread(idafpol,iodapol,xscm,7,index)
c       call dawrit(idafpol,iodapol,xscm,7,index+1,navpol)
c 1-----
c     endif
c
      return
      end
      subroutine checdis(n,x,tooclos,xscm)
c------
c      checks whether the n classical particles with coordinates x
c      are not too close to the boundary surface and on the right
c      side of the boundary
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      dimension x(3,n)
      logical tooclos
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
      dimension xscm(*)
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      dimension sc(3)
c
      data zero /0.0d00/
c
c-----  begin
c
      ilxs = 81
      nbem = nbem1
      if (itwosur .eq. 1) then
        ilxs = 84
        nbem = nbem2
      endif
c
      ixs = 1
c
c-----  read boundary points
c
      ixs = igmem_alloc(3*nbem)
      call daread(idafdrf,iodadrf,xscm(ixs),3*nbem,ilxs)
c
c-----  loop over classical particles to be checked
c
      do 100, i = 1, n
c 1-----
c  -----  loop over boundary elements
c
        do 200, j = 1, nbem
c   2-----
c    -----  calculate distance vector from classical particle
c           to boundary element
c
          call distab(x(1,i),xscm(ixs+(j-1)*3),sc,dist)
c
          if (dist .lt. excld) then
c     3-----
            if (imcout .ge. 5) then
c       4-----
              write(iwr,1001)
 1001         format(' classical atom too close to boundary, ',
     1        'move rejected')
c       4-----
            endif
c
            tooclos = .true.
            return
c     3-----
          endif
c
c    -----  the particle might be at the wrong side of the boundary
c           to find this out calculate dot product between distance
c           vector and normal vector: this must be positive
c
c         call dotprod(sc,xscm(ixs+(j-1)*3),dot,3)
          dot = ddot(3,sc,1,xscm(ixs+(j-1)*3),1)
c
          if (dot .lt. zero) then
c      3-----
            if (imcout .ge. 5) then
c       4-----
              write(iwr,1002)
 1002         format(' classical atom wrong side of boundary, ',
     1        'move rejected')
c       4-----
            endif
c
            tooclos = .true.
            return
c     3-----
          endif
c   2-----  next boundary element
  200   continue
c 1-----  next classical particle
  100 continue
c
      call gmem_free(ixs)
c
      return
      end
      subroutine getpesc(ipesout,ipes,lpes,lpesc,iopt,
     2                   ieffpol,isodis,ithole,afact,istrt,npt,lcor,
     3                   xn,poln,coord,kpol)
c------
c      gets new coordinates for classical group moved in a pes scan
c      and calculates (if necessary) a new polarisability tensor
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      dimension xn(3,npt)
      dimension poln(6), coord(lcor)
      dimension kpol(npt)
c
c      dimension polatt(6,npt)
      common/scrtch/clas(3000),polat(3000),polatt(6,3000)
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      dimension centpol(3)
c
      data zero /0.0d00/
c
c-----  begin
c
c-----  read all coordinates
c
      if (npt.gt.3000) call caserr
     + ('getpesc: npt>3000 increase dimension of polat, polatt >3000')
      call daread (idafdrf,iodadrf,coord,lpes,98)
c
c-----  calculate starting address for this particular configuration
c
      nstrt = (ipes-1)*lpesc
c
c-----  loop over all group members (including centre)
c
      do 1000, i = 1, npt
c 1-----
c  -----  read new coordinates into xn
c  -----  no group centres recorded for options 0 and 2
c         they are to be computed!!
c
        if ((iopt .ne. 1) .and. (i .eq. 1)) goto 1000
c
        do 100, k = 1, 3
c   2-----
          nstrt = nstrt + 1
          xn(k,i) = coord(nstrt)
c   2-----
  100   continue
c 1-----
 1000 continue
c
c-----  recalculate the polarisability tensor and
c       group centre from the given new coordinates
c
      if (iopt .ne. 1) then
c 1-----
        call clear (centpol,3)
        totpol = zero
c
        do 2000, i = 2, npt
c   2-----
c    -----  calculate centre of polarisability
c
          do 1100, k = 1, 3
            centpol(k) = centpol(k) + xn(k,i)*alfext(istrt+i-1)
 1100     continue
c
          totpol = totpol + alfext(istrt+i-1)
c
c    -----  initialise kpol, array pointing to correct atomic
c           polarisability
c
          kpol(i-1) = i
c   2-----
 2000   continue
c
        do 2200, k = 1, 3
          centpol(k) = centpol(k)/totpol
          xn(k,1) = centpol(k)
 2200   continue
c
        np3 = 3*(npt-1)
        nnp3 = (np3*(np3+1))/2
c
        call clear (coord,nnp3)
c
        call drfamat
     1   (coord,1,.false.,npt-1,xn,kpol,alfext(istrt+1),
     2    ipesout,afact,ithole)
c
c  -----  invert polarisability coupling matrix
c
        call linv3p(coord,auxx,1,np3,ier)
c
        if (ipesout .ge. 3)
     1    call hatout(coord,np3,np3,3,'group-amat')
c
        if ((ieffpol .eq. 1) .and. (isodis .eq. 0)) then
c   2-----
c    -----  make the group polarisability matrix -amat- block diagonal
c           to get effective atom polarisabilities for use in
c           calculation of dispersion energy
c
          np = npt - 1
          call effpol(np,np3,ia,coord,coord(nnp3+1),
     1            polat,polatt)
c
          do 3100, i = 2, npt
c     3-----
            atpol(istrt+i-1) = polat(i-1)
            do 3200, k = 1, 6
              atpolt(istrt+i-1,k) = polatt(i-1,k)
 3200       continue
c     3-----
 3100     continue
c   2-----
        endif
c
c  -----  contract polarizability matrix to 3x3 point polarizability
c
        call drfreda(coord,poln,npt-1)
        if (ipesout .ge. 2) call hatout(poln,3,3,3,'grouppol')
c 1-----
      endif
c
      return
      end
      subroutine mcblout(nstat,isamp,iblock,nmovt,iacc,uav,uav2,
     1                   uavq,uavq2,uavc,uavc2,uavce,uavce2,
     2                   uavcd,uavcd2,uavcr,uavcr2,uavi,uavi2,
     3                   uave,uave2,uavd,uavd2,uavr,uavr2,uavp,uavp2,
     4                   uavn,uavn2,uavnp,uavnp2,uavdf,uavdf2,
     5                   uavcg,uavcg2,uavceg,uavceg2,
     6                   uavcdg,uavcdg2,uavcrg,uavcrg2,uavig,uavig2,
     7                   uaveg,uaveg2,uavdg,uavdg2,uavrg,uavrg2,string)
c------
c       output of means etc during mc run
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      dimension uav(nstat), uav2(nstat), uavq(nstat), uavq2(nstat),
     1  uavi(nstat), uavi2(nstat), uave(nstat), uave2(nstat),
     2  uavd(nstat), uavd2(nstat), uavr(nstat), uavr2(nstat),
     3  uavp(nstat), uavp2(nstat), uavn(nstat), uavn2(nstat),
     4  uavnp(nstat), uavnp2(nstat)
c
      dimension uavdf(nstat*(nstat+1)/2),uavdf2(nstat*(nstat+1)/2)
c
      dimension uavcg(mxgrpar), uavcg2(mxgrpar), uavceg(mxgrpar),
     1 uavceg2(mxgrpar), uavcdg(mxgrpar), uavcdg2(mxgrpar)
c
      dimension uavig(mxgran,mxst),uavig2(mxgran,mxst),
     1 uaveg(mxgran,mxst),uaveg2(mxgran,mxst),uavdg(mxgran,mxst),
     2 uavdg2(mxgran,mxst),uavrg(mxgran,mxst),uavrg2(mxgran,mxst)
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      character*40 string
      character*40 string2(4)
c
      data zero /0.0d00/
c
      data string2(1), string2(2), string2(3), string2(4)
     1 /' for second state ', ' for third state ',
     2  ' for fourth state ', ' for fifth state '/
c
c-----  begin
c
      rmsc2 = uavc2 - uavc**2
      if (rmsc2 .le. zero) rmsc2 = zero
      rmsc = sqrt(rmsc2)
c
      rmsce2 = uavce2 - uavce**2
      if (rmsce2 .le. zero) rmsce2 = zero
      rmsce = sqrt(rmsce2)
c
      rmscd2 = uavcd2 - uavcd**2
      if (rmscd2 .le. zero) rmscd2 = zero
      rmscd = sqrt(rmscd2)
c
      rmscr2 = uavcr2 - uavcr**2
      if (rmscr2 .le. zero) rmscr2 = zero
      rmscr = sqrt(rmscr2)
c
      rmsi2 = uavi2(1) - uavi(1)**2
      if (rmsi2 .le. zero) rmsi2 = zero
      rmsi = sqrt(rmsi2)
c
      rmse2 = uave2(1) - uave(1)**2
      if (rmse2 .le. zero) rmse2 = zero
      rmse = sqrt(rmse2)
c
      rmsd2 = uavd2(1) - uavd(1)**2
      if (rmsd2 .le. zero) rmsd2 = zero
      rmsd = sqrt(rmsd2)
c
      rmsr2 = uavr2(1) - uavr(1)**2
      if (rmsr2 .le. zero) rmsr2 = zero
      rmsr = sqrt(rmsr2)
c
      rmsp2 = uavp2(1) - uavp(1)**2
      if (rmsp2 .le. zero) rmsp2 = zero
      rmsp = sqrt(rmsp2)
c
      rmsn2 = uavn2(1) - uavn(1)**2
      if (rmsn2 .le. zero) rmsn2 = zero
      rmsn = sqrt(rmsn2)
c
      rmsnp2 = uavnp2(1) - uavnp(1)**2
      if (rmsnp2 .le. zero) rmsnp2 = zero
      rmsnp = sqrt(rmsnp2)
c
      write(iwr,1) string
    1 format(a40)
c
      rms2 = uav2(1) - uav(1)**2
      if (rms2 .le. zero) rms2 = zero
      rms = sqrt(rms2)
c
      write(iwr,2) isamp,iblock,nmovt,iacc,uav(1),rms
    2 format(2(3x,i4),2i10,2e20.10)
c
      rmsq2 = uavq2(1) - uavq(1)**2
      if (rmsq2 .le. zero) rmsq2 = zero
      rmsq = sqrt(rmsq2)
c
      write(iwr,3) uavq(1),rmsq
    3 format('  energy of qm system             =',2e20.10)
c
      write(iwr,4) uavc,rmsc
    4 format('  total classical energy          =',2e20.10)
c
      if (iclinte .eq. 1) then
        write(iwr,5) uavce,rmsce
    5   format('  electrostatic classical energy  =',2e20.10)
      endif
c
      if (iclintd .eq. 1) then
        write(iwr,6) uavcd,rmscd
    6   format('  classical model dispersion ene  =',2e20.10)
      endif
c
      if (iclintr .eq. 1) then
        write(iwr,7) uavcr,rmscr
    7   format('  classical model repulsion ene   =',2e20.10)
      endif
c
      write(iwr,8) uavi(1), rmsi
    8 format('  total interaction energy        =',2e20.10)
c
      write(iwr,9) uave(1),rmse
    9 format('  electrostatic interaction energy=',2e20.10)
c
      if (gamdrf .ne. zero) then
        write(iwr,10) uavd(1),rmsd
   10   format('  approximate dispersion energy   =',2e20.10)
      endif
c
      write(iwr,11) uavr(1),rmsr
   11 format('  model repulsion energy          =',2e20.10)
c
      write(iwr,12) uavp(1),rmsp
   12 format('  equilibrium polarisation energy =',2e20.10)
c
      if (neqrf .eq. 1) then
c 1-----
        write(iwr,13) uavn(1),rmsn
   13   format('  non-equilibrium rf energy       =',2e20.10)
        write(iwr,14) uavnp(1),rmsnp
   14   format('  non-equilibrium pol energy cost =',2e20.10)
c 1-----
      endif
c
      do 100, igr = 1, ngran
c 1-----
c 1-----
  100 continue
c
      do 200, ist = 2, imcst
c 1-----
        rms2 = uav2(ist) - uav(ist)**2
        if (rms2 .le. zero) rms2 = zero
        rms = sqrt(rms2)
c
        rmsq2 = uavq2(ist) - uavq(ist)**2
        if (rmsq2 .le. zero) rmsq2 = zero
        rmsq = sqrt(rmsq2)
c
        rmsi2 = uavi2(ist) - uavi(ist)**2
        if (rmsi2 .le. zero) rmsi2 = zero
        rmsi = sqrt(rmsi2)
c
        rmse2 = uave2(ist) - uave(ist)**2
        if (rmse2 .le. zero) rmse2 = zero
        rmse = sqrt(rmse2)
c
        rmsd2 = uavd2(ist) - uavd(ist)**2
        if (rmsd2 .le. zero) rmsd2 = zero
        rmsd = sqrt(rmsd2)
c
        rmsr2 = uavr2(ist) - uavr(ist)**2
        if (rmsr2 .le. zero) rmsr2 = zero
        rmsr = sqrt(rmsr2)
c
        rmsp2 = uavp2(ist) - uavp(ist)**2
        if (rmsp2 .le. zero) rmsp2 = zero
        rmsp = sqrt(rmsp2)
c
        rmsn2 = uavn2(ist) - uavn(ist)**2
        if (rmsn2 .le. zero) rmsn2 = zero
        rmsn = sqrt(rmsn2)
c
        rmsnp2 = uavnp2(ist) - uavnp(ist)**2
        if (rmsnp2 .le. zero) rmsnp2 = zero
        rmsnp = sqrt(rmsnp2)
c
        write(iwr,15) string2(ist-1)
   15   format(a40)
c
        write(iwr,2) isamp,iblock,nmovt,iacc,uav(ist),rms
        write(iwr,3) uavq(ist),rmsq
c
        write(iwr,8) uavi(ist),rmsi
c
        write(iwr,9) uave(ist),rmse
c
        if (gammc(ist) .ne. zero) then
          write(iwr,10) uavd(ist),rmsd
        endif
c
        write(iwr,11) uavr(ist),rmsr
c
        write(iwr,12) uavp(ist),rmsp
c
        if (neqrf .eq. 1) then
          write(iwr,13) uavn(ist),rmsn
          write(iwr,14) uavnp(ist),rmsnp
        endif
c
        do 120, kst = 1, ist
c   2-----
          indxdf = ia(ist) + kst
          rmsdf2 = uavdf2(indxdf) - uavdf(indxdf)**2
          if (rmsdf2 .le. zero) rmsdf2 = zero
          rmsdf = sqrt(rmsdf2)
c
          write(iwr,16) uavdf(indxdf),rmsdf
   16   format('  state difference energy         =',2e20.10)
c   2-----
  120   continue
c 1-----
  200 continue
c
      return
      end
      subroutine mcout
c------
c      prints monte carlo energy results
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
      real*8 ureft, uold, eqm, eqmol, eclas, eclasol, eclase, ecleol
      real*8 eclasd, ecldol, eclasr, eclrol, eint, eintol, eelst, eelstol
      real*8 edisp, edisol, erep, erepol, extnol, epol, epolol
      real*8 eneq, eneqol, eneqp, eneqpol, edifol
      common /mcener/ ureft(mxst),uold(mxst),eqm(mxst),eqmol(mxst),
     +  eclas,eclasol,eclase,ecleol,eclasd,ecldol,eclasr,eclrol,
     +  eint(mxst),eintol(mxst),eelst(mxst),eelstol(mxst),edisp(mxst),
     +  edisol(mxst),erep(mxst),erepol(mxst),extnol(mxst),epol(mxst),
     +  epolol(mxst),eneq(mxst),eneqol(mxst),eneqp(mxst),eneqpol(mxst),
     +  edifol(mxst*(mxst+1)/2)
c
      real*8 utav, utav2, utavq, utavq2, utavc
      real*8 utavc2, utavce, utavce2, utavcd, utavcd2, utavcr, utavcr2
      real*8 utavi, utavi2, utave, utave2, utavd, utavd2, utavr, utavr2
      real*8 utavp, utavp2, utavn, utavn2, utavnp, utavnp2
      real*8 utavdf, utavdf2, q1tot
      common /mcener2/ utav(mxst),utav2(mxst),utavq(mxst),utavq2(mxst),
     +  utavc,utavc2,utavce,utavce2,utavcd,utavcd2,utavcr,utavcr2,
     +  utavi(mxst),utavi2(mxst),utave(mxst),utave2(mxst),utavd(mxst),
     +  utavd2(mxst),utavr(mxst),utavr2(mxst),utavp(mxst),utavp2(mxst),
     +  utavn(mxst),utavn2(mxst),utavnp(mxst),utavnp2(mxst),
     +  utavdf(mxst*(mxst+1)/2),utavdf2(mxst*(mxst+1)/2),q1tot(mxst)
c
      integer itotacc
      common /mcres/ itotacc
c
      real*8 eclasg, eclsolg, eclaseg, ecleolg, eclasdg, ecldolg
      real*8 eclasrg, eclrolg, eintg, eintolg, eelstg, elstolg
      real*8 edispg, edisolg, erepg, erepolg, extnolg
      common/mcener3/eclasg(mxgrpar),eclsolg(mxgrpar),eclaseg(mxgrpar),
     +  ecleolg(mxgrpar),eclasdg(mxgrpar),ecldolg(mxgrpar),
     +  eclasrg(mxgrpar),eclrolg(mxgrpar),eintg(mxgran,mxst),
     +  eintolg(mxgran,mxst),eelstg(mxgran,mxst),elstolg(mxgran,mxst),
     +  edispg(mxgran,mxst),edisolg(mxgran,mxst),erepg(mxgran,mxst),
     +  erepolg(mxgran,mxst),extnolg(mxgran,mxst)
c
      real*8 utavcg, utavcg2, utavceg, utvceg2, utavcdg, utvcdg2
      real*8 utavcrg, utvcrg2, utavig, utavig2, utaveg, utaveg2
      real*8 utavdg, utavdg2, utavrg, utavrg2
      common/mcener4/utavcg(mxgrpar),utavcg2(mxgrpar),utavceg(mxgrpar),
     +  utvceg2(mxgrpar),utavcdg(mxgrpar),utvcdg2(mxgrpar),
     +  utavcrg(mxgrpar),utvcrg2(mxgrpar),utavig(mxgran,mxst),
     +  utavig2(mxgran,mxst),utaveg(mxgran,mxst),utaveg2(mxgran,mxst),
     +  utavdg(mxgran,mxst),utavdg2(mxgran,mxst),utavrg(mxgran,mxst),
     +  utavrg2(mxgran,mxst)
c
c
      dimension a1(mxst)
c
      data zero,one  /0.0d00, 1.0d00/
c
c-----  begin
c
      ist = 1
c
      rms2 = utav2(ist) - utav(ist)**2
      if (rms2 .le. zero) rms2 = zero
      rms = sqrt(rms2)
c
      rmsq2 = utavq2(ist) - utavq(ist)**2
      if (rmsq2 .le. zero) rmsq2 = zero
      rmsq = sqrt(rmsq2)
c
      rmsc2 = utavc2 - utavc**2
      if (rmsc2 .le. zero) rmsc2 = zero
      rmsc = sqrt(rmsc2)
c
      rmsce2 = utavce2 - utavce**2
      if (rmsce2 .le. zero) rmsce2 = zero
      rmsce = sqrt(rmsce2)
c
      rmscd2 = utavcd2 - utavcd**2
      if (rmscd2 .le. zero) rmscd2 = zero
      rmscd = sqrt(rmscd2)
c
      rmscr2 = utavcr2 - utavcr**2
      if (rmscr2 .le. zero) rmscr2 = zero
      rmscr = sqrt(rmscr2)
c
      rmsi2 = utavi2(ist) - utavi(ist)**2
      if (rmsi2 .le. zero) rmsi2 = zero
      rmsi = sqrt(rmsi2)
c
      rmse2 = utave2(ist) - utave(ist)**2
      if (rmse2 .le. zero) rmse2 = zero
      rmse = sqrt(rmse2)
c
      rmsd2 = utavd2(ist) - utavd(ist)**2
      if (rmsd2 .le. zero) rmsd2 = zero
      rmsd = sqrt(rmsd2)
c
      rmsr2 = utavr2(ist) - utavr(ist)**2
      if (rmsr2 .le. zero) rmsr2 = zero
      rmsr = sqrt(rmsr2)
c
      rmsp2 = utavp2(ist) - utavp(ist)**2
      if (rmsp2 .le. zero) rmsp2 = zero
      rmsp = sqrt(rmsp2)
c
      rmsn2 = utavn2(ist) - utavn(ist)**2
      if (rmsn2 .le. zero) rmsn2 = zero
      rmsn = sqrt(rmsn2)
c
      rmsnp2 = utavnp2(ist) - utavnp(ist)**2
      if (rmsnp2 .le. zero) rmsnp2 = zero
      rmsnp = sqrt(rmsnp2)
c
      a1(ist) = (one/onekt)*log(q1tot(ist))
c
      write(iwr,1001) ngrpol,ngrpmc, nsamp, nblock, nmoves
 1001 format(//,
     1' =========================================================== ',/,
     2' =                monte carlo results:                     = ',/,
     3' = number of groups                   : ',i10,'         = ',/,
     4' = number of groups in monte carlo run: ',i10,'         = ',/,
     5' = number of samples                  : ',i10,'         = ',/,
     6' = number of blocks in each sample    : ',i10,'         = ',/,
     7' = number of moves in each block      : ',i10,'         = ',/,
     8' =========================================================== ')
c
      nmovt = nsamp*nblock*nmoves
      totacc = dble(itotacc)
      totmov = dble(nmovt)
      peracc = totacc/totmov
c
      write(iwr,1002) nmovt, itotacc, peracc*100.d00
 1002 format(/,
     1' total number of attempted monte carlo moves: ',i10,/,
     2' total number of accepted moves             : ',i10,/,
     3' acceptance percentage                      : ',f10.5)
c
      if (ncheck .ne. 0) then
c 1-----
        pi = 4.0*atan(1.0)
        write (iwr,101) darot/pi, dtrans
  101   format (/, ' final value of darot = ',f10.5,
     2   ' ; dtrans = ',f10.5)
c 1-----
      endif
c
      write(iwr,1003) utav(ist),rms
 1003 format(/,
     1' average total energy (a.u.)                : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      write(iwr,1004) utavq(ist),rmsq
 1004 format(/,
     1' average quantum system energy (a.u.)       : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      write(iwr,1005) utavc,rmsc
 1005 format(/,
     1' average total classical energy (a.u.)      : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      if (iclinte .eq. 1) then
        write(iwr,1006) utavce,rmsce
 1006   format(/,
     1' average electrostatic classical ener (a.u.): ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
      endif
c
      if (iclintd .eq. 1) then
        write(iwr,1007) utavcd,rmscd
 1007   format(/,
     1' average model dispersion clas ener (a.u.)  : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
      endif
c
      if (iclintr .eq. 1) then
        write(iwr,1008) utavcr,rmscr
 1008   format(/,
     1' average model repulsion clas ener (a.u.)   : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
      endif
c
      do 100, igr = 1, ngran
c 1-----
        do 20, jgr = 1, igr
c   2-----
          indxcl = ia(igr) + jgr
c
          write(iwr,1009) jgr, igr
 1009     format(/,
     1' interaction between classical analysis groups ',i2,' and ', i2)
c
          rmsc2 = utavcg2(indxcl) - utavcg(indxcl)**2
          if (rmsc2 .le. zero) rmsc2 = zero
          rmsc = sqrt(rmsc2)
          write(iwr,1005) utavcg(indxcl), rmsc
c
          if (iclinte .eq. 1) then
            rmsce2 = utvceg2(indxcl) - utavceg(indxcl)**2
            if (rmsce2 .le. zero) rmsce2 = zero
            rmsce = sqrt(rmsce2)
            write(iwr,1006) utavceg(indxcl),rmsce
          endif
c
          if (iclintd .eq. 1) then
            rmscd2 = utvcdg2(indxcl) - utavcdg(indxcl)**2
            if (rmscd2 .le. zero) rmscd2 = zero
            rmscd = sqrt(rmscd2)
            write(iwr,1007) utavcdg(indxcl), rmscd
          endif
c
          if (iclintr .eq. 1) then
            rmscr2 = utvcrg2(indxcl) - utavcrg(indxcl)**2
            if (rmscr2 .le. zero) rmscr2 = zero
            rmscr = sqrt(rmscr2)
            write(iwr,1008) utavcrg(indxcl), rmscr
          endif
c   2-----  next jgr
   20   continue
c 1-----  next igr
  100 continue
c
      write (iwr,1100)
 1100 format (/,
     1' total qm - classical interaction energies:')
c
      write(iwr,1101) utavi(ist),rmsi
 1101 format(/,
     1' average interaction energy (a.u.)          : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      write(iwr,1102) utave(ist),rmse
 1102 format(/,
     1' average electrost interaction ene (a.u.)   : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      if (gamdrf .ne. zero) then
        write(iwr,1103) utavd(ist),rmsd
 1103   format(/,
     1' average approximate dispersion  (a.u.)     : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
      endif
c
      write(iwr,1104) utavr(ist),rmsr
 1104 format(/,
     1' average model repulsion energy (a.u.)      : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      write(iwr,1105) utavp(ist),rmsp
 1105 format(/,
     1' average polarisation energy (a.u.)         : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
      if (neqrf .eq. 1) then
        write(iwr,1106) utavn(ist),rmsn
 1106   format(/,
     1' average non-equilibrium rf energy (a.u.)   : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c
        write(iwr,1107) utavnp(ist),rmsnp
 1107   format(/,
     1' average non-equilibrium energy cost (a.u.) : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
      endif
c
      if (gamdrf .ne. zero) write (iwr,1200)
 1200 format (/,
     1' warning: no separate analysis possible for dispersion',/,
     2' dispersion not included in interaction with analysis groups')
c
      do 200, igr = 1, ngran
c 1-----
        write(iwr,1201) igr
 1201   format(/,
     1' interaction with classical analysis group: ',i2)
c
        rmsi2 = utavig2(igr,ist) - utavig(igr,ist)**2
        if (rmsi2 .le. zero) rmsi2 = zero
        rmsi = sqrt(rmsi2)
        write(iwr,1101) utavig(igr,ist), rmsi
c
        rmse2 = utaveg2(igr,ist) - utaveg(igr,ist)**2
        if (rmse2 .le. zero) rmse2 = zero
        rmse = sqrt(rmse2)
        write(iwr,1102) utaveg(igr,ist), rmse
c
        rmsr2 = utavrg2(igr,ist) - utavrg(igr,ist)**2
        if (rmsr2 .le. zero) rmsr2 = zero
        rmsr = sqrt(rmsr2)
        write(iwr,1104) utavrg(igr,ist), rmsr
c 1-----
  200 continue
c
      write(iwr,1301) ureft(ist)
 1301 format(/,
     1' reference total energy                     : ',e20.10)
c
      write(iwr,1302) q1tot(ist), a1(ist)
 1302 format(/,
     1' partition function                         : ',e20.10,/,
     2' (helmholtz) free energy                    : ',e20.10)
c
      write(iwr,1303) (a1(ist) + ureft(ist))
 1303 format(/,
     1' total (helmholtz) free energy              : ',e20.10)
c
      do 300, ist = 2, imcst
c 1-----
        write(iwr,1401) ist
 1401   format(/,
     1' analysis for state number: ',i2)
c
        rms2 = utav2(ist) - utav(ist)**2
        if (rms2 .le. zero) rms2 = zero
        rms = sqrt(rms2)
c
        rmsq2 = utavq2(ist) - utavq(ist)**2
        if (rmsq2 .le. zero) rmsq2 = zero
        rmsq = sqrt(rmsq2)
c
        rmsi2 = utavi2(ist) - utavi(ist)**2
        if (rmsi2 .le. zero) rmsi2 = zero
        rmsi = sqrt(rmsi2)
c
        rmse2 = utave2(ist) - utave(ist)**2
        if (rmse2 .le. zero) rmse2 = zero
        rmse = sqrt(rmse2)
c
        rmsd2 = utavd2(ist) - utavd(ist)**2
        if (rmsd2 .le. zero) rmsd2 = zero
        rmsd = sqrt(rmsd2)
c
        rmsr2 = utavr2(ist) - utavr(ist)**2
        if (rmsr2 .le. zero) rmsr2 = zero
        rmsr = sqrt(rmsr2)
c
        rmsp2 = utavp2(ist) - utavp(ist)**2
        if (rmsp2 .le. zero) rmsp2 = zero
        rmsp = sqrt(rmsp2)
c
        rmsn2 = utavn2(ist) - utavn(ist)**2
        if (rmsn2 .le. zero) rmsn2 = zero
        rmsn = sqrt(rmsn2)
c
        rmsnp2 = utavnp2(ist) - utavnp(ist)**2
        if (rmsnp2 .le. zero) rmsnp2 = zero
        rmsnp = sqrt(rmsnp2)
c
        a1(ist) = (one/onekt)*log(q1tot(ist))
c
        write(iwr,1003) utav(ist),rms
c
        write(iwr,1004) utavq(ist),rmsq
c
        write(iwr,1100)
c
        write(iwr,1101) utavi(ist),rmsi
c
        write(iwr,1102) utave(ist),rmse
c
        if (gamdrf .ne. zero) then
          write(iwr,1103) utavd(ist),rmsd
        endif
c
        write(iwr,1104) utavr(ist),rmsr
c
        write(iwr,1105) utavp(ist),rmsp
c
        if (neqrf .eq. 1) then
          write(iwr,1106) utavn(ist),rmsn
          write(iwr,1107) utavnp(ist),rmsnp
        endif
c
        do 220, kst = 1, ist
c   2-----
          indxdf = ia(ist) + kst
          rmsdf2 = utavdf2(indxdf) - utavdf(indxdf)**2
          if (rmsdf2 .le. zero) rmsdf2 = zero
          rmsdf = sqrt(rmsdf2)
          write(iwr,1501) utavdf(indxdf), rmsdf
 1501     format(/,
     1' average state difference energy (a.u.)     : ',e20.10,/,
     2' r.m.s. deviation                           : ',e20.10)
c   2-----
  220   continue
c
        if (gamdrf .ne. zero) write (iwr,1200)
c
        do 240, igr = 1, ngran
c   2-----
          write(iwr,1201) igr
c
          rmsi2 = utavig2(igr,ist) - utavig(igr,ist)**2
          if (rmsi2 .le. zero) rmsi2 = zero
          rmsi = sqrt(rmsi2)
          write(iwr,1101) utavig(igr,ist), rmsi
c
          rmse2 = utaveg2(igr,ist) - utaveg(igr,ist)**2
          if (rmse2 .le. zero) rmse2 = zero
          rmse = sqrt(rmse2)
          write(iwr,1102) utaveg(igr,ist), rmse
c
          rmsr2 = utavrg2(igr,ist) - utavrg(igr,ist)**2
          if (rmsr2 .le. zero) rmsr2 = zero
          rmsr = sqrt(rmsr2)
          write(iwr,1104) utavrg(igr,ist), rmsr
c   2-----
  240   continue
c
        write(iwr,1301) ureft(ist)
c
        write(iwr,1302) q1tot(ist), a1(ist)
c
        write(iwr,1303) (a1(ist) + ureft(ist))
c 1-----
  300 continue
c
      return
      end
      subroutine mattrns(a,b,c,n,m,l)
c---------------
c
c       this subroutine performs a general matrix multiplication
c       of an nxm matrix a with a mxl matrix b to give an nxl matrix c
c
c                        c = a.b
c
c---------------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      integer i,j,k,l,m,n
c
      dimension a(n,m), b(m,l), c(n,l)
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
c
      do 100, i = 1, n
        do 100, j = 1, l
          sum = zero
          do 200, k = 1, m
            sum = sum + a(i,k)*b(k,j)
  200     continue
        c(i,j) = sum
  100 continue
      return
      end
      subroutine ranset(iseed)
c***********************************************************************
c     ranset sets the seed in /rndgen/ to be used by ranget
c***********************************************************************
      integer seed
      common /rndgen/ seed
      seed=(iseed/2)*2+1
      return
      end
      function ranget()
c***********************************************************************
c     ranget calls random number routine with seed stored in /rndgen/
c     here the routine random is used.
c***********************************************************************
      integer seed
      real*8 x, ranget
      common /rndgen/ seed
      call random(x,seed)
      ranget=x
      return
      end
      subroutine random (rand,ig)
c
cccccc r. geurtsen, groningen, wfvg, july 1987 ccccccccccccccccccccccccc
c      modified by: ton rullmann, utrecht, february 1990               c
c                                                                      c
c     subroutine random (rand,ig)                                      c
c                                                                      c
comment   random generates a random number rand, using a linear        c
c     congruential method. the recursion formula                       c
c                                                                      c
c         irand = mod(irand * b + 1, a)                                c
c                                                                      c
c     is used with  b = 31415821  and  a = 100000000. the last digit   c
c     from the random integer irand is chopped of, and the number      c
c     is scaled to a real value rand between 0 and 1, including 0 but  c
c     excluding 1.                                                     c
c                                                                      c
c     rand = delivered with random number between 0 and 1              c
c     ig = random number generator seed, is delivered with random      c
c          integer                                                     c
c                                                                      c
c modified: on subsequent entries irand is again calculated from ig,   c
c     just as on the first entry. for a continued series the user must c
c     provide the output value of irand as input to the next call.     c
c modified: print error message if out of range (include file needed)  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
      parameter (m=100000000, m1=10000, mult=31415821)
      data       irand /0/, new /0/
c
c
c*****tr: next line commented, so irand is recalculated from ig
c      if (new .ne. 0) goto 7
      new = 1
      irand = mod (iabs(ig),m)
    7 continue
c
c*****multiply irand by mult, but take into account that overflow must
c*****be discarded, and do not generate an error.
c
      irandh = irand / m1
      irandl = mod(irand, m1)
      multh = mult / m1
      multl = mod(mult, m1)
c
      irand = mod(irandh*multl + irandl*multh, m1) * m1 + irandl*multl
      irand = mod(irand + 1, m)
c
c*****convert irand to a real random number between 0 and 1.
c
      r = dble(irand / 10) * 10 / dble(m)
      if ((r .le. 0.d0) .or. (r .gt. 1.d0)) then
         r = 0.d0
         write (6,*)
     +     ' *** random: random number out of range, set to 0.'
      endif
      rand = r
      ig = irand
c
      return
      end
      subroutine drfaupd(amat,ndim,np,
     1                   xpts,mpol,polar,afact,ithole)
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c                                                                      c
c                                                                      c
c     this routine updates that part of the relay matrix to do with    c
c     the coupling between polarizabilities                            c
c                                                                      c
c     original literature (formulated for polarizabilities only!)      c
c     thole & van duijnen: theor.chim.acta (1980) 55 307 (form. 8,9)   c
c                                                                      c
c     thole: chem.phys.(1981) 59 341  (form.5)                         c
c                                                                      c
c                                                                      c
c     routines called:       distab,drftpq,hatout,linv3p               c
c                                                                      c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      dimension amat(ndim,ndim), xpts(3,1)
      dimension mpol(*), polar(*)
c
c-----  for communication with hondo8
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
c-----  local variables
c
      dimension p(3), q(3), pq(3)
      dimension bpq(3,3)
      dimension thspoln(6)
c
      data one,sixth /1.0d0,.16666666666666666666667d00/
c
c-----  begin
c
c-----  set pointer to modified group polarizability
c
      ip = inpt(1)
c
c-----  starting address in -a- for point -i-
c
      if = (ip-1)*3
c
c-----  invert group polarizability
c
c     call copyv (thispol,thspoln,6)
      call dcopy(6,thispol,1,thspoln,1)
      call linv3p(thspoln,auxx,1,3,ier)
      if (imcout .eq. 5) call hatout(thspoln,3,3,3,'grpol-1')
c
c-----  copy -thspoln- into -amat-
c
      kl = 0
      do 40, k = 1, 3
        do 50, l = 1, k
          kl = kl + 1
          amat(if+k,if+l) = thspoln(kl)
          amat(if+l,if+k) = thspoln(kl)
   50   continue
   40 continue
c
c-----  if group is translated calulate new t(p;q)
c
      if (notrans .eq. 0) then
c 1-----
        p(1) = xnpts(1,1)
        p(2) = xnpts(2,1)
        p(3) = xnpts(3,1)
c
        alfai = polar(ip)
c
c  -----  loop over polarizable points (to calculate t(i;j))
c
        do 100, jp = 1, np
c   2-----
c    -----  skip updated polarisability
c
          if (ip .eq. jp) goto 100
c
          jpol = mpol(jp)
          q(1) = xpts(1,jpol)
          q(2) = xpts(2,jpol)
          q(3) = xpts(3,jpol)
c
c    -----  staring address in -a- for point -j-
c
          jf = (jp-1)*3
c
c    -----  zero sub array
c
c    -----  off diagonal blocks are (modified) dipole interactions
c
          alfaj = polar(jp)
          call distab(p,q,pq,dist)
          dmind1 = one/dist
          dmin3 = dmind1**3
          v = dist/((alfai*alfaj)**sixth)
          call drftpq(pq,dmind1,dmin3,ithole,afact,v,bpq)
c
c    -----  copy -bpq- into -a-
c
          do 200, k = 1, 3
            do 300, l = 1, 3
              amat(if+k,jf+l) = bpq(k,l)
              amat(jf+l,if+k) = bpq(l,k)
  300       continue
  200     continue
c   2-----  next polarisability
  100   continue
c 1-----
      endif
c
      return
      end
      subroutine drfdfupd(ieps,b,xsurf,xnorm,area)
c------
c             routine updates coupling matrix elements between
c             classical polarizabilities and boundary surface elements
c             according to a move in a monte carlo simulation
c
c     ndima: dimension of polarizability matrix (if present)
c     ndimb: dimension of boundary problem
c            =   nbem if kappa.eq.0.0
c            = 2*nbem if kappa.ne.0.0
c     ndim: dimension of complete linear problem
c
c     b(ndim,ndim) is the relay-inp matrix, of which the parts
c     to do with the changed polarizability are updated
c
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  dummy arrays
c
      dimension b(ndim,ndim),xsurf(3,nbem),xnorm(3,nbem),area(nbem)
c
c-----  common blocks
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
c-----  local arrays:
c
      dimension tip(3,3),fip(3),xip(3),tipni(3),tipkip(3)
c
      logical kapnoz
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
c
c-----  begin
c
      if (ieps .eq. 0) then
        eps = eps1
        kappa = kappa1
      else
        eps = eps2
        kappa = kappa2
      endif
c
      pi = four*atan(one)
      kappas = kappa**2
      kapnoz = kappa .ne. zero
      epsfact = one/(two*pi*(one+eps))
      expkd = one
      fkapone = one
c
c-----  get index of the modified group polarizability
c
      np = inpt(1)
c
c-----  loop over boundary elements
c
      do 100, ni = 1, nbem
c 1-----
c  -----  calculate distance vector from -ipol- to surface element -ni-
c         -ipol- is the modified group polarizability
c
c  -----  xip = (p-i)
c
        call distab(xsurf(1,ni),xnpts(1,1),xip,dist)
        dmind1 = one/dist
        dmin2 = dmind1**2
        dmin3 = dmin2*dmind1
c
c  -----  fip = f(i;p): field  in -ipol- due to charge in -ni-
c
        do 30, k = 1, 3
          fip(k) = xip(k)*dmin3
   30   continue
c
c  -----  tip = t(i;p): field gradient of charge in -ni- at -ipol-
c
        call drftpq(xip,dmind1,dmin3,0,one,one,tip)
c
c  -----  calculate tipni = t(i;p) .n(i)
c
        call matvec(tip,xnorm(1,ni),tipni,3,.false.)
c
c  -----  calculate factors (only if finite ionic strength)
c
        if (kapnoz) then
          expkd = exp(-kappa*dist)
          fkapone = one + kappa*dist
        endif
c
        do 40, k = 1, 3
c   2-----
c    -----  potential of unit dipole in -ipol- at -ni-
c           = f(p;i) = -f(i;p)
c           scaled with 1/(2pi(1+eps)) to satisfy coupling equations
c
c         the negative of this matrix element should be stored in (i,p)
c
c         plus sign is a result of the use of f(i;p), whereas f(p;i)
c         is needed
c
          b(ndima+ni,(np-1)*3+k) = epsfact*fip(k)
c
c    -----  field of unit dipole density in -ni- at -ipol-
c           = - del(p) k(i;p)*s(i)
c           = - [(eps*(1+kd)*exp(-kd) -1)*t(i;p).n(i)
c
c              - eps*kappa**2*exp(-kd)*f(i;p).n(i)*(p-i)]*s(i)
c
c   *note* the second term is added only if kappa .ne. zero !!
c
c
c         the negative of this matrix element should be stored in (p;i)
c
          b((np-1)*3+k,ndima+ni) = (eps*fkapone*expkd - one)*
     1                              tipni(k)*area(ni)
c   2-----  next component
   40   continue
c
c  -----  only for poisson-boltzmann (i.e. finite ionic strength)
c
        if (kapnoz) then
c   2-----
c    -----  calculate fipni = f(i;p).n(i)
c
c         fipni = adotb(fip,xnorm(1,ni),3)
          fipni = ddot(3,fip,1,xnorm(1,ni),1)
c
c    -----  matrix elements
c
          do 50, k = 1, 3
c     3-----
c      -----  update del(p) k(i;p) *s(i)
c
            b((np-1)*3+k,ndima+ni) = b((np-1)*3+k,ndima+ni) -
     2           eps*kappa2*expkd*
     3           fipni*xip(k)*area(ni)
c
c      -----  minus field of unit dipole in -ipol- at -ni-
c             contracted with normal vector at -ni-
c             = t(p;i) . n(i)
c
c           the negative of this element should be stored in (i+nbem,p)
c           scaled by eps/(2pi(1+eps)) to satisfy the coupling equations
c
            b(ndima+nbem+ni,(np-1)*3+k) = - eps*epsfact*tipni(k)
c
c      -----  field of unit charge density in -ni- at -ipol-
c             = -del(p) l(i;p) * s(i)
c             = - [(1+kd)*exp(-kd) - 1]* f(i;p) * s(i)
c
c           the negative of this element should be stored in (p;i+nbem)
c
            b((np-1)*3+k,ndima+ni+nbem) =
     1               (fkapone*expkd - one)*fip(k)*area(ni)
c     3-----  next component
   50     continue
c   2-----
        endif
c 1-----  next boundary element
  100 continue
c
      return
      end
      subroutine zfpupd(zfp,vrp,eclas,eelst,edisp,erep,egrcls,edum,
     1           ngrpair)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c          routine forms vectors and scalars describing the inter-    c
c          actions between the external points.                       c
c          since the fields and potentials of the point charges       c
c          will be affected by polarization effects, they are stored  c
c          as the last column of the matrix "wt" which will be        c
c          contracted with 'a' (the polarizability)                   c
c                                                                     c
c       clas is the electrostatic energy between the points           c
c                                                                     c
c                                                                     c
c      interactions between points belonging to the same group        c
c      (i.e. having the same group names) are excluded: this is       c
c      essentially a quantum mechanical interaction, and should be    c
c      treated as such                                                c
c                                                                     c
c      interactions between "overlapping" distributions, as           c
c      defined by their polarizabilities, should be reduced           c
c      consistently with the reduction in routine drftpq              c
c                                                                     c
c     routines called: distab,clear                                   c
c                                                                     c
c                                                                     c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  dummy arguments
c
      dimension zfp(ndim,ngran), vrp(ndim,ngran)
      dimension eelst(ngrpair), edisp(ngrpair), erep(ngrpair),
     1          egrcls(ngrpair), edum(5*ngrpair)
c
c-----  common blocks
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
      real*8 thispol
      integer istrt
      common  /mcpol/ thispol(6),istrt
c
c
c-----  local variables
c
      logical updpoli, updpolj
c
c      dimension zfpi(3,ngran)
      dimension zfpi(3,3000)
      dimension p(3), q(3), pq(3)
c
      dimension b(3,3), b2(3,3), polten(3,3)
c
      character*16 namei,namej
      logical poli,polj
c
      data zero,one,two,three,four/0.0d00,1.0d00,2.d0,3.d0,4.d0/
      data pt5,sixth /0.50d00,.1666666666666666666667d00/
      data pt3, pt75, onept5 /
     +     0.333333333333333333333333d0,7.5d-01,1.5d00/
c
      data repcut /1.0d02/
      data small /1.0d-03/
c
c-----  begin
c
      if (ngran.gt.3000) call caserr
     + ('mcsub.f: zfpupd: ngran>3000 increase dimension of zfpi')
      do 50, i = 1, ngrpair
c 1-----
        eelst(i) = zero
        edisp(i) = zero
        erep(i) = zero
        egrcls(i) = zero
c 1-----
   50 continue
c
c-----  count polarizable external points and updated points
c
      imp = 1
      iupd = 1
c
c-----  loop over external points
c
      do 100, ii = 1, nxtpts
c 1-----
        zi = chrg(ii)
        namei = nxcent(ii)
        ipol = mpol(imp)
        poli = ipol.eq.ii
        igrani = igranl(ii)
        alfi = alfext(ii)
        if (alfi .eq. zero) alfi = one
c
        if (poli) then
c   2-----
          call clear(zfpi,3*ngran)
          irow = (imp-1)*3
c   2-----
        endif
c
c  -----  get (updated) coordinate vector
c
        if ((ii .ge. istrt) .and. (ii .lt. istrt+nnpts)) then
c   2-----
          p(1) = xnpts(1,iupd)
          p(2) = xnpts(2,iupd)
          p(3) = xnpts(3,iupd)
          iupd = iupd + 1
          updpoli = .true.
c   2-----
        else
c   2-----
          p(1) = xpts(1,ii)
          p(2) = xpts(2,ii)
          p(3) = xpts(3,ii)
          updpoli = .false.
c   2-----
        endif
c
c  -----  count polarizable points and updated points
c
        jmp = 1
        jupd = 1
c
c  -----  loop over external points
c
        do 90, jj = 1, nxtpts
c   2-----
          jpol = mpol(jmp)
          polj = jpol .eq. jj
          namej = nxcent(jj)
          igranj = igranl(jj)
c
c    -----  get (updated) coordinate vector
c
          if ((jj .ge. istrt) .and. (jj .lt. istrt+nnpts)) then
c     3-----
            q(1) = xnpts(1,jupd)
            q(2) = xnpts(2,jupd)
            q(3) = xnpts(3,jupd)
            jupd = jupd + 1
            updpolj = .true.
c     3-----
          else
c     3-----
            q(1) = xpts(1,jj)
            q(2) = xpts(2,jj)
            q(3) = xpts(3,jj)
            updpolj = .false.
c     3-----
          endif
c
          if (jj .eq. ii) goto 40
c
c    -----  exclude interactions between members of the same group
c
cxxx      if (namei(ngrnam:) .eq. namej(ngrnam:)) goto 40
          if (namei(7:(6+ngrnam)) .eq. namej(7:(6+ngrnam))) goto 40
c
c    -----  set index for analysis of interactions
c
          indxen = ia(max(igrani,igranj)) + min(igrani,igranj)
c
          zj = chrg(jj)
c
c    -----  distance vector between ii and jj in -pq-
c
          call distab(q,p,pq,dist)
          if (dist .lt. small) then
             write(iwr,95) namei, namej, dist
   95        format(/,1x,'WARNING: non-excluded interaction ',
     1   'between ', a16, ' and ', a16, ' at', e15.8,
     2   ' bohr distance: skipped')
             goto 40
          endif
          alfj = alfext(jj)
          if (alfj .eq. zero) alfj = one
          dmind1 = one/dist
c
          v = one
          factp = one
          factf = one
c
c    -----  scale interaction between point charges according to thole
c           (optional)
c
          if (modxza .ne. 0) then
c     3-----
            s = (alfi*alfj)**sixth
            v = dist/s
            if (ithole .eq. 1) then
c       4-----
              if (v .le. afact) then
                av = v/afact
                factp = av**4 - two*av**3 + two*av
                factf = four*av**3 - three*av**4
              endif
            else if (ithole .eq. 2) then
              au = afact*v
              factp = (one - (pt5*au + one)*exp(-au))
              factf = (one - (pt5*au**2 + au + one)*exp(-au))
c       4-----
            endif
c     3-----
          endif
c
          if (poli) then
c     3-----
c      -----  add field of charge jj at polarizability imp
c             to total field of charges at polarizability imp
c
            dmin3 = dmind1**3
c
            do 20, k = 1, 3
              zfpi(k,igranj) = zfpi(k,igranj) + zj*dmin3*pq(k)*factf
   20       continue
c     3-----
          endif
c
          if (iclintd .eq. 1) then
c     3-----
            if (poli.and.polj) then
c       4-----

c        -----  calculate approximate dispersion between
c               classical entities via slater-kirkwood formula
c
              if ((igrppol .eq. 0) .and.
     1           (namei(1:5) .eq. 'group') .and.
     2           (namej(1:5) .eq. 'group')) goto 25
c
              fsk = sqrt(polar(imp)/vale(imp)) +
     1              sqrt(polar(jmp)/vale(jmp))
c
c        -----  isotropic or non-isotropic dispersion may be used
c
              if (isodis .eq. 1) then
c         5-----
c          -----  isotropic polarisabilities
c
                dmin6 = (dmin3*factf)**2
                discon = onept5*polar(imp)*polar(jmp)*dmin6
c         5-----
              else
c         5-----
c          -----  non-isotropic polarisabilities
c          -----  calculate interaction tensor bq
c
                call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c          -----  multiply polarisability tensors: bq = aj bq ai
c
                call tenmul(b,b,b2,3)
c
                if (updpoli) then
                  call hexpand(thispol,polten,3,3)
                else
                  call hexpand(grpol(1,imp),polten,3,3)
                endif
c
                call tenmul(b2,polten,b,3)
c
                if (updpolj) then
                  call hexpand(thispol,polten,3,3)
                else
                  call hexpand(grpol(1,jmp),polten,3,3)
                endif
c
                call tenmul(polten,b,b2,3)
c
c          -----  calculate trace of product tensor bq
c
                call trace(b2,discon,3)
                discon = discon / four
c         5-----
              endif
c
              edisp(indxen) = edisp(indxen) - discon / fsk
 25           continue
c       4-----
            else if (igrppol.eq.0) then
c       4-----
c           use atom polarisabilities

              dmin3 = dmind1**3
c
              if ((atpol(ii) .ne. zero) 
     1    .and. (atpol(jj) .ne. zero)) then
c         5-----
                fsk = sqrt(atpol(ii)/valat(ii)) +
     1                sqrt(atpol(jj)/valat(jj))
c
c          -----  isotropic or non-isotropic dispersion may be used
c
                if (isodis .eq. 1) then
c           6-----
c            -----  isotropic polarisabilities
c
                  dmin6 = (dmin3*factf)**2
                  discon = onept5*atpol(ii)*atpol(jj)*dmin6
c           6-----
                else
c           6-----
c            -----  non-isotropic polarisabilities
c            -----  calculate interaction tensor bq
c
                  call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c            -----  multiply polarisability tensors: bq = aj bq ai
c
                  call tenmul(b,b,b2,3)
                  call hexpand(atpolt(1,ii),polten,3,3)
                  call tenmul(b2,polten,b,3)
                  call hexpand(atpolt(1,jj),polten,3,3)
                  call tenmul(polten,b,b2,3)
c
c            -----  calculate trace of product tensor bq
c
                  call trace(b2,discon,3)
                  discon = discon / four
c           6-----
                endif
c
                edisp(indxen) = edisp(indxen) - discon / fsk
c         5-----
              endif
c       4-----
            endif
c     3-----
          endif
c
c    -----  calculate electrostatic interaction
c           between external charge groups
c
          eelst(indxen) = eelst(indxen) + zi*zj*dmind1*factp
c
          if ((iclintr .eq. 1) .and. (dist .le. repcut)) then
c     3-----
c      -----  calculate approximate repulsion
c             between external charge groups, with empirical r12-term
c             from charmm
c
            if ((namei(:2) .ne. 'xx') .and.
     1          (namei(:2) .ne. 'qq') .and.
     2          (namei(:2) .ne. 'e ') .and.
     3          (namei(:2) .ne. ' ') .and.
     4          (namei(:2) .ne. 'gr') .and.
     5          (namej(:2) .ne. 'xx') .and.
     6          (namej(:2) .ne. 'qq') .and.
     7          (namej(:2) .ne. 'e ') .and.
     8          (namej(:2) .ne. ' ') .and.
     9          (namej(:2) .ne. 'gr')) then
c       4-----
c             alfi = alfa(namei,0,ier)
              alfi = alfext(ii)
              call drfnval(namei,nvali,znuc)
              aoverni = sqrt(alfi/nvali)
c
c             alfj = alfa(namej,0,ier)
              alfj = alfext(jj)
              call drfnval(namej,nvalj,znuc)
              aovernj = sqrt(alfj/nvalj)
c
              ri = radext(ii)*rfact
              rj = radext(jj)*rfact
c
c          -----  account for h-bonding
c
              if ((ihbond .eq. 1)
     1            .and. (dist .le. hbondl)) then
c         5-----
                if(namej(:2).eq.'h'.and.
     1            (namei(:2).eq.'n'.or.
     2             namei(:2).eq.'o'.or.
     3             namei(:2).eq.'f')) rj=hbondr
                if(namei(:2).eq.'h'.and.
     1            (namej(:2).eq.'n'.or.
     2             namej(:2).eq.'o'.or.
     3             namej(:2).eq.'f')) ri=hbondr
c         5-----
              endif
c
c        -----  calculate model repulsion
c
              fac1 = aoverni + aovernj
              fac2 = pt75*alfi*alfj
              fac3 = (ri+rj)**6
              fac4 = dmind1**12
              erep(indxen) = erep(indxen) + fac2*fac3*fac4/fac1
c       4-----
            endif
c     3-----
          endif
c
   40     if(polj) jmp=jmp+1
c   2-----
   90   continue
c
c  -----  store zfpi in zfp
c
        if (poli) then
c   2-----
c    -----  field at polarisability imp due to analysis groups
c
          do 60, l = 1, ngran
            do 60 k = 1, 3
              zfp(irow+k,l) = zfp(irow+k,l) + zfpi(k,l)
              vrp(irow+k,l) = vrp(irow+k,l) - zfpi(k,l)
   60     continue
c
          imp=imp+1
c   2-----
        endif
c 1-----
  100 continue
c
c-----  correct energies for double counting and calculate
c       total classical energy
c
      eclas = zero
      do 200, i = 1, ngrpair
c 1-----
        eelst(i) = pt5*eelst(i)
        if (iclinte .eq. 0) eelst(i) = zero
        edisp(i) = pt5*edisp(i)
        erep(i) = pt5*erep(i)
c
        egrcls(i) = eelst(i) + edisp(i) + erep(i)
c
        eclas = eclas + egrcls(i)
c 1-----
  200 continue
c
      return
      end
      subroutine zfbeupd(ieps,ineq,zfbe,vrbe,xsurf,xnorm,area)
c------
c             updates potentials and fields on boundary from
c             classical point charges
c
c     ndima: dimension of polarization matrix (if present)
c     ndimb: dimension of boundary matrix.
c            =   nbem if kappa.eq.0.0
c            = 2*nbem if kappa.ne.0.0
c     ndim: dimension of complete relay matrix
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  dummy arrays
c
      dimension zfbe(ndim,ngran),vrbe(ndim,ngran)
      dimension area(nbem)
      dimension xsurf(3,nbem),xnorm(3,nbem)
c
c-----  common blocks
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
      integer neq2eps, momsav
      common /neqpar1/ neq2eps,momsav
c
      real*8 kapneq1, kapneq2, epsneq1, epsneq2
      common /neqpar2/ kapneq1,kapneq2,epsneq1,epsneq2
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
c
c-----  local arrays
c
      dimension pi(3), pin(3)
c
c-----  local variables
c
      logical kapnoz
c
      data zero,one,two,four /0.0d00,1.0d00,2.0d00,4.0d00/
c
c-----  begin
c
c-----  set some variables
c
      if (ineq .eq. 0) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
          epsilon = eps1
          kappa = kappa1
        else
          epsilon = eps2
          kappa = kappa2
c   2-----
        endif
c 1-----
      else
c 1-----
        if (ieps .eq. 0) then
c   2-----
          epsilon = epsneq1
          kappa = kapneq1
        else
          epsilon = epsneq2
          kappa = kapneq2
c   2-----
        endif
c 1-----
      endif
c
      pie = four*atan(one)
      epsfact = one/(two*pie*(one+epsilon))
      kapnoz = kappa .ne. zero
      expkd = one
      fkapone = one
      expkdn = one
      fkaponn = one
c
c-----  loop over boundary elements -ni-
c
      do 100, ni = 1, nbem
c 1-----
c  -----  loop over updated external charges -np-
c
        do 200, nn = 2, nnpts
c   2-----
c    -----  pointer to updated point charge
c
          np = inpt(nn)
          igranp = igranl(np)
c
c    -----  pi:  vector from -np- to -ni- = (i-p)
c           dist:  length of pi
c    -----  pin:  updated vector from -np- to -ni- = (i-p)
c           distn:  length of pin
c
          call distab(xpts(1,np),xsurf(1,ni),pi,dist)
          call distab(xnpts(1,nn),xsurf(1,ni),pin,distn)
          dmind1 = one/dist
          dmin3 = (dmind1**2)*dmind1
          dmin1n = one/distn
          dmin3n = (dmin1n**2)*dmin1n
c
          if (kapnoz) then
            expkd = exp(-(kappa*dist))
            fkapone = one + (kappa*dist)
            expkdn = exp(-(kappa*distn))
            fkaponn = one+(kappa*distn)
          endif
c
c    -----  fpini:  field of unit (positive) charge
c                   in -np- at -ni-, contracted with
c                   normal vector in -ni- = (i-p).n(i)/dist**3
c                   = f(p;i) . n(i)
c
c         fpini = dmin3*adotb(pi,xnorm(1,ni),3)
c         fpinni = dmin3n*adotb(pin,xnorm(1,ni),3)
          fpini = dmin3*ddot(3,pi,1,xnorm(1,ni),1)
          fpinni = dmin3n*ddot(3,pin,1,xnorm(1,ni),1)
c
c    -----  potential of source charge in -np- at -ni-
c           scaled with 1/(2pi(1+eps)), as input for coupling
c           equations for w(i)
c
          zfbe(ndima+ni,igranp) = zfbe(ndima+ni,igranp)
     1    - epsfact*(chrg(np)*dmind1 - chrg(np)*dmin1n)
c
c    -----  reaction potential energy operator
c           for source charges, containing
c           contribution of unit surface dipoles w(i),
c           already multiplied by source charges
c
c         sum(p) k(i;p)s(i)*q(p)
c
c           the energy is evaluated in subr drfomga,
c           by contracting with the induced dipole density array
c
c    -----  the minus sign is a result of the use of the
c           inverted field: needed is f(i;p), fpini=f(p;i) . n(i)
c
          vrbe(ndima+ni,igranp) = vrbe(ndima+ni,igranp) +
     1         (epsilon*fkapone*expkd - one)*fpini*area(ni)*chrg(np)
     2       - (epsilon*fkaponn*expkdn - one)*fpinni*area(ni)*chrg(np)
c
c    -----  if non-zero ionic strength, a set of equations
c           for z(i) is added
c
          if (kapnoz) then
c     3-----
c      -----  minus field of source charge in -np- at -ni-,
c             contracted with normal vector at -ni-,
c             scaled with eps/(2pi(1+eps)), as input for
c             coupling equations for z(i)
c
            zfbe(ndima+nbem+ni,igranp) = zfbe(ndima+nbem+ni,igranp) +
     1             epsilon*epsfact*chrg(np)*fpini
     2           - epsilon*epsfact*chrg(np)*fpinni
c
c      -----  reaction potential energy operator
c             for source charges, containing
c             contribution of unit surface charges z(i),
c             already multiplied by source charges
c
c         sum(p) l(i;p)s(i)*q(p)
c
c             the energy is evaluated in subr drfomga,
c             by contacting with the induced charge density array
c
            vrbe(ndima+nbem+ni,igranp) = vrbe(ndima+nbem+ni,igranp) -
     1            (one - expkd)*dmind1*area(ni)*chrg(np)
     2          + (one - expkdn)*dmin1n*area(ni)*chrg(np)
c     3-----
          endif
c   2-----  next charge
  200   continue
c 1-----  next boundary element
  100 continue
c
      return
      end
      subroutine wtvrupd(xexp,wt,vr)
c------
c       calculation of (expanded) source and reaction fields
c       of/on (formal) qm particles (nuclei, electrons)
c       or representations (dp charges, mulliken charges
c       and dipoles) thereof, at polarizabilities,
c       boundary elements and external charges.
c
c       update only when translation is allowed
c
c       --------  p.th. van duijnen, ibm-kingston 1985, and
c                 groningen, dec. 1991.
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
c-----  dummy arrays
c
      dimension xexp(3,nexp),wt(nwtr,nwtc),
     1          vr(nwtr,nwtc)
c
c-----  common blocks
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
      character*8 bflab
      common/mollab/mwanam(128),mwbnam(128),bflab(1024)
c
      character*8 zruntp, zguess, zconf, zstate
      character*8 zcom,title,ztagg,zsymm
      character*8 zorb,zpseud
      character*8 anam
      character*10 bnam
      common/runlab/zcom(19),title(10),anam(maxat),bnam(maxorb),
     +           ztagg(maxat),zsymm(7),zorb(maxorb),zpseud(maxat)
c
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
c
c-----  local variables
c
      character*16 namj
c
      dimension p(3), q(3), pq(3), w(3)
      dimension b(3,3)
      dimension pn(3), pqn(3), wn(3)
      dimension bn(3,3)
c
      character*8 errmsg(3)
c
      data errmsg/'program',' stop in','-wtvrupd'/
      data third/.33333333333333333333333333333d00/
      data sixth/.16666666666666666666666666667d00/
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,pt5,pt75/0.0d00,0.5d00,0.75d00/
      data one,onept5 /1.0d00,1.5d00/
c
c-----  begin
c
c
c-----  get assets of moved polarisable point
c
      np = mpol(inpt(1))
c
      p(1) = xpts(1,np)
      p(2) = xpts(2,np)
      p(3) = xpts(3,np)
c
      pn(1) = xnpts(1,1)
      pn(2) = xnpts(2,1)
      pn(3) = xnpts(3,1)
c
      if = (inpt(1)-1)*3
c
c-----  loop over the expansion centers
c
      do 500, j = 1, nexp
c 1-----
        if (j .le. nat) then
c   2-----
c    -----  nucleus at expanson centre
c
c    -----  nuclear charge -za-
c
          za = czan(j)
c
c    -----  check whether the qm atom is "ambiguous"
c           if it is, this means it is treated both quantum
c           mechanically and classically.
c           replace the atom polarisability by the classical
c           polarisability.
c
          if (ncutpt(j) .ne. 0) then
c     3-----
            namj = nxcent(ncutpt(j))
            alfj = polar(ncutpt(j))
          else
            namj = anam(j)
c
c      -----  polarisability of atom corresponding to expansion centre
c
c           alfj = alfa(namj,0,ier)
            alfj = alfat(j)
c     3-----
          endif
c
c    -----  distance criterion for scaling interactions
c
          rj = (alfj**third)*afact
c   2-----
        else
c   2-----
c    -----  non-nulcei are given a polarizability 1.
c
          alfj = one
c   2-----
        endif
c
c  -----  pointer to arrays w.r.t. expansion center -j-
c
        jf = (j-1)*4
        nzf = jf*3
c
c  -----  position vector -j- into -q-
c
        q(1) = xexp(1,j)
        q(2) = xexp(2,j)
        q(3) = xexp(3,j)
c
c  -----  calculate old and new distance vectors between -q- and -p-
c
        call distab(q,p,pq,dist)
        call distab(q,pn,pqn,distn)
        dmind1 = one/dist
        dmin3 = dmind1*(dmind1**2)
        dmin1n = one/distn
        dmin3n = dmin1n*(dmin1n**2)
c
        factp = one
        factf = one
        v = one
c
        factpn = one
        factfn = one
        vn = one
c
c  -----  account for penetration effects (optional)
c
        if (modxza .ne. 0) then
c   2-----
          s = (alfj*alfj)**sixth
          v = dist/s
          vn = distn/s
c
          if (ithole .eq. 1) then
            if (v .le. afact) then
              av = v/afact
              factp = av**4 - two*av**3 + two*av
              factf = four*av**3 - three*av**4
            endif
          else if (ithole .eq. 2) then
            au = afact*v
            factp = (one - (pt5*au + one)*exp(-au))
            factf = (one - (pt5*au**2 + au + one)*exp(-au))
c     3-----
          endif
c
          if (ithole .eq. 1) then
            if (vn .le. afact) then
              avn = vn/afact
              factpn = avn**4 - two*avn**3 + two*avn
              factfn = four*avn**3 - three*avn**4
            endif
          else if (ithole .eq. 2) then
            au = afact*vn
            factpn = (one - (pt5*au + one)*exp(-au))
            factfn = (one - (pt5*au**2 + au + one)*exp(-au))
          endif
c   2-----
        endif
c
        factp = factp*dmind1
        factf = factf*dmin3
        factpn = factpn*dmin1n
        factfn = factfn*dmin3n
c
c  -----  calculate field gradient tensors
c
        call clear(b,9)
        call clear(bn,9)
c
c  -----  b(3,3) = t(p;q) : field gradient of charge in -p- at -q-
c         note that the interaction is scaled by v
c
        call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
        call drftpq(pqn,dmin1n,dmin3n,ithole,afact,vn,bn)
c
c  -----  calculate -w- = t(p;q) . q , this is part of the
c         taylor expansion of the field
c
        call matvec(b,q,w,3,.false.)
        call matvec(bn,q,wn,3,.false.)
c
c  -----  -wt- and -vr- matrices for expanded field and potential
c         of charge (distribution) in expansion centra
c         at the polarizable points (-wt-) and vice versa (-vr-)
c
c         depending on the type of sources (charges, dipoles,
c         charge distributions), specified in rfin by ifldin
c         and "recipients", specified by ifldout, the matrices
c         are constructed partly or completely.
c
        do 290, k = 1, 3
c   2-----
c    -----  scale the distance vector
c
          pqk = pq(k)*factf
          pqkn = pqn(k)*factfn
c
          do 280, l = 1, 3
c     3-----
c      -----  copy -b- into -wt-  and/or -vr-
c
c             this expansion can be used to calculate the field and
c             reaction field of a unit dipole in -q- (or -j-, the
c             expansion centre) at polarizable point -p-.
c             therefore, it is always calculated except when only
c             distributed monopoles are used to expand the source
c             and reaction fields of the quantum motif
c
c      -----  note: drftpq gives - d/dq f(q;p) !!
c
            if (ifldin .gt. 1) then
              wt(if+k,jf+l) = wt(if+k,jf+l) - bn(k,l) + b(k,l)
            endif
            if (ifldout .gt. 1) then
              vr(if+k,jf+l) = vr(if+k,jf+l) + bn(k,l) - b(k,l)
            endif
c     3-----
  280     continue
c
c    -----  if the source/reaction field is not expanded (i.e. only
c           distributed mono- /dipoles are used), only the potential
c           part is stored in -wt-/-vr-
c
          if (ifldin .gt. 2) then
            wt(if+k,jf+4) = wt(if+k,jf+4) + (wn(k)+pqkn)
     1                                    - (w(k)+pqk)
          else
            wt(if+k,jf+4) = wt(if+k,jf+4) + pqkn - pqk
          endif
c
          if (ifldout .gt. 2) then
            vr(if+k,jf+4) = vr(if+k,jf+4) - (wn(k) + pqkn)
     1                                    + (w(k) + pqk)
          else
            vr(if+k,jf+4) = vr(if+k,jf+4) - pqkn + pqk
          endif
c
c    -----  form -zfn-: the fields in the polarizable points
c           due to the internal nuclei
c           -zfn- is in fact the sum of the nuclear fields
c           and potentials
c
          if (j .le. nat) then
            wt(if+k,nexp4+ngran+1) = wt(if+k,nexp4+ngran+1)
     1                              + (pqkn - pqk)*za
            vr(if+k,nexp4+ngran+1) = vr(if+k,nexp4+ngran+1)
     1                              - (pqkn - pqk)*za
          endif
c   2-----
  290   continue
c 1-----  next expansion centre
  500 continue
c
      return
      end
      subroutine expupd(xexp,extnucn,epseun,
     *  extncgn,repmdgn,istat)
c------
c      updates the (expanded) field due to
c      external charges at the expansion centra of the
c      quantum mechanical system
c
c      the interaction energy between the external charges
c      and the nuclei (extnucn) is also updated
c
c      penetration (overlap) effects
c      resulting from close contacts
c      can be accounted for (modxza flag)
c
c      --------  p.th. van duijnen, ibm-kingston 1985, and
c               groningen, dec. 1991.
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
c
      dimension xexp(3,nexp)
c
c
      dimension extncgn(mxgran,mxst) ,repmdgn(mxgran,mxst)
c
c-----  common blocks
c
c
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
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
      character*8 bflab
      common/mollab/mwanam(128),mwbnam(128),bflab(1024)
c
      character*8 zruntp, zguess, zconf, zstate
      character*8 zcom,title,ztagg,zsymm
      character*8 zorb,zpseud
      character*8 anam
      character*10 bnam
      common/runlab/zcom(19),title(10),anam(maxat),bnam(maxorb),
     +           ztagg(maxat),zsymm(7),zorb(maxorb),zpseud(maxat)
c
c
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep
      integer neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga
      integer idipcal,iclinte,iclintd,ieffpol,isodis,iclintr
      integer iextdip,itsolv,itermax,isolsav,imomsav,idisadd
      integer ianal,ineqex,ithole,irevdis,ihbond,igrppol
      integer maxci_drf
      common/drfpar1/ibem,nodiscr,ixrelay,ialfso,neqsta,neqdis,neqrep,
     +               neqrf,ixamat,ixbmat,ixzfp,ixexp,ixwtvr,ixomga,
     +               idipcal,iclinte,iclintd,ieffpol,isodis,iclintr,
     +               iextdip,itsolv,itermax,isolsav,imomsav,idisadd,
     +               ianal,ineqex,ithole,irevdis,ihbond,igrppol,
     +               maxci_drf
c
      integer nxtpts, npol, nexp, natint, nshint, namb, nspec, ngran
      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
c
      integer ndima, ndimb, nbem, ndim, nexp4, npol3, nwtr, nwtc, nzfa
      integer lomga, nomga, nchd, ngrnam, ngrnam2, nzfp, nzfn
      integer nneq, nneqrf, maxneq, neqdim, nqdim, nqcls
      common/drfpar3/ndima,ndimb,nbem,ndim,nexp4,npol3,nwtr,nwtc,nzfa,
     +               lomga,nomga,nchd,ngrnam,ngrnam2,nzfp,nzfn,
     +               nneq,nneqrf,maxneq,neqdim,nqdim,nqcls
c
      integer ifldin, ifldout, idrfout, modxza, irepopt
      integer iexpza, icmexp, iadexp, nodpe, igetden, iarfcal, iqmclr
      common/drfpar4/ifldin,ifldout,idrfout,modxza,irepopt,
     +               iexpza,icmexp,iadexp,nodpe,igetden,iarfcal,
     +               iqmclr
c
      real*8 gamdrf, dstmin, dstmax, hbondl, hbondr, afact, rfact
      real*8 cvgrel, agrpe, agrpm, agrpc, scffact, acur
      real*8 conci_drf
      common/drfpar5/gamdrf,dstmin,dstmax,hbondl,hbondr,afact,rfact,
     +               cvgrel,agrpe,agrpm,agrpc,scffact,acur,conci_drf
c
      integer iind, isur1, isur2, ilwt, ilvr, illur, ilindx
      common/drfpar6/iind,isur1,isur2,ilwt,ilvr,illur,ilindx
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 zfa
      common/drfzfa1/zfa(4*(mxat +1),mxgran)
      real*8 zpa
      common/drfzfa2/zpa(mxat +1,mxgran)
      real*8 dzfa
      common/derzfa1/dzfa(12*(mxat+1),mxgran)
      real*8 dpseu
      common/derzfa2/dpseu(3,(mxat+1),mxgran)
c
c
      real*8 xpts
      common /extpts/ xpts(3,mxpts)
c
      character*16 nxcent
      common /namext/ nxcent(mxpts)
c
      real*8 chrg
      common /extcha/ chrg(mxpts)
c
      real*8 polar
      common /extpol/ polar(mxpol)
c
      real*8 atpol
      common /extpol2/ atpol(mxpts)
c
      real*8 atpolt
      common /extpolt/ atpolt(6,mxpts)
c
      real*8 alfext
      common /extalf/ alfext(mxpts)
c
      integer mpol
      common /ipol/ mpol(mxpol)
c
      integer ncutpt
      common /drfcut/ ncutpt(mxpts)
c
      integer nspecl
      common /drfspe/ nspecl(mxspe)
c
      integer nvalel
      common /drfval/ nvalel(mxpts)
c
      real*8 vale
      common /clasval/ vale(mxpts)
c
      real*8 valat
      common /atval/ valat(mxpts)
c
      real*8 polars
      common /claspol/ polars(mxpts)
c
      real*8 alfat
      common /atalf/ alfat(mxat)
c
c
      integer ngrpol,igrpol
      common /polgrp/ ngrpol,igrpol(mxgrp+1)
c
      real*8 grpol
      common /grpols/ grpol(6,mxgrp+1)
c
      integer nvalgrp
      common /grpnval/ nvalgrp(mxgrp+1)
c
      integer igranl
      common /grpanl/ igranl(mxpts)
c
      integer imemgrp
      common /grpmem/ imemgrp(mxpts)
c
      integer neqgrp
      common /grpneq/ neqgrp(mxgran)
c
      character*80 grlabel
      common /grlab/ grlabel(mxgrp)
c
c
c
      integer itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
      common /drfbem1/ itwoeps,ioponly,itwosur,ixnsurf,ibemout,iuniout
c
      integer leveli, levelo, nbem1, nbem2, nraw
      common /drfbem2/ leveli,levelo,nbem1,nbem2,nraw
c
      character*32 solnam
      real*8 radmax, eps1, eps2, spherad, rprobe, rprobej
      real*8 solrad, swidth, sdist
      real*8 kappa
      real*8 kappa1, kappa2, kappas
      common /drfbem3/ radmax,eps1,eps2,kappa1,kappa2,spherad,rprobe,
     +                 rprobej,solnam,solrad,swidth,sdist
c
      real*8 radat, cgrav
      common /drfbem4/ radat(maxat),cgrav(3)
c
c
      real*8 radext
      common /drfrad/ radext(mxpts)
c
c
c
      integer nsamp, nblock, nmoves, ncheck
      common /mcrunp1/ nsamp,nblock,nmoves,ncheck
c
      real*8 darot, dtrans, temp, onekt, excld, delmax, enmin
      common /mcrunp2/ darot,dtrans,temp,onekt,excld,delmax,enmin
c
      logical mcupdt,secstat,excess
      integer notrans, norot, imcout, iseed, imcst, iactst
      common /mcrunp3/ notrans,norot,imcout,iseed,imcst,iactst,
     +               mcupdt,secstat,excess
c
      character*16 outfor
      common /mcrunp4/ outfor
c
      real*8 ratmin, ratmax, amxrot, amxtrn, amnrot, amntrn
      common /mcrunp5/ ratmin,ratmax,amxrot,amxtrn,amnrot,amntrn
c
      real*8 gammc
      common /mcrunp6/ gammc(5)
c
c
      real*8 xnpts
      common /mcmod1/ xnpts(3,mxnpts)
c
      integer ngrpmc, nnpts, inpt
      common /mcmod2/ ngrpmc,nnpts,inpt(mxnpts)
c
      integer ngrppt
      common /mcmod3/ ngrppt(mxnpts)
c
      integer igrpst
      common /mcmod4/ igrpst(mxgrp+1)
c
      integer igrlast
      common /mcmod5/ igrlast(mxgrp+1)
c
      integer ibitmc
      common /mcmod6/ ibitmc(mxgrp+1)
c
c
c-----  local variables
c
      logical ipol
c
      dimension p(3), q(3), pq(3), w(3)
      dimension pqn(3), wn(3)
      dimension b(3,3), bn(3,3)
c
      character*16 namj, nami
c
      character*8 errmsg(3)
c
      data errmsg/'program',' stop in','-expupd-'/
      data third/.33333333333333333333333333333d00/
      data sixth/.16666666666666666666666666667d00/
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,one,pt5,pt75/0.0d00,1.0d00,0.5d00,0.75d00/
      data small /1.0d-03/
c
c-----  begin
c
c-----  update zfa, dzfa (derivative) and pseudo repulsion
c
      pseufac=pt75
c
c-----  loop over the expansion centra
c
      do 500, j = 1, nexp
c 1-----
        if (j .le. nat) then
c   2-----
          za = czan(j)
          namj = anam(j)
c         alfj = alfa(namj,0,ier)
          alfj = alfat(j)
c
c    -----  # of  valence electrons
c
          call drfnval(namj,nvalj,znuc)
c
c    -----  distance criterium for scaling interactions
c
          if (alfj .ne. zero) then
            aovernj = sqrt(alfj/nvalj)
          endif
c   2-----
        else
c   2-----
c    -----  non-nulcei are given a polarizability 1.
c
          alfj = one
c   2-----
        endif
c
c  -----  pointers to arrays w.r.t. expansion center -j-
c
        jf = (j-1)*4
        nzf = jf*3
c
c  -----  position vector -j- into -q-
c
        q(1) = xexp(1,j)
        q(2) = xexp(2,j)
        q(3) = xexp(3,j)
c
c  -----  pointer to external points in -xtpts-
c         with special properties
c         imp: polarizable external points
c
        imp=1
c
c  -----  loop over updated external points
c
        do 300, in = 2, nnpts
c   2-----
c
c    -----  check whether external point is polarizable
c
          ipol = ii .eq. mpol(imp)
c
c    -----  check whether external point is an ambiguous atom.
c           The interaction with the qm system is taken care
c           of by inclusion into the qm system.
c
          if (ncutpt(ii).ne.0) goto 290
c
c    -----  get position of old point
c
          ii = inpt(in)
          igrani = igranl(ii)
c
c    -----  position vector of external point -ii- in -p-
c
          p(1)=xpts(1,ii)
          p(2)=xpts(2,ii)
          p(3)=xpts(3,ii)
c
          call distab(q,p,pq,dist)
          call distab(q,xnpts(1,in),pqn,distn)
c
c
          if ((dist .lt. small) .or. (distn .lt. small)) then
              write(iwr,95) anam(j), nxcent(ii), dist
  95          format(/,1x,'WARNING: interaction ',
     1   'between ', a16, ' and ', a16, ' at', e15.8,
     2   ' bohr distance: skipped')
              goto 300
          endif
c
          alfi = alfext(ii)
          if (alfi .eq. zero) alfi = one
          dmind1 = one/dist
          dmin3 = dmind1*(dmind1**2)
          factp = one
          factf = one
          v = one
          dmin1n = one/distn
          dmin3n = dmin1n*(dmin1n**2)
          factpn = one
          factfn = one
          vn = one
c
c    -----  account for penetration effects (optional)
c
          if (modxza .ne. 0) then
c     3-----
c      -----  set distance criterium for penetration effects
c             depending on polarizabilities of -j- and -ii-
c
            s = (alfj*alfi)**sixth
            v = dist/s
            vn = distn/s
c
c      -----  check distance criterium
c             modify interactions (thole's cone model)
c
            if (ithole .eq. 1) then
              if (v .le. afact) then
                av = v/afact
                factp = av**4 - two*av**3 + two*av
                factf = four*av**3 - three*av**4
              endif
            else if (ithole .eq. 2) then
              au = afact*v
              factp = (one - (pt5*au + one)*exp(-au))
              factf = (one - (pt5*au**2 + au + one)*exp(-au))
            endif
            factp = factp*dmind1
            factf = factf*dmin3
c
            if (ithole .eq. 1) then
              if (vn .le. afact) then
                avn = vn/afact
                factpn = avn**4 - two*avn**3 + two*avn
                factfn = four*avn**3 - three*avn**4
              endif
            else if (ithole .eq. 2) then
              au = afact*vn
              factpn = (one - (pt5*au + one)*exp(-au))
              factfn = (one - (pt5*au**2 + au + one)*exp(-au))
            endif
            factpn = factpn*dmin1n
            factfn = factfn*dmin3n
c     3-----
          endif
c
c    -----  update coulomb interaction with internal nuclei
c
          if (j .le. nat) then
            extnuca = za*chrg(ii)*(factpn - factp)
            extnucn = extnucn + extnuca
            extncgn(igrani,istat) = extncgn(igrani,istat) + extnuca
          endif
c
c    -----  update -zfa-
c
c             -zfa- : sum(p) {f(p;q) + v(p;q) - f(p;q).q} * q(p)
c
c             the field operator and potential of all external charges
c             in p (-p-) at expansion centre q (-j- = -q-)
c
          sum = zero
          do 210, k = 1, 3
            zfa(jf+k,igrani) = zfa(jf+k,igrani)+chrg(ii)*
     1                  (pqn(k)*factfn - pq(k)*factf)
            sum = sum - q(k)*(pqn(k)*factfn -  pq(k)*factf)
  210     continue
          zfa(jf+4,igrani) = zfa(jf+4,igrani)
     1                      + chrg(ii)*(sum+(factpn - factp))
c
c    -----  update -zpa-
c
          zpa(j,igrani) = zpa(j,igrani) + chrg(ii)*(factpn - factp)
c
c    -----  classical field gradient tensor into -b-
c
          call clear(b,9)
          call clear(bn,9)
c
c    -----  b(3,3) = t(p;q) : field gradient of charge in -p- at -q-
c           the gradient is scaled according to v
c
          call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
          call drftpq(pqn,dmin1n,dmin3n,ithole,afact,vn,bn)
c
c    -----  multiply position vector of -j- with -b- to give -w-
c
          call matvec(b,q,w,3,.false.)
          call matvec(bn,q,wn,3,.false.)
c
c    -----  update derivative of -zfa-
c
          do 214, k = 1, 3
            do 213, l = 1, 3
              dzfa(nzf+(k-1)*4+l,igrani) = dzfa(nzf+(k-1)*4+l,igrani)
     1                - (bn(k,l) - b(k,l))*chrg(ii)
  213       continue
            dzfa(nzf+(k-1)*4+4,igrani) = dzfa(nzf+(k-1)*4+4,igrani)
     1                  +(wn(k) - w(k))*chrg(ii)
  214     continue
c
c    -----  non bonding repulsion taken from charmm:
c           j.comp.chem. 4 (1983) 213
c
c    -----  this is called the pseudo-repulsion throughout
c           the program
c
          if ((field(:4) .ne. ' ') .and.
c     3-----
     1          (j .le. nat) .and.
     1          (iqmclr .eq. 1) .and. 
     1          (namj(:2) .ne. 'bq') .and.
     1       (nxcent(ii)(:2) .ne. 'qq') .and.
     2       (nxcent(ii)(:2) .ne. 'e ') .and.
     3       (nxcent(ii)(:2) .ne. '  ')) then
c
c      -----  this might be a close atom: apply distance
c             criterium for pauli repulsion
c
            if (dist .lt. dstmin) then
c       4-----
c        -----  calculate pauli repulsion and
c               derivative, using some model potential
c
c               the necessary ingredients:
c
c        -----  # valence electrons
c
              nami=nxcent(ii)
c
c        -----  check if the point represents a group
c
              if (nami(:5).eq.'group') then
c           
c        -----  non-bonded interactions evaluated only between
c               internal atoms and individual members of a group:
c               skip the group-representing point
c
                goto 290
              else
c
c        -----  get the number of valence electrons of the
c               external point
c
                call drfnval(nxcent(ii),nvali,znuc)
c
c        -----  check if it was already found
c               close to another internal atom
c
                do 220, nsp = 1, nspec
                  if (nspecl(nsp) .eq. ii) goto 230
  220           continue
c
c       -----  this external point has not been found
c              close to an internal atom:
c              add it to the list
c
                nspec = nspec+1
                nvalel(nspec)=nvali
c               polars(nspec) = alfa(nxcent(ii),0,ier)
                polars(nspec) = alfext(ii)
                nspecl(nspec)=ii
  230           continue
c
c        -----  1/r12 term for selected atoms
c
c        -----  get polarizability of external atom
c
                if (ipol) then
                  alfi = polar(imp)
                else
c                 alfi = alfa(nami,0,ier)
                  alfi = alfext(ii)
                endif

                if (alfi .ne. zero) then
                  aoverni = sqrt(alfi/nvali)
                endif
c
c        -----  calculate equilibrium lj distance from
c               "size" of internal and external atom
c
c        -----  size of external atom: ri
c
                ri = radext(ii)*rfact
c
c        -----  size of external atom: rjx
c
                rjx = radat(j)*rfact
c
c        -----  in case of a h-bond, the sizes need to be
c               modified
c
c        -----  check h-bond criterium
c
                if ((ihbond .eq. 1)
     1            .and. (dist .le. hbondl)) then
c         5-----
                  if(namj(:2).eq.'h'.and.
     1                  (nami(:2).eq.'n'.or.
     2                   nami(:2).eq.'o'.or.
     3                   nami(:2).eq.'f')) rjx=hbondr
                  if(nami(:2).eq.'h'.and.
     1                  (namj(:2).eq.'n'.or.
     2                   namj(:2).eq.'o'.or.
     3                   namj(:2).eq.'f')) ri=hbondr
c         5-----
                endif
c
c        -----  calculate the 1/r12 lj term for the close contact
c
                if ((alfj .ne. zero)
     1            .and. (alfi .ne. zero)) then
c         5-----
                  pseu = pseufac*alfj*alfi/(aovernj+aoverni)
                  pseu1 = (ri+rjx)**6
                  epseun = epseun - pseu*pseu1*(dmind1**12)
                  repmdgn(igrani,istat) = repmdgn(igrani,istat)
     1          - pseu*pseu1*(dmind1**12)
c
c          -----  and derivative
c
                  do 270 k=1,3
                    dpseu(k,j,igrani) = dpseu(k,j,igrani) -
     1                     twelve*pq(k)*pseu*pseu1*(dmind1**14)
  270             continue
c         5-----
                endif
c       4-----
              endif
c
              if (distn .lt. dstmin) then
c       4-----
c        -----  calculate pauli repulsion and
c               derivative, using some model potential
c
c               the necessary ingredients:
c
c        -----  # valence electrons
c
                nami = nxcent(ii)
                call drfnval(nxcent(ii),nvali,znuc)
c
c        -----  check if it was already found
c               close to another internal atom
c
                do 320, nsp = 1, nspec
                  if (nspecl(nsp) .eq. ii) goto 330
  320           continue
c
c        -----  this external point has not been found
c               close to an internal atom:
c               add it to the list
c
                nspec = nspec + 1
                nvalel(nspec) = nvali
                nspecl(nspec) = ii
  330           continue
c
c        -----  1/r12 term for selected atoms
c
c               alfi = alfa(nami,0,ier)
                alfi = alfext(ii)
                if (alfi .ne. zero) then
                  aoverni = sqrt(alfi/nvali)
                endif
c
c        -----  calculate equilibrium lj distance from
c               "size" of internal and external atom
c
c        -----  size of external atom: ri
c
                rin = radext(ii)*rfact
c
c        -----  size of external atom: rjx
c
                rjxn = radat(j)*rfact
c
c        -----  in case of a h-bond, the sizes need to be
c               modified
c
c        -----  check h-bond criterium
c
                if ((ihbond .eq. 1) 
     1            .and. (distn .le. hbondl)) then
                  if(namj(:2).eq.'h'.and.
     1                  (nami(:2).eq.'n'.or.
     2                   nami(:2).eq.'o'.or.
     3                   nami(:2).eq.'f')) rjxn=hbondr
                  if(nami(:2).eq.'h'.and.
     1                  (namj(:2).eq.'n'.or.
     2                   namj(:2).eq.'o'.or.
     3                   namj(:2).eq.'f')) rin=hbondr
                endif
c
c        -----  calculate the 1/r12 lj term for the close contact
c
                if ((alfj .ne. zero) 
     1      .and. (alfi .ne. zero)) then
c         5-----
                  pseu = pseufac*alfj*alfi/(aovernj+aoverni)
                  pseu1 = (rin+rjxn)**6
                  epseun = epseun + pseu*pseu1*(dmin1n**12)
                  repmdgn(igrani,istat) = repmdgn(igrani,istat)
     1          +   pseu*pseu1*(dmin1n**12)
c
c          -----  and derivative
c
                  do 280, k = 1, 3
                    dpseu(k,j,igrani) = dpseu(k,j,igrani) +
     1                 twelve*pqn(k)*pseu*pseu1*(dmin1n**14)
  280             continue
c         5-----
                endif
c       4-----
              endif
c     3-----
            endif
c
          endif
c
  290     if (ipol) imp=imp+1
c
c   2-----  next new point
c
  300   continue
c
c 1-----  next expansion centre
c
  500 continue
c
      if (imcout .eq. 5) then
        do 700, igr = 1, ngran
          call hatout(dzfa(1,igr),4,3*nexp,5,'dzfa-upd')
          call hatout(zfa(1,igr),4,nexp,5,'zfa-upd')
          call hatout(zpa(1,igr),1,nexp,5,'zpa-upd')
  700   continue
      endif
c
      return
      end
