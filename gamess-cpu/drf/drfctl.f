      subroutine rfcal(xscm)
c------
c      this routine drives the calculation of all
c      quantities needed for the evaluation of the reaction field
c
c      these are (in order of appearance:)
c
c      - the second order perturbation polarizability of the
c      quantum mechanically treated system (alfsop)
c
c      - the coordinates, normal vectros and areas of the
c      boundary elements, if a boundary is specified (drfsurf)
c
c
c      - relay-matrix, if the coupled linear response equations
c      are solved through lu-decomposition of this matrix
c
c      - wt and vr matrices, containing the expanded
c      source and reaction fields, if the iterative method
c      is not used
c
c      - omega, that is the matrix containing (expanded) interactions
c      between the quantum mechanical system and its environment
c      omega is either constructed through iterative solution of
c      the linear response equations (itsolv = .true.),
c      or through lu-decomposing the so-called relay-matrix and
c      contracting this with source and reaction fields/potentials
c      (wt and vr). in the iterative method wt and vr are calculated
c      repeatedly on the fly.
c
c  *   if electrons and nuclei are treated separately for
c  *   either source or reaction field
c  *
c  *   - the assignation of overlap, first and second moment
c  *   integrals to expansion centra (drfass)
c
c  *  only when the source fields of the electrons and nuclei
c  *  are treated separately  a n d  the reaction fields are
c  *  coupled back to the electrons and nuclei separately
c  *
c  *   - the one-electron reaction field integrals (drfone)
c
c
c      it is called from standv
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
c
c-----  general common blocks
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
c
      integer ia
      common /ijpair/ ia(3*mxpts)
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
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
      dimension xscm(*)
c
c-----  local variables
c
      dimension cm(3)
      dimension optneq(7)
      dimension ibemdim(2),xbemdim(2)
      equivalence (ibemdim(1),xbemdim(1))
c
      data zero /0.0d00/
c
c-----  begin
c
      if (.not. mcupdt) write(iwr,11)
   11 format(/1x,104('-')/)
c
      if ((neqsta .eq. 1) .or. (neqrf .eq. 1)) then
c 1-----
        call daread(idafinf,iodainf,dimm,1,ineqix(nneq))
        neqdim = int(dimm)
c
        call daread(idafcst,iodacst,optneq,3,ineqix(nneq))
        ustaneq(iactst) = optneq(1)
        if (neqdis .eq. 1) ustaneq(iactst) = ustaneq(iactst) + optneq(2)
        if (neqrep .eq. 1) ustaneq(iactst) = ustaneq(iactst) + optneq(3)
c 1-----
      endif
c
c-----  default value for nonequilibrium potential: as potential at
c       expansion centra
c
      momsav = 0
c
      if (neqrf .eq. 1) then
c 1-----
c  -----  read parameters for non-equilibrium reaction field
c
        call daread(idafpol,iodapol,optneq,7,ineqix(nneq))
        upolneq(iactst) = optneq(1)
        epsneq1 = optneq(2)
        kapneq1 = optneq(3)
        epsneq2 = optneq(4)
        kapneq2 = optneq(5)
        neq2eps = int(optneq(6))
        momsav = int(optneq(7))
c 1-----
      endif
c
      uneqnuc(iactst) = zero
      ndimb1 = 0
      ndimb2 = 0
      nchd = (num*(num+1))/2
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
c
c-----  save vacuum one electron hamiltonian on drf da file also
c
      if (.not. mcupdt .and. .not. secstat) then
c 1-----
        if ((idrfout .ge. 1) .and. (nneq .le. 1)) then
c   2-----
          write(iwr,1001)
 1001     format(//'------  expanded reaction field --------'//)
c   2-----
        endif
c
c  -----  in mc calculation, save one-electron hamiltonian on
c       record -11- of unit -31- the first time only
c
c       need = loc10 + nchd
c       call setc(need)
	i10 = igmem_alloc(nchd)
        call daread(idafh   ,ioda   ,xscm(i10),nchd,11)
        call dawrit(idafdrf,iodadrf,xscm(i10),nchd,11,navdrf)
	call gmem_free(i10)
c       call setc(loadcm)
c 1-----
      endif
c
c-----  calculate polarizability of qm system if required,
c       by second order perturbation theory
c
      if (ialfso .ne. 0) then
        ifast = 0
        call alfsop(xscm)
        ialfso = 0
      endif
c
      ipass = 0
c
c-----  define boundary surface if required
c
      ieps = 0
   90 continue
c
      if (ibem .ne. 0) then
        if (ixnsurf .ne. 0) then
          call drfsurf(ieps,xscm)
        else if (field(5:) .ne. ' ') then
          m2 = lenint(2)
          call daread(idafdrf,iodadrf,xbemdim,m2,49)
          nbem1 = ibemdim(1)
          nbem2 = ibemdim(2)
        endif
      endif
c
      if ((itwosur .eq. 1) .and. (ieps. eq. 0)) then
        ieps = 1
        goto 90
      endif
c
      npol3 = npol*3
      ndima = npol3
      nnpol3 = npol3*(npol3+1)/2
c
      ieps = 0
      ineq = 0
c
c-----  save dimensions of boundary problem
c
      if (ibem .ne. 0) then
        ibemdim(1) = nbem1
        ibemdim(2) = nbem2
        m2 = lenint(2)
        call dawrit(idafdrf,iodadrf,xbemdim,m2,49,navdrf)
      endif
c
  100 continue
c
      ipass = ipass + 1
c
      nbem = nbem1
      if (ineq .eq. 0) then
        kappa = kappa1
      else
        kappa = kapneq1
      endif
c
      if (ieps .eq. 1) then
c 1-----
        if (itwosur .eq. 1) then
          nbem = nbem2
        endif
c
        if (ineq .eq. 0) then
          kappa = kappa2
        else
          kappa = kapneq2
        endif
c 1-----
      endif
c
      ndimb = nbem
      if (kappa .ne. zero) ndimb = 2*ndimb
c
      nwtr = npol3 + ndimb
      ndim = ndima + ndimb
c
c-----  calculate relay-matrix and lu decompose it, if
c       it isn't done yet, and if required
c
      if (ixrelay .ne. 0) call drfrelay(ieps,xscm)
c
c-----  calculate source and reaction field matrices
c       -wt- and -vr- if required
c
      if (ixwtvr .ne. 0) then
c 1-----
c  -----  set memory partitioning
c
	ixw = igmem_alloc(nwtc*nwtr)
	ixv = igmem_alloc(nwtc*nwtr)
	ixe = igmem_alloc(3*nexp)
c       ixw = i10
c       ixv = ixw + nwtc*nwtr
c       ixe = ixv + nwtc*nwtr
c       ixs = ixe + 3*nexp
c
        if (ibem .ne. 0) then
c   2-----
c    -----  dielectric surroundings included
c
	ixs = igmem_alloc(3*nbem)
	ixn = igmem_alloc(3*nbem)
	ixa = igmem_alloc(nbem)
c         ixn = ixs + 3*nbem
c         ixa = ixn + 3*nbem
c         last = ixa + nbem + 1
c       else
cixs = igmem_alloc(1)
cixn = igmem_alloc(1)
cixa = igmem_alloc(1)
c         ixn = ixs + 1
c         ixa = ixn + 1
c         last = ixa + 1
c   2-----
        endif
c
c       need = loc10 + last - i10
c       call setc(need)
        call drfwt(xscm(ixe),xscm(ixw),xscm(ixv),
     1           xscm(ixs),xscm(ixn),xscm(ixa),ieps,ineq,momsav)
c       call setc(loadcm)
        if (ibem .ne. 0) then
	call gmem_free(ixa)
	call gmem_free(ixn)
	call gmem_free(ixs)
	endif
	call gmem_free(ixe)
	call gmem_free(ixv)
	call gmem_free(ixw)
c 1-----
      endif
c
c-----  calculate the -omega- matrix
c
      if (ixomga .ne. 0) call drfomg(ieps,xscm)
c
c-----  further dimensioning
c
      ixvr = igmem_alloc(max(nwtc*nwtr,neqdim))
	ixo = igmem_alloc(nchd)
	ixdx = igmem_alloc(nchd)
	ixdy = igmem_alloc(nchd)
	ixdz = igmem_alloc(nchd)
c     ixvr = i10
c     ixo = ixvr + max(nwtc*nwtr,neqdim)
c     ixdx = ixo + nchd
c     ixdy = ixdx + nchd
c     ixdz = ixdy + nchd
c     ixrxx = ixdz + nchd
c
      if ((neqsta .eq. 1) .and. (ipass .eq. 1)) then
c 1-----
c  -----  integrals for external electrostatic field from
c         a previous calculation
c
	ixone = igmem_alloc(nchd)
	ixoneh = igmem_alloc(nchd)
	ixie = igmem_alloc(nchd)
c       ixone = ixdz + nchd
c       ixoneh = ixone + nchd
c       ixie = ixoneh
c       last = ixie + nchd + 1
c       need = loc10+ last - i10
c       call setc(need)
c
        call neqon2(1,ieps,xscm(ixvr),xscm(ixo),xscm(ixdx),
     1              xscm(ixdy),xscm(ixdz),xscm(ixone),xscm(ixoneh),
     2              xscm(ixie))
c       call setc(loadcm)
	call gmem_free(ixie)
	call gmem_free(ixoneh)
	call gmem_free(ixone)
c 1-----
      endif
c
c-----  if non-equilibrium reaction field, calculate integrals
c       only for static dielectric contribution
c
      if (ineq .eq. 1) then
c 1-----
	ixone = igmem_alloc(nchd)
	ixoneh = igmem_alloc(nchd)
	ixind = igmem_alloc(nwtr)
	ixie = igmem_alloc(nchd)
c       ixone = ixdz + nchd
c       ixoneh = ixone + nchd
c       ixind = ixoneh + nchd
c       ixie = ixind + nwtr
c       last = ixie + nchd + 1
c       need = loc10+ last - i10
c       call setc(need)
c
        call neqone(ieps,momsav,xscm(ixvr),xscm(ixo),xscm(ixdx),
     1              xscm(ixdy),xscm(ixdz),xscm(ixone),xscm(ixoneh),
     2              xscm(ixie),xscm(ixind),xscm(ixvr))
c       call setc(loadcm)
	call gmem_free(ixie)
	call gmem_free(ixind)
	call gmem_free(ixoneh)
	call gmem_free(ixone)
      call gmem_free(ixdz)
      call gmem_free(ixdy)
      call gmem_free(ixdx)
      call gmem_free(ixo)
      call gmem_free(ixvr)
        goto 200
c 1-----
      endif
c
c-----  calculate one-electron reaction field integrals
c
c-----  set memory partitioning
c
c      ovlap at ixo
c      dipx  at ixdx
c      dipy  at ixdy
c      dipz  at ixdz
c      rxx   at ixrxx
c      ryy   at ixryy
c      rzz   at ixrzz
c      rxy   at ixrxy
c      rxz   at ixrxz
c      ryz   at ixryz
c      omega at ixomg
c      xexp  at ixe
c      velq  at ixvel
c      veln  at ixvsne
c      vsnua at ixvsn
c      vself at ixvs
c      vext  at ixvex
c      relay at ixr
c      aind  at ixind
c      iexp  at ixie
c      ijbit at ixij
c      indx  at ixindx
c
c     if (field(5:) .eq. 'scf') then
c 1-----
	ixrxx = igmem_alloc(nchd)
	ixryy = igmem_alloc(nchd)
	ixrzz = igmem_alloc(nchd)
	ixrxy = igmem_alloc(nchd)
	ixrxz = igmem_alloc(nchd)
	ixryz = igmem_alloc(nchd)
	ixomg = igmem_alloc(nomga)
	ixvel = igmem_alloc(ngran*nchd)
	ixvsne = igmem_alloc(nchd)
	ixvsn = igmem_alloc(nchd)
	ixvs = igmem_alloc(nchd)
	ixvex = igmem_alloc(ngran*nchd)
c       ixryy = ixrxx + nchd
c       ixrzz = ixryy + nchd
c       ixrxy = ixrzz + nchd
c       ixrxz = ixrxy + nchd
c       ixryz = ixrxz + nchd
c       ixomg = ixryz + nchd
c       ixe = ixomg + nomga
c       ixvel = ixe + 3*nexp
c       ixvsne = ixvel + ngran*nchd
c       ixvsn = ixvsne + nchd
c       ixvs = ixvsn + nchd
c       ixvex = ixvs + nchd
c       ixr = ixvex + ngran*nchd
c 1-----
c     else
cixrxx = igmem_alloc(1)
cixryy = igmem_alloc(1)
cixrzz = igmem_alloc(1)
cixrxy = igmem_alloc(1)
cixrxz = igmem_alloc(1)
cixryz = igmem_alloc(1)
cixomg = igmem_alloc(1)
cixvel = igmem_alloc(1)
cixvsne = igmem_alloc(1)
cixvsn = igmem_alloc(1)
cixvs = igmem_alloc(1)
cixvex = igmem_alloc(1)
c 1-----
c       ixryy = ixrxx + 1
c       ixrzz = ixryy + 1
c       ixrxy = ixrzz + 1
c       ixrxz = ixrxy + 1
c       ixryz = ixrxz + 1
c       ixomg = ixryz + 1
c       ixe = ixomg + 1
c       ixvel = ixe + 3*nexp
c
c       if (iextdip .eq. 0) then
c   2-----
c         ixvsne = ixvel + 1
c       else
c         ixvsne = ixvel + nchd
c   2-----
c       endif
c
c       ixvsn = ixvsne + 1
c       ixvs = ixvsn + 1
c       ixvex = ixvs + 1
c       ixr = ixvex + ngran*nchd
c 1-----
c     endif
c
      ixe = igmem_alloc(3*nexp)
      ixr = igmem_alloc(ndim*ndim)
c     if (iextdip .eq. 1) then
c 1-----
c       ixind = igmem_alloc(ndim)
c       ixie  = igmem_alloc(nchd)
c       ixij  = igmem_alloc(nchd)
c       ixindx  = igmem_alloc(ndim)
c       ixind = ixr + ndim*ndim
c       ixie = ixind + ndim
c       ixij = ixie + nchd
c       ixindx = ixij + 1
c       last = ixindx + ndim + 1
c     else
c       ixind = igmem_alloc(1)
c       ixie  = igmem_alloc(nchd)
c       ixij  = igmem_alloc(1)
c       ixindx  = igmem_alloc(1)
c       ixind = ixr + 1
c       ixie = ixind + 1
c       ixij = ixie + nchd
c       ixindx = ixij + 1
c       last = ixindx + 1
c 1-----
c     endif
        ixind = igmem_alloc(ndim)
        ixie  = igmem_alloc(nchd)
        ixij  = igmem_alloc(nchd)
        ixindx  = igmem_alloc(ndim)
c
c     need = loc10+ last - i10
c     call setc(need)
c     call clear(xscm(ixvr),last-i10+1)
c
      call drfone(ieps,xscm(ixo),xscm(ixdx),xscm(ixdy),xscm(ixdz),
     1            xscm(ixrxx),xscm(ixryy),xscm(ixrzz),
     2            xscm(ixrxy),xscm(ixrxz),xscm(ixryz),xscm(ixomg),
     3            xscm(ixe),xscm(ixie),xscm(ixvel),xscm(ixvsne),
     4            xscm(ixvsn),xscm(ixvs),xscm(ixvex),xscm(ixr),
     5            xscm(ixvr),xscm(ixind),xscm(ixvr),
     6            xscm(ixij),xscm(ixindx))
c     call setc(loadcm)
      call gmem_free(ixindx)
      call gmem_free(ixij)
      call gmem_free(ixie)
      call gmem_free(ixind)
      call gmem_free(ixr)
      call gmem_free(ixe)
      call gmem_free(ixvex)
      call gmem_free(ixvs)
      call gmem_free(ixvsn)
      call gmem_free(ixvsne)
      call gmem_free(ixvel)
      call gmem_free(ixomg)
      call gmem_free(ixryz)
      call gmem_free(ixrxz)
      call gmem_free(ixrxy)
      call gmem_free(ixrzz)
      call gmem_free(ixryy)
      call gmem_free(ixrxx)
      call gmem_free(ixdz)
      call gmem_free(ixdy)
      call gmem_free(ixdx)
      call gmem_free(ixo)
      call gmem_free(ixvr)
c
c-----  repeat calculations from relay on, if two components
c       of the dielectric response are treated separately
c
      if ((itwoeps .eq. 1) .and. (ieps .eq. 0)) then
        ieps = 1
        goto 100
      endif
c
      ixrelay = 0
      ixomga = 0
c
      if (neqrf .eq. 1) then
        ineq = 1
        ieps = 0
        goto 100
      endif
c
  200 continue
c
      if ((neq2eps .eq. 1) .and. (ieps .eq. 0)) then
        ieps = 1
        goto 100
      endif
c
      ixzfp = 0
      ixwtvr = 0
c
c      call cmem(ngotcm)
c      if(ngotcm.ne.loadcm) call setc(loadcm)
      if (.not. mcupdt) write(iwr,11)
c
      return
      end
      subroutine drfwt(xexp,wt,vr,
     1           xsurf,xnorm,area,ieps,ineq,momsav)
c------
c       calculation of (expanded) source and reaction fields
c       of/on (formal) qm particles (nuclei, electrons)
c       or representations (dp charges, mulliken charges
c       and dipoles) thereof, at polarizabilities,
c       boundary elements and external charges.
c
c       --------  p.th. van duijnen, ibm-kingston 1985, and
c                 groningen, dec. 1991.
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
c
c-----  dummy arrays
c
      dimension xsurf(3,nbem),xnorm(3,nbem),area(nbem)
      dimension xexp(3,nexp),wt(nwtr,nwtc),
     1          vr(nwtr,nwtc)
c
c-----  common blocks
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
c      dimension clas((ngran*(ngran+1)*5)+1)
      dimension clas(3000)
c
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,one,pt5,pt75/0.0d00,1.0d00,0.5d00,0.75d00/
c
c-----  begin
c
      if ((ngran*(ngran+1)*5)+1.gt.3000) call caserr(
     + 'drfctl.f - drfwt: increase dimension of clas >3000')
      ilxs = 81
      ilxn = 82
      ilar = 83
      ilcl = 21
c
      if (ineq .eq. 0) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
          ilzf = 22
          ilzv = 23
          ilwt = 26
          ilvr = 28
c   2-----
        else
c   2-----
          ilzf = 24
          ilzv = 25
          ilwt = 30
          ilvr = 32
          if (itwosur .eq. 1) then
c     3-----
            ilxs = 84
            ilxn = 85
            ilar = 86
c     3-----
          endif
c   2-----
        endif
c 1-----
      else if (momsav .eq. 1) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
          ilvr = 34
        else
          ilvr = 36
          if (itwosur .eq. 1) then
c     3-----
            ilxs = 84
            ilxn = 85
            ilar = 86
c     3-----
          endif
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  the nonequilibrium rf has been saved as an expanded potential,
c         not as a set of induced moments: calculation of vr is not
c         necessary
c
        return
c 1-----
      endif
c
c-----  read expansion centers.
c       = # internal  atoms
c       (+ center of internal nuclear charge)
c
      call daread(idafdrf,iodadrf,xexp,3*nexp,1)
c
c-----  get surface points, normal vectors and areas of
c       boundary elements, if present (in xsurf, xnorm and area -
c       these are calculated in the subroutine drfsurf
c
      if ((field(5:) .ne. ' ') .and. (ibem .ne. 0)) then
c 1-----
        call daread(idafdrf,iodadrf,xsurf,3*nbem,ilxs)
        call daread(idafdrf,iodadrf,xnorm,3*nbem,ilxn)
        call daread(idafdrf,iodadrf,area ,  nbem,ilar)
c 1-----
      endif
c
c-----  calculate zfp: fields in polarizable points due to the charges
c       in other parts of the external system.
c
c       these fields and potentials will give rise to induced
c       charges and dipoles.
c       therefore, -zfp- is made part of -wt-, so that all transformatio
c       with -relay- may be done in one go.
c       -zfp- are the (4*nexp + 1)'th to (4*nexp + ngran)'th columns of
c       (elements 1 to npol3), where ngran is the number of analysis
c       groups for the classical system
c
c       in the same way, the reaction field due to unit induced
c       dipoles and charges at the external charges is calculated
c       and made part of -vr-
c
c       also the classical energy (the interaction energy between
c       the external charges without reaction field contributions
c       is calculated (clas(1))
c
      ngrpair = ngran*(ngran+1)/2
c
      if (ixzfp .ne. 0) then
c 1-----
        call clear(wt(1,nexp4+1),ngran*ndim)
        call clear(vr(1,nexp4+1),ngran*ndim)
        call clear (clas,ngrpair*10+1)
c
        call drfzfp(wt(1,nexp4+1),vr(1,nexp4+1),
     1          clas(1),clas(2),clas(2+ngrpair),clas(2+2*ngrpair),
     2          clas(2+3*ngrpair),clas(2+4*ngrpair),ngrpair)
c
        if ((field(5:) .ne. ' ') .and. (ibem .ne. 0)) then
c   2-----
c    -----  calculate zfbe: fields in boundary elements due to the charg
c           in other parts of the external system.
c
c         these fields and potentials will give rise to induced
c         dipoles (and charges if kappa .ne. 0)
c
c         -zfbe- are the (4*nexp + 1)'th to (4*nexp + ngran)'th columns
c         (elements npol3+1 to npol3+ndimb), where ngran is the number
c         of analysis groups for the classical system
c
c         in the same way, the reaction field due to unit induced
c         dipoles and charges at the external charges is calculated
c         and made part of -vr-
c
          call drfzfbe(ieps,ineq,wt(1,nexp4+1),vr(1,nexp4+1),
     1         xsurf,xnorm,area)
c   2-----
        endif
c
c  -----  save -zfp- and -zfbe- on disk
c
        if (ineq .eq. 0) then
c   2-----
          call dawrit(idafdrf,iodadrf,wt(1,nexp4+1),
     1                ngran*ndim,ilzf,navdrf)
          call dawrit(idafdrf,iodadrf,vr(1,nexp4+1),
     1                ngran*ndim,ilzv,navdrf)
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  read the information from disk and update if required
c
        call daread(idafdrf,iodadrf,clas,ngrpair*10+1,ilcl)
c
        if (mcupdt) then
c   2-----
c    -----  read wt and vr from disk
c
          call daread(idafdrf,iodadrf,wt,nwtr*nwtc,ilwt)
          call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvr)
c
c    -----  recalculate zfp and zfbe + vr according to mc move
c
          do 100, igr = 1, ngran
            call clear(wt(1,nexp4+igr),ndima)
            call clear(vr(1,nexp4+igr),ndima)
  100     continue
c
          call zfpupd(wt(1,nexp4+1),vr(1,nexp4+1),clas(1),clas(2),
     1         clas(2+ngrpair),clas(2+2*ngrpair),clas(2+3*ngrpair),
     2         clas(2+4*ngrpair),ngrpair)
          clasn = clas(1)
c
          if (ibem .ne. 0) then
c     3-----
            call zfbeupd(ieps,ineq,wt(1,nexp4+1),vr(1,nexp4+1),
     1                   xsurf,xnorm,area)
c     3-----
          endif
c   2-----
        else
c   2-----
c    -----  read -zfp- and -zfbe- from disk
c
          call daread(idafdrf,iodadrf,wt(1,nexp4+1),ngran*ndim,ilzf)
          call daread(idafdrf,iodadrf,vr(1,nexp4+1),ngran*ndim,ilzv)
c   2-----
        endif
c 1-----
      endif
c
      if (mcupdt) then
c 1-----
        if (imcout .eq. 5) then
c   2-----
          clase = zero
          do 200, k = 1, ngrpair
            clase = clase + clas(1+k)
  200     continue
c
          write(iwr,9900) clase
 9900     format(
     1    /' - - - updated external coulombic energy   = ',f15.8 )
c   2-----
        endif
c 1-----
      else
c 1-----
        if (idrfout .eq. 3) then
c   2-----
          clase = zero
          do 300, k = 1, ngrpair
            clase = clase + clas(1+k)
  300     continue
c
          write(iwr,9901) clase
 9901     format(
     1/' - - -  external coulombic energy   = ',f15.8 )
c   2-----
        endif
c 1-----
      endif
c
      if (idrfout .eq. 3) then
        call hatout(wt(1,nexp4+1),ndim,ngran,2,'zfp')
        call hatout(vr(1,nexp4+1),ndim,ngran,2,'vrp')
      endif
c
c-----  calculate required parts of source and reaction fields
c       due to qm
c
      if ((field(5:) .ne. '   ') .or. (neqrf .ne. 0)
     1    .or. (iextdip .ne. 0)) then
c 1-----
        if (.not. mcupdt) then
c   2-----
          call wtvrcal(xexp,wt,vr,xsurf,xnorm,area,ieps,ineq)
        else if (notrans .ne. 1) then
          call wtvrupd(xexp,wt,vr)
c   2-----
        endif
c 1-----
      endif
c
      if (idrfout .eq. 3) then
        call hatout(wt,nwtr,nwtc,2,'wt')
        call hatout(vr,nwtr,nwtc,2,'vr')
      endif
c
c-----   save (updated) wt and vr matrices and classical energies
c
      if (mcupdt) then
        ilwt = ilwt + 1
        ilvr = ilvr + 1
        ilcl = 91
      endif
c
      if (ineq .eq. 0)
     1  call dawrit(idafdrf,iodadrf,wt,nwtr*nwtc,ilwt,navdrf)
      call dawrit(idafdrf,iodadrf,vr,nwtr*nwtc,ilvr,navdrf)
      call dawrit(idafdrf,iodadrf,clas,ngrpair*10+1,ilcl,navdrf)
c
      return
      end
      subroutine drfout
c------
c       prints the output of the reaction field contribution
c       analysis
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c      --------  a.h. de vries, groningen 1993 -----
c
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
c
c-----  common blocks
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
c-----  local arrays
c
c      dimension clas((ngran*(ngran+1)*5)+1)
      dimension clas(3000)
c
c
c-----  data statements
c
      data zero, pt25 /0.0d00, 0.25d00/
c
      if ((ngran*(ngran+1)*5)+1.gt.3000) call caserr
     + ('drfctl.f - drfout: increase dimension of clas >3000')
      iw = 6
      ist = iactst
c
      ngrpair = (ngran*(ngran+1))/2
      if (mcupdt .and. (nopes .eq. 0)) then
        call daread(idafdrf,iodadrf,clas,ngrpair*10+1,91)
      else
        call daread(idafdrf,iodadrf,clas,ngrpair*10+1,21)
      endif
c
c-----  print energy analysis if required
c
      if (.not. mcupdt .or. imcout .eq. 5 .or. ipesout .ge. 1) then
c 1-----
        write(iwr,8000) field(:4),field(5:),version
 8000   format(
     1//1x,t20,' analysis of reaction field contributions '
     2/1x ,t20,' static field: ',a4,' reaction field: ',a4,
     3/1x ,t20,' version:  ',a24/)
c
c  -----  repeat header
c
        write (iw,8010) (title(i), i = 1, 10)
 8010   format(10a8)
c
c  -----  options
cxxx  distinction not implemented in this version
c       write (iw,8020)
c8020   format(/,' the source field is brought about by:')
c
c       if (ifldin .eq. 1) write (iw,8030)
c8030   format(' molecule represented by dipole preserving charges')
c
c       if (ifldin .eq. 2) write (iw,8040)
c8040   format(' molecule represented by mulliken charges & dipoles')
c
c       if (ifldin .eq. 3) write (iw,8050)
c8050   format(' molecule after expansion of electronic contributions')
c
c       if (ifldin .eq. 4) write (iw,8060)
c8060   format(' nuclei and electrons separately')
c
c       write (iw,8070)
c8070   format(/,' the reaction field is coupled back to:')
c
c       if (ifldout .eq. 1) write (iw,8030)
c       if (ifldout .eq. 2) write (iw,8040)
c       if (ifldout .eq. 3) write (iw,8050)
c       if (ifldout .eq. 4) write (iw,8060)
c
c       if (idisadd .eq. 0) then
c   2-----
c         write (iw,8080)
c8080     format (/,' the field due to the external charges is ',
c    1              'treated separately ')
c       else
c         write (iw,8090)
c8090     format (/,' the field due to the external charges is ',
c    1              'added to the molecular field ')
c   2-----
c       endif
cxxx
        if (itwoeps .eq. 1) write (iw,8100)
 8100   format(/,' static and optical dielectric constants included')
c
        if (neqrf .eq. 1) write (iw,8110) textneq(ineqix(nneq)-1)
 8110   format(/,' non-equilibrium rf defined by job with header:',
     1         /,a80)
c
        write(iwr,9000)
 9000   format(/,' --- quantum system ---',/)
c
        ic = 0
c
c  -----  if nuclei and electrons are treated separately,
c         print interaction energies
c
        if (ifldout .ge. 3) then
c   2-----
c    -----  print nuclear repulsion energy
c
          ic = ic + 1
          write(iwr,9110) ic, unucrep
 9110     format (' -',
     1 i3,' - vacuum nuclear repulsion energy             =',f20.12)
c
c    -----  one-electron energy (from vacuum hamiltonian) of the qm syst
c           with density obtained from scf calculation with or without
c           static and/or reactionfield
c
          ic = ic + 1
          write(iwr,9120) ic, eoneel(ist)
 9120     format (' -',
     1 i3,' - one-electron energy                         =',f20.12)
c
c    -----  kinetic energy of the qm electrons
c
          ic = ic + 1
          write(iwr,9130) ic, ekin(ist)
 9130     format (' -',
     1 i3,' - electronic kinetic energy                   =',f20.12)
c
c    -----  nuclear attraction energy of qm system
c
          ic = ic + 1
          write(iwr,9140) ic, enua(ist)
 9140     format (' -',
     1 i3,' - nuclear attraction energy                   =',f20.12)
c
c    -----  two-electron energy
c
          ic = ic + 1
          write(iwr,9150) ic, etwoel(ist)
 9150     format (' -',
     1 i3,' - two-electron energy                         =',f20.12)
c   2-----
        endif
c
c  -----  total scf energy of qm system
c
        ic = ic + 1
        write(iwr,9160) ic, uscf(ist)
 9160   format (' -',
     1 i3,' - scf energy                                  =',f20.12)
c
        if (field(5:) .ne. ' ') then
c   2-----
          if (ifldout .eq. 4) then
c     3-----
c      -----  screening of the nuclear repulsion energy,
c             i.e. the interaction between the nuclei and
c             their own reaction field
c             only if the source field of the qm system is
c             treated separately for nuclei and electrons
c
            if (ifldin .eq. 4) then
c       4-----
              ic = ic + 1
              write(iwr,9200) ic, snucnuc(ist)
 9200         format (' -',
     1 i3,' - screening of the nuclear repulsion energy   =',f20.12)
c
c        -----  electronic contributions
c
              if (gamdrf .ne. zero) then
c         5-----
c          -----  interaction of electrons with reaction field
c                 induced by the electron itself
c
                ic = ic + 1
                write(iwr,9210) ic, selel(ist)
 9210           format (' -',
     1 i3,' - screening of the electronic self energy     =',f20.12)
c         5-----
              endif
c
c        -----  modification of electron-nuclear interaction
c               energy by reaction field induced by nuclei
c               and electrons
c
              ic = ic + 1
              write(iwr,9220) ic, snua(ist)
 9220         format (' -',
     1 i3,' - screening of nuclear attraction energy      =',f20.12)
c
c        -----  interaction of electrons with reaction field
c               induced by all electrons
c
              ic = ic + 1
              write(iwr,9230) ic, stwoel(ist)
 9230         format (' -',
     1 i3,' - screening of two-electron energy            =',f20.12)
c       4-----
            else
c       4-----
              if (idisadd .eq. 0) then
c         5-----
c          -----  interaction between nuclei and dipoles
c                 and/or boundary charges (-dipoles)
c                 induced by the molecular field
c
                ic = ic + 1
                write(iwr,9240) ic, smolnuc(ist)
 9240           format (' -',
     1 i3,' - nuclei/ molecular induced reaction field    =',f20.12)
c
c          -----  interaction of electrons with reaction field
c                 induced by molecular field
c
                ic = ic + 1
                write(iwr,9250) ic, smolel(ist)
 9250           format (' -',
     1 i3,' - electrons/ molecular induced reaction field =',f20.12)
c         5-----
              else
c         5-----
                ic = ic + 1
                write(iwr,9260) ic, stotnuc(ist)
 9260           format (' -',
     1 i3,' - nuclei/ total induced reaction field        =',f20.12)
c
                ic = ic + 1
                write(iwr,9270) ic, stotel(ist)
 9270           format (' -',
     1 i3,' - electrons/ total induced reaction field     =',f20.12)
c         5-----
              endif
c       4-----
            endif
c     3-----
          else
c     3-----
            if (ifldin .eq. 4) then
c       4-----
c        -----  separate nulear and electronic source field
c
              ic = ic + 1
              write(iwr,9300) ic, snucmol(ist)
 9300         format (' -',
     1 i3,' - molecule/ nuclear induced reaction field    =',f20.12)
c
              ic = ic + 1
              write(iwr,9310) ic, selmol(ist)
 9310         format (' -',
     1 i3,' - molecule/ electronic induced reaction field =',f20.12)
c       4-----
            else
c       4-----
              if (idisadd .eq. 0) then
c         5-----
                ic = ic + 1
                write(iwr,9320) ic, smolmol(ist)
 9320           format (' -',
     1 i3,' - molecule/ molecular induced reaction field  =',f20.12)
c         5-----
              else
c         5-----
                ic = ic + 1
                write(iwr,9330) ic, stotmol(ist)
 9330           format (' -',
     1 i3,' - molecule/ total induced reaction field      =',f20.12)
c         5-----
              endif
c       4-----
            endif
c     3-----
          endif
c
          if (idisadd .eq. 0) then
c     3-----
            ic = ic + 1
            write(iwr,9340) ic, suqm(ist)
 9340       format (' -',
     1 i3,' - molecular reaction field stabilisation      =',f20.12)
c
            ic = ic + 1
            write(iwr,9350) ic, upolqm(ist)
 9350       format (' -',
     1 i3,' - equilibrium molecular polarisation energy   =',f20.12)
c     3-----
          endif
c
          if (itwoeps .eq. 1) then
c     3-----
c      -----  separate optic and static contributions to rf
c
            write(iwr,9390)
 9390       format(/,' --- optic dielectric contributions ---')
c
            if (ifldout .eq. 4) then
c       4-----
c        -----  screening of the nuclear repulsion energy,
c               i.e. the interaction between the nuclei and
c               their own reaction field
c               only if the source field of the qm system is
c               treated separately for nuclei and electrons
c
              if (ifldin .eq. 4) then
c         5-----
                ic = ic + 1
                write(iwr,9200) ic, snucno(ist)
c
c          -----  electronic contributions
c
                if (gamdrf .ne. zero) then
c           6-----
c            -----  interaction of electrons with reaction field
c                   induced by the electron itself
c
                  ic = ic + 1
                  write(iwr,9210) ic, selel(ist)
c           6-----
                endif
c
c          -----  modification of electron-nuclear interaction
c                 energy by reaction field induced by nuclei
c                 and electrons
c
                ic = ic + 1
                write(iwr,9220) ic, snuao(ist)
c
c          -----  interaction of electrons with reaction field
c                 induced by all electrons
c
                ic = ic + 1
                write(iwr,9230) ic, stwoelo(ist)
c         5-----
              else
c         5-----
                if (idisadd .eq. 0) then
c           6-----
c            -----  interaction between nuclei and dipoles
c                   and/or boundary charges (-dipoles)
c                   induced by the molecular field
c
                  ic = ic + 1
                  write(iwr,9240) ic, smolno(ist)
c
c            -----  interaction of electrons with reaction field
c                   induced by molecular field
c
                  ic = ic + 1
                  write(iwr,9250) ic, smolelo(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9260) ic, stotno(ist)
c
                  ic = ic + 1
                  write(iwr,9270) ic, stotelo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            else
c       4-----
              if (ifldin .eq. 4) then
c         5-----
c          -----  separate nulear and electronic source field
c
                ic = ic + 1
                write(iwr,9300) ic, snucmo(ist)
c
                ic = ic + 1
                write(iwr,9310) ic, selmolo(ist)
c         5-----
              else
c         5-----
                if (idisadd .eq. 0) then
c           6-----
                  ic = ic + 1
                  write(iwr,9320) ic, smolmo(ist)
                else
                  ic = ic + 1
                  write(iwr,9330) ic, stotmo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            endif
c
            if (idisadd .eq. 0) then
c       4-----
              ic = ic + 1
              write(iwr,9340) ic, suqmo(ist)
c
              ic = ic + 1
              write(iwr,9350) ic, upolqmo(ist)
c       4-----
            endif
c
            write(iwr,9395)
 9395       format(/,' --- static dielectric contributions ---')
c
            if (ifldout .eq. 4) then
c       4-----
c        -----  screening of the nuclear repulsion energy,
c               i.e. the interaction between the nuclei and
c               their own reaction field
c               only if the source field of the qm system is
c               treated separately for nuclei and electrons
c
              if (ifldin .eq. 4) then
c         5-----
                ic = ic + 1
                write(iwr,9200) ic, snucnuc(ist) - snucno(ist)
c
c          -----  electronic contributions
c
c          -----  modification of electron-nuclear interaction
c                 energy by reaction field induced by nuclei
c                 and electrons
c
                ic = ic + 1
                write(iwr,9220) ic, snua(ist) - snuao(ist)
c
c          -----  interaction of electrons with reaction field
c                 induced by all electrons
c
                ic = ic + 1
                write(iwr,9230) ic, stwoel(ist) - stwoelo(ist)
c         5-----
              else
c         5-----
                if (idisadd .eq. 0) then
c           6-----
c            -----  interaction between nuclei and dipoles
c                   and/or boundary charges (-dipoles)
c                   induced by the molecular field
c
                  ic = ic + 1
                  write(iwr,9240) ic, smolnuc(ist) - smolno(ist)
c
c            -----  interaction of electrons with reaction field
c                   induced by molecular field
c
                  ic = ic + 1
                  write(iwr,9250) ic, smolel(ist) - smolelo(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9260) ic, stotnuc(ist) -  stotno(ist)
c
                  ic = ic + 1
                  write(iwr,9270) ic, stotel(ist) - stotelo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            else
c       4-----
              if (ifldin .eq. 4) then
c         5-----
c          -----  separate nulear and electronic source field
c
                ic = ic + 1
                write(iwr,9300) ic, snucmol(ist) - snucmo(ist)
c
                ic = ic + 1
                write(iwr,9310) ic, selmol(ist) - selmolo(ist)
c         5-----
              else
c         5-----
                if (idisadd .eq. 0) then
c           6-----
                  ic = ic + 1
                  write(iwr,9320) ic, smolmol(ist) - smolmo(ist)
                else
                  ic = ic + 1
                  write(iwr,9330) ic, stotmol(ist) - stotmo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            endif
c
            if (idisadd .eq. 0) then
c       4-----
              ic = ic + 1
              write(iwr,9340) ic, suqm(ist) - suqmo(ist)
c
              ic = ic + 1
              write(iwr,9350) ic, upolqm(ist) - upolqmo(ist)
c       4-----
            endif
c     3-----
          endif
c   2-----
        endif
c
        if (neqsta .eq. 1) then
c   2-----
c    -----  interaction of qm system with non-equilibrium external stati
c           field
c
          if (ifldout .eq. 4) then
c     3-----
            ic = ic + 1
            write(iwr,9360) ic, ustanuc(ist)
 9360       format (' -',
     1 i3,' - nonequilibrium static field/nuclear energy  =',f20.12)
c
            ic = ic + 1
            write(iwr,9370) ic, ustael(ist)
 9370       format (' -',
     1 i3,' - nonequilibrium static field/electronic int  =',f20.12)
c     3-----
          endif
c
          ic = ic + 1
          write(iwr,9380) ic, ustaqm(ist)
 9380     format (' -',
     1 i3,' - nonequilibrium static field/molecular int.n =',f20.12)
c   2-----
        endif
c
        if (neqrf .eq. 1) then
c   2-----
c    -----  interaction of qm system with non-equilibrium rf
c
          if (ifldout .eq. 4) then
c     3-----
            ic = ic + 1
            write(iwr,9365) ic, uneqnuc(ist)
 9365       format (' -',
     1 i3,' - non-equilibrium rf / nuclear int.n energy   =',f20.12)
c
            ic = ic + 1
            write(iwr,9375) ic, uneqel(ist)
 9375       format (' -',
     1 i3,' - non-equilibrium rf /electronic int.n energy =',f20.12)
c     3-----
          endif
c
          ic = ic + 1
          write(iwr,9385) ic, uneqqm(ist)
 9385     format (' -',
     1 i3,' - non-equilibrium rf / molecular int.n energy =',f20.12)
c   2-----
        endif
c
c  -----  classical (interaction) energy analysis
c
        if (nodiscr .eq. 0) then
c   2-----
          do 100, igr = 1, ngran
c     3-----
            do 200, jgr = 1, igr
c       4-----
              indxcl = ia(igr) + jgr
c
              write(iwr,9400) igr , jgr
 9400         format(//' --- classical system numbers: ',i2,
     1        ' and ',i2, ' ---',/)
c
              if (iclinte .eq. 1) then
c         5-----
                ic = ic + 1
                write(iwr,9410) ic, clas(1+indxcl)
 9410           format (' -',
     1 i3,' - vacuum classical electrostatic interaction  =',f20.12)
c         5-----
              endif
c
              if (iclintd .eq. 1) then
c         5-----
                ic = ic + 1
                write(iwr,9420) ic, clas(1+ngrpair+indxcl)
 9420           format (' -',
     1 i3,' - classical dispersion energy estimate        =',f20.12)
c         5-----
              endif
c
              if (iclintr .eq. 1) then
c         5-----
                ic = ic + 1
                write(iwr,9430) ic, clas(1+2*ngrpair+indxcl)
 9430           format (' -',
     1 i3,' - classical repulsion energy estimate         =',f20.12)
c         5-----
              endif
c
              ic = ic + 1
              write(iwr,9450) ic, clas(1+3*ngrpair+indxcl)
 9450         format (' -',
     1 i3,' - vacuum classical interaction energy         =',f20.12)
c
              if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0)) then
c         5-----
                if (idisadd .eq. 0) then
c           6-----
                  ic = ic + 1
                  write(iwr,9460) ic, suclasg(igr,jgr)
 9460             format (' -',
     1 i3,' - screening of classical electrostatic energy =',f20.12)
c
                  if (igr .eq. jgr) then
c             7-----
                    ic = ic + 1
                    write(iwr,9470) ic, upolclg(igr)
 9470               format (' -',
     1 i3,' - equilibrium classical polarisation energy   =',f20.12)
c             7-----
                  endif
c
                  if (itwoeps .eq. 1) then
c             7-----
                    write(iwr,9390)
c
                    ic = ic + 1
                    write(iwr,9460) ic, suclsog(igr,jgr)
c
                    if (igr .eq. jgr) then
                      ic = ic + 1
                      write(iwr,9470) ic, uplclog(igr)
                    endif
c
                    write(iwr,9395)
c
                    ic = ic + 1
                    write(iwr,9460) ic,
     1                   suclasg(igr,jgr) - suclsog(igr,jgr)
c
                    if (igr .eq. jgr) then
                      ic = ic + 1
                      write(iwr,9470) ic,
     1                     upolclg(igr) - uplclog(igr)
                    endif
c             7-----
                  endif
c           6-----
                endif
c         5-----
              endif
c       4-----  next jgr
  200       continue
c
            if (neqsta .eq. 1) then
c       4-----
              ic = ic + 1
              write(iwr,9480) ic, ustaclg(igr)
 9480       format (' -',
     1 i3,' - nonequilibrium static field/classical int.n =',f20.12)
c       4-----
            endif
c
            if (neqrf .eq. 1) then
c       4-----
              ic = ic + 1
              write(iwr,9485) ic, uneqclg(igr)
 9485       format (' -',
     1 i3,' - non-equilibrium rf / classical int.n energy =',f20.12)
c       4-----
            endif
c     3-----  next igr
  100     continue
c
c    -----  collect contributions
c
          write(iwr,9490)
 9490     format(/,
     1   ' --- total interaction energies classical system ---')
c
          if (iclinte .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9410) ic, uclase
c     3-----
          endif
c
          if (iclintd .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9420) ic, uclasd
c     3-----
          endif
c
          if (iclintr .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9430) ic, uclasr
c     3-----
          endif
c
          ic = ic + 1
          write(iwr,9495) ic, uclas
 9495     format (' -',
     1 i3,' - total vacuum classical interaction energy   =',f20.12)
c
          if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0)) then
c     3-----
            if (idisadd .eq. 0) then
c       4-----
              ic = ic + 1
              write(iwr,9460) ic, suclas
c
              ic = ic + 1
              write(iwr,9470) ic, upolcl
c
              if (itwoeps .eq. 1) then
c         5-----
                write(iwr,9390)
c
                ic = ic + 1
                write(iwr,9460) ic, suclaso
c
                ic = ic + 1
                write(iwr,9470) ic, upolclo
c
                write(iwr,9395)
c
                ic = ic + 1
                write(iwr,9460) ic, suclas - suclaso
c
                ic = ic + 1
                write(iwr,9470) ic, upolcl - upolclo
c         5-----
              endif
c       4-----
            endif
c     3-----
          endif
c
          if (neqsta .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9480) ic, ustacl
c     3-----
          endif
c
          if (neqrf .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9485) ic, uneqcl
c     3-----
          endif
c   2-----
        endif
c
        write(iwr,9500)
 9500   format(//,' --- interaction quantum mechanical - classical',
     1             ' system ---',/)
c
        if (nodiscr .eq. 0) then
c   2-----
          do 300, igr = 1, ngran
c     3-----
            write(iwr,9505) igr
 9505       format(/,' --- classical group number: ',i2,' ---',/)
c
            if (ifldout .eq. 4) then
c       4-----
              if (field(:4) .ne. ' ') then
                ic = ic + 1
                write(iwr,9510) ic, extnucg(igr,ist)
 9510           format (' -',
     1 i3,' - nuclei/ external charge interaction         =',f20.12)
c
c      -----  interaction energy between electrons and external charges
c
                ic = ic + 1
                write(iwr,9520) ic, extelg(igr,ist)
 9520           format (' -',
     1 i3,' - electrons/ external charge interaction      =',f20.12)
c
              endif
c
              if (field(5:) .ne. ' ') then
c         5-----
                if (idisadd .eq. 0) then
c           6-----
c            -----  screening of the interaction between nuclei
c                   and external charges
c
                  if (iextdip .eq. 0) then
c             7-----
                    ic = ic + 1
                    write(iwr,9530) ic, sxtnucg(igr,ist)
 9530               format (' -',
     1 i3,' - screening of the ext. charge/nuclear int.n  =',f20.12)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9540) ic, sextelg(igr,ist)
 9540               format (' -',
     1 i3,' - screening of the ext. charge/electron int.n =',f20.12)
c             7-----
                  endif
c
                  if (ifldin .eq. 4) then
c             7-----
c              -----  screening of the interaction between nuclei
c                     and external charges
c
                    ic = ic + 1
                    write(iwr,9550) ic, snucxtg(igr,ist)
 9550               format (' -',
     1 i3,' - screening of the nuclei/ ext. charge int.n  =',f20.12)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9560) ic, selextg(igr,ist)
 9560               format (' -',
     1 i3,' - screening of the electron/ext. charge int.n =',f20.12)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolxtg(igr,ist)
 9570               format (' -',
     1 i3,' - screening of the molecule/ext. charge int.n =',f20.12)
c             7-----
                  endif
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotxtg(igr,ist)
 9580             format (' -',
     1 i3,' - screening total system/ext. charge int.n    =',f20.12)
c           6-----
                endif
c
                if (itwoeps .eq. 1) then
c           6-----
                  write(iwr,9390)
c
                  if (idisadd .eq. 0) then
c             7-----
c              -----  screening of the interaction between nuclei
c                     and external charges
c
                    if (iextdip .eq. 0) then
c               8-----
                      ic = ic + 1
                      write(iwr,9530) ic, sextnog(igr,ist)
c
c                -----  screening of the interaction energy between
c                       electrons and external charges
c
                      ic = ic + 1
                      write(iwr,9540) ic, sxtelog(igr,ist)
c               8-----
                    endif
c
                    if (ifldin .eq. 4) then
c               8-----
c                -----  screening of the interaction between nuclei
c                       and external charges
c
                      ic = ic + 1
                      write(iwr,9550) ic, sncexog(igr,ist)
c
c                -----  screening of the interaction energy between
c                       electrons and external charges
c
                      ic = ic + 1
                      write(iwr,9560) ic, selxtog(igr,ist)
c               8-----
                    else
c               8-----
                      ic = ic + 1
                      write(iwr,9570) ic, smlexog(igr,ist)
c               8-----
                    endif
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9580) ic, sttexog(igr,ist)
c             7-----
                  endif
c
                  write(iwr,9395)
c
                  if (idisadd .eq. 0) then
c             7-----
c              -----  screening of the interaction between nuclei
c                     and external charges
c
                    if (iextdip .eq. 0) then
c               8-----
                      ic = ic + 1
                      write(iwr,9530) ic,
     1                sxtnucg(igr,ist) - sextnog(igr,ist)
c
c                -----  screening of the interaction energy between
c                       electrons and external charges
c
                      ic = ic + 1
                      write(iwr,9540) ic,
     1                sextelg(igr,ist) - sxtelog(igr,ist)
c               8-----
                    endif
c
                    if (ifldin .eq. 4) then
c               8-----
c                -----  screening of the interaction between nuclei
c                       and external charges
c
                      ic = ic + 1
                      write(iwr,9550) ic,
     1                snucxtg(igr,ist) - sncexog(igr,ist)
c
c                -----  screening of the interaction energy between
c                       electrons and external charges
c
                      ic = ic + 1
                      write(iwr,9560) ic,
     1                selextg(igr,ist) - selxtog(igr,ist)
c               8-----
                    else
c               8-----
                      ic = ic + 1
                      write(iwr,9570) ic,
     1                smolxtg(igr,ist) - smlexog(igr,ist)
c               8-----
                    endif
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9580) ic,
     1              stotxtg(igr,ist) - sttexog(igr,ist)
c             7-----
                  endif
c           6-----
                endif
c         5-----
              endif
c       4-----
            else
c       4-----
              if (field(5:) .ne. ' ') then
c         5-----
c          -----  screening of the interaction between nuclei
c                 and external charges
c
                if (idisadd .eq. 0) then
c           6-----
                  if (ifldin .eq. 4) then
c             7-----
                    ic = ic + 1
                    write(iwr,9550) ic, snucxtg(igr,ist)
c
                    ic = ic + 1
                    write(iwr,9560) ic, selextg(igr,ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolxtg(igr,ist)
c             7-----
                  endif
c
                  ic = ic + 1
                  write(iwr,9590) ic, sxtmolg(igr,ist)
 9590             format (' -',
     1 i3,' - screening of the ext. charge/molecule int.n =',f20.12)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotxtg(igr,ist)
c           6-----
                endif
c
                if (itwoeps .eq. 1) then
c           6-----
                  write(iwr,9390)
c
                  if (idisadd .eq. 0) then
c             7-----
                    if (ifldin .eq. 4) then
c               8-----
                      ic = ic + 1
                      write(iwr,9550) ic, sncexog(igr,ist)
c
                      ic = ic + 1
                      write(iwr,9560) ic, selxtog(igr,ist)
c               8-----
                    else
c               8-----
                      ic = ic + 1
                      write(iwr,9570) ic, smlexog(igr,ist)
c               8-----
                    endif
c
                    ic = ic + 1
                    write(iwr,9590) ic, sextmog(igr,ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9580) ic, sttexog(igr,ist)
c             7-----
                  endif
c
                  write(iwr,9395)
c
                  if (idisadd .eq. 0) then
c             7-----
                    if (ifldin .eq. 4) then
c               8-----
                      ic = ic + 1
                      write(iwr,9550) ic,
     1                snucxtg(igr,ist) - sncexog(igr,ist)
c
                      ic = ic + 1
                      write(iwr,9560) ic,
     1                selextg(igr,ist) - selxtog(igr,ist)
c               8-----
                    else
c               8-----
                      ic = ic + 1
                      write(iwr,9570) ic,
     1                smolxtg(igr,ist) - smlexog(igr,ist)
c               8-----
                    endif
c
                    ic = ic + 1
                    write(iwr,9590) ic,
     1              sxtmolg(igr,ist) - sextmog(igr,ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9580) ic,
     1              stotxtg(igr,ist) - sttexog(igr,ist)
c             7-----
                  endif
c           6-----
                endif
c         5-----
              endif
c       4-----
            endif
c
c      -----  molecular contributions
c
            if (field(:4) .ne. ' ') then
              ic = ic + 1
              write(iwr,9600) ic, uelstg(igr,ist)
 9600       format (/,' -',
     1 i3,' - electrostatic qm-classical interaction      =',f20.12)
c
            endif
c
            if ((field(5:) .ne. ' ')
     1       .and. (idisadd .eq. 0)) then
c       4-----
              ic = ic + 1
              write(iwr,9610) ic, suintg(igr,ist)
 9610         format (' -',
     1 i3,' - screening electrostatic qm-classical int.n  =',f20.12)
c
              if (itwoeps .eq. 1) then
c         5-----
                write(iwr,9390)
c
                ic = ic + 1
                write(iwr,9610) ic, suintog(igr,ist)
c
                write(iwr,9395)
c
                ic = ic + 1
                write(iwr,9610) ic, suintg(igr,ist) - suintog(igr,ist)
c         5-----
              endif
c       4-----
            endif
c  
            if (field(:4) .ne. ' ') then
              ic = ic + 1
              write(iwr,9720) ic, repmodg(igr,ist)
 9720       format (' -',
     1 i3,' - model repulsion energy                      =',f20.12)
            endif
c     3-----  next classical group
  300     continue
c
c    -----  collect contributions
c
          write(iwr,9725)
 9725     format(/,' --- total quantum mechanical ',
     1'- classical interaction ---',/)
c
          if (ifldout .eq. 4) then
c     3-----
            if (field(:4) .ne. ' ') then
c       4-----
              ic = ic + 1
              write(iwr,9510) ic, extnuc(ist)
c
c      -----  interaction energy between electrons and external charges
c
              ic = ic + 1
              write(iwr,9520) ic, extel(ist)
c       4-----
            endif
c
            if (field(5:) .ne. ' ') then
c       4-----
              if (idisadd .eq. 0) then
c         5-----
c          -----  screening of the interaction between nuclei
c                 and external charges
c
                if (iextdip .eq. 0) then
c           6-----
                  ic = ic + 1
                  write(iwr,9530) ic, sextnuc(ist)
c
c            -----  screening of the interaction energy between
c                   electrons and external charges
c
                  ic = ic + 1
                  write(iwr,9540) ic, sextel(ist)
c           6-----
                endif
c
                if (ifldin .eq. 4) then
c           6-----
c            -----  screening of the interaction between nuclei
c                   and external charges
c
                  ic = ic + 1
                  write(iwr,9550) ic, snucext(ist)
c
c            -----  screening of the interaction energy between
c                   electrons and external charges
c
                  ic = ic + 1
                  write(iwr,9560) ic, selext(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9570) ic, smolext(ist)
c           6-----
                endif
c         5-----
              else
c         5-----
                ic = ic + 1
                write(iwr,9580) ic, stotext(ist)
c         5-----
              endif
c
              if (itwoeps .eq. 1) then
c         5-----
                write(iwr,9390)
c
                if (idisadd .eq. 0) then
c           6-----
c            -----  screening of the interaction between nuclei
c                   and external charges
c
                  if (iextdip .eq. 0) then
c             7-----
                    ic = ic + 1
                    write(iwr,9530) ic, sextno(ist)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9540) ic, sextelo(ist)
c             7-----
                  endif
c
                  if (ifldin .eq. 4) then
c             7-----
c              -----  screening of the interaction between nuclei
c                     and external charges
c
                    ic = ic + 1
                    write(iwr,9550) ic, snucexo(ist)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9560) ic, selexto(ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolexo(ist)
c             7-----
                  endif
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotexo(ist)
c           6-----
                endif
c
                write(iwr,9395)
c
                if (idisadd .eq. 0) then
c           6-----
c            -----  screening of the interaction between nuclei
c                   and external charges
c
                  if (iextdip .eq. 0) then
c             7-----
                    ic = ic + 1
                    write(iwr,9530) ic, sextnuc(ist) - sextno(ist)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9540) ic, sextel(ist) - sextelo(ist)
                  endif
c             7-----
                  if (ifldin .eq. 4) then
c             7-----
c              -----  screening of the interaction between nuclei
c                     and external charges
c
                    ic = ic + 1
                    write(iwr,9550) ic, snucext(ist) - snucexo(ist)
c
c              -----  screening of the interaction energy between
c                     electrons and external charges
c
                    ic = ic + 1
                    write(iwr,9560) ic, selext(ist) - selexto(ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolext(ist) - smolexo(ist)
c             7-----
                  endif
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotext(ist) - stotexo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            endif
c     3-----
          else
c     3-----
            if (field(5:) .ne. ' ') then
c       4-----
c        -----  screening of the interaction between nuclei
c               and external charges
c
              if (idisadd .eq. 0) then
c         5-----
                if (ifldin .eq. 4) then
c           6-----
                  ic = ic + 1
                  write(iwr,9550) ic, snucext(ist)
c
                  ic = ic + 1
                  write(iwr,9560) ic, selext(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9570) ic, smolext(ist)
c           6-----
                endif
c
                ic = ic + 1
                write(iwr,9590) ic, sextmol(ist)
c         5-----
              else
c         5-----
                ic = ic + 1
                write(iwr,9580) ic, stotext(ist)
c         5-----
              endif
c
              if (itwoeps .eq. 1) then
c         5-----
c
                write(iwr,9390)
c
                if (idisadd .eq. 0) then
c           6-----
                  if (ifldin .eq. 4) then
c             7-----
                    ic = ic + 1
                    write(iwr,9550) ic, snucexo(ist)
c
                    ic = ic + 1
                    write(iwr,9560) ic, selexto(ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolexo(ist)
c             7-----
                  endif
c
                  ic = ic + 1
                  write(iwr,9590) ic, sextmo(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotexo(ist)
c           6-----
                endif
c
                write(iwr,9395)
c
                if (idisadd .eq. 0) then
c           6-----
                  if (ifldin .eq. 4) then
c             7-----
                    ic = ic + 1
                    write(iwr,9550) ic, snucext(ist) - snucexo(ist)
c
                    ic = ic + 1
                    write(iwr,9560) ic, selext(ist) - selexto(ist)
c             7-----
                  else
c             7-----
                    ic = ic + 1
                    write(iwr,9570) ic, smolext(ist) - smolexo(ist)
c             7-----
                  endif
c
                  ic = ic + 1
                  write(iwr,9590) ic, sextmol(ist) - sextmo(ist)
c           6-----
                else
c           6-----
                  ic = ic + 1
                  write(iwr,9580) ic, stotext(ist) - stotexo(ist)
c           6-----
                endif
c         5-----
              endif
c       4-----
            endif
c     3-----
          endif
c
c  -----  molecular contributions
c
          if (field(:4) .ne. ' ') then
            ic = ic + 1
            write(iwr,9600) ic, uelst(ist)
          endif
c
          if ((field(5:) .ne. ' ')
     1       .and. (idisadd .eq. 0)) then
c     3-----
            ic = ic + 1
            write(iwr,9610) ic, suint(ist)
c
            if (itwoeps .eq. 1) then
c       4-----
              write(iwr,9390)
c
              ic = ic + 1
              write(iwr,9610) ic, suinto(ist)
c
              write(iwr,9395)
c
              ic = ic + 1
              write(iwr,9610) ic, suint(ist) - suinto(ist)
c       4-----
            endif
c     3-----
          endif
c   2-----
        endif
c
c  -----  other interaction energies
c
        if (gamdrf .ne. zero) then
c   2-----
c    -----  estimate of the dispersion energy
c
          ic = ic + 1
          write(iwr,9700) ic, udisp(ist)
 9700     format (/,' -',
     1 i3,' - estimate of the dispersion energy           =',f20.12)
c
          if (irevdis .eq. 1) then
c     3-----
c      -----  estimate of the "reverse" dispersion energy
c
            ic = ic + 1
            write(iwr,9710) ic, rdisp(ist)
 9710       format (' -',
     1 i3,' - estimate of the "reverse" dispersion energy =',f20.12)
c     3-----
          endif
c   2-----
        endif
c
c  -----  estimated short-range qm-classical interaction
c
        if (field(:4) .ne. ' ') then
          ic = ic + 1
          write(iwr,9720) ic, repmod(ist)
        endif
c
c  -----  summary of energy contributions
c
        write(iwr,9800)
 9800   format(//' --- summary of energy contributions ---',/)
c
c  -----  energy of quantum mechanically treated system
c
        ic = 1
        write(iwr,9810) ic, uqm(ist)
 9810   format (' -',
     1 i3,' - total energy of quantum system              =',f20.12)
c
c  -----  energy of classically treated system
c
        if (nodiscr .eq. 0) then
c   2-----
          ic = ic + 1
          write(iwr,9820) ic, uclas + suclas
 9820     format (' -',
     1 i3,' - total energy of classical system            =',f20.12)
c   2-----
        endif
c
c  -----  interaction energy
c
        totint = uint(ist) + udisp(ist) + rdisp(ist)
c
        if (field(5:) .ne. ' ') then
c   2-----
          if (idisadd .eq. 0) then
            totint = totint + smolmol(ist) + suint(ist)
          else
            totint = totint + stabtot(ist)
          endif
c   2-----
        endif
c
        ic = ic + 1
        write(iwr,9830) ic, totint
 9830   format (' -',
     1 i3,' - quant. mech. /classical interaction energy  =',f20.12)
c
c  -----  polarization energy
c
        if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0)) then
c   2-----
          if (itwoeps .eq. 1) then
c     3-----
            ic = ic + 1
            write(iwr,9840) ic, upoleqo(ist)
 9840       format (' -',
     1 i3,' - equilibrium polarisation energy (optic)     =',f20.12)
c
            ic = ic + 1
            write(iwr,9850) ic, upoleq(ist) - upoleqo(ist)
 9850       format (' -',
     1 i3,' - equilibrium polarisation energy (static)    =',f20.12)
c     3-----
          endif
c
          ic = ic + 1
          write(iwr,9860) ic, upoleq(ist)
 9860     format (' -',
     1 i3,' - equilibrium polarisation energy             =',f20.12)
c
          if (isolsav .eq. 1) then
c     3-----
cnotimp     ic = ic + 1
c           write(iwr,9870) ic, ucstst(ist)
c9870       format (' -',
c    1 i3,' - equilibrium external static energy cost     =',f20.12)
cnotimp
cnotinf     ic = ic + 1
c           write(iwr,9875) ic, ucstpl(ist)
c9875       format (' -',
c    1 i3,' - equilibrium rf energy cost                  =',f20.12)
c     3-----
          endif
c   2-----
        endif
c
c  -----  nonequilibrium static field contributions
c
        if (neqsta .eq. 1) then
c   2-----
          ic = ic + 1
          write(iwr,9880) ic, usta(ist)
 9880     format (' -',
     1 i3,' - total int.n with nonequilibrium static fld  =',f20.12)
c
          ic = ic + 1
          write(iwr,9890) ic, ustaneq(ist)
 9890     format (' -',
     1 i3,' - nonequilibrium static field energy cost     =',f20.12)
c   2-----
        endif
c
c  -----  non-equilibrium rf contributions
c
        if (neqrf .eq. 1) then
c   2-----
          ic = ic + 1
          write(iwr,9885) ic, uneq(ist)
 9885     format (' -',
     1 i3,' - total interaction with non-equilibrium rf   =',f20.12)
c
          ic = ic + 1
          write(iwr,9895) ic, upolneq(ist)
 9895     format (' -',
     1 i3,' - non-equilibrium polarisation energy cost    =',f20.12)
c   2-----
        endif
c
c  -----  ensemble total energy
c
        ic = ic + 1
        write(iwr,9910) ic, uens(ist)
 9910   format (' -',
     1 i3,' - configuration total energy                  =',f20.12)
c 1-----
      endif
c
      write(iwr,9920)
 9920 format(/104('-')/)
c
      return
      end
      subroutine drfrelay(ieps,xscm)
      implicit real*8  (a-h,o-z),integer  (i-n)
c------
c      p.th. van duijnen, groningen, dec.1991.
c      a.h. de vries, groningen, july 1993.
c
c      form the relay matrix
c
c      the relay matrix consists of maximally 9 blocks, depending
c      on the type of environment(s) present. its general structure is:
c
c       (  a(-1) + t(q;p))  | del(p) k(i;p) s(i) |  del(p) l(i;p) s(i)
c                           |                    |
c       --------------------|--------------------|--------------------
c                           |                    |
c         -ef* f(q;i)       | 1 - ef*k(j;i) s(j) |  - ef*l(j;i) s(j)
c                           |                    |
c       --------------------|--------------------|--------------------
c                           |                    |
c -eps*ef*del(i) f(q;i).n(i)|   -ef*m(j;i) s(j)  |  1 - ef*n(j;i) s(j)
c                           |                    |
c
c       in which:-a(-1) + t(q;p) couple polarizabilities in q and p,
c                -del(p) k(i;p) is the field in polarizable point p
c                 due to a surface dipole layer element in i
c                 (only with dielectric)
c                -del(p) l(i;p) is the field in polarizable point p
c                 due to a surface charge layer element in i (only with
c                 finite ionic strength)
c                -k(j;i), l(j;i), m(j;i) and n(j;i) couple the surface
c                 dipoles and charges at j and i
c                 (l,m and n only with finite ionic strength)
c                -f(q;i) is the potential of a induced dipole in q at
c                 surface element i (only with dielectric)
c                -del(i) f(q;i).n(i) is the field of an induced dipole
c                 in q at surface element i, contracted with the surface
c                 normal vector at i (only with finite ionic strength)
c
c                -ef: factor: 1/(2pi*(1+eps)
c                -s(i): area of surface element i (constant element appr
c
c------
c
c-----  common blocks
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
c
      integer maxblnk
      common /freein/ maxblnk
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
      dimension xscm(*)
c
c----  local arrays
c
      character *8 errmsg(3)
      logical group
      character*80 card,nami*16
      dimension axyz(6)
      data errmsg/'program','stop in','-drfrela'/
      data zero, three /0.0d0, 3.0d00/
c
c-----  begin
c
c-----  set some dimensions
c
c     a is (npol3,npol3)
c
c     b is (nbem,nbem) or (if kappa.ne.0.0) (2*nbem,2*nbem)
c
c     f is (npol3,nbem) or (if kappa.ne.0.0) (npol3,2*nbem)
c
c     d is (nbem,npol3) or (if kappa.ne.0.0) (2*nbem,npol3)
c
      group = ngrpol .ne. 0
c
      if ( .not. mcupdt) then
        if (idrfout .ge. 1) then
          if (npol .ne. 0) write(iwr,9910) ndima
          if (nbem .ne. 0) write(iwr,9911) ndimb
          write(iwr,9912) ndim
 9910   format(
     1  /' dimension of polarizability problem, ndima= ',i6)
 9911   format(
     1  /' dimension of boundary problem,       ndimb= ',i6)
 9912   format(
     1  /' dimension of total linear  problem,   ndim= ',i6/)
        endif
      endif
c
c-----  calculate addresses for -relay-, and if dielectric:
c       -xsurf-, -xnorm-, and -area-
c
c-----  -relay- at ixr
c
c     call setscm(ixr)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(ixr))
      ixr = igmem_alloc(ndim*ndim)
      call clear(xscm(ixr),ndim*ndim)
c
c-----  coordinates of surface elements  at ixs
c-----  normal vectors at ixn
c-----  areas at ixa
c
      ixs = igmem_alloc(3*nbem)
      ixn = igmem_alloc(3*nbem)
      ixa = igmem_alloc(nbem)
c     ixs = ixr + ndim*ndim
c     ixn = ixs + 3*nbem
c     ixa = ixn + 3*nbem
c
c     last = max(ixa+nbem+1,ixs+ndimb*ndimb)
c     need = loc10 + last - ixr
c     call setc(need)
c     call clear(xscm(ixr),last-ixr+1)
c
      if (ieps .eq. 0) then
        ilr = 41
        ilru = 42
        illur = 43
        ilindx = 44
        ilxs = 81
        ilxn = 82
        ilar = 83
c
        if ((.not. mcupdt) .or. (imcout .eq. 5)) then
          if ((ibem .ne. 0) .and. (idrfout .ge. 1))
     1      write (iwr,9920)
 9920         format(/,'  -- total dielectric constant --')
        endif
      else
        ilr = 45
        ilru = 46
        illur = 47
        ilindx = 48
        ilxs = 81
	ilxn = 82
	ilar = 83
	if (itwosur .eq. 1) then
          ilxs = 84
          ilxn = 85
          ilar = 86
        endif
c
        if ((.not. mcupdt) .or. (imcout .eq. 5)) then
          if ((ibem .ne. 0) .and. (idrfout .ge. 1)) 
     1      write (iwr,9921)
 9921     format(/,'  -- optical dielectric constant --')
        endif
      endif
c
      if (ibem .ne. 0) then
c
c  -----  get surface points, normal vectors and areas
c
        call daread(idafdrf,iodadrf,xscm(ixs),3*nbem,ilxs)
        call daread(idafdrf,iodadrf,xscm(ixn),3*nbem,ilxn)
        call daread(idafdrf,iodadrf,xscm(ixa),  nbem,ilar)
      endif
c
c-----  read or form polarizability matrix if needed
c
      if ( .not. mcupdt) then
c 1-----
        if (npol .ne. 0) then
c   2-----
          nnpol3 = (npol3*(npol3+1))/2
c
          if (ixamat .eq. 0) then
c     3-----
c      -----  read -amat- from dafile
c
            call daread(idafdrf,iodadrf,xscm(ixr),nnpol3,38)
c     3-----
          else if (ixamat .eq. 1) then
c     3-----
c      -----  form -a- matrix
c
            call drfamat(xscm(ixr),1,group,npol,xpts,mpol,polar,
     1                   idrfout,afact,ithole)
c     3-----
          else if (ixamat .eq. 2) then
c     3-----
c      -----  read a-matrix from input
c
            rewind ir
   10       read(ir,'(a)',end=15) card
            if(card.ne.' $relay') goto 10
            goto 20
   15       write(iwr,*) ' no data group -$amat- found'
            call hnderr(3,errmsg)
   20       continue
            do 40, n = 1, npol
c       4-----
c        -----  find position in -a-
c
              if = (n-1)*3
              aver = zero
   25         read(ir,'(a)') card
              if (card.ne.' $end') then
c         5-----
                call freerd(card,nami,axyz,6,maxblnk)
                if (nami .eq. 'group') then
c           6-----
c            -----  put group polarizability (six components)
c                   into -relay-
c
                  do 30, k = 1, 3
                    do 30, l = 1, k
                      kl = ia(k) + l
c                     xscm(ia(if+k)+if+l) = axyz(kl)
                      xscm(ixr+ia(if+k)+if+l-1) = axyz(kl)
                      if (k .eq. l) aver = aver + axyz(kl)
   30             continue
c           6-----
                else
c           6-----
c            -----  read atom diagonal polarizabilities
c
                  do 35, k = 1, 3
c                   xscm(ia(if+k)+if+k) = axyz(k)
                    xscm(ixr+ia(if+k)+if+k-1) = axyz(k)
                    aver = aver + axyz(k)
   35             continue
c           6-----
                endif
                polar(n) = aver/three
c         5-----
              endif
c       4-----  next polarisability
   40       continue
c     3-----
          endif
c
c    -----  save lower triangle of -a- matrix
c
          call dawrit(idafdrf,iodadrf,xscm(ixr),nnpol3,38,navdrf)
          if (idrfout .eq. 3)
     1       call hatout(xscm(ixr),npol3,npol3,3,'a-mat')
cahv      ixamat = 3
c
c    -----  expand -a- to square form
c
          call hexpand(xscm(ixr),xscm(ixr),ndima,ndim)
c
c    -----  if boundary surface, construct -f-,and -d- matrices
c           these are the matrices that couple the polarisabilities
c           to the boundary elements
c
          if (ibem .ne. 0) then
c     3-----
            call drfdfmat(ieps,xscm(ixr),xscm(ixs),
     1                    xscm(ixn),xscm(ixa))
c     3-----
          endif
c   2-----
        endif
c
c  -----  construct -b- matrix, if needed
c         this matrix couples the boundary elements
c
        if (ibem .ne. 0) then
c   2-----
          if(ixbmat .eq. 1) then
c     3-----
            call drfbmat(ieps,xscm(ixr),xscm(ixs),
     1                   xscm(ixn),xscm(ixa))
	    ixwrit = igmem_alloc(ndimb*ndimb)
            call bwrit(xscm(ixr),xscm(ixwrit),npol3,ndimb,ieps)
	    call gmem_free(ixwrit)
c     3-----
          else
c     3-----
            call bread(xscm(ixr),xscm(ixs),npol3,ndimb,ieps)
c     3-----
          endif
c   2-----
        endif
c
c  -----  save relay-input matrix on dadrf file (31)
c
        call dawrit(idafdrf,iodadrf,xscm(ixr),ndim*ndim,ilr,navdrf)
c 1-----
      else
c 1-----
c  -----  (in mc-loop:) read previous relay-input matrix
c
        call daread(idafdrf,iodadrf,xscm(ixr),ndim*ndim,ilr)
c
c  -----  update relevant parts of relay-input matrix
c
        call drfaupd(xscm(ixr),ndim,npol,xpts,mpol,
     1               polar,afact,ithole)
c
c  -----  coupling between polarizabilities and boundary elements
c         only necessary if translation is included
c
        if ((notrans. eq. 0) .and. (ibem .ne. 0))
     1    call drfdfupd(ieps,xscm(ixr),xscm(ixs),
     2                  xscm(ixn),xscm(ixa))
c
c  -----  save updated relay-input matrix on dadrf file
c
        call dawrit(idafdrf,iodadrf,xscm(ixr),ndim*ndim,ilru,navdrf)
c 1-----
      endif
c
c-----  print relay-input matrix if required
c
      if (idrfout .eq. 3)
     1 call hatout(xscm(ixr),ndim,ndim,2,'relay-inp')
c
c-----  repartition memory
c       xsurf,xnorm and area are overwritten
c
      call gmem_free(ixa)
      call gmem_free(ixn)
      call gmem_free(ixs)
c       relay  at ixr
c       vv     at ixv
c       indx   at ixi
c
      ixv = igmem_alloc(ndim)
      ixi = igmem_alloc(ndim)
c     ixv = ixr + ndim*ndim
c     ixi = ixv + ndim
c     last = ixi + ndim + 1
c
c     need = loc10 + last - ixr
c     call setc(need)
c
c-----  lu-decompose the relay-input matrix; the result is put
c       back into relay
c-----  idgt is a check on the reliability of the lu-decomposition
c
      idgt =  4

      call ludatf2(xscm(ixr),ndim,ndim,
     1             idgt,d,d2,xscm(ixi),xscm(ixv),
     2             wa,ier)
c
c    ludatf should probably get different arguments for the input
c    and output matrices, even if the output should replace the 
c    input. Therefore xscm(ixr) is added.
c
c    GAMESS-UK routine ludatf (tdaf.m) requires extra argument
c 
c     call ludatf(xscm(ixr),xscm(ixr),ndim,ndim,
c    1             idgt,d,d2,xscm(ixi),xscm(ixv),
c    2             wa,ier)
caleko
c
      if (idrfout .eq. 3)
     1 call hatout(xscm(ixr),ndim,ndim,2,'relay-lu')
c
c-----  write the lu-decomposed relay matrix
c       and array -indx- on da-file 31
c
      call dawrit(idafdrf,iodadrf,xscm(ixr),ndim*ndim,illur,navdrf)
      call dawrit(idafdrf,iodadrf,xscm(ixi),ndim,ilindx,navdrf)
c
c     reset core memory
c
      call gmem_free(ixi)
      call gmem_free(ixv)
      call gmem_free(ixr)
c     call setc(loadcm)
c
      if (idrfout .ge. 1) then
        if ((.not.mcupdt) .or. (imcout .eq. 5)) 
     1   write(iwr,*) 'end of drfrelay'
      endif
      return
      end
      subroutine drfomg(ieps,xscm)
c------
c      drives the calculation of the -omega- matrix,
c      containing all formal interactions between the
c      expansion centra
c
c      -omega- is evaluated through solving the coupled
c      linear equations for unit charges and  dipoles at
c      the expansion centra, either through matrix inversion
c      (lu-decomposition) of the -relay- matrix, or
c      using an iterative (direct) method, where the interactions
c      are calculated on the fly
c------
      implicit real*8  (a-h,o-z),integer  (i-n)
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
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
      dimension xscm(*)
c
      data zero /0.0d00/
c
c-----  begin
c
      if (ieps .eq. 0) then
        ilomg = 50
        kappa = kappa1
      else
        ilomg = 51
        kappa = kappa2
      endif
c
c-----  iterative procedure or matrix inversion
c
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
c
      if (itsolv .eq. 1) then
c
c  ----- set memory partitioning for iterative procedure
c
c          xexp at ixe
c
c        ixe = 1
c       ixe = i10
c
c  -----  omega at ixo
c
c       ixo = ixe + 3*nexp
c
c  -----  dip at ixd; dip0 at ixd0
c
c       ixd = ixo + nomga
c
c  -----  polarizabilities
c
c       if (npol .ne. 0) then
c         ixd0 = ixd + 3*npol
c
c    -----  omg at ixom; omg0 at ixom0; zet at ixz; zet0 at ixz0
c
c         ixom = ixd0 + 3*npol
c
c    -----  dielectric
c
c         if (ibem .ne. 0) then
c           ixom0 = ixom + nbem
c           ixz = ixom0 + nbem
c
c    -----  xsurf at ixs; xnorm at ixn; area at ixa
c
c    -----  finite ionic strength
c
c           if (kappa .ne. zero) then
c             ixz0 = ixz + nbem
c             ixs = ixz0 + nbem
c           else
c             ixz0 = ixz + 1
c             ixs = ixz0 + 1
c           endif
c
c           ixn = ixs + 3*nbem
c           ixa = ixn + 3*nbem
c           last = ixa + nbem
c         else
c
c    -----  no dielectric
c
c           ixom0 = ixom + 1
c           ixz = ixom0 + 1
c
c    -----  xsurf at ixs; xnorm at ixn; area at ixa
c
c           ixz0 = ixz + 1
c           ixs = ixz0 + 1
c           ixn = ixs + 1
c           ixa = ixn + 1
c           last = ixa + 1
c         endif
c
c       else
c
c    -----  no polarizabilities
c
c         ixd0 = ixd + 1
c
c    -----  omg at ixom; omg0 at ixom0; zet at ixz; zet0 at ixz0
c
c         ixom = ixd0 + 1
c
c    -----  dielectric
c
c         if (ibem .ne. 0) then
c           ixom0 = ixom + nbem
c           ixz = ixom0 + nbem
c
c      -----  xsurf at ixs; xnorm at ixn; area at ixa
c
c
c      -----  finite ionic strength
c
c           if (kappa .ne. zero) then
c             ixz0 = ixz + nbem
c             ixs = ixz0 + nbem
c           else
c             ixz0 = ixz + 1
c             ixs = ixz0 + 1
c           endif
c
c           ixn = ixs + 3*nbem
c           ixa = ixn + 3*nbem
c           last = ixa + nbem
c         else
c
c      -----  no dielectric
c
c           ixom0 = ixom + 1
c           ixz = ixom0 + 1
c
c      -----  xsurf at ixs; xnorm at ixn; area at ixa
c
c           ixz0 = ixz + 1
c           ixs = ixz0 + 1
c           ixn = ixs + 1
c           ixa = ixn + 1
c           last = ixa + 1
c         endif
c       endif
c       need = loc10 + last - i10
c       call setc(need)
c
c       call omgait(ieps,xscm(ixe),xscm(ixo),xscm(ixd),xscm(ixd0),
c    1              xscm(ixom),xscm(ixom0),xscm(ixz),xscm(ixz0),
c    2              xscm(ixs),xscm(ixn),xscm(ixa))
c
      else
c
c  -----  construct the -omega- matrix, by solving the
c         linear response equations for the expanded
c         electrons, the external charges and the nuclei
c
c  -----  set memory partitioning
c
c        ixr = 1
	ixr = igmem_alloc(ndim*ndim)
	ixw = igmem_alloc(nwtr*nwtc)
	ixv = igmem_alloc(nwtr*nwtc)
	ixo = igmem_alloc(nomga)
	ixi = igmem_alloc(ndim)
c       ixr = i10
c       ixw = ixr + ndim**2
c       ixv = ixw + nwtr*nwtc
c       ixo = ixv + nwtr*nwtc
c       ixi = ixo + nomga
c       last = ixi + ndim + 1
c
c       need = loc10 + last - i10
c       call setc(need)
c
        call drfomga(ieps,xscm(ixr),xscm(ixw),xscm(ixv),
     1            xscm(ixo),xscm(ixi))
c
c-----  write -omega- on da file, record -50- or -51-
c
      call dawrit(idafdrf,iodadrf,xscm(ixo),nomga,ilomg,navdrf)
c
	call gmem_free(ixi)
	call gmem_free(ixo)
	call gmem_free(ixv)
	call gmem_free(ixw)
	call gmem_free(ixr)
      endif
c
c     reset core memory
c
c     call setc(loadcm)
c
      return
      end
         
         
         
