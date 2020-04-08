      subroutine rfin
cafc
c  this routine reads input in order to be able to perform
c  direct reaction field calculations as developed by
c  p. van duijnen et al., lab. org. & mol. inorg. chem.,
c  nijenborgh 4, 9747 ag groningen, the netherlands
c
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
c  commonblocks for drf
cafc
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
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer ibasis
      common /defpar4/ ibasis
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
      integer nprint, itol, icut, normf, normp, ifill
      common /restar/ nprint,itol,icut,normf,normp,ifill(11)
      real*8 anorm
      common /basnrm/ anorm(1024)
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer natmax, nshmax, nummax, ngsmax
      common /maxint/ natmax,nshmax,nummax,ngsmax
c
      integer maxpol, maxpnt, maxgrp, maxbem
      common /maxext/ maxpol,maxpnt,maxgrp,maxbem
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
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c  PRP/ROTA bevat de dimensionering van het array temprt(10)
c
      real*8 dipol, dmx, dmy, dmz
      common /mmtdip/ dipol,dmx,dmy,dmz
c
      real*8 dplxx, dplyy, dplzz, dplxy, dplxz, dplyz
      common /datdpl/ dplxx,dplyy,dplzz,dplxy,dplxz,dplyz
c
      real*8 fplxx, fplyy, fplzz, fplxy, fplxz, fplyz
      common /datfpl/ fplxx,fplyy,fplzz,fplxy,fplxz,fplyz
c
c
      logical expat
      integer icnexp, iexpas, iexpcn, nexpc, nexpx, nex, icch
      common /expanc2/ icnexp,iexpas,iexpcn,nexpc,
     +                 nexpx,nex,icch,expat
c
      character*16 expnam
      common /xnames/ expnam(maxex)
c
      real*8 ufact
      common /expun/ ufact
c
c
      integer nacal
      common /rota1/ nacal
c
      real*8 temprt, symfac
      common /rota2/ temprt(10), symfac
c
c
      character*80 rftex
      common /texneq/ rftex
c
      character *80 texrf
      common /neqtex/ texrf(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
c Parameter settings for the Molecular Surface program mscon
c
c maxatm     maximum number of atoms
c maxtyp     maximum number of atom types
c maxnbr     maximum number of neighbors an atom may have
c maxsph     maximum number of surface points on a sphere
c maxcir     maximum number of surface points on a circle
c maxarc     maximum number of surface points on an arc
c maxppp     maximum number of surface points per probe
c maxyon     maximum number of yon probes
c maxvic     maximum number of victim probes
c maxeat     maximum number of eaters of a probe's surface
c maxcub     maximum number of cubes in one direction
c mxconpts   maximum number of surface points
c
c maxatm must be greater than or equal to maxyon
c because they share the same cubing arrays
c
      parameter (maxatm=2000)
      parameter (maxrfa2=12000)
      parameter (maxtyp=100)
      parameter (maxnbr=200)
      parameter (maxsph=1000)
      parameter (maxcir=1000)
      parameter (maxarc=1000)
      parameter (maxppp=1000)
      parameter (maxyon=1000)
      parameter (maxvic=6000)
      parameter (maxeat=1000)
      parameter (maxcub=40)
      parameter (mxconpts=4000)

c
      character*8 scftyp
      common /scfopt2/ scftyp
c
cxxx  real*8 co,rtype,molnum,d,rp
      real*8 co,rtype,d,rp
      integer*2 molnum
      integer nmaxcon,nmincon,natom
      logical isurdens,iminpoint
      common/coninp/co(3,maxrfa2+1),rtype(maxrfa2+1),d,rp
      common/coninp2/molnum(maxrfa2+1)
      common/coninp3/natom,nmaxcon,nmincon,isurdens,iminpoint
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
      common/work/jrec,jump
      common/nmlbem/iradat,iradex
c
      integer keepij, ibeen
      common /assdrf/ keepij,ibeen
c
c
      real*8 dstgrp
      integer igroup, iunit
      common /drf_in/ dstgrp,igroup,iunit
c
      integer indstg, indstmn, indstmx, insp, inrp, inpro, inprj
      integer inbndr, inbndl, insw, insd, insphr, insr
      common /drf_un/ indstg,indstmn,indstmx,
     + insp,inrp,inpro,inprj,inbndr,inbndl,insw,insd,
     + insphr,insr
c
c
      integer irfpas, maxitrf, istci
      common /rfci1/ irfpas,maxitrf,istci
      real*8 rfcvg
      common /rfci2/ rfcvg
      logical extra
      common /rfci3/ extra
c
c     common/drf_in/dstgrp,igroup,iunit
c     common/drf_un/indstg,indstmn,indstmx,
c    1 inbndr,inbndl,ind,inrp
      common/expanx/iasexp,ifitx,iradial
      integer runflg
      common/nonos/nosopa,nodppop,nodppe,runflg
      common/toprdrf/nodrf,ifitc,iscfpol
      common/nottwi/obeen,obeen2,obeen3,obeen4
      logical obeen, obeen2, obeen3, obeen4
c
      character*8 a8, inopt1, inopt2, inopt3
      character*4 inopt11
      character*10 rfsurf
c
cxxx
      dimension qmcharge(128)
      dimension afct(2)
cxxx
      write(iwr,1000)
 1000 format(/1x,104('-')/)

      call initid
      if (obeen ) then
         obeen = .false.
         obeen2 = .false.
         obeen3 = .false.
         obeen4 = .false.
      endif
c
c  default-settings of reaction field variables
      afct(1)   = 1.662d0
      afct(2)   = 2.089d0
      afact   = 2.089d0
      agrpc   = 0.99d0
      agrpe   = 1.0d-3
      agrpm   = 0.01d0
      cvgrel  = 1.0d-4
      indstg  = 0
      dstgrp  = 15.0d0
      indstmx = 0
      dstmax  = 1.0d3
      indstmn = 0
      dstmin  = 1.0d2
      eps1    = 1.0d0
      eps2    = 1.0d0
      field   = '        '
      gamdrf  = 0.0d0
      inbndl  = 0
      hbondl  = 4.5d0
      inbndr  = 0
      hbondr  = 1.0d0
      iadexp  = 0
      ialfso  = 0
      ianal   = 1
      iasexp  = 1
      iradial = 0
      ibeen   = 0
      keepij  = 0
      ibem    = 0
      ibemout = 0
      iclintd = 1
      iclinte = 1
      iclintr = 1
      icmexp  = 0
      icnexp  = 0
      idipcal = 0
      idisadd = 0
      idrfout = 0
      ieffpol = 0
      iexpas  = 1
      iexpcn  = 0
      iexpza  = 1
      iextdip = 0
      ifitc   = 0
      ifitx   = 0
      ifldin  = 4
      ifldout = 4
      igroup  = 0
      igrppol = 0
      ihbond  = 0
      imomsav = 0
      ineqex  = 0
      intdrf  = 0
      ioponly = 0
      iradat  = 0
      iradex  = 0
      irepopt = 0
      iqmclr  = 1
      irevdis = 0
      iscfpol = 1
      isodis  = 1
      isolsav = 0
      itermax = 15
      ithole  = 2
      itsolv  = 0
      itwoeps = 0
      itwosur = 0
      iuniout = 0
      iunit   = 0
      ixamat  = 1
      ixbmat  = 1
      ixnsurf = 1
      ixrelay = 1
      ixwtvr  = 1
      ixzfp   = 1
      kappa1  = 0.0d0
      kappa2  = 0.0d0
      leveli  = 0
      levelo  = 0
      maxblnk = 1
      maxitrf = 10
      modxza  = 1
      nacal   = 0
      neqdis  = 1
      neqrep  = 1
      neqrf   = 0
      neqsta  = 0
      ngrnam  = 2
      ngrnam2 = 4
      nodiscr = 0
      nodpe   = 0
      igetden = 0
      iarfcal = 0
      nodppe  = 1
      nodppop = 1
      nodrf   = 0
      nopes   = 1
      nosopa  = 1
      rfact   = 1.0d0
      rfcvg   = 1.0d-5
      runflg  = 0
      solnam  = 'x'
      spherad = 1.0d0
      symfac  = 1.0d0
      call vclr(temprt, 1, 10)
      temprt(1) = 298.15d0
cxxx from incon
      insp = 0
      d=1.0d0
      inrp = 0
      rp=1.0d0
      inpro = 0
      inprj = 0
      insphr = 0
      insw = 0 
      insr = 0
      insd = 0
cxxx  nmaxcon=0
cxxx  isurdens=.true.
      nmaxcon=200
      isurdens=.false.
c     iminpnt=.false.
c
      radmax = 0.0d0
c
      maxci_drf = 50
      conci_drf = -1.0d0
c
      write(iwr,6666) version
c
c  read values for drf-variables as specified by user
30    call input
      if (idrfout .ge. 1) call outrec
31    call inpa(a8)
          if (a8(1:3).eq.'end') then
               goto 1100
            else if (a8(1:8).eq.'external') then
caleko
c               call extern
c There is a gamess entry called extern in casb.f
c Therefore, the name of the drf-routine is changed
c to externd
c
caleko
               call externd
            else if (a8(1:5).eq.'field') then
c  defines coupling between quantum and classical system
  20           if (jrec.ge.jump) goto 21
               call inpa4(inopt11)
               if (inopt11(1:4) .eq. 'stat') then 
                 call inpa4(field(1:4))
                 if (field(1:4) .eq. 'none') then
                   field(1:4) = '    '
                 else if ((field(1:4) .ne. 'pert') 
     1           .and. (field(1:4) .ne. 'scf ')) then
       call caserr('illegal keyword following stat')
                 endif
               else if (inopt11(1:4) .eq. 'reac') then 
                 call inpa4(field(5:8))
                 if (field(5:8) .eq. 'none') then
                   field(5:8) = '    '
                 else if ((field(5:8) .ne. 'pert') 
     1           .and. (field(5:8) .ne. 'scf ')) then
       call caserr('illegal keyword following reac')
                 endif
               else if (inopt11(1:4) .eq. 'reac') then 
                 call inpa4(field(5:8))
               else
               call caserr('unrecognised field field key')
               endif
               goto 20
  21           continue
cobsolete   else if (a8(1:5).eq.'nodrf') then
c  do not react
cobsolete      call inpoi(nodrf)
cnotimp     else if (a8(1:6).eq.'runflg') then
c  runtype control
cnotimp        call inpi(runflg)
cobsolete   else if (a8(1:7).eq.'maxblnk') then
c  maximum number of blanks in input
cobsolete      call inpi(maxblnk)
cnotimp     else if (a8(1:6).eq.'nosopa') then
c  do not perform sop calculation of polarizability
cnotimp        call inpoi(nosopa)
cnothere    else if (a8(1:7).eq.'nodppop') then
c  do not perform dipole preserving population analysis
cnothere       call inpoi(nodppop)
cnotimp     else if (a8(1:6).eq.'nodppe') then
c  perform standard dipole preserving population analysis
cnotimp        call inpoi(nodppe)
cnotimp     else if (a8(1:6).eq.'icnexp') then
c  defines expansion centra
cnotimp        call inpi(icnexp)
cnotimp     else if (a8(1:6).eq.'iexpas') then
c  defines division overlap overlap distribution
cnotimp        call inpi(iexpas)
cnotimp     else if (a8(1:5).eq.'ifitc') then
c  defines points for least squares potential fit
cnotimp        call inpi(ifitc)
cnotimp     else if (a8(1:6).eq.'iexpcn') then
c  defines dipole preserving charge centres
cnotimp        call inpi(iexpcn)
cnothere    else if (a8(1:5).eq.'nacal') then
c  calculate rotational free energy at temperatures temprt
cnothere       call inpoi(nacal)
cnothere    else if (a8(1:6).eq.'temprt') then
c  specification of temperatures for rotational free energy calculation
cnothere       do 40 i=1,10
cnothere       call inpf(temprt(i))
cnothere       if (temprt(i) .lt. 0.0d0) goto 30
cnothere 40    continue
cnothere    else if (a8(1:6).eq.'symfac') then
c  symmetryfactor for rotational free enrgy calculation
cnothere       call inpf(symfac)
cnotimp     else if (a8(1:7).eq.'iextdip') then
c  add classically induced rf to static external hamiltonian
cnotimp        call inpoi(iextdip)
            else if (a8(1:8).eq.'drftwoel') then
               call inpa(inopt1)
               if (inopt1(1:4) .eq. 'disk') then
                  if (ifill(1) .eq. 1) then
                    intdrf = 1
                  else
                    intdrf = 2
                  endif
               else if (inopt1(1:6) .eq. 'direct') then
                 intdrf = 0
               else
               call caserr('unrecognised drftwoel field key')
               endif
cxxx        else if (a8(1:6).eq.'intdrf') then
c  specification method of two-electron integral calculation
cxxx           call inpi(intdrf)
            else if (a8(1:7).eq.'inclpol') then
               call inpa(inopt1)
               if (inopt1(1:3) .eq. 'off') then
                 iscfpol=0
               else if (inopt1(1:3) .eq. 'on') then
                 iscfpol = 1
               else
               call caserr('unrecognised inclpol field key')
               endif
cxxx        else if (a8(1:7).eq.'iscfpol') then
c  optimise with polarisation costs paid in advance
cxxx           call inpoi(iscfpol)
cnotimp     else if (a8(1:7).eq.'ixrelay') then
c  compute relay matrix
cnotimp        call inpoi(ixrelay)
cnotimp     else if (a8(1:6).eq.'ixamat') then
c  specification amat treatment
cnotimp        call inpi(ixamat)
cnotimp     else if (a8(1:6).eq.'ixwtvr') then
c  compute wt and vr matrix
cnotimp        call inpoi(ixwtvr)
cnotimp     else if (a8(1:5).eq.'ixzfp') then
c  compute zfp matrix
cnotimp        call inpoi(ixzfp)
cxxx        else if (a8(1:6).eq.'modxza') then
c  use damping at short distances
cxxx           call inpoi(modxza)
            else if (a8(1:8).eq.'expandcc') then
cxxx        else if (a8(1:6).eq.'icmexp') then
c  use centre of nuclear charge as expansion centre
cxxx           call inpoi(icmexp)
               call inpa(inopt1)
               if (inopt1(1:2) .eq. 'on') then
                 icmexp=1
               else if (inopt1(1:3) .eq. 'off') then
                 icmexp = 0
               else
               call caserr('unrecognised expandcc field key')
               endif
cxxx        else if (a8(1:6).eq.'iadexp') then
c  definition of additional expansion centra
cxxx           call inpi(iadexp)
cnotimp     else if (a8(1:6).eq.'iexpza') then
c  use expanded integrals
cnotimp        call inpoi(iexpza)
cxxx        else if (a8(1:6).eq.'igroup') then
c  use group polarisabilities
cxxx           call inpoi(igroup)
            else if (a8(1:8).eq.'grouping') then
c  use group polarisabilities
               call inpa(inopt1)
               if (inopt1(1:3) .eq. 'off') then
                 igroup = 0
               else if (inopt1(1:2) .eq. 'on') then
                 igroup = 1
               else
               call caserr('unrecognised grouping field key')
               endif
               call inpi(ngrnam)
               call inpi(ngrnamb)
               if (ngrnam .lt. 0) then
                 ngrnam = 2
               else if (ngrnam .gt. 8) then
                 write(iwr,1010)
 1010            format('WARNING: NGRNAM too large,',
     1           ' reduced to 8')
                 ngrnam = 8
               endif
               if (ngrnamb .lt. 0) then
                 ngrnam2 = ngrnam + 2
               else if (ngrnamb .gt. 10-ngrnam) then
                 write(iwr,1020)
 1020            format('WARNING: NGRNAM2 too large,',
     1           ' reduced to fit length for names')
                 ngrnam2 = 10
               else
                 ngrnam2 = ngrnam + ngrnamb
               endif
            else if (a8(1:6).eq.'dstgrp') then
c  minimal distance between group and internal atoms
               call inpf(dstgrp)
               indstg = 1
            else if (a8(1:5).eq.'agrpe') then
c  maximum dif between group and internal atom idip interaction energy
               call inpf(agrpe)
            else if (a8(1:5).eq.'agrpm') then
c  maximum dif between group and internal atom idip length
               call inpf(agrpm)
            else if (a8(1:5).eq.'agrpc') then
c  minimum coincidence between group and internal atom idip direction
               call inpf(agrpc)
cxxx        else if (a8(1:6).eq.'ngrnam') then
c  defines character starting position of groupname
cxxx           call inpi(ngrnam)
            else if (a8(1:6).eq.'gamdrf') then
c  defines scaling of dispersion between quantum and classical system
               call inpf(gamdrf)
cnotimp     else if (a8(1:6).eq.'ialfso') then
c  calculate sop quantum polarisability
cnotimp        call inpoi(ialfso)
            else if (a8(1:5).eq.'units') then
cxxx        else if (a8(1:5).eq.'iunit') then
c  use angstroms (otherwise bohrs)
cxxx           call inpoi(iunit)
               call inpa4(inopt11)
               if (inopt11 .eq. 'angs') iunit = 1
            else if (a8(1:6).eq.'dstmin') then
c  minimum distance between external and quantum atoms
               call inpf(dstmin)
               indstmn = 1
            else if (a8(1:6).eq.'dstmax') then
c  maximum distance between external and quantum atoms
               call inpf(dstmax)
               indstmx = 1
cnotimp     else if (a8(1:5).eq.'ianal') then
c  perform analysis of rf energies
cnotimp        call inpoi(ianal)
cxxx        else if (a8(1:7).eq.'idrfout') then
c  specification amount of drf output
cxxx           call inpi(idrfout)
            else if (a8(1:6).eq.'drfout') then
               call inpa(inopt1)
               if (inopt1(1:4) .eq. 'some') then
                 idrfout = 1
               else if (inopt1(1:4) .eq. 'more') then
                 idrfout = 2
               else if (inopt1(1:8) .eq. 'matrices') then
                 idrfout = 3
               else if (inopt1(1:5) .eq. 'oneel') then
                 idrfout = 4
               else if (inopt1(1:5) .eq. 'twoel') then
                 idrfout = 5
               else if (inopt1(1:8) .eq. 'standard') then
                 idrfout = 0
               else
               call caserr('unrecognised drfout field key')
               endif
               if (idrfout .ge. 1) idipcal = 1
cxxx        else if (a8(1:6).eq.'ihbond') then
            else if (a8(1:5).eq.'hbond') then
c  defines use of h-bond repulsion radius
cxxx           call inpi(ihbond)
cxxx        else if (a8(1:6).eq.'hbondl') then
c  defines h-bond distance
cxxx           call inpf(hbondl)
cxxx        else if (a8(1:6).eq.'hbondr') then
c  defines h-atom radius in h-bonds
cxxx           call inpf(hbondr)
cxxx           ihbond = 1
  50           if (jrec.ge.jump) goto 51
               call inpa4(inopt11)
               if (inopt11(1:4) .eq. 'radi') then 
                  call inpf(hbondr)
                  inbndr = 1
               else if (inopt11(1:4) .eq. 'dist') then 
                  call inpf(hbondl)
                  inbndl = 1
               else if (inopt11(1:3) .eq. 'off') then 
                  ihbond = 0         
               else if (inopt11(1:2) .eq. 'on') then 
                  ihbond = 1         
               else
               call caserr('unrecognised hbond field key')
               endif
               goto 50
  51           continue
cxxx        else if (a8(1:5).eq.'afact') then
c  defines effective radius of interacting polarisabilities
cxxx           call inpf(afact)
cxxx        else if (a8(1:4).eq.'ibem') then
c  defines way of solving boundary elements equations
cxxx           call inpi(ibem)
cxxx        else if (a8(1:7).eq.'nodiscr') then
c  no discrete environment present
cxxx           call inpoi(nodiscr)
cxxx        else if (a8(1:7).eq.'idipcal') then
c  calculate dipole and quadrupole moments
cxxx           call inpoi(idipcal)
cxxx now collected under clasclas
cxxx        else if (a8(1:7).eq.'iclintd') then
c  compute classical dispersion energy
cxxx           call inpoi(iclintd)
cxxx        else if (a8(1:7).eq.'iclinte') then
c  compute classical electrostatic energy
cxxx           call inpoi(iclinte)
cxxx        else if (a8(1:7).eq.'iclintr') then
c  compute classical repulsion energy
cxxx           call inpoi(iclintr)
            else if (a8(1:8).eq.'clasclas') then
  55           if (jrec.ge.jump) goto 56
               call inpa(inopt1)
               if (inopt1(1:6) .eq. 'noelst') then 
                 iclinte = 0
               else if (inopt1(1:4) .eq. 'elst') then 
                 iclinte = 1       
               else if (inopt1(1:6) .eq. 'nodisp') then 
                 iclintd = 0       
               else if (inopt1(1:4) .eq. 'disp') then 
                 iclintd = 1       
               else if (inopt1(1:6) .eq. 'norepn') then 
                 iclintr = 0       
               else if (inopt1(1:4) .eq. 'repn') then 
                 iclintr = 1       
               else
               call caserr('unrecognised clasclas field key')
               endif
               goto 55
  56           continue
cnotimp     else if (a8(1:7).eq.'irepopt') then
c  use sop based on overlap: not operative !
cnotimp        call inpoi(irepopt)
cnotimp     else if (a8(1:6).eq.'itsolv') then
c  solve rf equations iteratively: not operative !
cnotimp        call inpoi(itsolv)
cnotimp     else if (a8(1:6).eq.'cvgrel') then
c  convergence criterion for it solution of rf equations
cnotimp        call inpf(cvgrel)
cnotimp     else if (a8(1:7).eq.'itermax') then
c  maximum number of iterations in it solution of rf equations
cnotimp        call inpi(itermax)
cnotimp     else if (a8(1:6).eq.'ifldin') then
c  defines type of source field
cnotimp        call inpi(ifldin)
cnotimp     else if (a8(1:7).eq.'ifldout') then
c  defines type of reaction field
cnotimp        call inpi(ifldout)
cnotimp     else if (a8(1:5).eq.'nodpe') then
c  use standard dipole preserving analysis
cnotimp        call inpoi(nodpe)
cnotimp     else if (a8(1:5).eq.'neqrf') then
c  add non-equilibrium polarisation energy
cnotimp        call inpoi(neqrf)
            else if (a8(1:5).eq.'neqrf') then
c  add non-equilibrium polarisation energy
               neqrf = 1
c              do ii = 1, 80
c                rftex(ii:ii) = char1(ii+6:ii+6)
c              enddo
               call inpal(rftex)
cnotimp     else if (a8(1:7).eq.'isolsav') then
c  save reaction field
cnotimp        call inpoi(isolsav)
            else if (a8(1:7).eq.'fldsave') then
c  save reaction field by induced moments
               isolsav = 1
               imomsav = 1
cnotimp     else if (a8(1:7).eq.'idisadd') then
c  add external source field
cnotimp        call inpoi(idisadd)
cxxx        else if (a8(1:6).eq.'ithole') then
c  specification functiontype for damping potentials
cxxx           call inpi(ithole)
            else if (a8(1:7).eq.'damping') then
               inafact = 0
  60           if (jrec.ge.jump) goto 61
               call inpa4(inopt11)
               if (inopt11(1:4) .eq. 'afct') then 
                  inafact = 1 
                  call inpf(afact) 
               else if (inopt11(1:4) .eq. 'cone') then 
                  modxza = 1 
                  ithole = 1 
               else if (inopt11(1:4) .eq. 'expo') then 
                  modxza = 1 
                  ithole = 2 
               else if (inopt11(1:3) .eq. 'off') then 
                  modxza = 0 
                  ithole = 0 
               else
               call caserr('unrecognised damping field key')
               endif
               goto 60
  61           continue
               if ((modxza .eq. 1) .and. (inafact .eq. 0)) 
     1            afact = afct(ithole)
            else if (a8(1:8).eq.'clasdisp') then
cxxx        else if (a8(1:6).eq.'isodis') then
c  use isotropic polarisabilities
cxxx           call inpoi(isodis)
  70           if (jrec.ge.jump) goto 71
               call inpa(inopt1)
               if (inopt1(1:6) .eq. 'noniso') then 
                  isodis = 0 
               else if (inopt1(1:6) .eq. 'isodis') then 
                  isodis = 1 
               else if (inopt1(1:6) .eq. 'effpol') then 
                  ieffpol = 1 
               else if (inopt1(1:6) .eq. 'orgpol') then 
                  ieffpol = 0 
               else if (inopt1(1:8) .eq. 'grouppol') then 
                  igrppol = 1 
               else if (inopt1(1:7) .eq. 'atompol') then 
                  igrppol = 0 
               else
               call caserr('unrecognised clasdisp field key')
               endif
               goto 70
  71           continue
            else if (a8(1:6).eq.'assign') then
cxxx        else if (a8(1:6).eq.'iasexp') then
c  defines division overlap distribution between expansion centra
cxxx           call inpi(iasexp)
               inassign = 0
  75           if (jrec.ge.jump) goto 76
               call inpa(inopt1)
               if (inopt1(1:8) .eq. 'equatocc') then
                  if (inassign .eq. 1) then
                    write(iwr,77)
  77                format('WARNING---multiple options for assign ',
     1     'directive: check input')
                    call caserr('incompatible assign field keys')
                  endif
                  inassign = 1
                  iasexp = 0
               else if (inopt1(1:8) .eq. 'olaptocc') then
                  if (inassign .eq. 1) then
                    write(iwr,77)
                    call caserr('incompatible assign field keys')
                  endif
                  inassign = 1
                  iasexp = -1
               else if (inopt1(1:6) .eq. 'midpts') then
                  if (inassign .eq. 1) then
                    write(iwr,77)
                    call caserr('incompatible assign field keys')
                  endif
                  inassign = 1
                  iasexp = -3
               else if (inopt1(1:8) .eq. 'alldistr') then
                  if (inassign .eq. 1) then
                    write(iwr,77)
                    call caserr('incompatible assign field keys')
                  endif
                  inassign = 1
                  iasexp = -7
                  iadexp = 5
               else if (inopt1(1:8) .eq. 'distance') then
                  if (inassign .eq. 1) then
                    write(iwr,77)
                    call caserr('incompatible assign field keys')
                  endif
                  inassign = 1
                  iasexp = 1
               else if (inopt1(1:8) .eq. 'radial') then
                  iradial = 1
               else if (inopt1(1:8) .eq. 'keep') then
                  keepij = 1
               else
               call caserr('unrecognised assign field key')
               endif
               goto 75
   76          continue
cnotimp     else if (a8(1:5).eq.'ifitx') then
c  defines measurement points for least squares potential fit
cnotimp        call inpi(ifitx)
            else if (a8(1:8).eq.'dielectr') then
c  specification dielectric continuum solvent
               call indiel(rfsurf)
cxxx next directives now read in/set in indiel
cxxx        else if (a8(1:6).eq.'solnam') then
c  specification dielectric continuum solvent
cxxx           call inpa(solnam)
cxxx        else if (a8(1:4).eq.'eps1') then
c  defines total dielectric constant
cxxx           call inpf(eps1)
cxxx        else if (a8(1:4).eq.'eps2') then
c  defines optic dielectric constant
cxxx           call inpf(eps2)
cxxx        else if (a8(1:6).eq.'kappa1') then
c  defines ionic strength of dielectric continuum
cxxx           call inpf(kappa1)
cxxx        else if (a8(1:6).eq.'kappa2') then
c  defines ionic strength of continuum with optic dielectric constant
cxxx           call inpf(kappa2)
cxxx        else if (a8(1:7).eq.'itwoeps') then
c  use optic dielectric constant
cxxx           call inpoi(itwoeps)
cxxx        else if (a8(1:7).eq.'ioponly') then
c  use optic dielectric constant only
cxxx           call inpoi(ioponly)
cxxx        else if (a8(1:7).eq.'itwosur') then
c  use different surfaces for total and optic dielectric constant
cxxx           call inpoi(itwosur)
cnotimp     else if (a8(1:7).eq.'ixnsurf') then
c  compute boundary surface
cnotimp        call inpoi(ixnsurf)
cnotimp     else if (a8(1:6).eq.'ixbmat') then
c  compute bmat
cnotimp        call inpoi(ixbmat)
cxxx        else if (a8(1:6).eq.'leveli') then
c  defines number of surface points per atom
cxxx           call inpi(leveli)
cxxx        else if (a8(1:6).eq.'levelo') then
c  defines number of surface triangles
cxxx           call inpi(levelo)
            else if (a8(1:6).eq.'qmradi') then
cxxx        else if (a8(1:6).eq.'iradat') then
c  defines atomic radii quantum atoms
cxxx           call inpi(iradat)
               call inpa(inopt1)
               if (inopt1(1:5) .eq. 'table') then 
                  iradat = 0 
               else if (inopt1(1:7) .eq. 'conepol') then 
                  iradat = 1 
               else if (inopt1(1:7) .eq. 'expopol') then 
                  iradat = 2 
               else
               call caserr('unrecognised qmradi field key')
               endif
            else if (a8(1:8).eq.'clasradi') then
cxxx        else if (a8(1:6).eq.'iradex') then
c  defines atomic radii external atoms
cxxx           call inpi(iradex)
               call inpa(inopt1)
               if (inopt1(1:5) .eq. 'table') then 
                  iradex = 0 
               else if (inopt1(1:7) .eq. 'conepol') then 
                  iradex = 1 
               else if (inopt1(1:7) .eq. 'expopol') then 
                  iradex = 2 
               else if (inopt1(1:7) .eq. 'userpol') then 
                  iradex = 3 
               else
               call caserr('unrecognised clasradi field key')
               endif
            else if (a8(1:7).eq.'qmclrep') then
               call inpa(inopt1)
               if (inopt1(1:3) .eq. 'off') then 
                  iqmclr = 0 
               else if (inopt1(1:2) .eq. 'on') then 
                  iqmclr = 1 
               else
               call caserr('unrecognised qmclrep field key')
               endif
cxxx        else if (a8(1:7).eq.'spherad') then
c  defines radius of sphere around centre of mass
cxxx           call inpf(spherad)
cxxx        else if (a8(1:7).eq.'ibemout') then
c  specification amount of bem output
cxxx           call inpi(ibemout)
cxxx        else if (a8(1:7).eq.'ieffpol') then
c  use effective atomic polarisabilities
cxxx           call inpoi(ieffpol)
cnotimp     else if (a8(1:7).eq.'imomsav') then
c  save induced moments
cnotimp        call inpoi(imomsav)
cnotimp     else if (a8(1:6).eq.'ineqex') then
c  include non-equilibrium component for excited states
cnotimp        call inpoi(ineqex)
cnotimp     else if (a8(1:7).eq.'irevdis') then
c  calculate reverse dispersion
cnotimp        call inpoi(irevdis)
cnotimp     else if (a8(1:7).eq.'iuniout') then
c  specification units to be used in output
cnotimp        call inpi(iuniout)
cnotimp     else if (a8(1:7).eq.'maxitrf') then
c  maximum number of ci cycles
cnotimp        call inpi(maxitrf)
cnotimp     else if (a8(1:6).eq.'neqdis') then
c  add non-equilibrium dispersion energy
cnotimp        call inpoi(neqdis)
cnotimp     else if (a8(1:6).eq.'neqrep') then
c  add non-equilibrium repulsion energy
cnotimp        call inpoi(neqrep)
cnotimp     else if (a8(1:6).eq.'neqsta') then
c  add non-equilibrium electrostatic potential
cnotimp        call inpoi(neqsta)
cnotimp     else if (a8(1:5).eq.'nopes') then
c  do not perform potential energy surface scan
cnotimp        call inpoi(nopes)
cnotimp     else if (a8(1:5).eq.'rfact') then
c  defines scaling of atomic radii in charmm model repulsion
cnotimp        call inpf(rfact)
cnotimp     else if (a8(1:5).eq.'rfcvg') then
c  convergence criterion in ci calculation
cnotimp        call inpf(rfcvg)
caleko
csindiel    else if (a8(1:8).eq.'connolly') then
c  read connolly subdirectives
csindiel       call incon 
            else if (a8(1:8).eq.'threshci') then
              call inpf(conci_drf)
            else if (a8(1:8).eq.'maxcycci') then
              call inpi(maxci_drf)
            else if (a8(1:6).eq.'densrf') then
              igetden = 1
            else if (a8(1:7).eq.'average') then
              iarfcal = 1
            else if (a8(1:6).eq.'montec') then
c  variable specification for monte carlo run
               call inmontec
            else
               call caserr('unrecognised reaction field key')
          endif
      if (jrec .ge. jump) then
        call input
        if (idrfout .ge. 1) call outrec
      endif
      go to 31
1100  continue
c
c  define internal radii
      do 400, n = 1, nat
c
c  radii extended with frecer's model (070296)
        if (anam(n) .eq. 'e' .or. anam(n) .eq. 'z') then
          radi = 1.0d-5
        elseif ((iradat .ge. 0 .and. iradat .le. 2) .or.
     +          (iradat .eq. 4)) then
          radi = radius(anam(n),iradat,qmcharge(n),afct)
        else
          call caserr('illegal option of iradat')
        endif
        if (radi .eq. 0.0d0) radi = afact
        radat(n) = radi
        alfat(n) = alfa(anam(n),ithole,ier)
        if (radi .gt. radmax) radmax = radi
        if (idrfout .eq. 3) write(iwr,*) anam(n), radat(n), alfat(n)
400   continue

caleko
c
c    Gamess-coordinates and vdWaals radii are stored to be
c    used in the generation of the Connolly surface.
c
      if (ibem .eq. 5) then
        do 121 n=1,nat
          co(1,n) = cnew(n,1)
          co(2,n) = cnew(n,2)
          co(3,n) = cnew(n,3)
          rtype(n) = radat(n)
          molnum(n) = 0
121     continue
        natom=nxtpts+nat 
      endif
c

      if (idrfout .eq. 3) write(iwr,683) radmax
683   format(/'radmax = ',f10.6/)
c
c  error checking and accounting for some locals
      if ((field(1:8) .eq. 'pert    ') .or.
     +    (field(1:8) .eq. 'scf     ')) then
        ixrelay = 0
        ixomga  = 0
      endif
c
      if ((field(5:8) .eq. 'pert') .or.
     + (field(5:8) .eq. 'scf ')) then
        ixomga = 1
        if ((field(5:8) .eq. 'scf ') .and.
     +     (ifldin .le. 3)) then
        write(iwr,7008)
7008    format(/
     +10x,'****************************************************'/
     +10x,'* expanded field required if reaction field is to **'/
     +10x,'*  be included in scf-procedure: change ifldin    **'/
     +10x,'****************************************************')
          call caserr('expanded field required: change ifldin')
        endif
        if ((field(5:8) .eq. 'scf ') .and.
     +     (ifldout .le. 3)) then
        write(iwr,7009)
7009    format(/
     +10x,'****************************************************'/
     +10x,'* expanded field required if reaction field is to **'/
     +10x,'*  be included in scf-procedure: change ifldout   **'/
     +10x,'****************************************************')
          call caserr('expanded field required: change ifldout')
        endif
      endif
c
      if ((field(1:1) .eq. ' ') .and. (runtyp .eq. 'montec')) then
        call caserr
     +  ('if you want to monte carlo, provide external field')
      endif
c
cxxx  if ((nodiscr .gt. 0) .and. (nxtpts .gt. 0)) then
cxxx    print *,'you seem to be adding a discrete part,
cxxx +   while simultaneously specifying there is none present'
cxxx  endif
      if (nxtpts .eq. 0) nodiscr = 1
      if ((ibem .eq. 0) .and. (nodiscr .eq. 1)) then
        nodrf = 1
        field = ' '
        write(iwr,8000)
8000    format(//1x,'** BEWARE: no classical atoms and no dielectric'/
     +          /1x,'** DRF module switched off *** '/)
      endif
c
c  following checks are relevant only if there is a
c  continuum/discrete boundary present
      if (ibem .ne. 0) then
c
c  if relay matrix is supposed to be read in (instead of calculated)
c  the same goes for the amat, bmat and boundary surface matrix
c  these matrices are expected to be on tape 31
        if (ixrelay.eq.0) then
           ixamat = 0
           ixbmat = 0
           ixnsurf = 0
        endif
c  compute number of boundary elements, which determines dimension
        nbem = 4**levelo*60
c
c  get total and optical dielectric constant from solvset
        if (solnam .ne. 'x') call solvset
      endif
c
c  if saving of reaction field due to charge distribution
c  is requested by user and he is using monte carlo, make
c  sure average potential at expansion centra is saved,
c  instead of induced moments at classical groups
      if ((isolsav .eq. 1) .and. (runtyp .eq. 'montec')) then
        imomsav = 0
      endif
c
c  cannot use two dielectric constants when including non-
c  equilibrium component for excited states is requested
      if ((itwoeps .eq. 1) .and. (ineqex .eq. 1)) then
        write(iwr,7007)
7007    format(/
     +10x,'****************************************************'/
     +10x,'* cannot include optical dieletric constant and   **'/
     +10x,'* non-equilibrium component for excited states    **'/
     +10x,'* at the same time                                **'/
     +10x,'****************************************************')
        call caserr
     +  ('treatment of optical dieletric constant constrained')
      endif
c
      if ((idisadd .eq. 1) .and. (ifldin .eq. 4)) then
        write(iwr,7006)
7006    format(/
     +10x,'*******************************************************'/
     +10x,'* cannot calculate source field of nuclei and electrons'/
     +10x,'* seperately if external source field is to be added'/
     +10x,'*******************************************************')
        call caserr
     +  ('cannot calculate source field of nuclei and electrons')
      endif
c
      if (iextdip .eq. 1) then
        if (field(1:3) .ne. 'scf') then
        write(iwr,7005)
7005    format(/
     +  10x,'**************************************************'/
     +  10x,'* cannot include externally induced dipoles if ***'/
     +  10x,'* static field is not included in scf-procedure **'/
     +  10x,'**************************************************')
          call caserr
     +    ('cannot include externally induced dipoles')
        else
          ixrelay = 1
        endif
        if (field(4:7) .eq. 'scf') then
        write(iwr,7004)
7004    format(/
     +  10x,'**************************************************'/
     +  10x,'* cannot include externally induced dipoles if ***'/
     +  10x,'* reaction field is included in scf-procedure  ***'/
     +  10x,'**************************************************')
          call caserr(
     +    'cannot include externally induced dipoles')
        endif
      endif
c
caleko
c
      if (field(5:7) .eq. 'scf') then
      if (intdrf.eq.0) then
        if (scftyp .eq. 'gvb') then
        write(iwr,7003)
7003    format(/
     +  10x,'**************************************************'/
     +  10x,'Option not implemented. The 2-e (D)RF integrals   '/
     +  10x,'cannot be calculated on the fly for a GVB wave    '/
     +  10x,'function, use:                                    '/
     +  10x,'   ***     SUPER ON and DRFTWOEL DISK    ***      '/
     +  10x,'**************************************************')
        call caserr('Option not implemented for GVB wavefunction')
        endif
c     else if ((intdrf.eq.1) .and. (ifill(1).eq.1)) then
c       continue
      else if ((intdrf.eq.2) .and. (ifill(1).eq.0)) then
        continue
      else
        write(iwr,7002)
7002    format(/
     +  10x,'**************************************************'/
     +  10x,'Option not implemented. The drf-contributions can '/
     +  10x,'only be calculated on the fly (DRFTWOEL DIRECT) or'/
     +  10x,'added to the gamess 2-e integrals with            '/
     +  10x,'   ***     SUPER ON and DRFTWOEL DISK    ***      '/
     +  10x,'**************************************************')
        call caserr('Option not implemented')
      endif
      endif
c 
c     if ((field(5:7) .eq. 'scf') .and. (intdrf .eq. 1) .and.
c    +    (gamdrf .ne. 1.0d0)) then
c       call caserr('cannot add reaction field contributions to
c    +   two-electron integral file: gamdrf has got to be 1.0')
c     endif
c
      if (iexpza .eq. 0) then
        write(iwr,7001)
7001    format(/
     +  10x,'****************************************************'/
     +  10x,'exact evalution of external static field integrals'/
     +  10x,'currently not possible: if due to a charge, make the'/
     +  10x,'centre ambigous'/
     +  10x,'****************************************************')
        call caserr
     +  ('exact evalution of external static field not possible')
      endif
c
      if ((irevdis .eq. 1) .and. (gamdrf .eq. 0.0d0)) then
        call caserr
     +  ('cannot calculate reverse dispersion if gamdrf is 0.0')
      endif
c
      if (gamdrf .ne. 0.0d0) then
        if (field(5:5) .eq. ' ') then
        write(iwr,7000)
7000    format(/
     +  10x,'******************************************************'/
     +  10x,'if gamdrf is not 0.0 a reaction field must be present:'/
     +  10x,'             change FIELD directive'/
     +  10x,'******************************************************')
          call caserr('gamdrf not 0.0, change FIELD directive')
        endif
        ixomga = 1
      endif
c
      if ((ixamat .eq. 0) .and. (iodadrf(1,38) .eq. 0)) then
c       call caserr('cannot read amat for she does not seem present')
        write(iwr,5555)
 5555   format(' no amat present on file: calculated anew if required')
        ixamat = 1
      endif
      if ((ixbmat .eq. 0) .and. (iodadrf(1,49) .eq. 0)) then
c       call caserr('cannot read bmat for she does not seem present')
        write(iwr,5566)
 5566   format(' no bmat present on file: calculated anew if required')
        ixbmat = 1
      endif
c
c  conversion of field from lower case to upper case (still) necessary
      if (field(1:8) .eq. 'pert    ') field(1:8) = 'pert    '
      if (field(1:8) .eq. 'pertpert') field(1:8) = 'pertpert'
      if (field(1:8) .eq. 'scf     ') field(1:8) = 'scf     '
      if (field(1:8) .eq. 'scf pert') field(1:8) = 'scf pert'
      if (field(1:8) .eq. 'scf scf ') field(1:8) = 'scf scf '
c
c  open da-files for non-equilibrium reaction fields and
c  equilibrium polarisation energies
      if ((isolsav .eq. 1) .or. (neqsta .eq. 1) .or.
     +    (neqrf .eq. 1)) then
        idafinf = 141
        call daopen(idafinf,iodainf,navinf,255)
        idafsta = 142
        call daopen(idafsta,iodasta,navsta,255)
        idafcst = 143
        call daopen(idafcst,iodacst,navcst,255)
        idafind = 144
        call daopen(idafind,iodaind,navind,255)
        idafino = 145
        call daopen(idafino,iodaino,navino,255)
        idafpol = 146
        call daopen(idafpol,iodapol,navpol,255)
        idafdis = 147
        call daopen(idafdis,iodadis,navdis,255)
        idafrep = 148
        call daopen(idafrep,iodarep,navrep,255)
        idafrqm = 149
        call daopen(idafrqm,iodarqm,navrqm,255)
        call neqread
      endif
c
c
      if (idrfout .ge. 2) call prdrfv
cxxx
cxxx  if (idrfout .ge. 1) then
        write(iwr,7701)
 7701   format(/,1x,'Overview of Direct Reaction Field Embedding')
        if (nodiscr .eq. 0) then
          write(iwr,7702) nxtpts, npol
 7702     format(/,5x,i4,' classical atoms present, of which',
     1              i4,' polarizable')
          if ((iclintr .eq. 1) .and. (ihbond .eq. 1)) write(iwr,7710)
 7710     format(5x,'CHARMM force-field H-bond repulsion')
        endif
        if (ibem .ne. 0) then
          write(iwr,7752) rfsurf
 7752     format(/,5x,'Dielectric continuum present',/,
     1    9x,'Boundary surface is ',a10, 'surface')
          if (ioponly .eq. 1) then
            write(iwr,7753) eps1, kappa1
 7753       format(9x,'Optic dielectric: epsilon= ',f8.4,
     1                 ' kappa= ',f8.4)
          else 
            write(iwr,7754) eps1, kappa1
 7754       format(9x,'Total dielectric: epsilon= ',f8.4,
     1                 ' kappa= ',f8.4)
            if (itwoeps .eq. 1)
     1      write(iwr,7753) eps2, kappa2
          endif
        endif
c
        if (iscfpol .eq. 0) then
          write(iwr,7801)
 7801     format(/,3x,'Polarization costs not in energy ',
     +    'functional: will be paid afterwards')
          scffact = 1.0d0
        else
          write(iwr,7802)
 7802     format(/,3x,'Polarization costs included ',
     +    'in energy functional')
        endif
        if (gamdrf .ne. 0.0d0) then
          write(iwr,7811) gamdrf
 7811   format(3x,'Direct RF option: QM/classical dispersion',
     1  ' scaling factor is',f8.4)
        else
          write(iwr,7812) 
 7812   format(3x,'Average RF option: no QM/classical ',
     1   'dispersion estimate included')
        endif
cxxx  endif
      write(iwr,1000)
cxxx
6666  format(/
     +/1x,t20,'*******************************************',
     +/1x,t20,'*  gamess reaction field extension        *',
     +/1x,t20,'*          version ', a12,'           *',
     +/1x,t20,'*******************************************'/)
      return
      end
c
      subroutine inpoi(logovar)
cafc
c  this routine belongs to the inp* family
c  it reads logicals disguised as integers,
c  with values assigned according to the
c  following rules : if the variable is
c  mentioned only or stated to be "on" or
c  "1" the value .true. will be assigned,
c  if it is stated to be "off" or "0" the
c  value .false. will be assigned
c
      implicit real*8  (a-h,o-z),integer  (i-n)
c
      common/work/jrec,jump
c
      integer     logovar
      character*4 charac
c
      logovar = 1
      call inpa4(charac)
      if (charac(1:3) .eq. 'off' .or. charac(1:1) .eq. '0') then
          logovar = 0
      else if (charac(1:2) .ne. 'on' .and. charac(1:1) .ne. '1'
     +  .and. charac(1:1) .ne. ' ') then
          jrec = jrec - 1
      endif
      return
      end
c
      subroutine initid
cafc
c  this routine performs some initialising
c
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
c  commonblocks for drf
cafc
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
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer ibasis
      common /defpar4/ ibasis
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
      integer nprint, itol, icut, normf, normp, ifill
      common /restar/ nprint,itol,icut,normf,normp,ifill(11)
      real*8 anorm
      common /basnrm/ anorm(1024)
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer natmax, nshmax, nummax, ngsmax
      common /maxint/ natmax,nshmax,nummax,ngsmax
c
      integer maxpol, maxpnt, maxgrp, maxbem
      common /maxext/ maxpol,maxpnt,maxgrp,maxbem
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
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c  PRP/ROTA bevat de dimensionering van het array temprt(10)
c
      real*8 dipol, dmx, dmy, dmz
      common /mmtdip/ dipol,dmx,dmy,dmz
c
      real*8 dplxx, dplyy, dplzz, dplxy, dplxz, dplyz
      common /datdpl/ dplxx,dplyy,dplzz,dplxy,dplxz,dplyz
c
      real*8 fplxx, fplyy, fplzz, fplxy, fplxz, fplyz
      common /datfpl/ fplxx,fplyy,fplzz,fplxy,fplxz,fplyz
c
c
      logical expat
      integer icnexp, iexpas, iexpcn, nexpc, nexpx, nex, icch
      common /expanc2/ icnexp,iexpas,iexpcn,nexpc,
     +                 nexpx,nex,icch,expat
c
      character*16 expnam
      common /xnames/ expnam(maxex)
c
      real*8 ufact
      common /expun/ ufact
c
c
      integer nacal
      common /rota1/ nacal
c
      real*8 temprt, symfac
      common /rota2/ temprt(10), symfac
c
c
      character*80 rftex
      common /texneq/ rftex
c
      character *80 texrf
      common /neqtex/ texrf(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
c
      character *24 ver0
      data ver0/'3.0.2 090696'/
c
      nxtpts = 0
      ngran  = 1
c
      nbem   = 0
      nbem1  = 0
      nbem2  = 0
      nummax = mxnum
      nshmax = mxsh
      natmax = mxat
      ngsmax = mxgs
      maxpol = mxpol
      maxpnt = mxpts
      maxgrp = mxgrp
      maxbem = 1000
      maxneq = mxneq
      mxpair = max(nummax,3*maxpnt)
      numint = num
      natint = nat
      nshint = nshell
      ngauss = kstart(nshell) + kng(nshell) - 1
      loc    = kloc(nshell) + kmax(nshell) - kmin(nshell)
c
      scffact = 0.5d0
      nneq    = 1
      nneqrf  = 0
      neq2eps = 0
      ixomga  = 0
      ixbw    = 1
      ixexp   = 1
      isur1   = 52
      isur2   = 53
      call vclr(rdisp, 1, mxst)
      call vclr(upolneq, 1, mxst)
c
      version = ver0
      call inithondo(0)
c
      idafdrf = 131
      call daopen(idafdrf,iodadrf,navdrf,100)
c
c  open da-files for non-equilibrium reaction fields and
c  equilibrium polarisation energies
c     if ((isolsav .eq. 1) .or. (neqsta .eq. 1) .or.
c    +    (neqrf .eq. 1)) then
c       idafind = 44
c       call daopen(idafind,iodaind,navind,255)
c       idafino = 45
c       call daopen(idafino,iodaino,navino,255)
c       idafpol = 46
c       call daopen(idafpol,iodapol,navpol,255)
c       call neqread
c     endif
      return
      end
c
      subroutine prdrfv
cafc
c  this routine prints out the values of all drf-variables
c  the print_out is ordened in groups based on type of the
c  variable and the original namelist to whom it belonged
c
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
c  commonblocks for drf
cafc
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
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer ibasis
      common /defpar4/ ibasis
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
      integer nprint, itol, icut, normf, normp, ifill
      common /restar/ nprint,itol,icut,normf,normp,ifill(11)
      real*8 anorm
      common /basnrm/ anorm(1024)
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer natmax, nshmax, nummax, ngsmax
      common /maxint/ natmax,nshmax,nummax,ngsmax
c
      integer maxpol, maxpnt, maxgrp, maxbem
      common /maxext/ maxpol,maxpnt,maxgrp,maxbem
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
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c  PRP/ROTA bevat de dimensionering van het array temprt(10)
c
      real*8 dipol, dmx, dmy, dmz
      common /mmtdip/ dipol,dmx,dmy,dmz
c
      real*8 dplxx, dplyy, dplzz, dplxy, dplxz, dplyz
      common /datdpl/ dplxx,dplyy,dplzz,dplxy,dplxz,dplyz
c
      real*8 fplxx, fplyy, fplzz, fplxy, fplxz, fplyz
      common /datfpl/ fplxx,fplyy,fplzz,fplxy,fplxz,fplyz
c
c
      logical expat
      integer icnexp, iexpas, iexpcn, nexpc, nexpx, nex, icch
      common /expanc2/ icnexp,iexpas,iexpcn,nexpc,
     +                 nexpx,nex,icch,expat
c
      character*16 expnam
      common /xnames/ expnam(maxex)
c
      real*8 ufact
      common /expun/ ufact
c
c
      integer nacal
      common /rota1/ nacal
c
      real*8 temprt, symfac
      common /rota2/ temprt(10), symfac
c
c
      character*80 rftex
      common /texneq/ rftex
c
      character *80 texrf
      common /neqtex/ texrf(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
c
      integer irfpas, maxitrf, istci
      common /rfci1/ irfpas,maxitrf,istci
      real*8 rfcvg
      common /rfci2/ rfcvg
      logical extra
      common /rfci3/ extra
c
c
      common/nmlbem/iradat,iradex
c
      real*8 dstgrp
      integer igroup, iunit
      common /drf_in/ dstgrp,igroup,iunit
c
      integer indstg, indstmn, indstmx, insp, inrp, inpro, inprj
      integer inbndr, inbndl, insw, insd, insphr, insr
      common /drf_un/ indstg,indstmn,indstmx,
     + insp,inrp,inpro,inprj,inbndr,inbndl,insw,insd,
     + insphr,insr
c
c     common/drf_in/dstgrp,igroup,iunit
      common/expanx/iasexp,ifitx,iradial
      integer runflg
      common/nonos/nosopa,nodppop,nodppe,runflg
      common/toprdrf/nodrf,ifitc,iscfpol
c
      write(iwr,12)
12    format(//'drf variables and their values :'/)
c
c cntrl
      write(iwr, 1) nodrf, runflg, maxblnk
1     format('nodrf = ',i1,'  runflg = ',i2,'  maxblnk = ',i2/)
c
c prp
      write(iwr, 2) nacal, nodppe, nodppop, nosopa, symfac, temprt(1)
2     format('nacal = ',i1,'  nodppe = ',i1,'  nodppop = ',i1,
     +'  nosopa = ',i1,'  symfac = ',f4.2,'  temprt(1) = ',f6.2/)
c
c expanc
      write(iwr,3) icnexp, iexpas, iexpcn, ifitc
3     format('icnexp = ',i2,'  iexpas = ',i2,'  iexpcn = ',i2,
     +'  ifitc = ',i2/)
c
c expanx
      write(iwr,4) iasexp, ifitx
4     format('iasexp = ',i2,'  ifitx = ',i2/)
c
c bem
      write(iwr,5) ibemout, ioponly, itwoeps, itwosur, iuniout, ixbmat,
     +ixnsurf, leveli, levelo
5     format('ibemout = ',i2,'  ioponly = ',i1,'  itwoeps = ',i1,
     +'  itwosur = ',i1,'  iuniout = ',i2,'  ixbmat = ',i1,
     +'  ixnsurf = ',i1,'  leveli = ',i2,'  levelo = ',i2/)
      write(iwr,6) solnam
6     format('solnam = ',a15/)
      write(iwr,7) eps1, eps2, kappa1, kappa2, spherad
7     format('eps1 = ',f7.3,'  eps2 = ',f7.3,'  kappa1 = ',f7.3,
     +'  kappa2 = ',f7.3,'  spherad = ',f7.3/)
c
c drf
      write(iwr,8) ialfso, ianal, iclintd, iclinte, iclintr, icmexp,
     +idipcal, idisadd, ieffpol, iexpza, iextdip, igroup, imomsav,
     +ineqex, irepopt, irevdis, iscfpol, isodis, isolsav, itsolv,
     +iunit, ixrelay, ixwtvr, ixzfp, modxza, neqdis, neqrep, neqrf,
     +neqsta, nodiscr, nodpe, nopes
8     format('ialfso = ',i1,'  ianal = ',i1,'  iclintd = ',i1,
     +'  iclinte = ',i1,'  iclintr = ',i1,'  icmexp = ',i1,
     +'  idipcal = ',i1,'  idisadd = ',i1,'  ieffpol = ',i1,/
     +'iexpza = ',i1,'  iextdip = ',i1,'  igroup = ',i1,
     +'  imomsav = ',i1,'  ineqex = ',i1,'  irepopt = ',i1,
     +'  irevdis = ',i1,
     +'  iscfpol = ',i1,'  isodis = ',i1,/'isolsav = ',i1,
     +'  itsolv = ',i1,'  iunit = ',i1,'  ixrelay = ',i1,'  ixwtvr = '
     +,i1,'  ixzfp = ',i1,'  modxza = ',i1,'  neqdis = ',i1,
     +'  neqrep = ',i1,'  neqrf = ',i1,/'neqsta = ',i1,'  nodiscr = '
     +,i1,'  nodpe = ',i1,'  nopes = ',i1/)
      write(iwr,9) iadexp, ibem, idrfout, ifldin, ifldout, intdrf,
     +iradex, iradat, itermax, ithole, ixamat, maxitrf, ngrnam
9     format('iadexp = ',i2,'  ibem = ',i2,'  idrfout = ',i2,
     +'  ifldin = ',i2,'  ifldout = ',i2,'  intdrf = ',i2,'  iradex = '
     +,i2,'  iradat = ',i2,'  itermax = ',i2,/'ithole = ',i2,
     +'  ixamat = ',i2,'  maxitrf = ',i2,'  ngrnam = ',i2/)
      write(iwr,10) field
10    format('field = ',a8/)
      write(iwr,11) afact, agrpc, agrpe, agrpm, cvgrel, dstgrp, dstmax,
     +dstmin, gamdrf, hbondl, hbondr, rfact, rfcvg
11    format('afact = ',f7.3,'  agrpc = ',f7.3,'  agrpe = ',f7.3,
     +'  agrpm = ',f7.3,'  cvgrel = ',f7.3,'  dstgrp = ',f7.3,
     +'  dstmax = ',f8.3,/'dstmin = ',f7.3,'  gamdrf = ',f7.3,
     +'  hbondl = ',f7.3,'  hbondr = ',f7.3,'  rfact = ',f7.3,
     +'  rfcvg = ',f7.3/)
      return
      end
c
      subroutine externd
cafc
c
c  this routine reads information about the classical external atoms
c
c  every card contains the info about an atom formatted as follows :
c  chem symbol - group name - charge - x - y - z - alfa - (radius)
c
c  note for programmers :
c  due to historical reasons the chem symbol and group name are
c  stored in the char*16 variable namei, which puts the chem
c  symbol in position 1:2 and the group name in position ngrnam:end
c
c  in addition the following directives can be used :
c
c group --> force grouping of polarisabilities of atoms in next
c           group, regardless of current settings
c
c qmcharge --> read qmcharges as calculated in a previous job
c              expected only when iradex is set to 4
c              it has to be followed by a card for every atom
c              specifying its name (as above) and charge
c
c tape n --> read info from tape n instead of normal input
c            specification file
c
c amb x y n --> next atom is treated amibuously, which means
c               that it remains classical for the classical
c               atoms, while becoming part of the qm system also
c
c   x -->	
c   '    ': use the standard basis set for the original qm atoms
c            if this means no standard, further definition is required
c            as in $basis on the cards following the name, charge etc.
c
c   ' sto': sto-3g basis                          (no ecp)
c   ' n21': 3-21g basis                           (no ecp)
c   ' n31': 4-31g basis                           (no ecp)
c   ' dzv': dunning's (9,5) [3,2] basis           (no ecp)
c   ' dzp': ibid. with polarisation function      (no ecp)
c   ' tzv': triple zeta basis                     (no ecp)
c   ' tzp': ibid. with polarisation function      (no ecp)
c   'ecdz': compact ecp with double zeta basis
c   'ecdp': compact ecp with double zeta basis and polarisation function
c   'ectz': compact ecp with triple zeta basis
c   'ectp': compact ecp with triple zeta basis and polarisation function
c   'mini': huzinaga's mini-3 basis               (no ecp)
c   'midi': huzinaga's midi-3 basis               (no ecp)
c
c   y -->
c    '   ': no ecp unless implicit in the basis set defined
c    'ecp': the centre will get an ecp
c           (1) if the basis set defined (global or local) equals
c               dzv, dzp, tzv or tzp, the basis will automatically
c               be changed to ecdz, ecdp, ectz or ectp, respectively
c           (2) if not, ecp input is expected on $ecp
c               (in the order of input atoms)
c
c   n -->
c           determines the number of electrons to be
c           added to the qm system for this ambiguous atom
c           if ecp's are defined, the correct number of core electrons
c           is deducted automatically.
c
c           example: add a c-atom with an ecp and 1 valence electron
c
c           amb ecdz 3
c           c 123     0.0  3.35  6.17  11.0
c
c           note: in case of there being ambiguous atoms
c                 -modxza- is (changed to) 1
c           this is necessary to saturate "dangling" bonds
c
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
c Parameter settings for the Molecular Surface program mscon
c
c maxatm     maximum number of atoms
c maxtyp     maximum number of atom types
c maxnbr     maximum number of neighbors an atom may have
c maxsph     maximum number of surface points on a sphere
c maxcir     maximum number of surface points on a circle
c maxarc     maximum number of surface points on an arc
c maxppp     maximum number of surface points per probe
c maxyon     maximum number of yon probes
c maxvic     maximum number of victim probes
c maxeat     maximum number of eaters of a probe's surface
c maxcub     maximum number of cubes in one direction
c mxconpts   maximum number of surface points
c
c maxatm must be greater than or equal to maxyon
c because they share the same cubing arrays
c
      parameter (maxatm=2000)
      parameter (maxrfa2=12000)
      parameter (maxtyp=100)
      parameter (maxnbr=200)
      parameter (maxsph=1000)
      parameter (maxcir=1000)
      parameter (maxarc=1000)
      parameter (maxppp=1000)
      parameter (maxyon=1000)
      parameter (maxvic=6000)
      parameter (maxeat=1000)
      parameter (maxcub=40)
      parameter (mxconpts=4000)

c  commonblocks for drf
cafc
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
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer ibasis
      common /defpar4/ ibasis
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
c
      integer nprint, itol, icut, normf, normp, ifill
      common /restar/ nprint,itol,icut,normf,normp,ifill(11)
      real*8 anorm
      common /basnrm/ anorm(1024)
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
      integer ir, iwr, ipnch,ipadiofil
      common /iofile/ ir,iwr,ipnch,ipadiofil(20)
c
      integer list
      common /hlistng/ list
c
c
      integer natmax, nshmax, nummax, ngsmax
      common /maxint/ natmax,nshmax,nummax,ngsmax
c
      integer maxpol, maxpnt, maxgrp, maxbem
      common /maxext/ maxpol,maxpnt,maxgrp,maxbem
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
c
      integer nopes, lpes, lpesc, intpesc, numpes, npesgrp, ipesout
      common /pescl/ nopes,lpes,lpesc,intpesc,numpes,npesgrp,ipesout
c
c
c  PRP/ROTA bevat de dimensionering van het array temprt(10)
c
      real*8 dipol, dmx, dmy, dmz
      common /mmtdip/ dipol,dmx,dmy,dmz
c
      real*8 dplxx, dplyy, dplzz, dplxy, dplxz, dplyz
      common /datdpl/ dplxx,dplyy,dplzz,dplxy,dplxz,dplyz
c
      real*8 fplxx, fplyy, fplzz, fplxy, fplxz, fplyz
      common /datfpl/ fplxx,fplyy,fplzz,fplxy,fplxz,fplyz
c
c
      logical expat
      integer icnexp, iexpas, iexpcn, nexpc, nexpx, nex, icch
      common /expanc2/ icnexp,iexpas,iexpcn,nexpc,
     +                 nexpx,nex,icch,expat
c
      character*16 expnam
      common /xnames/ expnam(maxex)
c
      real*8 ufact
      common /expun/ ufact
c
c
      integer nacal
      common /rota1/ nacal
c
      real*8 temprt, symfac
      common /rota2/ temprt(10), symfac
c
c
      character*80 rftex
      common /texneq/ rftex
c
      character *80 texrf
      common /neqtex/ texrf(mxneq+1)
c
      integer index, ineqix
      common /indxneq/ index,ineqix(mxneq+1)
c
c
      real*8 trmatr,trvec
      common /trnsf/ trmatr(3,3),trvec(3)
c
cxxx  real*8 co,rtype,molnum,d,rp
      real*8 co,rtype,d,rp
      integer*2 molnum
      integer nmaxcon,nmincon,natom
      logical isurdens,iminpoint
      common/coninp/co(3,maxrfa2+1),rtype(maxrfa2+1),d,rp
      common/coninp2/molnum(maxrfa2+1)
      common/coninp3/natom,nmaxcon,nmincon,isurdens,iminpoint
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
      common/work/jrec,jump
      common/nmlbem/iradat,iradex
c
      real*8 dstgrp
      integer igroup, iunit
      common /drf_in/ dstgrp,igroup,iunit
c
      integer indstg, indstmn, indstmx, insp, inrp, inpro, inprj
      integer inbndr, inbndl, insw, insd, insphr, insr
      common /drf_un/ indstg,indstmn,indstmx,
     + insp,inrp,inpro,inprj,inbndr,inbndl,insw,insd,
     + insphr,insr
c
c     common/drf_in/dstgrp,igroup,iunit
c     common/drf_un/indstg,indstmn,indstmx,
c    1 inbndr,inbndl,ind,inrp
c
c  locals
c  additional array for qm charges from a previous job
c  necessary only when iradex == 4
      dimension qmcharge(128)
c
      logical group, fgroup, mcgroup, nqgroup
      character*16 namei
      character*80 card, bufstrg
      dimension afct(2), cxyz(6), qq(3)
      character*10 namgrp
      logical amb, far, newgrp, newanal, okqmch
      character*3 ecpwrd
      integer jnear
      dimension mgran(mxgran)
      logical oreanal
      character*4 a4
      character*8 a8
c
      character*4 nbs0, nbs1, nbs2, nbs3, nbs4, nbs5, nbs6, nbs7,
     +            nbs8, nbs9, nbs10, nbs11, nbs12, nbs13
                                                        
c  initialisation of various variables
c
c     data  nbas0, nbas1, nbas2, nbas3 /'    ',' sto',' n31',' n21'/
c     data  nbas4, nbas5, nbas6, nbas7 /' dzv',' dzp',' tzv',' tzp'/
c     data  nbas8, nbas9,nbas10,nbas11 /'ecdz','ecdp','ectz','ectp'/
c     data nbas12,nbas13               /'mini','midi'/
c
      data  nbs0, nbs1, nbs2, nbs3 /'    ',' sto',' n31',' n21'/
      data  nbs4, nbs5, nbs6, nbs7 /' dzv',' dzp',' tzv',' tzp'/
      data  nbs8, nbs9,nbs10,nbs11 /'ecdz','ecdp','ectz','ectp'/
      data nbs12,nbs13             /'mini','midi'/
c
      data ecpwrd /'ecp'/
c
c     data keyecp, ecpbl /' ecp', '    '/
c
      data bohr /0.529177249d0/
      data afct /1.662d0, 2.089d0/
c
      card = ' '
      fgroup = .false.
      group = igroup .ne. 0
      igran = 1
      mcgroup = .true.
      namei = ' '
      nqgroup = .false.
      ngrp = 0
      ngroup = 0
      ngrpol = 0
      ngrpmc = 0
      npolg = 0
      npoints = 0
      ndisc = 0
      npolin = 0
      npol = 0
      namb = 0
      okqmch = .false.
cxxx  radmax = 0.0d0
      maxbl = maxblnk
      amb = .false.
      ztot = 0.0d0
      ier = 0
      ner = 0
      fact1 = 1.0d0
caleko
c     These variables should only be converted to bohrs if 
c     they have been input in angstroms. If you convert their
c     default values you enlarge them by 1.89 while not 
c     necessary.
c
      if (iunit .gt. 0) fact1 = 1.0d0/bohr
cxxx
c  if user has specified default values, they
c  should still be converted!
c      if (dstmin .ne. 1.0d2) dstmin = dstmin*fact1
c      if (dstmax .ne. 1.0d3) dstmax = dstmax*fact1
c      if (hbondl .ne. 4.5) hbondl = hbondl*fact1
c      if (hbondr .ne. 1.0) hbondr = hbondr*fact1

       if (indstg .eq. 1) dstgrp = dstgrp*fact1
       if (indstmn .eq. 1) dstmin = dstmin*fact1
       if (indstmx .eq. 1) dstmax = dstmax*fact1
       if (inbndl .eq. 1) hbondl = hbondl*fact1
       if (inbndr .eq. 1) hbondr = hbondr*fact1
caleko

      fact3 = fact1**3
      irdext = ir
c
c  initialize ambiguous point arrays
      do 50, i = 1, natmax
        nambpt(i) = 0
50    continue
      do 55, i = 1, maxpnt
        ncutpt(i) = 0
        imemgrp(i) = 0
55    continue
      do 56, i = 1, mxgran
        neqgrp(i) = 1
 56   continue
      call vclr(znamb, 1, natmax)
      call vclr(atpol, 1, maxpnt)
      call vclr(atpolt, 1, 6*maxpnt)
      call vclr(cxyz, 1, 6)
c
c  code
      if (afact .eq. 0.0d0) afact = afct(ithole)
c
c
c  if boundary elements only, skip processing external points
cxxx  if (nodiscr .gt. 0) goto 360
c

      if (idrfout .ge. 1) write(iwr,665)
c
10    call input
      if (idrfout .ge. 1) call outrec
      call inpa4(a4)
      card(1:4) = a4(1:4)
c
c
      if (a4(1:3) .eq. 'amb') then
        amb = .true.
c
caleko
c      From this point on, everything that tries to convey information
c      to the qm system about the ambiguous atoms is taken out.
c      The qm system ALREADY knows.
c
 
c  nr of electrons to be added to the qm system for this amb atom
c        call inpi(nambel)
c        ambel = 1.0d0 * nambel
c
c  basis defined for this amb atom ?
c        call inpa4(a4)
c        if (a4(1:1) .eq. ' ') goto 10
c        if (a4(1:4) .eq. nbs0) then
c
c  employ the same basis as for the qm defined atoms
c  if the basis for the qm atoms is given on input
c  the same is requested for the amiguous atoms
c  format: see $basis
c  the ecp's are expected on input as required for qm
c  treated atoms: see $ecp
c          ibasx = ibasis
c          iecpx = ecpbl
c          call inpa4(a4)
c          if (a4(1:3) .eq. ecpwrd) then
c            iecpx = keyecp
c            if (ibasx .eq. nbas4) ibasx = nbas8
c            if (ibasx .eq. nbas5) ibasx = nbas9
c            if (ibasx .eq. nbas6) ibasx = nbas10
c            if (ibasx .eq. nbas7) ibasx = nbas11
c          endif
c        else if (a4(1:4) .eq. nbs1) then
c          ibasx = nbas1
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs2) then
c          ibasx = nbas2
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs3) then
c          ibasx = nbas3
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs4) then
c          ibasx = nbas4
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs5) then
c          ibasx = nbas5
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs6) then
c          ibasx = nbas6
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs7) then
c          ibasx = nbas7
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs8) then
c          ibasx = nbas8
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs9) then
c          ibasx = nbas9
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs10) then
c          ibasx = nbas10
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs11) then
c          ibasx = nbas11
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs12) then
c          ibasx = nbas12
c          iecpx = ecpbl
c        else if (a4(1:4) .eq. nbs13) then
c          ibasx = nbas13
c          iecpx = ecpbl
c        else
c          call caserr('basis for ambigous atom not valid')
c        endif
      goto 10
      endif
c
      if (a4(1:4) .eq. 'grou') then
        fgroup = .true.
        group = .true.
caleko
c       There is no checking on variable igroup, therefore there is always
c       a grouping of the classical group's polarizabilities.
c       In this way you should be able to define a group and still use 
c       individual polarisabilities.
cxxx    The point here is that grouping is FORCED, regardless 
c       of other settings, therefore fgroup and group
c       are set true
cxxx    if (igroup.eq.0) then
cxxx      fgroup = .false.
cxxx      group = .false.
cxxx    endif
c
caleko

6767    call inpa4(a4)
        if (a4(1:4) .ne. ' ') then
          if (a4(1:4) .eq. 'nomc') mcgroup = .false.
          if (a4(1:4) .eq. 'anal') call inpi(igran)
          if (a4(1:4) .eq. 'neqs') nqgroup = .true.
          goto 6767
        endif
        goto 10
      endif
c
      if (a4(1:4) .eq. 'tape') then
        call inpi(irdext)
        write(iwr,666) irdext
        goto 10
      endif
c
      if (a4(1:4) .eq. 'qmch') then
        okqmch = .true.
        do 96, n = 1, nat
          call input
          if (idrfout .ge. 1) call outrec
          call inpa4(namei(1:2))
cxxx      call inpa(namei(ngrnam:))
cxxx      call inpa(namei(7:))
          call inpan(namei(7:))
          call inpf(cxyz(1))
          qmcharge(n) = cxyz(1)
96      continue
        goto 10
      endif
c
      if ((a4(1:3) .eq. 'end') .and. (nxtpts .eq. 0)) goto 300
c  reading section
      namei(1:2) = a4(1:2)
cxxx  call inpa(namei(7:))
      call inpan(namei(7:))
      call inpf(cxyz(1))
      call inpf(cxyz(2))
      call inpf(cxyz(3))
      call inpf(cxyz(4))
      call inpf(cxyz(5))
      if (jrec .lt. jump) then
        call inpf(cxyz(6))
      else
        cxyz(6) = 0.0d00
      endif
c
      charge = cxyz(1)
      alf = cxyz(5)
      radi = cxyz(6)
c
      if ((namei(1:2) .eq. '  ') .and. (charge .eq. 0.0d0)) then
        ndisc = ndisc + 1
        if (idrfout .ne. 0) then
          if (ndisc .eq. 1) write (iwr,667) namei, (cxyz(i),i=2,4)
        endif
        goto 10
      endif
c
c  convert to bohr
      cxyz(2) = fact1*cxyz(2)
      cxyz(3) = fact1*cxyz(3)
      cxyz(4) = fact1*cxyz(4)
      alf = alf*fact3
      radi = radi*fact1
c
caleko
c     The drf-coordinates are transformed in the same way
c     as the qm-system coordinates.
c
     
* NOSYMM NOT DEFINED ...
c     if (nosymm .eq. 0) then 
        call transform(cxyz(2),trvec,trmatr)
c     endif
c
caleko

c  check whether this is an acceptable point
      if (nxtpts .eq. 0) then
cxxx    namgrp = namei((7+ngrnam):(7+ngrnam2))
        namgrp = namei(7:(6+ngrnam2))
      endif
cxxx  if (nxtpts .eq. 0) namgrp = namei(7:)
      newgrp = 
     +  (namei((7+ngrnam):(6+ngrnam2)) 
     +   .ne. namgrp(ngrnam+1:ngrnam2)) 
     1         .or. (card(1:3) .eq. 'end')
     2         .or. fgroup
      if (newgrp) then
c  process completed group
c
c  save number of points making up this group
        if (namgrp(ngrnam+1:ngrnam2) .eq. ' ') 
     1    namgrp = namgrp(1:ngrnam)//'gr'
        dnear = 1.0d6
        do 100, j = 1, nat
c
c  exclude amiguous atoms (their field is not calculated
c  at this group anyway)
          if ((nambpt(j) .ge. (nxtpts - ngroup + 1)) .and.
     +        (nambpt(j) .le. nxtpts)) goto 100
          do 90, k =nxtpts-ngroup+1, nxtpts
            call distab(c(1,j),xpts(1,k),qq,dist)
            if (dist .lt. dnear) then
              dnear = dist
              jnear = j
            endif
90        continue
100     continue
c
c  check
        if (dnear .ge. dstmax) then
c
c  complete group is too far away: remove all points from list
          write(iwr,668) ngroup, namgrp
          nxtpts = nxtpts - ngroup
          npol = npol - npolg
          ndisc = ndisc + ngroup
        else if ((group .and. dnear .gt. dstgrp) .or. fgroup) then
c
c  process completed group
          if (igran .gt. 1) then
            newanal = .true.
            do 200, i = 1, nxtpts - 1
              if (igran .eq. igranl(i)) newanal = .false.
200         continue
            if (newanal) ngran = ngran + 1
          endif
          if (.not. nqgroup) then
            neqgrp(igran) = 0
          endif
c
c    -----  check members of group
c
          if (ngroup .gt. mxnpts) then
            write(iwr,*) 
     +      'too many atoms defining a group: the maximum is', 
     +       mxnpts
            call caserr
     +   ('too many atoms in a group: change mxnpts in rfin/sizesrf')
          endif
          call drfgrp(idrfout,nxtpts-ngroup+1,nxtpts,ngroup+1,jnear,
     +                 npol-npolg+1,npol,ngrp,igran,ngrnam,
     +                 namgrp,fgroup,
     +                 idipcal,afact,ithole,agrpc,agrpm,agrpe,
     +                 ieffpol)
          ngrppt(ngrp) = ngroup
          igrpst(ngrp) = nxtpts - ngroup
          igrlast(ngrp) = nxtpts
          write(grlabel(ngrp),669) bufstrg
         if (mcgroup) then
            if (ngroup .ge. mxnpts) then
              write (iwr,670) mxnpts, ngrp
              call caserr
     +   ('reduce group size or enlarge mxnpts in ctlnew.f and dimpar')
            endif
            ngrpmc = ngrpmc + 1
            ibitmc(ngrpmc) = ngrp
          endif
        else if (group .and. (dnear .le. dstgrp)) then
          write(iwr,6681) ngroup, namgrp, dnear, dstgrp
        endif
c
c  reset counters for next group
        ngroup = 0
        npolg = 0
cxxx    namgrp = namei((7+ngrnam):(7+ngrnam2))
cxxx    namgrp = namei(7:)
        namgrp = namei(7:(6+ngrnam2))
        fgroup = .false.
        if (igroup .eq. 0) group = .false.
        mcgroup = .true.
        igran = 1
        nqgroup = .false.
      endif
c
c  store points
      if (card(1:3) .eq. 'end') goto 300
      nxtpts = nxtpts + 1
      ngroup = ngroup + 1
      npoints = npoints + 1
      if (nxtpts .gt. maxpnt) then
        write(iwr,671) nxtpts, maxpnt
        call caserr('try again')
      endif
      if (amb) then
        namb = namb + 1
        if ((nat + 1) .gt. natmax) then
          write(iwr,672)
          call caserr('try again')
        endif
c        call addamb(ambel,charge,cxyz(2),namei,ibasx,iecpx)
        nambpt(nat) = nxtpts
        ncutpt(nxtpts) = nat
        amb = .false.
      endif
      if (alf .le. 0.0d0) alf = alfa(namei,ithole,ier)
      if ((namei(1:2) .ne. 'qq') .and. (radi .eq. 0.0d0)) then
        if (iradex .eq. 3) then
          radi = afact*alf**(1.0d0/3.0d0)
        else
c
c  new computation of radius (070296)
          radi = radius(namei,iradex,charge,afct)
        endif
        if (radi .eq. 0.0d0) radi = afact
      endif
      if ((namei(1:2) .ne. 'qq') .and. (namei(1:2) .ne. '  ')) then
        npol = npol + 1
        npolg = npolg + 1
        if (npol .gt. maxpol) then
          write(iwr,673) npol, maxpol
          call caserr('try again')
        endif
        mpol(npol) = nxtpts
        call drfnval(namei,nval,zat)
        if (ier .ne. 0) ner = ner + 1
        polar(npol) = alf
        vale(npol) = 1.0d0*nval
        npolin = npolin + 1
cxxx  else if (namei(3:) .eq. ' ') then
cxxx    write(namei(3:),674) nxtpts-1,nxcent(nxtpts-1)(ngrnam:)
      endif
cxxx  if (namei(1:2) .eq. 'qq') namei(1:2) = 'qq'
      xpts(1,nxtpts) = cxyz(2)
      xpts(2,nxtpts) = cxyz(3)
      xpts(3,nxtpts) = cxyz(4)
      ztot = ztot + charge
      chrg(nxtpts) = charge
      nxcent(nxtpts) = namei
      radext(nxtpts) = radi
      alfext(nxtpts) = alf
      if (radi .gt. radmax) radmax = radi
      igranl(nxtpts) = igran
caleko
c
c    Drf-cordinates and vdWaals radii are stored to be
c    used in the generation of the Connolly surface.
c
      co(1,nxtpts+nat) = xpts(1,nxtpts)
      co(2,nxtpts+nat) = xpts(2,nxtpts)
      co(3,nxtpts+nat) = xpts(3,nxtpts)
      rtype(nxtpts+nat) = radext(nxtpts)
cxxx  molnum(nxtpts+nat) = ngroup+1
      molnum(nxtpts+nat) = 0
      goto 10
c
300   continue
cxxx  NOW done in rfin, after all input has been processed
caleko
c
c    Gamess-coordinates and vdWaals radii are stored to be
c    used in the generation of the Connolly surface.
c
c     do 121 n=1,nat
c       co(1,n) = cnew(n,1)
c       co(2,n) = cnew(n,2)
c       co(3,n) = cnew(n,3)
c       rtype(n) = radius(anam(n),0,0.0,0.0)
c       molnum(n) = 1
c121  continue
c     natom=nxtpts+nat 
c
cxxx
c     check analysis groups
c
      oreanal = .false.
      do i = 1, mxgran
        mgran(i) = 0
      enddo
      do i = 1, nxtpts
        mgran(igranl(i)) = 1
        if (igranl(i) .gt. ngran) then
          oreanal = .true.
        endif
      enddo
      if (.not. oreanal .and. (mgran(1) .eq. 0))
     1  oreanal = .true.
      if (oreanal) then
        write(iwr,6750)
 6750   format(/,1x,'*****WARNING***** Analysis groups ',
     1  'renumbered!',/)
        ngran = 0
        do i = 1, mxgran
          if (mgran(i) .eq. 1) then
            ngran = ngran + 1
            mgran(i) = ngran
          endif
        enddo
        do i = 1, nxtpts
          igranl(i) = mgran(igranl(i))
        enddo
      endif
cxxx
c  output section
      write(iwr,675) npoints,npolin,nxtpts,npol-ngrp,ngrp
cxxx  write(iwr,675) npoints,npolin,nxtpts,npol-ngrp,ngrp,namb
      ngrpol = ngrp
      if (namb .ne. 0) then
c
c  an abiguous atom is by definition "close"
c  therefore modify the potentials and fields
        modxza = 1
        write(iwr,676)
      endif
cxxx  if (idrfout .ne. 0) then
        write(iwr,677)
        imp = 1
        do 350, ii = 1, nxtpts
c
c  -----  polarizability used in classical dispersion
c
          poldis = 0.0d0
          if (iclintd .eq. 1) then
            if (igrppol .eq. 1) then
              if (nxcent(ii)(:5) .eq. 'group') then
                poldis = polar(imp)
              else if (imemgrp(ii) .ne. 0) then
                poldis = 0.0d0
              else
                if (ieffpol .eq. 1) then
                  poldis = atpol(ii)
                else if (ii .eq. mpol(imp)) then
                  poldis = polar(imp)
                endif
              endif
            else
              if (nxcent(ii)(:5) .eq. 'group') then
                poldis = 0.0d0
              else if (ieffpol .eq. 1) then
                poldis = atpol(ii)
              else if (ieffpol .eq. 0) then
                if (ii .eq. mpol(imp)) then
                  poldis = polar(imp)
                else
                  poldis = atpol(ii)
                endif
              endif
            endif
          endif
c
c
c
          if (ii .eq. mpol(imp)) then
            write(iwr,678) nxcent(ii),(xpts(k,ii),k=1,3),chrg(ii),
cxxx +        radext(ii),polar(imp),atpol(ii),igranl(ii)
     +        radext(ii),polar(imp),poldis,igranl(ii)
            imp = imp + 1
          else
            write(iwr,678) nxcent(ii),(xpts(k,ii),k=1,3),chrg(ii),
cxxx +        radext(ii),0.0d0,atpol(ii),igranl(ii)
     +        radext(ii),0.0d0,poldis,igranl(ii)
          endif
350     continue
cxxx  endif
      if (ner .ne. 0) write(iwr,679) ner
      write(iwr,680) ztot

      if ((idrfout .ge. 1) .and. (nprint .ne. -5)) call drfdst

c
c  skip to here if there is no discrete environment present
360   continue
c
c  check number of analysis groups requested
      if (ngran .gt. mxgran) then
        write(iwr,681) ngran
        call caserr('try again')
      endif
c
c  check if internal polarisabilities are already present
c  for calculation of reverse dispersion or because it is required
      if (irevdis .eq. 1) ialfso = 1
      if (ialfso .ne. 0) then
        if (ioda(1,213) .ne. 0) ialfso = 0
      endif
c
c  check if atomic charges for qm atoms are read (070296)
      if (iradex .eq. 4 .and. .not. (nat .eq. 1 .and.
     +  (anam(1) .eq. 'z' .or. anam(1) .eq. 'e'))) then
        if (.not. okqmch) then
          write(iwr,682)
          call caserr('try again')
        endif
      endif
cahv
      if (nxtpts .eq. 0) nodiscr = 1
cahv
c
c  define internal radii
cxxx  do 400, n = 1, nat
c
c  radii extended with frecer's model (070296)
c       if (anam(n) .eq. 'e' .or. anam(n) .eq. 'z') then
c         radi = 1.0d-5
c       elseif ((iradat .ge. 0 .and. iradat .le. 2) .or.
c    +          (iradat .eq. 4)) then
c         radi = radius(anam(n),iradat,qmcharge(n),afct)
c       else
c         call caserr('illegal option of iradat')
c       endif
c       if (radi .eq. 0.0d0) radi = afact
c       radat(n) = radi
c       alfat(n) = alfa(anam(n),ithole,ier)
c       if (radi .gt. radmax) radmax = radi
c       if (idrfout .eq. 3) write(iwr,*) anam(n), radat(n), alfat(n)
c00   continue
c     if (idrfout .eq. 3) write(iwr,683) radmax
c
c  process classical group pes input if required
      if (nopes .eq. 0) call pesinpc
c
665   format(/'reading external points')
666   format('/going to read from tape# ',i2)
667   format(/'following atoms discarded (no name, charge zero)'/,
     +  a16,3f10.6)
668   format(/i2,' atoms of group ', a14, ' discarded (too far)'/)
6681  format(/i2,' atoms of group ', a14, ' NOT grouped',
     1 ' (too close to QM system: ',f8.4,' cf DSTGRP= ',f8.4,' )'/)
669   format(a)
670   format(/'number of mc group members exceeds maximum of ',i3,
     +  ' atoms for group ',i3)
671   format(/i10,' too many points in rfin: raise maxpnts = ',i10/)
672   format(/'too many qm atoms; increase mxat in hondo or ',
     +  ' decrease number of ambiguous atoms in $external'/)
673   format(/i10,' too many points in rfin: raise maxpol = ',i10/)
674   format(i3,a)
675   format(/i4,' points found on input ',
     +  'of which ',i4,' polarisable'/
     +  '   constructed'/
     +  '   ',i4,' (charged) points'/
     +  '   ',i4,' point polarisabilities and'/
     +  '   ',i4,' group polarisabilities'/
     +  '   '/)
cxxx +  '   ',i4,' ambiguous'/)
676   format(/'rfin changed: modxza = 1'/)
677   format(/,' classical system specification ',/,
     + 'name',t21,'x',t31,'y',t41,'z',t51,'charge',
     +  t58,'   radius alfa (b**3) disalf (b**3) analysis'/)
678   format (1x,a16,4f10.6,3(1x,f8.3,1x),7x,i2)
679   format(/i4,' atoms without polarisability'/)
680   format(/'total external charge ',e20.12/)
681   format(/'number of analysis groups requested = ',i3,
     +  ' is too large: increase mxgran in hondo and drf/dimpar'/)
682   format(/'when iradex = 4 and there are qm atoms present you
     +  need a qmcharge block with name and charge for every atom'/)
c83   format(/'radmax = ',f10.6/)
c
      return
      end
chvd
c     subroutines daopen, daclos, daread, and dawrite have been moved
c     to GAMESS-UK/m4/gamess_hondo.m as they are now shared with the
c     HONDO NMR code. This move was necessary to avoid name clashes
c     between the nmr and drf libraries.
chvd
