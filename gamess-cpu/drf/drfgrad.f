      SUBROUTINE DRFGRAD(xscm)
      IMPLICIT real*8 (A-H,O-Z)
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
      common/junk/de(3,maxat)
      common/hermit/h(45)
      common/wermit/w(45)
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
      real*8 drfde
      common /rfgrad/ drfde(3,mxpts)
c
C
      dimension xscm(*)
C
C     DRF-gradient contributions
C
C-----  BEGIN
C
      IEPS = 0
      CALL DRFGRAD2(IEPS,XSCM)
      CALL EXTGRAD(IEPS,XSCM)
      call blkegr(nxtpts,drfde)
      RETURN
      END
      SUBROUTINE DRFGRAD2(IEPS,xscm)
      IMPLICIT real*8 (A-H,O-Z)
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
      common/junk/de(3,maxat)
      common/hermit/h(45)
      common/wermit/w(45)
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
      real*8 drfde
      common /rfgrad/ drfde(3,mxpts)
c
C
      dimension xscm(*)
C
C     DRF-gradient contributions
C
      DIMENSION DRFD(3,maxat)
C
      DATA ZERO /0.D00/
C
C-----  BEGIN
C
C
C     -----  initialize DRF-derivative contributions
C
      CALL CLEAR(DRFD,3*NEXP)
C
C-----  Set DA-record numbers for reading RF contribution and matrices
C
      ILCL = 21
C
      IF (IEPS .EQ. 0) THEN
        ILWT = 26
        ILVR = 28
        ILLUR = 43
        ILINDX = 44
        ILIJ = 56
        ILOMG = 50
      ELSE
        ILWT = 30
        ILVR = 32
        ILLUR = 47
        ILINDX = 48
        ILIJ = 57
        ILOMG = 51
      ENDIF
C
      IF (MCUPDT) THEN
        ILWT = ILWT + 1
        ILVR = ILVR + 1
        ILCL = 91
      ENDIF
cahv
c     NOT YET IMPLEMENTED! By treating QM fragment classically
c     a faster gradient may be achieved - accuracy???
C
C-----  Calculate Mulliken charges and dipoles and Dipole Preserving
C       charges if either is required for the analysis
C
c     IF (((IFLDIN .LE. 2) .OR. (IFLDOUT .LE. 2))
c    1    .AND. (IEPS .EQ. 0)) THEN
C 1-----
c       IF (NODPE .EQ. 1) THEN
C   2-----
c         ILMC = 121
c         ILMD = 122
c         ILDC = 124
c         CALL DPPOP(xscm)
c       ELSE
c         ILMC = 126
c         ILMD = 127
c         ILDC = 129
c         CALL DPPOPE(xscm)
C   2-----
c       ENDIF
C 1-----
c     ENDIF
C
C     -------  nuclear DRF contributions
C
      CALL CLEAR(DRFDE,3*NEXP)
c
      ixexp = igmem_alloc(3*nexp)
C
C-----  Read expansion centres
C
      CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXEXP),3*NEXP,1)
c
      IF (FIELD(5:) .NE. ' ') THEN
C 1-----
        ixomg = igmem_alloc(nomga)
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXOMG),NWTC*NWTC,ILOMG)
      ELSE
        ixomg = igmem_alloc(1)
      ENDIF
c
c-----  gradient contribution 
c     due to nuclei
c
      CALL DRFGRADNUC(DRFDE,CZAN,xscm(ixexp),NAT,xscm(ixomg),NWTC)
      if (idrfout .ge. 3) 
     1  call hatout(drfde,3,nexp,2,'nucgrad')
C
C     ------- add to gradient
C
      DO I = 1 , NAT
        DO J = 1 , 3
          DRFD(J,I) = DRFD(J,I) + DRFDE(J,I)
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(drfd,3,nexp,2,'gradnuc')
c
c-----  start allocating necessary memory
c
      ixie = igmem_alloc(nchd)
      ixb = igmem_alloc(nchd)
      ixd = igmem_alloc(nchd)
      ixdb = igmem_alloc(nchd)
      ixdh = igmem_alloc(nchd)
C
C-----  Read density
C
      CALL DENSRD(XSCM(IXD),XSCM(IXDB),XSCM(IXDH))
      call gmem_free(ixdh)
C
C-----  Read overlap and first moment integrals
C
      ixol = igmem_alloc(nchd)
      ixdx = igmem_alloc(nchd)
      ixdy = igmem_alloc(nchd)
      ixdz = igmem_alloc(nchd)
      CALL DAREAD(IDAFh,IODA,XSCM(IXOL),NCHD,12)
      CALL DAREAD(IDAFh,IODA,XSCM(IXDX),NCHD,53)
      CALL DAREAD(IDAFh,IODA,XSCM(IXDY),NCHD,54)
      CALL DAREAD(IDAFh,IODA,XSCM(IXDZ),NCHD,55)
      if (field(5:) .ne. ' ') then
        ixrxx = igmem_alloc(nchd)
        ixrxy = igmem_alloc(nchd)
        ixryx = igmem_alloc(nchd)
        ixrxz = igmem_alloc(nchd)
        ixrzx = igmem_alloc(nchd)
        ixryy = igmem_alloc(nchd)
        ixryz = igmem_alloc(nchd)
        ixrzy = igmem_alloc(nchd)
        ixrzz = igmem_alloc(nchd)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrxx),NCHD,61)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXryy),NCHD,62)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrzz),NCHD,63)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrxy),NCHD,64)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrxz),NCHD,65)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXryz),NCHD,66)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXryx),NCHD,64)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrzx),NCHD,65)
        CALL DAREAD(IDAFdrf,IODAdrf,XSCM(IXrzy),NCHD,66)
      else
        ixrxx = igmem_alloc(1)
        ixrxy = igmem_alloc(1)
        ixryx = igmem_alloc(1)
        ixrxz = igmem_alloc(1)
        ixrzx = igmem_alloc(1)
        ixryy = igmem_alloc(1)
        ixryz = igmem_alloc(1)
        ixrzy = igmem_alloc(1)
        ixrzz = igmem_alloc(1)
      endif
C
      IF (FIELD(5:) .NE. ' ') THEN
C 1-----
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXOMG),NWTC*NWTC,ILOMG)
        ixwt = igmem_alloc(nwtr*nwtc)
        ixvr = igmem_alloc(nwtr*nwtc)
        ixr = igmem_alloc(ndim*ndim)
        ixi = igmem_alloc(ndim)
        ixind = igmem_alloc(ndim)
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXWT),NWTR*NWTC,ILWT)
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXVR),NWTR*NWTC,ILVR)
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXR),NDIM*NDIM,ILLUR)
        CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXI),NDIM,ILINDX)
        CALL DAREAD(IDAFIND,IODAIND,XSCM(IXIND),NDIM,INDEX)
C  -----
        IF (FIELD(5:) .EQ. 'scf') THEN
          CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXB),NCHD,ILIJ)
        ELSE
          call intinit(xscm(ixb),nchd,1)
        ENDIF
C 1-----
      ELSE
C 1-----
        ixwt = igmem_alloc(1)
        ixvr = igmem_alloc(1)
        ixr = igmem_alloc(1)
        ixi = igmem_alloc(1)
        ixind = igmem_alloc(1)
C 1-----
      ENDIF
C
C-----  Read assignment vector
C
      CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXIE),NCHD,2)
      CALL CLEAR(DRFDE,3*NEXP)
C
      if (field(5:) .ne. ' ') then
        ixdwt = igmem_alloc(nwtr*nwtc)
        ixdvr = igmem_alloc(nwtr*nwtc)
        ixdipix = igmem_alloc(nwtc*ndim)
      else
        ixdwt = igmem_alloc(1)
        ixdvr = igmem_alloc(1)
        ixdipix = igmem_alloc(1)
      endif
      if (nbem .ne. 0) then
        ixsur = igmem_alloc(3*nbem)
        ixnor = igmem_alloc(3*nbem)
        ixar = igmem_alloc(nbem)
      else
        ixsur = igmem_alloc(3)
        ixnor = igmem_alloc(3)
        ixar = igmem_alloc(1)
      endif
C
C     -------  gradient due to DRF integrals w.r.t. expansion centers
C
      CALL DRFGRADX(DRFDE,XSCM(IXD),XSCM(IXDB),XSCM(IXV),
     1     XSCM(IXOL),XSCM(IXDX),XSCM(IXDY),XSCM(IXDZ),
     2     XSCM(IXRXX),XSCM(IXRYY),XSCM(IXRZZ),
     3     XSCM(IXRXY),XSCM(IXRXZ),XSCM(IXRYZ),
     4     XSCM(IXSUR),XSCM(IXNOR),XSCM(IXAR),XSCM(IXEXP),
     4     XSCM(IXWT),XSCM(IXVR),XSCM(IXDWT),XSCM(IXDVR),
     5     XSCM(IXOMG),XSCM(IXR),XSCM(IXIND),XSCM(IXDIPix),
     6     XSCM(IXIE),XSCM(IXB),XSCM(IXI),
     7     IEPS,INEQ)
C
      if (idrfout .ge. 3) 
     1  call hatout(drfde,3,nexp,2,'gradx')
C
C     ------- add to gradient
C
      DO I = 1 , NEXP
        DO J = 1 , 3
          DRFD(J,I) = DRFD(J,I) + DRFDE(J,I)
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(drfd,3,nexp,2,'xgrad')
C
      call gmem_free(ixar)
      call gmem_free(ixnor)
      call gmem_free(ixsur)
      call gmem_free(ixdipix)
      call gmem_free(ixdvr)
      call gmem_free(ixdwt)
      call gmem_free(ixind)
      call gmem_free(ixi)
      call gmem_free(ixr)
      call gmem_free(ixvr)
      call gmem_free(ixwt)
c
C     -------  gradient due to DRF integrals w.r.t. nuclei
C
      CALL CLEAR(DRFDE,3*NEXP)
C
	   if (field(5:) .ne. ' ')
     1  CALL DAREAD(IDAFDRF,IODADRF,XSCM(IXOMG),NWTC*NWTC,ILOMG)
      ixdsx = igmem_alloc(nchd)
      ixdsy = igmem_alloc(nchd)
      ixdsz = igmem_alloc(nchd)
      ixdxx = igmem_alloc(nchd)
      ixdyx = igmem_alloc(nchd)
      ixdzx = igmem_alloc(nchd)
      ixdxy = igmem_alloc(nchd)
      ixdyy = igmem_alloc(nchd)
      ixdzy = igmem_alloc(nchd)
      ixdxz = igmem_alloc(nchd)
      ixdyz = igmem_alloc(nchd)
      ixdzz = igmem_alloc(nchd)
      ixdxxx = igmem_alloc(nchd)
      ixdxxy = igmem_alloc(nchd)
      ixdxxz = igmem_alloc(nchd)
      ixdyyx = igmem_alloc(nchd)
      ixdyyy = igmem_alloc(nchd)
      ixdyyz = igmem_alloc(nchd)
      ixdzzx = igmem_alloc(nchd)
      ixdzzy = igmem_alloc(nchd)
      ixdzzz = igmem_alloc(nchd)
      ixdxyx = igmem_alloc(nchd)
      ixdxyy = igmem_alloc(nchd)
      ixdxyz = igmem_alloc(nchd)
      ixdxzx = igmem_alloc(nchd)
      ixdxzy = igmem_alloc(nchd)
      ixdxzz = igmem_alloc(nchd)
      ixdyzx = igmem_alloc(nchd)
      ixdyzy = igmem_alloc(nchd)
      ixdyzz = igmem_alloc(nchd)
      if (field(5:) .ne. ' ') then
        ixv = igmem_alloc(num*num)
      else
        ixv = igmem_alloc(1)
      endif
      do iatm = 1, nat 
         call clear(xscm(ixdsx),nchd)
         call clear(xscm(ixdsy),nchd)
         call clear(xscm(ixdsz),nchd)
c
c  -----  derivatives of overlap matrix to
c         atom iatm
c
         call drfsint(iatm,
     1      xscm(ixdsx),xscm(ixdsy),xscm(ixdsz))
	 if (idrfout .ge. 3) then
           call hatout(xscm(ixdsx),num,num,3,'DSX')
	   call hatout(xscm(ixdsy),num,num,3,'DSY')
	   call hatout(xscm(ixdsz),num,num,3,'DSZ')
	 endif
         call clear(xscm(ixdxx),nchd)
         call clear(xscm(ixdyx),nchd)
         call clear(xscm(ixdzx),nchd)
         call clear(xscm(ixdxy),nchd)
         call clear(xscm(ixdyy),nchd)
         call clear(xscm(ixdzy),nchd)
         call clear(xscm(ixdxz),nchd)
         call clear(xscm(ixdyz),nchd)
         call clear(xscm(ixdzz),nchd)
c
c  -----  derivatives of dipole matrix to
c         atom iatm
c
         call drfdint(iatm,xscm(ixdxx),xscm(ixdxy),xscm(ixdxz),
     1    xscm(ixdyx),xscm(ixdyy),xscm(ixdyz),
     1    xscm(ixdzx),xscm(ixdzy),xscm(ixdzz))
	 if (idrfout .ge. 3) then
           call hatout(xscm(ixdxx),num,num,3,'DSXX')
           call hatout(xscm(ixdyx),num,num,3,'DSXY')
           call hatout(xscm(ixdzx),num,num,3,'DSXZ')
           call hatout(xscm(ixdxy),num,num,3,'DSYX')
           call hatout(xscm(ixdyy),num,num,3,'DSYY')
           call hatout(xscm(ixdzy),num,num,3,'DSYZ')
           call hatout(xscm(ixdxz),num,num,3,'DSZX')
           call hatout(xscm(ixdyz),num,num,3,'DSZY')
           call hatout(xscm(ixdzz),num,num,3,'DSZZ')
         endif
         if ((field(5:) .ne. ' ') .and. (gamdrf .ne. 0.0)) then
           call clear(xscm(ixdxxx),nchd)
           call clear(xscm(ixdxxy),nchd)
           call clear(xscm(ixdxxz),nchd)
           call clear(xscm(ixdyyx),nchd)
           call clear(xscm(ixdyyy),nchd)
           call clear(xscm(ixdyyz),nchd)
           call clear(xscm(ixdzzx),nchd)
           call clear(xscm(ixdzzy),nchd)
           call clear(xscm(ixdzzz),nchd)
           call clear(xscm(ixdxyx),nchd)
           call clear(xscm(ixdxyy),nchd)
           call clear(xscm(ixdxyz),nchd)
           call clear(xscm(ixdxzx),nchd)
           call clear(xscm(ixdxzy),nchd)
           call clear(xscm(ixdxzz),nchd)
           call clear(xscm(ixdyzx),nchd)
           call clear(xscm(ixdyzy),nchd)
           call clear(xscm(ixdyzz),nchd)
c
c    -----  derivatives of second moments
c           wrt atom iatm
c
         call qmdrfder(iatm,xscm(ixdxxx),xscm(ixdxxy),xscm(ixdxxz),
     1   xscm(ixdyyx),xscm(ixdyyy),xscm(ixdyyz),
     1   xscm(ixdzzx),xscm(ixdzzy),xscm(ixdzzz),
     1   xscm(ixdxyx),xscm(ixdxyy),xscm(ixdxyz),
     1   xscm(ixdxzx),xscm(ixdxzy),xscm(ixdxzz),
     1   xscm(ixdyzx),xscm(ixdyzy),xscm(ixdyzz))
           if (idrfout .ge. 3) then
             call hatout(xscm(ixdxxx),num,num,3,'DSXXX')
             call hatout(xscm(ixdxxy),num,num,3,'DSXXY')
             call hatout(xscm(ixdxxz),num,num,3,'DSXXZ')
             call hatout(xscm(ixdyyx),num,num,3,'DSYYX')
             call hatout(xscm(ixdyyy),num,num,3,'DSYYY')
             call hatout(xscm(ixdyyz),num,num,3,'DSYYZ')
             call hatout(xscm(ixdzzx),num,num,3,'DSZZX')
             call hatout(xscm(ixdzzy),num,num,3,'DSZZY')
             call hatout(xscm(ixdzzz),num,num,3,'DSZZZ')
             call hatout(xscm(ixdxyx),num,num,3,'DSXYX')
             call hatout(xscm(ixdxyy),num,num,3,'DSXYY')
             call hatout(xscm(ixdxyz),num,num,3,'DSXYZ')
             call hatout(xscm(ixdxzx),num,num,3,'DSXZX')
             call hatout(xscm(ixdxzy),num,num,3,'DSXZY')
             call hatout(xscm(ixdxzz),num,num,3,'DSXZZ')
             call hatout(xscm(ixdyzx),num,num,3,'DSYZX')
             call hatout(xscm(ixdyzy),num,num,3,'DSYZY')
             call hatout(xscm(ixdyzz),num,num,3,'DSYZZ')
           endif
         endif
	 do kxyz = 1, 3
	   if (kxyz .eq. 1) then
             ixdsd = ixdsx
             ixdxd = ixdxx
             ixdyd = ixdyx
             ixdzd = ixdzx
             ixdxxd = ixdxxx 
             ixdyyd = ixdyyx 
             ixdzzd = ixdzzx 
             ixdxyd = ixdxyx 
             ixdxzd = ixdxzx 
             ixdyzd = ixdyzx 
	   else if (kxyz .eq. 2) then
             ixdsd = ixdsy
             ixdxd = ixdxy
             ixdyd = ixdyy
             ixdzd = ixdzy
             ixdxxd = ixdxxy 
             ixdyyd = ixdyyy 
             ixdzzd = ixdzzy 
             ixdxyd = ixdxyy 
             ixdxzd = ixdxzy 
             ixdyzd = ixdyzy 
	   else
             ixdsd = ixdsz
             ixdxd = ixdxz
             ixdyd = ixdyz
             ixdzd = ixdzz
             ixdxxd = ixdxxz 
             ixdyyd = ixdyyz 
             ixdzzd = ixdzzz 
             ixdxyd = ixdxyz 
             ixdxzd = ixdxzz 
             ixdyzd = ixdyzz 
	   endif
c
c-----  gradient contributions from 
c       derivatives of overlap, first (and second) moments
c
           call momgrad(iatm,kxyz,drfde,xscm(ixd),xscm(ixdb),
     1     xscm(ixol),xscm(ixdx),xscm(ixdy),xscm(ixdz),
     1     xscm(ixdsd),xscm(ixdxd),xscm(ixdyd),xscm(ixdzd),
     1     xscm(ixomg),xscm(ixie),xscm(ixb))
	   if (field(5:) .ne. ' ') then
             call drfmom2(iatm,kxyz,drfde,xscm(ixd),xscm(ixdb),
     1     xscm(ixv),
     1     xscm(ixol),xscm(ixdx),xscm(ixdy),xscm(ixdz),
     2     XSCM(ixdxxd),XSCM(ixdyyd),XSCM(ixdzzd),
     3     XSCM(ixdxyd),XSCM(ixdxzd),XSCM(ixdyzd),
     4     xscm(ixdxyd),xscm(ixdxzd),xscm(ixdyzd),
     1     xscm(ixdsd),xscm(ixdxd),xscm(ixdyd),xscm(ixdzd),
     1     xscm(ixomg),xscm(ixie),xscm(ixb))
	   endif
	 enddo
      enddo
C
      if (idrfout .ge. 3) 
     1  call hatout(drfde,3,nexp,2,'gradx')
C
C     ------- add to gradient
C
      DO I = 1 , NEXP
        DO J = 1 , 3
          DRFD(J,I) = DRFD(J,I) + DRFDE(J,I)
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(drfd,3,nexp,2,'gradnuc')
C
C
C-----  Distribute force on centre of nuclear 
C       charge 
C
      IF (nexp .gt. nat) THEN
        CALL DISFOR(DRFD,CZAN,NEXP,NAT)
      ENDIF
      if (idrfout .ge. 3) 
     1  call hatout(drfd,3,nat,2,'drfgtot')
C
C-----  Add RF gradient to QM gradient
C
      call gmem_free(ixv)
      call gmem_free(ixdyzz)
      call gmem_free(ixdyzy)
      call gmem_free(ixdyzx)
      call gmem_free(ixdxzz)
      call gmem_free(ixdxzy)
      call gmem_free(ixdxzx)
      call gmem_free(ixdxyz)
      call gmem_free(ixdxyy)
      call gmem_free(ixdxyx)
      call gmem_free(ixdzzz)
      call gmem_free(ixdzzy)
      call gmem_free(ixdzzx)
      call gmem_free(ixdyyz)
      call gmem_free(ixdyyy)
      call gmem_free(ixdyyx)
      call gmem_free(ixdxxz)
      call gmem_free(ixdxxy)
      call gmem_free(ixdxxx)
      call gmem_free(ixdzz)
      call gmem_free(ixdyz)
      call gmem_free(ixdxz)
      call gmem_free(ixdzy)
      call gmem_free(ixdyy)
      call gmem_free(ixdxy)
      call gmem_free(ixdzx)
      call gmem_free(ixdyx)
      call gmem_free(ixdxx)
      call gmem_free(ixdsz)
      call gmem_free(ixdsy)
      call gmem_free(ixdsx)
c
      call gmem_free(ixrzz)
      call gmem_free(ixrzy)
      call gmem_free(ixryz)
      call gmem_free(ixryy)
      call gmem_free(ixrzx)
      call gmem_free(ixrxz)
      call gmem_free(ixryx)
      call gmem_free(ixrxy)
      call gmem_free(ixrxx)
      call gmem_free(ixdz)
      call gmem_free(ixdy)
      call gmem_free(ixdx)
      call gmem_free(ixol)
      call gmem_free(ixdb)
      call gmem_free(ixd)
      call gmem_free(ixb)
      call gmem_free(ixie)
      call gmem_free(ixomg)
      call gmem_free(ixexp)
c
      if (idrfout .ge. 3) 
     1  call hatout(de,3,nat,2,'ghon')
      DO I = 1 , NAT
        DO J = 1 , 3
          DE(J,I) = DE(J,I) + DRFD(J,I)
        ENDDO
      ENDDO
C
      if (idrfout .ge. 3) 
     1  call hatout(DE,3,nat,2,'gradtot')
c
      RETURN
      END
      subroutine intinit(iarray,n,nvalue)
      integer iaraay, n, nvalue
      dimension iarray(n)
      do i = 1, n
        iarray(i) = nvalue
      enddo
      return
      end
      SUBROUTINE DRFGRADNUC(drfde,zan,xexp,nat,omega,nomg)
C
C     -----  DRF gradient: nuclear contributions
C     
      IMPLICIT real*8 (A-H,O-Z)
C
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
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
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
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      dimension drfde(3,nat), xexp(3,nat)
      dimension zan(nat)
      dimension omega(nomg,nomg)
C
      logical odebug
C
      DATA ZERO, PT5, ONE /0.D00, 5.D-01, 1.D00/
C
C-----  BEGIN
C
      odebug = .false.
C
      nzfp = nexp*4
C     
      DO N = 1, NAT
C 1-----
        JF = (N-1)*4
C     
C  -----  external charges
C     
        DNUEQX = ZERO
        DNUEQY = ZERO
        DNUEQZ = ZERO
        DO IGR = 1,NGRAN
          DNUEQX =ZFA(JF+1,IGR)* ZAN(N) + DNUEQX
          DNUEQY =ZFA(JF+2,IGR)* ZAN(N) + DNUEQY
          DNUEQZ =ZFA(JF+3,IGR)* ZAN(N) + DNUEQZ
        END DO
        DRFDE(1,N)=DRFDE(1,N) + DNUEQX
        DRFDE(2,N)=DRFDE(2,N) + DNUEQY
        DRFDE(3,N)=DRFDE(3,N) + DNUEQZ
C
        IF (FIELD(5:) .NE. ' ') THEN
C   2-----
C    -----  external dipoles
C
          DNUEDX = ZERO
          DNUEDY = ZERO
          DNUEDZ = ZERO
          DO IGR = 1, NGRAN     
            DNUEDX = (OMEGA(NZFP+IGR,JF+1) 
     1              + OMEGA(JF+1,NZFP+IGR))* ZAN(N)
            DNUEDY = (OMEGA(NZFP+IGR,JF+2) 
     1              + OMEGA(JF+2,NZFP+IGR))* ZAN(N)
            DNUEDZ = (OMEGA(NZFP+IGR,JF+3) 
     1              + OMEGA(JF+3,NZFP+IGR))* ZAN(N)
          END DO
C     
          DRFDE(1,N) = DRFDE(1,N) + pt5*DNUEDX
          DRFDE(2,N) = DRFDE(2,N) + pt5*DNUEDY
          DRFDE(3,N) = DRFDE(3,N) + pt5*DNUEDZ
          if (odebug) CALL hATOUT(drfde(1,n),3,1,2,'drfde')
C     
C    -----  screening nuclear repulsion
C     
          DO M = 1, NAT
C     3-----
            q(1) = XEXP(1,M)
            q(2) = XEXP(2,M)
            q(3) = XEXP(3,M)
            q(4) = ONE
            CALL DRFOAB(N,M,NWTC,OMEGA)
            CALL MATVEC(OMGAB,q,P,4,.FALSE.)
            DSNREPX = ZAN(N)*ZAN(M)*P(1)
            DSNREPY = ZAN(N)*ZAN(M)*P(2)
            DSNREPZ = ZAN(N)*ZAN(M)*P(3)
            CALL DRFOAB(M,N,NWTC,OMEGA)
            CALL MATVEC(OMGAB,q,P,4,.TRUE.)
            DSNREPX = ZAN(N)*ZAN(M)*P(1) + DSNREPX
            DSNREPY = ZAN(N)*ZAN(M)*P(2) + DSNREPY
            DSNREPZ = ZAN(N)*ZAN(M)*P(3) + DSNREPZ
C     
            DRFDE(1,N) = DRFDE(1,N) + PT5*DSNREPX
            DRFDE(2,N) = DRFDE(2,N) + PT5*DSNREPY
            DRFDE(3,N) = DRFDE(3,N) + PT5*DSNREPZ
            if (odebug) CALL hATOUT(drfde(1,n),3,1,2,'drfde')
C     3-----
          END DO
C   2-----
        ENDIF
C     
C     ------  model repulsion
C     
        DO IGR = 1, NGRAN
          DRFDE(1,N) = DRFDE(1,N) + DPSEU(1,N,IGR)
          DRFDE(2,N) = DRFDE(2,N) + DPSEU(2,N,IGR)
          DRFDE(3,N) = DRFDE(3,N) + DPSEU(3,N,IGR)
        ENDDO
        if (odebug) CALL hATOUT(drfde(1,n),3,1,2,'drfde')
C 1-----        
      END DO
      RETURN
      END
      SUBROUTINE DRFGRADX(DRFDE,DA,DB,V,SS,DX,DY,DZ,
     1 RXX,RYY,RZZ,RXY,RXZ,RYZ,XSURF,XNORM,AREA,
     2 XEXP,WT,VR,DWT,DVR,OMEGA,RELAY,DIPIND,DIPIX,
     3 IEXPC,IJBIT,INDX,
     4 IEPS,INEQ)
C------
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION DRFDE(3,128)
      DIMENSION DA(1),SS(1),DB(1),DX(1),DY(1),DZ(1),
     1 RXX(1),RYY(1),RZZ(1),RXY(1),RXZ(1),RYZ(1),
     2 DOMEGA(1),IEXPC(1),IJBIT(1)
      DIMENSION V(NUM,NUM)
      DIMENSION XEXP(3,NEXP)
      DIMENSION RELAY(NDIM,NDIM)
      DIMENSION DIPIND(1)
      DIMENSION INDX(NDIM)
      DIMENSION WT(NWTR,NWTC), VR(NWTR,NWTC)
      DIMENSION DWT(NWTR,NWTC), DVR(NWTR,NWTC)
      DIMENSION OMEGA(NWTC,NWTC)
      DIMENSION XSURF(3,NDIMB),XNORM(3,NDIMB),AREA(NDIMB)
      DIMENSION DIPIX(NWTC,NDIM)
C
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
      logical odebug
C
C------ DRF gradients w.r.t. expansion centers
C
      odebug = .false.
C
C------ loop over x,y,z components
C
      DO KXYZ = 1, 3
C 1-----
C  -----  derivatives of sources AND receptors
C
        IF (FIELD(5:) .NE. ' ') THEN
c   2-----
          CALL DWTVRCAL(KXYZ,XEXP,DWT,DVR,
     1       XSURF,XNORM,AREA,IEPS,INEQ)
          if (odebug) then
            call hatout(dwt,nwtr,nwtc,2,'DWT')
            call hatout(dvr,nwtr,nwtc,2,'DVR')
          endif
C
C    -----  Construct formal interactions with 
C           derivative of source
C
          CALL DRFdOMGA(IEPS,RELAY,DWT,VR,
     1    OMEGA,INDX,DIPIX)
          if (odebug) call hatout(omega,nwtc,nwtc,2,'DOM')
C
C    -----  Electrons are source
C
          CALL DRFGRAX2(KXYZ,DRFDE,XEXP,DA,DB,V,SS,DX,DY,DZ,
     1         RXX,RYY,RZZ,RXY,RXZ,RYZ,OMEGA,IEXPC,IJBIT)
C
C
C    -----  Construct formal interactions with 
C           derivative of recipient
C
          CALL DRFdOMGA(IEPS,RELAY,WT,DVR,
     1         OMEGA,INDX,DIPIX)
c
C    -----  Electrons are recipients
c
          CALL DRFGRAX2B(KXYZ,DRFDE,XEXP,DA,DB,V,SS,DX,DY,DZ,
     1         RXX,RYY,RZZ,RXY,RXZ,RYZ,OMEGA,IEXPC,IJBIT)
C
C    -----  Electrons are source
C
C    -----  Induced dipoles coming back on electrons
C
        ENDIF
C
C  -----  Calculate electronic contribution to 
C         derivative from interaction 
C         with external charges
C
        CALL DRFGRAX(KXYZ,DRFDE,DA,SS,DX,DY,DZ,
     1       IEXPC,IJBIT)
C 1-----
      END DO
      RETURN
      END	
      SUBROUTINE DWTVRCAL(KXYZ,XEXP,WT,VR,
     1           XSURF,XNORM,AREA,IEPS,INEQ)
C------
C       Calculation of derivative of
C       (expanded) source and reaction fields
C       of/on (formal) QM particles (nuclei, electrons)
C       or representations (DP charges, Mulliken charges
C       and dipoles) thereof, at polarizabilities,
C       boundary elements and external charges.
C
C       --------  P.Th. van Duijnen, IBM-KINGSTON 1985, and
C                 Groningen, Dec. 1991.
C
C       --------  Adapted from WTVRCAL
C                 KXYZ gives the index of the 
C                 cartesian derivative
C                 1=X; 2=Y; 3=Z
C                 A.H. de Vries, Daresbury Lab 1997, 
C                 and Groningen, 1997
C------
      IMPLICIT real*8 (A-H,O-Z)
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
C
C-----  Dummy arrays
C
      DIMENSION XSURF(3,NBEM),XNORM(3,NBEM),AREA(NBEM)
      DIMENSION XEXP(3,NEXP),WT(NWTR,NWTC),
     1          VR(NWTR,NWTC)
C
C-----  Common blocks
C
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
C
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
      integer neq2eps, momsav
      common /neqpar1/ neq2eps,momsav
c
      real*8 kapneq1, kapneq2, epsneq1, epsneq2
      common /neqpar2/ kapneq1,kapneq2,epsneq1,epsneq2
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
C
C-----  Local variables
C
      LOGICAL KAPNOZ
C
      CHARACTER*16 NAMJ
C
      DIMENSION P(3), q(3)
      DIMENSION PQ(3),W(3),B(3,3),BQ(3),DELKIQ(3),DELLIQ(3)
      DIMENSION QI(3), T(3,3), TQ(3), TNI(3)
      real*8    FQINI, EXPKD, FKAPONE, DIST, DMIND1, DMIN3
      real*8    XKIQ, XLIQ
      DIMENSION O(3,3)
C
      CHARACTER*8 ERRMSG(3)
C
      DATA ERRMSG/'PROGRAM','STOP IN','-DWTVRCA'/
      DATA THIRD/.33333333333333333333333333333D00/
      DATA SIXTH/.16666666666666666666666666667D00/
      DATA TWO,THREE,FOUR,TWELVE/2.0D00,3.0D00,4.0D00,12.0D00/
      DATA ZERO,PT5,PT75 /0.0D00,0.5D00,0.75D00/
      DATA ONE,ONEPT5 /1.0D00,1.5D00/
C
C-----  BEGIN
C
      IF (IBEM .NE. 0) THEN
C 1-----
        IF (INEQ .EQ. 0) THEN
C   2-----
          IF (IEPS .EQ. 0) THEN
C     3-----
            EPS = EPS1
            KAPPA = KAPPA1
          ELSE
            EPS = EPS2
            KAPPA = KAPPA2
C     3-----
          ENDIF
C   2-----
        ELSE
C   2-----
          IF (IEPS .EQ. 0) THEN
C     3-----
            EPS = EPSNEQ1
            KAPPA = KAPNEQ1
          ELSE
            EPS = EPSNEQ2
            KAPPA = KAPNEQ2
C     3-----
          ENDIF
C   2-----
        ENDIF
C
C  -----  Set logical  KAPNOZ
C         for non-zero ionic strength
C
        KAPNOZ = KAPPA .NE. ZERO
C
C  -----  Initialize some important factors
C
        EXPKD = ONE
        FKAPONE = ONE
        PI=FOUR*ATAN(ONE)
        EPSFACT= ONE/(TWO*PI*(ONE+EPS))
        KAPPAS = KAPPA**2
C 1-----
      ENDIF
C
      IF (MODXZA .EQ. 0) ITHOLE = 0
      ITHOL = ITHOLE
C
C-----  zfp and vrp must be cleared,
c       (cf wtvrcal)
C       since their values do not
C       depend on the co-ordinates of 
C       the expansion centres
C
      CALL CLEAR(WT(1,1),(NEXP4+NGRAN+1)*NDIM)
      CALL CLEAR(VR(1,1),(NEXP4+NGRAN+1)*NDIM)
C
C-----  Loop over the expansion centra
C
      DO J = 1, NEXP
C 1-----
        IF (J .LE. NAT) THEN
C   2-----
C    -----  Nucleus at expansion centre
C
C    -----  Skip ambiguous atoms
C
          IF (NAMBPT(J) .NE. 0) GOTO 500
C
C    -----  Nuclear charge -ZA-
C
          ZA = CZAN(J)
          NAMJ = ANAM(J)
C
C    -----  Polarisability of atom corresponding to expansion centre
C
c         ALFJ = ALFA(NAMJ,0,IER)
          alfj = alfat(j)
C   2-----
        ELSE
C   2-----
C    ----- Non-nulcei are given a polarizability 1.
C
          ALFJ = ONE
C   2-----
        ENDIF
C
C  -----  Pointer to arrays w.r.t. expansion centre -J-
C
        JF = (J-1)*4
C
C  -----  Position vector -J- into -qq(ivoff+)-
C
        q(1) = XEXP(1,J)
        q(2) = XEXP(2,J)
        q(3) = XEXP(3,J)
C
C  -----  Loop over polarizable points
C
        DO 300, II = 1, NPOL
C   2-----
          NP = MPOL(II)
C
C    -----  Skip ambiguous points. If polarisability is required for
C           these points w.r.t. the QM system, add basis functions
C           to represent the polarisability
C
          IF (NCUTPT(NP) .NE. 0) GOTO 300
C
C    -----  Position vector into -P-
C
          P(1) = XPTS(1,NP)
          P(2) = XPTS(2,NP)
          P(3) = XPTS(3,NP)
C
C    -----  Pointer in field arrays w.r.t. polarizable point
C
          IF = (II-1)*3
C
C    -----  Calculate distance vector between -qq(ivoff+)- and -P-
C    -----  PQ contains (p-q)
C
          CALL DISTAB(q,P,PQ,DIST)
C
C    -----  Skip very close polarisabilities (as may occur for
C           ambiguous atoms)
C
          IF (DIST .LE. 1.0D-03) GOTO 300
C
          DMIND1 = ONE/DIST
          DMIN3 = DMIND1*(DMIND1**2)
          FACTP = ONE
          FACTF = ONE
          V = ONE
c
          alfi = alfext(np)
          if (alfi .eq. zero) alfi = one
C
C    -----  Account for penetration effects (optional)
C
          IF (MODXZA .NE. 0) THEN
C     3-----
            S = (ALFJ*ALFi)**SIXTH
            V = DIST/S
            IF (ITHOLE .EQ. 1) THEN
C       4-----
              IF (V .LE. AFACT) THEN
                AV = V/AFACT
                FACTP = AV**4 - TWO*AV**3 + TWO*AV
                FACTF = FOUR*AV**3 - THREE*AV**4
              ENDIF
C       4-----
            ELSE IF (ITHOLE .EQ. 2) THEN
C       4-----
              AU = AFACT*V
              FACTP = (ONE - (PT5*AU + ONE)*EXP(-AU))
c              FACTP = (ONE - (PT5*AU - ONE)*EXP(-AU))
              FACTF = (ONE - (PT5*AU**2 + AU + ONE)*EXP(-AU))
C       4-----
            ENDIF
C     3-----
          ENDIF
          FACTP = FACTP*DMIND1
          FACTF = FACTF*DMIN3
C
C    -----  B(3,3) = t(p;q) : field gradient of charge in -P- at -qq(ivoff+)-
C           Note that the interaction may be scaled
C
          CALL DRFTPQ(PQ,DMIND1,DMIN3,ITHOL,AFACT,V,B)
C
C    -----  O(3,3) = o(KXYZ;p;q) : KXYZ derivative of
C           field gradient of charge in -P- at -qq(ivoff+)-
C           o = -d/drj t(i;j) = d/dri t(i;j)
C           Note that the interaction may be scaled
C
          CALL DRFOPQ(PQ,DMIND1,O,KXYZ,ITHOL,AFACT,V)
C
C    -----  Calculate -W- = o(KXYZ;p;q) . q , this is part of the
C           Taylor expansion of the field
C
          CALL MATVEC(O,q,W,3,.FALSE.)
C
C    -----  -WT- and -VR- matrices for expanded field and potential
C           of charge (distribution) in expansion centra
C           at the polarizable points (-WT-) and vice versa (-VR-)
C
C           Depending on the type of sources (charges, dipoles,
C           charge distributions), specified in DRFINP by IFLDIN
C           and "recipients", specified by IFLDOUT, the matrices
C           are constructed partly or completely.
C
          DO 290, K = 1, 3
C     3-----
C      -----  derivative of the distance vector
C
            PQK = B(KXYZ,K)
            DO 280, L = 1, 3
C       4-----
C        -----  Copy -O- into -WT-  and/or -VR-
C
C               This expansion can be used to calculate the field and
C               reaction field of a unit dipole in -qq(ivoff+)- (or -J-, the
C               expansion centre) at polarizable point -P-.
C               Therefore, it is always calculated except when only
C               distributed monopoles are used to expand the source
C               AND reaction fields of the quantum motif
C
C        -----  Note: DRFOPQ gives d/dQ t(qq(ivoff+);P) !!
C
              IF (IFLDIN .GT. 1) THEN
                WT(IF+K,JF+L) = -O(K,L)
              ENDIF
              IF (IFLDOUT .GT. 1) THEN
                VR(IF+K,JF+L) = O(K,L)
              ENDIF
C       4-----
  280       CONTINUE
C
C      -----  If the source/reaction field is not expanded (i.e. only
C             distributed mono- /dipoles are used), only the potential
C             part is stored in -WT-/-VR-
C
            IF (IFLDIN .GT. 2) THEN
c             WT(IF+K,JF+4) = (W(K)-PQK)
c             WT(IF+K,JF+4) = (W(K)+PQK)
              WT(IF+K,JF+4) = W(K)
            ELSE
              WT(IF+K,JF+4) = -PQK
            ENDIF
C
            IF (IFLDOUT .GT. 2) THEN
c             VR(IF+K,JF+4) = - (W(K) - PQK)
c             VR(IF+K,JF+4) = - (W(K) + PQK)
              VR(IF+K,JF+4) = -W(K)
            ELSE
              VR(IF+K,JF+4) = PQK
            ENDIF
C
C      -----  Form -ZFN-: the fields in the polarizable points
C             due to the INTERNAL nuclei
C             -ZFN- is in fact the sum of the nuclear fields
C             and potentials
C
            IF (J .LE. NAT) THEN
              WT(IF+K,NEXP4+NGRAN+1) = WT(IF+K,NEXP4+NGRAN+1) - PQK*ZA
              VR(IF+K,NEXP4+NGRAN+1) = VR(IF+K,NEXP4+NGRAN+1) + PQK*ZA
            ENDIF
C     3-----
  290     CONTINUE
C
C    -----  End of polarizable points
C   2-----
  300   CONTINUE
C
C  -----  Fields (and fields dot normal) at boundary elements
C         due to charges, dipoles and charge distributions
C         (expanded) in the expansion centra
C
C  -----  Loop over boundary elements
C         
C         NOTE: DERIVATIVES NOT DONE YET!!!!
cahv - set to 0
C
        DO 400, NI = 1, NBEM
C   2-----
C    -----  QI: Vector from expansion centre -J- to
C           boundary element -NI- = (i-q)
C           DIST: Length of QI
C
          CALL DISTAB(q,XSURF(1,NI),QI,DIST)
          DMIND1 = ONE/DIST
          DMIN2 = DMIND1**2
          DMIN3 = DMIND1*DMIN2
          IF (KAPNOZ) THEN
            EXPKD = EXP(-(KAPPA*DIST))
            FKAPONE = ONE + (KAPPA*DIST)
          ENDIF
C
C    -----  FQINI:  Field of unit (positive) charge
C                   in -J- at -NI-, contracted with
C                   normal vector in -NI- = (i-q).n(i)/dist**3
C                   = f(q;i) . n(i)
C
C           It is also the negative of the potential of a dipole
C           in the direction of n(i) at the expansion centre
C
          FQINI = DMIN3*ddot(3,QI,1,XNORM(1,NI),1)
c         FQINI = DMIN3*ADOTB(QI,XNORM(1,NI),3)
C
C    -----  REACTION potential energy operator
C           at expansion centre, first term in the expansion
C           containing contribution of unit dipoles w(i)
C
C                XKIQ = K(i;q) S(i)
C
C        = (eps*(1+kappa*dist)*exp(-kappa*dist) - 1) f(i;q).n(i) S(i)
C
C    -----  The minus sign is a result of the use of the
C           inverted field: needed is f(i;q), FQINI=f(q;i).n(i)
C
          XKIQ = - (EPS*FKAPONE*EXPKD - ONE)*FQINI*AREA(NI)
C
          IF (KAPNOZ) THEN
C
C      -----  REACTION potential energy operator
C             at expansion centre, first term in the expansion,
C             containing contribution of unit charges z(i)
C
C                XLIQ = L(i;q) S(i)
C
C         = (1 - exp(-kappa*dist) ) V(i;q) S(i)
C
            XLIQ = (ONE - EXPKD)*DMIND1*AREA(NI)
          ENDIF
C
C    -----  Check if expansion centre coincides with
C           position of nucleus; exclude ambiguous atoms
C
          IF ((J .LE. NAT) .AND. (NAMBPT(J) .EQ. 0)) THEN
C     3-----
C      -----  POTENTIAL of all source nuclear charges in -J- at -NI-
C             scaled with 1/(2pi(1+eps)), as input for coupling
C             equations for w(i)
C
C                  V(q;i) = sum(q) qq(ivoff+q) / dist
C
C      -----  This is the boundary element part of -ZFN-
C
            IF (IFLDIN .GT. 2) THEN
cahv          WT(NPOL3+NI,NEXP4+NGRAN+1) =
cahv 1            WT(NPOL3+NI,NEXP4+NGRAN+1) + EPSFACT*ZA*DMIND1
              WT(NPOL3+NI,NEXP4+NGRAN+1) = 0.0
            ENDIF
C
C      -----  The reaction potential energy operator at the
C             expansion centre -J-, multiplied by the source
C             nuclear charge
C
C      -----  The interaction energy is evaluated through contracting
C             this array with the induced dipole density array
C
            IF (IFLDOUT .GT. 2) THEN
cahv          VR(NPOL3+NI,NEXP4+NGRAN+1) =
cahv 1                    VR(NPOL3+NI,NEXP4+NGRAN+1) + XKIQ*ZA
              VR(NPOL3+NI,NEXP4+NGRAN+1) =  0.0
            ENDIF
C
            IF (KAPNOZ) THEN
C       4-----
C        -----  MINUS FIELD of all source nuclear charges in -J- at -NI-
C               contracted with normal vector at -NI-,
C               scaled with eps/(2pi(1+eps)), as input for
C               coupling equations for z(i)
C
C                    f(q;i) = qq(ivoff+q) (i-q)/ dist**3
C
              IF (IFLDIN .GT. 2) THEN
cahv            WT(NPOL3+NI+NBEM,NEXP4+NGRAN+1) =
cahv 1          WT(NPOL3+NI+NBEM,NEXP4+NGRAN+1) - EPS*EPSFACT*ZA*FQINI
                WT(NPOL3+NI+NBEM,NEXP4+NGRAN+1) = 0.0
              ENDIF
C
C        -----  The reaction potential energy operator at the
C               expansion centre -J-, multiplied by the source
C               nuclear charge
C
C        -----  The interaction energy is evaluated through contracting
C               this array with the induced charge density array
C
              IF (IFLDOUT .GT. 2) THEN
cahv            VR(NPOL3+NI+NBEM,NEXP4+NGRAN+1) =
cahv 1          VR(NPOL3+NI+NBEM,NEXP4+NGRAN+1) + XLIQ*ZA
                VR(NPOL3+NI+NBEM,NEXP4+NGRAN+1) = 0.0
              ENDIF
C       4-----
            ENDIF
C     3-----
          ENDIF
C
C    -----  Expansion of source potential and field
C           of surface charge distribution
C           and reaction potentials at and around the expansion centra
C
C             GENERAL:
C
C         Expansion in Taylor series in x around q:
C
C   For source potential: V(x;i) = V(q;i) + del(x) V(x;i) (x=q) .(x-q)
C         = V(q;i) - f(q;i).q + f(q;i).x
C
C   For source field:     f(x;i) = f(q;i) + del(x) f(x;i) (x=q) .(x-q)
C         = f(q;i) + t(q;i).q - t(q;i).x
C
C   For reaction potential due to induced dipoles (operator)
C                         K(i;x) = K(i;q) + del(x) K(i;x) (x=q) .(x-q)
C   = K(i;q) - [(eps*(1+kd)*exp(-kd) -1)*t(i;q).n(i) -
C                eps*(kappa**2)*exp(-kd)*f(i;q).n(i)*(q-i)].q
C            + [(eps*(1+kd)*exp(-kd) -1)*t(i;q).n(i) -
C                eps*(kappa**2)*exp(-kd)*f(i;q).n(i)*(q-i)].x
C
C   For reaction potential due to induced charges (operator)
C                         L(i;x) = L(i;q) + del(x) L(i;x) (x=q) .(x-q)
C   = L(i;q) - [(1+kd)*exp(-kd) -1)*f(i;q)].q
C            + [(1+kd)*exp(-kd) -1)*f(i;q)].x
C
C-----
C      SOURCE POTENTIAL IN -NI- (scaled with 1/(2pi(1+eps))
C
          IF (IFLDIN .GT. 2) THEN
C-----
C      V(q;i) - f(q;i).q
C-----
            WT(NPOL3+NI,JF+4) = EPSFACT*(DMIND1-
     +          DMIN3*ddot(3,QI,1,q,1))
c    +          DMIN3*ADOTB(QI,q,3))
          ELSE
            WT(NPOL3+NI,JF+4) = EPSFACT*DMIND1
          ENDIF
cahv
          WT(NPOL3+NI,JF+4) = 0.0           
cahv
C-----
C      the operator f(q;i) (.x)
C-----
          IF (IFLDIN .GT. 1) THEN
            DO 320, K = 1, 3
              WT(NPOL3+NI,JF+K) = EPSFACT*DMIN3*QI(K)
cahv
          WT(NPOL3+NI,JF+K) = 0.0                  
cahv
  320       CONTINUE
          ENDIF
C
C    -----  Calculate del(i) f(q;i) = t(q;i) = t(i;q) = T
C
          CALL DRFTPQ(QI,DMIND1,DMIN3,0,ONE,ONE,T)
C
C    -----  Calculate TNI = t(q;i) . n(i)
C
          CALL MATVEC(T,XNORM(1,NI),TNI,3,.FALSE.)
C
          IF (KAPNOZ) THEN
C     3-----
C      MINUS SOURCE FIELD IN -NI- contracted with normal vector in -NI-
C                           and scaled with (eps/(2pi(1+eps))
C-----
C
C      -----  Calculate contraction of t(q;i) with origin of
C             expansion q: TQ = t(q;i).q
C
            CALL MATVEC(T,q,TQ,3,.FALSE.)
C-----
C      - [f(q;i).n(i) + [t(q;i).q].n(i)]
C-----
            IF (IFLDIN .GT. 2) THEN
              WT(NPOL3+NI+NBEM,JF+4) = - EPS*EPSFACT*
     1        (FQINI + ddot(3,TQ,1,XNORM(1,NI),1))
c    1        (FQINI + ADOTB(TQ,XNORM(1,NI),3))
            ELSE
              WT(NPOL3+NI+NBEM,JF+4) = - EPS*EPSFACT*FQINI
            ENDIF
cahv
            WT(NPOL3+NI+NBEM,JF+4) = 0.0                 
cahv
C-----
C      The operator - [-t(q;i) (.x)] .n(i)
C                 =   [ t(q;i).n(i)] (.x)
C-----
            IF (IFLDIN .GT. 1) THEN
              DO 330, K = 1, 3
                WT(NPOL3+NI+NBEM,JF+K) = EPS*EPSFACT*TNI(K)
cahv
            WT(NPOL3+NI+NBEM,JF+K) = 0.0                
cahv
  330         CONTINUE
            ENDIF
C     3-----
          ENDIF
C-----
C      REACTION POTENTIAL DUE TO DIPOLE IN -NI- AT -X-
C-----
C
C    -----  Calculate del(x) K(i;x) (x=q) *S(i) = DELKIQ
C
C      The minus sign for the second term remains, since both
C      field dot normal(FQINI) and vector (QI) are from -J- to -NI-
C      whereas the reverse of vectors is needed
C
          IF (IFLDOUT .GT. 1) THEN
            DO 340, K = 1, 3
              DELKIQ(K)=  ((EPS*FKAPONE*EXPKD - ONE)*TNI(K) -
     1                  EPS*KAPPAS*EXPKD*FQINI*QI(K))*AREA(NI)
C-----
C      The operator del(x) K(i;x) (x=q) *S(i) (.x)
C-----
              VR(NPOL3+NI,JF+K)= DELKIQ(K)
cahv
              VR(NPOL3+NI,JF+K)= 0.0       
cahv
  340       CONTINUE
          ENDIF
C
          IF (IFLDOUT .GT. 2) THEN
C-----
C      [K(i;q) - del(x) K(i;x) (x=q) .q] S(i)
C-----
            VR(NPOL3+NI,JF+4)=  XKIQ - ddot(3,DELKIQ,1,q,1)
c           VR(NPOL3+NI,JF+4)=  XKIQ - ADOTB(DELKIQ,q,3)
          ELSE
C-----
C       K(i;q)
C-----
            VR(NPOL3+NI,JF+4)=  XKIQ
          ENDIF
cahv
            VR(NPOL3+NI,JF+4)=  0.0  
cahv
C
          IF (KAPNOZ) THEN
C     3-----
C      REACTION POTENTIAL DUE TO CHARGE IN -NI- AT -X-
C-----
C      Calculate del(x) L(i;x) (x=q) *S(i) = DELLIQ
C
C      Minus sign is a result of using f(q;i) = (i-q)/dist**3,
C      whereas f(i;q) is needed
C
            IF (IFLDOUT .GT. 1) THEN
              DO 350, K = 1, 3
                DELLIQ(K)=-(ONE-FKAPONE*EXPKD)*QI(K)*DMIN3*AREA(NI)
C-----
C      The operator [del(x) L(i;x) (x=q)] *S(i) (.x)
C-----
                VR(NPOL3+NI+NBEM,JF+K) = DELLIQ(K)
cahv
                VR(NPOL3+NI+NBEM,JF+K) = 0.00     
cahv
  350         CONTINUE
            ENDIF
C
            IF (IFLDOUT .GT. 2) THEN
C-----
C      [L(i;q) - del(x) L(i;x) (x=q) . q] *S(i)
C-----
              VR(NPOL3+NI+NBEM,JF+4)=  XLIQ - ddot(3,DELLIQ,1,q,1)
c             VR(NPOL3+NI+NBEM,JF+4)=  XLIQ - ADOTB(DELLIQ,q,3)
            ELSE
C-----
C       L(i;q)
C-----
              VR(NPOL3+NI+NBEM,JF+4) = XLIQ
            ENDIF
cahv
              VR(NPOL3+NI+NBEM,JF+4) = 0.0  
cahv
C     3-----
          ENDIF
C   2-----
  400   CONTINUE
 500    CONTINUE
C 1-----
      END DO
C
      RETURN
      END
      SUBROUTINE DRFGRAX(KXYZ,GRAD,D,OL,DX,DY,DZ,
     1           IEXPC,IJBIT)
C------
C      Electronic contribution to RF gradient 
C      on expansion centres
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION GRAD(3,NEXP)
      DIMENSION D(NCHD), OL(NCHD), DX(NCHD), DY(NCHD), DZ(NCHD)
      DIMENSION IEXPC(NCHD), IJBIT(NCHD)
C
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
C
      DIMENSION P(4)
C
      DATA ONE,TWO /1.D00, 2.D00/
C
C-----  BEGIN
C
      IJ = 0
C
C-----  Loop over charge distributions
C
      DO I = 1, NUM
        DO J = 1, I
C   1-----
          IJ = IJ + 1
C
          FAC = TWO
          IF (I.EQ.J) FAC = ONE
C
C    -----  Get expansion centre
C           of IJ-th charge distribution
C
          IJEXP = IEXPC(IJ)
          JF = (IJEXP-1)*4
          NZF = JF*3
C
C    -----  Collect overlap, dipole
C
          P(1) = FAC*D(IJ)*DX(IJ)
          P(2) = FAC*D(IJ)*DY(IJ)
          P(3) = FAC*D(IJ)*DZ(IJ)
          P(4) = FAC*D(IJ)*OL(IJ)
C
C    -----  Interaction with external charges
C
          DO IGR = 1, NGRAN
            DO L = 1, 4
C
C        -----  NOTE: --sign because coupling is to 
C               electrons
C
              GRAD(KXYZ,IJEXP) = GRAD(KXYZ,IJEXP)
     1             - DZFA(NZF+(KXYZ-1)*4+L,IGR)*P(L)
            ENDDO
          ENDDO
C   1-----
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(grad,3,nexp,2,'grax')
C-----
      RETURN
      END
      SUBROUTINE DRFGRAX2(KXYZ,GRAD,XEXP,DA,DB,V,OL,DX,DY,DZ,
     1           RXX,RYY,RZZ,RXY,RXZ,RYZ,OMEGA,
     1           IEXPC,IJBIT)
C------
C      Electronic contribution to RF gradient 
C      on expansion centres
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION GRAD(3,NEXP), XEXP(3,NEXP)
      DIMENSION OL(NCHD), DX(NCHD), DY(NCHD), DZ(NCHD)
      DIMENSION RXX(NCHD),RYY(NCHD),RZZ(NCHD),
     1          RXY(NCHD),RXZ(NCHD),RYZ(NCHD)
      DIMENSION OMEGA(NWTC,NWTC)
      DIMENSION IEXPC(NCHD), IJBIT(NCHD)
      DIMENSION V(NUM,NUM)
      DIMENSION DA(NCHD), DB(NCHD)
C
      LOGICAL UHF,ROHF,RGVB,ROGVB
      LOGICAL CORE,OPEN,PAIR
C
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
C
      integer ihlp
      common /ihelp/ ihlp(maxorb)
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
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
C
c
      real*8 alphij, betaij
      integer nopen, ncorb, nopset, npairs
      common /alpbet/ alphij(325),betaij(325),nopen(10),ncorb,
     +                nopset,npairs
c
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopenn
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopenn,nope,noe(10)
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
C
C
      DIMENSION OP(4), OQ(4)
      DIMENSION RR(4,4)
C
      DATA THRESH2 /1.0D-06/
C
      DATA ZERO,ONE,TWO,PT5,PT25/0.0D00,1.0D00,2.0D00,0.5D00,0.25D00/
C
C-----  BEGIN
C
C
      UHF = SCFTYP .EQ. 'UHF'
      ROHF = (SCFTYP.EQ. 'RHF') .AND. (NOPSET.GT.0)
      RGVB = SCFTYP .EQ. 'GVB'
      ROGVB = ROHF .OR. RGVB
C
      CORE = NCORB.NE.0
      OPEN = NOPSET.NE.0
      PAIR = NPAIRS.NE.0
C
C-----  Set up mapping of FOCK matrices
C
      IF (ROGVB) THEN
C 1-----
        NONE = 1
        IF (.NOT. CORE) NONE = 0
        IF (CORE) THEN
C   2-----
          DO 10, I = 1, NCORB
            IHLP(I) = 1
   10     CONTINUE
C   2-----
        ENDIF
        NOP = 0
C
        IF (OPEN) THEN
C   2-----
          DO 30, ISET = 1, NOPSET
C     3-----
            IOP = NOPEN(ISET)
            DO 20, I = 1, IOP
              IHLP(NCORB+NOP+I) = NONE + ISET
   20       CONTINUE
            NOP = NOP + IOP
C     3-----
   30     CONTINUE
C   2-----
        ENDIF
C
C
        IF (PAIR) THEN
C   2-----
          NGEM = 2*NPAIRS
          DO 40, IGEM = 1, NGEM
            IHLP(NCORB+NOP+IGEM) = NONE + NOPSET + IGEM
   40     CONTINUE
C   2-----
        ENDIF
C
        NORB = NCORB + NOP + 2*NPAIRS
        NCO1 = NCORB + 1
C
C  -----  Get vectors and constuct density
C
        CALL DAREAD(IDAFh,IODA,V,NUM*NUM,15)
        IF (CORE) THEN
C   2-----
          DO 70, I = 1, NUM
C     3-----
            DO 60, J = 1, I
C       4-----
              DUM = ZERO
              DO 50, K = 1, NCORB
                DUM = DUM + V(I,K)*V(J,K)
   50         CONTINUE
              IJ = IA(I) + J
              DB(IJ) = DUM
C       4-----
   60       CONTINUE
C     3-----
   70     CONTINUE
C   2-----
        ENDIF
C 1-----
      ENDIF
C
      NZFP = NEXP*4
      NZFN = NEXP*4 + NGRAN + 1
C
      IJ = 0
C
C-----  Loop over charge distributions
C
      DO I = 1, NUM
        DO J = 1, I
C   1-----
          IJ = IJ + 1
C
          IF (IJBIT(IJ) .EQ. 0) GOTO 100
C
          FAC = TWO
          IF (I.EQ.J) FAC = ONE
          FACTIJ = FAC
C
C    -----  Get expansion centre
C           of IJ-th charge distribution
C
          IJEXP = IEXPC(IJ)
          JF = (IJEXP-1)*4
C
C    -----  Collect overlap, dipole
C
          P(1) = DX(IJ)
          P(2) = DY(IJ)
          P(3) = DZ(IJ)
          P(4) = OL(IJ)
C
C    -----  Interaction with externally induced dipoles
C
          DO IGR = 1, NGRAN
	    contr = 0.0d0
            DO L = 1, 4
C       2-----
C        -----  NOTE: --sign because coupling is to 
C               electrons
C
C        -----  Electrons are source
C
              contr = contr
     1             - pt5*FACTIJ*DA(IJ)
     2             *omega(nzfp+igr,jf+l)*p(l)
	      contr = contr
     1             - pt5*FACTIJ*DA(IJ)
     2             *omega(jf+l,nzfp+igr)*p(l)
C       2-----
            ENDDO
            grad(kxyz,ijexp) = grad(kxyz,ijexp) + contr
          ENDDO
C
C    -----  Screening of the nuclear attraction
C
C      -----  Electrons are source
C
	    do ll = 1, nat
	      q(1) = xexp(1,ll)
	      q(2) = xexp(2,ll)
	      q(3) = xexp(3,ll)
	      q(4) = one
              CALL DRFOAB(LL,IJEXP,NWTC,OMEGA)
              CALL MATVEC(OMGAB,P,OP,4,.FALSE.)
              contr = - pt5*czan(ll)*DA(IJ)*FACTIJ*
     +        ddot(4,OP,1,q,1)
c    +        ADOTB(OP,q,4)
            GRAD(KXYZ,ijexp) = GRAD(KXYZ,ijexp) + contr
	    enddo
C
C      -----  Electrons are recipients
C
	    do ll = 1, nat
	      q(1) = xexp(1,ll)
	      q(2) = xexp(2,ll)
	      q(3) = xexp(3,ll)
	      q(4) = one
              CALL DRFOAB(IJEXP,LL,NWTC,OMEGA)
              CALL MATVEC(OMGAB,q,OQ,4,.FALSE.)
            contr = - pt5*czan(ll)*DA(IJ)*FACTIJ*
     +      ddot(4,OQ,1,P,1)
c    +      ADOTB(OQ,P,4)
            GRAD(KXYZ,ijexp) = GRAD(KXYZ,ijexp) + contr
	    enddo
C
C    -----  Self energy; dispersion
C
          IF (GAMDRF .NE. ZERO) THEN
C     2-----
            call drfoab(ijexp,ijexp,nwtc,omega)
            rr(1,1) = rxx(ij)
            rr(1,2) = rxy(ij)
            rr(2,1) = rxy(ij)
            rr(1,3) = rxz(ij)
            rr(3,1) = rxz(ij)
            rr(1,4) = dx(ij)
            rr(4,1) = dx(ij)
            rr(2,2) = ryy(ij)
            rr(2,3) = ryz(ij)
            rr(3,2) = ryz(ij)
            rr(2,4) = dy(ij)
            rr(4,2) = dy(ij)
            rr(3,3) = rzz(ij)
            rr(3,4) = dz(ij)
            rr(4,3) = dz(ij)
            rr(4,4) = ol(ij)
            contr = ddot(16,rr,1,omgab,1)*da(ij)*factij*pt5
c           contr = adotb(rr,omgab,16)*da(ij)*factij*pt5
            grad(kxyz,ijexp) = grad(kxyz,ijexp) + contr*gamdrf
C     2-----
          ENDIF
C
C    -----  Two-electron interaction
C
	  kl = 0
          DO K = 1, num
            DO L = 1, K
C       2-----
	      kl = kl + 1
              IF (IJBIT(KL) .EQ. 0) GOTO 150
              KLEXP = IEXPC(KL)
C
C        -----  Collect overlap, dipole
C
              q(1) = DX(KL)
              q(2) = DY(KL)
              q(3) = DZ(KL)
              q(4) = OL(KL)
C
c             CALL DRFOAB(IJEXP,KLEXP,NWTC,OMEGA)
C
c             CALL MATVEC(OMGAB,q,OQ,4,.FALSE.)
C
              CALL DRFOAB(KLEXP,IJEXP,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,P,OP,4,.FALSE.)
C
C        -----  Calculate addresses
C
              ik = ia(max(i,k)) + min(i,k)
              jk = ia(max(j,k)) + min(j,k)
              il = ia(max(i,l)) + min(i,l)
              jl = ia(max(j,l)) + min(j,l)
C
              II = I
              JJ = J
              KK = K
              LL = L
C
C        -----  Calculate density factors
C
              CALL DRFDAB(DA,DB,V,DCOUL,DEXCH,NORB,ROHF,RGVB,UHF)
C
              IF ((ABS(DCOUL) .GT. THRESH2) .OR.
     1            (ABS(DEXCH) .GT. THRESH2)) THEN
C         5-----
                FACTKL = TWO
                IF (L .EQ. K) FACTKL = ONE
                FACT = FACTIJ*FACTKL
		fact=fact*pt5
C
C          -----  Factor 0.25 for double counting twice
C
                jexp = klexp
		contr = 
     1             PT25*(DCOUL+DEXCH)*FACT*
     2             ddot(4,OP,1,q,1)
c    2             ADOTB(OP,q,4)
		grad(kxyz,jexp) = grad(kxyz,jexp) + contr
C
C         5-----
              ENDIF
C
 150          CONTINUE
C       2-----
            ENDDO
          ENDDO
C
 100      CONTINUE
C   1-----
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(grad,3,nexp,2,'grax2')
C-----
      RETURN
      END
      SUBROUTINE DRFGRAX2B(KXYZ,GRAD,XEXP,DA,DB,V,OL,DX,DY,DZ,
     1           RXX,RYY,RZZ,RXY,RXZ,RYZ,OMEGA,
     1           IEXPC,IJBIT)
C------
C      Electronic contribution to RF gradient 
C      on expansion centres
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION GRAD(3,NEXP), XEXP(3,NEXP)
      DIMENSION OL(NCHD), DX(NCHD), DY(NCHD), DZ(NCHD)
      DIMENSION RXX(NCHD),RYY(NCHD),RZZ(NCHD),
     1          RXY(NCHD),RXZ(NCHD),RYZ(NCHD)
      DIMENSION OMEGA(NWTC,NWTC)
      DIMENSION IEXPC(NCHD), IJBIT(NCHD)
      DIMENSION V(NUM,NUM)
      DIMENSION DA(NCHD), DB(NCHD)
C
      LOGICAL UHF,ROHF,RGVB,ROGVB
      LOGICAL CORE,OPEN,PAIR
C
C
C
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
C
      integer ihlp
      common /ihelp/ ihlp(maxorb)
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
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
C
c
      real*8 alphij, betaij
      integer nopen, ncorb, nopset, npairs
      common /alpbet/ alphij(325),betaij(325),nopen(10),ncorb,
     +                nopset,npairs
c
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopenn
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopenn,nope,noe(10)
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
C
      DIMENSION OP(4), OQ(4)
      DIMENSION RR(4,4)
C
      DATA THRESH2 /1.0D-06/
C
      DATA ZERO,ONE,TWO,PT5,PT25/0.0D00,1.0D00,2.0D00,0.5D00,0.25D00/
C
C-----  BEGIN
C
C
      UHF = SCFTYP .EQ. 'UHF'
      ROHF = (SCFTYP.EQ. 'RHF') .AND. (NOPSET.GT.0)
      RGVB = SCFTYP .EQ. 'GVB'
      ROGVB = ROHF .OR. RGVB
C
      CORE = NCORB.NE.0
      OPEN = NOPSET.NE.0
      PAIR = NPAIRS.NE.0
C
C-----  Set up mapping of FOCK matrices
C
      IF (ROGVB) THEN
C 1-----
        NONE = 1
        IF (.NOT. CORE) NONE = 0
        IF (CORE) THEN
C   2-----
          DO 10, I = 1, NCORB
            IHLP(I) = 1
   10     CONTINUE
C   2-----
        ENDIF
        NOP = 0
C
        IF (OPEN) THEN
C   2-----
          DO 30, ISET = 1, NOPSET
C     3-----
            IOP = NOPEN(ISET)
            DO 20, I = 1, IOP
              IHLP(NCORB+NOP+I) = NONE + ISET
   20       CONTINUE
            NOP = NOP + IOP
C     3-----
   30     CONTINUE
C   2-----
        ENDIF
C
C
        IF (PAIR) THEN
C   2-----
          NGEM = 2*NPAIRS
          DO 40, IGEM = 1, NGEM
            IHLP(NCORB+NOP+IGEM) = NONE + NOPSET + IGEM
   40     CONTINUE
C   2-----
        ENDIF
C
        NORB = NCORB + NOP + 2*NPAIRS
        NCO1 = NCORB + 1
C
C  -----  Get vectors and constuct density
C
        CALL DAREAD(IDAFh,IODA,V,NUM*NUM,15)
        IF (CORE) THEN
C   2-----
          DO 70, I = 1, NUM
C     3-----
            DO 60, J = 1, I
C       4-----
              DUM = ZERO
              DO 50, K = 1, NCORB
                DUM = DUM + V(I,K)*V(J,K)
   50         CONTINUE
              IJ = IA(I) + J
              DB(IJ) = DUM
C       4-----
   60       CONTINUE
C     3-----
   70     CONTINUE
C   2-----
        ENDIF
C 1-----
      ENDIF
C
      NZFP = NEXP*4
      NZFN = NEXP*4 + NGRAN + 1
C
      IJ = 0
C
C-----  Loop over charge distributions
C
      DO I = 1, NUM
        DO J = 1, I
C   1-----
          IJ = IJ + 1
C
          IF (IJBIT(IJ) .EQ. 0) GOTO 100
C
          FAC = TWO
          IF (I.EQ.J) FAC = ONE
          FACTIJ = FAC
C
C    -----  Get expansion centre
C           of IJ-th charge distribution
C
          IJEXP = IEXPC(IJ)
          JF = (IJEXP-1)*4
C
C    -----  Collect overlap, dipole
C
          P(1) = DX(IJ)
          P(2) = DY(IJ)
          P(3) = DZ(IJ)
          P(4) = OL(IJ)
C
C    -----  Interaction with externally induced dipoles
C
          DO IGR = 1, NGRAN
            DO L = 1, 4
C       2-----
C        -----  NOTE: --sign because coupling is to 
C               electrons
C
C        -----  Electrons are recipients
C
              GRAD(KXYZ,IJEXP) = GRAD(KXYZ,IJEXP)
     1             - pt5*FACTIJ*DA(IJ)
     2             * omega(jf+l,nzfp+igr)*p(l)
              GRAD(KXYZ,IJEXP) = GRAD(KXYZ,IJEXP)
     1             - pt5*FACTIJ*DA(IJ)
     2             * omega(nzfp+igr,jf+l)*p(l)
C       2-----
            ENDDO
          ENDDO
C
C    -----  Screening of the nuclear attraction
C
C      -----  Electrons are source
C
	    do ll = 1, nat
	      q(1) = xexp(1,ll)
	      q(2) = xexp(2,ll)
	      q(3) = xexp(3,ll)
	      q(4) = one
              CALL DRFOAB(LL,IJEXP,NWTC,OMEGA)
              CALL MATVEC(OMGAB,P,OP,4,.FALSE.)
            GRAD(KXYZ,ijexp) = GRAD(KXYZ,ijexp)
     1           - pt5*czan(ll)*DA(IJ)*FACTIJ*
     2             ddot(4,OP,1,q,1)
c    2             ADOTB(OP,q,4)
	    enddo
C
C      -----  Electrons are recipients
C
	    do ll = 1, nat
	      q(1) = xexp(1,ll)
	      q(2) = xexp(2,ll)
	      q(3) = xexp(3,ll)
	      q(4) = one
              CALL DRFOAB(IJEXP,LL,NWTC,OMEGA)
              CALL MATVEC(OMGAB,q,OQ,4,.FALSE.)
            GRAD(KXYZ,ijexp) = GRAD(KXYZ,ijexp)
     1           - pt5*czan(ll)*DA(IJ)*FACTIJ*
     2             ddot(4,OQ,1,P,1)
c    2             ADOTB(OQ,P,4)
	    enddo
C
C
C    -----  Self energy; dispersion
C
          IF (GAMDRF .NE. ZERO) THEN
C     2-----
            call drfoab(ijexp,ijexp,nwtc,omega)
            rr(1,1) = rxx(ij)
            rr(1,2) = rxy(ij)
            rr(2,1) = rxy(ij)
            rr(1,3) = rxz(ij)
            rr(3,1) = rxz(ij)
            rr(1,4) = dx(ij)
            rr(4,1) = dx(ij)
            rr(2,2) = ryy(ij)
            rr(2,3) = ryz(ij)
            rr(3,2) = ryz(ij)
            rr(2,4) = dy(ij)
            rr(4,2) = dy(ij)
            rr(3,3) = rzz(ij)
            rr(3,4) = dz(ij)
            rr(4,3) = dz(ij)
            rr(4,4) = ol(ij)
            contr = ddot(16,rr,1,omgab,1)*da(ij)*factij*pt5
c           contr = adotb(rr,omgab,16)*da(ij)*factij*pt5
            grad(kxyz,ijexp) = grad(kxyz,ijexp) + contr*gamdrf
C     2-----
          ENDIF
C
C    -----  Two-electron interaction
C
	  kl = 0
          DO K = 1, num
            DO L = 1, K
C       2-----
	      kl = kl + 1
              IF (IJBIT(KL) .EQ. 0) GOTO 150
              KLEXP = IEXPC(KL)
C
C        -----  Collect overlap, dipole
C
              q(1) = DX(KL)
              q(2) = DY(KL)
              q(3) = DZ(KL)
              q(4) = OL(KL)
C
              CALL DRFOAB(IJEXP,KLEXP,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,q,OQ,4,.FALSE.)
C
C        -----  Calculate addresses
C
              ik = ia(max(i,k)) + min(i,k)
              jk = ia(max(j,k)) + min(j,k)
              il = ia(max(i,l)) + min(i,l)
              jl = ia(max(j,l)) + min(j,l)
              II = I
              JJ = J
              KK = K
              LL = L
C
C        -----  Calculate density factors
C
              CALL DRFDAB(DA,DB,V,DCOUL,DEXCH,NORB,ROHF,RGVB,UHF)
C
              IF ((ABS(DCOUL) .GT. THRESH2) .OR.
     1            (ABS(DEXCH) .GT. THRESH2)) THEN
C         5-----
                FACTKL = TWO
                IF (L .EQ. K) FACTKL = ONE
                FACT = FACTIJ*FACTKL
		fact = fact*pt5
C
C          -----  Factor 0.25 for double counting twice
C
c               VAL = VAL*FACT
c               SCOUL = SCOUL + DCOUL*VAL
c               SEXCH = SEXCH + DEXCH*VAL
C
                jexp = klexp
C
                GRAD(KXYZ,JEXP) = GRAD(KXYZ,JEXP) +
     1             PT25*(DCOUL+DEXCH)*FACT*ddot(4,OQ,1,P,1)
c    1             PT25*(DCOUL+DEXCH)*FACT*ADOTB(OQ,P,4)
C         5-----
              ENDIF
C
 150          CONTINUE
C       2-----
            ENDDO
          ENDDO
C
 100      CONTINUE
C   1-----
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(grad,3,nexp,2,'grax2b')
C-----
      RETURN
      END
      subroutine drfdomga(ieps,relay,wt,vr,omega,indx,dipix)
c------
c      this routine constructs the matrix omega, containing the
c      reaction field contributions to the total energy for
c      the expanded electrons (formally for one electron per
c      expansion centre), the nuclei and the external charges.
c
c      it is formed by contracting the reaction potential -vr-
c      with the (formal) solutions of the relay equations, i.e.
c      the solutions of ri=s, with r the relay matrix, and s the
c      -wt- matrix.
c
c      depending on the type of source and recipients (i.e.
c      mono-/dipole expansion, molecular, separate nuclei and
c      electrons), parts or all of omega will be used for analysis.
c      also, the matrix elements of -wt- and -vr- differ, depending
c      on the type of expansion, resulting in differing omega
c      matrix elements. however, this does not influence the
c      calculation of omega. omega might contain mainly zero's.
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
      dimension relay(ndim,ndim), wt(nwtr,nwtc), vr(nwtr,nwtc),
     1 omega(nwtc,nwtc)
c
      integer indx(ndim)
      dimension dipix(nwtc,ndim)
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
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
c-----  end of common blocks
c
c-----  data
c
      data zero,one,pt5,pt75 /0.0d00,1.0d00,0.5d00,0.75d00/
c
c-----  begin
c
      call clear(omega,nomga)
c
c-----  solve linear eqs. and form -omega-
c
      do 100, i = 1, nwtc
        call luelmf(relay,wt(1,i),indx,ndim,ndim,auxx)
        if (idrfout.eq.3 .or. imcout .eq. 5) then
          call hatout(auxx,1,ndim ,22,'auxx_uit')
        endif
        do k = 1, ndim
          dipix(i,k) = auxx(k)
        enddo
c
c  -----  evaluation of energy contributions after solving
c         ri=s
c
c  -----  first for electrons (expansion, expansion centre)
c         then for external charges
c
        do 200, j = 1, nwtc
          omega(i,j) = ddot(ndim,auxx,1,vr(1,j),1)
c         omega(i,j) = adotb(auxx,vr(1,j),ndim)
  200   continue
  100 continue
c
      if (idrfout.eq.3 .or. imcout .eq. 5)
     1        call hatout(omega,nwtc,nwtc,2,'domega')
      return
      end
      SUBROUTINE MOMGRAD(IATM,KXYZ,GRAD,D,DB,OL,DX,DY,DZ,
     1           DS,DXD,DYD,DZD,OMEGA,
     1           IEXPC,IJBIT)
C------
C      Electronic contribution to RF gradient 
C      on expansion centres
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION GRAD(3,NEXP)
      DIMENSION DB(NCHD)
      DIMENSION D(NCHD), OL(NCHD), DX(NCHD), DY(NCHD), DZ(NCHD)
      DIMENSION DS(NCHD), DXD(NCHD), DYD(NCHD), DZD(NCHD)
      DIMENSION OMEGA(NWTC,NWTC)
      DIMENSION IEXPC(NCHD), IJBIT(NCHD)
C
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
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
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
C
      DIMENSION P(4)
C
      DATA ONE,TWO /1.D00, 2.D00/
C
C-----  BEGIN
C
      IJ = 0
C
C-----  Loop over charge distributions
C
      DO I = 1, NUM
        DO J = 1, I
C   1-----
          IJ = IJ + 1
C
          FAC = TWO
c         IF (I.EQ.J) FAC = ONE
C
C    -----  Get expansion centre
C           of IJ-th charge distribution
C
          IJEXP = IEXPC(IJ)
          JF = (IJEXP-1)*4
          NZF = JF*3
C
C    -----  Collect overlap, dipole
C
          P(1) = FAC*D(IJ)*DXD(IJ)
          P(2) = FAC*D(IJ)*DYD(IJ)
          P(3) = FAC*D(IJ)*DZD(IJ)
          P(4) = FAC*D(IJ)*DS(IJ)
C
C    -----  Interaction with external charges
C
          DO IGR = 1, NGRAN
C
C        -----  NOTE: --sign because coupling is to 
C               electrons
C
            GRAD(KXYZ,IATM) = GRAD(KXYZ,IATM)
     1             - ddot(4,p,1,zfa(jf+1,igr),1)
c    1             - adotb(p,zfa(jf+1,igr),4)
          ENDDO
C   1-----
        ENDDO
      ENDDO
      if (idrfout .ge. 3) 
     1  call hatout(grad,3,nexp,2,'grad')
C-----
      RETURN
      END
      SUBROUTINE DRFMOM2(iatm,KXYZ,GRAD,DA,DB,V,OL,DX,DY,DZ,
     1           RXX,RYY,RZZ,RXY,RXZ,RYZ,
     1           RYX,RZX,RZY,
     1           ds,dxd,dyd,dzd,OMEGA,
     1           IEXPC,IJBIT)
C------
C      Electronic contribution to RF gradient 
C      on expansion centres
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      DIMENSION GRAD(3,NEXP)
      DIMENSION OL(NCHD), DX(NCHD), DY(NCHD), DZ(NCHD)
      DIMENSION ds(NCHD), DXd(NCHD), DYd(NCHD), DZd(NCHD)
      DIMENSION RXX(NCHD),RYY(NCHD),RZZ(NCHD),
     1          RXY(NCHD),RXZ(NCHD),RYZ(NCHD),
     2          RYX(NCHD),RZX(NCHD),RZY(NCHD)
      DIMENSION OMEGA(NWTC,NWTC)
      DIMENSION IEXPC(NCHD), IJBIT(NCHD)
      DIMENSION V(NUM,NUM)
      DIMENSION DA(NCHD), DB(NCHD)
C
      LOGICAL UHF,ROHF,RGVB,ROGVB
      LOGICAL CORE,OPEN,PAIR
C
C
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
C
      integer ihlp
      common /ihelp/ ihlp(maxorb)
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
      real*8 p, q, omgab
      common /drfexp/ p(4),q(4),omgab(4,4)
c
c
      real*8 val
      integer i, j, k, l
      common /drfint/ i,j,k,l,val
c
c
      integer ii, jj, kk, ll, ij, kl, ik, jl, il, jk
      common /drfindx/ ii,jj,kk,ll,ij,kl,ik,jl,il,jk
c
C
c
      real*8 alphij, betaij
      integer nopen, ncorb, nopset, npairs
      common /alpbet/ alphij(325),betaij(325),nopen(10),ncorb,
     +                nopset,npairs
c
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopenn
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopenn,nope,noe(10)
c
c
      character*8 scftyp
      common /scfopt2/ scftyp
c
C
      DIMENSION OP(4), OQ(4)
      DIMENSION RR(4,4)
      dimension OsP(4), OdQ(4), sP(4), dQ(4), dN(4), OdN(4)
C
      DATA THRESH2 /1.0D-20/
C
      DATA ZERO,ONE,TWO,PT5,PT25/0.0D00,1.0D00,2.0D00,0.5D00,0.25D00/
C
C-----  BEGIN
C
      if (idrfout .ge. 3) then
	 call hatout(ds,num,num,3,'ds')
	 call hatout(dxd,num,num,3,'dxd')
	 call hatout(dyd,num,num,3,'dyd')
	 call hatout(dzd,num,num,3,'dzd')
      endif
C
      UHF = SCFTYP .EQ. 'UHF'
      ROHF = (SCFTYP.EQ. 'RHF') .AND. (NOPSET.GT.0)
      RGVB = SCFTYP .EQ. 'GVB'
      ROGVB = ROHF .OR. RGVB
C
      CORE = NCORB.NE.0
      OPEN = NOPSET.NE.0
      PAIR = NPAIRS.NE.0
C
C-----  Set up mapping of FOCK matrices
C
      IF (ROGVB) THEN
C 1-----
        NONE = 1
        IF (.NOT. CORE) NONE = 0
        IF (CORE) THEN
C   2-----
          DO 10, I = 1, NCORB
            IHLP(I) = 1
   10     CONTINUE
C   2-----
        ENDIF
        NOP = 0
C
        IF (OPEN) THEN
C   2-----
          DO 30, ISET = 1, NOPSET
C     3-----
            IOP = NOPEN(ISET)
            DO 20, I = 1, IOP
              IHLP(NCORB+NOP+I) = NONE + ISET
   20       CONTINUE
            NOP = NOP + IOP
C     3-----
   30     CONTINUE
C   2-----
        ENDIF
C
C
        IF (PAIR) THEN
C   2-----
          NGEM = 2*NPAIRS
          DO 40, IGEM = 1, NGEM
            IHLP(NCORB+NOP+IGEM) = NONE + NOPSET + IGEM
   40     CONTINUE
C   2-----
        ENDIF
C
        NORB = NCORB + NOP + 2*NPAIRS
        NCO1 = NCORB + 1
C
C  -----  Get vectors and constuct density
C
        CALL DAREAD(IDAF,IODA,V,NUM*NUM,15)
        IF (CORE) THEN
C   2-----
          DO 70, I = 1, NUM
C     3-----
            DO 60, J = 1, I
C       4-----
              DUM = ZERO
              DO 50, K = 1, NCORB
                DUM = DUM + V(I,K)*V(J,K)
   50         CONTINUE
              IJ = IA(I) + J
              DB(IJ) = DUM
C       4-----
   60       CONTINUE
C     3-----
   70     CONTINUE
C   2-----
        ENDIF
C 1-----
      ENDIF
C
      NZFP = NEXP*4
      NZFN = NEXP*4 + NGRAN + 1
C
      dN(1) = zero
      dN(2) = zero
      dN(3) = zero
      dN(4) = zero
      dN(kxyz) = one
c
      IJ = 0
C
C-----  Loop over charge distributions
C
      DO I = 1, NUM
        DO J = 1, I
C   1-----
          IJ = IJ + 1
C
c         IF (IJBIT(IJ) .EQ. 0) GOTO 100
C
          FACTdij = TWO
          FACTIJ = TWO
          IF (I.EQ.J) FACTIJ = ONE
C
C    -----  Get expansion centre
C           of IJ-th charge distribution
C
          IJEXP = IEXPC(IJ)
          JF = (IJEXP-1)*4
C
C    -----  Collect overlap, dipole
C
          P(1) = DXd(IJ)
          P(2) = DYd(IJ)
          P(3) = DZd(IJ)
          P(4) = ds(IJ)
C
          sP(1) = DX(IJ)
          sP(2) = DY(IJ)
          sP(3) = DZ(IJ)
          sP(4) = OL(IJ)
C    -----  Interaction with externally induced dipoles
C
          DO IGR = 1, NGRAN
            DO L = 1, 4
C       2-----
C        -----  NOTE: --sign because coupling is to 
C               electrons
C
C        -----  Electrons are source
C
              GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm)
     1     - pt5*FACTdIJ*DA(IJ)*OMEGA(NZFP+IGR,JF+L)*P(L)
     1     - pt5*FACTdIJ*DA(IJ)*OMEGA(jf+l,NZFP+IGR)*P(L)
C       2-----
            ENDDO
          ENDDO
C
C    -----  Screening of the nuclear attraction
C
          DO L = 1, 4
C     2-----
C      -----  NOTE: --sign because coupling is to 
C             electrons
C
C      -----  Electrons are source
C
            GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm)
     1     - pt5*FACTdIJ*DA(IJ)*OMEGA(NZFN,JF+L)*P(L)
     1     - pt5*FACTdIJ*DA(IJ)*OMEGA(jf+l,NZFN)*P(L)
C     2-----
          ENDDO
C
              CALL DRFOAB(IJEXP,Iatm,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,dN,OdN,4,.FALSE.)
C
              CALL DRFOAB(Iatm,IJEXP,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,sP,OsP,4,.FALSE.)
C
            GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm)
     1     - pt5*FACTIJ*czan(iatm)*DA(IJ)*ddot(4,OdN,1,sP,1)
     1     - pt5*FACTIJ*czan(iatm)*DA(IJ)*ddot(4,OsP,1,dN,1)
c    1     - pt5*FACTIJ*czan(iatm)*DA(IJ)*adotb(OdN,sP,4)
c    1     - pt5*FACTIJ*czan(iatm)*DA(IJ)*adotb(OsP,dN,4)
C
C    -----  Self energy; dispersion
C
          IF (GAMDRF .NE. ZERO) THEN
C     2-----
            call drfoab(ijexp,ijexp,nwtc,omega)
            rr(1,1) = two*rxx(ij)
            rr(1,2) = two*rxy(ij)
            rr(2,1) = two*ryx(ij)
            rr(1,3) = two*rxz(ij)
            rr(3,1) = two*rzx(ij)
            rr(1,4) = two*dxd(ij)
            rr(4,1) = two*dxd(ij)
            rr(2,2) = two*ryy(ij)
            rr(2,3) = two*ryz(ij)
            rr(3,2) = two*rzy(ij)
            rr(2,4) = two*dyd(ij)
            rr(4,2) = two*dyd(ij)
            rr(3,3) = two*rzz(ij)
            rr(3,4) = two*dzd(ij)
            rr(4,3) = two*dzd(ij)
            rr(4,4) = two*ds(ij)
            contr = ddot(16,rr,1,omgab,1)*da(ij)*factdij*pt25
c           contr = adotb(rr,omgab,16)*da(ij)*factdij*pt25
            grad(kxyz,iatm) = grad(kxyz,iatm) + contr*gamdrf
C     2-----
          ENDIF
C
C    -----  Two-electron interaction
C
          DO K = 1, I 
            LMAX = K
            IF (K .EQ. I) LMAX = J
            DO L = 1, LMAX
C       2-----
              KL = IA(max(K,l)) + min(L,k)
              KLEXP = IEXPC(KL)
C
C        -----  Collect overlap, dipole
C
              q(1) = DX(KL)
              q(2) = DY(KL)
              q(3) = DZ(KL)
              q(4) = OL(KL)
C
              dQ(1) = DXd(KL)
              dQ(2) = DYd(KL)
              dQ(3) = DZd(KL)
              dQ(4) = ds(KL)
C
              CALL DRFOAB(IJEXP,KLEXP,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,q,OQ,4,.FALSE.)
C
              CALL MATVEC(OMGAB,dQ,OdQ,4,.FALSE.)
C
              CALL DRFOAB(KLEXP,IJEXP,NWTC,OMEGA)
C
              CALL MATVEC(OMGAB,P,OP,4,.FALSE.)
C
              CALL MATVEC(OMGAB,sP,OsP,4,.FALSE.)
C
C        -----  Calculate addresses
C
              ik = ia(max(i,k)) + min(i,k)
              jk = ia(max(j,k)) + min(j,k)
              il = ia(max(i,l)) + min(i,l)
              jl = ia(max(j,l)) + min(j,l)
C
              II = I
              JJ = J
              KK = K
              LL = L
C
C        -----  Calculate density factors
C
              CALL DRFDAB(DA,DB,V,DCOUL,DEXCH,NORB,ROHF,RGVB,UHF)
C
              IF ((ABS(DCOUL) .GT. THRESH2) .OR.
     1            (ABS(DEXCH) .GT. THRESH2)) THEN
C         5-----
                FACTKL = TWO
                IF (L .EQ. K) FACTKL = one 
                factdkl = two
		fact = one
                IF ((I .EQ. K) .AND. (J .EQ. L)) FACT = PT5*FACT
		fact = fact*pt5
C
C          -----  Factor 0.25 for double counting twice
C
                GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm) +
     1             (DCOUL+DEXCH)*factdij*factkl*FACT*PT25*
     2             ddot(4,OP,1,q,1)
c    2             ADOTB(OP,q,4)
                GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm) +
     1             (DCOUL+DEXCH)*factdij*factkl*FACT*PT25*
     2             ddot(4,OQ,1,P,1)
c    2             ADOTB(OQ,P,4)
                GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm) +
     1        (DCOUL+DEXCH)*factij*factdkl*FACT*PT25*
     2             ddot(4,OsP,1,dQ,1)
c    2             ADOTB(OsP,dQ,4)
                GRAD(KXYZ,Iatm) = GRAD(KXYZ,Iatm) +
     1        (DCOUL+DEXCH)*factij*factdkl*FACT*PT25*
     2             ddot(4,OdQ,1,sP,1)
c    2             ADOTB(OdQ,sP,4)
C         5-----
              ENDIF
C
 150          CONTINUE
C       2-----
            ENDDO
          ENDDO
C
 100      CONTINUE
C   1-----
        ENDDO
      ENDDO
      if (idrfout .ge. 3)
     1  call hatout(grad,3,nexp,2,'grax2')
C-----
      RETURN
      END
      SUBROUTINE DISFOR(GRAD,q,NCEN,NQ)
C------
C      Distribute force on centre of charge
C      to charge-carriers only
C------
      IMPLICIT real*8 (A-H,O-Z)
C
      integer ncen, nq
      DIMENSION GRAD(3,NCEN)
      DIMENSION q(NQ)
C     
      QTOT = 0.D00
      DO I = 1, NQ
        QTOT = QTOT + q(I)
      ENDDO
C
      IF (QTOT .EQ. 0.) THEN
        STOP 'TOTAL NUCLEAR CHARGE ZERO'
      ENDIF
C
      DO I = 1, NQ
        DO J = 1, 3
          GRAD(J,I) = GRAD(J,I) + (q(I)/QTOT)*GRAD(J,NCEN)
        ENDDO
      ENDDO
C
      RETURN
      END
