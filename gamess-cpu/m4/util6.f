c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util6.m,v $
c  $State: Exp $
c  
      subroutine ijconr(s,t,nx,mapnr)
      implicit real*8  (a-h,o-z)
c     condense lower triangle to n.r. members
c
      dimension s(*),t(*),mapnr(*)
c
c     nx = length of triangle
c
      do 20 i = 1 , nx
         j = mapnr(i)
         if (j.ne.0) then
c
c     s = input : t =output
c
            t(j) = s(i)
         end if
 20   continue
      return
      end
      subroutine grhfbl(scftyp)
c     set up grhf block for closed and oscf types
      implicit real*8  (a-h,o-z)
      character *8 scftyp,oscf,grhf,open
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
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
      dimension i11(11)
      data i11/0,11,22,33,44,55,66,77,88,99,110/
      data oscf,grhf,open/'oscf','grhf','open'/
      if (scftyp.eq.grhf .or. scftyp.eq.open) return
c
c
      nact = num
      do 20 i = 1 , num
         iactiv(i) = i
 20   continue
      nbshel(1) = na
      ilfshl(1) = 0
      ilfshl(2) = na
      do 40 i = 1 , 11
         fjk(i) = 0.0d0
         fcan(i) = 0.0d0
         do 30 j = 1 , 11
            erga(i11(i)+j) = 0.0d0
            cana(i11(i)+j) = 0.0d0
            ergb(i11(i)+j) = 0.0d0
            canb(i11(i)+j) = 0.0d0
 30      continue
 40   continue
      fjk(1) = 2.0d0
      erga(i11(1)+1) = 2.0d0
      ergb(i11(1)+1) = -1.0d0
      cana(i11(1)+1) = 2.0d0
      cana(i11(2)+1) = 2.0d0
      canb(i11(1)+1) = -1.0d0
      canb(i11(2)+1) = -1.0d0
      fcan(1) = 2.0d0
      fcan(2) = 2.0d0
      if (scftyp.eq.oscf) then
         njk = 2
         njk1 = 3
         fjk(2) = 1.0d0
         fjk(3) = 0.0d0
         nbshel(2) = nb - na
         nbshel(3) = num - nb
         ilfshl(3) = nb
         erga(i11(1)+2) = 1.0d0
         erga(i11(2)+1) = 1.0d0
         erga(i11(2)+2) = 0.5d0
         erga(i11(1)+3) = 0.0d0
         erga(i11(3)+1) = 0.0d0
         erga(i11(2)+3) = 0.0d0
         erga(i11(3)+2) = 0.0d0
         erga(i11(3)+3) = 0.0d0
         ergb(i11(1)+2) = -0.5d0
         ergb(i11(2)+1) = -0.5d0
         ergb(i11(2)+2) = -0.5d0
         ergb(i11(3)+1) = 0.0d0
         ergb(i11(1)+3) = 0.0d0
         ergb(i11(2)+3) = 0.0d0
         ergb(i11(3)+2) = 0.0d0
         ergb(i11(3)+3) = 0.0d0
         return
      else
c
c     closed
c
         njk = 1
         njk1 = 2
         fjk(2) = 0.0d0
         nbshel(2) = num - na
         erga(i11(1)+2) = 0.0d0
         erga(i11(2)+1) = 0.0d0
         erga(i11(2)+2) = 0.0d0
         ergb(i11(1)+2) = 0.0d0
         ergb(i11(2)+1) = 0.0d0
         ergb(i11(2)+2) = 0.0d0
         return
      end if
      end
      subroutine dr2lag(dd,v,e,ndim)
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
c     extracts the lagrangian for various types of wavefunction
c
      dimension dd(*),e(*),v(ndim,*)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
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
      logical exist
      character *8 open,oscf,grhf
      data zer,one,two /0.0d0,1.0d0,2.0d0/
      data open/'open'/
      data oscf /'oscf'/
      data grhf/'grhf'/
c
      if (lci .or. lmcscf .or. cigr .or. mcgr) then
         mtyp = 0
         call secloc(isecll,exist,iblok)
         if (exist) then
            call rdedx(dd,nx,iblok,ifild)
         else
            call caserr('lagrangian not found')
         end if
         return
      else if (scftyp.eq.grhf .or. scftyp.eq.oscf) then
         call secloc(isect(42),exist,iblok)
         if (exist) then
            call rdedx(dd,lds(isect(42)),iblok,ifild)
         else
            call caserr('lagrangian not found')
         end if
         return
      else
         if (mp2 .or. mp3) then
            mtyp = 0
            call secloc(isecll,exist,iblok)
            if (exist) then
               if (odebug(30). and. nprint.ne.-5)  then
                if (mp2 ) write (iwr,6010) isecll
                if (mp3 ) write (iwr,6020) isecll
               endif
               call rdedx(dd,nx,iblok,ifild)
            else
               if (mp2) call caserr('mp2 lagrangian not found')
               if (mp3) call caserr('mp3 lagrangian not found')
            end if
         else
            call vclr(dd,1,nx)
         end if
         occ = one
         if (scftyp.ne.open) occ = two
         mtyp = 0
         call secloc(isect(8),exist,iblok)
         if (exist) then
            iblvec = iblok + mvadd
            call rdedx(v,num*ncoorb,iblvec,ifild)
         else
            call caserr('vectors not found')
         end if
         mtyp = 0
         call secloc(isect(9),exist,iblok)
         if (exist) then
            call rdedx(e,lds(isect(9)),iblok,ifild)
         else
            call caserr('eigenvalues not found')
         end if
         ij = 0
         do 40 i = 1 , num
            do 30 j = 1 , i
               ij = ij + 1
               dum = zer
               do 20 k = 1 , na
                  dum = dum - e(k)*v(i,k)*v(j,k)
 20            continue
               dd(ij) = dd(ij) + dum*occ
 30         continue
 40      continue
         if (scftyp.ne.open) return
         if (nb.eq.0) return
         mtyp = 0
         call secget(isect(11),mtyp,iblok)
         iblvec = iblok + mvadd
         call rdedx(v,num*ncoorb,iblvec,ifild)
         mtyp = 0
         call secget(isect(12),mtyp,iblok)
         call rdedx(e,lds(isect(12)),iblok,ifild)
         ij = 0
         do 70 i = 1 , num
            do 60 j = 1 , i
               ij = ij + 1
               dum = zer
               do 50 k = 1 , nb
                  dum = dum - e(k)*v(i,k)*v(j,k)
 50            continue
               dd(ij) = dd(ij) + dum*occ
 60         continue
 70      continue
         return
      end if
 6010 format (/1x,'mp2 lagrangian at section ',i4)
 6020 format (/1x,'mp3 lagrangian at section ',i4)
      end
      subroutine nrmapo(mnr,nocca,noccb,nsa,iky,mapie)
      implicit real*8  (a-h,o-z)
      dimension mnr(*),mapie(*),iky(*)
c
c     nsa is total number of active m.o.'s
c     nocca is number of doubly occupied
c     noccb is total occupied
c
c     maps from lower triangle to non-redundant pairs
c
      nmax = mapie(nsa)
      lenn = iky(nmax+1)
      do 20 i = 1 , lenn
         mnr(i) = 0
 20   continue
      nsoc = noccb - nocca
      nvirta = nsa - noccb
      ntpls1 = noccb + 1
      ndpls1 = nocca + 1
      if (nocca.ne.0) then
         do 40 j = 1 , nocca
            do 30 i = ndpls1 , nsa
               it = (i-ndpls1)*nocca + j
               ij = iky(mapie(i)) + mapie(j)
               mnr(ij) = it
 30         continue
 40      continue
      end if
      if (noccb.ne.nocca) then
         do 60 j = ndpls1 , noccb
            do 50 i = ntpls1 , nsa
               it = nvirta*nocca + (i-nsoc-1)*nsoc + j - nocca
               ij = iky(mapie(i)) + mapie(j)
               mnr(ij) = it
 50         continue
 60      continue
      end if
      return
      end
      subroutine qhq1(a,q,ilifq,ncore,h,iky,nbasis)
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
      dimension a(*),q(*),h(*),ilifq(*),iky(*)
      common/small/p(maxorb),y(maxorb)
c...   a=q(transpose) * h * q
c...   a and h stored in triangle form
      m1 = ilifq(1)
      nb = ilifq(2) - m1
      m1 = m1 + 1
      do 30 j = 1 , ncore
         map = ilifq(j)
         do 20 i = 1 , nbasis
            m = iky(i) + 1
            ilen = i
            p(i) = ddot(ilen,q(map+1),1,h(m),1)
            iless1 = i - 1
            call daxpy(iless1,q(map+i),h(m),1,p,1)
 20      continue
         call mxma(q(m1),nb,1,p,1,nbasis,a(iky(j)+1),1,j,j,nbasis,1)
 30   continue
      return
      end
      subroutine rdedv (v,n,nv,ibl,idev)
      implicit real*8  (a-h,o-z)
      dimension v(n,nv)
      call search(ibl,idev)
      call find(idev)
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      call get(v(i,iv),nw)
         i = i + nw
         k = k - nw
         if (k.gt.0 .or. iv.lt.nv) call find(idev)
         if (k.lt.0) go to 40
         if (k.ne.0) go to 20
 30   continue
      return
 40   call caserr('error in rdedv')
      end
      subroutine rdedvs(v,n,nv,idev)
      implicit real*8  (a-h,o-z)
      dimension v(n,nv)
      call find(idev)
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      call get(v(i,iv),nw)
         i = i + nw
         k = k - nw
         if (k.gt.0 .or. iv.lt.nv) call find(idev)
         if (k.lt.0) go to 40
         if (k.ne.0) go to 20
 30   continue
      return
 40   call caserr('error in rdedvs')
      end
      subroutine squars(t,sq,n)
      implicit real*8  (a-h,o-z)
c     triangle to square
      dimension sq(*),t(*)
      ij = 0
      ii = 0
      do 30 i = 1 , n
         jj = 0
         do 20 j = 1 , i
            ij = ij + 1
            sq(ii+j) = -t(ij)
            sq(jj+i) = t(ij)
            jj = jj + n
 20      continue
         ii = ii + n
 30   continue
      return
      end
      subroutine stvin2
c
c     ----- gauss-hermite quadrature using minimum point formula -----
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
      common/junk/xyzsp(3,maxat),
     * xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h(28)
      common/wermit/w(28)
      dimension min(7),max(7)
      data min /1,2,4,7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero,one /0.0d0,1.0d0/
      xint = zero
      yint = zero
      zint = zero
      npts = (ni+nj)/2
      imin = min(npts)
      imax = max(npts)
      do 160 i = imin , imax
         px = one
         py = one
         pz = one
         dum = h(i)/t
         ptx = dum + x0
         pty = dum + y0
         ptz = dum + z0
         ax = ptx - xi
         ay = pty - yi
         az = ptz - zi
         bx = ptx - xj
         by = pty - yj
         bz = ptz - zj
         go to (80,70,60,50,40,30,20) , ni
 20      px = ax
         py = ay
         pz = az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      px = px*ax
         py = py*ay
         pz = pz*az
 70      px = px*ax
         py = py*ay
         pz = pz*az
 80      go to (150,140,130,120,110,100,90) , nj
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     px = px*bx
         py = py*by
         pz = pz*bz
 140     px = px*bx
         py = py*by
         pz = pz*bz
 150     dum = w(i)
         xint = xint + dum*px
         yint = yint + dum*py
         zint = zint + dum*pz
 160  continue
      return
      end
      subroutine clenms(messge)
c
c
      implicit real*8  (a-h,o-z)
      character messge*(*)
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
c
      l = len(messge)
      write (iwr,'(/1x,a)') messge(1:l)
      call clenup
      stop
      end
      subroutine wrt3z(iblk,num3,len)
      implicit real*8  (a-h,o-z)
      common/blkin/q(511)
      call search(iblk,num3)
      do 20 i = 1 , len
         call vclr(q,1,511)
         call put(q,511,num3)
 20   continue
      call search(iblk,num3)
      return
      end
      subroutine vwrt3s(v,n,nv,idev)
      implicit real*8  (a-h,o-z)
      dimension v(n,nv)
      data ilen/511/
      do 30 iv = 1 , nv
         i = 1
         k = n
 20      nw = min(k,ilen)
         call put(v(i,iv),nw,idev)
         i = i + nw
         k = k - nw
         if (k.gt.0) go to 20
 30   continue
      return
      end
      subroutine daxpyi(n,sa,sx,indx,sy)
      implicit real*8  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension sx(*),sy(*),indx(*)
      if(n.le.0.or.sa.eq.0.d0) return
      m = n - (n/4)*4
      if( m .eq. 0 ) go to  4
      do  3 i = 1,m
        sy(indx(i)) = sy(indx(i)) + sa*sx(i)
    3 continue
      if( n .lt. 4 ) return
    4 mp1 = m + 1
      do  8 i = mp1,n,4
        sy(indx(i)) = sy(indx(i)) + sa*sx(i)
        sy(indx(i+1)) = sy(indx(i+1)) + sa*sx(i + 1)
        sy(indx(i+2)) = sy(indx(i+2)) + sa*sx(i + 2)
        sy(indx(i+3)) = sy(indx(i+3)) + sa*sx(i + 3)
    8 continue
      return
      end
      subroutine dlstgt(n,x,incx,scalar,nindx,indx)
      implicit real*8  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).le.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstne(n,x,incx,scalar,nindx,indx)
      implicit real*8  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).eq.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstge(n,x,incx,scalar,nindx,indx)
      implicit real*8  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).lt.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine dlstlt(n,x,incx,scalar,nindx,indx)
      implicit real*8  (a-h,o-z),integer    (i-n)
      dimension x(*),indx(*)
      nindx=0
      ij=1
      do 1 loop=1,n
      if(x(ij).ge.scalar) go to 1
      nindx=nindx+1
      indx(nindx)=ij
 1    ij=ij+incx
      return
      end
      subroutine ver_util6(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util6.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
