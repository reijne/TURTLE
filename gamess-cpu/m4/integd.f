c 
c  $Author: mrdj $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integd.m,v $
c  $State: Exp $
c  
      subroutine dmder(q)
c
c---------------------------------------------------------------
c     dipole and quadrupole derivatives
c----------------------------------------------------------------
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
      dimension q(*)
      common/small/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
c
c     get dipole and quadrupole for general interest
c
      call dipmom(q)
c
      call dmdint(q,dipi,dipn,dipd)
      call dmdwfn(q,dipd)
      return
      end
      subroutine dipmom(q)
c
c-------------------------------------------------------------
c     dipole moments
c-------------------------------------------------------------
c
      implicit real*8  (a-h,o-z)
      logical mpir
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
      logical lmos,lnucl
      common/small/y(maxorb)
      common/prpsec/isecd1,itypd1,isecd2,itypd2,lmos,isecv,itypv
     1   ,isecd3,itypd3,lnucl
      common/maxlen/maxq
      dimension q(*),origin(3)
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
      dimension de(3),dn(3),dt(3)
      character *8 grhf,closed,dipmp
      data grhf/'grhf'/
      data dipmp/'mpdipder'/
      data zero,fac,closed/0.0d0,2.5417701d0,'closed'/
c       if(mp2.and.scftyp.eq.closed) call mp2nat
      if (ich.ne.0 .and. oprn(6)) write (iwr,6010)
      origin(1) = gx
      origin(2) = gy
      origin(3) = gz
      mpir = runtyp.eq.dipmp
c     ----- calculate dipole moment integrals -----
      if(oprn(21).and.oprn(40)) write (iwr,6020) origin
      i1 = nx + 1
      i2 = nx + i1
      i3 = nx + i2
      if (isecd1.eq.0) isecd1 = isect(7)
      if (lmcscf .or. lci .or. mpir .or. mp2 .or. mp3) then
         call onepdm(q(1),q(i1))
         go to 30
      else
         call secget(isecd1,itypd1,iblok)
         call rdedx(q(1),nx,iblok,ifild)
      end if
      if (scftyp.ne.grhf .and. scftyp.ne.closed) then
         if (isecd2.eq.0) isecd2 = isect(10)
         call secget(isecd2,itypd2,iblok)
         call rdedx(q(i1),nx,iblok,ifild)
         do 20 i = 1 , nx
            q(i) = q(i) + q(i1+i-1)
 20      continue
      end if
 30   if (lmos) then
         call secget(isecv,itypv,iblok)
         call rdedx(q(i2),num*ncoorb,iblok+mvadd,ifild)
         call demoao(q(1),q(i1),q(i2),y,num,num,num)
         call dcopy(nx,q(i1),1,q(1),1)
      end if
      call dmints(q(i1),q(i2),q(i3))
c     ----- electronic contribution to dipole moment -----
      de(1) = -tracep(q(1),q(i1),num)
      de(2) = -tracep(q(1),q(i2),num)
      de(3) = -tracep(q(1),q(i3),num)
c     ----- nuclear contribution
      do 50 j = 1 , 3
         dn(j) = zero
         if (lnucl) then
            do 40 i = 1 , nat
               dn(j) = dn(j) + czan(i)*(c(j,i)-origin(j))
 40         continue
         end if
 50   continue
      do 60 j = 1 , 3
         dt(j) = de(j) + dn(j)
 60   continue
      dipol = dsqrt(dt(1)*dt(1)+dt(2)*dt(2)+dt(3)*dt(3))
      if(oprn(21).and.oprn(40)) then
       write (iwr,6030)
       write (iwr,6040) de , dn , dt , dipol
      endif
      dipol = fac*dipol
      do 70 i = 1 , 3
         de(i) = de(i)*fac
         dn(i) = dn(i)*fac
         dt(i) = dt(i)*fac
 70   continue
      if(oprn(21).and.oprn(40)) then
       write (iwr,6050)
       write (iwr,6040) de , dn , dt , dipol
      endif
      return
 6010 format (//10x,'warning ---- dipole moment for charged species',
     +        ' is origin dependent')
 6020 format (/1x,'origin for dipole moment evaluation',3f12.6)
 6030 format (/
     +      40x,'***********************'/
     +      40x,'*    dipole moment    *'/
     +      40x,'*   (atomic units )   *'/
     +      40x,'***********************'/)
 6040 format (/30x,'x',15x,'y',15x,'z'//10x,'electronic',3(2x,f14.8)
     +        /10x,'nuclear',3x,3(2x,f14.8)/10x,'total',5x,3(2x,f14.8)
     +        //10x,'dipole moment',f15.8)
 6050 format (/
     +      40x,'***********************'/
     +      40x,'*    dipole moment    *'/
     +      40x,'*      (debyes)       *'/
     +      40x,'***********************'/)
      end
      subroutine dmdrot(r,dipd,nat)
c
c      rotational sum rules for dipole derivatives
c
      implicit real*8  (a-h,o-z)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension r(3,nat),dipd(3,3,nat),rot(9)
      character*2 head
      dimension head(9)
      data head/'xx','xy','xz','yx','yy','yz','zx','zy','zz'/
      n = 0
      do 40 ialp = 1 , 3
         do 30 ibet = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 20 i = 1 , nat
               rot(n) = rot(n) + dipd(ibet,ialp,i)
 20         continue
 30      continue
 40   continue
      if (oprn(40) .and. oprn(25) ) then
       write (iwr,6010) (head(i),i=1,9)
       write (iwr,6030) (rot(i),i=1,9)
      endif
      n = 0
      do 90 ibet = 1 , 3
         do 80 idel = 1 , 3
            n = n + 1
            rot(n) = 0.0d0
            do 70 ialp = 1 , 3
               do 60 igam = 1 , 3
                  skew = etijk(ialp,igam,idel)
                  if (skew.ne.0.0d0) then
                     do 50 i = 1 , nat
                        rot(n) = rot(n) + skew*r(igam,i)
     +                           *dipd(ibet,ialp,i)
 50                  continue
                  end if
 60            continue
 70         continue
 80      continue
 90   continue
      if (oprn(40) .and. oprn(25) ) then
       write (iwr,6020) (head(i),i=1,9)
       write (iwr,6030) (rot(i),i=1,9)
      endif
      return
 6010 format (/
     +20x,'**********************************'/
     +20x,'* dipole translational sum rules *'/
     +20x,'**********************************'/
     +   /1x,9(10x,a2))
 6020 format (//
     +20x,'*******************************'/
     +20x,'* dipole rotational sum rules *'/
     +20x,'*******************************'/
     +   /1x,9(10x,a2))
 6030 format (1x,9f12.8)
      end
      subroutine dmdsym(dipd,iso,ict,natdim,nshels)
c----------------------------------------------------------
c     symmetrizes dipole moment derivatives
c----------------------------------------------------------
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
      common/scrtch/
     + pptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     + dp(3,3,maxat),qd(6,3,maxat)
      common/junkp/ptr(3,144),djunkp(6,288)
c
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
      dimension dipd(3,3,maxat)
      dimension iso(nshels,*),ict(natdim,*)
      if (nt.eq.1) return
      call rdedx(ptr,nw196(1),ibl196(1),ifild)
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      zero = 0.0d0
c     one = 1.0d0
      do 40 nd = 1 , 3
         do 30 nc = 1 , 3
            do 20 n = 1 , nat
               dp(nd,nc,n) = zero
 20         continue
 30      continue
 40   continue
c
c     ----- set transformation table: atoms versus symmetry operations.
c
      do 70 ii = 1 , nshell
         ic = katom(ii)
         do 60 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 60      continue
 70   continue
c
c
      do 130 n = 1 , nat
         do 120 nd = 1 , 3
            do 110 nc = 1 , 3
               do 100 nop = 1 , nt
                  icnu = ict(n,nop)
                  nn = 3*(nop-1)
                  do 90 ndd = 1 , 3
                     do 80 ncd = 1 , 3
                        dp(nd,nc,n) = dp(nd,nc,n) + dipd(ndd,ncd,icnu)
     +                                *ptr(ndd,nn+nd)*ptr(ncd,nn+nc)
 80                  continue
 90               continue
 100           continue
 110        continue
 120     continue
 130  continue
      dum = dble(nt)
      do 160 n = 1 , nat
         do 150 i = 1 , 3
            do 140 j = 1 , 3
               dipd(i,j,n) = dp(i,j,n)/dum
 140        continue
 150     continue
 160  continue
      return
      end
      subroutine dmdint(dd,dipi,dipn,dipd)
c
c     dipole moment derivatives
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
      logical out,norm
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     4   ijx(100),ijy(100),ijz(100)
      dimension dd(*)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension comp(3),dipd(3,3,maxat),dipi(3,3,maxat),
     &  dipn(3,3,maxat)
c
      character*5 comp
      data comp/'d/dx','d/dy','d/dz'/
      data ndim/5/
      data one /1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
      call onepdm(dd,dd(nx+1))
      tol = 2.30258d0*itol
      out = odebug(15)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
      do 50 n = 1 , nat
         do 30 i = 1 , 3
            do 20 j = 1 , 3
               dipn(i,j,n) = 0.0d0
               dipi(i,j,n) = 0.0d0
               dipd(i,j,n) = 0.0d0
 20         continue
 30      continue
         do 40 i = 1 , 3
            dipn(i,i,n) = czan(n)
 40      continue
 50   continue
c     ----- ishell
      do 130 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 120 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +           ndim,1,1)
            dxx = 0.0d0
            dyy = 0.0d0
            dzz = 0.0d0
            dxy = 0.0d0
            dxz = 0.0d0
            dyx = 0.0d0
            dyz = 0.0d0
            dzx = 0.0d0
            dzy = 0.0d0
c     ----- i primitive
            do 110 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 100 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                           csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 70 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 60 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
 60                     continue
 70                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 90 i = mini , maxi
                        in = loci + i
                        do 80 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
                           dum = dd(nn)*dij(n)
                           dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
                           dxx = dxx + dum*xdipd(nnx)*yin(ny)*zin(nz)
                           dyy = dyy + dum*xin(nnx)*ydipd(ny)*zin(nz)
                           dzz = dzz + dum*xin(nnx)*yin(ny)*zdipd(nz)
                           dxy = dxy + dum*xdip(nnx)*yd(ny)*zin(nz)
                           dyx = dyx + dum*xd(nnx)*ydip(ny)*zin(nz)
                           dxz = dxz + dum*xdip(nnx)*yin(ny)*zd(nz)
                           dzx = dzx + dum*xd(nnx)*yin(ny)*zdip(nz)
                           dyz = dyz + dum*xin(nnx)*ydip(ny)*zd(nz)
                           dzy = dzy + dum*xin(nnx)*yd(ny)*zdip(nz)
 80                     continue
 90                  continue
                  end if
 100           continue
 110        continue
            dipi(1,1,iat) = dipi(1,1,iat) - dxx
            dipi(2,2,iat) = dipi(2,2,iat) - dyy
            dipi(3,3,iat) = dipi(3,3,iat) - dzz
            dipi(1,2,iat) = dipi(1,2,iat) - dxy
            dipi(1,3,iat) = dipi(1,3,iat) - dxz
            dipi(2,1,iat) = dipi(2,1,iat) - dyx
            dipi(3,1,iat) = dipi(3,1,iat) - dzx
            dipi(2,3,iat) = dipi(2,3,iat) - dyz
            dipi(3,2,iat) = dipi(3,2,iat) - dzy
 120     continue
 130  continue
c
      if (out) then
         write (iwr,6010)
         do 150 n = 1 , nat
            write (iwr,6020)
            do 140 nc = 1 , 3
               write (iwr,6030) zaname(n) , comp(nc) ,
     +                         (dipi(nn,nc,n),nn=1,3)
 140        continue
 150     continue
      end if
      do 180 i = 1 , 3
         do 170 j = 1 , 3
            do 160 k = 1 , nat
               dipd(i,j,k) = dipd(i,j,k) + dipn(i,j,k) + dipi(i,j,k)
 160        continue
 170     continue
 180  continue
      return
 6010 format (//35x,'integral derivative contribution'//30x,'x',15x,'y',
     +        15x,'z',/)
 6020 format (//)
 6030 format (5x,a8,5x,a5,3f16.8)
      end
      subroutine dmdwfn(dd,dipd)
c     subroutine dmdwfn(dd,dipi,dipn,dipd)
c-----------------------------------------------------------------
c      wavefunction derivative contribution to dipole derivative
c-------------------------------------------------------------------
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),orb(3,3,maxat)
      dimension dd(*)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      dimension comp(3),dipd(3,3,maxat)
c     dimension dipi(3,3,maxat),dipn(3,3,maxat)
c
      character*5 comp
c     character *8 closed
c     data closed/'closed'/
      data comp/'d/dx','d/dy','d/dz'/
c
      ioff = 1
      i1 = nx + 1
      i2 = i1 + nx
      i3 = i2 + nx
      i4 = i3 + nx
c     i5 = i4 + nx
      nat3 = nat*3
      lennew = iky(ncoorb+1)
      if (lmeth2) then
         write (iwr,*) 'dipole derivatives from alternative formula'
         call dmder0(dd(1),dd(i1),dd(i2),dd(i3),dd(i4),dipd,lennew,nat3
     +  )
      else
         icomp = 21
         do 20 i = 1 , 3
            ioff = ioff + nx
            icomp = icomp + 1
            call vclr(dd(ioff),1,lennew)
            if (ione(i+3).eq.0) then
               write (iwr,6010) i
            else
               call secget(isect(icomp),icomp,iblok)
               call rdedx(dd(ioff),lennew,iblok,ifild)
            end if
 20      continue
         iblls = lensec(lennew)
         iposs = iochf(15)
         do 40 n = 1 , nat
            do 30 nc = 1 , 3
               call rdedx(dd,lennew,iposs,ifockf)
               iposs = iposs + iblls
               dumx = tracep(dd(i1),dd,ncoorb)
               dumy = tracep(dd(i2),dd,ncoorb)
               dumz = tracep(dd(i3),dd,ncoorb)
               dipd(1,nc,n) = dipd(1,nc,n) - dumx
               dipd(2,nc,n) = dipd(2,nc,n) - dumy
               dipd(3,nc,n) = dipd(3,nc,n) - dumz
 30         continue
 40      continue
      end if
      call dmdsym(dipd,dd,dd(nw196(5)+1),nat,nshell)
      if (oprn(40)) then
       if (opunch(3)) write (ipu,6060)
       if (oprn(25)) write (iwr,6020)
       do 60 n = 1 , nat
       if (oprn(25)) write (iwr,6030)
          do 50 nc = 1 , 3
             if (opunch(3)) write (ipu,6050) (dipd(nn,nc,n),nn=1,3)
             if(oprn(25))
     +  write (iwr,6040) zaname(n) , comp(nc) , (dipd(nn,nc,n),nn=1,3)
 50       continue
 60    continue
       call dmdrot(c,dipd,nat)
      endif
      return
 6010 format (//10x,'dipole integrals missing for component',i6)
 6020 format (/
     + 20x,'****************************'/
     + 20x,'*  scf dipole derivatives  *'/
     + 20x,'*      (atomic units)      *'/
     + 20x,'****************************'//
     + 30x,'x',15x,'y',15x,'z')
 6030 format (/)
 6040 format (5x,a8,5x,a5,3f16.8)
 6050 format (1x,3e20.12)
 6060 format (1x,'scf dipole derivatives')
      end
      subroutine dmder0(da,db,dc,fa,e,dipd,ltri,nat3)
c--------------------------------------------------------------
c     another way of doing dipole derivatives
c     (used for testing only)
c---------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension e(*)
      dimension da(ltri),db(ltri),dc(ltri),fa(ltri),dipd(3,nat3)
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
      m = 0
      call secget(isect(31),m,iblok)
      call search(iblok,ifild)
      call reads(da,ltri,ifild)
      call reads(db,ltri,ifild)
      call reads(dc,ltri,ifild)
c
      call search(iochf(13),ifockf)
      do 20 n = 1 , nat3
         call reads(fa,ltri,ifockf)
         dipd(1,n) = dipd(1,n) - tracep(fa,da,ncoorb)
         dipd(2,n) = dipd(2,n) - tracep(fa,db,ncoorb)
         dipd(3,n) = dipd(3,n) - tracep(fa,dc,ncoorb)
 20   continue
c
      call secget(isect(9),9,iblok)
      call rdedx(e,ncoorb,iblok,ifild)
      do 80 k = 1 , 3
         call rdedx(da,ltri,iochf(k+1),ifockf)
         m = 0
         call secget(isect(30+k),m,iblok)
         call rdedx(db,ltri,iblok,ifild)
         call vclr(dc,1,ltri)
         ij = 0
         do 40 i = 1 , nocc
            do 30 j = 1 , i
               ij = ij + 1
               dc(ij) = 2.0d0*da(ij)
 30         continue
 40      continue
         do 60 i = nocc + 1 , ncoorb
            do 50 j = 1 , nocc
               ij = iky(i) + j
               dc(ij) = db(ij)*e(j)
 50         continue
 60      continue
         call search(iochf(14),ifockf)
         do 70 n = 1 , nat3
            call reads(fa,ltri,ifockf)
            dipd(k,n) = dipd(k,n) + tracep(dc,fa,ncoorb)
 70      continue
 80   continue
      return
      end
      subroutine qmmo(q)
c------------------------------------------------------------------
c     quadrupole integrals into m.o. basis
c-----------------------------------------------------------------
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
      common/maxlen/maxq
      dimension q(*)
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
      data half,thrhal/0.5d0,1.5d0/
      nword = nx
      i1 = nx + 1
      i2 = i1 + nx
      i3 = i2 + nx
      i4 = i3 + nx
      i5 = i4 + nx
      i6 = i5 + nx
      if (i6.gt.maxq) call caserr('insufficient core')
      len = lensec(nword)
      ibl = iblks
      call qmints(q(1),q(i1),q(i2),q(i3),q(i4),q(i5))
      do 40 j = 1 , nword
         do 20 k = 4 , 6
            l = (k-1)*nx
            q(l+j) = q(l+j)*thrhal
 20      continue
         av = q(j) + q(nx+j) + q(nx+nx+j)
         do 30 k = 1 , 3
            l = (k-1)*nx
            q(j+l) = q(j+l)*thrhal - av*half
 30      continue
 40   continue
      do 60 i = 1 , 6
         l = (i-1)*nx
         do 50 j = 1 , nword
            q(j) = q(j+l)
 50      continue
         call wrt3(q(1),nword,ibl,ifils)
         ibl = ibl + len
 60   continue
      return
      end
      subroutine qmints(xxs,yys,zzs,xys,xzs,yzs)
c--------------------------------------------------------------
c     quadrupole moment integrals (a.o. basis)
c-------------------------------------------------------------
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
      logical iandj,out,norm,double
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
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
      common/small/dij(225),sxx(225),syy(225),szz(225),
     +sxy(225),sxz(225),syz(225),xin(50),yin(50),zin(50),
     +xxin(25),yyin(25),zzin(25),ijx(225),ijy(225),ijz(225)
      dimension xxs(*),yys(*),zzs(*),xys(*),xzs(*),yzs(*)
      data zero,one/0.0d0,1.0d0/
      origx = gx
      origy = gy
      origz = gz
      tol = 2.30258d0*itol
      out = odebug(19)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 110 ii = 1 , nshell
         i = katom(ii)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , ii
            j = katom(jj)
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            iandj = ii.eq.jj
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,
     +           5,1,1)
            do 20 i = 1 , ij
               sxx(i) = zero
               syy(i) = zero
               sxy(i) = zero
               sxz(i) = zero
               syz(i) = zero
               szz(i) = zero
 20         continue
c     ----- i primitive
            jgmax = j2
            do 70 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c
c     ----- j primtive
c
               if (iandj) jgmax = ig
               do 60 jg = j1 , jgmax
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     double = iandj .and. ig.ne.jg
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                  csj,cpj,cdj,cfj,cgj,
     +                  mini,maxi,minj,maxj,iandj,double,norm)
c     ----- second moment integrals -----
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     in = -5
                     do 40 i = 1 , lit
                        in = in + 5
                        ni = i
                        do 30 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xin(jn+25) = xintx*tinv
                           yin(jn+25) = yinty*tinv
                           zin(jn+25) = zintz*tinv
                           call qmsint()
                           xxin(jn) = xintx*tinv
                           yyin(jn) = yinty*tinv
                           zzin(jn) = zintz*tinv
 30                     continue
 40                  continue
                     do 50 i = 1 , ij
                        nnx = ijx(i)
                        ny = ijy(i)
                        nz = ijz(i)
                        dd = dij(i)
                        sxy(i) = sxy(i) + dd*xin(nnx+25)*yin(ny+25)
     +                           *zin(nz)
                        sxz(i) = sxz(i) + dd*xin(nnx+25)*yin(ny)
     +                           *zin(nz+25)
                        syz(i) = syz(i) + dd*xin(nnx)*yin(ny+25)
     +                           *zin(nz+25)
                        sxx(i) = sxx(i) + dd*xxin(nnx)*yin(ny)*zin(nz)
                        syy(i) = syy(i) + dd*xin(nnx)*yyin(ny)*zin(nz)
                        szz(i) = szz(i) + dd*xin(nnx)*yin(ny)*zzin(nz)
 50                  continue
                  end if
 60            continue
 70         continue
c     ----- set up second moment matrices -----
            max = maxj
            nn = 0
            do 90 i = mini , maxi
               li = loci + i
               in = (li*(li-1))/2
               if (iandj) max = i
               do 80 j = minj , max
                  lj = locj + j
                  jn = lj + in
                  nn = nn + 1
                  xxs(jn) = sxx(nn)
                  yys(jn) = syy(nn)
                  zzs(jn) = szz(nn)
                  xys(jn) = sxy(nn)
                  xzs(jn) = sxz(nn)
                  yzs(jn) = syz(nn)
 80            continue
 90         continue
 100     continue
 110  continue
      if (out) then
         write (iwr,6010)
         call prtris(xxs,num,iwr)
         write (iwr,6020)
         call prtris(yys,num,iwr)
         write (iwr,6030)
         call prtris(zzs,num,iwr)
         write (iwr,6040)
         call prtris(xys,num,iwr)
         write (iwr,6050)
         call prtris(xzs,num,iwr)
         write (iwr,6060)
         call prtris(yzs,num,iwr)
      end if
      return
 6010 format (/,10x,25('-'),/,10x,'xx-second moment integrals',/,10x,
     +        25('-'))
 6020 format (/,10x,25('-'),/,10x,'yy-second moment integrals',/,10x,
     +        25('-'))
 6030 format (/,10x,25('-'),/,10x,'zz-second moment integrals',/,10x,
     +        25('-'))
 6040 format (//' xy-integrals'//)
 6050 format (//' xz-integrals'//)
 6060 format (//' yz-integrals'//)
      end
      subroutine qmderi(dd)
c----------------------------------------------------------------
c     quadrupole derivative integrals
c----------------------------------------------------------------
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
      logical out,norm
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
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
      common/small/dipd(3,3,maxat),qudd(6,3,maxat)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     4    qxx(25),qyy(25),qzz(25),qxxd(25),qyyd(25),qzzd(25),
     5    ijx(100),ijy(100),ijz(100)
      dimension dd(*)
       character *8 comp,grhf,closed
      dimension comp(3)
      dimension derivs(27*maxat)
      equivalence (derivs(1),dipd(1,1,1))
      data comp/'d/dx','d/dy','d/dz'/
      data grhf/'grhf'/
      data closed/'closed'/
      data ndim/5/
      data two,thalf/2.0d0,1.5d0/
      data zer,one /0.0d0,1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
      call secget(isect(7),7,iblok)
      call rdedx(dd,nx,iblok,ifild)
      if (scftyp.ne.grhf .and. scftyp.ne.closed) then
         call secget(isect(10),10,iblok)
         call rdedx(dd(nx+1),nx,iblok,ifild)
         do 20 i = 1 , nx
            dd(i) = dd(i+nx) + dd(i)
 20      continue
      end if
      tol = 2.30258d0*itol
      out = odebug(19)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
      do 30 n = 1 , nat
         zz = czan(n)
         zx = (c(1,n)-gx)*zz
         zy = (c(2,n)-gy)*zz
         zz = (c(3,n)-gz)*zz
         qudd(4,1,n) = thalf*zy
         qudd(4,2,n) = thalf*zx
         qudd(4,3,n) = zer
         qudd(5,1,n) = thalf*zz
         qudd(5,2,n) = zer
         qudd(5,3,n) = thalf*zx
         qudd(6,1,n) = zer
         qudd(6,2,n) = zz*thalf
         qudd(6,3,n) = zy*thalf
         qudd(1,1,n) = two*zx
         qudd(1,2,n) = -zy
         qudd(1,3,n) = -zz
         qudd(2,1,n) = -zx
         qudd(2,2,n) = two*zy
         qudd(2,3,n) = -zz
         qudd(3,1,n) = -zx
         qudd(3,2,n) = -zy
         qudd(3,3,n) = two*zz
 30   continue
c     ----- ishell
      do 110 ii = 1 , nshell
         iat = katom(ii)
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +           ndim,1,1)
            qxxx = 0.0d0
            qxxy = 0.0d0
            qxxz = 0.0d0
            qyyx = 0.0d0
            qyyy = 0.0d0
            qyyz = 0.0d0
            qzzx = 0.0d0
            qzzy = 0.0d0
            qzzz = 0.0d0
            qxyx = 0.0d0
            qxyy = 0.0d0
            qxyz = 0.0d0
            qxzx = 0.0d0
            qxzy = 0.0d0
            qxzz = 0.0d0
            qyzx = 0.0d0
            qyzy = 0.0d0
            qyzz = 0.0d0
c     ----- i primitive
            do 90 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 80 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                  csj,cpj,cdj,cfj,cgj,
     +                  mini,maxi,minj,maxj,.false.,.false.,
     +                  norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 50 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 40 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
                           call qmsint()
                           qxx(jn) = xintx*tinv
                           qyy(jn) = yinty*tinv
                           qzz(jn) = zintz*tinv
 40                     continue
 50                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
                     call oneld(qxx,qyy,qzz,qxxd,qyyd,qzzd,ai,lit,ljt,1,
     +                          ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 70 i = mini , maxi
                        in = loci + i
                        do 60 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
                           dum = dd(nn)*dij(n)
                           dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
c
c
                           qxxx = qxxx + dum*qxxd(nnx)*yin(ny)*zin(nz)
                           qxxy = qxxy + dum*qxx(nnx)*yd(ny)*zin(nz)
                           qxxz = qxxz + dum*qxx(nnx)*yin(ny)*zd(nz)
c
c
                           qyyx = qyyx + dum*xd(nnx)*qyy(ny)*zin(nz)
                           qyyy = qyyy + dum*xin(nnx)*qyyd(ny)*zin(nz)
                           qyyz = qyyz + dum*xin(nnx)*qyy(ny)*zd(nz)
c
c
                           qzzx = qzzx + dum*xd(nnx)*yin(ny)*qzz(nz)
                           qzzy = qzzy + dum*xin(nnx)*yd(ny)*qzz(nz)
                           qzzz = qzzz + dum*xin(nnx)*yin(ny)*qzzd(nz)
c
c
                           qxyx = qxyx + dum*xdipd(nnx)*ydip(ny)*zin(nz)
                           qxyy = qxyy + dum*xdip(nnx)*ydipd(ny)*zin(nz)
                           qxyz = qxyz + dum*xdip(nnx)*ydip(ny)*zd(nz)
c
c
                           qxzx = qxzx + dum*xdipd(nnx)*yin(ny)*zdip(nz)
                           qxzy = qxzy + dum*xdip(nnx)*yd(ny)*zdip(nz)
                           qxzz = qxzz + dum*xdip(nnx)*yin(ny)*zdipd(nz)
c
c
                           qyzx = qyzx + dum*xd(nnx)*ydip(ny)*zdip(nz)
                           qyzy = qyzy + dum*xin(nnx)*ydipd(ny)*zdip(nz)
                           qyzz = qyzz + dum*xin(nnx)*ydip(ny)*zdipd(nz)
 60                     continue
 70                  continue
                  end if
 80            continue
 90         continue
            qudd(1,1,iat) = qudd(1,1,iat) - 0.5d0*(qxxx+qxxx-qyyx-qzzx)
            qudd(2,2,iat) = qudd(2,2,iat) - 0.5d0*(qyyy+qyyy-qxxy-qzzy)
            qudd(3,3,iat) = qudd(3,3,iat) - 0.5d0*(qzzz+qzzz-qxxz-qyyz)
            qudd(1,2,iat) = qudd(1,2,iat) - 0.5d0*(qxxy+qxxy-qzzy-qyyy)
            qudd(1,3,iat) = qudd(1,3,iat) - 0.5d0*(qxxz+qxxz-qyyz-qzzz)
            qudd(2,1,iat) = qudd(2,1,iat) - 0.5d0*(qyyx+qyyx-qxxx-qzzx)
            qudd(3,1,iat) = qudd(3,1,iat) - 0.5d0*(qzzx+qzzx-qxxx-qyyx)
            qudd(2,3,iat) = qudd(2,3,iat) - 0.5d0*(qyyz+qyyz-qxxz-qzzz)
            qudd(3,2,iat) = qudd(3,2,iat) - 0.5d0*(qzzy+qzzy-qxxy-qyyy)
            qudd(4,1,iat) = qudd(4,1,iat) - 1.5d0*qxyx
            qudd(4,2,iat) = qudd(4,2,iat) - 1.5d0*qxyy
            qudd(4,3,iat) = qudd(4,3,iat) - 1.5d0*qxyz
            qudd(5,1,iat) = qudd(5,1,iat) - 1.5d0*qxzx
            qudd(5,2,iat) = qudd(5,2,iat) - 1.5d0*qxzy
            qudd(5,3,iat) = qudd(5,3,iat) - 1.5d0*qxzz
            qudd(6,1,iat) = qudd(6,1,iat) - 1.5d0*qyzx
            qudd(6,2,iat) = qudd(6,2,iat) - 1.5d0*qyzy
            qudd(6,3,iat) = qudd(6,3,iat) - 1.5d0*qyzz
 100     continue
 110  continue
c
      if (out) then
         write (iwr,6010)
         do 130 n = 1 , nat
            write (iwr,6020)
            do 120 nc = 1 , 3
               write (iwr,6030) zaname(n) , comp(nc) ,
     +                         (qudd(nn,nc,n),nn=1,6)
 120        continue
 130     continue
      end if
c     term from differentiating orbitals
c     nplus1 = nocca + 1
      iposs = iochf(15)
      icomp = 24
      ioff = 1
      i1 = nx + 1
      i2 = i1 + nx
      i3 = i2 + nx
      lennew = iky(ncoorb+1)
      iblll = lensec(lennew)
      do 140 i = 1 , 3
         ioff = ioff + nx
         icomp = icomp + 1
         call vclr(dd(ioff),1,lennew)
         if (itwo(i).eq.0) then
            write (iwr,6040) i
         else
            call secget(isect(icomp),icomp,iblok)
            call rdedx(dd(ioff),lennew,iblok,ifild)
         end if
 140  continue
      do 160 n = 1 , nat
         do 150 nc = 1 , 3
            call rdedx(dd,lennew,iposs,ifockf)
            iposs = iposs + iblll
            dumx = tracep(dd(i1),dd,ncoorb)
            dumy = tracep(dd(i2),dd,ncoorb)
            dumz = tracep(dd(i3),dd,ncoorb)
            qudd(1,nc,n) = qudd(1,nc,n) - dumx
            qudd(2,nc,n) = qudd(2,nc,n) - dumy
            qudd(3,nc,n) = qudd(3,nc,n) - dumz
 150     continue
 160  continue
      icomp = 27
      ioff = 1
      do 170 i = 4 , 6
         icomp = icomp + 1
         ioff = ioff + nx
         call vclr(dd(ioff),1,lennew)
         if (itwo(i).eq.0) then
            write (iwr,6040) i
         else
            call secget(isect(icomp),icomp,iblok)
            call rdedx(dd(ioff),lennew,iblok,ifild)
         end if
 170  continue
      iposs = iochf(15)
      do 190 n = 1 , nat
         do 180 nc = 1 , 3
            call rdedx(dd,lennew,iposs,ifockf)
            iposs = iposs + iblll
            qudd(4,nc,n) = qudd(4,nc,n) - tracep(dd,dd(i1),ncoorb)
            qudd(5,nc,n) = qudd(5,nc,n) - tracep(dd,dd(i2),ncoorb)
            qudd(6,nc,n) = qudd(6,nc,n) - tracep(dd,dd(i3),ncoorb)
 180     continue
 190  continue
      call qmdsym(qudd,dd,dd(nw196(5)+1),nat,nshell)
      if (oprn(40)) then
       if (oprn(25)) write (iwr,6010)
       do 210 n = 1 , nat
          if (oprn(25)) write (iwr,6020)
          do 200 nc = 1 , 3
             if (oprn(25)) write (iwr,6030) zaname(n) , comp(nc) , 
     +      (qudd(nn,nc,n),nn=1,6)
 200      continue
 210   continue
      endif
      call timit(3)
      lenblk = lensec(27*maxat)
      call secput(isect(50),50,lenblk,iblok)
      call revind
      lds(isect(50)) = 27*maxat
      call wrt3(derivs,lds(isect(50)),iblok,ifild)
      call revise
      call clredx
      return
 6010 format (/
     +   20x,'*************************************'/
     +   20x,'*   quadrupole moment derivatives   *'/
     +   20x,'*           (atomic units)          *'/
     +   20x,'*************************************'//
     + 29x,'xx',14x,'yy',14x,'zz',14x,'xy',14x,'xz',14x,'yz'/)
 6020 format (/)
 6030 format (5x,a8,2x,a8,6f16.8)
 6040 format (//10x,'quadrupole integrals missing for component',i6)
      end
      subroutine qmsint
c------------------------------------------------------------
c     quadrupole moment subsidiary integrals
c-------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj,gx,gy,gz
      common/hermit/h(21)
      common/wermit/w(21)
      dimension min(6),max(6)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data zero /0.0d0/
      xintx = zero
      yinty = zero
      zintz = zero
      npts = (ni+nj)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 140 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
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
         go to (60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      go to (130,120,110,100,90,80,70) , nj
 70      px = px*bx
         py = py*by
         pz = pz*bz
 80      px = px*bx
         py = py*by
         pz = pz*bz
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
 130     xintx = xintx + px*(ptx-gx)*(ptx-gx)
         yinty = yinty + py*(pty-gy)*(pty-gy)
         zintz = zintz + pz*(ptz-gz)*(ptz-gz)
 140  continue
      return
      end
      subroutine qmdsym(qudd,iso,ict,natdim,nshels)
c----------------------------------------------------------------
c     symmetrises quadrupole derivatives
c---------------------------------------------------------------
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
      common/scrtch/
     + pptr(3,144),ddtr(6,288),ftr(10,480),gtr(15,720),
     + dp(3,3,maxat),qd(6,3,maxat)
      common/junkp/ptr(3,144),dtr(6,288)
c
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
      dimension qudd(6,3,maxat)
      dimension iso(nshels,*),ict(natdim,*)
      if (nt.eq.1) return
      call rdedx(ptr,nw196(1),ibl196(1),ifild)
      call rdedx(dtr,nw196(2),ibl196(2),ifild)
      nav = lenwrd()
      call readi(iso,nw196(5)*nav,ibl196(5),ifild)
      zero = 0.0d0
c     one = 1.0d0
      do 40 nd = 1 , 6
         do 30 nc = 1 , 3
            do 20 n = 1 , nat
               qd(nd,nc,n) = zero
 20         continue
 30      continue
 40   continue
c
c     ----- set transformation table: atoms versus symmetry operations.
c
      do 70 ii = 1 , nshell
         ic = katom(ii)
         do 60 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 60      continue
 70   continue
c
c
      do 130 n = 1 , nat
         do 120 nd = 1 , 6
            do 110 nc = 1 , 3
               do 100 nop = 1 , nt
                  icnu = ict(n,nop)
                  nnn = 6*(nop-1)
                  nn = 3*(nop-1)
                  do 90 ndd = 1 , 6
                     do 80 ncd = 1 , 3
                        qd(nd,nc,n) = qd(nd,nc,n) + qudd(ndd,ncd,icnu)
     +                                *dtr(ndd,nnn+nd)*ptr(ncd,nn+nc)
 80                  continue
 90               continue
 100           continue
 110        continue
 120     continue
 130  continue
      dum = dble(nt)
      do 160 n = 1 , nat
         do 150 i = 1 , 6
            do 140 j = 1 , 3
               qudd(i,j,n) = qd(i,j,n)/dum
 140        continue
 150     continue
 160  continue
      return
      end
      subroutine dmints(xs,ys,zs)
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
      logical iandj,out,norm,double
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
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,t,x0,y0,z0,
     1 xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
      common/small/dij(225),xin(125),yin(125),zin(125),
     +       sx(225),sy(225),sz(225),ijx(225),ijy(225),ijz(225)
      dimension xs(*),ys(*),zs(*)
      data zero,one/0.0d0,1.0d0/
      origx = gx
      origy = gy
      origz = gz
      tol = 2.30258d0*itol
      out = odebug(19)
      norm = normf.ne.1 .or. normp.ne.1
c     ----- ishell
      do 110 ii = 1 , nshell
         i = katom(ii)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , ii
            j = katom(jj)
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
            iandj = ii.eq.jj
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,
     +           5,1,1)
            do 20 i = 1 , ij
               sx(i) = zero
               sy(i) = zero
               sz(i) = zero
 20         continue
c     ----- i primitive
            jgmax = j2
            do 70 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               if (iandj) jgmax = ig
               do 60 jg = j1 , jgmax
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     double = iandj .and. ig.ne.jg
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                  csj,cpj,cdj,cfj,cgj,
     +                  mini,maxi,minj,maxj,iandj,double,norm)
c     ----- dipole moment integrals -----
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     in = -5
                     do 40 i = 1 , lit
                        in = in + 5
                        ni = i
                        do 30 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xin(jn+25) = xintx*tinv
                           yin(jn+25) = yinty*tinv
                           zin(jn+25) = zintz*tinv
 30                     continue
 40                  continue
                     do 50 i = 1 , ij
                        nnx = ijx(i)
                        ny = ijy(i)
                        nz = ijz(i)
                        sx(i) = sx(i) + dij(i)
     +                          *(xin(nnx+25)*yin(ny)*zin(nz))
                        sy(i) = sy(i) + dij(i)
     +                          *(xin(nnx)*yin(ny+25)*zin(nz))
                        sz(i) = sz(i) + dij(i)
     +                          *(xin(nnx)*yin(ny)*zin(nz+25))
 50                  continue
                  end if
 60            continue
 70         continue
c     ----- set up dipole moment matrices -----
            max = maxj
            nn = 0
            do 90 i = mini , maxi
               li = loci + i
               in = (li*(li-1))/2
               if (iandj) max = i
               do 80 j = minj , max
                  lj = locj + j
                  jn = lj + in
                  nn = nn + 1
                  xs(jn) = sx(nn)
                  ys(jn) = sy(nn)
                  zs(jn) = sz(nn)
 80            continue
 90         continue
 100     continue
 110  continue
      if (out) then
         write (iwr,6010)
         call prtris(xs,num,iwr)
         write (iwr,6020)
         call prtris(ys,num,iwr)
         write (iwr,6030)
         call prtris(zs,num,iwr)
      end if
      return
 6010 format (/,10x,25('-'),/,10x,'x-dipole moment integrals',/,10x,
     +        25('-'))
 6020 format (/,10x,25('-'),/,10x,'y-dipole moment integrals',/,10x,
     +        25('-'))
 6030 format (/,10x,25('-'),/,10x,'z-dipole moment integrals',/,10x,
     +        25('-'))
      end
      subroutine ver_integd(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integd.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
