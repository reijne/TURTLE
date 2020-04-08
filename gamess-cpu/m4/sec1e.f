c 
c  $Author: mrdj $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/sec1e.m,v $
c  $State: Exp $
c  
      subroutine dr2frm(der2,dd,natom,nat3)
c---------------------------------------------------------------
c     transfers partial answer into force constant matrix
c---------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension der2(nat3,*),dd(12,12),natom(4)
      do 50 i = 1 , 4
         n = natom(i)
         if (n.ne.0) then
            ii = (i-1)*3
            nn = (n-1)*3
            do 40 j = 1 , 4
               m = natom(j)
               if (m.ne.0) then
                  jj = (j-1)*3
                  mm = (m-1)*3
                  do 30 k = 1 , 3
                     do 20 l = 1 , 3
                        dum = dd(ii+k,jj+l)
                        der2(nn+k,mm+l) = der2(nn+k,mm+l) + dum
 20                  continue
 30               continue
               end if
 40         continue
         end if
 50   continue
      return
      end
      subroutine dr2sym(der2,dd,iso,ict,natom,nat3,nshels)
c-----------------------------------------------------------------
c     symmetrisation of second derivatives
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension der2(nat3,*),dd(nat3,*)
      dimension iso(nshels,*),ict(natom,*)
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/blkin/ptr(3,144)
      dimension t(3,3)
      data zero,one/0.0d0,1.0d0/
      if (nt.eq.1) return
      call rdedx(ptr,nw196(1),ibl196(1),ifild)
      do 30 nnb = 1 , nat3
         do 20 nna = 1 , nat3
            dd(nna,nnb) = zero
 20      continue
 30   continue
c     ----- set transformation table: atoms versus symmetry operations.
      do 60 ii = 1 , nshell
         ic = katom(ii)
         do 50 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 50      continue
 60   continue
c
c     loop over atoms
c
      do 110 natoma = 1 , natom
         do 100 natomb = 1 , natom
            n1 = (natoma-1)*3
            n2 = (natomb-1)*3
c     loop over symmetry operations
            do 90 jop = 1 , nt
               nan = ict(natoma,jop)
               nbn = ict(natomb,jop)
c     nuclei nan and nbn are equivalent to natoma and natomb
c
               n3 = (nan-1)*3
               n4 = (nbn-1)*3
c     transformation the 3*3 block corresponding to atoms nan,nbn
c     and add it to 3*3 block (natoma,natomb)
               np = (jop-1)*3
               do 70 k = 1 , 3
                  u1 = der2(n3+k,n4+1)
                  u2 = der2(n3+k,n4+2)
                  u3 = der2(n3+k,n4+3)
                  t(k,1) = u1*ptr(1,np+1) + u2*ptr(2,np+1)
     +                     + u3*ptr(3,np+1)
                  t(k,2) = u1*ptr(1,np+2) + u2*ptr(2,np+2)
     +                     + u3*ptr(3,np+2)
                  t(k,3) = u1*ptr(1,np+3) + u2*ptr(2,np+3)
     +                     + u3*ptr(3,np+3)
 70            continue
               do 80 k = 1 , 3
                  u1 = t(1,k)
                  u2 = t(2,k)
                  u3 = t(3,k)
                  dd(n1+1,n2+k) = u1*ptr(1,np+1) + u2*ptr(2,np+1)
     +                            + u3*ptr(3,np+1) + dd(n1+1,n2+k)
                  dd(n1+2,n2+k) = u1*ptr(1,np+2) + u2*ptr(2,np+2)
     +                            + u3*ptr(3,np+2) + dd(n1+2,n2+k)
                  dd(n1+3,n2+k) = u1*ptr(1,np+3) + u2*ptr(2,np+3)
     +                            + u3*ptr(3,np+3) + dd(n1+3,n2+k)
 80            continue
 90         continue
 100     continue
 110  continue
      dum = one/dble(nt)
      do 130 nnb = 1 , nat3
         do 120 nna = 1 , nat3
            der2(nna,nnb) = dd(nna,nnb)*dum
 120     continue
 130  continue
      return
      end
      subroutine dr2nc0(dd,odebug,iw,iblk,ifil)
c------------------------------------------------------------------
c     second derivatives of nuclear energy
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      logical odebug(*)
      dimension dd(*)
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
      real*8 clat, zlat
      integer nlat, maxlat
      parameter (maxlat=216)
      common /lattic/ clat(3,maxlat),zlat(maxlat),nlat
c
c
      ndim = max(nat,nlat)
      i10 = 1
      i20 = nat*ndim + i10
      i30 = nat*ndim + i20
c     nat3 = nat*3
c     nlen = nat3*nat3
c     last = nlen + nlen + i30
      call dr2nc1(dd(i10),dd(i20),dd(i30),nat,odebug,iw,iblk,ifil)
      return
      end
      subroutine dr2nc1(r2,r5,der2,natom,odebug,iw,iblk,ifil)
c------------------------------------------------------------------
c     second derivatives of nuclear energy
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      logical odebug
      dimension odebug(40)
      dimension r2(natom,*),r5(natom,*),der2(*)
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
      real*8 clat, zlat
      integer nlat, maxlat
      parameter (maxlat=216)
      common /lattic/ clat(3,maxlat),zlat(maxlat),nlat
c
c
c core-core repulsion parameters (see integb_lib for input)
c
      real*8 eta, d
      integer ichgat, indx, jndx
      logical indi
      integer mxpairs
      parameter(mxpairs=10)
      common/repp/eta(mxpairs),d(mxpairs),ichgat(750),npairs,
     +            indx(mxpairs),jndx(mxpairs),indi
      common/small/icol(maxorb)
c
      logical out
      data zero,one,three/0.0d0,1.0d0,3.0d0/
c
      if (indi) call 
     +  caserr('core-core repulsion second derivatives not implemented')
      out = odebug(6) .or. odebug(7)
      nat3 = natom*3
      nlen = nat3*nat3
      call vclr(der2,1,nlen)
      do 20 i = 1 , nat3
         icol(i) = (i-1)*nat3
 20   continue
c
      do 50 k = 2 , natom
         k1 = k - 1
         do 40 l = 1 , k1
            rkl = zero
            do 30 i = 1 , 3
               rkl = rkl + (c(i,k)-c(i,l))**2
 30         continue
            r2(k,l) = rkl
            r2(l,k) = rkl
            r5(k,l) = czan(k)*czan(l)/(rkl**2.5d0)
            r5(l,k) = r5(k,l)
 40      continue
 50   continue
c
      do 90 i = 1 , natom
         i0 = (i-1)*3
         do 80 k = 1 , 3
            do 70 l = 1 , 3
               delt = zero
               if (k.eq.l) delt = one
               do 60 j = 1 , natom
                  if (i.ne.j) then
                     rkl = (c(k,i)-c(k,j))*(c(l,i)-c(l,j))
                     der2(i0+k+icol(i0+l)) = der2(i0+k+icol(i0+l))
     +                  + r5(i,j)*(three*rkl-delt*r2(i,j))
                     j0 = (j-1)*3
                     der2(i0+k+icol(j0+l)) = r5(i,j)
     +                  *(delt*r2(i,j)-three*rkl)
                  end if
 60            continue
 70         continue
 80      continue
 90   continue
      if (nlat.ne.0) then
c
         do 120 k = 1 , natom
            do 110 l = 1 , nlat
               rkl = zero
               do 100 i = 1 , 3
                  rkl = rkl + (c(i,k)-clat(i,l))**2
 100           continue
               r2(k,l) = rkl
               r5(k,l) = czan(k)*zlat(l)/(rkl**2.5d0)
 110        continue
 120     continue
         do 160 i = 1 , natom
            i0 = (i-1)*3
            do 150 k = 1 , 3
               do 140 l = 1 , 3
                  delt = zero
                  if (k.eq.l) delt = one
                  do 130 j = 1 , nlat
                     rkl = (c(k,i)-clat(k,j))*(c(l,i)-clat(l,j))
                     der2(i0+k+icol(i0+l)) = der2(i0+k+icol(i0+l))
     +                  + r5(i,j)*(three*rkl-delt*r2(i,j))
 130              continue
 140           continue
 150        continue
 160     continue
      end if
c
      call vdwaals_hessian(der2,3*natom)
c
c      this now accumulates on section 46
c      it contains the mp2 second derivative
c
      call rdedx(der2(nlen+1),nlen,iblk,ifil)
      do 170 ms = 1 , nlen
         der2(nlen+ms) = der2(nlen+ms) + der2(ms)
 170  continue
      call wrt3(der2(nlen+1),nlen,iblk,ifil)
      if (.not.out) return
      write (iw,6010)
      call prnder(der2,nat3,iw)
      return
 6010 format (//10x,'nuclear contribution to second derivatives'//)
      end
      subroutine dr2ovl(dd,iso,nshels,isec46)
c------------------------------------------------------------------
c     second derivatives of overlap
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension dd(*),iso(*)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      common/junk/spjnk(3,maxat),
     1 xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/blkin/cx,cy,cz,
     2xin(36),yin(36),zin(36),xad(36),yad(36),zad(36),
     3xadd(36),yadd(36),zadd(36),
     4dij(100),ijx(100),ijy(100),ijz(100)
c
      dimension bder(3,3),mi(48)
      logical allout,out,norm
      data zero,one /0.0d0,1.0d0/
      data ndim/6/
c
c     second derivatives of overlap matrix
c
      nat3 = nat*3
      nlen = nat3*nat3
      ioff = lenrel(nw196(5))
      ioffs = nw196(5) + lenint(nat*nt)
      i10 = ioffs + 1
      lab0 = nlen + i10
      lab1 = lab0 + nx
      lab2 = lab1 + num*num
c     lasta = lab2 + num*num
c     lastb = lab0 + nlen
c     last = max(lasta,lastb)
      call rdedx(dd(1),nw196(5),ibl196(5),ifild)
      call dr2lag(dd(lab0),dd(lab1),dd(lab2),num)
      call vclr(dd(i10),1,nlen)
      tol = 2.30258d0*itol
      allout = odebug(7)
      out = allout .or. odebug(6)
      norm = normf.ne.1 .or. normp.ne.1
c
c     ----- ishell
c
      do 300 ii = 1 , nshell
         iat = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii+iliso(it))
            if (id.gt.ii) go to 300
            mi(it) = id
 30      continue
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit2 = lit + 2
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
         do 290 jj = 1 , ii
            jat = katom(jj)
            if (iat.ne.jat) then
               n2 = 0
               do 50 it = 1 , nt
                  jd = iso(jj+iliso(it))
                  if (jd.gt.ii) go to 290
                  id = mi(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 290
                  if (id.eq.ii .and. jd.eq.jj) n2 = n2 + 1
 50            continue
               q2 = dble(nt)/dble(n2)
               xj = c(1,jat)
               yj = c(2,jat)
               zj = c(3,jat)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions
c
               call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +                    ndim,1,1)
               xxdd = zero
               yydd = zero
               zzdd = zero
               xydd = zero
               xzdd = zero
               yzdd = zero
c
c     ----- i primitive
c
               do 260 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*xi
                  ayi = ai*yi
                  azi = ai*zi
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
c
c     ----- j primitive
c
                  do 250 jg = j1 , j2
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
                        ax = (axi+aj*xj)*aainv
                        ay = (ayi+aj*yj)*aainv
                        az = (azi+aj*zj)*aainv
c
c     ----- density factor
c
                        nn = 0
                        do 210 i = mini , maxi
                           in = iky(loci+i) + lab0 - 1
                           go to (60,70,120,120,80,120,120,100,120,120,
     +                            90,120,120,110,120,120,120,120,120,
     +                            100) , i
 60                        dum1 = csi
                           go to 120
 70                        dum1 = cpi
                           go to 120
 80                        dum1 = cdi
                           go to 120
 90                        dum1 = cfi
                           go to 120
 100                       if (norm) dum1 = dum1*root3
                           go to 120
 110                       if (norm) dum1 = dum1*root5
 120                       do 200 j = minj , maxj
                              go to (130,140,190,190,150,190,190,170,
     +                               190,190,160,190,190,180,190,190,
     +                               190,190,190,170) , j
 130                          dum2 = dum1*csj
                              go to 190
 140                          dum2 = dum1*cpj
                              go to 190
 150                          dum2 = dum1*cdj
                              go to 190
 160                          dum2 = dum1*cfj
                              go to 190
 170                          if (norm) dum2 = dum2*root3
                              go to 190
 180                          if (norm) dum2 = dum2*root5
 190                          nn = nn + 1
                              dum = dd(in+locj+j)
                              dij(nn) = dum2*q2*(dum+dum)
 200                       continue
 210                    continue
c
c     ----- overlap
c
                        t = dsqrt(aa)
                        tinv = one/t
                        x0 = ax
                        y0 = ay
                        z0 = az
                        in = -ndim
                        do 230 i = 1 , lit2
                           in = in + ndim
                           ni = i
                           do 220 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call stvin2()
                              xin(jn) = xint*tinv
                              yin(jn) = yint*tinv
                              zin(jn) = zint*tinv
 220                       continue
 230                    continue
                        call oneld(xin,yin,zin,xad,yad,zad,ai,lit+1,ljt,
     +                             1,ndim)
                        call oneld(xad,yad,zad,xadd,yadd,zadd,ai,lit,
     +                             ljt,1,ndim)
                        do 240 i = 1 , ij
                           nnx = ijx(i)
                           nny = ijy(i)
                           nnz = ijz(i)
                           d1 = dij(i)
                           xxdd = xxdd + d1*xadd(nnx)*yin(nny)*zin(nnz)
                           yydd = yydd + d1*xin(nnx)*yadd(nny)*zin(nnz)
                           zzdd = zzdd + d1*xin(nnx)*yin(nny)*zadd(nnz)
                           xydd = xydd + d1*xad(nnx)*yad(nny)*zin(nnz)
                           xzdd = xzdd + d1*xad(nnx)*yin(nny)*zad(nnz)
                           yzdd = yzdd + d1*xin(nnx)*yad(nny)*zad(nnz)
 240                    continue
                     end if
c
c     ----- end of primitive loops -----
c
 250              continue
 260           continue
c
c     ----- form integrals over derivatives -----
c
c
c      add contributions to hessian matrix
c
               i0 = (iat-1)*3
               j0 = (jat-1)*3
               bder(1,1) = xxdd
               bder(2,2) = yydd
               bder(3,3) = zzdd
               bder(1,2) = xydd
               bder(1,3) = xzdd
               bder(2,3) = yzdd
               bder(2,1) = bder(1,2)
               bder(3,1) = bder(1,3)
               bder(3,2) = bder(2,3)
               if (allout) write (iwr,6010) ii , jj , bder
               do 280 k = 1 , 3
                  do 270 l = 1 , 3
                     ici = (i0+l-1)*nat3 + i10 - 1
                     icj = (j0+l-1)*nat3 + i10 - 1
                     dd(i0+k+ici) = dd(i0+k+ici) + bder(k,l)
                     dd(j0+k+icj) = dd(j0+k+icj) + bder(k,l)
                     dd(i0+k+icj) = dd(i0+k+icj) - bder(k,l)
                     dd(j0+k+ici) = dd(j0+k+ici) - bder(k,l)
 270              continue
 280           continue
            end if
 290     continue
 300  continue
      call rdedx(dd(lab0),nlen,isec46,ifild)
      do 310 i = 1 , nlen
         dd(lab0-1+i) = dd(lab0-1+i) + dd(i10-1+i)
 310  continue
      call wrt3(dd(lab0),nlen,isec46,ifild)
      if (.not.out) return
      call dr2sym(dd(i10),dd(lab0),iso(1),iso(ioff+1),nat,nat3,
     +            nshels)
c
      write (iwr,6020)
      call prnder(dd(i10),nat3,iwr)
      return
 6010 format (5x,2i5,5x,3f15.7/20x,3f15.7/20x,3f15.7)
 6020 format (//5x,'contribution from second derivatives of',
     +        ' overlap matrix'//)
      end
      subroutine dr2ke(dd,iso,nshels,isec46)
c------------------------------------------------------------------
c     second derivatives of kinetic energy
c------------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension dd(*),iso(*)
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      common/junk/spjnk(3,maxat),
     1 xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/blkin/cx,cy,cz,
     2xin(36),yin(36),zin(36),xad(36),yad(36),zad(36),
     3xadd(36),yadd(36),zadd(36),
     4delx(36),dely(36),delz(36),delxa(36),delya(36),delza(36),
     5delxaa(36),delyaa(36),delzaa(36),
     6dij(100),ijx(100),ijy(100),ijz(100)
c
      logical allout,out,norm
      dimension bder(3,3),mi(48)
      data two,half/2.0d0,0.5d0/
      data zero,one /0.0d0,1.0d0/
      data ndim/6/
c
c     second derivatives of kinetic energy matrix
c
      nat3 = nat*3
      nlen = nat3*nat3
      ioff = lenrel(nw196(5))
      ioffs = nw196(5) + lenint(nat*nt)
      i10 = ioffs + 1
      lab0 = nlen + i10
      lab1 = lab0 + nx
c     lasta = lab1 + nx
c     lastb = lab0 + nlen
c     last = max(lasta,lastb)
      call rdedx(dd(1),nw196(5),ibl196(5),ifild)
      call onepdm(dd(lab0),dd(lab1))
      call vclr(dd(i10),1,nlen)
      tol = 2.30258d0*itol
      allout = odebug(7)
      out = allout .or. odebug(6)
      norm = normf.ne.1 .or. normp.ne.1
c
c     ----- ishell
c
      do 300 ii = 1 , nshell
         iat = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii+iliso(it))
            if (id.gt.ii) go to 300
            mi(it) = id
 30      continue
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit2 = lit + 2
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
         do 290 jj = 1 , ii
            jat = katom(jj)
            if (iat.ne.jat) then
               n2 = 0
               do 50 it = 1 , nt
                  jd = iso(jj+iliso(it))
                  if (jd.gt.ii) go to 290
                  id = mi(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 290
                  if (id.eq.ii .and. jd.eq.jj) n2 = n2 + 1
 50            continue
               q2 = dble(nt)/dble(n2)
               xj = c(1,jat)
               yj = c(2,jat)
               zj = c(3,jat)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions
c
               call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,.false.,
     +                    ndim,1,1)
               xxdd = zero
               yydd = zero
               zzdd = zero
               xydd = zero
               xzdd = zero
               yzdd = zero
c
c     ----- i primitive
c
               do 260 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*xi
                  ayi = ai*yi
                  azi = ai*zi
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
c
c     ----- j primtive
c
                  do 250 jg = j1 , j2
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
                        ax = (axi+aj*xj)*aainv
                        ay = (ayi+aj*yj)*aainv
                        az = (azi+aj*zj)*aainv
c
c     ----- density factor
c
                        nn = 0
                        do 210 i = mini , maxi
                           in = iky(loci+i) + lab0 - 1
                           go to (60,70,120,120,80,120,120,100,120,120,
     +                            90,120,120,110,120,120,120,120,120,
     +                            100) , i
 60                        dum1 = csi
                           go to 120
 70                        dum1 = cpi
                           go to 120
 80                        dum1 = cdi
                           go to 120
 90                        dum1 = cfi
                           go to 120
 100                       if (norm) dum1 = dum1*root3
                           go to 120
 110                       if (norm) dum1 = dum1*root5
 120                       do 200 j = minj , maxj
                              go to (130,140,190,190,150,190,190,170,
     +                               190,190,160,190,190,180,190,190,
     +                               190,190,190,170) , j
 130                          dum2 = dum1*csj
                              go to 190
 140                          dum2 = dum1*cpj
                              go to 190
 150                          dum2 = dum1*cdj
                              go to 190
 160                          dum2 = dum1*cfj
                              go to 190
 170                          if (norm) dum2 = dum2*root3
                              go to 190
 180                          if (norm) dum2 = dum2*root5
 190                          nn = nn + 1
                              dum = dd(in+locj+j)
                              dij(nn) = dum2*q2*(dum+dum)
 200                       continue
 210                    continue
c
c     ----- overlap
c
                        t = dsqrt(aa)
                        tinv = one/t
                        t1 = -two*aj*aj*tinv
                        t2 = aj*tinv
                        t3 = -half*tinv
                        x0 = ax
                        y0 = ay
                        z0 = az
                        in = -ndim
                        do 230 i = 1 , lit2
                           in = in + ndim
                           ni = i
                           do 220 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call stvin2()
c     overlap integrals
                              xin(jn) = xint*tinv
                              yin(jn) = yint*tinv
                              zin(jn) = zint*tinv
c     elements of del-squared
                              dum = dble(j+j-1)*t2
                              delx(jn) = dum*xint
                              dely(jn) = dum*yint
                              delz(jn) = dum*zint
                              nj = j + 2
                              call stvin2()
                              delx(jn) = delx(jn) + xint*t1
                              dely(jn) = dely(jn) + yint*t1
                              delz(jn) = delz(jn) + zint*t1
                              if (j.gt.2) then
                                 nj = j - 2
                                 call stvin2()
                                 dum = dble((j-1)*(j-2))*t3
                                 delx(jn) = delx(jn) + xint*dum
                                 dely(jn) = dely(jn) + yint*dum
                                 delz(jn) = delz(jn) + zint*dum
                              end if
 220                       continue
 230                    continue
                        call oneld(xin,yin,zin,xad,yad,zad,ai,lit+1,ljt,
     +                             1,ndim)
                        call oneld(xad,yad,zad,xadd,yadd,zadd,ai,lit,
     +                             ljt,1,ndim)
                        call oneld(delx,dely,delz,delxa,delya,delza,ai,
     +                             lit+1,ljt,1,ndim)
                        call oneld(delxa,delya,delza,delxaa,delyaa,
     +                             delzaa,ai,lit,ljt,1,ndim)
                        do 240 i = 1 , ij
                           nnx = ijx(i)
                           nny = ijy(i)
                           nnz = ijz(i)
                           d1 = dij(i)
                           xxdd = xxdd +
     +                            d1*(delxaa(nnx)*yin(nny)*zin(nnz)+
     +                            xadd(nnx)*dely(nny)*zin(nnz)+
     +                            xadd(nnx)*yin(nny)*delz(nnz))
                           yydd = yydd +
     +                            d1*(delx(nnx)*yadd(nny)*zin(nnz)+
     +                            xin(nnx)*delyaa(nny)*zin(nnz)+
     +                            xin(nnx)*yadd(nny)*delz(nnz))
                           zzdd = zzdd +
     +                            d1*(delx(nnx)*yin(nny)*zadd(nnz)+
     +                            xin(nnx)*dely(nny)*zadd(nnz)+
     +                            xin(nnx)*yin(nny)*delzaa(nnz))
                           xydd = xydd +
     +                            d1*(delxa(nnx)*yad(nny)*zin(nnz)+
     +                            xad(nnx)*delya(nny)*zin(nnz)+
     +                            xad(nnx)*yad(nny)*delz(nnz))
                           xzdd = xzdd +
     +                            d1*(delxa(nnx)*yin(nny)*zad(nnz)+
     +                            xad(nnx)*dely(nny)*zad(nnz)+
     +                            xad(nnx)*yin(nny)*delza(nnz))
                           yzdd = yzdd +
     +                            d1*(delx(nnx)*yad(nny)*zad(nnz)+
     +                            xin(nnx)*delya(nny)*zad(nnz)+
     +                            xin(nnx)*yad(nny)*delza(nnz))
 240                    continue
                     end if
c
c     ----- end of primitive loops -----
c
 250              continue
 260           continue
c
c
c
c      add contributions to hessian matrix
c
               i0 = (iat-1)*3
               j0 = (jat-1)*3
               bder(1,1) = xxdd
               bder(2,2) = yydd
               bder(3,3) = zzdd
               bder(1,2) = xydd
               bder(1,3) = xzdd
               bder(2,3) = yzdd
               bder(2,1) = bder(1,2)
               bder(3,1) = bder(1,3)
               bder(3,2) = bder(2,3)
               if (allout) write (iwr,6010) ii , jj , bder
               do 280 k = 1 , 3
                  do 270 l = 1 , 3
                     ici = (i0+l-1)*nat3 + i10 - 1
                     icj = (j0+l-1)*nat3 + i10 - 1
                     dd(i0+k+ici) = dd(i0+k+ici) + bder(k,l)
                     dd(j0+k+icj) = dd(j0+k+icj) + bder(k,l)
                     dd(i0+k+icj) = dd(i0+k+icj) - bder(k,l)
                     dd(j0+k+ici) = dd(j0+k+ici) - bder(k,l)
 270              continue
 280           continue
            end if
 290     continue
 300  continue
      call rdedx(dd(lab0),nlen,isec46,ifild)
      do 310 i = 1 , nlen
         dd(lab0-1+i) = dd(lab0-1+i) + dd(i10-1+i)
 310  continue
      call wrt3(dd(lab0),nlen,isec46,ifild)
      if (.not.out) return
      call dr2sym(dd(i10),dd(lab0),iso,iso(ioff+1),nat,nat3,
     +            nshels)
c
      write (iwr,6020)
      call prnder(dd(i10),nat3,iwr)
      return
 6010 format (5x,2i5,5x,3f15.7/20x,3f15.7/20x,3f15.7)
 6020 format (//5x,'contribution from second derivatives of',
     +        ' kinetic energy matrix'//)
      end
      subroutine dr2pe(dd,iso,nshels,isec46)
c------------------------------------------------------------------
c     second derivatives of potential energy
c-----------------------------------------------------------------
      implicit real*8  (a-h,o-z)
      dimension dd(*),iso(*)
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      common/junk/spjunk(3,maxat),
     1 xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/blkin/cx,cy,cz,
     2xaxa,yaya,zaza,xaya,xaza,yaza,
     3xbxb,ybyb,zbzb,xbyb,xbzb,ybzb,
     4xaxb,yayb,zazb,xayb,yaxb,xazb,zaxb,
     5yazb,zayb,bder(12,12),
     6xin(36),yin(36),zin(36),xad(36),yad(36),zad(36),
     7xadd(36),yadd(36),zadd(36),
     8xbd(36),ybd(36),zbd(36),xbdd(36),ybdd(36),zbdd(36),
     9xabd(36),yabd(36),zabd(36),
     adij(100),ijx(100),ijy(100),ijz(100)
c
      logical allout,out,norm,nota,notb,iandj,double
      dimension xxxx(165),natom(4),mi(48)
      equivalence (xxxx(1),xaxa)
      data ten/10.0d0/
      data pi212/1.1283791670955d0/
      data zero,one /0.0d0,1.0d0/
      data ndim/6/
c
c     second derivatives of potential energy matrix
c
      tol10 = one/(ten**itol)
      nat3 = nat*3
      nlen = nat3*nat3
      ioff = lenrel(nw196(5))
      ioffs = nw196(5) + lenint(nat*nt)
      i10 = ioffs + 1
      lab0 = nlen + i10
      lab1 = lab0 + nx
c     lasta = lab1 + nx
c     lastb = lab0 + nlen
c     last = max(lasta,lastb)
      call rdedx(dd(1),nw196(5),ibl196(5),ifild)
      call onepdm(dd(lab0),dd(lab1))
      call vclr(dd(i10),1,nlen)
      tol = 2.30258d0*itol
      allout = odebug(7)
      out = allout .or. odebug(6)
      norm = normf.ne.1 .or. normp.ne.1
      natom(4) = 0
c
c     ----- ishell
c
      do 420 ii = 1 , nshell
         iat = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii+iliso(it))
            if (id.gt.ii) go to 420
            mi(it) = id
 30      continue
         natom(2) = iat
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit2 = lit + 2
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
         do 410 jj = 1 , ii
            n2 = 0
            do 50 it = 1 , nt
               jd = iso(jj+iliso(it))
               if (jd.gt.ii) go to 410
               id = mi(it)
               if (id.lt.jd) then
                  nd = id
                  id = jd
                  jd = nd
               end if
               if (id.eq.ii .and. jd.gt.jj) go to 410
               if (id.eq.ii .and. jd.eq.jj) n2 = n2 + 1
 50         continue
            q2 = dble(nt)/dble(n2)
            jat = katom(jj)
            natom(1) = jat
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            iandj = ii.eq.jj
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            ljt2 = ljt + 2
            nroots = (lit2+ljt)/2
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions
c
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,iandj,ndim,
     +           1,1)
c
c     loop over nuclei
c
            do 400 ic = 1 , nat
               nota = iat.eq.ic
               notb = jat.eq.ic
               if (.not.(nota .and. notb)) then
                  natom(3) = ic
                  do 60 i = 1 , 165
                     xxxx(i) = zero
 60               continue
                  znuc = -czan(ic)
                  cx = c(1,ic)
                  cy = c(2,ic)
                  cz = c(3,ic)
c
c     ----- i primitive
c
                  jgmax = j2
                  do 290 ig = i1 , i2
                     ai = ex(ig)
                     arri = ai*rr
                     axi = ai*xi
                     ayi = ai*yi
                     azi = ai*zi
                     csi = cs(ig)
                     cpi = cp(ig)
                     cdi = cd(ig)
                     cfi = cf(ig)
c
c     ----- j primitive
c
                     if (iandj) jgmax = ig
                     do 280 jg = j1 , jgmax
                        aj = ex(jg)
                        aa = ai + aj
                        aainv = one/aa
                        dum = aj*arri*aainv
                        if (dum.le.tol) then
                           fac = dexp(-dum)*pi212*aainv
                           csj = cs(jg)*fac
                           cpj = cp(jg)*fac
                           cdj = cd(jg)*fac
                           cfj = cf(jg)*fac
                           ax = (axi+aj*xj)*aainv
                           ay = (ayi+aj*yj)*aainv
                           az = (azi+aj*zj)*aainv
c
c     ----- density factor
c
                           nn = 0
                           double = iandj .and. ig.ne.jg
                           jmax = maxj
                           do 220 i = mini , maxi
                              ni = iky(loci+i) + lab0 - 1
                              go to (70,80,130,130,90,130,130,110,130,
     +                               130,100,130,130,120,130,130,130,
     +                               130,130,110) , i
 70                           dum1 = csi
                              go to 130
 80                           dum1 = cpi
                              go to 130
 90                           dum1 = cdi
                              go to 130
 100                          dum1 = cfi
                              go to 130
 110                          if (norm) dum1 = dum1*root3
                              go to 130
 120                          if (norm) dum1 = dum1*root5
 130                          if (iandj) jmax = i
                              do 210 j = minj , jmax
                                 nj = locj + j
                                 go to (140,150,200,200,160,200,200,180,
     +                                  200,200,170,200,200,190,200,200,
     +                                  200,200,200,180) , j
 140                             dum2 = dum1*csj
                                 if (double) then
                                    if (i.gt.1) then
                                       dum2 = dum2 + csi*cpj
                                    else
                                       dum2 = dum2 + dum2
                                    end if
                                 end if
                                 go to 200
 150                             dum2 = dum1*cpj
                                 if (double) dum2 = dum2 + dum2
                                 go to 200
 160                             dum2 = dum1*cdj
                                 if (double) dum2 = dum2 + dum2
                                 go to 200
 170                             dum2 = dum1*cfj
                                 if (double) dum2 = dum2 + dum2
                                 go to 200
 180                             if (norm) dum2 = dum2*root3
                                 go to 200
 190                             if (norm) dum2 = dum2*root5
 200                             nn = nn + 1
                                 dum = dd(ni+nj)
                                 if ((loci+i).ne.(locj+j)) dum = dum +
     +                               dum
                                 dij(nn) = dum2*q2*dum
 210                          continue
 220                       continue
                           aax = aa*ax
                           aay = aa*ay
                           aaz = aa*az
                           pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                           if (nroots.le.3) call rt123
                           if (nroots.eq.4) call roots4
                           if (nroots.eq.5) call roots5
                           if (nroots.gt.5) call rootss
c
c     loop over roots
c
                           do 270 k = 1 , nroots
                              uu = aa*u(k)
                              ww = w(k)*znuc
                              tt = aa + uu
                              t = dsqrt(tt)
                              tinv = one/tt
                              x0 = (aax+uu*cx)*tinv
                              y0 = (aay+uu*cy)*tinv
                              z0 = (aaz+uu*cz)*tinv
                              in = -ndim
                              do 240 i = 1 , lit2
                                 in = in + ndim
                                 ni = i
                                 do 230 j = 1 , ljt2
                                    jn = in + j
                                    nj = j
                                    call stvin2()
                                    xin(jn) = xint
                                    yin(jn) = yint
                                    zin(jn) = zint*ww
 230                             continue
 240                          continue
c
c      differentiate subsiduary functions
c
                              if (.not.(nota)) then
                                 call oneld(xin,yin,zin,xad,yad,zad,ai,
     +                              lit+1,ljt+1,1,ndim)
                                 call oneld(xad,yad,zad,xadd,yadd,zadd,
     +                              ai,lit,ljt,1,ndim)
                                 if (notb) go to 250
                                 call oneld(xad,yad,zad,xabd,yabd,zabd,
     +                              aj,lit,ljt,2,ndim)
                              end if
                              call oneld(xin,yin,zin,xbd,ybd,zbd,aj,
     +                           lit+1,ljt+1,2,ndim)
                              call oneld(xbd,ybd,zbd,xbdd,ybdd,zbdd,aj,
     +                           lit,ljt,2,ndim)
c
c     assemble
c
 250                          do 260 i = 1 , ij
                                 nnx = ijx(i)
                                 nny = ijy(i)
                                 nnz = ijz(i)
                                 d1 = dij(i)
                                 if (dabs(d1).lt.tol10) go to 260
                                 if (.not.(nota)) then
                                    xaxa = xaxa + d1*xadd(nnx)*yin(nny)
     +                                 *zin(nnz)
                                    yaya = yaya + d1*xin(nnx)*yadd(nny)
     +                                 *zin(nnz)
                                    zaza = zaza + d1*xin(nnx)*yin(nny)
     +                                 *zadd(nnz)
                                    xaya = xaya + d1*xad(nnx)*yad(nny)
     +                                 *zin(nnz)
                                    xaza = xaza + d1*xad(nnx)*yin(nny)
     +                                 *zad(nnz)
                                    yaza = yaza + d1*xin(nnx)*yad(nny)
     +                                 *zad(nnz)
                                    if (notb) go to 260
                                    xaxb = xaxb + d1*xabd(nnx)*yin(nny)
     +                                 *zin(nnz)
                                    yayb = yayb + d1*xin(nnx)*yabd(nny)
     +                                 *zin(nnz)
                                    zazb = zazb + d1*xin(nnx)*yin(nny)
     +                                 *zabd(nnz)
                                    xayb = xayb + d1*xad(nnx)*ybd(nny)
     +                                 *zin(nnz)
                                    xazb = xazb + d1*xad(nnx)*yin(nny)
     +                                 *zbd(nnz)
                                    yazb = yazb + d1*xin(nnx)*yad(nny)
     +                                 *zbd(nnz)
                                    yaxb = yaxb + d1*xbd(nnx)*yad(nny)
     +                                 *zin(nnz)
                                    zaxb = zaxb + d1*xbd(nnx)*yin(nny)
     +                                 *zad(nnz)
                                    zayb = zayb + d1*xin(nnx)*ybd(nny)
     +                                 *zad(nnz)
                                 end if
                                 xbxb = xbxb + d1*xbdd(nnx)*yin(nny)
     +                                  *zin(nnz)
                                 ybyb = ybyb + d1*xin(nnx)*ybdd(nny)
     +                                  *zin(nnz)
                                 zbzb = zbzb + d1*xin(nnx)*yin(nny)
     +                                  *zbdd(nnz)
                                 xbyb = xbyb + d1*xbd(nnx)*ybd(nny)
     +                                  *zin(nnz)
                                 xbzb = xbzb + d1*xbd(nnx)*yin(nny)
     +                                  *zbd(nnz)
                                 ybzb = ybzb + d1*xin(nnx)*ybd(nny)
     +                                  *zbd(nnz)
 260                          continue
 270                       continue
                        end if
c     end of loop over roots
c
 280                 continue
 290              continue
c
c     end of loops over primitives
c     now have second derivatives <a''| 1/c  | b>,
c     <a | 1/c  |b''>  and <a' | 1/c | b' >
c     derivatives w.r.t. nuclear coordinate c obtained
c     using translational invariance
c
                  bder(1,1) = xbxb
                  bder(2,2) = ybyb
                  bder(3,3) = zbzb
                  bder(1,2) = xbyb
                  bder(1,3) = xbzb
                  bder(2,3) = ybzb
                  bder(2,1) = bder(1,2)
                  bder(3,1) = bder(1,3)
                  bder(3,2) = bder(2,3)
                  bder(4,4) = xaxa
                  bder(5,5) = yaya
                  bder(6,6) = zaza
                  bder(4,5) = xaya
                  bder(4,6) = xaza
                  bder(5,6) = yaza
                  bder(5,4) = bder(4,5)
                  bder(6,4) = bder(4,6)
                  bder(6,5) = bder(5,6)
                  bder(4,1) = xaxb
                  bder(5,2) = yayb
                  bder(6,3) = zazb
                  bder(5,1) = yaxb
                  bder(6,1) = zaxb
                  bder(6,2) = zayb
                  bder(4,2) = xayb
                  bder(4,3) = xazb
                  bder(5,3) = yazb
                  do 310 i = 4 , 6
                     do 300 j = 1 , 3
                        bder(j,i) = bder(i,j)
 300                 continue
 310              continue
                  if (iat.eq.jat) then
                     do 330 i = 4 , 6
                        do 320 j = 1 , 6
                           bder(i-3,j) = bder(i-3,j) + bder(i,j)
                           bder(i,j) = zero
 320                    continue
 330                 continue
                     do 350 i = 1 , 3
                        do 340 j = 4 , 6
                           bder(i,j-3) = bder(i,j-3) + bder(i,j)
                           bder(i,j) = zero
 340                    continue
 350                 continue
                  end if
                  if (allout) write (iwr,6010) ii , jj ,
     +                               ((bder(i,j),j=1,6),i=1,6)
                  do 370 k = 7 , 9
                     do 360 l = 1 , 6
                        bder(k,l) = -(bder(k-3,l)+bder(k-6,l))
                        bder(l,k) = bder(k,l)
 360                 continue
 370              continue
                  do 390 k = 7 , 9
                     do 380 l = 7 , k
                        bder(k,l) = -(bder(k-3,l)+bder(k-6,l))
                        bder(l,k) = bder(k,l)
 380                 continue
 390              continue
                  call dr2frm(dd(i10),bder,natom,nat3)
               end if
c     end of loop over nuclei
 400        continue
c
 410     continue
 420  continue
c     end of loops over shells
      call rdedx(dd(lab0),nlen,isec46,ifild)
      do 430 i = 1 , nlen
         dd(lab0-1+i) = dd(lab0-1+i) + dd(i10-1+i)
 430  continue
      call wrt3(dd(lab0),nlen,isec46,ifild)
      if (.not.out) return
      call dr2sym(dd(i10),dd(lab0),iso,iso(ioff+1),nat,nat3,
     +            nshels)
c
      write (iwr,6020)
      call prnder(dd(i10),nat3,iwr)
      call rdedx(dd(i10),nlen,isec46,ifild)
      call dr2sym(dd(i10),dd(lab0),iso,iso(ioff+1),nat,nat3,
     +            nshels)
      write (iwr,6030)
      call prnder(dd(i10),nat3,iwr)
      return
 6010 format (5x,2i5,5x,6f15.7,5(/20x,6f15.7))
 6020 format (//5x,'contribution from second derivatives of',
     +        ' potential energy matrix'//)
 6030 format (//5x,'total before two-electron contribution')
      end
      subroutine ver_sec1e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/sec1e.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
