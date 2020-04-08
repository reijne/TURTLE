      subroutine dppop(xscm)
c------
c      sets memory partitioning for standard mulliken and
c      dipole preserving charge analysis with equal partitioning
c      of overlap distribution contributions between originating
c      centres
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
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      character *8 errmsg(3)
      dimension xscm(*)
      data errmsg /'program', 'stop in', '-dppop-'/
c
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
      lwor = igmem_max_memory()
      i10 = igmem_alloc(lwor)
      i20=i10+nx
      i30=i20+nx+1
      i40=i30+nat*4
      last = i40 + 2*(nat*(nat+1))
      need = last - i10
      if (need .gt. lwor) then
        write(iwr,1001)
 1001   format(/,'Insufficient memory allocated for DPC analysis')
        call hnderr(3,errmsg)
      endif
      call dippop(xscm(i10),xscm(i20),xscm(i30),xscm(i40))
      call gmem_free(i10)
c
      return
      end
      subroutine dippop(a,b,q,d)
c------
c      perform dipole preserving analysis, based on equal partitioning
c      of overlap distributions contributions between originating
c      centres
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      dimension a(nx), b(nx)
      dimension q(nat,4), d(nat*(nat+1)/2,4)
c
c-----  common blocks
c
cxxxINCLUDE(comdrf/iofil)
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
      integer icent, itxyz
      common /nmorb/ icent(mxnum),itxyz(mxnum)
c
cxxxINCLUDE(comdrf/opt)
c
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
cxxxINCLUDE(comdrf/scfopt)
      character *8 title2,scftyp
      common/restrz/title2(10),scftyp
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
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
      common /restar/ nprint,ipadrestar(702)
c
c-----  local variables
c
      dimension sum(4)
      dimension imat(4)
      logical out, mcout
c
c-----  data statements
c
c     data imat /53,54,55,12/
      data imat/4,5,6,1/
      dimension potnuc(10),ostf(6)
      logical ostf
      data one,two /1.0d00,2.0d00/
c
c-----  begin
c
      out = nprint .eq. 3
c
      mcout = ((.not. mcupdt) .or. (imcout .eq. 5))
c
      if (mcout) then
        write(iwr,1001)
 1001   format(//t20,22('-')/t20,' dp analysis'/
     1 t20,22('-')/)
      endif
c
c-----  assign overlap distributions
c
      call lookupd
c
c-----  read density
c
      call rdedx(b,nx,ibl3pa,idaf)
cxxx  call daread(idafh,ioda,b,nx,16)
cxxx  if(scftyp.eq.'uhf') then
      if((scftyp.eq.'uhf') .or. (scftyp.eq.'rohf')
     1 .or. (scftyp .eq. 'gvb').or.(scftyp.eq.'grhf')) then
        call rdedx(a,nx,ibl3pb,idaf)
cxxx    call daread(idafh,ioda,a,nx,20)
        do 10 l=1,nx
          b(l)=b(l)+a(l)
   10 continue
      endif
c
c-----  loop over x,y,z dipole moments and overlap of distributions
c
      do 150 l=1,4
c
c  -----  read x, y, z, s
c
        do i=1,6
            ostf(i)=.false.
        enddo
        ostf(imat(l))=.true.
        call getmat(a,a,a,a,a,a,potnuc,num,ostf,ionsec)
c       call daread(idafh,ioda,a,nx,imat(l))
c
c  -----  set up loop over charge distributions
c
        ij=0
        do 140 i=1,num
          do 130 j=1,i
            if(i.eq.j) then
              t=one
            else
              t=two
            endif
            ij=ij+1
            a(ij)=a(ij)*b(ij)*t
  130     continue
  140   continue
c
c  -----  contract dipole and charges into mulliken dipoles
c         and charges
c
c       call condns(num,nat,a,d(1,l),q(1,l),icent,1.,sum(l))
        call condns(num,nat,a,d(1,l),q(1,l),icent,one,sum(l))
c
  150 continue
c
      if (out) call hatout(q,nat,4,5,'dip_q')
c
c-----  rearrange dipole moments for tdpop
c
      nn=0
      do 160 i=1,nat
        do 170 k=1,3
          nn=nn+1
          a(nn)=q(i,k)
  170   continue
  160 continue
c
c-----  calculate mulliken charges and dipoles and determine dp charges
c
      call tdpop(a,q(1,4),auxx)
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,260) sum(4)
  260   format(//t20,'total number of electrons'  ,f10.3/)
      endif
c
      return
      end
      subroutine lookupd
c------
c      looks up atomic centres of origin of overlap distributions
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
c-----  common blocks
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
      integer icent, itxyz
      common /nmorb/ icent(mxnum),itxyz(mxnum)
c
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
      do 50 ii=1,nshell
        i=katom(ii)
        mini=kmin(ii)
        maxi=kmax(ii)
        loci=kloc(ii)-mini
        do 30 m=mini,maxi
          li=loci+m
          itxyz(li)=m
          icent(li)=i
   30   continue
   50 continue
c
      return
      end
      subroutine condns(num,nat,d,c,q,icent,oc,sum)
c------
c  * * * condenses (inter) orbital distributions 'd' to (inter) atomic
c  * * * distributions 'c' and one centre distributions 'q'
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
      dimension d(*),c(*),q(*),icent(*)
c
c  * * *  dimensions:
c  * * *  d  num*(num+1)/2
c  * * *  c  nat*(nat+1)/2
c  * * *  q  nat
c
      data zero,pt5/0.0d00,0.5d00/
c
      nn=nat*(nat+1)/2
      call clear(c,nn)
      call clear(q,nat)
c
      k=0
      sum=zero
      do 20 i=1,num
        ic=icent(i)
        im=ic*(ic-1)/2
        do 10 j=1,i
          k=k+1
          jc=icent(j)
          jm=im+jc
          del=oc*d(k)
          sum=sum+del
          c(jm)=c(jm)+del
          del=pt5*del
          q(ic)=q(ic)+del
          q(jc)=q(jc)+del
   10   continue
   20 continue
c
      return
      end
      subroutine tdpop(p,q,w)
c------
c  * * *  population analysis on occupied orbitals preserving the
c  * * *  dipole moment
c
c         input parameters
c         nat: number of internal nuclei
c         q(nat) gross mulliken electronic charges
c         p(3,nat) gross mulliken electronic dipoles
c         czan(nat) nuclear charges
c
c         output parameters
c         q(nat) first: mulliken gross atomic charges
c         q(nat) second: effective point charges conserving molecular
c         charge and dipole moment
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
      dimension p(3,nat)
      dimension q(nat), w(*)
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
      integer icent, itxyz
      common /nmorb/ icent(mxnum),itxyz(mxnum)
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
      real*8 auxx
      common /aux/ auxx(3*mxpts)
c
C
      integer ihlp
      common /ihelp/ ihlp(maxorb)
c
c
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer nchrp, ichrp
      common /elpinf/ nchrp(21), ichrp(21)
      integer ilmc, ilmd, ildc
      common /analid/ ilmc,ilmd,ildc
c
c
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
c-----  local variables
c
      character*16 names(1000)
c
c-----  data statements
c
      data zero,small/0.0d00,1.d-4/
c
c  * * *  mulliken analysis
c  * * *  reduce electronic charge with nuclear charge
c
      call clear(w,6)
      ctot=zero
      do  6 i=1,nat
cxxx    names(i)=anam(i)//bnam(i)
        names(i)=anam(i)
c
c  * * *  only intrinsic dipoles are redistributed
c
        do 4 l=1,3
          p(l,i)=q(i)*c(l,i)-p(l,i)
    4   continue
        q(i)=czan(i)-q(i)
        ctot=ctot+q(i)
        do 5 l=1,3
          w(l)=w(l)+q(i)*c(l,i)
          w(l+3)=w(l+3)+p(l,i)
    5   continue
    6 continue
c
c    -----  save mulliken charges and dipoles
c
c  * * *  output mulliken charges and dipoles
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,10)
   10   format(//t20,22('-')/t20,'gross mulliken charges'/
     1 t20,22('-')/)
      endif
c
      call printq(2,nat,q,c,names)
c
      if (abs(ctot).ge.small) then
        do 15, l = 1, 3
          w(l) = w(l)/ctot
   15   continue
        if ((.not. mcupdt) .or. (imcout .eq. 5)) then
          write(iwr,16)
   16     format(//t20,26('-')/t20,'center of mulliken charges'/
     1  t20,26('-')/)
          call hatout(w,1,3,2,' ')
          write(iwr,17)
   17     format(//t20,25('-')/t20,'charge free dipole moment'/
     1  t20,25('-')/)
          call hatout(w(4),1,3,2,' ')
        endif
      endif
c
c-----  write expansion centra on -da31-, record -101-
c
      nchrp(2) = nat
      ichrp(2) = 101
      nchrp(3) = nat
      ichrp(3) = 101
      nchrp(4) = nat
      ichrp(4) = 101
cxxx  call dawrit (idafdrf,iodadrf,c,3*nat,101,navdrf)
c
c * * * punch mulliken charges
c
cxxx  write(ipnch,*) title
cxxx  write(ipnch,*)'  name, charge, x, y, z'
cxxx  write(ipnch,*)' $mul-charges'
cxxx  do 20  i=1,nat
cxxx    write(ipnch,650) anam(i),bnam(i),q(i),(c(k,i),k=1,3)
cxx20 continue
cxxx  write(ipnch,*) ' $end'
c
c-----  write mulliken charges on -da31-, record -111-
c
c     call dawrit(idafdrf,iodadrf,q,nat,111,navdrf)
cxxx  call dawrit(idafdrf,iodadrf,q,nat,ilmc,navdrf)
c
c * * * punch mulliken dipoles
c
cxxx  write(ipnch,*)'  name, dipx,dipy,dipz, x, y, z'
cxxx  write(ipnch,*)' $mul-dipoles'
cxxx  do 30  i=1,nat
cxxx    write(ipnch,650) anam(i),bnam(i),(p(l,i),l=1,3),(c(k,i),k=1,3)
cxx30 continue
cxxx  write(ipnch,*) ' $end'
c
c-----  write mulliken dipoles on -da31-, record -112-
c
c     call dawrit(idafdrf,iodadrf,p,3*nat,112,navdrf)
cxxx  call dawrit(idafdrf,iodadrf,p,3*nat,ilmd,navdrf)
c
c-----  calculate dp charges
c
      call dpopan(nat,nat,c,c,q,p,w,ihlp)
c
c  * * *  printed output of td charges
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,500)
  500   format(
     1//t20,36('-'),/t20,'charges preserving the dipole moment',
     2/t20,36('-'))
      endif
c
      call printq(2,nat,q,c,names)
c
c-----  punch charges on internal points
c
cxxx  if(nprint.eq.-5) rewind ipnch
cxxx  write(ipnch,550)
cx550 format('----  dipole conserving charges  ------'/)
cxxx  write(ipnch,*)'  name, charge, x, y, z'
cxxx  write(ipnch,*)' $dp-charges'
cxxx  do 600 i=1,nat
cxxx    write(ipnch,650) anam(i),bnam(i),q(i),(c(k,i),k=1,3)
cx600 continue
cxxx  write(ipnch,*) ' $end'
  650 format(a8,a2,6f10.6)
c
c-----  write dipole preserving charges on -da31-, record -114-
c
      nchrp(5) = nat
      ichrp(5) = 101
c     call dawrit(idafdrf,iodadrf,q,nat,114,navdrf)
cxxx  call dawrit(idafdrf,iodadrf,q,nat,ildc,navdrf)
c
      return
      end
      subroutine dpopan(ndip,npnts,cp,cq,q,p,w,isign)
c------
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c
c     ndip = number of dipole moments
c     npnts= number of expansion points
c     cp   = coordinates of dipoles
c     cq   = coordinates of expansion centers
c     q    = charges(input&output)
c     p    = dipoles
c     w    = work space
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
      dimension cp(3,ndip),q(ndip),p(3,ndip),w(*)
      dimension cq(3,npnts)
      dimension isign(npnts)
c
      dimension ri(3),ar(3),br(3),cr(3),dr(3),b(3,3),xav(3)
      dimension bvec(3,3),ia(3),ib(3)
c
      data ia/0,1,3/
      data ib/0,3,6/
      data zero,one/0.d00,1.d00/
      data thresh,smallw/1.0d-05,1.0d-10/
c
      weight(x)= exp(-x**2)
c
c-----  begin
c
c-----  loop over the muliken charges and dipoles
c
      do 400 i=1,ndip
c
c  -----  position vector in ar
c
        ar(1)=cp(1,i)
        ar(2)=cp(2,i)
        ar(3)=cp(3,i)
c
c  -----  bl is used for scaling the weight function
c
        bl=10000000.0d0
        do 15 l=1,npnts
          call distab(ar,cq(1,l),br,bll)
          if(bll.lt.bl.and.bll.gt.thresh) bl=bll
   15   continue
c
c  -----  loop over expansion points for finding averages etc.
c
        nsign=0
        sum=zero
        call clear(xav,3)
        call clear(b,9)
c
        do 35 l=1,npnts
c
c    -----  position vector in ri
c
          do 20 k=1,3
            ri(k)=cq(k,l)
   20     continue
c
c    -----  distance to dipole(i) into cr
c
          call distab(ri,ar,cr,al)
c
c    -----  weighting function
c
          wl=weight(al/bl)
c
c         if (wl .gt. small) then
c
            nsign=nsign+1
            isign(nsign)=l
            w(nsign)=wl
            sum=sum+wl
            call accum (ri,xav,wl,3)
c
c    -----  update matrix b
c
            do 30 im=1,3
              do 30 jm=1,3
                b(im,jm)=b(im,jm)+ri(im)*ri(jm)*wl
   30       continue
c
c         endif
c
   35   continue
c
c  -----  average coordinates in xav
c
        fact=one/sum
        do 36 im=1,3
          xav(im)=xav(im)*fact
   36   continue
c
c  -----  complete matrix b
c
        do 38 im=1,3
          do 38 jm=1,3
            b(im,jm)=b(im,jm)*fact-xav(im)*xav(jm)
   38   continue
c
c  -----  diagonalize matrix b
c
        call clear(cr,3)
        call shrink(b,b,3)
cxxx    call diagiv(b,bvec,cr,ia,3,3,3)
        call jacobi(b,ia,3,bvec,ib,3,cr,2,2,
     * 1.0d-08)
cxxx    call jacobi(a,iky,newb,q,ilifq,nrow,e,iop1,iop2,
cxxx * thresh)
c
c  -----  charge-free dipole in dr
c
        do 140 l=1,3
          dr(l)=-p(l,i)
  140   continue
c
c  -----  transform dipole vector
c
        call matvec(bvec,dr,ar,3,.true.)
c
c  -----  form u(dag)(p-q<x>)
c
        do 150 im=1,3
          dr(im)=ar(im)/(cr(im)+thresh*(cr(3)+thresh))
  150   continue
c
c  -----  transform result to br
c
        call matvec(bvec,dr,br,3,.false.)
c
c  -----  loop  over significant points and compute partial charges
c
        do 160 l=1,nsign
          ll=isign(l)
c
c    -----  position of point l relative to centre
c
          do 155 k=1,3
            ri(k)=cq(k,ll)-xav(k)
  155     continue
c
c    -----  charge on atom l
c
          qad=w(ll)*fact*ddot(3,br,1,ri,1)
c         qad=w(ll)*fact*adotb(br,ri,3)
          q(ll)=q(ll)-qad
  160   continue
  400 continue
c
      return
      end
      subroutine dppope(xscm)
c------
c      drives the mulliken and dp charge analysis with
c      charges distributions assigned on the basis of
c      distance criteria (in stead of dividing contributions
c      between source centra as in dppop)
c
c      namelist $expanc
c
c      icnexp: electronic expansion centra option
c
c           0: employ only atomic centra as expansion centra (default)
c           1: add centre of charge as expansion centre to the
c              atomic centra
c           2: add centre of charge and centra to be read from
c              $expc as expansion centra to atomic centra
c           3: add centra to be read from $expc as expansion centra
c              to atomic centra
c           4: use centre of charge as sole expansion centre
c           5: use centre of charge and centra to be read from
c              $expc as expansion centra
c           6: use only centra read from $expc as expansion centra
c           7: use each centre of overlap distribution as expansion cent
c              (this option implies the use of atomic centra, since
c               one-centre s-s overlap distributions are located
c               at the atom centre)
c
c          -n: define extra expansion centra if overlap distributions
c              are about equally distant from defined expansion centra
c
c              note: newly defined expansion centra are added to the lis
c              and available straight away
c
c      iexpas: option for assignment of overlap distributions
c
c          -3: all overlap distributions are assigned to the
c              geometric middle of the originating centres
c          -2: assign according to least-squares fit, but in case
c              of ambiguity to centre of charge
c          -1: assign all two-centre overlap distributions to the centre
c              of charge
c           0: if no unambiguous assignment can be made, assign
c              distribution to centre of charge
c           1: if no unambiguous assignment can be made, let it go
c              astray (default)
c           2: assign according to least squares fit to potential
c              due to overlap charge and dipole
c              only with icnexp < 4 !!!
c              because one-centre distributions are assigned to that cen
c
c      ifitc : option for centra at which the potential of overlap
c              distributions is to be calculated for fitting procedure
c              (with iexpas=2 only)
c
c           0: use the same centra as defined by icnexp(above) (default)
c              only with icnexp >= 0!!!!
c           1: use atomic centra only
c           2: use atomic centra + centre of charge
c           3: use atomic centra + centre of charge + centra read
c              from input ($fitc)
c           4: use atomic centra + centra read from input ($fitc)
c           5: use centra read from input ($fitc) only
c
c
c      iexpcn: dp expansion centra option
c
c           0: employ only atomic centra as expansion centra (default)
c           1: add centre of charge as expansion centre to the
c              atomic centra
c           2: add centre of charge and centra to be read from
c              $expx as expansion centra to atomic centra
c           3: add centra to be read from $expx as expansion centra
c              to atomic centra
c           4: use centre of charge as sole expansion centre
c           5: use centre of charge and centra to be read from
c              $expx as expansion centra
c           6: use only centra read from $expx as expansion centra
c           7: use each centre of overlap distribution as expansion cent
c              (this option implies the use of atomic centra, since
c               one-centre s-s overlap distributions are located
c               at the atom centre)
c           !! allowed only with icnexp = 7 !!
c
c      iunit: unit of length for extra expansion centres found on input
c           0: bohr (default)
c           1: angstrom
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
C
      integer ihlp
      common /ihelp/ ihlp(maxorb)
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
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
c
      integer idafdrf, navdrf, iodadrf
      common /drfdaf/ idafdrf,navdrf,iodadrf(2,255)
c
c
c
      integer nchrp, ichrp
      common /elpinf/ nchrp(21), ichrp(21)
      integer ilmc, ilmd, ildc
      common /analid/ ilmc,ilmd,ildc
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
cahv new call to asgnrf
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
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480)
      common /restar/ nprint
c
      dimension xscm(*)
c
c-----  local variables
c
      dimension xex(3,maxex)
      character*80 text
      character *8 errmsg(3)
      logical out, mcout
      logical oradial
      character*10 duma
      dimension cm(3)
c
      logical defnew, each, asscen
c
      namelist /expanc2/icnexp,iexpas,ifitc,iexpcn,iunit
c
c      data icnexp, iexpas, ifitc, iexpcn, iunit  /0,1,0,0,0/
      data  ifitc, iunit  /0,0/
      data errmsg /'program', 'stop in', '-dppope-'/
c      data expat /.false./
      data oradial /.true./
      data bohr /.529177249d00/
      data zero, one /0.0d00,1.0d00/
c
c...  comdrf/expan
c
      icnexp = 0
      iexpas = 1
      iexpcn = 0
      expat  = .false.
c
c-----  begin
c
      noprt = 1
c
      out =noprt.eq.1.and.nprint.eq.3
c
c
      mcout = ((.not. mcupdt) .or. (imcout .eq. 5))
c
      if (mcout) then
        write(iwr,1000)
 1000   format(//t20,22('-')/t20,'alternative dp analysis'/
     1 t20,22('-')/)
      endif
c
      rewind ir
      read(ir,expanc2,end=1)
c
c-----  check validity of options
c
    1 continue
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,expanc2)
      endif
c
      if ((icnexp .le. 0 .or. icnexp .eq. 3 .or.
     1     icnexp .eq. 6) .and. (iexpas .eq. 0 .or.
     2     iexpas .eq. -1 .or. iexpas .eq. -2)) then
        write(iwr,1001) icnexp, iexpas
 1001   format(/,' input parameters icnexp= ',i2,
     1         ' and iexpas= ',i2,' conflict',/,
     2         ' please change them on input: $expanc ')
        call hnderr(3,errmsg)
      endif
      if (iexpcn .eq. 7 .and. icnexp .ne. 7) then
        write(iwr,1011) iexpcn, icnexp
 1011   format(/,' input parameters iexpcn= ',i2,
     1         ' and icnexp= ',i2,' conflict',/,
     2         ' please change them on input: $expanc ')
        call hnderr(3,errmsg)
      endif
      if (iexpas .eq. 2) then
        if (ifitc .eq. 0 .and. icnexp .lt. 0) then
          write(iwr,1021) icnexp, ifitc
 1021     format(/,' input parameters icnexp= ',i2,
     1         ' and ifitc= ',i2,' conflict',/,
     2         ' please change them on input: $expanc ')
          call hnderr(3,errmsg)
        endif
        if (icnexp .gt. 3) then
          write(iwr,1022) icnexp, iexpas
 1022     format(/,' input parameters icnexp= ',i2,
     1         ' and iexpas= ',i2,' conflict',/,
     2         ' please change them on input: $expanc ')
          call hnderr(3,errmsg)
        endif
      endif
c
      ufact = one
      if (iunit .eq. 1) ufact = one/bohr
c
c
c-----  define electronic expansion centra
c
      icnexa = abs(icnexp)
      defnew = icnexp .lt. 0
      each = icnexa .eq. 7 .or. iexpas .eq. -3
      asscen = iexpas .le. 0
c
c-----  initialize
c
      nexpc=0
      nex = 0
c
c-----  check if atomic centra are required as expansion centra
c
      if (icnexa .le. 3 .or. each) then
        nexpc = nat
        expat = .true.
      endif
c
c-----  check if centre of charge is required as expansion centre
c
      if (icnexa .eq. 1 .or. icnexa .eq. 2 .or.
     1      icnexa .eq. 4 .or. icnexa .eq. 5) then
        nex = nex + 1
        expnam(nex) = ' cent '
        call drfcm(xex(1,nex))
        icch = nexpc + 1
      endif
c
c-----  read expansion centra from input if required
c       (data block beginning with $expc)
c       and store them temporarily in xex
c       currently, a maximum of 50 expansion centra to be read in is all
c
      if (icnexa .eq. 2 .or. icnexa .eq. 3 .or.
     1      icnexa .eq. 5 .or. icnexa .eq. 6) then
        rewind(ir)
   10   read(ir,1100,end=99) text
 1100   format(a80)
        if (text(:6) .ne. ' $expc') goto 10
   20   read(ir,1100) text
        if (text(:5) .eq. ' $end') goto 30
        nex = nex + 1
c
        if (nex .gt. maxex) then
          write(iwr,1200) maxex
 1200     format(/' number of extra expansion centra on $expc ',
     1             'greater than ',i4,/,
     2             ' increase maxex near common blocks /expanc/')
          call hnderr(3,errmsg)
        endif
c
        read(text,1300) expnam(nex),(xex(k,nex),k=1,3)
 1300   format(a10,3f20.10)
        do 25, k = 1, 3
          xex(k,nex) = xex(k,nex)*ufact
   25   continue
        goto 20
   99   write(iwr,*)' list $expc not found on input file, ',
     1               'please check'
        call hnderr(3,errmsg)
   30   continue
      endif
c
      nexpc = nexpc + nex
c
c-----  check assignment method and define points in which potential is
c       to be calculated
c
      nfitc = 0
      nexf = 0
      if (iexpas .eq. 2) then
        if (ifitc .eq. 0) then
          nfitc = nexpc
        endif
        if (ifitc .ge. 1 .and. ifitc .le. 4) then
          nfitc = nat
        endif
        if (ifitc .eq. 2 .or. ifitc .eq. 3) then
          nexf = nexf + 1
        endif
        if (ifitc .ge. 3 .and. ifitc .le. 5) then
c
c    -----  dummy read points where potential is to be read
c
          rewind(ir)
   11     read(ir,1100,end=199) text
          if (text(:6) .ne. ' $fitc') goto 11
   21     read(ir,1100) text
          if (text(:5) .eq. ' $end') goto 31
          nexf = nexf + 1
c
          if (nexf .gt. maxex) then
            write(iwr,1201) maxex
 1201       format(/' number of extra expansion centra on $fitc ',
     1             'greater than ',i4,/,
     2             ' increase maxex near common blocks /expanc/')
            call hnderr(3,errmsg)
          endif
c
          read(text,1300) duma,dum,dum,dum
          goto 21
  199     write(iwr,*)' list $fitc not found on input file, ',
     1               'please check'
          call hnderr(3,errmsg)
   31     continue
        endif
        nfitc = nfitc + nexf
      endif
c
c-----  set memory partitioning
c
c        xscm
c        ovlap at ixo
c        dipx  at ixdx
c        dipy  at ixdy
c        dipz  at ixdz
c        fexp  at ixfx
c        cexp  at ixc
c        iexpc ar ixie
c
c        note: integer array before some real arrays # expansion centra
c              may be enlarged in asgnrf!!!
c
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
      lwor = igmem_max_memory()
      i10 = igmem_alloc(lwor)
      ixo = i10
      ixdx = ixo + nx
      ixdy = ixdx + nx
      ixdz = ixdy + nx
      ixfx = ixdz + nx
      ixel = ixfx + 3*nfitc + 1
cahv
      ixis = ixel + nfitc*225 + 1
cahv
      ixie = ixis + nw196(5)
      ixdel = ixie + 2*nx + 1
      ixc = ixdel + nexpc*225
      last = ixie + 3*nexpc + 1
c
c     need = loc10 + last - i10
c     need = loc10 + last - i10
c     call setc(need)
      call clear(xscm(ixo),last-ixo+1)
c
c-----  fill in expansion centra in xscm
c
      ixex = ixc-1
      if (icnexa .le. 3 .or. each) then
c
c  -----  atomic centra
c
        do 100, i = 1, nat
          do 200, k= 1, 3
            ixex = ixex + 1
            xscm(ixex) = c(k,i)
  200     continue
  100   continue
      endif
c
c-----  extra expansion centra
c
      do 300, i = 1, nex
        do 400, k = 1, 3
          ixex = ixex + 1
          xscm(ixex) = xex(k,i)
  400   continue
  300 continue
c
c-----  define centra at which potential is to be calculated
c
      if (iexpas .eq. 2) then
        if (ifitc .eq. 0) then
c
c    -----  copy expansion centra
c
          do 410, i = 1, 3*nexpc
            xscm(ixfx+i-1) = xscm(ixc+i-1)
  410     continue
        else
          ixf = ixfx - 1
          if (ifitc .ge. 1 .and. ifitc .le. 4) then
c
c      -----  copy nuclear coordinates
c
            do 420, i = 1, nat
              do 430, k = 1, 3
                ixf = ixf + 1
                xscm(ixf) = c(k,i)
  430         continue
  420       continue
          endif
c
          if (ifitc .eq. 2 .or. ifitc .eq. 3) then
c
c      -----  add centre of charge to centra at which potential
c             to be calculated
c
            call drfcm(xscm(ixf+1))
            ixf = ixf + 3
          endif
          if (ifitc .ge. 3 .and. ifitc .le. 5) then
c
c      -----  read extra centra at which the potential is to be evaluate
c
            rewind(ir)
   12       read(ir,1100) text
            if (text(:6) .ne. ' $fitc') goto 12
   22       read(ir,1100) text
            if (text(:5) .eq. ' $end') goto 32
c
            read(text,1300) duma,(xscm(ixf+k),k=1,3)
            do 23, k = 1, 3
              xscm(ixf+k) = xscm(ixf+k)*ufact
   23       continue
            ixf = ixf + 3
            goto 22
   32       continue
c
          endif
        endif
      endif
c
c-----  transform expansion centra with centre of charge
c
      call drfcm(cm)
      ic = ixc-1
      do 450, i = 1, nexpc
        do 460, k = 1, 3
          ic = ic + 1
          xscm(ic) = xscm(ic) - cm(k)
  460   continue
  450 continue
c
c-----  assign expansion centra according to requested criteria
c
cahv      call asgnrf(nfitc,nexpc,icch,iexpas,defnew,each,asscen,
cahv     1     xscm(ixo),xscm(ixdx),xscm(ixdy),xscm(ixdz),
cahv     1     xscm(ixfx),xscm(ixel),xscm(ixdel),
cahv     2     xscm(ixc),xscm(ixie))
c
c----- read in transformation matrices for s,p,d,f basis functions.
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      call rdedx(dtr,nw196(2),ibl196(2),idaf)
      call rdedx(ftr,nw196(3),ibl196(3),idaf)
c
c----- read in symmetry array - iso
c
      call rdedx(xscm(ixis),nw196(5),ibl196(5),idaf)
c
      call asgnrf(nfitc,nexpc,icch,iexpas,defnew,each,asscen,
     1 oradial,nshell,xscm(ixo),xscm(ixdx),xscm(ixdy),xscm(ixdz),
     1     xscm(ixfx),xscm(ixel),xscm(ixdel),
     2     xscm(ixc),xscm(ixie),xscm(ixis))
c
      if (out)  call imatou(xscm(ixie),num,num,3,'iexpc')
c
c-----  electronic expansion centra at ixce
c
      ixce = i10
      ixex = ixce
c
      if (.not. expat) then
c
c  -----  add the nuclei as expansion centra (for mulliken analysis)
c
        nexpc = nexpc + nat
        ixex = ixex - 1
        do 500, i = 1, nat
          do 510, k = 1, 3
            ixex = ixex + 1
            xscm(ixex) = c(k,i)
  510     continue
  500   continue
c
c  -----  reassign the overlap distributions
c
        do 520, i = 1, nx
          xscm(ixie+i-1) = xscm(ixie+i-1) + nat
  520   continue
c
      endif
c
c-----  rewrite expansion centra at xscm(ixc), overwriting overlap etc
c
      do 530, i = 1, 3*nexpc
        xscm(ixex+i-1) = xscm(ixc+i-1)
  530 continue
c
c-----  count overall expansion centra
c
c-----  initialize
c
      nexpx=0
      nex = 0
c
c-----  check if atomic centra are required as expansion centra
c
      if (iexpcn .le. 3) then
        nexpx = nat
      endif
c
c-----  check if centre of charge is required as expansion centre
c
      if (iexpcn .eq. 1 .or. iexpcn .eq. 2 .or.
     1      iexpcn .eq. 4 .or. iexpcn .eq. 5) then
        nex = nex + 1
      endif
c
c-----  count expansion centra from input if required
c       (data block beginning with $expx)
c       and store them temporarily in xex
c       currently, a maximum of 50 expansion centra to be read in is all
c
      if (iexpcn .eq. 2 .or. iexpcn .eq. 3 .or.
     1      iexpcn .eq. 5 .or. iexpcn .eq. 6) then
        rewind(ir)
  110   read(ir,1100,end=98) text
        if (text(:6) .ne. ' $expx') goto 110
  120   read(ir,1100) text
        if (text(:5) .eq. ' $end') goto 130
        nex = nex + 1
c
        if (nex .gt. maxex) then
          write(iwr,1210) maxex
 1210     format(/' number of extra expansion centra on $expx ',
     1             'greater than ',i4,/,
     2             ' increase maxex near common blocks /expanc/')
          call hnderr(3,errmsg)
        endif
c
        read(text,1300) duma,dum,dum,dum
        goto 120
   98   write(iwr,*)' list $expx not found on input file, ',
     1               'please check'
        call hnderr(3,errmsg)
  130   continue
      endif
c
      nexpx = nexpx + nex
c
c-----  overall expansion centra at ixcx
c
      ixcx = ixce + 3*nexpc
c
c-----  further memory partitioning
c
c       ovlap,dipx,dipy,dipz at ixo
c       density              at ixd(respectively)
c       charges and dipoles  at ixq
c       workspace            at ixw
c       dp charges           at ixdp
c       workspace            at ixis
c
c
      ixo = ixcx + 3*nexpx
      ixd = ixo + nx
      ixq = ixd + nx
      ixw = ixq + 4*nexpc
      ixdp = ixw + max(nexpx,6)
      ixie = ixdp + nexpx
      ixis = ixie + nx
      last = ixis + nexpx + 1
c     need = loc10 + last - i10
c     call setc(need)
c
      nnew = 3*nx + 4*nexpc + max(nexpx,6) + 5*nexpx
      call clear(xscm(ixcx),nnew)
c
      call daread(idafh,ioda,xscm(ixie),nx,46)
c
c-----  call mulliken and dp charge calc. driver
c
      call dippope(nexpc,nexpx,xscm(ixce),xscm(ixcx),xscm(ixie),
     1             xscm(ixo),xscm(ixd),xscm(ixq),xscm(ixw),
     2             xscm(ixdp),xscm(ixis))
c     call setc(loadcm)
      call gmem_free(i10)
c
      return
      end
      subroutine dippope(nexpc,nexpx,cexp,xexp,iexp,a,b,q,w,qdp,is)
c------
c      calculates dp charges
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
c-----  dummy variables
c
      dimension cexp(3,nexpc), xexp(3,nexpx), q(nexpc,4)
      dimension iexp(nx), is(nexpx)
      dimension a(nx), b(nx), w(*)
      dimension qdp(nexpx)
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
      common /restar/ nprint
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
      character*8 scftyp
      common /scfopt2/ scftyp
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
c-----  local variables
c
      dimension imat(4)
      dimension sum(4)
      logical out
c
      data one,two /1.0d00,2.0d00/
      data imat /53,54,55,12/
c
c-----  begin
c
      noprt = 1
c
      out =noprt.eq.1.and.nprint.eq.3
c
c-----  get density in array b
c
      call daread(idafh,ioda,b,nx,16)
      if(scftyp.eq.'uhf') then
c
c  -----  add beta density if uhf wave function
c
        call daread(idafh,ioda,a,nx,20)
        do 10 l=1,nx
          b(l)=b(l)+a(l)
   10   continue
      endif
c
c-----  loop over x,y,z dipoles and overlap of charge distributions
c
      do 150 l=1,4
c
c  -----  read x,y,z or s into array a
c
        call daread(idafh,ioda,a,nx,imat(l))
c
c  -----  loop over charge distributions
c
        ij = 0
        do 140 i=1,num
          do 130 j=1,i
            if(i.eq.j) then
              t=one
            else
              t=two
            endif
            ij=ij+1
            a(ij)=a(ij)*b(ij)*t
  130     continue
  140   continue
c
c  -----  condense contributions of charge distributions
c
c       call condnse(nexpc,a,q(1,l),iexp,1.,sum(l))
        call condnse(nexpc,a,q(1,l),iexp,one,sum(l))
  150 continue
      if (out) call hatout(q,nexpc,4,5,'dip_q_exp')
c
c-----  rearrange dipole moments for tdpop
c
      nn=0
      do 160 i=1,nexpc
        do 160 k=1,3
          nn=nn+1
  160     a(nn)=q(i,k)
c
c-----  calculate dp charges
c
      call tdpope(cexp,xexp,a,q(1,4),w,qdp,is)
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,260) sum(4)
  260   format(//t20,'total number of electrons'  ,f10.3/)
      endif
c
      return
      end
      subroutine condnse(nexp,d,q,iexp,oc,sum)
c------
c  * * * condenses (inter) orbital distributions 'd' to
c  * * * one centre distributions 'q' in -nexp- expansion centers
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
c-----  dummy variables
c
      dimension d(nx), q(nexp)
      dimension iexp(nx)
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
c-----  data statements
c
      data zero,pt5/0.0d00,0.5d00/
c
      sum=zero
      do 10, ij = 1, nx
        del=oc*d(ij)
        sum=sum+del
        q(iexp(ij))=q(iexp(ij))+del
   10 continue
c
      return
      end
      subroutine tdpope(cexp,xexp,p,q,w,qdp,isign)
c------
c         population analysis on occupied orbitals preserving the
c         dipole moment
c
c         input parameters
c         nexpc: number of expansion centers
c         p(3,nexpc) gross mulliken electronic dipoles
c         czan(nat) nuclear charges
c         q(nexpc) : mulliken gross atomic charges
c
c         output parameters
c         qdp(nexpx) : effective point charges conserving molecular
c         charge and dipole moment
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
      parameter (nnam=1024)
c
c-----  dummy variables
c
      dimension cexp(3,nexpc), p(3,nexpc)
      dimension xexp(3,nexpx)
      dimension q(nexpc), w(*)
      dimension qdp(nexpx)
      dimension isign(nexpx)
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
      integer nchrp, ichrp
      common /elpinf/ nchrp(21), ichrp(21)
      integer ilmc, ilmd, ildc
      common /analid/ ilmc,ilmd,ildc
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
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
c
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
c-----  local variables
c
      character*16 names(nnam)
      character*80 text
      character *8 errmsg(3)
c
c-----  data statements
c
      data zero,small/0.0d00,1.d-4/
      data errmsg /'program', 'stop in', '-tdpope-'/
c
c-----  begin
c
c-----  check on length of names array
c
      if (nexpc .gt. nnam) then
        write(iwr,1000) nexpc
 1000   format(/,'  too many expansion centra defined for overlap',
     1  ' distributions: ',/,' -----  enlarge nnam in tdpope ',
     2  ' to at least ',i4)
        call hnderr(3,errmsg)
      endif
c
c-----  mulliken analysis
c
      call clear(w,6)
      ctot=zero
c
      do  6 i=1,nexpc
        if(i.le.nat) then
          names(i)=anam(i)//bnam(i)
          zz=czan(i)
        else if (icnexp .lt. 0 .or. icnexp .eq. 7
     1           .or. iexpas .eq. -3) then
          write(names(i),1001) ' cent ', i
 1001     format(a8,i4)
          zz=zero
        else
          names(i)=expnam(i-nat)
          zz=zero
        endif
c
c  -----  only "charge free" dipoles are redistributed
c
        do 4 l=1,3
          p(l,i)=q(i)*cexp(l,i)-p(l,i)
    4   continue
        q(i)=zz-q(i)
        ctot=ctot+q(i)
        do 5 l=1,3
          w(l)=w(l)+q(i)*cexp(l,i)
          w(l+3)=w(l+3)+p(l,i)
    5   continue
    6 continue
c
c-----  save mulliken charges and dipoles
c
      call dawrit (idafdrf,iodadrf,cexp,3*nexpc,102,navdrf)
c
c-----  output mulliken charges and dipoles
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,10)
   10   format(//t20,22('-')/t20,'gross mulliken charges'/
     1 t20,22('-')/)
      endif
c
      call printq(2,nexpc,q,cexp,names)
c
      if(abs(ctot).ge.small) then
        do 15 l=1,3
          w(l)=w(l)/ctot
   15   continue
        if ((.not. mcupdt) .or. (imcout .eq. 5)) then
          write(iwr,16)
   16     format(//t20,26('-')/t20,'center of mulliken charges'/
     1  t20,26('-')/)
          call hatout(w,1,3,2,' ')
          write(iwr,17)
   17     format(//t20,25('-')/t20,'charge free dipole moment'/
     1  t20,25('-')/)
          call hatout(w(4),1,3,2,' ')
        endif
      endif
c
c-----  punch mulliken charges on ipnch
c
      write(ipnch,*) title
      write(ipnch,*) 'expanded'
      write(ipnch,*)'  name, charge, x, y, z'
      write(ipnch,*)' $mul-charges'
      do 20  i=1,nexpc
        write(ipnch,650) names(i),q(i),(cexp(k,i),k=1,3)
   20 continue
      write(ipnch,*) ' $end'
c
c-----  write mulliken charges on -da31-, record -116-
c
      nchrp(7) = nexpc
      ichrp(7) = 102
      nchrp(8) = nexpc
      ichrp(8) = 102
      nchrp(9) = nexpc
      ichrp(9) = 102
c     call dawrit(idafdrf,iodadrf,q,nexpc,116,navdrf)
      call dawrit(idafdrf,iodadrf,q,nexpc,ilmc,navdrf)
c
c-----  punch mulliken dipoles
c
      write(ipnch,*)'  name, dipx,dipy,dipz, x, y, z'
      write(ipnch,*)' $mul-dipoles'
      do 30  i=1,nexpc
        write(ipnch,650) names(i),(p(l,i),l=1,3),(cexp(k,i),k=1,3)
   30 continue
      write(ipnch,*) ' $end'
c
c-----  write mulliken dipoles on -da31-, record -117-
c
c     call dawrit(idafdrf,iodadrf,p,3*nexpc,117,navdrf)
      call dawrit(idafdrf,iodadrf,p,3*nexpc,ilmd,navdrf)
c
c-----  define overall expansion centra
c
      iexpx = 0
c
c-----  check if atomic centra are required as expansion centra
c
      if (iexpcn .le. 3) then
c
c  -----  atomic centra
c
        do 2100, i = 1, nat
          do 2200, k= 1, 3
            xexp(k,i) = c(k,i)
 2200     continue
 2100   continue
c
        iexpx = iexpx + nat
      endif
c
c-----  check if centre of charge is required as expansion centre
c
      if (iexpcn .eq. 1 .or. iexpcn .eq. 2 .or.
     1      iexpcn .eq. 4 .or. iexpcn .eq. 5) then
c
        call drfcm(xexp(1,iexpx+1))
        iexpx = iexpx + 1
      endif
c
c-----  read expansion centra from input if required
c       (data block beginning with $expx)
c
      if (iexpcn .eq. 2 .or. iexpcn .eq. 3 .or.
     1      iexpcn .eq. 5 .or. iexpcn .eq. 6) then
        rewind(ir)
  110   read(ir,1100,end=98) text
 1100   format(a80)
        if (text(:6) .ne. ' $expx') goto 110
c
  120   read(ir,1100) text
        if (text(:5) .eq. ' $end') goto 130
c
        iexpx = iexpx + 1
        read(text,1300) names(iexpx), (xexp(k,iexpx), k= 1,3)
 1300   format(a10,3f20.10)
        do 125, k = 1, 3
          xexp(k,iexpx) = xexp(k,iexpx)
  125   continue
        goto 120
c
   98   write(iwr,*)' list $expx not found on input file, ',
     1               'please check'
        call hnderr(3,errmsg)
  130   continue
      endif
c
      if (iexpcn .eq. 7 .and. icnexp .eq. 7) then
c
c  -----  expansion centra for overlap distributions are
c         also used to expand dp charges
c
        do 2400, i = 1, nexpc
          do 2500, k = 1, 3
            xexp(k,i) = cexp(k,i)
 2500     continue
 2400   continue
      endif
c
      call dpopane(nexpc,nexpx,cexp,xexp,q,p,w,qdp,isign)
c
c-----  printed output of td charges
c
      if ((.not. mcupdt) .or. (imcout .eq. 5)) then
        write(iwr,500)
  500   format(
     1//t20,36('-'),/t20,'charges preserving the dipole moment',
     2/t20,36('-'))
      endif
c
      call printq(2,nexpx,qdp,xexp,names)
c
c-----  punch charges on internal points
c
      write(ipnch,550)
  550 format('----  dipole conserving charges  ------'/)
      write(ipnch,*)'  name, charge, x, y, z'
      write(ipnch,*)' $dp-charges'
      do 600 i=1,nexpx
        write(ipnch,650) names(i),qdp(i),(xexp(k,i),k=1,3)
  600 continue
      write(ipnch,*) ' $end'
  650 format(a10,6f10.6)
c
c-----  write dipole preserving charges on -da31-, record -119-
c
      nchrp(10) = nexpx
      ichrp(10) = 103
      call dawrit (idafdrf,iodadrf,xexp,3*nexpx,103,navdrf)
c     call dawrit(idafdrf,iodadrf,qdp,nexpx,119,navdrf)
      call dawrit(idafdrf,iodadrf,qdp,nexpx,ildc,navdrf)
c
      return
      end
      subroutine dpopane(ndip,npnts,cp,cq,q,p,w,qdp,isign)
c------
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c
c     ndip = number of dipole moments
c     npnts= number of expansion points
c     cp   = coordinates of dipoles
c     cq   = coordinates of expansion centers
c     q    = charges(input)
c     p    = dipoles
c     w    = work space
c     qdp  = charges(output)
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
      dimension cp(3,ndip),q(ndip),p(3,ndip),w(*)
      dimension cq(3,npnts)
      dimension isign(npnts)
      dimension qdp(npnts)
c
      common/iofile/ir,iw,ipnch
c
      dimension ri(3),ar(3),br(3),cr(3),dr(3),b(3,3),xav(3)
      dimension bvec(3,3),ia(3)
      dimension ib(3)
c
      data ia/0,1,3/
      data zero,one/0.d00,1.d00/
      data thresh,smallw/1.0d-05,1.0d-10/
c
      weight(x)= exp(-x**2)
c
c-----  begin
c
c-----  initialize dp charges
c
      call clear(qdp,npnts)
c
c-----  loop over the muliken charges and dipoles
c
      do 400 i=1,ndip
c
c  -----  position vector in ar
c
        ar(1)=cp(1,i)
        ar(2)=cp(2,i)
        ar(3)=cp(3,i)
c
c  -----  bl is used for scaling the weight function
c
        bl=10000000.0d0
        do 15 l=1,npnts
          call distab(ar,cq(1,l),br,bll)
          if(bll.lt.bl.and.bll.gt.thresh) bl=bll
   15   continue
c
c  -----  loop over expansion points for finding averages etc.
c
        nsign=0
        sum=zero
        call clear(xav,3)
        call clear(b,9)
c
        do 35 l=1,npnts
c
c    -----  position vector in ri
c
          do 20 k=1,3
            ri(k)=cq(k,l)
   20     continue
c
c    -----  distance to dipole(i) into cr
c
          call distab(ri,ar,cr,al)
c
c    -----  weighting function
c
          wl=weight(al/bl)
c
c         if (wl .gt. small) then
c
            nsign=nsign+1
            isign(nsign)=l
            w(nsign)=wl
            sum=sum+wl
            call accum (ri,xav,wl,3)
c
c    -----  update matrix b
c
            do 30 im=1,3
              do 30 jm=1,3
                b(im,jm)=b(im,jm)+ri(im)*ri(jm)*wl
   30       continue
c
c         endif
c
   35   continue
c
c       write(iwr,*) 'nsign = ',nsign,(isign(ii),ii=1,nsign)
c       write(iwr,*) 'w = ',(w(ii),ii=1,nsign)
c
c  -----  average coordinates in xav
c
        fact=one/sum
        do 36 im=1,3
          xav(im)=xav(im)*fact
   36   continue
c
c  -----  complete matrix b
c
        do 38 im=1,3
          do 38 jm=1,3
            b(im,jm)=b(im,jm)*fact-xav(im)*xav(jm)
   38   continue
c
c  -----  diagonalize matrix b
c
        call clear(cr,3)
        call shrink(b,b,3)
c       call diagiv(b,bvec,cr,ia,3,3,3)
        call gldiag(3,3,3,b,ib,cr,bvec,ia,2)
c
c  -----  first redistribute the charge and update charge-free
c         dipole resulting from shifting the charge
c
        do 60 l=1,nsign
          ll=isign(l)
c
c    -----  position of point l relative to centre
c
c         do 55 k=1,3
c           ri(k)=cq(k,ll)-xav(k)
c  55     continue
          do 55 k=1,3
            ri(k)=cq(k,ll)-cp(k,i)
   55     continue
c
c    -----  charge on atom l
c
          qad=w(ll)*fact*q(i)
          qdp(ll)=qdp(ll)+qad
c
c    -----  correction to dipole at -i-
c
          do 56, k = 1, 3
            p(k,i) = p(k,i) - qad*ri(k)
   56     continue
c
   60   continue
c
c  -----  charge-free dipole in dr
c
        do 140 l=1,3
          dr(l)=-p(l,i)
  140   continue
c
c  -----  transform dipole vector
c
        call matvec(bvec,dr,ar,3,.true.)
c
c  -----  form u(dag)(p-q<x>)
c
        do 150 im=1,3
          dr(im)=ar(im)/(cr(im)+thresh*(cr(3)+thresh))
  150   continue
c
c  -----  transform result to br
c
        call matvec(bvec,dr,br,3,.false.)
c
c  -----  loop  over significant points and compute partial charges
c
        do 160 l=1,nsign
          ll=isign(l)
c
c    -----  position of point l relative to centre
c
          do 155 k=1,3
            ri(k)=cq(k,ll)-xav(k)
  155     continue
c
c    -----  charge on atom l
c
          qad=w(ll)*fact*ddot(3,br,1,ri,1)
c         qad=w(ll)*fact*adotb(br,ri,3)
          qdp(ll)=qdp(ll)-qad
  160   continue
  400 continue
      return
      end
c      subroutine copyv(a,b,n)
cc------
cc      copies vector a of length n into b
cc------
c      implicit real*8  (a-h,o-z),integer  (i-n)
cINCLUDE(../m4/common/sizes)
cINCLUDE(comdrf/sizesrf)
c      dimension a(n), b(n)
c      do 100, i = 1, n
c        b(i) = a(i)
c  100 continue
c      return
c      end
      subroutine asgnrf(nfit,nexp,icch,iexpas,defnew,each,asscen,
     1           oradial,nshels,ovl,
     1           dipx,dipy,dipz,xexf,elpot,delpot,xexp,iexpc,
     2 iso)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
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
caleko
      common/nottwi/obeen,obeen2,obeen3,obeen4
caleko
      common/junk/s(225),g(225),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,
     *tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      common/blkin/dxyz(4),gg(225),ft(225),fx(225),dij(225),
     * pin(225),qin(225),rin(225),
     * ijx(225),ijy(225),ijz(225)
c mechanics
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
      common/g80nb/natmod,nnbg,mapmod(maxat),ipnonb(maxat+1),
     +             inonb(maxat*6)
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
c ***** omit specified charges from attraction terms ***
c---------
c  omtchg - omit 1e- integrals involving the specified charges
c
c   integrals of the form  < bf(a) | q(b) | bf (a) > are deleted
c
c   where a and b are centres indicated as follows
c
c  ztagcc - list of centre names a
c  ztagcq - for each a in ztagcc, the list of centres b
c
c  omtslf - omit self-energy of point charge array (applied to all
c           atoms with centre label of form bq* )
c      CHECK
c      parameter (mxombq=20,mxomat=50)
c---------

c
c obqf  - special treatment of forces .. skip force contributions
c         on non-bq centres
c         for use in qm/mm optimisations of surrounding
c         assumption is that bq's in this case don't have basis
c         functions or ecps
c
c omitted charge specifications
c

      integer mxombq, mxomat
      parameter (mxombq=50,mxomat=20)
      character*8 ztagcc, ztagcq
      logical omtchg, omtslf, obqf
      common /chgcc/ ztagcc(mxomat),ztagcq(mxomat,mxombq),omtchg,
     &     omtslf,obqf
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c ******
cdrf
c     drf extension
c     ====================================================================
c         note: in hondo, o,x,y and z are real*8
c NCLUDE(comdrf/opt)
caleko
c nexp was not defined here; let's see if this helps
c      common/drfpar2/nxtpts,npol,nexp,natint,nshint,namb,nspec,ngran
caleko
          common/hefcpar/edumm(5,1000),nedumm,iefc
          common/hfldpar/fldxyz(3),ifld
          character*4 keyfld, keyefc, iefc, ifld
c NCLUDE(comdrf/dafil)
         integer idafh,navh
         common/hdafile/idafh,navh,ioda(2,1000)
         common/c_of_m/pcm,qcm,rcm
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      integer isingl, nbits
      common /machin1/ isingl,nbits
c
      integer ixddaf
      integer ndar10, ndar20, ndar31, ndar41, ndar42, ndar43, ndar44
      integer ndar45, ndar46, ndar47, ndar48, ndar49
      common /dafindx/ ixddaf(8192),ndar10,ndar20,
     +                 ndar31,ndar41,ndar42,ndar43,ndar44,ndar45,
     +                 ndar46,ndar47,ndar48,ndar49
c
      integer lrec10, lrec20, lrec31, lrec41, lrec42, lrec43
      integer lrec44, lrec45, lrec46, lrec47, lrec48, lrec49
      common /dafrec/ lrec10,lrec20,lrec31,lrec41,lrec42,lrec43,
     +                lrec44,lrec45,lrec46,lrec47,lrec48,lrec49
c
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
      character*8 scftyp
      common /scfopt2/ scftyp
c
c
      integer nambpt
      common /drfamb/ nambpt(mxat)
c
      real*8 znamb
      common /drfamb2/ znamb(mxat)
c
cahv
c
c-----  dummy variables
c
      real*8 xexp
      dimension xexp(3,*)
      dimension iexpc(nx)
      real*8 ovl 
      dimension ovl(nx)
      dimension dipx(nx), dipy(nx), dipz(nx)
      real*8 xexf
      dimension xexf(3,nfit)
      dimension delpot(nexp,225)
      dimension elpot(nfit,225)
      logical defnew, each, asscen
      logical onec
      logical norm,double
      logical block,blocks,blocki
      logical some,out
      logical repeat, onetwo
      dimension cm(3)
      logical oradial
cdrf  ===================  end drfexts ===============================
cahv      dimension q(*),
      dimension iso(nshels,*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension m0(48)
c
      data tols, told1, told2, told3 /1.d-02, 1.d-04, 1.d-20, 1.d-03/
      data  m51/51/
c mechanics
      data keyfld, keyefc /' fld',' efc'/
      data dzero,pt5,done,two,three,five,seven /0.0d0,0.5d0,1.0d0,
     + 2.0d0,3.0d0,5.0d0,7.0d0/
      data rnine/9.0d0/
      data eleven /11.0d0/
      data pi212 /1.1283791670955d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
cahv
      noprt = 1
      out =noprt.eq.1.and.nprint.eq.3
      num2=(num*(num+1))/2
c
c-----  define expansion centra
c
      call drfcm(cm)
      if (out) call hatout(cm  ,3,1   ,2,'cm  ')
c
c-----  read overlap integrals from da10
c
      call daread(idafh,ioda,ovl, num2,12)
c
c-----  calculate dipole moments relative to -cm-
c
cmw put the centre-of-mass information into the common
c   and call the gamess dipint
c   note that dipint was adapted for cm().
c   a common c_of_m is used to pass cm()
      pcm=cm(1)
      qcm=cm(2)
      rcm=cm(3)
cmw      call dipint(dipx,dipy,dipz,cm,num)
c     call dipint(dipx,dipy,dipz)
cahv  call dipxyz(dipx,dipy,dipz,num)
      call dipxyz(dipx,dipy,dipz)
c
      if (abs(iexpas) .ne. 2) then
c
c----- calculate radial overlap of charge distributions -----
c           ( this is accomplished by setting -lit-
c           and -ljt- equal to -1- regardless of
c           -ktype(ii)- and -ktype(jj)-. we are interested
c           in the -s- component only )
c
        if(out) write(iwr,9999)
 9999 format(//' radial overlap of charge distributions '
     1 /,1x,31(1h-))
c
        onetwo = iexpas .eq. -3
c
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = natmod*3
      l5 = natmod
      outv = oprint(59)
      nav = lenwrd()
      if (nprint.eq.-5) outv = .false.
c
c     ----- set pointers for partitioning of core -----
c
      i10 = 0
      i20 = i10 + l2
      i30 = i20 + l2
cdrf  i301 space for hamilonian in fldint
      i301 = i30 + l2
cdrf  last = i30 + l2
      last = i301 + l2
cIF(parallel)
      if (lpseud.le.1) then
         i40 = i30 + l2
         last = i40 + l2 + l1
      else
         i40 = i30 + l2
         i50 = i40 + l2
         i60 = i50 + l3
         i70 = i60 + 225*num
         i80 = i70 + 225*20
         i90 = i80 + (nat*nt+nav-1)/nav
         i100 = i90 + nw196(4)
         last = i100 + nw196(5)
      end if
c
c crystal field correction
c
c
c ***** allocate memory for charge exclusion
c
      if (omtchg) then
         i110 = last
         last = i110 + nat * nat
      endif
      call izero(l2,iexpc,1)
c
      length = last - i10
 
cahv      i10 = igmem_alloc(length)
      i10 = 1
      i20 = i10 + l2
      i30 = i20 + l2
      if (lpseud.le.1) then
         i40 = i30 + l2
         last = i40 + l2 + l1
      else
         i40 = i30 + l2
         i50 = i40 + l2
         i60 = i50 + l3
         i70 = i60 + 225*num
         i80 = i70 + 225*20
         last = i80 + (nat*nt+nav-1)/nav
      end if
c
c crystal field correction
c
      if(ocryst)then
         ic10 = last
         last = ic10 + l2
      endif
c ...
c mechanics
c ...
      if (oaminc) then
         ia40 = last
         ia50 = ia40 + l4
         ia60 = ia50 + l5
         last = ia60 + l5
      end if
c ***** allocate memory for charge exclusion
c<<<<<< intega.m
c      if (omtchg) then
c       i110 = last
c       last = i110 + nat * nat
c      endif
c======
      if (omtchg) then
         i110 = last
         last = i110 + nat * nat
      endif
c
c     ----- calculate -s- matrix -----
c
c     - s- at x(i10)
c
         cpu = cpulft(1)
cahv
         ncall = 0
cahv
         if (ncall.eq.0) call cpuwal(begin,ebegin)
cahv         if (ncall.eq.0 .and. outv) write (iwr,6030) cpu
         tol = rln10*itol
         out = nprint.eq.3
         if (out) then
           do ii = 1,6
            oprn(30+ii) = .true.
           enddo
          else
           do ii = 1,6
            if(oprn(30+ii)) out = .true.
           enddo
         endif
         onorm = normf.ne.1 .or. normp.ne.1
         ndum = l2 + l2 + l2
c         call vclr(q(i10),1,ndum)
c
c     ----- ishell
c
cjmht psh said this was an event tracer for tcgmsg and so is redundant
cjmht_IF(parallel)
cjmht         call pg_evbgin('par.')
cjmht_ENDIF
         do 440 ii = 1 , nshell
            i = katom(ii)
	    iat = i
	    iatom = i
c
c     ----- eliminate ishell -----
c
            do 450 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 440
               m0(it) = id
450         continue
            icent = i
            pi = c(1,i)
            qi = c(2,i)
            ri = c(3,i)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            lit = ktype(ii)
            mini = kmin(ii)
            maxi = kmax(ii)
            loci = kloc(ii) - mini
c
c     ----- jshell
c
            do 430 jj = 1 , ii
               j = katom(jj)
	       jat = j
	       onec = iat .eq. jat
               n2 = 0
               do 470 it = 1 , nt
                  jd = iso(jj,it)
                  if (jd.gt.ii) go to 430
                  id = m0(it)
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  if (id.eq.ii .and. jd.gt.jj) go to 430
                  if (id.eq.ii.and.jd.eq.jj) then
                      n2 = n2 + 1
                  end if
470            continue
               q2 = dble(nt)/dble(n2)
               jcent = j
               pj = c(1,j)
               qj = c(2,j)
               rj = c(3,j)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               nroots = (lit+ljt-2)/2 + 1
               rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
               oiandj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
               ij = 0
               max = maxj
               do 30 i = mini , maxi
                  nnx = ix(i)
                  nny = iy(i)
                  nnz = iz(i)
                  if (oiandj) max = i
                  do 20 j = minj , max
                     ij = ij + 1
                     ijx(ij) = nnx + jx(j)
                     ijy(ij) = nny + jy(j)
                     ijz(ij) = nnz + jz(j)
                     if (j.le.1) then
                        ft(ij) = three
                     else if (j.le.4) then
                        ft(ij) = five
                     else if (j.le.10) then
                        ft(ij) = seven
                     else if (j.gt.20) then
                        ft(ij) = eleven
                     else
                        ft(ij) = rnine
                     end if
 20               continue
 30            continue
               do 40 i = 1 , ij
                  s(i) = dzero
 40            continue
c
c     ----- i primitive
c
               jgmax = j2
               do 400 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*pi
                  ayi = ai*qi
                  azi = ai*ri
                  csi = cs(ig)
                  if (oradial) then
                    cpi = cs(ig)
                    cdi = cs(ig)
                    cfi = cs(ig)
                    cgi = cs(ig)
                  else
                    cpi = cp(ig)
                    cdi = cd(ig)
                    cfi = cf(ig)
                    cgi = cg(ig)
                  endif
c
c     ----- j primtive
c
                  if (oiandj) jgmax = ig
                  do 390 jg = j1 , jgmax
                     aj = ex(jg)
                     aa = ai + aj
                     aa1 = done/aa
                     dum = aj*arri*aa1
                     if (dum.le.tol) then
                        fac = dexp(-dum)
                        csj = cs(jg)
                        if (oradial) then
                          cpj = cs(jg)
                          cdj = cs(jg)
                          cfj = cs(jg)
                          cgj = cs(jg)
                        else
                          cpj = cp(jg)
                          cdj = cd(jg)
                          cfj = cf(jg)
                          cgj = cg(jg)
                        endif
                        ax = (axi+aj*pj)*aa1
                        ay = (ayi+aj*qj)*aa1
                        az = (azi+aj*rj)*aa1
                        odoub = oiandj .and. ig.ne.jg
c
c     ----- density factor
c
                        max = maxj
                        nn = 0
                        do 220 i = mini , maxi
                           go to (50,60,120,120,
     +                            70,120,120,80,120,120,
     +                            90,120,120,100,120,120,120,120,120,
     +                            110,
     +                            112,120,120,114,120,120,120,120,120,
     +                            116,120,120,118,120,120), i
 50                        dum1 = csi*fac
                           go to 120
 60                        dum1 = cpi*fac
                           go to 120
 70                        dum1 = cdi*fac
                           go to 120
 80                        if (onorm) dum1 = dum1*sqrt3
                           go to 120
 90                        dum1 = cfi*fac
                           go to 120
 100                       if (onorm) dum1 = dum1*sqrt5
                           go to 120
 110                       if (onorm) dum1 = dum1*sqrt3
                           go to 120
 112                       dum1 = cgi*fac
                           go to 120
 114                       if (onorm) dum1 = dum1*sqrt7
                           go to 120
 116                       if (onorm) dum1 = dum1*sqrt5/sqrt3
                           go to 120
 118                       if (onorm) dum1 = dum1*sqrt3
 120                       if (oiandj) max = i
                           do 210 j = minj , max
                              go to (130,140,200,200,
     +                               150,200,200,160,200,200,
     +                               170,200,200,180,200,200,
     +                               200,200,200,190,
     +                               192,200,200,194,200,200,200,200,
     +                               200,196,200,200,198,200,200),j
 130                          dum2 = dum1*csj
                              if (odoub) then
                                 if (i.gt.1) then
                                    dum2 = dum2 + csi*cpj*fac
                                 else
                                    dum2 = dum2 + dum2
                                 end if
                              end if
                              go to 200
 140                          dum2 = dum1*cpj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 150                          dum2 = dum1*cdj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 160                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 170                          dum2 = dum1*cfj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 180                          if (onorm) dum2 = dum2*sqrt5
                              go to 200
 190                          if (onorm) dum2 = dum2*sqrt3
                              go to 200
 192                          dum2 = dum1*cgj
                              if (odoub) dum2 = dum2 + dum2
                              go to 200
 194                          if (onorm) dum2 = dum2*sqrt7
                              go to 200
 196                          if (onorm) dum2 = dum2*sqrt5/sqrt3
                              go to 200
 198                          if (onorm) dum2 = dum2*sqrt3
 200                          nn = nn + 1
                              dij(nn) = dum2
 210                       continue
 220                    continue
c
c     ----- overlap and kinetic energy
c
                        t = dsqrt(aa1)
                        t1 = -two*aj*aj*t
                        t2 = -pt5*t
                        p0 = ax
                        q0 = ay
                        r0 = az
                        in = -5
                        do 240 i = 1 , lit
                           in = in + 5
                           if (oradial) then
                             ni = 1
                           else
                             ni = i
                           endif
                           do 230 j = 1 , ljt
                              jn = in + j
                              if (oradial) then
                                nj = 1
                              else
                                nj = j
                              endif
                              call stvint
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
 230                       continue
 240                    continue
                        do 250 i = 1 , ij
                         nnx = ijx(i)
                         nny = ijy(i)
                         nnz = ijz(i)
                         pyz = qin(nny)*rin(nnz)
                         dum = pyz*pin(nnx)
                         s(i) = s(i) + dij(i)*dum
 250                    continue
c
                     end if
c ...
c ...
 390              continue
 400           continue
c
c
c     ----- set up overlap matrix
c
               max = maxj
               nn = 0
               do 420 i = mini , maxi
                  li = loci + i
                  in = (li*(li-1))/2
                  if (oiandj) max = i
                  do 410 j = minj , max
                     lj = locj + j
                     jn = lj + in
                     nn = nn + 1
                     if (onec .and. oradial) then
                       iexp=iatom
                       goto 7660
                     endif
                     if (iexpas .eq. -1) then
                       if (onec) then
                         iexp = iatom
                       else
                         iexp = icch
                       endif
                       goto 7660
                     endif
c	     write(iwr,*) nn, jn, s(nn), ovl(jn)
                     if(abs(s(nn)) .le. tols .or. onetwo) then
c
c      -----  (overlap is small), centre is placed at middle of
c             originating atoms
c
                       ctrx = pt5*(c(1,iat)+c(1,jat)) - cm(1)
                       ctry = pt5*(c(2,iat)+c(2,jat)) - cm(2)
                       ctrz = pt5*(c(3,iat)+c(3,jat)) - cm(3)
c
                     else
c
c      -----  centre is located at dipole/overlap
c
                       dum=done/s(nn)
                       sign=done
                       dovl= abs(ovl(jn))
                       if(dovl.gt.dzero) sign=dovl/ovl(jn)
                       ctrx=dipx(jn)*dum*sign
                       ctry=dipy(jn)*dum*sign
                       ctrz=dipz(jn)*dum*sign
c
                     endif
c
c	     write (iwr,*) jn, ctrx, ctry, ctrz
                     if (each) then
c
c      -----  assign to each overlap distribution its own centre
c
c      -----  first check if overlap distribution coincides with
c             an already defined centre
c
                       do 7800, iex = 1, nexp
              dist2 = (xexp(1,iex)-ctrx)**2  
     2               +(xexp(2,iex)-ctry)**2
     1               +(xexp(3,iex)-ctrz)**2
                         if (dist2 .le. told3) then
                           iexp = iex
                           goto 7810
                         endif
 7800                  continue
c
c      -----  add new expansion centre
c
                       nexp = nexp + 1
                       xexp(1,nexp) = ctrx
                       xexp(2,nexp) = ctry
                       xexp(3,nexp) = ctrz
                       iexp = nexp
c
 7810                  continue
c
                     else
c
c      -----  find closest expansion center
c
                       nclos = 0
            dist=sqrt((ctrx-xexp(1,1))**2+
     1       (ctry-xexp(2,1))**2+
     2       (ctrz-xexp(3,1))**2)
                       rmin=dist
                       iexp=1
c
                       do 7630 npnt=2,nexp
              dist=sqrt((ctrx-xexp(1,npnt))**2+
     1       (ctry-xexp(2,npnt))**2+
     2       (ctrz-xexp(3,npnt))**2)
                         diff = dist - rmin
                         if (abs(diff) .le. told1) 
     1             nclos = nclos + 1
                         if(diff .gt. told1)
     1             goto 7630
                         rlast=rmin
                         rmin=dist
                         iexp=npnt
                         if (rmin-rlast .lt. -told1)
     1               nclos = 0
 7630                  continue
                       repeat = nclos .ne. 0
                       if (repeat) then
c           if (repeat .and. .not. onec) then
c
c        -----  two or more expansion centra are about equally
c               distant from the centre of the overlap distribution
c               create a new one if defnew = .true.
c
                         if (defnew) then
                           nexp = nexp + 1
                           xexp(1,nexp) = ctrx
                           xexp(2,nexp) = ctry
                           xexp(3,nexp) = ctrz
                           iexp = nexp
                         else if (iexpas .eq. 0) then
c
c          -----  assign ambiguous distribution to centre of charge
c
                           iexp = icch
                         endif
                       endif
c
c      -----  one-centre overlap distributions are assigned to the
c             origin of the basis functions
c
                     endif
c
                     if (nambpt(katom(ii)) .ne. 0) then
                       iexp = katom(ii)
                     endif
                     if (nambpt(katom(jj)) .ne. 0) then
                       iexp = katom(jj)
                     endif
c
 7660                iexpc(jn)=iexp
                     if (out) then
               write(6,9966) li,lj,dum,ctrx,ctry,ctrz,iexp
 9966          format(2i4,4e15.6,i4)
                     endif
c
 410              continue
 420           continue
 430        continue
 440     continue
c
       else
         write(iwr,*) 'ASSIGN: option not implemented'
         call caserr('ASSIGN: option not implemented')
       endif
c
c-----  transform -xexp- to global origin
c
      do 9100 i=1,nexp
        do 9100 l=1,3
          xexp(l,i)=xexp(l,i)+cm(l)
 9100 continue
cahv
      return
      end
