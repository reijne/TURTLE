c 
c  $Author: jmht $
c  $Date: 2012-01-23 16:07:50 +0100 (Mon, 23 Jan 2012) $
c  $Locker:  $
c  $Revision: 6250 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd1.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   sorted6 = newmrd1    =
c ******************************************************
c ******************************************************
c
      subroutine sortdriver(q)
      implicit real*8 (a-h,o-z)
      common/symchk/ crtsym,excit,nrep,lsymd,lsadap
      logical lsymd,lsadap
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
      logical osortmr,osort2,odontgened6
      integer mapeirem,ikyrem,ninnexrem,multrem
      integer newlz
      common/mrdcisort/osortmr,osort2,mapeirem(201),
     +       ikyrem(maxorb),
     +       ninnexrem,multrem(8,8),newlz(maxorb),
     +       odontgened6
c
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
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
      dimension q(*)
c     initally allocate buffer arrays for fort31
c     generation in sorted6
c
      nx1 = nx + 1
      iston    = 1
      isacoul  = iston + 2000
      isaexc   = isacoul  + nx1
      isacoulo = isaexc   + nx1
      isaexco  = isacoulo + nx1
      ibob     = isaexco  + nx1
      need     = ibob     + nx1
c
      iston = igmem_alloc(need)
      isacoul  = iston + 2000
      isaexc   = isacoul  + nx1
      isacoulo = isaexc   + nx1
      isaexco  = isacoulo + nx1
      ibob     = isaexco  + nx1
c
c     now allocate remaining memory
c
      isorco = igmem_alloc_all(maxint)
      call setparmssym
      if (lsadap) call symsap2(q(isorco))
      if (lsymd) call sym12
c
      lfile=m6file
      do i=1,lfile
       lotape(i) = m6tape(i)
       liblk(i)  = m6blk(i)
       llblk(i)  = liblk(i) - m6last(i)
      enddo
      call sorted6(q(iston), q(isacoul), q(isaexc),
     +   q(isacoulo), q(isaexco), q(ibob), nx1,
     +   q(isorco),q(isorco),maxint)
      call gmem_free(isorco)
      call gmem_free(iston)
      return
      end

      subroutine setparmssym
      implicit real*8 (a-h,o-z)
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical osortmr,osort2,odontgened6
      integer mapeirem,ikyrem,ninnexrem,multrem
      integer newlz
      common/mrdcisort/osortmr,osort2,mapeirem(201),
     +       ikyrem(maxorb),
     +       ninnexrem,multrem(8,8),newlz(maxorb),
     +       odontgened6
c
      common/table/nirr,mult(8,8),ifrep,nirr8,jfrep
      logical modei
      integer iro, iky, mapei, mapie
      common/mappr/mapei(maxorb),modei,mapie(maxorb),
     +             iro(maxorb),iky(maxorb)
      logical lsymd,lsadap
      common/symchk/crtsym,excit,nrep,lsymd,lsadap
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
c
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/craypk/nact,nacta(maxorb),nprina,
     *              mact,macta(maxorb),mprina, ihamp,iseca
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,core,ncolo,nbas,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
      dimension msym(36)
c
      data isecb/470/
      data msym/
     *1,
     *2,1,
     *3,4,1,
     *4,3,2,1,
     *5,6,7,8,1,
     *6,5,8,7,2,1,
     *7,8,5,6,3,4,1,
     *8,7,6,5,4,3,2,1/
c
      nav = lenwrd()
      call secget(isecb,1005,iblka)
      call readi(nact,mach(14)*nav,iblka,idaf)
      ninnexrem=nact
      isec3=jdsec
      call isquar(multrem,msym,8,8)
      do ir=1,maxorb
         mapei(ir)=ir
         mapie(ir)=ir
         kr=ir*(ir-1)/2
         iky(ir)=kr
      enddo

      oprint(28)=.false.
      nsect2 = m6file
      do 20 i = 1 , nsect2
        notap2(i) = m6tape(i)
        iblk2(i) = m6blk(i)
        lblk2(i) = m6last(i)
 20   continue
      do ir=1,8
        do jr=1,8
           mult(ir,jr)=multrem(ir,jr)
        enddo
      enddo
       ninnex=ninnexrem
       crtsym=0.1d-4
      call setsto(maxorb,0,iro)
cRH
      return
      end

      subroutine sorted6(ston, sacoul, saexc,
     +  sacoulold, saexcold, bob,len1e,q,iq,maxint)
      implicit real*8 (a-h,o-z)
      integer *4 iq
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      logical osortmr,osort2,odontgened6
      integer mapeirem,ikyrem,ninnexrem,multrem
      integer newlz
      common/mrdcisort/osortmr,osort2,mapeirem(201),
     +       ikyrem(maxorb),
     +       ninnexrem,multrem(8,8),newlz(maxorb),
     +       odontgened6
c
c
      logical debugh, debugl
      common /mrdcid/ debugh, debugl
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
      common/bloksort/ibb,numw
      common/blkin/g(510),nnn
      common/craypk/i205(1360)
      integer  ntape
      common /ftap/ ntape
c
c core 
c
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
c
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)

      dimension q(*),iq(*)
c
c fort31 arrays
c
      real*8 ston, sacoul, saexc, sacoulold, saexcold
      real*8 bob
      dimension ston(2000), sacoul(len1e), saexc(len1e)
      dimension sacoulold(len1e),saexcold(len1e)
      dimension bob(len1e)
c
c fort31 
c
      real*8 core,etot,ecore
      integer jsym
      dimension jsym(36)
      integer nirr(8)
      integer idump(4),isyml(4)
      dimension idump22(4)
      character*3 name
      integer symlabeq
      integer symlabel(4)
      integer id(4)
c
      real*8 q2
      integer nit, iflop, isort, noi
      integer isym, nitlabels, imo, newl
c old mrdci
      integer idk
      common /junk/q2(512),
     +             iflop(667), isort(667), noi(667), nit(667),
     +             isym(8,maxorb), imo(maxorb), newl(maxorb),
     +             nitlabels(667), idk(maxorb+1)
c
      common/table/ nirrep,mult(8,8),ifrep,nirr8,jfrep
      logical modei
      integer iro, iky, mapei, mapie
      common/mappr/mapei(maxorb),modei,mapie(maxorb),
     +             iro(maxorb),iky(maxorb)
      common/lsort/iijjnk(4*maxorb),
     *pop(maxorb),potn,coren,ncolo,nbas,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nval,mapaie(maxorb),ilifa(maxorb),
     *iqsec,jdsec
c
      external fget
c
      lind(i,j)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
c
c data init.
c
      symcheck=0.0d0
      top=cpulft(1)
      write(iwr,1100)ntape,top
1100  format(//104('-')//
     *' **** sort of integrals and generation of fort.',i2,
     +       ' called at', f13.3,' secs')
      if (osort2) then 
         write(iwr,1301)ntape
1301  format(' **** generation of fort.',i2,' for MRDCI program')
      else 
         write(iwr,1302)ntape
1302  format(
     * ' **** generation of fort.',i2,
     + ' for new (C/C++) MRDCI program')
      endif

      do i=1,8
         do j=1,i
          jsym(lind(i,j))=mult(i,j)
         enddo
      enddo

      do i=1,8
         nirr(i)=0
         do j=1,maxorb
            if (iro(j).eq.i) nirr(i)=nirr(i)+1
         enddo
      enddo

      nit(1)=0
      do ikl=1,667
         nit(ikl)=0
      enddo
*----------------------------------------------------------------------*
*     Construct IDK                                                    *
*----------------------------------------------------------------------*
*
      idk(1) = 0
      Do i=2,maxorb
         idk(i) = idk(i-1) + i - 1
      End Do
      call vclr(sacoulold,1,len1e)
      call vclr(saexcold,1,len1e)
      call vclr(bob,1,len1e)
      ibb=0
c
      write(iwr,1101)nirrep,(in,nirr(in),in=1,nirrep)
1101  format(/
     * '   number of irreducible representations  ',i10//
     * '   irrep.    no. of mos'/
     * '   ===================='/
     * 8(I7,I12,/))

      nmo=0
      do i=1,nirrep
        nmo=nmo+nirr(i)
      enddo
      do ii=1,nmo
            imo(ii)=iro(ii)
      enddo
      ehf=0.0d0
      etot=ehf
c
c  generate new mo labels
c
      ionl=0
      do ii=1,nirrep
         ind=0
         do iii=1,nmo
            if (imo(iii).eq.ii) then
               ind=ind+1
               isym(ii,ind)=iii
               newl(iii)=ind
               newlz(iii)=newl(iii)+ionl
            endif
         enddo
         ionl=ionl+nirr(ii)
      enddo
c
      write(iwr,1200)
1200  format(/
     *'  MO    symmetry  new label '/
     *' ========================== ')

      do ikljo=1,nmo
         write(iwr,1201)ikljo,imo(ikljo),newlz(ikljo)
      enddo
1201  format(i4,i12,i11)
c
c generate nitlabels
c
      init=1
      nobr=0
      maxl=0
      nx=0
      do i=1,nirrep
         do j=1,i
            do k=1,i
               if (i.eq.k) then
                  do l=1,j
                     ij4 = jsym(lind(i,j))
                     kl4 = jsym(lind(k,l))
                     if (jsym(lind(ij4,kl4)).eq.1) then
                         idump(1)=i
                         idump(2)=j
                         idump(3)=k
                         idump(4)=l
                         init=init+1
                         numi=numberofint(i,j,k,l,nirr)
                         nx=nx+numi
                         nit(init)=nit(init-1)+numi
                         nobr=nobr+1
                         call packidump(igab,idump)
                         nitlabels(nobr)=igab
                         noi(nobr)=0
                         iflop(nobr)=numi
                         isort(nobr)=0
                         if (numi.gt.maxl) maxl=numi
                     else
                         numi=0
                         init=init+1
                         nit(init)=nit(init-1)+numi
                     endif
                  enddo
               endif
               if (i.ne.k) then
                  do l=1,k
                     ij4 = jsym(lind(i,j))
                     kl4 = jsym(lind(k,l))
                     if (jsym(lind(ij4,kl4)).eq.1) then
                         idump(1)=i
                         idump(2)=j
                         idump(3)=k
                         idump(4)=l
                         init=init+1
                         numi=numberofint(i,j,k,l,nirr)
                         nx=nx+numi
                         nit(init)=nit(init-1)+numi
                         nobr=nobr+1
                         call packidump(igab,idump)
                         nitlabels(nobr)=igab
                         iflop(nobr)=numi
                         noi(nobr)=0
                         isort(nobr)=0
                         if (numi.gt.maxl) maxl=numi
                     else
                         numi=0
                         init=init+1
                         nit(init)=nit(init-1)+numi
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
c
c no. of buckets is nobr
c largest block contains maxl integrals
c 
c     iz=maxl+(maxl/2)+mod(maxl,2)
c
c  calculate bucketsize and blocksize numw
c
      nobr2=nobr
      if (nobr2.le.1) nobr2=2
      ibsi=maxint/nobr2
      maxinb=ibsi*2/3-1
      numw=ibsi
c
      iav=maxint-numw-1
      npass=maxl/iav+1
      if(debugl) then
       write(iwr,*)'nobr2,maxint,ibsi,maxinb,numw,iav',
     +              nobr2,maxint,ibsi,maxinb,numw,iav
      endif
      write(iwr,1102)npass
1102  format (/'  number of sort passes ',i4/)
c
      name="sor"
      call opencc(name,3,2,ier)
      if (ier.ne.0)call caserr("during open sor file")
      iblok=1
      nint=1
      call setsto(1360,0,i205)
c
      do 20 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
 10   lbl=lbl+1
      call get(g(1),mw)
      if(mw.eq.0) go to 20
      if(lbl.ne.0) call find(iunit)
      int4=1
      call unpack(g(num2e+1),lab816,i205,numlab)
      do ik1 = 1, nnn
      idump22(2) = i205(int4)
      idump22(1) = i205(int4+1)
      idump22(4) = i205(int4+2)
      idump22(3) = i205(int4+3)
      int4 = int4 + 4
      if(debugh) then
       write(iwr,90000) ik1,idump22,g(ik1)
90000 format(1x,'sorted6: nword,i,j,k,l,g = ', i3,
     +   2x,4i5,3x,f20.9)
      endif
         is1=max(imo(idump22(1)),imo(idump22(2)))
         is2=min(imo(idump22(1)),imo(idump22(2)))
         is3=max(imo(idump22(3)),imo(idump22(4)))
         is4=min(imo(idump22(3)),imo(idump22(4)))
         if (is1.gt.is3) then
             isyml(1)=max(is1,is2)
             isyml(2)=min(is1,is2)
             isyml(3)=max(is3,is4)
             isyml(4)=min(is3,is4)
         else if ((is1.eq.is3).and.(is2.gt.is4))then
             isyml(1)=max(is1,is2)
             isyml(2)=min(is1,is2)
             isyml(3)=max(is3,is4)
             isyml(4)=min(is3,is4)
         else if ((is1.eq.is3).and.(is2.lt.is4))then
             isyml(1)=max(is3,is4)
             isyml(2)=min(is3,is4)
             isyml(3)=max(is1,is2)
             isyml(4)=min(is1,is2)
         else
             isyml(1)=max(is3,is4)
             isyml(2)=min(is3,is4)
             isyml(3)=max(is1,is2)
             isyml(4)=min(is1,is2)
         endif
         call packidump(isymz,isyml)
         do ik2=1,nobr
             igab=nitlabels(ik2)
             ip=ik2
             if (igab.eq.isymz) goto 15
         enddo
         if (dabs(g(ik1)).gt.symcheck)symcheck=dabs(g(ik1))
         goto 3
15       iad=ibsi*(ip-1)
         noi(ip)=noi(ip)+1
         iad=iad+noi(ip)
c
c
         symlabeq=isymz
         do ikl2=1,4
          symlabel(ikl2) = 0
         enddo
         call upackidump(symlabeq,symlabel)
c
         isyml(1)=imo(idump22(1))
         isyml(2)=imo(idump22(2))
         isyml(3)=imo(idump22(3))
         isyml(4)=imo(idump22(4))
         int=1
         do ikl2=1,4
            do ikl3=1,4
               if (symlabel(ikl2).eq.isyml(ikl3)) then
                  id(int)=newl(idump22(ikl3))
                  isyml(ikl3)=0
                  int=int+1
                  goto 5
               endif
            enddo
   5        continue
         enddo
c
         ix=symlabel(1)
         jx=symlabel(2)
         kx=symlabel(3)
         lx=symlabel(4)
         if ((ix.eq.jx).and.(kx.eq.lx).and.(ix.eq.kx)) then
                  call case1o(id,index)
                  if (osort2) then
                  if ( id(1).eq.id(2) .and. id(3).eq.id(4) )then
                     iaa=0
                     do icc=1,ix-1
                        iaa=iaa+nirr(icc)
                     enddo
                     ibab=iaa
                     iaa=iaa+id(1)
                     ibab=ibab+id(3)
                     sacoulold(idk(iaa)+ibab)=g(ik1)
                     if ((idk(iaa)+ibab).gt.len1e) 
     1                    call caserr('dimensioning problem')
                  endif
                  if ( id(1).eq.id(3) .and. id(2).eq.id(4) )then
                     iaa=0
                     do icc=1,ix-1
                        iaa=iaa+nirr(icc)
                     enddo
                     ibab=iaa
                     iaa=iaa+id(1)
                     ibab=ibab+id(2)
                     saexcold(idk(iaa)+ibab)=g(ik1)
                     if ((idk(iaa)+ibab).gt.len1e) 
     1                    call caserr('dimensioning problem')
                  endif
                  endif
             else if ((ix.eq.kx).and.(jx.eq.lx).and..not.(ix.eq.jx))
     1       then
                  call case3o(id,nirr(jx),index)
                  if (osort2) then
                  if ( id(1).eq.id(3) .and. id(2).eq.id(4) )then
                     iaa=0
                     ibab=0
                     do icc=1,ix-1
                        iaa=iaa+nirr(icc)
                     enddo
                     do icc=1,jx-1
                        ibab=ibab+nirr(icc)
                     enddo
                     iaa=iaa+id(1)
                     ibab=ibab+id(2)
                     saexcold(idk(iaa)+ibab)=g(ik1)
                     if ((idk(iaa)+ibab).gt.len1e) 
     1                    call caserr('dimensioning problem')
                  endif
                  endif
             else if ((ix.eq.jx).and.(kx.eq.lx).and..not.(ix.eq.kx))
     1       then
                  call case2o(id,nirr(kx),index)
                  if (osort2) then
                  if ( id(1).eq.id(2) .and. id(3).eq.id(4) )then
                     iaa=0
                     ibab=0
                     do icc=1,ix-1
                        iaa=iaa+nirr(icc)
                     enddo
                     do icc=1,kx-1
                        ibab=ibab+nirr(icc)
                     enddo
                     iaa=iaa+id(1)
                     ibab=ibab+id(3)
                     sacoulold(idk(iaa)+ibab)=g(ik1)
                     if ((idk(iaa)+ibab).gt.len1e) 
     1                    call caserr('dimensioning problem')
                  endif
                  endif
             else
                  call case4o(id,nirr(jx),nirr(kx),
     1            nirr(lx),index)
         endif
         q(iad)=g(ik1)
         iad=2*maxinb+2*(ibsi*(ip-1))+noi(ip)
         iq(iad)=index
         if(noi(ip).ge.maxinb) then
             iad=ibsi*(ip-1)+1
             call putsort(q(iad),q(iad),
     1                    noi(ip),isort(ip))
             noi(ip)=0
         endif
3        continue
      enddo
      iblok=iblok+1
      goto 10
20    continue
c
c write remaining blocks on sort file
c
      do ip=1,nobr
         if (noi(ip).ne.0) then
            iad=ibsi*(ip-1)+1
            call putsort(q(iad),q(iad),
     1                   noi(ip),isort(ip))
            noi(ip)=0
         endif
      enddo
c
c now write first and second record of fort.31
c
      write(iwr,1103)symcheck
1103  format(/' *** Symmetry check : Largest forbidden integral : ',
     1e10.5)
      call secini(ibl3d,idaf)
      call secget(isec3,1004,jblkk)
      write(iwr,1105)isec3,ibl3d,yed(idaf)
1105  format(//,
     *' transformed 1-electron integrals restored from section',i4,
     *' of dumpfile starting at block',i8,' of ',a4)
      call rdedx(pop,mach(15),jblkk,idaf)
      jblkk=jblkk+lensec(mach(15))
      repel=potn
      core=coren
      n1e=0
c
      call firstsort(nirrep,nirr,jsym,repel,etot,ntape)
      call secrec(nit,nx,nirr,nirrep,ntape)
c
c now go back over sort file, sort integrals and write to fortran
c unit ntape
c
      call vclr(q,1,maxint)
      iad=maxint-numw
      if (iav.gt.iad)call caserr('core error maxl')
      iadl=2*maxinb+2*iad
      ioff=0
      iston=0
      nrecs = 2
      if (nobr.eq.1) goto 115
      do ip=1,nobr
         do ikz=1,npass
             iorig=isort(ip)
 60          if (iorig.eq.0) goto 65
             call getsort(q(iad),q(iad),numi,iorig)
             do i=1,numi
                if (((iq(iadl+i-2)-(ikz-1)*iav).le.iav).and.
     1              ((iq(iadl+i-2)-(ikz-1)*iav).gt.0)) then
                   q(ioff+iq(iadl+i-2)-(ikz-1)*iav)=q(iad+i-1)
                endif
             enddo
             if (iorig.ne.0)goto 60
 65          if (iflop(ip).lt.iav) then
                 nqint=iflop(ip)
             else
                 nqint=iav
             endif
             iflop(ip)=iflop(ip)-iav
             if (nqint.lt.0) nqint=0
             do ilo=1,nqint
                iston=iston+1
                ston(iston)=q(ioff+ilo)
                if (iston.eq.2000) then
                   write(ntape)ston
                   if(debugh) then
                    write(iwr,*)'STON: iston = ', iston
                    write(iwr,90001)ston
90001               format(1x,7f12.8)
                   endif
                   nrecs = nrecs + 1
                   iston=0
                endif
             enddo
             call vclr(q(ioff+1),1,nqint)
         enddo
      enddo
      if (iston.ne.0) then
        nrecs = nrecs + 1
        write(ntape)ston
        if(debugh) then
          write(iwr,*)'STON: iston = ', iston
          write(iwr,90001)(ston(loop),loop=1,iston)
        endif
      endif
c
c ... sort integrals for no symmetry
c
      goto 116
115   continue
      call vclr(ston,1,2000)
      ntotalblocks=isort(1)
      iston=0
      nston=1
      ithisblock=0
 120  continue
      ithisblock=ithisblock+1
      iorig=ithisblock
      call getsort(q(iad),q(iad),numi,iorig)
      do i=1,numi
           if ((iq(iadl+i-2)-(nston-1)*2000) .gt.2000) then
               nston=nston+1
               write(ntape)ston
               nrecs = nrecs + 1
               if(debugh) then
                write(iwr,*)'STON:'
                write(iwr,90001)ston
               endif
               call vclr(ston,1,2000)
               iston=iq(iadl+i-2)-(nston-1)*2000
               ston(iq(iadl+i-2)-(nston-1)*2000)=q(iad+i-1)
           else
               iston=iq(iadl+i-2)-(nston-1)*2000
               ston(iq(iadl+i-2)-(nston-1)*2000)=q(iad+i-1)
           endif
      enddo
      if (ithisblock.ne.ntotalblocks)goto 120
      write(ntape)ston
      if(debugh) then
       write(iwr,*)'STON:'
       write(iwr,90001)ston
      endif
      nston=nston+1
      nrecs = nrecs + 1
116   continue
c
      write(iwr,61000) ntape, nrecs
61000 format(/1x,
     + 'total number of 2e-integral records written to fort.',i2,
     + ' = ',i5)
c
c 2 electron integral sort complete, now do 1 electron integral sort
c
      n1e=nmo*(nmo+1)/2
      call vclr(q,1,maxint)
      call rdedx(q,n1e,jblkk,idaf)
c
c  sort 1e int.
c
      ns=1
      do n=1,nirrep
         do ii=1,nirr(n)
            do jj=1,ii
               ilx=isym(n,ii)
               jlx=isym(n,jj)
               if (lind(ilx,jlx).gt.n1e) then
                   print*,'lind,n1e::',lind(ilx,jlx),n1e
                   call caserr('1e sort')
               endif
               bob(ns)=q(lind(ilx,jlx))
               if(debugh) then
                write(iwr,90002) bob(ns),ilx,jlx,ns
90002           format(' bob: ilx,jlx, ns = ',f12.7,2x,3i8)
               endif
               ns=ns+1
            enddo
          enddo
      enddo
c
      write(ntape)bob
      if (osort2) then
          write(ntape)sacoulold
          write(ntape)saexcold
      else
          write(ntape)sacoul
          write(ntape)saexc
      endif
c     ecore=etot+core
      ecore=core
      write(ntape)ecore
      call closecc(2)
      call delcc("sor",3,ier)
      if (ier.ne.0) call caserr("during delcc of sor file")
      if (ionsv4.eq.0) goto 11
c
c x,y,z integrals
c
      mtape=32
      m2=2
      call secget(ionsv4,m2,iblkx)
      lenl=lensec(nbas4*(nbas4+1)/2)
      iblkx=iblkx+1+3*lenl
      write(mtape)nirrep
      write(mtape)(nirr(n),n=1,nirrep)
c     do ijk=1,nirrep
c        write(mtape)(mult(ijk,kll2),kll2=1,nirrep)
c     enddo
      do ijk=1,3
      call vclr(q,1,maxint)
      call vclr(bob,1,len1e)
      call rdedx(q,n1e,iblkx,idaf)
c     
c  sort 1e int.
c     
c find out the irrep of this operator
      iop=-1
      dmaxint=-1000.0d0
      do n=1,nirrep
        do n2=1,n
           do ii=1,nirr(n)
              do jj=1,nirr(n2)
                 ilx=isym(n,ii)
                 jlx=isym(n2,jj)
                 if (dabs(q(lind(ilx,jlx))).gt.dmaxint) then
                    dmaxint=dabs(q(lind(ilx,jlx)))
                    iop=mult(n,n2)
                    ijp=ieor(n-1,n2-1)
                    iopx=2**ijp
                 endif
              enddo
           enddo
        enddo
      enddo
      ns=1
      do n=1,nirrep
       do n2=1,n
         if (mult(n,n2).eq.iop) then
           do ii=1,nirr(n2)
             if (n.eq.n2) then
              jjend=ii
             else
              jjend=nirr(n)
             endif
             do jj=1,jjend
               ilx=isym(n2,ii)
               jlx=isym(n,jj)
               if (lind(ilx,jlx).gt.n1e) then
                   print*,'lind,n1e::',lind(ilx,jlx),n1e
                   call caserr('1e sort')
               endif
               bob(ns)=q(lind(ilx,jlx))
               if(debugh) then
                write(iwr,90003) bob(ns),ilx,jlx,ns
90003           format(' bob: ilx,jlx, ns = ',f12.7,2x,3i8)
               endif
               ns=ns+1
             enddo
           enddo
         endif
       enddo
      enddo
c
       ns=ns-1
      write(mtape)ns,iopx
      write(mtape)(bob(ijkl),ijkl=1,ns)
      iblkx=iblkx+lenl
      enddo
c end x,y,z integrals
      close(mtape)
  11  continue
      top=cpulft(1)
      write(iwr,1104)top
1104  format(//,
     *' **** end of sort of integrals and generation of fort.31 at',
     * f13.3,' secs'//104('-'))
      end

      subroutine getsort(q,iq,numi,iblok)
      implicit real*8 (a-h,o-z)
      integer *4 iq
      common/bloksort/ibb,numw
c
      logical debugh, debugl
      common /mrdcid/ debugh, debugl
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension q(*),iq(*)
      call srchccs(2,iblok,ierror,numw)
      if (ierror.ne.0) call caserr('during search sor')
      call getccs(2,q,ierror,numw)
      if (ierror.ne.0) call caserr('during getccs')
      call upack2(q(numw),numi,iblok)
      if(debugh) then
       write(iwr,*)'getsort: num,iblok = ', numi, iblok
       write(iwr,10)(q(i),i=1,numi)
10     format(1x,6f12.7)
      endif
      return
      end

      subroutine putsort(q,iq,numi,iblok)
      implicit real*8 (a-h,o-z)
      integer *4 iq
      common/bloksort/ibb,numw
c
      logical debugh, debugl
      common /mrdcid/ debugh, debugl
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension q(*),iq(*)
      ibb=ibb+1
      q(numw) = pack2(numi,iblok)
      if(debugh) then
       write(iwr,*)'putsort: num,iblok = ', numi, iblok
      endif
      call putccs(2,q,ierror,numw)
      if (ierror.ne.0)call caserr('write error during putccs')
      iblok=ibb
      return
      end

      subroutine case1o(ijkl,ind)
      implicit real*8 (a-h,o-z)
      integer i,j,k,l
      integer ijkl(4)
c
      ia(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      ib(k,l)=max(k,l)*(max(k,l)-1)/2+min(k,l)
      index1(i,j,k,l)=max(ia(i,j),ib(k,l))*
     1 (max(ia(i,j),ib(k,l))-1)/2+min(ia(i,j),ib(k,l))
c
      if (ijkl(1).lt.0.or.ijkl(2).lt.0.or.ijkl(3).lt.0.or.ijkl(4).lt.0)
     1   call caserr("in sorted someone less than 0")
      ind=index1(ijkl(1),ijkl(2),ijkl(3),ijkl(4))
      return
      end


      subroutine case2o(ijkl,kx,ind)
      implicit real*8 (a-h,o-z)
      integer kx
      integer i,j,k,l,inj,inda,indb
      integer ijkl(4)
c
      inda(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      indb(k,l)=max(k,l)*(max(k,l)-1)/2+min(k,l)
      index2(i,j,k,l,inj)=(inda(i,j)-1)*inj+indb(k,l)
c
c      met inj:
c      inj=nj*(nj+1)/2
c      nj=number of mos in irrep j (B)
c
      inj=kx*(kx+1)/2
      if (ijkl(1).lt.0.or.ijkl(2).lt.0.or.ijkl(3).lt.0.or.ijkl(4).lt.0)
     1   call caserr("in sorted someone less than 0")
      ind=index2(ijkl(1),ijkl(2),ijkl(3),ijkl(4),inj)
      return
      end

      subroutine case3o(ijkl,jx,ind)
      implicit real*8 (a-h,o-z)
      integer jx,inb
      integer i,j,k,l
      integer ijkl(4)
c
      ia3(i,j,inb)=(i-1)*inb+j
      ib3(k,l,inb)=(k-1)*inb+l
      index3(i,j,k,l,inb)=max(ia3(i,j,inb),ib3(k,l,inb))*
     1 (max(ia3(i,j,inb),ib3(k,l,inb))-1)/2+min(ia3(i,j,inb),
     2  ib3(k,l,inb))
c
c     met inb:
c     inb = number of mos in irrep b
      if (ijkl(1).lt.0.or.ijkl(2).lt.0.or.ijkl(3).lt.0.or.ijkl(4).lt.0)
     1   call caserr("in sorted someone less than 0")
      ind=index3(ijkl(1),ijkl(2),ijkl(3),ijkl(4),jx)
      return
      end

      subroutine case4o(ijkl,jx,kx,lx,ind)
      implicit real*8 (a-h,o-z)
      integer jx,kx,lx
      integer i,j,k,l
      integer nj,nk,nl
      integer ijkl(4)
c
      index4(i,j,k,l,nj,nk,nl)=(i-1)*nj*nk*nl+
     1 (j-1)*nk*nl+(k-1)*nl+l
c
c      nj=number of mos in irrep (B)
c      nk=number of mos in irrep (C)
c      nl=number of mos in irrep (D)
c
      if (ijkl(1).lt.0.or.ijkl(2).lt.0.or.ijkl(3).lt.0.or.ijkl(4).lt.0)
     1   call caserr("in sorted someone less than 0")
      ind=index4(ijkl(1),ijkl(2),ijkl(3),ijkl(4),jx,kx,lx)
      return
      end

      function numberofint(i,j,k,l,nirr)
      implicit real*8 (a-h,o-z)
      integer nirr(8)
      if ((i.eq.j).and.(k.eq.l).and.(j.eq.k)) then
c case 1
          ix=nirr(i)
          nz=ix*(ix+1)/2
          nn=nz*(nz+1)/2
      else if ((i.eq.j).and.(k.eq.l).and..not.(j.eq.k)) then
c case 2
          ix=nirr(i)
          kx=nirr(k)
          nn=(ix*(ix+1)/2)*(kx*(kx+1)/2)
      else if ((i.eq.k).and.(j.eq.l).and..not.(i.eq.j))then
c case 3
          ix=nirr(i)
          jx=nirr(j)
          nn=ix*jx*(ix*jx+1)/2
      else
c case 4
c     nb. this code results in i*2 truncation on the alpha
c         nn=nirr(i)*nirr(j)*nirr(k)*nirr(l)
          ix=nirr(i)
          jx=nirr(j)
          kx=nirr(k)
          lx=nirr(l)
          nn = ix * jx * kx * lx
      endif
      numberofint=nn
      return
      end

      subroutine secrec(nit,nx,lj,nocc,ntape)
      implicit real*8 (a-h,o-z)
      integer nit(666), ntape
      integer lg(8),iball(8),itill(8),lj(8)
      integer mcompl
      common/scra/mcompl(100)
      do ik=1,100
         mcompl(ik)=0
      enddo
      IX=0
      JX=0
      KX=0
      DO  I=1,NOCC
        LG(I)  =IX
        ITILL(I)=JX
        IBALL(I)=KX
        IX=IX+LJ(I)*(LJ(I)+1)/2
        JX=JX+LJ(I)*LJ(I)
        KX=KX+LJ(I)
      End Do
      write(ntape)nit,nx,lg,iball,itill,mcompl
      return
      end

      subroutine firstsort(n,ljl,jabl,repel,etot,ntape)
      implicit real*8 (a-h,o-z)
      integer ntape
      integer ljl(8),jabl(36)
c
c     note that the following are not set, but pad out ntape
c
      integer mjl,kjl,njl,ntill,nball
      integer lsyml,isyml,ncompl
      common/junk3/mjl(8),kjl(8),njl(8),ntill(8),nball(9),
     + isyml(8),lsyml(800),ncompl(100)
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
      logical osortmr,osort2,odontgened6
      integer mapeirem,ikyrem,ninnexrem,multrem
      integer newlz
      common/mrdcisort/osortmr,osort2,mapeirem(201),
     +       ikyrem(maxorb),
     +       ninnexrem,multrem(8,8),newlz(maxorb),
     +       odontgened6
c
      iorbs=0
      knu=0
      icmo=0
      nt=0
      nod=2000
      jod=2000
      nsel=0
      if (osort2) then
         nod=1000
         nsel=n
         ntill(1)=0
         do ii=1,n
            icmo=icmo+ljl(ii)
            isyml(ii)=ii
         enddo
         icimo=0
         do ilif=1,8
            icimo=icimo+ljl(ilif)
            kjl(ilif)=0
         enddo
         do ii=2,n
            ntill(ii)=ntill(ii-1)+ljl(ii-1)
         enddo
      endif
      write (ntape) n,nod,jod,nt,icmo,mjl,kjl,ljl,njl,nsel,
     +              ntill,nball,isyml,jabl,iorbs,knu,lsyml,
     +              ncompl,repel,etot,newlz
      return
      end
      subroutine packidump(igab,idump)
      implicit real*8  (a-h,o-z), integer (i-n)
      integer igab,idump
      dimension idump(4)
c
      igab =   ior(ishft(idump(1),24),
     +         ior(ishft(idump(2),16),
     +         ior(ishft(idump(3), 8),idump(4))))
c
      return
      end
      subroutine upackidump(igab,idump)
      implicit real*8  (a-h,o-z), integer (i-n)
      integer igab,idump
      dimension idump(4)
      data m8/z'ff'/
c
      idump(1) = iand(ishft(igab,-24),m8)
      idump(2) = iand(ishft(igab,-16),m8)
      idump(3) = iand(ishft(igab, -8),m8)
      idump(4) = iand(      igab,    m8)
c
      return
      end

      subroutine symsap2(cr)
c
      implicit real*8  (a-h,o-z),integer  (i-n)
c...
c...   pick up symmetries from gamess adapt/orbitals
c...
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
       parameter (mxorb1=maxorb+1)
      dimension cr(*)
c...    common/craypk/ used to receive section 490
        common/craypk/ nira,mula(8,8),isymao(maxorb),isymmo(maxorb)
c...
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c
      character *8 rname, sname, text11, text22, anumt, bnumt
      common/cic/rname(8),sname,text11(14),text22(14),anumt(8),bnumt(8)
c
        logical lsymd,lsadap
        common/symchk/ crtsym,excit,nrep,lsymd,lsadap
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      logical modei
      integer iro, iky, mapei, mapie
      common/mappr/mapei(maxorb),modei,mapie(maxorb),
     +             iro(maxorb),iky(maxorb)
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c..     next common blocks for reading of vectors
      logical otrad
      common/lsort/ilifx(maxorb),ilifs(maxorb),ilifq(maxorb),
     * lelim(maxorb),
     *comman(maxorb),pottt(2),ncolo,nbas,newb,ncore,
     *mapcie(maxorb),ilifc(maxorb),nsact,mapaie(maxorb),ilifa(maxorb)
     *,iqsec
      common/scra  /ilifd(maxorb),ntrad(maxorb),itran(mxorb3),
     * ctran(mxorb3),otrad,itrad,deig(maxorb),
     *docc(mxorb1),nbasis,newbas,ncolq,jeig,jocc,jtemp,
     *potnuc,dx,dy,dz,sone(508)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
c
c
c...  now get vector info
c
      call secget(isec3,1004,jblkk)
      call rdedx(comman,mach(15),jblkk,idaf)
c...
c...  now get eigen vectors from dumpfile
c...
      len1=lensec(mach(8))
      nbsq=nbas*nbas
      call secget(iqsec,3,jblkk)
      nav = lenwrd()
      call readi(ilifd,mach(9)*nav,jblkk+len1+1,idaf)
      iblkrl = iposun(idaf)
      call rdedx(cr,nbsq,iblkrl,idaf)
c
c...     if vectors were not symmetry adapted switch to symdet
c
      if (otrad) then
         lsadap = .false.
         lsymd = .true.
         write(iwr,560)
         return
      end if
c..
      call secget(isect(490),51,iblock)
      call readi(nira,mach(13)*nav,iblock,idaf)
c..
      do i=1,8
       do j=1,8
        mult(i,j) = mula(i,j)
       enddo
      enddo
      nirr = nira
cjvl  if nirr not 1,2,4,8, make it so 
      if (nirr.eq.3) nirr = 4
      if (nirr.eq.5.or.nirr.eq.6.or.nirr.eq.7) nirr = 8
cjvl
      nrep = nirr
c..   assign labels
      sname = anumt(nirr)
      do i=1,nirr
       rname(i) = anumt(i)
      enddo
c..
c..        build iro  (reordered)
c..
      aax = 0.0d0
      iix = 0
c
      do loop=1,ninnex
         jc=mapaie(mapie(loop))
         ilif = (jc-1)*nbas
c...     find position of absolute max. element
         ix = 1
         ax = dabs(cr(ilif+1))
         do i=2,newb
            if (dabs(cr(ilif+i)).gt.ax) then
               ax = dabs(cr(ilif+i))
               ix = i
            end if
         enddo
c
c...     symmetry of mo is symmetry of ao with largest coef
c
         iro(loop) = isymao(ix)
c
c..      check if orbital pure symmetry (i.e. find worst case)
c
         do i=1,newb
            if (isymao(i).ne.iro(loop)) then
               if (dabs(cr(ilif+i)).gt.aax) then
                  aax = dabs(cr(ilif+i))
                  iix = loop
               end if
            end if
         enddo
      enddo
c
      if(.not.oprint(28)) then
       write(iwr,559)
       do j=1,nirr
        write(iwr,557) (mult(i,j),i=1,nirr)
       enddo
      endif
      if (nirr.ne.nira) write(iwr,556) nira,nirr
c..
       if (aax.gt.1.0d-5) then
         write(iwr,562) aax,iix
      end if
c
      return
c
560   format(/' ** no adapt info .. switch to symdet **')
557   format(5x,8i3)
556   format(/' ** # representations extended from',i6,' to',i6,' ** ')
559   format(/1x,'** symmetry taken from adapt ** ',//
     1       ,7x,'multiplication table'/,
     2        7x,'====================')
562   format(/' *** largest deviation from pure symmetry of ',e12.5,
     *       ' detected in orbital ',i4,' ***')
      end
      subroutine defaults_mrdci(cr)
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
c
c     this routine is designed to:
c  1. construct indexing arrays for the core and
c     discarded MOs required by the analysis MRDCI routines
c  2. construct a reasonable
c     set of defaults for the MRDCI module given a set
c     of input orbitals from the dumpfile, plus CORE
c     and (ACTIVE) data
c
      parameter (mxorb1=maxorb+1)
      parameter (maxroot=50)
      parameter (maxshl=50)
      parameter (maxshl1=maxshl+1)
      dimension cr(*)
c...    common/craypk/ used to receive section 490
      common/craypk/nira,mula(8,8),isymao(maxorb),isymmo(maxorb),
     +              nspace,
c       and section isect(470)
     +              nact,nacta(maxorb),nprina,
     +              mact,macta(maxorb),mprina, ihamp,iseca
c...
c
      integer nirr,mult,ifrep,nirr8,jfrep
      common /table/ nirr,mult(8,8),ifrep,nirr8,jfrep
c
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
        logical lsymd,lsadap
        common/symchk/ crtsym,excit,nrep,lsymd,lsadap
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      integer isec3, nsect2, notap2, iblk2, lblk2
      common /file12/ isec3,nsect2,notap2(20),iblk2(20),lblk2(20)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      logical debugh, debugl
      common /mrdcid/ debugh, debugl
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
      parameter (maxref = 256)
c
      real*8 edavit,cdavit,extrapit,eigvalr,cradd, 
     +     weighb, rootdel, ethreshit
      integer mbuenk,nbuenk,mxroots
      logical ifbuen, ordel, odave, oweight,odavit
      common /comrjb2/edavit(maxroot),cdavit(maxroot),
     +                odavit(maxroot),
     +                extrapit(maxroot),ethreshit(maxroot),
     +                eigvalr(maxref),cradd(3),weighb,
     +                rootdel,ifbuen,mbuenk,mxroots,
     +                ordel,odave,nbuenk,oweight
c
      character*8 zjob,zdate,ztime,zprog,ztype,zspace,ztext
      common/junkc/zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
      common/atmol3/mina(2),mouta,moutb
c..  next common blocks for reading of vectors
      common/lsort/ncore,mapcie(maxorb),ilifc(maxorb),
     + nsact,mapaie(maxorb),ilifa(maxorb),iqsec,itmp1,
c    local arrays ..
     + iro(maxorb), iroc(maxorb), irod(maxorb), mapd(maxorb), 
     + numb(maxorb),ilecj(maxorb),nvmo(8),
     + nj(8),njd(8),njs(8),nsorto(8),lj(8),njpair(8),
     + njps(40),njp(40),eigval(maxorb),
     + eigvald(maxorb,8),eigvals(maxorb,8)
      logical otrad
      common/scra/ilifd(maxorb),ntrad(maxorb),itran(mxorb3),
     * ctran(mxorb3),otrad,itrad,
     + deig(maxorb),docc(mxorb1),nbasis,newbas,ncol,jeig,jocc,jtemp
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      logical ospecc, ospece, ospecs
      common/trann/ic,i7,mj(8),kj(8),isecv,
     +             icore,ick(maxorb),izer,iaus(maxorb),kjunk,
     +             ospece, ospecc, ospecs
      common/multic/radius(40),irad(57+mcfzc),itype(maxorb),isecn
c
c     mrdci common blocks
c
      logical oconf,debugs
      common /parkin/egey1,trash,tdel,
     +              ical0,nele,nko,mxex,nmulp,ispacep,nprin,
     1              nft31,nstarv,ncorci,nform,oconf,debugs,
     +              lsng,maxci,ipt0,isym(maxroot)
c
      logical debugp,ospecp,oprop1e,noscfp,doscf
      common/cprop1e/istate,ipig(maxroot),ipmos,ipaos,iaopr(11),
     +               imopr(11),jkonp(maxshl1*maxroot),debugp,ospecp
     +               ,oprop1e,noscfp,doscf
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp
c
      logical out
c
      data m29,isecb/29,470/
      character *8 zcas
      data zcas/'mcscf'/
c
c     clear core and discarded arrays for new mrdci
c
      nav = lenwrd()
      out = nprint.ne.-5
c
      write(iwr,8912)
8912  format(/1x,'pre-processing mrd-ci data')
c
      nsomos = 0
      ic = 0
      i7 = 0
      do loop = 1, 8
       mj(loop) = 0
       kj(loop) = 0
       lj(loop) = 0
       njd(loop) = 0
       njs(loop) = 0
       njpair(loop)=0
      enddo
      icore = 0
      izer = 0
      call setsto(maxorb,0,ick)
      call setsto(maxorb,0,iaus)
c
c     first, restore header blocks from vectors section
c
      isecv=mouta
      if(moutb.ne.0)isecv=moutb
      if(zscftp.eq.zcas) iscftp = 1
      if(iscftp.eq.1.and.isecn.gt.0)isecv=isecn
c
      call secget(isecv,3,jblkk)
      call rdchr(zjob,m29,jblkk,idaf)
      call reads(deig,mach(8),idaf)
      len1=lensec(mach(8))
      len2=lensec(mach(9))
      iblki = jblkk+len1+1
      iblkv=  iblki+len2
      if(nbasis.ne.num) call caserr(
     *'restored eigenvectors are in incorrect format')
      nbsq=nbasis*nbasis
      if(out)
     +   write(iwr,89004)isecv,ztype,ztime,zdate,zjob,ztext
89004 format(/' vectors restored from section ',i4//
     *' header block information : '/
     * 1x,a7,'vectors created at ',a8,' on ',a8,
     *' in the job ',a8/
     *' with the title: ',10a8)
      if(iscftp.ne.1) then
       call readi(ilifd,mach(9)*nav,iblki,idaf)
      endif
c     now restore the vectors
      call rdedx(cr,nbasis*nbasis,iblkv,idaf)
c
      if(iscftp.eq.1) then
         call comharm(cr,'vectors',ilifq)
         ncol =newbas0
         nbasis = newbas0
      else
c     allow for spherical harmonics
         ncol = newbas0
      end if
c
c     now restore =active= and =core= specifications
c
      call secget(isecb,1005,iblka)
      call readi(nact,mach(14)*nav,iblka,idaf)
c
c     active list
c     allow for spherical harmonics
      if (nbasis.ne.newbas0) then
c     redefine active list if not specified
       if (nact.gt.ncol) then
         nact=ncol
         call wrt3i(nact,mach(14)*nav,iblka,idaf)
       endif
      endif
      nsact = nact
c
      if(nsact.lt.1.or.nsact.gt.ncol) then
       write(iwr,4001) nsact, ncol
 4001  format(/' **** invalid number of active orbitals'/
     +         ' nact, ncol = ', 2i5)
       call caserr('invalid number of active orbitals')
      endif
      do i=2,nsact
       if((locat1(nacta,i-1,nacta(i))).ne.0)
     *  call caserr(
     * 'invalid orbital specified in active list')
      enddo
      do i=1,nsact
       j=nacta(i)
       mapaie(i)=j
       ilifa(i)=ilifq(j)
      enddo
c
      if(debugl) then
       write(iwr,1005)
 1005  format(/
     * ' the following orbitals included in active list'/
     * ' =============================================='/
     * 4(' function e/i label',4x)/1x,88('-'))
       write(iwr,1006)(mapaie(k),k,k=1,nsact)
 1006  format(4(i9,i10,4x))
      endif
c
c     frozen list
c
      ncore=mact
      if(ncore.le.0) then
       if(debugl) write(iwr,2003)
 2003  format(/
     * ' no functions specified in frozen core list')
       go to 3000
      else
       if(ncore.gt.ncol)call caserr(
     *'invalid number of functions in core list')
       do i=1,ncore
       j=macta(i)
       mapcie(i)=j
       ilifc(i)=ilifq(j)
       enddo
       do i=2,ncore
       if((locat1(mapcie,i-1,mapcie(i))).ne.0) call caserr(
     * 'invalid function specified in frozen core list')
       enddo
       if (debugl) then
        write(iwr,2008)
 2008   format(/
     *  ' the following mos included in frozen core list'/
     *  ' ----------------------------------------------'/
     *  4(' function e/i label',4x))
        write(iwr,1006)(mapcie(k),k,k=1,ncore)
       endif
      endif
c
c...  if vectors were not symmetry adapted this analysis will
c     fail; we could switch to symdet, but the code will need
c     reworking here .. for the moment caserr
c
 3000 if (otrad) then
         lsadap = .false.
         lsymd = .true.
         write(iwr,560)
 560     format(/' ** no adapt info .. switch to symdet **')
         call caserr('must have symmetry adapted info. available')
      end if
c..
c     restore symmetry data 
c..
      call secget(isect(490),51,iblka)
      call readi(nira,mach(13)*nav,iblka,idaf)
c
      if(nira.le.0.or.nira.gt.8) call caserr(
     +     'symmetry data has been corrupted')
c..
        do i=1,8
         do j=1,8
         mult(i,j) = mula(i,j)
         enddo
        enddo
        nirr = nira
cjvl     if nirr not 1,2,4,8, make it so 
        if (nirr.eq.3) nirr = 4
        if (nirr.eq.5.or.nirr.eq.6.or.nirr.eq.7) nirr = 8
cjvl
        nrep = nirr
c..
c..        build iro  (reordered)
c..
      aax = 0.0d0
      iix = 0
c
      do loop=1,nsact
         jc=mapaie(loop)
         ilif = (jc-1)*nbasis
c...     find position of absolute max. element
         ix = 1
         ax = dabs(cr(ilif+1))
         do i=2,newbas
            if (dabs(cr(ilif+i)).gt.ax) then
               ax = dabs(cr(ilif+i))
               ix = i
            end if
         enddo
c
c...     symmetry of mo is symmetry of ao with largest coef
c
         iro(loop) = isymao(ix)
         eigval(loop) = deig(jc)
c
c..      check if orbital pure symmetry (i.e. find worst case)
c
         do i=1,newbas
            if (isymao(i).ne.iro(loop)) then
               if (dabs(cr(ilif+i)).gt.aax) then
                  aax = dabs(cr(ilif+i))
                  iix = loop
               end if
            end if
         enddo
      enddo
c..
      if(debugl) then
       write(iwr,559)
559    format(/'   ** symmetry taken from adapt ** ',//
     1       ,7x,'multiplication table'/,
     2        7x,'====================')
       do j=1,nirr
        write(iwr,557) (mult(i,j),i=1,nirr)
557    format(5x,8i3)
       enddo
      endif
c
      if (nirr.ne.nira) write(iwr,556) nira,nirr
556   format(/' ** # representations extended from',i6,' to',i6,' ** ')
c..
       if (aax.gt.1.0d-5) then
         write(iwr,562) aax,iix
562   format(/' *** largest deviation from pure symmetry of ',e12.5,
     *       ' detected in orbital ',i4,' ***')
      end if
c
c     characterise active mos
c
      do i = 1,nirr
       do j = 1,nsact
        if(iro(j).eq.i) then
         lj(i) = lj(i) + 1
        endif
       enddo
      enddo
c
      iix = 0
      do i =1,nirr
       iix = iix + lj(i)
      enddo
      if(iix.ne.nsact) call caserr(
     +  'problem in assigning symmetry of active m.o.s')
c
      if (debugl) then
        write(iwr,*) 'lj (active) = ', lj
        write(iwr,*) 'ncore = ', ncore
        write(iwr,*) 'mapcie = ', (mapcie(loop),loop=1,ncore)
      endif
c
c     now deduce symmetry of core and discarded m.o.s
c
      if(ncore.gt.0) then
       ic = 1
       if(out) then
        write(iwr,8082)
 8082   format(/
     +  ' *********************************'/
     +  ' ** frozen orbital specification *'/
     +  ' *********************************')
       endif
c     first handle core functions
c
      icore = ncore
      iix = 0
      aax = 0.0d0
c
      do loop=1,ncore
         jc=mapcie(loop)
         ilif = (jc-1)*nbasis
c...     find position of absolute max. element
         ix = 1
         ax = dabs(cr(ilif+1))
         do i=2,newbas
            if (dabs(cr(ilif+i)).gt.ax) then
               ax = dabs(cr(ilif+i))
               ix = i
            end if
         enddo
c
c...     symmetry of mo is symmetry of ao with largest coef
c
         iroc(loop) = isymao(ix)
c
c..      check if orbital pure symmetry (i.e. find worst case)
c
         do i=1,newbas
            if (isymao(i).ne.iroc(loop)) then
               if (dabs(cr(ilif+i)).gt.aax) then
                  aax = dabs(cr(ilif+i))
                  iix = loop
               end if
            end if
         enddo
      enddo
      do i = 1,nirr
       do j = 1,ncore
        if(iroc(j).eq.i) then
         mj(i) = mj(i) + 1
        endif
       enddo
      enddo
      nsd = 0
      do i=1,nirr
       mjj = mj(i)
       do j = 1, mjj
        nsd = nsd + 1
        ick(nsd) = j
       enddo
      enddo
       if (debugl) then
        write(iwr,*)'core, mj = ', mj
        write(iwr,*)'core, ick = ', (ick(loop),loop=1,ncore)
       endif
      if (aax.gt.1.0d-5) then
       write(iwr,563) aax,iix
563    format(/' *** largest deviation from pure symmetry of ',
     =        e12.5,' detected in core orbital ',i4,' ***')
      end if
c
      k=0
      do i=1,nirr
       if(out.and.mj(i).ne.0) then
        write(iwr,8072)i
 8072   format(/' orbital symmetry representation no.',i3)
        iij=k+1
        k=k+mj(i)
        write(iwr,8076)mj(i)
 8076   format(/' no. of frozen orbitals ',i3)
        write(iwr,8077)
        write(iwr,8078)(ick(j),j=iij,k)
 8077   format(/' orbital sequence nos.')
 8078   format(/1x,20i4)
       endif
      enddo
c
      endif
c
c     now deduce symmetry of discarded mos
c
      if (debugl) then
       write(iwr,*) 'nsact = ', nsact
       write(iwr,*) 'mapaie = ', (mapaie(loop),loop=1,nsact)
       write(iwr,*) 'ncol = ', ncol
      endif
      ndisc = newbas0 - (nsact+ncore)
      if (ndisc.gt.0) then
c
c     in most general case 
c     can only handle this by process of elimination (?)
c     as no track is kept of discarded functions
c
       nsd = 0
       do loop = 1, ncol
c
c     core functions
c
       do i = 1, ncore
       if(loop.eq.mapcie(i)) go to 600
       enddo
c
c     active functions
c
       do i = 1,nsact
       if(loop.eq.mapaie(i)) go to 600
       enddo
c
c     found an m.o. that is in neither list
c
       nsd = nsd + 1
       mapd(nsd) = loop
600    continue
       enddo
c
       if(ndisc.ne.nsd) then
        write(6,*) 'ndisc, nsd = ', ndisc, nsd
        write(6,*) 'mapd = ', (mapd(i),i=1,nsd)
        call caserr('error in discarded m.o. logic')
       endif
c
c      check that the discarded mos are not in fact stupid
c      null vectors that are added in the vectors/enter phase
c      note that it would seem that o/p'ing NOS from the
c      table-ci code is screwed if this section is specified
c      using enter ...
c
       do loop = 1,nsd
        jc = mapd(loop)
        ilif = (jc-1)*nbasis
        ix = 1
        ax = dabs(cr(ilif+1))
        do i=2,newbas
           if (dabs(cr(ilif+i)).gt.ax) then
              ax = dabs(cr(ilif+i))
              ix = i
           end if
        enddo
        if (ax.gt.1.0d-5) go to 700
       enddo
c
c      null vectors they are .. reset ndisc
c      this should not occurr given use of newbas0
c
       ncol =  nsact+ncore
       ndisc = 0
      endif
c
700    if (ndisc.gt.0) then
       i7 = 1
       if(out) then
        write(iwr,8083)
 8083   format(/
     +  ' ***********************************'/
     +  ' * discarded orbital specification *'/
     +  ' ***********************************'/)
       endif
       iix = 0
       aax = 0.0d0
c
       do loop=1,ndisc
          jc=mapd(loop)
          ilif = (jc-1)*nbasis
c...      find position of absolute max. element
          ix = 1
          ax = dabs(cr(ilif+1))
          do i=2,newbas
             if (dabs(cr(ilif+i)).gt.ax) then
                ax = dabs(cr(ilif+i))
                ix = i
             end if
          enddo
c
c...     symmetry of mo is symmetry of ao with largest coef
c
          irod(loop) = isymao(ix)
c
c..      check if orbital pure symmetry (i.e. find worst case)
c
          do i=1,newbas
             if (isymao(i).ne.irod(loop)) then
                if (dabs(cr(ilif+i)).gt.aax) then
                   aax = dabs(cr(ilif+i))
                   iix = loop
                end if
             end if
          enddo
       enddo
c
       do i = 1,nirr
        do j = 1,ndisc
         if(irod(j).eq.i) then
          kj(i) = kj(i) + 1
         endif
        enddo
       enddo
       nsd = 0
       do i=1,nirr
        ix = lj(i) + mj(i)
        mjj = kj(i)
        do j = 1, mjj
         nsd = nsd + 1
         iaus(nsd) = j + ix
        enddo
       enddo
c
       if (debugl) then
         write(iwr,*)'discarded, kj = ', kj
         write(iwr,*)'discarded, iaus = ', (iaus(loop),loop=1,ndisc)
       endif
       if (aax.gt.1.0d-5) then
        write(iwr,564) aax,iix
564     format(/' *** largest deviation from pure symmetry of ',
     =        e12.5,' detected in discarded orbital ',i4,' ***')
       end if
       l=0
       do i=1,nirr
       if(out.and.kj(i).ne.0) then
        write(iwr,8072)i
        iij=l+1
        l=l+kj(i)
        write(iwr,8081)kj(i)
 8081   format(/' no. of discarded orbitals ',i3)
        write(iwr,8077)
        write(iwr,8078)(iaus(j),j=iij,l)
       endif
       enddo
c
      endif
c
c     now check for consistency of input data. This is driven
c     through 3 variables, ospecc (conf data), ospece (cntrl or
c     nele data) and ospecs (symmetry). If these have not
c     been specified by input, construct defaults.
c
      if (ospece) then
c
c     check consistency of input data
c
       nel = ne - ncore - ncore
       if(nele.eq.nel) then
        nint = ne - nopen - npair - npair
        nint = nint/2 - ncore
         ndomos = nint
         if (nseto.ne.0) then
           do i = 1 , nseto
             nsomos = nsomos + no(i)
           enddo
         end if
         nint = ndomos + nsomos + npair + npair
       else
        if(.not.ospecc) then
         call caserr(
     + 'invalid no. of electrons for default configuration generator')
        endif
       endif
      else
c
c     construct default no. of active electrons
c
       nint = ne - nopen - npair - npair
       nint = nint/2 - ncore
       ndomos = nint 
       if (nseto.ne.0) then
        do i = 1 , nseto
           nsomos = nsomos + no(i)
        enddo
       end if
       nint = ndomos + nsomos + npair + npair
       nele = ne - ncore - ncore
       if (nele.ne.(ndomos*2+nsomos)) then
         write(iwr,*) ' no. of. DOMOS = ', ndomos
         write(iwr,*) ' no. of. SOMOS = ', nsomos
         write(iwr,*) ' nele          = ', nele  
         write(iwr,*) ' ncore         = ', ncore  
         write(iwr,*) ' npair         = ', npair  
CMR      Note: ndomos != (ne-nopen-2*npair-2*ncore)/2 for quartet and higher
CMR      Perhaps ne-nsomos-2*npair-2*ncore would be better? will force-fix for now
         ndomos=(nele-nsomos)/2
       endif
         
c
      endif
c
      if (debugl) then
       write(iwr,*) ' no. of. DOMOS = ', ndomos
       write(iwr,*) ' no. of. SOMOS = ', nsomos
      endif
c
      if(nele.lt.2.or.nele.gt.200) call caserr(
     + 'invalid number of active electrons')
c
c... determine degeneracy for nterm30
cjvl this ia not really right in case of symmetry
cjvl but if you know it so well you can specify conf in the input
c
      nterm30d = 0
      lloop = 0
      do loop=1,nsact
         if (lloop-1.gt.loop) cycle
         do lloop=loop,nsact 
            if (dabs(eigval(loop)-eigval(lloop)).gt.0.01d0) exit
         end do
         nterm30d = max(nterm30d,lloop-loop)
      end do
c
c     must consider required symmetry for reference functions
c
      if (.not.ospecs) then
c
c     determine ispacep based on input orbitals (open shell)
c
       if (nsomos.gt.0) then
        ispacep = iro(ndomos+1)
        do i = ndomos+2, ndomos+nsomos
         ispacep = mult(ispacep,iro(i))
        enddo
       else
        ispacep = 1
       endif
      endif
c
      write(iwr,8084) ispacep
 8084 format(/1x,'**** wavefunction symmetry = ', i2)
      if(ispacep.lt.1.or.ispacep.gt.8)call caserr(
     + 'invalid symmetry requested for CI wavefunction')
c
c     print attributes of both complete and
c     active orbital space
c
      write(iwr,1101)nrep
1101  format(/
     * '                ORBITAL SPECIFICATIONS                '/
     * ' ====================================================='/
     * '     number of irreducible representations  ',i6       /
     * '                  COMPLETE ORBITAL space              '/
     * ' irrep.   no. of mos     frozen    active    discarded'/
     * ' =====================================================')
      do i=1,nirr
      ntotal =  mj(i) +  lj(i) + kj(i)
      if(ntotal.gt.0) then
      write(6,1102)  i, ntotal, mj(i), lj(i), kj(i)
1102  format(i7,4i10)
      endif
      enddo
c     now print active orbital space
      ik = 0
c
      write(iwr,1103)
1103  format(
     * ' ====================================================='/
     * '                   ACTIVE ORBITAL space               '/
     * ' irrep.   no. of mos     occupied              virtual'/
     * ' =====================================================')
      do i = 1, nrep
      ik0 = ik
      if(lj(i).gt.0) then
        do j = 1, ndomos+nsomos
         if(iro(j).eq.i) then
          njd(i) = njd(i) + 1
         endif
        enddo
       if(njd(i).gt.0) then
       write(iwr,1104) i, lj(i), ik0+1, ik0+njd(i), 
     +                 ik0+njd(i)+1, ik0+lj(i)
1104   format(i7,5x,i5,7x, i3,' to ',i3,10x, i3,' to ',i3)
       else
       write(iwr,1106) i, lj(i), ik0+1, ik0+lj(i)
1106   format(i7,5x,i5,7x, 20x, i3,' to ',i3)
       endif
       ik = ik + lj(i)
      endif
      enddo
      write(iwr,1105)
1105  format(
     + ' ====================================================='/)
c
      do i=1,8
       njd(i) = 0
      enddo
c

      if (ospecc) then
c     check input configurations
       if(nko.le.0.or.nko.gt.256) call caserr(
     + 'invalid number of reference functions')
c
c     delay checks for direct-style specification until
c     selection processing
c
      if (.not.oconf) then
       irdc=11
       nshl = 80
       call rewftn(irdc)
       do i=1,nko
c       read the configurations
        read(irdc,*)np
        if(np.lt.0.or.np.gt.9) call caserr(
     +  'invalid no. of open shell mo.s. in reference function')
        j=nmulp+np
        if(j-2*(j/2).eq.0.or. np.gt.nele. or.
     +    np.gt.nmulp+3) then
         call caserr(
     +   'open shell structure inconsistent with multiplicity')
        endif
        mm = (nele+np)/2
        read(irdc,*)(numb(j),j=1,mm)
        do j = 1, mm
         if(numb(j).le.0.or.numb(j).gt.nsact) call caserr(
     +      'invalid orbital index in reference function')
        enddo
       enddo
       call rewftn(irdc)
c
       call chk_mrdci_data(nsact,nrep)
c
      endif
c
      else
c
c     construct default root set comprising the SCF
c     configuration, plus the lowest double and single
c     excitation from the highest DOMO of each symmetry
      nko = 1
c     first, categorise the MOS by MRDCI numbering
c
c     first domos
      ik = 0
      ij0 = nsomos
c
      do i = 1, nrep
      ik0 = ik
       do j = 1, ndomos
       if(iro(j).eq.i) then
        ij0 = ij0 + 1
        ik0 = ik0 + 1
        njd(i) = njd(i) + 1
        eigvald(njd(i),i) = eigval(j)
        numb(ij0) = ik0
       endif
       enddo
      ik = ik + lj(i)
      enddo
      if(debugl) then
       do j = 1,nrep
       njdj = njd(j)
       write(iwr,2546) j,njdj
       if(njdj.gt.0) write(iwr,2547) (eigvald(loop,j),loop=1,njdj)
2546   format(1x,'*** symmetry = ',i1,' no. of. DOMOS = ', i3)
2547   format(1x,'eigenvalues = '/ 1x,8f12.7)
       enddo
      endif
c     now somos
      if(nsomos.gt.0) then
       ik = 0
       ij0 = 0
       do i = 1, nrep
       ik0 = ik
        do j = ndomos+1,ndomos+nsomos
        if(iro(j).eq.i) then
         ij0 = ij0 + 1
         ik0 = ik0 + 1
         ik1 = ik0 + njd(i)
         njs(i) = njs(i) + 1
         eigvals(njs(i),i) = eigval(j)
         numb(ij0) = ik1
        endif
        enddo
       ik = ik + lj(i)
       enddo
      endif
      if(debugl) then
       do j = 1,nrep
       njds = njs(j)
       write(iwr,2548)j,njds
2548   format(1x,'*** symmetry = ',i1,' no. of. SOMOS= ',i3)
       if(njds.gt.0) write(iwr,2547) (eigvals(loop,j),loop=1,njds)
       enddo
      endif
c
c     now GVB pair orbitals (if any)
c     njpair holds no of pair orbitals of each irrep
c     njps holds symmetry of GVB orbital
c     njp  nolds numbering of FBO orbital within irrep manifold
c
      if(npair.gt.0) then
       npair2 = npair + npair
       ij0 = ndomos+nsomos
        do i = 1, npair2
         isymm = iro(i+ij0)
         ik0 = 0
         do j = 1, isymm-1
         ik0 = ik0 + lj(j)
         enddo
         njps(i) = isymm
         njpair(isymm) = njpair(isymm) + 1
         njp(i) = ik0 + njd(isymm)+njs(isymm)+njpair(isymm)
        enddo
        if (debugl) then
         write(iwr,*) 'njpair = ', (njpair(loop),loop=1,8)
         write(iwr,*) 'njps = ', (njps(loop),loop=1,npair2)
         write(iwr,*) 'njp = ', (njp(loop),loop=1,npair2)
        endif
c
      endif
c
      write(iwr,102) 
102   format(/1x,'default reference functions'/
     +        1x,'===========================')
c
c
      if(npair.eq.0) then
       if (debugl) write(6,*) 'reference = ', (numb(i),i=1,nint)
       ilecj(1) = nsomos
       ij = 1
       do i = 1, ndomos + nsomos
        ij = ij + 1
        ilecj(ij) = numb(i)
       enddo
       write(11,101)ilecj(1),(ilecj(j),j=2,ij)
       write(iwr,103)ilecj(1),(ilecj(j),j=2,ij)
101    format (i4,',',/,256i4)
103    format(/1x,i2,2x,16(16i4/))
      endif
c
      if(npair.gt.0) then
c
c     must be careful with GVB pair calculations
c     as the nint MOs of numb include ALL the
c     pair orbitals. We must re-order the numb array.
c     Generate the (npair+npair) reference configurations
c     below
c
      ilecj(1) = nsomos
      do i = 1,nsomos
       ilecj(i+1) = numb(i)
      enddo
c
      nko = 1
c     first construct the leading term
      ik = 0
      ij0 = 0
      do i = 1, nrep
      ik0 = ik
       do j = 1, ndomos
        if(iro(j).eq.i) then
         ij0 = ij0 + 1
         ik0 = ik0 + 1
         numb(ij0) = ik0
        endif
       enddo
       k1 = 1
       do kk = 1, npair
        isymm = njps(k1)
        kk1 = k1
        if(isymm.eq.i) then
         ij0 = ij0 + 1
         numb(ij0) = njp(kk1)
        endif
        k1 = k1 + 2
       enddo
       ik = ik + lj(i)
      enddo
      do i = 1, ndomos+npair
      ilecj(nsomos+i+1) = numb(i)
      enddo
      ij = ndomos+nsomos+npair+1
      write(11,101)ilecj(1),(ilecj(j),j=2,ij)
      write(iwr,103)ilecj(1),(ilecj(j),j=2,ij)
      if(.not.ospecp) then
c     load up SCF config into cprop1e for default props
       call icopy(ij,ilecj,1,jkonp,1)
      endif
c    now add the double excitation from each pair to reference set
c

      do k = 1, npair
        nko = nko + 1
        ik = 0
        ij0 = 0
        do i = 1, nrep
        ik0 = ik
         do j = 1, ndomos
          if(iro(j).eq.i) then
           ij0 = ij0 + 1
           ik0 = ik0 + 1
           numb(ij0) = ik0
          endif
         enddo
         k1 = 1
         do kk = 1, npair
          isym1 = njps(k1)
          isym2 = njps(k1+1)
          if (k.eq.kk) then
c
c     double excitation of kth pair
c
            kk1 = k1 + 1
            isymm = isym2
          else
            kk1 = k1
            isymm = isym1
          endif
          if(isymm.eq.i) then
           ij0 = ij0 + 1
           numb(ij0) = njp(kk1)
          endif
          k1 = k1 + 2
         enddo
        ik = ik + lj(i)
        enddo
      do i = 1, ndomos+npair
      ilecj(nsomos+i+1) = numb(i)
      enddo
      ij = ndomos+nsomos+npair+1
      write(11,101)ilecj(1),(ilecj(j),j=2,ij)
      write(iwr,103)ilecj(1),(ilecj(j),j=2,ij)
c
      enddo
c
      return
c
      endif
c
c     load up SCF config into cprop1e for default props
c     analysis
c
      if(.not.ospecp) then
       call icopy(ij,ilecj,1,jkonp,1)
      endif
c
c     for closed and high-spin open shell calculations
c     add lowest single & double excitation from DOMO
c     manifold to low lying  VMOs (1 per symmetry)
c     nvmo will specify the no. of unoccupied MOS of
c     each irrep
c
c     for GVB-pair calculations, just add the two
c     configurations from each pair
c
      do i = 1, nrep
       nj(i) = njd(i) + njs(i)
       nvmo(i) = lj(i) - nj(i)
      enddo
c
c
c     treat each irrep in turn
c
c     1st - check that the symmetry of the input MOs, and hence
c     the generated reference configurations, is consistent
c     with the requested state symmetry of CI wavefunction

       ispacd = 1
       if (nsomos.gt.0) then
        ispacd = iro(ndomos+1)
        do i = ndomos+2, ndomos+nsomos
         ispacd = mult(ispacd,iro(i))
        enddo
       else
        ispacd = 1
       endif
       isymtop = iro(ndomos)
       isymtop1= iro(ndomos-1)

      if(ispacd.lt.1.or.ispacd.gt.8.or.ispacd.ne.ispacep) then
       write(6,3566)ispacd,ispacep
3566   format(//
     + 5x,'*********************************************************'/
     + 5x,'* input MO symmetry inconsistent with CI state symmetry *'/
     + 5x,'* symmetry of input orbital set = ', i8,'              *'/
     + 5x,'* requested CI state symmetry   = ', i8,'              *'/
     + 5x,'*********************************************************'/
     +  )
       call caserr(
     + 'input MO symmetry inconsistent with CI state symmetry')
      endif
c
      do k = 1, nrep
c
c     are there appropriate VMOs and a DOMO?
c
       if(nvmo(k).ne.0.and.njd(k).ne.0) then
c
       nterm1 = 1
       if(nvmo(k).lt.2) then
        nterm2 = 1
        nterm30 = 1
       else
        nterm2 = 2
c       nterm30 = 2
        nterm30 = max(2,nterm30d)
        if (nterm30.gt.2) write(iwr,*) '# extra refs set to',nterm30,
     +                  ' due to degeneracies // check orbital order'
        if (ifbuen) then
c       this is set artificially high to avoid crucial singles
c       being omitted from the reference set
c        if(nvmo(k).ge.6) nterm30 = 6
         if(nsomos.eq.0.and.nvmo(k).ge.3) nterm2 = 3
        endif
       endif
       
c
       do nnvmo = 1, nterm1
        nko = nko + 1
        ik = 0
        ij0 = nsomos
        do i = 1, nrep
         ik0 = ik
         do j = 1, ndomos
          if(iro(j).eq.i) then
           ij0 = ij0 + 1
           ik0 = ik0 + 1
           numb(ij0) = ik0
          endif
         enddo
         if(ik0.ne.ik.and.i.eq.k)
     +      numb(ij0) = ik0 + njs(i) + nnvmo
         ik = ik + lj(i)
        enddo
c
c      add the root somos
c
        ik = 0
        ij0 = 0
        do i = 1, nrep
         ik0 = ik
         do j = ndomos+1,ndomos+nsomos
          if(iro(j).eq.i) then
           ij0 = ij0 + 1
           ik0 = ik0 + 1
           ik1 = ik0 + njd(i)
           numb(ij0) = ik1
          endif
         enddo
         ik = ik + lj(i)
        enddo

        if (debugl) write(iwr,*) 'reference = ',
     +             (numb(i),i=1,nint)
c
        ilecj(1) = nsomos
        ij = 1
        do i = 1, ndomos+nsomos
         ij = ij + 1
         ilecj(ij) = numb(i)
        enddo
        write(11,101)ilecj(1),(ilecj(j),j=2,ij)
        write(iwr,103)ilecj(1),(ilecj(j),j=2,ij)
       enddo
c
c      and the corresponding single excitations from the DOMO
c      to VMO manifold
c      Increase default for top DOMO symmetry
c
       if(ifbuen.and.nsomos.eq.0.and.
     +    (k.eq.isymtop.or.k.eq.isymtop1).and.
     +    nvmo(k).ge.5) then
        nterm2u = 5
       else
        nterm2u = nterm2
       endif
       do nnvmo = 1, nterm2u
        nko = nko + 1
        ik = 0
        ij0 = nsomos
        do i = 1, nrep
         ik0 = ik
         do j = 1, ndomos
          if(iro(j).eq.i) then
           ij0 = ij0 + 1
           ik0 = ik0 + 1
           numb(ij0) = ik0
          endif
         enddo
         if(ik0.ne.ik.and.i.eq.k) then
         nsing1 = numb(ij0)
         nsing2 = ik0 + nnvmo + njs(i)
         ij0 = ij0 -1
         endif
         ik = ik + lj(i)
        enddo
c       add the root SOMOS
        ik = 0
        ij0 = 0
        do i = 1, nrep
        ik0 = ik
         do j = ndomos+1,ndomos+nsomos
         if(iro(j).eq.i) then
          ij0 = ij0 + 1
          ik0 = ik0 + 1
          ik1 = ik0 + njd(i)
          numb(ij0) = ik1
         endif
         enddo
        ik = ik + lj(i)
        enddo
c
         if(debugl) write(iwr,*) 'reference = ',
     +          nsing1, nsing2, (numb(i),i=1,nint-1)
         do i = 1, ndomos + nsomos - 1
          ilecj(i) = numb(i)
         enddo
c
         m2 = 2 + nsomos
c
c      circumvent (pathetic?) ordering imposed by MRDCI on
c      somo ordering
c
         nsorto(1) = nsing1
         nsorto(2) = nsing2
         do i = 1, nsomos
          nsorto(i+2) = ilecj(i)
         enddo
         do min = 1 , m2
            jm = min
            im = nsorto(min)
            do  j = min , m2
               if (nsorto(j).lt.im) then
                  im = nsorto(j)
                  jm = j
               end if
            enddo
            if (jm.ne.min) then
               itmp = nsorto(jm)
               nsorto(jm) = nsorto(min)
               nsorto(min) = itmp
            end if
         enddo
c
         write(11,101)m2,(nsorto(i),i=1,m2),
     +                 (ilecj(j),j=nsomos+1,ndomos+nsomos-1)
         write(iwr,103)m2,(nsorto(i),i=1,m2),
     +              (ilecj(j),j=nsomos+1,ndomos+nsomos-1)
c
        enddo
c
c      and the corresponding single excitations from the SOMOs
c      to VMO manifold. The number to be added (nterm3) depends on 
c      the nature of the SOMO (lower valence or excited)
c
       if(njs(k).ne.0) then
c
c      decide on nterm3
c
c      emin = 0.0d0
c      do i = 1,nrep
c      if (i.ne.k) then
c       do j = 1,njs(i)
c        emin = min(emin,eigvals(j,i)
c       enddo
c      enddo
c      if(abs(eigvals(njs(k))-emin).lt.0.2d0) then
c
       nterm3 = nterm30
       if(ifbuen.and.abs(eigvals(njs(k),k)).lt.0.3d0.and.
     +    nvmo(k).ge.6) then
c       excited
        nterm3 = 6
       endif
       if(debugl) write(6,*) 'symmetry, nterm3, eigvals = ', 
     +            k,nterm3, eigvals(njs(k),k)
c
        do nnvmo = 1, nterm3
         nko = nko + 1
         ik = 0
         ij0 = nsomos
         do i = 1, nrep
          ik0 = ik
          do j = 1, ndomos
           if(iro(j).eq.i) then
            ij0 = ij0 + 1
            ik0 = ik0 + 1
            numb(ij0) = ik0
           endif
          enddo
          ik = ik + lj(i)
         enddo
c        add the root SOMOS plus the single excitations from the 
c        top SOMO
         ik = 0
         ij0 = 0
         do i = 1, nrep
         ik0 = ik
          do j = ndomos+1,ndomos+nsomos
          if(iro(j).eq.i) then
           ij0 = ij0 + 1
           ik0 = ik0 + 1
           ik1 = ik0 + njd(i)
           numb(ij0) = ik1
          endif
          enddo
          if(i.eq.k) then
           nsing1 = numb(ij0)
           numb(ij0) = ik1 + nnvmo
           nsing2 = numb(ij0)
          endif
         ik = ik + lj(i)
         enddo
c
         if(debugl) write(iwr,*) 'reference = ',
     +           nsing1, nsing2, (numb(i),i=1,nint-1)
         do i = 1, ndomos + nsomos
          ilecj(i) = numb(i)
         enddo
c
c      impose ordering expected by MRDCI on somo ordering
c
         do i = 1, nsomos
          nsorto(i) = ilecj(i)
         enddo
         do min = 1 , nsomos
           jm = min
           im = nsorto(min)
           do  j = min , nsomos
              if (nsorto(j).lt.im) then
                 im = nsorto(j)
                 jm = j
              end if
           enddo
           if (jm.ne.min) then
              itmp = nsorto(jm)
              nsorto(jm) = nsorto(min)
              nsorto(min) = itmp
           end if
         enddo
c
        write(11,101)nsomos,(nsorto(i),i=1,nsomos),
     +                 (ilecj(j),j=nsomos+1,ndomos+nsomos)
        write(iwr,103)nsomos,(nsorto(i),i=1,nsomos),
     +              (ilecj(j),j=nsomos+1,ndomos+nsomos)
c
        enddo

       endif
*****
c
       endif
c
      enddo
c
      endif
      return
      end
      subroutine chk_mrdci_data(nbox,nsym)
c
c     This routine scans the input data for the semi-direct 
c     MRDCI module, particularly the input configuration list.
c     This was originally dealt with at selection time i.e. too late.
c     This routine is designed to remove redundant configurations
c     from the input configuration list read from nf11,
c     rather than flag an error condition.
c     Mains with up to 9 open shells are accepted.
c     configurations up to seventh sk are generated.
c     Resulting configurations having more then 9 open shells will not
c     be taken into account any further.
c
      implicit real*8 (a-h,o-z)
c
      integer nbox,nsym
c
c                                                                        c      
c------------------------------------------------------------------------c      
c  include - file fuer parkwa2.f                                         c      
c------------------------------------------------------------------------c      
c                                                                               
c --- festlegung der deks                                                       
c --- festlegung der parameterinhalte                                           
      integer nteint
*  laenge der table-information (meist auf /a/ abgelegt)                        
      integer ndtab
      parameter (ndtab = 7 829)
*  maximale dimension der ci
      integer maxci
      parameter (maxci= 800 001)
      integer nlca
*  laenge des commons a : nlca
      parameter (nlca =  39 218)
      integer nt44
      parameter (nt44 =  44 000)
*  laenge verschiedener vektoren : noch unbekannt !!!                           
      integer n3zt
      parameter (n3zt   = 30 000)                                               
*  number of mains                                                             
      integer maxref
      parameter (maxref = 256)                                                  
*  number of roots
      parameter (mxroot = 50)
* maxshl haengt mit der anzahl der elektronen zusammen !                        
* maxele maximum number of electrons
      integer maxele, maxshl, ndcon
      parameter (maxele = 90 )                                                  
      parameter (maxshl = 50 )                                                  
      parameter (ndcon  = maxshl*maxref )                                       
*  anzahl der offenen schalen                                                   
      integer nopmax
      parameter (nopmax = 9)                                                    
      integer n12, n36
      parameter (n12 = 12 )                                                     
      parameter (n36 = 36 )                                                     
*  anzahl der superkategorien                                                   
      integer ndk5, iswhm
      parameter (ndk5 = 5 )                                                     
      parameter (iswhm= ndk5)                                                   
      integer nn3
      parameter (nn3  = 3 )                                                     
*  d2h - punktgruppe                                                            
      integer maxsym
      parameter (maxsym = 8)                                                    
*  laenge des iot - feldes                                                      
      integer niot
      parameter (niot = 100 000)                                                
*  null und andere numerische parameter                                         
      integer nzero
      parameter (nzero = 0 )                                                    
*  laenge des labelfiles                                                        
      integer nnid, nnid8
      parameter (nnid   = 2 000)                                                
      parameter (nnid8  =   250)                                                
*  laenge der integralbloecke vom ft31 : stoney                                 
      integer ndimh
      parameter (ndimh  = 1 000)                                                
*  laenge des nit-feldes                                                        
      integer nitmax
      parameter (nitmax =  666)                                                 
c------------------------------------------------------------------------c      
c ---- end of  include file                                              c      
c------------------------------------------------------------------------c      
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
       integer nhead
       common /ftap/ nston, mtape, mdisk,
     .               ideli, ltype, linf,  ntab,  kfile,
     .               kclab, ntype, mstvt, nf01, nf62,
     +               nhead, nf99, mtapev, nf11
c
      integer nytl
      common /cny/  nytl(ndk5)
c
      integer nconf
      integer jab, isym, mj, iab, kj, jcon, icuf, iduf
      integer kcon, ilee, ntil, nbal, jkon, nop
      integer isw, nyzl, nj, ijoe, ndub
      parameter (nmmo=256)
      common /aeng/ nconf(ndk5),jab(n36),isym(maxsym),
     .              mj(maxsym),iab(nmmo),kj(maxsym),
     .              jcon(nmmo),icuf(n12),iduf(maxshl),
     .              kcon(ndimh),ilee(maxsym+1),
     .              ntil(maxsym),nbal(maxsym+1),
     .              jkon(ndcon),nop(maxref),
     .              isw(ndk5),nyzl(maxref),nj(maxsym),
     .              ijoe(maxref),ndub(maxref)
c
      logical oconf,debugs
      common/parkin/egey1,trash,tdel,
     +       ical,nele,nko,mxex,nmul,ispace,nprin,
     +       nft31,nstarv,ncorci,nform,
     +       oconf,debugs,
     +       lsng,maxciin,ipt0,isymin(mxroot)
c
      integer jkonr(maxshl)
c
      integer isaf
      dimension isaf(0:9)
c
      parameter (ndkneu=ndk5+2,nopneu=nopmax+4)
c
      integer nytln,nconfn,iswn
      common /momain/ nytln(ndkneu),nconfn(ndkneu),
     +      iswn(ndkneu)
c
      nops   = 5
      nbxa=nbox-1
      nbxb=nbox-2
      nnnx=(nele+nmul-3)/2
      do i=1,ndkneu
       nnnx=nnnx+1
       nytln(i)=nnnx
      enddo
      do kdoo=1,ndk5
        nytl(kdoo)=nytln(kdoo)
      enddo
      write(iwr,8497) nbox,nele,nko,mxex,nsym,nmul,ispace,nprin
      if(mxex.gt.4.or.mxex.lt.0) then
       write (iwr,768)
       call caserr('excitation class not allowed')
      endif
      if(nmul.lt.1) then
       write(iwr,784)
       call caserr('multiplicity out of bounds')
      endif
      do i=1,ndk5
        isw(i)=0
      enddo
      do i=1,ndkneu
       iswn(i)=0
      enddo
      if (nmul.eq.1) then
        icount=1
        do kdoo=2,ndkneu,2
          icount=icount+1
          iswn(kdoo)=icount
        enddo
       else
        icount=0
        do kdoo=nmul-1,ndkneu,2
          icount=icount+1
          iswn(kdoo)=icount
        enddo
      endif
      if(nele.lt.2.or.nele.gt.maxele) then
       write(iwr,780)
       call caserr('requested number of electrons troublesome')
      endif
      if(ispace.lt.1.or.ispace.gt.8) then
       write(iwr,786)
      call caserr('symmetry parameter out of bounds')
      endif
      if(nko.gt.maxref) then
       write(iwr,760) nko,maxref
       call caserr('too many reference configurations')
      endif
c
      im=0
      call rewftn(nf11)
c --- read the configurations (mains) from input
c --- stream nf11
      call rewftn(nf11)
c     process list of mains, removing duplicate entries
c
      nkorig = nko
      nko = 0
10    nt = im +    1
      im = im + maxshl
c ---- reading the configurations
      do j=nt,im
       jkon(j) = 0
      enddo
      read (nf11,*,end=20,err=20)np
      mm = (nele+np)/2
cvp reference configurations with up to nopmax open shells are allowed
      if(np.gt.nopmax) then
       write(iwr,776)
       call caserr('too many open shells in mains for dimensions')
      endif
      read (nf11,*)(jkonr(j),j=1,mm)
      nko = nko + 1
      do j=nt,im
       jkon(j) = jkonr(j-nt+1)
      enddo
      nyt=1
      if(np.eq.0) go to  94
      nyt = iswn(np)
 94   nyt=nytln(nyt)
      nyzl(nko)=nyt
      if(nko.eq.1) go to 503
      i1=nko-1
      mx=-maxshl
      do 504 j=1,i1
      mx=mx+maxshl
      if(nop(j).ne.np) go to 504
      nv=mx
      la=nt-1
      do k=1,nyt
       la=la+1
       nv=nv+1
       if(jkon(la).ne.jkon(nv)) go to 504
      enddo
      write(iwr,778)
      write(iwr,777) nop(i), (jkon(nt-1+loop),loop=1,nyt)
      nko = nko - 1
      do loop = nt, im
       jkon(loop) = 0
      enddo
      im = im - maxshl
      go to 10
 504  continue
c
 503  na=nmul+np
c
c     I dont understand the "np.gt.nmul+3" condition below given the 
c     current setting for nopmax e.g. 6 open shells when nmul=1 
c     will trigger it.
c     I've tried editing it out, but causes problems with 
c     "not all mains generated"
c
      if(na-2*(na/2).eq.0 .or. np.gt.nele .or.np.gt.nmul+3) then
       write(iwr,772)
       call caserr(
     + 'open shell structure inconsistent with multiplicity')
      endif
      nop(nko)=np
      ndub(nko)=(nele-np)/2
      if(np.lt.2) go to 500
      nz=jkon(nt)
      if(nz.lt.1.or.nz.gt.nbox) go to 773
      do j=2,np
       nt=nt+1
       nv=jkon(nt)
       if(nv.le.nz) go to 769
       if(nv.gt.nbox) go to 773
       nz=nv
      enddo
      if(np.eq.nele) go to 10
518   mx=np+im-maxshl+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      lb=im-maxshl
      do 506 j=1,np
      lb=lb+1
      if(jkon(lb)-nv) 506,769,507
506   continue
507   jm=np+2
      if(jm.gt.nyt) go to 10
      kp=mx
      do 508 j=jm,nyt
      kg=mx
      mx=mx+1
      nv=jkon(mx)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      do 509 k=kp,kg
      nz=jkon(k)
      if(nv.le.nz) go to 769
 509   continue
      kg=im-maxshl
      do 510 k=1,np
      kg=kg+1
      if(jkon(kg)-nv) 510,769,508
 510  continue
 508  continue
      go to 10
 500  if(np.eq.0) go to 511
      nv=jkon(nt)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      go to 518
 511  nz=jkon(nt)
      if(nz.lt.1.or.nz.gt.nbox) go to 773
      if(nyt.eq.1) go to 10
      kp=nt
      do 512 j=2,nyt
      kg=nt
      nt=nt+1
      nv=jkon(nt)
      if(nv.lt.1.or.nv.gt.nbox) go to 773
      do 512 k=kp,kg
      nz=jkon(k)
      if(nv.le.nz) go to 769
 512  continue
      go to 10
  20  continue
c
      write(iwr,8498) nko
c
c     check requested size of zero order space
c
      do  j=(nmul-1),nopmax,2
       isaf(j)=numsaf(nmul,j)
      enddo
      k=0
      do j=1,nko
        k= k + isaf(nop(j))
      enddo
      write(iwr,553) nko,k
      if (k.gt.maxref)
     + call caserr('too many SAFs from specified mains')

c     end of reading and testing the reference configurations

      write(iwr,520)
      kg = -maxshl
      do i=1,nko
         kg=kg+maxshl
         np=nop(i)
         kb=kg+1
         nv=kg+nyzl(i)
         write(iwr,522)i, np,(jkon(j),j=kb,nv)
      enddo
c
c     revise nf11 id identical mains detected
      if(nkorig.ne.nko) then
       call rewftn(nf11)
       im = 0
       do i=1,nko
        nt = im +    1
        im = im + maxshl
        mm = (nele+nop(i))/2
        write(nf11,101) nop(i), (jkon(j),j=nt,nt+mm-1)
       enddo
      endif
      call rewftn(nf11)
c
      return
c
c     error section
c   
 769  write(iwr,770) i,nv,nz,jkon(kg),kg,jkon(lb),lb
 770  format(5x,'pauli was right or maybe permutation error',7i5)
      call caserr('possible permutation error')
 773  write(iwr,774)
 774  format(5x,'orbital numbering is weird in mains')
      call caserr(
     + 'invalid orbital index in reference function')
      return
 101  format (i4,',',/,256i4)
 522  format(5x,i3,8x,i3,9x,24i4)
 520  format(/10x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 10x,
     +'numbers of open shells and corresponding main configurations'/
     + 10x,
     +'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + / 5x, 'input     number of       main'
     + / 5x, 'index     open shells     configuration'/)
 553  format(/1x,'The set of',i4,' main configurations leads to',
     +i4,' spin adapted functions.')
8498  format(/1x,'no. of processed mains ',i3)
8497  format(/1x,'number of orbitals ',i7/
     &        1x,'no. of electrons   ',i7/
     &        1x,'no. of input mains ',i7/
     &        1x,'excitations        ',i7/
     &        1x,'no. of irreps.     ',i7/
     &        1x,'multiplicity       ',i7/
     &        1x,'ispace             ',i7/
     &        1x,'nprint             ',i7)
 768  format(5x, 'excitation class not allowed')
 784  format(5x,'multiplicity out of bounds')
 780  format(5x,'requested number of electrons troublesome')
 786  format(5x,'ispace parameter out of bounds')
 760  format(5x,'too many mains',2i6)
 776  format(5x,'too many open shells in mains for dimensions')
 778  format(/10x,
     +  '**             WARNING                        **'/10x,
     +  '* two of the main configurations are identical *'/10x,
     +  '* duplicate removed and number of mains reset  *'/)
 777  format(10x,'**',i5,2x,50i3)
 772  format(5x,
     + 'open shell structure in mains inconsistent with multiplicity')
      end
      subroutine ver_newmrd1(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd1.m,v $
     +     "/
      data revision /"$Revision: 6250 $"/
      data date /"$Date: 2012-01-23 16:07:50 +0100 (Mon, 23 Jan 2012) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
