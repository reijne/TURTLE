      subroutine adap2e(core)
c
c ----- generate full 2-e integral file and symmetry adapt
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      character*8 fnm
      character*6 snm
      data fnm,snm/"master.m","adap2e"/
c
      dimension core(*)
      data ono / .false./
c
      inull = igmem_null()
      dum = 0.0d0
      idum = 0
      ntsave = nt
c
c ----- reset symmetry to c1
c
      nt = 1
      iofsym = 1
      if (.not.(opass1)) then
         if (irest.lt.1) call standv(0,core)
c
c ----- 1e-integrals
c
         if (tim.ge.timlim) go to 20
      end if
c
c ----- 2e-integrals
c
      if (opass2) then
         if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         nopk = 1
         if (irest.le.1) then
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "prefac",IGMEM_DEBUG)
            call rdmake(core(iprefa))
            call jandk(zscftp,core,core(inull),core(inull),core(inull),
     +                 core(inull),core(inull),core(iprefa),core(inull))
            call gmem_free_inf(iprefa,fnm,snm,"prefac")
         endif
         if (tim.ge.timlim) go to 20
      end if
      if (irest.le.2) call adapti(core,ono)
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 20   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
      return
      end
      subroutine analis(core)
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
      integer idpa
      common/dppopc/idpa

      common/dfcalc/cdata(5,mxcalc),ictype(mxcalc),icsect(mxcalc),
     &     icgrid(mxcalc),icgsec(mxcalc),ndata(mxcalc),
     &     icstat(mxcalc),icdata(5,mxcalc),ncalc

      common/dfgrid/geom(12,mxgrid),gdata(5,mxgrid),igtype(mxgrid),
     &     igsect(mxgrid),nptot(mxgrid),npt(3,mxgrid),
     &     igdata(5,mxgrid),igstat(mxgrid),ngdata(mxgrid),
     &     ngrid
      character*10 charwall
      dimension core(*)
c
      call cpuwal(begin,ebegin)
      write (iwr,6010) begin,charwall()
      nav = lenwrd()
      newbas = num
      nbfnd = num
      ied3 = idaf
      nwl = 2*(maxorb+1)/nav
      nwm = (maxorb+6)/nav 
      nw1e2 = 400 + 204/nav
c     nwp2 = 2*maxorb + 82 + 60/nav
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      if(iplot.eq.0)ibl7la = 1
      if (idma.gt.0) ibl7la = ibl7la + lensec(maxat) + lensec(nwdma)
      if (iprop.gt.0) ibl7la = ibl7la + lensec(nw1e2)+ lensec(831)
      if (local.gt.0) ibl7la = ibl7la + lensec(nwl)
      if (imull.gt.0) ibl7la = ibl7la + lensec(nwm) + lensec(num)
c
      ofcheck = (iprop.gt.0)

      if (iplot.gt.0) then
c
c F-function limitation applies only to 
c mo, density and their gradients and 
c atomic density
c
         do icalc = 1,ncalc
            if(ictype(icalc) .eq. 2 .or.
     &         ictype(icalc) .eq. 3 .or.
     &         ictype(icalc) .eq. 5 .or.
     &         ictype(icalc) .eq. 6) then
               ofcheck = .true.
            endif
         enddo

      end if
c
c need to rationalise this with later f/g checks
c
      if (ofcheck .and. ispchk().ge.4) then
         write (iwr,6020)
         go to 30
      end if
c
c     ibl7la is 1st vacant block for use in =analysis= mode
c
      if (local.ge.1) then
c
c ----- boys or pipek/mezey localisation
c
         call boyloc(core)
         if (tim.ge.timlim) go to 30
      end if
c
c ----- allocate core storage (all available memory except in 
c ----- nbo analysis where core may be need for direct-SCF)
c
      if (inbo.ge.1) then
c
c -----  natural bond order analysis
c
         call runnbo(core)
         go to 30
      else
c
        ibase = igmem_alloc_all(lwordp)
c
        if (iplot.ge.1) then
c
c ----- graphical analysis (new code)
c
           call graph(core(ibase),ofcheck)
           if (tim.ge.timlim) go to 20
        end if
c
c  ---- potential fitted atomic charges
c
        if(opotf(0))call pdc1(core(ibase))
c
c ----- distr mult analysis
c
        if (idma.ge.1) then
           call dma(core(ibase))
           if (tim.ge.timlim) go to 20
        end if
        if (iprop.ge.1) then
c
c ----- 1e-properties
c
         call propg
         call propin(core(ibase))
         if (tim.ge.timlim) go to 20
        end if
        if (imull.ge.1) then
c
c ----- mulliken analysis
c
         call analm(core(ibase),lwordp)
        end if
      endif
c
c ----- reset core allocation
c
20    call gmem_free(ibase)
      lwordp=lwordpold
        if (idpa.ge.1) then
            call dppop(core)
        endif
c
30    cpu = cpulft(1)
      write (iwr,6030) cpu,charwall()
      call timana(14)
c
      return
 6010 format (/1x,104('-')//40x,21('*')/40x,'wavefunction analysis'/40x,
     +        21('*')//' commence analysis at ',f12.2,' seconds',a10,
     +        ' wall'/)
 6020 format (/1x,
     + '**** analysis option precluded by f- or g-functions in basis')
 6030 format (//' end of wavefunction analysis at ',f12.2,' seconds',
     +        a10,' wall',//1x,104('-')/)
      end
      subroutine hfscf(core)
c
c     ----- calculate hf energy of molecule -----
c                  rhf , uhf or casscf method
c     1. compute 1e-integrals
c     2. compute 2e-integrals
c     3. carry out scf procedure
c
c     note- initial density matrix is always assumed.
c
c     irest = 0     normal start and normal running condition.
c     irest = 1     2e-integral restart ( 1e + mo*s saved)
c     irest = 2     scf restart ( 1e + mo*s saved; 2e saved)
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      character *8 title,guess,conf
      common/restrz/title(12),guess,conf(37)
      common/restrl/ociopt(2),omp2,ogr(12),
     +omcscf,oforce,oci,ocart(32)
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      real*8 ptspac, delect, facnrm
      integer iconv,inkblk,itotbq,iscf,ispace,msecp
      logical opssc
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
      integer ncore,nact,nels
      common /fccwfn2 / ncore,nact,nels
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)

cdrf
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
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
      integer idafh,navh
      common/hdafile/idafh,navh,ioda(2,1000)
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
      character*8 fieldold
      logical oarf, oarfcvg
      common/enrghf/enhh,etothh,ehfhh
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
c
      real*8 timcheck
      character*8 fnm
      character*5 snm
      data fnm,snm/"master.m","hfscf"/
c
      dimension core(*),zcas(3)
      data szero /0.0d0/
      data ono/.false./
      data zuhf,zcas/'uhf','casscf','mcscf','vb'/
c     data zgrhf,zgvb,/'grhf','gvb'/
      data m16,m10/16,10/
c
c     ----- reset symmetry to c1 if this is a casscf run
c           or a gvb run with npair .ne. 0 ,or user requested
c
      inull = igmem_null()
      ntsave = nt
      if (npair.ne.0 .or. zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2) .or.
     +    zscftp.eq.zcas(3) .or. offsym) then
       nt = 1
       iofsym = 1
      endif
      enrgy = szero
      if (oprint(48)) then
         i10 = igmem_alloc_inf(nw196(5),fnm,snm,"i10",IGMEM_DEDUG)
         itmp = lenint(nat*nt)
         i11 = igmem_alloc_inf(itmp,fnm,snm,"i11",IGMEM_DEDUG)
         call symtab(core(i10),nshell,core(i11),nat)
         call gmem_free_inf(i11,fnm,snm,"i11")
         call gmem_free_inf(i10,fnm,snm,"i10")
      endif
c
      if (.not.(opass1)) then
c
c     ----- 1e- integrals -----
c
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
c
cjmht Need to be using the same time here, as otherwise different
c     processors end up going through different paths in the code.
c     Copy tim into timcheck so that we don't overwrite everyone
c     else's time with roots' as this could screw timing analysis.
c     Bit dirty, but hard-coding type as 7123 and length as 8 as 
c     that's what everyone else seems to have done. Baaaaaa...
c
            timcheck=tim
            if (timcheck.ge.timlim) go to 20
         endif
      end if
      ionsv = 1
c
c     ----- 2e-integrals -----
c
      if (omp2 .or. mp3 .or. oci .or. omcscf .or. omas) nopk = 1
      if (.not.(opass2 .or. odscf)) then
c        open = zscftp.eq.zuhf .or. zscftp.eq.zgvb .or.
c    +        zscftp.eq.zgrhf .or. na.ne.nb
c        opandk = nopk.ne.1 .and. open
         call start_time_period(TP_INTEG)
         if (irest.le.1) then
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "iprefa",IGMEM_NORMAL)
            call rdmake(core(iprefa))
            call jandk(zscftp,core,core(inull),core(inull),
     +                        core(inull),core(inull),core(inull),
     +                        core(iprefa),core(inull))
            call gmem_free_inf(iprefa,fnm,snm,"iprefa")
         endif
         call end_time_period(TP_INTEG)
         if (tim.ge.timlim) go to 20
      else
       if (opass2) then
          if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
       end if
      end if
c
      if (zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
     &    .or.zscftp.eq.zcas(3)) then
      oarf = (field(5:) .eq. 'scf ' .or. field(:4) .eq. 'scf ')
      else
      oarf = .false.
      endif
      fieldold = field
      enadd = 0.d00
      if (oarf) ocifp = .true.
      narfpts = 0
99999 continue
      oarfcvg = .false.
      if (oarf) then
            call arfupd(core,enadd,fieldold)
            field = ' '
         write(iwr,901) narfpts+1
901   format(/,1x,78('='),/,10x,'Average Reaction field MCSCF/CI ',
     1       'Calculation cycle',i3,/,1x,78('='),/,1x)
      end if
c
c     ----- symmetry adaption
c
      if (zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)) then
         if (.not.opass4) then
            if (irest.le.2) call adapti(core,ono)
            if (tim.ge.timlim) go to 20
         end if
         nt = ntsave
      end if
c
c     ----- scf -----
c
      if (.not.opass5) then
         if (irest .le. 4) then
            if (opssc) then
               call scf2(core)
            else if(omorok) then
                 if(odscf) then
                  call morokdscf(core)
                 else
                  call morokscf(core)
                 endif
            else if(zscftp.eq.zcas(3)) then
                 call vbstart(core,2)
c                goto 1000
            else
               call scf(core)
               if ((odscf. and. guess.eq.'atoms').or.oso) go to 1000
            end if

            call putdev(core,mouta,7,1)
            if (zscftp.eq.zuhf) call putdev(core,moutb,10,2)
c
            call secget(isect(494),m16,iblk34)
            call rdedx(potnuc,m10,iblk34,idaf)
            enuc = potnuc(1)
            ehfock = potnuc(2)
            etotal = potnuc(3)
c
            call secget(isect(13),13,iblok)
            call wrt3(enuc,lds(isect(13)),iblok,idaf)
         end if
      end if
c
cahv see comments in scf.m  0206
      if (fieldold .ne. ' ') then
        enrgy = etotal
c..     restore/generate density matrix
        call dipmat(zscftp,core)   
        if (zscftp.eq.zcas(1).or.zscftp.eq.zcas(2)
     *      .or.zscftp.eq.zcas(3)) then
           l2 = num*(num+1)/2
           i10 = igmem_alloc_inf(l2,fnm,snm,"i10-l2",IGMEM_NORMAL)
           call rdedx(core(i10),l2,ibl7pa,num8)
           call dawrit(idafh,ioda,core(i10),l2,16,navh)
           call gmem_free_inf(i10,fnm,snm,"i10-l2")
        endif
      endif
      if (oarf) then
c
         call arfcvg(core,enrgy,enadd,oarfcvg,narfpts)
c
      end if
      if (fieldold .ne. ' ') then
         etothh = enrgy
      endif
      narfpts = narfpts + 1
      if (oarf .and. .not. oarfcvg) goto 99999
      if (oarf) field = fieldold
c
1000  opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
 20   nt = ntsave
c
      return
      end
      subroutine initiala
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
c symmetry assignment
      common/atmol3o/osim1,osim2
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
c
c
c****** direct-scf common blocks
c
      logical oentry
      integer idblk, numdir, inddir
      common /wrtd/ idblk(30),numdir,inddir,oentry
c
c****** gen2 contains the integral test data
      common/gen2/ripftd(21,2),rijkld(21,2),rjpftd(21,2),
     + ripfal(21,2),rijkal(21,2),
     + prefac(21),tinner(21),rifall,riffew,ricfb1,ricfb2
c
c****** gen contains most of the stuff about the options
      common/gen/maxitd,iterd,concrt,conv,energy,thresh,dcmax,
     -      olog1,olog2,olog5,olog9,olog6,odbg(10),ogmp2,mp2fo,mp2fv,
     -      oread(3),iaccur,ocvary,ogrdt(5),
     -      noptd,npoint,timeg,ogmull,iswap,dshift,irhess
      logical ofull,odiisn
      common/gen3/ofull(100),odiisn
c
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
c
      real*8  dx0, dxm, dele, drms
      integer ntpr, ntrun, imcyc, ncyc, ntmin, nspacb
      common /modmin/ ntpr,ntrun,imcyc,ncyc,ntmin,nspacb,
     +                dx0,dxm,dele,drms
c
c
      real*8 bbeta, box
      integer ntx, ntxo, nrc, nrcx, npm, nrp, nsm, nram, ntb
      integer nspaca
      common /moda/ ntx,ntxo,nrc,nrcx,npm,nrp,nsm,nram,ntb,nspaca,
     +  bbeta,box(3)
c
c
      integer ntp, ntc, ntcc, ntf, ntid, ntn, ntnb, nsnb
      integer nrep, nrepc, ipro
      common /modf/ ntp,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,nrep,
     +              nrepc,ipro
c
c
      real*8 tempi, tempo, tautp, t, tims
      integer nstlim, mdpr, init, nspacc
      common /modmd/ nstlim,mdpr,init,nspacc,
     +               tempi,tempo,tautp,t,tims
c
c
      integer natom, nres, nbonh, nbona, ntheth, ntheta, nphih
      integer nphia, nnb, ntypes, nrcp, nconp, mema, memb, maxnb
      integer mbona, mtheta, mphia, natc, npro, numbnd, numang
      integer numphi, natyp, nphb, imodf, jmodf
      common /modfov/ natom,nres,nbonh,nbona,ntheth,ntheta,nphih,
     +                nphia,nnb,ntypes,nrcp,nconp,mema,memb,maxnb,
     +                mbona,mtheta,mphia,natc,npro,numbnd,numang,
     +                numphi,natyp,nphb,imodf(25),jmodf(50)
c
c
      integer nword, isw, kprint, ifzmat
      common /mod2d/ nword,isw,kprint,ifzmat
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
      common/nbterm/cut,scnb,scee,dielc,idiel,nbuck
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
c
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      common /zjorgs/ystat(3)
      common /jorgs/maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     *outjor,outdeb,maxsta
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
      real*8 tr, trx, try, trz, prmoms, aprev
      integer indmx, jaxis, igroup, igrp80
      logical oprmom
      common /molsym/ tr(3,3),trx,try,trz,indmx,jaxis,igroup,igrp80,
     +                prmoms(3),aprev(maxat3,3),oprmom
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      real*8 slder, alimin
      common/trntim/slder,alimin
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
c The oblock array determines what is written out to the punchfile:
c
c Below sourced from server.m
c oblock(1)   - coords - single point or last point
c keyword: coor
      integer PUN_COORD
      parameter (PUN_COORD=1)
c oblock(2)   - coordinates at each optimisation step
c keyword: opti
      integer PUN_OPT_COORD
      parameter (PUN_OPT_COORD=2)
c oblock(3)   - starting coordinates
c keyword: init
      integer PUN_START_COORD
      parameter (PUN_START_COORD=3)
c oblock(4)   - atom connectivity
c keyword: conn
      integer PUN_ATOM_CONN
      parameter (PUN_ATOM_CONN=4)
c oblock(5)   - molecular orbital occupations
c keyword: occu
      integer PUN_MO_OCC
      parameter (PUN_MO_OCC=5)
c oblock(6)   - eigenvectors
c keyword: vect
      integer PUN_EIGV
      parameter (PUN_EIGV=6)
c oblock(7)   - run type
c keyword: type
      integer PUN_RUNT
      parameter (PUN_RUNT=7)
c oblock(8)   - gvb pair data
c keyword: gvb
      integer PUN_GVBPAIR
      parameter (PUN_GVBPAIR=8)
c oblock(9)   - force constant matrix
c keyword: forc
      integer PUN_FCMAT
      parameter (PUN_FCMAT=9)
c oblock(10)  - vibrational frequencies
c keyword: vibr
      integer PUN_VIBFREQ
      parameter (PUN_VIBFREQ=10)
c oblock(11)  - basis set
c keyword: basi
      integer PUN_BASIS
      parameter (PUN_BASIS=11)
c oblock(12)  - orbital energies
c keyword: eige
      integer PUN_ORB_ENER
      parameter (PUN_ORB_ENER=12)
c oblock(13)  - total scf energies
c keyword: scfe
      integer PUN_SCF_ENER
      parameter (PUN_SCF_ENER=13)
c oblock(14)  - job title
c keyword: titl
      integer PUN_JOB_TITLE
      parameter (PUN_JOB_TITLE=14)
c oblock(15)  - overlap matrix
c keyword: over
      integer PUN_OVLP_MAT
      parameter (PUN_OVLP_MAT=15)
c oblock(16)  - grid data from dumpfile section(s)...
c keyword: grid
      integer PUN_GRID_DATA
      parameter (PUN_GRID_DATA=16)
c oblock(17)  - mulliken analysis
c keyword: mull
      integer PUN_MULLIKEN
      parameter (PUN_MULLIKEN=17)
c oblock(18)  - lowdin analysis
c keyword: lowd
      integer PUN_LOWDIN
      parameter (PUN_LOWDIN=18)
c oblock(19)  - spin densities
c keyword: spin
      integer PUN_SPIN_DENS
      parameter (PUN_SPIN_DENS=19)
c oblock(20)  - scf type
c keyword: leve
      integer PUN_SCF_TYPE
      parameter (PUN_SCF_TYPE=20)
c oblock(21)  - normal coordinates for vibrational modes
c keyword: norm
      integer PUN_VIB_MODES
      parameter (PUN_VIB_MODES=21)
c oblock(22)  - two electron integrals
c keyword: twoe
      integer PUN_2E_INT
      parameter (PUN_2E_INT=22)
c oblock(23)  - gradients
c keyword: grad
      integer PUN_GRADIENTS
      parameter (PUN_GRADIENTS=23)
c oblock(24)  - transformation matrix
c keyword: tran
      integer PUN_TRANS_MAT
      parameter (PUN_TRANS_MAT=24)
c oblock(25)  - potential derived charges
c keyword: pdc
      integer PUN_POT_DERV_CHG
      parameter (PUN_POT_DERV_CHG=25)
c oblock(26)  - grid symmetry array
c keyword: gsym
      integer PUN_GRID_SYMM
      parameter (PUN_GRID_SYMM=26)
c oblock(27)  - cartesian hessian matrix
c keyword: secd
      integer PUN_HESSIAN
      parameter (PUN_HESSIAN=27)
c oblock(28)  - distributed multipole analysis
c keyword: dma
      integer PUN_MPOLE_ANAL
      parameter (PUN_MPOLE_ANAL=28)
c oblock(29)  - density matrix
c keyword: dens
      integer PUN_DENS_MAT
      parameter (PUN_DENS_MAT=29)
c oblock(30)  - total energy
c keyword: ener
      integer PUN_TOT_ENER
      parameter (PUN_TOT_ENER=30)
c oblock(31)  - internal coordinate hessian matrix
c keyword: hess
      integer PUN_ZMAT_HESS
      parameter (PUN_ZMAT_HESS=31)
c oblock(32)  - dipole moment
c keyword: dipo
      integer PUN_DIPOLE
      parameter (PUN_DIPOLE=32)
c oblock(33)  - infrared intensities 
c keyword: infr
      integer PUN_INFRARED
      parameter (PUN_INFRARED=33)
c oblock(34)  - raman intensities
c keyword: rama
      integer PUN_RAMAN
      parameter (PUN_RAMAN=34)
c CURRENTLY ONLY USED FOR XML
c oblock(35)  - metadata (user, compilation date etc)
c keyword: meta
      integer PUN_METADATA
      parameter (PUN_METADATA=35)
c CURRENTLY ONLY USED FOR XML
c oblock(36)  - input parameters
c keyword: inpa
      integer PUN_INPUT_PARAM
      parameter (PUN_INPUT_PARAM=36)


      integer LENPUN
      parameter(LENPUN=60)

      logical oblock, opunold
      integer nfblck, nblsec, iblsec, nblgri, iblgri
      integer itwoe, ntwoe, iblfmt
      common/blocks/nfblck,nblsec,iblsec(10),nblgri,iblgri(10),
     +              itwoe,ntwoe,iblfmt(LENPUN),oblock(LENPUN),opunold
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
      integer ianz, iz, lbl, lalpha, lbeta, nz, nvar
      real*8 bl, alpha, beta
      common /czmat/ ianz(maxnz),iz(maxnz,4),bl(maxnz),
     +               alpha(maxnz),beta(maxnz),lbl(maxnz),
     +               lalpha(maxnz),lbeta(maxnz),nz,nvar
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      real*8 cicoef, f, aaaaa, bbbb
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),aaaaa(325),bbbb(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      real*8 fieldx, fieldy, fieldz
      logical ofield
      common /gms_field/  fieldx, fieldy, fieldz, ofield
c
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
c
c
      real*8 func0
      integer ncoord, npts, nserch, iupdat, icode, iupcod
      common /seerch/ func0,ncoord,npts,nserch,iupdat,icode,iupcod
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
      common/l/llock
      common/scfopt/maxcyc,mconv,nconv,npunch,
     * accdi1,accdi2,odiis,icoupl,ifolow,irotb,
     * dmpcut,acccc(6),iaccv(2),
     * rshift(4),iextn,junks,dampoo(9),idiisf,
     * optester,odynamic,omaxcyc,onoor,otestd
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
      common /copt / rcvopt
      common/foropt/vibsiz,nvib
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
c
      real*8 alphas,chargat,spinat
      integer nswapa, nswapb, nswap, next,isecat
      integer natconf,nnatconf,iatcon,iatstates,natdiff
      parameter (nnatconf=20)
      logical odiff,oexc,oground,oatdo,oalway,uhfatom,nonrel,forcat
      character*8 zatconf,atmode,zatdiff
      character*60 string_ch
      character*2 zatstate
c
      common /datgue/ alphas(maxorb),nswapa,nswapb,nswap(80),
     +                next(maxorb),oground,oatdo,isecat,oalway,
     +                odiff,oexc
      common /datgue2/ chargat(nnatconf),spinat(2,4,nnatconf),
     +                 natconf,uhfatom,iatcon,nonrel(nnatconf),
     +                 forcat(nnatconf),iatstates(nnatconf),
     +                 natdiff
      common /datgue3/ string_ch(2,nnatconf),zatconf(nnatconf),atmode,
     +                 zatstate(nnatconf),zatdiff(nnatconf)
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      common/saveco/uuuu,inttyp(5)
c
c     ------ scrf common blocks
c
      real*8 aradix,aradiy,aradiz, dielec
      logical oscrf
      common/scrf/aradix,aradiy,aradiz,dielec,oscrf
      common/dplrfc/tele(3),tnuc(3),tmol(3),dtot
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
c
c ----- Pisa solvation common blocks
c
      real*8 ptspac, delect, facnrm
      integer iconv,inkblk,itotbq,iscf,ispace,msecp
      logical opssc
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
c
c ----- crystal field
c
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst

c
c ----- omitcharge common block
c
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
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
c
      real*8 cutomp
      common /mpcut/ cutomp
c
c
c ... fortran data streams 
c
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      real*8 sco,pco,dco,fco,gco
      integer ilopri,ngroup,iatom
      integer mng,mpergr
      parameter (mng=10,mpergr=30)
      common /groups/ ilopri,ngroup,iatom(mpergr,mng),
     +                sco(mng), pco(mng), dco(mng), fco(mng),
     +                gco(mng)
c
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
c
      integer mp2gcycl
      logical opg_grad, opg_alloc, opg_alloc_sch, opg_alloc_dmat
      logical opg_grad_com
      logical opg_pertbns_sep
      common/mp2grd/mp2gcycl, opg_grad, opg_alloc, opg_alloc_sch, 
     +              opg_alloc_dmat, opg_grad_com,
     +              opg_pertbns_sep
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
c ax bx cx: limits of the bracket
c fa fb fc: energies at ax bx cx
c bta: exponets ratio
c nexp: number of even tempered expanded exponents
c outexp: outermost core exponent
c juniq:  number of symmetry unique atoms
c tnprim: number of primitives per atom
c
       real*8 tod_ax,tod_bx,tod_cx,tod_fa,tod_fb,tod_fc,tod_bta
       integer nexp,outexp,ipoint,juniq,tnprim
       logical plusmin,optbas
       common/coptbas/tod_ax,tod_bx,tod_cx,tod_fa,tod_fb,
     +              tod_fc,tod_bta,plusmin
     +              ,ipoint,optbas,nexp,outexp,juniq,tnprim
c
c    indx needs to be rest for taskfarming
      common/indxsv/indx
      integer idpa
      common/dppopc/idpa

c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
c
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
      integer ncore,nact,nels
      common /fccwfn2 / ncore,nact,nels
c
c    initialisation of structchk for bond distance checking
c
      real*8 scalel, vbond1, vbond2
      logical opoint1,ochkdst,ofatal
      integer mxbnd, nbond1, nbond2, mbond1, mbond2
      parameter (mxbnd = 1000)
c
      common/structchk/ vbond1(mxbnd),mbond1(2,mxbnd),
     +                  vbond2(mxbnd),mbond2(2,mxbnd),
     +                  scalel, nbond1, nbond2, opoint1,
     +                  ochkdst, ofatal
c
c ...
      data mit/200/
      data tenm2/1.0d-03/
      data zc1/'c1'/
c
c     set up initial default options
c
      nav = lenwrd()
c
c ... /coptbas/
c
      optbas=.false.
c
c
c ... /timez/ and /trntim/
c
      safrun = 2.0d0
      facrun = 0.85d0
      timlim = (timlim-safrun)*facrun
      safe = 10.0d0
      dumtim = 1800.0d0
c
      alimin = timlim
c
c .. /fpinfo/
c
      ocifp = .false.
c
c .. /drfopt/
c
      oreact = .false.
c
c .. /chgcc/
c
      omtchg = .false.
      omtslf = .true.
      obqf   = .false.
      do 55 loop =1,mxomat
      ztagcc(loop)=' '
      do loop2 = 1,mxombq
         ztagcq(loop,loop2) = ' '
      enddo
 55   continue
c
c ... /statis/
c
      do 20 loop = 1 , 50
         timsec(loop) = 0.0d0
         walsec(loop) = 0.0d0
 20   continue
      call cpuwal(begin,ebegin)
c
c ... /root/
c
      do 30 loop = 1 , 60
         dji(loop) = 1.0d0/(loop+loop-1)
 30   continue
c
c ... /mapper/
c
      do 40 i = 1 , maxorb
         mapie(i) = i
         i4096(i) = i*256
         k = i*(i-1)/2
         iky(i) = k
         ikyp(i) = k + i
 40   continue
c
c ... /datgue/
c
      nswapa = 0
      nswapb = 0
c
c ... /data1/
c
      mulla = 0
      omulla = .false.
      omullo = .false.
      omullp = .false.
      omullg = .false.
      omullu = .false.
      oactrn = .false.
      call setsto(maxorb,0,mullac)
      nacta = 0
      nactb = 0
      newpro = 0
      nocp = 0
      mcen = 0
      norbtp = 0
      norbc = 0
      norbcp = 0
      ihamcp = 0
      ihamc = isect(466)
      norb2 = -1
      norb3 = -1
      norbg = 0
c
c ... /cndx41/
c
      ionsv = 0
c
c ... /atmol3/
c
      call vclr(symtag,1,9)
      do 95 loop = 1,50
95    odirec(loop) = .false.
      mina = 0
      minb = 0
      mouta = 0
      moutb = 0
      gapa1 = 0.0d0
      gapa2 = 0.0d0
      gapb1 = 0.0d0
      gapb2 = 0.0d0
      scaleg = 0.0d0
      ibrk = 999
      lock = 2
      nsymo = 0
      oming = .false.
      oextg = .false.
      ovcfre = .false.
      osemi = .false.
      ostopm = .false.
      olevd = .false.
      isunor = 0
      osym_op = .false.
      osym_dens = .false.
c
c ... /runlab/
c
      zsymm = ' '
      call setstc(10,' ',ztitle)
      zgroup = zc1
      zconf = 'fix'
      zscftp = 'rhf'
      zruntp = zrunt(1)
c
c ... /sector/
c
      call secini(0,0)
c
c ... /infoa/
c
      ich = 0
      ne = 0
      mul = 1
      call setsto(maxat,0,ipseud)
      lpseud = 0
c
c ... /czmat/
c
      nz = 0
      nvar = 0
c
c ... /natorb/
c
      ispacg = 0
      isping = 0
c
c ... /mp2grd/
c
      opg_pertbns_sep = .false.
c
c ... /restar/
c
      call setsto(mxtask,-1,itask)
      local = 0
      intg76 = 2
      normf = 0
      normp = 0
      nprint = 0
      itol = 20
      icut = 9
      nopk = 0
      nopkr = 0
      iofsym = 0
      iofrst = 0
      ist = 1
      jst = 1
      kst = 1
      lst = 1
      nrec = 1
      nindmx = 0
      omaxb = .false.
      imaxb_ic = 0
      irest = 0
c
c ... /molsym/
c
      do 50 loop = 1 , 3
         prmoms(loop) = -1.0d0
 50   continue
      jaxis = 0
c ...
c     mechanics option. /modj/ etc
c ...
c     flag if hessian is to be generated by mechanics
      oamhes = .false.
      oaminc = .false.
      omapp = .false.
      omodel = .false.
      ifree = -1
c topology file format. unformatted ( binary ).
      kform = 1
c not constant pressure for dynamics.
      opress = .false.
c constant temperature for dynamics.
      otemp = .true.
c modify charges of ab initio atoms using input from mapp directive.
      ocharg = .false.
c regeneration of z-matrix variables.
      osubs = .false.
c
c ... /moda/
c
c final & restart coordinates file format. formatted.
      ntxo = 0
c read from unit 21 ( formatted).
      ntx = 0
      nrc = 0
      nrcx = 0
      npm = 1
      nram = 0
      nsm = 0
c periodicity with truncated octahedran (rectangular) constant
c  temperature.
      ntb = 0
c angle between x and z axes for monoclinic box.
      bbeta = 90.0d0
c box dimensions for dynamcis.
      box(1) = 0.0d0
      box(2) = 0.0d0
      box(3) = 0.0d0
c
c ... /modmin/
c
      ntrun = 2
      ntpr = 5
c number of mm optimisation cycles.
      imcyc = 9999
c
c ... /modf/
c
      ntp = 0
      ntnb = 1
      ntid = 0
      ipro = 1
c complete interaction is used in force evaluation.
      ntf = 1
c nb interactions are calculated using residue-based cutoff.
c                   the nb pairs are stored as atom pairs.
      ntn = 3
c update non - bonded pair list every nsnb steps.
      nsnb = 100
c
c ... /mod2d/
c
      isw = 0
c
c ... /modfov/
c
      nrcp = nrc
      npro = ipro
c
c ... /modmd/
c
c nstlim: this may cause problems later.
c may be necessary to have 2 values.
      nstlim = imcyc
c dump intermediate output during mm opt every cycle.
      mdpr = ntpr
c temperature for dynamics . initial velocities
      tempi = 50.0d0
c temperature for dynamics if constant temperature.
      tempo = 500.0d0
c coupling constant for maintaining constant temperature.
      tautp = 0.4d0
c initial velocities in maxwellian distribution.
      init = 0
c dynamics time step ( expressed in picoseconds).
      tims = 0.002d0
c
c ... /nbterm/
c
c cutoff distance for nb interactions.(angstroms may be wrong units)
      cut = 15.0d0
c distance - dependent dielectric function.
      idiel = 0
c dielectric constant if idiel .ne.0 maybe change this later.
      dielc = 1.0d0
c scale factor for 1-4 van der waals interactions.
      scnb = 2.0d0
c scale factor for 1-4 electrostatic interactions.
      scee = 2.0d0
c constant pressure if set.
c     press = 101325.0d0
c
c ... /infob/
c
      ozmat = .false.
c
c ... /restrj/
c
      ltask = 0
      lcoord = 0
      opass1 = .false.
      opass2 = .false.
      opass3 = .false.
      opass4 = .false.
      opass5 = .false.
      opass6 = .false.
      offsym = .false.
      orege = .false.
      orgall = .false.
      opass8 = .false.
      opass9 = .false.
      opass10 = .false.
      opass11 = .false.
      odisc = .false.
      oirc = .false.
      orpa = .false.
      odirpa = .false.
      omclr = .false.
      iplot = 0
      idma = 0
cdrf
      idpa = 0
cdrf
      iprop = 0
      imull = 0
      omrdci = .false.
      ofulci = .false.
      omech = .false.
      occsd = .false.
      occsdt = .false.
      oqcisd = .false.
      oqcisdt= .false.
      opark = .false.
      oclunp = .false.
      odiesel = .false.
      osortmr = .false.
      osort2  = .false.

      odlfind = .false.
c
c ... /sortp/
c

      osortp = .false.

      oschw  = .true.
c
c ... /dm/ -- default section numbers
c
      igrad = 0
      ifout1 = 0
      ifout2 = 0
      ifmola = 0
      if2 = 0
      if4 = 0
      iflagr = 0
      isecdm = isect(481)
      isecda = isect(488)
      iseclm = isect(480)
      isecla = isect(487)
      iseclg = isect(500)
c
c  /xfield/
c
      isecfx = isect(460)
      ocryst = .false.
c
c ----- setup constants /phycon/
c
      call phyfil
      ifau = 0
      nosymm = 0
      iseczz = isect(420)
c
c ... /field/
c
      fieldx = 0.0d0
      fieldy = 0.0d0
      fieldz = 0.0d0
      ofield = .false.
c
c ... /zjorgs/
c
      ystat(1) = 'mini'
      ystat(2) = 'sadd'
      ystat(3) = 'bran'
c
c ... /jorgs/
c
      maxsta = 3
      maxj = 40
      irecal = 0
      opowel = .false.
      obfgs = .true.
      obfgx = .false.
      onrfo = .true.
      ocut = .false.
      outjor = .true.
      outdeb = .false.
c
c ... /cntl1/
c
      eigmin = 1.0d-3
      eigmax = 25.0d0
      accin(1) = 0.3d0
      accin(2) = 0.003d0
      accin(3) = 0.2d0
      accin(4) = 0.1d0
      accin(5) = 0.15d0
      accin(6) = 0.10d0
      isadle = 1
      jm = mit
      mnls = mit
      iterch = mit
      ofcm = .false.
      nftotl = 0
      lintyp = 0
      lqstyp = 4
c
c ... /tol/
c
      toler = 1.0d-5
      toler2 = 2.5d-7
      isymtl = 0
c
c ... /iofile/
c
      idaf = 4
      num8 = 8
      ibl3d = 1
      numlib = 1
      iblkpl = 1
      numdis = 0
      ibldis = 1
c
c ... /dump3/
c
c ----- setup default section numbers
c
      ionsec = isect(492)
      ions2 = isect(482)
      isecqm = isect(467)
c
c ... /cndx40/
c
      npas41 = 1
      npas42 = 1
      iacc4 = 10
c
c ... /mpcut/
c
      cutomp = 1.0d-10
c 
c ... /ftape/
c
      nf1 = 1
      nf22 = 22
      nf48 = 48
      nf52 = 62
      nf58 = 58
      nf41 = 41
      nf2 = 2
      nf3 = 3
      nf4 = 4
      nf8 = 8
      nf9 = 9
      nf10 = 10
      nf31 = 31
      nf32 = 32
      nf33 = 33
      nf34 = 34
      nf35 = 35
      nf36 = 36
      nf38 = 38
      nf39 = 39
      nf42 = 42
c
c ... /prints/
c
      ciprnt = 0.01d0
      do 60 loop = 1 , 60
         oprint(loop) = .false.
 60   continue
c
c ... symmetry input
c
      osim1=.false.
      osim2=.false.
c
c ... symmetry assignment /fsymas/
c
      otsym=.true.
      degecr=5.0d-6
      omydbg=.false.
      oingam=.true.
      ostart=.true.
      iscsym=isect(408)
c
c ... lowdin analysis /groups/
c
      ilopri=0
      ngroup=0
c
c ... /blocks/
c
      do 65 loop = 1 , LENPUN
         oblock(loop) = .false.
         iblfmt(loop) = 0
 65   continue
      opunold = .false.
      nblsec = 0
      nblgri = 0
      do 66 loop = 1 , 10
       iblsec(loop) = 0
       iblgri(loop) = 0
 66   continue
      itwoe = 1
      ntwoe = 0
c
c ... /morokuma/
c
      imode = 0
      omorok = .false.
c
c ... /scfwfn/
c
      nseto = 0
      ncores = 0
      npair = 0
      nset = 0
      nopen = 0
      nope = 0
      do 70 loop = 1 , 10
         no(loop) = 0
 70   continue
      old = .false.
c
c ... /seerch/
c
      func0 = 0.0d0
      nserch = 0
      iupcod = 0
c
c ... /l/
c
      llock = 0
c
c ... /files/
c
      n2file = 1
      n4file = 1
      n6file = 1
      n5file = 1
      n9file = 1
      ntfile = 1
      n1file = 1
      n11fil = 1
      n12fil = 1
      n13fil = 1
      do 80 loop = 1 , 20
         n2blk(loop) = 1
         n4blk(loop) = 1
         n6blk(loop) = 1
         n5blk(loop) = 1
         n9blk(loop) = 1
         ntblk(loop) = 1
         n1blk(loop) = 1
         n11bl(loop) = 1
         n12bl(loop) = 1
         n13bl(loop) = 1
 80   continue
      do 90 loop = 1 , 20
         n2last(loop) = 0
         n4last(loop) = 0
         n6last(loop) = 0
         n5last(loop) = 0
         n9last(loop) = 0
         ntlast(loop) = 0
         n1last(loop) = 0
         n11lst(loop) = 0
         n12lst(loop) = 0
         n13lst(loop) = 0
 90   continue
      do 100 loop = 1 , 20
         n2tape(loop) = 3
         n4tape(loop) = 5
         n6tape(loop) = 7
         n5tape(loop) = 6
         n9tape(loop) = 10
         nttape(loop) = 11
         n1tape(loop) = 2
         n11tap(loop) = 12
         n12tap(loop) = 13
         n13tap(loop) = 14
 100  continue
      do 110 loop = 1 , maxlfn
         maxbl(loop) = 0
         minbl(loop) = 0
         maxset(loop) = 0
 110  continue
c
c ... /scfopt/
c
      odiis = .true.
      accdi1 = 0.1d0
      accdi2 = 0.0d0
      dmpcut = 0.0d0
      maxcyc = 50
      mconv = 5
      nconv = 5
      ifolow = 0
      npunch = 0
      idiisf = 3
      optester = .false.
      odynamic = .false.
      omaxcyc  = .false.
      onoor = .false.
      otestd = .false.
c
c ... /tran/
c
      otran = .false.
      otri = .true.
c
c ... /copt/
c
      rcvopt = 1.0d-3
c
c ... /foropt/
c
      vibsiz = tenm2
      nvib = 1
c
c ... /cslosc/
c
      odscf = .false.
      odnew = .false.
      oimag = .false.
      dlntol = -dlog(1.0d-8)
      tolitr(1) = -dlog(1.0d-8)
c
c  these settings were originally employed but in certain
c  cases led to convergence problems at 'low' testers ..
c  it may be that the tighter rotated axis routines will have
c  cured this problem, but for the time being we'll play safe.
c     tolitr(1)=-alog(1.0e-6)
c     tolitr(2)=-alog(1.0e-7)
c
      tolitr(2) = -dlog(1.0d-8)
      tolitr(3) = -dlog(1.0d-8)
      itrtol(1) = 3
      itrtol(2) = 6
      itrtol(3) = 999
c     deltol = dlog(0.05d0)
c ... in the interests of safety, not efficiency ...
      deltol = dlog(0.00000001d0)
      delfac = dlog(0.1d0)
      odelta = .false.
c
c ... /unocas/
c
      iuno = 0
      iunopp = 1
      iunosp = -1
      iunspp = 1
      oanil = .false.
c
c ... /fccwfn/ ormas common
c
      omas = .false.
      omasd = .false.
      nspace = 0
      ncore  = 0
      nact   = 0
      nels   = 0
c
c ... /scra7/
c
      ibl7la = 1
c
c ... /work/
c
      jrec = 0
c
c ... /saveco/
c
      do 160 loop=1,5
 160  inttyp(loop) = -1
c
c ... /wrtd/
c
      oentry = .false.
c
c ... /scrf/ and /psscrf/ (pisa common blocks)
c
      ptspac = 0.0d0
      delect = 0.0d0
      opssc = .false.
      iconv = 5
      itotbq = 0
      iscf = 0
      msecp = isect(409)
c
      oscrf = .false.
      aradix = 0.0d0
      aradiy = 0.0d0
      aradiz = 0.0d0
      dielec = 0.0d0
c
c   /gen/ ,/gen2/ and /gen3/  (direct-scf)
c
      iaccur = 4
      ocvary = .true.
      olog1 = .true.
      olog2 = .true.
      olog9 = .false.
      olog6 = .false.
      olog5 = .false.
      do 130 i = 1 , 10
         odbg(i) = .false.
 130  continue
      do 140 i = 1 , 21
         prefac(i) = 1.0d-8
         tinner(i) = 1.0d-8
 140  continue
      do 170 j = 1 , 2
         do 180 i = 1 , 21
            ripfal(i,j) = 1.0d0
            rijkal(i,j) = 1.0d0
            ripftd(i,j) = 1.0d0
            rijkld(i,j) = 1.0d0
 180     continue
 170  continue
      rifall = 1.0d0
      riffew = 1.0d0
      iswap = 0
      timeg = 0.0d0
      ogmp2 = .false.
      do 150 loop = 1 , 5
         ogrdt(loop) = .false.
 150  continue
      ogmull = .false.
      noptd = 4
      mp2fo = 0
      mp2fv = 0
      do 190 i = 1,100
 190  ofull(i) = .false.
      odiisn = .true.
      oread(1) = .false.
      oread(2) = .false.
      oread(3) = .false.
c
c   /gjs/
c
      symm_diag = .true.
c
c   /gvalue/
c
      gx = 0.0d0
      gy = 0.0d0
      gz = 0.0d0
      orgin = .false.
      do i = 1, 3
       ggx(i) = 0.0d0
      enddo
c ... taskfarming  intialisations
      indx = 0
      indx_xml = 0
c
c   /structchk/ for bond distance checking
c
      ochkdst = .false.
      opoint1 = .true.
      ofatal = .false.
      scalel = 1.3d0
c
c     /common_llvmo/
c
      call llvmo_defaults
      return
      end
      subroutine initialb
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
c symmetry assignment
      common/atmol3o/osim1,osim2
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
c
c
c
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
c
      real*8  dx0, dxm, dele, drms
      integer ntpr, ntrun, imcyc, ncyc, ntmin, nspacb
      common /modmin/ ntpr,ntrun,imcyc,ncyc,ntmin,nspacb,
     +                dx0,dxm,dele,drms
c
c
      real*8 bbeta, box
      integer ntx, ntxo, nrc, nrcx, npm, nrp, nsm, nram, ntb
      integer nspaca
      common /moda/ ntx,ntxo,nrc,nrcx,npm,nrp,nsm,nram,ntb,nspaca,
     +  bbeta,box(3)
c
c
      integer ntp, ntc, ntcc, ntf, ntid, ntn, ntnb, nsnb
      integer nrep, nrepc, ipro
      common /modf/ ntp,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,nrep,
     +              nrepc,ipro
c
c
      real*8 tempi, tempo, tautp, t, tims
      integer nstlim, mdpr, init, nspacc
      common /modmd/ nstlim,mdpr,init,nspacc,
     +               tempi,tempo,tautp,t,tims
c
c
      integer natom, nres, nbonh, nbona, ntheth, ntheta, nphih
      integer nphia, nnb, ntypes, nrcp, nconp, mema, memb, maxnb
      integer mbona, mtheta, mphia, natc, npro, numbnd, numang
      integer numphi, natyp, nphb, imodf, jmodf
      common /modfov/ natom,nres,nbonh,nbona,ntheth,ntheta,nphih,
     +                nphia,nnb,ntypes,nrcp,nconp,mema,memb,maxnb,
     +                mbona,mtheta,mphia,natc,npro,numbnd,numang,
     +                numphi,natyp,nphb,imodf(25),jmodf(50)
c
c
      integer nword, isw, kprint, ifzmat
      common /mod2d/ nword,isw,kprint,ifzmat
c
      common/nbterm/cut,scnb,scee,dielc,idiel,nbuck
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
c
      common /zjorgs/ystat(3)
      common /jorgs/maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     *outjor,outdeb,maxsta
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c The oblock array determines what is written out to the punchfile:
c
c Below sourced from server.m
c oblock(1)   - coords - single point or last point
c keyword: coor
      integer PUN_COORD
      parameter (PUN_COORD=1)
c oblock(2)   - coordinates at each optimisation step
c keyword: opti
      integer PUN_OPT_COORD
      parameter (PUN_OPT_COORD=2)
c oblock(3)   - starting coordinates
c keyword: init
      integer PUN_START_COORD
      parameter (PUN_START_COORD=3)
c oblock(4)   - atom connectivity
c keyword: conn
      integer PUN_ATOM_CONN
      parameter (PUN_ATOM_CONN=4)
c oblock(5)   - molecular orbital occupations
c keyword: occu
      integer PUN_MO_OCC
      parameter (PUN_MO_OCC=5)
c oblock(6)   - eigenvectors
c keyword: vect
      integer PUN_EIGV
      parameter (PUN_EIGV=6)
c oblock(7)   - run type
c keyword: type
      integer PUN_RUNT
      parameter (PUN_RUNT=7)
c oblock(8)   - gvb pair data
c keyword: gvb
      integer PUN_GVBPAIR
      parameter (PUN_GVBPAIR=8)
c oblock(9)   - force constant matrix
c keyword: forc
      integer PUN_FCMAT
      parameter (PUN_FCMAT=9)
c oblock(10)  - vibrational frequencies
c keyword: vibr
      integer PUN_VIBFREQ
      parameter (PUN_VIBFREQ=10)
c oblock(11)  - basis set
c keyword: basi
      integer PUN_BASIS
      parameter (PUN_BASIS=11)
c oblock(12)  - orbital energies
c keyword: eige
      integer PUN_ORB_ENER
      parameter (PUN_ORB_ENER=12)
c oblock(13)  - total scf energies
c keyword: scfe
      integer PUN_SCF_ENER
      parameter (PUN_SCF_ENER=13)
c oblock(14)  - job title
c keyword: titl
      integer PUN_JOB_TITLE
      parameter (PUN_JOB_TITLE=14)
c oblock(15)  - overlap matrix
c keyword: over
      integer PUN_OVLP_MAT
      parameter (PUN_OVLP_MAT=15)
c oblock(16)  - grid data from dumpfile section(s)...
c keyword: grid
      integer PUN_GRID_DATA
      parameter (PUN_GRID_DATA=16)
c oblock(17)  - mulliken analysis
c keyword: mull
      integer PUN_MULLIKEN
      parameter (PUN_MULLIKEN=17)
c oblock(18)  - lowdin analysis
c keyword: lowd
      integer PUN_LOWDIN
      parameter (PUN_LOWDIN=18)
c oblock(19)  - spin densities
c keyword: spin
      integer PUN_SPIN_DENS
      parameter (PUN_SPIN_DENS=19)
c oblock(20)  - scf type
c keyword: leve
      integer PUN_SCF_TYPE
      parameter (PUN_SCF_TYPE=20)
c oblock(21)  - normal coordinates for vibrational modes
c keyword: norm
      integer PUN_VIB_MODES
      parameter (PUN_VIB_MODES=21)
c oblock(22)  - two electron integrals
c keyword: twoe
      integer PUN_2E_INT
      parameter (PUN_2E_INT=22)
c oblock(23)  - gradients
c keyword: grad
      integer PUN_GRADIENTS
      parameter (PUN_GRADIENTS=23)
c oblock(24)  - transformation matrix
c keyword: tran
      integer PUN_TRANS_MAT
      parameter (PUN_TRANS_MAT=24)
c oblock(25)  - potential derived charges
c keyword: pdc
      integer PUN_POT_DERV_CHG
      parameter (PUN_POT_DERV_CHG=25)
c oblock(26)  - grid symmetry array
c keyword: gsym
      integer PUN_GRID_SYMM
      parameter (PUN_GRID_SYMM=26)
c oblock(27)  - cartesian hessian matrix
c keyword: secd
      integer PUN_HESSIAN
      parameter (PUN_HESSIAN=27)
c oblock(28)  - distributed multipole analysis
c keyword: dma
      integer PUN_MPOLE_ANAL
      parameter (PUN_MPOLE_ANAL=28)
c oblock(29)  - density matrix
c keyword: dens
      integer PUN_DENS_MAT
      parameter (PUN_DENS_MAT=29)
c oblock(30)  - total energy
c keyword: ener
      integer PUN_TOT_ENER
      parameter (PUN_TOT_ENER=30)
c oblock(31)  - internal coordinate hessian matrix
c keyword: hess
      integer PUN_ZMAT_HESS
      parameter (PUN_ZMAT_HESS=31)
c oblock(32)  - dipole moment
c keyword: dipo
      integer PUN_DIPOLE
      parameter (PUN_DIPOLE=32)
c oblock(33)  - infrared intensities 
c keyword: infr
      integer PUN_INFRARED
      parameter (PUN_INFRARED=33)
c oblock(34)  - raman intensities
c keyword: rama
      integer PUN_RAMAN
      parameter (PUN_RAMAN=34)
c CURRENTLY ONLY USED FOR XML
c oblock(35)  - metadata (user, compilation date etc)
c keyword: meta
      integer PUN_METADATA
      parameter (PUN_METADATA=35)
c CURRENTLY ONLY USED FOR XML
c oblock(36)  - input parameters
c keyword: inpa
      integer PUN_INPUT_PARAM
      parameter (PUN_INPUT_PARAM=36)


      integer LENPUN
      parameter(LENPUN=60)

      logical oblock, opunold
      integer nfblck, nblsec, iblsec, nblgri, iblgri
      integer itwoe, ntwoe, iblfmt
      common/blocks/nfblck,nblsec,iblsec(10),nblgri,iblgri(10),
     +              itwoe,ntwoe,iblfmt(LENPUN),oblock(LENPUN),opunold
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      real*8 fieldx, fieldy, fieldz
      logical ofield
      common /gms_field/  fieldx, fieldy, fieldz, ofield
c
c
      real*8 func0
      integer ncoord, npts, nserch, iupdat, icode, iupcod
      common /seerch/ func0,ncoord,npts,nserch,iupdat,icode,iupcod
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
      common/l/llock
      common/scfopt/maxcyc,mconv,nconv,npunch,
     * accdi1,accdi2,odiis,icoupl,ifolow,irotb,
     * dmpcut,acccc(6),iaccv(2),
     * rshift(4),iextn,junks,dampoo(9),idiisf,
     * optester,odynamic,omaxcyc,onoor,otestd
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common /copt / rcvopt
      common/foropt/vibsiz,nvib
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
      real*8 alphas,chargat,spinat
      integer nswapa, nswapb, nswap, next,isecat
      integer natconf,nnatconf,iatcon,iatstates,natdiff
      parameter (nnatconf=20)
      logical odiff,oexc,oground,oatdo,oalway,uhfatom,nonrel,forcat
      character*8 zatconf,atmode,zatdiff
      character*60 string_ch
      character*2 zatstate
c
      common /datgue/ alphas(maxorb),nswapa,nswapb,nswap(80),
     +                next(maxorb),oground,oatdo,isecat,oalway,
     +                odiff,oexc
      common /datgue2/ chargat(nnatconf),spinat(2,4,nnatconf),
     +                 natconf,uhfatom,iatcon,nonrel(nnatconf),
     +                 forcat(nnatconf),iatstates(nnatconf),
     +                 natdiff
      common /datgue3/ string_ch(2,nnatconf),zatconf(nnatconf),atmode,
     +                 zatstate(nnatconf),zatdiff(nnatconf)
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      common/saveco/uuuu,inttyp(5)
      common/gen2/ripftd(21,2),rijkld(21,2),rjpftd(21,2),
     + ripfal(21,2),rijkal(21,2),
     + prefac(21),tinner(21),rifall,riffew,ricfb1,ricfb2
c
c ----- omitcharge common block
c
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
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
c
      real*8 cutomp
      common /mpcut/ cutomp
c
      common /multic/radius,trust1,tfac1,trust2,tfac2,sparse,
     + convmc,econv
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c    indx needs to be rest for taskfarming
      common/indxsv/indx
c
      real*8 scalel, vbond1, vbond2
      logical opoint1,ochkdst,ofatal
      integer mxbnd, nbond1, nbond2, mbond1, mbond2
      parameter (mxbnd = 1000)
c
      common/structchk/ vbond1(mxbnd),mbond1(2,mxbnd),
     +                  vbond2(mxbnd),mbond2(2,mxbnd),
     +                  scalel, nbond1, nbond2, opoint1,
     +                  ochkdst, ofatal
c
      integer idpa
      common/dppopc/idpa

c
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
      integer ncore,nact,nels
      common /fccwfn2 / ncore,nact,nels
c
      data mit/200/
      data tenm2/1.0d-03/
c
c...  get rid of lingering  aimpac stuff
c
      call aimpac_request("off",0)
c
c     set up initial default options
c
      call cpuwal(begin,ebegin)
c
c ... /mapper/
c
      do 20 i = 1 , maxorb
         mapie(i) = i
         i4096(i) = i*256
         k = i*(i-1)/2
         iky(i) = k
         ikyp(i) = k + i
 20   continue
c
c ... /datgue/
c
      nswapa = 0
      nswapb = 0
c
c ... /data1/
c
      mulla = 0
      omulla = .false.
      omullo = .false.
      omullp = .false.
      omullg = .false.
      omullu = .false.
      call setsto(maxorb,0,mullac)
      nacta = 0
      nactb = 0
      newpro = 0
      nocp = 0
      mcen = 0
      norbtp = 0
      norbc = 0
      norbcp = 0
      ihamcp = 0
      ihamc = isect(466)
      norb2 = -1
      norb3 = -1
      norbg = 0
c
c ... /cndx41/
c
      ionsv = 0
c
c ... /atmol3/
c
      mina = 0
      minb = 0
c     mouta=0
c     moutb=0
      gapa1 = 0.0d0
      gapa2 = 0.0d0
      gapb1 = 0.0d0
      gapb2 = 0.0d0
      scaleg = 0.0d0
      ibrk = 999
      lock = 2
      nsymo = 0
      ovcfre = .false.
      osemi = .false.
      ostopm = .false.
      olevd = .false.
      isunor = 0
      osym_op = .false.
      osym_dens = .false.
c
c ... /runlab/
c
      zscftp = 'rhf'
      zruntp = zrunt(1)
c
c ... /infoa/
c
      call setsto(maxat,0,ipseud)
c
c ... /restar/
c
      local = 0
      normf = 0
      normp = 0
      nprint = 0
      itol = 20
      icut = 9
      nopk = 0
      ist = 1
      jst = 1
      kst = 1
      lst = 1
      nrec = 1
      nindmx = 0
      omaxb = .false.
      imaxb_ic = 0
      irest = 0
c ...
c     mechanics option. /modj/ etc
c ...
c     flag if hessian is to be generated by mechanics
      oamhes = .false.
      oaminc = .false.
      omapp = .false.
      omodel = .false.
      ifree = -1
c topology file format. unformatted ( binary ).
      kform = 1
c not constant pressure for dynamics.
      opress = .false.
c constant temperature for dynamics.
      otemp = .true.
c modify charges of ab initio atoms using input from mapp directive.
      ocharg = .false.
c regeneration of z-matrix variables.
      osubs = .false.
c
c ... /moda/
c
c final & restart coordinates file format. formatted.
      ntxo = 0
c read from unit 21 ( formatted).
      ntx = 0
      nrc = 0
      nrcx = 0
      npm = 1
      nram = 0
      nsm = 0
c periodicity with truncated octahedran (rectangular) constant
c  temperature.
      ntb = 0
c angle between x and z axes for monoclinic box.
      bbeta = 90.0d0
c box dimensions for dynamcis.
      box(1) = 0.0d0
      box(2) = 0.0d0
      box(3) = 0.0d0
c
c ... /modmin/
c
      ntrun = 2
      ntpr = 5
c number of mm optimisation cycles.
      imcyc = 9999
c
c ... /modf/
c
      ntp = 0
      ntnb = 1
      ntid = 0
      ipro = 1
c complete interaction is used in force evaluation.
      ntf = 1
c nb interactions are calculated using residue-based cutoff.
c                   the nb pairs are stored as atom pairs.
      ntn = 3
c update non - bonded pair list every nsnb steps.
      nsnb = 100
c
c ... /mod2d/
c
      isw = 0
c
c ... /modfov/
c
      nrcp = nrc
      npro = ipro
c
c ... /modmd/
c
c nstlim: this may cause problems later.
c may be necessary to have 2 values.
      nstlim = imcyc
c dump intermediate output during mm opt every cycle.
      mdpr = ntpr
c temperature for dynamics . initial velocities
      tempi = 50.0d0
c temperature for dynamics if constant temperature.
      tempo = 500.0d0
c coupling constant for maintaining constant temperature.
      tautp = 0.4d0
c initial velocities in maxwellian distribution.
      init = 0
c dynamics time step ( expressed in picoseconds).
      tims = 0.002d0
c
c ... /nbterm/
c
c cutoff distance for nb interactions.(angstroms may be wrong units)
      cut = 15.0d0
c distance - dependent dielectric function.
      idiel = 0
c dielectric constant if idiel .ne.0 maybe change this later.
      dielc = 1.0d0
c scale factor for 1-4 van der waals interactions.
      scnb = 2.0d0
c scale factor for 1-4 electrostatic interactions.
      scee = 2.0d0
c constant pressure if set.
c     press = 101325.0d0
c
c ... /restrj/
c
      ltask = 0
      lcoord = 1
      opass1 = .false.
      opass2 = .false.
      opass3 = .false.
      opass4 = .false.
      opass5 = .false.
      opass6 = .false.
      offsym = .false.
      opass8 = .false.
      opass9 = .false.
      opass10 = .false.
      opass11 = .false.
      odisc  = .false.
      orege = .false.
      orgall = .false.
      oirc = .false.
      orpa = .false.
      odirpa = .false.
      omclr = .false.
      iplot = 0
      idma = 0
cdrf
      idpa = 0
cdrf
      inbo = 0
      iprop = 0
      imull = 0
      omrdci = .false.
      ofulci = .false.
      omech = .false.
      occsd = .false.
      occsdt = .false.
      oqcisd = .false.
      oqcisdt= .false.
      opark = .false.
      oclunp = .false.
c
c ... /field/
c
      fieldx = 0.0d0
      fieldy = 0.0d0
      fieldz = 0.0d0
      ofield = .false.
c
c ... /zjorgs/
c
      ystat(1) = 'mini'
      ystat(2) = 'sadd'
      ystat(3) = 'bran'
c
c .. /fpinfo/
c
      ocifp = .false.
c
c .. /drfopt/
c
      oreact = .false.
c
c .. /chgcc/
c
      omtchg = .false.
      omtslf = .true.
      obqf   = .false.
      do 55 loop =1,mxomat
      ztagcc(loop)=' '
      do loop2 = 1,mxombq
         ztagcq(loop,loop2) = ' '
      enddo
 55   continue
c
c ... /jorgs/
c
      maxsta = 3
      maxj = 40
      irecal = 0
      opowel = .false.
      obfgs = .true.
      obfgx = .false.
      onrfo = .true.
      ocut = .false.
      outjor = .true.
      outdeb = .false.
c
c ... /cntl1/
c
      eigmin = 1.0d-3
      eigmax = 25.0d0
      accin(1) = 0.3d0
      accin(2) = 0.003d0
      accin(3) = 0.2d0
      accin(4) = 0.1d0
      accin(5) = 0.15d0
      accin(6) = 0.10d0
      isadle = 1
      jm = mit
      mnls = mit
      iterch = mit
      ofcm = .false.
      nftotl = 0
      lintyp = 0
      lqstyp = 4
c
c ... /tol/
c
      toler = 1.0d-5
      toler2 = 2.5d-7
      isymtl = 0
c
c ... /iofile/
c
      numlib = 1
      iblkpl = 1
      numdis = 0
      ibldis = 1
c
c ... /cndx40/
c
      npas41 = 1
      npas42 = 1
      iacc4 = 10
c
c ... /mpcut/
c
      cutomp = 1.0d-10
c
c ... /prints/
c
      ciprnt = 0.01d0
      do 30 loop = 1 , 60
         oprint(loop) = .false.
 30   continue
c
c ... suppress printing of distance matrix
c ... on second and subsequent passes
      oprint(21) = .true.
c
c ... symmetry input
c
      osim1=.false.
      osim2=.false.
c
c ... symmetry assignment /fsymas/
c
      otsym=.true.
      degecr=5.0d-6
      oingam=.true.
      omydbg=.false.
      ostart=.true.
      iscsym=isect(408)
c
c
c  ... /blocks/
c
      do 65 loop = 1 , 60
         oblock(loop) = .false.
         iblfmt(loop) = 0
 65   continue
      opunold = .false.
      nblsec = 0
      nblgri = 0
      do 66 loop = 1 , 10
         iblsec(loop) = 0
         iblgri(loop) = 0
 66   continue
      itwoe = 1
      ntwoe = 0
c
c ... /scfwfn/
c
      nseto = 0
      ncores = 0
      npair = 0
      nset = 0
      nopen = 0
      nope = 0
      do 40 loop = 1 , 10
         no(loop) = 0
 40   continue
      old = .false.
c
c ... /fccwfn/ ormas common
c
      omas = .false.
      omasd = .false.
      nspace = 0
      ncore  = 0
      nact   = 0
      nels   = 0
c
c ... /seerch/
c
      func0 = 0.0d0
      nserch = 0
      iupcod = 0
c
c ... /l/
c
      llock = 0
c
c ... /scfopt/
c
      odiis = .true.
      accdi1 = 0.1d0
      accdi2 = 0.0d0
      dmpcut = 0.0d0
      maxcyc = 50
      mconv = 5
      nconv = 5
      ifolow = 0
      npunch = 0
      idiisf = 3
      optester = .false.
      odynamic = .false.
      omaxcyc  = .false.
      onoor = .false.
      otestd = .false.
c
c ... /copt/
c
      rcvopt = 1.0d-3
c
c ... /foropt/
c
      vibsiz = tenm2
      nvib = 1
c
c ... /cslosc/
c
      odscf = .false.
      odnew = .false.
      oimag = .false.
      dlntol = -dlog(1.0d-8)
      tolitr(1) = -dlog(1.0d-8)
c
c  these settings were originally employed but in certain
c  cases led to convergence problems at 'low' testers ..
c  it may be that the tighter rotated axis routines will have
c  cured this problem, but for the time being we'll play safe.
c     tolitr(1)=-alog(1.0e-6)
c     tolitr(2)=-alog(1.0e-7)
c
      tolitr(2) = -dlog(1.0d-8)
      tolitr(3) = -dlog(1.0d-8)
      itrtol(1) = 3
      itrtol(2) = 6
      itrtol(3) = 999
c     deltol = dlog(0.05d0)
c ... in the interests of safety, not efficiency ...
      deltol = dlog(0.00000001d0)
      delfac = dlog(0.1d0)
      odelta = .false.
c
c ... /unocas/
c
      iuno = 0
      iunopp = 1
      iunosp = -1
      iunspp = 1
      oanil = .false.
c
c ... /scra7/
c
      ibl7la = 1
c
c ... /work/
c
      jrec = 0
c
c ... /saveco/
c
      do 160 loop=1,5
 160  inttyp(loop) = -1
c
c ... /gen2/
c
      do 170 j = 1 , 2
         do 180 i = 1 , 21
            ripfal(i,j) = 1.0d0
            rijkal(i,j) = 1.0d0
            ripftd(i,j) = 1.0d0
            rijkld(i,j) = 1.0d0
 180     continue
 170  continue
      rifall = 1.0d0
      riffew = 1.0d0
c
c reset potential derived charges
      odumf = opotf(2)
c
c ... /multic/ reset econv
c
      econv = 0.0d0
c
c ... taskfarming  intialisations
      indx = 0
      indx_xml = 0
c
c
c   /structchk/ for bond distance checking
c
      opoint1 = .true.
      scalel = 1.3d0
      ofatal = .false.
c
c     /common_llvmo/
c
      call llvmo_defaults
      return
      end
      subroutine integ(q)
c
c     driving routine for integrals
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
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
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
      integer ncore,nact,nels
      common /fccwfn2 / ncore,nact,nels
c
      dimension q(*)
      character*8 fnm
      character*5 snm
      data fnm,snm/"master.m","integ"/
      inull = igmem_null()
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,q)
            call putstv(q)
            if (tim.ge.timlim) return
         end if
      end if
      ionsv = 1
      if(offsym) iofsym = 1
      opass1 = .false.
      if (opass2) then
         if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         call start_time_period(TP_INTEG)
         if (irest.le.1 ) then
            if (mp2 .or. mp3 .or. lci .or. lmcscf .or. omas) nopk = 1
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "prefac",IGMEM_DEBUG)
            call rdmake(q(iprefa))
            call jandk(zscftp,q,q(inull),q(inull),q(inull),q(inull),
     +                 q(inull),q(iprefa),q(inull))
            call gmem_free_inf(iprefa,fnm,snm,"prefac")
            if (tim.ge.timlim) return
         call end_time_period(TP_INTEG)
         end if
      end if
      opass2 = .false.
      return
      end
      subroutine main2e(core)
c
c ----- generate 1-e and 2-e integrals only
c ----- all integrals to mainfile
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      character*8 fnm
      character*6 snm
      data fnm,snm/"master.m","main2e"/
c
      dimension core(*)
      data ono / .false./
c
      dum = 0.0d0
      idum = 0
      ntsave = nt
      if (offsym) iofsym = 1
      if (.not.(opass1)) then
         if (irest.lt.1) call standv(0,core)
c
c ----- 1e-integrals
c
         if (tim.ge.timlim) go to 20
      end if
c
c ----- 2e-integrals
c
      if (opass2) then
         if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if (irest.le.1) then
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "prefac",IGMEM_DEBUG)
            call rdmake(core(iprefa))
            call jandk(zscftp,core,core(inull),core(inull),core(inull),
     +                 core(inull),core(inull),core(iprefa),core(inull))
            call gmem_free_inf(iprefa,fnm,snm,"prefac")
         endif
         if (tim.ge.timlim) go to 20
      end if
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 20   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
      return
      end
      subroutine driver (core)
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
cdrf------------------------
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
      common/nottwi/obeen,obeen2,obeen3,obeen4
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
      real*8 enci, ec, eci
      integer nstat
      common /enrgci/ enci,ec,eci(10),nstat
c
      integer istate
      common /enrgc2/ istate(10)
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
      integer idpa
      common/dppopc/idpa

cdrf-------------------------
c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
c
      real*8  dx0, dxm, dele, drms
      integer ntpr, ntrun, imcyc, ncyc, ntmin, nspacb
      common /modmin/ ntpr,ntrun,imcyc,ncyc,ntmin,nspacb,
     +                dx0,dxm,dele,drms
c
c ...
      common/restrl/ociopt(2),omp2,ogr(13),
     + oforce,oci,ocart(32)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odumpm, ooprnt, opun7, omax, ologf
      integer maxbrd, logf, iofsp
      common /utilc/ odumpm,ooprnt,opun7,omax,maxbrd,ologf,
     +               logf,iofsp(4)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      common/restrz/ztttt(11),zzzzzz,zzspac(38)
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
      character *8 zdate,ztime,zaccno,zanam
      common/jinfo/zdate,ztime,zaccno,zanam
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
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
      integer mxbltyp
      parameter(mxbltyp = 50)

      logical oblur
      common/blurrr/oblur

      integer nbltyp, bltab(maxat), nutab(maxat)

      real*8 blexpo(maxat), blwght(maxat)
      real*8 blexpo2(mxbltyp), blwght2(mxbltyp)

      common/blurr2/blexpo, blwght,
     &     blexpo2, blwght2,
     &     bltab,nbltyp,nutab

      character*8 ztagbl
      common/blurrc/ztagbl(mxbltyp)

      character *8 fnm
      character *6 snm
      integer igmem_alloc_all_inf
      dimension core(*),zcas(3)
cdrf------------------------
      dimension dumm(mxgran), dumx(mxgrpar)
cdrf------------------------
      data fnm/"master.m"/
      data snm/"driver"/
      data zprop/'prop'/
      data zgame/' gamess'/ , zcas/'casscf','mcscf','vb'/
c
c     ----- set up job -----
c
      call gmem_null_initialise(core)
      call cpuwal(tstart,estart)
c
      eci(1) = 0.d00
c
      zcom(1) = zanam
      zcom(2) = zdate
      zcom(3) = ztime
      zcom(4) = zgame
      zcom(5) = ' '
      zcom(6) = zaccno
c
      call inpa(ztext)
      ytext = ztext(1:4)
c
      if (ytext.eq.'util' .or. ztext.eq.'serv') then
         call gserv(core,oquit)
         if (oquit) go to 280
      else if (ztext.eq.'servec') then
         call servec(core,'driver')
         call input
         if (.not.oswit) then
            call inpa4(ytext)
            if (ytext.eq.'stop') oswit = .true.
         end if
         if (oswit) then
            oquit = .true.
            go to 280
         end if
      end if

      if (ytext.eq.'mopa') then
         oldm = odumpm
         odumpm = .true.
         call mopac(core,oquit)
         odumpm = oldm
         if (oquit) go to 280
         call input
      end if
      jrec = 0
      oresti = .false.
      oterm = .false.
      ofirt = .true.
      oprop = .false.
c
c     ----- start- read in basis set and get initial mo*s -----
c
 20   call start(core,oresti,ostop,ostopr)
c
c punch requests table in output file
c
      call blktab(ostopr)
c
c write initial data to punchfile
c
      call blk11(core)
c
c
      if (ostop) then
c
c  end of this run
c
         if (ostopr) then
c
c user requested stop ... punch any available information
c on the dumpfile
c
            call blk(core)
 
            write (iwr,7000)
 7000       format (/' program termination requested'/)
 
         endif
 
         call closda

         go to 280
c ...
c     is this either a molecular mechanics, a molecular dynamics
c     or a mopac run ?
c ...
      else if (omodel) then
         call model(core,core,0)
c ...
c     ntmin = 3 : dynamics. molecular dynamics
c     calculation on the atoms surrounding the 'ab initio' atoms.
c ...
         if (ntmin.eq.3) then
            call model(core,core,8)
            go to 230
         end if
         if (ntmin.ne.0) go to 230
      end if
      if (omopac) then
         call mopac(core,oquit)
         go to 230
      end if
c
c      ------ check system is small enough to treat
c
      call chksiz(.false.)
c
      call vdwaals_print
c
c
c initialise CCP1 DFT module
c
      idum = CD_init(core,iwr)
      if (opg_root().and.CD_active()) then
         call CD_print_joboptions(iwr)
      endif

      if(oblur)then
c
c set up the external gaussian distribution
c
         if(nbltyp .ne. 0)then
c
c input via atom labels (blur directive)
c            
            do i = 1, nat
               blexpo(i) = -1.0d0
               blwght(i) =  0.0d0
               do j = 1, nbltyp
                  if(zaname(i) .eq. ztagbl(j))then
                     blexpo(i) = blexpo2(j)
                     blwght(i) = blwght2(j)
                  endif
               enddo
            enddo
         else
c
c input via blur keyword on geometry block
c
            do i = 1, nat
               write(6,*)'atom',i,'blur', blexpo(i), blwght(i)
            enddo
         endif
         idum = gden_init(core,iwr,.false.)
         call gden_add_gauss(nat,iwr,.false.)
      endif
c
c     ----- get initial m.o.s
c
      if (zruntp.ne.zrunt(10)) then
 
        call mogues(core)
        if (ostopm) go to 260
c
c     === attempt to allocate memory for integral storage
c
        if(omem(3).and.ofirt) call ed2mem(num,nat)
         ofirt = .false.
      endif
      irun = locatc(zrunt,mxtask,zruntp)
c
      osymm = .true.
      omcscf = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2) 
     &         .or. zscftp.eq.zcas(3)
      oprop = .false.
      
      if (irun.eq.0) go to 230
      go to ( 30, 40, 40, 50, 50, 60, 60, 70, 80, 90,
     +        40, 40, 50, 50,100,420, 40, 55,310, 50,
     +       110,120,130,150,140,200,160,170,180,230,
     +       210,220,225,215,230,230,230,230,230,230,
     +       230,300,227,230,230,230,230,230,230,230) , irun
c
c ==========================================================
c    if irun.eq.1     * calculate hf/casscf/mp2/mp3/masscf energy
c    if irun.eq.2     * 4-index transformation
c    if irun.eq.3     * direct-ci/full-ci/mrd-ci calculation
c    if irun.eq.4     * saddle point minimization
c    if irun.eq.5     * optimize molecular geometry
c    if irun.eq.6     * optimize molecular geometry (optxyz)
c    if irun.eq.7     * single point gradient evaluation
c    if irun.eq.8     * calculate force constant matrix
c    if irun.eq.9     * symmetry adapted transformation
c    if irun.eq.10    * wavefunction analysis
c    if irun.eq.11    * green's function (ovgf) calculation
c    if irun.eq.12    * green's function (tda) calculation
c    if irun.eq.15    * integral evaluation only
c    if irun.eq.16    * NMR chemical shifts
c    if irun.eq.17    * sac-ci calculation (no longer in use)
c    if irun.eq.18    * intrinsic reaction coordinate calc.
c    if irun.eq.19    * RPA with gradients calculation
c    if irun.eq.20    * RPA geometry optimisation (optimise)
c    if irun.eq.21    * 4-index transformation (revised)
c    if irun.eq.22    * polarizability
c    if irun.eq.23    * fockder
c    if irun.eq.24    * analytic second-derivatives (hessian)
c    if irun.eq.25    * dipole-moment derivatives (dipder)
c    if irun.eq.26    * polarisability derivatives (polder)
c    if irun.eq.27    * coupled hartree-fock magnetizability
c    if irun.eq.28    * intensity
c    if irun.eq.29    * test
c    if irun.eq.31    * hyperpolarizability
c    if irun.eq.32    * lazzer
c    if irun.eq.33    * raman intensities
c    if irun.eq.34    * infra-red intensities
c    if irun.eq.41    * mopac
c    if irun.eq.42    * response (mclr etc)
c    if irun.eq.43    * montec (drf etc)
c ==========================================================
c  
c     calculate hf/casscf/mp2/mp3/masscf energy
c
 30   call scfgvb(core)
      odscf1 = odscf
      oprop = iprop.lt.0
      go to 230
c
c     4-index transformation
c     direct-ci/full-ci/mrd-ci calculation
c
 40   call tranci(core)
      oprop = iprop.lt.0
      go to 230
c
c     hondo  4-index transformation
c
 110  call trndrv(core)
      oprn(4) = .false.
      osymm = .false.
      go to 230
c
c     saddle point minimization / optimize molecular geometry
c
 50   continue
      if (.not. odlfind) then
        call minit(core)
      else
c
c     DL-FIND, various optimisation methods
c
        call dlfind_gamess(core)
      endif
        oprop = iprop.lt.0
      go to 230
c
c     intrinsic reaction coordinate calculation
c
 55   call irc(core,iwr)
      go to 230
c
c     optimize molecular geometry (optxyz)
c     single point gradient evaluation
c
 60   continue
      if (.not. odlfind) then
        call optx(core)
      else
c
c     DL-FIND, various optimisation methods
c
        call dlfind_gamess(core)
      endif
        oprop = iprop.lt.0
      go to 230
c
c     calculate force constant matrix
c
 70   oforce = .true.
      call forcx(core)
      go to 230
c
c     symmetry adapted transformation
c
 80   call adap2e(core)
      go to 230
c
c     wavefunction analysis
c
 90   call analis(core)
      go to 230
c
c     integral evaluation only
c
 100  call main2e(core)
      go to 230
c
c     calculate polarizabilities
c
 120  if(mp3.or.oci.or.omcscf)
     + call caserr2('polarisability - only scf or mp2 allowed')
      call mpdder(core)
      go to 230
c
c     solve chf equations for nuclear displacements
c
 130  oprn(4) = .false.
      if (osymm) call indxsy(core)
      oprn(4) = .false.
      osymm = .false.
      call trnfkd(core)
      call chfndr(core,core)
      go to 230
c
c      calculate dipole derivatives
c
 140  if (omp2) then
         call mpdder(core)
      else if (mp3.or.oci.or.omcscf) then
         call caserr2
     +   ('correlated dipole derivatives - only mp2 allowed')
      else
         call scfdd(core)
         call anairr(core)
      end if
      go to 230
c
c     second derivatives
c
 150  continue
      if (omp2) then
         call mp2dd(core)
      else if( mp3.or.oci.or.omcscf) then
         call caserr2
     +   ('correlated second derivatives - only mp2 allowed')
      else
         call scfdd(core)
         call anairr(core)
      end if
      go to 230
c
c    magnetizability
c
 160  if(omp2.or.mp3.or.oci.or.omcscf) then
       call caserr2
     + ('no magnetic properties for correlated wavefunctions')
      endif
      call scfdd(core)
      call chfidr(core,core)
      oprn(40) = .true.
      call anairr(core)
      go to 230
c
c     intensities
c
 170  continue
      call anairr(core)
      go to 230
c
c      test  ---- user supplied routine for testing
c
 180  call test(core)
      go to 230
c
c     calculate polarisability derivatives
c
 200  if(omp2.or.mp3.or.oci.or.omcscf) then
       call caserr2
     +   ('polarisability derivatives - only scf allowed')
      endif
      call hfpder(core)
      go to 230
c
c     calculate hyperpolarisabilities
c
 210  if(omp2.or.mp3.or.oci.or.omcscf) then
       call caserr2
     +  ('hyperpolarisability - only scf allowed')
      endif
      call hfhpol(core)
      go to 230
c
 220  call hfdder(core)
      go to 230
c
c     calculate raman intensities
c
 225  if(omp2.or.mp3.or.oci.or.omcscf) then
       call caserr2('only closed-shell scf for raman intensities')
      endif
      call raman(core)
      go to 230
c
227   continue
      call mcx(core)
      go to 230
c
c      calculate infra red intensities
c
 215  if (omp2) then
         zruntp = 'hessian'
         zzzzzz = 'hessian'
         call mp2dd(core)
         zruntp = 'infrared'
         zzzzzz = 'infrared'
         iscftp = 99
         mprest = 3
         opass6 = .true.
         call mpdder(core)
         oprn(40) = .true.
         call anairr(core)
      else if (mp3.or.oci.or.omcscf) then
         call caserr2
     +   ('correlated infra-red intensities - only mp2 allowed')
      else
         call scfdd(core)
         oprn(40) = .true.
         call anairr(core)
      end if
      go to 230
c
c     NMR chemical shifts
c
 420  continue
      call gamerr(
     &'NMR functionality not supported by this build.',
     &ERR_NO_CODE, ERR_SWITCHED_OUT, ERR_SYNC, ERR_NO_SYS)
      go to 230
c
c     linear response
c
 300  continue
      call respon(core)
      goto 230
c
c     rpa+gradients
c
 310  call gamerr(
     &'RPA gradients functionality not supported by this build.',
     &ERR_NO_CODE, ERR_SWITCHED_OUT, ERR_SYNC, ERR_NO_SYS)
      goto 230
c
c     run special program options
c
c     ----- property options -----
c
 230  if (oso) then
         write(iwr,*) ' properties for spin-orbit zora are experimental'
      else
         if (omodel) then
           call model(core,core,9)
         endif
         if (zruntp.eq.zprop.and..not.opass10) then
           call hfprop(zscftp,core)
         endif
         if (oprop) then
           call hfprop2(zscftp,core)
           iprop = 0
         endif
      end if
cdrf-----------------------------------------------------
      if (oreact) then
       if (field.ne.'        ') then
        ist=iactst
c
c      -----  analyse rf contributions
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
     5 uelstg(1,ist),suintg(1,ist),uintg(1,ist),repmodg(1,ist),core)
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
     5 core)
c
        endif
c
        call drfout
       endif
      if (idpa .gt. 0) then
        call dppop(core)
c       call dppope(core)
      endif
      endif
cdrf----------------------------------------------------
c
c     --- formatted file for graphics ---
c
c     this call deals with results that have been stored on the
c     dumpfile
c
      call blk(core)

c
c     ---  Bader interface
c
      call aimpac(core)


      irun = locatc(zrunt,mxtask,zruntp)
      if (irun.eq.0) go to 260
      go to (260,250,240,260,260,260,260,260,250,260,250,250,260,260,
     +       260,250,250,260,260,260,260,260,260,260,260,260,260,260,
     +       260,260,260,260,260,260,260,260,260,260,260,260,260,250,
     +       260,260,260,260,260,260,260,260) , irun
 240  if (.not.ofulci) call closbf2(-1)
 250  if (.not.odirpa) call closbf(-1)
      go to 270
 260  if (zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
     &    .or. zscftp.eq.zcas(3)) go to 250
 270  oresti = .true.
      itask(mtask) = irest
      call revise
c
c  flush incore files to disk (GA ed3 only at present)
      call ioupdate
c
      go to 20
 280  continue
      call gmem_null_finalise(core)
      return
      end
      subroutine pacini
c
c    initialise values of variables
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
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      logical jaco
      common/jacobc/jaco
c
c      common/funcon/bohr(8)
      real*8 bohr,charge,veloc,hbar,pi,amu,avo,debye,hartre
      common/funcon/bohr,charge,veloc,hbar,pi,amu,avo,debye,hartre
c
c
      common/jspec/nspj(maxorb+1)
c
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      character *8 qnames
      common/qercmx/qnames(50)
      logical oqskip,omixed
      common/qercom/frqq(31),nfrqq,npoleb,iq6(36),
     + iqsec(50),oqskip(50),iqang(50),npb,ijunq(4),omixed,ifileb,iblkb
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/enuc,etotal,ehfock,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              mcyc,ischm,lock,maxit,nconv,npunch,lokcyc
      common/gesblk/alphs(maxorb),iswop(2,100),iextra(100),nswop,
     1              nextra,ifilva,iblka1,iblka2,ifilvb,iblkb1,
     2              iblkb2
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
c
      logical lmos,lnucl
      common/prpsec/isecd1,itypd1,isecd2,itypd2,lmos,isecv,
     1  itypv,isecd3,itypd3,lnucl
c
      real*8 clat, zlat
      integer nlat, maxlat
      parameter (maxlat=216)
      common /lattic/ clat(3,maxlat),zlat(maxlat),nlat
c
c
      common/junke/maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
c
      common/crnams/ iangc(40),iangs(24)
      character *8 pnamc,ac6,pnams
      character *4 bfnam,atnam
      common/crnamx/ac6(6),pnamc(40),pnams(24),bfnam(20),atnam(37)
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common /picon/ pito52,pidiv4,root3,root5,root53,root7
c
c
      character *8 hcore
c
      data zero,one,pt5/0.0d0,1.0d0,0.5d0/
      data hcore /'hcore'/
c
      call paczer
      root3 = dsqrt(3.0d0)
      root5 = dsqrt(5.0d0)
      root53= root5/root3
      root7 = dsqrt(7.0d0)
      mppas = 1
      mcflag = 0
      nmp = 0
      do 20 i = 1 , ntypr
         opunch(i) = .false.
         odebug(i) = .false.
         oprn(i) = .false.
 20   continue
      oprn(1) = .true.
      oprn(2) = .true.
      oprn(4) = .true.
      oprn(5) = .true.
      oprn(6) = .true.
      oprn(7) = .true.
      oprn(10) = .true.
c
      jaco = .true.
      pi     = dacos(-1.0d0)
      bohr   = toang(1)
      charge = toang(3)              ! was: 4.803250d-10
      veloc  = toang(9)              !      2.997924562d10
      hbar   = toang(4)/(2*pi)*1.0d7 !      1.0545887d-27
      amu    = toang(2)*1.0d3        !      1.660531d-24
      avo    = toang(5)              !      6.022169d23
      debye = charge*bohr*1.0d10
      hartre = toang(8)
c
c      default /prpsec/
c
      isecd1 = 207
      isecd2 = 210
      isecd3 = 241
      isecv = 208
      lmos = .false.
      lnucl = .true.
      itypd1 = 0
      itypd2 = 0
      itypd3 = 0
      itypv = 0
c
c       zero restart blocks
c
      do 30 i = 1 , 12
         ione(i) = 1
 30   continue
c
c     /cigrad/ block
c
      isecdd = 301
      isecll = 302
      iscigr = 303
      isecmo = 307
      isecnd = 305
      isecsy = 306
      ifil2d = n11tap(1)
      iblk2d = n11bl(1)
      cigr = .false.
      mpgr = .false.
      mcgr = .false.
      cicx = .false.
      cicv = .false.
      umpgr = .false.
c
      do 40 i = 1 , 20
         mpstrm(i) = i + 20
 40   continue
      ncore = 0
      nvr = 0
      do 50 ix = 1 , maxorb
         nspj(ix) = 0
 50   continue
c
c     defaults in options block
c
c     itol=18
      itoli = 18
      nruns = 1
c     icut=10
c     nopk=1
      icuti = 10
c     irest=-1
      iresti = 0
      norder = 2
c     nprint=-1
c     ist=1
c     jst=1
c     kst=1
c     lst=1
c     intloc=1
c     nrec=1
c     mul=1
      nmol = 1
      nmul(1) = 1
      nmul(2) = 1
      iochf(1) = n12bl(1)
      iblk(1) = n2blk(1)
      lblk(1) = n2last(1)
c     iblkd=1
c     iblks=1
      iblkf = n12bl(1)
      ifockf = n12tap(1)
c
      ifill = 1
      iblkl = 1
c
      notape(1) = n2tape(1)
c     ifild=4
c     ifils=8
c
      jjfile = n2file
      nnfile = n4file
      mmfile = n6file
      nofile(1) = n4tape(1)
      nufile(1) = n6tape(1)
      jblk(1) = n4blk(1)
      kblk(1) = n6blk(1)
c
      irfile = jjfile
      irunit = notape(1)
      irblok = iblk(1)
      irbl = irblok - lblk(1) - 1
c
      ngpts = 999999
      nps1 = 1
      nps2 = 1
      guess = hcore
c     nopt=4
      iorder = 2
c     nvib=1
c     vibsiz=0.001e0
      ropt = 0.01d0
      scftyp = ' '
      runtyp = ' '
      llibry = .true.
      minvec = 1
c
      bfgs = .true.
c
      nswop = 0
      nextra = 0
      nlat = 0
      ifilva = 4
      iblka1 = 0
      iblka2 = 0
      ifilvb = 0
      iblkb1 = 0
      iblkb2 = 0
c
c     scf block default options
c
      lock = 0
      lokcyc = 999
      nconv = 5
      npunch = 0
      maxit = 50
      mcyc = 5
      ischm = 3
      sh1(1) = pt5
      sh1(2) = pt5
      sh2(1) = 0.2d0
      sh2(2) = 0.2d0
      ldam12 = .false.
      ldam23 = .false.
      ldam13 = .false.
      d13 = one
      d23 = one
      d12 = one
      lfield = .false.
      fixed = .false.
c     iscftp = 0
      gap1(1) = one
      gap1(2) = one
      gap2(1) = pt5
      gap2(2) = pt5
      canna = 2.0d0
      cannb = 0.0d0
      cannc = -2.0d0
      fx = zero
      fy = zero
      fz = zero
      ldiis = .true.
      lex = .true.
c
c     property/polarisability defaults
c
      ipol = 0
      ldens = .true.
      ldiag = .false.
      npstar = 0
      iconvv = 8
      npole = 0
      nfreq = 0
      np = 0
      npa = 0
      do 60 i = 1 , 9
         ipang(i) = iangc(i)
         pnames(i) = pnamc(i)
 60   continue
      do 70 i = 1 , 50
         ipsec(i) = 0
         opskip(i) = .true.
 70   continue
      ogen = .false.
      do 80 i = 1 , 9
         if (ione(i+3).ne.0) then
            np = np + 1
            ipsec(i) = isect(i+21)
         end if
 80   continue
      ospher = .false.
      omixed = .false.
      npa = np
      do 90 i = 1 , 30
         oc6(i) = .false.
 90   continue
      do 100 i = 6 , 10
         oc6(i) = .true.
 100  continue
      nc6min = 6
      nc6max = 10
c
      npass1 = 1
      npass2 = 1
      iacc = 10
      iontrn = 0
c
c     /lattic/ block
c
      nlat = 0
c
      return
      end
      subroutine paczer
c
c  initialise
c
      implicit real*8  (a-h,o-z)
      logical logg
      character *8 ztit
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
      common/restrr/realg(30)
      common/restrl/logg(50)
      common/restri/ispacg(2095)
      common/restrz/ztit(50)
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
      character *8 blank
      data blank/'        '/
      do 20 i = 1 , 50
         ztit(i) = blank
 20   continue
      do 30 i = 1 , 30
         realg(i) = 0.0d0
 30   continue
      do 40 i = 1 , 50
         logg(i) = .false.
 40   continue
      call setsto(18,0,irest2)
      call setsto(len_cndx41,0,ncoorb)
      call setsto(63,0,ispacg)
      return
      end
      subroutine inittxt
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      character *4 ylabel
      character *8 ztype, zlist, znumb, zhead
      character *8 zheadr, zunit, ztitle, ztitl, zspin
      character *8 zsta, zhz, zgauss, zcont, zconta, zmul
      common/junkc/ylabel(26),ztype(35),zlist(300),znumb(100),zhead(90),
     * zheadr(20),zunit(33),ztitle(100),ztitl(120),zspin(30),
     * zsta,zhz,zgauss,zcont(10),zconta(maxat),zmul(2,maxorb)
c
      character *8 zzsta,zzhz,zzgauss
      character *8 zzhead, zzheadr, zzunit, zztitle, zztitl
      character *8 zzspin, zztype
      character *4 yylabel
      dimension zzhead(90),zzheadr(20),zzunit(33),zztype(35)
      dimension zztitle(100),zztitl(120),zzspin(30)
      dimension yylabel(26)
c
      data zzhead/'      ','    po','tentia','l     ','      ',
     $          '      ','magnet','ic shi','elding','      ',
     $          '      ','  elec','tric f','ield  ','      ',
     $          '   ele','ctric ','field ','gradie','nts   ',
     $          '      ','  dipo','le mom','ents  ','      ',
     $          '      ','quadru','pole m','oments','      ',
     $          ' diama','gnetic',' susce','ptibil','ities ',
     $          '      ','  seco','nd mom','ents  ','      ',
     $          '      ',' octup','ole mo','ments ','      ',
     $          '    th','ird mo','ments ','(part ','a)    ',
     $          '    th','ird mo','ments ','(part ','b)    ',
     $          '     h','exadec','apole ','moment','s     ',
     $          '    fo','urth m','oments',' (even',')     ',
     $          '     f','ourth ','moment','s (odd',')     ',
     $          '      ','     o','verlap','      ','      ',
     $          '     p','lanar ','densit','y     ','      ',
     $          '      ','line d','ensity','      ','      ',
     $          '     c','harge ','densit','y     ','      '/
      data zzheadr/
     *'  isot','ropic ','coupli','ng con','stants',
     *'anisot','ropic ','coupli','ng con','stants',
     *'     a','ngular',' momen','tum (l',')     ',
     *'      ','     l','/(r**3',')     ','      '/
      data zztitle/' 1/r   ',9*'      ',
     $           'sg(xx)','sg(yy)','sg(zz)','sg(xy)','sg(xz)','sg(yz)',
     1           ' 1/r  ','xx/rrr','yy/rrr','zz/rrr',
     $           '  ex  ','  ey  ','  ez  ',7*'      ',
     $           'fg(xx)','fg(yy)','fg(zz)','fg(xy)','fg(xz)','fg(yz)',
     1           4*'      ',
     $           '  x   ','  y   ','  z   ',7*'      ',
     $           'qm(xx)','qm(yy)','qm(zz)','qm(xy)','qm(xz)','qm(yz)',
     1           ' rsqd ',3*'      ',
     $           'ch(xx)','ch(yy)','ch(zz)','ch(xy)','ch(xz)','ch(yz)',
     1           ' rsqd ',3*'      ',
     $           '  xx  ','  yy  ','  zz  ','  xy  ','  xz  ','  yz  ',
     1           '  rr  ',3*'      ',
     $           'o(xxx)','o(xyy)','o(xzz)','o(yxx)','o(yyy)','o(yzz)',
     1           'o(zxx)','o(zyy)','o(zzz)','o(xyz)',
     $           ' xxx  ',' xyy  ',' xzz  ',' yxx  ',' yyy  ',' yzz  ',
     1           ' zxx  ',' zyy  ',' zzz  ',' xyz  '/
      data zztitl/ ' xrr  ',' yrr  ',' zrr  ',7*'      ',
     $           'hxxxx ','hyyyy ','hzzzz ',' xxrr ',' yyrr ',' zzrr ',
     1           ' rrrr ',' xyrr ',' xzrr ',' yzrr ',
     $           ' xxxx ',' yyyy ',' zzzz ',' xxyy ',' xxzz ',' yyzz ',
     1           ' xxrr ',' yyrr ',' zzrr ',' rrrr ',
     $           ' xyxx ',' xyyy ',' xyzz ',' xzxx ',' xzyy ',' xzzz ',
     1           ' yzxx ',' yzyy ',' yzzz ','      ',
     $           '  1   ',9*'      ',
     $           ' dxy  ',' dxz  ',' dyz  ',7*'      ',
     $           '  dx  ','  dy  ','  dz  ',7*'      ',
     $           '  df  ',9*'      ',
     *'  a   ',9*'      ',
     *' t(xx)',' t(yy)',' t(zz)',' t(xy)',' t(xz)',' t(yz)',
     *4*'      ',
     *'  lx  ','  ly  ','  lz  ',7*'      ',
     *'lx/rrr','ly/rrr','lz/rrr',7*'      '/
      data zzunit/
     *'esu/cm',' ',' ',
     *'ppm',' ',' ',
     *'10**8 dy','n/esu',' ',
     *'10**16 e','su/cm**3',' ',
     *'debye',' ' ,' ' ,
     *'10**-26 ','esu.cm**',' ',
     *'10**-6  ','emu/mol',' ',
     *'10**-16 ','cm**2',' ',
     *'10**-34 ','esu.cm**','3',
     *'10**-24 ','cm**3',' ',
     *'10**-24 ','cm**3',' ' /
      data zzsta,zzhz,zzgauss/'****','mhz','gauss'/
      data zzspin/
     *'h1  ','he3 ','li7 ','be9 ','b11 ','c13 ','n14 ',
     *'o17 ','f19 ','ne21','na23','mg25','al27','si29',
     *'p31 ','s33 ','cl35','ar39','k39 ','ca43','sc45',
     *'ti47','v51 ','cr53','mn55','fe57','co59','ni61',
     *'cu63','zn67'/
      data zztype/'s  ',
     +        'x','y','z',
     +        'xx','yy','zz','xy','xz','yz',
     +        'xxx','yyy','zzz','xxy','xxz','xyy','yyz','xzz',
     +        'yzz','xyz',
     +        'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz',
     +        'zzzx','zzzy','xxyy','xxzz','yyzz','xxyz','yyxz',
     +        'zzxy' /
      data yylabel /'1s','2s','2p','2sp','3s','3p','3sp',
     +     '3d','4s','4p','4sp','4d','5s','5p','5d',
     +     's','p','d','f','g','k','l','m','5sp','sp','****'/
c
      do i = 1,26
       ylabel(i) = yylabel(i)
       zheadr(i) = zzheadr(i)
      enddo
c
      do i = 1,35
       ztype(i) = zztype(i)
      enddo
c
      do i = 1,90
       zhead(i) = zzhead(i)
      enddo
c
      do i = 1,33
       zunit(i) = zzunit(i)
      enddo
c
      do i = 1,100
       ztitle(i) = zztitle(i)
      enddo
      do i = 1,120
       ztitl (i) = zztitl (i)
      enddo
c
      do i = 1, 30
       zspin(i) = zzspin(i)
      enddo
c
      zsta = zzsta
      zhz  = zzhz
      zgauss = zzgauss
c
      return
      end
      subroutine scfgvb(core)
      implicit none
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
c...  polarisation-disp common
c
      integer iunit_a,iblk_a,isec_a,iunit_b,iblk_b,isec_b
      integer nocc_a,nmo_a,nocc_b,nmo_b
      logical odisp
      common/disp/ odisp,iunit_a,iblk_a,isec_a,iunit_b,iblk_b,isec_b,
     1             nocc_a,nmo_a,nocc_b,nmo_b
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
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
c
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
      integer ncore,nact,nels
      common /fccwfn2 / ncore,nact,nels
c
      dimension core(*)
c
      real*8 dum, etotal, core
      integer idum,ierror
      logical opass

      idum = 0
      dum = 0.0d0
      opass = opass5

      if (mp2 .or. mp3) then
         if(odscf)then
c
c direct mp2 energy (murray code)
c
            if (CD_active()) then
               ierror = CD_jfit_init1()
               if (ierror.ne.0) then
                  write(iwr,600)ierror
                  call caserr2('Out of memory in incore coulomb fit')
               endif
            endif
            call hfscf(core)
            call drcdrv(core,'mp2')
            if (CD_active()) then
               ierror = CD_jfit_clean1()
               if (ierror.ne.0) then
                  call caserr2(
     +                 'Memory problem detected in CD_jfit_clean1')
               endif
            endif
         else
            call emp23(core,etotal)
         endif
      else if (omas) then
        if(odscf)then
c
c direct masscf energy
c
          call masscf(core,etotal)
         else
          call masscf(core,etotal)
         endif
      else
         if (CD_active()) then
            ierror = CD_jfit_init1()
            if (ierror.ne.0) then
               write(iwr,600)ierror
               call caserr2('Out of memory in incore coulomb fit')
            endif
         endif
 600     format('*** Need ',i10,' more words to store fitting ',
     +          'coefficients and Schwarz tables in core')
         call hfscf(core)
         if (CD_active()) then
            ierror = CD_jfit_clean1()
            if (ierror.ne.0) then
               call caserr2('Memory problem detected in CD_jfit_clean1')
            endif
         endif
      endif

      call ioupdate

      itask(mtask) = irest
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
      if (.not.(tim.ge.timlim .or. opass.or.odisp)) then

cgdf  07.04.08  1st, output den1, then fix hfprop
          if (omas) then
             write(iwr,*)' '
             write(iwr,*)'(masscf properties unavailable)'
             write(iwr,*)' '
             return
          end if

         if (.not.opass10) call hfprop(zscftp,core)
      end if
      return
      end
      subroutine scfrun(q)
c
c      driving routine for scf
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/enuc,etotal,ehfock,sh1(2),sh2(2),gap1(2),gap2(2),
     &              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     &              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     &              mcyc,ischm,lockc,maxit,nconv,npunch,lokcyc
      common/blkin/potnuc(10)
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
      dimension q(*)
      data m16,m10/16,10/
c
c
      if (.not.opass5) then
       if (irest.le.4) then
        call scf(q)
        if (odscf. and. guess.eq.'atoms') go to 1000
       endif
      end if
      call putdev(q,mouta,7,1)
      if (scftyp.eq.'open') call putdev(q,moutb,10,2)
c
      call secget(isect(494),m16,iblk34)
      call rdedx(potnuc,m10,iblk34,ifild)
      enuc = potnuc(1)
      ehfock = potnuc(2)
      etotal = potnuc(3)
c
      call secget(isect(13),13,iblok)
      call wrt3(enuc,lds(isect(13)),iblok,ifild)
1000  opass1 = .false.
      opass2 = .false.
      opass3 = .false.
      opass4 = .false.
      opass5 = .false.
      return
      end

      subroutine start(core,oresti,ostop,ostopr)
c
c     read in basis set + options
c
c     irest = 0     normal start + normal running condition.
c     irest = 1     2e-integral restart ( 1e +mo*s saved)
c     irest = 2     integral symmetry adaption restart
c     irest = 3     scf restart ( 1e + mo*s saved; 2e saved)
c     irest = 4     2-particle density transformation restart
c     irest = 5     1e-gradient restart ( mo*s saved; no gradient saved)
c     irest = 6     2e-gradient restart ( mo*s saved; gradient saved)
c     irest = 7     1e-properties restart
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      integer idpa
      common/dppopc/idpa

c
      dimension core(*)
      dimension ztestl(2),zstats(18)
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
      character*2 zelem(0:103)
      integer nelem
      parameter(nelem=103)
      common/periodic/zelem
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      real*8 fieldx, fieldy, fieldz
      logical ofield
      common /gms_field/  fieldx, fieldy, fieldz, ofield
c
      common/jorgs/maxj,irecal,opowel,obfgs,obfgx,onrfo,ocut,
     *outjor,outdeb,maxsta
c
c
      real*8 tchg
      integer ifree, kform
      logical opress,otemp,omapp,omodel,ocharg,oaminc,oamhes,osubs
      common /modj/ kform,opress,otemp,omapp,omodel,ocharg,
     +              oaminc,oamhes,tchg(maxat),osubs,ifree
c
c
      real*8 bbeta, box
      integer ntx, ntxo, nrc, nrcx, npm, nrp, nsm, nram, ntb
      integer nspaca
      common /moda/ ntx,ntxo,nrc,nrcx,npm,nrp,nsm,nram,ntb,nspaca,
     +  bbeta,box(3)
c
      common/memla/i02,i04,i06,i08,i10,i12,i14,i16,i18,i20,
     +             i22,i24,i26,i28,i30,i32,i34,i36,i38,i40,
     +             i42,i44,i46,i48,i50,i52,i54,i56,i58,i60,
     +             i62,i64,i66,i68,i70,i72,i74,i76,i78,i80
      common/memlb/l05,l10,l15,l20,l25,l30,l35,l40,l45,l50,
     +             l55,l60,l65,l70,l75,l80,l85,l90,l95,l100,
     +             k02,k04,k06,k08,k10,k12,k14,k16,k18,k20,
     +             k22,k24,k26,k28,k30,k32,k34,k36,k38,k40
      common/g80nb/natmod,nnbg,mapmod(maxat),ipnonb(maxat+1),
     +             inonb(maxat*6)
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
      real*8 func0
      integer ncoord, npts, nserch, iupdat, icode, iupcod
      common /seerch/ func0,ncoord,npts,nserch,iupdat,icode,iupcod
c
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/avstat/nwt,ispr(3),weight(5)
c
      integer n20cas
      integer ifsimu, isimul, iam, noci, ipople
      parameter (n20cas = 100)
      common /simul/ ifsimu,isimul(n20cas),iam(n20cas),noci(n20cas),
     +               ipople
c
c
      integer mode, nrever, ifhes, isort, iaugm, isp, ibuffn
      common /ctrl/ mode,nrever,ifhes,isort(n20cas),iaugm,
     +              isp,ibuffn
c
      common/testop/oto(20)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/blksiz2/ns2
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
c
c
      real*8 fudgit, fmax, swnr, swsimu
      integer nhesnr
      logical ojust, osuped
      common /ciconv/ fudgit(10),fmax,ojust,osuped,swnr,swsimu,nhesnr
c
      common/degen /ifdeg,idegen(100),ophase
      common/exc/kkk,kkk1,czl(5),idom(5),iocca(5,15),ioccb(5,15),irestu
     * ,ifstr
      common/cipar /clevc,accci,timeci,nitci,nfudge,ncfl,ncfu,nprinc
      common/blbpar/accblb,clevs,timblb,nitblb(2)
c
      real*8 cccnv
      logical oconv, oconvf
      integer icanon
c
      common /gms_finish/ cccnv,oconv,icanon,oconvf
C       common/prints/oprint(20),
C      * odist, obass, oadap, osadd, ovect, oanal, ohess,
C      * onopr(13),
C      * opmull, oplowd, opmos , opsadd, opvect, opadap, opdiis,
C      * opsymm, opfock, opq   , opao  , opmo  , opdump, opgf  ,
C      * optda , opdist, opdire, oprin(3),
C      * ciprnt, oecho
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
c
      common/l     /llock
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
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
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      common/restrr/
     + gxx,gyy,gzz,rspace(21),tiny,tits(2),scalec,roptc,roptcc
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
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
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
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
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
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
      real*8 alphas,chargat,spinat
      integer nswapa, nswapb, nswap, next,isecat
      integer natconf,nnatconf,iatcon,iatstates,natdiff
      parameter (nnatconf=20)
      logical odiff,oexc,oground,oatdo,oalway,uhfatom,nonrel,forcat
      character*8 zatconf,atmode,zatdiff
      character*60 string_ch
      character*2 zatstate
c
      common /datgue/ alphas(maxorb),nswapa,nswapb,nswap(80),
     +                next(maxorb),oground,oatdo,isecat,oalway,
     +                odiff,oexc
      common /datgue2/ chargat(nnatconf),spinat(2,4,nnatconf),
     +                 natconf,uhfatom,iatcon,nonrel(nnatconf),
     +                 forcat(nnatconf),iatstates(nnatconf),
     +                 natdiff
      common /datgue3/ string_ch(2,nnatconf),zatconf(nnatconf),atmode,
     +                 zatstate(nnatconf),zatdiff(nnatconf)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nnacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
      common/junkc/ylab(26),ztype(35),zlist(100),zcentg(100),
     * ztagc(100),zhead(496),ztita(10),zatagg(maxat),zmul(maxorb)
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      character * 132 char1,char2,char1c
      common/workc/char1,char2,char1c
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
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common /copt / rcvopt
c
      real*8 tr, trx, try, trz, prmoms, aprev
      integer indmx, jaxis, igroup, igrp80
      logical oprmom
      common /molsym/ tr(3,3),trx,try,trz,indmx,jaxis,igroup,igrp80,
     +                prmoms(3),aprev(maxat3,3),oprmom
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
c The oblock array determines what is written out to the punchfile:
c
c Below sourced from server.m
c oblock(1)   - coords - single point or last point
c keyword: coor
      integer PUN_COORD
      parameter (PUN_COORD=1)
c oblock(2)   - coordinates at each optimisation step
c keyword: opti
      integer PUN_OPT_COORD
      parameter (PUN_OPT_COORD=2)
c oblock(3)   - starting coordinates
c keyword: init
      integer PUN_START_COORD
      parameter (PUN_START_COORD=3)
c oblock(4)   - atom connectivity
c keyword: conn
      integer PUN_ATOM_CONN
      parameter (PUN_ATOM_CONN=4)
c oblock(5)   - molecular orbital occupations
c keyword: occu
      integer PUN_MO_OCC
      parameter (PUN_MO_OCC=5)
c oblock(6)   - eigenvectors
c keyword: vect
      integer PUN_EIGV
      parameter (PUN_EIGV=6)
c oblock(7)   - run type
c keyword: type
      integer PUN_RUNT
      parameter (PUN_RUNT=7)
c oblock(8)   - gvb pair data
c keyword: gvb
      integer PUN_GVBPAIR
      parameter (PUN_GVBPAIR=8)
c oblock(9)   - force constant matrix
c keyword: forc
      integer PUN_FCMAT
      parameter (PUN_FCMAT=9)
c oblock(10)  - vibrational frequencies
c keyword: vibr
      integer PUN_VIBFREQ
      parameter (PUN_VIBFREQ=10)
c oblock(11)  - basis set
c keyword: basi
      integer PUN_BASIS
      parameter (PUN_BASIS=11)
c oblock(12)  - orbital energies
c keyword: eige
      integer PUN_ORB_ENER
      parameter (PUN_ORB_ENER=12)
c oblock(13)  - total scf energies
c keyword: scfe
      integer PUN_SCF_ENER
      parameter (PUN_SCF_ENER=13)
c oblock(14)  - job title
c keyword: titl
      integer PUN_JOB_TITLE
      parameter (PUN_JOB_TITLE=14)
c oblock(15)  - overlap matrix
c keyword: over
      integer PUN_OVLP_MAT
      parameter (PUN_OVLP_MAT=15)
c oblock(16)  - grid data from dumpfile section(s)...
c keyword: grid
      integer PUN_GRID_DATA
      parameter (PUN_GRID_DATA=16)
c oblock(17)  - mulliken analysis
c keyword: mull
      integer PUN_MULLIKEN
      parameter (PUN_MULLIKEN=17)
c oblock(18)  - lowdin analysis
c keyword: lowd
      integer PUN_LOWDIN
      parameter (PUN_LOWDIN=18)
c oblock(19)  - spin densities
c keyword: spin
      integer PUN_SPIN_DENS
      parameter (PUN_SPIN_DENS=19)
c oblock(20)  - scf type
c keyword: leve
      integer PUN_SCF_TYPE
      parameter (PUN_SCF_TYPE=20)
c oblock(21)  - normal coordinates for vibrational modes
c keyword: norm
      integer PUN_VIB_MODES
      parameter (PUN_VIB_MODES=21)
c oblock(22)  - two electron integrals
c keyword: twoe
      integer PUN_2E_INT
      parameter (PUN_2E_INT=22)
c oblock(23)  - gradients
c keyword: grad
      integer PUN_GRADIENTS
      parameter (PUN_GRADIENTS=23)
c oblock(24)  - transformation matrix
c keyword: tran
      integer PUN_TRANS_MAT
      parameter (PUN_TRANS_MAT=24)
c oblock(25)  - potential derived charges
c keyword: pdc
      integer PUN_POT_DERV_CHG
      parameter (PUN_POT_DERV_CHG=25)
c oblock(26)  - grid symmetry array
c keyword: gsym
      integer PUN_GRID_SYMM
      parameter (PUN_GRID_SYMM=26)
c oblock(27)  - cartesian hessian matrix
c keyword: secd
      integer PUN_HESSIAN
      parameter (PUN_HESSIAN=27)
c oblock(28)  - distributed multipole analysis
c keyword: dma
      integer PUN_MPOLE_ANAL
      parameter (PUN_MPOLE_ANAL=28)
c oblock(29)  - density matrix
c keyword: dens
      integer PUN_DENS_MAT
      parameter (PUN_DENS_MAT=29)
c oblock(30)  - total energy
c keyword: ener
      integer PUN_TOT_ENER
      parameter (PUN_TOT_ENER=30)
c oblock(31)  - internal coordinate hessian matrix
c keyword: hess
      integer PUN_ZMAT_HESS
      parameter (PUN_ZMAT_HESS=31)
c oblock(32)  - dipole moment
c keyword: dipo
      integer PUN_DIPOLE
      parameter (PUN_DIPOLE=32)
c oblock(33)  - infrared intensities 
c keyword: infr
      integer PUN_INFRARED
      parameter (PUN_INFRARED=33)
c oblock(34)  - raman intensities
c keyword: rama
      integer PUN_RAMAN
      parameter (PUN_RAMAN=34)
c CURRENTLY ONLY USED FOR XML
c oblock(35)  - metadata (user, compilation date etc)
c keyword: meta
      integer PUN_METADATA
      parameter (PUN_METADATA=35)
c CURRENTLY ONLY USED FOR XML
c oblock(36)  - input parameters
c keyword: inpa
      integer PUN_INPUT_PARAM
      parameter (PUN_INPUT_PARAM=36)


      integer LENPUN
      parameter(LENPUN=60)

      logical oblock, opunold
      integer nfblck, nblsec, iblsec, nblgri, iblgri
      integer itwoe, ntwoe, iblfmt
      common/blocks/nfblck,nblsec,iblsec(10),nblgri,iblgri(10),
     +              itwoe,ntwoe,iblfmt(LENPUN),oblock(LENPUN),opunold
c

c
      real*8  vg, gxp, gyp, gzp, d2, valmo, occa, occb
      real*8  rij, timat, denat, encoat, eltot
      integer idenfu,nangpr, ncorxx, ndeg, nomo
      logical odenfu, osymm, osc
      common /dnfnw/ odenfu, osymm, osc, idenfu, nangpr, ncorxx,
     +               vg(200),gxp(200),gyp(200),gzp(200),d2(200),
     +               valmo(200),occa(200),occb(200),
     +               rij(1275),timat(50),denat(50),encoat(50),ndeg(50),
     +               eltot,nomo
c
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
c
      common/dnfdbg/denptx,denpty,denptz,odenpt
c
c      following common block for rfcontinuum scf calcs
c      contains input data
c
      real*8 aradix,aradiy,aradiz, dielec
      logical oscrf
      common/scrf/aradix,aradiy,aradiz,dielec,oscrf
c
c ----- common blocks for psscrf calculation ----
c
      real*8 ptspac, delect, facnrm
      integer iconv,inkblk,itotbq,iscf,ispace,msecp
      logical opssc
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
      common/surfgn/numbq(maxat),irejec(maxat),cavsze(maxat)
     * ,nocent,itmax,damper
c
c ----- common blocks for xfield correction term
c
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c
      character *8  all,filb
      logical lres
c
c
c user defined orientation
c
      common/oldtr/trold(3,3),trano(3),igpold,naxold,osetor
c
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      logical okillv
      integer isecvv, itypvv
      common /vectrn/ isecvv,itypvv,okillv
c
      character *8 groupy
      common/indsyx/groupy
      common/smaljr/ncorea(8),nacta(8),nvia(8),numsym(8)
      logical jaco
      common/jacobc/jaco
c
c
      character *8 pnames
      common/tdhfx/pnames(50)
c
c
      real*8 freq, w0
      integer nfreq, npole, ic6, nc6, iblkhi, ifreq, nc6min, nc6max
      integer ipsec, ipang, npa, ispa
      logical oc6, opskip, ogen, ospher
      common /tdhf/ freq(30),w0,nfreq,npole,ic6,nc6,iblkhi,ifreq,
     +              oc6(30),nc6min,nc6max,ipsec(50),opskip(50),
     +              ipang(50),npa,ogen,ospher,ispa
c
c
      character *8 qnames
      common/qercmx/qnames(50)
      logical oqskip,omixed
      common/qercom/frqq(31),nfrqq,npoleb,iq6(36),
     + iqsec(50),oqskip(50),iqang(50),npb,ijunq(4),omixed,ifileb,iblkb
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
      real*8 fjk, erga, ergb, cana, canb, damgen, shfgen, fcan
      integer nact, iactiv, nbshel, ilfshl, njk, njk1, nspace
      common /ghfblk/ nact,iactiv(maxorb),nbshel(11),ilfshl(11),
     + njk,njk1,nspace,fjk(11),erga(121),ergb(121),cana(121),canb(121),
     + damgen(121),shfgen(121),fcan(11)
c
c
      logical lmos,lnucl
      common/prpsec/isecd1,itypd1,isecd2,itypd2,lmos,isecv,
     1  itypv,isecd3,itypd3,lnucl
c
      logical lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis
      common/scfblk/enuc,etotal,ehfock,sh1(2),sh2(2),gap1(2),gap2(2),
     1              d12,d13,d23,canna,cannb,cannc,fx,fy,fz,
     2              lfield,fixed,lex,ldam12,ldam13,ldam23,ldiis,
     3              mcyc,ischm,lockc,maxit,nconvc,mpunch,lokcyc
c
      character *8 groupx
      common/symtrx/groupx
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      character *8 grp
      dimension grp(8)
      character *8 zpunch
      dimension zpunch(ntypr)
      character *4 prn
      dimension prn(ntypr)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      common/grest2/lres,if2v,ib2v1,ib2v2,irdum(3)
c
      common/small/iint(maxorb)
      common/jspec/nspj(maxorb+1)
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
      character *8 buff,buff1
      character *8 words
      character *8 grhf
c
      character *8 pnamc,ac6,pnams
      character *4 bfnam,atnam
      common/crnamx/ac6(6),pnamc(40),pnams(24),bfnam(20),atnam(37)
c
c...  counterpoise common
c
      character *8 ztag00
      integer maxghost,nghost
      parameter (maxghost=100)
      common/ghost/ ztag00(maxghost)
      common/ghostn/ nghost
c
c...  polarisation-disp common
c
      integer iunit_a,iblk_a,isec_a,iunit_b,iblk_b,isec_b
      integer nocc_a,nmo_a,nocc_b,nmo_b
      logical odisp
      common/disp/ odisp,iunit_a,iblk_a,isec_a,iunit_b,iblk_b,isec_b,
     1             nocc_a,nmo_a,nocc_b,nmo_b
c
c****** gen contains most of the stuff about the direct-scf options
      common/gen/maxitd,iterd,concrt,convd,energd,thrshd,dcmax
     +,olog1,olog2,olog5,olog9,olog6,odbg(10),omp2d,mp2fo,mp2fv
     +,oreadm,owritm,orestd,iaccur,ocvary,ogrdt,ogopt,ogdma
     +,ogdold
     +,ogscf,noptd,npoint,timeg,ogmull,iswapd,dshift,irhess
      common/gen2/ripftd(21,2),rijkld(21,2),rjpftd(21,2),
     + ripfal(21,2),rijkal(21,2),
     + prefac(21),tinner(21),rifall,riffew,ricfb1,ricfb2
      logical ofull,odiisn
      common/gen3/ofull(100),odiisn
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
c
c symmetry assignment
c
      common/atmol3o/osim1,osim2
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
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
c
      real*8 sco,pco,dco,fco,gco
      integer ilopri,ngroup,iatom
      integer mng,mpergr
      parameter (mng=10,mpergr=30)
      common /groups/ ilopri,ngroup,iatom(mpergr,mng),
     +                sco(mng), pco(mng), dco(mng), fco(mng),
     +                gco(mng)
c
c
      real*8 cutomp
      common /mpcut/ cutomp
c
c
      character *8 stop
      character *8 yes,on
c
c ---- for iwhere array (mfile memory etc)
c
      integer isel,iselr ,iselw,irep,ichek,ipos,
     *  nam,keep,istat,
     *  len,iop,iwhere,isync,
     *  iposf,keepf,
     *  istatf,lrecf,nrecf
      logical oform
      logical oedpr, oftpr, ogafsmode
      integer ioupd
      logical oadd
      common/disc/isel,iselr ,iselw,irep,ichek,ipos(maxlfn),
     *  nam(maxlfn),keep(maxlfn),istat(maxlfn),
     *  len(maxlfn),iop(maxlfn),iwhere(maxlfn),isync(maxlfn),
     *  iposf(maxfrt),keepf(maxfrt),
     *  istatf(maxfrt),lrecf(maxfrt),nrecf(maxfrt)
     * ,oform(maxfrt)
     *  ,oedpr(maxlfn),oftpr(maxfrt), ogafsmode,
     *  ioupd,oadd(maxlfn)
c
c ---- for morokuma analysis
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

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
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
       logical odft
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
c.... for harmonic
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c...   common for overruling local criteria (e.g. dependency, diag thresholds)
      integer  i_depen_check, i_global_diag, n_depen_check
      logical  o_depen_print
      common/overrule/ i_depen_check,i_global_diag, o_depen_print, 
     1                 n_depen_check
c
c     i_depen_check : criterion for throwing out vectors
c     n_depen_check : number  of vectors to be thrown out 
c     o_depen_print : if true +> extra print 

c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      logical mpassi
      common/integmp/mpassi
c
caleko
c ---- for transformation drf-coordinates

c
      real*8 trmatr,trvec
      common /trnsf/ trmatr(3,3),trvec(3)
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
cdrf------------------------
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
      common/nottwi/obeen,obeen2,obeen3,obeen4
cdrf------------------------
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
      dimension zscf(12),zchar(21)
      parameter (nvec=9)
      dimension zvect(nvec+1)
      dimension yout(20),yuptyp(3),ymull(7)
      dimension ydiis(8)
      dimension words(120),irndg(50),irndc(50)
      integer atomss(maxat), natomss
      integer grd_type(maxat)
      equivalence (words(1),zdd3(1))
      equivalence (zopen, words(8))
      equivalence (grhf,words(13))
      equivalence (zoscf,words(12))
c
      integer mxbltyp
      parameter(mxbltyp = 50)

      logical oblur
      common/blurrr/oblur

      integer nbltyp, bltab(maxat), nutab(maxat)

      real*8 blexpo(maxat), blwght(maxat)
      real*8 blexpo2(mxbltyp), blwght2(mxbltyp)

      common/blurr2/blexpo, blwght,
     &     blexpo2, blwght2,
     &     bltab,nbltyp,nutab

      character*8 ztagbl
      common/blurrc/ztagbl(mxbltyp)
c 
      integer icoord,ncons_co,nimage_co,iopt_co,mem_co
      integer upd_co,rec_co,fd_co,maxc_co,dump_co
      integer task_co, po_pop_size_co, po_distribution_co
      integer po_maxcycle_co, po_init_pop_size_co, po_reset_co
      integer po_nsave_co, ntasks_co 
      real*8    delta_co,nebk_co,time_co,fri0_co,frif_co,frip_co
      real*8    soft_co,tol_co,maxs_co
      real*8    temperature_co, po_radius_co, po_contraction_co
      real*8    po_tolerance_r_co, po_tolerance_g_co
      real*8    po_mutation_rate_co, po_death_rate_co, po_scalefac_co
      logical odlfind,rst_co
      logical tdlf_farm_co
c real vars
      common/dlfind/delta_co,
     + time_co,fri0_co,frif_co,frip_co,
     +     nebk_co,
     + soft_co,tol_co,
     + temperature_co, po_radius_co, 
     + po_contraction_co, po_tolerance_r_co, po_tolerance_g_co,
     + po_mutation_rate_co, po_death_rate_co, po_scalefac_co,
c integer/logical vars
     + odlfind,maxc_co,icoord,ncons_co,nimage_co,
     +     iopt_co,mem_co,upd_co,rec_co,fd_co,
     +     maxs_co,dump_co,rst_co,
     + task_co, po_pop_size_co, 
     + po_distribution_co, po_maxcycle_co, po_init_pop_size_co, 
     + po_reset_co, po_nsave_co, ntasks_co, tdlf_farm_co
c character vars
      character(64) geom2
      character*256 geomfile
      common/dlfindc/geom2,geomfile
      integer maxfreeze,nfreeze,ifreeze
      parameter (maxfreeze = 1000)
      common/dlfindfreez/nfreeze,ifreeze(maxfreeze)
      character*8 zopti,zfrozen
      common/runopt/zopti(maxat),zfrozen
      logical o_var_high,odumdum
      common /csubst/ values(maxvar),intvec(maxvar),fpvec(maxvar),
     1                cmin10(maxvar),cmin20(maxvar),o_var_high,odumdum
c
c     Variables controlling the Fermi-Dirac smearing
c     Search for the subroutine fermi_smear...
c
c     osmear: should Fermi-Dirac smearing be used?
c     esmear: smearing parameter controlling the width of the Fermi 
c             cut-off function (Hartree)
c
c     esmear_start: the initial smearing parameter
c     esmear_final: the final smearing parameter
c     esmear_scale: the smear parameter follows the "tester" scaled
c                   with this factor from the start to the final value.
c
      logical osmear
      real*8    esmear_start, esmear_final, esmear_scale
c
      common/fermidirac/esmear_start,esmear_final,esmear_scale,
     &                  osmear
c
c
c     Constants for converting between atomic units and other units
c
      real*8 bohr_to_ang,   ang_to_bohr
      real*8 hartree_to_eV, eV_to_hartree
      parameter (bohr_to_ang   = 0.529177249d0)
      parameter (ang_to_bohr   = 1.0d0/bohr_to_ang)
      parameter (hartree_to_eV = 27.21165d0)
      parameter (eV_to_hartree = 1.0d0/hartree_to_eV)
     

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
c     MASSCF common fccwfn 
c
      integer nspacem
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd,fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspacem,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
c
      dimension ysupe(5),isup(3),yint(12),yau(3)
      dimension ynoprd(21),yprind(21)
      dimension ymult(9)
      dimension yblckd(60)
      dimension zgrp(19)
      character*8 zdumz(5)
      integer     idumz(5)
      character*10 charwall

      external all_data

      data yblckd /
     * 'coor','opti','init','conn','occu','vect','type','gvb',
     * 'forc','vibr','basi','eige','scfe','titl','over','grid',
     * 'mull','lowd','spin','leve','norm','twoe','grad','tran',
     * 'pdc ','gsym','secd','dma ','dens','ener','hess','dipo',
     * 'infr','rama', 26*'****'/
      data zchar/'ssss','psss','psps','ppss','ppps','pppp'
     &,'dsss','dsps','dspp','dsds','dpss','dpps','dppp'
     &,'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
      data ymult/
     *'sing','doub','trip','quar','quin','sext',
     *'sept','octe','none'/
      data yint / 'hond','high','rys',
     +            'g76','g80','low','rota',
     +            'game','tigh','255','mpas',' ' /
      data ynames/'ed7'/
      data zblank/' '/,yblank/' '/
      data ogeom/.false./
      data ysupe/'on','off','forc','noso','sort'/
      data yone,ytwo,yada/'one','two','adap' /
      data yscf,y4in,yhf,ygrad,ychf/
     +     'scf','tran','hf','grad','chf'/
      data isup/0,1,-1/
      data yiang,yiangs/'angs','ang '/
      data yau/'au','a.u.','aus'/
      data zc1/'c1'/
      data ynoprd/
     * 'dist','basi','adap','hist','vect','anal','hess',
     * 'ci'  ,'mrdc','ccsd','bq',9*'****','syma'/
      data yprind /
     * 'mull','lowd','scf','opti','gues','adap',
     * 'diis','symm','fock','vect','aopr','mopr','dump','lmo','stat',
     * 'dist','dscf','tran','inte','smea','syma' /
      data zend/'end'/
      data m1,m21/1,21/
      data yold/'old'/
      data yonee/'onee'/
      data prn/'star','dihe','alls','irre','mode',
     1         '    ','occv','allv','hell','grad',
     2         'mp2e','chfr','chfs','allo','mp2g',
     3         '    ','    ','    ','    ','geom',
     3         'dipo','quad','forc','chan','deri',
     4         'pola','    ','iner','    ','    ',
     +         's'   ,'t'   ,'t+v' ,'x'   ,'y'   ,
     +         'z'   , 4*' '/
c
      data all/'all'/
      data stop/'end'/
     2    ,yes /'yes'/, zno /'no'/, on /'on'/
      data zcore/'core'/
      data zclose/'closed'/
      data grp/'c1','ci','cs','cn','cnv','dn','cnh','dnh'/
      data zgrp /'c1','cs','ci','cn','s2n','cnh',
     + 'cnv' ,'dn'  ,'dnh' ,'dnd' ,'cinfv','dinfh','t'   ,
     + 'th'  ,'td'  ,'o'   ,'oh'  ,'i'   ,'ih'  /
      data zpunch/'vcd','fcm','dipder','polder','polariza',
     1 'vectors','geometry','dipole','gradient','optgeom',
     2 'polarisa','trans','raman','infrared','fullci','casscf',
     3 'mcscf', 23*' '/
      data noptx/120/
      data nopty/70/
      data irndg /
     * 1, 7, 7, 5, 8, 0, 0, 0, 0, 0,
     * 0,10,21,22, 0, 0,27, 4,15,23,
     *24,25,26, 0, 1,10, 0,10,28, 5,
     * 1, 0, 0, 0,22,27, 0,30,32,31,
     * 0,33,34,7*0 /
      data irndc/
     * 1, 0, 0,18, 4, 4, 2, 5, 0, 0,
     * 0, 0, 0, 0,19, 0, 0, 0, 0, 0,
     *13,14,20,21,22,23,17,29,33,38,
     *40,39,42,43, 0, 0, 0, 0, 0, 0,
     * 10*0 /
     *
      data ygvb/'gvb'/
      data ydiis/'off','prin','onse','test','dyna','leve','etes',' '/
      data yout /'cide','vect','cico','diag','supd',
     *           '----','----','all' ,'----','----',
     *           'nove','grad','----','----','conf',
     *           'lagr','cicy','hami','mo'  ,'ao'  /
      data yuptyp/'bfgs','ms','powe'/
      data oyes,ono/.true.,.false./
      data ymull/'atom','orbi','basi','grou','prin','unpa',' '/
      data zscf/'rhf','uhf','gvb','grhf','casscf','mcscf','multi',
     &          'direct','mp2','mp3','vb','masscf'/
      data zoimag/'intmag'/,zdelta/'delta'/
      data zatext,zbtext/'a','b'/
      data orgvb/.false./
      data zbeta/'beta'/
      data zvect/'minguess','extguess','hcore','getq','extra',
     *'alphas','cards','nogen','atoms','mosaved'/
      data yend/'end'/
      data m13,m53/13,53/
      data m0,m999,m2/0,999,2/
      data zold/'old'/
      data zstats/'a','d','delt','e','gamm','p','pi',
     *'s','sig-','sig+','t','delta','del','sigma-','sigm',
     *'sigma+','sigp','average' /
      data yrest/'rest'/
       data mit/200/
      data olevel/.false./
      data zoff/'off'/
c
c initialise defaults in CCPDFT input processor
      idum = CD_defaults(iwr)
cDEBUG
c     if (.not.opg_root()) ipu=99
cDEBUG
      o_var_high = .false.
      odft = .false.
      odftgrad = .false.
      oatmdft = .false.
      odenscfdft = .false.

      if (oresti) then
         if (oterm) then
            ostop = .true.
            go to 2790
         end if
         call initialb
         orest = .true.
      else
         orest = .false.
         call initiala
         ynamed = 'ed3'
      end if



      indi = .false.
      lres = .false.
      iprint = 0
      key = 1
      ocidu = ono
      onconv = ono
      oactiv = ono
      odci = ono
      oinv = ono
      jdeg = 1
      jrun = 1
      ostop = ono
      ostopr = ono
      onel = ono
      ocha = ono
      orun = ono
      open = ono
      oscf = ono
      oirc = ono
      ogvb = ono
      ogrhf = ono
      ouhf = ono
      omcscf = ono
      ovb = ono
      osmear = ono
      ifasti = 0
      ifastd = 1
      moutaa = 0
      moutbb = 0
      nav = lenwrd()
      mpassi = .false.
c
c     ----- initialise settings
c
      call pacini
      call inittxt
c
c     ----- atomic mass tables
c
      imassv = 1
      call mass_init
c


      oblur = .false.
      nbltyp = 0
      nghost = 0
      odisp = .false.
      do i=1,maxat
         zopti(i) = 'yes'
      end do
      zfrozen = 'no'
c
c     ----- read molecular data -----
c
      obasis = .true.
      osafe = .false.
      ofirst = .true.
      opass7 = .false.
      osetor = .false.
      oarch = .false.
      ognore = .false.
      osupe = .false.
c
c     check for any fatal errors from initial data pre-screening
c     in preinit
c
      if (okillv) then
       write(iwr,6011)
       call caserr2('fatal error detected in data pre-screening')
      endif
      if (oresti) then
         obas = .true.
c
         cpu = cpulft(1)
         write (iwr,6010) cpu,charwall()
         call restre
         local = 0
         runtyp = zblank
         scftyp = zblank
         call setlab
      else
         obas = .false.
         go to 30
      end if
c
c ----- now commence directive input and processing
c
      oservec = .false.
 20   call input
      if (oswit) then
         if (oservec) then
            ostop = .true.
            ostopr = .true.
            return
         end if
c
c ----- end of data input detected
c       terminate job if on 2nd pass
c

         if (oresti) then
            ostop = .true.
            call timana(1)
            go to 2790
         else
            oterm = .true.
            if (.not.(ofirst)) go to 2730
            ofirst = .false.
            if (.not.(obas)) go to 560
            go to 570
         end if
      end if
      oservec = .false.
 30   call inpa(ztext)
      ytext = ytrunc(ztext)
      i = locatc(ydd,limit1,ytext)
      if (i.gt.0) then
         go to (  140,  150,  230,  250,  260,
     +            200,  170,  270,  270,  290,
     +             20,  300,  300,  310,  310,
     +            320,  320,  330,  330,  370,
     +            340,  390,  350,  360,  400,
     +            450,  410,  460,  470,  490,
     +            540,  540,  550,  550,  550,
     +            550,  500,  430,  480,  440,
     +            530,  130,  130,   50,   40,
     +             60,   70,  120,  380,10010,
     +          10020,10030,10040,10050,10060,
     +            435,10070,10080,10100,10110,
     +          10120,10140,10150,10160,10165,
     +          10166,10170,10175,10177,10180,
     +          10190,10200,10210,10220,10230,
     +          10101,10400,10410,10420,10520) , i
      else if (ztext.eq.zold) then
c
c
c    ***** old **** directive facilitates using original
c                   gvb gamess data
c
         old = oyes
         go to 20
      else
         i = locatc(ydd(101),limit2,ytext)
         if (i.le.0) then
            buff = ztext
            i = locatc(words,noptx,buff)
            if (i.le.0) then

               call gamerr(
     &       'unrecognised DIRECTIVE or faulty DIRECTIVE ordering',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)

c               call caserr2(
c     +       'unrecognised DIRECTIVE or faulty DIRECTIVE ordering')

            else
               go to (200,640,700,690,710, 20,730,740,750,760,770,740,
     +                740,780,790,740,800,810,820,830,840,850,860,870,
     +                880,900,910,920,930,940,950,960,970,980,990,1030,
     +                1050,740,1070,1090,1110,1120,1060,1130,  20,1150,
     +                1160,1170,1180,1190,1200,1210,1220,1230,1250,1260,
     +                1270,1280,1290,1300,1310,1320,20,1000,1010,660,
     +                1240,1020,1040,20,1350,1350,1350,1350,1350,1350,
     +                1350,1350,1350,1350,1350,1350,1400,1340,20,20,
     +                1360,1350,1350,1350,1350,1350,1350,1390,1350,1350,
     +                1350,1350,1350,1350,1350,1350,1350,1350,1340,1360,
     +                1350,1350,1350,1350,1350,20,20,20,20,20,20,20,20,
     +                20) , i
            end if
         else
            if (((i.gt.60.and.i.le.96) .or. i.eq.39) .and. ofirst) then
               ofirst = .false.
               if (.not.(obas)) go to 560
               go to 570
            end if
            go to (2390,2400,2410,2420,2470,2440,2480,2490,2310,2320,
     +             2330,2340,2350,2360,2370,2380,2250,1870,1890,1900,
     +             1940,1950,1960,1980,2030,2040,2050,2060,2070,2010,
     +             1850,2090,2110,2120,2130,1830,2160,2170,2190,2200,
     +             1840,1820,1780,1790,1810,2630,1420,1430,1440,1450,
     +             1460,1470,1480,1490,1500,1510,1520,2650,2680,2710,
     +             2720,2430,2180,1730,1730,1730,1720,1670,1675,1690,
     +             1710,1660,1660,1610,1600,1600,1590,1570,1580,1530,
     +             1540,1550,1560,2460,2500,2485,1725,1735,1745,1746,
     +             1595,1747,1748,1749,10130,10130,10240,10250,10260,
     +            10300,10270,10275), i
         end if
      end if
c
c    #################################################
c     class 1 directives
c    #################################################
c
c     ----- field  (input x,y and z values for chf calculation)
c
 40   ofield = .true.
      call inpf(fieldx)
      call inpf(fieldy)
      call inpf(fieldz)
      go to 20
c
c     ----- format high (controls field width on vector print)
c
 50   call inpa4(ytext)
      call inpa(ztext)
      if (ytext.eq.yint(2)) then
         if (ztext(1:4).eq.'vari') then
            o_var_high = .true.
         else
            oprint(20) = .true.
         end if
      end if
      go to 20
c
c ...  combined <oneelec>. are partial charges
c      to be included in 1-electron part of hamiltonian?
 60   call inpa4(ytext)
      if (ytext.eq.yonee) oaminc = .true.
      go to 20
c ... mapping  ( ab initio -> mechanics mapping vector.
c     { igam imec chg } )
c ... read in ab initio -> mechanics mapping vector and
c     charges on atoms. spread over several lines,
c     terminated by an end directive.
 70   nrp = 0
 80   call inpa(ztext)
      if (ztext.ne.zend) then
         if (ztext.eq.zblank) then
            call input
            go to 80
         end if
         jrec = jrec - 1
         nrp = nrp + 1
         call inpi(natg)
 90      call inpa4(ytext)
         if (ytext.ne.yblank) then
            jrec = jrec - 1
            call inpi(mapmod(natg))
         else
            call input
            go to 90
         end if
 100     call inpa4(ytext)
         if (ytext.ne.yblank) then
            jrec = jrec - 1
            call inpf(tchg(natg))
         else
            call input
            go to 100
         end if
         go to 80
      else if (.not.((ozmat .or. ogeom) .and. (nrp.ne.nat))) then
         omapp = (nrp.ge.2)
         omodel = omapp
         go to 20
      end if
 110  write (iwr,6020) nrp , nat
      call caserr2('inconsistency in number of atoms specified.')
c
c     ----- freeze  -----
c
c     default  freeze all    .. all dummy atoms frozen
c                     off    .. no freezing of dummy elements
c                   first    .. freeze first real atom, others free
c     otherwise              .. freeze nuclei
c
 120  call inpa(ztext)
      if (ztext.eq.'off') ifree = 0
      if (ztext(1:4).eq.'firs') ifree = 1
      if (ztext.eq.'off'.or.ztext(1:4).eq.'firs'
     1.or.ztext.eq.'all') go to 20
c
c*** freeze nuclei 
c
      if (.not.ogeom.or..not.obas) 
     1   call caserr('freeze nuclei after geometry and basis')
      nfixat = 0
      zfrozen = 'yes'
      write(iwr,122)
121   if (ztext.eq.'end') then
         if (nfixat.gt.0) write(iwr,123)(zdumz(i),idumz(i),i=1,nfixat)
         nfixat = 0
         write(iwr,124) 
122   format(/8x,'***** freeze *****',
     1       /' ***** the following atoms *****')
123   format(5x,5(a8,1x,'(',i5,' atoms)'))
124   format( ' ***** This options is *not* guaranteed (jvl) *****',/)
         go to 20
      else
         nfixat = nfixat + 1
         zdumz(nfixat) = ztext
         idumz(nfixat) = 0
         wght = 1.0d100
         do j=1,nat 
            if (ztext.eq.ztag(j)) then
               zopti(j) = 'no'
               idumz(nfixat) = idumz(nfixat) + 1
               call mass_editvec(1,j,-2,wght)
            end if  
         end do  
         if (idumz(nfixat).eq.0) call caserr('unknown center in freeze')
         if (nfixat.eq.5) then
            write(iwr,123)(zdumz(i),idumz(i),i=1,nfixat)
            nfixat = 0
         end if
         call inpa(ztext)
         go to 121
      end if
c
c     ----- pseudo  -----
c
 130  if (.not.obas) call caserr2(
     +   'pseudo or ecp directive appears before basis directive')
      i10 = igmem_alloc(maxat)
      call inpseu(core(i10))
      call gmem_free(i10)
      go to 30
c
c     ----- safety ----
c
 140  osafe = .true.
      call inpf(safe)
      if (jump.eq.3) call inpf(dumtim)
      go to 20
c
c     ----- minblock ----
c
 150  if (jump.ne.3) then
         call caserr2('invalid syntax of minblock directive')
      end if
      call inpa4(ytext)
      i = locatc(yed,maxlfn,ytext)
      if (i.ne.0) then
         call inpi(minbl(i))
      else
         call caserr2('invalid ddname parameter')
      end if
      go to 20
c
c     ----- noprint ----
c
 170  if (jump.ne.1) then
 180     call inpa4(ytext)
         if (ytext.eq.yblank) go to 620
         i = locatc(ynoprd,21,ytext)
         if(i.gt.20) then
         if(i.eq.21) otsym=.false.
         else
         if (i.ne.0) then
            oprint(i+20) = .true.
         else
            i = locatc(prn,ntypr,ytext)
            if( i.ne.0) oprn(i) = .false.
         end if
         endif
         go to 180
      else
         do 190 i = 1 , 20
            oprint(i+20) = .true.
 190     continue
         do 195 i = 1 , ntypr
            oprn(i) = .false.
 195     continue
         otsym = .false.
         go to 20
      end if
 620  if (.not.oprn(7)) oprn(8) = .false.
      go to 20
c
c     ----- iprint / output
c
 200  if (jump.ne.1) then
 210     call inpa4(ytext)
c
         if (ytext.eq.'anal') then
            call inpi(ilopri)
            goto 20
         endif
c
         if (ytext.eq.yblank) go to 610
         i = locatc(yprind,21,ytext)
         if (i.gt.20) then
            if (i.eq.21) omydbg=.true.
         else
            if (i.ne.0) then
               oprint(i+40) = .true.
            else
               i = locatc(prn,ntypr,ytext)
               if (i.ne.0) oprn(i) = .true.
            end if
         endif
         go to 210
      else
         do 220 i = 1 , 20
            oprint(i+40) = .true.
 220     continue
         do 225 i = 1 , ntypr
            oprn(i) = .true.
 225     continue
         go to 20
      end if
 610  if (oprn(8)) oprn(7) = .false.
      if (oprn(7)) oprn(8) = .false.
c
c ----- bypass
c
 230  if (jump.ne.1) then
 240     call inpa4(ytext)
         if (ytext.ne.yblank) then
            if (ytext.eq.yone) opass1 = .true.
            if (ytext.eq.ytwo) opass2 = .true.
            if (ytext.eq.yold) opass3 = .true.
            if (ytext.eq.yada) opass4 = .true.
            if (ytext.eq.yscf) opass5 = .true.
            if (ytext.eq.y4in) opass6 = .true.
            if (ytext.eq.'4ind') opass11 = .true.
            if (ytext.eq.yhf)  opass7 = .true.
            if (ytext.eq.ygrad) opass8 = .true.
            if (ytext.eq.ychf) opass9 = .true.
            if (ytext.eq.'igno') ognore = .true.
            if (ytext.eq.'anal'.or.ytext.eq.'prop') opass10 = .true.
            go to 240
         else if (opass6.and..not.o4ind) then
            opass5 = .true.
            opass4 = .true.
         else if (opass5) then
            opass4 = .true.
         else if (.not.(opass4)) then
            if (opass7) opass5 = .true.
            go to 20
         end if
      end if
      opass2 = .true.
      opass1 = .true.
      go to 20
c
c ----- nosymmetry
c
 250  if (ozmat) then
         call caserr2(
     +   'nosym directive must precede zmatrix or geometry')
      else
         nosymm = 1
      end if
      go to 20
c
c ----- zmat
c
 260  if (ogeom) then

         call gamerr(
     &        'unrecognised DIRECTIVE or faulty DIRECTIVE ordering',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)

c         call caserr(
c     +   'unrecognised directive or faulty directive ordering')

      else
         call inpa4(ytext)
         if (ytext.eq.yiang .or. ytext.eq.yiangs) then
           ifau = 1
           if(jump.eq.3) call inpi(iseczz)
         else if (ytext.eq.'mopa') then
          irold = ird
          write(iwr,6170)
          call archiv ('zmat',iwr)
          ifau = 1
          oarch = .true.
         else if (locatc(yau,3,ytext).gt.0) then
           ifau = 0
           if(jump.eq.3) call inpi(iseczz)
         else
           if(jump.eq.2) then
            jrec=jrec-1
            call inpi(iseczz)
           end if
         end if
         if(iseczz.le.0) call caserr2('invalid z-matrix section')
         call rdgeom
         call symm(iwr,core)
         if (oarch) ird = irold
         ozmat = .true.
         ogeom = .true.
         oarch = .false.
         if (.not.((omapp) .and. (nrp.ne.nat))) go to 20
         go to 110
      end if
c
c ----- angstrom
c
 270  ifau = 1
      go to 20
c
c ----- linewidth
c
c280  call inpi(iwid)
c     call inpwid(iwid)
c     go to 20
c
c ----- adapt off/ions2
c
 290  if (jump.ne.1) then
         call inpa4(ytext)
         if (ytext.ne.ysupe(2)) then
            jrec = jrec - 1
            call inpi(ions2)
            if (ions2.le.0) ions2 = isect(482)
            otri = .false.
         else
            otran = .true.
         end if
      end if
      go to 20
c
c     ----- mfile/main
c
 300  call inpa4(ytext)
      if(ytext.eq.'memo'.or.ytext.eq.'in-c') then
       omem(3) = .true.
       if (jump.gt.2) call filmem(n2file,n2tape,maxset(3))
      else
       jrec = jrec-1
       call filein(n2file,n2tape)
      endif
      go to 20
c
c     ----- sfile/seco
c
 310  call filein(n4file,n4tape)
      go to 20
c
c     ----- ffile/tfile
c
 320  call filein(n6file,n6tape)
      go to 20
c
c     ----- dmfi/twop
c
 330  call filein(n11fil,n11tap)
      go to 20
c
c     ----- afile
c
 340  call filein(n1file,n1tape)
      go to 20
c
c     ----- loop
c
 350  call filein(n9file,n9tape)
      go to 20
c
c     ----- rloop
c
 360  call filein(ntfile,nttape)
      go to 20
c
c     ----- cifile
c
 370  call filein(n5file,n5tape)
      go to 20
c
c     ----- hfile
c
 380  call filein(n12fil,n12tap)
      go to 20
c
c     ---- punchfile -----
c
10010 continue
      call inpa4(ytext)
      iform = 0
      if(ytext.eq.'old ')then
         opunold = .true.
      else if(ytext.eq.'all ')then
         do i=1,60
            oblock(i) = .true.
            iblfmt(i) = iform
         enddo
      else if(ytext.eq.'high')then
         iform = 1
      else if(ytext.eq.'efor')then
         iform = 2
      else
         jrec = jrec - 1
      endif
      if((jump.eq.1).or.
     &   (iform.ne.0.and.jump.eq.2))then
c
c -- default set of blocks for simple visualisation
c     geom, conn, occu, vect, runt, scft, basi
c
         oblock(1)=.true.
         oblock(4)=.true.
         oblock(5)=.true.
         oblock(6)=.true.
         oblock(7)=.true.
         oblock(11)=.true.
         oblock(20)=.true.
c
         iblfmt(1)=iform
         iblfmt(4)=iform
         iblfmt(5)=iform
         iblfmt(6)=iform
         iblfmt(7)=iform
         iblfmt(11)=iform
         iblfmt(20)=iform
c
         go to 20
c
      endif
 7203 call inpa4(ytext)
      if(ytext.eq.yblank)go to 20
      if(ytext.eq.'high')then
         iform = 1
         goto 7203
      else if(ytext.eq.'efor')then
         iform = 2
         goto 7203
      endif
      i=locatc(yblckd,60,ytext)
      if(i)7205,7205,7204
 7204 oblock(i)=.true.
      iblfmt(i)=iform
c read optional integer values
      if(i.eq.5.or.i.eq.6.or.i.eq.12.or.i.eq.13)then
c
c --- vectors sections
c
 7206    call inpa4(ytext)
        if(ytext.eq.yblank)go to 20
	        call intval(ytext,ksec,onumb)
        if(onumb)then
           if(nblsec.eq.10)then
              write(iwr,7208)ksec,ytext
 7208         format(1x,'warning : numeric argument ',i3,
     &  ' on punchfile keyword ',a4,' was ignored',/,1x,
     &  'maximum number of vectors sections which may be specified ',
     &  'is 10')
           else
              nblsec=nblsec+1
              iblsec(nblsec)=ksec
           endif
           goto 7206
        else
           jrec = jrec - 1
        endif
      else if(i.eq.16)then
c
c --- grid sections
c
 7207    call inpa4(ytext)
        if(ytext.eq.yblank)go to 20
        call intval(ytext,ksec,onumb)
        if(onumb)then
           if(nblgri.eq.10)then
              write(iwr,7209)ksec,ytext
 7209         format(1x,'warning : numeric argument ',i3,
     &  ' on punchfile keyword ',a4,' was ignored',/,1x,
     &  'maximum number of grid sections which may be specified ',
     &  'is 10')
           else
              nblgri=nblgri+1
              iblgri(nblgri)=ksec
           endif
           goto 7207
        else
           jrec = jrec - 1
        endif
      else if(i.eq.22)then
c
c --- blocks of integrals
c
        call inpa4(ytext)
        if(ytext.eq.yblank)go to 20
        call intval(ytext,ksec,onumb)
        if(onumb)then
           itwoe = ksec
        else
           jrec = jrec - 1
           goto 7203
        endif
        call inpa4(ytext)
        if(ytext.eq.yblank)go to 20
        call intval(ytext,ksec,onumb)
        if(onumb)then
           ntwoe = ksec
        else
           jrec = jrec - 1
        endif
      endif
      go to 7203
 7205 call caserr2('unrecognised punchfile keyword '//ytext)
c
c     specifications for density functional calculation of the
c     correlation energy
c
10020 odenfu = .true.
      odenpt = .false.
      idenfu = 0
      nangpr = 2
10021 call inpa(ztext)
      ytext=ytrunc(ztext)
      if(ytext.eq.'symm')then
         osymm=.true.
      else if(ytext.eq.'nosy')then
c (default)
         osymm=.false.
      else if(ytext.eq.'sc  ')then
         osc=.true.
      else if(ytext.eq.'cs1 ')then
         idenfu = 1
      else if(ytext.eq.'cs2 ')then
         idenfu = 2
      else if(ytext.eq.'vwn1')then
         idenfu = 3
      else if(ytext.eq.'vwn2')then
         idenfu = 4
      else if(ytext.eq.'slat')then
         idenfu = 5
      else if(ytext.eq.'p1')then
         idenfu = 6
      else if(ytext.eq.'p2')then
         idenfu = 7
      else if(ytext.eq.'poin')then
         call inpf(denptx)
         call inpf(denpty)
         call inpf(denptz)
         odenpt = .true.
      else if(ytext.eq.'core')then
         call inpi(ncorxx)
      else if(ytext.eq.'grid')then
         call inpi(nangpr)
      else if(ytext.eq.'    ')then
         if(osc.and.osymm)call caserr2
     &   ('symm and sc keywords on denfun directive are incompatible')
         if(idenfu.eq.0)call caserr2
     &        ('no density functional type specified')
         goto 20
      else
         call caserr2('bad keyword on denf directive')
      endif
      goto 10021
c
10030 oscrf = .true.
      call inpf(aradix)
      call inpf(aradiy)
      call inpf(aradiz)
      call inpf(dielec)
      go to 20
c
c    ------ psscrf (input dielectric, point charge density,
c           no. of centres, radii of spheres, maxit, damper)
c
10040 opssc = .true.
      iscf = 0
      call inpf(delect)
      call inpf(ptspac)
      call inpi(nocent)
      do 10041 id=1,nocent
      if(jrec.ge.jump) call input
      call inpf(cavsze(id))
10041 continue
      call inpi(itmax)
      if (itmax.eq.0)then
c     ----- defaults dvsg -----
       itmax=4
       damper=0.7d0
       go to 20
      endif
      call inpf(damper)
      go to 20
c
c orient 
c
10050 osetor = .true.
      call inpa(ztext)
      igpold = locatc(zgrp,19,ztext)
      call inpi(naxold)
      if(igpold.le.0.or.jump.ne.3)call
     +     caserr2('missing or invalid pt group on orient directive')
      do 10052 ii = 1,3
         call input
         if(jump.ne.4)call 
     +        caserr2('insufficient data on orient matrix record')
         do 10051 jj = 1,3
            call inpf(trold(jj,ii))
10051    continue
         call inpf(trano(ii))
10052 continue
      goto 20
c
c  symtol
c
10060 call inpi(isymtl)
      if(isymtl.lt.1.or.isymtl.gt.20)call 
     +     caserr2('invalid symtol value')
      goto 20
c
c  degecr
c
10070 call inpf(degecr)
      otsym= (degecr.ge.0d0) 
      write(iwr,6185) degecr
      goto 20
c
c  collect atoms into several groups (for lowdin analysis)
c
10080 call inpi (ngroup)
      if (ngroup .gt. mng) 
     + call caserr2('too many groups of atoms')
      do 10082 jgroup=1,ngroup
         call setsto(mpergr,0,iatom(1,jgroup))
         call input
         if (jump .gt. mpergr) 
     +   call caserr2('too many atoms in this group')
         do 10081 j=1,jump
            call inpi(iatom(j,jgroup))
10081    continue
10082 continue
      goto 20
c
c versions
c
10100 continue
      call prtver(iwr)
      call prsize(iwr)
      call wrtkeys(iwr)
      call inpa4(ytest)
      if (ytest.eq."cont") then 
         goto 20
      elseif (ytest.ne." ".and.ytest.ne."stop") then
         call caserr2(
     +   "invalid versions sub-directive: use VERSIONS [STOP|CONTINUE]")
      endif
      call whtps
      call shut
      call clenup2
      call pg_end(0)
      goto 20
c
c keys
c
10101 continue
      call wrtkeys(iwr)
      goto 20
c
c   -------- xfield (a flag to indicate that bq atoms
c            are to be treated as a model for a crystal
c            field)
c
10110 continue
      ocryst = .true.
      goto 20
c
c morokuma - perform morokuma energy decomposition
c            analysis
c
c routes code through morokscf driver, forces nopk=1
c
10120 continue
      call inpa(ztext)
      ytext=ytrunc(ztext)
      if(ytext.eq.'frag') then
         call inpi(imode)
         call inpa(ztext)
         if(imode .ne. 1 .and. imode .ne. 2)
     &        call caserr2('bad morokuma fragment index')
         fragname(imode) = ztext
         omorok = .true.
      else 
         call inpa(ztext)
         fragname(1) = ztext
         call inpa(ztext)
         fragname(2) = ztext
         call inpi(imode)
         if(imode.eq.0)imode=3
         omorok = .true.
      endif
      nopk = 1
      goto 20
c
c     vdwaals: van der Waals energy "corrections"
c
10300 continue
      call vdwaals_input
      goto 20
C **********************************************************************
C *Bolt in:     ccpdft                                                 *
C *             implement cdft directives                              *
C **********************************************************************
10130 continue
      idum = CD_request()
c
c
      if (.not.odft) then
c
c     establish a default set of attributes at the first occurrence of
c     a cdft directives. 
c
        call dft_defaults(atomss,grd_type,natomss,ierr0,ierr1)
c
c       do not attempt to re-set these defaults downstream 
c
        odft = .true.
c
c       for runtypes hessian and force turn the gradients of the 
c       quadrature on.
c
        if (jrun.eq.8.or.jrun.eq.24) ierror = CD_gradquad(.true.)
c
      endif
c
      call dft_input(atomss,grd_type,natomss,odftgrad,ierr0,ierr1,
     &               ierror)
c
      go to 20 
C **********************************************************************
C *End ccpdft bolt in                                                  *
C **********************************************************************
c
c...   react
c
c
c
c  ***************************************************************
c
c      react, i.e. reaction field option is being processed
c
c  **************************************************************
c
c
10140 oreact = .true.
      call rfin                               
      goto 20
c
c...   harmonic
c
10150  oharm = .true.
       if (obas) call caserr2(
     +   'harmonic directive appears after basis directive')
10151  call inpa4(ytext)
       if (ytext.eq.'off') oharm = .false.
       if (ytext.eq.'on')  oharm = .true.
       if (ytext.eq.'prin') opharm = .true.
       if (ytext.eq.yblank) then
          if (.not.oharm) write(iwr,10153) 
10153  format(' ** Cartesian Gaussians are used **')
          go to 20
       end if
       go to 10151
c
c...   zora
c
10160  ozora = .true.
c Jens to remove calls to Zora
       if (.not.orest .and. obas) 
     +    call caserr2(
     +    'zora directive presented AFTER basis directive')
       call zora_input
       go to 30

c     
c ...  chm
c
c     This directive is just used to control
c     diagnostic output for the charmm interface
c     for a non-charmm build it is ignored
c
10175  continue
      goto 20
c
c... newscf
c
10170  continue
       call gamerr(
     & 'newscf option is not supported by this build', 
     & ERR_NO_CODE, ERR_SWITCHED_OUT, ERR_SYNC, ERR_NO_SYS)
      go to 20
c
c ----- bqbq
c
10165 continue
      omtslf = .false.
      go to 20
c
c ---- blur
c
c  allows specific tags to be associated with
c  a gaussian charge distribution
c
10166 continue
      oblur = .true.
10167 call input
      call inpa(ztext)
      if(ztext.eq.zend)goto 10168
      if(jump.le.1)
     &     call caserr2('bad blur specification')
      nbltyp = nbltyp + 1
      if(nbltyp .gt. mxbltyp) 
     &     call caserr2('too many blurred charge types')
      ztagbl(nbltyp) = ztext
      call inpf(blexpo2(nbltyp))
      call inpf(blwght2(nbltyp))
      goto 10167
10168 continue
      go to 20
C
C*** Calling CRESTR
C
10180 continue
      oscf = .true.
      mp2 = .false.
      mp3 = .false.
      lmcscf = .false.
      lci = .false.
      mcgr = .false.
      cigr = .false.
      zscftp = 'vb      '
      nopk = 1
      otri = ono
      if (intg76.ne.0) intg76 = 3
      call vbstart(core,0)
      goto 20
C
C*** Calling vbin for input TURTLE
C
10190 continue
      if (opssc) then
          call mod1ed(c,czan,1)
      endif
      call vbstart(core,1)
      goto 20
c
c  --- dependency -- 
c      exp(-i_depen_check) is used as criterion for dependent basis
c
10200 call inpa4(ytest)
      if (ytest.eq.'crit') then
         call inpi(i_depen_check)
         i_depen_check = max(i_depen_check,0)
         write(iwr,10201) i_depen_check
10201    format(/' ** dependent vectors eliminated at 10**(-',i3,
     1           ') (0=-20) **')
      else if (ytest.eq.'numb'.or.ytest.eq.'quan') then
         call inpi(n_depen_check)
         write(iwr,10202) n_depen_check
10202    format(/' **',i5,' possibly dependent vectors will be'
     1          ,'  eliminated **')
      else if (ytest.eq.'prin') then
         o_depen_print = .true.
      else
         go to 20
      end if
      go to 10200
c
c*** calling servec for general vector manipulations .....
c
10210 if (ztext.ne.'servec'.and.ztext.ne.'server')  
     1   call caserr2('servec and server should be called by'//
     2               ' their full name')
      if (ztext.eq.'servec') then
         call servec(core,'start')
         oservec = .true.
         orest = .true.
      else if (ztext.eq.'server') then
         call gserv(core,oquit)
         if (oquit) go to 2790
      end if
      go to 20
c
c*** ghost : denote atoms as ghost atoms for easier counterpoise
c   
10220  call inpa(ztext)
       if (ztext.eq.' ') then
          call input
          go to 10220
       else if (ztext.eq.'end') then
          go to 20
       else
          nghost = nghost + 1
          if (nghost.gt.maxghost) call caserr2('too many ghosts')
          ztag00(nghost) = ztext
          go to 10220
      end if
c
c*** disp  polarisation dispersion energy
c
10230  odisp = .true.
       call inpa4(ytest)
       call inpi(iblk_a)
       call inpi(isec_a)
       call filec(ytest,iblk_a,iunit_a,irep)
       if (irep.ne.0) call caserr2('error in disp directive')
       call inpa4(ytest)
       call inpi(iblk_b)
       call inpi(isec_b)
       call filec(ytest,iblk_b,iunit_b,irep)
       if (irep.ne.0) call caserr2('error in disp directive')
       otran = .true.
       if (nosymm.ne.1) call caserr2('disp only with nosym')
       call outrec
       go to 20
c
c ----- bqforce
c
10177 continue
      obqf = .true.
      go to 20
c
c ----- smear (fermi-dirac smearing)
c
c     general form of input is:
c
c       smear <start temp> [<end temp>] [au|hartree|eV] 
c             [scale <scale factor>]
c
c     the end temperature and the units are optional. 
c     the default end temperature is 0
c     the default unit is Hartree 
c     the default scale factor is 1.
c
10240 continue
      osmear = oyes
      esmear_final = 0.0d0
      esmear_scale = 1.0d0
      call inpf(esmear_start)
      if (esmear_start.lt.0.0d0) then
         call caserr2("the smearing parameter must be positive")
      endif
      call inpa(ytest)
      if (ytest.eq.'au'.or.ytest.eq.' '.or.ytest.eq.'hart') then
      else if (ytest.eq.'ev') then
         esmear_start = esmear_start*eV_to_hartree
      else if (ytest.eq.'scal') then
         call inpf(esmear_scale)
         if (esmear_scale.le.0.0d0) then
           call caserr2("the scale factor must be positive")
         endif
      else
        jrec = jrec - 1
        call inpf(esmear_final)
        if (esmear_final.lt.0.0d0) then
           call caserr2("the smearing parameter must be positive")
        endif
        call inpa(ytest)
        if (ytest.eq.'au'.or.ytest.eq.' '.or.ytest.eq.'hart') then
        else if (ytest.eq.'ev') then
          esmear_start = esmear_start*eV_to_hartree
          esmear_final = esmear_final*eV_to_hartree
        else if (ytest.eq.'scal') then
          call inpf(esmear_scale)
          if (esmear_scale.le.0.0d0) then
            call caserr2("the scale factor must be positive")
          endif
        else
          write(iwr,*)'Unknown units: ',ytest
          call caserr2("Unknown units specified")
        endif
      endif
      call inpa(ytest)
      if (ytest.eq.'scal') then
         call inpf(esmear_scale)
         if (esmear_scale.le.0.0d0) then
           call caserr2("the scale factor must be positive")
         endif
         call inpa(ytest)
      endif
      if (ytest.ne.' ') then
         call caserr2("unknown fermi-smearing subdirective")
      endif
      if (esmear_start.lt.esmear_final) then
         if (opg_root()) then
          write(iwr,*)"The starting temperature you specified is ",
     &              esmear_start," (Hartree)"
          write(iwr,*)"The final    temperature you specified is ",
     &              esmear_final," (Hartree)"
          write(iwr,*)"However the starting smearing temperature must ",
     &              "be equal to or higher than the final temperature"
         endif
         call caserr2(
     &      "The starting smearing lower than final temperature")
      endif
      go to 20
c
c ----- average (state averaging in llvmo)
c
c     general form of input is:
c
c       average on|off|<tolerance>
c
c     the state averaging is "on" by default and
c     the default tolerance currently is 1.0d-4 Hartree
c
10250 continue
      call inpa(ytest)
      if (ytest.eq."on".or.ytest.eq."off") then
         call llvmo_average(ytest)
      else
         jrec = jrec - 1
         call inpf(rtest)
         call llvmo_average("on")
         call llvmo_tolerance(rtest)
      endif
      go to 20
c
c ----- diesel (Direct individually selecting MR-CI code)
c
c     this will start a part of the diesel C++ code
c     the first line contains the program to start
c     the rest is input for that program
c     the reading of the input starts at end
10260 continue
      write (iwr,*) "Reading the input for DIESEL"
      call dieselinp
      odiesel = oyes
      omrdci = ono
      call mrdcin2
c Reset opark to false, it has been set to true by mrdcin2
      opark = ono
      go to 20

c
c ----- dl-find geometry optimiser
c
10270 continue
      odlfind = .true.
c     Defaults: in general, set to negative values to use the defaults 
c      from dl-find.f90
      icoord=0
      ncons_co=0
      iopt_co=-1
      mem_co=-1
      time_co=-1.D0
      fri0_co=-1.D0
      frif_co=-1.D0
      frip_co=-1.D0
      nimage_co=-1
      nebk_co=-1.D0
      delta_co=-1.D0
      upd_co=-1
      rec_co=-1
      fd_co=2
      soft_co=1.D20
      tol_co=0.003D0
      maxc_co=60
      maxs_co=-1.D0
      geom2=''
      geomfile=''
      dump_co=0
      rst_co=.false.
c JMC new variables
      task_co=-1
      temperature_co=-1.D0
      po_pop_size_co=-1
      po_radius_co=-1.D0
      po_contraction_co=-1.D0
      po_tolerance_r_co=-1.D0
      po_tolerance_g_co=-1.D0
      po_distribution_co=-1
      po_maxcycle_co=-1
      po_init_pop_size_co=-1
      po_reset_co=-1
      po_mutation_rate_co=-1.D0
      po_death_rate_co=-1.D0
      po_scalefac_co=-1.D0
      po_nsave_co=-1
      ntasks_co=1
      nfreeze=0
c     tdlf_farm_co=

c     loop over input lines
10271 call input
      call inpa4(ytext)
      if(ytext.eq.yend)go to 20
      if(ytext.eq.'test')then
        call inpi(ntest)
      else if(ytext.eq.'coor')then
        call inpa4(ytext)
        if(ytext.eq.'mass')then
          icoord=icoord+5
        else if(ytext.eq.'hdlc')then
          icoord=icoord+1
        else if(ytext.eq.'htc')then
          icoord=icoord+2
        else if(ytext.eq.'dlc')then
          icoord=icoord+3
        else if(ytext.eq.'tc')then
          icoord=icoord+4
        else if(ytext.eq.'cart')then
        else
          write (iwr,*) "DL-FIND: argument to coordinates: ",ytext,
     &         " not recognised"
          call gamerr(
     &         'DL-FIND DIRECTIVE not recognised',
     &         ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
        endif
c     NEB: nimage, nebk
      else if(ytext.eq.'neb')then
        icoord=icoord+120
        call inpi(nimage_co)
        call inpf(nebk_co)
        if(nebk_co.eq.0.D0) nebk_co=-1.D0
c     NEBMODE
      else if(ytext.eq.'nebm')then
        call inpa4(ytext)
        if(ytext.eq.'free') then
          icoord=icoord-20
        else if(ytext.eq.'perp') then
          icoord=icoord-10
        else if(ytext.eq.'cart') then
          icoord=icoord+30
        end if
        if(ytext.ne.'cart') then
          call inpa4(ytext)
          if (ytext.eq.'cart') icoord=icoord+30
        end if
c     Dimer method: maxrot, nointerpolation, tolrot, delta
      else if(ytext.eq.'dime')then
        icoord=icoord+220
        call inpa4(ytext)
        if(ytext.eq.'noin')then
          coord=icoord-10
        end if
      else if(ytext.eq.'delt')then
        call inpf(delta_co)
        if(delta_co.eq.0.D0) delta_co=-1.D0
c     Optimiser
      else if(ytext.eq.'opti')then
        call inpa4(ytext)
        if(ytext.eq.'lbfg')then
          iopt_co=3
          call inpi(mem_co)
          if(mem_co.eq.0) mem_co=-1
        else if(ytext.eq.'prfo')then
          iopt_co=10
        else if(ytext.eq.'cg')then
          iopt_co=2
        else if(ytext.eq.'sd')then
          iopt_co=0
        else if(ytext.eq.'dyn')then
          iopt_co=30
          call inpf(time_co)
          call inpf(fri0_co)
          call inpf(frif_co)
          call inpf(frip_co)
c JMC parallel opts input
        else if(ytext.eq.'stoc')then
          iopt_co=51
          call inpi(po_pop_size_co)
          call inpf(po_radius_co)
          call inpf(po_contraction_co)
          call inpf(po_tolerance_r_co)
          call inpf(po_tolerance_g_co)
          call inpi(po_distribution_co)
          call inpi(po_maxcycle_co)
          call inpf(po_scalefac_co)
        else if(ytext.eq.'gene')then
          iopt_co=52
          call inpi(po_pop_size_co)
          call inpi(po_init_pop_size_co)
          call inpf(po_radius_co)
          call inpf(po_tolerance_g_co)
          call inpi(po_maxcycle_co)
          call inpi(po_reset_co)
          call inpf(po_mutation_rate_co)
          call inpf(po_death_rate_co)
          call inpi(po_nsave_co)
        end if
c     UPDATE
      else if(ytext.eq.'upda')then
        call inpa4(ytext)
        if(ytext.eq.'bofi')then
          upd_co=2
        else if(ytext.eq.'powe')then
          upd_co=1
        else if(ytext.eq.'none')then
          upd_co=0
        end if
        call inpi(rec_co)
        if(rec_co.eq.0) rec_co=-1
        call inpi(fd_co)
        if(fd_co.eq.0) fd_co=-1
        call inpf(soft_co)
        if(soft_co.eq.0.D0) soft_co=1.D20
      else if(ytext.eq.'xtol') then
        call inpf(tol_co)
        if(tol_co.eq.0.D0) tol_co=0.003D0
      else if(ytext.eq.'minm') then
        call inpi(maxc_co)
        if(maxc_co.eq.0) maxc_co=60
      else if(ytext.eq.'step') then
        call inpf(maxs_co)
        if(maxs_co.eq.0.D0) maxs_co=-1.D0
c     2nd set of input coordinates
      else if(ytext.eq.'geom') then
        call inpanpc(geomfile)
        geom2=geomfile(1:64)
c     dump
      else if(ytext.eq.'dump') then
        call inpi(dump_co)
c     restart
      else if(ytext.eq.'rest') then
        rst_co=.true.
c     Constraints
      else if(ytext.eq.'cons')then
10272   call input
        call inpa4(ytext)
        if(ytext.eq.yend)go to 10273
        if (ytext.eq.'froz') then
           nfreeze=nfreeze+1
           call inpi(itest)
           ifreeze(nfreeze)=itest
        endif
        ncons_co=ncons_co+1
        go to 10272
10273   continue
c JMC taskfarming input
c (not used in this version, so ntasks_co=1 is fine as a default [just
c needs to be > 0].  Gamess itself will make the task farms and DLF
c will get the info via a call to dlf_get_taskfarm, as tdlf_farm is .false.)
      else if(ytext.eq.'task')then
        call inpi(ntasks_co)
      else
        write (iwr,*) "DL-FIND keyword ",ytext," not recognised"
        call gamerr(
     &       'DL-FIND DIRECTIVE not recognised',
     &       ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
      endif
c     end of loop over input lines
      go to 10271
      go to 20
c
c ----- masscf data input
c
10275 omas = oyes
      call masscf_input
      go to 30
c
10400 continue
10410 continue
      call gamerr('this executable was not build with xml options',
     &        ERR_NO_CODE, ERR_SWITCHED_OUT, ERR_SYNC, ERR_NO_SYS)
c
c 
c      copq repeat
c
10420 call set_copq_again('init')
      go to 20
c
c     cpwa - cputime vs walltime
c
10520 call inpf(cpwa)
      write(iwr,10521) cpwa
10521 format(' ******** allowed relation between cpu and wall',g12.3)
      go to 20
c
c     ----- maxblock edx iblk
c
 390  if (jump.ne.3) then
         call caserr2('invalid syntax of maxblock directive')
      else
         call inpa4(ytext)
         i = locatc(yed,maxlfn,ytext)
         if (i.ne.0) then
            call inpi(maxbl(i))
         else
            call caserr2('invalid ddname parameter')
         end if
      end if
      go to 20

c
c ----- dumpfile
c
 400  call inpa4(ytext)
      if(ytext.eq.yblank) go to 20
      if(ytext.eq.'memo'.or.ytext.eq.'in-c') then
       omem(4) = .true.
       ynamed = 'ed3'
       call inpi(idum)
       if(idum.gt.0) maxset(4) = idum
       go to 400
      endif
      if(ytext.eq.'scra') then
       omem(8) = .true.
       call inpi(idum)
       if(idum.gt.0) maxset(8) = idum
       go to 400
      endif
      if (ytext.eq.'mini') then
c   .. minimise sections o/p to dumpfile
       iacsct(205) = -1
       iacsct(206) = -1
       iacsct(207) = -1
       iacsct(208) = -1
       iacsct(209) = -1
       iacsct(420) = -1
       go to 400
      endif
      ynamed = ytext
      if (jump.ge.3) then
         call inpi(ibl3d)
         call inpa4(ytext)
         if(ytext.eq.'size') then
          call inpi(maxb)
          if(maxb.le.0) maxb = 99999
         else
          jrec = jrec - 1
          call inpi(ionsec)
          if (ionsec.le.0) ionsec = isect(492)
         endif
      end if
      go to 20
c
c ----- title -----
c
 410  call input
      k = 1
      do 420 i = 1 , 10
         ztitle(i) = char1(k:k+7)
         k = k + 8
 420  continue
      write (iwr,6030) ztitle
      go to 20
c
c ----- charge ----
c
 430  if (ogeom.and..not.oresti) then
         call caserr2(
     +  'charge directive must precede geometry specification')
      else if (jump.ne.2) then
         call caserr2('error on charge directive')
      else
         if(onel)call caserr2(
     +        'electrons and charge directives are incompatible')
         call inpi(ich)
         ocha = .true.
      end if
      go to 20
c
c ----- elec ----
c
 435  if (ogeom) then
         call caserr2(
     +   'electrons directive must precede geometry specification')
      else if (jump.ne.2) then
         call caserr2('error on electrons directive')
      else
         if(ocha)call caserr2(
     +        'charge and electrons directives are incompatible')
         call inpi(ne)
         onel = .true.
      end if
      go to 20
c
c ----- mult ----
c
 440  if (ogeom.and..not.oresti) then
         call caserr2(
     +  'multiplicity directive must precede geometry specification')
      else
         call inpa4(ytext)
         mul = locatc(ymult,9,ytext)
         if (mul.le.0) then
            jrec = jrec - 1
            call inpi(mul)
            if (mul.le.0) call caserr2(
     +      'invalid multiplicity specified')
         end if
      end if
      go to 20
c
c ----- integral -----
c
 450  call inpa4(ytext)
      i = locatc(yint,12,ytext)
      if (i.eq.0.or.i.eq.12) then
         i = 5
         go to 20
      end if
      if (i.le.3) then
       intg76 = 0
      else if(i.eq.11) then
       mpassi = .true.
      else if (i.ge.4.and.i.le.7) then
       if (i.eq.6) ifastd = 0
       intg76 = 2
      else if (i.eq.8) then
       ifasti = 1
       ifastd = 2
      else if (i.eq.10) then
       o255i = .false.
      else
       ifasti = 2
       ifastd = 2
      endif
      go to 450
c
c ----- print -----
c
 460  call inpi(nprint)
      go to 20
c
c ----- accu ----
c
 470  call inpa4(ytext)
      if (ytext.eq.'dire') then
         call inpi(iaccur)
c**** accuracy of integral evaluation as determined by accuracy
c**** of incomplete gamma function calculation
c**** (don't vary this accuracy)
         ocvary = .false.
      else if (ytext.eq.'diag') then
c***  set accuracy of all diagonalisations
         call inpi(i_global_diag)
         if (i_global_diag.lt.-990) call caserr2('accuracy diag error')
         call outrec
      else
         jrec = jrec - 1
         call inpi(itol)
         call inpi(icut)
      end if
      go to 20
c
c ----- norm ----
c
 480  call inpi(normf)
      call inpi(normp)
      go to 20
c
c ----- supe ----
c
 490  if (ztext.eq.'superci') go to 1830
      osupe = .true.
      offsym = .false.
      iofsym = 0
 495  call inpa4(ytext)
      if (ytext.eq.yblank) go to 20
      i = locatc(ysupe,3,ytext)
      if (i.ne.0) then
         nopk = isup(i)
      else
         if (ytext.eq.ysupe(4)) then
            osortp = .false.
            go to 495
         end if
         if (ytext.eq.ysupe(5)) then
            osortp = .true.
            go to 495
         end if
c NB schwarz directive
         if (ytext.eq.'nosc') then
            oschw  = .false.
            go to 495
         end if
         if (ytext.ne.ydd(4))
     +       call caserr2('parameters of super directive invalid')
         offsym = .true.
      end if
      go to 495
c
c ----- symmetry -----
c
 500  if (jrec.lt.jump) then
         call inpf(toler)
         call inpf(toler2)
         write(iwr,*)
     +  ' symmetry tolerance thresholds changed to ',toler,toler2
      end if
      osim1=.false.
      osim2=.false.
      if (ozmat) then
         call caserr2(
     +   'symmetry directive must precede geometry specification')
      else
         call input
         call inpa(zgroup)
         call inpi(jaxis)
         call input
         call inpa(ztext)
         if (ztext.ne.zend) then
            jrec = 0
            osim1=.true.
            do 510 i = 1 , 6
               call inpf(symtag(i))
 510        continue
            call input
            call inpa(ztext)
            if (ztext.ne.zend) then
               jrec = 0
               osim2=.true.
               do 520 i = 7 , 9
                  call inpf(symtag(i))
 520           continue
               call inpa(zsymm)
               call input
               call inpa(ztext)
               if (ztext.ne.zend) call caserr2(
     +             'parameters of symmetry directive invalid')
            end if
         end if
      end if
      go to 20
c
c ----- restore task ----
c
 530  if (ztext.eq.'restrict') go to 950
      if (ztext.eq.'restore') then
         call rstart(core,orest,ynamed)
         obas = .true.
         ogeom = .true.
         orun = .true.
      else
         call restore(orest,ynamed)
      end if
c
c     reset multiplicity and charge to defaults
c
      mul = 1
      ich = 0
c
      scftyp = zblank
      runtyp = zblank
      go to 20
c
c ----- atoms/geom  -
c
 540  if (ogeom) then
         call caserr2( 'geometry has already been specified')
      else
         call inpa4(ytext)
         if (ytext.eq.'mopa') then
          irold = ird
          write(iwr,6180)
          call archiv ('geom',iwr)
          oarch = .true.
         endif
         call readat(core,oarch)
         if (oarch) ird = irold
         if(iseczz.le.0) call caserr2('invalid z-matrix section')
         if (zgroup.eq.zc1) call symm(iwr,core)
         ogeom = .true.
         oarch = .false.
         if (.not.((omapp) .and. (nrp.ne.nat))) go to 20
         go to 110
      end if
c
c ----- basis/gtos/stos/gauss
c
 550  if (obas) then
         call caserr2('basis directive presented twice')
      endif
      obas = .true.
      if (.not.orest. and. (orun.or.oscf.or.oinv)) call caserr2(
     +'basis directive presented AFTER runtype, scftype or vectors')
      if (.not.orest. and. odft) 
     +call caserr2('basis directive presented AFTER dft directive')
 560  continue   



c     Read the geometry from an xyz file
CMR:  was  call dlf_geom(maxat,ztag,cat,nat_)
      if (lgt(geomfile(1:1),' ').and.
     >    lle(geomfile(1:1),'~')) then
        call read_xyz(geomfile,maxat,ztag,cat,nat)
        if(nat.ne.0) then
          do iat=1,nat
            do j=1,3
              cnew(iat,j)= cat(j,iat)
              c(j,iat)   = cat(j,iat)
            enddo
            call lcase(ztag(iat))
            zaname(iat) = ztag(iat)

            iznum=jsubst(ztag(iat))
            czan(iat) = iznum
            czann(iat) = iznum
            symz(iat)  = czan(iat)
            imass(iat) = idint(czan(iat)+0.01d0)

          enddo
          ogeom = .true.
          nonsym = nat
c         call symm(iwr, core)
        endif
      endif

      call filec(ynames,m1,num8,irep)
      if (irep.ne.0) call caserr2(
     +             'scratchfile has been incorrectly assigned')
      if (.not.orest) then
         call filec(ynamed,ibl3d,idaf,irep)
         if (irep.ne.0) call caserr2(
     +              'dumpfile has been incorrectly assigned')
         numdu = idaf
         iblkdu = ibl3d
         orevis(1) = .false.
         len501 = lensec(ldsect(isect(501))) + lensec(lds(isect(501)))
         call secput(isect(501),m21,len501,ibl3rs)
         call openda
         irest = 0
         write (iwr,6040)
      end if

      write (iwr,6050) ynamed , iblkdu, ynames , m1
c
c     generate the transformation matrices.
c

      call ptgrp(1)
c
c     read the basis set for the unique centres,
c     and generate the molecular basis set.
c
      if (.not.ogeom) then
         call caserr2(
     +  'attempt to define basis before geometry specified')
      else
c
c     don't need a basis set for a mechanics or a dynamics run .
c ...
         if (jump.gt.1) then
            call inpa4(ytext)
            if (ytext.eq.ysupe(2)) then
               obasis = .false.
               go to 20
            else
               jrec = jrec - 1
            end if
         end if

         ilen = 2*maxat*48/nav + 2*mxshel*48/nav + 10*mxprim + maxat 
         i10 = igmem_alloc(ilen)
         i11 = i10 + maxat*48/nav
         i12 = i11 + mxprim
         i13 = i12 + mxprim
         i14 = i13 + mxprim
         i15 = i14 + mxprim
         i132 = i15 + mxprim
         i16 = i132 + maxat
         i17 = i16 + mxprim
         i18 = i17 + mxprim
         i19 = i18 + mxprim
         i20 = i19 + mxprim
         i25 = i20 + mxprim
         i30 = i25 + mxshel*48/nav
         i40 = i30 + mxshel*48/nav
         last = i40 + maxat*48/nav

         length = last - i10
         if(length .ne. ilen)call caserr2('mem size error')
         call atoms(core,core(i11),core(i12),core(i13),core(i14),
     +        core(i15),ztita,core(i132),zatagg,
     +        core(i16),core(i17),core(i18),core(i19),core(i20),
     +        core(i10),core(i25),core(i30),core(i40),
     +        maxat,mxshel,mxprim,obas,onel,oshell)

         call gmem_free(i10)
c
         call setlab
         call putzm(0,iseczz)

c
c Import geometry to DFT module
c It is done here (i.e. before we know if we are 
c really running a DFT calculation) because we 
c will need it either 
c a) when a dft directive is encountered
c b) if any centres are blurred bqs (e.g. from
c    charmm) at the point the charmm geometry
c    is imported
c
c     write(6,*)'import geom'

      do i = 1, nat
         atomss(i)=isubst(zaname(i))
      enddo
      ierr0=CD_import_geom(nat,ne,c,atomss)
      if(ierr0.gt.0)then
         call caserr2('cd_import_geom failed')
      endif
         if (obas .and. .not.oterm) go to 20

      end if
 570  if (ifau.ne.0) toler = toler/toang(1)

      if (lcoord.le.0) then
c
c     ----- write atomic cooordinates in the first record of daf -----
c
         idum = 0
         call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
      else
c
c     ----- read actual coordinates -----
c
         call rdrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
      end if
c
      call revise
      call revind
      iblkd = iblkdu
      ifild = idaf
      ndump4 = idaf
c...  ncoorb seems to be meant as # orbitals left after reduction
c...  mainly used in mp2/mp3 (beware of resetting this variable)
      ncoorb = num
c...  norbt is max # active orbitals
c     norbt = num
      norbt = newbas0
      nsa4= num
      nsb = num
      lenb4 = ikyp(num)
      do 580 loop = 1 , maxorb
         ilifq(loop) = (loop-1)*num
         ilifm(loop) = ilifq(loop)
 580  continue
      do 590 loop = 1 , 10
         title(loop) = ztitle(loop)
 590  continue
      groupy = zgroup
      groupx = zgroup
      jrec = 0
      if (.not.(oterm)) go to 30
      go to 2730
c     ----- now process run option directives
c
c ###############################################################
c   class-3 directives
c ###############################################################
c
c     cards
c
 640  call inpa(buff)
      do 650 i = 1 , ntypr
         if (buff.eq.zpunch(i)) opunch(i) = .true.
 650  continue
      if (opunch(11)) opunch(5) = .true.
      if (jrec.ge.jump) go to 20
      go to 640
c
c    debug directive
c
c
 660  call inpa(buff)
      if (buff.eq.'all') then
       do loop = 1, 40
        odebug(loop) = .true.
       enddo
      else if (buff.eq.'direct') then
c
c**** debug directive gives extra output for direct-scf
c
c    1   regurgitation of input data by direct scf code
c    2   integral test data after everry iteration
c    3   switches off symmetry after data has been read in
c    4   prints out 1-electron integrals
c    5   prints out the guess generated density matrix
c    6   prints out fock and density matrix at each iteration
c    7   prints out integral statistics
c    8   print out separate contributions to the gradient
c
         do 670 loop = 1 , jump - 2
            call inpi(k)
            if (k.lt.9) odbg(k) = .true.
 670     continue
      else
c    1   vcd
c    2   outa
c    3   densd/denpol
c    4   fd2c/fd2d/fd1
c    5   atoms/ptgrp/trmat ( print=1 )
c    6   chfcon ( print=-7 )
c    7   ddjk
c    8   construction of rhs
c    9   extra scf output
c   10   transformation table in symmos
c   11   4-index
c   12   jandk
c   13   jkder
c   14   mogues
c   15  dpdint
c   16   whatsl
c   17   stvder
c   18   print all 1-electron gradient (print=-3)
c   19   standv ( print=3 )
c   20 fockd2/hsymd
c   21 fockabd/polder
c   30 section/file handling
c   31 timing statistics
c   32 full-ci
c   33 ccsd
c   34 mcscf
c   38 uhfmp2
c   39 memory allocation
c   40 direct-ci
c
         jrec = jrec - 1
         do 680 loop = 1 , jump - 1
            call inpi(k)
            if (k.le.ntypr) odebug(k) = .true.
 680     continue
      end if
      go to 20
c
c     gtens
 690  lgten = .true.
      go to 20
c
c  vcd
c
 700  lvcd = .true.
      go to 20
c
c     field
c
 710  call inpa(buff)
      if (buff.eq.zoff) lfield = .false.
      if (buff.ne.zoff) then
         jrec = jrec - 1
         call inpf(fx)
         call inpf(fy)
         call inpf(fz)
         lfield = .true.
      end if
      go to 20
c
c  imptprint
c
 730  call inpi(i)
      if (i.eq.0) go to 20
      if (i.lt.1 .or. i.gt.15) call caserr2('error in impt input')
      call inpi(iopp(i))
      go to 730
c
c     open,oscf,grhf
c
 740  if (buff.eq.zscf(6)) then
         lmcscf = .true.
         lmcdat = .true.
         call inpa(buff1)
         if (buff1.eq.yes .or. buff1.eq.on .or. buff1.eq.zblank)
     +       lmcscf = .true.
         if (buff1.eq.zno .or. buff1.eq.zoff) lmcscf = .false.
         write (iwr,*) ' lmcscf in start ' , lmcscf
      else
         if (scftyp.eq.buff) go to 20
         scftyp = buff
         call reinit
         write (iwr,*) 'setting up new scf type '
      end if
      if (scftyp.eq.words(16)) scftyp = zopen
      oss = .false.
      go to 20
c
c     ncoorb
c
 750  call inpi(ncoorb)
      go to 20
c
c     keepfock
c
 760  fkder = zdd4(20)
      call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) fkder = zblank
      go to 20
c
c  jacobi
c
 770  jaco = .true.
 775  call inpa4(ytest)
      if (ytest.eq.'off' .or. ytest.eq.'no') then
       jaco = .false.
      else if (ytest.eq.'nosy') then
       symm_diag = .false.
      else if (ytest.eq.' ') then
       go to 20
      endif
      go to 775
c
c     contract
 780  lcontr = .true.
      go to 20
c
c     weights 
c
 790  call input
      if (nat.lt.1) call
     & gamerr('weights directive must follow geometry specification',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)

      call inpa(buff)

      if (buff.eq.stop) then

        go to 20

      elseif (buff(1:4).eq.'isot')then
c
c add a new isotope (will also edit the existing table)
c
         call inpa4(ybuf)
         iz = isubst(ybuf)
         call inpf(wght)
         call mass_newiso(iz,wght)

      elseif (buff(1:4).eq.'vect')then
c
c  read the atomic weight vector
c
         call inpa4(ybuf)
         if(ybuf .eq. 'crea')then
            imassv = mass_rdvec(-1)
         else
            jrec = jrec - 1
            idumf = mass_rdvec(1)
         endif
      
      elseif (buff(1:4).eq.'subs')then
c
c modify substitution pattern for mass vector 1
c
         call inpa(ztext)

         if(ztext(1:4) .eq. 'vect')then
c
c optional mass vector specification
c
            call inpi(imassv)
            call inpa(ztext)
         endif
c
c check an integer or a label has been entered
c
         call intval(ztext,iat,onumb)

         call inpa4(ybuf)

         if(ybuf .eq. 'isot')then
            call inpi(iso)
            if(onumb)then
               call mass_editvec(imassv,iat,iso,dum)
            else
               do iat = 1,nat
                  if(ztext .eq. ztag(iat))
     &                 call mass_editvec(imassv,iat,iso,dum)
               enddo
            endif

         else if(ybuf .eq. 'mass')then

            if(imassv .ne. 1)call 
     &gamerr('user defined atomic masses only valid for mass vector 1',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)

            call inpf(wght)

            if(onumb)then
               call mass_editvec(1,iat,-2,wght)
            else
               do iat = 1,nat
                  if(ztext .eq. ztag(iat))
     &                 call mass_editvec(1,iat,-2,wght)
               enddo
            endif
         else
            call 
     &  gamerr('bad substitute keyword in weights specification',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
         endif
      else
c
c original weights style input
c
c
c  <iat> <mass> 
c  <element> <iat> <imass>
c
         call intval(buff,iat,onumb)

         if(onumb) then
           call inpf(wwz)
            if(imassv .ne. 1)call 
     &gamerr('user defined atomic masses only valid for mass vector 1',
     &        ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)

           call mass_editvec(1,iat,-2,wwz)
         else
             call inpi(iat)
             i1 = isubst(buff)
             i2 = isubst(ztag(iat))
             if(i1 .ne. i2)then
                 write(iwr,*)'atom=',iat,'z=',i2,'weights z=',i1
                 call
     &  gamerr('weights directive invalid - atom is of wrong type',
     &  ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
             endif
             call inpi(imas)
             call mass_editvec(imassv,iat,imas,0.0d0)
         endif
      end if

      goto 790
c
c     canonical 
c
 800  call inpf(canna)
      call inpf(cannb)
      call inpf(cannc)
      go to 20
c
c     open shell singlet
c
 810  scftyp = zoscf
      oss = .true.
      mul = 1
      go to 20
c
c     sections
c
 820  call secion
      go to 20
c
c      density
c
 830  call inpi(isecd1)
      call inpi(itypd1)
      if (jump.ne.3) then
         call inpi(isecd2)
         call inpi(itypd2)
         if (jump.ne.5) then
            call inpi(isecd3)
            call inpi(itypd3)
         end if
      end if
      go to 20
c
c     mos
c
 840  lmos = .true.
      call inpi(isecv)
      call inpi(itypv)
      go to 20
c     order
 850  call inpi(norder)
      go to 20
c     epstein nesbet denominators
 860  ldenom = .true.
      go to 20
c     ignore intermolecular overlap
 870  ignore = .true.
      go to 20
c
c     salvage dumpfile directory block ( not possible with
c     some versions of i/o routines )
c
 880  write (iwr,6060)
 890  call input
      call inpa(buff)
      if (buff.eq.stop) then
         orevis(1) = .false.
         call clenup
         go to 20
      else
         jrec = jrec - 1
         call inpi(iseci)
         call inpi(icla)
         call inpi(iposi)
         iposi = iposi - ibl3d
         call inpi(length)
         iblkla = iblkla + length
         call inpi(lds(iseci))
         call qsector(iseci,iposi,iclass,length,'put')
         go to 890
      end if
c
c    ----- shells  ---- read in grhf data
c
 900  if (num.eq.0) call caserr2(
     +       ' shells directive must follow atoms specification ')
      call grhfin(.true.,.false.,iwr)
      go to 20
c
c    orbcore
c
 910  call inpi(ncore)
cj    cicv=ncore.gt.0.or.nvr.gt.0
cj    ncact=ncoorb-ncore-nvr
cj    nupact=ncore+ncact
      call inpi(nvr)
      go to 20
c
c     gauge
c
 920  call inpf(gxx)
      call inpf(gyy)
      call inpf(gzz)
      orgin = .true.
      ggx(1) = gxx
      ggx(2) = gyy
      ggx(3) = gzz
      go to 20
c
c     copy nblocks from edn to edm
c
 930  call copyfi
      go to 20
c
c     reorder
c
 940  call reor(core,jump,num,iint)
      go to 20
c
c     restrict
c
 950  call inpi(iscftp)
      if (jump.eq.1) iscftp = 1
      if (iscftp.eq.0) iscftp = 99
      go to 20
c
c     reset certain restart parameters
 960  call inpi(irest2)
      go to 20
 970  call inpi(irest3)
      go to 20
 980  call inpi(irest4)
      go to 20
 990  call inpi(irest5)
      go to 20
 1000 call inpi(irestp)
      go to 20
 1010 call inpi(irest1)
      go to 20
 1020 call inpi(irest6)
      if (irest6.eq.0) iochf(18) = 0
      go to 20
c
c     mp2 <on / off >
c
 1030 call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) then
         mp2 = .false.
         call reinit
         write (iwr,*) ' mp2 setting cleared'
         go to 20
      end if
      if (.not.(mp2)) then
         mp2 = .true.
         call reinit
      end if
      go to 20
c
c    symmetry order mos for ci/mcscf calculations
c
 1040 ordmo = .true.
      go to 20
c
c
c     seczero -  wipe out section on dumpfile
c
 1050 call inpi(jsect)
      call qsector(jsect,0,0,0,'put')
      orevis(1) = .false.
      call revind
      go to 20
c
c     mp3 <on / off >
c
 1060 call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) then
         mp3 = .false.
         call reinit
         write (iwr,*) ' mp3 setting cleared'
         go to 20
      end if
      if (.not.(mp3)) then
         mp3 = .true.
         call reinit
      end if
      go to 20
cj
cj   spdoc - doubly occupied orbitals which change occupation in the
cj   multi-reference ci configuration list
cj
 1070 cicx = .true.
 1080 call inpi(ix)
      if (ix.eq.0) go to 20
cj    nspj(ix)=1
      nspj(ix) = ix
      write (iwr,*) ' special core orbital ' , ix
      go to 1080
cj
cj    spuoc - unoccupied orbitals which change occupation in the
cj    multi-reference ci configuration list
cj
 1090 cicx = .true.
 1100 call inpi(ix)
      if (ix.eq.0) go to 20
cj    nspj(ix)=2
      nspj(ix) = num + ix
      write (iwr,*) ' special virtual orbital ' , ix
      go to 1100
c   helfeygr
c
 1110 hfgr = .true.
      call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) then
         hfgr = .false.
      end if
      go to 20
c
c     to reset value of mprest
c
 1120 call inpi(mprest)
      go to 20
c
c     pwftol
c
 1130 call inpf(pwftol)
      go to 20
c
c     fcm - force constant from external file
c
 1150 call fcmin('input')
      go to 20
ci
ci      pump2
ci
 1160 call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) then
         pump2 = .false.
         call reinit
         write (iwr,*) ' pump2 setting cleared'
         go to 20
      end if
      if (.not.(pump2)) then
         pump2 = .true.
         mp2 = .true.
         scftyp = zopen
         call reinit
      end if
      go to 20
c
c     ci
c
 1170 if (jump.eq.2) then
         call inpa(buff)
         if (buff.eq.zoff .or. buff.eq.zno) then
            lci = .false.
            cigr = .false.
            cicv = .false.
            cicx = .false.
            ncore = 0
            nvr = 0
            nupact = 0
            ncact = 0
            write (iwr,*) ' ci setting cleared '
            call reinit
         end if
      end if
      if (.not.(lci)) then
         lci = .true.
         call reinit
         if (lci) write (iwr,*) ' lci is true '
      end if
      go to 20
c
c     icflag
c
 1180 call inpi(icflag)
      go to 20
cj
cj   frozen specifies the core and virtual orbitals for the ci part of
cj   an mcscf/ci calculation
cj
cj    format :    frozen
cj                core 2 symm 3
cj                core 1 symm 4   ( for 2fzc3 1fzc4)
cj                virtual 1 symm 4 ( fzv4)
cj                end
cj
 1190 call input
      write (iwr,*) ' in frozen '
      call inpa(buff)
      if (buff.eq.zcore) then
         call inpi(ino1)
         write (iwr,*) ' ino1 = ' , ino1
         call inpa(buff)
         call inpi(ino2)
         ncorea(ino2) = ino1
         go to 1190
      else
         if (buff.eq.stop) go to 20
         call inpi(ino1)
         call inpa(buff)
         call inpi(ino2)
         nvia(ino2) = ino1
         go to 1190
      end if
c
c     minfc - minimum value of eigenvalue of f c m in optimisation
c
 1200 call inpf(tiny)
      go to 20
cj
cj   adaptint - to adapt integrals before calculation of scf
cj
 1210 ladp = .true.
      call inpa(buff)
      if (buff.eq.zoff .or. buff.eq.zno) ladp = .false.
      go to 20
c
c     set logical cpf
c
 1220 if (jump.eq.2) then
         call inpa(buff)
         if (buff.eq.zoff .or. buff.eq.zno) then
            lci = .false.
            lcpf = .false.
            cigr = .false.
            cicv = .false.
            cicx = .false.
            ncore = 0
            nvr = 0
            nupact = 0
            ncact = 0
            write (iwr,*) ' cpf setting cleared '
            call reinit
         end if
      end if
      if (.not.(lcpf)) then
         lcpf = .true.
         lci = .true.
         call reinit
         write (iwr,*) ' lcpf is true '
      end if
      go to 20
cj
cj   optorb - to optimise orbitals for scf based correlation methods
cj
 1230 loptor = .true.
      call inpa(buff)
      if (buff.eq.zno .or. buff.eq.zoff) loptor = .false.
      go to 20
c
c     mp2 passes - this saves disc space on the sort file ( short term)
c
 1240 call inpi(mppas)
      go to 20
c
c     mpflag
c
 1250 call inpi(mpflag)
      go to 20
c
c    nonuke
c
 1260 lnucl = .false.
      go to 20
c
c
c     frequency
c
 1270 if (ipol.eq.0) call caserr2(
     + 'frequency directive should follow polarisability')
      call frecky
      go to 20
c
c     onelec
c
 1280 call onelec(ione)
      go to 20
c
c     cutoff
c
 1290 call inpi(iacc)
      icut = iacc
      cutomp = 10.0d0**(-iacc)
      go to 20
c
c    chfconv
c
 1300 call inpi(iconvv)
      go to 20
c
c    perturb
c
 1310 if (ipol.eq.0) call caserr2(
     + ' perturbations directive should follow polarisability')
c
      call genpo
      go to 20
c
c           mixed   c6   coefficients
c
 1320 if (ipol.eq.0) call caserr2(
     +  'mixed directive should follow polarisability')
      call inpa(filb)
      call inpi(iblb)
      do 1330 i = 1 , maxlfn
         if (filb(1:4).eq.yed(i)) ifilb = i
 1330 continue
      write (iwr,6070) filb , iblb
      lenw = 31 + lenint(192)
      lenwx = 50
      call rdchr(qnames,lenwx,iblb,ifilb)
      call reads(frqq,lenw,ifilb)
      iblkb = iblb
      ifileb = ifilb
      irestp = 4
      omixed = .true.
      go to 20
c
c     polarisa / polariza
c
 1340 ipol = 1
      runtyp = buff
      key = i - nopty
      go to 20
c
c
c     hfscf, gradient, etc.
c
 1350 runtyp = buff
      key = i - nopty
      go to 20
c
c    magnetise
c
c
 1360 write (iwr,6080)
      runtyp = buff
      key = i - nopty
      npstar = 0
      np = 6
      npa = 6
      do 1370 i = 7 , 50
         opskip(i) = .true.
         ipsec(i) = 0
 1370 continue
      do 1380 i = 1 , 6
         opskip(i) = .false.
         pnames(i) = pnamc(i+34)
         ipsec(i) = 286 + i
 1380 continue
      ipsec(4) = 284
      ipsec(5) = 285
      ipsec(6) = 286
      npole = 0
      nfreq = 0
      go to 20
c -----------------------------------------------------
c
c
c      adpat, transform and dispersion
c
 1390 ladapt = .true.
      if (nopk.eq.0) go to 1410
      runtyp = buff
c     lword = maxq
      lsym = .true.
      key = i - nopty
      write (iwr,6090)
      ncoorb = num
      isecvv = isect(51)
      itypvv = 0
      go to 20
 1400 write (iwr,6100)
      if (nopk.ne.0) then
         runtyp = buff
         call inpa(buff)
         if (buff.eq.'disc') odisc = .true.
c        lword = maxq
         key = i - nopty
         isecvv = isect(8)
         itypvv = 8
         go to 20
      end if
c
 1410 write (iwr,6110)
      call caserr2('use supermatrix off')
      go to 20
c ###############################################################
c   class-2 directives
c ###############################################################
c ...
c     simons & jorgensen optimisation directives.
c ...
c ...
c    maxjorg : maximum number of steps allowed
c ...
 1420 call inpi(maxj)
      if (maxj.gt.mit) then
         write (iwr,6120) maxj
         maxj = mit
      end if
      go to 20
c ...
c     recalc : recalculate hessian every irecal points
c     default is to update hessian only.
c ...
 1430 if (jump.lt.2) then
         irecal = 10
      else
         call inpa(ztest)
         if (ztest.ne.zoff) then
            jrec = jrec - 1
            call inpi(irecal)
         else
            irecal = -1
         end if
      end if
      go to 20
c ...
c    powell update. should be default for transition state (check )
c ...
 1440 opowel = .true.
      obfgs = .false.
      obfgx = .false.
      go to 20
c ...
c     bfgs update of hessian. should be default for minimum (check )
c ...
 1450 obfgs = .true.
      opowel = .false.
      obfgx = .false.
      go to 20
c ...
c     bfgs update with safegaurds to ensure retention
c          of positive definite hessian.
c ...
 1460 obfgx = .true.
      opowel = .false.
      obfgs = .false.
      go to 20
c ...
c     rfo on  : take rfo or p-rfo step only ( default)
c     rfo off : take rfo or p-rfo step for "wrong" hessian
c               otherwise take newton-raphson step.
c ...
 1470 if (jump.eq.2) then
         call inpa(ztext)
         if (ztext.eq.zoff) onrfo = .false.
      end if
      go to 20
c ...
c     relaxation of certain cutoffs to control optimisation.
c ...
 1480 ocut = .true.
      go to 20
c ...
c     print option. controls extra printing in optimisation.
c     default is on.
c ...
 1490 if (jump.eq.2) then
         call inpa(ztext)
         if (ztext.eq.zoff) outjor = .false.
      end if
      go to 20
c ...
c     debug printing. default is off.
c ...
 1500 if (jump.eq.2) then
         call inpa(ztext)
         if (ztext.ne.zoff) outdeb = .true.
      end if
      go to 20
c ...
c     minimum allowable magnitude of the hessian eigenvalues.
c ...
 1510 if (jump.gt.2) then
         call inpf(egm)
         call inpi(igm)
         eigmin = egm*10.0d0**iabs(igm)
      else if (jump.eq.2) then
         call inpf(eigmin)
      end if
      go to 20
c ...
c     maximum allowable magnitude of the hessian eigenvalues.
c ...
 1520 if (jump.gt.2) then
         call inpf(egm)
         call inpi(igm)
         eigmax = egm*10.0d0**iabs(igm)
      else if (jump.eq.2) then
         call inpf(eigmax)
      end if
      go to 20
c ...
c     end of simons & jorgensen - specific optimsation directives.
c ...
c
c ----- model
c
 1530 omech = oyes
      odnew = ono
      call mechin
      go to 30
c
c ----- rpa (no longer used)
c
 1540 go to 30
c
c ----- sac-ci (no longer used)
c
 1550 go to 30
c
c ----- mopac
c
 1560 omopac = oyes
      call mopin(core)
      go to 30
c
c ----- mrdci
c
 1570 call inpa4(ytext)
      if(ytext.eq.'new'.or.ytext.eq.'dire') then
       omrdci = ono
       opark = .true.
       call mrdcin2
      else
       omrdci = oyes
c...    for unpack
       oclunp = .true.
       call mrdcin
      endif
      go to 30
c
c ----- fullci
c
 1580 ofulci = oyes
      call fulcin
      go to 30
c
c ----- ccsd
c
 1745 occsd = oyes
      if(ztext.eq.'ccsd(t)') occsdt = oyes
      call ccsdin
      go to 30
c
c ----- qcisd
c
 1746 oqcisd = oyes
      if(ztext.eq.'qcisd(t)') oqcisdt = oyes
      call ccsdin
      go to 30
c
c ----- savef
c
 1747 continue
      call savef
      go to 20
c
c ----- hmat
c
 1748 continue
      call write_fcm
      go to 20
c
c ----- nbo
c
 1749 continue
      inbo = 1
 1751 call input 
      call inpa4(ztext)    
      ytext = ytrunc(ztext)
      k=locatc(ydd(101),limit2,ytext)
      if (k.eq.0) goto 1751
      jrec=jrec-1
      go to 30
c
c ----- dpa
c
 1595 continue  
      idpa = 1
      go to 30
c
c ----- dma
c
 1590 call dmain
      idma = 1
      go to 30
c
c ----- multi/mcscf
c
 1600 call multin(core)
      go to 30
c
c ----- mulliken
c
 1610 imull = 1
 1620 call inpa4(ytest)
      j = locatc(ymull,7,ytest)
      if (j.ne.0) then
         go to (1630,1640,1640,1645,1650,1656,1655) , j
      else
         jrec = jrec - 1
         call active(mulla,mullac,ono,ono)
         go to 1655
      end if
 1630 omulla = oyes
      go to 1620
 1640 omullo = oyes
      go to 1620
 1645 omullg = oyes
      go to 1620
 1650 omullp = oyes
      go to 1620
 1656 omullu = oyes
      go to 1620
 1655 if (omullg) then
        call mulinp(zmul,num)
      endif
      go to 20
c
c ----- i.p. / ip
c
 1660 call inip
      go to 30
c
c ----- onelec / core
c
 1670 jump1 = jump - 1
      if (jump1.gt.0) then
         do 1680 i = 1 , jump1
            call inpa4(ytest)
            if (ytest.eq.yout(2).or.ytest.eq.'prin') then
               norbcp = 1
            else if (ytest.ne.yout(18)) then
               jrec = jrec - 1
               call inpi(ihamc)
            else
               ihamcp = 1
            end if
 1680    continue
      end if
      call active(norbc,norbca,oyes,oyes)
      go to 20
c
c ---- irc calculation
c
 1675 oirc = oyes
      call ircdat(iwr)
      go to 30
c
c ----- active
c
 1690 jump1 = jump - 1
      oactrn = oyes
      if (jump1.gt.0) then
         do 1700 i = 1 , jump1
            call inpa4(ytest)
            if (ytest.ne.yout(2).and.ytest.ne.'prin') then
               jrec = jrec - 1
               call inpi(ionsv4)
            else
               norbtp = 1
            end if
 1700    continue
      end if
      call active(norbt,norbta,ono,oyes)
      go to 20
c
c ----- direct
c
 1710 call ciin(core,.false.)
      odci = .true.
      go to 30
c
c ----- graphics
c       allocate core for user defined grids on cards
c
 1720 iplot = iplot + 1
      if(iplot.gt.1)call caserr2('only one graphics directive allowed')
c
c     allocate available memory
      ibse = igmem_alloc_all(lwrd)
c
      call plotin(core(ibse),lwrd)
c     free memory
      call gmem_free(ibse)
c
      go to 30
c
c ----- potfit (atomic charges fitted to electrostatic grids)
c 
 1725 if(opotf(0))call caserr2('only one potfit directive allowed')
      odumf=opotf(1)
      jump1 = jump - 1
      i = 0
 1723 continue
      call inpa4(ytext)
      i = i+1
      call intval(ytext,j,onumb)
      if(onumb)then
         call potfsc(j)
      else
         if(ytext.eq.'dipo')then
            call inpf(dipx)
            call inpf(dipy)
            call inpf(dipz)
            call potfdp(dipx,dipy,dipz)
            i = i + 3
         else if (ytext.eq.'symm')then
            call potfsy
         else if (ytext.eq.'char')then
            call inpf(char)
            call potfch(char)
            i = i + 1
         else if (ytext.eq.'cuto')then
            call inpf(rcut)
            call potfrc(rcut)
            i = i + 1
         endif
      endif
      if(i.lt.jump1)goto 1723
      goto 30
c
c ----- omitcharge
c
1734  call caserr2('omitcharge data line invalid')
1735  noncc = 0
      omtchg = .true.
c
c blank out arrays
c
      do ii = 1,mxomat
         ztagcc(ii) = ' '
         do jj=1,mxombq
            ztagcq(ii,jj) = ' '
         enddo
      enddo
c
1736  call input
      call inpa(ztext)
      if(ztext.eq.zend)go to 1737
      if(jump.le.1) go to 1734
      noncc = noncc + 1
      ztagcc(noncc) = ztext
      jump1 = jump -1
      if (jump1.gt.mxombq) then
          write(iwr,*)'too many centres omited for', ztagcc(noncc)
          write(iwr,*)'limit is ',mxombq
          go to 1734
      endif
      do 1738 i = 1, jump1
1738  call inpa (ztagcq(noncc,i))
      go to 1736
1737  if(noncc.le.0.or.noncc.gt.mxomat) call caserr2(
     + 'omitcharge data invalid')
      go to 20
c
c
c ----- localize -
c
 1730 local = 1
      jump1 = jump - 1
      ideff = index(char1,'defa')
      if (jump.le.1) then
         call active(nnacta,louta,ono,oyes)
         go to 20
      else if (jump1.le.3) then
         call inpa4(ytext)
         if(ytext.eq.'over'.or.ytext.eq.'pipe') then
          local = 2
          jump1 = jump1 - 1
         else if (ytext.eq.'boys') then
          local = 1
          jump1 = jump1 - 1
         else
          jrec = jrec -1
         endif
         go to 1750
      end if
 1740 if(ztext.eq.'default') then
       go to 20
      else
       call caserr2('invalid parameters specified in local directive')
      endif
 1750 if(jump1.le.0.and.ideff.eq.0) then
       call active(nnacta,louta,ono,oyes)
       go to 20
      endif
      do 1760 i = 1 , jump1
         call inpa(ztestl(i))
 1760 continue
      do 1770 i = 1 , jump1
         ztext = ztestl(i)
         if (ztext.ne.zatext) then
            if (ztext.ne.zbtext) go to 1740
            call active(nactb,loutb,ono,oyes)
         else
            call active(nnacta,louta,ono,oyes)
         end if
 1770 continue
      go to 20
c
c ----- nuclidic
c
 1780 call gamerr('nuclidic directive obsolete - use weights',
     &     ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
c      call input
c      call inpa(ztest)
c      if (ztest.eq.zend) go to 20
c      mcen = mcen + 1
c      jrec = jrec - 1
c      call inpa(zcentg(mcen))
c      call inpf(cenmas(mcen))
c      go to 1780
c
c ----- property
c
 1790 call inpa4(ytest)
      if(ytest.eq.'atom') then
       iprop = -1
       go to 20
      endif
      iprop = 1
 1800 call input
      call inpa(ztest)
      if (ztest.eq.zend) go to 20
      jrec = jrec - 1
      newpro = newpro + 1
      if (newpro.gt.npmx) call caserr2(
     +  'invalid number of properties requested')
      call inpi(itemp(newpro))
      call inpa(ztagc(newpro))
      call inpi(jtemp(newpro))
      go to 1800
c
c ----- centres
c
 1810 call input
      call inpa(ztest)
      if (ztest.eq.zend) go to 20
      nocp = nocp + 1
      if (nocp.gt.ncmx) call caserr2(
     +   'invalid number of additional centres specified')
      jrec = jrec - 1
      call inpf(vlist(nocp,1))
      call inpf(vlist(nocp,2))
      call inpf(vlist(nocp,3))
      call inpa(zlist(nocp))
      go to 1810
c
c     ----- update ----
c
 1820 call inpa4(ytest)
      iupcod = locatc(yuptyp,3,ytest)
      go to 20
c
c     ----- superci ---
c     ----- now allow dynamic switching of method in CASSCF
c
 1830 call inpa4(ytest)
      if (ytest.ne.'dyna') then
         osuped = .false.
         jrec = jrec -1
         call inlist(isort,n20cas,m0,m999)
      else
c...     initialised in mainci.f
 1831    call inpa4(ytest)
         osuped = .true.
         if (ytest.eq.'off') then
            osuped = .false.
         else if (ytest.eq.'nr') then
            call inpf(swnr)
         else if (ytest.eq.'simu') then
            call inpf(swsimu)
         else if (ytest.eq.'hess') then
            call inpi(nhesnr)
         else if (ytest.eq.'end'.or.ytest.eq.' ') then
      go to 20
         end if
         go to 1831
      end if
      go to 20
c
c     ----- pople ----
c
 1840 ipople = 1
      go to 20
c
c ----- jkop ---
c
 1850 if (jump.ge.2) then
         do 1860 j = 2 , jump
            call inpi(itest)
            oto(itest) = oyes
 1860    continue
      end if
      go to 20
c
c
c    =orbsym=
c
 1870 nsymo = jump - 1
      do 1880 i = 1 , nsymo
         call inpi(iorbsy(i))
 1880 continue
      go to 20
c
c... =noci=
c
 1890 call inlist(noci,n20cas,m1,m999)
      go to 20
c
c... =state=
c
 1900 call inpi(kkk)
      if (kkk.le.0 .or. kkk.gt.5) call caserr2(
     +    'invalid no. of states specified in state directive.')
      if (oactiv .and. kkk.gt.1)
     +    call caserr2('state directive must preceed config.')
      i = jump - 2
      if (i.ne.0) then
         nwt = i
         if (nwt.ne.kkk) call caserr2(
     +                      'invalid no. of state weightings presented.'
     +                      )
         do 1910 i = 1 , nwt
            call inpf(weight(i))
 1910    continue
         tot = 0.0d0
         do 1920 i = 1 , nwt
            tot = tot + weight(i)
 1920    continue
         do 1930 i = 1 , nwt
            weight(i) = weight(i)/tot
 1930    continue
      end if
      go to 20
c
c... =simul=
c
 1940 ifsimu = 1
      call setsto(n20cas,0,isimul)
      call inlist(isimul,n20cas,m1,m999)
      go to 20
c
c... =sort=
c
 1950 call inpi(nsz)
      if (nsz.lt.1 .or. nsz.gt.10) nsz = 10
      ns2 = nsz
      go to 20
c
c.... =augment=
c
 1960 idef = 0
      if (jump.lt.2) idef = 1
      call setsto(n20cas,idef,iam)
      if (jump.ge.2) then
         do 1970 i = 2 , jump
            call inpi(j)
            iam(j) = 1
 1970    continue
      end if
      go to 20
 1980 j = 1
c
c ... =onepdm=
c
 1990 if (j.ge.jump) go to 20
      call inpa4(ytext)
      if (ytext.eq.yout(19)) then
         call inpi(isecdm)
         if (isecdm.gt.0 .or. isecdm.le.204) then
            j = j + 2
            go to 1990
         else
            call caserr2('invalid dumpfile section specified')
         end if
      end if
 2000 if (ytext.ne.yout(20)) then
         j = j + 2
         go to 1990
      else
         call inpi(isecda)
         if (isecda.le.0 .or. isecda.gt.204) then
            call caserr2('invalid dumpfile section specified')
            go to 2000
         else
            j = j + 2
            go to 1990
         end if
      end if
c
c ... =lagr=
c
 2010 j = 1
 2020 if (j.eq.jump) go to 20
      call inpa4(ytext)
      if (ytext.eq.yout(19)) then
         call inpi(iseclm)
         if (iseclm.le.0 .or. iseclm.gt.204) then
            call caserr2('invalid dumpfile section specified')
            go to 2000
         else
            j = j + 2
            go to 2020
         end if
      else if (ytext.eq.yout(20)) then
         call inpi(isecla)
         if (isecla.le.0 .or. isecla.gt.204) then
            call caserr2('invalid dumpfile section specified')
            go to 2000
         else
            j = j + 2
            go to 2020
         end if
      else if (ytext.ne.ygvb) then
         j = j + 2
         go to 2020
      else
         call inpi(iseclg)
         if (iseclg.le.0 .or. iseclg.gt.204) then
            call caserr2('invalid dumpfile section specified')
            go to 2000
         else
            j = j + 2
            go to 2020
         end if
      end if
c
c... =hessian= directive
c
 2030 call inlist(isort,n20cas,m2,m999)
      go to 20
c
c... =civec= directive
c
 2040 call inpf(ciprnt)
      ocidu = oyes
      go to 20
c
c... =canon= cirective
c
 2050 call inpi(icanon)
      go to 20
c
c... =newton= directive
c
 2060 call inlist(isort,n20cas,m1,m2)
      go to 20
c
c... =linear= directive
c
 2070 call input
      ifdeg = 1
 2080 call inpi(i)
      call inpi(j)
      ilifc(jdeg) = i
      ilifc(jdeg+1) = j
      jdeg = jdeg + 2
      idegen(i) = j
      idegen(j) = i
      call input
      call inpa4(ytest)
      if (ytest.eq.yend) go to 20
      ophase = oyes
      jrec = jrec - 1
      go to 2080
c
c... =ciacc= directive
c
 2090 do 2100 i = 1 , 10
         call inpf(fudgit(i))
         if (fudgit(i).lt.1.0d-10) fudgit(i) = fudgit(i-1)
 2100 continue
      if (fudgit(1).lt.1.0d-10)
     +    call caserr2('first ci accuracy too small or zero.')
      ojust = ono
      go to 20
c
c... =supacc= directive
c
 2110 call inpf(accblb)
      go to 20
c
c... =tracc=
c
 2120 call inpi(iacc4)
      go to 20
c
c... =print=
c
 2130 k = jump - 1
      if (k.gt.0) then
         do 2150 i = 1 , k
            call inpa4(ytest)
            do 2140 j = 1 , 16
               oprint(j) = oprint(j) .or. ytest.eq.yout(j)
 2140       continue
            if (ytest.eq.yout(17)) call inpi(nprinc)
 2150    continue
      end if
      write (iwr,6130)
      go to 20
 2160 call inpi(npas41)
      call inpi(npas42)
      if (npas41.lt.1 .or. npas42.lt.1)
     +    call caserr2('npass1 or npass2 must be positive')
      odisc = .true.
      go to 20
c
c... =cipass= directive
c
 2170 call inpi(ipassc)
      write (iwr,6140) ipassc
      go to 20
c
c
c ----- config --  print (pass ipass) nodrt/bypass nosort
c
 2180 if (oactiv .or. zscftp.ne.zscf(5)) then
         call caserr2(
     +    'faulty ordering of CONFIG directive' )
      else
         oactiv = oyes
         call incas1
      end if
      go to 20
c
c... =stop= directive
c
 2190 ostop = oyes
      ostopr = oyes
      go to 2720
c
c    ***** open **** read in open shell data
c
 2200 open = .true.
      if (ouhf) then
c     OPEN will be ignored in UHF calculation
       open = .false.
       go to 20
      endif
      if (.not.ogrhf) then
         ogvb = oyes
         zscftp = zscf(3)
      end if
      if (ifolow.eq.0) lock = 4
      jumper = jump - 1
      do 2210 i = 1 , jumper
         call inpa(zstate)
         j = locatc(zstats,18,zstate)
         if (j.ne.0) go to 2220
         if (i.eq.jumper) zstate = zblank
 2210 continue
      go to 2230
 2220 jumper = jumper - 1
 2230 nset = jumper/2
      if ((nset*2).ne.jumper)
     +    call caserr2('invalid parameters specified in open directive')
      nseto = nset
      jrec = 1
      do 2240 i = 1 , nset
         call inpi(no(i))
         call inpi(noe(i))
         nopen = nopen + no(i)
         nope = nope + noe(i)
         if (i.ne.1) then
            if (no(i).eq.2 .and. no(i-1).eq.2
     +           .and. zstate.ne.zstats(18) ) nseto = nseto + 2
         end if
 2240 continue
      go to 20
c
c     ----- diis (p.pulay, chem.phys.letts 1981)
c
 2250 odiis = oyes
 2260 call inpa4(ytest)
      i = locatc(ydiis,8,ytest)
      if (i.ne.0) then
         go to (2270,2280,2290,2291,2301,2302,2303,2300) , i
 2270    odiis = ono
         odiisn = ono
         go to 20
 2280    opdiis = oyes
         go to 2260
 2290    call inpf(accdi1)
         if (accdi1.le.0.0d0) accdi1 = 0.01d0
      else
         i = locatc(yed,maxlfn,ytest)
         if (i.ne.0) then
            numdis = i
            call inpi(ibldis)
            if (ibldis.le.0) ibldis = 1
         end if
      end if
      go to 2260
 2291    call inpi(idiisf)
         if(idiisf.le.0) idiisf = 3
      go to 2260
 2301    odynamic = .true.
      go to 2260
 2302    olevd = .true.
      go to 2260
 2303    otestd = .true.
      go to 2260
 2300 mconv = 5
      go to 20
c
c     ----- value ----
c
 2310 call inpf(accin(1))
      go to 20
c
c     ----- xtol ----
c
 2320 call inpf(accin(2))
      rcvopt = accin(2)
      go to 20
c
c     ----- step ----
c
 2330 call inpf(accin(3))
      go to 20
c
c     ----- lsearch ----
c
 2340 call inpi(lintyp)
      call inpi(lqstyp)
      go to 20
c
c     ----- tolmax ----
c
 2350 call inpf(accin(4))
      go to 20
c
c     ----- tolstp ----
c
 2360 call inpf(accin(5))
      go to 20
c
c     ----- tanstp ----
c
 2370 call inpf(accin(6))
      go to 20
c
c     ----- minmax ----
c
 2380 call inpa4(ytest)
      if (ytest.eq.'revi') then
         call inpi(iterch)
         if (iterch.le.0 .or. iterch.gt.mit) iterch = mit
      else
         jrec = jrec - 1
      end if
      call inpi(jm)
      call inpi(mnls)
      if (jm.le.0 .or. jm.gt.mit) jm = mit
      if (mnls.le.0 .or. mnls.gt.mit) mnls = mit
      go to 20
c
c     ----- maxcyc ----
c
 2390 call inpi(maxcyc)
      call inpa4(ytest)
      if(ytest.eq.'igno') omaxcyc = .true.
      go to 20
c
c     ----- thre ----
c
 2400 call inpi(nconv)
      onconv = oyes
      if(.not.oreact) then
        if (nconv.gt.9) nconv = 9
      endif
      if(jump.gt.2) then
       call inpa4(ytext)
       if (ytext.eq.'prin') then
        optester = .true.
       else
        jrec = jrec -1
        call inpi(icnv)
        if(icnv.le.0.or.icnv.gt.9) icnv = 5
       endif
      endif
      go to 20
c
c     ----- npunch ----
c
 2410 call inpi(npunch)
      go to 20
c
c     ----- damp ----
c
 2420 call inpf(dmpcut)
      go to 20
c
c     ----- runtyp ----
c
 2430 orun = .true.
      call rntype(jrun)
      if (jrun.eq.42) then
       go to 30
      else
       go to 20
      endif
c
c ... swap
c
 2440 call inpa(ztext)
      iswap = 0
      ibase = 0
      if (ztext.eq.zbeta) ibase = 40
 2450 call input
      call inpa(ztest)
      if (ztest.eq.zend) then
         if (iswap.gt.20) call caserr2(
     +     'invalid number of orbital pairs in swap directive')
         if (ztext.ne.zbeta) nswapa = iswap
         if (ztext.eq.zbeta) nswapb = iswap
         go to 20
      else
         jrec = jrec - 1
         call inpi(moa)
         call inpi(mob)
         if (moa.le.0 .or. mob.le.0 .or. moa.gt.maxorb .or.
     +       mob.gt.maxorb) call caserr2(
     +           'invalid orbital specified in swap directive')
         iswap = iswap + 1
         nswap(ibase+1) = moa
         nswap(ibase+2) = mob
         ibase = ibase + 2
         go to 2450
      end if
c
c     ----- vectors ----
c
 2460 call invect(orest)
      oinv = oyes
      go to 20
c
c     ----- lock ----
c
 2470 call inpa(ztext)
      if (ztext.ne.zoff) then
         lock = 4
         llock = 1
         ifolow = 0
         write (iwr,6150)
      else
         ifolow = 1
         lock = 2
         llock = -1
      end if
      go to 20
c
c     ----- conv ----
c
 2480 call inpi(mconv)
      if (mconv.lt.0 .or. mconv.gt.15)
     +        call caserr2('invalid conv option requested')
      odiis = ono
      odiisn = ono
      go to 20
c
c     ----- natorb (for unocas) -----
c
 2485 call inpa(ztext)
      if (ztext.eq.'spin') then
         call inpi(iunosp)
 2486    call inpa(ztext)
         if (ztext.eq.zblank) go to 20
         if (ztext(1:4).eq.'prin') iunspp = 2
         if (ztext(1:4).eq.'anni') oanil = .true.
         go to 2486
      else
         jrec = jrec - 1
         call inpi(iuno)
 2487    call inpa(ztext)
         if (ztext.eq.zblank) go to 20
         if (ztext(1:4).eq.'prin') iunopp = 2
         if (ztext(1:4).eq.'anni') oanil = .true.
         if (ztext(1:4).eq.'cano') call incanon 
         go to 2487
      end if
      go to 20
c
c     ----- level ----
c
 2490 if(.not.oscf. and. .not.open) then
c     neither scftyp nor open have been specified yet
c     must drive off defaults
c
       call inpf(gapa1)
       olevel = oyes
       if (na.eq.nb) then
          call inpi(ibrk)
          call inpf(gapa2)
        else
          call inpf(gapb1)
          call inpi(ibrk)
          call inpf(gapa2)
          call inpf(gapb2)
       endif
       call inpa(ztext)
       if (ztext(1:4).eq.'diis') olevd = .true.
       if (ibrk.le.0) ibrk = 999
      else
       call inpf(gapa1)
       olevel = oyes
       if (zscftp.ne.zscf(5) .and. zscftp.ne.zscf(6)
     &     .and. zscftp.ne.zscf(11)) then
          if (zscftp.ne.zscf(1)) call inpf(gapb1)
          call inpi(ibrk)
          call inpf(gapa2)
          if (zscftp.ne.zscf(1)) call inpf(gapb2)
          call inpa(ztext)
          if (ztext(1:4).eq.'diis') olevd = .true.
          if (ibrk.le.0) ibrk = 999
       end if
      end if
      go to 20
c
c     ----- scftyp ----
c
 2500 oscf = .true.
      mp2 = .false.
      mp3 = .false.
      lmcscf = .false.
      lci = .false.
      mcgr = .false.
      cigr = .false.
      omas = .false.
      call inpa(ztest)
 2510 if (odscf .and. (ztest.ne.zscf(1).and.ztest.ne.zscf(2)
     +          .and.  ztest.ne.zscf(3).and.ztest.ne.zscf(4)))
     +    call caserr2('direct scf only for rhf, uhf and gvb')
      if (mp2 .and. (ztest.ne.zscf(1) .and. ztest.ne.zscf(2)))
     +    call caserr2('mp2 only for closed shell rhf and uhf')
      if (mp3 .and. (ztest.ne.zscf(1) .and. ztest.ne.zscf(2)))
     +    call caserr2('mp3 only for closed shell rhf and uhf')
      ii = locatc(zscf,12,ztest)
      if (ii.eq.0) then
         call caserr2('unrecognised scftype option requested')
      end if
      zscftp = zscf(ii)
      go to (20,2535,2530,2540,2585,2570,2560,2590,2610,2620,
     +       2575,2595) , ii
 2535 zscftp = zscf(2)
      ouhf = oyes
      go to 20
 2530 zscftp = zscf(3)
      ogvb = oyes
      go to 2550
 2540 zscftp = zscf(4)
      ogrhf = oyes
 2550 if (ifolow.eq.0) lock = 4
      if (.not.old) then
         call inpi(npair)
         if (npair.ne.0 .and. ogrhf) then
            call caserr2('unrecognised scftype option requested')
         end if
      end if
      call inpa4(ytest)
      if (ytest.eq.yrest) orgvb = oyes
      go to 20
 2560 zscftp = zscf(6)
 2570 lmcscf = .true.
 2585 omcscf = .true.
      go to 2580
 2575 ovb = .true.
 2580 nopk = 1
      otri = ono
c
cjvl      see adapt
cjvl      if(otran)call caserr2('symmetry adaption must be used with mcsc
c
c ----- temp. fix until shell ordering in g80//g76 sorted out
c
      if (intg76.ne.0) intg76 = 3
      go to 20
c
c *** direct scf
c
 2590 odscf = oyes
      nopk = 1
      odnew = .false.
      if (jump.gt.2) then
         do 2600 loop = 3 , jump
            call inpa(ztest)
            if (ztest.eq.'new') then
               odnew = oyes
 2603          call inpi(iterdd)
               call inpi(jterdd)
               if(iterdd.gt.0) then
                do 2601 i = iterdd,jterdd
 2601           ofull(i) = .true.
                go to 2603
               else
                go to 2602
               endif
            end if
            if (ztest.eq.'old') then
               odnew = ono
               otran = ono
               go to 2600
            end if
            if (ztest.eq.zscf(9)) then
               mp2 = oyes
               omp2d = oyes
               odnew = oyes
               call inpi(mp2fo)
               call inpi(mp2fv)
               if (intg76.ne.0) intg76 = 3
               call reinit
            end if
            if (ztest.eq.zscf(12)) then
               omas = oyes
               omasd = oyes
               if (intg76.ne.0) intg76 = 3
               call reinit
            end if
            if (ztest.eq.zscf(2).or.ztest.eq.zscf(3).or.
     +          ztest.eq.zscf(4)) then
             go to 2604
            endif
 2600    continue
      else
c     reset ztest for 'scftype direct' open shell to gvb
        if(mul.ne.1) then
         ztest = zscf(3)
        else
         ztest = zscf(1)
        endif
        go to 2604
      end if
 2602 ztest = zscf(1)
 2604  if (odnew) then
        if (omech. or .lpseud.ne.0) then
         if (omp2d) then
            call caserr2(
     +     'Cannot do ECPs + Direct MP2 in serial version')
         endif
         odnew = ono
         otran = ono
        endif
      endif
      go to 2510
c
c masscf
c
 2595 omas = .true.
      zscftp = zscf(12)
      nopk = 1
      if (intg76.ne.0) intg76 = 3
      call reinit
      if(jump.gt.2) then
       call inpa(ztest)
      else
c
cgdf  no unrestricted form of masscf 02.04.08
c     nb. zscf(1)='rhf', zscf='uhf'
c
c      if (na.eq.nb) then
        ztest = zscf(1)
c       else
c       ztest = zscf(2)
c      endif
      endif
      go to 2510
c
c mp2
c
 2610 mp2 = .true.
 2810 nopk = 1
      if (intg76.ne.0) intg76 = 3
      call reinit
      if(jump.gt.2) then
       call inpa(ztest)
      else
       if (na.eq.nb) then
        ztest = zscf(1)
        else
        ztest = zscf(2)
       endif
      endif
      go to 2510
c
c mp3
c
 2620 mp3 = .true.
      go to 2810
c
c *** select   intmag     delta <deltol>   <delfac>
c ***   <itrtol(i)> <tolitr(i)>
c *** end
c
 2630 if (jump.eq.1) then
c
c    apply more aggressive DSCF defaults
c    with delta-density matrix (see initiala)
c    problems on first cycle with non-shell structured
c    basis sets .. increase accuracy for these cases
c    (oshell = .false. from atoms)
c
        dlntol = -dlog(1.0d-8)
        if (oshell) then
         tolitr(1) = -dlog(1.0d-6)
        else
         tolitr(1) = -dlog(1.0d-7)
        endif
        tolitr(2) = -dlog(1.0d-7)
        tolitr(3) = -dlog(1.0d-8)
        itrtol(1) = 3
        itrtol(2) = 5
        itrtol(3) = 999
        deltol = dlog(0.05d0)
        delfac = dlog(0.1d0)
        go to 20
      endif
      call inpa(ztest)
      if (ztest.eq.zoimag) then
         oimag = oyes
         call inpa(ztest)
      end if
      if (ztest.ne.zblank) then
         if (ztest.eq.zdelta) then
            call inpf(deltol)
            deltol = dlog(deltol)
            call inpa(ztest)
            if (ztest.ne.zblank) then
               jrec = jrec - 1
               call inpf(delfac)
               delfac = dlog(delfac)
            end if
            iii = 1
            go to 2640
         end if
         call caserr2('unknown select option.')
      end if
      iii = 1
 2640 call input
      call inpa(ztest)
      if (ztest.eq.zend) go to 20
      if (iii.gt.3) call caserr2('only allowed 3 select values.')
      jrec = jrec - 1
      call inpi(itrtol(iii))
      call inpf(tolitr(iii))
      tolitr(iii) = -dlog(dmax1(tolitr(iii),1.0d-20))
c ***
      iii = iii + 1
      go to 2640
c
c     =prefac= pre-factor integral test thresholds
c
 2650 call input
      call inpa(ztest)
      ipref = 0
      if (ztest.eq.zend) go to 20
      call inpf(pfthr)
      call inpi(iexpf)
      pfthr = pfthr*10.0d0**iexpf
      if (ztest.eq.all) then
         do 2660 i = 1 , 21
            prefac(i) = pfthr
 2660    continue
         go to 2650
      end if
      do 2670 i = 1 , 21
         if (ztest.eq.zchar(i)) ipref = i
 2670 continue
      prefac(ipref) = pfthr
      if (ipref.eq.0) call caserr2('error in prefac directive')
      go to 2650
c
c     accuracy
c     =inner= integral test thresholds
c
 2680 call input
      call inpa(ztest)
      ipref = 0
      if (ztest.eq.zend) go to 20
      call inpf(pfthr)
      call inpi(iexpf)
      pfthr = pfthr*10.0d0**iexpf
      if (ztest.eq.all) then
         do 2690 i = 1 , 21
            tinner(i) = pfthr
 2690    continue
         go to 2680
      end if
      do 2700 i = 1 , 21
         if (ztest.eq.zchar(i)) ipref = i
 2700 continue
      tinner(ipref) = pfthr
      if (ipref.eq.0) call caserr2('error in inner directive')
      go to 2680
c
c ----- =lesstore=
c ----- use io to store fewer n*n matrices
c
 2710 olog1 = .true.
      go to 20
c
c ----- enter
c
 2720 call inpi(moutaa)
      call inpi(moutbb)
c
 2730 maxit = maxcyc
      maxitd = maxit
      nps1 = npas41
      nps2 = npas42
c
c = runtype = specification
c
      if (runtyp.eq.zblank) then
c        no hondo specification .. load up gamess
         if (orun) then
            jrun = locatc(zrunt,50,zruntp)
            if (jrun.eq.0) go to 2800
            irun = irndc(jrun)
            if (irun.ne.0) then
               runtyp = zdd4(irun)
            else
               runtyp = zdd4(1)
            end if
         else
            runtyp = zdd4(1)
         end if
c        no gamess specification  .. load up hondo
      else if (.not.orun) then
         irun = locatc(zdd4,50,runtyp)
         if (irun.eq.0) go to 2800
         jrun = irndg(irun)
         if (jrun.eq.0) go to 2800
         zruntp = zrunt(jrun)
      end if
c
c     set runtype=ci if any of the following have been
c     specified: direct, fullci, mrdci, ccsd, ccsd(t)
c     and runtype is NOT specified
c
      if (.not.ocifp) then
       if(odci.or.ofulci.or.omrdci.or.occsd.or.oqcisd.or.
     +    opark) then
        if(.not.orun) then
         jrun = 3
         zruntp = zrunt(jrun)
        endif
       endif
      endif
c
      if (na.ne.nb .and. (.not.open) .and. (.not.ouhf)
     +   .and. (.not.omcscf) .and. (.not.ovb)) then
c
c ---- setup default open-shell specifications
c ---  open-shell gvb/grhf module
c
         ogvb = oyes
         zscftp = zscf(3)
         if (ifolow.eq.0) lock = 4
         zstate = zblank
         nset = na - nb
         nseto = nset
         do loop = 1 , nset
            no(loop) = 1
            noe(loop) = 1
            nopen = nopen + 1
            nope = nope + 1
         enddo
      end if
c
c = scftype = specification
c
      if (scftyp.eq.zblank) then
c        no hondo specification .. load up gamess
         isscf = locatc(zscf,8,zscftp)
         if (isscf.eq.1) scftyp = zclose
         if (isscf.eq.2) scftyp = 'open'
         if (isscf.eq.3) scftyp = 'grhf'
         if (isscf.eq.4) scftyp = 'grhf'
         if (isscf.eq.6) scftyp = 'mcscf'
c        no gamess specification  .. load up hondo
      else if (scftyp.eq.zclose) then
         zscftp = zscf(1)
      else if (scftyp.eq.'open') then
         zscftp = zscf(2)
      else if (scftyp.eq.'oscf') then
         zscftp = zscf(3)
      else if (scftyp.eq.'grhf') then
         zscftp = zscf(4)
      else if (scftyp.eq.'mcscf') then
         zscftp = zscf(6)
      else
         call caserr2('invalid scftype')
      end if
c
c     specify default vectors options if omitted
c

      call defvect(oinv,moutaa,moutbb,orest,oresti)
c
c     optimise=optimize
      if (runtyp.eq.zdd4(30)) runtyp = zdd4(4)
c     polarisa=polariza
      if (runtyp.eq.zdd4(35)) runtyp = zdd4(14)
c     magnetisa=magnetiza
      if (runtyp.eq.zdd4(36)) runtyp = zdd4(17)
c     energy=hfscf
      if (runtyp.eq.zdd4(31)) runtyp = zdd4(1)
      limpt = limpt .or. (key.ge.7 .and. key.le.11)
      lskip = lskip .or. (key.ge.8 .and. key.le.10)
      if (limpt) nopk = 1
      if (runtyp.eq.zdd4(14) .and. fkder.eq.zdd4(20)) then
         ldiag = .true.
         fkder = zblank
      end if
      if (lci .and. runtyp.eq.zdd4(1)) then
         lcontr = .true.
      end if
cj  lopti gets set in master - this keeps ifil2d as ed4
cj
      if (runtyp.ne.zdd4(4) .and. runtyp.ne.zdd4(30)) lopti = .false.
      if (runtyp.ne.zdd4(5)) lforce = .false.
      if (lci .and. (.not.lopti) .and. (.not.lforce)) then
         ncoorb = num
         ncact = ncoorb - ncore - nvr
         nupact = ncact + ncore
         ntot = ncoorb
         ifil2d = 6
         iblk2d = 1
         ifilh = 1
         iblkh = 1
         isecdd = 301
         isecll = 302
         if (runtyp.eq.zdd4(2) .or. runtyp.eq.zdd4(4) .or.
     +       runtyp.eq.zdd4(5)) cigr = .true.
         cicv = cigr .and. (ncore.gt.0 .or. nvr.gt.0)
      end if
c
      if (lmcscf .and. (runtyp.eq.zdd4(2) .or. runtyp.eq.zdd4(4) .or.
     +    runtyp.eq.zdd4(5).or.runtyp.eq.zdd4(18))) mcgr = .true.
      if (pump2 .and. runtyp.ne.zdd4(1))
     +    call caserr2('spin-projected ump2 energies only')
      if (scftyp.eq.zclose .and. mul.ne.1) then
         call caserr2('scf type closed, but multiplicity not 1')
      end if
c
c     error trap for mp2 and mp3
c
      if ((mp2 .or. mp3) .and. (scftyp.ne.zclose .and. scftyp.ne.zopen))
     +    call caserr2('mp2 and mp3 only for closed shell and uhf ')
      if (mp2 .or. mp3 .or. lci .or. omas) then
         do 2750 i = 1 , 8
            if (groupx.eq.grp(i)) then
               igrp = i
               go to 2770
            end if
 2750    continue
 2760    write (iwr,6160)
         call caserr2('modify molecular symmetry and point group')
 2770    if (jaxis.ne.2 .and. igrp.gt.3) go to 2760
      end if
      if (runtyp.eq.zdd4(21)) then
         if (scftyp.eq.zopen) call caserr(
     +  'analytic second derivatives not available:use numerical method'
     +  )
         if (oss) call caserr(
     +    ' analytic second derivatives not available, see grhf option')
      end if
      if (.not.onconv) then
         if (odscf .and. odnew) then
            if (jrun.gt.3) nconv = 6
            go to 2780
         end if
         if (omas .or. mp2 .or. mp3) nconv = max(7,nconv)
         if (lci .or. lmcscf) nconv = max(6,nconv)
         if (jrun.gt.3.and.jrun.ne.42) then
          if (zscftp.eq.zscf(5)) then
            nconv = max(6,nconv)
          else
            nconv = max(7,nconv)
          endif
         endif
      end if
      if ((mp2 .or. mp3 .or. lci .or. lmcscf .or.
     +    zscftp.eq.zscf(5).or.omas). and. nconv.gt.7 ) then
       icut = icut + (nconv-7)
       iacc = iacc + (nconv-7)
       cutomp = 10.0d0**(-icut)
      endif
 2780 nconvc = nconv
c
c ... if scftype gvb specified, switch to grhf
c ... for the specified runtypes
c
      if(zscftp.eq.zscf(3)) then
        if(jrun.gt.20.and.jrun.le.34) then
           ogrhf = oyes
           zscftp = zscf(4)
           if(npair.ne.0)
     +     call caserr2(
     +       'pair-GVB option not available for specified RUNTYPE')
        endif
      endif
c
      if(jrun.gt.20.and.jrun.le.34.and.lpseud.eq.2) call caserr2(
     + 'ecp calculations not available for specified RUNTYPE')
c
      l384 = lenint(3*maxorb)
      lds(isect(44)) = l384
      m44 = 44
      len44 = lensec(lds(isect(44)))
      call secput(isect(44),m44,len44,iblok)
      call wrt3i(ilifq,lds(isect(44))*nav,iblok,ifild)
c
      l100 = 70 + lenint(60)
      lds(isect(103)) = l100
      mtype = 103
      len103 = lensec(lds(isect(103)))
      call secput(isect(103),mtype,len103,iblok) 
      call wrt3(cigr,lds(isect(103)),iblok,ifild)
c
      l27 = 20 + lenint(14)
      lds(isect(13)) = l27
      len13 = lensec(lds(isect(13)))
      m13 = 13
      call secput(isect(13),m13,len13,iblok)
      call wrt3(en,lds(isect(13)),iblok,ifild)
c
      if (zguess.ne.zvect(10)) guess = zguess
c
c  decide if format of integral storage given number of
c  basis functions (this controls all packing/unpacking
c  routines, using the parameters in atmblk) needs revision.
c
      if (num.gt. 255 .or. .not. o255i) then
c
c  supress supermatrix format as default ..
c        if(.not.osupe. and . nopk.ne.1) nopk = 1
c  now supress supermatrix format regardless for high bfn count
c  or simulation through o255i
c
         nopk = 1
         osortp = .false.
         num2e   = 255
         num2ep  = 255
         num2ejk = 170
         o255i   = .false.
         lab1632 = 32
         lab816  = 16
         numlab  = 1020
         numlabp = 510
         numlabjk = 340
         do loop = 1, maxorb
           i4096(loop) = loop * 65536
         enddo
c
c     now exclude runtypes/scftype that will not allow for more
c     than 255 basis functions (this list can be reduced with
c     effort)
c
         if (jrun.eq.1. or. jrun.eq.2. or. jrun.eq.4. or. jrun.eq.5. 
     +   or. jrun.eq.6. or. jrun.eq.7. or. jrun.eq.8. or. jrun.eq.10. 
     +   or. jrun.eq.11.or. jrun.eq.12.
     +   or. jrun.eq.13.
     +   or.jrun.eq.14. or. jrun.eq.15. or. jrun.eq.16. or. jrun.eq.18.
     +   or.jrun.eq.21.
     +   or.jrun.eq.42) then
c     check that the scftype is valid for more than 255 functions
cjvl1            if(zscftp.eq.zscf(11)) then
cjvl1             call caserr2(
cjvl1     +       'SCFTYPE not possible with more than 255 functions')
cjvl1            endif
         else
c        call caserr2(
c    +   'RUNTYPE not possible with more than 255 functions')
         endif
      endif
c
c    now exclude runtypes that will require 2-electron integral files
c    and will abort/fail if direct-SCF processing has been requested
c
      if (jrun.eq.1. or. (jrun.ge.4. and. jrun.le.8). 
     +               or.jrun.eq.10. or. jrun.eq.13. or.jrun.eq.14. 
     +               or.(jrun.ge.16. and. jrun.le.20).
     +or.jrun.eq.29. or. jrun.eq.41. or. jrun.eq.42 ) then
c     fine ... direct-scf is allowed for these
      else
       if(odscf.or.odnew) call caserr2(
     +     'RUNTYPE not compatible with direct-SCF processing')
      endif
c
c
c     Check that the effective number of electrons (i.e. after
c     considering the molecule, the charge, AND any effective core
c     potentials) matches the spin multiplicity.
c
      if(mod(ne,2).eq.mod(mul,2)) then
        call caserr2('inconsistent multiplicity and no. of electrons')
      endif
c
      call start2(core,jrun,jdeg,olevel,oactiv,ocidu,ogvb,
     +            ogrhf,orgvb,odci,obasis,osafe,odft,odftgrad,nav)
c
      odirec(jrun) = .true.
c
      l888 = 748 + lenint(maxorb+26)
      lds(isect(53)) = l888
      len53 = lensec(lds(isect(53)))
      call secput(isect(53),m53,len53,iblok)
      call wrt3i(nact,lds(isect(53))*nav,iblok,ifild)
c
      ifils = num8
      iblks = ibl7la
c     final check on validity of runtype with basis
      call start3(jrun,oreact)
c
 2790 continue
      return
 2800 call caserr2('invalid runtype option')
      return
 6011 format(//10x,
     +'***************************************************************'/
     +10x,
     +'** This job will be terminated prior to dumpfile corruption  **'/
     +10x,
     +'** An attempt has been made to restore vectors from a section *'/
     +10x,
     +'** in a startup job. Please specify the restart directive or  *'/
     +10x,
     +'** redefine the vectors directive e.g. vectors atoms         **'/
     +10x,
     +'***************************************************************'/
     + )
 6010 format (/1x,
     +    '************************************************************'
     +    ,'**********',/1x,
     +    '************************************************************'
     +    ,'**********',/1x,
     +    '*                                                           '
     +    ,'         *'/1x,'*     job  option  complete at  ',f12.2,
     +    ' seconds',a10,' wall  *'/1x,
     +    '*                                                           '
     +    ,'         *'/1x,
     +    '************************************************************'
     +    ,'**********',/1x,
     +    '************************************************************'
     +    ,'**********',/)
 6020 format (' the numer of atoms specified in the mapp directive (',
     +        i4,') is'/
     + ' not equal to the number specified by the zmat/geom directive ('
     + ,i4,').')
 6030 format (//10x,90('*')/10x,'*',88x,'*'/10x,'*',4x,10a8,4x,'*'/10x,
     +        '*',88x,'*'/10x,90('*')//)
 6040 format (1x,21('-')/' this is a startup job'/1x,21('-'))
 6050 format ('    dump file on ',a4,' starting at block',
     +        i6/' scratch file on ',a4,' starting at block',i6)
 6060 format (///1x,'recreating dumpfile directory block'///)
 6070 format (///1x,'mixed dispersion coefficents to be calculated'//1x,
     +        'information for second molecule from file',a5,
     +        ' at block',i5)
 6080 format (//1x,'coupled hartree-fock magnetizability calculation'//)
 6090 format (1x,'symmetry transformation section')
 6100 format (1x,'integral transformation section')
 6110 format (//1x,'4-index transformation section does not work',
     +        ' with the integrals in supermatrix form')
 6120 format (' the maximum number of steps has been reduced from ',i3,
     +        ' to 200')
 6130 format (/)
 6140 format (' ci coefficient passing from cycle',i4)
 6150 format (/' configuration locking requested'/)
 6160 format(
     +10x,'*****************************************************'/
     +10x,'* The molecular point group prohibits use of either *'/
     +10x,'* the requested SCFTYPE or RUNTYPE. Reduce the      *'/
     +10x,'* molecular symmetry by modifying the nuclear TAGs  *'/
     +10x,'*     (see chapter 2 of the User Manual)            *'/
     +10x,'*****************************************************'/)
 6170 format (1x,'restore z-matrix from mopac archive file')
 6180 format (1x,'restore geometry from mopac archive file')
 6185 format(/1x,
     +'The value of degeneracy threshold for symmetry assignment was'/
     +1x,' changed to ',f16.8)
      end
c
      subroutine zora_input
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      character*14 zzstring
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
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      real*8 xocc_z,xoccav_z
      integer nd_z,no_z,nocc_z
      parameter (nocc_z=1024)
      common/occnr_z/xocc_z(nocc_z),xoccav_z,nd_z,no_z
c
c used for spin-orbit ZORA calculation.
c
c      xocc: contains occupation numbers for kramers degenerate orbitals
c      xoccav_z : fractional occupation open shell (not implemented yet)
c      nd_z : number of orbitals in closed shells
c      no_z : number of orbitals in open shell
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
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
       lshift_z = 0
10161  call inpa4(ytext)
       ystring = ' '
       if (ytext.eq.'atom') then
          nat_z = 2
          icoul_z = 2
          oscalz = .false.
          oscalatz = .false.
       else if (ytext.eq.'mole') then
          nat_z = 0
          icoul_z = 1
          oscalz = .true.
          oscalatz = .false.
       else if (ytext.eq.'rest') then
          irest_z = 2
          call inpa4(ytext)
          numdu_z = locatc(yed,maxlfn,ytext)
          if (numdu_z.ne.0) then
             call inpi(ibldu_z)
          else
             jrec = jrec - 1
          end if        
       else if (ytext.eq.'limi') then
          call inpa4(ystring)
          if (ystring.eq.'s') onlys = .true.
          if (ystring.eq.'p') onlyp = .true.
          if (ystring.eq.'d') onlyd = .true.
          if (ystring.eq.'f') onlyf = .true.
       else if (ytext.eq.'c   ') then
          call inpf(cspeed)
       else if (ytext.eq.'prin') then
          opzora = .true.
       else if (ytext.eq.'pri2') then
          op2zora = .true.
       else if (ytext.eq.'scal') then
          call inpa4(ytext)
          if (ytext.eq.'on') then
               oscalz = .true.
               if (nat_z.eq.2) then
                   oscalz = .false.
                   oscalatz = .true.
               end if
          end if 
          if (ytext.eq.'off') then
             oscalz = .false.
             oscalatz = .false.
          end if
          if (ytext.eq.'iora') then
             oscalz = .false.
             oioraz = .true.
          end if
          if (ytext.eq.'atom') then
             oscalz = .false.
             oscalatz = .true.
          endif
          if (ytext.eq.'mole') then
             if (nat_z.ne.0) 
     1   call caserr('molecular scaling for atomic zora too difficult')
             oscalz = .true.
             oscalatz = .false.
          endif
       else if (ytext.eq.'cont') then
          ocontr = .true.  
       else if (ytext.eq.'exch') then
          call inpa4(ytext)
          if (ytext.eq.'on') isexch_z = 1
          if (ytext.eq.'off') isexch_z = 0
       else if (ytext.eq.'coul') then
          icoul_z = 5
          call inpa4(ytext)
          if (ytext.eq.'full') icoul_z = 1
          if (ytext.eq.'atom') icoul_z = 2
          if (ytext.eq.'none') icoul_z = 3
          if (ytext.eq.'smal') icoul_z = 4
          if (icoul_z.gt.4) call caserr2('unimplemented coul')
       else if (ytext.eq.'nite') then
          call inpi(niter_z)
       else if (ytext.eq.'noco') then
          onoext_z = .true.
       else if (ytext .eq. 'spin'.or.ytext.eq.'so') then
            print *,'spin orbit zora calculation'
            oso = .true.
            maxcyc = 1000      !maximum iterations increased -XJ
            nosymm = 1
            zscftp = 'uhf'
c..      nd and no are the # of fully and fractional occupieds
c..      not operational
c           call inpi(nd_z)
c           call inpi(no_z)
            nd_z = ne
            no_z = 0
       else if (ytext.eq.'inte') then
          call inpa4(ystring)
          if (ystring.eq.'auto') then
             int_zora = 1
          else if (ystring.eq.'none'.or.ystring.eq.'off') then
              int_zora = 0
          else if (ystring.eq.'manu') then
               int_zora = 2
          else
              call caserr2('unknown internal zora option')
          end if
       else if (ytext.eq.'get') then
          call inpa4(ystring)
          if (ystring.eq.'vect') then
             call inpi(is_z)
             if (is_z.le.0.or.is_z.gt.204) call caserr2('get_z section')
          else if (ystring.eq.'dens') then
             call inpi(is_z)
             if (is_z.le.0) call caserr2('getd_z section')
             is_z = -is_z
          else
             call caserr2('illegal zora get option')
          end if
       else if (ytext.eq.'atom') then
             nat_z = 2
             icoul_z = 2
       else if (ytext.eq.'atoi') then
             call inpi(nat_z)
             nat_z = max(nat_z,1)
             icoul_z = 2
       else if (ytext.eq.'crit') then
             call inpi(iiii)
             critat_z = 10.0d0**(-iiii)
       else if (ytext.eq.'none') then
             is_z = 0
             nat_z = 0
             icoul_z = 1
c...         above may change depending on defaults
       else if (ytext.eq.'gaug') then
          call inpi(igauge_z)
       else if (ytext.eq.'off') then
          ozora = .false.
       else if (ytext.eq.'shif') then
          call inpa4(ytext)
          if (ytext.eq.'game') lshift_z = 0
          if (ytext.eq.'henk') lshift_z = 1
          if (ytext.eq.'off') lshift_z = -1
       else if (ytext.eq.' ') then
          call input
          call inpa(ztext)
          if (ztext(1:4).eq.'zora') then
             go to 10161
          else
             jrec = jrec - 1
          end if
          zzstring = '       scaled '
          if (.not.oscalz) then
             if (oioraz) then
                 zzstring = ' iora metric  '
             else if (oscalatz) then
                 zzstring = 'atomic scaled '
                 if (nat_z.le.0.or..not.(icoul_z.eq.2.or.icoul_z.eq.3))
     1           call caserr2('scale atom only for get atom')
             else 
                 zzstring = '              '
             end if 
          end if  
          if (.not.ozora) zzstring ='  *no* '
          zzstr = 'full '
          if (icoul_z.eq.2) zzstr = 'atomic '
          if (icoul_z.eq.3) zzstr = '  no '
          if (icoul_z.eq.4) zzstr = '  small '
          write(iwr,'(a,/,a,a,a,9x,a,/,a,f14.7,27x,a)') 
     1 ' *************************************************************',
     2 ' ** ',zzstring,' ZORA Relativistic Approximation ','**',
     3 ' ** speed of light ',cspeed,'**'
          if (oioraz) write(iwr,'(a,/a,/a)')
     1 ' ** experimental IORA using just large-component density    **',
     2 ' ** after K.G.Dyall,E.van Lenthe, JCP 111(4),1366-1372(1999)**',
     3 ' ** but no eq. (53); renormalised density used throughout   **'
          if (ystring.ne.' ') write(iwr,'(a,a,a)')
     1 ' ** Internal projection only up to ',ystring,'-functions **'
          if (nat_z.eq.0) then
             write(iwr,'(a,a,a,28x,a)') ' ** ',zzstr,
     1                                  ' coulomb matrix used','**'
          else
             write(iwr,'(a,a,a,27x,a)') ' ** ','atomic ',
     1                                  ' coulomb and zora used','**'
          end if
          if (isexch_z.eq.1)
     1    write(iwr,'(a,a,a,28x,a)') ' ** ','        ',
     1                               ' exchange included  ','**'
          if (is_z.eq.0.and.nat_z.eq.0) then
             write(iwr,'(a,i5,a)') ' ** coulomb matrix updated every ',
     1             niter_z,' iteration            **'
          else if (is_z.gt.0) then
             write(iwr,'(2a,i5,a)')' ** coulomb matrix from vectors at',
     1             ' section ',is_z,'            **'
          else if (is_z.lt.0) then
             write(iwr,'(2a,i5,a)')' ** coulomb matrix from density at',
     1             ' section ',-is_z,'           **'
          else if (nat_z.gt.0) then
             if (irest_z.eq.0) then
                write(iwr,'(2a,i4,a,/,a,20x,a,1pe7.1,a)') 
     1                ' ** coulomb matrix from start',
     1                'density (',nat_z,' atom iterations) **',' **',
     1                'converge atomic energy till ',critat_z,'  **'
             else
              if (numdu_z.eq.0) then
                write(iwr,'(2a     )') ' ** ZORA corrections from cur',
     1                'rent dumpfile                  **'
              else
                write(iwr,'(3a,i5,a)') ' ** ZORA corrections from ',
     1          yed(numdu_z),' at block',ibldu_z,'                **'
              end if
             end if
          else if (irest_z.eq.0) then
             call caserr('Restore only for atomic ZORA')
          end if
          if (igauge_z.ne.0) write(iwr,'(a,i10,23x,a)')
     1       ' ** GAUGE shifted by      ',igauge_z,' **'
          write(iwr,'(a)') 
     1 ' *************************************************************'
          if (int_zora.eq.2) call zora_int
          go to 30
       else
          call caserr('unrecognised zora option')
       end if
       go to 10161
c
  30   continue
c
       return
       end
      subroutine dft_defaults(atomss,grd_type,natomss,ierr0,ierr1)
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer atomss(*), natomss,ierr0,ierr1,ierror
      integer grd_type(*)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert

c     Bfcheck is used to test whether an atom/bq-centre has a basis set.
      logical bfcheck(maxat)
c
c     establish a default set of attributes at the first occurrence of
c     a cdft directives. At present these are set as:
c     B-LYP, medium quadrature, explicit J
c     N.B. consistency check will be required here
c     when overriding these
c
c     grd_type(i) == -1 Means there is NO grid on this centre!!!

c
c     drive the dft module through AO+MO code by default
c     (although AO+MO should be best we do not have the timing
c      right yet, so use AO by default).
c
      ks_bas     = KS_AO
      ks_grd_bas = KS_GRD_AO
      ks_rhs_bas = KS_RHS_MO
      ks_lhs_bas = KS_LHS_MO
      ks_dt_bas  = KS_DT_MO
      ks_hes_bas = KS_HES_AO
      ks_dx_bas  = KS_DX_AO
c
c     symmetrise the XCorr matrix
c
      osym_op = .true.
      ierr1   = 0
c
c       explicit coulomb
c
      ierr0=CD_4c2eon()
      if(ierr0.ne.0) ierr1=ierr1+1
c
c       B-LYP functional
c
      ierr0=CD_set_functional('blyp')
      if(ierr0.ne.0) ierr1=ierr1+1
      do i = 1, nat
         atomss(i)=isubst(zaname(i))
         grd_type(i)  = -1
      enddo

c the geometry is now imported from start
c hvd: but do it again here for now as at this point we know the 
c      correct number of electrons in case ECPs are used...
c      a neater solution would be better...
      ierr0=CD_import_geom(nat,ne,c,atomss)

      if(ierr0.ne.0) ierr1=ierr1+1

      do i=1, nat
         bfcheck(i) = .false.
      enddo
      do i=1, nshell
         bfcheck(katom(i)) = .true.
      enddo

      do 6666 i = 1, nat
         if( .not. bfcheck(i))goto 6666
         if (zaname(i)(1:2).eq.'bq') then
            write(iwr,*)'Default (Carbon) grid assigned to BQ centre',i
            write(iwr,*)'You may wish to provide grid parameters'
            do j = 1, i-1
               if(atomss(j).eq.atomss(i).and.bfcheck(j)) then
                  grd_type(i)=grd_type(j)
                  goto 6666
               endif
            enddo
            grd_type(i)=CD_create_grid()
            goto 6666
         endif
         do j = 1, i-1
            if(atomss(j).eq.atomss(i)) then
               grd_type(i)=grd_type(j)
               goto 6666
            endif
         enddo
         grd_type(i)=CD_create_grid()
 6666 continue
c
c       Created initial grid types.
c       Variable atomss no longer needed and reused for input processing
c
      do i=1,nat
         if (grd_type(i).eq.-1) then
            ierr0=CD_assign_grid(i,0)
         else
            ierr0=CD_assign_grid(i,grd_type(i))
         endif
         if(ierr0.ne.0) ierr1=ierr1+1
      enddo
c
      if (ierr1.ne.0) then
        call caserr2('problems with default dft options')
      endif
c
      return
      end
      subroutine dft_input(atomss,grd_type,natomss,odftgrad,ierr0,ierr1,
     +                     ierror)
      implicit none
c
c     process one or more cdft/dftdirectives
c
c     implicit real*8 (a-h,p-w),integer (i-n),logical (o)
c     implicit character *8 (z),character *1 (x)
c     implicit character *4 (y)
c
      integer atomss(*), natomss,ierr0,ierr1,ierror
      integer grd_type(*)
      logical oscreen
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
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
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c... zorac is included for nat_z only
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
c
c     COMMON/DRIVE_DFT/ contains options that specify how GAMESS-UK
c     should drive the CCP1 DFT module.
c
      integer KS_AO,     KS_MO,     KS_AOMO
      integer KS_GRD_AO, KS_GRD_MO, KS_GRD_AOMO
      integer KS_RHS_AO, KS_RHS_MO, KS_RHS_AOMO
      integer KS_LHS_AO, KS_LHS_MO, KS_LHS_AOMO
      integer KS_DT_AO,  KS_DT_MO,  KS_DT_AOMO
      integer KS_HES_AO, KS_HES_MO, KS_HES_AOMO
      integer KS_DX_AO,  KS_DX_MO,  KS_DX_AOMO
      parameter (KS_AO    =01,KS_MO    =02,KS_AOMO    =03)
      parameter (KS_GRD_AO=04,KS_GRD_MO=05,KS_GRD_AOMO=06)
      parameter (KS_RHS_AO=07,KS_RHS_MO=08,KS_RHS_AOMO=09)
      parameter (KS_LHS_AO=10,KS_LHS_MO=11,KS_LHS_AOMO=12)
      parameter (KS_DT_AO =13,KS_DT_MO =14,KS_DT_AOMO =15)
      parameter (KS_HES_AO=16,KS_HES_MO=17,KS_HES_AOMO=18)
      parameter (KS_DX_AO =19,KS_DX_MO =20,KS_DX_AOMO =21)
c
c     if ks_bas.eq.KS_AO           call CD_energy_ao
c     if ks_bas.eq.KS_MO           call CD_energy_mo
c     if ks_bas.eq.KS_AOMO         call CD_energy
c
c     if ks_grd_bas.eq.KS_GRD_AO   call CD_force_ao
c     if ks_grd_bas.eq.KS_GRD_MO   call CD_force_mo
c     if ks_grd_bas.eq.KS_GRD_AOMO call CD_force
c
c     if ks_rhs_bas.eq.KS_RHS_AO   call CD_chf_rhs_ao
c     if ks_rhs_bas.eq.KS_RHS_MO   call CD_chf_rhs_mo
c     if ks_rhs_bas.eq.KS_RHS_AOMO call CD_chf_rhs
c
c     if ks_lhs_bas.eq.KS_LHS_AO   call CD_chf_lhs_ao
c     if ks_lhs_bas.eq.KS_LHS_MO   call CD_chf_lhs_mo
c     if ks_lhs_bas.eq.KS_LHS_AOMO call CD_chf_lhs
c
c     if ks_dt_bas.eq.KS_DT_AO     call CD_dksm_ao
c     if ks_dt_bas.eq.KS_DT_MO     call CD_dksm_mo
c     if ks_dt_bas.eq.KS_DT_AOMO   call CD_dksm
c
c     if ks_hes_bas.eq.KS_HES_AO   call CD_hess_ao
c     if ks_hes_bas.eq.KS_HES_MO   call CD_hess_mo
c     if ks_hes_bas.eq.KS_HES_AOMO call CD_hess
c
c     if ks_dx_bas.eq.KS_DX_AO     call CD_dksm_exp_ao
c     if ks_dx_bas.eq.KS_DX_MO     call CD_dksm_exp_mo
c     if ks_dx_bas.eq.KS_DX_AOMO   call CD_dksm_exp
c
      integer ks_bas,     ks_grd_bas, ks_hes_bas
      integer ks_rhs_bas, ks_lhs_bas, ks_dt_bas
      integer ks_dx_bas
      logical ogeompert
c
      common/drive_dft/ks_bas, ks_grd_bas, 
     &       ks_rhs_bas, ks_lhs_bas, ks_dt_bas,
     &       ks_hes_bas, ks_dx_bas,
     &       ogeompert
cDEBUG
      integer batch_counter
      logical screen_save_sw, big_mem_sw
      common/gezeik/batch_counter,screen_save_sw,big_mem_sw
cDEBUG
c
      character*8 ztext, zblank
      character*4 ytext, ytrunc, yblank
      logical odftgrad, onumb, oisint, omore
      integer ncent, icut_dft, itol_dft, icutd_dft, itold_dft
      integer nangzn, nphizn, i, nzone, imemo, itmp, nrad, ii, irow
      integer igen
      real*8 rbndzn, p(5), radm, fac
      real*8 dft_exptol, dft_psitol, dft_rhotol, dft_dentol, dft_wghtol,
     &     dft_wrntol
      dimension nangzn(20), nphizn(20), rbndzn(19)
      data yblank/' '/
      data zblank/' '/
c
10131 call inpa(ztext)
      ytext=ytrunc(ztext)
c
      ierror=0
      if(ytext.eq.'    ') then
        goto 20
      else if(ytext.eq.'4c2e') then
         ierror=CD_4c2eon()
      else if(ytext.eq.'atom') then
c
c        dft atom [on|off|<#it>]
c
         call inpa(ztext)
         if (ztext.ne.zblank) then
            if (ztext.eq.'on') then
               oatmdft=.true.
            else if (ztext.eq.'off') then
               oatmdft=.false.
            else
               oatmdft=.true.
               call intval(ztext,itmp,oisint)
               if (oisint) then
                  nat_z = itmp
               else
                  jrec = jrec - 1
               endif
            endif
         else
            oatmdft=.true.
         endif
         odenscfdft=odenscfdft.or.oatmdft
      else if(ytext.eq.'dens') then
c
c        dft denscf [on|off]
c
         call inpa(ztext)
         if (ztext.ne.zblank) then
            if (ztext.eq.'on') then
               odenscfdft=.true.
            else if (ztext.eq.'off') then
               odenscfdft=.false.
            else
               odenscfdft=.true.
               jrec = jrec - 1
            endif
         else
            odenscfdft=.true.
         endif
      else if(ztext.eq.'jfit') then
         ierror=CD_jfiton(.false.)
         if (ks_bas.eq.KS_AOMO) ks_bas = KS_AO
         call inpa(ztext)
         if (ztext.eq.'memory') then
            ierror=CD_jfiton(.true.)
            if (ierror.ne.0) then 
               call caserr2('Could not obtain memory for Coulomb fit!')
            endif
            call inpa(ztext)
            if (ztext.ne.yblank) then
               call intval(ztext,imemo,onumb)
               if (onumb) then
                  call store_spare(imemo)
               else
                  jrec = jrec - 1
               endif
            endif
         else
            jrec = jrec - 1
         endif
      else if(ztext.eq.'jfitg') then
         ierror=CD_jfitgon()
         if (ks_grd_bas.eq.KS_GRD_AOMO) ks_grd_bas = KS_GRD_AO
      else if(ytext.eq.'jmul') then
         ierror=CD_jmulton()
      else if(ytext.eq.'xcfi') then
         ierror=CD_xcfiton()
c
c...  attempt to set functional (if this fails then
c...  try the various aliases)
c
      else if(0.eq.CD_set_functional(ztext)) then
         call vdwaals_functional(ztext)
c
c...  aliases for exchange functionals 
c
      else if(ytext.eq.'hfex') then
         ierror=CD_set_functional('hf_x')
      else if(ytext.eq.'beck') then
         ierror=CD_set_functional('b88_x')
      else if(ytext.eq.'lda') then
         ierror=CD_set_functional('lda_x')
      else if(ytext.eq.'xa-f') then
         ierror=CD_set_functional('ft97a_x')
      else if(ytext.eq.'x-ft'.or.ytext.eq.'xb-f') then
         ierror=CD_set_functional('ft97b_x')
      else if(ytext.eq.'b3') then
         ierror=CD_set_functional('b3_x')
      else if(ytext.eq.'xpbe') then
         ierror=CD_set_functional('pbe_x')
c
c...  aliases for correlation functionals
c
      else if(ytext.eq.'vwnr') then
         ierror=CD_set_functional('vwn5rpa')
      else if(ytext.eq.'c-ft') then
         ierror=CD_set_functional('ft97_c')
      else if(ytext.eq.'noco') then
         ierror=CD_set_functional('null_c')
c
c...  aliases for exchange-correlation functionals
c
      else if(ytext.eq.'b-ly'.or.ytext.eq.'blyp') then
         ierror=CD_set_functional('blyp')
         call vdwaals_functional('blyp')
      else if(ytext.eq.'s-vw'.or.ytext.eq.'svwn') then
         ierror=CD_set_functional('svwn')
      else if(ytext.eq.'bp86'.or.ytext.eq.'b-p8') then
         ierror=CD_set_functional('bp86')
         call vdwaals_functional('bp86')
c
c...  end of functionals
c
      else if(ytext.eq.'angp') then
c
c        valid inputs:
c
c           angprune 
c
c           angprune on/off/auto
c
c           angprune label c1 c2 on/off/auto
c
c           angprune label c1 o1 on/off/auto label c2 h2 on/off/auto
c
         ierror=0
         call inpa4(ytext)
         if (ytext.eq.' ') then
            ierror=CD_MHL_ang_prune(0,'on')
         else if (ytext.eq.'auto') then 
            ierror=CD_auto_ang_prune(0,'on')
         else if (ytext.eq.'labe'.or.ytext.eq.'elem') then
 101        call read_types(nat,zaname,atomss,natomss,ytext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
            call inpa4(ytext)
            if (ytext.eq.'auto') then
               do i=1,natomss
                  ierror=ierror
     +                  +CD_auto_ang_prune(grd_type(atomss(i)),'on')
               enddo
            else
               do i=1,natomss
                  ierror=ierror
     +                  +CD_MHL_ang_prune(grd_type(atomss(i)),ytext)
               enddo
            endif
            call inpa4(ytext)
            if (ytext.eq.'labe'.or.ytext.eq.'elem') then
               goto 101
            else if (ytext.ne.' ') then
               jrec=jrec-1
            endif
         else
            ierror=CD_MHL_ang_prune(0,ytext)
         endif
         if (ierror.ne.0) then
            write(iwr,*)'Invalid suboption with angprune:',ytext
            call caserr2('Invalid suboption with angprune')
         endif
      else if(ytext.eq.'grad') then
         call inpa4(ytext)
         odftgrad = .true.
         if(ytext.eq.'on'.or.ytext.eq.'yes') then
            ierror=CD_gradquad(.true.)
         else if(ytext.eq.'off'.or.ytext.eq.'no') then
            ierror=CD_gradquad(.false.)
         else
            write(iwr,*)'Invalid suboption with gradients of ',
     +                  'quadrature:',ytext
            call caserr2('Invalid option with gradients of quadrature')
         endif
      else if(ytext.eq.'sort') then
c
c        valid inputs:
c
c           sort on 
c
c           sort off
c
         call inpa4(ytext)
         if(ytext.eq.'on'.or.ytext.eq.'yes') then
            ierror=CD_sortpoints(.true.)
         else if(ytext.eq.'off'.or.ytext.eq.'no') then
            ierror=CD_sortpoints(.false.)
         else
            write(iwr,*)'Invalid suboption with sorting of grid ',
     +                  'points:',ytext
            call caserr2('Invalid option with sort grid points')
         endif
      else if(ytext.eq.'lebe') then
c
c        valid inputs:
c
c           lebedev 302
c
c           lebedev 194 0.1 302 0.5 434
c
c           lebedev label c1 c2 194 0.1 302 0.5 434 element cl 590
c
c           lebedev row 1 194
c
c           lebedev rows 194  266  302  266  302  302  266
c
 3002    natomss=-1
         call inpa4(ytext)
         if (ytext.eq.'labe'.or.ytext.eq.'elem') then
            call read_types(nat,zaname,atomss,natomss,ytext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else if (ytext(1:4).eq.'rows') then
            do irow = 1, 7
               call inpi(nangzn(1))
               ierror=CD_ang_npoints_row(irow,nangzn(1))
            enddo
            go to 10131
         else if (ytext(1:3).eq.'row') then
            call inpi(irow)
         else
            jrec=jrec-1
         endif
         call inpi(nangzn(1))
         nzone = 1
 3000    call inpa(ztext)
         if ((ztext(1:1).lt.'a'.or.'z'.lt.ztext(1:1))
     +       .and.ztext(1:1).ne.' ') then
            jrec = jrec - 1
            call inpf(rbndzn(nzone))
            nzone = nzone + 1
            call inpi(nangzn(nzone))
            goto 3000
         endif
         jrec = jrec - 1
         if (ytext.eq.'row') then
            ierror=CD_ang_npoints_row(irow,nangzn(1))
            if (ierror.eq.1) then
               call caserr('Row of periodic table out of range')
            else if (ierror.eq.2) then
               call caserr(
     &              'Number of angular grid points must exceed 0')
            endif
         else if (natomss.eq.-1) then
            ierror=CD_lebedevon(0,nzone,nangzn,rbndzn)
            if (ierror.eq.1) then
               call caserr2('Grid type out of range')
            else if (ierror.eq.2) then
               call caserr2('invalid number of radial zones')
            else if (ierror.eq.3.or.ierror.eq.4) then
               call caserr2(
     +              'zone boundaries should be ordered increasingly')
            else if (ierror.eq.5) then
               call caserr2('No lebedev grid for the specified size')
            endif
         else
            do i = 1, natomss
               ierror=CD_lebedevon(grd_type(atomss(i)),
     +                             nzone,nangzn,rbndzn)
               if (ierror.eq.1) then
                  call caserr2('Grid type out of range')
               else if (ierror.eq.2) then
                  call caserr2('invalid number of radial zones')
               else if (ierror.eq.3.or.ierror.eq.4) then
                  call caserr2(
     +               'zone boundaries should be ordered increasingly')
               else if (ierror.eq.5) then
                  call caserr2(
     +               'No lebedev grid for the specified size')
               endif
            enddo
         endif
         if (ztext(1:4).eq.'labe'.or.ztext(1:4).eq.'elem'.or.
     +       ztext(1:3).eq.'row') goto 3002
      else if(ytext.eq.'gaus') then
c
c        valid inputs:
c
c           gausslegendre 302
c
c           gausslegendre 194 0.1 302 0.5 434
c
c           gausslegendre label c1 c2 194 0.1 302 0.5 434 element cl 590
c
c           gausslegendre row 1 194
c
 3012    natomss=-1
         call inpa4(ytext)
         if (ytext.eq.'labe'.or.ytext.eq.'elem') then
            call read_types(nat,zaname,atomss,natomss,ytext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else if (ytext.eq.'row') then
            call inpi(irow)
         else
            jrec=jrec-1
         endif
         nzone = 1
         call inpi(nangzn(1))
 3010    call inpa(ztext)
         if ((ztext(1:1).lt.'a'.or.'z'.lt.ztext(1:1))
     +       .and.ztext(1:1).ne.' ') then
            jrec = jrec - 1
            call inpf(rbndzn(nzone))
            nzone = nzone + 1
            call inpi(nangzn(nzone))
            goto 3010
         endif
         jrec = jrec - 1
         do i = 1, nzone
            nphizn(i)=nangzn(i)*2
         enddo
         if (ytext.eq.'row') then
            ierror=CD_ang_npoints_row(irow,nangzn(1))
         else if (natomss.eq.-1) then
            ierror=CD_gausslon(0,nzone,nangzn,nphizn,rbndzn)
            if (ierror.eq.1) then
               call caserr2('Grid type out of range')
            else if (ierror.eq.2) then
               call caserr2('invalid number of radial zones')
            else if (ierror.eq.3.or.ierror.eq.4) then
               call caserr2(
     +              'zone boundaries should be ordered increasingly')
            endif
         else 
            do i = 1, natomss
               ierror=CD_gausslon(grd_type(atomss(i)),
     +                            nzone,nangzn,nphizn,rbndzn)
               if (ierror.eq.1) then
                  call caserr2('Grid type out of range')
               else if (ierror.eq.2) then
                  call caserr2('invalid number of radial zones')
               else if (ierror.eq.3.or.ierror.eq.4) then
                  call caserr2(
     +                 'zone boundaries should be ordered increasingly')
               endif
            enddo
         endif
         if (ztext(1:4).eq.'labe'.or.ztext(1:4).eq.'elem'.or.
     +       ztext(1:3).eq.'row') goto 3012
      else if(ytext.eq.'eule' .or. ytext .eq.'log') then
c
c        valid input:
c
c           log 45 3.0
c
c           log label c1 h1 45 3.0
c
c           log element c h 45 3.0 label o1 20 3.0
c
c           log row 1 45 
c
c           log rows 20   35   35   35   35   40   35
c
c           euler 45
c
c           euler label c1 h1 45
c
c           euler element c h 45 label o1 20 
c
c           euler row 1 45 
c
         call inpa(ztext)
 3015    continue
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else if (ztext.eq.'rows') then
            do irow = 1, 7
               call inpi(nrad)
               ierror=ierror+CD_rad_npoints_row(irow,nrad)
            enddo
            go to 10131
         else if (ztext.eq.'row') then
            call inpi(irow)
         else
            jrec=jrec-1
         endif
         call inpi(nrad)
         if (ytext.eq.'log') then
            if (ztext.ne.'row') then
               call inpf(radm)
            endif
            ierror=0
            if (ztext.eq.'label'.or.ztext.eq.'element') then
               do i=1,natomss
                  ierror=CD_logon(grd_type(atomss(i)),nrad,radm)
                  if (ierror.eq.1) call caserr2('grid out of range')
                  if (ierror.eq.2) call caserr2('invalid num points')
                  if (ierror.eq.3) call caserr2('exponent missing')
                  if (ierror.ne.0) call caserr2('CD_logon failed')
               enddo
            else if (ztext.eq.'row') then
               ierror=ierror+CD_rad_npoints_row(irow,nrad)
            else
               ierror=CD_logon(0,nrad,radm)
               if (ierror.eq.1) call caserr2('grid out of range')
               if (ierror.eq.2) call caserr2('invalid num points')
               if (ierror.eq.3) call caserr2('exponent missing')
               if (ierror.ne.0) call caserr2('CD_logon failed')
            endif
            if (ierror.ne.0) then
               call caserr2('CD_logon failed')
            endif
         else
            ierror=0
            if (ztext.eq.'label'.or.ztext.eq.'element') then
               do i=1,natomss
                  ierror=ierror+CD_euleron(grd_type(atomss(i)),nrad)
               enddo
            else if (ztext.eq.'row') then
               ierror=ierror+CD_rad_npoints_row(irow,nrad)
            else
               ierror=ierror+CD_euleron(0,nrad)
            endif
            if (ierror.ne.0) then
               call caserr2('CD_euleron failed')
            endif
         endif
         call inpa(ztext) 
         call intval(ztext,nrad,onumb)
         if (ztext.eq.'label'.or.ztext.eq.'element'.or.
     &       ztext.eq.'row'.or.onumb) then
            goto 3015
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'scal') then
c
c        valid inputs:
c
c           scale 3.0
c
c           scale label c1 h1 3.0 label c2 o1 4.0
c
c           scale 4.0 label c1 h1 3.0
c
c           scale 4.0 element c1 h1 3.0
c
         call inpa(ztext)
 3022    natomss=-1
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else 
            jrec = jrec-1
         endif
         call inpf(fac)
         ierror=0
         if (natomss.eq.-1) then
            ierror=ierror+CD_gridscale(0,fac)
         else
            do i=1,natomss
               ierror=ierror+CD_gridscale(grd_type(atomss(i)),fac)
            enddo
         endif
         if (ierror.ne.0) call caserr2('CD_gridscale failed!!!')
         call inpa(ztext)
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            goto 3022
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'radi') then
c
c        valid inputs:
c
c           radii 3.0
c
c           radii label c1 h1 3.0 label c2 o1 4.0
c
c           radii 4.0 label c1 h1 3.0
c
c           radii 4.0 element c1 h1 3.0
c
c        purpose: change the "atom sizes" used for the radial grid 
c                 scale factors
c
         call inpa(ztext)
 3024    natomss=-1
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else 
            jrec = jrec-1
         endif
         call inpf(fac)
         ierror=0
         if (natomss.eq.-1) then
            ierror=ierror+CD_gridatomradius(0,fac)
         else
            do i=1,natomss
               ierror=ierror+CD_gridatomradius(grd_type(atomss(i)),fac)
            enddo
         endif
         if (ierror.ne.0) call caserr('CD_gridatomradius failed!!!')
         call inpa(ztext)
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            goto 3024
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'gsca') then
c
c        valid inputs:
c
c           gscale mk   - recommended by Mura-Knowles
c
c           gscale gam1 - 3.3 * Bragg-Slater radius
c
c           gscale gam2 - optimized parameters
c
c        purpose: change radial grid scale scheme
c
         call inpa4(ytext)
         ierror = CD_radscale_scheme(ytext)
c
      else if(ytext.eq.'wrad') then
c
c        valid inputs:
c
c           wradii 3.0
c
c           wradii label c1 h1 3.0 label c2 o1 4.0
c
c           wradii 4.0 label c1 h1 3.0
c
c           wradii 4.0 element c1 h1 3.0
c
c        purpose: change the atomic size adjustment radii for the 
c                 weighting schemes
c
         call inpa(ztext)
 3025    natomss=-1
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else 
            jrec = jrec-1
         endif
         call inpf(fac)
         ierror=0
         if (natomss.eq.-1) then
            ierror=ierror+CD_weightatomradius(0,fac)
         else
            do i=1,natomss
               ierror=ierror+CD_weightatomradius(grd_type(atomss(i)),
     &                                           fac)
            enddo
         endif
         if (ierror.ne.0) call caserr('CD_weightatomradius failed!!!')
         call inpa(ztext)
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            goto 3025
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'prad') then
c
c        valid inputs:
c
c           pradii 3.0
c
c           pradii label c1 h1 3.0 label c2 o1 4.0
c
c           pradii 4.0 label c1 h1 3.0
c
c           pradii 4.0 element c1 h1 3.0
c
c        purpose: change the atomic radii for the 
c                 angular grid pruning schemes
c
         call inpa(ztext)
 3026    natomss=-1
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else 
            jrec = jrec-1
         endif
         call inpf(fac)
         ierror=0
         if (natomss.eq.-1) then
            ierror=ierror+CD_pruneatomradius(0,fac)
         else
            do i=1,natomss
               ierror=ierror+CD_pruneatomradius(grd_type(atomss(i)),
     &                                           fac)
            enddo
         endif
         if (ierror.ne.0) call caserr('CD_pruneatomradius failed!!!')
         call inpa(ztext)
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            goto 3026
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'srad') then
c
c        valid inputs:
c
c           sradii 3.0
c
c           sradii label c1 h1 3.0 label c2 o1 4.0
c
c           sradii 4.0 label c1 h1 3.0
c
c           sradii 4.0 element c1 h1 3.0
c
         call inpa(ztext)
 3027    natomss=-1
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            call read_types(nat,zaname,atomss,natomss,ztext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,iwr)
         else 
            jrec = jrec-1
         endif
         call inpf(fac)
         ierror=0
         if (natomss.eq.-1) then
            ierror=ierror+CD_screenatomradius(0,fac)
         else
            do i=1,natomss
               ierror=ierror+CD_screenatomradius(grd_type(atomss(i)),
     &                                           fac)
            enddo
         endif
         if (ierror.ne.0) call caserr2('CD_screenatomradius failed!!!')
         call inpa(ztext)
         if (ztext.eq.'label'.or.ztext.eq.'element') then
            goto 3027
         else if (ztext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'grid') then
c
c        valid inputs
c
c           grid label c1 h2
c
c           grid element c h
c
         call read_types(nat,zaname,atomss,natomss,ytext)
         call new_gridlist(nat,grd_type,natomss,atomss,zaname,
     +                     iwr)
c
      else if(ytext.eq.'quad') then
c
c        valid inputs
c
c           quadrature old
c
c           quadrature generation 1
c
c           quadrature medium
c
c           quadrature generation 2 medium
c
c           quadrature label c1 h2 high
c
c           quadrature low label c2 c3 veryhigh
c
c           quadrature medium ignore
c
 4011    call inpa4(ytext)
         if(ytext.eq.'old') then
c
c           use to reset to original dft defaults for integration grids
c
            ierror = CD_defaults_old()
         else if (ytext.eq.'gene') then
c
c           get the generation number for the grids
c
            call inpi(igen)
            ierror = CD_generation(igen)
            if (ierror.ne.0) call caserr2('CD_generation failed!')
            goto 4011
         else
 4012       if (ytext.eq.'labe'.or.ytext.eq.'elem') then
               call read_types(nat,zaname,atomss,natomss,ytext)
               call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,
     +                              iwr)
               call inpa4(ytext)
               do i = 1, natomss
                  ierror=CD_accuracy(grd_type(atomss(i)),ytext)
                  write(iwr,4013)ytext
                  if (ierror.ne.0) call caserr2('CD_accuracy failed!')
               enddo
            else
               ierror=CD_accuracy(0,ytext)
               if (ierror.ne.0) then
                  write(iwr,4013)ytext
                  call caserr2('CD_accuracy failed!')
               endif
            endif
            call inpa4(ytext)
            if (ytext.eq.'labe'.or.ytext.eq.'elem') then
               goto 4012
            else if (ytext.eq.'igno') then
               ierror = CD_set_ignore_accuracy(.true.)
            else if (ytext.ne.' ') then
               jrec=jrec-1
            endif
 4013       format(" Unrecognised predefined grid: ",a4,/,
     +             " Expected one of 'low', 'medium', 'high', ",
     +             "'veryhigh', 'ref' or 'sg1'")
         endif
      else if(ytext.eq.'warn') then
c
c        valid inputs
c
c           warn 1.0d-6
c
c           warn label c1 h2 1.0d-8
c
c           warn 1.0d-4 label c2 c3 1.0d-8
c
         dft_wrntol = -1.0d0
         call inpa4(ytext)
 4014    if (ytext.eq.'labe'.or.ytext.eq.'elem') then
            call read_types(nat,zaname,atomss,natomss,ytext)
            call reduce_gridlist(nat,grd_type,natomss,atomss,zaname,
     +                           iwr)
            call inpf(dft_wrntol)
            do i = 1, natomss
               ierror=CD_warn(grd_type(atomss(i)),dft_wrntol)
               if (ierror.ne.0) call caserr2('CD_warn failed!')
            enddo
         else
            jrec=jrec-1
            call inpf(dft_wrntol)
            ierror=CD_warn(0,dft_wrntol)
            if (ierror.ne.0) call caserr2('CD_warn failed!')
         endif
         call inpa4(ytext)
         if (ytext.eq.'labe'.or.ytext.eq.'elem') then
            goto 4014
         else if (ytext.ne.' ') then
            jrec=jrec-1
         endif
      else if(ytext.eq.'abor') then
         ierror=CD_abort()
      else if(ytext.eq.'jbas') then
	call store_jbasis
      else if(ytext.eq.'aoba') then
cDEBUG
        write(*,*)'*** read AOBAS'
cDEBUG
 4020   call inpa(ztext)
        if (ztext.eq.'energy') then
           ks_bas = KS_AO
           goto 4020
        else if (ztext.eq.'force') then
           ks_grd_bas = KS_GRD_AO
           goto 4020
        else if (ztext.eq.'rhs') then
           ks_rhs_bas = KS_RHS_AO
           goto 4020
        else if (ztext.eq.'lhs') then
           ks_lhs_bas = KS_LHS_AO
           goto 4020
        else if (ztext.eq.'dks') then
           ks_dt_bas = KS_DT_AO
           goto 4020
        else if (ztext.eq.'hessian') then
           ks_hes_bas = KS_HES_AO
           goto 4020
        else if (ztext.eq.'dksx') then
           ks_dx_bas = KS_DX_AO
           goto 4020
        else if (ztext.ne.' ') then
           jrec=jrec-1
        endif
      else if(ytext.eq.'moba') then
cDEBUG
        write(*,*)'*** read MOBAS'
cDEBUG
 4022   call inpa(ztext)
        if (ztext.eq.'energy') then
           ks_bas = KS_MO
           goto 4022
        else if (ztext.eq.'force') then
           ks_grd_bas = KS_GRD_MO
           goto 4022
        else if (ztext.eq.'rhs') then
           ks_rhs_bas = KS_RHS_MO
           goto 4022
        else if (ztext.eq.'lhs') then
           ks_lhs_bas = KS_LHS_MO
           goto 4022
        else if (ztext.eq.'dks') then
           ks_dt_bas = KS_DT_MO
           goto 4022
        else if (ztext.eq.'hessian') then
           ks_hes_bas = KS_HES_MO
           goto 4022
        else if (ztext.eq.'dksx') then
           ks_dx_bas = KS_DX_MO
           goto 4022
        else if (ztext.ne.' ') then
           jrec=jrec-1
        endif
      else if(ytext.eq.'aomo') then
cDEBUG
        write(*,*)'*** read AOMOBAS'
cDEBUG
 4024   call inpa(ztext)
        if (ztext.eq.'energy') then
           ks_bas = KS_AOMO
           goto 4024
        else if (ztext.eq.'force') then
           ks_grd_bas = KS_GRD_AOMO
           goto 4024
        else if (ztext.ne.' ') then
           jrec=jrec-1
        endif
      else if(ytext.eq.'weig') then
	call inpa(ztext)
	ierror = CD_set_weight(ztext)
      else if(ytext.eq.'over') then
c
c  tolerance for shell overlap test
c
         call inpi(ii)
         ierror = CD_over(ii)
      else if(ytext.eq.'pene') then
c
c  tolerance for multpole penetration test
c
         call inpi(ii)
         ierror = CD_pener(ii)
      else if(ytext.eq.'pole') then
c
c  tolerance for multpole penetration test
c
         call inpi(ii)
         ierror = CD_pole(ii)
      else if(ytext.eq.'schw') then
c
c  tolerance for Schwarz inequality test (3 centre 2 electron Jfit)
c
         call inpi(ii)
         ierror = CD_schwarz(ii)
      else if(ytext.eq.'scre') then
c
c enable screening + set the tolerances
c
         dft_exptol = -1.0d0
         dft_psitol = -1.0d0
         dft_rhotol = -1.0d0
         dft_dentol = -1.0d0
         dft_wghtol = -1.0d0

         oscreen = .true.
         omore = .true.
         do while (omore) 
            call inpa4(ytext)
            if(ytext.eq.'p')then
               call inpf(dft_dentol)
            else if(ytext.eq.'rho')then
               call inpf(dft_rhotol)
            else if(ytext.eq.'weig')then
               call inpf(dft_wghtol)
            else if(ytext.eq.'psi')then
               call inpa(ztext)
 3032          natomss=-1
               if (ztext.eq.'label'.or.ztext.eq.'element') then
                  call read_types(nat,zaname,atomss,natomss,ztext)
                  call reduce_gridlist(nat,grd_type,natomss,atomss,
     +                                 zaname,iwr)
               else
                  jrec=jrec-1
               endif
               call inpf(dft_psitol)
               ierror=0
               if (natomss.eq.-1) then
                  ierror=CD_psitol(0,dft_psitol)
               else
                  do i = 1, natomss
                     ierror=ierror
     +                     +CD_psitol(grd_type(atomss(i)),dft_psitol)
                  enddo
               endif
               if (ierror.ne.0) call caserr2('CD_psitol failed')
               call inpa(ztext)
               if (ztext.eq.'label'.or.ztext.eq.'element') then
                  goto 3032
               else if (ztext.ne.' ') then
                  jrec=jrec-1
               endif
            else if(ytext.eq.'conv') then
               ierror = CD_conv_prune_on()
            else if(ytext.eq.'off') then
               oscreen = .false.
            else
               omore = .false.
               jrec = jrec - 1
            endif
         enddo

         ierror=CD_screen(oscreen,dft_dentol,dft_rhotol,dft_wghtol)

      else if (ytext.eq.'nosy') then
         osym_op = .false.

      else if(ytext.eq.'intt') then
c
c set tolerances for coulomb integral and derivative codes
c
         call inpa4(ytext)

         if(ytext.eq.'2cen')then
            ncent = 2
         else if(ytext.eq.'3cen')then
            ncent = 3
         elseif(ytext.eq.'4cen')then
            ncent = 4
         else if(ytext.eq.'all ')then
            ncent = 0
         else
            write(iwr,10133)
10133       format(1x,'cdft directive inttol must be followed by',
     &           'one of 2centre, 3centre, 4centre or all')
            call gamerr('bad cdft inttol directive',
     &      ERR_NO_CODE, ERR_INCOMPREHENSIBLE, ERR_SYNC, ERR_NO_SYS)
         endif

         icut_dft  = 0
         itol_dft  = 0
         icutd_dft = 0
         itold_dft = 0

         omore = .true.
         do while (omore) 
            call inpa4(ytext)
            if(ytext.eq.'cut ')then
               call inpi(icut_dft)
            else if(ytext.eq.'cutd')then
               call inpi(icutd_dft)
            else if(ytext.eq.'tol ')then
               call inpi(itol_dft)
            else if(ytext.eq.'told')then
               call inpi(itold_dft)
            else
               omore = .false.
               jrec = jrec - 1
            endif
         enddo

         ierror=CD_inttol(ncent,itol_dft,icut_dft,itold_dft,icutd_dft)

      else if(ytext.eq.'disp') then
         if (jrec .lt. jump) then
           call inpa4(ytext)
           if (.not.(ytext.eq."on".or.ytext.eq."off".or.
     +               ytext.eq."aver".or.ytext.eq."geom")) then
             jrec = jrec - 1
             ytext = "on"
           endif
         else
           ytext = "on"
         endif
         call vdwaals_basic_control(ytext)
      else if(ytext.eq.'debu') then
10132    call inpa(ztext)
         if(ztext .eq. '        ')then
            goto 20
         else
            ierror=CD_debug(ztext)
            if(ierror.eq.-1)then
c
c bad print spec (ambiguous)
c
               call caserr2('bad ccpdft print spec')
            else if(ierror.eq.1)then
c
c try and process as alternative key
c
               jrec = jrec -1
               goto 10131
            else
c
c  that was OK, check for another print key
c
               goto 10132
            endif
         endif
         if(ierror.ne.0.or.ierr0.ne.0.or.ierr1.ne.0)
     +      call caserr2('problem processing ccpdft key')
         goto 10132
      else
         write(iwr,*)'problem with keyword ',ytext
         call caserr2('bad ccpdft keyword '//ytext)
      endif
      goto 10131
c
 20   continue
      return
      end
      subroutine read_types(nat,zaname,atomss,natomss,ytyp)
      implicit none
      integer nat, natomss
      character*8 zaname(nat)
      character*4 ytyp
      integer atomss(nat)
c
c     Assuming we have hit the string "label" or "element" we want to 
c     read a number of atom-labels. Each label we read is translated 
c     into 1 or more indeces in the array zaname (1 label may occur for
c     various atoms).
c
c     It is assumed that the string "label" is followed by 1 or more
c     strings each of which exactly match at least 1 string zaname.
c
c     It is assumed that the string "element" is followed by 1 or more
c     strings that are names of chemical elements, each of these 
c     elements must be represented in zaname at least once.
c
c     The routine terminates if we hit a string that is not found in 
c     zaname, and is not equal "label" and not equal "element". 
c     On return natomss will be the number of of atoms identified.
c
      integer i, j, ielm, isubst
      logical found
      character*8 label
c
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
c
      natomss=0
      label = ytyp
 10   if (label(1:4).eq."grid") then
         call inpa(label)
         goto 10
      else if (label(1:4).eq."labe") then
 20      call inpa(label)
         found = .false.
         do i = 1, nat
            if(label.eq.zaname(i)) then
               do j = 1, natomss
                  if(atomss(j).eq.i) then
                     call caserr2('Atom specified more than once!!')
                  endif
               enddo
               natomss=natomss+1
               atomss(natomss)=i
               found = .true.
            endif
         enddo
         if (found) then
            goto 20
         else
            goto 10
         endif
      else if (label(1:4).eq."elem") then
 30      call inpa(label)
         ielm = isubst(label)
         found = .false.
         do i = 1, nat
            if(ielm.eq.isubst(zaname(i))) then
               do j = 1, natomss
                  if(atomss(j).eq.i) then
                     call caserr2('Atom specified more than once!!')
                  endif
               enddo
               natomss=natomss+1
               atomss(natomss)=i
               found = .true.
            endif
         enddo
         if (found) then
            goto 30
         else
            goto 10
cDEBUG
c           if (ielm.eq.-10) then
c              goto 10
c           else 
c              return
c           endif
cDEBUG
         endif
      endif
c
      if (label.ne.' ') then
         jrec=jrec-1
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine reduce_gridlist(nat,grd_type,natomss,atomss,
     &                           zaname,iwr)
      implicit none
      integer nat, natomss
      integer grd_type(nat), atomss(natomss)
      character*8 zaname(nat)
      integer iwr
c
c     The arguments:
c
c     - nat: the number of atoms in the geometry.
c     - grd_type: contains the mapping that associates each atom with
c            1 grid type.
c     - natomss: the number of atoms referenced in atomss
c     - atomss: the atoms for which the grid specification will be
c            changed.
c
c     The task:
c
c     The subroutine performs the task that will be explained in the
c     example below.
c
c     Imagine we have a methane molecule then the initial situation
c     might be:
c
c     atom number 1  2  3  4  5
c     atom label  c  h1 h2 h3 h4
c     grd_type    1  2  2  2  2
c
c     imagine that we want to change the grid definition for h1 and h2
c     (i.e. natomss=2, atomss/2,3/).
c     To be able to do that without changing the grid definition of 
c     h3 and h4, we need to create a new grid type for atoms h1 and h2.
c     I.e. the above situation will be transformed to:
c
c     atom number 1  2  3  4  5
c     atom label  c  h1 h2 h3 h4
c     grd_type    1  3  3  2  2
c
c     To guarantee that any changes applied earlier to all hydrogen
c     atoms are carried over grid type 3 initially has to be an exact
c     copy of grid type 2.
c
c     However if in the initial situation atomss was to reference all
c     atoms of type 2 (i.e. natomss=4, atomss/2,3,4,5/) then no new
c     grid types need to be created.
c
c     Note: if grd_type(i).eq.-1 then this centre does not have a
c           grid! Any attempt to mess with this centre is an error!!!
c
      integer min1, max1, i, j, ni, nj, newtype, ierror
      integer CD_clone_grid, CD_assign_grid
c
      min1 = nat
      max1 = 0
      do i = 1, natomss
         min1 = min(grd_type(atomss(i)),min1)
         max1 = max(grd_type(atomss(i)),max1)
      enddo
c
c     min1: The minimum value for the grid type present in atomss
c           (min1=1 in the initial situation in the example)
c     max1: The maximum value for the grid type present in atomss
c           (max1=2 in the initial situation in the example)
c
      if (min1.eq.-1) then
10       format(5x,'*** centre',i6,' labelled ',a8)
         write(iwr,*)
         write(iwr,*)'ERROR: Attempt to change grid parameters for ',
     &             'gridless centres!!!'
         write(iwr,*)'The offending centres are:'
         do i = 1, natomss
            if (grd_type(atomss(i)).eq.-1) then
               write(iwr,10)atomss(i),zaname(atomss(i))
            endif
         enddo
         write(iwr,*)
         call caserr2('Cannot change grid of gridless centre!!!')
      endif

      do i = min1, max1
         nj = 0
         ni = 0
         do j = 1, natomss
            if (grd_type(atomss(j)).eq.i) then
               nj=nj+1
            endif
         enddo
c
c        nj: the number of times atomss refers to an atom with grid
c            type i
c
         if (nj.gt.0) then
            do j = 1, nat
               if (grd_type(j).eq.i) then
                  ni=ni+1
               endif
            enddo
         endif
c
c        ni: the number of times grid type i occurs in grd_type
c
         if (nj.lt.ni) then
c
c           the set of atoms with grid type i that atomss refers to 
c           is smaller than the set of all atoms with grid type i,
c           therefore we need to create a new type for the atoms 
c           referenced in atomss
c
            newtype = CD_clone_grid(i)
            do j = 1, natomss
               if (grd_type(atomss(j)).eq.i) then
                  grd_type(atomss(j))=newtype
                  ierror=CD_assign_grid(atomss(j),newtype)
                  if (ierror.ne.0) then
                     write(iwr,*)'CD_assign_grid: ',atomss(j),newtype
                     call caserr2('CD_assign_grid failed!!!')
                  endif
               endif
            enddo
         else if (nj.eq.ni) then
c
c           atomss references all atoms of grid type i so no need to
c           do anything
c
         else if (nj.gt.ni) then
c
c           atomss references more atoms of grid type i than there 
c           exist, clearly something is wrong!
c
            call caserr2('Reduce_gridlist messed up!!!')
         endif
      enddo
c     
      end
c
c-----------------------------------------------------------------------
c
      subroutine new_gridlist(nat,grd_type,natomss,atomss,
     &                        zaname,iwr)
      implicit none
      integer nat, natomss
      integer grd_type(nat), atomss(natomss)
      character*8 zaname(nat)
      integer iwr
c
c     The arguments:
c
c     - nat: the number of atoms in the geometry.
c     - grd_type: contains the mapping that associates each atom with
c            1 grid type.
c     - natomss: the number of atoms referenced in atomss
c     - atomss: the atoms for which the grid specification will be
c            changed.
c
c     The task:
c
c     The subroutine will assign grid types to specified atoms 
c     that do not have one yet.
c
c     Note: if grd_type(i).ne.-1 then this centre does have a grid
c           already! Any attempt to mess with this centre is an error!!!
c
      logical oerror
      integer i, j, ierror
      integer CD_create_grid, CD_assign_grid, isubst
c
      oerror = .false.
      do i = 1, natomss
         oerror = oerror.or.(grd_type(atomss(i)).ne.-1) 
      enddo
      if (oerror) then
10       format(5x,'*** centre',i6,' labelled ',a8)
         write(iwr,*)
         write(iwr,*)'ERROR: Attempt to create grid parameters for ',
     &             'centres that have a grid already!!!'
         write(iwr,*)'The offending centres are:'
         do i = 1, natomss
            if (grd_type(atomss(i)).ne.-1) then
               write(iwr,10)atomss(i),zaname(atomss(i))
            endif
         enddo
         write(iwr,*)
         call caserr('Cannot create grid for gridded centre!!!')
      endif
c
c     Now create new grid type by element and store them in grd_type
c
      do 20 i = 1, natomss
         do j = 1, i-1
            if (isubst(zaname(atomss(j))).eq.isubst(zaname(atomss(i))))
     &      then
c
c              atomss(i) is of the same element as atomss(j) so
c              reuse the gridtype of atomss(j)
c
               grd_type(atomss(i)) = grd_type(atomss(j))
               goto 20
            endif
         enddo
         grd_type(atomss(i)) = CD_create_grid()
 20   continue
c
c     Associate the specified atoms with the new grid types
c
      do i = 1, natomss
         ierror=CD_assign_grid(atomss(i),grd_type(atomss(i)))
         if (ierror.ne.0) then
            write(iwr,*)'CD_assign_grid: ',atomss(i),grd_type(atomss(i))
            call caserr('CD_assign_grid failed!!!')
         endif
      enddo
c     
      end
c
c-----------------------------------------------------------------------
c
      subroutine store_spare(imemory)
      implicit none
c
c     This subroutine and entry point are meant to store the amount of
c     spare memory for the incore coulomb fit. 
c
c     The problem is that for the incore coulomb fit one has to allocate
c     the memory to store the 2- and 3-centre integrals quite high up
c     in the call tree. To establish how much memory we can use for this
c     storage we estimate the amount of memory needed below this point.
c     At the moment however we cannot determine the exact amount needed.
c     Therefore there is a chance that we get it wrong and the 
c     calculation runs out of memory. 
c
c     In the above case where a calculation runs out of memory we allow
c     for the user to set aside some "spare memory". I.e. the user can 
c     take away some memory from the integral storage, in the hope that
c     the calculation will then be able to complete succesfully. 
c
c     The current 2 routines are meant to store this amount of spare
c     memory in words for each node.
c
      integer imemory
c
      integer ispare_memory
      save    ispare_memory
      data ispare_memory/0/
c
      if (imemory.lt.0) then
         call caserr2('The amount of spare memory must be non-negative')
      endif
      ispare_memory = imemory
      return
c
      entry retrieve_spare(imemory)
      imemory = ispare_memory
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine defvect(oinv,moutaa,moutbb,orest,oresti)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     this routine attempts to apply a sensible and consistent
c     set of defaults for sections specified on both the vectors
c     and enter directive. These are based on the following
c     settings
c                           Alpha (1st section)     Beta (2nd section)
c   ==================================================================
c     closed shell SCF mos  isect(301) = 1                -
c     UHF SCF mos           isect(302) = 2          isect(303) = 3
c     open shell GVB mos    isect(304) = 4          isect(305) = 5
c     CASSCF MOs            isect(306) = 6          isect(307) = 7
c     MCSCF/MULTI Mos       isect(308) = 8          isect(309) = 9
c   ==================================================================
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      dimension zscf(12)
      data zscf/'rhf','uhf','gvb','grhf','casscf','mcscf','multi',
     +          'direct','mp2','mp3','vb','masscf'/
      data zguessg, zguessv/ 'atoms', 'mosaved'/
c
c     has the vectors section been omitted i.e. defaults to apply
c
      if (oinv) go to 12345
c
       if (orest) then
c
        if (oresti) then
c
c     2nd pass through data input ? 
c     (orest = .true., oresti=.true.)
c
         zguess = zguessv
         mina = mouta
         if (zscftp.eq.zscf(2)) then
            if (moutb.le.0) then
               minb = mouta
            else
               minb = moutb
            end if
         end if
         if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4)) then
            if (moutb.le.0) then
               minb = mouta
            else
               minb = moutb
            end if
         end if
         if (zscftp.eq.zscf(5) .or. zscftp.eq.zscf(6).or.
     1      zscftp.eq.zscf(11)) then
            minb = 0
         end if
c
c     is this the first step in a restart job? 
c     (orest = .true., oresti=.false.)
c
        else
c
         zguess = zguessv
         if(zscftp.eq.zscf(1)) then
          mina = isect(301)
         else if(zscftp.eq.zscf(2)) then
          mina = isect(302)
         else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
c ... use canonicalised mo section for gvb calcs
          mina = isect(305)
         else if(zscftp.eq.zscf(5)) then
          mina = isect(306)
         else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7).or.
     1   zscftp.eq.zscf(11)) then
          mina = isect(308)
         else
         endif
c ... does section section exist ?
         if (mina.gt.0) then
            call sectst(mina,i)
         else
            i = 0
         end if
         if (i.eq.0) then
c ... no default section from previous scf of same type
c ... try the closed shell section, if this is absent, then
c ... revert to atomic guess as last resort
          mina = isect(301) 
          call sectst(mina,i)
          if (i.eq.0) zguess = zguessg
          go to 12345
         endif
c ... now deal with minb
         if (zscftp.eq.zscf(2) .or. zscftp.eq.zscf(3) .or.
     +       zscftp.eq.zscf(4) .or. zscftp.eq.zscf(5) .or.
     +       zscftp.eq.zscf(6) .or. zscftp.eq.zscf(7) .or.
     +       zscftp.eq.zscf(11)) then
            if(zscftp.eq.zscf(2)) then
             minb = isect(303)
            else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
c ... use canonicalised mo section for gvb calcs
             minb = isect(305)
c ... no 2nd section for casscf or mcscf calculations
            else if(zscftp.eq.zscf(5)) then
             minb = 0
            else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7)) then
             minb = 0
            else
            endif
            if(minb.ne.0) then
              call sectst(minb,i)
c ... if beta section does not exist, set minb = mina
              if (i.eq.0) minb = mina
            endif
         endif
c
       endif
c
c      startip job (orest and oresti = .false.)
c
c ... default options now replaced by 'atoms'
c
       else if (oming) then
            zguess = zguessg
       else if (oextg) then
            zguess = zguessg
       else if (zscftp.eq.zscf(1)) then
            zguess = zguessg
       else
            zguess = zguessg
c
       end if
c
c     now deal with sections formally assigned to =enter=
c     specify default output vectors sections if =enter= omitted
c
12345 continue
c
c     identify case where second section may be required
c
      osecond = zscftp.eq.zscf(2) .or. zscftp.eq.zscf(3) .or.
     +          zscftp.eq.zscf(4) .or. zscftp.eq.zscf(5) .or.
     +          zscftp.eq.zscf(6) .or. zscftp.eq.zscf(7) .or.
     +          zscftp.eq.zscf(11) .or. zscftp.eq.zscf(1)
c
      if (oresti) then
c
c     note that some of the logic here is probably unnecessary
c     i.e. the same logic would apply whether oresti is true
c     or false, but repeat for the moment ...
c
         if (moutaa.le.0) then
c
c        runtype analyse is a special case; set mouta = mina
c        in default as only one set of MOS will be analysed
c
            if(zscftp.eq.zscf(1)) then
            mouta = isect(301)
            else if(zscftp.eq.zscf(2)) then
            mouta = isect(302)
            else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
            mouta = isect(304)
            else if(zscftp.eq.zscf(5)) then
            mouta = isect(306)
            else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7).or.
     1      zscftp.eq.zscf(11)) then
            mouta = isect(308)
            else
            endif
         else
            mouta = moutaa
         end if
         if (osecond) then
            if (moutbb.gt.0) then
               moutb = moutbb
            else if(zscftp.eq.zscf(2)) then
               moutb = isect(303)
            else if(zscftp.eq.zscf(5)) then
               moutb = isect(307)
            else if(zscftp.eq.zscf(11)) then
               moutb = mouta
            else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
                if(moutaa.gt.0) then
                 moutb = mouta
               else
                 moutb = isect(305)
               endif
            else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7).or.
     1        zscftp.eq.zscf(11)) then
               moutb = isect(309)
            else
               moutb = moutaa
            end if
         end if
c
c        runtype analyse is a special case; set mouta = mina
c        if no output vectors specified (but need to allow
c        for analysis modules producing new vectors e.g. lmo)
c
         if(zruntp.eq.'analyse'.and.oinv) then
          mouta = mina
          if(osecond) moutb = mina
         endif
      else
         if (moutaa.le.0) then
c
c        runtype analyse is a special case; set mouta = mina
c        in default as only one set of MOS will be analysed
c
            if(zruntp.eq.'analyse') then
            mouta = mina
            else if(zscftp.eq.zscf(1)) then
            mouta = isect(301)
            else if(zscftp.eq.zscf(2)) then
            mouta = isect(302)
            else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
            mouta = isect(304)
            else if(zscftp.eq.zscf(5)) then
            mouta = isect(306)
            else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7).or.
     1       zscftp.eq.zscf(11)) then
            mouta = isect(308)
            else
            endif
         else
            mouta = moutaa
         end if
        
         if (osecond) then
            if (moutbb.gt.0) then
               moutb = moutbb
            else if (zruntp.eq.'analyse') then
               moutb = mouta
            else if(zscftp.eq.zscf(2)) then
               moutb = isect(303)
            else if(zscftp.eq.zscf(5)) then
               moutb = isect(307)
            else if(zscftp.eq.zscf(11)) then
               moutb = mouta
            else if(zscftp.eq.zscf(3).or.zscftp.eq.zscf(4)) then
                if(moutaa.gt.0) then
                 moutb = mouta
               else
                 moutb = isect(305)
               endif
            else if(zscftp.eq.zscf(6).or.zscftp.eq.zscf(7).or.
     1         zscftp.eq.zscf(11)) then
               moutb = isect(309)
            else
               moutb = moutaa
            end if
         end if
      end if
c
      return
      end
      subroutine start2(core,jrun,jdeg,
     * olevel,oactiv,ocidu,ogvb,ogrhf,orgvb,
     * odci,obasis,osafe,odft,odftgrad,nav)
c
c    now check and print input options
c    then initiate computation
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *35 titdf
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
cdrf
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
c
       integer idafh, navh, ioda
       common /hdafile/ idafh,navh,ioda(2,1000)
c
      common/nottwi/obeen,obeen2,obeen3,obeen4
      logical vector
      character*8 scftp
      common/scfopt2/scftp
cdrf
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      dimension core(*)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
      real*8  dx0, dxm, dele, drms
      integer ntpr, ntrun, imcyc, ncyc, ntmin, nspacb
      common /modmin/ ntpr,ntrun,imcyc,ncyc,ntmin,nspacb,
     +                dx0,dxm,dele,drms
c
c
      integer nword, isw, kprint, ifzmat
      common /mod2d/ nword,isw,kprint,ifzmat
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
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/blksiz2/ns2
c
      real*8 cccnv
      logical oconv, oconvf
      integer icanon
c
      common /gms_finish/ cccnv,oconv,icanon,oconvf

      common/prints/oprint(20),
     * odist, obass, oadap, osadd, ovect, oanal, ohess,
     * ocinpr, omrdcinpr, occsdnpr, onobqp, onopr(9),
     * opmull, oplowd, opmos , opsadd, opvect, opadap, opdiis,
     * opsymm, opfock, opq   , opao  , opmo  , opdump, opgf  ,
     * optda , opdist, opdire, opcrp0, opinte, opcrp1        ,
     * ciprnt

      common/restrl/ociopt,ocifor,omp2,ohfg(12),omcscf,oforce,oci,
     +              ogr(32),oatmdft
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
      integer n2file, n2tape, n2blk, n2last
      integer n4file, n4tape, n4blk, n4last
      integer n6file, n6tape, n6blk, n6last
      integer n5file, n5tape, n5blk, n5last
      integer n9file, n9tape, n9blk, n9last
      integer ntfile, nttape, ntblk, ntlast
      integer n1file, n1tape, n1blk, n1last
      integer n11fil, n11tap, n11bl, n11lst
      integer n12fil, n12tap, n12bl, n12lst
      integer n13fil, n13tap, n13bl, n13lst
      integer maxbl, minbl, maxset
      common/gms_files/
     +      n2file,n2tape(20),n2blk(20),n2last(20),
     +      n4file,n4tape(20),n4blk(20),n4last(20),
     +      n6file,n6tape(20),n6blk(20),n6last(20),
     +      n5file,n5tape(20),n5blk(20),n5last(20),
     +      n9file,n9tape(20),n9blk(20),n9last(20),
     +      ntfile,nttape(20),ntblk(20),ntlast(20),
     +      n1file,n1tape(20),n1blk(20),n1last(20),
     +      n11fil,n11tap(20),n11bl(20),n11lst(20),
     +      n12fil,n12tap(20),n12bl(20),n12lst(20),
     +      n13fil,n13tap(20),n13bl(20),n13lst(20),
     +       maxbl(maxlfn),minbl(maxlfn),maxset(maxlfn)
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
      common/data1/vlist(100,3),cenmas(100),
     * newpro,nocp,mcen,mspace,itemp(100),jtemp(100)
     * ,nacta,nactb,louta(maxorb),loutb(maxorb)
     * ,norbt,norbta(maxorb),norbtp,
     *  norbc,norbca(maxorb),norbcp,ihamcp,ihamc
     * ,norb2,nseco(maxorb),norb3,nthird(maxorb),norbg,norbga(maxorb),
     *  norba,norbaa(maxorb),
     *  omulla,omullo,omullg,omullp,omullu,mullac(maxorb),mulla,oactrn
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
      common/junkc/ylab(26),ztype(35),zlist(100),zcentg(100),
     * ztagc(100),znumb(100),zhead(5,22),zunit(286),
     * ztita(10),zatagg(maxat),zmul(maxorb)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/scfopt/maxcyc,mconv,nconv,npunch,
     * accdi1,accdi2,odiis,icoupl,ifolow,irotb,
     * dmpcut,acccc(6),iaccv(2),
     * rshift(4),iextn,junks,dampoo(9),idiisf,
     * optester,odynamic,omaxcyc,onoor,otestd
      common /copt / rcvopt
      common/foropt/vibsiz,nvib
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      real*8 acc, stol, stepmx, tolmax, tolstp, tanstp
      real*8 accin, eigmin, eigmax, grms
      integer jm, mnls, nls, isadle, ifcm, iblfcm, isecfcm
      integer nentry, iterch, nftotl, ismax, lintyp, lqstyp
      logical ofcm

      common /cntl1/ acc,stol,stepmx,tolmax,tolstp,tanstp,
     +               jm,mnls,nls,isadle,imnter,idumcntl1,
     +               accin(6),nftotl,ismax,lintyp,lqstyp,
     +               eigmin,eigmax,grms(2),nentry,iterch,
     +               ifcm,iblfcm,isecfcm,ofcm
c
c
      real*8 fieldx, fieldy, fieldz
      logical ofield
      common /gms_field/  fieldx, fieldy, fieldz, ofield
c
c
      integer nint, next, ninnex, nwm1, nw, nwp1, nintwo
      integer nint21, nintp1
      common /orb/ nint,next,ninnex,nwm1,nw,nwp1,nintwo,
     +             nint21,nintp1
c
c
      real*8  vg, gxp, gyp, gzp, d2, valmo, occa, occb
      real*8  rij, timat, denat, encoat, eltot
      integer idenfu,nangpr, ncorxx, ndeg, nomo
      logical odenfu, osymm, osc
      common /dnfnw/ odenfu, osymm, osc, idenfu, nangpr, ncorxx,
     +               vg(200),gxp(200),gyp(200),gzp(200),d2(200),
     +               valmo(200),occa(200),occb(200),
     +               rij(1275),timat(50),denat(50),encoat(50),ndeg(50),
     +               eltot,nomo
c
c
c      following common block for rfcontinuum scf calcs
c      contains input data
c
      real*8 aradix,aradiy,aradiz, dielec
      logical oscrf
      common/scrf/aradix,aradiy,aradiz,dielec,oscrf
c
c ----- common blocks for psscrf calculation ----
c
      real*8 ptspac, delect, facnrm
      integer iconv,inkblk,itotbq,iscf,ispace,msecp
      logical opssc
      common/psscrf/ptspac,delect,iconv,opssc,
     +              facnrm,inkblk,itotbq,iscf,ispace,msecp
      common/surfgn/numbq(maxat),irejec(maxat),cavsze(maxat)
     * ,nocent,itmax,damper
c
c - xfield..
c
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
c
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
c   /gen/ and /gen2/ contain information about the direct-scf options
c
      common/gen/maxitd,iterd,concrt,convd,energd,thrshd,dcmax
     +,olog1,olog2,olog5,olog9,olog6,odbg(10),omp2d,mp2fo,mp2fv
     +,oreadm,owritm,orestd,iaccur,ocvary,ogrdt,ogopt,ogdma
     +,ogdold
     +,ogscf,noptd,npoint,timeg,ogmull,iswapd,dshift,irhess
      common/gen2/ripftd(21,2),rijkld(21,2),rjpftd(21,2),
     + ripfal(21,2),rijkal(21,2),
     + prefac(21),tinner(21),rifall,riffew,ricfb1,ricfb2
c
      logical oentry
      integer idblk, numdir, inddir
      common /wrtd/ idblk(30),numdir,inddir,oentry
c
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c
c     Variables controlling the Fermi-Dirac smearing
c     Search for the subroutine fermi_smear...
c
c     osmear: should Fermi-Dirac smearing be used?
c     esmear: smearing parameter controlling the width of the Fermi 
c             cut-off function (Hartree)
c
c     esmear_start: the initial smearing parameter
c     esmear_final: the final smearing parameter
c     esmear_scale: the smear parameter follows the "tester" scaled
c                   with this factor from the start to the final value.
c
      logical osmear
      real*8    esmear_start, esmear_final, esmear_scale
c
      common/fermidirac/esmear_start,esmear_final,esmear_scale,
     &                  osmear
c
c
c     Constants for converting between atomic units and other units
c
      real*8 bohr_to_ang,   ang_to_bohr
      real*8 hartree_to_eV, eV_to_hartree
      parameter (bohr_to_ang   = 0.529177249d0)
      parameter (ang_to_bohr   = 1.0d0/bohr_to_ang)
      parameter (hartree_to_eV = 27.21165d0)
      parameter (eV_to_hartree = 1.0d0/hartree_to_eV)
     

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
c
      real*8 scalel, vbond1, vbond2
      logical opoint1,ochkdst,ofatal
      integer mxbnd, nbond1, nbond2, mbond1, mbond2
      parameter (mxbnd = 1000)
c
      common/structchk/ vbond1(mxbnd),mbond1(2,mxbnd),
     +                  vbond2(mxbnd),mbond2(2,mxbnd),
     +                  scalel, nbond1, nbond2, opoint1,
     +                  ochkdst, ofatal
c
c
c
c     MASSCF common fccwfn 
c
      integer nspace
      integer msta,mnum,mini,maxi
      integer iami,iama,ibmi,ibma,idim
      integer lbst,nref0
      logical omas,omasd
      logical fdirct,qcorr
      real *8 c0sq
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq,
     *                omas, omasd
c
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
      integer ncore_f,nact_f,nels_f
      common /fccwfn2 / ncore_f,nact_f,nels_f
      parameter (mxnoro=20, mxfrz=20)
      logical fcore,osrso
      common /masinp/ toleng,acurcy_mas,damp_mas,mxpn,mxit,maxmas,
     *       norot,nrotab(2,mxnoro),nfrz,mofrz(mxfrz),fcore,osrso
      character*4 tmini, tmaxi, tmnum, tmsta
      data tmini, tmaxi, tmnum, tmsta /
     +      'mini', 'maxi', 'mnum', 'msta' /
c
      dimension mapcie(maxorb),icorec(maxorb)
      equivalence (mapcie(1),louta(1)),(icorec(1),loutb(1))
      dimension zscf(12),oflg1(20),oflg2(20)
      parameter (nvec=10)
      dimension zvect(nvec+1)
      dimension ztda(113)
      equivalence (ztda(1),ztita(1))
      equivalence (oflg1(1),odist),(oflg2(1),opmull)
c
      data ono/.false./
      data zscf/'rhf','uhf','gvb','grhf','casscf','mcscf','multi',
     &          'direct','mp2','mp3','vb','masscf'/
      data zatext,zbtext/'a','b'/
      data zvect/'minguess','extguess','hcore','getq','extra',
     *'alphas','cards','nogen','atoms','mcards', 'mosaved'/

      data m4,m23/4,23/
      data m3,m9/3,9/
      data m1,m2/1,2/
      data zsort,zpsort/'sort','psort'/
      data mit/200/
      data maxbit /64/
c
c  Print out nuclear masses if not the defaults
c
      if(omass_mod())call mass_print
c
c
      if(CD_active())then
c
c make some changes for the DFT cases
c
      idum = CD_set_2e()

      if( CD_HF_exchange() .or. CD_HF_coulomb() )then
c
c for DFT, coulomb and exchange operators must be separable
c so disable use of supermatrix 
c         
         if(nopk .ne. 1)then
            if(opg_root())then
               write(iwr,*)'Supermatrix option switched off for DFT'
            endif
            nopk = 1
         endif
      else
c
c We are doing a DFT calc without the need for 2-electron integrals,
c switch to Direct mode. This avoids extra logic in the SCF drivers
c and enables existing atomic SCF guess to work 
c (eventually should use DFT code in guess as well)
c
         if(.not.odscf)then
            if(opg_root())write(iwr,*)
     &           '2-electron integrals not required for requested DFT ',
     &           'options, switching to direct SCF mode'
            odscf = .true.
         endif
      endif
      idum = CD_reset_2e()
      endif
c
chvd  if (jrun.ge.4 .and. jrun.le.9) then
c     the above line is stupid. The DFT code always checks whether
c     derivatives of the weights are required for the current term.
c
c   if default gradquad setting (.false.) and directive not invoked, 
c   then reset gradquad to true if 1st or 2nd row transition metal 
c   involved in gradient evaluation
c
         otm = .false.
         do i = 1, nat
         nucz = isubst(zaname(i))
         if ( (nucz.ge.21.and.nucz.le.30) .or.
     +        (nucz.ge.39.and.nucz.le.48) )  otm = .true.
         enddo
c
         if (.not.odftgrad.and.otm) then 
           ierror=CD_gradquad(.true.)
         endif
chvd  endif
c     eventually make dft atoms the default 
c     if dft has been requested
c     if (odft) then
c      if(zguess.eq.'atoms') then
c       oatmdft = .true.
c      endif
c     endif

      if (obasis) then
c
c     don't need this block of code for a mechanics or mopac run.
c     normally obasis should be true if this part of the code
c     is reached successfully. thus, for a mechanics or mopac run obasis
c     should be false as no basis set is required.
c
         if (itol.eq.0) itol = 20
         if (icut.eq.0) icut = 9
         if (opsymm) write (iwr,6010) toler , toler2
         if (osafe) write (iwr,6040) safe , dumtim
         if (otran) write (iwr,6030)
         if (ofield) write (iwr,6020) fieldx , fieldy , fieldz
         if(odenfu) then
          if(idenfu.eq.1) titdf='Lee-Yang-Parr Formula              '
          if(idenfu.eq.2) titdf='Lee-Yang-Parr Formula (2nd order)  '
          if(idenfu.eq.3) titdf='Vosko-Wilk-Nusair Formula          '
          if(idenfu.eq.4) titdf='Vosko-Wilk-Nusair Formula + SIC    '
          nang=72
          if(nangpr.eq.2) nang=128
          if(nangpr.eq.3) nang=194
          write(iwr,7010) titdf,nang
          if(osymm) write(iwr,7020 )
          if(ncorxx.ne.0) then
           icores=ncorxx*2
           write(iwr,7030) icores
          end if
         end if
c      write out scrf input data
       if(oscrf) write(iwr,7060)dielec,aradix,aradiy,aradiz
c
c  -----flag xfield setting
       if(ocryst) write(iwr,7065)
c
c   ----- write out psscrf data -----
       if (opssc) then
          write(iwr,7070)delect,ptspac,nocent
          do 360 ii=1,nocent
           write(iwr,7080)ii,cavsze(ii)
360       continue
          write(iwr,7090)itmax,damper,msecp
       end if
      end if
c
      if (osafe .and. obasis) then
         write (iwr,6040) safe , dumtim
      else if (osafe) then
         write (iwr,6050) safe
      end if
      if (ogvb .or. old .or. ogrhf) then
         call gvbin
         if (ogrhf) call gvbrhf(.true.)
      end if
      if (orgvb) then
         call secget(isect(504),m23,k)
         k = k + 2
         call rdedx(cicoef,mach(5),k,idaf)
      end if
      acurcy = 10.0d0**(-nconv)
      concrt = acurcy
      convd = acurcy/10.0d0
      cccnv = acurcy
      accdi2 = acurcy*acurcy
      if (omrdci .and. oharm) then
       call caserr2(
     +  'harmonic basis not functional in mrdci calculations')
      endif
      if ((jrun.eq.3 .or. (jrun.eq.5.and.ocifp))
     +         .and. .not.ofulci .and. .not.omrdci
     +         .and. .not.occsd  .and. .not.oqcisd
     +         .and. .not.opark  .and. .not.odiesel ) then
         if (.not.odci) then
            call ciin(core,.true.)
            odci=.true.
         end if
      end if
      write (iwr,6060) zruntp
      if (jrun.eq.9 .or. jrun.eq.15 .or. jrun.eq.41 .or. (ntmin.eq.1)
     +    .or. (ntmin.eq.2)) go to 50
      zcom(5) = zscftp
      if (zguess.ne.zvect(nvec+1)) go to 20
      if (mina.le.0 .or. mina.gt.508) then
         call caserr2('invalid dumpfile section specified')
c
      else if (zscftp.ne.zscf(5) .and. zscftp.ne.zscf(6) .and.
     &         zscftp.ne.zscf(11) ) then
         if (zscftp.ne.zscf(2)) go to 20
         if (minb.eq.0) minb = mina
      else
         if (minb.lt.0) minb = 0
      end if
      if (minb.lt.0 .or. minb.gt.508) then
         call caserr2('invalid dumpfile section specified')
      end if
 20   if (mouta.le.0 .or. mouta.gt.508) then
         call caserr2('invalid dumpfile section specified')
      else if (zscftp.ne.zscf(2)) then
         if (zscftp.ne.zscf(5) .or. zscftp.ne.zscf(11)) go to 30
      end if
      if (moutb.le.0 .or. moutb.gt.508) then
         call caserr2('invalid dumpfile section specified')
      end if
 30   write (iwr,6070) zscftp , zguess
      if (zguess.eq.zvect(1) .or. zguess.eq.zvect(2)) then
         if (scaleg.ne.0.0d0) then
            write (iwr,6080) zguess , scaleg
         end if
      end if
      if (zguess.eq.zvect(10)) then
         write (iwr,6090) zatext , mina
         if (zscftp.eq.zscf(2)) write (iwr,6090) zbtext , minb
      end if
      write (iwr,6100) zatext , mouta
      if (zscftp.eq.zscf(2)) then
         write (iwr,6100) zbtext , moutb
         if (mouta.eq.moutb.and.jrun.ne.10) call caserr2(
     +    'identical section specified for -a- and -b- uhf orbitals')
      endif
      if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4)) then
         if (moutb.le.0) moutb = mouta
         write (iwr,6110) moutb
         if (moutb.gt.508) then
            call caserr2('invalid dumpfile section specified')
         end if
      end if
      if (omas) then
       write(iwr,3310) ncore_f,nact_f,nels_f,
     +                 toleng,acurcy_mas,damp_mas,maxmas,
     +                 norot,nfrz,mxpn,mxit,fcore,osrso
       if (osrso) write(iwr,3315)
 3310 format (/1x,
     +        '**************************************************** '/
     +     1x,'*             MASSCF Specifications                * '/
     +     1x,'**************************************************** '/
     +     1x,'** ncore = ', i2, '  ** nact = ', i2, ' ** nels = ', i2/
     +     1x,'* energy convergence tolerance = ', e8.2,/,
     +     1x,'* newton-raphson and davidson convergence = ', e8.2,/,
     +     1x,'* newton-raphson damping factor = ', e8.2,/,
     +     1x,'* maximum number of masscf iterations = ', 1i3,/,
     +     1x,'* number of frozen rotations = ', 1i3,/,
     +     1x,'* number of frozen orbitals = ', 1i3,/,
     +     1x,'* maximum number of davidson expansions = ', 1i3,/,
     +     1x,'* maximum number of davidson iterations = ', 1i3,/,
     +     1x,'* option to freeze all core orbitals = ', 1l1,/,
     +     1x,'* perform single-ref second order calc = ', 1l1)
 3315  format(
     +     1x,'*      a single scf iteration will be computed ',/,
     +     1x,'*      to initialize the second-order method ')
       write(iwr,3320) zscftp, nspace
 3320  format(1x,'* SCFTYPE = ',a8/
     +     1x,'* number of active spaces = ', i3)
       write(iwr,3330) tmini,(mini(loop),loop = 1, nspace)
       write(iwr,3330) tmaxi,(maxi(loop),loop = 1, nspace)
       write(iwr,3330) tmnum,(mnum(loop),loop = 1, nspace)
       write(iwr,3330) tmsta,(msta(loop),loop = 1, nspace)
 3330  format (1x,'* ',a4,5x,10i4)
      endif
      if (oharm) write(iwr,6071)
      if (iprop.lt.0) write(iwr,6072)
      if (npunch.ne.0) then
         write (iwr,6120) npunch
      end if
      if (odscf) then
         if (odnew) then
c
c     using dumpfile for storage of most num*num arrays.
c
            numdir = idaf
            lend2 = lensec(num*num)
            lendir = 30*lend2
            call secput(isect(507),507,lendir,iblkdr)
            idblk(1) = iblkdr
            do 40 i = 2 , 30
               idblk(i) = idblk(i-1) + lend2
 40         continue
            if (omp2d) then
               write (iwr,6130)
            else
               write (iwr,6140)
            end if
         else
            write (iwr,6150) dexp(-tolitr(1)),itrtol(1),
     +                       dexp(-tolitr(2)),itrtol(2),
     +                       dexp(-tolitr(3)),itrtol(3),
     +                      dexp(deltol) , dexp(delfac) , oimag
c is delta density ever to be invoked? if not, we can
c reduce i/o activity around section 471 ..
c
            odelta = deltol .gt. dlog(acurcy)
         end if
      end if
 50   if (ltask.ne.0) then
         if (ltask.ne.jrun) call caserr2(
     +     'tasks nominated by restart and runtype differ')
      end if
      mtask = jrun
c ...
c skip the rest of this routine if this a molecular mechanics
c or mopac run.
c ...
      if (omopac) then
         write (iwr,6160)
         go to 340
      end if
      if ((ntmin.eq.1) .or. (ntmin.eq.2)) then
         write (iwr,6170)
         go to 340
      else if (ntmin.eq.3) then
         write (iwr,6180)
         go to 340
      else
         write (iwr,6190)
      end if
c
      if (jrun.eq.5 .and. dabs(accin(1)-0.3d0).le.1.0d-7) then
        accin(1) = 0.6d0
      endif
      go to ( 60, 60, 60, 60, 60, 60, 60, 60,130, 70,
     +        70, 70, 60, 70,130, 70, 70, 70, 70, 70,
     +        60, 60, 60, 60, 60, 60, 60, 60, 60, 60,
     +        60, 60, 60, 60, 60, 60, 60, 60, 60, 60,
     +        60, 60, 60, 60, 60, 60, 60, 60, 60, 60) , jrun
 60   jconv = 0
      write (iwr,6200) nconv
      if (mconv.ne.5) then
         write (iwr,6210) mconv
      end if
      if (omaxcyc) write(iwr,6211) maxcyc
      if (osym_op) write(iwr,6215)
      if (dmpcut.ne.0.0d0) write (iwr,6220) dmpcut
      if (lock.eq.4) then
         write (iwr,6230)
      end if
      if (mod(mconv,2).eq.0) then
         write (iwr,6240)
         jconv = jconv + 1
      end if
      if (mod(mconv,8).gt.3) then
         call levsin(olevel,odft)
         jconv = jconv + 1
      end if
      if (mod(mconv,4).gt.1) then
         write (iwr,6320)
         jconv = jconv + 1
      end if
      if (mconv.ge.7) then
         if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4)) then
            write (iwr,6330)
            jconv = jconv + 1
         end if
      end if
      if (jconv.eq.0) write (iwr,6340)
c
      if (odiis) then
       if (odynamic) then
        write (iwr,6351) accdi2
       else 
        write (iwr,6350) accdi1 , accdi2
       endif
      endif
      if (numdis.gt.0) then
         call filec(yed(numdis),ibldis,itemp,irep)
         if (irep.ne.0) call caserr2(
     +          'scf control file has invalid attributes')
         write (iwr,6360) yed(numdis) , ibldis
      end if
      if (opdiis) write (iwr,6370)
      if (osmear) then
         write (iwr,6376)esmear_start,esmear_start*hartree_to_eV,
     &                   esmear_final,esmear_final*hartree_to_eV,
     &                   esmear_scale
      endif
      write (iwr,6380)
      if (oanil) write(iwr,6381)
c
 70   if (otran) opadap = ono
c ... flag printing requests
      do 80 loop = 1 , 20
         if (oflg1(loop) .or. oflg2(loop)) go to 90
 80   continue
      go to 100
 90   write (iwr,6390)
      if (opmull) write (iwr,6400)
      if (oplowd) write (iwr,6410)
      if (opmos) write (iwr,6420)
      if (opvect) write (iwr,6430)
      if (opadap) write (iwr,6440)
      if (opsadd) write (iwr,6450)
      if (opsymm) write (iwr,6460)
      if (opfock) write (iwr,6470)
      if (opq) write (iwr,6480)
      if (opao) write (iwr,6490)
      if (opmo) write (iwr,6500)
      if (opdump) write (iwr,6510)
      if (opdist) write (iwr,6520)
      if (opdire) write (iwr,6530)
      if (opinte) write (iwr,6540)
c
c     ----- noprin specifications
c
      if (odist) write (iwr,6550)
      if (obass) write (iwr,6560)
      if (oadap) write (iwr,6570)
      if (osadd) write (iwr,6580)
      if (ovect) write (iwr,6590)
      if (oanal) write (iwr,6600)
      if (ohess) write (iwr,6610)
      if (onobqp) write (iwr,6615)
      write (iwr,6190)
c
 100  if (zruntp.eq.zrunt(4) .or. zruntp.eq.zrunt(5) .or.
     +    zruntp.eq.zrunt(6) .or. zruntp.eq.zrunt(7) .or.
     +    zruntp.eq.zrunt(8)) then
         if (odscf .and. odnew) then
            if (omp2d) call caserr2(
     +     'direct-mp2 not available in gradient calculations')
            ogrdt = .true.
c***** do gradient calculations
            do 110 i = 1 , 21
               tinner(i) = tinner(i)*1.0d-2
 110        continue
            do 120 i = 1 , 21
               prefac(i) = prefac(i)*1.0d-2
 120        continue
            if (zruntp.ne.zrunt(7)) ogopt = .true.
         end if
         if (zruntp.ne.zrunt(7)) then
            write (iwr,6620)
            if (zruntp.eq.zrunt(4) .or. zruntp.eq.zrunt(5)) then
               write (iwr,6630)
               if (ifcm.gt.0) write (iwr,6640) yed(ifcm) , iblfcm
               if (ochkdst) write(iwr,6632)
               if (jm.le.0 .or. jm.gt.mit .or. mnls.le.0 .or.
     +             mnls.gt.mit) call caserr2(
     +             'invalid number of calculations specified in minmax')
            else if (zruntp.eq.zrunt(8)) then
               write (iwr,6650) nvib , vibsiz
            else
               write (iwr,6660) rcvopt
               if (ochkdst) write(iwr,6632)
            end if
            write (iwr,6670)
         end if
      end if
c
 130  if (.not.odscf .and. jrun.ne.10)
     +    call filchk(n2file,n2blk,n2last,n2tape,m1)
      if (opass3) then
c
c ----- reset internal file monitoring for requested mainfile
c
         m2file = n2file
         do 140 i = 1 , n2file
            m2tape(i) = n2tape(i)
            m2blk (i) = n2blk(i)
            m2last(i) = n2last(i)
 140     continue
      end if
      if (jrun.eq.10) go to 270
      if (zscftp.eq.zscf(5)) then
         call incas2(oactiv,ocidu,jdeg,jrun)
         go to 260
      else if (zscftp.eq.zscf(6)) then
         call incas3(jrun)
         go to 260
      else
c
         if (omp2 .or. mp3) then
c
c     now set-up =active= and =core= specifications (mp2 / mp3)
c
c     active list
c
            if (oactrn) then
c     trap mp2 gradients with frozen or discarded orbitals
c     this is a bug that needs fixing ...
             if (jrun.ge.4.and.jrun.le.8.or.jrun.eq.13.or.
     +       jrun.eq.14. or.jrun.eq.18.or.jrun.eq.24) call caserr2(
     +      'frozen/discarded MOs invalidate mp2 gradients')
            else
c...         if this bug is fixed harmonic may use reduced # orbitals
             if (jrun.ge.4.and.jrun.le.8.or.jrun.eq.13.or.
     +       jrun.eq.14. or.jrun.eq.18.or.jrun.eq.24) norbt = num
             do 150 i = 1 , norbt
             norbta(i) = i
 150         continue
            end if
            nsa4= norbt
            nsb = nsa4
            ncoorb = norbt
c...        above statement seems to fix active problem in mp3
            if (nsa4.lt.1 .or. nsa4.gt.num) then
               call caserr2('invalid number of active orbitals')
            else
               do 160 i = 2 , nsa4
                  if ((locat1(norbta,i-1,norbta(i))).ne.0) call caserr2(
     +                'invalid orbital specified in active list')
 160           continue
               do 170 i = 1 , nsa4
                  j = norbta(i)
                  mapie(i) = j
                  ilifm(i) = ilifq(j)
 170           continue
            end if
            write (iwr,6680)
            write (iwr,6690) (mapie(k),k,k=1,nsa4)
c
c     frozen list
c
            ncore = norbc
            if (ncore.lt.1) then
               write (iwr,6700)
            else
               if (ncore.gt.num) call caserr2(
     +             'invalid number of functions in core list')
               do 180 i = 1 , ncore
                  j = norbca(i)
                  mapcie(i) = j
                  icorec(i) = ilifq(j)
 180           continue
               if (ncore.gt.1) then
                  do 190 i = 2 , ncore
                     if ((locat1(mapcie,i-1,mapcie(i))).ne.0)
     +                call caserr2(
     +                'invalid function specified in frozen core list')
 190              continue
               end if
               write (iwr,6710)
               write (iwr,6690) (mapcie(k),k,k=1,ncore)
            end if
            if (omp2 .and. omp2d) then
               mp2fo = ncore + 1
               mp2fv = ncore + nsa4
            end if
            if ((omp2 .and. .not.omp2d) .or. mp3) then
               call filchk(n4file,n4blk,n4last,n4tape,m2)
               call filchk(n6file,n6blk,n6last,n6tape,m3)
               acc1 = 10.0d0**(-iacc4)
               write (iwr,6720) acc1
               call setbfa
               write (iwr,6730) zsort , nsz
               go to 250
            end if
         end if
c
         go to (270,220,200,270,265,270,270,270,220,270,
     +          220,220,270,270,270,220,220,270,270,270,
     +          210,220,220,220,220,220,220,220,220,220,
     +          220,220,220,220,270,270,270,270,270,270,
     +          270,220,270,270,270,270,270,270,270,270) , jrun
      end if
c
 265  if (.not. ocifp) go to 270
 200  if (omrdci) go to 270
      if (occsd.or.oqcisd) then
       if (zscftp.ne.zscf(1)) call caserr2(
     + 'ccsd only available for closed shell scf wavefunctions')
      endif
c
 210  ncore = norbc
 220  if(odirpa) go to 270
      call filchk(n4file,n4blk,n4last,n4tape,m2)
      call filchk(n6file,n6blk,n6last,n6tape,m3)
      acc1 = 10.0d0**(-iacc4)
      write (iwr,6720) acc1
c
      call setbfa
      write (iwr,6730) zsort , nsz
c
 230  if (jrun.eq.3.or.ocifp) then
         if (omrdci) go to 270
         if (.not.ofulci .and. .not. occsd.and. .not. oqcisd
     +                   .and. .not. opark) then
            write (iwr,6730) zpsort , ns2
            call setbfb
            call filchk(n5file,n5blk,n5last,n5tape,m9)
         end if
      end if
      if (jrun.eq.9) go to 330
c
      if (.not.(oactrn)) then
c....    initialise active taking acount of frozen/harmonic
c...     reset norbt to num for some  runtypes (harmonic)
         if ((jrun.ge.24.and.jrun.le.27).or.jrun.eq.31.or.
     1        jrun.eq.33.or.jrun.eq.34.or.jrun.eq.22) norbt = num
c
         norbt = norbt - ncore
         do i = 1 , norbt
            norbta(i) = i + ncore
         enddo
      else
c...     adjust active (take out cores)
         nnorbc = 0
 245     do i=1,norbt
            do j=1,norbc
               if (norbta(i).eq.norbca(j)) then
                  nnorbc = nnorbc + 1
                  norbt = norbt - 1
                  do k=i,norbt
                     norbta(k) = norbta(k+1)
                  end do
                  go to 245
               end if
            end do
         end do
         if (nnorbc.gt.0) then
           write(iwr,246) nnorbc,norbt
 246       format(' *** in eliminating',i6,' cores # actives reduced to'
     +            ,i6,' ***')
           if (nnorbc.ne.norbc) call caserr2('some but not all cores')
         end if
      end if
 250  if(odci) then
       if(ninnex.gt.norbt)
     + call caserr2(
     + 'inconsistent orbital specification in direct-ci')
      endif
      write (iwr,6740) ihamc
      nw1 = 117
      len = lensec(mach(14)) + lensec(mach(16)) + lensec(nw1)
      call secput(isect(470),1005,len,isecbl)
      call wrt3i(norbt,mach(14)*nav,isecbl,idaf)
      call wrt3is(norb2,mach(16)*nav,idaf)
      call wrtcs(ztda,nw1,idaf)
      go to 270
 260  if (jrun.eq.2 .or. jrun.eq.3 .or.ocifp) go to 230
 270  if (jrun.eq.10) then
         write (iwr,6750)
         if (idma.ne.0) then
            write (iwr,6760)
         end if
         if (iplot.ne.0) then
            write (iwr,6770)
         end if
         if (iprop.ne.0) then
            write (iwr,6780)
            nword = 400 + 204/nav
            call wrt3(vlist,nword,ibl7la,num8)
            nword = 831
            call wrtcs(ztype,nword,num8)
            ibl7la = iposun(num8)
         end if
         if (local.ne.0) then
            if(local.eq.1) write (iwr,6790)
            if(local.eq.2) write (iwr,6795)
            if (nacta.le.0) go to 310
            if (nacta.le.num) go to 290
 280        call caserr2(
     +      'error detected in list of orbitals to be localised')
 290        do 300 i = 1 , nacta
               if (louta(i).gt.num) go to 280
 300        continue
            write (iwr,6800) zatext
            write (iwr,6810) (louta(i),i=1,nacta)
 310        if (nactb.gt.0) then
               if (nactb.gt.num) go to 280
               do 320 i = 1 , nactb
                  if (loutb(i).gt.num) go to 280
 320           continue
               write (iwr,6800) zbtext
               write (iwr,6810) (loutb(i),i=1,nactb)
            end if
            nword = 2*(maxorb+1)/nav
            call wrt3i(nacta,nword*nav,ibl7la,num8)
            ibl7la = iposun(num8)
         end if
         if (imull.gt.0) then
            write (iwr,6820)
            nword = (maxorb+6)
            call wrt3i(omulla,nword,ibl7la,num8)
            call wrtcs(zmul,num,num8)
            ibl7la = iposun(num8)
         end if
         if(opotf(0))write(iwr,6825)
         if(inbo .gt. 0)write(iwr,6826)
         if(local.eq. 0)write(iwr,6827)
      end if
c
c     ----- save /scfwfn/ on direct access -----
c
      call secput(isect(504),m23,m4,iblkw)
      call wrt3(cicoef,mach(5),iblkw,idaf)
      iblkw = iblkw + 2
      call wrt3(cicoef,mach(5),iblkw,idaf)
c
c     ----- fix integral block count
c
330   nintmx = num2e
      if (omp2 .or. mp3 .or. oci .or. omcscf .or. ocifp. or.
     +    zscftp.eq.zscf(5). or. zscftp.eq.zscf(6).
     +    or. jrun.eq.9 .or. oso) nopk = 1
      opk = nopk.ne.1
      opandk = (opk .and. zscftp.eq.zscf(2)) .or. 
     +         (opk .and. na.ne.nb).or.
     +         (opk .and. zscftp.eq.zscf(3)) .or.
     +         (opk .and. zscftp.eq.zscf(4)) .or.
     +         (opk .and. nopk.eq.-1)
      if (opandk .and. nopk.eq.0) nopk = -1
c
      if (.not.opass5) then
         if (na.ne.nb .and. zscftp.eq.zscf(1)) call caserr2(
     +       'closed-shell scf requested for open-shell system')
      end if
      if (opandk) nintmx = num2ejk
c
      if(nat.eq.1.or.zgroup.eq.'c1'.or.npair.ne.0.or.
     +   zscftp.eq.zscf(5). or. zscftp.eq.zscf(6).or.
     +   zruntp.eq.zrunt(8)) 
     +      otsym = .false.
c
      if(nat.gt.5. and. oschw) then
       oschw = .true.
      else
       oschw = .false.
      end if
      if (nprint.ne.0) write (iwr,6830) nprint
      if (obasis) then
         if (.not.odnew) then
            icut0 = iabs(icut)
            write (iwr,6840) itol , icut0 , ist , jst , kst , lst
            if (intg76.eq.0) write (iwr,6850)
            if (nopk.eq.0) write (iwr,6860)
            if (nopk.eq.-1) write (iwr,6870)
            if (nopk.eq.1) write (iwr,6880)
            if (oschw)  write(iwr,7100)
            if (omem(3)) write(iwr,7110)
         end if
         if (offsym) write (iwr,6890)
         if (.not.odnew) write (iwr,6900)
         if (orege) then
            opass1 = .false.
            opass2 = .false.
            opass3 = .false.
            opass4 = .false.
           if (orgall) then
              opass5 = .false.
              opass6 = .false.
              opass8 = .false.
              opass9 = .false.
           end if
         end if
         if (opass3) opass2 = .true.
         if (opass1 .or. opass2 .or. opass3 .or. opass4 .or. opass5 .or.
     +       opass6) write (iwr,6910)
         if (opass1) write (iwr,'(1x,a)') '1e-integral evaluation'
         if (opass2) write (iwr,'(1x,a)') '2e-integral evaluation'
         if (opass3) write (iwr,'(1x,a)') 'old mainfile to be used'
         if (opass4) write (iwr,'(1x,a)') 
     *               'symmetry adapted integral evaluation'
         if (opass5) write (iwr,'(1x,a)') 'scf processing'
         if (opass6) write (iwr,'(1x,a)') 'integral transformation'
         if (opass8) write (iwr,'(1x,a)') 'gradient evaluation'
         if (opass9) write (iwr,'(1x,a)') 'chf processing'
         if (opass10) write (iwr,'(1x,a)') 'wavefunction analysis'
         if (opass11) write (iwr,'(1x,a)') '4index transformation'
         if (opass1 .or. opass2 .or. opass3 .or. opass4 .or. opass5 .or.
     +       opass6 .or. opass8 .or. opass9 .or. opass10 .or. opass11 ) 
     +       write (iwr,6980)
      end if
      write (iwr,6990)
      if (.not.(odist)) then
c     suppress distance checking until optimisation itself
         otmp = ochkdst
         ochkdst = .false.
         call intr(core)
         ochkdst = otmp
      end if
c
c     add checking and dumpfile section for ecp calculations
c
 340  if(lpseud.eq.1) then
       call inpseu2(core)
      endif
      if (.not.odnew .or. jrun.ne.10) then
         if (intg76.gt.0) then
c
c     calculate error function table for rotated-axis
c     integrals
c
            call sectst(isect(502),i)
            if (i.eq.0) call gamgen
         end if
c
         otest1 = jrun.ge.4 .and. jrun.le.8
         otest2 = jrun.ge.23 .and. jrun.le.26
         if (otest1 .or. otest2) then
c
c     calculate error function table for rotated-axis
c     derivative integrals
c
            if (intg76.gt.1) then
               call sectst(isect(503),i)
               if (i.eq.0) call gamgn2
            end if
         end if
      end if
        if (oreact) then
         if (field .eq. ' ' .and. .not. obeen2) then
c
c -- for Dipole Preserving Charge analysis
c    in the absence of a DRF calculation
c
           isingl=lenwrd()
           nbits =maxbit/isingl
           vector=.false.
           idafh=110
           ndar=999
           call daopen(idafh,ioda,navh,ndar)
           scftp = zscftp  
         endif
        endif
      call timana(1)
      call timit(0)
      return
 6010 format (/1x,'accuracy factors for symmetry operations',2(2x,2e9.2)
     +        )
 6020 format (/1x,'input x,y,z values'/1x,'for chf calculation',4x,
     +        3(1x,f11.7))
 6030 format (' **** symmetry adaption suppressed')
 6040 format (/' time safety factor      ',f8.2,
     +        ' seconds'/' integral dumping period ',f8.2,' seconds'/)
 6050 format (/' time safety factor      ',f8.2,' seconds'/)
 6060 format (/1x,
     +        '**************************************************** '/1x
     +        ,'*            JOB OPTIONS in EFFECT                 * '/1
     +        x,'**************************************************** '/
     +        1x,'* RUN TYPE                                ',a8,' *')
 6070 format (1x,'* SCF TYPE                                ',a8,
     +        ' *'/1x,'* Molecular orbital starting point        ',a8,
     +        ' *')
 6080 format (1x,'* Scale factor for ',a8,' routines ',f8.2,'      *')
 6090 format (1x,'* Restore ',a1,'-vectors from section            ',i6,
     +        ' *')
 6100 format (1x,'* Route   ',a1,'-vectors to section               ',
     +        i5,' *')
 6110 format (1x,'* Route canonicalised orbitals to section      ',i3,
     +        ' *')
 6071 format (1x,'* Use Spherical Harmonic Gaussians                 *')
 6072 format (1x,'* Compute default one-electron properties          *')
 6120 format (1x,'* Punch control                           ',i8,' *'/)
 6130 format (1x,'**************************************************** '
     +        /1x,
     +        '* Direct-mp2 *** compute integrals each iteration ** ')
 6140 format (1x,'**************************************************** '
     +        /1x,
     +        '* Direct-scf *** compute integrals each iteration ** ')
 6150 format (1x,'**************************************************** '
     +        /1x,
     +        '* Direct-scf *** compute integrals each iteration ** '/1x
     +        ,'**************************************************** '/1
     +        x,'* prefactor tolerance ',e10.2,' to iter ',i4,
     +        '     **'/1x,'*                     ',e10.2,' to iter ',
     +        i4,'     **'/1x,'*                     ',e10.2,
     +        ' to iter ',i4,'     **'/1x,
     +        '* delta density used when delta <=  ',e10.2,'    **'/1x,
     +        '* delta factor for tolerance        ',e10.2,'    **'/1x,
     +        '* integral magnitudes analysed (t/f)',l3,'           **')
 6160 format (/1x,
     +        '* Semi-empirical mopac option requested            * '/1x
     +        ,'**************************************************** '/)
 6170 format (/1x,
     +        '* Molecular mechanics option requested             * '/1x
     +        ,'**************************************************** '/)
 6180 format (/1x,
     +        '* Molecular dynamics option requested              * '/1x
     +        ,'**************************************************** '/)
 6190 format (1x,'**************************************************** '
     +        /)
 6200 format (1x,
     +        '****************************************************** '/
     +        1x,
     +        '*        SCF CONVERGENCE CONTROLS in EFFECT          * '/
     +        1x,
     +        '****************************************************** '/
     +        1x,'* Wavefunction convergence                    1.0e-',
     +        i1,' *')
 6210 format (1x,'* Convergence index                     ',i12,' *')
 6211 format (1x,
     +        '* Flag SCF convergence after',i4,
     +                                        ' iterations          *')
 6215 format (1x,
     +        '* Symmetrize Kohn-Sham matrix                        *')
 6220 format (1x,'* Damping cut-off tolerance         ',f14.2,'   *')
 6230 format (1x,
     +        '* Configuration locking requested                    *')
 6240 format (1x,
     +        '* Extrapolation                                      *')
 6320 format (1x,
     +        '* Damping (davidson)                                 *')
 6330 format (1x,
     +        '* Restrict occupied orbital mixing                   *')
 6340 format (1x,
     +        '* No damping,level-shifting, extrapolation or        *'/1
     +        x,'* Occupied m.o. mixing restriction requested         *'
     +        )
 6350 format (1x,'* Commence diis treatment at threshold of ',f10.5,
     +        ' *'/1x,'* Finish diis when residuum less than ',e14.6,
     +        ' *')
 6351 format (1x,
     +        '* Commence diis treatment at dynamic threshold       *'/
     +     1x,'* Finish diis when residuum less than ',e14.6,' *')
 6360 format (1x,'* Diis control information to ',a4,' at block ',i9,
     +        '*')
 6370 format (1x,
     +        '* diis output requested                              *')
 6376 format (1x,
     +        '* Fermi-Dirac smearing                               *'/
     +     1x,'* start smearing at    ',e10.2,' (',f9.2,' eV)     *'/
     +     1x,'* finish smearing at   ',e10.2,' (',f9.2,' eV)     *'/
     +     1x,'* scale factor         ',e10.2,'                    *')
 6380 format (1x,
     +        '****************************************************** ')
 6381 format (1x,
     +        '* perform UHF annihilation analysis                  *')
 6390 format (/1x,
     +        '**************************************************** '/1x
     +        ,'*            PRINT OPTIONS in EFFECT               * '/1
     +        x,'**************************************************** ')
 6400 format (1x,'* complete output of mulliken analysis             * '
     +        )
 6410 format (1x,'* complete output of lowdin analysis               * '
     +        )
 6420 format (1x,'* complete scf print throughout optimization       * '
     +        )
 6430 format (1x,'* output of initial vectors                        * '
     +        )
 6440 format (1x,'* output of mos in symmetry adapted basis          * '
     +        )
 6450 format (1x,'* intermediate optimization output                 * '
     +        )
 6460 format (1x,'* intermediate symmetry output                     * '
     +        )
 6470 format (1x,'* output diagonal fock elements                    * '
     +        )
 6480 format (1x,'* full vector output on completion of optimisation * '
     +        )
 6490 format (1x,'* output of property integrals (ao basis)          * '
     +        )
 6500 format (1x,'* output of property integrals (mo basis)          * '
     +        )
 6510 format (1x,'* output of dumpfile summary                       * '
     +        )
 6520 format (1x,'* output of distance matrix                        * '
     +        )
 6530 format (1x,'* intermediate direct-scf output                   * '
     +        )
 6540 format (1x,'* intermediate integral output                     * '
     +        )
 6550 format (1x,'* suppress printing of distance matrix             * '
     +        )
 6560 format (1x,'* suppress printing of basis functions             * '
     +        )
 6570 format (1x,'* suppress printing of sabf specification          * '
     +        )
 6580 format (1x,'* suppress printing of optimization history file   * '
     +        )
 6590 format (1x,'* suppress printing of non-canonicalised vectors   * '
     +        )
 6600 format (1x,'* suppress printing of wave function analysis      * '
     +        )
 6610 format (1x,'* suppress printing of force constant matrix       * '
     +        )
 6615 format (1x,'* suppress printing of bq coordinates and forces   * '
     +        )
 6620 format (/1x,
     +     '********************************************************** '
     +     /1x,
     +     '*          optimisation options in effect                * '
     +     /1x,
     +     '********************************************************** '
     +     )
 6630 format (1x,
     +      '* optimization procedure invoked                         *'
     +      )
 6632 format (1x,
     +      '* bond distance checking requested                       *'
     +      )
 6640 format (1x,'* restore hessian from dumpfile on ',a4,' at block',
     +        i8,' *')
 6650 format (1x,
     +     '* force constant matrix calculation                      * '
     +     /1x,'* use ',i1,'-point differencing formula',
     +     '                       * '/1x,
     +     '* step size in differencing formula ',f8.3,13x,'* ')
 6660 format (1x,'* convergence threshold in optimization run ',f8.3,4x,
     +        ' *')
 6670 format (1x,
     +     '********************************************************** '
     +     )
 6680 format (/' the following orbitals included in active list'/
     +        ' =============================================='/5
     +        (' function e/i label',4x)/1x,110('-'))
 6690 format (5(i9,i11,4x))
 6700 format (/' no functions specified in frozen core list')
 6710 format (/' the following mos included in frozen core list'/
     +        ' ----------------------------------------------'/5
     +        (' function e/i label',4x))
 6720 format (/1x,'4-index accuracy factor',e12.1)
 6730 format (' blocking factor for ',a5,' file = ',i3)
 6740 format (/1x,'route transformed core hamiltonian to section',i4)
 6750 format (/1x,
     +     '********************************************************** '
     +     /1x,
     +     '*              analysis options in effect                * '
     +     /1x,
     +     '********************************************************** '
     +     )
 6760 format (1x,
     +   '* perform distributed multipole analysis                 *')
 6770 format (1x,
     +   '* perform graphical analysis                             *')
 6780 format (1x,
     +   '* evaluate one-electron properties                       *')
 6790 format (/1x,
     +   '* localize molecular orbitals by boys algorithm          *')
 6795 format (/1x,
     +   '* localize molecular orbitals by Pipek-Mezey algorithm   *')
 6800 format (/1x,'* list of active  -',a1,'-  mos')
 6810 format (/20i4)
 6820 format (1x,
     +   '* perform extended mulliken analysis                     *')
 6825 format (1x,
     +   '* evaluate potential derived charges                     *')
 6826 format (1x,
     +   '* perform natural bond orbital analysis                  *')
 6827 format (1x,
     +   '********************************************************** '
     +     )
 6830 format (/' output print option                       ',i5/)
 6840 format (/1x,
     +        '**************************************************** '/1x
     +        ,'* 2-electron integral options                      * '/1
     +        x,'**************************************************** '/
     +        1x,'* prefactor tolerance for integrals        1.0e-',i2,
     +        ' *'/1x,'* integral cutoff                          1.0e-'
     +        ,i2,' *'/1x,'* starting shells                 ',4i4,' *')
 6850 format (1x,'* high accuracy integral generator requested       *')
 6860 format (1x,'* generate integrals in p-supermatrix form         *')
 6870 format (1x,'* generate integrals in p-k supermatrix form       *')
 6880 format (1x,'* generate integrals in non-supermatrix form       *')
 6890 format (1x,'* suppress symmetry in 2e-integral generator       *')
 7100 format (1x,'* apply schwarz inequality screening               *')
 7110 format (1x,'* store 2e-integrals in memory                     *')
 6900 format (1x,'****************************************************'/
     +        )
 6910 format (/1x,36('=')/' bypass'/1x,36('='))
 6980 format (1x,36('=')//)
 6990 format (/1x,104('-'))
 7030 format(  5x,'- ',i3,' core electrons not considered')
 7020 format(  5x,'- symmetry equivalences among atoms',/)
 7010 format(/,5x,'Correlation energy associated with HF density',/,
     +          5x,'evaluated using',/,
     +          5x,'- ',a35,/,
     +          5x,'- ','points in angular integration',i4)
 7060 format(/1x,'values entered for scrf calculation'/
     +        1x,'solvent dielectric',7x,f10.4/
     +        1x,'cavity dimensions x,y,z',1x,
     + 3(1x,f10.4) )
 7065 format(/1x,'bq atoms will be treated as crystal field model')
 7070 format(/1x,54('*')/
     + 1x,'**** polarised surface reaction field calculation ****'/
     + 1x,'****     solvent dielectric =',f8.3,'             ****'/
     + 1x,'**** surface point charges spaced at ',f8.3,'degs ****'/
     + 1x,'**** no. of user defined centres for surface = ',i3,' ***')
 7080 format(1x,'******   radius of sphere ',i3,' = ',f7.4,
     + ' angstrom ******')
 7090 format(
     + ' **** using ',i3,' iterations damped by factor ',f5.4,'  ****'/
     + ' **** point charges dumped to section       ',i5,  '  ****'/
     +   1x,54('*')/)
      end
      subroutine levsin(olevel,odft)
c
c    now check level shifters and set-up defaults
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      common/blbpar/accblb,clevs,timblb,nitblb(2)
c
      dimension zscf(6)
      data zscf/'rhf','uhf','gvb','grhf','casscf','mcscf'/
      data done, pt3, pt5, pt05, pt01 / 
     +     1.0d0, 0.3d0, 0.5d0, 0.05d0, 0.01d0 /

        write (iwr,6250)
         if (olevel) then
c  level shifters have been input - do not override
            write (iwr,6260)
         else
c
c  set up default values
c  now modify original settings to reflect presence of 1st row
c  transition metal, and subsequent reqirement for higher values
c
         otm = .false.
         do i = 1, nat
         nucz = isubst(zaname(i))
         if (nucz.ge.21.and.nucz.le.30) otm = .true.
         enddo
c
          if (zscftp.eq.zscf(1)) then
c  retain higher level shifter throughout for dft calculations
             gapa1 = done
             if(.not.odft) then
                ibrk = 5
               gapa2 = pt3
             endif
          else if (zscftp.eq.zscf(2)) then
c
c ----- default level shifters for uhf
c
             gapa1 = done
             gapb1 = done
             if(.not.odft) then
c  retain higher level shifter throughout for dft calculations
               ibrk = 5
               gapa2 = pt3
               gapb2 = pt3
             endif
          else if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4)) then
c
c ----- default level shifters for grhf/gvb
c
             gapa1 = pt05
             gapb1 = done
             ibrk = 5
             gapa2 = pt01
             gapb2 = pt5
c
c ----- casscf/ mcscf
c
          else
             gapa1 = clevs
          end if
c 
          if (otm) then
           gapa1 = gapa1 * 3.0d0
           gapb1 = gapb1 * 3.0d0
           gapa2 = gapa1
           gapb2 = gapb1
          endif
c
         end if
c
c
         if (zscftp.eq.zscf(5) .or. zscftp.eq.zscf(6)) then
            write (iwr,6270) gapa1
            clevs = gapa1
         else
            if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4))
     +          write (iwr,6280)
            write (iwr,6290) gapa1 , ibrk , gapa2
            if (zscftp.ne.zscf(1)) then
               if (zscftp.eq.zscf(2)) write (iwr,6310)
               if (zscftp.eq.zscf(3) .or. zscftp.eq.zscf(4))
     +             write (iwr,6300)
               write (iwr,6290) gapb1 , ibrk , gapb2
            end if
         end if
c
 6250 format (1x,
     +        '* Level shifting                                     *')
 6260 format (1x,
     +        '* Input level shifters specified                     *')
 6270 format (1x,'* Super-ci level shifter = ',f8.2,
     +        '                  *')
 6280 format (1x,
     +        '* Occupied-occupied orbital mixing                   *')
 6290 format (1x,'* Level shifter = ',f8.3,' to cycle',i3,' then ',f8.3,
     +        ' *')

 6300 format (1x,
     +        '* Occupied-virtual orbital mixing                    *')
 6310 format (1x,
     +        '* Beta hamiltonian :-                                *'/
     +     1x,'* Occupied-occupied orbital mixing                   *')

         return
         end
c
c-----------------------------------------------------------------------
c
      subroutine start3(jrun,oreact)
c
c     Now check the options requested in the input file against the
c     implemented capabilities.
c
c     implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
c     implicit character *8 (z),character *1 (x)
c     implicit character *4 (y)
      implicit none
c
c     In the type declarations below please respect the traditional
c     implicit declaration conventions...
c
      integer jrun
      logical oreact
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      real*8  vg, gxp, gyp, gzp, d2, valmo, occa, occb
      real*8  rij, timat, denat, encoat, eltot
      integer idenfu,nangpr, ncorxx, ndeg, nomo
      logical odenfu, osymm, osc
      common /dnfnw/ odenfu, osymm, osc, idenfu, nangpr, ncorxx,
     +               vg(200),gxp(200),gyp(200),gzp(200),d2(200),
     +               valmo(200),occa(200),occb(200),
     +               rij(1275),timat(50),denat(50),encoat(50),ndeg(50),
     +               eltot,nomo
c
c
      real*8  dx0, dxm, dele, drms
      integer ntpr, ntrun, imcyc, ncyc, ntmin, nspacb
      common /modmin/ ntpr,ntrun,imcyc,ncyc,ntmin,nspacb,
     +                dx0,dxm,dele,drms
c
c
      integer  ispchk
      external ispchk
      integer  ibasis
c
c     finally, check that the calculation requested is possible
c     given the basis set functions (s,p,d,f,g), and store
c     logical pointers in restrl specifying the latter.
c     
c     first, ensure that section 491 is present
c
      if (omech) then
c
        if (ntmin.eq.1  .or. ntmin.eq.2 .or. ntmin.eq.3) return
c
      endif
c
      ibasis = ispchk()
c
c     runtypes excluded by f-functions in basis
c
c
c disabled to get potentials
c
      if (ofbas) then
        if (iprop.gt.0) call caserr2(
     +    'property option not available with f-functions in the basis')
c       if (iplot.gt.0 ) then
c         additional check on potentials and density required here
c         now conducted in analis
c       endif
        if (odenfu) call caserr2(
     +    'old dft option not available with g-functions')
      endif
c
c     runtypes excluded by g-functions in basis
c
      if (ogbas) then
        if (oreact) call caserr2(
     +    'drf option not available with g-functions in the basis')
        if (odenfu) call caserr2(
     +    'old dft option not available with g-functions')

        if (odscf) then
          if (odnew) call caserr2(
     +      'direct calculations not possible with g-functions')
        end if
c
c       now screen specified runtype
c
        go to (60,60,60,60,60,60,60,60,60,70,
     +         60,60,60,60,60,60,60,60,60,60,
     +         80,80,80,80,80,80,80,80,60,80,
     +         80,80,80,80,60,60,60,60,60,60,
     +         60,90,50,60,60,60,60,60,60,60) , jrun
c
 70     if (iprop.gt.0 .or. idma.gt.0) then
          call caserr2(
     +    'analysis option not available with g-functions in the basis')
        endif
        go to  60
c
 80     call caserr2(
     +    'runtype option not available with g-functions in the basis')
        go to 60
c
 50     call caserr2(
     +    'drf option not available with g-functions in the basis')
        go to 60
c
 90     if (odirpa) call caserr2(
     +    'direct-RPA not available with g-functions in the basis')
c
 60     continue
c
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine test(q)
c
c     routine provided for user-supplied purposes
c     invoked by runtype test; includes core array (q) 
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*10 charwall
      dimension q(*)
      write(iwr,10) cpulft(1) ,charwall()
10    format(/1x,'end of test at ',f8.2,' seconds',a10,' wall')
      return
      end
      subroutine tranci(core)
c
c   1. compute 1&2 electron integrals
c   2. carryout scf/casscf/mcscf calculation
c   3. perform integral transformation
c   4. perform ci calculation
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
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
      real*8 enci, ec, eci
      integer nstat
      common /enrgci/ enci,ec,eci(10),nstat
c
      integer istate
      common /enrgc2/ istate(10)
c
c
      integer irfpas, maxitrf, istci
      common /rfci1/ irfpas,maxitrf,istci
      real*8 rfcvg
      common /rfci2/ rfcvg
      logical extra
      common /rfci3/ extra
c
      real*8 fzero, f1, ffp, alphf
      integer k, nvarf, modef, nstep, indfx, lamba
      logical ocnvrg, ocifp
      common /fpinfo/ fzero,f1(4),ffp,alphf,k,nvarf,ocnvrg,
     +                modef,nstep,indfx,lamba,ocifp
      logical oarf, oarfcvg
      character*8 fieldold
cdrf
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
      dimension core(*),zci(3),zcas(3)
c
      character*8 fnm
      character*6 snm
      data fnm,snm/"master.m","tranci"/
c
      data oyes,ono / .true.,.false./
      data zci,zcas  / 'ci','gf','tda', 
     + 'casscf','mcscf','vb'/
      data m10,m16/10,16/
c
      inull = igmem_null()
      oarf = field(5:) .eq. 'scf '
      fieldold = field
      enadd = 0.d00
      if (oarf) ocifp = .true.
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
     &         .or. zscftp.eq.zcas(3)
      nt = 1
      iofsym = 1
      nopk = 1
      enrgy = 0.0d0
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
            if (tim.ge.timlim) go to 80
         end if
      end if
      ionsv = 1
c
      if (opass2) then
         if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if (irest.le.1) then
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "prefac",IGMEM_DEBUG)
            call rdmake(core(iprefa))
            call jandk(zscftp,core,core(inull),core(inull),core(inull),
     +                 core(inull),core(inull),core(iprefa),core(inull))
            call gmem_free_inf(iprefa,fnm,snm,"prefac")
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6010)
               go to 80
            end if
         end if
      end if
c
      narfpts = 0
 888  continue
      oarfcvg = .false.
      if (oarf) then
         if (narfpts.gt.0) then
            call arfupd(core,enadd,fieldold)
            field = ' '
         end if
         write(iwr,901) narfpts+1
901   format(/,1x,78('='),/,10x,'Average Reaction field MCSCF/CI ',
     1       'Calculation cycle',i3,/,1x,78('='),/,1x)
      end if
c
      if (omccas) then
         if (.not.opass4) then
            if (irest.le.2) then
               call adapti(core,ono)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6020)
                  go to 80
               end if
            end if
         end if
      nt = ntsave
      end if
c
      if (.not.opass5) then
         if (irest.le.4) then
            call scf(core)
            if(.not.omccas) then
               call putdev(core,mouta,7,1)
c
               call secget(isect(494),m16,iblk34)
               call rdedx(potnuc,m10,iblk34,idaf)
               enuc = potnuc(1)
               ehfock = potnuc(2)
               etotal = potnuc(3)
c
               call secget(isect(13),13,iblok)
               call wrt3(enuc,lds(isect(13)),iblok,idaf)
            endif
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6030)
               go to 80
            end if
         end if
      end if
c
      if (oarf.and.narfpts.eq.0)
     1     call arfupd(core,enadd,fieldold)
c
      if (omrdci) then
         call mrdcim(core,enrgy)
      else
         if(opark) then
c consistency check for mrdci data and defaults
          call parkinc(core)
         endif
         if (.not.opass6) then
            if (irest.le.8) then
               call adapti(core,oyes)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6040)
                  go to 80
               end if
            end if
         end if
         i = locatc(zci,3,zruntp)
         if (i.eq.0) go to 80
         go to (20,30,40) , i
 20      if (ofulci) then
            call fci(core,enrgy)
         else if (occsd.or.oqcisd) then
            call titandrv(core,enrgy)
c rwah
         else if (opark) then
            call mrdcim2(core,enrgy)
         else if (odiesel) then
c Here we will create the calls to the diesel code
            call diesel
         else
            call direct(core,enrgy)
         end if
         go to 70
 30      call gf(core)
         go to 70
 40      call tda(core)
         go to 70
      end if
 70   continue
c
      if (oarf) call arfcvg(core,enrgy,enadd,oarfcvg,narfpts)
      if (fieldold .ne. ' ') then
        istci = iactst
      else 
        istci = 1
      endif
      eci(istci) = enrgy
      narfpts = narfpts + 1
      if (oarf .and. .not. oarfcvg) goto 888
      if (oarf) field = fieldold
c
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
      subroutine respon(core)
c
c   1. compute 1&2 electron integrals
c   2. carryout scf (rpa,dirrpa) or mcscf (mclr) calculation
c   3. perform integral transformation
c   4. perform response calculation
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      character *8 title,guess,conf
      common/restrz/title(12),guess,conf(37)
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
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
c
      dimension core(*),zcas(3)
c
      character*8 fnm
      character*6 snm
      data fnm,snm/"master.m","respon"/
c
      data oyes,ono / .true.,.false./
      data zcas  / 'casscf','mcscf','vb'/
      data m10,m16/10,16/
c
      inull = igmem_null()
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2) 
     &         .or. zscftp.eq.zcas(3)
      if (.not.odirpa) then
         nt = 1
         iofsym = 1
         nopk = 1
      endif
      enrgy = 0.0d0
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
            if (tim.ge.timlim) go to 80
         end if
      end if
      ionsv = 1
c
      if (opass2) then
         if(.not.ognore.and..not.odirpa) 
     &        call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if(.not.odscf) then
            if (irest.le.1) then
               iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                                  "prefac",IGMEM_DEBUG)
               call rdmake(core(iprefa))
               call jandk(zscftp,core,core(inull),core(inull),
     +                    core(inull),core(inull),core(inull),
     +                    core(iprefa),core(inull))
               call gmem_free_inf(iprefa,fnm,snm,"prefac")
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6010)
                  go to 80
               end if
            end if
         end if
      end if
c
      if (omccas) then
         if (.not.opass4) then
            if (irest.le.2) then
               call adapti(core,ono)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6020)
                  go to 80
               end if
            end if
         end if
      nt = ntsave
      end if
c
      if (.not.opass5) then
         if (irest.le.4) then
            call scf(core)
            if (odscf. and. guess.eq.'atoms') go to 1000
            if(.not.omccas) then
               call putdev(core,mouta,7,1)
c
               call secget(isect(494),m16,iblk34)
               call rdedx(potnuc,m10,iblk34,idaf)
               enuc = potnuc(1)
               ehfock = potnuc(2)
               etotal = potnuc(3)
c
               call secget(isect(13),13,iblok)
               call wrt3(enuc,lds(isect(13)),iblok,idaf)
            endif
1000        if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6030)
               go to 80
            end if
         end if
      end if
c
      if (orpa) then
         if (.not.opass6) then
            if (irest.le.8) then
               call adapti(core,oyes)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6040)
                  go to 80
               end if
            end if
         end if
      end if
c
c    linear response
c
      if (omclr) then 
c
c     allocate core (all available memory)
c
         ibase = igmem_alloc_all(lword)
         call lrmain(core(ibase),core(ibase),lword)
         call gmem_free(ibase)
      else if (orpa .or. odirpa) then
         call rpdriv(core,core)
      else
         call caserr2 ('unrecognised response type')
      end if
c
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
      subroutine citran(core)
c
c   1. compute 1&2 electron integrals
c   2. carryout scf/casscf/mcscf calculation
c   3. perform integral transformation
c   4. perform ci calculation
c
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      logical oarf, oarfcvg
cdrf
c
      character*8 field
      character*24 version
      common /drfopt/ field,version
      logical oreact
      integer intdrf
      common /drfopti/ oreact,intdrf
c
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
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
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
c
      dimension core(*),zcas(3)
c
      character*8 fnm
      character*6 snm
      data fnm,snm/"master.m","citran"/
c
      data oyes,ono / .true.,.false./
      data zcas  /  'casscf','mcscf','vb'/
      data m10,m16/10,16/
c
      enadd = 0.d00
      oarf = field(5:) .ne. ' '
      inull = igmem_null()
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
     &         .or. zscftp.eq.zcas(3)
      nt = 1
      iofsym = 1
      nopk = 1
      enrgy = 0.0d0
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
            if (tim.ge.timlim) go to 80
         end if
      end if
      ionsv = 1
c
      if (opass2) then
         if (.not.ognore)call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if (irest.le.1) then
            iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                               "prefac",IGMEM_DEBUG)
            call rdmake(core(iprefa))
            call jandk(zscftp,core,core(inull),core(inull),core(inull),
     +                 core(inull),core(inull),core(iprefa),core(inull))
            call gmem_free_inf(iprefa,fnm,snm,"prefac")
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6010)
               go to 80
            end if
         end if
      end if
c
      if (omccas) then
         if (.not.opass4) then
            if (irest.le.2) then
               call adapti(core,ono)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6020)
                  go to 80
               end if
            end if
         end if
      nt = ntsave
      end if
c
      if (.not.opass5) then
         if (irest.le.4) then
            call scf(core)
            if(.not.omccas) then
               call putdev(core,mouta,7,1)
c
               call secget(isect(494),m16,iblk34)
               call rdedx(potnuc,m10,iblk34,idaf)
               enuc = potnuc(1)
               ehfock = potnuc(2)
               etotal = potnuc(3)
c
               call secget(isect(13),13,iblok)
               call wrt3(enuc,lds(isect(13)),iblok,idaf)
            endif
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6030)
               go to 80
            end if
         end if
      end if
c
      narfpts = 0
 888  continue
      oarfcvg = .false.
      if (oarf) call arfupd(core,enadd,fieldold)
      if (omrdci) then
         call mrdcim(core,enrgy)
      else
         if(opark) then
c consistency check for mrdci data and defaults
          call parkinc(core)
         endif
         if (.not.opass6) then
            if (irest.le.8) then
               call adapti(core,oyes)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6040)
                  go to 80
               end if
            end if
         end if
         if (ofulci) then
            call fci(core,enrgy)
         else if (occsd.or.oqcisd) then
            call titandrv(core,enrgy)
         else if (opark) then
            call mrdcim2(core,enrgy)
         else
            call direct(core,enrgy)
         end if
      end if
      if (oarf) call arfcvg(core,enrgy,enadd,oarfcvg,narfpts)
      narfpts = narfpts + 1
      if (oarf .and. .not. oarfcvg) goto 888
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
c
c  Check validity to integral file and print diagnostics
c
      subroutine checkinteg(nopk,nopkr,iofsym,iofrst)
      implicit none
      integer nopk,nopkr,iofsym,iofrst
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      if(nopk.ne.nopkr)then
         write(iwr,*)' invalid format for integral file in bypass'
         write(iwr,*)' nopk:',nopk,' nopkr:',nopkr
         write(iwr,*)
     +   ' If you know it is safe, you can defeat the check'
         write(iwr,*)' with the ignore keyword on the bypass directive'
      endif
      if(iofsym.ne.iofrst)then
         write(iwr,*)' inconsistent use of skeletonisation in bypass'
         write(iwr,*)' iofsym:',iofsym,' iofrst:',iofrst
         write(iwr,*)
     +   ' If you know it is safe, you can defeat the check'
         write(iwr,*)' with the ignore keyword on the bypass directive'
      endif
      if(nopk.ne.nopkr) call caserr2(
     +     'invalid format for integral file in bypass')
      if(iofsym.ne.iofrst) call caserr2(
     +     'inconsistent use of skeletonisation in bypass')
      end
      subroutine ver_master(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/master.m,v $
     +     "/
      data revision /"$Revision: 6270 $"/
      data date /"$Date: 2012-12-03 14:27:59 +0100 (Mon, 03 Dec 2012) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
