**==aolim.f
      subroutine aolim
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
      common/junk/limlow(maxat),limsup(maxat)
c
      do 20 i = 1 , nat
         limlow(i) = 0
         limsup(i) = 0
 20   continue
      lat = katom(1)
      limlow(lat) = 1
      do 30 i = 1 , nshell
         iat = katom(i)
         if (lat.ne.iat) then
            limsup(lat) = kloc(i) - 1
            lat = iat
            limlow(lat) = kloc(i)
         end if
 30   continue
      limsup(iat) = num
c
      return
      end
**==atpop.f
      subroutine atpop(a,b,nat)
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
      common/junk/limlow(maxat),limsup(maxat)
      dimension a(*),b(*)
c
      data dzero /0.0d0/
c
      ind(i,j) = ((max(i,j)*(max(i,j)-1))/2+min(i,j))
c
      j = nat*(nat+1)/2
      call vclr(b,1,j)
      do 50 i = 1 , nat
         i1 = limlow(i)
         i2 = limsup(i)
         do 40 j = 1 , nat
            j1 = limlow(j)
            j2 = limsup(j)
            dum = dzero
            if (i1.ne.0 .and. j1.ne.0) then
               do 30 k = i1 , i2
                  do 20 l = j1 , j2
                     kl = ind(k,l)
                     dum = dum + a(kl)
 20               continue
 30            continue
            end if
            ij = ind(i,j)
            b(ij) = dum
 40      continue
 50   continue
      return
      end
**==boyloc.f
      subroutine boyloc(q)
c
c     boyloc calculates the localized molecular orbitals
c     by either :
c     1.  method of foster and boys. this implementation is qcpe 
c         program 354,  modified to run as part of gamess.
c     2.  population method of j. pipek and p. g. mezey,
c         j. chem. phys. 90, 4916 (1989),
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/miscop/nacta,nactb,moout(maxorb),mooutb(maxorb)
     * , eiga(maxorb),frocca(maxorb),eigb(maxorb),froccb(maxorb)
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
      common/junkc/zcom1(29),zcomm(19),zbuff(10)
      common/blkcore/corev(512),charge(4),cpad(6)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension q(*),ytype(2),o1e(6)
      logical once
c
      data ytype/'-a-','-b-'/
      data itvec/3/
c
c     allocate space for localisation
c
      once = .false.
c
c     restore active list
c
      nav = lenwrd()
      nw1e2 = 400 + 204/nav
      nwdma = 5 + 4*maxat + (13+maxat+nav)/nav
      nwp2 = 2*maxorb + 82 + 60/nav
      i = 1
      if (idma.gt.0) i = i + lensec(maxat) + lensec(nwdma)
      if (iplot.gt.0) i = i + iplot*(lensec(nwp2)+lensec(140))
      if (iprop.gt.0) i = i + lensec(nw1e2) + lensec(816)
      nw = 2*(maxorb+1)
      call readi(nacta,nw,i,num8)
      ibl7qa = ibl7la
      lenb = num*(num+1)/2
      l3 = num*num
      n1 = max(nacta,nactb)
      if (n1.le.0 ) then
c
c     construct default list of active mos, allowing
c     for core mos as a function of scftype
c
        odef = .true.
        call lmopsi(zscftp,mcore,mdoc,mact,mocc)
      else
        odef = .false.
      endif
      n1 = max(nacta,nactb)
      if (n1.lt.0 .or. n1.gt.newbas) go to 30
      n2 = n1*(n1+1)/2
      n3 = n1*n1
c
c     allocate space. 
c
      i10 = 1
c
      if(local.eq.1) then
c
c     the arrays used in boyorb for each address are
c     i10  cc,cl
c     i20  b,dmao,iord
c     i30  iir
c     i40  ia
c     i50  rij
c     i60  ri
c     i70  moout
c     i80  qpix
c     i90  qpjx
c     i100 norot
c
1      i20 = i10 + l3
       i30 = i20 + max(newbas*n1,lenb)
       i50 = i30 + newbas
       i60 = i50 + max(3*n2,lenb)
       i70 = i60 + newbas + newbas + newbas
       i80 = i70 + newbas
       i90 = i80 + newbas
       i100= i90 + newbas
       if (.not.once) then
          once = .true.
          need = i100 + n3
          i10 = igmem_alloc_inf(need,'anala.m','boyloc',
     1                         'i10',IGMEM_DEBUG )
          write (iwr,6010)
          go to 1
       end if
c
      else if (local.eq.2) then
c
2      ltran = i10  + newbas*newbas
       lvout = ltran + n1*n1
       ltri  = lvout + newbas*newbas
       lq    = ltri  + lenb
       lwrk  = lq    + l3
       liord = lwrk  + newbas
       liir  = liord + n1
       lia   = liir  + n1
       lrij  = lia   + newbas
       lqpix  = lrij  + n2*nat
       lqpjx = lqpix + newbas
       lnorot= lqpjx + newbas
       if (.not.once) then
          once = .true.
          need  = lnorot+ n3
          i10 = igmem_alloc_inf(need,'anala.m','boyloc',
     1                         'i10',IGMEM_DEBUG )
          write (iwr,6040)
          go to 2
       end if
       else
      call caserr2('invalid localisation index detected')
      endif
c     if (lwordp.lt.need) then
c        write (iwr,6020) lwordp , need
c        call caserr2(
c    +    'insufficient memory for localised orbital analysis')
c     end if
c
c     read in the full vector matrices and write them out on ed7
c
      if (nacta.eq.0) go to 60
      if (mina.gt.0) go to 40
 20   call caserr2('invalid section specified for input vectors')
 30   call caserr2(
     +     'error detected in list of orbitals to be localised')
 40   call secget(mina,itvec,iblkv)
      call getqp(zcomm,zbuff,eiga,frocca,nbas,newb,ncola,ieiga,idiffa,
     +           maxorb,iblkv)
      do 50 i = 1 , nacta
         if (moout(i).gt.ncola) go to 30
 50   continue
      write (iwr,6030) ytype(1) , mina , ibl3d , yed(ied3) ,
     +                (zcomm(7-i),i=1,6) , zbuff
      call rdedx(q(i10),l3,iblkv,idaf)
c
c      conduct analysis in non-adapted basis set
c
      if (.not.otran) call tdown(q(i10),ilifq,q(i10),ilifq,ncola)
      call wrt3(q(i10),l3,ibl7qa,num8)
c
 60   if (nactb.ne.0) then
c
c      restore b-set vectors
c
         if (minb.le.0) go to 20
         call secget(minb,itvec,iblkv)
         call getqp(zcomm,zbuff,eigb,froccb,nbas,newb,ncolb,ieigb,
     +              idiffb,maxorb,iblkv)
         write (iwr,6030) ytype(2) , minb , ibl3d , yed(ied3) ,
     +                   (zcomm(7-i),i=1,6) , zbuff
         do 70 i = 1 , nactb
            if (mooutb(i).gt.ncolb) go to 20
 70      continue
         ibl7qb = iposun(num8)
         call rdedx(q(i10),l3,iblkv,idaf)
         if (.not.otran) call tdown(q(i10),ilifq,q(i10),ilifq,ncolb)
         call wrt3(q(i10),l3,ibl7qb,num8)
      end if
c
c     reset iftran regardless of -a- and -b- vectors
c     this now involves resetting contents of /tran/
c     to reflect non-adapted ctrans list
c
      otran = .false.
      loop = 0
      do i = 1 , num
         ilifc(i) = loop
         ntran(i) = 1
         ctran(i) = 1.0d0
         itran(i) = i
         loop = loop + 1
      enddo
      do i = 1,6
         o1e(i) = .false.
      enddo
c
      if(local.eq.1) then
c
c ----- now restore dipole ints for Boys and save on scratch
c
       do i =4,6
         o1e(i) = .true.
       enddo
       call getmat(q(i10),q(i10),q(i10),q(i10),q(i20),q(i50),
     +             charge,newbas,o1e,ionsec)
       ibl7x = iposun(num8)
       call wrt3(q(i10),lenb,ibl7x,num8)
       ibl7y = iposun(num8)
       call wrt3(q(i20),lenb,ibl7y,num8)
       ibl7z = iposun(num8)
       call wrt3(q(i50),lenb,ibl7z,num8)
       ibl7la = iposun(num8)
c
c     call the Foster-Boys localization module
c
       call boyorb(zscftp,q(i10),q(i10),q(i20),q(i20),q(i20),q(i30),
     +             q(i50),q(i60),q(i80),q(i90),q(i100),
     +             newbas,lenb,l3,n1,n2,odef,nprint,q)
c
      else if(local.eq.2) then 
c
c ----- now restore overlap ints for localisation and save on scratch
c
       o1e(1) = .true.
       call getmat(q(i10),q(i10),q(i10),q(i10),q(i10),q(i10),
     +             charge,newbas,o1e,ionsec)
       ibl7x = iposun(num8)
       call wrt3(q(i10),lenb,ibl7x,num8)
       ibl7la = iposun(num8)
       lsao=ltri
       lmap=lia
c
c     call the Population localization module
c
       call poporb(zscftp,q(i10),q(ltran),q(lvout),
     +             q(liord),q(liir),q(lmap),q(lrij),
     +             q(lqpix),q(lqpjx),q(lnorot),q(lsao),
     +             newbas,lenb,l3,n1,n2,nat,odef,nprint)
      else
       call caserr2('invalid localisation index detected')
      endif
c
      call gmem_free_inf(i10,'anala.m','boyloc','i10')
c
      return
 6010 format (/40x,31('=')/
     +         40x,'Foster-Boys Localisation Module'/
     +         40x,'Rev.Mod.Phys. 32, p300 (1960)'/
     +         40x,31('=')/)
 6040 format (/40x,34('=')/
     +         40x,'MOS Localised by Population method'/
     *         40x,'J. Chem. Phys. 90, 4916 (1989)'/
     +         40x,34('=')/)
 6020 format (/10x,'core store analysis'/10x,19('-')/
     +         10x,'main core available = ',i8,' words'/
     +         10x,'main core required  = ',i8,' words')
 6030 format (//1x,a3,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',
     +        a4//' header block information : '/
     +        ' vectors created under account ',a8/1x,a7,
     +        'vectors created by ',a8,' program at ',a8,' on ',a8,
     +        ' in the job ',a8/' with the title : ',10a8)
      end
**==boyorb.f
      subroutine boyorb(zzzzzz,cc,cl,b,dmao,iord,iir,rij,ri,
     + qpix,qpjx,norot,l1,l2,l3,n1,n2,odef,nprint,q)
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
      dimension cc(*),cl(l1,n1),dmao(l2),iir(l1),qpix(l1),
     +qpjx(l1),b(l1,n1),iord(l1),ri(l1,3),rij(n2,3),ytype(2)
      dimension norot(n1,n1),ibl(3)
      dimension q(*)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/miscop/nacta,nactb,moout(maxorb),mooutb(maxorb)
     * , eiga(maxorb),frocca(maxorb),eigb(maxorb),froccb(maxorb)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      character*10 charwall
      parameter (maxboys=1500)
c...   ** does not cost a thing (jvl 2006) **
      data ytype/'-a-','-b-'/
      data m1/1/
c
c ----
c
      cpu = cpulft(1)
      write (iwr,6130) cpu,charwall()
      nsav = mouta
      mi = 0
      norb = nacta
      iblock = ibl7qa
      ibl(1) = ibl7x
      ibl(2) = ibl7y
      ibl(3) = ibl7z
 20   mi = mi + 1
      if (norb.eq.0) go to 330
      call gnorot(zzzzzz,mi,norot,n1,norb,odef)
      nredo = 0
 30   call rdedx(cc,l3,iblock,num8)
c
c     calculate (i/mx/j) for mo's i nad j and mx=mx,my,mz
c     ri are diagonal elements
c     rij are the off diagonal
      do 80 kk = 1 , 3
         call rdedx(dmao,l2,ibl(kk),num8)
         mm = 0
         do 70 i = 1 , norb
            iim = ilifq(moout(i))
            do 60 j = 1 , i
               if (i.ne.j) mm = mm + 1
               jjm = ilifq(moout(j))
               sum = 0.0d0
               do 50 k = 1 , num
                  ccc = cc(iim+k)
                  do 40 l = 1 , num
                     kl = iky(max(k,l)) + min(k,l)
                     sum = sum + ccc*cc(l+jjm)*dmao(kl)
 40               continue
 50            continue
               if (i.eq.j) ri(i,kk) = sum
               if (i.ne.j) rij(mm,kk) = sum
 60         continue
 70      continue
 80   continue
      if (nredo.eq.0) qq = randd(1.0d0,cc,l1)
      write (iwr,6150) ytype(mi)
      write (iwr,6100)
      do 90 i = 1 , norb
         write (iwr,6110) moout(i) , (ri(i,k),k=1,3)
 90   continue
      do 110 i = 1 , norb
         do 100 j = 1 , i
            cl(i,j) = 0.00d0
            cl(j,i) = 0.00d0
 100     continue
         cl(i,i) = 1.00d0
 110  continue
      iter = 0
      shift = datan(1.00d0)
 120  change = 0.00d0
      iter = iter + 1
      do 130 i = 1 , norb
         iir(i) = i
 130  continue
      nnn = norb
      do 140 i = 1 , norb
         qq = randd(0.00d0,cc,l1)
         iii = idint(qq*dble(nnn)+1.00d0)
         iord(i) = iir(iii)
         iir(iii) = iir(nnn)
         nnn = nnn - 1
 140  continue
c     for each pair of orbitals a 2X2 unitary transformation
c     is performed. the transformation is  
c         psi prime(i)=   cos(t)*psi(i)+sin(t)*psi(j)  and    
c         psi prime (j)= -sin(t)*psi(i)+cos(t)*psi(j).
c     the boys method requires that  t be such as to maximize
c     the sum of the squares of the one center molecular orbital
c     dipole moment integrals.
      do 220 iii = 1 , norb
         i = iord(iii)
         if(norot(i,i).eq.1) go to 220
         im = i - 1
         jm = 1
         ijm = im*(im-1)/2 + 1
         rm = 0.00d0
         tm = 0.00d0
         sm = 0.00d0
         cm = 1.00d0
         do 180 j = 1 , norb
            if (i.lt.j) then
               ij = (j-1)*(j-2)/2 + i
            else if (i.eq.j) then
               go to 180
            else
               ij = im*(im-1)/2 + j
            end if
            if(norot(i,j).eq.1) go to 180
            t = 0.00d0
            tx = 0.00d0
            do 150 kk = 1 , 3
               t = t + 4.00d0*rij(ij,kk)**2 - ri(i,kk)**2 - ri(j,kk)
     +             **2 + 2.00d0*ri(i,kk)*ri(j,kk)
               tx = tx + rij(ij,kk)*(ri(j,kk)-ri(i,kk))
 150        continue
            if ((dabs(t).le.1.0d-10).and.(dabs(tx).le.1.0d-10))
     +          go to 180
            tx = 4.00d0*tx
            t = datan2(tx,t)
            t = t/4.00d0
            sign = 1.00d0
            if (t.gt.0.00d0) sign = -1.00d0
            t = t + sign*shift
            itim = 0
 160        s = dsin(t)
            itim = itim + 1
            co = dcos(t)
            rin = 0.00d0
            do 170 kk = 1 , 3
               qpi = co*co*ri(i,kk) + s*s*ri(j,kk) + 
     +          2.00d0*co*s*rij(ij,kk)
               qpj = co*co*ri(j,kk) + s*s*ri(i,kk) - 
     +          2.00d0*co*s*rij(ij,kk)
               rin = rin + qpi*qpi + qpj*qpj - ri(i,kk)**2 - ri(j,kk)**2
 170        continue
            ttest = dabs(t) - shift
            if ((dabs(t).gt.1.0d-08).and.(dabs(ttest).gt.1.0d-08)) then
               if (rin.lt.-1.0d-08) then
                  if (itim.le.1) then
                     sign = 1.00d0
                     if (t.gt.0.00d0) sign = -1.00d0
                     t = t + shift*sign
                     go to 160
                  else
                     write (iwr,6030) i , j
                     write (iwr,6040) t , s , co , rin
                     go to 230
                  end if
               end if
            end if
            if (rin.gt.rm) then
               ijm = ij
               rm = rin
               jm = j
               sm = s
               cm = co
               tm = t
            end if
 180     continue
         t = tm
         rin = rm
         s = sm
         co = cm
         j = jm
         ij = ijm
         if(norot(i,j).eq.1) go to 220
         change = change + t*t
         do 210 kk = 1 , 3
            qpi = co*co*ri(i,kk) + s*s*ri(j,kk) + 2.00d0*co*s*rij(ij,kk)
            qpj = co*co*ri(j,kk) + s*s*ri(i,kk) - 2.00d0*co*s*rij(ij,kk)
            qpij = (co*co-s*s)*rij(ij,kk) + co*s*(ri(j,kk)-ri(i,kk))
            do 190 k = 1 , norb
               if (i.lt.k) then
                  ik = (k-1)*(k-2)/2 + i
               else if (i.eq.k) then
                  go to 190
               else
                  ik = (i-1)*(i-2)/2 + k
               end if
               if (j.lt.k) then
                  jk = (k-1)*(k-2)/2 + j
               else if (j.eq.k) then
                  go to 190
               else
                  jk = (j-1)*(j-2)/2 + k
               end if
               qpix(k) = co*rij(ik,kk) + s*rij(jk,kk)
               qpjx(k) = co*rij(jk,kk) - s*rij(ik,kk)
 190        continue
            do 200 k = 1 , norb
               if (i.lt.k) then
                  ik = (k-1)*(k-2)/2 + i
               else if (i.eq.k) then
                  go to 200
               else
                  ik = (i-1)*(i-2)/2 + k
               end if
               if (j.lt.k) then
                  jk = (k-1)*(k-2)/2 + j
               else if (j.eq.k) then
                  go to 200
               else
                  jk = (j-1)*(j-2)/2 + k
               end if
               rij(ik,kk) = qpix(k)
               rij(jk,kk) = qpjx(k)
 200        continue
            rin = rin + qpi + qpj - ri(i,kk) - ri(j,kk)
            ri(i,kk) = qpi
            ri(j,kk) = qpj
            rij(ij,kk) = qpij
 210     continue
      call drot(norb,cl(1,i),1,cl(1,j),1,co,s)
 220  continue
c      if convergence has not been reached start another series
c      of two center rotations.
      change = dsqrt(2.00d0*change/(norb*(norb-1)))
      if (iter.le.maxboys) then
         if (iter.ge.1 .and. nprint.eq.(-7)) write (iwr,6050) iter ,
     +       change
         if (change.ge.1.0d-10) go to 120
      end if
 230  write (iwr,6060) iter,change
      if (iter.ge.maxboys .or. change.gt.1.d-6) then
         nredo = nredo + 1
         if (nredo.ne.3) then
            write (iwr,6070)
            go to 30
         end if
         write(iwr,6071)
      else
         write(iwr,6072)
      end if
      if (mi.eq.1) write (iwr,6080)
      if (mi.eq.2) write (iwr,6090)
      do 240 i = 1 , norb
         call dcopy(norb,cl(1,i),1,b(1,i),1)
 240  continue
      call rdedx(cc,l3,iblock,num8)
      write(10)norb
      do i=1,norb
         write(10)(b(j,i),j=1,norb)
      enddo
      do 260 i = 1 , norb
         call vclr(ri(1,1),1,num)
         do 250 j = 1 , norb
         iij = ilifq(moout(j))
         call daxpy(num,b(j,i),cc(iij+1),1,ri(1,1),1)
 250     continue
         call dcopy(num,ri(1,1),1,b(1,i),1)
 260  continue
c
c ----- now load up lmos into cc array
c
      do 270 i = 1 , norb
         ii = ilifq(moout(i))
         call dcopy(num,b(1,i),1,cc(ii+1),1)
 270  continue
c
c...   if we want harmonic do it now
c
      call prsql(cc,l1,l1,l1)
      if (oharm) then
         write(6,'(a)') 'Generating Harmonic local orbitals'
         if (mi.eq.1) then
           call harmonic
           i10 = igmem_alloc_inf(l3,'anala.m','boyorb',
     1                           'i10',IGMEM_DEBUG )
           call anorm(q(i10),q)
           call gmem_free_inf(i10,'anala.m','boyorb','i10')
         end if
         call tback(cc,ilifq,cc,ilifq,newbas0)
      end if
c
c ----- now save onto dumpfile
c
      if (mi.eq.1) call putq(zcom,ztitle,eiga,frocca,l1,l1,l1,m1,m1,cc,
     +                       nsav,imu)
      if (mi.eq.2) call putq(zcom,ztitle,eigb,froccb,l1,l1,l1,m1,m1,cc,
     +                       nsav,imu)
c ----
      call prsql(cc,l1,l1,l1)
      if (nprint.eq.7) then
         call pusql(cc,l1,l1,l1)
      end if
      write (iwr,6010) ytype(mi) , nsav
      do 310 imu = 1 , 3
         call rdedx(dmao,l2,ibl(imu),num8)
         do 300 i = 1 , norb
            sum = 0.0d0
            iim = ilifq(moout(i))
            do 290 k = 1 , num
               iik = iim + k
               do 280 l = 1 , num
                  kl = iky(max(k,l)) + min(k,l)
                  sum = sum + cc(iik)*cc(l+iim)*dmao(kl)
 280           continue
 290        continue
            ri(i,imu) = sum
 300     continue
 310  continue
      write (iwr,6140)
      write (iwr,6100)
      do 320 i = 1 , norb
         write (iwr,6110) moout(i) , (ri(i,k),k=1,3)
 320  continue
      write (iwr,6120) norb , num
c
 330  if (mi.ne.1 .or. nactb.eq.0) then
         cpu = cpulft(1)
         write (iwr,6020) cpu,charwall()
         return
      else
         iblock = ibl7qb
         nsav = moutb
         norb = nactb
         do 340 i = 1 , norb
            moout(i) = mooutb(i)
 340     continue
         go to 20
      end if
 6010 format (///1x,a3,' lmo output to section',i4,' of the dumpfile')
 6020 format (/' end of localization at ',f10.2,' seconds',a10,' wall',
     * //1x,104('-'))
 6030 format (/' ***** ',
     +        'no rotation increases integrals  --- program aborted'///
     +        10x,'i= ',i3,5x,'j= ',i3)
 6040 format (5x,'theta = ',g20.10/5x,'sin(theta)= ',f10.7,
     +        '   cos(theta)= ',f10.7/5x,
     +        'total change to this point = ',g20.10)
 6050 format (//10x,'iteration number',i3,'   change =',f20.10)
 6060 format (/25x,'program stopped after',i3,' iterations and change',
     1        1pe11.3/)
 6070 format (//10x,'*** localization has been unsucessful ***'
     +        //10x,'program will restart with new random number'
     +         /10x,'and rotation sequence for orbitals')
 6071 format(10x,'*** localisation failed ***')
 6072 format(10x,'*** localisation converged ***')
 6080 format (//50x,14('*')/50x,'alpha spin lmo'/50x,14('*')/)
 6090 format (//50x,14('*')/50x,' beta spin lmo'/50x,14('*')/)
 6100 format (//'  contributions to the electric dipole moment'/2x,
     +        43('-')//' m.o.',14x,'x',15x,'y',15x,'z'/1x,51('-')/)
 6110 format (1x,i3,5x,3(5x,f10.5))
 6120 format (//10x,'this case had ',i3,' active mos. and ',i3,
     +        ' basis functions'/)
 6130 format (/' commence boys localization at ',f8.2,' seconds',
     * a10,' wall')
 6140 format (//' **** final m.o.s'/)
 6150 format (/' **** starting m.o.s  ----',3x,a3,' vectors ')
      end
**==denhf.f
      subroutine denhf(zscftp,da,db,l,iblka,iblkb)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension da(*),db(*)
c
      data zrhf,zcas,zmcscf/'rhf','casscf','mcscf'/
c
      call rdedx(da,l,iblka,num8)
      if (zscftp.ne.zrhf.and.zscftp.ne.zcas.and.
     +    zscftp.ne.zmcscf) then
         call rdedx(db,l,iblkb,num8)
      end if
      return
      end
**==denint.f
      subroutine denint
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junk/pint,qint,rint,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj
c
      data done /1.0d0/
c
      pint = done
      qint = done
      rint = done
      ptxi = p0 - pi
      ptyi = q0 - qi
      ptzi = r0 - ri
      ptxj = p0 - pj
      ptyj = q0 - qj
      ptzj = r0 - rj
      go to (60,50,40,30,20,10,5 ) , ni
  5   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 10   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 20   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 30   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 40   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 50   pint = pint*ptxi
      qint = qint*ptyi
      rint = rint*ptzi
 60   continue
      go to (110,100,90,80,75,70,65) , nj
 65   pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 70   pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 75   pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 80   pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 90   pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 100  pint = pint*ptxj
      qint = qint*ptyj
      rint = rint*ptzj
 110  return
      end
**==dipmat.f
      subroutine dipmat(zscftp,q)
c
c     ----- restore the dipole matrix elements -----
c           density, eigenvectors and eigenvalues
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      common/restrl/ ociopt(2),omp2
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
      common/blkcore/corev(512),array(5),cpad(5)
      common/multic/radius(40),irad(57+mcfzc),itype(maxorb),isecn
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      dimension q(*),o1e(6)
c
      data zrhf/'rhf'/
      data zuhf/'uhf'/
      data zgvb/'gvb'/
      data zgrhf/'grhf'/
      data zcas /'casscf'/
      data zmcscf /'mcscf'/
      data zvb /'vb'/
      data m1,m19,m990/1,19,990/
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      len = lensec(l1)
c
c     ----- set pointers for partitioning of core -----
c
      omcscf = zscftp.eq.zmcscf.or.zscftp.eq.zvb
      if (zscftp.eq.zvb) isecn=mouta
      if (omcscf) then
        if (isecn.eq.0) then
           write (iwr,6037)
           return
        end if
        i10 = igmem_alloc_inf(5*l2,'anala.m','dipmat',
     1                        'i10',IGMEM_DEBUG )
        i20 = i10 + l2
        i30 = i20 + l2
        i40 = i30 + l2
        i50 = i40 + l2
      else
        i10 = igmem_alloc_inf(4*l2,'anala.m','dipmat',
     1                        'i10',IGMEM_DEBUG )
        i20 = i10 + l2
        i30 = i20 + l2
        i40 = i30 + l2
      endif

      if (nprint.eq.3) write (iwr,6070) i10 , i20 , i30 , i40
c
      ocas = zscftp.eq.zcas
      if (ocas) then
         call secloc(isecda,oexist,iblko)
         if (.not.oexist) then
            write (iwr,6035)
            write (iwr,6040)
         end if
      end if
      if (omp2) then
         call secloc(isect(101),oexist,iblko)
         if (oexist) then
            write (iwr,6010)
         else
            call secloc(isect(45),oexist,iblko)
            if (oexist) then
               write (iwr,6020)
            else
               write (iwr,6030)
               write (iwr,6040)
            end if
         end if
      end if
      if (mp3) then
         call secloc(isect(101),oexist,iblko)
         if (oexist) then
            write (iwr,6050)
         else
            write (iwr,6060)
            write (iwr,6040)
         end if
      end if
c
c     ---- first move vectors etc to ed7
c     ---- and transform from sabf to ao basis
c
c     retrieve and store a-density matrix
c
c     suppress header print in getq
c
      ndum = nprint
      nprint = -5
      ibl7pa = ibl7la
      if (omcscf) then
c
c     compute spinfree density matrix from NOS (see below)
c
      iposv = isecn
      else if (ocas) then
       call secget(isecda,m990,iblok)
       call rdedx(q(i10),l2,iblok,idaf)
       iposv = moutb
      else
       iposv = mouta
       call secget(isect(497),m19,iblok)
       ioffsa = iblok
       ioffsb = iblok + len + lensec(l2)
       call rdedx(q(i10),l2,ioffsa,idaf)
c      store
       if (zscftp.eq.zrhf) then
          if (omp2 .or. mp3) then
             call secloc(isecdd,oexist,iblok)
             if (.not.oexist .and. omp2)
     +           call secloc(isect(45),oexist,iblok)
             if (oexist) then
                call rdedx(q(i20),l2,iblok,idaf)
                do 20 ms = 1 , l2
                   q(i10-1+ms) = q(i10-1+ms) + q(i20-1+ms)
 20             continue
             end if
          end if
       end if
      endif
      if (omcscf) then
c
c     retrieve MCSCF NOs
c
       call getq(q(i10),q(i30),q(i40),nbasis,newbas,m1,m1,m1,
     +           iposv,zscftp)
       if (.not.otran)
     +   call tdown(q(i10),ilifq,q(i10),ilifq,l1)
       call dmtx(q(i50),q(i10),q(i40),iky,newbas,newbas,newbas)
       call wrt3(q(i50),l2,ibl7pa,num8)
       ibl7qa = iposun(num8)
       call wrt3(q(i10),l3,ibl7qa,num8)
      else
       call wrt3(q(i10),l2,ibl7pa,num8)
c
c     retrieve a-vectors
c
       call getq(q(i10),q(i30),q(i40),nbasis,newbas,m1,m1,m1,
     +           iposv,zscftp)
       call tdown(q(i10),ilifq,q(i10),ilifq,l1)
       ibl7qa = iposun(num8)
       call wrt3(q(i10),l3,ibl7qa,num8)
      endif
      ibl7ea = iposun(num8)
      if (ocas.or.omcscf.or.(ozora.and.oso)) then
c
c     store occupation numbers
c
       call wrt3(q(i40),l1,ibl7ea,num8)
      else
c
c     store a-eigenvalues
c
       call wrt3(q(i30),l1,ibl7ea,num8)
      end if
      if (.not.(ocas.or.omcscf).and.zscftp.ne.zrhf) then
c
c     retrieve and store b-density matrix
c
         call rdedx(q(i10),l2,ioffsb,idaf)
         ibl7pb = iposun(num8)
         call wrt3(q(i10),l2,ibl7pb,num8)
         if (zscftp.ne.zgvb .and. zscftp.ne.zgrhf) then
c
            if ((omp2 .or. mp3) .and. zscftp.eq.zuhf) then
               call secloc(isecdd,oexist,iblok)
               if (oexist) then
c ... note : adding ump2 correction to b-term
                  call rdedx(q(i20),l2,iblok,idaf)
                  do 30 ms = 1 , l2
                     q(i10-1+ms) = q(i10-1+ms) + q(i20-1+ms)
 30               continue
               end if
            end if
c
c     retrieve b-vectors
c
            call getq(q(i10),q(i30),q(i40),nbasis,newbas,m1,m1,m1,moutb,
     +                zscftp)
            call tdown(q(i10),ilifq,q(i10),ilifq,l1)
            ibl7qb = iposun(num8)
            call wrt3(q(i10),l3,ibl7qb,num8)
c
c     store b-eigenvalues
c
            ibl7eb = iposun(num8)
            if (.not.(ozora.and.oso)) then
               call wrt3(q(i30),l1,ibl7eb,num8)
            else
c...           for so store occupations
               call wrt3(q(i40),l1,ibl7eb,num8)
            end if
         end if
      end if
      nprint = ndum
      ibl7x = iposun(num8)
c
c     ----- restore dipole moment integrals -----
c
      do loop = 1,6
       o1e(loop) = .true.
      enddo
      o1e(2) = .false.
      o1e(3) = .false.
      call getmat(q(i40),q(i10),q(i10),q(i10),q(i20),q(i30),
     +            array,l1,o1e,ionsec)
c
c     ----- save x,y,z on ed7      -----
c
      call wrt3(q(i10),l2,ibl7x,num8)
      ibl7y = iposun(num8)
      call wrt3(q(i20),l2,ibl7y,num8)
      ibl7z = iposun(num8)
      call wrt3(q(i30),l2,ibl7z,num8)
      ibl7la = iposun(num8)
      call wrt3(q(i40),l2,ibl7s,num8)

      call gmem_free_inf(i10,'anala.m','dipmat','i10')

      return
 6010 format (/1x,'mp2 properties calculated as energy derivatives')
 6020 format (/1x,'mp2 properties calculated as expectation values')
 6030 format (/1x,'mp2 density matrix not found')
 6035 format (/1x,'casscf density matrix not found')
 6037 format (/1x,'mcscf NOS not found')
 6040 format (/1x,'analysing scf density matrix')
 6050 format (/1x,'mp3 properties calculated as energy derivatives')
 6060 format (/1x,'mp3 density matrix not found')
 6070 format (1x,'core assignment'/1x,'i10,  i20,  i30, i40  =',4i10)
      end
**==dipmo.f
      subroutine dipmo(zscftp,dipmu,v,rmu,nconf,l1,l2,l3,nprint)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
      dimension dipmu(l2,*),v(l1,*),rmu(3,*),nconf(*)
c
      data dzero,two/0.0d0,2.0d0/
      data zgvb/'gvb'/
      data fac/2.54158059d0/
c
c     ----- read in the mo's -----
c
      call rdedx(v,l3,ibl7qa,num8)
c
c     ----- calculate the number of occupied orbitals -----
c
      norb = nco + npair + npair
      if (nseto.ne.0) then
         do 20 i = 1 , nseto
            norb = norb + no(i)
 20      continue
      end if
      ic = 0
c
c     ----- set up map for orbitals to fock operators -----
c
      if (nco.ne.0) call setsto(nco,1,nconf)
      if (nseto.ne.0) then
         nbase = ncores
         ic = 0
         do 40 i = 1 , nseto
            nop = no(i)
            do 30 j = 1 , nop
               nconf(ic+nco+j) = nbase + 1
 30         continue
            ic = ic + nop
            nbase = nbase + 1
 40      continue
      end if
      if (npair.ne.0) then
         np2 = npair + npair
         do 50 i = 1 , np2
            nconf(i+nco+ic) = ncores + nseto + i
 50      continue
      end if
      call vclr(rmu,1,3*l1)
c
c     ----- loop over the norb orbitals -----
c
      do 90 i = 1 , norb
         do 80 j = 1 , 3
            kl = 1
            dum = 0.0d0
            do 60 k = 1 , l1
               if (v(k,i).ne.dzero) then
                  dum = dum + ddot(k,dipmu(kl,j),1,v(1,i),1)*v(k,i)
               end if
               kl = kl + k
 60         continue
            dum = two*dum
            kl = 1
            do 70 k = 1 , l1
               dum = dum - v(k,i)*dipmu(kl,j)*v(k,i)
               kl = kl + k + 1
 70         continue
            iham = nconf(i)
            rmu(j,i) = -two*f(iham)*dum*fac
 80      continue
 90   continue
      write (iwr,6010)
      do 100 i = 1 , norb
         write (iwr,6020) i , (rmu(j,i),j=1,3)
 100  continue
      return
 6010 format (/1x,80('=')//15x,32('-'),/15x,
     +        'molecular orbital dipole moments'/15x,32('-'),//,10x,
     +        '   mo   ',4x,'x',11x,'y',11x,'z',//)
 6020 format (10x,i5,3f12.6)
      end
**==dipole.f
      subroutine dipole(zscf,q)
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension q(*),xx(3),dmq(3),dmxe(3),dxt(3)
c
      data zrhf,zgvb,zgrhf/'rhf','gvb','grhf'/
      data zcas / 'casscf' /, zmcscf / 'mcscf'/
      data fac /2.54158059d0/
      data xx /'x','y','z'/
      data m3/3/
      data small/1.0d-16/
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c     ----- set pointers for partitioning of core -----
c
c
c measure required size
c
      i10 = 0
      i20 = i10 + l2
      i30 = i20 + l2
      i40 = i30 + l2
      i50 = i40 + l2
      i51 = i40 + l3
      i60 = i50 + l2
      i61 = i51 + l1
      i71 = i61 + l1 + l1 + l1
      i81 = i71 + l1
      last = i50
      if (zscf.ne.zrhf) last = i60
      if (zscf.eq.zgvb .or. zscf.eq.zgrhf) last = i81
      length = last - i10

      i10 = igmem_alloc_inf(length,'anala.m','dipole',
     1                      'i10',IGMEM_DEBUG )
      i20 = i10 + l2
      i30 = i20 + l2
      i40 = i30 + l2
      i50 = i40 + l2
      i51 = i40 + l3
      i60 = i50 + l2
      i61 = i51 + l1
      i71 = i61 + l1 + l1 + l1
      i81 = i71 + l1

      if (nprint.eq.3) write (iwr,6040) i10 , i20 , i30 , i40 , i50 ,
     +                                  i51 , i60 , i61 , i71 , i81 ,
     +                                  last
c
c     -----  read in dipole moment integrals -----
c
      call rdedx(q(i10),l2,ibl7x,num8)
      call rdedx(q(i20),l2,ibl7y,num8)
      call rdedx(q(i30),l2,ibl7z,num8)
c
c     ----- electronic contribution to dipole moment -----
c
      if (nprint.ne.(-5)) then
c
c     ----- first, on an mo basis if gvb or grhf-----
c
        if (zscftp.eq.zgvb.and.npair.ne.0) then
         call dipmo(zscf,q(i10),q(i40),q(i61),q(i71),l1,l2,l3,nprint)
        endif
      endif
c
c     ----- now the total electric dipole moment -----
c
      if (zscf.ne.zrhf.and.zscf.ne.zcas.and.
     +    zscf.ne.zmcscf) then
        call denhf(zscf,q(i40),q(i50),l2,ibl7pa,ibl7pb)
        dmxe(1) = -tracep(q(i40),q(i10),l1)
        dmxe(2) = -tracep(q(i40),q(i20),l1)
        dmxe(3) = -tracep(q(i40),q(i30),l1)
        dmxe(1) = -tracep(q(i50),q(i10),l1) + dmxe(1)
        dmxe(2) = -tracep(q(i50),q(i20),l1) + dmxe(2)
        dmxe(3) = -tracep(q(i50),q(i30),l1) + dmxe(3)
      else
        call denhf(zscf,q(i40),q(i40),l2,ibl7pa,ibl7pb)
        dmxe(1) = -tracep(q(i40),q(i10),l1)
        dmxe(2) = -tracep(q(i40),q(i20),l1)
        dmxe(3) = -tracep(q(i40),q(i30),l1)
      endif
c
c     ----- nuclear contribution (excluding bq's when xfield
c           option is on)
c
      do   1 i=1,3
         dmq(i)=0.0d0
 1    continue
      do   2 i = 1,nat
         if( .not. ocryst .or. (zaname(i)(1:2) .ne. 'bq')) then
            dmq(1) = dmq(1)+czan(i)*c(1,i)
            dmq(2) = dmq(2)+czan(i)*c(2,i)
            dmq(3) = dmq(3)+czan(i)*c(3,i)
         endif
 2    continue

      do 20 i = 1 , 3
         dxt(i) = dmq(i) + dmxe(i)
 20   continue
      dtot = ddot(m3,dxt,1,dxt,1)
      if (dtot.gt.small) dtot = dsqrt(dtot)
      dipol = dtot*fac
      if (nprint.ne.-5) then
         if(ocryst)then
            write (iwr,6008)
         else
            write (iwr,6009)
         endif
         write (iwr,6010)
         do 30 i = 1 , 3
            write (iwr,6020) xx(i) , dmq(i) , dmxe(i) , dxt(i)
 30      continue
         write (iwr,6030) dtot , dipol
      end if
c
c  write to punchfile if requested
c
      call blkdip(dxt)
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,'anala.m','dipole','i10')

      return
 6008 format (/17x,'dipole moments excluding bq centres')
 6009 format (/17x,'dipole moments')
 6010 format (//11x,'nuclear',6x,'electronic',11x,
     +        'total'/)
 6020 format (1x,a1,3f16.7)
 6030 format (/' total dipole moment = ',f16.7,
     +        ' (a.u.)'/'                       ',f16.7,' (debye)')
 6040 format (' core assignement ',/,' i10, i20, i30, i40, i50, i51,',
     +        ' i60, i61, i71, i81 = ',10i8,/,' last = ',i8)
      end
**==getqp.f
      subroutine getqp(zcomm,ztit,eig,deff,norb,norbn,ncolu,ieig,
     *ideff,nomx,iblk)
c...     read eigenvectors
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
      parameter (mxorb1 = maxorb+1)
      dimension zcomm(*),ztit(*),eig(*),deff(*)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
      common/blkorbs/value(maxorb),pop(mxorb1),
     *nbasis,newb,ncol,ivalue,ipop,ift
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/junkc/zcom(19),ztitle(10)
c
      data m29/29/
c
      call rdchr(zcom,m29,iblk,numdu)
      call reads(value,mach(8),numdu)
      if (nbasis.ne.nbfnd .or. ncol.le.0 .or. ncol.gt.nomx) then
       write(iwr,40) nbasis, ncol, nbfnd,nomx
 40    format(/1x,'*** retrieved nbasis,ncol = ',2i6/
     +         1x,'*** checkpoint values     = ',2i6)
       call caserr2('retrieved eigenvectors are invalid')
      end if
      norb = nbasis
      norbn = newb
      ncolu = ncol
      ieig = ivalue
      ideff = ipop
      do 20 i = 1 , 10
         ztit(i) = ztitle(i)
 20   continue
      do 30 i = 1 , 19
         zcomm(i) = zcom(i)
 30   continue
      call dcopy(ncolu,value,1,eig,1)
      if (ideff.ge.0) then
         call dcopy(ncolu,pop,1,deff,1)
      end if
      nav = lenwrd()
      call readis(ilifc,mach(9)*nav,numdu)
      iblk = iblk + 1 + lensec(mach(8)) + lensec(mach(9))
      return
      end
**==gnorot.f
      subroutine gnorot(zscftp,ipass,norot,m1,nact,odef)
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
c
      dimension norot(m1,*)
      dimension numcor(103)
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
      data numcor/2*0,   
     *            2*1,                  6*1,  
     *            2*5,                  6*5,
     *            2*9,          10*9,   6*14,
     *            2*18,         10*18,  6*23,
     *            2*27,  14*27, 10*34,  6*39,
     *            2*43,  14*43,    50/
      data zrhf,zuhf,zrohf/'rhf','uhf','rohf' /
c
c
c     ----- set up defaults for the wavefunction -----
c
c     -mcore- is the number of chemical core (nonvalence) orbitals,
c     -mdoc- is the number of doubly occupied valence orbitals,
c     -mact- is the number of partially occupied orbitals.
c     for uhf, -mcore- is still the chemical core, but -mdoc- is the
c     number of singly occupied valence orbitals, of alpha or beta
c     spin depending on the value for -ipass-.
c
      if(.not.odef) then
c
      do 160 i=1,nact
      do 160 j=1,nact
         norot(j,i) = 0
  160 continue
c
      else
c    
      mcore= 0
      mdoc = 0
c
c     first count the chemical core orbitals.
c
      do 100 i=1,nat
         iz = int(czan(i)+0.01d0) 
         if(iz.le.0) go to 100
         if(iz.gt.103) call caserr2(
     +    'z.gt.103 in establishing core count')
         mcore = mcore + numcor(iz) 
  100 continue
c
c     assign valence orbitals according to the scf type.
c
      if(zscftp.eq.zrhf) then
         mdoc = na - mcore
         nact = mdoc
      end if
      if(zscftp.eq.zuhf) then
c
       if(ipass.eq.1) then
         mdoc = na - mcore
         nact = mdoc
       else
         mdoc = nb - mcore
         nact = mdoc
       endif
c
      end if
      if(zscftp.eq.zrohf) then
         mdoc = nb - mcore
      end if
c
c     ----- now set up defaults for freezing orbitals -----
c
      do 120 i=1,nact
      do 125 j=1,nact
         norot(j,i) = 0
  125 continue
  120 continue
c
c     set up an invariant localization, by restricting any
c     rotations mixing doubly and partially occupied orbitals.
c     leave out any row and column involving a frozen orbital.
c
      ndbocc = mcore + mdoc
      ij=0
       do 280 i=1,nact
       do 270 j=1,i
          ij = ij+1
          if(i.gt.ndbocc  .and.  j.le.ndbocc) norot(i,j)=1
  270  continue
  280  continue
c
      endif
c
      return
      end
**==grossc.f
      subroutine grossc(a,b,n)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
c
      data dzero /0.0d0/
c
      ind(i,j) = ((max(i,j)*(max(i,j)-1))/2+min(i,j))
c
      do 30 i = 1 , n
         dum = dzero
         do 20 j = 1 , n
            ij = ind(i,j)
            dum = dum + a(ij)
 20      continue
         b(i) = dum
 30   continue
c
      return
      end
**==hfprop.f
      subroutine hfprop(zzzzzz,core)
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
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
      real*8 aoccc
      integer isetc, iync, ns_canon, isec_canon
      logical opr_canon,o_canon,oset_canon
      common/canon_nat/aoccc(10),isetc(10),iync(10),
     1                 ns_canon,isec_canon,
     2                 opr_canon,o_canon,oset_canon
c
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/restrl/ociopt(2),omp2
      common/multic/radius(40),irad(57+mcfzc),itype(maxorb),isecn
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
      character*10 charwall
      dimension core(*)
c
      data zscf,zprop/'scf','prop'/
      data zmcscf,zcas,zvb/'mcscf','casscf','vb'/
c
      if (zzzzzz.eq.zvb) zzzzzz = zcas
      call cpuwal(begin,ebegin)
      if (nprint.ne.-5) then
         write (iwr,6010) begin,charwall()
      end if


c
c     ----- restore the dipole matrix elements -----
c
      call dipmat(zzzzzz,core)
      if(zzzzzz.eq.zmcscf.and.isecn.eq.0) return
c
      if (ilopri .gt. 0) call mylowd(core)
c
c     ----- calculate the d transformation matrix and then transform
c           all affected matrices -----
c
      if (zruntp.eq.zscf .or. zruntp.eq.zprop) then
        call transd(zzzzzz,core,nprint)
      endif
c
c     ----- dipole moment -----
c
      call dipole(zzzzzz,core)
      if (nprint.ne.(-5)) then
         if (.not.oprint(26)) then
c
c     ----- symmetrically transform vectors and density matrix -----
c
            call lowdin(zzzzzz,core,nprint)
c
c     ----- mulliken and lowdin population analysis -----
c
            call mulken(zzzzzz,core,nprint)
c
c     ----- atomic spin density -----
c
            if (zzzzzz.ne.zcas.and.zzzzzz.ne.zmcscf) 
     +      call spind(zzzzzz,core)
c
            if (omp2. or. mp3) then
c
c     ----- MP2 Natural orbitals -----
c
             if (iuno.ge.0) call mp2nat(core,iuno,iunopp,0,zzzzzz)
             if (iunosp.ge.0) call mp2nat(core,iunosp,iunspp,1,zzzzzz)
c
            else if (zzzzzz.eq.'uhf') then
c
c     ----- UHF Natural orbitals -----
c
             if (iuno.ge.0) call uhfnat(core,iuno,iunopp,0)
             if (iunosp.ge.0) call uhfnat(core,iunosp,iunspp,1)
            else
            endif
c
            if (o_canon) then
c
c ----- canonicalisation of natural orbitals
c
c
               iv = igmem_alloc_inf(num*num,'anala.m','hfprop',
     1                              'iv',IGMEM_DEBUG )
               call excanon(core(iv),num,iuno,core)
               call gmem_free_inf(iv,'anala.m','hfprop','iv')
c
            end if
c           
            cpu = cpulft(1)
            write (iwr,6020) cpu,charwall()
         end if
      end if
      call timana(14)
      return
 6010 format (/1x,80('=')//40x,21('*')/40x,'wavefunction analysis'/40x,
     +        21('*')//' commence analysis at ',f9.2,' seconds'
     +        ,a10,' wall')
 6020 format (//' end of wavefunction analysis at ',f8.2,
     +        ' seconds.',a10,' wall'//1x,80('=')//)
      end
**==hfprop2.f
      subroutine hfprop2(zzzzzz,core)
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
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
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
      common/restrl/ociopt(2),omp2
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
      character*10 charwall
      dimension core(*)
c
c
      call cpuwal(begin,ebegin)
      if (nprint.ne.-5) then
         write (iwr,6010) begin,charwall()
      end if
      if (nprint.ne.(-5)) then
         if (.not.oprint(26)) then
            call hfprop3(zzzzzz,core)
            cpu = cpulft(1)
            write (iwr,6020) cpu,charwall()
         end if
      end if
      call timana(14)
      return
 6010 format (/1x,80('=')//40x,17('*')/40x,'property analysis'/40x,
     +        17('*')//' commence analysis at ',f9.2,' seconds'
     +        ,a10,' wall')
 6020 format (//' end of property analysis at ',f8.2,
     +        ' seconds',a10,' wall'//1x,80('=')//)
      end
**==latpop.f
      subroutine latpop(dd,da,nat)
c
c     - compress from orbitals to atoms for lowdin analysis -----
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
      common/junk/limlow(maxat),limsup(maxat)
      dimension dd(*),da(*)
c
      data dzero /0.0d0/
c
      do 30 i = 1 , nat
         da(i) = dzero
         i1 = limlow(i)
         i2 = limsup(i)
         if (i1.ne.0) then
            do 20 j = i1 , i2
               da(i) = da(i) + dd(j)
 20         continue
         end if
 30   continue
c
      return
      end
**==lmopsi.f
      subroutine lmopsi(zscftp,mcore,mdoc,mact,mocc)
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
      dimension numcor(103)
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
      common/miscop/nacta,nactb,moouta(maxorb),mooutb(maxorb)
     * , eiga(maxorb),frocca(maxorb),eigb(maxorb),froccb(maxorb)
c
      data numcor/2*0,   
     *            2*1,                  6*1,  
     *            2*5,                  6*5,
     *            2*9,          10*9,   6*14,
     *            2*18,         10*18,  6*23,
     *            2*27,  14*27, 10*34,  6*39,
     *            2*43,  14*43,    50/
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb' /
      data zatext,zbtext/'a','b'/
c
c
c     ----- set up defaults for the wavefunction -----
c
c     -mcore- is the number of chemical core (nonvalence) orbitals,
c     -mdoc- is the number of doubly occupied valence orbitals,
c     -mact- is the number of partially occupied orbitals.
c
c     for uhf, -mcore- is still the chemical core, but -mdoc- is the
c     number of singly occupied valence orbitals, of alpha or beta
c     spin 
c
      mcore= 0
      mdoc = 0
      mact = 0
c
c     first count the chemical core orbitals.
c
      do 100 i=1,nat
         iz = int(czan(i)+0.01d0) 
         if(iz.le.0) go to 100
         if(iz.gt.103) call caserr2(
     +    'z.gt.103 in establishing core count')
         mcore = mcore + numcor(iz) 
  100 continue
c
c     assign valence orbitals according to the scf type.
c
      if(zscftp.eq.zrhf) then
         mdoc = na - mcore
         mact = 0
         nacta = mdoc
         do 130 i = 1,nacta
130      moouta(i) = mcore + i
      end if
      if(zscftp.eq.zuhf) then
c
         mdoc = na - mcore
         nacta = mdoc
         do 140 i = 1,nacta
140      moouta(i) = mcore + i
c
         mdoc = nb - mcore
         nactb = mdoc
         do 150 i = 1,nactb
150      mooutb(i) = mcore + i
c
         mact = 0
      end if
      if(zscftp.eq.zgrhf.or.zscftp.eq.zgvb) then
         mdoc = nb - mcore
         mact = na - nb
      end if
c
      mocc = mcore + mdoc + mact
c
      if (nacta.gt.0) then
       if (nacta.le.num) go to 290
 280   call caserr2(
     + 'error detected in default list of orbitals to be localised')
 290   do 300 i = 1 , nacta
          if (moouta(i).gt.num) go to 280
 300   continue
       write (iwr,6800) zatext
       write (iwr,6810) (moouta(i),i=1,nacta)
      endif
      if (nactb.gt.0) then
         if (nactb.le.num) go to 340
 350   call caserr2(
     + 'error detected in default list of orbitals to be localised')
 340     do 320 i = 1 , nactb
            if (mooutb(i).gt.num) go to 350
 320     continue
         write (iwr,6800) zbtext
         write (iwr,6810) (mooutb(i),i=1,nactb)
      end if
      return
c
 6800 format (/1x,'* list of active  -',a1,'-  mos')
 6810 format (/20i4)
c
      end
**==lowdin.f
      subroutine lowdin(zscftp,q,nprint)
c
c     ----- lowdin population analysis -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical otri
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension q(*),ibl(6)
c
      data zrhf,zuhf,zcas,zmcscf/'rhf','uhf','casscf','mcscf'/
c
      l1 = num
      l2 = num*(num+1)/2
      l3 = num*num
c
c     ----- set pointers for partioning of core -----
c
      otri = .false.
      ntri = 6
      length = l1 +  ntri*l2
      i10 = igmem_alloc_inf(length,'anala.m','lowdin',
     1                         'i10',IGMEM_DEBUG )
      i20 = i10 + l2
      i30 = i20 + l2
      i31 = i30 + l3
      i40 = i30 + l2
      if (otri) then
       i50 = i10
       i60 = i50 + l2
       i70 = i40 + l2
      else
       i50 = i40 + l2
       i60 = i50 + l2
       i70 = i60 + l2
      endif
      i80 = i70 + l1

      if (nprint.eq.3) write (iwr,6010) i10 , i20 , i30 , i31 ,
     +                                  i50 , i60 , i70 , i80
c
c
c     ----- compute the symmetric transformation matrix
c
c     -s - at q(i10)     overlap matrix
c     -sv- at q(i30)     eigenvectors of s
c     -se- at q(i31)     eigenvalues of s
c     -sh- at q(i10)     symmetric transformation matrix
c
      ibl7so = ibl7la
      call rdedx(q(i10),l2,ibl7s,num8)
      call symtrn(q(i10),q(i30),q(i31),iky,q(i70),l1,l2,l1,'lowdin')
      call wrt3(q(i10),l2,ibl7so,num8)
c
c ----- define block positions on ed7 for orthogonal p,q,e
c
      lenq = lensec(l3)
      lenp = lensec(l2)
      lenv = lensec(l1)
      ibl7o(1) = iposun(num8)
      ibl7o(2) = ibl7o(1) + lenq
      ibl7o(3) = ibl7o(2) + lenp
      ibl7o(4) = ibl7o(3) + lenv
      ibl7o(5) = ibl7o(4) + lenq
      ibl7o(6) = ibl7o(5) + lenp
c
c     ----- transform the vectors and the density matrix
c           to the symmetric orthogonal basis set
c
c     -sh- at q(i30)     symmetric transformation matrix
c     -v - at q(i10)
c     -vo- at q(i10)     vectors in orthogonal basis set
c     -vo2-at q(i30)     orbital population matrix
c                        vo2(i,j)= vo(i,j)**2*sign(vo(i,j))
c     -do- at q(i30)     density matrix in the orthogonal basis
c     -dod-at q(i60)     diagonal of do
c
      ndav = 1
      ibl(1) = ibl7qa
      ibl(2) = ibl7pa
      ibl(3) = ibl7ea
      ibl(4) = ibl7qb
      ibl(5) = ibl7pb
      ibl(6) = ibl7eb
c

      it = 0
 20   if (it.ne.1 .or. zscftp.eq.zuhf) then
c     ----- calculate orbital populations -----
         call rdedx(q(i30),l2,ibl7so,num8)
         call rdedx(q(i10),l3,ibl(ndav),num8)
         call symtrv(q(i30),q(i10),q(i30),q(i50),l1,l1,l1)
         call wrt3(q(i30),l3,ibl7o(ndav),num8)
         if (oprint(42)) then
            if (zscftp.ne.zrhf .and. it.eq.0) write (iwr,6020)
            if (zscftp.ne.zrhf .and. it.eq.1) write (iwr,6030)
            if (zscftp.eq.zrhf) write (iwr,6040)
         end if
         call rdedx(q(i70),l1,ibl(ndav+2),num8)
         if (oprint(42)) call prev(q(i50),q(i70),l1,l1,l1)
      end if
      call rdedx(q(i10),l2,ibl7so,num8)
      call rdedx(q(i20),l2,ibl(ndav+1),num8)
      call symtrd(q(i10),q(i20),q(i30),q(i50),q(i10),q(i70),l1)
      call wrt3(q(i10),l2,ibl7o(ndav+1),num8)
      call wrt3(q(i70),l1,ibl7o(ndav+2),num8)
      it = it + 1
      if (it.ge.2 .or. zscftp.eq.zrhf.or.zscftp.eq.zcas.or.
     +                 zscftp.eq.zmcscf) then
         ibl7la = iposun(num8)

         call gmem_free_inf(i10,'anala.m','lowdin','i10')
         return
      else
         ndav = ndav + 3
         go to 20
      end if
 6010 format (' core assignement'/' i10, i20, i30, i31, i50, i60',
     +        ' i70, i80 = '/7i8/' last = ',i8)
 6020 format (//1x,80('=')//10x,58('-')/10x,'lowdin orbital (vo(i,j)',
     +        '**2)*','sign(vo(i,j))   alpha orbitals'/10x,58('-'))
 6030 format (//1x,80('=')//10x,57('-')/10x,'lowdin orbital (vo(i,j)',
     +        '**2)*','sign(vo(i,j))   beta orbitals'/10x,57('-'))
 6040 format (//1x,80('=')//10x,41('-')/10x,'lowdin orbital (vo(i,j)',
     +        '**2)*','sign(vo(i,j))'/10x,41('-'))
 6050 format(/1x,' **** 4-triangle GA-based code invoked *****'/)
c
      end
**==mopop.f
      subroutine mopop(s,v,t,ia,zscftp,l1,l3)
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
      dimension s(*),v(l1,*),t(*),ia(*)
c
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
      data zcas,zmcscf/'casscf','mcscf'/
      data zalph,zbet,zboth/'alpha','beta','both'/
      data dzero,done,two/0.0d0,1.0d0,2.0d0/
c
      norb = na
      if (zscftp.eq.zgrhf .or. zscftp.eq.zgvb) then
c
c     ----- if grhf or gvb, calculate norb -----
c
         norb = nco + npair + npair
         if (nseto.ne.0) then
            do 20 i = 1 , nseto
               norb = norb + no(i)
 20         continue
         end if
      end if
      if(zscftp.eq.zcas.or.zscftp.eq.zmcscf.or.(ozora.and.oso))norb = l1
      ipass = 1
      fact = done
      if (zscftp.eq.zrhf) fact = two
      mblq = ibl7qa
      mble = ibl7ea
      if (oprint(41)) write (iwr,6010)
 30   zlab = zboth
      if (ipass.eq.1 .and. zscftp.ne.zrhf) zlab = zalph
      if (ipass.eq.2 .and. zscftp.ne.zrhf) zlab = zbet
      if (oprint(41)) write (iwr,6020) zlab
      call rdedx(v,l3,mblq,num8)
      do 60 i = 1 , norb
         do 50 j = 1 , l1
            sum = dzero
            do 40 k = 1 , l1
               jk = ia(max(j,k)) + min(j,k)
               sum = sum + v(k,i)*s(jk)*v(j,i)
 40         continue
            t(j) = sum*fact
 50      continue
         call dcopy(l1,t,1,v(1,i),1)
 60   continue
c
      call rdedx(t,l1,mble,num8)
c
      if (zscftp.eq.zcas.or.zscftp.eq.zmcscf.or.(ozora.and.oso)) then
c
c     ----- if casscf weight by occupation number -----
c
      do 65 i = 1,norb
 65   call dscal(l1,t(i),v(1,i),1)
c
      endif
c
      if (zscftp.eq.zgrhf .or. zscftp.eq.zgvb) then
c
c     ----- if grhf or gvb, weight by occupation number -----
c
         iorb = nco
         iham = ncores
         if (nco.gt.0) then
            sum = two*f(iham)
            ncol1 = nco*l1
            call dscal(ncol1,sum,v(1,1),1)
         end if
         if (nseto.gt.0) then
           do 70 i = 1 , nseto
           nop = no(i)
           nopl1 = nop*l1
           iham = iham + 1
           sum = two*f(iham)
           call dscal(nopl1,sum,v(1,iorb+1),1)
           iorb = iorb + nop
 70        continue
         end if
         if (npair.gt.0) then
            np2 = npair + npair
            do 80 i = 1 , np2
               iham = iham + 1
               iorb = iorb + 1
               sum = two*f(iham)
               call dscal(l1,sum,v(1,iorb),1)
 80         continue
         end if
c
c     ----- if grhf or gvb, zero out the virtual eigenvalues -----
c
         norbp1 = norb + 1
         if (norbp1.le.l1) then
            call vclr(t(norbp1),1,l1-norbp1+1)
         end if
      end if
      if (oprint(41)) call prev(v,t,norb,l1,l1)
      if (zscftp.ne.zuhf .or. ipass.eq.2) then
         return
      else
         ipass = 2
         mblq = ibl7qb
         mble = ibl7eb
         norb = na
         go to 30
      end if
 6010 format (/1x,80('=')//10x,44('-')/10x,
     +        'mulliken populations over ','molecular orbitals'/10x,
     +        44('-')/)
 6020 format (//10x,a8,' orbitals'/)
c
      end
**==mp2nat.f
      subroutine mp2nat(q,nsav,iprinv,ispin,zzzzzz)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c..   get spinfree natural orbitals from mp2 function
c..   ispin = 0 : spinfree natural orbitals
c
c..
      dimension q(*)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      common/restrl/ociopt(2),omp2
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      dimension ztype(2)
c
      data ztype/'spinfree','  spin  '/
      data zrhf/'rhf'/
      data zuhf/'uhf'/
      data zgvb/'gvb'/
      data zgrhf/'grhf'/
      data m19/19/
c..
c..    core partitioning
c..
       l1 = num
       l2 = num*(num+1)/2
       l3 = num*num
       len = lensec(l1)
c
       i10 = igmem_alloc_inf(3*l3+2*l1,'anala.m','mp2nat',
     1                       'i10',IGMEM_DEBUG )
       i20 = i10 + l3
       i30 = i20 + l3
       i50 = i30 + l3
       i60 = i50 + l1
c      i50 is scratch for minv (crayxmp) : 2*l1
c      i50,i60 scratch for minvrt (others) l1 each
c      last = i60 + l1
c
c     ----- get mp2 density matrix -----
c
      if (omp2) then
         call secloc(isect(101),oexist,iblko)
         if (oexist) then
            write (iwr,6010)
         else
            call secloc(isect(45),oexist,iblko)
            if (oexist) then
               write (iwr,6020)
            else
               write (iwr,6030)
               write (iwr,6040)
               go to 10
            end if
         end if
      end if
      if (mp3) then
         call secloc(isect(101),oexist,iblko)
         if (oexist) then
            write (iwr,6050)
         else
            write (iwr,6060)
            write (iwr,6040)
            go to 10
         end if
      end if
c
c     -da- at q(i10)  and  -db- at q(i20)
c
c    then total at (i10)
c    ispin 0 : +   // ispin 1 : -
c
      call secget(isect(497),m19,iblok)
      ioffsa = iblok
      ioffsb = iblok + len + lensec(l2)
      call rdedx(q(i10),l2,ioffsa,idaf)
c
      if (zzzzzz.eq.zrhf) then
         call secloc(isecdd,oexist,iblok)
         if (.not.oexist .and. omp2)
     +       call secloc(isect(45),oexist,iblok)
         if (oexist) then
            call rdedx(q(i20),l2,iblok,idaf)
            do ms = 1 , l2
               q(i10-1+ms) = q(i10-1+ms) + q(i20-1+ms)
            enddo
         end if
      end if
c
      if (zzzzzz.ne.zrhf) then
c
c     retrieve and store b-density matrix
c
         call rdedx(q(i20),l2,ioffsb,idaf)
         if (zzzzzz.ne.zgvb .and. zzzzzz.ne.zgrhf) then
c
            if (zzzzzz.eq.zuhf) then
               call secloc(isecdd,oexist,iblok)
               if (oexist) then
c ... note : adding ump2 correction to b-term
                  call rdedx(q(i30),l2,iblok,idaf)
                  do ms = 1 , l2
                     q(i20-1+ms) = q(i20-1+ms) + q(i30-1+ms)
                  enddo
               end if
            end if
c
         end if
         if (ispin.eq.0) then
            call vadd(q(i10),1,q(i20),1,q(i10),1,l2)
         else
            call vsub(q(i10),1,q(i20),1,q(i10),1,l2)
         end if
      end if
c
c..   transform density matrix to orthonormal basis (scf MOs)
c
      call rdedx(q(i20),l3,ibl3qa,idaf)
      call tdown(q(i20),ilifq,q(i20),ilifq,l1)
      call vclr(q(i50),1,2*l1)
c
c      invert transformation matrix t -> t(inv)
c
      d = 0.0d0
      call minvrt(q(i20),l1,d,q(i50),q(i60))
c
c      transpose it  t(inv) -> t(inv(dagger))
c
      call dagger(l1,l1,q(i20),l1,q(i30),l1)
c
c      then make t(inv).d.t(inv(dagger))
c      i30: t / i20: output d / i10: input d
c
      call mult2(q(i30),q(i20),q(i10),l1,l1,l1)
c
c..    diagonalise orthonormalised d-matrix (sort in decreasing order)
c..    result : occupations in i50 ; vectors in i30
c
      call gldiag(l1,l1,l1,q(i20),q(i60),q(i50),q(i30),iky,3)
c
c...  back-transform the eigenvectors (with symmetry adapted vectors)
c
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call tfsqc(q(i30),q(i10),q(i20),l1,l1,l1)
c
c...   save ** mo's and eigenvalues = occupations ** on section nsav
c
      if (nsav.ne.0) then
         if (iprinv.gt.0) write(iwr,9027) ztype(ispin+1),nsav
         call putq(zcom,ztitle,q(i50),q(i50),l1,l1,l1,1,1,
     *             q(i30),nsav,ibl3qa)
      end if
c
c..     transform back to original basis and print if requested
c
      if (iprinv.gt.1) then
         call tdown(q(i30),ilifq,q(i30),ilifq,l1)
         write(iwr,9028) ztype(ispin+1)
         call prev(q(i30),q(i50),l1,l1,l1)
      endif
      if (iprinv.eq.1) then
         write(iwr,9029) ztype(ispin+1)
         write(iwr,9030) (q(i50+i-1),i=1,l1)
         write(iwr,9030)
      end if
c
  10  continue
      call gmem_free_inf(i10,'anala.m','mp2nat','i10')
c
 6010 format (/1x,'mp2 properties calculated as energy derivatives')
 6020 format (/1x,'mp2 properties calculated as expectation values')
 6030 format (/1x,'mp2 density matrix not found')
 6040 format (/1x,'**** mp2 natural orbitals not generated')
 6050 format (/1x,'mp3 properties calculated as energy derivatives')
 6060 format (/1x,'mp3 density matrix not found')
 9027 format (/' ** ',a8,' MP2 natural orbitals saved in section'
     *               ,i4,' **')
 9028 format (/30x,29('-')/
     +         30x,a8,' mp2 natural orbitals '/
     +         30x,29('-')//)
 9029 format  (//10x,'----- ',a8,' MP2 natorb occupations -----')
 9030 format  (/10x,8f14.7)
      return
      end
**==mulken.f
      subroutine mulken(zzzzzz,q,nprint)
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
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension q(*)
c

      data zalpha,zbeta,ztotal,zspin /'   alpha','    beta','     all',
     +     '    spin'/
      data zrhf,zcas,zmcscf/'rhf','casscf', 'mcscf'/
      data zskip/'skippop'/
c
      out = nprint.eq.3
      ohfcas = zzzzzz.eq.zrhf.or.zzzzzz.eq.zcas.or.
     +         zzzzzz.eq.zmcscf
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c     ----- set pointers for partitioning of core -----
c
      length = 4*nat + 3*l2 + (nat*(nat+1))/2 + 4*l1
      i10 = igmem_alloc_inf(length,'anala.m','mulken',
     1                      'i10',IGMEM_DEBUG )

      i30 = i10 + l2
      i40 = i30 + (nat*(nat+1))/2
      i50 = i40 + l2
      i60 = i50 + l1
      i70 = i60 + nat
      i80 = i70 + l1
      i90 = i80 + nat
      i41 = i90
      i51 = i41 + l2
      i61 = i51 + l1
      i71 = i61 + nat
      i81 = i71 + l1
      i91 = i81 + nat

      if (nprint.eq.3) write (iwr,6090) i10 , i30 , i40 , i50 , i60 ,
     +                                  i70 , i41 , i51 , i61 , i71 ,
     +                                  i91
c
c     ----- determine number of basis functions per atom -----
c
      call aolim
c
c     ----- read in overlap matrix -----
c     -s- at q(i10)
      call rdedx(q(i10),l2,ibl7s,num8)
c
c     ----- calculate mo populations -----
c
      call mopop(q(i10),q(i40),q(i40+l3),iky,zzzzzz,l1,l3)
      if (zruntp.ne.zskip) then
c
c     ----- read in density matrices -----
c     -da- at q(i40)
c     -db- at q(i41)
c
         call denhf(zzzzzz,q(i40),q(i41),l2,ibl7pa,ibl7pb)
         ipass = 1
         iden = i40
         igoc = i50
         igac = i60
         iblock = ibl7o(3)
         iloc = i70
         ilac = i80
         olast = .false.
         zelec = zalpha
         if (ohfcas) zelec = ztotal
c
c     ----- do a mulliken population analysis ----
c           calculate overlap population
c
         call vmul(q(iden),1,q(i10),1,q(iden),1,l2)
 20      write (iwr,6010) zelec
         if (out) call prtril(q(iden),l1)
c
c     ----- calculate total gross population in ao*s ----
c
         call grossc(q(iden),q(igoc),l1)
         call rdedx(q(iloc),l1,iblock,num8)
c     combine alpha and beta lowdin populations if all electrons
         if (olast) then
            call vadd(q(iloc),1,q(i71),1,q(iloc),1,l1)
         end if
         write (iwr,6020)
         do 30 i = 1 , l1
            write (iwr,6030) i , zbflab(i) , q(i-1+igoc) , q(i-1+iloc)
 30      continue
c
c --- punchfile
c
         call blkao(q(igoc),q(iloc),zelec)

c
c     ----- compress from orbitals to atoms -----
c
         call atpop(q(iden),q(i30),nat)
         write (iwr,6040)
         if (out) call prtril(q(i30),nat)
c
c     ----- compress from orbitals to atoms for lowdin analysis -----
c
         call latpop(q(iloc),q(ilac),nat)
c
c     ----- calculate total gross population on atoms -----
c
         call grossc(q(i30),q(igac),nat)
         write (iwr,6050)
         do 40 i = 1 , nat
            write (iwr,6080) i , zaname(i) , czan(i) , q(i-1+igac) ,
     +                      q(i-1+ilac)
 40      continue

c
c --- punchfile
c
         call blkat(q(igac),q(ilac),zelec)
         if (.not.ohfcas) then
            if (.not.olast) then
               if (ipass.eq.2) then
c
c     ----- calculate orbital and atomic spin densities -----
c
                  call vsub(q(i50),1,q(i51),1,q(i50),1,l1)
                  i70m1 = i70 - 1
                  i71m1 = i71 - 1
                  do 50 i = 1 , l1
                     dum = q(i+i70m1) + q(i+i71m1)
                     q(i+i71m1) = q(i+i70m1) - q(i+i71m1)
                     q(i+i70m1) = dum
 50               continue
                  call vsub(q(i60),1,q(i61),1,q(i60),1,l1)
                  i80m1 = i80 - 1
                  i81m1 = i81 - 1
                  do 60 i = 1 , nat
                     dum = q(i+i80m1) + q(i+i81m1)
                     q(i+i81m1) = q(i+i80m1) - q(i+i81m1)
                     q(i+i80m1) = dum
 60               continue
                  write (iwr,6060)
                  do 70 i = 1 , num
                     write (iwr,6030) i , zbflab(i) , q(i-1+i50) ,
     +                               q(i-1+i71)
 70               continue
                  write (iwr,6070)
                  do 80 i = 1 , nat
                  write (iwr,6080) i , zaname(i) , czan(i) , q(i-1+i60)
     +                               , q(i-1+i81)
 80               continue
                  call blkat(q(i60),q(i81),zspin)
c
c     ----- do all electrons -----
c
                  ipass = 1
                  olast = .true.
                  zelec = ztotal
                  call vadd(q(i40),1,q(i41),1,q(i40),1,l2)
                  iden = i40
                  igoc = i50
                  igac = i60
                  iloc = i70
                  ilac = i80
c     restore the alpha lowdin orbital populations
c
                  call rdedx(q(i71),l1,ibl7o(3),num8)
               else
                  ipass = 2
                  iden = i41
                  iblock = ibl7o(6)
                  igoc = i51
                  igac = i61
                  iloc = i71
                  ilac = i81
                  zelec = zbeta
                   call vmul(q(iden),1,q(i10),1,q(iden),1,l2)
               end if
               go to 20
            end if
         end if
      end if
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,'anala.m','mulken','i10')

      return
 6010 format (//1x,80('=')//10x,39('-')/10x,
     +        'mulliken and lowdin population ','analyses',10x,a8,
     +        ' electrons'/10x,39('-'))
 6020 format (/10x,'----- total gross population in aos ------'/)
 6030 format (10x,i5,2x,a10,2f12.5)
 6040 format (/10x,'----- condensed to atoms -----'/)
 6050 format (/10x,'----- total gross population on atoms ----'/)
 6060 format (/10x,'----- aos spin population ------'/)
 6070 format (/10x,'----- atomic spin population -----'/)
 6080 format (10x,i5,2x,a8,2x,f6.1,2f12.5)
 6090 format (' core assignement'/' i10, i30, i40, i50,',
     +        ' i60, i70, i41, i51, i61, i71 = '/10i8/' last = ',i8)
c
      end
**==mylowd.f
      subroutine mylowd(q)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *8 obnam
      character*3 char3i,symbol(14)
      character*6 symtex(14)
      character*4 bfnam(35)
c...dimension 100 eventually has to be increased if the value 
c   of nmg is altered
      character*100 line
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
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,
     +             newbas,ncol,ivalue,ioccup,ispa,inx(maxorb),
     +             ipoint(maxorb),isit(maxorb),itogrp(maxat),
     +             ilabel(maxorb),iok
      common/craypk/mmmm(65),isymao(maxorb),isymmo(maxorb)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension q(*)
      dimension obnam(29),norb(8)
      character*8 tlabel
      common/junkc/ zlabel(maxorb),tlabel(maxorb)
c...carsten1
      data itex / 29 /
      data ev / 27.21165d0 /
      data symbol/'a''','a"','a1','b1','b2','a2',
     +            'ag','b1u','b2u','b3g','b3u','b2g','b1g','au'/
      data symtex/'a''','a''''','a_1','b_1','b_2','a_2','a_g',
     +            'b_{1u}','b_{2u}','b_{3g}','b_{3u}','b_{2g}',
     +            'b_{1g}','a_u'/
      data bfnam  /'s','x','y','z',
     +            'xx','yy','zz','xy','xz','yz',
     +            'xxx','yyy','zzz','xxy','xxz','xyy','yyz','xzz',
     +            'yzz','xyz',
     +            'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz',
     +            'zzzx','zzzy','xxyy','xxzz','yyzz','xxyz','yyxz',
     +            'zzxy' /
c...
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c     allocate memory
      nreq = l3 * 3 + l2 + 2*l1
      i500 = igmem_alloc_inf(nreq,'anala.m','mylowd',
     1                       'i500',IGMEM_DEBUG )
      i501 = i500 + l2 
      i502 = i501 + l3
      i503 = i502 + l3
      i504 = i503 + l3
      i505 = i504 + l1
      nirr = mmmm(1)
      nirs = nirr - 2
      if (nirr.eq.1) nirs = 0
c
      call setsto(8,0,norb)
      do 300 i=1,num
      k = isymmo(i)
      l = norb(k) + 1
      ipoint(i) = l
300   norb(k) = l
c...retrieve information from symmetry assignment
      call resyas
c...read in overlap matrix
      call rdedx(q(i500),l2,ibl7s,num8)
c...square overlap matrix
      call square(q(i501),q(i500),num,num)
c...diagonalize overlap matrix
      call tred2(num,num,q(i501),q(i504),q(i505),q(i502))
      call tql2(num,num,q(i504),q(i505),q(i502),ierr)
      if (ierr.ne.0) call caserr2(
     + 'problems with calculation of s**(+1/2)')
c...compute s**(+1/2)
      call vsqrt(q(i504),1,q(i504),1,num)
      i502x = i502
      i503x = i503
      i504x = i504
c...here do v*d**(+1/2)
      do 500 i = 1 , num
      call vsmul(q(i502x),1,q(i504x),q(i503x),1,num)
         i502x = i502x + num
         i503x = i503x + num
         i504x = i504x + 1
500   continue
c...here do [v*d**(+1/2)]*v(t)
      call mxmaa(q(i503),1,num,q(i502),num,1,q(i501),1,
     +                     num,num,num,num)
c...read in coefficient matrix
      call rdedx (q(i502),l3,ibl7qa,num8)
c...compute s**(+1/2) * c
      call mxmaa(q(i501),1,num,q(i502),1,num,q(i503),1,
     +                     num,num,num,num)
c...some preparation for output
      if (ngroup .gt. 0) then
         n = 0
c...find out where the AOs are centered
         do 100 ii = 1 , nshell
            iat = katom(ii)
            mini = kmin(ii)
            maxi = kmax(ii)
            do 110 i = mini , maxi
               n = n + 1
               isit(n) = iat
110         continue
100      continue
c...find out to which group the individual atoms belong
         do 120 jgroup = 1 , ngroup
            do 130 ia = 1 , mpergr
               if (iatom(ia,jgroup) .eq. 0) go to 120
               itogrp(iatom(ia,jgroup)) = jgroup
130         continue
120      continue
c...find out if d,f or g functions are present
         imaxbf = 0
         do 140 m = 1 , num
            i = locatc(bfnam,35,zbflab(m)(7:10))
            if (i .eq. 0) call caserr2(
     +      'unrecognized basis function label in mylowd')
            imaxbf = max(i,imaxbf)
140      continue
      end if
c...output of results
      iloout = ilopri / 10
      ilotex = ilopri - 10*iloout
      if (iloout .gt. 3) iloout = 3
      if (ilotex .gt. 2) ilotex = 2
      open (itex,file='lowdin.tex',status='unknown')
      pageth = 50.0d0
      thresh = 0.8d0
c...header of tex file
      write (itex,1010)
      if (ilotex .eq. 2) then
         write (itex,1020) thresh
      else
         if (ngroup .eq. 0) then
            write (itex,1030)
         else
            write (itex,1040) ('r|',jgroup=1,ngroup)
            write (itex,1041) ('&',jgroup=1,ngroup)
            write (itex,1042) 
     +     (' & \\multicolumn{1}{c|}{',jgroup,'}',jgroup=1,ngroup)
            write (itex,1043) ('&',jgroup=1,ngroup)
            write (itex,1044)
          end if
      end if
      rlines = 3.6d0
c...read in orbital energies
      call secget(mouta,3,iblko)
      call rdchr(obnam,29,iblko,idaf)
      call reads(evalue,mach(8),idaf)
c...print lowdin orbitals if desired
      if (iloout .eq. 3) then
         write (iwr,10)
         call prev(q(i503),evalue,num,num,num)
      end if
c...print summary of lowdin orbitals
      if (iloout .eq. 1) then
         write (iwr,700)
      else         
         write (iwr,710)
      end if
      i503x = i503 - 1
      write (iwr,715)
      n = 1
c...loop over orbitals
10000 if (n .gt. num) go to 900
      scont = 0.0d0
      pcont = 0.0d0
      dcont = 0.0d0
      fcont = 0.0d0
      gcont = 0.0d0
      if (ngroup .gt. 0) then
         call vclr(sco,1,ngroup)
         call vclr(pco,1,ngroup)
         call vclr(dco,1,ngroup)
         call vclr(fco,1,ngroup)
         call vclr(gco,1,ngroup)
      end if
c...evaluate s,p,d,f,g contributions
      do 718 m = 1 , num
         i = locatc(bfnam,35,zbflab(m)(7:10))
         if (i .eq. 0) call caserr2(
     +      'unrecognized basis function label in mylowd')
         c2 = q(i503x+m)**2
         if (i .eq. 1) then
            scont = scont + c2
            if (ngroup .gt. 0) sco(itogrp(isit(m)))
     +                        = sco(itogrp(isit(m))) + c2
         else if (i .le. 4) then 
            pcont = pcont + c2
            if (ngroup .gt. 0) pco(itogrp(isit(m)))
     +                        = pco(itogrp(isit(m))) + c2
         else if (i .le. 10) then
            dcont = dcont + c2
            if (ngroup .gt. 0) dco(itogrp(isit(m)))
     +                        = dco(itogrp(isit(m))) + c2
         else if (i .le. 20) then
            fcont = fcont + c2
            if (ngroup .gt. 0) fco(itogrp(isit(m)))
     +                        = fco(itogrp(isit(m))) + c2
         else
            gcont = gcont + c2
            if (ngroup .gt. 0) gco(itogrp(isit(m)))
     +                        = gco(itogrp(isit(m))) + c2
         end if
718   continue
      scont = scont * 100.0d0
      pcont = pcont * 100.0d0
      dcont = dcont * 100.0d0
      fcont = fcont * 100.0d0
      gcont = gcont * 100.0d0
c...determine multiplicity
      imult = 1
      if (iok .ne. 2) then
         if (zlabel(n)(1:1) .eq. 'e') imult = 2
         if (zlabel(n)(1:1) .eq. 't') imult = 3
      end if
      if (iloout.eq.1 .and. ilotex.eq.1) go to 890
c...sort coefficients
      do 720 m = 1 , num
         inx(m) = m
720   continue
      do 740 i = 1 , num
         k = i
         if (num .gt. 1) then
            cmax = dabs(q(i503x+i))
            do 730 j = i+1 , num
               if (dabs(q(i503x+j)) .gt. cmax) then
                  k = j
                  cmax = dabs(q(i503x+j))
               end if
730         continue
         end if
         t = q(i503x+i)
         q(i503x+i) = q(i503x+k)
         q(i503x+k) = t
         it = inx(i)
         inx(i) = inx(k)
         inx(k) = it
740   continue
      total = 0.0d0
      do 750 m = 1 , num
         total = total + q(i503x+m)**2
         if (total .ge. thresh) go to 760
750   continue
760   continue
c...output of most important contributions
      n2 = min(m,4)
      if (iloout .eq. 1) go to 7705
c...standard output, first line
      if (iok .ne. 2) then
         write (iwr,769) n,ilabel(n),zlabel(n)(1:4),
     +        evalue(n)*ev,(q(i503x+i),zbflab(inx(i)),inx(i),i=1,n2)
      else
         write (iwr,770) n,ipoint(n),symbol(nirs+isymmo(n)),
     +        evalue(n)*ev,(q(i503x+i),zbflab(inx(i)),inx(i),i=1,n2)
      end if
7705  if (ilotex .eq. 1) go to 7745
c...tex file, first line
      if (iok .ne. 2) then
       write (itex,773) n,ilabel(n),tlabel(n),evalue(n)*ev
      else
        write (itex,783) n,ipoint(n),symtex(nirs+isymmo(n)),evalue(n)*ev
      end if
      do 775 i = 1 , n2
         write (itex,774) q(i503x+i),zbflab(inx(i))
775   continue
      do 777 i = 1 , 4-n2
         write (itex,776)
777   continue
      write (itex,778)
      rlines = rlines + 1.0d0
7745  n1 = 5
      n2 = 8
780   if (m .lt. n1) go to 890
      n2 = min(n2,m)
      if (iloout .eq. 1) go to 7905
c...standard output, subsequent lines
      write (iwr,790) (q(i503x+i),zbflab(inx(i)),inx(i),i=n1,n2)
7905  if (ilotex .eq. 1) go to 800
c...tex file, subsequent lines
      write (itex,791)
      do 793 i = n1 , n2
         write (itex,792) q(i503x+i),zbflab(inx(i))
793   continue
      do 795 i = 1 , 4-(n2-n1+1)
         write (itex,776)
795   continue
      write (itex,778)
      rlines = rlines + 1.0d0
800   n1 = n1 + 4
      n2 = n2 + 4
      go to 780
890   i503x = i503x + imult*num
c...output of s,p,d,f contributions
c
c...standard output
      if (ngroup .eq. 0) then
         if (iloout .eq. 1) then
            if (iok .ne. 2) then
               write (iwr,2001) n,ilabel(n),zlabel(n)(1:4),
     +                         evalue(n)*ev,
     +                         scont,pcont,dcont,fcont,gcont
            else
               write (iwr,2002) n,ipoint(n),symbol(nirs+isymmo(n)),
     +                         evalue(n)*ev,
     +                         scont,pcont,dcont,fcont,gcont
            end if
         else
            write (iwr,2200) scont,pcont,dcont,fcont,gcont
         end if
      else
         iline = 0
         do 2205 jgroup = 1 , ngroup
            line(iline+1:iline+4) = '|   '
            line(iline+5:iline+7) = char3i(jgroup)
            line(iline+8:iline+8) = ' '
            iline = iline + 8
2205     continue
         line(iline+1:iline+9) = '| total |'
         iline = iline + 9
         if (iloout .eq. 1) then
            if (iok .ne. 2) then
               write (iwr,2210) n,ilabel(n),zlabel(n)(1:4),
     +                    evalue(n)*ev,(line(jline:jline),jline=1,iline)
            else
               write (iwr,2220) n,ipoint(n),symbol(nirs+isymmo(n)),
     +                    evalue(n)*ev,(line(jline:jline),jline=1,iline)
            end if
         else
            write (iwr,2230) (line(jline:jline),jline=1,iline)
         end if
         write (iwr,2240) ('-',i=1,ngroup*8+14)
c...94 should be changed to nmg*8+14 if the value of nmg is enlarged
         write (iwr,2206) 's',(sco(jgroup)*100.0d0,jgroup=1,ngroup),
     +             scont
         write (iwr,2240) ('-',i=1,ngroup*8+14)
         if (imaxbf .gt. 1) then
            write (iwr,2206) 'p',(pco(jgroup)*100.0d0,jgroup=1,ngroup)
     +            ,pcont
            write (iwr,2240) ('-',i=1,ngroup*8+14)
         end if
         if (imaxbf .gt. 4) then
            write (iwr,2206) 'd',(dco(jgroup)*100.0d0,jgroup=1,ngroup)
     +            ,dcont
            write (iwr,2240) ('-',i=1,ngroup*8+14)
         end if
         if (imaxbf .gt. 10) then
            write (iwr,2206) 'f',(fco(jgroup)*100.0d0,jgroup=1,ngroup)
     +            ,fcont
            write (iwr,2240) ('-',i=1,ngroup*8+14)
         end if
         if (imaxbf .gt. 15) then
            write (iwr,2206) 'g',(gco(jgroup)*100.0d0,jgroup=1,ngroup)
     +            ,gcont
            write (iwr,2240) ('-',i=1,ngroup*8+14)
         end if
         write (iwr,2207)
     +      ((sco(jgroup)+pco(jgroup)+dco(jgroup)+fco(jgroup)+
     +        gco(jgroup))*100.0d0,
     +       jgroup=1,ngroup),100.0d0
         write (iwr,2240) ('-',i=1,ngroup*8+14)
      end if
c...tex file
      if (ngroup .eq. 0) then
         if (ilotex .eq. 1) then
            if (iok .ne. 2) then
               write (itex,3001) n,ilabel(n),tlabel(n),
     +                           evalue(n)*ev,
     +                           scont,pcont,dcont,fcont,gcont
            else
               write (itex,3002) n,ipoint(n),symtex(nirs+isymmo(n)),
     +                           evalue(n)*ev,
     +                           scont,pcont,dcont,fcont,gcont
            end if
         else
            write (itex,3200) scont,pcont,dcont,fcont,gcont
         end if
         rlines = rlines + 1.0d0
      else
         if (ilotex .eq. 1) then
            if (iok .ne. 2) then
               write (itex,3210) n,ilabel(n),tlabel(n),evalue(n)*ev,
     +                 ('& $',sco(jgroup)*100.0d0,'$',jgroup=1,ngroup),
     +                 '& $',scont,'$'
            else
               write (itex,3220) n,ipoint(n),symtex(nirs+isymmo(n)),
     +                 evalue(n)*ev,
     +                 ('& $',sco(jgroup)*100.0d0,'$',jgroup=1,ngroup),
     +                 '& $',scont,'$'
            end if
            write (itex,778)
            rlines = rlines + 1.0d0
            if (imaxbf .gt. 1) then
               write (itex,3206) 'p',('& $',pco(jgroup)*100.0d0,'$'
     +               ,jgroup=1,ngroup),
     +                                '& $',pcont,'$'
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 4) then
               write (itex,3206) 'd',('& $',dco(jgroup)*100.0d0,'$'
     +               ,jgroup=1,ngroup),
     +                                '& $',dcont,'$'
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 10) then
               write (itex,3206) 'f',('& $',fco(jgroup)*100.0d0,'$'
     +               ,jgroup=1,ngroup),
     +                                '& $',fcont,'$'
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 15) then
               write (itex,3206) 'g',('& $',gco(jgroup)*100.0d0,'$'
     +               ,jgroup=1,ngroup),
     +                                '& $',fcont,'$'
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            write (itex,3205) ngroup+5
            write (itex,3206) 'all',
     +         ('& $',(sco(jgroup)+pco(jgroup)+dco(jgroup)
     +            +fco(jgroup)+gco(jgroup))*100.0d0,'$',
     +          jgroup=1,ngroup),'& $',100.0d0,'$'
            write (itex,778)
            rlines = rlines + 1.3d0
         else
            write (itex,3230) ('r|',jgroup=1,ngroup)
            write (itex,3231) (jgroup,jgroup=1,ngroup)
            write (itex,3232)
            write (itex,3209) 's',('&',sco(jgroup)*100.0d0,
     +                         jgroup=1,ngroup),'&',scont
            write (itex,778)
            rlines = rlines + 2.3d0
            if (imaxbf .gt. 1) then
               write (itex,3209) 'p',('&',pco(jgroup)*100.0d0,
     +                         jgroup=1,ngroup),'&',pcont
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 4) then
               write (itex,3209) 'd',('&',dco(jgroup)*100.0d0,
     +                         jgroup=1,ngroup),'&',dcont
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 10) then
               write (itex,3209) 'f',('&',fco(jgroup)*100.0d0,
     +                         jgroup=1,ngroup),'&',fcont
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            if (imaxbf .gt. 15) then
               write (itex,3209) 'g',('&',gco(jgroup)*100.0d0,
     +                         jgroup=1,ngroup),'&',gcont
               write (itex,778)
               rlines = rlines + 1.0d0
            end if
            write (itex,891)
            write (itex,3209) 'all',
     +         ('&',(sco(jgroup)+pco(jgroup)+dco(jgroup)
     +          +fco(jgroup)+gco(jgroup))*100.0d0,
     +          jgroup=1,ngroup),'&',100.0d0
            write (itex,778)
            write (itex,3240)
            rlines = rlines + 1.6d0
         end if
      end if
      write (iwr,715)
      write (itex,891)
      rlines = rlines + 0.3d0
      if (rlines.ge.pageth .and. n.lt.num) then
         write (itex,2140)
         if (ilotex .eq. 2) then
            write (itex,1020) thresh
         else
            if (ngroup .eq. 0) then
               write (itex,1030)
            else
               write (itex,1040) ('r|',jgroup=1,ngroup)
               write (itex,1041) ('&',jgroup=1,ngroup)
               write (itex,1042) (jgroup,jgroup=1,ngroup)
               write (itex,1043) ('&',jgroup=1,ngroup)
               write (itex,1044)
            end if
         end if
         rlines = 3.6d0
      end if
      n = n + imult
      go to 10000
900   continue
c...close tex file
      write (itex,3010)
      close (itex)
c...reset core memory
      call gmem_free_inf(i500,'anala.m','mylowd','i500')
      return
10    format (//1x,120('-')//
     +     /40x,'*******************************************'
     +     /40x,'                                      1/2  '
     +     /40x,'    lowdin orthonormalization:   G = S   *C'
     +     /40x,'*******************************************'
     +    //)
700   format (//1x,130('-')//
     +      40x,'*********************************'/
     +      40x,' Analysis of orbital characters  '/
     +      40x,'*********************************'//
     +      4x,'MO',5x,'energy (eV)')
710   format (//1x,130('-')//
     +         40x,'*********************************'/
     +         40x,'Summary of most important weights'/
     +         40x,'& analysis of orbital characters '/
     +         40x,'*********************************'//
     +         4x,'MO',5x,
     +            'energy (eV)   Most important AO contributions')
715   format (115('-'))
769   format (i3,1x,i3,a4,f9.2,1x,4(f9.4,' (',a10,'/',i3,')'))
770   format (i3,1x,i3,a3,1x,f9.2,1x,4(f9.4,' (',a10,'/',i3,')'))
773   format ('$',i3,'$ & $',i3,a8,'$ & $',f9.2,'$')
783   format ('$',i3,'$ & $',i3,a6,'$ & $',f9.2,'$')
774   format (' & $',f9.4,'$ & (',a10,')')
776   format ('&&')
790   format (21x,4(f9.4,' (',a10,'/',i3,')'))
791   format ('&&')
792   format ('& $',f9.4,'$ & (',a10,')')
2001  format (i3,1x,i3,a4,f9.2,4x,
     +        'orbital character: ',f5.1,' % s, ',
     +        f5.1,' % p, ',f5.1,' % d, ',f5.1,' % f ',
     +        f5.1,' % g')
2002  format (i3,1x,i3,a3,1x,f9.2,4x,
     +        'orbital character: ',f5.1,' % s, ',
     +        f5.1,' % p, ',f5.1,' % d, ',f5.1,' % f ',
     +        f5.1,' % g')
2200  format (24x,'orbital character: ',f5.1,' % s, ',
     +     f5.1,' % p, ',f5.1,' % d, ',f5.1,' % f, ',f5.1,' % g')
2210  format (i3,1x,i3,a4,f9.2,17x,100a1)
2220  format (i3,1x,i3,a3,1x,f9.2,17x,100a1)
2230  format (37x,100a1)
2240  format (32x,94a1)
c...format 11(...) should be changed to <nmg+1>(...) 
c...if the value of nmg is enlarged
2206  format (34x,a1,'  |',11(f6.1,1x,'|'))
2207  format (33x,'all |',11(f6.1,1x,'|'))
3210  format (i3,' & $',i3,a8,' $ & $',f9.2,'$ & s ',11(a3,f5.1,a1))
3220  format (i3,' & $',i3,a6,' $ & $',f9.2,'$ & s ',11(a3,f5.1,a1))
3206  format ('&&&',a3,11(a3,f5.1,a1))
3231  format ('r|}'/'&',10(i2,' & '))
3209  format (a3,11(a1,f5.1))
1010  format ('\\documentstyle[11pt]{article}'/
     +        '\\pagestyle{plain}'/
     +        '\\topmargin=-1.8cm'/
     +        '\\oddsidemargin=-1.5cm'/
     +        '\\evensidemargin=-0.5cm'/
     +        '\\textwidth=17cm'/
     +        '\\textheight=26cm'/
     +        '\\parindent=0cm'/
     +        '\\begin{document}')
1020  format('\\begin{tabular}{|r|l|r|rlrlrlrl|} \\hline ',
     + '\\multicolumn{2}{|c|}{}&&& &&& &&& \\\\'/
     + '\\multicolumn{2}{|c|}{MO} & \\multicolumn{1}{c|}{energy (eV)}',
     + ' & \\multicolumn{8}{l|}{Most important AO contributions',
     + ' ($\\sum\\limits_i c^2_i\\,\\geq\\,',f4.2,'$)}',
     + '\\\\ \\multicolumn{2}{|c|}{}&&& &&& &&& \\\\ \\hline')
1030  format ('\\begin{tabular}{|r|l|r|r|r|r|r|} \\hline ',
     +  '\\multicolumn{2}{|c|}{} &&&&& \\\\'/
     +  '\\multicolumn{2}{|c|}{MO} & \\multicolumn{1}{c|}{energy (eV)}',
     +  ' & \\multicolumn{1}{c|}{s} & \\multicolumn{1}{c|}{p}',
     +  ' & \\multicolumn{1}{c|}{d} & \\multicolumn{1}{c|}{f} \\\\',
     +  '\\multicolumn{2}{|c|}{}&&&&& \\\\ \\hline')
c...format 10a2 should be changed to <nmg>a2 if the value of 
c...nmg is enlarged
1040  format ('\\begin{tabular}{|r|l|r|r|',10a2)
c...format 10a1 should be changed to <nmg>a1 if the value of 
c...nmg is enlarged
1041  format ('r|} \\hline \\multicolumn{2}{|c|}{} &&&',10a1)
c...format 10(...) should be changed to <nmg>(...) if the value of 
c...nmg is enlarged
1042  format (' \\\\'/
     +'\\multicolumn{2}{|c|}{MO} & \\multicolumn{1}{c|}{energy (eV)} & '
     +,10(a23,i2,a1))
c...format 10a1 should be changed to <nmg>a1 if the value of 
c...nmg is enlarged
1043  format ('& total \\\\'/'\\multicolumn{2}{|c|}{} &&& ',10a1)
1044  format (' \\\\ \\hline')
778   format ('\\\\')
3001  format (i3,
     +     ' & $',i3,a8,'$ & $',f9.2,'$',4(' & ',f5.1,'\\,\\%'),'\\\\')
3002  format (i3,
     +     ' & $',i3,a6,'$ & $',f9.2,'$',4(' & ',f5.1,'\\,\\%'),'\\\\')
3200  format (
     + '&&& \\multicolumn{8}{l|}{orbital character: ',f5.1,' \\% s, ',
     +        f5.1,' \\% p, ',f5.1,' \\% d, ',f5.1,' \\% f} \\\\')
3205  format ('\\cline{4-',i2,'}')
3230  format ('&&& \\multicolumn{8}{l|}{\\begin{tabular}{|r|',10a2)
3232  format (' total \\\\ \\hline')
3240  format ('\\hline'/'\\end{tabular}} \\\\ \\hline')
891   format ('\\hline')
2140  format('\\end{tabular}'/'\\newpage')
3010  format ('\\end{tabular}'/'\\end{document}')
      end
**==resyas.f
      subroutine resyas
c
c--------------------------------------------------------------------
c   Restores information from symmetry assignment
c   and converts orbital labels to TeX format
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
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
      parameter(isymtp=99)
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
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,
     +             newbas,ncol,ivalue,ioccup,ispa,inx(maxorb),
     +             ipoint(maxorb),isit(maxorb),itogrp(maxat),
     +             ilabel(maxorb),iok
      character*8 zlabel,tlabel
      common/junkc/zlabel(maxorb),tlabel(maxorb)
c
      isymsc = isect(499)
      call qsector(isymsc,m,ityp,nblk,'get')
      if (m .eq. 0) then
         write (iwr,1) isymsc
1        format (//'WARNING: Retrieval of information about symmetry ',
     +       'assignment failed : Section ',i3,' is not present.'//)
         iok = 2
      else
         call secget (isymsc,isymtp,iblnum)
         call readi (iok,1,iblnum,idaf)
         call readis (ilabel,num,idaf)
         call rdchrs (zlabel,num,idaf)
         if (iok .ne. 2) then ! convert zlabel into tex format
            do 250 m = 1 , num
               tlabel(m)(1:1) = zlabel(m)(1:1)
               if (zlabel(m)(2:2) .eq. ' ') go to 250
               if (zlabel(m)(2:2) .eq. '''') then
                  tlabel(m) = zlabel(m)
               else
                  tlabel(m)(2:3) = '_{'
                  if (zlabel(m)(2:2).ne.' ' .and. 
     +                   zlabel(m)(2:2).ne.'''') then
                     tlabel(m)(4:4) = zlabel(m)(2:2)
                     i = 3
                  end if
                  if (zlabel(m)(3:3).ne.' ' .and. 
     +                   zlabel(m)(3:3).ne.'''') then
                     tlabel(m)(5:5) = zlabel(m)(3:3)
                     i = 4
                  end if
                  tlabel(m)(i+2:i+2) = '}'
                  tlabel(m)(i+3:i+3) = zlabel(m)(i:i)
                  tlabel(m)(i+4:i+4) = zlabel(m)(i+1:i+1)
               end if
250         continue
         end if
      end if
      return
      end
**==poporb.f
      subroutine poporb(zzzzzz,vecin,tran,vecout,iord,iir,
     *                  map,rij,qpix,qpjx,norot,sao,
     *                  l1,l2,l3,n1,n2,natoms,odef,nprint)
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
c
      dimension vecin(*),tran(n1,n1),sao(l2),norot(n1,n1),
     +          vecout(*),rij(n2,natoms),
     +          iir(l1),map(l1),iord(l1),qpix(l1),qpjx(l1),
     +          ytype(2)
c
c
      common /junk / limlow(maxat),limsup(maxat)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common /miscop/ nacta,nactb,moout(maxorb),mooutb(maxorb)
     * , eiga(maxorb),frocca(maxorb),eigb(maxorb),froccb(maxorb)
c
      character*10 charwall
      parameter (dzero=0.0d+00, done=1.0d+00, two=2.0d+00, four=4.0d+00,
     *           tenm3=1.0d-03, tenm8=1.0d-08, tenm10=1.0d-10)
c
      data zuhf/'uhf'/
      data ytype/'-a-','-b-'/
      data m1/1/
c
      ind(i,j) = ((max(i,j)*(max(i,j)-1))/2+min(i,j))
c
c     poporb calculates the localized molecular orbitals by the
c     population method of j. pipek and p. g. mezey,
c     j. chem. phys. 90, 4916 (1989),
c
      cpu = cpulft(1)
      write (iwr,6130) cpu,charwall()
      nsav = mouta
      mi = 0
      norb = nacta
      iblock = ibl7qa
c
      nbasis= l1
c     nocc  = nbasis * norb
      maxit = max(maxcyc,1500)
      cvgloc = 10.0d0**(-nconv)
c
      call aolim
      do 20 k=1,nbasis
      do 20 i=1,natoms
      if(limlow(i).le.k.and.limsup(i).ge.k) map(k)=i
   20 continue
c
2000  mi = mi + 1
      if (norb.eq.0) go to 330
      call gnorot(zzzzzz,mi,norot,n1,norb,odef)
      nredo = 0
      call rdedx(vecin,l3,iblock,num8)
c
      nredo = 0
  110 continue
      nredo=nredo+1
c
c          construct initial atomic populations
c
      call vclr(rij,1,n2*natoms)
      call rdedx(sao,l2,ibl7x,num8)
      ij = 0
      do 280 i = 1,norb
      iim = ilifq(moout(i))
         do 280 j = 1,i
            ij = ij+1
            jjm = ilifq(moout(j))
            do 260 k = 1,nbasis
               kk = map(k)
               do 260 l = 1,nbasis
                  ll = map(l)
                  kl = ind(k,l)
                  sum = vecin(k+iim)*vecin(l+jjm)*sao(kl)/two
                  rij(ij,kk) = rij(ij,kk) + sum
                  rij(ij,ll) = rij(ij,ll) + sum
  260       continue
  280 continue
c
c          compute initial localization sum
c
      sumrr = dzero
      do 320 i=1,norb
         if(norot(i,i).eq.1) go to 320
         ii = ind(i,i)
         do 310 k=1,natoms
            sumrr = sumrr + rij(ii,k)**2
  310    continue
  320 continue
      write(iwr,9010) sumrr
      call popqat(rij,vecout,n2,norb,natoms,zzzzzz,iwr)
c
c          seed the random function, initialize
c          the localization transformation, etc.
c
      donept = done
      if (nredo.eq.1) dd = rndlmo(donept,vecin,l1)
c
      call vclr(tran,1,norb*norb)
      do 340 i = 1,norb
         tran(i,i) = done
  340 continue
      iter = 0
      shift = atan(donept)
c
c          begin localization cycles
c
  360 continue
      change = dzero
      iter = iter+1
      do 380 i = 1,norb
         iir(i) = i
  380 continue
      nnn = norb
      do 400 i = 1,norb
         dd = rndlmo(change,vecin,l1)
         iii = int(dd*dble(nnn)+done)
         iord(i) = iir(iii)
         iir(iii) = iir(nnn)
         nnn = nnn-1
  400 continue
c
c        for each pair of orbitals a two dimensional unitary
c        transformation is performed. the transformation is
c
c           psi'(i) =  cos(t)*psi(i) + sin(t)*psi(j)  and
c           psi'(j) = -sin(t)*psi(i) + cos(t)*psi(j).
c
c        localization requires that t be such as to maximize
c        the sum of the squares of the atomic populations.
c
      do 920 iii = 1,norb
      i = iord(iii)
      if(norot(i,i).eq.1) go to 920
      ii = ind(i,i)
      jm = 1
      rm = dzero
      tm = dzero
      sm = dzero
      cm = done
      do 580 j = 1,norb
      if(i.eq.j) go to 580
      if(norot(i,j).eq.1) go to 580
      ij = ind(i,j)
      jj = ind(j,j)
      t = dzero
      tx = dzero
      do 480 kk = 1,natoms
         t = t + four*rij(ij,kk)**2 - rij(ii,kk)**2 - rij(jj,kk)**2
     *         + two*rij(ii,kk)*rij(jj,kk)
         tx = tx + rij(ij,kk)*(rij(jj,kk) - rij(ii,kk))
  480 continue
      if ((dabs(t) .le. tenm10) .and. (dabs(tx) .le. tenm10)) go to 580
      tx = four*tx
      t = datan2(tx,t)/four
      sign = done
      if (t .gt. dzero) sign = -done
      t = t+sign*shift
      itim = 0
  500 itim = itim+1
      s = dsin(t)
      cc = dcos(t)
      rin = dzero
      do 520 kk = 1,natoms
         qpi = cc*cc*rij(ii,kk)+s*s*rij(jj,kk)+two*cc*s*rij(ij,kk)
         qpj = cc*cc*rij(jj,kk)+s*s*rij(ii,kk)-two*cc*s*rij(ij,kk)
         rin = rin+qpi*qpi+qpj*qpj-rij(ii,kk)**2-rij(jj,kk)**2
  520 continue
      ttest = dabs(t)-shift
      if ((dabs(t) .le. tenm8) .or. (dabs(ttest) .le. tenm8)) go to 560
      if (rin .ge. -tenm8) go to 560
      if (itim .le. 1) go to 540
         write (iwr,9020) i,j
         write (iwr,9030) t,s,cc,rin
      return
c
  540 sign = done
      if (t .gt. dzero) sign = -done
      t = t+shift*sign
      go to 500
c
  560 if (rin .le. rm) go to 580
      rm = rin
      tm = t
      sm = s
      cm = cc
      jm = j
  580 continue
c
      rin = rm
      t = tm
      s = sm
      cc = cm
      j = jm
      ij = ind(i,j)
      jj = ind(j,j)
      if(norot(i,j).eq.1) go to 920
c
c        accumulate the 2x2 rotation
c
      change = change+t*t
      call drot(norb,tran(1,i),1,tran(1,j),1,cc,s)
c
c        update the atomic populations
c
      do 880 kk = 1,natoms
         qpi = cc*cc*rij(ii,kk)+s*s*rij(jj,kk)+two*cc*s*rij(ij,kk)
         qpj = cc*cc*rij(jj,kk)+s*s*rij(ii,kk)-two*cc*s*rij(ij,kk)
         qpij = (cc*cc-s*s)*rij(ij,kk)+cc*s*(rij(jj,kk)-rij(ii,kk))
         do 720 k = 1,norb
            if (i.eq.k.or.j.eq.k) goto 720
            ik = ind(i,k)
            jk = ind(j,k)
            qpix(k) = cc*rij(ik,kk)+s*rij(jk,kk)
            qpjx(k) = cc*rij(jk,kk)-s*rij(ik,kk)
            rij(ik,kk) = qpix(k)
            rij(jk,kk) = qpjx(k)
  720    continue
         rin = rin+qpi+qpj-rij(ii,kk)-rij(jj,kk)
         rij(ii,kk) = qpi
         rij(jj,kk) = qpj
         rij(ij,kk) = qpij
  880 continue
  920 continue
c
c          test for convergence of localization procedure
c
      change = sqrt(two*change/(norb*(norb-1)))
      if(nprint.eq.-7)  write(iwr,9050) iter,change
      if(mod(iter,10).eq.0) write(iwr,9050) iter,change
      if(iter.lt.maxit  .and.  change.gt.tenm3*cvgloc) go to 360
      if(change.le.cvgloc) go to 1000
         if(nredo.le.2) write(iwr,9060)
         if(nredo.le.2) go to 110
            write(iwr,9070)
            go to 1010
c
c          finished with localization cycles
c
 1000 continue
      write (iwr,9080) iter,cpulft(1),charwall()
c
c        transform to final orbitals, copy virtual space
c
1010  call rdedx(vecout,l3,iblock,num8)
      write(10)norb
      do i=1,norb
         write(10)(tran(j,i),j=1,norb)
      enddo
      do 930 i = 1 , norb
         call vclr(rij(1,1),1,nbasis)
         do 940 j = 1 , norb
         iij = ilifq(moout(j))+1
         sc = tran(j,i)
         call daxpy(nbasis,sc,vecout(iij),1,rij(1,1),1)
 940     continue
         call dcopy(nbasis,rij(1,1),1,vecin(1+ilifq(i)),1)
 930  continue
c
c ----- now load up lmos into vecin array
c
      do 270 i = 1 , norb
         ii = ilifq(moout(i))
         call dcopy(nbasis,vecin(1+ilifq(i)),1,vecout(ii+1),1)
 270  continue
c
c ----- now save onto dumpfile
c
      if (mi.eq.1) call putq(zcom,ztitle,eiga,frocca,l1,l1,l1,m1,m1,
     +             vecout,nsav,imu)
      if (mi.eq.2) call putq(zcom,ztitle,eigb,froccb,l1,l1,l1,m1,m1,
     +             vecout,nsav,imu)
c ----
      write (iwr,6010) ytype(mi) , nsav
c
c          print/punch final localized orbitals
c
      write(iwr,9090)
      if (mi.eq.1) write (iwr,6080)
      if (mi.eq.2) write (iwr,6090)
      call prsql(vecout,l1,l1,l1)
c
      if(oprint(54)) then
         write(iwr,9120)
         call prsq(tran,norb,norb,n1)
      end if
c
      write (iwr,9130) norb,nbasis
c
      if (nprint.eq.7) then
       if(mi.eq.1) write(ipu,8000)
       call pusql(vecout,norb,l1,l1)
       if (zzzzzz.ne.zuhf  .or.  mi.eq.2)write(ipu,8010)
      end if
c
c          construct final atomic populations
c
      call vclr(rij,1,n2*natoms)
      call rdedx(sao,l2,ibl7x,num8)
      ij = 0
      do 1280 i = 1,norb
      iim = ilifq(moout(i))
         do 1280 j = 1,i
         jjm = ilifq(moout(j))
            ij = ij+1
            do 1260 k = 1,nbasis
               kk = map(k)
               iik = iim + k
               do 1260 l = 1,nbasis
                  ll = map(l)
                  kl = ind(k,l)
                  sum = vecout(iik)*vecout(l+jjm)*sao(kl)/two
                  rij(ij,kk) = rij(ij,kk) + sum
                  rij(ij,ll) = rij(ij,ll) + sum
 1260       continue
 1280 continue
c
c               compute final localization sum (skip frozen mo-s)
c
      sumrr = dzero
      do 1380 i=1,norb
         if(norot(i,i).eq.1) go to 1380
         ii = ind(i,i)
         do 1370 k = 1,natoms
            sumrr = sumrr + rij(ii,k)**2
 1370    continue
 1380 continue
      write(iwr,9140) sumrr
      call popqat(rij,vecin,n2,norb,natoms,zzzzzz,iwr)
 330  if (mi.ne.1 .or. nactb.eq.0) then
         cpu = cpulft(1)
         write (iwr,6020) cpu,charwall()
      else
         iblock = ibl7qb
         nsav = moutb
         norb = nactb
         do 350 i = 1 , norb
            moout(i) = mooutb(i)
 350     continue
         go to 2000
      end if
      return
c
 8000 format('pipek-mezey population localized orbitals'/' $vec')
 8010 format(' $end')
 6010 format (/1x,a3,' lmo output to section',i4,' of the dumpfile')
 6020 format (/' end of localization at ',f10.2,' seconds'
     +        ,a10,' wall'//1x,80('='))
 6130 format (/1x,
     +  'commence Pipek-Mezey localization at ',f8.2,' seconds'
     +  ,a10,' wall')
c
 9010 format(/10x,'the initial localization sum is',f14.6)
 9020 format(1x,'no rotation increases atomic populations',
     *           ' --- localization aborted'/
     *           10x,'i=',i3,5x,'j=',i3)
 9030 format(5x,8htheta = ,g20.10/5x,12hsin(theta)= ,f10.7,
     +     15h   cos(theta)= ,f10.7/5x,29htotal change to this point = ,
     +     g20.10)
 9050 format(10x,'iteration',i4,'   orbital change=',g20.10)
 9060 format (//10x,'*** localization has been unsucessful ***'
     +        //10x,'program will restart with new random number'
     +         /10x,'and rotation sequence for orbitals')
 9070 format(/10x,'+++++++++++++++++++++++++++++++++++++++++++++++'/
     +        10x,'+ localization fails --- localization aborted +'/
     +        10x,'+++++++++++++++++++++++++++++++++++++++++++++++')
 9080 format(/1x,'localization converged in',i3,
     +           ' iterations at ',f9.2,' seconds',a10,' wall')
 9090 format(/10x,'Pipek-Mezey population localized orbitals are'/)
 6080 format(50x,14('*')/50x,'alpha spin lmo'/50x,14('*')/)
 6090 format(50x,14('*')/50x,' beta spin lmo'/50x,14('*')/)
 9120 format(/10x,'the transformation from initial mo-s (rows)',
     *             ' to localized mo-s (columns) is')
 9130 format(/10x,'this localization had',i3,' m.o. s  and',i3,
     +           ' basis functions')
 9140 format(10x,'the final localization sum is',f16.6/)
      end
**==popqat.f
      subroutine popqat(rij,pop,n2,norb,natoms,zscftp,iw)
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
      dimension rij(n2,1),pop(natoms,norb)
      parameter (done=1.0d+00, two=2.0d+00)
c
c     data zrhf/'rhf'/
      data zuhf/'uhf'/
c
c     compute mulliken gross atomic populations
c
      ind(i,j) = ((max(i,j)*(max(i,j)-1))/2+min(i,j))
      docc = two
      if (zscftp.eq.zuhf) docc = done
      do 10 i=1,norb
      ii = ind(i,i)
      do 10 j=1,natoms
   10 pop(j,i) = docc*rij(ii,j)
      write(iw,1000)
 1000 format(/10x,'mulliken gross atomic populations'/
     +        10x,'*********************************'/)
      call prsq(pop,norb,natoms,natoms)
      write(iw,1001)
 1001 format(/)
      return
      end
**==randd.f
      function randd(qq,d,l1)
c...    semi-random number generator
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
      common /saveco/ u,idadap,ifirst,ijato,klato,ijshlo
      dimension d(l1,l1)
      if (qq.eq.0.0d0) then
         u = (3.14160d0+u)**5
         pp = dble(idint(u))
         u = u - pp
         randd = u
         return
      else
         n = iabs(na-num) + 1
         m = n + 5
         pp = d(n,m)*datan(1.00d0)
         u = (3.14160d0+pp)**5
         pp = dble(idint(u))
         u = u - pp
         randd = u
         return
      end if
      end
**==rndlmo.f
      function rndlmo(dxx,d,l1)
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
      dimension d(l1,l1),u(1)
c
      save u
c
      parameter (dzero=0.0d+00, done=1.0d+00)
c
      pi = acos(-done)
      if (dxx .ne. dzero) then
         n = abs(na-num)+1
         m = n+5
         dxy = d(n,m)*atan(done)
         u(1) = (pi+dxy)**5
         dxy = dble(int(u(1)))
         u(1) = u(1)-dxy
         rndlmo = u(1)
         return
      endif
c
      u(1) = (pi+u(1))**5
      dxy = dble(int(u(1)))
      u(1) = u(1)-dxy
      rndlmo = u(1)
      return
      end
**==spind.f
      subroutine spind(zzzzzz,q)
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
      common/junk/pint,qint,rint,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,
     + ss(225),dij(225),pin(32),qin(32),rin(32),
     + ijx(225), ijy(225), ijz(225)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension q(*)
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      data dzero,two /0.0d0,2.0d0/
      data conpi /3.14159265d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     *          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     *          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     *          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     *          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     *          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     *          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     *          1,16, 1, 6, 1,11,11, 1, 6, 6,
     *          1,21, 1, 6, 1,16,16, 1, 6,11,
     *          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     *          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     *          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     *          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     *          1, 1,16, 1, 6, 1, 6,11,11, 6,
     *          1, 1,21, 1, 6, 1, 6,16,16, 1,
     *         11,11, 6, 6,11/
      data zrhf/'rhf'/
      data zskip/'skippop'/
c
      if (zzzzzz.eq.zrhf) return
      if (zruntp.eq.zskip) return
      if (na.eq.nb) return
      l1 = num
      l2 = (num*(num+1))/2
c
c     ----- set pointers for partitioning of core -----
c     ----- get core memory -----
c
      i10 = igmem_alloc_inf(3*l2+nat,'anala.m','spind',
     1                      'i10',IGMEM_DEBUG )
      i20 = i10 + l2
      i30 = i20 + l2
      i40 = i30 + l2
c     last = i40 + nat
c
c     ----- get density matrices -----
c
c     -da- at x(i10)
c     -db- at x(i20)
c
      call denhf(zzzzzz,q(i10),q(i20),l2,ibl7pa,ibl7pb)
c
      write (iwr,6020)
      tol = rln10*itol
      out = nprint.eq.3
      onorm = normf.ne.1 .or. normp.ne.1
c
c     ----- loop over atoms -----
c
      do 320 iat = 1 , nat
         p0 = c(1,iat)
         q0 = c(2,iat)
         r0 = c(3,iat)
c
c     ----- ishell
c
         do 310 ii = 1 , nshell
            i = katom(ii)
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
            do 300 jj = 1 , ii
               j = katom(jj)
               pj = c(1,j)
               qj = c(2,j)
               rj = c(3,j)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
               oianj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
               ij = 0
               max = maxj
               do 30 i = mini , maxi
                  nnx = ix(i)
                  nny = iy(i)
                  nnz = iz(i)
                  if (oianj) max = i
                  do 20 j = minj , max
                     ij = ij + 1
                     ijx(ij) = nnx + jx(j)
                     ijy(ij) = nny + jy(j)
                     ijz(ij) = nnz + jz(j)
 20               continue
 30            continue
               do 40 i = 1 , ij
                  ss(i) = dzero
 40            continue
c
c     ----- i primitive
c
               jgmax = j2
               do 270 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*pi
                  ayi = ai*qi
                  azi = ai*ri
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
                  cgi = cg(ig)
c
c     ----- j primtive
c
                  if (oianj) jgmax = ig
                  do 260 jg = j1 , jgmax
                     aj = ex(jg)
                     aa = ai + aj
                     csj = cs(jg)
                     cpj = cp(jg)
                     cdj = cd(jg)
                     cfj = cf(jg)
                     cgj = cg(jg)
                     ax = (axi+aj*pj)/aa
                     ay = (ayi+aj*qj)/aa
                     az = (azi+aj*rj)/aa
                     dum = aj*arri/aa +
     +                     aa*((p0-ax)**2+(q0-ay)**2+(r0-az)**2)
                     if (dum.le.tol) then
                        fac = dexp(-dum)
c
c     ----- density factor
c
                        odoub = oianj .and. ig.ne.jg
                        max = maxj
                        nn = 0
                        do 220 i = mini , maxi
                         go to (50,60,120,120,
     +                          70,120,120,80,120,120,
     +                          90,120,120,100,120,120,120,120,120,110,
     +                         112,120,120,114,120,120,120,120,120,116,
     +                         120,120,118,120,120), i
c
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
 120                       if (oianj) max = i
                           do 210 j = minj , max
                              go to (130,140,200,200,
     +                         150,200,200,160,200,200,
     +                         170,200,200,180,200,200,200,200,200,190,
     +                         192,200,200,194,200,200,200,200,200,196,
     +                         200,200,198,200,200),j
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
c     ----- density integrals -----
c
                        in = -5
                        do 240 i = 1 , lit
                           in = in + 5
                           ni = i
                           do 230 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call denint
                              pin(jn) = pint
                              qin(jn) = qint
                              rin(jn) = rint
 230                       continue
 240                    continue
                        do 250 i = 1 , ij
                           nnx = ijx(i)
                           nny = ijy(i)
                           nnz = ijz(i)
                           ss(i) = ss(i) + dij(i)*pin(nnx)*qin(nny)
     +                             *rin(nnz)
 250                    continue
                     end if
 260              continue
 270           continue
c
c     ----- set up density integral matrix -----
c
               max = maxj
               nn = 0
               do 290 i = mini , maxi
                  li = loci + i
                  in = (li*(li-1))/2
                  if (oianj) max = i
                  do 280 j = minj , max
                     lj = locj + j
                     jn = lj + in
                     nn = nn + 1
                     q(jn-1+i30) = ss(nn)
 280              continue
 290           continue
 300        continue
 310     continue
         if (out) then
            write (iwr,6010) iat
            call prtril(q(i30),l1)
         end if
c
c      original constants removed (see below)
c
       q(iat-1+i40) = tracep(q(i10),q(i30),l1)-tracep(q(i20),q(i30),l1)
c
 320  continue
c
c - punchfile output
c
      call blkspd(q(i40))
c
c     generate the original multiplicative factor on expectation value
c
      dum = (conpi/dble(na-nb)) / two
      do 330 i = 1 , nat
         if (czan(i).gt.dzero) then
         dum1 = dum * q(i-1+i40)
          write (iwr,6030) i , zaname(i) , czan(i) , q(i-1+i40), dum1
         endif
 330  continue
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,'anala.m','spind','i10')
      return
 6010 format (/10x,22('-')/10x,'spin density integrals',10x,' atom ',
     +        i5/10x,22('-'))
 6020 format (//1x,80('=')//10x,19('-')/10x,'atomic spin density'/10x,
     +        19('-')/
     +        50x,'( * pi/(2*(nalpha-nbeta))'/)
 6030 format (10x,i5,2x,a8,2x,f6.1,f12.5, 5x, '(',f12.5,' )' )
c
      end
**==symtrd.f
      subroutine symtrd(sh,d,shsq,dsq,dt,dod,n)
c
c     ----- transform the density matrix to the symmetrically
c           orthogonal basis set -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension sh(*),d(*),dt(n,*),dod(*)
      dimension shsq(n,*),dsq(n,*)
c
      call square(shsq,sh,n,n)
      call square(dsq,d,n,n)
      call mxma(shsq,1,n,dsq,1,n,dt,1,n,n,n,n)
c
      call mxma(dt,1,n,shsq,1,n,dsq,1,n,n,n,n)
      ij = 0
      do i = 1 , n
      dod(i) = dsq(i,i)
         do j = 1 , i
         ij = ij + 1
         sh(ij) = dsq(j,i)
         enddo
      enddo
c
      return
      end
**==symtrn.f
      subroutine symtrn(s,sv,se,ia,b,n,n2,ndim,string)
c
c     ----- calculate the symmetric s **(+1/2) matrix
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension s(*),sv(ndim,*),se(*),ia(*),b(*)
      character *(*) string
c
c     ----- diagonalize the overlap matrix -----
c
      call gldiag(n,n,ndim,s,b,se,sv,ia,2)
c
c     ----- take the square root of the eigenvalues
c     ----- check for too small or < 0
c
      if (se(1).lt.1.0d-14) then
         write(iwr,*) (se(i),i=1,n)
         write(iwr,*) ' symtrn called from ',string
         call caserr('basis dependent (symtrn)') 
      endif
c
      call vsqrt(se,1,se,1,n)
      call vclr(s,1,n2)
c
c     ----- form s **(+1/2) -----
c
      ij = 1
      do 30 i = 1 , n
         do 20 k = 1 , n
            fact = sv(i,k)*se(k)
            if (fact.ne.0.0d0) then
               call daxpy(i,fact,sv(1,k),1,s(ij),1)
            end if
 20      continue
         ij = ij + i
 30   continue
c
      return
      end
**==symtrv.f
      subroutine symtrv(sh,v,vo,vo2,m,n,ndim)
c
c     ----- transform the orbitals to the symmetrically
c           orthogonalized basis set -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension sh(*),v(ndim,*),vo(ndim,*),vo2(ndim,*)
c
      call square(vo2,sh,m,n)
      call mxma(vo2,1,ndim,v,1,ndim,vo,1,ndim,m,n,n)
      do 50 j = 1 , m
         do 40 i = 1 , n
            vt = vo(i,j)
            vo2(i,j) = vt*vt
            if (vt.lt.0.0d0) vo2(i,j) = -vo2(i,j)
 40      continue
 50   continue
c
      return
      end
**==transd.f
      subroutine transd(zzzzzz,q,nprint)
c
c     ----- transform the d basis functions from xx,yy, and zz to
c           rr, xx-yy, 3zz-rr -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical otri
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
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      dimension dtr(3,3),dtrinv(3,3)
      dimension q(*)
c
      data zrhf,zcas,zmcscf/'rhf','casscf','mcscf'/
      data dtr/
     + 0.4472135954999579d0, 0.4472135954999579d0, 0.4472135954999579d0,
     + 0.8660254037844386d0,-0.8660254037844386d0, 0.0000000000000000d0,
     +-0.5000000000000000d0,-0.5000000000000000d0, 1.0000000000000000d0/
      data dtrinv/
     + 0.7453559924999299d0, 0.7453559924999299d0, 0.7453559924999299d0,
     + 0.5773502691896258d0,-0.5773502691896258d0, 0.0000000000000000d0,
     +-0.3333333333333333d0,-0.3333333333333333d0, 0.6666666666666667d0/
c
c
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c     ----- set pointers for partitioning of memory -----
c
c     -------- -------- -------- -------- ------- -------
c     i10      i20      i30      i40      i50     i60
c     -------- -------- -------- -------- ------- -------
c     l2       l2       l2       l2       l2      l2
c     -------------- -- -------------- -- ------------ --
c                    i21
c     -------------- -- 
c     l3             l1 
c
      otri = .false.
      ntri = 6
      i10 = igmem_alloc_inf(ntri*l2,'anala.m','transd',
     1                      'i10',IGMEM_DEBUG )
      i20 = i10 + l2
      i21 = i10 + l3
      i30 = i20 + l2
      i40 = i30 + l2
      if (otri) then
       i50 = i10
      else
       i50 = i40 + l2
      endif
      i60 = i50 + l2
      last = i60 + l2

      if (nprint.eq.5) write (iwr,6030) i10 , i20 , i21 , i30 , 
     +                                  i40 , i50,  i60, last
c
c     ----- transformation matrix at q(i10) -----
c
      call trinit(q(i10),dtr,l1,odbas,ifbf)
c
      if (odbas) then
c
c     ----- transform the overlap matrix -----
c
c     -s - at q(i30)
c     -st- at q(i50)
c
         call rdedx(q(i30),l2,ibl7s,num8)
c
c     ----- check if the overlap is already orthogonal -----
c
         locodi = iky(ifbf+1) + ifbf
         if (q(i30-1+locodi).gt.0.001d0) then
            call mult2(q(i10),q(i50),q(i30),l1,l1,l1)
            call wrt3(q(i50),l2,ibl7s,num8)
            if (nprint.eq.5) call prtril(q(i50),l1)
c
c     ----- transform the dipole matrix elements -----
c
c     xyz  at q(i30)
c     xyzt at q(i50)
c
            do i = 1 , 3
               if(i.eq.1) then
                mbld = ibl7x
               else if(i.eq.2) then
                mbld = ibl7y
               else
                mbld = ibl7z
               endif
               call rdedx(q(i30),l2,mbld,num8)
               call trinit(q(i10),dtr,l1,odbas,ifbf)
               call mult2(q(i10),q(i50),q(i30),l1,l1,l1)
               call wrt3(q(i50),l2,mbld,num8)
            enddo
c
c     ----- transform the density matrix -----
c
c     -d - at x(i30)
c     -dt- at x(i40)
c
            ipass = 1
            mblq = ibl7qa
            mblp = ibl7pa
            mble = ibl7ea
c
 110        continue
            call rdedx(q(i30),l2,mblp,num8)
c
c      ----- calculate the inverse of the transformation matrix -----
c
            call trinit(q(i10),dtrinv,l1,odbas,ifbf)
            call mult2(q(i10),q(i50),q(i30),l1,l1,l1)
            call wrt3(q(i50),l2,mblp,num8)
c
c     ----- transform the eigenvectors -----
c
c     -v - at x(i30)
c     -vt- at x(i30)
c
            if (ipass.eq.1) then
               call rdedx(q(i30),l3,mblq,num8)
               call trinit(q(i10),dtrinv,l1,odbas,ifbf)
               call vclr(q(i50),1,l3)
               call mxmb(q(i10),l1,1,q(i30),1,l1,q(i50),1,l1,l1,l1,l1)
               call wrt3(q(i50),l3,mblq,num8)
               if (nprint.eq.5) call prev(q(i50),q(i10),l1,l1,l1)
            end if
c
c     ----- print actual value of energy -----
c
c     -v- at x(i10)
c     -d- at x(i30)
c     -e- at x(i21)
c
            if (nprint.ne.5 .and. nprint.ne.-5) then
               if (oprint(43)) then
                  if (ipass.eq.1) call rdedx(q(i10),l3,mblq,num8)
                  call rdedx(q(i30),l2,mblp,num8)
                  if (ipass.eq.1) call rdedx(q(i21),l1,mble,num8)
                  if (ipass.eq.1) write (iwr,6010)
                  if (ipass.eq.1) call prev(q(i10),q(i21),l1,l1,l1)
                  write (iwr,6020)
                  call prtril(q(i30),l1)
               end if
            end if
c
c     ----- reset core memory -----
c
            if (ipass.ne.2) then
               if (zzzzzz.ne.zrhf.and.zzzzzz.ne.zcas.and.
     +             zzzzzz.ne.zmcscf) then
                  ipass = 2
                    mblq = ibl7qb
                    mblp = ibl7pb
                    mble = ibl7eb
                  go to 110
               end if
            end if
         end if
      end if
      call gmem_free_inf(i10,'anala.m','transd','i10')
      return
 6010 format (/1x,80('=')//10x,38('-')/10x,
     +        'eigenvectors (d functions transformed)'/10x,38('-'))
 6020 format (//1x,80('=')//10x,40('-')/10x,
     +        'density matrix (d functions transformed)'/10x,40('-'))
 6030 format (' core assignement',/' i10, i20, i21,',
     +        '  i30, i40 ,i50, i60  = ',/,7i8,/,' last = ',i8)
 6040 format(/1x,' **** 4-triangle GA-based code invoked *****'/)
      end
**==trinit.f
      subroutine trinit(trns,dtr,l1,odbas,ifbf)
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
      dimension dtr(3,3),zdlab(3)
      dimension trns(l1,*)
c
      data zdlab/'     rr ','   xx-yy','   zz-rr'/
c
c
c     ----- initialize the transformation matrix to the unit matrix ----
c
      odbas = .false.
      call vclr(trns,1,num*num)
      do 20 i = 1 , num
         trns(i,i) = 1.0d0
 20   continue
      ifbf = 0
c
c     ----- check each shell for d basis functions -----
c
      do 60 i = 1 , nshell
         if (ktype(i).eq.3) then
            odbas = .true.
            ifbf = kloc(i)
            icol = ifbf
            do 40 j = 1 , 3
               irow = ifbf
               do 30 k = 1 , 3
                  trns(irow,icol) = dtr(k,j)
                  irow = irow + 1
 30            continue
               icol = icol + 1
 40         continue
            do 50 j = 1 , 3
               zbflab(ifbf+j-1) = zdlab(j)
 50         continue
         end if
 60   continue
c
      return
      end
**==uhfnat.f
      subroutine uhfnat(q,nsav,iprinv,ispin)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c..   get spinfree natural orbitals from uhf function
c..   may also be useful in starting up casscf (cf. Pulay - UNOCAS)
c..   ispin = 0 : spinfree natural orbitals
c..   ispin = 1 : spin natural orbitals
c     oanil  = .true. : generate annihilated density matrices
c
c      density matrices need to be transformed inverse
c..
c..     ** called from hfprop **
c..      jvl Utrecht 1991
c
      dimension q(*)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      dimension ztype(2)
c
      data ztype/'spinfree','  spin  '/
c
      if (ispin.eq.0) zcom(5) = 'uhfno'
      if (ispin.eq.1) zcom(5) = 'uhfsno'
c..
c..    core partitioning
c..
       l1 = num
       l2 = num*(num+1)/2
       l3 = num*num
c
c      first compute core requirements
c
       i10 = 0
       i20 = i10 + l2 + l2
       i30 = i20 + l2 + l2
       i50 = i30 + l2 + l2
       i60 = i50 + l1
c       i50 is scratch for minv (crayxmp) : 2*l1
c       i50,i60 scratch for minvrt (others) l1 each
       last = i60 + l1
c
c      for annihilation
c
       if(oanil) then
        i21 = i10 + l2
        i31 = i21 + l2
        i41 = i31 + l2
        i51 = i41 + l2
        i61 = i51 + l2
        lasta = i61 + l2
        last = max (last,lasta)
       endif
c
c..    now get core and assign as above
c
       i10 = igmem_alloc_inf(last,'anala.m','uhfnat',
     1                       'i10',IGMEM_DEBUG )
c
       i20 = i10 + l2 + l2
       i30 = i20 + l2 + l2
       i50 = i30 + l2 + l2
       i60 = i50 + l1
c
c      for annihilation
c
       if(oanil) then
        i21 = i10 + l2
        i31 = i21 + l2
        i41 = i31 + l2
        i51 = i41 + l2
        i61 = i51 + l2
       endif
c
c     ----- get uhf density matrices -----
c
c     -da- at q(i10)
c     -db- at q(i20)
c     then total at q(i10)
c    ispin 0 : +   // ispin 1 : -
c
      if(oanil) then
        call anil(q(i10),q(i21),q(i31),q(i51),q(i61))
        call dcopy(l2,q(i21),1,q(i20),1)
      else
        call rdedx(q(i10),l2,ibl3pa,idaf)
        call rdedx(q(i20),l2,ibl3pb,idaf)
      endif
      if (ispin.eq.0) then
         call vadd(q(i10),1,q(i20),1,q(i10),1,l2)
      else
         call vsub(q(i10),1,q(i20),1,q(i10),1,l2)
      end if
c
c..   transform density matrix to orthonormal basis (UHF A-vectors)
c
      call rdedx(q(i20),l3,ibl3qa,idaf)
      call tdown(q(i20),ilifq,q(i20),ilifq,l1)
      call vclr(q(i50),1,2*l1)
c
c      invert transformation matrix t -> t(inv)
c
      d = 0.0d0
      call minvrt(q(i20),l1,d,q(i50),q(i60))
c
c      transpose it  t(inv) -> t(inv(dagger))
c
      call dagger(l1,l1,q(i20),l1,q(i30),l1)
c
c      then make t(inv).d.t(inv(dagger))
c      i30: t / i20: output d / i10: input d
c
      call mult2(q(i30),q(i20),q(i10),l1,l1,l1)
c
c..    diagonalise orthonormalised d-matrix (sort in decreasing order)
c..    result : occupations in i50 ; vectors in i30
c
      call gldiag(l1,l1,l1,q(i20),q(i60),q(i50),q(i30),iky,3)
c
c...  back-transform the eigenvectors (with symmetry adapted vectors)
c
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call tfsqc(q(i30),q(i10),q(i20),l1,l1,l1)
c
c...   save ** mo's and eigenvalues = occupations ** on section nsav
c
      if (nsav.ne.0) then
         if (iprinv.gt.0) then
           if (oanil) then
            write(iwr,9027) ztype(ispin+1),nsav
           else
            write(iwr,9028) ztype(ispin+1),nsav
           endif
         endif
         call putq(zcom,ztitle,q(i50),q(i50),l1,l1,l1,1,1,
     *             q(i30),nsav,ibl3qa)
      end if
c
c..     transform back to original basis and print if requested
c
      if (iprinv.gt.1) then
         call tdown(q(i30),ilifq,q(i30),ilifq,l1)
         if (oanil) then
          write(iwr,9029) ztype(ispin+1)
         else
          write(iwr,9030) ztype(ispin+1)
         endif
         call prev(q(i30),q(i50),l1,l1,l1)
      endif
      if (iprinv.eq.1) then
         if (oanil) then
          write(iwr,9031) ztype(ispin+1)
         else
          write(iwr,9032) ztype(ispin+1)
         endif
         write(iwr,9040) (q(i50+i-1),i=1,l1)
         write(iwr,9040)
      end if
c
      call gmem_free_inf(i10,'anala.m','uhfnat','i10')
c
      return
9027  format(/' ** annihilated ',a8,
     +        ' UHF natural orbitals saved in section',i4,' **')
9028  format(/' ** ',a8,' UHF natural orbitals saved in section'
     +       ,i4,' **')
9029  format(//30x,40('-')/30x,'annihilated ',a8,
     +      ' uhf natural orbitals'/30x,40('-')/)
9030  format(//30x,29('-')/30x,a8,' uhf natural orbitals'/
     +       30x,29('-')/)
9031  format(//10x,'----- annihilated ',a8,
     +        ' UHF natural orbital occupations -----')
9032  format(//10x,'----- ',a8,' UHF natural orbital occupations -----')
9040  format(/10x,7f14.7)
      end
      subroutine ver_anala(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/anala.m,v $
     +     "/
      data revision /
     +     "$Revision: 6321 $"
     +      /
      data date /
     +     "$Date: 2015-03-15 01:39:24 +0100 (Sun, 15 Mar 2015) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
      subroutine anil(q1,q2,q3,q5,q6)
      implicit real*8 (a-h,o-z)
      logical o1e
      character *4 annil
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
      dimension q1(*),q2(*),q3(*),q5(*),q6(*)
      dimension coa(8),cob(8),mmi(8),o1e(6)
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
      common/blkcore/corev(512),charge(4),cpad(6)
      common/junk/coa1,coa2,coa3,coa4,coa5,coa6,coa7,coa8,
     +            cob1,cob2,cob3,cob4,cob5,cob6,cob7,cob8,
     +            i1,i2,i3,i4,i5,i6,i7,i8
      equivalence (coa(1),coa1),(cob(1),cob1),(mmi(1),i1)
      data stwo,quart,half,one,two,thr,four,six,ten/
     +     -2.0d0,0.25d0, 0.5d0 ,1.0d0,2.0d0,3.0d0,
     +      4.0d0, 6.0d0,10.0d0/
      data annil/'anil'/
      zcom(4)=annil
      l3 = num*num
c
c     restore s matrix
c
      do loop =1,6
       o1e(loop) = .false.
      enddo
      o1e(1) = .true.
      call getmat(q6,q6,q6,q6,q6,q6,charge,num,o1e,ionsec)
c     restore alpha vectors and transform to AO basis
      call rdedx(q3,l3,ibl3qa,idaf)
      call tdown(q3,ilifq,q3,ilifq,num)
      call densum(q1,q3,ilifq,na,num)
      call wrt3(q1,nx,ibl7la,num8)
c
      i1=ibl7la
      i2=iposun(num8)
c
c     restore beta vectors and tdown
      call rdedx(q3,l3,ibl3qb,idaf)
      call tdown(q3,ilifq,q3,ilifq,num)
c
      call densum(q5,q3,ilifq,nb,num)
      call wrt3(q5,nx,i2,num8)
      ja=iposun(num8)
c ... ps....ja
      call anilm1(q1,q6,q3,iky,ilifq,num)
      call wrt3(q3,l3,ja,num8)
c ... qs....jb
      jb=iposun(num8)
      call anilm1(q5,q6,q1,iky,ilifq,num)
      call wrt3(q1,l3,jb,num8)
c ... tr(psqs)
      tab=0.0d0
      do 1 i=1,num
      ii=ilifq(i)
      do 1 j=1,num
 1    tab=tab+q3(ii+j)*q1(ilifq(j)+i)
c ... psq
      call anilm2(q3,q5,q1,iky,ilifq,num)
      ij=0
      do 2 i=1,num
      ii=ilifq(i)
      do 2 j=1,i
      ij=ij+1
 2    q3(ij)=q1(ii+j)+q1(ilifq(j)+i)
c ... (psq+qsp)....i3
      i3=iposun(num8)
      call wrt3(q3,nx,i3,num8)
c ... qs*psq = qspsq .... i7
      i7=iposun(num8)
      call rdedx(q3,l3,jb,num8)
      call anilm3(q3,q1,q5,ilifq,num)
      call wrt3(q5,nx,i7,num8)
c ... ps*qspsq = psqspsq
c ... (psqspsq+qspsqsp) ... i5
      i5=iposun(num8)
      call rdedx(q3,l3,ja,num8)
      call anilm2(q3,q5,q1,iky,ilifq,num)
      ij=0
      do 3 i=1,num
      ii=ilifq(i)
      do 3 j=1,i
      ij=ij+1
 3    q3(ij)=q1(ii+j)+q1(i+ilifq(j))
      call wrt3(q3,nx,i5,num8)
c ... qs * psqspsq = qspsqspsq ... i6
      i6=iposun(num8)
      call rdedx(q3,l3,jb,num8)
      call anilm3(q3,q1,q5,ilifq,num)
      call wrt3(q5,nx,i6,num8)
      i4=iposun(num8)
c ... tr(ps,qspsqspsq*s)
      call rdedx(q3,l3,ja,num8)
      call anilm1(q5,q6,q1,iky,ilifq,num)
      tab3=0.0d0
      do 4 i=1,num
      ii=ilifq(i)
      do 4 j=1,num
 4    tab3=tab3+q3(i+ilifq(j))*q1(ii+j)
c ... qs *p = qsp
      call rdedx(q3,l3,jb,num8)
      call rdedx(q5,nx,i1,num8)
       call anilm2(q3,q5,q1,iky,ilifq,num)
c ... ps * qsp = psqsp .... i4
      call rdedx(q3,l3,ja,num8)
      call anilm3(q3,q1,q5,ilifq,num)
      call wrt3(q5,nx,i4,num8)
      i8=iposun(num8)
c ... tr(psqsp.s.qs)
      call anilm1(q5,q6,q1,iky,ilifq,num)
      call rdedx(q3,l3,jb,num8)
      tabab=0.0d0
      do 6 i=1,num
      ii=ilifq(i)
      do 6 j=1,num
 6    tabab=tabab+q1(i+ilifq(j))*q3(j+ii)
c ... qs * psqsp = qspsqsp
      call anilm2(q3,q5,q1,iky,ilifq,num)
c ... ps * qspsqsp = psqspsqsp ... i8
      call rdedx(q3,l3,ja,num8)
      call anilm3(q3,q1,q5,ilifq,num)
      call wrt3(q5,nx,i8,num8)
      p=na
      q=nb
      fn=p+q
      const=q+q-p-two
      tab2=tab+tab
      coa1=const*(const-tab2)+q*(p-one)+(thr-fn+tab2)*tab
     +    -tabab-tabab
      cob1=q-tab
      coa2=p-tab
      cob2=coa1+q-p
      tab4=tab2+tab2
      coa3=tab2+one-p-const
      cob3=tab2+one-q-const
      coa4=fn-tab4-thr+const+const
      cob7=coa4
      coa5=stwo
      cob5=stwo
      coa6=0.0d0
      cob8=0.0d0
      coa7=one
      cob4=one
      coa8=four
      cob6=coa8
      cc=fn*(fn+two)*quart-p*q
      s2sd=cc-tab
      s4sd=cc*cc+p*q+two*(tab*tab-tabab)
     +       -(two*cc+fn-two)*tab
      s6sd=cc*cc*cc + cc*p*q + p*q*(two*cc+fn-two)
     +    -tab*(thr*cc*cc+thr*cc*(fn-two)+(fn-two)**2
     +    +p*q+four*(p-one)*(q-one))+two*(thr*cc+thr*fn-ten)
     +    *(tab**2-tabab)-six*(tab**3-thr*tab*tabab+two*tab3)
      sdash=half*(p-q)
      sda1=sdash+one
      sda2=sdash+two
      s2aa=(s6sd-two*sda1*sda2*s4sd+
     +    (sda1*sda2)**2*s2sd)/
     +    (s4sd-two*sda1*sda2*s2sd+
     +    (sda1*sda2)**2)
      write(iwr,7)s2sd,s2aa
7     format(/' uhf (after annihilation) analysis'/1x,32('-')//
     + 20x,'<s**2> before annihilation   =',f16.7/
     + 20x,'<s**2> after  annihilation   =',f16.7)
      anorm=const*const-two*const*tab+p*q-fn*tab
     +    +tab2+two*tab*tab-two*tabab
      anorm=one/anorm
      call vclr(q1,1,nx)
      call vclr(q2,1,nx)
      do 8 i=1,8
      coal=coa(i)
      cobl=cob(i)
      call rdedx(q3,nx,mmi(i),num8)
      do 8 j=1,nx
      cc=q3(j)
      q1(j)=q1(j)+coal*cc
 8    q2(j)=q2(j)+cobl*cc
      do 9 i=1,nx
      q1(i)=q1(i)*anorm
 9    q2(i)=q2(i)*anorm
      return
      end
      subroutine anilm1(a,b,c,map,ilif,nbasis)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*),map(*),ilif(*)
c ... c=a*b  (a and b are triangles)-no overwrite
      do 1 j=1,nbasis
      jj=map(j)
      iii=ilif(j)
      do 1 i=1,nbasis
      ii=map(i)
      ij=iii+i
      x=0.0d0
      do 7 k=1,nbasis
      kk=map(k)
      if(i.ge.k)go to 3
      ik=i+kk
      go to 4
 3    ik=k+ii
 4    if(j.ge.k)go to 6
      jk=j+kk
      go to 7
  6   jk=k+jj
 7    x=x+a(ik)*b(jk)
 1    c(ij)=x
      return
      end
      subroutine anilm2(a,b,c,map,ilif,nbasis)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*),map(*),ilif(*)
c ... c=a*b (b is a triangle)   - no overwrite
      do 8 j=1,nbasis
      jj=map(j)
      ii=ilif(j)
      do 8 i=1,nbasis
      x=0.0d0
      ij=ii+i
      do 9 k=1,nbasis
      if(j.ge.k)go to 11
      jk=j+map(k)
      go to 9
 11   jk=k+jj
 9    x=x+a(i+ilif(k))*b(jk)
 8    c(ij)=x
      return
      end
      subroutine anilm3(a,b,c,ilif,nbasis)
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*),ilif(*)
c ... c=a*b  c is a triangle -- no overwrite
      ij=0
      do 13 i=1,nbasis
      do 13 j=1,i
      jj=ilif(j)
      ij=ij+1
      x=0.0d0
      do 14 k=1,nbasis
 14   x=x+a(ilif(k)+i)*b(jj+k)
 13   c(ij)=x
      return
      end
      subroutine densum(p,a,ilif,nocc,nbasis)
      implicit real*8 (a-h,o-z)
      dimension a(*),p(*),ilif(*)
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
      call vclr(p,1,iky(nbasis+1))
      do 32 k=1,nocc
      m=ilif(k)
      ij=1
      do 32 i=1,nbasis
      call daxpy(i,a(m+i),a(m+1),1,p(ij),1)
32    ij=ij+i
      return
      end
c****  excanon
      subroutine excanon(v,l1,iuno,q)
c
      implicit real*8 (a-h,p-w),integer   (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...  subroutine to perform a canonicalisation of natural orbitals
c...  by subset as specified in the canon directive
c...  extra routines/commons are :
c...  incanon
c...  common/canon_nat/ in incanon, excanon mains
c
c...   reads vectors from iuno (specified on natorb) 
c...   write result to isec_canon (on canon subdirective)
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
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      real*8 aoccc
      integer isetc, iync, ns_canon, isec_canon
      logical opr_canon,o_canon,oset_canon
      common/canon_nat/aoccc(10),isetc(10),iync(10),
     1                 ns_canon,isec_canon,
     2                 opr_canon,o_canon,oset_canon
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/junk/eig(maxorb),frocc(maxorb)
      common/junkc/zcomm(19),ztita(10)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)       
      common/scra/iso(mxshel,48)
c
      dimension v(l1,l1)
      dimension q(*)
c
c
      l2 = l1*(l1+1)/2
      l3 = l1*l1
c
      idens = igmem_alloc_inf(4*l2+2*l2,'anala.m','excanon',
     1                        'idens',IGMEM_DEBUG)
      idmat = idens
      ifock = idens + l2
      itri = ifock + l2
      iexch = itri
      itri2 = itri + l2
      isq = itri2 + l2
      irdmat = igmem_null()
      iprefa = igmem_null()
      if (odscf) then
         irdmat = iexch + l2
         iprefa = irdmat + ikyp(nshell)
      endif
      write(iwr,606) iuno,isec_canon
606   format(/' C A N O N I C A L I S A T I O N of NATURAL ORBITALS',
     1       /' vectors read from section',i4,' and written to',i4)
c
c.... get vectors
c
      iposv = iuno
      if (iposv.le.0 .or. iposv.gt.350)
     +    call caserr2('invalid section specified for eigenvectors')
c
c...  read natural orbitals // suppress print
c
      ndum = nprint
      nprint = -5
      call getq(v,eig,frocc,nbas,newbas,1,ieig,ipop,iposv,'nat-')
      nprint = ndum
c
      if (ipop.ne.1)call caserr2('vectors without occup. in excanon')
      zcomm(5) = 'canonic'
      if (nbas.ne.l1) call caserr2('vector error in canon')
c...  unadapted in isq
      call tdown(q(isq),ilifq,v,ilifq,newbas)
c
c...  build density matrix
c...  print the occupations
c
      n = 0
      nn = 0
      do i=1,l1
         if (dabs(frocc(i)-2.0d0).gt.1.0d-13) go to 1
         n = n + 1
      end do
1     do i=l1,1,-1
         if (dabs(frocc(i)).gt.1.0d-13) go to 2
         nn = nn + 1
      end do
2     write(iwr,6) n,l1-n-nn,nn
      if (l1-n-nn.gt.0) write(iwr,7) (frocc(i),i=n+1,l1-nn)
6     format(' # Doubly occupied',i5,' # Variably occupied',i5,
     1       ' # Empty orbitals',i5)
7     format(' Occupations ',5f12.8,/,(13x,5f12.8))
c
      if (.not.oset_canon) then
c
c...  translate occ to set
c
       do i=2,l1
        if (frocc(i).gt.frocc(i-1)) call caserr2(
     +                          'unordered nos - canon')
       end do
c
10     if (ns_canon.eq.0) go to 30
       if (frocc(1).lt.aoccc(1)) then
          ns_canon = ns_canon - 1
          do i=1,ns_canon
             aoccc(i) = aoccc(i+1)
             iync(i) = iync(i+1)
          end do
          go to 10
       end if
       do is=1,ns_canon
         do i=1,l1
           if (frocc(i).lt.aoccc(is)) go to 20
         end do
20       isetc(is) = i-1
       end do
      end if  
c
30    continue
c
c...  ** empty sets are just skipped
c
c...  build density matrix
c
      call vclr(q(idens),1,l2)
      do i=1,l1
         pop = frocc(i)
         if (pop.gt.1.0d-13) then
            kl = 0
            ii = (i-1)*l1+isq-1
            do k=1,l1
               do l=1,k
                  q(idens+kl) = q(idens+kl) + pop*q(k+ii)*q(l+ii)
                  kl = kl + 1
               end do
            end do
         end if
      end do

c
c     ----- read in transformation matrices for s,p,d,f and g basis
c           functions.
c
      nav = lenwrd()
      if (nt.gt.1) then
       call rdedx(ptr,nw196(1),ibl196(1),idaf)
       if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
       if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
       if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
       call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      endif

c
c.... build either fock-matrix or coulomb matrix
c
      zsave = zscftp
      zscftp = 'rhf'
c
c
c Modify fock builder options
c
      idum = CD_set_2e()

c
      if (odscf) then
         call rdmake(q(iprefa))
         call mkrdmt('canon',q(irdmat),q(idmat),l2,nprint)
         dlntol = tolitr(1)
         call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),0)
      else
         call hstar(q(idens),q(ifock),q(iexch),nopk)
      end if

c
c restore fock builder options
c
      idum = CD_reset_2e()
c
c...  add density functional contribution
c
      idum = CD_rks()
      if (CD_active()) then
         idum= CD_energy_ao(flop,q(ifock),flop,q(idens),flop,
     1                   energy,q,q,.false.,0.0d0,iwr
     3                   )
      endif
c
c.... symmetrize skeleton fock matrix
c
      call symh(q(ifock),q(itri),iky,0,0)
c
c     ---- read in core hamiltonian matrix, one electron ints in q(i20)
c     - h - at q(i10)
c
      call rdedx(q(itri),l2,ibl7f,num8)
c
      if (ozora) then
            call zora(q,q(idmat),q(itri),'read')
      end if
c
      call vadd(q(ifock),1,q(itri),1,q(itri2),1,l2)
c
c...  adapt
c
      call tranp(q(itri2),q(idens))
c
      zscftp = zsave
c
c...  we have a Fock matrix in adapted basis and vectors in the basis
c...  (fock is kept safe in q(idens)
c
c...  canonicalise in sets
c
      isetc(ns_canon+1) = l1
      ib = 1
      call vclr(eig,1,l1)
      call vclr(frocc,1,l1)
c
      do is=0,ns_canon
       call dcopy(l2,q(idens),1,q(ifock),1)
       ie = isetc(is+1)
       ns = ie-ib+1
       if (ns.eq.0.or.iync(is+1).eq.0) go to 40
c...    ** itri2/ifock effectively square (in mult2)
       call mult2(v(1,ib),q(itri2),q(ifock),ns,ns,l1)
       call jacobi(q(itri2),iky,ns,q(isq),ilifq,l1,eig(ib),2,2,1.0d-10)
c...    ** v,isq,itri square
       call tfsqc(q(isq),v(1,ib),q(itri),ns,l1,l1)
       call dcopy(ns*l1,q(isq),1,v(1,ib),1)
40     ib = ie+1
      end do
c
c...  print and save to disk
c
      write(iwr,603) 
603   format(/' vectors are canonicalised in sets',/,
     1       ' set from to    orbital energies')
      ib = 1
      do is=0,ns_canon
       ie = isetc(is+1)
       if (ib.gt.ie.or.iync(is+1).eq.0) then
          if (ib.le.ie) write(iwr,602) is+1,ib,ie
       else
          write(iwr,601) is+1,ib,ie,
     1                   (eig(i),i=ib,ie)
       end if
601    format(i3,1x,2i4,(t13,8f12.5))
602    format(i3,1x,2i4,(t16,'*not* canonicalised'))
       ib = ie + 1
      end do
c
      call putq(zcomm,ztita,eig,frocc,l1,l1,l1,1,0,v,isec_canon,iblk)
      write(iwr,604) isec_canon
604   format(/' canonicalised vectors saved in section',i5)
c
      if (opr_canon) then
         write(iwr,605)
605      format(/' canonicalised vectors  ..... ')
         call tdown(v,ilifq,v,ilifq,l1)
         call prev(v,eig,l1,l1,l1)
      end if
c
      call gmem_free_inf(idens,'anala.m','excanon','idens')
c
      return 
      end
