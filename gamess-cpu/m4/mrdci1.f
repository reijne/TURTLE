c 
c  $AuTor: wab $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mrdci1.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   Table-ci (adapt module) =
c ******************************************************
c ******************************************************
      subroutine amrdci(x,lword)
      implicit real*8  (a-h,o-z), integer (i-n)
      logical orout,route
      character *4 fdump
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
      common /blkin/ pp(4),y(507)
      common/craypk/igam,orout,isec33,numsp
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
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common /linkmr/ imap(510),ipr,lbuff,igame
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
      common/lsort /atomo(4,maxat),poto,nnsho(maxat)
      common/b/sum(2),itape(2),ntape,mtape,ltape,mscr
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      common/restri/nfils(63),lds(508),isect(508),ldsect(508)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      dimension iq(6),x(*)
c
      data fdump/'dump'/
      data m1,m2,m511/1,2,511/
c
c
      if(.not.oprint(29)) write(iwr,35)
      igame=igam
      route=orout
      if(route) then
      m1file=n1file
      do 1020 loop=1,m1file
      m1tape(loop)=n1tape(loop)
      m1blk (loop)=n1blk(loop)
 1020 m1last(loop)=n1last(loop)
      endif
      lfile=m2file
      do 1000 loop=1,lfile
      lotape(loop)=m2tape(loop)
      liblk (loop)=m2blk(loop)
 1000 llblk (loop)=m2blk(loop)-m2last(loop)
      isec3=isec33
      if(isec3.gt.350.and.isec3.ne.482)call caserr(
     *'invalid section specified for 1-electron integrals')
      if(isec3.ne.0)ions2=isec3
      if(.not.oprint(29)) then
        write(iwr,3)fdump,yed(idaf),ibl3d
        write(iwr,32)
        call filprn(n2file,n2blk,n2last,n2tape)
        if(route)go to 104
        write(iwr,105)
        go to 106
 104    write(iwr,33)ions2
        write(iwr,36)
        call filprn(m1file,m1blk,m1last,m1tape)
        if(ions2.eq.ionsec)write(iwr,34)
 106    write(iwr,18)lword
      endif
c        1st tackle the dfile
      call secget(isect(491),m1,iblkv)
c
c ----- decide if this section has originated from an
c       atmol or gamess run - modified for arbitrary mxprim
c
      call search(iblkv,idaf)
      nblk1 = lensec(mxprim) - 1
      nw1 = mxprim - nblk1 * 511
      if (nblk1. gt. 0) then
       do 3000 loop = 1, nblk1
       call find(idaf)
       call get(pp(1),nw)
 3000  continue
      endif
      call find(idaf)
      call get(pp(1),nw)
      if(nw.eq.nw1) then
          call gamrd(iblkv)
      else
          call atmrd(iblkv)
      endif
c
c          restore potnuc
c
      call secget(ionsec,m2,iblk33)
      call rdedx(pp,m511,iblk33,idaf)
      poto=pp(1)
c
      call rewftn(nf2)
      lenbas=num*(num+1)/2
      iq(1)=1
      do 26 i=2,4
 26   iq(i)=iq(i-1)+lenbas
      iq(5) = iq(4) + lenbas
      iq(6) = iq(5) + 255 * mxprms
      last =  iq(6) + 255 * mxprms
      if(last.gt.lword)call caserr(
     * 'insufficient memory available')
c
      call face1(x(1),lword,iwr)
      call amrdm(x(1),x(iq(1)),x(iq(2)),x(iq(3)),x(iq(4)),
     +           x(iq(5)),x(iq(6)),lword)
      call face2(x(iq(1)),x(iq(2)),x(iq(3)),x(iq(4)),route)
c     now rewind data sets (Pentium problem)
      call rewftn(ntape)
      call rewftn(mtape)
      if (ltape.ne.ntape) then
       call rewftn(ltape)
      endif
      call rewftn(mscr)
      return
 18   format(' main core available = ',i8,' words'/)
 34   format(/
     *' **** note : original 1-electron integrals to be overwritten'
     */)
 36   format(/1x,'symmetry adapted integral files'/
     *1x,32('*'))
 33   format(/' *** sabf routing specified'//
     * ' one-electron integrals to section ',i3,' of dumpfile')
 105  format(/' *** no sabf routing specified'/)
 32   format(/1x,'2-electron integral files'/1x,25('*'))
 3    format(/1x,a4,'file on ',a4,' at block',i6 )
 35   format(//1x,104('=')//40x,35('*')/
     * 40x,'Table-Ci  --  symmetry adapt module'/40x,35('*')/)
      end
      subroutine atmrd(iblkv)
      implicit real*8  (a-h,o-z), integer (i-n)
      character *8 ztitle,tagg
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
      common/junk/ishp(mxgaus),ityp(mxgaus),igp(mxgaus),exxp(2,mxgaus)
      common /lsort / atomo(4,maxat),poto,nsho(maxat)
      common /scra  / zet(10,160),anorm(10,160),
     + cx(100),cy(100),cz(100),z(100),itype(160),
     + itypee(160),
     + icentr(160),nbase(160),ntrm(160),title(2),ititl(8),
     + norb,non,ibl(3),ngrps,ncentg(2),cgx(36),idddd,norbd,
     + mmbase(160)
      common/junkc/ztitle(10),tagg(100)
      common/linkmr/new(255),newei(255),ipr,lbuff
      dimension ngrpp(4)
c
      data ngrpp/1,3,6,10/
c
      nav = lenwrd()
c     nump=mxgrps*(20+6/nav)+maxat*4+38+18/nav
c     m110=10+maxat
      nump=160*(20+6/nav)+100*4+38+18/nav
      m110=110
      call rdedx(zet,nump,iblkv,idaf)
      call rdchrs(ztitle,m110,idaf)
      if(non.le.0.or.non.gt.maxat) go to 6000
      if(ngrps.le.0.or.ngrps.gt.mxgrps) go to 6000
      if(norb.gt.0.and.norb.le.maxorb) go to 6010
6000  call caserr(
     *'parameter error detected on atmol dumpfile')
6010  nat=non
      num=norb
      if(idddd.ne.0) call caserr(
     *'spherical harmonic option not permitted')
c
      if(oprint(31))write(iwr,1)
 1    format(
     *' section 191 created by atmol integrals program *****'
     * )
c
c       primitive normalisation (xx)
c
      pi=3.14159265359d0
      pip=1.0d0/(2.0d0**0.25d0*pi**0.625d0)
c
      do 20 i=1,non
      atomo(1,i)=z(i)
      atomo(2,i)=cx(i)
      atomo(3,i)=cy(i)
      atomo(4,i)=cz(i)
   20 nsho(i)=0
c
      do 21 i=1,ngrps
      l=icentr(i)
   21 nsho(l)=nsho(l)+1
      ishell=0
      ipr=0
      iij=0
      do 22 i=1,non
      do 23 j=1,ngrps
      if(icentr(j).ne.i) go to 23
      ishell=ishell+1
      m=itype(j)
      iik=ngrpp(m)
      iil=nbase(j)
      do 25 l=1,iik
      iij=iij+1
      new(iil+l)=iij
 25   newei(iij)=iil+l
      nt1=ntrm(j)
      do 24 k=1,nt1
      ipr=ipr+1
      if(ipr.gt.mxgaus) go to 6000
      ishp(ipr)=ishell
      ityp(ipr)=m
      igp(ipr)=ipr
      coef=anorm(k,j)*pip
      zes=zet(k,j)
      zet2=zes+zes
      cc=(pi/zet2)**1.5d0
      go to(800,801,802,803), m
 801  cc=cc/(zet2+zet2)
      go to 800
 802  cc=cc/(4.0d0*zet2*zet2)
      go to 800
 803  call caserr('cartesian f-functions not valid in atmol')
 800  coef=coef*dsqrt(cc)
      exxp(1,ipr)=zes
      exxp(2,ipr)=coef
   24 continue
   23 continue
   22 continue
c
      return
      end
      subroutine gamrd(iblkv)
      implicit real*8  (a-h,o-z),integer (i-n)
      character *8 title,tagg
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/junk/ishp(mxgaus),ityp(mxgaus),igp(mxgaus),exxp(2,mxgaus)
      common/lsort/atomo(4,maxat),poto,nsho(maxat)
      common/junkc/title(10),tagg(maxat)
      common/blkin/h(2)
      common/scra/ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     * cf(mxprim),cg(mxprim),z(maxat),
     * kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     * kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell,non,norb,nspace
      common/linkmr/new(255),newei(255),iprim,lbuff,igame
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/restri/nfils(63),lds(508),isect(508),ldsect(508)
      data m15/15/
c
      if(oprint(31))write(iwr,1)
 1    format(
     *' section 191 created by gamess program *****'/)
      igame=1
      m110=10+maxat
      m1420=mxprim*5+maxat
      call rdedx(ex,mxprim,iblkv,idaf)
      call reads(cs,m1420,idaf)
      call rdchrs(title,m110,idaf)
      nav = lenwrd()
      call readis(kstart,mach(2)*nav,idaf)
c
      if(non.le.0.or.non.gt.maxat)go to 1000
      if(nshell.le.0.or.nshell.gt.mxshel)go to 1000
      if(norb.gt.0.and.norb.le.maxorb)go to 1010
 1000 call caserr(
     *'parameter error detected in adapt pre-processor')
 1010 nat=non
      mach4 = 12*nat + 3 + 4/lenwrd()
      call secget(isect(493),m15,iblkv)
      call rdedx(h,mach4,iblkv,idaf)
c
      itemp = 9*nat + 1
      do 20 i=1,non
      atomo(1,i)=z(i)
      atomo(2,i)=h(itemp )
      atomo(3,i)=h(itemp + 1)
      atomo(4,i)=h(itemp + 2)
      nsho(i)=0
 20   itemp = itemp + 3
c
      iprim=0
      num=0
      ishell=0
c
      do 22 iat=1,non
      do 23 ii=1,nshell
      i=katom(ii)
      if(i.ne.iat)go to 23
      is=kstart(ii)
      ipri=kng(ii)
      if=is+ipri-1
      mini=kmin(ii)
      maxi=kmax(ii)
      kk=kloc(ii)-mini
c
      do 25 iorb=mini,maxi
      li=kk+iorb
      num=num+1
      new(num)=li
 25   newei(li)=num
c
      ishell=ishell+1
      j=maxi-mini
      if(j.eq.0)go to 26
      if(j-3)27,26,28
c ----- s
 26   m=1
      do 2008 ig=is,if
      iprim=iprim+1
      ishp(iprim)=ishell
      ityp(iprim)=m
      igp(iprim)=iprim
      exxp(1,iprim)=ex(ig)
 2008 exxp(2,iprim)=cs(ig)
      if(j.eq.0)go to 24
c ----- sp
      ishell=ishell+1
      nsho(iat)=nsho(iat)+1
c ----- p
 27   m=2
      do 2009 ig=is,if
      iprim=iprim+1
      ishp(iprim)=ishell
      ityp(iprim)=m
      igp(iprim)=iprim
      exxp(1,iprim)=ex(ig)
 2009 exxp(2,iprim)=cp(ig)
      go to 24
c ----- d, f and g
 28   if(j.eq.5) then
      m=3
      do 2010 ig=is,if
      iprim=iprim+1
      ishp(iprim)=ishell
      ityp(iprim)=m
      igp(iprim)=iprim
      exxp(1,iprim)=ex(ig)
 2010 exxp(2,iprim)=cd(ig)
c
      else if(j.eq.9) then
c
      m=4
      do 2011 ig=is,if
      iprim=iprim+1
      ishp(iprim)=ishell
      ityp(iprim)=m
      igp(iprim)=iprim
      exxp(1,iprim)=ex(ig)
 2011 exxp(2,iprim)=cf(ig)
      else
c     g functions
      m = 5
      do 2012 ig=is,if
      iprim=iprim+1
      ishp(iprim)=ishell
      ityp(iprim)=m
      igp(iprim)=iprim
      exxp(1,iprim)=ex(ig)
 2012 exxp(2,iprim)=cg(ig)
c
      endif
c
 24   nsho(iat)=nsho(iat)+1
      if(iprim.gt.mxgaus)go to 1000
 23   continue
 22   continue
c
      if(ishell.le.0.or.ishell.gt.mxgrps)go to 1000
      if(num.ne.norb)go to 1000
      return
      end
      subroutine amrdm(t,sm,hm,sa,ha,bbuff,cbbuff,lword)
      implicit real*8 (a-h,o-z), integer (i-n)
      character *1 dash
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/junk/ixa(3400),
     *bb(2550),cbb(2550),x(255),y(255),q(255),
     *nnum(256),i2(256),i3(256),i4(256),nord(256),ilife(256),
     *nco(100),nrw(100),ityp(256),icon(256),
     *iord(256),ncont(100)
      common/lsort /natt(100),ise(7),ip(700),ibas(100),
     *npf(7),npz(7),nel(3),
     *nbon(100),nfun(768),npar(768),mpar(768)
      common/scra /sexp(1000),coe(1000),
     *cm(100),cx(100),cy(100),cz(100),
     *z(100),xx(100),yy(100),zz(100)
     *,mcomp(256),ksym(2040),nfil(256),
     *imix(3),ibuk(3),inuk(3),nzer(8),nblz(8),ncal(8),
     *nzil(8),nsog(256),msym(2040),kcomp(256),
     *nwil(8),ngal(8),kj(8)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,is,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/linkmr/neway(255),newya(255)
      dimension guf(mxcrec),c(400)
      dimension h(7401)
      dimension t(*),sm(*),hm(*),sa(*),ha(*)
      dimension bbuff(*),cbbuff(*)
      equivalence (bb(1),h(1)),(c(1),cm(1))
c
      data dash/'-'/
      data mx2550 /2550/
      jtape=iwr
      itape=ird
      ntape=nf2
      ltape=nf22
      ipup=ltape
      mscr=nf1
      mtape=nf3
      jpup=mtape
      nrac=10
      cpu=cpulft(1)
      if(.not.oprint(29)) write(jtape,8999)cpu
 8999 format(/' commence symmetry adaption at ',f8.2,
     *' secs.')
      norb=255
      namx=100
      nsmx=700
      ifrk=500
      ifrk1=ifrk+1
      if2=1000
      ik=(lword-500)/3
      ik=ik+ik
      irs=(3*ik+500)/2
c
      read(ntape) natoms,nbas,nhfunc,(iord(i),i=1,
     1 nhfunc),repel,newya
      read(ntape) (iord(i),i=1,natoms),
     1(z(i),i=1,natoms),(xx(i),i=1,natoms),(yy(i),i=1,natoms),
     2(zz(i),i=1,natoms)
      nc2=2*natoms*(natoms-1)
      nc3=nc2*(natoms-2)/2
      nrec2=nc3*(natoms-3)/8+nc3+nc2+natoms
      nbasr=2
      i1=0
      i7=0
      do 110 i=1,nbas
      read(ntape)nraw,ncon
      nbasr=nbasr+1
      if(nraw.eq.0.or.ncon.eq.0)go to 111
      read(ntape)(sexp(i1+j),j=1,nraw),(coe(i1+j),j=1,nraw),
     1 (ityp(i7+j),j=1,ncon),(icon(i7+j),j=1,ncon)
      do j = i7+1, i7+ncon
        if (icon(j).gt.mxprms) then
          write(*,*)'*** nbas   = ',nbas
          write(*,*)'*** icon   = ',(icon(k),k=1,i7+ncon)
          write(*,*)'*** mxprms = ',mxprms
          write(*,*)'*** The number of primitive functions in the '
          write(*,*)'*** contraction <icon> exceeds parameter <mxprms>'
          call caserr(
     $       "Number of prim. functions <icon> exceeds <mxprms>")
        endif
      enddo
      nbasr=nbasr+1
 111  nco(i)=ncon
      nrw(i)=nraw
      i1=i1+nraw
      i7=i7+ncon
 110  continue
      iorbs=0
      do 114 k=1,natoms
      ncc=iord(k)
      iorbs=iorbs+nco(ncc)
      ncont(k)=nco(ncc)
 114  continue
      do 832 i=1,natoms
      natt(i)=0
      chg=z(i)
      ns1=ncont(i)
      if(dabs(chg).gt.1.0d-5)go to 831
      if(ns1.eq.0)natt(i)=3
      if(ns1.ne.0)natt(i)=1
      go to 832
 831  if(ns1.eq.0)natt(i)=2
 832  continue
      if(iorbs.gt.norb) call caserr(
     *'internal orbital count corrupted .. call for help .. 279')
      do 73 i=1,nbas
  73  nord(i)=0
      do 72 i=1,natoms
      isal=iord(i)
  72  nord(isal)=nord(isal)+1
      iexp=0
      ican=0
      ibat=0
      do 80 i=1,nbas
      isal=nord(i)
      ncw=nco(i)
      if(ncw.eq.0)go to 80
      do 82 j=1,ncw
      ican=ican+1
      ithn=icon(ican)
      ittp=ityp(ican)
      ibat=ibat+1
      nyx=nix(ittp)
      nyy=niy(ittp)
      nyz=niz(ittp)
      do 71 k=1,ithn
      iexp=iexp+1
      coff=coe(iexp)
      expf=sexp(iexp)
      jbat=ibat-ncw
      do 71 l=1,isal
      jbat=jbat+ncw
      jbati = (jbat-1)*mxprms + k
      bbuff(jbati)=expf
  71  cbbuff(jbati)=coff
      jbat=ibat-ncw
      do 82 k=1,isal
      jbat=jbat+ncw
      nnum(jbat)=ithn
      i2(jbat)=nyx
      i3(jbat)=nyy
  82  i4(jbat)=nyz
      ibat=jbat
 80   continue
c
c     ilife indexing array for bb and cbb
c
      jbat = 0
      do 770 i=1,iorbs
      ilife(i)=jbat
      jbat=jbat + nnum(i)
770   continue
      if(jbat.gt.mx2550) then
       call caserr('dimensioning problem in amrdm')
      endif
      kkk=0
      do 771 i=1,iorbs
      jbat=ilife(i)
      isal=(i-1)*mxprms
      jk=nnum(i)
      do 772 j=1,jk
      jbat = jbat + 1
      bb(jbat) = bbuff(isal+j)
      cbb(jbat) = cbbuff(isal+j)
772   continue
771   continue
c     
      jbat=0
      do 83 i=1,nbas
      isal=nord(i)
      ncv=nco(i)
      if(ncv.eq.0)go to 83
      jk=0
      do 408 k=1,natoms
      jp=iord(k)
      if (jp.ne.i) go to 408
      xt=xx(k)
      yt=yy(k)
      zt=zz(k)
      do 406 j=1,ncv
      jbat=jbat+1
      x(jbat)=xt
      y(jbat)=yt
 406  q(jbat)=zt
      jk=jk+1
      if (jk.eq.isal) go to 83
 408  continue
      call caserr(
     *'internal orbital count corrupted .. call for help .. 280')
  83  continue
      do 21 i=1,nbas
 21      ibas(i)=0
         ibs=0
         do 22 i=1,natoms
         k=iord(i)
         kkk=natt(i)
         if(kkk-1)22,2040,2041
2040     ibs=1
         go to 2042
2041     ibs=2
2042     ibas(k)=kkk
 22     continue
      do 7i=1,nsmx
7     ip(i)=0
      icx=-namx
      do 3 i=1,7
      icx=icx+namx
      ix=isx(i)
      iy=isy(i)
      iz=isz(i)
      do 9 j=1,nbas
      if (ibas(j).ne.0) go to 9
      isal=nord(j)
      isl=0
      idx=icx
      do 8 k=1,natoms
      idx=idx+1
      if(natt(k).ne.0.or.ip(idx).ne.0) go to 8
      if (iord(k).ne.j) go to 8
      nfl=0
      if (ix.eq.0) go to 10
      if (dabs(xx(k)).lt.1.0d-4) go to 11
      gx=-xx(k)
      nfl=1
      go to 12
11    gx=0.0d0
      go to 12
10    gx=xx(k)
12    if (iy.eq.0) go to 13
      if (dabs(yy(k)).lt.1.0d-4) go to 14
      gy=-yy(k)
      nfl=1
      go to 15
14    gy=0.0d0
      go to 15
13    gy=yy(k)
15    if(iz.eq.0) go to 16
      if (dabs(zz(k)).lt.1.0d-4) go to 17
      gz=-zz(k)
      nfl=1
      go to 18
17    gz=0.0d0
      go to 18
16    gz=zz(k)
18    if (nfl.eq.1) go to 19
      ip(idx)=k
      isl=isl+1
      if(isl.eq.isal) go to 9
      go to 8
19    k1=k+1
      if (k1.gt.natoms) go to 36
      iex=idx
      do 20 l=k1,natoms
      iex=iex+1
      if (natt(l).ne.0.or.ip(iex).ne.0) go to 20
      if (iord(l).ne.j) go to 20
      if (dabs(gx-xx(l)).gt.1.0d-4) go to 20
      if (dabs(gy-yy(l)).gt.1.0d-4) go to 20
      if (dabs(gz-zz(l)).gt.1.0d-4) go to 20
      ip(idx)=l
      ip(iex)=k
      isl=isl+2
      if(isl.eq.isal) go to 9
      go to 8
20    continue
36    ise(i)=0
      go to 3
8     continue
9     continue
      if(ibs.eq.0) go to 23
      do 24 j=1,nbas
      if (ibas(j).ne.1) go to 24
      isal=nord(j)
      isl=0
      idx=icx
      do 25 k=1,natoms
      idx=idx+1
      if (natt(k).ne.1.or.ip(idx).ne.0) go to 25
      if (iord(k).ne.j) go to 25
      nfl=0
      if (ix.eq.0) go to 26
      if (dabs(xx(k)).lt.1.0d-4) go to 27
      gx=-xx(k)
      nfl=1
      go to 28
27    gx=0.0d0
      go to 28
26    gx=xx(k)
28    if(iy.eq.0) go to 29
      if (dabs(yy(k)).lt.1.0d-4) go to 30
      gy=-yy(k)
      nfl=1
      go to 31
30    gy=0.0d0
      go to 31
29    gy=yy(k)
31    if(iz.eq.0) go to 6
      if(dabs(zz(k)).lt.1.0d-4) go to 32
      gz=-zz(k)
      nfl=1
      go to 33
32    gz=0.0d0
      go to 33
6     gz=zz(k)
33    if(nfl.eq.1) go to 34
      ip(idx)=k
      isl=isl+1
      if(isl.eq.isal) go to 24
      go to 25
34    k1=k+1
      if(k1.gt.natoms) go to 37
      iex=idx
      do 35 l=k1,natoms
      iex=iex+1
      if(natt(l).ne.1.or.ip(iex).ne.0) go to 35
      if (iord(l).ne.j) go to 35
      if(dabs(gx-xx(l)).gt.1.0d-4) go to 35
      if (dabs(gy-yy(l)).gt.1.0d-4) go to 35
      if (dabs(gz-zz(l)).gt.1.0d-4) go to 35
      ip(idx)=l
      ip(iex)=k
      isl=isl+2
      if(isl.eq.isal) go to 24
      go to 25
35    continue
37    ise(i)=0
      go to 3
25    continue
24    continue
      do 700 j=1,nbas
      if (ibas(j).ne.2) go to 700
      isal=nord(j)
      isl=0
      idx=icx
      do 701 k=1,natoms
      idx=idx+1
      if (natt(k).ne.2.or.ip(idx).ne.0) go to 701
      if (iord(k).ne.j) go to 701
      nfl=0
      if (ix.eq.0) go to 702
      if (dabs(xx(k)).lt.1.0d-4) go to 703
      gx=-xx(k)
      nfl=1
      go to 704
703   gx=0.0d0
      go to 704
702   gx=xx(k)
704   if(iy.eq.0) go to 705
      if (dabs(yy(k)).lt.1.0d-4) go to 706
      gy=-yy(k)
      nfl=1
      go to 707
706   gy=0.0d0
      go to 707
705   gy=yy(k)
707   if(iz.eq.0) go to 708
      if(dabs(zz(k)).lt.1.0d-4) go to 709
      gz=-zz(k)
      nfl=1
      go to 710
709   gz=0.0d0
      go to 710
708   gz=zz(k)
710   if(nfl.eq.1) go to 711
      ip(idx)=k
      isl=isl+1
      if(isl.eq.isal) go to 700
      go to 701
711   k1=k+1
      if(k1.gt.natoms) go to 712
      iex=idx
      do 713 l=k1,natoms
      iex=iex+1
      if(natt(l).ne.2.or.ip(iex).ne.0) go to 713
      if (iord(l).ne.j) go to 713
      if(dabs(gx-xx(l)).gt.1.0d-4) go to 713
      if (dabs(gy-yy(l)).gt.1.0d-4) go to 713
      if (dabs(gz-zz(l)).gt.1.0d-4) go to 713
      ip(idx)=l
      ip(iex)=k
      isl=isl+2
      if(isl.eq.isal) go to 700
      go to 701
713   continue
712   call caserr(
     *'internal orbital count corrupted .. call for help .. 281')
701   continue
700   continue
23    ise(i)=1
3     continue
      icx=-namx
      ix=0
      do 38 i=1,7
      icx=icx+namx
      if (ise(i).eq.0) go to 38
      iy=0
      ix=ix+1
      idx=icx
      do 39 j=1,natoms
      idx=idx+1
      if (ip(idx).eq.j) go to 39
      iy=iy+ncont(j)
39    continue
      npf(i)=iy
38    continue
      isymx=ix
      if (ix.gt.0 .and. natoms.ne.1 ) go to 40
      if(.not.oprint(29)) write(jtape,41)
 41   format(/1x,74('=')/1x,
     *'no symmetry elements have been recognized in the current ',
     1'coordinate system' /1x,74('=')/)
      isymx=0
      ix = 0
      nsel=1
      nj(1)=iorbs
      nbal(1)=0
      ntil(1)=0
      do 4000 i=1,iorbs
      ncomp(i)=1
 4000 lsym(i)=i
      jsym(1)=1
      knu=0
      do 4001 i=1,natoms
      knu=knu+1
      cm(knu)=z(i)
      cx(knu)=xx(i)
      cy(knu)=yy(i)
 4001 cz(knu)=zz(i)
      write(ltape)natoms,nbas,nhfunc,
     *(iord(i),i=1,nhfunc),repel,newya
      write(ltape)nsel,nj,iorbs,knu,c,jsym
      write(ltape)lsym,ncomp,ntil,nbal
      do 4005 i=1,natoms
 4005 nord(i)=iord(i)
      write(ltape)repel,h
      go to 4002
c
c ----- open sort file
c
40    call setbfa
      if (ix.eq.1) jx=1
      if (ix.eq.3) jx=2
      if (ix.eq.7) jx=3
      nsel=ix+1
      kx=0
      do 43 i=1,7
43    npz(i)=0
      do 45 j=1,jx
      lx=10000
      do 44 i=1,7
      if (ise(i).eq.0) go to 44
      if (npz(i).eq.1) go to 44
      nf=npf(i)
      if (nf.ge.lx) go to 44
      lx=nf
      ji=i
      if (nf.eq.0) go to 46
44    continue
      npz(ji)=1
      nel(j)=ji
      go to 45
46    kx=kx+1
      npz(ji)=1
      nel(j)=ji
45    continue
      if (jx.lt.3) go to 47
      ia=nel(1)
      ib=nel(2)
      kk=nel(3)
      igx=isx(ia)+isx(ib)+isx(kk)
      if (igx.gt.1) igx=igx-2
      if (igx.ne.0) go to 47
      igx=isy(ia)+isy(ib)+isy(kk)
      if (igx.gt.1) igx=igx-2
      if (igx.ne.0) go to 47
      igx=isz(ia)+isz(ib)+isz(kk)
      if (igx.gt.1) igx=igx-2
      if (igx.ne.0) go to 47
      lx=10000
      if(npf(kk).eq.0) kx=kx-1
      do 48 i=1,7
      if (npz(i).eq.1) go to 48
      nel(3)=i
      go to 49
48    continue
49    if(npf(i).eq.0) kx=kx+1
47    nbon(1)=0
      if (natoms.eq.1) go to 50
      do 51 i=2,natoms
      il=i-1
51    nbon(i)=nbon(il)+ncont(il)
50    icx=-iorbs
      do 52 i=1,jx
      icx=icx+iorbs
      ji=nel(i)
      idx=(ji-1)*namx
      ix=isx(ji)
      iy=isy(ji)
      iz=isz(ji)
      do 53 j=1,natoms
53    ibas(j)=0
      ican=0
      do 54 j=1,nbas
      isal=nord(j)
      ncw=nco(j)
      if(ncw.eq.0)go to 54
      icbn=ican
      ist=0
      do 55 k=1,natoms
      if (iord(k).ne.j) go to 55
      if(ibas(k).ne.0) goto 55
      ipx=idx+k
      ipx=ip(ipx)
      ka=nbon(k)
       nb=ka+icx
      ican=icbn
      if (ipx.ne.k) go to 56
      ist=ist+1
      do 57 l=1,ncw
      ican=ican+1
      nb=nb+1
         ka=ka+1
      ittp=ityp(ican)
      isum=0
      if (ix.eq.1) isum=isum+nix(ittp)
      if (iy.eq.1) isum=isum+niy(ittp)
      if (iz.eq.1) isum=isum+niz(ittp)
      if (isum-(isum/2)*2.eq.0) go to 58
      nfun(nb)=-ka
      npar(nb)=1
      go to 57
58    nfun(nb)=ka
      npar(nb)=0
57    continue
      if (ist.eq.isal) go to 54
      go to 55
 56       im=nbon(ipx)
      mb=im+icx
      ibas(ipx)=1
      ist=ist+2
      do 59 l=1,ncw
      nb=nb+1
       ka=ka+1
       im=im+1
      mb=mb+1
      npar(nb)=-1
      npar(mb)=-1
      ican=ican+1
      ittp=ityp(ican)
      isum=0
      if(ix.eq.1) isum=isum+nix(ittp)
      if (iy.eq.1) isum=isum+niy(ittp)
      if(iz.eq.1) isum=isum+niz(ittp)
      if (isum-(isum/2)*2.eq.0) go to 74
      nfun(nb)=-im
      nfun(mb)=-ka
      go to 59
74    nfun(nb)=im
      nfun(mb)=ka
59    continue
      if (ist.eq.isal) go to 54
55    continue
54    continue
52    continue
      if (kx.eq.jx) go to 84
      ix=0
      iy=0
      lx=0
      do 60 i=1,iorbs
60    nfil(i)=0
      ly=-iorbs
      do 61 i=1,iorbs
      ly=ly+1
      if(nfil(i).ne.0) go to 61
      my=ly
      do 62 j=1,jx
      my=my+iorbs
62    imix(j)=npar(my)
      iz=0
      jz=0
      do 63 j=1,jx
      if (imix(j).ge.0) go to 64
      iz=iz+1
      ibuk(iz)=j
      go to 63
64    jz=jz+1
      inuk(jz)=j
63    continue
      if (iz.gt.0) go to 65
      lx=lx+1
      mcomp(lx)=1
      ix=ix+1
      ksym(ix)=i
      do 75 j=1,jx
      iy=iy+1
75    mpar(iy)=imix(j)
      go to 61
65    ih=ix
      ix=ix+1
      ksym(ix)=i
      k=ibuk(1)
      ig=iorbs*(k-1)+i
      ix=ix+1
      mq=nfun(ig)
      ksym(ix)=mq
      if (mq.lt.0) mq=-mq
      nfil(mq)=1
      if (iz.eq.1) go to 66
      jg=iorbs*(ibuk(2)-1)
      ix=ix+1
      mq=nfun(jg+i)
      ksym(ix)=mq
      if (mq.lt.0) mq=-mq
      nfil(mq)=1
      ix=ix+1
      ig=ksym(ih+2)
      if (ig.lt.0) go to 68
      mq=nfun(jg+ig)
      ksym(ix)=mq
      if(mq.lt.0)mq=-mq
      nfil(mq)=1
      go to 69
68    mq=-nfun(jg-ig)
      ksym(ix)=mq
      if (mq.lt.0) mq=-mq
      nfil(mq)=1
69    if (iz.eq.2) go to 67
      jg=iorbs+iorbs
      ix=ix+1
      mq=nfun(jg+i)
      ksym(ix)=mq
      if(mq.lt.0) mq=-mq
      nfil(mq)=1
      do 86 j=1,3
      ix=ix+1
      ig=ksym(ih+j)
      if(ig.lt.0) go to 85
      mq=nfun(jg+ig)
      ksym(ix)=mq
      if(mq.lt.0)mq=-mq
      go to 86
85    mq=-nfun(jg-ig)
      ksym(ix)=mq
      if(mq.lt.0) mq=-mq
86    nfil(mq)=1
      do 91 j=1,8
      lx=lx+1
91    mcomp(lx)=8
      do 92 j=1,24
      iy=iy+1
92    mpar(iy)=irrep(j)
      ix=ix-7
      do 93 j=1,7
      ix=ix+8
93    ksym(ix)=i
      nq=0
      ih=ih+1
      do 94 j=1,7
      ih=ih+1
      mq=ksym(ih)
      kh=ih
      do 94 k=1,7
      nq=nq+1
      kh=kh+8
      if (idh(nq).lt.0) go to 95
      ksym(kh)=mq
      go to 94
95    ksym(kh)=-mq
94    continue
      ix=kh
      go to 61
67    do 96 j=1,4
      lx=lx+1
96    mcomp(lx)=4
      ix=ix-3
      do 97 j=1,3
      ix=ix+4
97    ksym(ix)=i
      nq=0
      ih=ih+1
      do 98 j=1,3
      ih=ih+1
      mq=ksym(ih)
      kh=ih
      do 98 k=1,3
      nq=nq+1
      kh=kh+4
      if (jdh(nq).lt.0) go to 99
      ksym(kh)=mq
      go to 98
99    ksym(kh)=-mq
98    continue
      ix=kh
      if (jz.eq.0) go to 120
      lq=inuk(1)
      mq=imix(lq)
      ky=iy+lq-3
      do 121 j=1,4
      ky=ky+3
121   mpar(ky)=mq
      ky=iy+ibuk(1)
      mpar(ky)=0
      mpar(ky+6)=0
      mpar(ky+3)=1
      mpar(ky+9)=1
      ky=iy+ibuk(2)
      mpar(ky)=0
      mpar(ky+3)=0
      mpar(ky+6)=1
      mpar(ky+9)=1
      iy=iy+12
      go to 61
120   mpar(iy+1)=0
      mpar(iy+2)=0
      mpar(iy+3)=1
      mpar(iy+4)=0
      mpar(iy+5)=0
      mpar(iy+6)=1
      mpar(iy+7)=1
      mpar(iy+8)=1
      iy=iy+8
      go to 61
66    ksym(ix+1)=ksym(ix-1)
      ksym(ix+2)=-ksym(ix)
      ix=ix+2
      do 122 j=1,2
      lx=lx+1
122   mcomp(lx)=2
      if (jz.lt.2) go to 123
      lq=inuk(1)
      mq=imix(lq)
      ky=iy+lq-3
      do 124 j=1,2
      ky=ky+3
124   mpar(ky)=mq
      lq=inuk(2)
      mq=imix(lq)
      ky=iy+lq-3
      do 125 j=1,2
      ky=ky+3
125   mpar(ky)=mq
      lq=ibuk(1)+iy
      mpar(lq)=0
      mpar(lq+3)=1
      iy=iy+6
      go to 61
123   if(jz.eq.0) go to 126
      lq=inuk(1)
      mq=imix(lq)
      ky=iy+lq-2
      do 127 j=1,2
      ky=ky+2
127   mpar(ky)=mq
      lq=ibuk(1)+iy
      mpar(lq)=0
      mpar(lq+2)=1
      iy=iy+4
      go to 61
126   mpar(iy+1)=0
      mpar(iy+2)=1
      iy=iy+2
61    continue
      go to 130
84    do 131 i=1,nsel
131   mj(i)=0
      kx=jx
      ly=-iorbs
      do 132 i=1,iorbs
      ly=ly+1
      my=ly
      isum=1
      ny=1
      do 133 j=1,jx
      my=my+iorbs
      lz=npar(my)
      if (lz.eq.1) isum=isum+ny
133   ny=ny+ny
      nc=mj(isum)+1
      mj(isum)=nc
      nir(i)=isum
132   loc(i)=nc
      nzer(1)=0
      do 134 i=2,nsel
      i9=i-1
134   nzer(i)=nzer(i9)+mj(i9)
      jv=0
      do 145 i=1,nsel
      js=mj(i)
      nj(i)=js
      ntil(i)=jv
      nbal(i)=jv
145   jv=jv+js
      do 135 i=1,iorbs
      ni=nir(i)
      nz=nzer(ni)+1
      nzer(ni)=nz
135   lsym(nz)=i
      do 136 i=1,iorbs
136   ncomp(i)=1
      ix=iorbs
      go to 220
130   do 140 i=1,nsel
140   mj(i)=0
      if (kx.eq.0) go to 141
      ly=-iorbs
      do 142 i=1,iorbs
      ly=ly+1
      my=ly
      isum=1
      ny=1
      do 143 j=1,kx
      my=my+iorbs
      lg=npar(my)
      if (lg.eq.1) isum=isum+ny
143   ny=ny+ny
      nc=mj(isum)+1
      mj(isum)=nc
      nir(i)=isum
142   loc(i)=nc
      ltape=jpup
141   do 144 i=1,nsel
      nj(i)=0
144   nblz(i)=0
      i9=0
      do 146 i=1,iorbs
      ny=1
      isum=1
      do 147 j=1,jx
      i9=i9+1
      if (mpar(i9).eq.1) isum=isum+ny
147   ny=ny+ny
      nsog(i)=isum
      nj(isum)=nj(isum)+1
146   nblz(isum)=nblz(isum)+mcomp(i)
      nzil(1)=0
      ncal(1)=0
      do 148 i=2,nsel
      i9=i-1
      ncal(i)=ncal(i9)+nblz(i9)
148   nzil(i)=nzil(i9)+nj(i9)
      nork=ncal(nsel)+nblz(nsel)
      do 149 i=1,nsel
      ntil(i)=nzil(i)
149   nbal(i)=ncal(i)
      lz=0
      nbal(nsel+1)=nork
      do 150 i=1,iorbs
      isum=nsog(i)
      mc=mcomp(i)
      nz=nzil(isum)+1
      nt=ncal(isum)
      ncomp(nz)=mc
      nzil(isum)=nz
      do 151 j=1,mc
      lz=lz+1
      nt=nt+1
151   lsym(nt)=ksym(lz)
150   ncal(isum)=nt
220   continue
      knu=0
      do 1000 i=1,natoms
      knu=knu+1
      cm(knu)=z(i)
      cx(knu)=xx(i)
      cy(knu)=yy(i)
      cz(knu)=zz(i)
1000  continue
      nb=-1
      write(ltape)nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel
     * ,newya
      write(ltape)nsel,mj,nj,iorbs,knu,c,jsym,msym
      write (ltape)nir,loc,lsym,ncomp,ntil,nbal
      do 4006 i=1,natoms
 4006 nord(i)=iord(i)
      write (ltape) repel,h
4002  nrecy=4
      do 9980 ipass=1,3
      nrec1=nbasr
      call getshm(sa,ha,natoms,ncont,iky,g,900,ntape,nrec1,jtape)
      iz=0
      ih=0
      ig=0
      do 202 i=1,nsel
      mjx=nj(i)
      if (mjx.eq.0) go to 202
      jz=iz
      mh=ih
      do 203 j=1,mjx
      jz=jz+1
      kz=iz
      mc=ncomp(jz)
      nh=ih
      do 215 k=1,j
      kz=kz+1
      ig=ig+1
      nc=ncomp(kz)
      hum=0.0d0
      sum=0.0d0
      lh=mh
      do 204 jj=1,mc
      lh=lh+1
      md=lsym(lh)
      if (md.lt.0) go to 205
      isig=0
      go to 206
205   md=-md
      isig=1
206   me=iky(md)
      kh=nh
      do 204 kk=1,nc
      kh=kh+1
      nd=lsym(kh)
      if (nd.lt.0) go to 207
      if (nd.gt.md) go to 208
      jg=me+nd
210   sam=sa(jg)
      ham=ha(jg)
      if (isig.eq.0) go to 213
      sam=-sam
      ham=-ham
      go to 213
208   jg=iky(nd)+md
      go to 210
207   nd=-nd
      if (nd.gt.md) go to 211
      jg=me+nd
212   sam=sa(jg)
      ham=ha(jg)
      if (isig.eq.1) go to 213
      sam=-sam
      ham=-ham
      go to 213
211   jg=iky(nd)+md
      go to 212
213   sum=sum+sam
204   hum=hum+ham
      sm(ig)=sum
      hm(ig)=hum
215   nh=kh
203   mh=kh
      ih=kh
      iz=kz
202   continue
      nblk=(ig-1)/mxcrc2+1
      nrecy=nblk+nrecy
      i9=0
      do 216 i=1,nblk
      lg=mxcrc2
      do 217 j=1,mxcrc2
      lg=lg+1
      i9=i9+1
      guf(j)=hm(i9)
217   guf(lg)=sm(i9)
216   write(ltape)guf
9980  continue
c
      if(oprint(31))write(jtape,6500)(dash,i=1,129)
 6500 format(/1x,129a1)
      call amrdm2(nork,ipup,isymx,t)
      return
      end
      subroutine amrdm2(nork,ipup,isymx,t)
      implicit real*8 (a-h,o-z), integer (i-n)
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
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/junk/ixa(3400),
     *bb(2550),cbb(2550),x(255),y(255),q(255),
     *nnum(256),i2(256),i3(256),i4(256),nord(256)
     *,nco(100),nrw(100),ityp(256),icon(256),
     *iord(256),ncont(100)
      common/lsort /natt(100),ise(7),ip(700),ibas(100),
     *npf(7),npz(7),nel(3),
     *nbon(100),nfun(768),npar(768),mpar(768)
      common/scra /sexp(1000),coe(1000),
     *cm(100),cx(100),cy(100),cz(100),
     *z(100),xx(100),yy(100),zz(100)
     *,mcomp(256),ksym(2040),nfil(256),
     *imix(3),ibuk(3),inuk(3),nzer(8),nblz(8),ncal(8),
     *nzil(8),nsog(256),msym(2040),kcomp(256),
     *nwil(8),ngal(8),kj(8)
c
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,is,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/linkmr/neway(255),newya(255)
      dimension guf(mxcrec),c(400),iwrit(50)
      dimension h(7401)
      dimension t(*)
      equivalence (bb(1),h(1)),(c(1),cm(1))
      if(isymx.gt.0)go to 4003
      call amrd00
      go to 8997
4003  if(oprint(31))write(jtape,6501)
 6501 format(/40x,
     *'symmetry adapted basis functions'/40x,32('-'))
      if (kx.eq.0) go to 300
      call amrd1(t(1),t(1),t(1))
      if(jx.ne.kx) go to 400
       n=0
       nc=1
       k=0
       do 861 i=1,nsel
      if(oprint(31))write(jtape,6502)i
 6502 format(/
     *' irreducible representation no. ',i2/1x,33('-'))
        njx=nj(i)
       if(njx.eq.0) go to 6503
       n=n+1
      if(oprint(31))write(jtape,6504)
 6504 format(/
     *' sequence     no. of    lcbf'/
     *' no. of sabf  terms    (atmol numbering)'/)
       do 868 j=1,njx
       k=k+1
       kkkk=lsym(k)
       l=newya(iabs(kkkk))
       if(kkkk.lt.0)l=-l
 868    if(oprint(31))write(jtape,6505)j,nc,l
 6505 format(4x,i4,8x,i2,7x,8i5)
       nj(n)=njx
       irper(n)=i
      if(oprint(31))write(jtape,6506)njx,n
 6506 format(/
     *' *** total no. of sabfs = ',i3/
     */'    revised representation no. = ',i3/)
      go to 861
 6503 if(oprint(31))write(jtape,6507)
 6507 format(/
     *' *** there are no sabf of this representation *** ')
  861  continue
       go to 8997
300   call amrd0(t(1),t(1))
      go to 600
400   ltape=ipup
      n=0
      do 401 i=1,nsel
      if(nj(i).ne.0) go to 402
      ircor(i)=0
      go to 401
402   nb=nbal(i)+1
      nb=lsym(nb)
      n=n+1
      if (nb.lt.0) nb=-nb
      ircor(i)=nir(nb)
401   continue
      nm=0
      do 403 i=1,ntel
      do 404 j=1,nsel
      if (ircor(j).ne.i) go to 404
      nm=nm+1
      irper(nm)=j
      if (nm.eq.n) go to 405
404   continue
403   continue
405   lx=0
      mx=0
      do 407 i=1,n
      nb=irper(i)
      nc=nbal(nb)
      nd=ntil(nb)
      mjx=nj(nb)
      kj(i)=mjx
      nwil(i)=mx
      ngal(i)=lx
      do 407 j=1,mjx
      nd=nd+1
      mx=mx+1
      kc=ncomp(nd)
      kcomp(mx)=kc
      do 407 k=1,kc
      lx=lx+1
      nc=nc+1
407   msym(lx)=lsym(nc)
      do 417 i=1,nork
      lx=msym(i)
      if(lx.lt.0) go to 418
      lsym(i)=loc(lx)
      go to 417
418   lsym(i)=-loc(-lx)
417   continue
      do 409 i=1,iorbs
409   ncomp(i)=kcomp(i)
      do 410 i=1,n
      nbal(i)=ngal(i)
      nj(i)=kj(i)
410   ntil(i)=nwil(i)
      nbal(n+1)=nork
      mx=0
      do 218 i=1,n
      nb=irper(i)
      mjx=nj(i)
      do 218 j=1,mjx
      mx=mx+1
      nir(mx)=nb
218   loc(mx)=j
      lx=0
      mx=0
      iw=1
       mi=1
      do 411 i=1,ntel
      mjx=mj(i)
      mtil(i)=mx
      mbal(i)=lx
      if (mjx.eq.0) go to 411
      mx=mx+mjx
      do 412 it=iw,n
      il=irper(it)
      if (ircor(il).eq.i) go to 412
      mi=it
      lx=nbal(it)
      go to 411
412   continue
411   iw=mi
      call rewftn(mtape)
      read(mtape)nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
          nb=-4
      write(ltape)nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
       read(mtape)nsel,nwil,nwil,iorbs,knu,c
      write(ltape)nsel,mj,nj,iorbs,knu,c,msym
      read(mtape)
      write(ltape)nir,loc,lsym,ncomp,ntil,nbal
      read(mtape) repel,h
      write(ltape) repel,h
      nblk=nrecy-4
      do 413 i=1,nblk
      read(mtape) guf
413   write(ltape) guf
      read(mtape) jx
      write(ltape) jx
      do 420 ii=1,jx
      read(mtape) i,j,k,l
      write(ltape) i,j,k,l
      if (i.ne.j) go to 421
      if (i.eq.k) go to 430
      go to 440
421   if (i.eq.k) go to 450
      go to 460
430   call ad1iii(i,t(1),t(1))
      go to 420
440   call ad1iik(i,k,t(1),t(1))
      go to 420
450   call ad1iji(i,j,t(1),t(1))
      go to 420
460   call ad1ijk(i,j,k,l,t(1),t(1))
420   continue
      do 470 i=1,nsel
470   inper(i)=0
      do 471 i=1,n
      j=irper(i)
471   inper(j)=i
      do 472 i=1,ntel
472   mver(i)=0
      do 473 i=1,n
      j=irper(i)
      j=ircor(j)
473   mver(j)=mver(j)+1
      mvom(1)=0
      do 474 i=2,ntel
474   mvom(i)=mvom(i-1)+mver(i-1)
      go to 800
600   n=0
      do 601 i=1,nsel
      if (nj(i).ne.0) go to 602
      ircor(i)=0
      go to 601
602   n=n+1
      ntil(n)=ntil(i)
      nbal(n)=nbal(i)
      nj(n)=nj(i)
      ircor(i)=1
601   continue
      do 603 i=1,n
603   irper(i)=i
      mvom(1)=0
      mver(1)=n
       mj(1)=iorbs
      do 604 i=1,nsel
604   inper(i)=0
      do 605 i=1,n
      j=irper(i)
605   inper(j)=i
      mtil(1)=0
      mbal(1)=0
      iu=0
      do 607 i=1,n
      nb=irper(i)
      nx=nj(i)
      do 607 j=1,nx
      iu=iu+1
      nir(iu)=nb
607   loc(iu)=j
800   call rewftn(ltape)
      call rewftn(mtape)
      do 650 i=1,n
      ij(i)=nj(i)-1
      ix=irper(i)
      iy=ircor(ix)
650   npal(i)=ntil(i)-mtil(iy)+1
      read (ltape) nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
      nb=-2
      write(mtape)nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
      read(ltape)nsel,nwil,nwil,iorbs,knu,c,msym
       if(kx.gt.0) go to 163
         do 164 i=1,nork
 164      msym(i)=lsym(i)
 163  write(mtape)nsel,n,ntel,mj,nj,iorbs,knu,c,nrecy,msym
      read(ltape)
      write(mtape)nir,loc,mver,mvom,lsym,ncomp,nbal,ntil,mbal,mtil,irper
     1,inper,ircor,jsym,npal,ij,mvil
      read(ltape) repel,h
      write(mtape) repel,h
      nblk=nrecy-4
      do 801 i=1,nblk
      read(ltape)guf
801   write(mtape)guf
      read(ltape)jx
      write(mtape) jx
      do 820 ii=1,jx
      read(ltape)i,j,k,l
      if (i.ne.j) go to 821
      if (i.eq.k) go to 830
      go to 840
821   if (i.eq.k) go to 850
      go to 860
830   call ad2iii(i,t(1),t(1))
      go to 820
840   call ad2iik(i,k,t(1),t(1))
      go to 820
850   call ad2iji(i,j,t(1),t(1))
      go to 820
860   call ad2ijk(i,j,k,l,t(1),t(1))
820   continue
      call rewftn(mtape)
      call rewftn(ltape)
      read(mtape) nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
      nb=-3
      write(ltape)nb,nbas,nhfunc,(iord(i),i=1,nhfunc),repel,newya
      read(mtape)nsel,n,ntel,mj,nj,iorbs,knu,c,nrecy,msym
      write(ltape)nsel,n,ntel,mj,nj,iorbs,knu,c,nrecy,msym,irs,ifrk1
      read(mtape)
      write(ltape)nir,loc,mver,mvom,lsym,ncomp,nbal,ntil,mtil,irper,inp
     1er,ircor,jsym,npal,ij,mvil
      read(mtape) repel,h
      write(ltape) repel,h
      do 900 i=1,nblk
      read(mtape) guf
900   write(ltape) guf
      read(mtape) jx
      iz=0
      do 863 i=1,nsel
      iw=inper(i)
      if(iw.eq.0) go to 863
       iz=iz+1
       mj(iz)=nj(iw)
       irper(iz)=i
 863  continue
      iz=0
      do 864 i=1,nsel
      iw=inper(i)
      if(oprint(31))write(jtape,6502)i
      if(iw.eq.0)go to 6600
      kg=ntil(iw)
      iz=iz+1
      njx=mj(iz)
      it=nbal(iw)
      if(oprint(31))write(jtape,6504)
       do 865 j=1,njx
        kg=kg+1
        nc=ncomp(kg)
        do 869 k=1,nc
        it=it+1
 869    lsym(k)=msym(it)
       do 6601 k=1,nc
       kkkk=lsym(k)
       if(kkkk.lt.0)go to 6602
       iwrit(k)=newya(kkkk)
       go to 6601
 6602  iwrit(k)=-newya(-kkkk)
 6601  continue
 865   if(oprint(31))write(jtape,6505) j,nc,(iwrit(l),l=1,nc)
      if(oprint(31))write(jtape,6506)njx,iz
      go to 864
 6600 if(oprint(31))write(jtape,6507)
 864     continue
      call amrd2(t(1),t(1),t(1))
 8997 cpu=cpulft(1)
      if(.not.oprint(29)) write(jtape,8998)cpu
 8998 format(/
     *' end of symmetry adaption at ',f8.2,' secs.'/)
      return
      end
      subroutine ad2iii(ma,t,is)
      implicit real*8 (a-h,o-z), integer (i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5119)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(gin(1),gout(1)),
     2  (p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      im=mtil(ma)
      write(mtape)im,im,im,im
      mjx=mj(ma)
      kxl=mtil(ma)
      njx=mver(ma)
      nzl=mvom(ma)
      ina=iky(mjx+1)
      im=ik/ina
      if (im.gt.ina) im=ina
      il=(ina-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak= ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do  1  i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
    1 nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if (iq.eq.1.or.ma.eq.1) go to 700
      call rewftn(mscr)
  701 read (ltape) f
      write(mscr) f
      if (n2.ne.0) go to 701
      ntape=mscr
      call rewftn(ntape)
      go to 702
  700 ntape=ltape
c
c    when ma=1 , in either single- or multi-pass mode
c    no copy is taken i.e. have to rewind and span nrecy
c    records when in multi pass mode .
c
  702 do  2  i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
   11 read (ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
        j=intin(int4  )
      idx=iky(i)+j-ifm
      if (idx.lt.1.or.idx.gt.imax) go to 3
      k=intin(int4+3)
      iln=(idx-1)/im+1
      kz=idx-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=q(jj)
      l=intin(int4+2)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad2iii: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm)=(kz-1)*ina+iky(k)+l
      if=nfg(iln)+1
      if (if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
    3 int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
      idx=kxl+isb
      jdx=kxl+jsb
       i=isb
 200    i=i+1
      j=jsb
c     oijkl(1)=i
      idx=idx+1
      it=nir(idx)
      it=iky(it)
      kdx=jdx
 201      j=j+1
      if(imc.lt.im) goto 40
      if (ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c  31 t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
   40 lv=nzl
      kdx=kdx+1
      kt=nir(kdx)
      lh=inper(kt)
      kt=jsym(kt+it)
      k5=iky(kt)
      imc=imc+1
c     oijkl(2)=j
      do 80 km=1,njx
      lv=lv+1
      kd=npal(lv)
      if (kd.gt.i) go to 30
      kr=irper(lv)
      if (kr.gt.kt) go to 81
      kr=k5+kr
      go to 82
   81 kr=iky(kr)+kt
   82 kr=jsym(kr)
      inr=inper(kr)
      if(inr.eq.0.or.inr.gt.lv) go to 80
      lp=npal(inr)
      kw=lp+ij(inr)
      mmx=nbal(inr)
      mmm=ntil(inr)
      lx=nbal(lv)
      kk=ntil(lv)
       lw=kd+ij(lv)
      if (lw.gt.i) lw=i
      do 41 k=kd,lw
      if (k.lt.i) go to 880
      mlim=j
      go to 881
  880 mlim=k
  881 kk=kk+1
      mcp=ncomp(kk)
      do 39 ll=1,mjx
   39 v(ll)=0.0d0
c     oijkl(3)=k
      do 43 ll=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if (ld.lt.0) go to 44
      if (ld.eq.1) go to 850
      lda=ld-1
      iw=iwa+iky(ld)
      do 70 l=1,lda
      iw=iw+1
      v(l)=v(l)+t(iw)
   70 continue
  850 iw=iwa+ld
      do 851 l=ld,mjx
      iz=iw+iky(l)
      v(l)=v(l)+t(iz)
  851 continue
      go to 43
   44 ld=-ld
      if (ld.eq.1) go to 852
      lda=ld-1
      iw=iwa+iky(ld)
      do 853 l=1,lda
      iw=iw+1
      v(l)=v(l)-t(iw)
  853 continue
  852 iw=iwa+ld
      do 854 l=ld,mjx
      iz=iw+iky(l)
      v(l)=v(l)-t(iz)
  854 continue
   43  continue
      mc=kw
      if (mc.gt.mlim) mc=mlim
      mx=mmx
      mm=mmm
      do 55 l=lp,mc
      mm=mm+1
      ncp=ncomp(mm)
      sum=0.0d0
      do 53 ll=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if (ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
   54 sum=sum-v(-ld)
   53 continue
       if(dabs(sum).lt.1.0d-12) go to 55
c     oijkl(4)=l
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(i)+j
      intout(ig4+1)=i4096(k)+l
      if (ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2iii labels'
       write(6,92233) (pli(loop),loop=1,250)
92233  format(5(1x,z16))
      endif
c
      ig4=-1
      write (mtape) p
      ig=0
   55 continue
   41 continue
   80 continue
   30 iwa=iwa+ina
        if( j.lt.i) go to 201
      jsb=0
           jdx=kxl
      if(i.lt.mjx) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      if(oprint(31)) then
       write(6,*)' ad2iii labels'
       write(6,92233) (pli(loop),loop=1,250)
      endif
      write(mtape) p
      return
   50 isb=i-1
      jsb=j-1
      call rewftn(ntape)
      if (ma.gt.1) go to 762
      do 57 i=1,nrecy
   57 read (ntape)
  762 do 58 i=1,il
      nfg(i)=0
   58 nrc(i)=maxb
    2 continue
      return
      end
      subroutine ad2iik(ma,mb,t,is)
      implicit real*8  (a-h,o-z), integer (i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5119)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     2(gin(1),gout(1)),(p(1),pq(1))
     3,(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      mjx=mtil(ma)
      mjy=mtil(mb)
      write(mtape) mjx,mjx,mjy,mjy
      mjx=mj(mb)
      mjy=mj(ma)
      ina=iky(mjx+1)
      inb=iky(mjy+1)
      im=ik/ina
      kxl=mtil(ma)
      njx=mver(mb)
       nzl=mvom(mb)
      if(im.gt.inb) im=inb
      il=(inb-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
       iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
 1    nfg(i)=0
      ig=0
       isb=0
      ig4=-1
       jsb=0
      if(iq.eq.1) go to 700
      call rewftn(mscr)
 701  read(ltape) f
      write(mscr) f
      if(n2.ne.0) goto 701
       ntape=mscr
        call rewftn(ntape)
      goto 702
 700   ntape=ltape
 702  do 2 i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
       irec=0
       ilc=0
  11  read(ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      j=intin(int4  )
      idx=iky(i)+j-ifm
      if(idx.lt.1.or.idx.gt.imax) goto 3
      k=intin(int4+3)
      iln=(idx-1)/im+1
      kz=idx-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=q(jj)
      l=intin(int4+2)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad2iik: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm)=(kz-1)*ina+iky(k)+l
      if=nfg(iln)+1
      if (if.eq.ib) goto 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
       goto 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
 3    int4=int4+4
      goto 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
       imc=im
      idx=kxl+isb
      jdx=kxl+jsb
        i=isb
 200   i=i+1
      j=jsb
c     oijkl(1)=i
      idx=idx+1
      it=nir(idx)
      it=iky(it)
      kdx=jdx
 201    j=j+1
      if(imc.lt.im) goto 40
      if(ilc.eq.il) goto 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c31   t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
 40    lv=nzl
c     oijkl(2)=j
      kdx=kdx+1
      kt=nir(kdx)
      kt=jsym(kt+it)
      k5=iky(kt)
      imc=imc+1
      do 80 km=1,njx
      lv=lv+1
      kd=npal(lv)
      kr=irper(lv)
      if(kr.gt.kt) goto 81
       kr=k5+kr
       goto 82
 81   kr=iky(kr)+kt
  82   kr=jsym(kr)
      inr=inper(kr)
      if(inr.eq.0.or.inr.gt.lv) goto 80
      lp=npal(inr)
      kw=lp+ij(inr)
      mmx=nbal(inr)
      mmm=ntil(inr)
      lx=nbal(lv)
      kk=ntil(lv)
      lw=kd+ij(lv)
      do 41 k=kd,lw
       kk=kk+1
       mcp=ncomp(kk)
      do 39 ll=1,mjx
 39    v(ll)=0.0d0
c     oijkl(3)=k
      do 43 ll=1,mcp
       lx=lx+1
      ld=lsym(lx)
       if(ld.lt.0) go to 44
       if(ld.eq.1) goto 850
      lda=ld-1
      iw=iwa+iky(ld)
      do 70 l=1,lda
      iw=iw+1
      v(l)=v(l)+t(iw)
 70     continue
 850     iw=iwa+ld
       do 851 l=ld,mjx
      iz=iw+iky(l)
       v(l)=v(l)+t(iz)
 851     continue
         goto 43
 44   ld=-ld
         if(ld.eq.1) go to 852
       lda=ld-1
      iw=iwa+iky(ld)
      do 853 l=1,lda
      iw=iw+1
       v(l)=v(l)-t(iw)
 853       continue
 852     iw=iwa+ld
      do 854 l=ld,mjx
      iz=iw+iky(l)
      v(l)=v(l)-t(iz)
 854    continue
 43     continue
      mc=kw
      if(mc.gt.k) mc=k
      mx=mmx
      mm=mmm
      do 55 l=lp,mc
      mm=mm+1
      ncp=ncomp(mm)
      sum=0
      do 53 ll=1,ncp
      mx=mx+1
      ld=lsym(mx)
       if(ld.lt.0) goto 54
      sum=sum+v(ld)
       goto 53
 54      sum=sum-v(-ld)
 53     continue
        if(dabs(sum).lt.1.0d-12) go to 55
c     oijkl(4)=l
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(i)+j
      intout(ig4+1)=i4096(k)+l
      if(ig.lt.ifrk) goto 55
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2iik labels'
       write(6,92233) (pli(loop),loop=1,250)
92233  format(5(1x,z16))
      endif
c
      write(mtape) p
       ig=0
       ig4=-1
  55   continue
  41   continue
  80   continue
        iwa=iwa+ina
       if(j.lt.i) go to 201
       jsb=0
              jdx=kxl
       if(i.lt.mjy) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2iik labels'
       write(6,92233) (pli(loop),loop=1,250)
      endif
c
       write(mtape) p
        return
 50     isb=i-1
       jsb=j-1
      call rewftn(ntape)
       do 58 i=1,il
      nfg(i)=0
 58    nrc(i)=maxb
 2      continue
      return
      end
      subroutine ad2iji(ma,mb,t,is)
      implicit real*8  (a-h,o-z), integer (i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(gin(1),gout(1)),
     2(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      mjx=mtil(ma)
      mjy=mtil(mb)
      write(mtape) mjx,mjy,mjx,mjy
      mjx=mj(ma)
      mjy=mj(mb)
      kxl=mtil(ma)
      kyl=mtil(mb)
      ifm=-mjy
      do 900 i=1,mjx
      ifm=ifm + mjy
  900 jdeks(i) = ifm
      ina=mjx*mjy
      im=ik/ina
      if (im.gt.ina) im=ina
      njx=mver(ma)
      nzl=mvom(ma)
      il=(ina-1)/im + 1
      iq=(il-1)/if2+1
      il=(il-1)/iq + 1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix = -ib
      iy = - ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i) = iy
      mli(i) = ix
      nrc(i) = maxb
    1 nfg(i) = 0
      ig = 0
      ig4 = -1
      isb = 0
      jsb = 0
      if(iq.eq.1) go to 700
      call rewftn(mscr)
  701 read(ltape) f
      write(mscr) f
      if(n2.ne.0) go to 701
      ntape = mscr
      call rewftn(ntape)
      go to 702
  700 ntape = ltape
  702 do 2 i8 = 1,iq
      call setsto(2000,0,intin)
      ifm = ifm + imax
      irec = 0
      ilc = 0
   11 read (ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      j=intin(int4  )
      idx=jdeks(i) + j - ifm
      if(idx.lt.1.or.idx.gt.imax) go to 3
      k = intin(int4+3)
      iln = (idx-1)/im + 1
      kz = idx - (iln-1)*im
      mq=mlr(iln) + 1
      mm = mli(iln)  + 1
      t(mq) = q(jj)
      l=intin(int4+2)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad2iji: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm) = (kz-1)*ina + jdeks(k) + l
      if = nfg(iln) + 1
      if(if.eq.ib) go to 10
      nfg(iln) = if
      mlr(iln) = mq
      mli(iln) = mm
      go to 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
    3 int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
      idx = kxl + isb
      jdx = kyl + jsb
       i=isb
 200    i=i+1
       j=jsb
c     oijkl(1) = i
      idx = idx + 1
      it = nir(idx)
      i5 = iky(it)
      kdx = jdx
 201      j=j+1
      if(imc.lt.im) goto 40
      if(ilc.eq.il) go to 50
      ilc = ilc + 1
      imc=0
      iwa=0
      call vclr(t,1,icak)
c     do 31 kk=1,icak
c  31 t(kk)=0.0
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
   40 lv = nzl
cc    oijkl(2) = j
      kdx = kdx + 1
      kt = nir(kdx)
      lh = inper(kt)
      if(kt.gt.it) go to 940
      kt = jsym(i5+kt)
      go to 941
  940 kt = iky(kt) + it
      kt = jsym(kt)
  941 k5 = iky(kt)
      imc=imc+1
      do 80 km=1,njx
      lv = lv + 1
      kd = npal(lv)
      if(kd.gt.i) go to 30
      kr = irper(lv)
      if(kr.gt.kt) go to 81
      kr = k5 + kr
      go to 82
   81 kr = iky(kr) + kt
   82 kr = jsym(kr)
      inr = inper(kr)
      if (inr.eq.0) go to 80
      lp = npal(inr)
      kw = lp + ij(inr)
      mmx = nbal(inr)
      mmm = ntil(inr)
      lx = nbal(lv)
      kk = ntil(lv)
      lw = kd + ij(lv)
      if (lw.gt.i) lw = i
      do 41 k=kd,lw
      mlim = kw
      if(k.lt.i) go to 880
      if(kw.gt.j) mlim = j
  880 kk = kk + 1
      mcp = ncomp(kk)
      do 39 ll=1,mjy
   39 v(ll) = 0.0d0
c     oijkl(3) = k
      do 43 ll=1,mcp
      lx = lx + 1
      ld = lsym(lx)
      if(ld.lt.0) go to 44
      iw = iwa + jdeks(ld)
      do 70 l=1,mjy
      iw = iw + 1
      v(l) = v(l) + t(iw)
   70 continue
      go to 43
   44 iw = iwa + jdeks(-ld)
      do 71 l=1,mjy
      iw = iw + 1
      v(l) = v(l) - t(iw)
   71 continue
   43 continue
      mx = mmx
      mm = mmm
      do 55 l=lp,mlim
      mm = mm + 1
      ncp = ncomp(mm)
      sum = 0.0d0
      do 53 ll=1,ncp
      mx = mx + 1
      ld = lsym(mx)
      if(ld.lt.0) go to 54
      sum = sum + v(ld)
      go to 53
   54 sum = sum - v(-ld)
   53 continue
        if(dabs(sum).lt.1.0d-12) go to 55
c     oijkl(4) = l
      ig = ig + 1
      pq(ig) = sum
      ig4=ig4+2
      intout(ig4  )=i4096(i)+j
      intout(ig4+1)=i4096(k)+l
      if(ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2iji labels'
       write(6,92233) (pli(loop),loop=1,250)
92233  format(5(1x,z16))
      endif
c
      ig4=-1
      write (mtape) p
      ig = 0
   55 continue
   41 continue
   80 continue
   30 iwa = iwa + ina
       if(j.lt.mjy) go to 201
      jsb=0
           jdx=kyl
       if(i.lt.mjx) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2iji labels'
       write(6,92233) (pli(loop),loop=1,250)
      endif
c
      write (mtape) p
      return
   50 isb = i-1
      jsb = j-1
      call rewftn(ntape)
      do 58 i=1,il
      nfg(i) = 0
   58 nrc(i) = maxb
    2 continue
      return
      end
      subroutine amrd2(t,is,ixa)
      implicit real*8  (a-h,o-z), integer (i-n)
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
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkk,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,jss,kss,lss,lr,jr,nrec1
      dimension f(751)
      dimension q(500),fm(251),is(*),ixa(*),t(*)
      dimension nfil(666),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),
     1(fm(251),n2),(g(1),f(1)),(gin(1),gout(1))
      integer *1 oijkl,oic,ojs,oks,ols
      integer *4 ijkl,ic,js,ks,ls
      dimension oijkl(4),oic(4),ojs(4),oks(4),ols(4)
      equivalence (oijkl(1),ijkl), (oic(1),ic)
      equivalence (ojs(1),js),(oks(1),ks),(ols(1),ls)
      data ijkl,ic,js,ks,ls /5*0/
      data maxb/9999999/
      jv=0
      ix=0
      il=0
      nmel=iky(nsel+1)
      nmel=iky (nmel+1)
      do 81 i=1,nmel
81    nfil(i)=0
      call setsto(2000,0,intin)
c
      nav = lenwrd()
c
      do 2 i=1,nsel
      mi=inper(i)
      do 2 j=1,i
      il=il+1
      ist=jsym(il)
      jj=inper(j)
      if (mi.eq.0)jj=0
      jl=0
      do 2 k=1,i
      mk=inper(k)
      if(jj.eq.0) mk=0
      mlim=k
      if (i.eq.k) mlim=j
      do 2 l=1,mlim
      ix=ix+1
      jl=jl+1
      if (mk.eq.0)  goto2
      if (inper(l).eq.0) goto 2
      if (ist.ne.jsym(jl)) goto 2
      jv=jv+1
      nfil(ix)=jv
2     continue
      ib=(3*jv)/2
      ib=(irs/ib)-1
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
        write(ltape) jv,ib
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      do 1 i=1,jv
       ix=ix+itx
       iy=iy+ity
       mlr(i)=iy
       mli(i)=ix
       nfg(i)=0
1     nrc(i)=maxb
      irec=0
      do 4 ii=1,jx
      read  (mtape) ia,ja,ka,la
6     read  (mtape) f
      call unpack(fm,8,intin,2000)
      int4=1
      do 5 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0) go to 4
      j=intin(int4  )
      l=intin(int4+2)
      k=intin(int4+3)
      i=i+ia
      j=j+ja
      k=k+ka
      l=l+la
      ir=nir(i)
      ic=loc(i)
      jr=nir(j)
      js=loc(j)
      kr=nir(k)
      ks=loc(k)
      lr=nir(l)
      ls=loc(l)
      if (jr.le.ir) goto 10
      nz=ic
      ic=js
      js=nz
      nz=ir
      ir=jr
      jr=nz
10    if (lr.le.kr) goto 11
      nz=ks
      ks=ls
      ls=nz
      nz=kr
      kr=lr
      lr=nz
11    idx=iky(ir) + jr
      jdx=iky(kr) + lr
      if (jdx-idx) 12,64,65
65    idx = iky(jdx) + idx
      nz=ic
      ic=ks
      ks=nz
      nz=js
      js=ls
      ls=nz
      goto 13
64    if(ic.ge.ks) go to 12
      nz=ic
      ic=ks
      ks=nz
      nz=js
      js=ls
      ls=nz
12    idx=iky(idx) + jdx
13    idx=nfil(idx)
      mq=mlr(idx)+1
      mm=mli(idx)+1
      if=nfg(idx)+1
      oijkl(2) = oic(1)
      oijkl(1) = ojs(1)
      oijkl(4) = oks(1)
      oijkl(3) = ols(1)
      t(mq)=q(jj)
      is(mm)=ijkl
      if( if.eq.ib) go to 14
      nfg(idx)=if
      mlr(idx)=mq
      mli(idx)=mm
      go  to 5
 14   call stopbk
      nwbnwb=ib
      lnklnk=nrc(idx)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(idx)=md
      mli(idx)=me
      nrc(idx)=irec
      nfg(idx)=0
      irec=irec+nsz
5     int4=int4+4
      goto 6
4     continue
      do 16 i=1,jv
      nz=nfg(i)
      if(nz.eq.0)go to 16
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(i)
      mq=mlr(i)-nz
      mm=mli(i)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      nrc(i)=irec
      irec=irec+nsz
 16   continue
c
      call stopbk
c
      idx=0
      do 20 ii=1,nsel
      do 20 jj=1,ii
      do 20 kk=1,ii
      mlim=kk
      if (ii.eq.kk) mlim=jj
      do 20 ll=1,mlim
      idx=idx+1
      jdx=nfil(idx)
      if (jdx.eq.0) goto 20
      lnklnk=nrc(jdx)
      write(ltape) ii,jj,kk,ll
      icl=0
      go to 22
 23   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      do 24 i=1,nwbnwb
      icl=icl+1
      q(icl)=gin(i)
      intout(icl)=ixa(i)
      if(icl.lt.ifrk)go to 24
      n2=icl
      call pack(fm,32,intout,500)
      write(ltape)f
      icl=0
 24   continue
 22   if(lnklnk.ne.maxb) go to 23
      intout(icl+1)=0
      call pack(fm,32,intout,500)
      n2=0
      write(ltape)f
 20   continue
      return
      end
      subroutine amrd1(t,is,ixa)
      implicit real*8  (a-h,o-z), integer (i-n)
      integer *2 icc
      integer *4 ijkl
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/blkin/g(510),mword
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/linkmr/new(255)
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      common/craypk/intin(1360)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/aplus/gsp(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jb,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkk,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension q(500),fm(251),is(*),ixa(*),t(*)
      dimension f(751),gin(5118)
      dimension nfil(666),intout(500)
      dimension icc(2)
      equivalence (ijkl,icc(1))
      equivalence (gout(1),gin(1)),(gsp(1),q(1)),
     1(gsp(501),fm(1)),
     2(fm(251),n2),(gsp(1),f(1)),(intin(1),intout(1))
c
      external fget
      data maxb/9999999/
      data ijkl/0/
c
      ntel=8
      if(kx.eq.2) ntel=4
      if(kx.eq.1) ntel=2
c
      nav = lenwrd()
c
      nmel=iky(ntel+1)
      nmel=iky(nmel+1)
      ix=0
      jx=0
      il=0
      do 1 i=1,nmel
    1  nfil(i)=0
      do 2 i=1,ntel
       mi=mj(i)
      do 2 j=1,i
      il=il+1
      ist=jsym(il)
      jj=mj(j)
      if(mi.eq.0) jj=0
      jl=0
      do 2 k=1,i
      mk=mj(k)
      if(jj.eq.0) mk=0
      mlim=k
      if(i.eq.k) mlim=j
      do 2 l=1,mlim
      ix=ix+1
      jl=jl+1
      if(mk.eq.0) go to 2
      if(mj(l).eq.0) go to 2
      if(ist.ne.jsym(jl)) go to 2
      jx=jx+1
      nfil(ix)=jx
 2    continue
      ib=(3*jx)/2
      ib=(irs/ib)-1
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      do 3 i=1,jx
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nfg(i)=0
 3    nrc(i)=maxb
      irec=0
      call setsto(1360,0,intin)
      do 4 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)120,4,120
 120   num=lotape(ifile)
      call search(liblk(ifile),num)
 710  call fget(g,mmmm,num)
      if(mmmm)720,4,720
720   call unpack(g(num2e+1),lab816,intin,numlab)
      int4=1
      do 5 int=1,mword
      j=new(intin(int4  ))
      i=new(intin(int4+1))
      l=new(intin(int4+2))
      k=new(intin(int4+3))
      gg=g(int)
      if(oprint(31)) then
       write(6,5566) int,i,j,k,l,gg
5566   format(1x,'amrd1: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      if(i.ge.j)go to 700
      mfg=i
      i=j
      j=mfg
 700  if(k.ge.l)go to 701
      mfg=k
      k=l
      l=mfg
 701  ijj=j+i4096(i)
      kl=l+i4096(k)
      if(ijj.ge.kl)go to 702
      mfg=k
      k=i
      i=mfg
      mfg=l
      l=j
      j=mfg
 702  ir=nir(i)
      ic=loc(i)
      jr=nir(j)
      js=loc(j)
      kr=nir(k)
      ks=loc(k)
      lr=nir(l)
      ls=loc(l)
      if(jr-ir)10,60,61
61    nz=ic
      ic=js
      js=nz
      nz=ir
      ir=jr
       jr=nz
      goto 10
 60   if(js.le.ic) goto 10
      nz=ic
      ic=js
      js=nz
10    if(lr-kr)11,62,63
63    nz=ks
      ks=ls
      ls=nz
      nz=kr
      kr=lr
      lr=nz
      goto 11
62    if(ls.le.ks) goto 11
      nz=ks
      ks=ls
      ls=nz
11    idx=iky(ir)+jr
      jdx=iky(kr)+lr
      if(jdx-idx) 12,64,65
 65   idx=iky(jdx)+idx
      nz=ic
      ic=ks
      ks=nz
      nz=js
      js=ls
      ls=nz
      go to 13
64    if(ic-ks) 67,68,12
67    nz=ic
      ic=ks
      ks=nz
      nz=js
      js=ls
      ls=nz
      goto 12
68    if(js.ge.ls) goto 12
      nz=js
      js=ls
      ls=nz
 12   idx=iky(idx)+jdx
 13   idx=nfil(idx)
      mq=mlr(idx)+1
      mm=mli(idx)+1
      if=nfg(idx)+1
      icc(1)=js+i4096(ic)
      icc(2)=ls+i4096(ks)
      t(mq)=gg
      is(mm)=ijkl
      if(if.eq.ib) go to 14
      nfg(idx)=if
      mlr(idx)=mq
      mli(idx)=mm
      go to 5
 14   call stopbk
      nwbnwb=ib
      lnklnk=nrc(idx)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(idx)=md
      mli(idx)=me
      nrc(idx)=irec
      nfg(idx)=0
      irec=irec+nsz
 5    int4=int4+4
      lbl=lbl+1
      if(lbl)710,4,710
 4    continue
      do 16 i=1,jx
      nz=nfg(i)
      if(nz.eq.0)go to 16
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(i)
      mq=mlr(i)-nz
      mm=mli(i)-nz
      call dcopy(nz,t(mq+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      nrc(i)=irec
      irec=irec+nsz
 16   continue
c
      call stopbk
c
      write(ltape) jx
      idx=0
      do 20 ii=1,ntel
      do 20 jj=1,ii
      do 20 kk=1,ii
      mlim=kk
      if(ii.eq.kk) mlim=jj
      do 20 ll=1,mlim
      idx=idx+1
       jdx=nfil(idx)
      if(jdx.eq.0) go to 20
      lnklnk=nrc(jdx)
      write(ltape) ii,jj,kk,ll
      icl=0
      go to 22
 23   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      do 24 i=1,nwbnwb
      icl=icl+1
      q(icl)=gin(i)
      intout(icl)=ixa(i)
      if(icl.lt.ifrk)go to 24
      n2=icl
      call pack(fm,32,intout,500)
      write(ltape)f
      icl=0
 24   continue
 22   if(lnklnk.ne.maxb) go to 23
      intout(icl+1)=0
      call pack(fm,32,intout,500)
      n2=0
      write(ltape)f
 20   continue
      return
      end
      subroutine getshm (s,h,natoms,ncon,if,buff,lbuff,oneu,nrec1,iwr)
      implicit real*8  (a-h,o-z),integer(i-n)
      logical odum
      integer oneu
      dimension buff(lbuff)
      dimension s(*),h(*),ncon(*),if(*)
      int=lbuff
      ipend=0
      do 19 ia=1,natoms
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 19
      do 18 ip=ipbeg,ipend
      it2=if(ip)
      do 17 iq=ipbeg,ip
      ix=it2+iq
      int=int+1
      if(int.le.lbuff) go to 16
      read(oneu,end=120,err=120) odum,buff
      nrec1=nrec1+1
      int=1
   16 s(ix)=buff(int)
   17 continue
   18 continue
   19 continue
      if(natoms.eq.1) go to 30
      ipend=ncon(1)
      do 29 ia=2,natoms
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 29
      iam1=ia-1
      iqend=0
      do 28 ib=1,iam1
      iqbeg=iqend+1
      iqend=iqend+ncon(ib)
      if(ncon(ib).eq.0)go to 28
      do 27 ip=ipbeg,ipend
      it2=if(ip)
      do 26 iq=iqbeg,iqend
      ix=it2+iq
      int=int+1
      if(int.le.lbuff) go to 25
      read(oneu,end=120,err=120) odum,buff
      nrec1=nrec1+1
      int=1
   25 s(ix)=buff(int)
   26 continue
   27 continue
   28 continue
   29 continue
   30 ipend=0
      do 39 ia=1,natoms
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 39
      do 38 ip=ipbeg,ipend
      it2=if(ip)
      do 37 iq=ipbeg,ip
      ix=it2+iq
      int=int+1
      if(int.le.lbuff) go to 36
      read(oneu,end=120,err=120) odum,buff
      nrec1=nrec1+1
      int=1
   36 h(ix)=buff(int)
   37 continue
   38 continue
   39 continue
      if(natoms.eq.1) go to 50
      ipend=ncon(1)
      do 49 ia=2,natoms
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      iam1=ia-1
      iqend=0
      if(ncon(ia).eq.0)go to 49
      do 48 ib=1,iam1
      iqbeg=iqend+1
      iqend=iqend+ncon(ib)
      if(ncon(ib).eq.0)go to 48
      do 47 ip=ipbeg,ipend
      it2=if(ip)
      do 46 iq=iqbeg,iqend
      ix=it2+iq
      int=int+1
      if(int.le.lbuff) go to 45
      read(oneu,end=120,err=120) odum,buff
      nrec1=nrec1+1
      int=1
   45 h(ix)=buff(int)
   46 continue
   47 continue
   48 continue
   49 continue
   50 if (int.eq.lbuff) read (oneu)
      if (int.eq.lbuff) nrec1=nrec1+1
      return
  120 write(iwr,130)
  130 format (//' ** error detected on reading one-electron integrals')
      call caserr(
     * 'error reading 1-electron integrals')
      return
      end
      subroutine readin(idxd,bufd,n,ni,lu)
      implicit real*8  (a-h,o-z),integer(i-n)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension bufd(n)
      dimension idxd(ni)
      read(lu,end=101,err=102) idxd,bufd
      return
  101 write(iwr,601) n,lu
      go to 103
  102 write(iwr,602) n,lu
 103  call caserr(
     * 'i/o error in reading integrals')
      return
  601 format(//' end of data trying to read',i6,' integrals; unit',i3)
  602 format(//' i/o error trying to read',i6,' integrals on unit',i3)
      end
      subroutine ad1iji(ma,mb,t,is)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/craypk/intin(2000),intout(2000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251),
     2is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(pli(251),m2),(gin(1),gout(1)),
     2(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      llx=mbal(ma)
      lly=mbal(mb)
      mjx=mj(ma)
      mjy=mj(mb)
      ifm=-mjy
      do 900 i=1,mjx
      ifm=ifm+mjy
  900 jdeks(i)=ifm
      kx=mtil(ma)
      ky=mtil(mb)
      ina=mjx*mjy
      im=ik/ina
      if (im.gt.ina) im=ina
      if (im.eq.0) call caserr(
     *'insufficient main memory for program to continue')
      il=(ina-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
    1 nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if (iq.eq.1) go to 700
      call rewftn(mscr)
  701 read (mtape) f
      write (mscr) f
      if (n2.ne.0) go to 701
      ntape=mscr
      call rewftn(ntape)
      go to 702
  700 ntape=mtape
  702 do 2 i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
   11 read (ntape) f
      call unpack(fm,8,intin,2000)
      int4=1
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      gg=q(jj)
      j=intin(int4  )
      l=intin(int4+2)
      k=intin(int4+3)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,gg
5566   format(1x,'ad1iji: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      idx=jdeks(i)+j
      kdx=jdeks(k)+l
      kt=kdx-ifm
      if (kt.lt.1.or.kt.gt.imax) go to 9
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+idx
      if=nfg(iln)+1
      if (if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 9
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
    9 if (idx.eq.kdx) go to 3
      kt=idx-ifm
      if (kt.lt.1.or.kt.gt.imax) go to 3
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+kdx
      if=nfg(iln)+1
      if (if.eq.ib) go to 12
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 3
 12   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5577 loop=1,ib
c5577 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
    3 int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
       k=isb
 200     k=k+1
       l=jsb
c     oijkl(3)=k
 201      l=l+1
      if(imc.lt.im) goto 40
      if (ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
      call vclr(t,1,icak)
c     do 31 kk=1,icak
c  31 t(kk)=0.0
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
   40 lx=llx
      kk=kx
      imc=imc+1
c     oijkl(4)=l
      do 41 ii=1,mjx
      kk=kk+1
      mcp=ncomp(kk)
      do 39 i=1,mjy
   39 v(i)=0.0d0
c     oijkl(1)=ii
      do 43 i=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if (ld.lt.0) go to 44
      iw=iwa+jdeks(ld)
      do 70 j=1,mjy
      iw=iw+1
      v(j)=v(j)+t(iw)
   70 continue
      go to 43
   44 iw=iwa+jdeks(-ld)
      do 71 j=1,mjy
      iw=iw+1
      v(j)=v(j)-t(iw)
   71 continue
   43 continue
      mx=lly
      mm=ky
      do 55 jj=1,mjy
      mm=mm+1
      ncp=ncomp(mm)
      sum=0
      do 53 i=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if (ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
   54 sum=sum-v(-ld)
   53 continue
      if (dabs(sum).lt.1.0d-13) go to 55
c     oijkl(2)=jj
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(ii)+jj
      intout(ig4+1)=i4096(k)+l
      if (ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
      m2=ig
      write (ltape) p
      ig=0
      ig4=-1
   55 continue
   41 continue
      iwa=iwa+ina
        if(l.lt.mjy) go to 201
      jsb=0
      if(k.lt.mjx) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      m2=0
      write (ltape) p
      return
   50 isb=k-1
      jsb=l-1
      call rewftn(ntape)
      do 58 i=1,il
      nfg(i)=0
   58 nrc(i)=maxb
    2 continue
       return
      end
       subroutine ad1ijk(le,lb,lc,lf,t,is)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(pli(251),m2),(gin(1),gout(1))
     2,(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      llx = mbal(le)
      lly=mbal(lb)
      mjx=mj(le)
      mjy=mj(lb)
      mja=mj(lc)
      mjb=mj(lf)
      kx=mtil(le)
      ky=mtil(lb)
      ina=mjx*mjy
      inb=mja*mjb
      ifm=-mjy
      do 901 i=1,mjx
      ifm=ifm+mjy
901   kdeks(i)=ifm
      ifm=-mjb
      do 900 i=1,mja
      ifm=ifm+mjb
900   jdeks(i)=ifm
      im=ik/ina
      if (im.gt.inb) im=inb
      if (im.eq.0) call caserr(
     *'insufficient main memory for program to continue')
      il=(inb-1)/im + 1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
1     nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if (iq .eq. 1 ) go to 700
      call rewftn(mscr)
701    read (mtape) f
      write (mscr) f
      if (n2 .ne. 0) goto 701
      ntape=mscr
      call rewftn(ntape)
      goto 702
700   ntape=mtape
702    do 2 i8=1 ,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
11    read (ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      k=intin(int4+3)
      l=intin(int4+2)
      kdx=jdeks(k) + l - ifm
      if (kdx .lt.1  .or. kdx .gt. imax) goto 3
      iln=(kdx-1)/im +1
      kz=kdx-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=q(jj)
      j=intin(int4  )
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad1ijk: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm)=(kz-1)*ina+kdeks(i)+j
      if=nfg(iln)+1
      if (if .eq. ib) goto 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      goto 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
3      int4=int4+4
      goto 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
       k=isb
 200      k=k+1
       l=jsb
c     oijkl(3)=k
 201     l=l+1
      if(imc.lt.im) goto 40
      if (ilc .eq. il) goto 50
      ilc=ilc+1
      imc=0
      iwa=0
      call vclr(t,1,icak)
c     do 31 kk=1,icak
c  31   t(kk)=0.0
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
40     lx=llx
      kk=kx
      imc=imc+1
c     oijkl(4)=l
      do 41 ii=1,mjx
      kk=kk+1
      mcp=ncomp(kk)
      do 39 i=1,mjy
39     v(i)=0.0d0
c     oijkl(1)=ii
      do 43 i=1,mcp
      lx=lx+1
      lt=lsym(lx)
      if (lt .lt. 0) goto 44
      iw=iwa+kdeks(lt)
      do 70 j=1,mjy
      iw=iw+1
      v(j)=v(j)+t(iw)
70    continue
      goto 43
44    iw=iwa+kdeks(-lt)
      do 71 j=1,mjy
      iw=iw+1
      v(j)=v(j)-t(iw)
71    continue
43    continue
      mx=lly
      mm=ky
      do 55 j=1,mjy
      mm=mm+1
      ncp=ncomp(mm)
      sum=0.0d0
      do 53  i=1,ncp
      mx=mx+1
      lt=lsym(mx)
      if (lt.lt.0) goto 54
      sum=sum+v(lt)
      goto 53
54    sum=sum-v(-lt)
53    continue
      if (dabs(sum) .lt. 1.0d-13) goto 55
c     oijkl(2)=j
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4)  =i4096(ii)+j
      intout(ig4+1)=i4096(k)+l
      if (ig .lt. ifrk) goto 55
      call pack(pli,16,intout,1000)
      m2=ig
      write (ltape) p
      ig=0
      ig4=-1
55    continue
41    continue
      iwa=iwa+ina
       if(l.lt.mjb) go to 201
      jsb=0
       if(k.lt.mja) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      m2=0
      write(ltape) p
      return
50    isb=k-1
      jsb=l-1
      call rewftn(ntape)
      do 58 i=1,il
      nfg(i)=0
58    nrc(i)=maxb
2     continue
      return
      end
      subroutine amrd0(t,is)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      common/blkin/g(510),mword
      common/linkmr/new(255)
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      common/aplus/gsp(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2400)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension pq(500),pli(251),is(*),gin(5118)
     1,t(*),intout(2000)
      equivalence (p(1),pq(1)),(p(501),pli(1)),(gin(1),gout(1)),
     2(pli(251),m2)
      equivalence(intout(1),intin(1))
      external fget
      data maxb/9999999/
c
      nav = lenwrd()
c
      ina=iky(iorbs+1)
      im=ik/ina
      if(im.gt.ina)im=ina
      if (im.eq.0) call caserr(
     *'insufficient main memory for program to continue')
      il=(ina-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      jx=1
      write(ltape) jx
      write(ltape)jx,jx,jx,jx
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
1     nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      do 2 i8=1,iq
      ifm=ifm+imax
      ilc=0
      irec=0
      call setsto(1360,0,intin)
      do 4 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)120,4,120
 120   num=lotape(ifile)
      call search(liblk(ifile),num)
 710  call fget(g,mmmm,num)
      if(oprint(31)) then
        write(6,*)'amrd0: mmmm,mword = ', mmmm,mword
      endif
      if(mmmm)720,4,720
720   call unpack(g(num2e+1),lab816,intin,numlab)
      int4=1
      do 5 int=1,mword
      gg=g(int)
      j=new(intin(int4  ))
      i=new(intin(int4+1))
      l=new(intin(int4+2))
      k=new(intin(int4+3))
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,gg
5566   format(1x,'amrd0: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      if(i.ge.j)go to 700
      mfg=i
      i=j
      j=mfg
 700  if(k.ge.l)go to 701
      mfg=k
      k=l
      l=mfg
 701  ijj=j+i4096(i)
      kl=l+i4096(k)
      if(ijj.ge.kl)go to 702
      mfg=k
      k=i
      i=mfg
      mfg=l
      l=j
      j=mfg
 702  idx=min(i,j)+iky(max(i,j))
      kdx=min(k,l)+iky(max(k,l))
      kt=kdx-ifm
      if (kt.lt.1.or.kt.gt.imax) go to 9
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+idx
      if=nfg(iln)+1
      if(if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 9
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
      call dcopy(ib,t(md+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
9     if (idx.eq.kdx) go to 5
      kt=idx-ifm
      if (kt.lt.1.or.kt  .gt.imax) go to 5
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+kdx
      if=nfg(iln)+1
      if (if.eq.ib) go to 12
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 5
 12   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
      call dcopy(ib,t(md+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
 5    int4=int4+4
      lbl=lbl+1
      if(lbl)710,4,710
 4    continue
      do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
      call dcopy(nz,t(mq+1),1,gout(1),1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
       k=isb
 200    k=k+1
      l=jsb
c     oijkl(3)=k
 201      l=l+1
      if(imc.lt.im) goto 40
      if (ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c  31 t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
40    lx=0
      imc=imc+1
c     oijkl(4)=l
      do 41 ii=1,iorbs
      mcp=ncomp(ii)
c     do 39 i=1,iorbs
c39   v(i)=0.0e0
      call vclr(v,1,iorbs)
c     oijkl(1)=ii
      do 43 i=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if(ld.lt.0) go to 44
      if(ld.eq.1) go to 850
      lda=ld-1
      iw=iwa+iky(ld)
      do 70 j=1,lda
      iw=iw+1
      v(j)=v(j)+t(iw)
70    continue
850   iw=iwa+ld
      do 851 j=ld,iorbs
      iz=iw+iky(j)
      v(j)=v(j)+t(iz)
851   continue
      go to 43
44    ld=-ld
      if (ld.eq.1) go to 852
      lda=ld-1
      iw=iwa+iky(ld)
      do 853 j=1,lda
      iw=iw+1
      v(j)=v(j)-t(iw)
853   continue
852   iw=iwa+ld
      do 854 j=ld,iorbs
      iz=iw+iky(j)
      v(j)=v(j)-t(iz)
854   continue
43    continue
      mx=0
      do 55 jj=1,ii
      ncp=ncomp(jj)
      sum=0.0d0
      do 53 i=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if(ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
54    sum=sum-v(-ld)
53    continue
      if(dabs(sum).lt.1.0d-13) go to 55
c     oijkl(2)=jj
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(ii)+jj
      intout(ig4+1)=i4096(k)+l
      if (ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
      m2=ig
      write(ltape) p
      ig=0
      ig4=-1
55    continue
41    continue
      iwa=iwa+ina
       if(l.lt.k) go to 201
      jsb=0
        if(k.lt.iorbs) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      m2=0
      write(ltape ) p
      return
50    isb=k-1
      jsb=l-1
      do 58 i=1,il
      nfg(i)=0
58    nrc(i)=maxb
2     continue
      return
      end
      subroutine ad1iii(ma,t,is)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
       dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(pli(251),m2),(gin(1),gout(1)),
     2(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      llx=mbal(ma)
      mjx=mj(ma)
      kx=mtil(ma)
      ina=iky(mjx+1)
      im=ik/ina
      if(im.gt.ina)im=ina
      if (im.eq.0) call caserr(
     *'insufficient main memory for program to continue')
      il=(ina-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
1     nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if(iq.eq.1.or.ma.eq.1) go to 700
      call rewftn(mscr)
701   read(mtape) f
      write(mscr) f
      if (n2.ne.0) go to 701
      ntape=mscr
      call rewftn(ntape)
      go to 702
700   ntape=mtape
702   do 2i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
11    read(ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      j=intin(int4  )
      l=intin(int4+2)
      k=intin(int4+3)
      gg=q(jj)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,gg
5566   format(1x,'ad1iii: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      idx=iky(i)+j
      kdx=iky(k)+l
      kt=kdx-ifm
      if (kt.lt.1.or.kt.gt.imax) go to 9
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+idx
      if=nfg(iln)+1
      if(if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 9
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
9     if(idx.eq.kdx) go to 3
      kt=idx-ifm
      if (kt.lt.1.or.kt.gt.imax) go to 3
      iln=(kt-1)/im+1
      kz=kt-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=gg
      is(mm)=(kz-1)*ina+kdx
      if =nfg(iln)+1
      if(if.eq.ib) go to 12
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 3
 12   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5577 loop=1,ib
c5577 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
3     int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
       k=isb
 200      k=k+1
      l=jsb
c     oijkl(3)=k
 201      l=l+1
       if(imc.lt.im) goto 40
      if(ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c  31 t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
40    lx=llx
      kk=kx
      imc=imc+1
c     oijkl(4)=l
      do 41 ii=1,mjx
      kk=kk+1
      mcp=ncomp(kk)
      do 39 i=1,mjx
39    v(i)=0.0d0
c     oijkl(1)=ii
      do 43 i=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if(ld.lt.0) go to 44
      if (ld.eq.1) go to 850
      lda=ld-1
      iw=iwa+iky(ld)
      do 70 j=1,lda
      iw=iw+1
      v(j)=v(j)+t(iw)
70    continue
850   iw=iwa+ld
      do 851 j=ld,mjx
      iz=iw+iky(j)
      v(j)=v(j)+t(iz)
851   continue
      go to 43
44    ld=-ld
      if(ld.eq.1) go to 852
      lda=ld-1
      iw=iwa+iky(ld)
      do 853 j=1,lda
      iw=iw+1
      v(j)=v(j)-t(iw)
853   continue
852   iw=iwa+ld
      do 854 j=ld,mjx
      iz=iw+iky(j)
      v(j)=v(j)-t(iz)
854   continue
43    continue
      mx=llx
      mm=kx
      do 55 jj=1,ii
      mm=mm+1
      ncp=ncomp(mm)
      sum=0.0d0
      do 53 i=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if(ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
54    sum=sum-v(-ld)
53    continue
      if(dabs(sum).lt.1.0d-13) go to 55
c     oijkl(2)=jj
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(ii)+jj
      intout(ig4+1)=i4096(k)+l
      if(ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
      m2=ig
      write(ltape)p
      ig=0
      ig4=-1
55    continue
41    continue
      iwa=iwa+ina
       if(l.lt.k) go to 201
        jsb=0
        if(k.lt.mjx) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      m2=0
      write(ltape) p
      return
50    isb=k-1
      jsb=l-1
      call rewftn(ntape)
      if(ma.gt.1) go to 762
      do 57 i=1,nrecy
57    read(ntape)
762   do 58 i=1,il
      nfg(i)=0
58    nrc(i)=maxb
2     continue
      return
      end
      subroutine ad1iik(ma,mb,t,is)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251)
     2,is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(pli(251),m2),(gin(1),gout(1)),
     2(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
       data maxb/9999999/
c
      nav = lenwrd()
c
      llx=mbal(ma)
      mjx=mj(ma)
      mjy=mj(mb)
      kx=mtil(ma)
      ina=iky(mjx+1)
      inb=iky(mjy+1)
      im=ik/ina
      if(im.gt.inb)im=inb
      if (im.eq.0) call caserr(
     *'insufficient main memory for program to continue')
      il=(inb-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak=ina*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
1     nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if(iq.eq.1) go to 700
      call rewftn(mscr)
701   read(mtape) f
      write(mscr) f
      if(n2.ne.0) go to 701
      ntape=mscr
      call rewftn(ntape)
      go to 702
700   ntape=mtape
702   do 2 i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
11    read(ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      j=intin(int4  )
      l=intin(int4+2)
      k=intin(int4+3)
      kdx=iky(k)+l-ifm
      if(kdx.lt.1.or.kdx.gt.imax) go to 3
      iln=(kdx-1)/im+1
      kz=kdx-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=q(jj)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad1iik: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm)=(kz-1)*ina+iky(i)+j
      if=nfg(iln)+1
      if(if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
3     int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
        k=isb
 200     k=k+1
      l=jsb
c     oijkl(3)=k
 201      l=l+1
      if(imc.lt.im) goto 40
      if (ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c31   t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
40    lx=llx
      imc=imc+1
c     oijkl(4)=l
      kk=kx
      do 41 ii=1,mjx
      kk=kk+1
      mcp=ncomp(kk)
      do 39 i=1,mjx
39    v(i)=0.0d0
c     oijkl(1)=ii
      do 43 i=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if (ld.lt.0) go to 44
      if (ld.eq.1) go to 850
      lda=ld-1
      iw=iwa+iky(ld)
      do 70 j=1,lda
      iw=iw+1
      v(j)=v(j)+t(iw)
70    continue
850   iw=iwa+ld
      do 851 j=ld,mjx
      iz=iw+iky(j)
      v(j)=v(j)+t(iz)
851   continue
      go to 43
44    ld=-ld
      if(ld.eq.1) go to 852
      lda=ld-1
      iw=iwa+iky(ld)
      do 853 j=1,lda
      iw=iw+1
      v(j)=v(j)-t(iw)
853   continue
852   iw=iwa+ld
      do 854 j=ld,mjx
      iz=iw+iky(j)
      v(j)=v(j)-t(iz)
854   continue
43    continue
      mx=llx
      mm=kx
      do 55 jj=1,ii
      mm=mm+1
      ncp=ncomp(mm)
      sum=0.0d0
      do 53 i=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if(ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
54    sum=sum-v(-ld)
53    continue
      if(dabs(sum).lt.1.0d-13) go to 55
c     oijkl(2)=jj
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(ii)+jj
      intout(ig4+1)=i4096(k)+l
      if(ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
      m2=ig
      write(ltape) p
      ig=0
      ig4=-1
55    continue
41    continue
      iwa=iwa+ina
      if(l.lt.k) go to 201
       jsb=0
        if(k.lt.mjy) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
      m2=0
      write(ltape) p
      return
50    isb=k-1
      jsb=l-1
      call rewftn(ntape)
      do 58 i=1,il
      nfg(i)=0
58    nrc(i)=maxb
2     continue
      return
      end
      subroutine readsh(s,h,ao,mst)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      common /blkin/ pp(4),y(512)
      common/linkmr/new(255)
c
      dimension s(*),h(*),ao(*),mst(2)
c
      data mmat,m2,m511/ 2,2,511/
c
      call secget(ionsec,m2,iblk33)
      call rdedx(pp,m511,iblk33,idaf)
      lenb = lensec(nx)
      do imat=1,mmat
      j=mst(imat)
      ibl = iblk33 + 1 + (j-1) * lenb
      call rdedx(ao,nx,ibl,idaf)
c
      if (j.eq.3.and.ozora) then
          if (nwcor_z.ne.nx) call caserr('zora confusion in readsh')
c....     no dynamic memory here ; mrdcim must have called zora
          call rdedx(h,nwcor_z,ibcor_z,num8)
          call vadd(h,1,ao,1,ao,1,nwcor_z)
          if (oscalz) then
c...  add scaling correction to 1-electron matrix
             call rdedx(h,nwcor_z,ibscalao_z,num8)
             call vadd(h,1,ao,1,ao,1,nwcor_z)
          end if
      end if
c
      ij = 0
        do ia = 1,num
         do ja = 1,ia
         ij = ij + 1
         yy=ao(ij)
         ih=new(ia)
         jh=new(ja)
         iinn=min(ih,jh)+iky(max(ih,jh))
         h(iinn)=yy
         if(imat.eq.1) s(iinn)=yy
         enddo
        enddo
      enddo
      return
      end
      subroutine face1(s,ns,iw)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      dimension s(*),mst(6)
      common /linkmr/ map(510),ipr,lbuff
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
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      data mst / 1,3,2,4,5,6 /
c
      i1=1
      i2=i1+nat
      i3=i2+nat
      i4=i3+3*mxgaus
      i5=i4+num
      i6=i5+num
      i7=i6+4*nat
      i8=i7+nat
      i9=i8+num
      i10=i9+num
      icx=i10+(num+1)/2
      if (icx.le.ns) go to 1000
1001  call caserr('insufficient memory in adapt pre-processor')
1000  call basis(s(i6) ,s(i2) ,s(i3) ,s(i4) ,s(i5), s(i1),
c                atom   nsh    exx    ntyp   nnum   ncon
     *           s(i7) ,s(i8) ,s(i9),nat,num,mxgaus,iw)
c                nraw   nord  ncent
      lbuff=900
      i1=i2
      i2=lbuff+i1
      i3=nx+i2
      i4=max(nx,na)+i3
      last=nx+i4
      if (last.gt.ns) go to 1001
      ik=1
      do 200 k=1,3
      call readsh(s(i2),s(i3),s(i4),mst(ik))
c                 s     h     dummy
      call putshm(s(i2),s(i3),s(i1),s(1),lbuff)
c                 s     h     buff   ncon
200   ik=ik+2
      call rewftn(nf2)
      return
      end
      subroutine wrt(i,j,a,k,l)
      implicit real*8  (a-h,o-z),integer(i-n)
      dimension a(k),j(l)
      write (i) j,a
      return
      end
      subroutine basis(atom,nsh,exx,ntyp,nnum,ncon,nraw,nord,nhelp,
     +                 natoms,ior,mxe,iw)
      implicit real*8  (a-h,o-z),integer(i-n)
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
      common/junk/ishp(mxgaus),ityp(mxgaus),igp(mxgaus),exxp(2,mxgaus)
      common/lsort/atomo(4,maxat),poto,nsho(maxat)
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
      common/linkmr/news(255),newei(255),ipr,lbuff,igame
c
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      common/junkc/title(10),tagg(maxat)
      dimension ioff(5),mtyp(15)
      dimension atom(4,natoms),nsh(natoms),exx(3,mxe),
     *ntyp(ior),nnum(ior),ncon(natoms),nraw(natoms),nhelp(ior),
     *nord(natoms)
      dimension mal(10),iotypa(35),icar(6),iotypg(35)
      character *8 title,tagg
      character *1 icar,dash,ishh
      data mal /1,3,6,10,15,1,4,6,10,15/
      data iotypa/1,2,3,4,6,7,9,5,8,10,
     *11,16,20,12,13,15,17,18,19,14,21,22,23,24,25,26,27,28,29,30,31,32,
     *33,34,35/
      data iotypg/ 1, 2, 3, 4, 5, 8,10, 6, 7, 9,
     1            11,16,20,12,13,15,17,18,19,14,
     2            21,22,23,24,25,26,27,28,29,30,31,32,33,34,35/
      data icar /'s','p','d','f','g','0'/,
     *ioff /0,1,4,10,20/
      data dash/'-'/
c
      potnuc=poto
      do 9999 i=1,nat
      nsh(i)=nsho(i)
      ncon(i)=0
      nraw(i)=0
      do 9999 j=1,4
 9999 atom(j,i)=atomo(j,i)
      iat=1
      norbs=0
      ksh=nsh(iat)
      if(ksh.eq.0)call caserr(
     *'indexing problem in adaptation 1 .. call for help ..259')
      iprim=0
      jsh=0
      kshel=1
      jorbs=0
      kprim=1
      iexx=1
      ioo=0
      i=0
2     continue
      i=i+1
      if(i.eq.ipr+1) then
        ig=ig+1
        ioo=-1
        goto 7
      endif 
      ish=ishp(i)
      ity=ityp(i)
      ig=igp(i)
      do 9997 j=1,2
 9997 exx(j,iexx)=exxp(j,ig)
      im=mal(ity)-1
      if (im.eq.0) goto 3
      do 500 j=1,im
      exx(1,iexx+j)=exx(1,iexx)
500   exx(2,iexx+j)=exx(2,iexx)
    3 iexx=iexx+im+1
      if (iexx.gt.mxgaus) call caserr(
     *'primitive allocation exceeded in adaptation')
      if (ish.eq.kshel) goto 503
      ioo=0
7     kshel=ish
      jsh=jsh+1
      do 600 j=1,kmal
      iprim=ig-kprim+iprim
      nnum(jorbs+j)=ig-kprim
600   ntyp(jorbs+j)=mtyp(j)
      kprim=ig
      jorbs=jorbs+kmal
      if (jsh.lt.ksh) goto 503
      nraw(iat)=iprim
      iprim=0
      ncon(iat)=jorbs-norbs
      norbs=jorbs
      iat=iat+1
      ksh=nsh(iat)
      jsh=0
503   if (ioo) 5,8,2
8     ioo=1
      if(ity .gt. 5) ity=ity-5
      kmal=im+1
      nhelp(jorbs+1)=ity
      if (im.eq.0) goto 11
      do 10 j=1,im
10    nhelp(jorbs+j+1)=6
11    kstrt=ioff(ity)
      do 703 j=1,kmal
      jjj=iotypa(kstrt+j)
      if(igame.eq.1)jjj=iotypg(kstrt+j)
 703  mtyp(j)=jjj
      goto 2
5     continue
      nbas=1
      nord(1)=1
      lnb=ncon(1)
      lpr=nraw(1)
      if (nat.eq.1) goto 99
      do 12 i=2,nat
      ig=0
      kmal=0
      ish=ncon(i)
      zi=atom(1,i)
      iprim=nraw(i)
      iat=i-1
      do 13 j=1,iat
      zij=zi-atom(1,j)
      jsh=ncon(j)
      kprim=nraw(j)
      if (jsh.ne.ish) goto 14
      if (kprim.ne.iprim) goto 14
      if(dabs(zij).gt.0.003d0)go to 14
      if(ish.eq.0)go to 2027
      do 16 k=1,ish
      if (nnum(lnb+k).ne.nnum(kmal+k)) goto 14
      if (ntyp(lnb+k).ne.ntyp(kmal+k)) goto 14
16    continue
      do 17 k=1,kprim
      do 15 l=1,2
      if (dabs(exx(l,lpr+k)-exx(l,ig+k)).gt.1.d-5) goto 14
15    continue
17    continue
2027  nord(i)=nord(j)
      goto 18
14    ig=ig+kprim
      kmal=kmal+jsh
13    continue
      nbas=nbas+1
      nord(i)=nbas
18    lnb=lnb+ish
12    lpr=lpr+iprim
99    if(.not.oprint(29)) write(iw,9168)title
 9168 format(//20x,90('*')/
     *20x,'*',88x,'*'/
     *20x,'*',88x,'*'/
     *20x,'*',4x,10a8,4x,'*'/
     *20x,'*',88x,'*'/
     *20x,90('*')/)
      if(oprint(31))write(iw,100) num,nat,
     *(i,(atom(j,i),j=1,4),nord(i),i=1,nat)
100   format(/
     *' no. of basis functions      ',i6/
     *' no. of nuclei               ',i6//
     *'  centre      charge     x          y          z       basis'/
     *2x,58('-')
     */(1x,i2,7x,4(f8.2,3x),i2))
      call rewftn(nf2)
      k=1
      write (nf2) nat,nbas,k,k,potnuc,newei
      write (nf2) (nord(i),i=1,nat),
     *((atom(i,j),j=1,nat),i=1,4)
      lnb=0
      lpr=0
      nbas=0
      if(oprint(31))write(iw,7001)(dash,i=1,129)
 7001 format(/1x,129a1)
      if(oprint(32))write(iw,7000)
 7000 format(/40x,28('-')/
     *40x,'basis function specification'/40x,28('-'))
      do 31 i=1,nat
      nco=ncon(i)
      nra=nraw(i)
      if (nord(i).le.nbas) goto 19
      jpr=lpr
      nbas=nbas+1
      if(oprint(32))write(iw,101) nbas,nco,nra
101   format(/5x,' **** atom type no.',i3,' no. of gtos: ',i3,
     *' no. of primitives: ',i3)
      if(nco.eq.0)go to 2037
      if(oprint(32))write(iw,104) (nnum(lnb+j),j=1,nco)
104   format(/' gto contractions: ',50i2)
      if(oprint(32))write(iw,6666)
6666  format(/)
      do 34 j=1,nco
      jsh=nnum(lnb+j)
      ksh=mal(nhelp(lnb+j))
      if (nhelp(lnb+j).eq.6) goto 34
      ishh=icar(nhelp(lnb+j))
      call trmrsp(exx(1,jpr+1),ksh,jsh)
      do 9 ity=1,jsh
9     if(oprint(32))write (iw,103) ishh,exx(2,ity+jpr),exx(1,ity+jpr)
      jpr=jpr+jsh*ksh
34    continue
103   format(1x,a1,3x,2f14.6)
 2037 write(nf2)nra,nco
      if(nra.ne.0.and.nco.ne.0)
     *write (nf2) ((exx(k,lpr+j),j=1,nra),k=1,2),
     *(ntyp(lnb+j),j=1,nco),(nnum(lnb+j),j=1,nco)
19    lpr=lpr+nra
31    lnb=lnb+nco
      return
      end
      subroutine putshm(s,h,buff,ncon,lbu)
      implicit real*8  (a-h,o-z), integer (i-n)
      logical otrue,ofalse
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
      common /linkmr/ map(510),ipr,lbuff
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
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
c
      dimension buff(lbu),s(*),h(*),ncon(*)
      otrue=.true.
      ofalse=.false.
      ir=0
      ipend=0
      do 19 ia=1,nat
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 19
      do 18 jp=ipbeg,ipend
      it2=iky(jp)
      do 17 jq=ipbeg,jp
      ix=it2+jq
      ir=ir+1
      if (ir .le. lbuff) goto 16
      ir=1
      write (nf2) ofalse,buff
   16 buff(ir)=s(ix)
   17 continue
   18 continue
   19 continue
      if(nat.eq.1) go to 30
      ipend=ncon(1)
      do 29 ia=2,nat
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 29
      iam1=ia-1
      iqend=0
      do 28 ib=1,iam1
      iqbeg=iqend+1
      iqend=iqend+ncon(ib)
      if(ncon(ib).eq.0)go to 28
      do 27 jp=ipbeg,ipend
      it2=iky(jp)
      do 26 jq=iqbeg,iqend
      ix=it2+jq
      ir=ir+1
      if (ir .le. lbuff) goto 25
      ir=1
      write (nf2) ofalse,buff
   25 buff(ir)=s(ix)
   26 continue
   27 continue
   28 continue
   29 continue
   30 ipend=0
      do 39 ia=1,nat
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 39
      do 38 jp=ipbeg,ipend
      it2=iky(jp)
      do 37 jq=ipbeg,jp
      ix=it2+jq
      ir=ir+1
      if (ir .le. lbuff) goto 36
      ir=1
      write (nf2) ofalse,buff
   36 buff(ir)=h(ix)
   37 continue
   38 continue
   39 continue
      if(nat.eq.1) goto 1
      ipend=ncon(1)
      do 49 ia=2,nat
      ipbeg=ipend+1
      ipend=ipend+ncon(ia)
      if(ncon(ia).eq.0)go to 49
      iam1=ia-1
      iqend=0
      do 48 ib=1,iam1
      iqbeg=iqend+1
      iqend=iqend+ncon(ib)
      if(ncon(ib).eq.0)go to 48
      do 47 jp=ipbeg,ipend
      it2=iky(jp)
      do 46 jq=iqbeg,iqend
      ix=it2+jq
      ir=ir+1
      if (ir .le. lbuff) goto 45
      ir=1
      write (nf2) ofalse,buff
   45 buff(ir)=h(ix)
   46 continue
   47 continue
   48 continue
   49 continue
1     write (nf2) otrue,buff
      return
      end
      subroutine trmrsp(e,i1,i2)
       implicit real*8  (a-h,o-z), integer (i-n)
      dimension e(3,2)
      if (i1.eq.1) goto 999
      if (i2.eq.1) goto 999
      nz=i1*i2
      m=1
      n=1
      do 1 i=1,nz
      e(3,i)=e(1,n)
      n=n+i1
      if (n.le.nz) goto 11
      m=m+1
      n=m
11    continue
1     continue
      do 2 i=1,nz
2     e(1,i)=e(3,i)
      m=1
      n=1
      do 3 i=1,nz
      e(3,i)=e(2,n)
      n=n+i1
      if (n.le.nz) goto 12
      m=m+1
      n=m
12    continue
3     continue
      do 4 i=1,nz
4     e(2,i)=e(3,i)
999   continue
      return
      end
      subroutine face2(sout,vout,s,tv,route)
      implicit real*8  (a-h,o-z), integer (i-n)
      logical route
      character *1 dash
      character *4 dump,file
      character *8 title
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
c ... fortran data streams 
c
      integer nf1, nf2, nf3, nf4, nf8, nf9, nf10, nf22
      integer nf31, nf32, nf33, nf34, nf35, nf36, nf38, nf39, nf42
      integer nf48, nf52, nf58, nf41
      common/ftape/nf1,nf2,nf3,nf4,nf8,nf9,nf10,nf22,
     +             nf31,nf32,nf33,nf34,nf35,nf36,nf38,nf39,nf42,
     +             nf48,nf52,nf58,nf41
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
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
      common/junk/e(7401)
      common/junkc/title(10)
      common /blkin/ nirr,mult(8,8),ir(maxorb),
     * ispaa(maxorb),ispbb(2,maxorb)
      common/craypk/intin(2000),intout(680),nword
      common /scra / potnuc,dx,dy,dz,ss(510),
     + coord(400),iord(256),msym(2048),
     + ncomp(256),nbal(9),ntil(8),mtil(4),irper(8),
     + inper(8),lsym(2040),
     + nj(8),nk(8),nst(256),nbi(8),nci(8)
      common/restri/nfils(63),lds(508),isect(508),ldsect(508)
      dimension guf(mxcrec),sm(mxcrc2),hm(mxcrc2)
      dimension fn(251),a(751)
      dimension sout(*),vout(*),s(*),tv(*)
      dimension gout(511),im(2)
c
      equivalence (e(1),guf(1)),
     +            (guf(1),hm(1)),(guf(mxcrc2+1),sm(1)),
     +            (gout(1),nirr),
     +            (guf(1),a(1)),(a(501),fn(1)),(fn(251),n2),
     +            (im(1),ss(1)),
     +            (gout(511),mword)
c
      data m51/51/
      data m0,m1,m2/0,1,2/
      data dump,file /'dump','main'/
      data dash/'-'/
c
      m511 = 511
      cpu=cpulft(1)
      if(.not.oprint(29)) write(iwr,7006)(dash,i=1,129)
      if(route)go to 666
      if(.not.oprint(29)) write(iwr,667)
      go to 668
 666  if(.not.oprint(29)) write(iwr,14)cpu
      if(.not.oprint(29)) write(iwr,1)dump,yed(idaf),ibl3d
      num2=m1tape(1)
      iblk2=m1blk(1)
      if(ibl3d) 4,4,5
    4 call caserr('invalid starting block for direct access file')
    5 if(.not.oprint(29)) write(iwr,1) file,yed(num2),iblk2
      if(iblk2) 4,4,7
c
c            now set up data from amrdm
c
  7   ntape=nf22
      ifrk=500
      call rewftn(ntape)
      read(ntape) nb,ndum,nhfunc,(iord(i),i=1,nhfunc),anucr
      if(nb.le.0)go to 6
      read(ntape)nsel,nj,iorbs,knu,coord,jsym
      read(ntape)lsym,ncomp,ntil,nbal
      read(ntape)anucr,e
      go to 66
 6    if(nb.eq.-1) go to 8
      read(ntape) nsel,nir,ndum,nk,nj,iorbs,knu,coord,ndum,lsym
      read(ntape) iord,nst,msym,ncomp,nbal,ntil,mtil,irper,inper
      read(ntape) anucr,e
      go to 9
c
    8 read(ntape) nsel,nj,nk,iorbs,knu,coord
      read(ntape) iord,lsym,ncomp,ntil,nbal
      read(ntape) anucr,e
 66   nir=0
      do 10 i=1,nsel
      ns=nj(i)
      if(ns.gt.0) go to 11
      inper(i)=0
      go to 10
   11 nir=nir+1
      inper(i)=nir
      nj(nir)=nj(i)
      ntil(nir)=ntil(i)
      nbal(nir)=nbal(i)
   10 continue
c
    9 if(oprint(31))write(iwr,12) title,knu,iorbs,nir
      if(oprint(31)) then
        write(iwr,8002)nsel,nir,iorbs,knu,iord
        write(iwr,8002)nj,nk,lsym,msym
        write(iwr,8002)ncomp,nbal,ntil,mtil
        write(iwr,8002)irper,inper
      endif
c
c ----- load up symmetry section for casscf etc
c
      nirr=nsel
      it=0
      do 80 i=1,8
      do 80 j=1,i
      it=it+1
      mult(i,j)=jsym(it)
      if(i.ne.j)mult(j,i)=mult(i,j)
 80   continue
c
      it=0
      do 81 i=1,nsel
      iw=inper(i)
      if(iw.eq.0)go to 81
      ns=nj(iw)
      if(ns.le.0)go to 81
      do 82 j=1,ns
      it=it+1
 82   ir(it)=i
 81   continue
c
      if(it.ne.iorbs)call caserr(
     *'indexing problem in adaptation .. call for help ..271')
c
      call secput(isect(490),m51,lensec(mach(13)),iblk33)
      call wrt3(gout,mach(13),iblk33,idaf)
      ig=0
      nums1=0
      it=0
      do 13 j=1,nsel
      iw=inper(j)
      if(iw.eq.0) go to 13
      it=it+1
      kg=ntil(iw)
      ns=nj(iw)
      nst(it)=ns
      nbi(it)=nums1
      nci(it)=ig
      nums1=nums1+iky(ns+1)
      ig=ig+ns
   13 continue
c
c  now restore 1e integrals from amrdm for storage on dumpfile
c  3 pass process to minimise store in amrdm
c    1.s,t+v ....2.t,x ....3.y,z
c  insert integrals into section ionsec on dumpfile
      nblk=(nums1-1)/mxcrc2+1
c
      call secget(ionsec,m2,iblk33)
c
c     decide if sabf 1e-ints section has been specified
c     create if necessary
c
      lenb = lensec(nx)
      if(ions2.ne.ionsec) then
       i = 1 * 6*lenb
       call secput(ions2,m2,i,iblk3a)
      else
       iblk3a=iblk33
      endif
c
c              input 1st block for potnuc
c              check against amrdm here
c
      call rdedx(potnuc,m511,iblk33,idaf)
      if(dabs(potnuc-anucr).gt.1.0d-5) call caserr(
     *'inconsistency detected in one-electron integral section')
c
      if(iblk33.ne.iblk3a) then
       call wrt3(potnuc,m511,iblk3a,idaf)
      endif
c
      ibl3a1 = iblk3a + 1
      do 22 ipass=1,3
      kg=0
      do 18 i=1,nblk
      read(ntape)guf
      do 18 j=1,mxcrc2
      kg=kg+1
      s(kg)=sm(j)
      tv(kg)=hm(j)
      if(kg.eq.nums1)go to 181
 18   continue
c
181   continue
      call vclr(sout,1,nx)
      call vclr(vout,1,nx)
      ig=0
      do 20 i=1,nir
      mjx=nst(i)
      nbase=nci(i)
      do 21 j=1,mjx
      inti=nbase+j
      iii=iky(inti)
      do 21 k=1,j
      intj=nbase+k
      ig=ig+1
      ij=intj+iii
      sout(ij)=s(ig)
   21 vout(ij)=tv(ig)
   20 continue
c
      if(ipass.eq.1) then
        ibl1 = ibl3a1 
        ibl2 = ibl3a1 + 2*lenb
      else if(ipass.eq.2) then
        ibl1 = ibl3a1 + lenb
        ibl2 = ibl3a1 + 3*lenb
      else
        ibl1 = ibl3a1 + 4*lenb
        ibl2 = ibl1 + lenb
      endif
c
c     now output ao to sec. isec3
c
      call wrt3(sout,nx,ibl1,idaf)
      call wrt3(vout,nx,ibl2,idaf)
c
 22   continue
      cpu=cpulft(1)
      if(.not.oprint(29)) write(iwr,33)cpu,ions2
c
c            now for the 2el integrals
c
      call search(iblk2,num2)
      nword=-1
      mword=0
      call setsto(2000,0,intin)
      ijint=0
c
      read(ntape) kg
      do 400 ix=1,kg
      read(ntape) i,j,k,ndum
c     write(iwr,8000)i,j,k,l,n2
      if(i.eq.j) go to 42
      if(i.eq.k) go to 41
  420 read(ntape) a
      if(n2.ne.0) go to 420
      go to 400
c...
c...  (ia jb/ka lb)
 41   iqx=nci(i)
      jqx=nci(j)
   53 read(ntape) a
      call unpack(fn,8,intin,2000)
      int4=1
      isg=0
   54 isg=isg+1
      inti=intin(int4+1)
      if(inti.eq.0)go to 400
      iq=inti+iqx
      jq=intin(int4  )+jqx
      lq=intin(int4+2)+iqx
      kq=intin(int4+3)+jqx
      ijint=ijint+1
      mword=mword+1
      nword=nword+2
      gout(mword)=a(isg)
      loop=min(iq,jq)+i4096(max(iq,jq))
      moop=min(kq,lq)+i4096(max(kq,lq))
      intout(nword  )=max(loop,moop)
      intout(nword+1)=min(loop,moop)
      if(mword.ne.340) go to 700
      call pack(gout(num2ep+1),lab1632,intout,numlabp)
      call put(gout(1),m511,num2)
      mword=0
      nword=-1
      iblk2=iblk2+1
700   int4=int4+4
      if(isg.lt.ifrk) go to 54
      go to 53
c
c...  (ia ja / kb lb)
c
 42   if(i.eq.k) go to 43
      jqx=nci(k)
      iqx=nci(i)
   56 read(ntape) a
      call unpack(fn,8,intin,2000)
      int4=1
      isg=0
   57 isg=isg+1
      inti=intin(int4+1)
      if(inti.eq.0)go to 400
      iq=inti+iqx
      jq=intin(int4  )+iqx
      lq=intin(int4+2)+jqx
      kq=intin(int4+3)+jqx
      ijint=ijint+1
      mword=mword+1
      nword=nword+2
      gout(mword)=a(isg)
      loop=min(iq,jq)+i4096(max(iq,jq))
      moop=min(kq,lq)+i4096(max(kq,lq))
      intout(nword  )=max(loop,moop)
      intout(nword+1)=min(loop,moop)
      if(mword.ne.340) go to 800
      call pack(gout(num2ep+1),lab1632,intout,numlabp)
      call put(gout(1),m511,num2)
      mword=0
      nword=-1
      iblk2=iblk2+1
800   int4=int4+4
      if(isg.lt.ifrk) go to 57
      go to 56
c
c     (iaja / kala)
c
 43   iqx=nci(i)
   51 read(ntape) a
      call unpack(fn,8,intin,2000)
      int4=1
      isg=0
   58 isg=isg+1
      inti=intin(int4+1)
      if(inti.eq.0)go to 400
      iq=inti    +iqx
      jq=intin(int4  )+iqx
      lq=intin(int4+2)+iqx
      kq=intin(int4+3)+iqx
      ijint=ijint+1
      mword=mword+1
      nword=nword+2
      gout(mword)=a(isg)
      loop=min(iq,jq)+i4096(max(iq,jq))
      moop=min(kq,lq)+i4096(max(kq,lq))
      intout(nword  )=max(loop,moop)
      intout(nword+1)=min(loop,moop)
      if(mword.ne.340) go to 900
      call pack(gout(num2ep+1),lab1632,intout,numlabp)
      call put(gout(1),m511,num2)
      mword=0
      nword=-1
      iblk2=iblk2+1
900   int4=int4+4
      if(isg.lt.ifrk) go to 58
      go to 51
c
 400  continue
c
c
      if(mword.lt.m1) go to 9999
      call pack(gout(num2ep+1),lab1632,intout,numlabp)
      call put(gout,m511,num2)
9999  call put(gout,m0,num2)
c
      cpu=cpulft(1)
      if(.not.oprint(29)) write(iwr,60)cpu,ijint
668   if(.not.oprint(29)) then
        write(iwr,61)
        write(iwr,7006)(dash,i=1,129)
       endif
c
      if(oprint(31))call secsum
      call revind
      if(oprint(31))call whtps
      call clredx
      return
 61   format(/5x,'****** symmetry adaption complete *****'/)
 60   format(' mainfile output complete at ',f8.2,' secs.'//
     *         ' number of 2-e integrals =  ',i10)
c8000 format(10x,'******',4i5,i10)
 7006 format(1x,129a1)
 667  format(/' *** no output of symmetry adapted integrals',
     *' to direct access files'/)
 14   format(/
     *' *** output of symmetry adapted integrals to direct access files'
     */' *** commencing at ',f8.2,' secs.')
 1    format(/ 5x,a4,'file on ',a4,' commencing at block',i5)
 12   format(/' case : ',10a8//
     *        ' no. of nuclei                    :',i3/
     *        ' no. of basis functions           :',i3/
     *        ' no. of irreducible reps          :',i3/)
 8002 format(5x,20i5)
   33 format(' dumpfile output complete at ',f8.2,' secs.'/
     *        ' 1-e integrals to section',i4)
c
      end
      subroutine ad2ijk(i,j,k,l,t,is)
      implicit real*8  (a-h,o-z), integer (i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/aplus/g(900),p(751),v(255),idxd(300),
     *nir(256),loc(256),ncomp(256),jdeks(256),kdeks(256),
     *lsym(2040),mj(8),nj(8),ntil(8),nbal(9),irper(8),
     *ircor(8),mtil(4),mbal(4),inper(8),mver(4),
     *mvom(4),npal(8),ij(8),mvil(4)
c
      integer irrep, nix, niy, niz, isx, isy, idh, jdh, isz, jsym
      common /blockc/ irrep(24),nix(35),niy(35),niz(35),isx(7),isy(7),
     +                idh(49),jdh(9),isz(7),jsym(36)
c
      common/three/mli(1000),mlr(1000),nrc(1000),nfg(1000)
      common/bufb/nwbnwb,lnklnk,gout(5118)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      common/stak/btri,mlow(2),irec
      common/junk/ixa(3400)
      common/craypk/intin(2000),intout(2000)
      common/b/ sum,gg,itape,jtape,ntape,mtape,ltape,mscr,nrac,nrecy,
     1nblk,ifrk,ifrk1,ik,if2,im,mjx,kxl,njx,nzl,ina,il,iq,ib,icak,
     2ib2,itx,nbyt,ity,ix,iy,imax,ifm,ig,isb,jsb,jx,ilc,jj,idx,
     3iln,kz,mq,mm,if,md,iba,nz,jdx,iwa,it,kdx,kt,lh,k5,jm,ixq,lv,km,
     4kd,kr,inr,lp,kw,mmx,mmm,lx,kk,lw,mlim,mcp,ll,ld,lda,iw,iz,mc,mx,
     5ncp,kx,nrec2,ntel,irs,ijkl,nmel,ist,mi,jl,mk,ibl,jn,in1,in,ncas,
     6nlop,ia,ja,ka,la,kp,lenth,int,ii,icl,iorbs,llx,mjy,inb,i8,lly,ky,
     7mja,mjb,lt,kyl,i5,nsel,n,jv,js,ks,ls,lr,jr,nrec1,imc
      dimension f(751)
      dimension q(500),fm(251),pq(500),pli(251),
     2is(*),t(*),gin(5118)
      equivalence (g(1),q(1)),(g(501),fm(1)),(fm(251),n2),
     1(gin(1),gout(1)),
     2(p(1),pq(1)),(p(501),pli(1))
     3,(g(1),f(1))
      data maxb/9999999/
c
      nav = lenwrd()
c
      mjx=mtil(i)
      mjy=mtil(j)
      kxl=mtil(k)
      kyl=mtil(l)
      write(mtape) mjx,mjy,kxl,kyl
      mjx=mj(i)
      mjy=mj(j)
      kxl=mtil(i)
      kyl=mtil(j)
      ifm=-mjy
c     novector
      do 901 ll=1,mjx
      ifm=ifm+mjy
  901 kdeks(ll)=ifm
      mja=mj(k)
      mjb=mj(l)
      ifm=-mjb
      do 900 ll=1,mja
      ifm=ifm+mjb
  900 jdeks(ll)=ifm
      njx=mver(k)
      nzl=mvom(k)
      ina=mjx*mjy
      inb=mja*mjb
      im=ik/inb
      if(im.gt.ina) im=ina
      il=(ina-1)/im+1
      iq=(il-1)/if2+1
      il=(il-1)/iq+1
      ib=ik/il
      icak= inb*im
      itx=ib/nav
      if(itx*nav.ne.ib) ib=ib-1
      if (ib.gt.nsz340) ib=nsz340
      itx=(nav+1)*ib
      ity=itx/nav
      ix=-ib
      iy=-ity
      imax=im*il
      ifm=-imax
      do 1 i=1,il
      ix=ix+itx
      iy=iy+ity
      mlr(i)=iy
      mli(i)=ix
      nrc(i)=maxb
    1 nfg(i)=0
      ig=0
      ig4=-1
      isb=0
      jsb=0
      if(iq.eq.1) go to 700
      call rewftn(mscr)
  701 read(ltape) f
      write(mscr) f
      if(n2.ne.0) go to 701
      ntape=mscr
      call rewftn(ntape)
      go to 702
  700 ntape=ltape
  702 do 2 i8=1,iq
      call setsto(2000,0,intin)
      ifm=ifm+imax
      irec=0
      ilc=0
  11  read(ntape) f
      int4=1
      call unpack(fm,8,intin,2000)
      do 3 jj=1,ifrk
      i=intin(int4+1)
      if(i.eq.0)go to 705
      j=intin(int4  )
      idx=kdeks(i)+j-ifm
      if(idx.lt.1.or.idx.gt.imax) go to 3
      k=intin(int4+3)
      iln=(idx-1)/im+1
      kz=idx-(iln-1)*im
      mq=mlr(iln)+1
      mm=mli(iln)+1
      t(mq)=q(jj)
      l=intin(int4+2)
      if(oprint(31)) then
       write(6,5566) jj,i,j,k,l,q(jj)
5566   format(1x,'ad2ijk: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      is(mm)=(kz-1)*inb+jdeks(k)+l
      if=nfg(iln)+1
      if(if.eq.ib) go to 10
      nfg(iln)=if
      mlr(iln)=mq
      mli(iln)=mm
      go to 3
 10   call stopbk
      nwbnwb=ib
      lnklnk=nrc(iln)
      md=mq-ib
      me=mm-ib
c     do 5566 loop=1,ib
c5566 gout(loop)=t(md+loop)
      call dcopy(ib,t(md+1),1,gout,1)
      call pack(gout(nsz341),32,is(me+1),nsz340)
      call sttout
      mlr(iln)=md
      mli(iln)=me
      nrc(iln)=irec
      nfg(iln)=0
      irec=irec+nsz
   3  int4=int4+4
      go to 11
 705  do 13 jj=1,il
      nz=nfg(jj)
      if(nz.eq.0)go to 13
      call stopbk
      nwbnwb=nz
      lnklnk=nrc(jj)
      mq=mlr(jj)-nz
      mm=mli(jj)-nz
c     do 5588 loop=1,nz
c5588 gout(loop)=t(mq+loop)
      call dcopy(nz,t(mq+1),1,gout,1)
      call pack(gout(nsz341),32,is(mm+1),nsz340)
      call sttout
      mlr(jj)=mq
      mli(jj)=mm
      nrc(jj)=irec
      irec=irec+nsz
 13   continue
c
      call stopbk
c
      imc=im
      idx=kxl+isb
      jdx=kyl+jsb
         i=isb
 200      i=i+1
        j=jsb
c     oijkl(1)=i
      idx=idx+1
      it=nir(idx)
      i5=iky(it)
      kdx=jdx
 201     j=j+1
      if(imc.lt.im)  goto 40
      if(ilc.eq.il) go to 50
      ilc=ilc+1
      imc=0
      iwa=0
c     do 31 kk=1,icak
c 31  t(kk)=0.0
      call vclr(t,1,icak)
      lnklnk=nrc(ilc)
      go to 32
 33   iblok=lnklnk
      call rdbak(iblok)
      call stopbk
      call unpack(gin(nsz341),32,ixa,nsz340)
      call dsctr(nwbnwb,gin,ixa,t)
 32   if(lnklnk.ne.maxb)go to 33
  40  lv=nzl
c     oijkl(2)=j
      kdx=kdx+1
      kt=nir(kdx)
      if(kt.gt.it) go to 940
      kt=jsym(i5+kt)
      go to 941
  940 kt=iky(kt)+it
      kt=jsym(kt)
  941 k5=iky(kt)
      imc=imc+1
      do 80 km=1,njx
      lv=lv+1
      kd=npal(lv)
      kr=irper(lv)
      if(kr.gt.kt) go to 81
      kr=kr+k5
      go to 82
  81  kr=iky(kr)+kt
  82  kr=jsym(kr)
      inr=inper(kr)
      if(inr.eq.0) go to 80
      lp=npal(inr)
      kw=lp+ij(inr)
      mmx=nbal(inr)
      mmm=ntil(inr)
      lx=nbal(lv)
      kk=ntil(lv)
      lw=kd+ij(lv)
      do 41 k=kd,lw
      kk=kk+1
      mcp=ncomp(kk)
      do 39 ll=1,mjb
  39  v(ll)=0.0d0
c     oijkl(3)=k
      do 43 ll=1,mcp
      lx=lx+1
      ld=lsym(lx)
      if(ld.lt.0) go to 44
      iw=iwa+jdeks(ld)
      do 70 l=1,mjb
      iw=iw+1
      v(l)=v(l)+t(iw)
  70  continue
      go to 43
  44  iw=iwa+jdeks(-ld)
      do 71 l=1,mjb
      iw=iw+1
      v(l)=v(l)-t(iw)
  71  continue
  43  continue
      mx=mmx
      mm=mmm
      do 55 l=lp,kw
      mm=mm+1
      ncp=ncomp(mm)
      sum=0.0d0
      do 53 ll=1,ncp
      mx=mx+1
      ld=lsym(mx)
      if(ld.lt.0) go to 54
      sum=sum+v(ld)
      go to 53
  54  sum=sum-v(-ld)
  53  continue
       if(dabs(sum).lt.1.0d-12) go to 55
c     oijkl(4)=l
      ig=ig+1
      pq(ig)=sum
      ig4=ig4+2
      intout(ig4  )=i4096(i)+j
      intout(ig4+1)=i4096(k)+l
      if(ig.lt.ifrk) go to 55
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2ijk labels'
       write(6,92233) (pli(loop),loop=1,250)
92233  format(5(1x,z16))
      endif
c
      write(mtape) p
      ig=0
      ig4=-1
  55  continue
  41  continue
  80  continue
      iwa=iwa+inb
         if(j.lt.mjy) go to 201
       jsb=0
        jdx=kyl
      if(i.lt.mjx) go to 200
      ig4=ig4+2
      intout(ig4)=0
      intout(ig4+1)=0
      call pack(pli,16,intout,1000)
c
      if(oprint(31)) then
       write(6,*)' ad2ijk labels'
       write(6,92233) (pli(loop),loop=1,250)
      endif
c
      write(mtape) p
      return
  50  isb=i-1
      jsb=j-1
      call rewftn(ntape)
      do 58 i=1,il
      nfg(i)=0
  58  nrc(i)=maxb
   2  continue
      return
      end
      subroutine amrd00
      implicit real*8  (a-h,o-z), integer (i-n)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      common/linkmr/new(255)
      common/b/sum(2),itape(4),ltape,mscr(4),ifrk
      common/blkin/gin(510),mword
      common/aplus/g(500),gm(250),n2
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      dimension f(751)
      common/craypk/intin(1360),intout(1000)
      equivalence (f(1),g(1))
      external fget
c
      jx=1
      write(ltape)jx
      write(ltape)jx,jx,jx,jx
      icl=0
      ic4=-1
      call setsto(1360,0,intin)
      do 4 ifile=1,lfile
      lbl=llblk(ifile)
      if(lbl)120,4,120
 120   num=lotape(ifile)
      call search(liblk(ifile),num)
 710  call fget(gin,mmmm,num)
      if(oprint(31)) then
        write(6,*)'amrd00: mmmm,mword = ', mmmm,mword
      endif
      if(mmmm)10,4,10
10    call unpack(gin(num2e+1),lab816,intin,numlab)
      int4=1
      do 3 iword=1,mword
      j=new(intin(int4  ))
      i=new(intin(int4+1))
      l=new(intin(int4+2))
      k=new(intin(int4+3))
      ic4=ic4+2
      icl=icl+1
      g(icl)=gin(iword)
      if(oprint(31)) then
       write(6,5566) iword,i,j,k,l,g(icl)
5566   format(1x,'amrd00: i,j,k,l,val = ',i4,2x,4i4,5x,f20.10)
      endif
      ij=min(i,j)+i4096(max(i,j))
      kl=min(k,l)+i4096(max(k,l))
      intout(ic4  )=max(ij,kl)
      intout(ic4+1)=min(ij,kl)
      if(icl.lt.ifrk)go to 3
      n2=icl
      call pack(gm,16,intout,1000)
      write(ltape)f
      icl=0
      ic4=-1
 3    int4=int4+4
      lbl=lbl+1
      if(lbl)710,4,710
 4    continue
      ic4=ic4+2
      intout(ic4  )=0
      intout(ic4+1)=0
      call pack(gm,16,intout,1000)
      n2=0
      write(ltape)f
      return
      end
      subroutine ver_mrdci1(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mrdci1.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
