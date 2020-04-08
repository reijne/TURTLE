c 
c  $Author: jvl $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util2.m,v $
c  $State: Exp $
c  
c     deck=util2
c ******************************************************
c ******************************************************
c             =   util2  =
c ******************************************************
c ******************************************************
      subroutine blocki
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/craypk/integ(680),ipad(680)
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
      common/blkin/gout(510),mword
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      if(omaxb) return 
      mword=icount-1
      if(mword.eq.0)go to 41
      ochek = .false.
c
c...    for incore scf go on to see how much space is needed ...
      if(imaxb_ic.eq.1) go to 50
c
      if(nopk.ne.1) then
      if(.not.opank)then
         call pack(gout(num2ep+1),lab1632,integ,numlabp)
      else
         call pack(gout(num2ejk+num2ejk+1),lab1632,
     +             integ,numlabjk)
      endif
c
      else
c
      labss = lab816 + lab816
      numlabi = numlab / 2
      call pack(gout(num2e+1),labss,integ,numlabi)
      endif
      if(iposun(mainp).ne.iblkmp)callsearch(iblkmp,mainp)
      call put(gout(1),m511,mainp)
c ======================================================
c = debug print for integrals (nprint=4)
c ======================================================
      if(.not.out) go to 240
      if(nopk.eq.1) go to 230
      if(opank) then
      call pr2ejk(mainp,iblkmp)
      else
      call pr2ep (mainp,iblkmp)
      endif
      go to 240
 230  call pr2ei (mainp,iblkmp)
 240  continue
c ============================================================
 50   icount=1
      ic4 = 1
      nrec=nrec+1
      iblkmp=iblkmp+1
      mblp=mblp+1
      if(mblp)41,40,41
 40    m2last(m2file)=iblkmp
       mfilep=mfilep+1
      if(mfilep.le.n2file)go to 1
      if (omem(mainp)) then
         imaxb_ic = 1
c...    memory segment is the only one ??
         mfilep = mfilep - 1
      else
         omaxb = .true.
      end if
c
      if (omaxb) then
      write(iwr,2)yed(mainp)
 2    format(//
     *' *** sorry ***'//
     *' size of mainfile section ',a4,' has been exceeded'//
     *' restart with additional section allocated'//)
      end if
c
      go to 41
  1   mainp=n2tape(mfilep)
      iblkmp=n2blk(mfilep)
      mblp=iblkmp-n2last(mfilep)
      m2file=m2file+1
      m2tape(m2file)=mainp
      m2blk(m2file)=iblkmp
      m2last(m2file)=-1
  41  return
      end
      subroutine qout(g)
c
c ..  qout is now restricted to conventional scf
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
c
      common/craypk/integ(680),ipad(680)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      common/blkin/gout(510),nword
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      dimension g(*)
c
c
c     ----- pack the 4 indices of two integrals into one word
c     ----- write label + integral on mainfile
c
      ijn = 0
      jmax = maxj
      do 2600 i = mini,maxi
      if (oianj) jmax = i
      i1 = loci + i
      ipack = i4096(i1)
      do 2400 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      i2 = locj + j
      if(i1-i2)10, 20, 20
 10   lab1 = i4096(i2) + i1
      go to 30
 20   lab1 = ipack + i2
 30   lmax = maxl
      kln = 0
      do 2200 k = mink,maxk
      if (okanl) lmax = k
      i3 = lock + k
      kpack = i4096(i3)
      do 2000 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 2400
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 2000
      i4 = locl + l
      if(i3-i4) 40, 50, 50
 40   lab2 = i4096(i4) + i3
      go to 60
 50   lab2 = kpack + i4
 60   gout(icount) = val
      integ(ic4  )=max(lab1,lab2)
      integ(ic4+1)=min(lab1,lab2)
      ic4=ic4+2
      icount = icount+1
      if (icount .le. nintmx) go to 2000
      call blocki
      if(omaxb)go to 261
 2000 continue
 2200 continue
 2400 continue
 2600 continue
 261  return
      end
      subroutine dbuild(fock,dmat,g,
     &  fac1, fac2, facex, ocoul, oexch)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c

      dimension g(*),fock(*),dmat(*)
c
c     ----- DSCF fock builder
c
      ijn = 0
      jmax = maxj

      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
c
c Coulomb term
c
        if(ocoul)then
           val2=val+val
           val4=val2+val2
           fock(itr12) = facij*val4*dmat(itr34) + fock(itr12)
           if(itr12 .ne. itr34)then
              fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
           endif
        endif
c
c exchange
c
        if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
         val = val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2

         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
  200 continue
  220 continue
  240 continue
  260 continue
      end

      subroutine dir_build_uhf(fock,ak,p,q,g,
     &   fac1, fac2, facex, ocoul, oexch)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      dimension g(*),fock(*),ak(*),p(*),q(*)
c
c     ----- Direct open-shell UHF Fock builder 
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
c
c Coulomb term
c
        if(ocoul)then

           gik=val
           val2=gik+gik
           val4=val2+val2

           fock(itr12) = facij*val4*p(itr34) + fock(itr12)
           if(itr12 .ne. itr34)then
              fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
           endif

        endif
c
c exchange
c
        if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
           gik=val*facex
           val2=gik+gik
           val4=val2+val2
           gil=gik
           if(i1.eq.i3.or.i2.eq.i4)gik=val2
           if(i2.eq.i3)gil=val2
           if(i2.ge.i3)goto 1
           itr23=ikyk+i2
           if(i2.ge.i4)goto 1
           itr24=iky(i4)+i2
 1         ajk=fock(itr23)-gil*p(itr14)
           bjk=ak(itr23)+gil*q(itr14)
           ail=fock(itr14)-gil*p(itr23)
           bil=ak(itr14)+gil*q(itr23)
           aik=fock(itr13)-gik*p(itr24)
           bik=ak(itr13)+gik*q(itr24)
           fock(itr24)=fock(itr24)-gik*p(itr13)
           ak(itr24)=ak(itr24)+gik*q(itr13)
           fock(itr23)=ajk
           ak(itr23)=bjk
           fock(itr14)=ail
           ak(itr14)=bil
           fock(itr13)=aik
           ak(itr13)=bik

      endif
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open(coul,exch,dens,g)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      dimension g(*),coul(*),exch(*),dens(*)
c
c     ----- Direct open-shell SCF J & K builder (nshell=1)
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        gik=val
        val2=val+val
        ikyi=iky(i1)
        ikyj=iky(i2)
        ikyk=iky(i3)
        itr13=ikyi+i3
        itr14=ikyi+i4
        itr12=ikyi+i2
        itr23=ikyj+i3
        itr24=ikyj+i4
        itr34=ikyk+i4
        gil=val
        if(i1.eq.i3.or.i2.eq.i4)gik=val2
        if(i2.eq.i3)gil=val2
        if(i2.ge.i3)goto 280
        itr23=ikyk+i2
        if(i2.ge.i4)goto 280
        itr24=iky(i4)+i2
  280   bij=val2*dens(itr34)+coul(itr12)
        coul(itr34)=val2*dens(itr12)+coul(itr34)
        coul(itr12)=bij
        bjk=exch(itr23)+gil*dens(itr14)
        bil=exch(itr14)+gil*dens(itr23)
        bik=exch(itr13)+gik*dens(itr24)
        exch(itr24)=exch(itr24)+gik*dens(itr13)
        exch(itr23)=bjk
        exch(itr14)=bil
        exch(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open2(l2,coul,exch,dens,g)
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
      logical ospbas, onocnt, ochek, opdbas, opfbas, opgbas
      integer kad, istd
      common /ijlab/ kad(mxshel),ospbas,onocnt,istd,ochek,
     +               opdbas,opfbas,opgbas
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      dimension g(*),coul(*),exch(*),dens(*)
c
c     ----- Direct open-shell SCF J & K builder (nshell > 1)
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=val+val
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      gil=val
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)go to 280
      itr23=ikyk+i2
      if(i2.ge.i4)go to 280
      itr24=iky(i4)+i2
  280 continue
      do 300 iiii=1,nsheld
      bij=val2*dens(itr34)+coul(itr12)
      coul(itr34)=val2*dens(itr12)+coul(itr34)
      coul(itr12)=bij
      bjk=exch(itr23)+gil*dens(itr14)
      bil=exch(itr14)+gil*dens(itr23)
      bik=exch(itr13)+gik*dens(itr24)
      exch(itr24)=exch(itr24)+gik*dens(itr13)
      exch(itr23)=bjk
      exch(itr14)=bil
      exch(itr13)=bik
      itr12=itr12+l2
      itr34=itr34+l2
      itr13=itr13+l2
      itr14=itr14+l2
      itr23=itr23+l2
      itr24=itr24+l2
  300 continue
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine pkfile(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,
     +                   g,q)
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
c ***
      dimension q(*),g(*)
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c ***
      common/craypk/integ(680),ipad(680)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      common/blkin/gout(510),nword
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer igt, jgt, kgt, lgt
      integer ibt, jbt, kbt, lbt, lenddd
      common /flips/ igt(3),jgt(3),kgt(3),lgt(3),
     +               ibt(mxprms),jbt(mxprms),kbt(mxprms),lbt(mxprms),
     +               lenddd
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
      data dzero,pt5 /0.0d0,0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
c
c     type = 1 for (ii ii ii ii)
c            2     (ii jj jj jj)
c            3     (ii ii kk kk)
c            4     (ii ii ii ll)
c            5     (ii jj kk kk)
c            6     (ii jj jj ll)
c            7     (ii ii kk ll)
c            8     (ii jj kk ll)
c
      if (ii-jj) 100, 100, 240
  100 if (jj-kk) 120, 120, 180
  120 if (kk-ll) 140, 140, 160
  140 ntyp = 1
      go to 380
  160 ntyp = 4
      go to 380
  180 if (kk-ll) 200, 200, 220
  200 ntyp = 3
      go to 380
  220 ntyp = 7
      go to 380
  240 if (jj-kk) 260, 260,320
  260 if (kk-ll) 280,280,300
  280 ntyp = 2
      go to 380
  300 ntyp = 6
      go to 380
  320 if (kk-ll) 340,340,360
  340 ntyp = 5
      go to 380
  360 ntyp = 8
  380 continue
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 50626
      norg3 = 101251
      if (oskpa .and. .not. onpsym) norg2 = 1
      if (oskpc .and. .not. onpsym) norg3 = 50626
      if (oskpb .and. .not. onpsym) norg3 = 1
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      do 1060 i = mini,maxi
      if (oianj) jmax = i
      do 1040 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 1020 k = mink,kmaxx
      if (okanl) lmax = k
      do 1000 l = minl,lmax
      ia = i-mini+1
      ja = j-minj+1
      ka = k-mink+1
      la = l-minl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1 = ibt(ia)+jbt(ja)+kbt(ka)+lbt(la)+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2 = ibt(ia)+jbt(ka)+kbt(ja)+lbt(la)+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      go to (400,420,400,440,420,420,440,420),ntyp
  400 if (ia .eq. ja) go to 440
  420 n3 = ibt(ia)+jbt(la)+kbt(ja)+lbt(ka)+norg3
      go to 460
  440 n3 = ibt(ja)+jbt(ka)+kbt(ia)+lbt(la)+norg3
  460 continue
c
c     ----- form first linear combination -----
c
      jump = 1
      i1 = loci+i
      i2 = locj+j
      i3 = lock+k
      i4 = locl+l
      if (i2 .eq. i3) jump = 2
      if ((i2 .eq. i4) .or. (i1 .eq. i3)) jump = 3
      g1 = g(n1)
      g2 = g(n2)
      g3 = g(n3)
      go to (480,500,520,540),ind
  480 valk = g2+g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  500 valk = g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  520 valk = g2
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  540 valk = dzero
      valp = (g1+g1)+(g1+g1)
  560 continue
      g(n1) = valp
c     nn = n1
      go to 820
c
c     ----- form second linear combination -----
c
  580 go to (600,620,640,660),ind
  600 valk = g3+g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  620 valk = g1+g3
      valp = -valk
      go to 680
  640 valk = g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  660 valk = g1
      valp = -valk
  680 continue
      g(n2) = valp
c     nn = n2
      jump = 2
      if ((i1 .eq. i2) .or. (i3 .eq. i4)) jump = 3
      n = i2
      i2 = i3
      i3 = n
      go to 820
c
c     ----- form third linear combination -----
c
  700 go to (720,740,760,780),ind
  720 valk = g1+g2
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  740 valk = g1
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  760 valk = g1+g2
      valp = -valk
      go to 800
  780 valk = g1
      valp = -valk
  800 continue
      g(n3) = valp
c     nn = n3
      i2 = locl+l
      i3 = locj+j
      i4 = lock+k
      jump = 3
  820 continue
c
c     ----- store integral and indices -----
c
      if (opank) go to 880
c
c     ----- -p- supermatrix only. -----
c
      if ( dabs(valp) .lt. cutoff) go to 980
c     valp0 = valp
      if (i1 .eq. i3 .and. i2 .eq. i4) valp = valp*pt5
      integ(ic4  )=iky(i1)+i2
      integ(ic4+1)=iky(i3)+i4
      ic4=ic4+2
      gout(icount) = valp
      icount = icount+1
      if (icount .le. nintmx) go to 980
c ***
      if (osortp) then
        call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
      else
        call blocki
      endif
c ***
      if(omaxb)go to 1061
      go to 980
  880 continue
c
c     ----- -p- and -k- supermatrices -----
c
      if (dabs(valp).lt.cutoff .and. dabs(valk).lt.cutoff) go to 980
c     valp0 = valp
c     valk0 = valk
      if (i1 .ne. i3 .or. i2 .ne. i4) go to 900
      valp = valp*pt5
      valk = valk*pt5
  900 continue
      integ(ic4  )=iky(i1)+i2
      integ(ic4+1)=iky(i3)+i4
      ic4 = ic4 + 2
      gout(icount) = valp
      gout(num2ejk+icount) = valk
      icount = icount+1
      if (icount .le. nintmx) go to 980
      call blocki
      if(omaxb)go to 1061
  980 go to (580,700,1000),jump
 1000 continue
 1020 continue
 1040 continue
 1060 continue
 1061 return
      end
      subroutine pfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,
     +                  g,q)
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
c ***
      dimension g(*),q(*)
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c ***
      common/craypk/integ(680),ipad(680)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      dimension jump(256),i1(256),i2(256),i3(256),i4(256)
     *,n1(256),n2(256),n3(256)
     *,o12(256),o34(256)
     *,valp1(256),i2out(256),i3out(256)
     *,valp2(256),valp3(256)
      common/blkin/goutx(510),nword
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 257
      norg3 = 513
      if ( .not. oskpa) go to 100
      if ( .not. onpsym) norg2 = 1
      ib2 = ib1
      jb2 = jb1
      kb2 = kb1
      lb2 = lb1
  100 if ( .not. oskpb) go to 140
      if ( .not. onpsym) norg3 = norg1
      if (ii .ne. kk) go to 120
      if (jj .eq. ll) go to 120
      ib3 = kb1
      jb3 = lb1
      kb3 = jb1
      lb3 = ib1
      go to 180
  120 ib3 = ib1
      jb3 = jb1
      kb3 = lb1
      lb3 = kb1
      go to 180
  140 continue
      if ( .not. oskpc) go to 180
      if ( .not. onpsym) norg3 = norg2
      if (ii .ne. jj) go to 160
      if (kk .eq. ll) go to 160
      ib3 = kb2
      jb3 = lb2
      kb3 = ib2
      lb3 = jb2
      go to 180
  160 ib3 = ib2
      jb3 = jb2
      kb3 = kb2
      lb3 = lb2
  180 continue
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      mmijkl=0
      do 780 i = mini,maxi
      if (oianj) jmax = i
      do 780 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 780 k = mink,kmaxx
      if (okanl) lmax = k
      do 780 l = minl,lmax
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(mmijkl) = ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2(mmijkl) = ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n3(mmijkl) = ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
      do 781 i=1,mmijkl
      jump(i) = 1
      if (i2(i) .eq. i3(i)) jump(i) = 2
      if (i2(i) .eq. i4(i) .or. i1(i) .eq. i3(i)) jump(i) = 3
      o34(i)=i3(i).eq.i4(i)
      o12(i)=i1(i).eq.i2(i)
 781   continue
c
      loop2=0
      loop3=0
c
      go to (200,220,240,260),ind
c
200   do 280 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      g3=g(n3(i))
      valp=4.0d0*g1-g2-g3
      valp1(i)=valp
      if(jump(i)-2)2001,2002,280
 2001 valp=4.0d0*g2-g1-g3
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 280
 2002 valp=4.0d0*g3-g1-g2
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 280  continue
      go to 3000
c
220   do 380 i=1,mmijkl
      g1=g(n1(i))
      g3=g(n3(i))
      valp=4.0d0*g1-g3
      valp1(i)=valp
      if(jump(i)-2)3001,3002,380
 3001 valp=-g1-g3
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 380
 3002 valp=4.0d0*g3-g1
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 380  continue
      go to 3000
c
240   do 480 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      valp=4.0d0*g1-g2
      valp1(i)=valp
      if(jump(i)-2)4001,4002,480
 4001 valp=4.0d0*g2-g1
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 480
 4002 valp=-g1-g2
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 480  continue
      go to 3000
c
260   do 580 i=1,mmijkl
      g1=g(n1(i))
      valp=4.0d0*g1
      valp1(i)=valp
      if(jump(i)-2)5001,5002,580
 5001 valp=-g1
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 580
 5002 valp=-g1
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 580  continue
c
c     ----- store integral and indices -----
c
c     ----- -p- supermatrix only. -----
c
3000  do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff) go to 560
      if (i1(i) .eq. i3(i) .and. i2(i) .eq. i4(i))
     * valp1(i)=valp1(i)*pt5
      integ(ic4  )=iky(i1(i))+i2(i)
      integ(ic4+1)=iky(i3(i))+i4(i)
      ic4 = ic4 + 2
      goutx(icount) = valp1(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 560   continue
c
      do 570 i=1,loop2
      if ( dabs(valp2(i)) .lt. cutoff) go to 570
      iii=i2out(i)
      if (i1(iii) .eq. i2(iii) .and. i3(iii) .eq. i4(iii))
     * valp2(i)=valp2(i)*pt5
      integ(ic4  )=iky(i1(iii))+i3(iii)
      integ(ic4+1)=iky(i2(iii))+i4(iii)
      ic4 = ic4 + 2
      goutx(icount) = valp2(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 570   continue
c
      do 590 i=1,loop3
      if ( dabs(valp3(i)) .lt. cutoff) go to 590
      iii=i3out(i)
      if (i1(iii) .eq. i2(iii) .and. i4(iii) .eq. i3(iii))
     * valp3(i)=valp3(i)*pt5
      integ(ic4  )=iky(i1(iii))+i4(iii)
      integ(ic4+1)=iky(i2(iii))+i3(iii)
      ic4 = ic4 + 2
      goutx(icount) = valp3(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 590   continue
  782 return
      end
      subroutine pfi70s(ii,jj,kk,ll,g,q)
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
      dimension g(*),q(*)
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c ***
      common/craypk/integ(680),ipad(680)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
      dimension i1(256),i2(256),i3(256),i4(256),
     *            n1(3,256),valpg(3,256),
     *            valp1(256),valp2(256),valp3(256)
      common/blkin/goutx(510),nword
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c     data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
c     ind = 1
      norg1 = 1
      norg2 = 257
      norg3 = 513
      mmijkl=0
      do 780 i = mini,maxi
      do 780 j = minj,maxj
      do 780 k = mink,maxk
      do 780 l = minl,maxl
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(1,mmijkl)=ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n1(2,mmijkl)=ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n1(3,mmijkl)=ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
      call dgthr(mmijkl*3,g,valpg,n1)
c
      do 280 i=1,mmijkl
      valp1(i)=4.0d0*valpg(1,i)-valpg(2,i)-valpg(3,i)
      valp2(i)=4.0d0*valpg(2,i)-valpg(1,i)-valpg(3,i)
      valp3(i)=4.0d0*valpg(3,i)-valpg(1,i)-valpg(2,i)
 280  continue
c
c     ----- store integral and indices -----
c
c     ----- -p- supermatrix only. -----
c
      do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff) go to 560
      integ(ic4  )=iky(i1(i))+i2(i)
      integ(ic4+1)=iky(i3(i))+i4(i)
      ic4 = ic4 + 2
      goutx(icount) = valp1(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 560   continue
c
      do 570 i=1,mmijkl
      if ( dabs(valp2(i)) .lt. cutoff) go to 570
      integ(ic4  )=iky(i1(i))+i3(i)
      integ(ic4+1)=iky(i2(i))+i4(i)
      ic4 = ic4 + 2
      goutx(icount) = valp2(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 570   continue
c
      do 590 i=1,mmijkl
      if ( dabs(valp3(i)) .lt. cutoff) go to 590
      integ(ic4  )=iky(i1(i))+i4(i)
      integ(ic4+1)=iky(i2(i))+i3(i)
      ic4 = ic4 + 2
      goutx(icount) = valp3(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 590   continue
  782 return
      end
      subroutine pkfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,g)
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
      dimension g(*)
      common/craypk/integ(680),ipad(680)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
      dimension jump(256),i1(256),i2(256),i3(256),i4(256)
     *, n1(256),n2(256),n3(256)
     *, o12(256),o34(256)
     *, valp1(256),i2out(256),i3out(256)
     *, valp2(256),valp3(256),valk1(256),valk2(256),valk3(256)
      common/blkin/gout(510),nword
c     logical opk,opank,oskips
c     common/pkfil/opk,opank,oskips
      logical opk,opank
      common/pkfil/opk,opank
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 257
      norg3 = 513
      if ( .not. oskpa) go to 100
      if ( .not. onpsym) norg2 = 1
      ib2 = ib1
      jb2 = jb1
      kb2 = kb1
      lb2 = lb1
  100 if ( .not. oskpb) go to 140
      if ( .not. onpsym) norg3 = norg1
      if (ii .ne. kk) go to 120
      if (jj .eq. ll) go to 120
      ib3 = kb1
      jb3 = lb1
      kb3 = jb1
      lb3 = ib1
      go to 180
  120 ib3 = ib1
      jb3 = jb1
      kb3 = lb1
      lb3 = kb1
      go to 180
  140 continue
      if ( .not. oskpc) go to 180
      if ( .not. onpsym) norg3 = norg2
      if (ii .ne. jj) go to 160
      if (kk .eq. ll) go to 160
      ib3 = kb2
      jb3 = lb2
      kb3 = ib2
      lb3 = jb2
      go to 180
  160 ib3 = ib2
      jb3 = jb2
      kb3 = kb2
      lb3 = lb2
  180 continue
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      mmijkl=0
      do 780 i = mini,maxi
      if (oianj) jmax = i
      do 780 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 780 k = mink,kmaxx
      if (okanl) lmax = k
      do 780 l = minl,lmax
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(mmijkl) = ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2(mmijkl) = ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n3(mmijkl) = ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
      do 781 i=1,mmijkl
      jump(i) = 1
      if (i2(i) .eq. i3(i)) jump(i) = 2
      if (i2(i) .eq. i4(i) .or. i1(i) .eq. i3(i)) jump(i) = 3
      o34(i)=i3(i).eq.i4(i)
      o12(i)=i1(i).eq.i2(i)
 781   continue
c
      loop2=0
      loop3=0
c
      go to (200,220,240,260),ind
c
200   do 280 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      g3=g(n3(i))
      valk=g2+g3
      valp=4.0d0*g1-valk
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)2001,2002,280
 2001 valk=g1+g3
      valp=4.0d0*g2-valk
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 280
 2002 valk=g1+g2
      valp=4.0d0*g3-valk
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 280  continue
      go to 3000
c
220   do 380 i=1,mmijkl
      g1=g(n1(i))
      g3=g(n3(i))
      valk=g3
      valp=4.0d0*g1-valk
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)3001,3002,380
 3001 valk=g1+g3
      valp=-valk
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 380
 3002 valk=g1
      valp=4.0d0*g3-valk
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 380  continue
      go to 3000
c
240   do 480 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      valk=g2
      valp=4.0d0*g1-valk
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)4001,4002,480
 4001 valk=g1
      valp=4.0d0*g2-valk
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 480
 4002 valk=g1+g2
      valp=-valk
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 480  continue
      go to 3000
c
260   do 580 i=1,mmijkl
      g1=g(n1(i))
      valp=4.0d0*g1
      valp1(i)=valp
      valk1(i)=0.0d0
      if(jump(i)-2)5001,5002,580
 5001 valk=g1
      valp=-valk
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 580
 5002 valk=g1
      valp=-valk
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 580  continue
c
c     ----- store integral and indices -----
c
c     ----- -p+k- supermatrix only. -----
c
3000  do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff.and.  dabs(valk1(i)).lt.cutoff)
     * go to 560
      if (i1(i) .ne. i3(i) .or. i2(i) .ne. i4(i))
     * go to 5601
      valp1(i)=valp1(i)*pt5
      valk1(i)=valk1(i)*pt5
5601  continue
      integ(ic4  )=iky(i1(i))+i2(i)
      integ(ic4+1)=iky(i3(i))+i4(i)
      ic4 = ic4 + 2
      gout(icount) = valp1(i)
      gout(num2ejk+icount) = valk1(i)
      icount = icount+1
      if (icount .le. nintmx) go to 560
c
      call blocki
      if(omaxb)go to 782
 560   continue
c
      do 570 i=1,loop2
      if ( dabs(valp2(i)) .lt. cutoff. and.  dabs(valk2(i)).lt.cutoff)
     *go to 570
      iii=i2out(i)
      if (i1(iii) .ne. i2(iii) .or. i3(iii) .ne. i4(iii)) go to 5701
       valp2(i)=valp2(i)*pt5
       valk2(i)=valk2(i)*pt5
 5701 continue
      integ(ic4  )=iky(i1(iii))+i3(iii)
      integ(ic4+1)=iky(i2(iii))+i4(iii)
      ic4 = ic4 + 2
      gout(icount) = valp2(i)
      gout(num2ejk+icount) = valk2(i)
      icount = icount+1
      if (icount .le. nintmx) go to 570
c
      call blocki
      if(omaxb)go to 782
 570   continue
c
      do 590 i=1,loop3
      if ( dabs(valp3(i)) .lt. cutoff. and . dabs(valk3(i)).lt.cutoff)
     * go to 590
      iii=i3out(i)
      if (i1(iii) .ne. i2(iii) .or. i4(iii) .ne. i3(iii))
     * go to 5901
       valp3(i)=valp3(i)*pt5
       valk3(i)=valk3(i)*pt5
 5901 continue
      integ(ic4  )=iky(i1(iii))+i4(iii)
      integ(ic4+1)=iky(i2(iii))+i3(iii)
      ic4 = ic4 + 2
      gout(icount) = valp3(i)
      gout(num2ejk+icount) = valk3(i)
      icount = icount+1
      if (icount .le. nintmx) go to 590
c
      call blocki
      if(omaxb)go to 782
 590   continue
  782 return
      end
      subroutine zsortp(gbuf,klbuf,iptbuf,ijbas)
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
c ***
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
      dimension gbuf(lengbf,*),klbuf(lengbf,*),iptbuf(*)
      dimension gx(340),ijx(340),klx(340)
c ***
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
      common/craypk/intij(340),intkl(340),ipad(680)
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
      common/blkin/goutx(510),nword
c ***
c ***
c *** copy goutx, intij, intkl into local work space,
c *** sort contents into buffers writing out as usual
c ***
      mword = icount - 1
      icount = 1
      ic4 = 1
      if(mword.le.0) return
      call dcopy(mword,goutx,1,gx,1)
      call icopy(mword,intij(1),2,ijx,1)
      call icopy(mword,intij(2),2,klx,1)
      do 10 i=1,mword
          ibuf = ijx(i) - ijbas
          if(ibuf.le.0 .or. ibuf.gt.ngbf)
     &      call caserr('ibuf value is invalid.')
          iptbuf(ibuf) = iptbuf(ibuf) + 1
          gbuf (iptbuf(ibuf),ibuf) = gx(i)
          klbuf(iptbuf(ibuf),ibuf) = klx(i)
          if(iptbuf(ibuf).eq.lengbf) call clrgbf(ibuf,gbuf,klbuf,iptbuf)
10    continue
      return
      end
      subroutine clrgbf(ibuf,gbuf,klbuf,iptbuf)
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
c ***
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
      dimension gbuf(lengbf,*),klbuf(lengbf,*),iptbuf(*)
c ***
      common/craypk/intij(680),ipad(680)
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
      common/blkin/goutx(340)
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c ***
c *** empty out buffer ibuf with usual blocki call sequence
c ***
c *** will not work for lengbf other than 340 as /blkin/ and
c *** /crypak/ will be destroyed by zsortp.
c ***
      mword = iptbuf(ibuf)
      if (mword.le.0) return
      ipt = 1
      ijval = ijbase + ibuf
10    nl = nintmx - icount + 1
      nleft = min(nl,mword)
c *** put smallest value at front to minimise precision loss
      i = idamin(nleft,gbuf(ipt,ibuf),1)
      i=i+ipt-1
      gx = gbuf(ipt,ibuf)
      gbuf(ipt,ibuf) = gbuf(i,ibuf)
      gbuf(i,ibuf) = gx
      m = klbuf(ipt,ibuf)
      klbuf(ipt,ibuf) = klbuf(i,ibuf)
      klbuf(i,ibuf) = m
      iqq = 2*icount-1
      do 11 i=1,nleft
          intij(iqq)=ijval
          intij(iqq+1) = klbuf(i+(ipt-1),ibuf)
11        iqq = iqq + 2
      call dcopy(nleft,gbuf(ipt,ibuf),1,goutx(icount),1)
      call rinsrt(goutx(icount),nleft)
      icount = icount + nleft
      ic4 = 2*icount-1
      if (icount .gt. nintmx) call blocki
      mword = mword - nleft
      ipt = ipt + nleft
      if(mword .gt. 0) goto 10
      iptbuf(ibuf)=0
      return
      end
      subroutine rinsrt(r,m)
      implicit real*8  (a-h,o-z)
c ***
c *** insert last len=9 bits of m into r at position ipos=55.
c ***
      integer r,m,mtemp 
      mtemp = iand(511,m) 
      r=iand(-512,r) 
      r=ior(r,mtemp) 
      return
      end
      subroutine dbutci(ista,jsta,ksta,lsta)
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
      common/restar/nprint,itol,icut,normf,normp,nopk,
     *                     nrest,nrec,intloc,ist,jst,kst,lst
     *   ,nspfil(434), n5file,n5tape(20),n5blk(20),n5last(20)
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen
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
c
      real*8 de
      common/grad2/de(3,maxat)
c
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
c
      common/blkin/qqqq(510),mword
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
      write (iwr,9008)
c
c     ----- set starting parameters -----
c
      outd = nprint .eq. -4
      icutd=icut-1
      icutd=max(icutd,8)
      itold=itol
      itold=max(itold,16)
      dcut=done/(ten**icutd)
      tol1=e*(itold+2)
      tol2=e*itold
      tol3=done/ten**(itold+2)
      tol4=done/ten**itold
      onorm=normf.ne.1.or.normp.ne.1
c
c     ----- read in 1e-gradient -----
c
      call rdgrd(de,nrest,ista,jsta,ksta,lsta)
c
c     ----- position 2-particle density matrix file
c
      lfile=n5file
      do 1 i=1,lfile
      lotape(i)=n5tape(i)
      liblk  (i)=n5blk (i)
 1    llblk  (i)=liblk(i)-n5last(i)
c
      ifiled=1
      mword=0
      iwordd=0
      idmmc=lotape(1)
      iblkmm =liblk(1)
      lblmm  =llblk(1)
      call search(iblkmm,idmmc)
c
      return
 9008 format(/1x,22('-')/1x,'gradient of the energy'/1x,22('-')/)
      end
      subroutine dabci(ii,jj,kk,ll,q4,oeof,abdens)
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
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
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
      common/blkin/gin(510),mword
      common/craypk/ijkl(4,340)
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
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
      common/tgrad/dgout(9)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      dimension abdens(*)
c
      data eight,four,pt5/8.0d0,4.0d0,0.5d0/
      call vclr(abdens,1,lendd)
      if(oeof)return
      iword=iwordd
      oijdif=ii.ne.jj
      okldif=kk.ne.ll
      oikdif=ii.ne.kk.or.jj.ne.ll
      mini=kloc(ii)
      minj=kloc(jj)
      mink=kloc(kk)
      minl=kloc(ll)
      maxi=kloc(ii+1)-1
      maxj=kloc(jj+1)-1
      maxk=kloc(kk+1)-1
      maxl=kloc(ll+1)-1
1     if(iword.ne.mword)goto10
      if(lblmm.ne.0)go to 5
      ifiled=ifiled+1
      if(ifiled.gt.lfile)go to 888
      idmmc=lotape(ifiled)
      call search(liblk(ifiled),idmmc)
      lblmm=llblk(ifiled)
 5    callfind(idmmc)
      callget(gin,m)
      lblmm=lblmm+1
      if(m)2,888,2
2     iword=0
      call unpack(gin(num2e+1),lab816,ijkl,numlab)
10    iword=iword+1
      j = ijkl(1,iword)
      i = ijkl(2,iword)
      l = ijkl(3,iword)
      k = ijkl(4,iword)
      if(i.lt.mini)goto1
      if(i.gt.maxi)goto999
      if(j.lt.minj)goto1
      if(j.gt.maxj)goto999
      if(k.lt.mink)goto1
      if(k.gt.maxk)goto999
      if(l.lt.minl)goto1
      if(l.gt.maxl)goto999
      d4=eight
      if(i.eq.j)d4=four
      if(k.eq.l)d4=d4*pt5
      nn=(i-mini)*inc5+(j-minj)*inc4+(k-mink)*inc3+l-minl+inc2
      buff=gin(iword)*q4*d4
      abdens(nn)=buff
      if(oijdif)goto20
      nn=(j-minj)*inc5+(i-mini)*inc4+(k-mink)*inc3+l-minl+inc2
      abdens(nn)=buff
      if(okldif)goto20
      nn=(j-minj)*inc5+(i-mini)*inc4+(l-minl)*inc3+k-mink+inc2
      abdens(nn)=buff
20    if(okldif)goto30
      nn=(i-mini)*inc5+(j-minj)*inc4+(l-minl)*inc3+k-mink+inc2
      abdens(nn)=buff
30    if(oikdif)goto1
      nn=(k-mink)*inc5+(l-minl)*inc4+(i-mini)*inc3+j-minj+inc2
      abdens(nn)=buff
      if(oijdif)goto40
      nn=(k-mink)*inc5+(l-minl)*inc4+(j-minj)*inc3+i-mini+inc2
      abdens(nn)=buff
      if(okldif)goto40
      nn=(l-minl)*inc5+(k-mink)*inc4+(j-minj)*inc3+i-mini+inc2
      abdens(nn)=buff
40    if(okldif)goto1
      nn=(l-minl)*inc5+(k-mink)*inc4+(i-mini)*inc3+j-minj+inc2
      abdens(nn)=buff
      goto1
999   iwordd=iword-1
      return
888   oeof=.true.
      return
      end
      function chrint(i)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c converts an integer to a left justified character string
c without leading zeros, and with minus sign (if negative)
c   works for up to 15 decimal digit integers
      character*16 c,chrint
      character*1 digit(0:9)
      data c/' '/
      data digit/'0','1','2','3','4','5','6','7','8','9'/
      num=iabs(i)
      l=16
      do 1 j=16,2,-1
      new=num/10
      n=num-10*new
      if(n.gt.0) l=j
      c(j:j)=digit(n)
1     num=new
      if(i.lt.0) then
          l=l-1
          c(l:l)='-'
      endif
      chrint=c(l:16)
      return
      end
      function lstchr(a)
      character*(*) a
c routine to return index of last non-blank character in character var
      n=len(a)
      do 10 i=n,1,-1
      if( a(i:i).ne.' ' ) go to 20
10    continue
      i=0
20    lstchr=i
      return
      end
      subroutine wrtc(z,nword,iblk,num3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension z(*)
      common/junkcc/dchar(511)
      call search(iblk,num3)
      l=511
      j=1
      k=nword
 1    k=k-511
      if(k)10495,2,2
10495 l=k+511
c2    call cmoved(z(j),dchar,l)
 2    do 3 loop=1,l
 3    read(z(j+loop-1),'(a8)') dchar(loop)
      call put(dchar,l,num3)
      j=j+511
      if(k)10500,10500,1
10500 return
      end
      subroutine rdchr(z,nword,iblk,num3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junkcc/dchar(511)
      dimension z(*)
      call search(iblk,num3)
      j=1
10491 call find(num3)
      call get(dchar,l)
c     call dmovec(dchar,z(j),l)
      do 3 loop=1,l
 3    write(z(j+loop-1),'(a8)') dchar(loop)
      j=j+l
      if(j.le.nword)go to 10491
      if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present ', nword, j-1
       call caserr('invalid number of words in rdchr')
      endif
      return
      end
      subroutine wrtcs(z,nword,num3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension z(*)
      common/junkcc/dchar(511)
      l=511
      j=1
      k=nword
 1    k=k-511
      if(k)10495,2,2
10495 l=k+511
c2    call cmoved(z(j),dchar,l)
 2    do 3 loop=1,l
 3    read(z(j+loop-1),'(a8)') dchar(loop)
      call put(dchar,l,num3)
      j=j+511
      if(k)10500,10500,1
10500 return
      end
      subroutine rdchrs(z,nword,num3)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/junkcc/dchar(511)
      dimension z(*)
      j=1
10491 call find(num3)
      call get(dchar,l)
c     call dmovec(dchar,z(j),l)
      do 3 loop=1,l
 3    write(z(j+loop-1),'(a8)') dchar(loop)
      j=j+l
      if(j.le.nword)go to 10491
      if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present ', nword, j-1
       call caserr('invalid number of words in rdchrs')
      endif
      return
      end
      subroutine setstc(nf,ztext,zlab)
      character *(*) ztext,zlab
      dimension zlab(*)
      if(nf.eq.0)go to 2
      do 1 i=1,nf
 1    zlab(i)=ztext
 2    return
      end

      function eq(i,j)
c alliant compiler only compares lowest 32 bits of i8 quantities
c so have to use eq if top 32 bits may be significant
      logical eq
      integer i(2),j(2)
      eq = (i(1).eq.j(1)) .and. (i(2).eq.j(2))
      end

      function locatz(label,nf,itext)
      implicit real*8  (a-h,p-z),integer    (i-n)
      integer label(2,*),itext(2)
      logical eq
      do 1 i=1,nf
      if(eq(label(1,i),itext))go to 2
 1    continue
      locatz=0
      return
 2    locatz=i
      return
      end
      subroutine hstarg(trans,coul,exch,dens,nopk,
     + nclf,nscmf,iscmf,l1,l2,nprint)
c
c     ----- -hstarg- forms multiple fock matrices 
c                    doing only one integral pass -----
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
      common/craypk/integ(680),ipad(680)
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      common/blkin/gout(510),nint
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension coul(*),exch(*),dens(*)
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
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/vinteg/q(1)
      dimension trans(l1,*)
      data pt5,two/0.5d0,2.0d0/
c
      nshell=nclf+nscmf
c     load up the density matrices
      ioff = 1
      if (nclf .gt. 0) then
c      closed shell part
       call dengvb(trans,dens,1,iky,l1)
       ioff=ioff+l2
      endif
c     open shell part
      if (nscmf .gt. 0) then
       ilow = iscmf
       ihi = ilow + nscmf - 1
       do 200 i = ilow,ihi
       call dengvb(trans,dens(ioff),i,iky,l1)
 200   ioff = ioff + l2
      endif
c
c     ----- divide the diagonal density elements by 2
c
      ioff=0
      do 300 i = 1,nshell
      if (nprint.eq.5) then
       write(iwr,*)' density matrix, shell = ',i
       call prtri(dens(ioff+1),l1)
      endif
      do 1 j = 1,l1
      nij = ikyp(j) + ioff
   1  dens(nij)=dens(nij)*pt5
 300  ioff=ioff+l2
c
      if (nopk .ne. 1) go to 380
c
c      integrals in non-supermatrix form
c
      ndzero=nshell*l2
      call vclr(coul,1,ndzero)
      call vclr(exch,1,ndzero)
c
      do 2000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 2000
c
       if (o255i) then
       call sgmatn(l2,nshell,coul,exch,dens,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       else
       call sgmatn_255(l2,nshell,coul,exch,dens,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
1050   imemo = imemo + 512
      else
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 2004  jblock=jblock+1
       call get(gout,nw)
       if(nw)2001,2000,2001
 2001  if(jblock.ne.0)call find(main)
c
       if(o255i) then
          call sgmat(l2,nshell,coul,exch,dens)
       else
          call sgmat_255(l2,nshell,coul,exch,dens)
       endif
c
       if(jblock)2004,2000,2004
      endif
 2000 continue
c
      if(nclf.gt.0) then
       if (nprint.eq.5) then
        write(iwr,*)' closed shell J matrix'
        call prtri(coul,l1)
        write(iwr,*)' closed shell K matrix'
        call prtri(exch,l1)
       endif
c
       call dscal(l2,two,coul,1)
       call vsub(coul,1,exch,1,coul,1,l2)
c
      endif
      return
c
c     ----- integrals in supermatrix form -----
c
  380 continue
      labss = num2ejk + num2ejk + 1
      ndzero=nshell*l2
      call vclr(coul,1,ndzero)
      call vclr(exch,1,ndzero)
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
c     memory resident integrals
c
       imemo = iqqoff(main) + 1
       do 2020 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 1000
c
       call jksupm(coul,coul(l2+1),exch(l2+1),dens,dens(l2+1),nscmf,
     + l2, q(imemo), q(imemo+num2ejk), q(imemo+num2ejk+num2ejk), 
     +     q(imemo+510) )
2020   imemo = imemo + 512
      else
c
c     file resident integrals
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gout,nw)
       if(nw)1001,1000,1001
1001   if(jblock.ne.0)call find(main)
c
c ----------build the fock matrices--------------
c
       iword=1
       if (nclf.ne.0.and.nscmf.eq.0) then
c
c---------closed shells only------------
c
c     unpack a buffer of labels
c
        call unpack(gout(labss),lab1632,integ,numlabjk)
        do 480 m = 1,nint
        ij = integ(iword  )
        kl = integ(iword+1)
        goutpm = gout(m)
        coul(ij) = coul(ij) + goutpm*dens(kl)
        coul(kl) = coul(kl) + goutpm*dens(ij)
        iword=iword+2
 480    continue
c      
       else if (nclf.eq.0.and.nscmf.ne.0) then
c
c---------open shells only----------------
c
c     unpack a buffer of labels
c
       call unpack(gout(labss),lab1632,integ,numlabjk)
        do 500 m = 1,nint
        ij = integ(iword  )
        kl = integ(iword+1)
        goutpm = gout(m)
        goutkm = gout(num2ejk+m)
        do 501 i = 1,nscmf
        coul(ij) = coul(ij) + goutpm*dens(kl)
        exch(ij) = exch(ij) + goutkm*dens(kl)
        coul(kl) = coul(kl) + goutpm*dens(ij)
        exch(kl) = exch(kl) + goutkm*dens(ij)
        ij = ij + l2
 501    kl = kl + l2
        iword=iword+2
 500    continue
c
       else
c
c---------both closed + open shells------------
c
       call jksupe(coul,coul(l2+1),exch(l2+1),dens,dens(l2+1),nscmf,l2)
      endif
c
c-----------check next buffer-------------------
c
       if(jblock) 1004,1000,1004
      endif
 1000 continue
c
      if(nscmf.gt.0) then
       ioff = 1
       if(nclf.ne.0) ioff = ioff + l2
       do 720 ifo = 1,nscmf
       call vadd(coul(ioff),1,exch(ioff),1,coul(ioff),1,l2)
       call dscal(l2,pt5,coul(ioff),1)
 720   ioff=ioff+l2
      endif
      return
      end
      subroutine hstar(d,f,exch,nopk)
c
c     ----- -hstar- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order-
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
c     ----- ***********************************  ------
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension d(*),f(*),exch(*)
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
      logical oexch, ocoul, onodft
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
      common/blkin/gin(510),nint
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
      common/vinteg/q(1)
c
      data pt5/0.5d0/
c
      call vclr(f,1,nx)
c
      do 1 m=1,num
      nij=ikyp(m)
   1  d(nij)=d(nij)*pt5
      if (nopk.ne.1) go to 280
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
c    strategy a function of whether modified dft fock matrices are
c    required - if so, abandon the vectorised/cal strategy and
c    fall back to single-block fortran processing
c
      odft = CD_active()
      onodft = .not. odft
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
      ovect = o255i .and. .not. odft
      call vclr(exch,1,nx)
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
       call upack2(q(imemo+511),ibll,nw)  
       if(nw.eq.0) go to 1000
c
       if (o255i) then
          call sgmtmm(f,d,q(imemo),q(imemo+num2e),q(imemo+510) )
       else
          call sgmtmm_255(f,d,q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
c
 1050  imemo = imemo + 512
      else
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gin,nw)
       if(nw.eq.0)go to 1005
       if(jblock.ne.0)call find(main)
c
       if (ovect) then
          call sgmata(f,d)
       else if(.not.o255i) then
          call sgmata_255(f,d)
       else
          call sgmata_dft(f,d,facex,ocoul,oexch)
       endif
c
       if(jblock)1004,1000,1004
 1005  llblk(ifile)=liblk(ifile)-iposun(main)+1
      endif
 1000 continue
c
      go to 360
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c
 280  continue
c
       do 2000 ifile=1,lfile
       main=lotape(ifile)
       if(omem(main)) then
        imemo = iqqoff(main) + 1
        if(nopk.ne.-1.and.osortp) then
         do 2010 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsuppm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2010     imemo = imemo + 512
c
        else if(nopk.ne.-1.and..not.osortp) then
c
         do 2020 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsupm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2020     imemo = imemo + 512
c
        else 
c
         do 2030 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsupam(f,d,q(imemo),q(imemo+num2ejk+num2ejk),
     +               q(imemo+510) )
2030     imemo = imemo + 512
        endif
c
       else
c
        call search(liblk(ifile),main)
        call find(main)
        jblock=llblk(ifile)
 2004   jblock=jblock+1
        call get(gin,nw)
        if(nw)2001,2005,2001
 2001   if(jblock.ne.0)call find(main)
c
        if(nopk.eq.-1)call gsupa(f,d)
c
        if(nopk.ne.-1) then
          if(osortp) then
            call gsupp(f,d)
          else
            if(o255i) then
               call gsup(f,d)
            else
               call gsup_255(f,d)
            endif
          endif
        endif
c
        if(jblock)2004,2000,2004
 2005   llblk(ifile)=liblk(ifile)-iposun(main)+1
       endif
 2000  continue
c
c ----
c
 360  continue
      do 3   m = 1,num
      nij = ikyp(m)
   3  d(nij) = d(nij)+d(nij)
      call dscal(nx,pt5,f,1)
      return
      end
      subroutine hstaru(da,db,fa,fb,nopk)
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
      common/craypk/integ(680),ipad(680)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/blkin/gout(510),nint
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      logical oexch, ocoul, onodft
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
      common/vinteg/q(1)
      dimension da(*),fa(*),db(*),fb(*)
      equivalence (gg1,gg),(ikyj,jl),(ikyk,jk),(ikyi,il)
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
      data pt5/0.5d0/
c     data two/2.0d0/
      do 100 m = 1,nx
      duma=da(m)
      dumb=db(m)
      da(m)=duma+dumb
 100  db(m)=duma-dumb
      call vclr(fa,1,nx)
      call vclr(fb,1,nx)
      do 1   m=1,num
      nij=ikyp(m)
      da(nij)=da(nij)*pt5
   1  db(nij)=db(nij)*pt5
      if(nopk.ne.1)goto280
c
c     ----- integrals are not in supermatrix form (nopk=.true.)
c
      odft = CD_active()
      onodft = .not. odft
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 1000
c
       if (o255i) then
       call proc2m(fa,da,fb,db,
     +      facex,ocoul,oexch,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       else
       call proc2m_255(fa,da,fb,db,
     +      facex,ocoul,oexch,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
1050   imemo = imemo + 512
      else
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gout,nw)
       if(nw)1001,1000,1001
1001   if(jblock.ne.0)call find(main)
c
       if (o255i) then
          call proc2f(fa,da,fb,db,facex,ocoul,oexch)
       else
          call proc2_255 (fa,da,fb,db,facex,ocoul,oexch)
       endif
c
       if(jblock)1004,1000,1004
c
      endif
 1000 continue
      do 700 i=1,nx
      duma=fa(i)
      dumb=fb(i)
      fa(i)=duma-dumb
 700  fb(i)=duma+dumb
c
      go to 5
c
c     ----- integrals are in a supermatrix form (nopk=.false.)
c
 280  continue
c
      do 2000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
       imemo = iqqoff(main) + 1
       do 2050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 2000
c
       call gsupum(fa,fb,da,db,
     +     q(imemo),q(imemo+num2ejk),q(imemo+num2ejk+num2ejk),
     +     q(imemo+510) )
2050   imemo = imemo + 512
      else
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 2004  jblock=jblock+1
       call get(gout,nw)
       if(nw)2001,2000,2001
2001   if(jblock.ne.0)call find(main)
       call unpack(gout(num2ejk+num2ejk+1),lab1632,
     +             integ,numlabjk)
       iword=1
       do 360 m=1,nint
       nij1=integ(iword  )
       nkl1=integ(iword+1)
       valp1 = gout(m)
       valk1 = gout(num2ejk+m)
       dump = valp1*da(nkl1)
       dumk = valk1*db(nkl1)
       fa(nij1) = fa(nij1)+dump-dumk
       fb(nij1) = fb(nij1)+dump+dumk
       dump = valp1*da(nij1)
       dumk = valk1*db(nij1)
       fa(nkl1) = fa(nkl1)+dump-dumk
       fb(nkl1) = fb(nkl1)+dump+dumk
       iword=iword+2
 360   continue
c
       if(jblock)2004,2000,2004
      endif
c
 2000 continue
    5 continue
      do   6 m = 1,num
      nij = ikyp(m)
      da(nij) = da(nij)+da(nij)
    6 db(nij) = db(nij)+db(nij)
      do 420 m = 1,nx
      duma = da(m)
      dumb = db(m)
      da(m) = (duma+dumb)*pt5
 420  db(m) = (duma-dumb)*pt5
c
      call dscal(nx,pt5,fa,1)
      call dscal(nx,pt5,fb,1)
c
      return
      end
      subroutine gsupam(h,p,g,integ,mword)
      implicit real*8  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer *2 integ
      dimension h(*),p(*),g(*),integ(2,*),mword(*)
      do 1 iw=1,mword(1)
      ij=integ(1,iw)
      kl=integ(2,iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
      return
      end
      subroutine gsupm(f,d,g,integ,mword)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *2 integ
      dimension g(*),integ(2,*),mword(*)
      dimension f(*),d(*)
      do 1 iw=1,mword(1)
      ij=integ(1,iw)
      kl=integ(2,iw)
      f(ij)=d(kl)*g(iw)+f(ij)
      f(kl)=d(ij)*g(iw)+f(kl)
    1 continue
      return
      end

      subroutine gsuppm(h,p,g,ijkl,mword)
      implicit real*8  (a-h,o-z)
      integer*4 i8temp(2),fiv11
      integer *2 ijkl
      dimension h(*),p(*),g(*),ijkl(2,*),mword(2)
      equivalence (i8temp,r8temp)
      data fiv11/511/

      if(mword(1).le.0) return
c
c
      iwlo = 1
 10   continue
      r8temp=g(iwlo)
c
      nij=and(fiv11,i8temp(1))  
      ij = ijkl(1,iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword(1)) call caserr('invalid length in gsuppm')
      if(ij.ne.ijkl(1,iwhi)) call caserr('ij mismatch in gsuppm')
c
c fortran source
c
       temp=0.0d0
       pij = p(ij)
       do 20 iw = iwlo,iwhi
           kl = ijkl(2,iw)

           h(kl) = h(kl) + pij * g(iw)
 20        temp = temp + p(kl) * g(iw)
       h(ij) = h(ij) + temp
c
      iwlo = iwhi + 1
      if(iwlo.le.mword(1)) goto 10
      return
      end


      subroutine gsupum(fa,fb,da,db,goutp,goutk,integ,nint)
      implicit real*8  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer *2 integ
      dimension nint(*)
      dimension goutp(*),goutk(*)
      dimension fa(*),fb(*),da(*),db(*),integ(*)
c
       iword=1
       do 360 m=1,nint(1)
       nij1=integ(iword  )
       nkl1=integ(iword+1)
       valp1=goutp(m)
       valk1 = goutk(m)
       dump = valp1*da(nkl1)
       dumk = valk1*db(nkl1)
       fa(nij1) = fa(nij1)+dump-dumk
       fb(nij1) = fb(nij1)+dump+dumk
       dump = valp1*da(nij1)
       dumk = valk1*db(nij1)
       fa(nkl1) = fa(nkl1)+dump-dumk
       fb(nkl1) = fb(nkl1)+dump+dumk
       iword=iword+2
 360   continue
       return
       end



      subroutine jksupm(hscm0,hscm1,hscm2,dscm0,dscm1,nscmf,l2,
     + goutp, goutk, integ, nint)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer *2 integ
      dimension hscm0(*),hscm1(*),hscm2(*),dscm0(*),dscm1(*)
      dimension goutp(*),goutk(*)
      dimension nint(2),integ(2,*)
c
c---------both closed + open shells------------
c
      if(nscmf.eq.1) then
c
       do   10 m = 1,nint(1)
       ij = integ(1,m)
       kl = integ(2,m)
       goutpm = goutp(m)
       goutkm = goutk(m)
       hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
       hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
       hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
       hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
       hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
       hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
   10  continue
      else
       do   1 m = 1,nint(1)
       ij = integ(1,m)
       kl = integ(2,m)
       goutpm = goutp(m)
       goutkm = goutk(m)
       hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
       hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
       do 2   i = 1,nscmf
       hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
       hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
       hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
       hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
       ij = ij + l2
    2  kl = kl + l2
    1  continue
      endif
      return
      end

      subroutine proc2m(fock,p,ak,q,facex,ocoul,oexch,
     +                  gg,integ,mword)

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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension ak(*),q(*),p(*),fock(*),gg(*)

      dimension mword(*)
      logical *1 integ
      dimension integ(*)
      logical ocoul,oexch
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

c  variant uses logical*1 on RHS of unpack
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
c
      iword=1
      do 6000 iw=1,mword(1)
       jlog(1) = integ(iword  )
       ilog(1) = integ(iword+1)
       llog(1) = integ(iword+2)
       klog(1) = integ(iword+3)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
 6000 iword=iword+4
      return
      end

      subroutine proc2m_255(fock,p,ak,q,facex,ocoul,
     +                  oexch,
     +                  gg,integ,mword)

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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension ak(*),q(*),p(*),fock(*),gg(*)

      dimension mword(*)
      integer *2 integ
      dimension integ(*)
      logical ocoul,oexch
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

c  variant uses logical*1 on RHS of unpack
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
c
      iword=1
      do 6000 iw=1,mword(1)
       jlog(1) = integ(iword  )
       ilog(1) = integ(iword+1)
       llog(1) = integ(iword+2)
       klog(1) = integ(iword+3)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
 6000 iword=iword+4
      return
      end

      subroutine sgmatn(l2,nshell,b,c,p,gg,intin,mword)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension b(*),c(*),p(*),gg(*)
      dimension mword(*)
      logical *1 intin
      dimension intin(*)
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
      int4=1
      if(nshell.gt.1) then
c
      do 6 iw=1,mword(1)
       jlog(1)=intin(int4  )
       ilog(1)=intin(int4+1)
       llog(1)=intin(int4+2)
       klog(1)=intin(int4+3)
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
    6 int4=int4+4
      else
      do 66 iw=1,mword(1)
      jlog(1)=intin(int4  )
      ilog(1)=intin(int4+1)
      llog(1)=intin(int4+2)
      klog(1)=intin(int4+3)
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
   66 int4=int4+4
      endif
      return
      end
      subroutine sgmatn_255(l2,nshell,b,c,p,gg,intin,mword)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension b(*),c(*),p(*),gg(*)
      dimension mword(*)
      integer *2 intin
      dimension intin(*)
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
      int4=1
      if(nshell.gt.1) then
c
      do 6 iw=1,mword(1)
       jlog(1)=intin(int4  )
       ilog(1)=intin(int4+1)
       llog(1)=intin(int4+2)
       klog(1)=intin(int4+3)
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
    6 int4=int4+4
      else
      do 66 iw=1,mword(1)
      jlog(1)=intin(int4  )
      ilog(1)=intin(int4+1)
      llog(1)=intin(int4+2)
      klog(1)=intin(int4+3)
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
   66 int4=int4+4
      endif
      return
      end
      subroutine sgmtmm(fock,p,gg,integ,mword)
      implicit real*8  (a-h,o-z)
      dimension p(*),fock(*),gg(*)
      dimension mword(*)
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
      logical *1 integ 
      dimension integ(*)
      logical oexch, ocoul, onodft
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
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
c
      iword = 1
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
      do 6000 iw=1,mword(1)
      jlog(1) = integ(iword  )
      ilog(1) = integ(iword+1)
      llog(1) = integ(iword+2)
      klog(1) = integ(iword+3)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
 6000 iword=iword+4
      return
      end
      subroutine sgmtmm_255(fock,p,gg,integ,mword)
      implicit real*8  (a-h,o-z)
      dimension p(*),fock(*),gg(*)
      dimension mword(*)
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
      integer *2 integ 
      dimension integ(*)
      logical oexch, ocoul, onodft
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
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
c
      iword = 1
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
      do 6000 iw=1,mword(1)
      jlog(1) = integ(iword  )
      ilog(1) = integ(iword+1)
      llog(1) = integ(iword+2)
      klog(1) = integ(iword+3)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
      if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
 6000 iword=iword+4
      return
      end
      subroutine ver_util2(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util2.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
