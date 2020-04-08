c 
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mrdci4.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   Table-ci (analysis module) =
c ******************************************************
c ******************************************************
      subroutine nmrdci(q,lword)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
      common/blkorbs/value(maxorb),occ(maxorb+1),iorbs,
     *newbas,ksum,ivalue,iocc,ipad
      common/ftape/ispac(12),ktape,ifox,idum4(6)
      common/scrtch/nconf(276)
c
      dimension q(*)
c
      write(iwr,3)yed(idaf),ibl3d
 3    format(/1x,104('=')//
     *40x,36('*')/
     *40x,'Table-ci  --  natural orbital module'/
     *40x,36('*')//
     *1x,'dumpfile on ',a4,' at block',i6/)
       if(ibl3d)12,12,13
 12   call caserr('invalid starting block for dumpfile')
 13   write(iwr,18)lword
 18   format(/' main core available = ',i8,' words'//)
      if(num.le.0.or.num.gt.maxorb) call caserr(
     * 'invalid number of basis functions')
      ivalue=-1
      iocc=1
      call rewftn(ifox)
      read(ifox)iwod,vnuc,zero,imo,m,nconf,ksum,iorbs
      if(iorbs.ne.num.or.imo.gt.ksum.or.ksum.gt.iorbs)
     * call caserr('inconsistent parameters on dumpfile')
      call rewftn(ifox)
      newbas=iorbs
      len1=ksum*iorbs
      nx=imo*(imo+1)/2
      len3=imo*imo
      i1=1
      i2=i1+len1
      i3=i2+len1
      i4=i3+nx
      i5=i4+len3
      if(i5.gt.lword)call caserr('insufficient memory available')
      lwor= lword-i5-1
      write(iwr,19)lwor
 19   format(' main core not used = ',i8,' words '/)
      call nmrd0(q(i1),q(i2),q(i3),q(i4),q(i5),len1,nx,len3)
      call clredx
      call rewftn(ifox)
      call rewftn(ktape)
c
      return
      end
      subroutine eigene(a,r,n,mv)
      implicit real*8  (a-h,o-z), integer (i-n)
      dimension a(*),r(*)
      if(mv-1) 10,25,10
  10  iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0d0
      if(i-j) 20,15,20
  15  r(ij)=1.0d0
  20  continue
  25  anorm=0.0d0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
  30  ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
  35  continue
      if(anorm) 165,165,40
  40  anorm=1.414d0*dsqrt(anorm)
      anrmx=anorm*1.0d-6/dble(n)
      ind=0
      thr=anorm
  45  thr=thr/dble(n)
  50  l=1
  55  m=l+1
  60  mq=((m*m)-m)/2
      lq=((l*l)-l)/2
      lm=l+mq
      if(dabs(a(lm))-thr) 130,65,65
  65  ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
  70  y=-y
  75  if(y .gt. 1.0d0) y=1.0d0
      if(y .lt. -1.0d0) y=-1.0d0
      sinx=y/dsqrt(2.0d0*(1.0d0+(dsqrt(1.0d0-y*y))))
      sinx2=sinx*sinx
      cosx=dsqrt(1.0d0-sinx2)
      cosx2=cosx*cosx
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
  80  if(i-m) 85,115,90
  85  im=i+mq
      go to 95
  90  im=m+iq
  95  if(i-l) 100,105,105
 100  il=i+lq
      go to 110
 105  il=l+iq
 110  x=a(il)*cosx-a(im)*sinx
      a(im)=(a(il)*sinx)+(a(im)*cosx)
      a(il)=x
 115  if(mv-1) 120,125,120
 120  ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=(r(ilr)*sinx)+(r(imr)*cosx)
      r(ilr)=x
 125  continue
      x=2.0d0*a(lm)*sincs
      y= (a(ll)*cosx2)+(a(mm)*sinx2)-x
      x=(a(ll)*sinx2)+(a(mm)*cosx2)+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
 130  if(m-n) 135,140,135
 135  m=m+1
      go to 60
 140  if(l-n+1) 145,150,145
 145  l=l+1
      go to 55
 150  if(ind-1) 160,155,160
 155  ind=0
      go to 50
 160  if(thr-anrmx) 165,165,45
 165  iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
 170  x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
 175  do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
 180  r(imr)=x
 185  continue
      return
      end
      subroutine nmrd0(qat,y,gx,w,q,len1,len2,len3)
      implicit real*8  (a-h,o-z), integer (i-n)
      character *1 dash,star
      character *4 tagg
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
      common/blkorbs/vvv(maxorb),occat(maxorb+1),nspabc(6)
      common/craypk/itag(mxroot),isec(mxroot),jsec(mxroot),nwi
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scrtch/imat(45000),jmat(3780),f(mxcrec),
     *irm(mxcrec),h(100),ibal(8),itil(8),
     *mcomp(256),e(7401),c(400),nzer(256),
     *tc(255),rr(255),nfil(8),ndog(20),
     *g(mxcrc2),gj(mxcrc2),jrm(126),isecc(18),iper(100),ihog(48),
     *fj(mxcrec)
      common/bufb/nconf(5),lsym(2040),ncomp(256),
     *imap(504),jmap(140),ikan(mxcrc2),jkan(mxcrc2),
     *mper(maxorb),intp(maxorb),intn(maxorb),jcon(maxorb),
     +nda(20),mj(20),ided(20),nj(8),kj(8),ntil(8),nbal(9),
     +idep(8),lj(8),nstr(5),nytl(5),nplu(5),ndub(5),
     +jtest(mxnshl),ktest(mxnshl),ibug(4)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/ftape/
     * lscr,ispaci(11),ktape,ifox,idum4(6)
      common/lsort/ical,iswh,m,mu,m2,imo
c
      dimension qat(len1),w(len3),y(len1),gx(len2),q(len3)
      dimension cf(22500),icf(22500)
      dimension tagg(2),ilifs(255),newya(255)
c
      equivalence (cf(1),imat(1))
      equivalence (imat(22501),icf(1))
c
      data two/2.0d0/
      data tagg/'sabf','a.o.'/
      data dash,star/'-','*'/
c
 9500 format(/' *** commence natural orbital analysis at ',f8.2,
     * ' secs.')
      read(ifox)iwod,vnuc,zero,imo,m,nconf,newya,
     + nytl,nplu,ndub,iswh,ksum,iorbs,
     + knu,cf,icf,c,ibal,itil,mcomp,kj,lj,n,ifrk,e,lsym,nsel,nj,
     + ntil,nbal,ncomp
c
c     if (oprint(32)) then
c      call print_ifox(y,len1,iky,
c    + iwod,vnuc,zero,imo,m,nconf,newya,
c    + nytl,nplu,ndub,iswh,ksum,iorbs,
c    + knu,cf,icf,c,ibal,itil,mcomp,kj,lj,n,
c    + ifrk,e,lsym,nsel,nj,ntil,nbal,ncomp)
c     endif
c
c     specify default settings ....
c     density matrix o/p to ft42
c     natural orbital generation
c     abelian point groups only
c
      lwater=1
      iprp=1
      ndeg=0
c
      if(nwi.le.0.or.nwi.gt.mxroot)call caserr(
     *'invalid number of orbital sets requested')
      do 7002 i=1,nwi
      if(itag(i).le.0.or.itag(i).gt.1000)call caserr(
     *'invalid orbital set requested')
 7002 continue
      write(iwr,8000)nwi,(itag(i),i=1,nwi)
 8000 format(/
     *' *** natural orbital analysis requested for',i3,
     *' ci vectors'//5x,
     *'with following locations on the ci vector file (ft36)'
     *,20i3)
      if(oprint(32))
     *write(iwr,8501)
 8501 format(/
     *' print of density matrix and n.o.s in mo basis requested'/)
      do 20006 i=1,nwi
      if(isec(i).gt.0)write(iwr,7004)tagg(1),i,isec(i)
      if(jsec(i).gt.0)write(iwr,7004)tagg(2),i,jsec(i)
      if(isec(i).gt.350. or .jsec(i).gt.350)
     * call caserr('invalid dumpfile section specified for nos.')
20006 continue
 7004 format(/
     *' ** route n.o.s ( ',a4,' basis) for state',i3,
     *' to section',i4,' of dumpfile')
      write(iwr,7053)imo,m
 7053 format(//
     *' no. of active orbitals ',i3/
     *' no. of active electrons',i3/)
      write(iwr,3)(dash,i=1,40)
 3    format(
     *' irreducible       no. of'/
     *' representation    active orbitals'/1x,40a1)
      do 7060 i=1,n
 7060 write(iwr,7061)i,lj(i)
 7061 format(7x,i5,10x,i3)
      write(iwr,7071)(dash,i=1,40)
 7071 format(1x,40a1)
       jx=0
       do 95 i=1,iorbs
       nzer(i)=jx
 95   jx=jx+ncomp(i)
      icore=ksum-imo
      icorp=icore-1
      im9=iky(imo+1)
      if(im9.ne.len2)call caserr(
     *'inconsistency detected in contents of table-ci interfaces')
      iimo=0
      do 7000 i=1,ksum
      ilifs(i)=iimo
 7000 iimo=iimo+iorbs
      ncol=ksum*iorbs
      jx=0
       mb=0
       ibuk=0
       mg=icore*iorbs
       do 303 i=1,n
       kp=kj(i)
       lp=kp-lj(i)
       kx=ntil(i)
       do 303 j=1,kp
       if(j.gt.lp) go to 360
       jz=ibuk
       id=ibuk+1
       ibuk=ibuk+iorbs
       jw=ibuk
       go to 361
 360   jz=mg
       id=mg+1
       mg=mg+iorbs
       jw=mg
361    do 304 k=id,jw
 304   y(k)=0
       mb=mb+1
       ll=mcomp(mb)
       do 305 k=1,ll
       jx=jx+1
       vm=cf(jx)
       nzh=icf(jx)+kx
       lx=nzer(nzh)
       nzh=ncomp(nzh)
       do 305 l=1,nzh
       lx=lx+1
       ibj=lsym(lx)
       if(ibj.gt.0) go to 306
       y(jz-ibj)=-vm
       go to 305
 306  y(jz+ibj)=vm
 305   continue
 303   continue
       if(iprp.eq.0) go to 308
       call rewftn(ktape)
       write(ktape)iorbs,imo,icore,n,lj,knu,m,y,e,c,newya
       if(lwater.eq.0) go to 20
308   if(ndeg.eq.0) go to 6
      call input
      do 7005 i=1,n
7005  call inpi(idep(i))
      write(iwr,2) (idep(i),i=1,n)
 2    format(5x,25i5)
      ib=icorp*iorbs
      kx=0
      ideg=0
      mg=0
      do 7 i=1,n
      ldj=lj(i)
      id=idep(i)
      ideg=ideg+id
      if (id.eq.1) go to 9
      call input
      do 7006 j=1,id
 7006 call inpi(ibug(j))
      do 7007 j=1,ldj
 7007 call inpi(iper(j))
      jx=0
      lx=kx
      kb=0
      do 8 j=1,id
      mg=mg+1
      ibj=ibug(j)
      mj(mg)=ibj
      nda(mg)=i
      do 8 k=1,ibj
      jx=jx+1
      jy=iper(jx)
      kx=kx+1
      jz=jy+lx
      mper(kx)=jz
      jb=ib+jy*iorbs
      do 8 l=1,iorbs
      jb=jb+1
      kb=kb+1
8     qat(kb)=y(jb)
      mb=ib+iorbs
      do 11 j=1,kb
      mb=mb+1
11    y(mb)=qat(j)
      go to 7
9     do 10 j=1,ldj
      kx=kx+1
10    mper(kx)=kx
      mg=mg+1
      mj(mg)=ldj
      nda(mg)=i
7     ib=ib+ldj*iorbs
      write(iwr,2) (mper(i),i=1,imo)
      write(iwr,2) (mj(i),i=1,ideg)
      call input
      do 7008 i=1,ideg
 7008 call inpi(ided(i))
      write(iwr,2) (ided(i),i=1,ideg)
      go to 20
6     do 14 i=1,imo
14    mper(i)= i
      do 15 i=1,n
      nda(i)=i
15    mj(i)=lj(i)
20    mu=m-2
      m2=m+2
      nzh=0
      cpu=cpulft(1)
      write(iwr,9500)cpu
      do 1000 i=1,nwi
      write(iwr,9000)(star,j=1,129)
 9000 format(/1x,129a1)
      write(iwr,8502)i
 8502 format(/40x,
     *'natural orbital analysis for state no.',i3/40x,41('-')/)
      call vclr(qat,1,ncol)
      iimo=0
      nid=itag(i)
30    nzh=nzh+1
      if (nzh.ne.nid) go to 31
      call rewftn(lscr)
      do 23 j=1,iswh
      if (nconf(j).eq.0) go to 23
      read (ifox) nhb,imax,ndt,kml,imap,ihog
      write (lscr)nhb,imax,ndt,kml,imap,ihog
      do 610 k=1,nhb
      read (ifox) jkan,fj,g,h
  610 write (lscr)jkan,fj,g,h
   23 continue
      call rewftn(lscr)
      go to 32
   31 do 21 j=1,iswh
      if (nconf(j).eq.0) go to 21
      read (ifox) nhb
      do 22 k=1,nhb
   22 read (ifox)
   21 continue
      go to 30
   32 call vclr(q,1,im9)
      call nmrd1(q)
      if(.not.oprint(32))go to 8503
      write(iwr,511) i
 511  format(/20x,
     *'first order ci density matrix for state no.',i2,
     *'  (active mo basis)'/)
      call writel(q,imo)
      write(iwr,9000)(dash,j=1,129)
8503  if (iprp.eq.0) go to 501
      write(ktape) im9,q
      if(lwater.eq.0) go to 1000
 501  if(oprint(32))
     *write(iwr,451)i
 451  format(/20x,
     *' *** natural orbitals (active mo basis) for state',
     *i2,' ****'//)
      call rewftn(lscr)
      do 315 j=1,n
 315   nfil(j)=0
      if(icore.eq.0) go to 450
      lw=-iorbs
      do 322 k=1,n
      ix=kj(k)-lj(k)
      if(ix.eq.0) go to 322
      iz=ntil(k)
      lx=nj(k)
      ky=nbal(k)+1
      nfil(k)=nfil(k)+ix
      do 324 l=1,ix
      lw=lw+iorbs
      jz=iz
      ls=ky
      do 323 mb=1,lx
      jz=jz+1
      ks=lsym(ls)+lw
      ls=ls+ncomp(jz)
 323   rr(mb)=y(ks)
       iimo=iimo+1
       iii=iz+ilifs(iimo)
       call dcopy(lx,rr,1,qat(iii+1),1)
 324   continue
       occat(iimo)=two
 322  continue
450   if(ndeg.ne.0) go to 502
      do 503 j=1,n
503   ided(j)=j
      ideg=n
502   ix=0
      imom=0
      orb=0
      do 504 j=1,ideg
      ix=ix+1
      mjk=mj(j)
      if (ided(j).eq.ix) go to 505
      ix=ix-1
      imom=imom+mjk
      go to 504
505   jmom=imom
      ndog(1)=j
      imom=imom+mjk
      iz=(jmom+icorp)*iorbs
      nstr(1)=iz
      nbas=1
      lx=jmom
      idx=0
      do 506 l=1,mjk
      lx=lx+1
      kp=mper(lx)
      kp=iky(kp)
      kx=jmom
      do 506 k=1,l
      kx=kx+1
      lp=mper(kx)+kp
      idx=idx+1
506   gx(idx)=q(lp)
      if (j.eq.ideg) go to 507
      ja=j+1
      jmom=imom
      do 508 jb=ja,ideg
      if (ided(jb).eq.ix) go to 509
      jmom=jmom+mj(jb)
      go to 508
509   iz=(jmom+icorp)*iorbs
      nbas=nbas+1
      ndog(nbas)=jb
      nstr(nbas)=iz
      lx=jmom
      idx=0
      do 510 l=1,mjk
      lx=lx+1
      kp=mper(lx)
      kp=iky(kp)
      kx=jmom
      do 510 k=1,l
      kx=kx+1
      lp=mper(kx)+kp
      idx=idx+1
510   gx(idx)=gx(idx)+q(lp)
      jmom=jmom+mjk
508   continue
507   vm=0.0d0
      if (mjk.eq.1) go to 550
      call eigene(gx,w,mjk,0)
      ida=0
      idb=0
      do 520 k=1,mjk
      idb=idb+k
      orc=gx(idb)
      vm=vm+orc
      if(oprint(32) )write(iwr,521)(star,l=1,129)
     *,ix,k,orc
 521  format(/1x,129a1/
     * ' irrep. no. ',i1,3x,
     *' n.o. sequence no.',i3,5x,
     *' occupation ',f14.8//)
      idc=ida+1
      ida=ida+mjk
      if(oprint(32) )write(iwr,522) (w(l),l=idc,ida)
522   format(/5x,10f12.8)
      do 526 l=1,nbas
      iz=nstr(l)
      do 524 la=1,iorbs
      iz=iz+1
      jz=iz
      cz=0
      do 525 lb=idc,ida
      jz=jz+iorbs
525   cz=cz+w(lb)*y(jz)
524   tc(la)=cz
      ig=ndog(l)
      ia=nda(ig)
      id=nfil(ia)+1
      write(lscr)orc,(tc(la),la=1,iorbs)
      nfil(ia)=id
      lb=nj(ia)
      iw=ntil(ia)
      ig=nbal(ia)+1
      do 351 la=1,lb
      iw=iw+1
      ks=lsym(ig)
      ig=ig+ncomp(iw)
351    rr(la)=tc(ks)
c
      iimo=iimo+1
      iii=ilifs(iimo)+ntil(ia)
      call dcopy(lb,rr,1,qat(iii+1),1)
      occat(iimo)=orc
c
 526  continue
520   continue
      write(iwr,527)ix,vm
 527  format(/' **** irrep. no. ',i1,
     *10x,'**** sum of occupation numbers',f14.8/)
      orb=orb+vm
      if(oprint(32))write(iwr,9000)(dash,l=1,129)
      go to 504
550   orc=gx(1)
      orb=orb+orc
      if(oprint(32))write(iwr,521)(star,l=1,129), ix,mjk,orc
      do 551 ld=1,nbas
      la=nstr(ld)+iorbs
      lb=la+1
      la=la+iorbs
      ig=ndog(ld)
      ia=nda(ig)
      id=nfil(ia)+1
      nfil(ia)=id
      write(lscr)orc,(y(ll),ll=lb,la)
      la=lb-1
          lb=nj(ia)
          iw=ntil(ia)
       ig=nbal(ia)+1
           do 352 ir=1,lb
        iw=iw+1
               ks=lsym(ig)+la
        ig=ig+ncomp(iw)
352   rr(ir)=y(ks)
      iimo=iimo+1
      iii=ilifs(iimo)+ntil(ia)
      call dcopy(lb,rr,1,qat(iii+1),1)
      occat(iimo)=orc
551   continue
      write(iwr,527)ix,orc
504   continue
      write(iwr,553) orb
553   format(/2x,'total active electron sum',f20.8)
      if(isec(i).le.0)go to 9004
      call putqno(qat,zero,tagg(1),isec(i),1)
9004  call rewftn(lscr)
c
c   now output n.o.s in ao representation
c   o/p to lineprinter comprises active orbitals ** only **
c   symmetry ordered, and occupation ordered within each irrep
c   o/p to dumpfile , section jsec(i), includes core orbitals
c
      ibuk=0
      jcore=0
      mg=icore
      nact=0
      do 9003 ir=1,n
      kp=kj(ir)
      lp=kp-lj(ir)
      do 9003 ld=1,kp
      if(ld.gt.lp)go to 9001
      jcore=jcore+1
      iimo=ilifs(jcore)
      occat(jcore)=two
      do 9002 k=1,iorbs
 9002 qat(newya(k)+iimo)=y(ibuk+k)
      ibuk=ibuk+iorbs
      go to 9003
 9001 nact=nact+1
      mgg=ilifs(mg+nact)
      read(lscr)orc,(qat(newya(k)+mgg),k=1,iorbs)
      occat(mg+nact)=orc
 9003 continue
      write(iwr,9000)(dash,j=1,129)
      write(iwr,9600)i
 9600 format(/40x,'natural orbitals for state',i3,
     *'  (a.o. basis)'/40x,43('-')/)
      if(oprint(31))call prev(qat(icore*iorbs+1),occat(icore+1),
     *  nact,iorbs,iorbs)
      write(iwr,9601)
      mg1=mg+1
      mg2=mg+nact
      write(iwr,9602)(occat(k),k=mg1,mg2)
 9602 format(/10x,8f14.7)
 9601 format(/50x,'occupation numbers'/50x,18('-'))
      if(jsec(i).le.0)go to 1000
      call putqno(qat,zero,tagg(2),jsec(i),2)
1000  continue
      cpu=cpulft(1)
      write(iwr,8001)cpu
 8001 format(/1x,104('*')//
     *' *** end of natural orbital analysis at ',f8.2,
     * ' secs.'/)
      return
      end
      subroutine nmrd1(q)
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
      common/blkorbs/vvv(maxorb),occat(maxorb+1),nspabc(6)
      common/craypk/itag(mxroot),isec(mxroot),jsec(mxroot),nwi
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scrtch/imat(45000),jmat(3780),f(mxcrec),
     *irm(mxcrec),h(100),ibal(8),itil(8),
     *mcomp(256),e(7401),c(400),nzer(256),
     *tc(255),rr(255),nfil(8),ndog(20),
     *g(mxcrc2),gj(mxcrc2),jrm(126),isecc(18),iper(100),ihog(48),
     +fj(mxcrec)
      common/bufb/nconf(5),lsym(2040),ncomp(256),
     +imap(504),jmap(140),ikan(mxcrc2),jkan(mxcrc2),
     +mper(maxorb),intp(maxorb),intn(maxorb),jcon(maxorb),
     +nda(20),mj(20),ided(20),nj(8),kj(8),ntil(8),nbal(9),idep(8),lj(8),
     +nstr(5),nytl(5),nplu(5),ndub(5),jtest(mxnshl),ktest(mxnshl),
     +ibug(4)
      common/ftape/
     * lscr,ispaci(11),ktape,ifox,idum4(6)
      common/lsort/ical,iswh,m,mu,m2,imo
c
      dimension cf(22500),q(*),icf(22500)
c
      equivalence (cf(1),imat(1))
      equivalence (imat(22501),icf(1))
c
      ical=1
      do 43 j=1,iswh
      nc=nconf(j)
      if (nc.eq.0) go to 43
      read (lscr) jhb,jmax,nd,kml,imap,ihog
      nl=nytl(j)
      nmns=j-1
      nqns=j-2
      nps=nplu(j)
      ndb=ndub(j)
      ndbq=ndb+ndb
      nod=nmns+nps
      nod2=nod-m2
      jcl=0
      idel=nd*m
      do 100 na=1,jhb
      jcl=jcl+jmax
      if(jcl.gt.nc) jmax=jmax+nc-jcl
      read(lscr) ikan,fj,gj,h
      nxj=0
      k4=0
      k5=0
      do 44 k=1,jmax
      do 50 l=1,nl
      nxj=nxj+1
   50 jtest(l)=ikan(nxj)
      if(nmns.gt.0) go to 56
      irm(k)=1
      if (nps.eq.0) go to 57
      do 58 l=1,nps
      k4=k4+1
58    imat(k4)=jtest(l)
      if (ndb.eq.0) go to 44
57    kx=nps
      do 59 l=1,ndb
      kx=kx+1
      k4=k4+1
      lk=jtest(kx)
      imat(k4)=lk
      k4=k4+1
59    imat(k4)=-lk
      go to 44
56    jq=1
      lx=k4
      do 60 l=1,nd
      k5=k5+1
      irm(k5)=1
      jz=jq+nqns
      jy=imap(jq)
      do 61 kk=1,nod
      lx=lx+1
      if (kk.ne.jy)  go to 62
      if (jz.eq.jq) go to 63
      jq=jq+1
      jy=imap(jq)
63    imat(lx)=-jtest(kk)
      go to 61
62    imat(lx)=jtest(kk)
61    continue
      jq=jz+1
60    lx=lx+ndbq
      if (ndb.eq.0) go to 101
      kx=nod2+k4
      kk=nod
      do 64 l=1,ndb
      kk=kk+1
      lk=jtest(kk)
      kx=kx+2
      mx=kx
      do 64 ll=1,nd
      mx=mx+m
      imat(mx+1)=lk
64    imat(mx+2) =-lk
      k4=mx+2
      go to 44
 101  k4=lx
 44   continue
      if(ical.gt.1) go to 102
      ical=2
      go to 103
 102  call rewftn(lscr)
      if(j.eq.1) go to 104
      jw=j-1
      if(j.eq.2) go to 105
      j8=j-2
      do 106 k=1,j8
      mc=nconf(k)
      if(mc.eq.0) go to 106
      read(lscr) nhb
      do 107 l=1,nhb
 107  read(lscr)
 106  continue
 105  mc=nconf(jw)
      if(mc.eq.0) go to 108
      read(lscr) nhb,imax,ndj,kmj,jmap
      nlj=nytl(jw)
      jmns=jw-1
      jqns=jw-2
      jps=nplu(jw)
      jdb=ndub(jw)
      jdbq=jdb+jdb
      nodj=jps+jmns
      nodj2=nodj-m2
      icl=imax
      do 72 iw=1,mc
      if (icl.lt.imax) go to 37
      read (lscr) jkan,f
      icl=1
      nx=0
      if=0
      go to 38
   37 icl=icl+1
      if=if+ndj
   38 do 66 kw=1,nlj
      nx=nx+1
   66 ktest(kw) =jkan(nx)
      call setsto(imo,0,jcon)
      if(nodj.eq.0) go to 110
      do 111 l=1,nodj
      nt=ktest(l)
 111  jcon(nt)=1
      if(jdb.eq.0) go to 70
 110  kk=nodj
      do 113 l=1,jdb
      kk=kk+1
      nt=ktest(kk)
 113  jcon(nt)=2
70    if (jmns.gt.0) go to 80
      jrm(1)=1
      if (jps.eq.0) go to 81
      do 82 kw=1,jps
82    jmat(kw)=ktest(kw)
      if (jdb.eq.0) go to 94
81    lx=jps
      kx=jps
      do 84 kw=1,jdb
      kx=kx+1
      lx=lx+1
      lk=ktest(kx)
      jmat(lx)=lk
      lx=lx+1
84    jmat(lx)=-lk
      go to 94
80    lx=0
      jq=1
      do 85 kw=1,ndj
      jrm(kw)=1
      jz=jq+jqns
      jy=jmap(jq)
      do 86 lw=1,nodj
      lx=lx+1
      if (lw.ne.jy) go to 87
      if (jz.eq.jq) go to 88
      jq=jq+1
      jy=jmap(jq)
88    jmat(lx)=-ktest(lw)
      go to 86
87    jmat(lx)=ktest(lw)
86    continue
      jq=jz+1
85    lx=lx+jdbq
      if (jdb.eq.0) go to 94
      kx=nodj2
      kk=nodj
      do 89 kw=1,jdb
      kk=kk+1
      lk=ktest(kk)
      kx=kx+2
      mx=kx
      do 89 lw=1,ndj
      mx=mx+m
      jmat(mx+1)=lk
89    jmat(mx+2)=-lk
94    nxj=0
      inis=-nd
      igi=-kml
      do 114 l=1,jmax
      inis=inis+nd
      igi=igi+kml
      do 115 ll=1,nl
      nxj=nxj+1
 115  jtest(ll)=ikan(nxj)
      nix=0
      if(nod.eq.0) go to 116
      do 117 ll=1,nod
      jt=jtest(ll)
      if(jcon(jt).gt.0) go to 117
      if(nix.eq.1) go to 114
      nix=1
 117  continue
      if(ndb.eq.0) go to 119
 116  kk=nod
      do 120 ll=1,ndb
      kk=kk+1
      jt=jtest(kk)
      jb=jcon(jt)
      if(jb.eq.2) go to 120
      if(jb.eq.0.or.nix.eq.1) go to 114
      nix=1
 120  continue
 119  if(ndj.gt.kml) go to 150
      inj=-m
      ls=if
      do 151 lw=1,ndj
      ls=ls+1
      orb=f(ls)
      inj=inj+m
      do 152 kw=1,imo
      intp(kw)=0
152   intn(kw)=0
      ink=inj
      do 153 kw=1,m
      ink=ink+1
      lx=jmat(ink)
      if (lx.lt.0) go to 154
      intp(lx)=kw
      go to 153
154   intn(-lx)=kw
153   continue
      ks=igi
      do 155 kw=1,kml
      nix=0
      ks=ks+1
      ir=ihog(kw)+inis
      ini=(ir-1)*m
      ink=ini
      do 156 mw=1,m
      ink=ink+1
      lo=imat(ink)
      if (lo.lt.0) go to 157
      if (intp(lo).gt.0) go to 156
160   if (nix.eq.1) go to 155
      nix=1
      iodd=mw
      go to 156
157   if(intn(-lo).eq.0) go to 160
156   continue
      if (iodd.eq.1) go to 161
      kx=1
      ky=iodd-1
      ink=ini
169   do 162 mw=kx,ky
      ink=ink+1
      lx=imat(ink)
      if (lx.lt.0) go to 163
      lo=intp(lx)
      go to 164
163   lo=intn(-lx)
164   if (lo.eq.mw) go to 162
      ix=inj+lo
      iy=inj+mw
      ni=jmat(iy)
      jmat(ix)=ni
      jmat(iy)=lx
      jrm(lw)=-jrm(lw)
      if (lx.lt.0) go to 165
      intp(lx)=mw
      go to 166
165   intn(-lx)=mw
166   if (ni.lt.0) go to 167
      intp(ni)=lo
      go to 162
167   intn(-ni)=lo
162   continue
      if (ky.gt.mu) go to 170
161   kx=iodd+1
      ky=m
      ink=ini+iodd
      go to 169
170   vm=orb*gj(ks)
      if (irm(ir).ne.jrm(lw))vm=-vm
      ink=iodd+ini
      ia=imat(ink)
      ink=iodd+inj
      ib=jmat(ink)
      if (ia.gt.0) go to 171
      ia=-ia
      ib=-ib
171   ia=min(ia,ib)+iky(max(ia,ib))
      q(ia)=q(ia)+vm
155   continue
151   continue
      go to 114
150   ks=igi
      do 180 lw=1,kml
      ks=ks+1
      ir=ihog(lw)+inis
      ini=(ir-1)*m
      orb=gj(ks)
      do 181 kw=1,imo
      intp(kw)=0
181   intn(kw)=0
      ink=ini
      do 182 kw=1,m
      ink=ink+1
      lx=imat(ink)
      if (lx.lt.0) go to 183
      intp(lx)=kw
      go to 182
183   intn(-lx)=kw
182   continue
      inj=-m
      ls=if
      do 184 kw=1,ndj
      ls=ls+1
      inj=inj+m
      ink=inj
      nix=0
      do 185 mw=1,m
      ink=ink+1
      lo=jmat(ink)
      if (lo.lt.0) go to 186
      if (intp(lo).gt.0) go to 185
187   if(nix.eq.1) go to 184
      nix=1
      iodd=mw
      go to 185
186   if(intn(-lo).eq.0) go to 187
185   continue
      if (iodd.eq.1) go to 188
      kx=1
      ky=iodd-1
      ink=inj
189   do 190 mw=kx,ky
      ink=ink+1
      lx=jmat(ink)
      if (lx.lt.0) go to 191
      lo=intp(lx)
      go to 192
191   lo=intn(-lx)
192   if (lo.eq.mw) go to 190
      ix=ini+lo
      iy=ini+mw
      ni=imat(iy)
      imat(ix)=ni
      imat(iy)=lx
      irm(ir)=-irm(ir)
      if (lx.lt.0) go to 193
      intp(lx)=mw
      go to 194
193   intn(-lx)=mw
194   if (ni.lt.0) go to 195
      intp(ni)=lo
      go to 190
195   intn(-ni)=lo
190   continue
      if (ky.gt.mu) go to 196
188   kx=iodd+1
      ky=m
      ink=inj+iodd
      go to 189
  196 vm=orb*f(ls)
      if (irm(ir).ne.jrm(kw)) vm=-vm
      ink=iodd+ini
      ia=imat(ink)
      ink=iodd+inj
      ib=jmat(ink)
      if (ia.gt.0) go to 197
      ia=-ia
      ib=-ib
197   ia = min(ia,ib)+iky(max(ia,ib))
      q(ia)=q(ia)+vm
184   continue
180   continue
114   continue
   72 continue
 108  if(na.eq.1) go to 121
 104  nb=na-1
      read(lscr) nhb,imax
      do 122 jw=1,nb
      read(lscr) jkan,f,g
      nx=0
      ig=-kml
      do 122 iw=1,imax
      do 203 kw=1,nl
      nx=nx+1
  203 ktest(kw)=jkan(nx)
      ig=ig+kml
      call setsto(imo,0,jcon)
      if(nod.eq.0) go to 124
      do 125 l=1,nod
      nt=ktest(l)
 125  jcon(nt)=1
      if (ndb.eq.0) go to 207
124   kk=nod
      do 126 l=1,ndb
      kk=kk+1
      nt=ktest(kk)
126   jcon(nt)=2
207   if (nmns.gt.0) go to 210
      jrm(1)=1
      if (nps.eq.0) go to 211
      do 212 kw=1,nps
212   jmat(kw)=ktest(kw)
      if (ndb.eq.0) go to 220
211   lx=nps
      kx=nps
      do 213 kw=1,ndb
      kx=kx+1
      lx=lx+1
      lk=ktest(kx)
      jmat(lx)=lk
      lx=lx+1
213   jmat(lx)=-lk
      go to 220
210   lx=0
      do 214 kw=1,kml
      jrm(kw)=1
      jq=(ihog(kw)-1)*nmns+1
      jz=jq+nqns
      jy=imap(jq)
      do 215 lw=1,nod
      lx=lx+1
      if (lw.ne.jy) go to 216
      if (jz.eq.jq) go to 217
      jq=jq+1
      jy=imap(jq)
217   jmat(lx)=-ktest(lw)
      go to 215
216   jmat(lx)=ktest(lw)
215   continue
214   lx=lx+ndbq
      if (ndb.eq.0) go to 220
      kx=nod2
      kk=nod
      do 218 kw=1,ndb
      kk=kk+1
      lk=ktest(kk)
      kx=kx+2
      mx=kx
      do 218 lw=1,kml
      mx=mx+m
      jmat(mx+1)=lk
218   jmat(mx+2) =-lk
220   nxj=0
      inis=-nd
      do 127 l=1,jmax
      inis=inis+nd
      inx=(inis-1)*m
      do 128 ll=1,nl
      nxj=nxj+1
 128  jtest(ll)=ikan(nxj)
      nix=0
      if(nod.eq.0) go to 129
      do 130 ll=1,nod
      jt=jtest(ll)
      if(jcon(jt).gt.0) go to 130
      if(nix.eq.1) go to 127
      nix=1
 130  continue
      if(ndb.eq.0) go to 131
 129  kk=nod
      do 132 ll=1,ndb
      kk=kk+1
      jt=jtest(kk)
      jb=jcon(jt)
      if(jb.eq.2) go to 132
      if(jb.eq.0.or.nix.eq.1) go to 127
      nix=1
 132  continue
 131  ks=ig
      inj=-m
      do 222 kw=1,kml
      ks=ks+1
      orb=g(ks)
      inj=inj+m
      do 223 lw=1,imo
      intp(lw)=0
223   intn(lw)=0
      ink=inj
      do 224 lw=1,m
      ink=ink+1
      lx=jmat(ink)
      if (lx.lt.0) go to 225
      intp(lx)=lw
      go to 224
225   intn(-lx)=lw
224   continue
      ini=inx
      iny=inis
      do 226 lw=1,nd
      ini=ini+m
      ink=ini
      iny=iny+1
      nix=0
      do 227 mw=1,m
      ink=ink+1
      lo=imat(ink)
      if (lo.lt.0) go to 228
      if (intp(lo).gt.0) go to 227
229   if (nix.eq.1) go to 226
      nix=1
      iodd=mw
      go to 227
228   if (intn(-lo).eq.0) go to 229
227   continue
      if (iodd.eq.1) go to 230
      kx=1
      ky=iodd-1
      ink=ini
233   do 231 mw=kx,ky
      ink=ink+1
      lx=imat(ink)
      if (lx.lt.0) go to 232
      lo=intp(lx)
      go to 234
232   lo=intn(-lx)
234   if(lo.eq.mw) go to 231
      ix=inj+lo
      iy=inj+mw
      ni=jmat(iy)
      jmat(ix)=ni
      jmat(iy)=lx
      jrm(kw)=-jrm(kw)
      if (lx.lt.0) go to 235
      intp(lx)=mw
      go to 236
235   intn(-lx)=mw
236   if (ni.lt.0) go to 237
      intp(ni)=lo
      go to 231
237   intn(-ni)=lo
231   continue
      if (ky.gt.mu) go to 238
230   kx=iodd+1
      ky=m
      ink=ini+iodd
      go to 233
238   vm=orb*fj(iny)
      if (irm(iny).ne.jrm(kw)) vm=-vm
      ink=iodd+ini
      ia=imat(ink)
      ink=iodd+inj
      ib=jmat(ink)
      if (ia.gt.0) go to 240
      ia=-ia
      ib=-ib
240   ia = min(ia,ib) + iky(max(ia,ib))
      q(ia)=q(ia)+vm
226   continue
222   continue
127   continue
122   continue
      go to 133
121   read(lscr)
 133  read(lscr)
 103  nx=0
      ini=-idel
      inis=-nd
      in2=-kml
      do 134 l=1,jmax
      in2=in2+kml
      inis=inis+nd
      ini=ini+idel
      call setsto(imo,0,jcon)
      if(nod.eq.0) go to 840
      do 141 ll=1,nod
      nx=nx+1
      nt=ikan(nx)
 141  jcon(nt)=1
      if(ndb.eq.0) go to 142
 840  do 143 ll=1,ndb
      nx=nx+1
      nt=ikan(nx)
 143  jcon(nt)=2
 142  nx=nx-nl
      orb=h(l)
      if(nod.eq.0) go to 135
      do 136 ll=1,nod
      nx=nx+1
      nt=ikan(nx)+1
      idx=iky(nt)
 136  q(idx)=q(idx)+orb
      if(ndb.eq.0) go to 137
 135  orb=orb+orb
      do 138 ll=1,ndb
      nx=nx+1
      nt=ikan(nx)+1
      idx=iky(nt)
 138  q(idx)=q(idx)+orb
 137  if(l.eq.1) go to 134
      l1=l-1
      inj=-idel-m
      injs=-nd
      nxj=0
      do 139 ll=1,l1
      inj=inj+idel
      injs=injs+nd
      do 243 kw=1,nl
      nxj=nxj+1
 243  ktest(kw)=ikan(nxj)
      nix=0
      if(nod.eq.0) go to 244
      do 245 kw=1,nod
      jt=ktest(kw)
      if(jcon(jt).gt.0) go to 245
      if(nix.eq.1) go to 139
      nix=1
 245  continue
      if(ndb.eq.0) go to 246
 244  kk=nod
      do 247 kw=1,ndb
      kk=kk+1
      jt=ktest(kk)
      jc=jcon(jt)
      if(jc.eq.2) go to 247
      if(jc.eq.0.or.nix.eq.1) go to 139
      nix=1
 247  continue
 246  ink=in2
      do 248 kw=1,kml
      ihq=ihog(kw)
      ink=ink+1
      orb=gj(ink)
      in3=(ihq-1)*m+ini
      in7=ihq+inis
      do 249 lw=1,imo
      intp(lw)=0
 249  intn(lw)=0
      in8=in3
      do 250 lw=1,m
      in8=in8+1
      lx=imat(in8)
      if(lx.lt.0) go to 251
      intp(lx)=lw
      go to 250
 251  intn(-lx)=lw
 250  continue
      in4=inj
      in6=injs
      do 252 lw=1,nd
      in4=in4+m
      in5=in4
      nix=0
      in6=in6+1
      do 253 mw=1,m
      in5=in5+1
      lo=imat(in5)
      if(lo.lt.0) go to 254
      if(intp(lo).gt.0) go to 253
 255  if(nix.eq.1) go to 252
      nix=1
      iodd=mw
      go to 253
 254  if(intn(-lo).eq.0) go to 255
 253  continue
      if(iodd.eq.1) go to 256
      kx=1
      ky=iodd-1
      in5=in4
 257  do 258 mw=kx,ky
      in5=in5+1
      lx=imat(in5)
      if(lx.lt.0) go to 259
      lo=intp(lx)
      go to 260
 259  lo=intn(-lx)
 260  if(lo.eq.mw) go to 258
      ix=in3+lo
      iy=in3+mw
      ni=imat(iy)
      imat(ix)=ni
      imat(iy)=lx
      irm(in7)=-irm(in7)
      if(lx.lt.0) go to 261
      intp(lx)=mw
      go to 262
 261  intn(-lx)=mw
 262  if(ni.lt.0) go to 263
      intp(ni)=lo
      go to 258
 263  intn(-ni)=lo
 258  continue
      if(ky.gt.mu) go to 264
 256  kx=iodd+1
      ky=m
      in5=in4+iodd
      go to 257
 264  vm=orb*fj(in6)
      if(irm(in6).ne.irm(in7)) vm=-vm
      ia=iodd+in4
      ia=imat(ia)
      ib=iodd+in3
      ib=imat(ib)
      if(ia.gt.0) go to 265
      ia=-ia
      ib=-ib
 265  ia=min(ia,ib)+iky(max(ia,ib))
      q(ia)=q(ia)+vm
 252  continue
 248  continue
 139  continue
 134  continue
 100  continue
43    continue
      return
      end
      subroutine putqno(q,etot,tbas,mpos,mtype)
      implicit real*8  (a-h,o-z), integer (i-n)
      character *(*) tbas
      character *8 text,type,title,com
      logical otran,otri
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
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     * ctran(mxorb3),otran,otri
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/blkorbs/value(maxorb),occ(maxorb+1),
     * nbasis,newbas,ncol,ivalue,iocc,ipad
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character *8 zdate,ztime,zaccno,zanam
      common/jinfo/zdate,ztime,zaccno,zanam
      common/junkc/com(19),title(10)
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
      dimension q(*),type(2)
      data text/' mrd-ci'/
      data type/'no-sabf','no-aos'/
      data m29/29/
c
      otran=.true.
      if(mpos.gt.350)call caserr(
     * 'invalid section specified for natural orbital output')
      nbsq=newbas*ncol
      j=1+lensec(mach(8))+lensec(mach(9))+lensec(nbsq)
      call secput(mpos,3,j,iblk)
      do 20 i=1,10
 20   title(i)=ztitle(i)
      write(iwr,1)text,tbas,mpos,ibl3d,yed(idaf)
 1    format(//a7,
     *' natural orbitals ( ',a4,' basis) stored in section',i4,
     *' of dumpfile starting at block',i6,
     *' of ',a4)
      occ(maxorb+1)=etot
      com(1)=zanam
      com(2)=zdate
      com(3)=ztime
      com(4)=text(2:7)
      com(5)=type(mtype)
      call wrtc(com,m29,iblk,idaf)
      call wrt3s(value,mach(8),idaf)
      nav = lenwrd()
      call wrt3is(ilifc(1),mach(9)*nav,idaf)
      call wrt3s(q,nbsq,idaf)
      call clredx
      call revind
      return
      end
      subroutine pmrd1(q,f,s,scr1,scr2,lj,ideks,iorbs,icore,knu,
     1imo,ksum,n,mfg)
      implicit real*8  (a-h,o-z), integer (i-n)
      dimension q(*),f(*),s(*),scr1(*),scr2(*),lj(*),ideks(*),mfg(*)
      imfg=0
      do 1 i=1,icore
      do 1 j=1,iorbs
      imfg=imfg+1
 1    scr1(imfg)=2.0d0*f(imfg)
      ig=0
      do 4 i=1,n
      ljn=lj(i)
      if=ig
      ig=ig+ljn
      iscr1=imfg-iorbs
      do 4 j=1,ljn
      lorie=if+j
      matt=ideks(lorie)
      do 4 kay=1,iorbs
      imfg=imfg+1
      tim=0.0d0
      kate=kay+iscr1
      do 5 lamb=1,ljn
      ianita=if+lamb
      if (j.lt.lamb) go to 6
      kent=matt+ianita
      go to 7
 6    kent=ideks(ianita)+lorie
 7    kate=kate+iorbs
 5    tim=tim+f(kate)*q(kent)
 4    scr1(imfg)=tim
      mike=0
      do 8 i=1,iorbs
      do 8 j=1,i
      tim=0.0d0
      mike=mike+1
      martyn=-iorbs
      do 9 kay=1,ksum
      martyn=martyn+iorbs
      iscr1=martyn+i
      kate=martyn+j
 9    tim=tim+scr1(iscr1)*f(kate)
 8    scr2(mike)=tim
      mike=0
      do 10 i=1,iorbs
      do 10 j=1,i
      mike=mike+1
      tom=s(mike)*scr2(mike)
      if (i.ne.j) tom=tom+tom
      scr1(mike)=tom
 10   continue
      knu2=ideks(knu+1)
      call vclr(scr2,1,knu2)
      ig=0
      do 30 i=1,iorbs
      kay=mfg(i)
      mike=ideks(kay)
      do 30 j=1,i
      imfg=mfg(j)
      if (kay.lt.imfg) go to 31
      kate=mike+imfg
      go to 32
  31  kate=ideks(imfg)+kay
  32  ig=ig+1
      scr2(kate)=scr2(kate)+scr1(ig)
  30  continue
      return
      end
      subroutine nmruc(chg,cx,cy,cz,knu)
      implicit real*8  (a-h,o-z), integer (i-n)
      common/scrtch/dinx,diny,dinz,qnxx,qnyy,qnzz,qnxy,qnxz,qnyz
      dimension chg(*),cx(*),cy(*),cz(*)
      dinx=0.0d0
      diny=0.0d0
      dinz=0.0d0
      qnxx=0.0d0
      qnyy=0.0d0
      qnzz=0.0d0
      qnxy=0.0d0
      qnxz=0.0d0
      qnyz=0.0d0
      do 51 i=1,knu
      za=chg(i)
      zx=cx(i)
      zy=cy(i)
      zz=cz(i)
      dinx=za*(zx) + dinx
      diny=za*(zy) + diny
      dinz=za*(zz) + dinz
      qnxx=qnxx+ za*(zx)*(zx)
      qnyy=qnyy+ za*(zy)*(zy)
      qnzz=qnzz+ za*(zz)*(zz)
      qnxy=qnxy+ za*(zx)*(zy)
      qnxz=qnxz+ za*(zx)*(zz)
      qnyz=qnyz+ za*(zy)*(zz)
   51 continue
      return
      end
      subroutine pmrdci(x,lword)
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/ftape/nf2(12),ifox,idum4(7)
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
      dimension lj(8),iq(15),x(*)
c
      write(iwr,3)yed(idaf),ibl3d
 3    format(/1x,104('=')//
     *40x,44('*')/
     *40x,'Table-Ci  --  one-electron properties module'/
     *40x,44('*')//
     *1x,'dumpfile on ',a4,' at block',i6)
      if(ibl3d)12,12,13
 12   call caserr('invalid starting block for dumpfile')
 13   write(iwr,18)lword
 18   format(/' main core available = ',i8,' words')
      if(num.le.0.or.num.gt.maxorb) call caserr(
     *'invalid number of basis functions')
      call rewftn(ifox)
      read(ifox)iorbs,imo,icore,nsym,lj
      if(iorbs.ne.num.or.imo.gt.iorbs)
     * call caserr('inconsistent parameters on dumpfile')
      nx=iorbs*(iorbs+1)/2
      lenmo=imo*(imo+1)/2
      lensym=icore
      do 4 i=1,nsym
 4    lensym=lensym+lj(i)*(lj(i)+1)/2
      lenact=(imo+icore)*iorbs
      iq(1)=1
      do 5000 i=2,15
5000  iq(i)=iq(i-1)+nx
      if(nx*15.gt.lword)call caserr('insufficient memory available')
      lwor= lword-nx*15-1
      write(iwr,19)lwor
 19   format(/' main core not used = ',i8,' words')
      call pmrd0(x(iq(1)),x(iq(2)),x(iq(3)),x(iq(4)),x(iq(5)),
     * x(iq(6)),x(iq(7)),x(iq(8)),x(iq(9)),x(iq(10)),x(iq(11)),
     * x(iq(3)),x(iq(4)),x(iq(6)),x(iq(7)),x(iq(10)),
     + x(iq(12)),x(iq(14)),
     * nx,lenact,lensym,lenmo,iwr)
      call clredx
      call rewftn(ifox)
c
      return
      end
      subroutine pmrd0(ovl,dipx,dipy,dipz,
     * qdxx,qdyy,qdzz,qdxy,qdxz,qdyz,tvl,tigl,
     * temp,pig,q,f,scr1,scr2,len,lenact,lensym,lenmo,iwr)
      implicit real*8  (a-h,o-z), integer (i-n)
      logical iaopr,imopr
      character *2 stype
      character *1 dash
      character *4 imos,iaos,mopr,is
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/scrtch/dinx,diny,dinz,qnxx,qnyy,qnzz,qnxy,qnxz,qnyz,
     *chg(100),cx(100),cy(100),cz(100),
     *cc(2550),cbb(2550),x(255),y(255),z(255),
     *ic(256),ii(256),ij(256),ik(256),ii4(256),ilife(256),
     *zeta(2550),d(2550),anorm(255),pvec(3),
     *h(11),s(11),t(11),lcomp(256),ipigg(10),lj(8),
     *jcon(maxorb),jkon(20)
      common/aaaa/icom(21)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/ftape/ktape,ntape,irdd(10),iput,idum4(7)
      common/craypk/istate,ipig(mxroot),ipmos,ipaos,
     * iaopr(11),imopr(11),jkonp(21*mxroot)
c
      dimension ovl(len),dipx(len),dipy(len),dipz(len)
      dimension a(7401)
      dimension qdxx(len),qdyy(len),qdzz(len),qdxy(len),
     *qdxz(len),qdyz(len),tigl(len),tvl(len),f(lenact),
     *temp(lenact),q(lenmo),pig(lensym),ch(400),
     +scr1(len),scr2(len)
      dimension mfg(256)
c
      dimension t1(5),t2(5),t3(5),zxf(12),zyf(12),zzf(12),alf(12),
     1 bf(5),af(7),dif(3),dif2(3)
      dimension newya(255),new(255),
     * stype(10),inumb(10),mopr(11),is(3)
      dimension ia(4),ja(4),ka(4)
      character*10 charwall
c
      equivalence (a(1),cc(1)),(ch(1),chg(1))
      data inumb/0,9,3,1,12,10,4,18,6,2/
      data stype/
     *'s','x','y','z','xy','xz','yz','xx','yy','zz'/
      data dash/'-'/
      data mopr/
     *'s','x','y','z','xx','yy','zz','xy','xz','yz','t'/
      data d5/1.0d-5/
      data is,iaos,imos/
     *'s','p','d','a.o.','m.o.'/
      ind(i,j)=max(new(i),new(j))*
     + (max(new(i),new(j))-1)/2 + min(new(i),new(j))
      cpu=cpulft(1)
      write(iwr,7002)cpu ,charwall()
 7002 format(/' *** commence property evaluation at ',
     *f8.2,' seconds',a10,' wall')
      nshl=20
      thrs=140.0d0
      fac=2.0d0**1.5d0
      call rewftn(iput)
      read (iput)iorbs,imo,icore,n,lj,knu,mlec,f,a,ch,newya
c
      call prpind(cx,cy,cz,ii4,newya,new,iorbs,knu,iwr)
c
      if(istate.le.0.or.istate.gt.mxroot)call caserr(
     *'invalid number of ci vectors requested for analysis')
      do 7000 i=1,istate
      if(ipig(i).le.0)call caserr(
     *'invalid ci vector specified for analysis')
 7000 continue
      write(iwr,9000)istate,
     * iorbs,imo,icore,mlec,ztitle,(ipig(i),i=1,istate)
 9000 format(/
     *' no. of states to be treated         ',i3/
     *' no. of a.o. basis functions         ',i3/
     *' no. of active orbitals              ',i3/
     *' no. of frozen orbitals              ',i3/
     *' no. of correlated electrons         ',i3//
     *' *** case : ',10a8//
     *' locations on ft42 of associated     '/
     *' ci density matrices                 ',20i3/)
c
c     ao integral print specification
c
      if(oprint(32))write(iwr,9100)
 9100 format(
     *' print out of basis functions requested'/)
      if(ipaos.ne.0)go to 9004
      write(iwr,9005)iaos
 9005 format(
     *' no printing of integrals in ',a4,' basis required')
      go to 9006
 9004 write(iwr,9007)iaos
 9007 format(' print following integrals in ',a4,' basis : ')
      do 9008 j=1,11
      if(.not.iaopr(j))
     *write(iwr,9010)mopr(j)
 9010 format(/10x,a2,' - matrix')
 9008 continue
c
c    mo integral print specification
c
 9006 if(ipmos.ne.0)go to 9011
      write(iwr,9005)imos
      go to 9012
 9011 write(iwr,9007)imos
      do 9013 j=1,11
      if(.not.imopr(j))
     *write(iwr,9010)mopr(j)
 9013 continue
c
 9012 write(iwr,22)
 22   format(/29x,'geometry (a.u.)'/29x,15('-')/
     *8x,'x',15x,'y',15x,'z',17x,'charge'/1x,63('='))
      do 9001 m=1,knu
 9001 write(iwr,9002)cx(m),cy(m),cz(m),chg(m)
 9002 format(4f16.7)
      write(iwr,9003)knu
 9003 format(/' no. of nuclei = ',i3)
      if(.not.oprint(32))go to 9029
      write(iwr,9020)iorbs
 9020 format(/1x,104('-')//40x,19('-')/40x,'molecular basis set'/
     *40x,19('-')//
     *' no. of gtos = ',i3/
     */7x,'ordered list of basis functions'/)
      do 9021 i=1,iorbs
      kk=new(i)
      x0=x(kk)
      y0=y(kk)
      z0=z(kk)
      do 9022 j=1,knu
      if( (dabs(x0-cx(j)).lt.d5) .
     * and. (dabs(y0-cy(j)).lt.d5).
     * and. (dabs(z0-cz(j)).lt.d5) ) go to 9023
 9022 continue
      call caserr(
     *'attempt to cite basis function on unknown centre')
 9023 icen=j
      j=ii(kk)*9+ij(kk)*3+ik(kk)
      do 9024 k=1,10
      if(j.eq.inumb(k))go to 9025
 9024 continue
      call caserr(
     *'basis function of unknown type detected')
 9025 if(k.eq.1)j=1
      if(k.gt.1.and.k.le.4)j=2
      if(k.ge.5)j=3
      l=ic(kk)
      write(iwr,9026)i,icen,is(j),stype(k)
      do 9028 j=1,l
      kkj = j + ilife(kk)
 9028 write(iwr,9002)cbb(kkj),cc(kkj)
 9021 continue
 9029 write(iwr,9009)(dash,j=1,129)
 9009 format(/1x,129a1)
 9026 format(/' orbital centre   type sub-type'/
     *i8,i7,4x,a1,7x,a2//
     *11x,'ctran',12x,'zeta')
      do 421 i=1,iorbs
      kk=i
      x0=x(kk)
      y0=y(kk)
      z0=z(kk)
      do 422 j=1,knu
      if( (dabs(x0-cx(j)).lt.d5) .
     * and. (dabs(y0-cy(j)).lt.d5).
     * and. (dabs(z0-cz(j)).lt.d5) ) go to 423
 422  continue
      call caserr('parameter error in property routines')
 423  icen=j
      mfg(i)=j
 421  continue
c     write (iwr,5243) (mfg(i),i=1,22)
c5243 format(5x,18i8)
      call nmruc(chg,cx,cy,cz,knu)
      debye=2.541587d0
      h(1)=0.0d0
      h(2)=dinx
      h(3)=diny
      h(4)=dinz
      h(5)=qnxx
      h(6)=qnyy
      h(7)=qnzz
      h(8)=qnxy
      h(9)=qnxz
      h(10)=qnyz
      h(11)=0.0d0
      call rewftn(ktape)
      write (ktape) f
      call rewftn(ktape)
      t1(1)=4.0d0
      t2(1)=2.0d0
      t3(1)=0.5d0
      j=1
      do 128 i=2,5
      j=j+2
      t1(i)=(t1(i-1)*4)/j
      t3(i)=t3(i-1)*j*0.5d0
  128 t2(i)=dsqrt(t1(i))
      mm=0
      do 101 k=1,iorbs
      lcomp(k)=ic(k)
      lco=lcomp(k)
      kkk = ilife(k)
      kkkk = kkk
      do 30 l=1,lco
      zeta(l+kkkk)=cc(l+kkk)
  30  d(l+kkkk)=cbb(l+kkk)
      lco=ic(k)
      mi2=ii(k)
      t1f=fac
      if(mi2.gt.0) t1f=t1f*t2(mi2)
      mj2=ij(k)
      if(mj2.gt.0) t1f=t1f*t2(mj2)
      mk2=ik(k)
      if(mk2.gt.0) t1f=t1f*t2(mk2)
      mmmax=mi2
      if(mmmax.lt.mj2) mmmax=mj2
      if(mmmax.lt.mk2) mmmax=mk2
      mi3=mi2+1
      mj3=mj2+1
      mk3=mk2+1
      fun=mi2+mj2+mk2
      farx=0.5d0*fun+0.75d0
       gun=fun*0.5d0
      mi4=mi2*(mi2-1)/2
      mj4=mj2*(mj2-1)/2
      mk4=mk2*(mk2-1)/2
      x0=x(k)
      y0=y(k)
      z0=z(k)
      sum=0.0d0
      do 550 l=1,lco
      a1=zeta(l+kkkk)
      ga=d(l+kkkk)*(a1**(0.75d0+gun))
      do 550 l2=1,lco
      b=zeta(l2+kkkk)
      gb=dsqrt(b)/(a1+b)
  550 sum=sum+ga*d(l2+kkkk)*(gb**(1.5d0+fun))
      sum=sum*fac*(2.0d0**fun)
      anorm(k)=1.0d0/dsqrt(sum)
      t1f=t1f*anorm(k)
      do 2 l=1,k
      mm=mm+1
      llll = ilife(l)
      ni3=ii(l)
      nj3=ij(l)
      nk3=ik(l)
      suma=0.0d0
      sumb=0.0d0
      sumc=0.0d0
      sumd=0.0d0
      sume=0.0d0
      sumf=0.0d0
      sumg=0.0d0
      sumh=0.0d0
      sumi=0.0d0
      sumj=0.0d0
      sumk=0.0d0
      t2f=t1f*anorm(l)
      if(ni3.gt.0) t2f=t2f*t2(ni3)
      if(nj3.gt.0) t2f=t2f*t2(nj3)
      if(nk3.gt.0) t2f=t2f*t2(nk3)
      fbrx=0.5d0*(ni3+nj3+nk3)+0.75d0
      nax=ni3
      if(nax.lt.nj3) nax=nj3
      if(nax.lt.nk3) nax=nk3
      nax=nax+2
      x5=x(l)
      y5=y(l)
      z5=z(l)
      ni4=ni3+1
      ni5=ni4+1
      nia=ni5+1
      nib=ni3-1
      if(ni3.lt.3) go to 143
      ni8=(ni3-3)*(ni3-2)/2-1
  143 ni6=ni3*(ni3-1)/2
      ni7=ni5*(ni5-1)/2-1
      nj4=nj3+1
      nj5=nj4+1
      nja=nj5+1
      njb=nj3-1
      if(nj3.lt.3) go to 144
      nj8=(nj3-3)*(nj3-2)/2-1
  144 nj6=nj3*(nj3-1)/2
      nj7=nj5*(nj5-1)/2-1
      nk4=nk3+1
      nk5=nk4+1
      nka=nk5+1
      nkb=nk3-1
      if(nk3.lt.3) go to 145
      nk8=(nk3-3)*(nk3-2)/2-1
  145 nk6=nk3*(nk3-1)/2
      nk7=nk5*(nk5-1)/2-1
      dif(1)=x(l)-x0
      dif(2)=y(l)-y0
      dif(3)=z(l)-z0
      difsum=0.0d0
      do 10 i1=1,3
      dif2(i1)=dif(i1)*dif(i1)
  10  difsum=difsum+dif2(i1)
      lc1=ic(l)
      nxx=mi2+ni3+1
      nxy=mj2+nj3+1
      nxz=mk2+nk3+1
      mxx=nxx+1
      mxy=nxy+1
      mxz=nxz+1
      nox=mxx
      if(nox.lt.mxy) nox=mxy
      if(nox.lt.mxz) nox=mxz
      bax=dif(1)
      bay=dif(2)
      baz=dif(3)
      zxf(1)=bax
      zyf(1)=bay
      zzf(1)=baz
      do 130 i=2,mxx
  130 zxf(i)=zxf(i-1)*bax
      do 132 i=2,mxy
  132 zyf(i)=zyf(i-1)*bay
      do 134 i=2,mxz
  134 zzf(i)=zzf(i-1)*baz
      do 3 l1=1,lco
      aa=zeta(l1+kkkk)
      qx=aa*x0
      qy=aa*y0
      qz=aa*z0
      tuffy=d(l1+kkkk)*(aa**farx)
      am=-aa
      af(1)=am
      aab2=difsum*aa
      do 136 i=2,nax
  136 af(i)=am*af(i-1)
      do 3 l2=1,lc1
      bb=zeta(l2+llll)
      bb3=-2.0d0*bb*bb
      alp=aa+bb
      alg=1.0d0/alp
      arg=aab2*bb*alg
      if(arg.gt.thrs) go to 3
      tuf=tuffy*dexp(-arg)
      tuf=tuf*d(l2+llll)*(bb**fbrx)
      if(mmmax.eq.0) go to 141
      bf(1)=bb
      if(mmmax.eq.1) go to 141
      do 907 i=2,mmmax
  907 bf(i)=bf(i-1)*bb
  141 tuf=tuf*(alg**1.5d0)
      alf(1)=alg
      do 151 i=2,nox
  151 alf(i)=alf(i-1)*alg
      sxma=0.0d0
      sxmb=0.0d0
      sxmc=0.0d0
      syma=0.0d0
      symb=0.0d0
      symc=0.0d0
      szma=0.0d0
      szmb=0.0d0
      szmc=0.0d0
      sxmd=0.0d0
      sxme=0.0d0
      symd=0.0d0
      syme=0.0d0
      szmd=0.0d0
      szme=0.0d0
      px=(qx+bb*x5)*alg
      py=(qy+bb*y5)*alg
      pz=(qz+bb*z5)*alg
      mi5=mi4
      do 140 i=1,mi3
      w4=1.0d0
      if(i.eq.1) go to 142
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  142 ni9=ni6
      do 146 j=1,ni4
      w5=w4
      if(j.eq.1) go to 147
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  147 i9=i+j-2
      if(i9.gt.0) w5=w5*zxf(i9)
      i8=nxx-i9
      i7=i8/2
      i4=i7+i9
      if(i8-i7*2.eq.0) go to 153
      w6=w5*alf(i4+1)
      if(i4.gt.0) w5=w5*alf(i4)
      w6=w6*t3(i7+1)
      if(i7.gt.0) w5=w5*t3(i7)
      sxma=sxma+w5
      sxmc=sxmc+w6
      go to 146
 153  w5=w5*alf(i4)*t3(i7)
      sxmb=sxmb+w5
 146   continue
       ni9=ni7
       do 300 j=1,nia
       ni9=ni9+1
       i9=i+j-2
       i8=mxx-i9
       i7=i8/2
       if(i8-i7*2.ne.0) go to 300
       w5=w4
       if(i7.eq.0) go to 301
       w5=w5*t3(i7)*alf(i7+i9)
       go to 302
 301   if(i9.gt.0) w5=w5*alf(i9)*zxf(i9)
       go to 303
 302    if(i9.gt.0) w5=w5*zxf(i9)
 303   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       sxmd=sxmd+w5
 300   continue
       if(ni3.lt.2) go to 140
       ni9=ni8
       do 304 j=1,nib
       ni9=ni9+1
       i9=i+j-2
       i8=mxx-i9-4
       i7=i8/2
       if(i8-i7*2.ne.0) go to 304
       w5=w4
       if(i7.eq.0) go to 305
       w5=w5*t3(i7)*alf(i7+i9)
       go to 306
305    if(i9.gt.0) w5=w5*alf(i9)*zxf(i9)
       go to 307
 306    if(i9.gt.0) w5=w5*zxf(i9)
 307   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       sxme=sxme+w5
 304   continue
 140   continue
      mi5=mj4
      do 160 i=1,mj3
      w4=1.0d0
      if(i.eq.1) go to 162
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  162 ni9=nj6
      do 166 j=1,nj4
      w5=w4
      if(j.eq.1) go to 167
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  167 i9=i+j-2
      if(i9.gt.0) w5=w5*zyf(i9)
      i8=nxy-i9
      i7=i8/2
      i4=i7+i9
      if(i8-i7*2.eq.0) go to 173
      w6=w5*alf(i4+1)
      if(i4.gt.0) w5=w5*alf(i4)
      w6=w6*t3(i7+1)
      if(i7.gt.0) w5=w5*t3(i7)
      syma=syma+w5
      symc=symc+w6
      go to 166
 173  w5=w5*alf(i4)*t3(i7)
      symb=symb+w5
 166  continue
       ni9=nj7
       do 320 j=1,nja
       ni9=ni9+1
       i9=i+j-2
       i8=mxy-i9
       i7=i8/2
       if(i8-i7*2.ne.0) go to 320
       w5=w4
       if(i7.eq.0) go to 321
       w5=w5*t3(i7)*alf(i7+i9)
       go to 322
 321   if(i9.gt.0) w5=w5*alf(i9)*zyf(i9)
       go to 323
 322    if(i9.gt.0) w5=w5*zyf(i9)
 323   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       symd=symd+w5
 320   continue
       if(nj3.lt.2) go to 160
       ni9=nj8
       do 324 j=1,njb
       ni9=ni9+1
       i9=i+j-2
       i8=mxy-i9-4
       i7=i8/2
       if(i8-i7*2.ne.0) go to 324
       w5=w4
       if(i7.eq.0) go to 325
       w5=w5*t3(i7)*alf(i7+i9)
       go to 326
325    if(i9.gt.0) w5=w5*alf(i9)*zyf(i9)
       go to 327
 326    if(i9.gt.0) w5=w5*zyf(i9)
 327   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       syme=syme+w5
 324   continue
 160   continue
      mi5=mk4
      do 180 i=1,mk3
      w4=1.0d0
      if(i.eq.1) go to 182
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  182 ni9=nk6
      do 186 j=1,nk4
      w5=w4
      if(j.eq.1) go to 187
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  187 i9=i+j-2
      if(i9.gt.0) w5=w5*zzf(i9)
      i8=nxz-i9
      i7=i8/2
      i4=i7+i9
      if(i8-i7*2.eq.0) go to 193
      w6=w5*alf(i4+1)
      if(i4.gt.0) w5=w5*alf(i4)
      w6=w6*t3(i7+1)
      if(i7.gt.0) w5=w5*t3(i7)
      szma=szma+w5
      szmc=szmc+w6
      go to 186
 193  w5=w5*alf(i4)*t3(i7)
      szmb=szmb+w5
 186   continue
       ni9=nk7
       do 340 j=1,nka
       ni9=ni9+1
       i9=i+j-2
       i8=mxz-i9
       i7=i8/2
       if(i8-i7*2.ne.0) go to 340
       w5=w4
       if(i7.eq.0) go to 341
       w5=w5*t3(i7)*alf(i7+i9)
       go to 342
 341   if(i9.gt.0) w5=w5*alf(i9)*zzf(i9)
       go to 343
 342    if(i9.gt.0) w5=w5*zzf(i9)
 343   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       szmd=szmd+w5
 340   continue
       if(nk3.lt.2) go to 180
       ni9=nk8
       do 344 j=1,nkb
       ni9=ni9+1
       i9=i+j-2
       i8=mxz-i9-4
       i7=i8/2
       if(i8-i7*2.ne.0) go to 344
       w5=w4
       if(i7.eq.0) go to 345
       w5=w5*t3(i7)*alf(i7+i9)
       go to 346
345    if(i9.gt.0) w5=w5*alf(i9)*zzf(i9)
       go to 347
 346    if(i9.gt.0) w5=w5*zzf(i9)
 347   if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
       szme=szme+w5
 344   continue
 180   continue
       sumg=sumg+tuf*sxma*syma*szma
       suma=suma+tuf*(sxmb+px*sxma)*syma*szma
       sumb=sumb+tuf*(symb+py*syma)*szma*sxma
       sumc=sumc+tuf*(szmb+pz*szma)*syma*sxma
       sumd=sumd+tuf*(sxmc+2.0d0*px*sxmb+px*px*sxma)*syma*szma
       sume=sume+tuf*(symc+2.0d0*py*symb+py*py*syma)*sxma*szma
       sumf=sumf+tuf*(szmc+2.0d0*pz*szmb+pz*pz*szma)*sxma*syma
       sumh=sumh+tuf*(sxmb+px*sxma)*(symb+py*syma)*szma
       sumi=sumi+tuf*(sxmb+px*sxma)*(szmb+pz*szma)*syma
       sumj=sumj+tuf*(symb+py*syma)*(szmb+pz*szma)*sxma
       w5=bb3*sxmd
       w5=w5+(2.0d0*ni3+1.0d0)*bb*sxma
       w5=w5-ni6*sxme
       w4=w5*syma*szma
       w5=bb3*symd
       w5=w5+(2.0d0*nj3+1.0d0)*bb*syma
       w5=w5-nj6*syme
       w4=w5*sxma*szma+w4
       w5=bb3*szmd
       w5=w5+(2.0d0*nk3+1.0d0)*bb*szma
       w5=w5-nk6*szme
       w4=w4+w5*sxma*syma
       sumk=sumk+tuf*w4
 3    continue
       ovl(mm)=t2f*sumg
       dipx(mm)=-t2f*suma
       dipy(mm)=-t2f*sumb
       dipz(mm)=-t2f*sumc
       qdxx(mm)=-t2f*sumd
       qdyy(mm)=-t2f*sume
       qdzz(mm)=-t2f*sumf
       qdxy(mm)=-t2f*sumh
       qdxz(mm)=-t2f*sumi
       qdyz(mm)=-t2f*sumj
       tvl(mm)=t2f*sumk
 2     continue
 101  continue
      call rewftn(ntape)
      write     (ntape) ovl
      write     (ntape) dipx
      write     (ntape) dipy
      write     (ntape) dipz
      write     (ntape) qdxx
      write     (ntape) qdyy
      write     (ntape) qdzz
      write     (ntape) qdxy
      write     (ntape) qdxz
      write     (ntape) qdyz
      write     (ntape) tvl
c
c    now print ao integrals if requested
c
      j=0
      do 9030 i=1,11
      if(iaopr(i))go to 9030
      write(iwr,9031)mopr(i)
 9031 format(//40x,'a.o. integrals ***   ',a2,'- matrix ***'/
     *40x,36('-')//)
      call writex(ovl(j+1),new,iorbs,iwr)
      write(iwr,9009)(dash,l=1,129)
 9030 j=j+len
c
      call rewftn(ntape)
      read (ktape) f
      call rewftn(ktape)
      ksum=icore+imo
      jmax=ksum*iorbs
      do 41 i=1,11
      read (ntape) tigl
      call vclr(temp,1,jmax)
      jdx=0
      do 42 j=1,iorbs
      do 42 k=1,j
      jdx=jdx+1
      big=tigl(ind(j,k))
      if (dabs(big).lt.1.0d-14) go to 42
      kdx=-iorbs
      do 43 l=1,ksum
      kdx=kdx+iorbs
      ig=kdx+j
      jg=kdx+k
      temp(ig)=temp(ig)+big*f(jg)
      if (j.eq.k) go to 43
      temp(jg)=temp(jg)+big*f(ig)
   43 continue
   42 continue
      kdx=0
      if (icore.eq.0) go to 44
      do 46 j=1,icore
      sum=0.0d0
      do 45 k=1,iorbs
      kdx=kdx+1
   45 sum=sum+f(kdx)*temp(kdx)
   46 pig(j)=sum
   44 imark=icore
      do 47 j=1,n
      ljn=lj(j)
      mdx=kdx
      do 48 k=1,ljn
      ldx=kdx
      do 50 l=1,k
      imark=imark+1
      sum=0.0d0
      mdy=mdx
      do 49 m=1,iorbs
      mdy=mdy+1
      ldx=ldx+1
   49 sum=sum+f(mdy)*temp(ldx)
   50 pig(imark)=sum
   48 mdx=mdy
   47 kdx=mdx
      if(imopr(i))go to 41
      write(iwr,9009)(dash,j=1,129)
      write (iwr,81)mopr(i)
 81   format(//40x,'m.o. integrals ***   ',a2,'- matrix ***'/
     *40x,36('-')//)
      if (icore.eq.0) go to 77
      write(iwr,9040)
 9040 format(/10x,'frozen core orbitals'/
     *'   orbital   integral'//)
      do 78 jj=1,icore
   78 write (iwr,79) jj,pig(jj)
   79 format (i10,f20.8)
   77 kdx=icore
      write(iwr,9041)
 9041 format(//10x,'integrals over active orbitals'//
     *4(' irrep.  i  j',19x)//)
      l=0
      do 80 jj=1,n
      ljn=lj(jj)
      do 80 j=1,ljn
      do 80 k=1,j
      kdx=kdx+1
      l=l+1
      ia(l)=jj
      ja(l)=j
      ka(l)=k
      temp(l)=pig(kdx)
      if(l.lt.4)go to 80
      write (iwr,9042) (ia(l),ja(l),ka(l),temp(l),l=1,4)
 9042 format(4(i4,3x,2i3,1x,f12.6,6x ))
      l=0
 80   continue
      if(l.ne.0)write(iwr,9042)
     *(ia(j),ja(j),ka(j),j=1,l)
   41 write (ktape) pig
      ich=0
      ijp=1
      do 53 i=1,istate
      write(iwr,9050)i
 9050 format(/1x,104('-')//
     *40x,37('-')/
     *40x,'molecular properties for state no.',i3/
     *40x,37('-')/)
      call setsto(imo,0,jcon)
      np=jkonp(ijp)
      do 7001 j=1,nshl
7001  jkon(j)=jkonp(ijp+j)
      ijp=ijp+21
      if(np.eq.0) go to 32
      do 34 j=1,np
      idx=jkon(j)
 34   jcon(idx)=1
 32   ib=(mlec-np)/2+np
      np1=np+1
      do 35 j=np1,ib
      idx=jkon(j)
 35   jcon(idx)=2
      write(iwr,1)
 1    format(/
     *' *** corresponding single configuration for this state :')
      if(np.eq.0)go to 9060
      write(iwr,9051)np,(jkon(j),j=1,np)
 9051 format(/' sequence nos. of ',i2,
     *' open shell orbitals'//
     *10x,6i3)
      go to 9052
 9060 write(iwr,9054)
 9054 format(/' no open shell orbitals')
 9052 j=ib-np
      write(iwr,9055)j,(jkon(j),j=np1,ib)
 9055 format(/
     *' sequence nos. of ',i2,' doubly occupied orbitals'//
     *10x,30i3)
   54 ich=ich+1
      if (ich.eq.ipig(i)) go to 55
      read (iput)
      go to 54
   55 read (iput) k,(q(j),j=1,k)
      call rewftn(ntape)
      read (ntape) tigl
      call rewftn(ntape)
      call pmrd1(q,f,tigl,scr1,scr2,lj,iky,iorbs,icore,
     +           knu,imo,ksum,n,mfg)
      idx=0
      do 61 j=1,imo
      do 61 k=1,j
      idx=idx+1
      if (j.eq.k) go to 61
      q(idx)=q(idx)+q(idx)
   61 continue
      call rewftn(ktape)
      do 70 j=1,11
      read (ktape) pig
      sum=0.0d0
      if (icore.eq.0) go to 72
      do 71 k=1,icore
   71 sum=sum+2.0d0*pig(k)
   72 scf=sum
      kdx=icore
      ipr=0
      do 73 k=1,n
      ljn=lj(k)
      lpr=ipr
      do 74 l=1,ljn
      lpr=lpr+1
      lrr=lpr
      mpr=ipr
      jrr=iky(lrr)
      do 74 m=1,l
      kdx=kdx+1
      mpr=mpr+1
      mrr=mpr+jrr
      ant=pig(kdx)
      sum=sum+ant*q(mrr)
      if (l.ne.m) go to 74
      if(jcon(lrr).eq.0) goto 74
      scf=scf+ant*jcon(lrr)
   74 continue
   73 ipr=lpr
      s(j)=scf
   70 t(j)=sum
      write (iwr,76) ich
   76 format(/30x,' molecular properties for state no.',i3,' on output t
     1ape'//)
      write            (iwr,201)
  201 format(25x,'el.scf',7x,'nuclear',10x,'el.ci',6x,'total scf',7x,'to
     1tal ci'/10x,85('=')/)
      ovs=s(1)
      ovc=t(1)
      zero=0.0d0
      write            (iwr,202) ovs,zero,ovc,ovs,ovc
  202 format(10x,'overlap',3x,5f15.8/)
      dinx=h(2)
      disx=s(2)
      dicx=t(2)
      f1=dinx+disx
      f2=dinx+dicx
      write            (iwr,203)disx,dinx,dicx,f1,f2
  203 format(10x,'dipole(x)',1x,5f15.8/)
      diny=h(3)
      disy=s(3)
      dicy=t(3)
      f1=diny+disy
      f2=diny+dicy
      write            (iwr,204)disy,diny,dicy,f1,f2
  204 format(10x,'dipole(y)',1x,5f15.8/)
      dinz=h(4)
      disz=s(4)
      dicz=t(4)
      f1=dinz+disz
      f2=dinz+dicz
      write            (iwr,205)disz,dinz,dicz,f1,f2
  205 format(10x,'dipole(z)',1x,5f15.8/)
      dinxd=dinx*debye
      disxd=disx*debye
      dicxd=dicx*debye
      dinyd=diny*debye
      disyd=disy*debye
      dicyd=dicy*debye
      dinzd=dinz*debye
      diszd=disz*debye
      diczd=dicz*debye
      write(iwr,2200)
 2200 format(/10x,'dipole moment in debye'/)
      f1=dinxd+disxd
      f2=dinxd+dicxd
      write(iwr,203)  disxd,dinxd,dicxd,f1,f2
      f1=dinyd+disyd
      f2=dinyd+dicyd
      write(iwr,204) disyd,dinyd,dicyd,f1,f2
      f1= dinzd+diszd
      f2=dinzd+diczd
      write(iwr,205) diszd,dinzd,diczd,f1,f2
      qnxx=h(5)
      qusxx=s(5)
      qcxx=t(5)
      f1=qnxx+qusxx
      f2=qnxx+qcxx
      write            (iwr,206) qusxx,qnxx,qcxx,f1,f2
  206 format(//10x,'quad(xx)',2x,5f15.8/)
      qnyy=h(6)
      qusyy=s(6)
      qcyy=t(6)
      f1=qnyy+qusyy
      f2=qnyy+qcyy
      write            (iwr,207) qusyy,qnyy,qcyy,f1,f2
  207 format(10x,'quad(yy)',2x,5f15.8/)
      qnzz=h(7)
      quszz=s(7)
      qczz=t(7)
      f1=qnzz+quszz
      f2=qnzz+qczz
      write            (iwr,208) quszz,qnzz,qczz,f1,f2
  208 format(10x,'quad(zz)',2x,5f15.8/)
      qnxy=h(8)
      qusxy=s(8)
      qcxy=t(8)
      f1=qnxy+qusxy
      f2=qnxy+qcxy
      write            (iwr,209) qusxy,qnxy,qcxy,f1,f2
  209 format(10x,'quad(xy)',2x,5f15.8/)
      qnxz=h(9)
      qusxz=s(9)
      qcxz=t(9)
      f1=qnxz+qusxz
      f2=qnxz+qcxz
      write            (iwr,210) qusxz,qnxz,qcxz,f1,f2
  210 format(10x,'quad(xz)',2x,5f15.8/)
      qnyz=h(10)
      qusyz=s(10)
      qcyz=t(10)
      f1=qnyz+qusyz
      f2=qnyz+qcyz
      write            (iwr,211) qusyz,qnyz,qcyz,f1,f2
  211 format(10x,'quad(yz)',2x,5f15.8/)
      tvs=s(11)
      tvc=t(11)
      write            (iwr,222) tvs,zero,tvc,tvs,tvc
  222 format(10x,'kinetic e.',5f15.8/)
   53 continue
      write(iwr,9009)(dash,j=1,129)
      call rewftn(ktape)
      call rewftn(ntape)
      cpu=cpulft(1)
      write(iwr,8001)cpu ,charwall()
 8001 format(//
     *' *** end of Table-Ci properties calculation at ',
     *f8.2,' seconds',a10,' wall'/)
      return
      end
      subroutine prpind(cx,cy,cz,iord,newya,new,iorbs,knu,jtape)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      dimension newya(*),new(*),cx(*),cy(*),cz(*),iord(*)
c
      data d5 /1.0d-5/
c
      nc=0
      nhfunc = 0
      do 555 loop=1,knu
      if(iord(loop).gt.nhfunc) nhfunc=iord(loop)
555   continue
c
      do 442 loop=1,nhfunc
      do 552 iat=1,non
      if(iord(iat).ne.loop)go to 552
      x0=c(1,iat)
      y0=c(2,iat)
      z0=c(3,iat)
      do 9027 j=1,non
      if( (dabs(x0-cx(j)).lt.d5) .
     * and. (dabs(y0-cy(j)).lt.d5).
     * and. (dabs(z0-cz(j)).lt.d5) ) go to 8023
 9027 continue
      call caserr(
     *'prpind: attempt to site basis function on unknown centre')
 8023 continue
      do 553 iii=1,nshell
      i=katom(iii)
      if(i.ne.j)go to 553
      mini=kmin(iii)
      maxi=kmax(iii)
      kk=kloc(iii)-mini
c
      do 25 iorb=mini,maxi
      li=kk+iorb
      nc=nc+1
 25   newya(nc)=li
c
 553  continue
 552  continue
 442  continue
      do 8998 i=1,iorbs
8998  new(newya(i)) = i
      if(oprint(32)) then
       do 8997 i=1,iorbs
 8997  write(jtape,*) 'prpind:  i,newya,new' , i, newya(i), new(i)
      endif
      return
      end
      subroutine moment(q,lword)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      common/scrtch/nconf(276)
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
      common/ftape/ispaci(13),ifox,idum4(6)
c
      dimension q(*)
c
c ifox=36 default
c iwod=1000 (block size)
c vnuc is nuclear repulsion
c zero is total scf energy
c imo is no. of active mos
c m is no. of active electrons
c ksum is active+frozen mos
c iorbs is total mos (active+frozen+discarded)
c nconf(5) are number of selected configs (not safs) in last extr pass
c 
      write(iwr,3)yed(idaf),ibl3d
 3    format(/1x,104('=')//
     *40x,38('*')/
     *40x,'Table-Ci  --  transition moment module'/
     *40x,38('*')//
     *1x,'dumpfile on ',a4,' at block',i6 )
      if(ibl3d)12,12,13
 12   call caserr('invalid starting block for dumpfile')
 13   write(iwr,18)lword
 18   format(/' main core available = ',i8,' words')
      if(num.le.0.or.num.gt.maxorb) call caserr(
     *'invalid number of basis functions')
      call rewftn(ifox)
      read(ifox)iwod,vnuc,zero,imo,m,nconf,ksum,iorbs
      if(iorbs.ne.num.or.imo.gt.ksum.or.ksum.gt.iorbs)
     * call caserr('invalid parameters on dumpfile')
      call rewftn(ifox)
      l22500=22500
      lentot=ksum*iorbs
      lensq=imo*imo
      nx=iorbs*(iorbs+1)/2
      i1=1
      i2=i1+lentot
      i3=i2+nx
      i4=i3+lentot
      i5=i4+l22500
      if(i5.gt.lword)call caserr('insufficient memory available')
      lwor= lword-i5-1
      write(iwr,8001)lwor
 8001 format(/' *** main core not used = ',i8,' words ')
      call tm(q(1),q(i2),q(i3),q(1),q(i4),
     * l22500,lentot,nx,lensq)
      call clredx
      call rewftn(ifox)
c
      return
      end
      subroutine tm(y,pig,dbq,q,df,lenf,lentot,lenbas,lensq)
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
      character *1 dash,is
      character *2 stype
      logical onew
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/aaaa/icom(21)
c
c 256 = maxorb, here is sometimes also 255
c 10  = max primitives in one contracted gaussian !!
c ispa4 was introduced to map the mrdci order to gamess scf 
c order, set in the symmetry adaption .
c
      common/scrtch/bc(2550),cbb(2550),x(255),p(255),r(255),
     *ic(256),ii(256),ij(256),ik(256),ispa4(256),ilife(256),
     *jmatsp(2855),
     *ovl(500),dixf(500),diyf(500),dizf(500),xnab(500),ynab(500),
     *znab(500),spc(19000),
     *icf(22500),ibal(8),itil(8),mcomp(256),kj(8),lsym(2040),
     *nj(8),ntil(8),nbal(9),ncomp(256),nzer(256),xdum(1300),
     *nconf(5),jconf(5),nytl(5),ndub(5),jtest(mxnshl),
     *nplu(5),itest(mxnshl),jytl(5),jdub(5),jplu(5),lj(8),imap(504),
     *ihog(48),jkan(mxcrc2),jmap(504),jhog(48),ikan(mxcrc2),
     +jcon(maxorb),imat(3780),irm(126),jrm(mxcrec),
     +intp(maxorb),intn(maxorb),
     +zeta(2550),d(2550),xorb(255),yorb(255),
     *zzorb(255),anorm(255),psep(6),pvec(3),dif(3),dif2(3),
     *wc(100),wx(100),wy(100),wz(100),t1(5),t2(5),t3(5),
     *zxf(11),zyf(11),zzf(11),alf(11),bf(5),af(6),
     *f(mxcrec),g(mxcrc2),fj(mxcrec),gj(mxcrc2),ci(126),cj(126),
     +di(48),dj(48)
c
      common/ftape/iaos,imos,iscr,idum4(17),itmfil
c iaos are one el integ in ao basis written in tm2
c imos are used only in this routine
c itmfil is the unit to permit subsequent TM analysis
      common/craypk/iput,idx,jput,jdx,nstate
      common/lsort/fac,acore,enee,ilifq(maxorb),iorbs,iswh,m,m2,mu,id
     *,ni8,nj8,nk8,imo,ig
      dimension y(lentot),dbq(lentot),pig(lenbas),q(lensq)
      dimension cf(22500),df(lenf),bloc(9100),ovm(1300),
     * xm(1300),ym(1300),zm(1300),xam(1300),yam(1300),zam(1300),
     * e(7401),
     * jmat(30000),dihf(3500),newya(255)
c
      character*10 charwall
      dimension new(maxorb),stype(10),inumb(10),is(3)
      dimension icen(maxorb)
      equivalence (e(1),bc(1),bloc(1),jmat(1))
      equivalence (cf(1),ovl(1)),(dihf(1),ovl(1)),
     * (bloc(1),ovm(1)),
     * (bloc(1301),xm(1)),(bloc(2601),ym(1)),(bloc(3901),
     *zm(1)),(bloc(5201),xam(1)),(bloc(6501),yam(1)),
     * (bloc(7801),zam(1))
c
      data d5,dash/1.0d-5,'-'/
      data inumb/0,9,3,1,12,10,4,18,6,2/
      data stype/
     *'s','x','y','z','xy','xz','yz','xx','yy','zz'/
      data is/'s','p','d'/
c
      ind(i,j)=max(new(i),new(j))*
     + (max(new(i),new(j))-1)/2 + min(new(i),new(j))
c
      cpu=cpulft(1)
      write(iwr,9000)cpu ,charwall()
 9000 format(/' *** commence moment evaluation at ',f8.2,
     * ' seconds',a10,' wall')
      istont=0
      if(iput.le.0.or.jput.le.0)call caserr(
     *'invalid fortran stream specified for ci vector')
      call rewftn(iscr)
      call rewftn(iaos)
      call rewftn(imos)
c
      irbmax=255
      mmax=30
      ksmax=255
c jd= block length on imos file
      jd=1300
      id=500
      fack=8066.0d0
      ev=27.211d0
      amult=2.02592d-6
      t1(1)=4.0d0
      t2(1)=2.0d0
      t3(1)=0.5d0
      j=1
      do 128 i=2,5
      j=j+2
      t1(i)=(t1(i-1)*4)/j
      t3(i)=t3(i-1)*j*0.5d0
  128 t2(i)=dsqrt(t1(i))
      fac=2.0d0**1.5d0
      read(iput)iwod,vnuc,zero,imo,m,nconf,newya,
     *nytl,nplu,ndub,iswh,ksum,
     1iorbs,knu,cf,icf,wc,wx,wy,wz,ibal,itil,mcomp,kj,lj,n,ifrk,e,lsym,
     2 nsel,nj,ntil,nbal,ncomp
      iflw=imo*imo
      icore=ksum-imo
      write(iwr,36)iput,idx,ztitle,jput,jdx,nstate,
     * iorbs,imo,icore,m
 36   format(/
     *' location of ground state ci vector on ft',i2,7x,i3/
     *' case : ',10a8//
     *' excited state vector(s) reside on ft',i2/
     *' location of first excited state vector           ',i3/
     *' number of excited state vectors to be treated    ',i3/
     *' number of contracted a.o. basis functions        ',i3/
     *' number of active orbitals                        ',i3/
     *' number of frozen orbitals                        ',i3/
     *' number of correlated electrons                   ',i3)
      if(oprint(32))write(iwr,9004)
 9004 format(/
     *' print out of basis orbitals requested')
      do 900 i=1,iorbs
      new(i)=0
900   continue
c     do 31 i=1,iorbs
c31   new(newya(i))=i
      call prpind(wx,wy,wz,ispa4,newya,new,iorbs,knu,iwr)
      onew=.false.
      do 910 i=1,iorbs
      if(new(i).ne.i) onew=.true.
910   continue
      if(onew) then
       if(oprint(32)) then
       write(iwr,*)' tm: array "new" is NOT identical permutation'
       write(iwr,32)(new(i),i=1,iorbs)
32     format(20i5)
       endif
      else
      if(oprint(32))write(iwr,*)
     +              ' tm: array "new" is identical permutation'
      endif
      if(oprint(32)) write(iwr,9997)
 9997 format(/29x,'geometry (a.u.)'//
     *8x,'x',15x,'y',15x,'z',17x,'charge'/)
      if(oprint(32)) then
       do 9001 j=1,knu
 9001  write(iwr,9002)wx(j),wy(j),wz(j),wc(j)
      endif
 9002 format(4f16.7)
      if(oprint(32)) write(iwr,9003)knu
 9003 format(/' no. of nuclei = ',i3)
      if(oprint(32)) write(iwr,9020)iorbs
 9020 format(/' *** basis function specification'//
     *' no. of gtos = ',i3//
     *7x,'ordered list of basis functions'/)
      do 9021 i=1,iorbs
      kk=new(i)
      x0=x(kk)
      y0=p(kk)
      z0=r(kk)
      do 9022 j=1,knu
      if( (dabs(x0-wx(j)).lt.d5) .
     * and. (dabs(y0-wy(j)).lt.d5).
     * and. (dabs(z0-wz(j)).lt.d5) ) go to 9023
 9022 continue
      call caserr(
     *'attempt to cite basis function on unknown centre')
 9023 icen(kk)=j
      j=ii(kk)*9+ij(kk)*3+ik(kk)
      do 9024 k=1,10
      if(j.eq.inumb(k))go to 9025
 9024 continue
      call caserr(
     *'basis function of unknown type detected')
 9025 if(k.eq.1)j=1
      if(k.gt.1.and.k.le.4)j=2
      if(k.ge.5)j=3
      l=ic(kk)
      if(oprint(32)) write(iwr,9026)i,icen(kk),is(j),stype(k)
      if(oprint(32)) then
      kkkk = ilife(kk)
      do 9028 j=1,l
 9028 write(iwr,9002)cbb(kkkk+j),bc(kkkk+j)
      endif
 9021 continue
c
c save important information for tmatrix conversion 
c to cartesian representation
c
      if(itmfil.ne.0) then
      call rewftn(itmfil)
      write(itmfil) nstate,iorbs,imo,icore,m,knu
c  geometry
      write(itmfil) (wc(j),j=1,knu)
      write(itmfil) (wx(j),j=1,knu)
      write(itmfil) (wy(j),j=1,knu)
      write(itmfil) (wz(j),j=1,knu)
c
c  basis
c     new is permutation of 1..iorbs
c  we will better permute the molcao coefficients with newya
c
c     write(itmfil) (new(j),j=1,iorbs)
      write(itmfil) (j,j=1,iorbs)
c centers
      write(itmfil) (icen(j),j=1,iorbs)
c x,y,z exponents
      write(itmfil) (ii(j),j=1,iorbs)
      write(itmfil) (ij(j),j=1,iorbs)
      write(itmfil) (ik(j),j=1,iorbs)
c number of primitives
      write(itmfil) (ic(j),j=1,iorbs)
c contraction coefficients and exponents
      do 8000, j=1,iorbs
      kkkk = ilife(j)
      do 8001, k=1,ic(j)
      write(itmfil) cbb(kkkk+k),bc(kkkk+k)
8001  continue
8000  continue
c
      endif
c
      write(iwr,9009)(dash,j=1,129)
 9009 format(/1x,129a1)
 9026 format(/
     *' orbital centre   type sub-type'/
     *i8,i7,4x,a1,7x,a2//
     *11x,'ctran',12x,'zeta')
      jy=0
      do 450 i=1,iorbs
      nzer(i)=jy
450     jy=jy+ncomp(i)
       jy=0
      ma=0
       ibuk=0
        jcl=icore*iorbs
        do 451 i=1,n
        k2=kj(i)
        k4=k2-lj(i)
         kx=ntil(i)
                  do 451 j=1,k2
        if(j.gt.k4) go to 452
        j8=ibuk
        k5=j8+1
        ibuk=ibuk+iorbs
        jt=ibuk
        go to 453
 452    j8=jcl
        k5=jcl+1
         jcl=jcl+iorbs
         jt=jcl
453    do 454 k=k5,jt
454    y(k)=0.0d0
       ma=ma+1
       ll=mcomp(ma)
       do 455 k=1,ll
       jy=jy+1
        vm=cf(jy)
            nxx=icf(jy)+kx
        lx=nzer(nxx)
        nxx=ncomp(nxx)
        do 455 l=1,nxx
        lx=lx+1
        ifj=lsym(lx)
        if(ifj.gt.0) go to 456
        y(j8-ifj)=-vm
        go to 455
456     y(j8+ifj)=vm
455    continue
451    continue
c
c  save y(ao,mo) coeffs - coeffs of transformed mos to aos
c  core and discarded orbitals are skipped and the rest is sorted 
c  according to irreps
c  in AO the order is like in the scf printout , or like here
c  in increasing order with indirection through the array newya()
c
c the permutation new is usually identity, but does not have to be
c
c  arrays cf and icf contain only those coefficients which are finite
c  due to the symmetry and they are kept only with one sign
c  you can have a look at the transformation above
c  the coefficients refer to normalised contracted gaussians
c
c when omitting the permutation new, molcao coefs are in the 
c following order:
c mo index: mrdci numbering(as in configuration specification)
c ao index: scf numbering permuted through newya
c
c the tmat plotting will be shielded from it by permuting molcao 
c coefs when saving them to tmfile 
c (by inverse permutation the the one discussed above)
c
      if(itmfil.ne.0) then
      write(iwr,*)' molcao coefficients written to tmfile',itmfil
      if(oprint(32)) then 
       write(iwr,*)' 1st index is active mo no. in mrdci convention'
       write(iwr,*)' 2nd index is ao in the usual scf ordering'
      endif
      ibase=icore*iorbs
      do 8004,i=1,imo
      write(itmfil)(y(ibase+newya(j)),j=1,iorbs)
      if(oprint(32)) then
         do 9091,j=1,iorbs
9091     write(iwr,*)'i,j,newya(j),y(newya(j),i) ',
     +    i,j,newya(j),y(ibase+newya(j))
      endif
      ibase=ibase+iorbs
8004  continue
      endif
c
      zbx=0.0d0
      zby=0.0d0
      zbz=0.0d0
      do 901 i=1,knu
      c=wc(i)
      zbx=zbx+c*wx(i)
      zby=zby+c*wy(i)
  901 zbz=zbz+c*wz(i)
      write (iwr,902) zbx,zby,zbz
  902 format(/'  nuclear dipole components : ',3f15.8)
      if (iorbs.gt.irbmax) go to 750
      if (m.gt.mmax) go to 750
      if (ksum.gt.ksmax) go to 750
      j=0
c set column pointers, one column after another
      do 148 i=1,imo
      ilifq(i)=j
 148  j=j+imo
      if (iput.eq.jput) go to 40
c
      call rewftn(jput)
      read (jput) iwod,vnuc,zero,jmo,m2,jconf,newya,jytl,
     * jplu,jdub,iswh,jsum,
     1jorbs,knu,df
      jcore=jsum-jmo
      if (jmo.ne.imo) go to 751
      if (jcore.ne.icore) go to 752
      if (jorbs.ne.iorbs) go to 753
      if (jsum.ne.ksum) go to 754
      if (m2.ne.m) go to 755
      if (istont.ne.0) go to 40
      do 35 i=1,jy
      if (dabs(cf(i)-df(i)).gt.1.0d-7) go to 757
35    continue
40    call tm2
c necessary integrals calculated
      fcorx=0.0d0
      fcory=0.0d0
      fcorz=0.0d0
       iaa=iorbs*(iorbs+1)/2
       iorq=imo*(imo+1)/2
       iblk=(iaa-1)/id
       ires=iaa-iblk*id
       iblk=iblk+1
       iors=ksum*iorbs
       ixy=icore*iorbs
       call rewftn(imos)
       iy=-id
       do 812 ma=1,7
       call rewftn(iaos)
       iy=iy+id
       ig=0
       iend=id
       do 813 i=1,iblk
       read(iaos) dihf
       if(i.eq.iblk) iend=ires
       jy=iy
       do 813 j=1,iend
       ig=ig+1
       jy=jy+1
 813   pig(ig)=dihf(jy)
      call vclr(dbq,1,iors)
       do 815 i=1,iorbs
       do 815 j=1,i
       jab=ind(i,j)
       big=pig(jab)
       if(dabs(big).lt.1.0d-14) go to 815
       kab=-iorbs
       do 816 k=1,ksum
       kab=kab+iorbs
       ig=kab+i
       jg=kab+j
       dbq(jg)=dbq(jg)+big *y(ig)
       if(i.eq.j) go to 816
       if(ma.gt.4) go to 216
       dbq(ig)=dbq(ig)+big*y(jg)
        go to 816
216    dbq(ig)=dbq(ig)-big*y(jg)
816   continue
815   continue
       kab=0
       if(icore.eq.0) go to 817
c no core orbitals
       if(ma.eq.1.or.ma.gt.4) go to 822
       if(ma-3) 818,819,820
818    do 821 i=1,icore
      fcorx=fcorx+ddot(iorbs,y(kab+1),1,dbq(kab+1),1)
 821  kab=kab+iorbs
       go to 817
822    kab=kab+ixy
       go to 817
819    do 823 i=1,icore
      fcory=fcory+ddot(iorbs,y(kab+1),1,dbq(kab+1),1)
 823  kab=kab+iorbs
       go to 817
820    do 824 i=1,icore
      fcorz=fcorz+ddot(iorbs,y(kab+1),1,dbq(kab+1),1)
 824  kab=kab+iorbs
817    jg=0
       mdx=kab
       do 825 i=1,imo
       ldx=kab
       do 826 j=1,i
       mdy=mdx
       jg=jg+1
       sum=0.0d0
       do 827 k=1,iorbs
       mdy=mdy+1
       ldx=ldx+1
 827   sum=sum+y(mdy)*dbq(ldx)
       xdum(jg)=sum
c here is the multiplication with molcao coefs
       if(jg.lt.jd) go to 826
       write(imos) xdum
       jg=0
 826   continue
 825   mdx=mdy
       if(jg.gt.0) write(imos) xdum
 812   continue
       iblk=(iorq-1)/jd
       ires=iorq-iblk*jd
       iblk=iblk+1
       fcorx=fcorx+fcorx
       fcory=fcory+fcory
       fcorz=fcorz+fcorz
       write(iwr,603) fcorx,fcory,fcorz
603    format(/5x,'core dipole components : ',3f15.8)
c
c subtract the nuclear component (which should be zero because of 
c standard orientation with [0,0,0] in the charge centre
c but even if the molecule is shifted, it must be multiplied with the
c overlap, which is zero
c
       fcorx=fcorx-zbx
       fcory=fcory-zby
       fcorz=fcorz-zbz
      if (idx.eq.1) go to 50
      idy=idx-1
      do 45 i=1,idy
      do 51 j=1,iswh
      nc=nconf(j)
      if (nc.eq.0) go to 51
      read(iput) nhb
      do 52 k=1,nhb
52    read(iput)
51    continue
45    continue
   50 do 57 j=1,iswh
      if (nconf(j).eq.0) go to 57
      read(iput)nhb,imax,ndt,kml,imap,ihog,eneg
      write(iscr)nhb,imax,ndt,kml,imap,ihog
      do 58 k=1,nhb
      read(iput) jkan,f,g
58    write(iscr) jkan,f,g
   57 continue
      if (iput.ne.jput) go to 60
      do 66 i=1,iswh
      jconf(i)=nconf(i)
      jytl(i)=nytl(i)
      jdub(i)=ndub(i)
66    jplu(i)=nplu(i)
 60   if(iput.eq.jput)jdx=jdx-idx
      if (jdx.eq.1) go to 64
      jdy=jdx-1
      do 65 i=1,jdy
      do 69 j=1,iswh
      nc=jconf(j)
      if (nc.eq.0) go to 69
      read(jput) nhb
      do 70 k=1,nhb
70    read(jput)
69    continue
65    continue
c
64    continue
c
c  here is the main loop for all states treated
c
      do 200 i=1,nstate
      write(iwr,9030)i
 9030 format(/1x,104('-')//40x,43('-')/
     *40x,'moment calculation for excited state no.',i3/
     *40x,43('-')/)
      acore=0.0d0
      call vclr(q,1,iflw)
      call tm3(q,lensq)
      write(iwr,499) acore
499   format(/2x,'overlap between upper and lower state is',f20.16)
      if(oprint(32)) then
      write(iwr,299)
 299  format(/
     *40x,'transition density matrix'/40x,25('-'))
c          in mrd-ci frozen, discarded and reordered mo basis
      call writem(q,ilifq,imo,imo)
      endif
c
c save tmat for this state
c
      if(itmfil.ne.0) then
c write in columns
      do 8003,j=1,imo
      jpoint=ilifq(j)
      write(itmfil) (q(k+jpoint),k=1,imo)
8003  continue
      endif
c
      del=enee-eneg
      delv=del*ev
      delc=delv*fack
      write(iwr,501) eneg,enee,del,delv,delc,45.564d0/del
 501  format(//' transition energy data'/1x,22('-')//
     *' total ci energy of ground state   ',f16.6,' hartree'/
     *' total ci energy of excited state  ',f16.6,' hartree'//
     *' excitation energy                 ',f16.6,' hartree'/
     *' excitation energy                 ',f16.6,' e.v.'/
     *' excitation energy                 ',f16.6,' cm-1'/
     *' excitation wavelength             ',f16.6,' nm'/)
      x1g=acore*fcorx
      x2g=acore*fcory
      x3g=acore*fcorz
c the core and nuclear contributions are multiplied by the overlap
c nuclear was added to the core somewhere above
      y1g=0.0d0
      y2g=0.0d0
      icam=0
      iend=jd
      y3g=0.0d0
      ig=jd
      write(iwr,9009)(dash,j=1,60)
      if(oprint(32))then
        finthr=0.01d0
      else
        finthr=0.1d0
      endif
      write(iwr,9031) finthr
 9031 format(//' finite contributions to oscillator strengths '/
     + ' (threshold is |d(k|l)| >= ',f6.3,' )'/1x,80('-')///
     *6x,'i',4x,'j',8x,'d(i,j)',14x,'x',14x,'y',14x,'z',
     *9x,'del<x>',9x,'del<y>',9x,'del<z>',14x,'s'/5x,127('-')//)
c
c use tmat and integrals to calculate contributions
c loops only for 1 triangle, but checks both 
c jk and kj element if k.ne.j
c
      do 502 j=1,imo
      iaf=ilifq(j)
      do 502 k=1,j
      iag=iaf+k
c row k column j -th element
      cp=q(iag)
      if(ig.lt.jd) go to 830
      icam=icam+1
       call rewftn(imos)
      if(icam.eq.iblk) iend=ires
      ig=0
      kz=-jd
      do 832 lc=1,7
      kz=kz+jd
      do 831 l=1,iblk
      if(l.eq.icam) go to 833
      read(imos)
      go to 831
833   read(imos)xdum
      lz=kz
      do 834 ld=1,iend
      lz=lz+1
834   bloc(lz)=xdum(ld)
831   continue
832    continue
830   ig=ig+1
      if (dabs(cp).gt.finthr)
     +  write(iwr,505)j,k,cp,xm(ig),ym(ig),zm(ig),
     +                  xam(ig),yam(ig),zam(ig),ovm(ig)
  505 format(2x,2i5,8f15.6)
      if (j.eq.k) go to 503
      jaf=ilifq(k)+j
      cq=q(jaf)
      if (dabs(cq).gt.finthr)
     +  write(iwr,505) k,j,cq,xm(ig),ym(ig),zm(ig),xam(ig),
     +                   yam(ig),zam(ig),ovm(ig)
      cplu=cp+cq
      cmin=cp-cq
      x1g=x1g+xm(ig)*cplu
      x2g=x2g+ym(ig)*cplu
      x3g=x3g+zm(ig)*cplu
      y1g=y1g+xam(ig)*cmin
      y2g=y2g+yam(ig)*cmin
      y3g=y3g+zam(ig)*cmin
      go to 502
503   x1g=x1g+xm(ig)*cp
      x2g=x2g+ym(ig)*cp
      x3g=x3g+zm(ig)*cp
502   continue
      delc=delc*delc*delc*amult
      x1h=x1g*x1g
      y1h=x2g*x2g
      z1h=x3g*x3g
      x1q=delc*x1h
      y1q=delc*y1h
      z1q=delc*z1h
      dbt=x1h+y1h+z1h
      dct=dbt*delc
      write (iwr,2010)
 2010 format(/8x,3h(x),12x,3h(y),12x,3h(z))
      write(iwr,700)x1g,x2g,x3g
700   format(/2x,7f15.6)
      write (iwr,2015)
 2015 format(/8x,4h(x)2,11x,4h(y)2,11x,4h(z)2,12x,3hsum )
      write(iwr,700) x1h,y1h,z1h,dbt
      write(iwr,701)
 701  format(/
     *' dipole transition probability and life-time (tau)'
     */1x,49('-'))
      tau=1.0d0/dct
      write(iwr,702)x1q,y1q,z1q,dct,tau
 702  format(/
     *'      dx       ',5x,e15.6/
     *'      dy       ',5x,e15.6/
     *'      dz       ',5x,e15.6//
     *' dtotal (sec-1)',5x,e15.6//
     *' tau (sec)     ',5x,e15.6/)
      fun=del/1.5d0
      x1h=fun*x1h
      y1h=fun*y1h
      z1h=fun*z1h
      ftt=fun*dbt
      write (iwr,2016)
 2016 format(' oscillator strengths -- dipole length operator'/
     *1x,46('-')//
     *9x,'fx',13x,'fy',13x,'fz',13x,'f(r)')
      write(iwr,700) x1h,y1h,z1h,ftt
      write (iwr,2060)
 2060 format(/8x,'delx',11x,'dely',11x,'delz')
      write(iwr,700) y1g,y2g,y3g
      y1g=y1g*y1g
      y2g=y2g*y2g
      y3g=y3g*y3g
      ftt=y1g+y2g+y3g
      write (iwr,2061)
 2061 format(/8x,5hdelx2,10x,5hdely2,10x,5hdelz2,11x,3hsum)
      write(iwr,700) y1g,y2g,y3g,ftt
      write (iwr,2062)
 2062 format(//' oscillator strengths -- dipole velocity operator'
     +/1x,48('-')//
     *8x,'fx(del)',8x,'fy(del)',8x,'fz(del)',11x,'f(del)')
      fun=1.0d0/(1.5d0*del)
      y1g=fun*y1g
      y2g=fun*y2g
      y3g=fun*y3g
      ftt=fun*ftt
      write(iwr,700) y1g,y2g,y3g,ftt
200   continue
      call rewftn(iscr)
      call rewftn(iaos)
      call rewftn(imos)
      if(iput.ne.jput) call rewftn(jput)
      cpu=cpulft(1)
      write(iwr,8002)cpu ,charwall()
 8002 format(/' *** end of moment calculation at ',f8.2,
     *' seconds',a10,' wall'/)
      return
750   write(iwr,758)iorbs,m,ksum
 758  format(//
     *' dimensioning error has occurred'//
     *' no. of basis functions                   ',i3//
     =' no. of correlated electrons              ',i3//
     *' no. of core+active orbitals              ',i3//)
      call caserr('dimensioning error has occurred')
751   write(iwr,759) jmo,imo
 759  format(//
     *' inconsistent number of active orbitals'/
     *' no. of active orbitals in excited state ',i3//
     *' no. of active orbitals in  ground state ',i3//)
      call caserr('inconsistent number of active orbitals')
752   write(iwr,760) icore,jcore
 760  format(//
     *' inconsistent number of frozen core orbitals'/
     *' no. of frozen orbitals in ground state ',i3//
     *' no. of frozen orbitals in excited state',i3//)
      call caserr('inconsistent no. of frozen core orbitals')
753   write(iwr,761) iorbs,jorbs
 761  format(//
     *' inconsistent number of a.o. basis functions'//
     *' number of a.o.s in ground state ',i3//
     *' number of a.o.s in excited state',i3//)
      call caserr('inconsistent number of a.o. basis functions')
754   write(iwr,762)jsum,ksum
 762  format(//
     *' inconsistent number of core+active orbitals'//
     *' number in ground state ',i3//
     *' number in excited state',i3//)
      call caserr('inconsistent no. of core+active orbitals')
755   write(iwr,763)m2,m
 763  format(//
     *' inconsistent number of correlated electrons'//
     *' number of active electrons in excited state',i3//
     *' number of active electrons in  ground state',i3//)
      call caserr('inconsistent number of correlated electrons')
757   write(iwr,765) i,cf(i),df(i)
 765  format(//
     *' inconsistency in eigen vectors for ground and excited state'//
     *' element',i6,/
     *' coefficient in ground state vector ',f20.8//
     *' coefficient in excited state vector',f20.8//)
      call caserr('inconsistent ground and excited state vectors')
      return
      end
      subroutine tm2
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
      common/aaaa/icom(21)
      common/scrtch/bc(2550),cbb(2550),x(255),p(255),r(255),
     *ic(256),ii(256),ij(256),ik(256),ispa4(256),ilife(256),
     *jmatsp(2855),
     *ovl(500),dixf(500),diyf(500),dizf(500),xnab(500),ynab(500),
     *znab(500),spc(19000),
     *icf(22500),ibal(8),itil(8),mcomp(256),kj(8),lsym(2040),
     *nj(8),ntil(8),nbal(9),ncomp(256),nzer(256),xdum(1300),
     *nconf(5),jconf(5),nytl(5),ndub(5),jtest(mxnshl),
     *nplu(5),itest(mxnshl),jytl(5),jdub(5),jplu(5),lj(8),imap(504),
     *ihog(48),jkan(mxcrc2),jmap(504),jhog(48),ikan(mxcrc2),
     +jcon(maxorb),imat(3780),irm(126),jrm(mxcrec),
     +intp(maxorb),intn(maxorb),
     +zeta(2550),d(2550),xorb(255),yorb(255),
     *zorb(255),anorm(255),psep(6),pvec(3),dif(3),dif2(3),
     *wc(100),wx(100),wy(100),wz(100),t1(5),t2(5),t3(5),
     *zxf(11),zyf(11),zzf(11),alf(11),bf(5),af(6),
     *f(mxcrec),g(mxcrc2),fj(mxcrec),gj(mxcrc2),ci(126),cj(126),
     *di(48),dj(48)
c
      common/ftape/iaos,imos,iscr,idum4(17)
      common/craypk/iput,idx,jput,jdx,nstate
      common/lsort/fac,acore,enee,ilifq(maxorb),iorbs,iswh,m,m2,mu,id
     *,ni8,nj8,nk8,imo,ig
      dimension cf(22500),bloc(9100),ovm(1300),
     * xm(1300),ym(1300),zm(1300),xam(1300),yam(1300),zam(1300),
     * e(7401),
     * jmat(30000),dihf(3500)
c
      equivalence (e(1),bc(1),bloc(1),jmat(1))
      equivalence (cf(1),ovl(1)),(dihf(1),ovl(1)),
     * (bloc(1),ovm(1)),
     * (bloc(1301),xm(1)),(bloc(2601),ym(1)),(bloc(3901),
     *zm(1)),(bloc(5201),xam(1)),(bloc(6501),yam(1)),
     * (bloc(7801),zam(1))
c
      thrs=140.0d0
      mm=0
      m2=m+2
      mu=m-2
      do 101 k=1,iorbs
      lco=ic(k)
      mi2=ii(k)
      t1f=fac
      if(mi2.gt.0) t1f=t1f*t2(mi2)
      mj2=ij(k)
      if(mj2.gt.0) t1f=t1f*t2(mj2)
      mk2=ik(k)
      if(mk2.gt.0) t1f=t1f*t2(mk2)
      mmmax=mi2
      if(mmmax.lt.mj2) mmmax=mj2
      if(mmmax.lt.mk2) mmmax=mk2
      mi3=mi2+1
      mj3=mj2+1
      mk3=mk2+1
      fun=mi2+mj2+mk2
      farx=0.5d0*fun+0.75d0
       gun=fun*0.5d0
      mi4=mi2*(mi2-1)/2
      mj4=mj2*(mj2-1)/2
      mk4=mk2*(mk2-1)/2
      x0=x(k)
      y0=p(k)
      z0=r(k)
      kkk = ilife(k)
      kkkk = kkk
      do 4 l=1,lco
      zeta(l+kkkk)=bc(kkk+l)
   4  d(l+kkkk)   =cbb(kkk+l)
      xorb(k)=x0
      yorb(k)=y0
      zorb(k)=z0
      sum=0.0d0
      do 550 l=1,lco
      a1=zeta(l+kkkk)
      ga=d(l+kkkk)*(a1**(0.75d0+gun))
      do 550 l2=1,lco
      b=zeta(l2+kkkk)
      gb=dsqrt(b)/(a1+b)
  550 sum=sum+ga*d(l2+kkkk)*(gb**(1.5d0+fun))
      sum=sum*fac*(2.0d0**fun)
      anorm(k)=1.0d0/dsqrt(sum)
      t1f=t1f*anorm(k)
      do 102 l=1,k
      mm=mm+1
      llll = ilife(l)
      ni3=ii(l)
      nj3=ij(l)
      nk3=ik(l)
      suma=0.0d0
      sumb=0.0d0
      sumc=0.0d0
      sumd=0.0d0
      sume=0.0d0
      sumf=0.0d0
      sumg=0.0d0
      t2f=t1f*anorm(l)
      if(ni3.gt.0) t2f=t2f*t2(ni3)
      if(nj3.gt.0) t2f=t2f*t2(nj3)
      if(nk3.gt.0) t2f=t2f*t2(nk3)
      fbrx=0.5d0*(ni3+nj3+nk3)+0.75d0
      nax=ni3
      if(nax.lt.nj3) nax=nj3
      if(nax.lt.nk3) nax=nk3
      nax=nax+1
      x5=xorb(l)
      y5=yorb(l)
      z5=zorb(l)
      ni4=ni3+1
      ni5=ni4+1
      if(ni3.lt.2) go to 143
      ni8=(ni3-1)*(ni3-2)/2-1
  143 ni6=ni3*(ni3-1)/2
      ni7=ni4*(ni4-1)/2-1
      nj4=nj3+1
      nj5=nj4+1
      if(nj3.lt.2) go to 144
      nj8=(nj3-1)*(nj3-2)/2-1
  144 nj6=nj3*(nj3-1)/2
      nj7=nj4*(nj4-1)/2-1
      nk4=nk3+1
      nk5=nk4+1
      if(nk3.lt.2) go to 145
      nk8=(nk3-1)*(nk3-2)/2-1
  145 nk6=nk3*(nk3-1)/2
      nk7=nk4*(nk4-1)/2-1
      dif(1)=xorb(l)-x0
      dif(2)=yorb(l)-y0
      dif(3)=zorb(l)-z0
      difsum=0.0d0
      do 10 i1=1,3
      dif2(i1)=dif(i1)*dif(i1)
  10  difsum=difsum+dif2(i1)
      lc1=ic(l)
      nxx=mi2+ni3+1
      nxy=mj2+nj3+1
      nxz=mk2+nk3+1
      nox=nxx
      if(nox.lt.nxy) nox=nxy
      if(nox.lt.nxz) nox=nxz
      bax=dif(1)
      bay=dif(2)
      baz=dif(3)
      zxf(1)=bax
      zyf(1)=bay
      zzf(1)=baz
      if(nxx.eq.1) go to 129
      do 130 i=2,nxx
  130 zxf(i)=zxf(i-1)*bax
  129 if(nxy.eq.1) go to 131
      do 132 i=2,nxy
  132 zyf(i)=zyf(i-1)*bay
  131 if(nxz.eq.1) go to 133
      do 134 i=2,nxz
  134 zzf(i)=zzf(i-1)*baz
  133 do 3 l1=1,lco
      aa=zeta(l1+kkkk)
      qx=aa*x0
      qy=aa*y0
      qz=aa*z0
      tuffy=d(l1+kkkk)*(aa**farx)
      am=-aa
      af(1)=am
      aab2=difsum*aa
      if(nax.eq.1) go to 135
      do 136 i=2,nax
  136 af(i)=am*af(i-1)
  135 do 3 l2=1,lc1
      bb=zeta(l2+llll)
      bb2=bb+bb
      alp=aa+bb
      alg=1.0d0/alp
      arg=aab2*bb*alg
      if(arg.gt.thrs) go to 3
      tuf=tuffy*dexp(-arg)
      tuf=tuf*d(l2+llll)*(bb**fbrx)
      if(mmmax.eq.0) go to 141
      bf(1)=bb
      if(mmmax.eq.1) go to 141
      do 907 i=2,mmmax
  907 bf(i)=bf(i-1)*bb
  141 tuf=tuf*(alg**1.5d0)
      alf(1)=alg
      if(nox.eq.1) go to 150
      do 151 i=2,nox
  151 alf(i)=alf(i-1)*alg
  150 sxma=0.0d0
      sxmb=0.0d0
      sxmc=0.0d0
      sxmd=0.0d0
      syma=0.0d0
      symb=0.0d0
      symc=0.0d0
      symd=0.0d0
      szma=0.0d0
      szmb=0.0d0
      szmc=0.0d0
      szmd=0.0d0
      px=(qx+bb*x5)*alg
      py=(qy+bb*y5)*alg
      pz=(qz+bb*z5)*alg
      mi5=mi4
      do 140 i=1,mi3
      w4=1.0d0
      if(i.eq.1) go to 142
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  142 ni9=ni6
      do 146 j=1,ni4
      w5=w4
      if(j.eq.1) go to 147
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  147 i9=i+j-2
      if(i9.gt.0) w5=w5*zxf(i9)
      i8=nxx-i9
      i7=i8/2
      i4=i7+i9
      if(i4.gt.0) w5=w5*alf(i4)
      if(i8-i7*2.eq.0) go to 153
      if(i7.gt.0) w5=w5*t3(i7)
      sxma=sxma+w5
      go to 146
  153 sxmb=sxmb+w5*t3(i7)
  146 continue
c     next comes fx(m+1)  ***********
      ni9=ni7
      do 154 j=1,ni5
      ni9=ni9+1
      i9=i+j-2
      i8=nxx-i9
      i7=i8/2
      if(i8-i7*2.ne.0) go to 154
      w5=w4
      if(i7.eq.0) go to 800
      w5=w5*t3(i7)*alf(i7+i9)
      go to 801
  800 if(i9.gt.0) w5=w5*alf(i9)*zxf(i9)
      go to 802
  801 if(i9.gt.0) w5=w5*zxf(i9)
  802 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      sxmc=sxmc+w5
 154  continue
      if(ni3.eq.0) go to 140
      ni9=ni8
      do 157 j=1,ni3
      ni9=ni9+1
      i9=i+j-2
      i8=nxx-i9-2
      i7=i8/2
      if(i8-i7*2.ne.0) go to 157
      w5=w4
      if(i7.eq.0) go to 158
      w5=w5*t3(i7)*alf(i7+i9)
      go to 159
  158 if(i9.gt.0) w5=w5*alf(i9)*zxf(i9)
      go to 803
  159 if(i9.gt.0) w5=w5*zxf(i9)
  803 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      sxmd=sxmd+w5
  157 continue
  140 continue
      mi5=mj4
      do 160 i=1,mj3
      w4=1.0d0
      if(i.eq.1) go to 162
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  162 ni9=nj6
      do 166 j=1,nj4
      w5=w4
      if(j.eq.1) go to 167
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  167 i9=i+j-2
      if(i9.gt.0) w5=w5*zyf(i9)
      i8=nxy-i9
      i7=i8/2
      i4=i7+i9
      if(i4.gt.0) w5=w5*alf(i4)
      if(i8-i7*2.eq.0) go to 173
      if(i7.gt.0) w5=w5*t3(i7)
      syma=syma+w5
      go to 166
  173 symb=symb+w5*t3(i7)
  166 continue
      ni9=nj7
      do 174 j=1,nj5
      ni9=ni9+1
      i9=i+j-2
      i8=nxy-i9
      i7=i8/2
      if(i8-i7*2.ne.0) go to 174
      w5=w4
      if(i7.eq.0) go to 804
      w5=w5*t3(i7)*alf(i7+i9)
      go to 805
  804 if(i9.gt.0) w5=w5*alf(i9)*zyf(i9)
      go to 806
  805 if(i9.gt.0) w5=w5*zyf(i9)
  806 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      symc=symc+w5
  174 continue
      if(nj3.eq.0) go to 160
      ni9=nj8
      do 177 j=1,nj3
      ni9=ni9+1
      i9=i+j-2
      i8=nxy-i9-2
      i7=i8/2
      if(i8-i7*2.ne.0) go to 177
      w5=w4
      if(i7.eq.0) go to 178
      w5=w5*t3(i7)*alf(i7+i9)
      go to 179
  178 if(i9.gt.0) w5=w5*alf(i9)*zyf(i9)
      go to 807
  179 if(i9.gt.0) w5=w5*zyf(i9)
  807 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      symd=symd+w5
  177 continue
  160 continue
      mi5=mk4
      do 180 i=1,mk3
      w4=1.0d0
      if(i.eq.1) go to 182
      mi5=mi5+1
      w4=icom(mi5)*bf(i-1)
  182 ni9=nk6
      do 186 j=1,nk4
      w5=w4
      if(j.eq.1) go to 187
      ni9=ni9+1
      w5=w5*icom(ni9)*af(j-1)
  187 i9=i+j-2
      if(i9.gt.0) w5=w5*zzf(i9)
      i8=nxz-i9
      i7=i8/2
      i4=i7+i9
      if(i4.gt.0) w5=w5*alf(i4)
      if(i8-i7*2.eq.0) go to 193
      if(i7.gt.0) w5=w5*t3(i7)
      szma=szma+w5
      go to 186
  193 szmb=szmb+w5*t3(i7)
  186 continue
      ni9=nk7
      do 194 j=1,nk5
      ni9=ni9+1
      i9=i+j-2
      i8=nxz-i9
      i7=i8/2
      if(i8-i7*2.ne.0) go to 194
      w5=w4
      if(i7.eq.0) go to 808
      w5=w5*t3(i7)*alf(i7+i9)
      go to 809
  808 if(i9.gt.0) w5=w5*alf(i9)*zzf(i9)
      go to 810
  809 if(i9.gt.0) w5=w5*zzf(i9)
  810 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      szmc=szmc+w5
  194 continue
      if(nk3.eq.0) go to 180
      ni9=nk8
      do 197 j=1,nk3
      ni9=ni9+1
      i9=i+j-2
      i8=nxz-i9-2
      i7=i8/2
      if(i8-i7*2.ne.0) go to 197
      w5=w4
      if(i7.eq.0) go to 198
      w5=w5*t3(i7)*alf(i7+i9)
      go to 199
  198 if(i9.gt.0) w5=w5*alf(i9)*zzf(i9)
      go to 811
  199 if(i9.gt.0) w5=w5*zzf(i9)
  811 if(j.gt.1) w5=w5*icom(ni9)*af(j-1)
      szmd=szmd+w5
  197 continue
  180 continue
      sumg=sumg+tuf*sxma*syma*szma
      suma=suma+tuf*(sxmb+px*sxma)*syma*szma
      sumb=sumb+tuf*sxma*(symb+py*syma)*szma
      sumc=sumc+tuf*sxma*syma*(szmb+pz*szma)
      sumd=sumd+tuf*(ni3*sxmd-bb2*sxmc)*syma*szma
      sume=sume+tuf*(nj3*symd-bb2*symc)*sxma*szma
      sumf=sumf+tuf*(nk3*szmd-bb2*szmc)*sxma*syma
   3  continue
      ovl(mm)=t2f*sumg
      dixf(mm)=t2f*suma
      diyf(mm)=t2f*sumb
      dizf(mm)=t2f*sumc
      xnab(mm)=-t2f*sumd
      ynab(mm)=-t2f*sume
      znab(mm)=-t2f*sumf
       if(mm.lt.id) go to 102
       write(iaos) dihf
      mm=0
 102   continue
  101 continue
           write(iaos) dihf
      return
      end
      subroutine tm3(q,lensq)
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
      common/aaaa/icom(21)
      common/scrtch/bc(2550),cbb(2550),x(255),p(255),r(255),
     *ic(256),ii(256),ij(256),ik(256),ispa4(256),ilife(256),
     *jmatsp(2855),
     *ovl(500),dixf(500),diyf(500),dizf(500),xnab(500),ynab(500),
     *znab(500),spc(19000),
     *icf(22500),ibal(8),itil(8),mcomp(256),kj(8),lsym(2040),
     *nj(8),ntil(8),nbal(9),ncomp(256),nzer(256),xdum(1300),
     *nconf(5),jconf(5),nytl(5),ndub(5),jtest(mxnshl),
     *nplu(5),itest(mxnshl),jytl(5),jdub(5),jplu(5),lj(8),imap(504),
     *ihog(48),jkan(mxcrc2),jmap(504),jhog(48),ikan(mxcrc2),
     +jcon(maxorb),imat(3780),irm(126),jrm(mxcrec),
     +intp(maxorb),intn(maxorb),
     +zeta(2550),d(2550),xorb(255),yorb(255),
     *zorb(255),anorm(255),psep(6),pvec(3),dif(3),dif2(3),
     *wc(100),wx(100),wy(100),wz(100),t1(5),t2(5),t3(5),
     *zxf(11),zyf(11),zzf(11),alf(11),bf(5),af(6),
     *f(mxcrec),g(mxcrc2),fj(mxcrec),gj(mxcrc2),ci(126),cj(126),
     +di(48),dj(48)
c
      common/ftape/iaos,imos,iscr,idum4(17)
      common/craypk/iput,idx,jput,jdx,nstate
      common/lsort/fac,acore,enee,ilifq(maxorb),iorbs,iswh,m,m2,mu,id
     *,ni8,nj8,nk8,imo,ig
      dimension q(lensq)
      dimension cf(22500),bloc(9100),ovm(1300),
     * xm(1300),ym(1300),zm(1300),xam(1300),yam(1300),zam(1300),
     * e(7401),
     * jmat(30000),dihf(3500)
c
      equivalence (e(1),bc(1),bloc(1),jmat(1))
      equivalence (cf(1),ovl(1)),(dihf(1), ovl(1)),
     * (bloc(1),ovm(1)),
     * (bloc(1301),xm(1)),(bloc(2601),ym(1)),(bloc(3901),
     *zm(1)),(bloc(5201),xam(1)),(bloc(6501),yam(1)),
     * (bloc(7801),zam(1))
c
      do 201 j=1,iswh
      jc=jconf(j)
      if (jc.eq.0) go to 201
      read(jput)jhb,jmax,ndj,kmj,jmap,jhog,enee
      nlj=jytl(j)
      jmns=j-1
      jqns=j-2
      jps=jplu(j)
      jdb=jdub(j)
      jdbq=jdb+jdb
      nodj=jmns+jps
      nodj2=nodj-m2
      jcl=0
      idel=ndj*m
      do 202 k=1,jhb
      jcl=jcl+jmax
      if(jcl.gt.jc) jmax=jc+jmax-jcl
      read(jput) jkan,fj,gj
      nxj=0
      call rewftn(iscr)
      k5=0
      k4=0
      do 111 k6=1,jmax
      do 112 l=1,nlj
       nxj=nxj+1
 112  jtest(l)=jkan(nxj)
      if(jmns .gt. 0) go to 213
      jrm(k6)=1
      if(jps.eq.0) go to 214
      do 215 l=1,jps
      k4=k4+1
 215  jmat(k4)=jtest(l)
      if(jdb.eq.0) go to 111
 214  kx=jps
      do 217 l=1,jdb
      kx=kx+1
      k4=k4+1
      lk=jtest(kx)
      jmat(k4)=lk
      k4=k4+1
 217  jmat(k4)=-lk
      go to 111
213   jq=1
      lx=k4
      do 218 l=1,ndj
      k5=k5+1
      jrm(k5)=1
      jz=jq+jqns
      jy=jmap(jq)
      do 219 kk=1,nodj
      lx=lx+1
      if (kk.ne.jy) go to 220
      if (jz.eq.jq) go to 221
      jq=jq+1
      jy=jmap(jq)
221   jmat(lx)=-jtest(kk)
      go to 219
220   jmat(lx)=jtest(kk)
219   continue
      jq=jz+1
218   lx=lx+jdbq
      if (jdb.eq.0) go to 113
      kx=nodj2+k4
      kk=nodj
      do 222 l=1,jdb
      kk=kk+1
      lk=jtest(kk)
      kx=kx+2
      mx=kx
      do 222 ll=1,ndj
      mx=mx+m
      jmat(mx+1)=lk
222   jmat(mx+2) =-lk
      k4=mx+2
      go to 111
  113 k4=lx
 111  continue
      if (j.eq.1) go to 300
      j9=j-1
      if(j.eq.2) go to 114
      j8=j-2
      do 115 l=1,j8
      nc=nconf(l)
      if(nc.eq.0) go to 115
      read(iscr)nhb
      do 116 ll=1,nhb
 116  read(iscr)
 115  continue
 114  nc=nconf(j9)
      if(nc.eq.0) go to 300
      read(iscr)nhb,imax,ndt,kml,imap
      nl=nytl(j9)
      nmns=j9-1
      nqns=j9-2
      nps=nplu(j9)
      ndb=ndub(j9)
      ndbq=ndb+ndb
      nod=nmns+nps
      nodm2=nod-m2
      icl=imax
      do 224 k9=1,nc
      if (icl.lt.imax) go to 225
      read(iscr) ikan,f
      icl=1
      nx=0
      if=0
      go to 226
225   icl=icl+1
226   do 227 l9=1,nl
      nx=nx+1
227   itest(l9)=ikan(nx)
      do 208 l=1,imo
 208  jcon(l)=0
      if(nod.eq.0) go to 209
      do 210 l=1,nod
      nt=itest(l)
 210  jcon(nt)=1
      if(ndb.eq.0) go to 211
 209  kk=nod
      do 212 l=1,ndb
      kk=kk+1
      nt=itest(kk)
 212  jcon(nt)=2
 211  do 233 kw=1,ndt
      if=if+1
 233  ci(kw)=f(if)
      if(nmns.gt.0) go to 235
      irm(1)=1
      if (nps.eq.0) go to 236
      do 237 kw=1,nps
237   imat(kw)=itest(kw)
      if(ndb.eq.0) go to 240
236   lx=nps
      kx=nps
      do 238 kw=1,ndb
      kx=kx+1
      lx=lx+1
      lk=itest(kx)
      imat(lx)=lk
      lx=lx+1
238   imat(lx)=-lk
      go to 240
235   lx=0
      jq=1
      do 239 kw=1,ndt
      irm(kw)=1
      jz=jq+nqns
      jy=imap(jq)
      do 241 lw=1,nod
      lx=lx+1
      if (lw.ne.jy) go to 242
      if (jz.eq.jq) go to 243
      jq=jq+1
      jy=imap(jq)
243   imat(lx)=-itest(lw)
      go to 241
242   imat(lx)=itest(lw)
241   continue
      jq=jz+1
239   lx=lx+ndbq
      if (ndb.eq.0) go to 240
      kx=nodm2
      kk=nod
      do 244 kw=1,ndb
      kk=kk+1
      lk=itest(kk)
      kx=kx+2
      mx=kx
      do 244 lw=1,ndt
      mx=mx+m
      imat(mx+1)=lk
244   imat(mx+2)=-lk
 240  nxj=0
      injs=-ndj
      igj=0
      do 117 l=1,jmax
      injs=injs+ndj
      do 118 ll=1,nlj
      nxj=nxj+1
 118  jtest(ll)=jkan(nxj)
      nix=0
      if(nodj.eq.0) go to 228
      do 229 kw=1,nodj
      jt=jtest(kw)
      if(jcon(jt).gt.0) go to 229
      if(nix.eq.1) go to 230
      nix=1
  229 continue
      if(jdb.eq.0) go to 231
 228  kk=nodj
      do 232 kw=1,jdb
      kk=kk+1
      jt=jtest(kk)
      jb=jcon(jt)
      if(jb.eq.2) go to 232
      if(jb.eq.0.or.nix.eq.1) go to 230
      nix=1
 232  continue
 231  do 206     ll=1,kmj
      igj=igj+1
 206  dj(ll)=gj(igj)
      go to 119
  230 igj=igj+kmj
      go to 117
 119  if(ndt.gt.kmj) go to 280
      ini=-m
      do 245 lw=1,ndt
      orb=ci(lw)
      ini=ini+m
      do 246 kw=1,imo
      intp(kw)=0
246   intn(kw)=0
      ink=ini
      do 247 kw=1,m
      ink=ink+1
      lx=imat(ink)
      if(lx.lt.0) go to 248
      intp(lx)=kw
      go to 247
248   intn(-lx)=kw
247   continue
      do 249 kw=1,kmj
      nix=0
      jr=jhog(kw)+injs
      inj=(jr-1)*m
      ink=inj
      do 250 mw=1,m
      ink=ink+1
      la=jmat(ink)
      if (la.lt.0) go to 251
      if (intp(la).gt.0) go to 250
252   if (nix.eq.1) go to 249
      nix=1
      iodd=mw
      go to 250
251   if (intn(-la).eq.0) go to 252
250   continue
      if (iodd.eq.1) go to 253
      kx=1
      ky=iodd-1
      ink=inj
261   do 254 mw=kx,ky
      ink=ink+1
      lx=jmat(ink)
      if(lx.lt.0) go to 255
      la=intp(lx)
      go to 256
255   la=intn(-lx)
256   if (la.eq.mw) go to 254
      ix=ini+la
      iy=ini+mw
      ni=imat(iy)
      imat(ix)=ni
      imat(iy)=lx
      irm(lw)=-irm(lw)
      if (lx.lt.0) go to 257
      intp(lx)=mw
      go to 258
257   intn(-lx)=mw
258   if(ni.lt.0) go to 259
      intp(ni)=la
      go to 254
259   intn(-ni)=la
254   continue
      if (ky.gt.mu) go to 260
253   kx=iodd+1
      ky=m
      ink=inj+iodd
      go to 261
260   vm=orb*dj(kw)
      if (irm(lw).ne.jrm(jr))vm=-vm
      ink=ini+iodd
      ia=imat(ink)
      ink=inj+iodd
      ja=jmat(ink)
      if(ia.gt.0) go to 262
      ia=-ia
      ja=-ja
262   ia=ilifq(ia)+ja
      q(ia)=q(ia)+vm
249   continue
245   continue
      go to 117
280   do 281 lw=1,kmj
      jr=jhog(lw)+injs
      inj=(jr-1)*m
      orb=dj(lw)
      do 282 kw=1,imo
      intp(kw)=0
282   intn(kw)=0
      ink=inj
      do 283 kw=1,m
      ink=ink+1
      lx=jmat(ink)
      if (lx.lt.0) go to 284
      intp(lx)=kw
      go to 283
284   intn(-lx)=kw
283   continue
      ini=-m
      do 285 kw=1,ndt
      ini=ini+m
      ink=ini
      nix=0
      do 286 mw=1,m
      ink=ink+1
      la=imat(ink)
      if (la.lt.0) go to 287
      if (intp(la).gt.0) go to 286
288   if (nix.eq.1) go to 285
      nix=1
      iodd=mw
      go to 286
287   if (intn(-la).eq.0) go to 288
286   continue
      if (iodd.eq.1) go to 289
      kx=1
      ky=iodd-1
      ink=ini
290   do 291 mw=kx,ky
      ink=ink+1
      lx=imat(ink)
      if (lx.lt.0) go to 292
      la=intp(lx)
      go to 293
292   la=intn(-lx)
293   if(la.eq.mw) go to 291
      ix=inj+la
      iy=inj +mw
      ni=jmat(iy)
      jmat(ix)=ni
      jmat(iy)=lx
      jrm(jr)=-jrm(jr)
      if (lx.lt.0) go to 294
      intp(lx)=mw
      go to 295
294   intn(-lx)=mw
295   if (ni.lt.0) go to 296
      intp(ni)=la
      go to 291
296   intn(-ni)=la
291   continue
      if (ky.gt.mu) go to 297
289   kx=iodd+1
      ky=m
      ink=ini+iodd
      go to 290
297   vm=orb*ci(kw)
      if (irm(kw).ne.jrm (jr))vm=-vm
      ink=iodd+ini
      ia=imat(ink)
      ink=iodd+inj
      ja=jmat(ink)
      if (ia.gt.0) go to 298
      ia=-ia
      ja=-ja
298   ia=ilifq(ia)+ja
      q(ia)=q(ia)+vm
285   continue
281   continue
 117  continue
 224  continue
 300  jpig=j+1
      if(j.eq.iswh) jpig=j
      do 323 j9=j,jpig
      nc=nconf(j9)
      if(nc.eq.0) go to 323
      read(iscr)nhb,imax,ndt,kml,imap,ihog
      nl=nytl(j9)
      nmns=j9-1
      nqns=j9-2
      nps=nplu(j9)
      ndb=ndub(j9)
      ndbq=ndb+ndb
      nod=nmns+nps
      nodm2=nod-m2
      icl=imax
      do 324 k9=1,nc
      if (icl.lt.imax) go to 325
      read(iscr) ikan,f,g
      icl=1
      nx=0
      ig=0
      go to 326
325   icl=icl+1
326   do 327 l9=1,nl
      nx=nx+1
327   itest(l9)=ikan(nx)
      do 120 l=1,imo
 120  jcon(l)=0
      if(nod.eq.0) go to 123
      do 121 l=1,nod
      nt=itest(l)
 121  jcon(nt)=1
      if(ndb.eq.0) go to 122
 123  kk=nod
      do 124 l=1,ndb
      kk=kk+1
      nt=itest(kk)
 124  jcon(nt)=2
  122  do 125 kw=1,kml
      ig=ig+1
 125  di(kw)=g(ig)
      if(nmns.gt.0) go to 335
      irm(1)=1
      if (nps.eq.0) go to 336
      do 337 kw=1,nps
337   imat(kw)=itest(kw)
      if(ndb.eq.0) go to 340
336   lx=nps
      kx=nps
      do 338 kw=1,ndb
      kx=kx+1
      lx=lx+1
      lk=itest(kx)
      imat(lx)=lk
      lx=lx+1
338   imat(lx)=-lk
      go to 340
335   lx=0
      do 339 kw=1,kml
      irm(kw)=1
      jq=(ihog(kw)-1)*nmns+1
      jz=jq+nqns
      jy=imap(jq)
      do 341 lw=1,nod
      lx=lx+1
      if (lw.ne.jy) go to 342
      if (jz.eq.jq) go to 343
      jq=jq+1
      jy=imap(jq)
343   imat(lx)=-itest(lw)
      go to 341
342   imat(lx)=itest(lw)
341   continue
339   lx=lx+ndbq
      if (ndb.eq.0) go to 340
      kx=nodm2
      kk=nod
      do 344 kw=1,ndb
      kk=kk+1
      lk=itest(kk)
      kx=kx+2
      mx=kx
      do 344 lw=1,kml
      mx=mx+m
      imat(mx+1)=lk
344   imat(mx+2)=-lk
 340  nxj=0
      ifj=0
      jsig=-m-idel
      do 126 l=1,jmax
      jsig=jsig+idel
      do 127 ll=1,nlj
      nxj=nxj+1
 127  jtest(ll)=jkan(nxj)
      nix=0
      if(nodj.eq.0) go to 328
      do 329 kw=1,nodj
      jt=jtest(kw)
      if(jcon(jt).gt.0) go to 329
      if(nix.eq.1) go to 330
      nix=1
 329  continue
      if(jdb.eq.0) go to 331
 328  kk=nodj
      do 332 kw=1,jdb
      kk=kk+1
      jt=jtest(kk)
      jb=jcon(jt)
      if(jb.eq.2) go to 332
      if(jb.eq.0.or.nix.eq.1) go to 330
      nix=1
 332  continue
 331  injs=ifj
      do 333 kw=1,ndj
      ifj=ifj+1
 333  cj(kw)=fj(ifj)
      go to 334
 330  ifj=ifj+ndj
      go to 126
334   if (nix.eq.1) go to 434
      vm=0.0d0
      do 400 kw=1,kml
      jr=jhog(kw)
400   vm=vm+cj(jr)*di(kw)
      acore=acore+vm
      if (nod.eq.0) go to 402
      do 401 kw=1,nod
      jt=itest(kw)
      ia=ilifq(jt)+jt
401   q(ia)=q(ia)+vm
      if(ndb.eq.0) go to 126
402   kk=nod
      vm=vm+vm
      do 403 kw=1,ndb
      kk=kk+1
      jt=itest(kk)
      ia=ilifq(jt)+jt
403   q(ia)=q(ia)+vm
      go to 126
 434  if(ndj.gt.kml) go to 380
      inj=jsig
      do 345 lw=1,ndj
      orb=cj(lw)
      inj=inj+m
      do 346 kw=1,imo
      intp(kw)=0
346   intn(kw)=0
      ink=inj
      do 347 kw=1,m
      ink=ink+1
      lx=jmat(ink)
      if(lx.lt.0) go to 348
      intp(lx)=kw
      go to 347
348   intn(-lx)=kw
347   continue
      ini=-m
      do 349 kw=1,kml
      nix=0
      ini=ini+m
      ink=ini
      do 350 mw=1,m
      ink=ink+1
      la=imat(ink)
      if (la.lt.0) go to 351
      if (intp(la).gt.0) go to 350
352   if(nix.eq.1) go to 349
      nix=1
      iodd=mw
      go to 350
351   if(intn(-la).eq.0) go to 352
350   continue
      if (iodd.eq.1) go to 353
      kx=1
      ky=iodd-1
      ink=ini
361   do 354 mw=kx,ky
      ink=ink+1
      lx=imat(ink)
      if(lx.lt.0) go to 355
      la=intp(lx)
      go to 356
355   la=intn(-lx)
356   if(la.eq.mw) go to 354
      ix=inj+la
      iy=inj+mw
      ni=jmat(iy)
      jmat(ix)=ni
      jmat(iy)=lx
       jrm(lw+injs)=-jrm(lw+injs)
      if(lx.lt.0) go to 357
      intp(lx)=mw
      go to 358
357   intn(-lx)=mw
358   if(ni.lt.0) go to 359
      intp(ni)=la
      go to 354
359   intn(-ni)=la
354   continue
      if(ky.gt.mu) go to 360
353   kx=iodd+1
      ky=m
      ink=ini+iodd
      go to 361
360   vm=orb*di(kw)
       if(irm(kw).ne.jrm(lw+injs)) vm=-vm
      ink=ini+iodd
      ia=imat(ink)
      ink=inj+iodd
      ja=jmat(ink)
      if(ia.gt.0) go to 362
      ia=-ia
      ja=-ja
362   ia=ilifq(ia)+ja
      q(ia)=q(ia)+vm
349   continue
345   continue
      go to 126
380   ini=-m
      do 381 lw=1,kml
      orb=di(lw)
      ini=ini+m
      do 382 kw=1,imo
      intp(kw)=0
382   intn(kw)=0
      ink=ini
      do 383 kw=1,m
      ink=ink+1
      lx=imat(ink)
      if(lx.lt.0) go to 384
      intp(lx)=kw
      go to 383
384   intn(-lx)=kw
383   continue
      inj=jsig
      do 385 kw=1,ndj
      inj=inj+m
      ink=inj
      nix=0
      do 386 mw=1,m
      ink=ink+1
      la=jmat(ink)
      if(la.lt.0) go to 387
      if (intp(la).gt.0) go to 386
388   if (nix.eq.1) go to 385
      nix=1
      iodd=mw
      go to 386
387   if (intn(-la).eq.0) go to 388
386   continue
      if (iodd.eq.1) go to 389
      kx=1
      ky=iodd-1
      ink=inj
390   do 391 mw=kx,ky
      ink=ink+1
      lx=jmat(ink)
      if (lx.lt.0) go to 392
      la=intp(lx)
      go to 393
392   la=intn(-lx)
393   if(la.eq.mw) go to 391
      ix=ini+la
      iy=ini+mw
      ni=imat(iy)
      imat(ix)=ni
      imat(iy)=lx
      irm(lw)=-irm(lw)
      if (lx.lt.0) go to 394
      intp(lx)=mw
      go to 395
394   intn(-lx)=mw
395   if(ni.lt.0) go to 396
      intp(ni)=la
      go to 391
396   intn(-ni)=la
391   continue
      if (ky.gt.mu) go to 397
389   kx=iodd+1
      ky=m
      ink=inj+iodd
      go to 390
397   vm=orb*cj(kw)
      if(irm(lw).ne.jrm(kw+injs)) vm=-vm
      ink=iodd+ini
      ia=imat(ink)
      ink=iodd+inj
      ja=jmat(ink)
      if(ia.gt.0) go to 398
      ia=-ia
      ja=-ja
398   ia=ilifq(ia)+ja
      q(ia)=q(ia)+vm
385   continue
381   continue
 126   continue
324   continue
323   continue
202   continue
201   continue
       return
       end
      subroutine writex(p,map,newbas,jtape)
      implicit real*8  (a-h,o-z), integer (i-n)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      dimension p(*),map(*)
      ind(i,j)=max(map(i),map(j))*(max(map(i),map(j))-1)/2
     +        +min(map(i),map(j))
      m5=12
      if(oprint(20)) m5=8
      m=1
      n=m5
   6  if(newbas.lt.m)return
      if(n.gt.newbas)n=newbas
      if(.not.oprint(20))write(jtape,200)(i,i=m,n)
      if(oprint(20))write(jtape,100)(i,i=m,n)
100   format(//3x,8i14)
200   format(//12i9)
      write(jtape,101)
101   format(/)
      do 1 j=1,newbas
      if(oprint(20))write(jtape,102)j,(p(ind(i,j)),i=m,n)
      if(.not.oprint(20))write(jtape,202)j,(p(ind(i,j)),i=m,n)
 1    continue
102   format(7x,i3,8f14.7)
202   format(1x,i3,12f9.4)
      m=m+m5
      n=n+m5
      goto 6
      end
      subroutine ver_mrdci4(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mrdci4.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
