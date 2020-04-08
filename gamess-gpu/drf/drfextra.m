c ======================================================================
c  this is needed to get drf things to work!
c  for now numerous values are hard-set
c
      subroutine inithondo(mode)
      implicit REAL  (a-h,o-z),integer  (i-n)
c changed:
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/darw)
INCLUDE(comdrf/dafil)
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/bas)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/mem)
INCLUDE(comdrf/mollab)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/scfopt)
INCLUDE(comdrf/runpar)
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/alpbet)
INCLUDE(comdrf/drfzfa)
c gamess:
INCLUDE(../m4/common/cslosc)
c
      REAL cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopenn
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopenn,nope,noe(10)
c
c the data statement is taken from 'master'
c the commonblock runlab in gen/mollab contains
c 'zsymm(7)' in which the third is actually 'zscftp'
      character*8 zrhf,zuhf,zgrhf,zgvb
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
c      from 'ctlnew.f'
      data maxbit /64/
c     ----------------------------------------------
c      in gen/scm the xscm is set :
c     ----------------------------------------------
c      data maxcor /8000/
c
c     comdrf/mem
c
      maxcor = 8000
c---------------------------------------------------mode0
      if (mode .eq. 0) then
         nzfa = (mxat+1)*mxgran
c        if ((zsymm(3).eq.zrhf).or.
c    1       (zsymm(3).eq.zuhf).or.
c    2       (zsymm(3).eq.zgrhf).or.
c    3       (zsymm(3).eq.zgvb)) then
         if (zsymm(4) .eq. 'ci') then
           wfntyp='        '
         else
           wfntyp='scf     '
         endif
         scftyp = zsymm(3)
         if ( odscf) hndtyp = 'direct'
c number of states is:
         iactst=1
         if (itol.le.0) itol=15
c
c-----   set hardware characteristic parameters
c           isingl = 1 for 64-bit word hardware
c           isingl = 2 for 32-bit word hardware
c           vector =.false. for scalar hardware
c           vector =.true.  for vector hardware
         isingl=lenwrd()
         nbits =maxbit/isingl
c        vector=.false.
c
c-----   define file numbers
c           ir: input data are read from this file
c           iw: output is written to this file
c           ipnch: vectors and other information is written to this file
c           ijk: integral file 1
c           ipk: integral file 2
         ipnch=7
         ijk=8
         ipk=9
c-----  correctly open the defined files
cafc dit is de plaats waar de lege vectors-, integrals1- en integrals2-
cafc file worden aangemaakt, skipping this doesn't seem to hurt
cafc
c         open (unit=ipnch, access='sequential', form='formatted',
c     +         status='unknown', file='vectors')
c         open (unit=ijk, access='sequential', form='unformatted',
c     +         status='unknown', file='integrals1')
c         open (unit=ipk, access='sequential', form='unformatted',
c     +         status='unknown', file='integrals2')
c-----  open idafh: general purpose da-file
         idafh=110
_IF(cray)
         ndar=512
_ELSE
         ndar=999
_ENDIF
         call daopen(idafh,ioda,navh,ndar)
c
cold         do 99 i=0,3000
         do 99 i=1,3*mxpts
            ia(i)=(i*(i-1))/2
   99    continue
      else
c----------------------------------------------mode=1
cmw         nprint=3
cxxx     print *,"inithond: array ia() hardcoded 3000..."
         call bfnorm
cxxx     print *,"inithond: anorm build"
         call bfnlab
cxxx     print *,"inithond: bflab build"
cxxx
         do i = 1, 325
           alphij(i) = alpha(i)
           betaij(i) = beta(i)
         enddo
         do i = 1, 10
           nopen(i) = no(i)
         enddo
         nopset = nseto
         npairs = npair
         if (scftyp .eq. 'gvb') then
           ncorb = nco
         else
           ncorb = nb
         endif
cxxx
      endif
c-------------------------------
cxxx  print *,"inithondo:ready"
      return
      end
c=======================================================================
c      subroutine dipcal(n,x0,x,q)
c=====================================================================
c      subroutine addamb(elnum,charge,cxyz,name,ibas,iecpx)
c      subroutine drfdst
c some functions like  and distab
c===================================================================
c      subroutine fdagal(f,al,fa,ia,n)
c----------------------------
c      function amas(name)
c=======================================================================
c      subroutine drfreda(a,redpol,n)
c======================================================================
c      subroutine shrink(b,a,n)
c=======================================================================
      subroutine bfnorm
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
INCLUDE(comdrf/bas)
      data zero,pt5,one /0.0d+00,0.5d+00,1.0d+00/
      data pt75 /0.75d+00/
      data pt187 /1.875d+00/
      data pt6562 /6.5625d+00/
      data pi32 /5.56832799683170d+00/
      data toln /1.0d-10/
c
c     ----- ishell
c
      do 100 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      snrm=zero
      pnrm=zero
      dnrm=zero
      fnrm=zero
cmw      gnrm=zero
      do 8 i=mini,maxi
      go to (3,4,8,8,5,8,8,8,8,8,6,8,8,8,8,8,8,8,8,8,
     1       7,8,8,8,8,8,8,8,8,8,8,8,8,8,8),i
    3 snrm=one
      go to 8
    4 pnrm=one
      go to 8
    5 dnrm=one
      go to 8
    6 fnrm=one
      go to 8
    7 continue
cmw    gnrm=one
    8 continue
      if(normf.eq.1) go to 25
c
c     ----- calculate normalization factor      -----
c
      snrmi=zero
      pnrmi=zero
      dnrmi=zero
      fnrmi=zero
cmw      gnrmi=zero
      do 20 ig=i1,i2
      do 20 jg=i1,ig
      ee=ex(ig)+ex(jg)
      dum=ee* sqrt(ee)
      dums=cs(ig)*cs(jg)/dum
      dump=pt5*cp(ig)*cp(jg)/(ee*dum)
      dumd=pt75*cd(ig)*cd(jg)/(ee**2*dum)
      dumf=pt187*cf(ig)*cf(jg)/(ee**3*dum)
cmw      dumg=pt6562*cg(ig)*cg(jg)/(ee**4*dum)
      if(jg.eq.ig) go to 10
      dums=dums+dums
      dump=dump+dump
      dumd=dumd+dumd
      dumf=dumf+dumf
cmw      dumg=dumg+dumg
   10 snrmi=snrmi+dums
      pnrmi=pnrmi+dump
      dnrmi=dnrmi+dumd
      fnrmi=fnrmi+dumf
cmw      gnrmi=gnrmi+dumg
   20 continue
      if(snrmi.gt.toln) snrm=one/ sqrt(snrmi*pi32)
      if(pnrmi.gt.toln) pnrm=one/ sqrt(pnrmi*pi32)
      if(dnrmi.gt.toln) dnrm=one/ sqrt(dnrmi*pi32)
      if(fnrmi.gt.toln) fnrm=one/ sqrt(fnrmi*pi32)
cmw      if(gnrmi.gt.toln) gnrm=one/ sqrt(gnrmi*pi32)
c
   25 continue
c
      do 80 i=mini,maxi
      li=loci+i
      go to (30,40,40,40,50,50,50,50,50,50,
     1       60,60,60,60,60,60,60,60,60,60,
     2       70,70,70,70,70,70,70,70,70,70,70,70,70,70,70),i
   30 anorm(li)=snrm
      go to 80
   40 anorm(li)=pnrm
      go to 80
   50 anorm(li)=dnrm
      go to 80
   60 anorm(li)=fnrm
      go to 80
   70 continue
cmw   anorm(li)=gnrm
   80 continue
c
      do 90 ig=i1,i2
      cs(ig)=cs(ig)*snrm
      cp(ig)=cp(ig)*pnrm
      cd(ig)=cd(ig)*dnrm
      cf(ig)=cf(ig)*fnrm
cmw      cg(ig)=cg(ig)*gnrm
   90 continue
c
  100 continue
      return
      end
c==========================================================
      subroutine bfnlab
      implicit REAL  (a-h,o-z),integer  (i-n)
      character*2 lab
      character*2 bfnam,atnam1,atnam2,number
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/mollab)
      dimension number(128),bfnam(35),atnam1(57),atnam2(57)
      dimension lab(4)
      character*8 bflabl
      equivalence (bflabl,lab(1))
      data bfnam  /' s',' x',' y',' z','xx','yy','zz','xy','xz','yz',
     1             'f ','f ','f ','f ','f ','f ','f ','f ','f ','f ',
     2             'g ','g ','g ','g ','g ','g ','g ','g ','g ','g ',
     3             'g ','g ','g ','g ','g '/
      data atnam1 /'  ',' h',' l',' b','  ','  ','  ','  ','  ',' n',
     1             ' n',' m',' a',' s','  ','  ',' c',' a','  ',' c',
     2             ' s',' t','  ',' c',' m',' f',' c',' n',' c',' z',
     3             ' g',' g',' a',' s',' b',' k',' r',' s','  ',' z',
     4             ' n',' m',' t',' r',' r',' p',' a',' c',' i',' s',
     5             ' s',' t','  ',' x',' c',' b','  '/
      data atnam2 /'h ','e ','i ','e ','b ','c ','n ','o ','f ','e ',
     1             'a ','g ','l ','i ','p ','s ','l ','r ','k ','a ',
     2             'c ','i ','v ','r ','n ','e ','o ','i ','u ','n ',
     3             'a ','e ','s ','e ','r ','r ','b ','r ','y ','r ',
     4             'b ','o ','c ','u ','h ','d ','g ','d ','n ','n ',
     5             'b ','e ','i ','e ','s ','a ','  '/
      data number /' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',
     1             '11','12','13','14','15','16','17','18','19','20',
     2             '21','22','23','24','25','26','27','28','29','30',
     3             '31','32','33','34','35','36','37','38','39','40',
     4             '41','42','43','44','45','46','47','47','49','50',
     5             '51','52','53','54','55','56','57','58','59','60',
     6             '61','62','63','64','65','66','67','68','69','70',
     7             '71','72','73','74','75','76','77','78','79','80',
     8             '81','82','83','84','85','86','87','88','89','90',
     9             '91','92','93','94','95','96','97','97','99','00',
     1             '01','02','03','04','05','06','07','08','09','10',
     2             '11','12','13','14','15','16','17','18','19','20',
     3             '21','22','23','24','25','26','27','28'/
c
      n=0
      do 100 ii=1,nshell
      iat=katom(ii)
      j = nint(czan(iat))
      if(j.lt.1) j=0
      if(j.eq.0.or.j.gt.56) j=57
      mini=kmin(ii)
      maxi=kmax(ii)
      do 100 i=mini,maxi
      n=n+1
      lab(1)=number(iat)
      lab(2)=atnam1(  j)
      lab(3)=atnam2(  j)
      lab(4)= bfnam(  i)
      bflab(n)=bflabl
  100 continue
      return
      end
c --------------------------------------------------
      subroutine freerd(line,name,x,nmax,maxblnk)
c------
c   routine reads (atom)name, charge and coordinates in
c   free format from a 80 character line line.
c   separators are commas (always) and
c   in name a row of blanks.ge.maxblnk
c
c   it is prudent to end name by a comma!
c   if name is longer than 16 characters, it is truncated
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
c
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/iofil)
c
      dimension x(nmax)
      character*8 fform,aform
      character*16 name
      character*80 line
      character *1 char,comma,blank
c
      logical alfa,separ
      data comma/','/
      data blank/' '/
      data fform,aform/'(e20.0)','(a16)'/
c
      name = blank
      nblanks=0
c
   10 format(a)
      do 11, i = 1, nmax
        x(i) = 0.0d0
   11 continue
      alfa = .true.
      separ = .false.
      nx = 0
      nchar = 0
      nc1 = 1
c
c	* * * find # significant characters in next block	
c
   20 do 90, nc = nc1, 80
c 1-----
        read(line(nc:nc),'(a)') char
        if (char .eq. blank) then
c   2-----
          if (nchar .ne. 0) then
c     3-----
            if (alfa) then
c       4-----
c        -----  nblanks.ge.maxblnk in alfa is a separator,
c               so count and check	
c
              nblanks = nblanks + 1
              if (nblanks .le. maxblnk) then
c         5-----
                n2 = nc
                nchar = nchar + 1
                if (nchar .gt. 16) then
                  separ = .true.
                endif
c         5-----
              else
c         5-----
                separ =.true.
c         5-----
              endif
c       4-----
            else
c       4-----
c        -----  a blank is a separator in numerical info	
c
              separ=.true.
c       4-----
            endif
c     3-----
          endif
c   2-----
        else
c   2-----
          nblanks=0
c
c    -----  a comma is always a separator	
c
          if (char .eq. comma) then
c     3-----
            separ = .true.
c     3-----
          else
c     3-----
            if (nchar .eq. 0) n1 = nc
            nchar = nchar + 1
            n2 = nc
            if (nc .eq. 80) separ = .true.
c     3-----
          endif
c   2-----
        endif
c
c  -----  we found a separator: process info in this block	
c
        if (separ) then
c   2-----
          if (alfa) then
c     3-----
            if (nchar .gt. 16) nchar = 16
            write(aform(3:4),'(i2)') nchar
            read(line(n1:n2),aform) name
            alfa = .false.
c     3-----
          else
c     3-----
            nx = nx + 1
            if (nx .gt. nmax) return
            if (nchar .ne. 0) then
c       4-----
              write(fform(3:4),'(i2)') nchar
              read(line(n1:n2),fform) x(nx)
c       4-----
            endif
c     3-----
          endif
          nchar = 0
          separ = .false.
c
c    -----  find possible trailing blanks and/or comma
c
          do 30, ntb = nc+1, 80
c     3-----
            read(line(ntb:ntb),'(a)') char
            if (char .ne. blank) then
c       4-----
              if (char .eq. comma) then
c         5-----
                next = ntb + 1
                goto 95
              else
                goto 40
c         5-----
              endif
c       4-----
            endif
c     3-----
   30     continue
   40     next = ntb
          goto 95
c   2-----
        endif
c 1-----
   90 continue
   95 continue
      nc1 = next
      if (next .lt. 80) goto 20
      return
      end
c --------------------------------------
c     subroutine hnderr now handled through gamess_hondo.m
c --------------------------------------
c     subroutine lcpend
c     implicit REAL  (a-h,o-z),integer  (i-n)
c     character*10 fdate,day
c     character*8 hour
cINCLUDE(../m4/common/sizes)
cINCLUDE(comdrf/sizesrf)
cINCLUDE(comdrf/iofil)
cibm      call date(day)
c     call time(hour)
ccyb      hour=time()
ccyb      day=date()
c     day=fdate()
c     write(iwr,9999) hour,day
c     return
c9999 format(/,' job ended at ',a8,' on ',a10)
c     end
c --------------------------------------------------------------------
c
c    gauss-hermite calc
c
c --------------------------------------------------------------------
c     subroutine sxyz now handled through gamess_hondo.m
c     where sxyzh [analg.m] is called.
c ----------------------------------------------------
      subroutine qdpints(xxs,yys,zzs,xys,xzs,yzs,cm,l1)
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical iandj
      logical norm,double
      logical out
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
INCLUDE(comdrf/bas)
INCLUDE(comdrf/ijpair)
cmw      common/output/nprint
cmw      common/iofile/ir,iw
cmw      common/nshel/
cmw     1ex(2048),cs(2048),cp(2048),cd(2048),cf(2048),cg(2048),
cmw     1kstart(512),katom(512),ktype(512),kng(512),kloc(512),kmin(512),
cmw     2kmax(512),nshell
cmw      common/baspar/normf,normp,itol
      common/xyzqdp/xint,yint,zint,xintx,yinty,zintz,
     1              xintxx,yintyy,zintzz,t,xc,yc,zc,
     2              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
cmw      common/ijpair/ia(1)
      dimension xxs(*),yys(*),zzs(*),xys(*),xzs(*),yzs(*)
      dimension sxx(225),syy(225),szz(225),sxy(225),sxz(225),syz(225)
      dimension dij(225)
      dimension   xin(5,5),  yin(5,5),  zin(5,5)
      dimension  xxin(5,5), yyin(5,5), zzin(5,5)
      dimension xxxin(5,5),yyyin(5,5),zzzin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension cm(3)
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      tol=rln10*itol
      out=nprint.eq.6
      norm=normf.ne.1.or.normp.ne.1
c
      xc=cm(1)
      yc=cm(2)
      zc=cm(3)
c
c     ----- ishell -----
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
      sxx(ij)=zero
      syy(ij)=zero
      szz(ij)=zero
      sxy(ij)=zero
      sxz(ij)=zero
      syz(ij)=zero
  100 continue
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
cmw      cgi=cg(ig)
c
c     ----- j primitive -----
c
      jgmax=j2
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
cmw      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      double=iandj.and.ig.ne.jg
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
cmw  180 dum1=cgi*fac
  180 go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      jmax=maxj
      if(iandj) jmax=i
      do 360 j=minj,jmax
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      if(.not.double) go to 350
      if(i.gt.1) go to 240
      dum2=dum2+dum2
      go to 350
  240 dum2=dum2+csi*cpj*fac
      go to 350
  250 dum2=dum1*cpj
      if(double) dum2=dum2+dum2
      go to 350
  260 dum2=dum1*cdj
      if(double) dum2=dum2+dum2
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      if(double) dum2=dum2+dum2
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
cmw  310 dum2=dum1*cgj
  310 if(double) dum2=dum2+dum2
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 dij(ij)=dum2
c
c     ----- quadrupole moment integrals -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljt
      nj=j
      do 370 i=1,lit
      ni=i
      call qdpxyzs
        xin(i,j)=xint*t1
        yin(i,j)=yint*t1
        zin(i,j)=zint*t1
       xxin(i,j)=xintx*t1
       yyin(i,j)=yinty*t1
       zzin(i,j)=zintz*t1
      xxxin(i,j)=xintxx*t1
      yyyin(i,j)=yintyy*t1
      zzzin(i,j)=zintzz*t1
  370 continue
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 390 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      ij=ij+1
      sxx(ij)=sxx(ij)+dij(ij)*(xxxin(ix,jx)*  yin(iy,jy)*  zin(iz,jz))
      syy(ij)=syy(ij)+dij(ij)*(  xin(ix,jx)*yyyin(iy,jy)*  zin(iz,jz))
      szz(ij)=szz(ij)+dij(ij)*(  xin(ix,jx)*  yin(iy,jy)*zzzin(iz,jz))
      sxy(ij)=sxy(ij)+dij(ij)*( xxin(ix,jx)* yyin(iy,jy)*  zin(iz,jz))
      sxz(ij)=sxz(ij)+dij(ij)*( xxin(ix,jx)*  yin(iy,jy)* zzin(iz,jz))
      syz(ij)=syz(ij)+dij(ij)*(  xin(ix,jx)* yyin(iy,jy)* zzin(iz,jz))
  390 continue
 6000 continue
 7000 continue
c
c     ----- set up quadrupole moment matrices -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=ia(loci+i)+(locj+j)
      xxs(nn)=sxx(ij)
      yys(nn)=syy(ij)
      zzs(nn)=szz(ij)
      xys(nn)=sxy(ij)
      xzs(nn)=sxz(ij)
      yzs(nn)=syz(ij)
 7500 continue
c
 8000 continue
 9000 continue
c
      if(.not.out) go to 9100
      write(iwr,9999)
      call prtrl(xxs,l1)
      write(iwr,9998)
      call prtrl(yys,l1)
      write(iwr,9997)
      call prtrl(zzs,l1)
      write(iwr,9996)
      call prtrl(xys,l1)
      write(iwr,9995)
      call prtrl(xzs,l1)
      write(iwr,9994)
      call prtrl(yzs,l1)
 9100 continue
c
      return
 9999 format(/,10x,26(1h-),/,10x,'xx-second moment integrals',/,
     1         10x,26(1h-))
 9998 format(/,10x,26(1h-),/,10x,'yy-second moment integrals',/,
     1         10x,26(1h-))
 9997 format(/,10x,26(1h-),/,10x,'zz-second moment integrals',/,
     1         10x,26(1h-))
 9996 format(/,10x,26(1h-),/,10x,'xy-second moment integrals',/,
     1         10x,26(1h-))
 9995 format(/,10x,26(1h-),/,10x,'xz-second moment integrals',/,
     1         10x,26(1h-))
 9994 format(/,10x,26(1h-),/,10x,'yz-second moment integrals',/,
     1         10x,26(1h-))
      end
      subroutine qdpxyzs
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      common/xyzqdp/xint,yint,zint,xintx,yinty,zintz,
     1 xintxx,yintyy,zintzz,t,xc,yc,zc,x0,y0,z0,
     2 xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h1,h2(2),h3(3),h4(4),h5(5),h6(6),h7(7)
      common/wermit/w1,w2(2),w3(3),w4(4),w5(5),w6(6),w7(7)
      dimension h(28),w(28),min(7),max(7)
      equivalence (h(1),h1),(w(1),w1)
      data min /1,2,4, 7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      xintx=zero
      yinty=zero
      zintz=zero
      xintxx=zero
      yintyy=zero
      zintzz=zero
      npts=(ni+nj-2+2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
      xintx=xintx+px*(ptx-xc)
      yinty=yinty+py*(pty-yc)
      zintz=zintz+pz*(ptz-zc)
      xintxx=xintxx+px*(ptx-xc)*(ptx-xc)
      yintyy=yintyy+py*(pty-yc)*(pty-yc)
      zintzz=zintzz+pz*(ptz-zc)*(ptz-zc)
   16 continue
      return
      end
c ------------------------------------------from scfnew.f
      subroutine scfinp
      implicit REAL  (a-h,o-z),integer  (i-n)
      integer gvbflg,uhfflg
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/dafil)
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/bas)
INCLUDE(../m4/common/infoa)
      common/sqfile/ijk,ipk
      common/scfcvg/tolscf,itermx,micro,npnch,nselct
      common/scfprt/imoprt
      common/scfout/mprint
      common/dmping/oshft,vshft,dmp,dmp0,ixtrp
caleko 
c 
INCLUDE(comdrf/scfopt)
caleko
           
c
INCLUDE(comdrf/runpar)
c
      dimension cicoef(2,12),f(25),alpha(325),beta(325)
      dimension no(10)
      character *8 errmsg(3)

caleko
      dimension ci(2,12),fii(25)
INCLUDE(comdrf/alpbet)
caleko
c
      logical singlet
c
      character*8 rhf, uhf, gvb
c
      namelist /scf/ damp,damp0,iextrp,vshift,oshift,
     1 maxit,micit,npunch,npflg,moordr,acurcy,uhfflg,gvbflg,
     2 nco,nseto,npair,no,cicoef,f,alpha,beta,nspin,moprt,nchar,
     3 noavg,nopel
c
      namelist /scf2/ damp,damp0,iextrp,vshift,oshift,
     1 maxit,micit,npunch,npflg,moordr,acurcy,uhfflg,gvbflg,
     2 nco,nseto,npair,no,cicoef,f,alpha,beta,nspin,moprt,nchar,
     3 noavg,nopel
c
      namelist /scf3/ damp,damp0,iextrp,vshift,oshift,
     1 maxit,micit,npunch,npflg,moordr,acurcy,uhfflg,gvbflg,
     2 nco,nseto,npair,no,cicoef,f,alpha,beta,nspin,moprt,nchar,
     3 noavg,nopel
c
      namelist /scf4/ damp,damp0,iextrp,vshift,oshift,
     1 maxit,micit,npunch,npflg,moordr,acurcy,uhfflg,gvbflg,
     2 nco,nseto,npair,no,cicoef,f,alpha,beta,nspin,moprt,nchar,
     3 noavg,nopel
c
      namelist /scf5/ damp,damp0,iextrp,vshift,oshift,
     1 maxit,micit,npunch,npflg,moordr,acurcy,uhfflg,gvbflg,
     2 nco,nseto,npair,no,cicoef,f,alpha,beta,nspin,moprt,nchar,
     3 noavg,nopel
c
      data errmsg        /'program ','stop in ','-scfinp-'/
      data rhf,uhf,gvb   /'rhf     ','uhf     ','gvb     '/
      data uhfflg,gvbflg /0,0/
      data singlet /.false./
c
      data damp,damp0    /1.0d+00,0.0d+00/
      data vshift,oshift /0.0d+00,0.0d+00/
      data maxit         /30/
      data micit         /15/
      data iextrp        /5/
      data moordr        /0/
      data npunch        /0/
      data npflg         /0/
      data acurcy        /1.0d-05/
      data nspin         /1/
      data nchar         /0/
      data nco           /-1/
      data nseto         /0/
      data no            /10*0/
      data npair         /0/
      data mxseto        /10/
      data mxnham        /25/
      data mxpair        /12/
c
      data noavg         /1/
      data nopel         /0/
c
      data f      / 1.0d+00, 0.5d+00, 0.5d+00,
     1          22*0.0d+00/
      data alpha  / 2.0d+00, 1.0d+00, 0.5d+00, 1.0d+00,0.5d+00, 0.5d+00,
     1         319*0.0d+00/
      data beta   /-1.0d+00,-0.5d+00,-0.5d+00,-0.5d+00,0.5d+00,-0.5d+00,
     1         319*0.0d+00/
      data cicoef /9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,
     1             9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,
     2             9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,
     3             9.0d-01,-2.0d-01,9.0d-01,-2.0d-01,9.0d-01,-2.0d-01/
      data zero, pt5, one, two    /0.d00,0.5d+00,1.d00,2.d00/
c      ----- read namelist $scf -----
c
c     maxit  = maximum number of iterations.
c
c     npunch = punch out option (used in conjunction with -nprint- )
c
c     npflg  = print flag ( also in conjunction with -nprint- )
c
c     moprt  = number of orbitals printed at the end of scf
c
c     moordr = flag for selection method of the new mo's
c              during scf.
c
c     acurcy = convergence criterion on density matrix elements.
c
c     damp = damping factor for density matrix
c            d(+1) = ( d( 0) + damp * d(-1) )/( 1. +damp )
c
c     damp0 = cutoff parameter for extrapolation while
c             damping. the value of cutoff is @
c             cutoff = 0.5 + damp0
c             damp0 should have a value between -.5 and +.5
c             the larger cutoff, the more extrapolation
c
c     iextrp = interval frequency for extrapolation while
c              damping.
c
c     vshift = level shifting parameter for virtual orbitals.
c
c              the level shifting technique takes precedence
c              over the damping technique, i.e. if vshift is
c              not equal to zero, then -damp- is set equal to zero.
c
c     oshift = level shifting parameter for open shells.
c
c     micit  = maximum number of micro iterations.
c
c     uhfflg = flag for uhf wavefunction
c
c     gvbflg = flag for gvb wavefunction
c
c     nco    = number of closed shells.
c
c     nseto  = number of open shell classes.
c
c     npair  = number of gvb pairs.
c
c     no     = number of open shells in each class.
c
c     f     = f term for gvb energy expression
c
c     alpha = alpha coupling coefficients
c
c     beta   = beta coupling coefficients
c
c     nspin  = spin multiplicity.
c
c     nchar  = total charge
c
c
c     noavg  = 1: do not calculate average state orbitals (default)
c              0: calculate average state orbitals
c              note: forces gvbflg = 1
c              and expects nopel to be given
c
c     nopel  = number of electrons in open shells for calcualtion
c              of average state orbitals
c
      moprt=num
      rewind ir
c
c-----  read scf input for actual state to be calculated
c
      goto (81,82,83,84,85) iactst
   81   read(ir,scf,end=5)
        goto 5
   82   read(ir,scf2,end=5)
        goto 5
   83   read(ir,scf3,end=5)
        goto 5
   84   read(ir,scf4,end=5)
        goto 5
   85   read(ir,scf5,end=5)
        goto 5
c
    5 continue
      if(moprt.lt.0.or.moprt.gt.num) moprt=num
      imoprt=moprt
      npnch =npunch
      mprint=npflg
      scftyp=rhf
c
      if(uhfflg.ne.0) scftyp=uhf
      if(gvbflg.ne.0.or.npair.ne.0) scftyp=gvb
      if(scftyp.eq.rhf.and.nseto.gt.1) scftyp=gvb
c
      if (noavg .eq. 0) then
        if (nopel .ne. 0) then
          write(iwr,1001)
 1001     format
     1 (/,'  note: average state calculation: via gvb option',/,
     2    '        ignore value of spin multiplicity')
          scftyp = gvb
          nspin = 1
        else
          write(iwr,1011)
 1011     format
     1 (/,'  number of open shell electrons not given in $scf',
     2 ' while noavg = 0; check input')
          call hnderr(3,errmsg)
        endif
      endif
c
      if(scftyp.eq.gvb) gvbflg=1
      if(scftyp.ne.gvb) npair=0
      nselct=moordr
      tolscf=acurcy
      itermx=maxit
      oshft=oshift
      vshft=vshift
      dmp  =damp
      dmp0 =damp0
      ixtrp=iextrp
      micro=micit
c
c     ----- set wavefunction data -----
c
      mul=nspin
      na=(ne+(nspin-1)-nchar)/2
      nb=(ne-(nspin-1)-nchar)/2
      ich=ne-(na+nb)
c
      if(scftyp.eq.uhf) goto 100
c
      singlet = nspin .eq. 1
c
      if (nco .lt. 0) then
c 1-----
        ncorb=nb
c
c  -----  skip this business for average state calculation
c
        if ((noavg .eq. 1) .and. (na .eq. nb)) then
c   2-----
c * * *  singlet: rhf or gvb (default states for two open shell electron
c
          if ((nseto.ne.0) .or.(npair.ne.0)) then
c     3-----
c * * *  singlet rohf-gvb for two shells only
c
            singlet = .true.
            if(nseto.gt.1) nseto=1
            ncorb=ncorb-nseto-npair
            if(nseto.ne.0) then
c       4-----
              nseto=2
              no(1)=1
              no(2)=1
c       4-----
            endif
c     3-----
          endif
c   2-----
        else if ((noavg .eq. 1) .and. (na .ne. nb)) then
c   2-----
c * * *  highest multiplicity
c    -----  note!!!!!!!!!!!!!!
c           for open-shell rhf, a closed shell is always
c           supposed to be there. provisions are made to define
c           one even if nco = 1.
c           this is not so for gvb, so provisions are made here.
c
          singlet = .false.
c
          if(scftyp.eq.rhf) then
c     3-----
            nseto=1
            no(1)=na-nb
c     3-----
          else
c     3-----
c      -----  for gvb functions, make the maximum multiplicity
c             more general
c
            nseto=na - nb
            do 400, i = 1, nseto
              no(i)=1
  400       continue
c     3-----
          endif
c   2-----
        endif
c
c  -----  set number of closed shells
c
        nco=ncorb
c 1-----
      endif
c
c-----  set proper coupling constants for gvb wf
c
      if (gvbflg .eq. 1) then
c 1-----
        if (noavg .eq. 0) then
c   2-----
          call gvbpar(ne,nco,nopel,nseto,no,f,alpha,beta)
c   2-----
        else if ((nco .eq. 0) .and. (f(1) .eq. one)) then
c   2-----
c    -----  re-initialise occupation and couplings for
c           a gvb wavefunction with nco = 1,
c           unless the user has shown by fractionally
c           occupying the first shell that he knows there are
c           no closed shells.
c
          ij = 0
          do 500, i = 1, nseto
            f(i) = pt5
            do 410, j = 1, i
              ij = ij + 1
              alpha(ij) = pt5
              beta(ij) = -pt5
              if (singlet .and. (i .ne. j)) beta(ij) = pt5
  410       continue
  500     continue
c   2-----
        endif
c 1-----
      endif
c
      if(nseto.le.mxseto) go to 10
      write(iwr,9999)
      call hnderr(3,errmsg)
   10 continue
      nham=nseto+2*npair
      if(nco.gt.0) nham=nham+1
      if(nham.le.mxnham) go to 15
      write(iwr,9997)
      call hnderr(3,errmsg)
   15 continue
      if(npair.le.mxpair) go to 20
      write(iwr,9996)
      call hnderr(3,errmsg)
   20 continue
      ncorb=nco
      nopset=nseto
      npairs=npair
      do 25 i=1,nseto
   25 nopen(i)=no(i)
c
c     ----- define coupling coefficients -----
c
      nham=nseto
      if(nco.gt.0) nham=nham+1
      ij=0
      do 30 iham=1,nham
      fii(iham)=f(iham)
      do 30 jham=1,iham
      ij=ij+1
      alphij(ij)=alpha(ij)
      betaij(ij)= beta(ij)
   30 continue
      if(gvbflg.ne.0.and.nspin.eq.3.and.nseto.eq.2
     1 .and. noavg .eq. 1 .and. (nco .gt. 0)) betaij(5)=-pt5
c
c     ----- more wavefunction data -----
c
      nb=nco
      na=nco
      if(scftyp.ne.gvb.or.npair.eq.0) go to 70
      do 60 ipair=1,npair
      dum=cicoef(1,ipair)**2+cicoef(2,ipair)**2
      dum= sqrt(dum)
      cicoef(1,ipair)=cicoef(1,ipair)/dum
      cicoef(2,ipair)=cicoef(2,ipair)/dum
      ci(1,ipair)=cicoef(1,ipair)
      ci(2,ipair)=cicoef(2,ipair)
   60 continue
   70 continue
c
      na=nco
      if(scftyp.eq.gvb) na=na+2*npair
      if(nseto.eq.0) go to 90
      do 80 iseto=1,nseto
   80 na=na+no(iseto)
   90 continue
      nb=na
c
  100 continue
c
      return
 9999 format(' too many sets of open shells (max=10). stop.')
 9998 format(' no namelist $scf found. stop.')
 9997 format(47h too many hamiltonian operators (max=25). stop.)
 9996 format(39h too many orbital pairs (max=12). stop.)
      end
      subroutine gvbpar(ne,nco,nopel,nseto,no,f,alpha,beta)
c------
c      defines the coupling parameters for an average state calculation
c      see: motecc, hondo (chapter 6)
c           mcweeny
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
c
      dimension no(10)
      dimension f(25), alpha(325), beta(325)
c
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ijpair)
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-gvbpar-'/
c
      data zero, one, two /0.0d00, 1.0d00, 2.0d00/
      data pt5 /0.5d00/
c
c-----  count number of open shell orbitals
c
      noporb = 0
      nset = nseto
      do 100, i = 1, nset
c 1-----
c  -----  split set of open shells if more than one open shell is given
c         this is necessary since there seem to be problems with having
c         more than one open shell orbital in the same set (see hondo
c         manual)
c
        if (no(i) .gt. 1) then
c   2-----
          do 90, k = 2, no(i)
c     3-----
            nseto = nseto + 1
            if (nseto .gt. 10) then
              write(iwr,1001)
 1001         format
     1 (/,'  too many open shells have been defined, sorry!')
              call hnderr(3,errmsg)
            endif
c
            no(nseto) = 1
            noporb = noporb + 1
c     3-----
   90     continue
c
          no(nset) = 1
c   2-----
        endif
c
        noporb = noporb + 1
c 1-----
  100 continue
c
c-----  calculate occupation of the shells and coupling constants
c       between shells
c
c-----  closed shell
c
      nco = nco - nopel/2
c
      ncfs = 0
c
      if (nco .gt. 0) then
        ncfs = ncfs + 1
        f(ncfs) = one
        alpha(ncfs) = two
        beta (ncfs) = -one
      endif
c
      elnum = dble(nopel)
      orbnum = dble(noporb)
      occfac = elnum/(two*orbnum)
      coupfac = (elnum-one)/((two*orbnum)-one)
c
      do 200, i = (ncfs + 1), (ncfs + nseto)
c 1-----
        f(i) = occfac
c
c  -----  closed-open
c
        if (nco .gt. 0) then
          alpha(ia(i) + 1) = two*occfac
          beta (ia(i) + 1) = -occfac
        endif
c
c  -----  open-open
c
        do 190, k = (ncfs+1), i
c   2-----
          alpha(ia(i) + k) = two*occfac*coupfac
          beta (ia(i) + k) = -occfac*coupfac
c   2-----
  190   continue
c 1-----
  200 continue
c
      return
      end
c      function dottri(a,b,n)
c      implicit REAL  (a-h,o-z),integer  (i-n)
cc
cc     ----- trace of product of 2 symmetric matrices -----
cc           -a- and -b- stored linearly
cc
c      dimension a(*),b(*)
c      data zero /0.0d+00/
c      dum=zero
c      k=0
c      do 20 i=1,n
c      do 10 j=1,i
c      k=k+1
c      ab=a(k)*b(k)
c   10 dum=dum+(ab+ab)
c   20 dum=dum- ab
c      dottri=dum
c      return
c      end
c ----------------------------------------------------------from ctl.f
      subroutine tftr(h,f,q,t,ia,m,n,ndim)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      dimension h(*),f(*),q(ndim,*),t(*),ia(*)
      data small /1.0d-11/
      data zero  /0.0d+00/
      ij = 0
      do 180 j = 1,m
      ik = 0
      do 140 i = 1,n
      dum = zero
      qij = q(i,j)
      max = i-1
      if (max .eq. 0) go to 120
      do 100 k = 1,max
      ik = ik+1
      t(k) = t(k)+ f(ik)*qij
      dum  = dum + f(ik)*q(k,j)
  100 continue
  120 ik = ik+1
      t(i) = dum + f(ik)*qij
  140 continue
      do 160 i = 1,j
      ij = ij+1
      hij = ddot(n,q(1,i),1,t,1)
      if( abs(hij).lt.small) hij=zero
      h(ij)=hij
  160 continue
  180 continue
      return
      end
c -----------------------------------from ctl.f
      subroutine prtrl(d,n)
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c     ----- print out a triangular matrix with labels
c
INCLUDE(comdrf/dafil)
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/mollab)
      dimension d(*),dd(10)
      if(list.eq.0) mmax=10
      if(list.eq.1) mmax=7
      if(list.eq.2) mmax=7
      imax = 0
  100 imin = imax+1
      imax = imax+mmax
      if (imax .gt. n) imax = n
      write (iwr,9008)
      if(list.eq.0) write (iwr,9028) (bflab(i),i = imin,imax)
      if(list.eq.1) write (iwr,9128) (bflab(i),i = imin,imax)
      if(list.eq.2) write (iwr,9228) (bflab(i),i = imin,imax)
      write (iwr,9008)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      ii = max( i, j)
      jj = min( i, j)
      ij = (ii*(ii-1))/2 + jj
  140 dd(k) = d(ij)
      if(list.eq.0) write (iwr,9048) j,bflab(j),(dd(i),i = 1,k)
      if(list.eq.1) write (iwr,9148) j,bflab(j),(dd(i),i = 1,k)
      if(list.eq.2) write (iwr,9248) j,bflab(j),(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(15x,10(2x,a8,1x))
 9048 format(i5,2x,a8,10f11.5)
 9128 format(15x,7(4x,a8,3x))
 9148 format(i5,2x,a8,7f15.10)
 9228 format(15x,7(4x,a8,3x))
 9248 format(i5,2x,a8,7e15.8)
      end
