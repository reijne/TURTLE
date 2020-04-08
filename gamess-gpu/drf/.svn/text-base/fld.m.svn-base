cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hnd8 : fld fortran
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***********************************************************************
c*                                                                     *
c*  copyright i.b.m. 1989,1990                                         *
c*  all rights reserved.                                               *
c*                                                                     *
c*  permission is hereby granted to use but not to reproduce or        *
c*  distribute this material. the use is restricted to research        *
c*  purposes only. this material may not be reproduced for             *
c*  commercial purposes or included in a commercial product            *
c*  without the specific written permission of i.b.m. .                *
c*                                                                     *
c*  any publication based upon results obtained from this material     *
c*  must include the following attribution:                            *
c*                                                                     *
c*  this work is based, in part, on results from the motecc(tm)        *
c*  package.  motecc is a trademark of i.b.m. .                        *
c*                                                                     *
c*                                                                     *
c***********************************************************************
      subroutine fldinp
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      common/iofile/ir,iwr,ip
cmw      common/fldpar/fieldx,fieldy,fieldz,ifield
      common/symtry/invt(48),nt,ntmax,ntwd,nosym
cmw      common/runopt/runtyp
INCLUDE(comdrf/opt)
      character*4 keyblk,keyfld
      character*8 ramir
c
      namelist /fld/ fldx,fldy,fldz
c
cmw      data keyblk /4h    /
      data keyblk /'    '/
cmw      data keyfld /4h fld/
      data keyfld /' fld'/
cmw      data ramir  /8hir-raman/
      data ramir /'ir-raman'/
c
      data fldx,fldy,fldz /0.0d+00,0.0d+00,0.0d+00/
      data zero   /0.0d+00/
c
c     ----- read namelist -$fld- -----
c
c     ifield=keyblk
c
c     ----- if -irramx- job , force no field -----
c
      if(runtyp.eq.ramir) go to 10
c
      rewind ir
      read(ir,fld,end=10)
c     ifield=keyfld
      write(iwr,9999) fldx,fldy,fldz
   10 continue
      fieldx=fldx
      fieldy=fldy
      fieldz=fldz
c
c     ----- print warning about symmetry -----
c
      if((fieldx.ne.zero.or.fieldy.ne.zero.or.fieldz.ne.zero).and.
     1   (nosym.eq.0)) write(iwr,9998)
      return
c
 9999 format(/,10x,29(1h-),/,10x,'dipole field strengths (a.u.)',
     1       /,10x,29(1h-),
     1       /,' fldx = ',f10.6,' fldy = ',f10.6,' fldz = ',f10.6)
 9998 format(/,1x,81(1h.),
     1       /,1x,'warning : an external electric field reduces the ',
     1            'symmetry of the electron density.',
     1       /,1x,'you may have to lower or turn off symmetry in ',
     1            '$cntrl- .'
     1       /,1x,81(1h.))
      end
_IF()
      subroutine fldint(h0,xfld,yfld,zfld)
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      logical some
      common/iofile/ir,iwr,ip
      common/sqfile/ijk,ipk
cmw      common/dafile/idaf,nav,ioda(255)
INCLUDE(comdrf/dafil)
INCLUDE(../m4/common/infoa)
cmw      common/fldpar/fldx,fldy,fldz,ifld
         common/restar/nprint
INCLUDE(comdrf/opt)
      common /c_of_m/pcm,qcm,rcm
      dimension h0(*),xfld(*),yfld(*),zfld(*)
      dimension cm(3)
      data zero /0.0d+00/
c
      some=nprint.ne.-5
      if(some) write(iwr,9999)
      do 10 i=1,3
   10 cm(i)=zero
      l1= num
      l2=(num*(num+1))/2
c
c     ----- calculate dipole moment integrals -----
c
cmw put the centre-of-mass information into the common
c   and call the gamess dipint
c   note that dipint was adapted for cm().
c   a common c_of_m is used to pass cm()
      pcm=cm(1)
      qcm=cm(2)
      rcm=cm(3)
cmw      call dipint(xfld,yfld,zfld,cm,l1)
cahv  call dipxyz(dipx,dipy,dipz,num)
      call dipxyz(dipx,dipy,dipz)
c
c     ----- add -hfld- to -h0- -----
c
      call daread(idafh,ioda,h0,l2,11)
      do 20 i=1,l2
   20 h0(i)=h0(i)-(-xfld(i)*fldx-yfld(i)*fldy-zfld(i)*fldz)
      call dawrit(idafh,ioda,h0,l2,11,nav)
c
      if(some) write(iwr,9998)
 9999 format(/,10x,15(1h-),/,10x,'-fld- integrals',/,10x,15(1h-))
 9998 format(' ...... end of -fld- integrals ......')
      return
      end
      subroutine fldder
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
cmw      common/memory/maxcor,maxlcm
INCLUDE(comdrf/mem)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/scfopt)
INCLUDE(../m4/common/infoa)
      character *8 scf, uhf
      common/scm/x(1)
      character *8 errmsg(3)
      data errmsg  /'program ','stop in ','-fldder-'/
      data scf,uhf /'scf     ','uhf     '/
c
c     ----- get core memory -----
c
      call cmem(loadcm)
      locx1=loccm(x(1))
      ngotmx=maxcor-locx1
c
      i10=1
      i20=i10+(num*(num+1))/2
      i30=i20+(num*(num+1))/2
c
      last=i20+1
      if(wfntyp.eq.scf.and.scftyp.eq.uhf) last=i30
      last=last-1
      if(last.le.ngotmx) go to 10
      need=locx1+last
      call needcm(last,need)
      call hnderr(3,errmsg)
   10 continue
c
c     ----- get density matrix -----
c
cnot      call denddm(x(i10),x(i20),num)
c
c     ----- get gradient of field interaction -----
c
      call flddin(x(i10))
c
c     ----- reset core memory -----
c
      call cmem(ngotcm)
      if(ngotcm.ne.loadcm) call setc(loadcm)
      return
      end
      subroutine flddin(dab)
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      logical iandj,out,some,norm
      logical iskip
      common/times/timlim,ti,tx,tim,ti0
INCLUDE(comdrf/opt)
c
c     common/rys/xx,u(5),w(5),nroots
c
      common/iofile/ir,iwr,ip
      common/sqfile/ijk,ipk
cmw      common/dafile/idaf,nav,ioda(255)
INCLUDE(comdrf/dafil)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
cmw      common/baspar/normf,normp,itol
INCLUDE(comdrf/bas)
cmw      common/fldpar/fldx,fldy,fldz,ifld
      common/xyzdip/xint0,yint0,zint0,xintx,yinty,zintz,t,xc,yc,zc,
     1              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
cmw      common/ijpair/ia(1)
INCLUDE(comdrf/ijpair)
      common/grad12/de(3,128)
      dimension dij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension xs(6,5),ys(6,5),zs(6,5)
      dimension xx(6,5),yy(6,5),zz(6,5)
      dimension dxsdi(5,5),dysdi(5,5),dzsdi(5,5)
      dimension dxxdi(5,5),dyydi(5,5),dzzdi(5,5)
      dimension dnam(3)
      dimension dab(*)
      data dnam   /'e*x','e*y','e*z'/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
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
      tol =rln10*itol
      out =nprint.eq.-3
      some=nprint.ne.-5
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- calculate field derivatives -----
c
      if(out) write(iwr,9999)
      nder=1
c
      xc=zero
      yc=zero
      zc=zero
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
      lit=ktype(ii)+1
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
      litder=lit+nder
      iat=i
c
c     ----- jshell -----
c
      do 8000 jj=1,nshell
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
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      dum=ai+ai
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
cmw      cgi=cg(ig)
c
c     ----- j primitive -----
c
      do 6000 jg=j1,j2
      aj=ex(jg)
      aa =ai+aj
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
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
cmw  310 dum2=dum1*cgj
  310 go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      nn=ia(max(loci+i,locj+j))
     1  +   min(loci+i,locj+j)
      den=dab(nn)
      den=den+den
      ij=ij+1
  360 dij(ij)=dum2*den
c
c     ----- -fld- derivatives -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljt
      nj=j
      do 370 i=1,litder
      ni=i
      call dipxyz2
      xs(i,j)=xint0*t1
      ys(i,j)=yint0*t1
      zs(i,j)=zint0*t1
      xx(i,j)=xintx*t1
      yy(i,j)=yinty*t1
      zz(i,j)=zintz*t1
  370 continue
c
cnot      call deri(dxsdi,dysdi,dzsdi,xs,ys,zs,lit,ljt,ai)
cnot      call deri(dxxdi,dyydi,dzzdi,xx,yy,zz,lit,ljt,ai)
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      do 380 j=minj,maxj
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      dflxdx=dxxdi(ix,jx)*   ys(iy,jy)*   zs(iz,jz)
      dflydx=dxsdi(ix,jx)*   yy(iy,jy)*   zs(iz,jz)
      dflzdx=dxsdi(ix,jx)*   ys(iy,jy)*   zz(iz,jz)
      dflxdy=   xx(ix,jx)*dysdi(iy,jy)*   zs(iz,jz)
      dflydy=   xs(ix,jx)*dyydi(iy,jy)*   zs(iz,jz)
      dflzdy=   xs(ix,jx)*dysdi(iy,jy)*   zz(iz,jz)
      dflxdz=   xx(ix,jx)*   ys(iy,jy)*dzsdi(iz,jz)
      dflydz=   xs(ix,jx)*   yy(iy,jy)*dzsdi(iz,jz)
      dflzdz=   xs(ix,jx)*   ys(iy,jy)*dzzdi(iz,jz)
      ij=ij+1
      de(1,iat)=de(1,iat)-dij(ij)*(-dflxdx*fldx-dflydx*fldy-dflzdx*fldz)
      de(2,iat)=de(2,iat)-dij(ij)*(-dflxdy*fldx-dflydy*fldy-dflzdy*fldz)
      de(3,iat)=de(3,iat)-dij(ij)*(-dflxdz*fldx-dflydz*fldy-dflzdz*fldz)
  380 continue
  390 continue
c
c     ----- end of primitive loops -----
c
 6000 continue
 7000 continue
c
      if(.not.out) go to 8000
      write(iwr,9994) ii,jj
      mmax=0
 7500 min=mmax+1
      mmax=mmax+8
      if(mmax.gt.nat) mmax=nat
      write(iwr,9998)
      write(iwr,9997) (i,i=min,mmax)
      write(iwr,9998)
      do 7600 i=1,3
 7600 write(iwr,9996) dnam(i),(de(i,j),j=min,mmax)
      if(mmax.lt.nat) go to 7500
c
 8000 continue
 9000 continue
c
      if(.not.out) go to 9300
      mmax=0
 9100 min=mmax+1
      mmax=mmax+8
      if(mmax.gt.nat) mmax=nat
      write(iwr,9998)
      write(iwr,9997) (i,i=min,mmax)
      write(iwr,9998)
      do 9200 i=1,3
 9200 write(iwr,9993) dnam(i),(de(i,j),j=min,mmax)
      if(mmax.lt.nat) go to 9100
      write(iwr,9995)
c
 9300 continue
      return
 9999 format(/,10x,31(1h-),/,10x,' -fld- contribution to gradient',
     1       /,10x,31(1h-))
 9998 format(/)
 9997 format(5x,'atom',8(6x,i2,7x))
 9996 format(7x,a3,8e15.7)
 9995 format(/,' ...... end of -fld- gradient ...... ',/)
 9994 format(' shells ii, jj ',2i5)
 9993 format(7x,a3,8f15.7)
      end
      subroutine efcinp
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      common/iofile/ir,iwr,ip
cmw      common/efcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
caleko
      character*8 efclab
      common/hefcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
      character *4 iefc
caleko
cmw      common/runopt/runtyp
INCLUDE(comdrf/opt)
      common/scm/bond(1000),alpha(1000),beta(1000),sign(1000),
     1 icnox(1000),iefcon(3,1000),inefc(1000)
      dimension title(10)
      character *8 errmsg(3)
      character*8 word,efcnam
      character*8 blank,endwrd,efcwrd,ramir
      character*4 keyblk,keyefc
c
      data errmsg /'program ','stop in ','-efcinp-'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data maxefc /1000/
      data blank  /'        '/
      data endwrd /' $end   '/
      data efcwrd /' $efc   '/
      data ramir  /'ir-raman'/
      data keyblk /'    '/
      data keyefc /' efc'/
c
      iefc=keyblk
      nefc=0
      do 10 k=1,maxefc
      efclab(k) = blank
      efcz(k)=zero
      do 10 l=1,3
   10 efcc(l,k)=zero
c
c     ----- if -irramx- job , force no point charges -----
c
      if(runtyp.eq.ramir) go to 100
c
c     ----- search for  fractional charges -----
c
      rewind ir
   20 read(ir,9999,end=100) word
      if(word.ne.efcwrd) go to 20
      iefc=keyefc
c
      read(ir,9999) title
      write(iwr,9998) title
      read(ir,9997) iunit
      nefcin=0
   30 continue
      nefc1=nefc+1
      call efcrd(efcnam,zefc,x,y,z,efcc,nefc,maxefc,iunit,
     1 bond,alpha,beta,sign,iconx,iefcon,inefc,nefcin)
      if(efcnam.eq.endwrd) go to 100
      nefc=nefc+1
      if(nefc.gt.maxefc) go to 110
      efclab(nefc)=efcnam
      efcz(nefc)=zefc
      efcc(1,nefc)=x
      efcc(2,nefc)=y
      efcc(3,nefc)=z
      write(iwr,9996) nefc,efclab(nefc),efcz(nefc),
     1 efcc(1,nefc),efcc(2,nefc),efcc(3,nefc)
      go to 30
c
  100 continue
      if(nefc.eq.0) iefc=keyblk
c      if(nefc.eq.0) write(iwr,9994)
      return
c
  110 write(iwr,9995) maxefc
      call hnderr(3,errmsg)
      return
 9999 format(10a8)
 9998 format(/,10x,13(1h-),/,10x,'point charges',/,10x,13(1h-),
     1       /,1x,10a8,/)
 9997 format(i5)
 9996 format(i4,2x,a8,' efc= ',f9.3,' x= ',f10.5,' y= ',f10.5,
     1                                           ' z= ',f10.5)
 9995 format(' too many point charges. maxefc = ',i5)
 9994 format(' no point charges .')
      end
      subroutine efcrd(atom,znuc,x,y,z,c,nat,natmax,iunit,
     1 bond,alpha,beta,sign,iconx,iatcon,inatom,natin)
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      common/iofile/ir,iwr,ip
      dimension bond(*),alpha(*),beta(*),sign(*),iconx(*),
     1 iatcon(3,*),inatom(*)
      dimension c(3,*)
      dimension cc(3,1004),xyz(3),xyz0(3),t(3,3)
      dimension nconx(7),npt(1005),ipt(3)
      dimension npt001(100),npt101(100),npt201(100),
     1          npt301(100),npt401(100),npt501(100),
     1          npt601(100),npt701(100),npt801(100),
     1          npt901(100),np1001(  5)
      character *8 errmsg(3)
      equivalence (npt(  1),npt001(1)),(npt(101),npt101(1)),
     1            (npt(201),npt201(1)),(npt(301),npt301(1)),
     1            (npt(401),npt401(1)),(npt(501),npt501(1)),
     1            (npt(601),npt601(1)),(npt(701),npt701(1)),
     1            (npt(801),npt801(1)),(npt(901),npt901(1)),
     1            (npt(1001),np1001(1))
      data errmsg /'program ','stop in ','- efcrd-'/
      data endwrd /' $end   '/
      data pi2    /6.28318530717958d+00/
      data degree /360.0d+00/
      data unit   /0.529177249d+00/
      data nconx  /'  lc',' pcc','npcc','ccpa',' ptc',' tct','    '/
      data npt001 /
     1'   1','   2','   3','   4','   5','   6','   7','   8','   9',
     1'  10','  11','  12','  13','  14','  15','  16','  17','  18',
     1'  19','  20','  21','  22','  23','  24','  25','  26','  27',
     1'  28','  29','  30','  31','  32','  33','  34','  35','  36',
     1'  37','  38','  39','  40','  41','  42','  43','  44','  45',
     1'  46','  47','  48','  49','  50',
     1'  51','  52','  53','  54','  55','  56','  57','  58','  59',
     1'  60','  61','  62','  63','  64','  65','  66','  67','  68',
     1'  69','  70','  71','  72','  73','  74','  75','  76','  77',
     1'  78','  79','  80','  81','  82','  83','  84','  85','  86',
     1'  87','  88','  89','  90','  91','  92','  93','  94','  95',
     1'  96','  97','  98','  99',' 100'/
      data npt101 /
     1' 101',' 102',' 103',' 104',' 105',' 106',' 107',' 108',' 109',
     1' 110',' 111',' 112',' 113',' 114',' 115',' 116',' 117',' 118',
     1' 119',' 120',' 121',' 122',' 123',' 124',' 125',' 126',' 127',
     1' 128',' 129',' 130',' 131',' 132',' 133',' 134',' 135',' 136',
     1' 137',' 138',' 139',' 140',' 141',' 142',' 143',' 144',' 145',
     1' 146',' 147',' 148',' 149',' 150',
     1' 151',' 152',' 153',' 154',' 155',' 156',' 157',' 158',' 159',
     1' 160',' 161',' 162',' 163',' 164',' 165',' 166',' 167',' 168',
     1' 169',' 170',' 171',' 172',' 173',' 174',' 175',' 176',' 177',
     1' 178',' 179',' 180',' 181',' 182',' 183',' 184',' 185',' 186',
     1' 187',' 188',' 189',' 190',' 191',' 192',' 193',' 194',' 195',
     1' 196',' 197',' 198',' 199',' 200'/
      data npt201 /
     2' 201',' 202',' 203',' 204',' 205',' 206',' 207',' 208',' 209',
     2' 210',' 211',' 212',' 213',' 214',' 215',' 216',' 217',' 218',
     2' 219',' 220',' 221',' 222',' 223',' 224',' 225',' 226',' 227',
     2' 228',' 229',' 230',' 231',' 232',' 233',' 234',' 235',' 236',
     2' 237',' 238',' 239',' 240',' 241',' 242',' 243',' 244',' 245',
     2' 246',' 247',' 248',' 249',' 250',
     2' 251',' 252',' 253',' 254',' 255',' 256',' 257',' 258',' 259',
     2' 260',' 261',' 262',' 263',' 264',' 265',' 266',' 267',' 268',
     2' 269',' 270',' 271',' 272',' 273',' 274',' 275',' 276',' 277',
     2' 278',' 279',' 280',' 281',' 282',' 283',' 284',' 285',' 286',
     2' 287',' 288',' 289',' 290',' 291',' 292',' 293',' 294',' 295',
     2' 296',' 297',' 298',' 299',' 300'/
      data npt301 /
     3' 301',' 302',' 303',' 304',' 305',' 306',' 307',' 308',' 309',
     3' 310',' 311',' 312',' 313',' 314',' 315',' 316',' 317',' 318',
     3' 319',' 320',' 321',' 322',' 323',' 324',' 325',' 326',' 327',
     3' 328',' 329',' 330',' 331',' 332',' 333',' 334',' 335',' 336',
     3' 337',' 338',' 339',' 340',' 341',' 342',' 343',' 344',' 345',
     3' 346',' 347',' 348',' 349',' 350',
     3' 351',' 352',' 353',' 354',' 355',' 356',' 357',' 358',' 359',
     3' 360',' 361',' 362',' 363',' 364',' 365',' 366',' 367',' 368',
     3' 369',' 370',' 371',' 372',' 373',' 374',' 375',' 376',' 377',
     3' 378',' 379',' 380',' 381',' 382',' 383',' 384',' 385',' 386',
     3' 387',' 388',' 389',' 390',' 391',' 392',' 393',' 394',' 395',
     3' 396',' 397',' 398',' 399',' 400'/
      data npt401 /
     4' 401',' 402',' 403',' 404',' 405',' 406',' 407',' 408',' 409',
     4' 410',' 411',' 412',' 413',' 414',' 415',' 416',' 417',' 418',
     4' 419',' 420',' 421',' 422',' 423',' 424',' 425',' 426',' 427',
     4' 428',' 429',' 430',' 431',' 432',' 433',' 434',' 435',' 436',
     4' 437',' 438',' 439',' 440',' 441',' 442',' 443',' 444',' 445',
     4' 446',' 447',' 448',' 449',' 450',
     4' 451',' 452',' 453',' 454',' 455',' 456',' 457',' 458',' 459',
     4' 460',' 461',' 462',' 463',' 464',' 465',' 466',' 467',' 468',
     4' 469',' 470',' 471',' 472',' 473',' 474',' 475',' 476',' 477',
     4' 478',' 479',' 480',' 481',' 482',' 483',' 484',' 485',' 486',
     4' 487',' 488',' 489',' 490',' 491',' 492',' 493',' 494',' 495',
     4' 496',' 497',' 498',' 499',' 500'/
      data npt501 /
     5' 501',' 502',' 503',' 504',' 505',' 506',' 507',' 508',' 509',
     5' 510',' 511',' 512',' 513',' 514',' 515',' 516',' 517',' 518',
     5' 519',' 520',' 521',' 522',' 523',' 524',' 525',' 526',' 527',
     5' 528',' 529',' 530',' 531',' 532',' 533',' 534',' 535',' 536',
     5' 537',' 538',' 539',' 540',' 541',' 542',' 543',' 544',' 545',
     5' 546',' 547',' 548',' 549',' 550',
     5' 551',' 552',' 553',' 554',' 555',' 556',' 557',' 558',' 559',
     5' 560',' 561',' 562',' 563',' 564',' 565',' 566',' 567',' 568',
     5' 569',' 570',' 571',' 572',' 573',' 574',' 575',' 576',' 577',
     5' 578',' 579',' 580',' 581',' 582',' 583',' 584',' 585',' 586',
     5' 587',' 588',' 589',' 590',' 591',' 592',' 593',' 594',' 595',
     5' 596',' 597',' 598',' 599',' 600'/
      data npt601 /
     6' 601',' 602',' 603',' 604',' 605',' 606',' 607',' 608',' 609',
     6' 610',' 611',' 612',' 613',' 614',' 615',' 616',' 617',' 618',
     6' 619',' 620',' 621',' 622',' 623',' 624',' 625',' 626',' 627',
     6' 628',' 629',' 630',' 631',' 632',' 633',' 634',' 635',' 636',
     6' 637',' 638',' 639',' 640',' 641',' 642',' 643',' 644',' 645',
     6' 646',' 647',' 648',' 649',' 650',
     6' 651',' 652',' 653',' 654',' 655',' 656',' 657',' 658',' 659',
     6' 660',' 661',' 662',' 663',' 664',' 665',' 666',' 667',' 668',
     6' 669',' 670',' 671',' 672',' 673',' 674',' 675',' 676',' 677',
     6' 678',' 679',' 680',' 681',' 682',' 683',' 684',' 685',' 686',
     6' 687',' 688',' 689',' 690',' 691',' 692',' 693',' 694',' 695',
     6' 696',' 697',' 698',' 699',' 700'/
      data npt701 /
     7' 701',' 702',' 703',' 704',' 705',' 706',' 707',' 708',' 709',
     7' 710',' 711',' 712',' 713',' 714',' 715',' 716',' 717',' 718',
     7' 719',' 720',' 721',' 722',' 723',' 724',' 725',' 726',' 727',
     7' 728',' 729',' 730',' 731',' 732',' 733',' 734',' 735',' 736',
     7' 737',' 738',' 739',' 740',' 741',' 742',' 743',' 744',' 745',
     7' 746',' 747',' 748',' 749',' 750',
     7' 751',' 752',' 753',' 754',' 755',' 756',' 757',' 758',' 759',
     7' 760',' 761',' 762',' 763',' 764',' 765',' 766',' 767',' 768',
     7' 769',' 770',' 771',' 772',' 773',' 774',' 775',' 776',' 777',
     7' 778',' 779',' 780',' 781',' 782',' 783',' 784',' 785',' 786',
     7' 787',' 788',' 789',' 790',' 791',' 792',' 793',' 794',' 795',
     7' 796',' 797',' 798',' 799',' 800'/
      data npt801 /
     8' 801',' 802',' 803',' 804',' 805',' 806',' 807',' 808',' 809',
     8' 810',' 811',' 812',' 813',' 814',' 815',' 816',' 817',' 818',
     8' 819',' 820',' 821',' 822',' 823',' 824',' 825',' 826',' 827',
     8' 828',' 829',' 830',' 831',' 832',' 833',' 834',' 835',' 836',
     8' 837',' 838',' 839',' 840',' 841',' 842',' 843',' 844',' 845',
     8' 846',' 847',' 848',' 849',' 850',
     8' 851',' 852',' 853',' 854',' 855',' 856',' 857',' 858',' 859',
     8' 860',' 861',' 862',' 863',' 864',' 865',' 866',' 867',' 868',
     8' 869',' 870',' 871',' 872',' 873',' 874',' 875',' 876',' 877',
     8' 878',' 879',' 880',' 881',' 882',' 883',' 884',' 885',' 886',
     8' 887',' 888',' 889',' 890',' 891',' 892',' 893',' 894',' 895',
     8' 896',' 897',' 898',' 899',' 900'/
      data npt901 /
     9' 901',' 902',' 903',' 904',' 905',' 906',' 907',' 908',' 909',
     9' 910',' 911',' 912',' 913',' 914',' 915',' 916',' 917',' 918',
     9' 919',' 920',' 921',' 922',' 923',' 924',' 925',' 926',' 927',
     9' 928',' 929',' 930',' 931',' 932',' 933',' 934',' 935',' 936',
     9' 937',' 938',' 939',' 940',' 941',' 942',' 943',' 944',' 945',
     9' 946',' 947',' 948',' 949',' 950',
     9' 951',' 952',' 953',' 954',' 955',' 956',' 957',' 958',' 959',
     9' 960',' 961',' 962',' 963',' 964',' 965',' 966',' 967',' 968',
     9' 969',' 970',' 971',' 972',' 973',' 974',' 975',' 976',' 977',
     9' 978',' 979',' 980',' 981',' 982',' 983',' 984',' 985',' 986',
     9' 987',' 988',' 989',' 990',' 991',' 992',' 993',' 994',' 995',
     9' 996',' 997',' 998',' 999','1000'/
      data np1001 /'   o','   i','   j','   k','    '/
      data zero,one,two /0.0d+00,1.0d+00,2.0d+00/
      data iplus,iminus /'+','-'/
c
      if(iunit.eq.0) go to 1000
c
c     ----- read in cartesian coordinates -----
c
      read(ir,9999) atom,znuc,x,y,z
      if(iunit.lt.0) return
      x=x/unit
      y=y/unit
      z=z/unit
      return
c
c     ----- read in internal coordinates.
c           convert to cartesian coordinates. -----
c
c     bond lengths must be given in angstroms,
c     angles must be given in degrees.
c
c
c     default values_
c
c     bond length = 0.0 angstrom
c     alpha = 0.0 degree
c     beta  = 0.0 degree
c     sign  = +
c     connection type = lc = linear connection
c     connection points = 1. origin of master frame
c                         2. unit point on x-axis of master frame
c                         3. unit point on y-axis of master frame
c
c     possible connections are_
c          1. linear connnection                     -  lc-
c          2. planar central connection              - pcc-
c          3. non-planar central connection          -npcc-
c          4. central connection with polar angle    -ccpa-
c          5. planar terminal connection             - ptc-
c          6. terminal connection with torsion       - tct-
c
 1000 continue
      read(ir,9998) atom,znuc,jconx,r,alph,bet,jsign,
     1 (ipt(i),i=1,3)
      if(atom.eq.endwrd) return
      natin=natin+1
      bond(natin)=r/unit
      alpha(natin)=alph*pi2/degree
      beta(natin)=bet*pi2/degree
      jmax=natmax+5
      do 1200 i=1,3
      iat=natmax+6
      do 1100 j=1,jmax
      if(ipt(i).eq.npt(j)) iat=j
 1100 continue
      if(iat.eq.natmax+5.and.i.eq.1) iat=natmax+1
      if(iat.eq.natmax+5.and.i.eq.2) iat=natmax+2
      if(iat.eq.natmax+5.and.i.eq.3) iat=natmax+3
      if(iat.ge.natmax+1.and.iat.le.natmax+4) go to 1150
      if(iat.lt.natin) go to 1150
      write(iwr,9997) natin
      call hnderr(3,errmsg)
 1150 continue
      jat=iat
      if(iat.le.natmax) jat=inatom(iat)
 1200 iatcon(i,natin)=jat
      jmax=natmax+4
      do 1300 j=1,jmax
      do 1300 i=1,3
 1300 cc(i,j)=zero
      cc(1,natmax+2)=one
      cc(2,natmax+3)=one
      cc(3,natmax+4)=one
      if(nat.eq.0) go to 1500
      do 1400 j=1,nat
      do 1400 i=1,3
 1400 cc(i,j)=c(i,j)
 1500 continue
      if(jsign.ne.iminus) jsign=iplus
      sign(natin)=one
      if(jsign.eq.iminus) sign(natin)=-one
      kconx=8
      do 1600 k=1,7
 1600 if(jconx.eq.nconx(k)) kconx=k
      if(kconx.le.7) go to 1700
      write(iwr,9996) natin
      call hnderr(3,errmsg)
 1700 if(kconx.eq.7) kconx=1
      iconx(natin)=kconx
c
c     ----- calculate cartesian coordinates -----
c
      iat1=iatcon(1,natin)
      iat2=iatcon(2,natin)
cnot      rab=rij(cc,iat1,iat2)
      do 1800 i=1,3
 1800 t(i,1)=(cc(i,iat2)-cc(i,iat1))/rab
      if(kconx.gt.1) go to 2000
c
c     ----- linear connection -----
c
      dum=-sign(natin)*bond(natin)
      do 1900 i=1,3
 1900 xyz(i)=dum*t(i,1)
      go to 3000
c
c     ----- tri-atomic connection -----
c
 2000 continue
      iat3=iatcon(3,natin)
cnot      rac=rij(cc,iat1,iat3)
cnot      rbc=rij(cc,iat2,iat3)
c
c     ----- define local frame -----
c
      do 2100 i=1,3
 2100 t(i,2)=cc(i,iat3)-cc(i,iat1)
      dot=t(1,1)*t(1,2)+t(2,1)*t(2,2)+t(3,1)*t(3,2)
      do 2200 i=1,3
 2200 t(i,2)=t(i,2)-dot*t(i,1)
      dot=t(1,2)*t(1,2)+t(2,2)*t(2,2)+t(3,2)*t(3,2)
      dot= sqrt(dot)
      do 2300 i=1,3
 2300 t(i,2)=t(i,2)/dot
      t(1,3)=t(2,1)*t(3,2)-t(3,1)*t(2,2)
      t(2,3)=t(3,1)*t(1,2)-t(1,1)*t(3,2)
      t(3,3)=t(1,1)*t(2,2)-t(2,1)*t(1,2)
c
c     ----- define polar coodinates -----
c
      alph=alpha(natin)
      phi=alph
      go to (2400,2500,2600,2700,2500,2700),kconx
 2400 continue
      call hnderr(3,errmsg)
 2500 theta=zero
      go to 2800
 2600 bet=beta(natin)
      gam= acos((rab*rab+rac*rac-rbc*rbc)/(two*rab*rac))
      theta=
     1  acos(( cos(bet)- cos(alph)* cos(gam))/( sin(alph)* sin(gam)))
      go to 2850
 2700 bet=beta(natin)
      theta=bet
      go to 2850
 2800 continue
      phi=phi*sign(natin)
 2850 continue
      theta=theta*sign(natin)
      r=bond(natin)
c
c     ----- get cartesian coordinates in local frame -----
c
      xyz0(1)=r* cos(phi)
      xyz0(2)=r* sin(phi)* cos(theta)
      xyz0(3)=r* sin(phi)* sin(theta)
c
c     ----- get cartesian coordinates in master frame -----
c
      do 2900 i=1,3
      xyz(i)=zero
      do 2900 j=1,3
      xyz(i)=xyz(i)+t(i,j)*xyz0(j)
 2900 continue
 3000 continue
      x=xyz(1)+cc(1,iat1)
      y=xyz(2)+cc(2,iat1)
      z=xyz(3)+cc(3,iat1)
      return
 9999 format(a8,2x,f5.3,3f20.12)
 9998 format(a8,2x,f5.3,1x,a4,3f10.5,4x,a1,3(1x,a4))
 9997 format(' error in the internal coordinates input of atom ',i5)
 9996 format(' error in connection type for atom ',i5)
      end
      subroutine efcprt
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      common/iofile/ir,iw,ip
INCLUDE(../m4/common/infoa)
cmw      common/mollab/anam(128),bnam(128),bflab(1024)
INCLUDE(comdrf/mollab)
cmw      common/efcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
caleko
      common/hefcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
      character *4 iefc
caleko
INCLUDE(comdrf/opt)
      character*4 keyefc
      data keyefc /' efc'/
c
      if(iefc.ne.keyefc) return
c
      return
      end
      subroutine efcint(hefc,h0)
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      logical iandj
      logical norm,double
      logical some,out
      common/times/timlim,ti,tx,tim,ti0
INCLUDE(comdrf/opt)
      common/iofile/ir,iwr,ip
cmw      common/dafile/idaf,nav,ioda(255)
INCLUDE(comdrf/dafil)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
cmw      common/baspar/normf,normp,itol
INCLUDE(comdrf/bas)
cmw      common/efcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc

caleko
      common/hefcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
      character *4 iefc
caleko
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/ssgg/s(225),g(225),ek(225)
cmw      common/ijpair/ia(1)
INCLUDE(comdrf/ijpair)
      common/rys/xx,u(5),w(5),nroots
      dimension hefc(*),h0(*)
      dimension dij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension xv(5,5,5),yv(5,5,5),zv(5,5,5)
      data rln10 /2.30258d+00/
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data pi212 /1.1283791670955d+00/
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
      out =nprint.eq.3
      some=nprint.ne.-5
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- calculate -hefc- matrix -----
c
      if(some) write(iwr,9999)
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
      nroots=(lit+ljt-2)/2+1
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
      g(ij)=zero
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
c     ----- -efc- interaction -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ic=1,nefc
      znuc=-efcz(ic)
      cx=efcc(1,ic)
      cy=efcc(2,ic)
      cz=efcc(3,ic)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      tt=one/(aa+uu)
      t = sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i
      call sxyz
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint*ww
  410 continue
  420 continue
c
      ij=0
      do 450 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 440 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      dum=zero
      do 430 iroot=1,nroots
  430 dum=dum+xv(ix,jx,iroot)*yv(iy,jy,iroot)*zv(iz,jz,iroot)
      dum=dum*(aa1*pi212)
      ij=ij+1
      g(ij)=g(ij)+dum*dij(ij)
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
c     ----- set up -hefc- matrix -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=ia(loci+i)+(locj+j)
      hefc(nn)=g(ij)
 7500 continue
c
 8000 continue
 9000 continue
c
c     ----- add -hefc- to -h0- -----
c
      if(out) call prtrl(hefc,num)
      num2=(num*(num+1))/2
      call daread(idafh,ioda,h0,num2,11)
      do 9100 i=1,num2
 9100 h0(i)=h0(i)+hefc(i)
      call dawrit(idafh,ioda,h0,num2,11,nav)
c
      if(some) write(iwr,9400)
 9400 format(' ...... end of -efc- integrals ......')
      return
 9999 format(/,10x,15(1h-),/,10x,'-efc- integrals',
     1       /,10x,15(1h-))
      end
      subroutine efcder
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
cmw      common/memory/maxcor,maxlcm
INCLUDE(comdrf/mem)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/scfopt)
      character *8 scf, uhf
      common/scm/x(1)
      character *8 errmsg(3)
      data errmsg  /'program ','stop in ','-efcder-'/
      data scf,uhf /'scf     ','uhf     '/
c
c     ----- get core memory -----
c
      call cmem(loadcm)
      locx1=loccm(x(1))
      ngotmx=maxcor-locx1
c
      i10=1
      i20=i10+(num*(num+1))/2
      i30=i20+(num*(num+1))/2
c
      last=i20+1
      if(wfntyp.eq.scf.and.scftyp.eq.uhf) last=i30
      last=last-1
      if(last.le.ngotmx) go to 10
      need=locx1+last
      call needcm(last,need)
      call hnderr(3,errmsg)
   10 continue
c
c     ----- get density matrix -----
c
cnot      call denddm(x(i10),x(i20),num)
c
c     ----- get gradient of -efc- potential -----
c
      call efdint(x(i10))
c
c     ----- reset core memory -----
c
      call cmem(ngotcm)
      if(ngotcm.ne.loadcm) call setc(loadcm)
      return
      end
      subroutine efdint(dab)
      implicit REAL (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      logical norm
      logical out
cmw      common/output/nprint
INCLUDE(comdrf/opt)
      common/iofile/ir,iwr,ip
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
cmw      common/baspar/normf,normp,itol
INCLUDE(comdrf/bas)
cmw      common/efcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc

caleko
      common/hefcpar/efcc(3,1000),efcz(1000),efclab(1000),nefc,iefc
      character *4 iefc
caleko  
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
     1                              ,cx,cy,cz
cmw      common/ijpair/ia(1)
INCLUDE(comdrf/ijpair)
      common/grad12/de(3,128)
      common/rys/xx,u(5),w(5),nroots
      dimension dab(*)
      dimension dij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension  xv  (6,5,5), yv  (6,5,5), zv  (6,5,5)
      dimension dxvdi(5,5,5),dyvdi(5,5,5),dzvdi(5,5,5)
      dimension dnam(3)
      character *8 errmsg(3)
      data errmsg /'program ','stop in ','-efdint-'/
      data dnam   /3he'x,3he'y,3he'z/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
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
      out =nprint.eq.-3
      norm=normf.ne.1.or.normp.ne.1
c
c     ----- calculate -efc- derivative term -----
c
      if(out) write(iwr,9999)
      nder=1
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
      litder=lit+nder
      iat   =i
c
c     ----- jshell -----
c
      do 8000 jj=1,nshell
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
      ljtmod=ljt+2
c
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(iwr,9997) maxrys,lit,ljt,nroots
         call hnderr(3,errmsg)
      endif
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
      do 6000 jg=j1,j2
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
      do 360 j=minj,maxj
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      go to 350
  250 dum2=dum1*cpj
      go to 350
  260 dum2=dum1*cdj
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
cmw  310 dum2=dum1*cgj
  310 go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      nn=ia(max(loci+i,locj+j))
     1  +   min(loci+i,locj+j)
      den=dab(nn)
      den=den+den
      ij=ij+1
  360 dij(ij)=dum2*den
c
c     ----- -efc- derivatives -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ic=1,nefc
      znuc=-efcz(ic)
      cx=efcc(1,ic)
      cy=efcc(2,ic)
      cz=efcc(3,ic)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      tt=one/(aa+uu)
      t = sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
c
      do 410 j=1,ljt
      nj=j
      do 410 i=1,litder
      ni=i
cnot      call dsxyz
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint*ww
  410 continue
c
cnot      call deri(dxvdi(1,1,iroot),dyvdi(1,1,iroot),dzvdi(1,1,iroot),
cnot     1           xv  (1,1,iroot), yv  (1,1,iroot), zv  (1,1,iroot),
cnot     2          lit,ljt,ai)
c
  420 continue
c
      ij=0
      do 450 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      do 440 j=minj,maxj
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      dumx=zero
      dumy=zero
      dumz=zero
      do 430 iroot=1,nroots
      dumx=dumx+dxvdi(ix,jx,iroot)* yv  (iy,jy,iroot)* zv  (iz,jz,iroot)
      dumy=dumy+ xv  (ix,jx,iroot)*dyvdi(iy,jy,iroot)* zv  (iz,jz,iroot)
  430 dumz=dumz+ xv  (ix,jx,iroot)* yv  (iy,jy,iroot)*dzvdi(iz,jz,iroot)
      ij=ij+1
      de(1,iat)=de(1,iat)+dumx*(aa1*pi212*dij(ij))
      de(2,iat)=de(2,iat)+dumy*(aa1*pi212*dij(ij))
      de(3,iat)=de(3,iat)+dumz*(aa1*pi212*dij(ij))
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
      if(.not.out) go to 8000
      write(iwr,9993) ii,jj
      mmax=0
 7500 mmin=mmax+1
      mmax=mmax+8
      if(mmax.gt.nat) mmax=nat
      write(iwr,9996)
      write(iwr,9995) (i,i=mmin,mmax)
      write(iwr,9996)
      do 7600 i=1,3
 7600 write(iwr,9994) dnam(i),(de(i,j),j=mmin,mmax)
      if(mmax.lt.nat) go to 7500
c
 8000 continue
 9000 continue
c
      if(.not.out) go to 9300
      mmax=0
 9100 mmin=mmax+1
      mmax=mmax+8
      if(mmax.gt.nat) mmax=nat
      write(iwr,9996)
      write(iwr,9995) (i,i=mmin,mmax)
      write(iwr,9996)
      do 9200 i=1,3
 9200 write(iwr,9992) dnam(i),(de(i,j),j=mmin,mmax)
      if(mmax.lt.nat) go to 9100
      write(iwr,9998)
 9300 continue
c
      return
 9999 format(/,10x,21(1h-),/,10x,'-efc- derivative term',
     1       /,10x,21(1h-))
 9998 format(' ...... end of -efc- derivative term ......')
 9997 format(' in -efdint- , the rys quadrature in not implemented',
     1       ' beyond -nroots- = ',i3,/,
     2       ' lit,ljt,nroots= ',3i3)
 9996 format(/)
 9995 format(5x,'atom',8(6x,i3,6x))
 9994 format(6x,a3,8e15.7)
 9993 format(' shells ii, jj ',2i5)
 9992 format(6x,a3,8f15.7)
      end
_ENDIF
