c
c     additional hondo routines needed for the nmr
c
c-----------------------------------------------------------------------
c
      subroutine advfil(ijk,ipk,pandk)
      implicit REAL (a-h,o-z)
      logical pandk
      logical pack2e
      common/pckopt/nhex,ntupl,pack2e
      common/intopt/nopk,nok,nosqur
      common/intfil/nintmx
      if(pack2e) go to 10
      read(ijk)
      if(pandk) read(ipk)
      return
   10 call caserr("advfil: should not get here")
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine diaaxs(a,vec,eig,ia,nvec,n,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- routine to substitute diagiv for diagonalization -----
c           of symmetric 3x3 matrix a  in triangular form
c
      character*8 errmsg
      common/iofile/ir,iw
      dimension a(6),h(3,3),vec(3,3),eig(3),ia(1),big(6)
      dimension errmsg(3)
      dimension iky(3),ilifq(3)
      data errmsg /'program ','stop in ','-diaaxs-'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data conv   /1.0d-10/
      data maxit  /50/
      data iky    /0,1,3/
      data ilifq  /0,3,6/
c
      if(n.eq.3.and.ndim.eq.3) go to 10
         write(iw,9999)
         call hnderr(3,errmsg)
   10 continue
c
      do 30 i=1,3
         do 20 j=1,3
            vec(j,i)=zero
   20    continue
         vec(i,i)=one
   30 continue
      call jacobi(a,iky,n,vec,ilifq,n,eig,1,0,1.0d-10)
c
c     ----- check for right handedness, correct if not -----
c
      test =   vec(1,3)*( vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2) )
     1       + vec(2,3)*( vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2) )
     2       + vec(3,3)*( vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2) )
      if(test.gt.zero) return
      if( abs(eig(1)-eig(2)).gt.conv) go to 60
         t = eig(1)
         eig(1) = eig(2)
         eig(2) = t
         do 50 i=1,3
            t = vec(i,1)
            vec(i,1) = vec(i,2)
            vec(i,2) = t
   50    continue
         return
   60 if( abs(eig(2)-eig(3)).gt.conv) go to 80
         t = eig(2)
         eig(2) = eig(3)
         eig(3) = t
         do 70 i=1,3
            t = vec(i,2)
            vec(i,2) = vec(i,3)
            vec(i,3) = t
   70    continue
         return
   80 do 90 i=1,3
         vec(i,3) = - vec(i,3)
   90 continue
      return
 9999 format(/,' -diaaxs- diagonalization only set up for 3x3 matrix')
      end
c
c-----------------------------------------------------------------------
c
      subroutine intout(i1,i2,i3,i4,q4,nn,val,ncall)
      implicit REAL (a-h,o-z)
      logical out
      common/intprt/q(3),v(3),jc,n1(3),j1(3),j2(3),j3(3),j4(3),out
      common/iofile/ir,iw
      if(ncall.ne.0) go to 10
      jc=jc+1
      j1(jc)=i1
      j2(jc)=i2
      j3(jc)=i3
      j4(jc)=i4
      q(jc)=q4
      n1(jc)=nn
      v(jc)=val
      if(jc.lt.3)  return
   10 continue
      if(jc.eq.0) return
      write(iw,9999) (j1(m),j2(m),j3(m),j4(m),q(m),n1(m),v(m),m=1,jc)
      jc=0
 9999 format(3(4i3,f6.3,i5,e20.12))
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine jkbcdf(b00,b01,b10,c00,d00,f00,dij,dkl,
     1                  abv,cv,rwv,numg,nroots)
      implicit REAL (a-h,o-z)
      logical nmaxs,nmaxp,mmaxs,mmaxp
      common/shlgnm/nmaxs,nmaxp,mmaxs,mmaxp
      dimension b00(numg,nroots,3),b01(numg,nroots,3),b10(numg,nroots,3)
      dimension c00(numg,nroots,3)
      dimension d00(numg,nroots,3)
      dimension f00(numg,nroots,3),dij(numg,nroots,3),dkl(numg,nroots,3)
      dimension abv(5,1),cv(18,1)
      dimension rwv(2,numg,nroots)
      data pt5,one /0.5d+00,1.0d+00/
c
      do 40 nr=1,nroots
      do 30 ng=1,numg
      aa =abv(1,ng)
      bb =abv(2,ng)
      rho=abv(3,ng)
      qab=abv(4,ng)
      uu =rho*rwv(1,ng,nr)
      ww =    rwv(2,ng,nr)
      aauu=aa+uu
      bbuu=bb+uu
      f00(ng,nr,1)=ww*qab
      f00(ng,nr,2)=one
      f00(ng,nr,3)=one
      dum2=pt5/(aa*bb+uu*(aa+bb))
      audum=aauu*dum2
      budum=bbuu*dum2
       udum=  uu*dum2
      b00(ng,nr,1)= udum
      b00(ng,nr,2)= udum
      b00(ng,nr,3)= udum
      b01(ng,nr,1)=audum
      b01(ng,nr,2)=audum
      b01(ng,nr,3)=audum
      b10(ng,nr,1)=budum
      b10(ng,nr,2)=budum
      b10(ng,nr,3)=budum
       udum= udum+ udum
      if(mmaxs) go to 10
      audum=audum+audum
      d00(ng,nr,1)= udum*cv( 1,ng) + audum*cv( 2,ng)
      d00(ng,nr,2)= udum*cv( 3,ng) + audum*cv( 4,ng)
      d00(ng,nr,3)= udum*cv( 5,ng) + audum*cv( 6,ng)
   10 if(nmaxs) go to 20
      budum=budum+budum
      c00(ng,nr,1)= udum*cv( 8,ng) + budum*cv( 7,ng)
      c00(ng,nr,2)= udum*cv(10,ng) + budum*cv( 9,ng)
      c00(ng,nr,3)= udum*cv(12,ng) + budum*cv(11,ng)
   20 continue
c
   30 continue
   40 continue
c
      do 60 nr=1,nroots
      do 50 ng=1,numg
      dij(ng,nr,1)=cv(13,ng)
      dij(ng,nr,2)=cv(14,ng)
      dij(ng,nr,3)=cv(15,ng)
      dkl(ng,nr,1)=cv(16,ng)
      dkl(ng,nr,2)=cv(17,ng)
      dkl(ng,nr,3)=cv(18,ng)
   50 continue
   60 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine jkwrys(rwv,abv,numg)
      implicit REAL (a-h,o-z)
      common/hfk/xx,u(9),w(9),nroots
      common/root/yy,t(12),v(12),mroots
      dimension rwv(2,numg,1),abv(5,1)
c
      if(nroots.gt.5) go to 100
c
c     ----- nroots .le. 5 -----
c
      do 20 ng=1,numg
      yy=abv(5,ng)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call root4
      if(nroots.eq.5) call root5
      do 10 nr=1,nroots
      rwv(1,ng,nr)=t(nr)
      rwv(2,ng,nr)=v(nr)
   10 continue
   20 continue
      return
c
  100 continue
c
c     ----- nroots .gt. 5 -----
c
      do 120 ng=1,numg
      xx=abv(5,ng)
      call droot
      do 110 nr=1,nroots
      rwv(1,ng,nr)=u(nr)
      rwv(2,ng,nr)=w(nr)
  110 continue
  120 continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine oehnd(nijg)
      implicit REAL (a-h,o-z)
      parameter (mxgsh =30,mxgsh2=mxgsh*mxgsh)
      parameter (mxprim=MXPRIM)
      parameter (mxshel=MXSHEL)
      parameter (maxorb=MAXORB)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shlord/modshl(mxshel),invshl(mxshel)
      common/ijchrg/acharg(11,mxgsh2),dasi(mxgsh2),dasj(mxgsh2),nij
      common/intmem/iclb,inijg,ibuf,igint,iijklg,igijkl,ignkl,
     1 ignm,idij,idkl,ib00,ib01,ib10,ic00,id00,if00,
     2 isj,isk,isl,idijsi,idijsj,idklsk,idklsl,iabv,icv,irw,ichrg
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/vcore)
      dimension nijg(2,1)
c
c     ----- one-electron charge distribution -----
c
      nij0=0
      do 9000 ii=1,nshell
      iimod=modshl(ii)
      do 8000 jj=1,ii
      jjmod=modshl(jj)
      ish=iimod
      jsh=jjmod
      call oeshel(ish,jsh)
      call oechrg
      iijj=iky(max(iimod,jjmod))+min(iimod,jjmod)
      nijg(1,iijj)=nij0
      nijg(2,iijj)=nij
      if(nij.eq.0) go to 7000
      call oewrit(Q(ichrg),nij0,acharg,dasi,dasj,nij)
 7000 continue
      nij0=nij0+nij
 8000 continue
 9000 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine oeshel(ish,jsh)
      implicit REAL (a-h,o-z)
      parameter (mxatom=MAXAT)
      parameter (mxprim=MXPRIM)
      parameter (mxshel=MXSHEL)
      parameter (mxgsh =30,mxgsh2=mxgsh*mxgsh)
      logical    iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical    spi,spj,spk,spl,spij,spkl,spijkl
      logical    expndi,expndk
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(mxatom),
     1                                          c(3,mxatom)
      common/nshel/
     1ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),cf(mxprim),cg(mxprim),
     1kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     2kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/shltyp/spi,spj,spk,spl,spij,spkl,spijkl
      common/shlxpn/expndi,expndk
      common/shlpar/lit,ljt,lkt,llt,loci,locj,lock,locl,
     1              mini,minj,mink,minl,maxi,maxj,maxk,maxl
      common/shlnum/numi,numj,numk,numl,ij,kl,ijkl
      common/shlinf/ga(mxgsh),cca(mxgsh),ccas(mxgsh),
     1              gb(mxgsh),ccb(mxgsh),ccbs(mxgsh),
     1              gc(mxgsh),ccc(mxgsh),cccs(mxgsh),
     1              gd(mxgsh),ccd(mxgsh),ccds(mxgsh),
     1              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     1              nga,ngb,ngc,ngd
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      dimension cspdfg(mxprim,5)
      equivalence (cspdfg(1,1),cs(1))
c
      iieqjj=ish.eq.jsh
c
c     ----- ishell -----
c
      i=katom(ish)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ish)
      i2=i1+kng(ish)-1
      lit=ktype(ish)
      mini=kmin(ish)
      maxi=kmax(ish)
      numi=maxi-mini+1
      loci=kloc(ish)-mini
      spi=lit.eq.2.and.mini.eq.1
      nga=0
      do 10 i=i1,i2
      nga=nga+1
      ga(nga)=ex(i)
      cca(nga)=cspdfg(i,lit)
      if(spi) ccas(nga)=cspdfg(i,1)/cspdfg(i,2)
   10 continue
c
c     ----- jshell -----
c
      j=katom(jsh)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jsh)
      j2=j1+kng(jsh)-1
      ljt=ktype(jsh)
      minj=kmin(jsh)
      maxj=kmax(jsh)
      numj=maxj-minj+1
      locj=kloc(jsh)-minj
      spj=ljt.eq.2.and.minj.eq.1
      ngb=0
      do 20 j=j1,j2
      ngb=ngb+1
      gb(ngb)=ex(j)
      ccb(ngb)=cspdfg(j,ljt)
      if(spj) ccbs(ngb)=cspdfg(j,1)/cspdfg(j,2)
   20 continue
      rri=((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
      spij=spi.or.spj
      expndi=lit.ge.ljt
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine oechrg
      implicit REAL (a-h,o-z)
      parameter (mxgsh =30,mxgsh2=mxgsh*mxgsh)
      logical iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      logical spi,spj,spk,spl,spij,spkl,spijkl
      logical expndi,expndk
      common/shltol/rtol,dtol
      common/shlequ/iieqjj,kkeqll,ijeqkl,ijgtkl,ijltkl
      common/shltyp/spi,spj,spk,spl,spij,spkl,spijkl
      common/shlxpn/expndi,expndk
      common/shlinf/ga(mxgsh),cca(mxgsh),ccas(mxgsh),
     1              gb(mxgsh),ccb(mxgsh),ccbs(mxgsh),
     1              gc(mxgsh),ccc(mxgsh),cccs(mxgsh),
     1              gd(mxgsh),ccd(mxgsh),ccds(mxgsh),
     1              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     1              nga,ngb,ngc,ngd
      common/ijchrg/acharg(11,mxgsh2),dasi(mxgsh2),dasj(mxgsh2),nij
      data one /1.0d+00/
      data pt5 /0.5d+00/
c
c     ----- -ij- charge distribution -----
c
      xc=xi
      yc=yi
      zc=zi
      dxij=xi-xj
      dyij=yi-yj
      dzij=zi-zj
      if(expndi) go to 10
      xc=xj
      yc=yj
      zc=zj
      dxij=xj-xi
      dyij=yj-yi
      dzij=zj-zi
   10 continue
c
c     ----- - i- primitive           -----
c
      nij=0
      do 40 ia=1,nga
      ai=ga(ia)
      arri=ai*rri
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      cci=cca(ia)
c
c     ----- - j- primitive           -----
c
      jbmax=ngb
      if(iieqjj) jbmax=ia
      do 30 jb=1,jbmax
      aj=gb(jb)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.rtol) go to 30
      daexpa=cci*ccb(jb)* exp(-dum)*aa1
      dum= abs(daexpa)
      if(dum.le.dtol) go to 30
c
      nij=nij+1
      acharg( 1,nij)= daexpa
      acharg( 2,nij)= aa
      acharg( 3,nij)=(axi+aj*xj)*aa1
      acharg( 4,nij)=(ayi+aj*yj)*aa1
      acharg( 5,nij)=(azi+aj*zj)*aa1
      acharg( 6,nij)= xc
      acharg( 7,nij)= yc
      acharg( 8,nij)= zc
      acharg( 9,nij)= dxij
      acharg(10,nij)= dyij
      acharg(11,nij)= dzij
      if(spi) dasi(nij)=ccas(ia)
      if(spj) dasj(nij)=ccbs(jb)
c
      if(.not.iieqjj) go to 30
      if(jb.eq.ia)    go to 20
      acharg(1,nij)=acharg(1,nij)+acharg(1,nij)
   20 continue
      if(.not.spi)    go to 30
      dasi(nij)= ccas(ia)*ccbs(jb)
      dasj(nij)=(ccas(ia)+ccbs(jb))*pt5
   30 continue
   40 continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine oewrit(chrg,nij0,acharg,dasi,dasj,nij)
      implicit REAL (a-h,o-z)
      dimension acharg(11,1),dasi(1),dasj(1)
      dimension chrg(13,1)
c
      do 20 ij=1,nij
      do 10 i=1,11
   10 chrg( i,ij+nij0)=acharg(i,ij)
      chrg(12,ij+nij0)=dasi(ij)
      chrg(13,ij+nij0)=dasj(ij)
   20 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine oeread(chrg,nij0,acharg,dasi,dasj,nij)
      implicit REAL (a-h,o-z)
      dimension acharg(11,1),dasi(1),dasj(1)
      dimension chrg(13,1)
c
      do 20 ij=1,nij
      do 10 i=1,11
   10 acharg(i,ij)=chrg( i,ij+nij0)
      dasi(ij)    =chrg(12,ij+nij0)
      dasj(ij)    =chrg(13,ij+nij0)
   20 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pkread(ip,ik,xp,xk,ix,nxx,nmax)
      implicit REAL (a-h,o-z)
      character*8 errmsg
      logical     pack2e
      common/iofile/ir,iw
      common/pcklab/labsiz
      common/pckopt/nhex,ntupl,pack2e
      common/intopt/nopk,nok,nosqur
      common/intfil/nintmx
      dimension xp(nmax),xk(nmax),ix(labsiz*nmax)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-pkread-'/
c
      if(pack2e) go to 10
      read(ip,end=20,err=20) xp,ix,nx
      read(ik,end=20,err=20) xk,ix,nx
      nxx=nx
      return
   10 call caserr("pkread: should not get here")
      return
c
c     ----- error ... -----
c
   20 continue
      write(iw,*) 'error while reading -p and k- . stop.'
      call hnderr(3,errmsg)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine tfsq(h,f,q1,t,m,n,ndim)
      implicit REAL (a-h,o-z)
c
c     ----- h(m,m) = q1(dagger)(n,m) * f(n,n) * q1(n,m) -----
c
      character*8 errmsg
      common/iofile/ir,iw
      common/memvec/iadvec,maxvec
      common/zeropt/zertol
INCLUDE(../m4/common/vcore)
      dimension h(ndim,*),f(ndim,*),q1(ndim,*),t(*)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-tfsq  -'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
c
      i10=1+iadvec
      i20=i10+n*n
      i30=i20+n*m
      i40=i30+m*m
      last=i40-1
      need=last-iadvec
      if(need.gt.maxvec) then
         write(iw,9999) need
         call hnderr(3,errmsg)
      endif
c
      do j=1,n
         do i=1,n
            dum = f(i,j)
            Q(i+(j-1)*n+i10-1)=dum
         enddo
      enddo
c
      call dgemm('n','n',n,m,n,one,Q(i10),n,q1,ndim,zero,
     +           Q(i20),n)
c
      call dgemm('t','n',m,m,n,one,q1,ndim,Q(i20),n,zero,
     +           Q(i30),m)
c
      do j=1,m
         do i=1,m
            dum=Q(i+(j-1)*m+i30-1)
            if( abs(dum).lt.zertol) dum=zero
            h(i,j)=dum
         enddo
      enddo
      return
 9999 format(' in -tfsq- , not enough scratch space for vector.',
     1       ' needed = ',i10)
      end

