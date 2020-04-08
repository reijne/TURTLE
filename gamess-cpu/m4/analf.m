c *****************************************************
c *****************************************************
c             = analf  =
c             = includes solvation code =
c ******************************************************
c             =   potentials generator  =
c             =   dvsg jan. 1991        =
c ******************************************************
c ******************************************************
       subroutine clrefg(efgmol)
       implicit REAL  (a-h,p-w), integer(i-n), logical(o)
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/psscrf)
       dimension efgmol(2,*)
c
       ist=nat+1
       ifn=nat+itotbq
       do 1 i=1,2
          do 2 j=ist,ifn
                efgmol(i,j)=0.0d0
2          continue
1       continue
c
       return
      end
      function denuc(n,cz,c,numbq,iwr)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension cz(*),c(3,*)
      data dzero /0.0d0/
      denuc = dzero
      if (n .eq. 1) return
      do 120 i = 2,n
      ni = i-1
      do 120 j = 1,ni
      rr = dzero
      do 100 k = 1,3
  100 rr = rr+(c(k,i)-c(k,j))**2
  120 denuc = denuc+cz(i)*cz(j)/dsqrt(rr)
c     write(iwr,160)denuc
c 160 format(1x,'    nuclear energy = ',f20.10)
      enuc1=denuc
      do 130 i=1,n
       do 140 j=n+1,n+numbq
      rr = dzero
      do 150 k = 1,3
  150 rr = rr+(c(k,i)-c(k,j))**2
      denuc = denuc+cz(i)*cz(j)/dsqrt(rr)
  140 continue
  130 continue
      enuc2=denuc-enuc1
c     write(iwr,170)enuc2
c 170 format(1x,'chg-nuclear energy = ',f20.10)
      return
      end
      subroutine dpogrd(ngrid,dens,ga,coord,cnrmpt,
_IF(parallel)
c***   **MPP**
     + efgsum,efgmol,iwr)
c***   **MPP**
_ENDIF
_IFN(parallel)
     + efgmol,iwr)
_ENDIF
      implicit REAL  (a-h,o-z)
      logical iandj,norm,double
INCLUDE(common/sizes)
c
c     common for psscrf
INCLUDE(common/psscrf)
c
INCLUDE(common/restar)
INCLUDE(common/root)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
      common/junk/s(225),g(225),
     + xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,
     + tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
      common/blkin/dxyz(4),ggg(225),ft(225),dij(225),
     + xin(125),yin(125),zin(125),
     + ijx(225),ijy(225),ijz(225)
      dimension dens(*),ga(*)
      dimension coord(3,*),cnrmpt(3,*),efgmol(2,*)
INCLUDE(common/parallel)
_IF(parallel)
c***   **MPP**
      dimension efgsum(2,*)
_ENDIF
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      data zero,one /0.0d+00,1.0d+00/
c     data pt5,two,three/0.5d+00,2.0d+00,3.0d+00/
      data pi212 /1.1283791670955d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data rln10 /2.30258d+00/
c
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
c
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
      ist=nat+1
      ifn=nat+ngrid
      np36=ngrid * 225 * 2
_IF(parallel)
c***   **MPP**
      call pg_dlbreset
      next = ipg_dlbtask()
      ifn2 = ifn + ifn
      call vclr(efgmol,1,ifn2)
c***   **MPP**
_ENDIF
c
c     ----- ishell
c
_IF(parallel)
      do 9000 ii=nshell,1,-1
_ELSE
      do 9000 ii=1,nshell
_ENDIF
      i=katom(ii)
      xi=coord(1,i)
      yi=coord(2,i)
      zi=coord(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell
c
_IF(parallel)
      do 8000 jj=ii,1,-1
c***   **MPP**
      icount_dlb = icount_dlb + 1
      if (icount_dlb.eq.next) then
c***   **MPP**
_ELSE
      do 8000 jj=1,ii
_ENDIF
      j=katom(jj)
      xj=coord(1,j)
      yj=coord(2,j)
      zj=coord(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
      nroots=(lit+ljt-2)/2+1
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      iandj=ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
      ij=0
      max=maxj
      do 50 i=mini,maxi
      mx=ix(i)
      my=iy(i)
      mz=iz(i)
      if(iandj) max=i
      do 50 j=minj,max
      ij=ij+1
      ijx(ij)=mx+jx(j)
      ijy(ij)=my+jy(j)
      ijz(ij)=mz+jz(j)
   50 continue
c
      call vclr(ga,1,np36)
c
c     ----- i primitive
c
      jgmax=j2
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
      cgi=cg(ig)
c
c     ----- j primitive
c
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac=dexp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor
c
      double=iandj.and.ig.ne.jg
      max=maxj
      nn=0
      do 170 i = mini,maxi
      go to (200,220,280,280,
     +       240,280,280,260,280,280,
     +       210,280,280,215,280,280,280,280,280,217,
     +       222,280,280,224,280,280,280,280,280,
     +       226,280,280,228,280,280), i
  200 dum1 = csi*fac
      go to 280
  220 dum1 = cpi*fac
      go to 280
  240 dum1 = cdi*fac
      go to 280
  260 if (norm) dum1 = dum1*sqrt3
      go to 280
  210 dum1=cfi*fac
      go to 280
  215 if (norm) dum1=dum1*sqrt5
      go to 280
  217 if (norm) dum1=dum1*sqrt3
      go to 280
  222 dum1=cgi*fac
      go to 280
  224 if (norm) dum1 = dum1*sqrt7
      go to 280
  226 if (norm) dum1 = dum1*sqrt5/sqrt3
      go to 280
  228 if (norm) dum1 = dum1*sqrt3
  280 if (iandj) max = i
      do 170 j = minj,max
      go to (300,340,400,400,
     +       360,400,400,380,400,400,
     +       310,400,400,315,400,400,400,400,400,317,
     +       322,400,400,324,400,400,400,400,
     +       400,326,400,400,328,400,400),j
  300 dum2 = dum1*csj
      if ( .not. double) go to 400
      if (i .gt. 1) go to 320
      dum2 = dum2+dum2
      go to 400
  320 dum2 = dum2+csi*cpj*fac
      go to 400
  340 dum2 = dum1*cpj
      if (double) dum2 = dum2+dum2
      go to 400
  360 dum2 = dum1*cdj
      if (double) dum2 = dum2+dum2
      go to 400
  380 if (norm) dum2 = dum2*sqrt3
      go to 400
  310 dum2=dum1*cfj
      if (double) dum2 = dum2+dum2
      go to 400
  315 if (norm) dum2 = dum2*sqrt5
      go to 400
  317 if (norm) dum2 = dum2*sqrt3
      go to 400
  322 dum2=dum1*cgj
      if (double) dum2 = dum2 + dum2
      go to 400
  324 if (norm) dum2 = dum2*sqrt7
      go to 400
  326 if (norm) dum2 = dum2*sqrt5/sqrt3
      go to 400
  328 if (norm) dum2 = dum2*sqrt3
  400 nn = nn+1
  170 dij(nn) = dum2
c
c     ----- nuclear attraction
c
      dum=pi212*aa1
      do 500 i=1,ij
  500 dij(i)=dij(i)*dum
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
c
c ----- altered to do 2 sets of potentials
c       at given points for psscrf
c
      ipt=0
      do 450 idg=1,2
      do 450 igg=ist,ifn
      znuc=one
      ipt = ipt+1
        if(idg.eq.1) then
          cx=coord(1,igg)
          cy=coord(2,igg)
          cz=coord(3,igg)
        else
          cx=cnrmpt(1,igg)
          cy=cnrmpt(2,igg)
          cz=cnrmpt(3,igg)
        endif
      pp=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3)call rt123
      if(nroots.eq.4)call roots4
      if(nroots.eq.5)call roots5
      mm=0
      do 420 k=1,nroots
      uu=aa*u(k)
      ww=w(k)*znuc
      tt=one/(aa+uu)
      t= dsqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      in=-5+mm
      do 410 i=1,lit
      in=in+5
      ni=i
      do 410 j=1,ljt
      jn=in+j
      nj=j
      call stvint
      xin(jn)=xint
      yin(jn)=yint
      zin(jn)=zint*ww
  410 continue
  420 mm=mm+25
      do 440 i=1,ij
      mx=ijx(i)
      my=ijy(i)
      mz=ijz(i)
      dum=zero
      mm=0
      do 430 k=1,nroots
      dum=dum+xin(mx+mm)*yin(my+mm)*zin(mz+mm)
  430 mm=mm+25
  440 g(i)= dum*dij(i)
      index=225*(ipt-1)
      max=maxj
      nn2=0
      do  460 i = mini,maxi
      if (iandj) max = i
      do  460 j = minj, max
      index=index+1
      nn2 = nn2 +1
  460 ga(index) = ga(index) + g(nn2)
  450 continue
 6000 continue
 7000 continue
c
c     ----- set up  efgmol matrix
c
      ipt=0
      do 470 idg=1,2
      do 470 igg=ist,ifn
      ipt = ipt + 1
      index=225*(ipt-1)
      vt=zero
      max = maxj
      do  480 i=mini,maxi
      li=loci+i
      in=iky(li)
      if (iandj) max=i
      do  480 j =minj,max
      lj= locj + j
      jn = lj + in
      index=index+1
      dw = dens(jn) * ga(index)
      if(li.ne.lj) dw= dw + dw
      vt = vt +dw
  480 continue
      efgmol(idg,igg)=efgmol(idg,igg)-vt
  470 continue
_IF(parallel)
      next = ipg_dlbtask()
      endif
_ENDIF
 8000 continue
 9000 continue
_IF(parallel)
c***   ***node-MPP***
      call pg_dgop(903,efgmol,ifn2,'+')
      call vadd(efgmol,1,efgsum,1,efgsum,1,ifn2)
      call pg_dlbpush
c***   ***node-MPP***
_ENDIF
      ss=cpulft(1)
      write(iwr,3001)ss
 3001 format(' potential surface complete at ',f8.2,' secs')
      return
      end
      subroutine dponuc(coord,chg,cnrmpt,efgmol)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
c     common for psscrf
c
INCLUDE(common/psscrf)
INCLUDE(common/infoa)
c
      dimension coord(3,*),chg(*),cnrmpt(3,*),efgmol(2,*)
c
      data dzero/0.0d0/
c
c ... setup nuclear contribution to potential grid
c ... changed to use surface point coordinates
c
      ist=nat+1
      ifn=nat+itotbq
c
      do 3 idg=1,2
         do 1 i=ist,ifn
           if(idg.eq.1) then
             pppp=coord(1,i)
             qqqq=coord(2,i)
             rrrr=coord(3,i)
             top=dzero
           else
             pppp=cnrmpt(1,i)
             qqqq=cnrmpt(2,i)
             rrrr=cnrmpt(3,i)
             top=dzero
           endif
           do 2 ii=1,nat
             ctj=pppp-coord(1,ii)
             pow1=qqqq-coord(2,ii)
             pow2=rrrr-coord(3,ii)
             ctj=ctj*ctj+pow1*pow1+pow2*pow2
             top=top+chg(ii)/dsqrt(ctj)
2          continue
1        efgmol(idg,i)=top
3     continue
c
      return
      end
      subroutine dtgrid(real8,coord,atmchg,cnrmpt,efgmol)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension real8(*),coord(3,*),atmchg(*),cnrmpt(3,*),efgmol(2,*)
c
c     common for psscrf
c
INCLUDE(common/psscrf)
INCLUDE(common/mapper)
INCLUDE(common/timez)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/restar)
      character*10 charwall
c
c     data m2/2/
c
c     oend = .true.
c
      i10=1
      i20=i10+nx
_IF(parallel)
c***   **MPP**
      i30=i20+225*itotbq*2
      non2=(nat+itotbq)*2
      last=i30+non2
c***   **MPP**
_ENDIF
      write(iwr,51004) cpulft(1),charwall()
51004 format(/' **** commence potential grid generation at ',f8.2,
     + ' seconds',a10,' wall')
c
c  ... commence grid generation
c
      call dponuc(coord,atmchg,cnrmpt,efgmol)
c
      call rdedx(real8(i10),nx,ibl3pa,idaf)
      call dpogrd(itotbq,real8(i10),real8(i20),
_IF(parallel)
c***   **MPP**
     + coord,cnrmpt,efgmol,real8(i30),iwr)
c***   **MPP**
_ELSE
     + coord,cnrmpt,efgmol,iwr)
_ENDIF
c
c     if(.not.oend) then
c      cpu=cpulft(1)
c      call clredx
c      irest=10
c      tim=timlim+0.5d0
c      return
c     endif
c
      return
      end
      subroutine mod1e(q,coord,atmchg)
      implicit REAL       (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/scra7)
INCLUDE(common/tran)
INCLUDE(common/statis)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/atmol3)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
c
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/cslosc)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
      common/blkcore/corev(512),charge(10)
c
      dimension q(*),coord(3,*),atmchg(*),o1e(6)
      character*7 fnm
      character*5 snm
      data fnm,snm/"analf.m","mod1e"/
c
      data m19/19/
c
      call cpuwal(begin,ebegin)
c
c ----- first calculate s,t,f, store in section 492 and on ed7
c
c read in symmetry array - iso
c
      i10 = igmem_alloc_inf(nw196(5),fnm,snm,'i10',IGMEM_DEBUG)
      call rdedx(q(i10),nw196(5),ibl196(5),idaf)
      call nucatt(q,coord,atmchg,q(i10),nshell)
      call gmem_free_inf(i10,fnm,snm,'i10')
c
      nav = lenwrd()
      l2 = num*(num+1)/2
c     l3 = num*num
c
c      determine overall memory required
c
      i10 = 0
      i20 = i10 + l2
      i30 = i20 + l2
      max1 = i30 + l2
      i21 = i20 + mxgaus
      i22 = i21 + mxgaus
      i23 = i22 + mxprms*num
      i24 = i23 + mxprms*num
      i25 = i24 + (7*maxat+1)/nav
      i26 = i25 + (3*num+1)/nav
      i27 = i26 + (3*num+1)/nav
      i28 = i27 + (3*num+1)/nav
      i29 = i28 + (8*num+1)/nav
      i31 = i29 + (8*num+1)/nav
      max2 = i31 + (num+1)/nav
      max1 = max(max1,max2)
      length = max1
c
c     now allocate required memory and determine
c     pointers
c
      i10 = igmem_alloc_inf(length,fnm,snm,"i10",IGMEM_DEBUG)
      i20 = i10 + l2
      i30 = i20 + l2
c
c -----  .. compute storage requirements for adaption
c
      max1 = i30 + l2
c
c...   everybody gets a dummy symmetry adaption at least
c...   adapt will remember this
c
      i21 = i20 + mxgaus
      i22 = i21 + mxgaus
      i23 = i22 + mxprms*num
      i24 = i23 + mxprms*num
      i25 = i24 + (7*maxat+1)/nav
      i26 = i25 + (3*num+1)/nav
      i27 = i26 + (3*num+1)/nav
      i28 = i27 + (3*num+1)/nav
      i29 = i28 + (8*num+1)/nav
      i31 = i29 + (8*num+1)/nav
      max2 = i31 + (num+1)/nav
c
c ----- restore the s,t,t+v integrals
c
      do loop = 1,3
        o1e(loop) = .true.
        o1e(loop+3) = .false.
      enddo
      call getmat(q(i10),q(i20),q(i30),q(i10),q(i10),q(i10),
     +            charge,num,o1e,ionsec)
      ibl7s = ibl7la
      call wrt3(q(i10),l2,ibl7s,num8)
      ibl7t = iposun(num8)
      call wrt3(q(i20),l2,ibl7t,num8)
      ibl7f = iposun(num8)
      call wrt3(q(i30),l2,ibl7f,num8)
      ibl7st = iposun(num8)
      call adapt(q(i20),q(i21),q(i22),q(i23),q(i24),q(i25),q(i26),q(i27)
     +           ,q(i28),q(i29),q(i31),num)
c
c ----- now transform 1e-integrals to sabf and output to ed7
c
      if (odscf .and. odnew) then
         call anorm(q(i10),q)
c ----- restore s-matrix
         call rdedx(q(i10),l2,ibl7s,num8)
      end if
c
      call tranp(q(i10),q(i20))
      call wrt3(q(i20),l2,ibl7st,num8)
      ibl7tt = iposun(num8)
      call rdedx(q(i10),l2,ibl7t,num8)
      call tranp(q(i10),q(i20))
      call wrt3(q(i20),l2,ibl7tt,num8)
      ibl7ft = iposun(num8)
      call rdedx(q(i10),l2,ibl7f,num8)
      call tranp(q(i10),q(i20))
c
      call wrt3(q(i20),l2,ibl7ft,num8)
c
c ----- now decide on blocks for density matrix storage
c ----- and eigen value storage on dumpfile .. create sufficient
c ----- space for both -a- and -b- set vectors
c ----- dumped to section 497 of dumpfile
c
      ibl7la = iposun(num8)
      length = lensec(l2)
      lene = lensec(num)
      len2 = 2*(lene+length)
      call secput(isect(497),m19,len2,ibl3pa)
      ibl3ea = ibl3pa + length
      ibl3pb = ibl3ea + 1
      ibl3eb = ibl3pb + length
c
c ----- reset core
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
c ----- now act on vectors specification
c
      call revind
c
c     ----- read saved mo*s -----
c
      if(mina.le.0) mina = mouta
      zguess = 'mosaved'
      if (zscftp.ne.'vb') call mofile(q,zscftp,0)
c
      call timana(2)
      return
      end
      subroutine mod1ed(coord,atmchg,idump)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/timez)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/psscrf)
c
      dimension coord(3,*),atmchg(*)
      data m16/16/
c
      non = nat + itotbq
      nav = lenwrd()
      mword1 = 3 + (6/nav)
      mword2 = 3 * non
      mword3 = non
      lenbb = lensec(mword1) + lensec(mword2) + lensec(mword3)
      call secput(msecp,m16,lenbb,iblk16)
      if(idump.le.0) then
       call rdedx(ptspac,mword1,iblk16,idaf)
       call reads(coord,mword2,idaf)
       call reads(atmchg,mword3,idaf)
      else
       call wrt3(ptspac,mword1,iblk16,idaf)
       call wrt3s(coord,mword2,idaf)
       call wrt3s(atmchg,mword3,idaf)
      endif
      return
      end
      subroutine nucatt(q,coord,atmchg,iso,nshels)
      implicit REAL       (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/segm)
INCLUDE(common/timez)
INCLUDE(common/symtry)
INCLUDE(common/prints)
INCLUDE(common/restar)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
INCLUDE(common/statis)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
      common/junk/s(225),g(225),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,
     *tol,ii,jj,lit,ljt,mini,minj,maxi,maxj,iandk
INCLUDE(common/root)
      common/blkin/dxyz(4),gg(225),ft(225),dij(225),
     + pin(125),qin(125),rin(125),
     + ijx(225),ijy(225),ijz(225)
c mechanics
INCLUDE(common/modj)
INCLUDE(common/g80nb)
c mechanics
INCLUDE(common/psscrf)
INCLUDE(common/parallel)
      dimension q(*),coord(3,*),atmchg(*),iso(nshels,*)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      character*7 fnm
      character*6 snm
      data fnm,snm/"analf.m","nucatt"/
      data  m51/51/
      data dzero,pt5,done,two,three,five,seven /0.0d0,0.5d0,1.0d0,
     + 2.0d0,3.0d0,5.0d0,7.0d0/
      data rnine,eleven/9.0d0,11.0d0/
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
     +         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     +         21, 1, 1,16,16, 6, 1, 6, 1,11,
     +         11, 1,11, 6, 6/
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
_IF(parallel)
c***   **MPP**
      call pg_dlbreset
      next = ipg_dlbtask()
c***   **MPP**
_ENDIF
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      natoms = nat + itotbq
      l4 = natmod*3
      l5 = natmod
      ono = .false.
c     outv = oprint(59)
c     if (nprint.eq.-5) outv = .false.
      nav = lenwrd()
c
c     ----- determine required core -----
c
      i10  = 0
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
c ...
c mechanics
c ...
      if (oaminc) then
         ia40 = last
         ia50 = ia40 + l4
         ia60 = ia50 + l5
         last = ia60 + l5 + 1
      end if
c
c     ----- get core memory -----
c
      length = last - i10
      i10 = igmem_alloc_inf(length,fnm,snm,'i10',IGMEM_DEBUG)
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
c ...
c mechanics
c ...
      if (oaminc) then
         ia40 = last
         ia50 = ia40 + l4
         ia60 = ia50 + l5
         last = ia60 + l5 + 1
      end if
c
c
      dxyz(1) = enucf(natoms,atmchg,coord)
      do 450 i = 2 , 4
         dxyz(i) = 0.0d0
 450  continue
      do 460 i = 1 , natoms
         dxyz(2) = dxyz(2) + atmchg(i)*coord(1,i)
         dxyz(3) = dxyz(3) + atmchg(i)*coord(2,i)
         dxyz(4) = dxyz(4) + atmchg(i)*coord(3,i)
 460  continue
c
         if (oaminc) then
c     retrieve mechanics coordinates. xmod
            call secget(isect(472),m51,iblk172)
            call rdedx(q(ia40),l4,iblk172,idaf)
c                      retrieve nbact array. nbact
            call secget(isect(473),m51,iblk173)
            idum = (natmod-1) / nav + 1
            call rdedx(q(ia60),idum,iblk173,idaf)
c                      retrieve charges. chgmod
            call secget(isect(474),m51,iblk174)
            call rdedx(q(ia50),l5,iblk174,idaf)
         end if
c
c     ----- calculate -s- and -h0- matrices -----
c
c     - s- at x(i10)
c     -h0- at x(i20)
c
c        cpu = cpulft(1)
c
         tol = rln10*itol
         out = nprint.eq.3
         onorm = normf.ne.1 .or. normp.ne.1
         ndum = l2 + l2
         call vclr(q(i10),1,ndum)
c
c     ----- ishell
c
_IF(parallel)
c***   **MPP**
         do 440 ii = nshell,1,-1
_ELSE
         do 440 ii = 1 , nshell
_ENDIF
            i = katom(ii)
            icent = i
            pi = coord(1,i)
            qi = coord(2,i)
            ri = coord(3,i)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            lit = ktype(ii)
            mini = kmin(ii)
            maxi = kmax(ii)
            loci = kloc(ii) - mini
c
c     ----- jshell
c
_IF(parallel)
c***   **MPP**
            do 430 jj = ii, 1, -1
            icount_dlb = icount_dlb + 1
            if (icount_dlb.eq.next) then
c***   **MPP**
_ELSE
            do 430 jj = 1 , ii
_ENDIF
               j = katom(jj)
               jcent = j
               pj = coord(1,j)
               qj = coord(2,j)
               rj = coord(3,j)
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
                  mx = ix(i)
                  my = iy(i)
                  mz = iz(i)
                  if (oiandj) max = i
                  do 20 j = minj , max
                     ij = ij + 1
                     ijx(ij) = mx + jx(j)
                     ijy(ij) = my + jy(j)
                     ijz(ij) = mz + jz(j)
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
                  gg(i) = dzero
                  g(i) = dzero
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
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
                  cgi = cg(ig)
c
c     ----- j primitive
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
                        cpj = cp(jg)
                        cdj = cd(jg)
                        cfj = cf(jg)
                        cgj = cg(jg)
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
                           ni = i
                           do 230 j = 1 , ljt
                              jn = in + j
                              nj = j
                              call stvint
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
                              nj = j + 2
                              call stvint
                              pin(jn+25) = pint*t1
                              qin(jn+25) = qint*t1
                              rin(jn+25) = rint*t1
                              nj = j - 2
                              if (nj.gt.0) then
                                 call stvint
                              else
                                 pint = dzero
                                 qint = dzero
                                 rint = dzero
                              end if
                              n = (j-1)*(j-2)
                              dum = dfloat(n)*t2
                              pin(jn+50) = pint*dum
                              qin(jn+50) = qint*dum
                              rin(jn+50) = rint*dum
 230                       continue
 240                    continue
                        do 250 i = 1 , ij
                           mx = ijx(i)
                           my = ijy(i)
                           mz = ijz(i)
                           pyz = qin(my)*rin(mz)
                           dum = pyz*pin(mx)
                           dum1 = (pin(mx+25)+pin(mx+50))*pyz 
     +                     + (qin(my+25)+qin(my+50))*pin(mx)*rin(mz) 
     +                     + (rin(mz+25)+rin(mz+50))*pin(mx)*qin(my)
                           s(i) = s(i) + dij(i)*dum
                           g(i) = g(i) + dij(i)*(dum*aj*ft(i)+dum1)
                           gg(i) = gg(i) + dij(i)*(dum*aj*ft(i)+dum1)
 250                    continue
c
c     ----- nuclear attraction
c
                        dum = pi212*aa1
c                       facinv = aa/(fac*pi212)
                        do 260 i = 1 , ij
                           dij(i) = dij(i)*dum
 260                    continue
                        aax = aa*ax
                        aay = aa*ay
                        aaz = aa*az
                        do 320 ic = 1 , natoms
                           pnuc = -atmchg(ic)
                           cx = coord(1,ic)
                           cy = coord(2,ic)
                           cz = coord(3,ic)
                           pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
                           if (nroots.le.3) call rt123
                           if (nroots.eq.4) call roots4
                           if (nroots.eq.5) call roots5
                           mm = 0
                           do 290 k = 1 , nroots
                              uu = aa*u(k)
                              ww = w(k)*pnuc
                              tt = done/(aa+uu)
                              t = dsqrt(tt)
                              p0 = (aax+uu*cx)*tt
                              q0 = (aay+uu*cy)*tt
                              r0 = (aaz+uu*cz)*tt
                              in = -5 + mm
                              do 280 i = 1 , lit
                                 in = in + 5
                                 ni = i
                                 do 270 j = 1 , ljt
                                    jn = in + j
                                    nj = j
                                    call stvint
                                    pin(jn) = pint
                                    qin(jn) = qint
                                    rin(jn) = rint*ww
 270                             continue
 280                          continue
                              mm = mm + 25
 290                       continue
                           do 310 i = 1 , ij
                              mx = ijx(i)
                              my = ijy(i)
                              mz = ijz(i)
                              dum = dzero
                              mm = 0
                              do 300 k = 1 , nroots
                                 dum = dum + pin(mx+mm)*qin(my+mm)
     +                                 *rin(mz+mm)
                                 mm = mm + 25
 300                          continue
                              g(i) = g(i) + dum*dij(i)
 310                       continue
 320                    continue
                        if (oaminc) then
c ...
c    skip this code if dealing with a dummy atom.
c ...
                           ic3 = -4
                           do 380 icmod = 1 , natmod
                              ic3 = ic3 + 3
                              call skip80(icent,jcent,icmod,q(ia60),
     +                           oskmod)
                              if (.not.(oskmod)) then
                                 ia40p = ia40 + ic3
                                 pnuc = -q(ia50+icmod-1)
                                 cx = q(ia40p+1)
                                 cy = q(ia40p+2)
                                 cz = q(ia40p+3)
                                 pp = aa*((ax-cx)**2+(ay-cy)**2+(az-cz)
     +                                **2)
                                 if (nroots.le.3) call rt123
                                 if (nroots.eq.4) call roots4
                                 if (nroots.eq.5) call roots5
                                 mm = 0
                                 do 350 kmod = 1 , nroots
                                    uu = aa*u(kmod)
                                    ww = w(kmod)*pnuc
                                    tt = done/(aa+uu)
                                    t = dsqrt(tt)
                                    p0 = (aax+uu*cx)*tt
                                    q0 = (aay+uu*cy)*tt
                                    r0 = (aaz+uu*cz)*tt
                                    in = -5 + mm
                                    do 340 imod = 1 , lit
                                       in = in + 5
                                       ni = imod
                                       do 330 jmod = 1 , ljt
                                         jn = in + jmod
                                         nj = jmod
                                         call stvint
                                         pin(jn) = pint
                                         qin(jn) = qint
                                         rin(jn) = rint*ww
 330                                   continue
 340                                continue
                                    mm = mm + 25
 350                             continue
                                 do 370 imod = 1 , ij
                                    mx = ijx(imod)
                                    my = ijy(imod)
                                    mz = ijz(imod)
                                    dum = 0.0d0
                                    mm = 0
                                    do 360 kmod = 1 , nroots
                                       dum = dum + pin(mx+mm)*qin(my+mm)
     +                                    *rin(mz+mm)
                                       mm = mm + 25
 360                                continue
                                    g(imod) = g(imod) + dum*dij(imod)
 370                             continue
                              end if
 380                       continue
                        end if
                     end if
c ...
c ...
 390              continue
 400           continue
c
c
c     - s- at x(i10)
c     -h0- at x(i20)
c     -t-  at x(i30)
c
c
c     ----- set up overlap and h-core matrices
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
                     q(jn-1+i10) = s(nn)
                     q(jn-1+i20) = g(nn)
                     q(jn-1+i30) = gg(nn)
 410              continue
 420           continue
_IF(parallel)
              next = ipg_dlbtask()
              endif
_ENDIF
 430        continue
 440     continue
c
_IF(parallel)
c***   integrals have been calculated partly on each node
c***   now gather them and send them to each other
c
      call pg_dgop(101,q(i10),l2,'+')
      call pg_dgop(102,q(i20),l2,'+')
      call pg_dgop(103,q(i30),l2,'+')
      call pg_dlbpush
c
c***   **MPP**
_ENDIF
         if (lpseud.eq.1) then
c
c ---- if ecp pseudopotentials required, call ecpint
c
            call vclr(q(i40),1,l2)
            call ecpint(q,q(i40),q(i20),iso,nshell,ono,out)
         endif
         if (lpseud.eq.2) then
c
c ---- if nonlocal pseudopotentials required, call xpsnlc
c
            call xpsnlc(q(i20),q(i50),q(i60),q(i40),q(i70),q(i80),
     +           iso,l2,non,num,nshell,ono)
         end if
c
c
c     ----- output 1-electron integrals to dumpfile
c
         call sec192(q(i10),q(i20),q(i30),l2)
c
c     ----- now compute dipole moment integrals
c
_IF(parallel)
         call vclr(q(i10),1,l2*3)
_ENDIF
         call dipxyz(q(i10),q(i20),q(i30))
c
c
_IF(parallel)
c
c***   ***MPP***
c
c      combine and send to each other
c
      call pg_dgop(104,q(i10),l2,'+')
      call pg_dgop(105,q(i20),l2,'+')
      call pg_dgop(106,q(i30),l2,'+')
c
c***   **MPP**
_ENDIF
c     ---- load to ionsec
c
         call secdip(q(i10),q(i20),q(i30),l2)
c
c     ---- print if requested
c
_IF(parallel)
c***  **MPP**
      if (out.and.opg_root()) call prt1e(q(i10),ionsec)
c***   **MPP**
_ELSE
      if (out) call prt1e(q(i10),ionsec)
_ENDIF
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,fnm,snm,'i10')
      return
      end
      subroutine potgen(q,coord,atmchg,cnrmpt,efgmol,iw)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c     common for psscrf
INCLUDE(common/psscrf)
INCLUDE(common/infoa)
INCLUDE(common/timez)
c
      character*10 charwall
      dimension q(*),coord(3,*),atmchg(*),cnrmpt(3,*),efgmol(2,*)
c
c     driver for psscrf analysis
c
      write(iw,89003) cpulft(1),charwall()
89003 format(/1x,'entered psscrf potentials module at ',f8.2,
     +  ' seconds',a10,' wall')
c     generate potential surface
      call dtgrid (q,coord,atmchg,cnrmpt,efgmol)
c     if(tim.ge.timlim)go to 3
c
      call clredx
      top=cpulft(1)
      write(iw,24)top,charwall()
 24   format(/1x,
     * 'end of potential surface calculation at ',f8.2,' seconds',
     * a10,' wall')
c
      return
      end
      subroutine potini(core,coord,atmchg,deltas,cnrmpt,efgmol,esp)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
INCLUDE(common/iofile)
INCLUDE(common/statis)
INCLUDE(common/psscrf)
INCLUDE(common/runlab)
INCLUDE(common/timez)
INCLUDE(common/limy)
INCLUDE(common/scra7)
INCLUDE(common/infoa)
c
      dimension core(*),coord(3,*),atmchg(*),deltas(2,*),cnrmpt(3,*),
     + efgmol(2,*),esp(*)
      character*7 fnm
      character*6 snm
      data fnm,snm/"analf.m","potini"/
c
c ----- allocate core storage, assuming all available core
c ----- is to be allocated
c
      l2 = num*(num+1)/2
      lwordp = igmem_max_memory() 
     +       - memreq_pg_dgop(max(l2,2*(nat+itotbq)),'+')
     +       - 2*igmem_overhead()
      ibase = igmem_alloc_inf(lwordp,fnm,snm,'ibase',IGMEM_DEBUG)
c
c -----  calculate esp surface
c
      call potgen(core(ibase),coord,atmchg,cnrmpt,efgmol,iwr)
c
      call splzbq(core(ibase),coord,atmchg,deltas,cnrmpt,efgmol,
     +            esp,iwr)
c
      call timana(14)
c
c ----- reset core allocation
c
      call gmem_free_inf(ibase,fnm,snm,'ibase')
c
      return
      end
      subroutine psanin
c
c     this is data input routine for graphical analysis
c     primitive checking only (now a null routine)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/direc)
INCLUDE(common/machin)
INCLUDE(common/scra7)
c
INCLUDE(common/psscrf)
c
c     default settings
c
      inkblk=ibl7la
c
      return
      end
       subroutine splzbq(ptchg,coord,atmchg,deltas,cnrmpt,efgmol,
     +                   esp,iw)
        implicit REAL  (a-h,p-w),integer(i-n),logical(o)
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/prints)
INCLUDE(common/psscrf)
INCLUDE(common/vdwrad)
       common/surfgn/inumbq(maxat),irejec(maxat),cavszd(maxat),
     +       nocent,maxit,damper
c
       dimension ptchg(*),esp(*)
       dimension coord(3,*),atmchg(*),deltas(2,*),cnrmpt(3,*),
     +           efgmol(2,*)
c
        data pi/3.14159265359d0/
        data p1,p2,p4,ps1/1.0d0,2.0d0,4.0d0,-1.0d0/
c       data p32/1.5d0/
c
       thresh = 1.0d-5
       totold=0.0d0
       shpcor=0.0d0
       spzcor=0.0d0
       corfac=(delect-p1)/(p4*pi*delect)
       ist=nat+1
       ifn=nat+itotbq
c
       do 1 i=ist,ifn
          esp(i)=efgmol(1,i)
          efgmol(1,i)=(efgmol(1,i)-efgmol(2,i))/facnrm
          efgmol(1,i)=corfac*efgmol(1,i)*deltas(1,i)
          ptchg(i)=efgmol(1,i)
          atmchg(i)=0.0d0
1         continue
c
       idaz=nat
       write(iw,1000)
1000   format(1x,40('*')/
     +        1x,'*****  initial charge distribution *****')
       do 21 i=1,nat
       sum=0.0d0
       do 22 j=1,inumbq(i)
           idaz=idaz+1
           sum=sum+ptchg(idaz)
22       continue
       write(iw,1001)i,sum
1001       format(1x,'***** charge on sphere ',i2,' = ',f8.6,'****')
21       continue
       write(iw,1002)
1002       format(1x,40('*')/
     +            1x,10('*'),'      Iterating     ',10('*'))
       do 2 i=1,maxit
          do 3 ii=ist,ifn
              atmchg(ii)=damper*ptchg(ii)+(p1-damper)*atmchg(ii)
3          continue
          do 4 j=ist,ifn
          ptesp1=0.0d0
          ptesp2=0.0d0
          ptchg(j)=0.0d0
              do 5 k=ist,ifn
                 ptefg=0.0d0
                 if(k.ne.j) then
              distjk= (coord(1,j)-coord(1,k))**2
     +               +(coord(2,j)-coord(2,k))**2
     +               +(coord(3,j)-coord(3,k))**2
              distjk=dsqrt(distjk)
              ptesp1=p1/distjk
c
              distjk= (cnrmpt(1,j)-coord(1,k))**2
     +               +(cnrmpt(2,j)-coord(2,k))**2
     +               +(cnrmpt(3,j)-coord(3,k))**2
              distjk=dsqrt(distjk)
              ptesp2=p1/distjk
c
              ptefg=atmchg(k)*(ptesp2-ptesp1)
         if(((ptefg.lt.0.0d0).and.(atmchg(k).gt.0.0d0)).or.
     +     ((ptefg.gt.0.0d0).and.(atmchg(k).lt.0.0d0)))
     +     ptefg=ps1*ptefg
c
                ptchg(j)=ptchg(j)+ptefg
                  endif
5               continue
c
              shpcor=(deltas(1,j)/facnrm)*ptchg(j)
              spzcor=dsqrt(deltas(1,j)/(p4*pi*(deltas(2,j)**2)))
              ptchg(j)=efgmol(1,j)-corfac*(shpcor
     +              -p2*pi*atmchg(j)*(p1-spzcor))
c
                daz1=dabs(ptchg(j)/atmchg(j))
                if(daz1.gt.5.0d0)ptchg(j)=atmchg(j)
4          continue
c
          totchg=0.0d0
          sumneg=0.0d0
          sumpos=0.0d0
          do 6 j=ist,ifn
              if(ptchg(j).lt.0.0d0)sumneg=sumneg+ptchg(j)
              if(ptchg(j).gt.0.0d0)sumpos=sumpos+ptchg(j)
              totchg=totchg+ptchg(j)
6          continue
           write(iw,1003)i,totchg
1003       format(1x,'**** cycle',i2,' total charge =',f8.6,
     *            ' ****')
           if(dabs(totchg-totold).le.thresh)goto 999
           totold=totchg
c
2       continue
c
999       write(iw,1004)
1004       format(40('*'))
          write(iw,101)sumpos
        write(iw,102)sumneg
101       format(1x,'*** total positive charge = ',f9.6,' ***')
102       format(1x,'*** total negative charge = ',f9.6,' ***')
        totchg=totchg+ich
        if ((sumpos.lt.thresh).or.(sumneg.gt.(ps1*thresh)))
     *  totchg=p2*totchg
       do 7 i=ist,ifn
           atmchg(i)=ptchg(i)
          if (atmchg(i).gt.0.0d0)then
              atmchg(i)=atmchg(i)*(p1-totchg/(p2*sumpos))
          else
              atmchg(i)=atmchg(i)*(p1-totchg/(p2*sumneg))
          endif
7       continue
       totchg=0.0d0
       do 9 i=ist,ifn
          totchg=totchg+atmchg(i)
9       continue
       write(iw,1006)totchg,itotbq
1006   format(1x,
     +     '***** total charge after compensation =',f10.6,'*****'/
     +  1x,'***** number of surface points        =',i10,'*****')
c
       if(oprint(54)) then
        write(iw,1005)itotbq
1005    format(1x,40('*')/
     +  1x,'**** charges on',i5,' surface points ****'/1x,40('*'))
        do 8 i=ist,ifn
        write(iw,103)coord(1,i),coord(2,i),coord(3,i),atmchg(i)
103     format(' bq  ',4(f10.6))
8       continue
       endif
c
       return
      end
      subroutine surgen(q,c,czan,delta,cnrmpt,cavsze,nocycl,iw)
c
      implicit REAL     (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/segm)
      dimension q(*),c(3,*),czan(*),delta(2,*),cnrmpt(3,*),
     +          cavsze(*)
      character*7 fnm
      character*6 snm
      data fnm,snm/"analf.m","surgen"/
c
c ===
c === this is a conservative estimate of memory requirements
c === in surfgen .. avoids overwriting of c,deltas and cnrmpt arrays
c
      lwor9=maxat*20
      lword= 8 * lwor9
      i10 = igmem_alloc_inf(lword,fnm,snm,'i10',IGMEM_DEBUG)
c
      i20 = i10 + 3*lwor9
      i30 = i20 + 3*lwor9
      call surgn2(c,czan,delta,cnrmpt,cavsze,
     +  q(i10),q(i20),q(i30),lwor9,nocycl,iw)
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
      return
      end
      subroutine surgn2(cc,cczan,deltas,cnrmpt,cavsze,
     +                  qc,qnrmpt,qdelt,lwor9,nocycl,iw)
c
      implicit REAL  (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z), character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/runlab)
INCLUDE(common/phycon)
c -----   common blocks for psscrf -----
INCLUDE(common/psscrf)
INCLUDE(common/vdwrad)
      common/surfgn/numbq(maxat),irejec(maxat),cavszd(maxat),
     *  nocent,itmax,damper
c -----   common blocks for psscrf -----
c
      dimension qc(3,*),qnrmpt(3,*),qdelt(2,*)
      dimension cc(3,*),cczan(*),cnrmpt(3,*),deltas(2,*),
     + cavsze(*)
c
      data pi/3.14159265359d0/
      data reject/-9999.0d0/
      data cutoff/0.01d0/
c
      non=nat
c
      call flushn(iw)
      write(iw,1005)
1005  format(1x,40('*')/
     + 1x,'****** Surface generation routine ******')
      if(nocycl.eq.1)then
      if(nocent.eq.0)then
      write(iw,1004)
1004  format(1x,'*** using default Van der Waals radii **')
      do 2 i=1,non
         if (radvdw(imass(i)).le.0.001d0) then
            write(iwr,*) '*** no VDW radius for element',imass(i)
            write(iwr,*) '*** 4.0 is used; adapt the program or use own'
            radvdw(imass(i)) = 4.0d0
         end if
         cavsze(i)=radvdw(imass(i))
2      continue
      else
      write(iw,1006)
1006  format(1x,'****** using user supplied radii *******')
      do 21 i=1,nocent
       cavsze(i)=cavszd(i)/toang(1)
21    continue
      endif
      endif
c
       iangno=360.0d0/ptspac
       facnrm=0.0d0
c
      do 3 iptyp=1,2
       if(iptyp.eq.2) facnrm=0.01d0
       icnt=non
       do 4 i=1,non
          radius=cavsze(i)-facnrm
          numbq(i)=0
        do 5 j=1,(iangno/2)
         theta=j*ptspac-(ptspac/2.0d0)
         theta=(theta/180.0d0)*pi
         do 6 k=0,iangno/4
              fi=k*ptspac
              fi=(fi/180.0d0)*pi
c
              goto(111,222),iptyp
c
c    surface points generated
c
111           continue
c
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) 'extend core 8'
                 go to 1009
              end if
              qc(1,icnt)=radius*dcos(fi)*dsin(theta)
              qc(2,icnt)=radius*dsin(fi)*dsin(theta)
              qc(3,icnt)=radius*dcos(theta)
c
              qdelt(1,icnt)=dabs(2.0d0*(cavsze(i)**2)*dsin(theta)*
     +         dsin((ptspac*pi)/360.0d0)*ptspac*(pi/180.0d0))
              qdelt(2,icnt)=cavsze(i)
c
              numbq(i)=numbq(i)+1
c
c    as symmetrical about z-axis:
c
              ix=idint(qc(1,icnt)*100.0d0)
              if(ix.eq.0) goto 61
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) 'extend core 7'
                 go to 1009
              end if
              qc(1,icnt)=-1.0d0*qc(1,icnt-1)
              qc(2,icnt)=qc(2,icnt-1)
              qc(3,icnt)=qc(3,icnt-1)
c
              qdelt(1,icnt)=dabs(2.0d0*(cavsze(i)**2)*dsin(theta)*
     +         dsin((ptspac*pi)/360.0d0)*ptspac*(pi/180.0d0))
              qdelt(2,icnt)=cavsze(i)
c
              numbq(i)=numbq(i)+1
c
c    and symmetrical about x axis:
c
61            iy=idint(qc(2,icnt)*100.0d0)
              if(iy.eq.0) goto 6
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core 6'
                 go to 1009
              end if
              qc(1,icnt)=-1.0d0*qc(1,icnt-1)
              qc(2,icnt)=-1.0d0*qc(2,icnt-1)
              qc(3,icnt)=qc(3,icnt-1)
c
              qdelt(1,icnt)=dabs(2.0d0*(cavsze(i)**2)*dsin(theta)*
     +         dsin((ptspac*pi)/360.0d0)*ptspac*(pi/180.0d0))
              qdelt(2,icnt)=cavsze(i)
c
              numbq(i)=numbq(i)+1
c
              if(ix.eq.0) goto 6
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core 5'
                 go to 1009
              end if
              qc(1,icnt)=-1.0d0*qc(1,icnt-1)
              qc(2,icnt)=qc(2,icnt-1)
              qc(3,icnt)=qc(3,icnt-3)
c
              qdelt(1,icnt)=dabs(2.0d0*(cavsze(i)**2)*dsin(theta)*
     +         dsin((ptspac*pi)/360.0d0)*ptspac*(pi/180.0d0))
              qdelt(2,icnt)=cavsze(i)
c
              numbq(i)=numbq(i)+1
              goto 6
c    points on normal to surface points generated
c
c
222           continue
c
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core 4'
                 go to 1009
              end if             
              qnrmpt(1,icnt)=radius*dcos(fi)*dsin(theta)
              qnrmpt(2,icnt)=radius*dsin(fi)*dsin(theta)
              qnrmpt(3,icnt)=radius*dcos(theta)
              numbq(i)=numbq(i)+1
c
c    as symmetrical about z-axis:
c
              ix=idint(qnrmpt(1,icnt)*100.0d0)
              if(ix.eq.0) goto 62
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core 3'
                 go to 1009
              end if
              qnrmpt(1,icnt)=-1.d0*qnrmpt(1,icnt-1)
              qnrmpt(2,icnt)=qnrmpt(2,icnt-1)
              qnrmpt(3,icnt)=qnrmpt(3,icnt-1)
c
              numbq(i)=numbq(i)+1
c
c    and symmetrical about x axis:
c
62            iy=idint(qnrmpt(2,icnt)*100.d0)
              if(iy.eq.0) goto 6
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core 2'
                 go to 1009
              end if
              qnrmpt(1,icnt)=-1.0d0*qnrmpt(1,icnt-1)
              qnrmpt(2,icnt)=-1.0d0*qnrmpt(2,icnt-1)
              qnrmpt(3,icnt)=qnrmpt(3,icnt-1)
              numbq(i)=numbq(i)+1
c
              if(ix.eq.0) goto 6
              icnt=icnt+1
              if(icnt.gt.lwor9) then
                 write(iwr,*) ' extend core  1'
                 go to 1009
              end if
              qnrmpt(1,icnt)=-1.0d0*qnrmpt(1,icnt-1)
              qnrmpt(2,icnt)=qnrmpt(2,icnt-1)
              qnrmpt(3,icnt)=qnrmpt(3,icnt-1)
              numbq(i)=numbq(i)+1
c
6             continue
5        continue
4      continue
3     continue
c
c     surface points now generated
c     held in array of type (xcord,atom)
c                           (ycord,atom)
c                           (zcord,atom)
c     number of surface points coords per atom=numbq(atom)
c
c     now convert them into coordinate system of atoms
c
      icnt=non
      do 7 i=1,non
         do 8 j=1, numbq(i)
            icnt=icnt+1
            if(icnt.gt.lwor9) then
               write(iwr,*) ' extend core do 8'
               go to 1009
            end if
            qc(1,icnt)  =qc(1,icnt)  +c(1,i)
            qc(2,icnt)  =qc(2,icnt)  +c(2,i)
            qc(3,icnt)  =qc(3,icnt)  +c(3,i)
c
            qnrmpt(1,icnt)=qnrmpt(1,icnt)+c(1,i)
            qnrmpt(2,icnt)=qnrmpt(2,icnt)+c(2,i)
            qnrmpt(3,icnt)=qnrmpt(3,icnt)+c(3,i)
8        continue
7      continue
c
c     now test to see if they overlap with each other
c
c
      do 1112 i=1,non
         irejec(i)=0
1112  continue
      do 9 i=1,(non-1)
           ist=(i-1)*numbq(i) + non
           do 10 j=(i+1),non
                jst=(j-1)*numbq(i) + non
                do 11 k=1,numbq(i)
                      icnt=ist+k
c     dist=distance of kth surface point on atom i
c           from atom j
                     sqdst= (qc(1,icnt)-c(1,j))**2
     *                     +(qc(2,icnt)-c(2,j))**2
     *                     +(qc(3,icnt)-c(3,j))**2
                     dist=dsqrt(sqdst)
                     if (dist.lt.(cavsze(j)+cutoff)) then
                        qc(1,icnt)  =reject
                        qc(2,icnt)  =reject
                        qc(3,icnt)  =reject
                        irejec(i)=irejec(i)+1
                     endif
11              continue
c
                do 12 l=1,numbq(j)
                     jcnt=jst+l
c     dist=distance of lth surface point on atom j
c           from atom i
                     sqdst= (qc(1,jcnt)-c(1,i))**2
     *                     +(qc(2,jcnt)-c(2,i))**2
     *                     +(qc(3,jcnt)-c(3,i))**2
                     dist=dsqrt(sqdst)
                     if (dist.lt.(cavsze(i)+cutoff)) then
                          qc(1,jcnt)=reject
                          qc(2,jcnt)=reject
                          qc(3,jcnt)=reject
                          irejec(j)=irejec(j)+1
                     endif
12              continue
10       continue
9     continue
c
c     now remove bad surface points and contract arrays accordingly
c
      icnt=non
      idas=non
      itotbq=0
      do 13 i=1,non
      cc(1,i)=c(1,i)
      cc(2,i)=c(2,i)
      cc(3,i)=c(3,i)
      cczan(i)=czan(i)
            do 14 j=1,numbq(i)
                  icnt=icnt+1
                  if (qc(1,icnt).ne.reject) then
                     idas=idas+1
                     if(idas.gt.maxat) then
                        write(iwr,*) 'idas =',idas,' maxat =',maxat,
     1                               ' adapt sizes'
                        go to 1009
                     end if
                     cc(1,idas)=qc(1,icnt)
                     cc(2,idas)=qc(2,icnt)
                     cc(3,idas)=qc(3,icnt)
                     deltas(1,idas)=qdelt(1,icnt)
                     deltas(2,idas)=qdelt(2,icnt)
c
                     cnrmpt(1,idas)=qnrmpt(1,icnt)
                     cnrmpt(2,idas)=qnrmpt(2,icnt)
                     cnrmpt(3,idas)=qnrmpt(3,icnt)
c
                  endif
14           continue
             numbq(i)=idas-non-itotbq
             itotbq=idas-non
13    continue
c
c
c     finished
c
       write(iw,1002)itotbq
1002   format(1x,40('*')/
     + 1x,'**** number of surface points = ',i4,' ***')
       do 1000 i=1,non
       write(iw,1001)i,numbq(i)
1001   format(1x,'**** ','no. of points on atom ',i2,' = ',i3,' ****')
1000       continue
       write(iw,1003)
1003   format(1x,40('*')/)
c
c ----- add surface points to gamess deck as bq coords
c
      non = nat
      do 15 i=1,itotbq
         non=non+1
c        zaname(nat)='bq'
c        imass(nat)=0
c        amas(nat)=0.0d0
         czan(non)=0.0d0
15    continue
       call flushn(iw)
      if(nat.le.maxat) return
      write(iwr,*) 'nat =',nat,' maxat =',maxat,' adapt sizes'
1009  call caserr(
     + 'excessive number of points in surface generator')
c
      return
      end
       subroutine tofren(etot,etot1,atmchg,esp)
       implicit REAL  (a-h,p-w),integer(i-n),logical(o)
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/psscrf)
INCLUDE(common/vdwrad)
       dimension atmchg(*),esp(*)
c
       data pt5/0.5d0/
c
c       convert hf energy to free energy of solvation
c
       ist=nat+1
       ifn=nat+itotbq
       sum=0.0d0
       do 1 i=ist,ifn
          esp(i)=atmchg(i)*esp(i)
          sum=sum+esp(i)
1       continue
c
       etot1=etot-pt5*sum
c
       return
      end
      subroutine ver_analf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/analf.m,v $
     +     "/
      data revision /"$Revision: 6219 $"/
      data date /"$Date: 2010-12-20 17:30:03 +0100 (Mon, 20 Dec 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
