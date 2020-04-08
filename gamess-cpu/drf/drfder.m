      subroutine drfsint(iatom,dsx,dsy,dsz)
c
c this is an old routine (gamess 91) with new includes ....
c sf 02-97
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
_IF(parallel)
INCLUDE(../m4/common/nodinf)
_ENDIF
INCLUDE(../m4/common/timez)
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
c
c     common/nshel_z/ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
c    +               cf(mxprim),cg(mxprim),
c    +       kstart(mxshel),katom(mxshel),ktype(mxshel),
c    +       kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
c    +               nshell,non,numorb,ndumm
c

      common/intdip/pint,qint,rint,pintd,qintd,rintd,
     1t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,cx,cy,cz
c     common/junk/desp(3,maxat),
c    *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj
c    + ,cx,cy,cz,
c    +  pin(25),qin(25),rin(25),pd(25),qd(25),rd(25),dij(100),
c    +  ijx(100),ijy(100),ijz(100),sx(100),sy(100),sz(100)
      dimension pin(25),qin(25),rin(25),
     +  pd(25),qd(25),rd(25),dij(100),
     +  ijx(100),ijy(100),ijz(100),sx(100),sy(100),sz(100)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/segm)
      common/blkin/gout(510),nword
INCLUDE(../m4/common/misc)
c
c
c
c     dimension sz1(*),sx1(*),sy1(*)
      dimension dsx(*), dsy(*), dsz(*)
      data sqrt3 /1.73205080756888d0/
      data dzero,done /0.0d0,1.0d0/
      data rln10 /2.30258d0/
      data zgrhf,zgvb/'grhf','gvb'/
      data m511/511/
c
c     ----- calculate derivatives of the overlap matrix -----
c
      tol = rln10*itol
      out = nprint.eq. - 3 .or. nprint.eq. - 10
      onorm = normf.ne.1 .or. normp.ne.1
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      nword = 1
c
c     ----- ishell
c
      do 110 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 110
         pi = c(1,iat)
         qi = c(2,iat)
         ri = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit1 = lit + 1
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
c        do 100 jj = 1 , ii
         do 100 jj = 1 , nshell
            jat = katom(jj)
c           if (iat.ne.jat) then
c  sf on all atoms!
c
               pj = c(1,jat)
               qj = c(2,jat)
               rj = c(3,jat)
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
               call idxadrf(ijx,ijy,ijz,ij,
     +           mini,maxi,minj,maxj,.false.,5,
     +                    1,1)
               do 20 i = 1 , ij
                  sx(i) = dzero
                  sy(i) = dzero
                  sz(i) = dzero
 20            continue
c
c     ----- i primitive
c
               do 70 ig = i1 , i2
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
                  do 60 jg = j1 , j2
                     aj = ex(jg)
                     aa = ai + aj
                     aa1 = done/aa
                     dum = aj*arri*aa1
                     if (dum.le.tol) then
                        fac = dexp(-dum)
                        csj = cs(jg)*fac
                        cpj = cp(jg)*fac
                        cdj = cd(jg)*fac
                        cfj = cf(jg)*fac
                        cgj = cg(jg)*fac
                        ax = (axi+aj*pj)*aa1
                        ay = (ayi+aj*qj)*aa1
                        az = (azi+aj*rj)*aa1
c
c     ----- density factor
c
                        call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                                  csj,cpj,cdj,cfj,cgj,
     +                              mini,maxi,minj,maxj,.false.,.false.,
     +                              onorm)
c
c     ----- overlap
c
                        t = dsqrt(aa1)
                        p0 = ax
                        q0 = ay
                        r0 = az
                        in = -5

                        do 40 i = 1 , lit1
                           in = in + 5
                           ni = i
                           do 30 j = 1 , ljt
                              jn = in + j
                              nj = j
c                             call vint
                              call vintdrf
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
 30                        continue
 40                     continue
                        call oneld(pin,qin,rin,pd,qd,rd,ai,lit,ljt,1,5)

c
                        do 50 i = 1 , ij
                           mx = ijx(i)
                           my = ijy(i)
                           mz = ijz(i)
                         sx(i) = sx(i) + dij(i)*pd(mx)*qin(my)*rin(mz)
                         sy(i) = sy(i) + dij(i)*pin(mx)*qd(my)*rin(mz)
                         sz(i) = sz(i) + dij(i)*pin(mx)*qin(my)*rd(mz)
 50                     continue
                     end if
c
c     ----- end of primitive loops -----
c
 60               continue
 70            continue
c
c     ----- calculate derivatives of overlap matrix -----
c
               n = 0
               do 90 i = mini , maxi
                  in = loci + i
                  do 80 j = minj , maxj
                     n = n + 1
                     jn = locj + j
c                    if (jn.le.in) then
c                       nn = iky(in) + jn
                        nn = iky(max(in,jn)) + min(in,jn)
                        dsx(nn) = dsx(nn) + sx(n)
                        dsy(nn) = dsy(nn) + sy(n)
                        dsz(nn) = dsz(nn) + sz(n)
c                       sx1(nn)=sx(n)
c                       sy1(nn)=sy(n)
c                       sz1(nn)=sz(n)
c                    end if
 80               continue
 90            continue
c           end if
c idem
 100     continue
 110  continue
c
      return
 6010 format (/40x,'lagrangian weighted density'/40x,27('*'))
 6020 format (//5x,'i',4x,'j',15x,'sx',18x,'sy',18x,'sz',18x,'lij')
 6030 format (1x,2i5,5x,3f20.8,i4)
      end
_IF()
      subroutine sdummy(dsx,dsy,dsz)
c
      implicit REAL  (a-h,o-z)
INCLUDE(../m4/common/sizes)
      logical out,norm
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
INCLUDE(../m4/common/common)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/mapper)
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     4   ijx(100),ijy(100),ijz(100)
c     dimension dd(*)
c
INCLUDE(../m4/common/prnprn)
c     dimension comp(3),dipd(3,3,maxat),dipi(3,3,maxat),
c    &  dipn(3,3,maxat)
      dimension comp(3)
      dimension dsx(nx),dsy(nx),dsz(nx)
      dimension sx(100),sy(100),sz(100)
c
      character*5 comp
      data comp/'d/dx','d/dy','d/dz'/
      data ndim/5/
      data dzero,one /0.0d0,1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
c     call onepdm(dd,dd(nx+1))
      tol = 2.30258d0*itol
      out = odebug(15)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
c     do 50 n = 1 , nat
c        do 30 i = 1 , 3
c           do 20 j = 1 , 3
c              dipn(i,j,n) = 0.0d0
c              dipi(i,j,n) = 0.0d0
c              dipd(i,j,n) = 0.0d0
c20         continue
c30      continue
c        do 40 i = 1 , 3
c           dipn(i,i,n) = czan(n)
c40      continue
c50   continue
c     ----- ishell
      do 130 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 130
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit1 = lit + 1
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 120 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,
     +                  .false.,ndim,1,1)
            dxx = 0.0d0
            dyy = 0.0d0
            dzz = 0.0d0
            dxy = 0.0d0
            dxz = 0.0d0
            dyx = 0.0d0
            dyz = 0.0d0
            dzx = 0.0d0
            dzy = 0.0d0
            do i = 1, ij
              sx(i) = dzero
              sy(i) = dzero
              sz(i) = dzero
            enddo
c     ----- i primitive
            do 110 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 100 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     in = -ndim
                     do 70 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 60 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call vintdrf
                           xin(jn) = xint0*t
                           yin(jn) = yint0*t
                           zin(jn) = zint0*t
c                          xdip(jn) = xintx*tinv
c                          ydip(jn) = yinty*tinv
c                          zdip(jn) = zintz*tinv
 60                     continue
 70                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
c                    call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
c    +                          ljt,1,ndim)
c     ----- calculate derivatives of dipole matrix -----
                        do 50 i = 1 , ij
                           mx = ijx(i)
                           my = ijy(i)
                           mz = ijz(i)
                           sx(i) = sx(i) 
     1    + dij(i)*xd(mx)*yin(my)*zin(mz)
                           sy(i) = sy(i) 
     1    + dij(i)*xin(mx)*yd(my)*zin(mz)
                           sz(i) = sz(i) 
     1    + dij(i)*xin(mx)*yin(my)*zd(mz)
 50                     continue
                     n = 0
                     do 90 i = mini , maxi
                        in = loci + i
                        do 80 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
c                          dum = dd(nn)*dij(n)
                           dum = dij(n)
                           dum = dum + dum
c                          nnx = ijx(n)
c                          ny = ijy(n)
c                          nz = ijz(n)
                           dsx(n) = dsx(n) + dum * sx(nn)
                           dsy(n) = dsy(n) + dum * sy(nn)
                           dsz(n) = dsz(n) + dum * sz(nn)
c                          dsx(nn) = dum*xd(nnx)*yin(ny)*zin(nz)
c                          dsy(nn) = dum*xin(nnx)*yd(ny)*zin(nz)
c                          dsz(nn) = dum*xin(nnx)*yin(ny)*zd(nz)
c                          dxx = dxx + dum*xdipd(nnx)*yin(ny)*zin(nz)
c                          dyy = dyy + dum*xin(nnx)*ydipd(ny)*zin(nz)
c                          dzz = dzz + dum*xin(nnx)*yin(ny)*zdipd(nz)
c                          dxy = dxy + dum*xdip(nnx)*yd(ny)*zin(nz)
c                          dyx = dyx + dum*xd(nnx)*ydip(ny)*zin(nz)
c                          dxz = dxz + dum*xdip(nnx)*yin(ny)*zd(nz)
c                          dzx = dzx + dum*xd(nnx)*yin(ny)*zdip(nz)
c                          dyz = dyz + dum*xin(nnx)*ydip(ny)*zd(nz)
c                          dzy = dzy + dum*xin(nnx)*yd(ny)*zdip(nz)
 80                     continue
 90                  continue
                  end if
 100           continue
 110        continue
c           dipi(1,1,iat) = dipi(1,1,iat) - dxx
c           dipi(2,2,iat) = dipi(2,2,iat) - dyy
c           dipi(3,3,iat) = dipi(3,3,iat) - dzz
c           dipi(1,2,iat) = dipi(1,2,iat) - dxy
c           dipi(1,3,iat) = dipi(1,3,iat) - dxz
c           dipi(2,1,iat) = dipi(2,1,iat) - dyx
c           dipi(3,1,iat) = dipi(3,1,iat) - dzx
c           dipi(2,3,iat) = dipi(2,3,iat) - dyz
c           dipi(3,2,iat) = dipi(3,2,iat) - dzy
 120     continue
 130  continue
c
c     if (out) then
c        write (iwr,6010)
c        do 150 n = 1 , nat
c           write (iwr,6020)
c           do 140 nc = 1 , 3
c              write (iwr,6030) zaname(n) , comp(nc) ,
c    +                         (dipi(nn,nc,n),nn=1,3)
c140        continue
c150     continue
c     end if
c     do 180 i = 1 , 3
c        do 170 j = 1 , 3
c           do 160 k = 1 , nat
c              dipd(i,j,k) = dipd(i,j,k) + dipn(i,j,k) + dipi(i,j,k)
c160        continue
c170     continue
c180  continue
      return
 6010 format (//35x,'integral derivative contribution'//30x,'x',15x,'y',
     +        15x,'z',/)
 6020 format (//)
 6030 format (5x,a8,5x,a5,3f16.8)
      end
_ENDIF
      subroutine drfdint(iatom,dipxx,dipxy,dipxz,
     1 dipyx,dipyy,dipyz,dipzx,dipzy,dipzz)
c
c     dipole moment derivatives
c
      implicit REAL  (a-h,o-z)
INCLUDE(../m4/common/sizes)
      logical out,norm
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
INCLUDE(../m4/common/common)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/mapper)
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     3    scrtchpad(150),
     4    ijx(100),ijy(100),ijz(100)
     
c     dimension dd(*)
c
INCLUDE(../m4/common/prnprn)
c     dimension comp(3),dipd(3,3,maxat),dipi(3,3,maxat),
c    &  dipn(3,3,maxat)
      dimension comp(3)
      dimension dipxx(nx),dipyx(nx),dipzx(nx)
      dimension dipxy(nx),dipyy(nx),dipzy(nx)
      dimension dipxz(nx),dipyz(nx),dipzz(nx)
c
      character*5 comp
      data comp/'d/dx','d/dy','d/dz'/
      data ndim/5/
      data one /1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
c     call onepdm(dd,dd(nx+1))
      tol = 2.30258d0*itol
      out = odebug(15)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
c     do 50 n = 1 , nat
c        do 30 i = 1 , 3
c           do 20 j = 1 , 3
c              dipn(i,j,n) = 0.0d0
c              dipi(i,j,n) = 0.0d0
c              dipd(i,j,n) = 0.0d0
c20         continue
c30      continue
c        do 40 i = 1 , 3
c           dipn(i,i,n) = czan(n)
c40      continue
c50   continue
c     ----- ishell
      do 130 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 130
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 120 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,
     +           .false.,ndim,1,1)
            dxx = 0.0d0
            dyy = 0.0d0
            dzz = 0.0d0
            dxy = 0.0d0
            dxz = 0.0d0
            dyx = 0.0d0
            dyz = 0.0d0
            dzx = 0.0d0
            dzy = 0.0d0
c     ----- i primitive
            do 110 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 100 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 70 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 60 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
 60                     continue
 70                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 90 i = mini , maxi
                        in = loci + i
                        do 80 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
c                          dum = dd(nn)*dij(n)
                           dum = dij(n)
c                          dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
                           dipxx(nn) = dipxx(nn) 
     1   + dum*xdipd(nnx)*yin(ny)*zin(nz)
                           dipyy(nn) = dipyy(nn) 
     1   + dum*xin(nnx)*ydipd(ny)*zin(nz)
                           dipzz(nn) = dipzz(nn) 
     1   + dum*xin(nnx)*yin(ny)*zdipd(nz)
                           dipxy(nn) = dipxy(nn) 
     1   + dum*xdip(nnx)*yd(ny)*zin(nz)
                           dipyx(nn) = dipyx(nn) 
     1   + dum*xd(nnx)*ydip(ny)*zin(nz)
                           dipxz(nn) = dipxz(nn) 
     1   + dum*xdip(nnx)*yin(ny)*zd(nz)
                           dipzx(nn) = dipzx(nn) 
     1   + dum*xd(nnx)*yin(ny)*zdip(nz)
                           dipyz(nn) = dipyz(nn) 
     1   + dum*xin(nnx)*ydip(ny)*zd(nz)
                           dipzy(nn) = dipzy(nn) 
     1   + dum*xin(nnx)*yd(ny)*zdip(nz)
c                          dxx = dxx + dum*xdipd(nnx)*yin(ny)*zin(nz)
c                          dyy = dyy + dum*xin(nnx)*ydipd(ny)*zin(nz)
c                          dzz = dzz + dum*xin(nnx)*yin(ny)*zdipd(nz)
c                          dxy = dxy + dum*xdip(nnx)*yd(ny)*zin(nz)
c                          dyx = dyx + dum*xd(nnx)*ydip(ny)*zin(nz)
c                          dxz = dxz + dum*xdip(nnx)*yin(ny)*zd(nz)
c                          dzx = dzx + dum*xd(nnx)*yin(ny)*zdip(nz)
c                          dyz = dyz + dum*xin(nnx)*ydip(ny)*zd(nz)
c                          dzy = dzy + dum*xin(nnx)*yd(ny)*zdip(nz)
 80                     continue
 90                  continue
                  end if
 100           continue
 110        continue
c           dipi(1,1,iat) = dipi(1,1,iat) - dxx
c           dipi(2,2,iat) = dipi(2,2,iat) - dyy
c           dipi(3,3,iat) = dipi(3,3,iat) - dzz
c           dipi(1,2,iat) = dipi(1,2,iat) - dxy
c           dipi(1,3,iat) = dipi(1,3,iat) - dxz
c           dipi(2,1,iat) = dipi(2,1,iat) - dyx
c           dipi(3,1,iat) = dipi(3,1,iat) - dzx
c           dipi(2,3,iat) = dipi(2,3,iat) - dyz
c           dipi(3,2,iat) = dipi(3,2,iat) - dzy
 120     continue
 130  continue
c
c     if (out) then
c        write (iwr,6010)
c        do 150 n = 1 , nat
c           write (iwr,6020)
c           do 140 nc = 1 , 3
c              write (iwr,6030) zaname(n) , comp(nc) ,
c    +                         (dipi(nn,nc,n),nn=1,3)
c140        continue
c150     continue
c     end if
c     do 180 i = 1 , 3
c        do 170 j = 1 , 3
c           do 160 k = 1 , nat
c              dipd(i,j,k) = dipd(i,j,k) + dipn(i,j,k) + dipi(i,j,k)
c160        continue
c170     continue
c180  continue
      return
 6010 format (//35x,'integral derivative contribution'//30x,'x',15x,'y',
     +        15x,'z',/)
 6020 format (//)
 6030 format (5x,a8,5x,a5,3f16.8)
      end
      subroutine vintdrf
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
      common/intdip/pint,qint,rint,pintd,qintd,rintd,
     1t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,cx,cy,cz
INCLUDE(../m4/common/hermit)
INCLUDE(../m4/common/wermit)
      dimension min(6),max(6)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data dzero /0.0d0/
c
      pint = dzero
      qint = dzero
      rint = dzero
      npts = (ni+nj-2+1)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 140 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)*t
         ptx = dum + p0
         pty = dum + q0
         ptz = dum + r0
         ax = ptx - pi
         ay = pty - qi
         az = ptz - ri
         bx = ptx - pj
         by = pty - qj
         bz = ptz - rj
         go to (60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      go to (130,120,110,100,90,80,70) , nj
 70      px = px*bx
         py = py*by
         pz = pz*bz
 80      px = px*bx
         py = py*by
         pz = pz*bz
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     pint = pint + px
         qint = qint + py
         rint = rint + pz
 140  continue
      return
      end
      subroutine idxadrf(ijx,ijy,ijz,ij,mini,maxi,
     &minj,maxj,iandj,inc1,inc2,inc3)
      implicit REAL  (a-h,o-z)
      logical iandj
      dimension ijx(100),ijy(100),ijz(100)
INCLUDE(../m4/common/inxblk)
      dimension jx(20),jy(20),jz(20)
c
_IF1(ct)cdir$ novector
_IF1(a)cvd$  novector
_IF1(x)c$dir scalar
      do 20 j = minj , maxj
         jx(j) = ix(j)*inc2
         jy(j) = iy(j)*inc2
         jz(j) = iz(j)*inc2
 20   continue
      ij = 0
      jmax = maxj
      do 40 i = mini , maxi
         nx = ix(i)*inc1 + inc3
         ny = iy(i)*inc1 + inc3
         nz = iz(i)*inc1 + inc3
c        if (iandj) jmax = i
         do 30 j = minj , jmax
            ij = ij + 1
            ijx(ij) = nx + jx(j)
            ijy(ij) = ny + jy(j)
            ijz(ij) = nz + jz(j)
 30      continue
 40   continue
_IF1(ct)cdir$ vector
      return
      end
      subroutine qmdrfder(iatom,quxxx,quxxy,quxxz,
     1 quyyx,quyyy,quyyz,quzzx,quzzy,quzzz,
     2 quxyx,quxyy,quxyz,quxzx,quxzy,quxzz,
     3 quyzx,quyzy,quyzz)
c----------------------------------------------------------------
c     quadrupole derivative integrals
c----------------------------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(../m4/common/sizes)
      logical out,norm
INCLUDE(../m4/common/prnprn)
INCLUDE(../m4/common/symtry)
INCLUDE(../m4/common/timez)
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
INCLUDE(../m4/common/common)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/cndx41)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
c
      common/small/dipd(3,3,maxat),qudd(6,3,maxat)
INCLUDE(../m4/common/mapper)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     4    qxx(25),qyy(25),qzz(25),qxxd(25),qyyd(25),qzzd(25),
     5    ijx(100),ijy(100),ijz(100)
c     dimension dd(*)
      dimension quxxx(nx),quxxy(nx),quxxz(nx)
      dimension quyyx(nx),quyyy(nx),quyyz(nx)
      dimension quzzx(nx),quzzy(nx),quzzz(nx)
      dimension quxyx(nx),quxyy(nx),quxyz(nx)
      dimension quxzx(nx),quxzy(nx),quxzz(nx)
      dimension quyzx(nx),quyzy(nx),quyzz(nx)
       character *8 comp,grhf,closed
      dimension comp(3)
      dimension derivs(27*maxat)
      equivalence (derivs(1),dipd(1,1,1))
      data comp/'d/dx','d/dy','d/dz'/
      data grhf/'grhf'/
      data closed/'closed'/
      data ndim/5/
      data two,thalf/2.0d0,1.5d0/
      data zer,one /0.0d0,1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
c     call secget(isect(7),7,iblok)
c     call rdedx(dd,nx,iblok,ifild)
c     if (scftyp.ne.grhf .and. scftyp.ne.closed) then
c        call secget(isect(10),10,iblok)
c        call rdedx(dd(nx+1),nx,iblok,ifild)
c        do 20 i = 1 , nx
c           dd(i) = dd(i+nx) + dd(i)
c20      continue
c     end if
      tol = 2.30258d0*itol
      out = odebug(19)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
c     do 30 n = 1 , nat
c        zz = czan(n)
c        zx = (c(1,n)-gx)*zz
c        zy = (c(2,n)-gy)*zz
c        zz = (c(3,n)-gz)*zz
c        qudd(4,1,n) = thalf*zy
c        qudd(4,2,n) = thalf*zx
c        qudd(4,3,n) = zer
c        qudd(5,1,n) = thalf*zz
c        qudd(5,2,n) = zer
c        qudd(5,3,n) = thalf*zx
c        qudd(6,1,n) = zer
c        qudd(6,2,n) = zz*thalf
c        qudd(6,3,n) = zy*thalf
c        qudd(1,1,n) = two*zx
c        qudd(1,2,n) = -zy
c        qudd(1,3,n) = -zz
c        qudd(2,1,n) = -zx
c        qudd(2,2,n) = two*zy
c        qudd(2,3,n) = -zz
c        qudd(3,1,n) = -zx
c        qudd(3,2,n) = -zy
c        qudd(3,3,n) = two*zz
c30   continue
c     ----- ishell
      do 110 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 110
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,
     +                  .false.,ndim,1,1)
            qxxx = 0.0d0
            qxxy = 0.0d0
            qxxz = 0.0d0
            qyyx = 0.0d0
            qyyy = 0.0d0
            qyyz = 0.0d0
            qzzx = 0.0d0
            qzzy = 0.0d0
            qzzz = 0.0d0
            qxyx = 0.0d0
            qxyy = 0.0d0
            qxyz = 0.0d0
            qxzx = 0.0d0
            qxzy = 0.0d0
            qxzz = 0.0d0
            qyzx = 0.0d0
            qyzy = 0.0d0
            qyzz = 0.0d0
c     ----- i primitive
            do 90 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 80 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 50 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 40 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
                           call qmsint()
                           qxx(jn) = xintx*tinv
                           qyy(jn) = yinty*tinv
                           qzz(jn) = zintz*tinv
 40                     continue
 50                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
                     call oneld(qxx,qyy,qzz,qxxd,qyyd,qzzd,ai,lit,ljt,1,
     +                          ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 70 i = mini , maxi
                        in = loci + i
                        do 60 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
c                          dum = dd(nn)*dij(n)
                           dum = dij(n)
c                          dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
c
                           quxxx(nn) = quxxx(nn)
     1                      + dum*qxxd(nnx)*yin(ny)*zin(nz)
c
c                          qxxx = qxxx + dum*qxxd(nnx)*yin(ny)*zin(nz)
                           quxxy(nn) = quxxy(nn)
     1                      + dum*qxx(nnx)*yd(ny)*zin(nz)
c                          qxxy = qxxy + dum*qxx(nnx)*yd(ny)*zin(nz)
                           quxxz(nn) = quxxz(nn)
     1                      + dum*qxx(nnx)*yin(ny)*zd(nz)
c                          qxxz = qxxz + dum*qxx(nnx)*yin(ny)*zd(nz)
c
c
                           quyyx(nn) = quyyx(nn)
     1                      + dum*xd(nnx)*qyy(ny)*zin(nz)
c                          qyyx = qyyx + dum*xd(nnx)*qyy(ny)*zin(nz)
                           quyyy(nn) = quyyy(nn)
     1                      + dum*xin(nnx)*qyyd(ny)*zin(nz)
c                          qyyy = qyyy + dum*xin(nnx)*qyyd(ny)*zin(nz)
                           quyyz(nn) = quyyz(nn)
     1                      + dum*xin(nnx)*qyy(ny)*zd(nz)
c                          qyyz = qyyz + dum*xin(nnx)*qyy(ny)*zd(nz)
c
c
                           quzzx(nn) = quzzx(nn)
     1                      + dum*xd(nnx)*yin(ny)*qzz(nz)
c                          qzzx = qzzx + dum*xd(nnx)*yin(ny)*qzz(nz)
                           quzzy(nn) = quzzy(nn)
     1                      + dum*xin(nnx)*yd(ny)*qzz(nz)
c                          qzzy = qzzy + dum*xin(nnx)*yd(ny)*qzz(nz)
                           quzzz(nn) = quzzz(nn)
     1                      + dum*xin(nnx)*yin(ny)*qzzd(nz)
c                          qzzz = qzzz + dum*xin(nnx)*yin(ny)*qzzd(nz)
c
c
                           quxyx(nn) = quxyx(nn)
     1                      + dum*xdipd(nnx)*ydip(ny)*zin(nz)
c                          qxyx = qxyx + dum*xdipd(nnx)*ydip(ny)*zin(nz)
                           quxyy(nn) = quxyy(nn)
     1                      + dum*xdip(nnx)*ydipd(ny)*zin(nz)
c                          qxyy = qxyy + dum*xdip(nnx)*ydipd(ny)*zin(nz)
                           quxyz(nn) = quxyz(nn)
     1                      + dum*xdip(nnx)*ydip(ny)*zd(nz)
c                          qxyz = qxyz + dum*xdip(nnx)*ydip(ny)*zd(nz)
c
c
                           quxzx(nn) = quxzx(nn)
     1                      + dum*xdipd(nnx)*yin(ny)*zdip(nz)
c                          qxzx = qxzx + dum*xdipd(nnx)*yin(ny)*zdip(nz)
                           quxzy(nn) = quxzy(nn)
     1                      + dum*xdip(nnx)*yd(ny)*zdip(nz)
c                          qxzy = qxzy + dum*xdip(nnx)*yd(ny)*zdip(nz)
                           quxzz(nn) = quxzz(nn)
     1                      + dum*xdip(nnx)*yin(ny)*zdipd(nz)
c                          qxzz = qxzz + dum*xdip(nnx)*yin(ny)*zdipd(nz)
c
c
                           quyzx(nn) = quyzx(nn)
     1                      + dum*xd(nnx)*ydip(ny)*zdip(nz)
c                          qyzx = qyzx + dum*xd(nnx)*ydip(ny)*zdip(nz)
                           quyzy(nn) = quyzy(nn)
     1                      + dum*xin(nnx)*ydipd(ny)*zdip(nz)
c                          qyzy = qyzy + dum*xin(nnx)*ydipd(ny)*zdip(nz)
                           quyzz(nn) = quyzz(nn)
     1                      + dum*xin(nnx)*ydip(ny)*zdipd(nz)
c                          qyzz = qyzz + dum*xin(nnx)*ydip(ny)*zdipd(nz)
 60                     continue
 70                  continue
                  end if
 80            continue
 90         continue
c           qudd(1,1,iat) = qudd(1,1,iat) - 0.5d0*(qxxx+qxxx-qyyx-qzzx)
c           qudd(2,2,iat) = qudd(2,2,iat) - 0.5d0*(qyyy+qyyy-qxxy-qzzy)
c           qudd(3,3,iat) = qudd(3,3,iat) - 0.5d0*(qzzz+qzzz-qxxz-qyyz)
c           qudd(1,2,iat) = qudd(1,2,iat) - 0.5d0*(qxxy+qxxy-qzzy-qyyy)
c           qudd(1,3,iat) = qudd(1,3,iat) - 0.5d0*(qxxz+qxxz-qyyz-qzzz)
c           qudd(2,1,iat) = qudd(2,1,iat) - 0.5d0*(qyyx+qyyx-qxxx-qzzx)
c           qudd(3,1,iat) = qudd(3,1,iat) - 0.5d0*(qzzx+qzzx-qxxx-qyyx)
c           qudd(2,3,iat) = qudd(2,3,iat) - 0.5d0*(qyyz+qyyz-qxxz-qzzz)
c           qudd(3,2,iat) = qudd(3,2,iat) - 0.5d0*(qzzy+qzzy-qxxy-qyyy)
c           qudd(4,1,iat) = qudd(4,1,iat) - 1.5d0*qxyx
c           qudd(4,2,iat) = qudd(4,2,iat) - 1.5d0*qxyy
c           qudd(4,3,iat) = qudd(4,3,iat) - 1.5d0*qxyz
c           qudd(5,1,iat) = qudd(5,1,iat) - 1.5d0*qxzx
c           qudd(5,2,iat) = qudd(5,2,iat) - 1.5d0*qxzy
c           qudd(5,3,iat) = qudd(5,3,iat) - 1.5d0*qxzz
c           qudd(6,1,iat) = qudd(6,1,iat) - 1.5d0*qyzx
c           qudd(6,2,iat) = qudd(6,2,iat) - 1.5d0*qyzy
c           qudd(6,3,iat) = qudd(6,3,iat) - 1.5d0*qyzz
 100     continue
 110  continue
c
c     if (out) then
c        write (iwr,6010)
c        do 130 n = 1 , nat
c           write (iwr,6020)
c           do 120 nc = 1 , 3
c              write (iwr,6030) zaname(n) , comp(nc) ,
c    +                         (qudd(nn,nc,n),nn=1,6)
c120        continue
c130     continue
c     end if
c     term from differentiating orbitals
c     nplus1 = nocca + 1
c     iposs = iochf(15)
c     icomp = 24
c     ioff = 1
c     i1 = nx + 1
c     i2 = i1 + nx
c     i3 = i2 + nx
c     lennew = iky(ncoorb+1)
c     iblll = lensec(lennew)
c     do 140 i = 1 , 3
c        ioff = ioff + nx
c        icomp = icomp + 1
c        call vclr(dd(ioff),1,lennew)
c        if (itwo(i).eq.0) then
c           write (iwr,6040) i
c        else
c           call secget(isect(icomp),icomp,iblok)
c           call rdedx(dd(ioff),lennew,iblok,ifild)
c        end if
c140  continue
c     do 160 n = 1 , nat
c        do 150 nc = 1 , 3
cIF(secd_parallel)
c          call fetch(dd,lennew,'pdens',(n-1)*3+nc)
cELSE
c           call rdedx(dd,lennew,iposs,ifockf)
cENDIF
c           iposs = iposs + iblll
c           dumx = tracep(dd(i1),dd,ncoorb)
c           dumy = tracep(dd(i2),dd,ncoorb)
c           dumz = tracep(dd(i3),dd,ncoorb)
c           qudd(1,nc,n) = qudd(1,nc,n) - dumx
c           qudd(2,nc,n) = qudd(2,nc,n) - dumy
c           qudd(3,nc,n) = qudd(3,nc,n) - dumz
c150     continue
c160  continue
c     icomp = 27
c     ioff = 1
c     do 170 i = 4 , 6
c        icomp = icomp + 1
c        ioff = ioff + nx
c        call vclr(dd(ioff),1,lennew)
c        if (itwo(i).eq.0) then
c           write (iwr,6040) i
c        else
c           call secget(isect(icomp),icomp,iblok)
c           call rdedx(dd(ioff),lennew,iblok,ifild)
c        end if
c170  continue
c     iposs = iochf(15)
c     do 190 n = 1 , nat
c        do 180 nc = 1 , 3
cIF(secd_parallel)
csecd
c           call fetch(dd,lennew,'pdens',(n-1)*3+nc)
cELSE
c           call rdedx(dd,lennew,iposs,ifockf)
cENDIF
c           iposs = iposs + iblll
c           qudd(4,nc,n) = qudd(4,nc,n) - tracep(dd,dd(i1),ncoorb)
c           qudd(5,nc,n) = qudd(5,nc,n) - tracep(dd,dd(i2),ncoorb)
c           qudd(6,nc,n) = qudd(6,nc,n) - tracep(dd,dd(i3),ncoorb)
c180     continue
c190  continue
c     call qmdsym(qudd,dd,dd(nw196(5)+1),nat,nshell)
c     if (oprn(40)) then
c      if (oprn(25)) write (iwr,6010)
c      do 210 n = 1 , nat
c         if (oprn(25)) write (iwr,6020)
c         do 200 nc = 1 , 3
c            if (oprn(25)) write (iwr,6030) zaname(n) , comp(nc) , 
c    +      (qudd(nn,nc,n),nn=1,6)
c200      continue
c210   continue
c     endif
      call timit(3)
c     lenblk = lensec(27*maxat)
c     call secput(isect(50),50,lenblk,iblok)
c     call revind
c     lds(isect(50)) = 27*maxat
c     call wrt3(derivs,lds(isect(50)),iblok,ifild)
c     call revise
c     call clredx
      return
 6010 format (/
     +   20x,'*************************************'/
     +   20x,'*   quadrupole moment derivatives   *'/
     +   20x,'*           (atomic units)          *'/
     +   20x,'*************************************'//
     + 29x,'xx',14x,'yy',14x,'zz',14x,'xy',14x,'xz',14x,'yz'/)
 6020 format (/)
 6030 format (5x,a8,2x,a8,6f16.8)
 6040 format (//10x,'quadrupole integrals missing for component',i6)
      end
