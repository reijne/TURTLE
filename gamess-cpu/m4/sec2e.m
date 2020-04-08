c 
c  $Author: mrdj $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/sec2e.m,v $
c  $State: Exp $
c  
      subroutine dr2den(abdens,da,db,dc,ia,ii,jj,kk,ll,q4,scftyp,oss,
     &  mp2)
c
c     density matrices for scf wavefunctions
c
      implicit REAL  (a-h,o-z)
      logical mp2,oss
      character *8 scftyp
      dimension abdens(*),da(*),db(*),dc(*),ia(*)
c
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/ghfblk)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/incrdd)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
      logical open, ijeq, kleq
      character *8 closed,grhf
      data closed/'closed'/
      data grhf/'grhf'/
      data half,four /0.5d0,4.0d0/
c
      open = scftyp.ne.closed
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
_IF(ccpdft)
      hf_wght = CD_HF_exchange_weight()
_ENDIF
      if (cigr) then
c
c    additions for ci core
c
         if (ncore.eq.0) return
         ni = 1
         do 50 i = mini , maxi
            nj = ni
            i1 = loci + i
            ii1 = ia(i1)
            do 40 j = minj , maxj
               nk = nj
               j1 = locj + j
               jj1 = ia(j1)
               mij = ii1 + j1
               if (j1.gt.i1) mij = jj1 + i1
               ijeq = i1.eq.j1
               do 30 k = mink , maxk
                  nl = nk
                  k1 = lock + k
                  kk1 = ia(k1)
                  mik = ii1 + k1
                  if (k1.gt.i1) mik = kk1 + i1
                  mjk = jj1 + k1
                  if (k1.gt.j1) mjk = kk1 + j1
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                  do 20 l = minl , maxl
                     nn = nl
                     l1 = locl + l
                     ll1 = ia(l1)
                     mkl = kk1 + l1
                     if (l1.gt.k1) mkl = ll1 + k1
                     mjl = jj1 + l1
                     if (l1.gt.j1) mjl = ll1 + j1
                     mil = ii1 + l1
                     if (l1.gt.i1) mil = ll1 + i1
                     kleq = k1.eq.l1
c
c     db ----- core density
c     dc ----- ci density
c
                     dfaca = db(mij)*db(mkl)*four - db(mik)*db(mjl)
     +                       - db(mjk)*db(mil)
                     dfacb = (db(mij)*dc(mkl)+db(mkl)*dc(mij))
     +                       *four - db(mil)*dc(mjk) - db(mjk)*dc(mil)
     +                       - db(mik)*dc(mjl) - db(mjl)*dc(mik)
                     dfac = dfaca*four + dfacb + dfacb
                     if (ijeq) dfac = dfac*half
                     if (kleq) dfac = dfac*half
                     if (i1.eq.k1 .and. j1.eq.l1) dfac = dfac*half
                     abdens(nn) = abdens(nn) + dfac*q4
                     nl = nl + inc2
 20               continue
                  nk = nk + inc3
 30            continue
               nj = nj + inc4
 40         continue
            ni = ni + inc5
 50      continue
         return
      else if (scftyp.eq.grhf) then
c
c     general case
c
         ni = 1
         do 110 i = mini , maxi
            nj = ni
            i1 = loci + i
            ii1 = ia(i1)
            do 100 j = minj , maxj
               nk = nj
               j1 = locj + j
               jj1 = ia(j1)
               mij = ii1 + j1
               if (j1.gt.i1) mij = jj1 + i1
               ijeq = i1.eq.j1
               do 90 k = mink , maxk
                  nl = nk
                  k1 = lock + k
                  kk1 = ia(k1)
                  mik = ii1 + k1
                  if (k1.gt.i1) mik = kk1 + i1
                  mjk = jj1 + k1
                  if (k1.gt.j1) mjk = kk1 + j1
                  do 80 l = minl , maxl
                     nn = nl
                     l1 = locl + l
                     ll1 = ia(l1)
                     mkl = kk1 + l1
                     if (l1.gt.k1) mkl = ll1 + k1
                     mjl = jj1 + l1
                     if (l1.gt.j1) mjl = ll1 + j1
                     mil = ii1 + l1
                     if (l1.gt.i1) mil = ll1 + i1
                     kleq = k1.eq.l1
                     dfac = 0.0d0
                     ioff = 0
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                     do 70 is = 1 , njk
                        isi = (is-1)*11
                        joff = 0
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                        do 60 js = 1 , njk
                           dfac = dfac + 4.0d0*erga(isi+js)
     +                            *(da(ioff+mij)*da(joff+mkl)
     +                            +da(ioff+mkl)*da(joff+mij))
     +                            + 2.0d0*ergb(isi+js)
     +                            *(da(ioff+mik)*da(joff+mjl)
     +                            +da(ioff+mjl)*da(joff+mik)
     +                            +da(ioff+mil)*da(joff+mjk)
     +                            +da(ioff+mjk)*da(joff+mil))
                           joff = joff + nx
 60                     continue
                        ioff = ioff + nx
 70                  continue
                     if (ijeq) dfac = dfac*half
                     if (kleq) dfac = dfac*half
                     if (i1.eq.k1 .and. j1.eq.l1) dfac = dfac*half
                     abdens(nn) = dfac*q4
                     nl = nl + inc2
 80               continue
                  nk = nk + inc3
 90            continue
               nj = nj + inc4
 100        continue
            ni = ni + inc5
 110     continue
         return
      else
         ni = 1
         do 150 i = mini , maxi
            nj = ni
            i1 = loci + i
            ii1 = ia(i1)
            do 140 j = minj , maxj
               nk = nj
               j1 = locj + j
               jj1 = ia(j1)
               mij = ii1 + j1
               if (j1.gt.i1) mij = jj1 + i1
               dij = da(mij)*four
               ijeq = i1.eq.j1
               do 130 k = mink , maxk
                  nl = nk
                  k1 = lock + k
                  kk1 = ia(k1)
                  mik = ii1 + k1
                  if (k1.gt.i1) mik = kk1 + i1
                  mjk = jj1 + k1
                  if (k1.gt.j1) mjk = kk1 + j1
_IF(ccpdft)
                  dik = da(mik)*hf_wght
                  djk = da(mjk)*hf_wght
_ELSE
                  dik = da(mik)
                  djk = da(mjk)
_ENDIF
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                  do 120 l = minl , maxl
                     nn = nl
                     l1 = locl + l
                     ll1 = ia(l1)
                     mkl = kk1 + l1
                     if (l1.gt.k1) mkl = ll1 + k1
                     mjl = jj1 + l1
                     if (l1.gt.j1) mjl = ll1 + j1
                     mil = ii1 + l1
                     if (l1.gt.i1) mil = ll1 + i1
                     kleq = k1.eq.l1
                     dfac = dij*da(mkl) - dik*da(mjl) - djk*da(mil)
_IF(ccpdft)
                     if (open) dfac = dfac - (db(mik)*db(mjl) + db(mjk)
     +                                *db(mil))*hf_wght
                     if (oss) dfac = dfac - (dc(mik)*dc(mjl) + dc(mjk)
     +                               *dc(mil))*hf_wght
     +                               + 3.0d0*(db(mik)*dc(mjl)+dc(mik)
     +                               *db(mjl)+db(mjk)*dc(mil)+dc(mjk)
     +                               *db(mil))
_ELSE
                     if (open) dfac = dfac - db(mik)*db(mjl) - db(mjk)
     +                                *db(mil)
                     if (oss) dfac = dfac - dc(mik)*dc(mjl) - dc(mjk)
     +                               *dc(mil)
     +                               + 3.0d0*(db(mik)*dc(mjl)+dc(mik)
     +                               *db(mjl)+db(mjk)*dc(mil)+dc(mjk)
     +                               *db(mil))
_ENDIF
                     if (ijeq) dfac = dfac*half
                     if (kleq) dfac = dfac*half
                     if (i1.eq.k1 .and. j1.eq.l1) dfac = dfac*half
                     if (mp2) then
                        abdens(nn) = abdens(nn) + dfac*q4
                     else
                        abdens(nn) = dfac*q4
                     end if
                     nl = nl + inc2
 120              continue
                  nk = nk + inc3
 130           continue
               nj = nj + inc4
 140        continue
            ni = ni + inc5
 150     continue
         return
      end if
      end

      subroutine dr2fre(core,fc,a,vec,e,rm,coord,g,nat3,iw,last)
c--------------------------------------------------------------------
c     analyse second derivative matrix to give frequencies
c     modification of hondo routine
c--------------------------------------------------------------------
c
      implicit REAL  (a-h,o-z)
      dimension fc(*),a(*),vec(nat3,*),e(*),rm(*),coord(*)
      dimension g(*)
      dimension core(*)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/prnprn)
INCLUDE(common/mapper)
INCLUDE(common/phycon)
INCLUDE(common/runopt)
INCLUDE(common/maxlen)
c
      dimension del(3*maxat),cnew(3*maxat)
      character*1 clab
      dimension clab(3)
      data clab/'x','y','z'/
      data eps/1.0d-10/
      data done/1.0d0/
      data ten5,ten18,four/1.0d5,1.0d18,4.0d0/
c
c     Compute the conversion factor from physical constants to 
c     ensure consistency with other quantities in the code.
c
c!old data tfact /2.6436411815482d+07/
c
c     avog:    avogadro constant, in mol**-1
c     slight:  speed of light, in cm/sec
c     hartre:  joules per hartree
c     bohr:   bohrs per angstrom
c
      avog = toang(5)
      slight = toang(9)
      hartre = toang(8)
      bohr  = toang(1)
c
c     factor:  converts frequencies to wavenumbers.
c     conver:  converts force constants in atomic units to
c              mdyn/angstrom**2
c
      pi = dacos(-done) 
      factor = (ten5/(four*pi*pi))*(avog/slight/slight)
      conver = (ten18*hartre)/(bohr*bohr)
      tfact = factor*conver
c
c
      nnc = (nat3*(nat3+1))/2
*
* Allocate memory for call to thermochemsitry
*
      i10 = 1
      i20 = i10 + nnc
      i30 = i20 + nat3*nat3
      i40 = i30 + nat3*nat3
      iend = i40 + nat3
      if ( last + iend .gt. maxq ) then
            write (iw,6020) maxq , iend+last
            call caserr('dr2fre: insufficient core available')
      endif
*
      do 30 i = 1 , nat3
         do 20 j = 1 , i
            fc(iky(i)+j) = vec(i,j)
 20      continue
 30   continue
c
c     ----- diagonalize cartesian force matrix -----
c
      call dcopy(nnc,fc,1,a,1)
      call ngdiag(a,vec,e,iky,nat3,nat3,1,1.0d-20)
c
c     work out quadratic corrections to geometry
c
      do 50 i = 1 , nat3
         del(i) = 0.0d0
         cnew(i) = 0.0d0
 50   continue
      do 80 j = 1 , nat3
         if (dabs(g(j)).gt.eps) then
            do 70 k = 1 , nat3
               s = -g(j)*vec(j,k)
c
c    ignore any modes with freq le 10 wave no.
c    note : if not at stationary point this may predict
c    rotation
c         : if very flat surface some modes may be
c    omitted
c
               if (dsqrt(dabs(e(k)*tfact)).gt.10.0d0) then
                  s = s/e(k)
                  do 60 i = 1 , nat3
                     del(i) = del(i) + s*vec(i,k)
 60               continue
               end if
 70         continue
         end if
 80   continue
      do 90 i = 1 , nat3
         cnew(i) = coord(i) + del(i)
 90   continue
      if (oprn(24)) then
       write (iw,6070)
       n = 0
       do 110 i = 1 , nat
          do 100 j = 1 , 3
             n = n + 1
             write (iw,6060) zaname(i) , clab(j) , coord(n) , 
     +                       g(n) ,  del(n) ,  cnew(n)
 100      continue
        write (iw,6010)
 110   continue
      endif
c
c-weights
c loop over vectors of atomic masses
c
      nmv  = mass_numvec()
      do imv = 1, nmv

         call mass_prvec(imv)

c     
c     ----- create full force matrix with mass weighted transformation -
c
         call dr2mas(imv,rm)

         do 130 i = 1 , nat3
            do 120 j = 1 , i
               ij = iky(i) + j
               a(ij) = rm(i)*fc(ij)*rm(j)
 120        continue
 130     continue

c
c     ----- get normal modes and frequencies -----
c
         call ngdiag(a,vec,e,iky,nat3,nat3,1,1.0d-20)
c
c     ----- convert frequencies to inverse cm -----
c
         do 140 i = 1 , nat3
            e(i) = dsign(1.0d0,e(i))*dsqrt(dabs(tfact*e(i)))
 140     continue
c
c     ----- construct mass weighted displacement vectors
c           and renormalize them -----
c
         do 160 i = 1 , nat3
            do 150 j = 1 , nat3
               vec(j,i) = vec(j,i)*rm(j)
 150        continue
 160     continue
         if (oprn(5)) then
            write (iw,6080)
            mmax = 0
 170        mmin = mmax + 1
            mmax = mmax + 8
            if (mmax.gt.nat3) mmax = nat3
            write (iw,6010)
            write (iw,6020) (j,j=mmin,mmax)
            write (iw,6010)
            write (iw,6030) (e(j),j=mmin,mmax)
            write (iw,6010)
            do 180 iat = 1 , nat
               i0 = 3*(iat-1)
               write (iw,6040) iat , zaname(iat) , clab(1) ,
     +              (vec(i0+1,j),j=mmin,mmax)
               write (iw,6050) clab(2) , (vec(i0+2,j),j=mmin,mmax)
               write (iw,6050) clab(3) , (vec(i0+3,j),j=mmin,mmax)
 180        continue
            if (mmax.lt.nat3) go to 170
         else
            write (iw,6090) (e(j),j=1,nat3)
         end if
*
* Compute thermo chemistry for this force constant matrix
*
          call dcopy(nnc,fc,1,core(i10),1)
          call vibfrq(nat,mul,iky,coord,nat3,core(i10),rm,
     +            core(i20),core(i30),
     +            vec,core(i40),a,e,toang,iw)
*
      enddo
      return
 6010 format (/)
 6020 format (18x,8(3x,i3,4x))
 6030 format (18x,8f10.6)
 6040 format (i3,3x,a8,3x,a1,8f10.6)
 6050 format (17x,a1,8f10.6)
 6060 format (10x,a8,3x,a1,4f20.8)
 6070 format (//
     + 10x,'========================='/
     + 10x,'predicted geometry change'/
     + 10x,'========================='//
     + 10x,'atom',12x,'present position',12x,'gradient',8x,
     + 'displacement',8x,'new position'//)
 6080 format (/
     +20x,'**********************************************'/
     +20x,'*    harmonic frequencies and normal modes   *'/
     +20x,'* translation and rotation not projected out *'/
     +20x,'**********************************************'/)
 6090 format (/
     + 20x,'************************************'/
     + 20x,'* unprojected harmonic frequencies *'/
     + 20x,'************************************'//
     + (1x,8f12.3))
      end
      subroutine dr2as1(x,y,z,xa,ya,za,xdd,ydd,zdd,ix,iy,iz,
     + abdens,ncdim)
c----------------------------
c     assembly routine case 1
c----------------------------
      implicit REAL  (a-h,o-z)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*)
      dimension xa(ncdim,*),ya(ncdim,*),za(ncdim,*)
      dimension xdd(ncdim,*),ydd(ncdim,*),zdd(ncdim,*)
      dimension abdens(*),ix(*),iy(*),iz(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/root)
INCLUDE(common/ddshln)
INCLUDE(common/incrdd)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225),dd(144)
c
      logical unroll
c
cvd$r assoc
c
c     i1 = inc1
      g1 = 0.0d0
      g2 = 0.0d0
      g3 = 0.0d0
      g5 = 0.0d0
      g6 = 0.0d0
      g9 = 0.0d0
      unroll = ncontr.le.6 .and. ijkld2.ge.16
      if (.not.unroll) then
         do 30 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
c
c
            ss1 = 0.0d0
            ss2 = 0.0d0
            ss3 = 0.0d0
            ss5 = 0.0d0
            ss6 = 0.0d0
            ss9 = 0.0d0
            do 20 nr = 1 , ncontr
               ss1 = ss1 + xdd(nr,mx)*y(nr,my)*z(nr,mz)
               ss2 = ss2 + xa(nr,mx)*ya(nr,my)*z(nr,mz)
               ss3 = ss3 + xa(nr,mx)*y(nr,my)*za(nr,mz)
               ss5 = ss5 + x(nr,mx)*ydd(nr,my)*z(nr,mz)
               ss6 = ss6 + x(nr,mx)*ya(nr,my)*za(nr,mz)
               ss9 = ss9 + x(nr,mx)*y(nr,my)*zdd(nr,mz)
 20         continue
            g1 = ab*ss1 + g1
            g2 = ab*ss2 + g2
            g3 = ab*ss3 + g3
            g5 = ab*ss5 + g5
            g6 = ab*ss6 + g6
            g9 = ab*ss9 + g9
 30      continue
      else
         go to (40,60,80,100,120,140) , ncontr
 40      do 50 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*xdd(1,mx)*y(1,my)*z(1,mz)
            g2 = g2 + ab*xa(1,mx)*ya(1,my)*z(1,mz)
            g3 = g3 + ab*xa(1,mx)*y(1,my)*za(1,mz)
            g5 = g5 + ab*x(1,mx)*ydd(1,my)*z(1,mz)
            g6 = g6 + ab*x(1,mx)*ya(1,my)*za(1,mz)
            g9 = g9 + ab*x(1,mx)*y(1,my)*zdd(1,mz)
 50      continue
         go to 160
 60      do 70 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*(xdd(1,mx)*y(1,my)*z(1,mz)+xdd(2,mx)*y(2,my)
     +           *z(2,mz))
            g2 = g2 + ab*(xa(1,mx)*ya(1,my)*z(1,mz)+xa(2,mx)*ya(2,my)
     +           *z(2,mz))
            g3 = g3 + ab*(xa(1,mx)*y(1,my)*za(1,mz)+xa(2,mx)*y(2,my)
     +           *za(2,mz))
            g5 = g5 + ab*(x(1,mx)*ydd(1,my)*z(1,mz)+x(2,mx)*ydd(2,my)
     +           *z(2,mz))
            g6 = g6 + ab*(x(1,mx)*ya(1,my)*za(1,mz)+x(2,mx)*ya(2,my)
     +           *za(2,mz))
            g9 = g9 + ab*(x(1,mx)*y(1,my)*zdd(1,mz)+x(2,mx)*y(2,my)
     +           *zdd(2,mz))
 70      continue
         go to 160
 80      do 90 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*(xdd(1,mx)*y(1,my)*z(1,mz)+xdd(2,mx)*y(2,my)
     +           *z(2,mz)+xdd(3,mx)*y(3,my)*z(3,mz))
            g2 = g2 + ab*(xa(1,mx)*ya(1,my)*z(1,mz)+xa(2,mx)*ya(2,my)
     +           *z(2,mz)+xa(3,mx)*ya(3,my)*z(3,mz))
            g3 = g3 + ab*(xa(1,mx)*y(1,my)*za(1,mz)+xa(2,mx)*y(2,my)
     +           *za(2,mz)+xa(3,mx)*y(3,my)*za(3,mz))
            g5 = g5 + ab*(x(1,mx)*ydd(1,my)*z(1,mz)+x(2,mx)*ydd(2,my)
     +           *z(2,mz)+x(3,mx)*ydd(3,my)*z(3,mz))
            g6 = g6 + ab*(x(1,mx)*ya(1,my)*za(1,mz)+x(2,mx)*ya(2,my)
     +           *za(2,mz)+x(3,mx)*ya(3,my)*za(3,mz))
            g9 = g9 + ab*(x(1,mx)*y(1,my)*zdd(1,mz)+x(2,mx)*y(2,my)
     +           *zdd(2,mz)+x(3,mx)*y(3,my)*zdd(3,mz))
 90      continue
         go to 160
 100     do 110 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*(xdd(1,mx)*y(1,my)*z(1,mz)+xdd(2,mx)*y(2,my)
     +           *z(2,mz)+xdd(3,mx)*y(3,my)*z(3,mz)+xdd(4,mx)*y(4,my)
     +           *z(4,mz))
            g2 = g2 + ab*(xa(1,mx)*ya(1,my)*z(1,mz)+xa(2,mx)*ya(2,my)
     +           *z(2,mz)+xa(3,mx)*ya(3,my)*z(3,mz)+xa(4,mx)*ya(4,my)
     +           *z(4,mz))
            g3 = g3 + ab*(xa(1,mx)*y(1,my)*za(1,mz)+xa(2,mx)*y(2,my)
     +           *za(2,mz)+xa(3,mx)*y(3,my)*za(3,mz)+xa(4,mx)*y(4,my)
     +           *za(4,mz))
            g5 = g5 + ab*(x(1,mx)*ydd(1,my)*z(1,mz)+x(2,mx)*ydd(2,my)
     +           *z(2,mz)+x(3,mx)*ydd(3,my)*z(3,mz)+x(4,mx)*ydd(4,my)
     +           *z(4,mz))
            g6 = g6 + ab*(x(1,mx)*ya(1,my)*za(1,mz)+x(2,mx)*ya(2,my)
     +           *za(2,mz)+x(3,mx)*ya(3,my)*za(3,mz)+x(4,mx)*ya(4,my)
     +           *za(4,mz))
            g9 = g9 + ab*(x(1,mx)*y(1,my)*zdd(1,mz)+x(2,mx)*y(2,my)
     +           *zdd(2,mz)+x(3,mx)*y(3,my)*zdd(3,mz)+x(4,mx)*y(4,my)
     +           *zdd(4,mz))
 110     continue
         go to 160
 120     do 130 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*(xdd(1,mx)*y(1,my)*z(1,mz)+xdd(2,mx)*y(2,my)
     +           *z(2,mz)+xdd(3,mx)*y(3,my)*z(3,mz)+xdd(4,mx)*y(4,my)
     +           *z(4,mz)+xdd(5,mx)*y(5,my)*z(5,mz))
            g2 = g2 + ab*(xa(1,mx)*ya(1,my)*z(1,mz)+xa(2,mx)*ya(2,my)
     +           *z(2,mz)+xa(3,mx)*ya(3,my)*z(3,mz)+xa(4,mx)*ya(4,my)
     +           *z(4,mz)+xa(5,mx)*ya(5,my)*z(5,mz))
            g3 = g3 + ab*(xa(1,mx)*y(1,my)*za(1,mz)+xa(2,mx)*y(2,my)
     +           *za(2,mz)+xa(3,mx)*y(3,my)*za(3,mz)+xa(4,mx)*y(4,my)
     +           *za(4,mz)+xa(5,mx)*y(5,my)*za(5,mz))
            g5 = g5 + ab*(x(1,mx)*ydd(1,my)*z(1,mz)+x(2,mx)*ydd(2,my)
     +           *z(2,mz)+x(3,mx)*ydd(3,my)*z(3,mz)+x(4,mx)*ydd(4,my)
     +           *z(4,mz)+x(5,mx)*ydd(5,my)*z(5,mz))
            g6 = g6 + ab*(x(1,mx)*ya(1,my)*za(1,mz)+x(2,mx)*ya(2,my)
     +           *za(2,mz)+x(3,mx)*ya(3,my)*za(3,mz)+x(4,mx)*ya(4,my)
     +           *za(4,mz)+x(5,mx)*ya(5,my)*za(5,mz))
            g9 = g9 + ab*(x(1,mx)*y(1,my)*zdd(1,mz)+x(2,mx)*y(2,my)
     +           *zdd(2,mz)+x(3,mx)*y(3,my)*zdd(3,mz)+x(4,mx)*y(4,my)
     +           *zdd(4,mz)+x(5,mx)*y(5,my)*zdd(5,mz))
 130     continue
         go to 160
 140     do 150 n = 1 , ijkld2
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            ab = abdens(n)
            g1 = g1 + ab*(xdd(1,mx)*y(1,my)*z(1,mz)+xdd(2,mx)*y(2,my)
     +           *z(2,mz)+xdd(3,mx)*y(3,my)*z(3,mz)+xdd(4,mx)*y(4,my)
     +           *z(4,mz)+xdd(5,mx)*y(5,my)*z(5,mz)+xdd(6,mx)*y(6,my)
     +           *z(6,mz))
            g2 = g2 + ab*(xa(1,mx)*ya(1,my)*z(1,mz)+xa(2,mx)*ya(2,my)
     +           *z(2,mz)+xa(3,mx)*ya(3,my)*z(3,mz)+xa(4,mx)*ya(4,my)
     +           *z(4,mz)+xa(5,mx)*ya(5,my)*z(5,mz)+xa(6,mx)*ya(6,my)
     +           *z(6,mz))
            g3 = g3 + ab*(xa(1,mx)*y(1,my)*za(1,mz)+xa(2,mx)*y(2,my)
     +           *za(2,mz)+xa(3,mx)*y(3,my)*za(3,mz)+xa(4,mx)*y(4,my)
     +           *za(4,mz)+xa(5,mx)*y(5,my)*za(5,mz)+xa(6,mx)*y(6,my)
     +           *za(6,mz))
            g5 = g5 + ab*(x(1,mx)*ydd(1,my)*z(1,mz)+x(2,mx)*ydd(2,my)
     +           *z(2,mz)+x(3,mx)*ydd(3,my)*z(3,mz)+x(4,mx)*ydd(4,my)
     +           *z(4,mz)+x(5,mx)*ydd(5,my)*z(5,mz)+x(6,mx)*ydd(6,my)
     +           *z(6,mz))
            g6 = g6 + ab*(x(1,mx)*ya(1,my)*za(1,mz)+x(2,mx)*ya(2,my)
     +           *za(2,mz)+x(3,mx)*ya(3,my)*za(3,mz)+x(4,mx)*ya(4,my)
     +           *za(4,mz)+x(5,mx)*ya(5,my)*za(5,mz)+x(6,mx)*ya(6,my)
     +           *za(6,mz))
            g9 = g9 + ab*(x(1,mx)*y(1,my)*zdd(1,mz)+x(2,mx)*y(2,my)
     +           *zdd(2,mz)+x(3,mx)*y(3,my)*zdd(3,mz)+x(4,mx)*y(4,my)
     +           *zdd(4,mz)+x(5,mx)*y(5,my)*zdd(5,mz)+x(6,mx)*y(6,my)
     +           *zdd(6,mz))
 150     continue
      end if
 160  dd(ioff+1) = g1 + dd(ioff+1)
      dd(ioff+2) = g2 + dd(ioff+2)
      dd(ioff+3) = g3 + dd(ioff+3)
      dd(ioff+14) = g5 + dd(ioff+14)
      dd(ioff+15) = g6 + dd(ioff+15)
      dd(ioff+27) = g9 + dd(ioff+27)
      return
      end
      subroutine dr2as2(x,y,z,xa,ya,za,xb,yb,zb,xdd,ydd,zdd,ix,iy,
     1 iz,abdens,diag,ncdim)
c----------------------------
c     assembly routine case 2
c----------------------------
      implicit REAL  (a-h,o-z)
      logical diag
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*)
      dimension xa(ncdim,*),ya(ncdim,*),za(ncdim,*)
      dimension xb(ncdim,*),yb(ncdim,*), zb(ncdim,*)
      dimension xdd(ncdim,*),ydd(ncdim,*),zdd(ncdim,*)
      dimension ix(*),iy(*),iz(*),abdens(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/root)
INCLUDE(common/ddshln)
INCLUDE(common/incrdd)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     + aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     + aaa(9*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225),dd(144)
c
cvd$r assoc
c     i1 = inc1
      g1 = 0.0d0
      g2 = 0.0d0
      g3 = 0.0d0
      g5 = 0.0d0
      g6 = 0.0d0
      g9 = 0.0d0
      g4 = 0.0d0
      g7 = 0.0d0
      g8 = 0.0d0
      do 40 n = 1 , ijkld2
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         ab = abdens(n)
c------------------------------
         ss1 = 0.0d0
         ss2 = 0.0d0
         ss3 = 0.0d0
         ss5 = 0.0d0
         ss6 = 0.0d0
         ss9 = 0.0d0
         if (diag) then
            do 20 nr = 1 , ncontr
               ss1 = ss1 + xdd(nr,mx)*y(nr,my)*z(nr,mz)
               ss2 = ss2 + xa(nr,mx)*yb(nr,my)*z(nr,mz)
               ss3 = ss3 + xa(nr,mx)*y(nr,my)*zb(nr,mz)
               ss5 = ss5 + x(nr,mx)*ydd(nr,my)*z(nr,mz)
               ss6 = ss6 + x(nr,mx)*ya(nr,my)*zb(nr,mz)
               ss9 = ss9 + x(nr,mx)*y(nr,my)*zdd(nr,mz)
 20         continue
         else
            ss4 = 0.0d0
            ss7 = 0.0d0
            ss8 = 0.0d0
            do 30 nr = 1 , ncontr
               ss1 = ss1 + xdd(nr,mx)*y(nr,my)*z(nr,mz)
               ss2 = ss2 + xa(nr,mx)*yb(nr,my)*z(nr,mz)
               ss3 = ss3 + xa(nr,mx)*y(nr,my)*zb(nr,mz)
               ss5 = ss5 + x(nr,mx)*ydd(nr,my)*z(nr,mz)
               ss6 = ss6 + x(nr,mx)*ya(nr,my)*zb(nr,mz)
               ss9 = ss9 + x(nr,mx)*y(nr,my)*zdd(nr,mz)
               ss4 = ss4 + xb(nr,mx)*ya(nr,my)*z(nr,mz)
               ss7 = ss7 + xb(nr,mx)*y(nr,my)*za(nr,mz)
               ss8 = ss8 + x(nr,mx)*yb(nr,my)*za(nr,mz)
 30         continue
         end if
c-------------------------------
         g1 = ab*ss1 + g1
         g2 = ab*ss2 + g2
         g3 = ab*ss3 + g3
         g5 = ab*ss5 + g5
         g6 = ab*ss6 + g6
         g9 = ab*ss9 + g9
         if (.not.(diag)) then
            g4 = ab*ss4 + g4
            g7 = ab*ss7 + g7
            g8 = ab*ss8 + g8
         end if
 40   continue
      dd(ioff+1) = g1 + dd(ioff+1)
      dd(ioff+2) = g2 + dd(ioff+2)
      dd(ioff+3) = g3 + dd(ioff+3)
      dd(ioff+14) = g5 + dd(ioff+14)
      dd(ioff+15) = g6 + dd(ioff+15)
      dd(ioff+27) = g9 + dd(ioff+27)
      if (diag) return
      dd(ioff+13) = g4 + dd(ioff+13)
      dd(ioff+25) = g7 + dd(ioff+25)
      dd(ioff+26) = g8 + dd(ioff+26)
      return
      end
      subroutine dr2as3(x,y,z,xa,ya,za,xdd,ydd,zdd,ix,iy,iz,
     1 abdens, noform, ncdim)
c--------------------
c     assembly case 3
c--------------------
      implicit REAL  (a-h,o-z)
      logical noform
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*)
      dimension xa(ncdim,*),ya(ncdim,*),za(ncdim,*)
      dimension xdd(ncdim,*),ydd(ncdim,*),zdd(ncdim,*)
      dimension noform(*),ix(*),iy(*),iz(*), abdens(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/root)
      parameter (ncmax=65)
INCLUDE(common/ddshln)
INCLUDE(common/incrdd)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     + aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     + aaa(9*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225),dd(144)
c
cvd$r assoc
c
c     i1 = inc1
      g1 = 0.0d0
      g2 = 0.0d0
      g3 = 0.0d0
      g5 = 0.0d0
      g6 = 0.0d0
      g9 = 0.0d0
      n = 0
      nn = 0
      do 40 i = 1 , ijd2
         mmax = ik(i)
         do 30 k = 1 , mmax
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               dum = abdens(n)
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               ss1 = 0.0d0
               ss2 = 0.0d0
               ss3 = 0.0d0
               ss5 = 0.0d0
               ss6 = 0.0d0
               ss9 = 0.0d0
               do 20 nr = 1 , ncontr
                  d1 = ddij(nr,i)*ddkl(nr,k)
                  ss1 = ss1 + d1*xdd(nr,mx)*y(nr,my)*z(nr,mz)
                  ss2 = ss2 + d1*xa(nr,mx)*ya(nr,my)*z(nr,mz)
                  ss3 = ss3 + d1*xa(nr,mx)*y(nr,my)*za(nr,mz)
                  ss5 = ss5 + d1*x(nr,mx)*ydd(nr,my)*z(nr,mz)
                  ss6 = ss6 + d1*x(nr,mx)*ya(nr,my)*za(nr,mz)
                  ss9 = ss9 + d1*x(nr,mx)*y(nr,my)*zdd(nr,mz)
 20            continue
               g1 = dum*ss1 + g1
               g2 = dum*ss2 + g2
               g3 = dum*ss3 + g3
               g5 = dum*ss5 + g5
               g6 = dum*ss6 + g6
               g9 = dum*ss9 + g9
            end if
 30      continue
 40   continue
      dd(ioff+1) = g1 + dd(ioff+1)
      dd(ioff+2) = g2 + dd(ioff+2)
      dd(ioff+3) = g3 + dd(ioff+3)
      dd(ioff+14) = g5 + dd(ioff+14)
      dd(ioff+15) = g6 + dd(ioff+15)
      dd(ioff+27) = g9 + dd(ioff+27)
      return
      end
      subroutine dr2as4(x,y,z,xa,ya,za,xb,yb,zb,xdd,ydd,zdd,
     1 ix,iy,iz,abdens,noform,diag,ncdim)
c----------------------------
c     assembly routine case 4
c----------------------------
      implicit REAL  (a-h,o-z)
      logical diag,noform
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*)
      dimension xa(ncdim,*),ya(ncdim,*),za(ncdim,*)
      dimension xb(ncdim,*),yb(ncdim,*),zb(ncdim,*)
      dimension xdd(ncdim,*),ydd(ncdim,*),zdd(ncdim,*)
      dimension ix(*),iy(*),iz(*),abdens(*),noform(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
      parameter (ncmax=65)
INCLUDE(common/root)
INCLUDE(common/ddshln)
INCLUDE(common/incrdd)
      common/small/t1(ncmax),t2(ncmax),t3(ncmax),t4(ncmax),t5(ncmax),
     +   t6(ncmax),t7(ncmax),t8(ncmax),t9(ncmax)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225),dd(144)
c
cvd$r assoc
c
c     i1 = inc1
      g1 = 0.0d0
      g2 = 0.0d0
      g3 = 0.0d0
      g5 = 0.0d0
      g6 = 0.0d0
      g9 = 0.0d0
      n = 0
      nn = 0
c
c
      if (diag) then
c
c
         do 40 i = 1 , ijd2
            mmax = ik(i)
            do 30 k = 1 , mmax
               nn = nn + 1
               if (.not.(noform(nn))) then
                  n = n + 1
                  dum = abdens(n)
                  mx = ix(n)
                  my = iy(n)
                  mz = iz(n)
                  s1 = 0.0d0
                  s2 = 0.0d0
                  s3 = 0.0d0
                  s5 = 0.0d0
                  s6 = 0.0d0
                  s9 = 0.0d0
                  do 20 nr = 1 , ncontr
                     d1 = ddij(nr,i)*ddkl(nr,k)
                     s1 = s1 + d1*xdd(nr,mx)*y(nr,my)*z(nr,mz)
                     s2 = s2 + d1*xa(nr,mx)*yb(nr,my)*z(nr,mz)
                     s3 = s3 + d1*xa(nr,mx)*y(nr,my)*zb(nr,mz)
                     s5 = s5 + d1*x(nr,mx)*ydd(nr,my)*z(nr,mz)
                     s6 = s6 + d1*x(nr,mx)*ya(nr,my)*zb(nr,mz)
                     s9 = s9 + d1*x(nr,mx)*y(nr,my)*zdd(nr,mz)
 20               continue
                  g1 = dum*s1 + g1
                  g2 = dum*s2 + g2
                  g3 = dum*s3 + g3
                  g5 = dum*s5 + g5
                  g6 = dum*s6 + g6
                  g9 = dum*s9 + g9
               end if
 30         continue
 40      continue
         dd(ioff+1) = g1 + dd(ioff+1)
         dd(ioff+2) = g2 + dd(ioff+2)
         dd(ioff+3) = g3 + dd(ioff+3)
         dd(ioff+14) = g5 + dd(ioff+14)
         dd(ioff+15) = g6 + dd(ioff+15)
         dd(ioff+27) = g9 + dd(ioff+27)
         return
      else
c
         g4 = 0.0d0
         g7 = 0.0d0
         g8 = 0.0d0
c
         do 70 i = 1 , ijd2
            do 60 k = 1 , ik(i)
               nn = nn + 1
               if (.not.(noform(nn))) then
                  n = n + 1
                  dum = abdens(n)
                  mx = ix(n)
                  my = iy(n)
                  mz = iz(n)
                  s1 = 0.0d0
                  s2 = 0.0d0
                  s3 = 0.0d0
                  s4 = 0.0d0
                  s5 = 0.0d0
                  s6 = 0.0d0
                  s7 = 0.0d0
                  s8 = 0.0d0
                  s9 = 0.0d0
                  do 50 nr = 1 , ncontr
                     d1 = (ddij(nr,i)*ddkl(nr,k))
                     s1 = s1 + d1*xdd(nr,mx)*y(nr,my)*z(nr,mz)
                     s2 = s2 + d1*xa(nr,mx)*yb(nr,my)*z(nr,mz)
                     s3 = s3 + d1*xa(nr,mx)*y(nr,my)*zb(nr,mz)
                     s5 = s5 + d1*x(nr,mx)*ydd(nr,my)*z(nr,mz)
                     s6 = s6 + d1*x(nr,mx)*ya(nr,my)*zb(nr,mz)
                     s9 = s9 + d1*x(nr,mx)*y(nr,my)*zdd(nr,mz)
                     s4 = s4 + d1*xb(nr,mx)*ya(nr,my)*z(nr,mz)
                     s7 = s7 + d1*xb(nr,mx)*y(nr,my)*za(nr,mz)
                     s8 = s8 + d1*x(nr,mx)*yb(nr,my)*za(nr,mz)
 50               continue
                  g1 = dum*s1 + g1
                  g2 = dum*s2 + g2
                  g3 = dum*s3 + g3
                  g5 = dum*s5 + g5
                  g6 = dum*s6 + g6
                  g9 = dum*s9 + g9
                  g4 = dum*s4 + g4
                  g7 = dum*s7 + g7
                  g8 = dum*s8 + g8
               end if
 60         continue
 70      continue
         dd(ioff+1) = g1 + dd(ioff+1)
         dd(ioff+2) = g2 + dd(ioff+2)
         dd(ioff+3) = g3 + dd(ioff+3)
         dd(ioff+14) = g5 + dd(ioff+14)
         dd(ioff+15) = g6 + dd(ioff+15)
         dd(ioff+27) = g9 + dd(ioff+27)
         dd(ioff+13) = g4 + dd(ioff+13)
         dd(ioff+25) = g7 + dd(ioff+25)
         dd(ioff+26) = g8 + dd(ioff+26)
         return
      end if
c
      end
      subroutine dr2end(der2,dd,dacc,index,ii,jj,nshell,isec46)
c-----------------------
c     close down routine
c-----------------------
c
      implicit REAL  (a-h,o-z)
      dimension der2(*),dd(*),dacc(*)
c
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/timez)
      character*10 charwall
c
      nlen = lds(isect(46))
_IF(parallel)
c
c      write(6,*)'dr2end',ipg_nodeid(), ii,jj, index
c
c index = 0 timout
c index = 1 finish
c index = 2 accumulate
c
c#accum a cycle of 2e- drv
c
c     add contents of present cycle
      call vadd(dacc,1,der2,1,dacc,1,nlen)
c
c     clear buffer
      call vclr(der2,1,nlen)
c
      if(index.eq.0)then
c
         call caserr('out of time')
c
      elseif (index.eq.1)then
c
c sum accum buffer from each processor
c
c         write(6,*)'gop',ipg_nodeid()
         call pg_dgop(130, dacc, nlen,'+')
c
c  accumulate with disk contents
c
         call rdedx(dd,nlen,isec46,ifild)
         call vadd(dd,1,dacc,1,dd,1,nlen)
         call wrt3(dd,nlen,isec46,ifild)
c
c read result into der2
c
         call rdedx(der2,nlen,isec46,ifild)
         ist = 1
         jst = 1
         kst = 1
         lst = 1
         nrec = 1
         intloc = 1
         irest = 14
         call revise
         call clredx
         write (iwr,6010) cpulft(1)
      elseif (index.eq.2)then
c
c         write(6,*)'accum',ipg_nodeid()
c
      endif
      return
_ELSE

c     read what has been calculated so far
      call rdedx(dd,nlen,isec46,ifild)
c
c#accum a cycle of 2e- drv
c
c     add contents of present cycle
      call vadd(der2,1,dd,1,der2,1,nlen)

c     write out again
      call wrt3(der2,nlen,isec46,ifild)

c     clear buffer
      call vclr(der2,1,nlen)
      if (index.ne.2) then
         if (index.eq.1) then
c
c     all finished
c
            call rdedx(der2,nlen,isec46,ifild)
            ist = 1
            jst = 1
            kst = 1
            lst = 1
            nrec = 1
            intloc = 1
            irest = 14
            call revise
            call clredx
            write (iwr,6010) cpulft(1) ,charwall()
            return
         else
c     ----- dump restart data -----
            if (jj.eq.nshell) return
            ist = ii
            jst = jj + 1
            if (jj.eq.ii) ist = ii + 1
            if (jj.eq.ii) jst = 1
            kst = 1
            lst = 1
            nrec = 1
            intloc = 1
            irest = 13
            call revise
            write (iwr,6020)
            write (iwr,6030) ist , jst , kst , lst
            call clenms('stop')
         end if
      end if
c
c   precautionary dump
c
      if (jj.eq.nshell) return
      istr = ist
      jstr = jst
      kstr = kst
      lstr = lst
      ist = ii
      jst = jj + 1
      if (jj.eq.ii) ist = ii + 1
      if (jj.eq.ii) jst = 1
      kst = 1
      lst = 1
      call revise
      call clredx
      ist = istr
      jst = jstr
      kst = kstr
      lst = lstr
      return
_ENDIF

 6010 format (/' end of 2-electron derivative integral evaluation at ',
     +        f8.2,' seconds',a10,' wall')
 6020 format (//10x,'insufficient time to complete evaluation of',
     +        ' two-electron contibution')
 6030 format (//10x,'restart at shells i=',i4,'  j=',i4,'  k=',i4,
     +        '  l=',i4)
      end
      subroutine dr2gen(qq,iqq,noform)
c-----------------------------------------------
c     contains all the loops over the primitives
c-----------------------------------------------
      implicit REAL  (a-h,o-z)
      logical noform
      dimension qq(*)
      dimension iqq(*),noform(*)
c
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/common)
INCLUDE(common/maxlen)
INCLUDE(common/iofile)
INCLUDE(common/root)
INCLUDE(common/ddmisc)
INCLUDE(common/prnprn)
INCLUDE(common/ddshln)
INCLUDE(common/ddshli)
      parameter (ncmax=65)
INCLUDE(common/setdd)
INCLUDE(common/incrdd)
      common/bufb/axak(mxp2),ayak(mxp2),azak(mxp2),axai(mxp2),
     1  ayai(mxp2),azai(mxp2),abv(mxp2),aandbv(mxp2),rhov(mxp2),
     2  xxv(mxp2),c1xv(mxp2),c2xv(mxp2),c3xv(mxp2),c4xv(mxp2),
     3  c1yv(mxp2),c2yv(mxp2),c3yv(mxp2),c4yv(mxp2),
     4  c1zv(mxp2),c2zv(mxp2),c3zv(mxp2),c4zv(mxp2),expev(mxp2)
c
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     +   aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +   aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +   dd(4*mxp2),ijden(225),ik(225),
     +   ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +   dij(225),dkl(225),ijgt(225),klgt(225),ddtemp(144)
c
      logical spij,spkl
      logical double
      logical sptru
      data one/1.0d0/
      data pi252/34.986836655250d0/
c
      call dr2set(nimax,njmax,nkmax,nlmax,lit,ljt,lkt,llt,skip)
      kln2 = 1
      kln1 = nlmax + 1
      ijn2 = kln1*(nkmax+1)
      ijn1 = ijn2*(njmax+1)
      inc1 = ijn1*(nimax+1)
c     if((inc1/4)*4.eq.inc1)inc1=inc1+1
c
c     offsets for subsiduary integrals storage
c
      ic0 = iabd
      ic1 = isub
      isubsp = maxqq - isub - 1
      ncmmm = isubsp/(12*inc1)
      ncnnn = ncmax - 1
      if (ncmmm.gt.ncnnn) ncmmm = ncnnn
      ncmmm = (ncmmm/nroots)*nroots
c
      ic2 = ic1 + inc1*ncmmm
      ic3 = ic2 + inc1*ncmmm
      ic4 = ic3 + inc1*ncmmm
      ic5 = ic4 + inc1*ncmmm
      ic6 = ic5 + inc1*ncmmm
      ic7 = ic6 + inc1*ncmmm
      ic8 = ic7 + inc1*ncmmm
      ic9 = ic8 + inc1*ncmmm
      ic10 = ic9 + inc1*ncmmm
      ic11 = ic10 + inc1*ncmmm
      ic12 = ic11 + inc1*ncmmm
c
      if (nimax.lt.njmax) then
         is = nimax
         nimax = njmax
         njmax = is
         ij1x = ijn2
         ij2x = ijn1
         xc = xj
         yc = yj
         zc = zj
         dxij = xj - xi
         dyij = yj - yi
         dzij = zj - zi
      else
         ij1x = ijn1
         ij2x = ijn2
         xc = xi
         yc = yi
         zc = zi
         dxij = xi - xj
         dyij = yi - yj
         dzij = zi - zj
      end if
      if (nkmax.lt.nlmax) then
         is = nlmax
         nlmax = nkmax
         nkmax = is
         kl1x = kln2
         kl2x = kln1
         xd = xl
         yd = yl
         zd = zl
         dxkl = xl - xk
         dykl = yl - yk
         dzkl = zl - zk
      else
         xd = xk
         yd = yk
         zd = zk
         dxkl = xk - xl
         dykl = yk - yl
         dzkl = zk - zl
         kl1x = kln1
         kl2x = kln2
      end if
      nmax = nimax + njmax
      mmax = nkmax + nlmax
      maxn = nmax + 1
      do 20 i = 1 , maxn
         n = i - 1
         if (n.le.nimax) iorg(i) = ij1x*n + 1
         if (n.gt.nimax) iorg(i) = ij1x*nimax + ij2x*(n-nimax) + 1
 20   continue
      maxm = mmax + 1
      do 30 k = 1 , maxm
         n = k - 1
         if (n.le.nkmax) korg(k) = kl1x*n
         if (n.gt.nkmax) korg(k) = kl1x*nkmax + kl2x*(n-nkmax)
 30   continue
      ncontr = 0
c
c     indexing
c
      call indexa(ijx,ijy,ijz,ijd2,mini,maxi,minj,maxj,iandj,
     +            ijn1,ijn2,1)
      call indexa(klx,kly,klz,kld2,mink,maxk,minl,maxl,kandl,
     +            kln1,kln2,0)
      nn = 0
      ijkld2 = 0
      do 50 i = 1 , ijd2
         maxik = ik(i)
         do 40 k = 1 , maxik
            nn = nn + 1
            if (.not.(noform(nn))) then
               ijkld2 = ijkld2 + 1
               iqq(ixi+ijkld2-1) = ijx(i) + klx(k)
               iqq(iyi+ijkld2-1) = ijy(i) + kly(k)
               iqq(izi+ijkld2-1) = ijz(i) + klz(k)
            end if
 40      continue
 50   continue
c
      do 60 n = 1 , nij
         axak(n) = aa(n)*(x1(n)-xd)
         ayak(n) = aa(n)*(y1(n)-yd)
         azak(n) = aa(n)*(z1(n)-zd)
         axai(n) = aa(n)*(x1(n)-xc)
         ayai(n) = aa(n)*(y1(n)-yc)
         azai(n) = aa(n)*(z1(n)-zc)
 60   continue
c
      abmax = 0.0d0
      do 70 i = 1 , ijkld2
         ab = dabs(qq(iabd+i-1))
         if (ab.gt.abmax) abmax = ab
 70   continue
      spkl = (mink.eq.1 .and. maxk.eq.4) .or.
     +       (minl.eq.1 .and. maxl.eq.4)
      spij = (mini.eq.1 .and. maxi.eq.4) .or.
     +       (minj.eq.1 .and. maxj.eq.4)
      sptru = spkl .or. spij
c
c
c     k primitive
c
      maxlg = ngd
      do 300 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*xk
         akyk = ak*yk
         akzk = ak*zk
         csk = csc(kg)*pi252
         cpk = cpc(kg)*pi252
         cdk = cdc(kg)*pi252
         cfk = cfc(kg)*pi252
c
c     ----- l primitive
c
         if (kandl) maxlg = kg
         do 290 lg = 1 , maxlg
            al = dg(lg)
            b = ak + al
            binv = one/b
            bbrrk = al*brrk*binv
            if ((bbrrk+rsmall).le.tol(1)) then
               exkl = dexp(-bbrrk)
               csl = csd(lg)*binv
               cpl = cpd(lg)*binv
               cdl = cdd(lg)*binv
               cfl = cfd(lg)*binv
               xb = (akxk+al*xl)*binv
               yb = (akyk+al*yl)*binv
               zb = (akzk+al*zl)*binv
               bxbk = b*(xb-xd)
               bybk = b*(yb-yd)
               bzbk = b*(zb-zd)
               bxbi = b*(xb-xc)
               bybi = b*(yb-yc)
               bzbi = b*(zb-zc)
c
c     ----- density factor
c
               double = kandl .and. kg.ne.lg
               n = 0
               lmax = maxl
               do 190 k = mink , maxk
                  if (kandl) lmax = k
                  go to (80,90,120,120,100,120,120,120,120,120,110,120,
     +                   120,120,120,120,120,120,120,120) , k
 80               dum1 = csk
                  go to 120
 90               dum1 = cpk
                  go to 120
 100              dum1 = cdk
                  go to 120
 110              dum1 = cfk
 120              do 180 l = minl , lmax
                     go to (130,140,170,170,150,170,170,170,170,170,160,
     +                      170,170,170,170,170,170,170,170,170) , l
 130                 dum2 = dum1*csl
                     if (double) then
                        if (k.gt.1) then
                           dum2 = dum2 + csk*cpl
                        else
                           dum2 = dum2 + dum2
                        end if
                     end if
                     go to 170
 140                 dum2 = dum1*cpl
                     if (double) dum2 = dum2 + dum2
                     go to 170
 150                 dum2 = dum1*cdl
                     if (double) dum2 = dum2 + dum2
                     go to 170
 160                 dum2 = dum1*cfl
                     if (double) dum2 = dum2 + dum2
 170                 n = n + 1
                     dkl(n) = dum2
 180              continue
 190           continue
               dkld = dkl(1)
               if (sptru) then
                  do 200 k = 1 , n
                     dkl(k) = dkl(k)/dkld
 200              continue
               end if
               dkld = dkld*exkl
               if (dabs(dkld*abmax).ge.tol(3)) then
c
c     ----- pair of i,j primitives
                  do 210 n = 1 , nij
                     abv(n) = aa(n)*b
                     aandbv(n) = aa(n) + b
                     expev(n) = exij(n)/dsqrt(aa(n)+b)
                     rhov(n) = abv(n)/aandbv(n)
                     xxv(n) = rhov(n)
     +                        *((x1(n)-xb)**2+(y1(n)-yb)**2+(z1(n)-zb)
     +                        **2)
                     c1xv(n) = bxbk + axak(n)
                     c2xv(n) = bxbk*aa(n)
                     c3xv(n) = bxbi + axai(n)
                     c4xv(n) = b*axai(n)
                     c1yv(n) = bybk + ayak(n)
                     c2yv(n) = bybk*aa(n)
                     c3yv(n) = bybi + ayai(n)
                     c4yv(n) = b*ayai(n)
                     c1zv(n) = bzbk + azak(n)
                     c2zv(n) = bzbk*aa(n)
                     c3zv(n) = bzbi + azai(n)
                     c4zv(n) = b*azai(n)
 210              continue
c
                  n = 0
                  nn = 0
                  jgmax = ngb
                  do 280 ig = 1 , nga
                     ai = ag(ig)
                     if (iandj) jgmax = ig
                     do 270 jg = 1 , jgmax
                        n = n + 1
                        if ((bbrrk+r(n)).lt.tol(2)) then
                           aj = bg(jg)
                           dijd = dd(nn+1)
                           if (sptru) then
                              dddd = one/dijd
                              do 220 i = 1 , ijd2
                                 dij(i) = dd(ijden(i)+nn)*dddd
 220                          continue
                           end if
                           expe = dkld*dijd*expev(n)
                           if (dabs(expe*abmax).ge.tol(4)) then
                              pp = xxv(n)
c
c     ----- roots and weights for quadrature
c
                              if (nroots.le.3) call rt123
                              if (nroots.eq.4) call roots4
                              if (nroots.eq.5) call roots5
                              if (nroots.gt.5) call rootss
c                             mm = 0
***                           max = nmax + 1
c
c     compute two-electron  integrals for each root
c
                              nnn0 = ncontr
                              do 230 m = 1 , nroots
                               ncontr = ncontr + 1
                               u2 = u(m)*rhov(n)
                               f00(ncontr) = expe*w(m)
                               dum2 = 0.5d0/(abv(n)+u2*aandbv(n))
                               dum = dum2 + dum2
                               bp01(ncontr) = (aa(n)+u2)*dum2
                               b00(ncontr) = u2*dum2
                               b10(ncontr) = (b+u2)*dum2
                               xcp00(ncontr) = (u2*c1xv(n)+c2xv(n))*dum
                               xc00(ncontr) = (u2*c3xv(n)+c4xv(n))*dum
                               ycp00(ncontr) = (u2*c1yv(n)+c2yv(n))*dum
                               yc00(ncontr) = (u2*c3yv(n)+c4yv(n))*dum
                               zcp00(ncontr) = (u2*c1zv(n)+c2zv(n))*dum
                               zc00(ncontr) = (u2*c3zv(n)+c4zv(n))*dum
                               aei(ncontr) = ai
                               aej(ncontr) = aj
                               aek(ncontr) = ak
                               ael(ncontr) = al
 230                          continue
                              if (sptru) then
                                 ncontr = nnn0
                                 do 260 m = 1 , nroots
                                    ncontr = ncontr + 1
                                    do 240 iii = 1 , ijd2
                                       ddij(ncontr,iii) = dij(iii)
 240                                continue
                                    do 250 iii = 1 , kld2
                                       ddkl(ncontr,iii) = dkl(iii)
 250                                continue
 260                             continue
                              end if
                              if (ncontr.ge.ncmmm) then
c
c     ----- form (i,j//k,l) integrals
c
                       call dr2xyz(qq(ic1),qq(ic2),qq(ic3),ncmmm)
                       if (.not.(skip(1))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aei,ljt,lkt,llt,lit+1,ijn2,kln1,kln2,ijn1,ncmmm)
                          call dr2pri(qq(ic4),qq(ic5),qq(ic6),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          ioff = noff(1)
                          if (.not.sptru)
     +                        call dr2as1(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
                          if (sptru)
     +                        call dr2as3(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
                       end if
                       if (.not.(skip(2))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  aej,lit+1,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          ioff = noff(2)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(2),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(2),ncmmm)
                       end if
                       if (.not.(skip(3))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  aek,lit+1,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          ioff = noff(3)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(3),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(3),ncmmm)
                       end if
                       if (.not.(skip(4))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  ael,lit+1,ljt,lkt,llt,ijn1,ijn2,kln1,kln2,
     +  ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
                          ioff = noff(4)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(4),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(4),ncmmm)
                       end if
                       if (.not.(skip(5))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aej,lit,lkt,llt,ljt+1,ijn1,kln1,kln2,ijn2,ncmmm)
                          call dr2pri(qq(ic4),qq(ic5),qq(ic6),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          ioff = noff(5)
                          if (.not.sptru)
     +                        call dr2as1(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
                          if (sptru)
     +                        call dr2as3(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
                       end if
                       if (.not.(skip(6))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  aek,lit,ljt+1,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          ioff = noff(6)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(6),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(6),ncmmm)
                       end if
                       if (.not.(skip(7))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  ael,lit,ljt+1,lkt,llt,ijn1,ijn2,kln1,kln2,ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
                          ioff = noff(7)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(7),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(7),ncmmm)
                       end if
                       if (.not.(skip(8))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aek,lit,ljt,llt,lkt+1,ijn1,ijn2,kln2,kln1,ncmmm)
                          call dr2pri(qq(ic4),qq(ic5),qq(ic6),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
                          ioff = noff(8)
                          if (.not.sptru)
     +                        call dr2as1(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
                          if (sptru)
     +                        call dr2as3(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
                       end if
                       if (.not.(skip(9))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic7),qq(ic8),qq(ic9),
     +  ael,lit,ljt,lkt+1,llt,ijn1,ijn2,kln1,kln2,ncmmm)
                          call dr2pri(qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
                          ioff = noff(9)
                          if (.not.sptru)
     +                        call dr2as2(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(9),ncmmm)
                          if (sptru)
     +                        call dr2as4(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(9),ncmmm)
                       end if
                       if (.not.(skip(10))) then
                          call dr2pri(qq(ic1),qq(ic2),qq(ic3),
     +  qq(ic4),qq(ic5),qq(ic6),
     +  ael,lit,ljt,lkt,llt+1,ijn1,ijn2,kln1,kln2,ncmmm)
                          call dr2pri(qq(ic4),qq(ic5),qq(ic6),
     +  qq(ic10),qq(ic11),qq(ic12),
     +  ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,kln2,ncmmm)
                          ioff = noff(10)
                          if (.not.sptru)
     +                        call dr2as1(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
                          if (sptru)
     +                        call dr2as3(qq(ic1),qq(ic2),
     +  qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +  iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
                       end if
c
                       ncontr = 0
                    end if
                 end if
              end if
c
c     end of loops over primitives
c
                        nn = nn + 4
 270                 continue
 280              continue
               end if
            end if
 290     continue
 300  continue
      if (ncontr.ne.0) then
         call dr2xyz(qq(ic1),qq(ic2),qq(ic3),ncmmm)
         if (.not.(skip(1))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aei,ljt,lkt,llt,lit+1,ijn2,kln1,kln2,ijn1,ncmmm)
            call dr2pri(qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),
     +      qq(ic12),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            ioff = noff(1)
            if (.not.sptru) call dr2as1(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
            if (sptru) call dr2as3(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
         end if
         if (.not.(skip(2))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      aej,lit+1,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            ioff = noff(2)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +        qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(2),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +        qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(2),ncmmm)
         end if
         if (.not.(skip(3))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      aek,lit+1,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            ioff = noff(3)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +        qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(3),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +        qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(3),ncmmm)
         end if
         if (.not.(skip(4))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      ael,lit+1,ljt,lkt,llt,ijn1,ijn2,kln1,kln2,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,ijn1,ncmmm)
            ioff = noff(4)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +       qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +       qq(ic10),qq(ic11),qq(ic12),
     +       iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(4),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +        qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(4),ncmmm)
         end if
         if (.not.(skip(5))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aej,lit,lkt,llt,ljt+1,ijn1,kln1,kln2,ijn2,ncmmm)
            call dr2pri(qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            ioff = noff(5)
            if (.not.sptru) call dr2as1(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
            if (sptru) call dr2as3(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
         end if
         if (.not.(skip(6))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      aek,lit,ljt+1,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            ioff = noff(6)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +        qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(6),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +         qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +         qq(ic11),qq(ic12),
     +         iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(6),ncmmm)
         end if
         if (.not.(skip(7))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      ael,lit,ljt+1,lkt,llt,ijn1,ijn2,kln1,kln2,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,ijn2,ncmmm)
            ioff = noff(7)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +        qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(7),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +        qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(7),ncmmm)
         end if
         if (.not.(skip(8))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aek,lit,ljt,llt,lkt+1,ijn1,ijn2,kln2,kln1,ncmmm)
            call dr2pri(qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
            ioff = noff(8)
            if (.not.sptru) call dr2as1(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
            if (sptru) call dr2as3(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
         end if
         if (.not.(skip(9))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic7),qq(ic8),
     +                  qq(ic9),
     +      ael,lit,ljt,lkt+1,llt,ijn1,ijn2,kln1,kln2,ncmmm)
            call dr2pri(qq(ic7),qq(ic8),qq(ic9),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,kln1,ncmmm)
            ioff = noff(9)
            if (.not.sptru) call dr2as2(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),
     +        qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),diag(9),ncmmm)
            if (sptru) call dr2as4(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic7),qq(ic8),qq(ic9),qq(ic10),
     +        qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,diag(9),ncmmm)
         end if
         if (.not.(skip(10))) then
            call dr2pri(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                  qq(ic6),
     +      ael,lit,ljt,lkt,llt+1,ijn1,ijn2,kln1,kln2,ncmmm)
            call dr2pri(qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),
     +                  qq(ic12),
     +      ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,kln2,ncmmm)
            ioff = noff(10)
            if (.not.sptru) call dr2as1(qq(ic1),qq(ic2),qq(ic3),
     +        qq(ic4),qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),ncmmm)
            if (sptru) call dr2as3(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +        qq(ic5),qq(ic6),qq(ic10),qq(ic11),qq(ic12),
     +        iqq(ixi),iqq(iyi),iqq(izi),qq(ic0),noform,ncmmm)
         end if
      end if
      ncontr = 0
c
c
c     use translational invariance to find missing triangles
c
c
      call dr2trn(ddtemp,ncen,natom,iwr,odebug(7))
      return
      end

      subroutine dr2int(qq,iqq,iso,nshels,noform,isec46)
c-----------------------------------------------
c  second derivatives of integrals
c  driving routine for two-electron contribution
c  to second derivatives
c-----------------------------------------------
      implicit REAL  (a-h,o-z)
      logical noform
      dimension noform(*),qq(*),iqq(*),iso(*)
c
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/ghfblk)
INCLUDE(common/grad2)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/maxlen)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/timez)
INCLUDE(common/mapper)
INCLUDE(common/ddmisc)
INCLUDE(common/ddshln)
INCLUDE(common/incrdd)
INCLUDE(common/prnprn)
INCLUDE(common/runopt)
c
      logical dpres,gpres,lab,labc,labcd,dtru
      common/blkin/zzzzz(510),mword,nlenx,kworx,kworxx
      common/craypk/labs(1360)
c
c     ncmax is (one more than) maximum length of vector loop
c     in integral assembly
c
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225),ddtemp(144)
_IF(parallel)
INCLUDE(common/parallel)
_ENDIF
c
c     noform onwards to end of block is dynamically allocated
c     maximum length of block is maxq
c
c
      dimension m0(48),m1(48),m2(48),m3(48)
      character*10 charwall

      character *8 oscf,grhf
      data oscf,grhf/'oscf','grhf'/
      data zero/0.0d0/
c
c
      nav = lenwrd()
      maxqq = maxq

      kloc(nshell+1) = num + 1
      call setsto(1360,0,labs)
      call timit(3)
      tim0 = tim
      tim1 = tim
c     tim2 = tim
      dtmax = zero
      nat3 = nat*3
      nlen = nat3*nat3
      if (scftyp.eq.grhf) then
         call secget(isect(53),53,iblok)
         call readi(nact,lds(isect(53))*nav,iblok,ifild)
      end if
      if (irest.lt.14) then
c
c     space requirements
c
         dpres = .false.
         fpres = .false.
         gpres = .false.
         do 20 i = 1 , nshell
            if (ktype(i).eq.3) dpres = .true.
            if (ktype(i).eq.4) fpres = .true.
            if (ktype(i).eq.5) gpres = .true.
 20      continue
c
c     initial core for indexing arrays
c
         lnddm = 257
         if (dpres) lnddm = 1297
         if (fpres) lnddm = 10001
         if (gpres) lnddm = 50625
         nreq = lenint(4)*lnddm + nw196(5) + lenint(nt*nat)
         iof = lenrel(nw196(5))
         ioff1 = iof + nt*nat
c
c     id0 is where second-derivative matrix is accumulated
c
_IF(parallel)
         id00 = nreq + 1
         id0 = id00 + nlen
_ELSE
         id0 = nreq + 1
         id00 = id0
_ENDIF
         id1 = id0 + nlen
c
c     density matrices start at id1
c
         ndens = 1
         if (scftyp.eq.oscf) ndens = 2
         if (oss) ndens = 3
         if (scftyp.eq.grhf) ndens = njk
         nreq = nreq + nlen + nx*ndens
_IF(parallel)
         nreq = nreq + nlen
_ENDIF
c
         id2 = id1 + nx
         if (scftyp.eq.grhf) id2 = id1 + nx*njk
         id3 = id2 + nx
c
c     2-pdm at iabd
c
         iabd = nreq + 1
c
c     allow space for derivative subsidiary integrals
c
         nreq = iabd + lnddm + 4860
         if (dpres) nreq = iabd + lnddm + 18504
         if (fpres) nreq = iabd + lnddm + 67500
c        g-requirements to be determined ...
         if (gpres) nreq = iabd + lnddm + 67500
         write (iwr,6010) cpulft(1) ,charwall()
         write (iwr,6020) maxqq , nreq
         if (nreq.gt.maxqq) then
            call caserr('seq2e: insufficient core available')
         end if
         call vclr(qq(id0),1,nlen)
         call dr2ini(qq(id1),qq(id2),qq(id3),qq(iabd))
         mword = 0
         kworx = 999
         iword = 0
         call rdedx(qq(1),nw196(5),ibl196(5),ifild)
         if (mp2) call search(iblk2d,ifil2d)
c
c          start start start start
c
_IF(parallel)
c
c
c   zero common segment
c
      nlen = lds(isect(46))
      call vclr(qq(id00),1,nlen)
      call pg_dlbreset
      mcount=0
      next = ipg_dlbtask()
_ENDIF

         if (ist.le.nshell) then
c     ----- ishell -----
_IF(parallel)
            do 180 ii = nshell , ist , -1
_ELSE
            do 180 ii = ist , nshell
_ENDIF
               do 40 it = 1 , nt
                  id = iso(ii+iliso(it))
                  if (id.gt.ii) go to 180
                  m0(it) = id
 40            continue
               iceni = katom(ii)
c     ----- jshell -----
               j0 = jst
_IF(parallel)
               do 170 jj = ii , j0 , -1
_ELSE
               do 170 jj = j0 , ii
_ENDIF
                  jst = 1
                  do 60 it = 1 , nt
                     id = m0(it)
                     jd = iso(jj+iliso(it))
                     if (jd.gt.ii) go to 170
                     if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (id.eq.ii .and. jd.gt.jj) go to 170
                     m1(it) = id
                     m2(it) = jd
 60               continue
                  lab = katom(jj).eq.iceni
c
c     store information about the pair (ij)
c
                  call dr2shl(1,ii,jj,kk,ll)
                  call dr2prm()
                  if (nij.ne.0) then
_IF(parallel)
                  mcount = mcount + 1
cdb
c                  write(6,*)'node=',ipg_nodeid(),'mcount=',mcount,
c     &                 'next=',next
                  if(mcount .eq. next) then
_ENDIF
c
c     ----- kshell -----
                     k0 = kst
                     do 160 kk = k0 , ii
                        kst = 1
                        do 80 it = 1 , nt
                           kd = iso(kk+iliso(it))
                           if (kd.gt.ii) go to 160
                           m3(it) = kd
 80                     continue
                        labc = lab .and. katom(kk).eq.iceni
c     ----- lshell -----


                        l0 = lst
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                        do 150 ll = l0 , maxll
                           lst = 1
                           labcd = labc .and. katom(ll).eq.iceni
                           if (.not.(labcd)) then
                              n4 = 0
                              do 100 it = 1 , nt
                                 ld = iso(ll+iliso(it))
                                 if (ld.gt.ii) go to 150
                                 kd = m3(it)
                                 if (kd.lt.ld) then
                                    nd = kd
                                    kd = ld
                                    ld = nd
                                 end if
                                 id = m1(it)
                                 jd = m2(it)
                                 if (id.eq.ii .or. kd.eq.ii) then
                                    if (kd.ge.id) then
                                       if (kd.ne.id .or. ld.gt.jd) then
                                         nd = id
                                         id = kd
                                         kd = nd
                                         nd = jd
                                         jd = ld
                                         ld = nd
                                       end if
                                    end if
                                    if (jd.ge.jj) then
                                       if (jd.gt.jj) go to 150
                                       if (kd.ge.kk) then
                                         if (kd.gt.kk) go to 150
                                         if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 150
                                         n4 = n4 + 1
                                         end if
                                       end if
                                    end if
                                 end if
 100                          continue
c     ----- calculate q4 factor for this group of shells -----
                              qq4 = dfloat(nt)/dfloat(n4)
c     ----- check for redundant combinations -----
                              call dr2skp(ii,jj,kk,ll,iwr,odebug(7))
                              if (npass.ne.0) then
                                 call vclr(ddtemp,1,144)
                                 call dr2shl(2,ii,jj,kk,ll)
c     ----- form products of density matrix elements -----
                                 if (mp2)
     +                               call dr2pdm(qq(iabd),ii,jj,kk,ll,
     +  qq4)
                                 call dr2den(qq(iabd),qq(id1),qq(id2),
     +  qq(id3),iky,ii,jj,kk,ll,qq4,scftyp,oss,mp2)
c     ----mess about with the density matix by eliminating
c     zero elements
c     absorb normalisation of d-functions into density matrix
                                 dtru = lit.ge.3 .or. ljt.ge.3 .or.
     +                                  lkt.ge.3 .or. llt.ge.3
                                 if (dtru .and. norm)
     +                               call dr2ndf(qq(iabd))
                                 nn = ioff1
                                 do 120 i = 1 , ijd2
                                    mix = ik(i)
                                    do 110 k = 1 , mix
                                       nn = nn + 1
                                       n = ijgt(i) + klgt(k)
                                       noform(nn) = .false.
                                       if (dabs(qq(iabd-1+n)).lt.cutoff)
     +                                    noform(nn) = .true.
 110                                continue
 120                             continue
                        call dcopy(lendd,qq(iabd),1,qq(iabd+lendd),1)
                                 nn = ioff1
                                 ijkld2 = 0
                                 do 140 i = 1 , ijd2
                                    mix = ik(i)
                                    do 130 k = 1 , mix
                                       nn = nn + 1
                                       if (.not.(noform(nn))) then
                                         n = ijgt(i) + klgt(k)
                                         ijkld2 = ijkld2 + 1
                                         qq(iabd-1+ijkld2)
     +                                      = qq(iabd+lendd-1+n)
                                       end if
 130                                continue
 140                             continue
                                 if (ijkld2.ne.0) then
                                    call dr2gen(qq,iqq(ioff1+1),
     +  noform(ioff1+1))
                                    call dr2frm(qq(id0),ddtemp,natom,
     +  nat*3)
                                 end if
                              end if
                           end if
 150                    continue
 160                 continue
                     call timit(3)
                     tim3 = tim - tim1
                     tim1 = tim
                     if (tim3.gt.dtmax) dtmax = tim3
                     if ((timlim-tim).le.(dtmax*1.8d0+10.0d0)) then
                        call dr2end(qq(id0),qq(iabd),qq(id00),
     +                              0,ii,jj,nshell,isec46)
                     end if
                     call dr2end(qq(id0),qq(iabd),qq(id00),
     +                           2,ii,jj,nshell,isec46)
_IF(parallel)
c                  write(6,*)'get task',ipg_nodeid()
                  next = ipg_dlbtask()
c                  write(6,*)' done ',next, ipg_nodeid()

                  endif
_ENDIF
                  end if
 170           continue
               call timit(3)
c              dtim = tim - tim2
c              tim2 = tim
 180        continue
         end if

c     ----- end of 'shell' loops -----

         call dr2end(qq(id0),qq(iabd),qq(id00),
     +               1,ii,jj,nshell,isec46)
         call dr2sym(qq(id0),qq(iabd),iso(1),iso(iof+1),nat,nat3,
     +               nshels)
_IF(parallel)
         call pg_dlbpush()
_ENDIF
c
c... for frozen atoms in hessian  - make hessian big
c
      big = 1000000.0d0
      do j=1,nat 
         if (zopti(j).eq.'no') then
            do k=1,nat*3
               do jc=(j-1)*3+1,j*3
                   qq(k+(jc-1)*nat*3+id0-1) = 0.0d0
                   qq(jc+(k-1)*nat*3+id0-1) = 0.0d0
                   if (k.eq.jc) qq(k+(jc-1)*nat*3+id0-1) = big
               end do  
            end do  
         end if  
      end do  
c
c#accum ??? rewrite after symmetrisation???
         call wrt3(qq(id0),nlen,isec46,ifild)
         if (oprn(23)) then
            write (iwr,6030)
            call prnder(qq(id0),nat3,iwr)
         endif
      end if
      ncc = nat3*(nat3+1)/2
      ii0 = 1
      ii1 = nlen + ii0
      ii2 = ncc + ii1
      ii3 = ncc + ii2
      ii4 = ii3 + nat3
      last = ii4 + nat3
      call rdedx(qq(ii0),nlen,isec46,ifild)
      call secget(isect(14),14,isec14)
      call rdedx(de,maxat*3,isec14,ifild)
      call revise
      call dr2fre(qq(last),qq(ii1),qq(ii2),qq(ii0),qq(ii3),qq(ii4),
     +     c,de,nat3,iwr,last)
      ti = tim0
      call timit(3)
      return
 6010 format (/1x,35('*')/1x,'2-electron 2nd-derivative integrals'/1x,
     +        35('*')//
     +        ' commence 2-electron derivative integral evaluation at ',
     +        f8.2,' seconds',a10,' wall')
 6020 format (/1x,'core available  ',i7,' words'/1x,'core required   ',
     +        i7,' words'//)
 6030 format (/20x,'*****************************************'/
     +        20x,'* total cartesian 2nd derivative matrix *'/
     +        20x,'*****************************************')
      end
c
c-----------------------------------------------------------------------
c
      subroutine dr2dft(q,iq,isec46)
      implicit none
c------------------------------------------------------------------
c     include the 2nd derivative contributions for DFT calculations
c
c     Essentially we only call the appropriate routine from the DFT
c     module. But we need to mess around a bit with the GAMESS-UK
c     data structures to make sure that everything gets added on to
c     the right things to make it all work in the end...
c------------------------------------------------------------------
c
c     Parameters
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
c     In variables
c
INCLUDE(common/restri)
INCLUDE(common/restar)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/atmblk)
INCLUDE(common/drive_dft)
c
      integer isec46
c
c     Work space
c
      REAL q(*)
      integer iq(*)
c
c     Local variables
c
      integer nat3
      integer nlen, l2
      integer ierror, i, inull
      integer ida, idb, id0, id1
      integer isecd1
      character *7 fnm
      character *6 snm
      data fnm/'sec2e.m'/
      data snm/'dr2dft'/
c
c     Functions
c
      integer igmem_alloc_inf
      integer igmem_null
INCLUDE(common/ccpdft.hf77)
c
c     Code
c
      if (.not.CD_active()) return
      inull = igmem_null()
c
      call gmem_summarize(fnm,snm,'before_dft_hess',IGMEM_DEBUG)
c
      nat3 = 3*nat
      nlen = nat3*nat3
      l2   = num*(num+1)/2
      id0  = igmem_alloc_inf(nlen,fnm,snm,'dft_hess',IGMEM_DEBUG)
      ida  = igmem_alloc_inf(num*num,fnm,snm,'alpha-vectors',
     &                       IGMEM_NORMAL)
c
      call vclr(q(id0),1,nlen)
      if (ks_hes_bas.eq.KS_HES_AO) then
         call secget(isect(7),7,isecd1)
         call rdedx(q(ida),lds(isect(7)),isecd1,ifild)
         ierror = CD_hess_ao(iq,q,q(ida),q(inull),q(id0),.false.,iwr)
      else if (ks_hes_bas.eq.KS_HES_MO) then
         call secget(isect(8),8,isecd1)
         isecd1 = isecd1 + mvadd
         call rdedx(q(ida),num*ncoorb,isecd1,ifild)
         ierror = CD_hess_mo(iq,q,ncoorb,nocc,0,q(ida),q(inull),
     &                       q(id0),.false.,iwr)
      else
         write(iwr,*)'ks_hes_bas = ',ks_hes_bas
         call caserr('dr2dft: invalid value of ks_hes_bas!')
      endif
      call gmem_free_inf(ida,fnm,snm,'alpha-vectors')
c
c        ierror = CD_hess_ao(iq,q,q(ida),q(inull),q(id0),.false.,iwr)
c
c     Add the DFT contribution to the stuff stored on isec46
c
      id1 = igmem_alloc_inf(nlen,fnm,snm,'tmp_hess',IGMEM_DEBUG)
      call rdedx(q(id1),nlen,isec46,ifild)
      do i = 0, nlen-1
         q(id1+i) = q(id1+i) + q(id0+i)
      enddo
      call wrt3(q(id1),nlen,isec46,ifild)
      call gmem_free_inf(id1,fnm,snm,'tmp_hess')
c
      call gmem_free_inf(id0,fnm,snm,'dft_hess')
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dr2ndf(abdens)
c----------------------------------------------
c     insert normalisation of d and f-functions
c----------------------------------------------
      implicit REAL  (a-h,o-z)
      dimension abdens(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/ddmisc)
INCLUDE(common/ddshln)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225),dd(144)
      data one,root3,root5/1.0d0,1.7320508076d0,2.2360679775d0/
      n = 0
      jmax = maxj
      dum1 = one
      do 30 i = mini , maxi
         if (i.eq.8) dum1 = root3
         if (i.eq.14) dum1 = root5
         if (i.eq.20) dum1 = dum1*root3
         dum2 = dum1
         if (iandj) jmax = i
         do 20 j = minj , jmax
            if (j.eq.8) dum2 = dum2*root3
            if (j.eq.14) dum2 = dum2*root5
            if (j.eq.20) dum2 = dum2*root3
            n = n + 1
            dij(n) = dum2
 20      continue
 30   continue
      n = 0
      dum1 = one
      lmax = maxl
      do 50 k = mink , maxk
         if (k.eq.8) dum1 = root3
         if (k.eq.14) dum1 = root5
         if (k.eq.20) dum1 = dum1*root3
         dum2 = dum1
         if (kandl) lmax = k
         do 40 l = minl , lmax
            if (l.eq.8) dum2 = dum2*root3
            if (l.eq.14) dum2 = dum2*root5
            if (l.eq.20) dum2 = dum2*root3
            n = n + 1
            dkl(n) = dum2
 40      continue
 50   continue
      do 70 i = 1 , ijd2
         d1 = dij(i)
         n1 = ijgt(i)
         maxik = ik(i)
         do 60 k = 1 , maxik
            n = n1 + klgt(k)
            abdens(n) = abdens(n)*d1*dkl(k)
 60      continue
 70   continue
      return
      end
      subroutine dr2prm
c---------------------------------------------
c     information regarding sets of primitives
c---------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/ddmisc) 
INCLUDE(common/ddshln)
INCLUDE(common/ddshli)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     + aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     + a(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),dij(4*mxp2),
     + ijden(225)
      data one/1.0d0/
      jmax = maxj
      n = 0
      nn = 0
      do 70 i = mini , maxi
         go to (20,20,30,30,20,30,30,30,30,30,20,30,30,30,30,30,30,30,
     +          30,30) , i
 20      nm = nn
 30      nn = nm
         if (iandj) jmax = i
         do 60 j = minj , jmax
            go to (40,40,50,50,40,50,50,50,50,50,40,50,50,50,50,50,50,
     +             50,50,50) , j
 40         nn = nn + 1
 50         n = n + 1
            ijden(n) = nn
 60      continue
 70   continue
c     ----- i primitive
      nij = 0
      jbmax = ngb
      do 230 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*xi
         ayi = ai*yi
         azi = ai*zi
         csi = csa(ia)
         cpi = cpa(ia)
         cdi = cda(ia)
         cfi = cfa(ia)
c     ----- j primitive
         if (iandj) jbmax = ia
         do 220 jb = 1 , jbmax
            aj = bg(jb)
            aa = ai + aj
            ainv = one/aa
            dum = aj*arri*ainv
            csj = csb(jb)*ainv
            cpj = cpb(jb)*ainv
            cdj = cdb(jb)*ainv
            cfj = cfb(jb)*ainv
            nm = (nij+nij) + (nij+nij)
            nn = nm
            nij = nij + 1
            r(nij) = dum
            a(nij) = aa
            x1(nij) = (axi+aj*xj)*ainv
            y1(nij) = (ayi+aj*yj)*ainv
            z1(nij) = (azi+aj*zj)*ainv
c     ----- density factor
            do 190 i = mini , maxi
               if (iandj) jmax = i
               go to (80,90,190,190,100,190,190,190,190,190,110,190,190,
     +                190,190,190,190,190,190,190) , i
 80            dum1 = csi
               go to 120
 90            dum1 = cpi
               go to 120
 100           dum1 = cdi
               go to 120
 110           dum1 = cfi
 120           do 180 j = minj , jmax
                  go to (130,140,180,180,150,180,180,180,180,180,160,
     +                   180,180,180,180,180,180,180,180,180) , j
 130              dum2 = dum1*csj
                  go to 170
 140              dum2 = dum1*cpj
                  go to 170
 150              dum2 = dum1*cdj
                  go to 170
 160              dum2 = dum1*cfj
 170              nn = nn + 1
                  dij(nn) = dum2
 180           continue
 190        continue
            if (.not.iandj) go to 220
            if (ia.eq.jb) go to 220
            go to (210,200,210,210) , lit
 200        if (mini.ne.2) then
               dij(nm+2) = dij(nm+2) + csi*cpj
               dij(nm+3) = dij(nm+3) + dij(nm+3)
            end if
 210        dij(nm+1) = dij(nm+1) + dij(nm+1)
 220     continue
 230  continue
      if (nij.eq.0) return
      rsmall = r(1)
      do 240 n = 1 , nij
         exij(n) = dexp(-r(n))
 240  continue
      do 250 n = 1 , nij
         if (rsmall.gt.r(n)) rsmall = r(n)
 250  continue
      if (rsmall.ge.tol(1)) nij = 0
      return
      end
      subroutine dr2mas(iv,rm)
c--------------------------------------------------------------------
c     mass weighting of force constant routine
c     basically the hondo routine
c-weights
c     modified to deal with multiple isotopic substitution vectors
c     arg iv specifies which vector we want
c-------------------------------------------------------------------
c     ----- construct the -g- matrix.
c
      implicit REAL  (a-h,o-z)
      dimension rm(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
c
      data one /1.0d0/
c
      k = 0
      do 40 iat = 1 , nat
         dum = one/dsqrt(amass_get(iv,iat))
         do 30 j = 1 , 3
            k = k + 1
            rm(k) = dum
 30      continue
 40   continue
      return
      end
      subroutine dr2set(ni,nj,nk,nl,lit,ljt,lkt,llt,skip)
c--------------------------------
c     set which of shells to skip
c--------------------------------
      implicit REAL  (a-h,o-z)
      logical skip(10)
c
      do 20 i = 1 , 10
         skip(i) = .not.skip(i)
 20   continue
c     (1)=i'i'    (2)=i'j'    (3)=i'k'    (4)=i'l'
c     (5)=j'j'    (6)=j'k'    (7)=j'l'
c     (8)=k'k'    (9)=k'l'
c     (10)=l'l'
      ni = lit - 1
      nj = ljt - 1
      nk = lkt - 1
      nl = llt - 1
      if (skip(2) .or. skip(3) .or. skip(4)) ni = lit
      if (skip(2) .or. skip(6) .or. skip(7)) nj = ljt
      if (skip(3) .or. skip(6) .or. skip(9)) nk = lkt
      if (skip(4) .or. skip(7) .or. skip(9)) nl = llt
      if (skip(1)) ni = lit + 1
      if (skip(5)) nj = ljt + 1
      if (skip(8)) nk = lkt + 1
      if (skip(10)) nl = llt + 1
      do 30 i = 1 , 10
         skip(i) = .not.skip(i)
 30   continue
      return
      end
      subroutine dr2shl(nelec,ish,jsh,ksh,lsh)
c------------------------------------------
c     information regarding pairs of shells
c------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/ddmisc)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/root)
INCLUDE(common/ddshln)
INCLUDE(common/ddshli)
INCLUDE(common/incrdd)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),
     + aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     + aaa(9*mxp2),ijden(225),
     + ik(225),ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225),ddtemp(144)
c
      if (nelec.eq.2) then
         kandl = ksh.eq.lsh
         same = ish.eq.ksh .and. jsh.eq.lsh
         k = katom(ksh)
         xk = c(1,k)
         yk = c(2,k)
         zk = c(3,k)
         k1 = kstart(ksh)
         k2 = k1 + kng(ksh) - 1
         lkt = ktype(ksh)
         mink = kmin(ksh)
         maxk = kmax(ksh)
         lock = kloc(ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex(k)
            csc(ngc) = cs(k)
            cpc(ngc) = cp(k)
            cdc(ngc) = cd(k)
            cfc(ngc) = cf(k)
 20      continue
         l = katom(lsh)
         xl = c(1,l)
         yl = c(2,l)
         zl = c(3,l)
         l1 = kstart(lsh)
         l2 = l1 + kng(lsh) - 1
         llt = ktype(lsh)
         minl = kmin(lsh)
         maxl = kmax(lsh)
         locl = kloc(lsh) - minl
         ngd = 0
         do 30 l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex(l)
            csd(ngd) = cs(l)
            cpd(ngd) = cp(l)
            cdd(ngd) = cd(l)
            cfd(ngd) = cf(l)
 30      continue
         nroots = (lit+ljt+lkt+llt)/2
         rrk = ((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
         inc2 = 1
         inc3 = inc2*(maxl-minl+1)
         inc4 = inc3*(maxk-mink+1)
         inc5 = inc4*(maxj-minj+1)
         lendd = inc5*(maxi-mini+1)
         if (mod(lendd,8).eq.0) lendd = lendd + 1
         ijd2 = 0
         jmax = maxj
         do 50 i = mini , maxi
            if (iandj) jmax = i
            ittt = inc5*(i-mini) + 1
            do 40 j = minj , jmax
               ijd2 = ijd2 + 1
               ijgt(ijd2) = ittt
               ittt = ittt + inc4
 40         continue
 50      continue
         kld2 = 0
         lmax = maxl
         do 70 k = mink , maxk
            if (kandl) lmax = k
            ittt = inc3*(k-mink)
            do 60 l = minl , lmax
               kld2 = kld2 + 1
               klgt(kld2) = ittt
               ittt = ittt + inc2
 60         continue
 70      continue
         do 80 i = 1 , ijd2
            ik(i) = kld2
            if (same) ik(i) = i
 80      continue
         ijkld2 = ijd2*kld2
         if (same) ijkld2 = ijd2*(ijd2+1)/2
         ixi = lendd + 1
         iyi = ixi + lendd
         izi = iyi + lendd
         isub = iabd + lendd
         return
      else
         iandj = ish.eq.jsh
         i = katom(ish)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ish)
         i2 = i1 + kng(ish) - 1
         lit = ktype(ish)
         mini = kmin(ish)
         maxi = kmax(ish)
         loci = kloc(ish) - mini
         nga = 0
         do 90 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex(i)
            csa(nga) = cs(i)
            cpa(nga) = cp(i)
            cda(nga) = cd(i)
            cfa(nga) = cf(i)
 90      continue
         j = katom(jsh)
         xj = c(1,j)
         yj = c(2,j)
         zj = c(3,j)
         j1 = kstart(jsh)
         j2 = j1 + kng(jsh) - 1
         ljt = ktype(jsh)
         minj = kmin(jsh)
         maxj = kmax(jsh)
         locj = kloc(jsh) - minj
         ngb = 0
         do 100 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex(j)
            csb(ngb) = cs(j)
            cpb(ngb) = cp(j)
            cdb(ngb) = cd(j)
            cfb(ngb) = cf(j)
 100     continue
         rri = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
         return
      end if
      end
      subroutine dr2skp(ii,jj,kk,ll,iw,outdb)
c-------------------------------------------------------------------
c     determine which shells omitted due to translational
c     invariance
c-------------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      logical outdb
c
INCLUDE(common/sizes)
INCLUDE(common/ddmisc)
INCLUDE(common/nshel)
c
      logical test
      logical inej,inek,inel,jnek,jnel,knel
c
      data test/.false./
      do 20 i = 1 , 10
         skip(i) = .true.
         diag(i) = .true.
         noff(i) = 0
 20   continue
      ncen = 0
      npass = 0
      do 30 i = 1 , 4
         natom(i) = 0
 30   continue
      lit = ktype(ii)
      ljt = ktype(jj)
      lkt = ktype(kk)
      llt = ktype(ll)
      iat = katom(ii)
      jat = katom(jj)
      kat = katom(kk)
      lat = katom(ll)
      inej = iat.ne.jat
      inek = iat.ne.kat
      inel = iat.ne.lat
      jnek = jat.ne.kat
      jnel = jat.ne.lat
      knel = kat.ne.lat
      if (inej) then
         if (inek) then
            if (jnek) then
               if (.not.(jnel)) go to 50
               if (.not.(knel)) go to 60
c     iat # jat # kat # lat  -- omit one centre
               if (test) then
c     alternative for 4-center integrals
                  skip(1) = .false.
                  skip(2) = .false.
                  skip(3) = .false.
                  skip(5) = .false.
                  skip(6) = .false.
                  skip(8) = .false.
                  diag(2) = .false.
                  diag(3) = .false.
                  diag(6) = .false.
                  natom(1) = iat
                  natom(2) = jat
                  natom(3) = kat
                  natom(4) = lat
                  ncen = 4
                  npass = 6
                  noff(1) = 0
                  noff(2) = 3
                  noff(3) = 6
                  noff(5) = 39
                  noff(6) = 42
                  noff(8) = 78
               else
                  skip(2) = .false.
                  skip(3) = .false.
                  diag(3) = .false.
                  skip(4) = .false.
                  diag(4) = .false.
                  skip(6) = .false.
                  skip(7) = .false.
                  diag(7) = .false.
                  skip(9) = .false.
                  ncen = 4
                  natom(1) = iat
                  natom(2) = jat
                  natom(3) = kat
                  natom(4) = lat
                  noff(2) = 3
                  noff(3) = 6
                  noff(4) = 9
                  noff(6) = 42
                  noff(7) = 45
                  noff(9) = 81
                  npass = 6
               end if
               go to 70
            else
               if (jnel) then
c     ----- jat = kat # iat # lat -----
                  skip(1) = .false.
                  skip(4) = .false.
                  diag(4) = .false.
                  skip(10) = .false.
                  natom(1) = iat
                  natom(2) = lat
                  natom(3) = jat
                  ncen = 3
                  noff(4) = 3
                  noff(10) = 39
                  npass = 3
               else
c     jat=kat=lat    (i'j/kl)
                  skip(1) = .false.
                  natom(1) = iat
                  natom(2) = jat
                  ncen = 2
                  npass = 1
               end if
               go to 70
            end if
         else if (inel) then
            if (jnel) go to 40
c     iat=kat   jat=lat
            n1 = (lit+1)*(lkt+1)*ljt*llt
            n2 = lit*lkt*(ljt+1)*(llt+1)
            if (n1.ge.n2) go to 40
            go to 50
         else
c     iat=kat=lat   derivative (ij''/kl)
            skip(5) = .false.
            natom(1) = jat
            natom(2) = iat
            ncen = 2
            npass = 1
            go to 70
         end if
      else if (inek) then
         if (.not.(knel)) then
c     iat=jat  ,  kat=lat   differentiate one pair
            n1 = lit*ljt*(lkt+1)*(llt+1)
            n2 = (lit+1)*(ljt+1)*lkt*llt
            if (n2.lt.n1) go to 60
         end if
      else
         if (inel) then
c     iat=jat=kat   derivative (ij/kl'')
            skip(10) = .false.
            natom(1) = lat
            natom(2) = iat
            ncen = 2
            npass = 1
         end if
         go to 70
      end if
c     iat=jat   derivatives (ij/k''l) (ij/k'l') (ij/kl'')
      skip(8) = .false.
      skip(9) = .false.
      diag(9) = .false.
      skip(10) = .false.
      noff(9) = 3
      noff(10) = 39
      natom(1) = kat
      natom(2) = lat
      natom(3) = iat
      ncen = 3
      npass = 3
      go to 70
c     iat=kat  do (ij''/kl) (ij'/kl') and (ij/kl'')
 40   skip(5) = .false.
      skip(7) = .false.
      diag(7) = .false.
      skip(10) = .false.
      natom(1) = jat
      natom(2) = lat
      natom(3) = iat
      ncen = 3
      noff(7) = 3
      noff(10) = 39
      npass = 3
      go to 70
c      jat=lat   (i''j/kl) (i'j/k'l) and (ij/k''l)
 50   skip(1) = .false.
      skip(3) = .false.
      diag(3) = .false.
      skip(8) = .false.
      natom(1) = iat
      natom(2) = kat
      natom(3) = jat
      ncen = 3
      noff(3) = 3
      noff(8) = 39
      npass = 3
      go to 70
c     kat=lat   derivatives (i'j/kl) and (ij'/kl)
 60   skip(1) = .false.
      skip(2) = .false.
      diag(2) = .false.
      skip(5) = .false.
      natom(1) = iat
      natom(2) = jat
      natom(3) = kat
      ncen = 3
      noff(2) = 3
      noff(5) = 39
      npass = 3
 70   if (.not.outdb) return
      write (iw,6010) ii , jj , kk , ll , natom , noff , skip , diag
      return
 6010 format (5x,4i5,5x,4i5,5x,10i5/5x,10l4,5x,10l4)
      end
      subroutine dr2ini(da,db,dd,qq)
c----------------------------
c     initialisation routine
c----------------------------
      implicit REAL  (a-h,o-z)
      dimension da(*),db(*),dd(*),qq(*)
c
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
INCLUDE(common/ddmisc)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/ghfblk)
INCLUDE(common/infoa)
INCLUDE(common/prnprn)
c
      character *8 grhf,closed,open
      data grhf/'grhf'/
      data open,one,ten,e /'open',1.0d0,1.0d1,2.30258d0/
      data closed/'closed'/
      data half/0.5d0/
c
      if (scftyp.eq.grhf) then
c
c     general case
c
         m = 0
         call secget(isect(8),m,iblok)
         call rdedx(qq,num*ncoorb,iblok+mvadd,ifild)
         m = 0
         do 50 is = 1 , njk
            call vclr(da(m+1),1,nx)
            nsi = nbshel(is)
            ils = ilfshl(is)
            do 40 ni = 1 , nsi
               nc = iactiv(ils+ni)
               ncol = (nc-1)*num
               ij = 0
               do 30 i = 1 , num
                  qqin = qq(i+ncol)
                  do 20 j = 1 , i
                     ij = ij + 1
                     da(ij+m) = da(ij+m) + qqin*qq(j+ncol)
 20               continue
 30            continue
 40         continue
            m = m + nx
 50      continue
c
c
      else
c
c     ordinary scf density matrices
c
         call secget(isect(7),7,isecd1)
         call rdedx(da,lds(isect(7)),isecd1,ifild)
         if (scftyp.ne.closed) then
            call secget(isect(10),10,isecd2)
            call rdedx(db,lds(isect(10)),isecd2,ifild)
            if (oss) call secget(isect(41),41,isecd3)
            if (oss) call rdedx(dd,lds(isect(41)),isecd3,ifild)
            do 60 i = 1 , nx
               duma = da(i)
               dumb = db(i)
               da(i) = duma + dumb
               if (scftyp.eq.open) db(i) = duma - dumb
               if (oss) then
                  dumd = dd(i)
                  db(i) = half*(dumb+dumd)
                  dd(i) = half*(dumb-dumd)
               end if
 60         continue
         end if
      end if
c
c
      out = odebug(7)
      icutd = icut - 3
      itold = itol - 1
      itold = max(itol,16)
      cutoff = one/(ten**icutd)
      tol(1) = e*(itold+2)
      tol(2) = e*itold
      tol(3) = one/ten**(itold+2)
      tol(4) = one/ten**itold
      norm = normf.ne.1 .or. normp.ne.1
      return
c
      end
      subroutine dr2pri(x,y,z,xd,yd,zd,b,m1,m2,m3,m4,i1,i2,k1,k2,
     +    ncdim)
c------------------------------------------------------------------
c     calculate derivatives of primitive integrals from
c     undifferentiated set
c-----------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),
     +         xd(ncdim,*),yd(ncdim,*),zd(ncdim,*)
      dimension b(*)
c
INCLUDE(common/ddshln)
c
      parameter (ncmax=65)
      dimension a(ncmax)
c
      do 20 i = 1 , ncontr
         a(i) = b(i) + b(i)
 20   continue
      n1 = 1
c
c
      do 120 i = 1 , m1
         n2 = n1
         do 110 j = 1 , m2
            n3 = n2
            do 100 k = 1 , m3
               n4 = n3
               do 90 l = 1 , m4
                  go to (30,50,70,70,70,70,70) , l
c
c
 30               do 40 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2)
 40               continue
                  n4 = n4 + k2
                  go to 90
c
c
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 50               do 60 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) - x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) - y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) - z(nr,n4-k2)
 60               continue
                  n4 = n4 + k2
                  go to 90
c
c
c
 70               fac = -dfloat(l-1)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 80 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) + fac*x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) + fac*y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) + fac*z(nr,n4-k2)
 80               continue
c
                  n4 = n4 + k2
 90            continue
               n3 = n3 + k1
 100        continue
            n2 = n2 + i2
 110     continue
         n1 = n1 + i1
 120  continue
      return
      end
      subroutine dr2trn(dd,ncen,natom,iw,out)
c-------------------------------------------------------
c     find terms missing due to translational invariance
c-------------------------------------------------------
      implicit REAL  (a-h,o-z)
      logical out
      dimension dd(12,12),natom(4)
c
      logical test
      data test/.false./
c
      go to (250,260,320,20) , ncen
 20   if (test) then
c     alternative way for 4-centers
         do 40 i = 1 , 9
            do 30 j = 1 , i
               dd(j,i) = dd(i,j)
 30         continue
 40      continue
         do 60 i = 10 , 12
            do 50 j = 1 , 9
               dd(i,j) = -(dd(i-9,j)+dd(i-6,j)+dd(i-3,j))
               dd(j,i) = dd(i,j)
 50         continue
 60      continue
         do 80 i = 10 , 12
            do 70 j = 10 , 12
               dd(i,j) = -(dd(i-9,j)+dd(i-6,j)+dd(i-3,j))
 70         continue
 80      continue
         if (out) call prsqm(dd,12,12,12,iw)
         return
      else
         do 100 i = 4 , 12
            i3 = i - 3
            do 90 j = 1 , i3
               dd(j,i) = dd(i,j)
 90         continue
 100     continue
         do 120 i = 1 , 3
            do 110 j = 1 , i
               dd(i,j) = -(dd(i+3,j)+dd(i+6,j)+dd(i+9,j))
               dd(j,i) = dd(i,j)
 110        continue
 120     continue
         do 140 i = 4 , 5
            i2 = i - 2
            do 130 j = i2 , 3
               dd(i,j) = -(dd(i-3,j)+dd(i+3,j)+dd(i+6,j))
               dd(j,i) = dd(i,j)
 130        continue
 140     continue
         do 160 i = 4 , 6
            do 150 j = 4 , i
               dd(i,j) = -(dd(i-3,j)+dd(i+3,j)+dd(i+6,j))
               dd(j,i) = dd(i,j)
 150        continue
 160     continue
         do 180 i = 7 , 8
            i2 = i - 2
            do 170 j = i2 , 6
               dd(i,j) = -(dd(i-6,j)+dd(i-3,j)+dd(i+3,j))
               dd(j,i) = dd(i,j)
 170        continue
 180     continue
         do 200 i = 7 , 9
            do 190 j = 7 , i
               dd(i,j) = -(dd(i-6,j)+dd(i-3,j)+dd(i+3,j))
               dd(j,i) = dd(i,j)
 190        continue
 200     continue
         do 220 i = 10 , 11
            i2 = i - 2
            do 210 j = i2 , 9
               dd(i,j) = -(dd(i-9,j)+dd(i-6,j)+dd(i-3,j))
               dd(j,i) = dd(i,j)
 210        continue
 220     continue
         do 240 i = 10 , 12
            do 230 j = 10 , i
               dd(i,j) = -(dd(i-9,j)+dd(i-6,j)+dd(i-3,j))
               dd(j,i) = dd(i,j)
 230        continue
 240     continue
         if (out) call prsqm(dd,12,12,12,iw)
      end if
 250  return
c     two centres
 260  do 280 i = 1 , 3
         do 270 j = 1 , i
            dd(j,i) = dd(i,j)
 270     continue
 280  continue
      do 300 i = 1 , 3
         do 290 j = 1 , 3
            dd(i+3,j+3) = dd(i,j)
            dd(i+3,j) = -dd(i,j)
            dd(i,j+3) = -dd(i,j)
 290     continue
 300  continue
 310  if (out) call prsqm(dd,6,6,12,iw)
      return
c     three centres
 320  do 340 i = 1 , 6
         do 330 j = 1 , i
            dd(j,i) = dd(i,j)
 330     continue
 340  continue
      do 360 i = 7 , 9
         do 350 j = 1 , 6
            dd(i,j) = -(dd(i-3,j)+dd(i-6,j))
            dd(j,i) = dd(i,j)
 350     continue
 360  continue
      do 380 i = 7 , 9
         do 370 j = 7 , 9
            dd(i,j) = -(dd(i-6,j)+dd(i-3,j))
 370     continue
 380  continue
      if (natom(1).eq.natom(2)) then
         do 400 i = 1 , 3
            do 390 j = 1 , 9
               dd(i,j) = dd(i,j) + dd(i+3,j)
               dd(i+3,j) = dd(i+6,j)
               dd(i+6,j) = 0.0d0
 390        continue
 400     continue
         do 420 i = 1 , 6
            do 410 j = 1 , 3
               dd(i,j) = dd(i,j) + dd(i,j+3)
               dd(i,j+3) = dd(i,j+6)
               dd(i,j+6) = 0.0d0
 410        continue
 420     continue
         natom(2) = natom(3)
         natom(3) = natom(4)
         natom(4) = 0
         ncen = ncen - 1
         go to 310
      else
         if (out) call prsqm(dd,9,9,12,iw)
         return
      end if
      end
      subroutine dr2xyz(x,y,z,ncdim)
c------------------------------------------------
c   primitive integrals
c------------------------------------------------
      implicit REAL  (a-h,o-z)
      dimension x(*),y(*),z(*)
c
INCLUDE(common/ddshln)
      parameter (ncmax=65)
INCLUDE(common/setdd)
      common/small/ca(ncmax),cb(ncmax)
c
      logical n0,n1,m0,m1
      dimension i(12),k(12)
      data zero,one /0.0d+00,1.0d+00/
c
c
      do 20 n = 1 , nmax + 1
         i(n) = (iorg(n)-1)*ncdim + 1
 20   continue
      do 30 n = 1 , mmax + 1
         k(n) = korg(n)*ncdim
 30   continue
      ij1 = ij1x*ncdim
      ij2 = ij2x*ncdim
      kl1 = kl1x*ncdim
      kl2 = kl2x*ncdim
      ink = 1
c
      n0 = nmax.eq.0
      n1 = nmax.le.1
      m0 = mmax.eq.0
      m1 = mmax.le.1
      if (n0) then
         i1 = i(1)
         ia = 0
         do 40 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 40      continue
         if (m0) go to 670
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
         do 50 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 50      continue
         if (m1) go to 670
c     ----- i(0,nk) -----
c
         do 60 nc = 1 , ncontr
            ca(nc) = zero
 60      continue
         i3 = i1
         i4 = i1 + k2
         do 90 nk = 2 , mmax
            do 70 nc = 1 , ncontr
               ca(nc) = ca(nc) + bp01(nc)
 70         continue
            i5 = i1 + k(nk+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 80 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
               ia = ia + ink
 80         continue
            i3 = i4
            i4 = i5
 90      continue
         if (nlmax.eq.0) go to 670
c     ----- i(0,0,nk,nl) -----
         i5 = k(mmax+1)
         min = nkmax
         go to 610
      else if (m0) then
         i1 = i(1)
         ia = 0
         do 100 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 100     continue
         if (n0) go to 600
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
         do 110 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 110     continue
         if (n1) go to 600
c     ----- i(ni,0) -----
c
         do 120 nc = 1 , ncontr
            ca(nc) = zero
 120     continue
         i3 = i1
         i4 = i2
         do 150 ni = 2 , nmax
            do 130 nc = 1 , ncontr
               ca(nc) = ca(nc) + b10(nc)
 130        continue
            i5 = i(ni+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 140 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
               ia = ia + ink
 140        continue
            i3 = i4
            i4 = i5
 150     continue
         if (njmax.eq.0) go to 600
c     ----- i(ni,nj,0,0) -----
         i5 = i(nmax+1)
         min = nimax
         go to 540
      else
c     ----- i(0,0) -----
         i1 = i(1)
         ia = 0
c
         do 160 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 160     continue
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
         do 170 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 170     continue
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
         do 180 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 180     continue
c     ----- i(1,1) -----
         i3 = i2 + k2
         ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 190 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)*x(i2+ia) + b00(nc)
            y(i3+ia) = ycp00(nc)*y(i2+ia) + b00(nc)
            z(i3+ia) = zcp00(nc)*z(i2+ia) + b00(nc)*f00(nc)
            ia = ia + ink
 190     continue
         if (.not.(n1)) then
            do 200 nc = 1 , ncontr
               cb(nc) = b00(nc)
               ca(nc) = zero
 200        continue
            i3 = i1
            i4 = i2
            do 250 n = 2 , nmax
c     ----- i(n,0) -----
               i5 = i(n+1)
               do 210 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b10(nc)
 210           continue
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 220 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
                  ia = ia + ink
 220           continue
               do 230 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 230           continue
c     ----- i(n,1) -----
               i3 = i5 + k2
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 240 nc = 1 , ncontr
                  x(i3+ia) = xcp00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = ycp00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zcp00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 240           continue
               i3 = i4
               i4 = i5
 250        continue
         end if
         if (.not.(m1)) then
            do 260 nc = 1 , ncontr
               ca(nc) = zero
               cb(nc) = b00(nc)
 260        continue
            i3 = i1
            i4 = i1 + k2
            do 310 m = 2 , mmax
               do 270 nc = 1 , ncontr
                  ca(nc) = ca(nc) + bp01(nc)
 270           continue
c     ----- i(0,m) -----
               i5 = i1 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 280 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
                  ia = ia + ink
 280           continue
               do 290 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 290           continue
c     ----- i(1,m) -----
               i3 = i2 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 300 nc = 1 , ncontr
                  x(i3+ia) = xc00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = yc00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zc00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 300           continue
               i3 = i4
               i4 = i5
 310        continue
         end if
         if (.not.(n1 .or. m1)) then
c     ----- i(n,m) -----
c
            do 320 nc = 1 , ncontr
               ca(nc) = b00(nc)
 320        continue
            k3 = k2
            do 370 m = 2 , mmax
               k4 = k(m+1)
c
               do 330 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b00(nc)
                  cb(nc) = b10(nc)
 330           continue
               i3 = i1
               i4 = i2
               do 360 n = 2 , nmax
                  i5 = i(n+1)
                  ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 340 nc = 1 , ncontr
                     x(i5+k4+ia) = cb(nc)*x(i3+k4+ia) + xc00(nc)
     +                             *x(i4+k4+ia) + ca(nc)*x(i4+k3+ia)
                     y(i5+k4+ia) = cb(nc)*y(i3+k4+ia) + yc00(nc)
     +                             *y(i4+k4+ia) + ca(nc)*y(i4+k3+ia)
                     z(i5+k4+ia) = cb(nc)*z(i3+k4+ia) + zc00(nc)
     +                             *z(i4+k4+ia) + ca(nc)*z(i4+k3+ia)
                     ia = ia + ink
 340              continue
                  do 350 nc = 1 , ncontr
                     cb(nc) = cb(nc) + b10(nc)
 350              continue
                  i3 = i4
                  i4 = i5
 360           continue
               k3 = k4
 370        continue
         end if
         if (njmax.eq.0) go to 450
c     ----- i(ni,nj,m) -----
         m = 0
         i5 = i(nmax+1)
      end if
 380  min = nimax
      km = k(m+1)
 390  n = nmax
      i3 = i5 + km
 400  i4 = i(n) + km
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 410 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 410  continue
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 400
      min = min + 1
      if (min.lt.nmax) go to 390
      if (nimax.ne.0) then
         i3 = ij2 + km + i1
         do 440 nj = 1 , njmax
            i4 = i3
            do 430 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 420 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 420           continue
               i4 = i4 + ij1
 430        continue
            i3 = i3 + ij2
 440     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 380
 450  if (nlmax.eq.0) go to 530
c
c     ----- i(ni,nj,nk,nl) -----
c
      i5 = k(mmax+1)
      ia = i1
      ni = 0
 460  nj = 0
      ib = ia
      min = nkmax
 470  m = mmax
      i3 = ib + i5
 480  i4 = ib + k(m)
      ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 490 nc = 1 , ncontr
         x(i3+ic) = x(i3+ic) + dxkl*x(i4+ic)
         y(i3+ic) = y(i3+ic) + dykl*y(i4+ic)
         z(i3+ic) = z(i3+ic) + dzkl*z(i4+ic)
         ic = ic + ink
 490  continue
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 480
      min = min + 1
      if (min.lt.mmax) go to 470
      if (nkmax.ne.0) then
         i3 = ib + kl2
         do 520 nl = 1 , nlmax
            i4 = i3
            do 510 nk = 1 , nkmax
               ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 500 nc = 1 , ncontr
                  x(i4+ic) = x(i4+ic+kl1-kl2) + dxkl*x(i4+ic-kl2)
                  y(i4+ic) = y(i4+ic+kl1-kl2) + dykl*y(i4+ic-kl2)
                  z(i4+ic) = z(i4+ic+kl1-kl2) + dzkl*z(i4+ic-kl2)
                  ic = ic + ink
 500           continue
               i4 = i4 + kl1
 510        continue
            i3 = i3 + kl2
 520     continue
      end if
      nj = nj + 1
      ib = ib + ij2
      if (nj.le.njmax) then
         min = nkmax
         go to 470
      else
         ni = ni + 1
         ia = ia + ij1
         if (ni.le.nimax) go to 460
      end if
 530  return
 540  ni = nmax
      i3 = i5
 550  i4 = i(ni)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 560 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 560  continue
      i3 = i4
      ni = ni - 1
      if (ni.gt.min) go to 550
      min = min + 1
      if (min.lt.nmax) go to 540
      if (nimax.ne.0) then
         i3 = ij2 + i1
         do 590 nj = 1 , njmax
            i4 = i3
            do 580 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 570 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 570           continue
               i4 = i4 + ij1
 580        continue
            i3 = i3 + ij2
 590     continue
      end if
 600  return
 610  nk = mmax
      i3 = i1 + i5
 620  i4 = i1 + k(nk)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 630 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxkl*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dykl*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzkl*z(i4+ia)
         ia = ia + ink
 630  continue
      i3 = i4
      nk = nk - 1
      if (nk.gt.min) go to 620
      min = min + 1
      if (min.lt.mmax) go to 610
      if (nkmax.ne.0) then
         i3 = i1 + kl2
         do 660 nl = 1 , nlmax
            i4 = i3
            do 650 nk = 1 , nkmax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 640 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+kl1-kl2) + dxkl*x(i4+ia-kl2)
                  y(i4+ia) = y(i4+ia+kl1-kl2) + dykl*y(i4+ia-kl2)
                  z(i4+ia) = z(i4+ia+kl1-kl2) + dzkl*z(i4+ia-kl2)
                  ia = ia + ink
 640           continue
               i4 = i4 + kl1
 650        continue
            i3 = i3 + kl2
 660     continue
      end if
 670  return
      end
      subroutine dr2pdm(abdens,ii,jj,kk,ll,q4)
c-------------------------------------------------------------------
c     will read 2 - pdm for correlated cases ( so far mp2 only)
c-------------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      dimension abdens(*)
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
INCLUDE(common/cigrad)
INCLUDE(common/nshel)
INCLUDE(common/incrdd)
c
c     i/o blocks
c
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
      common/blkin/gin(510),mword,nlenx,kworx,kworxx
c
      logical ijdiff,kldiff,ikdiff
      data four,pt5/4.0d0,0.5d0/
c
      call vclr(abdens,1,lendd)
      if (kworx.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kworx)
         iword = 0
c
c     packing scheme depends on block size and number of bytes
c     per label
c
_IFN1(iv)        call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)         call upak8v(gin(num2e+1),i205)
         if (kworx.le.0) go to 30
      end if
      iword = iword + 1
_IFN1(iv)      m = (iword+iword) + (iword+iword)
_IFN1(iv)      i = labs(m-3)
_IFN1(iv)      j = labs(m-2)
_IFN1(iv)      k = labs(m-1)
_IFN1(iv)      l = labs(m)
_IF1(iv)      i = i205(iword)
_IF1(iv)      j = j205(iword)
_IF1(iv)      k = k205(iword)
_IF1(iv)      l = l205(iword)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = four
                  if (k.eq.l) d4 = d4*pt5
                  if (i.eq.k .and. j.eq.l) d4 = d4*pt5
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dr2pun(title,c,grad,fcm,nat,nat3,ip,ndim)
c------------------------------------------------------
c     punch out force constant matrix to fortran unit 7
c------------------------------------------------------
      implicit REAL  (a-h,o-z)
      character *8 title
      dimension title(10),c(3,nat),grad(3,nat),fcm(ndim,ndim)
c
      write (ip,6010) (title(i),i=1,9)
      write (ip,6020)
      write (ip,6030) (c(1,i),c(2,i),c(3,i),i=1,nat)
      write (ip,6040)
      write (ip,6030) (grad(1,i),grad(2,i),grad(3,i),i=1,nat)
      write (ip,6050)
      do 20 i = 1 , nat
         write (ip,6030) (fcm(j,i*3-2),fcm(j,i*3-1),fcm(j,i*3),j=1,nat3)
 20   continue
      return
 6010 format (1x,9a8)
 6020 format (1x,'geometry')
 6030 format (1x,3e20.12)
 6040 format (1x,'gradient')
 6050 format (1x,'cartesian second derivatives')
      end
      subroutine ver_sec2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/sec2e.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
