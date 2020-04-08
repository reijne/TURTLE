      subroutine sp0011(gout)
c        *****  special fast routine for -p- loop for 0011 *****
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
      common/miscg/mab,mcd,ngangb
c
      real*8 ap, bp, cq, dq, px, py, pz, qx, qy, qz, rpq, rpqsq
      real*8 pq1, pq2, pq3
      real*8 c11, c12, c13, c21, c22, c23, c31, c32, c33
      common /pqgeom/ ap,bp,cq,dq,px,py,pz,qx,qy,qz,rpq,rpqsq,
     +                pq1,pq2,pq3,
     +                c11,c12,c13,c21,c22,c23,c31,c32,c33
      common/ginf/ga,gb,gc,gd,sa,sb,sc,sd,pa,pb,pc,pd,gab,gcd
c
      real*8 gp, ep, dp00p, dp01p, dp10p, dp11p, app, bpp
      common /pgeom/ gp(mxprms*mxprms),ep(mxprms*mxprms),
     +               dp00p(mxprms*mxprms),dp01p(mxprms*mxprms),
     +               dp10p(mxprms*mxprms),dp11p(mxprms*mxprms),
     +               app(mxprms*mxprms),bpp(mxprms*mxprms)
c
c
      real*8 acx, acy, acz, acy2, cosg, sing
      common/qgeom/acx,acy,acz,acy2,cosg,sing
c
      common/astore/qq,theta,n
      common/inttab/
     +  a0(333),b0(333),c0(333),abc1,
     +  a1(333),b1(333),c1(333),abc2,
     +  a2(333),b2(333),c2(333),abc3,
     +  a3(333),b3(333),c3(333),abc4,
     +  a4(333),b4(333),c4(333),abc5,
     +  a5(333),b5(333),c5(333),abc6
c
      real*8 auxvar, var
      integer ifasti, ifastd
c
c     ifasti = 0 default settings in filmax
c            = 1 use GAMESS-US settings
c            = 2 tighten accuracy
c     ifastd = 0 loosen default settings (old values)
c            = 1 use revised default settings
c            = 2 tighten accuracy in g-80 routines
c
      common /auxvar/ auxvar,var(2),ifasti,ifastd
c
c
c     cmaxa/b/c/d = max of s/p contraction coeffs
c     ismlp = holds info on how to compute gamma function ?
c     ismlq = ?
c     isml  = ?
c     error1 = screening for prefactors
c     error2 = screening for prefactors
c
      real*8 cmax ,cmaxa, cmaxb, cmaxc, cmaxd, error1, error2
      integer ismlp, ismlq, isml
      common /maxc/ cmax(mxprim),
     +              cmaxa(mxprms),cmaxb(mxprms),cmaxc(mxprms),
     +              cmaxd(mxprms),error1,error2,
     +              ismlp(mxprms*mxprms),ismlq,isml
c
      real*8 ag,  csa, cpa
      real*8 bg,  csb, cpb
      real*8 cgg, csc, cpc
      real*8 dg,  csd, cpd
      integer nga, la, ngb, lb, ngc, lc, ngd, ld
      common/shllfo/nga,la,ag(mxprms) , csa(mxprms),cpa(mxprms),
     +              ngb,lb,bg(mxprms) , csb(mxprms),cpb(mxprms),
     +              ngc,lc,cgg(mxprms), csc(mxprms),cpc(mxprms),
     +              ngd,ld,dg(mxprms) , csd(mxprms),cpd(mxprms)
c
c
c     ax, ... = coords of a, b, c, d
c     rab, ... = distances
c     p/q = rotation matrix elements
c
      real*8 ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd
      real*8 rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13
      real*8 q21,q22,q23,q31,q32,q33
      common /geom/ ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,rab,rabsq,rcd,
     +          rcdsq,p11,p12,p13,p21,p22,p23,p31,p32,p33,q11,q12,q13,
     +          q21,q22,q23,q31,q32,q33
c
c
      dimension gout(*)
c
      data dzero,done/0.0d0,1.0d0/
      data sixty,tenm12/60.0d0,1.0d-12/
c
      do 940 k = 1,ngc
      gc = cgg(k)
      do 940 l = 1,ngd
      gd = dg(l)
      gcd = gc+gd
      ecd = done/gcd
      cq = gd*ecd*rcd
      dq = cq-rcd
      qqq = cq*dq*gcd
      if (qqq+sixty) 940,500,500
  500 v =  dexp(qqq)*ecd
  520 qqtest = cmaxc(k)*cmaxd(l)*v
      if (qqtest-error1) 560,560,540
  540 ismlq = 0
      go to 600
  560 if (qqtest-error2) 940,940,580
  580 ismlq = 1
  600 sc = csc(k)
      sd = csd(l)
      pc = cpc(k)
      pd = cpd(l)
      dq00 = sc*sd*v
      dq01 = sc*pd*v
      dq10 = pc*sd*v
      dq11 = pc*pd*v
      aqx = acx+sing*cq
      aqz = acz+cosg*cq
      qperp2 = aqx*aqx+acy2
      qperp = dsqrt(qperp2)
      if (qperp-tenm12) 640,640,620
  620 cosp = -aqx/qperp
      sinp = -acy/qperp
      go to 660
  640 cosp = done
      sinp = 0.0d0
 660  h0000 = 0.d0
      h0001 = 0.d0
      h0003 = 0.d0
      h0011 = 0.d0
      h0013 = 0.d0
      h0033 = 0.d0
      do 180 i = 1,ngangb
      isml = ismlq+ismlp(i)
      if (isml .ge. 2) go to 180
      auxvar = var(isml+1)
      pqab = aqz-app(i)
      pqab2 = pqab*pqab
      g = 1.d0/(ep(i)+ecd)
      p = g*(pqab2+qperp2)
      g = g*ecd
      if (p .le. auxvar) go to 140
      f0 = dp00p(i)*dsqrt(0.7853981625d0/(p*(gp(i)+gcd)))
      gtx = g/p
      f1 = 0.5d0*f0*gtx
      f2 = 1.5d0*f1*gtx
      go to 160
  140 q = dp00p(i)/dsqrt(gp(i)+gcd)
      gy = g*q
      ggy = g*gy
      qq = p*12.5d0
      n =  idint(qq)
      theta = qq- dble(n)
      theta2 = theta*(theta-1.d0)
      theta3 = theta2*(theta-2.d0)
      theta4 = theta2*(theta+1.d0)
      f0 = (a0(n+1)+theta*b0(n+1)-theta3*c0(n+1)+theta4*c0(n+2))*q
      f1 = (a1(n+1)+theta*b1(n+1)-theta3*c1(n+1)+theta4*c1(n+2))*gy
      f2 = (a2(n+1)+theta*b2(n+1)-theta3*c2(n+1)+theta4*c2(n+2))*ggy
  160 h0000 = h0000+f0
      h0001 = h0001+f1
      h0003 = h0003-f1*pqab
      h0011 = h0011+f2
      h0013 = h0013-f2*pqab
      h0033 = h0033+f2*pqab2
  180 continue
      h0022 = 0.5d0*ecd*(h0000-h0001)
      h0001 = h0001*qperp
      h0011 = h0011*qperp2+h0022
      h0013 = h0013*qperp
      h0033 = h0033+h0022
      if(sinp)120,100,120
 100  if(cosp)1000,120,920
 120  v44 = cosp*cosp
      v77 = v44
      v47 = done-v44
      v74 = v47
      v54 = cosp*sinp
      v57 = -v54
      g0011 = v44*h0011+v47*h0022
      g0012 = v54*h0011+v57*h0022
      g0022 = v74*h0011+v77*h0022
      g0013 = cosp*h0013
      g0023 = sinp*h0013
      g0033 = h0033
      g0001 = cosp*h0001
      g0002 = sinp*h0001
      g0003 = h0003
      g0000 = h0000
      go to 2000
  920 g0000 = h0000
      g0001 = h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
      go to 2000
1000  g0000 = h0000
      g0001 = -h0001
      g0002 = dzero
      g0003 = h0003
      g0011 = h0011
      g0012 = dzero
      g0013 = -h0013
      g0022 = h0022
      g0023 = dzero
      g0033 = h0033
 2000 continue
      r13 = cq*sing
      r33 = cq*cosg
      r14 = dq*sing
      r34 = dq*cosg
      g0010 = g0001
      g0020 = g0002
      g0021 = g0012
      g0030 = g0003
      g0031 = g0013
      g0032 = g0023
      if (rcdsq) 220,220,200
 200  g0010 = g0010+r13*g0000
      g0011 = g0011+r13*g0001
      g0012 = g0012+r13*g0002
      g0013 = g0013+r13*g0003
      g0030 = g0030+r33*g0000
      g0031 = g0031+r33*g0001
      g0032 = g0032+r33*g0002
      g0033 = g0033+r33*g0003
      g0001 = g0001+r14*g0000
      g0011 = g0011+r14*g0010
      g0021 = g0021+r14*g0020
      g0031 = g0031+r14*g0030
      g0003 = g0003+r34*g0000
      g0013 = g0013+r34*g0010
      g0023 = g0023+r34*g0020
      g0033 = g0033+r34*g0030
220   gout( 1) = gout( 1)+g0000*dq00
      gout( 2) = gout( 2)+g0001*dq01
      gout( 3) = gout( 3)+g0002*dq01
      gout( 4) = gout( 4)+g0003*dq01
      gout( 5) = gout( 5)+g0010*dq10
      gout( 6) = gout( 6)+g0011*dq11
      gout( 7) = gout( 7)+g0012*dq11
      gout( 8) = gout( 8)+g0013*dq11
      gout( 9) = gout( 9)+g0020*dq10
      gout( 10) = gout( 10)+g0021*dq11
      gout( 11) = gout( 11)+g0022*dq11
      gout( 12) = gout( 12)+g0023*dq11
      gout( 13) = gout( 13)+g0030*dq10
      gout( 14) = gout( 14)+g0031*dq11
      gout( 15) = gout( 15)+g0032*dq11
      gout( 16) = gout( 16)+g0033*dq11
 940  continue
      ind = 0
      do 700 l = 1,4
      ind = ind+1
      i1 = 4+ind
      i2 = 8+ind
      i3 = 12+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
      gout(i3 ) = p13*t1+p23*t2+p33*t3
  700 continue
      ind = -3
      do 720 k = 1,4
      ind = ind+4
      i1 = 1+ind
      i2 = 2+ind
      i3 = 3+ind
      t1 = gout(i1)
      t2 = gout(i2)
      t3 = gout(i3)
      gout(i3 ) = p13*t1+p23*t2+p33*t3
      gout(i1 ) = p11*t1+p21*t2+p31*t3
      gout(i2 ) = p12*t1+p22*t2+p32*t3
  720 continue
      return
      end
