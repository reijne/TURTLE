      subroutine drfamat(a,nfirst,group,np,
     1                   xpts,mpol,polar,idrfout,afact,ithole)
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c                                                                      c
c     this routine forms that part of the relay matrix to do with      c
c     the coupling between polarizabilities                            c
c     the a(-1) + t matrix elements are stored in sub-diagonal form    c
c                                                                      c
c     original literature (formulated for polarizabilities only!)      c
c     thole & van duijnen: theor.chim.acta (1980) 55 307 (form. 8,9)   c
c                                                                      c
c     thole: chem.phys.(1981) 59 341  (form.5)                         c
c                                                                      c
c                                                                      c
c     output:                                                          c
c                                                                      c
c     subdiagonal of a (vector wise)                                   c
c                                                                      c
c     routines called: clear,distab,drftpq,hatout,linv3p               c
c                                                                      c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/auxdrf)
c
INCLUDE(comdrf/grinf)
c
c-----  dummy and local variables
c
      logical group,nxtgrp
c
      dimension bpq(3,3),pq(3),b(6),p(3),q(3)
      dimension a(*),xpts(3,*),mpol(*),polar(*)
      data one,sixth /1.0d0,.16666666666666666666667d00/
c
c-----  begin
c
      indx(i,j)=ia(i)+j
      igr = 1
c
c-----  loop over polarizable points
c
      do 1000, i = 1, np
        ip=nfirst+i-1
c
c  -----  get polarizability of point -i-
c
        alfai=polar(ip)
c
c  -----  starting address in -a- for point -i-
c
        if=(i-1)*3
c
c  -----  get coordinates of point -i-
c
        ipol=mpol(ip)
        nxtgrp=ip.eq.igrpol(igr)
        if(group.and.nxtgrp) then
c
c    -----  get group polarizability
c
          do 30 kl=1,6
   30       b(kl)=grpol(kl,igr)
          if(idrfout.eq.3) call hatout(b,3,3,3,'grpol')
c
c    -----  invert group polarizability
c
          call linv3p(b,auxx,1,3,ier)
          if(idrfout.eq.3) call hatout(b,3,3,3,'grpol-1')
c
c    -----  copy -b- into -a-
c
          kl=0
          do 40 k=1,3
            do 40 l=1,k
              kl=kl+1
              indxa=indx(if+k,if+l)
              a(indxa)=b(kl)
   40     continue
          igr=igr+1
        else
c
c    -----  fill diagonal elements with inverse polarizability
c
          do 50 k=1,3
            indxa=indx(if+k,if+k)
   50       a(indxa)=one/alfai
        endif
        p(1)=xpts(1,ipol)
        p(2)=xpts(2,ipol)
        p(3)=xpts(3,ipol)
c
c  -----  loop over polarizable points (to calculate t(i;j))
c
        do 500 j=1,i
          jp=nfirst+j-1
          if(i.eq.j) goto 500
          jpol=mpol(jp)
          q(1)=xpts(1,jpol)
          q(2)=xpts(2,jpol)
          q(3)=xpts(3,jpol)
c
c    -----  staring address in -a- for point -j-
c
          jf=(j-1)*3
c
c    -----  zero sub array
c
c    -----  off diagonalblocks are (modified) dipole interactions
c
          alfaj=polar(jp)
          call distab(p,q,pq,dist)
          dmind1=one/dist
          dmin3=dmind1**3
          v=dist/((alfai*alfaj)**sixth)
          call drftpq(pq,dmind1,dmin3,ithole,afact,v,bpq)
c
c    -----  copy -bpq- into -a-
c
          do 100 k=1,3
            do 100 l=1,3
            indxa=indx(if+l,jf+k)
  100       a(indxa)=bpq(k,l)
  500   continue
c
 1000 continue
      return
      end
      subroutine drfbmat(ieps,b,xsurf,xnorm,area)
c------
c             routine calcuates coupling matrix elements for
c             boundary surface
c
c     ndima: dimension of polarizability matrix -a- (if present)
c     ndimb: dimension of linear problem
c            =   nbem if kappa.eq.0.0
c            = 2*nbem if kappa.ne.0.0
c     ndim: dimension of complete linear problem
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension b(ndim,ndim),xsurf(3,nbem),xnorm(3,nbem)
      dimension area(nbem)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
c
INCLUDE(comdrf/drfpar)
c
INCLUDE(comdrf/drfbem)
c
c-----  local arrays:
c
      dimension tij(3,3),fij(3),xij(3),tijni(3)
c
      character *8 errmsg(3)
c
      logical kapnoz
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
c
c-----  threshold for warning boundary elements are too close
c
      data thresh /1.0d-04/
c
      data errmsg /'program', 'stop in','-drfbmat'/
c
c-----  begin
c
      pi = four*atan(one)
      if (ieps .eq. 0) then
        eps = eps1
        kappa = kappa1
      else
        eps = eps2
        kappa = kappa2
      endif
c
      kapnoz = kappa .ne. zero
      kappas = kappa**2
      epsfact = one/(two*pi*(one+eps))
      expkd = one
      fkapone = one
c
c-----  loop over boundary elements
c
      do 2000, ni = 1, nbem
c
c  -----  diagonal elements equal 1
c
c  -----  1 - k(i;i) = 1
c
        b(ndima+ni,ndima+ni) = one
        if (kapnoz) then
c
c    -----  1 - n(i;i) = 1
c
          b(ndima+nbem+ni,ndima+nbem+ni) = one
c
c    -----  l(i;i) = m(i;i) = 0
c
          b(ndima+ni,ndima+nbem+ni) = zero
          b(ndima+nbem+ni,ndima+ni) = zero
        endif
c
c  -----  loop over boundary elements
c
        do 1000, nj = 1, nbem
c
c    -----  skip diagonal elements
c
          if (ni .eq. nj) goto 1000
c
c    -----  potential and field operators between elements
c
c    -----  xij = (j - i)
c
          call distab(xsurf(1,ni),xsurf(1,nj),xij,dist)
c
          if (dist .lt. thresh) then
c     3-----
            write (iwr,2011) ni, nj
 2011       format (/,' error in drfbmat: distance between elements ',
     2      i4, ' and',i4,' almost zero, check surface!!!')
            call hnderr (3,errmsg)
c     3-----
          endif
c
          dmind1 = one/dist
          dmin3 = dmind1**2*dmind1
c
c    -----  fij : field at -nj- due to charge in -ni-
c
          do 20, k = 1, 3
            fij(k) = xij(k)*dmin3
   20     continue
c
c    -----  calculate factors (only if finite ionic strength)
c
          if (kapnoz) then
            expkd = exp(-kappa*dist)
            fkapone = one + kappa*dist
          endif
c
c    -----  k-coefficient
c
c           k(i;j): potential of unit dipole density in -ni- at -nj-
c                 = (eps*(1+kappa*dist)*exp(-kappa*dist) -1) *
c                   f(i;j) dot n(i) * s(i)
c
c           negative of this element should be stored in matrix element
c           scaled by   1/(2pi(1+eps)) to satisfy the coupling equations
c
          b(ndima+nj,ndima+ni) = - epsfact*(eps*fkapone*expkd - one)*
     1               ddot(3,fij,1,xnorm(1,ni),1)*
c    1               adotb(fij,xnorm(1,ni),3)*
     2                           area(ni)
c
          if (kapnoz) then
c
c      -----  l-coefficient
c
c             l(i;j): potential of unit charge (density) in -ni- at -nj-
c                 = (1 - exp(-kappa*dist) ) v(i;j) s(i)
c
c             the negative of this matrix element should be stored in (j
c             scaled by   1/(2pi(1+eps)) to satisfy the coupling equatio
c
            b(ndima+nj,ndima+ni+nbem) =
     1            - epsfact*(one - expkd)*dmind1*area(ni)
c
c      -----  m-coefficient
c
c         m(i;j): minus field of unit dipole (density) in -ni- at -nj-
c                 contracted with normal vector at -nj-
c
c         = ((1+kap*d) exp(-kap*d) -1) * (t(i;j) dot n(i)) dot  n(j)
c           - kappa**2 exp(-kap*d) (f(i;j) dot n(i) * (j-i) dot n(j))) *
c
c         the negative of this element should be stored in (j,i)
c         scaled by eps/(2pi(1+eps)) to satisfy the coupling equations
c
c
c      -----  calculate tensor tij = t(i;j)
c
            call drftpq(xij,dmind1,dmin3,0,one,one,tij)
c
c      -----  calculate tijni = t(i;j) . n(i)
c
            call matvec(tij,xnorm(1,ni),tijni,3,.false.)
c
c      -----  calculate tijninj = (t(i;j).n(i)) . n(j)
c
            tijninj = (one - expkd*fkapone)*
     +                 ddot(3,tijni,1,xnorm(1,nj),1)
c    +                 adotb(tijni,xnorm(1,nj),3)
c
c      -----  calculate fnixnj = kap**2*exp(-kd)*f(i;j).n(i) * (j-i).n(j
c
            fnixnj = kappas*expkd*
     +               ddot(3,fij,1,xnorm(1,ni),1)*
     +               ddot(3,xij,1,xnorm(1,nj),1)
c    +               adotb(fij,xnorm(1,ni),3)*
c    +               adotb(xij,xnorm(1,nj),3)
c
c      -----  m-coefficient
c
            b(ndima+nj+nbem,ndima+ni) = - eps*epsfact*
     1                                  (tijninj - fnixnj)*area(ni)
c
c      -----  n-coefficient
c
c          n(i;j): minus field of unit charge (density) in -ni-
c                  at -nj- contacted with normal vector at -nj-
c
c           =((1+kap*dist) exp(-kap*d) / eps -1) * f(i;j) dot n(j) *s(i)
c
c          the negative of this element should be stored in (j,i)
c          scaled by eps/(2pi(1+eps)) to satisfy the coupling equations
c
            b(ndima+nj+nbem,ndima+ni+nbem) = - eps*epsfact*
     1               ((one/eps)*fkapone*expkd - one)*
     +               ddot(3,fij,1,xnorm(1,nj),1)*area(ni)
c    2               adotb(fij,xnorm(1,nj),3)*area(ni)
          endif
 1000   continue
 2000 continue
c
      return
      end
      subroutine drfdfmat(ieps,b,xsurf,xnorm,area)
c------
c      routine calculates coupling matrix elements between
c      classical polarizabilities and boundary surface elements
c
c      ndima: dimension of polarizability matrix (if present)
c      ndimb: dimension of boundary problem
c             =   nbem if kappa.eq.0.0
c             = 2*nbem if kappa.ne.0.0
c      ndim: dimension of complete linear problem
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension b(ndim,ndim),xsurf(3,nbem),xnorm(3,nbem)
      dimension area(nbem)
c
c-----  common blocks
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
c
INCLUDE(comdrf/drfbem)
c
c-----  local arrays:
c
      dimension tip(3,3),fip(3),xip(3),tipni(3),tipkip(3)
c
      logical kapnoz
c
      data zero,one,two,four/0.0d00,1.0d00,2.0d00,4.0d00/
c
c-----  begin
c
      pi = four*atan(one)
      if (ieps .eq. 0) then
        eps = eps1
        kappa = kappa1
      else
        eps = eps2
        kappa = kappa2
      endif
c
      kappas = kappa**2
      kapnoz = kappa .ne. zero
      epsfact = one/(two*pi*(one+eps))
      expkd = one
      fkapone = one
c
c-----  loop over boundary elements
c
      do 2000, ni = 1, nbem
c
c  -----  loop over polarizable points
c
        do 1000, np = 1, npol
          ipol = mpol(np)
c
c    -----  distance vector from -ipol- to surface element -ni-
c
c    -----  xip = (p-i)
c
          call distab(xsurf(1,ni),xpts(1,ipol),xip,dist)
          dmind1 = one/dist
          dmin2 = dmind1**2
          dmin3 = dmin2*dmind1
c
c    -----  fip = f(i;p): field  in -ipol- due to charge in -ni-
c
          do 30, k = 1, 3
            fip(k) = xip(k)*dmin3
   30     continue
c
c    -----  tip = t(i;p): field gradient of charge in -ni- at -ipol-
c
          call drftpq(xip,dmind1,dmin3,0,one,one,tip)
c
c    -----  calculate tipni = t(i;p) .n(i)
c           n(i) is the surface normal vector at -ni-
c
          call matvec(tip,xnorm(1,ni),tipni,3,.false.)
c
c    -----  modify factors if finite ionic strength
c
          if (kapnoz) then
            expkd = exp(-kappa*dist)
            fkapone = one + kappa*dist
          endif
c
          do 40, k = 1, 3
c
c      -----  potential of unit dipole in -ipol- at -ni-
c             = f(p;i) = -f(i;p)
c             scaled with 1/(2pi(1+eps)) to satisfy coupling equations
c
c             the negative of this matrix element should be stored in (i
c
c             plus sign is a result of the use of f(i;p), whereas f(p;i)
c             is needed
c
            b(ndima+ni,(np-1)*3+k) =   epsfact*fip(k)
c
c      -----  field of unit dipole density in -ni- at -ipol-
c             = - del(p) k(i;p)*s(i)
c             = - [(eps*(1+kd)*exp(-kd) -1)*t(i;p).n(i)
c
c               - eps*kappa**2*exp(-kd)*f(i;p).n(i)*(p-i)]*s(i)
c
c   *note* the second term is added only if kappa .ne. zero !!
c
c
c             the negative of this matrix element should be stored in (p
c
            b((np-1)*3+k,ndima+ni) = (eps*fkapone*expkd - one)*
     1                               tipni(k)*area(ni)
c
   40     continue
c
c    -----  only for poisson-boltzmann (i.e. finite ionic strength)
c
          if (kapnoz) then
c
c      -----  calculate fipni = f(i;p).n(i)
c
            fipni = ddot(3,fip,1,xnorm(1,ni),1)
c           fipni = adotb(fip,xnorm(1,ni),3)
c
c      -----  matrix elements
c
            do 50, k = 1, 3
c
c        -----  update del(p) k(i;p) *s(i)
c
              b((np-1)*3+k,ndima+ni) = b((np-1)*3+k,ndima+ni) -
     2           eps*kappas*expkd*
     3           fipni*xip(k)*area(ni)
c
c        -----  minus field of unit dipole in -ipol- at -ni-
c               contracted with normal vector at -ni-
c               = t(p;i) . n(i)
c
c               the negative of this element should be stored in (i+nbem
c               scaled by eps/(2pi(1+eps)) to satisfy the coupling equat
c
              b(ndima+nbem+ni,(np-1)*3+k)= - eps*epsfact*tipni(k)
c
c        -----  field of unit charge density in -ni- at -ipol-
c               = -del(p) l(i;p) * s(i)
c               = - [(1+kd)*exp(-kd) - 1]* f(i;p) * s(i)
c
c               the negative of this element should be stored in (p;i+nb
c
              b((np-1)*3+k,ndima+ni+nbem) =
     1               (fkapone*expkd - one)*fip(k)*area(ni)
c
   50       continue
          endif
 1000   continue
 2000 continue
c
      return
      end
      subroutine bwrit(relay,bmat,npol3,ndimb,ieps)
c------
c      writes coupling between boundary elements to dafile
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c----- dummy variables
c
      dimension relay(npol3+ndimb,npol3+ndimb),bmat(ndimb,ndimb)
c
c-----  common blocks
c
INCLUDE(comdrf/drfdaf)
c
c-----  get correct index number
c
      if (ieps .eq. 0) then
        ilb = 39
      else
        ilb = 40
      endif
c
c-----  fill -bmat- array with proper elements from -relay-
c       that is, the lower rectangle
c
      do 100, i = 1, ndimb
        do 200, j = 1, ndimb
          bmat(i,j) = relay(npol3+i,npol3+j)
  200   continue
  100 continue
c
c-----  write -bmat- to dafile on record -39- (static) or -40- (optic)
c
      call dawrit(idafdrf,iodadrf,bmat,ndimb*ndimb,ilb,navdrf)
c
      return
      end
      subroutine bread(relay,bmat,npol3,ndimb,ieps)
c------
c      reads coupling between boundary elements to dafile
c      and puts them at the appropriate place in -relay-
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c----- dummy variables
c
      dimension relay(npol3+ndimb,npol3+ndimb),bmat(ndimb,ndimb)
c
c-----  common blocks
c
INCLUDE(comdrf/drfdaf)
c
c-----  get correct index
c
      if (ieps .eq. 0) then
        ilb = 39
      else
        ilb = 40
      endif
c
c-----  read -bmat- to dafile on record -39- (static) or -40- (optic)
c
      call daread(idafdrf,iodadrf,bmat,ndimb*ndimb,ilb)
c
c-----  fill -relay- array with proper elements from -bmat-
c       that is, the lower rectangle
c
      do 100, i = 1, ndimb
        do 200, j = 1, ndimb
          relay(npol3+i,npol3+j) = bmat(i,j)
  200   continue
  100 continue
c
      return
      end
      subroutine wtvrcal(xexp,wt,vr,
     1           xsurf,xnorm,area,ieps,ineq)
c------
c       calculation of (expanded) source and reaction fields
c       of/on (formal) qm particles (nuclei, electrons)
c       or representations (dp charges, mulliken charges
c       and dipoles) thereof, at polarizabilities,
c       boundary elements and external charges.
c
c       --------  p.th. van duijnen, ibm-kingston 1985, and
c                 groningen, dec. 1991.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension xsurf(3,nbem),xnorm(3,nbem),area(nbem)
      dimension xexp(3,nexp),wt(nwtr,nwtc),
     1          vr(nwtr,nwtc)
c
c-----  common blocks
c
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/mollab)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/grinf)
INCLUDE(comdrf/drfamb)
INCLUDE(comdrf/neqpar)
c
INCLUDE(comdrf/drfbem)
INCLUDE(comdrf/rad)
c
c-----  local variables
c
      logical kapnoz
c
      character*16 namj
c
      dimension p(3), q(3)
      dimension pq(3),w(3),b(3,3),bq(3),delkiq(3),delliq(3)
      dimension qi(3), t(3,3), tq(3), tni(3)
c
      character*8 errmsg(3)
c
      data errmsg/'program','stop in','-wtvrcal'/
      data third/.33333333333333333333333333333d00/
      data sixth/.16666666666666666666666666667d00/
      data two,three,four,twelve/2.0d00,3.0d00,4.0d00,12.0d00/
      data zero,pt5,pt75 /0.0d00,0.5d00,0.75d00/
      data one,onept5 /1.0d00,1.5d00/
c
c-----  begin
c
      if (ibem .ne. 0) then
c 1-----
        if (ineq .eq. 0) then
c   2-----
          if (ieps .eq. 0) then
c     3-----
            eps = eps1
            kappa = kappa1
          else
            eps = eps2
            kappa = kappa2
c     3-----
          endif
c   2-----
        else
c   2-----
          if (ieps .eq. 0) then
c     3-----
            eps = epsneq1
            kappa = kapneq1
          else
            eps = epsneq2
            kappa = kapneq2
c     3-----
          endif
c   2-----
        endif
c
c  -----  set logical  kapnoz
c         for non-zero ionic strength
c
        kapnoz = kappa .ne. zero
c
c  -----  initialize some important factors
c
        expkd = one
        fkapone = one
        pi=four*atan(one)
        epsfact= one/(two*pi*(one+eps))
        kappas = kappa**2
c 1-----
      endif
c
      call clear(wt(1,1),nexp4*ndim)
      call clear(vr(1,1),nexp4*ndim)
c
      call clear(wt(1,nexp4+ngran+1),ndim)
      call clear(vr(1,nexp4+ngran+1),ndim)
c
c-----  loop over the expansion centra
c
      do 500, j = 1, nexp
c 1-----
        if (j .le. nat) then
c   2-----
c    -----  nucleus at expanson centre
c
c    -----  skip ambiguous atoms
c
          if (nambpt(j) .ne. 0) goto 500
c
c    -----  nuclear charge -za-
c
          za = czan(j)
          namj = anam(j)
c
c    -----  polarisability of atom corresponding to expansion centre
c
          alfj = alfat(j)
c         alfj = alfa(namj,0,ier)
c   2-----
        else
c   2-----
c    ----- non-nulcei are given a polarizability 1.
c
          alfj = one
c   2-----
        endif
c
c  -----  pointer to arrays w.r.t. expansion centre -j-
c
        jf = (j-1)*4
c
c  -----  position vector -j- into -q-
c
        q(1) = xexp(1,j)
        q(2) = xexp(2,j)
        q(3) = xexp(3,j)
c
c  -----  loop over polarizable points
c
        do 300, ii = 1, npol
c   2-----
          np = mpol(ii)
c
c    -----  skip ambiguous points. if polarisability is required for
c           these points w.r.t. the qm system, add basis functions
c           to represent the polarisability
c
          if (ncutpt(np) .ne. 0) goto 300
c
c    -----  position vector into -p-
c
          p(1) = xpts(1,np)
          p(2) = xpts(2,np)
          p(3) = xpts(3,np)
c
c    -----  pointer in field arrays w.r.t. polarizable point
c
          if = (ii-1)*3
c
c    -----  calculate distance vector between -q- and -p-
c
          call distab(q,p,pq,dist)
c
c    -----  skip very close polarisabilities (as may occur for
c           ambiguous atoms)
c
          if (dist .le. 1.0d-03) goto 300
c
          dmind1 = one/dist
          dmin3 = dmind1*(dmind1**2)
          factp = one
          factf = one
          v = one
c
          alfi = alfext(np)
          if (alfi .eq. zero) alfi = one
c
c    -----  account for penetration effects (optional)
c
          if (modxza .ne. 0) then
c     3-----
            s = (alfj*alfi)**sixth
            v = dist/s
            if (ithole .eq. 1) then
c       4-----
              if (v .le. afact) then
                av = v/afact
                factp = av**4 - two*av**3 + two*av
                factf = four*av**3 - three*av**4
              endif
c       4-----
            else if (ithole .eq. 2) then
c       4-----
              au = afact*v
              factp = (one - (pt5*au + one)*exp(-au))
              factf = (one - (pt5*au**2 + au + one)*exp(-au))
c       4-----
            endif
c     3-----
          endif
          factp = factp*dmind1
          factf = factf*dmin3
c
c    -----  calculate field gradient tensor
c
          call clear(b,9)
c
c    -----  b(3,3) = t(p;q) : field gradient of charge in -p- at -q-
c           note that the interaction is scaled by v
c
          call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c    -----  calculate -w- = t(p;q) . q , this is part of the
c           taylor expansion of the field
c
          call matvec(b,q,w,3,.false.)
c
c    -----  -wt- and -vr- matrices for expanded field and potential
c           of charge (distribution) in expansion centra
c           at the polarizable points (-wt-) and vice versa (-vr-)
c
c           depending on the type of sources (charges, dipoles,
c           charge distributions), specified in rfin by ifldin
c           and "recipients", specified by ifldout, the matrices
c           are constructed partly or completely.
c
          do 290, k = 1, 3
c     3-----
c      -----  scale the distance vector
c
            pqk = pq(k)*factf
            do 280, l = 1, 3
c       4-----
c        -----  copy -b- into -wt-  and/or -vr-
c
c               this expansion can be used to calculate the field and
c               reaction field of a unit dipole in -q- (or -j-, the
c               expansion centre) at polarizable point -p-.
c               therefore, it is always calculated except when only
c               distributed monopoles are used to expand the source
c               and reaction fields of the quantum motif
c
c        -----  note: drftpq gives - d/dq f(q;p) !!
c
              if (ifldin .gt. 1) then
                wt(if+k,jf+l) = -b(k,l)
              endif
              if (ifldout .gt. 1) then
                vr(if+k,jf+l) = b(k,l)
              endif
c       4-----
  280       continue
c
c      -----  if the source/reaction field is not expanded (i.e. only
c             distributed mono- /dipoles are used), only the potential
c             part is stored in -wt-/-vr-
c
            if (ifldin .gt. 2) then
              wt(if+k,jf+4) = (w(k)+pqk)
            else
              wt(if+k,jf+4) = pqk
            endif
c
            if (ifldout .gt. 2) then
              vr(if+k,jf+4) = - (w(k) + pqk)
            else
              vr(if+k,jf+4) = -pqk
            endif
c
c      -----  form -zfn-: the fields in the polarizable points
c             due to the internal nuclei
c             -zfn- is in fact the sum of the nuclear fields
c             and potentials
c
            if (j .le. nat) then
              wt(if+k,nexp4+ngran+1) = wt(if+k,nexp4+ngran+1) + pqk*za
              vr(if+k,nexp4+ngran+1) = vr(if+k,nexp4+ngran+1) - pqk*za
            endif
c     3-----
  290     continue
c
c    -----  end of polarizable points
c   2-----
  300   continue
c
c  -----  fields (and fields dot normal) at boundary elements
c         due to charges, dipoles and charge distributions
c         (expanded) in the expansion centra
c
c  -----  loop over boundary elements
c
        do 400, ni = 1, nbem
c   2-----
c    -----  qi: vector from expansion centre -j- to
c           boundary element -ni- = (i-q)
c           dist: length of qi
c
          call distab(q,xsurf(1,ni),qi,dist)
          dmind1 = one/dist
          dmin2 = dmind1**2
          dmin3 = dmind1*dmin2
          if (kapnoz) then
            expkd = exp(-(kappa*dist))
            fkapone = one + (kappa*dist)
          endif
c
c    -----  fqini:  field of unit (positive) charge
c                   in -j- at -ni-, contracted with
c                   normal vector in -ni- = (i-q).n(i)/dist**3
c                   = f(q;i) . n(i)
c
c           it is also the negative of the potential of a dipole
c           in the direction of n(i) at the expansion centre
c
          fqini = dmin3*ddot(3,qi,1,xnorm(1,ni),1)
c         fqini = dmin3*adotb(qi,xnorm(1,ni),3)
c
c    -----  reaction potential energy operator
c           at expansion centre, first term in the expansion
c           containing contribution of unit dipoles w(i)
c
c                xkiq = k(i;q) s(i)
c
c        = (eps*(1+kappa*dist)*exp(-kappa*dist) - 1) f(i;q).n(i) s(i)
c
c    -----  the minus sign is a result of the use of the
c           inverted field: needed is f(i;q), fqini=f(q;i).n(i)
c
          xkiq = - (eps*fkapone*expkd - one)*fqini*area(ni)
c
          if (kapnoz) then
c
c      -----  reaction potential energy operator
c             at expansion centre, first term in the expansion,
c             containing contribution of unit charges z(i)
c
c                xliq = l(i;q) s(i)
c
c         = (1 - exp(-kappa*dist) ) v(i;q) s(i)
c
            xliq = (one - expkd)*dmind1*area(ni)
          endif
c
c    -----  check if expansion centre coincides with
c           position of nucleus; exclude ambiguous atoms
c
          if ((j .le. nat) .and. (nambpt(j) .eq. 0)) then
c     3-----
c      -----  potential of all source nuclear charges in -j- at -ni-
c             scaled with 1/(2pi(1+eps)), as input for coupling
c             equations for w(i)
c
c                  v(q;i) = sum(q) q(q) / dist
c
c      -----  this is the boundary element part of -zfn-
c
            if (ifldin .gt. 2) then
              wt(npol3+ni,nexp4+ngran+1) =
     1            wt(npol3+ni,nexp4+ngran+1) + epsfact*za*dmind1
            endif
c
c      -----  the reaction potential energy operator at the
c             expansion centre -j-, multiplied by the source
c             nuclear charge
c
c      -----  the interaction energy is evaluated through contracting
c             this array with the induced dipole density array
c
            if (ifldout .gt. 2) then
              vr(npol3+ni,nexp4+ngran+1) =
     1                    vr(npol3+ni,nexp4+ngran+1) + xkiq*za
            endif
c
            if (kapnoz) then
c       4-----
c        -----  minus field of all source nuclear charges in -j- at -ni-
c               contracted with normal vector at -ni-,
c               scaled with eps/(2pi(1+eps)), as input for
c               coupling equations for z(i)
c
c                    f(q;i) = q(q) (i-q)/ dist**3
c
              if (ifldin .gt. 2) then
                wt(npol3+ni+nbem,nexp4+ngran+1) =
     1          wt(npol3+ni+nbem,nexp4+ngran+1) - eps*epsfact*za*fqini
              endif
c
c        -----  the reaction potential energy operator at the
c               expansion centre -j-, multiplied by the source
c               nuclear charge
c
c        -----  the interaction energy is evaluated through contracting
c               this array with the induced charge density array
c
              if (ifldout .gt. 2) then
                vr(npol3+ni+nbem,nexp4+ngran+1) =
     1          vr(npol3+ni+nbem,nexp4+ngran+1) + xliq*za
              endif
c       4-----
            endif
c     3-----
          endif
c
c    -----  expansion of source potential and field
c           of surface charge distribution
c           and reaction potentials at and around the expansion centra
c
c             general:
c
c         expansion in taylor series in x around q:
c
c   for source potential: v(x;i) = v(q;i) + del(x) v(x;i) (x=q) .(x-q)
c         = v(q;i) - f(q;i).q + f(q;i).x
c
c   for source field:     f(x;i) = f(q;i) + del(x) f(x;i) (x=q) .(x-q)
c         = f(q;i) + t(q;i).q - t(q;i).x
c
c   for reaction potential due to induced dipoles (operator)
c                         k(i;x) = k(i;q) + del(x) k(i;x) (x=q) .(x-q)
c   = k(i;q) - [(eps*(1+kd)*exp(-kd) -1)*t(i;q).n(i) -
c                eps*(kappa**2)*exp(-kd)*f(i;q).n(i)*(q-i)].q
c            + [(eps*(1+kd)*exp(-kd) -1)*t(i;q).n(i) -
c                eps*(kappa**2)*exp(-kd)*f(i;q).n(i)*(q-i)].x
c
c   for reaction potential due to induced charges (operator)
c                         l(i;x) = l(i;q) + del(x) l(i;x) (x=q) .(x-q)
c   = l(i;q) - [(1+kd)*exp(-kd) -1)*f(i;q)].q
c            + [(1+kd)*exp(-kd) -1)*f(i;q)].x
c
c-----
c      source potential in -ni- (scaled with 1/(2pi(1+eps))
c
          if (ifldin .gt. 2) then
c-----
c      v(q;i) - f(q;i).q
c-----
            wt(npol3+ni,jf+4) = epsfact*(dmind1-dmin3
     +      *ddot(3,qi,1,q,1))
c    +      *adotb(qi,q,3))
          else
            wt(npol3+ni,jf+4) = epsfact*dmind1
          endif
c-----
c      the operator f(q;i) (.x)
c-----
          if (ifldin .gt. 1) then
            do 320, k = 1, 3
              wt(npol3+ni,jf+k) = epsfact*dmin3*qi(k)
  320       continue
          endif
c
c    -----  calculate del(i) f(q;i) = t(q;i) = t(i;q) = t
c
          call drftpq(qi,dmind1,dmin3,0,one,one,t)
c
c    -----  calculate tni = t(q;i) . n(i)
c
          call matvec(t,xnorm(1,ni),tni,3,.false.)
c
          if (kapnoz) then
c     3-----
c      minus source field in -ni- contracted with normal vector in -ni-
c                           and scaled with (eps/(2pi(1+eps))
c-----
c
c      -----  calculate contraction of t(q;i) with origin of
c             expansion q: tq = t(q;i).q
c
            call matvec(t,q,tq,3,.false.)
c-----
c      - [f(q;i).n(i) + [t(q;i).q].n(i)]
c-----
            if (ifldin .gt. 2) then
              wt(npol3+ni+nbem,jf+4) = - eps*epsfact*
     1        (fqini + ddot(3,tq,1,xnorm(1,ni),1))
c    1        (fqini + adotb(tq,xnorm(1,ni),3))
            else
              wt(npol3+ni+nbem,jf+4) = - eps*epsfact*fqini
            endif
c-----
c      the operator - [-t(q;i) (.x)] .n(i)
c                 =   [ t(q;i).n(i)] (.x)
c-----
            if (ifldin .gt. 1) then
              do 330, k = 1, 3
                wt(npol3+ni+nbem,jf+k) = eps*epsfact*tni(k)
  330         continue
            endif
c     3-----
          endif
c-----
c      reaction potential due to dipole in -ni- at -x-
c-----
c
c    -----  calculate del(x) k(i;x) (x=q) *s(i) = delkiq
c
c      the minus sign for the second term remains, since both
c      field dot normal(fqini) and vector (qi) are from -j- to -ni-
c      whereas the reverse of vectors is needed
c
          if (ifldout .gt. 1) then
            do 340, k = 1, 3
              delkiq(k)=  ((eps*fkapone*expkd - one)*tni(k) -
     1                  eps*kappas*expkd*fqini*qi(k))*area(ni)
c-----
c      the operator del(x) k(i;x) (x=q) *s(i) (.x)
c-----
              vr(npol3+ni,jf+k)= delkiq(k)
  340       continue
          endif
c
          if (ifldout .gt. 2) then
c-----
c      [k(i;q) - del(x) k(i;x) (x=q) .q] s(i)
c-----
            vr(npol3+ni,jf+4)=  xkiq - ddot(3,delkiq,1,q,1)
c           vr(npol3+ni,jf+4)=  xkiq - adotb(delkiq,q,3)
          else
c-----
c       k(i;q)
c-----
            vr(npol3+ni,jf+4)=  xkiq
          endif
c
          if (kapnoz) then
c     3-----
c      reaction potential due to charge in -ni- at -x-
c-----
c      calculate del(x) l(i;x) (x=q) *s(i) = delliq
c
c      minus sign is a result of using f(q;i) = (i-q)/dist**3,
c      whereas f(i;q) is needed
c
            if (ifldout .gt. 1) then
              do 350, k = 1, 3
                delliq(k)=-(one-fkapone*expkd)*qi(k)*dmin3*area(ni)
c-----
c      the operator [del(x) l(i;x) (x=q)] *s(i) (.x)
c-----
                vr(npol3+ni+nbem,jf+k) = delliq(k)
  350         continue
            endif
c
            if (ifldout .gt. 2) then
c-----
c      [l(i;q) - del(x) l(i;x) (x=q) . q] *s(i)
c-----
              vr(npol3+ni+nbem,jf+4)=  
     +        xliq - ddot(3,delliq,1,q,1)
c             vr(npol3+ni+nbem,jf+4)=  xliq - adotb(delliq,q,3)
            else
c-----
c       l(i;q)
c-----
              vr(npol3+ni+nbem,jf+4) = xliq
            endif
c     3-----
          endif
c   2-----
  400   continue
c 1-----
  500 continue
c
      return
      end
      subroutine dipsf(npol3,ndimb,nexp,wt,dip,sf)
c------
c      this routine calculates the field at the
c      polarizabilities and /or boundary charges
c      due to dipoles in the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension wt(npol3+ndimb,nexp*4+2),dip(3,nexp)
      dimension sf(npol3+ndimb)
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate pointer
c
        ip = (i-1)*4
c
c  -----  loop over x,y,z components
c
        do 200, k = 1, 3
c
c    -----  loop over polarizabilities (x,y,z) and boundary elements
c
          do 300, j = 1, npol3+ndimb
c
c      -----  calculate total field at polarizability/boundary element
c
            sf(j) = sf(j) + dip(k,i)*wt(j,ip+k)
  300     continue
  200   continue
c
  100 continue
      return
      end
      subroutine qsf(npol3,ndimb,nexp,wt,q,sf)
c------
c      this routine calculates the field at the
c      polarizabilities and /or boundary charges
c      due to charges in the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension wt(npol3+ndimb,nexp*4+2)
      dimension q(nexp), sf(npol3+ndimb)
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate pointer
c
        ip = (i-1)*4+4
c
c  -----  loop over polarizabilities (x,y,z) and boundary elements
c
        do 200, j = 1, npol3+ndimb
c
c    -----  calculate total field at polarizability/boundary element
c
          sf(j) = sf(j) + q(i)*wt(j,ip)
  200   continue
c
  100 continue
      return
      end
      subroutine elsf(npol3,ndimb,nexp,num,wt,d,s,rx,ry,rz,sf,iexpc)
c------
c      this routine calculates the field at the
c      polarizabilities and /or boundary charges
c      due to the electronic density
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension wt(npol3+ndimb,nexp*4+2)
      dimension d(num*(num+1)/2), rx(num*(num+1)/2),ry(num*(num+1)/2),
     1          rz(num*(num+1)/2),s(num*(num+1)/2)
      dimension iexpc(num*(num+1)/2)
      dimension sf(npol3+ndimb)
c
c----- local arrays
c
      dimension f(4)
c
      data one, two /1.0d00, 2.0d00/
c
c-----  effective loop over charge distributions
c
      ij = 0
      do 100, i = 1, num
        do 200, j = 1, i
          ij = ij + 1
c
          fact = two
          if (j .eq. i) fact = one
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
c    -----  premultiply density with dipole and overlap moments
c           of the charge distribution
c
          f(1) = fact*d(ij)*rx(ij)
          f(2) = fact*d(ij)*ry(ij)
          f(3) = fact*d(ij)*rz(ij)
          f(4) = fact*d(ij)*s(ij)
c
c    -----  execute expansion
c
          do 300, k = 1, 4
            do 400, l = 1, npol3+ndimb
              sf(l) = sf(l) - f(k)*wt(l,ip+k)
  400       continue
  300     continue
c
  200   continue
  100 continue
      return
      end
      subroutine diprf(npol3,ndimb,nexp,vr,dip,dipind,ediprf)
c------
c      this routine calculates the interaction of induced dipoles
c      and /or boundary charges
c      with dipoles at the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vr(npol3+ndimb,nexp*4+2),dip(3,nexp)
      dimension dipind(npol3+ndimb)
c
      data zero /0.0d00/
c
c-----  initialize energy contribution
c
      ediprf = zero
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate pointer
c
        ip = (i-1)*4
c
c  -----  loop over polarizabilities (x,y,z) and boundary elements
c
        do 200, j = 1, npol3+ndimb
c
c    -----  loop over x,y,z components
c
          do 300, k = 1, 3
c
c      -----  calculate interaction energy induced moments
c             with dipoles at expansion centra
c
            ediprf = ediprf + dip(k,i)*vr(j,ip+k)*dipind(j)
  300     continue
  200   continue
c
  100 continue
      return
      end
      subroutine qrf(npol3,ndimb,nexp,vr,q,dipind,eqrf)
c------
c      this routine calculates the field at the
c      polarizabilities and /or boundary charges
c      due to charges in the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vr(npol3+ndimb,nexp*4+2)
      dimension q(nexp), dipind(npol3+ndimb)
c
      data zero /0.0d00/
c
c-----  initialize energy contrbution
c
      eqrf = zero
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate pointer
c
        ip = (i-1)*4+4
c
c  -----  loop over polarizabilities (x,y,z) and boundary elements
c
        do 200, j = 1, npol3+ndimb
c
c    -----  calculate interaction energy of charges with rf
c
          eqrf = eqrf + q(i)*vr(j,ip)*dipind(j)
  200   continue
c
  100 continue
      return
      end
      subroutine elrf(npol3,ndimb,nexp,num,vr,d,rx,ry,rz,s,iexpc,
     1                dipind,eelmol)
c------
c      this routine calculates the interaction energy between
c      the induced dipoles and charges and the electronic
c      density
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vr(npol3+ndimb,nexp*4+2)
      dimension d(num*(num+1)/2), rx(num*(num+1)/2),ry(num*(num+1)/2),
     1          rz(num*(num+1)/2),s(num*(num+1)/2)
      dimension iexpc(num*(num+1)/2)
      dimension dipind(npol3+ndimb)
c
c----- local arrays
c
      dimension f(4)
c
      data zero, one, two /0.0d00, 1.0d00, 2.0d00/
c
c-----  initialize interaction energy
c
      eelmol = zero
c
c-----  effective loop over charge distributions
c
      ij = 0
      do 100, i = 1, num
        do 200, j = 1, i
          ij = ij + 1
c
          fact = two
          if (j .eq. i) fact = one
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
c    -----  premultiply density with dipole and overlap moments
c           of the charge distribution
c
          f(1) = fact*d(ij)*rx(ij)
          f(2) = fact*d(ij)*ry(ij)
          f(3) = fact*d(ij)*rz(ij)
          f(4) = fact*d(ij)*s(ij)
c
c    -----  execute expansion
c
          do 300, k = 1, 4
            do 400, l = 1, npol3+ndimb
              eelmol = eelmol - f(k)*vr(l,ip+k)*dipind(l)
  400       continue
  300     continue
c
  200   continue
  100 continue
      return
      end
      subroutine drfomga(ieps,relay,wt,vr,omega,indx)
c------
c      this routine constructs the matrix omega, containing the
c      reaction field contributions to the total energy for
c      the expanded electrons (formally for one electron per
c      expansion centre), the nuclei and the external charges.
c
c      it is formed by contracting the reaction potential -vr-
c      with the (formal) solutions of the relay equations, i.e.
c      the solutions of ri=s, with r the relay matrix, and s the
c      -wt- matrix.
c
c      depending on the type of source and recipients (i.e.
c      mono-/dipole expansion, molecular, separate nuclei and
c      electrons), parts or all of omega will be used for analysis.
c      also, the matrix elements of -wt- and -vr- differ, depending
c      on the type of expansion, resulting in differing omega
c      matrix elements. however, this does not influence the
c      calculation of omega. omega might contain mainly zero's.
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension relay(ndim,ndim), wt(nwtr,nwtc), vr(nwtr,nwtc),
     1 omega(nwtc,nwtc)
c
      dimension indx(ndim)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/auxdrf)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
c
INCLUDE(comdrf/runpar)
c
c-----  end of common blocks
c
c-----  data
c
      data zero,one,pt5,pt75 /0.0d00,1.0d00,0.5d00,0.75d00/
c
c-----  begin
c
      if (ieps .eq. 0) then
        ilr = 43
        ili = 44
        ilwt = 26
        ilvr = 28
      else
        ilr = 47
        ili = 48
        ilwt = 30
        ilvr = 32
      endif
c
      if (mcupdt) then
        ilwt = ilwt + 1
        ilvr = ilvr + 1
      endif
c
c-----  read the lu-decomposed relay matrix from record -ilr-
c       and array -indx- from record -ili-
c
      call daread(idafdrf,iodadrf,relay,ndim*ndim,ilr)
      call daread(idafdrf,iodadrf,indx,ndim,ili)
c
      if (idrfout.eq.3 .or. imcout .eq. 5) then
        call hatout(relay,ndim,ndim,2,'relay-mat')
        call imatou(indx,ndim,1,2,'relay-indx')
      endif
c
c-----  read -wt- and -vr- from disk
c
      call daread(idafdrf,iodadrf,wt,nwtr*nwtc,ilwt)
      call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvr)
c
      call clear(omega,nomga)
c
c-----  solve linear eqs. and form -omega-
c
      do 100, i = 1, nwtc
        call luelmf(relay,wt(1,i),indx,ndim,ndim,auxx)
        if (idrfout.eq.3 .or. imcout .eq. 5) then
          call hatout(auxx,1,ndim ,22,'auxx_uit')
        endif
c
c  -----  evaluation of energy contributions after solving
c         ri=s
c
c  -----  first for electrons (expansion, expansion centre)
c         then for external charges
c
        do 200, j = 1, nwtc
          omega(i,j) = ddot(ndim,auxx,1,vr(1,j),1)
c         omega(i,j) = adotb(auxx,vr(1,j),ndim)
  200   continue
  100 continue
c
      if (idrfout.eq.3 .or. imcout .eq. 5)
     1        call hatout(omega,nwtc,nwtc,2,'omega')
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c-----  form -wt*a- for use in derivative
c
c
c     do 700 i=1,nexp4
c       if(ibem.eq.0) then
c         call fdagal(wt(1,i),a,auxx,ia,npol3)
c       else
c         andere transform
c       endif
c       do 650 k=1,npol3+ndimb
c         wt(k,i)=auxx(k)
c 650   continue
c 700 continue
c
c
c     if(idrfout.eq.5) call hatout(wt  ,(npol3+ndmib),nexp4,2,'wt*a ')
c
c-----  save  -zfp*a- separately on disk.
c
c     call dawrit(idafdrf,iodadrf,auxx,npol3+ndimb,5,navdrf)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      return
      end
      subroutine drfoab(iexpa,iexpb,nc,omega)
c------
c
c       routine fetches required part of the -omega- matrix.
c       the required part is a 4x4 matrix. if a diagonal
c       block is required, the lower triangle of -omega-
c       is expanded to square, symmetric form
c
c       iexpa,iexpb label the expansion centra
c
c       the required part is put into -omgab-, which
c       is transported through the common block drfexp
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension omega(nc,nc)
c
c-----  common blocks
c
INCLUDE(comdrf/drfexp)
c
c-----  begin
c
      ia = iexpa
      ib = iexpb
c
c-----  set pointers
c
      ik = (ia-1)*4
      jk = (ib-1)*4
c
      do 20, k = 1, 4
        do 30, l = 1, 4
          omgab(k,l) = omega(ik+k,jk+l)
   30   continue
   20 continue
      return
      end
         
