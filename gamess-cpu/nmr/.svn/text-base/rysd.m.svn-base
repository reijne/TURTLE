      subroutine droot
      implicit REAL  (a-h,o-z)
      REAL xx,uf,wf
      common/hfk/xx,uf(9),wf(9),nroots
      common/rysff/ff(19)
      common/rysrw/r(9,9),w(9,9)
      dimension c(10,10),s(10,10),a(10),rt(10)
      data pt5,zero,one,four /0.5d+00,0.0d+00,1.0d+00,4.0d+00/
c
c     this version uses christoffel formula for weights.
c     ith root of the jth rys polynomial is returned in r(i,j) with
c     the corresponding weight factor in w(i,j).   j=1,2,...,n
c
      n=nroots
      x=xx
      if(n.lt.2) n=2
      n1=n+1
      nn=n+n
      call rysfun(x,nn)
      do 10 i=1,n1
      do 10 j=1,n1
   10 s(i,j)=ff(i+j-1)
      call ryssmt(c,s,n1)
      do 20 i=1,n
      do 20 j=1,i
      w(i,j)= zero
   20 r(i,j)= zero
      wsum=ff(1)
      w(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum= sqrt(c(2,3)**2-four *c(1,3)*c(3,3))
      r(1,2)=   pt5 *(-c(2,3)-dum)/c(3,3)
      r(2,2)=   pt5 *(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 70
      do 25 i=3,n1
   25 rt(i)=  one
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 60 k=3,n
      k1=k+1
      do 30 i=1,k1
   30 a(i)=c(i,k1)
      call rysnod(a,rt,k)
      do 50 i=1,k
   50 r(i,k)=rt(i)
   60 continue
   70 do 150 k=2,n
      jmax=k-1
      do 150 i=1,k
      root=r(i,k)
      dum=  one  /ff(1)
      do 110 j=1,jmax
      j1=j+1
      poly=c(j1,j1)
      do 100 m=1,j
  100 poly=poly*root+c(j1-m,j1)
  110 dum=dum+poly*poly
  150 w(i,k)=  one  /dum
      do 160 k=1,nroots
      dum=r(k,nroots)
      uf(k)=dum/(  one  -dum)
  160 wf(k)=w(k,nroots)
      return
      end
      subroutine ryssmt(c,s,n)
      implicit REAL  (a-h,o-z)
      dimension c(10,10),s(10,10),v(10),y(10)
      data zero,one /0.0d+00,1.0d+00/
c
c     routine returns an n by n triangular matrix c such that
c     c(transpose)sc=i,  where i is an n by n identity matrix.
c
      do 10 i=1,n
      do 10 j=1,i
   10 c(i,j)= zero
      do 100 j=1,n
      kmax=j-1
      fac=s(j,j)
      if(kmax.eq.0) go to 60
      do 20 k=1,kmax
      v(k)= zero
   20 y(k)=s(k,j)
      do 50 k=1,kmax
      dot= zero
      do 30 i=1,k
   30 dot=c(i,k)*y(i)+dot
      do 40 i=1,k
   40 v(i)=v(i)-dot*c(i,k)
   50 fac=fac-dot*dot
   60 fac=one/ sqrt(fac)
      c(j,j)=fac
      if(kmax.eq.0) go to 100
      do 70 k=1,kmax
   70 c(k,j)=fac*v(k)
  100 continue
      return
      end
      subroutine rysfun(x,n)
      implicit REAL  (a-h,o-z)
      common/rysff/ff(19)
      data pt5,one,two /0.5d+00,1.0d+00,2.0d+00/
      tol=1.0d-14
      xx=x+x
      facmin=xx
      e=0.5409855304296342219319112d-78
      if(facmin.lt.2*180.2160d+00) e= exp(-x)
      if(facmin.gt.   80.0000d+00) go to 100
      term=one
      sum =one
      fac=n
      fac=fac+pt5
   10 fac=fac+one
      term=term*x/fac
      sum=sum+term
      if(fac.le.facmin) go to 10
      t=term
      s=sum
      if(t.gt.s*tol) go to 10
      fac=n+n+1
      ff(n+1)=sum*e/fac
      m=n-1
      fac=m+m+1
   20 if(m.lt.0) return
      ff(m+1)=(e+xx*ff(m+2))/fac
      m=m-1
      fac=fac-two
      go to 20
c
c     use asymptotic expansion for large arguments.
c
  100 a= sqrt(.7853981633974483096156608d+00/x)
      tmax=a*tol/e
      term=one/xx
      sum=term
      fac=one
  110 fac=fac-two
      term=fac*term/xx
      sum=term+sum
      t=term
      if( abs(t).gt.tmax) go to 110
      ff(1)=a-e*sum
      fac=-one
      m=0
  120 if(m.eq.n) return
      m=m+1
      fac=fac+two
      ff(m+1)=(fac*ff(m)-e)/xx
      go to 120
      end
      subroutine rysnod(a,rt,k)
      implicit REAL  (a-h,o-z)
      character*8 errmsg
      dimension a(10),rt(10)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-rysnod-'/
      data zero /0.0d+00/
c
c     routine returns rt(i) the ith root of a polynomial of order
c     k whose mth coefficient is stored in a(m+1).  it is assumed that
c     the initial values in rt bracket the final values.
c
      tol=1.0d-11
      k1=k+1
      r2= zero
      p2=a(1)
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt. zero) go to 20
      write(6,15) m,k
   15 format(//,' root number ',i3,' was not found for polynomial',
     1 ' of order ',i3,//)
      call hnderr(3,errmsg)
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if( abs(delta).lt.tol) go to 90
      dr=0.0625d+00*dr
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
      p5=p5*r5+a(k1-i)
   40 p6=p6*r6+a(k1-i)
   45 prod=p5*p6
      if(prod.lt. zero) go to 30
      prod=p3*p5
      if(prod.gt. zero) go to 60
      r5=0.25d+00*r3+0.75d+00*r5
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
   60 r6=0.25d+00*r4+0.75d+00*r6
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
      return
      end
