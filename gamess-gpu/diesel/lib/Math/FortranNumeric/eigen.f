*###### bibliotheksroutinen fuer mrdci-program
*
* v0.0a (21.8.1992)
*
*
c
      subroutine eigen(a,r,n,mv)
c
c primitive jacobi-diagonalisierung
c
c
      integer ij,j,i,n
      dimension a(n*(n+1)/2),r(n*n)
      double precision a,r,anorm,anrmx,thr,x,y,sinx,sinx2,cosx,cosx2,
     1 sincs
      if(mv-1) 10,25,10
  10  iq=-n
      do j=1,n
       iq = iq + n
       do i=1,n
        ij = iq + i
        r(ij)=0.0d0
        if(i.eq.j) r(ij)=1.0d0
       end do
      end do
  25  anorm=0.0d0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
  30  ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
  35  continue
      if(anorm) 165,165,40
  40  anorm=1.414d0*dsqrt(anorm)
      anrmx=anorm*1.0d-12/float(n)
      ind=0
      thr=anorm
  45  thr=thr/float(n)
  50  l=1
  55  m=l+1
  60  mq=((m*m)-m)/2
      lq=((l*l)-l)/2
      lm=l+mq
  62  if(dabs(a(lm))-thr) 130,65,65
  65  ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
  68  y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
  70  y=-y
  75  if(y .gt. 1.0d0) y=1.0d0
      if(y .lt. -1.0d0) y=-1.0d0
      sinx=y/dsqrt(2.0d0*(1.0d0+(dsqrt(1.0d0-y*y))))
      sinx2=sinx*sinx
  78  cosx=dsqrt(1.0d0-sinx2)
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
