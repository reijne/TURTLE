c 
c  $Author: wab $
c  $Date: 2009-12-17 11:37:45 +0100 (Thu, 17 Dec 2009) $
c  $Locker:  $
c  $Revision: 6115 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util3.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   util3   =
c ******************************************************
c ******************************************************
c     deck=util3
_IF(vax)
      subroutine fmove(a,b,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      n5=(n/5)*5
      if(n5.le.0)goto 2
      do 1 i=1,n5,5
      b(i  ) = a(i  )
      b(i+1) = a(i+1)
      b(i+2) = a(i+2)
      b(i+3) = a(i+3)
 1    b(i+4) = a(i+4)
 2    if(n.le.n5)return
      n5=n5+1
      do 3 i=n5,n
 3    b(i)=a(i)
      return
      end
      function vecsum(a,b,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      vecsum=0.0d0
      n5=(n/5)*5
      if(n5.le.0)goto 2
      do 1 i=1,n5,5
    1 vecsum=a(i)*b(i)+a(i+1)*b(i+1)
     *           +a(i+2)*b(i+2)+a(i+3)*b(i+3)+a(i+4)*b(i+4)+vecsum
    2  if(n.le.n5)return
      n5=n5+1
      do 3 i=n5,n
    3 vecsum=a(i)*b(i)+vecsum
      return
      end
      subroutine mult2(q,a,h,ncore,ndum,nbasis)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension a(*),q(*),h(*)
      dimension p(maxorb),y(maxorb)
INCLUDE(common/mapper)
c...   a=q(transpose) * h * q
c...   a and h stored in triangle form
      do 1 j=1,ncore
      call dcopy(nbasis,q(ilifq(j)+1),1,y(1),1)
      do 2 i=1,nbasis
      bot=y(i)
      m=iky(i)
      top=h(m+i)*bot
      n=i-1
      if(n)2,2,4
4     do 5 k=1,n
       gut=h(m+k)
      top=y(k)*gut+top
5     p(k)=bot*gut+p(k)
2     p(i)=top
      m=iky(j)
      do 1 i=1,j
      iiii=ilifq(i)
      sum=0.0d0
      do 11 k=1,nbasis
 11   sum=sum+p(k)*q(iiii+k)
  1   a(m+i)=sum
      return
      end
      subroutine szero(c,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension c(*)
      if(n.gt.0)then
      do 1 loop=1,n
    1 c(loop)=0.0d0
      endif
      return
      end
      subroutine gmake(a,p,iky)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer * 2 i,j,k,l
      dimension a(*),p(*),iky(*)
      logical *1 iii(2),jjj(2),kkk(2),lll(2),mij,mkl
      common/blkin/g(340),mij(680),mkl(680),mword
      data i,j,k,l/4*0/
      equivalence(iii(1),i),(jjj(1),j),(kkk(1),k),(lll(1),l)
      equivalence (gg1,gg),(ikyj,jl),(ikyk,jk),(ikyi,il)
c ...
c ... computes symmetric g matrix
c ... g=2j<p>-k<p>
c ... p is symmetric density matrix
c ...
      do 60 iw=1,mword
      gg=g(iw)
      gg2=gg
      gg3=gg
      gg4=gg
      gg5=gg
      gg6=gg
      gg7=gg
      gg8=gg+gg
      gg9=gg8
      il=iw+iw
c
c vax ... i*2 i*256+j equivalenced to 2 logical*1 s
c vax ... has first logical=j, second=i.
c
      iii(1)=mij(il)
      jjj(1)=mij(il-1)
      kkk(1)=mkl(il)
      lll(1)=mkl(il-1)
      ikyi=iky(i)
      ij=ikyi+j
      ik=ikyi+k
      il=ikyi+l
      ikyk=iky(k)
      kl=ikyk+l
      ikyj=iky(j)
      if(i.eq.j)then
      gg3=0.0
      gg4=0.0
      gg6=0.0
      gg7=0.0
      gg9=gg
      endif
      if(ij.eq.kl)then
      gg5=0.0
      gg6=0.0
      gg7=0.0
      gg9=0.0
      endif
      if(k.eq.l)then
      gg2=0.0
      gg4=0.0
      gg7=0.0
      gg8=gg
      endif
      if(i.eq.k)gg1=gg1+gg5
      if(j-k)4,5,6
    5 gg3=gg3+gg6
      go to 7
    4 gg3=gg6
    7 jk=ikyk+j
      if(j-l)8,9,10
    9 gg4=gg4+gg7
      go to 10
    8 gg4=gg7
      jl=iky(l)+j
      go to 11
    6 jk=ikyj+k
   10 jl=ikyj+l
   11 a(ij)=(gg8+gg8)*p(kl)+a(ij)
      a(kl)=(gg9+gg9)*p(ij)+a(kl)
      a(ik)=a(ik)-gg1*p(jl)
      a(jl)=a(jl)-gg4*p(ik)
      a(il)=a(il)-gg2*p(jk)
      a(jk)=a(jk)-gg3*p(il)
   60 continue
      return
      end
      subroutine proc2(a,p,b,q,iky)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *2 i,j,k,l
      dimension a(*),p(*),b(*),q(*),iky(*)
      logical *1 iii(2),jjj(2),kkk(2),lll(2),mij,mkl
      common/blkin/g(340),mij(680),mkl(680),mword
      data i,j,k,l/4*0/
      equivalence(iii(1),i),(jjj(1),j),(kkk(1),k),(lll(1),l)
c ... hamiltonian builder for rhf open shell .. spatially restricted
c ... define following symmetric density matrices
c ...  r1 .. closed shell
c ...  r2 .. open shell
c ... p=r1+r2/2
c ... setup ..
c ...     2j<p>-k<p> in a .... p=p
c ...      k<r2>/2   in b ....q=r2/2
      do 60 iw=1,mword
      gg=g(iw)
      gg1=gg
      gg2=gg
      gg3=gg
      gg4=gg
      gg5=gg
      gg6=gg
      gg7=gg
      gg8=gg+gg
      gg9=gg8
      il=iw+iw
      iii(1)=mij(il)
      jjj(1)=mij(il-1)
      kkk(1)=mkl(il)
      lll(1)=mkl(il-1)
      ikyi=iky(i)
      ij=ikyi+j
      ik=ikyi+k
      il=ikyi+l
      ikyk=iky(k)
      kl=ikyk+l
      ikyj=iky(j)
      if(i.eq.j)then
      gg3=0.0
      gg4=0.0
      gg6=0.0
      gg7=0.0
      gg9=gg
      endif
      if(ij.eq.kl)then
      gg5=0.0
      gg6=0.0
      gg7=0.0
      gg9=0.0
      endif
      if(k.eq.l)then
      gg2=0.0
      gg4=0.0
      gg7=0.0
      gg8=gg
      endif
      if(i.eq.k)gg1=gg1+gg5
      if(j-k)4,5,6
    5 gg3=gg3+gg6
      go to 7
    4 gg3=gg6
    7 jk=ikyk+j
      if(j-l)8,9,10
    9 gg4=gg4+gg7
      go to 10
    8 gg4=gg7
      jl=iky(l)+j
      go to 11
    6 jk=ikyj+k
   10 jl=ikyj+l
   11 a(ij)=(gg8+gg8)*p(kl)+a(ij)
      a(kl)=(gg9+gg9)*p(ij)+a(kl)
      a(ik)=a(ik)-gg1*p(jl)
      a(jl)=a(jl)-gg4*p(ik)
      a(il)=a(il)-gg2*p(jk)
      a(jk)=a(jk)-gg3*p(il)
      b(ik)=b(ik)-gg1*q(jl)
      b(jl)=b(jl)-gg4*q(ik)
      b(il)=b(il)-gg2*q(jk)
   60 b(jk)=b(jk)-gg3*q(il)
      return
      end
_ENDIF
_IF1()      subroutine ifmove(ia,ib,n)
_IF1()      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      dimension ia(*),ib(*)
_IF1()      if(n.le.0)return
_IF1()      do 1 loop=1,n
_IF1()    1 ib(loop)=ia(loop)
_IF1()      return
_IF1()      end
_IF(ibm,vax,1s)
      subroutine subvec(r,a,b,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*),r(*)
      if(n.gt.0) then
      do 1 loop=1,n
    1 r(loop)=a(loop)-b(loop)
      endif
      return
      end
      subroutine gtriad(n,scalar,r,a,b)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),b(*)
      if(n.gt.0)then
      do 1 loop=1,n
    1 r(loop)=a(loop)+scalar*b(loop)
      endif
      return
      end
      subroutine scaler(n,scalar,r,b)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),b(*)
      if(n.gt.0)then
      do 1 loop=1,n
    1 r(loop)=scalar*b(loop)
      endif
      return
      end
      subroutine vvtv(n,r,a,b)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      dimension a(*),b(*),r(*)
      if(n.gt.0) then
      do 1 loop=1,n
    1  r(loop)=a(loop)*b(loop)
      endif
      return
      end
_ENDIF
_IFN(cray,apollo,hp700)
      subroutine setsto(nf,itext,ilab)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ilab(*)
      if(nf.gt.0)then
      do 1 i=1,nf
   1  ilab(i)=itext
      endif
      end
_ENDIF
_IF(apollo,hp700)
       subroutine `setsto(nf,itext,ilab)'
       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
       implicit character *8 (z),character *1 (x)
       implicit character *4 (y)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
       dimension ilab(*)
       if(nf.gt.0) then
           call vec_$iinit(ilab,nf,itext)
       endif
       return
      end
_ENDIF
_IFN1(c)       subroutine setstl(nf,otext,olab)
_IFN1(c)       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IFN1(c)       implicit character *8 (z),character *1 (x)
_IFN1(c)       implicit character *4 (y)
_IFN1(c)       dimension olab(*)
_IFN1(c)       if(nf.gt.0)then
_IFN1(c)       do 1 loop=1,nf
_IFN1(c)    1  olab(loop)=otext
_IFN1(c)       endif
_IFN1(c)       return
_IFN1(c)       end
      subroutine triangle(r,a,n)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c.... transform symmetric square r to lower triangle a
c.... a and r may be the same
c
      dimension r(*),a(*)
c
      ka = 0
      do i=1,n
         kr = (i-1)*n
         do j=1,i
            ka = ka + 1
            a(ka) = r(kr+j)
         end do
      end do
c
      return
      end
_IFN(cray)
      subroutine square(r,a,mrowr,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*)
c... convert triangle a to square r
      iii=0
      jjj=1
      k=1
      do 4 i=1,n
      jj=jjj
_IF1(t)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 3 j=1,i
      r(iii+j)=a(k)
      r(jj)=a(k)
      k=k+1
    3 jj=jj+mrowr
      iii=iii+mrowr
    4 jjj=jjj+1
      return
      end
      subroutine sqtrip(r,a,nk)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(nk,*),a(*)
      do 1 i=1,nk
    1 r(i,i)=0.0d0
      m=1
      do 2 j=2,nk
      j1=j-1
_IF1(t)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 2 i=1,j1
      r(i,j)=a(m)
      r(j,i)=-a(m)
    2 m=m+1
      return
      end
_ENDIF
      subroutine trsqsq(sq,tr,n)
      implicit REAL  (a-h,o-z)
      dimension sq(n,n),tr(*)
      ij = 0
      do 30 i = 1 , n
         do 20 j = 1 , i
            ij = ij + 1
            tr(ij) = (sq(i,j)+sq(j,i))*0.5d0
 20      continue
 30   continue
      return
      end
_IF1(civ)      subroutine vsav(n,scalar,r,b)
_IF1(civ)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(civ)      implicit character *8 (z),character *1 (x)
_IF1(civ)      implicit character *4 (y)
_IF1(civ)      dimension b(*),r(*)
_IF1(civ)      if(n.gt.0) then
_IF1(civ)      do 1 loop=1,n
_IF1(civ)    1  r(loop)=b(loop)+scalar
_IF1(civ)      endif
_IF1(civ)      return
_IF1(civ)      end
_IFN(j90)
      function locat1(label,nf,itext)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(f)      word label,itext
      dimension label(*)
      do 1 loop=1,nf
      if(label(loop).eq.itext)go to 2
    1 continue
      locat1=0
      return
    2 locat1=loop
      return
      end
_ENDIF
_IF(cray)
      subroutine minvrt(a,n,d,l,m)
      implicit none
      REAL d, a(*)
      integer n, l(*), m(*)
c
c call cray routine
c   NB scratch space m is not required
c..   if singularities are not handled correctly cf minvrt below (jvl)
c
      call minv(a,n,n,l,d,1.0d-50,0,1)
      end
_ELSE
      subroutine minvrt(a,n,d,l,m)
c
c...  routine made to go on if singularity is detected 
c...  this is generally the best atitude / do set d to 0.0 then
c...  J.H. van Lenthe / R.W.A. Havenith  1998
c
c...  If on entry d.eq.-1.0d0 then don't calculate the determinant
c...  of the matrix otherwise do calculate the determinant.
c...  I assume that the determinant is computed to provide a test
c...  whether the matrix was singular. Perhaps someone wants
c...  the true determinant for other purposes, so we keep the option
c...  to calculate it. However, there are cases where computing the
c...  true determinant causes floating point overflows. So we allow
c...  its calculation to be switched off (e.g. in ZORA).
c...  Huub van Dam 1999.
c
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),l(*),m(*)
      data done,dzero/1.0d0,0.0d0/
      odeterminant = .not.(d.eq.-1.0d0)
      d=1.0d0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
      if( dabs(biga)- dabs(a(ij))) 15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
   45 if(biga) 48,46,48
c...   singularity ; fix by setting biga to 1.0
   46 biga = 1.0d0
      d=dzero
c
   48 do 55 i=1,n
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
      kj=k-n
      do 75 j=1,n
      kj=kj+n
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
      if (odeterminant) d=d*biga
      a(kk)=done/biga
   80 continue
      k=n
   90 k=(k-1)
      if(k) 92 ,92 ,93
   93 i=l(k)
      if(i-k) 94 ,94 ,95
   95 jq=n*(k-1)
      jr=n*(i-1)
      do 96  j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
   96 a(ji) =hold
   94 j=m(k)
      if(j-k)  90, 90, 91
   91 ki=k-n
      do 97  i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
   97 a(ji) =hold
      go to  90
   92 return
      end
_ENDIF
_IFN(vax)
c     ****f* transform/mult2
c
c     NAME
c
c       mult2 - wrapper to decide whether to parallel or serial
c               implementation
c
c     SYNOPSIS
c
c       mult2(q,r,f,na,nb,nbasis)
c
c     FUNCTION
c
c       Both a serial and parallel implementation of mult2 are
c       available. The parallel version is a sensible choice only
c       if there are 2 or more processors available to split the
c       work across and if there is enough work to do.
c
c     SOURCE
c
      subroutine mult2(q,r,f,na,nb,nbasis)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),r(*),f(*)
c
INCLUDE(common/parcntl)
c
      integer  ipg_nnodes
      external ipg_nnodes
c
      nproc = ipg_nnodes()
      if (nproc.gt.1.and.nbasis.gt.idpmult2) then
        call pmult2(q,r,f,na,nb,nbasis)
      else
        call smult2(q,r,f,na,nb,nbasis)
      endif
c
      end
c     ******
c
c     ****f* transform/pmult2
c
c     NAME
c
c       pmult2 - implements a matrix transformation in parallel
c
c     SYNOPSIS
c
c       pmult2(q,r,f,na,nb,nbasis)
c
c     ARGUMENTS
c
c       Q - An nbasis x na matrix holding the transformation matrix
c           which remains unchanged.
c       R - An na x nbasis sized memory segment, on return it holds
c           the upper triangle of the transformed matrix.
c       F - An max(na x na, nbasis x (nbasis + 1) / 2) memory segment,
c           on entry it holds the upper triangle of the matrix to
c           transform, on return the contents are destroyed.
c       na - The dimension of the space in which R will be represented
c       nb - If nb == na: R includes the diagonal,
c            If nb == na-1: R excludes the diagonal.
c       nbasis - The dimension of the space in which F is represented
c                on entry.
c
c     FUNCTION
c
c       The routine computes the matrix transformation
c
c         |latex \begin{equation}
c         |latex R = Q^T F Q
c         |latex \end{equation}
c         |html  R = Q<sup>T</sup> F Q
c
c       in parallel using the BLAS matrix-matrix multiplication
c       routines.
c
c       For the parallelisation the fact that all the matrices are
c       replicated data objects is exploited. Therefore the only 
c       concerns are the distribution of work and the re-replication
c       of the results.
c
c       The re-replication of results involves a 
c       log2(nproc)*O(nbasis**2) cost which is incurred for both 
c       matrix-matrix multiplications. However this is deemed 
c       acceptable as each matrix-matrix multiplication involves
c       O(nbasis**3) compute cost. 
c
c       The work is partitioned by splitting one dimension across the
c       processors. Therefore the work per processor scales as
c       O(nbasis**2) even in the most favourable case. This is
c       considered acceptable as the cost of the re-replication 
c       unavoidably has a similar scaling. Hence a more fine grained
c       partitioning of the work would not result in a better
c       scalability.
c
c       The first matrix-matrix multiplication involves a triangular
c       matrix. As the lower columns are much shorter than the higher
c       columns this could cause load balancing problems. However,
c       realize that the combined length of the first and last column
c       is the same as the combined length of the second and last but
c       one column, etc. Therefore by combining the work of two columns
c       a constant task size can be obtained. This helps in load
c       balancing the work even if static load balancing is used.
c       The code to exploit this is available below. However in
c       practice it does not seem as effective as the simpler version.
c       This may be due to cache misses which will reduce as the
c       matrix-matrix multiplication progresses in the current 
c       algorithm. With the more advanced algorithm the constant
c       alternation between a short vector and a long vector could
c       cause a high level of cache misses to persist throughout
c       the matrix-matrix multiplication.
c
c     SOURCE
c
      subroutine pmult2(q,r,f,na,nb,nbasis)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(parallel)
      common/scftim/tdiag(5),tmult,tdum(4)
_ENDIF
INCLUDE(common/timeperiods)
INCLUDE(common/parallel)
      dimension q(*),r(*),f(*)
c
      integer  ipg_dlbtask
      external ipg_dlbtask
c...
      nbb=nbasis
_IF(parallel)
      dumtim=dclock()
_ENDIF
      call start_time_period(TP_MULT2)
c...
c...  halve diagonal of f (only nbasis operations do not parallelise)
c...
      m=0
      do loop=1,nbasis
        m=m+loop
        f(m)=f(m)*0.5d0
      enddo
c...
c...  r = q(dagger) * f(scaled)
c...  lower triangle of  f  effectively zero
c...
      call vclr(r,1,nbasis*na)
      ll = nbasis*na+1
      mm = nbasis*(nbasis+1)/2+1
      next = ipg_dlbtask()
      do loop = nbasis, 1, -1
        mm = mm - loop
        ll = ll - na
        icount_dlb = icount_dlb + 1
        if (icount_dlb.eq.next) then
          call mxmb(q,nbb,1,f(mm),1,1,r(ll),1,1,na,loop,1)
          next = ipg_dlbtask()
        endif
      enddo
c
c     HvD Begin
c
c     The following section of code is a variant of the above that
c     was intended to demonstrate better load balancing properties.
c     The mxmb calls are structured such that every task given to
c     a processor is formally of equal size to all other tasks.
c     In practice however this performs worse, I assume do to 
c     cache thrashing. When I tested it on 4 processors the CPU
c     efficiency for a whole valinomycin calculation was reduced
c     from 98% with the code above to 96%, or stated otherwise
c     we lost about 90 seconds on 3900 by using the code below.
c     So for now we stick with the code above.
c
c     next = ipg_dlbtask()
c     nhalf = (nbasis+1)/2
c     do lbot = 1, nhalf
c       icount_dlb = icount_dlb + 1
c       if (icount_dlb.eq.next) then
c         ltop = nbasis+1-lbot
c         mmb = lbot*(lbot-1)/2+1
c         mmt = ltop*(ltop-1)/2+1
c         llb = (lbot-1)*na+1
c         llt = (ltop-1)*na+1
c
c...      The two calls to mxmb combine to constant task sizes.
c...      This helps with the load balancing even if static load
c...      balancing is used.
c
c         call mxmb(q,nbb,1,f(mmb),1,1,r(llb),1,1,na,lbot,1)
c         if (lbot.ne.ltop) then
c           call mxmb(q,nbb,1,f(mmt),1,1,r(llt),1,1,na,ltop,1)
c         endif
c         next = ipg_dlbtask()
c       endif
c     enddo
c
c     HvD End
c
      call pg_dgop(3010,r,nbasis*na,'+')
c...
c...  f = r * q
c...
      call vclr(f,1,na*na)
      do loop = 1, na
        icount_dlb = icount_dlb + 1
        if (icount_dlb.eq.next) then
          iiq = (loop-1)*nbb+1
          iif = (loop-1)*na+1
          call mxmb(r,1,na,q(iiq),1,nbb,f(iif),1,na,na,nbasis,1)
          next = ipg_dlbtask()
        endif
      enddo
      call pg_dgop(3020,f,na*na,'+')
c...
c...  r = f + f(dagger) (O(na**2) operations, do not parallelise)
c...
      m=0
      ll=0
      do loop=1,na
        mm=loop
        limit=loop
        if(loop.eq.na)limit=nb
        do moop=1,limit
          m=m+1
          r(m)=f(ll+moop)+f(mm)
          mm=mm+na
        enddo
        ll=ll+na
      enddo
      call end_time_period(TP_MULT2)
_IF(parallel)
      tmult=tmult+(dclock()-dumtim)
_ENDIF
      call pg_dlbpush
      return
      end
c     ******
c
c     ****f* transform/smult2
c
c     NAME
c
c       smult2 - implements a matrix transformation in serial
c
c     SYNOPSIS
c
c       smult2(q,r,f,na,nb,nbasis)
c
c     ARGUMENTS
c
c       Q - An nbasis x na matrix holding the transformation matrix
c           which remains unchanged.
c       R - An na x nbasis sized memory segment, on return it holds
c           the upper triangle of the transformed matrix.
c       F - An max(na x na, nbasis x (nbasis + 1) / 2) memory segment,
c           on entry it holds the upper triangle of the matrix to
c           transform, on return the contents are destroyed.
c       na - The dimension of the space in which R will be represented
c       nb - See na.
c       nbasis - The dimension of the space in which F is represented
c                on entry.
c
c     FUNCTION
c
c       The routine computes the matrix transformation
c
c         |latex \begin{equation}
c         |latex R = Q^T F Q
c         |latex \end{equation}
c         |html  R = Q<sup>T</sup> F Q
c
c       in serial using the BLAS matrix-matrix multiplication routines.
c
c     SOURCE
c
      subroutine smult2(q,r,f,na,nb,nbasis)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF1(f)      common/iofile/iread(42),omax,max
_IF(parallel)
      common/scftim/tdiag(5),tmult,tdum(4)
_ENDIF
INCLUDE(common/timeperiods)
      dimension q(*),r(*),f(*)
c...
_IF1(f)      if(omax) then
_IF1(f)        call maxqhq(q,r,f,na,nb,nbasis)
_IF1(f)        return
_IF1(f)      endif
      nbb=nbasis
_IF(parallel)
      dumtim=dclock()
_ENDIF
      call start_time_period(TP_MULT2)
c...
c...
c...  halve diagonal of f
c...
      m=0
      do loop=1,nbasis
        m=m+loop
        f(m)=f(m)*0.5d0
      enddo
c...
c...  r = q(dagger) * f(scaled)
c...  lower triangle of  f  effectively zero
c...
      call vclr(r,1,nbasis*na)
      call mxmt(q,nbb,1,f,r,na,nbasis)
c...
c...  f = r * q
c...
      call vclr(f,1,na*na)
      call mxmb(r,1,na,q,1,nbb,f,1,na,na,nbasis,na)
c...
c...  r = f + f(dagger)
c...
      m=0
      ll=0
      do loop=1,na
        mm=loop
        limit=loop
        if(loop.eq.na)limit=nb
        do moop=1,limit
          m=m+1
          r(m)=f(ll+moop)+f(mm)
          mm=mm+na
        enddo
        ll=ll+na
      enddo
      call end_time_period(TP_MULT2)
_IF(parallel)
      tmult=tmult+(dclock()-dumtim)
_ENDIF
      return
      end
c     ******
_ENDIF
_IF(meiko)
       subroutine mult2(q,a,h,ncore,ndum,nbasis)
c     this is the pure-scalar version of mult2
c     parallelised over outer ncore loop
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical oipsci
      parameter (maxorb=256)
      dimension a(*),q(*),h(*)
      dimension p(maxorb),y(maxorb)
INCLUDE(common/mapper)
      common/scftim/tdiag(5),tmult,tmult2(4)
      common/bufb/temp(13000)
c...   a=q(transpose) * h * q
c...   a and h stored in triangle form
      dumtim=dclock()
      l2=nbasis*(nbasis+1)/2
      call vclr(a,1,l2)
      iflop=iipsci()
      do 1 j=1,ncore
        if (oipsci()) go to 1
      call dcopy(nbasis,q(ilifq(j)+1),1,y(1),1)
      do 2 i=1,nbasis
      bot=y(i)
      m=iky(i)
      top=h(m+i)*bot
      n=i-1
      if(n)2,2,4
 4    do 5 k=1,n
       gut=h(m+k)
      top=y(k)*gut+top
 5    p(k)=bot*gut+p(k)
 2    p(i)=top
      m=iky(j)
      do 6 i=1,j
      iiii=ilifq(i)
      sum=0.0d0
      do 11 k=1,nbasis
 11   sum=sum+p(k)*q(iiii+k)
  6   a(m+i)=sum
  1   continue
      call gdsum(a,l2,temp)
      iflop=iipsci()
      tmult=tmult+(dclock()-dumtim)
      return
      end
_ENDIF
_IF(cyber205)
      subroutine mult2(q,r,f,nact,mq,nbas)
c..
c..    mult2  r = q(dagger)/f(triangle)/q
c..    nbas = nbasis  / nact = nactiv /  mq not used
c..    r (nbas*nact) and f(triangle/nbas*nbas) must be big enough
c..
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension r(*),f(*),q(*)
      common/scr205/t(2)
       if (nact.ne.mq) call caserr2(' bad call to mult2')
c..    half diagonal of triangle f
      call vxvxts(nbas,0.5,f,iky(2),f,iky(2))
c..    transpose q
      call dagger(nbas,nact,q,nbas,t,nact)
c..    r = q(dagger) * f(lower-trangle)
      call szero(r,nbas*nact)
      call mxmm(t,nact,f,r,nact,nact,nbas)
c..    f = r*q
      call mymo(r,nact,q,1,nbas,f,nact,nact,nbas,nact)
c..    symmetrize f to triangle r
      call symm1(r,f,nact)
c..
      return
      end
c      subroutine mult2(q,r,f,na,nb,nbasis)
cc...
cc...  mult2 r = q(dagger).f.q   (gutil 1985 version)
cc...     ** this version works as well **
cc...
c      dimension q(*),r(*),f(*)
cc..
c      nbb = nbasis
cc...    half f-diagonal
c      m = 0
c      do 1 loop=1,nbasis
c         m = m+loop
c1     f(m) = f(m)*0.5e0
cc...
cc...  r=q(dagger) * f(scaled)
cc...  lower triangle of f effectively zero
cc...
c      call szero(r,nbasis*na)
c      call mxmt(q,nbb,1,f,r,na,nbasis)
cc...
cc...  f = r*q
cc...
c      call szero(f,na*na)
c      call mxmb(r,1,na,q,1,nbb,f,1,na,na,nbasis,na)
cc...
cc...  r = f + f(dagger)
cc...
c      m = 0
c      ll = 0
c      do 4 loop =1,na
c         mm = loop
c         limit = loop
c         if (loop.eq.na) limit = nb
c         do 5 moop=1,limit
c            m = m + 1
c            r(m) = f(ll+moop)+f(mm)
c5        mm = mm + na
c4     ll = ll + na
cc...
c      return
c      end
_ENDIF
      subroutine vtamv (amat,vec,d1,nsa)
      implicit REAL  (a-h,o-z)
c
c     amat = vec(transpose)*amat*vec
c
      dimension amat(nsa*nsa),vec(nsa*nsa),d1(nsa*nsa)
c
      n2 = nsa*nsa
c
      call vclr(d1,1,n2)
      call mxmb(vec,1,nsa,amat,1,nsa,d1,1,nsa,nsa,nsa,nsa)
c
      call vclr(amat,1,n2)
      call mxmb(d1,1,nsa,vec,nsa,1,amat,1,nsa,nsa,nsa,nsa)
c
      return
      end
      subroutine vtamvu(amat,vec,d,ans,nsa)
      implicit REAL  (a-h,o-z)
      dimension amat(nsa*nsa),vec(nsa*nsa),d(nsa*nsa),ans(nsa*nsa)
c
      n2 = nsa*nsa
      call vclr(d,1,n2)
      call mxmb(vec,1,nsa,amat,1,nsa,d,1,nsa,nsa,nsa,nsa)
c
      call mxmb(vec,1,nsa,d,nsa,1,ans,nsa,1,nsa,nsa,nsa)
      return
      end
      subroutine cpuwal(cpu,elapse)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      cpu=cpulft(1)
      call walltime(elapse)
      return
      end
      subroutine timana(mode)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/statis)
c
      cpusec=cpulft(1)
      call walltime(wall)
      timsec(mode)=cpusec-begin +timsec(mode)
      walsec(mode)=wall  -ebegin+walsec(mode)
c
      return
      end
c
c  Time period logging
c
      block data init_time_periods
      implicit none
INCLUDE(common/timeperiods)
      data ttotw/maxtp*0.0d0/
      data ttotc/maxtp*0.0d0/
      data taggc/maxtp*0.0d0/
      data ttotu/maxtp*0.0d0/
      data ttots/maxtp*0.0d0/
      data ntpc/maxtp*0/
      data lab(TP_ENTIRE              )/'entire'/
      data lab(TP_SCF                 )/'scf'/
      data lab(TP_ORFOG               )/'orfog'/
      data lab(TP_DHSTAR              )/'dhstar comp'/
      data lab(TP_RDMAT               )/'rdmat'/
      data lab(TP_DIIS                )/'diis'/
      data lab(TP_DIAG                )/'diag'/
      data lab(TP_APRDMP2             )/'aprdmp2'/
      data lab(TP_DHSTAR_GOP          )/'dhstar gop'/
      data lab(TP_APRM1234            )/'aprm1234'/
      data lab(TP_APRQ34              )/'aprq34'/
      data lab(TP_APRQ1               )/'aprq1'/
      data lab(TP_APRQ2               )/'aprq2'/
      data lab(TP_APRQ2D              )/'aprq2d'/
      data lab(TP_APRQ34D             )/'aprq34d'/
      data lab(TP_APRMP2E             )/'aprmp2e'/
      data lab(TP_MP1PDM              )/'mp1pdm'/
      data lab(TP_MP1PDM_1            )/'mp1pdm 1'/
      data lab(TP_MP1PDM_2            )/'mp1pdm 2'/
      data lab(TP_APR1PDM             )/'apr1pdm'/
      data lab(TP_MP2HESS             )/'mp2hess'/
      data lab(TP_MP2CHF              )/'mp2chf'/
      data lab(TP_MP2MAKEW            )/'mp2makew'/
      data lab(TP_MP2DS               )/'mp2 ds'/
      data lab(TP_MP2BACK_P           )/'mp2back_p'/
      data lab(TP_MP2BACK_F           )/'mp2back_f'/
      data lab(TP_MP2BACKTRAN_2       )/'mp2back_2'/
      data lab(TP_MP2MCDAB            )/'mp2mcdab'/
      data lab(TP_DGOP                )/'dgop'/
      data lab(TP_BCAST               )/'bcast'/
      data lab(TP_NXTVAL              )/'nxtval'/
      data lab(TP_GENRAL              )/'genral(tr)'/
      data lab(TP_GENRAL_1PDM         )/'genral(1pdm)'/
      data lab(TP_GA_PUT_Q2D          )/'apr q2d put'/
      data lab(TP_GA_ACC_Q2D          )/'apr q2d acc'/
      data lab(TP_MKT2AO              )/'mkt2ao'/
      data lab(TP_APRL2               )/'aprl2'/
      data lab(TP_GA_GET_L2           )/'aprl2 get'/
      data lab(TP_APRL34              )/'aprl34'/
      data lab(TP_DGENRL              )/'dgenrl'/
      data lab(TP_MP2                 )/'mp2'/
      data lab(TP_JKDER               )/'jkder'/
      data lab(TP_APRL1234            )/'aprl1234'/
      data lab(TP_APRL1               )/'aprl1'/
      data lab(TP_APRDMP2_I           )/'aprdmp2 init'/
      data lab(TP_JKDER_GET           )/'jkder get'/
      data lab(TP_MP1PDM_3            )/'mp1mdm 3'/
      data lab(TP_MP1PDM_4            )/'mp1mdm 4'/
      data lab(TP_MKT2MO              )/'mkt2mo'/
      data lab(TP_2D_AOINTS           )/'AO integs'/
      data lab(TP_2D_SCF              )/'SCF'/
      data lab(TP_2D_HFGRDN           )/'HF grad'/
      data lab(TP_2D_INDX2T           )/'index2'/
      data lab(TP_2D_MOINTS           )/'MO ints'/
      data lab(TP_2D_TRNFKD           )/'trnfkd'/
      data lab(TP_2D_CHFNDR           )/'CPHF'/
      data lab(TP_2D_QMDER            )/'Quad Mom'/
      data lab(TP_2D_DMDER            )/'Dip Mom'/
      data lab(TP_2D_2D               )/'2dd ints'/
      data lab(TP_2D_CHF              )/'CHF'/
      data lab(TP_2D_NUC              )/'nuclear'/
      data lab(TP_2D_OVL              )/'overlap'/
      data lab(TP_2D_KE               )/'KE'/
      data lab(TP_2D_PE               )/'PE'/
      data lab(TP_2D_2E               )/'2-elec'/
      data lab(TP_2D_TOTAL            )/'hessian'/
      data lab(TP_2D_CHFDRV           )/'chfdrv'/
      data lab(TP_2D_PDENS            )/'pdens'/
      data lab(TP_2D_PFOCK            )/'pfock'/
      data lab(TP_2D_CHFHESS          )/'hessian'/
      data lab(TP_2D_CHFRHS           )/'rhs'/
      data lab(TP_2D_SYMMRHS          )/'symm rhs'/
      data lab(TP_2D_SYMMU            )/'symm soln'/
      data lab(TP_2D_PFOCK_OOOO       )/'oooo'/
      data lab(TP_2D_PFOCK_VOOO       )/'vooo'/
      data lab(TP_2D_PFOCK_VVOO       )/'vvoo'/
      data lab(TP_2D_PFOCK_VOVO       )/'vovo'/
      data lab(TP_2D_PFOCK_SUM        )/'gsum'/
      data lab(TP_2D_AOGEN            )/'AO genral'/
      data lab(TP_2D_AOUT             )/'AO out'/
      data lab(TP_TEST1               )/'test 1'/
      data lab(TP_TEST2               )/'test 2'/
      data lab(TP_TEST3               )/'test 3'/
      data lab(TP_TEST4               )/'test 4'/
      data lab(TP_TEST5               )/'test 5'/
      data lab(TP_TEST6               )/'test 6'/
      data lab(TP_HFGRAD              )/'hfgrad'/
      data lab(TP_GAMULT2             )/'GA mult2'/
      data lab(TP_GAORTHOG            )/'GA orthog'/
      data lab(TP_PDIAG               )/'para diag'/
      data lab(TP_MULT2               )/'mult2'/
      data lab(TP_INTEG               )/'integ'/
      data lab(TP_IOFM1               )/'find (host)'/
      data lab(TP_IOF0                )/'find (node)'/
      data lab(TP_IOF1                )/'find (virt)'/
      data lab(TP_IOF2                )/'find (node0)'/
      data lab(TP_IOF3                )/'find (conc)'/
      data lab(TP_IOF4                )/'find (mem)'/
      data lab(TP_IOF5                )/'find (mem0)'/
      data lab(TP_IOF6                )/'find (ga)'/
      data lab(TP_IOF7                )/'find (ga0)'/
      data lab(TP_IOGM1               )/'get (host)'/
      data lab(TP_IOG0                )/'get (node)'/
      data lab(TP_IOG1                )/'get (virt)'/
      data lab(TP_IOG2                )/'get (node0)'/
      data lab(TP_IOG3                )/'get (conc)'/
      data lab(TP_IOG4                )/'get (mem)'/
      data lab(TP_IOG5                )/'get (mem0)'/
      data lab(TP_IOG6                )/'read (ga)'/
      data lab(TP_IOG7                )/'read (ga0)'/
      data lab(TP_IOPM1               )/'put (host)'/
      data lab(TP_IOP0                )/'put (node)'/
      data lab(TP_IOP1                )/'put (virt)'/
      data lab(TP_IOP2                )/'put (node0)'/
      data lab(TP_IOP3                )/'put (conc)'/
      data lab(TP_IOP4                )/'put (mem)'/
      data lab(TP_IOP5                )/'put (mem0)'/
      data lab(TP_IOP6                )/'wrt (ga)'/
      data lab(TP_IOP7                )/'wrt (ga0)'/
      data lab(TP_IOOM1               )/'open (host)'/
      data lab(TP_IOO0                )/'open (node)'/
      data lab(TP_IOO1                )/'open (virt)'/
      data lab(TP_IOO2                )/'open (node0)'/
      data lab(TP_IOO3                )/'open (conc)'/
      data lab(TP_IOO4                )/'alloc'/
      data lab(TP_IOO5                )/'alloc0'/
      data lab(TP_IOO6                )/'open (ga)'/
      data lab(TP_IOO7                )/'open (ga0)'/
      data lab(TP_IO_GAFILE_READ      )/'ga load'/
      data lab(TP_IO_GAFILE_DUMP      )/'ga dump'/
      data lab(TP_PEIGS               )/'peigs'/
      data lab(TP_PDSYEV              )/'pdsyev'/
      data lab(TP_PDSYEVX             )/'pdsyevx'/
      data lab(TP_PDSYEVD             )/'pdsyevd'/
      data lab(TP_PDSYEVR             )/'pdsyevr'/
      data lab(TP_DFT_JFIT            )/'jfit'/
      data lab(TP_DFT_JFIT_VFORM      )/'form v'/
      data lab(TP_DFT_JFIT_TR         )/'form tr'/
      data lab(TP_DFT_JFIT_NR         )/'form nr'/
      data lab(TP_DFT_JFIT_COEF       )/'fit coef'/
      data lab(TP_DFT_JFIT_KSMAT      )/'form ks'/
      data lab(TP_DFT_JFIT_ENERGY     )/'energy'/
      data lab(TP_DFT_EXQUAD          )/'xc quad'/
      data lab(TP_DFT_EXQUAD_INTRO    )/'intro'/
      data lab(TP_DFT_EXQUAD_INTEG    )/'integ'/
      data lab(TP_DFT_EXQUAD_DGOP     )/'dgop'/
      data lab(TP_DFT_EXQUADF         )/'xc quad f'/
      data lab(TP_DFT_EXQUADF_INTRO   )/'intro'/
      data lab(TP_DFT_EXQUADF_INTEG   )/'integ'/
      data lab(TP_DFT_EXQUADF_DGOP    )/'dgop'/
      data lab(TP_DFT_EXQUADLHS       )/'xc quad lhs'/
      data lab(TP_DFT_EXQUADLHS_INTRO )/'intro'/
      data lab(TP_DFT_EXQUADLHS_INTEG )/'integ'/
      data lab(TP_DFT_EXQUADLHS_DGOP  )/'dgop'/
      data lab(TP_DFT_EXQUADHES       )/'xc quad hes'/
      data lab(TP_DFT_EXQUADHES_INTRO )/'intro'/
      data lab(TP_DFT_EXQUADHES_INTEG )/'integ'/
      data lab(TP_DFT_EXQUADHES_DGOP  )/'dgop'/
      data lab(TP_DFT_EXQUADRHS       )/'xc quad rhs'/
      data lab(TP_DFT_EXQUADRHS_INTRO )/'intro'/
      data lab(TP_DFT_EXQUADRHS_INTEG )/'integ'/
      data lab(TP_DFT_EXQUADRHS_DGOP  )/'dgop'/
      data lab(TP_DFT_EXQUADDKSX      )/'xc quad dksx'/
      data lab(TP_DFT_EXQUADDKSX_INTRO)/'intro'/
      data lab(TP_DFT_EXQUADDKSX_INTEG)/'integ'/
      data lab(TP_DFT_EXQUADDKSX_DGOP )/'dgop'/
      data lab(TP_DFT_EXQUADDKS       )/'xc quad dks'/
      data lab(TP_DFT_EXQUADDKS_INTRO )/'intro'/
      data lab(TP_DFT_EXQUADDKS_INTEG )/'integ'/
      data lab(TP_DFT_EXQUADDKS_DGOP  )/'dgop'/
      data lab(TP_DFT_JMULT           )/'jmult'/
      data lab(TP_DFT_JMULT_INTRO     )/'intro'/
      data lab(TP_DFT_JMULT_SB        )/'sb'/
      data lab(TP_DFT_JMULT_FOCK      )/'fock'/
      data lab(TP_DFT_JMULT_CJAT0     )/'cjat0'/
      data lab(TP_DFT_JFITG           )/'jfit gradient'/
      data lab(TP_DFT_JFITG_VFORM     )/'form v'/
      data lab(TP_DFT_JFITG_TR        )/'form tr'/
      data lab(TP_DFT_JFITG_NR        )/'form nr'/
      data lab(TP_DFT_JFITG_COEF      )/'fit coef'/
      data lab(TP_DFT_JFITG_2C        )/'2c deriv ints'/
      data lab(TP_DFT_JFITG_3C        )/'3c deriv ints'/
      data lab(TP_DFT_JFIT_TR_INIT    )/'fit tr/3c2e'/
      data lab(TP_DFT_JFIT_INV        )/'fit invert'/
      data lab(TP_VB                  )/'vb'/
      data lab(TP_VB_STRUC            )/'vb struct.'/
      data lab(TP_VB_ME               )/'vb matrix el.'/
      data lab(TP_VB_DIAG             )/'vb diag.'/
      data lab(TP_VB_TRAN             )/'vb int.tran.'/
      data lab(TP_VB_VIRT             )/'vb virt.space'/
      data lab(TP_VB_LADM             )/'vb lagr.dens.'/
      data lab(TP_VB_DTRAN            )/'vb dens. tran.'/
      data lab(TP_ISEND               )/'isend'/
      data lab(TP_IRECV               )/'irecv'/
      data lab(TP_WAIT                )/'comm. wait'/
      data lab(TP_STVECP              )/'stvecp'/
      data lab(TP_TVDER               )/'tvder'/
      data lab(TP_SDER                )/'sder'/
      data lab(TP_SGRAD               )/'sgrad'/
      data lab(TP_HELFEY              )/'helfey'/
c
c     F90 timeperiods 
c
      data lab(TP_F90_START           )/'start'/                !206
      data lab(TP_F90_SCF             )/'scf'/
      data lab(TP_F90_BUILD           )/'Fock build'/
      data lab(TP_F90_DIIS            )/'diis'/
      data lab(TP_F90_SIMIL           )/'simil'/
      data lab(TP_F90_DIAG            )/'diag'/
      data lab(TP_F90_BACK            )/'back'/
      data lab(TP_F90_ASSIGN          )/'assign'/
      data lab(TP_F90_ORTHOG          )/'orthog'/
      data lab(TP_F90_MAKE_DENS       )/'make dens'/
      data lab(TP_F90_END             )/'end'/
      data lab(TP_F90_LEV_SHIFT       )/'lev shift'/
      data lab(TP_F90_TESTER_EVAL     )/'tester eval'/
      data lab(TP_F90_DELTA_EVAL      )/'delta eval' /       !219
      data lab(TP_F90_TDOWN           )/'tdown' /            !220
      data lab(TP_F90_RDMAT           )/'rdmat'/             !221
      data lab(TP_F90_INTS            )/'ints'/              !222
c
c     F90 DENSCF timeperiods
c
      data lab(TP_NEWSCF              )/'newscf'/            !223
      data lab(TP_DENSCF              )/'denscf'/            !224
      data lab(TP_DENSCF_BUILD        )/'Fock build'/        !225
      data lab(TP_DENSCF_RDMAT        )/'rdmat'/             !226
      data lab(TP_DENSCF_INTS         )/'ints'/              !227
      data lab(TP_DENSCF_DIAG_S       )/'diag Overlap'/      !228
      data lab(TP_DENSCF_SIMIL        )/'simil'/             !229
      data lab(TP_DENSCF_DIAG         )/'diag'/              !230
      data lab(TP_DENSCF_BACK         )/'back'/              !231
      data lab(TP_DENSCF_MAKE_DENS    )/'make dens '/        !232
      data lab(TP_DENSCF_TDOWN        )/'tdown'/             !233
      data lab(TP_DRHFCL_GA           )/'drhfcl_ga'/         !234
c
c     response theory
c
      data lab(TP_RESPONSE            )/'response'/
      data lab(TP_RPA                 )/'RPA'/
      data lab(TP_TDA                 )/'TDA'/
      data lab(TP_RPANAL              )/'analysis'/
      data lab(TP_RPA_MO2AO           )/'tr mo2ao'/
      data lab(TP_RPA_INT             )/'2e ints'/
      data lab(TP_RPA_CNTRCT          )/'contraction'/
      data lab(TP_RPA_AO2MO           )/'tr ao2mo'/
      data lab(TP_TDA_MO2AO           )/'tr mo2ao'/
      data lab(TP_TDA_INT             )/'2e ints'/
      data lab(TP_TDA_CNTRCT          )/'contraction'/
      data lab(TP_TDA_AO2MO           )/'tr ao2mo'/
c
      data itpdepth/0/
      data itpstack(0)/-1/
      data parent/maxtp*-2/
_IF(vampir)
      data vtclass/maxtp*-2/
      data vtstate/maxtp*-2/
_ENDIF(vampir)
      end

      subroutine reset_time_periods

INCLUDE(common/timeperiods)
INCLUDE(common/iofile)
      logical opg_root
      do i = 1,maxtp
         ttotw(i) = 0.0d0
         ttotc(i) = 0.0d0
         ttots(i) = 0.0d0
         ttotu(i) = 0.0d0
         taggc(i) = 0.0d0
         ntpc(i) = 0
      enddo
c
c     rebuild parent/child relationships
c    
      do i=1,maxtp
         parent(i)=-2
      enddo

      parent(TP_ENTIRE) = -1

      parent(TP_INTEG) = TP_ENTIRE

      parent(TP_SCF) = TP_ENTIRE
      parent(TP_ORFOG) = TP_SCF
c$$$      parent(TP_DHSTAR) = TP_SCF
c$$$      parent(TP_DHSTAR_GOP) = TP_SCF
      parent(TP_RDMAT) = TP_SCF
      parent(TP_DIIS) = TP_SCF
c$$$      parent(TP_DIAG) = TP_SCF
      parent(TP_TEST1) = TP_SCF
      parent(TP_TEST2) = TP_SCF
      parent(TP_TEST3) = TP_SCF
      parent(TP_TEST4) = TP_SCF
      parent(TP_TEST5) = TP_SCF
      parent(TP_TEST6) = TP_SCF

      parent(TP_HFGRAD) = TP_ENTIRE
c
c  mp2 code
c
      parent(TP_MP2) = TP_ENTIRE
      parent(TP_APRDMP2) = TP_MP2
      parent(TP_APRDMP2_I) = TP_APRDMP2
      parent(TP_APRM1234) = TP_APRDMP2
      parent(TP_GENRAL) = TP_APRM1234
      parent(TP_APRQ1) = TP_APRM1234
      parent(TP_APRQ2) = TP_APRM1234
      parent(TP_APRQ2D) = TP_APRM1234
      parent(TP_GA_PUT_Q2D) = TP_APRQ2D
      parent(TP_GA_ACC_Q2D) = TP_APRQ2D
      parent(TP_APRQ34) = TP_APRM1234
      parent(TP_APRQ34D) = TP_APRM1234
      parent(TP_APRMP2E) = TP_APRDMP2
c
c 1pdm
c
      parent(TP_MP1PDM) =TP_MP2
      parent(TP_MP1PDM_1) = TP_MP1PDM
      parent(TP_MP1PDM_2) = TP_MP1PDM
      parent(TP_MP1PDM_3) = TP_MP1PDM
      parent(TP_MP1PDM_4) = TP_MP1PDM

c
c apr1pdm
c
      parent(TP_APR1PDM) =TP_MP2
      parent(TP_MKT2AO) = TP_APR1PDM

      parent(TP_APRL1234) = TP_APR1PDM
      parent(TP_GENRAL_1PDM) = TP_APRL1234
      parent(TP_APRL1) = TP_APRL1234
      parent(TP_APRL2) = TP_APRL1234
      parent(TP_GA_GET_L2) = TP_APRL2
      parent(TP_APRL34) = TP_APRL1234
c gdf:  added tp_mkt2mo 14.3.95
      parent(TP_MKT2MO) = TP_APR1PDM

      parent(TP_MP2HESS) = TP_MP2
      parent(TP_MP2CHF) = TP_MP2
      parent(TP_MP2MAKEW) = TP_MP2
      parent(TP_MP2DS) = TP_MP2

c add jkder
      parent(TP_JKDER) = TP_MP2
      parent(TP_JKDER_GET) = TP_JKDER
      parent(TP_MP2BACK_P) = TP_JKDER
      parent(TP_MP2BACK_F) = TP_JKDER
      parent(TP_MP2BACKTRAN_2) = TP_JKDER
      parent(TP_MP2MCDAB) = TP_JKDER
      parent(TP_DGENRL) = TP_JKDER

c 2nd
      parent(TP_2D_TOTAL) =TP_ENTIRE

      parent(TP_2D_AOINTS) =TP_2D_TOTAL
      parent(TP_2D_AOGEN) = TP_2D_AOINTS
      parent(TP_2D_AOUT) = TP_2D_AOINTS
      parent(TP_2D_SCF) =TP_2D_TOTAL
      parent(TP_2D_HFGRDN) =TP_2D_TOTAL
      parent(TP_2D_INDX2T) =TP_2D_TOTAL
      parent(TP_2D_MOINTS) =TP_2D_TOTAL
      parent(TP_2D_TRNFKD) =TP_2D_TOTAL
      parent(TP_2D_CHFNDR) =TP_2D_TOTAL
      parent(TP_2D_QMDER) =TP_2D_TOTAL
      parent(TP_2D_DMDER) =TP_2D_TOTAL
      parent(TP_2D_2D) =TP_2D_TOTAL

      parent(TP_2D_CHF) =TP_2D_2D
      parent(TP_2D_NUC) =TP_2D_2D
      parent(TP_2D_OVL) =TP_2D_2D
      parent(TP_2D_KE) =TP_2D_2D
      parent(TP_2D_PE) =TP_2D_2D
      parent(TP_2D_2E) =TP_2D_2D

csecd
      parent(TP_2D_CHFHESS) = TP_2D_CHFNDR
      parent(TP_2D_CHFRHS) = TP_2D_CHFNDR
      parent(TP_2D_CHFDRV) = TP_2D_CHFNDR
ccc      parent(TP_MP2CHF) = TP_2D_CHFNDR
      parent(TP_2D_PFOCK) = TP_2D_CHFNDR

      parent(TP_2D_PFOCK_OOOO) = TP_2D_PFOCK
      parent(TP_2D_PFOCK_VOOO) = TP_2D_PFOCK
      parent(TP_2D_PFOCK_VVOO) = TP_2D_PFOCK
      parent(TP_2D_PFOCK_VOVO) = TP_2D_PFOCK
      parent(TP_2D_PFOCK_SUM) = TP_2D_PFOCK

      parent(TP_2D_PDENS) = TP_2D_CHFNDR
      parent(TP_2D_SYMMRHS) = TP_2D_CHFNDR
      parent(TP_2D_SYMMU) = TP_2D_CHFNDR
c
c  CCP1 DFT code timers
c
      parent(TP_DFT_JFIT) = TP_SCF
      parent(TP_DFT_JMULT) = TP_SCF
      parent(TP_DFT_EXQUAD) = TP_SCF
      parent(TP_DFT_EXQUADF) = TP_HFGRAD
      parent(TP_DFT_EXQUADHES) = TP_2D_2D
      parent(TP_DFT_EXQUADLHS) = TP_2D_CHFDRV
      parent(TP_DFT_EXQUADRHS) = TP_2D_CHFRHS
      parent(TP_DFT_EXQUADDKS) = TP_2D_PFOCK
      parent(TP_DFT_EXQUADDKSX) = TP_2D_TOTAL
      parent(TP_DFT_JFITG) = TP_HFGRAD

      parent(TP_DFT_JFIT_VFORM) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_INV) = TP_DFT_JFIT_VFORM

      parent(TP_DFT_JFIT_TR) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_TR_INIT) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_NR) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_COEF) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_KSMAT) = TP_DFT_JFIT
      parent(TP_DFT_JFIT_ENERGY) = TP_DFT_JFIT

      parent(TP_DFT_JFITG_VFORM) = TP_DFT_JFITG
      parent(TP_DFT_JFITG_TR) = TP_DFT_JFITG
      parent(TP_DFT_JFITG_NR) = TP_DFT_JFITG
      parent(TP_DFT_JFITG_COEF) = TP_DFT_JFITG
      parent(TP_DFT_JFITG_2C) = TP_DFT_JFITG
      parent(TP_DFT_JFITG_3C) = TP_DFT_JFITG

      parent(TP_DFT_EXQUAD_INTRO) = TP_DFT_EXQUAD
      parent(TP_DFT_EXQUAD_INTEG) = TP_DFT_EXQUAD
      parent(TP_DFT_EXQUAD_DGOP)  = TP_DFT_EXQUAD

      parent(TP_DFT_EXQUADF_INTRO) = TP_DFT_EXQUADF
      parent(TP_DFT_EXQUADF_INTEG) = TP_DFT_EXQUADF
      parent(TP_DFT_EXQUADF_DGOP)  = TP_DFT_EXQUADF

      parent(TP_DFT_EXQUADHES_INTRO) = TP_DFT_EXQUADHES
      parent(TP_DFT_EXQUADHES_INTEG) = TP_DFT_EXQUADHES
      parent(TP_DFT_EXQUADHES_DGOP)  = TP_DFT_EXQUADHES

      parent(TP_DFT_EXQUADLHS_INTRO) = TP_DFT_EXQUADLHS
      parent(TP_DFT_EXQUADLHS_INTEG) = TP_DFT_EXQUADLHS
      parent(TP_DFT_EXQUADLHS_DGOP)  = TP_DFT_EXQUADLHS

      parent(TP_DFT_EXQUADRHS_INTRO) = TP_DFT_EXQUADRHS
      parent(TP_DFT_EXQUADRHS_INTEG) = TP_DFT_EXQUADRHS
      parent(TP_DFT_EXQUADRHS_DGOP)  = TP_DFT_EXQUADRHS

      parent(TP_DFT_EXQUADDKS_INTRO) = TP_DFT_EXQUADDKS
      parent(TP_DFT_EXQUADDKS_INTEG) = TP_DFT_EXQUADDKS
      parent(TP_DFT_EXQUADDKS_DGOP)  = TP_DFT_EXQUADDKS

      parent(TP_DFT_EXQUADDKSX_INTRO) = TP_DFT_EXQUADDKSX
      parent(TP_DFT_EXQUADDKSX_INTEG) = TP_DFT_EXQUADDKSX
      parent(TP_DFT_EXQUADDKSX_DGOP)  = TP_DFT_EXQUADDKSX

      parent(TP_DFT_JMULT_INTRO) = TP_DFT_JMULT
      parent(TP_DFT_JMULT_SB) = TP_DFT_JMULT
      parent(TP_DFT_JMULT_FOCK) = TP_DFT_JMULT
      parent(TP_DFT_JMULT_CJAT0) = TP_DFT_JMULT_FOCK


c one electron int der timings
      parent(TP_STVECP) = TP_HFGRAD
      parent(TP_TVDER) = TP_HFGRAD
      parent(TP_SDER) = TP_HFGRAD
      parent(TP_SGRAD) = TP_SDER
      parent(TP_HELFEY) = TP_HFGRAD

*IJB newscf_f90 stuff
      parent( TP_NEWSCF ) = TP_SCF

      parent( TP_F90_START ) = TP_NEWSCF
      parent( TP_F90_SCF   ) = TP_NEWSCF
      parent( TP_F90_END   ) = TP_NEWSCF

      parent( TP_F90_BUILD       ) = TP_F90_SCF
      parent( TP_F90_DIIS        ) = TP_F90_SCF
      parent( TP_F90_SIMIL       ) = TP_F90_SCF
      parent( TP_F90_DIAG        ) = TP_F90_SCF
      parent( TP_F90_BACK        ) = TP_F90_SCF
      parent( TP_F90_DIAG        ) = TP_F90_SCF
      parent( TP_F90_ASSIGN      ) = TP_F90_SCF
      parent( TP_F90_ORTHOG      ) = TP_F90_SCF
      parent( TP_F90_MAKE_DENS   ) = TP_F90_SCF
      parent( TP_F90_LEV_SHIFT   ) = TP_F90_SCF
      parent( TP_F90_TESTER_EVAL ) = TP_F90_SCF
      parent( TP_F90_DELTA_EVAL  ) = TP_F90_SCF
      parent( TP_F90_TDOWN       ) = TP_F90_SCF

      parent( TP_F90_RDMAT ) = TP_F90_BUILD
      parent( TP_F90_INTS  ) = TP_F90_BUILD


      parent( TP_DENSCF ) = TP_SCF

      parent( TP_DENSCF_BUILD     ) = TP_DENSCF
      parent( TP_DENSCF_DIAG_S    ) = TP_DENSCF
      parent( TP_DENSCF_SIMIL     ) = TP_DENSCF
      parent( TP_DENSCF_DIAG      ) = TP_DENSCF
      parent( TP_DENSCF_BACK      ) = TP_DENSCF
      parent( TP_DENSCF_MAKE_DENS ) = TP_DENSCF
      parent( TP_DENSCF_TDOWN     ) = TP_DENSCF

      parent( TP_DENSCF_RDMAT ) = TP_DENSCF_BUILD
      parent( TP_DENSCF_INTS  ) = TP_DENSCF_BUILD

      parent( TP_DRHFCL_GA ) = TP_SCF

      parent( TP_RESPONSE   ) = TP_ENTIRE
      parent( TP_RPA        ) = TP_RESPONSE
      parent( TP_TDA        ) = TP_RESPONSE
      parent( TP_RPANAL     ) = TP_RESPONSE
      parent( TP_RPA_MO2AO  ) = TP_RPA
      parent( TP_RPA_INT    ) = TP_RPA
      parent( TP_RPA_CNTRCT ) = TP_RPA
      parent( TP_RPA_AO2MO  ) = TP_RPA
      parent( TP_TDA_MO2AO  ) = TP_TDA
      parent( TP_TDA_INT    ) = TP_TDA
      parent( TP_TDA_CNTRCT ) = TP_TDA
      parent( TP_TDA_AO2MO  ) = TP_TDA
_IF(vampir)
      call vtclassdef("GAMESS-UK",vtclassgam,ierr)
      if (ierr.ne.0) then
         write(iwr,*)"vtclassdef ",i,lab(i),"GAMESS-UK",ierr
         call caserr2("reset_time_periods: vtclassdef failed")
      endif
      do i = 1, maxtp
         call vtclassdef(lab(i),vtclass(i),ierr)
         if (ierr.ne.0) then
            write(iwr,*)"vtclassdef ",i,lab(i),ierr
            call caserr2("reset_time_periods: vtclassdef failed")
         endif
         call vtfuncdef(lab(i),vtclass(i),vtstate(i),ierr)
         if (ierr.ne.0) then
            write(iwr,*)"vtfuncdef ",i,lab(i),ierr
            call caserr2("reset_time_periods: vtfuncdef failed")
         endif
      enddo
_ENDIF(vampir)
      end

      subroutine start_time_period(itab)
      implicit none
INCLUDE(common/timeperiods)
INCLUDE(common/iofile)
      REAL wall, cpusec
      REAL user, sys
      REAL buff(3)
      logical opg_root

      integer itab
      integer ierr

      if(itab.le.0.or.itab.gt.maxtp)then
         write(6,*)itab
         call caserr('bad time period')
      endif

      if(ntpc(itab).lt.0)then
          write(6,*)itab
c         call caserr('same time period started more than once')
      endif
      ntpc(itab) = -ntpc(itab) - 1

_IF(dynamic_tp)
c
c     Update the stack of time periods
c
      itpdepth = itpdepth+1
      if (itpdepth.gt.maxtpdepth) then
         call caserr('time period stack overflow')
      endif
      itpstack(itpdepth) = itab
c
c     Assign parent for this time period. The basic assumption is
c     that every time period has one and only one parent (except 
c     TP_ENTIRE of course which has no parent). If this assumption is
c     violated the only thing we can do is to orphane the time period
c     so that it will be shown on a line of its own. Short of 
c     dynamically mapping out the whole tree this is the best we can do.
c
      if (parent(itab).eq.itpstack(itpdepth-1)) then
c
c        The parent has been initialised already.
c
      else if (parent(itab).eq.-2) then
c
c        The parent has not been initialised yet, so set it.
c
         parent(itab) = itpstack(itpdepth-1)
c
      else
c
c        This time period has more than one parent,
c        recover by orphaning it.
c
         parent(itab) = -1
c
      endif 
_ENDIF(dynamic_tp)
_IF(vampir)
      call vtbegin(vtstate(itab),ierr)
      if (ierr.ne.0) then
         write(iwr,*)"start_time_period: itab,ierr = ",itab,ierr
         call caserr2("start_time_period: vtbegin failed")
      endif
_ENDIF(vampir)
c
c CPU times
c
      call gms_cputime(buff)
      user=buff(1)
      sys=buff(2)
      cpusec=buff(3)
c
c WALL times
c
_IF(ipsc)
      wall = cpusec
_ELSE
      call walltime(wall)
_ENDIF
c
c set start times
      tsc(itab)=cpusec
      tsw(itab)=wall  
      tsu(itab)=user
      tss(itab)=sys

c      write(6,*)'start_time_period ',lab(itab), wall

c      if(opg_root() )then
c         write(6,*)'start_time_period ',lab(itab),wall
c      endif

      return
      end

      subroutine end_time_period(itab)
      implicit none

      integer itab
INCLUDE(common/timeperiods)
INCLUDE(common/iofile)
      REAL  wall, cpusec, user, sys, buff(3)
      logical opg_root
      integer ierr
_IF(dynamic_tp)
c
c     pop this time period off the stack
c
      itpdepth = itpdepth-1
      if (itpdepth.lt.0) then
         call caserr('time period stack underflow')
      endif
_ENDIF(dynamic_tp)
c
c increment counter, store as -ve to indicate
c period in progress
c

      if(itab.le.0.or.itab.gt.maxtp)then
         write(6,*)itab
         call caserr('bad time period')
      endif
c
_IF(vampir)
      call vtend(vtstate(itab),ierr)
      if (ierr.ne.0) then
         write(iwr,*)"end_time_period: itab,ierr = ",itab,ierr
         call caserr2("end_time_period: vtend failed")
      endif
_ENDIF(vampir)
c
c  ntpc(itab) should be negative if we have started
c  this timeperiod
c
      if(ntpc(itab).ge.0)then
         if(opg_root())then
            write(iwr,100)ntpc(itab),itab, lab(itab)
 100        format(1x,'Warning - internal problem with timing data:',
     +           i3,i3,a20,/,'Program will continue, but timing data',
     +           ' may be affected')
         endif
         return
      endif
c
c CPU times
c
      call gms_cputime(buff)
      user=buff(1)
      sys=buff(2)
      cpusec=buff(3)
c
c WALL times
c
_IF(ipsc)
      wall = cpusec
_ELSE
      call walltime(wall)
_ENDIF
c     
c accumulate times


c      if(opg_root())then
c        write(6,*)'end_time_period ',lab(itab), wall- tsw(itab)
c         write(6,*)'end time ',s,nt(itab),
c     &        cpusec - tsc(itab), wall - tsw(itab)
c     &    ,   user - tsu(itab), sys - tss(itab)
c      endif

      ttotc(itab) =  ttotc(itab) + cpusec - tsc(itab)
      ttotw(itab) = ttotw(itab) + wall - tsw(itab)
      ttotu(itab) = ttotu(itab) + user - tsu(itab)
      ttots(itab) = ttots(itab) + sys - tss(itab)

      ntpc(itab) = -ntpc(itab)

      return
      end

      subroutine list_time_periods(all_nodes,aggregate)
      implicit none
INCLUDE(common/timeperiods)
INCLUDE(common/parcntl)
_IF(parallel)
      integer me, nnode, ipg_nodeid, ipg_nnodes
      integer ifrom,len
      integer liw
      integer lenwrd
c     integer itp
c     integer ibuff(maxlablen)
      integer inode
_ENDIF
      integer i
      logical opg_root
      logical all_nodes
      logical aggregate
      REAL  wall, cpusec, user, sys, aggcpu, buff(3)

c
c For loading as library (chemshell, charmm etc)
c
_IF(chemshell,charmm)
      external init_time_periods
_ENDIF

_IF(parallel)
      me = ipg_nodeid()
      nnode = ipg_nnodes()
_ENDIF

_IF(parallel)
      liw = 8 / lenwrd()
_ENDIF

      call gms_cputime(buff)
      user=buff(1)
      sys=buff(2)
      cpusec=buff(3)
      call walltime(wall)
c
c  root node times
c

c 
c generate aggregate times over all nodes
c
      if(aggregate)then
         aggcpu = cpusec
         call pg_dgop(301,aggcpu,1,'+')
         do i = 1,maxtp
            taggc(i) = ttotc(i)
         enddo
         call pg_dgop(301,taggc,maxtp,'+')
      endif

      if(opg_root() .and. iparapr .gt. 0)then
_IF(parallel)
         inode=0
_ENDIF
ccc         write(6,*)'---------------interval timings-----------------'
         do i = 1,maxtp
c
c disable this all info should now be in the timing tree
            if(ntpc(i).gt.0)then
c               write(6,100)lab(i),ttotc(i),ttotw(i)
c     +              ,ttotu(i),ttots(i),ntpc(i)
c 100           format(1x,a10,4f14.3,i7)
            else if(ntpc(i).lt.0)then
c               write(6,101)lab(i),ttotc(i) + cpusec - tsc(i),
c     +              ttotw(i) + wall - tsw(i),
c     +              ttotu(i) + user - tsu(i),
c     +              ttots(i) + sys -  tss(i),
c     +              -ntpc(i)
c 101        format(1x,a10,4f14.3,i7,' *')
           endif
         enddo

         call timetree(0,cpusec,wall,user,sys,aggcpu,aggregate)

      endif
_IF(parallel)

      if(iparapr.eq.2 .and. all_nodes)then
         do inode = 1,nnode - 1
            call pg_synch(1)
            if(me.eq.0)then
c collect and print
c timsec array - comment out those parts that are now hardwired

c               call pg_rcv(101,ntp,liw,len,inode,ifrom,1)
               call pg_rcv(102,ntpc,maxtp*liw,len,inode,ifrom,1)
c               do itp = 1,ntp
c                  call pg_rcv(106+itp,ibuff,maxlablen*liw,
c     &                 len,inode,ifrom,1)
c                  call itoch(ibuff,lab(itp))
c               enddo
               call pg_rcv(104,ttotc,maxtp*8,len,inode,ifrom,1)
               call pg_rcv(105,ttotw,maxtp*8,len,inode,ifrom,1)
               call pg_rcv(106,ttots,maxtp*8,len,inode,ifrom,1)
               call pg_rcv(107,ttotu,maxtp*8,len,inode,ifrom,1)

c        write(6,*)'-----------------node ',inode,' ------------------'
               do i = 1,maxtp

c                  if(ntpc(i).ne.0)
c     +                 write(6,100)lab(i),ttotc(i),ttotw(i)
c     +                          ,ttotu(i),ttots(i),ntpc(i)
               enddo

               call timetree(inode,0.0d0,0.0d0,0.0d0,0.0d0,
     &              0.0d0,.false.)

            else if (me.eq.inode) then
c               call pg_snd(101,maxtp,liw,0,1)
               call pg_snd(102,ntpc,maxtp*liw,0,1)
c               do itp=1,maxtp
c                  call chtoi(ibuff,lab(itp))
c                  call pg_snd(106+itp,ibuff,liw*maxlablen,0,1)
c               enddo
               call pg_snd(104,ttotc,maxtp*8,0,1)
               call pg_snd(105,ttotw,maxtp*8,0,1)
               call pg_snd(106,ttots,maxtp*8,0,1)
               call pg_snd(107,ttotu,maxtp*8,0,1)
            else
c     sit this out
            endif
         enddo
      endif
_ENDIF

cc      call  wrttasks

      return
      end

      subroutine timetree(node, cpu, wall, user, sys, 
     &     aggcpu, aggregate)
      implicit none
      REAL cpu, wall, user, sys, aggcpu
      logical aggregate

INCLUDE(common/timeperiods)
INCLUDE(common/iofile)

      integer node
      integer i, i2, i3, i4, i5, i6
      integer iflag2, iflag3, iflag4, iflag5, iflag6

      REAL timec, timeac, timew, timeu, times
      REAL timec2, timeac2, timew2, timeu2, times2

      REAL sumc2, sumac2, sumw2, sumu2, sums2
      REAL sumc3, sumac3, sumw3, sumu3, sums3
      REAL sumc4, sumac4, sumw4, sumu4, sums4
      REAL sumc5, sumac5, sumw5, sumu5, sums5
      REAL sumc6, sumac6, sumw6, sumu6, sums6

      character*20 lother

      lother = 'other'

      write(iwr,198)node
 198  format(1x,'Timing analysis for node ',i3,/,
     &       1x,'============================')
c
c   sub-period analysis 
c
 199  format(1x,
     &     15x,37x,"        CPU                 Aggr.    WALL  eff.")
 200  format(1x,a15,15x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
 2001 format(1x,a15,15x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,i3,
     &     f9.2,1x,' *')
 201  format(1x,3x,a15,12x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
 2011 format(1x,3x,a15,12x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,i3,
     &     f9.2,1x,' *')
 202  format(1x,6x,a15,9x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
 203  format(1x,9x,a15,6x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
 204  format(1x,12x,a15,3x,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
 205  format(1x,15x,a15,i7,f9.2,' (',f9.2,' u+',f9.2,' s)',f12.2,
     &     f9.2,1x,i3)
      write(iwr,199)

      do i = 1,maxtp     
c
c establish unconnected intervals
            if(parent(i) .gt. 0 .and. ntpc(parent(i)).eq. 0)
     &           parent(i) = -2
c
c list hierarchical times
            if(parent(i) .eq. -1.and.ntpc(i) .ne. 0)then
               if(ntpc(i).gt.0)then
c
c     complete
                  timec = ttotc(i)
                  if(aggregate)then
                     timeac = taggc(i)
                  else
                     timeac = 0.0d0
                  endif
                  timeu = ttotu(i)
                  times = ttots(i)
                  timew = ttotw(i)
                  if(dabs(timew).gt.0.0d0) then
                    write(iwr,200)lab(i),ntpc(i),timec,timeu,times,
     +                    timeac,timew,nint(100.0d0*timec/timew)
                  endif
               else
c
c     incomplete
                  timec = ttotc(i) + cpu - tsc(i)
                  timeac = 0.0
                  timeu = ttotu(i) + user - tsu(i)
                  times = ttots(i) + sys - tss(i)
                  timew = ttotw(i) + wall - tsw(i)

                  if(dabs(timew).gt.0.0d0) then
                    write(iwr,2001)lab(i),-ntpc(i),timec,timeu,times,
     +                    timeac,timew,nint(100.0d0*timec/timew)
                  endif
               endif

               sumc2 = 0.0d0
               sumac2 = 0.0d0
               sumw2 = 0.0d0
               sumu2 = 0.0d0
               sums2 = 0.0d0
               iflag2 = 0
               do i2 = 1,maxtp     
                  if(parent(i2) .eq. i.and.ntpc(i2) .ne. 0)then

                     if(ntpc(i2).gt.0)then
c
c     complete
                        timec2 = ttotc(i2)
                        timeac2 = taggc(i2)
                        timeu2 = ttotu(i2)
                        times2 = ttots(i2)
                        timew2 = ttotw(i2)
                        if(dabs(timew2).gt.0.0d0) then
                          write(iwr,201)lab(i2),ntpc(i2),timec2,
     +                    timeu2,times2,timeac2,timew2,
     &                    nint(100.0d0*timec2/timew2)
                        endif

                     else
c
c     incomplete
                        timec2 = ttotc(i2) + cpu - tsc(i2)
                        timeac2 = 0.0d0
                        timeu2 = ttotu(i2) + user - tsu(i2)
                        times2 = ttots(i2) + sys - tss(i2)
                        timew2 = ttotw(i2) + wall - tsw(i2)

                        if(dabs(timew2).gt.0.0d0) then
                          write(iwr,2011)lab(i2),-ntpc(i2),timec2,
     +                    timeu2,times2,timeac2,timew2,
     &                    nint(100.0d0*timec2/timew2)
                        endif

                     endif

                     iflag2 = 1
                     sumc2 = sumc2 + timec2
                     sumac2 = sumac2 + timeac2
                     sumw2 = sumw2 + timew2
                     sumu2 = sumu2 + timeu2
                     sums2 = sums2 + times2


               sumc3 = 0.0d0
               sumac3 = 0.0d0
               sumw3 = 0.0d0
               sumu3 = 0.0d0
               sums3 = 0.0d0
               iflag3 = 0
               do i3 = 1,maxtp     
                  if(parent(i3) .eq. i2.and.ntpc(i3) .ne. 0)then
                     iflag3 = 1
                     sumc3 = sumc3 + ttotc(i3)
                     sumac3 = sumac3 + taggc(i3)
                     sumw3 = sumw3 + ttotw(i3)
                     sumu3 = sumu3 + ttotu(i3)
                     sums3 = sums3 + ttots(i3)
                     write(iwr,202)lab(i3),ntpc(i3),
     &                    ttotc(i3),ttotu(i3),ttots(i3),
     &                    taggc(i3),
     &                    ttotw(i3),
     &            nint(100.0d0*ttotc(i3)/(max(1.0d-10,ttotw(i3))))
               sumc4 = 0.0d0
               sumac4 = 0.0d0
               sumw4 = 0.0d0
               sumu4 = 0.0d0
               sums4 = 0.0d0
               iflag4 = 0
               do i4 = 1,maxtp     
                  if(parent(i4) .eq. i3.and.ntpc(i4) .ne. 0)then
                     iflag4 = 1
                     sumc4 = sumc4 + ttotc(i4)
                     sumac4 = sumac4 + taggc(i4)
                     sumw4 = sumw4 + ttotw(i4)
                     sumu4 = sumu4 + ttotu(i4)
                     sums4 = sums4 + ttots(i4)
                     if(dabs(ttotw(i4)).gt.0.0d0) then
                       write(iwr,203)lab(i4),ntpc(i4),
     &                 ttotc(i4),ttotu(i4),ttots(i4),taggc(i4),
     &                 ttotw(i4),nint(100.0d0*ttotc(i4)/ttotw(i4))
                     endif
               sumc5 = 0.0d0
               sumac5 = 0.0d0
               sumw5 = 0.0d0
               sumu5 = 0.0d0
               sums5 = 0.0d0
               iflag5 = 0
               do i5 = 1,maxtp     
                  if(parent(i5) .eq. i4.and.ntpc(i5) .ne. 0)then
                     iflag5 = 1
                     sumc5 = sumc5 + ttotc(i5)
                     sumac5 = sumac5 + taggc(i5)
                     sumw5 = sumw5 + ttotw(i5)
                     sumu5 = sumu5 + ttotu(i5)
                     sums5 = sums5 + ttots(i5)
                     if(dabs(ttotw(i5)).gt.0.0d0) then
                       write(iwr,204)lab(i5),ntpc(i5),
     &                 ttotc(i5),ttotu(i5),ttots(i5),taggc(i5),
     &                ttotw(i5),nint(100.0d0*ttotc(i5)/ttotw(i5))
                     endif
               sumc6 = 0.0d0
               sumac6 = 0.0d0
               sumw6 = 0.0d0
               sumu6 = 0.0d0
               sums6 = 0.0d0
               iflag6 = 0
               do i6 = 1,maxtp     
                  if(parent(i6) .eq. i5.and.ntpc(i6) .ne. 0)then
                     iflag6 = 1
                     sumc6 = sumc6 + ttotc(i6)
                     sumac6 = sumac6 + taggc(i6)
                     sumw6 = sumw6 + ttotw(i6)
                     sumu6 = sumu6 + ttotu(i6)
                     sums6 = sums6 + ttots(i6)
                     if(dabs(ttotw(i6)).gt.0.0d0) then
                       write(iwr,205)lab(i6),ntpc(i6),
     &                 ttotc(i6),ttotu(i6),ttots(i6),taggc(i6),
     &                  ttotw(i6),nint(100.0d0*ttotc(i6)/ttotw(i6))
                     endif
                  endif
               enddo
               if(iflag6.ne.0.and.dabs(ttotw(i5)-sumw6).gt.0.0d0)
     &              write(iwr,205)lother,0,ttotc(i5) - sumc6, 
     &              ttotu(i5) - sumu6,
     &              ttots(i5) - sums6, 
     &              taggc(i5) - sumac6,
     &              ttotw(i5) - sumw6,
     &          nint(100.0d0*(ttotc(i5)-sumc6)/(ttotw(i5)-sumw6))

                  endif
               enddo
               if(iflag5.ne.0.and.dabs(ttotw(i4)-sumw5).gt.0.0d0)
     &              write(iwr,204)lother,0,ttotc(i4) - sumc5,  
     &              ttotu(i4) - sumu5, 
     &              ttots(i4) - sums5,
     &              taggc(i4) - sumac5,
     &              ttotw(i4) - sumw5,
     &          nint(100.0d0*(ttotc(i4)-sumc5)/(ttotw(i4)-sumw5))
                  endif
               enddo
               if(iflag4.ne.0.and.dabs(ttotw(i3)-sumw4).gt.0.0d0)
     &              write(iwr,203)lother,0,ttotc(i3) - sumc4, 
     &              ttotu(i3) - sumu4, 
     &              ttots(i3) - sums4,
     &              taggc(i3) - sumac4,
     &              ttotw(i3) - sumw4,
     &          nint(100.0d0*(ttotc(i3)-sumc4)/(ttotw(i3)-sumw4))
                  endif
               enddo
               if(iflag3.ne.0.and.dabs(timew2-sumw3).gt.0.0d0)
     &              write(iwr,202)lother,0,timec2 - sumc3,
     &              timeu2 - sumu3,
     &              times2 - sums3,
     &              timeac2 - sumac3,
     &              timew2 - sumw3,
     &          nint(100.0d0*(timec2-sumc3)/(timew2-sumw3))
                  endif
               enddo
               if(iflag2.ne.0.and.dabs(timew-sumw2).gt.0.0d0)
     &              write(iwr,201)lother,0,timec - sumc2, 
     &              timeu - sumu2, 
     &              times - sums2,
     &              timeac - sumac2,
     &              timew - sumw2,
     &          nint(100.0d0*(timec-sumc2)/(timew-sumw2))
            endif
         enddo
c
c orphan time periods
         do i = 1,maxtp
            if(ntpc(i).ne.0.and.parent(i).eq.-2)then
               timec = ttotc(i)
               timeac = taggc(i)
               timeu = ttotu(i)
               times = ttots(i)
               timew = ttotw(i)
               if(dabs(timew).gt.0.0d0) then
                  write(iwr,200)lab(i),ntpc(i),timec,timeu,times,
     &                 timeac,timew,nint(100.0d0*timec/timew)
               endif
               
            endif
         enddo
c
_IF(ga)
_IFN(chemshell,charmm)
      call ga_print_stats
_ENDIF
_ENDIF
c
      end
      subroutine timout
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *7 ztext,ztext1,ztext2
INCLUDE(common/statis)
INCLUDE(common/iofile)
INCLUDE(common/parcntl)
INCLUDE(common/timeperiods)
_IF(parallel)
      common/scftim/tdiag,tmxmb,torfog,tfock,
     +    tdiis,tmult,tspac(4)
_ENDIF
      dimension ztext(200)
      dimension ztext1(100),ztext2(100)
      equivalence (ztext1(1),ztext(1)),(ztext2(1),ztext(101))
      data ztext1/
     + 'input  ','       ','       ','       ',
     + 'vector ','generat','ion    ','       ',
     + '1e-inte','gral ev','aluatio','n      ',
     + '2e-inte','gral ev','aluatio','n      ',
     + 'scf    ','       ','       ','       ',
     + '1e-inte','gral de','rivativ','es     ',
     + '2e-inte','gral de','rivativ','es     ',
     + 'symmetr','y adapt','ion    ','       ',
     + 'casscf ','       ','       ','       ',
     + 'mc-scf ','       ','       ','       ',
     + 'integra','l trans','formati','on     ',
     + 'direct-','ci     ','       ','       ',
     + 'table-c','i      ','       ','       ',
     + 'wave-fu','nction ','analysi','s      ',
     + 'greens ','functio','n i.p. ','       ',
     + 't.d.a. ','i.p.   ','       ','       ',
     + 'utility',' functi','on s   ','       ',
     + 'full-ci','       ','       ','       ',
     + 'model e','nergy +',' force ','       ',
     + 'model m','inimisa','tion   ','       ',
     + 'model d','ynamics','       ','       ',
     + 'coupled',' hartre','e-fock ','eqtns. ',
     + '1e-inte','gral 2n','d deriv','atives ',
     + '2e-inte','gral 2n','d deriv','atives ',
     + 'mp2 ene','rgy    ','       ','       ' /
      data ztext2/
     + 'rpa    ','       ','       ','       ',
     + 'mopac  ','       ','       ','       ',
     + 'intrins','ic reac','tion co','ord.   ',
     + 'ccsd   ','       ','       ','       ',
     + 'respons','e      ','       ','       ',
     + 'ExCorr.',' energy','       ','       ',
     + 'ExCorr.',' gradie','nt     ','       ',
     + 68*'       ',
     + 'other  ','       ','       ','       '/
_IF(parallel)
      dimension tbuff(2,51)
_ENDIF
c
      begin=cpulft(1)
      call walltime(ebegin)
      top=0.0d0
      etop=0.0d0
      do 4 i=1,49
         etop=etop+walsec(i)
         top=top+timsec(i)
 4    continue
      timsec(50)=begin-top-tstart
      walsec(50)=ebegin-etop-estart
      tempt=begin-tstart
      tempe=ebegin-estart
      if(tempt.gt.0.0d0.and.tempe.gt.0.0d0) then
         t1=100.0d0/tempt
         e1=100.0d0/tempe
         write(iwr,1)
         jj=0
         do 2 i=1,50
            if(timsec(i).gt.1.0d-3.or.walsec(i).gt.1.0d-3)then
               top=timsec(i)*t1
               etop=walsec(i)*e1
               write(iwr,3)(ztext(jj+k),k=1,4),timsec(i),top,
     &              walsec(i),etop
 3             format(1x,4a7,1x,f13.2,3x,f7.2,2x,f14.2,3x,f7.2)
            endif
            jj=jj+4
 2       continue
         write(iwr,5)tempt,tempe
 5       format(1x,78('*')/
     +        1x,'Total',24x,f13.2,12x,f14.2/1x,78('*'))
      else
         write(6,*)'times.. skip',tempt,tempe
      endif
_IF(parallel)
c
c  stats for nodes other than 0
      if(iparapr.eq.2)then
         mynode = ipg_nodeid()
         nnode = ipg_nnodes()
         do inode = 1,nnode - 1

            call pg_synch(1)

            if(mynode.eq.0)then
c     collect and print
c     timsec array
               write(6,*)'****** node ',inode,' *******'
               call pg_rcv(100,tbuff,102*8,len,inode,ifrom,1)
               do i=1,50
                  timsec(i)=tbuff(1,i)
                  walsec(i)=tbuff(2,i)
               enddo
               tempt= tbuff(1,51)
               tempe= tbuff(2,51)
               if(tempt.gt.0.0d0.and.tempe.gt.0.0d0) then
                  t1=100.0d0/tempt
                  e1=100.0d0/tempe
c     write(iwr,1)
                  jj=0
                  do 22 i=1,50
                     if(timsec(i).gt.1.0d-3.or.walsec(i).gt.1.0d-3)then
                        top=timsec(i)*t1
                        etop=walsec(i)*e1
                        write(iwr,3)(ztext(jj+k),k=1,4),timsec(i),top,
     &                       walsec(i),etop
                     endif
                     jj=jj+4
 22               continue
                  write(iwr,5)tempt,tempe
               endif
            else if (mynode.eq.inode) then
               do i=1,50
                  tbuff(1,i)=timsec(i)
                  tbuff(2,i)=walsec(i)
               enddo
               tbuff(1,51) = tempt
               tbuff(2,51) = tempe
               call pg_snd(100,tbuff,102*8,0,1)
            else
c     sit this out
            endif
         enddo
      endif
c
c now IPSC time data
c
      write(iwr,100) tdiag,torfog,tfock,
     +              tdiis,tmult
100   format(/
     + 1x,'matrix operations (seconds)'/
     + 1x,'==========================='/
     + 1x,'diagonalisation (jacobi) ',1x,f8.2/
     + 1x,'orthogonalisation (orfog)',1x,f8.2/
     + 1x,'fock-build (dbuild)      ',1x,f8.2/
     + 1x,'DIIS (mult2 etc)         ',1x,f8.2/
     + 1x,'q*hq (mult2)             ',1x,f8.2/
     + 1x,'==========================='//)
_ENDIF
c
c  now the time_period statistics..
c
      call end_time_period(TP_ENTIRE)
_IF(parallel,timings)
      call list_time_periods(.true., .true.)
_ENDIF
      return
 1       format(//
     *        1x,78('*')/
     *        1x,'gamess timing analysis'/
     *1x,'task',25x,'cpu (seconds)',3x,'percent',2x,'wall (seconds)',
     *    3x,'percent'/
     *    1x,78('*'))

      end
_IF()
c
c - old version - to check .f..
c 
      subroutine timoutxx
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *7 ztext,ztext1,ztext2
INCLUDE(common/statis)
INCLUDE(common/iofile)
_IF1(h)      common/scftim/tdiag,tmxmb,torfog,tfock,
_IF1(h)     +              tdiis,tmult,tspac(4)
      dimension ztext(200)
      dimension ztext1(100),ztext2(100)
      equivalence (ztext1(1),ztext(1)),(ztext2(1),ztext(101))
      data ztext1/
     *'input  ','       ','       ','       ','vector ','generat',
     *'ion    ','       ','1e-inte','gral ev','aluatio','n      ',
     *'2e-inte','gral ev','aluatio','n      ','scf    ','       ',
     *'       ','       ','1e-inte','gral de','rivativ','es     ',
     *'2e-inte','gral de','rivativ','es     ','symmetr','y adapt',
     *'ion    ','       ','casscf ','       ','       ','       ',
     *'mc-scf ','       ','       ','       ','integra','l trans',
     *'formati','on     ','direct-','ci     ','       ','       ',
     *'table-c','i      ','       ','       ','wave-fu','nction ',
     *'analysi','s      ','greens ','functio','n i.p. ','       ',
     *'t.d.a. ','i.p.   ','       ','       ','utility',' functi',
     *'on s   ','       ','full-ci','       ','       ','       ',
     *'model e','nergy +',' force ','       ','model m','inimisa',
     *'tion   ','       ','model d','ynamics','       ','       ',
     *'coupled',' hartre','e-fock ','eqtns. ','1e-inte','gral 2n',
     3'd deriv','atives ','2e-inte','gral 2n','d deriv','atives ',
     *'mp2 ene','rgy    ','       ','       '/
      data ztext2/
     *'rpa    ','       ','       ','       ','mopac  ','       ',
     *'       ','       ','intrins','ic reac','tion co','ord.   ',
     *'ccsd   ','       ','       ','       ','respons','e      ',
     *'       ','       ',
     *76*'       ',
     *'other  ','       ','       ','       '/
c
      begin=cpulft(1)
      call walltime(ebegin)
      top=0.0d0
      etop=0.0d0
      do 4 i=1,49
      etop=etop+walsec(i)
 4    top=top+timsec(i)
      timsec(50)=begin-top-tstart
      walsec(50)=ebegin-etop-estart
      tempt=begin-tstart
      tempe=ebegin-estart
      if(tempt.gt.0.0d0.and.tempe.gt.0.0d0) then
      t1=100.0d0/tempt
      e1=100.0d0/tempe
c
      write(iwr,1)
 1    format(//
     *1x,78('*')/
     *1x,'gamess timing analysis'/
     *1x,'task',25x,'cpu (seconds)',3x,'percent',2x,'wall (seconds)',
     *3x,'percent'/
     *1x,78('*'))
      jj=0
      do 2 i=1,50
      if(timsec(i).lt.1.0d-3.and.walsec(i).lt.1.0d-3)go to 2
      top=timsec(i)*t1
      etop=walsec(i)*e1
      write(iwr,3)(ztext(jj+k),k=1,4),timsec(i),top,walsec(i),etop
 3    format(1x,4a7,1x,f13.2,3x,f7.2,2x,f14.2,3x,f7.2)
 2    jj=jj+4
      write(iwr,5)tempt,tempe
 5    format(1x,78('*')/
     +       1x,'Total',24x,f13.2,12x,f14.2/1x,78('*'))
      endif
_IF1(h)      write(iwr,100) tdiag,torfog,tfock,
_IF1(h)     +              tdiis,tmult
_IF1(h)100   format(/
_IF1(h)     + 1x,'matrix operations (seconds)'/
_IF1(h)     + 1x,'==========================='/
_IF1(h)     + 1x,'diagonalisation (jacobi) ',1x,f8.2/
_IF1(h)     + 1x,'orthogonalisation (orfog)',1x,f8.2/
_IF1(h)     + 1x,'fock-build (dbuild)      ',1x,f8.2/
_IF1(h)     + 1x,'DIIS (mult2 etc)         ',1x,f8.2/
_IF1(h)     + 1x,'q*hq (mult2)             ',1x,f8.2/
_IF1(h)     + 1x,'==========================='//)
      return
      end
_ENDIF
      subroutine timrem(tleft)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/timez)
INCLUDE(common/iofile)
      t1=cpulft(1)
      tleft=timlim-t1
      call flushn(iwr)
      return
      end
      subroutine flushn(iunit)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
c
c...  flush sequential fortran stream on unix machines
c
INCLUDE(common/work)
_IFN(ibm,cray,ipsc,rs6000,hp700,hpux11,vax)
_IF(altix,linux)
c...  just flush output (flush is something else)
      if (iunit.eq.6.and.oflush) call flushout
_ELSE
      if (oflush) then
         call flush(iunit)
      endif
_ENDIF
_ENDIF
      return
      end
      subroutine filchk(nfile,iblk,lblk,notape,ilab)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/files)
INCLUDE(common/discc)
      common/disc/isel(5),ipos(maxlfn)
      dimension iblk(*),lblk(*),notape(*),ztext(40),itext(10)
      data xdash/'-'/
      data itext/28,18,28,28,17,22,25,28,14,9/
      data ztext/
     *'2-electr','on ao in','tegral f','iles    ',
     *'secondar','y mainfi','le      ','        ',
     *'2-electr','on mo in','tegral f','iles    ',
     *'2-partic','le ao de','nsity ma','trix    ',
     *'loop for','mula tap','e       ','        ',
     *'reordere','d formul','a tape  ','        ',
     *'symmetry',' adapted',' mainfil','e       ',
     *'2-partic','le mo de','nsity ma','trix    ',
     *'direct-c','i file  ','        ','        ',
     *'Fock fil','e       ','        ','        '/
c
c ----- print and check file specifications
c
      ii=(ilab-1)*4
      ifield=itext(ilab)
      do 1 i=1,nfile
      yjjj=yed(notape(i))
      if(minbl(notape(i)).ne.0)iblk(i)=minbl(notape(i))
      call filec(yjjj,iblk(i),iunit,irep)
      if(irep.ge.0) go to 2
      if(irep+1)8,9,9
c
c      ----- invalid ddname ---
c
 9    call caserr2('invalid ddname parameter')
c
c      ----- invalid starting block ----
c
 8    write(iwr,10)(ztext(ii+j),j=1,4),yjjj,iblk(i)
 10   format(//1x,4a8//
     *' *** invalid starting block specified on ',a4,
     *' : ',i6)
      call caserr2('illegal search of data set')
 2    if(irep.eq.0) go to 7
      isel(1)=iunit
      write(iwr,3)(ztext(ii+j),j=1,4),yjjj
 3    format(//1x,4a8,
     *' : section ',a4,' not allocated')
_IF1(iv)      call fault
      call caserr2('data set not assigned')
 7    if(maxbl(iunit).gt.0)lblk(i)=maxbl(iunit)
 1    continue
      write(iwr,4)(ztext(ii+j),j=1,4)
      write(iwr,5)(xdash,i=1,ifield)
 5    format(1x,40a1)
 4    format(//1x,4a8)
      call filprn(nfile,iblk,lblk,notape)
      return
      end
_IFN(ibm)
      subroutine filec(yname,iblk,iunit,irep)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/discc)
c
      maxb=99999
      irep=0
      do  1 iunit=1,maxlfn
      if(yname.eq.yed(iunit))go to 2
 1    continue
      irep=-1
      return
 2    if(iblk.lt.1.or.iblk.gt.maxb)then
      irep=-2
      endif
      return
      end
_ENDIF
      subroutine getmat(s,t,f,p,q,r,charge,nbasis,
     *ostf,isec1e)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension s(*),t(*),f(*),p(*),q(*),r(*),charge(*)
      dimension ostf(6)
      common/blkin/potnuc(4),pint(508)
INCLUDE(common/dump3)
      data m2/2/
      call secinf(idum,numdu,idum,idum)
c
      call secget(isec1e,m2,iblk34)
      lentri=nbasis*(nbasis+1)/2
      lenb = lensec(lentri)
c
      call rdedx(potnuc, 511, iblk34, numdu)
      mword = 4 + 2/lenwrd()
      call dcopy(mword,potnuc,1,charge,1)
      ibl3s = iblk34 + 1
c
      if (ostf(1)) then
        call rdedx(s,lentri,ibl3s,numdu)
      endif
      ibl3t = ibl3s + lenb
c
      if (ostf(2)) then
        call rdedx(t,lentri,ibl3t,numdu)
      endif
      ibl3f = ibl3t + lenb
c
      if (ostf(3)) then
        call rdedx(f,lentri,ibl3f,numdu)
      endif
      ibl3x = ibl3f + lenb
c
      if(ostf(4)) then
        call rdedx(p,lentri,ibl3x,numdu)
        do i = 1,lentri
         p(i) = -p(i)
        enddo
      endif
      ibl3y = ibl3x + lenb
c
      if(ostf(5)) then
        call rdedx(q,lentri,ibl3y,numdu)
        do i = 1,lentri
         q(i) = -q(i)
        enddo
      endif
      ibl3z = ibl3y + lenb
c
      if(ostf(6)) then
        call rdedx(r,lentri,ibl3z,numdu)
        do i = 1,lentri
         r(i) = -r(i)
        enddo
      endif
c
      return
      end

      subroutine upack3(packed,i,j,k)
      implicit none
      integer i,j,k
_IF(i8)
      integer packed
_IF(oldpack3)
      i=shiftr(packed,32)
      j=shiftr(packed,16).and.177777b
      k=packed.and.177777b
_ELSE
      i=ishft(packed,-40)
      j=iand(ishft(packed,-20),z'FFFFF')
      k=iand(packed,z'FFFFF')
_ENDIF
_ELSEIF(GIGA_DUMP)
c     use integer *8 quantity for packing
c     and unpacking in c
      integer * 8 packed
c
      call upack3cc(packed,i,j,k)
_ELSE
      REAL packed
      call upack3cc(packed,i,j,k)
_ENDIF
      return
      end


_IF(bits8)
      subroutine pack3(i,j,k,packed)
      implicit none
      integer i,j,k
      REAL packed
      call pack3cc(packed,i,j,k)
_ELSE
      function pack3(i,j,k)
      implicit none
      integer i,j,k
      REAL pack3
_IF(GIGA_DUMP)
c     use integer *8 quantity for packing and unpacking in c
      integer * 8 packed
      REAL dummy
      equivalence  (dummy,packed)
      call pack3cc(packed,i,j,k)
      pack3=dummy
_ELSEIF(i8)
      integer packed
      REAL dummy
      equivalence  (dummy,packed)
      packed=ior(ior(ishft(i,40),ishft(j,20)),k)
      pack3=dummy
_ELSE
      call pack3cc(pack3,i,j,k)
_ENDIF
_ENDIF
      return
      end

_IFN(upck-equiv)
_IF1(iv)      subroutine sgmat(l2,nshell,iky,b,c,p)
_IFN1(iv)      subroutine sgmat(l2,nshell,b,c,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
      dimension b(*),c(*),p(*)
_IF1(iv)      dimension iky(*)
      common/blkin/gg(510),mword
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IFN1(iv)      common/craypk/intin(1360)
_IFN(ibm,vax)
INCLUDE(common/mapper)
_ENDIF
_IF1(iv)      call upak8v(gg(num2e+1),i205)
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
_IFN1(iv)      call unpack(gg(num2e+1),lab816,intin,numlab)
_ELSE
      subroutine sgmat(l2,nshell,b,c,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      logical *1 intin,iii,jjj,kkk,lll
      dimension b(*),c(*),p(*)
_IF(i8)
      dimension iii(8),jjj(8),kkk(8),lll(8)
_ELSE
      dimension iii(4),jjj(4),kkk(4),lll(4)
_ENDIF
      common/blkin/gg(340),intin(1360),mword
INCLUDE(common/mapper)
      equivalence (iii(1),i),(jjj(1),j),(kkk(1),k),
     +  (lll(1),l)
      data i,j,k,l/4*0/
_ENDIF
_IFN1(iv)      int4=1
      if(nshell.gt.1) then
      do 6 iw=1,mword
_IF(upck-equiv)
_IF(littleendian)
      jjj(1) = intin(int4 )
      iii(1) = intin(int4+1)
      lll(1) = intin(int4+2)
      kkk(1) = intin(int4+3)
_ELSEIF(i8)
      iii(8) =intin(int4  )
      jjj(8) =intin(int4+1)
      kkk(8) =intin(int4+2)
      lll(8) =intin(int4+3)
_ELSE
      iii(4) =intin(int4  )
      jjj(4) =intin(int4+1)
      kkk(4) =intin(int4+2)
      lll(4) =intin(int4+3)
_ENDIF
_ELSEIF(ibm,vax)
      i=i205(iw)
      j=j205(iw)
      k=k205(iw)
      l=l205(iw)
_ELSEIF(littleendian)
      j=intin(int4  )
      i=intin(int4+1)
      l=intin(int4+2)
      k=intin(int4+3)
_ELSE
      i=intin(int4  )
      j=intin(int4+1)
      k=intin(int4+2)
      l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
_IF1(x)c$dir scalar
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
_IF1(iv)    6 continue
_IFN1(iv)    6 int4=int4+4
      else
      do 66 iw=1,mword
_IF(upck-equiv)
_IF(littleendian)
      jjj(1) = intin(int4 )
      iii(1) = intin(int4+1)
      lll(1) = intin(int4+2)
      kkk(1) = intin(int4+3)
_ELSEIF(i8)
      iii(8) =intin(int4  )
      jjj(8) =intin(int4+1)
      kkk(8) =intin(int4+2)
      lll(8) =intin(int4+3)
_ELSE
      iii(4) =intin(int4  )
      jjj(4) =intin(int4+1)
      kkk(4) =intin(int4+2)
      lll(4) =intin(int4+3)
_ENDIF
_ELSEIF(ibm,vax)
      i=i205(iw)
      j=j205(iw)
      k=k205(iw)
      l=l205(iw)
_ELSEIF(littleendian)
      j=intin(int4  )
      i=intin(int4+1)
      l=intin(int4+2)
      k=intin(int4+3)
_ELSE
      i=intin(int4  )
      j=intin(int4+1)
      k=intin(int4+2)
      l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
_IF1(iv)   66 continue
_IFN1(iv)   66 int4=int4+4
      endif
      return
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      end
_IFN(ibm,vax)
      subroutine sgmat_255(l2,nshell,b,c,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
      dimension b(*),c(*),p(*)
      common/blkin/gg(510),mword
      common/craypk/intin(1360)
INCLUDE(common/mapper)
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
      call unpack(gg(num2e+1),lab816,intin,numlab)
      int4=1
      if(nshell.gt.1) then
      do 6 iw=1,mword
_IF(littleendian)
       j=intin(int4  )
       i=intin(int4+1)
       l=intin(int4+2)
       k=intin(int4+3)
_ELSE
       i=intin(int4  )
       j=intin(int4+1)
       k=intin(int4+2)
       l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
_IF1(x)c$dir scalar
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
    6 int4=int4+4
      else
      do 66 iw=1,mword
_IF(littleendian)
       j=intin(int4  )
       i=intin(int4+1)
       l=intin(int4+2)
       k=intin(int4+3)
_ELSE
       i=intin(int4  )
       j=intin(int4+1)
       k=intin(int4+2)
       l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
   66 int4=int4+4
      endif
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end
_ENDIF
_IFN(ibm)
      subroutine jksupe(hscm0,hscm1,hscm2,dscm0,dscm1,nscmf,l2)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension hscm0(*),hscm1(*),hscm2(*),dscm0(*),dscm1(*)
      common/blkin/gout(510),nint
INCLUDE(common/atmblk)
_IF1(v)      common/craypk/ij205(204),kl205(204)
_IFN1(v)      common/craypk/integ(1360)
c
_IF1(v)      call upak6v(gout(num2ejk+num2ejk+1),ij205)
_IFN1(v)      call unpack(gout(num2ejk+num2ejk+1),lab1632,
_IFN1(v)     +            integ,numlabjk)
c
c---------both closed + open shells------------
c
_IFN1(v)      iword=1
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
      do   1 m = 1,nint
_IF1(v)      ij = ij205(m)
_IF1(v)      kl = kl205(m)
_IFN1(v)      ij = integ(iword  )
_IFN1(v)      kl = integ(iword+1)
      goutpm = gout(m)
      goutkm = gout(num2ejk+m)
      hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
      hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
_IF1(x)c$dir scalar
      do 2   i = 1,nscmf
      hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
      hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
      hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
      hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
      ij = ij + l2
    2 kl = kl + l2
_IF1(v)    1 continue
_IFN1(v)    1 iword=iword+2
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end
      subroutine gsup(h,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IF(cray,ksr)
      common/craypk/integ(1360)
      common/blkin /g(510),mword
      dimension h(*),p(*)
      call unpack(g(num2ep+1),lab1632,integ,numlabp)
      iword=1
      do 1 iw=1,mword
      ij=integ(iword)
      kl=integ(iword+1)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 iword=iword+2
_ELSEIF(vax)
      integer * 2 iij,kkl
      common/blkin/g(340),iij(340),kkl(340),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=iij(iw)
      kl=kkl(iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
_ELSE
      integer *2 integ
      common/blkin /g(340),integ(2,340),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=integ(1,iw)
      kl=integ(2,iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
_ENDIF
      return
      end
      subroutine gsup_255(h,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IF(cray,ksr)
      common/craypk/integ(1360)
      common/blkin /g(510),mword
      dimension h(*),p(*)
      call unpack(g(num2ep+1),lab1632,integ,numlabp)
      iword=1
      do 1 iw=1,mword
      ij=integ(iword)
      kl=integ(iword+1)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 iword=iword+2
_ELSEIF(vax)
      common/blkin/g(255),iij(255),kkl(255),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=iij(iw)
      kl=kkl(iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
_ELSE
_IF(i8)
      integer *4 integ
_ENDIF
      common/blkin /g(255),integ(2,255),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=integ(1,iw)
      kl=integ(2,iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
_ENDIF
      return
      end
_IF(oldpack)
_IF1(i)@process directive('*vdir:')
      subroutine gsupp(h,p)
      implicit REAL  (a-h,o-z)
_IF1(atpgbdrhs)      integer i8temp(2),fiv11
_IF1(x)      integer *8 fiv11,i8temp
_IF1(k)      integer fiv11,i8temp
_IF1(xpgbdrhsk)      integer *2 ijkl
_IF1(iv)      integer fiv11
_IF1(iv)      integer *2 ij205,kl205
_IF1(iv)      dimension i8temp(2)
      dimension h(*),p(*)
_IF1(xpgbdrhsk)      common/blkin/g(340),ijkl(2,340),mword
_IFN1(ivxpgbdrhsk)      common/blkin/g(340),gij(170),mword
_IFN1(ivxpgbdrhsk)      common/craypk/ij205(340),kl205(340)
_IF1(iv)      common/blkin/g(340),ij205(340),kl205(340),mword
_IFN1(ivxpgbdrhsk)      dimension iij(170),iballs(680)
_IF1(iv)      equivalence (r8temp,i8temp(1))
_IFN1(ivxpgbdrhsk)      equivalence (iij(1),gij(1))
_IFN1(civ)      equivalence (i8temp,r8temp)
_IFN1(c)      data fiv11/511/
      if(mword.le.0) return
c
_IF1(cat)c
_IF1(cat)c need kl205 with unit stride to use spdot/spaxpy
_IF1(cat)      call unpack(iij,16,iballs,680)
_IF1(cat)      iword = 1
_IF1(cat)      do 19 i = 1,mword
_IF1(cat)          ij205(i) = iballs(iword)
_IF1(cat)          kl205(i) = iballs(iword+1)
_IF1(cat)19        iword = iword + 2
c
      iwlo = 1
 10   continue
_IF1(c)      nij = and(511,g(iwlo))
_IFN1(c)      r8temp=g(iwlo)
_IF1(x)      nij=iand(fiv11,i8temp)
_IF1(k)      nij=and(fiv11,i8temp)
_IF1(bivt)      nij=iand(fiv11,i8temp(2))
_IF1(ad)      nij=iand(fiv11,i8temp(1))
_IF1(pgsr)      nij=and(fiv11,i8temp(2))
_IF1(h)      nij=and(fiv11,i8temp(1))
_IFN1(xpgbdrhsk)      ij = ij205(iwlo)
_IF1(xpgbdrhsk)      ij = ijkl(1,iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword) call caserr2('invalid length in gsupp')
_IFN1(xpgbdrhsk)      if(ij.ne.ij205(iwhi)) call caserr2('ij mismatch in gsupp')
_IF1(xpgbdrhsk)      if(ij.ne.ijkl(1,iwhi)) call caserr2('ij mismatch in gsupp')
_IF1(c)c
_IF1(c)c cray cal calls ... no faster than cft on xmp. need on 1s.
_IF1(c)c
_IF1(c)      h(ij) = h(ij) + spdot(nij,p,kl205(iwlo),g(iwlo))
_IF1(c)      call spaxpy(nij,p(ij),g(iwlo),h,kl205(iwlo))
c
_IFN1(c)c
_IFN1(c)c fortran source
_IFN1(c)c
_IFN1(c)       temp=0.0d0
_IFN1(c)       pij = p(ij)
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
_IFN1(c)       do 20 iw = iwlo,iwhi
_IFN1(cxpgbdrhsk)           kl = kl205(iw)
_IF1(xpgbdrhsk)           kl = ijkl(2,iw)
_IFN1(c)           h(kl) = h(kl) + pij * g(iw)
_IFN1(c) 20        temp = temp + p(kl) * g(iw)
_IFN1(c)       h(ij) = h(ij) + temp
c
      iwlo = iwhi + 1
      if(iwlo.le.mword) goto 10
      return
      end
_ELSE
_IF1(i)@process directive('*vdir:')
      subroutine gsupp(h,p)
      implicit REAL  (a-h,o-z)
_IF(convex)
      integer *8 fiv11,i8temp
      integer *2 ijkl
_ELSEIF(ksr,i8)
      integer fiv11,i8temp
      integer *2 ijkl
_ELSEIF(ibm,vax)
      integer fiv11
      integer *2 ij205,kl205
      dimension i8temp(2)
_ELSEIF(cray)
_ELSE
      integer i8temp(2),fiv11
      integer *2 ijkl
_ENDIF
      dimension h(*),p(*)
_IF(cray,titan,alliant)
      common/blkin/g(340),gij(170),mword
      common/craypk/ij205(340),kl205(340)
      dimension iij(170),iballs(680)
      equivalence (iij(1),gij(1))
_ELSEIF(ibm,vax)
      common/blkin/g(340),ij205(340),kl205(340),mword
      equivalence (r8temp,i8temp(1))
_ELSE
      common/blkin/g(340),ijkl(2,340),mword
_ENDIF

_IFN1(civ)      equivalence (i8temp,r8temp)
_IFN1(c)      data fiv11/511/
      if(mword.le.0) return
c
_IF(cray,alliant,titan)
c
c need kl205 with unit stride to use spdot/spaxpy
      call unpack(iij,16,iballs,680)
      iword = 1
      do 19 i = 1,mword
          ij205(i) = iballs(iword)
          kl205(i) = iballs(iword+1)
19        iword = iword + 2
_ENDIF
c
      iwlo = 1
 10   continue
_IF(cray)
      nij = and(511,g(iwlo))
_ELSE
      r8temp=g(iwlo)
_IF(convex,ksr,i8)
      nij=IAND64(fiv11,i8temp)
_ELSEIF(littleendian)
      nij=IAND32(fiv11,i8temp(1))
_ELSE
      nij=IAND32(fiv11,i8temp(2))
_ENDIF
_ENDIF
_IF(cray,alliant,titan,ibm)
      ij = ij205(iwlo)
_ELSE
      ij = ijkl(1,iwlo)
_ENDIF
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword) call caserr2('invalid length in gsupp')

_IF(cray,alliant,titan,ibm)
      if(ij.ne.ij205(iwhi)) call caserr2('ij mismatch in gsupp')
_ELSE
      if(ij.ne.ijkl(1,iwhi)) call caserr2('ij mismatch in gsupp')
_ENDIF
_IF(cray)
c
c cray cal calls ... no faster than cft on xmp. need on 1s.
c
      h(ij) = h(ij) + spdot(nij,p,kl205(iwlo),g(iwlo))
      call spaxpy(nij,p(ij),g(iwlo),h,kl205(iwlo))
c
_ELSE
c
c fortran source
c
       temp=0.0d0
       pij = p(ij)
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
       do 20 iw = iwlo,iwhi
_IF(alliant,titan,ibm)
           kl = kl205(iw)
_ELSE
           kl = ijkl(2,iw)
_ENDIF
           h(kl) = h(kl) + pij * g(iw)
 20        temp = temp + p(kl) * g(iw)
       h(ij) = h(ij) + temp
c
_ENDIF
      iwlo = iwhi + 1
      if(iwlo.le.mword) goto 10
      return
      end
_ENDIF
_ENDIF
_IFN(ibm)
      subroutine gsupa(h,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IF1(v)      integer *2 iij,kkl
_IFN1(v)      common/blkin /g(510),mword
_IF1(v)      common/blkin/g(408),iij(204),kkl(204),mword
_IFN1(v)      common/craypk/integ(1360)
      dimension h(*),p(*)
_IFN1(v)      call unpack(g(num2ejk+num2ejk+1),lab1632,integ,numlabjk)
_IFN1(v)      iword=1
      do 1 iw=1,mword
_IFN1(v)      ij=integ(iword)
_IFN1(v)      kl=integ(iword+1)
_IF1(v)      ij=iij(iw)
_IF1(v)      kl=kkl(iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
_IFN1(v)    1 iword=iword+2
_IF1(v)    1 continue
      return
      end
_ENDIF
      function seccpu()
      implicit REAL  (a-h,o-z)
INCLUDE(common/iofile)
_IFN1(f)      seccpu=cpulft(1)
_IF1(f)      icode=sys$gettime(ak,wall)
_IF1(f)      seccpu=ak
      call flushn(iwr)
      return
      end

_IF(bits8)
c     As a subroutine
      subroutine pack2(nw,m,i4)
      implicit none
      integer nw,m
      integer*4 i4(2)
_IF(littleendian)
      i4(2)=nw
      i4(1)=m
_ELSE
      i4(1)=nw
      i4(2)=m
_ENDIF
      return
      end
_ELSE
c     As a function
      function pack2(nw,m)
      implicit none
      integer nw,m
      integer*4 i4(2)
_IF(i8drct)
      integer*8 pack, pack2
_ELSE
      REAL pack, pack2      
_ENDIF(i8drct)
      equivalence (i4(1),pack)

_IF(littleendian)
      i4(2)=nw
      i4(1)=m
_ELSE
      i4(1)=nw
      i4(2)=m
_ENDIF
      pack2 = pack
      return
      end
_ENDIF(bits8)


_IF(bits8)
      subroutine upack2(i4,left,iright)
      implicit none
      integer left,iright
      integer*4 i4(2)
_IF(littleendian)
      left=i4(2)
      iright=i4(1)
_ELSE
      left=i4(1)
      iright=i4(2)
_ENDIF
      return
      end
_ELSE
      subroutine upack2(b,left,iright)
      implicit none
      integer left,iright
      integer*4 i4(2)
_IF(i8drct)
      integer*8 b,unpack
_ELSE
      REAL b,unpack
_ENDIF(i8drct)      
      equivalence (i4(1),unpack)
      unpack=b
_IF(littleendian)
      left=i4(2)
      iright=i4(1)
_ELSE
      left=i4(1)
      iright=i4(2)
_ENDIF
      return
      end
_ENDIF

_IF1(ivf)       subroutine upak8v(ix,intij)
_IF1(ivf)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)      implicit character *8 (z),character *1 (x)
_IF1(ivf)      implicit character *4 (y)

_IF1(i)      logical * 1 mij,mkl
_IF1(i)      logical*  1 i205,j205,k205,l205
_IF1(i)      dimension ix(*),intij(*)
_IF1(i)      common/blkin/gout(340),mij(680),mkl(680),nint
_IF1(i)      common/craypk/i205(1360),j205(1360),k205(1360),l205(1360)
_IF1(i)      i4=4
_IF1(i)      do 1 i=1,nint
_IF1(i)      il=i+i
_IF1(i)      i205(i4)=mij(il-1)
_IF1(i)      j205(i4)=mij(il  )
_IF1(i)      k205(i4)=mkl(il-1)
_IF1(i)      l205(i4)=mkl(il  )
_IF1(i)    1 i4=i4+4

_IF1(v)      logical * 1 mij,mkl
_IF1(v)      dimension ix(*),intij(*)
_IF1(v)      common/blkin/gout(340),mij(680),mkl(680),nint
_IF1(v)      common/craypk/i205(340),j205(340),k205(340),l205(340)
_IF1(v)      do 1 i=1,nint
_IF1(v)      il=i+i
_IF1(v)      j205(i)=mij(il-1)
_IF1(v)      i205(i)=mij(il  )
_IF1(v)      l205(i)=mkl(il-1)
_IF1(v) 1    k205(i)=mkl(il  )
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      ii=1
_IF1(f)      do 1 i=1,85
_IF1(f)      intij(ii  )=extract(ix(i),0,8)
_IF1(f)      intij(ii+1)=extract(ix(i),16,8)
_IF1(f)      intij(ii+2)=extract(ix(i),32,8)
_IF1(f)      intij(ii+3)=extract(ix(i),48,8)
_IF1(f)      jj=ii+340
_IF1(f)      intij(jj  )=extract(ix(i),8,8)
_IF1(f)      intij(jj+1)=extract(ix(i),24,8)
_IF1(f)      intij(jj+2)=extract(ix(i),40,8)
_IF1(f)      intij(jj+3)=extract(ix(i),56,8)
_IF1(f)      kk=jj+340
_IF1(f)      j=i+85
_IF1(f)      intij(kk  )=extract(ix(j),0,8)
_IF1(f)      intij(kk+1)=extract(ix(j),16,8)
_IF1(f)      intij(kk+2)=extract(ix(j),32,8)
_IF1(f)      intij(kk+3)=extract(ix(j),48,8)
_IF1(f)      ll=kk+340
_IF1(f)      intij(ll  )=extract(ix(j),8,8)
_IF1(f)      intij(ll+1)=extract(ix(j),24,8)
_IF1(f)      intij(ll+2)=extract(ix(j),40,8)
_IF1(f)      intij(ll+3)=extract(ix(j),56,8)
_IF1(f)   1  ii=ii+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(ivf)       subroutine pak4v(intij,ix)
_IF1(ivf)       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)       implicit character *8 (z),character *1 (x)
_IF1(ivf)       implicit character *4 (y)
_IF1(iv)      integer*2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      j=341
_IF1(iv)      do 1 i=1,340
_IF1(iv)      ix(i  )= intij(i)
_IF1(iv)      ix(j  )= intij(j)
_IF1(iv)    1 j=j+1
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      int=1
_IF1(f)      do 1 i=1,170
_IF1(f)      ix(i)=insert(
_IF1(f)     *             insert(intij(int+3),32,16,intij(int+2)),0,32,
_IF1(f)     *             insert(intij(int+1),32,16,intij(int)) )
_IF1(f)   1  int=int+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(iv)      subroutine upak4v(ix,intij)
_IF1(iv)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(iv)      implicit character *8 (z),character *1 (x)
_IF1(iv)      implicit character *4 (y)
_IF1(iv)      integer * 2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      j=341
_IF1(iv)      do 1 i=1,340
_IF1(iv)      intij(i)=ix(i)
_IF1(iv)      intij(j)=ix(j)
_IF1(iv)    1 j=j+1
_IF1(iv)      return
_IF1(iv)      end
_IF1(ivf)       subroutine pak6v(intij,ix)
_IF1(ivf)       implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(ivf)       implicit character *8 (z),character *1 (x)
_IF1(ivf)       implicit character *4 (y)
_IF1(iv)      integer*2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      common/blkin/gout(510),mword
_IF1(iv)      j=205
_IF1(iv)      do 1 i=1,mword
_IF1(iv)      ix(i  )= intij(i)
_IF1(iv)      ix(j  )= intij(j)
_IF1(iv)    1 j=j+1
_IF1(f)      word ix
_IF1(f)      dimension ix(*),intij(*)
_IF1(f)      int=1
_IF1(f)      do 1 i=1,102
_IF1(f)      ix(i)=insert(
_IF1(f)     *             insert(intij(int+3),32,16,intij(int+2)),0,32,
_IF1(f)     *             insert(intij(int+1),32,16,intij(int)) )
_IF1(f)    1 int=int+4
_IF1(ivf)      return
_IF1(ivf)      end
_IF1(iv)      subroutine upak6v(ix,intij)
_IF1(iv)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(iv)      implicit character *8 (z),character *1 (x)
_IF1(iv)      implicit character *4 (y)
_IF1(iv)      integer * 2 ix
_IF1(iv)      dimension intij(*),ix(*)
_IF1(iv)      common/blkin/gout(510),mword
_IF1(iv)      j=205
_IF1(iv)      do 1 i=1,mword
_IF1(iv)      intij(i)=ix(i)
_IF1(iv)      intij(j)=ix(j)
_IF1(iv)    1 j=j+1
_IF1(iv)      return
_IF1(iv)      end
_IFN(cray,ibm,vax)
      subroutine sgmata(fock,p)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(upck-equiv)
      logical *1 integ,iii,jjj,kkk,lll
      dimension p(*),fock(*)
_IF(i8)
      dimension iii(8),jjj(8),kkk(8),lll(8)
_ELSE
      dimension iii(4),jjj(4),kkk(4),lll(4)
_ENDIF
      common/blkin/gg(340),integ(1360),mword
      equivalence (iii(1),i),(jjj(1),j),(kkk(1),k),
     +            (lll(1),l)
      data i,j,k,l/4*0/
_ELSE
      dimension p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
c     this was set a while back by HVD in response to
c     apparent failure under Linux. I can find no trace
c     of this (if required, izero should be invoked in 
c     the level above sgmata ...)
c
c     call izero(1360,integ,1)
      call unpack(gg(num2e+1),lab816,integ,numlab)
_ENDIF
c
      iword = 1
      do 6000 iw=1,mword
_IF(upck-equiv)
_IFN(littleendian)
_IF(i8)
      iii(8) = integ(iword )
      jjj(8) = integ(iword+1)
      kkk(8) = integ(iword+2)
      lll(8) = integ(iword+3)
_ELSE
      iii(4) = integ(iword )
      jjj(4) = integ(iword+1)
      kkk(4) = integ(iword+2)
      lll(4) = integ(iword+3)
_ENDIF
_ELSE
      jjj(1) = integ(iword )
      iii(1) = integ(iword+1)
      lll(1) = integ(iword+2)
      kkk(1) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
      j = integ(iword )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
 6000 iword=iword+4
      return
      end
      subroutine sgmata_255(fock,p)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
      dimension p(*),fock(*)
_IFN(ipsc,t3d)
      integer *2 integ
      common/blkin/gg(255),integ(1020),mword
_ELSE
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
_ENDIF
c
_IF(t3d,ipsc)
      call unpack(gg(num2e+1),lab816,integ,numlab)
_ENDIF
_IF(ccpdft)
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
c
      iword = 1
      do 6000 iw=1,mword
_IF(ipsc,littleendian)
      j = integ(iword )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex

        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
_ELSE
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
_ENDIF
 6000 iword=iword+4
      return
      end
_IF(ccpdft)
      subroutine sgmata_dft(fock,p,facex,ocoul,oexch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      logical oexch, ocoul
_IF(upck-equiv)
      logical *1 integ,iii,jjj,kkk,lll
      dimension p(*),fock(*)
_IF(i8)
      dimension iii(8),jjj(8),kkk(8),lll(8)
_ELSE
      dimension iii(4),jjj(4),kkk(4),lll(4)
_ENDIF
      common/blkin/gg(340),integ(1360),mword
      equivalence (iii(1),i),(jjj(1),j),(kkk(1),k),
     +            (lll(1),l)
      data i,j,k,l/4*0/
_ELSE
      dimension p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
_ENDIF
c
      iword = 1
      do 6000 iw=1,mword
_IF(upck-equiv)
_IFN(littleendian)
_IF(i8)
      iii(8) = integ(iword )
      jjj(8) = integ(iword+1)
      kkk(8) = integ(iword+2)
      lll(8) = integ(iword+3)
_ELSE
      iii(4) = integ(iword )
      jjj(4) = integ(iword+1)
      kkk(4) = integ(iword+2)
      lll(4) = integ(iword+3)
_ENDIF
_ELSE
      jjj(1) = integ(iword )
      iii(1) = integ(iword+1)
      lll(1) = integ(iword+2)
      kkk(1) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
      j = integ(iword )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      gik=gg(iw)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
c
c coulomb term
c
      if (ocoul) then
        ij=ikyi+j
        kl=ikyk+l
        g2=gik+gik
        g4=g2+g2
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
      if (oexch) then
c
c DFT exchange with scaling factor
c
        ik=ikyi+k
        il=ikyi+l
        jk=ikyj+k
        jl=ikyj+l
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
 6000 iword=iword+4
      return
      end
_ENDIF
_ENDIF
_IF(cray)
_IF(ccpdft)
      subroutine sgmata_dft(fock,p,facex,ocoul,oexch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      logical oexch, ocoul
      dimension p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
cdir$ list
cdir$ novector
      iword = 1
      do 6000 iw=1,mword
      i = integ(iword )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
      gik=gg(iw)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
c
c coulomb term
c
      if (ocoul) then
        ij=ikyi+j
        kl=ikyk+l
        g2=gik+gik
        g4=g2+g2
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
      if(oexch) then
c
c DFT exchange with scaling factor
c
        ik=ikyi+k
        il=ikyi+l
        jk=ikyj+k
        jl=ikyj+l
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
 6000 iword=iword+4
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end
_ENDIF
      subroutine sgmata_255(fock,p)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
      dimension p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
_IF(ccpdft)
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
c
cdir$ list
cdir$ novector
      iword = 1
      do 6000 iw=1,mword
      i = integ(iword )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
c
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex

        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
_ELSE
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
_ENDIF
 6000 iword=iword+4
cdir$ vector
cdir$ nolist
      return
      end
_ENDIF
_IF(ibm)
      subroutine sgmata(fock,p)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
      integer * 2 i,j,k,l
      logical *1 iii(2),jjj(2),kkk(2),lll(2),mij,mkl
      dimension p(*),fock(*)
      common/blkin/g(340),mij(680),mkl(680),mword
c
      equivalence(iii(1),i),(jjj(1),j),(kkk(1),k),(lll(1),l)
      data i,j,k,l/4*0/
c
_IF(ccpdft)
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
      do 6000 iw=1,mword
      il=iw+iw
      iii(2)=mij(il-1)
      jjj(2)=mij(il)
      kkk(2)=mkl(il-1)
      lll(2)=mkl(il)
      gik=g(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
c
c coulomb term
c
      if (ocoul) then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
 
      else if(oexch) then
c
c DFT exchange with scaling factor
c
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
_ELSE
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
_ENDIF
 6000 continue
      return
      end
_ENDIF
_IFN(cray)
_IF(ccpdft)
      subroutine proc2f(fock,p,ak,q,
     +                 facex,ocoul,oexch)
_ELSE
      subroutine proc2f(fock,p,ak,q)
_ENDIF
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(ccpdft)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
_ENDIF
_IF(upck-equiv)
      logical *1 integ,iii,jjj,kkk,lll
      dimension ak(*),q(*),p(*),fock(*)
_IF(i8)
      dimension iii(8),jjj(8),kkk(8),lll(8)
_ELSE
      dimension iii(4),jjj(4),kkk(4),lll(4)
_ENDIF
      common/blkin/gg(340),integ(1360),mword
      equivalence (iii(1),i),(jjj(1),j),(kkk(1),k),
     + (lll(1),l)
      data i,j,k,l/4*0/
c
_ELSE
      dimension ak(*),q(*),p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
_ENDIF
      iword=1
      do 6000 iw=1,mword
_IF(upck-equiv)
_IFN(littleendian)
_IF(i8)
      iii(8) = integ(iword  )
      jjj(8) = integ(iword+1)
      kkk(8) = integ(iword+2)
      lll(8) = integ(iword+3)
_ELSE
      iii(4) = integ(iword  )
      jjj(4) = integ(iword+1)
      kkk(4) = integ(iword+2)
      lll(4) = integ(iword+3)
_ENDIF
_ELSE
      jjj(1) = integ(iword  )
      iii(1) = integ(iword+1)
      lll(1) = integ(iword+2)
      kkk(1) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
_ELSE
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
_ENDIF
 6000 iword=iword+4
      return
      end
_IF(ccpdft)
      subroutine proc2_255(fock,p,ak,q,
     +                 facex,ocoul,oexch)
_ELSE
      subroutine proc2_255(fock,p,ak,q)
_ENDIF
      implicit REAL  (a-h,o-z)
      dimension ak(*),q(*),p(*),fock(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(ccpdft)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
_ENDIF
_IFN(ipsc,t3d)
      integer *2 integ
      common/blkin/gg(255),integ(1020),mword
_ELSE
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
_ENDIF
c
_IF(ipsc,t3d)
      call unpack(gg(num2e+1),lab816,integ,numlab)
_ENDIF
      iword=1
      do 6000 iw=1,mword
_IF(ipsc,littleendian)
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
_ELSE
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
_ENDIF
 6000 iword=iword+4
      return
      end
_ENDIF
_IF(cray)
_IF(ccpdft)
      subroutine proc2f(fock,p,ak,q,facex,ocoul,oexch)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
      dimension ak(*),q(*),p(*),fock(*)
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
cdir$ list
cdir$ novector
      iword=1
      do 6000 iw=1,mword
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      gik=gg(iw)
c
c coulomb term
c
      if (ocoul) then
        ij=ikyi+j
        kl=ikyk+l
        g2=gik+gik
        g4=g2+g2
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if (oexch) then
c
c full term for HF or weighted term for b3lyp etc
c
        ik=ikyi+k
        il=ikyi+l
        jk=ikyj+k
        jl=ikyj+l
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
 6000 iword=iword+4
cdir$ vector
cdir$ nolist
      return
      end
_ENDIF
_IF(ccpdft)
      subroutine proc2_255(fock,p,ak,q,
     +                 facex,ocoul,oexch)
_ELSE
      subroutine proc2_255(fock,p,ak,q)
_ENDIF
      implicit REAL  (a-h,o-z)
      dimension ak(*),q(*),p(*),fock(*)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IF(ccpdft)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
_ENDIF
      common/blkin/gg(510),mword
      common/craypk/integ(1360)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
cdir$ list
cdir$ novector
      iword=1
      do 6000 iw=1,mword
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
_ELSE
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
_ENDIF
 6000 iword=iword+4
cdir$ vector
cdir$ nolist
      return
      end
_ENDIF
_IF(alliant)
      recursive subroutine sgmtap(fock,p,l2,gg,mij,integ,mword)
cvd$r noconcur
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension p(l2),fock(l2),gg(340),integ(1360)
      byte mij(1360)
c
cvd$  vector
      do 10 loop=1,mword*4
         integ(loop)=mij(loop)
 10   continue
c
      iword = 1
      do 6000 iw=1,mword
      j = integ(iword )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
 6000 iword=iword+4
      return
      end
_ENDIF
_IF1()      subroutine triad(n,scalar,a,b)
_IF1()      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      dimension a(*),b(*)
_IF1()      if(n.gt.0)then
_IF1()      do 1 loop=1,n
_IF1()   1  a(loop)=a(loop)+scalar*b(loop)
_IF1()      endif
_IF1()      return
_IF1()      end
_IF1(vatpgbdrsk)      subroutine getcor(lword,mfree)
_IF1(vatpgbdrsk)      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1(vatpgbdrsk)      implicit character *8 (z),character *1 (x)
_IF1(vatpgbdrsk)      implicit character *4 (y)
_IF1(v)      common/vcore/core(500000)
_IF1(atpgbdrsk)      if(lword.eq.0) lword=500000
_IF1(v)      lword=500000
_IF1(vatpgbdrsk)      mfree=40000
_IF1(vatpgbdrsk)      return
_IF1(vatpgbdrsk)      end
_IF1()c       OLD convex memory allocation
_IF1()      subroutine getcor(lword,mfree)
_IF1()      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      common/vcore/icore(1)
_IF1()      if(lword.eq.0) lword=500000
_IF1()      itpr=malloc((lword+1000)*8)
_IF1()c ... above addition of 1000 words to cure apparent
_IF1()c ... corruption of top of allocated core (guga-ci/casscf)
_IF1()      if ( itpr .eq. 0 ) then
_IF1()       write(6,*) ' insufficient space on the runtime stack '
_IF1()       stop
_IF1()      endif
_IF1()      icore(1)=itpr
_IF1()      mfree=40000
_IF1()      return
_IF1()      end
_IF1(x)      subroutine getcor(lword,mfree)
_IF1(x)      implicit REAL (a-h,p-w),integer (i-n),logical  (o)
_IF1(x)      implicit character *8 (z),character *1 (x)
_IF1(x)      implicit character *4 (y)
_IF1(x)      common/vcore/icore(1)
_IF1(x)      logical first
_IF1(x)      data first /.true./
_IF1(x)      save first
_IF1(x)      ioldptr=icore(1)
_IF1(x)      if (first) goto 1
_IF1(x)      call dalloc (ioldptr, ier)
_IF1(x)      if( ier .lt. 0 ) then
_IF1(x)      write(6,*) ' invalid memory address for deallocation'
_IF1(x)      stop
_IF1(x)      endif
_IF1(x) 1    if(lword.eq.0) lword=500000
_IF1(x)      ll = (lword+1000)*8
_IF1(x)      call nalloc (ll ,iptr, ier)
_IF1(x)c ... above addition of 1000 words to cure apparent
_IF1(x)c ... corruption of top of allocated core (guga-ci/casscf)
_IF1(x)      if ( ier .lt. 0 ) then
_IF1(x)       write(6,*) ' insufficient memory available for allocation '
_IF1(x)       stop
_IF1(x)      endif
_IF1(x)      icore(1)=iptr
_IF1(x)      mfree=40000
_IF1(x)      first = .false.
_IF1(x)      return
_IF1(x)      end
_IF1()      subroutine clock(ztime)
_IF1()      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
_IF1()      implicit character *8 (z),character *1 (x)
_IF1()      implicit character *4 (y)
_IF1()      data zblank/' '/
_IF1()      ztime=zblank
_IF1()      return
_IF1()      entry date(zdate)
_IF1()      zdate=zblank
_IF1()      return
_IF1()      entry tmleft(isec)
_IF1()      isec=0
_IF1()      return
_IF1()      entry userid(zanam)
_IF1()      zanam=zblank
_IF1()      return
_IF1()      end
_IF(apollo,hp700)
      function `dsum(n,sx,incx)'
      implicit REAL  (a-h,o-z)
      dimension sx(*)
      j=1
      if(n.gt.0) then
       if(incx.eq.1) then
        `dsum'=`vec_$dsum(sx,n)'
       else
        `dsum'=vec_$dsum_i(sx,incx,n)
       endif
      else
       `dsum'= 0.0d0
      endif
      return
      end
      subroutine `dfill(n,sa,sx,incx)'
      implicit REAL  (a-h,o-z)
      dimension sx(*)
      if(n.le.0) return
      if(incx.eq.1) then
       call vec_$dinit(sx,n,sa)
      else
       ii=1
       do 10 i=1,n
       sx(ii) = sa
10     ii = ii + incx
      endif
      return
      end
_ENDIF
_IFN(cray,sgi)
_IF(i8,_NOT(sgi))
      function smach (jj)
      implicit REAL  (a-h,o-z)
_IF(ibm,vax)
      dimension consts(3)
      data consts /z3410000000000000,
     +             z0010000000000000,
     +             z7fffffffffffffff /
      smach = consts(jj)
_ELSE
      data tol/0.0d0/
      if(jj.eq.1) smach=x02aaf(tol)
      if(jj.eq.2) smach=x02agf(tol)
      if(jj.eq.3) smach=x02acf(tol)
_ENDIF
      return
      end
_ENDIF
_ENDIF
_IF(ibm)
_IFN(assem)
      subroutine jksupe(hscm0,hscm1,hscm2,dscm0,dscm1,nscmf,l2)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension hscm0(*),hscm1(*),hscm2(*),dscm0(*),dscm1(*)
      common/blkin/gout(510),nint
      common/craypk/ij205(204),kl205(204)
c
      call upak6v(gout(num2ejk+num2ejk+1),ij205)
c
c---------both closed + open shells------------
c
      do   1 m = 1,nint
      ij = ij205(m)
      kl = kl205(m)
      goutpm = gout(m)
      goutkm = gout(num2ejk+m)
      hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
      hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
      do 2   i = 1,nscmf
      hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
      hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
      hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
      hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
      ij = ij + l2
    2 kl = kl + l2
    1 continue
      return
      end
      subroutine gsup(h,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer * 2 iij,kkl
      common/blkin/g(340),iij(340),kkl(340),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=iij(iw)
       kl=kkl(iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
      return
      end
      subroutine gsupa(h,p)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *2 iij,kkl
      common/blkin/g(408),iij(204),kkl(204),mword
      dimension h(*),p(*)
      do 1 iw=1,mword
      ij=iij(iw)
      kl=kkl(iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
      return
      end
      subroutine gmake(a,p,iky)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer * 2 i,j,k,l
      dimension a(*),p(*),iky(*)
      logical *1 iii(2),jjj(2),kkk(2),lll(2),mij,mkl
      common/blkin/g(340),mij(680),mkl(680),mword
      data i,j,k,l/4*0/
      equivalence(iii(1),i),(jjj(1),j),(kkk(1),k),(lll(1),l)
      equivalence (gg1,gg),(ikyj,jl),(ikyk,jk),(ikyi,il)
c
c ... computes symmetric g matrix
c ... g=2j<p>-k<p>
c ... p is symmetric density matrix
c ...
      do 60 iw=1,mword
      gg=g(iw)
      gg2=gg
      gg3=gg
      gg4=gg
      gg5=gg
      gg6=gg
      gg7=gg
      gg8=gg+gg
      gg9=gg8
      il=iw+iw
      iii(2)=mij(il-1)
      jjj(2)=mij(il)
      kkk(2)=mkl(il-1)
      lll(2)=mkl(il)
      ikyi=iky(i)
      ij=ikyi+j
      ik=ikyi+k
      il=ikyi+l
      ikyk=iky(k)
      kl=ikyk+l
      ikyj=iky(j)
      if(i.eq.j)then
      gg3=0.0d0
      gg4=0.0d0
      gg6=0.0d0
      gg7=0.0d0
      gg9=gg
      endif
      if(ij.eq.kl)then
      gg5=0.0d0
      gg6=0.0d0
      gg7=0.0d0
      gg9=0.0d0
      endif
      if(k.eq.l)then
      gg2=0.0d0
      gg4=0.0d0
      gg7=0.0d0
      gg8=gg
      endif
      if(i.eq.k)gg1=gg1+gg5
      if(j-k)4,5,6
    5 gg3=gg3+gg6
      go to 7
    4 gg3=gg6
    7 jk=ikyk+j
      if(j-l)8,9,10
    9 gg4=gg4+gg7
      go to 10
    8 gg4=gg7
      jl=iky(l)+j
      go to 11
    6 jk=ikyj+k
   10 jl=ikyj+l
   11 a(ij)=(gg8+gg8)*p(kl)+a(ij)
      a(kl)=(gg9+gg9)*p(ij)+a(kl)
      a(ik)=a(ik)-gg1*p(jl)
      a(jl)=a(jl)-gg4*p(ik)
      a(il)=a(il)-gg2*p(jk)
      a(jk)=a(jk)-gg3*p(il)
   60 continue
      return
      end
      subroutine proc2(a,p,b,q,iky)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer *2 i,j,k,l
      dimension a(*),p(*),b(*),q(*)
      dimension iky(*)
      logical *1 iii(2),jjj(2),kkk(2),lll(2),mij,mkl
      common/blkin/g(340),mij(680),mkl(680),mword
      data i,j,k,l/4*0/
      equivalence(iii(1),i),(jjj(1),j),(kkk(1),k),(lll(1),l)
c
c ... hamiltonian builder for rhf open shell .. spatially restricted
c ... define following symmetric density matrices
c ...  r1 .. closed shell
c ...  r2 .. open shell
c ... p=r1+r2/2
c ... setup ..
c ...     2j<p>-k<p> in a .... p=p
c ...      k<r2>/2   in b ....q=r2/2
c
      do 60 iw=1,mword
      gg=g(iw)
      gg1=gg
      gg2=gg
      gg3=gg
      gg4=gg
      gg5=gg
      gg6=gg
      gg7=gg
      gg8=gg+gg
      gg9=gg8
      il=iw+iw
      iii(2)=mij(il-1)
      jjj(2)=mij(il)
      kkk(2)=mkl(il-1)
      lll(2)=mkl(il)
      ikyi=iky(i)
      ij=ikyi+j
      ik=ikyi+k
      il=ikyi+l
      ikyk=iky(k)
      kl=ikyk+l
      ikyj=iky(j)
      if(i.eq.j)then
      gg3=0.0d0
      gg4=0.0d0
      gg6=0.0d0
      gg7=0.0d0
      gg9=gg
      endif
      if(ij.eq.kl)then
      gg5=0.0d0
      gg6=0.0d0
      gg7=0.0d0
      gg9=0.0d0
      endif
      if(k.eq.l)then
      gg2=0.0d0
      gg4=0.0d0
      gg7=0.0d0
      gg8=gg
      endif
      if(i.eq.k)gg1=gg1+gg5
      if(j-k)4,5,6
    5 gg3=gg3+gg6
      go to 7
    4 gg3=gg6
    7 jk=ikyk+j
      if(j-l)8,9,10
    9 gg4=gg4+gg7
      go to 10
    8 gg4=gg7
      jl=iky(l)+j
      go to 11
    6 jk=ikyj+k
   10 jl=ikyj+l
   11 a(ij)=(gg8+gg8)*p(kl)+a(ij)
      a(kl)=(gg9+gg9)*p(ij)+a(kl)
      a(ik)=a(ik)-gg1*p(jl)
      a(jl)=a(jl)-gg4*p(ik)
      a(il)=a(il)-gg2*p(jk)
      a(jk)=a(jk)-gg3*p(il)
      b(ik)=b(ik)+gg1*q(jl)
      b(jl)=b(jl)+gg4*q(ik)
      b(il)=b(il)+gg2*q(jk)
   60 b(jk)=b(jk)+gg3*q(il)
      return
      end
      subroutine szero(c,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension c(*)
      if(n.gt.0)then
      do 1 loop=1,n
    1 c(loop)=0.0d0
      endif
      return
      end
_ENDIF
_ENDIF
_IFN(cray)
      subroutine setstz(ntimes,int,conf)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IF(ibm,vax)
      dimension conf(*)
      flag=pad(int)
      do 1 i=1,ntimes
   1  conf(i)=flag
_ELSEIF(convex,i8drct)
      integer *8 conf
      dimension conf(*)
      do 1 i=1,ntimes
   1  conf(i)=int
_ELSEIF(ksr)
      integer conf
      dimension conf(*)
      do 1 i=1,ntimes
   1  conf(i)=int
_ELSE
      integer*4 conf(2,*),int
      do 1 i=1,ntimes
_IF(littleendian)
        conf(1,i) = int
        conf(2,i) = 0
_ELSE
        conf(1,i) = 0
        conf(2,i) = int
_ENDIF
   1  continue
_ENDIF
      return
      end
_ENDIF
_IF1()   c logic for following routine inverted 
_IF1()   c_IIF(alliant,titan,apollo,sun,sgi,ipsc,rs6000,dec,hp700,ksr)
_IFN(cray,ibm,convex) 
      subroutine vsqrt(a,incri,b,incrj,n)
      implicit REAL  (a-h,o-z)
      dimension a(*),b(*)
      loopi=1
      loopj=1
      do 10 loop=1,n
         b(loopj)= dsqrt(a(loopi))
         loopj=loopj+incrj
         loopi=loopi+incri
 10   continue
      return
      end
      subroutine vmul(a,ia,b,ib,c,ic,n)
      implicit REAL  (a-h,o-z)
      dimension c(*),a(*),b(*)
       if(ia.eq.1. and. ib.eq.1. and .ic.eq.1) then
         do 10 loop=1,n
10       c(loop)=a(loop)*b(loop)
      else
         loopc=1
         loopa=1
         loopb=1
         do 20 loop=1,n
            c(loopc)=a(loopa)*b(loopb)
            loopa=loopa+ia
            loopb=loopb+ib
            loopc=loopc+ic
 20      continue
      endif
      return
      end
      subroutine vfill(scaler,c,ic,n)
      implicit REAL  (a-h,o-z)
      dimension c(*)
      loopc=1
      do 10 loop=1,n
         c(loopc)=scaler
         loopc=loopc+ic
 10   continue
      return
      end
      subroutine vma(a,ia,b,ib,c,ic,d,id,n)
      implicit REAL  (a-h,o-z)
      dimension d(*),c(*),a(*),b(*)
      loopc=1
      loopa=1
      loopb=1
      loopd=1
      do 10 loop=1,n
         d(loopd)=a(loopa)*b(loopb) +  c(loop)
         loopa=loopa+ia
         loopb=loopb+ib
         loopc=loopc+ic
         loopd=loopd+id
 10   continue
      return
      end
      subroutine vneg(a,incri,b,incrj,n)
      implicit REAL  (a-h,o-z)
      dimension a(*),b(*)
      loopi=1
      loopj=1
      do 10 loop=1,n
         b(loopj)=-a(loopi)
         loopj=loopj+incrj
         loopi=loopi+incri
 10   continue
c
      return
      end
      subroutine rmtran(a,mrowa,r,mrowr,ncol,nrow)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character  (z),character *1 (x)
      implicit character *4 (y)
      dimension a(mrowa,*),r(mrowr,*)
      if(ncol.le.nrow) then
         do 10 j=1,ncol
            do 20 i=1,nrow
               r(i,j)=a(j,i)
 20         continue
 10      continue
      else
         do 30 j=1,nrow
            do 40 i=1,ncol
               r(j,i)=a(i,j)
 40         continue
 30      continue
      endif
      return
      end
      subroutine maxmgv(x,ix,valmax,imax,n)
      implicit REAL (a-h,o-z)
      dimension x(*)
c
      valmax = 0.0d0
      imax = 0
      ipt = 1
      do 10 i = 1,n
         ax = dabs(x(ipt))
         if (ax.ge.valmax) then
            valmax = ax
            imax = i
         endif
         ipt = ipt + ix
 10   continue
c
      return
      end
_ENDIF
_IFN(ibm,convex) 
      subroutine vdiv(a,ia,b,ib,c,ic,n)
      implicit REAL  (a-h,o-z)
      dimension c(*),a(*),b(*)
      loopc=1
      loopa=1
      loopb=1
      do 1 loop=1,n
         c(loopc)=a(loopa)/b(loopb)
         loopa=loopa+ia
         loopb=loopb+ib
         loopc=loopc+ic
 1    continue
      return
      end
      subroutine vsmsb(b,ib,scalar,a,ia,r,ir,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character  (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),b(*)
      loopa = 1
      loopb = 1
      loopr = 1
      do 10 loop=1,n
         r(loopr)=scalar*b(loopb)-a(loopa)
         loopa = loopa + ia
         loopr = loopr + ir
         loopb = loopb + ib
 10   continue
      return
      end
_ENDIF
_IF1()   c logic for following routine inverted 
_IF1()   c_IIF(titan,sun,sgi,ipsc,rs6000,dec,ksr)
_IFN(cray,hp700,apollo,ibm,convex)
      subroutine vmov(a,ia,b,ib,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),b(*)
      if(ia.eq.1. and. ib.eq.1) then
        do 20 loop=1,n
20      b(loop)=a(loop)
      else
        inca = 1
        incb = 1
          do 10 loop=1,n
           b(incb)=a(inca)
          inca=inca+ia
10        incb=incb+ib
      endif
      return
      end
      subroutine vadd(x,incx,y,incy,z,incz,n)
      implicit REAL  (a-h,o-z)
      dimension x(*),y(*),z(*)
      if(incx.eq.1. and .incy.eq.1.  and  .incz.eq.1) then
        do 10 loop=1,n
10      z(loop) = x(loop)+y(loop)
      else
        ix =1
        iy =1
        iz =1
        do 20 loop=1,n
           z(iz) = x(ix) + y(iy)
           ix = ix + incx
           iy = iy + incy
           iz = iz + incz
20      continue
      endif
      return
      end
      subroutine vsmul (a,ia,scaler,c,ic,n)
      implicit REAL (a-h,o-z)
      dimension a(*),c(*)
      if(ia.eq.1. and .ic.eq.1) then
         do 10 loop=1,n
10       c (loop) = a (loop) * scaler
      else
         j=1
         k=1
         do 20 i=1,n
            c (j) = a (k) * scaler
            k = k +ia
            j = j + ic
20       continue
      endif
      return
      end
      subroutine vsma(b,ib,scalar,a,ia,r,ir,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character  (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),b(*)
      loopr=1
      loopa=1
      loopb=1
      do 10 i=1,n
         r(loopr)=a(loopa)+scalar*b(loopb)
         loopr = loopr + ir
         loopa = loopa + ia
         loopb = loopb + ib
 10   continue
      return
      end
      subroutine vsadd(a,ia,scaler,c,ic,n)
      implicit REAL  (a-h,o-z)
      dimension a(*),c(*)
      loopa=1
      loopc=1
      do 10 loop=1,n
         c(loopc)=a(loopa)+scaler
         loopc=loopc+ic
         loopa=loopa+ia
 10   continue
      return
      end
_ENDIF
_IF1()   c logic for following routine inverted 
_IF1()   c_IIF(titan,sun,sgi,ipsc,rs6000,dec,hp700,ksr)
_IFN(cray,ibm,apollo,convex)
      subroutine vclr(a,incr,n)
      implicit REAL  (a-h,o-z)
      dimension a(*)
      if(incr.eq.1) then
        do 10 loop=1,n
10      a(loop) = 0.0d0
      else
        loopi=1
        do 20 loop=1,n
           a(loopi)=0.0d0
           loopi=loopi+incr
20      continue
      endif
      return
      end
      subroutine vsub(x,incx,y,incy,z,incz,n)
      implicit REAL  (a-h,o-z)
      dimension x(*),y(*),z(*)
      if(incx.eq.1. and .incy.eq.1.  and  .incz.eq.1) then
        do 10 loop=1,n
10      z(loop) = x(loop)-y(loop)
      else
        ix =1
        iy =1
        iz =1
        do 20 loop=1,n
           z(iz) = x(ix) - y(iy)
           ix = ix + incx
           iy = iy + incy
           iz = iz + incz
20      continue
      endif
      return
      end
_ENDIF
_IF(apollo,hp700)
      subroutine vmov(a,ia,b,ib,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension a(*),b(*)
      if(n.le.0) return
      if(ia.eq.1. and. ib.eq.1) then
        call `vec_$dcopy(a,b,n)'
      else
        call vec_$dcopy_i(a,ia,b,ib,n)
      endif
      return
      end
_IFN(hp700)
      subroutine vclr(a,incr,n)
      implicit REAL  (a-h,o-z)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension a(*)
      if(n.le.0) return
       if(incr.eq.1) then
         call vec_$dzero(a,n)
       else
         call vec_$dzero_i(a,incr,n)
       endif
      return
      end
_ENDIF
      subroutine vadd(x,incx,y,incy,z,incz,n)
      implicit REAL  (a-h,o-z)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension x(*),y(*),z(*)
      if(n.le.0) return
      if(incx.eq.1. and. incy.eq.1. and. incz.eq.1) then
         call vec_$dadd_vector(x,y,n,z)
      else
         call vec_$dadd_vector_i(x,incx,y,incy,n,z,incz)
      endif
      return
      end
_IFN(hp700)
      subroutine vsub(x,incx,y,incy,z,incz,n)
      implicit REAL  (a-h,o-z)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension x(*),y(*),z(*)
      if(n.le.0) return
      if(incx.eq.1. and. incy.eq.1. and. incz.eq.1) then
         call vec_$dsub(x,y,n,z)
      else
         call vec_$dsub_i(x,incx,y,incy,n,z,incz)
      endif
      return
      end
_ENDIF
      subroutine vsmul (a,ia,scaler,c,ic,n)
      implicit REAL (a-h,o-z)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension a(*),c(*)
      if(n.le.0) return
      if(ia.eq.1. and. ic.eq.1) then
_IF(apollo)
         call vec_$dmult_constant(a,n,scaler,c)
      else
         call vec_$dmult_constant_i(a,ia,n,scaler,c,ic)
_ELSE
         do i=1,n
            c(i) = scaler * a(i)
         enddo
      else
         do i=0,n-1
            c(i*ic+1) = scaler * a(i*ia+1)
         enddo
_ENDIF
      endif
      return
      end
      subroutine vsma(b,ib,scalar,a,ia,r,ir,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character  (z),character *1 (x)
      implicit character *4 (y)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension a(*),r(*),b(*)
      if(n.le.0) return
      if(ib.eq.1. and. ia.eq.1. and. ir.eq.1) then
      call vec_$dmult_add(a,b,n,scalar,r)
      else
      call vec_$dmult_add_i(a,ia,b,ib,n,scalar,r,ir)
      endif
      return
      end
      subroutine vsadd(a,ia,scaler,c,ic,n)
      implicit REAL  (a-h,o-z)
_IFN1(b)%include '/sys/ins/base.ins.ftn'
_IFN1(b)%include '/sys/ins/vec.ins.ftn'
      dimension a(*),c(*)
      if(n.le.0) return
      if(ia.eq.1. and. ic.eq.1) then
         call vec_$dadd_constant(a,n,scaler,c)
      else
         call vec_$dadd_constant_i(a,ia,n,scaler,c,ic)
      endif
      return
      end
_ENDIF
_IF(cray)
      subroutine vsadd(a,ia,scaler,c,ic,n)
      implicit none
      REAL a(*),c(*), scaler
      integer ia, ic, n, loopa, loopc, loop
      if(n.le.0) return
      if(ia.eq.1. and. ic.eq.1) then
         call vsav(n,scaler,c,a)
      else
         loopa=1
         loopc=1
         do 10 loop=1,n
            c(loopc)=a(loopa)+scaler
            loopc=loopc+ic
            loopa=loopa+ia
 10      continue
      endif
      return
      end
_ENDIF
_IF(convex,titan)
      subroutine jandk2(l2,nshell,njk,h3,p,iq,jq,kq,lq,gq,pg,
     * ind,nint)
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension h3(l2*nshell*2*njk),p(l2*nshell),gq(nint),
     &          pg(6*nshell*nint)
      integer iq(nint),jq(nint),kq(nint),lq(nint),ind(6*nshell*nint)
c
c build j and k matrices for nshell density matrices.
c on xmp all code for nshell<3 should vectorize and jkscat can use
c vector scatter to handle overwrite logic with a vector
c length of nshell*njk*6 <= 64 (128 on c2). crucial assumption in logic.
c
c for nshell<3 njk copies of the j and k matrices are made to
c enhance vectorisation. for larger values of nshell njk is ignored.
c
c l2 = size of triangle
c nshell = no. of shells
c njk = no. of copies to enhance vectorisation/parallelisation
c h3 = h3(l2,nshell,itype,njk), itype = 1 is k, 2 is j
c p = p(l2,nshell) input density matrices
c iq,jq,kq,lq = integral indices
c gq = integrals
c pg,ind = (6*nint*njk*nshell) workspace, only for nshell<3
c nint = no. of input integrals
c
      if (nint.le.0) return
_IF1(c)      if(nshell*njk*6.gt.64)call caserr2(' nshell*njk too large',0)
_IF1(x)      if(nshell*njk*6.gt.128)call caserr2(' nshell*njk too large',0)
_IF1(t)      if(nshell*njk*6.gt.32)call caserr2(' nshell*njk too large',0)
      mskip = l2*nshell*2
      icoff = l2*nshell
      if (nshell.eq.1) then
        ic4 = 1
_IF1(c)cdir$ align
_IF1(x)c$dir no_recurrence
_IF1(tc)cdir$ ivdep
        do 51 iww = 1,nint
            i = iq(iww)
            j = jq(iww)
            k = kq(iww)
            l = lq(iww)
            ikyi = iky(i)
            ikyj = iky(j)
            ikyk = iky(k)
            ikyl = iky(l)
            ij = ikyi + j
            ik = ikyi + k
            il = ikyi + l
            jk = max(ikyj + k,ikyk + j)
            jl = max(ikyj + l,ikyl + j)
            kl = ikyk + l
c
            g = gq(iww)
            g2 = g + g
_IF1()c            xik = cvmgt(g2,g,(i.eq.k .or. j.eq.l))
_IF1()c            xil = cvmgt(g2,g,(j.eq.k))
            xik = g
            if(i.eq.k .or. j.eq.l) xik = g2
            xil = g
            if(j.eq.k) xil = g2
c
            moff = mod(iww,njk) * mskip
            pg(ic4  ) = xik * p(jl)
            pg(ic4+1) = xil * p(jk)
            pg(ic4+2) = xil * p(il)
            pg(ic4+3) = xik * p(ik)
            pg(ic4+4) =   g * p(kl)
            pg(ic4+5) =   g * p(ij)
c
            ind(ic4  ) = ik + moff
            ind(ic4+1) = il + moff
            ind(ic4+2) = jk + moff
            ind(ic4+3) = jl + moff
            ind(ic4+4) = ij + moff + icoff
            ind(ic4+5) = kl + moff + icoff
            ic4 = ic4 + 6
51      continue
      else if(nshell.eq.2) then
c loop is messy to force vectorisation with cft 1.16
        ic4 = 1
        ic42 = 7
_IF1()cdir$ align,ivdep
c$dir no_recurrence
_IF1(ct)cdir$ ivdep
        do 52 iww = 1,nint
            ikyi = iky(iq(iww))
            ij = ikyi + jq(iww)
            ik = ikyi + kq(iww)
            il = ikyi + lq(iww)
            jk = max(iky(jq(iww)) + kq(iww),iky(kq(iww)) + jq(iww))
            jl = max(iky(jq(iww)) + lq(iww),iky(lq(iww)) + jq(iww))
            kl = iky(kq(iww)) + lq(iww)
c
            g = gq(iww)
            g2 = g+g
_IF1()c            xik = cvmgt((2.*gq(iww)),gq(iww),
_IF1()c     &            (iq(iww).eq.kq(iww) .or. jq(iww).eq.lq(iww)))
_IF1()c            xil = cvmgt((2.*gq(iww)),gq(iww),(jq(iww).eq.kq(iww)))
            xik = g
            if(iq(iww).eq.kq(iww) .or. jq(iww).eq.lq(iww))
     &        xik = g2
            xil = g
            if(jq(iww).eq.kq(iww)) xil = g2
c
            moff = mod(iww,njk) * mskip
            pg(ic4  ) = xik * p(jl)
            pg(ic4+1) = xil * p(jk)
            pg(ic4+2) = xil * p(il)
            pg(ic4+3) = xik * p(ik)
            pg(ic4+4) =   g * p(kl)
            pg(ic4+5) =   g * p(ij)
c
            ind(ic4  ) = ik + moff
            ind(ic4+1) = il + moff
            ind(ic4+2) = jk + moff
            ind(ic4+3) = jl + moff
            ind(ic4+4) = ij + moff + icoff
            ind(ic4+5) = kl + moff + icoff
            ic4 = ic4 + 12
c
            ij = ij + l2
            ik = ik + l2
            il = il + l2
            jk = jk + l2
            jl = jl + l2
            kl = kl + l2
c
            pg(ic42  ) = xik * p(jl)
            pg(ic42+1) = xil * p(jk)
            pg(ic42+2) = xil * p(il)
            pg(ic42+3) = xik * p(ik)
            pg(ic42+4) =   g * p(kl)
            pg(ic42+5) =   g * p(ij)
c
            ind(ic42  ) = ik + moff
            ind(ic42+1) = il + moff
            ind(ic42+2) = jk + moff
            ind(ic42+3) = jl + moff
            ind(ic42+4) = ij + moff + icoff
            ind(ic42+5) = kl + moff + icoff
            ic42 = ic42 + 12
52      continue
      else
        ic4 = 1
_IF1()cdir$ align
        do 53 iww = 1,nint
            i = iq(iww)
            j = jq(iww)
            k = kq(iww)
            l = lq(iww)
            ikyi = iky(i)
            ikyj = iky(j)
            ikyk = iky(k)
            ikyl = iky(l)
            ij = ikyi + j
            ik = ikyi + k
            il = ikyi + l
            jk = max(ikyj + k,ikyk + j)
            jl = max(ikyj + l,ikyl + j)
            kl = ikyk + l
c
            g = gq(iww)
            g2 = g + g
            xik = g
            if(i.eq.k .or. j.eq.l) xik = g2
            xil = g
            if(j.eq.k) xil = g2
c$dir scalar
          do 15 ishell = 1,nshell
              cij = h3(ij+icoff) + g*p(kl)
              h3(kl+icoff) = h3(kl+icoff) + g*p(ij)
              h3(ij+icoff) = cij
              hik = h3(ik) + xik * p(jl)
              hil = h3(il) + xil * p(jk)
              hjk = h3(jk) + xil * p(il)
              h3(jl) = h3(jl) + xik * p(ik)
              h3(jk) = hjk
              h3(il) = hil
              h3(ik) = hik
              ij = ij + l2
              ik = ik + l2
              il = il + l2
              jk = jk + l2
              jl = jl + l2
              kl = kl + l2
15        continue
53      continue
      endif
c
      if (nshell.le.2) then
        njkl = mod(nint,njk)
        call jkscat(h3,pg,ind,nint*nshell*6,njk*nshell*6,njkl*nshell*6)
      endif
c
      return
      end
      subroutine jkscat(h,pg,ind,n6,njk6,njkl6)
      implicit REAL  (a-h,o-z)
      dimension h(*),pg(*),ind(*)
      if (n6.le.0) return
_IF(titan)
c
c scalar or short vector code ... use of vector scatter
c to perform overwrite logic gives factor of 6 on vl.
c
      njk = njk6/6
      njkl = njkl6/6
      if(njkl.eq.0) njkl = njk
      mjk = njk
      ic = n6/6
      nbloc = (ic-1)/njk + 1
      ic4 = 1
      do 10 ibloc = 1,nbloc
          if(ibloc.eq.nbloc) mjk = njkl
_IF1()cdir$ ivdep,shortloop,ivdmo
_IF1(t)cdir$ ivdep
         do 20 m = 1,mjk
             hik = h(ind(ic4  )) + pg(ic4)
             hil = h(ind(ic4+1)) + pg(ic4+1)
             hjk = h(ind(ic4+2)) + pg(ic4+2)
             hjl = h(ind(ic4+3)) + pg(ic4+3)
             h(ind(ic4+3)) = hjl
             h(ind(ic4+2)) = hjk
             h(ind(ic4+1)) = hil
             h(ind(ic4  )) = hik
             hij = h(ind(ic4+4)) + pg(ic4+4)
             hkl = h(ind(ic4+5)) + pg(ic4+5)
             h(ind(ic4+4)) = hij
             h(ind(ic4+5)) = hkl
             ic4 = ic4+6
20        continue
10    continue
_ENDIF
_IF(convex)
c
c following is cray/convex or vector register specific
_IF1()cc equivalent to the above. use cal version.
      nbloc = (n6-1)/njk6
      njk6m1 = njk6 - 1
      ic4 = 1
      do 10 ibloc = 1,nbloc
_IF1()cdir$ ivdep,shortloop,ivdmo
_IF1(x)c$dir no_recurrence
_IF1(ct)cdir$ ivdep
      do 20 m = ic4,ic4+njk6m1
          h(ind(m)) = h(ind(m)) + pg(m)
20        continue
      ic4 = ic4+njk6
10    continue
      ntemp = njkl6
      if(ntemp.eq.0) ntemp=njk6
_IF1()cdir$ ivdep,shortloop,ivdmo
_IF1(x)c$dir no_recurrence
_IF1(tc)cdir$ ivdep
      do 30 m = ic4,ic4+ntemp-1
          h(ind(m)) = h(ind(m)) + pg(m)
30    continue
_ENDIF
      return
      end
      subroutine upak8z(nword,iii,i,j,k,l)
      implicit REAL  (a-h,o-z)
      logical *1 iii,i,j,k,l
      dimension iii(*),i(4,*),j(4,*),k(4,*),l(4,*)
      kk=1
      do 1 loop=1,nword
      i(4,loop)=iii(kk  )
      j(4,loop)=iii(kk+1)
      k(4,loop)=iii(kk+2)
      l(4,loop)=iii(kk+3)
  1   kk=kk+4
      return
      end
_ENDIF
      subroutine dagger(ncol,nrow,a,mrowa,r,mrowr)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(mrowa,*),r(mrowr,*)
      if(ncol.gt.nrow)go to 3
      do 1 j=1,ncol
      do 1 i=1,nrow
   1  r(i,j)=a(j,i)
      return
   3  do 2 j=1,nrow
      do 2 i=1,ncol
   2  r(j,i)=a(i,j)
      return
      end
      subroutine symm1(r,a,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(n,*)
      m=1
      do 1 i=1,n
      do 1 j=1,i
      r(m)=a(i,j)+a(j,i)
   1  m=m+1
      return
      end
      function nocat1(label,nf,itext)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension label(*)
      do 1 loop=1,nf
      if(label(loop).ne.itext)go to 2
    1 continue
      nocat1=0
      return
    2 nocat1=loop
      return
      end
      subroutine ibasgn(n,ibegin,iadd,iarr)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension iarr(*)
      j=ibegin
      do 1 i=1,n
      iarr(i)=j
    1 j=j+iadd
      return
      end
_IF(f77)
_IFN(rs6000,sun,sgi,dec,cray,t3d)
      function rand(is)
c should work on any 32 bit machine. actual generator
c is only fair and not suitable for detailed work.
      integer is,imult,imod,is1,is2,iss2
      data imult/16807/,imod/2147483647/,scale/4.656612875d-10/
c
c is = mod(is*16807,2**31-1).
      if(is.le.0) is = 1
      is2 = mod(is,32768)
      is1 = (is-is2)/32768
      iss2 = is2 * imult
      is2 = mod(iss2,32768)
      is1 = mod(is1*imult+(iss2-is2)/32768,65536)
      is = mod(is1*32768+is2,imod)
      c = scale * dfloat(is)
      rand=c
      return
      end
_ENDIF
_ENDIF
_IF(c90)
      subroutine szero(c,n)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension c(*)
      if(n.gt.0)then
        do 1 loop=1,n
    1   c(loop)=0.0d0
      endif
      return
      end
      subroutine scatt(n,scalar,v,mapv)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c  CALL SCATT(N,SCALAR,V,MAPV)
c  V(MAPV) = SCALAR FOR N ELEMENTS
c  invoked by DIAG on cray
c
      dimension v(*), mapv(*)
c
      do 10 loop=1,n
 10   v(mapv(loop)) = scalar
c
      return
      end
_ENDIF
_IF(cray)
      subroutine `vadd(x,incx,y,incy,z,incz,n)'
      implicit REAL  (a-h,o-z)
      dimension x(*),y(*),z(*)
      if(incx.eq.1. and .incy.eq.1.  and  .incz.eq.1) then
        do loop=1,n
          z(loop) = x(loop)+y(loop)
        enddo
      else
        ix =1
        iy =1
        iz =1
        do  loop=1,n
           z(iz) = x(ix) + y(iy)
           ix = ix + incx
           iy = iy + incy
           iz = iz + incz
        enddo
      endif
      return
      end
      subroutine vneg(a,incri,b,incrj,n)
      implicit REAL  (a-h,o-z)
      dimension a(*),b(*)
      loopi=1
      loopj=1
      do loop=1,n
         b(loopj)=-a(loopi)
         loopj=loopj+incrj
         loopi=loopi+incri
      enddo
      return
      end
_ENDIF
_IF(t3d)
      subroutine sgthr(n,a,r,map)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),map(*)
      if(n.gt.0)then
       do 1 loop=1,n
    1  r(loop)=a(map(loop))
      endif
      return
      end
      subroutine ssctr(n,a,map,r)
      implicit REAL  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),map(*)
      if(n.gt.0) then
       do 1 loop=1,n
    1  r(map(loop))=a(loop)
      endif
      return
      end
_ENDIF
      subroutine izero(nel, iarray, incx)
      implicit REAL  (a-h,o-z),integer  (i-n)
      dimension iarray(*)
      if(nel.le.0)return
      if(incx.eq.1) then
       do 1 loop = 1,nel
   1   iarray(loop)=0
      else
       ii = 1
       do 2 loop = 1,nel
       iarray(ii) = 0
   2   ii = ii + incx
      endif
      return
      end
_IFN(cray,t3d)
      subroutine zzero(array,nel)
      implicit REAL  (a-h,o-z),integer  (i-n)
_IF(convex,i8drct)
      integer *8 array,z0
_ELSEIF(ksr,i8)
      integer array,z0
_ELSE
      REAL array,z0
_ENDIF
      dimension  array(*)
_IF(ksr,convex,i8,i8drct)
      data z0/0/
_ELSEIF(bits8)
      call pad8(0,z0)
_ELSE
      z0=pad(0)
_ENDIF
      if(nel.le.0)return
      do 1 loop =1,nel
   1   array(loop)=z0
      return
      end
_ENDIF
c
c  print out of RCS version numbers and dates
c     file name sorted in alphabetical order per module
c
      subroutine prtver(iwr)
      integer iwr
c
c     In the version routines the strings are declared and used
c     as follows:
c
c        declaration             usage          length
c        -----------             -----          ------
c        character*80 source     source(9:)     72
c        character*30 revision   revision(11:)  20
c        character*60 date       data(7:)       54
c
c     so the declarations below should match the length above...
c
      character*20 r 
      character*54 d 
      character*72 s

      write(iwr,102)
      write(iwr,101)
      write(iwr,102)
      call ver_anala(s,r,d)
      if(s.ne." ")write(iwr,100)"anala",r,s,d
      call ver_analb(s,r,d)
      if(s.ne." ")write(iwr,100)"analb",r,s,d
      call ver_analc(s,r,d)
      if(s.ne." ")write(iwr,100)"analc",r,s,d
      call ver_anald(s,r,d)
      if(s.ne." ")write(iwr,100)"anald",r,s,d
      call ver_anale(s,r,d)
      if(s.ne." ")write(iwr,100)"anale",r,s,d
      call ver_analf(s,r,d)
      if(s.ne." ")write(iwr,100)"analf",r,s,d
      call ver_analg(s,r,d)
      if(s.ne." ")write(iwr,100)"analg",r,s,d
      call ver_basis(s,r,d)
      if(s.ne." ")write(iwr,100)"basis",r,s,d
      call ver_basis1(s,r,d)
      if(s.ne." ")write(iwr,100)"basis1",r,s,d
      call ver_basis2(s,r,d)
      if(s.ne." ")write(iwr,100)"basis2",r,s,d
      call ver_basis3(s,r,d)
      if(s.ne." ")write(iwr,100)"basis3",r,s,d
      s = " "
      r = " "
      d = " "
      call ver_c(s,len(s),r,len(r),d,len(d))
      if(s.ne." ")write(iwr,100)"c",r,s,d
      call ver_casa(s,r,d)
      if(s.ne." ")write(iwr,100)"casa",r,s,d
      call ver_casb(s,r,d)
      if(s.ne." ")write(iwr,100)"casb",r,s,d
      call ver_ccsd(s,r,d)
      if(s.ne." ")write(iwr,100)"ccsd",r,s,d
      call ver_cphf(s,r,d)
      if(s.ne." ")write(iwr,100)"cphf",r,s,d
      call ver_derdrv(s,r,d)
      if(s.ne." ")write(iwr,100)"derdrv",r,s,d
      call ver_dft(s,r,d)
      if(s.ne." ")write(iwr,100)"dft",r,s,d
      call ver_dircta(s,r,d)
      if(s.ne." ")write(iwr,100)"dircta",r,s,d
      call ver_dirctb(s,r,d)
      if(s.ne." ")write(iwr,100)"dirctb",r,s,d
      call ver_dirctc(s,r,d)
      if(s.ne." ")write(iwr,100)"dirctc",r,s,d
      call ver_dirctd(s,r,d)
      if(s.ne." ")write(iwr,100)"dirctd",r,s,d
      call ver_direct(s,r,d)
      if(s.ne." ")write(iwr,100)"direct",r,s,d
      call ver_dirrpa(s,r,d)
      if(s.ne." ")write(iwr,100)"dirrpa",r,s,d
      call ver_drv1e(s,r,d)
      if(s.ne." ")write(iwr,100)"drv1e",r,s,d
      call ver_drv2e(s,r,d)
      if(s.ne." ")write(iwr,100)"drv2e",r,s,d
      call ver_drv80(s,r,d)
      if(s.ne." ")write(iwr,100)"drv80",r,s,d
      call ver_drvmp(s,r,d)
      if(s.ne." ")write(iwr,100)"drvmp",r,s,d
      call ver_fullci(s,r,d)
      if(s.ne." ")write(iwr,100)"fullci",r,s,d
_IF(ga)
      call ver_ga(s,r,d)
      if(s.ne." ")write(iwr,100)"ga",r,s,d
_ENDIF
      call ver_gff(s,r,d)
      if(s.ne." ")write(iwr,100)"gff",r,s,d
      call ver_guess(s,r,d)
      if(s.ne." ")write(iwr,100)"guess",r,s,d
      call ver_index4(s,r,d)
      if(s.ne." ")write(iwr,100)"index4",r,s,d
      call ver_inpci(s,r,d)
      if(s.ne." ")write(iwr,100)"inpci",r,s,d
      call ver_input(s,r,d)
      if(s.ne." ")write(iwr,100)"input",r,s,d
      call ver_intega(s,r,d)
      if(s.ne." ")write(iwr,100)"intega",r,s,d
      call ver_integb(s,r,d)
      if(s.ne." ")write(iwr,100)"integb",r,s,d
      call ver_integb_lib(s,r,d)
      if(s.ne." ")write(iwr,100)"integb_lib",r,s,d
      call ver_integb_nl(s,r,d)
      if(s.ne." ")write(iwr,100)"integb_nl",r,s,d
      call ver_integc(s,r,d)
      if(s.ne." ")write(iwr,100)"integc",r,s,d
      call ver_integd(s,r,d)
      if(s.ne." ")write(iwr,100)"integd",r,s,d
      call ver_intege(s,r,d)
      if(s.ne." ")write(iwr,100)"intege",r,s,d
      call ver_integs(s,r,d)
      if(s.ne." ")write(iwr,100)"integs",r,s,d
      call ver_machci(s,r,d)
      if(s.ne." ")write(iwr,100)"machci",r,s,d
      call ver_machscf(s,r,d)
      if(s.ne." ")write(iwr,100)"machscf",r,s,d
      call ver_mainci(s,r,d)
      if(s.ne." ")write(iwr,100)"mainci",r,s,d
      call ver_mains(s,r,d)
      if(s.ne." ")write(iwr,100)"mains",r,s,d
      call ver_mass(s,r,d)
      if(s.ne." ")write(iwr,100)"mass",r,s,d
      call ver_master(s,r,d)
      if(s.ne." ")write(iwr,100)"master",r,s,d
_IF(never)
      call ver_mcdab_ga(s,r,d)
      if(s.ne." ")write(iwr,100)"mcdab_ga",r,s,d
_ENDIF
      call ver_mcscfa(s,r,d)
      if(s.ne." ")write(iwr,100)"mcscfa",r,s,d
      call ver_mclr(s,r,d)
      if(s.ne." ")write(iwr,100)"mclr",r,s,d
      call ver_mcscfb(s,r,d)
      if(s.ne." ")write(iwr,100)"mcscfb",r,s,d
      call ver_mcscfc(s,r,d)
      if(s.ne." ")write(iwr,100)"mcscfc",r,s,d
      call ver_model(s,r,d)
      if(s.ne." ")write(iwr,100)"model",r,s,d
      call ver_morokuma(s,r,d)
      if(s.ne." ")write(iwr,100)"morokuma",r,s,d
_IF(mp2_parallel)
      call ver_mp2_parallel(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2_parallel",r,s,d
_ENDIF
_IF(masscf)
      call ver_masscf(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2_parallel",r,s,d
_ENDIF
      call ver_mp2a(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2a",r,s,d
      call ver_mp2b(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2b",r,s,d
      call ver_mp2c(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2c",r,s,d
      call ver_mp2d(s,r,d)
      if(s.ne." ")write(iwr,100)"mp2d",r,s,d
      call ver_mp3(s,r,d)
      if(s.ne." ")write(iwr,100)"mp3",r,s,d
      call ver_mrdci1(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci1",r,s,d
      call ver_mrdci2(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci2",r,s,d
      call ver_mrdci3(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci3",r,s,d
      call ver_mrdci4(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci4",r,s,d
      call ver_mrdci5(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci5",r,s,d
      call ver_mrdci6(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci6",r,s,d
      call ver_mrdci7(s,r,d)
      if(s.ne." ")write(iwr,100)"mrdci7",r,s,d
cjmht      call ver_nag(s,r,d)
cjmht      if(s.ne." ")write(iwr,100)"nag",r,s,d
_IF(mrdci)
      call ver_newmrd1(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd1",r,s,d
      call ver_newmrd2(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd2",r,s,d
      call ver_newmrd3(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd3",r,s,d
      call ver_newmrd4(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd4",r,s,d
      call ver_newmrd5(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd5",r,s,d
      call ver_newmrd6(s,r,d)
      if(s.ne." ")write(iwr,100)"newmrd6",r,s,d
_ENDIF
_IF(nbo)
      call ver_nbo(s,r,d)
      if(s.ne." ")write(iwr,100)"nbo",r,s,d
_ENDIF
      call ver_nvccsd(s,r,d)
      if(s.ne." ")write(iwr,100)"nvccsd",r,s,d
      call ver_optef(s,r,d)
      if(s.ne." ")write(iwr,100)"optef",r,s,d
      call ver_optim(s,r,d)
      if(s.ne." ")write(iwr,100)"optim",r,s,d
_IF(never)
      call ver_pack(s,r,d)
      if(s.ne." ")write(iwr,100)"pack",r,s,d
_ENDIF
      call ver_parallel(s,r,d)
      if(s.ne." ")write(iwr,100)"parallel",r,s,d
_IF(peigs)
      s = " "
      r = " "
      d = " "
      call ver_peigs_interface(s,len(s),r,len(r),d,len(d))
      if(s.ne." ")write(iwr,100)"peigs_interface",r,s,d
_ENDIF
_IF(never)
      call ver_pdiag(s,r,d)
      if(s.ne." ")write(iwr,100)"pdiag",r,s,d
      call ver_plot(s,r,d)
      if(s.ne." ")write(iwr,100)"plot",r,s,d
      call ver_ptest(s,r,d)
      if(s.ne." ")write(iwr,100)"ptest",r,s,d
_ENDIF
      call ver_rpa(s,r,d)
      if(s.ne." ")write(iwr,100)"rpa",r,s,d
_IF(rpagrad)
      call ver_rpagrad(s,r,d)
      if(s.ne." ")write(iwr,100)"rpagrad",r,s,d
_ENDIF
      call ver_scf(s,r,d)
      if(s.ne." ")write(iwr,100)"scf",r,s,d
_IF(secd_parallel)
      call ver_secd_parallel(s,r,d)
      if(s.ne." ")write(iwr,100)"secd_parallel",r,s,d
      call ver_secdchf(s,r,d)
      if(s.ne." ")write(iwr,100)"secdchf",r,s,d
_ENDIF
      call ver_tsort(s,r,d)
      if(s.ne." ")write(iwr,100)"tsort",r,s,d
      s = " "
      r = " "
      d = " "
      call ver_tsortc(s,len(s),r,len(r),d,len(d))
      if(s.ne." ")write(iwr,100)"tsortc",r,s,d
      call ver_sec1e(s,r,d)
      if(s.ne." ")write(iwr,100)"sec1e",r,s,d
      call ver_sec2e(s,r,d)
      if(s.ne." ")write(iwr,100)"sec2e",r,s,d
      call ver_secmp2(s,r,d)
      if(s.ne." ")write(iwr,100)"secmp2",r,s,d
      call ver_server(s,r,d)
      if(s.ne." ")write(iwr,100)"server",r,s,d
      call ver_tdaf(s,r,d)
      if(s.ne." ")write(iwr,100)"tdaf",r,s,d
      call ver_tran4(s,r,d)
      if(s.ne." ")write(iwr,100)"tran4",r,s,d
      call ver_util1(s,r,d)
      if(s.ne." ")write(iwr,100)"util1",r,s,d
      call ver_util2(s,r,d)
      if(s.ne." ")write(iwr,100)"util2",r,s,d
      call ver_util3(s,r,d)
      if(s.ne." ")write(iwr,100)"util3",r,s,d
      call ver_util4(s,r,d)
      if(s.ne." ")write(iwr,100)"util4",r,s,d
      call ver_util5(s,r,d)
      if(s.ne." ")write(iwr,100)"util5",r,s,d
      call ver_util6(s,r,d)
      if(s.ne." ")write(iwr,100)"util6",r,s,d
      call ver_util7(s,r,d)
      if(s.ne." ")write(iwr,100)"util7",r,s,d
      call ver_util8(s,r,d)
      if(s.ne." ")write(iwr,100)"util8",r,s,d
_IF(vdw)
      call ver_inter_vdwaals(s,r,d)
      if(s.ne." ")write(iwr,100)"vdwaals interface",r,s,d
      call ver_vdwaals(s,r,d)
      if(s.ne." ")write(iwr,100)"vdwaals",r,s,d
_ENDIF
      call ver_zora(s,r,d)
      if(s.ne." ")write(iwr,100)"zora",r,s,d
_IF(ccpdft)
      write(iwr,*)
      write(iwr,*)'CCP1 DFT Module'
      write(iwr,*)
      call ver_dft_basis(s,r,d)
      if(s.ne." ")write(iwr,100)"basis",r,s,d
      call ver_dft_chf(s,r,d)
      if(s.ne." ")write(iwr,100)"chf",r,s,d
      call ver_dft_coulomb(s,r,d)
      if(s.ne." ")write(iwr,100)"coulomb",r,s,d
      call ver_dft_deriv2e(s,r,d)
      if(s.ne." ")write(iwr,100)"deriv2e",r,s,d
      call ver_dft_exp_dksm_hess(s,r,d)
      if(s.ne." ")write(iwr,100)"exp_dksm_hess",r,s,d
_IF(ga)
      call ver_dft_gadft(s,r,d)
      if(s.ne." ")write(iwr,100)"gadft",r,s,d
_ENDIF
      call ver_dft_gamess(s,r,d)
      if(s.ne." ")write(iwr,100)"gamess",r,s,d
      call ver_dft_gden(s,r,d)
      if(s.ne." ")write(iwr,100)"gden",r,s,d
      call ver_dft_global(s,r,d)
      if(s.ne." ")write(iwr,100)"global",r,s,d
      call ver_dft_integ_data(s,r,d)
      if(s.ne." ")write(iwr,100)"integ_data",r,s,d
      call ver_dft_integ_te2c_norm(s,r,d)
      if(s.ne." ")write(iwr,100)"integ_te2c_norm",r,s,d
      call ver_dft_integ2e(s,r,d)
      if(s.ne." ")write(iwr,100)"integ2e",r,s,d
      call ver_dft_interface(s,r,d)
      if(s.ne." ")write(iwr,100)"interface",r,s,d
      call ver_dft_intpack(s,r,d)
      if(s.ne." ")write(iwr,100)"intpack",r,s,d
      call ver_dft_jfitg(s,r,d)
      if(s.ne." ")write(iwr,100)"jfitg",r,s,d
      call ver_dft_lebedevlaikov(s,r,d)
      if(s.ne." ")write(iwr,100)"lebedev-laikov",r,s,d
      call ver_dft_matform(s,r,d)
      if(s.ne." ")write(iwr,100)"matform",r,s,d
      call ver_dft_matpack(s,r,d)
      if(s.ne." ")write(iwr,100)"matpack",r,s,d
_IF(never)
      call ver_dft_mpole1(s,r,d)
      if(s.ne." ")write(iwr,100)"mpole1",r,s,d
      call ver_dft_mpole2(s,r,d)
      if(s.ne." ")write(iwr,100)"mpole2",r,s,d
_ENDIF
      call ver_dft_readbasis(s,r,d)
      if(s.ne." ")write(iwr,100)"readbasis",r,s,d
      call ver_dft_xc(s,r,d)
      if(s.ne." ")write(iwr,100)"xc",r,s,d
      call ver_dft_xc_lib(s,r,d)
      if(s.ne." ")write(iwr,100)"xc_lib",r,s,d
      call ver_dft_weights(s,r,d)
      if(s.ne." ")write(iwr,100)"weights",r,s,d
_ENDIF

      write(iwr,102)
100   format(1x,a12,2x,a14,2x,a50,2x,a20)
101   format(1x,'module',8x,'version',9x,'source',47x,'date')
102   format(1x,102('='))
      end
c


      subroutine chksiz(olog)

      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)

INCLUDE(common/sizes)
INCLUDE(common/iofile)      
INCLUDE(common/infoa)      
INCLUDE(common/nshel)      
INCLUDE(common/czmat)

      opr = olog .and. opg_root()

      ok = .true.

      if(opr)write(iwr,*)'Checking dimensions:'

      otst = maxat .ge. nat
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'atoms','maxat',maxat,nat,otst

      otst = maxorb .ge. num
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'orbitals','maxorb',maxorb,num,otst

      otst = mxshel .ge. nshell
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'shells','mxshel',mxshel,nshell,otst

      i = kstart(nshell) + kng(nshell) - 1
      otst = mxprim .ge. i
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'primitives','mxprim',mxprim,i,otst

      otst = maxnz .ge. nz
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'z-atoms','maxnz',maxnz,nz,otst

      otst = maxvar .ge. nvar
      ok = ok .and. otst
      if(opr.or..not.otst)
     &     write(iwr,100)'z-variables','maxvar',maxvar,nvar,otst

 100  format(1x,a20,a8,i8,i8,l1)

      if( .not. ok ) call caserr2('code dimensions too small')

      return 
      end

_IF(charmm)
c
c  Pass gamess-uk results back to charmm Note there are a few
c  different versions of this routine. This is for versions 4&5 of
c  the interface see www.cse.clrc.ac.uk/ccg/software/chmguk
c  for details.
c
_IF(charmm_if4)
c for version 4
      subroutine gms2chm(gtot,dx,dy,dz,map,block,natom)
_ELSE
c from version 5 we recover a quantum derived charge
      subroutine gms2chm(gtot,dx,dy,dz,qmchg,map,block,natom)
_ENDIF

      implicit none
INCLUDE(common/sizes)
INCLUDE(common/funct)
INCLUDE(common/iofile)
INCLUDE(common/chmgms)
INCLUDE(common/infoa)

_IF(charmm_if4)
      REAL gtot, DX(*), DY(*), DZ(*), block
_ELSE
      REAL gtot, DX(*), DY(*), DZ(*), block, qmchg(*)
_ENDIF
      integer map(*), natom

      integer i, ii
      REAL fac
      integer ipg_nnodes
c
c These factors have been taken from charmm v 26
c Correction for 1/nodes is empirical
c
      REAL bohrr
      parameter (bohrr = 0.529177249d0)

      REAL tokcal
      parameter (tokcal = 627.5095d0)
c
      fac = ipg_nnodes()
      fac = block / fac
      fac = fac * tokcal / bohrr
      do i = 1,nat
         ii = map(i)
         dx(ii) = dx(ii) + fac*egrad(3*(i-1)+1)
         dy(ii) = dy(ii) + fac*egrad(3*(i-1)+2)
         dz(ii) = dz(ii) + fac*egrad(3*(i-1)+3)
      enddo
      fac = ipg_nnodes()
      fac = block / fac
      fac = fac * tokcal
      gtot = gtot + fac*(enrgy - eoff)

_IF(charmm_if4)
c no charges passed in this interface
_ELSE
      do i=1,natqm
         ii = map(i)
         qmchg(ii) = qatch(i)
      enddo
_ENDIF

      return
      end
_ENDIF
      subroutine onepdm(da,db)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c     extracts the one-particle density matrix for various types
c     of wavefunction
INCLUDE(common/cigrad)
INCLUDE(common/infoa)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/iofile)
INCLUDE(common/atmblk)
c
      logical exist
      dimension da(*),db(*)
      character *8 closed,open,grhf
      data closed/'closed'/
      data  grhf/'grhf'/
      data open /'open'/
      if (lci .or. lmcscf .or. cigr .or. mcgr) then
c        mtyp = 0
         write (iwr,*) ' core orbitals ' , ncore
         call secloc(isecdd,exist,iblok)
         if (exist) then
            call rdedx(da,nx,iblok,ifild)
         else
            call caserr2('ci density matrix not found')
         end if
         if (ncore.eq.0) return
         isec = isect(8)
         if (mcgr) isec = isecmo
c        mtyp = 0
         call secloc(isec,exist,iblok)
         if (exist) then
            call rdedx(db,num*ncoorb,iblok+mvadd,ifild)
         else
            call caserr2('vectors not found')
         end if
         ij = 0
         do 40 i = 1 , num
            do 30 j = 1 , i
               ij = ij + 1
               do 20 k = 1 , ncore
                  kk = (k-1)*num
                  da(ij) = da(ij) + 2.0d0*db(kk+i)*db(kk+j)
 20            continue
 30         continue
 40      continue
         return
      else
c        mtyp = 0
         call secloc(isect(7),exist,iblok)
         if (exist) then
            call rdedx(da,lds(isect(7)),iblok,ifild)
         else
            call caserr2('closed shell density matrix not found')
         end if
         if (.not.((mp2 .or. mp3) .and. scftyp.eq.closed)) then
            if (scftyp.eq.grhf .or. scftyp.eq.closed) return
c           mtyp = 0
            call secloc(isect(10),exist,iblok)
            if (exist) then
               call rdedx(db,lds(isect(10)),iblok,ifild)
            else
               call caserr2('open shell density matrix not found')
            end if
            do 50 i = 1 , nx
               da(i) = da(i) + db(i)
 50         continue
            if (.not.((mp2 .or. mp3) .and. scftyp.eq.open)) then
               return
            end if
         end if
      end if
c
c     mtyp = 0
      call secloc(isecdd,exist,iblok)
      if (.not.exist .and. scftyp.eq.closed .and. mp2)
     +    call secloc(isect(45),exist,iblok)
      if (.not.exist) return
      call rdedx(db,nx,iblok,ifild)
c
      do 60 ms = 1 , nx
         da(ms) = da(ms) + db(ms)
 60   continue
c
      return
      end
      subroutine prnder(der2,ndim,iw)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
      character *1 clab
      dimension der2(ndim,ndim),clab(3)
      data clab/'x','y','z'/
      max = 0
 20   min = max + 1
      max = max + 3
      if (max.gt.nat) max = nat
      write (iw,6020) (n,n=min,max)
      write (iw,6030) (zaname(n),n=min,max)
      write (iw,6040) ((clab(m),m=1,3),n=min,max)
      write (iw,6010)
      j0 = 3*(min-1) + 1
      j1 = 3*max
      do 30 iat = 1 , nat
         i0 = 3*(iat-1)
         write (iw,6050) iat , zaname(iat) , clab(1) ,
     +                   (der2(i0+1,j),j=j0,j1)
         write (iw,6060) clab(2) , (der2(i0+2,j),j=j0,j1)
         write (iw,6060) clab(3) , (der2(i0+3,j),j=j0,j1)
 30   continue
      if (max.lt.nat) go to 20
      return
 6010 format (/)
 6020 format (//17x,3(10x,i3,14x))
 6030 format (/17x,3(9x,a8,10x))
 6040 format (/17x,3(4x,a1,8x,a1,8x,a1,4x))
 6050 format (i3,3x,a8,2x,a1,2x,9f9.5)
 6060 format (16x,a1,2x,9f9.5)
      end

      subroutine coperm(ip,flop,n,iwr)
c
c     complete a permutation ip   (flop is scratch space)
c
c
      implicit REAL  (a-h,o-z),integer  (i-n)
      logical flop
      dimension ip(n),flop(n)
c
      do 10 i=1,n
10    flop(i) = .false.
c
c     check what numbers are already there
c
      do 20 i=1,n
         if (ip(i).le.0) go to 20
         k = i
         if (flop(ip(i))) go to 100
         flop(ip(i)) = .true.
20    continue
c
c...  fill in missing numbers
c
      k = 0
      do 40 i=1,n
         if (ip(i).gt.0) go to 40
30       k = k + 1
         if (flop(k)) go to 30
         ip(i) = k
40    continue
c
      return
c
c...  error a number appears twice
100   write(iwr,1) k,(ip(i),i=1,n)
1     format(' ****** error ******* the number ',i5,' appears twice',
     1 ' in the following permutation ',/,
     2 (15x,20(i4),1x))
c
      call caserr2('permutation error in coperm')
      return
c
      end
cjmht etijk needed by both anald and integd, but for base build with drf
c     these files are not needed together so this function now separate
      function etijk(i,j,k)
      implicit REAL  (a-h,o-z)
      dimension e(3,3,3)
      data e/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,
     &   0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,
     &  0.0d0,-1.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
      etijk = e(i,j,k)
      return
      end
cjmht - this previously in util8. It is currently only used
c       by drf/mscon and mopac7
      subroutine imatx(a)
      implicit none
      REAL a
      integer i,j
c------
c load identity matrix
c------
      dimension a(3,3)
      do 100 i = 1,3
         do 50 j = 1,3
            a(i,j) = 0.0d0
 50       continue
         a(i,i) = 1.0d0
100   continue
      return
      end
      subroutine ver_util3(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util3.m,v $
     +     "/
      data revision /"$Revision: 6115 $"/
      data date /"$Date: 2009-12-17 11:37:45 +0100 (Thu, 17 Dec 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
