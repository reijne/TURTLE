      subroutine gwmix(ir,ic,ig,nblock,ialfa,w1,supg,nelec,n1,val)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
INCLUDE(common/tractlt)
      dimension w1(n1)
      dimension ig(5,nblock)
      dimension ir(nelec),ic(nelec),supg(0:n2int)
      integer intpo
!$acc routine (intpo)
c
      val=0.0d0
c     print *,'remco ig::',nblock,ialfa
c     do i=1,nblock
c       print *,(ig(k,i),k=1,5)
c     enddo
!$acc data copyin(ir,ic) present(supg)
!$acc& copyin(ig(5,nblock),w1(1:n1),nblock) 
!$acc& copy(val)
!$acc parallel loop private(msta,mend) reduction(+:val) present (supg) 
!$acc& gang vector_length(16)
      do 51 m=1,nblock
         msta = ig(3,m)
         mend = ig(3,m) + ig(1,m) - 1
!$acc loop
         do 41 k=msta+1,mend
            do 31 l=msta,k-1
!$acc loop
               do 21 i=msta+1,mend
!$acc loop seq reduction(+:val)
                  do 11 j=msta,i-1
                     val=val+w1(i)*w1(j)*(supg(intpo(ir(i),ic(k),ir(j),
     1            ic(l)))-supg(intpo(ir(i),ic(l),ir(j),ic(k))))
11                continue
21             continue
31         continue
41       continue
51    continue
c
!$acc parallel loop reduction(+:val) gang vector_length(16) 
!$acc& present (supg,ir,ic,ig,w1) 
      do 60 m=1,nblock-1
!$acc loop seq
               do 30 n=m+1,nblock
!$acc loop 
         do 50 l=ig(4,m),ig(4,m+1)-1
!$acc loop 
            do 40 j=ig(3,m),ig(3,m+1)-1
!$acc loop 
                  do 20 k=ig(4,n),ig(4,n)+ig(2,n)-1
!$acc loop seq reduction(+:val)
                     do 10 i=ig(3,n),ig(3,n)+ig(1,n)-1
                         val=val+w1(i)*w1(j)*
     1                      supg(intpo(ir(i),ic(k),ir(j),ic(l)))
10                   continue
20                continue
40          continue
50       continue
30             continue
60    continue
c
!$acc parallel loop reduction(+:val) present (supg,ir,ic,ig,w1) 
!$acc& gang vector_length(16) 
      do 120 m=1,ialfa-1
               do 90 n=m+1,ialfa
         do 110 l=ig(4,m),ig(4,m+1)-1
            do 100 j=ig(3,m),ig(3,m+1)-1
                  do 80 k=ig(4,n),ig(4,n)+ig(2,n)-1
!$acc loop seq  reduction(+:val)
                     do 70 i=ig(3,n),ig(3,n)+ig(1,n)-1
                         val=val-w1(i)*w1(j)*
     1                      supg(intpo(ir(i),ic(l),ir(j),ic(k)))
70                   continue
80                continue
100         continue
110      continue
90             continue
120   continue
!$acc parallel loop reduction(+:val) present (supg,ir,ic,ig,w1)
!$acc& gang vector_length(16) 
      do 180 m=ialfa+1,nblock-1
               do 150 n=m+1,nblock
         do 170 l=ig(4,m),ig(4,m+1)-1
            do 160 j=ig(3,m),ig(3,m+1)-1
                  do 140 k=ig(4,n),ig(4,n)+ig(2,n)-1
!$acc loop  seq  reduction(+:val)
                     do 130 i=ig(3,n),ig(3,n)+ig(1,n)-1
                         val=val-w1(i)*w1(j)*
     1                      supg(intpo(ir(i),ic(l),ir(j),ic(k)))
130                  continue
140               continue
160         continue
170      continue
150            continue
180   continue
cc!$acc update self (val)
cc!$acc end kernels
!$acc end data
      print *,'remco val::',val

      return
      end
   
       function intpo(i,j,k,l)
!$acc routine
       intpo=0
        return
       end
      subroutine delintdev(supg)
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(common/tractlt)
      dimension supg(0:n2int)
!$acc exit data delete (supg)
      return
      end
    
      subroutine putintdev(supg)
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(common/tractlt)
      dimension supg(0:n2int)
!$acc enter data copyin(supg)
c     return
      end
       
