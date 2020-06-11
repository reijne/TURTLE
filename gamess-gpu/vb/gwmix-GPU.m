      subroutine gwmix(ir,ic,ig,nblock,ialfa,w1,supg,nelec,n1,sum)
c     Combines cikjl, pikjl, wmix, gmix, (gather, ddot and subvec)
c     from matre3.
c     ======================NOLIFT + ACC=========================
      implicit REAL  (a-h,o-z) , integer   (i-n)
      integer count_0, count_1, count_rate, count_max
c
INCLUDE(common/tractlt)
      dimension w1(n1)
      dimension ig(5,nblock)
      dimension ir(nelec),ic(nelec),supg(*)
!$acc routine(intpos)

      sum=0.0d0
!$acc data copyin(ir,ic) present(supg), copy(sum)
!$acc& copyin(ig(5,nblock), w1(1:n1), nblock)
      do 51 m=1,nblock
          do 41 k=ig(3,m)+1, ig(3,m)+ig(1,m)-1
!$acc parallel loop reduction(+:sum)
          do 31 l=ig(3,m), k-1
!$acc loop reduction(+:sum)
              do 21 i=ig(3,m)+1, ig(3,m)+ig(1,m)-1
!$acc loop reduction(+:sum) private(ii, jj, kk, ll, ig5min)
                do 11 j=ig(3,m), i-1
                  !cikjl: calculate the loop vars
                  ii=i-ig(3,m)
                  jj=j-ig(3,m)
                  kk=k-ig(3,m)
                  ll=l-ig(3,m)

                  !jacobi ratio theorem to calculate 2nd-o-cofacs
                  soo=w1((kk)*ig(1,m)+ii+ig(5,m))*
     &                w1((ll)*ig(1,m)+jj+ig(5,m))-
     &                w1((ll)*ig(1,m)+ii+ig(5,m))*
     &                w1((kk)*ig(1,m)+jj+ig(5,m))
                  !(subvec ipos with ipose) * 2nd-o-cofactor, add to total
                  ! matrix element sum
                  sum=sum+soo*(supg(intpos(ir(i),ic(k),ir(j),ic(l)))
     1                -supg(intpos(ir(i),ic(l),ir(j),ic(k))))
11                continue
21             continue
31         continue
!$acc end parallel loop
41       continue
51    continue

      do 690 m=1, nblock-1
        do 590 l=ig(4,m), ig(4,m+1)-1
          do 490 j=ig(3,m), ig(3,m+1)-1
!$acc parallel loop reduction(+:sum)
            do 390 n=m+1, nblock
!$acc loop reduction(+:sum)
              do 290 k=ig(4,n), ig(4,n)+ig(2,n)-1
!$acc loop reduction(+:sum) private(jl, ik, soo)
                do 190 i=ig(3,n), ig(3,n)+ig(1,n)-1
                  ! wmix: get the first order cofactor from w1
                  jl=(l-ig(4,m))*(ig(3,m+1)-ig(3,m))+j-ig(3,m)+ig(5,m)
                  ! get first order cofactor from w1
                  ik=(k-ig(4,n))*(ig(3,n)+ig(1,n)-ig(3,n))+i
     1                -ig(3,n)+ig(5,n)
                  soo = w1(ik) * w1(jl)
                  ! calc 2nd order cofactor * integral calc
                  ! add to total
                  if (m>=ialfa+1.or.m<=ialfa-1.and.n<=ialfa) then
                    sum=sum-soo*supg(intpos(ir(i),ic(l),ir(j),ic(k)))
                  end if
                    sum=sum+soo*supg(intpos(ir(i),ic(k),ir(j),ic(l)))
190               continue
290             continue
390           continue
!$acc end parallel loop
490         continue
590       continue
690     continue

!$acc end data
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
      return
      end
        
      function intposx(i,j,k,l)
      implicit REAL (a-h,o-z)
!$acc routine
      intposx=0
      return
      end