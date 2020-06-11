      subroutine gwmix(ir,ic,ig,nblock,ialfa,w1,supg,nelec,n1,sum)
c     Combines cikjl, pikjl, wmix, gmix, (gather, ddot and subvec)
c     from matre3.
c     ======================CUSTOM 2nd LOOPS======================
      implicit REAL  (a-h,o-z) , integer   (i-n)
!$acc routine(intpos)
c
INCLUDE(common/tractlt)
      dimension w1(n1)
      dimension ig(5,nblock)
      dimension ir(nelec),ic(nelec),supg(*)
      integer jl, end, le
      integer count_0, count_1, count_rate, count_max
c  supg(0:n2int)
c       
      sum=0.0d0
      ! if (ig(4, 2) == 22) then
      !   print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      ! print *, nblock, ialfa
      ! call system_clock(count_0, count_rate, count_max)
      ! end if
c     enddo
!$acc data copyin(ir,ic) present(supg), copy(sum)
!$acc& copyin(ig(5,nblock), w1(1:n1), nblock)
      do 51 m=1,nblock
          !msta = ig(3,m)
          !mend = ig(3,m) + ig(1,m) - 1
          do 41 k=ig(3,m)+1, ig(3,m)+ig(1,m)-1
!$acc parallel loop reduction(+:sum)
          do 31 l=ig(3,m), k-1
!$acc loop reduction(+:sum)
              do 21 i=ig(3,m)+1, ig(3,m)+ig(1,m)-1
!$acc loop reduction(+:sum) private(ii, jj, kk, ll, ig5min)
                do 11 j=ig(3,m), i-1
                  !cikjl: calculate the loop vars
                  ! ig3mp = ig(3,m)+1 
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
                  !it=it+1 !why is this still here? leftover
11                continue
21             continue
31         continue
!$acc end parallel loop
41       continue
51    continue
      ! sum_af = sum
      if (nblock==2 .and. ig(3,1)==1 .and. ig(4,1)) then
        ! print *, "custom"
        le = ig(4,2)-1
        end = (ig(3,2)-1)*(le)
!$acc parallel loop reduction(+:sum) copyin(end, le) private(l, j)
        do 1 oe=1, end
          l = modulo(oe-1, le)+1
          if (l==0) then
            l = le
          end if
          j = ((oe-1)/le) + 1
!$acc loop reduction(+:sum) private(j, l)
              do 296 k=ig(4,2), ig(4,2)+ig(2,2)-1
!$acc loop reduction(+:sum) private(jl, ik, soo, j, l)
                do 196 i=ig(3,2), ig(3,2)+ig(1,2)-1
                  ! wmix: get the first order cofactor from w1
                  jl=(l-ig(4,1))*(ig(3,2)-ig(3,1))+j-ig(3,1)+ig(5,1)
                  ! get first order cofactor from w1
                  ik=(k-ig(4,2))*(ig(3,2)+ig(1,2)-ig(3,2))+i
     1                -ig(3,2)+ig(5,2)
                  soo = w1(ik) * w1(jl)
                  ! calc 2nd order cofactor * integral calc
                  ! add to total
                  if (1>=ialfa+1.or.1<=ialfa-1.and.2<=ialfa) then
                    sum=sum-soo*supg(intpos(ir(i),ic(l),ir(j),ic(k)))
                  end if
                    sum=sum+soo*supg(intpos(ir(i),ic(k),ir(j),ic(l)))
196               continue
296             continue
1         continue
!$acc end parallel loop
      else
      do 690 m=1, nblock-1
        ! print *, "m:", m
        do 490 j=ig(3,m), ig(3,m+1)-1
          do 590 l=ig(4,m), ig(4,m+1)-1
!$acc parallel loop reduction(+:sum)
            do 390 n=m+1, nblock
              ! print *, "n:", n
!$acc loop reduction(+:sum)
              do 290 k=ig(4,n), ig(4,n)+ig(2,n)-1
!$acc loop reduction(+:sum) private(jl, ik, soo)
                do 190 i=ig(3,n), ig(3,n)+ig(1,n)-1
                  ! print *, l*100000000+j*1000000+k*10000+i
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
590       continue
490         continue
690     continue
      end if
!$acc end data
      ! print *, sum
      ! stop
      ! if (ig(4, 2) == 22) then
      !   call system_clock(count_1, count_rate, count_max)
      !   print *, count_1*1.0/count_rate
      !   stop
      ! end if
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
c     line 54
!       do 60 m=1,nblock-1
!         print *, "first loop structure"
!         do 50 l=ig(4,m),ig(4,m+1)-1
!             do 40 j=ig(3,m),ig(3,m+1)-1
!               ! wmix: get the first order cofactor from w1
!               jl=(l-ig(4,m))*(ig(3,m+1)-ig(3,m))+j-ig(3,m)+ig(5,m)
!               scalar=w1(jl)
!               do 30 n=m+1,nblock
!                 do 20 k=ig(4,n),ig(4,n)+ig(2,n)-1
!                   do 10 i=ig(3,n),ig(3,n)+ig(1,n)-1
!                     ! get first order cofactor from w1
!                     ik=(k-ig(4,n))*(ig(3,n)+ig(1,n)-ig(3,n))+i
!      1                  -ig(3,n)+ig(5,n)
!                     ! calc 2nd order cofactor * integral calc
!                     ! add to total
!                     val=val+w1(ik)*scalar*
!      1                  supg(intpos(ir(i),ic(k),ir(j),ic(l)))
! 10                continue
! 20              continue
! 30            continue
! 40          continue
! 50       continue
! 60    continue
! c
!       do 120 m=1,ialfa-1
!         print *, "second loop structure"
!         do 110 l=ig(4,m),ig(4,m+1)-1
!             do 100 j=ig(3,m),ig(3,m+1)-1
!               ! wmix: get the first order cofactor from w1
!               jl=(l-ig(4,m))*(ig(3,m+1)-ig(3,m))+j-ig(3,m)+ig(5,m)
!               scalar=w1(jl)
!               do 90 n=m+1,ialfa
!                 do 80 k=ig(4,n),ig(4,n)+ig(2,n)-1
!                   do 70 i=ig(3,n),ig(3,n)+ig(1,n)-1
!                       ! get first order cofactor from w1
!                       ik=(k-ig(4,n))*(ig(3,n)+ig(1,n)-ig(3,n))+i
!      2                    -ig(3,n)+ig(5,n)
!                       ! calc 2nd order cofactor * integral calc
!                       ! add to total
!                       val=val-w1(ik)*scalar*
!      2                    supg(intpos(ir(i),ic(l),ir(j),ic(k)))
! 70                   continue
! 80                continue
! 90             continue
! 100         continue
! 110      continue
! 120   continue
!       do 180 m=ialfa+1,nblock-1
!         print *, "third loop structure"
!          do 170 l=ig(4,m),ig(4,m+1)-1
!             do 160 j=ig(3,m),ig(3,m+1)-1
!               ! wmix: get the first order cofactor from w1
!               jl=(l-ig(4,m))*(ig(3,m+1)-ig(3,m))+j-ig(3,m)+ig(5,m)
!               scalar=w1(jl)
!                 do 150 n=m+1,nblock
!                   do 140 k=ig(4,n),ig(4,n)+ig(2,n)-1
!                      do 130 i=ig(3,n),ig(3,n)+ig(1,n)-1
!                         ! get first order cofactor from w1
!                         ik=(k-ig(4,n))*(ig(3,n)+ig(1,n)-ig(3,n))+i
!      3                     -ig(3,n)+ig(5,n)
!                         ! calc 2nd order cofactor * integral calc
!                         ! add to total
!                         val=val-w1(ik)*scalar*
!      3                      supg(intpos(ir(i),ic(l),ir(j),ic(k)))
! 130                  continue
! 140               continue
! 150            continue
! 160         continue
! 170      continue
! 180   continue
