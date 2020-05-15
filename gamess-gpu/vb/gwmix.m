      subroutine gwmix(ir,ic,ig,nblock,ialfa,w1,supg,nelec,n1,val)
c     Combines pikjl, wmix, gather and gmix
c     
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
INCLUDE(common/tractlt)
      dimension w1(n1)
      dimension ig(5,nblock)
      dimension ir(nelec),ic(nelec),supg(*)
!$acc routine (intpos)
c  supg(0:n2int)
c
      val=0.0d0
c     print *,'remco ig::',nblock,ialfa
c     do i=1,nblock
c       print *,(ig(k,i),k=1,5)
c     enddo

!$acc data copyin(ir,ic) present(supg)
!$acc& copyin(ig(5,nblock),w1(1:n1),nblock)
!$acc& copy(val)
!$acc kernels present (supg)
! parallel loops
! reduction val

      do 51 m=1,nblock
         msta = ig(3,m)
         mend = ig(3,m) + ig(1,m) - 1
         do 41 k=msta+1,mend
            do 31 l=msta,k-1
              do 21 i=msta+1,mend
                do 11 j=msta,i-1
                  !cikjl: calculate the loop vars
                  ii=i-msta+1
                  jj=j-msta+1
                  kk=k-msta+1
                  ll=l-msta+1
                  !jacobi ratio theorem to calculate 2nd-o-cofacs
                  scal=w1((kk-1)*ig(1,m)+ii+ig(5,m)-1)*
     &                 w1((ll-1)*ig(1,m)+jj+ig(5,m)-1)-
     &                 w1((ll-1)*ig(1,m)+ii+ig(5,m)-1)*
     &                 w1((kk-1)*ig(1,m)+jj+ig(5,m)-1)
                  !(subvec ipos with ipose) * 2nd-o-cofactor, add to total
                  ! matrix element val
                  val=val+scal*(supg(intpos(ir(i),ic(k),ir(j),
     1            ic(l)))-supg(intpos(ir(i),ic(l),ir(j),ic(k))))
                  !it=it+1 !why is this still here? leftover
11                continue
21             continue
31         continue
41       continue
51    continue
c     it works yeey
      do 690 m=1,nblock-1
        do 590 l=ig(4,m),ig(4,m+1)-1
          do 490 j=ig(3,m),ig(3,m+1)-1
            ! wmix: get the first order cofactor from w1
            jl=(l-ig(4,m))*(ig(3,m+1)-ig(3,m))+j-ig(3,m)+ig(5,m)
               scalar=w1(jl)
            do 390 n=m+1,nblock
              do 290 k=ig(4,n),ig(4,n)+ig(2,n)-1
                do 190 i=ig(3,n),ig(3,n)+ig(1,n)-1
                  ! get first order cofactor from w1
                  ik=(k-ig(4,n))*(ig(3,n)+ig(1,n)-ig(3,n))+i
     1                -ig(3,n)+ig(5,n)
                  ! calc 2nd order cofactor * integral calc
                  ! add to total
                  if (m >= ialfa+1) then
                    ! print *, "third"
                    val=val - w1(ik)*scalar *
     3                  supg(intpos(ir(i),ic(l),ir(j),ic(k)))
                  end if

                  if (m <= ialfa-1 .and. n <= ialfa) then
                    ! print *, "second"
                    val=val - w1(ik)*scalar *
     2                  supg(intpos(ir(i),ic(l),ir(j),ic(k)))
                  end if
                    val=val+w1(ik)*scalar*
     1                  supg(intpos(ir(i),ic(k),ir(j),ic(l)))
190               continue
290             continue
390           continue
490         continue
590       continue
690     continue
      ! print *, "ending val concat", val
      ! end if
!$acc end kernels
!$acc end data
      ! stop
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