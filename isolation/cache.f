      program mixer
        implicit none       
        integer :: m, l, j, n, k, i, x
        integer :: mst, lst, jst, nst, kst, ist
        integer :: mop, lop, jop, nop, kop, iop
        integer :: it, nblock
        integer :: ikjl, total
        integer cache(194481, 2, 4)
        integer ipos(194481)
        integer ipose(194481)
        logical equal

        cache(1, 1, 1) = -1
        total = 1000000

        nblock = 2
        mst = 1
        lst = 1
        jst = 1
        nst = 2
        kst = 22
        ist = 22

        mop = mst + nblock - 1
        lop = lst + 20 
        jop = jst + 20 
        nop = nst + nblock - 1 
        kop = kst + 20 
        iop = ist + 20 
        do total=1, total
          if (cache(1, 1, 1) == -1) then 
            it = 0
            do 60 m=1, nblock - 1
              do 50 l=lst, lop 
                do 40 j=jst, jop 
                  do 30 n=m+1, nblock          
                    do 20 k=kst, kop  
                      do 10 i=ist, iop 
                        it = it + 1
                        ikjl = intpos(i, k, j, l)
                        cache(it, 1, 1) = i
                        cache(it, 1, 2) = k
                        cache(it, 1, 3) = j
                        cache(it, 1, 4) = l !(1144, _)
10                continue
20              continue
30            continue
40          continue
50        continue
60      continue
        it = 0
        do 160 m=1, nblock - 1
          do 150 l=lst, lop 
            do 140 j=jst, jop 
              do 130 n=m+1, nblock          
                do 120 k=kst, kop  
                  do 110 i=ist, iop 
                    it = it + 1
                    ikjl = intpos(i, l, j, k) 
                    cache(it, 2, 1) = i
                    cache(it, 2, 2) = k
                    cache(it, 2, 3) = j
                    cache(it, 2, 4) = l ! (1144, 2255)
110                continue
120              continue
130            continue
140          continue
150        continue
160      continue
          else
!$acc parallel loop
            do x=1, it
              i = cache(x, 1, 1)
              k = cache(x, 1, 2)
              j = cache(x, 1, 3)
              l = cache(x, 1, 4)
              ikjl = intpos(i, k, j, l)
              ipos(x) = i+j+k+l 
              i = cache(x, 2, 1)
              k = cache(x, 2, 2)
              j = cache(x, 1, 3)
              l = cache(x, 1, 4)
              ikjl = intpos(i, l, j, k)
              ipose(x) = i+j+k+l
            end do
!$acc end parallel loop
          end if
        end do

      print *, "Made it through", total
      do i= 1, 20
        print *, ipos(i), ipose(i)
      end do
      ! print *, 1000000*i+10000*k+100*j+l
      contains
      integer function intpos(ll, jj, kk, ii)
!$acc routine
      integer :: ll, jj, kk, ii, ijkl
          ijkl = ii + kk
          ijkl = ijkl + 2 * 50 + ll
          ijkl = ijkl + kk * ii * 20
          if (l <= kk) then
            ijkl = ijkl + 22 * kk * ii * ll
          else 
            ijkl = ijkl + 21 * kk * ii * ll
          end if
          ijkl = jj + ii * ll + kk
          return
      end function
      end program mixer