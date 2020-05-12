      program mixer
        implicit none       
        integer :: m, l, j, n, k, i
        integer :: mst, lst, jst, nst, kst, ist
        integer :: mop, lop, jop, nop, kop, iop
        integer :: it, nblock
        integer :: ikjl, total
        integer ipos(194481)
        integer ipose(194481)
        logical equal

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
        ! print *, mst, lst, jst, nst, kst, ist
        ! print *, mop, lop, jop, nop, kop, iop

        do total=1, total
          it = 0
          do 60 m=1, nblock - 1
            do 50 l=lst, lop 
              do 40 j=jst, jop 
                do 30 n=m+1, nblock          
                  do 20 k=kst, kop  
                    do 10 i=ist, iop 
                      it = it + 1
                      ikjl = intpos(i, k, j, l)
                      ipos(it) = i+j+l+k
10                  continue
20                continue
30              continue
40            continue
50          continue
60        continue
          it = 0
          do 160 m=1, nblock - 1
            do 150 l=lst, lop 
              do 140 j=jst, jop 
                do 130 n=m+1, nblock          
                  do 120 k=kst, kop  
                    do 110 i=ist, iop 
                      it = it + 1
                      ikjl = intpos(i, l, j, k)
                      ipose(it) = i+j+l+k
110                  continue
120                continue
130              continue
140            continue
150          continue
160        continue
        end do

        print *, "Made it through", total 
        do i= 1, 20
          print *, ipos(i), ipose(i)
        end do

        contains 
        integer function intpos(ll, jj, kk, ii)
        integer :: ll, jj, kk, ii
        integer :: ijkl
            ijkl = ii + kk
            ijkl = ijkl + 2 * 50 + ll
            ijkl = ijkl + kk * ii * 20
            if (ll <= kk) then
              ijkl = ijkl + 22 * kk * ii * ll
            else 
              ijkl = ijkl + 21 * kk * ii * ll
            end if
            ijkl = ijkl + jj + ii * ll + kk
            return
        end function
      end program mixer