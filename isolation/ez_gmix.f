      program mixer
        implicit none       
        integer :: m, l, j, n, k, i
        integer :: mst, lst, jst, nst, kst, ist
        integer :: mdif, ldif, jdif, ndif, kdif, idif
        integer :: mit, lit, jit, nit, kit, iit
        integer :: index, it, nblock, iar
        integer :: mop, lop, jop, nop, kop, iop
        
        logical equal
        ! common /posit/ iky(3)     
        nblock = 256

        ! Start values of the loops
        ! mst = 1
        ! lst = fake_ig(4, mst) ! Do all these inits in the inner most loop to make sure we get the right values
        ! jst = fake_ig(3, mst)
        ! nst = mst + 1
        ! kst = fake_ig(4, nst)
        ! ist = fake_ig(3, nst)
        mst = 1
        lst = 3
        jst = 3
        nst = 3
        kst = 4
        ist = 4


        ! Stop conditions
        ! mop = nblock-1
        ! lop = fake_ig(4,m+1)-1
        ! jop = fake_ig(3,m+1)-1
        ! nop = nblock
        ! kop = fake_ig(4,n)+fake_ig(2,n)-1
        ! iop = fake_ig(3,n)+fake_ig(1,n)-1
        mop = mst + 2
        lop = lst + 1 
        jop = jst + 1 
        nop = nst + 1
        kop = kst + 3 
        iop = ist + 3 
        

        ! ! Number of iterations each loop does : stop - start + 1

        ! print *, "diffs:", mdif, ldif, jdif, ndif, kdif, idif

        ! Current iteration counter for every loop
        ! mit = m - mst
        ! lit = l - lst
        ! jit = j - jst
        ! nit = n - nst
        ! kit = k - kst
        ! iit = i - ist
        print *, "Start the loops" 
        it = 0
        do 60 m=mst, mop
          do 50 l=lst, lop 
            do 40 j=jst, jop !// 1 2//s20
              do 30 n=nst, nop          !1 2 3    
                do 20 k=kst, kop  !1 2 3 4 5 6
!$acc parallel loop
                  do 10 i=ist, iop !1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
                    it = it + 1

                    mdif = mop + 1 - mst
                    ldif = lop + 1 - lst
                    jdif = jop + 1 - jst
                    ndif = nop + 1 - nst
                    kdif = kop + 1 - kst
                    idif = iop + 1 - ist

                    ! current its
                    mit = m - mst
                    lit = l - lst
                    jit = j - jst
                    nit = n - nst
                    kit = k - kst
                    iit = i - ist

                    ! print *, n, nst, nit, nop
                    ! print *, "iteration variables"
                    ! print *, m, l, j, n, k, i
                    ! print *, "start values"
                    ! print *, mst, lst, jst, nst, kst, ist
                    ! print *, "It counters"
                    ! print *, mit, lit, jit, nit, kit, iit
                    ! index = (iit + 1) + kit * idif + nit * kdif * idif + jit * ndif * kdif * idif + lit * jdif * ndif * kdif * idif + mit * ldif * jdif * ndif * kdif * idif
                    ! index = (iit + 1) +kit*idif+nit*kdif*idif
                    ! index = index + jit * ndif * kdif * idif
                    ! index = index + lit * jdif * ndif *kdif*idif
                    ! index = index +  mit*ldif*jdif*ndif*kdif*idif

                    index = (iit + 1) + kit * idif
                    index = index + nit * kdif * idif  
                    index = index + jit * ndif * kdif * idif
                    index = index + lit * jdif * ndif * kdif * idif
                    iar = mit * ldif * jdif * ndif * kdif * idif
                    index = index + iar
                    print *, "it = ", it
                    print *, "index = ", index
                    print *, ""
                    nst = 2
10                continue
!$acc end parallel loop
20              continue
30            continue
40          continue
50        continue
60      continue

        print *, "Made it tkankerhrough" 
        ! print *, "it = ", it
        print *, "index = ", index
        contains
        integer function fake_ig(val, val2) result(out)
          implicit none
          integer, intent(in) :: val, val2
          out = 20 / val * val2
        end function fake_ig 
      end program mixer