      SUBROUTINE CSORT(i,j,k,l)
C                                                                               
C DIESE SUBROUTINE Sortiert die indexe der INTEGRALE                           
C in canonische Ordnung
C
C     i,j,k,l : Indexe des Integrals
C                                                                               
      implicit none
C 
      integer*4 i,j,k,l
      integer*4 i0,j0,k0,l0
      integer*4 mij,mkl      
      mij = max(i,j)
      mkl = max(k,l)
       If (i.lt.j) then
         write(6,*) 'Sort needed'
         i0 = i
         j0 = j
         i  = j0
         j  = i0
       End If
       If (k.lt.l) then
         write(6,*) 'Sort needed'
         k0 = k
         l0 = l
         k  = l0
         l  = k0
       End If
       If (i.lt.k) then
         write(6,*) 'Sort needed'
         i0 = i
         j0 = j
         k0 = k
         l0 = l
         i  = k0
         j  = l0
         k  = i0
         l  = j0
       Else If (i.eq.k) then
        If (j.lt.l) then
          write(6,*) 'Sort needed'
            i0 = i
            j0 = j
            k0 = k
            l0 = l
            i  = k0
            j  = l0
            k  = i0
            l  = j0
         End If
        End If
C
      RETURN                                                                    
      END                                                                       
