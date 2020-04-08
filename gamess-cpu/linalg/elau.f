      subroutine elau(hinv,l,d,a,e)
c
      double precision a(*)
      double precision d(l)
      double precision e(l)
      double precision f
      double precision g
      double precision half
      double precision hh
      double precision hinv
      double precision zero
c
      parameter (zero = 0.0d+00, half = 0.5d+00)
c
      jl = l
      e(1) = a(1) * d(1)
      jk = 2
      do 210 j = 2, jl
         f = d(j)
         g = zero
         jm1 = j - 1
c
         do 200 k = 1, jm1
            g = g + a(jk) * d(k)
            e(k) = e(k) + a(jk) * f
            jk = jk + 1
  200    continue
c
         e(j) = g + a(jk) * f
         jk = jk + 1
  210 continue
c
c        .......... form p ..........
c
      f = zero
      do 245 j = 1, l
         e(j) = e(j) * hinv
         f = f + e(j) * d(j)
  245 continue
c
c     .......... form q ..........
c
      hh = f * half * hinv
      do 250 j = 1, l
  250 e(j) = e(j) - hh * d(j)
c
      return
      end
