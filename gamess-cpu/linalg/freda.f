      subroutine freda(l,d,a,e)
c
      double precision a(*)
      double precision d(l)
      double precision e(l)
      double precision f
      double precision g
c
      jk = 1
c
c     .......... form reduced a ..........
c
      do 280 j = 1, l
         f = d(j)
         g = e(j)
c
         do 260 k = 1, j
            a(jk) = a(jk) - f * e(k) - g * d(k)
            jk = jk + 1
  260    continue
c
  280 continue
      return
      end
