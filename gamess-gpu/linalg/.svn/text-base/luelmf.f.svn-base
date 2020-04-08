      subroutine luelmf (a,b,ipvt,n,ia,x)
c
c-----------------------------------------------------------------------
c
c   purpose             - elimination part of solution of ax=b
c                           (full storage mode)
c
c   usage               - call luelmf (a,b,ipvt,n,ia,x)
c
c   arguments    a      - a = lu (the result computed in the
c                           routine ludatf) where l is a lower
c                           triangular matrix with ones on the main
c                           diagonal. u is upper triangular. l and u
c                           are stored as a single matrix a and the
c                           unit diagonal of l is not stored. (input)
c                b      - b is a vector of length n on the right hand
c                           side of the equation ax=b. (input)
c                ipvt   - the permutation matrix returned from the
c                           routine ludatf, stored as an n length
c                           vector. (input)
c                n      - order of a and number of rows in b. (input)
c                ia     - row dimension of a exactly as specified in
c                           the dimension statement in the calling
c                           program. (input)
c                x      - the result x. (output)
c
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through routine uhelp
c
c-----------------------------------------------------------------------
c
c
      dimension          a(ia,*),b(*),ipvt(*),x(*)
c      REAL  a,b,x,sum
      double precision  a,b,x,sum
c                                  first executable statement
c                                  solve ly = b for y
      do 5 i=1,n
    5 x(i) = b(i)
      iw = 0
      do 20 i=1,n
         ip = ipvt(i)
         sum = x(ip)
         x(ip) = x(i)
         if (iw .eq. 0) go to 15
         im1 = i-1
         do 10 j=iw,im1
            sum = sum-a(i,j)*x(j)
   10    continue
         go to 20
   15    if (sum .ne. 0.d0) iw = i
   20 x(i) = sum
c                                  solve ux = y for x
      do 30 ib=1,n
         i = n+1-ib
         ip1 = i+1
         sum = x(i)
         if (ip1 .gt. n) go to 30
         do 25 j=ip1,n
            sum = sum-a(i,j)*x(j)
   25   continue
   30 x(i) = sum/a(i,i)
      return
      end
