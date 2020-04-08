      function idmin(n,sx,incx)
c
c   this function finds the smallest index of a minimum element of a
c   real array sx whose n elements are stored sequentially with
c   spacing incx >= 1.  if n <= 0, the value zero is returned.
c   thus, if i = ismin(n,sx,1), then sx(i) is an element of array sx
c   of minimum value.
c
c   description of parameters
c
c    --input--
c        n  number of elements in input vector
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c
c    --output--
c    idmin  smallest index (zero if n .le. 0)
c
      implicit double precision  (a-h,o-z)
      dimension sx(*)
      idmin = 0
      if( n .lt. 1 ) return
      idmin = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smin = (sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if((sx(ix)).ge.smin) go to 5
         idmin = i
         smin = (sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smin = (sx(1))
      do 30 i = 2,n
         if((sx(i)).ge.smin) go to 30
         idmin = i
         smin = (sx(i))
   30 continue
      return
      end
