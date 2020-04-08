      function idamin(n,sx,incx)
      implicit double precision  (a-h,o-z)
c
c   this function finds the smallest index of an element of minimum
c   magnitude of a real array sx whose n elements are stored
c   sequentially with spacing incx >= 1.  if n <= 0, the value zero
c   is returned.  thus, if i = isamin(n,sx,1), then sx(i) is an
c   element of array sx of minimum absolute value.
c
c   description of parameters
c
c    --input--
c        n  number of elements in input vector
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c
c    --output--
c   isamin  smallest index (zero if n .le. 0)
c
      dimension sx(*)
      idamin = 0
      if( n .lt. 1 ) return
      idamin = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smin = dabs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(sx(ix)).ge.smin) go to 5
         idamin = i
         smin = dabs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smin = dabs(sx(1))
      do 30 i = 2,n
         if(dabs(sx(i)).ge.smin) go to 30
         idamin = i
         smin = dabs(sx(i))
   30 continue
      return
      end
