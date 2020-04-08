      subroutine bmove (a,b,n)
      implicit real*8 (a-h,p-w), integer (i-n), logical (o)
      dimension a(*),b(*)
      do 10 i=n,1,-1
10    b(i) = a(i)
      return
      end
