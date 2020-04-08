      subroutine dsctr(n,a,map,r)
      implicit double precision  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),r(*),map(*)
      if(n.gt.0) then
      do 1 loop=1,n
   1  r(map(loop))=a(loop)
      endif
      return
      end
