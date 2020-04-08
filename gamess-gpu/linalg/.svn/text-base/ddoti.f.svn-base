      function ddoti(n,sx,isy,sy)
c
c     returns the dot product of real sx and sy(isy)
c
      implicit double precision  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension sx(*),sy(*),isy(*)
      ddoti = 0.d0
      if(n.le.0)return
      do  10 loop = 1,n
         ddoti = ddoti + sx(loop)*sy(isy(loop))
   10 continue
      return
      end
