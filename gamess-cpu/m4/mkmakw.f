      subroutine mkmakw(y,z,w,ibly,iblz,iblw,a1,a2,nocca,
     1   nvirta,ncoorb,ifils,e,ifint1,ifint2,ifile1)
c
      implicit real*8  (a-h,o-z)
      dimension y(ncoorb,ncoorb),w(ncoorb,ncoorb),
     1  z(nocca*nvirta),a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),
     2 e(ncoorb)
c
      nsq = ncoorb*ncoorb
      ntri = ncoorb*(ncoorb+1)/2
      call rdedx(z,nocca*nvirta,iblz,ifils)
      call rdedx(y,nsq,ibly,ifile1)
      do 30 i = 1 , nocca
         do 20 ia = nocca + 1 , ncoorb
            iai = (ia-nocca-1)*nocca + i
            y(ia,i) = -z(iai)
 20      continue
 30   continue
      call wrt3(y,nsq,ibly,ifile1)
c
c
      call rdedx(a1,nsq,iblw,ifile1)
      do 50 is = 1 , ncoorb
         do 40 ir = 1 , ncoorb
            w(ir,is) = y(ir,is)*e(is)
 40      continue
 50   continue
      do 70 j = 1 , nocca
         do 60 i = 1 , nocca
            w(i,j) = w(i,j) - 0.5d0*a1(i,j)
 60      continue
 70   continue
      do 90 j = nocca + 1 , ncoorb
         do 80 i = nocca + 1 , ncoorb
            w(i,j) = w(i,j) - 0.5d0*a1(i,j)
 80      continue
 90   continue
      do 110 ia = nocca + 1 , ncoorb
         do 100 i = 1 , nocca
            w(i,ia) = w(i,ia) - a1(i,ia)
 100     continue
 110  continue
c
c
      do 150 j = 1 , nocca
         do 140 k = 1 , j
            call rdedz(a2,ntri,ifint1)
            call squr(a2,a1,ncoorb)
            call rdedz(a2,nsq,ifint2)
            do 130 iq = 1 , ncoorb
               do 120 ip = 1 , ncoorb
                  a1(ip,iq) = 4.0d0*a1(ip,iq) - a2(ip,iq) - a2(iq,ip)
 120           continue
 130        continue
            temp = ddot(ncoorb*ncoorb,y,1,a1,1)
            w(j,k) = w(j,k) + 0.5d0*temp
            if (j.ne.k) then
               w(k,j) = w(k,j) + 0.5d0*temp
            end if
 140     continue
 150  continue
c
      call wrt3(w,nsq,iblw,ifile1)
      return
      end
