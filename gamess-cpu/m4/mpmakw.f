      subroutine mpmakw(y,z,w,w2,ibly,iblz,iblw,iblw2,a1,a2,nocca,
     1   nvirta,ncoorb,ifils,e,ifint1,ifint2,ifile1)
c
      implicit real*8  (a-h,o-z)
      dimension y(ncoorb,ncoorb),w(ncoorb,ncoorb),
     1  z(nocca*nvirta),a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),
     2 e(ncoorb),w2(ncoorb,ncoorb)
c
      nsq = ncoorb*ncoorb
      call vclr(w2,1,nsq)
c
c
      call search(1,ifint1)
      call search(1,ifint2)
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
      do 170 ip = 1 , ncoorb
         do 160 iq = 1 , ip
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            yf = (y(ip,iq)+y(iq,ip))
            if (ip.eq.iq) yf = 0.5d0*yf
            yf1 = 0.5d0*yf
            yf2 = -yf
            if (yf.ne.0.0d0) then
               do 130 k = 1 , nocca
                  do 120 j = 1 , nocca
                     w(j,k) = w(j,k)
     +                        + yf1*(a1(j,k)*4.0d0-a2(j,k)-a2(k,j))
 120              continue
 130           continue
               do 150 k = 1 , ncoorb
                  do 140 j = 1 , ncoorb
                     w2(j,k) = w2(j,k)
     +                         + yf2*(a1(j,k)*4.0d0-a2(j,k)-a2(k,j))
 140              continue
 150           continue
            end if
 160     continue
 170  continue
c
      call wrt3(w,nsq,iblw,ifile1)
      call wrt3(w2,nsq,iblw2,ifile1)
      return
      end
