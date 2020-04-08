c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp2c.m,v $
c  $State: Exp $
c  
c===========================================================
      subroutine form1(ifort1,ifort2,q,nov,n2,ifort3)
      implicit real*8  (a-h,o-z)
      dimension q(*)
      m1 = 1
      call search(m1,ifort1)
      call search(m1,ifort2)
      call search(m1,ifort3)
      do 20 i = 1 , nov
         call rdedz(q,n2,ifort1)
         call wtedz(q,n2,ifort3)
         call rdedz(q,n2,ifort2)
         call wtedz(q,n2,ifort3)
 20   continue
      return
      end
      subroutine mpdag3(ifort1,ifort2,nov,n2,b,c,ifort3,cx,bx)
      implicit real*8  (a-h,o-z)
c
c     matrix transposition ( middle part of transformation )
c     in-core version
      dimension b(nov*n2),bx(nov*n2)
      dimension c(nov),cx(nov)
c
      m1 = 1
      call search(m1,ifort1)
      call search(m1,ifort3)
c
c     read the matrix into core
c
      ib = 1
      do 20 ixx = 1 , nov
         call rdedz(b(ib),n2,ifort1)
         call rdedz(bx(ib),n2,ifort3)
         ib = ib + n2
 20   continue
c
c
c     rewind file
c     ( can be overwriting original )
c
      m1 = 1
      call search(m1,ifort2)
c
c     write out transposed matrix
c     (could cause a few paging problems here ! on a virtual
c      machine )
c
      do 40 iel = 1 , n2
         ipos = iel
         do 30 iij = 1 , nov
            c(iij) = b(ipos)
            cx(iij) = bx(ipos)
            ipos = ipos + n2
 30      continue
         call wtedz(c,nov,ifort2)
         call wtedz(cx,nov,ifort2)
 40   continue
      return
      end
      subroutine umpprj(
     + epsa,epsb,s,euhf,na,nb,nall,spin,c,h,cb,hb,cbb,
     +hbb,cbbb,hbbb,ct,ht,ctt,htt,ctb,scr,scr1,c1,h1,c2,h2,c3,h3,
     +c4,h4,c5,h5,c6,h6,c7,h7,c8,h8,ta,tb,dc3,dh3,nva,nvb,
     +   ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,ifort7,ifort8,
     +   ifort9,qiof,iqiof,maxqq,temp,eump2)
      implicit real*8  (a-h,o-z)
      integer aa,ab,ba,bb
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension qiof(maxqq),iqiof(2)
      dimension h00(3),h01(3),s01(3),s00(4),f(3,3),ff(3),z(210)
      dimension c(nvb*nvb),h(nvb*nvb),cb(nvb*nvb),hb(nvb*nvb)
      dimension cbb(nvb*nvb),hbb(nvb*nvb),cbbb(nvb*nvb),hbbb(nvb*nvb)
      dimension ct(nvb*nvb),ht(nvb*nvb),ctt(nvb*nvb),htt(nvb*nvb)
      dimension ctb(nvb*nvb),temp(2)
      dimension scr(nall*nall),scr1(nall*nall)
      dimension c1(nall*nvb),h1(nall*nvb),c2(nall*nvb),h2(nall*nvb)
      dimension c3(na*nvb),h3(na*nvb),c4(na*nvb),h4(na*nvb)
      dimension c5(na*nvb),h5(na*nvb),c6(na*nvb),h6(na*nvb)
      dimension c7(na*nvb),h7(na*nvb),c8(na*nvb),h8(na*nvb)
      dimension s(nall*nall),ta(nall*nall),tb(nall*nall)
      dimension epsa(nall),epsb(nall),dc3(na*nvb),dh3(na*nvb)
      nvab = nva*nvb
      nvaa = nva*nva
c     nvbb = nvb*nvb
      nvbtri = nvb*(nvb+1)/2
      nvatri = nva*(nva+1)/2
      naa = na*na
      nbb = nb*nb
      nab = na*nb
      do 20 i = 1 , 210
         z(i) = 0.0d0
 20   continue
      call mxma(s,1,nall,s,nall,1,ta,1,nall,nall,nb,nall)
      call mxma(s,nall,1,s,1,nall,tb,1,nall,nall,na,nall)
      g = 0.25d0*(na+nb-2.0d0)
      do 30 i = 1 , na
         g = g - ta(i+(i-1)*nall)
 30   continue
      r2 = na*nb
      do 50 k = 1 , na
         r2 = r2 - 4.0d0*g*ta(k+(k-1)*nall)
         do 40 l = 1 , na
            r2 = r2 - 2.0d0*ta(k+(k-1)*nall)*ta(l+(l-1)*nall)
     +           - 2.0d0*ta(k+(l-1)*nall)*ta(k+(l-1)*nall)
 40      continue
 50   continue
c
c
      r1 = 0.0d0
      do 60 i = 1 , na
         r1 = r1 + ta(i+(i-1)*nall)
 60   continue
c............assume h(ka,lb,aa,bb) available on stream 1
      z1 = 0.0d0
      z2 = 0.0d0
      do 90 ki = 1 , na
         do 80 kk = 1 , na
            z1 = z1 + ta(ki+(kk-1)*nall)*ta(ki+(kk-1)*nall)
            do 70 kj = 1 , na
               z2 = z2 + ta(ki+(kk-1)*nall)*ta(kk+(kj-1)*nall)
     +              *ta(kj+(ki-1)*nall)
 70         continue
 80      continue
 90   continue
      s00(4) = 4.0d0*g*(r1*r1+4.0d0*z1-4.0d0*r1-r2) + r1*r2 -
     +         8.0d0*r1*r1 - 4.0d0*r1 + 8.0d0*z1*(r1-1) + 12.0d0*z2
      m1 = 1
      call search(m1,ifort1)
      call search(m1,ifort9)
      do 130 lb = 1 , nb
         do 120 ka = 1 , na
c.......read in h(ka,lb,aa,bb)to h(aa,bb)
            call rdedz(h,nvab,ifort1)
            do 110 aa = 1 , nva
               do 100 bb = 1 , nvb
                  c(aa+(bb-1)*nva) = -h(aa+(bb-1)*nva)
     +                               /(epsa(aa+na)+epsb(bb+nb)-epsa(ka)
     +                               -epsb(lb))
 100           continue
 110        continue
c..........write out c(ka,lb,aa,bb) to stream 1
            call wtedz(c,nvab,ifort9)
 120     continue
 130  continue
c.....sort c(ka,lb,aa,bb) to c(aa,bb,ka,lb)
c.....sort h(ka,lb,aa,bb)to h(aa,bb,ka,lb)
      call mpsrt3(ifort1,ifort6,nab,nvab,ifort9,qiof,iqiof,maxqq)
      call form1(ifort1,ifort9,qiof,nab,nvab,ifort7)
      call delfil(ifort1)
      call delfil(ifort9)
      call vclr(qiof,1,maxqq)
      call search(m1,ifort6)
      call search(m1,ifort2)
      call search(m1,ifort3)
      call search(m1,ifort4)
      do 150 bb = 1 , nvb
         do 140 aa = 1 , nva
c..... read in c to c(ka,lb)
c..... read in h to h(ka,lb)
            call rdedz(h,nab,ifort6)
            call rdedz(c,nab,ifort6)
c......form h(aa,bb,ka,ib)*t(ib,lb)=hb(aa,bb,ka,lb)
            call mxma(c(1),1,na,tb(1),1,nall,cb(1),1,na,na,nb,nb)
            call mxma(h(1),1,na,tb(1),1,nall,hb(1),1,na,na,nb,nb)
c.......form h(aa,bb,ia,lb)*t(ia,ka)=hbbb(aa,bb,ka,lb)
            call mxma(ta(1),nall,1,c(1),1,na,cbbb(1),1,na,na,na,nb)
            call mxma(ta(1),nall,1,h(1),1,na,hbbb(1),1,na,na,na,nb)
c.......form h(aa,bb,ia,ib)*s(ia,lb)*s(ka,ib)=hbb(aa,bb,ka,lb)
            call mxma(s(1),1,nall,c(1),na,1,scr(1),1,na,na,nb,na)
            call mxma(scr(1),1,na,s(1),1,nall,cbb(1),1,na,na,na,nb)
            call mxma(s(1),1,nall,h(1),na,1,scr(1),1,na,na,nb,na)
            call mxma(scr(1),1,na,s(1),1,nall,hbb(1),1,na,na,na,nb)
            call wtedz(cb,nab,ifort2)
            call wtedz(hb,nab,ifort2)
            call wtedz(cbb,nab,ifort3)
            call wtedz(hbb,nab,ifort3)
            call wtedz(cbbb,nab,ifort4)
            call wtedz(hbbb,nab,ifort4)
c....
c.....write out cb,hb on styream 2
c.....write out cbb,hbb on stream 3
c.....write out cbbb,hbbb on stream 4
 140     continue
 150  continue
      call delfil(ifort6)
c....sort cb(aa,bb,ka,lb) to cb(ka,lb,aa,bb)
      call mpsrt3(ifort2,ifort2,nvab,nab,ifort2,qiof,iqiof,maxqq)
      call mpsrt3(ifort3,ifort3,nvab,nab,ifort3,qiof,iqiof,maxqq)
      call mpsrt3(ifort4,ifort4,nvab,nab,ifort4,qiof,iqiof,maxqq)
c.....do same for hb,cbb,hbb,cbbb,hbbb
      call vclr(qiof,1,maxqq)
      call search(m1,ifort7)
      call search(m1,ifort2)
      call search(m1,ifort3)
      call search(m1,ifort4)
      do 340 lb = 1 , nb
         do 330 ka = 1 , na
            call rdedz(h,nvab,ifort7)
            call rdedz(c,nvab,ifort7)
            call rdedz(cb,nvab,ifort2)
            call rdedz(hb,nvab,ifort2)
            call rdedz(cbb,nvab,ifort3)
            call rdedz(hbb,nvab,ifort3)
            call rdedz(cbbb,nvab,ifort4)
            call rdedz(hbbb,nvab,ifort4)
c.......read in h(ka,lb,aa,bb) into h(aa,bb)
c.......do same for c,cb,hb,cbb,hbb,cbbb,hbbb
c.....on streams 1,2,3,4
c....form c1(pb,bb)=c(aa,bb)*s(aa,pb) and h1
            call mxma(s(na+1),nall,1,c(1),1,nva,c1(1),1,nall,nall,nva,
     +                nvb)
            call mxma(s(na+1),nall,1,h(1),1,nva,h1(1),1,nall,nall,nva,
     +                nvb)
c.....form c2(pa,aa)=c(aa,bb)*s(pa,bb) and h2
            call mxma(s(nall*nb+1),1,nall,c(1),nva,1,c2(1),1,nall,nall,
     +                nvb,nva)
            call mxma(s(nall*nb+1),1,nall,h(1),nva,1,h2(1),1,nall,nall,
     +                nvb,nva)
c.......form c3(lb,aa)=c(ka,lb,aa,bb)*s(ka,bb) and h3
            call mxmb(c(1),1,nva,s(ka+nb*nall),nall,1,c3(lb),nb,1,nva,
     +                nvb,1)
            call mxmb(h(1),1,nva,s(ka+nb*nall),nall,1,h3(lb),nb,1,nva,
     +                nvb,1)
c......form c4(ka,db)=c(ka,lb,aa,db)*s(aa,lb) and h4
            call mxmb(c(1),nva,1,s(na+1+(lb-1)*nall),1,1,c4(ka),na,1,
     +                nvb,nva,1)
            call mxmb(h(1),nva,1,s(na+1+(lb-1)*nall),1,1,h4(ka),na,1,
     +                nvb,nva,1)
c.......form c5(ka,ca)=c(ia,lb,ca,bb)*t(lb,bb) and h5
            call mxmb(c(1),1,nva,tb(lb+nb*nall),nall,1,c5(ka),na,1,nva,
     +                nvb,1)
            call mxmb(h(1),1,nva,tb(lb+nb*nall),nall,1,h5(ka),na,1,nva,
     +                nvb,1)
c......form c6(lb,db)=c(ka,lb,aa,db)*t(aa,ka) and h6
            call mxmb(c(1),nva,1,ta(na+1+(ka-1)*nall),1,1,c6(lb),nb,1,
     +                nvb,nva,1)
            call mxmb(h(1),nva,1,ta(na+1+(ka-1)*nall),1,1,h6(lb),nb,1,
     +                nvb,nva,1)
c.......form c7(ka,db)=c(ka,lb,aa,db)*s(aa,kb)*t(lb,kb) and h7
            call mxma(s(na+1),1,nall,tb(1+(lb-1)*nall),1,nall,scr(1),1,
     +                nva,nva,nb,1)
            call mxmb(c(1),nva,1,scr(1),1,1,c7(ka),na,1,nvb,nva,1)
            call mxmb(h(1),nva,1,scr(1),1,1,h7(ka),na,1,nvb,nva,1)
c.......form c8(lb,aa)=c(ka,lb,aa,bb)*s(la,bb)*t(ka,la) and h8
            call mxma(s(nb*nall+1),nall,1,ta(1+(ka-1)*nall),1,nall,
     +                scr(1),1,nvb,nvb,na,1)
            call mxmb(c(1),1,nva,scr(1),1,1,c8(lb),nb,1,nva,nvb,1)
            call mxmb(h(1),1,nva,scr(1),1,1,h8(lb),nb,1,nva,nvb,1)
c.....no 2 c2(la,aa)*hbbb(aa,bb)*s(la,bb)=c2(la,aa)*scr(la,aa)
            call mxma(hbbb(1),1,nva,s(nb*nall+1),nall,1,scr(1),na,1,nva,
     +                nvb,na)
            do 160 la = 1 , na
               z(2) = z(2) - ddot(nva,c2(la),nall,scr(la),na)
 160        continue
c.....no 4 c1(kb,db)*hbbb(ca,db)*s(ca,kb)=c1(kb,db)*scr(kb,db)
            call mxma(hbbb(1),nva,1,s(na+1),1,nall,scr(1),nall,1,nvb,
     +                nva,nall)
            do 170 kb = 1 , nb
               z(4) = z(4) - ddot(nvb,c1(kb),nall,scr(kb),nall)
 170        continue
c.....no 6 c1(bb,db)*hbbb(ca,bb)*s(ca,db)=c1(bb,db)*scr(db,bb)
            do 180 bb = 1 , nvb
               z(6) = z(6)
     +             - ddot(nvb,c1(nb+bb),nall,scr((bb-1)*nall+nb+1),1)
 180        continue
c.....no 3 c2(la,aa)*hb(aa,bb)*s(la,bb)=c2(la,aa)*scr(la,aa)
            call mxma(hb(1),1,nva,s(nb*nall+1),nall,1,scr(1),na,1,nva,
     +                nvb,na)
            do 190 la = 1 , na
               z(3) = z(3) - ddot(nva,c2(la),nall,scr(la),na)
 190        continue
c.....no 5 c1(kb,db)*hb(ca,db)*s(ca,kb)=c1(kb,db)*scr(kb,db)
            call mxma(hb(1),nva,1,s(na+1),1,nall,scr(1),nall,1,nvb,nva,
     +                nall)
            do 200 kb = 1 , nb
               z(5) = z(5) - ddot(nvb,c1(kb),nall,scr(kb),nall)
 200        continue
c.....no 7 c1(bb,db)*hb(ca,bb)*s(ca,db)=c1(bb,db)*scr(db,bb)
            do 210 bb = 1 , nvb
               z(7) = z(7)
     +             - ddot(nvb,c1(nb+bb),nall,scr((bb-1)*nall+nb+1),1)
 210        continue
c.....no 23 c2(la,aa)*hbb(aa,bb)*s(la,bb)=c2(la,aa)*scr(la,aa)
            call mxma(hbb(1),1,nva,s(nb*nall+1),nall,1,scr(1),na,1,nva,
     +                nvb,na)
            do 220 la = 1 , na
               z(23) = z(23) + ddot(nva,c2(la),nall,scr(la),na)
 220        continue
c.....no 29 c1(kb,db)*hbb(ca,db)*s(ca,kb)=c1(kb,db)*scr(kb,db)
            call mxma(hbb(1),nva,1,s(na+1),1,nall,scr(1),nall,1,nvb,nva,
     +                nall)
            do 230 kb = 1 , nb
               z(29) = z(29) + ddot(nvb,c1(kb),nall,scr(kb),nall)
 230        continue
c.....no 24 c1(bb,db)*hb(ca,bb)*s(ca,db)=c1(bb,db)*scr(db,bb)
            do 240 bb = 1 , nvb
               z(24) = z(24)
     +              + ddot(nvb,c1(nb+bb),nall,scr((bb-1)*nall+nb+1),1)
 240        continue
c......no 1 -c2(la,ca)*h2(ia,ca)*t(ia,la)=-scr(ia,la)*ta(ia,la)
            call mxma(h2(1),1,nall,c2(1),nall,1,scr(1),1,na,na,nva,na)
            do 250 ia = 1 , na
               z(1) = z(1) - ddot(na,scr(ia),na,ta(ia),nall)
               z(11) = z(11) - g*scr(ia+(ia-1)*na)
 250        continue
c.....no 11 -g*c2(la,ca)*h2(la,ca)
c.....no 8 -c1(kb,db)*h1(ib,db)*t(ib,kb)=-scr(ib,kb)*t(ib,kb)
            call mxma(h1(1),1,nall,c1(1),nall,1,scr(1),1,nb,nb,nvb,nb)
            do 260 kb = 1 , nb
               z(8) = z(8)
     +            - ddot(nb,scr((kb-1)*nb+1),1,tb((kb-1)*nall+1),1)
               z(13) = z(13) - g*scr(kb+(kb-1)*nb)
 260        continue
c.....no 13 -g*c1(kb,db)*h1(kb,db)
c...no 12 -g*c1(bb,db)*h1(db,bb)
            do 270 bb = 1 , nvb
               z(12) = z(12)
     +            - g*ddot(nvb,c1(nb+bb),nall,h1(nb+1+(bb-1)*nall),1)
 270        continue
c....no 21 c2(la,ca)*s(ca,kb)*h1(kb,bb)*s(la,bb)=scr(la,kb)*scr1(la,kb)
            call mxma(c2(1),1,nall,s(na+1),1,nall,scr(1),1,na,na,nva,nb)
            call mxma(s(nb*nall+1),1,nall,h1(1),nall,1,scr1(1),1,na,na,
     +                nvb,nb)
            z(21) = z(21) + ddot(na*nb,scr(1),1,scr1(1),1)
c....no 39 -c1(kb,db)*h1(db,bb)*t(kb,bb)
            call mxma(c1(1),1,nall,h1(nb+1),1,nall,scr(1),1,nb,nb,nvb,
     +                nvb)
            do 280 kb = 1 , nb
               z(39) = z(39) - ddot(nvb,scr(kb),nb,tb(kb+nb*nall),nall)
 280        continue
c...no 40 -c1(bb,db)*h1(kb,bb)*t(kb,db)
            call mxma(h1(1),1,nall,c1(nb+1),1,nall,scr(1),1,nb,nb,nvb,
     +                nvb)
            do 290 kb = 1 , nb
               z(40) = z(40) - ddot(nvb,scr(kb),nb,tb(kb+nb*nall),nall)
 290        continue
c...nos 41,42,43,44,45,46,47,48,49
            z(41) = z(41) + ddot(nvab,cb(1),1,hbbb(1),1)
            z(42) = z(42) - g*ddot(nvab,c(1),1,hbb(1),1)
            z(43) = z(43) - ddot(nvab,cbbb(1),1,hbb(1),1)
            z(44) = z(44) - ddot(nvab,hbbb(1),1,cbb(1),1)
            z(45) = z(45) + ddot(nvab,cbbb(1),1,hbbb(1),1)
            z(46) = z(46) + g*ddot(nvab,c(1),1,hbbb(1),1)
            z(47) = z(47) + ddot(nvab,cb(1),1,hb(1),1)
            z(48) = z(48) + g*ddot(nvab,c(1),1,hb(1),1)
            z(49) = z(49) + 0.25d0*r2*ddot(nvab,c(1),1,h(1),1)
c....nos 50,51,52,53,163
            z(163) = z(163) + ddot(nvab,c(1),1,h(1),1)
            z(50) = z(50) + r1*ddot(nvab,c(1),1,h(1),1)
            z(51) = z(51) + ddot(nvab,c(1),1,hbb(1),1)
            z(52) = z(52) - ddot(nvab,c(1),1,hb(1),1)
            z(202) = z(202) - ddot(nvab,h,1,cb,1)
            z(201) = z(201) - ddot(nvab,h,1,cbb,1)
            z(53) = z(53) - ddot(nvab,c(1),1,hbbb(1),1)
c....no 54 c1(db,bb)*h1(bb,db)
            do 300 bb = 1 , nvb
               z(54) = z(54) + ddot(nvb,c1((bb-1)*nall+nb+1),1,h1(bb+nb)
     +                 ,nall)
 300        continue
c....no 55 c2(ka,aa)*h2(ka,aa)
            do 310 ia = 1 , na
               z(55) = z(55) + ddot(nva,c2(ia),nall,h2(ia),nall)
 310        continue
c.....no 56 c1(mb,bb)*h1(mb,bb)
            do 320 jb = 1 , nb
               z(56) = z(56) + ddot(nvb,c1(jb),nall,h1(jb),nall)
 320        continue
 330     continue
 340  continue
      call delfil(ifort3)
c.....end of big read loop
c.....no 9 c4(ka,db)*h4(ia,db)*t(ia,ka)
      call mxma(c4(1),1,na,h4(1),na,1,scr(1),1,na,na,nvb,na)
      do 350 ka = 1 , na
         z(9) = z(9) + ddot(na,scr(ka),na,ta(1+(ka-1)*nall),1)
         z(14) = z(14) + g*scr(ka+(ka-1)*na)
 350  continue
c......no 14 g*c4(ka,db)*h4(ka,db)
c......no 10 c3(lb,ca)*h3(ib,ca)*t(lb,ib)
      call mxma(c3(1),1,nb,h3(1),nb,1,scr(1),1,nb,nb,nva,nb)
      do 360 lb = 1 , nb
         z(10) = z(10) + ddot(nb,scr(lb),nb,tb(1+(lb-1)*nall),1)
         z(15) = z(15) + g*scr(lb+(lb-1)*nb)
 360  continue
c........no 15 g*c3(lb,ca)*h3(lb,ca)
c........no 16 -c3(lb,ca)*h5(ia,ca)*s(ia,lb)
      call mxma(c3(1),1,nb,h5(1),na,1,scr(1),1,nb,nb,nva,na)
      do 370 ia = 1 , na
         z(16) = z(16) - ddot(nb,scr(1+(ia-1)*nb),1,s(ia),nall)
 370  continue
c........no 17 -h3(lb,ca)*c5(ia,ca)*s(ia,lb)
      call mxma(h3(1),1,nb,c5(1),na,1,scr(1),1,nb,nb,nva,na)
      do 380 ia = 1 , na
         z(17) = z(17) - ddot(nb,scr(1+(ia-1)*nb),1,s(ia),nall)
 380  continue
c......nos 18,19,22,28,35,36
      z(18) = ddot(nb*nva,c3(1),1,h8(1),1)
      z(19) = ddot(nb*nva,h3(1),1,c8(1),1)
      z(22) = ddot(na*nva,c5(1),1,h5(1),1)
      z(28) = ddot(nb*nvb,c6(1),1,h6(1),1)
      z(35) = ddot(na*nvb,c4(1),1,h7(1),1)
      z(36) = ddot(na*nvb,h4(1),1,c7(1),1)
c.....no 27 -c3(lb,aa)*h3(lb,ca)*t(aa,ca)
      call mxma(h3(1),1,nb,ta(1+na+na*nall),1,nall,scr(1),1,nb,nb,nva,
     +          nva)
      z(27) = -ddot(nb*nva,c3(1),1,scr(1),1)
c.....no 20 -c4(ka,db)*h4(ka,bb)*t(bb,db)
      call mxma(h4(1),1,na,tb(1+nb+nb*nall),1,nall,scr(1),1,na,na,nvb,
     +          nvb)
      z(20) = -ddot(na*nvb,c4(1),1,scr(1),1)
c......no 25 c3(lb,aa)*h6(lb,bb)*s(aa,bb)
      call mxma(h6(1),1,nb,s(na+1+nb*nall),nall,1,scr(1),1,nb,nb,nvb,
     +          nva)
      z(25) = ddot(nb*nva,c3(1),1,scr(1),1)
c......no 26 h3(lb,aa)*c6(lb,bb)*s(aa,bb)
      call mxma(c6(1),1,nb,s(na+1+nb*nall),nall,1,scr(1),1,nb,nb,nvb,
     +          nva)
      z(26) = ddot(nb*nva,h3(1),1,scr(1),1)
c.........no 37 c4(ka,db)*h5(ka,aa)*s(aa,db)
      call mxma(h5(1),1,na,s(1+na+nb*nall),1,nall,scr(1),1,na,na,nva,
     +          nvb)
      z(37) = ddot(na*nvb,c4(1),1,scr(1),1)
c.........no 38 h4(ka,db)*c5(ka,aa)*s(aa,db)
      call mxma(c5(1),1,na,s(1+na+nb*nall),1,nall,scr(1),1,na,na,nva,
     +          nvb)
      z(38) = ddot(na*nvb,h4(1),1,scr(1),1)
c....no 33 c6(lb,db)*h4(ia,db)*s(ia,lb)
      call mxma(c6(1),1,nb,h4(1),na,1,scr(1),1,nb,nb,nvb,na)
      do 390 ia = 1 , na
         z(33) = z(33) - ddot(nb,scr((ia-1)*nb+1),1,s(ia),nall)
 390  continue
c....no 34 h6(lb,db)*c4(ia,db)*s(ia,lb)
      call mxma(h6(1),1,nb,c4(1),na,1,scr(1),1,nb,nb,nvb,na)
      do 400 ia = 1 , na
         z(34) = z(34) - ddot(nb,scr((ia-1)*nb+1),1,s(ia),nall)
 400  continue
c......no 30 c3(lb,aa)*s(aa,bb)*h4(ia,bb)*s(ia,lb)
      call mxma(c3(1),1,nb,s(1+na+nb*nall),1,nall,scr(1),1,nb,nb,nva,
     +          nvb)
      call mxma(s(1),nall,1,h4(1),1,na,scr1(1),1,nb,nb,na,nvb)
      z(30) = -ddot(nb*nvb,scr(1),1,scr1(1),1)
c......no 31 h3(lb,aa)*s(aa,bb)*c4(ia,bb)*s(ia,lb)
      call mxma(h3(1),1,nb,s(1+na+nb*nall),1,nall,scr(1),1,nb,nb,nva,
     +          nvb)
      call mxma(s(1),nall,1,c4(1),1,na,scr1(1),1,nb,nb,na,nvb)
      z(31) = -ddot(nb*nvb,scr(1),1,scr1(1),1)
c.......no 32 (c3(lb,aa)*s(aa,lb))*(h3(kb,ca)*s(ca,kb))
      zz1 = 0.0d0
      zz2 = 0.0d0
      do 410 lb = 1 , nb
         zz1 = zz1 + ddot(nva,c3(lb),nb,s((lb-1)*nall+na+1),1)
         zz2 = zz2 + ddot(nva,h3(lb),nb,s((lb-1)*nall+na+1),1)
 410  continue
      z(32) = zz1*zz2
c.....no 57,58
      z(57) = -ddot(na*nvb,c4(1),1,h4(1),1)
      z(58) = -ddot(nb*nva,c3(1),1,h3(1),1)
c.....no 60,61
      do 420 jb = 1 , nb
         z(60) = z(60) + ddot(nva,h3(jb),nb,s((jb-1)*nall+na+1),1)
         z(61) = z(61) + ddot(nva,c3(jb),nb,s((jb-1)*nall+na+1),1)
 420  continue
c....no 70,74
      do 430 jb = 1 , nb
         z(70) = z(70) + ddot(nvb,h6(jb),nb,tb(jb+nb*nall),nall)
         z(74) = z(74) + ddot(nvb,c6(jb),nb,tb(jb+nb*nall),nall)
 430  continue
c....nos 71,75,73,77
      do 440 ia = 1 , na
         z(71) = z(71) - ddot(nvb,h7(ia),na,s(ia+nb*nall),nall)
         z(75) = z(75) - ddot(nvb,c7(ia),na,s(ia+nb*nall),nall)
         z(73) = z(73) - g*ddot(nvb,h4(ia),na,s(ia+nb*nall),nall)
         z(77) = z(77) - g*ddot(nvb,c4(ia),na,s(ia+nb*nall),nall)
 440  continue
c...no 72,76
      call mxma(s(1),1,nall,tb(1+nb*nall),1,nall,scr(1),1,na,na,nb,nvb)
      do 450 ia = 1 , na
         z(72) = z(72) - ddot(nvb,h4(ia),na,scr(ia),na)
         z(76) = z(76) - ddot(nvb,c4(ia),na,scr(ia),na)
 450  continue
c....can now lose streams with c,h,cbb,hbb on streams 1 and 3
c.....end of (ab|ab)
      call search(m1,ifort5)
      call search(m1,ifort9)
c.....assume h(ka,la,aa,ba) available for (aa|aa) case ,stream 5
      do 490 la = 1 , na
         do 480 ka = 1 , na
            call rdedz(h,nvatri,ifort5)
c...read in h(ka,la,aa,ba) into h(aa,ba)
            iab = 0
            do 470 aa = 1 , nva
               do 460 ba = 1 , aa
                  iab = iab + 1
                  c(iab) = -h(iab)
     +                     /(epsa(aa+na)+epsa(ba+na)-epsa(ka)-epsa(la))
 460           continue
 470        continue
c....write out c(aa,ba) on stream 5
            call wtedz(c,nvatri,ifort9)
 480     continue
 490  continue
c....sort c(ka,la,aa,ba) into c(aa,ba,ka,la), and also for h
      call mpsrt3(ifort5,ifort6,naa,nvatri,ifort9,qiof,iqiof,maxqq)
      call form1(ifort5,ifort9,qiof,naa,nvaa,ifort1)
      call delfil(ifort5)
      call delfil(ifort9)
      call vclr(qiof,1,maxqq)
      call search(m1,ifort6)
      call search(m1,ifort5)
      call search(m1,ifort9)
      call search(m1,ifort3)
      do 500 iab = 1 , nvatri
c....read c(aa,ba,ka,la) into c(ka,la)
c....read h(aa,ba,ka,la) into h(ka,la)... stream 5...
         call rdedz(h,naa,ifort6)
         call rdedz(c,naa,ifort6)
c....c(ia,ja)*s(ja,lb)=ct(ia,lb)
         call mxma(c(1),1,na,s(1),1,nall,ct(1),1,na,na,na,nb)
         call mxma(h(1),1,na,s(1),1,nall,ht(1),1,na,na,na,nb)
c......c(ia,ja)*t(ja,ma)*t(la,ma)=ctb(ia,la)
         call mxma(c(1),1,na,ta(1),1,nall,scr(1),1,na,na,na,na)
         call mxma(scr(1),1,na,ta(1),1,nall,ctb(1),1,na,na,na,na)
c.....c(ia,ja)*s(ia,mb)*s(ja,nb)=ctt(mb,nb) and htt
         call mxma(c(1),1,na,s(1),1,nall,scr(1),1,na,na,na,nb)
         call mxma(s(1),nall,1,scr(1),1,na,ctt(1),1,nb,nb,na,nb)
         call mxma(h(1),1,na,s(1),1,nall,scr(1),1,na,na,na,nb)
         call mxma(s(1),nall,1,scr(1),1,na,htt(1),1,nb,nb,na,nb)
         call wtedz(ct,nab,ifort5)
         call wtedz(ht,nab,ifort5)
         call wtedz(ctt,nbb,ifort9)
         call wtedz(htt,nbb,ifort9)
         call wtedz(ctb,naa,ifort3)
c....write out ct(ia,lb),ht(ia,lb),ctt(mb.lb),htt(mb,lb),ctb(ia,ma)
c.... on streams 6,7,3, respectively
 500  continue
      call delfil(ifort6)
      call mpsrt3(ifort5,ifort5,nvatri,nab,ifort5,qiof,iqiof,maxqq)
      call mpsrt3(ifort9,ifort9,nvatri,nbb,ifort9,qiof,iqiof,maxqq)
      call mpsrt0(ifort3,ifort3,nvatri,naa,qiof,iqiof,maxqq)
      call vclr(qiof,1,maxqq)
c.... sort ct(aa,ba,ia,lb) to ct(ia,lb,aa,ba)
c.... sort ht(aa,ba,ia,lb) to ht(ia,lb,aa,ba)
c.... sort ctt(aa,ba,mb,lb) to ctt(mb,lb,aa,ba)
c.... sort htt(aa,ba,mb,lb) to htt(mb,lb,aa,ba)
c.... sort ctb(aa,ba,ia,ma) to ctb(ia,ma,aa,ba)
      call search(m1,ifort3)
      call search(m1,ifort1)
      call search(m1,ifort9)
      call search(m1,ifort5)
      do 530 ja = 1 , na
         do 520 ia = 1 , na
c....... read in c(ia,ja,aa,ba) to c(nva,nva), and h on stream 5
c....... read in ctb(ia,ja,aa,ba) to ctb(nva,nva) on stream 3
            call rdedz(temp,nvatri,ifort1)
            call squars(temp,h,nva)
            call rdedz(temp,nvatri,ifort1)
            call squars(temp,c,nva)
            call rdedz(temp,nvatri,ifort3)
            call squars(temp,ctb,nva)
c....form c1(kb,aa)=c(ia,ja,aa,ba)*s(ba,kb)
c.....form dc3(ia,aa)=c(ia,ja,aa,ba)*t(ba,ja), and h1 and dh3
            call mxma(s(1+na),nall,1,c(1),nva,1,c1(1),1,nb,nb,nva,nva)
            call mxma(s(1+na),nall,1,h(1),nva,1,h1(1),1,nb,nb,nva,nva)
            call mxmb(c(1),1,nva,ta(1+na+(ja-1)*nall),1,1,dc3(ia),na,1,
     +                nva,nva,1)
            call mxmb(h(1),1,nva,ta(1+na+(ja-1)*nall),1,1,dh3(ia),na,1,
     +                nva,nva,1)
c......no 88,90,92,84,86,161
            z(161) = z(161) + 0.25d0*ddot(nva*nva,c(1),1,h(1),1)
            z(88) = z(88) + 0.25d0*r2*ddot(nva*nva,c(1),1,h(1),1)
            z(90) = z(90) + 0.25d0*r1*ddot(nva*nva,c(1),1,h(1),1)
            z(92) = z(92) + 0.5d0*ddot(nb*nva,c1(1),1,h1(1),1)
            z(84) = z(84) - 2.0d0*g*ddot(nb*nva,c1(1),1,h1(1),1)
            z(86) = z(86) + 2.0d0*ddot(nva*nva,ctb(1),1,h(1),1)
c......no 81
            call mxma(c1(1),1,nb,h1(1),nb,1,scr(1),1,nb,nb,nva,nb)
            do 510 kb = 1 , nb
               z(81) = z(81) - 2.0d0*ddot(nb,scr(kb),nb,tb(kb),nall)
 510        continue
c......no 83
            call mxma(h1(1),1,nb,ta(na+1+na*nall),1,nall,scr(1),1,nb,nb,
     +                nva,nva)
            z(83) = z(83) + ddot(nb*nva,c1(1),1,scr(1),1)
 520     continue
 530  continue
      call delfil(ifort1)
c....no 82
      z(82) = 4.0d0*ddot(na*nva,dc3(1),1,dh3(1),1)
c....nos 89,79
      do 540 ia = 1 , na
         z(89) = z(89) + 2.0d0*ddot(nva,dh3(ia),na,ta(ia+na*nall),nall)
         z(79) = z(79) + 2.0d0*ddot(nva,dc3(ia),na,ta(ia+na*nall),nall)
 540  continue
      do 560 jb = 1 , nb
         do 550 ib = 1 , nb
c.... read in ctt(nvab) and htt on strwam 7
            call rdedz(ctt,nvatri,ifort9)
            call rdedz(htt,nvatri,ifort9)
            z(85) = z(85) + 2.0d0*ddot(nvatri,ctt(1),1,htt(1),1)
 550     continue
 560  continue
c... now lose ctb on stream 3
c       call delfil(ifort3)
      do 580 lb = 1 , nb
         do 570 ka = 1 , na
c.... read in ct(ka,lb,aa,ba) to ct(aa,ba) on stream 6
            call rdedz(temp,nvatri,ifort5)
            call squars(temp,ct,nva)
c.... read in ht(ka,lb,aa,ba) to ht(aa,ba) on stream 6
            call rdedz(temp,nvatri,ifort5)
            call squars(temp,ht,nva)
c....nos 87,91
            z(87) = z(87) + 2.0d0*g*ddot(nva*nva,ct(1),1,ht(1),1)
            z(91) = z(91) - 0.5d0*ddot(nva*nva,ct(1),1,ht(1),1)
c.....no 80
            call mxma(ct(1),1,nva,ta(1+na+na*nall),1,nall,scr(1),1,nva,
     +                nva,nva,nva)
            z(80) = z(80) - 4.0d0*ddot(nva*nva,scr(1),1,ht(1),1)
 570     continue
 580  continue
c.....  now lose ctb
      call delfil(ifort3)
c.....end of (aa|aa) terms
c....now doing (ab|aa) terms
c......no 101,121
      z(101) = ddot(na*nva,c5(1),1,dh3(1),1)
      z(121) = ddot(na*nva,h5(1),1,dc3(1),1)
c...no 100,120
      call mxma(dh3(1),1,na,c3(1),nb,1,scr(1),1,na,na,nva,nb)
      do 590 ja = 1 , na
         z(100) = z(100) - ddot(nb,s(ja),nall,scr(ja),na)
 590  continue
      call mxma(dc3(1),1,na,h3(1),nb,1,scr(1),1,na,na,nva,nb)
      do 600 ja = 1 , na
         z(120) = z(120) - ddot(nb,s(ja),nall,scr(ja),na)
 600  continue
c....nos 104,124
      call mxma(dh3(1),1,na,s(1+na+nb*nall),1,nall,scr(1),1,na,na,nva,
     +          nvb)
      z(104) = ddot(na*nvb,c4(1),1,scr(1),1)
      call mxma(dc3(1),1,na,s(1+na+nb*nall),1,nall,scr(1),1,na,na,nva,
     +          nvb)
      z(124) = ddot(na*nvb,h4(1),1,scr(1),1)
      call search(m1,ifort5)
      call search(m1,ifort2)
      call search(m1,ifort4)
      call search(m1,ifort7)
      do 670 jb = 1 , nb
         do 660 ia = 1 , na
c.....read in ht(ia,jb,aa,ca) to ht(aa,ca)  and ct on stream 6
c....read in cb(ia,jb,aa,cb) to cb(aa,cb) and hb on stream 2
c....read in cbbb(ia,jb,aa,cb) to cbbb(aa,cb) and hbbb on stream 4
c....read in c(ia,jb,aa,bb) to c(aa,bb) and h on stream 1
c..........form c2(ia,jb,pa,aa)=c(aa,bb)*s(pa,bb) and h2
            call rdedz(temp,nvatri,ifort5)
            call squars(temp,ct,nva)
            call rdedz(temp,nvatri,ifort5)
            call squars(temp,ht,nva)
            call rdedz(cb,nvab,ifort2)
            call rdedz(hb,nvab,ifort2)
            call rdedz(cbbb,nvab,ifort4)
            call rdedz(hbbb,nvab,ifort4)
            call rdedz(h,nvab,ifort7)
            call rdedz(c,nvab,ifort7)
            call mxma(s(1+nb*nall),1,nall,c(1),nva,1,c2(1),1,nall,nall,
     +                nvb,nva)
            call mxma(s(1+nb*nall),1,nall,h(1),nva,1,h2(1),1,nall,nall,
     +                nvb,nva)
c.....no 102,122
            call mxma(c2(1),1,nall,ht(1),1,nva,scr(1),1,na,na,nva,nva)
            do 610 ka = 1 , na
               z(102) = z(102)-ddot(nva,scr(ka),na,ta(ka+na*nall),nall)
 610        continue
            call mxma(h2(1),1,nall,ct(1),1,nva,scr(1),1,na,na,nva,nva)
            do 620 ka = 1 , na
               z(122) = z(122)-ddot(nva,scr(ka),na,ta(ka+na*nall),nall)
 620        continue
c.......no 103,123,106,126,105,125
            call mxma(s(1+na+nb*nall),1,nall,cbbb(1),nva,1,scr(1),1,nva,
     +                nva,nvb,nva)
            z(103) = z(103) + ddot(nva*nva,scr(1),1,ht(1),1)
            call mxma(s(1+na+nb*nall),1,nall,hbbb(1),nva,1,scr(1),1,nva,
     +                nva,nvb,nva)
            z(123) = z(123) + ddot(nva*nva,scr(1),1,ct(1),1)
            call mxma(cb(1),1,nva,s(1+na+nb*nall),nall,1,scr(1),1,nva,
     +                nva,nvb,nva)
            z(106) = z(106) - ddot(nva*nva,scr(1),1,ht(1),1)
            call mxma(hb(1),1,nva,s(1+na+nb*nall),nall,1,scr(1),1,nva,
     +                nva,nvb,nva)
            z(126) = z(126) - ddot(nva*nva,scr(1),1,ct(1),1)
            call mxma(ht(1),1,nva,ta(1+na+na*nall),1,nall,scr(1),1,nva,
     +                nva,nva,nva)
            do 630 ba = 1 , nva
               z(105) = z(105) - ddot(nva,scr(ba),nva,c2(na+ba),nall)
 630        continue
            call mxma(ct(1),1,nva,ta(1+na+na*nall),1,nall,scr(1),1,nva,
     +                nva,nva,nva)
            do 640 ba = 1 , nva
               z(125) = z(125) - ddot(nva,scr(ba),nva,h2(na+ba),nall)
 640        continue
c......no 107
            do 650 ba = 1 , nva
               z(107) = z(107)
     +                  - g*ddot(nva,c2(ba+na),nall,ht(1+(ba-1)*nva),1)
               z(127) = z(127)
     +                  - g*ddot(nva,h2(ba+na),nall,ct(1+(ba-1)*nva),1)
               z(111) = z(111)
     +                  + ddot(nva,h2(ba+na),nall,ct(1+(ba-1)*nva),1)
               z(112) = z(112)
     +                  + ddot(nva,c2(ba+na),nall,ht(1+(ba-1)*nva),1)
 650        continue
 660     continue
 670  continue
c..... now lose ct,ht on stream 6
      call delfil(ifort5)
      call search(m1,ifort8)
      call search(m1,ifort1)
c.....end of (ab|aa) case
c.....assume h(kb,lb,ab,bb) available for (bb|bb) case
      do 710 lb = 1 , nb
         do 700 kb = 1 , nb
c... read in h(kb,lb,ab,bb) into h(ab,bb) on stream 8
            call rdedz(h,nvbtri,ifort8)
            iab = 0
            do 690 ab = 1 , nvb
               do 680 bb = 1 , ab
                  iab = iab + 1
                  c(iab) = -h(iab)
     +                     /(epsb(ab+nb)+epsb(bb+nb)-epsb(kb)-epsb(lb))
 680           continue
 690        continue
c.... write out c(ab,bb) on stream 8
            call wtedz(c,nvbtri,ifort1)
 700     continue
 710  continue
c.... sort c(kb,lb,ab,bb) into c(ab,bb,kb,lb), and also for h
      call mpsrt3(ifort8,ifort5,nbb,nvbtri,ifort1,qiof,iqiof,maxqq)
      call form1(ifort8,ifort1,qiof,nbb,nvbtri,ifort6)
      call vclr(qiof,1,maxqq)
      call delfil(ifort1)
      call delfil(ifort8)
      call search(m1,ifort5)
      call search(m1,ifort1)
      call search(m1,ifort8)
      call search(m1,ifort3)
      do 720 iab = 1 , nvbtri
c....read c(ab,bb,kb,lb) into c(kb,lb) on stream 8
c....read h(ab,bb,kb,lb) into h(kb,lb)...
         call rdedz(h,nbb,ifort5)
         call rdedz(c,nbb,ifort5)
c....c(ib,jb)*s(la,jb)=ct(la,ib)
         call mxma(c(1),1,nb,s(1),nall,1,ct(1),na,1,nb,nb,na)
         call mxma(h(1),1,nb,s(1),nall,1,ht(1),na,1,nb,nb,na)
c......c(ib,jb)*t(jb,mb)*t(lb,mb)=ctb(ib,lb)
         call mxma(c(1),1,nb,tb(1),1,nall,scr(1),1,nb,nb,nb,nb)
         call mxma(scr(1),1,nb,tb(1),1,nall,ctb(1),1,nb,nb,nb,nb)
c.....c(ib,jb)*s(ma,ib)*s(na,jb)=ctt(ma,na) and htt
         call mxma(c(1),1,nb,s(1),nall,1,scr(1),1,nb,nb,nb,na)
         call mxma(s(1),1,nall,scr(1),1,nb,ctt(1),1,na,na,nb,na)
         call mxma(h(1),1,nb,s(1),nall,1,scr(1),1,nb,nb,nb,na)
         call mxma(s(1),1,nall,scr(1),1,nb,htt(1),1,na,na,nb,na)
c.... write out ct(ia,lb),ht(ia,lb),ctt(mb.lb),htt(ma,la),ctb(ib,mb)
         call wtedz(ct,nab,ifort8)
         call wtedz(ht,nab,ifort8)
         call wtedz(ctt,naa,ifort1)
         call wtedz(htt,naa,ifort1)
         call wtedz(ctb,nbb,ifort3)
 720  continue
      call delfil(ifort5)
c.... sort ct(ab,bb,ia,lb) to ct(ia,lb,ab,bb) stream 6
c.... sort ht(ab,bb,ia,lb) to ht(ia,lb,ab,bb)
c.... sort ctt(ab,bb,ma,la) to ctt(ma,la,ab,bb) stream 9
c.... sort htt(ab,bb,ma,la) to htt(ma,la,ab,bb)
c.... sort ctb(ab,bb,ib,mb) to ctb(ib,mb,ab,bb) stream 3
      call mpsrt3(ifort8,ifort8,nvbtri,nab,ifort8,qiof,iqiof,maxqq)
      call mpsrt3(ifort1,ifort1,nvbtri,naa,ifort1,qiof,iqiof,maxqq)
      call mpsrt0(ifort3,ifort3,nvbtri,nbb,qiof,iqiof,maxqq)
      call vclr(qiof,1,maxqq)
      call search(m1,ifort3)
      call search(m1,ifort8)
      call search(m1,ifort1)
      call search(m1,ifort9)
      call search(m1,ifort6)
      call vclr(dc3,1,na*nvb)
      call vclr(dh3,1,na*nvb)
      do 790 jb = 1 , nb
         do 780 ib = 1 , nb
c....... read in c(ib,jb,ab,bb)  to c(nvb,nvb), and h on stream 8
c....... read in ctb(ib,jb,ab,bb) to ctb(nvb,nvb) on stream 3
c....... read in ctt and htt to ctt(nva,nva) stream 7
            call rdedz(temp,nvatri,ifort9)
            call squars(temp,ctt,nva)
            call rdedz(temp,nvatri,ifort9)
            call squars(temp,htt,nva)
            call rdedz(temp,nvbtri,ifort3)
            call squars(temp,ctb,nvb)
            call rdedz(temp,nvbtri,ifort6)
            call squars(temp,h,nvb)
            call rdedz(temp,nvbtri,ifort6)
            call squars(temp,c,nvb)
c....form c1(pa,ab)=c(ib,jb,ab,bb)*s(pa,bb)
c.....form dc3(ib,ab)=c(ib,jb,ab,bb)*t(bb,jb), and h1 and dh3
            call mxma(s(1+nb*nall),1,nall,c(1),nvb,1,c1(1),1,nall,nall,
     +                nvb,nvb)
            call mxma(s(1+nb*nall),1,nall,h(1),nvb,1,h1(1),1,nall,nall,
     +                nvb,nvb)
            call mxmb(c(1),1,nvb,tb(1+nb+(jb-1)*nall),1,1,dc3(ib),nb,1,
     +                nvb,nvb,1)
            call mxmb(h(1),1,nvb,tb(1+nb+(jb-1)*nall),1,1,dh3(ib),nb,1,
     +                nvb,nvb,1)
c......no 188,190,192,184,186,162
            z(162) = z(162) + 0.25d0*ddot(nvb*nvb,c(1),1,h(1),1)
            z(188) = z(188) + 0.25d0*r2*ddot(nvb*nvb,c(1),1,h(1),1)
            z(190) = z(190) + 0.25d0*r1*ddot(nvb*nvb,c(1),1,h(1),1)
            do 730 ja = 1 , na
               z(192)=z(192)+0.5d0*ddot(nvb,c1(ja),nall,h1(ja),nall)
               z(184)=z(184)-2.0d0*g*ddot(nvb,c1(ja),nall,h1(ja),nall)
 730        continue
            z(186) = z(186) + 2.0d0*ddot(nvb*nvb,ctb(1),1,h(1),1)
c......no 181
            call mxma(c1(1),1,nall,h1(1),nall,1,scr(1),1,na,na,nvb,na)
            do 740 ka = 1 , na
               z(181) = z(181) - 2.0d0*ddot(na,scr(ka),na,ta(ka),nall)
 740        continue
c......no 183
            call mxma(h1(1),1,nall,tb(nb+1+nb*nall),1,nall,scr(1),1,na,
     +                na,nvb,nvb)
            do 750 ja = 1 , na
               z(183) = z(183) + ddot(nvb,c1(ja),nall,scr(ja),na)
 750        continue
c............(aa|bb) terms
            call mxma(htt(1),1,nva,s(1+na+nb*nall),1,nall,scr(1),1,nva,
     +                nva,nva,nvb)
            do 760 aa = 1 , nva
               z(110) = z(110) - ddot(nvb,c1(aa+na),nall,scr(aa),nva)
 760        continue
            call mxma(ctt(1),1,nva,s(1+na+nb*nall),1,nall,scr(1),1,nva,
     +                nva,nva,nvb)
            do 770 aa = 1 , nva
               z(115) = z(115) - ddot(nvb,h1(aa+na),nall,scr(aa),nva)
 770        continue
 780     continue
 790  continue
c....no 182
      call delfil(ifort6)
      z(182) = 4.0d0*ddot(nb*nvb,dc3(1),1,dh3(1),1)
c....nos 189,179
      do 800 ib = 1 , nb
         z(189)=z(189)+2.0d0*ddot(nvb,dh3(ib),nb,tb(ib+nb*nall),nall)
         z(179)=z(179)+2.0d0*ddot(nvb,dc3(ib),nb,tb(ib+nb*nall),nall)
 800  continue
      do 820 ja = 1 , na
         do 810 ia = 1 , na
c.... read in ctt(nva,nvb) and htt on stream 9
            call rdedz(ctt,nvbtri,ifort1)
            call rdedz(htt,nvbtri,ifort1)
            z(185) = z(185) + 2.0d0*ddot(nvbtri,ctt(1),1,htt(1),1)
 810     continue
 820  continue
c... now lose ctb and ctt and htt streams 3 and 7
      call delfil(ifort3)
      call delfil(ifort9)
      call delfil(ifort1)
      do 840 lb = 1 , nb
         do 830 ka = 1 , na
c.... read in ct(ka,lb,ab,bb) to ct(ab,bb) on stream 6
c.... read in ht(ka,lb,ab,bb) to ht(ab,bb)
            call rdedz(temp,nvbtri,ifort8)
            call squars(temp,ct,nvb)
            call rdedz(temp,nvbtri,ifort8)
            call squars(temp,ht,nvb)
c....nos 187,191
            z(187) = z(187) + 2.0d0*g*ddot(nvb*nvb,ct(1),1,ht(1),1)
            z(191) = z(191) - 0.5d0*ddot(nvb*nvb,ct(1),1,ht(1),1)
c.....no 180
            call mxma(ct(1),1,nvb,tb(1+nb+nb*nall),1,nall,scr(1),1,nvb,
     +                nvb,nvb,nvb)
            z(180) = z(180) - 4.0d0*ddot(nvb*nvb,scr(1),1,ht(1),1)
 830     continue
 840  continue
c.....end of (bb|bb) terms
c....now doing (ab|bb) terms
c......no 131,141
      z(131) = ddot(nb*nvb,c6(1),1,dh3(1),1)
      z(141) = ddot(nb*nvb,h6(1),1,dc3(1),1)
c...no 130,140
      call mxma(dh3(1),1,nb,c4(1),na,1,scr(1),1,nb,nb,nvb,na)
      do 850 jb = 1 , nb
         z(130) = z(130) - ddot(na,s(1+(jb-1)*nall),1,scr(jb),nb)
 850  continue
      call mxma(dc3(1),1,nb,h4(1),na,1,scr(1),1,nb,nb,nvb,na)
      do 860 jb = 1 , nb
         z(140) = z(140) - ddot(na,s(1+(jb-1)*nall),1,scr(jb),nb)
 860  continue
c....nos 134,144
      call mxma(dh3(1),1,nb,s(1+na+nb*nall),nall,1,scr(1),1,nb,nb,nvb,
     +          nva)
      z(134) = ddot(nb*nva,c3(1),1,scr(1),1)
      call mxma(dc3(1),1,nb,s(1+na+nb*nall),nall,1,scr(1),1,nb,nb,nvb,
     +          nva)
      z(144) = ddot(nb*nva,h3(1),1,scr(1),1)
      call search(m1,ifort7)
      call search(m1,ifort2)
      call search(m1,ifort4)
      call search(m1,ifort8)
      do 930 jb = 1 , nb
         do 920 ia = 1 , na
c.....read in ht(ia,jb,ab,cb) to ht(ab,cb)  and ct on stream 6
c....read in cb(ia,jb,aa,cb) to cb(aa,cb) and hb on stream 2
c....read in cbbb(ia,jb,aa,cb) to cbbb(aa,cb) and hbbb on stream 4
c....read in c(ia,jb,aa,bb) to c(aa,bb) and h on stream 1
c..........form c1(ia,jb,pb,ab)=c(aa,ab)*s(aa,pb) and h1
            call rdedz(temp,nvbtri,ifort8)
            call squars(temp,ct,nvb)
            call rdedz(temp,nvbtri,ifort8)
            call squars(temp,ht,nvb)
            call rdedz(cb,nvab,ifort2)
            call rdedz(hb,nvab,ifort2)
            call rdedz(h,nvab,ifort7)
            call rdedz(c,nvab,ifort7)
            call rdedz(cbbb,nvab,ifort4)
            call rdedz(hbbb,nvab,ifort4)
            call mxma(s(1+na),nall,1,c(1),1,nva,c1(1),1,nall,nall,nva,
     +                nvb)
            call mxma(s(1+na),nall,1,h(1),1,nva,h1(1),1,nall,nall,nva,
     +                nvb)
c.....no 132
            call mxma(c1(1),1,nall,ht(1),1,nvb,scr(1),1,nb,nb,nvb,nvb)
            do 870 kb = 1 , nb
               z(132)=z(132)-ddot(nvb,scr(kb),nb,tb(kb+nb*nall),nall)
 870        continue
            call mxma(h1(1),1,nall,ct(1),1,nvb,scr(1),1,nb,nb,nvb,nvb)
            do 880 kb = 1 , nb
               z(142)=z(142)-ddot(nvb,scr(kb),nb,tb(kb+nb*nall),nall)
 880        continue
c......no 133,143,136,146,135,145
            call mxma(s(1+na+nb*nall),nall,1,cb(1),1,nva,scr(1),1,nvb,
     +                nvb,nva,nvb)
            z(133) = z(133) + ddot(nvb*nvb,scr(1),1,ht(1),1)
            call mxma(s(1+na+nb*nall),nall,1,hb(1),1,nva,scr(1),1,nvb,
     +                nvb,nva,nvb)
            z(143) = z(143) + ddot(nvb*nvb,scr(1),1,ct(1),1)
            call mxma(cbbb(1),nva,1,s(1+na+nb*nall),1,nall,scr(1),1,nvb,
     +                nvb,nva,nvb)
            z(136) = z(136) - ddot(nvb*nvb,scr(1),1,ht(1),1)
            call mxma(hbbb(1),nva,1,s(1+na+nb*nall),1,nall,scr(1),1,nvb,
     +                nvb,nva,nvb)
            z(146) = z(146) - ddot(nvb*nvb,scr(1),1,ct(1),1)
            call mxma(ht(1),1,nvb,tb(1+nb+nb*nall),1,nall,scr(1),1,nvb,
     +                nvb,nvb,nvb)
            do 890 bb = 1 , nvb
               z(135) = z(135) - ddot(nvb,scr(bb),nvb,c1(nb+bb),nall)
 890        continue
            call mxma(ct(1),1,nvb,tb(1+nb+nb*nall),1,nall,scr(1),1,nvb,
     +                nvb,nvb,nvb)
            do 900 bb = 1 , nvb
               z(145) = z(145) - ddot(nvb,scr(bb),nvb,h1(nb+bb),nall)
 900        continue
c......no 137
            do 910 bb = 1 , nvb
               z(137) = z(137)
     +                  - g*ddot(nvb,c1(bb+nb),nall,ht(1+(bb-1)*nvb),1)
               z(147) = z(147)
     +                  - g*ddot(nvb,h1(bb+nb),nall,ct(1+(bb-1)*nvb),1)
               z(113) = z(113)
     +                  + ddot(nvb,h1(bb+nb),nall,ct(1+(bb-1)*nvb),1)
               z(114) = z(114)
     +                  + ddot(nvb,c1(bb+nb),nall,ht(1+(bb-1)*nvb),1)
 910        continue
 920     continue
 930  continue
c..... now lose ct,ht on stream 6
c.... lose cb,hb on stream 2
c.... lose cbbb,hbbb on stream 4
      call delfil(ifort8)
      call delfil(ifort2)
      call delfil(ifort4)
      call delfil(ifort7)
      s00(1) = 1.0d0
      s00(2) = r1
      s00(3) = r2
      s01(1) = 0.0d0
      s01(2) = z(61)
      s01(3) = z(79) + z(179) + 4.0d0*(z(74)+z(75)+z(76)+z(77))
      h00(1) = euhf
      h00(2) = euhf*s00(2) + z(60)
      h00(3) = euhf*s00(3) + z(89) + z(189)
     +         + 4.0d0*(z(70)+z(71)+z(72)+z(73))
      h01(1) = z(161) + z(162) + z(163)
      h01(2) = z(90) + z(91) + z(92) + z(190) + z(191) + z(192) + z(111)
     +         + z(112) + z(113) + z(114) + euhf*s01(2)
      h01(3) = z(110) + z(115) + euhf*s01(3)
      do 940 i = 1 , 9
         h01(2) = h01(2) + z(49+i)
         h01(3) = h01(3) + z(79+i) + z(179+i)
 940  continue
      h013 = 0.0d0
      h014 = 0.0d0
      h015 = 0.0d0
      h016 = 0.0d0
      do 950 i = 1 , 8
         h013 = h013 + z(99+i)
         h014 = h014 + z(129+i)
         h015 = h015 + z(119+i)
         h016 = h016 + z(139+i)
         h01(3) = h01(3) + 4.0d0*(z(99+i)+z(129+i)+z(119+i)+z(139+i))
 950  continue
      h017 = 0.0d0
      do 960 i = 1 , 49
         h017 = h017 + z(i)
         h01(3) = h01(3) + 4.0d0*z(i)
 960  continue
      f(3,1) = (1-nb/(2*spin+2))*(1-nb/(2*(2*spin+3)))
      f(3,2) = (1-nb/(4*spin+6))/(2*spin+2) + (1-nb/(2*spin+2))
     +         /(4*spin+6)
      f(3,3) = 1/((2*spin+2)*(4*spin+6))
      f(2,1) = (1-nb/(2*spin+2))
      f(2,2) = 1/(2*spin+2)
      f(1,1) = 1.0d0
      write(iwr,6010)
      do 980 iproj = 1 , 3
         do 970 i = 1 , iproj
            ff(i) = f(iproj,i)
 970     continue
         znorm = ddot(iproj,ff(1),1,s00(1),1)
         puhf = ddot(iproj,ff(1),1,h00(1),1)/znorm
         fff = 0.0d0
         if (iproj.gt.1) fff = -puhf*ddot(iproj-1,ff(2),1,s01(2),1)
         pmp2 = (ddot(iproj,ff(1),1,h01(1),1)+fff)/znorm
         ssq = spin*(spin+1) + nb - ddot(iproj,ff(1),1,s00(2),1)/znorm
         eump2 = puhf + pmp2
         write (iwr,6020) iproj-1,puhf,pmp2,ssq,eump2
 980  continue
 6010 format (/10x,49('*')/
     +        10x,'projected uhf-mp2 calculation'/10x,49('*'))
 6020 format (
     +        10x,'projection number ',i3/
     +        10x,49('*')/
     +        10x,'projected uhf energy             ',f15.8/
     +        10x,'projected second order energy    ',f15.8/
     +        10x,'s-squared after projection       ',f15.8/
     +        10x,'mp2 total energy after projection',f15.8/
     +        10x,49('*'))
      return
      end
      subroutine mpsrt3(ifort1,ifort2,nrows,ncols,ifort3,q,iq,maxqq)
      implicit real*8  (a-h,o-z)
c
c     out of core matrix transposition = sorting routine
c     for middle of transformation. uses the usual bucket-sort
c     algorithm
c     sorts 2 sets together
c
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      dimension q(maxqq),iq(2)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
c
      n2 = ncols
      nov = nrows
c
c     can the matrix transposition be done in core ?
c
      iop = 0
      if (iop.ne.1) then
         if ((2*(nov+nov*n2)).le.maxqq) then
            i1 = nov + 1
            i2 = i1 + nov*n2
            i3 = i2 + nov
            i4 = i3 + nov*n2
            if(odebug(39))write (iwr,6020)
            call mpdag3(ifort1,ifort2,nov,n2,q(i1),q(1),ifort3,q(i2),
     +  q(i3))
            return
         end if
      end if
c
c    ibl5 = no of integrals in block of sortfile
c    only machine dependent feature should be structure of /bufb/
c
      ibl5 = nsz340
      iilen = nsz340*lenwrd()/2
      if(odebug(39)) write (iwr,6030)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
c
c       are going to sort matrix which is (n2,nov) into one
c       which is (nov,n2)
c
c       maxt is the number of columns which can be held in core
c       which is the number in each bucket
c
      maxt = maxqq/(nov*2)
      niqq = maxqq - n2 - lenint(n2+n2)
      nword = (niqq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
      nbuck = ((n2-1)/maxt) + 1
      nadd = min(maxt,n2)
      maxa = nbuck*(ibl5+ibl5i) + n2 + lenint(n2+n2)
c
c     nbuck is the number of buckets actually needed
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxqq , maxa
         call caserr('insufficient core')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q,1,maxa)
      i1 = lenrel(n2) + 1
      i2 = i1 + n2
      i4 = i2 + n2
      i3 = n2 + lenint(n2+n2) + 1
c     note that i3 and i4 are the same addresses really
      call setbfa
      call rdsrt3(q,iq(i1),iq(i2),q(i3),iq(i4),ifort1,ifort3)
c
c       read through the sort file to give final result
c
      maxa = nov*nadd*2
      call wtsrt3(q(1),maxa,ifort2)
c
      call closbf(0)
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
 6020 format(/1x,'doing in-core sort for mpsrt3')
 6030 format(/1x,'doing out-of-core sort for mpsrt3')
      end
      subroutine rdsrt3(a,ia,ib,aa,iaa,ifort,ifort2)
c
      implicit real*8  (a-h,o-z)
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/bufb/nkk,mkk,g(5119)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      dimension a(*),ia(*),ib(*),aa(*),iaa(*)
      data lastb/999999/
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c       the first n2 real words of core are used as input
c       workspace.
c
c     iilen = nsz340*lenwrd()/2
c
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i)
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      iblock = 0
      call vclr(g,1,nsz340+nsz170)
c
c
      m1 = 1
      call search(m1,ifort)
      call search(m1,ifort2)
c     ninb = no of elements ib bucket(coreload)
      ninb = nadd*nov
c
      do 70 i = 1 , nov
c
c     read column i of original matrix
c
         nrow = n2
         call rdeds(a,ia,nrow,ifort)
c
         do 30 j = 1 , nrow
c     iaddr is address of (i,j) in transposed matrix
            iaddr = i + (ia(j)-1)*nov
            ibuck = (iaddr-1)/ninb
            ia(j) = iaddr - ninb*ibuck
            ib(j) = ibuck + 1
 30      continue
c
         do 40 j = 1 , nrow
c
            ibuck = ib(j)
c
c     element goes in bucket ibuck with modified address
c
            nwb = nwbuck(ibuck) + 1
            aa(ibase(ibuck)+nwb) = a(j)
            iaa(ibasen(ibuck)+nwb) = ia(j)
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - empty
c
               call stopbk
               mkk = mark(ibuck)
               nkk = nwb
               call dcopy(ibl5,aa(ibase(ibuck)+1),1,g,1)
               call pack(g(nsz341),32,iaa(ibasen(ibuck)+1),ibl5)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 40      continue
         nrow = n2
         call rdeds(a,ia,nrow,ifort2)
c
         do 50 j = 1 , nrow
c     iaddr is address of (i,j) in transposed matrix
            iaddr = i + (ia(j)-1)*nov
            ibuck = (iaddr-1)/ninb
            ia(j) = iaddr - ninb*ibuck + ninb
            ib(j) = ibuck + 1
 50      continue
c
c
         do 60 j = 1 , nrow
c
c
            ibuck = ib(j)
c
c     element goes in bucket ibuck with modified address
c
            nwb = nwbuck(ibuck) + 1
            aa(ibase(ibuck)+nwb) = a(j)
            iaa(ibasen(ibuck)+nwb) = ia(j)
            nwbuck(ibuck) = nwb
            if (nwb.eq.ibl5) then
c
c     this block full - empty
c
               call stopbk
               mkk = mark(ibuck)
               nkk = nwb
               call dcopy(ibl5,aa(ibase(ibuck)+1),1,g,1)
               call pack(g(nsz341),32,iaa(ibasen(ibuck)+1),ibl5)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 60      continue
 70   continue
c
c
c      empty anything remaining in buckets
c
      do 80 ibuck = 1 , nbuck
         nwb = nwbuck(ibuck)
         if (nwb.ne.0) then
            call stopbk
            mkk = mark(ibuck)
            nkk = nwb
            call dcopy(ibl5,aa(ibase(ibuck)+1),1,g,1)
            call pack(g(nsz341),32,iaa(ibasen(ibuck)+1),ibl5)
            call sttout
            nwbuck(ibuck) = 0
            mark(ibuck) = iblock
            iblock = iblock + nsz
         end if
 80   continue
c
c
      call stopbk
      return
      end
      subroutine wtsrt3(q,maxq,ifort)
c
      implicit real*8  (a-h,o-z)
      dimension q(maxq)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n2,nbuck
c
      integer maxbuc, mark, ibase, ibasen
      parameter (maxbuc=1500)
      common /three/ mark(maxbuc),ibase(maxbuc),ibasen(maxbuc)
c
      common/junk/nwbuck(maxbuc)
      common/sortpk/labs(1)
      common/bufb/nkk,mkk,g(5119)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
c
      integer nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl
      integer nsz510,nsz341,nsz342,nsz680
      integer nsstat,nslen,nsmax,nsort, nstid
      logical oasync
      common/blksiz/nsz,nsz512,nsz340,nsz170,nsz85,nszij,nszkl,
     * nsz510,nsz341,nsz342,nsz680,nsstat,nslen,nsmax,nsort,
     * nstid(2),oasync
c
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
c     ibl5i = lenint(ibl5)
c
      m1 = 1
      call search(m1,ifort)
c
      min = 1
      max = nadd
      call stopbk
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         mkk = mark(i)
         call vclr(q,1,maxq)
 20      if (mkk.eq.lastb) then
c
c     coulumns min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),nov,ifort)
               call wtedz(q(j+nov*nadd),nov,ifort)
               j = j + nov
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.n2) max = n2
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
            call unpack(g(nsz341),32,labs,ibl5)
            call dsctr(nkk,g,labs,q)
            go to 20
         end if
 40   continue
      return
      end
      subroutine ver_mp2c(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp2c.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
