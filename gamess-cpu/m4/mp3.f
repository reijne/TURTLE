c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp3.m,v $
c  $State: Exp $
c  
      subroutine umpm3(y1,y2,zl,a1,a2,b,nocc,ncoorb,mn,
     1 istrma,istrmb,ibly1,ibly2,iblw,uhf)
      implicit real*8  (a-h,o-z)
      dimension y1(ncoorb,ncoorb),y2(ncoorb,ncoorb),zl(ncoorb,ncoorb),
     1 a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),b(mn)
      logical uhf
c
      ifile1 = 1
      nocc1 = nocc + 1
      nsq = ncoorb*ncoorb
      call rdedx(y1,nsq,ibly1,ifile1)
      if (.not.uhf) then
         call dcopy(nsq,y1,1,y2,1)
      else
         call rdedx(y2,nsq,ibly2,ifile1)
      end if
      call rdedx(zl,nsq,iblw,ifile1)
      call rewedz(istrma)
      call rewedz(istrmb)
c
      do 50 ip = 1 , ncoorb
         do 40 iq = 1 , ip
            call rdedz(a1,nsq,istrma)
            call rdedz(a2,nsq,istrmb)
            if (iq.le.nocc) then
               if (ip.gt.nocc) then
c
                  zl(ip,iq) = -2.0d0*zl(ip,iq)
                  do 30 j = 1 , ncoorb
                     do 20 k = 1 , ncoorb
                        zl(ip,iq) = zl(ip,iq) + y1(j,k)*a1(j,k)
     +                              + y2(j,k)*a2(j,k)*2.0d0
 20                  continue
 30               continue
               end if
            end if
 40      continue
 50   continue
c
      do 70 ia = nocc1 , ncoorb
         do 60 i = 1 , nocc
            kt = (ia-nocc-1)*nocc + i
            b(kt) = zl(ia,i) + zl(i,ia)*2.0d0
 60      continue
 70   continue
      call wrt3(zl,nsq,iblw,ifile1)
      return
      end
      subroutine cuhf3(eps,ea,eb,amaa,amab,amba,ambb,nocca,
     1 noccb,ncoorb,mna,mn,istmaa,istmab,istmba,istmbb,iblw,ifw,uhf)
      implicit real*8  (a-h,o-z)
      logical uhf
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labout(1360)
      common/blkin/g(510),nint,nxtr
      dimension amaa(ncoorb,ncoorb),amab(ncoorb,ncoorb),
     1 amba(ncoorb,ncoorb),ambb(ncoorb,ncoorb),ea(ncoorb),
     1 eb(ncoorb),eps(mn)
      logical lstop,skipp
      dimension skipp(100)
      data m0/0/
c
      call izero(1360,labout,1)
      nsq = ncoorb*ncoorb
      nint = 0
c
      iai = 0
      do 30 ia = nocca + 1 , ncoorb
         do 20 i = 1 , nocca
            iai = iai + 1
            eps(iai) = 1.0d0/(ea(ia)-ea(i))
 20      continue
 30   continue
      if (uhf) then
         iai = 0
         do 50 ia = noccb + 1 , ncoorb
            do 40 i = 1 , noccb
               iai = iai + 1
               eps(mna+iai) = 1.0d0/(eb(ia)-eb(i))
 40         continue
 50      continue
      end if
c
      call search(iblw,ifw)
      call rewedz(istmaa)
      call rewedz(istmab)
      if (uhf) then
         call rewedz(istmba)
         call rewedz(istmbb)
      end if
      do 130 ip = 1 , ncoorb
         do 120 iq = 1 , ip
            call rdedz(amaa,nsq,istmaa)
            call rdedz(amab,nsq,istmab)
            if (uhf) then
               call rdedz(amba,nsq,istmba)
               call rdedz(ambb,nsq,istmbb)
            end if
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
                  iai = (ip-nocca-1)*nocca + iq
                  ibj = 0
                  do 70 ib = nocca + 1 , ncoorb
                     do 60 j = 1 , nocca
                        ibj = ibj + 1
                        if (iai.ge.ibj) then
                           nint = nint + 1
                           g(nint) = amaa(ib,j)
                           nint2 = nint + nint
                           labout(nint2-1) = iai
                           labout(nint2) = ibj
                           if (nint.eq.nintmx) then
                           call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
                           call put(g,m511,ifw)
                           nint = 0
                           call izero(1360,labout,1)
                           end if
                        end if
 60                  continue
 70               continue
c
                  ibj = 0
                  do 90 ib = noccb + 1 , ncoorb
                     do 80 j = 1 , noccb
                        ibj = ibj + 1
                        if (.not.(.not.uhf .and. iai.lt.ibj)) then
                           nint = nint + 1
                           g(nint) = amab(ib,j)*2.0d0
                           nint2 = nint + nint
                           labout(nint2-1) = iai
                           labout(nint2) = mna + ibj
                           if (.not.uhf) labout(nint2) = ibj
                           if (nint.eq.nintmx) then
                            call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
                            call put(g,m511,ifw)
                            nint = 0
                            call izero(1360,labout,1)
                           end if
                        end if
 80                  continue
 90               continue
               end if
            end if
c
            if (uhf) then
               if (ip.gt.noccb) then
                  if (iq.le.noccb) then
                     iai = (ip-noccb-1)*noccb + iq
                     ibj = 0
                     do 110 ib = noccb + 1 , ncoorb
                        do 100 j = 1 , noccb
                           ibj = ibj + 1
                           if (iai.ge.ibj) then
                           nint = nint + 1
                           g(nint) = ambb(ib,j)
                           nint2 = nint + nint
                           labout(nint2-1) = mna + iai
                           labout(nint2) = mna + ibj
                           if (nint.eq.nintmx) then
                           call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
                           call put(g,m511,ifw)
                           nint = 0
                           call izero(1360,labout,1)
                           end if
                           end if
 100                    continue
 110                 continue
                  end if
               end if
            end if
 120     continue
 130  continue
      call pack(g(num2ep+1),lab1632,labout,numlabp)
      call put(g,m511,ifw)
      call put(g,m0,ifw)
c
      skipp(1) = .false.
      lstop = .false.
c
      call chfdrv(eps,lstop,skipp)
      return
      end
      subroutine umps3(veca,vecb,y,ysym,ifytr,ifdma,ifdmb,
     1                 nocca,noccb,ncoorb,iblkya,iblkyb,uhf)
      implicit real*8  (a-h,o-z)
      dimension y(ncoorb*ncoorb),ysym(ncoorb*ncoorb),veca(ncoorb*ncoorb)
     1 ,vecb(ncoorb*ncoorb)
      logical uhf
      nsq = ncoorb*ncoorb
      ifile1 = 1
c
      call rdedx(y,nsq,ifytr,ifile1)
      do 30 i = 1 , ncoorb
         do 20 j = 1 , ncoorb
            ij = (j-1)*ncoorb + i
            ji = (i-1)*ncoorb + j
            ysym(ij) = (y(ij)+y(ji))/2.0d0
 20      continue
 30   continue
c
      call vclr(y,1,nsq)
      call mxmb(veca,1,ncoorb,veca,ncoorb,1,y,1,ncoorb,ncoorb,nocca,
     +          ncoorb)
      call vclr(veca,1,nsq)
      call mxmb(vecb,1,ncoorb,vecb,ncoorb,1,veca,1,ncoorb,ncoorb,noccb,
     +          ncoorb)
c
      call wrt3(ysym,nsq,ifytr,ifile1)
      call wrt3(y,nsq,ifdma,ifile1)
      call wrt3(veca,nsq,ifdmb,ifile1)
c
      call rdedx(veca,nsq,iblkya,ifile1)
      if (.not.uhf) then
         call dcopy(nsq,veca,1,vecb,1)
      else
         call rdedx(vecb,nsq,iblkyb,ifile1)
      end if
      do 50 i = 1 , ncoorb
         do 40 j = 1 , ncoorb
            ij = (j-1)*ncoorb + i
            ji = (i-1)*ncoorb + j
            ysym(ij) = (veca(ij)+veca(ji))/2.0d0
            y(ij) = (vecb(ij)+vecb(ji))/2.0d0
 40      continue
 50   continue
      call wrt3(ysym,nsq,iblkya,ifile1)
      call wrt3(y,nsq,iblkyb,ifile1)
      return
      end
      subroutine umpmw3(e,y1,y2,w,a1,a2,nocc,ncoorb,
     1 istrma,istrmb,ibly1,ibly2,iblw,uhf)
      implicit real*8  (a-h,o-z)
      dimension y1(ncoorb,ncoorb),y2(ncoorb,ncoorb),
     1 a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),w(ncoorb,ncoorb),
     1 e(ncoorb)
c
      logical uhf
      nsq = ncoorb*ncoorb
      nocc1 = nocc + 1
      ifile1 = 1
      call rewedz(istrma)
      call rewedz(istrmb)
      call rdedx(y1,nsq,ibly1,ifile1)
      if (.not.uhf) then
         call dcopy(nsq,y1,1,y2,1)
      else
         call rdedx(y2,nsq,ibly2,ifile1)
      end if
      call rdedx(a1,nsq,iblw,ifile1)
c
      do 30 ir = 1 , nocc
         do 20 is = 1 , nocc
            w(ir,is) = -y1(ir,is)*(e(ir)+e(is))*0.5d0
c   minus sign as y has minus occ-occ and virt-occ elements
c     w(ir,is)=-y1(ir,is)*e(is)
            w(ir,is) = w(ir,is) + a1(ir,is)
 20      continue
 30   continue
c
      do 50 ir = nocc1 , ncoorb
         do 40 is = nocc1 , ncoorb
            w(ir,is) = -y1(ir,is)*(e(ir)+e(is))*0.5d0
            w(ir,is) = w(ir,is) + a1(ir,is)
 40      continue
 50   continue
c
      do 70 ia = nocc1 , ncoorb
         do 60 i = 1 , nocc
            w(ia,i) = -y1(ia,i)*e(i)
            w(i,ia) = a1(i,ia)*2.0d0
 60      continue
 70   continue
c
      do 110 ip = 1 , ncoorb
         do 100 iq = 1 , ip
            call rdedz(a1,nsq,istrma)
            call rdedz(a2,nsq,istrmb)
            zz = 1.0d0
            if (ip.eq.iq) zz = 0.5d0
            zz = -zz
            yfa = zz*(y1(ip,iq)+y1(iq,ip))
            yfb = zz*(y2(ip,iq)+y2(iq,ip))
            yf1 = 0.5d0*yfa
            yf2 = 0.5d0*yfb
            do 90 j = 1 , nocc
               do 80 k = 1 , nocc
                  w(j,k) = w(j,k) + yf1*a1(j,k)
                  w(j,k) = w(j,k) + yf2*a2(j,k)*2.0d0
 80            continue
 90         continue
 100     continue
 110  continue
      call wrt3(w,nsq,iblw,ifile1)
      return
      end
      subroutine umpmka(umm,amm,ammz,uds,ads,t1,t2,t3,a1,a2,e1,e2,
     +           y,zl,dum,zint,nocc,noccz,ndep,ndepz,mn,mnz,ncoorb,
     +           ifort1,ifort2,ifort3,istrma,ibly,iblw,uhf)
c
      implicit real*8  (a-h,o-z)
      dimension umm(ndep),amm(ndep),ammz(ndepz),
     +          ads(mn,mnz),uds(mn,mnz),t1(ncoorb,ncoorb),
     +          t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),a1(ncoorb,ncoorb),
     +          a2(ncoorb,ncoorb),zint(mn),e1(ncoorb),e2(ncoorb),
     +          y(ncoorb,ncoorb),zl(ncoorb,ncoorb),dum(ncoorb,ncoorb)
      logical uhf
c
c
      ifile1 = 1
      nocc1 = nocc + 1
      nvirt = ncoorb - nocc
c     nvirtz = ncoorb - noccz
      noccz1 = noccz + 1
      nsq = ncoorb*ncoorb
      ijm = (nocc+1)*nocc/2
      ijmz = (noccz+1)*noccz/2
      call vclr(y,1,nsq)
      call vclr(zl,1,nsq)
      call rewedz(ifort1)
      call rewedz(ifort2)
      if (uhf) call rewedz(ifort3)
      call rewedz(istrma)
c
      do 650 ip = 1 , ncoorb
         do 640 iq = 1 , ip
            call rdedz(t1,nsq,ifort1)
            call rdedz(t2,nsq,ifort2)
            if (.not.uhf) then
               call dcopy(nsq,t1,1,t3,1)
            else
               call rdedz(t3,nsq,ifort3)
            end if
c
            do 30 i = 1 , ncoorb
               do 20 j = 1 , ncoorb
                  a1(i,j) = t1(i,j) + t1(i,j) - t2(i,j) - t2(j,i)
                  a2(i,j) = t3(i,j) + t3(i,j)
 20            continue
 30         continue
            call wtedz(a1,nsq,istrma)
            call vclr(a1,1,nsq)
            call vclr(a2,1,nsq)
c
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ib = ip
                  j = iq
                  ebj = e1(ib) - e1(j)
c
                  do 50 ia = nocc1 , ncoorb
                     do 40 i = 1 , nocc
                        difen = ebj + e1(ia) - e1(i)
                        a2(i,ia) = (t1(i,ia)-t2(i,ia))/difen
c
c
                        mfact = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        ij = (i-1)*i/2 + j
                        if (j.gt.i) then
                           ij = (j-1)*j/2 + i
                           mfact = -mfact
                        end if
                        iabij = (iab-1)*ijm + ij
                        a1(i,ia) = umm(iabij)/difen*4.0d0*mfact
 40                  continue
 50               continue
c
c   construct   v   matrix
c
c  (ac\jb)tiajb   and   (ki\jb)tkajb
c
c
                  call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
c
c
                  call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,
     +                      1,ncoorb,nocc,nvirt,nocc)
c
c   construct  zl  matrix
c
                  call mxmb(a1(1,nocc1),1,ncoorb,t1(nocc1,1),1,ncoorb,
     +                      zl,ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t1,1,ncoorb,a1(1,nocc1),1,ncoorb,zl(1,nocc1)
     +                      ,1,ncoorb,ncoorb,nocc,nvirt)
c
c
c   add on ump2 contributions
c
                  do 70 ia = nocc1 , ncoorb
                     do 60 i = 1 , nocc
                        difen = ebj + e1(ia) - e1(i)
                        a2(i,ia) = -a2(i,ia)
                        a1(i,ia) = t1(i,ia)/difen
 60                  continue
 70               continue
c
                  call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,
     +                      1,ncoorb,nocc,nvirt,nocc)
c
                  call mxmb(a2(1,nocc1),1,ncoorb,t1(nocc1,1),1,ncoorb,
     +                      zl,ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t1,1,ncoorb,a2(1,nocc1),1,ncoorb,zl(1,nocc1)
     +                      ,1,ncoorb,ncoorb,nocc,nvirt)
c
c
c
c  (ij\kc)sjakc
c
                  ic = ip
                  j = iq
                  do 130 ib = nocc1 , ncoorb
                     mfact = 1
                     ibc = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                     if (ib.gt.ic) then
                        ibc = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                        mfact = -mfact
                     end if
                     do 100 ia = nocc1 , ncoorb
                        mfact1 = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact1 = -mfact1
                        end if
c
c
c
                        do 90 l = 1 , nocc
                           lj = (j-1)*j/2 + l
                           mfact3 = mfact1
                           if (l.gt.j) then
                              lj = (l-1)*l/2 + j
                              mfact3 = -mfact1
                           end if
                           iablj = (iab-1)*ijm + lj
c
                           a2(l,ia) = amm(iablj)*mfact3
c
                           do 80 k = 1 , nocc
                              kl = (k-1)*k/2 + l
                              mfact2 = mfact
                              if (l.gt.k) then
                                 kl = (l-1)*l/2 + k
                                 mfact2 = -mfact
                              end if
                              ibckl = (ibc-1)*ijm + kl
c
                              a1(k,l) = amm(ibckl)*mfact2
 80                        continue
 90                     continue
 100                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,ncoorb,1,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,nocc,nvirt)
c
c   (ic\kj) sjakc
c
                     do 120 k = 1 , nocc
                        do 110 l = 1 , nocc
                           a1(k,l) = -a1(k,l)
 110                    continue
 120                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,1,ncoorb,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,nocc,nvirt)
c
 130              continue
c
c
c  (ij\kc)sjakc
c
                  do 190 ib = noccz1 , ncoorb
                     do 160 ia = nocc1 , ncoorb
                        do 150 k = 1 , nocc
                           do 140 l = 1 , noccz
                              ikc = (ic-nocc1)*nocc + k
                              ilb = (ib-noccz1)*noccz + l
                              ija = (ia-nocc1)*nocc + j
                              a1(k,l) = ads(ikc,ilb)
                              a2(l,ia) = ads(ija,ilb)
 140                       continue
 150                    continue
 160                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,ncoorb,1,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,noccz)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,noccz,nvirt)
c
c  (ic\kj)sjakc
c
                     do 180 k = 1 , nocc
                        do 170 l = 1 , noccz
                           a1(k,l) = -a1(k,l)
 170                    continue
 180                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,1,ncoorb,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,noccz)
c
c
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,noccz,nvirt)
c
 190              continue
               end if
c
c   db term
c
               if (iq.gt.nocc) then
                  ib = ip
                  id = iq
c
c  (ic\bd)sabcd
c
                  call vclr(dum,1,nsq)
                  do 220 l = 1 , nocc
                     do 210 k = 1 , nocc
                        mfact = 1
                        mfact1 = 1
                        kl = (k-1)*k/2 + l
                        if (l.gt.k) then
                           kl = (l-1)*l/2 + k
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 200 ic = nocc1 , ncoorb
                           ia = ic
                           mfact2 = mfact
                           idc = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                           if (id.gt.ic) then
                              idc = (id-nocc1)*(id-nocc)/2 + ic - nocc
                              mfact2 = -mfact
                           end if
                           klcd = (idc-1)*ijm + kl
                           mfact3 = mfact1
                           iba = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iba = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact3 = -mfact1
                           end if
                           klba = (iba-1)*ijm + kl
c
                           a1(k,ic) = amm(klcd)*mfact2
                           a2(k,ia) = amm(klba)*mfact3*0.5d0
c
c
 200                    continue
 210                 continue
c
                     call mxmb(a1(1,nocc1),ncoorb,1,a2(1,nocc1),1,
     +                         ncoorb,dum(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     if (ib.ne.id)
     +                   call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,
     +                   ncoorb,dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,
     +                   nvirt)
c
 220              continue
c
                  call mxmb(t1(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,
     +                      ncoorb,zl(1,nocc1),1,ncoorb,ncoorb,nvirt,
     +                      nvirt)
c
c  now bc term
c
                  ib = ip
                  ic = iq
c
c  now bc term
c
                  ib = ip
                  ic = iq
c
c   (ba\kc)sibkc
c
                  do 290 id = nocc1 , ncoorb
                     mfact = 1
                     mfact1 = 1
                     icd = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                     if (id.gt.ic) then
                        icd = (id-nocc1)*(id-nocc)/2 + ic - nocc
                        mfact = -mfact
                     end if
                     ibd = (ib-nocc1)*(ib-nocc)/2 + id - nocc
                     if (id.gt.ib) then
                        ibd = (id-nocc1)*(id-nocc)/2 + ib - nocc
                        mfact1 = -mfact1
                     end if
                     do 240 k = 1 , nocc
                        do 230 j = 1 , nocc
                           i = k
                           kj = (k-1)*k/2 + j
                           mfact2 = mfact
                           if (j.gt.k) then
                              kj = (j-1)*j/2 + k
                              mfact2 = -mfact
                           end if
                           icdkj = (icd-1)*ijm + kj
                           ij = (i-1)*i/2 + j
                           mfact3 = mfact1
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact3 = -mfact1
                           end if
                           ibdij = (ibd-1)*ijm + ij
                           a1(k,j) = amm(icdkj)*mfact2
                           a2(i,j) = amm(ibdij)*mfact3
c
 230                    continue
 240                 continue
c
c  (bc\ka)sibkc
c
                     call vclr(dum,1,nsq)
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,dum,1,ncoorb,
     +                         nocc,nocc,nocc)
c
c
c
                     call mxmb(t2,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c   (ba\kc)sibkc
c
                     do 260 k = 1 , nocc
                        do 250 j = 1 , nocc
                           dum(k,j) = -dum(k,j)
 250                    continue
 260                 continue
c
                     call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call vclr(dum,1,nsq)
                     if (ib.ne.ic) then
                        call mxmb(a2,1,ncoorb,a1,ncoorb,1,dum,1,ncoorb,
     +                            nocc,nocc,nocc)
c
                        call mxmb(dum,ncoorb,1,t2,1,ncoorb,zl,ncoorb,1,
     +                            nocc,nocc,ncoorb)
c
                        do 280 k = 1 , nocc
                           do 270 j = 1 , nocc
                              dum(k,j) = -dum(k,j)
 270                       continue
 280                    continue
c
                        call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                            ncoorb,nocc,nocc)
c
c
                     end if
 290              continue
c
                  do 360 id = noccz1 , ncoorb
                     do 310 i = 1 , nocc
                        do 300 j = 1 , noccz1
                           k = i
                           ibi = (ib-nocc1)*nocc + i
                           ikc = (ic-nocc1)*nocc + k
                           ijd = (id-noccz1)*noccz + j
                           a1(k,j) = ads(ikc,ijd)
                           a2(i,j) = ads(ibi,ijd)
c
 300                    continue
 310                 continue
c
c
                     call vclr(dum,1,nsq)
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,dum,1,ncoorb,
     +                         nocc,noccz,nocc)
c
                     call mxmb(t2,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c  (bc\ka)sibkc
c
                     do 330 k = 1 , nocc
                        do 320 j = 1 , nocc
                           dum(k,j) = -dum(k,j)
 320                    continue
 330                 continue
c
                     call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c
                     if (ib.ne.ic) then
                        call vclr(dum,1,nsq)
                        call mxmb(a2,1,ncoorb,a1,ncoorb,1,dum,1,ncoorb,
     +                            nocc,noccz,nocc)
c
                        call mxmb(dum,ncoorb,1,t2,1,ncoorb,zl,ncoorb,1,
     +                            nocc,nocc,ncoorb)
c
c
                        do 350 k = 1 , nocc
                           do 340 j = 1 , nocc
                              dum(k,j) = -dum(k,j)
 340                       continue
 350                    continue
c
                        call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                            ncoorb,nocc,nocc)
c
                     end if
c
 360              continue
               end if
            end if
c
c
c  now lj term
c
c   (ka\lj)sikjl
c
            if (ip.le.nocc) then
               if (iq.le.nocc) then
                  l = ip
                  j = iq
c
                  call vclr(dum,1,nsq)
                  do 390 id = nocc1 , ncoorb
                     do 380 ic = nocc1 , ncoorb
                        mfact = 1
                        mfact1 = 1
                        icd = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                        if (id.gt.ic) then
                           icd = (id-nocc1)*(id-nocc)/2 + ic - nocc
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 370 k = 1 , nocc
                           i = k
                           kl = (k-1)*k/2 + l
                           mfact2 = mfact
                           if (l.gt.k) then
                              kl = (l-1)*l/2 + k
                              mfact2 = -mfact
                           end if
                           icdkl = (icd-1)*ijm + kl
c
                           ij = (i-1)*i/2 + j
                           mfact3 = mfact1
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact3 = -mfact1
                           end if
                           icdij = (icd-1)*ijm + ij
c
                           a1(k,ic) = amm(icdkl)*mfact2
                           a2(i,ic) = amm(icdij)*mfact3*0.5d0
c
 370                    continue
 380                 continue
c
                     call mxmb(a1(1,nocc1),1,ncoorb,a2(1,nocc1),ncoorb,
     +                         1,dum,1,ncoorb,nocc,nvirt,nocc)
c
                     if (l.ne.j)
     +                   call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),
     +                   ncoorb,1,dum,1,ncoorb,nocc,nvirt,nocc)
c
 390              continue
c
                  call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                      nocc,nocc)
               end if
            end if
c
c
c  now beta terms - firstly  jb term
c
            if (ip.le.noccz) go to 600
            if (iq.le.noccz) then
               ib = ip
               j = iq
               ebj = e2(ib) - e2(j)
               ijb = (ib-noccz1)*noccz + j
c
c   (ac\jb)ticjb  and  (ki\jb)tkajb
c
               do 410 ia = nocc1 , ncoorb
                  do 400 i = 1 , nocc
                     difen = ebj + e1(ia) - e1(i)
                     iia = (ia-nocc1)*nocc + i
                     a2(i,ia) = t3(i,ia)/difen*2.0d0
                     a1(i,ia) = uds(iia,ijb)/difen*4.0d0
 400              continue
 410           continue
c
               call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                   y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
               call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,1,
     +                   ncoorb,nocc,nvirt,nocc)
c
               call mxmb(a1(1,nocc1),1,ncoorb,t3(nocc1,1),1,ncoorb,zl,
     +                   ncoorb,1,nocc,nvirt,ncoorb)
c
               call mxmb(t3,1,ncoorb,a1(1,nocc1),1,ncoorb,zl(1,nocc1),1,
     +                   ncoorb,ncoorb,nocc,nvirt)
c
c add on ump2 contributions
c
               do 430 ia = nocc1 , ncoorb
                  do 420 i = 1 , nocc
                     iai = (ia-nocc1)*nocc + i
                     zint(iai) = a2(i,ia)*0.5d0
                     a2(i,ia) = -a2(i,ia)*0.5d0
 420              continue
 430           continue
c
               call mxmb(zint,nocc,1,a2(1,nocc1),1,ncoorb,y(nocc1,nocc1)
     +                   ,1,ncoorb,nvirt,nocc,nvirt)
c
               call mxmb(zint,1,nocc,a2(1,nocc1),ncoorb,1,y,1,ncoorb,
     +                   nocc,nvirt,nocc)
c
               call mxmb(a2(1,nocc1),1,ncoorb,t3(nocc1,1),1,ncoorb,zl,
     +                   ncoorb,1,nocc,nvirt,ncoorb)
c
c
c
               call mxmb(t3,1,ncoorb,a2(1,nocc1),1,ncoorb,zl(1,nocc1),1,
     +                   ncoorb,ncoorb,nocc,nvirt)
c
c  now kc term
c
c   (ij\kc)sjakc
c
c
               ic = ip
               k = iq
               ikc = (ic-noccz1)*noccz + k
c
               call vclr(zint,1,mn)
               do 490 ia = nocc1 , ncoorb
                  do 480 j = 1 , nocc
                     do 450 ib = nocc1 , ncoorb
                        mfact = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        do 440 l = 1 , nocc
                           ilb = (ib-nocc1)*nocc + l
                           jl = (j-1)*j/2 + l
                           mfact2 = mfact
                           if (l.gt.j) then
                              jl = (l-1)*l/2 + j
                              mfact2 = -mfact
                           end if
                           iabjl = (iab-1)*ijm + jl
                           ija = (ia-nocc1)*nocc + j
c
                           zint(ija) = zint(ija) + amm(iabjl)
     +                                 *ads(ilb,ikc)*mfact2
 440                    continue
 450                 continue
c
c
                     do 470 ib = noccz1 , ncoorb
                        mfact = 1
                        icb = (ic-noccz1)*(ic-noccz)/2 + ib - noccz
                        if (ib.gt.ic) then
                           icb = (ib-noccz1)*(ib-noccz)/2 + ic - noccz
                           mfact = -mfact
                        end if
                        do 460 l = 1 , noccz
                           ilb = (ib-noccz1)*noccz + l
                           ija = (ia-nocc1)*nocc + j
                           kl = (k-1)*k/2 + l
                           mfact2 = mfact
                           if (l.gt.k) then
                              kl = (l-1)*l/2 + k
                              mfact2 = -mfact
                           end if
                           icbkl = (icb-1)*ijmz + kl
c
                           zint(ija) = zint(ija) + ads(ija,ilb)
     +                                 *ammz(icbkl)*mfact2
 460                    continue
 470                 continue
 480              continue
 490           continue
c
               call mxmb(t3,1,ncoorb,zint,1,nocc,zl(1,nocc1),1,ncoorb,
     +                   ncoorb,nocc,nvirt)
c
               call mxmb(t3(1,nocc1),1,ncoorb,zint,nocc,1,zl,1,ncoorb,
     +                   ncoorb,nvirt,nocc)
c
            end if
c
c  now bd term
c
c  (ic\bd)sabcd
c
            if (iq.gt.noccz) then
               ib = ip
               id = iq
               call vclr(dum,1,nsq)
c
               do 570 l = 1 , noccz
                  ibl = (ib-noccz1)*noccz + l
                  idl = (id-noccz1)*noccz + l
                  ic = id
                  ibj = ibl
                  icj = idl
c
                  call mxmb(ads(1,ibl),nocc,1,ads(1,idl),1,nocc,
     +                      dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(ads(1,ibj),1,nocc,ads(1,icj),nocc,1,dum,1,
     +                      ncoorb,nocc,nvirt,nocc)
c
                  if (ib.ne.id) then
                     call mxmb(ads(1,idl),nocc,1,ads(1,ibl),1,nocc,
     +                         dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,
     +                         nvirt)
c
                     call mxmb(ads(1,icj),1,nocc,ads(1,ibj),nocc,1,dum,
     +                         1,ncoorb,nocc,nvirt,nocc)
                  end if
c
 570           continue
c
               call mxmb(t3(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,ncoorb,
     +                   zl(1,nocc1),1,ncoorb,ncoorb,nvirt,nvirt)
c
c   (bc/ka)sibkc
c
               do 590 i = 1 , nocc
                  do 580 j = 1 , nocc
                     dum(i,j) = -dum(i,j)
 580              continue
 590           continue
c
               call mxmb(t3,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                   nocc,nocc)
            end if
c
c  now  lj term
c
c   (ka\lj)sijkl
c
 600        if (ip.le.noccz) then
               if (iq.le.noccz) then
                  l = ip
                  j = iq
                  call vclr(dum,1,nsq)
                  do 610 id = noccz1 , ncoorb
                     ild = (id-noccz1)*noccz + l
                     ijd = (id-noccz1)*noccz + j
                     k = l
                     ib = id
                     ijb = ijd
                     ikb = ild
c
                     call mxmb(ads(1,ild),1,nocc,ads(1,ijd),nocc,1,dum,
     +                         1,ncoorb,nocc,nvirt,nocc)
c
                     call mxmb(ads(1,ijb),nocc,1,ads(1,ikb),1,nocc,
     +                         dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,
     +                         nvirt)
c
                     if (l.ne.j) then
                        call mxmb(ads(1,ijd),1,nocc,ads(1,ild),nocc,1,
     +                            dum,1,ncoorb,nocc,nvirt,nocc)
c
                        call mxmb(ads(1,ikb),nocc,1,ads(1,ijb),1,nocc,
     +                            dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,
     +                            nvirt)
                     end if
c
 610              continue
c
                  call mxmb(t3,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                      nocc,nocc)
c
c  (ic\kj)sjakc
c
                  do 630 ia = nocc1 , ncoorb
                     do 620 ib = nocc1 , ncoorb
                        dum(ia,ib) = -dum(ia,ib)
 620                 continue
 630              continue
                  call mxmb(t3(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,
     +                      ncoorb,zl(1,nocc1),1,ncoorb,ncoorb,nvirt,
     +                      nvirt)
               end if
            end if
c
c
 640     continue
 650  continue
c
c
c
      do 670 i = 1 , nocc
         do 660 j = 1 , nocc
            y(i,j) = -y(i,j)
 660     continue
 670  continue
c
      call wrt3(y,nsq,ibly,ifile1)
      call wrt3(zl,nsq,iblw,ifile1)
      return
      end
      subroutine umpmkb(umm,amm,ammz,uds,ads,t1,t2,t3,a1,a2,e1,e2,
     1               y,zl,dum,zint,nocc,noccz,ndep,ndepz,mn,mnz,ncoorb,
     1                ifort1,ifort2,ifort3,istrma,ibly,iblw)
c
      implicit real*8  (a-h,o-z)
      dimension umm(ndep),amm(ndep),ammz(ndepz),
     +          ads(mnz,mn),uds(mnz,mn),t1(ncoorb,ncoorb),
     +          t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),a1(ncoorb,ncoorb),
     +          a2(ncoorb,ncoorb),zint(mn),e1(ncoorb),e2(ncoorb),
     +          y(ncoorb,ncoorb),zl(ncoorb,ncoorb),dum(ncoorb,ncoorb)
c
      ifile1 = 1
      nocc1 = nocc + 1
      nvirt = ncoorb - nocc
c     nvirtz = ncoorb - noccz
      noccz1 = noccz + 1
      nsq = ncoorb*ncoorb
      ijm = (nocc+1)*nocc/2
      ijmz = (noccz+1)*noccz/2
      call vclr(y,1,nsq)
      call vclr(zl,1,nsq)
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      call rewedz(istrma)
c
      do 650 ip = 1 , ncoorb
         do 640 iq = 1 , ip
            call rdedz(t1,nsq,ifort1)
            call rdedz(t2,nsq,ifort2)
            call rdedz(t3,nsq,ifort3)
c
            do 30 i = 1 , ncoorb
               do 20 j = 1 , ncoorb
                  a1(i,j) = t1(i,j) + t1(i,j) - t2(i,j) - t2(j,i)
                  a2(i,j) = t3(i,j) + t3(i,j)
 20            continue
 30         continue
            call wtedz(a1,nsq,istrma)
            call vclr(a1,1,nsq)
            call vclr(a2,1,nsq)
c
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ib = ip
                  j = iq
                  ebj = e1(ib) - e1(j)
c
                  do 50 ia = nocc1 , ncoorb
                     do 40 i = 1 , nocc
                        difen = ebj + e1(ia) - e1(i)
                        a2(i,ia) = (t1(i,ia)-t2(i,ia))/difen
c
                        mfact = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        ij = (i-1)*i/2 + j
                        if (j.gt.i) then
                           ij = (j-1)*j/2 + i
                           mfact = -mfact
                        end if
                        iabij = (iab-1)*ijm + ij
                        a1(i,ia) = umm(iabij)/difen*4.0d0*mfact
 40                  continue
 50               continue
c
c   construct   v   matrix
c
c  (ac\jb)tiajb   and   (ki\jb)tkajb
c
c
                  call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,
     +                      1,ncoorb,nocc,nvirt,nocc)
c
c   construct  zl  matrix
c
                  call mxmb(a1(1,nocc1),1,ncoorb,t1(nocc1,1),1,ncoorb,
     +                      zl,ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t1,1,ncoorb,a1(1,nocc1),1,ncoorb,zl(1,nocc1)
     +                      ,1,ncoorb,ncoorb,nocc,nvirt)
c
c   add on ump2 contribution
c
                  do 70 ia = nocc1 , ncoorb
                     do 60 i = 1 , nocc
                        difen = ebj + e1(ia) - e1(i)
                        a2(i,ia) = -a2(i,ia)
                        a1(i,ia) = t1(i,ia)/difen
 60                  continue
 70               continue
c
                  call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,
     +                      1,ncoorb,nocc,nvirt,nocc)
c
                  call mxmb(a2(1,nocc1),1,ncoorb,t1(nocc1,1),1,ncoorb,
     +                      zl,ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t1,1,ncoorb,a2(1,nocc1),1,ncoorb,zl(1,nocc1)
     +                      ,1,ncoorb,ncoorb,nocc,nvirt)
c
c  (ij\kc)sjakc
c
                  ic = ip
                  j = iq
                  do 130 ib = nocc1 , ncoorb
                     mfact = 1
                     ibc = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                     if (ib.gt.ic) then
                        ibc = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                        mfact = -mfact
                     end if
                     do 100 ia = nocc1 , ncoorb
                        mfact1 = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact1 = -mfact1
                        end if
c
                        do 90 l = 1 , nocc
                           lj = (j-1)*j/2 + l
                           mfact3 = mfact1
                           if (l.gt.j) then
                              lj = (l-1)*l/2 + j
                              mfact3 = -mfact1
                           end if
                           iablj = (iab-1)*ijm + lj
c
                           a2(l,ia) = amm(iablj)*mfact3
c
                           do 80 k = 1 , nocc
                              kl = (k-1)*k/2 + l
                              mfact2 = mfact
                              if (l.gt.k) then
                                 kl = (l-1)*l/2 + k
                                 mfact2 = -mfact
                              end if
                              ibckl = (ibc-1)*ijm + kl
c
                              a1(k,l) = amm(ibckl)*mfact2
 80                        continue
 90                     continue
 100                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,ncoorb,1,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,nocc,nvirt)
c
c   (ic\kj) sjakc
c
                     do 120 k = 1 , nocc
                        do 110 l = 1 , nocc
                           a1(k,l) = -a1(k,l)
 110                    continue
 120                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,1,ncoorb,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,nocc,nvirt)
c
 130              continue
c
c  (ij\kc)sjakc
c
                  do 190 ib = noccz1 , ncoorb
                     do 160 ia = nocc1 , ncoorb
                        do 150 k = 1 , nocc
                           do 140 l = 1 , noccz
                              ikc = (ic-nocc1)*nocc + k
                              ilb = (ib-noccz1)*noccz + l
                              ija = (ia-nocc1)*nocc + j
                              a1(k,l) = ads(ilb,ikc)
                              a2(l,ia) = ads(ilb,ija)
 140                       continue
 150                    continue
 160                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,ncoorb,1,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,noccz)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,noccz,nvirt)
c
c  (ic\kj)sjakc
c
                     do 180 k = 1 , nocc
                        do 170 l = 1 , noccz
                           a1(k,l) = -a1(k,l)
 170                    continue
 180                 continue
c
                     call vclr(dum,1,nsq)
                     call mxmb(t2,1,ncoorb,a1,1,ncoorb,dum,1,ncoorb,
     +                         ncoorb,nocc,noccz)
c
                     call mxmb(dum,1,ncoorb,a2(1,nocc1),1,ncoorb,
     +                         zl(1,nocc1),1,ncoorb,ncoorb,noccz,nvirt)
c
 190              continue
               end if
c
c   db term
c
               if (iq.gt.nocc) then
                  ib = ip
                  id = iq
c
c  (ic\bd)sabcd
c
                  call vclr(dum,1,nsq)
                  do 220 l = 1 , nocc
                     do 210 k = 1 , nocc
                        mfact = 1
                        mfact1 = 1
                        kl = (k-1)*k/2 + l
                        if (l.gt.k) then
                           kl = (l-1)*l/2 + k
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 200 ic = nocc1 , ncoorb
                           ia = ic
                           mfact2 = mfact
                           idc = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                           if (id.gt.ic) then
                              idc = (id-nocc1)*(id-nocc)/2 + ic - nocc
                              mfact2 = -mfact
                           end if
                           klcd = (idc-1)*ijm + kl
                           mfact3 = mfact1
                           iba = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iba = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact3 = -mfact1
                           end if
                           klba = (iba-1)*ijm + kl
c
                           a1(k,ic) = amm(klcd)*mfact2
                           a2(k,ia) = amm(klba)*mfact3*0.5d0
c
 200                    continue
 210                 continue
c
                     call mxmb(a1(1,nocc1),ncoorb,1,a2(1,nocc1),1,
     +                         ncoorb,dum(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     if (ib.ne.id)
     +                   call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,
     +                   ncoorb,dum(nocc1,nocc1),1,ncoorb,nvirt,nocc,
     +                   nvirt)
c
 220              continue
c
                  call mxmb(t1(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,
     +                      ncoorb,zl(1,nocc1),1,ncoorb,ncoorb,nvirt,
     +                      nvirt)
c
c  now bc term
c
                  ib = ip
                  ic = iq
c
c  now bc term
c
                  ib = ip
                  ic = iq
c
c   (ba\kc)sibkc
c
                  do 290 id = nocc1 , ncoorb
                     mfact = 1
                     mfact1 = 1
                     icd = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                     if (id.gt.ic) then
                        icd = (id-nocc1)*(id-nocc)/2 + ic - nocc
                        mfact = -mfact
                     end if
                     ibd = (ib-nocc1)*(ib-nocc)/2 + id - nocc
                     if (id.gt.ib) then
                        ibd = (id-nocc1)*(id-nocc)/2 + ib - nocc
                        mfact1 = -mfact1
                     end if
                     do 240 k = 1 , nocc
                        do 230 j = 1 , nocc
                           i = k
                           kj = (k-1)*k/2 + j
                           mfact2 = mfact
                           if (j.gt.k) then
                              kj = (j-1)*j/2 + k
                              mfact2 = -mfact
                           end if
                           icdkj = (icd-1)*ijm + kj
                           ij = (i-1)*i/2 + j
                           mfact3 = mfact1
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact3 = -mfact1
                           end if
                           ibdij = (ibd-1)*ijm + ij
                           a1(k,j) = amm(icdkj)*mfact2
                           a2(i,j) = amm(ibdij)*mfact3
c
 230                    continue
 240                 continue
c
c  (bc\ka)sibkc
c
                     call vclr(dum,1,nsq)
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,dum,1,ncoorb,
     +                         nocc,nocc,nocc)
c
                     call mxmb(t2,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c   (ba\kc)sibkc
c
                     do 260 k = 1 , nocc
                        do 250 j = 1 , nocc
                           dum(k,j) = -dum(k,j)
 250                    continue
 260                 continue
c
                     call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
                     call vclr(dum,1,nsq)
                     if (ib.ne.ic) then
                        call mxmb(a2,1,ncoorb,a1,ncoorb,1,dum,1,ncoorb,
     +                            nocc,nocc,nocc)
c
                        call mxmb(dum,ncoorb,1,t2,1,ncoorb,zl,ncoorb,1,
     +                            nocc,nocc,ncoorb)
c
                        do 280 k = 1 , nocc
                           do 270 j = 1 , nocc
                              dum(k,j) = -dum(k,j)
 270                       continue
 280                    continue
c
                        call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                            ncoorb,nocc,nocc)
c
c
                     end if
 290              continue
c
                  do 360 id = noccz1 , ncoorb
                     do 310 i = 1 , nocc
c     do 21 k=1,nocc
                        do 300 j = 1 , noccz1
                           k = i
                           ibi = (ib-nocc1)*nocc + i
                           ikc = (ic-nocc1)*nocc + k
                           ijd = (id-noccz1)*noccz + j
                           a1(k,j) = ads(ijd,ikc)
                           a2(i,j) = ads(ijd,ibi)
c
c
 300                    continue
 310                 continue
c
                     call vclr(dum,1,nsq)
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,dum,1,ncoorb,
     +                         nocc,noccz,nocc)
c
                     call mxmb(t2,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c  (bc\ka)sibkc
c
                     do 330 k = 1 , nocc
                        do 320 j = 1 , nocc
                           dum(k,j) = -dum(k,j)
 320                    continue
 330                 continue
c
                     call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                         ncoorb,nocc,nocc)
c
c
                     if (ib.ne.ic) then
                        call vclr(dum,1,nsq)
                        call mxmb(a2,1,ncoorb,a1,ncoorb,1,dum,1,ncoorb,
     +                            nocc,noccz,nocc)
c
                        call mxmb(dum,ncoorb,1,t2,1,ncoorb,zl,ncoorb,1,
     +                            nocc,nocc,ncoorb)
c
                        do 350 k = 1 , nocc
                           do 340 j = 1 , nocc
                              dum(k,j) = -dum(k,j)
 340                       continue
 350                    continue
c
                        call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,
     +                            ncoorb,nocc,nocc)
c
                     end if
c
 360              continue
               end if
            end if
c
c  now lj term
c
c   (ka\lj)sikjl
c
            if (ip.le.nocc) then
               if (iq.le.nocc) then
                  l = ip
                  j = iq
c
                  call vclr(dum,1,nsq)
                  do 390 id = nocc1 , ncoorb
                     do 380 ic = nocc1 , ncoorb
                        mfact = 1
                        mfact1 = 1
                        icd = (ic-nocc1)*(ic-nocc)/2 + id - nocc
                        if (id.gt.ic) then
                           icd = (id-nocc1)*(id-nocc)/2 + ic - nocc
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 370 k = 1 , nocc
                           i = k
                           kl = (k-1)*k/2 + l
                           mfact2 = mfact
                           if (l.gt.k) then
                              kl = (l-1)*l/2 + k
                              mfact2 = -mfact
                           end if
                           icdkl = (icd-1)*ijm + kl
c
                           ij = (i-1)*i/2 + j
                           mfact3 = mfact1
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact3 = -mfact1
                           end if
                           icdij = (icd-1)*ijm + ij
c
                           a1(k,ic) = amm(icdkl)*mfact2
                           a2(i,ic) = amm(icdij)*mfact3*0.5d0
c
 370                    continue
 380                 continue
c
                     call mxmb(a1(1,nocc1),1,ncoorb,a2(1,nocc1),ncoorb,
     +                         1,dum,1,ncoorb,nocc,nvirt,nocc)
c
                     if (l.ne.j)
     +                   call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),
     +                   ncoorb,1,dum,1,ncoorb,nocc,nvirt,nocc)
c
 390              continue
c
                  call mxmb(t1,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                      nocc,nocc)
               end if
            end if
c
c  now beta terms - firstly  jb term
c
            if (ip.le.noccz) go to 600
            if (iq.le.noccz) then
               ib = ip
               j = iq
               ebj = e2(ib) - e2(j)
               ijb = (ib-noccz1)*noccz + j
c
c   (ac\jb)ticjb  and  (ki\jb)tkajb
c
               do 410 ia = nocc1 , ncoorb
                  do 400 i = 1 , nocc
                     difen = ebj + e1(ia) - e1(i)
                     iia = (ia-nocc1)*nocc + i
                     a2(i,ia) = t3(i,ia)/difen*2.0d0
                     a1(i,ia) = uds(ijb,iia)/difen*4.0d0
 400              continue
 410           continue
c
               call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,ncoorb,
     +                   y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
               call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),ncoorb,1,y,1,
     +                   ncoorb,nocc,nvirt,nocc)
c
               call mxmb(a1(1,nocc1),1,ncoorb,t3(nocc1,1),1,ncoorb,zl,
     +                   ncoorb,1,nocc,nvirt,ncoorb)
c
               call mxmb(t3,1,ncoorb,a1(1,nocc1),1,ncoorb,zl(1,nocc1),1,
     +                   ncoorb,ncoorb,nocc,nvirt)
c
c  add on ump2 contributions
c
               do 430 ia = nocc1 , ncoorb
                  do 420 i = 1 , nocc
                     iai = (ia-nocc1)*nocc + i
                     zint(iai) = a2(i,ia)*0.5d0
                     a2(i,ia) = -a2(i,ia)*0.5d0
 420              continue
 430           continue
c
               call mxmb(zint,nocc,1,a2(1,nocc1),1,ncoorb,y(nocc1,nocc1)
     +                   ,1,ncoorb,nvirt,nocc,nvirt)
c
               call mxmb(zint,1,nocc,a2(1,nocc1),ncoorb,1,y,1,ncoorb,
     +                   nocc,nvirt,nocc)
c
               call mxmb(a2(1,nocc1),1,ncoorb,t3(nocc1,1),1,ncoorb,zl,
     +                   ncoorb,1,nocc,nvirt,ncoorb)
c
               call mxmb(t3,1,ncoorb,a2(1,nocc1),1,ncoorb,zl(1,nocc1),1,
     +                   ncoorb,ncoorb,nocc,nvirt)
c
c  now kc term
c
c   (ij\kc)sjakc
c
               ic = ip
               k = iq
               ikc = (ic-noccz1)*noccz + k
c
               call vclr(zint,1,mn)
               do 490 ia = nocc1 , ncoorb
                  do 480 j = 1 , nocc
                     do 450 ib = nocc1 , ncoorb
                        mfact = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        do 440 l = 1 , nocc
                           ilb = (ib-nocc1)*nocc + l
                           jl = (j-1)*j/2 + l
                           mfact2 = mfact
                           if (l.gt.j) then
                              jl = (l-1)*l/2 + j
                              mfact2 = -mfact
                           end if
                           iabjl = (iab-1)*ijm + jl
                           ija = (ia-nocc1)*nocc + j
c
                           zint(ija) = zint(ija) + amm(iabjl)
     +                                 *ads(ikc,ilb)*mfact2
 440                    continue
 450                 continue
c
                     do 470 ib = noccz1 , ncoorb
                        mfact = 1
                        icb = (ic-noccz1)*(ic-noccz)/2 + ib - noccz
                        if (ib.gt.ic) then
                           icb = (ib-noccz1)*(ib-noccz)/2 + ic - noccz
                           mfact = -mfact
                        end if
                        do 460 l = 1 , noccz
                           ilb = (ib-noccz1)*noccz + l
                           ija = (ia-nocc1)*nocc + j
                           kl = (k-1)*k/2 + l
                           mfact2 = mfact
                           if (l.gt.k) then
                              kl = (l-1)*l/2 + k
                              mfact2 = -mfact
                           end if
                           icbkl = (icb-1)*ijmz + kl
c
                           zint(ija) = zint(ija) + ads(ilb,ija)
     +                                 *ammz(icbkl)*mfact2
 460                    continue
 470                 continue
 480              continue
 490           continue
c
               call mxmb(t3,1,ncoorb,zint,1,nocc,zl(1,nocc1),1,ncoorb,
     +                   ncoorb,nocc,nvirt)
c
               call mxmb(t3(1,nocc1),1,ncoorb,zint,nocc,1,zl,1,ncoorb,
     +                   ncoorb,nvirt,nocc)
            end if
c
c  now bd term
c
c  (ic\bd)sabcd
c
c
            if (iq.gt.noccz) then
               ib = ip
               id = iq
               call vclr(dum,1,nsq)
c
               do 570 l = 1 , noccz
                  ibl = (ib-noccz1)*noccz + l
                  idl = (id-noccz1)*noccz + l
                  ic = id
                  ibj = ibl
                  icj = idl
c
                  call mxmb(ads(ibl,1),mnz*nocc,mnz,ads(idl,1),mnz,
     +                      mnz*nocc,dum(nocc1,nocc1),1,ncoorb,nvirt,
     +                      nocc,nvirt)
c
                  call mxmb(ads(ibj,1),mnz,mnz*nocc,ads(icj,1),mnz*nocc,
     +                      mnz,dum,1,ncoorb,nocc,nvirt,nocc)
c
                  if (ib.ne.id) then
                     call mxmb(ads(idl,1),mnz*nocc,mnz,ads(ibl,1),mnz,
     +                         mnz*nocc,dum(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     call mxmb(ads(icj,1),mnz,mnz*nocc,ads(ibj,1),
     +                         mnz*nocc,mnz,dum,1,ncoorb,nocc,nvirt,
     +                         nocc)
                  end if
c
 570           continue
c
c
               call mxmb(t3(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,ncoorb,
     +                   zl(1,nocc1),1,ncoorb,ncoorb,nvirt,nvirt)
c
c   (bc/ka)sibkc
c
               do 590 i = 1 , nocc
                  do 580 j = 1 , nocc
                     dum(i,j) = -dum(i,j)
 580              continue
 590           continue
c
               call mxmb(t3,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                   nocc,nocc)
            end if
c
c  now  lj term
c
c   (ka\lj)sijkl
c
 600        if (ip.le.noccz) then
               if (iq.le.noccz) then
                  l = ip
                  j = iq
                  call vclr(dum,1,nsq)
                  do 610 id = noccz1 , ncoorb
                     ild = (id-noccz1)*noccz + l
                     ijd = (id-noccz1)*noccz + j
                     k = l
                     ib = id
                     ijb = ijd
                     ikb = ild
c
                     call mxmb(ads(ild,1),mnz,mnz*nocc,ads(ijd,1),
     +                         mnz*nocc,mnz,dum,1,ncoorb,nocc,nvirt,
     +                         nocc)
c
                     call mxmb(ads(ijb,1),mnz*nocc,mnz,ads(ikb,1),mnz,
     +                         mnz*nocc,dum(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     if (l.ne.j) then
                        call mxmb(ads(ijd,1),mnz,mnz*nocc,ads(ild,1),
     +                            mnz*nocc,mnz,dum,1,ncoorb,nocc,nvirt,
     +                            nocc)
c
                        call mxmb(ads(ikb,1),mnz*nocc,mnz,ads(ijb,1),
     +                            mnz,mnz*nocc,dum(nocc1,nocc1),1,
     +                            ncoorb,nvirt,nocc,nvirt)
                     end if
c
 610              continue
c
                  call mxmb(t3,1,ncoorb,dum,1,ncoorb,zl,1,ncoorb,ncoorb,
     +                      nocc,nocc)
c
c  (ic\kj)sjakc
c
                  do 630 ia = nocc1 , ncoorb
                     do 620 ib = nocc1 , ncoorb
                        dum(ia,ib) = -dum(ia,ib)
 620                 continue
 630              continue
                  call mxmb(t3(1,nocc1),1,ncoorb,dum(nocc1,nocc1),1,
     +                      ncoorb,zl(1,nocc1),1,ncoorb,ncoorb,nvirt,
     +                      nvirt)
               end if
            end if
c
 640     continue
 650  continue
c
c
      do 670 i = 1 , nocc
         do 660 j = 1 , nocc
            y(i,j) = -y(i,j)
 660     continue
 670  continue
c
      call wrt3(y,nsq,ibly,ifile1)
      call wrt3(zl,nsq,iblw,ifile1)
      return
      end
      subroutine squrh(t,sq,n)
      implicit real*8  (a-h,o-z)
      dimension sq(*),t(*)
      ij = 0
      ii = 0
      do 30 i = 1 , n
         jj = 0
         do 20 j = 1 , i
            ij = ij + 1
            zfact = 0.5d0
            if (i.eq.j) zfact = 1.0d0
            sq(ii+j) = t(ij)*zfact
            sq(jj+i) = t(ij)*zfact
            jj = jj + n
 20      continue
         ii = ii + n
 30   continue
      return
      end
      subroutine umps5(iso,nshels,ys,ya,yb,cssa,cssb,gtp,veca,vecb,
     +     dum,zinta,zint1a,zint2a,zintb,zint1b,zint2b,
     +     gtps,ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,
     +     nocca,noccb,ncoorb,nova,novb,nijas,nijbs,nabas,nabbs,ifytr,
     +     ifdma,ifdmb,iblw,ifw,iblkya,iblkyb,uhf)
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      logical uhf
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common/craypk/labout(1360)
      common/blk1/g(510),nint,nxtr
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      dimension m0(48),iso(nshels,*)
      logical lab,labc,labcd
      logical ijump,jjump
      dimension ys(ncoorb*ncoorb),ya(ncoorb*ncoorb),yb(ncoorb*ncoorb),
     +  cssa(ncoorb*ncoorb),cssb(ncoorb*ncoorb),
     +  veca(ncoorb,ncoorb),vecb(ncoorb,ncoorb),gtp(ncoorb,ncoorb),
     +  dum(ncoorb,ncoorb),zinta(nova),zint1a(nabas),zint2a(nijas),
     +  zintb(novb),zint1b(nabbs),zint2b(nijbs),gtps(*)
      data  cut/1.0d-10/
      data mzero/0/
c
      call izero(1360,labout,1)
      ntri = ncoorb*(ncoorb+1)/2
      nsq = ncoorb*ncoorb
      nocca1 = nocca + 1
      noccb1 = noccb + 1
      nija = (nocca+1)*nocca/2
      nijb = (noccb+1)*noccb/2
      nvirta = ncoorb - nocca
      nvirtb = ncoorb - noccb
      naba = (nvirta+1)*nvirta/2
      nabb = (nvirtb+1)*nvirtb/2
      fact = 2.0d0
      if (uhf) fact = 1.0d0
      ifile1 = 1
      call rdedx(ys,nsq,ifytr,ifile1)
      call rdedx(ya,nsq,iblkya,ifile1)
      call rdedx(cssa,nsq,ifdma,ifile1)
      call rdedx(yb,nsq,iblkyb,ifile1)
      call rdedx(cssb,nsq,ifdmb,ifile1)
c
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      if (uhf) then
         call rewedz(ifort4)
         call rewedz(ifort5)
         call rewedz(ifort6)
      end if
      call search(iblw,ifw)
      nint = 0
      do 160 ii = 1 , nshell
         do 30 it = 1 , nt
            id = iso(ii,it)
            ijump = id.gt.ii
            m0(it) = id
 30      continue
         iceni = katom(ii)
         do 150 jj = 1 , ii
            if (.not.(ijump)) then
               do 50 it = 1 , nt
                  id = m0(it)
                  jd = iso(jj,it)
                  jjump = jd.gt.ii
                  if (id.lt.jd) then
                     nd = id
                     id = jd
                     jd = nd
                  end if
                  jjump = id.eq.ii .and. jd.gt.jj
 50            continue
               lab = katom(jj).eq.iceni
            end if
            mini = kmin(ii)
            minj = kmin(jj)
            maxi = kmax(ii)
            maxj = kmax(jj)
            loci = kloc(ii) - mini
            locj = kloc(jj) - minj
            ntimes = (maxi-mini+1)*(maxj-minj+1)
            imax = loci + maxi
c
            do 80 itimes = 1 , ntimes
               ib1 = (itimes-1)*ntri
               call rdedz(zinta,nova,ifort1)
               call rdedz(zint1a,naba,ifort2)
               call rdedz(zint2a,nija,ifort3)
               if (uhf) then
                  call rdedz(zintb,novb,ifort4)
                  call rdedz(zint1b,nabb,ifort5)
                  call rdedz(zint2b,nijb,ifort6)
               end if
c
               call squrh(zint1a,dum,nvirta)
               call dcopy(nabas,dum,1,zint1a,1)
               call vclr(dum,1,nsq)
               call vclr(dum,1,nsq)
               call squrh(zint2a,dum,nocca)
               call dcopy(nijas,dum,1,zint2a,1)
c
               if (uhf) then
                  call squrh(zint1b,dum,nvirtb)
                  call dcopy(nabbs,dum,1,zint1b,1)
                  call squrh(zint2b,dum,noccb)
                  call dcopy(nijbs,dum,1,zint2b,1)
               end if
c
c
               if (.not.(ijump .or. jjump)) then
c
c                 iv1 = ncoorb*nocca + 1
c                 iv2 = ncoorb*noccb + 1
                  call vclr(dum,1,nsq)
                  call vclr(gtp,1,nsq)
c
                  call mxmb(veca(1,nocca1),1,ncoorb,zinta,nocca,1,dum,1,
     +                      ncoorb,imax,nvirta,nocca)
c
                  call mxmb(veca,1,ncoorb,dum,ncoorb,1,gtp,1,ncoorb,
     +                      imax,nocca,imax)
c
                  call vclr(dum,1,nsq)
                  call mxmb(veca(1,nocca1),1,ncoorb,zint1a,1,nvirta,dum,
     +                      1,ncoorb,imax,nvirta,nvirta)
c
                  call mxmb(veca(1,nocca1),1,ncoorb,dum,ncoorb,1,gtp,1,
     +                      ncoorb,imax,nvirta,imax)
c
                  call vclr(dum,1,nsq)
                  call mxmb(veca,1,ncoorb,zint2a,1,nocca,dum,1,ncoorb,
     +                      imax,nocca,nocca)
c
                  call mxmb(veca,1,ncoorb,dum,ncoorb,1,gtp,1,ncoorb,
     +                      imax,nocca,imax)
c
                  if (uhf) then
                     call vclr(dum,1,nsq)
                     call mxmb(vecb(1,noccb1),1,ncoorb,zintb,noccb,1,
     +                         dum,1,ncoorb,imax,nvirtb,noccb)
c
                     call mxmb(vecb,1,ncoorb,dum,ncoorb,1,gtp,1,ncoorb,
     +                         imax,noccb,imax)
c
                     call vclr(dum,1,nsq)
                     call mxmb(vecb(1,noccb1),1,ncoorb,zint1b,1,nvirtb,
     +                         dum,1,ncoorb,imax,nvirtb,nvirtb)
c
                     call mxmb(vecb(1,noccb1),1,ncoorb,dum,ncoorb,1,gtp,
     +                         1,ncoorb,imax,nvirtb,imax)
c
                     call vclr(dum,1,nsq)
                     call mxmb(vecb,1,ncoorb,zint2b,1,noccb,dum,1,
     +                         ncoorb,imax,noccb,noccb)
c
                     call mxmb(vecb,1,ncoorb,dum,ncoorb,1,gtp,1,ncoorb,
     +                         imax,noccb,imax)
                  end if
c
                  do 70 ms1 = 1 , imax
                     do 60 ms2 = 1 , ms1
                        ms12 = iky(ms1) + ms2
                        gtps(ib1+ms12) = (gtp(ms1,ms2)+gtp(ms2,ms1))
     +                     *0.5d0
 60                  continue
 70               continue
               end if
 80         continue
c
            if (.not.(ijump .or. jjump)) then
c
               do 140 kk = 1 , ii
                  labc = lab .and. katom(kk).eq.iceni
                  maxll = kk
                  if (kk.eq.ii) maxll = jj
                  do 130 ll = 1 , maxll
                     labcd = labc .and. katom(ll).eq.iceni
                     if (.not.(labcd)) then
                        mink = kmin(kk)
                        minl = kmin(ll)
                        maxk = kmax(kk)
                        maxl = kmax(ll)
                        lock = kloc(kk) - mink
                        locl = kloc(ll) - minl
                        icount = 1
                        do 120 i = mini , maxi
                           i1 = loci + i
                           do 110 j = minj , maxj
                              j1 = locj + j
                              do 100 k = mink , maxk
                                 k1 = lock + k
                                 do 90 l = minl , maxl
                                    l1 = locl + l
c
                                    if (k1.ge.l1) then
                                       ikl = iky(k1) + l1
                                    else
                                       ikl = iky(l1) + k1
                                    end if
                                    ipos = (icount-1)*ntri
                                    val = gtps(ipos+ikl)*fact
c
c                                   i11 = (i1-1)*ncoorb
                                    j11 = (j1-1)*ncoorb
                                    k11 = (k1-1)*ncoorb
                                    l11 = (l1-1)*ncoorb
c
                                    imn = j11 + i1
                                    ils = l11 + k1
                                    iml = k11 + i1
                                    ims = l11 + i1
                                    inl = k11 + j1
                                    ins = l11 + j1
c
                                    vv1 = ys(imn)*cssa(ils)
                                    vv2 = ys(ils)*cssa(imn)
                                    vv3 = ys(imn)*cssb(ils)
                                    vv4 = ys(ils)*cssb(imn)
                                    val = val + 0.5d0*(vv1+vv2+vv3+vv4)
c
                                    vv1 = ya(iml)*cssa(ins)
                                    vv2 = ya(ims)*cssa(inl)
                                    vv3 = ya(inl)*cssa(ims)
                                    vv4 = ya(ins)*cssa(iml)
                                    vv5 = yb(iml)*cssb(ins)
                                    vv6 = yb(ims)*cssb(inl)
                                    vv7 = yb(inl)*cssb(ims)
                                    vv8 = yb(ins)*cssb(iml)
                                    val = val -
     +                                 (vv1+vv2+vv3+vv4+vv5+vv6+vv7+vv8)
     +                                 /4.0d0
c
                                    if (dabs(val).gt.cut) then
                                      nint = nint + 1
                                      g(nint) = -val
                                      nint4 = nint + nint + nint + nint
                                      labout(nint4-3) = i1
                                      labout(nint4-2) = j1
                                      labout(nint4-1) = k1
                                      labout(nint4) = l1
                                      if (nint.eq.nintmx) then
                                      call pack(g(num2e+1),lab816,
     +                                          labout,numlab)
                                      call put(g,m511,ifw)
                                      nint = 0
                                    call izero(1360,labout,1)
                                       end if
                                    end if
c
 90                              continue
 100                          continue
                              icount = icount + 1
 110                       continue
 120                    continue
                     end if
 130              continue
 140           continue
            end if
 150     continue
 160  continue
      call pack(g(num2e+1),lab816,labout,numlab)
      call put(g,m511,ifw)
      call put(g,mzero,ifw)
      call clredx
      return
      end
      subroutine ump3sa(umm,amm,ammz,uds,ads,vec,vecz,dum,zint,zint1,
     1                  zint2,a1,a2,s1,s2,yint,yintz,e,ez,
     1      ncoorb,nocc,noccz,ndep,ndepz,mn,mnz,ifort1,ifort2,ifort3)
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension umm(ndep),amm(ndep),ammz(ndepz),uds(mn,mnz),ads(mn,mnz),
     1        vec(ncoorb,ncoorb),vecz(ncoorb,ncoorb),dum(ncoorb,ncoorb),
     1    zint(ncoorb*ncoorb),zint1(ncoorb*ncoorb),zint2(ncoorb*ncoorb),
     1      a1(ncoorb,ncoorb),
     1        a2(ncoorb,ncoorb),s1(ncoorb,ncoorb),s2(ncoorb,ncoorb),
     1   yint(mn),yintz(mnz),e(ncoorb),ez(ncoorb)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      nsq = ncoorb*ncoorb
      nocc1 = nocc + 1
      noccz1 = noccz + 1
      nvirt = ncoorb - nocc
      nvirtz = ncoorb - noccz
      ijm = (nocc+1)*nocc/2
      ijmz = (noccz+1)*noccz/2
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
c
      call vclr(a1,1,nsq)
      call vclr(yintz,1,mnz)
      call vclr(zint,1,nsq)
      call vclr(zint1,1,nsq)
      call vclr(zint2,1,nsq)
      do 650 ip = 1 , ncoorb
         do 640 iq = 1 , ip
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ib = ip
                  j = iq
                  ejb = e(ib) - e(j)
                  ijb = (ib-nocc1)*nocc + j
                  do 30 ia = nocc1 , ncoorb
                     mfact = 1
                     iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                     if (ib.gt.ia) then
                        iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                        mfact = -mfact
                     end if
                     do 20 i = 1 , nocc
                        mfact2 = mfact
                        ij = (i-1)*i/2 + j
                        if (j.gt.i) then
                           ij = (j-1)*j/2 + i
                           mfact2 = -mfact
                        end if
                        iabij = (iab-1)*ijm + ij
                        difen = ejb + e(ia) - e(i)
                        iai = (ia-nocc1)*nocc + i
c
                        yint(iai) = amm(iabij)*mfact2
                        a1(i,ia) = -umm(iabij)/difen*mfact2*4.0d0
 20                  continue
 30               continue
c
c
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(a1(1,nocc1),1,ncoorb,vec(1,nocc1),ncoorb,1,
     +                      dum,1,ncoorb,nocc,nvirt,ncoorb)
c
                  call vclr(zint,1,nsq)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
                  call vclr(dum,1,nsq)
c
c   ump2 2-pdm inside u3suba
c
                  call mxmb(yint,1,nocc,vec(1,nocc1),ncoorb,1,dum,1,
     +                      ncoorb,nocc,nvirt,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  do 50 ia = noccz1 , ncoorb
                     do 40 i = 1 , noccz
                        iai = (ia-noccz1)*noccz + i
                        difen = ejb + ez(ia) - ez(i)
                        yintz(iai) = -uds(ijb,iai)/difen*4.0d0
 40                  continue
 50               continue
c
                  call mxmb(yintz,1,noccz,vecz(1,noccz1),ncoorb,1,dum,1,
     +                      ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c
                  call vclr(dum,1,nsq)
c
c  ump2 2-pdm
c
                  call mxmb(ads(ijb,1),mn,mn*noccz,vecz(1,noccz1),
     +                      ncoorb,1,dum,1,ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c  next term
c
                  ia = ip
                  i = iq
                  iai = (ia-nocc1)*nocc + i
                  call vclr(yint,1,mn)
                  do 110 ic = nocc1 , ncoorb
                     do 100 k = 1 , nocc
                        ikc = (ic-nocc1)*nocc + k
                        do 70 ib = nocc1 , ncoorb
                           mfact = 1
                           mfact1 = 1
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact = -mfact
                           end if
                           icb = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                           if (ib.gt.ic) then
                              icb = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                              mfact1 = -mfact1
                           end if
                           do 60 j = 1 , nocc
                              ij = (i-1)*i/2 + j
                              mfact2 = mfact
                              if (j.gt.i) then
                                 ij = (j-1)*j/2 + i
                                 mfact2 = -mfact
                              end if
                              kj = (k-1)*k/2 + j
                              mfact3 = mfact1
                              if (j.gt.k) then
                                 kj = (j-1)*j/2 + k
                                 mfact3 = -mfact1
                              end if
                              iabij = (iab-1)*ijm + ij
                              icbkj = (icb-1)*ijm + kj
c
                              val1 = amm(iabij)*mfact2
                              val2 = amm(icbkj)*mfact3
                              yint(ikc) = yint(ikc) - val1*val2
 60                        continue
 70                     continue
c
                        do 90 ib = noccz1 , ncoorb
                           do 80 j = 1 , noccz
                              ijb = (ib-noccz1)*noccz + j
c
                              yint(ikc) = yint(ikc) - ads(iai,ijb)
     +                           *ads(ikc,ijb)
 80                        continue
 90                     continue
 100                 continue
 110              continue
c
                  call vclr(yintz,1,mnz)
                  do 170 ic = noccz1 , ncoorb
                     do 160 k = 1 , noccz
                        ikc = (ic-noccz1)*noccz + k
c
                        do 130 ib = nocc1 , ncoorb
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           mfact = 1
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact = -mfact
                           end if
                           do 120 j = 1 , nocc
                              ijb = (ib-nocc1)*nocc + j
                              ij = (i-1)*i/2 + j
                              mfact2 = mfact
                              if (j.gt.i) then
                                 ij = (j-1)*j/2 + i
                                 mfact2 = -mfact
                              end if
                              iabij = (iab-1)*ijm + ij
c
                              yintz(ikc) = yintz(ikc) - amm(iabij)
     +                           *ads(ijb,ikc)*mfact2
 120                       continue
 130                    continue
c
                        do 150 ib = noccz1 , ncoorb
                           mfact = 1
                           icb = (ic-noccz1)*(ic-noccz)/2 + ib - noccz
                           if (ib.gt.ic) then
                              icb = (ib-noccz1)*(ib-noccz)
     +                              /2 + ic - noccz
                              mfact = -mfact
                           end if
                           do 140 j = 1 , noccz
                              mfact2 = mfact
                              kj = (k-1)*k/2 + j
                              if (j.gt.k) then
                                 kj = (j-1)*j/2 + k
                                 mfact2 = -mfact
                              end if
                              icbkj = (icb-1)*ijmz + kj
                              ijb = (ib-noccz1)*noccz + j
c
                              yintz(ikc) = yintz(ikc) - ads(iai,ijb)
     +                           *ammz(icbkj)*mfact2
 140                       continue
 150                    continue
 160                 continue
 170              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(yint,1,nocc,vec(1,nocc1),ncoorb,1,dum,1,
     +                      ncoorb,nocc,nvirt,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(yintz,1,noccz,vecz(1,noccz1),ncoorb,1,dum,1,
     +                      ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
                  do 190 im = 1 , ncoorb
                     do 180 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        dum(im,in) = (zint(imn)+zint(inm))/2.0d0
 180                 continue
 190              continue
c
                  icount = 1
                  do 230 ii = 1 , nshell
                     do 220 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 210 i = mini , maxi
                           i1 = loci + i
                           do 200 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint(icount) = dum(i1,j1)
                              icount = icount + 1
 200                       continue
 210                    continue
 220                 continue
 230              continue
c
                  call wtedz(zint,nsq,ifort1)
               end if
c
               if (iq.gt.nocc) then
                  ib = ip
                  id = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 280 j = 1 , nocc
                     ijb = (ib-nocc1)*nocc + j
                     ijd = (id-nocc1)*nocc + j
                     do 250 ia = nocc1 , ncoorb
                        mfact = 1
                        mfact1 = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        iad = (ia-nocc1)*(ia-nocc)/2 + id - nocc
                        if (id.gt.ia) then
                           iad = (id-nocc1)*(id-nocc)/2 + ia - nocc
                           mfact1 = -mfact1
                        end if
                        do 240 i = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iadij = (iad-1)*ijm + ij
c
                           a1(i,ia) = amm(iabij)*mfact2
                           a2(i,ia) = -amm(iadij)*mfact3*0.25d0
 240                    continue
 250                 continue
c
                     do 270 ic = noccz1 , ncoorb
                        do 260 k = 1 , noccz
                           ikc = (ic-noccz1)*noccz + k
                           yintz(ikc) = -ads(ijb,ikc)*0.5d0
 260                    continue
 270                 continue
c
                     call mxmb(a1(1,nocc1),ncoorb,1,a2(1,nocc1),1,
     +                         ncoorb,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     call mxmb(yintz,noccz,1,ads(ijd,1),mn,mn*noccz,
     +                         s2(noccz1,noccz1),1,ncoorb,nvirtz,noccz,
     +                         nvirtz)
c
                     if (ib.ne.id) then
                        call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,
     +                            ncoorb,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                            nocc,nvirt)
c
                        call mxmb(ads(ijd,1),mn*noccz,mn,yintz,1,noccz,
     +                            s2(noccz1,noccz1),1,ncoorb,nvirtz,
     +                            noccz,nvirtz)
c
                     end if
 280              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1(nocc1,nocc1),1,ncoorb,vec(1,nocc1),
     +                      ncoorb,1,dum,1,ncoorb,nvirt,nvirt,ncoorb)
c
                  call vclr(zint1,1,nsq)
c
                  call mxmb(vec(1,nocc1),1,ncoorb,dum,1,ncoorb,zint1,1,
     +                      ncoorb,ncoorb,nvirt,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2(noccz1,noccz1),1,ncoorb,vecz(1,noccz1),
     +                      ncoorb,1,dum,1,ncoorb,nvirtz,nvirtz,ncoorb)
c
                  call mxmb(vecz(1,noccz1),1,ncoorb,dum,1,ncoorb,zint1,
     +                      1,ncoorb,ncoorb,nvirtz,ncoorb)
c
c  next term
c
                  ia = ip
                  ic = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 310 ib = nocc1 , ncoorb
                     mfact = 1
                     mfact1 = 1
                     iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                     if (ib.gt.ia) then
                        iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                        mfact = -mfact
                     end if
                     icb = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                     if (ib.gt.ic) then
                        icb = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                        mfact1 = -mfact1
                     end if
                     do 300 i = 1 , nocc
                        do 290 j = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           icbij = (icb-1)*ijm + ij
c
                           a1(i,j) = amm(iabij)*mfact2
                           a2(i,j) = amm(icbij)*mfact3*0.5d0
 290                    continue
 300                 continue
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,s1,1,ncoorb,nocc,
     +                         nocc,nocc)
c
                     if (ia.ne.ic) call mxmb(a2,1,ncoorb,a1,ncoorb,1,s1,
     +                   1,ncoorb,nocc,nocc,nocc)
 310              continue
c
                  do 340 ib = noccz1 , ncoorb
                     do 330 j = 1 , noccz
                        do 320 i = 1 , nocc
                           ijb = (ib-noccz1)*noccz + j
                           iai = (ia-nocc1)*nocc + i
                           ici = (ic-nocc1)*nocc + i
c
                           a1(i,j) = ads(iai,ijb)
                           a2(i,j) = ads(ici,ijb)*0.5d0
 320                    continue
 330                 continue
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,s1,1,ncoorb,nocc,
     +                         noccz,nocc)
c
                     if (ia.ne.ic) call mxmb(a2,1,ncoorb,a1,ncoorb,1,s1,
     +                   1,ncoorb,nocc,noccz,nocc)
 340              continue
c
                  do 370 j = 1 , nocc
                     ija = (ia-nocc1)*nocc + j
                     ijc = (ic-nocc1)*nocc + j
c
                     do 360 iia = noccz1 , ncoorb
                        do 350 iii = 1 , noccz
                           iiia = (iia-noccz1)*noccz + iii
                           yintz(iiia) = ads(ija,iiia)*0.5d0
 350                    continue
 360                 continue
c
                     call mxmb(yintz,1,noccz,ads(ijc,1),mn*noccz,mn,s2,
     +                         1,ncoorb,noccz,nvirtz,noccz)
c
                     if (ia.ne.ic)
     +                   call mxmb(ads(ijc,1),mn,mn*noccz,yintz,noccz,1,
     +                   s2,1,ncoorb,noccz,nvirtz,noccz)
 370              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1,1,ncoorb,vec,ncoorb,1,dum,1,ncoorb,nocc,
     +                      nocc,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint1,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2,1,ncoorb,vecz,ncoorb,1,dum,1,ncoorb,
     +                      noccz,noccz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint1,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
                  do 390 im = 1 , ncoorb
                     do 380 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        a1(im,in) = (zint1(imn)+zint1(inm))/2.0d0
 380                 continue
 390              continue
c
                  icount = 1
                  do 430 ii = 1 , nshell
                     do 420 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 410 i = mini , maxi
                           i1 = loci + i
                           do 400 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint1(icount) = a1(i1,j1)
                              icount = icount + 1
 400                       continue
 410                    continue
 420                 continue
 430              continue
                  call wtedz(zint1,nsq,ifort2)
               end if
            end if
c
            if (ip.le.nocc) then
               if (iq.le.nocc) then
                  l = ip
                  j = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 480 ib = nocc1 , ncoorb
                     ijb = (ib-nocc1)*nocc + j
                     ilb = (ib-nocc1)*nocc + l
                     do 450 ia = nocc1 , ncoorb
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        mfact = 1
                        mfact1 = 1
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 440 i = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                           end if
                           il = (i-1)*i/2 + l
                           if (l.gt.i) then
                              il = (l-1)*l/2 + i
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iabil = (iab-1)*ijm + il
c
                           a1(i,ia) = -amm(iabij)*mfact2*0.25d0
                           a2(i,ia) = amm(iabil)*mfact3
 440                    continue
 450                 continue
c
                     do 470 ic = noccz1 , ncoorb
                        do 460 k = 1 , noccz
                           ikc = (ic-noccz1)*noccz + k
                           yintz(ikc) = -ads(ijb,ikc)*0.5d0
 460                    continue
 470                 continue
c
                     call mxmb(a1(1,nocc1),1,ncoorb,a2(1,nocc1),ncoorb,
     +                         1,s1,1,ncoorb,nocc,nvirt,nocc)
c
                     call mxmb(yintz,1,noccz,ads(ilb,1),mn*noccz,mn,s2,
     +                         1,ncoorb,noccz,nvirtz,noccz)
c
                     if (l.ne.j) then
                        call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),
     +                            ncoorb,1,s1,1,ncoorb,nocc,nvirt,nocc)
c
                        call mxmb(ads(ilb,1),mn,mn*noccz,yintz,noccz,1,
     +                            s2,1,ncoorb,noccz,nvirtz,noccz)
                     end if
 480              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1,1,ncoorb,vec,ncoorb,1,dum,1,ncoorb,nocc,
     +                      nocc,ncoorb)
c
                  call vclr(zint2,1,nsq)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint2,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2,1,ncoorb,vecz,ncoorb,1,dum,1,ncoorb,
     +                      noccz,noccz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint2,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c   part of (ac/ki) sacki term
c
                  i = ip
                  k = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
                  do 510 j = 1 , nocc
                     mfact = 1
                     mfact1 = 1
                     ij = (i-1)*i/2 + j
                     if (j.gt.i) then
                        ij = (j-1)*j/2 + i
                        mfact = -mfact
                     end if
                     kj = (k-1)*k/2 + j
                     if (j.gt.k) then
                        kj = (j-1)*j/2 + k
                        mfact1 = -mfact1
                     end if
c
                     do 500 ia = nocc1 , ncoorb
                        do 490 ib = nocc1 , ncoorb
                           mfact2 = mfact
                           mfact3 = mfact1
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iabkj = (iab-1)*ijm + kj
c
                           a1(ia,ib) = amm(iabij)*mfact2
                           a2(ia,ib) = amm(iabkj)*mfact3*0.5d0
 490                    continue
 500                 continue
c
                     call mxmb(a1(nocc1,nocc1),1,ncoorb,a2(nocc1,nocc1),
     +                         ncoorb,1,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nvirt,nvirt)
c
                     if (i.ne.k) call mxmb(a2(nocc1,nocc1),1,ncoorb,
     +                   a1(nocc1,nocc1),ncoorb,1,s1(nocc1,nocc1),1,
     +                   ncoorb,nvirt,nvirt,nvirt)
c
 510              continue
c
                  do 540 j = 1 , noccz
                     do 530 ia = nocc1 , ncoorb
                        do 520 ib = noccz1 , ncoorb
                           ijb = (ib-noccz1)*noccz + j
                           iai = (ia-nocc1)*nocc + i
                           iak = (ia-nocc1)*nocc + k
c
                           a1(ia,ib) = ads(iai,ijb)
                           a2(ia,ib) = ads(iak,ijb)*0.5d0
 520                    continue
 530                 continue
c
                     call mxmb(a1(nocc1,noccz1),1,ncoorb,
     +                         a2(nocc1,noccz1),ncoorb,1,s1(nocc1,nocc1)
     +                         ,1,ncoorb,nvirt,nvirtz,nvirt)
c
                     if (i.ne.k) call mxmb(a2(nocc1,noccz1),1,ncoorb,
     +                   a1(nocc1,noccz1),ncoorb,1,s1(nocc1,nocc1),1,
     +                   ncoorb,nvirt,nvirtz,nvirt)
c
 540              continue
c
                  do 570 ib = nocc1 , ncoorb
                     ibi = (ib-nocc1)*nocc + i
                     ibk = (ib-nocc1)*nocc + k
c
                     do 560 iia = noccz1 , ncoorb
                        do 550 iii = 1 , noccz
                           iiia = (iia-noccz1)*noccz + iii
                           yintz(iiia) = ads(ibi,iiia)*0.5d0
 550                    continue
 560                 continue
c
                     call mxmb(yintz,noccz,1,ads(ibk,1),mn,mn*noccz,
     +                         s2(noccz1,noccz1),1,ncoorb,nvirtz,noccz,
     +                         nvirtz)
c
                     if (i.ne.k) call mxmb(ads(ibk,1),mn*noccz,mn,yintz,
     +                   1,noccz,s2(noccz1,noccz1),1,ncoorb,nvirtz,
     +                   noccz,nvirtz)
c
 570              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1(nocc1,nocc1),1,ncoorb,vec(1,nocc1),
     +                      ncoorb,1,dum,1,ncoorb,nvirt,nvirt,ncoorb)
c
                  call mxmb(vec(1,nocc1),1,ncoorb,dum,1,ncoorb,zint2,1,
     +                      ncoorb,ncoorb,nvirt,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2(noccz1,noccz1),1,ncoorb,vecz(1,noccz1),
     +                      ncoorb,1,dum,1,ncoorb,nvirtz,nvirtz,ncoorb)
c
                  call mxmb(vecz(1,noccz1),1,ncoorb,dum,1,ncoorb,zint2,
     +                      1,ncoorb,ncoorb,nvirtz,ncoorb)
c
                  do 590 im = 1 , ncoorb
                     do 580 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        a2(im,in) = (zint2(imn)+zint2(inm))/2.0d0
 580                 continue
 590              continue
c
                  icount = 1
                  do 630 ii = 1 , nshell
                     do 620 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 610 i = mini , maxi
                           i1 = loci + i
                           do 600 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint2(icount) = a2(i1,j1)
                              icount = icount + 1
 600                       continue
 610                    continue
 620                 continue
 630              continue
                  call wtedz(zint2,nsq,ifort3)
               end if
            end if
c
 640     continue
 650  continue
      return
      end
      subroutine ump3sb(umm,amm,ammz,uds,ads,vec,vecz,dum,zint,zint1,
     1                  zint2,a1,a2,s1,s2,yint,yintz,e,ez,
     1      ncoorb,nocc,noccz,ndep,ndepz,mn,mnz,ifort1,ifort2,ifort3)
c
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension umm(ndep),amm(ndep),ammz(ndepz),uds(mnz,mn),ads(mnz,mn),
     1        vec(ncoorb,ncoorb),vecz(ncoorb,ncoorb),dum(ncoorb,ncoorb),
     1    zint(ncoorb*ncoorb),zint1(ncoorb*ncoorb),zint2(ncoorb*ncoorb),
     1      a1(ncoorb,ncoorb),
     1        a2(ncoorb,ncoorb),s1(ncoorb,ncoorb),s2(ncoorb,ncoorb),
     1   yint(mn),yintz(mnz),e(ncoorb),ez(ncoorb)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      nsq = ncoorb*ncoorb
      nocc1 = nocc + 1
      noccz1 = noccz + 1
      nvirt = ncoorb - nocc
      nvirtz = ncoorb - noccz
      ijm = (nocc+1)*nocc/2
      ijmz = (noccz+1)*noccz/2
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
c
      call vclr(a1,1,nsq)
      call vclr(yintz,1,mnz)
      call vclr(zint,1,nsq)
      call vclr(zint1,1,nsq)
      call vclr(zint2,1,nsq)
      do 650 ip = 1 , ncoorb
         do 640 iq = 1 , ip
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ib = ip
                  j = iq
                  ejb = e(ib) - e(j)
                  ijb = (ib-nocc1)*nocc + j
                  do 30 ia = nocc1 , ncoorb
                     mfact = 1
                     iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                     if (ib.gt.ia) then
                        iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                        mfact = -mfact
                     end if
                     do 20 i = 1 , nocc
                        mfact2 = mfact
                        ij = (i-1)*i/2 + j
                        if (j.gt.i) then
                           ij = (j-1)*j/2 + i
                           mfact2 = -mfact
                        end if
                        iabij = (iab-1)*ijm + ij
                        difen = ejb + e(ia) - e(i)
                        iai = (ia-nocc1)*nocc + i
c
                        yint(iai) = amm(iabij)*mfact2
                        a1(i,ia) = -umm(iabij)/difen*mfact2*4.0d0
 20                  continue
 30               continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(a1(1,nocc1),1,ncoorb,vec(1,nocc1),ncoorb,1,
     +                      dum,1,ncoorb,nocc,nvirt,ncoorb)
c
                  call vclr(zint,1,nsq)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
c ump2 2-pdm inside u3subb
c
                  call mxmb(yint,1,nocc,vec(1,nocc1),ncoorb,1,dum,1,
     +                      ncoorb,nocc,nvirt,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  do 50 ia = noccz1 , ncoorb
                     do 40 i = 1 , noccz
                        iai = (ia-noccz1)*noccz + i
                        difen = ejb + ez(ia) - ez(i)
                        yintz(iai) = -uds(iai,ijb)/difen*4.0d0
 40                  continue
 50               continue
c
                  call mxmb(yintz,1,noccz,vecz(1,noccz1),ncoorb,1,dum,1,
     +                      ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c
                  call vclr(dum,1,nsq)
c
c   ump2 2-pdm
c
                  call mxmb(ads(1,ijb),1,noccz,vecz(1,noccz1),ncoorb,1,
     +                      dum,1,ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c  next term
c
                  ia = ip
                  i = iq
                  iai = (ia-nocc1)*nocc + i
                  call vclr(yint,1,mn)
                  do 110 ic = nocc1 , ncoorb
                     do 100 k = 1 , nocc
                        ikc = (ic-nocc1)*nocc + k
                        do 70 ib = nocc1 , ncoorb
                           mfact = 1
                           mfact1 = 1
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact = -mfact
                           end if
                           icb = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                           if (ib.gt.ic) then
                              icb = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                              mfact1 = -mfact1
                           end if
                           do 60 j = 1 , nocc
                              ij = (i-1)*i/2 + j
                              mfact2 = mfact
                              if (j.gt.i) then
                                 ij = (j-1)*j/2 + i
                                 mfact2 = -mfact
                              end if
                              kj = (k-1)*k/2 + j
                              mfact3 = mfact1
                              if (j.gt.k) then
                                 kj = (j-1)*j/2 + k
                                 mfact3 = -mfact1
                              end if
                              iabij = (iab-1)*ijm + ij
                              icbkj = (icb-1)*ijm + kj
c
                              val1 = amm(iabij)*mfact2
                              val2 = amm(icbkj)*mfact3
                              yint(ikc) = yint(ikc) - val1*val2
 60                        continue
 70                     continue
c
                        do 90 ib = noccz1 , ncoorb
                           do 80 j = 1 , noccz
                              ijb = (ib-noccz1)*noccz + j
c
                              yint(ikc) = yint(ikc) - ads(ijb,iai)
     +                           *ads(ijb,ikc)
 80                        continue
 90                     continue
 100                 continue
 110              continue
c
                  call vclr(yintz,1,mnz)
                  do 170 ic = noccz1 , ncoorb
                     do 160 k = 1 , noccz
                        ikc = (ic-noccz1)*noccz + k
c
                        do 130 ib = nocc1 , ncoorb
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           mfact = 1
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact = -mfact
                           end if
                           do 120 j = 1 , nocc
                              ijb = (ib-nocc1)*nocc + j
                              ij = (i-1)*i/2 + j
                              mfact2 = mfact
                              if (j.gt.i) then
                                 ij = (j-1)*j/2 + i
                                 mfact2 = -mfact
                              end if
                              iabij = (iab-1)*ijm + ij
c
                              yintz(ikc) = yintz(ikc) - amm(iabij)
     +                           *ads(ikc,ijb)*mfact2
 120                       continue
 130                    continue
c
                        do 150 ib = noccz1 , ncoorb
                           mfact = 1
                           icb = (ic-noccz1)*(ic-noccz)/2 + ib - noccz
                           if (ib.gt.ic) then
                              icb = (ib-noccz1)*(ib-noccz)
     +                              /2 + ic - noccz
                              mfact = -mfact
                           end if
                           do 140 j = 1 , noccz
                              mfact2 = mfact
                              kj = (k-1)*k/2 + j
                              if (j.gt.k) then
                                 kj = (j-1)*j/2 + k
                                 mfact2 = -mfact
                              end if
                              icbkj = (icb-1)*ijmz + kj
                              ijb = (ib-noccz1)*noccz + j
c
                              yintz(ikc) = yintz(ikc) - ads(ijb,iai)
     +                           *ammz(icbkj)*mfact2
 140                       continue
 150                    continue
 160                 continue
 170              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(yint,1,nocc,vec(1,nocc1),ncoorb,1,dum,1,
     +                      ncoorb,nocc,nvirt,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(yintz,1,noccz,vecz(1,noccz1),ncoorb,1,dum,1,
     +                      ncoorb,noccz,nvirtz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
                  do 190 im = 1 , ncoorb
                     do 180 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        dum(im,in) = (zint(imn)+zint(inm))/2.0d0
 180                 continue
 190              continue
c
                  icount = 1
                  do 230 ii = 1 , nshell
                     do 220 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 210 i = mini , maxi
                           i1 = loci + i
                           do 200 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint(icount) = dum(i1,j1)
                              icount = icount + 1
 200                       continue
 210                    continue
 220                 continue
 230              continue
c
                  call wtedz(zint,nsq,ifort1)
               end if
c
               if (iq.gt.nocc) then
                  ib = ip
                  id = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 280 j = 1 , nocc
                     ijb = (ib-nocc1)*nocc + j
                     ijd = (id-nocc1)*nocc + j
                     do 250 ia = nocc1 , ncoorb
                        mfact = 1
                        mfact1 = 1
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                        end if
                        iad = (ia-nocc1)*(ia-nocc)/2 + id - nocc
                        if (id.gt.ia) then
                           iad = (id-nocc1)*(id-nocc)/2 + ia - nocc
                           mfact1 = -mfact1
                        end if
                        do 240 i = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iadij = (iad-1)*ijm + ij
c
                           a1(i,ia) = amm(iabij)*mfact2
                           a2(i,ia) = -amm(iadij)*mfact3*0.25d0
 240                    continue
 250                 continue
c
                     do 270 ic = noccz1 , ncoorb
                        do 260 k = 1 , noccz
                           ikc = (ic-noccz1)*noccz + k
                           yintz(ikc) = -ads(ikc,ijb)*0.5d0
 260                    continue
 270                 continue
c
                     call mxmb(a1(1,nocc1),ncoorb,1,a2(1,nocc1),1,
     +                         ncoorb,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nocc,nvirt)
c
                     call mxmb(yintz,noccz,1,ads(1,ijd),1,noccz,
     +                         s2(noccz1,noccz1),1,ncoorb,nvirtz,noccz,
     +                         nvirtz)
c
                     if (ib.ne.id) then
                        call mxmb(a2(1,nocc1),ncoorb,1,a1(1,nocc1),1,
     +                            ncoorb,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                            nocc,nvirt)
c
                        call mxmb(ads(1,ijd),noccz,1,yintz,1,noccz,
     +                            s2(noccz1,noccz1),1,ncoorb,nvirtz,
     +                            noccz,nvirtz)
c
                     end if
 280              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1(nocc1,nocc1),1,ncoorb,vec(1,nocc1),
     +                      ncoorb,1,dum,1,ncoorb,nvirt,nvirt,ncoorb)
c
                  call vclr(zint1,1,nsq)
c
                  call mxmb(vec(1,nocc1),1,ncoorb,dum,1,ncoorb,zint1,1,
     +                      ncoorb,ncoorb,nvirt,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2(noccz1,noccz1),1,ncoorb,vecz(1,noccz1),
     +                      ncoorb,1,dum,1,ncoorb,nvirtz,nvirtz,ncoorb)
c
                  call mxmb(vecz(1,noccz1),1,ncoorb,dum,1,ncoorb,zint1,
     +                      1,ncoorb,ncoorb,nvirtz,ncoorb)
c
c  next term
c
                  ia = ip
                  ic = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 310 ib = nocc1 , ncoorb
                     mfact = 1
                     mfact1 = 1
                     iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                     if (ib.gt.ia) then
                        iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                        mfact = -mfact
                     end if
                     icb = (ic-nocc1)*(ic-nocc)/2 + ib - nocc
                     if (ib.gt.ic) then
                        icb = (ib-nocc1)*(ib-nocc)/2 + ic - nocc
                        mfact1 = -mfact1
                     end if
                     do 300 i = 1 , nocc
                        do 290 j = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           icbij = (icb-1)*ijm + ij
c
                           a1(i,j) = amm(iabij)*mfact2
                           a2(i,j) = amm(icbij)*mfact3*0.5d0
 290                    continue
 300                 continue
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,s1,1,ncoorb,nocc,
     +                         nocc,nocc)
c
                     if (ia.ne.ic) call mxmb(a2,1,ncoorb,a1,ncoorb,1,s1,
     +                   1,ncoorb,nocc,nocc,nocc)
 310              continue
c
                  do 340 ib = noccz1 , ncoorb
                     do 330 j = 1 , noccz
                        do 320 i = 1 , nocc
                           ijb = (ib-noccz1)*noccz + j
                           iai = (ia-nocc1)*nocc + i
                           ici = (ic-nocc1)*nocc + i
c
                           a1(i,j) = ads(ijb,iai)
                           a2(i,j) = ads(ijb,ici)*0.5d0
 320                    continue
 330                 continue
c
                     call mxmb(a1,1,ncoorb,a2,ncoorb,1,s1,1,ncoorb,nocc,
     +                         noccz,nocc)
c
                     if (ia.ne.ic) call mxmb(a2,1,ncoorb,a1,ncoorb,1,s1,
     +                   1,ncoorb,nocc,noccz,nocc)
 340              continue
c
                  do 370 j = 1 , nocc
                     ija = (ia-nocc1)*nocc + j
                     ijc = (ic-nocc1)*nocc + j
c
                     do 360 iia = noccz1 , ncoorb
                        do 350 iii = 1 , noccz
                           iiia = (iia-noccz1)*noccz + iii
                           yintz(iiia) = ads(iiia,ija)*0.5d0
 350                    continue
 360                 continue
c
                     call mxmb(yintz,1,noccz,ads(1,ijc),noccz,1,s2,1,
     +                         ncoorb,noccz,nvirtz,noccz)
c
                     if (ia.ne.ic) call mxmb(ads(1,ijc),1,noccz,yintz,
     +                   noccz,1,s2,1,ncoorb,noccz,nvirtz,noccz)
 370              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1,1,ncoorb,vec,ncoorb,1,dum,1,ncoorb,nocc,
     +                      nocc,ncoorb)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint1,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2,1,ncoorb,vecz,ncoorb,1,dum,1,ncoorb,
     +                      noccz,noccz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint1,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
                  do 390 im = 1 , ncoorb
                     do 380 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        a1(im,in) = (zint1(imn)+zint1(inm))/2.0d0
 380                 continue
 390              continue
c
                  icount = 1
                  do 430 ii = 1 , nshell
                     do 420 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 410 i = mini , maxi
                           i1 = loci + i
                           do 400 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint1(icount) = a1(i1,j1)
                              icount = icount + 1
 400                       continue
 410                    continue
 420                 continue
 430              continue
                  call wtedz(zint1,nsq,ifort2)
               end if
            end if
c
            if (ip.le.nocc) then
               if (iq.le.nocc) then
                  l = ip
                  j = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
c
                  do 480 ib = nocc1 , ncoorb
                     ijb = (ib-nocc1)*nocc + j
                     ilb = (ib-nocc1)*nocc + l
                     do 450 ia = nocc1 , ncoorb
                        iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                        mfact = 1
                        mfact1 = 1
                        if (ib.gt.ia) then
                           iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                           mfact = -mfact
                           mfact1 = -mfact1
                        end if
                        do 440 i = 1 , nocc
                           mfact2 = mfact
                           mfact3 = mfact1
                           ij = (i-1)*i/2 + j
                           if (j.gt.i) then
                              ij = (j-1)*j/2 + i
                              mfact2 = -mfact
                           end if
                           il = (i-1)*i/2 + l
                           if (l.gt.i) then
                              il = (l-1)*l/2 + i
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iabil = (iab-1)*ijm + il
c
                           a1(i,ia) = -amm(iabij)*mfact2*0.25d0
                           a2(i,ia) = amm(iabil)*mfact3
 440                    continue
 450                 continue
c
                     do 470 ic = noccz1 , ncoorb
                        do 460 k = 1 , noccz
                           ikc = (ic-noccz1)*noccz + k
                           yintz(ikc) = -ads(ikc,ijb)*0.5d0
 460                    continue
 470                 continue
c
                     call mxmb(a1(1,nocc1),1,ncoorb,a2(1,nocc1),ncoorb,
     +                         1,s1,1,ncoorb,nocc,nvirt,nocc)
c
                     call mxmb(yintz,1,noccz,ads(1,ilb),noccz,1,s2,1,
     +                         ncoorb,noccz,nvirtz,noccz)
c
                     if (l.ne.j) then
                        call mxmb(a2(1,nocc1),1,ncoorb,a1(1,nocc1),
     +                            ncoorb,1,s1,1,ncoorb,nocc,nvirt,nocc)
c
                        call mxmb(ads(1,ilb),1,noccz,yintz,noccz,1,s2,1,
     +                            ncoorb,noccz,nvirtz,noccz)
                     end if
 480              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1,1,ncoorb,vec,ncoorb,1,dum,1,ncoorb,nocc,
     +                      nocc,ncoorb)
c
                  call vclr(zint2,1,nsq)
c
                  call mxmb(vec,1,ncoorb,dum,1,ncoorb,zint2,1,ncoorb,
     +                      ncoorb,nocc,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2,1,ncoorb,vecz,ncoorb,1,dum,1,ncoorb,
     +                      noccz,noccz,ncoorb)
c
                  call mxmb(vecz,1,ncoorb,dum,1,ncoorb,zint2,1,ncoorb,
     +                      ncoorb,noccz,ncoorb)
c
c   part of (ac/ki) sacki term
c
                  i = ip
                  k = iq
                  call vclr(s1,1,nsq)
                  call vclr(s2,1,nsq)
                  do 510 j = 1 , nocc
                     mfact = 1
                     mfact1 = 1
                     ij = (i-1)*i/2 + j
                     if (j.gt.i) then
                        ij = (j-1)*j/2 + i
                        mfact = -mfact
                     end if
                     kj = (k-1)*k/2 + j
                     if (j.gt.k) then
                        kj = (j-1)*j/2 + k
                        mfact1 = -mfact1
                     end if
c
                     do 500 ia = nocc1 , ncoorb
                        do 490 ib = nocc1 , ncoorb
                           mfact2 = mfact
                           mfact3 = mfact1
                           iab = (ia-nocc1)*(ia-nocc)/2 + ib - nocc
                           if (ib.gt.ia) then
                              iab = (ib-nocc1)*(ib-nocc)/2 + ia - nocc
                              mfact2 = -mfact
                              mfact3 = -mfact1
                           end if
                           iabij = (iab-1)*ijm + ij
                           iabkj = (iab-1)*ijm + kj
c
                           a1(ia,ib) = amm(iabij)*mfact2
                           a2(ia,ib) = amm(iabkj)*mfact3*0.5d0
 490                    continue
 500                 continue
c
                     call mxmb(a1(nocc1,nocc1),1,ncoorb,a2(nocc1,nocc1),
     +                         ncoorb,1,s1(nocc1,nocc1),1,ncoorb,nvirt,
     +                         nvirt,nvirt)
c
                     if (i.ne.k) call mxmb(a2(nocc1,nocc1),1,ncoorb,
     +                   a1(nocc1,nocc1),ncoorb,1,s1(nocc1,nocc1),1,
     +                   ncoorb,nvirt,nvirt,nvirt)
c
 510              continue
c
                  do 540 j = 1 , noccz
                     do 530 ia = nocc1 , ncoorb
                        do 520 ib = noccz1 , ncoorb
                           ijb = (ib-noccz1)*noccz + j
                           iai = (ia-nocc1)*nocc + i
                           iak = (ia-nocc1)*nocc + k
c
                           a1(ia,ib) = ads(ijb,iai)
                           a2(ia,ib) = ads(ijb,iak)*0.5d0
 520                    continue
 530                 continue
c
                     call mxmb(a1(nocc1,noccz1),1,ncoorb,
     +                         a2(nocc1,noccz1),ncoorb,1,s1(nocc1,nocc1)
     +                         ,1,ncoorb,nvirt,nvirtz,nvirt)
c
                     if (i.ne.k) call mxmb(a2(nocc1,noccz1),1,ncoorb,
     +                   a1(nocc1,noccz1),ncoorb,1,s1(nocc1,nocc1),1,
     +                   ncoorb,nvirt,nvirtz,nvirt)
c
 540              continue
c
                  do 570 ib = nocc1 , ncoorb
                     ibi = (ib-nocc1)*nocc + i
                     ibk = (ib-nocc1)*nocc + k
c
                     do 560 iia = noccz1 , ncoorb
                        do 550 iii = 1 , noccz
                           iiia = (iia-noccz1)*noccz + iii
                           yintz(iiia) = ads(iiia,ibi)*0.5d0
 550                    continue
 560                 continue
c
c
                     call mxmb(yintz,noccz,1,ads(1,ibk),1,noccz,
     +                         s2(noccz1,noccz1),1,ncoorb,nvirtz,noccz,
     +                         nvirtz)
c
                     if (i.ne.k) call mxmb(ads(1,ibk),noccz,1,yintz,1,
     +                   noccz,s2(noccz1,noccz1),1,ncoorb,nvirtz,noccz,
     +                   nvirtz)
c
 570              continue
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s1(nocc1,nocc1),1,ncoorb,vec(1,nocc1),
     +                      ncoorb,1,dum,1,ncoorb,nvirt,nvirt,ncoorb)
c
                  call mxmb(vec(1,nocc1),1,ncoorb,dum,1,ncoorb,zint2,1,
     +                      ncoorb,ncoorb,nvirt,ncoorb)
c
                  call vclr(dum,1,nsq)
c
                  call mxmb(s2(noccz1,noccz1),1,ncoorb,vecz(1,noccz1),
     +                      ncoorb,1,dum,1,ncoorb,nvirtz,nvirtz,ncoorb)
c
                  call mxmb(vecz(1,noccz1),1,ncoorb,dum,1,ncoorb,zint2,
     +                      1,ncoorb,ncoorb,nvirtz,ncoorb)
c
                  do 590 im = 1 , ncoorb
                     do 580 in = 1 , ncoorb
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        a2(im,in) = (zint2(imn)+zint2(inm))/2.0d0
 580                 continue
 590              continue
c
                  icount = 1
                  do 630 ii = 1 , nshell
                     do 620 jj = 1 , ii
                        mini = kmin(ii)
                        minj = kmin(jj)
                        maxi = kmax(ii)
                        maxj = kmax(jj)
                        loci = kloc(ii) - mini
                        locj = kloc(jj) - minj
c
                        do 610 i = mini , maxi
                           i1 = loci + i
                           do 600 j = minj , maxj
                              j1 = locj + j
c
c                             ij1 = (j1-1)*ncoorb + i1
                              zint2(icount) = a2(i1,j1)
                              icount = icount + 1
 600                       continue
 610                    continue
 620                 continue
 630              continue
                  call wtedz(zint2,nsq,ifort3)
               end if
            end if
c
 640     continue
 650  continue
      return
      end
      subroutine ump3dm(q,iblw,ifw,ifytr,ifdma,ifdmb,iblkya,iblkyb,
     1            ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,iw)
c
      implicit real*8  (a-h,o-z)
      logical dpres,gpres,uhf
      dimension q(*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      character *8 open
      data open/'open'/
      uhf = scftyp.eq.open
c
      dpres = .false.
      fpres = .false.
      gpres = .false.
      do 20 i = 1 , nshell
         if (ktype(i).eq.3) dpres = .true.
         if (ktype(i).eq.4) fpres = .true.
         if (ktype(i).eq.5) gpres = .true.
         if (ktype(i).gt.5) call caserr('ump3dm: upto g-functions only')
 20   continue
      nsq = ncoorb*ncoorb
      nova = nocca*nvirta
      novb = noccb*nvirtb
      ntri = ncoorb*(ncoorb+1)/2
      mna = nocca*nvirta
      mnb = noccb*nvirtb
      mn = mna + mnb
      nija = nocca*(nocca+1)/2
      nijb = noccb*(noccb+1)/2
      naba = nvirta*(nvirta+1)/2
      nabb = nvirtb*(nvirtb+1)/2
      ntri = ncoorb*(ncoorb+1)/2
      nsq = ncoorb*ncoorb
c     ifile1 = 1
c     ndepa = naba*nija
c     ndepb = nabb*nijb
      nabas = nvirta*nvirta
      nabbs = nvirtb*nvirtb
      nijas = nocca*nocca
      nijbs = noccb*noccb
c
c
c
      ibase = igmem_alloc_all(maxq)
      call mpsrt0(ifort1,ifort1,nova,nsq,q(ibase),q(ibase),maxq)
      call mpsrt0(ifort2,ifort2,naba,nsq,q(ibase),q(ibase),maxq)
      call mpsrt0(ifort3,ifort3,nija,nsq,q(ibase),q(ibase),maxq)
      if (uhf) then
         call mpsrt0(ifort4,ifort4,novb,nsq,q(ibase),q(ibase),maxq)
         call mpsrt0(ifort5,ifort5,nabb,nsq,q(ibase),q(ibase),maxq)
         call mpsrt0(ifort6,ifort6,nijb,nsq,q(ibase),q(ibase),maxq)
      end if
      call gmem_free(ibase)
c
c
c     i1=1+nsq
      iso = 1
      i0 = iso + nw196(5)
      i1 = i0 + nsq
      i2 = i1 + nsq
      i3 = i2 + nsq
      i4 = i3 + nsq
      i5 = i4 + nsq
      i6 = i5 + nsq
      i7 = i6 + nsq
      i8 = i7 + nsq
      i9 = i8 + nsq
      i10 = i9 + nova
      i11 = i10 + nabas
      i12 = i11 + nijas
      i13 = i12 + novb
      i14 = i13 + nabbs
      i15 = i14 + nijbs
      itop = i15 + ntri*16
      if (dpres) itop = i15 + ntri*36
      if (fpres) itop = i15 + ntri*100
      if (gpres) itop = i15 + ntri*225
      maxq = igmem_max_memory()
      if (itop.gt.maxq) then
         write (iw,6010) maxq , itop
         call caserr(' not enough core ')
      end if
      iso = igmem_alloc(nw196(5))
      i0  = igmem_alloc(nsq)
      i1  = igmem_alloc(nsq)
      i2  = igmem_alloc(nsq)
      i3  = igmem_alloc(nsq)
      i4  = igmem_alloc(nsq)
      i5  = igmem_alloc(nsq)
      i6  = igmem_alloc(nsq)
      i7  = igmem_alloc(nsq)
      i8  = igmem_alloc(nsq)
      i9  = igmem_alloc(nova)
      i10 = igmem_alloc(nabas)
      i11 = igmem_alloc(nijas)
      i12 = igmem_alloc(novb)
      i13 = igmem_alloc(nabbs)
      i14 = igmem_alloc(nijbs)
      if (gpres) then
         i15 = igmem_alloc(ntri*225)
      else if (fpres) then
         i15 = igmem_alloc(ntri*100)
      elseif (dpres) then
         i15 = igmem_alloc(ntri*36)
      else
         i15 = igmem_alloc(ntri*16)
      endif
c
c     m8 = 8
c     m11 = 11
      itype = 0
      call secget(isect(8),itype,iblk8)
      iblk8 = iblk8 + mvadd
      call rdedx(q(i6),nsq,iblk8,ifild)
      if (.not.uhf) then
         call dcopy(nsq,q(i6),1,q(i7),1)
      else
         itype = 0
         call secget(isect(11),itype,iblk11)
         iblk11 = iblk11 + mvadd
         call rdedx(q(i7),nsq,iblk11,ifild)
      end if
c
c
      call umps3(q(i6),q(i7),q(i8),q(i5),ifytr,ifdma,ifdmb,nocca,noccb,
     +  ncoorb,iblkya,iblkyb,uhf)
c
c     m8 = 8
c     m11 = 11
      itype = 0
      call secget(isect(8),itype,iblk8)
      iblk8 = iblk8 + mvadd
      call rdedx(q(i6),nsq,iblk8,ifild)
      if (.not.uhf) then
         call dcopy(nsq,q(i6),1,q(i7),1)
      else
         itype = 0
         call secget(isect(11),itype,iblk11)
         iblk11 = iblk11 + mvadd
         call rdedx(q(i7),nsq,iblk11,ifild)
      end if
c
      call rdedx(q(iso),nw196(5),ibl196(5),ifild)
      call umps5(q(iso),nshell,
     +  q(i0),q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),
     +  q(i8),q(i9),q(i10),q(i11),q(i12),q(i13),q(i14),q(i15),ifort1,
     +  ifort2,ifort3,ifort4,ifort5,ifort6,nocca,noccb,ncoorb,nova,
     +  novb,nijas,nijbs,nabas,nabbs,ifytr,ifdma,ifdmb,iblw,ifw,iblkya,
     +  iblkyb,uhf)
c
      call gmem_free(i15)
      call gmem_free(i14)
      call gmem_free(i13)
      call gmem_free(i12)
      call gmem_free(i11)
      call gmem_free(i10)
      call gmem_free(i9)
      call gmem_free(i8)
      call gmem_free(i7)
      call gmem_free(i6)
      call gmem_free(i5)
      call gmem_free(i4)
      call gmem_free(i3)
      call gmem_free(i2)
      call gmem_free(i1)
      call gmem_free(i0)
      call gmem_free(iso)
      return
 6010 format (/' insufficient core for umps5 ',i9,' real words ',
     +        ' need ',i9,' real words ')
      end
      subroutine ump3pd(q,umma,amma,ummb,ammb,uds,ads,
     +           ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,iw)
c
      implicit real*8  (a-h,o-z)
      dimension q(*)
      dimension umma(*),amma(*),ummb(*),ammb(*),uds(*),ads(*)
c
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
      logical uhf
      character *8 open
      data open/'open'/
      uhf = scftyp.eq.open
c
c     nova = nocca*nvirta
c     novb = noccb*nvirtb
c     ntri = ncoorb*(ncoorb+1)/2
c     ifile1 = 1
      nsq = ncoorb*ncoorb
      mna = nocca*nvirta
      mnb = noccb*nvirtb
      mn = mna + mnb
      nija = nocca*(nocca+1)/2
      nijb = noccb*(noccb+1)/2
      naba = nvirta*(nvirta+1)/2
      nabb = nvirtb*(nvirtb+1)/2
      nsq = ncoorb*ncoorb
      ndepa = naba*nija
      ndepb = nabb*nijb
c
      i6 = 1
      i7 = i6 + nsq
      i8 = i7 + nsq
      i9 = i8 + nsq
      i10 = i9 + nsq
      i11 = i10 + nsq
      i12 = i11 + nsq
      i13 = i12 + nsq
      i14 = i13 + nsq
      i15 = i14 + nsq
      i16 = i15 + nsq
      i17 = i16 + mna
      i18 = i17 + mnb
      i19 = i18 + ncoorb
      itop = i19 + ncoorb
      maxq = igmem_max_memory()
      if (itop.gt.maxq) then
         write (iw,6010) maxq , itop
         call caserr(' not enough core ')
      end if
c
      i6  = igmem_alloc(nsq)
      i7  = igmem_alloc(nsq)
      i8  = igmem_alloc(nsq)
      i9  = igmem_alloc(nsq)
      i10 = igmem_alloc(nsq)
      i11 = igmem_alloc(nsq)
      i12 = igmem_alloc(nsq)
      i13 = igmem_alloc(nsq)
      i14 = igmem_alloc(nsq)
      i15 = igmem_alloc(nsq)
      i16 = igmem_alloc(mna)
      i17 = igmem_alloc(mnb)
      i18 = igmem_alloc(ncoorb)
      i19 = igmem_alloc(ncoorb)
c
      itype = 0
      call secget(isect(8),itype,iblk8)
      iblk8 = iblk8 + mvadd
      call rdedx(q(i6),nsq,iblk8,ifild)
      if (.not.uhf) then
         call dcopy(nsq,q(i6),1,q(i7),1)
      else
         itype = 0
         call secget(isect(11),itype,iblk11)
         iblk11 = iblk11 + mvadd
         call rdedx(q(i7),nsq,iblk11,ifild)
      end if
c
c
c
      m9 = 9
      m12 = 12
      call secget(isect(9),m9,isec9)
      call rdedx(q(i18),lds(isect(9)),isec9,ifild)
      if (.not.uhf) then
         call dcopy(ncoorb,q(i18),1,q(i19),1)
      else
         call secget(isect(12),m12,isec12)
         call rdedx(q(i19),lds(isect(12)),isec12,ifild)
      end if
c
c
c
      call ump3sa(umma,amma,ammb,uds,ads,q(i6),q(i7),q(i8),q(i9),
     +  q(i10),q(i11),q(i12),q(i13),q(i14),q(i15),q(i16),q(i17),q(i18),
     +  q(i19),ncoorb,nocca,noccb,ndepa,ndepb,mna,mnb,ifort1,ifort2,
     +  ifort3)
c
c
      if (.not.uhf) goto 999
c
c
c     m8 = 8
c     m11 = 11
      itype = 0
      call secget(isect(8),itype,iblk8)
      iblk8 = iblk8 + mvadd
      call rdedx(q(i6),nsq,iblk8,ifild)
      itype = 0
      call secget(isect(11),itype,iblk11)
      iblk11 = iblk11 + mvadd
      call rdedx(q(i7),nsq,iblk11,ifild)
c
c
c
      m9 = 9
      m12 = 12
      call secget(isect(9),m9,isec9)
      call rdedx(q(i18),lds(isect(9)),isec9,ifild)
      call secget(isect(12),m12,isec12)
      call rdedx(q(i19),lds(isect(12)),isec12,ifild)
c
c
      call ump3sb(ummb,ammb,amma,uds,ads,q(i7),q(i6),q(i8),q(i9)
     +  ,q(i10),q(i11),q(i12),q(i13),q(i14),q(i15),q(i17),q(i16),q(i19)
     +  ,q(i18),ncoorb,noccb,nocca,ndepb,ndepa,mnb,mna,ifort4,ifort5,
     +  ifort6)
c
 999  call gmem_free(i19)
      call gmem_free(i18)
      call gmem_free(i17)
      call gmem_free(i16)
      call gmem_free(i15)
      call gmem_free(i14)
      call gmem_free(i13)
      call gmem_free(i12)
      call gmem_free(i11)
      call gmem_free(i10)
      call gmem_free(i9)
      call gmem_free(i8)
      call gmem_free(i7)
      call gmem_free(i6)
c
      return
 6010 format (/' insufficient core for umps3a ',i9,' real words ',
     +        ' need ',i9,' real words ')
      end
      subroutine uhfmp3(q,iq)
      implicit real*8  (a-h,o-z)
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer master, indxi, indxj, nfiles, junits, jblkrs, jblkas
      integer nfilef, junitf, jblkrf, jblkaf, isecbl, iscftp
      integer lword4, ilow4, ncol4, nsa4, newb4, nbas4, lenb4
      integer ndump4, iblkq4, nblkq4, lena4, nbb4, ionsv4, isecv4
      logical oprin4, oindx4
      integer npas41, npas42, iacc4
c
      common /cndx40/ master,indxi,indxj,nfiles,junits,jblkrs,jblkas,
     + nfilef,junitf,jblkrf,jblkaf,isecbl,iscftp,
     + lword4,ilow4,ncol4,nsa4,newb4,nbas4,lenb4,ndump4,iblkq4,
     + nblkq4,lena4,nbb4,
     + oprin4(10),oindx4,ionsv4,npas41,npas42,iacc4,isecv4
c
       integer len_cndx40
       parameter (len_cndx40=37)
c      used: restre(util1),revise(util1),utyp21(server)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c
      dimension q(*),iq(*)
c
      integer ixoxxa, ixoxxb
      integer ispin, iblkzz, iblkz, npassm, intblk, mupblk
      common /uhfspn/ ispin,iblkzz,iblkz,npassm,intblk(4),mupblk(8),
     +      ixoxxa,ixoxxb
c
      common/scfblk/enucl,etot,ehf
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      logical uhf
      character *8 open,energy
      data open /'open'/
      data energy /'hfscf'/
c
      l1 = newbas1
c
      uhf = scftyp.eq.open
      ifort1 = 17
      ifort2 = 18
      ifort3 = 19
      ifort4 = 20
      ifort5 = 21
      ifort6 = 22
      ifort7 = 23
      ifort8 = 24
      ifort9 = 9
      ifor10 = 10
      ifor11 = 11
      ibl11 = 1
c
      if (.not.uhf) then
         npassm = 1
         mupblk(1) = 1
         nvirtb = nvirta
      end if
      nsq = ncoorb*l1
      ntri = (ncoorb+1)*ncoorb/2
      ntrica = (nocca+1)*nocca/2
      ntricb = (noccb+1)*noccb/2
      ntriva = (nvirta+1)*nvirta/2
      ntrivb = (nvirtb+1)*nvirtb/2
      nlatri = (ntrica+1)*ntrica/2
      nlbtri = (ntricb+1)*ntricb/2
      mnatri = ntrica*ntriva
      mnbtri = ntricb*ntrivb
      mna = nocca*nvirta
      mnb = noccb*nvirtb
c
c  set up array positions
c
      iq1 = 1 + mnatri
      iq2 = iq1 + mnatri
      iq3 = iq2 + mnbtri
      iq4 = iq3 + mnbtri
      iq5 = iq4 + mna*mnb
      iq6 = iq5 + mna*mnb
      iq7 = iq6 + nsq
      iq8 = iq7 + nsq
      iq9 = iq8 + nsq
      iq10 = iq9 + ncoorb
      iq11 = iq10 + ncoorb
      itop = iq11 + nlatri
      maxq = igmem_max_memory()
c
      if (ncoorb.ne.nvirta+nocca.or.ncoorb.ne.nvirtb+noccb)
     1   call caserr(' mp3 with orbital restrictions is not allowed ')
c
      if (mprest.gt.2) then
c
c
c  ump3 gradient
c
         iq0 = igmem_alloc(mnatri)
         iq1 = igmem_alloc(mnatri)
         iq2 = igmem_alloc(mnbtri)
         iq3 = igmem_alloc(mnbtri)
         iq4 = igmem_alloc(mna*mnb)
         iq5 = igmem_alloc(mna*mnb)
c
         call rdedx(q(iq0),mnatri,ibl11,ifor11)
         call reads(q(iq1),mnatri,ifor11)
         call reads(q(iq2),mnbtri,ifor11)
         call reads(q(iq3),mnbtri,ifor11)
         call reads(q(iq4),mna*mnb,ifor11)
         call reads(q(iq5),mna*mnb,ifor11)
         call ump3gr(q,q(iq0),q(iq1),q(iq2),q(iq3),q(iq4),q(iq5),
     +               ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,ifort7,
     +               ifort8,ifort9,ifor10,ifor11,ibl11,iwr)
c
      else
c
c
         iblzz = 1
         iblki = 1
         call mpsrtj(q,iq,ncoorb,mupblk(iblki),junitf,ifort1)
         ipss = 1
         call mpsrtk(q,iq,ncoorb,mupblk(iblki),junitf,ifort2,ipss)
         if (uhf) then
            iblki = 2
            iblki1 = npassm*2 + 1
            ixoxxa = 0
            ixoxxb = 0
            call mpsrt6(q,iq,ncoorb,iblki,iblki1,junitf,ifort5,iblzz,
     +  iblkz,ifilz,ixoxxa,ixoxxb)
            call mpsrt6(q,iq,ncoorb,iblki1,iblki,junitf,ifort6,iblkz,
     +  iblzz,ifilz,ixoxxa,ixoxxb)
            iblki = npassm*2 + 2
            call mpsrtj(q,iq,ncoorb,mupblk(iblki),junitf,ifort3)
            ipss = 1
            call mpsrtk(q,iq,ncoorb,mupblk(iblki),junitf,ifort4,ipss)
         end if
c
         if (itop.gt.maxq) then
            write (iwr,6000) maxq , itop
            call caserr(' insufficient core for mp3 calculation')
         end if
c
         iq0  = igmem_alloc(mnatri)
         iq1  = igmem_alloc(mnatri)
         iq2  = igmem_alloc(mnbtri)
         iq3  = igmem_alloc(mnbtri)
         iq4  = igmem_alloc(mna*mnb)
         iq5  = igmem_alloc(mna*mnb)
         iq6  = igmem_alloc(nsq)
         iq7  = igmem_alloc(nsq)
         iq8  = igmem_alloc(nsq)
         iq9  = igmem_alloc(ncoorb)
         iq10 = igmem_alloc(ncoorb)
         iq11 = igmem_alloc(nlatri)
c
         m9 = 9
         m12 = 12
         call secget(isect(9),m9,isec9)
         call rdedx(q(iq9),lds(isect(9)),isec9,ifild)
         if (.not.uhf) then
            call dcopy(ncoorb,q(iq9),1,q(iq10),1)
         else
            call secget(isect(12),m12,isec12)
            call rdedx(q(iq10),lds(isect(12)),isec12,ifild)
         end if
c
c
c
         write(iwr,6010)
         e2aa = 0.0d0
         e2ab = 0.0d0
         e2ba = 0.0d0
         e2bb = 0.0d0
         eaa3 = 0.0d0
         ebb3 = 0.0d0
         eab3 = 0.0d0
         eba3 = 0.0d0
c
         call vclr(q(iq4),1,mna*mnb)
         ispin = 0
         call umpe3a(eaa3,eba3,e2aa,e2ba,q(iq6),q(iq7),q(iq8),q(iq9),
     +  q(iq10),q(iq11),q(iq0),q(iq1),q(iq4),q(iq5),nocca,noccb,nlatri,
     +  mnatri,ncoorb,mna,mnb,ifort1,ifort2,ifort5,junitf,uhf)
c
         if (uhf) then
            ispin = 1
            call umpe3b(ebb3,eab3,e2bb,e2ab,q(iq6),q(iq7),q(iq8),
     +  q(iq10),q(iq9),q(iq11),q(iq2),q(iq3),q(iq4),q(iq5),noccb,nocca,
     +  nlbtri,mnbtri,ncoorb,mnb,mna,ifort3,ifort4,ifort6,junitf)
         end if
c
c  ump2 energy
c
         e2 = -0.25d0*(e2aa+e2ba)*2.0d0
         write (iwr,6050)e2aa,e2ba
         if (uhf) then
            e2 = -0.25d0*(e2aa+e2ab+e2ba+e2bb)
            write (iwr,6060)e2ab,e2bb
         end if
         write (iwr,6080) e2
c
         write (iwr,6030)eaa3,eab3,eba3,ebb3
         e3 = eaa3*2.0d0 + eba3
         if (uhf) e3 = eaa3 + ebb3 + eab3
c
         etot = enucl + ehf + e3 + e2
         m13 = 13
         length = lensec(lds(isect(13)))
         call secput(isect(13),m13,length,isec13)
         call wrt3(enucl,lds(isect(13)),isec13,ifild)
         write (iwr,6040) e3 , etot
c
         ntri = ncoorb*(ncoorb+1)/2
         nsq = ncoorb*l1
         m8 = 8
         m11 = 11
         isec = isect(8)
c        itype = 0
         call secget(isec,m8,ibl)
         ibl = ibl + mvadd
         call rdedx(q(iq6),nsq,ibl,ifild)
         if (.not.uhf) then
            call dcopy(nsq,q(iq6),1,q(iq7),1)
         else
            isec = isect(11)
c           itype = 0
            call secget(isec,m11,ibl)
            ibl = ibl + mvadd
            call rdedx(q(iq7),nsq,ibl,ifild)
         end if
         m5 = 5
         isec = isect(5)
         call secget(isec,m5,ibl)
         call rdedx(q(iq8),ntri,ibl,ifild)
c
c  s squared expectation value
c
         call umpsq(q(iq6),q(iq7),q(iq8),nocca,noccb,l1,ntri,ssq)
c
         write(iwr,6070)ssq
c
         call gmem_free(iq11)
         call gmem_free(iq10)
         call gmem_free(iq9)
         call gmem_free(iq8)
         call gmem_free(iq7)
         call gmem_free(iq6)
c
         if (runtyp.ne.energy) then
            call wrt3(q(iq0),mnatri,ibl11,ifor11)
            call wrt3s(q(iq1),mnatri,ifor11)
            call wrt3s(q(iq2),mnbtri,ifor11)
            call wrt3s(q(iq3),mnbtri,ifor11)
            call wrt3s(q(iq4),mna*mnb,ifor11)
            call wrt3s(q(iq5),mna*mnb,ifor11)
         end if
      end if

      call gmem_free(iq5)
      call gmem_free(iq4)
      call gmem_free(iq3)
      call gmem_free(iq2)
      call gmem_free(iq1)
      call gmem_free(iq0)
      return
 6000 format (/' insufficient core for e3 ',i9,' real words ',' need ',
     +        i9,' real words ')
 6010 format (/10x,47('*')/
     +        10x,'mp3 calculation'/10x,47('*'))
 6050 format (
     +        10x,'second order contributions      '/10x,47('*')/
     +        10x,'e2(aa)                          ',f15.8/
     +        10x,'e2(ba)                          ',f15.8)
 6060 format (
     +        10x,'e2(ab)                          ',f15.8/
     +        10x,'e2(bb)                          ',f15.8)
 6080 format (
     +        10x,'second order perturbation energy',f15.8)
 6030 format (10x,47('*')/
     +        10x,'third order contributions       '/10x,47('*')/
     +        10x,'e3(aa)                          ',f15.8/
     +        10x,'e3(ab)                          ',f15.8/
     +        10x,'e3(ba)                          ',f15.8/
     +        10x,'e3(bb)                          ',f15.8)
 6040 format (
     +        10x,'third order perturbation energy ',f15.8/10x,47('*')/
     +        10x,'total energy (mp3)              ',f15.8/
     +        10x,47('*'))
 6070 format (
     +        10x,'expectation value of S**2       ',f15.8/)
      end
      subroutine ump3gr(q,umma,amma,ummb,ammb,uds,ads,
     1    ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,
     1    ifort7,ifort8,ifort9,ifor10,ifor11,ibl11,iw)
c
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
      dimension q(*), umma(*), amma(*), ummb(*), ammb(*), uds(*), ads(*)
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 rcigrd
      integer isecdd, isecll, ifil2d, iblk2d, iword, mnnr, mnc
      integer mnv, mnx, iscigr, isecmo
      integer isecnd, isecsy, irlagr, iadfrc, nfc, intlgr
      integer ncepa, ispaer
      integer nd2mo, ncore, ncact, nvr, ifilh, iblkh, iblk1
      integer ibl222, ntot, nupact, ijr3
      logical cigr, cicv, mpgr, mcgr, cicx, umpgr
      logical lcisd, lcepa, lacpf, lnewci, lsingl
      common /cigrad/ cigr,isecdd,isecll,ifil2d,iblk2d,iword,cicv,
     +                mnnr,mnc,mnv,mnx,mpgr,mcgr,cicx,iscigr,isecmo,
     +                isecnd,isecsy,irlagr,iadfrc,nfc,intlgr,umpgr,
     +                lcisd,lcepa,lacpf,lnewci,lsingl,ncepa,ispaer(20),
     +                nd2mo,ncore,ncact,nvr,ifilh,iblkh,iblk1,
     +                ibl222,ntot,nupact,ijr3,rcigrd(70)
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      logical uhf
      character *8 open
      data open/'open'/
      uhf = scftyp.eq.open
      mna = nocca*nvirta
      mnb = noccb*nvirtb
      mn = mna + mnb
      nija = nocca*(nocca+1)/2
      nijb = noccb*(noccb+1)/2
      naba = nvirta*(nvirta+1)/2
      nabb = nvirtb*(nvirtb+1)/2
      ntri = ncoorb*(ncoorb+1)/2
      nsq = ncoorb*ncoorb
      ifile1 = 1
      ndepa = naba*nija
      ndepb = nabb*nijb
c
c  1) ya(ncoorb,ncoorb)
c  2) wa(ncoorb,ncoorb)
c  3) yb(ncoorb,ncoorb)
c  4) wb(ncoorb,ncoorb)
c  5) yta(ncoorb,ncoorb)
c  6) ytb(ncoorb,ncoorb)
c  7) ytrans(nbf,nbf)
c  8) wtrans(nbf,nbf)
c  9) dma(ncoorb,ncoorb)
c 10) dmb(ncoorb,ncoorb)
c
      mpblk(1) = 1
      mpblk(2) = mpblk(1) + lensec(nsq)
      mpblk(3) = mpblk(2) + lensec(nsq)
      mpblk(4) = mpblk(3) + lensec(nsq)
      mpblk(5) = mpblk(4) + lensec(nsq)
      mpblk(6) = mpblk(5) + lensec(num*num)
      mpblk(7) = mpblk(6) + lensec(num*num)
      mpblk(8) = mpblk(7) + lensec(num*num)
      mpblk(9) = mpblk(8) + lensec(num*num)
      mpblk(10) = mpblk(9) + lensec(nsq)
      mpblk(11) = mpblk(10) + lensec(nsq)
      call wrt3z(1,ifile1,mpblk(11))
c
      i6 = 1
      i7 = i6 + nsq
      i8 = i7 + nsq
      i9 = i8 + nsq
      i10 = i9 + nsq
      i11 = i10 + nsq
      i12 = i11 + ncoorb
      i13 = i12 + ncoorb
      i14 = i13 + nsq
      i15 = i14 + nsq
      i16 = i15 + nsq
      itop = i16 + mna
      maxq = igmem_max_memory()
      if (itop.gt.maxq) then
         write (iw,6010) maxq , itop
         call caserr(' not enough core ')
      end if
c
      i6   = igmem_alloc(nsq)
      i7   = igmem_alloc(nsq)
      i8   = igmem_alloc(nsq)
      i9   = igmem_alloc(nsq)
      i10  = igmem_alloc(nsq)
      i11  = igmem_alloc(ncoorb)
      i12  = igmem_alloc(ncoorb)
      i13  = igmem_alloc(nsq)
      i14  = igmem_alloc(nsq)
      i15  = igmem_alloc(nsq)
      i16  = igmem_alloc(mna)
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(q(i11),ncoorb,isec9,ifild)
      if (.not.uhf) then
         call dcopy(ncoorb,q(i11),1,q(i12),1)
      else
         m12 = 12
         call secget(isect(12),m12,isec12)
         call rdedx(q(i12),ncoorb,isec12,ifild)
      end if
c
c  rmp3 case
c
      if (.not.uhf) then
         call dcopy(ndepa,umma,1,ummb,1)
         call dcopy(ndepa,amma,1,ammb,1)
      end if
      call umpmka(umma,amma,ammb,uds,ads,q(i6),q(i7),q(i8),
     +  q(i9),q(i10),q(i11),q(i12),q(i13),q(i14),q(i15),q(i16),nocca,
     +  noccb,ndepa,ndepb,mna,mnb,ncoorb,ifort1,ifort2,ifort5,ifort7,
     +  mpblk(1),mpblk(2),uhf)
c
c
      if (uhf) call delfil(ifort1)
      call delfil(ifort2)
c
      if (uhf) then
         call umpmkb(ummb,ammb,amma,uds,ads,q(i6),q(i7),q(i8),
     +  q(i9),q(i10),q(i12),q(i11),q(i13),q(i14),q(i15),q(i16),noccb,
     +  nocca,ndepb,ndepa,mnb,mna,ncoorb,ifort3,ifort4,ifort6,ifort8,
     +  mpblk(3),mpblk(4))
         call delfil(ifort3)
         call delfil(ifort4)
      end if
c
      call gmem_free(i16)
      call gmem_free(i15)
      call gmem_free(i14)
c
      if (.not.uhf) ifort6 = ifort1
      call umpm3(q(i6),q(i7),q(i8),q(i9),q(i10),q(i13),nocca,ncoorb,mna
     +  ,ifort7,ifort6,mpblk(1),mpblk(3),mpblk(2),uhf)
c
      if (uhf) then
         call umpm3(q(i6),q(i7),q(i8),q(i9),q(i10),q(i13+mna),noccb,
     +  ncoorb,mnb,ifort8,ifort5,mpblk(3),mpblk(1),mpblk(4),uhf)
      end if
c
      myuk = mn
      if (.not.uhf) myuk = mna
      call wrt3(q(i13),myuk,iblks,ifils)
c
c   first half of tpdm backtransformation
c
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(q(i11),ncoorb,isec9,ifild)
      if (uhf) then
         m12 = 12
         call secget(isect(12),m12,isec12)
         call rdedx(q(i12),ncoorb,isec12,ifild)
      end if
c
c  solve cpuhf equations to get z matrix
c
      np = 1
      npstar = 0
      jblk(1) = 1
c
      if (.not.uhf) mn = mna
      call cuhf3(q(i13),q(i11),q(i12),q(i7),q(i8),q(i9),q(i10),
     +  nocca,noccb,ncoorb,mna,mn,ifort7,ifort6,ifort5,ifort8,
     +  jblk(1),nofile(1),uhf)
c
      iblkz = iblks + lensec(myuk)
      call rdedx(q(i9),myuk,iblkz,ifils)
c
      call umpzy3(q(i7),q(i8),q(i9),nocca,noccb,ncoorb,myuk,mpblk(1),
     +  mpblk(3),uhf)
c
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(q(i11),ncoorb,isec9,ifild)
      if (uhf) then
         m12 = 12
         call secget(isect(12),m12,isec12)
         call rdedx(q(i12),ncoorb,isec12,ifild)
      end if
c
      if (.not.uhf) ifort5 = ifort6
      call umpmw3(q(i11),q(i6),q(i7),q(i8),q(i9),q(i10),nocca,ncoorb,
     +  ifort7,ifort5,mpblk(1),mpblk(3),mpblk(2),uhf)
      call delfil(ifort7)
      call delfil(ifort5)
c
      if (uhf) then
         call umpmw3(q(i12),q(i6),q(i7),q(i8),q(i9),q(i10),noccb,ncoorb
     +  ,ifort8,ifort6,mpblk(3),mpblk(1),mpblk(4),uhf)
         call delfil(ifort8)
         call delfil(ifort6)
      end if
c
c ya, yb, wa, wb stored on ed0 in mo basis
c
      if (.not.uhf) then
         mpblk(10) = mpblk(7)
         mpblk(9) = mpblk(6)
         mpblk(6) = mpblk(8)
         mpblk(8) = mpblk(5)
         mpblk(7) = mpblk(4)
         mpblk(5) = mpblk(3)
      end if
c
      ityp = 0
      len = lensec(ntri)
      call secget(isect(8),ityp,iblok)
      iblok = iblok + mvadd
      call rdedx(q(i6),num*ncoorb,iblok,ifild)
      call vclr(q(i10),1,nsq)
      call rdedx(q(i8),nsq,mpblk(1),ifile1)
      call vtamvu(q(i8),q(i6),q(i9),q(i10),ncoorb)
c
c     iblkya = iblkz + lensec(mn)
c     iblkyb = iblkya + lensec(nsq)
c
      call wrt3(q(i10),nsq,mpblk(5),ifile1)
      call vclr(q(i13),1,nsq)
c
      if (uhf) then
         ityp = 0
         call secget(isect(11),ityp,iblok)
         iblok = iblok + mvadd
         call rdedx(q(i7),num*ncoorb,iblok,ifild)
         call rdedx(q(i8),nsq,mpblk(3),ifile1)
         call vtamvu(q(i8),q(i7),q(i9),q(i13),ncoorb)
         call wrt3(q(i13),nsq,mpblk(6),ifile1)
      end if
c
c  add ump2 + ump3 for ya yb ys
c
      if (.not.uhf)  then
        call dcopy(nsq,q(i10),1,q(i13),1)
      endif
c
c
      ij = 0
      do 30 i = 1 , ncoorb
         do 20 j = 1 , ncoorb
c     q(i14+ij) = -q(i14+ij) + q(i10+ij)
c     q(i15+ij) = -q(i15+ij) + q(i13+ij)
            q(i10+ij) = q(i10+ij) + q(i13+ij)
            ij = ij + 1
 20      continue
 30   continue
      call wrt3(q(i10),nsq,mpblk(7),ifile1)
c
c
      call trsqsq(q(i10),q(i9),ncoorb)
c
c  pick up mp2 one pdm
c
      fact = 1.0d0
      if (.not.uhf) fact = 2.0d0
      do 40 ii = 1 , ntri
         q(i9+ii-1) = -q(i9+ii-1)
 40   continue
c
c
      call secput(isecdd,ityp,len,isdd)
      call wrt3(q(i9),ntri,isdd,ifild)
      call rdedx(q(i8),nsq,mpblk(2),ifile1)
      call vclr(q(i10),1,nsq)
      call vtamvu(q(i8),q(i6),q(i9),q(i10),ncoorb)
      if (uhf) then
         call rdedx(q(i8),nsq,mpblk(4),ifile1)
         call vtamvu(q(i8),q(i7),q(i9),q(i10),ncoorb)
      end if
      call wrt3(q(i10),nsq,mpblk(8),ifile1)
      call trsqsq(q(i10),q(i9),ncoorb)
      do 50 ii = 1 , ntri
         q(i9+ii-1) = -q(i9+ii-1)*fact
 50   continue
c
      call secput(isecll,ityp,len,isll)
      call wrt3(q(i9),ntri,isll,ifild)
c
      call gmem_free(i13)
      call gmem_free(i12)
      call gmem_free(i11)
      call gmem_free(i10)
      call gmem_free(i9)
      call gmem_free(i8)
      call gmem_free(i7)
      call gmem_free(i6)
c
c  second half of tpdm backtransformation
c
c
      call rdedx(umma,ndepa,ibl11,ifor11)
      call reads(amma,ndepa,ifor11)
      call reads(ummb,ndepb,ifor11)
      call reads(ammb,ndepb,ifor11)
      call reads(uds,mna*mnb,ifor11)
      call reads(ads,mna*mnb,ifor11)
      if (.not.uhf) then
         call dcopy(ndepa,umma,1,ummb,1)
         call dcopy(ndepa,amma,1,ammb,1)
      end if
      call delfil(ifor11)
      call ump3pd(q,umma,amma,ummb,ammb,uds,ads,
     +     ifort1,ifort2,ifort3,ifort4,ifort9,ifor10,iw)
c
      call ump3dm(q,iblk2d,ifil2d,mpblk(7),mpblk(9),mpblk(10),mpblk(5),
     +  mpblk(6),ifort1,ifort2,ifort3,ifort4,ifort9,ifor10,iw)
      if(oprint(44)) call whtps
      call delfil(ifort3)
      call delfil(ifort2)
      call delfil(ifort1)
      if (uhf) then
         call delfil(ifor10)
         call delfil(ifort9)
         call delfil(ifort4)
      end if
c
c
      return
 6010 format (/' insufficient core for umpmk ',i9,' real words ',
     +        ' need ',i9,' real words ')
      end
      subroutine umpzy3(ya,yb,z,nocca,noccb,ncoorb,mn,iblya,iblyb,uhf)
c
      implicit real*8  (a-h,o-z)
      dimension ya(ncoorb,ncoorb),yb(ncoorb,ncoorb),z(mn)
      logical uhf
c
      nsq = ncoorb*ncoorb
      ifile1 = 1
      call rdedx(ya,nsq,iblya,ifile1)
      if (uhf) call rdedx(yb,nsq,iblyb,ifile1)
      iai = 0
      do 30 ia = nocca + 1 , ncoorb
         do 20 i = 1 , nocca
            iai = iai + 1
            ya(ia,i) = z(iai)
 20      continue
 30   continue
      call wrt3(ya,nsq,iblya,ifile1)
c
      if (uhf) then
         do 50 ia = noccb + 1 , ncoorb
            do 40 i = 1 , noccb
               iai = iai + 1
               yb(ia,i) = z(iai)
 40         continue
 50      continue
         call wrt3(yb,nsq,iblyb,ifile1)
      end if
      return
      end
      subroutine ver_mp3(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp3.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
