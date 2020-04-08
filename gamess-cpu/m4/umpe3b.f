      subroutine umpe3b(e3ss,e3ds,e2ss,e2ds,t1,t2,t3,e,e1,tocc,u,a,
     1           uab,aab,nocc,noccz,nltri,mltri,ncoorb,mn,mnz,
     1               ifort1,ifort2,ifort3,kunit)
c
      implicit real*8  (a-h,o-z)
      dimension t1(ncoorb,ncoorb),t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),
     1 e(ncoorb),e1(ncoorb),a(mltri),u(mltri),tocc(nltri),
     1 uab(mnz,mn),aab(mnz,mn)
c
      integer ixoxxa, ixoxxb
      integer ispin, iblkzz, iblkz, npassm, intblk, mupblk
      common /uhfspn/ ispin,iblkzz,iblkz,npassm,intblk(4),mupblk(8),
     +      ixoxxa,ixoxxb
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
      common/blkin/yy(510),nint
      common/craypk/labs(1360)
      nsq = ncoorb*ncoorb
c
      thresh = 1.0d-8
      noc1 = nocc + 1
      nocc1 = noc1
      noccz1 = noccz + 1
      call vclr(u,1,mltri)
      call setsto(1360,0,labs)
c
c
      call vclr(tocc,1,nltri)
      call vclr(u,1,mltri)
      do 40 ipass = 1 , npassm
         iblki = npassm*2 + ipass*2
         call search(mupblk(iblki),kunit)
 20      call find(kunit)
         call get(yy,nw)
         if (nw.ne.0) then
c
            call unpack(yy(num2e+1),lab816,labs,numlab)
c
            do 30 int = 1 , nint
               nint4 = int + int + int + int

               k = labs(nint4-2)
               i = labs(nint4-3)
               l = labs(nint4  )
               j = labs(nint4-1)
               if (k.le.nocc) then
c
                  ki = (k-1)*k/2 + i
                  lj = (l-1)*l/2 + j
                  kilj = (ki-1)*ki/2 + lj
                  tocc(kilj) = yy(int)
               end if
 30         continue
            go to 20
         end if
 40   continue
c
      ijm = (nocc+1)*nocc/2
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      do 370 ip = 1 , ncoorb
         do 360 iq = 1 , ip
            call rdedz(t1,nsq,ifort1)
            call rdedz(t2,nsq,ifort2)
c
            call rdedz(t3,nsq,ifort3)
c
            if (ip.gt.nocc) then
               if (iq.gt.nocc) then
                  ic = ip
                  id = iq
                  eac = e(ip) + e(iq)
c
                  n = 2
c
                  do 80 i = 1 , nocc
                     do 70 j = 1 , i
                        difen = eac - e(i) - e(j)
c
                        val = (t2(i,j)-t2(j,i))/difen
                        do 60 ia = noc1 , ncoorb
                           do 50 ib = noc1 , ia
c
                              ij = (i-1)*i/2 + j
                              iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                              iiajb = (iab-1)*ijm + ij
c
                              u(iiajb) = u(iiajb)
     +                           + (t2(ia,ib)-t2(ib,ia))*val*n*0.125d0
 50                        continue
 60                     continue
 70                  continue
 80               continue
c
                  ia = ip
                  ib = iq
                  eab = e(ia) + e(ib)
                  do 120 k = 1 , nocc
                     do 110 l = 1 , k
                        do 100 i = 1 , nocc
                           do 90 j = 1 , i
                              ki = (k-1)*k/2 + i
                              if (i.gt.k) ki = (i-1)*i/2 + k
                              lj = (l-1)*l/2 + j
                              if (j.gt.l) lj = (j-1)*j/2 + l
                              kilj = (ki-1)*ki/2 + lj
                              li = (l-1)*l/2 + i
                              if (i.gt.l) li = (i-1)*i/2 + l
                              kj = (k-1)*k/2 + j
                              if (j.gt.k) kj = (j-1)*j/2 + k
                              kjli = (kj-1)*kj/2 + li
                              if (li.gt.kj) kjli = (li-1)*li/2 + kj
c
                              ij = (i-1)*i/2 + j
                              iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                              iiajb = (iab-1)*ijm + ij
                              difen = eab - e(k) - e(l)
                              val = (t2(k,l)-t2(l,k))/difen
c
                              valx = val*(tocc(kilj)-tocc(kjli))
     +                               *0.125d0*2.0d0
                              if (dabs(valx).lt.thresh) then
                              end if
                              u(iiajb) = u(iiajb)
     +                           + val*(tocc(kilj)-tocc(kjli))
     +                           *0.125d0*2.0d0
 90                        continue
 100                    continue
 110                 continue
 120              continue
               end if
            end if
c
c   third term
c
c   (ac\bd)aicjd
c
            if (ip.gt.noccz) then
               if (iq.gt.noccz) then
                  ib = ip
                  id = iq
                  do 160 ia = nocc1 , ib
                     maxcc = ia
                     if (ib.eq.ia) maxcc = id
                     do 150 ic = nocc1 , maxcc
                        iac = (ia-1)*ia/2 + ic
                        ibd = (ib-1)*ib/2 + id
                        if (iac.eq.ibd) go to 170
                        do 140 j = 1 , noccz
                           do 130 i = 1 , nocc
                              ijd = (id-noccz1)*noccz + j
                              ijb = (ib-noccz1)*noccz + j
                              iic = (ic-nocc1)*nocc + i
                              iia = (ia-nocc1)*nocc + i
c
                              val = aab(ijd,iic)
                              val1 = aab(ijb,iic)
                              val2 = aab(ijd,iia)
                              val3 = aab(ijb,iia)
c
c
                              uab(ijb,iia) = uab(ijb,iia)
     +                           + val*t3(ia,ic)*0.25d0
c
                              if (ib.ne.id) uab(ijd,iia) = uab(ijd,iia)
     +                            + val1*t3(ia,ic)*0.25d0
c
                              if (ia.ne.ic) uab(ijb,iic) = uab(ijb,iic)
     +                            + val2*t3(ia,ic)*0.25d0
c
                              if ((ia.ne.ic) .and. (ib.ne.id))
     +                            uab(ijd,iic) = uab(ijd,iic)
     +                            + val3*t3(ia,ic)*0.25d0
c
 130                       continue
 140                    continue
 150                 continue
 160              continue
               end if
            end if
c
 170        if (ip.gt.noccz) then
               if (iq.le.noccz) then
                  ic = ip
                  k = iq
                  ekc1 = e1(ic) - e1(k)
c
                  do 210 ia = noc1 , ncoorb
                     do 200 i = 1 , nocc
c
c  ump2 enery
c
                        val = t3(i,ia)
                        val = val*val
                        difen = ekc1 + e(ia) - e(i)
                        ff = 1.0d0/difen
                        e2ds = val*ff*2.0d0 + e2ds
c
                        do 190 ib = noc1 , ia
                           do 180 j = 1 , i
                              ij = (i-1)*i/2 + j
                              iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                              iiajb = (iab-1)*ijm + ij
c
c
                              difen = ekc1 + e(ib) - e(j)
                              val1 = t3(j,ib)/difen
c
                              difen = ekc1 + e(ia) - e(j)
                              val2 = -t3(j,ia)/difen
c
                              difen = ekc1 + e(ib) - e(i)
                              val3 = -t3(i,ib)/difen
c
                              difen = ekc1 + e(ia) - e(i)
                              val4 = t3(i,ia)/difen
c
c
c
                              u(iiajb) = u(iiajb)
     +                           + (val1*t3(i,ia)+val2*t3(i,ib)
     +                           +val3*t3(j,ia)+val4*t3(j,ib))*0.25d0
 180                       continue
 190                    continue
 200                 continue
 210              continue
               end if
            end if
c
c
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ic = ip
                  k = iq
                  ekc = e(ic) - e(k)
                  ikc = (ic-nocc1)*nocc + k
c
                  do 270 ia = noc1 , ncoorb
                     do 260 i = 1 , nocc
                        iia = (ia-nocc1)*nocc + i
c
c  ump2 energy
c
                        val = t1(i,ia) - t2(i,ia)
                        val = val*val
                        difen = ekc + e(ia) - e(i)
                        ff = 1.0d0/difen
                        e2ss = val*ff + e2ss
c
                        do 230 ib = noc1 , ia
                           do 220 j = 1 , i
                              difen = ekc + e(ib) - e(j)
                              val1 = (t1(j,ib)-t2(j,ib))/difen
c
c
                              ij = (i-1)*i/2 + j
                              iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                              iiajb = (iab-1)*ijm + ij
c
                              difen = ekc + e(ia) - e(j)
                              val2 = -(t1(j,ia)-t2(j,ia))/difen
c
                              difen = ekc + e(ib) - e(i)
                              val3 = -(t1(i,ib)-t2(i,ib))/difen
c
                              difen = ekc + e(ia) - e(i)
                              val4 = (t1(i,ia)-t2(i,ia))/difen
c
                              u(iiajb) = u(iiajb)
     +                           + (val1*(t1(i,ia)-t2(ia,i))
     +                           +val2*(t1(i,ib)-t2(ib,i))
     +                           +val3*(t1(j,ia)-t2(ia,j))
     +                           +val4*(t1(j,ib)-t2(ib,j)))*0.25d0
 220                       continue
 230                    continue
c
c  (kb\\jc)aiakc and (ka\\ac)akcjb
c
                        difen = ekc + e(ia) - e(i)
                        do 250 ib = noccz1 , ncoorb
                           do 240 j = 1 , noccz
                              ijb = (ib-noccz1)*noccz + j
                              val = t1(i,ia) - t2(ia,i)
c
                              uab(ijb,iia) = uab(ijb,iia) + aab(ijb,ikc)
     +                           *val*0.25d0
c
                              val = (t1(i,ia)-t2(i,ia))/difen
c
                              difen1 = ekc + e1(ib) - e1(j)
c
                              uab(ijb,iia) = uab(ijb,iia) + aab(ijb,ikc)
     +                           *difen1*val*0.25d0
 240                       continue
 250                    continue
c
c
                        if (ic.le.ia) then
                           if (k.le.i) then
c
                              ik = (i-1)*i/2 + k
                              iac = (ia-noc1)*(ia-nocc)/2 + ic - nocc
                              iiakc = (iac-1)*ijm + ik
c
                              difen = ekc + e(ia) - e(i)
                              valx = (t1(i,ia)-t2(i,ia))/difen
                              a(iiakc) = (t1(i,ia)-t2(i,ia))/difen
                           end if
                        end if
 260                 continue
 270              continue
               end if
            end if
c
c  (ki\bc)akajc
c
            if (ip.le.noccz) then
               if (iq.le.noccz) then
                  k = ip
                  j = iq
c
                  do 310 ia = nocc1 , ncoorb
                     do 300 ic = nocc1 , ncoorb
                        do 290 ib = noccz1 , ncoorb
                           do 280 i = 1 , nocc
                              iic = (ic-nocc1)*nocc + i
                              ikb = (ib-noccz1)*noccz + k
                              iia = (ia-nocc1)*nocc + i
                              ijb = (ib-noccz1)*noccz + j
                              uab(ijb,iia) = uab(ijb,iia) - t3(ia,ic)
     +                           *aab(ikb,iic)*0.25d0
c
                              if (k.ne.j) uab(ikb,iia) = uab(ikb,iia)
     +                            - t3(ia,ic)*aab(ijb,iic)*0.25d0
 280                       continue
 290                    continue
 300                 continue
 310              continue
               end if
            end if
c
c   (ki\lj)akalb
c
            if (ip.le.noccz) then
               if (iq.le.noccz) then
                  l = ip
                  j = iq
                  maxkk = l
                  if (l.gt.nocc) maxkk = nocc
                  do 350 k = 1 , maxkk
                     maxii = k
                     if (k.eq.l) maxii = j
                     do 340 i = 1 , maxii
                        ki = (k-1)*k/2 + i
                        lj = (l-1)*l/2 + j
                        if (ki.eq.lj) go to 360
                        do 330 ib = noccz1 , ncoorb
                           do 320 ia = nocc1 , ncoorb
                              ilb = (ib-noccz1)*noccz + l
                              ijb = (ib-noccz1)*noccz + j
                              ika = (ia-nocc1)*nocc + k
                              iia = (ia-nocc1)*nocc + i
c
                              val = aab(ilb,ika)
                              val1 = aab(ijb,ika)
                              val2 = aab(ilb,iia)
                              val3 = aab(ijb,iia)
c
c
                              uab(ijb,iia) = uab(ijb,iia) + val*t3(k,i)
     +                           *0.25d0
c
                              if (l.ne.j) uab(ilb,iia) = uab(ilb,iia)
     +                            + val1*t3(k,i)*0.25d0
c
                              if (k.ne.i) uab(ijb,ika) = uab(ijb,ika)
     +                            + val2*t3(k,i)*0.25d0
c
                              if ((k.ne.i) .and. (l.ne.j)) uab(ilb,ika)
     +                            = uab(ilb,ika) + val3*t3(k,i)*0.25d0
c
 320                       continue
 330                    continue
 340                 continue
 350              continue
               end if
            end if
c
c
c
 360     continue
 370  continue
c
c
      do 430 ia = noc1 , ncoorb
         do 420 i = 1 , nocc
            iia = (ia-nocc1)*nocc + i
            do 390 ib = noc1 , ia
               do 380 j = 1 , i
                  ij = (i-1)*i/2 + j
                  iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                  iiajb = (iab-1)*ijm + ij
c
c
                  n = 4
                  if ((ia.eq.ib) .or. (j.eq.i)) n = 2
                  if ((ia.eq.ib) .and. (j.eq.i)) n = 1
c
                  e3ss = e3ss + u(iiajb)*a(iiajb)*n
 380           continue
 390        continue
c
cc
c
            do 410 ib = noccz1 , ncoorb
               do 400 j = 1 , noccz
                  ijb = (ib-noccz1)*noccz + j
                  e3ds = e3ds + uab(ijb,iia)*aab(ijb,iia)*4.0d0
 400           continue
 410        continue
 420     continue
 430  continue
c
c
      return
      end
