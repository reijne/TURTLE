      subroutine umpe3a(e3ss,e3ds,e2ss,e2ds,t1,t2,t3,e,e1,tocc,u,a,
     1           uab,aab,nocc,noccz,nltri,mltri,ncoorb,mn,mnz,
     1               ifort1,ifort2,ifort3,kunit,uhf)
c
      implicit real*8  (a-h,o-z)
      dimension t1(ncoorb,ncoorb),t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),
     1 e(ncoorb),e1(ncoorb),a(mltri),u(mltri),tocc(nltri),
     1 uab(mn,mnz),aab(mn,mnz)
      logical uhf
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
      fac = 1.0d0
      if (.not.uhf) fac = 0.5d0
c
      thresh = 1.0d-8
      noc1 = nocc + 1
      nocc1 = noc1
      noccz1 = noccz + 1
      call vclr(u,1,mltri)
      call setsto(1360,0,labs)
      if (.not.uhf) then
         iyuk = ifort3
         ifort3 = ifort1
      end if
      call rewedz(ifort3)
c
c  form aab
c
      do 50 ip = 1 , ncoorb
         do 40 iq = 1 , ip
            call rdedz(t3,nsq,ifort3)
            if (ip.gt.noccz) then
               if (iq.le.noccz) then
                  ib = ip
                  j = iq
                  ejb = e1(ib) - e1(j)
                  ijb = (ib-noccz1)*noccz + j
                  do 30 ia = nocc1 , ncoorb
                     do 20 i = 1 , nocc
                        iia = (ia-nocc1)*nocc + i
                        difen = ejb + e(ia) - e(i)
                        aab(iia,ijb) = t3(i,ia)/difen
 20                  continue
 30               continue
               end if
            end if
 40      continue
 50   continue
c
c
      call vclr(tocc,1,nltri)
      call vclr(u,1,mltri)
      do 80 ipass = 1 , npassm
         iblki = ipass*2 - 1
         call search(mupblk(iblki),kunit)
 60      call find(kunit)
         call get(yy,nw)
         if (nw.ne.0) then
c
            call unpack(yy(num2e+1),lab816,labs,numlab)
c
            do 70 int = 1 , nint
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
 70         continue
            go to 60
         end if
 80   continue
c
      ijm = (nocc+1)*nocc/2
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      do 400 ip = 1 , ncoorb
         do 390 iq = 1 , ip
            call rdedz(t1,nsq,ifort1)
            call rdedz(t2,nsq,ifort2)
c
            if (.not.uhf) then
               call dcopy(nsq,t1,1,t3,1)
            else
               call rdedz(t3,nsq,ifort3)
            end if
c
c
            if (ip.gt.nocc) then
               if (iq.gt.nocc) then
                  ic = ip
                  id = iq
                  eac = e(ip) + e(iq)
c
                  n = 2
c
                  do 120 i = 1 , nocc
                     do 110 j = 1 , i
                        difen = eac - e(i) - e(j)
c
                        val = (t2(i,j)-t2(j,i))/difen
c
                        do 100 ia = noc1 , ncoorb
                           do 90 ib = noc1 , ia
c
                              ij = (i-1)*i/2 + j
                              iab = (ia-noc1)*(ia-nocc)/2 + ib - nocc
                              iiajb = (iab-1)*ijm + ij
                              u(iiajb) = u(iiajb)
     +                           + (t2(ia,ib)-t2(ib,ia))*val*n*0.125d0
 90                        continue
 100                    continue
 110                 continue
 120              continue
c
                  ia = ip
                  ib = iq
                  eab = e(ia) + e(ib)
                  do 160 k = 1 , nocc
                     do 150 l = 1 , k
                        do 140 i = 1 , nocc
                           do 130 j = 1 , i
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
c
                              valx = val*(tocc(kilj)-tocc(kjli))
     +                               *0.125d0*2.0d0
                              if (dabs(valx).lt.thresh) then
                              end if
                              u(iiajb) = u(iiajb)
     +                           + val*(tocc(kilj)-tocc(kjli))
     +                           *0.125d0*2.0d0
 130                       continue
 140                    continue
 150                 continue
 160              continue
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
                  do 200 ic = nocc1 , ib
                     maxaa = ic
                     if (ib.eq.ic) maxaa = id
                     do 190 ia = nocc1 , maxaa
                        iac = (ic-1)*ic/2 + ia
                        ibd = (ib-1)*ib/2 + id
                        ffac = 1.0d0
                        if (ibd.eq.iac) ffac = fac
                        do 180 j = 1 , noccz
                           do 170 i = 1 , nocc
                              ijd = (id-noccz1)*noccz + j
                              ijb = (ib-noccz1)*noccz + j
                              iic = (ic-nocc1)*nocc + i
                              iia = (ia-nocc1)*nocc + i
                              val = aab(iic,ijd)
                              val1 = aab(iic,ijb)
                              val2 = aab(iia,ijd)
                              val3 = aab(iia,ijb)
c
                              uab(iia,ijb) = uab(iia,ijb)
     +                           + val*t3(ia,ic)*0.25d0*ffac
c
                              if (ib.ne.id) uab(iia,ijd) = uab(iia,ijd)
     +                            + val1*t3(ia,ic)*0.25d0*ffac
c
                              if (ia.ne.ic) uab(iic,ijb) = uab(iic,ijb)
     +                            + val2*t3(ia,ic)*0.25d0*ffac
c
                              if ((ia.ne.ic) .and. (ib.ne.id))
     +                            uab(iic,ijd) = uab(iic,ijd)
     +                            + val3*t3(ia,ic)*0.25d0*ffac
c
 170                       continue
 180                    continue
 190                 continue
 200              continue
               end if
            end if
c
            if (ip.gt.noccz) then
               if (iq.le.noccz) then
                  ic = ip
                  k = iq
                  ekc1 = e1(ic) - e1(k)
c
                  do 240 ia = noc1 , ncoorb
                     do 230 i = 1 , nocc
c
c  ump2 enery
c
                        val = t3(i,ia)
                        val = val*val
                        difen = ekc1 + e(ia) - e(i)
                        ff = 1.0d0/difen
                        e2ds = val*ff*2.0d0 + e2ds
c
                        do 220 ib = noc1 , ia
                           do 210 j = 1 , i
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
 210                       continue
 220                    continue
 230                 continue
 240              continue
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
                  do 300 ia = noc1 , ncoorb
                     do 290 i = 1 , nocc
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
                        do 260 ib = noc1 , ia
                           do 250 j = 1 , i
                              difen = ekc + e(ib) - e(j)
                              val1 = (t1(j,ib)-t2(j,ib))/difen
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
 250                       continue
 260                    continue
c
c  (kb\\jc)aiakc and (ka\\ac)akcjb
c
                        difen = ekc + e(ia) - e(i)
                        do 280 ib = noccz1 , ncoorb
                           do 270 j = 1 , noccz
                              ijb = (ib-noccz1)*noccz + j
                              val = t1(i,ia) - t2(ia,i)
c
                              uab(iia,ijb) = uab(iia,ijb) + aab(ikc,ijb)
     +                           *val*0.25d0
c
                              val = (t1(i,ia)-t2(i,ia))/difen
c
                              difen1 = ekc + e1(ib) - e1(j)
c
                              uab(iia,ijb) = uab(iia,ijb) + aab(ikc,ijb)
     +                           *difen1*val*0.25d0
 270                       continue
 280                    continue
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
 290                 continue
 300              continue
               end if
            end if
c
c  (kj\ac)aickb
c
            if (ip.le.noccz) then
               if (iq.le.noccz) then
                  k = ip
                  j = iq
c
                  do 340 ia = nocc1 , ncoorb
                     do 330 ic = nocc1 , ncoorb
                        do 320 ib = noccz1 , ncoorb
                           do 310 i = 1 , nocc
                              iic = (ic-nocc1)*nocc + i
                              ikb = (ib-noccz1)*noccz + k
                              iia = (ia-nocc1)*nocc + i
                              ijb = (ib-noccz1)*noccz + j
                              uab(iia,ijb) = uab(iia,ijb) - t3(ia,ic)
     +                           *aab(iic,ikb)*0.25d0
                              if (k.ne.j) uab(iia,ikb) = uab(iia,ikb)
     +                            - t3(ia,ic)*aab(iic,ijb)*0.25d0
 310                       continue
 320                    continue
 330                 continue
 340              continue
               end if
            end if
c
c
c   (ki\lj)akalb
c
            if (ip.le.noccz) then
               if (iq.le.noccz) then
                  l = ip
                  j = iq
                  do 380 k = 1 , l
                     maxii = k
                     if (k.eq.l) maxii = j
                     do 370 i = 1 , maxii
                        ki = (k-1)*k/2 + i
                        lj = (l-1)*l/2 + j
                        fffac = 1.0d0
                        if (ki.eq.lj) fffac = fac
                        do 360 ia = nocc1 , ncoorb
                           do 350 ib = noccz1 , ncoorb
                              ilb = (ib-noccz1)*noccz + l
                              ijb = (ib-noccz1)*noccz + j
                              ika = (ia-nocc1)*nocc + k
                              iia = (ia-nocc1)*nocc + i
                              val = aab(ika,ilb)
                              val1 = aab(ika,ijb)
                              val2 = aab(iia,ilb)
                              val3 = aab(iia,ijb)
c
                              uab(iia,ijb) = uab(iia,ijb) + val*t3(k,i)
     +                           *0.25d0*fffac
c
                              if (l.ne.j) uab(iia,ilb) = uab(iia,ilb)
     +                            + val1*t3(k,i)*0.25d0*fffac
c
                              if (k.ne.i) uab(ika,ijb) = uab(ika,ijb)
     +                            + val2*t3(k,i)*0.25d0*fffac
c
                              if ((k.ne.i) .and. (l.ne.j)) uab(ika,ilb)
     +                            = uab(ika,ilb) + val3*t3(k,i)
     +                              *0.25d0*fffac
c
 350                       continue
 360                    continue
 370                 continue
 380              continue
               end if
            end if
c
c
c
 390     continue
 400  continue
c
c
      do 460 ia = noc1 , ncoorb
         do 450 i = 1 , nocc
            iia = (ia-nocc1)*nocc + i
            do 420 ib = noc1 , ia
               do 410 j = 1 , i
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
 410           continue
 420        continue
c
            if (.not.uhf) then
               do 440 ib = noccz1 , ncoorb
                  do 430 j = 1 , noccz
                     ijb = (ib-noccz1)*noccz + j
                     if (ijb.ge.iia) then
                        uab(iia,ijb) = uab(iia,ijb) + uab(ijb,iia)
                     else
                        uab(iia,ijb) = uab(ijb,iia)
                     end if
                     e3ds = e3ds + uab(iia,ijb)*aab(iia,ijb)*4.0d0
 430              continue
 440           continue
            end if
 450     continue
 460  continue
c
      if (.not.uhf) ifort3 = iyuk
      return
      end
