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
      implicit REAL  (a-h,o-z)
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
      implicit REAL  (a-h,o-z)
      logical uhf
INCLUDE(common/restar)
INCLUDE(common/atmblk)
_IFN1(iv)      common/craypk/labout(1360)
_IF1(iv)      common/craypk/labij(340),labkl(340)
      common/blkin/g(510),nint,nxtr
      dimension amaa(ncoorb,ncoorb),amab(ncoorb,ncoorb),
     1 amba(ncoorb,ncoorb),ambb(ncoorb,ncoorb),ea(ncoorb),
     1 eb(ncoorb),eps(mn)
      logical lstop,skipp
      dimension skipp(100)
      data m0/0/
c
_IF1()c     write(iwrt,*) ' entered cpuhf1  '
_IFN1(civ)      call izero(1360,labout,1)
_IF1(iv)      call setsto(680,0,labij)
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
_IF(ibm,vax)
                           labij(nint) = iai
                           labkl(nint) = ibj
                           if (nint.eq.num2e) then
                           call pak4v(labij,g(num2ep+1))
_ELSE
                           nint2 = nint + nint
                           labout(nint2-1) = iai
                           labout(nint2) = ibj
                           if (nint.eq.nintmx) then
                           call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
_ENDIF
                           call put(g,m511,ifw)
                           nint = 0
_IFN1(civ)                           call izero(1360,labout,1)
_IF1(iv)                           call setsto(680,0,labij)
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
_IF(ibm,vax)
                           labij(nint) = iai
                           labkl(nint) = mna + ibj
                           if (.not.uhf) labkl(nint) = ibj
                           if (nint.eq.nintmx) then
                            call pak4v(labij,g(num2ep+1))
_ELSE
                           nint2 = nint + nint
                           labout(nint2-1) = iai
                           labout(nint2) = mna + ibj
                           if (.not.uhf) labout(nint2) = ibj
                           if (nint.eq.nintmx) then
                            call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
_ENDIF
                            call put(g,m511,ifw)
                            nint = 0
_IFN1(civ)                            call izero(1360,labout,1)
_IF1(iv)                            call setsto(680,0,labij)
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
_IF(ibm,vax)
                           labij(nint) = mna + iai
                           labkl(nint) = mna + ibj
                           if (nint.eq.nintmx) then
                           call pak4v(labij,g(num2ep+1))
_ELSE
                           nint2 = nint + nint
                           labout(nint2-1) = mna + iai
                           labout(nint2) = mna + ibj
                           if (nint.eq.nintmx) then
                           call pack(g(num2ep+1),lab1632,
     +                               labout,numlabp)
_ENDIF
                           call put(g,m511,ifw)
                           nint = 0
_IFN1(civ)                           call izero(1360,labout,1)
_IF1(iv)                           call setsto(680,0,labij)
                           end if
                           end if
 100                    continue
 110                 continue
                  end if
               end if
            end if
 120     continue
 130  continue
_IFN1(iv)      call pack(g(num2ep+1),lab1632,labout,numlabp)
_IF1(iv)      call pak4v(labij,g(num2ep+1))
      call put(g,m511,ifw)
      call put(g,m0,ifw)
c
      skipp(1) = .false.
      lstop = .false.
c
      call chfdrv(eps,lstop,skipp)
      return
      end
_EXTRACT(umpe3a,umpe3a,hp700,hp800,pclinux)
      subroutine umpe3a(e3ss,e3ds,e2ss,e2ds,t1,t2,t3,e,e1,tocc,u,a,
     1           uab,aab,nocc,noccz,nltri,mltri,ncoorb,mn,mnz,
     1               ifort1,ifort2,ifort3,kunit,uhf)
c
      implicit REAL  (a-h,o-z)
      dimension t1(ncoorb,ncoorb),t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),
     1 e(ncoorb),e1(ncoorb),a(mltri),u(mltri),tocc(nltri),
     1 uab(mn,mnz),aab(mn,mnz)
      logical uhf
INCLUDE(common/uhfspn)
INCLUDE(common/atmblk)
      common/blkin/yy(510),nint
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
      nsq = ncoorb*ncoorb
      fac = 1.0d0
      if (.not.uhf) fac = 0.5d0
c
      thresh = 1.0d-8
      noc1 = nocc + 1
      nocc1 = noc1
      noccz1 = noccz + 1
      call vclr(u,1,mltri)
_IFN1(iv)      call setsto(1360,0,labs)
_IF1(iv)      call setsto(1360,0,i205)
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
_IF(ibm,vax)
            call upak8v(yy(num2e+1),i205)
_ELSE
            call unpack(yy(num2e+1),lab816,labs,numlab)
_ENDIF
c
            do 70 int = 1 , nint
_IFN1(iv)               nint4 = int + int + int + int

_IF(ibm,vax)
               k = i205(int)
               i = j205(int)
               l = k205(int)
               j = l205(int)
_ELSEIF(littleendian)
               k = labs(nint4-2)
               i = labs(nint4-3)
               l = labs(nint4  )
               j = labs(nint4-1)
_ELSE
               k = labs(nint4-3)
               i = labs(nint4-2)
               l = labs(nint4-1)
               j = labs(nint4)
_ENDIF
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
_ENDEXTRACT
_EXTRACT(umpe3b,umpe3b,hp700,hp800,pclinux)
      subroutine umpe3b(e3ss,e3ds,e2ss,e2ds,t1,t2,t3,e,e1,tocc,u,a,
     1           uab,aab,nocc,noccz,nltri,mltri,ncoorb,mn,mnz,
     1               ifort1,ifort2,ifort3,kunit)
c
      implicit REAL  (a-h,o-z)
      dimension t1(ncoorb,ncoorb),t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),
     1 e(ncoorb),e1(ncoorb),a(mltri),u(mltri),tocc(nltri),
     1 uab(mnz,mn),aab(mnz,mn)
INCLUDE(common/uhfspn)
INCLUDE(common/atmblk)
      common/blkin/yy(510),nint
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
      nsq = ncoorb*ncoorb
c
      thresh = 1.0d-8
      noc1 = nocc + 1
      nocc1 = noc1
      noccz1 = noccz + 1
      call vclr(u,1,mltri)
_IFN1(iv)      call setsto(1360,0,labs)
_IF1(iv)      call setsto(1360,0,i205)
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
_IF(ibm,vax)
            call upak8v(yy(num2e+1),i205)
_ELSE
            call unpack(yy(num2e+1),lab816,labs,numlab)
_ENDIF
c
            do 30 int = 1 , nint
_IFN1(iv)               nint4 = int + int + int + int

_IF(ibm,vax)
               k = i205(int)
               i = j205(int)
               l = k205(int)
               j = l205(int)
_ELSEIF(littleendian)
               k = labs(nint4-2)
               i = labs(nint4-3)
               l = labs(nint4  )
               j = labs(nint4-1)
_ELSE
               k = labs(nint4-3)
               i = labs(nint4-2)
               l = labs(nint4-1)
               j = labs(nint4)
_ENDIF
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
_IF1()c                             write(iwrt,*) i,ia,k,ic,iiakc,valx
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
_ENDEXTRACT
      subroutine umps3(veca,vecb,y,ysym,ifytr,ifdma,ifdmb,
     1                 nocca,noccb,ncoorb,iblkya,iblkyb,uhf)
      implicit REAL  (a-h,o-z)
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
      implicit REAL  (a-h,o-z)
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
_EXTRACT(umpmka,hp700,hp800)
      subroutine umpmka(umm,amm,ammz,uds,ads,t1,t2,t3,a1,a2,e1,e2,
     +           y,zl,dum,zint,nocc,noccz,ndep,ndepz,mn,mnz,ncoorb,
     +           ifort1,ifort2,ifort3,istrma,ibly,iblw,uhf)
c
      implicit REAL  (a-h,o-z)
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
_IF1()               go to 560
c
_IF1()c
_IF1()c  next kc term
_IF1()c
_IF1()c   (ba\kc)sibkc
_IF1()c
_IF1()c
_IF1()               ikc = (ic-noccz1)*noccz + k
_IF1()c
_IF1()               call vclr(zint,1,mn)
_IF1()               do 550 ib = nocc1 , ncoorb
_IF1()                  do 540 i = 1 , nocc
_IF1()                     do 510 id = nocc1 , ncoorb
_IF1()                        mfact = 1
_IF1()                        ibd = (ib-nocc1)*(ib-nocc)/2 + id - nocc
_IF1()                        if (id.gt.ib) then
_IF1()                           ibd = (id-nocc1)*(id-nocc)/2 + ib - nocc
_IF1()                           mfact = -mfact
_IF1()                        end if
_IF1()                        do 500 j = 1 , nocc
_IF1()                           ijd = (id-nocc1)*nocc + j
_IF1()                           ij = (i-1)*i/2 + j
_IF1()                           mfact2 = mfact
_IF1()                           if (j.gt.i) then
_IF1()                              ij = (j-1)*j/2 + i
_IF1()                              mfact2 = -mfact
_IF1()                           end if
_IF1()                           ibdij = (ibd-1)*ijm + ij
_IF1()                           ibi = (ib-nocc1)*nocc + i
_IF1()c
_IF1()                           zint(ibi) = zint(ibi) + amm(ibdij)
_IF1()     +                                 *ads(ijd,ikc)*mfact2
_IF1() 500                    continue
_IF1() 510                 continue
_IF1()c
_IF1()                     do 530 id = noccz1 , ncoorb
_IF1()                        mfact = 1
_IF1()                        icd = (ic-noccz1)*(ic-noccz)/2 + id - noccz
_IF1()                        if (id.gt.ic) then
_IF1()                           icd = (id-noccz1)*(id-noccz)/2 + ic - noccz
_IF1()                           mfact = -mfact
_IF1()                        end if
_IF1()                        do 520 j = 1 , noccz
_IF1()                           ijd = (id-noccz1)*noccz + j
_IF1()                           ibi = (ib-nocc1)*nocc + i
_IF1()                           kj = (k-1)*k/2 + j
_IF1()                           mfact2 = mfact
_IF1()                           if (j.gt.k) then
_IF1()                              kj = (j-1)*j/2 + k
_IF1()                              mfact2 = -mfact
_IF1()                           end if
_IF1()                           icdkj = (icd-1)*ijmz + kj
_IF1()c
_IF1()                           zint(ibi) = zint(ibi) + ads(ibi,ijd)
_IF1()     +                                 *ammz(icdkj)*mfact2
_IF1() 520                    continue
_IF1() 530                 continue
_IF1() 540              continue
_IF1() 550           continue
_IF1()c
_IF1()               call mxmb(t3(1,nocc1),1,ncoorb,zint,nocc,1,zl,ncoorb,1,
_IF1()     +                   ncoorb,nvirt,nocc)
            end if
c
c  now bd term
c
c  (ic\bd)sabcd
c
_IF1()560        if (iq.gt.noccz) then
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
_ENDEXTRACT
      subroutine umpmkb(umm,amm,ammz,uds,ads,t1,t2,t3,a1,a2,e1,e2,
     1               y,zl,dum,zint,nocc,noccz,ndep,ndepz,mn,mnz,ncoorb,
     1                ifort1,ifort2,ifort3,istrma,ibly,iblw)
c
      implicit REAL  (a-h,o-z)
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
_IF1()               go to 560
_IF1()c
_IF1()c  next kc term
_IF1()c
_IF1()c   (ba\kc)sibkc
_IF1()c
_IF1()               ikc = (ic-noccz1)*noccz + k
_IF1()c
_IF1()               call vclr(zint,1,mn)
_IF1()               do 550 ib = nocc1 , ncoorb
_IF1()                  do 540 i = 1 , nocc
_IF1()                     do 510 id = nocc1 , ncoorb
_IF1()                        mfact = 1
_IF1()                        ibd = (ib-nocc1)*(ib-nocc)/2 + id - nocc
_IF1()                        if (id.gt.ib) then
_IF1()                           ibd = (id-nocc1)*(id-nocc)/2 + ib - nocc
_IF1()                           mfact = -mfact
_IF1()                        end if
_IF1()                        do 500 j = 1 , nocc
_IF1()                           ijd = (id-nocc1)*nocc + j
_IF1()                           ij = (i-1)*i/2 + j
_IF1()                           mfact2 = mfact
_IF1()                           if (j.gt.i) then
_IF1()                              ij = (j-1)*j/2 + i
_IF1()                              mfact2 = -mfact
_IF1()                           end if
_IF1()                           ibdij = (ibd-1)*ijm + ij
_IF1()                           ibi = (ib-nocc1)*nocc + i
_IF1()c
_IF1()                           zint(ibi) = zint(ibi) + amm(ibdij)
_IF1()     +                                 *ads(ikc,ijd)*mfact2
_IF1() 500                    continue
_IF1() 510                 continue
_IF1()c
_IF1()                     do 530 id = noccz1 , ncoorb
_IF1()                        mfact = 1
_IF1()                        icd = (ic-noccz1)*(ic-noccz)/2 + id - noccz
_IF1()                        if (id.gt.ic) then
_IF1()                           icd = (id-noccz1)*(id-noccz)/2 + ic - noccz
_IF1()                           mfact = -mfact
_IF1()                        end if
_IF1()                        do 520 j = 1 , noccz
_IF1()                           ijd = (id-noccz1)*noccz + j
_IF1()                           ibi = (ib-nocc1)*nocc + i
_IF1()                           kj = (k-1)*k/2 + j
_IF1()                           mfact2 = mfact
_IF1()                           if (j.gt.k) then
_IF1()                              kj = (j-1)*j/2 + k
_IF1()                              mfact2 = -mfact
_IF1()                           end if
_IF1()                           icdkj = (icd-1)*ijmz + kj
_IF1()c
_IF1()                           zint(ibi) = zint(ibi) + ads(ijd,ibi)
_IF1()     +                                 *ammz(icdkj)*mfact2
_IF1() 520                    continue
_IF1() 530                 continue
_IF1() 540              continue
_IF1() 550           continue
_IF1()c
_IF1()               call mxmb(t3(1,nocc1),1,ncoorb,zint,nocc,1,zl,ncoorb,1,
_IF1()     +                   ncoorb,nvirt,nocc)
            end if
c
c  now bd term
c
c  (ic\bd)sabcd
c
c
_IF1() 560        if (iq.gt.noccz) then
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
      implicit REAL  (a-h,o-z)
      dimension sq(*),t(*)
      ij = 0
      ii = 0
      do 30 i = 1 , n
         jj = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
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
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical uhf
INCLUDE(common/nshel)
INCLUDE(common/restar)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IFN1(iv)      common/craypk/labout(1360)
_IF1(iv)      common/craypk/labij(340),labkl(340)
      common/blk1/g(510),nint,nxtr
INCLUDE(common/symtry)
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
_IFN1(civ)      call izero(1360,labout,1)
_IF1(iv)      call setsto(680,0,labij)
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
_IF(ibm,vax)
                                      labij(nint) = j1 + i4096(i1)
                                      labkl(nint) = l1 + i4096(k1)
                                      if (nint.eq.nintmx) then
                                      call pak4v(labij,g(num2e+1))
_ELSE
                                      nint4 = nint + nint + nint + nint
                                      labout(nint4-3) = i1
                                      labout(nint4-2) = j1
                                      labout(nint4-1) = k1
                                      labout(nint4) = l1
                                      if (nint.eq.nintmx) then
                                      call pack(g(num2e+1),lab816,
     +                                          labout,numlab)
_ENDIF
                                      call put(g,m511,ifw)
                                      nint = 0
_IFN1(civ)                                    call izero(1360,labout,1)
_IF1(iv)                                        call setsto(680,0,labij)
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
_IFN1(iv)      call pack(g(num2e+1),lab816,labout,numlab)
_IF1(iv)      call pak4v(labij,g(num2e+1))
      call put(g,m511,ifw)
      call put(g,mzero,ifw)
      call clredx
      return
      end
_EXTRACT(ump3sa,hp700,hp800)
      subroutine ump3sa(umm,amm,ammz,uds,ads,vec,vecz,dum,zint,zint1,
     1                  zint2,a1,a2,s1,s2,yint,yintz,e,ez,
     1      ncoorb,nocc,noccz,ndep,ndepz,mn,mnz,ifort1,ifort2,ifort3)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension umm(ndep),amm(ndep),ammz(ndepz),uds(mn,mnz),ads(mn,mnz),
     1        vec(ncoorb,ncoorb),vecz(ncoorb,ncoorb),dum(ncoorb,ncoorb),
     1    zint(ncoorb*ncoorb),zint1(ncoorb*ncoorb),zint2(ncoorb*ncoorb),
     1      a1(ncoorb,ncoorb),
     1        a2(ncoorb,ncoorb),s1(ncoorb,ncoorb),s2(ncoorb,ncoorb),
     1   yint(mn),yintz(mnz),e(ncoorb),ez(ncoorb)
INCLUDE(common/nshel)
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
_ENDEXTRACT
      subroutine ump3sb(umm,amm,ammz,uds,ads,vec,vecz,dum,zint,zint1,
     1                  zint2,a1,a2,s1,s2,yint,yintz,e,ez,
     1      ncoorb,nocc,noccz,ndep,ndepz,mn,mnz,ifort1,ifort2,ifort3)
c
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension umm(ndep),amm(ndep),ammz(ndepz),uds(mnz,mn),ads(mnz,mn),
     1        vec(ncoorb,ncoorb),vecz(ncoorb,ncoorb),dum(ncoorb,ncoorb),
     1    zint(ncoorb*ncoorb),zint1(ncoorb*ncoorb),zint2(ncoorb*ncoorb),
     1      a1(ncoorb,ncoorb),
     1        a2(ncoorb,ncoorb),s1(ncoorb,ncoorb),s2(ncoorb,ncoorb),
     1   yint(mn),yintz(mnz),e(ncoorb),ez(ncoorb)
INCLUDE(common/nshel)
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
      implicit REAL  (a-h,o-z)
      logical dpres,gpres,uhf
      dimension q(*)
c
INCLUDE(common/sizes)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/nshel)
INCLUDE(common/atmblk)
INCLUDE(common/symtry)
c
      character *8 open
      data open/'open'/
      uhf = scftyp.eq.open
_IF1()c     write(iw,*) ' entered u3tpdm'
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
_EXTRACT(ump3pd,hp700,hp800)
      subroutine ump3pd(q,umma,amma,ummb,ammb,uds,ads,
     +           ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,iw)
c
      implicit REAL  (a-h,o-z)
      dimension q(*)
      dimension umma(*),amma(*),ummb(*),ammb(*),uds(*),ads(*)
c
INCLUDE(common/sizes)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/nshel)
INCLUDE(common/atmblk)
c
      logical uhf
      character *8 open
      data open/'open'/
      uhf = scftyp.eq.open
_IF1()c     write(iw,*) ' entered ump3pd'
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
_ENDEXTRACT
      subroutine uhfmp3(q,iq)
      implicit REAL  (a-h,o-z)
INCLUDE(common/cigrad)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/cndx40)
INCLUDE(common/sizes)
INCLUDE(common/harmon)
c
      dimension q(*),iq(*)
INCLUDE(common/uhfspn)
      common/scfblk/enucl,etot,ehf
INCLUDE(common/atmblk)
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
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension q(*), umma(*), amma(*), ummb(*), ammb(*), uds(*), ads(*)
c
INCLUDE(common/prints)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/atmblk)
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
      implicit REAL  (a-h,o-z)
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
