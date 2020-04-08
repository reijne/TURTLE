c*module ormas2  *deck fchcx1s
c     ------------------------------------------------------------------
      subroutine fchcx1s(si1,si2,index,nact,na,nb,ci,ab,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
      logical goparr,dskwrk,maswrk
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension si1(*),si2(*)
      dimension index((nact*(nact+1))/2+1)
      dimension ci(nci),ab(nci),q(nci),x(nx)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension iob(nact),lgmul(nsym,nsym),ktab(nsym)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast),lsbc(ibst)
      dimension iposa(na*(nact-na),nsym)
      dimension ipera(na*(nact-na),nsym)
      dimension iind1(na*(nact-na),nsym)
      dimension igroa(na*(nact-na),nsym)
      dimension immc(nsym)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(nsym+1,ibst+1)
c
      do ii=1,nci
          ab(ii) = 0.0d+00
      enddo
c
c  --- big loop over all alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_dlbnext(my_task)
      endif
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            if (goparr.and.kka.ne.my_task) goto 4899
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            do ii=1,nsym
               immc(ii)=0
            enddo
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) goto 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) goto 4870
c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
c
                     ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
              if (lbox2(ispa2).gt.iama(ispa2)) goto 4800
c
c  get group number
c
           call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
           nias = nast(iga2)
c
              call idpost(iacon2,na,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kasym = lsyma(kapos)
              kper1 = (-1)**jpera
              immc(kasym) = immc(kasym) + 1
              jspo = immc(kasym)
              iposa(jspo,kasym) = jpza2
              ipera(jspo,kasym) = kper1
              iind1(jspo,kasym) = ind
              igroa(jspo,kasym) = iga2
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 4800
c
c   determine the alpha string contribution to the matrix element.
c
              c = si1(ind)
c
              do 4712 ik=1,na
                 if (ik.eq.ia1) goto 4712
                 ion = iacon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 4712         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4690
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
c
c matrix element addition
c
                 d = c
                 do 4670 ik=1,nb
                    ion = ibcon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 4670            continue
c
                 t = d*kper1
                 ab(jci1) = ab(jci1) + t*ci(jci2)
                 ab(jci2) = ab(jci2) + t*ci(jci1)
 4680         continue
c
 4690         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4700         continue
c
c --  double alpha excitations start here  ---
c
 4800         continue
c
            if (ia1.eq.na) goto 4850
            if (jj.eq.nact) goto 4850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            iaes3 = 1
            do kk=1,ispa1-1
               iaes3 = iaes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispa1.
c
            do 3890 ispa3 = ispa1,nspace
               if (lbox1(ispa3).eq.0) goto 3887
               if (ispa3.eq.ispa1.and.ia1.eq.ieae) goto 3887
               ioc3 = lbox2(ispa3)
               if (ioc3.eq.0) goto 3887
c
               iaee3 = iaes3 + lbox2(ispa3)-1
               lbox3(ispa3) = lbox3(ispa3)-1
c
c  loop over electrons in ispa3, which are larger than ia1,
c  making sure it isn't the already excited electron.
c
               jsta3 = iaes3
               if (ispa3.eq.ispa1) then
                  jsta3=ia1
                  if (igap-1.eq.ia1) jsta3=jsta3+1
               endif
c
               do 3880 ia3=jsta3,iaee3
                  if (ia3.eq.igap-1) goto 3880
                  io3 = iacon2(ia3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igae3 = 0
                  do jik=1,ispa2-1
                     igae3 = igae3 + lbox2(jik)
                  enddo
                  do 3850 ispa4=ispa2,nspace
c
c  igas3, igae3 are electrons specifying ispa4 space electron limits.
c
                  lbox3(ispa4) = lbox3(ispa4) + 1
                  igas3 = igae3 + 1
                  igae3 = igae3 + lbox2(ispa4)
                  if (lbox3(ispa3).lt.iami(ispa3)) goto 3840
                  if (lbox3(ispa4).gt.iama(ispa4)) goto 3840
             if (ispa4.eq.ispa2.and.jj.eq.msta(ispa2+1)-1) goto 3840
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox3,iga3)
               nias3 = nast(iga3)
c
c  make gap information here.
c
                  igaa3 = max(igas3,igap)
c
                  if (lbox2(ispa4).eq.0) then
                     ista3 = msta(ispa4)
                     iend3 = msta(ispa4+1)-1
                  elseif (ispa4.eq.ispa2) then
                     ista3 = jj+1
                     iend3 = iacon2(igaa3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igae3) iend3=msta(ispa2+1)-1
c
                  else
                     ista3 = msta(ispa4)
                     iend3 = iacon2(igaa3)-1
                  endif
c
c  loop over gaps
c
                  do 3830 igap3=igaa3,igae3+1
c
                     do 3820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 3820
c
              call rede00(iacon2,ibcon1,na,ia3,igap3-1,jj3,jpera3)
              call idpost(ibcon1,na,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga3),ibcon2,jposa3)
              iper3 = (-1)**(jpera3+jpera)
              kapos3 = jposa3 + nias3
              jpza3 = lspa(kapos3)
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over beta strings of right symmetry and group.
c
              do 3700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga3).ne.1) goto 3700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci3 = jpza3 + ldisb(ksym,igb,iga3)
c
              do 3680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci3 = jci3 + 1
                 ab(jci1) = ab(jci1) + t*ci(jci3)
                 ab(jci3) = ab(jci3) + t*ci(jci1)
 3680         continue
c
 3700         continue
c
 3820                continue
c
                     ista3 = iacon2(igap3)+1
                     iend3 = iacon2(igap3+1)-1
                     if (igap3.eq.igae3) iend3=msta(ispa4+1)-1
 3830             continue
c
 3840             lbox3(ispa4) = lbox3(ispa4) - 1
 3850             continue
c
 3880          continue
c
               lbox3(ispa3) = lbox3(ispa3)+1
 3887          iaes3 = iaes3 + lbox3(ispa3)
 3890       continue
c
 4850                continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c
c  --- end of loop over single alpha excitations. ---
c      now to sort them by positions within symmetries.
c
            do ii=1,nsym
               call fccsrt3(igroa(1,ii),ipera(1,ii),iind1(1,ii),
     *                   iposa(1,ii),immc(ii))
            enddo
c
c  --- end of loop over pure alpha excitations.
c  now to loop over all simultaneous ab -> a'b' excitations.
c
      if (nspace.eq.1) goto 3400
c
c  ***** general case of more than one space !!!!! *******
c
c
      if (.not.fdirct) then
c
c  loop over single alpha excites within each symmetry.
c
      do 3000 isae=1,nsym
         kbsym=ktab(isae)
         do 2900 jsae=1,immc(isae)
            jposae=iposa(jsae,isae)
            jperae=ipera(jsae,isae)
            jindae=iind1(jsae,isae)
            jgroae=igroa(jsae,isae)
c
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym.and.iga.eq.jgroae) then
c
c  loop over all relevant betas
c
            labpos = jpza1
            labpos2 = jposae
            do 2813 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2813
               nibs = nbst(igb)
               do 2763 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  labpos2 = labpos2 + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2613 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2613
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  jcib2=jpza1+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
                  ab(labpos2) = ab(labpos2) + c*ci(jcib2)
                  ab(jcib2) = ab(jcib2) + c*ci(labpos2)
c
 2613          continue
c
 2763          continue
 2813       continue
c
       else
c
c  loop over all relevant (a-)beta dets.
c
            labpos = jpza1
            do 2800 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2800
               nibs = nbst(igb)
               do 2750 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2600 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2600
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2600          continue
c
 2750          continue
 2800       continue
c
c
c  loop over all relevant (a'-)beta dets.
c
            labpos = jposae
            do 2803 igb=1,itgb
               if (lgcom(igb,jgroae).ne.1) goto 2803
               nibs = nbst(igb)
               do 2753 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2603 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,iga).ne.1) goto 2603
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jpza1+ldisb(ksym,igb2,iga)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2603          continue
c
 2753          continue
 2803       continue
c
      endif
c
 2900    continue
 3000 continue
c
      else
c
c  ***** direct method below *****
c
      do 4000 isae=1,nsym
         kbsym = ktab(isae)
c
c  analyse excited a' groups for compatibility with b.
c
      ngrps = 0
      do 1976 ii=1,immc(isae)
         icgr = igroa(ii,isae)
         do jj=1,ngrps
            if (jb1gr(jj).eq.icgr) goto 1976
         enddo
         ngrps = ngrps + 1
         jb1gr(ngrps) = icgr
 1976 continue
c
c first type of betas, ksym -> kbsym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7400 igb=1,itgb
             if (lgcom(igb,iga).ne.1) goto 7399
             lab1 = jpza1 + ldisb(ksym,igb,iga)
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
             do 7300 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                labpos = lab1 + lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry kbsym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7290 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7290
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7285 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7280 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7270
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7270
c
c  get group number
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          nias = nbst(igb2)
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb2,iagrp).eq.1) goto 7555
             enddo
             goto 7270
 7555        continue
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7260 igap=igba,igbe+1
c
                     do 7250 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7250
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jb1p = lspb(iposb + nias)
c
            do 3900 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb2,jgroae).ne.1) goto 3900
c
               jposae=iposa(jsae,isae)
               jcib = jb1p + jposae + ldisb(kbsym,igb2,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
            ab(labpos) = ab(labpos) + c*ci(jcib)
            ab(jcib) = ab(jcib) + c*ci(labpos)
c
 3900       continue
c
 7250                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7260             continue
c
 7270              lbox3(ispb2) = lbox3(ispb2) - 1
 7280           continue
c
 7285           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7290           continue
c
 7300        continue
c
 7399        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7400     continue
c
c
c second type of betas, kbsym -> ksym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7801 igb=1,itgb
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb,iagrp).eq.1) goto 7655
             enddo
             goto 7799
 7655        continue
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
c
             do 7701 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                lab1 = lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry ksym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7691 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7691
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7685 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7681 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7671
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7671
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          if (lgcom(igb2,iga).ne.1) goto 7671
          jcibs = jpza1 + ldisb(ksym,igb2,iga)
          nias = nbst(igb2)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7660 igap=igba,igbe+1
c
                     do 7650 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7650
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jcib= jcibs + lspb(iposb + nias)
c
            do 4300 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb,jgroae).ne.1) goto 4300
c
               jposae=iposa(jsae,isae)
               labpos = lab1 + jposae + ldisb(kbsym,igb,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
            ab(labpos) = ab(labpos) + c*ci(jcib)
            ab(jcib) = ab(jcib) + c*ci(labpos)
c
 4300       continue
c
 7650                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7660             continue
c
 7671              lbox3(ispb2) = lbox3(ispb2) - 1
 7681           continue
c
 7685           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7691           continue
c
 7701        continue
c
 7799        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7801     continue
c
 4000 continue
c
      endif
c
c  ***** end of direct option ******
c
c  **** end of general case of more than one space ******
c
      goto 4898
c
 3400 continue
c
c ***** special case of one space !!!!! ******
c
c  loop over single alpha excites within each symmetry.
c
       do 2901 isae=1,nsym
          kbsym=ktab(isae)
          do 2801 jsae=1,immc(isae)
                jposae=iposa(jsae,isae)
                jperae=ipera(jsae,isae)
                jindae=iind1(jsae,isae)
                jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym) then
c
c  loop over all relevant betas
c
             labpos=jpza1
             labpos2=jposae
                do 2721 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    labpos2 = labpos2 + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2621 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib = jposae+jb1po(jbindx)
                  jcib2 = jpza1+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
                  ab(labpos2) = ab(labpos2) + c*ci(jcib2)
                  ab(jcib2) = ab(jcib2) + c*ci(labpos2)
c
 2621          continue
c
 2721          continue
c
      else
c
c  loop over all relevant (a-)beta dets.
c
                labpos = jpza1
                do 2751 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2601 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2601          continue
c
 2751          continue
c
c  loop over all relevant (a'-)beta dets.
c
               labpos = jposae
               do 2761 kkb=lsbs(kbsym,1),lsbs(kbsym+1,1)-1
                  labpos = labpos + 1
                  ibpos = lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2611 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jpza1+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2611          continue
c
 2761          continue
c
      endif
c
 2801       continue
 2901   continue
c
c  **** end of special case of one space *******
c
 4898       continue
            if (goparr) call ddi_dlbnext(my_task)
 4899       continue
            call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings. ---
c
c     reset counters and sum up ab
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_gsumf(2900,ab,nci)
      endif
c
c **
c  --- loop over all pure beta excitations.
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox5)
c
      do 8000 igb=1,itgb
c
         call resetde(lbox1,nspace,nb,msta,ibcon1)
c
c  kkb gives the actual position of the beta string ibcon1 in
c  the full beta string list.
c
         do 7900 kkb=nbst(igb)+1,nbst(igb+1)
            jpzb1 = lspb(kkb)
            kbsym = lsymb(kkb)
            ksym=ktab(kbsym)
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            iebs = nb+1
            do 7890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) goto 7890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1.
c  iebs, iebe are the electrons in space ispb1.
c
               do 7885 ib1=iebe,iebs,-1
                  io1 = ibcon1(ib1)
                  igbe = iebe - lbox1(ispb1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 7880 ispb2=ispb1,nspace
c
c  igbs, igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) goto 7870
c
c  make gap information here.
c
                  igbb = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igbb)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igbb)-1
                  endif
c
c  loop over gaps
c
                  do 7860 igap=igbb,igbe+1
c
                     do 7850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
                        ind = index(jj) + io1
c
              call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,jperb)
                  if (lbox2(ispb2).gt.ibma(ispb2)) goto 7800
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 7800
c
c  get group number
c
           call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox2,igb2)
           nibs = nbst(igb2)
c
              call idpost(ibcon2,nb,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon1,jposb)
              kbpos = jposb + nibs
              jpzb2 = lspb(kbpos)
              kper1 = (-1)**jperb
c
c   determine the beta string contribution to the matrix element.
c
              c = si1(ind)
c
              do 7712 ik=1,nb
                 if (ik.eq.ib1) goto 7712
                 ion = ibcon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 7712         continue
c
c  loop over alpha strings of the right group and symmetry.
c
              call resetco(lbox3,nspace,na,iama,iami,lbox4)
c
              do 7700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7690
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              call resetde(lbox3,nspace,na,msta,iacon1)
              ista1 = 1
              do 7680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 do iiz=ista1,iena1-1
                    call moveup2(lbox3,nspace,na,msta,iacon1)
                 enddo
                 ista1 = iena1
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
c
c matrix element addition
c
                 d = c
                 do 7670 ik=1,na
                    ion = iacon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 7670            continue
c
                 t = d*kper1
                 ab(jci1) = ab(jci1) + t*ci(jci2)
                 ab(jci2) = ab(jci2) + t*ci(jci1)
 7680         continue
c
 7690         call pushco(lbox3,nspace,na,iama,iami,lbox4,iend)
 7700         continue
c
c  --- double excitations start here
c
 7800         continue
c
           if (ib1.eq.nb) goto 7850
           if (jj.eq.nact) goto 7850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            ibes3 = 1
            do kk=1,ispb1-1
               ibes3 = ibes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispb1.
c
            do 6890 ispb3 = ispb1,nspace
               if (lbox1(ispb3).eq.0) goto 6887
               if (ispb3.eq.ispb1.and.ib1.eq.iebe) goto 6887
               ioc3 = lbox2(ispb3)
               if (ioc3.eq.0) goto 6887
c
               ibee3 = ibes3 + lbox2(ispb3)-1
               lbox3(ispb3) = lbox3(ispb3)-1
c
c  loop over electrons in ispb3, which are larger than ib1,
c  making sure it isn't the already excited electron.
c
               jstb3 = ibes3
               if (ispb3.eq.ispb1) then
                  jstb3=ib1
                  if (igap-1.eq.ib1) jstb3=jstb3+1
               endif
c
               do 6880 ib3=jstb3,ibee3
                  if (ib3.eq.igap-1) goto 6880
                  io3 = ibcon2(ib3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igbe3 = 0
                  do jik=1,ispb2-1
                     igbe3 = igbe3 + lbox2(jik)
                  enddo
                  do 6850 ispb4=ispb2,nspace
c
c  igbs3, igbe3 are electrons specifying ispb4 space electron limits.
c
                  lbox3(ispb4) = lbox3(ispb4) + 1
                  igbs3 = igbe3 + 1
                  igbe3 = igbe3 + lbox2(ispb4)
                  if (lbox3(ispb3).lt.ibmi(ispb3)) goto 6840
                  if (lbox3(ispb4).gt.ibma(ispb4)) goto 6840
             if (ispb4.eq.ispb2.and.jj.eq.msta(ispb2+1)-1) goto 6840
c
c  get group number
c
               call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb3)
               nibs3 = nbst(igb3)
c
c  make gap information here.
c
                  igbb3 = max(igbs3,igap)
c
                  if (lbox2(ispb4).eq.0) then
                     ista3 = msta(ispb4)
                     iend3 = msta(ispb4+1)-1
                  elseif (ispb4.eq.ispb2) then
                     ista3 = jj+1
                     iend3 = ibcon2(igbb3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igbe3) iend3=msta(ispb2+1)-1
c
                  else
                     ista3 = msta(ispb4)
                     iend3 = ibcon2(igbb3)-1
                  endif
c
c  loop over gaps
c
                  do 6830 igap3=igbb3,igbe3+1
c
                     do 6820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 6820
c
              call rede00(ibcon2,iacon1,nb,ib3,igap3-1,jj3,jperb3)
              call idpost(iacon1,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb3),iacon2,jposb3)
              iper3 = (-1)**(jperb3+jperb)
              kbpos3 = jposb3 + nibs3
              jpzb3 = lspb(kbpos3)
c
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over alpha strings of right symmetry and group.
c
              do 6700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb3,iga).ne.1) goto 6700
              jcib1 = jpzb1 + ldisb(kbsym,igb,iga)
              jcib3 = jpzb3 + ldisb(kbsym,igb3,iga)
              nias = nast(iga)
c
              do 6680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena3 = lsac(kka)
                 jcia = lspa(nias+iena3)
                 jci1 = jcia + jcib1
                 jci3 = jcia + jcib3
c
                 ab(jci1) = ab(jci1) + t*ci(jci3)
                 ab(jci3) = ab(jci3) + t*ci(jci1)
 6680         continue
c
 6700         continue
c
 6820                continue
c
                     ista3 = ibcon2(igap3)+1
                     iend3 = ibcon2(igap3+1)-1
                     if (igap3.eq.igbe3) iend3=msta(ispb4+1)-1
 6830             continue
c
 6840             lbox3(ispb4) = lbox3(ispb4) - 1
 6850             continue
c
 6880          continue
c
               lbox3(ispb3) = lbox3(ispb3)+1
 6887          ibes3 = ibes3 + lbox3(ispb3)
 6890       continue
c
c  --- end of loop over double beta excitations.
c
 7850                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7860             continue
c
 7870             lbox2(ispb2) = lbox2(ispb2) - 1
 7880          continue
c
 7885          continue
c
               lbox2(ispb1) = lbox2(ispb1)+1
 7890       continue
c
            call moveup2(lbox1,nspace,nb,msta,ibcon1)
 7900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox5,iend)
 8000 continue
c
c  --- end of loop over beta strings. ---
c
c   now for the diagonal contributions
c
      do 119 ijk=1,nci
         ab(ijk) = ab(ijk) + q(ijk)*ci(ijk)
  119 continue
c
c  all done, return
c
      return
      end
c
c*module ormas2  *deck fchcxys
c     ------------------------------------------------------------------
      subroutine fchcxys(si1,si2,index,nact,na,nb,ci,ab,nv,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
      logical goparr,dskwrk,maswrk
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension si1(*),si2(*)
      dimension index((nact*(nact+1))/2+1)
      dimension ci(nci,nv),ab(nci,nv),q(nci),x(nx)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension iob(nact),lgmul(nsym,nsym),ktab(nsym)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast),lsbc(ibst)
      dimension iposa(na*(nact-na),nsym)
      dimension ipera(na*(nact-na),nsym)
      dimension iind1(na*(nact-na),nsym)
      dimension igroa(na*(nact-na),nsym)
      dimension immc(nsym)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(nsym+1,ibst+1)
c
      do ii=1,nci
         do jj=1,nv
            ab(ii,jj) = 0.0d+00
         enddo
      enddo
c
c  --- big loop over all alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_dlbnext(my_task)
      endif
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            if (goparr.and.kka.ne.my_task) goto 4899
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            do ii=1,nsym
               immc(ii)=0
            enddo
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) goto 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) goto 4870
c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
c
                     ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
                  if (lbox2(ispa2).gt.iama(ispa2)) goto 4800
c
c  get group number
c
           call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
           nias = nast(iga2)
c
              call idpost(iacon2,na,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kasym = lsyma(kapos)
              kper1 = (-1)**jpera
              immc(kasym) = immc(kasym) + 1
              jspo = immc(kasym)
              iposa(jspo,kasym) = jpza2
              ipera(jspo,kasym) = kper1
              iind1(jspo,kasym) = ind
              igroa(jspo,kasym) = iga2
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 4800
c
c   determine the alpha string contribution to the matrix element.
c
              c = si1(ind)
c
              do 4712 ik=1,na
                 if (ik.eq.ia1) goto 4712
                 ion = iacon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 4712         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4690
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
c
c matrix element addition
c
                 d = c
                 do 4670 ik=1,nb
                    ion = ibcon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 4670            continue
c
                 t = d*kper1
c
                 do 44 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci2,kj)
                 ab(jci2,kj) = ab(jci2,kj) + t*ci(jci1,kj)
   44            continue
c
 4680         continue
c
 4690         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4700         continue
c
c --  double alpha excitations start here  ---
c
 4800         continue
c
            if (ia1.eq.na) goto 4850
            if (jj.eq.nact) goto 4850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            iaes3 = 1
            do kk=1,ispa1-1
               iaes3 = iaes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispa1.
c
            do 3890 ispa3 = ispa1,nspace
               if (lbox1(ispa3).eq.0) goto 3887
               if (ispa3.eq.ispa1.and.ia1.eq.ieae) goto 3887
               ioc3 = lbox2(ispa3)
               if (ioc3.eq.0) goto 3887
c
               iaee3 = iaes3 + lbox2(ispa3)-1
               lbox3(ispa3) = lbox3(ispa3)-1
c
c  loop over electrons in ispa3, which are larger than ia1,
c  making sure it isn't the already excited electron.
c
               jsta3 = iaes3
               if (ispa3.eq.ispa1) then
                  jsta3=ia1
                  if (igap-1.eq.ia1) jsta3=jsta3+1
               endif
c
               do 3880 ia3=jsta3,iaee3
                  if (ia3.eq.igap-1) goto 3880
                  io3 = iacon2(ia3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igae3 = 0
                  do jik=1,ispa2-1
                     igae3 = igae3 + lbox2(jik)
                  enddo
                  do 3850 ispa4=ispa2,nspace
c
c  igas3, igae3 are electrons specifying ispa4 space electron limits.
c
                  lbox3(ispa4) = lbox3(ispa4) + 1
                  igas3 = igae3 + 1
                  igae3 = igae3 + lbox2(ispa4)
                  if (lbox3(ispa3).lt.iami(ispa3)) goto 3840
                  if (lbox3(ispa4).gt.iama(ispa4)) goto 3840
             if (ispa4.eq.ispa2.and.jj.eq.msta(ispa2+1)-1) goto 3840
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox3,iga3)
               nias3 = nast(iga3)
c
c  make gap information here.
c
                  igaa3 = max(igas3,igap)
c
                  if (lbox2(ispa4).eq.0) then
                     ista3 = msta(ispa4)
                     iend3 = msta(ispa4+1)-1
                  elseif (ispa4.eq.ispa2) then
                     ista3 = jj+1
                     iend3 = iacon2(igaa3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igae3) iend3=msta(ispa2+1)-1
c
                  else
                     ista3 = msta(ispa4)
                     iend3 = iacon2(igaa3)-1
                  endif
c
c  loop over gaps
c
                  do 3830 igap3=igaa3,igae3+1
c
                     do 3820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 3820
c
              call rede00(iacon2,ibcon1,na,ia3,igap3-1,jj3,jpera3)
              call idpost(ibcon1,na,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga3),ibcon2,jposa3)
              iper3 = (-1)**(jpera3+jpera)
              kapos3 = jposa3 + nias3
              jpza3 = lspa(kapos3)
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over beta strings of right symmetry and group.
c
              do 3700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga3).ne.1) goto 3700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci3 = jpza3 + ldisb(ksym,igb,iga3)
c
              do 3680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci3 = jci3 + 1
c
                 do 55 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci3,kj)
                 ab(jci3,kj) = ab(jci3,kj) + t*ci(jci1,kj)
   55            continue
c
 3680         continue
c
 3700         continue
c
 3820                continue
c
                     ista3 = iacon2(igap3)+1
                     iend3 = iacon2(igap3+1)-1
                     if (igap3.eq.igae3) iend3=msta(ispa4+1)-1
 3830             continue
c
 3840             lbox3(ispa4) = lbox3(ispa4) - 1
 3850             continue
c
 3880          continue
c
               lbox3(ispa3) = lbox3(ispa3)+1
 3887          iaes3 = iaes3 + lbox3(ispa3)
 3890       continue
c
 4850                continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c
c  --- end of loop over single alpha excitations. ---
c      now to sort them by positions within symmetries.
c
            do ii=1,nsym
               call fccsrt3(igroa(1,ii),ipera(1,ii),iind1(1,ii),
     *                   iposa(1,ii),immc(ii))
            enddo
c
c  --- end of loop over pure alpha excitations.
c  now to loop over all simultaneous ab -> a'b' excitations.
c
      if (nspace.eq.1) goto 3400
c
c  ***** general case of more than one space !!!!! *******
c
      if (.not.fdirct) then
c
c  loop over single alpha excites within each symmetry.
c
      do 3000 isae=1,nsym
         kbsym=ktab(isae)
         do 2900 jsae=1,immc(isae)
            jposae=iposa(jsae,isae)
            jperae=ipera(jsae,isae)
            jindae=iind1(jsae,isae)
            jgroae=igroa(jsae,isae)
c
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym.and.iga.eq.jgroae) then
c
c  loop over all relevant betas
c
            labpos = jpza1
            labpos2 = jposae
            do 2813 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2813
               nibs = nbst(igb)
               do 2763 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  labpos2 = labpos2 + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2613 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2613
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  jcib2=jpza1+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  do 66 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
                  ab(labpos2,kj) = ab(labpos2,kj) + c*ci(jcib2,kj)
                  ab(jcib2,kj) = ab(jcib2,kj) + c*ci(labpos2,kj)
   66             continue
c
 2613          continue
c
 2763          continue
 2813       continue
c
       else
c
c  loop over all relevant (a-)beta dets.
c
            labpos = jpza1
            do 2800 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2800
               nibs = nbst(igb)
               do 2750 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2600 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2600
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  do 77 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   77             continue
c
 2600          continue
c
 2750          continue
 2800       continue
c
c
c  loop over all relevant (a'-)beta dets.
c
            labpos = jposae
            do 2803 igb=1,itgb
               if (lgcom(igb,jgroae).ne.1) goto 2803
               nibs = nbst(igb)
               do 2753 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2603 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,iga).ne.1) goto 2603
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jpza1+ldisb(ksym,igb2,iga)+jb1po(jbindx)
c
                  do 88 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   88             continue
c
 2603          continue
c
 2753          continue
 2803       continue
c
      endif
c
 2900    continue
 3000 continue
c
      else
c
c  ***** direct method below *****
c
      do 4000 isae=1,nsym
         kbsym = ktab(isae)
c
c  analyse excited a' groups for compatibility with b.
c
      ngrps = 0
      do 1976 ii=1,immc(isae)
         icgr = igroa(ii,isae)
         do jj=1,ngrps
            if (jb1gr(jj).eq.icgr) goto 1976
         enddo
         ngrps = ngrps + 1
         jb1gr(ngrps) = icgr
 1976 continue
c
c first type of betas, ksym -> kbsym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7400 igb=1,itgb
             if (lgcom(igb,iga).ne.1) goto 7399
             lab1 = jpza1 + ldisb(ksym,igb,iga)
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
             do 7300 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                labpos = lab1 + lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry kbsym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7290 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7290
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7285 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7280 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7270
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7270
c
c  get group number
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          nias = nbst(igb2)
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb2,iagrp).eq.1) goto 7555
             enddo
             goto 7270
 7555        continue
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7260 igap=igba,igbe+1
c
                     do 7250 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7250
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jb1p = lspb(iposb + nias)
c
            do 3900 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb2,jgroae).ne.1) goto 3900
c
               jposae=iposa(jsae,isae)
               jcib = jb1p + jposae + ldisb(kbsym,igb2,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
         do kj=1,nv
            ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
            ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
         enddo
c
 3900       continue
c
 7250                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7260             continue
c
 7270              lbox3(ispb2) = lbox3(ispb2) - 1
 7280           continue
c
 7285           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7290           continue
c
 7300        continue
c
 7399        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7400     continue
c
c
c second type of betas, kbsym -> ksym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7801 igb=1,itgb
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb,iagrp).eq.1) goto 7655
             enddo
             goto 7799
 7655        continue
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
c
             do 7701 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                lab1 = lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry ksym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7691 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7691
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7685 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7681 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7671
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7671
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          if (lgcom(igb2,iga).ne.1) goto 7671
          jcibs = jpza1 + ldisb(ksym,igb2,iga)
          nias = nbst(igb2)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7660 igap=igba,igbe+1
c
                     do 7650 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7650
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jcib= jcibs + lspb(iposb + nias)
c
            do 4300 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb,jgroae).ne.1) goto 4300
c
               jposae=iposa(jsae,isae)
               labpos = lab1 + jposae + ldisb(kbsym,igb,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
         do kj=1,nv
            ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
            ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
         enddo
c
 4300       continue
c
 7650                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7660             continue
c
 7671              lbox3(ispb2) = lbox3(ispb2) - 1
 7681           continue
c
 7685           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7691           continue
c
 7701        continue
c
 7799        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7801     continue
c
 4000 continue
c
      endif
c
c  ***** end of direct option ******
c
c  **** end of general case of more than one space ******
c
      goto 4898
c
 3400 continue
c
c ***** special case of one space !!!!! ******
c
c  loop over single alpha excites within each symmetry.
c
       do 2901 isae=1,nsym
          kbsym=ktab(isae)
          do 2801 jsae=1,immc(isae)
                jposae=iposa(jsae,isae)
                jperae=ipera(jsae,isae)
                jindae=iind1(jsae,isae)
                jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym) then
c
c  loop over all relevant betas
c
             labpos=jpza1
             labpos2=jposae
                do 2721 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    labpos2 = labpos2 + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2621 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib = jposae+jb1po(jbindx)
                  jcib2 = jpza1+jb1po(jbindx)
c
                  do 56 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
                  ab(labpos2,kj) = ab(labpos2,kj) + c*ci(jcib2,kj)
                  ab(jcib2,kj) = ab(jcib2,kj) + c*ci(labpos2,kj)
   56             continue
c
 2621          continue
c
 2721          continue
c
      else
c
c  loop over all relevant (a-)beta dets.
c
                labpos = jpza1
                do 2751 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2601 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jposae+jb1po(jbindx)
c
                  do 67 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   67             continue
c
 2601          continue
c
 2751          continue
c
c  loop over all relevant (a'-)beta dets.
c
               labpos = jposae
               do 2761 kkb=lsbs(kbsym,1),lsbs(kbsym+1,1)-1
                  labpos = labpos + 1
                  ibpos = lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2611 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*jperae*jb1pe(jbindx)
                  jcib=jpza1+jb1po(jbindx)
c
                  do 78 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   78             continue
c
 2611          continue
c
 2761          continue
c
      endif
c
 2801       continue
 2901   continue
c
c  **** end of special case of one space *******
c
 4898       continue
            if (goparr) call ddi_dlbnext(my_task)
 4899       continue
            call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings. ---
c
c     reset counters and sum up ab
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_gsumf(2900,ab,nci*nv)
      endif
c
c  --- loop over all pure beta excitations.
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox5)
c
      do 8000 igb=1,itgb
c
         call resetde(lbox1,nspace,nb,msta,ibcon1)
c
c  kkb gives the actual position of the beta string ibcon1 in
c  the full beta string list.
c
         do 7900 kkb=nbst(igb)+1,nbst(igb+1)
            jpzb1 = lspb(kkb)
            kbsym = lsymb(kkb)
            ksym=ktab(kbsym)
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            iebs = nb+1
            do 7890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) goto 7890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1.
c  iebs, iebe are the electrons in space ispb1.
c
               do 7885 ib1=iebe,iebs,-1
                  io1 = ibcon1(ib1)
                  igbe = iebe - lbox1(ispb1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 7880 ispb2=ispb1,nspace
c
c  igbs, igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) goto 7870
c
c  make gap information here.
c
                  igbb = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igbb)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igbb)-1
                  endif
c
c  loop over gaps
c
                  do 7860 igap=igbb,igbe+1
c
                     do 7850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
                        ind = index(jj) + io1
c
              call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,jperb)
                  if (lbox2(ispb2).gt.ibma(ispb2)) goto 7800
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 7800
c
c  get group number
c
           call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox2,igb2)
           nibs = nbst(igb2)
c
              call idpost(ibcon2,nb,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon1,jposb)
              kbpos = jposb + nibs
              jpzb2 = lspb(kbpos)
              kper1 = (-1)**jperb
c
c   determine the beta string contribution to the matrix element.
c
              c = si1(ind)
c
              do 7712 ik=1,nb
                 if (ik.eq.ib1) goto 7712
                 ion = ibcon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 7712         continue
c
c  loop over alpha strings of the right group and symmetry.
c
              call resetco(lbox3,nspace,na,iama,iami,lbox4)
c
              do 7700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7690
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              call resetde(lbox3,nspace,na,msta,iacon1)
              ista1 = 1
              do 7680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 do iiz=ista1,iena1-1
                    call moveup2(lbox3,nspace,na,msta,iacon1)
                 enddo
                 ista1 = iena1
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
c
c matrix element addition
c
                 d = c
                 do 7670 ik=1,na
                    ion = iacon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 7670            continue
c
                 t = d*kper1
c
                 do 87 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci2,kj)
                 ab(jci2,kj) = ab(jci2,kj) + t*ci(jci1,kj)
   87            continue
c
 7680         continue
c
 7690         call pushco(lbox3,nspace,na,iama,iami,lbox4,iend)
 7700         continue
c
c  --- double excitations start here
c
 7800         continue
c
           if (ib1.eq.nb) goto 7850
           if (jj.eq.nact) goto 7850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            ibes3 = 1
            do kk=1,ispb1-1
               ibes3 = ibes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispb1.
c
            do 6890 ispb3 = ispb1,nspace
               if (lbox1(ispb3).eq.0) goto 6887
               if (ispb3.eq.ispb1.and.ib1.eq.iebe) goto 6887
               ioc3 = lbox2(ispb3)
               if (ioc3.eq.0) goto 6887
c
               ibee3 = ibes3 + lbox2(ispb3)-1
               lbox3(ispb3) = lbox3(ispb3)-1
c
c  loop over electrons in ispb3, which are larger than ib1,
c  making sure it isn't the already excited electron.
c
               jstb3 = ibes3
               if (ispb3.eq.ispb1) then
                  jstb3=ib1
                  if (igap-1.eq.ib1) jstb3=jstb3+1
               endif
c
               do 6880 ib3=jstb3,ibee3
                  if (ib3.eq.igap-1) goto 6880
                  io3 = ibcon2(ib3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igbe3 = 0
                  do jik=1,ispb2-1
                     igbe3 = igbe3 + lbox2(jik)
                  enddo
                  do 6850 ispb4=ispb2,nspace
c
c  igbs3, igbe3 are electrons specifying ispb4 space electron limits.
c
                  lbox3(ispb4) = lbox3(ispb4) + 1
                  igbs3 = igbe3 + 1
                  igbe3 = igbe3 + lbox2(ispb4)
                  if (lbox3(ispb3).lt.ibmi(ispb3)) goto 6840
                  if (lbox3(ispb4).gt.ibma(ispb4)) goto 6840
             if (ispb4.eq.ispb2.and.jj.eq.msta(ispb2+1)-1) goto 6840
c
c  get group number
c
               call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb3)
               nibs3 = nbst(igb3)
c
c  make gap information here.
c
                  igbb3 = max(igbs3,igap)
c
                  if (lbox2(ispb4).eq.0) then
                     ista3 = msta(ispb4)
                     iend3 = msta(ispb4+1)-1
                  elseif (ispb4.eq.ispb2) then
                     ista3 = jj+1
                     iend3 = ibcon2(igbb3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igbe3) iend3=msta(ispb2+1)-1
c
                  else
                     ista3 = msta(ispb4)
                     iend3 = ibcon2(igbb3)-1
                  endif
c
c  loop over gaps
c
                  do 6830 igap3=igbb3,igbe3+1
c
                     do 6820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 6820
c
              call rede00(ibcon2,iacon1,nb,ib3,igap3-1,jj3,jperb3)
              call idpost(iacon1,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb3),iacon2,jposb3)
              iper3 = (-1)**(jperb3+jperb)
              kbpos3 = jposb3 + nibs3
              jpzb3 = lspb(kbpos3)
c
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over alpha strings of right symmetry and group.
c
              do 6700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb3,iga).ne.1) goto 6700
              jcib1 = jpzb1 + ldisb(kbsym,igb,iga)
              jcib3 = jpzb3 + ldisb(kbsym,igb3,iga)
              nias = nast(iga)
c
              do 6680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena3 = lsac(kka)
                 jcia = lspa(nias+iena3)
                 jci1 = jcia + jcib1
                 jci3 = jcia + jcib3
c
                 do 85 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci3,kj)
                 ab(jci3,kj) = ab(jci3,kj) + t*ci(jci1,kj)
   85            continue
c
 6680         continue
c
 6700         continue
c
 6820                continue
c
                     ista3 = ibcon2(igap3)+1
                     iend3 = ibcon2(igap3+1)-1
                     if (igap3.eq.igbe3) iend3=msta(ispb4+1)-1
 6830             continue
c
 6840             lbox3(ispb4) = lbox3(ispb4) - 1
 6850             continue
c
 6880          continue
c
               lbox3(ispb3) = lbox3(ispb3)+1
 6887          ibes3 = ibes3 + lbox3(ispb3)
 6890       continue
c
c  --- end of loop over double beta excitations.
c
 7850                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7860             continue
c
 7870             lbox2(ispb2) = lbox2(ispb2) - 1
 7880          continue
c
 7885          continue
c
               lbox2(ispb1) = lbox2(ispb1)+1
 7890       continue
c
            call moveup2(lbox1,nspace,nb,msta,ibcon1)
 7900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox5,iend)
 8000 continue
c
c  --- end of loop over beta strings. ---
c
c
c   now for the diagonal contributions
c
      do 119 ijk=1,nci
         do 118 kj=1,nv
           ab(ijk,kj) = ab(ijk,kj) + q(ijk)*ci(ijk,kj)
  118    continue
  119 continue
c
c  all done, return
c
      return
      end
c
c*module ormas2  *deck fchc01s
c     ------------------------------------------------------------------
      subroutine fchc01s(si1,si2,index,nact,na,nb,ci,ab,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,idim1,idim2,
     *           ispin,ihmcon)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
      logical goparr,dskwrk,maswrk
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension si1(*),si2(*)
      dimension index((nact*(nact+1))/2+1)
      dimension ci(nci),ab(nci),q(nci),x(nx)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension iob(nact),lgmul(nsym,nsym),ktab(nsym)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsbs(nsym+1,itgb)
      dimension lsbc(ibst)
      dimension iposa(na*(nact-na),nsym)
      dimension ipera(na*(nact-na),nsym)
      dimension iind1(na*(nact-na),nsym)
      dimension igroa(na*(nact-na),nsym)
      dimension immc(nsym)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(idim1,idim2)
      dimension ispin(*),ihmcon(1)
c
      do ii=1,nci
          ab(ii) = 0.0d+00
      enddo
c
c  --- big loop over all alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_dlbnext(my_task)
      endif
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            if (goparr.and.kka.ne.my_task) goto 4899
            ilima = kka-nast(iga)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            do ii=1,nsym
               immc(ii)=0
            enddo
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) goto 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) goto 4870

c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
c
                     ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
                  if (lbox2(ispa2).gt.iama(ispa2)) goto 4800
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
               nias = nast(iga2)
c
              call idpost(iacon2,na,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kasym = lsyma(kapos)
              kper1 = (-1)**jpera
              immc(kasym) = immc(kasym) + 1
              jspo = immc(kasym)
              iposa(jspo,kasym) = jpza2
              ipera(jspo,kasym) = kper1
              iind1(jspo,kasym) = ind
              igroa(jspo,kasym) = iga2
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 4800
c
c   determine the alpha string contribution to the matrix element.
c
              c = si1(ind)
c
              do 4712 ik=1,na
                 if (ik.eq.ia1) goto 4712
                 ion = iacon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 4712         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4690
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
c
c matrix element addition
c
                 d = c
                 do 4670 ik=1,nb
                    ion = ibcon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 4670            continue
c
                 t = d*kper1
                 ab(jci1) = ab(jci1) + t*ci(jci2)
                 ab(jci2) = ab(jci2) + t*ci(jci1)
 4680         continue
c
 4690         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4700         continue
c
c --  double alpha excitations start here  ---
c
 4800         continue
c
            if (ia1.eq.na) goto 4850
            if (jj.eq.nact) goto 4850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            iaes3 = 1
            do kk=1,ispa1-1
               iaes3 = iaes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispa1.
c
            do 3890 ispa3 = ispa1,nspace
               if (lbox1(ispa3).eq.0) goto 3887
               if (ispa3.eq.ispa1.and.ia1.eq.ieae) goto 3887
               ioc3 = lbox2(ispa3)
               if (ioc3.eq.0) goto 3887
c
               iaee3 = iaes3 + lbox2(ispa3)-1
               lbox3(ispa3) = lbox3(ispa3)-1
c
c  loop over electrons in ispa3, which are larger than ia1,
c  making sure it isn't the already excited electron.
c
               jsta3 = iaes3
               if (ispa3.eq.ispa1) then
                  jsta3=ia1
                  if (igap-1.eq.ia1) jsta3=jsta3+1
               endif
c
               do 3880 ia3=jsta3,iaee3
                  if (ia3.eq.igap-1) goto 3880
                  io3 = iacon2(ia3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igae3 = 0
                  do jik=1,ispa2-1
                     igae3 = igae3 + lbox2(jik)
                  enddo
                  do 3850 ispa4=ispa2,nspace
c
c  igas3, igae3 are electrons specifying ispa4 space electron limits.
c
                  lbox3(ispa4) = lbox3(ispa4) + 1
                  igas3 = igae3 + 1
                  igae3 = igae3 + lbox2(ispa4)
                  if (lbox3(ispa3).lt.iami(ispa3)) goto 3840
                  if (lbox3(ispa4).gt.iama(ispa4)) goto 3840
             if (ispa4.eq.ispa2.and.jj.eq.msta(ispa2+1)-1) goto 3840
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox3,iga3)
               nias3 = nast(iga3)
c
c  make gap information here.
c
                  igaa3 = max(igas3,igap)
c
                  if (lbox2(ispa4).eq.0) then
                     ista3 = msta(ispa4)
                     iend3 = msta(ispa4+1)-1
                  elseif (ispa4.eq.ispa2) then
                     ista3 = jj+1
                     iend3 = iacon2(igaa3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igae3) iend3=msta(ispa2+1)-1
c
                  else
                     ista3 = msta(ispa4)
                     iend3 = iacon2(igaa3)-1
                  endif
c
c  loop over gaps
c
                  do 3830 igap3=igaa3,igae3+1
c
                     do 3820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 3820
c
              call rede00(iacon2,ibcon1,na,ia3,igap3-1,jj3,jpera3)
              call idpost(ibcon1,na,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga3),ibcon2,jposa3)
              iper3 = (-1)**(jpera3+jpera)
              kapos3 = jposa3 + nias3
              jpza3 = lspa(kapos3)
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over beta strings of right symmetry and group.
c
              do 3700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga3).ne.1) goto 3700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci3 = jpza3 + ldisb(ksym,igb,iga3)
c
              do 3680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci3 = jci3 + 1
                 ab(jci1) = ab(jci1) + t*ci(jci3)
                 ab(jci3) = ab(jci3) + t*ci(jci1)
 3680         continue
c
 3700         continue
c
 3820                continue
c
                     ista3 = iacon2(igap3)+1
                     iend3 = iacon2(igap3+1)-1
                     if (igap3.eq.igae3) iend3=msta(ispa4+1)-1
 3830             continue
c
 3840             lbox3(ispa4) = lbox3(ispa4) - 1
 3850             continue
c
 3880          continue
c
               lbox3(ispa3) = lbox3(ispa3)+1
 3887          iaes3 = iaes3 + lbox3(ispa3)
 3890       continue
c
 4850                continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c
c  --- end of loop over single alpha excitations. ---
c      now to sort them by positions within symmetries.
c
            do ii=1,nsym
               call fccsrt3(igroa(1,ii),ipera(1,ii),iind1(1,ii),
     *                   iposa(1,ii),immc(ii))
            enddo
c
c  --- end of loop over pure alpha excitations.
c  now to loop over all simultaneous ab -> a'b' excitations.
c
      if (nspace.eq.1) goto 3400
c
c  ***** general case of more than one space !!!!! *******
c
      if (.not.fdirct) then
c
c  loop over single alpha excites within each symmetry.
c
      do 3000 isae=1,nsym
         kbsym=ktab(isae)
c
         do 2900 jsae=1,immc(isae)
            jposae=iposa(jsae,isae)
            jperae=ipera(jsae,isae)
            jindae=iind1(jsae,isae)
            jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym.and.iga.eq.jgroae) then
c
c  loop over all relevant betas
c
            do 2813 igb=iga,itgb
               if (lgcom(igb,iga).ne.1) goto 2813
               lab1 = jpza1 + ldisb(ksym,igb,iga)
               lab2 = jposae + ldisb(ksym,igb,iga)
               nibs = nbst(igb)
               do 2763 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2763
                  labpos = lab1 + lspb(ibpos)
                  labpos2 = lab2 + lspb(ibpos)
                  qjmoda=jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2613 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2613
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  jcib2=jpza1+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
                  ab(labpos2) = ab(labpos2) + c*ci(jcib2)
                  ab(jcib2) = ab(jcib2) + c*ci(labpos2)
c
 2613          continue
c
 2763          continue
 2813       continue
c
       else
c
c  loop over all relevant (a-)beta dets.
c
            do 2800 igb=iga,itgb
               if (lgcom(igb,iga).ne.1) goto 2800
               lab1 = jpza1 + ldisb(ksym,igb,iga)
               nibs = nbst(igb)
               do 2750 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2750
                  labpos = lab1 + lspb(ibpos)
                  qjmoda = jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2600 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2600
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2600          continue
c
 2750          continue
 2800       continue
c
c  loop over all relevant (a'-)beta dets.
c
            do 2803 igb=iga,itgb
               if (lgcom(igb,jgroae).ne.1) goto 2803
               lab1 = jposae + ldisb(kbsym,igb,jgroae)
               nibs = nbst(igb)
               do 2753 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2753
                  labpos = lab1 + lspb(ibpos)
                  qjmoda = jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2603 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,iga).ne.1) goto 2603
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jpza1+ldisb(ksym,igb2,iga)+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2603          continue
c
 2753          continue
 2803       continue
c
      endif
c
 2900    continue
 3000 continue
c
      else
c
c   ***** direct method below *****
c
      do 4000 isae=1,nsym
         kbsym = ktab(isae)
c
c  analyse excited a' groups for compatibility with b.
c
      ngrps = 0
      do 1976 ii=1,immc(isae)
         icgr = igroa(ii,isae)
         do jj=1,ngrps
            if (jb1gr(jj).eq.icgr) goto 1976
         enddo
         ngrps = ngrps + 1
         jb1gr(ngrps) = icgr
 1976 continue
c
c first type of betas, ksym -> kbsym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7400 igb=1,itgb
             if (igb.lt.iga) goto 7399
             if (lgcom(igb,iga).ne.1) goto 7399
             lab1 = jpza1 + ldisb(ksym,igb,iga)
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
             do 7300 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                if (ibpos.lt.kka) goto 7300
                labpos = lab1 + lspb(ibpos)
                qjmoda=1.0d+00
                if (ibpos.eq.kka) qjmoda=0.5d+00
c
c loop over single beta excites from ibpos which are of symmetry kbsym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7290 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7290
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7285 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7280 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7270
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7270
c
c  get group number
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          nias = nbst(igb2)
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb2,iagrp).eq.1) goto 7555
             enddo
             goto 7270
 7555        continue
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7260 igap=igba,igbe+1
c
                     do 7250 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7250
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)*qjmoda
            jb1p = lspb(iposb + nias)
c
            do 3900 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb2,jgroae).ne.1) goto 3900
c
               jposae=iposa(jsae,isae)
               jcib = jb1p + jposae + ldisb(kbsym,igb2,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
            ab(labpos) = ab(labpos) + c*ci(jcib)
            ab(jcib) = ab(jcib) + c*ci(labpos)
c
 3900       continue
c
 7250                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7260             continue
c
 7270              lbox3(ispb2) = lbox3(ispb2) - 1
 7280           continue
c
 7285           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7290           continue
c
 7300        continue
c
 7399        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7400     continue
c
c
c second type of betas, kbsym -> ksym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7800 igb=1,itgb
             if (igb.lt.iga) goto 7799
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb,iagrp).eq.1) goto 7655
             enddo
             goto 7799
 7655        continue
c
c insert check for a'b group compatibility.
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
c
             do 7700 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                if (ibpos.lt.kka) goto 7700
                lab1 = lspb(ibpos)
                qjmoda=1.0d+00
                if (ibpos.eq.kka) qjmoda=0.5d+00
c
c
c loop over single beta excites from ibpos which are of symmetry ksym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7690 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7690
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7685 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7680 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7670
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7670
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          if (lgcom(igb2,iga).ne.1) goto 7670
          jcibs = jpza1 + ldisb(ksym,igb2,iga)
          nias = nbst(igb2)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7660 igap=igba,igbe+1
c
                     do 7650 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7650
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)*qjmoda
            jcib= jcibs + lspb(iposb + nias)
c
            do 4300 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb,jgroae).ne.1) goto 4300
c
               jposae=iposa(jsae,isae)
               labpos = lab1 + jposae + ldisb(kbsym,igb,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
            ab(labpos) = ab(labpos) + c*ci(jcib)
            ab(jcib) = ab(jcib) + c*ci(labpos)
c
 4300       continue
c
 7650                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7660             continue
c
 7670              lbox3(ispb2) = lbox3(ispb2) - 1
 7680           continue
c
 7685           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7690           continue
c
 7700        continue
c
 7799        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7800     continue
c
 4000 continue
c
      endif
c
c ****** end of direct method *******
c
c  **** end of general case of more than one space ******
c
      goto 4898
c
 3400 continue
c
c ***** special case of one space !!!!! ******
c
c  loop over single alpha excites within each symmetry.
c
       do 2901 isae=1,nsym
          kbsym=ktab(isae)
          do 2801 jsae=1,immc(isae)
                jposae=iposa(jsae,isae)
                jperae=ipera(jsae,isae)
                jindae=iind1(jsae,isae)
                jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym) then
c
c  loop over all relevant betas
c
             labpos=jpza1
             labpos2=jposae
                do 2721 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    labpos2 = labpos2 + 1
                    ibpos =  lsbc(kkb)
                    if (ibpos.lt.ilima) goto 2721
                    qjmoda=jperae
                    if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2621 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib = jposae+jb1po(jbindx)
                  jcib2 = jpza1+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
                  ab(labpos2) = ab(labpos2) + c*ci(jcib2)
                  ab(jcib2) = ab(jcib2) + c*ci(labpos2)
c
 2621          continue
               qjmoda=jperae
c
 2721          continue
c
      else
c
c  loop over all relevant (a-)beta dets.
c
                labpos = jpza1
                do 2751 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    ibpos =  lsbc(kkb)
                    if (ibpos.lt.ilima) goto 2751
                    qjmoda=jperae
                    if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2601 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2601          continue
               qjmoda=jperae
c
 2751          continue
c
c  loop over all relevant (a'-)beta dets.
c
               labpos = jposae
               do 2761 kkb=lsbs(kbsym,1),lsbs(kbsym+1,1)-1
                  labpos = labpos + 1
                  ibpos = lsbc(kkb)
                  if (ibpos.lt.ilima) goto 2761
                  qjmoda=jperae
                  if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2611 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jpza1+jb1po(jbindx)
c
                  ab(labpos) = ab(labpos) + c*ci(jcib)
                  ab(jcib) = ab(jcib) + c*ci(labpos)
c
 2611          continue
               qjmoda=jperae
c
 2761          continue
c
      endif
c
 2801       continue
 2901   continue
c
c  **** end of special case of one space *******
c
 4898       continue
            if (goparr) call ddi_dlbnext(my_task)
 4899       continue
            call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_gsumf(2900,ab,nci)
      endif
c
c **
c
      is = (-1)**ispin(ihmcon(1))
      do 1111 iga=1,itga
         do 1122 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            ilima = kka-nast(iga)
c
            do 1133 igb=1,iga
            if (lgcom(igb,iga).ne.1) goto 1133
            icc1 = jpza1 + ldisb(ksym,igb,iga)
            icc2 = ldisb(jasym,iga,igb) + lspb(kka)
            nibs=nbst(igb)
            do 1144 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               ibpos = lsbc(kkb) + nibs
               if (ibpos.gt.kka) goto 1122
               ici2 = lspa(ibpos) + icc2
               if (ibpos.eq.kka) then
                  ab(ici2) = ab(ici2) + is*ab(ici2)
                  goto 1122
               endif
               ici1 = icc1+lspb(ibpos)
               qt = ab(ici1)
               ab(ici1) = ab(ici1) + is*ab(ici2)
               ab(ici2) = ab(ici2) + is*qt
 1144       continue
 1133       continue
 1122    continue
 1111 continue
c
c   now for the diagonal contributions
c
      do 119 ijk=1,nci
         ab(ijk) = ab(ijk) + q(ijk)*ci(ijk)
  119 continue
c
c  all done, return
c
      return
      end
c
c*module ormas2  *deck fchc0ys
c     ------------------------------------------------------------------
      subroutine fchc0ys(si1,si2,index,nact,na,nb,ci,ab,nv,q,nci,x,nx,
     *           iacon1,iacon2,ibcon1,ibcon2,
     *           lbox1,lbox2,lbox3,lbox4,lbox5,
     *           nsym,iob,lgmul,ktab,
     *           landet,lbndet,nast,nbst,lsyma,lgcom,
     *           lspa,lspb,ldisb,lsbs,lsbc,
     *           itga,itgb,iast,ibst,
     *           iposa,ipera,iind1,igroa,immc,
     *           nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,
     *           ispin,ihmcon)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical fdirct,qcorr
      logical goparr,dskwrk,maswrk
c
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension si1(*),si2(*)
      dimension index((nact*(nact+1))/2+1)
      dimension ci(nci,nv),ab(nci,nv),q(nci),x(nx)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension iob(nact),lgmul(nsym,nsym),ktab(nsym)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsbs(nsym+1,itgb)
      dimension lsbc(ibst)
      dimension iposa(na*(nact-na),nsym)
      dimension ipera(na*(nact-na),nsym)
      dimension iind1(na*(nact-na),nsym)
      dimension igroa(na*(nact-na),nsym)
      dimension immc(nsym)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(nsym+1,ibst+1)
      dimension ispin(*),ihmcon(1)
c
      do ii=1,nci
         do jj=1,nv
           ab(ii,jj) = 0.0d+00
         enddo
      enddo
c
c  --- big loop over all alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_dlbnext(my_task)
      endif
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            if (goparr.and.kka.ne.my_task) goto 4899
            ilima = kka-nast(iga)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            do ii=1,nsym
               immc(ii)=0
            enddo
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) goto 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) goto 4870
c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
c
                     ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
              if (lbox2(ispa2).gt.iama(ispa2)) goto 4800
c
c  get group number
c
           call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
           nias = nast(iga2)
c
              call idpost(iacon2,na,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kasym = lsyma(kapos)
              kper1 = (-1)**jpera
              immc(kasym) = immc(kasym) + 1
              jspo = immc(kasym)
              iposa(jspo,kasym) = jpza2
              ipera(jspo,kasym) = kper1
              iind1(jspo,kasym) = ind
              igroa(jspo,kasym) = iga2
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 4800
c
c   determine the alpha string contribution to the matrix element.
c
              c = si1(ind)
c
              do 4712 ik=1,na
                 if (ik.eq.ia1) goto 4712
                 ion = iacon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
                 c = c + si2(jj1) - si2(inx)
 4712         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4690
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
c
c matrix element addition
c
                 d = c
                 do 4670 ik=1,nb
                    ion = ibcon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    d = d + si2(jj1)
 4670            continue
c
                 t = d*kper1
c
                 do 44 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci2,kj)
                 ab(jci2,kj) = ab(jci2,kj) + t*ci(jci1,kj)
   44            continue
c
 4680         continue
c
 4690         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4700         continue
c
c --  double alpha excitations start here  ---
c
 4800         continue
c
            if (ia1.eq.na) goto 4850
            if (jj.eq.nact) goto 4850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            iaes3 = 1
            do kk=1,ispa1-1
               iaes3 = iaes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispa1.
c
            do 3890 ispa3 = ispa1,nspace
               if (lbox1(ispa3).eq.0) goto 3887
               if (ispa3.eq.ispa1.and.ia1.eq.ieae) goto 3887
               ioc3 = lbox2(ispa3)
               if (ioc3.eq.0) goto 3887
c
               iaee3 = iaes3 + lbox2(ispa3)-1
               lbox3(ispa3) = lbox3(ispa3)-1
c
c  loop over electrons in ispa3, which are larger than ia1,
c  making sure it isn't the already excited electron.
c
               jsta3 = iaes3
               if (ispa3.eq.ispa1) then
                  jsta3=ia1
                  if (igap-1.eq.ia1) jsta3=jsta3+1
               endif
c
               do 3880 ia3=jsta3,iaee3
                  if (ia3.eq.igap-1) goto 3880
                  io3 = iacon2(ia3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igae3 = 0
                  do jik=1,ispa2-1
                     igae3 = igae3 + lbox2(jik)
                  enddo
                  do 3850 ispa4=ispa2,nspace
c
c  igas3, igae3 are electrons specifying ispa4 space electron limits.
c
                  lbox3(ispa4) = lbox3(ispa4) + 1
                  igas3 = igae3 + 1
                  igae3 = igae3 + lbox2(ispa4)
                  if (lbox3(ispa3).lt.iami(ispa3)) goto 3840
                  if (lbox3(ispa4).gt.iama(ispa4)) goto 3840
             if (ispa4.eq.ispa2.and.jj.eq.msta(ispa2+1)-1) goto 3840
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox3,iga3)
               nias3 = nast(iga3)
c
c  make gap information here.
c
                  igaa3 = max(igas3,igap)
c
                  if (lbox2(ispa4).eq.0) then
                     ista3 = msta(ispa4)
                     iend3 = msta(ispa4+1)-1
                  elseif (ispa4.eq.ispa2) then
                     ista3 = jj+1
                     iend3 = iacon2(igaa3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igae3) iend3=msta(ispa2+1)-1
c
                  else
                     ista3 = msta(ispa4)
                     iend3 = iacon2(igaa3)-1
                  endif
c
c  loop over gaps
c
                  do 3830 igap3=igaa3,igae3+1
c
                     do 3820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 3820
c
              call rede00(iacon2,ibcon1,na,ia3,igap3-1,jj3,jpera3)
              call idpost(ibcon1,na,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga3),ibcon2,jposa3)
              iper3 = (-1)**(jpera3+jpera)
              kapos3 = jposa3 + nias3
              jpza3 = lspa(kapos3)
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
                 c = si2(inx) - si2(inx2)
                 t = c*iper3
c
c  loop over beta strings of right symmetry and group.
c
              do 3700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga3).ne.1) goto 3700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci3 = jpza3 + ldisb(ksym,igb,iga3)
c
              do 3680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci3 = jci3 + 1
c
                 do 55 kj=1,nv
                 ab(jci1,kj) = ab(jci1,kj) + t*ci(jci3,kj)
                 ab(jci3,kj) = ab(jci3,kj) + t*ci(jci1,kj)
   55            continue
c
 3680         continue
c
 3700         continue
c
 3820                continue
c
                     ista3 = iacon2(igap3)+1
                     iend3 = iacon2(igap3+1)-1
                     if (igap3.eq.igae3) iend3=msta(ispa4+1)-1
 3830             continue
c
 3840             lbox3(ispa4) = lbox3(ispa4) - 1
 3850             continue
c
 3880          continue
c
               lbox3(ispa3) = lbox3(ispa3)+1
 3887          iaes3 = iaes3 + lbox3(ispa3)
 3890       continue
c
 4850                continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c
c  --- end of loop over single alpha excitations. ---
c      now to sort them by positions within symmetries.
c
            do ii=1,nsym
               call fccsrt3(igroa(1,ii),ipera(1,ii),iind1(1,ii),
     *                   iposa(1,ii),immc(ii))
            enddo
c
c  --- end of loop over pure alpha excitations.
c  now to loop over all simultaneous ab -> a'b' excitations.
c
      if (nspace.eq.1) goto 3400
c
c  ***** general case of more than one space !!!!! *******
c
      if (.not.fdirct) then
c
c  loop over single alpha excites within each symmetry.
c
      do 3000 isae=1,nsym
         kbsym=ktab(isae)
         do 2900 jsae=1,immc(isae)
            jposae=iposa(jsae,isae)
            jperae=ipera(jsae,isae)
            jindae=iind1(jsae,isae)
            jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym.and.iga.eq.jgroae) then
c
c  loop over all relevant betas
c
            do 2813 igb=iga,itgb
               if (lgcom(igb,iga).ne.1) goto 2813
               lab1 = jpza1 + ldisb(ksym,igb,iga)
               lab2 = jposae + ldisb(ksym,igb,iga)
               nibs = nbst(igb)
               do 2763 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2763
                  labpos = lab1 + lspb(ibpos)
                  labpos2 = lab2 + lspb(ibpos)
                  qjmoda=jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2613 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2613
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  jcib2=jpza1+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  do 66 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
                  ab(labpos2,kj) = ab(labpos2,kj) + c*ci(jcib2,kj)
                  ab(jcib2,kj) = ab(jcib2,kj) + c*ci(labpos2,kj)
   66             continue
c
 2613          continue
c
 2763          continue
 2813       continue
c
       else
c
c  loop over all relevant (a-)beta dets.
c
            do 2800 igb=iga,itgb
               if (lgcom(igb,iga).ne.1) goto 2800
               lab1 = jpza1 + ldisb(ksym,igb,iga)
               nibs = nbst(igb)
               do 2750 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2750
                  labpos = lab1 + lspb(ibpos)
                  qjmoda = jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2600 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2600
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
c
                  do 77 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   77             continue
c
 2600          continue
c
 2750          continue
 2800       continue
c
c
c  loop over all relevant (a'-)beta dets.
c
            do 2803 igb=iga,itgb
               if (lgcom(igb,jgroae).ne.1) goto 2803
               lab1 = jposae + ldisb(kbsym,igb,jgroae)
               nibs = nbst(igb)
               do 2753 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                  ibpos = nibs + lsbc(kkb)
                  if (ibpos.lt.kka) goto 2753
                  labpos = lab1 + lspb(ibpos)
                  qjmoda = jperae
                  if (ibpos.eq.kka) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2603 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,iga).ne.1) goto 2603
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jpza1+ldisb(ksym,igb2,iga)+jb1po(jbindx)
c
                  do 88 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   88             continue
c
 2603          continue
c
 2753          continue
 2803       continue
c
      endif
c
 2900    continue
 3000 continue
c
      else
c
c  ***** direct method below *****
c
      do 4000 isae=1,nsym
         kbsym = ktab(isae)
c
c  analyse excited a' groups for compatibility with b.
c
      ngrps = 0
      do 1976 ii=1,immc(isae)
         icgr = igroa(ii,isae)
         do jj=1,ngrps
            if (jb1gr(jj).eq.icgr) goto 1976
         enddo
         ngrps = ngrps + 1
         jb1gr(ngrps) = icgr
 1976 continue
c
c first type of betas, ksym -> kbsym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7400 igb=1,itgb
             if (igb.lt.iga) goto 7399
             if (lgcom(igb,iga).ne.1) goto 7399
             lab1 = jpza1 + ldisb(ksym,igb,iga)
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
             do 7300 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                if (ibpos.lt.kka) goto 7300
                labpos = lab1 + lspb(ibpos)
                qjmoda=1.0d+00
                if (ibpos.eq.kka) qjmoda=0.5d+00
c
c loop over single beta excites from ibpos which are of symmetry kbsym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7290 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7290
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7285 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7280 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7270
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7270
c
c  get group number
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          nias = nbst(igb2)
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb2,iagrp).eq.1) goto 7555
             enddo
             goto 7270
 7555        continue
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7260 igap=igba,igbe+1
c
                     do 7250 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7250
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)*qjmoda
            jb1p = lspb(iposb + nias)
c
            do 3900 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb2,jgroae).ne.1) goto 3900
c
               jposae=iposa(jsae,isae)
               jcib = jb1p + jposae + ldisb(kbsym,igb2,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
         do kj=1,nv
            ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
            ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
         enddo
c
 3900       continue
c
 7250                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7260             continue
c
 7270              lbox3(ispb2) = lbox3(ispb2) - 1
 7280           continue
c
 7285           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7290           continue
c
 7300        continue
c
 7399        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7400     continue
c
c second type of betas, kbsym -> ksym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7800 igb=1,itgb
             if (igb.lt.iga) goto 7799
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb,iagrp).eq.1) goto 7655
             enddo
             goto 7799
 7655        continue
c
c insert check for a'b group compatibility.
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
c
             do 7700 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                if (ibpos.lt.kka) goto 7700
                lab1 = lspb(ibpos)
                qjmoda=1.0d+00
                if (ibpos.eq.kka) qjmoda=0.5d+00
c
c
c loop over single beta excites from ibpos which are of symmetry ksym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7690 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7690
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7685 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7680 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7670
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7670
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          if (lgcom(igb2,iga).ne.1) goto 7670
          jcibs = jpza1 + ldisb(ksym,igb2,iga)
          nias = nbst(igb2)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7660 igap=igba,igbe+1
c
                     do 7650 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7650
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)*qjmoda
            jcib= jcibs + lspb(iposb + nias)
c
            do 4300 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb,jgroae).ne.1) goto 4300
c
               jposae=iposa(jsae,isae)
               labpos = lab1 + jposae + ldisb(kbsym,igb,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               c = si2(ix)*qjper*jperae
c
         do kj=1,nv
            ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
            ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
         enddo
c
 4300       continue
c
 7650                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7660             continue
c
 7670              lbox3(ispb2) = lbox3(ispb2) - 1
 7680           continue
c
 7685           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7690           continue
c
 7700        continue
c
 7799        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7800     continue
c
 4000 continue
c
      endif
c
c ****** end of direct method *******
c
c  **** end of general case of more than one space ******
c
      goto 4898
c
 3400 continue
c
c ***** special case of one space !!!!! ******
c
c  loop over single alpha excites within each symmetry.
c
       do 2901 isae=1,nsym
          kbsym=ktab(isae)
          do 2801 jsae=1,immc(isae)
                jposae=iposa(jsae,isae)
                jperae=ipera(jsae,isae)
                jindae=iind1(jsae,isae)
                jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym) then
c
c  loop over all relevant betas
c
             labpos=jpza1
             labpos2=jposae
                do 2721 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    labpos2 = labpos2 + 1
                    ibpos =  lsbc(kkb)
                    if (ibpos.lt.ilima) goto 2721
                    qjmoda=jperae
                    if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2621 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib = jposae+jb1po(jbindx)
                  jcib2 = jpza1+jb1po(jbindx)
c
                  do 56 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
                  ab(labpos2,kj) = ab(labpos2,kj) + c*ci(jcib2,kj)
                  ab(jcib2,kj) = ab(jcib2,kj) + c*ci(labpos2,kj)
   56             continue
c
 2621          continue
               qjmoda=jperae
c
 2721          continue
c
      else
c
c  loop over all relevant (a-)beta dets.
c
                labpos = jpza1
                do 2751 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    ibpos =  lsbc(kkb)
                    if (ibpos.lt.ilima) goto 2751
                    qjmoda=jperae
                    if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2601 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jposae+jb1po(jbindx)
c
                  do 67 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   67             continue
c
 2601          continue
               qjmoda=jperae
c
 2751          continue
c
c  loop over all relevant (a'-)beta dets.
c
               labpos = jposae
               do 2761 kkb=lsbs(kbsym,1),lsbs(kbsym+1,1)-1
                  labpos = labpos + 1
                  ibpos = lsbc(kkb)
                  if (ibpos.lt.ilima) goto 2761
                  qjmoda=jperae
                  if (ibpos.eq.ilima) qjmoda=jperae/2.0d+00
c
c  loop over single beta excites from ibpos
c
               do 2611 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  c = si2(ix)*qjmoda*jb1pe(jbindx)
                  jcib=jpza1+jb1po(jbindx)
c
                  do 78 kj=1,nv
                  ab(labpos,kj) = ab(labpos,kj) + c*ci(jcib,kj)
                  ab(jcib,kj) = ab(jcib,kj) + c*ci(labpos,kj)
   78             continue
c
 2611          continue
               qjmoda=jperae
c
 2761          continue
c
      endif
c
 2801       continue
 2901   continue
c
c  **** end of special case of one space *******
c
 4898       continue
            if (goparr) call ddi_dlbnext(my_task)
 4899       continue
            call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings. ---
c
      if (goparr) then
         call ddi_dlbreset()
         call ddi_gsumf(2900,ab,nci*nv)
      endif
c
c **
c
      do 1111 iga=1,itga
         do 1122 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            ilima = kka-nast(iga)
c
            do 1133 igb=1,iga
            if (lgcom(igb,iga).ne.1) goto 1133
            icc1 = jpza1 + ldisb(ksym,igb,iga)
            icc2 = ldisb(jasym,iga,igb) + lspb(kka)
            nibs=nbst(igb)
            do 1144 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               ibpos = lsbc(kkb) + nibs
               if (ibpos.gt.kka) goto 1122
               ici2 = lspa(ibpos) + icc2
               if (ibpos.eq.kka) then
                  do 4444 kj=1,nv
                  is = (-1)**ispin(ihmcon(kj))
                  ab(ici2,kj) = ab(ici2,kj) + is*ab(ici2,kj)
 4444             continue
                  goto 1122
               endif
               ici1 = icc1+lspb(ibpos)
               do 3333 kj=1,nv
               is = (-1)**ispin(ihmcon(kj))
               qt = ab(ici1,kj)
               ab(ici1,kj) = ab(ici1,kj) + is*ab(ici2,kj)
               ab(ici2,kj) = ab(ici2,kj) + is*qt
 3333          continue
 1144       continue
 1133       continue
 1122    continue
 1111 continue
c
c   now for the diagonal contributions
c
      do 119 ijk=1,nci
         do 118 kj=1,nv
            ab(ijk,kj) = ab(ijk,kj) + q(ijk)*ci(ijk,kj)
  118    continue
  119 continue
c
c  all done, return
c
      return
      end
c*module ormas2  *deck masprt
      subroutine masprt(iw,some,vec,nast,itga,itgb,
     *             lsyma,iast,ibst,lgcom,lsbs,nsym,lsbc,
     *             lbox1,lbox2,lbox3,iacon1,ibcon1,ktab)
c
      implicit double precision(a-h,o-z)
c
      logical some,goparr,dskwrk,maswrk,svdskw,fdirct,qcorr
c
      parameter (mxrt=100)
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
      common /enrgys/ enucr,eelct,etot,stot,ssquar,ecore,escf,eerd,
     *                e1,e2,ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
c
      dimension vec(nci,nstate)
      dimension nast(itga),lsyma(iast)
      dimension lgcom(itgb,itga),lsbs(nsym+1,itgb),lsbc(ibst)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension iacon1(na),ibcon1(nb),ktab(nsym)
      character*250 name
c
      data check/8hcheck   /
c
c     ----- print the determinant based ci eigenvector -----
c
      svdskw = dskwrk
      dskwrk=.false.
c
      nsym = 2**igpdet
c
      if(some) write(iw,9140) grpdet
c
c        ----- print ci energies and eigenvectors -----
c        note that prtdet destroys the eigenvectors.
c
      if(nci.le.20) then
         iop=1
         numprt=nci
         if(some) write(iw,9150)
      else
         iop=2
         numprt=0
         if(some) write(iw,9160) prttol
      end if
c
       iunit88=88
      do 430 i=1,nstate
       write(name,'(A13,I3.3)')'civecs.ascii_',i
       if (iw.eq.6)open(unit=iunit88,file=name,form='formatted')
         if(some) then
            write(iw,9170) i,estate(i),spins(i),sz,stsym
            if(exetyp.ne.check) call maspri(iw,iop,numprt,
     *      vec(1,i),nast,itga,itgb,lsyma,iast,ibst,lgcom,
     *      lsbs,nsym,lsbc,lbox1,lbox2,lbox3,iacon1,ibcon1,ktab,iunit88)
c
         end if
        if (iw.eq.6)close(iunit88)
  430 continue
cc
      dskwrk=svdskw
      return
c
 9140 format(/1x,'ci eigenvectors will be labeled in group=',a8)
 9150 format(1x,'printing all non-zero ci coefficients')
 9160 format(1x,'printing ci coefficients larger than',f10.6)
 9170 format(/1x,'state',i4,'  energy= ',f20.10,'  s=',f6.2,
     *           '  sz=',f6.2,:,'  space sym=',a4/)
      end
c
c*module ormas2  *deck maspri
c     -----------------------------------------------------------
      subroutine maspri(iw,iop,num,
     *               ci,nast,itga,itgb,lsyma,
     *               iast,ibst,lgcom,lsbs,nsym,lsbc,
     *               lbox1,lbox2,lbox3,iacon1,ibcon1,ktab,iwp)
c     -----------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      parameter (mxrt=100)
c
      logical goparr,dskwrk,maswrk,fdirct,qcorr
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      character*500 cona,conb
      dimension ci(nci)
      dimension nast(itga+1),lsyma(iast)
      dimension lgcom(itgb,itga),lsbs(nsym+1,itgb),lsbc(ibst)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension iacon1(na),ibcon1(nb),ktab(nsym)
c
      do ii=1,500
         cona(ii:ii) = ' '
         conb(ii:ii) = ' '
      enddo
c
c  set up the table
c
      isiza = 0
      isizb = 0
      do ii=1,nspace
         isiza = isiza+(3*iama(ii))
         isizb = isizb+(3*ibma(ii))
      enddo
c
      if(isiza+1.gt.500  .or.  isizb+1.gt.500) then
         write(iw,*) 'too many orbitals to print ormas ci vector'
         return
      end if
c
      iapl = isiza/2 - 2
      ibpl = isizb/2 - 1
      cona(iapl:iapl+4) = 'alpha'
      conb(ibpl:ibpl+4) = 'beta '
      if(maswrk) write(iw,'(4a)') cona(1:isiza+1),'|',
     *                            conb(1:isizb+1),'| coefficient'
c
      do ii=1,isiza+1
         cona(ii:ii) = ' '
      enddo
      do ii=1,isizb+1
         conb(ii:ii) = ' '
      enddo
      iapl = 1
      ibpl = 1
      do ii=1,nspace
         if (ii.lt.10) write(cona(iapl+1:iapl+1),'(i1)') ii
         if (ii.lt.10) write(conb(ibpl+1:ibpl+1),'(i1)') ii
         if (ii.gt.10) write(cona(iapl+1:iapl+2),'(i2)') ii
         if (ii.gt.10) write(conb(ibpl+1:ibpl+2),'(i2)') ii
         iapl = iapl + iama(ii)*3
         ibpl = ibpl + ibma(ii)*3
      enddo
      if (maswrk) write(iw,'(4a)') cona(1:isiza+1),'|',
     *                  conb(1:isizb+1),'|'
      do ii=1,isiza+1
         cona(ii:ii) = '-'
      enddo
      do ii=1,isizb+1
         conb(ii:ii) = '-'
      enddo
      if(maswrk) write(iw,'(4a)') cona(1:isiza+1),'|',
     *                 conb(1:isizb+1),'|------------'
c
      if (iop.eq.1) then
c
      do 3000 kjk=1,num
c
         ici = 0
         ipos = -1
         pmax = 0.0d+00
c
         do 400 iga=1,itga
            do 410 kka=nast(iga)+1,nast(iga+1)
               jasym=lsyma(kka)
               ksym=ktab(jasym)
               do 500 igb=1,itgb
                  if (lgcom(igb,iga).ne.1) go to 500
                  do 510 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                     nend = lsbc(kkb)
                     ici=ici+1
                     if (abs(ci(ici)).gt.pmax) then
                        iza = iga
                        izb = igb
                        ika = kka - nast(iga)
                        ikb = nend
                        ipos = ici
                        pmax = abs(ci(ici))
                     endif
  510             continue
  500          continue
  410       continue
  400    continue
         if (ipos.eq.-1) go to 3000
c
c  make the determinant
c
         call resetco(lbox1,nspace,na,iama,iami,lbox3)
         do ii=1,iza-1
            call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
         enddo
         call resetde(lbox1,nspace,na,msta,iacon1)
         do ii=1,ika-1
            call moveup2(lbox1,nspace,na,msta,iacon1)
         enddo
c
         call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
         do ii=1,izb-1
            call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
         enddo
         call resetde(lbox2,nspace,nb,msta,ibcon1)
         do ii=1,ikb-1
            call moveup2(lbox2,nspace,nb,msta,ibcon1)
         enddo
c
c   now to print out the determinant
c
         do ii=1,500
            cona(ii:ii) = ' '
            conb(ii:ii) = ' '
         enddo
         ispa = 1
         ispb = 1
         ipla = 1
         iplb = 1
         iela = 1
         ielb = 1
         do ii=1,nspace
            ipla = ispa
            do jj=1,lbox1(ii)
               write(cona(ipla:ipla+2),'(i3)') iacon1(iela)
               ipla = ipla+3
               iela = iela+1
            enddo
            ispa = ispa + 3*iama(ii)
         enddo
         do ii=1,nspace
            iplb = ispb
            do jj=1,lbox2(ii)
               write(conb(iplb:iplb+2),'(i3)') ibcon1(ielb)
               iplb = iplb+3
               ielb = ielb + 1
            enddo
            ispb = ispb + 3*ibma(ii)
         enddo
         if (maswrk) write(iw,'(4a,f20.15)') cona(1:isiza+1),'|',
     *                  conb(1:isizb+1),'|  ',ci(ipos)
         if (maswrk) write(iwp,'(4a,f20.15)') cona(1:isiza+1),' ',
     *                  conb(1:isizb+1),'   ',ci(ipos)
         ci(ipos) = 0.0d+00
c
 3000 continue
c
      else
c
      do 4000 kjk=1,nci
c
         ici = 0
         ipos = -1
         pmax = 0.0d+00
c
         do 700 iga=1,itga
            do 710 kka=nast(iga)+1,nast(iga+1)
               jasym=lsyma(kka)
               ksym=ktab(jasym)
               do 800 igb=1,itgb
                  if (lgcom(igb,iga).ne.1) go to 800
                  do 810 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                     nend = lsbc(kkb)
                     ici=ici+1
                     if (abs(ci(ici)).gt.pmax) then
                        iza = iga
                        izb = igb
                        ika = kka - nast(iga)
                        ikb = nend
                        ipos = ici
                        pmax = abs(ci(ici))
                     endif
  810             continue
  800          continue
  710       continue
  700    continue
c
c  check if is bigger than crit
c
c        if (abs(ci(ipos)).lt.prttol) return
c
c  make the determinant
c
         call resetco(lbox1,nspace,na,iama,iami,lbox3)
         do ii=1,iza-1
            call pushco(lbox1,nspace,na,iama,iami,lbox3,ise)
         enddo
         call resetde(lbox1,nspace,na,msta,iacon1)
         do ii=1,ika-1
            call moveup2(lbox1,nspace,na,msta,iacon1)
         enddo
c
         call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
         do ii=1,izb-1
            call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,ise)
         enddo
         call resetde(lbox2,nspace,nb,msta,ibcon1)
         do ii=1,ikb-1
            call moveup2(lbox2,nspace,nb,msta,ibcon1)
         enddo
c
c   now to print out the determinant
c
         do ii=1,500
            cona(ii:ii) = ' '
            conb(ii:ii) = ' '
         enddo
         ispa = 1
         ispb = 1
         ipla = 1
         iplb = 1
         iela = 1
         ielb = 1
         do ii=1,nspace
            ipla = ispa
            do jj=1,lbox1(ii)
               write(cona(ipla:ipla+2),'(i3)') iacon1(iela)
               ipla = ipla+3
               iela = iela+1
            enddo
            ispa = ispa + 3*iama(ii)
         enddo
         do ii=1,nspace
            iplb = ispb
            do jj=1,lbox2(ii)
               write(conb(iplb:iplb+2),'(i3)') ibcon1(ielb)
               iplb = iplb+3
               ielb = ielb + 1
            enddo
            ispb = ispb + 3*ibma(ii)
         enddo
         if (abs(ci(ipos)).gt.prttol) then
         if (maswrk) write(iw,'(4a,f20.15)') cona(1:isiza+1),'|',
     *                  conb(1:isizb+1),'|  ',ci(ipos)
         endif
         if (maswrk) write(iwp,'(4a,f20.15)') cona(1:isiza+1),' ',
     *                  conb(1:isizb+1),'   ',ci(ipos)
         ci(ipos) = 0.0d+00
c
 4000 continue
c
      endif
c
      return
      end
c
c*module ormas2  *deck cicopy
c     -----------------------------------------------------------
      subroutine cicopy(ci,ab,ntot)
c     -----------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ci(ntot),ab(ntot)
      do ii=1,ntot
         ab(ii) = ci(ii)
      enddo
      return
      end
c
c*module ormas2  *deck masdm1
c     ------------------------------------------------------------------
      subroutine masdm1(iw,nprint,
     *           dm1,m2,nact,nci,na,nb,iroot,x,nx,ci,
     *           index,nsym,iob,
     *           lbox1,lbox2,lbox4,lbox5,
     *           ktab,iacon1,iacon2,ibcon1,ibcon2,
     *           landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *           lspa,lspb,ldisb,lsas,lsbs,lsac,
     *           itga,itgb,iast,ibst)
c     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical some,goparr,dskwrk,maswrk,fdirct,qcorr
c
      parameter (mxrt=100)
c
      common /enrgys/ enucr,eelct,etot,stot,ssquar,ecore,escf,eerd,
     *                e1,e2,ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
c
      dimension dm1(m2),x(nx),ci(nci)
      dimension index((nact*(nact+1))/2+1),iob(nact)
      dimension lbox1(nspace),lbox2(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension ktab(nsym)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast)
c
c     flag for 1st masscf iter
c
      logical first_fnr
      common /flag_fnr/ first_fnr
c
      some = maswrk  .and.  nprint.ne.-5 .and. first_fnr
c
      if (some) write(iw,9000) iroot
      if (some) write(iw,9010) iroot,estate(iroot)
c
c  generate the one particle density matrix for each state.
c
      do 13 i=1,m2
         dm1(i) = 0.0d+00
   13 continue
c
c ---  big loop over all alpha determinants ---
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) go to 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) go to 4870
                  if (lbox2(ispa2).gt.iama(ispa2)) go to 4870
c
c  get group number
c
           call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
           nias = nast(iga2)
c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip.
c
              if (is1.ne.is2) go to 4850
c
              ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
              call idpost(iacon2,na,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kper1 = (-1)**jpera
c
c  loop over beta strings of right symmetry and group
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
                 fc = ci(jci1)*ci(jci2)*kper1
                 dm1(ind) = dm1(ind) + fc
 4680         continue
c
 4700         continue
c
 4850          continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c   diagonal contributions here.
c
            do 67 ii=1,na
               i1 = iacon1(ii)
               ind1 = index(i1+1)
c
               do 53 igb=1,itgb
                  if (lgcom(igb,iga).ne.1) go to 53
                  jci1 = jpza1 + ldisb(ksym,igb,iga)
                  do 58 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                     jci1 = jci1 + 1
                     fc = ci(jci1)*ci(jci1)
                     dm1(ind1) = dm1(ind1) + fc
   58             continue
   53          continue
c
   67       continue
c
            call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings ----
c
c  --- big loop over beta -----
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox5)
c
      do 8000 igb=1,itgb
c
         call resetde(lbox1,nspace,nb,msta,ibcon1)
c
c  kkb gives the actual position of the beta string ibcon1 in
c  the full beta string list.
c
         do 7900 kkb=nbst(igb)+1,nbst(igb+1)
            jpzb1 = lspb(kkb)
            kbsym = lsymb(kkb)
            ksym=ktab(kbsym)
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            iebs = nb+1
            do 7890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) go to 7890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1.
c  iebs, iebe are the electrons in space ispb1.
c
               do 7885 ib1=iebe,iebs,-1
                  io1 = ibcon1(ib1)
                  igbe = iebe - lbox1(ispb1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 7880 ispb2=ispb1,nspace
c
c  igbs, igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) go to 7870
                  if (lbox2(ispb2).gt.ibma(ispb2)) go to 7870
c
c  get group number
c
           call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox2,igb2)
           nibs = nbst(igb2)
c
c  make gap information here.
c
                  igbb = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igbb)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igbb)-1
                  endif
c
c  loop over gaps
c
                  do 7860 igap=igbb,igbe+1
c
                     do 7850 jj=ista,iend
                        is2 = iob(jj)
                        ind = index(jj) + io1
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
                        if (is1.ne.is2) go to 7850
c
              call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,jperb)
c
              call idpost(ibcon2,nb,lbox2,nspace,msta,idim,x,nx,lbst,
     *                  lbndet(1,igb2),iacon1,jposb)
              kbpos = jposb + nibs
              jpzb2 = lspb(kbpos)
              kper1 = (-1)**jperb
c
c  loop over alpha strings of the right group and symmetry.
c
              do 7700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7700
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              do 7680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
c
                 fc = ci(jci1)*ci(jci2)*kper1
                 dm1(ind) = dm1(ind) + fc
c
 7680         continue
 7700         continue
c
 7850                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7860             continue
c
 7870             lbox2(ispb2) = lbox2(ispb2) - 1
 7880          continue
c
 7885          continue
c
               lbox2(ispb1) = lbox2(ispb1)+1
 7890       continue
c
c  remaining part of diagonal contributions.
c
            do 69 ii=1,nb
               i1 = ibcon1(ii)
               ind1 = index(i1+1)
c
              do 76 iga=1,itga
              if (lgcom(igb,iga).ne.1) go to 76
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 nias = nast(iga)
c
              do 78 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
c
                 fc = ci(jci1)*ci(jci1)
                 dm1(ind1) = dm1(ind1) + fc
c
   78         continue
   76         continue
c
   69       continue
c
            call moveup2(lbox1,nspace,nb,msta,ibcon1)
 7900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox5,iend)
 8000 continue
c
      return
c
 9000 format(/5x,27(1h-)/5x,'one particle density matrix'/5x,27(1h-)//
     *  1x,'density matrix will be saved for properties of state',i4)
 9010 format(/1x,'ci eigenstate',i4,' total energy =',f20.10)
c
      end
c
c*module ormas2  *deck masdm2
c     ----------------------------------------------------------------
      subroutine masdm2(iw,nprint,iwts,wstate,spins,ipures,s,k,grpdet,
     *     ncorsv,
     *     dm1,dm2,m2,m4,nact,nci,na,nb,ci,ab,
     *     y,nx,index,nsym,iob,
     *     lbox1,lbox2,lbox3,lbox4,lbox5,
     *     lgmul,ktab,iacon1,iacon2,ibcon1,ibcon2,
     *     landet,lbndet,nast,nbst,lsyma,lsymb,lgcom,
     *     lspa,lspb,ldisb,lsas,lsbs,lsac,lsbc,
     *     itga,itgb,iast,ibst,
     *     iposa,ipera,iind1,igroa,immc,
     *     nb1ex,jb1gr,jb1pe,jb1in,jb1po,jb1st,idim1,idim2,x)
c     ----------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      logical some,goparr,dskwrk,maswrk,fdirct,qcorr
c
      parameter (mxrt=100, mxatm=2000)
      parameter (zero=0.0d+00)
c
      common /enrgys/ enucr,eelct,etot,stot,ssquar,ecore,escf,eerd,
     *                e1,e2,ven,vee,epot,ekin,estate(mxrt),statn,edft(2)
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
      common /funct / e,egrad(3,mxatm)
      common /infoa_f / nat,ich,mul,num,nqmt,ne,ma,mb,
     *                zan(mxatm),c(3,mxatm),ian(mxatm)
      common /output_f/ nprintx,itol,icut,normf,normp,nopk
      common /par   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      common /pcklab/ labsiz
      common /runopt_f/ runtyp,exetyp,nevals,nglevl,nhlevl
c
c     flag for 1st masscf iter
c
      logical first_fnr
      common /flag_fnr/ first_fnr
c
c     gamess-uk memory management, replaces /fmcom/
c
c     common /fmcom / x(1)
      dimension x(*)
c
      dimension iwts(mxrt),wstate(mxrt),spins(mxrt)
      dimension dm1(m2),dm2(m4)
      dimension ci(nci,k),ab(nci,k),y(nx)
      dimension index((nact*(nact+1))/2+1),iob(nact)
      dimension lbox1(nspace),lbox2(nspace),lbox3(nspace)
      dimension lbox4(nspace),lbox5(nspace)
      dimension lgmul(nsym,nsym),ktab(nsym)
      dimension iacon1(na),iacon2(na),ibcon1(na),ibcon2(na)
      dimension landet(nspace,itga),lbndet(nspace,itgb)
      dimension nast(itga+1),nbst(itgb+1)
      dimension lsyma(iast),lsymb(ibst),lgcom(itgb,itga)
      dimension lspa(iast),lspb(ibst),ldisb(nsym,itgb,itga)
      dimension lsas(nsym+1,itga),lsbs(nsym+1,itgb)
      dimension lsac(iast),lsbc(ibst)
      dimension iposa(na*(nact-na),nsym)
      dimension ipera(na*(nact-na),nsym)
      dimension iind1(na*(nact-na),nsym)
      dimension igroa(na*(nact-na),nsym)
      dimension immc(nsym)
      dimension jb1gr(nb1ex),jb1pe(nb1ex),jb1in(nb1ex),jb1po(nb1ex)
      dimension jb1st(idim1,idim2)
c
c     ----  state-averaged 1e- and 2e- density matrix  ----
c
      l1 = num
      m1 = nact
      nocc1 = ncorsv + nact
      nocc2 = (nocc1*nocc1+nocc1)/2
      some = maswrk  .and.  nprint.ne.-5 .and. first_fnr
      if(some) write(iw,9310)
c
      mxstat=0
      mxnzw=0
c
      do 100 i=mxrt,1,-1
         if(wstate(i).gt.zero) then
            if(mxstat.eq.0) mxstat=i
            mxnzw=mxnzw+1
         end if
  100 continue
c
      if(mxstat.eq.0) then
         write(iw,*) 'oops, in -masdm2-, something happened to wstate'
         call abrt
      end if
c
      if (some) write(iw,9320) mxnzw
c
c     set state averaged energy, print root information
c
      e = zero
      nxtr=0
      do 310 ist=1,k
         if(ipures.eq.1  .and.  abs(spins(ist)-s).gt.0.03d+00) go to 310
         nxtr=nxtr+1
         if(wstate(nxtr).gt.zero) then
            e = e + wstate(nxtr) * estate(ist)
            if(some) write(iw,9340) ist,estate(ist),
     *                              wstate(nxtr),spins(ist)
         end if
         if(nxtr.gt.mxstat) go to 320
  310 continue
c
c     croak the job if we didn't calculate enough roots with the
c     desired spin multiplicity during the ci diagonalization.
c     if this happens on the 1st mcscf iter, we've already got
c     the ci expansions printed out, and should not repeat it.
c
      if(nxtr.lt.mxstat) then
         if(maswrk) write(iw,9350) nxtr,s,mxstat
         call abrt
      end if
c
  320 continue
c
c     copy ci coefficients for all states with non-zero weights into ab
c
      nxtw=1
      nxtr=0
      do 620 ist=1,k
         if (ipures.eq.1) then
            if (abs(spins(ist)-s).gt.0.03d+00) go to 620
            nxtr=nxtr+1
         else
            nxtr=ist
         endif
c
         if (nxtr.eq.iwts(nxtw)) then
            call cicopy(ci(1,ist),ab(1,nxtw),nci)
            nxtw = nxtw + 1
         endif
  620 continue
c
      do ii=1,m2
         dm1(ii) = 0.0d+00
      enddo
      do ii=1,m4
         dm2(ii) = 0.0d+00
      enddo
c
      nxtw = nxtw - 1
c
c  ---------  now to determine state averaged density matrices ---------
c
c
c  --- big loop over all alpha strings. ---
c
      call resetco(lbox1,nspace,na,iama,iami,lbox5)
c
      do 5000 iga=1,itga
c
         call resetde(lbox1,nspace,na,msta,iacon1)
c
c  kka gives the actual position of the alpha string iacon1 in
c  the full alpha string list.
c
         do 4900 kka=nast(iga)+1,nast(iga+1)
            jpza1 = lspa(kka)
            jasym = lsyma(kka)
            ksym=ktab(jasym)
            do ii=1,nsym
               immc(ii)=0
            enddo
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            ieas = na+1
            do 4890 ispa1=nspace,1,-1
               ioc1 = lbox1(ispa1)
               ieae = ieas - 1
               ieas = ieas - ioc1
               if (ioc1.eq.0) go to 4890
               lbox2(ispa1) = lbox2(ispa1)-1
c
c  loop electrons in space ispa1.
c  ieas, ieae are the electrons in space ispa1.
c
               do 4885 ia1=ieae,ieas,-1
                  io1 = iacon1(ia1)
                  igae = ieae - lbox1(ispa1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 4880 ispa2=ispa1,nspace
c
c  igas, igae are electrons specifying ispa2 space electron limits.
c
                  igas = igae + 1
                  igae = igae + lbox1(ispa2)
c
                  lbox2(ispa2) = lbox2(ispa2) + 1
                  if (lbox2(ispa1).lt.iami(ispa1)) go to 4870
c
c  make gap information here.
c
                  igaa = max(ia1+1,igas)
                  if (lbox1(ispa2).eq.0) then
                     ista = msta(ispa2)
                     iend = msta(ispa2+1)-1
                  elseif (ispa2.eq.ispa1) then
                     ista = io1+1
                     iend = iacon1(igaa)-1
                     if (ia1.eq.ieae) iend=msta(ispa1+1)-1
                  else
                     ista = msta(ispa2)
                     iend = iacon1(igaa)-1
                  endif
c
c  loop over gaps
c
                  do 4860 igap=igaa,igae+1
c
                     do 4850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
c
                     ind = index(jj) + io1
c
              call rede00(iacon1,iacon2,na,ia1,igap-1,jj,jpera)
              if (lbox2(ispa2).gt.iama(ispa2)) go to 4800
c
c  get group number
c
           call positco(lbox5,nspace,na,iama,iami,lbox4,lbox2,iga2)
           nias = nast(iga2)
c
              call idpost(iacon2,na,lbox2,nspace,msta,idim,y,nx,lbst,
     *                  landet(1,iga2),ibcon1,jposa)
              kapos = jposa + nias
              jpza2 = lspa(kapos)
              kasym = lsyma(kapos)
              kper1 = ((-1)**jpera)*2
              immc(kasym) = immc(kasym) + 1
              jspo = immc(kasym)
              iposa(jspo,kasym) = jpza2
              ipera(jspo,kasym) = kper1
              iind1(jspo,kasym) = ind
              igroa(jspo,kasym) = iga2
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) go to 4800
c
c  loop over appropriate beta dets and update the 1e dm1
c
              do 1700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 1700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              do 1680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
                 dm1(ind) = dm1(ind) + fc
 1680         continue
c
 1700         continue
c
c   determine the alpha string contribution to the matrix element.
c
              do 4712 ik=1,na
                 if (ik.eq.ia1) go to 4712
                 ion = iacon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
c
c  loop over appropriate beta dets and update the 1e dm2
c
              do 1705 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 1705
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              do 1685 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
                 dm2(jj1) = dm2(jj1) + fc
                 dm2(inx) = dm2(inx) - fc
 1685         continue
c
 1705         continue
c
 4712         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga2).ne.1) goto 4690
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci2 = jpza2 + ldisb(ksym,igb,iga2)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
                 jci2 = jci2 + 1
c
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
c
                 do 4670 ik=1,nb
                    ion = ibcon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    dm2(jj1) = dm2(jj1) + fc
 4670            continue
c
 4680         continue
c
 4690         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4700         continue
c
c --  double alpha excitations start here  ---
c
 4800         continue
c
            if (ia1.eq.na) goto 4850
            if (jj.eq.nact) goto 4850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            iaes3 = 1
            do kk=1,ispa1-1
               iaes3 = iaes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispa1.
c
            do 3890 ispa3 = ispa1,nspace
               if (lbox1(ispa3).eq.0) goto 3887
               if (ispa3.eq.ispa1.and.ia1.eq.ieae) goto 3887
               ioc3 = lbox2(ispa3)
               if (ioc3.eq.0) goto 3887
c
               iaee3 = iaes3 + lbox2(ispa3)-1
               lbox3(ispa3) = lbox3(ispa3)-1
c
c  loop over electrons in ispa3, which are larger than ia1,
c  making sure it isn't the already excited electron.
c
               jsta3 = iaes3
               if (ispa3.eq.ispa1) then
                  jsta3=ia1
                  if (igap-1.eq.ia1) jsta3=jsta3+1
               endif
c
               do 3880 ia3=jsta3,iaee3
                  if (ia3.eq.igap-1) goto 3880
                  io3 = iacon2(ia3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igae3 = 0
                  do jik=1,ispa2-1
                     igae3 = igae3 + lbox2(jik)
                  enddo
                  do 3850 ispa4=ispa2,nspace
c
c  igas3, igae3 are electrons specifying ispa4 space electron limits.
c
                  lbox3(ispa4) = lbox3(ispa4) + 1
                  igas3 = igae3 + 1
                  igae3 = igae3 + lbox2(ispa4)
                  if (lbox3(ispa3).lt.iami(ispa3)) goto 3840
                  if (lbox3(ispa4).gt.iama(ispa4)) goto 3840
             if (ispa4.eq.ispa2.and.jj.eq.msta(ispa2+1)-1) goto 3840
c
c  get group number
c
               call positco(lbox5,nspace,na,iama,iami,lbox4,lbox3,iga3)
               nias3 = nast(iga3)
c
c  make gap information here.
c
                  igaa3 = max(igas3,igap)
c
                  if (lbox2(ispa4).eq.0) then
                     ista3 = msta(ispa4)
                     iend3 = msta(ispa4+1)-1
                  elseif (ispa4.eq.ispa2) then
                     ista3 = jj+1
                     iend3 = iacon2(igaa3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igae3) iend3=msta(ispa2+1)-1
c
                  else
                     ista3 = msta(ispa4)
                     iend3 = iacon2(igaa3)-1
                  endif
c
c  loop over gaps
c
                  do 3830 igap3=igaa3,igae3+1
c
                     do 3820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 3820
c
              call rede00(iacon2,ibcon1,na,ia3,igap3-1,jj3,jpera3)
              call idpost(ibcon1,na,lbox3,nspace,msta,idim,y,nx,lbst,
     *                  landet(1,iga3),ibcon2,jposa3)
              iper3 = ((-1)**(jpera3+jpera))*2
              kapos3 = jposa3 + nias3
              jpza3 = lspa(kapos3)
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
c
c  loop over beta strings of right symmetry and group.
c
              do 3700 igb=1,itgb
              if (lgcom(igb,iga).ne.1.or.lgcom(igb,iga3).ne.1) goto 3700
              jci1 = jpza1 + ldisb(ksym,igb,iga)
              jci3 = jpza3 + ldisb(ksym,igb,iga3)
c
              do 3680 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
                 jci3 = jci3 + 1
c
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci3,kki)
                 enddo
                 fc = fc*iper3
                 dm2(inx) = dm2(inx) + fc
                 dm2(inx2) = dm2(inx2) - fc
 3680         continue
c
 3700         continue
c
 3820                continue
c
                     ista3 = iacon2(igap3)+1
                     iend3 = iacon2(igap3+1)-1
                     if (igap3.eq.igae3) iend3=msta(ispa4+1)-1
 3830             continue
c
 3840             lbox3(ispa4) = lbox3(ispa4) - 1
 3850             continue
c
 3880          continue
c
               lbox3(ispa3) = lbox3(ispa3)+1
 3887          iaes3 = iaes3 + lbox3(ispa3)
 3890       continue
c
 4850                continue
c
                  ista = iacon1(igap)+1
                  iend = iacon1(igap+1)-1
                  if (igap.eq.igae) iend=msta(ispa2+1)-1
 4860             continue
c
 4870             lbox2(ispa2) = lbox2(ispa2) - 1
 4880          continue
c
 4885          continue
c
               lbox2(ispa1) = lbox2(ispa1)+1
 4890       continue
c
c  diagonal elements here
c
            do 67 ii=1,na
               i1 = iacon1(ii)
               ind1 = index(i1+1)
c
c  loop over beta strings of right symmetry and group.
c
              do 3705 igb=1,itgb
              if (lgcom(igb,iga).ne.1) goto 3705
              jci1 = jpza1 + ldisb(ksym,igb,iga)
c
              do 3685 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
c
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci1,kki)
                 enddo
                 dm1(ind1) = dm1(ind1) + fc
 3685         continue
c
 3705         continue
c
              do 64 jj=ii+1,na
                 i2 = iacon1(jj)
                 ind2 = index(i2+1)
                 indm = ind2 - i2 + i1
                 j1 = index(indm+1)
                 j2 = index(ind2) + ind1
c
              do 3710 igb=1,itgb
              if (lgcom(igb,iga).ne.1) goto 3710
              jci1 = jpza1 + ldisb(ksym,igb,iga)
c
              do 3690 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 jci1 = jci1 + 1
c
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci1,kki)
                 enddo
                 dm2(j1) = dm2(j1) - fc
                 dm2(j2) = dm2(j2) + fc
 3690         continue
c
 3710         continue
c
   64         continue
c
c  loop over beta dets of the right group and symmetry.
c
              call resetco(lbox3,nspace,nb,ibma,ibmi,lbox4)
c
              do 4709 igb=1,itgb
              if (lgcom(igb,iga).ne.1) goto 4699
              jci1 = jpza1 + ldisb(ksym,igb,iga)
c
              call resetde(lbox3,nspace,nb,msta,ibcon1)
              istb = 1
              do 4689 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                 ienb = lsbc(kkb)
                 do iiz=istb,ienb-1
                    call moveup2(lbox3,nspace,nb,msta,ibcon1)
                 enddo
                 istb = ienb
                 jci1 = jci1 + 1
c
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci1,kki)
                 enddo
c
                 do 4679 ik=1,nb
                    i2 = ibcon1(ik)
                    ind2 = index(i2+1)
                    jmi = min(ind1,ind2)
                    jma = max(ind1,ind2)
                    j2 = index(jma) + jmi
                    dm2(j2) = dm2(j2) + fc
 4679            continue
c
 4689         continue
c
 4699         call pushco(lbox3,nspace,nb,ibma,ibmi,lbox4,iend)
 4709         continue
c
   67       continue
c
c  --- end of loop over single alpha excitations. ---
c      now to sort them by positions within symmetries.
c
            do ii=1,nsym
               call fccsrt3(igroa(1,ii),ipera(1,ii),iind1(1,ii),
     *                   iposa(1,ii),immc(ii))
            enddo
c
c  --- end of loop over pure alpha excitations.
c  now to loop over all simultaneous ab -> a'b' excitations.
c
      if (nspace.eq.1) goto 3400
c
c  ***** general case of more than one space !!!!! *******
c
      if (.not.fdirct) then
c
c  loop over single alpha excites within each symmetry.
c
      do 3000 isae=1,nsym
         kbsym=ktab(isae)
         do 2900 jsae=1,immc(isae)
            jposae=iposa(jsae,isae)
            jperae=ipera(jsae,isae)
            jindae=iind1(jsae,isae)
            jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym.and.iga.eq.jgroae) then
c
c  loop over all relevant betas
c
            labpos = jpza1
            labpos2 = jposae
            do 2813 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2813
               nibs = nbst(igb)
               do 2763 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  labpos2 = labpos2 + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2613 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2613
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  jcib2=jpza1+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  fc = 0.0d+00
                  fc1 = 0.0d+00
                 do kki=1,nxtw
           fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
           fc1 = fc1 + wstate(iwts(kki))*ab(labpos2,kki)*ab(jcib2,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  fc1 = fc1*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc + fc1
c
 2613          continue
c
 2763          continue
 2813       continue
c
       else
c
c  loop over all relevant (a-)beta dets.
c
            labpos = jpza1
            do 2800 igb=1,itgb
               if (lgcom(igb,iga).ne.1) goto 2800
               nibs = nbst(igb)
               do 2750 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2600 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,jgroae).ne.1) goto 2600
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib=jposae+ldisb(kbsym,igb2,jgroae)+jb1po(jbindx)
                  fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc
c
 2600          continue
c
 2750          continue
 2800       continue
c
c  loop over all relevant (a'-)beta dets.
c
            labpos = jposae
            do 2803 igb=1,itgb
               if (lgcom(igb,jgroae).ne.1) goto 2803
               nibs = nbst(igb)
               do 2753 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                  labpos = labpos + 1
                  ibpos = nibs + lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2603 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  igb2=jb1gr(jbindx)
                  if (lgcom(igb2,iga).ne.1) goto 2603
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib=jpza1+ldisb(ksym,igb2,iga)+jb1po(jbindx)
                  fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc
c
 2603          continue
c
 2753          continue
 2803       continue
c
      endif
c
 2900    continue
 3000 continue
c
      else
c
c  ***** direct method below *****
c
      do 4000 isae=1,nsym
         kbsym = ktab(isae)
c
c  analyse excited a' groups for compatibility with b.
c
      ngrps = 0
      do 1976 ii=1,immc(isae)
         icgr = igroa(ii,isae)
         do jj=1,ngrps
            if (jb1gr(jj).eq.icgr) goto 1976
         enddo
         ngrps = ngrps + 1
         jb1gr(ngrps) = icgr
 1976 continue
c
c first type of betas, ksym -> kbsym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7400 igb=1,itgb
             if (lgcom(igb,iga).ne.1) goto 7399
             lab1 = jpza1 + ldisb(ksym,igb,iga)
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
             do 7300 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                labpos = lab1 + lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry kbsym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7290 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7290
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7285 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7280 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7270
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7270
c
c  get group number
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          nias = nbst(igb2)
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb2,iagrp).eq.1) goto 7555
             enddo
             goto 7270
 7555        continue
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7260 igap=igba,igbe+1
c
                     do 7250 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7250
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,y,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jb1p = lspb(iposb + nias)
c
            do 3900 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb2,jgroae).ne.1) goto 3900
c
               jposae=iposa(jsae,isae)
               jcib = jb1p + jposae + ldisb(kbsym,igb2,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
c
             fc=0.0d+00
             do kki=1,nxtw
            fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
             enddo
             fc = fc*jperae*qjper
             dm2(ix) = dm2(ix) + fc
c
 3900       continue
c
 7250                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7260             continue
c
 7270              lbox3(ispb2) = lbox3(ispb2) - 1
 7280           continue
c
 7285           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7290           continue
c
 7300        continue
c
 7399        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7400     continue
c
c
c second type of betas, kbsym -> ksym
c
          call resetco(lbox2,nspace,nb,ibma,ibmi,lbox3)
          do 7801 igb=1,itgb
c
c check for a' group compatibility
c
             do ii=1,ngrps
                iagrp = jb1gr(ii)
                if (lgcom(igb,iagrp).eq.1) goto 7655
             enddo
             goto 7799
 7655        continue
c
             call resetde(lbox2,nspace,nb,msta,ibcon1)
             nibs = nbst(igb)
             istb = 1
c
             do 7701 kkb=lsbs(kbsym,igb),lsbs(kbsym+1,igb)-1
                ienb = lsbc(kkb)
                do iiz=istb,ienb-1
                   call moveup2(lbox2,nspace,nb,msta,ibcon1)
                enddo
                istb = ienb
c
                ibpos = nibs + lsbc(kkb)
                lab1 = lspb(ibpos)
c
c loop over single beta excites from ibpos which are of symmetry ksym
c
                mesym1 = lgmul(ksym,kbsym)
                do ii=1,nspace
                   lbox3(ii) = lbox2(ii)
                enddo
                iebs = nb + 1
                do 7691 ispb1=nspace,1,-1
                   ioc1 = lbox2(ispb1)
                   iebe = iebs - 1
                   iebs = iebs - ioc1
                   if (ioc1.eq.0) go to 7691
                   lbox3(ispb1) = lbox3(ispb1)-1
c
c  loop electrons in space ispb1
c  iebs, iebe are the electrons in space ispb1
c
                do 7686 ib1 = iebe,iebs,-1
                   io1 = ibcon1(ib1)
                   mesym2 = lgmul(mesym1,iob(io1))
                   igbe = iebe - lbox2(ispb1)
c
c  loop over possible spaces to excite into.
c
                do 7681 ispb2=ispb1,nspace
c
c  igbs,igbe are electrons specifying ispb2 space electron limits.
c
                   igbs = igbe + 1
                   igbe = igbe + lbox2(ispb2)
c
                   lbox3(ispb2) = lbox3(ispb2) + 1
                   if (lbox3(ispb1).lt.ibmi(ispb1)) go to 7671
                   if (lbox3(ispb2).gt.ibma(ispb2)) go to 7671
c
          call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb2)
          if (lgcom(igb2,iga).ne.1) goto 7671
          jcibs = jpza1 + ldisb(ksym,igb2,iga)
          nias = nbst(igb2)
c
c  make gap information here.
c
                  igba = max(ib1+1,igbs)
                  if (lbox2(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igba)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igba)-1
                  endif
c
c  loop over gaps
c
                  do 7660 igap=igba,igbe+1
c
                     do 7650 jj=ista,iend
c
c check to see if excited beta is of right symmetry (kbsym)
c
            if (iob(jj).ne.mesym2) goto 7650
c
            ind1 = index(jj) + io1
            call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,iper)
            call idpost(ibcon2,nb,lbox3,nspace,msta,idim,y,nx,lbst,
     *                  lbndet(1,igb2),iacon2,iposb)
            qjper = ((-1)**iper)
            jcib= jcibs + lspb(iposb + nias)
c
            do 4300 jsae=1,immc(isae)
               jgroae=igroa(jsae,isae)
               if (lgcom(igb,jgroae).ne.1) goto 4300
c
               jposae=iposa(jsae,isae)
               labpos = lab1 + jposae + ldisb(kbsym,igb,jgroae)
c
               jindae=iind1(jsae,isae)
               jma=max(jindae,ind1)
               jmi=min(jindae,ind1)
               ix=index(jma)+jmi
c
               jperae=ipera(jsae,isae)
               fc = 0.0d+00
               do kki=1,nxtw
            fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
               enddo
               fc = fc*jperae*qjper
               dm2(ix) = dm2(ix) + fc
c
 4300       continue
c
 7650                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7660             continue
c
 7671              lbox3(ispb2) = lbox3(ispb2) - 1
 7681           continue
c
 7686           continue
c
                   lbox3(ispb1) = lbox3(ispb1) + 1
 7691           continue
c
 7701        continue
c
 7799        call pushco(lbox2,nspace,nb,ibma,ibmi,lbox3,iend)
 7801     continue
c
 4000 continue
c
      endif
c
c  ***** end of direct option *****
c
c  **** end of general case of more than one space ******
c
      goto 4899
c
 3400 continue
c
c ***** special case of one space !!!!! ******
c
c  loop over single alpha excites within each symmetry.
c
       do 2901 isae=1,nsym
          kbsym=ktab(isae)
          do 2801 jsae=1,immc(isae)
                jposae=iposa(jsae,isae)
                jperae=ipera(jsae,isae)
                jindae=iind1(jsae,isae)
                jgroae=igroa(jsae,isae)
c
c  if isae.eq.jasym then special case.
c
       if (isae.eq.jasym) then
c
c  loop over all relevant betas
c
             labpos=jpza1
             labpos2=jposae
                do 2721 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    labpos2 = labpos2 + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2621 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib = jposae+jb1po(jbindx)
                  jcib2 = jpza1+jb1po(jbindx)
                  fc = 0.0d+00
                  fc1 = 0.0d+00
                 do kki=1,nxtw
          fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
          fc1 = fc1 + wstate(iwts(kki))*ab(labpos2,kki)*ab(jcib2,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  fc1 = fc1*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc + fc1
c
 2621          continue
c
 2721          continue
c
      else
c
c  loop over all relevant (a-)beta dets.
c
                labpos = jpza1
                do 2751 kkb=lsbs(ksym,1),lsbs(ksym+1,1)-1
                    labpos = labpos + 1
                    ibpos =  lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2601 jbindx=jb1st(kbsym,ibpos),jb1st(kbsym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib=jposae+jb1po(jbindx)
                  fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc
c
 2601          continue
c
 2751          continue
c
c  loop over all relevant (a'-)beta dets.
c
               labpos = jposae
               do 2761 kkb=lsbs(kbsym,1),lsbs(kbsym+1,1)-1
                  labpos = labpos + 1
                  ibpos = lsbc(kkb)
c
c  loop over single beta excites from ibpos
c
               do 2611 jbindx=jb1st(ksym,ibpos),jb1st(ksym+1,ibpos)-1
c
                  jma=max(jindae,jb1in(jbindx))
                  jmi=min(jindae,jb1in(jbindx))
                  ix=index(jma)+jmi
                  jcib=jpza1+jb1po(jbindx)
                  fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(labpos,kki)*ab(jcib,kki)
                 enddo
                  fc = fc*jperae*jb1pe(jbindx)
                  dm2(ix) = dm2(ix) + fc
c
 2611          continue
c
 2761          continue
c
      endif
c
 2801       continue
 2901   continue
c
c  **** end of special case of one space *******
c
 4899       call moveup2(lbox1,nspace,na,msta,iacon1)
 4900    continue
c
         call pushco(lbox1,nspace,na,iama,iami,lbox5,iend)
 5000 continue
c
c  --- end of loop over alpha strings. ---
c **
c  --- loop over all pure beta excitations.
c
      call resetco(lbox1,nspace,nb,ibma,ibmi,lbox5)
c
      do 8000 igb=1,itgb
c
         call resetde(lbox1,nspace,nb,msta,ibcon1)
c
c  kkb gives the actual position of the beta string ibcon1 in
c  the full beta string list.
c
         do 7900 kkb=nbst(igb)+1,nbst(igb+1)
            jpzb1 = lspb(kkb)
            kbsym = lsymb(kkb)
            ksym=ktab(kbsym)
c
            do ii=1,nspace
               lbox2(ii) = lbox1(ii)
            enddo
c
c  loop over spaces to excite electron from.
c
            iebs = nb+1
            do 7890 ispb1=nspace,1,-1
               ioc1 = lbox1(ispb1)
               iebe = iebs - 1
               iebs = iebs - ioc1
               if (ioc1.eq.0) goto 7890
               lbox2(ispb1) = lbox2(ispb1)-1
c
c  loop electrons in space ispb1.
c  iebs, iebe are the electrons in space ispb1.
c
               do 7885 ib1=iebe,iebs,-1
                  io1 = ibcon1(ib1)
                  igbe = iebe - lbox1(ispb1)
                  is1 = iob(io1)
c
c  loop over possible spaces to excite into.
c
               do 7880 ispb2=ispb1,nspace
c
c  igbs, igbe are electrons specifying ispb2 space electron limits.
c
                  igbs = igbe + 1
                  igbe = igbe + lbox1(ispb2)
c
                  lbox2(ispb2) = lbox2(ispb2) + 1
                  if (lbox2(ispb1).lt.ibmi(ispb1)) goto 7870
c
c  make gap information here.
c
                  igbb = max(ib1+1,igbs)
                  if (lbox1(ispb2).eq.0) then
                     ista = msta(ispb2)
                     iend = msta(ispb2+1)-1
                  elseif (ispb2.eq.ispb1) then
                     ista = io1+1
                     iend = ibcon1(igbb)-1
                     if (ib1.eq.iebe) iend=msta(ispb1+1)-1
                  else
                     ista = msta(ispb2)
                     iend = ibcon1(igbb)-1
                  endif
c
c  loop over gaps
c
                  do 7860 igap=igbb,igbe+1
c
                     do 7850 jj=ista,iend
                        is2 = iob(jj)
                        ip1 = lgmul(is1,is2)
                        ind = index(jj) + io1
c
              call rede00(ibcon1,ibcon2,nb,ib1,igap-1,jj,jperb)
                  if (lbox2(ispb2).gt.ibma(ispb2)) goto 7800
c
c   if deoccupied and newly occupied orbitals are of different symmetry,
c   skip to doubles.
c
              if (is1.ne.is2) goto 7800
c
c  get group number
c
           call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox2,igb2)
           nibs = nbst(igb2)
c
              call idpost(ibcon2,nb,lbox2,nspace,msta,idim,y,nx,lbst,
     *                  lbndet(1,igb2),iacon1,jposb)
              kbpos = jposb + nibs
              jpzb2 = lspb(kbpos)
              kper1 = ((-1)**jperb)*2
c
c  loop over alpha and update dm
c
              do 7705 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7705
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              do 7685 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
                 dm1(ind) = dm1(ind) + fc
 7685         continue
c
 7705         continue
c
              do 7712 ik=1,nb
                 if (ik.eq.ib1) goto 7712
                 ion = ibcon1(ik)
                 j1 = index(ion+1)
                 jma = max(j1,ind)
                 jmi = min(j1,ind)
                 jj1 = index(jma) + jmi
                 jma = max(ion,jj)
                 jmi = min(ion,jj)
                 j1 = index(jma)+jmi
                 jma = max(io1,ion)
                 jmi = min(io1,ion)
                 j2 = index(jma)+jmi
                 jma = max(j1,j2)
                 jmi = min(j1,j2)
                 inx = index(jma)+jmi
c
              do 7710 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7710
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              do 7695 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
                 dm2(jj1) = dm2(jj1) + fc
                 dm2(inx) = dm2(inx) - fc
 7695         continue
c
 7710         continue
c
 7712         continue
c
c  loop over alpha strings of the right group and symmetry.
c
              call resetco(lbox3,nspace,na,iama,iami,lbox4)
c
              do 7700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb2,iga).ne.1) goto 7690
                 jcib1 = ldisb(kbsym,igb,iga) + jpzb1
                 jcib2 = ldisb(kbsym,igb2,iga) + jpzb2
                 nias = nast(iga)
c
              call resetde(lbox3,nspace,na,msta,iacon1)
              ista1 = 1
              do 7680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena1 = lsac(kka)
                 do iiz=ista1,iena1-1
                    call moveup2(lbox3,nspace,na,msta,iacon1)
                 enddo
                 ista1 = iena1
                 jcia = lspa(nias+iena1)
                 jci1 = jcia + jcib1
                 jci2 = jcia + jcib2
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci2,kki)
                 enddo
                 fc = fc*kper1
c
                 do 7670 ik=1,na
                    ion = iacon1(ik)
                    j1 = index(ion+1)
                    jma = max(j1,ind)
                    jmi = min(j1,ind)
                    jj1 = index(jma) + jmi
                    dm2(jj1) = dm2(jj1) + fc
 7670            continue
c
 7680         continue
c
 7690         call pushco(lbox3,nspace,na,iama,iami,lbox4,iend)
 7700         continue
c
c  --- double excitations start here
c
 7800         continue
c
           if (ib1.eq.nb) goto 7850
           if (jj.eq.nact) goto 7850
c
            do ii=1,nspace
               lbox3(ii) = lbox2(ii)
            enddo
            ibes3 = 1
            do kk=1,ispb1-1
               ibes3 = ibes3 + lbox2(kk)
            enddo
c
c  loop over spaces to excite electrons from, must be ge than
c  space of first excitation, ispb1.
c
            do 6890 ispb3 = ispb1,nspace
               if (lbox1(ispb3).eq.0) goto 6887
               if (ispb3.eq.ispb1.and.ib1.eq.iebe) goto 6887
               ioc3 = lbox2(ispb3)
               if (ioc3.eq.0) goto 6887
c
               ibee3 = ibes3 + lbox2(ispb3)-1
               lbox3(ispb3) = lbox3(ispb3)-1
c
c  loop over electrons in ispb3, which are larger than ib1,
c  making sure it isn't the already excited electron.
c
               jstb3 = ibes3
               if (ispb3.eq.ispb1) then
                  jstb3=ib1
                  if (igap-1.eq.ib1) jstb3=jstb3+1
               endif
c
               do 6880 ib3=jstb3,ibee3
                  if (ib3.eq.igap-1) goto 6880
                  io3 = ibcon2(ib3)
                  is3 = iob(io3)
c
c  loop over spaces to excite electron into.  must be ge than
c  space first electron was excited into.
c
                  igbe3 = 0
                  do jik=1,ispb2-1
                     igbe3 = igbe3 + lbox2(jik)
                  enddo
                  do 6850 ispb4=ispb2,nspace
c
c  igbs3, igbe3 are electrons specifying ispb4 space electron limits.
c
                  lbox3(ispb4) = lbox3(ispb4) + 1
                  igbs3 = igbe3 + 1
                  igbe3 = igbe3 + lbox2(ispb4)
                  if (lbox3(ispb3).lt.ibmi(ispb3)) goto 6840
                  if (lbox3(ispb4).gt.ibma(ispb4)) goto 6840
             if (ispb4.eq.ispb2.and.jj.eq.msta(ispb2+1)-1) goto 6840
c
c  get group number
c
               call positco(lbox5,nspace,nb,ibma,ibmi,lbox4,lbox3,igb3)
               nibs3 = nbst(igb3)
c
c  make gap information here.
c
                  igbb3 = max(igbs3,igap)
c
                  if (lbox2(ispb4).eq.0) then
                     ista3 = msta(ispb4)
                     iend3 = msta(ispb4+1)-1
                  elseif (ispb4.eq.ispb2) then
                     ista3 = jj+1
                     iend3 = ibcon2(igbb3)-1
c
c  i am suspect about this next line, we'll see what happens.......
c
                     if (igap-1.eq.igbe3) iend3=msta(ispb2+1)-1
c
                  else
                     ista3 = msta(ispb4)
                     iend3 = ibcon2(igbb3)-1
                  endif
c
c  loop over gaps
c
                  do 6830 igap3=igbb3,igbe3+1
c
                     do 6820 jj3=ista3,iend3
                        is4 = iob(jj3)
                        ip2 = lgmul(is3,is4)
                        if (ip1.ne.ip2) goto 6820
c
              call rede00(ibcon2,iacon1,nb,ib3,igap3-1,jj3,jperb3)
              call idpost(iacon1,nb,lbox3,nspace,msta,idim,y,nx,lbst,
     *                  lbndet(1,igb3),iacon2,jposb3)
              iper3 = ((-1)**(jperb3+jperb))*2
              kbpos3 = jposb3 + nibs3
              jpzb3 = lspb(kbpos3)
c
                 jma=max(jj3,io3)
                 jmi=min(jj3,io3)
                 i2 = index(jma) + jmi
                 inx = index(i2) + ind
                 ii1 = index(jj3) + io1
                 jma=max(io3,jj)
                 jmi=min(io3,jj)
                 ii2 = index(jma) + jmi
                 jma=max(ii1,ii2)
                 jmi=min(ii1,ii2)
                 inx2 = index(jma) + jmi
c
c  loop over alpha strings of right symmetry and group.
c
              do 6700 iga=1,itga
              if (lgcom(igb,iga).ne.1.or.lgcom(igb3,iga).ne.1) goto 6700
              jcib1 = jpzb1 + ldisb(kbsym,igb,iga)
              jcib3 = jpzb3 + ldisb(kbsym,igb3,iga)
              nias = nast(iga)
c
              do 6680 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena3 = lsac(kka)
                 jcia = lspa(nias+iena3)
                 jci1 = jcia + jcib1
                 jci3 = jcia + jcib3
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci3,kki)
                 enddo
                 fc = fc*iper3
                 dm2(inx) = dm2(inx) + fc
                 dm2(inx2) = dm2(inx2) - fc
 6680         continue
c
 6700         continue
c
 6820                continue
c
                     ista3 = ibcon2(igap3)+1
                     iend3 = ibcon2(igap3+1)-1
                     if (igap3.eq.igbe3) iend3=msta(ispb4+1)-1
 6830             continue
c
 6840             lbox3(ispb4) = lbox3(ispb4) - 1
 6850             continue
c
 6880          continue
c
               lbox3(ispb3) = lbox3(ispb3)+1
 6887          ibes3 = ibes3 + lbox3(ispb3)
 6890       continue
c
c  --- end of loop over double beta excitations.
c
 7850                continue
c
                  ista = ibcon1(igap)+1
                  iend = ibcon1(igap+1)-1
                  if (igap.eq.igbe) iend=msta(ispb2+1)-1
 7860             continue
c
 7870             lbox2(ispb2) = lbox2(ispb2) - 1
 7880          continue
c
 7885          continue
c
               lbox2(ispb1) = lbox2(ispb1)+1
 7890       continue
c
c remaining diagonal contributions here
c
            do 69 ii=1,nb
               i1 = ibcon1(ii)
               ind1 = index(i1+1)
c
              do 6705 iga=1,itga
              if (lgcom(igb,iga).ne.1) goto 6705
              jcib1 = jpzb1 + ldisb(kbsym,igb,iga)
              nias = nast(iga)
c
              do 6685 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena3 = lsac(kka)
                 jcia = lspa(nias+iena3)
                 jci1 = jcia + jcib1
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci1,kki)
                 enddo
                 dm1(ind1) = dm1(ind1) + fc
 6685         continue
c
 6705         continue
c
              do 74 jj=ii+1,nb
                 i2 = ibcon1(jj)
                 ind2 = index(i2+1)
                 indm = ind2 - i2 + i1
                 j1 = index(indm+1)
                 jmi = min(ind1,ind2)
                 jma = max(ind1,ind2)
                 j2 = index(jma) + jmi
c
              do 6710 iga=1,itga
              if (lgcom(igb,iga).ne.1) goto 6710
              jcib1 = jpzb1 + ldisb(kbsym,igb,iga)
              nias = nast(iga)
c
              do 6690 kka=lsas(ksym,iga),lsas(ksym+1,iga)-1
                 iena3 = lsac(kka)
                 jcia = lspa(nias+iena3)
                 jci1 = jcia + jcib1
                 fc = 0.0d+00
                 do kki=1,nxtw
              fc = fc + wstate(iwts(kki))*ab(jci1,kki)*ab(jci1,kki)
                 enddo
                 dm2(j1) = dm2(j1) - fc
                 dm2(j2) = dm2(j2) +  fc
 6690         continue
c
 6710         continue
c
   74         continue
c
   69         continue
c
            call moveup2(lbox1,nspace,nb,msta,ibcon1)
 7900    continue
c
         call pushco(lbox1,nspace,nb,ibma,ibmi,lbox5,iend)
 8000 continue
c
c  --- end of loop over beta strings. ---
c
c ----    end of density matrix generation ----------
c
      call valfm(loadfm)
      lwrk   = loadfm + 1
      llabmo = lwrk   + nocc2
      llbabl = llabmo + l1
      llbirp = llbabl + m1
      lsyirp = llbirp + 12
      last   = lsyirp + 12
      need   = last   - loadfm - 1
      call getfm(need)
c
      call detgrp(grpdet,   iob   ,x(llbabl),ptgrp,x(llbirp),
     *            x(lsyirp),nsym,nirrp,l1,nact,ncorsv)
c
      cutoff = max(1.0d-11,10.0d+00**(-icut))
      if(some) write(iw,9370) x(lsyirp),ptgrp
c
      call wtdm12(dm1,dm2,x(llbabl),
     *            m1,m2,m4,x(lwrk),nocc2,ncorsv,cutoff)
c
      call retfm(need)
c
      if(some) write(iw,9390)
      if(some) call timit(1)
c
      return
c
 9310 format(/5x,55(1h-)/
     *   5x,' one and two particle ormas density matrix computation'/
     *   18x,'program written by joe ivanic'/
     *   5x,55(1h-))
 9320 format(/1x,'the densities are state averaged over',i4,' root(s)')
 9340 format(1x,'state=',i4,'   energy=',f20.10,'   weight=',f8.5,
     *           '   s=',f6.2)
 9350 format(/1x,'***** error *****'/
     *       1x,'this run found',i5,' ci eigenvectors with s=',f5.2,','/
     *       1x,'but you requested state averaging of',i5,' roots.'/
     *       1x,'please examine your choice of -nstate- input data.'/)
 9370 format(1x,'sieving the ',a4,
     *          ' symmetry nonzero density elements in group ',a8)
 9380 format(1x,i10,' nonzero dm2 elements written in',i8,
     *          ' records to file',i3)
 9390 format(1x,'..... done with 1 and 2 particle density matrix .....')
c
      end
c
c
c*module ormas2  *deck refwe
c     ----------------------------------------------------------------
      subroutine refwe(vec,nast,itga,itgb,
     *                 lsyma,iast,lgcom,lsbs,nsym,ktab)
c     ----------------------------------------------------------------
      implicit double precision(a-h,o-z)
c
      parameter (mxrt=100)
c
      logical fdirct,qcorr
c
      common /detwfn/ wstate(mxrt),spins(mxrt),crit,prttol,s,sz,
     *                grpdet,stsym,glist,
     *                nflgdm(mxrt),iwts(mxrt),ncorsv,ncor,nact,norb,
     *                na,nb,nstate,kst,iroot,ipures,maxw1,niter,
     *                maxp,nci,igpdet,kstsym
      common /fccwfn/ nspace,msta(51),mnum(51),mini(51),maxi(51),
     *                iami(51),iama(51),ibmi(51),ibma(51),idim(51),
     *                lbst(51),nref0,fdirct,qcorr,c0sq
c
      dimension vec(nci)
      dimension nast(itga),lsyma(iast)
      dimension lgcom(itgb,itga),lsbs(nsym+1,itgb),ktab(nsym)
c
      c0sq = 0.0d+00
      ici = 0
      icount = 0
c
      iga = 1
      do 410 kka=nast(iga)+1,nast(iga+1)
         jasym=lsyma(kka)
         ksym=ktab(jasym)
         do 500 igb=1,itgb
            if (lgcom(igb,iga).ne.1) goto 500
            do 510 kkb=lsbs(ksym,igb),lsbs(ksym+1,igb)-1
               ici=ici+1
               if (igb.eq.1) then
                  c0sq = c0sq + vec(ici)**2.0d+00
                  icount = icount + 1
               endif
  510       continue
  500    continue
  410 continue
c
      return
      end
c
c     layer to gamess-uk i/o
c
      subroutine put_mo_coeffs(coeffs,len)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/atmblk)
      integer len,iblok
      REAL coeffs(*)
      call secget(isect(8),8,iblok)
      iblok = iblok + mvadd
      call wrt3(coeffs,len,iblok,ifild)
      end 

      subroutine put_mo_coeffs2(coeffs,len)
INCLUDE(common/restar)
INCLUDE(common/restri)
INCLUDE(common/atmblk)
      integer len,iblok
      REAL coeffs(*)
      len2=lensec(len)+mvadd
      call secput(isect(8),8,len2,iblok)
      iblok = iblok + mvadd
      call wrt3(coeffs,len,iblok,ifild)
      end

