c     deck=guess
c ******************************************************
c ******************************************************
c             =   guess   =
c ******************************************************
c ******************************************************
_EXTRACT(adapt,mips4)
**==adapt.f
      subroutine adapt(sexp,coe,bb,cbb,ip,nfun,
     * npar,mpar,ksym,lsym,iwrit,newb)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *3 txtsym
INCLUDE(common/sizes)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/tran)
INCLUDE(common/harmon)
INCLUDE(common/nshel)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/restri)
INCLUDE(common/gjs)
      common/junk/p(maxorb),q(maxorb),r(maxorb),
     *i2(maxorb),i3(maxorb),i4(maxorb),nord(maxorb),
     *nco(maxat),nrw(maxat),ityp(maxorb),icon(maxorb),iord(maxat),
     *ncont(maxat),natt(maxat),ise(7),ibas( maxat),npf(7),npz(7),nel(3),
     *nbon(maxat),mcomp(maxorb),
     *nfil(maxorb),imix(3),ibuk(3),inuk(3),nzer(8),nblz(8),ncal(8),
     *nzil(8),nir(maxorb),loc(maxorb),ncomp(maxorb),
     *mj(8),nj(8),ntil(8),nbal(9),jtyp(maxorb)
     *,nnum(maxorb),nsog(maxorb),new(maxorb),newei(maxorb)
INCLUDE(common/blockc)
c
      dimension igtab(35),ittab(35),txtsym(9)
c
c magical constants here are 
c                      3=max. no. of generators of abel. group
c                      8=max group order - d2h, 7=8-1
c                     10=max no. of primitive GTO contracted together
c
      dimension sexp(mxgaus),coe(mxgaus),bb(newb,*),cbb(newb,*),
     + ip(7*maxat),nfun(3*newb),npar(3*newb),mpar(3*newb),
     + ksym(8*newb),lsym(8*newb),iwrit(newb)
c
      common /saveco/ u,idadap,ifirst,ijato,klato,ijshlo
c     save idadap
c
      data igtab/1,2,3,4,
     +           5,8,10,6,7,9,
     +          11,16,20,12,13,15,17,18,19,14,
     +          21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /
      data ittab/1,2,3,5,
     +           1,1,1,4,6,7,
     +           2,3,5,3,5,2,5,2,3,8,
     +           1,1,1,4,6,4,7,6,7,1,1,1,7,6,4/
      data m51/51/
      data txtsym/'C1','C2','C2v','Cs','Ci','D2h','D2','C2h','???'/
      data thresh/.001d0/
c
c...   if otran=.true. prepare a dummy adapt
c...   for force etc. it has to remember to do this every time
c...   it is being called    (mfg/jvl 1989)
c...   idadap is checked below together with the check on atomic calc.
c
c      data idadap/-1/
c...  this switches adapt
      if (idadap.eq.-1) then
         if (otran) then
            idadap = 1
            otran = .false.
         else
            idadap = 0
         end if
      end if
c
c ----- symmetry adaption generator ----
c
c ----- first reorder the shell structured basis list (/nshel/)
c       to the atom based list required by adaption code
c
c ----- first generate no. of contracted functions .. nco
c                and   no. of primitive  functions .. nrw
c -----      on each atom
c
c       generate exponent and coeficient list ordered (sexp,coe)
c       for primitives ... atom ordered
c
c ----- build up no. of primitives/gto in icon
c                      and gto type in    ityp
c
      if (oprint(46)) write (iwr,6010)
      do 20 i = 1 , nat
         nco(i) = 0
         nrw(i) = 0
 20   continue
c
      iprim = 0
      nc = 0
c
      do 140 iat = 1 , nat
c
         do 130 ii = 1 , nshell
            i = katom(ii)
            if (i.eq.iat) then
c    shell located at given atom, where it starts and 
c    number of components
               is = kstart(ii)
               ipri = kng(ii)
               if = is + ipri - 1
               mini = kmin(ii)
               maxi = kmax(ii)
               kk = kloc(ii) - mini
c
               do 120 iorb = mini , maxi
                  li = kk + iorb
                  nc = nc + 1
                  new(nc) = li
                  newei(li) = nc
c
c                 spdfg branching
                  go to (30,50,50,50,
     +                   70,70,70,70,70,70,
     +                   90,90,90,90,90,90,90,90,90,90,
     +                   95,95,95,95,95,95,95,95,95,95,
     +                   95,95,95,95,95) , iorb
c
 30               do 40 ig = is , if
                     iprim = iprim + 1
                     coe(iprim) = cs(ig)
                     sexp(iprim) = ex(ig)
 40               continue
                  go to 110
c
 50               do 60 ig = is , if
                     iprim = iprim + 1
                     coe(iprim) = cp(ig)
                     sexp(iprim) = ex(ig)
 60               continue
                  go to 110
c
 70               do 80 ig = is , if
                     iprim = iprim + 1
                     coe(iprim) = cd(ig)
                     sexp(iprim) = ex(ig)
 80               continue
                  go to 110
c
 90               do 100 ig = is , if
                     iprim = iprim + 1
                     coe(iprim) = cf(ig)
                     sexp(iprim) = ex(ig)
 100              continue
                  go to 110
c
 95               do 105 ig = is , if
                     iprim = iprim + 1
                     coe(iprim) = cg(ig)
                     sexp(iprim) = ex(ig)
 105              continue
 110              if (iprim.gt.mxgaus)then
                  write(iwr,6140) iprim,mxgaus
                go to 1100
                endif
                  nco(iat) = nco(iat) + 1
                  nrw(iat) = nrw(iat) + ipri
                  icon(nc) = ipri
                  jtyp(nc) = ittab(iorb)
                  ityp(nc) = igtab(iorb)
 120           continue
            end if
c
 130     continue
 140  continue
c
c
c ----- now classify atoms into type, based solely on charge + basis set
c ----- considerations, ignoring gamess symmetry classification
c
      nbas = 1
      iord(1) = 1
      lnb = nco(1)
      lpr = nrw(1)
c
c...   adapt is special for 1 atom or if adapt off was requested
c...   mfg 1989
c
      if (nat.eq.1 .or. idadap.eq.1) then
         do 150 loop = 1 , num
            ncomp(loop) = 1
 150     continue
         if (nat.eq.1) then
            do 160 loop = 1 , 8
               nj(loop) = 0
 160        continue
            if (idadap.eq.1) then
             nsel = 1
             nj(1) = num
             do 195 i = 1 , num
                lsym(i) = i
 195         continue
             go to 1010
            else
             nbas = 0
             nsel = 0
             do 180 joop = 1 , 8
                do 170 loop = 1 , num
                   if (jtyp(loop).eq.joop) then
                      nbas = nbas + 1
                      lsym(nbas) = loop
                      nj(joop) = nj(joop) + 1
                      nsel = max(nsel,joop)
                   end if
 170            continue
 180         continue
            endif
         else
            nsel = 1
            nj(1) = num
            do 190 i = 1 , num
               lsym(i) = i
 190        continue
         end if
         go to 1010
      end if
c
c...  end special adapt
c
c find all distinct (unique) atoms according to the basis
      do 250 i = 2 , nat
         ig = 0
         kmal = 0
c                       charge
         czani = czan(i)
         ish = nco(i)
         iprim = nrw(i)
         iat = i - 1
         do 230 j = 1 , iat
            czanij = czani - czan(j)
            jsh = nco(j)
            jprim = nrw(j)
            if (jsh.eq.ish) then
               if (jprim.eq.iprim) then
                  if (dabs(czanij).le.thresh) then
                     if (ish.ne.0) then
                        do 200 k = 1 , ish
                           if (icon(lnb+k).ne.icon(kmal+k)) go to 220
                           if (ityp(lnb+k).ne.ityp(kmal+k)) go to 220
 200                    continue
c
                        do 210 k = 1 , jprim
                           if (dabs(sexp(lpr+k)-sexp(ig+k)).gt.1.d-5)
     +                         go to 220
                           if (dabs(coe(lpr+k)-coe(ig+k)).gt.1.d-5)
     +                         go to 220
 210                    continue
                     end if
c                                      is equivalent
                     iord(i) = iord(j)
                     go to 240
                  end if
               end if
            end if
c                          was different
 220        ig = ig + jprim
c
            kmal = kmal + jsh
 230     continue
c
c        new distinct atom
         nbas = nbas + 1
         iord(i) = nbas
 240     lnb = lnb + ish
         lpr = lpr + iprim
 250  continue
c
c    iord is pointer from atom to its representative distinct atom
c    nbas is no. of distinct atoms found
c
c ----- now compress arrays down from atoms to distinct atoms(nbas)
c
c
      jj = 0
      iprim = 0
      lnb = 0
      lpr = 0
      nbas = 0
c
      do 280 i = 1 , nat
         ii = nco(i)
         nra = nrw(i)
         if (iord(i).gt.nbas) then
            nbas = nbas + 1
            nco(nbas) = ii
            nrw(nbas) = nra
c
            if (ii.ne.0) then
c
               do 260 j = 1 , ii
                  jj = jj + 1
                  icon(jj) = icon(lnb+j)
                  ityp(jj) = ityp(lnb+j)
 260           continue
               do 270 j = 1 , nra
c
                  iprim = iprim + 1
c
                  coe(iprim) = coe(lpr+j)
c
                  sexp(iprim) = sexp(lpr+j)
 270           continue
            end if
         end if
         lpr = lpr + nra
         lnb = lnb + ii
 280  continue
      namx = maxat
      nsmx = 7*maxat
      iorbs = 0
      do 290 k = 1 , nat
         ncc = iord(k)
         iorbs = iorbs + nco(ncc)
         ncont(k) = nco(ncc)
 290  continue
      if (iorbs.gt.num)then
      write(iwr,6150) iorbs, num
      go to 1100
      endif
      call setsto(nbas,0,nord)
      do 300 i = 1 , nat
         isal = iord(i)
         nord(isal) = nord(isal) + 1
 300  continue
c in nord are counts of references to distinct atoms
c the following prepares locations and xyz exponents of shells
      iexp = 0
      ican = 0
      ibat = 0
      do 350 i = 1 , nbas
         isal = nord(i)
         ncw = nco(i)
         if (ncw.ne.0) then
            do 340 j = 1 , ncw
               ican = ican + 1
               ithn = icon(ican)
               ittp = ityp(ican)
               ibat = ibat + 1
c              xyz exponents of cartesian gaussian
               nyx = nix(ittp)
               nyy = niy(ittp)
               nyz = niz(ittp)
               do 320 k = 1 , ithn
                  iexp = iexp + 1
                  coff = coe(iexp)
                  expf = sexp(iexp)
                  jbat = ibat - ncw
                  do 310 l = 1 , isal
                     jbat = jbat + ncw
                     bb(jbat,k) = expf
                     cbb(jbat,k) = coff
 310              continue
 320           continue
               jbat = ibat - ncw
               do 330 k = 1 , isal
                  jbat = jbat + ncw
                  nnum(jbat) = ithn
                  i2(jbat) = nyx
                  i3(jbat) = nyy
                  i4(jbat) = nyz
 330           continue
 340        continue
            ibat = jbat
         end if
 350  continue
      jbat = 0
      do 380 i = 1 , nbas
         isal = nord(i)
         ncv = nco(i)
         if (ncv.ne.0) then
            jk = 0
            do 370 k = 1 , nat
               jp = iord(k)
               if (jp.eq.i) then
                  pt = c(1,k)
                  qt = c(2,k)
                  rt = c(3,k)
                  do 360 j = 1 , ncv
                     jbat = jbat + 1
                     p(jbat) = pt
                     q(jbat) = qt
                     r(jbat) = rt
 360              continue
                  jk = jk + 1
                  if (jk.eq.isal) go to 380
               end if
 370        continue
            go to 1100
         end if
 380  continue
      do 390 i = 1 , nat
         natt(i) = nuct(i) - 1
 390  continue
      call setsto(nbas,0,ibas)
      ibs = 0
      do 400 i = 1 , nat
         k = iord(i)
         kkk = natt(i)
         if (kkk.lt.1) go to 400
         if (kkk.eq.1.or.kkk.eq.2) then
            ibs = 1
         else
            ibs = 2
         end if
         ibas(k) = kkk
 400  continue
      call setsto(nsmx,0,ip)
      icx = -namx
c
c loop over all possible group operations with exception of identity
c isx etc. contain x,y,z inversions under group operations
c the following code creates array ip containing transformations
c of atoms under group operations or zero, if that operation does
c not map atoms to ane another
c
      do 500 i = 1 , 7
         icx = icx + namx
         ix = isx(i)
         iy = isy(i)
         iz = isz(i)
         do 430 j = 1 , nbas
            if (ibas(j).eq.0) then
               isal = nord(j)
               isl = 0
               idx = icx
               do 420 k = 1 , nat
                  idx = idx + 1
                  if (natt(k).eq.0 .and. ip(idx).eq.0) then
                     if (iord(k).eq.j) then
                        nfl = 0
                        if (ix.eq.0) then
                           gx = c(1,k)
                        else if (dabs(c(1,k)).lt.1.0d-5) then
                           gx = 0.0d0
                        else
                           gx = -c(1,k)
                           nfl = 1
                        end if
                        if (iy.eq.0) then
                           gy = c(2,k)
                        else if (dabs(c(2,k)).lt.1.0d-5) then
                           gy = 0.0d0
                        else
                           gy = -c(2,k)
                           nfl = 1
                        end if
                        if (iz.eq.0) then
                           gz = c(3,k)
                        else if (dabs(c(3,k)).lt.1.0d-5) then
                           gz = 0.0d0
                        else
                           gz = -c(3,k)
                           nfl = 1
                        end if
                        if (nfl.eq.1) then
                           k1 = k + 1
                           if (k1.le.nat) then
                              iex = idx
                              do 410 l = k1 , nat
                                 iex = iex + 1
                                 if (natt(l).eq.0 .and. ip(iex).eq.0)
     +                               then
                                    if (iord(l).eq.j) then
                                       if (dabs(gx-c(1,l)).le.1.0d-5
     +)                                    then
                                         if (dabs(gy-c(2,l))
     +                                      .le.1.0d-5) then
                                         if (dabs(gz-c(3,l))
     +                                      .le.1.0d-5) then
                                         ip(idx) = l
                                         ip(iex) = k
                                         isl = isl + 2
                                         if (isl.ne.isal) go to 420
                                         go to 430
                                         end if
                                         end if
                                       end if
                                    end if
                                 end if
 410                          continue
                           end if
                           ise(i) = 0
                           go to 500
                        else
                           ip(idx) = k
                           isl = isl + 1
                           if (isl.eq.isal) go to 430
                        end if
                     end if
                  end if
 420           continue
            end if
 430     continue
         if (ibs.ne.0) then
            do iatt = 1, 3
            do 460 j = 1 , nbas
               if (ibas(j).eq.iatt) then
                  isal = nord(j)
                  isl = 0
                  idx = icx
                  do 450 k = 1 , nat
                     idx = idx + 1
                     if (natt(k).eq.iatt .and. ip(idx).eq.0) then
                        if (iord(k).eq.j) then
                           nfl = 0
                           if (ix.eq.0) then
                              gx = c(1,k)
                           else if (dabs(c(1,k)).lt.1.0d-5) then
                              gx = 0.0d0
                           else
                              gx = -c(1,k)
                              nfl = 1
                           end if
                           if (iy.eq.0) then
                              gy = c(2,k)
                           else if (dabs(c(2,k)).lt.1.0d-5) then
                              gy = 0.0d0
                           else
                              gy = -c(2,k)
                              nfl = 1
                           end if
                           if (iz.eq.0) then
                              gz = c(3,k)
                           else if (dabs(c(3,k)).lt.1.0d-5) then
                              gz = 0.0d0
                           else
                              gz = -c(3,k)
                              nfl = 1
                           end if
                           if (nfl.eq.1) then
                              k1 = k + 1
                              if (k1.le.nat) then
                                 iex = idx
                                 do 440 l = k1 , nat
                                    iex = iex + 1
                                    if (natt(l).eq.iatt .and. 
     &                                  ip(iex).eq.0) then
                                       if (iord(l).eq.j) then
                                         if (dabs(gx-c(1,l))
     +                                      .le.1.0d-5) then
                                         if (dabs(gy-c(2,l))
     +                                      .le.1.0d-5) then
                                         if (dabs(gz-c(3,l))
     +                                      .le.1.0d-5) then
                                         ip(idx) = l
                                         ip(iex) = k
                                         isl = isl + 2
                                         if (isl.ne.isal) go to 450
                                         go to 460
                                         end if
                                         end if
                                         end if
                                       end if
                                    end if
 440                             continue
                              end if
                              ise(i) = 0
                              go to 500
                           else
                              ip(idx) = k
                              isl = isl + 1
                              if (isl.eq.isal) go to 460
                           end if
                        end if
                     end if
 450              continue
               end if
 460        continue
            enddo ! iatt
         end if
         ise(i) = 1
 500  continue
c
c array ise contains 1 if that operation belongs to group, 0 otherwise
c
c begin AdM
c     get group order
      maxise = nise(ise)
1     if (nise(ise).gt.maxise) then
c (maxise.eq.2 or 4) reducing to C2v or C2
         if (maxise.eq.4) then
            ise(3)=0
            ise(5)=0
            ise(6)=0
            ise(7)=0
            write(iwr,6130)txtsym(3)
         elseif (maxise.eq.2) then
            ise(1)=0
            ise(2)=0
            ise(3)=0
            ise(7)=0
            ise(5)=0
            ise(6)=0
            write(iwr,6130)txtsym(2)
         elseif (maxise.eq.1) then
            ise(1)=0
            ise(2)=0
            ise(3)=0
            ise(7)=0
            ise(4)=0
            ise(5)=0
            ise(6)=0
            write(iwr,6130)txtsym(1)
         else
            call caserr2('Unable to reduce symmetry in adapt')
        endif
      endif
c end AdM
      icx = -namx
      ix = 0
      do 520 i = 1 , 7
         icx = icx + namx
         if (ise(i).ne.0) then
            iy = 0
            ix = ix + 1
            idx = icx
            do 510 j = 1 , nat
               idx = idx + 1
               if (ip(idx).ne.j) then
                  iy = iy + ncont(j)
               end if
 510        continue
            npf(i) = iy
         else
       npf(i)=0
c probably not necessary, but lets initialize all elements
         end if
 520  continue
      ifound=0
      oisd2=.false.
      if(ix.eq.0) then
      write(iwr,6160)txtsym(1)
      ifound=ifound+1
      endif
      if(ix.eq.1) then
      if(ise(1).eq.1.or.ise(2).eq.1.or.ise(3).eq.1) then
      write(iwr,6160)txtsym(4)
      ifound=ifound+1
      endif
      if(ise(4).eq.1.or.ise(5).eq.1.or.ise(6).eq.1) then
      write(iwr,6160)txtsym(2)
      ifound=ifound+1
      endif
      if(ise(7).eq.1)then
      write(iwr,6160)txtsym(5)
      ifound=ifound+1
      endif
      endif
      if(ix.eq.7)then
      write(iwr,6160)txtsym(6)
      ifound=ifound+1
      endif
      if(ix.eq.3) then
      if(ise(4)+ise(5)+ise(6).eq.3) then
      write(iwr,6160)txtsym(7)
      ifound=ifound+1
      oisd2=.true.
      endif
      if(ise(4)+ise(5)+ise(6).eq.1 .and. ise(1)+ise(2)+ise(3).eq.2)then
      write(iwr,6160)txtsym(3)
      ifound=ifound+1
      endif
      if(ise(4)+ise(5)+ise(6).eq.1 .and. ise(1)+ise(2)+ise(3).eq.1 
     +                             .and. ise(7).eq.1 )then
      write(iwr,6160)txtsym(8)
      ifound=ifound+1
      endif
      endif
      if(ifound.ne.1) then
      write(iwr,6160) txtsym(9)
      goto 1100
      endif
c 
c jx will be number of generators of given group
c
      if (ix.gt.0) then
c
c
         if (ix.eq.1) jx = 1
         if (ix.eq.3) jx = 2
         if (ix.eq.7) jx = 3
         nsel = ix + 1
         kx = 0
         do 530 i = 1 , 7
            npz(i) = 0
 530     continue
c
c the following selects the generators
c
         do 560 j = 1 , jx
            lx = 10000
            do 540 i = 1 , 7
               if (ise(i).ne.0) then
                  if (npz(i).ne.1) then
                     nf = npf(i)
                     if (nf.lt.lx) then
                        lx = nf
                        ji = i
                        if (nf.eq.0) go to 550
                     end if
                  end if
               end if
 540        continue
            npz(ji) = 1
            nel(j) = ji
            go to 560
 550        kx = kx + 1
            npz(ji) = 1
            nel(j) = ji
 560     continue
         if (jx.ge.3) then
            ia = nel(1)
            ib = nel(2)
            kk = nel(3)
            igx = isx(ia) + isx(ib) + isx(kk)
            if (igx.gt.1) igx = igx - 2
            if (igx.eq.0) then
               igx = isy(ia) + isy(ib) + isy(kk)
               if (igx.gt.1) igx = igx - 2
               if (igx.eq.0) then
                  igx = isz(ia) + isz(ib) + isz(kk)
                  if (igx.gt.1) igx = igx - 2
                  if (igx.eq.0) then
                     lx = 10000
                     if (npf(kk).eq.0) kx = kx - 1
                     do 570 i = 1 , 7
                        if (npz(i).ne.1) then
                           nel(3) = i
                           go to 580
                        end if
 570                 continue
 580                 if (npf(i).eq.0) kx = kx + 1
                  end if
               end if
            end if
         end if
         nbon(1) = 0
         if (nat.ne.1) then
            do 590 i = 2 , nat
               il = i - 1
               nbon(i) = nbon(il) + ncont(il)
 590        continue
         end if
         icx = -iorbs
c
c the following calculates transformation table of
c shells nfun and 'parity' array npar
c nfun((op-1)*nshells+ishell)=+/-result_shell(op)
c npar is -1 ->to another shell +/-, 0 to itsself, 1 to -itsself
c
c        loop for the generators
         do 640 i = 1 , jx
            icx = icx + iorbs
            ji = nel(i)
            idx = (ji-1)*namx
c           does operation invert that coordinate
            ix = isx(ji)
            iy = isy(ji)
            iz = isz(ji)
            call setsto(nat,0,ibas)
            ican = 0
c           loop over distinct atoms
            do 630 j = 1 , nbas
               isal = nord(j)
               ncw = nco(j)
               if (ncw.ne.0) then
                  icbn = ican
                  ist = 0
c             now investigate all atoms of that distinct type
                  do 620 k = 1 , nat
                     if (iord(k).eq.j) then
                        if (ibas(k).eq.0) then
c                          counter of orbs at previous atoms
                           ipx = idx + k
                           ipx = ip(ipx)
                           ka = nbon(k)
                           nbb = ka + icx
                           ican = icbn
                           if (ipx.ne.k) then
c                             has a counterpart
                              im = nbon(ipx)
                              mb = im + icx
c                        that shell has been found/processed
                              ibas(ipx) = 1
                              ist = ist + 2
                              do 600 l = 1 , ncw
                                 nbb = nbb + 1
                                 ka = ka + 1
                                 im = im + 1
                                 mb = mb + 1
                                 npar(nbb) = -1
                                 npar(mb) = -1
                                 ican = ican + 1
                                 ittp = ityp(ican)
c                           count parity with respect to the generator
                                 isum = 0
                                 if (ix.eq.1) isum = isum + nix(ittp)
                                 if (iy.eq.1) isum = isum + niy(ittp)
                                 if (iz.eq.1) isum = isum + niz(ittp)
                                 if (isum.eq.(isum/2)*2) then
                                    nfun(nbb) = im
                                    nfun(mb) = ka
                                 else
                                    nfun(nbb) = -im
                                    nfun(mb) = -ka
                                 end if
 600                          continue
                              if (ist.eq.isal) go to 630
                           else
c                             transforms to itsself
                              ist = ist + 1
                              do 610 l = 1 , ncw
                                 ican = ican + 1
                                 nbb = nbb + 1
                                 ka = ka + 1
                                 ittp = ityp(ican)
                                 isum = 0
                                 if (ix.eq.1) isum = isum + nix(ittp)
                                 if (iy.eq.1) isum = isum + niy(ittp)
                                 if (iz.eq.1) isum = isum + niz(ittp)
                                 if (isum.eq.(isum/2)*2) then
                                    nfun(nbb) = ka
                                    npar(nbb) = 0
                                 else
                                    nfun(nbb) = -ka
                                    npar(nbb) = 1
                                 end if
 610                          continue
                              if (ist.eq.isal) go to 630
                           end if
                        end if
                     end if
 620              continue
               end if
 630        continue
 640     continue
         if (kx.eq.jx) then
            do 650 i = 1 , nsel
               mj(i) = 0
 650        continue
            kx = jx
            ly = -iorbs
            do 670 i = 1 , iorbs
               ly = ly + 1
               my = ly
               isum = 1
               ny = 1
               do 660 j = 1 , jx
                  my = my + iorbs
                  lz = npar(my)
                  if (lz.eq.1) isum = isum + ny
                  ny = ny + ny
 660           continue
               nc = mj(isum) + 1
               mj(isum) = nc
               nir(i) = isum
               loc(i) = nc
 670        continue
            nzer(1) = 0
            do 680 i = 2 , nsel
               i9 = i - 1
               nzer(i) = nzer(i9) + mj(i9)
 680        continue
            jv = 0
            do 690 i = 1 , nsel
               js = mj(i)
               nj(i) = js
               ntil(i) = jv
               nbal(i) = jv
               jv = jv + js
 690        continue
            do 700 i = 1 , iorbs
               ni = nir(i)
               nz = nzer(ni) + 1
               nzer(ni) = nz
               lsym(nz) = i
 700        continue
            call setsto(iorbs,1,ncomp)
            ix = iorbs
         else
            ix = 0
            iy = 0
            lx = 0
            call setsto(iorbs,0,nfil)
            ly = -iorbs
c
c Here is the main part which generates adapted functions
c main input is nfun and npar, output ksym, mpar, mcomp,
c ksym are the linear combinations; mcomp are numbers of
c components of sabf (lengths of records in ksym)
c and mpar contains parity with respect to group generators,
c from which the irrep assignment is straightforward.
c The whole code is branched according to the number of generators,
c which map a given shell to another shell resp. itself
c conf. ary ibuk resp. inuk
c usually the number of components of given sabf is
c 2**(no. of generators mapping it to other shell),  but it does 
c not always have to be so, because two generators can 
c map one shell to +/- another shell in the D2 case
c
            do 890 i = 1 , iorbs
               ly = ly + 1
               if (nfil(i).eq.0) then
                  my = ly
                  do 710 j = 1 , jx
                     my = my + iorbs
                     imix(j) = npar(my)
 710              continue
                  iz = 0
                  jz = 0
                  do 720 j = 1 , jx
                     if (imix(j).ge.0) then
                        jz = jz + 1
                        inuk(jz) = j
                     else
                        iz = iz + 1
                        ibuk(iz) = j
                     end if
 720              continue
                  if (iz.gt.0) then
c                    at least one counterpart
                     ih = ix
                     ix = ix + 1
                     ksym(ix) = i
                     k = ibuk(1)
                     ig = iorbs*(k-1) + i
                     ix = ix + 1
                     mq = nfun(ig)
                     ksym(ix) = mq
                     mqsav=mq
                     if (mq.lt.0) mq = -mq
                     nfil(mq) = 1
                     if (iz.eq.1) then
c                       one counterpart code
                        ksym(ix+1) = ksym(ix-1)
                        ksym(ix+2) = -ksym(ix)
                        ix = ix + 2
                        do 730 j = 1 , 2
                           lx = lx + 1
                           mcomp(lx) = 2
 730                    continue
c                       branch according to the no. of 
c                       self mapping generators
                        if (jz.ge.2) then
                           lq = inuk(1)
                           mq = imix(lq)
                           ky = iy + lq - 3
                           do 740 j = 1 , 2
                              ky = ky + 3
                              mpar(ky) = mq
 740                       continue
                           lq = inuk(2)
                           mq = imix(lq)
                           ky = iy + lq - 3
                           do 750 j = 1 , 2
                              ky = ky + 3
                              mpar(ky) = mq
 750                       continue
                           lq = ibuk(1) + iy
                           mpar(lq) = 0
                           mpar(lq+3) = 1
                           iy = iy + 6
                        else if (jz.eq.0) then
                           mpar(iy+1) = 0
                           mpar(iy+2) = 1
                           iy = iy + 2
                        else
                           lq = inuk(1)
                           mq = imix(lq)
                           ky = iy + lq - 2
                           do 760 j = 1 , 2
                              ky = ky + 2
                              mpar(ky) = mq
 760                       continue
                           lq = ibuk(1) + iy
                           mpar(lq) = 0
                           mpar(lq+2) = 1
                           iy = iy + 4
                        end if
                     else
c                         2 or more counterparts code
                        jg = iorbs*(ibuk(2)-1)
                        ix = ix + 1
                        itmp = nfun(jg+i)
                      jtmp=itmp
c check if it not equal to mq from first operation
                        if(itmp.lt.0) itmp= -itmp
                      if(mq.eq.itmp) then
                         if(.not. oisd2) then
                            write(iwr,6190)
                            goto 1100
                           endif
c  this section here should fix the special case of td->d2 
c  octahedron+tetrahedron
c  and probably also some other rare cases
c  Jiri Pittner 1993
c                          first do the same as if iz.eq.1
c                           l.c. with opposite sign
c                          compensate the previous increment
                           ix=ix-1
                         ksym(ix+1) = ksym(ix-1)
                         ksym(ix+2) = -ksym(ix)
                         ix = ix + 2
c                           set up the component counter
                         do 735 j = 1 , 2
                            lx = lx + 1
                            mcomp(lx) = 2
735                        continue
c  remains to properly set up IR number
c  with respect to the first operation
                          mpar(iy+ibuk(1))=0
                          mpar(iy+ibuk(1)+2)=1
c  with respect to the second operation
                            if(jtmp.eq.mqsav) then
                            mpar(iy+ibuk(2))=0
                            mpar(iy+ibuk(2)+2)=1
                          else
                            mpar(iy+ibuk(2))=1
                            mpar(iy+ibuk(2)+2)=0
                          endif
                            iy=iy+4
c                           exit the current iz/jz branching
                         goto 885
                      endif
                        mq = nfun(jg+i)
                        ksym(ix) = mq
                        if (mq.lt.0) mq = -mq
                        nfil(mq) = 1
                        ix = ix + 1
                        ig = ksym(ih+2)
                        if (ig.lt.0) then
                           mq = -nfun(jg-ig)
                           ksym(ix) = mq
                           if (mq.lt.0) mq = -mq
                           nfil(mq) = 1
                        else
                           mq = nfun(jg+ig)
                           ksym(ix) = mq
                           if (mq.lt.0) mq = -mq
                           nfil(mq) = 1
                        end if
                        if (iz.eq.2) then
                           do 770 j = 1 , 4
                              lx = lx + 1
                              mcomp(lx) = 4
 770                       continue
                           ix = ix - 3
                           do 780 j = 1 , 3
                              ix = ix + 4
                              ksym(ix) = i
 780                       continue
                           nq = 0
                           ih = ih + 1
                           do 800 j = 1 , 3
                              ih = ih + 1
                              mq = ksym(ih)
                              kh = ih
                              do 790 k = 1 , 3
                                 nq = nq + 1
                                 kh = kh + 4
                                 if (jdh(nq).lt.0) then
                                    ksym(kh) = -mq
                                 else
                                    ksym(kh) = mq
                                 end if
 790                          continue
 800                       continue
                           ix = kh
                           if (jz.eq.0) then
                              mpar(iy+1) = 0
                              mpar(iy+2) = 0
                              mpar(iy+3) = 1
                              mpar(iy+4) = 0
                              mpar(iy+5) = 0
                              mpar(iy+6) = 1
                              mpar(iy+7) = 1
                              mpar(iy+8) = 1
                              iy = iy + 8
                           else
                              lq = inuk(1)
                              mq = imix(lq)
                              ky = iy + lq - 3
                              do 810 j = 1 , 4
                                 ky = ky + 3
                                 mpar(ky) = mq
 810                          continue
                              ky = iy + ibuk(1)
                              mpar(ky) = 0
                              mpar(ky+6) = 0
                              mpar(ky+3) = 1
                              mpar(ky+9) = 1
                              ky = iy + ibuk(2)
                              mpar(ky) = 0
                              mpar(ky+3) = 0
                              mpar(ky+6) = 1
                              mpar(ky+9) = 1
                              iy = iy + 12
                           end if
                        else
c                         3 counterparts
                           jg = iorbs + iorbs
                           ix = ix + 1
                           mq = nfun(jg+i)
                           ksym(ix) = mq
                           if (mq.lt.0) mq = -mq
                           nfil(mq) = 1
                           do 820 j = 1 , 3
                              ix = ix + 1
C AdM/JvL
                              ig = ksym(ih+j+1)
                              if (ig.lt.0) then
                                 mq = -nfun(jg-ig)
                                 ksym(ix) = mq
                                 if (mq.lt.0) mq = -mq
                              else
                                 mq = nfun(jg+ig)
                                 ksym(ix) = mq
                                 if (mq.lt.0) mq = -mq
                              end if
                              nfil(mq) = 1
 820                       continue
                           do 830 j = 1 , 8
                              lx = lx + 1
                              mcomp(lx) = 8
 830                       continue
                           do 840 j = 1 , 24
                              iy = iy + 1
                              mpar(iy) = irrep(j)
 840                       continue
                           ix = ix - 7
                           do 850 j = 1 , 7
                              ix = ix + 8
                              ksym(ix) = i
 850                       continue
                           nq = 0
                           ih = ih + 1
                           do 870 j = 1 , 7
                              ih = ih + 1
                              mq = ksym(ih)
                              kh = ih
                              do 860 k = 1 , 7
                                 nq = nq + 1
                                 kh = kh + 8
                                 if (idh(nq).lt.0) then
                                    ksym(kh) = -mq
                                 else
                                    ksym(kh) = mq
                                 end if
 860                          continue
 870                       continue
                           ix = kh
                        end if
c                             2 counterparts
                     end if
c                             1 counterpart
                  else
c                    no counterpart
                     lx = lx + 1
                     mcomp(lx) = 1
                     ix = ix + 1
                     ksym(ix) = i
                     do 880 j = 1 , jx
                        iy = iy + 1
                        mpar(iy) = imix(j)
 880                 continue
                  end if
c               this is label for the special case of D2
 885              continue
c                    not processed shells if
               end if
 890        continue
c
c the main job has been done; it remains to reorder the information
c reduce symmetry, if necessary, calc irrep numbers, and create
c arrays to be used in tdown etc.
c
        if ( lx .ne. iorbs ) then
         write(iwr,6170) iorbs, lx
         goto 1100
         endif
         ixsum=0
         do 895, i=1,lx
895        ixsum=ixsum+mcomp(i)
         if(ix.ne.ixsum) then
         write(iwr,6180) ix, ixsum 
         goto 1100
         endif
c                          group order
            do 900 i = 1 , nsel
               mj(i) = 0
 900        continue
            if (kx.ne.0) then
               ly = -iorbs
               do 920 i = 1 , iorbs
                  ly = ly + 1
                  my = ly
                  isum = 1
                  ny = 1
                  do 910 j = 1 , kx
                     my = my + iorbs
                     lg = npar(my)
                     if (lg.eq.1) isum = isum + ny
                     ny = ny + ny
 910              continue
                  nc = mj(isum) + 1
                  mj(isum) = nc
                  nir(i) = isum
                  loc(i) = nc
 920           continue
            end if
            do 930 i = 1 , nsel
               nj(i) = 0
               nblz(i) = 0
 930        continue
            i9 = 0
            do 950 i = 1 , iorbs
               ny = 1
               isum = 1
               do 940 j = 1 , jx
                  i9 = i9 + 1
                  if (mpar(i9).eq.1) isum = isum + ny
                  ny = ny + ny
 940           continue
               nsog(i) = isum
               nj(isum) = nj(isum) + 1
               nblz(isum) = nblz(isum) + mcomp(i)
 950        continue
            nzil(1) = 0
            ncal(1) = 0
            do 960 i = 2 , nsel
               i9 = i - 1
               ncal(i) = ncal(i9) + nblz(i9)
               nzil(i) = nzil(i9) + nj(i9)
 960        continue
            nork = ncal(nsel) + nblz(nsel)
            do 970 i = 1 , nsel
               ntil(i) = nzil(i)
               nbal(i) = ncal(i)
 970        continue
            lz = 0
            nbal(nsel+1) = nork
            do 990 i = 1 , iorbs
               isum = nsog(i)
               mc = mcomp(i)
               nz = nzil(isum) + 1
               nt = ncal(isum)
               ncomp(nz) = mc
               nzil(isum) = nz
               do 980 j = 1 , mc
                  lz = lz + 1
                  nt = nt + 1
                  lsym(nt) = ksym(lz)
 980           continue
               ncal(isum) = nt
 990        continue
         end if
      else
         write (iwr,6020)
c
         nsel = 1
         nj(1) = num
         do loop = 1 , num
            ncomp(loop) = 1
         enddo
         do loop = 1 , num
            lsym(loop) = loop
         enddo
      end if
c
 1010 otran = .false.
c
      if (oprint(46)) write (iwr,6030)
      n = 0
      jum = 0
      newbas = 0
      do 1040 i = 1 , nsel
         if (oprint(46)) write (iwr,6040) i
         njx = nj(i)
         if (njx.eq.0) then
            if (oprint(46)) write (iwr,6080)
         else
            n = n + 1
            if (oprint(46)) write (iwr,6050)
            do 1030 j = 1 , njx
               newbas = newbas + 1
               ilifc(newbas) = jum
               numb = ncomp(newbas)
               ntran(newbas) = numb
C AdM/JvL begin
               if (numb.lt.1) goto 1100
               if ((jum+numb).gt.mxorb3) then
                  if(idadap.eq.1)call caserr2('Array overflow in adapt')
                  maxise=nise(ise) / 2
                  if (maxise.eq.0) then
                    call caserr2('Impossible error in adapt')
                  endif 
                  goto 1
               endif
c AdM/JvL end
c
c create itran and ctran
c
               do 1020 kk = 1 , numb
                  jum = jum + 1
                  kkkk = lsym(jum)
                  cc = 1.0d0
                  if (kkkk.lt.0) then
                     cc = -cc
                     l = new(-kkkk)
                     iwrit(kk) = -l
                  else
                     l = new(kkkk)
                     iwrit(kk) = l
                  end if
                  ctran(jum) = cc
                  itran(jum) = l
 1020          continue
c
               if (oprint(46)) write (iwr,6060) newbas , numb ,
     +                                (iwrit(ll),ll=1,numb)
 1030       continue
            if (oprint(46)) write (iwr,6070) njx
         end if
 1040 continue
c
      if (newbas.ne.num) go to 1100
c
c ----- load up symmetry section 190
c
c     nsel is the highest symmetry encountered
c     set nirr correctly to corresponding group dimension
c
      if (nsel.gt.4) then
         nirr = 8
      else if (nsel.eq.3.or.nsel.eq.4) then
         nirr = 4
      else
         nirr = nsel
      end if
c
      itt = 0
      do 1070 i = 1 , 8
         do 1060 j = 1 , i
            itt = itt + 1
            mult(i,j) = jsym(itt)
         mult(j,i) = jsym(itt)
 1060    continue
 1070 continue
      itt = 0
      do 1090 i = 1 , nsel
         njx = nj(i)
         if (njx.ne.0) then
            do 1080 j = 1 , njx
               itt = itt + 1
               isymao(itt) = i
 1080       continue
         end if
 1090 continue
      i = maxorb - itt
      if (i.gt.0) call setsto(i,0,isymao(itt+1))

c...     HARMONIC OPTION
c
      if (oharm) call harmonic
c
      call secput(isect(490),m51,lensec(mach(13)),j)
      nav = lenwrd()
      call wrt3i(nirr,mach(13)*nav,j,idaf)
c
      if (oprint(46)) then
         write (iwr,6090)
      else
         if (oharm) then
            write (iwr,6101)
         else
            write (iwr,6100)
         end if
         do 1050 loop = 1 , nsel
            njx = nj(loop)
            if (njx.gt.0) then
               if (oharm) then
                  write (iwr,6110) loop,njx,nsym0(loop)
               else
                  write (iwr,6110) loop , njx
               end if
            end if
 1050    continue
         if (oharm) then
            write(iwr,6121)
         else
            write (iwr,6120)
         end if
      end if
c
      return
c
 1100 call caserr2('error detected in symmetry adaption algorithm')
      return
 6010 format (/40x,17('*')/40x,'symmetry adaption'/40x,17('*')/)
 6020 format (/1x,60('*')
     +        /' no symmetry elements have been recognised in the '/
     +        ' current coordinate system ..... '/1x,60('*')//)
 6030 format (40x,'symmetry adapted basis functions'/40x,32('-')/)
 6040 format (/' irreducible representation no. ',i2/1x,33('-'))
 6050 format (/' sequence     no. of    lcbf'/
     +        ' no. of sabf  terms    (gamess numbering)'/)
 6060 format (4x,i4,8x,i2,7x,8i5)
 6070 format (/' *** total no. of sabfs = ',i3)
 6080 format (/' **** there are no sabf of this representation ****'/)
 6090 format (//1x,104('-'))
 6100 format (/1x,30('=')/1x,'irrep  no. of symmetry adapted'/1x,
     +        '       basis functions'/1x,30('='))
 6101 format (/1x,50('=')/1x,'irrep  no. of symmetry adapted'
     +        ,'  no. harmonic s.a.',/1x,
     +        '       basis functions          basis functions'
     +        ,/1x,50('='))
 6110 format (1x,i3,i12,13x,i12)
 6120 format (1x,30('=')/)
 6121 format (1x,50('=')/)
c AdM/JvL
 6130 format (1x,'Array overflow in adapt. Symmetry reduced to ',a3)
 6140 format (/1x,'number of primitive orbitals exceeds mxgaus',2i5)
 6150 format (/1x,'inconsistency in number of contracted orbitals:',
     + 2i5)
 6160 format (/1x,'point group used for symmetry adaption: ',a3)
 6170 format (/1x,'adapt: ERROR inconsistent sabf count ',2i4)
 6180 format (/1x,'adapt : ix sum inconsistent',2i4)
 6190 format (/1x,'unexpected case in adapt')
      end
_ENDEXTRACT
      subroutine harmonic
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c
c...   create an ctrans that yields harmonic from cartesians
c...   to normalise (*not done*) - divide by sqrt(cstra) etc.
c...     **cnorm comments all these normalisations out **
c...   d-order xx,yy,zz,xy,xz,yz => 2zz-xx-yy,xx-yy,xy,xz,yz
c...   f-order xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz =>
c...           (3,0) -(2zz-3xx-3yy)z. (3,1) -(4zz-xx-yy)x,-(4zz-xx-yy)y,
c...           (3,2)  z(xx-yy),xyx,   (3,3)  x(xx-3yy),y(3xx-yy)
c...   g-order xxxx,yyyy,zzzz,xxxy,xxxz,yyyx,yyyz,zzzx,zzzy,
c...           xxyy,xxzz,yyzz,xxyz,yyxz,zzxy      =>
c...           (4,0) 3xxxx+6xxyy-24xxzz+3yyyy-24yyzz+8zzzz
c...           (4,1) 3xxxz+3yyxz-4zzzx, 3xxyz+3yyyz-4zzzy
c...           (4,2) xxxx-6xxzz-yyyy+6yyzz, xxxy+yyyx-6zzxy
c...           (4,3) xxxz-3yyxz, 3xxyz-yyyz
c...           (4,4) xxxx-6xxyy+yyyy, xxxy-yyyx
c...   normalising is not our concern at this moment
c...   see routine setlab
c...   !! a function may only be first in contraction once (cf trantran)
c
INCLUDE(common/harmon)
c
c      tran description : ......
c      otran : true : no contraction
c      then num times ... (num in /infoa/)
c      ntran(i) : # terms
c      ilifc(i) : 0 point (reference) in ctran/itran
c      itran/ctran(ilifc(k)+k) position and coefficient of kth term
c
c... the normal symmetry adaption only combines functions of same type
c... therefore a combination of harmonic and symmetry adaption is 
c... possible
c
INCLUDE(common/nshel)
c
      common/junk2/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +             ctran(mxorb3),
     +             ilifct(maxorb),ntrant(maxorb),itrant(mxorb3),
     +             ctrant(mxorb3)
c     dimension istra(1),iptra(3,3)
      dimension idtra(6,5),iftra(10,7),ifn(10)
      dimension igtra(15,9),ign(15)
cnorm      dimension cstra(1),cptra(3),cdtra(5),cftra(7),cgtra(9)
c                 s
c     data istra/ 1/
c                 x  y  z
cnorm      data cstra/1.0d0/
c     data iptra/ 1, 0, 0,
c    2            0, 1, 0,
c    3            0, 0, 1/
cnorm      data cptra/1.0d0,1.0d0,1.0d0/
c                xx yy zz xy xz yz
      data idtra/-1,-1, 2, 0, 0, 0,
     2            1,-1, 0, 0, 0, 0,
     3            0, 0, 0, 1, 0, 0,
     4            0, 0, 0, 0, 1, 0,
     5            0, 0, 0, 0, 0, 1/
cnorm      data cdtra/4.0d0,1.3333333333333333d0,1.0d0,1.0d0,1.0d0/
c
c..   f-contractions derived by Huub van Dam (1998)
c..   contractions for unnormalized cartesians; correct with dsqrt(ifn)
c..   The relevant portion of the normalisation is for a cartesian 
c..    x**L.y**M.z**N : L!!.M!!.N!! (L!!=prod(i=1,L)(2i-1) ) ,
c..    so L=1,2,3,4 : L!!=1,3,15,105 :=> factors in ifn/ign
c
c                 xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz  
      data iftra/  0,  0, -2,  0,  3,  0,  3,  0,  0,  0, 
     2             1,  0,  0,  0,  0,  1,  0, -4,  0,  0, 
     3             0,  1,  0,  1,  0,  0,  0,  0, -4,  0, 
     4             0,  0,  0,  0,  1,  0, -1,  0,  0,  0, 
     5             0,  0,  0,  0,  0,  0,  0,  0,  0,  1, 
     6             1,  0,  0,  0,  0, -3,  0,  0,  0,  0, 
     7             0, -1,  0,  3,  0,  0,  0,  0,  0,  0/
      data ifn /  15, 15, 15,  3,  3,  3,  3,  3,  3,  1/
cnorm      data cftra/60.0d0,40.0,40.0d0,4.0d0,1.0d0,24.0d0,24.0d0/
c
c..   g-contractions derived by Huub van Dam (1998)
c..   contractions for unnormalized cartesians; correct with dsqrt(ign)
c
c                xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy 
c                  xxyy xxzz yyzz xxyz yyxz zzxy
      data igtra/  3,   3,   8,   0,   0,  0,    0,   0,   0,
     1               6, -24, -24,   0,   0,   0,
     2             0,   0,   0,   0,   3,   0,   0,  -4,   0,
     2               0,   0,   0,   0,   3,   0,
     3             0,   0,   0,   0,   0,   0,   3,   0,  -4,
     3               0,   0,   0,   3,   0,   0,
     4             1,  -1,   0,   0,   0,   0,   0,   0,   0,
     4               0,  -6,   6,   0,   0,   0,
     5             0,   0,   0,   1,   0,   1,   0,   0,   0,
     5               0,   0,   0,   0,   0,  -6,
     6             0,   0,   0,   0,   1,   0,   0,   0,   0,
     6               0,   0,   0,   0,  -3,   0,
     7             0,   0,   0,   0,   0,   0,  -1,   0,   0,
     7               0,   0,   0,   3,   0,   0,
     8             1,   1,   0,   0,   0,   0,   0,   0,   0,
     8              -6,   0,   0,   0,   0,   0,
     9             0,   0,   0,   1,   0,  -1,   0,   0,   0,
     9               0,   0,   0,   0,   0,   0/
      data ign/   105, 105, 105,  15,  15,  15,  15,  15,  15,
     1                9,   9,   9,   3,   3,   3/
cnorm      data cgtra/6720d0,168.0d0,168.0d0,336.0d0,84.0d0,24.0d0,24.0d0,
cnorm     1           192.0d0,12.0d0/
c
c...   loop over the shells
c
       ibfh = 0
       numh = 0
       nnh = 0
c
       do ii=1,nshell
          mini = kmin(ii)
          maxi = kmax(ii)
c...    s 1..1  p 2..4  sp 1..4
          if (mini.le.2) then
             do ih =mini,maxi
                numh = numh + 1
                ntran(numh) = 1
                ilifc(numh) = nnh
                itran(nnh+1) = ibfh + 1
                ctran(nnh+1) = 1.0d0
                nnh = nnh + 1
                ibfh = ibfh + 1
             end do
c...    d 5..10
          else if (mini.eq.5) then
             do ih =1,5
                numh = numh + 1
                ntran(numh) = 0
                ilifc(numh) = nnh
                do i = mini,maxi
                   if (idtra(i-4,ih).ne.0) then
                      ntran(numh) = ntran(numh) + 1
                      itran(nnh+ntran(numh)) = ibfh + i-4
                      ctran(nnh+ntran(numh)) = idtra(i-4,ih)
                   end if
                end do 
                if (ih.eq.1) then
c
c...    to make trantran easier 
c...    in the 1st one, get order zz,xx,yy
c...    this makes sure a function is only found once as first in itranh
c
                   temp = ctran(nnh+1)
                   itemp = itran(nnh+1)
                   ctran(nnh+1) = ctran(nnh+3)
                   itran(nnh+1) = itran(nnh+3)
                   ctran(nnh+3) = ctran(nnh+2)
                   itran(nnh+3) = itran(nnh+2)
                   ctran(nnh+2) = temp
                   itran(nnh+2) = itemp
                end if
                nnh = nnh + ntran(numh)
             end do
cnorm
cnorm             do ih=1,5
cnorm                ihh = numh-5+ih
cnorm                do i=1,ntran(ihh)
cnorm                 ctran(ilifc(ihh)+i)=ctran(ilifc(ihh)+i)/sqrt(cdtra(ih))
cnorm                end do
cnorm             end  do
cnorm
             ibfh = ibfh + 6
c...    f 11..20
          else if (mini.eq.11) then
             do ih =1,7
                numh = numh + 1
                ntran(numh) = 0
                ilifc(numh) = nnh
                do i = mini,maxi
                   if (iftra(i-10,ih).ne.0) then
                      ntran(numh) = ntran(numh) + 1
                      itran(nnh+ntran(numh)) = ibfh + i-10
                      ctran(nnh+ntran(numh)) = iftra(i-10,ih) *
     1                                         dsqrt(1.0d0*ifn(i-10))
                   end if
                end do 
                if (ih.gt.5) then
c
c...    to make trantran easier 
c...    in the 6th and the 7th one, we interchange basis order
c...    this makes sure a function is only found once as first in itranh
c
                   temp = ctran(nnh+1)
                   itemp = itran(nnh+1)
                   ctran(nnh+1) = ctran(nnh+2)
                   itran(nnh+1) = itran(nnh+2)
                   ctran(nnh+2) = temp
                   itran(nnh+2) = itemp
                end if
c
                nnh = nnh + ntran(numh)
             end do
cnorm
cnorm             do ih=1,7
cnorm                ihh = numh-7+ih
cnorm                do i=1,ntran(ihh)
cnorm                 ctran(ilifc(ihh)+i)=ctran(ilifc(ihh)+i)/sqrt(cftra(ih))
cnorm                end do
cnorm             end do
cnorm
             ibfh = ibfh + 10
c...    g 21..35
          else if (mini.eq.21) then
             do ih=1,9
                numh = numh + 1
                ntran(numh) = 0
                ilifc(numh) = nnh
                do i = mini,maxi
                   if (igtra(i-20,ih).ne.0) then
                      ntran(numh) = ntran(numh) + 1
                      itran(nnh+ntran(numh)) = ibfh + i-20
                      ctran(nnh+ntran(numh)) = igtra(i-20,ih) *
     1                                         dsqrt(1.0d0*ign(i-20))
                   end if
                end do
c...       interchanging orders
                if (ih.eq.4.or.ih.ge.6) then
                   if (ih.eq.8) then
                      temp = ctran(nnh+1)
                      itemp = itran(nnh+1)
                      ctran(nnh+1) = ctran(nnh+3)
                      itran(nnh+1) = itran(nnh+3)
                      ctran(nnh+3) = temp
                      itran(nnh+3) = itemp
                   else
                      temp = ctran(nnh+1)
                      itemp = itran(nnh+1)
                      ctran(nnh+1) = ctran(nnh+2)
                      itran(nnh+1) = itran(nnh+2)
                      ctran(nnh+2) = temp
                      itran(nnh+2) = itemp
                   end if
                end if
                nnh = nnh + ntran(numh)
             end do
cnorm
cnorm             do ih=1,9
cnorm                ihh = numh-9+ih
cnorm                do i=1,ntran(ihh)
cnorm                 ctran(ilifc(ihh)+i)=ctran(ilifc(ihh)+i)/sqrt(cgtra(ih))
cnorm                end do
cnorm             end do
cnorm
             ibfh = ibfh + 15
          else
             call caserr2('higher than g not implemented in HARMONIC')
          end if
       end do
c
c...   newbash is # orbitals after harmonic
c...   newbas0 may be lowered later due to depencies
c
       if (newbash.ne.numh) call caserr2('harmonic foul up')
c
c...      transfrom the symmetry adaption to symmetry + harmonic
c...      and correct for non-normalized functions
c
       call trantran(newbash)
c
       return
       end
      subroutine preharm
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c
c...
c...   preparations of harmonic (determine newbas0 and set newbas1)
c...
c...  ******************************************************************
c...  ** To find harmonic related changes look for the variables      **
c...  ** newbas0  :  # of harmonic/reduced gaussians (new basis, l0)  **
c...  ** newbash  :  # of harmonic gaussians (bef. dependency anal.)  ** 
c...  ** newbas1  :  # of cartesian gaussians (original basis,num,l1) **
c...  **      J.H. van Lenthe, R.W.A. Havenith, Utrecht 1997-2006     **
c...  ******************************************************************
c
c...  This change is meant to work in conjunction with the ctrans/adapt
c...  If l0 is reduced to eliminate real dependency it is probably
c...  wise to adjust newbas0; And of course to check
c
c...  Module specific remarks ......
c     SCF : allow dependent vectors in qmat => l0 < l1
c           clear "missing" vectors (both orthog and other)
c     GVB : trying to keep all square => reset l0 to l1
c           so newbas0 occurs more often (often instead of l0)
c           
c
INCLUDE(common/sizes)
INCLUDE(common/harmon)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
c
       newbas1 = num
       if (.not.oharm) then
          newbas0 = num
          newbash = num
          return
       end if
c
c...   loop over the shells
c
       newbas0 = 0
       do ii=1,nshell
          mini = kmin(ii)
          if (mini.eq.1) then
             newbas0 = newbas0 + 1
             if (kmax(ii).eq.4) newbas0 = newbas0 + 3
c...    s 1..1  or  sp 1..4
          else if (mini.eq.2) then
c...    p 2..4 
             newbas0 = newbas0 + 3
          else if (mini.eq.5) then
c...    d 5..10
             newbas0 = newbas0 + 5
          else if (mini.eq.11) then
c...    f 11..20
             newbas0 = newbas0 + 7
          else if (mini.eq.21) then
c...    g 21..35
             newbas0 = newbas0 + 9
          else
             call caserr2('UNKNOWN SHELL in preharm')
          end if
       end do
       newbash = newbas0
c
       return
       end
       subroutine trantran(numh)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c      transform the adapt transformation in tran with the competing
c      transformation giving as formal paramters
c      this routine assumes that itranh contains a different first bf 
c      in each contraction
c      functions that are to be eliminated are made zero
c      let tranq sort this (numh is the real number left)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/tran)
INCLUDE(common/gjs)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
INCLUDE(common/scra7)
INCLUDE(common/restrj)
      common/junk2/ilifch(maxorb),ntranh(maxorb),itranh(mxorb3),
     +             ctranh(mxorb3),
     +             ilifct(maxorb),ntrant(maxorb),itrant(mxorb3),
     +             ctrant(mxorb3)
c
      if (otran) then
         do i=1,num
            ctrant(i) = 1.0d0
            itrant(i) = 1
            ntrant(i) = 1
            ilifct(i) = i-1
         end do
         otran = .false.
      end if
c
c...   save old ctrans for newmrd6 to allow production 
c...   of vectors in non-harmonic symmetry adapted basis
c
      if (opark) then
         ibl7ha = ibl7la
         call wrt3i(ilifc,2*maxorb+mxorb3,ibl7ha,num8)
         call wrt3s(ctran,mxorb3,num8)
         ibl7la = iposun(num8)
      end if
c
_IF(cray)
         call dcopy(maxorb,ilifc,1,ilifct,1)
         call dcopy(maxorb,ntran,1,ntrant,1)
         call dcopy(mxorb3,itran,1,itrant,1)
_ELSE
         call icopy(maxorb,ilifc,1,ilifct,1)
         call icopy(maxorb,ntran,1,ntrant,1)
         call icopy(mxorb3,itran,1,itrant,1)
_ENDIF
         call dcopy(mxorb3,ctran,1,ctrant,1)
c
      numm = 1
      ilifc(1) = 0
c
      do is=1,num
        ntran(numm) = 0
        ofound = .false.
        do it =1,ntrant(is)        
          ibf = itrant(ilifct(is)+it)
c  
c...     search for this basis function in itranh array
c
          do ih =1,numh
            if (ibf.eq.itranh(ilifch(ih)+1)) then
              ofound = .true.
              do i=1,ntranh(ih)
                ntran(numm) = ntran(numm) + 1
                itran(ilifc(numm)+ntran(numm)) = itranh(ilifch(ih)+i)
                ctran(ilifc(numm)+ntran(numm)) = ctrant(ilifct(is)+it)*
     1                                           ctranh(ilifch(ih)+i)
              end do
            end if
          end do            
        end do         
        if (.not.ofound) then
c...    not found => create 0-function (cf qmat)
               ntran(numm) = 1
               itran(ilifc(numm)+ntran(numm)) = 1
               ctran(ilifc(numm)+ntran(numm)) = 0.0d0
        end if
        ilifc(numm+1) = ilifc(numm) + ntran(numm)
        numm = numm + 1
      end do         
c
      numm = numm - 1
      if (numm.ne.num) call caserr2('error in trantran ')
c
      if (ilifc(numm+1).gt.mxorb3) 
     1    call caserr2('array overflow in trantran')
c
         if (nirr.eq.0) then
c...        symmetry was not initialised so fake c1
            nirr = 1
            do i=1,num
               isymao(i) = 1
            end do
         end if
c
c...   print adaption (if requested)
c
      if (opharm) then
         kb = 1
         write(iwr,600) numh,num
600      format(' ** HARMONIC option leaves ',i5,' out of',i5)
         do irr=1,nirr
            nn = 0
            do 1 ke=kb,num
               if (ctran(ilifc(ke)+1).ne.0.0d0) nn = nn + 1
1           if (isymao(ke).ne.irr) go to 2
            nn = nn + 1
            ke = num + 1
2           ke = ke - 1
            nn = nn - 1
            write(iwr,601) irr,nn
601         format(/' Irreducible representation ',i3,' of dim.',i5)
            do k=kb,ke
              il = ilifc(k)
              if (ctran(il+1).ne.0.0d0) then
                 write(iwr,602) k,ntran(k)
602              format(' function',i5,' no. of terms ',i5,
     1                  '      AO , Coefficient ..')
                 write(iwr,603) (itran(il+j),ctran(il+j),j=1,ntran(k))
603              format((t7,5(i4,f8.4)))
              end if
            end do
            kb = ke + 1
         end do
      end if
c
c...  count # irreps per symmetry for harmonics
c...  store irrep of deleted ao in ielimh
c
      do i=1,nirr
         nsym0(i) = 0
      end do
      do k=1,num
         if (ctran(ilifc(k)+1).eq.0.0d0) then
            ielimh(k) = isymao(k)
         else
            ielimh(k) = 0 
            nsym0(isymao(k)) = nsym0(isymao(k)) + 1
         end if
      end do
c
      return
      end
c
      subroutine stretch(q,nfrom,nto,ilifq,nvb0)
c
      implicit REAL (a-h,o-z)
c
c...  stretch vectors from nfrom long to nto long
c...  or reversely compress them from nfrom to nto
c...  the smallest number gives the number of vectors
c...  unless nvb0.ne.0 then do nvb0 vectors (for vb)
c
c
      dimension q(*),ilifq(*)
c
      if (nfrom.eq.nto) return
c
      if (nfrom.lt.nto) then
           nf = nfrom
           if (nvb0.ne.0) nf = nvb0
         do i=nf,1,-1
            kfrom = (i-1)*nfrom
            kto = (i-1)*nto
            call vclr(q(kto+nfrom+1),1,nto-nfrom)
            do k=nfrom,1,-1
               q(kto+k) = q(kfrom+k)
            end do
         end do
      else
           nt = nto
           if (nvb0.ne.0) nt = nvb0
         do i=1,nt
            kfrom = (i-1)*nfrom
            kto = (i-1)*nto
            do k=1,nto
               q(kto+k) = q(kfrom+k)
            end do
         end do
      end if
c
      do i=1,nto
         ilifq(i) = (i-1)*nto
      end do
c
      return
      end
      subroutine comharm(q,cmd,ilifq)
c
      implicit REAL (a-h,o-z)
c
c...  compress vectors and ctrans from normal to harmonic
c...  also compress 1-electron triangles (produced by uncompressed ctrans)
c...       then ilifq is supposed to be the symmetry array (cf qmat_symm)
c...  used for CASSCF, MCSCF
c...  cmd : ctrans  : adapt common/tran/
c...        vectors : compact vectors (take out rows)
c...        matrix  : compact matrix (take out columns)
c...        vbiii   : compact just iii vectors (for vb)
c...        tri     : compact triangle and symmetry array
c
INCLUDE(common/sizes)
INCLUDE(common/tran)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/restri)
c
      dimension q(*),ilifq(*)
      character*(*) cmd
c
      if (newbash.eq.newbas1) return
c
      if (cmd.eq.'ctrans') then
c
c...  adapt (smile) common/tran/
c...  need only change ilifc (addressing) and ntran
c
         k = 1
         do i=1,newbas1
            if (ielimh(i).ne.0) then
               nleft = newbas1-(i-k+1)
               do j=k,nleft
                  ilifc(j) = ilifc(j+1)
                  ntran(j) = ntran(j+1)
                  isymao(j) = isymao(j+1)
               end do
            else 
               k = k + 1
            end if
         end do
c
         k=k-1
         if (k.ne.newbash) call caserr2('error in comharm ')
c
c...     clear remaining part of isymao
         do i=newbash+1,newbas1
            isymao(i) = 0
         end do
c...     update symmetry section on disk
         call secput(isect(490),51,lensec(mach(13)),jjb)
         nav = lenwrd()
         call wrt3i(nirr,mach(13)*nav,jjb,idaf)
c
      else if (cmd.eq.'vectors'.or.cmd.eq.'matrix'
     1                         .or.cmd(1:2).eq.'vb') then
c
         if (cmd(1:2).eq.'vb') then
            read(cmd,'(2x,i3)') nvb0
         else
            nvb0 = newbash
         end if
c
c...  adapt the vectors 
c...  vectors need to be closely packed
c
         k = 1
         do i=1,newbas1
            if (ielimh(i).ne.0) then
               nleft = newbas1-(i-k+1)
               do l=1,nvb0
                  do j=k,nleft
                     q(ilifq(l)+j) = q(ilifq(l)+j+1)
                  end do
               end do
            else 
               k = k + 1
            end if
         end do
c
         k=k-1
         if (k.ne.newbash) call caserr2('error in comharm ')
c
         if (cmd.eq.'matrix') then
c
c...        compact the vectors by eliminating zero columns
c
            n = 0
            do i=1,newbas1
               if (ielimh(i).ne.0) then
                  n = n + 1
               else if (n.gt.0) then
                  do j=1,(newbas1-i+1)*newbas1
                     q((i-n-1)*newbas1+j) = q((i-1)*newbas1+j)
                  end do
               end if
            end do
            if (n.ne.newbas1-newbash) call caserr2('comharm matrix ')
         end if
c
         if (cmd(1:2).ne.'vb') nvb0 = 0
         call stretch(q,newbas1,newbash,ilifq,nvb0)
c
      else if (cmd.eq.'tri') then
c
c...     compact triangle + symmetry array
c
         nn = 0
         ij = 0
         nelim = 0
         do i=1,newbas1
            if (ielimh(i).ne.0) then
               nelim = nelim + 1
               nn = nn + i
            else
               ilifq(i-nelim) = ilifq(i)
               do j=1,i
                  nn = nn + 1
                  if (ielimh(j).eq.0) then
                     ij = ij + 1
                     q(ij) = q(nn)
                  end if
               end do
            end if
         end do
         if (nelim.ne.newbas1-newbash)callcaserr('comharmh tri foulup')
         call vclr(q(newbash*(newbash+1)/2+1),1,
     1             newbas1*(newbas1+1)/2-newbash*(newbash+1)/2)
      else
c
         call caserr2('wrong cmd in comharm')
c
      end if
c
      return
      end
      subroutine expharm(q,cmd,ilifq)
c
      implicit REAL (a-h,o-z)
c
c...  expand vectors and/or ctrans from harmonic to normal
c...  cmd : action required ...
c...    ctrans : adapt ctrans
c...    vectors : do vectors for casscf
c...              ilifq is only set here (by stretch)/mcscf
c
INCLUDE(common/sizes)
INCLUDE(common/tran)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/syminf)
INCLUDE(common/restri)
c
      dimension q(*),ilifq(*)
      character*(*) cmd
c    
      if (newbas0.eq.newbas1) return
c
c...  expand common/tran/
c...  need only change ilifc (addressing) and ntran
c...  and adapt the vectors as well
c...   *** reverse comharm ***
c
      if (cmd.eq.'vectors') then
c
         call exphvec(q,ielimh,ilifq,newbash,newbas1)
         return
      else if (cmd.eq.'ctrans') then
         nleft = newbas1 - newbas0
         do i=newbas1,1,-1
            if (ielimh(i).eq.0) then
               ilifc(i) = ilifc(i-nleft)
               ntran(i) = ntran(i-nleft)
               isymao(i) = isymao(i-nleft)
            else
               ilifc(i) = ilifc(i+1) - 1
               ntran(i) = 1
               isymao(i) = ielimh(i)
               nleft = nleft - 1
            end if
         end do
c
         if (nleft.ne.0) call caserr2('error in expharm ')
c
c...     update symmetry section on disk
         call secput(isect(490),51,lensec(mach(13)),j)
         nav = lenwrd()
         call wrt3i(nirr,mach(13)*nav,j,idaf)
c
      else
         call caserr2('unknown command in expharm')
      end if
c
      return
      end
      subroutine mcharm
c
      implicit REAL (a-h,o-z)
c
      dimension nsyma1(8)
c
c....  routine to set a few variables in mcscf commons
c....  called from inpci
c
INCLUDE(common/sizes)
INCLUDE(common/harmon)
INCLUDE(common/syminf)
INCLUDE(common/gjs)
INCLUDE(common/multic)
c
c...  set up bookkeeping
c
         if (nirrr.ne.nirr) call caserr2(
     +      'nirr gjs and multic .. mcharm')
c
         do i=1,nirr
            nsyma1(i) = 0
         end do
         do i=1,newbas1
               nsyma1(isymao(i)) = nsyma1(isymao(i)) + 1
         end do
c
         ns=1
         ks = newbas0 - nfreez
         do i=1,nirr
            do k=ks+1,ks+nsyma1(i)-nsymao(i)
               itype(k) = i
            end do
            ks = ks + nsyma1(i)-nsymao(i)
            nsymao(i) = nsyma1(i)
            nsymm(i)=nsymao(i)
            istart(i)=ns
            ns=ns+nsymao(i)
            mfin(i)=ns-1
         end do
         nbasao=ns-1
         nbasis=nbasao
         nblkq = nbasis*nbasis
c
      return
      end
      subroutine exphvec(q,ielim,ilifq,n0,n1)
c
      implicit REAL (a-h,o-z)
c
c...  expand vectors from harmonic to normal
c
      dimension q(*),ielim(*),ilifq(*)
c    
      if (n0.eq.n1) return
c
      call stretch(q,n0,n1,ilifq,0)
      nleft = n1 - n0
      do i=n1,1,-1
         if (ielim(i).eq.0) then
               do l=1,n0
                  q(ilifq(l)+i) = q(ilifq(l)+i-nleft)
               end do
         else
               do l=1,n0
                  q(ilifq(l)+i) = 0.0d0
               end do
            nleft = nleft - 1
         end if
      end do
c
      if (nleft.ne.0) call caserr2('error in exphvec')
c
c...   clear the non-existing vectors
c
      do i=n0+1,n1
         do l=1,n1
            q(ilifq(i)+l) = 0.0d0
         end do
      end do
c
      return
      end
      subroutine exphvs(q,isym)
c
c...  expand vectors for one symmetry
c...  we assume ctrans is already expanded
c...  ielimh is symmetry ordered
c...  used in mcgrad to expand lagrangian
c
      implicit REAL (a-h,o-z)
      dimension q(*)
INCLUDE(common/sizes)
INCLUDE(common/harmon)
INCLUDE(common/syminf)
      dimension ilifq(maxorb)
c
      kk = 1
      do i=1,isym-1
         kk = kk + nsymao(i)
      end do
      do i=1,nsymao(isym)
         ilifq(i) = (i-1)*nsymao(isym)
      end do
c
      call exphvec(q,ielimh(kk),ilifq,nsym0(isym),nsymao(isym))
c
      return
      end
      subroutine compv(v,l0,l1,q,cmd)
c
c.... complete the vector-set v to obtain a square matrix
c...  dimension of v : l1*l1
c...  # nonzero vectors : l0
c...  cmd  : 'add' : overwrite end of v
c...        'fill' : fill in already zeroed vectors (do not do e,p)
c
      implicit REAL (a-h,o-z)
      dimension v(l1,l1),q(*)
      character*(*) cmd
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      dimension ilifq(maxorb),iky(maxorb)
c
      if (cmd.ne.'add'.and.cmd.ne.'fill') 
     1     call caserr2('illegal cmd in compv')
      if (cmd.eq.'add') then
         do i=l0+1,l1
            do j=1,l1
               v(j,i) = 0.0d0
            end do
         end do
      end if
c
      iccc = igmem_alloc_inf(l1**2,'guess.m','compv','iccc',IGMEM_DEBUG)
      iev  = igmem_alloc_inf(l1**2,'guess.m','compv','iev',IGMEM_DEBUG)
      ieig = igmem_alloc_inf(l1,'guess.m','compv','ieig',IGMEM_DEBUG)
c
      do i=1,l1
         ilifq(i) = (i-1)*l1
         iky(i) = i*(i-1)/2
      end do
c
c...     generate q*qx
c
      call mxmg(v,1,l1,v,l1,1,q(iccc),1,l1,l1,l1,l1)
c
      do i=1,l1
       do j=1,i
         q(iccc-1+i*(i-1)/2+j) = q(iccc-1+(i-1)*l1+j)
       end do
      end do
      call jacobi(q(iccc),iky,l1,q(iev),ilifq,l1,
     1            q(ieig),2,2,1.0d-9)
c
c...     add 0 eigenvectors on the eliminated spots to v
c
      k = 0
      do 10 i=1,l1
         do j=1,l1
            if (v(j,i).ne.0.0d0) go to 10
         end do
         k = k + 1
c...         check eigenvalue (should be 0.0)
         if (q(ieig+k-1).gt.1.0d-10) 
     1      call caserr2('compv  non 0.0 eigenvalue ')
c...         check eigenvector which is to be replaced (should be 0.0)
         do j=1,l1
            v(j,i) = q(iev-1+ilifq(k)+j)
         end do
10    continue
c
      if (q(ieig+k).lt.1.0d-10) call caserr2('compv : excess 0')
      if (k.ne.l1-l0) call caserr2('compv : dim error')
c
      call gmem_free_inf(ieig,'guess.m','compv','ieig')
      call gmem_free_inf(iev,'guess.m','compv','iev')
      call gmem_free_inf(iccc,'guess.m','compv','iccc')
c
      return
      end
**==alphaz.f
      subroutine alphaz(t,f,scra,num,lennew)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/datgue)
INCLUDE(common/iofile)
      dimension t(*),f(*),scra(*)
      write (iwr,6010)
      call vsub(f,1,t,1,t,1,lennew)
      l = 0
      do 30 i = 1 , num
         m = i*(i+1)/2
         top = alphas(i) - f(m)
         test = t(m)
         alphas(i) = top
         scra(i) = test
         do 20 j = 1 , i
            l = l + 1
            f(l) = f(l) + t(l)*(top+alphas(j))/(test+scra(j))
 20      continue
 30   continue
      return
 6010 format (/1x,'initial guess orbitals generated from alphas'/)
      end
**==anorm.f
      subroutine anorm(sovl,q)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/scra7)
INCLUDE(common/infoa)
INCLUDE(common/tran)
INCLUDE(common/iofile)
INCLUDE(common/tranao)
INCLUDE(common/harmon)
      dimension  sovl(*),q(*)
c
      equivalence (num,nbasis)
c
      l2 = nbasis*(nbasis+1)/2
c
c     -----  generate inverse of ctrans for ao-> sabf tran.
c
      nbas2 = nbasis*nbasis
      maxllen = mxorb3 + mxorb3 + 1
      call vclr(sovl,1,nbas2)
      ii = 0
      do 30 i = 1 , num
         n = ntran(i)
         nsp = ilifc(i)
         do 20 j = 1 , n
            nsp = nsp + 1
            sovl(ii+itran(nsp)) = ctran(nsp)
 20      continue
         ii = ii + nbasis
 30   continue
c
c...  we have to complete the ctrans  to allow inversion
c...  this code might be used if we want to generate a 
c...  harmonic + the spherical components ctrans
c 
      if (oharm) call compv(sovl,newbas0,newbas1,q,'fill')
c
_IFN1(c)      tester=0.0d0
_IFN1(c)      call minvrt(sovl,nbasis,tester,loc,ncomp)
_IFN1(c)      if( dabs(tester).ge.1.0d-6)go to 50
_IF1(c)      call minv(sovl,nbasis,nbasis,loc,tester,1.0d-30,0,1)
_IF1(c)      if( dabs(tester).ge.1.0d-20)go to 50
      write (iwr,6010) tester
 40   call caserr2('error detected in symmetry adaption algorithm')
c
c    ----- now build up inverse ctrans in /tranao/
c
c
50    jum = 0
      ii = 0
      do 70 i = 1 , nbasis
         ilifd(i) = jum
         numb = 0
         do 60 j = 1 , num
            if (dabs(sovl(ii+j)).gt.1.0d-6) then
               numb = numb + 1
               jum = jum + 1
               ctrad(jum) = sovl(ii+j)
               itrad(jum) = j
            end if
 60      continue
         if (numb.lt.1 .or. jum.ge.maxllen) go to 40
         ntrad(i) = numb
         ii = ii + num
 70   continue
c
      return
 6010 format (//' ***** error *****'//
     +        ' ctrans matrix is singular : det = ',e20.13)
      end
**==aolim2.f
      subroutine aolim2
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
      common/junk/limlow(maxat),limsup(maxat)
      do 20 i = 1 , nat
         limlow(i) = 0
         limsup(i) = 0
 20   continue
      lat = katom(1)
      limlow(lat) = 1
      do 30 i = 1 , nshell
         iat = katom(i)
         if (lat.ne.iat) then
            limsup(lat) = kloc(i) - 1
            lat = iat
            limlow(lat) = kloc(i)
         end if
 30   continue
      limsup(iat) = num
      return
      end
**==atcond.f
      subroutine atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps,iecpt,og,
     1                  oprin,istate,iatdiff)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     for atom of nuclear charge zn find electron configuration
c     and k(l,l,0) coupling coefficient. simple aufbau is
c     assumed throughout the periodic system. the algorithm used
c     will work for zn less than 119.
c
c     to allow for pseudo-potentials the effective charge znps
c     is used to determine the number of closed shells filled
c     jvl  (daresbury 1988): corrected July 93
c     added possibility of exact groundstate for up to d-openshells
c.......................................................................
      dimension ncsh(*),nccup(*),nosh(*),ajmn(*)
      logical og,oprin
c
INCLUDE(common/sizes)
INCLUDE(common/datgue)
c
c...  Begin ECP-data
c
c     dimension nelecp(26), iecp(26)
      dimension isy(17,3)
c
c...  Different schemes to fill the core orbitals in the ECP
c...  (1=s, 2=p, 3=d, 4=f)
c
c....      LANL isymhw (version 6.2)
      data (isy(i,1),i=1,17)/ 1, 1,2, 1,2,3, 1,2,3 ,1,2,4,3, 1,2,4,3/
c....      CEP / SBKJC isymce
      data (isy(i,2),i=1,17)/ 1, 1,2, 1,2,3, 1,2,3,4, 1,2,3,4, 1,2,3/
c...       isymab
      data (isy(i,3),i=1,17)/ 1, 1,2, 1,3,2, 1,3,2,1, 4,3,2,1, 4,3,2/
c
c...  nelecp(i) gives the number of electrons in ECP i.
c...  iecp(i)   gives a filling scheme for ECP i.
c
c     data nelecp/ 
c    +     0,2,4,10,12,18,22,28,30,36,40,46,48,54,60,62,68,72,78,80,
c    +     86,92,94,100,104,110/
c     data iecp/
c    +     1,1,1,1 ,1 ,1 ,3 ,1 ,1 ,1 ,3 ,1 ,1 ,1 ,2 ,2 ,1 ,3 ,1 ,1 ,
c    +     1 ,2 ,2 ,1  ,3  ,1  /
c
c...  End ECP-Data
c.......................................................................
c
c     initialize .
c.......................................................................
      nsym = 0
      maxbas = 4
      do 20 i = 1 , maxbas+1
         ncsh(i) = 0
         nosh(i) = 0
         nccup(i) = 0
 20   continue
c.......................................................................
cmarcin      print *,'nz start: ',zn
      nz = zn + 0.1d0
      nzps = znps + 0.1d0
c
      idiff = nz - nzps
      ihw = 0
      if (idiff.ne.0) then
c
c        we have a pseudopotential
c
         ihw = iecpt
      endif
      nzps = idiff
c.......................................................................
c
c     fill up shells - all electron first
c.......................................................................
      nlast = 0
      do 50 i = 1 , maxbas
         do 40 j = 1 , 2
            ksym = i
            do 30 k = 1 , i
               nshell = nlast + 4*ksym - 2
               if (nz.lt.nshell) go to 60
               nsym = max(nsym,ksym)
               ncsh(ksym) = ncsh(ksym) + 1
               ksym = ksym - 1
               nlast = nshell
 30         continue
 40      continue
 50   continue
c.......................................................................
c
c     check if open shell atom. test for la and ac.
c.......................................................................
 60   if (nz.eq.57 .or. nz.eq.89) then
         ncsh(4) = ncsh(4) - ncsh(4)/2
         ksym = 3
      else if (nz.eq.79) then
c..      Au d10s1 instead of d9s2
         ncsh(3) = ncsh(3) + 1
         ncsh(1) = ncsh(1) - 1
         ksym = 1
         nlast = 78
      end if
c
c ... now consider pseudopotentials ...
c ... decrease ncsh according to nelcor (nzps)
c     
      if (nzps.gt.0) then
c     pick up correct library of ECPs
c       write(6,*) 'IECPT = ', iecpt
c     CEP / SBKJC
        if (iecpt.eq.1.or.iecpt.eq.2) then
         ihw = 2
         if(nz.eq.55.or.nz.eq.56) ihw = 1
c     LANL / HW (version 6.2)
        else if (iecpt.eq.3) then
         ihw = 1
c     LANL2DZ
        else if (iecpt.eq.4) then
         ihw = 1
         if (nz.ge.72.and.nz.le.80) ihw = 2
c     CRENBL
        else if (iecpt.eq.5) then
         ihw = 1
         if (nz.ge.72) ihw = 2
c     CRENBS
        else if (iecpt.eq.6) then
         ihw = 1
c     Stuttgart RLC
        else if (iecpt.eq.7) then
         ihw = 1
c     Stuttgart RSC
        else if (iecpt.eq.8) then
         ihw = 1
         if (nz.ge.72.and.nz.le.80) ihw = 2
         if (nz.ge.89.and.nz.le.103) ihw = 2
        else
c
c        ECP specified by CARDS directive. Have to guess which 
c        configuration to use.
c
         ihw = 3
 110     nelec = 0
         do i = 1,17
            nelec = nelec + 4*isy(i,ihw)-2
            if (nelec.eq.nzps) then
c              this filling scheme will do
               goto 100
            else if (nelec.gt.nzps) then
c              this filling scheme is toast try something else
               if (ihw.eq.3) ihw = 1
               if (ihw.eq.1) ihw = 2
               if (ihw.eq.2) then
c                 good luck !!!
                  goto 100
               endif
               goto 110
            endif
         enddo
 100     continue
        endif
         nelec = 0
         do i = 1,17
            nelec = nelec + 4*isy(i,ihw) - 2
            if (nelec.le.nzps) then
               isymm = isy(i,ihw)
               ncsh(isymm) = ncsh(isymm) - 1
            else
               go to 90
            endif
         enddo       
      endif
c
c     check for cases where pseudo-potential changes  state
c     (cu  pseudo incl d10 =>s1)
c
90    if (nz.eq.29.and.nzps.eq.28) then
         nlast = 28
         ksym = 1
      end if
c
      do 70 i = 1 , 24
         ajmn(i) = 0.0d0
 70   continue
      if (nz.ne.nlast) then
         nosh(ksym) = 1
         nsym = max(nsym,ksym)
         nccup(ksym) = nz - nlast
      end if
c
c...  modify (tweak) configuration according to hand given 
c...  configurations
c
      call atcont('conf',oprin)
c.......................................................................
c
c     set k(l,l,0). moved to stajmn
c.......................................................................
c        ind = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
c        t = nccup(ksym)
c        f = 4*ksym - 4
c        ajmn(ind) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
c..      if (og) also set K(l,l,2) and K(l,l,4) (and more) if needed
c
cmarcin      print *,'nz: ',nz,' nlast: ',nlast
cmarcin      print *,'ksym: ',ksym
cmarcin      print *,'nsym: ',nsym,'  nccup: ',nccup(ksym),' nosh: ',nosh(ksym)
cmarcin      print *,'istate before stajmn: ',istate
cmarcin      print *,'iatdiff before stajmn: ',iatdiff
      call stajmn(nosh,nccup,ajmn,nsym,og,istate,iatdiff)
c
      return
      end
      subroutine atcont(choice,oprin)
c
c...  tweak the atomic configuration according to input specification
c...  e.g s1d4p3  or  s0s1s9 (one shell down) or s2s1 (one shell up)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      character*(*) choice
c
INCLUDE(common/sizes)
INCLUDE(common/datgue)
INCLUDE(common/runlab)
INCLUDE(common/iofile)
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
c
      dimension ncsh_save(5),nosh_save(5),nccup_save(5)
      save ncsh_save,nosh_save,nccup_save
c
      character*65 string
      character*5 yspdf,text
      character*10 prepre
      character*7 config(2)
      dimension nspdf(4)
      data yspdf/'spdf '/
      data prepre/'+-12345678'/
      data nspdf/2,6,10,14/
      data config/'atomscf','density'/
c
      if (choice.eq.'conf') then
         ichoice = 1
c...   save config
         do i=1,5
            ncsh_save(i) = ncsh(i)
            nosh_save(i) = nosh(i)
            nccup_save(i) = nccup(i)
         end do
      else if (choice.eq.'dens') then
         ichoice = 2
         do i=1,5
            ncsh(i) = ncsh_save(i)
            nosh(i) = nosh_save(i)
            nccup(i) = nccup_save(i)
         end do
      else
         call caserr('wrong choice in atcont')
      end if
c
      iatcon = 0
      do iat=1,natconf
         if (zaname(iatom).eq.zatconf(iat)) go to 10
      end do
      if (oprin) then
c...     just print the configuration
         is = 5
         go to 100
      end if
c
      return
c
10    string = string_ch(ichoice,iat)
      iatcon = iat
      if (atmode.eq.'zora'.and..not.forcat(iat)) iatcon = -iat
      if (string.eq.' ') return
c
c...  as only the last shell or the one after may be  specified 
c...  the number may be forgotten
c
      ii = 1
15    ipre = index(prepre,string(ii:ii))
      if (ipre.gt.0) ii = ii + 1
      is = index(yspdf,string(ii:ii))
      iforce = 0
      if (ipre.gt.0) then
         ii = ii + 1
         if (ipre.eq.1) iforce = +1
         if (ipre.eq.2) iforce = -1
         if (ipre-2.lt.ncsh(is)) call caserr(' specify only outside')
         if (ipre-2.gt.ncsh(is)) iforce = +1
         if (ipre-2.eq.ncsh(is)) iforce = -1
      end if
c
      if (is.eq.0) call caserr(' error in atomic conf')  
100   if (is.eq.5) then
c...     ... end ...
c...     this print will do unusual configurations
         ll = 0
         if (ncsh(1).ge.1) call addstring(string,ll,'1s2')
         if (ncsh(1).ge.2) call addstring(string,ll,'.2s2')
         if (ncsh(2).ge.1) call addstring(string,ll,'.2p6')
         if (ncsh(1).ge.3) call addstring(string,ll,'.3s2')
         if (ncsh(2).ge.2) call addstring(string,ll,'.3p6')
         if (ncsh(1).ge.4) call addstring(string,ll,'.4s2')
         if (ncsh(3).ge.1) call addstring(string,ll,'.3d10')
         if (ncsh(2).ge.3) call addstring(string,ll,'.4p6')
         if (ncsh(1).ge.5) call addstring(string,ll,'.5s2')
         if (ncsh(3).ge.2) call addstring(string,ll,'.4d10')
         if (ncsh(2).ge.4) call addstring(string,ll,'.5p6')
         if (ncsh(1).ge.6) call addstring(string,ll,'.6s2')
         if (ncsh(4).ge.1) call addstring(string,ll,'.4f14')
         if (ncsh(3).ge.3) call addstring(string,ll,'.5d10')
         if (ncsh(2).ge.5) call addstring(string,ll,'.6p6')
         if (ncsh(1).ge.7) call addstring(string,ll,'.7s2')
         if (ncsh(4).ge.2) call addstring(string,ll,'.5f14')
         if (ncsh(3).ge.4) call addstring(string,ll,'.6d10')
         if (ncsh(2).ge.6) call addstring(string,ll,'.7p6')
         if (ncsh(1).ge.7) call addstring(string,ll,'.8s2')
         do i=1,4
            if (nosh(i).eq.1) then
               if (nccup(i).lt.10) then
                  write(text,'(a1,i1,a1,i1)') 
     +                 '.',ncsh(i)+i,yspdf(i:i),nccup(i)
                  call addstring(string,ll,text(1:4))
               else
                  write(text,'(a1,i1,a1,i2)') 
     +                 '.',ncsh(i)+i,yspdf(i:i),nccup(i)
                  call addstring(string,ll,text(1:5))
               end if
            end if
         end do
         write(iwr,600) zaname(iatom),config(ichoice),string
600      format(/5x,' atom ',a8,' - configuration for ',a7,' : ',a65)
         return
      end if
      ii = ii + 1
      nn = 1
      if (index(yspdf,string(ii+1:ii+1)).eq.0) nn = 2
      if (nn.eq.1) then
         read(string(ii:ii),'(i1)') n
         ii = ii + 1
      else
         read(string(ii:ii+1),'(i2)') n
         ii = ii + 2
      end if
      if (n.eq.nspdf(is)) then
         ncsh(is) = ncsh(is) + 1
         if (iforce.ne.+1) then
            nosh(is) = 0
            nccup(is) = 0
         end if
      else if (n.gt.nspdf(is)) then
         call caserr(' shell overflow in atomic conf ')
      else if (n.eq.0) then
         if (nosh(is).eq.1) then
            nosh(is) = 0
            nccup(is) = 0
         else
            ncsh(is) = ncsh(is) -1
         end if
      else
         if (nosh(is).eq.0.and.iforce.lt.+1) ncsh(is) = ncsh(is) - 1
         nosh(is) = 1
         nccup(is) = n
      end if
c
      go to 15
c
      end
      subroutine stajmn(nosh,nccup,ajmn,nsym,og,
     +                  istate,iatdiff)
c
c...  set additional coupling coefficents for atomscf (Faegri/Almlof)
c...  the program already sets the (l,l,0) coupling coefficient
c...  at position ind, (now incorporated here)
c...  only for groundstate / involves a bit of guessing
c...  we assume the presence of just one open shell ....
c...      now dxs1 also allowed
c...  the coefficients are fractions => nominator/denominator pairs
c...  if (og) try to set all couplings for real groundstate
c...  modified by Marcin, 2008: now, excited states for p-shell only
c...  interface for different X-shell excited states has been prepared
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/datgue)
c
      dimension nosh(nsym),nccup(nsym),ajmn(*)
      logical og
      REAL kpp2,kdd2,kdd4,ksd2,kff2,kff4,kff6
      integer kpp2t(3,5),kpp2n(3,5)
      integer kdd2t(9),kdd2n(9),kdd4t(9),kdd4n(9),ksd2t(9),ksd2n(9)
      integer kff2t(13),kff2n(13),kff4t(13),kff4n(13)
      integer kff6t(13),kff6n(13),iindx
      character*2 zpshell(2,5)
c
c     p-shell
c                 ground states  2P  3P  4S  3P  2P
      data (kpp2t(1,i), i=1,5) /  2, -1, -2, -1,  2/
      data (kpp2n(1,i), i=1,5) / 15, 15, 15, 60,375/
c             excited states #1  XX  1D  2D  1D  XX
      data (kpp2t(2,i), i=1,5) /  0, 13,  2, 13,  0/
      data (kpp2n(2,i), i=1,5) /  0, 75, 75,300,  0/
      data (zpshell(1,i), i=1,5) /'XX','1d','2d','1d','XX'/
c             excited states #2  XX  1S  2P  1S  XX
      data (kpp2t(3,i), i=1,5) /  0,  8,  2,  2,  0/
      data (kpp2n(3,i), i=1,5) /  0, 15, 15, 15,  0/
      data (zpshell(2,i), i=1,5) /'XX','1s','2p','1s','XX'/
c
c     d-shell
c                   2D    3F    4F    5D    6S    5D    4F    3F    2D
      data kdd2t/    2,  -26,  -58,   -1,   -2,   -1, -174,  -13,    2/
      data kdd2n/   35,  245,  735,   20,   35,   45,12005, 1960, 2835/
      data kdd4t/    2,    9,  -34,   -1,   -2,   -1,  -34,    9,    2/
      data kdd4n/   35,  245, 2205,   20,   35,   45,12005, 3920, 2835/
c
c     s-shell - d shell 
c                   3D    4F    5F    6D    7S    6D    5F    4F    3D
      data ksd2t/   -1,   -1,   -1,   -1,   -1,   -2,   -3,   -1,   -1/
      data ksd2n/    5,    5,    5,    5,    5,   15,   35,   20,   45/
c
c     f-shell
c     Thanks to Mariusz Klobukowski (Alberta, Canada, March 2001)
c
c            n       1        2       3        4        5       6     7 
c                   2F       3H      4I        5I       6H     7F    8S
c                    8        9      10        11       12      13
c                   7F       6H      5I        4I       3H      2F
      data kff2t/    4,     -23,   -256,     -17,    -344,    -34,   -4,
     +             -17,    -344,    -17,    -256,     -23,      4/
      data kff2n/  105,     315,   2835,     252,    7875,    945,  105,
     +             840,   25515,   1575,   38115,   11340,  17745/
      data kff4t/    2,     -53,   -722,     -74,   -1598,    -17,   -2,
     +             -17,   -1598,   -296,    -722,     -53,      2/
      data kff4n/   77,    2541,  22869,    2541,   63525,    693,   77,
     +            1232,  205821,  63525,  307461,   91476,  13013/
      data kff6t/  100,    3125,  17800,   -1325,   -2272,   -850, -100,
     +            -425,  -56800,    -53,   17800,    3125,    100/
      data kff6n/ 3003,   99099, 891891,  396396,   99099,  27027, 3003,
     +           24024, 8027019,  99099,11990979, 3567564, 507507/
c
      do ksym=1,nsym
        if (nosh(ksym).eq.1) then
c           set k(l,l,0)
          ind = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
          t = nccup(ksym)
          f = 4*ksym - 4
          ajmn(ind) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
c
          if (odiff) then
            if (og) then
              iindx = 1
            end if
            nel = nccup(nsym)
c
            if (ksym.eq.1) then
c
c...   s-shell
c
cmarcin              print *,'1 ksym: ',ksym,'  ind: ',ind
            else if (ksym.eq.2) then
c
c...   p-shell
c
              if (oexc) then
cmarcin                print *,'zatstate: ',zatstate(iatdiff)
cmarcin                print *,'zpshell 1: ',zpshell(1,nel)
cmarcin                print *,'zpshell 2: ',zpshell(2,nel)
c                do ij=1,5
                  if ( zatstate(iatdiff) .eq. zpshell(1,nel) ) then
                    iindx = 2
                  else if ( zatstate(iatdiff) .eq. zpshell(2,nel) ) then
                    iindx = 3
                  end if
c                end do
              end if
              kpp2 = (1.0d0*kpp2t(iindx,nel))/kpp2n(iindx,nel)
cmarcin              print *,'2 ksym: ',ksym,'  ind: ',ind
              ajmn(ind+1) = kpp2/2.0d0
c
c...   d-shell
c
            else if (ksym.eq.3) then
              kdd2 = (1.0d0*kdd2t(nel))/kdd2n(nel)
              kdd4 = (1.0d0*kdd4t(nel))/kdd4n(nel)
              ajmn(ind+1) = kdd2/2.0d0
              ajmn(ind+2) = kdd4/2.0d0
cmarcin              print *,'3 ksym: ',ksym,'  ind: ',ind
c....  d-shell - s-shell
              if (nosh(1).eq.1) then
cmarcin                 print *,'s-d- coupling'
                 ksd2 = (1.0d0*ksd2t(nel))/ksd2n(nel)
                 ajmn(5) = ksd2/2.0d0
c.....              (index found by hvd)
              end if
c
c...   f-shell
c
            else if (ksym.eq.4) then
              kff2 = (1.0d0*kff2t(nel))/kff2n(nel)
              kff4 = (1.0d0*kff4t(nel))/kff4n(nel)
              kff6 = (1.0d0*kff6t(nel))/kff6n(nel)
              ajmn(ind+1) = kff2/2.0d0
              ajmn(ind+2) = kff4/2.0d0
              ajmn(ind+3) = kff6/2.0d0
              if (nosh(1).eq.1) then
                 print *,' sorry, no sf-shells implemented '
                 print *,' reverting to only K(f,f,x) usage'
              end if
            end if
          end if
        end if
      end do
c
      return
      end
**==atomd.f
      subroutine atomd(og,oprin,iwr,znps,iecp,ic,isymax,hatom,
     + pcap,qcap,fc,fo,s,u,t,h,dc,dos,dt,dold,ss,
     + cvec, copn, smin, qmin, transf, cc , dalpha, dbeta, nbb, 
     + istate, iatdiff)
c
      implicit REAL  (a-h,o-z),integer   (i-n)
      parameter (nbig=1500, no=50)
      dimension ic(6,*),hatom(*)
      dimension pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      dimension h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      dimension cvec(*),copn(*),smin(nbb,*),qmin(nbb,*),transf(*),cc(*)
      dimension dalpha(*),dbeta(*)
c
      logical oprin,og,oprin_save,opg_root
c.......................................................................
c     atomic r h f code for gto-s. uses roothaan double diagonalization.
c	Roothaan, C.C.J. and P.S. Bagus, ATOMSCF. 
c                 Methods Comput. Phys. 2, 47, 1963. 2: p. 47
c       Rev.by B. Roos, C. Salez, A. Veillard, and E. Clementi, 
c                 IBM Research Report RJ518 (1968).
c       Rev. by K. Faegri for use in Disco
c       Rev. by J.H. van Lenthe, R. Zwaans and  H.J.J. van Dam 
c                 for use in GAMESS-UK and NWCHEM.
c.......................................................................
c
c      routine datoms .. tramad are a totally separate unit
c      they use commom junk to communicate
c      oprin,iwr are new parameters to control printing
c      hatom are the pseudo-corrections in the contracted basis from gam
c.......................................................................
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/datgue)
c
c     nsht      = total number of shells
c     n1(i)     = nbas(i) * (nbas(i) + 1 ) / 2
c     nbc(i)    = number of cont. orbitals in symmetry i
c     cont(i)   = contraction coeff. assosiated with primitive no. i
c     nstrt(i)  = number for first primitive in cont. no. i
c     nbct      = total number of cont. basis functions.
c.......................................................................
c     zn = znx
      if (dabs(zn).ge.1.d-8) then
c.......................................................................
c
c     distribute electrons according to aufbau
c.......................................................................
cmarcin         print *,'zn before atcond: ',zn
cmarcin         print *,'istate before atcond: ',istate
cmarcin         print *,'iatdiff before atconf: ',iatdiff
         call atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps,iecp,og,oprin,
     +               istate,iatdiff)
         oprin_save = oprin
         if (iatcon.gt.0) oprin = .true.
crz
c     fix for use of pseudopotentials
         znsave = zn
         zn = znps
c.......................................................................
c
c     move basis set information from transfer variables to working
c     variables. ** not necessary in present version
c.......................................................................
         nbct = 0
         ndim = 0
         nsqt = 0
         nbcdim = 0
         jcount = 0
         nsht = 0
         k = 0
         do 40 i = 1 , nsym
            n1(i) = nbas(i)*(nbas(i)+1)/2
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
            ndim = ndim + n1(i)
            nsht = nsht + ncsh(i) + nosh(i)
            nsqt = nsqt + nbas(i)**2
            nbct = nbct + nbc(i)
            do 20 j = 1 , nbas(i)
               jcount = jcount + 1
               nqn(jcount) = i
 20         continue
            do 30 j = 1 , n1(i)
               k = k + 1
               dold(k) = 0.0d0
 30         continue
 40      continue
         ns = 1
         nstrt(ns) = 1
         do 60 l = 1 , nsym
            iant = nbc(l)
            ksum = 0
            do 50 i = 1 , iant
               nstrt(ns+i) = nstrt(ns+i-1) + ic(l,i)
               ksum = ksum + ic(l,i)
 50         continue
            if (ksum.ne.nbas(l)) write (iwr,6020) l , ksum , nbas(l)
            ns = ns + iant
 60      continue
         do 70 i = 1 , nsqt
            cvec(i) = 0.d0
            copn(i) = 0.0d0
            cc(i) = 0.d0
 70      continue
         do 80 i = 1 , ndim
            dold(i) = 0.0d0
 80      continue
cjvl  few extra checks
         do 81 i=1,nsym
            if (nbc(i).lt.ncsh(i)+nosh(i)) then
               oprin = .true.
               ndim = 0
            end if
81       continue
cjvl      
         if (oprin) write (iwr,6030) zn
         if (oprin) write (iwr,6040) (nbas(i),i=1,nsym)
         if (oprin) write (iwr,6050) (nbc(i),i=1,nsym)
         if (oprin) write (iwr,6060) (ncsh(i),i=1,nsym)
         if (oprin) write (iwr,6070) (nosh(i),i=1,nsym)
         if (oprin) write (iwr,6080) (nccup(i),i=1,nsym)
         if (oprin) write(iwr,1700)(2*ajmn(i),i=1,20),(ajmn(i),i=21,24)
cjvl
         if (ndim.eq.0) then
            if (opg_root()) then
               write(iwr,6100)
            endif
            call caserr2('error in atom scf (see above)')
         endif
cjvl
         maxitr = 100
c..
c..     calculate 1-electron ints
c..
         call oeigd(fc,s,u,t,h)
c.......................................................................
c
c     copy overlap matrix to ss
c.......................................................................
         do 90 i = 1 , ndim
            ss(i) = s(i)
            fc(i) = s(i)
 90      continue
c.......................................................................
c
c     now transform ss to contracted basis, then set up transformation
c     matrix to o.n. contracted basis.
c.......................................................................
         call trafsd(nsym,nbas,ndim,ss,nbc,cont,nstrt,dc)
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)
c...
         nstep1 = 1
         nstep2 = 1
         do 100 i = 1 , nsym
crz      to surpress problems when there is one type of shell missing...
            if (nbc(i).ne.0) then
               call shalfd(fc(nstep1),transf(nstep2),nbc(i),nbb)
               call starcd(cc(nstep2),ss(nstep1),nbc(i),ncsh(i),
     +                     nosh(i))
               nstep1 = nstep1 + nbc(i)*(nbc(i)+1)/2
               nstep2 = nstep2 + nbc(i)**2
            end if
 100     continue
         nitscf = 0
         nconv = 0
         damp = .30d0
c        for Yb the damping is more and this converges in 61 cycles
         if (zn.eq.70.0) damp = .70d0
 110     nitscf = nitscf + 1
c.......................................................................
c
c     transform vectors and set up matrices in primitive basis,
c     then transform fock matrices to contracted basis.
c.......................................................................
         call tracd(cvec,cc,nsqt)
         call densid(dt,dold,dos,dalpha,dbeta,nsym,nosh,ncsh,
     +               nccup,cvec,damp,nconv,nbas,nitscf,tlarge,eps)
c
c... check for convergence on tlarge (max change of d-matrix)
c
         if (tlarge.le.1.0d-5) nconv = 1
         if (nitscf+20.ge.maxitr) then
            write (iwr,6010) nitscf , energ , tlarge
         end if
c
         call hamild(pcap,qcap,fc,fo,s,u,t,h,dos,dt,cvec,smin,qmin,nbb)
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)
         call trafsd(nsym,nbas,ndim,fo,nbc,cont,nstrt,dc)
c...
c...    now add the h-pseudo-contributions 
c...
         do 120 i = 1 , nbcdim
            fc(i) = fc(i) + hatom(i)
            fo(i) = fo(i) + hatom(i)
 120     continue
c.......................................................................
c
c     do double diagonalization by symmetries:
c         1. transform block to o.n.basis (contracted).
c         2. store o.n. transformation matrix in vector matrix.
c         3. diagonalize.
c         4. order eigenvectors by eigenvalue.
c         5. if necessary, merge open and closed vectors.
c.......................................................................
         nstep1 = 1
         nstep2 = 1
         knteps = 0
         do 160 i = 1 , nsym
            nbc1 = nbc(i)
            nbc2 = nbc1**2
            nbc3 = (nbc2+nbc1)/2
            if (ncsh(i).ne.0) then
               call tramad(fc(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,cc(nstep2),1)
               call jacod(fc(nstep1),cc(nstep2),nbc1,n1(i),nbc2,1,nbc1,
     +                    dc,dt,nbc1)
            end if
            if (nosh(i).ne.0) then
               call tramad(fo(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,copn(nstep2),1)
               call jacod(fo(nstep1),copn(nstep2),nbc1,n1(i),nbc2,1,
     +                    nbc1,dc,dt,nbc1)
            end if
            icount = nstep1
            do 130 j = 1 , nbc1
               dc(j) = fc(icount)
               dos(j) = fo(icount)
               icount = icount + 1 + j
 130        continue
            call orderd(cc(nstep2),nbc1,nbc1,idum,idum,dt,dc,100)
            if (nosh(i).gt.0) then
               call orderd(copn(nstep2),nbc1,nbc1,idum,idum,dt,dos,100)
               call cmergd(cc(nstep2),copn(nstep2),ncsh(i),nbc1,nosh(i))
            end if
            nstep1 = nstep1 + nbc1*(nbc1+1)/2
            nstep2 = nstep2 + nbc2
            if (nconv.gt.0) then
               if (ncsh(i).gt.0) then
                  do 140 j = 1 , ncsh(i)
                     knteps = knteps + 1
                     eps(knteps) = dc(j)
 140              continue
               end if
               if (nosh(i).gt.0) then
                  do 150 j = 1 , nosh(i)
                     knteps = knteps + 1
                     eps(knteps) = dos(ncsh(i)+j)
 150              continue
               end if
            end if
 160     continue
         if (nitscf.ge.maxitr) nconv = 1
         if (nconv.le.0) go to 110
c
c...  set configuration to determine density
c
         call atcont('dens',.false.)
c
         call densid(dt,dold,dos,dalpha,dbeta,nsym,nosh,ncsh,
     +               nccup,cc,damp,nconv,nbc,nitscf,tlarge,eps)
c
c...     correct energy for pseudopential or zora
c
         energ = energ + ddot(nbcdim,hatom,1,dt,1)
c        
         lm = 0
         do 190 i = 1 , nsym
            noff = lm
            nbci = nbc(i)
            do 180 l = 1 , nbci
               ll = noff + l*(l+1)/2
               do 170 m = 1 , l
                  mm = noff + m*(m+1)/2
                  lm = lm + 1
                  dt(lm) = dt(lm)*dsqrt(ss(ll)*ss(mm))
                  if (m.ne.l) dt(lm) = dt(lm)/2.0d0
                  if (uhfatom) then
c...              now dt - alpha / dos beta
                     dos(lm) = dos(lm)*dsqrt(ss(ll)*ss(mm))
                     if (m.ne.l) dos(lm) = dos(lm)/2.0d0
                  end if
 170           continue
 180        continue
 190     continue
         if (oprin) then
            call outpud(copn,cc,1,iwr)
         else
            call outpud(copn,cc,0,iwr)
         end if
      else
c.......................................................................
c
c     special section for handling the case of floating functions
c     on centers with no charge.
c.......................................................................
         nbcdim = 0
         do 200 i = 1 , isymax
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
 200     continue
         do 210 i = 1 , nbcdim
            dt(i) = 0.0d0
 210     continue
         nsym = isymax
         energ = 0.0d0
      end if
c
      oprin = oprin_save
c.......................................................................
c
c     output atomic density matrix to tape.
c.......................................................................
c        do 230 j=1,nurk
c        write(ntdens) nsym,isymax,(nxxbc(i),i=1,isymax),nbcdim,energ,
c    x        dt(i),i=1,nbcdim)
c 230    continue
c...     copy atomic density matrix to output array
c
      return
 6010 format (' it.',i4,'  energy',d19.10,'  div.',d13.5)
 6020 format ('-',' wrong contraction in symmetry   ',3i5)
 6030 format (/'      charge =',f10.6,//'      symmetry species',12x,
     +        's',5x,'p',5x,'d',5x,'f')
 6040 format (6x,'number of basis functions =',4(i2,4x))
 6050 format (6x,'number of cont. functions =',4(i2,4x))
 6060 format (6x,'number of closed shells   =',4(i2,4x))
 6070 format (6x,'number of open shells     =',4(i2,4x))
 6080 format (6x,'open shell occupation     =',4(i2,4x))
 6100 format (//,
     +     ' *** The error below is most likely caused by the number',
     +   /,' *** of electrons exceeding the capacity of the basis set.',
     +   /,' *** This is often due to:',
     +   /,' *** 1) specifying an ecp basis set but omitting to ',
     +     'specify the ECP itself',
     +   /,' *** 2) associating an unsuitable basis set with an atom',
     +   /,' *** Please make sure your input is correct.',//)
 1700 format(/'      vector coupling coefficients',/(6f16.8))
      end
**==cmergd.f
      subroutine cmergd(c1,c2,nc,nb,no)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     merge open and closed shell coefficient matrices.
c.......................................................................
      dimension c1(nb,nb), c2(nb,nb)
      nlast = nc + no
      nrow = nc + 1
      do 30 i = nrow , nlast
         do 20 j = 1 , nb
            c1(j,i) = c2(j,i)
 20      continue
 30   continue
      return
      end
**==complt.f
      subroutine complt(ar,b,aa,mdim,nbasis,nvec,north)
c
c    *******************************************************************
c
c          this subroutine normalizes, orthogonalizes, and completes
c     the mo coefficient matrix.  the real mo coefs are in ar
c     and the vectors used to complete the
c     coefs are in b.  aa, bb, and cc are scratch vectors.  mdim is
c     the dimension of the matrices, and nvec is the number of
c     vectors provided.
c
c    *******************************************************************
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ar(mdim,mdim),b(mdim,mdim)
      dimension aa(mdim)
      data big/1.0d10/
c     normalize and orthogonalize the vectors provided.
      if (nvec.ne.north) then
         if (north.le.1) call gnorm(ar(1,1),nbasis)
         lim = max(north,2)
c
         do 30 i = lim , nvec
            call vclr(aa,1,nbasis)
            im1 = i - 1
            do 20 j = 1 , im1
               dotr = ddot(nbasis,ar(1,i),1,ar(1,i),1)
               call daxpy(nbasis,dotr,ar(1,j),1,aa(1),1)
 20         continue
            call vsub(ar(1,i),1,aa,1,ar(1,i),1,nbasis)
      dotr=1.0d0/dnrm2(nbasis,ar(1,i),1)
      call dscal(nbasis,dotr,ar(1,i),1)
 30      continue
      end if
c
c     complete the remaining vectors.
c
      nstart = nvec + 1
      if (nstart.gt.nbasis) return
      do 70 ivec = nstart , nbasis
         dotmin = big
         jvec = ivec - 1
c
         do 50 jb = 1 , nbasis
            sumdot = 0.0d0
c
            do 40 j = 1 , jvec
               dotr = ddot(nbasis,ar(1,j),1,b(1,jb),1)
               sumdot = sumdot + dotr*dotr
 40         continue
            if (sumdot.le.dotmin) then
               dotmin = sumdot
               jbmin = jb
            end if
 50      continue
c
         call dcopy(nbasis,b(1,jbmin),1,ar(1,ivec),1)
c
         call vclr(aa,1,nbasis)
         do 60 j = 1 , jvec
            dotr = ddot(nbasis,ar(1,j),1,ar(1,ivec),1)
            call daxpy(nbasis,dotr,ar(1,j),1,aa,1)
 60      continue
c
         call vsub(ar(1,ivec),1,aa,1,ar(1,ivec),1,nbasis)
      dotr = 1.0d0/dnrm2(nbasis,ar(1,ivec),1)
      call dscal(nbasis,dotr,ar(1,ivec),1)
 70   continue
      return
      end
**==creded.f
      subroutine creded(d,db,dt,dbeta,iiloc,nbb,uhfatom)
c
      implicit REAL  (a-h,o-z),integer   (i-n)
c
c..   routine to merge atomic density matrix into molecular one
c..   is called straight after atom, so all info is still in
c..   common /junk/ : i.e. dt,nbc
c
      dimension d(*),db(*),dt(*),dbeta(*),iiloc(nbb,6)
      logical uhfatom
c
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
c
      dimension dmult(145)
c..
crz   order of d-s and f-s in GAMESS is completely different
crz   from MOLECULE !!
c
c     dmult consists of :
c     nothing for s;
c      p functions (1 line)
c      d functions (next 4 lines)
c      f functions (last lines)
c      -0.670820393 = -0.3*sqrt(5) (used for f)
c
      data dmult/ 1.0d0, 3*0.0d0, 1.0d0, 3*0.0d0, 1.0d0,
     x  1.0d0, -0.5d0, -0.5d0, 3*0.0d0,
     x -0.5d0,  1.0d0, -0.5d0, 3*0.0d0,
     x -0.5d0, -0.5d0,  1.0d0, 3*0.0d0,
     x  3*0.0d0, 1.0d0, 6*0.0d0, 1.0d0, 6*0.0d0, 1.0d0,
     x  1.0d0,   4*0.0d0, -0.670820393d0, 0.0d0,-0.670820393d0, 2*0.0d0,
     x  0.0d0,1.0d0,0.0d0,-0.670820393d0, 4*0.0d0,-0.670820393d0, 0.0d0,
     x2*0.0d0,1.0d0,0.0d0,-0.670820393d0, 0.0d0,-0.670820393d0, 3*0.0d0,
     x   0.0d0, -0.670820393d0,  0.0d0,  1.2d0, 4*0.0d0, -0.3d0,  0.0d0,
     x 2*0.0d0, -0.670820393d0,  0.0d0,  1.2d0,  0.0d0, -0.3d0, 3*0.0d0,
     x  -0.670820393d0,         4*0.0d0, 1.2d0,  0.0d0, -0.3d0, 2*0.0d0,
     x 2*0.0d0, -0.670820393d0,  0.0d0, -0.3d0,  0.0d0,  1.2d0, 3*0.0d0,
     x  -0.670820393d0,     4*0.0d0,    -0.3d0,  0.0d0,  1.2d0, 2*0.0d0,
     x   0.0d0, -0.670820393d0, 0.0d0, -0.3d0, 4*0.0d0,  1.2d0,   0.0d0,
     x 9*0.0d0,  1.0d0/
c.......................................................................
c
c         among the variables used are:
c
c              nsym       - highest l-quantum no. used in atomic calc.
c              dt         - atomic density matrix (for uhf alpha)
c              dbeta      - atomic density matrix (for uhf beta)
c              d          - area for final molecular density matrix. (a for uhf)
c              db         - area for final beta molecular density matrix.
c
c        jvl  daresbury 1988
c        added UHF 2 density matrices jvl 2004
c.......................................................................
c..
c..    triangle statement function
c..
      itrian(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c..
      lm = 0
      do 40 k = 1 , nsym
         nbci = nbc(k)
         if (k.eq.1) then
c.......................................................................
c
c       s-orbitals, simple distribution of matrix elements.
c
c.......................................................................
            do 30 l = 1 , nbci
               do 20 m = 1 , l
                  noff = itrian(iiloc(l,1),iiloc(m,1))
                  lm = lm + 1
                  d(noff) = dt(lm)
                  if (uhfatom) db(noff) = dbeta(lm)
 20            continue
 30         continue
         else if (k.le.4) then
            kdim = k*(k+1)/2
            kmone = k - 1
            kpoint = kmone*(kmone+1)*(((6*kmone+24)*kmone+26)*kmone+4)
     +               /120
            factor = 1.d0/(k+k-1)
            lm_save = lm
            call pdfded(kdim,d,dt,factor,dmult(kpoint),lm,nbci,
     +                  iiloc(1,k))
            if (uhfatom) then 
               lm = lm_save
               call pdfded(kdim,db,dbeta,factor,dmult(kpoint),lm,nbci,
     +                     iiloc(1,k))
            end if
         else
            call caserr('unexpected shell in creded')
         end if
 40   continue
c..
      return
      end
**==creorb.f
      subroutine creorb(c,ndim,eorb,oorb,cc,iiloc,nbb,iorb)
c
      implicit REAL  (a-h,o-z),integer   (i-n)
c
c..   routine to merge atomic orbital matrix into molecular one
c..   is called straight after atom, so all info is still in
c..   common /junk/ : i.e. dt,nbc
c..   we try to get them in orbital energy order(eps)
c..   but first closed shells then open shells
c     This routine works only for harmonic / adapt off
c     If you do want the orbitals in another representation, you might
c     consider a vectors getq in a following job
c
      dimension c(ndim,ndim),eorb(ndim),oorb(ndim),cc(*),iiloc(nbb,6)
c
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
c
INCLUDE(common/sizes)
INCLUDE(common/tran)
      dimension idim(4),ndone(4)
      data idim/1,3,5,7/
c..
      ncsht = 0
      do i=1,nsym
        ndone(i) = 0
        ncsht = ncsht + ncsh(i)
      end do
c..
      do ish=1,nsht
c...   look for lowest energy
c...   first do all closed shell's then the open shells
         e = 0.0d0
         kk = 0
         do i=1,nsym
            kk = kk + ndone(i)
            jend = ncsh(i)
            if (ish.gt.ncsht) jend = jend + nosh(i)
            do j=ndone(i)+1,jend
               kk = kk + 1
               if (eps(kk).lt.e) then
                  e = eps(kk)
                  ilow = i
                  jlow = j
               end if
            end do
         end do
c
         lm = 0
         do i=1,ilow - 1
            lm = lm + nbc(i)*nbc(i)
         end do
         lm = lm + ndone(ilow)*nbc(ilow)
c
         do m = 1 , nbc(ilow)
            lm = lm + 1
            iil = iiloc(m,ilow)
            do i=1,idim(ilow)
c....           skip those not in harmonic
10             if (ctran(ilifc(iil+i-1)+1).eq.0.0d0) then
                  iil = iil + 1
                  go to 10
               end if
c
               c(iil+i-1,iorb+i) = cc(lm)
c
            end do
         end do
c...  set orbital energies
         do i=1,idim(ilow)
            eorb(iorb+i) = e
         end do
c...  set occupations
         if (ndone(ilow)+1.le.ncsh(ilow)) then
            o = 2.0d0
         else
            if (nosh(ilow).ne.1) call caserr2('creorb nosh crash')
            o = nccup(ilow)
            o = o/(idim(ilow)*1.0d0)
         end if
         do i=1,idim(ilow)
            oorb(iorb+i) = o
         end do
c
         iorb = iorb + idim(ilow)
         ndone(ilow) = ndone(ilow) + 1
c..
      end do
c
      return
      end
**==datoms.f
      subroutine datoms(hbig,hatom,d,db,orb,eorb,oorb,oprin,
     + pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     + cvec, copn, smin, qmin, transf, cc, 
     + ic, iiloc, iisch, odone, dalpha,dbeta ,nbb,toteng,q)
c
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
c
c...   subroutine to coordinate atom-scf calls and d-matrix gathering
c...   for atomic startup
c...   hbig  : correction to  full h-matrix to supply to atom
c...   d     : full density matrix as return parameter
c...   alternatively the atomic orbitals may be combined to give
c...   then the orbitals and orbital occupations/energies 
c...   are written to disk and returned in orb/eorb/oorb
c...
c...   **note** data is transferred to atom directly via common/junk/
c
c...   for UHF (uhfatom) d and db are de alpha and betaq density matrices
c
INCLUDE(common/sizes)
      REAL orb,oorb
      parameter (nbig=1500,no=50)
      dimension hbig(*),hatom(*),d(*),orb(*),eorb(*),oorb(*)
      dimension pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      dimension h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      dimension cvec(*),copn(*),smin(nbb,*),qmin(nbb,*),transf(*),cc(*)
      dimension ic(6,nbb),iiloc(nbb,6),iisch(nbb,6),odone(*)
      dimension q(*), dalpha(*),dbeta(*) , db(*)
c
c...   commons where information comes from :
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/infob)
INCLUDE(common/datgue)
INCLUDE(common/runlab)
INCLUDE(common/harmon)
      common /saveco/ udum,idadap
INCLUDE(common/nshel)
      common/junk2/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +             ctran(mxorb3),
     +             ilifct(maxorb),ntrant(maxorb),itrant(mxorb3),
     +             ctrant(mxorb3)
      dimension cspd(mxprim,4)
      equivalence (cs(1),cspd(1,1))
c...   common where info goes to :
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
INCLUDE(common/phycon)
INCLUDE(common/restri)
c
c..    we need xy for d and xyz for f in gathering h-ints
c..    so set proper offsets in ioffhp (see do 140)
c
      dimension ioffhp(4)
      data ioffhp/0,0,3,9/
c..
      data pi32/5.56832799683170d0/
      data pt5,pt75,dzero/0.5d0,0.75d0,0.0d0/
      data pt187 /1.875d+00/
      data m100/100/
c..
c..    triangle statement function (not heavily used)
c..
      itrian(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c..
      if (nbb.gt.nbig) call caserr2('datoms nbig clash')
c
      nav = lenwrd()
      if (oatdo) then
         if(.not.oharm) call caserr2('atorbs requires harmonic')
         if (idadap.eq.0) 
     1      call caserr2('atorbs requires adapt off')
         call vclr(orb,1,num*num)
         call vclr(eorb,1,num)
         call vclr(oorb,1,num)
      endif
      call vclr(d,1,nx)
      if (uhfatom) call vclr(db,1,nx)
      toteng = dzero
c
      do 20 i = 1 , nat
         odone(i) = .false.
 20   continue
      norb = 0
      iatom = -1
c
c..   the full nuclear charges are in infob (ctranz) in the right order
c
c...  loop over atoms like adapt does
c
120   continue
cjvl      odone_one = .false.
      do 110 iat = 1 , nat
         if (odone(iat)) go to 110
c
c ... eliminate magic/ghost/point charges/dummies
c
         if (nuct(iat).ne.1) then
            odone(iat) = .true.
            go to 110
         end if
c..
c..   check if we have already treated this one or if it is of same
c..   type as the one we already did
c..   to complicated now to do only unique atoms; limit ourselves to 
c..   ton atoms close together, iteration disabled
c..
cjvl        if (odone_one) then
            if (osatod(iat,ic,iiloc,iisch,nbb)) then
               go to 100
            end if
cjvl           go to 110
cjvl        end if
cjvl        odone_one = .true.
c
c...  gather  shell / symmetry info
c
         iatom = iat
         do 30 i = 1 , 4
            nbc(i) = 0
 30      continue
c
c.. nbc  # shell-s / symmetry
c.. iisch  contains index of shell
c.. iiloc  contains position of starting ao of shell in "real" world
c
         do 50 ii = 1 , nshell
            i = katom(ii)
            if (i.eq.iat) then
               mini = kmin(ii)
               maxi = kmax(ii)
               kk = ktype(ii)
               do 40 iorb = mini , maxi
                  if (iorb.eq.1) then
c..  translate to 1 (s)
                     nbc(1) = nbc(1) + 1
                     iisch(nbc(1),1) = ii
                     iiloc(nbc(1),1) = kloc(ii)
                  else if (iorb.eq.2 .or. iorb.eq.5 .or. iorb.eq.11)
     +                     then
c..  translate to 2 (p) 3(d) or  4(f)
                     ispdf = kk
                     nbc(ispdf) = nbc(ispdf) + 1
                     iisch(nbc(ispdf),ispdf) = ii
                     iiloc(nbc(ispdf),ispdf) = kloc(ii) + iorb - mini
                  end if
 40            continue
            end if
 50      continue
c..
c..     we gathered symmetry/shell info ; now get the real thing
c..
         kkzc = 0
         kh = 0
         isymax = 0
         do 90 ispdf = 1 , 4
c..      nbas = total # primitives for this symmetry
            nbas(ispdf) = 0
            if (nbc(ispdf).gt.0) isymax = ispdf
            do 80 j = 1 , nbc(ispdf)
               ii = iisch(j,ispdf)
               is = kstart(ii)
               if = is + kng(ii) - 1
c..      ic = # number of primitives /contracted /symmetry
               ic(ispdf,j) = kng(ii)
               nbas(ispdf) = nbas(ispdf) + kng(ii)
               if (kkzc+kng(ii).gt.nbb) 
     1            call caserr2('nbb/nbig error in datoms')
c..      gather the primitives / watch the subtle use of 2-dim cspd
               do 60 k = is , if
                  kkzc = kkzc + 1
                  zeta(kkzc) = ex(k)
                  cont(kkzc) = cspd(k,ispdf)
c...     get contraction coeff-s as we are used to
                  ee = 2*zeta(kkzc)
                  fac = pi32/(ee*dsqrt(ee))
                  if (ispdf.eq.2) then
                     fac = pt5*fac/ee
                  else if (ispdf.eq.3) then
                     fac = pt75*fac/(ee*ee)
                  else if (ispdf.eq.4) then
                     fac = pt187*fac/(ee**3)
                  end if
                  cont(kkzc) = cont(kkzc)*dsqrt(fac)
 60            continue
c...     gather the correction (ecp,zora,dft) integrals 
c...     for the contracted ao-s
c...     they are added in at the right time (in atomd)
c...     use proper offset to use pure d or f functions (ioffhp)
               do 70 k = 1 , j
                  itria = itrian(iiloc(k,ispdf)+ioffhp(ispdf),
     *                           iiloc(j,ispdf)+ioffhp(ispdf))
                  kh = kh + 1
                  hatom(kh) = hbig(itria)
 70            continue
c...
 80         continue
 90      continue
c..
c..     all prepared call  atomd
c..     zeta,cont,nbas,nbc,nbas,ic,zn are passed via junk
c..     energ and the density matrix dt are received via junk
c..     note zn is the real nuclear charge / znps is the effective charg
c..
         zn = czanr(iat)
         istate = iatstates(iat)
         do ij=1,nat
           do ijj=1,natdiff
             if ( zaname(ij) .eq. zatdiff(ijj) ) then
               iatdiff = ijj
             end if
           end do
         end do
         znps = czan(iat)
         iecp = ipseud(iat)
c        write(6,*)' CALL ATOMD: iecp = ', iecp
c
         call atomd(oground,oprin,iwr,znps,iecp,ic,isymax,hatom,
     +       pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +       cvec, copn, smin, qmin, transf, cc, dalpha, dbeta, nbb, 
     +       istate,iatdiff)
c..
 100     continue
c..
c..      now add density-matrix to the molecular d-matrix
c..
         if (oatdo) then
c..       create combined orbitals (atorb)
            call creorb(orb,num,eorb,oorb,cc,iiloc,nbb,norb)
         endif
c..       normal mode; add density matices (atdens)
         call creded(d,db,dt,dos,iiloc,nbb,uhfatom)
c..
cmarcin         print *,'atomic energy of atom: ',energ
         toteng = toteng + energ
         odone(iat) = .true.
 110  continue
cjvl     if (odone_one) go to 120
c..
      if (oground) write(iwr,6011)
      if (oexc) then
        write(iwr,6012)
        do i=1,nat
          if ( iatstates(i) .eq. 1 ) then
            write(iwr,'(7x,a,a,a,a)') 
     +            'atom ',zatdiff(i),' in state ',zatstate(i)
          else
            write(iwr,'(7x,a,a,a)') 'atom ',zaname(i),' in ground state'
          end if
        end do
      end if
      if (oatdo) then
c
c..   complete vector set
c
        call compv(orb,norb,num,q,'add')
      end if
c..
      return
 6011 format (1x,'***** using real groundstate atoms *** ')
 6012 format (1x,'***** using real excitedstate atoms *** ')
      end
**==denat.f
      subroutine denat(q,mode)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c..   get starting vectors from superposition of atomic densities
c..   also used to calculate atomic zora corrections (jvl 2000)
c..   mode = 'start' : atomic startup
c..   mode = 'zora'  : zora corrections
c..   *note* making zora really atomic would make it more efficient
c..
c..   routines involved are (all present in this order) :
c..   denat  : general control routine
c..   datoms : translate orbital info, call atomscf , gather d-matrix
c..   atomd,tracd,trafsd,cmergd,oeigd,teigd,densid,denmad, ..
c..   .. hamild,outpud,shalfd,atcond,starcd,creded,pdfded, ..
c..   .. jacod,orderd,tramad
c..   oeigd has been changed to pick up the 1-elec ints
c..   atcond has been changed to allow for effective nuclear charge
c..   creded/pdfded have been adapted to yield directly a d-matrix
crz   pseudopotentials now picked up from ed7 (ibl7ec) (jvl,2001)
c..
c...    start all types after one cycle closed shell in denscf
c...    exception:UHF/UHFATOM; a- and b-density matrices are used in UHF
c
chvd
c     This code is called as follows:
c     In the case of a ZORA calculation (with or without DFT) denat
c     is called twice:
c     1) as DENAT(mode='zora')
c     2) as DENAT(mode='start')
c     In the case of no-ZORA and DFT atomic startup DENAT is called once
c     1) as DENAT(mode='zora') (but oatmdft is .true.)
c     In all other cases DENAT is called once 
c     1) as DENAT(mode='start) (but oatmdft is .false.)
c
c     Furthermore note the strange way in which the explicit 2-el terms
c     are modified. I.e. CD_set_2e is called within an IF statement 
c     whereas the CD_reset_2e is called unconditionally. The reason is 
c     that the first time datoms is called Hartree-Fock is the only 
c     thing that makes sense, in all subsequent calls the DFT terms 
c     should be taken into account. However I think it is a good idea to
c     tidy up after oneself, so CD_reset_2e is always called.
c
      dimension q(*)
      character*(*) mode
      character *7 fnm
      character *5 snm
      data fnm/"guess.m"/
      data snm/"denat"/
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/atmol3)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/sector)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/zorac)
INCLUDE(common/nshel)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
      integer idftbl
      save idftbl
      parameter(idftsc=422)
INCLUDE(common/datgue)
INCLUDE(common/restrl)
INCLUDE(common/restri)
INCLUDE(common/symtry)
c
c     allow for maximum of nbig (was 100) bfns on any given atom
c
c     nbig=1500
c     no = 50
      parameter (no=50)
c
      atmode = mode
      if (.not.ozora) atmode = 'flop'
      if (uhfatom.and.zscftp.ne.'uhf') 
     1    call caserr('atoms uhf only for uhf')
c...  calculate the maximum size of nbig 
c...  note in /junk/ it is still a parameter which must be bigger
c
      call datoms_pre(nbig)
c
      if (mode.eq.'start'.or.(mode.eq.'zora'.and..not.ozora)) then
         if (.not.oatdo) then
          if (oatmdft) then 
           write (iwr,6012)
           else
           write (iwr,6010)
          endif
         else
          write (iwr,6011)
         endif
      end if
c..
c..    core partitioning (like hcore + extra-s)
c..
c     out = nprint.eq.2
      oyes = .true.
      ono  = .not.oyes
      l1 = num
      l2 = num*(num+1)/2
      l3 = num*num
c
      ntr = nbig*(nbig+1)/2
      nsq = nbig * nbig
c
      i10 = igmem_alloc_inf(ntr,'guess.m','denat','i10',IGMEM_NORMAL)
      i15 = igmem_alloc(l2)
      i20 = igmem_alloc(l2)
      if (uhfatom) i20b = igmem_alloc(l2)
      if (oatdo) then 
c...   for atorbs 
         i21 = igmem_alloc(l3)
         i22 = igmem_alloc(l1)
         i23 = igmem_alloc(l1)
      end if
c...
      ipcap = igmem_alloc(ntr)
      iqcap = igmem_alloc(ntr)
      ifc = igmem_alloc(ntr)
      ifo = igmem_alloc(ntr)
      is = igmem_alloc(ntr)
      iu = igmem_alloc(ntr)
      it = igmem_alloc(ntr)
      ih = igmem_alloc(ntr)
      idc = igmem_alloc(ntr)
      idos = igmem_alloc(ntr)
      if (uhfatom) then
         idalpha  = igmem_alloc(ntr)
         idbeta  = igmem_alloc(ntr)
      end if
      idt = igmem_alloc(ntr)
      idold = igmem_alloc(ntr)
      iss = igmem_alloc(ntr)
      ic = igmem_alloc(nsq)
      icopn = igmem_alloc(nsq)
      ismin = igmem_alloc(nbig*no)
      iqmin = igmem_alloc(nbig*no)
      itransf = igmem_alloc(nsq)
      icc = igmem_alloc(nsq)
      iccc  = igmem_alloc(nbig*6)
      iiloc = igmem_alloc(nbig*6)
      iisch = igmem_alloc(nbig*6)
      iodon = igmem_alloc(maxat)
c
_IF(ccpdft)
      odft = CD_active().and.oatmdft
_ELSE
      odft = .false.
_ENDIF
      if (odft) then
         ikma = igmem_alloc(l2+1)
         if (mode.eq.'zora'.and.ozora) then
            call secput(isect(idftsc),2001,lensec(l2+1),idftbl)
         endif
      endif
      if (lpseud.gt.0) then
c
crz    pseudo corrections in i15
c..    **** not symmetry adapted ****
c      calculate 1-centre pseudopotential corrections in i15
c
        if (oint_zora) then
           iiso = igmem_alloc_inf(nwiso_z,fnm,snm,'iso',IGMEM_DEBUG)
           call rdedx(q(iiso),nwiso_z,ibiso_z,num8)
        else
           iiso = igmem_alloc_inf(nw196(5),fnm,snm,'iso',IGMEM_DEBUG)
           call rdedx(q(iiso),nw196(5),ibl196(5),idaf)
        end if
        iscr = igmem_alloc_inf(l2,fnm,snm,'scratch',IGMEM_NORMAL)
        if (lpseud.eq.1) then
           call vclr(q(iscr),1,l2)
           call ecpint(q,q(i15),q(iscr),q(iiso),nshell,oyes,ono)
        else if (lpseud.eq.2) then
           nav = lenwrd()
           call vclr(q(i15),1,l2)
           isc1 = igmem_alloc_inf(l3,fnm,snm,'isc1',IGMEM_NORMAL)
           isc2 = igmem_alloc_inf(255*num,fnm,snm,'isc2',IGMEM_DEBUG)
           isc3 = igmem_alloc_inf(255*20,fnm,snm,'isc3',IGMEM_DEBUG)
           isc4 = igmem_alloc_inf((nat*nt+nav-1)/nav,fnm,snm,'isc4',
     *                        IGMEM_DEBUG)
           call xpsnlc(q(i15),q(isc1),q(isc2),q(iscr),q(isc3),q(isc4),
     *                 q(iiso),l2,nat,num,nshell,oyes)
           call gmem_free_inf(isc4,fnm,snm,'isc4')
           call gmem_free_inf(isc3,fnm,snm,'isc3')
           call gmem_free_inf(isc2,fnm,snm,'isc2')
           call gmem_free_inf(isc1,fnm,snm,'isc1')
        else 
           call caserr2("denat: illegal lspeud value")
        endif
        call sym1e(q(i15),q(iscr),q(iiso),1,nshell)
        ibl7ec = ibl7la
        call wrt3(q(i15),l2,ibl7ec,num8)
        ibl7la = iposun(num8)
        call gmem_free_inf(iscr,fnm,snm,'scratch')
        call gmem_free_inf(iiso,fnm,snm,'iso')
      else
        call vclr(q(i15),1,l2)
      endif
c
c...  add 1 electron zora corrections
c
      if (ozora) then
         if (nat_z.gt.0.and.mode.eq.'start') then
            call zora(q,q,q(i15),'read')
         else
            call zora(q,q,q(i15),'start')
         end if
         if (q(i15).ne.q(i15)) 
     1      call caserr('ZORA gives NAN; dependency?')
      end if 
_IF(ccpdft)
      if (odft) then
         if (mode.eq.'start'.and.ozora) then
            call rdedx(q(ikma),l2+1,idftbl,numdu)
            call daxpy(l2,1.0d0,q(ikma),1,q(i15),1)
            edft = q(ikma+l2)
            ierror = CD_set_2e()
         else
            call vclr(q(ikma),1,l2+1)
            edft = 0.0d0
         endif
      endif
_ENDIF
c..
c..   now loop over the atoms / do atomic scf and gather d-matrix
c..
      iagain = 0
      toteng_save = 0.0d0
10    call datoms(q(i15),q(i10),q(i20),q(i20b),
     +    q(i21),q(i22),q(i23),oprint(45),
     +    q(ipcap), q(iqcap), q(ifc), q(ifo), q(is), q(iu), 
     +    q(it), q(ih), q(idc), q(idos), q(idt), q(idold), q(iss) ,
     +    q(ic), q(icopn), q(ismin), q(iqmin), q(itransf), q(icc) ,
     +    q(iccc), q(iiloc), q(iisch), q(iodon),q(idalpha),q(idbeta), 
     +    nbig,toteng,q) 
_IF(ccpdft)
      if (odft) then
         toteng = toteng - tracep(q(i20),q(ikma),l1) + edft
      endif
_ENDIF
      write(iwr,6050)toteng
c
      if (mode.eq.'zora') then
c
c..   if atomic zora is specified converge atomic (nat_z times)
c
c..   clear zora corrections and recalculate 
c
         iagain = iagain + 1
         if (ozora.and.odft) then
            write(iwr,6070)  iagain,abs(toteng-toteng_save)
         else if (ozora.and..not.odft) then
            write(iwr,6040) iagain,abs(toteng-toteng_save)
         else if (.not.ozora.and.odft) then
            write(iwr,6060) iagain
         else
c           do nothing
         endif
      end if
c
      if ((iagain.lt.nat_z.or.abs(toteng-toteng_save).gt.critat_z)
     1   .and.mode.eq.'zora') then
         toteng_save = toteng
c
         if (lpseud.gt.0) then
            call rdedx(q(i15),l2,ibl7ec,num8)
         else
            call vclr(q(i15),1,l2)
         endif
c
         if (ozora) then
            call zora(q,q(i20),q(i15),'force')
         endif
_IF(ccpdft)
c
c...     fiddle in the DFT contributions
c
         if(odft)then
            ierror = CD_set_2e()
c
c...        Store current spintyp setting and Coulomb fitting settings
c
            orks  = CD_is_rks()
            ojfit = CD_is_jfiton()
            omem  = CD_is_jfitmem()
c
c...        Make sure we are running closed shell DFT here
c...        and without any Coulomb fitting stuff
c
            ierror = CD_rks()
            ierror = CD_jfitoff()
            call vclr(q(ikma),1,l2)
_IF(debug_S)
            isma = igmem_alloc(l2)
_ENDIF
c
c           Suppress the check on the quadrature accuracy otherwise
c           the code will barf in case of any charged species.
c
            oignore = CD_ignore_accuracy()
            idum = CD_set_ignore_accuracy(.true.)
            idum = CD_energy_ao(c,q(ikma),dum,q(i20),dum,
     +             edft,q,q,.false.,dft_accu,iwr
_IF(debug_S)
     +            ,q(isma)
_ENDIF
     +             )
            idum = CD_set_ignore_accuracy(oignore)
_IF(debug_S)
            call gmem_free_inf(isma,"guess.m","denat",'isma')
_ENDIF
            q(ikma+l2) = edft
            if (ozora) then
               call wrt3(q(ikma),l2+1,idftbl,numdu)
            endif
            call daxpy(l2,1.0d0,q(ikma),1,q(i15),1)
c
c...        Restore original spintyp setting and Coulomb fit settings
c
            if (ojfit)     ierror=CD_jfiton(omem)
            if (.not.orks) ierror=CD_uks()
         endif
_ENDIF
c
         go to 10
c
      end if
_IF(ccpdft)
      ierror = CD_reset_2e()
_ENDIF
      if (odft) then
         call gmem_free_inf(ikma,"guess.m","denat",'ikma')
         if (mode.eq.'start'.and.ozora) then
            call secdrp(isect(idftsc))
         endif
      endif
c
c..   write d-matrix to disc before it is too late
c
      if (mode.eq.'start'.or.(mode.eq.'zora'.and..not.ozora)) then
         call wrt3(q(i20),l2,ibl3pa,idaf)
         if (uhfatom) call wrt3(q(i20b),l2,ibl3pb,idaf)
c
c...   perhaps also the orbitals to the vectors section
c
         if (oatdo) then
            moutt = mouta
            if (isecat.ne.-1) moutt = isecat
            zcom(5) = 'atorbs'
            if (oatdo) call putq(zcom,ztitle,q(i22),q(i23),l1,l1,l1,
     +                           1,1,q(i21),moutt,iblkq)
         end if
      end if
c
c..   print if requested
c
      if (oprint(45).and.
     +   (mode.eq.'start'.or.(mode.eq.'zora'.and..not.ozora)).and.
     +   .not.oatdo) then
         if (.not.uhfatom) then
            write (iwr,6020) 'total'
            call prtril(q(i20),l1)
         else
            write (iwr,6020) 'alpha'
            call prtril(q(i20),l1)
            write (iwr,6020) ' beta'
            call prtril(q(i20b),l1)
         end if
      end if
c
      if (oscalatz.and.mode.eq.'zora') then
         call atscf(q,i20)
         write(iwr,6030) 
      endif
c
      if (mode.eq.'zora') then
c...     dump atomic zora corrections to the dumpfile (for restarts)
         call zora(q,dum,dum,'dump')
      end if
c
c...  reset core
c
c     call gmem_free_inf(isma,"guess.m","denat",'isma')
      call gmem_free_set(i10,iodon)
c
      return
 6010 format (/' initial guess orbitals generated by ',
     +        'superposition of atomic densities'/)
 6011 format (/' initial guess orbitals are ',
     +        'atomic orbitals'/)
 6012 format (/' initial guess orbitals generated by ',
     +        'superposition of atomic dft densities'/)
 6020 format (//30x,28('-'),/,30x,'   initial guess ',a5,' density   ',
     +        /,30x,28('-')//)
 6030 format (1x,'***** adding ZORA scaling ***** ') 
 6040 format(1x,' ZORA Atomic iteration',i3,' delta-E ',1pe8.2)
 6050 format (1x,'***** total atomic energy ',f17.8,' *** ')
 6060 format(1x,' DFT Atomic iteration',i3)
 6070 format(1x,' ZORA-DFT Atomic iteration',i3,' delta-E ',1pe8.2)
      end
**==datoms_pre.f
      subroutine datoms_pre(nbig)
c
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
c
c...   subroutine to prepare  atom-scf calls for atomic startup
c...   by finding out the max. dimension to be encountered
c
INCLUDE(common/sizes)
c
c...   commons where information comes from :
INCLUDE(common/infoa)
INCLUDE(common/nshel)
c
c...  loop over atoms like datoms does
c...  we do all atoms (not unique ones like the real datoms)
c
      nbig = 0
      do 110 iat = 1 , nat
c
c ... eliminate magic/ghost/point charges/dummies
c
         if (nuct(iat).ne.1) go to 110
c
c...  gather  shell / symmetry info
c..   kng contains # primitives
c
         nbas = 0
         do 50 ii = 1 , nshell
            i = katom(ii)
            if (katom(ii).eq.iat) then
               mini = kmin(ii)
               maxi = kmax(ii)
               do 40 iorb = mini , maxi
                  if (iorb.eq.1.or.
     1                iorb.eq.2 .or. iorb.eq.5 .or. iorb.eq.11)
     2               nbas = nbas + kng(ii)
 40            continue
            end if
 50      continue
c..
      nbig = max(nbig,nbas)
c..
 110  continue
c..
      return
      end
**==denmad.f
      subroutine denmad(d,c,ns,nb,occ,nrow,occlst)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     make actual density matrices.
c     because of charge spec. highest shell
c.......................................................................
      dimension d(*), c(nb,nb)
      icount = 0
      do 40 i = 1 , nb
         do 30 j = 1 , i
            icount = icount + 1
            klast = nrow + ns - 1
            sum = 0.0d0
            do 20 k = nrow , klast-1
               sum = sum + c(i,k)*c(j,k)
 20         continue
            sum = sum*occ
            if (klast.ge.1) sum = sum + c(i,klast)*c(j,klast) * occlst
            if (i.ne.j) sum = 2*sum
            d(icount) = sum
 30      continue
 40   continue
      return
      end
**==atscf.f
      subroutine atscf(q,idens)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...   perform a atoms only 1 iteration scf using denat density
c...   used to create scaled ZORA corrections
c...   partially stolen from denascf (without ga)
c...   SF, JvL (1999)
c
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/zorac)
INCLUDE(common/iofile)
INCLUDE(common/mapper)
INCLUDE(common/machin)
INCLUDE(common/harmon)
INCLUDE(common/restri)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
c
      dimension q(*)
c
      l0 = num
      l1 = num
      l2 = l1*(l1+1)/2
      l3 = num*num
c
      i10 = igmem_alloc(6*l2+l1)
      ifock = i10
      i20 = i10 + l2
      i30 = i20 + l2
      i40 = i30 + l2
      i50 = i40 + l2
      i60 = i50 + l2
      i70 = i60 + l2
      last = i70 + l1
c
c...    make sure we produce 2J-K
c
      call swap_zora('small')
c
      isexch_zs = isexch_z
      isexch_z = 1
      oatint_z = .true.
      oatscf_z = .true.
c
      call standv(1,q)
      call coul_z(q,idens,ifock)
c
      isexch_z = isexch_zs
      oatint_z = .false.
c
      call swap_zora('nonsmall')
c
c...  read one electron ints
c
      call rdedx(q(i20),l2,ibshat_z+lensec(l2),num8)
c
      if (ozora) then
            call zora(q,q(idens),q(i20),'read')
      end if
c
      call vadd(q(i10),1,q(i20),1,q(i10),1,l2)
c
c      make orthormalizing vectors at q(i30)
c..    stolen from denscf (nov99)
c..
c..    first get canonical orthonormal vectors like hcore
c..    read s-matrix in i20, vectors end up at i30 
c..    nno symmetry here, we are dealing with atoms
c..    *note* s-matrix is blocked
c
c     ----- read in fock transformation matrix -----
c           transform hamiltonian matrix
c
c       at q(i10) fock matrix (for later)
c       at q(i20) overlapmatrix
c       at q(i30) orthonormalised vectors, output
c       at q(i70) and q(i50) scratch
c
c...  symmetry adapt s-ints (mainly for harmonic)
c
      call rdedx(q(i30),l2,ibshat_z,num8)
      call tranp(q(i30),q(i20))
c
      call jacobi(q(i20),iky,l1,q(i30),ilifq,l1,q(i70),2,2,1.0d-10)
c
c     ----- form canonical (normalizing) orthonormal orbitals -----
c
      l0 = 0
      do j=1,l1
         if (dabs(q(i70+j-1)).gt.1.0d-8) then
            l0 = l0 + 1
            call dscal(l1,1.0d0/dsqrt(q(i70+j-1)),q(i30+(j-1)*l1),1)
            if (j.ne.l0)
     +        call dcopy(l1,q(i30+(j-1)*l1),1,q(i30+(l0-1)*l1),1)
         end if
      end do
      if (l0.ne.l1) then
         if (l0.ne.newbas0) call caserr2('atom harmon problem in atscf')
         call vclr(q(i30+l0*l1),1,(l1-l0)*l1)
      end if
c
c     - q - at q(i30) orthonormalizing transformation vectors
c     - h - at q(i10)
c     - h'- at q(i50) transformed h matrix
c     (both h and h' require space of square)
c
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
      call vclr(q(i50),1,l2)
      call mult2(q(i30),q(i50),q(i10),l1,l1,l1)
c
c     ----- diagonalize new hamiltonian matrix -----
c
c      at q(i50) transformed fock matrix
c      at q(i10) eigenvectors
c      at q(i70) eigenvalues
c
      call jacobi(q(i50),iky,l1,q(i10),ilifq,l1,q(i70),2,2,1.0d-10)
c
c     ----- back-transform the eigenvectors -----
c
c     eigenvectors at q(i10)
c     the famous orthonormalising vectors from qmat at q(i30)
c     scratch area at q(i50)
c     the vectors are still used tdowned / do not read again
c
c...     complete transformation
c...     not really necessary / scaling can handle 0-vectors
c
c     call compv(q(i30),l0,l1,q,'add')
c
      call tfsqc(q(i10),q(i30),q(i50),l1,l1,l1)
c
c
c...  compute scaling matrix
c
      call zora(q,q(i30),q(i50),'scale')
c... build fock matrix (diagonal with orbital energies) 
      call vclr(q(i50),1,l2)
      do i=1,l1
         q(i50+iky(i+1)-1) = q(i70+i-1)
      end do
c
c... fock matrix at i50
c... scratch at i30
c... eigenvectors at i10
c... space for scaling matrix at i40 
c... scaling correction on ao basis at i40 (return)
c
      ehf = 0.0d0
      call scale_z(q,q(i50),q(i30),q(i10),q(i40),ehf,etot,ne,num,0)
c
c... now replace zora correction by scaled zora corrections
c... unscaled zora correcion + scaling correction
c
c... read zora corrections from ed7 
      call rdedx(q(i10),nwcor_z,ibcor_z,num8) 
c... add scaling correction
      call vadd(q(i10),1,q(i40),1,q(i10),1,l2)
c... replace unscaled zora correction
      call wrt3(q(i10),nwcor_z,ibcor_z,num8) 
c
      oatscf_z = .false.
c
      call gmem_free_inf(i10,"guess.m","atscf",'i10')
c
      return
      end
**==denscf.f
      subroutine denscf(q,zscf)
c
c     performs one cycle scf with the density matrix from
c     routine denat
c     suited for open shells and direct scf
c     works like rhfclm or drhfcl
c     corrections added for casscf,mcscf etc etc etc
c     renate zwaans 1991
c     modified to use 6, rather than 10, triangles (jun 1996)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*),nsymm(8),msymm(8)
c
INCLUDE(common/sizes)
c
INCLUDE(common/machin)
INCLUDE(common/tran)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
INCLUDE(common/gjs)
INCLUDE(common/atmol3)
INCLUDE(common/cslosc)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/diisd/sta(210),cca(20),ra(20),scalea(20),iposa(20),
     + nstora,mpa,odiisa,
     * nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     * mpb,odiisb,nspacb,nsss(30),igvbo(maxorb),igvb1(maxorb),igsp(4),
     * ekk(63),intci(150)
      common/scra/iso(mxshel,48)
INCLUDE(common/mapper)
INCLUDE(common/datgue)
INCLUDE(common/nshel)
INCLUDE(common/prints)
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/runlab)
INCLUDE(common/scfopt)
INCLUDE(common/scra7)
INCLUDE(common/segm)
INCLUDE(common/symtry)
INCLUDE(common/scfwfn)
INCLUDE(common/timeperiods)
INCLUDE(common/harmon)
INCLUDE(common/zorac)
_IF(ga)
INCLUDE(common/gadiis)
_ENDIF
_IF(parallel)
INCLUDE(common/parcntl)
INCLUDE(common/nodeio)
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/statis)
_ENDIF
INCLUDE(common/gmempara)
INCLUDE(common/restri)
c
      dimension zcas(2)

      Integer n_irreps
      Logical use_symmetry, force_serial

_IF(ga)
      character*1 xn
      data xn/'n'/
_ENDIF

      data zuhf/'uhf'/, zgvb/'gvb'/, zgrhf/'grhf'/
      data zcas/'casscf','mcscf'/
      data m51,small/51,1.0d-3/
_IF(parallel)
      data m3,m167 /3,167/
_ENDIF
      inull = igmem_null()
      call check_feature('denscf')
cjvl      call gmem_check_guards('in denscf') 
c
c..
c..   core partitioning (like hcore + extra-s)
c..   core for hstar seprate (i10/iexch/idmat)
c
      call start_time_period(TP_DENSCF)

      iter = 0
      nav = lenwrd()
      omcscf = zscf.eq.zcas(1) .or. zscf.eq.zcas(2)
      ouhf = zscf.eq.zuhf
      ogvb = zscf.eq.zgvb. or. zscf.eq.zgrhf
      out = nprint.eq.2
      l1 = num
      l2 = num*(num+1)/2
      lentri = l2
      l3 = num*num
      l4 = l2 + l2
      l0 = newbas0
c
      ifock  = inull
      irdmat = inull
      iprefa = inull
      iexch  = inull
      idmat  = inull
c
c.... for casscf/mcscf and harmonic reset section 490 etc
c.... because the adapt already compressed it
c.... this is only need for atomic startup in (mc/cas)scf
c.... and needs to be undone at end of denscf
c
      if (omcscf.and.oharm) call expharm(flop,'ctrans',flop)
c
_IF(parallel)
      oswed3(4) = .true.
      oswed3(8) = .true.
_ENDIF

      if (oatdo.and.isecat.eq.-1) then
c...     for atorbs most is not needed
c...     i10 orbs i50 backtranformed orbs / i40 occupations / i70 eps
         i10 = igmem_alloc(2*l3+2*l1)
         i50 = i10 + l3
         i70 = i50 + l3
         i40 = i70 + l1
         ndum = nprint
         nprint = -5
         call getq(q(i10),q(i70),q(i40),l1,l1,3,ieig,ipop,mouta,
     +             'natorb')
         nprint = ndum
c...     set noc1
         noc1 = 0
         do i=1,l1
            if (q(i40+i-1).gt.0.0d0) noc1 = noc1 + 1
         end do
         go to 10
      end if
c
c...  Begin constructing the core partitioning:
c
      if (odscf) then
         nss = 1
         i10 = igmem_alloc_inf(l2*6+l1,'guess.m','denscf',
     +                         'i10-ifock',IGMEM_NORMAL)
         ifock = i10
c
         i20 = i10 + l2
         iexch = i20
         i30 = i20 + l2
         idmat = i30
c
         i40 = i30 + l2
         irdmat = i40
         iprefa = irdmat + ikyp(nshell)
         i50 = i40 + l2
         i60 = i50 + l2
         i70 = i60 + l2
c
      else
c
_IF(unicos,convex,titan)
         loc10=loccm()
         lwor=nmaxly-loc10
         len=l4+2*l1
_IF(unicos)
         nss=max(2,min(10,(lwor-len)/l4))
_ELSEIF(convex)
         nss=max(2,min(18,(lwor-len)/l4))
_ELSEIF(titan)
         nss=max(2,min(5,(lwor-len)/l4))
_ENDIF
c
_ELSE
         nss = 2
_ENDIF
c
         i10 = igmem_alloc_inf(max(l2*6+l1,nss*14+l2),'guess.m',
     +                         'denscf','i10-ifock',IGMEM_NORMAL)
         iexch = i10 + l2
         idmat = i10 + nss*l4
c
         i20 = i10 + l2
         i30 = i20 + l2
         i40 = i30 + l2
         i50 = i40 + l2
         i60 = i50 + l2
         i70 = i60 + l2
c
      endif
c
_IF(ccpdft)
      ierror = CD_update_geom(c)
      if (CD_active().and.odenscfdft) then
         call retrieve_spare(imemspare)
         imemfree = igmem_max_memory() - imemspare
_IF(debug_S)
         imemfree  = imemfree - 2*l2 - 2*igmem_overhead()
_ENDIF
         imemreq  = CD_memreq_energy_ao(q,q,iwr)
         ierror = CD_jfit_init2(imemfree,imemreq,iwr)
         if (ierror.ne.0) then
            call caserr2('Out of memory in incore coulomb fit')
         endif
         if (CD_jfit_incore().and.ozora) then
            write(iwr,*)'*** WARNING: ZORA memory usage not accurately',
     +                  ' accounted for!!!'
            write(iwr,*)'*** WARNING: Calculation may run out of ',
     +                  'memory!!!'
         endif
      endif
_ENDIF
c
      omaster = .true.
_IF(parallel)
      if(ipiomode .eq. IO_NZ_S .and. .not. opg_root())
     &     omaster = .false.
_ENDIF

      call start_time_period(TP_DENSCF_BUILD)

      if (.not.odscf) then
c
c        read d-matrix from disc at q(idmat)
c
         call rdedx(q(idmat),l2,ibl3pa,idaf)
      else
         call start_time_period( TP_DENSCF_RDMAT )
         m171t=171
         len171=  lensec(l2)
         nshtri=ikyp(nshell)
         nshblk=lensec(nshtri)
         if (ouhf) then
            ntri = 8
            lent = len171 + len171
         else if (ogvb) then
            ntri = ncores + npair + npair + nseto
            lent = len171 * ntri
         else
            lent = len171
            ntri = 4
         endif
         iof171(1)= 0
         iof171(2)= iof171(1) + nshblk
         iof171(3)= iof171(2) + nshblk
         do  i=4,6
           iof171(i)=iof171(i-1)+lent
         enddo
         if(irest.eq.1.or.irest.eq.2) then
            call secget(isect(471),m171t,ibl171)
            write(iwr,2990)irest,ibl171
 2990       format(/1x,'restarting denscf in integrals :',
     *      'restart parameter =',i3/
     *      1x,'section 171 at block ',i5)
            call rdedx(q(iprefa),nshtri,ibl171,idaf)
            call reads(q(irdmat),nshtri,idaf)
            call rdedx(q(ifock),l2,ibl171+iof171(3),idaf)
            call rdedx(q(idmat),l2,ibl171+iof171(6),idaf)
            dlnmxd=-9999999.0d0
            do 1672 kkk=1,nshtri
               if(q(irdmat+kkk-1).le.dlnmxd) goto 1672
               dlnmxd=q(irdmat+kkk-1)
 1672       continue
_IFN(parallel)
            if(irest.eq.1) goto 6351
            if(irest.eq.2) goto 6352
_ENDIF
         endif
         if (.not.odelta) ntri = 0
         ln171=ntri*len171+nshblk*2
         call secput(isect(471),m171t,ln171,ibl171)
_IF(parallel)
         if(ipiomode.eq.IO_NZ_S)then
            oswed3(4)=.false.
            oswed3(8)=.false.
         endif
_ENDIF
         if(omaster.and.irest.ne.1) then
            call rdedx(q(idmat),l2,ibl3pa,idaf)
            call rdmake(q(iprefa))
            if(out) then
               write(iwr,2980) ibl171,ln171
 2980          format(/1x,
     &         'output section 171 to block ',i5/1x,
     &         'section length              ',i5/)
            endif
            call wrt3(q(iprefa),nshtri,ibl171,idaf)
            call mkrdmt('guess',q(irdmat),q(idmat),l2,nprint)
            call wrt3s(q(irdmat),nshtri,idaf)
            if(odelta) then
              call zer171(q(i10),l2,ntri,ibl171+iof171(3),len171,idaf)
              call wrt3(q(idmat),l2,ibl171+iof171(6),idaf)
            endif
            call clredx
         endif
_IF(parallel)
         if(ipiomode.eq.IO_NZ_S)then
            call pg_brdcst(7123,q(idmat),l2*8,0)
            call pg_brdcst(7122,q(irdmat),nshtri*8,0)
            call pg_brdcst(7124,dlnmxd,8,0)
            oswed3(4) = .true.
            oswed3(8) = .true.
         endif
_ENDIF(parallel)
      end if
c
 6351 continue
      if (odscf) then
c...     perform closed-shell direct scf
         if( out ) write(iwr,3030)
 3030    format(/5x,25('*'))
         dlntol = tolitr(1) - dmin1(dmax1(dlnmxd,delfac),0.0d0)
         if ( out ) write (iwr,6010) dlntol
c ...
c        we must reset "zscftp" in /runlab/ to ensure
c        the closed shell fock builder is called from intega etc
c
         zsave = zscftp
         zscftp = 'rhf'
         call end_time_period( TP_DENSCF_RDMAT )
      endif
c
c     ----- construct a skeleton fock matrix -----
c     h       at q(i10)
c     density at q(idmat)
c     exch    at q(iexch) (scratch)
c
c     Modify fock builder options
c
_IF(ccpdft)
      if (odenscfdft) then
         idum = CD_set_2e()
      endif
_ENDIF

      call start_time_period( TP_DENSCF_INTS )

      if (odscf) then

          if(.not.uhfatom) call dhstar(q,q(ifock),q(idmat),q(iprefa),
     +                                 q(irdmat),irest)
          zscftp = zsave
          if(irest.ne.0) go to 6352
c
      else
_IF(cray,convex,titan)
          if(.not.uhfatom) call hstar(q(idmat),q(i10),q(iexch),nopk,nss)
_ELSE
          if(.not.uhfatom) call hstar(q(idmat),q(i10),q(iexch),nopk)
_ENDIF
      end if

      call end_time_period( TP_DENSCF_INTS )

c
c restore fock builder options
c
_IF(ccpdft)
      if (odenscfdft) then
         idum = CD_reset_2e()
      endif
_ENDIF
c
_IF(parallel)
      if(ipiomode.eq.IO_NZ_S)then
         oswed3(4) = .false.
         oswed3(8) = .false.
      endif
_ENDIF(parallel)
      if(omaster) then
         if (out) then
            write (iwr,6060)
            call prtril(q(i10),l1)
         end if
c
c        ----- symmetrize skeleton fock matrix -----
c
c        -h- at q(i10)
c        scratch area at q(i20)
c
         if (omcscf) then
            ntsave = nt
            nt = 1
         endif
         call symh(q(i10),q(i20),iky,0,0)
c
         if (omcscf) nt = ntsave
         if (out) then
            write (iwr,6070)
            call prtril(q(i10),l1)
         end if
c
c        ---- read in core hamiltonian matrix, one electron ints in 
c             q(i20)
c        - h - at q(i10)
c
         call rdedx(q(i20),l2,ibl7f,num8)
c
         if (ozora) then
            call zora(q,q(idmat),q(i20),'read')
         end if
c
         call vadd(q(i10),1,q(i20),1,q(i10),1,nx)
c
_IF(ccpdft)
         if(CD_active().and.odenscfdft)then
c
c...        Store current spintyp setting
c
            orks = CD_is_rks()
            oignore = CD_ignore_accuracy()
c
c...        Make sure we are running closed shell DFT here
c
            ierror = CD_rks()
c
c...        Switch the accuracy check off to stop the code barfing at
c...        all charged species.
c
            ierror = CD_set_ignore_accuracy(.true.)
_IF(debug_S)
            isma = igmem_alloc(l2)
            ismb = igmem_alloc(l2)
            call vclr(q(isma),1,l2)
_ENDIF
            idum = CD_energy_ao(c,q(i10),dum,q(idmat),dum,
     +             edft,q,q,outon,dft_accu,iwr
_IF(debug_S)
     +             ,q(isma)
_ENDIF
     +             )
_IF(debug_S)
            call compare_S(q(isma),q(ismb),l2,num)
            call gmem_free_inf(ismb,"guess.m","denscf",'ismb')
            call gmem_free_inf(isma,"guess.m","denscf",'isma')
_ENDIF
            call symm_op(q(i10))
c
c...        Restore original spintyp setting
c
            ierror = CD_set_ignore_accuracy(oignore)
            if (.not.orks) ierror=CD_uks()
            ierror = CD_jfit_clean2()
            if (ierror.ne.0) then
               call caserr2(
     +            'Memory failure in rhfcld:CD_jfit_clean2')
            endif
         endif
_ENDIF
c
         if (out) then
            write (iwr,6050)
            call prtril(q(i10),l1)
         end if

         call end_time_period(TP_DENSCF_BUILD)
c
c        make orthormalizing vectors at q(i30)
c..
c..      first get canonical orthonormal vectors like hcore
c..      read s-matrix in i20, vectors end up at i30 and on idaf 
c..      (ibl3qs) this is symmetry adapted to keep scf happy
c
c        ----- read in fock transformation matrix -----
c              transform hamiltonian matrix
c
c        at q(i20) overlapmatrix
c        at q(i30) orthonormalised vectors, output
c        at q(i70) and q(i50) scratch
c
         call start_time_period( TP_DENSCF_DIAG_S )
         
         call rdedx(q(i20),l2,ibl7st,num8)
         call secget(isect(490),m51,iblk51)
         call readi(mmmm,mach(13)*nav,iblk51,idaf)
c
         n_irreps = 0
         Do i = 1, l1
            isymmos( i ) = isymaos( i )
            n_irreps = Max( n_irreps, isymaos( i ) )
         End Do
         use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
         call qmat_symm(q,q(i20),q(i30),q(i70),q(i50),iky,l0,l1,l3,l1,
     +                  out,isymmos, use_symmetry )
c
         call wrt3(q(i30),l3,ibl7la,num8)

         call end_time_period( TP_DENSCF_DIAG_S )

c
c        - q - at q(i30) orthonormalizing transformation vectors
c        - h - at q(i10)
c        - h'- at q(i50) transformed h matrix
c        (both h and h' require space of square)
c
         if (symm_diag) then
            call characterize_mo( l1, l0, q( i30 ), isymaos, 
     +                            n_irreps, isymmos, ierr )
            If( ierr .NE. 0 ) Then
               If( opg_root() ) Then
                  Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                          'symmetry of mo''s for use in the diag.'
                  Write( 6, * ) 'Ignoring symmetry in diag'
               end if
            end if
         else
            ierr = 999
         endif
         use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and. 
     +                  symm_diag
         call start_time_period( TP_DENSCF_TDOWN )
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
         call end_time_period( TP_DENSCF_TDOWN )
         call start_time_period( TP_DENSCF_SIMIL )
_IF(ga)
         if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
            call load_ga_from_square(ih_vec,q(i30),l1)
            call load_ga_from_triangle(ih_scr2,q(i10),l1)
            call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
            call load_triangle_from_ga(q(i50),ih_scr2, l1)
         else
            call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
         endif
_ELSEIF(_AND(newscf,scalapack))
         call mult2_f90(q(i30),q(i50),q(i10),l0,l0,l1)
_ELSE
         call mult2(q(i30),q(i50),q(i10),l0,l0,l1)
_ENDIF
         call end_time_period( TP_DENSCF_SIMIL )
c
c        ----- diagonalize new hamiltonian matrix -----
c
c         at q(i50) transformed fock matrix
c         at q(i10) eigenvectors
c         at q(i70) eigenvalues
c
         m2 = 2
c
c        WARNING: current setting of 1.0d-8 could cause potential 
c        problems in some high symmetry cases .. consider resetting 
c        diaacc to as low as 1.0d-11 (the lower bound at SCF 
c        convergence)
c
         diaacc = 1.0d-8

c        broadcast fock to all nodes (screened I/O and parallel diag only)
         if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &       call pg_brdcst(7127,q(i50),l2*8,0)
c
         call start_time_period( TP_DENSCF_DIAG )
         call start_time_period( TP_DIAG )
         call jacobi_symm(q(i50),iky,l0,q(i10),ilifq,l0,q(i70),
     +        m2,2,diaacc, isymmos, use_symmetry)
         call end_time_period( TP_DIAG )
         call end_time_period( TP_DENSCF_DIAG )
c
c        ----- back-transform the eigenvectors -----
c
c        eigenvectors at q(i10)
c        the famous orthonormalising vectors from qmat at q(i30)
c        scratch area at q(i50)
c
         call start_time_period( TP_DENSCF_BACK )
         call rdedx(q(i30),l3,ibl7la,num8)
_IF(ga)
         if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then
 
            call load_ga_from_square(ih_scr,q(i10),l1)
            call load_ga_from_square(ih_vec,q(i30),l1)
 
            call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &           ih_vec, ih_scr, 0.0d0, ih_scr2)
 
            call load_square_from_ga(q(i10),ih_scr2, l1)
 
         else
            call tfsqc(q(i10),q(i30),q(i50),l0,l1,l1)
         endif
_ELSEIF(_AND(newscf,scalapack))
         Do i = 1, l0
            Do j = 1, l0
               q( i10 + ( i - 1 ) * l0 + j - 1 ) = 
     +              q( i10 + ( i - 1 ) * l1 + j - 1 )
            End Do
         End Do
         call tfsqc_f90(l0,l1,q(i30),q(i10))
_ELSE
         call tfsqc(q(i10),q(i30),q(i50),l0,l1,l1)
_ENDIF
         call end_time_period( TP_DENSCF_BACK )

_IF(newscf)
c        Not used as newscf_f90 no longer uses DENSCF
c        but left as is for testing purposes
         call rdedx(q(i30),l3,ibl7la,num8)
_ENDIF
      endif
c
c...  re-entry oatdo
c     need i10 (orbitals) i70 (e's) i40 (occupations)
c
 10   continue
      if (omaster) then
c
c...     save mo-s in case of uhf
c
         if (ouhf) then
c...        save orbitals + orbital-energies
            call wrt3(q(i10),l3,ibl7la,num8)
            call dcopy(l1,q(i70),1,q(i30),1)
         end if
c
c...     select mo-s
c
         noc = na
         nsav = mouta
         mswap = nswapa
c
c...     swap orbitals if requested
c        align orbitals for mcscf / force     pjar
c
         if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
            write(iwr,*) 'mc swap call from denscf'
            call swap_mcf(q(i10),mswap,nswap)
         endif
         if (mswap.gt.0) call swap(q(i10),q(i70),num,mswap,nswap)
         if (oatdo.and.isecat.eq.-1) then
c..         swap occupations along and do not replace them for atorbs
            do i=1,mswap
               temp = q(i40-i+nswap(i))
               q(i40-i+nswap(i)) = q(i40-i+nswap(i+1))
               q(i40-i+nswap(i+1)) = temp
            end do
         else  
c
c..         set up occupations in i40 like hcore
c
            call llvmo(q(i70),q(i40),noc,noc1,l0)
            pop = 2.0d0
            if (ouhf) pop = 1.0d0
            call dscal(l0,pop,q(i40),1)
         end if
c
c...     save ** mo-s and eigenvalues / orbital occupations **
c...            (i10)      (i70)             (i40)
c
         call putq(zcom,ztitle,q(i70),q(i40),l1,l1,l0,1,1,q(i10),nsav,
     +             ibl3qa)
c
c ...    this analysis is required for casscf/mcscf incorporating atoms
c
         if (oprint(45).and..not.uhfatom) write(iwr,6080)
         if (.not.otran) then
            call secget(isect(490),m51,iblk51)
            call readi(mmmm,mach(13)*nav,iblk51,idaf)
            do 30 i=1,l0
               ii=i10+ilifq(i)
               ibig=idamax(l1,q(ii),1)
               isymmos(i)=isymaos(ibig)
               bigg = 0.0d0
               ibig = 0
               do 50 j=1,l1
                  if (isymmos(i).ne.isymaos(j).and.
     +                dabs(q(ii+j-1)).gt.bigg) then
                     bigg = dabs(q(ii+j-1))
                     ibig = j
                  end if
 50            continue
               if (bigg.gt.small) write(iwr,6090) i,ibig,bigg
 30         continue
c
c           final check    aos versus mos
c
            do 70 i=1,8
               msymm(i)=0
               nsymm(i)=0
 70         continue
c
c           build up ao array
c
            do 80 i=1,l1
               j=isymaos(i)
	       if (j .gt. 0) then
		   nsymm(j)=nsymm(j)+1
	       endif
 80         continue
c
c           build up mo array
c
            do 90 i=1,l1
               j=isymmos(i)
	       if (j .gt. 0) then
		   msymm(j)=msymm(j)+1
	       endif
 90         continue
c
c           now cross check
c
            do 170 i=1,8
               if(msymm(i).eq.nsymm(i))goto170
               if ((oharm.or.odepen).and.msymm(i).eq.nsym0(i)) go to 170
               write(iwr,6100)i,nsymm(i),msymm(i)
 170        continue
         endif
         if (oprint(45).and..not.uhfatom) then
            do 100 i=1,l0
               if (otran) then
                  write(iwr,6110)i,q(i70+i-1),q(i40+i-1)
               else
                  write(iwr,6120)i,isymmos(i),q(i70+i-1),q(i40+i-1)
               endif
 100        continue
            write(iwr,6130)
         endif
c
c        revised section 190 with mo symmetries
c
         if(.not.otran)call wrt3i(mmmm,mach(13)*nav,iblk51,idaf)
         call setsto(1360,0,mmmm)
c
c        = end of analysis =
c
c..      transform back to original basis and print if requested
c
         call start_time_period( TP_DENSCF_TDOWN )
         call tdown(q(i50),ilifq,q(i10),ilifq,l0)
         call end_time_period( TP_DENSCF_TDOWN )
c
         if (oprint(45).and..not.uhfatom) then
            if (ouhf) write (iwr,6030)
            write (iwr,6020)
            call prev(q(i50),q(i70),l0,l1,l1)
         end if
c
c...     now do density matrices
c...     i10 resulting d-matrix
c...     i50 orbitals
c...     i40 occupations
c...     i70 eigenvalues
c
         call start_time_period( TP_DENSCF_MAKE_DENS )
         if (.not.uhfatom) then
            call dmtx(q(i10),q(i50),q(i40),iky,noc1,l1,l1)
            call wrt3(q(i10),l2,ibl3pa,idaf)
         end if
         call wrt3(q(i70),l1,ibl3ea,idaf)
         call end_time_period( TP_DENSCF_MAKE_DENS )

         if (ouhf) then
c
c...        if uhf similar for beta orbitals
c
c...        restore orbitals / energies
            call rdedx(q(i10),l3,ibl7la,num8)
            call dcopy(l1,q(i30),1,q(i70),1)
c
c...        select mo-s
c
            noc = nb
            nsav = moutb
            mswap = nswapb
            ibase = 40
c
c...        swap orbitals if requested
c
            if (mswap.gt.0) call swap(q(i10),q(i70),num,mswap,
     +                                nswap(ibase+1))
c
c..         set up occupations in i40 like hcore
c..         no change for atorbs yet
c
            call llvmo(q(i70),q(i40),noc,noc1,l0)
c
c...        clear b-vectors etc. if nb=0
c
            if (nb.eq.0) call setz(q(i10),q(i10),q(i70),l1,l2,l3)
c
c...        save ** mo-s and eigenvalues / orbital occupations **
c...                (i10)      (i70)             (i40)
c
            call putq(zcom,ztitle,q(i70),q(i40),l1,l1,l0,1,1,q(i10),
     +                nsav,ibl3qb)
c
c..         transform back to original basis and print if requested
c
            call start_time_period( TP_DENSCF_TDOWN )
            call tdown(q(i50),ilifq,q(i10),ilifq,l0)
            call end_time_period( TP_DENSCF_TDOWN )
c
            if (oprint(45).and..not.uhfatom) then
               write (iwr,6040)
               write (iwr,6020)
               call prev(q(i50),q(i70),l0,l1,l1)
            end if
c
c...        now do density matrices
c...        i10 resulting d-matrix
c...        i50 orbitals
c...        i40 occupations
c...        i70 eigenvalues
c
            call start_time_period( TP_DENSCF_MAKE_DENS )
            if (.not.uhfatom) then
               call dmtx(q(i10),q(i50),q(i40),iky,noc1,l1,l1)
               call wrt3(q(i10),l2,ibl3pb,idaf)
            end if
            call wrt3(q(i70),l1,ibl3eb,idaf)
            call end_time_period( TP_DENSCF_MAKE_DENS )
         end if
      else
_IF(parallel)
c        ****  PROCESSING BY OTHER NODES 
_IF(diag_parallel)
c        receive fock from root node
         if (l1 .ge. idpdiag)then
            call pg_brdcst(7127,q(i50),l2*8,0)
            call start_time_period( TP_DENSCF_DIAG )
            call start_time_period( TP_DIAG )
            m2=2
c
c           WARNING: current setting of 1.0d-8 could cause potential 
c           problems in some high symmetry cases .. consider resetting 
c           diaacc to as low as 1.0d-11 (the lower bound at SCF 
c           convergence)
c
            diaacc = 1.0d-8
            call jacobi(q(i50),iky,l0,q(i10),ilifq,l1,q(i70),m2,2,
     +                  diaacc)
            call end_time_period( TP_DIAG )
            call end_time_period( TP_DENSCF_DIAG )
         endif
_ENDIF
c
c        now set up necessary block addresses for other nodes
c        qmat section
c
         len1=lensec(l1+1)
         lenv=len1+lensec(l2)
         lenq=lenv+lensec(l3)
         call secput(isecqm,m167,lenq,iblkv)
         ibl3qs=iblkv+lenv
c
c        vectors section
c
         len1=lensec(mach(8))
         lenv=lensec(l3)
         len2=lensec(mach(9))
         j=len2+lenv+len1+1
         call secput(mouta,m3,j,iblkt)
         ibl3qa=iblkt+len2+len1+1
         if(ouhf) then
            call secput(moutb,m3,j,iblkt)
            ibl3qb=iblkt+len2+len1+1
         endif 
_ENDIF
      endif
_IF(parallel)
c all nodes here
      if(ipiomode .eq. IO_NZ_S)then
         oswed3(4) = .true.
         oswed3(8) = .true.
      endif
_ENDIF
      irest = 0
c
c.... compress ctrans again
c
      if (omcscf.and.oharm) call comharm(flop,'ctrans',flop)
c
c...  reset core
c
 6352 continue
      call gmem_free_inf(i10,"guess.m","denscf",'i10')

      call end_time_period(TP_DENSCF)

c
      return
 6010 format (5x,'dlntol =  ',f15.10/5x,25('*'))
 6020 format (//30x,28('-'),/,30x,'   initial guess orbitals   ',/,30x,
     +        28('-')//)
 6030 format (//' ====== DENSCF UHF - Alpha =======')
 6040 format (//' ====== DENSCF UHF - Beta  =======')
 6050 format (20x,11('-')/20x,'fock matrix'/,20x,11('-'))
 6060 format (/20x,20('-')/20x,'skeleton fock matrix'/20x,20('-'))
 6070 format (/20x,23('-')/20x,'symmetrized fock matrix'/20x,
     +        23('-'))
 6080 format(/1x,47('=')/
     *' m.o.  irrep  orbital energy   orbital occupancy'/
     *1x,47('=')/)
 6090 format(' ** symmetry contamination of mo',i4,' at sabf',i4,
     *         ' of ',e12.5,'  **')
 6100 format(/
     *' *** warning -error in m.o. symmetry designation'//
     *4x,'irrep. ',i2,' no. of a.o.s =',i3/
     *14x,'no. of m.o.s =',i3/)
 6110 format(1x,i3,7x,f16.8,f20.7)
 6120 format(1x,i3,i7,f16.8,f20.7)
 6130 format(/1x,47('=')/)
      end
**==densid.f
      subroutine densid(dt,dold,dos,dalpha,dbeta,nsym,nosh,ncsh,
     x                  nccup,c,damp,nconv,nbas,nitscf,tlarge,eps)
c.......................................................................
c
c     driver routine for density matrix processing
c     on convergence dt and dold will (for uhfatom) contain a,b density
c.......................................................................
      implicit REAL  (a-h,o-z),integer   (i-n)
      dimension dt(*),dold(*),dos(*),ncsh(*),nosh(*),nccup(*),c(*),
     x          nbas(*),dalpha(*),dbeta(*),eps(*)
INCLUDE(common/sizes)   
INCLUDE(common/datgue)   
INCLUDE(common/iofile)   
      logical skipd,ochargat
      character*4 yspdf
      data yspdf/'spdf'/
c
      nstep1 = 1
      nstep2 = 1
      k = 0
      vamp1 = 1.0d0
      vamp2 = 0.0d0
      ochargat = .false.
      if (nconv.eq.0) then
         if (nitscf.gt.1) vamp1 = 1.0d0 - damp
         if (nitscf.gt.1) vamp2 = damp
      else if (iatcon.gt.0.and.chargat(max(1,iatcon)).ne.0.0d0) then
         ochargat = .true.
         write(iwr,600) chargat(iatcon)
600      format(1x,/,1x,'   for the density the charge is ',f12.6)
         if (uhfatom.and.
     +       spinat(1,i,iatcon)+spinat(2,i,iatcon).ne.0.0d0) 
     +       call caserr('charge and spin may not be mixed')
c
c...     figure out where to put the charge
c...     try  highest shell assume 1s2s2p3s3p3d..... order
c...     preference for open sherlls though
c...     if not enough electrons in a shell the user should adqpt conf
c...     too many electrons we will just do
c
         nshigh = 0
         do i=1,nsym
            if (ncsh(i)+nosh(i)+i-1.gt.nshigh) then
               nshigh = ncsh(i)+nosh(i)+i-1
            end if
         end do
c...     found highest shell,  now go for highest subshell,prefer open shells
         ishigh = 0
         do i=nsym,1,-1
            if (nosh(i).ne.0) then
               ishigh = i
               go to  10
            else if (ncsh(i)+i-1.eq.nshigh) then
               if (ishigh.eq.0) ishigh = i
            end if
         end do
c...     found highest subshell ; # electrons in it ...
10       if (nosh(ishigh).ne.0) then
            nelhigh = nccup(ishigh) 
         else
            nelhigh = 4*ishigh - 2
         end if
c...     check if we can do this
         if (nelhigh.lt.chargat(iatcon)) then
            write(iwr,605) 
     +       nshigh,yspdf(ishigh:ishigh),nelhigh,chargat(iatcon)
605          format(/,' *** Shell ',i1,a1,' has',i2,' electrons, so a',
     +               '  charge of',f11.6,'is not acceptable',
     +               ' - try conf ***')
            call caserr(' too much charge in atdens ')
         end if
c...     determine occupation for the shell involved
         occhigh = nelhigh - chargat(iatcon)
      end if
c
      skipd = .false.
      if (atmode.eq.'zora'.and.iatcon.ne.0.and.nconv.ne.0) then
         if (nonrel(abs(iatcon))) then
c...        skip non-relativistic atoms to save time
            skipd = .true.
         end if
      end if
c
      do 40 i = 1 , nsym
         occucl = 4*i - 2
         occuch = occucl
         occuop = nccup(i)
         if (ochargat) then
            if (ishigh.eq.i) then
c...           charge usually on open shells
               if (nosh(i).ne.0) then
                  occuop = occhigh
               else
                  occuch = occhigh
               end if
            end if
            if (occuop.ne.nccup(i)) 
     +      write(iwr,601)'  open ',yspdf(i:i),nccup(i)*1.0d0,occuop
601         format(4x,a7,a1,'-shell from ',f8.4,' to occupation ',f8.4)
            if (occuch.ne.occucl) write(iwr,601) 'closed ',yspdf(i:i),
     +                                            occucl,occuch
         end if
c
         nbas1 = nbas(i)
         do 30 m = 1 , nbas1
            do 20 n = 1 , m
               k = k + 1
               dt(k) = 0.0d0
               dos(k) = 0.0d0
               if (uhfatom) then
                  dalpha(k) = 0.0d0
                  dbeta(k) = 0.0d0
               end if
 20         continue
 30      continue
c
         if (skipd) go to 39
         if (ncsh(i).ne.0) call denmad(dt(nstep1),c(nstep2),ncsh(i),
     +                                 nbas1,occucl,1,occuch)
         if (nosh(i).ne.0) then
            if (atmode.ne.'zora'.and.uhfatom.and.iatcon.gt.0.and.
     +          nconv.ne.0.and.
     +          spinat(1,i,iatcon)+spinat(2,i,iatcon).ne.0.0d0) then
               write(iwr,602) yspdf(i:i),
     +                        spinat(1,i,iatcon),spinat(2,i,iatcon)
602            format(6x,a1,'-shell gets ',f8.4,' alpha and  ',f8.4,
     +               ' beta character')
               if (spinat(1,i,iatcon)+spinat(2,i,iatcon).ne.occuop)
     +            write(iwr,603) 
603               format(4x,'this is *not* the specified open shell ',
     +                      'occupation')
                  call denmad(dalpha(nstep1),c(nstep2),nosh(i),
     +                        nbas1,spinat(1,i,iatcon),ncsh(i)+1,
     +                        spinat(1,i,iatcon))
                  call denmad(dbeta(nstep1),c(nstep2),nosh(i),
     +                        nbas1,spinat(2,i,iatcon),ncsh(i)+1,
     +                        spinat(2,i,iatcon))
            else
                  call denmad(dos(nstep1),c(nstep2),nosh(i),
     +                        nbas1,occuop,ncsh(i)+1,occuop)
            end if
         end if
 39      nstep1 = nstep1 + nbas1*(nbas1+1)/2
         nstep2 = nstep2 + nbas1**2
 40   continue
c...
      if (nconv.ne.0) then
c...     density may be different => comparing agains old meaningless
         icount = 0
         do i = 1 , nsym
            do j = 1 , nbas(i)
               do k = 1 , j
                  icount = icount + 1
                  dt(icount) = (dt(icount)+dos(icount))
               end do
            end do
         end do
         if (uhfatom) then
c...        generate uhf matrices
            do i=1,icount
               dt(i) = dt(i)/2.0d0
               dos(i) = dt(i) + dbeta(i)
               dt(i) = dt(i) + dalpha(i)
            end do
         end if
c...
         return
      else 
c
         tlarge = 0.0d0
         icount = 0
         do 70 i = 1 , nsym
            do 60 j = 1 , nbas(i)
               do 50 k = 1 , j
                  icount = icount + 1
                  dt(icount) = (dt(icount)+dos(icount))
     +                         *vamp1 + dold(icount)*vamp2
                  ddiff = dabs(dt(icount)-dold(icount))
                  dold(icount) = dt(icount)
                  if (ddiff.gt.tlarge) then
c                    jmax = j
c                    kmax = k
                     tlarge = ddiff
                  end if
 50            continue
 60         continue
 70      continue
      end if
c
      return
c1000 format(' ',2e22.10,i7,10x,e8.2,2i4)
      end
**==extra.f
      subroutine extra(q,zscf,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/datgue)
INCLUDE(common/atmol3)
INCLUDE(common/runlab)
INCLUDE(common/discc)
INCLUDE(common/tran)
INCLUDE(common/scra7)
INCLUDE(common/harmon)
      dimension ond(maxorb),q(*)
      dimension zcas(2)
      equivalence (ond(1),alphas(1))
_IF(ga)
INCLUDE(common/parcntl)
INCLUDE(common/gadiis)
      character*1 xn
      data xn/'n'/
_ENDIF
c     data dzero /0.0d0/
      data done,two /1.0d0,2.0d0/
      data m0,m1,m3/0,1,3/
      data zuhf,zanam,zbname/'uhf',' -a-',' -b-'/
      data zcas/'casscf','mcscf'/
c
      write (iwr,6050)
      out = nprint.eq.2
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      ibl7qa = ibl7la
c
      ibl7qb = ibl7qa + lensec(l3)
_IF1(cu)      call setsto(l1,.false.,ond)
_IFN1(cu)      call setstl(l1,.false.,ond)
      if (mextra.ne.0) then
         write (iwr,6010)
         write (iwr,6020) (next(i),i=1,mextra)
         do 20 i = 1 , mextra
            if (ond(next(i)))
     +          call caserr2('lcbf specified twice in extra directive')
            ond(next(i)) = .true.
 20      continue
         if (.not.(otran)) then
c
c ----- now handle case when adapt specified
c
_IF1(cu)            call setsto(l1,.false.,ond)
_IFN1(cu)            call setstl(l1,.false.,ond)
c
c ----- original numbering refers to ao basis set
c ----- now convert to sabf numbering
c
            do 50 i = 1 , mextra
c
               j = next(i)
               do 40 k = 1 , l1
                  n = ntran(k)
                  nsp = ilifc(k)
                  do 30 l = 1 , n
                     nsp = nsp + 1
                     if (itran(nsp).eq.j) then
                        if (.not.(ond(k))) then
                           next(i) = k
                           ond(k) = .true.
                           go to 50
                        end if
                     end if
 30               continue
 40            continue
               call caserr2(
     +         'invalid function specified in extra directive')
 50         continue
            write (iwr,6030)
            write (iwr,6020) (next(i),i=1,mextra)
         end if
      else
         write (iwr,6040)
      end if

c
c     ----- set pointers for partitioning of core -----
c
      i10 = igmem_alloc(l2+l2+l2+l2+l1+l3)
      i20 = i10 + l2
      i21 = i10 + l3
      i30 = i20 + l2
      i40 = i30 + l2
      i41 = i30 + l3
      i60 = i40 + l2
      i70 = i60 + l1
      last = i70 + l3
c
c     ----- get canonical orbitals -----
c
c     -s- at x(i20)
c     -q- at x(i30)
c
      if (out) write (iwr,6080) i10 , i20 , i21 , i30 , i40 , i41 , 
     +                          i60 , i70 , last
      call rdedx(q(i20),l2,ibl7st,num8)
      call qmat(q,q(i20),q(i30),q(i41),q(i70),iky,l0,l1,l3,l1,out)
      ipass = 1
      mblq = ibl7qa
      mswap = nswapa
      ibase = 0
      call filec(yed(numg(1)),iblk3g(1),iunit,irep)
      if (irep.ne.0) then
         call caserr2('invalid data set specification')
      end if
 60   numg(1) = iunit
      call secini(iblk3g(1),numg(1))
c
      call getq(q(i10),q(i41),q(i41),nbas,newb,m3,ieig,ipop,isecg(1),
     +          zanam)
      ncc = num - nbas
      if (mextra.gt.ncc) call caserr2(
     +    'incorrect number of lcbf specified in extra directive')
 70   call vclr(q(i30),1,l3)
      do 90 i = 1 , num
         m = ilifq(i)
         n = (i-1)*nbas
         do 80 k = 1 , num
            if (.not.(ond(k))) then
               q(i30-1+m+k) = q(i10+n)
               n = n + 1
            end if
 80      continue
 90   continue
      if (mextra.ne.0) then
         do 100 i = 1 , mextra
            q(i30-1+ilifq(i+newb)+next(i)) = 1.0d0
 100     continue
      end if
c
c     ----- orthonormalize the orbitals -----
c
c     -s- at x(i10)
c     -q- at x(i10)
c     -v- at x(i30)
c
      call ortho1(q(i10),q(i10),q(i30),q(i41),iky,l0,l0,l1,l2,l3,l1)
      call rdedx(q(i10),l3,ibl3qs,idaf)
c
c     ----- back-transform the mo-s -----
c
c     -q- at x(i10)
c     -v- at x(i30)
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i10),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i60),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i60),l0,l1,l1)
_ENDIF
c
c     align orbitals for mcscf / force     pjar
c
      if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
         if (ibase.ne.0) call caserr2('extra error with mcscf / force')
         write(iwr,*) 'mc swap call from extra'
         call swap_mcf(q(i30),mswap,nswap)
      endif
      if (mswap.gt.0) call swap(q(i30),q(i41),num,mswap,nswap(ibase+1))
c
      call wrt3(q(i30),l3,mblq,num8)
      if (zscf.ne.zuhf) then
         call secini(ibl3d,idaf)
c
c     ----- rhf case -----
c
         if (oprint(45)) write (iwr,6060)
         pop = two
         call vclr(q(i41),1,l1)
         do 110 i = 1 , na
            if (i.gt.nb) pop = done
            q(i-1+i41) = pop
 110     continue
         call rdedx(q(i30),l3,ibl7qb,num8)
         call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),
     +             mouta,ibl3qa)
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
         if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
         call dmtx(q(i10),q(i30),q(i41),iky,na,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pa,idaf)
         call wrt3(q(i41),l1,ibl3ea,idaf)
         if (na.ne.nb.and.ibl3pb.ne.0) then
            if (nb.ne.0) call dmtx(q(i20),q(i30),q(i41),iky,nb,l1,l1)
            if (nb.eq.0) call setz(q(i30),q(i10),q(i41),l1,l2,l3)
            call vsub(q(i10),1,q(i20),1,q(i10),1,l2)
            call wrt3(q(i10),l2,ibl3pb,idaf)
         end if
         go to 120
      else if (nb.eq.0) then
         call setz(q(i30),q(i10),q(i21),l1,l2,l3)
         call wrt3(q(i30),l3,ibl7qb,num8)
      else if (ipass.ne.2) then
         ipass = 2
         mblq = ibl7qa
         mswap = nswapb
         ibase = 40
         call filec(yed(numg(2)),iblk3g(2),iunit,irep)
         if (irep.ne.0) then
            call caserr2('invalid data set specification')
            go to 60
         else
            numg(2) = iunit
            call secini(iblk3g(2),numg(2))
            call getq(q(i10),q(i41),q(i41),nbasb,newb,m3,ieig,ipop,
     +                isecg(2),zbname)
            if (nbasb.ne.nbas) call caserr2(
     +          'invalid eigenvectors retrieved by extra directive')
            go to 70
         end if
      end if
      call secini(ibl3d,idaf)
      if (oprint(45)) write (iwr,6060)
c
c     ----- uhf case -----
c
_IFN1(civ)      call vfill(done,q(i41),1,na)
_IF1(c)      call setsto(na,done,q(i41))
_IF1(iv)      call setstr(na,done,q(i41))
      if (na.lt.l1) call vclr(q(i41+na),1,l1-na)
      call rdedx(q(i30),l3,ibl7qb,num8)
      call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),mouta,
     +          ibl3qa)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
      call dmtx(q(i10),q(i30),q(i41),iky,na,l1,l1)
      if (out) call prtril(q(i10),l1)
      call wrt3(q(i10),l2,ibl3pa,idaf)
c
      call wrt3(q(i41),l1,ibl3ea,idaf)
      call vclr(q(i41),1,l1)
      if (nb.ne.0) then
_IFN1(civ)      call vfill(done,q(i41),1,nb)
_IF1(c)      call setsto(nb,done,q(i41))
_IF1(iv)      call setstr(nb,done,q(i41))
         if (oprint(45)) write (iwr,6070)
         call rdedx(q(i30),l3,ibl7qb,num8)
         call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),
     +             moutb,ibl3qb)
         call tdown(q(i30),ilifq,q(i30),ilifq,l1)
         if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
         call dmtx(q(i10),q(i30),q(i41),iky,nb,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pb,idaf)
c
         call wrt3(q(i41),l1,ibl3eb,idaf)
      end if
c
c     ----- reset core memory -----
c
 120  call gmem_free_inf(i10,"guess.m","extra",'i10')
      return
 6010 format (/' list of extra basis functions (ao numbering)')
 6020 format (/30i4)
 6030 format (//' list of extra sabf')
 6040 format (/' no extra basis functions specified')
 6050 format (/' initial guess orbitals read from foreign dumpfile',
     +        /' with basis set extension'//)
 6060 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')//)
 6070 format (//30x,28('-')/30x,'initial guess beta  orbitals'/30x,
     +        28('-')//)
 6080 format (' core assignement '/,' i10, i20, i21, i30, i40,',
     +        ' i41, i60, i70  = '/8i8/' last = ',i8)
      end
**==gesext.f
      subroutine gesext(q,zscf,zconff,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
      common/bufb/iulim(maxat),llim(maxat)
INCLUDE(common/scra7)
INCLUDE(common/mapper)
INCLUDE(common/machin)
INCLUDE(common/runlab)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
      common/blkin/eneg(11)
INCLUDE(common/restri)
_IF(ga)
INCLUDE(common/parcntl)
INCLUDE(common/gadiis)
_ENDIF
      Logical use_symmetry, force_serial
      Integer n_irreps
      dimension q(*)
      dimension iblkq(6)
      dimension lneg(29),ival1(10),ival2(10),icore(25)
      dimension zcas(2)
      equivalence (iblkq(1),ibl3qa)
_IF(ga)
      character*1 xn
      data xn/'n'/
_ENDIF
      data dzero,done,two /0.0d0,1.0d0,2.0d0/
      data m1,m51/1,51/
      data zuhf,zalter /'uhf'    ,'alter'  /
      data zcas/'casscf','mcscf'/
c
      write (iwr,6020)
      out = nprint.eq.2
      nav = lenwrd()
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      onelec = .false.
      ncore = 0
      nval = 0
c
c     ----- get core memory
c
      i10 = igmem_alloc(l3 + l1 + l3 + l1 +l3 + l1 + l1)
      i11 = i10 + l3
      i20 = i11 + l1
      i22 = i20 + l2
      i21 = i20 + l3
      i30 = i21 + l1
      i31 = i30 + l3
      i40 = i31 + l1
      last = i40 + l1
c
      if (out) write (iwr,6060) i10 , i11 , i20 , i21 , i30 , i31 , 
     +                          i40 , last
c
c     ----- get canonical orthonormal vectors -----
c
c     ----- set iulim and llim -----
c
      last = 0
      do i = 1 , nat
         if (nuct(i).le.3) then
            do j = 1 , nshell
               if (katom(j).eq.i) then
                  llim(i) = kloc(j)
                  if (last.gt.0) iulim(last) = llim(i) - 1
                  go to 30
               end if
            enddo
 30         last = i
         end if
      enddo
      iulim(last) = l1
c
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
       n_irreps = 0
       Do i = 1, l1
          isymmos( i ) = isymaos( i )
          n_irreps = Max( n_irreps, isymaos( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
c     ----- read in overlap matrix (sabf)
c
      call rdedx(q(i20),l2,ibl7st,num8)
c
c     ----- form the q matrix -----
c
c     call qmat(q,q(i10),q(i30),q(i31),q(i20),iky,l0,l1,l3,l1,out)
      call qmat_symm(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out,
     +              isymmos, use_symmetry)
c
      call wrt3(q(i30),l3,ibl7la,num8)
c
c     ----- now read in the overlap matrix and form huckel matrix
c
      call rdedx(q(i10),l2,ibl7s,num8)
c
      i101 = i10 - 1
c
c           construct huckel matrix -q(i10)-, in the overlap matrix
c                       (assumptions - h,he have valence s,
c                       alkali/alkali earths have valence s only
c                       main block has valence s,p
c                       transition metals have valence s,p,d)
c
      do n = 1 , nat
         if (nuct(n).le.3) then
            atnum = czan(n)
            nucz = nint(atnum)
            call hucklp(iwr,nucz,eneg,ncore,nval,onelec,.true.,icore,
     +                  ival1,ival2,lneg)
c
c                 zero off diagonal interactions with core
c                 basis functions, and set diagonal core energies
c
            if (ncore.ne.0) then
               i0 = llim(n) - 1
               do i = 1 , ncore
                  irow = i0 + icore(i)
                  do j = 1 , l1
                     ij = iky(max(j,irow)) + min(irow,j) + i101
                     q(ij) = dzero
                  enddo
                  ii = ikyp(irow) + i101
                  q(ii) = eneg(lneg(i))
               enddo
            end if
c
c       set diagonal elements of valence orbitals
c            some minimal and split valence bases define
c            p's on alkali, some don't. these elements are
c            always treated as split-valence for li,be ,
c            na,mg  otherwise minimal s, with any p's
c            or outer s handled as supplemental functions
c
            if (nval.ne.0) then
               i0 = llim(n) - 1
               fac1 = 0.5d00
               fac2 = 0.5d00
               if (onelec) fac1 = done
               if (nval.eq.1) then
                  irow = i0 + ival1(1)
                  ii = ikyp(irow) + i101
                  q(ii) = fac1*eneg(lneg(ncore+1))
                  if (.not.onelec) then
                     irow = i0 + ival2(1)
                     ii = ikyp(irow) + i101
                     q(ii) = fac2*eneg(lneg(ncore+1))
                  end if
               else
                  do i = 1 , nval
                     irow = i0 + ival1(i)
                     ii = ikyp(irow) + i101
                     q(ii) = fac1*eneg(lneg(ncore+i))
                     if (.not.(onelec)) then
                        irow = i0 + ival2(i)
                        ii = ikyp(irow) + i101
                        q(ii) = fac2*eneg(lneg(ncore+i))
                     end if
                  enddo
               end if
            end if
c
c                 bond,rydberg,or polarization functions
c
            if (.not.onelec) nval = nval + nval
            nextra = iulim(n) - llim(n) + 1 - ncore - nval
            if (nextra.lt.0) call caserr2('invalid indexing in gesext')
            if (nextra.ne.0) then
               small = -0.005d+00
               tiny = 0.0001d+00*small
               i0 = llim(n) + ncore + nval - 1
               ii = (i0*i0+i0)/2
               do i = 1 , nextra
                  ii = ii + i + i0
                  value = small + tiny*(10*nucz+i)
                  q(ii+i101) = value
               enddo
            end if
         end if
      enddo
c
c           generate off diagonal huckel elements
c
      fudge = 0.875d+00
      if (l1.ne.1) then
         ii = 1
         ij = 1
         do i = 2 , l1
            ii = ii + i
            hii = q(i101+ii)
            jj = 0
            do j = 1 , i - 1
               jj = jj + j
               ij = ij + 1
               sij = q(i101+ij)
               if (sij.ne.dzero) then
                  q(i101+ij) = fudge*sij*(hii+q(i101+jj))
               end if
            enddo
            ij = ij + 1
         enddo
      end if
c
c     eliminate s contaminant from d orbitals
c
      do n = 1 , nat
         atnum = czan(n)
         nucz = nint(atnum)
         if (nucz.gt.20 .and. nuct(n).le.3) then
            call hucklp(iwr,nucz,eneg,ncore,nval,onelec,.true.,icore,
     +                  ival1,ival2,lneg)
            i0 = llim(n)
            eval = eneg(6)
            if (nucz.le.36) then
               if (nucz.le.30) eval = fac1*eneg(6)
               i0 = i0 + 16
            else
               i0 = i0 + 20
            end if
            call trn5d(q(i10),l2,i0,eval)
            if (nucz.le.30) then
               i0 = i0 + 6
               eval = fac2*eneg(6)
               call trn5d(q(i10),l2,i0,eval)
            end if
            if (nucz.gt.38) then
               i0 = llim(n) + 26
               eval = eneg(9)
               if (nucz.le.48) eval = fac1*eneg(9)
               call trn5d(q(i10),l2,i0,eval)
               if (nucz.le.48) then
                  i0 = i0 + 6
                  eval = fac2*eneg(9)
                  call trn5d(q(i10),l2,i0,eval)
               end if
            end if
         end if
      enddo
c
      if (oprint(45)) then
         write (iwr,6010)
         call prtril(q(i10),l1)
      end if
c
c        - q - at q(i30) orthonormalizing transformation vectors
c        - h - at q(i10)
c        - h'- at q(i20) transformed h matrix
c        (both h and h' require space of square)
c
         if (symm_diag) then
            call characterize_mo( l1, l0, q( i30 ), isymaos,
     +                            n_irreps, isymmos, ierr )
            If( ierr .NE. 0 ) Then
               If( opg_root() ) Then
                  Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                          'symmetry of mo''s for use in the diag.'
                  Write( 6, * ) 'Ignoring symmetry in diag'
               end if
            end if
         else
            ierr = 999
         endif
         use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and.
     +                  symm_diag
c
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
c     ----- fock transformation h= q*h0*q -----
c
_IF(ga)
      if(l1.ge.idpmult2)then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i20),ih_scr2, l1)
      else
         call mult2(q(i30),q(i10),q(i20),l0,l0,l1)
      endif
_ELSEIF(_AND(newscf,scalapack))
      call mult2_f90(q(i30),q(i20),q(i10),l0,l0,l1)
_ELSE
      call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
_ENDIF
c
c     ----- diagonalize h matrix -----
c
c         at q(i20) transformed fock matrix
c         at q(i10) eigenvectors
c         at q(i40) eigenvalues
c
         m2 = 2
c
c        WARNING: current setting of 1.0d-8 could cause potential
c        problems in some high symmetry cases .. consider resetting
c        diaacc to as low as 1.0d-11 (the lower bound at SCF
c        convergence)
c
         diaacc = 1.0d-8

c        broadcast fock to all nodes (screened I/O and parallel diag only)

         if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &       call pg_brdcst(7127,q(i20),l2*8,0)
c
         call jacobi_symm(q(i20),iky,l0,q(i30),ilifq,l0,q(i40),
     +        m2,2,diaacc, isymmos, use_symmetry)
c
c        ----- back-transform the eigenvectors -----
c
c        eigenvectors at q(i30)
c        orthonormalising vectors from qmat at q(i10)
c        scratch area at q(i50)
c
         call rdedx(q(i10),l3,ibl7la,num8)
c
c     ----- back-transform the eigenvectors -----
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
      endif
_ELSEIF(_AND(newscf,scalapack))
         Do i = 1, l0
            Do j = 1, l0
               q( i30 + ( i - 1 ) * l0 + j - 1 ) =
     +              q( i30 + ( i - 1 ) * l1 + j - 1 )
            End Do
         End Do
         call tfsqc_f90(l0,l1,q(i10),q(i30))
_ELSE
      call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
_ENDIF
      iblkqq = ibl7la
      call wrt3(q(i30),l3,iblkqq,num8)
      iblkev = iposun(num8)
      call wrt3(q(i40),l1,iblkev,num8)
c
c     ----- get canonical orbitals -----
c
c     -s- at q(i20)
c     -q- at q(i30)
c
c     n0 = l0
      call rdedx(q(i20),l2,ibl7st,num8)
c     call qmat(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out)
      call qmat_symm(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out,
     +               isymmos, use_symmetry)
c
c     ----- orthonormalize the orbitals -----
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
      call rdedx(q(i30),l3,iblkqq,num8)
c load overlap
      call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i10), num, l0)
       else
          call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i20),q(i10),iky,
     *         ilifq,l0,l1,1)
       endif
_ELSE
          call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i20),q(i10),iky,
     *         ilifq,l0,l1,1)
_ENDIF
c
      call rdedx(q(i40),l1,iblkev,num8)
c
c     ----- select mo*s and occupation numbers -----
c           calculate density matrices
c
      if (zscf.eq.zuhf) call wrt3(q(i30),l3,iblkqq,num8)
      pop = two
      ndaf = 1
      if (zscf.eq.zuhf) pop = done
      ipass = 1
      noc = na
      nsav = mouta
      mswap = nswapa
      ibase = 0
      if (oprint(45)) write (iwr,6030)
 130  continue
c
c     align orbitals for mcscf / force     pjar
c
      if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
         if (ibase.ne.0) call caserr2('gesext error with mcscf / force')
         write(iwr,*) 'mc swap call from gesext'
         call swap_mcf(q(i30),mswap,nswap)
      endif
      if (mswap.gt.0) call swap(q(i30),q(i40),num,mswap,nswap(ibase+1))
      if (zconff.ne.zalter) call llvmo(q(i40),q(i21),noc,no,l0)
      call dscal(l0,pop,q(i21),1)
c
c     ----- save mo*s  + orbital energies -----
c
      if (oprint(45)) write (iwr,6040) (q(i-1+i21),i=1,l0)
      call putq(zcom,ztitle,q(i40),q(i21),l1,l1,l0,m1,m1,q(i30),nsav,
     +          iblkq(ndaf))
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
c
c ----- compute and save density matrix
c
      call dmtx(q(i10),q(i30),q(i21),iky,no,l1,l1)
      call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
      call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
      ndaf = 4
      if (zscf.eq.zuhf) then
         if (nb.eq.0) then
            nsav = moutb
            call setz(q(i30),q(i10),q(i40),l1,l2,l3)
            call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m1,m1,q(i30),
     +                nsav,iblkq(ndaf))
            call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
            call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
         else if (ipass.ne.2) then
            ipass = 2
            noc = nb
            nsav = moutb
            mswap = nswapb
c
            ibase = 40
c
            call rdedx(q(i30),l3,iblkqq,num8)
            if (oprint(45)) write (iwr,6050)
            go to 130
         end if
      end if
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,"guess.m","gesext",'i10')

      return
 6010 format (//30x,22('=')/30x,'extended huckel matrix'/30x,22('=')/)
 6020 format (/1x,'initial guess orbitals generated with extended',
     +        ' basis set option'/)
 6030 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')//)
 6040 format (15x,7f15.4)
 6050 format (//30x,28('-')/30x,'initial guess beta  orbitals'/30x,
     +        28('-')//)
 6060 format (' core assignement '/' i10, i11, i20, i21, i30, i31,',
     +        ' i40 = ',7i8/' last = ',i8)
      end
**==gesmin.f
      subroutine gesmin(q,zscf,zconff,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
      common/bufb/iulim(maxat),llim(maxat)
INCLUDE(common/scra7)
INCLUDE(common/mapper)
INCLUDE(common/machin)
INCLUDE(common/runlab)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
      common/blkin/eneg(11)
INCLUDE(common/restri)
_IF(ga)
INCLUDE(common/parcntl)
INCLUDE(common/gadiis)
_ENDIF
      Logical use_symmetry, force_serial
      Integer n_irreps
      dimension q(*)
      dimension zcas(2)
      dimension iblkq(6)
      dimension lneg(29),idum1(1)
      equivalence (iblkq(1),ibl3qa)
_IF(ga)
      character*1 xn
      data xn/'n'/
_ENDIF
c
      data dzero,done,two /0.0d0,1.0d0,2.0d0/
      data zuhf,zalter /'uhf'    ,'alter'  /
      data m1,m51/1,51/
      data zcas/'casscf','mcscf'/
c
      write (iwr,6020)
      out = nprint.eq.2
      nav = lenwrd()
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      onelec = .false.
      ncore = 0
      nval = 0
c
c     ----- get core memory
c
      i10 = igmem_alloc(l3 + l1 + l3 + l1 + l3 + l1 + l1)
      i11 = i10 + l3
      i20 = i11 + l1
      i22 = i20 + l2
      i21 = i20 + l3
      i30 = i21 + l1
      i31 = i30 + l3
      i40 = i31 + l1
      last = i40 + l1
c
c     ----- get canonical orthonormal vectors -----
c
      if (out) write (iwr,6060) i10 , i11 , i20 , i21 , i30 , i31 , 
     +                          i40 , last
c
c     ----- set iulim and llim -----
c
      last = 0
      do i = 1 , nat
         if (nuct(i).le.3) then
            do j = 1 , nshell
               if (katom(j).eq.i) then
                  llim(i) = kloc(j)
                  if (last.gt.0) iulim(last) = llim(i) - 1
                  go to 30
               end if
            enddo
 30         last = i
         end if
      enddo
      iulim(last) = l1
c
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
       n_irreps = 0
       Do i = 1, l1
          isymmos( i ) = isymaos( i )
          n_irreps = Max( n_irreps, isymaos( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
c     ----- read in overlap matrix (sabf)
c
      call rdedx(q(i20),l2,ibl7st,num8)
c
c     ----- form the q matrix -----
c
c     call qmat(q,q(i10),q(i30),q(i31),q(i20),iky,l0,l1,l3,l1,out)
      call qmat_symm(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out,
     +               isymmos, use_symmetry)
c
      call wrt3(q(i30),l3,ibl7la,num8)
c
c     ----- now read in the overlap matrix and form huckel matrix
c
      call rdedx(q(i10),l2,ibl7s,num8)
c
      i101 = i10 - 1
c
c           construct huckel matrix -q(i10)-, in the overlap matrix
c
c             (assumptions - h,he have valence s,
c             alkali/alkali earths have valence s only
c             main block has valence s,p
c             transition metals have valence d,s)
c
      do n = 1 , nat
         if (nuct(n).le.3) then
            atnum = czan(n)
            nucz = nint(atnum)
            call hucklp(iwr,nucz,eneg,ncore,nval,onelec,.false.,idum1,
     +                  idum1,idum1,lneg)
c
c                 zero off diagonal interactions with core
c                 basis functions, and set diagonal core energies
c
            if (ncore.ne.0) then
               i0 = llim(n) - 1
               do i = 1 , ncore
                  irow = i0 + i
                  do j = 1 , l1
                     ij = iky(max(j,irow)) + min(irow,j) + i101
                     q(ij) = dzero
                  enddo
                  ii = ikyp(irow) + i101
                  q(ii) = eneg(lneg(i))
               enddo
            end if
c
c                 set diagonal elements of valence orbitals
c                      some minimal valence bases define
c                      p-s on alkali, some do not. these elements are
c                      always treated as minimal s, with any p-s
c                      or outer s handled as supplemental functions
c
            if (nval.ne.0) then
               i0 = llim(n) + ncore - 1
               fac = done
               do i = 1 , nval
                  irow = i0 + i
                  ii = ikyp(irow) + i101
                  q(ii) = fac*eneg(lneg(ncore+i))
               enddo
            end if
c
c                 bond,rydberg,or polarization functions
c
            nextra = iulim(n) - llim(n) + 1 - ncore - nval
            if (nextra.lt.0) call caserr2('invalid indexing in gesext')
            if (nextra.ne.0) then
               small = -0.005d+00
               tiny = 0.0001d+00*small
               i0 = llim(n) + ncore + nval - 1
               ii = (i0*i0+i0)/2
               do i = 1 , nextra
                  ii = ii + i + i0
                  value = small + tiny*(10*nucz+i)
                  q(ii+i101) = value
               enddo
            end if
         end if
      enddo
c
c           generate off diagonal huckel elements
c
      fudge = 1.75d+00/2.0d+00
      if (l1.ne.1) then
         ii = 1
         ij = 1
         do i = 2 , l1
            ii = ii + i
            hii = q(i101+ii)
            jj = 0
            do j = 1 , i - 1
               jj = jj + j
               ij = ij + 1
               sij = q(i101+ij)
               if (sij.ne.dzero) then
                  q(i101+ij) = fudge*sij*(hii+q(i101+jj))
               end if
            enddo
            ij = ij + 1
         enddo
      end if
c
c           eliminate s contaminant from d orbitals
c
      do n = 1 , nat
         atnum = czan(n)
         nucz = nint(atnum)
         if (nucz.gt.20 .and. nuct(n).le.3) then
            call hucklp(iwr,nucz,eneg,ncore,nval,onelec,.false.,idum1,
     +                  idum1,idum1,lneg)
            i0 = llim(n) - 1 + 9
            eval = eneg(6)
            call trn5d(q(i10),l2,i0,eval)
            if (nucz.gt.38) then
               i0 = llim(n) - 1 + 19
               eval = eneg(9)
               call trn5d(q(i10),l2,i0,eval)
            end if
         end if
      enddo
c
      if (oprint(45)) then
         write (iwr,6010)
         call prtril(q(i10),l1)
      end if
c
c        - q - at q(i30) orthonormalizing transformation vectors
c        - h - at q(i10)
c        - h'- at q(i20) transformed h matrix
c        (both h and h' require space of square)
c
         if (symm_diag) then
            call characterize_mo( l1, l0, q( i30 ), isymaos,
     +                            n_irreps, isymmos, ierr )
            If( ierr .NE. 0 ) Then
               If( opg_root() ) Then
                  Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                          'symmetry of mo''s for use in the diag.'
                  Write( 6, * ) 'Ignoring symmetry in diag'
               end if
            end if
         else
            ierr = 999
         endif
         use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and.
     +                  symm_diag
c
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
c     ----- fock transformation h= q*h0*q -----
c
_IF(ga)
      if(l1.ge.idpmult2)then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i20),ih_scr2, l1)
      else
         call mult2(q(i30),q(i10),q(i20),l0,l0,l1)
      endif
_ELSEIF(_AND(newscf,scalapack))
      call mult2_f90(q(i30),q(i20),q(i10),l0,l0,l1)
_ELSE
      call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
_ENDIF
c
c     ----- diagonalize h matrix -----
c
c         at q(i20) transformed fock matrix
c         at q(i10) eigenvectors
c         at q(i40) eigenvalues
c
         m2 = 2
c
c        WARNING: current setting of 1.0d-8 could cause potential
c        problems in some high symmetry cases .. consider resetting
c        diaacc to as low as 1.0d-11 (the lower bound at SCF
c        convergence)
c
         diaacc = 1.0d-8
c
c        broadcast fock to all nodes (screened I/O and parallel diag only)
c
         if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &       call pg_brdcst(7127,q(i20),l2*8,0)
c
         call jacobi_symm(q(i20),iky,l0,q(i30),ilifq,l0,q(i40),
     +        m2,2,diaacc, isymmos, use_symmetry)
c
c        ----- back-transform the eigenvectors -----
c
c        eigenvectors at q(i30)
c        orthonormalising vectors from qmat at q(i10)
c        scratch area at q(i50)
c
         call rdedx(q(i10),l3,ibl7la,num8)
c
c     ----- back-transform the eigenvectors -----
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
      endif
_ELSEIF(_AND(newscf,scalapack))
         Do i = 1, l0
            Do j = 1, l0
               q( i30 + ( i - 1 ) * l0 + j - 1 ) =
     +              q( i30 + ( i - 1 ) * l1 + j - 1 )
            End Do
         End Do
         call tfsqc_f90(l0,l1,q(i10),q(i30))
_ELSE
      call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
_ENDIF
      iblkqq = ibl7la
      call wrt3(q(i30),l3,iblkqq,num8)
      iblkev = iposun(num8)
      call wrt3(q(i40),l1,iblkev,num8)
c
c     ----- get canonical orbitals -----
c
c     -s- at q(i20)
c     -q- at q(i30)
c
c     n0 = l0
      call rdedx(q(i20),l2,ibl7st,num8)
c     call qmat(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out)
      call qmat_symm(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out,
     +               isymmos, use_symmetry)
c
c     ----- orthonormalize the orbitals -----
c
c     -q- at q(i10)
c     -s- at q(i10)
c     -v- at q(i30)
c
      call rdedx(q(i30),l3,iblkqq,num8)
c load overlap
      call rdedx(q(i10),l2,ibl7st,num8)
_IF(ga)
       if(l1.ge.idporth .and. ipiomode .ne. IO_NZ_S)then
          call porth(q(i30),q(i10), num, l0)
       else
          call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i20),q(i10),iky,
     *         ilifq,l0,l1,1)
       endif
_ELSE
          call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
          call orfog(q(i30),q(i30),q(i20),q(i10),iky,
     *         ilifq,l0,l1,1)
_ENDIF
c
      call rdedx(q(i40),l1,iblkev,num8)
c
c     ----- select mo*s and occupation numbers -----
c           calculate density matrices
c
      if (zscf.eq.zuhf) call wrt3(q(i30),l3,iblkqq,num8)
      pop = two
      ndaf = 1
      if (zscf.eq.zuhf) pop = done
      ipass = 1
      noc = na
      nsav = mouta
      mswap = nswapa
      ibase = 0
      if (oprint(45)) write (iwr,6030)
 130  continue
c
c     align orbitals for mcscf / force     pjar
c
      if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
         if (ibase.ne.0) call caserr2('gesmin error with mcscf / force')
         write(iwr,*) 'mc swap call from gesmin'
         call swap_mcf(q(i30),mswap,nswap)
      endif
      if (mswap.gt.0) call swap(q(i30),q(i40),num,mswap,nswap(ibase+1))
      if (zconff.ne.zalter) call llvmo(q(i40),q(i21),noc,no,l0)
      call dscal(l0,pop,q(i21),1)
c
c     ----- save mo*s  + orbital energies -----
c
      if (oprint(45)) write (iwr,6040) (q(i-1+i21),i=1,l0)
      call putq(zcom,ztitle,q(i40),q(i21),l1,l1,l0,m1,m1,q(i30),nsav,
     +          iblkq(ndaf))
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
c
c ----- compute and save density matrix
c
      call dmtx(q(i10),q(i30),q(i21),iky,no,l1,l1)
      call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
      call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
      ndaf = 4
      if (zscf.eq.zuhf) then
         if (nb.eq.0) then
            nsav = moutb
            call setz(q(i30),q(i10),q(i40),l1,l2,l3)
            call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m1,m1,q(i30),
     +                nsav,iblkq(ndaf))
            call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
            call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
         else if (ipass.ne.2) then
            ipass = 2
            noc = nb
            nsav = moutb
            mswap = nswapb
c
            ibase = 40
c
            call rdedx(q(i30),l3,iblkqq,num8)
            if (oprint(45)) write (iwr,6050)
            go to 130
         end if
      end if
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,"guess.m","gesmin",'i10')

      return
 6010 format (//30x,22('=')/30x,'extended huckel matrix'/30x,22('=')/)
 6020 format (/' initial guess orbitals generated with minimal',
     +        ' basis set option'/)
 6030 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')//)
 6040 format (15x,7f15.4)
 6050 format (//30x,28('-')/30x,'initial guess beta  orbitals'/30x,
     +        28('-')//)
 6060 format (' core assignement '/' i10, i11, i20, i21, i30, i31,',
     +        ' i40 = ',7i8/' last = ',i8)
      end
**==getqq.f
      subroutine getqq(q,zscf,nprint)
c
c...  vectors getq option
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter ( mxorb1=maxorb+1)
INCLUDE(common/tran)
INCLUDE(common/machin)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/datgue)
INCLUDE(common/atmol3)
INCLUDE(common/runlab)
INCLUDE(common/harmon)
INCLUDE(common/discc)
INCLUDE(common/scra7)
      common/junk/ilifc2(maxorb),nterm2(maxorb),
     + it2(mxorb3),ct2(mxorb3),otran2,ispn2,
     + ex2(mxprim),cs2(mxprim),cp2(mxprim),cd2(mxprim),cf2(mxprim),
     + cg2(mxprim),
     + ccspac(maxat),
     + kstrt2(mxshel),ktom2(mxshel),ktype2(mxshel),kng2(mxshel),
     + kloc2(mxshel),kmin2(mxshel),kmax2(mxshel),nshel2,nspace(3)
     +,nprini(700),iovmat,iosvec,ioproj,iorthg,c2(3,maxat)
INCLUDE(common/restri)
      logical irspl
      common/junko/cspace(30),irspa(700),irspl(40),irspb(1590),
     *             irspc(336),irspd(8)
      dimension q(*)
      dimension zcas(2)
c
      data zcas/'casscf','mcscf'/
      data done,two /1.0d0,2.0d0/
c     data dzero/0.0d0/
      data m0,m1,m3,m21/0,1,3,21/
      data zuhf,zanam,zbname/'uhf',' -a-',' -b-'/
c
      write (iwr,6010)
      out = nprint.eq.2
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      nav = lenwrd()
      nw0 = maxorb + mxorb1 + 6/nav
      m1420 = 5*mxprim + maxat
      m110 = 10 + maxat
      len0 = 1 + lensec(nw0)
      lenv = lensec(l3)
      ibl7qa = ibl7la
      ibl7qb = ibl7qa + lenv
      ibl7qc = ibl7qb + lenv
c
c     ----- set pointers for partitioning of core -----
c
      i10 = igmem_alloc(l2 + l2 + l2 + l2 + l1 + l3)

      i20 = i10 + l2
      i21 = i10 + l3
      i30 = i20 + l2
      i40 = i30 + l2
      i41 = i30 + l3
      i60 = i40 + l2
      i70 = i60 + l1
      last = i70 + l3
c
c     ----- get canonical orbitals -----
c
c     -s- at x(i20)
c     -q- at x(i30)
c
      if (out) write (iwr,6040) i10 , i20 , i21 , i30 , i40 , i41 ,
     +                          i60, last
      call rdedx(q(i20),l2,ibl7st,num8)
      call qmat(q,q(i20),q(i30),q(i41),q(i70),iky,l0,l1,l3,l1,out)
c
c...  keep dumpfile in sync with memory (only and especially for redundancy)
c...  otherwise a later qmat/qmat_symm may sort it
c
      if (odepen) call revind
c
      ipass = 1
      mswap = nswapa
      ibase = 0
      mblq = ibl7qa
      call filec(yed(numg(1)),iblk3g(1),iunit,irep)
      if (irep.ne.0) then
         call caserr2('invalid data set specification')
      end if
 20   numg(1) = iunit
      call secini(iblk3g(1),numg(1))
      numab = numg(1)
c
c ---- read in vectors from secondary dumpfile
c
      call secget(isecg(1),m3,ioctr2)
      ioctr2 = ioctr2 + len0
      call getq(q(i30),q(i41),q(i41),i,j,m3,ieig,ipop,isecg(1),zanam)
      mbasis = i
      nbasis = num
c
c...  check for possible null vectors (vectors not occuring in harmonic
c...  and if so complete vector set
c...  (if from harmonic , the zero vectors are the last one (check this)
c
22    nonvec = 0
      do k = 1,mbasis
         kk = i30+(k-1)*mbasis-1
         do l=1,mbasis
            if (dabs(q(kk+l)).gt.1.0d-13) go to 25
         end do
c...     zero vector found  check if remainder is 0.0
         do l=k*mbasis,mbasis*mbasis-1
            if (q(i30+l).gt.1.0d-13) call caserr2('zero middlevec')
         end do
         go to 26
25       continue
      end do
      go to 27
c...   go to non contracted
26     continue
      call readi(ilifc2,mach(9)*nav,ioctr2,numab)
      call tdown2(q(i30),mbasis,k-1,ilifc2,nterm2,it2,ct2,otran2)
      otran2 = .true.
      do  i = 1 , mbasis
           nterm2(i) = 1
           ilifc2(i) = i - 1
           it2(i) = i
           ct2(i) = 1.0d0
       end do
      nonvec = mbasis - k +1
c...   this flag is used to avoid reading ilfc2 etc again
27    continue
c
c ---- read in geometry from the secondary dumpfile
c
 30   itmp = idaf
      call wrt3(q(i30),mbasis*mbasis,ibl7qc,num8)
      ibl7qd = iposun(num8)
      idaf = numab
c
c ----- reset coords to new ie c2=c
c
      call rdrec1(id,id,id,id,c2,c2,dum,dum,dum,c2,c2)
      idaf = itmp
      do 50 k = 1 , nat
         do 40 kk = 1 , 3
            c2(kk,k) = c(kk,k)
 40      continue
 50   continue
c
c
c ---- decide whether to force projection or not
c
      if (nonvec.eq.0) then
      call readi(ilifc2,mach(9)*nav,ioctr2,numab)
      end if
      ifprj = 0
      if (otran .and. .not.otran2) ifprj = 1
      if (otran2 .and. .not.otran) ifprj = 1
      if (mbasis.ne.nbasis) ifprj = 1
      if (ifprj.ne.1) then
         do 60 k = 1 , nat
            if (c2(1,k).ne.c(1,k)) go to 70
            if (c2(2,k).ne.c(2,k)) go to 70
            if (c2(3,k).ne.c(3,k)) go to 70
 60      continue
         ifprj = 0
      end if
      go to 80
 70   ifprj = 1
c...   always do projection for harmonic
      if (newbas1.ne.newbas0) ifprj = 1
 80   if (ifprj.ne.0) go to 120
c
c     check on =ctrans= list
c
      if (.not.(otran)) then
         do 100 i = 1 , nbasis
            n = ntran(i)
            if (n.ne.nterm2(i)) go to 110
            nsp = ilifc(i)
            if (nsp.ne.ilifc2(i)) go to 110
            do 90 j = 1 , n
               nsp = nsp + 1
               if (dabs(ctran(nsp)-ct2(nsp)).gt..001d0) go to 110
               if (itran(nsp).ne.it2(nsp)) go to 110
 90         continue
 100     continue
      end if
      go to 150
 110  ifprj = 1
c
c ---- set up shell information from secondary dumpfile
c
 120  mtype = 0
      call secget(isect(491),mtype,ibl3a)
      call rdedx(ex2,mxprim,ibl3a,numab)
      ibl3a = ibl3a + lensec(mxprim)
      call rdedx(cs2,m1420,ibl3a,numab)
      ibl3a = ibl3a + lensec(m110) + lensec(m1420)
      call readi(kstrt2,mach(2)*nav,ibl3a,numab)
c
c ---- set up restart block for secondary basis set
c
c ---- restore restart section on dumpfile
c
      call secget(isect(501),m21,ibl21)
      ibl21 = ibl21 + 1
      call rdedx(cspace,lds(isect(501)),ibl21,numab)
c
_IFN1(cuf)      call icopy(700,irspa,1,nprini,1)
_IF1(cuf)      call fmove(irspa,nprini,700)
c
c ---- allocate some more space for projection
c
      mdim = max(nbasis,mbasis)
c
c     following is attempt to reduce getqq store
c     requirements to 6-triangles.
c
c     call setscm(i81)
c     i90=i81 + mdim*mdim
c     i100=i90 +mdim*mdim
c     i110=i100+mdim
      i81 = i70
      i90 = i10
      i100 = i60
      i110 = i21
c     last=i110+mdim
c     len2=last-i81
c     call cmem(load2)
c     loc81=loccm()
c     need2=loc81+len2
c     call setc(need2)
      if (ipass.ne.2 .or. numg(1).ne.numg(2)) then
c
c ---- form the transformation matrix
c
         call rdedx(q(i30),l2,ibl7st,num8)
         call square(q(i81),q(i30),mdim,l1)
         call rootmt(q(i81),q(i90),q(i30),q(i110),q(i100),mdim,l1,1,
     #               newbas1-newbas0)
c
c ---- save transformation matrix
c
         iovmat = ibl7qd
         l4 = mdim*mdim
         call wrt3(q(i81),l4,iovmat,num8)
c
c ---- save vectors of s for completion of vectors if necessary
c
         iosvec = iposun(num8)
         call wrt3(q(i90),l4,iosvec,num8)
         ioproj = iposun(num8)
c
c ---- form the overlap matrix between the two basis sets
c
         call scross(q(i90),mdim)
c
c ---- read in ctrans for old basis
c
         if (nonvec.eq.0) then
          call readi(ilifc2,mach(9)*nav,ioctr2,numab)
         end if
c
c ---- transform overlap matrix to ctrans basis
c
         call tranpq(q(i90),q(i81),mdim,nbasis,mbasis)
c
c ---- form the projection matrix
c
         call rdedx(q(i90),l4,iovmat,num8)
         call mxmg(q(i90),1,mdim,q(i81),1,mdim,q(i30),1,mdim,nbasis,
     +             nbasis,mbasis)
c
c ---- save the projection matrix
c
         call wrt3(q(i30),l4,ioproj,num8)
         iorthg = iposun(num8)
      end if
c
c ---- now project out the vectors onto the current basis
c
      if (ifprj.ne.0) then
c
c ---- check that enough space exists and expand the vectors
c
*        if (length.le.l4) call caserr2(
*    +                      'not enough space for projection')
         call rdedx(q(i30),mbasis*mbasis,ibl7qc,num8)
         n1 = 0
         n2 = 0
         do 130 i = 1 , mbasis
            call dcopy(mbasis,q(i30+n1),1,q(i81+n2),1)
            n1 = n1 + mbasis
            n2 = n2 + mdim
 130     continue
         call projec(q(i81),q(i30),q(i90),q(i100),q(i110),mdim,nbasis,
     +               mbasis,nonvec,q)
c
c ---- realign the vector into x(i30)
c
         n1 = 0
         n2 = 0
         do 140 i = 1 , nbasis
            call dcopy(nbasis,q(i81+n2),1,q(i30+n1),1)
            n1 = n1 + nbasis
            n2 = n2 + mdim
 140     continue
      end if
c
c ---- return all extra core
c
c     call setc(load2)
c
c     ----- back-transform the mo-s -----
c
c     -q- at x(i10)
c     -v- at x(i30)
c
150   if (l0.ne.l1) then
         call rdedx(q(i10),l2,ibl7st,num8)
         call swap_qq(q(i30),q(i10),l0,l1)
      end if
      call ortho1(q(i10),q(i10),q(i30),q(i41),iky,l0,l0,l1,l2,l3,l1)
      call rdedx(q(i10),l3,ibl3qs,idaf)
      call tfsqc(q(i30),q(i10),q(i70),l0,l1,l1)
c
c     align orbitals for mcscf / force     pjar
c
      if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
         if (ibase.ne.0) call caserr2('getqq error with mcscf / force')
         write(iwr,*) 'mc swap call from getqq'
         call swap_mcf(q(i30),mswap,nswap)
      endif
      if (mswap.gt.0) call swap(q(i30),q(i41),num,mswap,nswap(ibase+1))
c
      call wrt3(q(i30),l3,mblq,num8)
      if (zscf.ne.zuhf) then
c
c     ----- rhf case -----
c
         call secini(ibl3d,idaf)
         if (oprint(45)) write (iwr,6020)
         call vclr(q(i41),1,l1)
_IFN1(civ)      call vfill(two,q(i41),1,nb)
_IF1(c)      call setsto(nb,two,q(i41))
_IF1(iv)      call setstr(nb,two,q(i41))
         if (na.gt.nb)
_IFN1(civ)     * call vfill(done,q(i41+nb),1,na-nb)
_IF1(c)     * call setsto(na-nb,done,q(i41+nb))
_IF1(iv)     * call setstr(na-nb,done,q(i41+nb))
         call rdedx(q(i30),l3,ibl7qa,num8)
         call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),
     +             mouta,ibl3qa)
         call tdown(q(i30),ilifq,q(i30),ilifq,l1)
         if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
         call dmtx(q(i10),q(i30),q(i41),iky,na,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pa,idaf)
         call wrt3(q(i41),l1,ibl3ea,idaf)
         if (na.ne.nb.and.ibl3pb.ne.0) then
            if (nb.ne.0) call dmtx(q(i20),q(i30),q(i41),iky,nb,l1,l1)
            if (nb.eq.0) call setz(q(i30),q(i10),q(i41),l1,l2,l3)
            call vsub(q(i10),1,q(i20),1,q(i10),1,l2)
            call wrt3(q(i10),l2,ibl3pb,idaf)
         end if
         go to 160
      else if (nb.eq.0) then
         call setz(q(i30),q(i10),q(i21),l1,l2,l3)
         call wrt3(q(i30),l3,ibl7qb,num8)
      else if (ipass.ne.2) then
         ipass = 2
         mblq = ibl7qb
         mswap = nswapb
         ibase = 40
         call filec(yed(numg(2)),iblk3g(2),iunit,irep)
         if (irep.ne.0) then
            call caserr2('invalid data set specification')
            go to 20
         else
            numg(2) = iunit
            call secini(iblk3g(2),numg(2))
            call secget(isecg(2),m3,ioctr2)
            ioctr2 = ioctr2 + len0
            call getq(q(i30),q(i41),q(i41),i,j,m3,ieig,ipop,isecg(2),
     +                zbname)
            numab = numg(2)
            mbasis = i
            go to 22
         end if
      end if
c
c     ----- uhf case -----
c
      call secini(ibl3d,idaf)
      if (oprint(45)) write (iwr,6020)
_IFN1(civ)      call vfill(done,q(i41),1,na)
_IF1(c)      call setsto(na,done,q(i41))
_IF1(iv)      call setstr(na,done,q(i41))
      if (na.lt.l1) call vclr(q(i41+na),1,l1-na)
      call rdedx(q(i30),l3,ibl7qa,num8)
      call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),mouta,
     +          ibl3qa)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
      call dmtx(q(i10),q(i30),q(i41),iky,na,l1,l1)
      if (out) call prtril(q(i10),l1)
      call wrt3(q(i10),l2,ibl3pa,idaf)
c
      call wrt3(q(i41),l1,ibl3ea,idaf)
      call vclr(q(i41),1,l1)
chvd  if (nb.ne.0) then
         if (oprint(45)) write (iwr,6030)
_IFN1(civ)      call vfill(done,q(i41),1,nb)
_IF1(c)      call setsto(nb,done,q(i41))
_IF1(iv)      call setstr(nb,done,q(i41))
         call rdedx(q(i30),l3,ibl7qb,num8)
         if (oprint(45)) call prev(q(i30),q(i41),l0,l1,l1)
         call putq(zcom,ztitle,q(i41),q(i41),l1,l1,l0,m0,m1,q(i30),
     +             moutb,ibl3qb)
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
         call dmtx(q(i10),q(i30),q(i41),iky,nb,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pb,idaf)
c
         call wrt3(q(i41),l1,ibl3eb,idaf)
chvd  end if
c
c     ----- reset core memory -----
c
 160  call gmem_free_inf(i10,"guess.m","getqq",'i10')
      return
 6010 format (/' initial guess orbitals read from foreign dumpfile'/)
 6020 format (//30x,28('-'),/,30x,'initial guess alpha orbitals',/,30x,
     +        28('-')//)
 6030 format (//30x,28('-'),/,30x,'initial guess beta  orbitals',/,30x,
     +        28('-')//)
 6040 format (' core assignement ',/,' i10, i20, i21, i30, i40,',
     +        ' i41, i60 = ',/,7i8,/,' last = ',i8)
      end
      subroutine swap_qq(v,s,l0,l1)
c
c     swap vectors to get 0.0 vectors at end (should be l1-l0)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension v(l1,l1),s(*)
c...   this function is little used
      itrian(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      ignore = 0
      do i=1,l0
         ss = 0.0d0
         do j=1,l1
            do k=1,l1
               ss = ss + v(j,i)*v(k,i)*s(itrian(j,k))
            end do
         end do
         if (ss.le.1.0d-12) then
            ignore = ignore + 1
c...     swap and clear the 0.0 vector properly
            call dswap(l1,v(1,i),1,v(1,l0+ignore),1)
            call vclr(v(1,l0+ignore),1,l1)
         end if
      end do
c
c     if (ignore.ne.l1-l0) call caserr2('swap_qq error')
c...   note that the vectors may be in proper order already
c
      return
      end
      subroutine tdown2(q,num,nnn,ilifc,ntran,itran,ctran,otran)
c
c...  simple tdown for use in getqq 
c...  ilifq and ilifn fixed to (i-1)*num
c...  input and output are the same
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
      dimension dum(maxorb)
c
      dimension q(*),ilifc(*),ntran(*),itran(*),ctran(*)
c
      if(otran) return
c
      do 731 i=1,nnn
      m=(i-1)*num
      call vclr(dum,1,num)
      do 733 j=1,num
      n=ntran(j)
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
      do 733 k=1,n
      l=ilifc(j)    +k
733   dum(itran(l))=ctran(l)*q(m+j)+dum(itran(l))
      call dcopy(num,dum(1),1,q((i-1)*num+1),1)
 731  continue
c 
60010 if (nnn.lt.num) then
c...    clear vectors that are extra
         do i=nnn+1,num
            call vclr(q((i-1)*num+1),1,num)
         end do
      end if
c
      return
      end
**==gnorm.f
      subroutine gnorm(ar,nbasis)
c
c          this subroutine normalizes the vector a.  the real and
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ar(nbasis)
      data done/1.0d0/
c
_IF1(u)      br=vecsum(ar(1),ar(1),nbasis)
_IF1(u)      br=done/dsqrt(br)
_IF1(u)      call scaler(nbasis,br,ar,ar)
_IFN1(u)      br=done/dnrm2(nbasis,ar,1)
_IFN1(u)      call dscal(nbasis,br,ar,1)
      return
      end
**==hamild.f
      subroutine hamild(pcap, qcap, fc, fo, s, u, t, h, dos, dt,
     +                  c, smin, qmin, nbb)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     construct fock matrices. this is a direct s c f procedure,
c     and the two-electron integrals are recalculated for every
c     iteration.
c.......................................................................
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension c(*), smin(nbb,*), qmin(nbb,*)
      dimension pcap(*), qcap(*), fc(*), fo(*), s(*)
      dimension u(*), t(*), h(*), dos(*), dt(*)
c
      call teigd (pcap, qcap, u, t, dt, dos)
c.......................................................................
c
c     compute smin and qmin
c.......................................................................
      nstep1 = 0
      nstep2 = 0
      nstep = 0
      do 50 i = 1 , nsym
         nsh = ncsh(i) + nosh(i)
         nbas1 = nbas(i)
         do 40 j = 1 , nsh
            naddr = nstep2 + (j-1)*nbas1
            j1 = nstep + j
            do 30 m = 1 , nbas1
               smin(m,j1) = 0.d0
               qmin(m,j1) = 0.d0
               do 20 n = 1 , nbas1
                  k = max(m,n)*(max(m,n)-1)/2 + min(m,n) + nstep1
                  smin(m,j1) = smin(m,j1) + s(k)*c(n+naddr)
                  qmin(m,j1) = qmin(m,j1) + qcap(k)*c(n+naddr)
 20            continue
 30         continue
 40      continue
         nstep = nstep + nsh
         nstep2 = nstep2 + nbas1**2
         nstep1 = nstep1 + n1(i)
 50   continue
c.......................................................................
c
c     compute fc and fo
c.......................................................................
      k = 1
      nstep = 0
      do 90 i = 1 , nsym
         occucl = 4*i - 2
         nosh1 = nosh(i)
         nsh = ncsh(i) + nosh1
         fact1 = occucl/(occucl-nccup(i))
         fact2 = nccup(i)/(occucl-nccup(i))
         nbas1 = nbas(i)
         do 80 m = 1 , nbas1
            do 70 n = 1 , m
               term1 = 0.d0
               term2 = 0.d0
               do 60 j = 1 , nsh
                  j1 = j + nstep
                  if (j.ne.nsh .or. nosh1.eq.0) then
                     term1 = term1 + smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  else
                     term2 = smin(m,j1)*qmin(n,j1) + qmin(m,j1)
     +                       *smin(n,j1)
                  end if
 60            continue
               fo(k) = fact1*term1 + pcap(k) + h(k) - qcap(k)
               fc(k) = fact2*term2 + pcap(k) + h(k)
               k = k + 1
 70         continue
 80      continue
         nstep = nstep + nsh
 90   continue
      return
      end
**==hcore.f
      subroutine hcore(q,zscf,oalpha,zconff,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/mapper)
INCLUDE(common/machin)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/runlab)
INCLUDE(common/scra7)
INCLUDE(common/symtry)
INCLUDE(common/harmon)
INCLUDE(common/gjs)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
_IF(ga)
INCLUDE(common/parcntl)
INCLUDE(common/gadiis)
_ENDIF
INCLUDE(common/restri)
      dimension iblkq(6),q(*)
      equivalence (iblkq(1),ibl3qa)
      Logical use_symmetry, force_serial
      Integer n_irreps
_IF(ga)
      character*1 xn
      data xn/'n'/
_ENDIF
      character*7 fnm
      character*5 snm 
      data fnm,snm/'guess.m','hcore'/

      data done,two /1.0d0,2.0d0/
      data zuhf,zalter /'uhf'    ,'alter'  /
      data m1,m51/1,51/
      if (oalpha) write (iwr,6010)
      out = nprint.eq.2
      nav = lenwrd()
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c     ----- set pointers for partitioning of core -----
c
c     The max is needed for qmat_symm because l3<l1+1 in the case 
c     when num=1.
      i10 = igmem_alloc_inf(max(l3,l1+1),fnm,snm,'i10',IGMEM_NORMAL)
      i11 = igmem_alloc_inf(l1,fnm,snm,'i11',IGMEM_DEBUG)
      i20 = igmem_alloc_inf(l3,fnm,snm,'i20',IGMEM_NORMAL)
      i21 = igmem_alloc_inf(l1,fnm,snm,'i21',IGMEM_DEBUG)
      i30 = igmem_alloc_inf(l3,fnm,snm,'i30',IGMEM_NORMAL)
      i31 = igmem_alloc_inf(l1,fnm,snm,'i31',IGMEM_DEBUG)
      i40 = igmem_alloc_inf(l1,fnm,snm,'i40',IGMEM_DEBUG)
c
c     ----- get canonical orthonormal vectors -----
c
      if (out) write (iwr,6050) i10 , i11 , i20 , i21 , i30 , i31 ,
     +                          i40, last
      call rdedx(q(i20),l2,ibl7st,num8)
c
      call secget(isect(490),m51,iblk51)
      call readi(mmmm,mach(13)*nav,iblk51,idaf)
       n_irreps = 0
       Do i = 1, l1
          isymmos( i ) = isymaos( i )
          n_irreps = Max( n_irreps, isymaos( i ) )
       End Do
       use_symmetry = n_irreps .GT. 1 .and. symm_diag
c
      call qmat_symm(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out,
     +               isymmos, use_symmetry)
c
      call wrt3(q(i30),l3,ibl7la,num8)
c
c     ----- read in bare nucleus hamiltonian matrix -----
c
      call rdedx(q(i10),l2,ibl7f,num8)
c
      if (.not.(oalpha)) then
c
c ... alphas directive has been invoked
c ... read in kinetic energy
c
         call rdedx(q(i20),l2,ibl7t,num8)
         call alphaz(q(i20),q(i10),q(i11),l1,l2)
      end if
c
c        - q - at q(i30) orthonormalizing transformation vectors
c        - h - at q(i10)
c        - h'- at q(i20) transformed h matrix
c        (both h and h' require space of square)
c
         if (symm_diag) then
            call characterize_mo( l1, l0, q( i30 ), isymaos,
     +                            n_irreps, isymmos, ierr )
            If( ierr .NE. 0 ) Then
               If( opg_root() ) Then
                  Write( 6, * ) 'WARNING: Unable to characterize ' //
     +                          'symmetry of mo''s for use in the diag.'
                  Write( 6, * ) 'Ignoring symmetry in diag'
               end if
            end if
         else
            ierr = 999
         endif
         use_symmetry = ierr .EQ. 0 .And. n_irreps .GT. 1 .and.
     +                  symm_diag
c
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
c
c     ----- fock transformation h= q*h0*q -----
c
_IF(ga)
      if(l1.ge.idpmult2.and.(ipiomode.ne.IO_NZ_S))then
         call load_ga_from_square(ih_vec,q(i30),l1)
         call load_ga_from_triangle(ih_scr2,q(i10),l1)
         call mult2_ga(ih_scr2, ih_vec, ih_scr,l1)
         call load_triangle_from_ga(q(i20),ih_scr2, l1)
      else
         call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
      endif
_ELSEIF(_AND(newscf,scalapack))
      call mult2_f90(q(i30),q(i20),q(i10),l0,l0,l1)
_ELSE
      call mult2(q(i30),q(i20),q(i10),l0,l0,l1)
_ENDIF
c
c     ----- diagonalize h matrix -----
c
c         at q(i20) transformed fock matrix
c         at q(i10) eigenvectors
c         at q(i40) eigenvalues
c
         m2 = 2
c
c        WARNING: current setting of 1.0d-8 could cause potential
c        problems in some high symmetry cases .. consider resetting
c        diaacc to as low as 1.0d-11 (the lower bound at SCF
c        convergence)
c
         diaacc = 1.0d-8

c        broadcast fock to all nodes (screened I/O and parallel diag only)
c
         if(ipiomode.eq.IO_NZ_S .and. odpdiag(l1))
     &       call pg_brdcst(7127,q(i20),l2*8,0)
c
         call jacobi_symm(q(i20),iky,l0,q(i30),ilifq,l0,q(i40),
     +        m2,2,diaacc, isymmos, use_symmetry)
c
c        ----- back-transform the eigenvectors -----
c
c        eigenvectors at q(i30)
c        orthonormalising vectors from qmat at q(i10)
c        scratch area at q(i50)
c
         call rdedx(q(i10),l3,ibl7la,num8)
c
c     ----- back-transform the eigenvectors -----
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
      endif
_ELSEIF(newscf)
         Do i = 1, l0
            Do j = 1, l0
               q( i30 + ( i - 1 ) * l0 + j - 1 ) =
     +              q( i30 + ( i - 1 ) * l1 + j - 1 )
            End Do
         End Do
         call tfsqc_f90(l0,l1,q(i10),q(i30))
_ELSE
      call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
_ENDIF
c
c     ----- select mo*s and occupation numbers
c           calculate density matrices.
c
      if (zscf.eq.zuhf) call wrt3(q(i30),l3,ibl7la,num8)
      pop = two
      ndaf = 1
      if (zscf.eq.zuhf) pop = done
      ipass = 1
      noc = na
      nsav = mouta
      mswap = nswapa
      ibase = 0
      if (oprint(45)) write (iwr,6020)
c
c ... swap orbitals if requested
c
 20   if (mswap.gt.0) call swap(q(i30),q(i40),num,mswap,nswap(ibase+1))
      if (zconff.ne.zalter) call llvmo(q(i40),q(i21),noc,no,l0)
      call dscal(l0,pop,q(i21),1)
      if (oprint(45)) write (iwr,6030) (q(i-1+i21),i=1,l0)
c
c     ----- save mo*s  + orbital energies -----
c
      call putq(zcom,ztitle,q(i40),q(i21),l1,l1,l0,m1,m1,q(i30),nsav,
     +          iblkq(ndaf))
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
c
c     ----- now density matrix ..
c
      call dmtx(q(i10),q(i30),q(i21),iky,no,l1,l1)
      call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
      call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
      ndaf = 4
      if (zscf.eq.zuhf) then
         if (nb.eq.0) then
            nsav = moutb
            call setz(q(i30),q(i10),q(i40),l1,l2,l3)
            call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m1,m1,q(i30),
     +                nsav,iblkq(ndaf))
            call wrt3(q(i10),l2,iblkq(ndaf+1),idaf)
            call wrt3(q(i40),l1,iblkq(ndaf+2),idaf)
         else if (ipass.ne.2) then
            ipass = 2
            noc = nb
            nsav = moutb
            mswap = nswapb
            ibase = 40
c
            call rdedx(q(i30),l3,ibl7la,num8)
            if (oprint(45)) write (iwr,6040)
            go to 20
         end if
      end if
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i40,fnm,snm,'i40')
      call gmem_free_inf(i31,fnm,snm,'i31')
      call gmem_free_inf(i30,fnm,snm,'i30')
      call gmem_free_inf(i21,fnm,snm,'i21')
      call gmem_free_inf(i20,fnm,snm,'i20')
      call gmem_free_inf(i11,fnm,snm,'i11')
      call gmem_free_inf(i10,fnm,snm,'i10')
      return
 6010 format (/' initial guess orbitals generated from bare nucleus',
     +        ' hamiltonian'/)
 6020 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')//)
 6030 format (15x,7f15.4)
 6040 format (/30x,28('-')/30x,'initial guess  beta orbitals'/30x,
     +        28('-')//)
 6050 format ('  core assignement'/' i10, i11, i20, i21, i30, i31,',
     +        '  i40 = ',7i8,/' last = ',i8)
      end
**==hucklp.f
      subroutine hucklp(iw,nucz,eneg,ncore,nval,onelec,odz,
     *                  icore,ival1,ival2,lneg)
      implicit REAL  (a-h,o-z)
      logical onelec,odz
      dimension eneg(*),icore(*),ival1(*),ival2(*),lneg(*)
c     dimension row1(2)
      dimension row2(3,8),row3(5,8),row4(8,8),row5(11,8)
      dimension tm1(7,10),tm2(10,10),mneg(29)
      dimension ico4(15),ival41(10),ival42(10)
      dimension ico5(25),ival51(10),ival52(10)
      dimension lrow(29),ltm(29)
c
      data mneg/1,2,3,3,3,4,5,5,5,6,6,6,6,6,6,7,8,8,8,
     *          9,9,9,9,9,9,10,11,11,11/
      data lrow/1,2,3,3,3,4,5,5,5,7,8,8,8,6,6,6,6,6,6,
     *          9,9,9,9,9,9,10,11,11,11/
      data ltm /1,2,3,3,3,4,5,5,5,7,8,8,8,6,6,6,6,6,6,
     *          10,11,11,11,9,9,9,9,9,9/
      data ico4/1,2,3,4,5,6,7,8,9,18,19,20,21,22,23/
      data ival41/10,11,12,13,18,19,20,21,22,23/
      data ival42/14,15,16,17,24,25,26,27,28,29/
c
      data ico5/1,2,3,4,5,6,7,8,9,10,11,12,13,22,23,24,25,26,27,
     *          28,29,30,31,32,33/
      data ival51/14,15,16,17,28,29,30,31,32,33/
      data ival52/18,19,20,21,34,35,36,37,38,39/
c
c     the following orbital energies are taken from
c     e.clementi,c.roetti at.nuc.data tables, vol 14.
c     the order is 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p
c     transition metal energies are from s**1 d**n states, except zn,cd
c
c     data row1/0.500d+00,0.918d+00/
      data row2/ 2.48d+00,0.196d+00,0.130d+00,
     *           4.73d+00,0.309d+00,0.190d+00,
     *           7.70d+00,0.495d+00,0.310d+00,
     *           11.3d+00,0.706d+00,0.433d+00,
     *           15.6d+00,0.945d+00,0.568d+00,
     *           20.7d+00,1.244d+00,0.632d+00,
     *           26.4d+00,1.573d+00,0.730d+00,
     *           32.8d+00,1.930d+00,0.850d+00/
      data row3/ 40.5d+00, 2.80d+00,1.52d+00,0.182d+00,0.120d+00,
     *           49.0d+00, 3.77d+00,2.28d+00,0.253d+00,0.160d+00,
     *           58.5d+00, 4.91d+00,3.22d+00,0.393d+00,0.210d+00,
     *           68.8d+00, 6.16d+00,4.26d+00,0.540d+00,0.297d+00,
     *           80.0d+00, 7.51d+00,5.40d+00,0.696d+00,0.392d+00,
     *           92.0d+00, 9.00d+00,6.68d+00,0.880d+00,0.437d+00,
     *          104.0d+00,10.6d+00,8.07d+00,1.073d+00,0.506d+00,
     *          118.6d+00,12.3d+00,9.57d+00,1.278d+00,0.591d+00/
      data row4/133.5d+00,14.5d+00,11.5d+00,1.75d+00,0.95d+00,
     *          0.000d+00,0.147d+00,0.000d+00,
     *          149.4d+00,16.8d+00,13.6d+00,2.24d+00,1.34d+00,
     *          0.000d+00,0.196d+00,0.000d+00,
     *          378.8d+00,48.2d+00,42.5d+00,6.40d+00,4.48d+00,
     *          1.193d+00,0.424d+00,0.208d+00,
     *          405.2d+00,52.1d+00,46.2d+00,7.19d+00,5.17d+00,
     *          1.635d+00,0.553d+00,0.287d+00,
     *          432.6d+00,56.3d+00,50.2d+00,8.03d+00,5.88d+00,
     *          2.113d+00,0.686d+00,0.369d+00,
     *          460.9d+00,60.7d+00,54.3d+00,8.93d+00,6.66d+00,
     *          2.650d+00,0.838d+00,0.403d+00,
     *          490.1d+00,65.2d+00,58.6d+00,9.87d+00,7.48d+00,
     *          3.220d+00,0.993d+00,0.457d+00,
     *          520.2d+00,69.9d+00,63.0d+00,10.8d+00,8.33d+00,
     *          3.825d+00,1.153d+00,0.524d+00/
      data tm1/165.8d+00,18.9d+00,15.5d+00,2.41d+00,1.44d+00,
     *          0.215d+00,0.195d+00,
     *         183.1d+00,21.2d+00,17.6d+00,2.70d+00,1.64d+00,
     *          0.273d+00,0.205d+00,
     *         201.3d+00,23.7d+00,19.8d+00,2.99d+00,1.84d+00,
     *         0.321d+00,0.214d+00,
     *         220.4d+00,26.2d+00,22.1d+00,3.29d+00,2.05d+00,
     *         0.373d+00,0.222d+00,
     *         240.4d+00,28.9d+00,24.6d+00,3.62d+00,2.30d+00,
     *         0.383d+00,0.227d+00,
     *         261.2d+00,31.7d+00,27.2d+00,3.96d+00,2.55d+00,
     *         0.406d+00,0.230d+00,
     *         282.9d+00,34.6d+00,29.9d+00,4.30d+00,2.80d+00,
     *         0.434d+00,0.233d+00,
     *         305.4d+00,37.7d+00,32.7d+00,4.65d+00,3.06d+00,
     *         0.457d+00,0.236d+00,
     *         328.8d+00,40.8d+00,35.6d+00,5.01d+00,3.32d+00,
     *         0.491d+00,0.238d+00,
     *         353.3d+00,44.4d+00,38.9d+00,5.63d+00,3.84d+00,
     *         0.783d+00,0.293d+00/
      data row5/
     *  551.0d+00, 75.0d+00, 67.9d+00,12.1d+00, 9.5d+00, 4.7d+00,
     *  1.52d+00,0.81d+00,0.000d+00,0.138d+00,0.000d+00,
     *  583.7d+00, 80.4d+00, 73.0d+00,13.5d+00,10.7d+00, 5.7d+00,
     *  1.90d+00,1.10d+00,0.000d+00,0.178d+00,0.000d+00,
     *  997.8d+00,149.4d+00,139.2d+00,29.6d+00,25.4d+00,17.6d+00,
     *  4.98d+00,3.51d+00,1.063d+00,0.372d+00,0.197d+00,
     * 1041.2d+00,157.0d+00,146.5d+00,31.6d+00,27.2d+00,19.2d+00,
     *  5.51d+00,3.97d+00,1.369d+00,0.476d+00,0.265d+00,
     * 1085.6d+00,164.8d+00,154.0d+00,33.6d+00,19.2d+00,20.8d+00,
     *  6.06d+00,4.45d+00,1.688d+00,0.582d+00,0.335d+00,
     * 1130.9d+00,172.8d+00,161.7d+00,35.8d+00,31.1d+00,22.5d+00,
     *  6.65d+00,4.95d+00,2.038d+00,0.701d+00,0.360d+00,
     * 1177.2d+00,180.9d+00,169.7d+00,37.9d+00,33.1d+00,24.3d+00,
     *  7.24d+00,5.47d+00,2.401d+00,0.821d+00,0.403d+00,
     * 1224.4d+00,189.3d+00,177.8d+00,40.2d+00,35.2d+00,26.1d+00,
     *  7.86d+00,6.01d+00,2.778d+00,0.944d+00,0.457d+00/
      data ((tm2(i,j),i=1,10),j=1,5) /
     *  616.6d+00, 85.7d+00, 78.1d+00,14.7d+00,11.8d+00, 6.5d+00,
     *  2.07d+00,1.22d+00,0.194d+00,0.192d+00,
     *  650.6d+00, 91.3d+00, 83.4d+00,15.9d+00,12.9d+00, 7.4d+00,
     *  2.30d+00,1.38d+00,0.249d+00,0.204d+00,
     *  685.4d+00, 97.0d+00, 88.8d+00,17.2d+00,14.0d+00, 8.3d+00,
     *  2.53d+00,1.55d+00,0.299d+00,0.214d+00,
     *  721.2d+00,102.9d+00, 94.5d+00,18.6d+00,15.3d+00, 9.3d+00,
     *  2.76d+00,1.72d+00,0.357d+00,0.222d+00,
     *  757.9d+00,108.9d+00,100.2d+00,20.0d+00,16.6d+00,10.3d+00,
     *  3.00d+00,1.91d+00,0.377d+00,0.222d+00/
      data ((tm2(i,j),i=1,10),j=6,10)/
     *  795.5d+00,115.2d+00,106.2d+00,21.4d+00,17.8d+00,11.3d+00,
     *  3.26d+00,2.10d+00,0.412d+00,0.222d+00,
     *  834.0d+00,121.6d+00,112.4d+00,22.9d+00,19.2d+00,12.4d+00,
     *  3.50d+00,2.29d+00,0.451d+00,0.220d+00,
     *  873.5d+00,128.1d+00,118.7d+00,24.4d+00,20.5d+00,13.5d+00,
     *  3.75d+00,2.48d+00,0.488d+00,0.220d+00,
     *  913.8d+00,134.9d+00,125.2d+00,25.9d+00,21.9d+00,14.7d+00,
     *  4.00d+00,2.68d+00,0.537d+00,0.220d+00,
     *  955.4d+00,142.1d+00,132.1d+00,27.7d+00,23.6d+00,16.1d+00,
     *  4.45d+00,3.05d+00,0.763d+00,0.265d+00/
c
c     ----- huckel parameterization -----
c
c     the valence atomic orbital energies are not good parameters.
c     the gaussian-xx program initial guess routines use numbers
c     which work better, and were obtained from experience.
c     the numbers used in gauss-xx are approximately related
c     to the atomic energies by the following factors.
c
c     eneg(1,1) = -0.537e+00
c     eneg(2,1) = -0.735e+00
c     do 300 i=3,10
c        eneg(i,2) = 2.0e+00 * eneg(i,2)
c        eneg(i,3) = 0.7e+00 * eneg(i,3)
c 300 continue
c
c
      onelec = .false.
c           bond functions
      if (nucz.le.0) then
         ncore = 0
         nval = 0
c            h
      else if (nucz.le.1) then
         ncore = 0
         nval = 1
         eneg(1) = -0.537d+00
         if (odz) then
            ival1(1) = 1
            ival2(1) = 2
         end if
c            he
      else if (nucz.le.2) then
         ncore = 0
         nval = 1
         eneg(1) = -0.735d+00
         if (odz) then
            ival1(1) = 1
            ival2(1) = 2
         end if
c           li-be
      else if (nucz.le.4) then
         ncore = 1
         if (odz) then
            nval = 4
            do 20 loop = 1 , nval
               ival1(loop) = ncore + loop
               ival2(loop) = ncore + loop + nval
 20         continue
            do 30 loop = 1 , ncore
               icore(loop) = loop
 30         continue
         else
            onelec = .true.
            nval = 1
         end if
         nuc = nucz - 2
         eneg(1) = -row2(1,nuc)
         eneg(2) = -row2(2,nuc)
         eneg(3) = -row2(3,nuc)
         eneg(2) = eneg(2) + eneg(2)
         eneg(3) = eneg(3)*0.70d+00
c            b-ne
      else if (nucz.le.10) then
         ncore = 1
         nval = 4
         nuc = nucz - 2
         eneg(1) = -row2(1,nuc)
         eneg(2) = -row2(2,nuc)
         eneg(3) = -row2(3,nuc)
c
         eneg(2) = eneg(2) + eneg(2)
         eneg(3) = eneg(3)*0.70d+00
         if (odz) then
            do 40 loop = 1 , nval
               ival1(loop) = ncore + loop
               ival2(loop) = ncore + loop + nval
 40         continue
            do 50 loop = 1 , ncore
               icore(loop) = loop
 50         continue
         end if
c           na-mg
c
c     for the second row, gauss-xx uses 2*2s, 1*2p, 2.2*3s, 2.3*3p
c     this does not predict the correct homo-lumo for sih2.
c
c     do 400 i=11,18
c        eneg(i,2) = 2.0e+00 * eneg(i,2)
c        eneg(i,4) = 4.0e+00 * eneg(i,4)
c        eneg(i,5) = 2.0e+00 * eneg(i,5)
c 400 continue
      else if (nucz.le.12) then
         ncore = 5
         nuc = nucz - 10
         if (odz) then
            nval = 4
            do 60 loop = 1 , nval
               ival1(loop) = ncore + loop
               ival2(loop) = ncore + loop + nval
 60         continue
            do 70 loop = 1 , ncore
               icore(loop) = loop
 70         continue
         else
            onelec = .true.
            nval = 1
         end if
         do 80 i = 1 , 5
            eneg(i) = -row3(i,nuc)
 80      continue
         eneg(2) = eneg(2) + eneg(2)
         eneg(4) = eneg(4)*4.0d0
         eneg(5) = eneg(5) + eneg(5)
c           al-ar
      else if (nucz.le.18) then
         ncore = 5
         nval = 4
         nuc = nucz - 10
         if (odz) then
            do 90 loop = 1 , nval
               ival1(loop) = ncore + loop
               ival2(loop) = ncore + loop + nval
 90         continue
            do 100 loop = 1 , ncore
               icore(loop) = loop
 100        continue
         end if
         do 110 i = 1 , 5
            eneg(i) = -row3(i,nuc)
 110     continue
c            k-ca
      else if (nucz.le.20) then
         ncore = 9
         nuc = nucz - 18
         do 120 i = 1 , 5
            eneg(i) = -row4(i,nuc)
 120     continue
         eneg(7) = -row4(7,nuc)
         if (odz) then
            nval = 4
            do 130 loop = 1 , nval
               ival1(loop) = ival41(loop)
               ival2(loop) = ival42(loop)
 130        continue
            do 140 loop = 1 , ncore
               icore(loop) = ico4(loop)
 140        continue
         else
            onelec = .true.
            nval = 1
         end if
c           sc-zn
      else if (nucz.le.30) then
c     for the 3rd row tm ,  use 2*3s, 1*3p, 2.0*3d, 2.0*4s
         ncore = 9
         nuc = nucz - 20
         do 150 i = 1 , 7
            eneg(i) = -tm1(i,nuc)
 150     continue
         if (odz) then
            nval = 10
            eneg(8) = eneg(7)*0.3333333d0
            do 160 loop = 1 , nval
               ival1(loop) = ival41(loop)
               ival2(loop) = ival42(loop)
 160        continue
            do 170 loop = 1 , ncore
               icore(loop) = ico4(loop)
 170        continue
         else
            nval = 7
         end if
         eneg(4) = eneg(4) + eneg(4)
         eneg(6) = eneg(6)*2.0d0
         eneg(7) = eneg(7) + eneg(7)
c           ga-kr
      else if (nucz.le.36) then
         ncore = 15
         nval = 4
         nuc = nucz - 28
         do 180 i = 1 , 8
            eneg(i) = -row4(i,nuc)
 180     continue
         if (odz) then
            do 190 loop = 1 , nval
               ival1(loop) = ival41(loop)
               ival2(loop) = ival42(loop)
 190        continue
            do 200 loop = 1 , ncore
               icore(loop) = ico4(loop)
 200        continue
         end if
c           rb-sr
      else if (nucz.le.38) then
         ncore = 19
         nuc = nucz - 36
         do 210 i = 1 , 8
            eneg(i) = -row5(i,nuc)
 210     continue
         eneg(10) = -row5(10,nuc)
         if (odz) then
            nval = 4
            do 220 loop = 1 , nval
               ival1(loop) = ival51(loop)
               ival2(loop) = ival52(loop)
 220        continue
            do 230 loop = 1 , ncore
               icore(loop) = ico5(loop)
 230        continue
         else
            onelec = .true.
            nval = 1
         end if
c            y-cd
      else if (nucz.le.48) then
         ncore = 19
         nuc = nucz - 38
         do 240 i = 1 , 10
            eneg(i) = -tm2(i,nuc)
 240     continue
         if (odz) then
            nval = 10
            eneg(11) = eneg(10)*0.3333333d0
            do 250 loop = 1 , nval
               ival1(loop) = ival51(loop)
               ival2(loop) = ival52(loop)
 250        continue
            do 260 loop = 1 , ncore
               icore(loop) = ico5(loop)
 260        continue
         else
            nval = 7
         end if
c           in-xe
      else if (nucz.le.54) then
         ncore = 25
         nval = 4
         nuc = nucz - 46
         do 270 i = 1 , 11
            eneg(i) = -row5(i,nuc)
 270     continue
         if (odz) then
            do 280 loop = 1 , nval
               ival1(loop) = ival51(loop)
               ival2(loop) = ival52(loop)
 280        continue
            do 290 loop = 1 , ncore
               icore(loop) = ico5(loop)
 290        continue
         end if
c           the rest of the elements...
      else
         write (iw,6010) nucz
         call caserr2('invalid nuclear charge for huckel guess')
      end if
      nsum = ncore + nval
      if (.not.(odz)) then
         do 300 loop = 1 , nsum
            lneg(loop) = mneg(loop)
 300     continue
         if (nucz.eq.19 .or. nucz.eq.20) lneg(10) = 7
         if (nucz.eq.37 .or. nucz.eq.38) lneg(20) = 10
      else if (nucz.le.30 .or. nucz.gt.48) then
         do 310 loop = 1 , nsum
            lneg(loop) = lrow(loop)
 310     continue
      else if (nucz.le.36) then
         do 320 loop = 1 , nsum
            lneg(loop) = mneg(loop)
 320     continue
      else
         do 330 loop = 1 , nsum
            lneg(loop) = ltm(loop)
 330     continue
      end if
      return
c
 6010 format (1x,'nucz=',i5,' is  .gt.54---huckel guess unavailable')
      end
**==jacod.f
      subroutine jacod(f,v,nb,nb1,nb2,nmin,nmax,big,jbig,maxao)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c      f is the matrix to be diagonalized. f is stored triangular
c      v is the array of eigenvectors. quadratic array, dimension nb*nb
c      big and jbig are temporary scratch areas of dimension nb
c      the rotations among the first nmin basis functions are not
c      accounted for.
c      the rotations among the last nb-nmax basis functions are not
c      accounted for.
c.......................................................................
      dimension big(nb),jbig(nb),f(nb1),v(nb2)
      data root2 /0.707106781186548d0/
c     data c1/1.d-12/
      data c2,c3,c4,c5,c6/1.d-12,4.d-16,2.d-16,8.d-9,3.d-6/
c...    ligen ignored
c     if(iligen.ne.1) goto 17171
c     call ligen(f,v,jbig,nb,maxao,big,1.e-13)
c     return
c       17171 continue
      if (nb.eq.1) then
         v(1) = 1.0d0
         return
      end if
      ii = 0
c.......................................................................
c      loop over rows (i) of triangular matrix
c.......................................................................
      do 30 i = 1 , nb
         big(i) = 0.00d0
         jbig(i) = 0
         if (i.ge.nmin .and. i.ne.1) then
            j = min(i-1,nmax)
c.......................................................................
c      loop over columns (k) of triangular matrix to determine
c      largest off-diagonal elements in row(i).
c.......................................................................
            do 20 k = 1 , j
               if (dabs(big(i)).lt.dabs(f(ii+k))) then
                  big(i) = f(ii+k)
                  jbig(i) = k
               end if
 20         continue
         end if
         ii = ii + i
 30   continue
 40   sd = 1.050d0
c.......................................................................
c      find smallest diagonal element and corresponding largest
c      off-diagonal element.
c.......................................................................
      jj = 0
      do 50 j = 1 , nb
         jj = jj + j
         sd = dmin1(sd,dabs(f(jj)))
 50   continue
      sd = dmax1(sd,c6)*c2
      t = 0.0d0
      i1 = max(2,nmin)
      ib = 1
      do 60 i = i1 , nb
         if (t.lt.dabs(big(i))) then
            t = dabs(big(i))
            ib = i
         end if
 60   continue
c.......................................................................
c      test for convergence, then determine rotation.
c.......................................................................
      ia = jbig(ib)
      if (t.lt.sd) then
         return
      else
         iaa = ia*(ia-1)/2
         ibb = ib*(ib-1)/2
         jaa = (ia-1)*maxao
         jbb = (ib-1)*maxao
         dif = f(iaa+ia) - f(ibb+ib)
         if (dabs(dif).gt.c3*t) then
            t2x2 = big(ib)/dif
            t2x25 = t2x2*t2x2
            if (t2x25.le.c4) then
               cx = 1.0d0
               sx = t2x2
            else if (t2x25.le.c5) then
               sx = t2x2*(1.0d0-1.5d0*t2x25)
               cx = 1.0d0 - 0.5d0*t2x25
            else if (t2x25.gt.c6) then
               t = 0.25d0/dsqrt(0.25d0+t2x25)
               cx = dsqrt(0.5d0+t)
               sx = dsign(dsqrt(0.5d0-t),t2x2)
            else
               cx = 1.0d0 + t2x25*(t2x25*1.375d0-0.5d0)
               sx = t2x2*(1.00d0+t2x25*(t2x25*3.875d0-1.5d0))
            end if
         else
            sx = root2
            cx = root2
         end if
         iear = iaa + 1
         iebr = ibb + 1
         do 90 ir = 1 , nb
            t = f(iear)*sx
            f(iear) = f(iear)*cx + f(iebr)*sx
            f(iebr) = t - f(iebr)*cx
            if (ir.lt.ia) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ir.eq.ia) then
               tt = f(iebr)
               ieaa = iear
               ieab = iebr
               f(iebr) = big(ib)
               iear = iear + ir - 1
               if (jbig(ir).ne.0) go to 70
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else
               t = f(iear)
               it = ia
               iear = iear + ir - 1
               if (ir.lt.ib) then
               else if (ir.eq.ib) then
                  f(ieaa) = f(ieaa)*cx + f(ieab)*sx
                  f(ieab) = tt*cx + f(iebr)*sx
                  f(iebr) = tt*sx - f(iebr)*cx
                  iebr = iebr + ir - 1
                  go to 70
               else
                  if (dabs(t).lt.dabs(f(iebr))) then
                     if (ib.le.nmax) then
                        t = f(iebr)
                        it = ib
                     end if
                  end if
                  iebr = iebr + ir - 1
               end if
            end if
            if (dabs(t).ge.dabs(big(ir))) then
               big(ir) = t
               jbig(ir) = it
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            else if (ia.ne.jbig(ir) .and. ib.ne.jbig(ir)) then
               iear = iear + 1
               iebr = iebr + 1
               go to 90
            end if
 70         kq = iear - ir - ia + 1
            big(ir) = 0.00d0
            ir1 = min(ir-1,nmax)
            do 80 i = 1 , ir1
               k = kq + i
               if (dabs(big(ir)).lt.dabs(f(k))) then
                  big(ir) = f(k)
                  jbig(ir) = i
               end if
 80         continue
            iear = iear + 1
            iebr = iebr + 1
 90      continue
         do 100 i = 1 , maxao
            t = v(jbb+i)*sx
            v(jbb+i) = v(jaa+i)*sx - v(jbb+i)*cx
            v(jaa+i) = v(jaa+i)*cx + t
 100     continue
         go to 40
      end if
      end
      subroutine fermi_smear(e,f,em,nfocc,noccmx,t,ne,nmo,iwr,out)
      implicit none
c
c     Description:
c     ------------
c
c     Calculate the occupation numbers for the orbitals using
c     Fermi-Dirac smearing:
c
c                 1
c        f_i = -------------------------, b = 1/t
c              1 + exp( b * [e_i - mu] )
c
c     mu is to be optimised such that sum_i f_i = Ne. 
c
c     Define Nd = sum_i f_i then the error is
c
c        d = (Ne - Nd) ** 2
c
c        error = sqrt(d)
c
c     The gradient of d with respect to mu is
c
c        du = -2 (Ne - Nd) dNdu
c
c                     b * exp( b * [e_i - mu] )
c        dNdu = sum_i -------------------------------
c                     ( 1 + exp( b * [e_i - mu] ))**2
c
c     The second derivative of d with respect to mu is
c
c        d2u = 2 dNdu ** 2 - 2 (Ne - Nd) d2Ndu
c
c                         ( b * exp( b * [e_i - mu] )) ** 2
c        d2dNdu = sum_i 2 ---------------------------------
c                         ( 1 + exp( b * [e_i - mu] )) ** 3
c
c                         b**2 * exp( b * [e_i - mu] )
c               - sum_i   ---------------------------------
c                         ( 1 + exp( b * [e_i - mu] )) ** 2
c
c     Now we can use a quasi newton update:
c
c         mu_i+1 = mu_i - du / d2dNdu
c
c     Start with mu = (E_lumo + E_homo)/2.
c
c
c     Also calculate the Mermin term to the total energy:
c
c        Em = 1/b * sum_i {(1-f_i)*ln(1-f_i) + f_i*ln(f_i)}
c
c     Adding this term to the total energy will give the correct
c     finite temperature energy. The gradient expressions should not
c     need any modification [1].
c
c     Comments:
c     ---------
c
c     1) The Fermi-Dirac smearing enforces strict aufbau ordering
c        of the orbitals through the occupations. Thus it cannot be
c        used together with options that may break the aufbau ordering
c        such as (large) level shifters, or locking.
c
c     2) The Fermi-Dirac smearing breaks the strict distinction 
c        between occupied and virtual orbitals. This has consequences
c        for the definition of the tester.
c
c     References:
c     -----------
c
c     [1] R.W. Warren, B.I. Dunlap
c         "Fractional occupation numbers and density functional
c          energy gradients with the linear combination of Gaussian-
c          type orbitals approach"
c         Chemical Physics Letters 262 (1996) 384-392.
c
c     Input
c
      integer nmo ! number of MO's
      integer ne  ! number of electrons
      REAL e(nmo) ! orbital energies
      REAL t      ! smearing parameter
      integer iwr ! unit number of standard output
      logical out ! .true. if output requested
c
c     Output
c
      REAL f(nmo)    ! occupation numbers
      REAL Em        ! the Mermin energy term
      integer nfocc  ! the number of fully occupied orbitals
      integer noccmx ! the number of (fully or partially) occupied 
                     ! orbitals
c
c     Functions
c
      logical opg_root
c
c     Local variables
c
      logical onotaufbau
      REAL mu
      REAL b,b2
      REAL Nd,dNdu,d2Ndu
      REAL z,z2,zp1,zp12,zp13
      REAL d,du,d2u
      REAL p,tt
      REAL Ehomo, Elumo
      integer i,j,k
      REAL tol
      parameter(tol=1.0d-8)
      integer maxit
      parameter(maxit=20)
c
      onotaufbau = .false.
      if (ne.le.0) then
         do i = 1, nmo
            f(i) = 0.0d0
         enddo
         return
      endif
c     write(*,*)'t=',t
      tt = max(t,tol)
      b  = 1.0d0/tt
      b2 = b*b
c
c     First we need to calculate the starting mu for
c     which we need the LUMO and HOMO energy levels
c     We need to do this because if level shifters are
c     used the energy levels may be in the "wrong" order.
c
      do i = 1, nmo
         f(i) = 0.0d0
      enddo
      do i = 1, Ne
         k = nmo
         do j = 1, nmo-1
            if (f(j).lt.0.5d0.and.e(j).lt.e(k)) then
               k=j
            endif
         enddo
         f(k) = 1.0d0
      enddo
      Elumo = e(nmo)
      Ehomo = e(1)
      do j = 2, nmo-1
         if (f(j).gt.0.5d0) then
            Ehomo = max(Ehomo,e(j))
         else
            Elumo = min(Elumo,e(j))
         endif
      enddo
      mu = 0.5d0 * (Ehomo + Elumo)
c
c     Now optimise mu to conserve the number of electrons
c
      if (out.and.opg_root()) then
         write(iwr,"(2(/1x,a4,1x,a5,10x,a12,8x,a6,1x,a9))")
     +        "iter","Nelec","Nelec_approx","Nerror","chem.pot.",
     +        "----","-----","------------","------","---------"
      endif
      j      = 0 
 10   j      = j+1
      Nd     = 0.0d0
      dNdu   = 0.0d0
      d2Ndu  = 0.0d0
      Em     = 0.0d0
      nfocc  = 0
      noccmx = 0
      do i = 1, nmo
         p = b*(e(i)-mu)
         if (p.gt.-dlog(tol)) then
            f(i)   = 0.0d0
         else if (p.lt.dlog(tol)) then
            f(i)   = 1.0d0
            Nd     = Nd+1.0d0
            nfocc  = nfocc+1
c           noccmx = noccmx+1
            noccmx = i
         else
            z      = dexp(p)
            z2     = z*z
            zp1    = 1.0d0+z
            zp12   = zp1*zp1
            zp13   = zp1*zp12
            Nd     = Nd + 1.0d0/zp1
            dNdu   = dNdu + b*z/zp12
            d2Ndu  = d2Ndu + 2.0d0*b2*z2/zp13 - b2*z/zp12
            f(i)   = 1.0d0/zp1
            Em     = Em + tt*((1.0d0-f(i))*dlog(1.0d0-f(i))+
     &                        f(i)*dlog(f(i)))
c           noccmx = noccmx+1
            noccmx = i
         endif
      enddo
      d   =  (Ne-Nd)
      du  = -2.0d0*(Ne-Nd)*dNdu
      d2u =  2.0d0*dNdu*dNdu-2.0d0*(Ne-Nd)*d2Ndu
      if (out.and.opg_root()) then
         write(iwr,"(1x,i4,1x,i5,f22.8,f14.8,f10.4)")j,Ne,Nd,d,mu
      endif
      if (dabs(d)/Ne.gt.tol) then
         if (j.lt.maxit) then
            mu = mu - du/d2u
            goto 10
         endif
      endif
c
c     Check that the orbitals energies are in aufbau
c     ordering. Otherwise the Fermi-Dirac smearing will mess
c     everything up as it insists to occuppy the orbitals 
c     according to aufbau.
c
c     Aufbau ordering is checked against the occupations.
c     Because interchanges of equally occupied orbitals
c     are insignificant and the orbital occupations are 
c     energy dependent this seems the best way to do it.
c
      do i=1,nmo-1
         if (f(i+1).gt.f(i)+tol) then
c           write(iwr,"(2(/1x,a3,a22))")
c    +                "mo","energy",
c    +                "--","------"
c           do j = 1, nmo
c              write(iwr,"(1x,i3,f22.8)")j,e(j)
c           enddo
            if (opg_root()) then
              write(iwr,"(1x,'WARNING: fermi_smear: aufbau violation',
     +                    ' detected')")
              write(iwr,"(1x,'orbitals ',i5,' and ',i5,' in non-aufbau',
     +                    ' order')")i,i+1
            endif
            onotaufbau = .true.
c           call caserr2("fermi_smear: aufbau violation detected")
         endif
      enddo
      if ((out.or.j.ge.maxit.or.onotaufbau).and.opg_root()) then
         write(iwr,"(/1x,'fully occupied orbitals: ',i5)")nfocc
         write(iwr,"( 1x,'all   occupied orbitals: ',i5)")noccmx
         write(iwr,"(2(/1x,a3,2a22))")
     +             "mo","energy","occupation",
     +             "--","------","----------"
         do i = 1, nmo
            write(iwr,"(1x,i3,2f22.8)")i,e(i),f(i)
         enddo
         write(iwr,*)
      endif
      if (j.ge.maxit) then
         call caserr2("fermi_smear: too many iterations")
      endif
      end
      subroutine llvmo_defaults
c
c     Set the defaults for handling the state averaging in llvmo
c
      implicit none
INCLUDE(common/llvmo)
      oaverage = .true.
      dtolavg  = 1.0d-4
      end
      subroutine llvmo_average(flag)
c
c     Turn state averaging on or off
c
      implicit none
      character*(*) flag
INCLUDE(common/iofile)
INCLUDE(common/llvmo)
      if (flag.eq."on") then
         oaverage = .true.
      else if (flag.eq."off") then
         oaverage = .false.
      else
         write(iwr,*)"invalid option   : ",flag
         write(iwr,*)"valid options are: on/off"
         call caserr("invalid option to AVERAGE directive")
      endif
      end
      subroutine llvmo_tolerance(dtol)
c
c     Set the orbital degeneracy tolerance for the state averaging
c
      implicit none
      REAL dtol
INCLUDE(common/iofile)
INCLUDE(common/llvmo)
      if (dtol.ge.0.0d0) then
         dtolavg = dtol
      else
         write(iwr,*)"invalid value   : ",dtol
         write(iwr,*)"valid values are: all non-negative values"
         call caserr("invalid value with AVERAGE directive")
      endif
      end
**==llvmo.f
      subroutine llvmo(e,q,noc,no,n)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),e(*)
INCLUDE(common/llvmo)
c     data dzero/0.0d0/
      data done/1.0d0/
c
c     ----- check for low lying vacant orbitals -----
c
      call vfill(done,q,1,noc)
      if (noc.lt.n) call vclr(q(noc+1),1,n-noc)
      no = noc
      if (.not.oaverage) then
c        if a single determinant wavefunction was requested explicitly
c        then we are done now.
         return
      endif
      if (noc.eq.n) return
      qo2 = e(noc)
      qo3 = e(noc+1)
      if (dabs(qo3-qo2).lt.dtolavg) then
c
c -----  how many orbitals are nearly degenerate ??
c
         nod = 0
         nou = 0
         do 20 i = 1 , noc
            if (dabs(e(i)-qo2).lt.dtolavg) nod = nod + 1
 20      continue
         nocp = noc + 1
         do 30 i = nocp , n
            if (dabs(e(i)-qo2).lt.dtolavg) nou = nou + 1
 30      continue
         ntot = nou + nod
         fac = dfloat(nod)/dfloat(ntot)
         m = noc - nod
         do 40 i = 1 , ntot
            m = m + 1
            q(m) = fac
 40      continue
         no = m
      end if
      return
      end
**==mofile.f
      subroutine mofile(q,zscf,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/runlab)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/prints)
INCLUDE(common/segm)
      dimension q(*)
      data zuhf /'uhf'    /
c     data dzero/0.0d0/
      data zanam,zbname/' -a-',' -b-'/
      data m1,m2/1,2/
      write (iwr,6010)
      out = nprint.eq.2
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c
c     ----- set pointers for partitioning of core -----
c
      i10 = igmem_alloc(l2+l3+3*l1)
      i20 = i10 + l2
      i30 = i20 + l3
      i40 = i30 + l1
      i50 = i40 + l1
      last = i50 + l1
c
CMR   debugging
CMR   write(iwr,*) 'CMR, mofile secsum says :'
CMR   call secsum()
CMR   debugging

      if (out) write (iwr,6040) i10 , i20 , i30 , i40 , last
c
      pop = 2.0d0
      if (zscf.eq.zuhf) pop = 1.0d0
      if (oprint(45)) write (iwr,6020)
_IFN1(civ)      call vfill(pop,q(i40),1,na)
_IF1(c)      call setsto(na,pop,q(i40))
_IF1(iv)      call setstr(na,pop,q(i40))
      if (na.lt.l1) call vclr(q(i40+na),1,l1-na)
      call getq(q(i20),q(i30),q(i50),i,j,m2,ieig,ipop,mina,zanam)
      if (nswapa.gt.0) call swap(q(i20),q(i30),l1,nswapa,nswap(1))
      call putq(zcom,ztitle,q(i30),q(i40),l1,l1,l1,m1,m1,q(i20),mouta,
     +          ibl3qa)
      call tdown(q(i20),ilifq,q(i20),ilifq,l1)
      if (oprint(45)) call prev(q(i20),q(i30),l1,l1,l1)
      call dmtx(q(i10),q(i20),q(i40),iky,na,l1,l1)
      call wrt3(q(i10),l2,ibl3pa,idaf)
c
      call wrt3(q(i30),l1,ibl3ea,idaf)
      if (zscf.eq.zuhf) then
c
c     ----- uhf case -- b orbitals
c
         if (nb.eq.0) then
            call setz(q(i20),q(i10),q(i30),l1,l2,l3)
         else
            if (oprint(45)) write (iwr,6030)
            call getq(q(i20),q(i30),q(i50),i,j,m2,ieig,ipop,minb,zbname)
            if (nswapb.gt.0) call swap(q(i20),q(i30),l1,nswapb,nswap(41)
     +                                 )
_IFN1(civ)            call vfill(pop,q(i40),1,nb)
_IF1(c)            call setsto(nb,pop,q(i40))
_IF1(iv)            call setstr(nb,pop,q(i40))
            if (nb.lt.l1) call vclr(q(i40+nb),1,l1-nb)
         end if
c
         call putq(zcom,ztitle,q(i30),q(i40),l1,l1,l1,m1,m1,q(i20),
     +             moutb,ibl3qb)
         call tdown(q(i20),ilifq,q(i20),ilifq,l1)
         if (oprint(45)) call prev(q(i20),q(i30),l1,l1,l1)
c
         call dmtx(q(i10),q(i20),q(i40),iky,nb,l1,l1)
         call wrt3(q(i10),l2,ibl3pb,idaf)
         call wrt3(q(i30),l1,ibl3eb,idaf)
      end if
c
c     ----- reset core memory -----
c
      call gmem_free_inf(i10,"guess.m","mofile",'i10')
c
      return
 6010 format (/,' molecular orbitals restored from dumpfile'/)
 6020 format (//30x,14('-'),/,30x,'alpha orbitals',/,30x,14('-')//)
 6030 format (//30x,14('-'),/,30x,'beta  orbitals',/,30x,14('-')//)
 6040 format (' core assignement '/,' i10, i20, i30, i40 =',4i8,/,
     +        ' last = ',i8)
      end
**==mogues.f
      subroutine mogues(q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/scra7)
INCLUDE(common/tran)
INCLUDE(common/statis)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/atmol3)
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/restrj)
INCLUDE(common/restri)
      character *8 title,guess
      common/restrz/title(12),guess
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/cslosc)
      common/blkcore/corev(512),charge(10)
INCLUDE(common/zorac)
INCLUDE(common/drfopt)
_IF(drf)
INCLUDE(../drf/comdrf/dafil)
INCLUDE(../drf/comdrf/drfpar)
      logical odrf
_ENDIF
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
      logical mpassi
      common/integmp/mpassi
c
      logical odft
      character*7 fnm
      character*6 snm
      data fnm,snm/'guess.m','mogues'/
c
      dimension q(*),o1e(6)
      data m19/19/
      data ztyp2,ztyp3,ztyp4,ztyp5,ztyp6 /'hcore',
     + 'minguess','extguess','cards' ,'mosaved'/
      data ztyp7,ztyp8,ztyp9,ztyp12,ztyp13/
     +'roread' ,'gvbread','lmoread',
     + 'nogen' ,'alphas'/
      data ztyp10,ztyp11/'getq','extra'/
      data ztyp16/'atoms'/
      data zuhf,zgrhf,zgvb/   'uhf'    ,'grhf'   ,'gvb'     /
c
c     ----- are 1e-routines to be called
c
      data oalpha/.true./
c
      data zap/'adapt'/
_IF(drf)
c
      odrf = field .ne. ' '
_ENDIF
_IF(ccpdft)
      odft = CD_active().and.oatmdft
_ELSE
      odft = .false.
_ENDIF
      call cpuwal(begin,ebegin)
      nav = lenwrd()
      if (.not.(opass1)) then
c
c ----- first restore s,t,f from section 192 and store on ed7
c
         if (itask(mtask).eq.-1) then
            call standv(1,q)
         end if
      end if
      l2 = num*(num+1)/2
c     l3 = num*num

_IF(ga)
c allocate  
      oscftp = zscftp.eq.zuhf
      call declare_diis_storage(num,oscftp)
_ENDIF

      i10=0
      if (mpassi) then
c     minimize overall memory usage
       i20 = i10 
       max1 = i20 + l2 + l2
      else
       i20 = i10 + l2
       i30 = i20 + l2
c ----- if adaption requested .. compute storage requirements
       max1 = i30 + l2
      endif
      i21 = i20 + mxgaus
      i22 = i21 + mxgaus
      i23 = i22 + mxprms*num
      i24 = i23 + mxprms*num
      i25 = i24 + (7*maxat+1)/nav
      i26 = i25 + (3*num+1)/nav
      i27 = i26 + (3*num+1)/nav
      i28 = i27 + (3*num+1)/nav
      i29 = i28 + (8*num+1)/nav
      i31 = i29 + (8*num+1)/nav
      max2 = i31 + (num+1)/nav
      max1 = max(max1,max2)
      last = max1
      length = last - i10

      i10 = igmem_alloc_inf(length,fnm,snm,'i10',IGMEM_DEBUG)
      if (mpassi) then
       i20 = i10
       max1 = i20 + l2 + l2
      else
       i20 = i10 + l2
       i30 = i20 + l2
c
c ----- if adaption requested .. compute storage requirements
c
       max1 = i30 + l2
      endif
c
c...   everybody gets a dummy symmetry adaption at least
c...   adapt will remember this
c
c
      i21 = i20 + mxgaus
      i22 = i21 + mxgaus
      i23 = i22 + mxprms*num
      i24 = i23 + mxprms*num
      i25 = i24 + (7*maxat+1)/nav
      i26 = i25 + (3*num+1)/nav
      i27 = i26 + (3*num+1)/nav
      i28 = i27 + (3*num+1)/nav
      i29 = i28 + (8*num+1)/nav
      i31 = i29 + (8*num+1)/nav
      max2 = i31 + (num+1)/nav
      max1 = max(max1,max2)
c
      if(mpassi) then
c
c      do this in 3 passes to reduce overall memory usage
c
       do mpass= 1,3
        do loop = 1,6
          o1e(loop) = .false.
        enddo
        o1e(mpass) = .true.
        call getmat(q(i10),q(i10),q(i10),q(i10),q(i10),q(i10),
     *              charge,num,o1e,ionsec)
        if (mpass.eq.1) then
         ibl7s = ibl7la
         call wrt3(q(i10),l2,ibl7s,num8)
        else if (mpass.eq.2) then
         ibl7t = iposun(num8)
         call wrt3(q(i10),l2,ibl7t,num8)
        else
         ibl7f = iposun(num8)
         call wrt3(q(i10),l2,ibl7f,num8)
         ibl7la = iposun(num8)
        endif
       enddo
      else
       do loop = 1,3
        o1e(loop) = .true.
        o1e(loop+3) = .false.
       enddo
       call getmat(q(i10),q(i20),q(i30),q(i10),q(i10),q(i10),
     +             charge,num,o1e,ionsec)

       ibl7s = ibl7la
       call wrt3(q(i10),l2,ibl7s,num8)
       ibl7t = iposun(num8)
       call wrt3(q(i20),l2,ibl7t,num8)
       ibl7f = iposun(num8)
       call wrt3(q(i30),l2,ibl7f,num8)
       ibl7la = iposun(num8)
      endif
c
      write (iwr,6010)
      call adapt(q(i20),q(i21),q(i22),q(i23),q(i24),q(i25),q(i26),q(i27)
     +           ,q(i28),q(i29),q(i31),num)
c
      if(mpassi) then
c     restore s matrix and i20 definition
       call rdedx(q(i10),l2,ibl7s,num8)
       i20 = i10 + l2
      endif
c
c ----- now transform 1e-integrals to sabf and output to ed7
c
      if (zguess.eq.ztyp12.or.(odscf.and.odnew)) then
         call anorm(q(i10),q)
         call rdedx(q(i10),l2,ibl7s,num8)
      end if
      call tranp(q(i10),q(i20))
      ibl7st = ibl7la
      call wrt3(q(i20),l2,ibl7st,num8)
_IF(ga)
c
c  check/create GA storage, load S, ST into GA storage
c
      call init_diis_storage(num,q(i10),q(i20))
_ENDIF
      ibl7tt = iposun(num8)
      call rdedx(q(i10),l2,ibl7t,num8)
      call tranp(q(i10),q(i20))
      call wrt3(q(i20),l2,ibl7tt,num8)
      ibl7ft = iposun(num8)
      call rdedx(q(i10),l2,ibl7f,num8)
      call tranp(q(i10),q(i20))
      call wrt3(q(i20),l2,ibl7ft,num8)
c
c ----- now decide on blocks for density matrix storage
c ----- and eigen value storage on dumpfile .. create sufficient
c ----- space for both -a- and -b- set vectors based on scftype
c ----- dumped to section 497 of dumpfile
c
      ibl7la = iposun(num8)
      length = lensec(l2)
      lene = lensec(num)
      len2 = lene + length
      oscftp = zscftp.eq.zuhf .or. zscftp.eq.zgrhf .or.
     +         zscftp.eq.zgvb 
      if(oscftp) len2 = len2 + len2
      if (ozora.and.oso) len2 = len2 + lene + length
      call secput(isect(497),m19,len2,ibl3pa)
      ibl3ea = ibl3pa + length
      if(oscftp) then
       ibl3pb = ibl3ea + lensec(num)
       ibl3eb = ibl3pb + length
      else
       ibl3pb = 0
      endif
c
c ----- reset core
c
      call gmem_free_inf(i10,fnm,snm,'i10')
c
c ----- now act on vectors specification
c
      call revind
_IF(drf)
      if (odrf) call hcore(q,zscftp,oalpha,zconf,nprint)
_ENDIF
c
c...  if an atomic zora is requested then calculate the zora
c...  corrections here, if not present --- so it is usable for all
c
      if (nat_z.gt.0.and.ozora.and.irest_z.eq.0) call denat(q,'zora')
c
c...   if vb, do not try to read anything here (will not understand)
      if (zruntp.ne.zap.and.zscftp.ne.'vb') then
         if (irest.gt.0 .and. zguess.ne.ztyp5 .and.
     +       zguess.ne.ztyp10 .and. zguess.ne.ztyp11 .and.
     +       zguess.ne.ztyp12) zguess = ztyp6
c
c
c     ----- start with bare nucleus hamiltonian or alphas matrix -----
c
         if (zguess.eq.ztyp2 .or. zguess.eq.ztyp13) then
            if (zguess.eq.ztyp13) oalpha = .false.
            call hcore(q,zscftp,oalpha,zconf,nprint)
c
c     ------  start with atomic densities  ------
c
         else if (zguess.eq.ztyp16) then
            if (nat_z.gt.0.and.odft.and..not.ozora) then
               call denat(q,'zora')
            else
               call denat(q,'start')
            endif
c
c     ----- minimal basis set initial guess routine -----
c
         else if (zguess.eq.ztyp3) then
            call gesmin(q,zscftp,zconf,nprint)
c
c     ----- extented basis set initial guess routine -----
c
         else if (zguess.eq.ztyp4) then
            call gesext(q,zscftp,zconf,nprint)
c
c     ----- read mo*s from cards -----
c
         else if (zguess.eq.ztyp5 .or. zguess.eq.ztyp7 .or.
     +            zguess.eq.ztyp8 .or. zguess.eq.ztyp9) then
            call readmo(q,zscftp,nprint)
c-ps
c     ----- read mo*s from cards (minimal memory version) -----
c
         else if (zguess.eq.'mcards') then
c
            call readmomm(q,zscftp,nprint)
c
c     ----- read saved mo*s -----
c
         else if (zguess.eq.ztyp6) then
            if (guess.ne.'atoms')
     +      call mofile(q,zscftp,nprint)
c
c     ----- generate additional natural orbitals for gvb pairs -----
c
         else if (zguess.eq.ztyp12) then
            call nogen(q)
c
c     ----- read mos from foreign dumpfile
c
         else if (zguess.eq.ztyp10) then
            call getqq(q,zscftp,nprint)
c
c     ----- read mos from foreign dumpfile
c           with basis set extension
c
         else if (zguess.eq.ztyp11) then
            call extra(q,zscftp,nprint)
c
c
c     ----- error in initial vectors option -----
c
         else
            write (iwr,6030) zguess
            call caserr2('unrecognised vectors option')
         end if
      end if
c
_IF(drf)
      if (odrf) then
      l3 = num*num
      i10 = igmem_alloc_inf(l3,fnm,snm,'i10',IGMEM_NORMAL)
c
      call rdedx(q(i10),l3,ibl3qs,idaf)
      call dawrit(idafh,ioda,q(i10),l3,13,navh)
      call rdedx(q(i10),l3,ibl3qa,idaf)
      call dawrit(idafh,ioda,q(i10),l3,15,navh)
c
c     ----- now density matrix ..
c
      call rdedx(q(i10),l2,ibl3pa,idaf)
      if (igetden.eq.1) then
        call daread(idafh,ioda,q(i10),l2,16)
      call wrt3(q(i10),l2,ibl3qs,idaf)
      else
        call dawrit(idafh,ioda,q(i10),l2,16,navh)
      endif
cahv  call wrt3(q(i40),l1,ibl3ea,idaf)
cahv  ndaf = 4
      if (zscftp.eq.zuhf) then
        call rdedx(q(i10),l3,ibl3qb,idaf)
        call dawrit(idafh,ioda,q(i10),l3,19,navh)
cahv     if (nb.eq.0) then
c           nsav = moutb
c           call setz(q(i30),q(i10),q(i40),l1,l2,l3)
c           call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m1,m1,q(i30),
c    +                nsav,iblkq(ndaf))
            call rdedx(q(i10),l2,ibl3pb,idaf)
            call dawrit(idafh,ioda,q(i10),l2,20,navh)
cahv        call wrt3(q(i40),l1,ibl3eb,idaf)
cahv     end if
      end if
c
      call gmem_free_inf(i10,fnm,snm,'i10')
      end if
_ENDIF
      write (iwr,6020)
      call timana(2)

      return
 6010 format (/40x,17('*')/40x,'vector generation'/40x,17('*')/)
 6020 format (/1x,104('-')/)
 6030 format (/1x,' error in initial vectors option, guess = ',a8)
      end
**==nise.f
      function nise(ise)
      integer ise(7) 
c identity is not represented in ise, so start count with n=1
      n=1
      do 1 i=1,7
 1      if (ise(i).eq.1) n=n+1 
      nise=n
      return
      end
**==nogen.f
      subroutine nogen(q)
c
c     ----- nogen will generate the the second natural orbital
c           of a gvb pair based on the first natural orbital.
c     the guess is made by placing a node between the
c     atom with the largest coefficients and the rest
c     of the mo-s, then orthogonalizing.  the rest of
c     the virtual space is filled up randomly -----
c
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/timez)
INCLUDE(common/runlab)
INCLUDE(common/segm)
      dimension q(*)
c
      data m0/0/
      l1 = num
      l2 = num*(num+1)/2
c
c     allocate space. the arrays used are
c     i10  v
c     i20  e
c     i30  ia
c     i40  s
c     i50  t
c     i60  moordr
c
c
      l3 = num*num
c
      i10 = igmem_alloc(l3 + l1 + l1 + l3 + l2 + l1 + 9*l1)
c
      i20 = i10 + l3
      i30 = i20 + l1
      i40 = i30 + l1
      i50 = i40 + l3
      i60 = i50 + l2
      i70 = i60 + l1
c
c     last = i70 + 9*l1
c
c     ----- read in all nogen input -----
c
c     ----- generate the natural orbitals -----
c
      call noin(nsplit,nmo,nswapa,q(i10),nswap,l1)
      call split(q(i10),nsplit,nmo,nat,l1)
c
c     ----- calculate the reduced density matrices -----
c
      call othogs(q(i10),q(i40),q(i50),q(i20),q(i60),iky,q(i70),nmo,
     +            nsplit,l1,l2,l3)
c
c     ----- store zero as eigenvalues for alpha and beta -----
c
      call redden(q(i10),q(i40),q(i50),q(i20),l1,l2)
c
c ----- transform vectors to sabf
c
      call vclr(q(i30),1,l1)
      call tback(q(i10),ilifq,q(i10),ilifq,l1)
      call putq(zcom,ztitle,q(i30),q(i30),l1,l1,l1,m0,m0,q(i10),mouta,
     +          ibl3qa)
      call wrt3(q(i30),l1,ibl3ea,idaf)
      call wrt3(q(i30),l1,ibl3eb,idaf)
c
      call gmem_free_inf(i10,"guess.m","nogen",'i10')
      return
      end
**==noin.f
      subroutine noin(numgen,nmo,norder,v,moordr,l1)
c
c     ----- read in all nogen input parameters -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/scfwfn)
      dimension v(l1,*)
      dimension moordr(*)
c
c     ----- generate numgen,nmo
c     numgen = number of new no-s to be generated
c     nmo   = number of mo-s to be read in
c     norder = array of new mo order -----
c
      data m2,zanam/2,' -a- '/
      numgen = npair
      j = nco
      if (nseto.ne.0) then
         do 20 i = 1 , nseto
            j = j + no(i)
 20      continue
      end if
c
      nmo = j + numgen
      if (nmo.eq.0) then
         if (numgen.ne.0 .or. norder.ne.0) then
            call caserr2('parameters of nogen directive are invalid')
         end if
      end if
 30   if ((numgen+nmo).gt.l1) then
         call caserr2('parameters of nogen directive are invalid')
         go to 30
      else
c
c        noc = nmo + numgen
c
         write (iwr,6010) nmo , numgen
         call getq(v,alphas,alphas,i,j,m2,ieig,ipop,mina,zanam)
         if (norder.gt.0) call swap(v,alphas,l1,norder,moordr)
c
c     ----- zero out the virtual space if any -----
c
         call tdown(v,ilifq,v,ilifq,nmo)
         nmo1 = nmo + 1
         i = (l1-nmo)*l1
         if (i.gt.0) call vclr(v(1,nmo1),1,i)
         return
      end if
 6010 format (/10x,26('*'),/,10x,'natural orbital generation',/,10x,
     +        26('*'),//,10x,'nmo',7x,i5,/,10x,'numgen',4x,i5,/)
      end
**==oeigd.f
      subroutine oeigd(fc, s, u, t, h)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     one-electron integrals. general.
c.......................................................................
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension fc(*),s(*),u(*),t(*),h(*)
      ufacl = dsqrt(8.0d0/3.1415926536d0)
      nstep = 0
      k = 0
      do 40 l = 1 , nsym
         nbas1 = nbas(l)
         expfac = dfloat(l) + 0.5d0
         do 30 i = 1 , nbas1
            do 20 j = 1 , i
               k = k + 1
               zp = zeta(nstep+i)
               zq = zeta(nstep+j)
               zpq = 0.5d0*(zp+zq)
               term1 = dsqrt(zpq)
               ppq = zp*zq/zpq
               rpq = dsqrt(ppq/zpq)
               s(k) = rpq**expfac
               fc(k) = s(k)
               u(k) = ufacl*s(k)*term1
               t(k) = expfac*s(k)*ppq
               h(k) = t(k) - zn*u(k)
 20         continue
 30      continue
         nstep = nstep + nbas1
         ufacl = ufacl*2.0d0*dfloat(l)/(2.0d0*dfloat(l)+1.0d0)
 40   continue
      return
      end
**==orderd.f
      subroutine orderd(amat,nrow,ncol,imemb,nblock,ind,vec,iflag)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c     this routine sorts a set of column vectors in amat(*,*) according
c     to increasing values in vec(*). care is taken that the output orde
c     of the vectors is as close to the input order as possible.
c.......................................................................
      dimension amat(nrow,ncol), ind(*), vec(*), imemb(*)
c.......................................................................
c
c     set tolerance for test, then determine for each element how
c     many smaller elements vec(*) contains.
c.......................................................................
      tol = 1.0d-10
      do 30 i = 1 , ncol
         test = vec(i) - tol
         indi = 1
         do 20 j = 1 , ncol
            if (vec(j).lt.test) indi = indi + 1
 20      continue
         ind(i) = indi
 30   continue
c.......................................................................
c
c     if desired,scan ind(*) to determine the number of different values
c     in vec(*) and the number of elements of each value.
c.......................................................................
      if (iflag.le.99) then
         do 40 i = 1 , ncol
            imemb(i) = 0
 40      continue
         do 50 i = 1 , ncol
            indi = ind(i)
            imemb(indi) = imemb(indi) + 1
 50      continue
         icount = 0
         do 60 i = 1 , ncol
            if (imemb(i).ne.0) then
               icount = icount + 1
               imemb(icount) = imemb(i)
            end if
 60      continue
         nblock = icount
      end if
c.......................................................................
c
c     establish order in degeneracies of the ordering vector.
c.......................................................................
      do 80 i = 2 , ncol
         itest = ind(i)
         im1 = i - 1
         do 70 j = 1 , im1
            if (ind(j).eq.itest) itest = itest + 1
 70      continue
         ind(i) = itest
 80   continue
c.......................................................................
c
c     ind(*) contains ordering indices for amat(*,*). sort following
c     input order as far as this is correct.
c.......................................................................
      ilow = 1
 90   do 100 i = 1 , ncol
         if (i.ne.ind(i)) go to 110
 100  continue
      go to 130
c.......................................................................
c
c     input order wrong. swap present vector into correct col.
c.......................................................................
 110  index = ind(i)
      do 120 j = 1 , nrow
         scra = amat(j,i)
         amat(j,i) = amat(j,index)
         amat(j,index) = scra
 120  continue
      scra = vec(i)
      vec(i) = vec(index)
      vec(index) = scra
      ind(i) = ind(index)
      ind(index) = index
      if (ind(i).ne.i) go to 110
c.......................................................................
c
c     amat(*,*) is ordered through i. go back to checking order.
c.......................................................................
      ilow = i + 1
      if (ilow.lt.ncol) go to 90
 130  return
      end
**==osatod.f
      function osatod(iat,ic,iiloc,iisch,nbb)
c
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
      logical osatod
c
c...   function checks if this atom is the same as one for which info
c...   is in /junk/ and from which we can use the density matrix
c...   for atomic startup
c...   check also on qtomic name
c...   if so osatod is true  and iiloc is updated
c...
c
INCLUDE(common/sizes)
c...   commons where information comes from :
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/runlab)
      dimension cspd(mxprim,4)
      equivalence (cspd(1,1),cs(1))
c
c...   common where info goes to :
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension ic(6,nbb),iiloc(nbb,6),iisch(nbb,6)
      dimension nbcn(5)
c..
      data pi32/5.56832799683170d0/
c     data dzero/0.0d0/
      data pt5,pt75/0.5d0,0.75d0/
      data pt187 /1.875d+00/
c..
      osatod = .false.
      if (iatom.le.0) return
      if (czan(iat).ne.zn) return
      if (zaname(iat).ne.zaname(iatom)) return
c
c...  gather  shell / symmetry info
c
      do 20 i = 1 , 4
         nbcn(i) = 0
 20   continue
c
c.. nbc  # shell-s / symmetry
c.. iisch  contains index of shell
c.. iiloc  contains position of starting ao of shell in "real" world
c
      do 40 ii = 1 , nshell
         i = katom(ii)
         if (i.eq.iat) then
            mini = kmin(ii)
            maxi = kmax(ii)
            kk = ktype(ii)
            do 30 iorb = mini , maxi
               if (iorb.eq.1) then
c..  translate to 1 (s)
                  nbcn(1) = nbcn(1) + 1
                  iisch(nbcn(1),1) = ii
                  iiloc(nbcn(1),1) = kloc(ii)
               else if (iorb.eq.2 .or. iorb.eq.5 .or. iorb.eq.11) then
c..  translate to 2 (p) 3(d) or  4(f)
                  ispdf = kk
                  nbcn(ispdf) = nbcn(ispdf) + 1
                  iisch(nbcn(ispdf),ispdf) = ii
                  iiloc(nbcn(ispdf),ispdf) = kloc(ii) + iorb - mini
               end if
 30         continue
         end if
 40   continue
c..     check nbcn
      do 50 i = 1 , 4
         if (nbc(i).ne.nbcn(i)) return
 50   continue
c..
c..     we gathered symmetry/shell info ; now check the real thing
c..
      kkzc = 0
c     kh = 0
c     isymax = 0
      do 80 ispdf = 1 , 4
c..      nbas = total # primitives for this symmetry
         nbasn = 0
c        if (nbc(ispdf).gt.0) isymax = ispdf
         do 70 j = 1 , nbc(ispdf)
            ii = iisch(j,ispdf)
            is = kstart(ii)
            if = is + kng(ii) - 1
c..      ic = # number of primitives /contracted /symmetry
            if (ic(ispdf,j).ne.kng(ii)) return
            nbasn = nbasn + kng(ii)
c..      check the primitives / watch the subtle use of 2-dim cspd
            do 60 k = is , if
               kkzc = kkzc + 1
               if (zeta(kkzc).ne.ex(k)) return
               contn = cspd(k,ispdf)
c...     get contraction coeff-s as we are used to
               ee = 2*zeta(kkzc)
               fac = pi32/(ee*dsqrt(ee))
               if (ispdf.eq.2) then
                  fac = pt5*fac/ee
               else if (ispdf.eq.3) then
                  fac = pt75*fac/(ee*ee)
               else if (ispdf.eq.4) then
                  fac = pt187*fac/(ee**3)
               end if
               contn = contn*dsqrt(fac)
               if (cont(kkzc).ne.contn) return
 60         continue
 70      continue
c...
         if (nbasn.ne.nbas(ispdf)) return
c...
 80   continue
c..
c..      all checks out
c..
      osatod = .true.
c..
      return
      end
**==othogs.f
      subroutine othogs(v,s,t,e,moordr,ia,b,m,nsplit,l1,l2,l3)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/prints)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
      dimension v(l1,*),s(l1,*),t(l2),e(l1),ia(l1),moordr(*)
      dimension b(*)
      data dzero,done /0.0d0,1.0d0/
      noc = m + nsplit
c
c     ----- make sure that the input orbitals are orthonormal -----
c
      call rdedx(t,l2,ibl7s,num8)
      call square(s,t,l1,l1)
      call schmds(v,s,t,m,l1)
c
c     ----- schmidt each split mo against the m old mo-s
c
c     first normalize all of the noc vectors
      if (nsplit.ne.0) then
         do 30 i = 1 , noc
            dum = dzero
            do 20 j = 1 , l1
               if (v(j,i).ne.dzero) dum = dum +
     +             ddot(l1,s(1,j),1,v(1,i),1)*v(j,i)
 20         continue
            dum = done/dsqrt(dum)
            call dscal(l1,dum,v(1,i),1)
 30      continue
         do 60 i = 1 , nsplit
            do 50 j = 1 , m
               dum = dzero
               do 40 k = 1 , l1
                  fact = v(k,i+m)
                  if (fact.ne.dzero) dum = dum +
     +                ddot(l1,s(1,k),1,v(1,j),1)*fact
 40            continue
               call daxpy(l1,-dum,v(1,j),1,v(1,i+m),1)
 50         continue
 60      continue
         do 80 i = 1 , nsplit
            dum = dzero
            do 70 k = 1 , l1
               fact = v(k,i+m)
               if (fact.ne.dzero) dum = dum +
     *          ddot(l1,s(1,k),1,v(1,i+m),1) * fact
 70         continue
            dum = dsqrt(done/dum)
c
c     ----- save the vectors in v - reuse v for scratch -----
c
            call dscal(l1,dum,v(1,i+m),1)
 80      continue
c
c     ----- calculate s over mo-s -----
c
         call wrt3(v,l3,ibl7la,num8)
c
c     ----- get s**(-1/2) -----
c
         call ttran(t,s,v,e,ia,m,nsplit,l1)
         m2 = nsplit*(nsplit+1)/2
c
c     ----- get the orthogonal v -----
c
         call symsqs(t,v,e,s,b,ia,nsplit,m2,nsplit)
         call rdedx(v,l3,ibl7la,num8)
         call vorth(v,s,e,m,nsplit,l1)
c
c     ----- reorder the no-s so that they conform to standard
c           format -----
c
         if (nsplit.ne.0) then
            ihi = m - nsplit
            if (ihi.ge.1) then
               do 90 i = 1 , ihi
                  moordr(i) = i
 90            continue
            end if
            ilo = m - nsplit + 1
            ii = ilo
            do 100 i = ilo , m
               moordr(ii) = i
               moordr(ii+1) = i + nsplit
               ii = ii + 2
 100        continue
         end if
         call reordr(v,moordr,noc,l1)
      end if
      if (oprint(45)) write (iwr,6010)
      if (oprint(45)) call prsql(v,m+nsplit,l1,l1)
      return
 6010 format (/10x,22('=')/10x,'gvb guess eigenvectors'/10x,22('='))
      end
**==outpud.f
      subroutine outpud(copn,cc,ntest,iwr)
      implicit REAL  (a-h,o-z),integer   (i-n)
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension cc(*),copn(*)
      if (ntest.gt.0) write (iwr,6010) nitscf , energ , cin , vir
      if (ntest.gt.0) write (iwr,6030)
      nbc1 = nbc(1)
      do 20 i = 2 , nsym
         nbc1 = max(nbc1,nbc(i))
 20   continue
      naddr = 0
      noddr = 0
      do 60 i = 1 , nsym
         do 50 j = 1 , ncsh(i) + nosh(i)
            do 30 k = 1 , nbc(i)
               naddr = naddr + 1
               noddr = noddr + 1
               copn(noddr) = cc(naddr)
 30         continue
            if (nbc(i).lt.nbc1) then
               do 40 k = nbc(i) + 1 , nbc1
                  noddr = noddr + 1
                  copn(noddr) = 0.0d0
 40            continue
            end if
 50      continue
         naddr = naddr + nbc(i)*(nbc(i)-ncsh(i)-nosh(i))
 60   continue
      jfirst = 1
 70   jlast = min(jfirst+7,nsht)
      if (ntest.gt.0) then
         write (iwr,6020) (eps(j),j=jfirst,jlast)
         write (iwr,6030)
         do 80 i = 1 , nbc1
            write (iwr,6040) (copn(i+(j-1)*nbc1),j=jfirst,jlast)
 80      continue
      end if
      if (jlast.eq.nsht) return
      jfirst = jfirst + 8
      if (ntest.gt.0) write (iwr,6030)
      go to 70
 6010 format (/,8x,'final scf results at iteration',i4,/,8x,
     +        'total hf energy',4x,'kinetic energy',4x,'virial theorem',
     +        /4x,3(e19.10),//,8x,'orbital energies and eigenvectors')
 6020 format (/,2x,8(2x,f12.5),1x)
 6030 format (' ')
 6040 format (2x,8f14.6)
      end
**==pdfded.f
      subroutine pdfded(kdim,d,dhelp,factor,dmult,lm,nbci,iiloc)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     routine for distributing atomic densities to the molecular
c     density matrix. note that the number of orbitals differ in
c     the molecular and atomic case for d and f orbitals. transformat-
c     ion matrices are provided in dmult(*,*). to work out these tables
c     real atomic orbitals are needed. for f functions these are:
c              sqrt(1/60)*(2zzz-3zyy-3zxx)
c              sqrt(1/40)*(4zzy-yyy-xxy)
c              sqrt(1/40)*(4zzx-xxx-xyy)
c                    xyz
c              sqrt(1/4)*(xxz-yyz)
c              sqrt(1/24)*(3xxy-yyy)
c              sqrt(1/24)*(3xyy-xxx)
c     normalization of primitives is given by (xyz:xyz)=1, (xxy:xxy)=3
c     (xxx:xxx)=15.
c
c.......................................................................
c..
      dimension d(*),dhelp(*),dmult(kdim,kdim),iiloc(nbci)
      itria2(i,j,na) = (max(i,j)+na-1)*(max(i,j)+na-2)/2 + min(i,j)
c..
c..
      do 50 l = 1 , nbci
         lmsave = lm
         do 40 na = 1 , kdim
            lm = lmsave
            do 30 m = 1 , l
               noff = itria2(iiloc(m),iiloc(l),na) - 1
               lm = lm + 1
               delem = dhelp(lm)*factor
               nbrang = kdim
               if (m.eq.l) nbrang = na
               do 20 nb = 1 , nbrang
                  noff = noff + 1
                  d(noff) = delem*dmult(na,nb)
 20            continue
 30         continue
 40      continue
 50   continue
c..
      return
      end
**==projec.f
      subroutine projec(a,b,res,aa,bb,mdim,nbasis,mbasis,nonvec,q)
c
c    *******************************************************************
c
c          this subroutine projects the mo coefficients in the array
c     a into the space of the current basis set.  the array b is
c     also used, but the projected mo coefs are returned in a.  the
c     routine also handles all necessary normalization and completion
c
c    *******************************************************************
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension a(mdim,mdim),b(mdim,mdim),aa(mdim),bb(mdim)
      dimension res(mdim,mdim)
      dimension q(*)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
      common/junk/ilifc2(maxorb),nterm2(maxorb),it2(mxorb3),
     * ctrn2(mxorb3),iftrn2(2),
     * ex2(6,mxprim),cspac(maxat),kstrt2(7,mxshel),kstr(4),nprini(700),
     *iovmat,iosvec,ioproj,iorthg,c2(3,maxat)
c
c     tread the projection matrix.
      l4 = mdim*mdim
      call rdedx(b,l4,ioproj,num8)
c
c     project the mo coefs.
      call mxmg(b,1,mdim,a,1,mdim,res,1,mdim,nbasis,mbasis,mbasis)
c
c     now, re-orthonormalize the projected molecular orbitals.
c
c     first, calculate the matrix q=c(transpose).c.    (remember that
c     these are the transformed coefficients.)
      call wrt3(res,l4,iorthg,num8)
c
      call dagger(nbasis,mbasis,res,mdim,b,mdim)
      call mxmg(b,1,mdim,res,1,mdim,a,1,mdim,mbasis,nbasis,mbasis)
c     now, take the inverse square root of the matrix q.
c...   assume that we will find at most the number of zero-s in the input
c...   set (nonvec) or due to harmonic (newbas1-newbas0)
      call rootmt(a,b,res,bb,aa,mdim,mbasis,1,
     1            max(nonvec,newbas1-newbas0))
c
c     now, right multiply the mo coefficients by sqrt(q).  the
c     resulting coefficients will be orthonormal.
      call rdedx(b,l4,iorthg,num8)
      call mxmg(b,1,mdim,a,1,mdim,res,1,mdim,nbasis,mbasis,mbasis)
      if (mbasis.lt.nbasis) then
         call rdedx(b,l4,iosvec,num8)
c
c     finally, complete the coefficient matrix if necessary.
         call complt(res,b,aa,mdim,nbasis,mbasis,mbasis)
      end if
c...   add vectors if null vectors were present
      if (nonvec.ne.0) call compv(res,nbasis-nonvec,nbasis,q,'fill')
c
      call rdedx(b,l4,iovmat,num8)
c
c     back transform the mo coefs.
c
      call mxmg(b,1,mdim,res,1,mdim,a,1,mdim,nbasis,nbasis,nbasis)
c
      return
      end
**==readmo.f
      subroutine readmo(q,zscf,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/scfwfn)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/runlab)
INCLUDE(common/scra7)
INCLUDE(common/harmon)
      dimension q(*)
      dimension zcas(2)
_IF(ga)
INCLUDE(common/parcntl)
INCLUDE(common/gadiis)
      character*1 xn
      data xn/'n'/
_ENDIF
      data dzero,done,two /0.0d0,1.0d0,2.0d0/
      data m0,m1/0,1/
      data zuhf,zgrhf,zgvb/   'uhf'    ,'grhf'   ,'gvb'     /
      data zcas/'casscf','mcscf'/
c
      if(ovcfre)then
         write (iwr,6011)
      else
         write (iwr,6010)
      endif
      if (osemi) then
       irold = ird
       write(iwr,6070)
       call archiv ('vect',iwr)
       call rmcore
       nvread = num
      else
       nvread = na
      endif
      out = nprint.eq.2
      l0 = newbas0
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      ibl7qa = ibl7la
      ibl7qb = ibl7qa + lensec(l3)
c
c     ----- set pointers for partitioning of core -----
c
      i10 = igmem_alloc(l3 + l1 + l3 + l1 + l3 + l1 + l1)
c
      i11 = i10 + l3
      i20 = i11 + l1
      i21 = i20 + l3
      i30 = i21 + l1
      i31 = i30 + l3
      i40 = i31 + l1
      last= i40 + l1
c
      if (out) write (iwr,6060) i10 , i11 , i20 , i21 , i30 , i31 ,
     +                          i40, last
c
c     ----- get canonical orbitals -----
c
c     -s- at x(i20)
c     -q- at x(i30)
c
      call rdedx(q(i20),l2,ibl7st,num8)
      call qmat(q,q(i20),q(i30),q(i31),q(i10),iky,l0,l1,l3,l1,out)
      pop = two
      if (zscf.eq.zuhf) pop = done
      ipass = 1
      mblq = ibl7qa
      noc = na
      mswap = nswapa
      ibase = 0
      if (zscf.eq.zgrhf .or. zscf.eq.zgvb) then
c
c     ----- if grhf or gvb, get the number of mo-s -----
c
         noc = nco + npair + npair
         if (nseto.ne.0) then
            do 20 i = 1 , nseto
               noc = noc + no(i)
 20         continue
         end if
      end if
 30   continue
      if (isunor.ne.0) then
c
c     -s- at x(i10)
c     form s**(-1/2) for transforming mopac vectors
c
          call rdedx(q(i10),l2,ibl7s,num8)
          call square(q(i20),q(i10),l1,l1)
          call rootmt(q(i20),q(i10),q(i30),q(i21),q(i31),l1,l1,1,0)
c
      endif
c
      if(ovcfre)then
c
c free format
c
         do 35 j = 1,nvread
            read (ird,*) (q(i30-1+(j-1)*l1+i),i=1,num)
 35      continue
      else
         do 50 j = 1 , nvread
            imax = 0
            ic = 0
 40         imin = imax + 1
            imax = imax + 5
            ic = ic + 1
            if (imax.gt.num) imax = num
            read (ird,6040) jj , icc , (q(i30-1+(j-1)*l1+i),i=imin,imax)
            if (jj.ne.j .or. icc.ne.ic) then
               write (iwr,6050) j , ic
               call caserr2('error in vector card input')
            end if
            if (imax.lt.num) go to 40
 50      continue
      endif
c
c...  transform orbitals back to sabf basis
c
      isovl = igmem_alloc(l3)
      call  anorm(q(isovl),q)
      call gmem_free_inf(isov1,"guess.m","readmo",'isov1')
      call tback(q(i30),ilifq,q(i30),ilifq,l0)
c
      if(isunor.ne.0)then
c
c     -s- at x(i10)
c     -q- at x(i30)
c
c     transform mopac vectors with s**(-1/2)
c
         call mxmg(q(i20),1,l1,q(i30),1,l1,q(i10),1,l1,l1,l1,l1)
         call dcopy(l3,q(i10),1,q(i30),1)
c
c assign pops, and write to section isunor
c
         do 55 i = 1 , l0
            pop = dzero
            if (i.le.na) pop = done
            if (i.le.nb) pop = two
            q(i-1+i40) = pop
 55      continue
         isect = isunor + ipass - 1
         if (oprint(45)) write (iwr,6020)
         if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
         write (iwr,7000)isect
         call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m0,m1,q(i30),
     +        isunor,iblkps)
      endif
c     ----- orthonormalize the orbitals -----
c
c     -s- at x(i10)
c     -q- at x(i10)
c     -v- at x(i30)
c
      call ortho1(q(i10),q(i10),q(i30),q(i40),iky,noc,l0,l1,l2,l3,l1)
      call rdedx(q(i10),l3,ibl3qs,idaf)
c
c     ----- back-transform the mo-s -----
c
c     -q- at x(i10)
c     -v- at x(i30)
c
_IF(ga)
      if(l1 .ge. idpmult2 .and. ipiomode .ne. IO_NZ_S)then

         call load_ga_from_square(ih_scr,q(i30),l1)
         call load_ga_from_square(ih_vec,q(i10),l1)

         call pg_dgemm(xn,xn,l1,l0,l0,1.0d0,
     &     ih_vec, ih_scr, 0.0d0, ih_scr2)

         call load_square_from_ga(q(i30),ih_scr2, l1)

      else
          call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
      endif
_ELSE
      call tfsqc(q(i30),q(i10),q(i20),l0,l1,l1)
_ENDIF
c
c     align orbitals for mcscf / force     pjar
c
      if (zscf.eq.zcas(2).and.zruntp.eq.'force') then
         if (ibase.ne.0) call caserr2('readmo error with mcscf / force')
         write(iwr,*) 'mc swap call from readmo'
         call swap_mcf(q(i30),mswap,nswap)
      endif
      if (mswap.gt.0) call swap(q(i30),q(i40),num,mswap,nswap(ibase+1))
c
      call wrt3(q(i30),l3,mblq,num8)
      if (zscf.ne.zuhf) then
c
c     ----- rhf case -----
c
         if (oprint(45)) write (iwr,6020)
         do 60 i = 1 , l0
            pop = dzero
            if (i.le.na) pop = done
            if (i.le.nb) pop = two
            q(i-1+i40) = pop
 60      continue
         call rdedx(q(i30),l3,ibl7qa,num8)
         call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m0,m1,q(i30),
     +             mouta,ibl3qa)
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
         if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
         call dmtx(q(i10),q(i30),q(i40),iky,na,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pa,idaf)
         call wrt3(q(i40),l1,ibl3ea,idaf)
         if (na.ne.nb) then
            if (nb.ne.0) call dmtx(q(i20),q(i30),q(i40),iky,nb,l1,l1)
            if (nb.eq.0) call setz(q(i30),q(i10),q(i40),l1,l2,l3)
            call vsub(q(i10),1,q(i20),1,q(i10),1,l2)
            call wrt3(q(i10),l2,ibl3pb,idaf)
         end if
         go to 70
      else if (nb.eq.0) then
         call setz(q(i30),q(i10),q(i40),l1,l2,l3)
         call wrt3(q(i30),l3,ibl7qb,num8)
      else if (ipass.ne.2) then
         ipass = 2
         noc = nb
         mblq = ibl7qb
         if (osemi) then
          nvread = num
         else
          nvread = na
         endif
         mswap = nswapb
         ibase = 40
         go to 30
      end if
c
c     ----- uhf case -----
c
      if (oprint(45)) write (iwr,6020)
_IFN1(civ)      call vfill(done,q(i40),1,na)
_IF1(c)      call setsto(na,done,q(i40))
_IF1(iv)      call setstr(na,done,q(i40))
      if (na.lt.l0) call vclr(q(i40+na),1,l0-na)
      call rdedx(q(i30),l3,ibl7qa,num8)
      if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
      call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m0,m1,q(i30),mouta,
     +          ibl3qa)
      call tdown(q(i30),ilifq,q(i30),ilifq,l0)
      call dmtx(q(i10),q(i30),q(i40),iky,na,l1,l1)
      if (out) call prtril(q(i10),l1)
      call wrt3(q(i10),l2,ibl3pa,idaf)
c
      call wrt3(q(i40),l1,ibl3ea,idaf)
      if (nb.ne.0) then
         if (oprint(45)) write (iwr,6030)
_IFN1(civ)         call vfill(done,q(i40),1,nb)
_IF1(iv)         call setstr(nb,done,q(i40))
_IF1(c)         call setsto(nb,done,q(i40))
         if (nb.lt.l0) call vclr(q(i40+nb),1,l0-nb)
         call rdedx(q(i30),l3,ibl7qb,num8)
         if (oprint(45)) call prev(q(i30),q(i40),l0,l1,l1)
         call putq(zcom,ztitle,q(i40),q(i40),l1,l1,l0,m0,m1,q(i30),
     +             moutb,ibl3qb)
         call tdown(q(i30),ilifq,q(i30),ilifq,l0)
         call dmtx(q(i10),q(i30),q(i40),iky,nb,l1,l1)
         if (out) call prtril(q(i10),l1)
         call wrt3(q(i10),l2,ibl3pb,idaf)
c
         call wrt3(q(i40),l1,ibl3eb,idaf)
      end if
c
c     ----- reset core memory -----
c
 70   call gmem_free_inf(i10,"guess.m","readmo",'i10')
c     mina = mouta
c     if (moutb.ne.0) minb = moutb
c     if ( isunor.ne.0) then
c      mina = isunor
c      if (moutb.ne.0) minb = isunor + 1
c     endif
      if (osemi) ird = irold
      return
 6010 format (/' initial guess orbitals read from cards'/)
 6011 format (/
     &     ' initial guess orbitals read from cards (free format)'/)
 6020 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')/)
 6030 format (//30x,28('-')/30x,'initial guess beta  orbitals'/30x,
     +        28('-')/)
 6040 format (i2,i3,5e15.8)
 6050 format (' in reading the mo coefficients from cards',/,
     +        ' the occupied orbitals are out of order ',2i5)
 6060 format (' core assignement '/' i10, i11, i20, i21, i30,',
     +        ' i31, i40 = '/7i8/' last = ',i8)
 6070 format (' restore vectors from mopac archive file')
 7000 format (/' vectors have been written to section      ',i4/
     +         ' of dumpfile after s**(-1/2) transformation')
      end
c
c  read mos from cards - this version runs in less than 6 triangles
c  so can be used to start up the parallel SCF code
c  RHF only, no orthogonalisation
c
      subroutine readmomm(q,zscf,nprint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/prints)
INCLUDE(common/infoa)
INCLUDE(common/segm)
INCLUDE(common/scfwfn)
INCLUDE(common/atmol3)
INCLUDE(common/datgue)
INCLUDE(common/runlab)
INCLUDE(common/scra7)
      dimension q(*)
      data dzero,done,two /0.0d0,1.0d0,2.0d0/
      data m0,m1/0,1/
      data zuhf/   'uhf'  /
c     data zgrhf,zgvb/  'grhf'   ,'gvb'     /
c
      write (iwr,6011) 

 6011 format (/
     &     ' initial guess orbitals read from cards (free format)'/)

      nvread = na

      out = nprint.eq.2
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      ibl7qa = ibl7la
      ibl7qb = ibl7qa + lensec(l3)

      if(zscftp .ne. 'rhf')call caserr2('invalid scf type in readmomm')
      l0=l1

      pop = two
      if (zscf.eq.zuhf) pop = done
c     noc = na

      ivec = igmem_alloc(l3)
      iden = igmem_alloc(l2)
      ipop = igmem_alloc(l1)
c
c read in free format
c
      do 35 j = 1,nvread
      write(6,*)'reading..',j
         read (ird,*) (q(ivec-1+(j-1)*l1+i),i=1,num)
 35   continue

c
c write to ed7
c
      call wrt3(q(ivec),l3,ibl7qa,num8)
c
c  define populations
c
      if (oprint(45)) write (iwr,6020)


      do 60 i = 1 , l0
         pop = dzero
         if (i.le.na) pop = done
         if (i.le.nb) pop = two
         q(ipop + i-1) = pop
 60   continue

      call putq(zcom,ztitle,q(ipop),q(ipop),l1,l1,l0,m0,
     &     m1,q(ivec), mouta,ibl3qa)

      call tdown(q(ivec),ilifq,q(ivec),ilifq,l0)
      if (oprint(45)) call prev(q(ivec),q(ipop),l0,l1,l1)
      call dmtx(q(iden),q(ivec),q(ipop),iky,na,l1,l1)
c                -d-     -v-    -p-

      if (out) call prtril(q(iden),l1)
      call wrt3(q(iden),l2,ibl3pa,idaf)
      call wrt3(q(ipop),l1,ibl3ea,idaf)
c
 70   continue

c     call setc(loadcm)
c     mina = mouta
c     if (moutb.ne.0) minb = moutb
c     if ( isunor.ne.0) then
c      mina = isunor
c      if (moutb.ne.0) minb = isunor + 1
c     endif
c      if (osemi) ird = irold
c
c     reset core
c
      call gmem_free_inf(ipop,"guess.m","readmomm",'ipop')
      call gmem_free_inf(iden,"guess.m","readmomm",'iden')
      call gmem_free_inf(ivec,"guess.m","readmomm",'ivec')
c
      return
 6010 format (/' initial guess orbitals read from cards'/)
 6020 format (//30x,28('-')/30x,'initial guess alpha orbitals'/30x,
     +        28('-')/)
 6030 format (//30x,28('-')/30x,'initial guess beta  orbitals'/30x,
     +        28('-')/)
 6040 format (i2,i3,5e15.8)
 6050 format (' in reading the mo coefficients from cards',/,
     +        ' the occupied orbitals are out of order ',2i5)
 6060 format (' core assignement '/' i10, i11, i20, i21, i30,',
     +        ' i31, i40 = '/7i8/' last = ',i8)
 6070 format (' restore vectors from mopac archive file')
 7000 format (/' vectors have been written to section      ',i4/
     +         ' of dumpfile after s**(-1/2) transformation')
      end




**==reordr.f
      subroutine reordr(v,moordr,nmo,l1)
c
c     ----- reorder a set of molecular orbitals -----
c
c
c     v       = initial vector set
c     moorder = reordering instruction. the i-th new mo
c               is the one currently stored in the moordr(i)
c               position
c     nmo     = number of molecular orbitals to reorder
c     l1      = dimension v
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(l1,*),moordr(*)
      do 30 i = 1 , nmo
         j = moordr(i)
_IF1(uv)         do 40 k = 1,l1
_IF1(uv)         tem = v(k,j)
_IF1(uv)         v(k,j) = v(k,i)
_IF1(uv)         v(k,i) = tem
_IF1(uv)   40    continue
_IFN1(uv)         call dswap(l1,v(1,j),1,v(1,i),1)
         i1 = i + 1
         do 20 k = i1 , nmo
            if (moordr(k).eq.i) moordr(k) = j
 20      continue
 30   continue
      return
      end
**==rmcore.f
      subroutine rmcore
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/runlab)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
c
      dimension ztit(103),nelcor(103)
      data maxtyp/103/
      data nelcor/
     $         0, 0, 2, 2, 2, 2, 2, 2, 2, 2,
     $         10, 10, 10, 10, 10, 10, 10, 10,
     $         18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 28,
     $                     28, 28, 28, 28, 28, 28,
     $ 67 * 0 /
      data ztit/
     $         'h ', 'he',
     $         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne',
     $         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $         'k ', 'ca',
     $                     'sc', 'ti', 'v ', 'cr', 'mn',
     $                     'fe', 'co', 'ni', 'cu', 'zn',
     $                     'ga', 'ge', 'as', 'se', 'br', 'kr',
     $ 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     $ 'in','sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     $ 'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     $ 'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po',
     $ 'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm',
     $ 'bk','cf','es','fm','md','no','lw'   /
c
      diff = 0.0d0
      do 10 iat = 1,nat
      zinp = zaname(iat)
      ztagg = ztit(jsubst(zinp))
      ktype = locatc(ztit,maxtyp,ztagg)
      if (ktype.le.0) then
         call caserr2('unidentified element in rmcore')
      end if
      fnel = dfloat(nelcor(ktype))
      czan(iat) = czan(iat) - fnel
      diff = diff + fnel
 10   continue
c
c  adjust number of electrons to take account of inner shells
c
      idiff = nint(diff)
      if (mod(idiff,2).ne.0)
     +    call caserr2('number of electrons is not even  --  rmcore')
      ne = ne - idiff
      na = na - idiff/2
      nb = nb - idiff/2
      write (iwr,6000) ne , na , nb
      return
 6000 format (/5x,'total no. of electrons =',i5/
     +          5x,'no. of occupied orbitals (alpha) =',i5/
     +          5x,'no. of occupied orbitals (beta)  =',i5/)
      end
**==rootmt.f
      subroutine rootmt(a,b,res,aa,bb,mdim,nbas,inv,nignore)
c
c          this routine takes the matrix currently in a and returns
c     either its square root or its inverse square root in a.  the
c     eigenvectors of the matrix initially in a are returned in b.
c     inv=0 means the square root of a is to be formed, and inv=1 means
c     that the inverse square root is to be formed.
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(mdim,mdim),b(mdim,mdim),aa(mdim),bb(mdim)
      dimension res(mdim,mdim)
      data done/1.0d0/
c
_IF1(f)      call eigrs(mdim,nbas,a,aa,bb,b,ierr)
_IFN1(f)      call f02abf(a,mdim,nbas,aa,b,mdim,bb,ifail)
c
      ignore = 0
      do 20 i = 1 , nbas
         if (aa(i).le.1.0d-12) then
             aa(i) = 1.0d0
             ignore = ignore + 1
         end if
         term = dsqrt(aa(i))
         if (inv.ne.0) term = done/term
c
_IF1(civu)         call scaler(nbas,term,res(1,i),b(1,i))
_IFN1(civu)         call vsmul(b(1,i),1,term,res(1,i),1,nbas)
 20   continue
c
      if (nignore.lt.ignore)
     +  call caserr2('old negative eigen value in projection matrix.')
c
      call mxmg(res,1,mdim,b,mdim,1,a,1,mdim,nbas,nbas,nbas)
      return
      end
**==schmds.f
      subroutine schmds(v,s,t,m,l1)
c
c     ----- schmidt orthonormalization routine -----
c
c     v = input vector set
c
c     s = overlap matrix
c
c     ia = triangular address array
c
c     t = scratch array of length l1
c
c     ilo,ihi = the low and high values for the range of mo-s
c               to be othogonalized
c
c     jlo,jhi = the low and high values for the range of mo-s
c               that ilo,ihi are to be orthogonalized against
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(l1,*),s(l1,*),t(*)
      data dzero,done,two/0.0d0,1.0d0,2.0d0/
      do 50 i = 1 , m
         do 20 l = 1 , l1
            t(l) = ddot(l1,v(1,i),1,s(1,l),1)
 20      continue
         do 30 j = 1 , i - 1
            dum = -ddot(l1,t,1,v(1,j),1)
            call daxpy(l1,dum,v(1,j),1,v(1,i),1)
 30      continue
         dum = dzero
         do 40 k = 1 , l1
            fact = two*v(k,i)
            if (fact.ne.dzero) then
               dum = dum + ddot(k,s(1,k),1,v(1,i),1)*fact
               dum = dum - v(k,i)*s(k,k)*v(k,i)
            end if
 40      continue
         dum = done/dsqrt(dum)
         call dscal(l1,dum,v(1,i),1)
 50   continue
      return
      end
**==scrint.f
      subroutine scrint
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      common/bufb/ssss(465),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,tolsp,nroosp,ni,nj
INCLUDE(common/hermit)
INCLUDE(common/wermit)
      dimension mmin(6),mmax(6)
      data mmin /1,2,4,7,11,16/
      data mmax /1,3,6,10,15,21/
      data dzero /0.0d0/
      pint = dzero
      qint = dzero
      rint = dzero
      npts = (ni+nj-2)/2 + 1
      imin = mmin(npts)
      imax = mmax(npts)
      do 140 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)*t
         ptx = dum + p0
         pty = dum + q0
         ptz = dum + r0
         ax = ptx - pi
         ay = pty - qi
         az = ptz - ri
         bx = ptx - pj
         by = pty - qj
         bz = ptz - rj
         go to (60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      go to (130,120,110,100,90,80,70) , nj
 70      px = px*bx
         py = py*by
         pz = pz*bz
 80      px = px*bx
         py = py*by
         pz = pz*bz
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     pint = pint + px
         qint = qint + py
         rint = rint + pz
 140  continue
      return
      end
**==scross.f
      subroutine scross(sout,mdim)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension sout(mdim,mdim)
INCLUDE(common/segm)
INCLUDE(common/timez)
INCLUDE(common/restar)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/zorac)
      common/bufb/s(225),g(225),pp,u(7),w(7),
     *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,
     *tol,nroots,ni,nj,ii,jj,lit,ljt,mini,minj,maxi,maxj,oiandj
     *,imap(mxprim)
      common/junk/ilifc3(maxorb),ntran3(maxorb),itran3(mxorb3),
     * ctran3(mxorb3),itrn3(2),
     *ex2(mxprim),cs2(mxprim),cp2(mxprim),cd2(mxprim),cf2(mxprim),
     *cg2(mxprim),cspace(maxat),
     +kstrt2(mxshel),ktom2(mxshel),ktype2(mxshel),kng2(mxshel),
     +kloc2(mxshel),kmin2(mxshel),kmax2(mxshel),nshel2,nspace(3),
     *nprin1,itol1,icut1,normf1,normp1
     *,nprini(695),iovmat,iosvec,ioproj,iorthg,c2(3,maxat)
      common/blkin/gg(225),ft(225),dij(225),
     * pin(225),qin(225),rin(225),
     * ijx(225),ijy(225),ijz(225)
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      data pi32/5.56832799683170d0/
      data pt75/0.75d0/
      data dzero,pt5,done,two,three,five,seven /0.0d0,0.5d0,1.0d0,
     + 2.0d0,3.0d0,5.0d0,7.0d0/
      data rnine/9.0d0/
      data eleven /11.0d0/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data rln10 /2.30258d0/
      data tol2/1.0d-10/
      data pt187,pt6562 /1.875d0,6.5625d0/
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     +          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     +          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     +          2, 0, 2, 1, 1/
      data ix / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     +          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     +          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     +          0, 2, 1, 2, 1/
      data iy / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     +          1,16, 1, 6, 1,11,11, 1, 6, 6,
     +          1,21, 1, 6, 1,16,16, 1, 6,11,
     +          1,11, 6,11, 6/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     +          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     +          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     +          2, 2, 1, 1, 2/
      data iz / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     +          1, 1,16, 1, 6, 1, 6,11,11, 6,
     +          1, 1,21, 1, 6, 1, 6,16,16, 1,
     +         11,11, 6, 6,11/
      tol = rln10*itol
      onorm = normf.ne.1 .or. normp.ne.1
      onorm1 = normf1.ne.1 .or. normp1.ne.1
      if ((onorm .and. .not.onorm1) .or. (onorm1 .and. .not.onorm).and.
     +    .not.oint_zora)
     +    call caserr2('basis function normalisation incompatible.')
      oiandj = .false.
      do 20 ishl2 = 1 , nshel2
         k1 = kstrt2(ishl2)
         imap(k1) = 0
 20   continue
      if (oint_zora) go to 71
      do 70 ishl2 = 1 , nshel2
         k1 = kstrt2(ishl2)
         k2 = k1 + kng2(ishl2) - 1
         if (imap(k1).eq.0) then
            imap(k1) = ishl2
            if (normp.ne.1) then
               do 30 ig = k1 , k2
                  ee = ex2(ig) + ex2(ig)
                  facs = pi32/(ee*dsqrt(ee))
                  facp = pt5*facs/ee
                  facd = pt75*facs/(ee*ee)
                  facf = pt187*facs/(ee**3)
                  facg = pt6562*facs/(ee**4)
                  cs2(ig) = cs2(ig)/dsqrt(facs)
                  cp2(ig) = cp2(ig)/dsqrt(facp)
                  cd2(ig) = cd2(ig)/dsqrt(facd)
                  cf2(ig) = cf2(ig)/dsqrt(facf)
                  cg2(ig) = cg2(ig)/dsqrt(facg)
 30            continue
            end if
            if (normf.ne.1) then
               facs = dzero
               facp = dzero
               facd = dzero
               facf = dzero
               facg = dzero
               do 50 ig = k1 , k2
                  do 40 jg = k1 , ig
                     ee = ex2(ig) + ex2(jg)
                     fac = ee*dsqrt(ee)
                     dums = cs2(ig)*cs2(jg)/fac
                     dump = pt5*cp2(ig)*cp2(jg)/(ee*fac)
                     dumd = pt75*cd2(ig)*cd2(jg)/(ee**2*fac)
                     dumf = pt187*cf2(ig)*cf2(jg)/(ee**3*fac)
                     dumg = pt6562*cg2(ig)*cg2(jg)/(ee**4*fac)
                     if (ig.ne.jg) then
                        dums = dums + dums
                        dump = dump + dump
                        dumd = dumd + dumd
                        dumf = dumf + dumf
                        dumg = dumg + dumg
                     end if
                     facs = facs + dums
                     facp = facp + dump
                     facd = facd + dumd
                     facf = facf + dumf
                     facg = facg + dumg
 40               continue
 50            continue
               do 60 ig = k1 , k2
                  if (facs.gt.tol2) cs2(ig) = cs2(ig)/dsqrt(facs*pi32)
                  if (facp.gt.tol2) cp2(ig) = cp2(ig)/dsqrt(facp*pi32)
                  if (facd.gt.tol2) cd2(ig) = cd2(ig)/dsqrt(facd*pi32)
                  if (facf.gt.tol2) cf2(ig) = cf2(ig)/dsqrt(facf*pi32)
                  if (facg.gt.tol2) cg2(ig) = cg2(ig)/dsqrt(facg*pi32)
 60            continue
            end if
         end if
 70   continue
 71   continue
      i = mdim*mdim
      call vclr(sout,1,i)
      do 370 ii = 1 , nshell
         i = katom(ii)
         pi = c(1,i)
         qi = c(2,i)
         ri = c(3,i)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
         do 360 jj = 1 , nshel2
            j = ktom2(jj)
            pj = c2(1,j)
            qj = c2(2,j)
            rj = c2(3,j)
            if (oatint_z) then
c...           atomic zora 1-center
               if (katom(ii).ne.ktom2(jj)) go to 360
               pi = 0.0d0
               qi = 0.0d0
               ri = 0.0d0
               pj = 0.0d0
               qj = 0.0d0
               rj = 0.0d0
            end if
            j1 = kstrt2(jj)
            j2 = j1 + kng2(jj) - 1
            ljt = ktype2(jj)
            minj = kmin2(jj)
            maxj = kmax2(jj)
            locj = kloc2(jj) - minj
            nroots = (lit+ljt-2)/2 + 1
            rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
            ij = 0
            mmax = maxj
            do 90 i = mini , maxi
               nnx = ix(i)
               nny = iy(i)
               nnz = iz(i)
               do 80 j = minj , mmax
                  ij = ij + 1
                  ijx(ij) = nnx + jx(j)
                  ijy(ij) = nny + jy(j)
                  ijz(ij) = nnz + jz(j)
                  if (j.le.1) then
                     ft(ij) = three
                  else if (j.le.4) then
                     ft(ij) = five
                  else if (j.le.10) then
                     ft(ij) = seven
                  else if (j.gt.20) then
                     ft(ij) = eleven
                  else
                     ft(ij) = rnine
                  end if
 80            continue
 90         continue
            do 100 i = 1 , ij
               s(i) = dzero
 100        continue
            jgmax = j2
            do 330 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*pi
               ayi = ai*qi
               azi = ai*ri
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
               do 320 jg = j1 , jgmax
                  aj = ex2(jg)
                  aa = ai + aj
                  aa1 = done/aa
                  dum = aj*arri*aa1
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs2(jg)
                     cpj = cp2(jg)
                     cdj = cd2(jg)
                     cfj = cf2(jg)
                     cgj = cg2(jg)
                     ax = (axi+aj*pj)*aa1
                     ay = (ayi+aj*qj)*aa1
                     az = (azi+aj*rj)*aa1
                     mmax = maxj
                     nn = 0
                     do 280 i = mini , maxi
                        go to (110,120,180,180,
     +                         130,180,180,140,180,180,
     +                         150,180,180,160,180,180,180,180,180,170,
     +                         172,180,180,174,180,180,180,180,180,
     +                         176,180,180,178,180,180), i
 110                    dum1 = csi*fac
                        go to 180
 120                    dum1 = cpi*fac
                        go to 180
 130                    dum1 = cdi*fac
                        go to 180
 140                    if (onorm) dum1 = dum1*sqrt3
                        go to 180
 150                    dum1 = cfi*fac
                        go to 180
 160                    if (onorm) dum1 = dum1*sqrt5
                        go to 180
 170                    if (onorm) dum1 = dum1*sqrt3
                        go to 180
 172                    dum1 = cgi*fac
                        go to 180
 174                    if (onorm) dum1 = dum1*sqrt7
                        go to 180
 176                    if (onorm) dum1 = dum1*sqrt5/sqrt3
                        go to 180
 178                    if (onorm) dum1 = dum1*sqrt3
 180                    do 270 j = minj , mmax
                           go to (190,200,260,260,
     +                      210,260,260,220,260,260,
     +                      230,260,260,240,260,260,260,260,260,250,
     +                      252,260,260,254,260,260,260,260,
     +                      260,256,260,260,258,260,260),j
 190                       dum2 = dum1*csj
                           go to 260
 200                       dum2 = dum1*cpj
                           go to 260
 210                       dum2 = dum1*cdj
                           go to 260
 220                       if (onorm) dum2 = dum2*sqrt3
                           go to 260
 230                       dum2 = dum1*cfj
                           go to 260
 240                       if (onorm) dum2 = dum2*sqrt5
                           go to 260
 250                       if (onorm) dum2 = dum2*sqrt3
                           go to 260
 252                       dum2 = dum1*cgj
                           go to 260
 254                       if (onorm) dum2 = dum2*sqrt7
                           go to 260
 256                       if (onorm) dum2 = dum2*sqrt5/sqrt3
                           go to 260
 258                       if (onorm) dum2 = dum2*sqrt3
 260                       nn = nn + 1
                           dij(nn) = dum2
 270                    continue
 280                 continue
                     t = dsqrt(aa1)
                     t1 = -two*aj*aj*t
                     t2 = -pt5*t
                     p0 = ax
                     q0 = ay
                     r0 = az
                     in = -5
                     do 300 i = 1 , lit
                        in = in + 5
                        ni = i
                        do 290 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call scrint
                           pin(jn) = pint*t
                           qin(jn) = qint*t
                           rin(jn) = rint*t
                           nj = j + 2
                           call scrint
                           pin(jn+25) = pint*t1
                           qin(jn+25) = qint*t1
                           rin(jn+25) = rint*t1
                           nj = j - 2
                           if (nj.gt.0) then
                              call scrint
                           else
                              pint = dzero
                              qint = dzero
                              rint = dzero
                           end if
                           n = (j-1)*(j-2)
                           dum = dfloat(n)*t2
                           pin(jn+50) = pint*dum
                           qin(jn+50) = qint*dum
                           rin(jn+50) = rint*dum
 290                    continue
 300                 continue
                     do 310 i = 1 , ij
                        nnx = ijx(i)
                        nny = ijy(i)
                        nnz = ijz(i)
                        cyz = qin(nny)*rin(nnz)
                        dum = cyz*pin(nnx)
                        dum1 = (pin(nnx+25)+pin(nnx+50))*cyz 
     +                       + (qin(nny+25)+qin(nny+50))*pin(nnx)
     +                         *rin(nnz) + (rin(nnz+25)+rin(nnz+50))
     +                         *pin(nnx)*qin(nny)
                        s(i) = s(i) + dij(i)*dum
 310                 continue
                  end if
 320           continue
 330        continue
            mmax = maxj
            nn = 0
            do 350 i = mini , maxi
               li = loci + i
               do 340 j = minj , mmax
                  lj = locj + j
                  nn = nn + 1
                  sout(li,lj) = s(nn)
 340           continue
 350        continue
 360     continue
 370  continue
      return
      end
**==shalfd.f
      subroutine shalfd(s,v,n,maxao)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     make s**(-1/2) for transformation to orthonormal basis.
c.......................................................................
      common/blkin/scr1(100),scr2(100)
INCLUDE(common/iofile)
      dimension s(*), v(n,n)
      icount = 0
      nnd = n*(n+1)/2
      nn = n*n
      if (n.gt.maxao) call caserr2('dimensioning error in shalfd')
      do 30 i = 1 , n
         do 20 j = 1 , i
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
 20      continue
         v(i,i) = 1.0d0
 30   continue
      call jacod(s,v,n,nnd,nn,1,n,scr1,scr2,n)
       icount = 0
       do i=1,n
          if (s(i*(i+1)/2).lt.1.0d-10) icount = icount + 1
       end do
       if (icount.gt.0) then
          write(iwr,33) icount,(s(i*(i+1)/2),i=1,n)
33        format(' *** ',i5,' eigenvalues of s-matrix under 1.0d-10',
     1           ' **dependency**',/'  eigenvalues',/,(1x,10(d10.3),1x))
          call caserr('dependency detected in atomscf')
       end if
      icount = 1
      do 40 i = 1 , n
         scr1(i) = 1.0d0/dsqrt(s(icount))
         icount = icount + i + 1
 40   continue
      icount = 0
      do 70 i = 1 , n
         do 60 j = 1 , i
            icount = icount + 1
            hlp = 0.0d0
            do 50 k = 1 , n
               hlp = hlp + v(i,k)*scr1(k)*v(j,k)
 50         continue
            s(icount) = hlp
 60      continue
 70   continue
      icount = 0
      do 90 i = 1 , n
         do 80 j = 1 , i
            icount = icount + 1
            v(j,i) = s(icount)
            v(i,j) = s(icount)
 80      continue
 90   continue
      return
      end
**==split.f
      subroutine split(v,nsplit,nmo,na,l1)
c
c     ----- split some of the occupied orbitals into pairs -----
c
c     nmo   = number of mo-s read in
c     nsplit = number to be split
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension v(l1,*)
      common/junk/ limlow(maxat),limsup(maxat)
      if (nsplit.eq.0) return
      call aolim2
      ii = 0
      ilo = nmo - nsplit + 1
      do 30 i = ilo , nmo
         ii = ii + 1
         iat = 1
         rmax = 0.0d0
         do 20 j = 1 , na
            il = limlow(j)
            if (il.ne.0) then
               iu = limsup(j)
               k = iu - il + 1
               rmaxs = ddot(k,v(il,i),1,v(il,i),1)
               if (rmax.le.rmaxs) then
                  iat = j
                  rmax = rmaxs
               end if
            end if
 20      continue
         il = limlow(iat)
         iu = limsup(iat)
         call dcopy(l1,v(1,i),1,v(1,ii+nmo),1)
         a = -1.0d0
         j = iu - il + 1
         call dscal(j,a,v(il,ii+nmo),1)
 30   continue
      return
      end
**==starcd.f
      subroutine starcd(c,ss,nbci,ncshi,noshi)
      implicit REAL  (a-h,o-z),integer   (i-n)
      dimension ss(*), c(nbci,nbci)
      common/blkin/a(100)
c.......................................................................
c
c     this routine provides a set of schmidt
c     orthogonalized start vectors.
c
c     inline indexing function for triangular matrices.
c
c.......................................................................
      index(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c.......................................................................
c
c     uneducated guess of trial vectors.
c
c.......................................................................
      nshti = ncshi + noshi
      if (nshti.eq.0) return
      if (nshti.gt.nbci) call caserr2(
     +  'more atomic shells than basis functions in atomic startup')
      ixx = nbci/nshti
      khi = 0
      do 30 j = 1 , nshti
         klo = khi + 1
         khi = khi + ixx
         do 20 k = klo , khi
            c(k,j) = 1.0d0
 20      continue
 30   continue
c.......................................................................
c
c     orthogonalize.
c
c.......................................................................
      do 130 j = 1 , nshti
c.......................................................................
c
c     take scalar products with preceding vectors of same symm.
c     first do a(*) = ss(*)*c(*,j).
c
c.......................................................................
         do 50 l = 1 , nbci
            sum = 0.0d0
            do 40 m = 1 , nbci
               indxlm = index(l,m)
               sum = sum + ss(indxlm)*c(m,j)
 40         continue
            a(l) = sum
 50      continue
c.......................................................................
c
c     if first vector of symmetry, normalize directly.
c
c.......................................................................
         if (j.ne.1) then
c.......................................................................
c
c     c(*,k)*a(*)
c
c.......................................................................
            jm1 = j - 1
            do 80 k = 1 , jm1
               sum = 0.0d0
               do 60 l = 1 , nbci
                  sum = sum + c(l,k)*a(l)
 60            continue
c.......................................................................
c
c     multiply c(*,k) by the scalar product and subtract from c(*,j)
c
c.......................................................................
               do 70 l = 1 , nbci
                  c(l,j) = c(l,j) - c(l,k)*sum
 70            continue
 80         continue
c.......................................................................
c
c     c(*,j) is orthogonal to all preceding vectors. normalize.
c
c.......................................................................
            do 100 l = 1 , nbci
               sum = 0.0d0
               do 90 m = 1 , nbci
                  indxlm = index(l,m)
                  sum = sum + ss(indxlm)*c(m,j)
 90            continue
               a(l) = sum
 100        continue
         end if
         sum = 0.0d0
         do 110 l = 1 , nbci
            sum = sum + c(l,j)*a(l)
 110     continue
         sum = dsqrt(sum)
         do 120 l = 1 , nbci
            c(l,j) = c(l,j)/sum
 120     continue
 130  continue
      return
      end
**==swap.f
      subroutine swap(q,e,nbas,nswap,iswap)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/iofile)
      dimension q(nbas,*),e(*),iswap(*)
      ii = 1
      do 20 imo = 1 , nswap
         j = iswap(ii)
         k = iswap(ii+1)
         write (iwr,6010) j , k
         temp = e(j)
         e(j) = e(k)
         e(k) = temp
_IF1(uv)         do 30 i=1,nbas
_IF1(uv)         temp=q(i,j)
_IF1(uv)         q(i,j)=q(i,k)
_IF1(uv)   30    q(i,k)=temp
_IFN1(uv)          call dswap(nbas,q(1,j),1,q(1,k),1)
         ii = ii + 2
 20   continue
      return
 6010 format (/' mos ',i4,' and',i4,' exchanged')
      end
**==symsqs.f
      subroutine symsqs(s,sv,se,sh,b,ia,n,n2,ndim)
c
c     ----- calculate the symmetric s **(+1/2) matrix
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension s(*),sv(ndim,*),se(*),ia(*),sh(ndim,*),b(*)
c
c     ----- diagonalize the overlap matrix -----
c
      data dzero,done/0.0d0,1.0d0/
      call gldiag(n,n,ndim,s,b,se,sv,ia,2)
c
c     ----- take the square root of the eigenvalues
c
c
c     ----- form s **(+1/2) -----
c
_IF1(civ)      do 40 i = 1,n
_IF1(civ)   40 se(i) = dsqrt(se(i))
_IF1(u)      call vsqrtv(n,se,se)
_IFN1(civu)      call vsqrt(se,1,se,1,n)
      call vclr(s,1,n2)
      do 30 k = 1 , n
         dum = done/se(k)
         ij = 1
         do 20 i = 1 , n
            fact = sv(i,k)*dum
            if (fact.ne.dzero) then
               call daxpy(i,fact,sv(1,k),1,s(ij),1)
            end if
            ij = ij + i
 20      continue
 30   continue
      call square(sh,s,n,n)
      return
      end
**==tback.f
      subroutine tback(qnew,ilifn,q,ilifq,nnn)
c
c...  create vectors in symmetry adapted basis from unadapted vectors
c...  inverse of tdown , anorm must have been called to fill tranao
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension qnew(*),q(*),ilifn(*),ilifq(*)
      dimension dum(maxorb)
INCLUDE(common/tran)
INCLUDE(common/tranao)
INCLUDE(common/infoa)
c
      if (otran) then
         do 20 i = 1 , nnn
            call dcopy(num,q(ilifq(i)+1),1,qnew(ilifn(i)+1),1)
 20      continue
         return
      else
         do 50 i = 1 , nnn
            m = ilifq(i)
            call vclr(dum,1,num)
            do 40 j = 1 , num
               n = ntrad(j)
               do 30 k = 1 , n
                  l = ilifd(j) + k
                  dum(itrad(l)) = ctrad(l)*q(m+j) + dum(itrad(l))
 30            continue
 40         continue
            call dcopy(num,dum(1),1,qnew(ilifn(i)+1),1)
 50      continue
         return
      end if
      end
**==teigd.f
      subroutine teigd(pcap,qcap,u,t,dt,dos)
      implicit REAL  (a-h,o-z),integer   (i-n)
      logical open,klnemn
c.......................................................................
c
c     two-electron integral routine for s,p,d, and f functions.
c.......................................................................
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension twopow(14), pifac(14)
      dimension pcap(*),qcap(*),u(*),t(*),dos(*),dt(*)
_IF(ccpdft)
      REAL facex
INCLUDE(common/ccpdft.hf77)
INCLUDE(common/restrl)
_ENDIF
c.......................................................................
c
c     angular factors for exchange integrals. obtained as sum of squares
c     of slater coefficients c(kappa;l1,m1;l2,m2) divided by
c     2*(2*l1+1)*(2*l2+1) for any kappa, l1, and l2.
c.......................................................................
      data ss0,sp1,pp0,pp2,sd2,pd1,pd3,dd0,dd2,
     x     sf3,pf2,pf4,df1,df3,df5,ff0,ff2,ff4,ff6
     x/ .50000000000d+00, .16666666667d+00, .16666666667d+00,
     x  .66666666667d-01, .10000000000d+00, .66666666667d-01,
     x  .42857142857d-01, .10000000000d+00, .28571428571d-01,
     x                    .71428571429d-01, .42857142857d-01,
     x  .31746031746d-01, .42857142857d-01, .19047619048d-01,
     x  .21645021645d-01, .71428571429d-01, .19047619048d-01,
     x  .12987012987d-01, .16650016650d-01/
      f0pol(a,b) = 3.0d0*(16.0d0*a**6+
     +             104.0d0*a**5*b+286.0d0*a**4*b**2+
     +             429.0d0*(a*b)**3+
     +             286.0d0*a**2*b**4+104.0d0*a*b**5+16.0d0*b**6)
      f2pol(a,b) = 8.0d0*a**4 + 52.0d0*a**3*b + 143.0d0*(a*b)**2 + 
     +            52.0d0*a*b**3 + 8.0d0*b**4
      df1pol(a,b) = 8.0d0*a**4 + 44.0d0*a**3*b + 99.0d0*(a*b)**2 + 
     +             44.0d0*a*b**3 +  8.0d0*b**4
      df3pol(a,b) = 2.0d0*a**2 + 11.0d0*a*b + 2.0d0*b**2
c.......................................................................
c     two-electron integral routine for lcgo atom scf.
c     restricted to principal quantum numbers 1,2 and 3 for respectively
c     s,p, and d orbitals.
c.......................................................................
c     pi = 3.14159265d0
c     on6 = 1.0d0/6.0d0
c     on15 = 1.0d0/15.d0
c     on35 = 2*on70
c     on70 = 1.0d0/70.d0
      pifac(1) = dsqrt(3.14159265d0)
      twopow(1) = 1
      do 20 i = 2 , 14
         twopow(i) = twopow(i-1)*0.5d0
         pifac(i) = pifac(i-1)*dfloat(2*i-1)
 20   continue
_IF(ccpdft)
      if (CD_active().and.oatmdft) then
         if (CD_HF_exchange()) then
            facex = CD_HF_exchange_weight()
         else
            facex = 0.0d0
         endif
      else
         facex = 1.0d0
      endif
_ENDIF
c.......................................................................
c
c     this part sets up the coefficients lambda,p,q and mu,r,s.
c.......................................................................
      j = 0
c     nstep1 = 0
      kmx = 0
      factkl = 1.0d0
      kl = 0
      pot = 0.0d0
      potn = 0.0d0
      cin = 0.0d0
      do 190 i = 1 , nsym
         prfac1 = twopow(i+1)*factkl
         kin = kmx + 1
         kmx = kin + nbas(i) - 1
         do 180 k = kin , kmx
            zp = zeta(k)
            do 170 l = kin , k
               kl = kl + 1
               pcap(kl) = 0.0d0
               qcap(kl) = 0.0d0
               zq = zeta(l)
               zpq = zp + zq
               prfac2 = prfac1*u(kl)
               xfac1 = prfac2*zpq**i*pifac(i)
c              nstep2 = 0
               mmx = 0
               factmn = 1.0d0
               mn = 0
               do 160 im = 1 , i
                  open = (nosh(i).ne.0 .and. nosh(im).ne.0)
                  prfac3 = prfac2*pifac(im)*twopow(im)*factmn
                  xfac2 = xfac1*factmn*twopow(im+1)
                  minx = mmx + 1
                  mmx = minx - 1 + nbas(im)
                  mmxp = mmx
                  if (im.eq.i) mmxp = k
                  do 150 m = minx , mmxp
                     zr = zeta(m)
                     zpr = zp + zr
                     zqr = zq + zr
                     nmx = l
                     if (m.lt.k) nmx = m
                     do 140 n = minx , nmx
                        mn = mn + 1
                        klnemn = (kl.ne.mn)
                        j = j + 1
c.......................................................................
c
c     j is the number label of the matrix elements to be calculated
c     i=lambda+1,k=p,l=q,im=mu+1,m=r,n=s
c.......................................................................
                        zs = zeta(n)
                        zrs = zr + zs
                        zqs = zq + zs
                        zps = zp + zs
                        zpqrs = zpq + zrs
                        zpqrs2 = 2.0d0*zpqrs**2
                        zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        xterm = (1.0d0/dsqrt(zpqrs))**(2*(i+im)-3)
                        prfac4 = prfac3*u(mn)*xterm
                        xfac3 = xfac2*u(mn)*xterm
                        xfac11 = xfac3*(zrs/zprzqs)**im
                        xfac21 = xfac3*(zrs/zpszqr)**im
                        xfsum = xfac11 + xfac21
                        ntest = i*(i-1)/2 + im
                        go to (30,40,50,60,70,80,90,100,110,120,140) ,
     +                         ntest
c.......................................................................
c
c     i=1,im=1,(ss)-loop. x0=j0(ss),y0=k0(ss)
c.......................................................................
 30                     x0 = prfac4
                        y0 = xfsum
_IF(ccpdft)
                        pj = x0 - y0*ss0*facex
                        qj = -ajmn(1)*y0*facex
_ELSE
                        pj = x0 - y0*ss0
                        qj = -ajmn(1)*y0
_ENDIF
                        go to 130
c.......................................................................
c
c     i=2,im=1,(sp)-loop. x0=j0(sp),y1=k1(sp)
c.......................................................................
 40                     x0 = prfac4*(3*zpq+2*zrs)
                        y1 = xfsum
_IF(ccpdft)
                        pj = x0 - y1*sp1*facex
                        qj = -ajmn(2)*y1*facex
_ELSE
                        pj = x0 - y1*sp1
                        qj = -ajmn(2)*y1
_ENDIF
                        go to 130
c.......................................................................
c
c     i=2,im=2,(pp)-loop. x0=j0(pp),y0=k0(pp),y2=k2(pp)
c.......................................................................
 50                     x0 = prfac4*(zpqrs2+zpq*zrs)
                        y0 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y0 = y0 + xfsum
                        y2 = xfsum*5
_IF(ccpdft)
                        pj = x0 - (y0*pp0 + y2*pp2)*facex
                        qj = -(ajmn(3)*y0+ajmn(4)*y2)*facex
_ELSE
                        pj = x0 - y0*pp0 - y2*pp2
                        qj = -(ajmn(3)*y0+ajmn(4)*y2)
_ENDIF
                        go to 130
c.......................................................................
c
c     i=3,im=1,(sd)-loop. x0=j0(sd),y2=k2(sd)
c.......................................................................
 60                     x0 = prfac4*(15.0d0*zpq**2+20.0d0*zpq*zrs
     +                               +8.0d0*zrs**2)
                        y2 = xfsum
_IF(ccpdft)
                        pj = x0 - y2*sd2*facex
                        qj = -ajmn(5)*y2*facex
_ELSE
                        pj = x0 - y2*sd2
                        qj = -ajmn(5)*y2
_ENDIF
                        go to 130
c.......................................................................
c
c     i=3,im=2,(pd)-loop. x0=j0(pd),y1=k1(pd),y3=k3(pd),x2=j2(pd)
c.......................................................................
 70                     x0 = prfac4*(10.0d0*zpq**3+
     +                       35.0d0*zpq**2*zrs+
     +                       28.0d0*zpq*zrs**2+8.0d0*zrs**3)
                        y1 = xfsum*zpqrs2
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
                        y1 = y1 + 3.0d0*xfsum
                        y3 = 7.0d0*xfsum
                        qj = 0.d0
_IF(ccpdft)
                        pj = x0 - (y1*pd1 + y3*pd3)*facex
_ELSE
                        pj = x0 - y1*pd1 - y3*pd3
_ENDIF
                        if (open) then
                           x2 = prfac4*5.0d0*zpq*zrs*
     +                                (7.0d0*zpq+2*zrs)
_IF(ccpdft)
                           qj = ajmn(21)*x2 - (ajmn(6)*y1+ajmn(7)*y3)
     +                                        *facex
_ELSE
                           qj = ajmn(21)*x2 - (ajmn(6)*y1+ajmn(7)*y3)
_ENDIF
                        end if
                        go to 130
c.......................................................................
c
c     i=3,im=3,(dd)-loop. x0=j0(dd),y0=k0(dd),y2=k2(dd),y4=k4(dd)
c.......................................................................
 80                     zprzqs = zpr*zqs
                        zpszqr = zps*zqr
                        x0 = prfac4*((zpqrs2+zpq*zrs)
     +                        *zpqrs2*2+7.0d0*(zpq*zrs)**2)
                        y01 = xfac11*((zpqrs2+zprzqs)
     +                        *zpqrs2*2+7.0d0*zprzqs**2)
                        y02 = xfac21*((zpqrs2+zpszqr)
     +                        *zpqrs2*2+7.0d0*zpszqr**2)
                        y0 = y01 + y02
                        xfac11 = xfac11*7.0d0*zprzqs
                        xfac21 = xfac21*7.0d0*zpszqr
                        y21 = xfac11*(zpqrs2+5.0d0*zprzqs)
                        y22 = xfac21*(zpqrs2+5.0d0*zpszqr)
                        y2 = y21 + y22
                        y4 = (xfac11*zprzqs+xfac21*zpszqr)*9.0d0
_IF(ccpdft)
                        pj = x0 - (y0*dd0 + (y2+y4)*dd2)*facex
                        qj = -ajmn(8)*y0 - ajmn(9)*y2 - ajmn(10)*y4
                        qj = qj*facex
_ELSE
                        pj = x0 - y0*dd0 - (y2+y4)*dd2
                        qj = -ajmn(8)*y0 - ajmn(9)*y2 - ajmn(10)*y4
_ENDIF
                        go to 130
c.......................................................................
c
c     i=4,im=1,(sf)-loop. x0=j0(sf),y3=k3(sf)
c.......................................................................
 90                     x0 = prfac4*(105.0d0*zpq**3+210.0d0*zpq**2*zrs+
     +                       168.0d0*zpq*zrs**2+48.0d0*zrs**3)
                        y3 = xfsum
_IF(ccpdft)
                        pj = x0 - y3*sf3*facex
                        qj = -ajmn(11)*y3*facex
_ELSE
                        pj = x0 - y3*sf3
                        qj = -ajmn(11)*y3
_ENDIF
                        go to 130
c.......................................................................
c
c     i=4,im=2,(pf)-loop. x0=j0(pf),y2=k2(pf),y4=k4(pf),x2=j2(pf)
c.......................................................................
 100                    x0 = prfac4*(70.0d0*zpq**4+315.0d0*zpq**3*zrs+
     +                       378.0d0*(zpq*zrs)**2+216.0d0*zpq*zrs**3+
     +                        48.0d0*zrs**4)
c.......................................................................
c                       y2 = xfsum*zpqrs2
c.......................................................................
                        y2 = xfac11*(2.0d0*zpr**2+
     +                               9.0d0*zpr*zqs+2*zqs**2)
     +                     + xfac21*(2.0d0*zps**2+
     +                               9.0d0*zps*zqr+2*zqr**2)
                        xfsum = xfac11*zpr*zqs + xfac21*zps*zqr
c.......................................................................
c                       y2 = y2 + 5*xfsum
c.......................................................................
                        y4 = 9.0d0*xfsum
                        qj = 0.0d0
_IF(ccpdft)
                        pj = x0 - (y2*pf2 + y4*pf4)*facex
_ELSE
                        pj = x0 - y2*pf2 - y4*pf4
_ENDIF
                        if (open) then
                           x2 = prfac4*5.0d0*zpq*zrs*(63.0d0*zpq**2+
     +                          26.0d0*zpq*zrs+ 8.0d0*zrs**2)
_IF(ccpdft)
                           qj = ajmn(22)*x2 - (ajmn(12)*y2+ajmn(13)*y4)
     +                                        *facex
_ELSE
                           qj = ajmn(22)*x2 - ajmn(12)*y2 - ajmn(13)*y4
_ENDIF
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=3,(df)-loop. x0=j0(df),y1=k1(df),y3=k3(df),y5=k5(df)
c                         x2=j2(df),x4=j4(df).
c.......................................................................
 110                    x0 = prfac4*(56.0d0*zpq**5+308.0d0*zpq**4*zrs+
     +                       693.0d0*zpq**3*zrs**2+
     +                       594.0d0*zpq**2*zrs**3+
     +                       264.0d0*zpq*zrs**4+48.0d0*zrs**5)
                        y1 = xfac11*df1pol(zpr,zqs)
     +                       + xfac21*df1pol(zps,zqr)
                        xfac11 = 9.0d0*zpr*zqs*xfac11
                        xfac21 = 9.0d0*zps*zqr*xfac21
                        y3 = xfac11*df3pol(zpr,zqs)
     +                       + xfac21*df3pol(zps,zqr)
                        y5 = 11.0d0*(xfac11*zpr*zqs+xfac21*zps*zqr)
_IF(ccpdft)
                        pj = x0 - (y1*df1 + y3*df3 + y5*df5)*facex
_ELSE
                        pj = x0 - y1*df1 - y3*df3 - y5*df5
_ENDIF
                        qj = 0.0d0
                        if (open) then
                           prfac4 = prfac4*7.0d0*zpq*zrs
                           x2 = prfac4*(18.0d0*zpq**3+99.0d0*zpq**2*zrs+
     +                          44.0d0*zpq*zrs**2+8.0d0*zrs**3)
                           x4 = 9.0d0*prfac4*(11.0d0*zpq+2*zrs)*zpq*zrs
_IF(ccpdft)
                           qj = x2*ajmn(23) + x4*ajmn(24) - (y1*ajmn(14)
     +                          + y3*ajmn(15) + y5*ajmn(16))*facex
_ELSE
                           qj = x2*ajmn(23) + x4*ajmn(24) - y1*ajmn(14)
     +                          - y3*ajmn(15) - y5*ajmn(16)
_ENDIF
                        end if
                        go to 130
c.......................................................................
c
c     i=4,im=4,(ff)-loop. x0=j0(ff),y0=k0(ff),y2=k2(ff),y4=k4(ff)
c              y6=k6(ff)
c.......................................................................
 120                    x0 = prfac4*f0pol(zpq,zrs)
                        y0 = xfac11*f0pol(zpr,zqs)
     +                       + xfac21*f0pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*9.0d0
                        xfac21 = xfac21*zps*zqr*9.0d0
                        y2 = xfac11*f2pol(zpr,zqs)
     +                       + xfac21*f2pol(zps,zqr)
                        xfac11 = xfac11*zpr*zqs*11
                        xfac21 = xfac21*zps*zqr*11
                        y4 = xfac11*(2.0d0*zpr**2+13.0d0*zpr*zqs+
     +                       2.0d0*zqs**2) + xfac21*(2.0d0*zps**2+
     +                      13.0d0*zps*zqr +2.0d0*zqr**2)
                        y6 = 13.0d0*(xfac11*zpr*zqs+xfac21*zps*zqr)
_IF(ccpdft)
                        pj = x0 - (y0*ff0 + y2*ff2 + y4*ff4 + y6*ff6)
     +                            *facex
_ELSE
                        pj = x0 - y0*ff0 - y2*ff2 - y4*ff4 - y6*ff6
_ENDIF
                        qj = -y0*ajmn(17) - y2*ajmn(18) - y4*ajmn(19)
     +                       - y6*ajmn(20)
_IF(ccpdft)
                        qj = qj*facex
_ENDIF
 130                    pcap(kl) = pcap(kl) + pj*dt(mn)
                        if (klnemn) pcap(mn) = pcap(mn) + pj*dt(kl)
                        term = dt(kl)*pj*dt(mn)
                        if (open) then
                           qcap(kl) = qcap(kl) + qj*dos(mn)
                           if (klnemn) qcap(mn) = qcap(mn) + qj*dos(kl)
                           term = term - dos(kl)*qj*dos(mn)
                        end if
                        if (.not.klnemn) term = term*0.5d0
                        pot = pot + term
 140                 continue
 150              continue
                  factmn = factmn/im
 160           continue
               potn = potn + u(kl)*dt(kl)
               cin = cin + t(kl)*dt(kl)
 170        continue
 180     continue
         factkl = factkl/i
 190  continue
c     istart = 1
      pot = pot - zn*potn
      energ = cin + pot
      vir = pot/cin
      return
      end
**==tracd.f
      subroutine tracd(c,cc,nsqt)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     transform vectors from contr. to primitive basis functions.
c.......................................................................
      parameter (nbig=1500, no=50)
      common /junk/nsym,nbas(5),ncsh(5),nosh(5),nccup(6),nsht,nitscf,
     + zn,n1(6),nqn(nbig),zeta(nbig),eps(no),cin,vir,energ,ajmn(24),
     + damp,nconv,nbc(5),cont(nbig),nbct,nstrt(nbig),ifcont,iatom
      dimension c(*),cc(*)
c
      do 20 i = 1 , nsqt
         c(i) = 0.0d0
 20   continue
c.......................................................................
c
c     cc(i,j)  contains vectors over cont. functions.
c              i runs over functions, j runs over orbitals
c.......................................................................
      nstep = 0
      nstep1 = 0
      nstep2 = 0
      do 60 l = 1 , nsym
         nbc1 = nbc(l)
         nbasl = nbas(l)
         nsh1 = ncsh(l) + nosh(l)
         ii = 1
         do 50 n = 1 , nbc1
            k1 = nstrt(nstep+n)
            k2 = nstrt(nstep+n+1) - 1
            do 40 k = k1 , k2
               do 30 j = 1 , nsh1
                  c(nstep2+ii+(j-1)*nbasl) = cc(nstep1+n+(j-1)*nbc1)
     +               *cont(k)
 30            continue
               ii = ii + 1
 40         continue
 50      continue
         nstep = nstep + nbc1
         nstep1 = nstep1 + nbc1**2
         nstep2 = nstep2 + nbasl**2
 60   continue
      return
      end
**==trafsd.f
      subroutine trafsd(nsym,nbas,ndim,a,nbc,cont,nstrt,ffc)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c     transform matrix a
c     from prim. to cont. basis.
c.......................................................................
      dimension ffc(*),nbas(*),a(*),nbc(*),cont(*),nstrt(*)
      do 20 i = 1 , ndim
         ffc(i) = a(i)
 20   continue
      mnstep = 0
      ijstep = 0
      k = 1
      ijbas = 0
      mnbas = 0
      do 70 l = 1 , nsym
         nbc1 = nbc(l)
         do 60 m = 1 , nbc1
            do 50 n = 1 , m
               m1 = nstrt(mnbas+m)
               n0 = nstrt(mnbas+n)
               m2 = nstrt(mnbas+m+1) - 1
               n2 = nstrt(mnbas+n+1) - 1
               sumc = 0.0d0
               do 40 i = m1 , m2
                  ij = ijstep + (i-ijbas-1)*(i-ijbas)/2
                  if (m.eq.n) then
                     n2 = i
                  end if
                  do 30 j = n0 , n2
                     jj = j - ijbas
                     if (m.ne.n) then
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     else if (i.ne.j) then
                        sumc = sumc + 2.0d0*cont(i)*cont(j)*ffc(ij+jj)
                     else
                        sumc = sumc + cont(i)*cont(j)*ffc(ij+jj)
                     end if
 30               continue
 40            continue
               a(k) = sumc
               k = k + 1
 50         continue
 60      continue
         mnstep = mnstep + nbc(l)*(nbc(l)+1)/2
         ijstep = ijstep + nbas(l)*(nbas(l)+1)/2
         mnbas = mnbas + nbc(l)
         ijbas = ijbas + nbas(l)
 70   continue
      return
      end
**==tramad.f
      subroutine tramad(a,b,c,mdima,mdim,scr)
      implicit REAL  (a-h,o-z),integer   (i-n)
c.......................................................................
c
c     this routine tranforms a to (b+)ab, where a is a lower
c     triangular symmetric matrix, and b is a square matrix.
c     the transformed matrix is returned in a.
c.......................................................................
      dimension a(mdima), b(mdim,mdim), c(mdima), scr(mdim)
      do 60 i = 1 , mdim
c.......................................................................
c
c     generate i-th column of ab.
c.......................................................................
         do 30 j = 1 , mdim
            sum = 0.0d0
            do 20 l = 1 , mdim
               maxjl = max(j,l)
               jl = maxjl*(maxjl-3)/2 + j + l
               sum = sum + a(jl)*b(l,i)
 20         continue
            scr(j) = sum
 30      continue
c.......................................................................
c
c     multiply this by rows of b+
c.......................................................................
         do 50 j = 1 , i
            sum = 0.0d0
            do 40 k = 1 , mdim
               sum = sum + scr(k)*b(k,j)
 40         continue
            c(i*(i-1)/2+j) = sum
 50      continue
 60   continue
c.......................................................................
c
c     transfer to a for return.
c.......................................................................
      do 70 i = 1 , mdima
         a(i) = c(i)
 70   continue
      return
      end
**==tranpq.f
      subroutine tranpq(p1,p2,mdim,nbasis,mbasis)
c
c ---- transforms a non-symmetric matrix to ctrans basis
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension p1(mdim,mdim),p2(mdim,mdim)
      common/blkin/ctx(mxorb3),itx(mxorb3)
INCLUDE(common/tran)
      common/junk/ilifc2(maxorb),nterm2(maxorb),
     * it2(mxorb3),ct2(mxorb3),otran2
      mdim2 = mdim*mdim
      if (otran) then
         do 20 i = 1 , nbasis
            ntran(i) = 1
            ilifc(i) = i - 1
            itran(i) = i
            ctran(i) = 1.0d0
 20      continue
      end if
      if (otran2) then
         do 30 i = 1 , mbasis
            nterm2(i) = 1
            ilifc2(i) = i - 1
            it2(i) = i
            ct2(i) = 1.0d0
 30      continue
      end if
c
c ---- zero result
c
c
c ---- loop over functions in basis set 1
c
      call vclr(p2,1,mdim2)
      do 80 i = 1 , nbasis
         nt1 = ntran(i)
         do 40 k = 1 , nt1
            l = ilifc(i) + k
            ctx(k) = ctran(l)
            itx(k) = itran(l)
 40      continue
c
c ---- loop over functions in basis set 2
c
         do 70 j = 1 , mbasis
            nt2 = nterm2(j)
            top = 0.0d0
            do 60 k = 1 , nt2
               kk = k + ilifc2(j)
               l2 = it2(kk)
               do 50 ii = 1 , nt1
                  l1 = itx(ii)
                  top = top + ctx(ii)*ct2(kk)*p1(l1,l2)
 50            continue
 60         continue
            p2(i,j) = top
 70      continue
 80   continue
      return
      end
**==trn5d.f
      subroutine trn5d(h,l2,i0,eval)
c
      implicit REAL  (a-h,o-z)
c
      dimension h(l2)
c
c     ----- eliminate the s contaminant from a d shell -----
c
c     on entry, eval is the energy assigned to the d xy,xz,yz orbs
c               i0 is the row in -h- immediately before the d xx orb
c
c       a 6 d basis must be fixed to generate 5 degenerate,
c       negative energy orbitals, and a s orbital with positive energy
c       the symmetric 3x3 matrix involving xx,yy,and zz equal to
c       2a/3+b
c       -a/3+b  2a/3+b
c       -a/3+b  -a/3+b  2a/3+b
c       has eigenvalues and eigenvectors of
c        a == xx-yy
c        a == xx+yy-2*zz
c       3b == xx+yy+zz
c       let a be the energy of the xy,xz,yz d orbitals, and b be +3
c       so the s contaminant gets a large positive energy.
c
c
c     ediag= 0.66666e+00*eval+3.0e+00
c     eoff =-0.33333e+00*eval+3.0e+00
      ediag = (4.0d0/9.0d0)*eval + 2.0d0
      eoff = -(2.0d0/9.0d0)*eval + 2.0d0
      ii = (i0*i0+i0)/2
      do 20 i = 1 , 3
         ii = ii + i + i0
         h(ii) = ediag
         if (i.ne.1) then
            ij = ii - 1
            h(ij) = eoff
            if (i.ne.2) then
               ij = ij - 1
               h(ij) = eoff
            end if
         end if
 20   continue
      return
      end
**==ttran.f
      subroutine ttran(h,f,q,t,ia,m,nsplit,n)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension h(*),f(n,*),q(n,*),t(*),ia(*)
      mp1 = m + 1
      do 40 jj = 1 , nsplit
         j = jj + m
         do 20 i = 1 , n
            t(i) = ddot(n,f(1,i),1,q(1,j),1)
 20      continue
         do 30 i = mp1 , j
            dum = ddot(n,q(1,i),1,t,1)
            h(ia(jj)+i-m) = dum
 30      continue
 40   continue
      return
      end
**==vorth.f
      subroutine vorth(v,s,vt,m,nsplit,l1)
c
c     ----- generate symmetrically orthogonalized vectors
c     using the transformation s**(1/2)  (in s) and
c     the vectors v -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(l1,*),s(nsplit,*),vt(l1)
      do 30 i = 1 , l1
         do 20 j = 1 , nsplit
            vt(j) = ddot(nsplit,v(i,m+1),l1,s(1,j),1)
 20      continue
         call dcopy(nsplit,vt,1,v(i,m+1),l1)
 30   continue
      return
      end
      subroutine ver_guess(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/guess.m,v $
     +     "/
      data revision /"$Revision: 6321 $"/
      data date /"$Date: 2015-03-15 01:39:24 +0100 (Sun, 15 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
      subroutine swap_mcf(vect,mswap,nswap)
      implicit integer (a-z)
      REAL  vect(*)
      integer nswap(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/mapper)
INCLUDE(common/machin)
INCLUDE(common/iofile)
      common/craypk/mmmm(65),isymaos(maxorb),isymmos(maxorb)
INCLUDE(common/restri)
      logical logsym
      character*8 wordt
      data m51/51/
      nav = lenwrd()
c
      call secget(isect(490),m51,iblk51)
_IFN1(civfuk)      call readi(mmmm,mach(13)*nav,iblk51,idaf)
_IF1(civfuk)      call rdedx(mmmm,mach(13),iblk51,idaf)
      do 35 i=1,num
      ii=1+ilifq(i)
      ibig=idamax(num,vect(ii),1)
      isymmos(i)=isymaos(ibig)
35    continue
c
      nswapin = mswap
      if (nswapin.eq.0) goto 1
      jbase = 80 - 2*nswapin
      do 5 isw = 1,2*nswapin
      nswap(jbase+isw) = nswap(isw)
5     continue
c
1     itag = 0
      isw = 0
      iorb2 = 0
20    itag = itag + 1
      iorb = iorb2
      wordt = zorb(itag)
      if (inporb(wordt,repcnt,clef,cod,sym,logsym).ne.0) goto 100
      if (.not.logsym) call caserr2('error 1 in swap_mcf')
      iorb1 = iorb + 1
      iorb2 = iorb + repcnt
      iorbr = iorb
50    jorb = iorbr
      iorbr = iorbr + 1
      if (iorbr.gt.iorb2) goto 20
30    jorb = jorb + 1
      jsym = isymmos(jorb)
      if (jsym.eq.0) call caserr2('error 2 in swap_mcf')
      if (jsym.ne.sym) goto 30
      if (jorb.eq.iorbr) goto 50
      do 40 korb = iorbr+1,jorb
      isw = isw + 1
      nswap(isw) = iorbr
      isw = isw + 1
      nswap(isw) = korb
40    continue
      isymmos(jorb) = isymmos(iorbr)
      isymmos(iorbr) = sym
      goto 50
c
100   if (nswapin.eq.0) goto 200
      do 110 jsw = 1,2*nswapin
      isw = isw + 1
      nswap(isw) = nswap(jbase+jsw)
110   continue
200   mswap = isw/2
      if (mswap.gt.40) call caserr2('mcscf/force : too many swaps')
      return
      end
