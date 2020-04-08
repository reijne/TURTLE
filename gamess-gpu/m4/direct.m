c 
c  $Author: jmht $
c  $Date: 2013-12-05 14:38:55 +0100 (Thu, 05 Dec 2013) $
c  $Locker:  $
c  $Revision: 6291 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/direct.m,v $
c  $State: Exp $
c  
      subroutine bodmp2(iscat2,zmem,zmemc,xijk,xija,xijt,xijp1,xijp2
     &,xijp3,nij,xklk,xkla,xklt,xklp1,xklp2,xklp3,igath,nkl
     &,orbmo,iocc,nxx,irpos,ispos
     &,ext1,ext2,ext3,i1,i2,j1,j2,k1,k2,l1,l2,prefac,xinner
     &,cindmp,iindmp,tindmp,intdoq,iwr)
c
c**** routine calculates integrals prior to mp2 calculation; little
c**** different in structure from body
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto),zmemc(icsize,memcto)
     &,xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,xkla(nkl),xklk(nkl),xklp1(nkl),xklp2(nkl),xklp3(nkl)
     &,ab(3),naind(3),ncind(3),pa(3),igath(nkl),xijt(nij),xklt(nkl)
     &,m1(48),m2(48),m3(48)
      dimension orbmo(nb,nb),iscat2(nkl)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/gen/maxit,iter,concrit,conv,energy,thresh,dmax
c
      common/tempit/zitest,zipass,zjtest,zjpass
      common/scra  /iso(mxshel,48),nt,nttttt
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/defunk/cobas(maxorb,3),nd
      logical ext1,ext2,ext3,ext4,ext5
c
_IF1(c)cdir$  list
_IF1(c)cdir$  novector
      if (nbreak.ne.1) then
         write (iwr,*) 'code assumes the above'
      end if
      ivlx = nvl(ix)
      jvlx = nvl(jx)
      kvlx = nvl(kx)
      lvlx = nvl(lx)
c     ipos = nposf(i1,ivlx,nd,np,ns)
c     jpos = nposf(j1,jvlx,nd,np,ns)
c     kpos = nposf(k1,kvlx,nd,np,ns)
c     lpos = nposf(l1,lvlx,nd,np,ns)
      memreq = ivlx*jvlx*kvlx*lvlx
      jmax = j2
      icorig = 0
      ipq = 0
c
c**** the main contracted loop over i and j shells
c
      do 400 ic = i1 , i2
         if (ext1) jmax = ic
         do 390 jc = j1 , jmax
c           dij = 4*1.0d0
c           codens = 1.4d0
            ext5 = ext1 .and. ic.eq.jc
            itprim = nc(ic)*nc(jc)
            if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
c***** symmetry check on jshell
            do 20 it = 2 , nt
               ids = iso(ic,it)
               jds = iso(jc,it)
               if (ids.lt.jds) then
                  nds = ids
                  ids = jds
                  jds = nds
               end if
               m1(it) = ids
               m2(it) = jds
 20         continue
c
            ab(1) = cobas(ic,1) - cobas(jc,1)
            ab(2) = cobas(ic,2) - cobas(jc,2)
            ab(3) = cobas(ic,3) - cobas(jc,3)
            rab = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
c           jpmax = nc(jc)
c
c*****  now the loop over the number of different kl
c*****  mini-batches in this particular maxi-batch
c
            do 380 nbrk = 1 , nbreak
               ext4 = .false.
               kmax = info(2,nbrk)
               lmax = info(4,nbrk)
               lmin = info(3,nbrk)
               kmin = info(1,nbrk)
               if (ext3) then
                  if (kmin.gt.ic) go to 380
                  if (kmax.ge.ic) then
                     if (kmin.eq.ic .and. lmin.gt.jc) go to 380
                     if (kmax.ne.ic .or. lmax.ge.jc) then
                        ext4 = .true.
                        kmax = ic
                        lmax = jc
                     end if
                  end if
               end if
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
               ncr = 0
               icount = icorig
               xijtm = 0.0d0
               do 30 ijprm = 1 , itprim
                  icount = icount + 1
                  xijtm = max(xijtm,xijt(icount))
 30            continue
c
               test1 = prefac/xijtm
               tout = tmax(nbrk)
               if (test1.lt.tout) then
                  test2 = xinner/xijtm
c                 tmax1 = tmax(nbrk)*10.0d0
                  jstart = info(5,nbrk)
                  if (nt.eq.1) then
                     call tstmp2(test2,xklt,
     +                           igath,zmem(1,4),kmin,kmax,
     +                           l1,l2,lmax,lmin,jstart,nkl,ipsize,ext2,
     +                           jcount,jc0,jc1)
                  else
                     call tst1sm(test2,xklt,
     +                           igath,zmem(1,4),ic,jc,kmin,kmax,
     +                           m1,m2,m3,xklk,l1,l2,lmax,lmin,jstart,
     +                           nkl,ipsize,ext2,jcount,jc0,jc1)
                  end if
                  zipass = zipass + jc0
                  zitest = zitest + jcount
                  if (jc0.ne.0) then
                     jstart = jstart + 1
c
c***** prepare all parameter arrays needed to form
c***** restricted list of intermediate primitive integrals
c
      call dgthr(jc0,xklp1(jstart),zmem(1,1),igath(1))
      call dgthr(jc0,xklp2(jstart),zmem(1,2),igath(1))
      call dgthr(jc0,xklp3(jstart),zmem(1,3),igath(1))
      call dgthr(jc0,xkla(jstart),zmem(1,7),igath(1))
      if (nt.eq.1)  then
      call dgthr(jc0,xklk(jstart),zmem(1,9),igath(1))
      endif
                     do 40 kl = 1 , jc0
                        zmem(kl,4) = zmem(kl,1) - zmem(kl,4)
                        zmem(kl,5) = zmem(kl,2) - zmem(kl,5)
                        zmem(kl,6) = zmem(kl,3) - zmem(kl,6)
 40                  continue
                     icount = icorig
c                    jcold = jc0
                     do 200 ijprm = 1 , itprim
                        icount = icount + 1
c
c***** the outer loop integral test
c
                        test1 = prefac/xijt(icount)
                        tout = tmax(nbrk)
                        if (test1.lt.tout) then
                           zjpass = zjpass + 1
c
                           pa(1) = xijp1(icount) - cobas(ic,1)
                           pa(2) = xijp2(icount) - cobas(ic,2)
                           pa(3) = xijp3(icount) - cobas(ic,3)
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c
c***** now do the loops for forming the intermediate
c***** primitive integrals
c
c**** call routine to form (ss|ss)m integrals + a few
c**** useful quantities.
c
c*** z(1-3) qi ::: z(4-6) qi-ci ::: z(7) c+d :::
c*** z(8) 1/(a+b+c+d) ::: z(9) xklk ::: z(10-12) wi-pi :::
c*** z(13-15) wi-qi ::: z(16-?) fm(t) :::
c
                           call sinter(zmem(1,1),zmem(1,2),zmem(1,3),
     +                                 zmem(1,13),zmem(1,7),zmem(1,8),
     +                                 zmem(1,9),zmem(1,10),zmem(1,11),
     +                                 zmem(1,12),zmem(1,16),
     +                                 xijp1(icount),xijp2(icount),
     +                                 xijp3(icount),xija(icount),
     +                                 xijk(icount),mx,jc0,ipsize,
     +                                 zmem(1,mx+16))
c
c*** loop over number of distinct integral classes
c*** needed to form the target integrals
c
                           memcur = 15 + mx
                           do 160 n = 1 , ninty
                              ityp = intgrl(1,n)
                              ktyp = intgrl(2,n)
                              ityp1 = ityp - 1
                              ktyp1 = ktyp - 1
                              ityp2 = ityp - 2
                              ktyp2 = ktyp - 2
                              ivl = nvl(ityp)
                              kvl = nvl(ktyp)
                              if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
                              if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
                              nrecur = intgrl(3,n)
c
                              mem1 = intgrl(4,n)
                              mem2 = intgrl(5,n)
                              mem3 = intgrl(6,n)
                              mem4 = intgrl(7,n)
                              mem5 = intgrl(8,n)
c
c now:
c loop over number of indices associated with i position
                              do 150 i = 1 , ivl
                                 nr = nrecur
                                 if (nr.lt.5) then
c***** reduction is at the i position
c
                                    index = ispdf1(i,1)
                                    indwp = index + 9
                                    natred = ispdf1(i,2)
                                    itemp = (natred-1)*kvl
                                    jmem1 = itemp + mem1
                                    jmem2 = itemp + mem2
                                    jmem5 = mem5 + (natred-1)*kvl1
                                    if (nr.ne.1 .and. nr.ne.3) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                                       nai = ixyz(index,natred,ityp2)
                                       if (nai.eq.0) then
                                         nr = nr - 1
                                         go to 50
                                       end if
                                       naind(2) = ixyz(2,natred,ityp2)
     +                                    + 1
                                       naind(3) = ixyz(3,natred,ityp2)
     +                                    + 1
                                       naind(index) = naind(index) - 1
                                       idoubr = iarray(naind(2),naind(3)
     +                                    )
                                       itemp2 = (idoubr-1)*kvl
                                       jmem3 = itemp2 + mem3
                                       jmem4 = itemp2 + mem4
                                    end if
                                 else
c
c**** reduction is at the k position
                                    itemp = (i-1)*kvl1
                                    jmem1 = itemp + mem1
                                    jmem2 = itemp + mem2
                                    itemp = (i-1)*kvl2
                                    jmem3 = itemp + mem3
                                    jmem4 = itemp + mem4
                                    if (nr.ne.5 .and. nr.ne.6) then
                                       naind(1) = ixyz(1,i,ityp1) + 1
                                       naind(2) = ixyz(2,i,ityp1) + 1
                                       naind(3) = ixyz(3,i,ityp1) + 1
                                    end if
                                 end if
c
c now:
c loop over number of indices associated with k position
c
 50                              do 140 k = 1 , kvl
                                    memcur = memcur + 1
                                    if (nr.lt.5) then
                                       imem1 = jmem1 + k
                                       imem2 = jmem2 + k
                                    else
                                       index = ispdf1(k,1)
                                       indqc = index + 3
                                       indwq = index + 12
                                       natred = ispdf1(k,2)
                                       imem1 = jmem1 + natred
                                       imem2 = jmem2 + natred
                                    end if
                                    go to (60,70,80,90,100,110,120,130)
     +                                 , nr
 60                                 call recur1(zmem(1,memcur),pa(index)
     +                                 ,zmem(1,imem1),zmem(1,indwp),
     +                                 zmem(1,imem2),jc0)
                                    go to 140
 70                                 imem3 = jmem3 + k
                                    imem4 = jmem4 + k
                                    call recur2(zmem(1,memcur),pa(index)
     +                                 ,zmem(1,imem1),zmem(1,indwp),
     +                                 zmem(1,imem2),xija(icount),nai,
     +                                 zmem(1,imem3),zmem(1,7),zmem(1,8)
     +                                 ,zmem(1,imem4),jc0)
                                    go to 140
 80                                 nci = ixyz(index,k,ktyp1)
                                    if (nci.eq.0) go to 60
                                    naind(2) = ixyz(2,k,ktyp1) + 1
                                    naind(3) = ixyz(3,k,ktyp1) + 1
                                    naind(index) = naind(index) - 1
                                    imem5 = jmem5 +
     +                                 iarray(naind(2),naind(3))
                                    call recur3(zmem(1,memcur),pa(index)
     +                                 ,zmem(1,imem1),zmem(1,indwp),
     +                                 zmem(1,imem2),zmem(1,8),nci,
     +                                 zmem(1,imem5),jc0)
                                    go to 140
 90                                 nci = ixyz(index,k,ktyp1)
                                    if (nci.eq.0) go to 70
                                    naind(2) = ixyz(2,k,ktyp1) + 1
                                    naind(3) = ixyz(3,k,ktyp1) + 1
                                    naind(index) = naind(index) - 1
                                    imem5 = jmem5 +
     +                                 iarray(naind(2),naind(3))
                                    imem3 = jmem3 + k
                                    jmem4 = jmem4 + k
                                    call recur4(zmem(1,memcur),pa(index)
     +                                 ,zmem(1,imem1),zmem(1,indwp),
     +                                 zmem(1,imem2),xija(icount),nai,
     +                                 zmem(1,imem3),zmem(1,7),zmem(1,8)
     +                                 ,zmem(1,imem4),nci,zmem(1,imem5),
     +                                 jc0)
                                    go to 140
 100                                call recur5(zmem(1,memcur),
     +                                 zmem(1,indqc),zmem(1,imem1),
     +                                 zmem(1,indwq),zmem(1,imem2),jc0)
                                    go to 140
 110                                nci = ixyz(index,natred,ktyp2)
                                    if (nci.eq.0) go to 100
                                    ncind(2) = ixyz(2,natred,ktyp2) + 1
                                    ncind(3) = ixyz(3,natred,ktyp2) + 1
                                    ncind(index) = ncind(index) - 1
                                    imem3 = jmem3 +
     +                                 iarray(ncind(2),ncind(3))
                                    imem4 = jmem4 +
     +                                 iarray(ncind(2),ncind(3))
                                    call recur6(zmem(1,memcur),
     +                                 zmem(1,indqc),zmem(1,imem1),
     +                                 zmem(1,indwq),zmem(1,imem2),
     +                                 zmem(1,7),nci,zmem(1,imem3),
     +                                 xija(icount),zmem(1,8),
     +                                 zmem(1,imem4),jc0)
                                    go to 140
 120                                if (naind(index).eq.1) go to 100
                                    naind(index) = naind(index) - 1
                                    imem5 = (iarray(naind(2),naind(3))
     +                                 -1)*kvl1 + natred + mem5
                                    call recur7(zmem(1,memcur),
     +                                 zmem(1,indqc),zmem(1,imem1),
     +                                 zmem(1,indwq),zmem(1,imem2),
     +                                 zmem(1,8),naind(index),
     +                                 zmem(1,imem5),jc0)
                                    naind(index) = naind(index) + 1
                                    go to 140
 130                                if (naind(index).eq.1) go to 110
                                    nci = ixyz(index,natred,ktyp2)
                                    if (nci.eq.0) go to 120
                                    ncind(2) = ixyz(2,natred,ktyp2) + 1
                                    ncind(3) = ixyz(3,natred,ktyp2) + 1
                                    ncind(index) = ncind(index) - 1
                                    imem3 = jmem3 +
     +                                 iarray(ncind(2),ncind(3))
                                    imem4 = jmem4 +
     +                                 iarray(ncind(2),ncind(3))
                                    naind(index) = naind(index) - 1
                                    imem5 = (iarray(naind(2),naind(3))
     +                                 -1)*kvl1 + natred + mem5
                                    call recur8(zmem(1,memcur),
     +                                 zmem(1,indqc),zmem(1,imem1),
     +                                 zmem(1,indwq),zmem(1,imem2),
     +                                 zmem(1,7),nci,zmem(1,imem3),
     +                                 xija(icount),zmem(1,8),
     +                                 zmem(1,imem4),naind(index),
     +                                 zmem(1,imem5),jc0)
                                    naind(index) = naind(index) + 1
 140                             continue
 150                          continue
 160                       continue
c
c***** now the primitive intermediate integrals are complete
c***** for one particular primitive pair of iprim,jprim.
c***** now the integrals are expanded to a complete list and
c***** and accumulated over all the primitives in one ic,jc
c***** contracted pairing to yield half-contracted integrals.
c
                           ncr = ncr + 1
                           do 190 n = 1 , ncon
                              imem1 = intcon(3,n)
                              imem2 = intcon(4,n)
                              memcur = 0
                   do 180 i = 1 , nvl(intcon(1,n))
                      do 170 k = 1 , nvl(intcon(2,n))
                         memcur = memcur + 1
                         imem2m = imem2 + memcur
                         imem1m = imem1 + memcur
                         if (ncr.eq.1) then
                            if (ext4) then
              call vclr(zmem(jc0,imem2m),1,jc0+itprim)
                            end if
         call dcopy(jc0,zmem(1,imem1m),1,zmem(1,imem2m),1)
                         else
         call vadd(zmem(1,imem2m),1
     +            ,zmem(1,imem1m),1
     +            ,zmem(1,imem2m),1
     +            ,jc0)
                         end if
 170                  continue
 180               continue
 190                       continue
                        end if
c
 200                 continue
c
c*****  now that the loops over the i and j primitives for the one
c*****  particular ic,jc contracted pairing have been closed, we
c*****  can contract the kl part of the integrals.
c*****  first we reopen the the kl integral loop
c
                     if (ncr.ne.0) then
_IF1(cu)          call gather(jc0,igath(1),iscat2(jstart),igath(1))
_IFN1(cu)          call igthr(jc0,iscat2(jstart),igath(1),igath(1))
                        call tst22(zmem,zmemc,igath(1),jc0,jc1,memcon)
c
c
c***** now we have all the contracted integrals. use the horizontal
c***** recursion relation (if necessary) to find the desired integral
c
                        ihzc = 0
                        if (jx.gt.1) then
c***** hrr must be applied on the ij position
                           jxtyp1 = jx - 1
                           ixtyp1 = ix - 1
                           do 270 kl = kx , klx
                              ihzc = ihzc + 1
                              mem1 = ihz(1,ihzc)
                              mem2 = ihz(2,ihzc)
                              mem3 = ihz(3,ihzc)
                              mem4 = ihz(4,ihzc)
                              mem5 = ihz(5,ihzc)
                              kvl = nvl(kl)
                              ijc = 0
                              do 260 i = 1 , ivlx
                                 naind(2) = ixyz(2,i,ixtyp1) + 1
                                 naind(3) = ixyz(3,i,ixtyp1) + 1
                                 do 250 j = 1 , jvlx
                                    ncind(2) = naind(2)
     +                                 + ixyz(2,j,jxtyp1)
                                    ncind(3) = naind(3)
     +                                 + ixyz(3,j,jxtyp1)
                                    ij = iarray(ncind(2),ncind(3))
                                    imem1 = ijc*kvl + mem1
                                    imem2 = (ij-1)*kvl + mem2
                                    ijc = ijc + 1
                                    if (rab.lt.1.0d-8) then
c********* when (ai-bi) is zero (a+b,0|c+d,0) = (a,b|c+d,0)
                                       do 210 k = 1 , kvl
         call dcopy(jc1,zmemc(1,imem2+k),1,zmemc(1,imem1+k),1)
 210                                   continue
                                       go to 250
                                    end if
                                    index1 = ispdf1(j,1)
                                    ncind(index1) = ncind(index1) - 1
                                    imem3a = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    if (jx.eq.2) then
*****  a p function is being added to the i/a position.
                                       do 220 k = 1 , kvl
                                         call hrz2(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),jc1)
 220                                   continue
                                       go to 250
                                    end if
                                    ired1 = ispdf1(j,2)
                                    index2 = ispdf1(ired1,1)
                                    ncind(index2) = ncind(index2) - 1
                                    imem4a = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index1) = ncind(index1) + 1
                                    imem3b = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    if (jx.eq.3) then
*****  a d function is being added to the i/a position
                                       do 230 k = 1 , kvl
                                         call hrz4(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),ab(index2)
     +                                      ,zmemc(1,imem3b+k),
     +                                      zmemc(1,imem4a+k),jc1)
 230                                   continue
                                       go to 250
                                    end if
                                    ired2 = ispdf1(ired1,2)
                                    index3 = ispdf1(ired2,1)
                                    ncind(index3) = ncind(index3) - 1
                                    imem4b = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index2) = ncind(index2) + 1
                                    imem3c = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    ncind(index1) = ncind(index1) - 1
                                    imem4c = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index2) = ncind(index2) - 1
                                    imem5 = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem5
                                    if (jx.eq.4) then
*****  an f function is being added to the i/a position
                                       do 240 k = 1 , kvl
                                         call hrz8(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),ab(index2)
     +                                      ,zmemc(1,imem3b+k),
     +                                      ab(index3),zmemc(1,imem3c+k)
     +                                      ,zmemc(1,imem4a+k),
     +                                      zmemc(1,imem4b+k),
     +                                      zmemc(1,imem4c+k),
     +                                      zmemc(1,imem5+k),jc1)
 240                                   continue
                                    end if
 250                             continue
 260                          continue
 270                       continue
c***** at this point the hrr on a and b was necessary and
c***** has been applied.
                        end if
                        if (lx.ne.1) then
c***** hrr is needed for c and d
c
c***** first form ci-di in mtemp position
                           jc1 = 0
                           mtemp1 = mtemp + 1
                           mtemp2 = mtemp + 2
                           do 290 kc = kmin , kmax
c*** a few logicals to determine extents of l loop
                              lmin1 = l1
                              lmax1 = l2
                              if (ext2) lmax1 = kc
                              if (kc.eq.kmax) lmax1 = lmax
                              if (kc.eq.kmin) lmin1 = lmin
                              do 280 lc = lmin1 , lmax1
                                 jc1 = jc1 + 1
                                 zmemc(jc1,mtemp) = cobas(kc,1)
     +                              - cobas(lc,1)
                                 zmemc(jc1,mtemp1) = cobas(kc,2)
     +                              - cobas(lc,2)
                                 zmemc(jc1,mtemp2) = cobas(kc,3)
     +                              - cobas(lc,3)
 280                          continue
 290                       continue
                           kxtyp1 = kx - 1
                           lxtyp1 = lx - 1
                           ihzc = ihzc + 1
                           klc = 0
                           mem1 = ihz(1,ihzc)
                           mem2 = ihz(2,ihzc)
                           mem3 = ihz(3,ihzc)
                           mem4 = ihz(4,ihzc)
                           mem5 = ihz(5,ihzc)
                           kl1 = kvlx*lvlx
                           kl2 = nvl(klx)
                           kl3 = nvl(klx-1)
                           kl4 = nvl(klx-2)
                           kl5 = nvl(klx-3)
                           do 370 k = 1 , kvlx
                              naind(2) = ixyz(2,k,kxtyp1) + 1
                              naind(3) = ixyz(3,k,kxtyp1) + 1
                              do 360 l = 1 , lvlx
                                 ncind(2) = ixyz(2,l,lxtyp1) + naind(2)
                                 ncind(3) = ixyz(3,l,lxtyp1) + naind(3)
                                 klc = klc + 1
                                 imem1 = mem1 + klc
                                 kl = iarray(ncind(2),ncind(3))
                                 imem2 = mem2 + kl
                                 index1 = ispdf1(l,1)
                                 indcd1 = mtemp - 1 + index1
                                 ncind(index1) = ncind(index1) - 1
                                 imem3a = mem3 +
     +                              iarray(ncind(2),ncind(3))
                                 ijc = -1
                                 if (lx.eq.2) then
*****  a p function is being added to the k/c position
                                    do 310 i = 1 , ivlx
                                       do 300 j = 1 , jvlx
                                         ijc = ijc + 1
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3 = imem3a + ijc*kl3
                                         call hrz2c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3),jc1)
 300                                   continue
 310                                continue
                                    go to 360
                                 end if
                                 ired1 = ispdf1(l,2)
                                 index2 = ispdf1(ired1,1)
                                 indcd2 = mtemp - 1 + index2
                                 ncind(index2) = ncind(index2) - 1
                                 imem4a = mem4 +
     +                              iarray(ncind(2),ncind(3))
                                 ncind(index1) = ncind(index1) + 1
                                 imem3b = mem3 +
     +                              iarray(ncind(2),ncind(3))
                                 if (lx.eq.3) then
c**** a d function is being added to the k/c position
                                    do 330 i = 1 , ivlx
                                       do 320 j = 1 , jvlx
                                         ijc = ijc + 1
                                         itemp = ijc*kl3
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3a = imem3a + itemp
                                         jmem3b = imem3b + itemp
                                         jmem4a = imem4a + ijc*kl4
                                         call hrz4c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3a),
     +                                      zmemc(1,indcd2),
     +                                      zmemc(1,jmem3b),
     +                                      zmemc(1,jmem4a),jc1)
 320                                   continue
 330                                continue
                                    go to 360
                                 end if
                                 ired2 = ispdf1(ired1,2)
                                 index3 = ispdf1(ired2,1)
                                 indcd3 = mtemp - 1 + index3
                                 ncind(index3) = ncind(index3) - 1
                                 imem4b = iarray(ncind(2),ncind(3))
     +                              + mem4
                                 ncind(index2) = ncind(index2) + 1
                                 imem3c = iarray(ncind(2),ncind(3))
     +                              + mem3
                                 ncind(index1) = ncind(index1) - 1
                                 imem4c = iarray(ncind(2),ncind(3))
     +                              + mem4
                                 ncind(index2) = ncind(index2) - 1
                                 imem5 = iarray(ncind(2),ncind(3))
     +                              + mem5
                                 if (lx.eq.4) then
c**** a d function is being added to the k/c position
                                    do 350 i = 1 , ivlx
                                       do 340 j = 1 , jvlx
                                         ijc = ijc + 1
                                         itemp = ijc*kl3
                                         itemp1 = ijc*kl4
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3a = imem3a + itemp
                                         jmem3b = imem3b + itemp
                                         jmem3c = imem3c + itemp
                                         jmem4a = imem4a + itemp1
                                         jmem4b = imem4b + itemp1
                                         jmem4c = imem4c + itemp1
                                         jmem5 = imem5 + ijc*kl5
                                         call hrz8c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3a),
     +                                      zmemc(1,indcd2),
     +                                      zmemc(1,jmem3b),
     +                                      zmemc(1,indcd3),
     +                                      zmemc(1,jmem3c),
     +                                      zmemc(1,jmem4a),
     +                                      zmemc(1,jmem4b),
     +                                      zmemc(1,jmem4c),
     +                                      zmemc(1,jmem5),jc1)
 340                                   continue
 350                                continue
                                 end if
 360                          continue
 370                       continue
                        end if
c***** the hrr has now been applied and the integrals lie in
c*****    zmemc(1,memans)
c***** send the integrals to four index transformation
                        kmaxt = k2 - k1 + 1
                        lmaxt = l2 - l1 + 1
                        intdoq = intdoq + 1
                        call indx4a(zmemc(1,memans+1),ivlx,jvlx,kvlx,
     +                              lvlx,kmaxt,lmaxt,ext2,orbmo,
     +                              iocc,nxx,ipq,irpos,ispos,memreq,
     +                              cindmp,iindmp,tindmp)
                     end if
                  end if
               end if
c
 380        continue
            icorig = icorig + itprim
            zjtest = zjtest + itprim
 390     continue
 400  continue
_IF1(c)cdir$  nolist
_IF1(c)cdir$  vector
      return
      end
      subroutine bodrv(iscat2,zmem,zmemc,xijk,xija,xijt,xijp1,xijp2
     &,xijp3,nij,xklk,xkla,xklt,xklp1,xklp2,xklp3,igath,nkl
     &,dens,codens,xklc,grad,grtemp,igath2
     &,ext1,ext2,ext3,i1,i2,j1,j2,k1,l1,l2,prefac,xinner)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c**** the fixed parameters in common/indexd/ are explained in the
c**** block data statement (in bdata.f?)
c**** the memory locations and requirements were established in
c**** cidrv and are stored in common/cider/ and cidriv
c
c**** routine calculates derivative integrals (see also body)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto),zmemc(icsize,memcto)
     &,xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,xkla(nkl),xklk(nkl),xklp1(nkl),xklp2(nkl),xklp3(nkl)
     &,ab(3),naind(3),ncind(3),pa(3),igath(nkl),xijt(nij),xklt(nkl)
     &,m0(48),m1(48),m2(48),m3(48)
     &,xklc(nkl),grad(na,3),grtemp(na,3)
      dimension codens(nb,nb),dens(nb,nb),iscat2(nkl),igath2(nkl)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/gen/maxit,iter,concrit,conv,energy,thresh,dmax
c
      common/tempit/zitest,zipass,zjtest,zjpass
      common/scra  /iso(mxshel,48),nt
c
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mtgrad
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/defunk/cobas(maxorb,3),nd
      logical ext1,ext2,ext3,ext4,ext5,ext7,centr1
c
_IF1(c)cdir$  list
_IF1(c)cdir$  novector
      ivlx = nvl(ix)
      jvlx = nvl(jx)
      kvlx = nvl(kx)
      lvlx = nvl(lx)
      ipos = nposf(i1,ivlx,nd,np,ns)
      jpos = nposf(j1,jvlx,nd,np,ns)
      kpos = nposf(k1,kvlx,nd,np,ns)
      lpos = nposf(l1,lvlx,nd,np,ns)
      memreq = ivlx*jvlx*kvlx*lvlx + 13
      jmax = j2
      icorig = 0
      centr1 = .false.
      if (i1.eq.j1 .and. i1.eq.k1 .and. i1.eq.l1) centr1 = .true.
c
c**** the main contracted loop over i and j shells
c
      do 690 ic = i1 , i2
         if (ext1) jmax = ic
c**** symetry test on ishell
         do 30 it = 2 , nt
            id = iso(ic,it)
            if (id.gt.ic) then
               do 20 jc = j1 , jmax
                  ext5 = ext1 .and. ic.eq.jc
                  itprim = nc(ic)*nc(jc)
                  if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
                  icorig = icorig + itprim
                  zjtest = zjtest + itprim
 20            continue
               go to 690
            end if
            m0(it) = id
 30      continue
         do 680 jc = j1 , jmax
            ext5 = ext1 .and. ic.eq.jc
            itprim = nc(ic)*nc(jc)
            if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
c***** symmetry check on jshell
            if (ext1) then
               do 40 it = 2 , nt
                  ids = m0(it)
                  jds = iso(jc,it)
                  if (jds.gt.ic) go to 670
                  if (ids.lt.jds) then
                     nds = ids
                     ids = jds
                     jds = nds
                  end if
                  if (ids.eq.ic .and. jds.gt.jc) go to 670
                  m1(it) = ids
                  m2(it) = jds
 40            continue
            else
               do 50 it = 2 , nt
                  ids = m0(it)
                  jds = iso(jc,it)
                  if (ids.eq.ic .and. jds.gt.jc) go to 670
                  m1(it) = ids
                  m2(it) = jds
 50            continue
            end if
c
            dij = 4*codens(ic,jc)
            ab(1) = cobas(ic,1) - cobas(jc,1)
            ab(2) = cobas(ic,2) - cobas(jc,2)
            ab(3) = cobas(ic,3) - cobas(jc,3)
            rab = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
            if (rab.ge.1.0d-8) centr1 = .false.
            if (.not.(centr1)) then
c              jpmax = nc(jc)
c
c*****  now the loop over the number of different kl
c*****  mini-batches in this particular maxi-batch
c
               do 660 nbrk = 1 , nbreak
                  ext4 = .false.
                  kmax = info(2,nbrk)
                  lmax = info(4,nbrk)
                  lmin = info(3,nbrk)
                  kmin = info(1,nbrk)
                  if (ext3) then
                     if (kmin.gt.ic) go to 660
                     if (kmax.ge.ic) then
                        if (kmin.eq.ic .and. lmin.gt.jc) go to 660
                        if (kmax.ne.ic .or. lmax.ge.jc) then
                           ext4 = .true.
                           kmax = ic
                           lmax = jc
                        end if
                     end if
                  end if
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
                  ncr = 0
                  icount = icorig
                  xijtm = 0.0d0
                  do 60 ijprm = 1 , itprim
                     icount = icount + 1
                     xijtm = max(xijtm,xijt(icount))
 60               continue
c
                  test1 = prefac/xijtm
                  tout = dmax*tmax(nbrk)
                  if (test1.lt.tout) then
                     test2 = xinner/xijtm
                     tmax1 = tmax(nbrk)*10.0d0
                     jstart = info(5,nbrk)
                     if (nt.eq.1) then
                        call tst1dr(codens,dij,dmaxij,tmax1,test2,
     +                              xklt(jstart+1),igath,zmem(1,4),ic,
     +                              jc,iscat2(jstart+1),kmin,
     +                              kmax,l1,l2,lmax,lmin,nkl,
     +                              ipsize,ext2,jcount,jc0,jc1)
                     else
                        call tst1sd(codens,dij,dmaxij,tmax1,test2,xklt,
     +                              igath,zmem(1,4),ic,jc,kmin,
     +                              kmax,m1,m2,m3,xklk,l1,l2,lmax,lmin,
     +                              jstart,nkl,ipsize,ext2,ext3,jcount,
     +                              jc0,jc1)
                     end if
                     if (jc0.ne.0) then
                        jstart = jstart + 1
c
c***** prepare all parameter arrays needed to form
c***** restricted list of intermediate primitive integrals
c
      call dgthr(jc0,xklp1(jstart),zmem(1,1),igath(1))
      call dgthr(jc0,xklp2(jstart),zmem(1,2),igath(1))
      call dgthr(jc0,xklp3(jstart),zmem(1,3),igath(1))
      call dgthr(jc0,xkla(jstart),zmem(1,7),igath(1))
      call dgthr(jc0,xklc(jstart),zmem(1,16),igath(1))
      if (nt.eq.1)  then
       call dgthr(jc0,xklk(jstart),zmem(1,9),igath(1))
      endif
                        do 70 kl = 1 , jc0
                           zmem(kl,4) = zmem(kl,1) - zmem(kl,4)
                           zmem(kl,5) = zmem(kl,2) - zmem(kl,5)
                           zmem(kl,6) = zmem(kl,3) - zmem(kl,6)
 70                     continue
                        icount = icorig
                        jcold = jc0
                        istar1 = nup(ic) - nc(ic) + 1
                        istar2 = nup(ic)
                        ijprim = 0
                        do 370 iprim = istar1 , istar2
                           alpha = pe(iprim)*2.0d0
                           do 360 jprim = 1 , nc(jc)
                              ijprim = ijprim + 1
                              if (ijprim.gt.itprim) go to 360
                              icount = icount + 1
c
c***** the outer loop integral test
c
                              test1 = prefac/xijt(icount)
                              tout = dmaxij*tmax(nbrk)
                              if (test1.ge.tout) go to 360
                              if (nbrk.eq.1) then
                                 zjpass = zjpass + 1
                              end if
c
                              pa(1) = xijp1(icount) - cobas(ic,1)
                              pa(2) = xijp2(icount) - cobas(ic,2)
                              pa(3) = xijp3(icount) - cobas(ic,3)
c
                              ext7 = .false.
                              if (ext4) then
                                 jc0 = jcold
                                 jctemp = jcold - itprim
                                 itxmax = itprim
                                 if (jctemp.lt.0) then
                                    jctemp = 0
                                    itxmax = jc0
                                 end if
_IF1(x)c$dir scalar
_IF1(t)cdir$ nextscalar
_IF1(c)ccdir$ nextscalar
                                 do 80 itx = 1 , itxmax
                                    jctemp = jctemp + 1
                                    itx1 = igath(jctemp) + jstart - 1
                                    if (itx1.ge.icount) then
                                       if (itx1.gt.icount) then
                                         jc0 = jctemp - 1
                                       else
                                         jc0 = jctemp
                                         ext7 = .true.
                                         zmem(jc0,9) = zmem(jc0,9)*0.5d0
                                       end if
                                       go to 90
                                    end if
 80                              continue
                              end if
 90                           zipass = zipass + jc0
                              zitest = zitest + jcount
c
c
c***** now do the loops for forming the intermediate
c***** primitive integrals
c
c**** call routine to form (ss|ss)m integrals + a few
c**** useful quantities.
c
c*** z(1-3) qi ::: z(4-6) qi-ci ::: z(7) c+d :::
c*** z(8) 1/(a+b+c+d) ::: z(9) xklk ::: z(10-12) wi-pi :::
c*** z(13-15) wi-qi ::: z(16) c ::: z(17-?) fm(t) :::
c
                              call sinter(zmem(1,1),zmem(1,2),zmem(1,3),
     +                           zmem(1,13),zmem(1,7),zmem(1,8),
     +                           zmem(1,9),zmem(1,10),zmem(1,11),
     +                           zmem(1,12),zmem(1,17),xijp1(icount),
     +                           xijp2(icount),xijp3(icount),
     +                           xija(icount),xijk(icount),mx+1,jc0,
     +                           ipsize,zmem(1,mx+18))
c
                              if (ext7) zmem(jc0,9) = zmem(jc0,9)*2.0d0
c
c*** loop over number of distinct integral classes
c*** needed to form the target integrals
c
                              memcur = 17 + mx
                              if (rab.lt.1.0d-8) ninty = ninty - 1
                              do 210 n = 1 , ninty
                                 ityp = intgrl(1,n)
                                 ktyp = intgrl(2,n)
                                 ityp1 = ityp - 1
                                 ktyp1 = ktyp - 1
                                 ityp2 = ityp - 2
                                 ktyp2 = ktyp - 2
                                 ivl = nvl(ityp)
                                 kvl = nvl(ktyp)
                                 if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
                                 if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
                                 nrecur = intgrl(3,n)
c
                                 mem1 = intgrl(4,n)
                                 mem2 = intgrl(5,n)
                                 mem3 = intgrl(6,n)
                                 mem4 = intgrl(7,n)
                                 mem5 = intgrl(8,n)
c
c now:
c loop over number of indices associated with i position
                                 do 200 i = 1 , ivl
                                    nr = nrecur
                                    if (nr.lt.5) then
c***** reduction is at the i position
c
                                       index = ispdf1(i,1)
                                       indwp = index + 9
                                       natred = ispdf1(i,2)
                                       itemp = (natred-1)*kvl
                                       jmem1 = itemp + mem1
                                       jmem2 = itemp + mem2
                                       jmem5 = mem5 + (natred-1)*kvl1
                                       if (nr.ne.1 .and. nr.ne.3) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                                         nai = ixyz(index,natred,ityp2)
                                         if (nai.eq.0) then
                                         nr = nr - 1
                                         go to 100
                                         end if
                                         naind(2) = ixyz(2,natred,ityp2)
     +                                      + 1
                                         naind(3) = ixyz(3,natred,ityp2)
     +                                      + 1
                                         naind(index) = naind(index) - 1
                                         idoubr = iarray(naind(2),
     +                                      naind(3))
                                         itemp2 = (idoubr-1)*kvl
                                         jmem3 = itemp2 + mem3
                                         jmem4 = itemp2 + mem4
                                       end if
                                    else
c
c**** reduction is at the k position
                                       itemp = (i-1)*kvl1
                                       jmem1 = itemp + mem1
                                       jmem2 = itemp + mem2
                                       itemp = (i-1)*kvl2
                                       jmem3 = itemp + mem3
                                       jmem4 = itemp + mem4
                                       if (nr.ne.5 .and. nr.ne.6) then
                                         naind(1) = ixyz(1,i,ityp1) + 1
                                         naind(2) = ixyz(2,i,ityp1) + 1
                                         naind(3) = ixyz(3,i,ityp1) + 1
                                       end if
                                    end if
c
c now:
c loop over number of indices associated with k position
c
 100                                do 190 k = 1 , kvl
                                       memcur = memcur + 1
                                       if (nr.lt.5) then
                                         imem1 = jmem1 + k
                                         imem2 = jmem2 + k
                                       else
                                         index = ispdf1(k,1)
                                         indqc = index + 3
                                         indwq = index + 12
                                         natred = ispdf1(k,2)
                                         imem1 = jmem1 + natred
                                         imem2 = jmem2 + natred
                                       end if
                                       go to (110,120,130,140,150,160,
     +                                    170,180) , nr
 110                                   call recur1(zmem(1,memcur),
     +                                    pa(index),zmem(1,imem1),
     +                                    zmem(1,indwp),zmem(1,imem2),
     +                                    jc0)
                                       go to 190
 120                                   imem3 = jmem3 + k
                                       imem4 = jmem4 + k
                                       call recur2(zmem(1,memcur),
     +                                    pa(index),zmem(1,imem1),
     +                                    zmem(1,indwp),zmem(1,imem2),
     +                                    xija(icount),nai,zmem(1,imem3)
     +                                    ,zmem(1,7),zmem(1,8),
     +                                    zmem(1,imem4),jc0)
                                       go to 190
 130                                   nci = ixyz(index,k,ktyp1)
                                       if (nci.eq.0) go to 110
                                       naind(2) = ixyz(2,k,ktyp1) + 1
                                       naind(3) = ixyz(3,k,ktyp1) + 1
                                       naind(index) = naind(index) - 1
                                       imem5 = jmem5 +
     +                                    iarray(naind(2),naind(3))
                                       call recur3(zmem(1,memcur),
     +                                    pa(index),zmem(1,imem1),
     +                                    zmem(1,indwp),zmem(1,imem2),
     +                                    zmem(1,8),nci,zmem(1,imem5),
     +                                    jc0)
                                       go to 190
 140                                   nci = ixyz(index,k,ktyp1)
                                       if (nci.eq.0) go to 120
                                       naind(2) = ixyz(2,k,ktyp1) + 1
                                       naind(3) = ixyz(3,k,ktyp1) + 1
                                       naind(index) = naind(index) - 1
                                       imem5 = jmem5 +
     +                                    iarray(naind(2),naind(3))
                                       imem3 = jmem3 + k
                                       imem4 = jmem4 + k
                                       call recur4(zmem(1,memcur),
     +                                    pa(index),zmem(1,imem1),
     +                                    zmem(1,indwp),zmem(1,imem2),
     +                                    xija(icount),nai,zmem(1,imem3)
     +                                    ,zmem(1,7),zmem(1,8),
     +                                    zmem(1,imem4),nci,
     +                                    zmem(1,imem5),jc0)
                                       go to 190
 150                                   call recur5(zmem(1,memcur),
     +                                    zmem(1,indqc),zmem(1,imem1),
     +                                    zmem(1,indwq),zmem(1,imem2),
     +                                    jc0)
                                       go to 190
 160                                   nci = ixyz(index,natred,ktyp2)
                                       if (nci.eq.0) go to 150
                                       ncind(2) = ixyz(2,natred,ktyp2)
     +                                    + 1
                                       ncind(3) = ixyz(3,natred,ktyp2)
     +                                    + 1
                                       ncind(index) = ncind(index) - 1
                                       imem3 = jmem3 +
     +                                    iarray(ncind(2),ncind(3))
                                       imem4 = jmem4 +
     +                                    iarray(ncind(2),ncind(3))
                                       call recur6(zmem(1,memcur),
     +                                    zmem(1,indqc),zmem(1,imem1),
     +                                    zmem(1,indwq),zmem(1,imem2),
     +                                    zmem(1,7),nci,zmem(1,imem3),
     +                                    xija(icount),zmem(1,8),
     +                                    zmem(1,imem4),jc0)
                                       go to 190
 170                                   if (naind(index).eq.1) go to 150
                                       naind(index) = naind(index) - 1
                                       imem5 =
     +                                    (iarray(naind(2),naind(3))-1)
     +                                    *kvl1 + natred + mem5
                                       call recur7(zmem(1,memcur),
     +                                    zmem(1,indqc),zmem(1,imem1),
     +                                    zmem(1,indwq),zmem(1,imem2),
     +                                    zmem(1,8),naind(index),
     +                                    zmem(1,imem5),jc0)
                                       naind(index) = naind(index) + 1
                                       go to 190
 180                                   if (naind(index).eq.1) go to 160
                                       nci = ixyz(index,natred,ktyp2)
                                       if (nci.eq.0) go to 170
                                       ncind(2) = ixyz(2,natred,ktyp2)
     +                                    + 1
                                       ncind(3) = ixyz(3,natred,ktyp2)
     +                                    + 1
                                       ncind(index) = ncind(index) - 1
                                       imem3 = jmem3 +
     +                                    iarray(ncind(2),ncind(3))
                                       imem4 = jmem4 +
     +                                    iarray(ncind(2),ncind(3))
                                       naind(index) = naind(index) - 1
                                       imem5 =
     +                                    (iarray(naind(2),naind(3))-1)
     +                                    *kvl1 + natred + mem5
                                       call recur8(zmem(1,memcur),
     +                                    zmem(1,indqc),zmem(1,imem1),
     +                                    zmem(1,indwq),zmem(1,imem2),
     +                                    zmem(1,7),nci,zmem(1,imem3),
     +                                    xija(icount),zmem(1,8),
     +                                    zmem(1,imem4),naind(index),
     +                                    zmem(1,imem5),jc0)
                                       naind(index) = naind(index) + 1
 190                                continue
 200                             continue
 210                          continue
c
c***** now the primitive intermediate integrals are complete
c***** for one particular primitive pair of iprim,jprim.
c***** next get half-contracted integrals.
c***** note there are 4 calls. see  hgp paper.
c
                              if (rab.lt.1.0d-8) ninty = ninty + 1
                              ncr = ncr + 1
                              do 240 n = 1 , nconx
                                 imem1 = iconx(3,n)
                                 imem2 = iconx(4,n)
                                 memcur = 0
                  do 230 i = 1 , nvl(iconx(1,n))
                     do 220 k = 1 , nvl(iconx(2,n))
                        memcur = memcur + 1
                        imem1m = imem1 + memcur
                        imem2m = imem2 + memcur
                        if (ncr.eq.1) then
                          if (ext4) then
         call vclr(zmem(jc0,imem2m),1,jc0+itprim)
                          end if
         call dcopy(jc0,zmem(1,imem1m),1,zmem(1,imem2m),1)
                        else
         call vadd(zmem(1,imem2m),1
     +            ,zmem(1,imem1m),1
     +            ,zmem(1,imem2m),1
     +            ,jc0)
                        end if
 220                 continue
 230                             continue
 240                          continue
c
                              do 270 n = 1 , nconc
                                 imem1 = iconc(3,n)
                                 imem2 = iconc(4,n)
                                 memcur = 0
                                 do 260 i = 1 , nvl(iconc(1,n))
                     do 250 k = 1 , nvl(iconc(2,n))
                        memcur = memcur + 1
                        imem1m = imem1 + memcur
                        imem2m = imem2 + memcur
                        if (ncr.eq.1) then
                          if (ext4) then
        call vclr(zmem(jc0,imem2m),1,jc0+itprim)
                          end if
        call dcopy(jc0,zmem(1,imem1m),1,zmem(1,imem2m),1)
                        else
         call vadd(zmem(1,imem2m),1
     +            ,zmem(1,imem1m),1
     +            ,zmem(1,imem2m),1
     +            ,jc0)
                        end if
 250                 continue
 260                             continue
 270                          continue
c
                              do 300 n = 1 , nconb
                                 imem1 = iconb(3,n)
                                 imem2 = iconb(4,n)
                                 memcur = 0
                     do 290 i = 1 , nvl(iconb(1,n))
                     do 280 k = 1 , nvl(iconb(2,n))
                        memcur = memcur + 1
                        imem2m = memcur + imem2
                        imem1m = memcur + imem1
                        if (ncr.eq.1) then
                          if (ext4) then
        call vclr(zmem(jc0,imem2m),1,jc0+itprim)
                          end if
        call dcopy(jc0,zmem(1,imem1m),1,zmem(1,imem2m),1)
                        else
         call vadd(zmem(1,imem2m),1
     +            ,zmem(1,imem1m),1
     +            ,zmem(1,imem2m),1
     +            ,jc0)
                        end if
 280                 continue
 290                             continue
 300                          continue
c
                              if (rab.ge.1.0d0-8) then
                                 do 350 n = 1 , ncona
                                    imem1 = icona(3,n)
                                    imem2 = icona(4,n)
                                    memcur = 0
                     do 340 i = 1 , nvl(icona(1,n))
                        do 330 k = 1 , nvl(icona(2,n))
                          memcur = memcur + 1
                          if (ncr.eq.1) then
                          if (ext4) then
                call vclr(zmem(jc0,imem2+memcur),1,jc0+itprim)
                          end if
                          do 310 jcount = 1 , jc0
                          zmem(jcount,imem2+memcur)
     +                       = alpha*zmem(jcount,
     +                       imem1+memcur)
 310                      continue
                          else
                          do 320 jcount = 1 , jc0
                          zmem(jcount,imem2+memcur)
     +                       = zmem(jcount,imem2+memcur)
     +                       + alpha*zmem(jcount,imem1+
     +                       memcur)
 320                      continue
                          end if
 330                    continue
 340                 continue
 350              continue
                              end if
c
 360                       continue
 370                    continue
                        if (ncr.ne.0) then
c
c**** multiply the half contracted c and d derivatives by the
c**** appropriate exponents
c
                           do 410 n = 1 , nconc
                              imem2 = iconc(4,n)
                              memcur = 0
                              do 400 i = 1 , nvl(iconc(1,n))
                                 do 390 k = 1 , nvl(iconc(2,n))
                                    memcur = memcur + 1
                                    do 380 jcount = 1 , jc0
                                       zmem(jcount,imem2+memcur)
     +                                    = zmem(jcount,imem2+memcur)
     +                                    *zmem(jcount,16)
 380                                continue
 390                             continue
 400                          continue
 410                       continue
c
                           do 420 jcount = 1 , jc0
                              zmem(jcount,16) = 2.0d0*zmem(jcount,7)
     +                           - zmem(jcount,16)
 420                       continue
c
                           do 460 n = 1 , nconb
                              imem2 = iconb(4,n)
                              memcur = 0
                              do 450 i = 1 , nvl(iconb(1,n))
                                 do 440 k = 1 , nvl(iconb(2,n))
                                    memcur = memcur + 1
                                    do 430 jcount = 1 , jc0
                                       zmem(jcount,imem2+memcur)
     +                                    = zmem(jcount,imem2+memcur)
     +                                    *zmem(jcount,16)
 430                                continue
 440                             continue
 450                          continue
 460                       continue
c
c**** here is a second tier integral test that operates on contracted
c**** derivative quantities rather than primitive quantities.
c
_IF1(cu)      call gather(jc0,igath(1),iscat2(jstart),igath(1))
_IFN1(cu)        call igthr(jc0,iscat2(jstart),igath(1),igath(1))
c
c**** form array that points to restricted list of contracted
c**** this array will go in igath
c
                           do 470 jcount = 1 , jc1
                              igath2(jcount) = 0
 470                       continue
                           do 480 jcount = 1 , jc0
                              igath2(igath(jcount)) = 1
 480                       continue
                           jc11 = 0
                           do 490 jcount = 1 , jc1
                              if (igath2(jcount).eq.1) then
                                 jc11 = jc11 + 1
                                 igath2(jcount) = jc11
                              end if
 490                       continue
                           do 500 jcount = 1 , jc0
                              igath(jcount) = igath2(igath(jcount))
 500                       continue
c
c*****  now that the loops over the i and j primitives for the one
c*****  particular ic,jc contracted pairing have been closed, we
c*****  can contract the kl part of the integrals.
                           call tst22(zmem,zmemc(1,memreq),igath(1),jc0,
     +                                jc11,memcon-memreq+1)
c
c***** form the gamma to multiply the derivative integrals
c
                           call drcdab(zmemc,ivlx,jvlx,kvlx,lvlx,ic,jc,
     +                                 kmin,kmax,lmin,lmax,l1,l2,ext2,
     +                                 dens,ipos,jpos,kpos,lpos,i1,j1,
     +                                 k1)
c
c***** now form ci-di in mtemp position
                           jc1 = 0
                           mtemp1 = mtemp + 1
                           mtemp2 = mtemp + 2
                           mtemp3 = mtemp + 3
                           do 520 kc = kmin , kmax
c*** a few logicals to determine extents of l loop
                              lmin1 = l1
                              lmax1 = l2
                              if (ext2) lmax1 = kc
                              if (kc.eq.kmax) lmax1 = lmax
                              if (kc.eq.kmin) lmin1 = lmin
                              do 510 lc = lmin1 , lmax1
                                 jc1 = jc1 + 1
                                 zmemc(jc1,mtemp1) = cobas(kc,1)
     +                              - cobas(lc,1)
                                 zmemc(jc1,mtemp2) = cobas(kc,2)
     +                              - cobas(lc,2)
                                 zmemc(jc1,mtemp3) = cobas(kc,3)
     +                              - cobas(lc,3)
 510                          continue
 520                       continue
c
c****form array that transforms complete list of contracted
c****to restricted list. put this array in igath.
c
                           jc11 = 0
                           do 530 jcount = 1 , jc1
                              if (igath2(jcount).ne.0) then
                                 jc11 = jc11 + 1
                                 igath(jc11) = jcount
                              end if
 530                       continue
c**** gather all the gammas and ci-di s  into restricted forms
c**** according to igath.
                           do 540 m = 1 , memreq - 9
c**** implicit in gather is igath(i).ge.i
            call dgthr(jc11,zmemc(1,m),zmemc(1,m),igath)
 540                       continue
c
c**** zero out temprorary gradient storage
                           do 560 icart = mtgrad + 1 , mtgrad + 9
                              do 550 jcount = 1 , jc1
                                 zmemc(jcount,icart) = 0.0d0
 550                          continue
 560                       continue
                           do 580 icart = 1 , 3
                              do 570 inat = 1 , na
                                 grtemp(inat,icart) = 0.0d0
 570                          continue
 580                       continue
c
                           ihzc = 0
                           if (rab.ge.1.0d-8) then
c*** find horizontal contributions to derivative integrals for a
                              call hzdrv(ix+1,jx,kx,lx,ix-1,jx,kx,lx,
     +                           ihzc,rab,zmemc,ab,jc11)
c*** put in contribution to gradient from derivative at position a
                              iatom = nwa(nup(ic))
                              mtgrad1 = mtgrad + 1
                              mtgrad2 = mtgrad + 2
                              mtgrad3 = mtgrad + 3
                              do 590 kl = 1 , jc11
                                 grtemp(iatom,1) = grtemp(iatom,1)
     +                              + zmemc(kl,mtgrad1)
                                 grtemp(iatom,2) = grtemp(iatom,2)
     +                              + zmemc(kl,mtgrad2)
                                 grtemp(iatom,3) = grtemp(iatom,3)
     +                              + zmemc(kl,mtgrad3)
 590                          continue
                           else
                              ihzc = ihzc + lx + 1
                              if (ix.gt.1) ihzc = ihzc + lx + 1
                           end if
c*** find horizontal contributions to derivative integrals for d
                           call hzdrv(ix,jx,kx,lx+1,ix,jx,kx,lx-1,ihzc,
     +                                rab,zmemc,ab,jc11)
c*** find horizontal contributions to derivative integrals for c
                           call hzdrv(ix,jx,kx+1,lx,ix,jx,kx-1,lx,ihzc,
     +                                rab,zmemc,ab,jc11)
                           jc1 = 0
                           jc11 = 0
                           do 610 kc = kmin , kmax
c*** a few logicals to determine extents of l loop
                              lmin1 = l1
                              lmax1 = l2
                              if (ext2) lmax1 = kc
                              if (kc.eq.kmax) lmax1 = lmax
                              if (kc.eq.kmin) lmin1 = lmin
                              katom = nwa(nup(kc))
                              do 600 lc = lmin1 , lmax1
                                 jc1 = jc1 + 1
                                 if (igath2(jc1).ne.0) then
                                    jc11 = jc11 + 1
                                    igath(jc11) = nwa(nup(lc))
                                    igath2(jc11) = katom
                                 end if
 600                          continue
 610                       continue
c**** put in contribution to gradient from derivative at position c
                           mtgrad1 = mtgrad + 7
                           mtgrad2 = mtgrad + 8
                           mtgrad3 = mtgrad + 9
                           do 620 kl = 1 , jc11
                              grtemp(igath2(kl),1)
     +                           = grtemp(igath2(kl),1)
     +                           + zmemc(kl,mtgrad1)
                              grtemp(igath2(kl),2)
     +                           = grtemp(igath2(kl),2)
     +                           + zmemc(kl,mtgrad2)
                              grtemp(igath2(kl),3)
     +                           = grtemp(igath2(kl),3)
     +                           + zmemc(kl,mtgrad3)
 620                       continue
c**** put in contribution to gradient from derivative at position d
                           mtgrad1 = mtgrad + 4
                           mtgrad2 = mtgrad + 5
                           mtgrad3 = mtgrad + 6
                           do 630 kl = 1 , jc11
                              grtemp(igath(kl),1) = grtemp(igath(kl),1)
     +                           + zmemc(kl,mtgrad1)
                              grtemp(igath(kl),2) = grtemp(igath(kl),2)
     +                           + zmemc(kl,mtgrad2)
                              grtemp(igath(kl),3) = grtemp(igath(kl),3)
     +                           + zmemc(kl,mtgrad3)
 630                       continue
c
c**** now infer the gradient at position b from the contents
c**** of grtemp. first zero out any contribution to be in grtemp
c
                           jatom = nwa(nup(jc))
                           grtemp(jatom,1) = 0.0d0
                           grtemp(jatom,2) = 0.0d0
                           grtemp(jatom,3) = 0.0d0
                           do 650 icart = 1 , 3
                              scalar = 0.0d0
                              do 640 inat = 1 , na
                                 grad(inat,icart) = grad(inat,icart)
     +                              + grtemp(inat,icart)
                                 scalar = scalar + grtemp(inat,icart)
 640                          continue
                              grad(jatom,icart) = grad(jatom,icart)
     +                           - scalar
 650                       continue
                        end if
                     end if
                  end if
c
 660           continue
            end if
 670        icorig = icorig + itprim
            zjtest = zjtest + itprim
 680     continue
 690  continue
_IF1(c)cdir$  nolist
_IF1(c)cdir$  vector
      return
      end
_IF(hpux11)
c HP compiler bug JAGae55357
c$HP$ OPTIMIZE LEVEL1
_ENDIF
      subroutine body(iscat2,zmem,zmemc,xijk,xija,xijt,xijp1,xijp2
     &,xijp3,nij,xklk,xkla,xklt,xklp1,xklp2,xklp3,igath,nkl
     &,fock,dens,codens
     &,ext1,ext2,ext3,i1,i2,j1,j2,k1,l1,l2,prefac,xinner)
c
c
c
c******** routine calculates two electron integrals
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto),zmemc(icsize,memcto)
     &,xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,xkla(nkl),xklk(nkl),xklp1(nkl),xklp2(nkl),xklp3(nkl)
     &,ab(3),naind(3),ncind(3),pa(3),igath(nkl),xijt(nij),xklt(nkl)
     &,m0(48),m1(48),m2(48),m3(48)
      dimension fock(nb,nb),codens(nb,nb),dens(nb,nb),iscat2(nkl)
      logical ofull,odiisn
_IF1()*     logical odo
      common/gen3/ofull(100),odiisn
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/gen/maxit,iter,concrit,conv,energy,thresh,dmax
c
      common/tempit/zitest,zipass,zjtest,zjpass
      common/scra  /iso(mxshel,48),nt
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/defunk/cobas(maxorb,3),nd
_IF(parallel)
c***   ***node-MPP***
INCLUDE(common/parallel)
      common/pardir/next,nbcnt(256)
c***   ***node-MPP***
_ENDIF
      logical ext1,ext2,ext3,ext4,ext5,ext7
c
      data factm /1.0d-02/
c
_IF1(c)cdir$  list
_IF1(c)cdir$  novector
      ivlx = nvl(ix)
      jvlx = nvl(jx)
      kvlx = nvl(kx)
      lvlx = nvl(lx)
      ipos = nposf(i1,ivlx,nd,np,ns)
      jpos = nposf(j1,jvlx,nd,np,ns)
      kpos = nposf(k1,kvlx,nd,np,ns)
      lpos = nposf(l1,lvlx,nd,np,ns)
      memreq = ivlx*jvlx*kvlx*lvlx
      jmax = j2
      icorig = 0
c
c**** the main contracted loop over i and j shells
c
      do 470 ic = i1 , i2
         if (ext1) jmax = ic
c**** symetry test on ishell
         do 30 it = 2 , nt
            id = iso(ic,it)
            if (id.gt.ic) then
               do 20 jc = j1 , jmax
                  ext5 = ext1 .and. ic.eq.jc
                  itprim = nc(ic)*nc(jc)
                  if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
                  icorig = icorig + itprim
                  zjtest = zjtest + itprim
 20            continue
               go to 470
            end if
            m0(it) = id
 30      continue
         do 460 jc = j1 , jmax
            ext5 = ext1 .and. ic.eq.jc
            itprim = nc(ic)*nc(jc)
            if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
c***** symmetry check on jshell
            if (ext1) then
               do 40 it = 2 , nt
                  ids = m0(it)
                  jds = iso(jc,it)
                  if (jds.gt.ic) go to 450
                  if (ids.lt.jds) then
                     nds = ids
                     ids = jds
                     jds = nds
                  end if
                  if (ids.eq.ic .and. jds.gt.jc) go to 450
                  m1(it) = ids
                  m2(it) = jds
 40            continue
            else
               do 50 it = 2 , nt
                  ids = m0(it)
                  jds = iso(jc,it)
                  if (ids.eq.ic .and. jds.gt.jc) go to 450
                  m1(it) = ids
                  m2(it) = jds
 50            continue
            end if
_IF(parallel)
            icount_dlb = icount_dlb + 1
            if(icount_dlb . eq. next) then
               mndp1 = ipg_nodeid() + 1
               nbcnt(mndp1) = nbcnt(mndp1) + 1
_ENDIF(parallel)
c
            dij = 4*codens(ic,jc)
            ab(1) = cobas(ic,1) - cobas(jc,1)
            ab(2) = cobas(ic,2) - cobas(jc,2)
            ab(3) = cobas(ic,3) - cobas(jc,3)
            rab = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
c           jpmax = nc(jc)
c
c*****  now the loop over the number of different kl
c*****  mini-batches in this particular maxi-batch
c
            do 440 nbrk = 1 , nbreak
               ext4 = .false.
               kmax = info(2,nbrk)
               lmax = info(4,nbrk)
               lmin = info(3,nbrk)
               kmin = info(1,nbrk)
               if (ext3) then
                  if (kmin.gt.ic) go to 440
                  if (kmax.ge.ic) then
                     if (kmin.eq.ic .and. lmin.gt.jc) go to 440
                     if (kmax.ne.ic .or. lmax.ge.jc) then
                        ext4 = .true.
                        kmax = ic
                        lmax = jc
                     end if
                  end if
               end if
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
               ncr = 0
               icount = icorig
               xijtm = 0.0d0
               do 60 ijprm = 1 , itprim
                  icount = icount + 1
                  xijtm = max(xijtm,xijt(icount))
 60            continue
c
               test1 = prefac/xijtm
               tout = dmax*tmax(nbrk)
               if (ofull(iter)) test1 = test1 * factm
_IF1()***** 
_IF1()*  THIS WORKS OK
_IF1()*              if (test1.lt.tout) then
_IF1()*****  TEST
_IF1()*              odo = test1.lt.tout .or. ofull(iter)
_IF1()*              if (odo) then
_IF1()*****  TEST
               if (test1.lt.tout) then
                  test2 = xinner/xijtm
                  tmax1 = tmax(nbrk)*10.0d0
                  jstart = info(5,nbrk)
c********* call to integral tests
                  if (nt.eq.1) then
                     call tst11(codens,dij,dmaxij,tmax1,test2,
     +                          xklt(jstart+1),igath,zmem(1,4),ic,jc,
     +                          iscat2(jstart+1),kmin,kmax,l1,l2,
     +                          lmax,lmin,nkl,ipsize,ext2,
     +                          jcount,jc0,jc1)
                  else
                     call tst1s(codens,dij,dmaxij,test2,
     +                          xklt(jstart+1),igath,zmem(1,4),ic,jc,
     +                          iscat2(jstart+1),kmin,kmax,m1,m2,
     +                          m3,l1,l2,lmax,lmin,nkl,ipsize,
     +                          ext2,ext3,jcount,jc0,jc1)
                  end if
                  if (jc0.ne.0) then
                     jstart = jstart + 1
c
c***** prepare all parameter arrays needed to form
c***** restricted list of intermediate primitive integrals
c
        call dgthr(jc0,xklp1(jstart),zmem(1,1),igath(1))
        call dgthr(jc0,xklp2(jstart),zmem(1,2),igath(1))
        call dgthr(jc0,xklp3(jstart),zmem(1,3),igath(1))
        call dgthr(jc0,xklk(jstart),zmem(1,9),igath(1))
                     if (nt.ne.1) then
                        do 70 kl = 1 , jc0
                           zmem(kl,9) = zmem(kl,9)*zmem(kl,7)
 70                     continue
                     end if
        call dgthr(jc0,xkla(jstart),zmem(1,7),igath(1))
                     do 80 kl = 1 , jc0
                        zmem(kl,4) = zmem(kl,1) - zmem(kl,4)
                        zmem(kl,5) = zmem(kl,2) - zmem(kl,5)
                        zmem(kl,6) = zmem(kl,3) - zmem(kl,6)
 80                  continue
                     icount = icorig
                     jcold = jc0
                     do 260 ijprm = 1 , itprim
                        icount = icount + 1
c
c***** the outer loop integral test
c
                        test1 = prefac/xijt(icount)
                        tout = dmaxij*tmax(nbrk)
                        if (test1.ge.tout) go to 260
                        if (nbrk.eq.1) zjpass = zjpass + 1
c
                        pa(1) = xijp1(icount) - cobas(ic,1)
                        pa(2) = xijp2(icount) - cobas(ic,2)
                        pa(3) = xijp3(icount) - cobas(ic,3)
c
c***** this next piece of code is not very important so don't
c***** worry if you don't understand it
c
                        ext7 = .false.
                        if (ext4) then
                           jc0 = jcold
                           jctemp = jcold - itprim
                           itxmax = itprim
                           if (jctemp.lt.0) then
                              jctemp = 0
                              itxmax = jc0
                           end if
_IF1(x)c$dir scalar
_IF1(t)cdir$ nextscalar
_IF1(c)ccdir$ nextscalar
                           do 90 itx = 1 , itxmax
                              jctemp = jctemp + 1
                              itx1 = igath(jctemp) + jstart - 1
                              if (itx1.ge.icount) then
                                 if (itx1.gt.icount) then
                                    jc0 = jctemp - 1
                                 else
                                    jc0 = jctemp
                                    ext7 = .true.
                                    zmem(jc0,9) = zmem(jc0,9)*0.5d0
                                 end if
                                 go to 100
                              end if
 90                        continue
                        end if
 100                    zipass = zipass + jc0
                        zitest = zitest + jcount
                        if (jc0.le.0) go to 265
c
c***** now do the loops for forming the intermediate
c***** primitive integrals
c
c**** call routine to form (ss|ss)m integrals + a few
c**** useful quantities.
c
c*** z(1-3) qi ::: z(4-6) qi-ci ::: z(7) c+d :::
c*** z(8) 1/(a+b+c+d) ::: z(9) xklk ::: z(10-12) wi-pi :::
c*** z(13-15) wi-qi ::: z(16-?) fm(t) :::
c
                        call sinter(zmem(1,1),zmem(1,2),zmem(1,3),
     +                              zmem(1,13),zmem(1,7),zmem(1,8),
     +                              zmem(1,9),zmem(1,10),zmem(1,11),
     +                              zmem(1,12),zmem(1,16),xijp1(icount),
     +                              xijp2(icount),xijp3(icount),
     +                              xija(icount),xijk(icount),mx,jc0,
     +                              ipsize,zmem(1,mx+16))
c
                        if (ext7) zmem(jc0,9) = zmem(jc0,9)*2.0d0
c
c*** loop over number of distinct integral classes
c*** needed to form the target integrals
c
                        memcur = 15 + mx
                        do 220 n = 1 , ninty
                           ityp = intgrl(1,n)
                           ktyp = intgrl(2,n)
                           ityp1 = ityp - 1
                           ktyp1 = ktyp - 1
                           ityp2 = ityp - 2
                           ktyp2 = ktyp - 2
                           ivl = nvl(ityp)
                           kvl = nvl(ktyp)
                           if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
                           if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
                           nrecur = intgrl(3,n)
c
                           mem1 = intgrl(4,n)
                           mem2 = intgrl(5,n)
                           mem3 = intgrl(6,n)
                           mem4 = intgrl(7,n)
                           mem5 = intgrl(8,n)
c
c now:
c loop over number of indices associated with i position
                           do 210 i = 1 , ivl
                              nr = nrecur
                              if (nr.lt.5) then
c***** reduction is at the i position
c
                                 index = ispdf1(i,1)
                                 indwp = index + 9
                                 natred = ispdf1(i,2)
                                 itemp = (natred-1)*kvl
                                 jmem1 = itemp + mem1
                                 jmem2 = itemp + mem2
                                 jmem5 = mem5 + (natred-1)*kvl1
                                 if (nr.ne.1 .and. nr.ne.3) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                                    nai = ixyz(index,natred,ityp2)
                                    if (nai.eq.0) then
                                       nr = nr - 1
                                       go to 110
                                    end if
                                    naind(2) = ixyz(2,natred,ityp2) + 1
                                    naind(3) = ixyz(3,natred,ityp2) + 1
                                    naind(index) = naind(index) - 1
                                    idoubr = iarray(naind(2),naind(3))
                                    itemp2 = (idoubr-1)*kvl
                                    jmem3 = itemp2 + mem3
                                    jmem4 = itemp2 + mem4
                                 end if
                              else
c
c**** reduction is at the k position
                                 itemp = (i-1)*kvl1
                                 jmem1 = itemp + mem1
                                 jmem2 = itemp + mem2
                                 itemp = (i-1)*kvl2
                                 jmem3 = itemp + mem3
                                 jmem4 = itemp + mem4
                                 if (nr.ne.5 .and. nr.ne.6) then
                                    naind(1) = ixyz(1,i,ityp1) + 1
                                    naind(2) = ixyz(2,i,ityp1) + 1
                                    naind(3) = ixyz(3,i,ityp1) + 1
                                 end if
                              end if
c
c now:
c loop over number of indices associated with k position
c
 110                          do 200 k = 1 , kvl
                                 memcur = memcur + 1
                                 if (nr.lt.5) then
                                    imem1 = jmem1 + k
                                    imem2 = jmem2 + k
                                 else
                                    index = ispdf1(k,1)
                                    indqc = index + 3
                                    indwq = index + 12
                                    natred = ispdf1(k,2)
                                    imem1 = jmem1 + natred
                                    imem2 = jmem2 + natred
                                 end if
                                 go to (120,130,140,150,160,170,180,190)
     +                                  , nr
cc*** question asked here is which of the vectorised vrr routines
c**** do we want
 120                             call recur1(zmem(1,memcur),pa(index),
     +                              zmem(1,imem1),zmem(1,indwp),
     +                              zmem(1,imem2),jc0)
                                 go to 200
 130                             imem3 = jmem3 + k
                                 imem4 = jmem4 + k
                                 call recur2(zmem(1,memcur),pa(index),
     +                              zmem(1,imem1),zmem(1,indwp),
     +                              zmem(1,imem2),xija(icount),nai,
     +                              zmem(1,imem3),zmem(1,7),zmem(1,8),
     +                              zmem(1,imem4),jc0)
                                 go to 200
 140                             nci = ixyz(index,k,ktyp1)
                                 if (nci.eq.0) go to 120
                                 naind(2) = ixyz(2,k,ktyp1) + 1
                                 naind(3) = ixyz(3,k,ktyp1) + 1
                                 naind(index) = naind(index) - 1
                                 imem5 = jmem5 +
     +                              iarray(naind(2),naind(3))
                                 call recur3(zmem(1,memcur),pa(index),
     +                              zmem(1,imem1),zmem(1,indwp),
     +                              zmem(1,imem2),zmem(1,8),nci,
     +                              zmem(1,imem5),jc0)
                                 go to 200
 150                             nci = ixyz(index,k,ktyp1)
                                 if (nci.eq.0) go to 130
                                 naind(2) = ixyz(2,k,ktyp1) + 1
                                 naind(3) = ixyz(3,k,ktyp1) + 1
                                 naind(index) = naind(index) - 1
                                 imem5 = jmem5 +
     +                              iarray(naind(2),naind(3))
                                 imem3 = jmem3 + k
                                 jmem4 = jmem4 + k
                                 call recur4(zmem(1,memcur),pa(index),
     +                              zmem(1,imem1),zmem(1,indwp),
     +                              zmem(1,imem2),xija(icount),nai,
     +                              zmem(1,imem3),zmem(1,7),zmem(1,8),
     +                              zmem(1,imem4),nci,zmem(1,imem5),jc0)
                                 go to 200
 160                             call recur5(zmem(1,memcur),
     +                              zmem(1,indqc),zmem(1,imem1),
     +                              zmem(1,indwq),zmem(1,imem2),jc0)
                                 go to 200
 170                             nci = ixyz(index,natred,ktyp2)
                                 if (nci.eq.0) go to 160
                                 ncind(2) = ixyz(2,natred,ktyp2) + 1
                                 ncind(3) = ixyz(3,natred,ktyp2) + 1
                                 ncind(index) = ncind(index) - 1
                                 imem3 = jmem3 +
     +                              iarray(ncind(2),ncind(3))
                                 imem4 = jmem4 +
     +                              iarray(ncind(2),ncind(3))
                                 call recur6(zmem(1,memcur),
     +                              zmem(1,indqc),zmem(1,imem1),
     +                              zmem(1,indwq),zmem(1,imem2),
     +                              zmem(1,7),nci,zmem(1,imem3),
     +                              xija(icount),zmem(1,8),zmem(1,imem4)
     +                              ,jc0)
                                 go to 200
 180                             if (naind(index).eq.1) go to 160
                                 naind(index) = naind(index) - 1
                                 imem5 = (iarray(naind(2),naind(3))-1)
     +                              *kvl1 + natred + mem5
                                 call recur7(zmem(1,memcur),
     +                              zmem(1,indqc),zmem(1,imem1),
     +                              zmem(1,indwq),zmem(1,imem2),
     +                              zmem(1,8),naind(index),zmem(1,imem5)
     +                              ,jc0)
                                 naind(index) = naind(index) + 1
                                 go to 200
 190                             if (naind(index).eq.1) go to 170
                                 nci = ixyz(index,natred,ktyp2)
                                 if (nci.eq.0) go to 180
                                 ncind(2) = ixyz(2,natred,ktyp2) + 1
                                 ncind(3) = ixyz(3,natred,ktyp2) + 1
                                 ncind(index) = ncind(index) - 1
                                 imem3 = jmem3 +
     +                              iarray(ncind(2),ncind(3))
                                 imem4 = jmem4 +
     +                              iarray(ncind(2),ncind(3))
                                 naind(index) = naind(index) - 1
                                 imem5 = (iarray(naind(2),naind(3))-1)
     +                              *kvl1 + natred + mem5
                                 call recur8(zmem(1,memcur),
     +                              zmem(1,indqc),zmem(1,imem1),
     +                              zmem(1,indwq),zmem(1,imem2),
     +                              zmem(1,7),nci,zmem(1,imem3),
     +                              xija(icount),zmem(1,8),zmem(1,imem4)
     +                              ,naind(index),zmem(1,imem5),jc0)
                                 naind(index) = naind(index) + 1
 200                          continue
 210                       continue
 220                    continue
c
c***** now the primitive intermediate integrals are complete
c***** for one particular primitive pair of iprim,jprim. they are
c***** contracted to yield half-contracted integrals.
c
                        ncr = ncr + 1
                        do 250 n = 1 , ncon
                           imem1 = intcon(3,n)
                           imem2 = intcon(4,n)
                           memcur = 0
            do 240 i = 1 , nvl(intcon(1,n))
               do 230 k = 1 , nvl(intcon(2,n))
                  memcur = memcur + 1
                  imem1m = memcur + imem1
                  imem2m = memcur + imem2
                  if (ncr.eq.1) then
                     if (ext4) then
          call vclr(zmem(jc0,imem2m),1,jc0+itprim)
                     end if
          call dcopy(jc0,zmem(1,imem1m),1,zmem(1,imem2m),1)
                  else
          call vadd(zmem(1,imem2m),1
     +             ,zmem(1,imem1m),1
     +             ,zmem(1,imem2m),1
     +             ,jc0)
                  end if
 230           continue
 240        continue
 250     continue
c
 260                 continue
 265                 continue
c
c*****  now that the loops over the i and j primitives for the one
c*****  particular ic,jc contracted pairing have been closed, we
c*****  can contract the kl part of the integrals.
c
                     if (jc0.gt.0.and.ncr.ne.0) then
_IF1(cu)      call gather(jc0,igath(1),iscat2(jstart),igath(1))
_IFN1(cu)       call igthr(jc0,iscat2(jstart),igath(1),igath(1))
                        call tst22(zmem,zmemc,igath(1),jc0,jc1,memcon)
c
c
c***** now we have all the contracted integrals. use the horizontal
c***** recursion relation (if necessary) to find the desired integral
c
                        ihzc = 0
                        if (jx.gt.1) then
c***** hrr must be applied on the ij position
                           jxtyp1 = jx - 1
                           ixtyp1 = ix - 1
                           do 330 kl = kx , klx
                              ihzc = ihzc + 1
                              mem1 = ihz(1,ihzc)
                              mem2 = ihz(2,ihzc)
                              mem3 = ihz(3,ihzc)
                              mem4 = ihz(4,ihzc)
                              mem5 = ihz(5,ihzc)
                              kvl = nvl(kl)
                              ijc = 0
                              do 320 i = 1 , ivlx
                                 naind(2) = ixyz(2,i,ixtyp1) + 1
                                 naind(3) = ixyz(3,i,ixtyp1) + 1
                                 do 310 j = 1 , jvlx
                                    ncind(2) = naind(2)
     +                                 + ixyz(2,j,jxtyp1)
                                    ncind(3) = naind(3)
     +                                 + ixyz(3,j,jxtyp1)
                                    ij = iarray(ncind(2),ncind(3))
                                    imem1 = ijc*kvl + mem1
                                    imem2 = (ij-1)*kvl + mem2
                                    ijc = ijc + 1
                                    if (rab.lt.1.0d-8) then
c********* when (ai-bi) is zero (a+b,0|c+d,0) = (a,b|c+d,0)
                                       do 270 k = 1 , kvl
         call dcopy(jc1,zmemc(1,imem2+k),1,zmemc(1,imem1+k),1)
 270                                   continue
                                       go to 310
                                    end if
                                    index1 = ispdf1(j,1)
                                    ncind(index1) = ncind(index1) - 1
                                    imem3a = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    if (jx.eq.2) then
*****  a p function is being added to the i/a position.
                                       do 280 k = 1 , kvl
                                         call hrz2(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),jc1)
 280                                   continue
                                       go to 310
                                    end if
                                    ired1 = ispdf1(j,2)
                                    index2 = ispdf1(ired1,1)
                                    ncind(index2) = ncind(index2) - 1
                                    imem4a = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index1) = ncind(index1) + 1
                                    imem3b = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    if (jx.eq.3) then
*****  a d function is being added to the i/a position
                                       do 290 k = 1 , kvl
                                         call hrz4(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),ab(index2)
     +                                      ,zmemc(1,imem3b+k),
     +                                      zmemc(1,imem4a+k),jc1)
 290                                   continue
                                       go to 310
                                    end if
                                    ired2 = ispdf1(ired1,2)
                                    index3 = ispdf1(ired2,1)
                                    ncind(index3) = ncind(index3) - 1
                                    imem4b = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index2) = ncind(index2) + 1
                                    imem3c = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem3
                                    ncind(index1) = ncind(index1) - 1
                                    imem4c = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem4
                                    ncind(index2) = ncind(index2) - 1
                                    imem5 = (iarray(ncind(2),ncind(3))
     +                                 -1)*kvl + mem5
                                    if (jx.eq.4) then
*****  an f function is being added to the i/a position
                                       do 300 k = 1 , kvl
                                         call hrz8(zmemc(1,imem1+k),
     +                                      zmemc(1,imem2+k),ab(index1),
     +                                      zmemc(1,imem3a+k),ab(index2)
     +                                      ,zmemc(1,imem3b+k),
     +                                      ab(index3),zmemc(1,imem3c+k)
     +                                      ,zmemc(1,imem4a+k),
     +                                      zmemc(1,imem4b+k),
     +                                      zmemc(1,imem4c+k),
     +                                      zmemc(1,imem5+k),jc1)
 300                                   continue
                                    end if
 310                             continue
 320                          continue
 330                       continue
c***** at this point the hrr on a and b was necessary and
c***** has been applied.
                        end if
                        if (lx.ne.1) then
c***** hrr is needed for c and d
c
c***** first form ci-di in mtemp position
                           jc1 = 0
                           mtemp1 = mtemp + 1
                           mtemp2 = mtemp + 2
                           do 350 kc = kmin , kmax
c*** a few logicals to determine extents of l loop
                              lmin1 = l1
                              lmax1 = l2
                              if (ext2) lmax1 = kc
                              if (kc.eq.kmax) lmax1 = lmax
                              if (kc.eq.kmin) lmin1 = lmin
                              do 340 lc = lmin1 , lmax1
                                 jc1 = jc1 + 1
                                 zmemc(jc1,mtemp) = cobas(kc,1)
     +                              - cobas(lc,1)
                                 zmemc(jc1,mtemp1) = cobas(kc,2)
     +                              - cobas(lc,2)
                                 zmemc(jc1,mtemp2) = cobas(kc,3)
     +                              - cobas(lc,3)
 340                          continue
 350                       continue
                           kxtyp1 = kx - 1
                           lxtyp1 = lx - 1
                           ihzc = ihzc + 1
                           klc = 0
                           mem1 = ihz(1,ihzc)
                           mem2 = ihz(2,ihzc)
                           mem3 = ihz(3,ihzc)
                           mem4 = ihz(4,ihzc)
                           mem5 = ihz(5,ihzc)
                           kl1 = kvlx*lvlx
                           kl2 = nvl(klx)
                           kl3 = nvl(klx-1)
                           kl4 = nvl(klx-2)
                           kl5 = nvl(klx-3)
                           do 430 k = 1 , kvlx
                              naind(2) = ixyz(2,k,kxtyp1) + 1
                              naind(3) = ixyz(3,k,kxtyp1) + 1
                              do 420 l = 1 , lvlx
                                 ncind(2) = ixyz(2,l,lxtyp1) + naind(2)
                                 ncind(3) = ixyz(3,l,lxtyp1) + naind(3)
                                 klc = klc + 1
                                 imem1 = mem1 + klc
                                 kl = iarray(ncind(2),ncind(3))
                                 imem2 = mem2 + kl
                                 index1 = ispdf1(l,1)
                                 indcd1 = mtemp - 1 + index1
                                 ncind(index1) = ncind(index1) - 1
                                 imem3a = mem3 +
     +                              iarray(ncind(2),ncind(3))
                                 ijc = -1
                                 if (lx.eq.2) then
*****  a p function is being added to the k/c position
                                    do 370 i = 1 , ivlx
                                       do 360 j = 1 , jvlx
                                         ijc = ijc + 1
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3 = imem3a + ijc*kl3
                                         call hrz2c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3),jc1)
 360                                   continue
 370                                continue
                                    go to 420
                                 end if
                                 ired1 = ispdf1(l,2)
                                 index2 = ispdf1(ired1,1)
                                 indcd2 = mtemp - 1 + index2
                                 ncind(index2) = ncind(index2) - 1
                                 imem4a = mem4 +
     +                              iarray(ncind(2),ncind(3))
                                 ncind(index1) = ncind(index1) + 1
                                 imem3b = mem3 +
     +                              iarray(ncind(2),ncind(3))
                                 if (lx.eq.3) then
c**** a d function is being added to the k/c position
                                    do 390 i = 1 , ivlx
                                       do 380 j = 1 , jvlx
                                         ijc = ijc + 1
                                         itemp = ijc*kl3
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3a = imem3a + itemp
                                         jmem3b = imem3b + itemp
                                         jmem4a = imem4a + ijc*kl4
                                         call hrz4c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3a),
     +                                      zmemc(1,indcd2),
     +                                      zmemc(1,jmem3b),
     +                                      zmemc(1,jmem4a),jc1)
 380                                   continue
 390                                continue
                                    go to 420
                                 end if
                                 ired2 = ispdf1(ired1,2)
                                 index3 = ispdf1(ired2,1)
                                 indcd3 = mtemp - 1 + index3
                                 ncind(index3) = ncind(index3) - 1
                                 imem4b = iarray(ncind(2),ncind(3))
     +                              + mem4
                                 ncind(index2) = ncind(index2) + 1
                                 imem3c = iarray(ncind(2),ncind(3))
     +                              + mem3
                                 ncind(index1) = ncind(index1) - 1
                                 imem4c = iarray(ncind(2),ncind(3))
     +                              + mem4
                                 ncind(index2) = ncind(index2) - 1
                                 imem5 = iarray(ncind(2),ncind(3))
     +                              + mem5
                                 if (lx.eq.4) then
c**** a d function is being added to the k/c position
                                    do 410 i = 1 , ivlx
                                       do 400 j = 1 , jvlx
                                         ijc = ijc + 1
                                         itemp = ijc*kl3
                                         itemp1 = ijc*kl4
                                         jmem1 = imem1 + ijc*kl1
                                         jmem2 = imem2 + ijc*kl2
                                         jmem3a = imem3a + itemp
                                         jmem3b = imem3b + itemp
                                         jmem3c = imem3c + itemp
                                         jmem4a = imem4a + itemp1
                                         jmem4b = imem4b + itemp1
                                         jmem4c = imem4c + itemp1
                                         jmem5 = imem5 + ijc*kl5
                                         call hrz8c(zmemc(1,jmem1),
     +                                      zmemc(1,jmem2),
     +                                      zmemc(1,indcd1),
     +                                      zmemc(1,jmem3a),
     +                                      zmemc(1,indcd2),
     +                                      zmemc(1,jmem3b),
     +                                      zmemc(1,indcd3),
     +                                      zmemc(1,jmem3c),
     +                                      zmemc(1,jmem4a),
     +                                      zmemc(1,jmem4b),
     +                                      zmemc(1,jmem4c),
     +                                      zmemc(1,jmem5),jc1)
 400                                   continue
 410                                continue
                                 end if
 420                          continue
 430                       continue
                        end if
c***** the hrr has now been applied and the integrals lie in
c*****    zmemc(1,memans)
c***** send the integrals to the fock build!
                        call fockb(zmemc(1,memans+1),ivlx,jvlx,kvlx,
     +                             lvlx,ic,jc,kmin,kmax,lmin,lmax,l1,l2,
     +                             ext2,fock,dens,zmemc(1,mtemp),memreq,
     +                             ipos,jpos,kpos,lpos,i1,j1,k1)
                     end if
                  end if
               end if
c
 440        continue
_IF(parallel)
         next = ipg_dlbtask()
         endif
_ENDIF
 450        icorig = icorig + itprim
            zjtest = zjtest + itprim
 460     continue
 470  continue
_IF1(c)cdir$  nolist
_IF1(c)cdir$  vector
      return
      end
_IF(hpux11)
c HP compiler bug JAGae55357
c$HP$ OPTIMIZE LEVEL3
_ENDIF
      subroutine tizer(iw,ifreem,ixz,jxz,kxz,lxz,ijxz,klxz,mxz)
c
      implicit REAL  (a-h,o-z)
_IF1()      character *2 buff1,buff2
      dimension iw(ijxz,klxz,mxz)
      common/iofile/iread,iwr
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
c
c**** routine decides for a given integral type which
c**** intermediate integrals are needed
c
      ix = ixz
      jx = jxz
      kx = kxz
      lx = lxz
      ijx = ijxz
      klx = klxz
      mx = mxz
      mtemp = 1
      do 40 m = 1 , mx
         do 30 ij = 1 , ijx
            do 20 kl = 1 , klx
               iw(ij,kl,m) = 0
 20         continue
 30      continue
 40   continue
      do 60 ij = ix , ijx
         do 50 kl = kx , klx
            iw(ij,kl,1) = -1
 50      continue
 60   continue
c
      maxm = 1
      ita = ijx + klx
 70   do 100 ij = 1 , ijx
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx
            if (ij.ge.itates) then
c
               do 90 kl = 1 , klx
                  ijkl = ij + kl
                  if (ijkl.eq.ita) then
c
                     do 80 m = 1 , maxm
                        if (iw(ij,kl,m).ne.0) then
                           if (kl.ne.1) then
                              call kcrs(iw,ijx,klx,mx,ij,kl,m)
                              go to 80
                           end if
                           call icrs(iw,ijx,klx,mx,ij,kl,m)
                        end if
 80                  continue
                  end if
 90            continue
            end if
         end if
 100  continue
c
      maxm = maxm + 1
      ita = ita - 1
      if (ita.gt.2) go to 70
c
c**** sort out memory requirements for the ssss(m) integrals
c
      mem = 15
      do 110 i = 1 , mx
         iw(1,1,i) = mem
         mem = mem + 1
 110  continue
      memche = mem + mx + 4
c
c**** find the memory locations of the child integrals
c**** and the other arrays needed to form an integral- also
c**** give a positon in a list to each integral
c
      ita = 3
      maxm = mx - 1
      ninty = 0
 120  do 160 ij = 1 , ijx
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx
            if (ij.ge.itates) then
c
               do 150 kl = 1 , klx
                  ijkl = kl + ij
                  if (ijkl.eq.ita) then
                     nadd = nvl(ij)*nvl(kl)
c
                     do 140 m = 1 , maxm
                        if (iw(ij,kl,m).eq.0) go to 140
                        ninty = ninty + 1
                        intgrl(1,ninty) = ij
                        intgrl(2,ninty) = kl
                        if (kl.ne.1) then
                           call kcrs2(iw,ijx,klx,mx,ij,kl,m,
     +                                intgrl(3,ninty))
                           go to 130
                        end if
                        call icrs2(iw,ijx,klx,mx,ij,kl,m,intgrl(3,ninty)
     +                             )
 130                    iw(ij,kl,m) = mem
                        mem = mem + nadd
 140                 continue
                  end if
 150           continue
            end if
         end if
 160  continue
      maxm = maxm - 1
      ita = ita + 1
      if (maxm.gt.0) go to 120
c
c***** now prepare memory for contractions part
c
      if (mem.lt.memche) mem = memche
      mempri = mem
c     mp1 = mem + 1
      memcon = 0
      ncon = 0
      do 180 ij = ix , ijx
         do 170 kl = kx , klx
            ncon = ncon + 1
            intcon(1,ncon) = ij
            intcon(2,ncon) = kl
            intcon(3,ncon) = iw(ij,kl,1)
            intcon(4,ncon) = mem
            intcon(5,ncon) = memcon
            iw(ij,kl,1) = memcon
            nadd = nvl(ij)*nvl(kl)
            memcon = memcon + nadd
            mem = mem + nadd
 170     continue
 180  continue
      mempto = mem
      ipsize = ifreem/(mempto+1)
      ictot1 = (mempri-1)*ipsize
      icsiz1 = ictot1/(memcon+1)
c
      jcount = 0
      nadd1 = nvl(ix)*nvl(jx)
      if (jx.le.1) then
c******** here hrr not applied to ij
         if (lx.eq.1) then
c******* no hrr applied at all
            memans = 0
            icsize = icsiz1
            go to 240
         end if
c******* hrr applied to kl only
         nadd = nvl(kx)*nvl(lx)*nadd1
         memans = memcon
         mtemp = memans + nadd + 1
c      if (memtes.gt.mempto) call caserr
c     &('bug3 in memory allocation')
         icsiz2 = ifreem/(mtemp+6)
         icsize = min(icsiz1,icsiz2)
         ihz(1,1) = memans
         do 190 kl = kx , klx
            indk = (klx-kl) + 2
            ihz(indk,1) = iw(ix,kl,1)
 190     continue
c********* then hrr must be used on ij part of integral
c********* and info needed for this is found here
      else if (.not.lx.eq.1) then
c******** here hrr must be applied to kl also
         memans = 0
         mem = nadd1*nvl(kx)*nvl(lx)
         if (mem.lt.memcon) call caserr('bug1 in cider')
         do 210 kl = kx , klx
            jcount = jcount + 1
            ihz(1,jcount) = mem
            mem = mem + nadd1*nvl(kl)
c
            do 200 ij = ix , ijx
               indi = (ijx-ij) + 2
               ihz(indi,jcount) = iw(ij,kl,1)
 200        continue
 210     continue
c
         mtemp = mem + 1
         icsiz2 = ifreem/(mem+7)
         icsize = min(icsiz1,icsiz2)
         jc1 = jcount + 1
         ihz(1,jc1) = 0
         do 220 jc = 1 , jcount
            indj = jcount - jc + 2
            ihz(indj,jc1) = ihz(1,jc)
 220     continue

      else
c******** no hrr is needed for kl part
         kl = klx
         ihz(1,1) = memcon
         do 230 ij = ix , ijx
            indi = (ijx-ij) + 2
            ihz(indi,1) = iw(ij,kl,1)
 230     continue
         memans = ihz(1,1)
         memtes = memcon + nadd1*nvl(klx) + 1
         icsiz2 = ifreem/memtes
         icsize = min(icsiz1,icsiz2)
      end if
c
_IF1()c     write(iwr,*) 'ninty is ',ninty
_IF1()c      write(iwr,104)
_IF1()c104   format(1x,'integral','  value of m','    start pt in memory')
 240  continue
_IF1()         do 270 ij = 1 , ijx
_IF1()         if (ij.eq.1) buff1 = 'ss'
_IF1()         if (ij.eq.2) buff1 = 'ps'
_IF1()         if (ij.eq.3) buff1 = 'ds'
_IF1()         if (ij.eq.4) buff1 = 'fs'
_IF1()         if (ij.eq.5) buff1 = 'gs'
_IF1()         do 260 kl = 1 , klx
_IF1()            if (kl.eq.1) buff2 = 'ss'
_IF1()            if (kl.eq.2) buff2 = 'ps'
_IF1()            if (kl.eq.3) buff2 = 'ds'
_IF1()            if (kl.eq.4) buff2 = 'fs'
_IF1()            if (kl.eq.5) buff2 = 'gs'
_IF1()            do 250 m = 1 , mx
_IF1()               if (iw(ij,kl,m).ne.0) then
_IF1()                  m1 = m - 1
_IF1()               end if
_IF1()c      write(iwr,105) buff1,buff2,m1,iw(ij,kl,m)
_IF1()c105   format(1x,'(',a2,'|',a2,')',6x,i5,5x,i10)
_IF1() 250        continue
_IF1() 260     continue
_IF1() 270  continue
_IF1()c      write(iwr,134)
_IF1()c134   format(1x,'integral',' recur eq.',15x,'parents position
_IF1()c     & in memory')
_IF1()      do 280 i = 1 , ninty
_IF1()         if (intgrl(1,i).eq.1) buff1 = 'ss'
_IF1()         if (intgrl(1,i).eq.2) buff1 = 'ps'
_IF1()         if (intgrl(1,i).eq.3) buff1 = 'ds'
_IF1()         if (intgrl(1,i).eq.4) buff1 = 'fs'
_IF1()         if (intgrl(1,i).eq.5) buff1 = 'gs'
_IF1()         if (intgrl(2,i).eq.1) buff2 = 'ss'
_IF1()         if (intgrl(2,i).eq.2) buff2 = 'ps'
_IF1()         if (intgrl(2,i).eq.3) buff2 = 'ds'
_IF1()         if (intgrl(2,i).eq.4) buff2 = 'fs'
_IF1()         if (intgrl(2,i).eq.5) buff2 = 'gs'
_IF1()c      write(iwr,135) buff1,buff2,(intgrl(j,i),j=3,8)
_IF1()c135   format(1x,'(',a2,'|',a2,')',2x,i5,7x,5i10)
_IF1() 280  continue
c
_IF1()c      write(iwr,*) 'ipsize is',ipsize
_IF1()c      write(iwr,*) 'icsize is',icsize
      memcto = ifreem/icsize
      return
      end
      subroutine tizex(iw,ifreem,ixz,jxz,kxz,lxz,ijxz,klxz,mxz)
c
      implicit REAL  (a-h,o-z)
_IF1()      character *2 buff1,buff2
      dimension iw(ijxz,klxz,mxz)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      common/cidex/intgrl(8,100),ihex(8)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize
c
c**** routine decides for a given exchange integral type which
c**** intermediate integrals are needed
c
      ix = ixz
      jx = jxz
      kx = kxz
      lx = lxz
      ijx = ijxz
      klx = klxz
      mx = mxz
      mtemp = 1
      do 40 m = 1 , mx
         do 30 ij = 1 , ijx
            do 20 kl = 1 , klx
               iw(ij,kl,m) = 0
 20         continue
 30      continue
 40   continue
      do 60 ij = ix , ijx
         do 50 kl = kx , klx
            iw(ij,kl,1) = -1
 50      continue
 60   continue
c
      maxm = 1
      ita = ijx + klx
 70   do 100 ij = 1 , ijx
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx
            if (ij.ge.itates) then
c
               do 90 kl = 1 , klx
                  ijkl = ij + kl
                  if (ijkl.eq.ita) then
c
                     do 80 m = 1 , maxm
                        if (iw(ij,kl,m).ne.0) then
                           if (kl.ne.1) then
                              call kcrs(iw,ijx,klx,mx,ij,kl,m)
                              go to 80
                           end if
                           call icrs(iw,ijx,klx,mx,ij,kl,m)
                        end if
 80                  continue
                  end if
 90            continue
            end if
         end if
 100  continue
c
      maxm = maxm + 1
      ita = ita - 1
      if (ita.gt.2) go to 70
c
c**** sort out memory requirements for the ssss(m) integrals
c
      mem = 11
      do 110 i = 1 , mx
         iw(1,1,i) = mem
         mem = mem + 1
 110  continue
      memche = mem + mx + 4
c
c**** find the memory locations of the child integrals
c**** and the other arrays needed to form an integral- also
c**** give a positon in a list to each integral
c
      ita = 3
      maxm = mx - 1
      ninty = 0
 120  do 160 ij = 1 , ijx
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx
            if (ij.ge.itates) then
c
               do 150 kl = 1 , klx
                  ijkl = kl + ij
                  if (ijkl.eq.ita) then
                     nadd = nvl(ij)*nvl(kl)
c
                     do 140 m = 1 , maxm
                        if (iw(ij,kl,m).eq.0) go to 140
                        ninty = ninty + 1
                        intgrl(1,ninty) = ij
                        intgrl(2,ninty) = kl
                        if (kl.ne.1) then
                           call kcrs2(iw,ijx,klx,mx,ij,kl,m,
     +                                intgrl(3,ninty))
                           go to 130
                        end if
                        call icrs2(iw,ijx,klx,mx,ij,kl,m,intgrl(3,ninty)
     +                             )
 130                    iw(ij,kl,m) = mem
                        mem = mem + nadd
 140                 continue
                  end if
 150           continue
            end if
         end if
 160  continue
      maxm = maxm - 1
      ita = ita + 1
      if (maxm.gt.0) go to 120
c
c***** now prepare memory for contractions part
c
      if (jx.eq.1) then
c**** (ds|ds)
         ihex(1) = iw(3,3,1)
      end if
      if (jx.eq.2) then
c**** (dp|dp)
         ihex(1) = iw(4,4,1)
         ihex(2) = iw(4,3,1)
         ihex(3) = iw(3,3,1)
      end if
      if (jx.eq.3) then
c**** (dd|dd)
         ihex(1) = iw(5,5,1)
         ihex(2) = iw(5,4,1)
         ihex(3) = iw(5,3,1)
         ihex(4) = iw(4,4,1)
         ihex(5) = iw(4,3,1)
         ihex(6) = iw(3,3,1)
      end if
      if (mem.lt.memche) mem = memche
      mempri = mem
      mempto = mem
      ipsize = ifreem/(mempto+1)
c
_IF1()c999   continue
_IF1()c     write(iwr,*) 'ninty is ',ninty
_IF1()c      write(iwr,104)
_IF1()c104   format(1x,'integral','  value of m','    start pt in memory')
_IF1()      do 190 ij = 1 , ijx
_IF1()         if (ij.eq.1) buff1 = 'ss'
_IF1()         if (ij.eq.2) buff1 = 'ps'
_IF1()         if (ij.eq.3) buff1 = 'ds'
_IF1()         if (ij.eq.4) buff1 = 'fs'
_IF1()         if (ij.eq.5) buff1 = 'gs'
_IF1()         do 180 kl = 1 , klx
_IF1()            if (kl.eq.1) buff2 = 'ss'
_IF1()            if (kl.eq.2) buff2 = 'ps'
_IF1()            if (kl.eq.3) buff2 = 'ds'
_IF1()            if (kl.eq.4) buff2 = 'fs'
_IF1()            if (kl.eq.5) buff2 = 'gs'
_IF1()            do 170 m = 1 , mx
_IF1()               if (iw(ij,kl,m).ne.0) then
_IF1()                  m1 = m - 1
_IF1()               end if
_IF1()c      write(iwr,105) buff1,buff2,m1,iw(ij,kl,m)
_IF1()c105   format(1x,'(',a2,'|',a2,')',6x,i5,5x,i10)
_IF1() 170        continue
_IF1() 180     continue
_IF1() 190  continue
_IF1()c      write(iwr,134)
_IF1()c134   format(1x,'integral',' recur eq.',15x,'parents position
_IF1()c     & in memory')
_IF1()      do 200 i = 1 , ninty
_IF1()         if (intgrl(1,i).eq.1) buff1 = 'ss'
_IF1()         if (intgrl(1,i).eq.2) buff1 = 'ps'
_IF1()         if (intgrl(1,i).eq.3) buff1 = 'ds'
_IF1()         if (intgrl(1,i).eq.4) buff1 = 'fs'
_IF1()         if (intgrl(1,i).eq.5) buff1 = 'gs'
_IF1()         if (intgrl(2,i).eq.1) buff2 = 'ss'
_IF1()         if (intgrl(2,i).eq.2) buff2 = 'ps'
_IF1()         if (intgrl(2,i).eq.3) buff2 = 'ds'
_IF1()         if (intgrl(2,i).eq.4) buff2 = 'fs'
_IF1()         if (intgrl(2,i).eq.5) buff2 = 'gs'
_IF1()c      write(iwr,135) buff1,buff2,(intgrl(j,i),j=3,8)
_IF1()c135   format(1x,'(',a2,'|',a2,')',2x,i5,7x,5i10)
_IF1() 200  continue
c
c      write(iwr,*) 'ipsize is',ipsize
c      write(iwr,*) 'icsize is',icsize
      return
      end
      subroutine cidonb(iw,iw2,ifreem,ixz,jxz,ijxz,mxz)
c
      implicit REAL  (a-h,o-z)
_IF1()      character *2 buff1,buffi
      dimension iw(ijxz+1,mxz+1),iw2(ijxz+1,4)
_IF1()      dimension buffi(6)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,ix,jx,kxun,lxun,mx,ijx,klxun,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mtgrad
c
c**** routine decides for a given derivative one electron integral type
c**** intermediate integrals are needed. (nuclear attraction 1-e integra
c
      do 40 i = 1 , 5
         do 20 j = 1 , 50
            ihz(i,j) = -1
 20      continue
         do 30 j = 1 , 20
            icona(i,j) = -1
            iconb(i,j) = -1
            iconc(i,j) = -1
            iconx(i,j) = -1
 30      continue
 40   continue
      ijxz1 = ijxz + 1
_IF1()      buffi(1) = 's'
_IF1()      buffi(2) = 'p'
_IF1()      buffi(3) = 'd'
_IF1()      buffi(4) = 'f'
_IF1()      buffi(5) = 'g'
_IF1()      buffi(6) = 'h'
      ix = ixz
      jx = jxz
      ijx = ijxz
      mx = mxz
      do 60 m = 1 , mx + 1
         do 50 ij = 1 , ijx + 1
            iw(ij,m) = 0
 50      continue
 60   continue
      do 80 m = 1 , 4
         do 70 ij = 1 , ijx + 1
            iw2(ij,m) = 0
 70      continue
 80   continue
c
c**** the four integral classes need are
c
c**(ix+1 | jx )
      do 90 ij = ix + 1 , ijx + 1
         iw2(ij,1) = -1
         iw(ij,1) = -1
 90   continue
c**(ix | kx+1)
      do 100 ij = ix , ijx + 1
         iw2(ij,2) = -1
         iw(ij,1) = -1
 100  continue
c**(ix-1 | jx )
      if (ix.gt.1) then
         do 110 ij = ix - 1 , ijx - 1
            iw2(ij,3) = -1
            iw(ij,1) = -1
 110     continue
      end if
c**(ix | jx-1 )
      if (jx.gt.1) then
         do 120 ij = ix , ijx - 1
            iw2(ij,4) = -1
            iw(ij,1) = -1
 120     continue
      end if
c
c**** find out all the intermediate integrals needed
c
      maxm = 1
      ita = ijx + 1
 130  do 150 ij = 1 , ijx + 1
c
         if (ij.eq.ita) then
c
            do 140 m = 1 , maxm
               if (iw(ij,m).ne.0) then
                  iw(ij-1,m) = -1
                  iw(ij-1,m+1) = -1
                  if (ij.gt.2) then
                     iw(ij-2,m) = -1
                     iw(ij-2,m+1) = -1
                  end if
               end if
 140        continue
         end if
 150  continue
c
      maxm = maxm + 1
      ita = ita - 1
      if (ita.gt.1) go to 130
c
c**** sort out memory requirements for the (s|a|s)(m) integrals
c
      mem = 7
      do 160 i = 1 , mx + 1
         iw(1,i) = mem
         mem = mem + 1
 160  continue
c
c**** find the memory locations of the child integrals
c**** and the other arrays needed to form an integral- also
c**** give a positon in a list to each integral
c
      ita = 2
      maxm = mx
      ninty = 0
 170  do 200 ij = 1 , ijx + 1
         if (ij.eq.ita) then
c
            do 190 m = 1 , maxm
               if (iw(ij,m).ne.0) then
                  ninty = ninty + 1
                  do 180 mq = 1 , 8
                     intgrl(mq,ninty) = -1
 180              continue
                  intgrl(1,ninty) = ij
                  intgrl(2,ninty) = iw(ij-1,m)
                  intgrl(3,ninty) = iw(ij-1,m+1)
                  if (ij.gt.2) then
                     intgrl(4,ninty) = iw(ij-2,m)
                     intgrl(5,ninty) = iw(ij-2,m+1)
                  end if
                  iw(ij,m) = mem
                  mem = mem + nvl(ij)
               end if
 190        continue
         end if
 200  continue
      maxm = maxm - 1
      ita = ita + 1
      if (maxm.gt.0) go to 170
_IF1()c      write(iwr,*) ' ninty is ',ninty
_IF1()c      write(iwr,104)
_IF1()c104   format(1x,'integral','  value of m','    start pt in memory')
_IF1()      do 220 ij = 1 , ijx + 1
_IF1()         buff1 = buffi(ij)
_IF1()         do 210 m = 1 , mx + 1
_IF1()            if (iw(ij,m).ne.0) then
_IF1()               m1 = m - 1
_IF1()            end if
_IF1()c      write(iwr,105) buff1,m1,iw(ij,m)
_IF1()c105   format(1x,'(',a2,'| s )',6x,i5,5x,i10)
_IF1() 210     continue
_IF1() 220  continue
_IF1()c      write(iwr,134)
_IF1()c134   format(1x,'integral',' recur eq.',15x,'parents position
_IF1()c     & in memory')
_IF1()      do 230 i = 1 , ninty
_IF1()         buff1 = buffi(intgrl(1,i))
_IF1()c      write(iwr,135) buff1,(intgrl(j,i),j=2,6)
_IF1()c135   format(1x,'(',a2,'| s )',2x,i5,7x,5i10)
_IF1() 230  continue
c
c**** workout the contraction and hrr for all six integral classes
c
c** start with 3 reduced integral classes
      nconx = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1708)
_IF1()c1708  format(1x)
      if (ix.gt.1) then
_IF1()c      write(iwr,1706) buffi(ix-1),buffi(jx)
_IF1()c1706  format(1x,'the target class is ','(',a1,'|a|',a1,')')
      end if
      if (jx.gt.1) then
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx-1)
      end if
_IF1()c      write(iwr,1704)
      do 240 ij = 1 , ijx
         if (iw2(ij,3).ne.0 .or. iw2(ij,4).ne.0) then
            nconx = nconx + 1
            iconx(1,nconx) = ij
            iconx(3,nconx) = iw(ij,1)
            iconx(4,nconx) = mem
            if (iw2(ij,3).ne.0) iw2(ij,3) = mem
            if (iw2(ij,4).ne.0) iw2(ij,4) = mem
_IF1()            buff1 = buffi(ij)
_IF1()c      write(iwr,1705) buff1,iconx(3,nconx),iconx(4,nconx)
            mem = mem + nvl(ij)
         end if
 240  continue
c
c**   (ix | a | jx+1) integral classes
c
      nconb = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx+1)
_IF1()c      write(iwr,1704)
      do 250 ij = ix , ijx + 1
         nconb = nconb + 1
         iconb(1,nconb) = ij
         iconb(3,nconb) = iw(ij,1)
         iconb(4,nconb) = mem
         iw2(ij,2) = mem
         mem = mem + nvl(ij)
_IF1()         buff1 = buffi(ij)
_IF1()c      write(iwr,1705) buff1,iconb(3,nconb),iconb(4,nconb)
 250  continue
c
c**   (ix+1 | jx ) integral classes
_IF1()c      write(iwr,1708)
_IF1()cc      write(iwr,1706) buffi(ix+1),buffi(jx)
_IF1()c      write(iwr,1704)
_IF1()c1704   format(1x,'integral',' uncontracted position',
_IF1()c     &' half-contracted position',
_IF1()c     &' contracted position')
      ncona = 0
      do 260 ij = ix + 1 , ijx + 1
         ncona = ncona + 1
         icona(1,ncona) = ij
         icona(3,ncona) = iw(ij,1)
         icona(4,ncona) = mem
         iw2(ij,1) = mem
         mem = mem + nvl(ij)
_IF1()         buff1 = buffi(ij)
_IF1()c      write(iwr,1705) buff1,icona(3,ncona),icona(4,ncona)
_IF1()c1705  format(1x,'(',a2,'| a | s )',6x,i5,25x,20x,i5)
 260  continue
c
      ncon = ncona + nconb + nconx
      mem = mem + nvl(ix)*nvl(jx)
      memold = mem
c
c**** terms giving d/da(ix jx |kx lx)
      jcount = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,*) 'the horizontal recursion relations'
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix+1),buffi(jx)
_IF1()c      if (ix.gt.1)
_IF1()c     &write(iwr,1706) buffi(ix-1),buffi(jx)
_IF1()c      write(iwr,1708)
      call hzcdrd(ix+1,jx,ix-1,jx,mem,jcount,iw2(1,1),iw2(1,3),ijxz1)
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx+1)
      memcmx = mem
      mem = memold
_IF1()c      if (jx.gt.1)
_IF1()c     &write(iwr,1706) buffi(ix),buffi(jx-1)
_IF1()c      write(iwr,1708)
      call hzcdrd(ix,jx+1,ix,jx-1,mem,jcount,iw2(1,2),iw2(1,4),ijxz1)
      memcmx = max(mem,memcmx)
_IF1()c      write(iwr,*) 'memcmx is ',memcmx
c
      mempto = memcmx + 1
      ipsize = ifreem/(mempto+1)
c
      return
      end
      subroutine cidone(iw,ifreem,ixz,kxz,ninty1,kinup,koveru,
     +                  kindow,koverd)
c
      implicit REAL  (a-h,o-z)
_IF1()       character *2 buff1,buff2
      dimension iw(ixz+1,kxz,2)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ixun,jxun,kxun,lxun,mxun,ijxun,klxun,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
c
c**** routine deals with the easier one electron derivative integrals
c**** routine decides for a given integral type which
c**** intermediate integrals are needed
c
      ixun = ixz + 1
      kxun = kxz
      ix = ixz + 1
      kx = kxz
      do 40 m = 1 , 2
         do 30 ij = 1 , ix
            do 20 kl = 1 , kx
               iw(ij,kl,m) = 0
 20         continue
 30      continue
 40   continue
      iw(ix,kx,1) = -1
      if (ix.gt.2) iw(ix-2,kx,1) = -1
      iw(ix,kx,2) = -1
      if (ix.gt.2) iw(ix-2,kx,2) = -1
c
c**** mark all necessary intermediates to form the kinetic
c**** integrals
c
c     maxm = 1
      ita = ix + kx
 50   do 70 ij = 1 , ix
         if (ij.lt.ita) then
            if ((ij+kx).ge.ita) then
c
               do 60 kl = 1 , kx
                  ijkl = ij + kl
                  if (ijkl.eq.ita) then
                     if (iw(ij,kl,2).ne.0) then
                        if (kl.ne.1) then
                           iw(ij,kl-1,2) = -1
                           iw(ij,kl,1) = -1
                           if (kl.gt.2) then
                              iw(ij,kl-2,2) = -1
                              iw(ij,kl-2,1) = -1
                           end if
                           if (ij.gt.1) iw(ij-1,kl-1,2) = -1
                        else
                           iw(ij-1,kl,2) = -1
                           iw(ij,kl,1) = -1
                           if (ij.gt.2) then
                              iw(ij-2,kl,2) = -1
                              iw(ij-2,kl,1) = -1
                           end if
                           if (kl.gt.1) iw(ij-1,kl-1,2) = -1
                        end if
                     end if
                  end if
 60            continue
            end if
         end if
 70   continue
c
      ita = ita - 1
      if (ita.gt.2) go to 50
c
c**** mark all necessary intermediates to form the overlap
c**** integrals
c
c     maxm = 1
      ita = ix + kx
 80   do 100 ij = 1 , ix
         if (ij.lt.ita) then
            if ((ij+kx).ge.ita) then
c
               do 90 kl = 1 , kx
                  ijkl = ij + kl
                  if (ijkl.eq.ita) then
                     if (iw(ij,kl,1).ne.0) then
                        if (kl.ne.1) then
                           iw(ij,kl-1,1) = -1
                           if (kl.gt.2) iw(ij,kl-2,1) = -1
                           if (ij.gt.1) iw(ij-1,kl-1,1) = -1
                        else
                           iw(ij-1,kl,1) = -1
                           if (ij.gt.2) iw(ij-2,kl,1) = -1
                           if (kl.gt.1) iw(ij-1,kl-1,1) = -1
                        end if
                     end if
                  end if
 90            continue
            end if
         end if
 100  continue
c
      ita = ita - 1
      if (ita.gt.2) go to 80
c
c**** sort out memory requirements for the (s|s) and (s|t|s) integrals
c
      iw(1,1,1) = 0
      iw(1,1,2) = 1
      mem = 10
c
c**** find the memory locations of the child overlap integrals
c**** and the other arrays needed to form the integral- also
c**** give a positon in a list to each integral
c
      ita = 3
      ninty = 0
 110  do 140 ij = 1 , ix
         if (ij.lt.ita) then
            if ((ij+kx).ge.ita) then
c
               do 130 kl = 1 , kx
                  ijkl = kl + ij
                  if (ijkl.eq.ita) then
                     if (iw(ij,kl,1).ne.0) then
                        nadd = nvl(ij)*nvl(kl)
c
                        ninty = ninty + 1
                        do 120 m = 1 , 8
                           intgrl(m,ninty) = -1
 120                    continue
                        intgrl(1,ninty) = ij
                        intgrl(2,ninty) = kl
                        if (kl.ne.1) then
                           intgrl(3,ninty) = iw(ij,kl-1,1)
                           if (kl.gt.2) intgrl(4,ninty) = iw(ij,kl-2,1)
                           if (ij.gt.1) intgrl(5,ninty)
     +                         = iw(ij-1,kl-1,1)
                        else
                           intgrl(3,ninty) = iw(ij-1,kl,1)
                           if (ij.gt.2) intgrl(4,ninty) = iw(ij-2,kl,1)
                        end if
                        iw(ij,kl,1) = mem
                        mem = mem + nadd
                     end if
                  end if
 130           continue
            end if
         end if
 140  continue
      ita = ita + 1
      if (ita.le.(ix+kx)) go to 110
c
c**** find the memory locations of the child kinetic integrals
c**** and the other arrays needed to form the integral- also
c**** give a positon in a list to each integral
c
      ita = 3
      ninty1 = ninty
 150  do 180 ij = 1 , ix
         if (ij.lt.ita) then
            if ((ij+kx).ge.ita) then
c
               do 170 kl = 1 , kx
                  ijkl = kl + ij
                  if (ijkl.eq.ita) then
                     if (iw(ij,kl,2).ne.0) then
                        nadd = nvl(ij)*nvl(kl)
c
                        ninty = ninty + 1
                        do 160 m = 1 , 8
                           intgrl(m,ninty) = -1
 160                    continue
                        intgrl(1,ninty) = ij
                        intgrl(2,ninty) = kl
                        if (kl.ne.1) then
                           intgrl(3,ninty) = iw(ij,kl-1,2)
                           intgrl(4,ninty) = iw(ij,kl,1)
                           if (kl.gt.2) then
                              intgrl(5,ninty) = iw(ij,kl-2,2)
                              intgrl(6,ninty) = iw(ij,kl-2,1)
                           end if
                           if (ij.gt.1) intgrl(7,ninty)
     +                         = iw(ij-1,kl-1,2)
                        else
                           intgrl(3,ninty) = iw(ij-1,kl,2)
                           intgrl(4,ninty) = iw(ij,kl,1)
                           if (ij.gt.2) then
                              intgrl(5,ninty) = iw(ij-2,kl,2)
                              intgrl(6,ninty) = iw(ij-2,kl,1)
                           end if
                        end if
                        iw(ij,kl,2) = mem
                        mem = mem + nadd
                     end if
                  end if
 170           continue
            end if
         end if
 180  continue
      ita = ita + 1
      if (ita.le.(ix+kx)) go to 150
      koveru = iw(ix,kx,1)
      kinup = iw(ix,kx,2)
      koverd = 0
      kindow = 0
      if (ix.gt.2) then
         koverd = iw(ix-2,kx,1)
         kindow = iw(ix-2,kx,2)
      end if
      mem = mem + (nvl(ixz)+nvl(ix))*nvl(kx)
      if (ix.gt.2) mem = mem + nvl(ix-2)*nvl(kx)
      mempto = mem + 5
      ipsize = ifreem/(mempto+1)
c
_IF1()c      write(iwr,*) 'ninty is ',ninty
_IF1()c      write(iwr,104)
_IF1()c104   format(1x,'integral','  value of m','    start pt in memory')
_IF1()      do 200 ij = 1 , ix
_IF1()         if (ij.eq.1) buff1 = 's'
_IF1()         if (ij.eq.2) buff1 = 'p'
_IF1()         if (ij.eq.3) buff1 = 'd'
_IF1()         if (ij.eq.4) buff1 = 'f'
_IF1()         do 190 kl = 1 , kx
_IF1()            if (kl.eq.1) buff2 = 's'
_IF1()            if (kl.eq.2) buff2 = 'p'
_IF1()            if (kl.eq.3) buff2 = 'd'
_IF1()            if (kl.eq.4) buff2 = 'f'
_IF1()            if (iw(ij,kl,1).eq.0) then
_IF1()            end if
_IF1()c      write(iwr,105) buff1,buff2,iw(ij,kl,1)
_IF1()c105   format(1x,'(',a2,'|',a2,')',6x,5x,i10)
_IF1() 190     continue
_IF1() 200  continue
_IF1()      do 220 ij = 1 , ix
_IF1()         if (ij.eq.1) buff1 = 's'
_IF1()         if (ij.eq.2) buff1 = 'p'
_IF1()         if (ij.eq.3) buff1 = 'd'
_IF1()         if (ij.eq.4) buff1 = 'f'
_IF1()         do 210 kl = 1 , kx
_IF1()            if (kl.eq.1) buff2 = 's'
_IF1()            if (kl.eq.2) buff2 = 'p'
_IF1()            if (kl.eq.3) buff2 = 'd'
_IF1()            if (kl.eq.4) buff2 = 'f'
_IF1()            if (iw(ij,kl,2).eq.0) then
_IF1()            end if
_IF1()c      write(iwr,905) buff1,buff2,iw(ij,kl,2)
_IF1()c905   format(1x,'(',a2,'| t |',a2,')',6x,5x,i10)
_IF1() 210     continue
_IF1() 220  continue
_IF1()c      write(iwr,134)
_IF1()c134   format(1x,'integral',' recur eq.',15x,'parents position
_IF1()c     & in memory')
_IF1()      do 230 i = 1 , ninty
_IF1()         if (intgrl(1,i).eq.1) buff1 = 's'
_IF1()         if (intgrl(1,i).eq.2) buff1 = 'p'
_IF1()         if (intgrl(1,i).eq.3) buff1 = 'd'
_IF1()         if (intgrl(1,i).eq.4) buff1 = 'f'
_IF1()         if (intgrl(2,i).eq.1) buff2 = 's'
_IF1()         if (intgrl(2,i).eq.2) buff2 = 'p'
_IF1()         if (intgrl(2,i).eq.3) buff2 = 'd'
_IF1()         if (intgrl(2,i).eq.4) buff2 = 'f'
_IF1()c      write(iwr,135) buff1,buff2,(intgrl(j,i),j=3,8)
_IF1()c135   format(1x,'(',a2,'||',a2,')',2x,i5,7x,5i10)
_IF1() 230  continue
c
c      write(iwr,*) 'ipsize is',ipsize
      return
      end
      subroutine cidrv(iw,iw2,ifreem,ixz,jxz,kxz,lxz,ijxz,klxz,mxz)
c
      implicit REAL  (a-h,o-z)
_IF1()      character *2 buffi,buffij,buff1,buff2
_IF1()      dimension buffij(6),buffi(4)
      dimension iw(ijxz+1,klxz+1,mxz+1),iw2(ijxz+1,klxz+1,6)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mtgrad
c
c**** routine decides for a given derivative integral type which
c**** intermediate integrals are needed
c
      do 40 i = 1 , 5
         do 20 j = 1 , 50
            ihz(i,j) = -1
 20      continue
         do 30 j = 1 , 20
            icona(i,j) = -1
            iconb(i,j) = -1
            iconc(i,j) = -1
            iconx(i,j) = -1
 30      continue
 40   continue
      ijxz1 = ijxz + 1
      klxz1 = klxz + 1
_IF1()      buffi(1) = 's'
_IF1()      buffi(2) = 'p'
_IF1()      buffi(3) = 'd'
_IF1()      buffi(4) = 'f'
_IF1()      buffij(1) = 'ss'
_IF1()      buffij(2) = 'ps'
_IF1()      buffij(3) = 'ds'
_IF1()      buffij(4) = 'fs'
_IF1()      buffij(5) = 'gs'
_IF1()      buffij(6) = 'hs'
      ix = ixz
      jx = jxz
      kx = kxz
      lx = lxz
      ijx = ijxz
      klx = klxz
      mx = mxz
      do 70 m = 1 , mx + 1
         do 60 ij = 1 , ijx + 1
            do 50 kl = 1 , klx + 1
               iw(ij,kl,m) = 0
 50         continue
 60      continue
 70   continue
      do 100 m = 1 , 6
         do 90 ij = 1 , ijx + 1
            do 80 kl = 1 , klx + 1
               iw2(ij,kl,m) = 0
 80         continue
 90      continue
 100  continue
c
c**** the six integral classes need are
c
c**(ix+1 jx | kx lx)
      do 120 ij = ix + 1 , ijx + 1
         do 110 kl = kx , klx
            iw2(ij,kl,1) = -1
            iw(ij,kl,1) = -1
 110     continue
 120  continue
c**(ix jx | kx lx+1)
      do 140 ij = ix , ijx
         do 130 kl = kx , klx + 1
            iw2(ij,kl,2) = -1
            iw(ij,kl,1) = -1
 130     continue
 140  continue
c**(ix jx | kx+1 lx)
      do 160 ij = ix , ijx
         do 150 kl = kx + 1 , klx + 1
            iw2(ij,kl,3) = -1
            iw(ij,kl,1) = -1
 150     continue
 160  continue
c**(ix-1 jx | kx lx )
      if (ix.gt.1) then
         do 180 ij = ix - 1 , ijx - 1
            do 170 kl = kx , klx
               iw2(ij,kl,4) = -1
               iw(ij,kl,1) = -1
 170        continue
 180     continue
      end if
c**(ix jx | kx lx-1 )
      if (lx.gt.1) then
         do 200 ij = ix , ijx
            do 190 kl = kx , klx - 1
               iw2(ij,kl,5) = -1
               iw(ij,kl,1) = -1
 190        continue
 200     continue
      end if
c**(ix jx | kx-1 lx )
      if (kx.gt.1) then
         do 220 ij = ix , ijx
            do 210 kl = kx - 1 , klx - 1
               iw2(ij,kl,6) = -1
               iw(ij,kl,1) = -1
 210        continue
 220     continue
      end if
c
c**** find out all the intermediate integrals needed
c
      maxm = 1
      ita = ijx + klx + 1
 230  do 260 ij = 1 , ijx + 1
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx - 1
            if (ij.ge.itates) then
c
               do 250 kl = 1 , klx + 1
                  ijkl = ij + kl
                  if (ijkl.eq.ita) then
c
                     do 240 m = 1 , maxm
                        if (iw(ij,kl,m).ne.0) then
                           if (kl.gt.klx) then
                              call kcrsd(iw,ijx,klx,mx,ij,kl,m)
                              go to 240
                           end if
                           if (ij.gt.ijx) then
                              call icrsd(iw,ijx,klx,mx,ij,kl,m)
                              go to 240
                           end if
                           if (kl.ne.1) then
                              call kcrsd(iw,ijx,klx,mx,ij,kl,m)
                              go to 240
                           end if
                           call icrsd(iw,ijx,klx,mx,ij,kl,m)
                        end if
 240                 continue
                  end if
 250           continue
            end if
         end if
 260  continue
c
      maxm = maxm + 1
      ita = ita - 1
      if (ita.gt.2) go to 230
c
c**** sort out memory requirements for the ssss(m) integrals
c
      mem = 16
      do 270 i = 1 , mx + 1
         iw(1,1,i) = mem
         mem = mem + 1
 270  continue
      memche = mem + mx + 4
c
c**** find the memory locations of the child integrals
c**** and the other arrays needed to form an integral- also
c**** give a positon in a list to each integral
c
      ita = 3
c      maxm=mx-1
      maxm = mx
      ninty = 0
 280  do 320 ij = 1 , ijx + 1
         itates = ita - 1
         if (ij.le.itates) then
            itates = ita - klx - 1
            if (ij.ge.itates) then
c
               do 310 kl = 1 , klx + 1
                  ijkl = kl + ij
                  if (ijkl.eq.ita) then
                     nadd = nvl(ij)*nvl(kl)
c
                     do 300 m = 1 , maxm
                        if (iw(ij,kl,m).eq.0) go to 300
                        ninty = ninty + 1
                        intgrl(1,ninty) = ij
                        intgrl(2,ninty) = kl
                        if (kl.gt.klx) then
                           call kcrsd2(iw,ijx,klx,mx,ij,kl,m,
     +                                 intgrl(3,ninty))
                           go to 290
                        end if
                        if (ij.gt.ijx) then
                           call icrsd2(iw,ijx,klx,mx,ij,kl,m,
     +                                 intgrl(3,ninty))
                           go to 290
                        end if
                        if (kl.ne.1) then
                           call kcrsd2(iw,ijx,klx,mx,ij,kl,m,
     +                                 intgrl(3,ninty))
                           go to 290
                        end if
                        call icrsd2(iw,ijx,klx,mx,ij,kl,m,
     +                              intgrl(3,ninty))
 290                    iw(ij,kl,m) = mem
                        mem = mem + nadd
 300                 continue
                  end if
 310           continue
            end if
         end if
 320  continue
      maxm = maxm - 1
      ita = ita + 1
      if (maxm.gt.0) go to 280
_IF1()c      write(iwr,*) ' ninty is ',ninty
_IF1()c      write(iwr,104)
_IF1()c104   format(1x,'integral','  value of m','    start pt in memory')
_IF1()      do 350 ij = 1 , ijx + 1
_IF1()         buff1 = buffij(ij)
_IF1()         do 340 kl = 1 , klx + 1
_IF1()            buff2 = buffij(kl)
_IF1()            do 330 m = 1 , mx + 1
_IF1()               if (iw(ij,kl,m).ne.0) then
_IF1()                  m1 = m - 1
_IF1()               end if
_IF1()c      write(iwr,105) buff1,buff2,m1,iw(ij,kl,m)
_IF1()c105   format(1x,'(',a2,'|',a2,')',6x,i5,5x,i10)
_IF1() 330        continue
_IF1() 340     continue
_IF1() 350  continue
_IF1()c      write(iwr,134)
_IF1()c134   format(1x,'integral',' recur eq.',15x,'parents position
_IF1()c     & in memory')
_IF1()      do 360 i = 1 , ninty
_IF1()         buff1 = buffij(intgrl(1,i))
_IF1()         buff2 = buffij(intgrl(2,i))
_IF1()c      write(iwr,135) buff1,buff2,(intgrl(j,i),j=3,8)
_IF1()c135   format(1x,'(',a2,'|',a2,')',2x,i5,7x,5i10)
_IF1() 360  continue
c
c***** now prepare memory for contractions part
c
      if (mem.lt.memche) mem = memche
      mempri = mem
c     mp1 = mem + 1
      memcon = nvl(ix)*nvl(jx)*nvl(kx)*nvl(lx)
      mtemp = memcon
      memcon = memcon + 3
      mtgrad = memcon
      memcon = memcon + 9
c
c**** workout the contraction and hrr for all six integral classes
c
c** start with 3 reduced integral classes
      nconx = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1708)
_IF1()c1708  format(1x)
_IF1()      if (ix.gt.1) then
_IF1()c      write(iwr,1706) buffi(ix-1),buffi(jx),buffi(kx),buffi(lx)
_IF1()c1706  format(1x,'the target class is ','(',a1,a1,'|',a1,a1,')')
_IF1()      end if
_IF1()      if (lx.gt.1) then
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx),buffi(lx-1)
_IF1()      end if
_IF1()      if (kx.gt.1) then
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx-1),buffi(lx)
_IF1()      end if
_IF1()c      write(iwr,1704)
      do 380 ij = 1 , ijx
         do 370 kl = 1 , klx
            if (iw2(ij,kl,4).ne.0 .or. iw2(ij,kl,5).ne.0 .or.
     +          iw2(ij,kl,6).ne.0) then
               nconx = nconx + 1
               iconx(1,nconx) = ij
               iconx(2,nconx) = kl
               iconx(3,nconx) = iw(ij,kl,1)
               iconx(4,nconx) = mem
               iconx(5,nconx) = memcon
               if (iw2(ij,kl,4).ne.0) iw2(ij,kl,4) = memcon
               if (iw2(ij,kl,5).ne.0) iw2(ij,kl,5) = memcon
               if (iw2(ij,kl,6).ne.0) iw2(ij,kl,6) = memcon
_IF1()               buff1 = buffij(ij)
_IF1()               buff2 = buffij(kl)
_IF1()c      write(iwr,1705) buff1,buff2,iconx(3,nconx),iconx(4,nconx)
_IF1()c     &,iconx(5,nconx)
               nadd = nvl(ij)*nvl(kl)
               memcon = memcon + nadd
               mem = mem + nadd
            end if
 370     continue
 380  continue
c
c**   (ix jx | kx+1 lx) integral classes
      nconc = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx+1),buffi(lx)
_IF1()c      write(iwr,1704)
      do 400 ij = ix , ijx
         do 390 kl = kx + 1 , klx + 1
            nconc = nconc + 1
            iconc(1,nconc) = ij
            iconc(2,nconc) = kl
            iconc(3,nconc) = iw(ij,kl,1)
            iconc(4,nconc) = mem
            iconc(5,nconc) = memcon
            iw2(ij,kl,3) = memcon
            nadd = nvl(ij)*nvl(kl)
            memcon = memcon + nadd
            mem = mem + nadd
_IF1()            buff1 = buffij(ij)
_IF1()            buff2 = buffij(kl)
_IF1()c      write(iwr,1705) buff1,buff2,iconc(3,nconc),iconc(4,nconc)
_IF1()c     &,iconc(5,nconc)
 390     continue
 400  continue
c
c**   (ix jx | kx lx+1) integral classes
c
      nconb = 0
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx),buffi(lx+1)
_IF1()c      write(iwr,1704)
      do 420 ij = ix , ijx
         do 410 kl = kx , klx + 1
            nconb = nconb + 1
            iconb(1,nconb) = ij
            iconb(2,nconb) = kl
            iconb(3,nconb) = iw(ij,kl,1)
            iconb(4,nconb) = mem
            iconb(5,nconb) = memcon
            iw2(ij,kl,2) = memcon
            nadd = nvl(ij)*nvl(kl)
            memcon = memcon + nadd
            mem = mem + nadd
_IF1()            buff1 = buffij(ij)
_IF1()            buff2 = buffij(kl)
_IF1()c      write(iwr,1705) buff1,buff2,iconb(3,nconb),iconb(4,nconb)
_IF1()c     &,iconb(5,nconb)
 410     continue
 420  continue
c
c**   (ix+1 jx | kx lx) integral classes
c
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix+1),buffi(jx),buffi(kx),buffi(lx)
_IF1()c      write(iwr,1704)
_IF1()c1704   format(1x,'integral',' uncontracted position',
_IF1()c     &' half-contracted position',
_IF1()c     &' contracted position')
      ncona = 0
      do 440 ij = ix + 1 , ijx + 1
         do 430 kl = kx , klx
            ncona = ncona + 1
            icona(1,ncona) = ij
            icona(2,ncona) = kl
            icona(3,ncona) = iw(ij,kl,1)
            icona(4,ncona) = mem
            icona(5,ncona) = memcon
            iw2(ij,kl,1) = memcon
            nadd = nvl(ij)*nvl(kl)
            memcon = memcon + nadd
            mem = mem + nadd
_IF1()            buff1 = buffij(ij)
_IF1()            buff2 = buffij(kl)
_IF1()c      write(iwr,1705) buff1,buff2,icona(3,ncona),icona(4,ncona)
_IF1()c     &,icona(5,ncona)
_IF1()c1705  format(1x,'(',a2,'|',a2,')',6x,i5,20x,i5,20x,i5)
 430     continue
 440  continue
c
      ncon = ncona + nconb + nconc + nconx
c
      mempto = mem + 1
      ipsize = ifreem/(mempto)
      ictot1 = (mempri-1)*ipsize
      icsiz1 = ictot1/(memcon+1)
c
c**** terms giving d/da(ix jx |kx lx)
      jcount = 0
      mem = memcon
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,*) 'the horizontal recursion relations'
_IF1()c      write(iwr,1708)
_IF1()c      write(iwr,1706) buffi(ix+1),buffi(jx),buffi(kx),buffi(lx)
_IF1()c      if (ix.gt.1)
_IF1()c     &write(iwr,1706) buffi(ix-1),buffi(jx),buffi(kx),buffi(lx)
_IF1()c      write(iwr,1708)
      call hzcdr(ix+1,jx,kx,lx,ix-1,jx,kx,lx,mem,jcount,iw2(1,1,1),
     +           iw2(1,1,4),ijxz1,klxz1)
      memcmx = mem
      mem = icona(5,1)
_IF1()c      write(iwr,*) 'memcmx is ',memcmx
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx),buffi(lx+1)
_IF1()c      if (lx.gt.1)
_IF1()c     &write(iwr,1706) buffi(ix),buffi(jx),buffi(kx),buffi(lx-1)
_IF1()c      write(iwr,1708)
      call hzcdr(ix,jx,kx,lx+1,ix,jx,kx,lx-1,mem,jcount,iw2(1,1,2),
     +           iw2(1,1,5),ijxz1,klxz1)
      memcmx = max(mem,memcmx)
      mem = iconb(5,1)
_IF1()c      write(iwr,*) 'memcmx is ',memcmx
_IF1()c      write(iwr,1706) buffi(ix),buffi(jx),buffi(kx+1),buffi(lx)
_IF1()c      if (kx.gt.1)
_IF1()c     &write(iwr,1706) buffi(ix),buffi(jx),buffi(kx-1),buffi(lx)
_IF1()c      write(iwr,1708)
      call hzcdr(ix,jx,kx+1,lx,ix,jx,kx-1,lx,mem,jcount,iw2(1,1,3),
     +           iw2(1,1,6),ijxz1,klxz1)
      memcmx = max(mem,memcmx)
c      write(iwr,*) 'memcmx is ',memcmx
c
      memcto = memcmx + 1
      icsiz2 = ifreem/(memcto)
      icsize = min(icsiz1,icsiz2)
c
      return
      end
      subroutine countd(i1,ii1,ii2,ixatbf,ixatp,ipar2,n10,na)
c
c******* routine splits the basis functions into batches for integral
c******* evaluation, according to how much memory is available.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension ii1(numspl),ii2(numspl),ixatbf(maxat),ixatp(maxat)
c
_IF1()c     write(iwr,*)'countd',i1,i2,ipar2,na
      n10 = 0
      icmax = 0
      ix1 = i1
      itest2 = ipar2
      ic1 = ixatp(1)
      jc1 = ixatbf(1)
_IF1()c     write(iwr,*)'ic1,jc1',ic1,jc1
      do 20 i = 1 , na
         ii = i + 1
         jc1 = jc1 + ixatbf(ii)
         ic1 = ic1 + ixatp(ii)
_IF1()c     write(iwr,*)'ic1,jc1,itest2',ic1,jc1,itest2
         if (ic1.gt.itest2 .or. i.eq.na) then
            ic1 = ic1 - ixatp(ii)
            if (ic1.gt.icmax) icmax = ic1
            n10 = n10 + 1
            ii1(n10) = ix1
            ix1 = ix1 + jc1 - ixatbf(ii)
            ii2(n10) = ix1 - 1
            ic1 = ixatp(ii)
            jc1 = ixatbf(ii)
         end if
 20   continue
      if (icmax.lt.ipar2) ipar2 = icmax
_IF1()c     write(iwr,*)'icmax,ipar2'
c
      return
      end
      subroutine dd1(xijk,xija,xijp1,xijp2,xijp3,i1,i2
     &,j1,j2,xt1,xt2,wx
     &,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ccc1,ccc2,ext1)
c
c
c........ bad code !
c
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (d/d) and (d/h/d) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
INCLUDE(common/sizes)
      logical ext1
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      dimension xija(n11),xijk(n11),ccc2(n11,3)
     &,s1(n11*6,6),h1(n11*6,6),sx(nb,nb),hx(nb,nb),ccc1(n11,3)
     &,lar(n11),f(9*n11),t(n11),xt1(n14*6,n13*6),xt2(n16*6,n15*6)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
     &,xijp1(n11),xijp2(n11),xijp3(n11)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11*36)
      icount = 0
c
      do 30 i = imin , imax
         if (ext1) jmax = i
         do 20 j = jmin , jmax
            icount = icount + 1
            ccc1(icount,1) = xijp1(icount) - coprm(i,1)
            ccc1(icount,2) = xijp2(icount) - coprm(i,2)
            ccc1(icount,3) = xijp3(icount) - coprm(i,3)
            ccc2(icount,1) = xijp1(icount) - coprm(j,1)
            ccc2(icount,2) = xijp2(icount) - coprm(j,2)
            ccc2(icount,3) = xijp3(icount) - coprm(j,3)
 20      continue
 30   continue
      do 60 k = 1 , na
c
         do 40 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 40      continue
c
         call fquik2(f,f,lar,t,5,n11)
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 50 icount = 1 , n11
            index1 = n11 + icount
            index2 = n11 + index1
            index3 = n11 + index2
            index4 = n11 + index3
            index5 = n11 + index4
            x7 = xijk(icount)*charge(k)/pi14
            f(icount) = f(icount)*x7
            f(index1) = f(index1)*x7
            f(index2) = f(index2)*x7
            f(index3) = f(index3)*x7
            f(index4) = f(index4)*x7
            xij2 = 1.0d0/(xija(icount)+xija(icount))
            x40 = (f(icount)-f(index1))*xij2
            x41 = (f(index1)-f(index2))*xij2
            x42 = (f(index2)-f(index3))*xij2
            x1a = xijp1(icount) - coord(1,k)
            x1b = xijp2(icount) - coord(2,k)
            x1c = xijp3(icount) - coord(3,k)
c
c     ***** xx1 has the k dependant factors of (pi/v/s)m=0 ***** 
c     ***** xx2 has the k dependant factors of (pi/v/s)m=1 ***** 
c     ***** xx3 has the k dependant factors of (pi/v/s)m=2 ***** 
c     ***** xx3b has (pi/v/s)m=3 *******************************  
c     ***** xx4 contains (1/2@)*( (pi/v/s)m=0 - (pi/v/s)m=1 )
c     ***** xx5 has (1/2@)*( (pi/v/s)m=1 - (pi/v/s)m=2 )
c
            x1x = ccc1(icount,1)*f(icount) - x1a*f(index1)
            x2x = ccc1(icount,1)*f(index1) - x1a*f(index2)
            x3x = ccc1(icount,1)*f(index2) - x1a*f(index3)
            x3bx = ccc1(icount,1)*f(index3) - x1a*f(index4)
            x4x = xij2*(x1x-x2x)
            x5x = xij2*(x2x-x3x)
            x1y = ccc1(icount,2)*f(icount) - x1b*f(index1)
            x2y = ccc1(icount,2)*f(index1) - x1b*f(index2)
            x3y = ccc1(icount,2)*f(index2) - x1b*f(index3)
            x3by = ccc1(icount,2)*f(index3) - x1b*f(index4)
            x4y = xij2*(x1y-x2y)
            x5y = xij2*(x2y-x3y)
            x1z = ccc1(icount,3)*f(icount) - x1c*f(index1)
            x2z = ccc1(icount,3)*f(index1) - x1c*f(index2)
            x3z = ccc1(icount,3)*f(index2) - x1c*f(index3)
            x3bz = ccc1(icount,3)*f(index3) - x1c*f(index4)
            x4z = xij2*(x1z-x2z)
            x5z = xij2*(x2z-x3z)
c
c***** xd1, xd2 and xd3 have (d|v|s)m=0, (d|v|s)m=1 and (d|v|s)m=2
c**** dxx
            xd1a = ccc1(icount,1)*x1x - x1a*x2x + x40
            xd2a = ccc1(icount,1)*x2x - x1a*x3x + x41
            xd3a = ccc1(icount,1)*x3x - x1a*x3bx + x42
c**** dxy
            xd1b = ccc1(icount,2)*x1x - x1b*x2x
            xd2b = ccc1(icount,2)*x2x - x1b*x3x
            xd3b = ccc1(icount,2)*x3x - x1b*x3bx
c**** dxz
            xd1c = ccc1(icount,3)*x1x - x1c*x2x
            xd2c = ccc1(icount,3)*x2x - x1c*x3x
            xd3c = ccc1(icount,3)*x3x - x1c*x3bx
c**** dyy
            xd1d = ccc1(icount,2)*x1y - x1b*x2y + x40
            xd2d = ccc1(icount,2)*x2y - x1b*x3y + x41
            xd3d = ccc1(icount,2)*x3y - x1b*x3by + x42
c**** dyz
            xd1e = ccc1(icount,3)*x1y - x1c*x2y
            xd2e = ccc1(icount,3)*x2y - x1c*x3y
            xd3e = ccc1(icount,3)*x3y - x1c*x3by
c**** dzz
            xd1f = ccc1(icount,3)*x1z - x1c*x2z + x40
            xd2f = ccc1(icount,3)*x2z - x1c*x3z + x41
            xd3f = ccc1(icount,3)*x3z - x1c*x3bz + x42
c
c**** form (1/2@)*( (pi|v|pj)m=0 - (pi|v|pj)m=1 ) in ppxx using
c**** a horizontal recursion relation
            ab1 = ccc2(icount,1) - ccc1(icount,1)
            ab2 = ccc2(icount,2) - ccc1(icount,2)
            ab3 = ccc2(icount,3) - ccc1(icount,3)
            ppxx = xij2*(xd1a-xd2a+ab1*(x1x-x2x))
            ppyy = xij2*(xd1d-xd2d+ab2*(x1y-x2y))
            ppzz = xij2*(xd1f-xd2f+ab3*(x1z-x2z))
            ppxy = xij2*(xd1b-xd2b+ab2*(x1x-x2x))
            ppyx = xij2*(xd1b-xd2b+ab1*(x1y-x2y))
            ppxz = xij2*(xd1c-xd2c+ab3*(x1x-x2x))
            ppzx = xij2*(xd1c-xd2c+ab1*(x1z-x2z))
            ppzy = xij2*(xd1e-xd2e+ab2*(x1z-x2z))
            ppyz = xij2*(xd1e-xd2e+ab3*(x1y-x2y))
c***** form (1/2@)*( (dij|v|s)m=0 - (dij|v|s)m=1 ) in ddxx
            ddxx = xij2*(xd1a-xd2a)
            ddxy = xij2*(xd1b-xd2b)
            ddxz = xij2*(xd1c-xd2c)
            ddyy = xij2*(xd1d-xd2d)
            ddyz = xij2*(xd1e-xd2e)
            ddzz = xij2*(xd1f-xd2f)
c
c***** (dxx|v|px)m=0 and (dxx|v|px)m=1 go in dxxpx and dxxpx1
c
c***** (dxx|v|px-z)
            dxxpx = ccc2(icount,1)*xd1a - x1a*xd2a + x4x + x4x
            dxxpy = ccc2(icount,2)*xd1a - x1b*xd2a
            dxxpz = ccc2(icount,3)*xd1a - x1c*xd2a
            dxxpx1 = ccc2(icount,1)*xd2a - x1a*xd3a + x5x + x5x
            dxxpy1 = ccc2(icount,2)*xd2a - x1b*xd3a
            dxxpz1 = ccc2(icount,3)*xd2a - x1c*xd3a
c***** (dxy|v|px-z)
            dxypx = ccc2(icount,1)*xd1b - x1a*xd2b + x4y
            dxypy = ccc2(icount,2)*xd1b - x1b*xd2b + x4x
            dxypz = ccc2(icount,3)*xd1b - x1c*xd2b
            dxypx1 = ccc2(icount,1)*xd2b - x1a*xd3b + x5y
            dxypy1 = ccc2(icount,2)*xd2b - x1b*xd3b + x5x
            dxypz1 = ccc2(icount,3)*xd2b - x1c*xd3b
c***** (dxz|v|px-z)
            dxzpx = ccc2(icount,1)*xd1c - x1a*xd2c + x4z
            dxzpy = ccc2(icount,2)*xd1c - x1b*xd2c
            dxzpz = ccc2(icount,3)*xd1c - x1c*xd2c + x4x
            dxzpx1 = ccc2(icount,1)*xd2c - x1a*xd3c + x5z
            dxzpy1 = ccc2(icount,2)*xd2c - x1b*xd3c
            dxzpz1 = ccc2(icount,3)*xd2c - x1c*xd3c + x5x
c***** (dyy|v|px-z)
            dyypx = ccc2(icount,1)*xd1d - x1a*xd2d
            dyypy = ccc2(icount,2)*xd1d - x1b*xd2d + x4y + x4y
            dyypz = ccc2(icount,3)*xd1d - x1c*xd2d
            dyypx1 = ccc2(icount,1)*xd2d - x1a*xd3d
            dyypy1 = ccc2(icount,2)*xd2d - x1b*xd3d + x5y + x5y
            dyypz1 = ccc2(icount,3)*xd2d - x1c*xd3d
c***** (dyz|v|px-z)
            dyzpx = ccc2(icount,1)*xd1e - x1a*xd2e
            dyzpy = ccc2(icount,2)*xd1e - x1b*xd2e + x4z
            dyzpz = ccc2(icount,3)*xd1e - x1c*xd2e + x4y
            dyzpx1 = ccc2(icount,1)*xd2e - x1a*xd3e
            dyzpy1 = ccc2(icount,2)*xd2e - x1b*xd3e + x5z
            dyzpz1 = ccc2(icount,3)*xd2e - x1c*xd3e + x5y
c***** (dzz|v|px-z)
            dzzpx = ccc2(icount,1)*xd1f - x1a*xd2f
            dzzpy = ccc2(icount,2)*xd1f - x1b*xd2f
            dzzpz = ccc2(icount,3)*xd1f - x1c*xd2f + x4z + x4z
            dzzpx1 = ccc2(icount,1)*xd2f - x1a*xd3f
            dzzpy1 = ccc2(icount,2)*xd2f - x1b*xd3f
            dzzpz1 = ccc2(icount,3)*xd2f - x1c*xd3f + x5z + x5z
c
c***** now the answer goes in the array h1
c
c***** (dxx|v|dkl)
            h1(icount,1) = h1(icount,1) + ccc2(icount,1)
     +                     *dxxpx - x1a*dxxpx1 + ddxx + ppxx + ppxx
            h1(index1,1) = h1(index1,1) + ccc2(icount,2)
     +                     *dxxpx - x1b*dxxpx1
            h1(index2,1) = h1(index2,1) + ccc2(icount,3)
     +                     *dxxpx - x1c*dxxpx1
            h1(index3,1) = h1(index3,1) + ccc2(icount,2)
     +                     *dxxpy - x1b*dxxpy1 + ddxx
            h1(index4,1) = h1(index4,1) + ccc2(icount,3)
     +                     *dxxpy - x1c*dxxpy1
            h1(index5,1) = h1(index5,1) + ccc2(icount,3)
     +                     *dxxpz - x1c*dxxpz1 + ddxx
c***** (dxy|v|dkl)
            h1(icount,2) = h1(icount,2) + ccc2(icount,1)
     +                     *dxypx - x1a*dxypx1 + ddxy + ppyx
            h1(index1,2) = h1(index1,2) + ccc2(icount,2)
     +                     *dxypx - x1b*dxypx1 + ppxx
            h1(index2,2) = h1(index2,2) + ccc2(icount,3)
     +                     *dxypx - x1c*dxypx1
            h1(index3,2) = h1(index3,2) + ccc2(icount,2)
     +                     *dxypy - x1b*dxypy1 + ddxy + ppxy
            h1(index4,2) = h1(index4,2) + ccc2(icount,3)
     +                     *dxypy - x1c*dxypy1
            h1(index5,2) = h1(index5,2) + ccc2(icount,3)
     +                     *dxypz - x1c*dxypz1 + ddxy
c***** (dxz|v|dkl)
            h1(icount,3) = h1(icount,3) + ccc2(icount,1)
     +                     *dxzpx - x1a*dxzpx1 + ddxz + ppzx
            h1(index1,3) = h1(index1,3) + ccc2(icount,2)
     +                     *dxzpx - x1b*dxzpx1
            h1(index2,3) = h1(index2,3) + ccc2(icount,3)
     +                     *dxzpx - x1c*dxzpx1 + ppxx
            h1(index3,3) = h1(index3,3) + ccc2(icount,2)
     +                     *dxzpy - x1b*dxzpy1 + ddxz
            h1(index4,3) = h1(index4,3) + ccc2(icount,3)
     +                     *dxzpy - x1c*dxzpy1 + ppxy
            h1(index5,3) = h1(index5,3) + ccc2(icount,3)
     +                     *dxzpz - x1c*dxzpz1 + ddxz + ppxz
c***** (dyy|v|dkl)
            h1(icount,4) = h1(icount,4) + ccc2(icount,1)
     +                     *dyypx - x1a*dyypx1 + ddyy
            h1(index1,4) = h1(index1,4) + ccc2(icount,2)
     +                     *dyypx - x1b*dyypx1 + ppyx + ppyx
            h1(index2,4) = h1(index2,4) + ccc2(icount,3)
     +                     *dyypx - x1c*dyypx1
            h1(index3,4) = h1(index3,4) + ccc2(icount,2)
     +                     *dyypy - x1b*dyypy1 + ddyy + ppyy + ppyy
            h1(index4,4) = h1(index4,4) + ccc2(icount,3)
     +                     *dyypy - x1c*dyypy1
            h1(index5,4) = h1(index5,4) + ccc2(icount,3)
     +                     *dyypz - x1c*dyypz1 + ddyy
c***** (dyz|v|dkl)
            h1(icount,5) = h1(icount,5) + ccc2(icount,1)
     +                     *dyzpx - x1a*dyzpx1 + ddyz
            h1(index1,5) = h1(index1,5) + ccc2(icount,2)
     +                     *dyzpx - x1b*dyzpx1 + ppzx
            h1(index2,5) = h1(index2,5) + ccc2(icount,3)
     +                     *dyzpx - x1c*dyzpx1 + ppyx
            h1(index3,5) = h1(index3,5) + ccc2(icount,2)
     +                     *dyzpy - x1b*dyzpy1 + ddyz + ppzy
            h1(index4,5) = h1(index4,5) + ccc2(icount,3)
     +                     *dyzpy - x1c*dyzpy1 + ppyy
            h1(index5,5) = h1(index5,5) + ccc2(icount,3)
     +                     *dyzpz - x1c*dyzpz1 + ddyz + ppyz
c***** (dzz|v|dkl)
            h1(icount,6) = h1(icount,6) + ccc2(icount,1)
     +                     *dzzpx - x1a*dzzpx1 + ddzz
            h1(index1,6) = h1(index1,6) + ccc2(icount,2)
     +                     *dzzpx - x1b*dzzpx1
            h1(index2,6) = h1(index2,6) + ccc2(icount,3)
     +                     *dzzpx - x1c*dzzpx1 + ppzx + ppzx
            h1(index3,6) = h1(index3,6) + ccc2(icount,2)
     +                     *dzzpy - x1b*dzzpy1 + ddzz
            h1(index4,6) = h1(index4,6) + ccc2(icount,3)
     +                     *dzzpy - x1c*dzzpy1 + ppzy + ppzy
            h1(index5,6) = h1(index5,6) + ccc2(icount,3)
     +                     *dzzpz - x1c*dzzpz1 + ddzz + ppzz + ppzz
 50      continue
 60   continue
c
      icount = 0
      do 80 i = imin , imax
         if (ext1) jmax = i
         do 70 j = jmin , jmax
            icount = icount + 1
            index1 = n11 + icount
            index2 = n11 + index1
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(index2) = pe(j)
            f(index1) = pe(i)
            f(icount) = pe(i)*pe(j)
 70      continue
 80   continue
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 90 icount = 1 , n11
         index1 = n11 + icount
         index2 = n11 + index1
         index3 = n11 + index2
         index4 = n11 + index3
         index5 = n11 + index4
c
c     ***** x4 is (s||s), x3 is (s|t|s)
c
         x1 = 1.0d0/xija(icount)
         x2 = f(icount)*x1
         x4 = xijk(icount)*pi14*dsqrt(x1)
         x1 = 0.5d0*x1
         x22 = 2.0d0*x2
         x3 = x2*(3.0d0-x22*t(icount))*x4
         x8 = x2*x4/f(index1)
         x9 = x2/f(index2)
c
c***** x1x-z is (pi||s) and xtx-z is (pi|v|s)
         x1x = ccc1(icount,1)*x4
         x1y = ccc1(icount,2)*x4
         x1z = ccc1(icount,3)*x4
c
         xtx = ccc1(icount,1)*x3 + x22*x1x
         xty = ccc1(icount,2)*x3 + x22*x1y
         xtz = ccc1(icount,3)*x3 + x22*x1z
c***** xdxx-xdzz is (dij||s) and xtxx-xtzz is (dij|v|s)
c***** xpxx-xpzz is (pi||pj) and xptxx-xptzz is (pi|v|pj)
         x6 = x1*x4
         xdxx = ccc1(icount,1)*x1x + x6
         xdxy = ccc1(icount,2)*x1x
         xdxz = ccc1(icount,3)*x1x
         xdyy = ccc1(icount,2)*x1y + x6
         xdyz = ccc1(icount,3)*x1y
         xdzz = ccc1(icount,3)*x1z + x6
         xpxx = ccc2(icount,1)*x1x + x6
         xpxy = ccc2(icount,2)*x1x
         xpxz = ccc2(icount,3)*x1x
         xpyy = ccc2(icount,2)*x1y + x6
         xpyz = ccc2(icount,3)*x1y
         xpzz = ccc2(icount,3)*x1z + x6
         xpzy = xpyz
         xpyx = xpxy
         xpzx = xpxz
c
         xt6 = x1*x3 - x8
         xtxx = ccc1(icount,1)*xtx + xdxx*x22 + xt6
         xtxy = ccc1(icount,2)*xtx + xdxy*x22
         xtxz = ccc1(icount,3)*xtx + xdxz*x22
         xtyy = ccc1(icount,2)*xty + xdyy*x22 + xt6
         xtyz = ccc1(icount,3)*xty + xdyz*x22
         xtzz = ccc1(icount,3)*xtz + xdzz*x22 + xt6
         xt6 = xt6 + x8
         xptxx = x1*(ccc2(icount,1)*xtx+xpxx*x22+xt6)
         xptxy = x1*(ccc2(icount,2)*xtx+xpxy*x22)
         xptxz = x1*(ccc2(icount,3)*xtx+xpxz*x22)
         xptyy = x1*(ccc2(icount,2)*xty+xpyy*x22+xt6)
         xptyz = x1*(ccc2(icount,3)*xty+xpyz*x22)
         xptzz = x1*(ccc2(icount,3)*xtz+xpzz*x22+xt6)
         xpxx = x1*xpxx
         xpxy = x1*xpxy
         xpxz = x1*xpxz
         xpyy = x1*xpyy
         xpyz = x1*xpyz
         xpzz = x1*xpzz
         xpzy = xpyz
         xpyx = xpxy
         xpzx = xpxz
         xptzy = xptyz
         xptyx = xptxy
         xptzx = xptxz
c
c***** now do the overlap integrals in s1 and the hcore integrals
c***** in h1
         x1x = x1x*x1
         x1y = x1y*x1
         x1z = x1z*x1
         xtx = xtx*x1
         xty = xty*x1
         xtz = xtz*x1
c***** (dxx||px-z) and (dxx|h|px-z)
         dxxpx = ccc2(icount,1)*xdxx + x1x + x1x
         dxxpy = ccc2(icount,2)*xdxx
         dxxpz = ccc2(icount,3)*xdxx
         txxpx = ccc2(icount,1)*xtxx + x22*dxxpx + xtx + xtx
         txxpy = ccc2(icount,2)*xtxx + x22*dxxpy
         txxpz = ccc2(icount,3)*xtxx + x22*dxxpz
c***** (dxy||px-z) and (dxy|h|px-z)
         dxypx = ccc2(icount,1)*xdxy + x1y
         dxypy = ccc2(icount,2)*xdxy + x1x
         dxypz = ccc2(icount,3)*xdxy
         txypx = ccc2(icount,1)*xtxy + x22*dxypx + xty
         txypy = ccc2(icount,2)*xtxy + x22*dxypy + xtx
         txypz = ccc2(icount,3)*xtxy + x22*dxypz
c***** (dxz||px-z) and (dxz|h|px-z)
         dxzpx = ccc2(icount,1)*xdxz + x1z
         dxzpy = ccc2(icount,2)*xdxz
         dxzpz = ccc2(icount,3)*xdxz + x1x
         txzpx = ccc2(icount,1)*xtxz + x22*dxzpx + xtz
         txzpy = ccc2(icount,2)*xtxz + x22*dxzpy
         txzpz = ccc2(icount,3)*xtxz + x22*dxzpz + xtx
c***** (dyy||px-z) and (dyy|h|px-z)
         dyypx = ccc2(icount,1)*xdyy
         dyypy = ccc2(icount,2)*xdyy + x1y + x1y
         dyypz = ccc2(icount,3)*xdyy
         tyypx = ccc2(icount,1)*xtyy + x22*dyypx
         tyypy = ccc2(icount,2)*xtyy + x22*dyypy + xty + xty
         tyypz = ccc2(icount,3)*xtyy + x22*dyypz
c***** (dyz||px-z) and (dyz|h|px-z)
         dyzpx = ccc2(icount,1)*xdyz
         dyzpy = ccc2(icount,2)*xdyz + x1z
         dyzpz = ccc2(icount,3)*xdyz + x1y
         tyzpx = ccc2(icount,1)*xtyz + x22*dyzpx
         tyzpy = ccc2(icount,2)*xtyz + x22*dyzpy + xtz
         tyzpz = ccc2(icount,3)*xtyz + x22*dyzpz + xty
c***** (dzz||px-z) and (dzz|h|px-z)
         dzzpx = ccc2(icount,1)*xdzz
         dzzpy = ccc2(icount,2)*xdzz
         dzzpz = ccc2(icount,3)*xdzz + x1z + x1z
         tzzpx = ccc2(icount,1)*xtzz + x22*dzzpx
         tzzpy = ccc2(icount,2)*xtzz + x22*dzzpy
         tzzpz = ccc2(icount,3)*xtzz + x22*dzzpz + xtz + xtz
c
c***** the answers are put in h1 and s1
c
         xtxx = xtxx*x1 - xdxx*x9
         xtxy = xtxy*x1 - xdxy*x9
         xtxz = xtxz*x1 - xdxz*x9
         xtyy = xtyy*x1 - xdyy*x9
         xtyz = xtyz*x1 - xdyz*x9
         xtzz = xtzz*x1 - xdzz*x9
         xdxx = xdxx*x1
         xdxy = xdxy*x1
         xdxz = xdxz*x1
         xdyy = xdyy*x1
         xdyz = xdyz*x1
         xdzz = xdzz*x1
c***** (dxx|v|dkl)
         s1(icount,1) = ccc2(icount,1)*dxxpx + xdxx + xpxx + xpxx
         s1(index1,1) = ccc2(icount,2)*dxxpx
         s1(index2,1) = ccc2(icount,3)*dxxpx
         s1(index3,1) = ccc2(icount,2)*dxxpy + xdxx
         s1(index4,1) = ccc2(icount,3)*dxxpy
         s1(index5,1) = ccc2(icount,3)*dxxpz + xdxx
         h1(icount,1) = ccc2(icount,1)*txxpx + xtxx + xptxx + xptxx +
     +                  x22*s1(icount,1) - h1(icount,1)
         h1(index1,1) = ccc2(icount,2)*txxpx + x22*s1(index1,1)
     +                  - h1(index1,1)
         h1(index2,1) = ccc2(icount,3)*txxpx + x22*s1(index2,1)
     +                  - h1(index2,1)
         h1(index3,1) = ccc2(icount,2)*txxpy + xtxx + x22*s1(index3,1)
     +                  - h1(index3,1)
         h1(index4,1) = ccc2(icount,3)*txxpy + x22*s1(index4,1)
     +                  - h1(index4,1)
         h1(index5,1) = ccc2(icount,3)*txxpz + xtxx + x22*s1(index5,1)
     +                  - h1(index5,1)
c***** (dxy|v|dkl)
         s1(icount,2) = ccc2(icount,1)*dxypx + xdxy + xpyx
         s1(index1,2) = ccc2(icount,2)*dxypx + xpxx
         s1(index2,2) = ccc2(icount,3)*dxypx
         s1(index3,2) = ccc2(icount,2)*dxypy + xdxy + xpxy
         s1(index4,2) = ccc2(icount,3)*dxypy
         s1(index5,2) = ccc2(icount,3)*dxypz + xdxy
         h1(icount,2) = ccc2(icount,1)*txypx + xtxy + xptyx +
     +                  x22*s1(icount,2) - h1(icount,2)
         h1(index1,2) = ccc2(icount,2)*txypx + xptxx + x22*s1(index1,2)
     +                  - h1(index1,2)
         h1(index2,2) = ccc2(icount,3)*txypx + x22*s1(index2,2)
     +                  - h1(index2,2)
         h1(index3,2) = ccc2(icount,2)*txypy + xtxy + xptxy +
     +                  x22*s1(index3,2) - h1(index3,2)
         h1(index4,2) = ccc2(icount,3)*txypy + x22*s1(index4,2)
     +                  - h1(index4,2)
         h1(index5,2) = ccc2(icount,3)*txypz + xtxy + x22*s1(index5,2)
     +                  - h1(index5,2)
c***** (dxz|v|dkl)
         s1(icount,3) = ccc2(icount,1)*dxzpx + xdxz + xpzx
         s1(index1,3) = ccc2(icount,2)*dxzpx
         s1(index2,3) = ccc2(icount,3)*dxzpx + xpxx
         s1(index3,3) = ccc2(icount,2)*dxzpy + xdxz
         s1(index4,3) = ccc2(icount,3)*dxzpy + xpxy
         s1(index5,3) = ccc2(icount,3)*dxzpz + xdxz + xpxz
         h1(icount,3) = ccc2(icount,1)*txzpx + xtxz + xptzx +
     +                  x22*s1(icount,3) - h1(icount,3)
         h1(index1,3) = ccc2(icount,2)*txzpx + x22*s1(index1,3)
     +                  - h1(index1,3)
         h1(index2,3) = ccc2(icount,3)*txzpx + xptxx + x22*s1(index2,3)
     +                  - h1(index2,3)
         h1(index3,3) = ccc2(icount,2)*txzpy + xtxz + x22*s1(index3,3)
     +                  - h1(index3,3)
         h1(index4,3) = ccc2(icount,3)*txzpy + xptxy + x22*s1(index4,3)
     +                  - h1(index4,3)
         h1(index5,3) = ccc2(icount,3)*txzpz + xtxz + xptxz +
     +                  x22*s1(index5,3) - h1(index5,3)
c***** (dyy|v|dkl)
         s1(icount,4) = ccc2(icount,1)*dyypx + xdyy
         s1(index1,4) = ccc2(icount,2)*dyypx + xpyx + xpyx
         s1(index2,4) = ccc2(icount,3)*dyypx
         s1(index3,4) = ccc2(icount,2)*dyypy + xdyy + xpyy + xpyy
         s1(index4,4) = ccc2(icount,3)*dyypy
         s1(index5,4) = ccc2(icount,3)*dyypz + xdyy
         h1(icount,4) = ccc2(icount,1)*tyypx + xtyy + x22*s1(icount,4)
     +                  - h1(icount,4)
         h1(index1,4) = ccc2(icount,2)*tyypx + xptyx + xptyx +
     +                  x22*s1(index1,4) - h1(index1,4)
         h1(index2,4) = ccc2(icount,3)*tyypx + x22*s1(index2,4)
     +                  - h1(index2,4)
         h1(index3,4) = ccc2(icount,2)*tyypy + xtyy + xptyy + xptyy +
     +                  x22*s1(index3,4) - h1(index3,4)
         h1(index4,4) = ccc2(icount,3)*tyypy + x22*s1(index4,4)
     +                  - h1(index4,4)
         h1(index5,4) = ccc2(icount,3)*tyypz + xtyy + x22*s1(index5,4)
     +                  - h1(index5,4)
c***** (dyz|v|dkl)
         s1(icount,5) = ccc2(icount,1)*dyzpx + xdyz
         s1(index1,5) = ccc2(icount,2)*dyzpx + xpzx
         s1(index2,5) = ccc2(icount,3)*dyzpx + xpyx
         s1(index3,5) = ccc2(icount,2)*dyzpy + xdyz + xpzy
         s1(index4,5) = ccc2(icount,3)*dyzpy + xpyy
         s1(index5,5) = ccc2(icount,3)*dyzpz + xdyz + xpyz
         h1(icount,5) = ccc2(icount,1)*tyzpx + xtyz + x22*s1(icount,5)
     +                  - h1(icount,5)
         h1(index1,5) = ccc2(icount,2)*tyzpx + xptzx + x22*s1(index1,5)
     +                  - h1(index1,5)
         h1(index2,5) = ccc2(icount,3)*tyzpx + xptyx + x22*s1(index2,5)
     +                  - h1(index2,5)
         h1(index3,5) = ccc2(icount,2)*tyzpy + xtyz + xptzy +
     +                  x22*s1(index3,5) - h1(index3,5)
         h1(index4,5) = ccc2(icount,3)*tyzpy + xptyy + x22*s1(index4,5)
     +                  - h1(index4,5)
         h1(index5,5) = ccc2(icount,3)*tyzpz + xtyz + xptyz +
     +                  x22*s1(index5,5) - h1(index5,5)
c***** (dzz|v|dkl)
         s1(icount,6) = ccc2(icount,1)*dzzpx + xdzz
         s1(index1,6) = ccc2(icount,2)*dzzpx
         s1(index2,6) = ccc2(icount,3)*dzzpx + xpzx + xpzx
         s1(index3,6) = ccc2(icount,2)*dzzpy + xdzz
         s1(index4,6) = ccc2(icount,3)*dzzpy + xpzy + xpzy
         s1(index5,6) = ccc2(icount,3)*dzzpz + xdzz + xpzz + xpzz
         h1(icount,6) = ccc2(icount,1)*tzzpx + xtzz + x22*s1(icount,6)
     +                  - h1(icount,6)
         h1(index1,6) = ccc2(icount,2)*tzzpx + x22*s1(index1,6)
     +                  - h1(index1,6)
         h1(index2,6) = ccc2(icount,3)*tzzpx + xptzx + xptzx +
     +                  x22*s1(index2,6) - h1(index2,6)
         h1(index3,6) = ccc2(icount,2)*tzzpy + xtzz + x22*s1(index3,6)
     +                  - h1(index3,6)
         h1(index4,6) = ccc2(icount,3)*tzzpy + xptzy + xptzy +
     +                  x22*s1(index4,6) - h1(index4,6)
         h1(index5,6) = ccc2(icount,3)*tzzpz + xtzz + xptzz + xptzz +
     +                  x22*s1(index5,6) - h1(index5,6)
 90   continue
c
c*****  put the result in the matrices to be transformed ***
c
      do 130 ic = 1 , 6
         icount = 0
         do 120 jc = 1 , 6
            jmax = n16
            do 110 i = 1 , n14
               ip = (i-1)*6
               if (ext1) jmax = i
               do 100 j = 1 , jmax
                  jp = (j-1)*6
                  icount = icount + 1
                  wx(ip+ic,jp+jc) = s1(icount,ic)
                  wy(ip+ic,jp+jc) = h1(icount,ic)
 100           continue
 110        continue
 120     continue
 130  continue
c
      if (ext1) then
         do 150 i = 1 , n14*6
            do 140 j = 1 , i
               wx(j,i) = wx(i,j)
               wy(j,i) = wy(i,j)
 140        continue
 150     continue
      end if
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16*6,wz,1,ipar2,n14*6,n16*6,n15*6)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wx,ipar2,1,n15*6,n14*6,n13*6)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16*6,wz,1,ipar2,n14*6,n16*6,n15*6)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wy,ipar2,1,n15*6,n14*6,n13*6)
c
      ix1 = 6*i1 - 3*np - ns*5 - 5
      ix2 = ix1 + n13*6 - 1
      jx1 = 6*j1 - 3*np - ns*5 - 5
      jx2 = jx1 + n15*6 - 1
      ic = 0
      do 170 i = ix1 , ix2
         ic = ic + 1
         jc = 0
         if (ext1) jx2 = i
         do 160 j = jx1 , jx2
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 160     continue
 170  continue
c
      return
      end
      subroutine def1dr(d,rdens,dens,grad,iwr)
c
c.... routine coordinates the calculation of the 1-electron derivative
c.... integrals
      implicit REAL  (a-h,o-z)
c     character *2 buff
c
c***** maxat is max number of atoms,maxorb is max num of basis
c***** functions,mxprim is max num of distinct primitives (ie px,py
c***** and pz count only once)
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
      dimension d(*)
      common/igarg/igarm(40),iparr(3)
      common/tempit/zitest,zipass,zjtest,zjpass
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
      dimension grad(na,9),dens(nb,nb),rdens(nb,nb)
c
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
      logical glog1,glog2,glog3,glog5,glog6,gdbg
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10)
c
      logical ext1
c
c     For memory allocation tracking
      character *8 fnm,snm
      data fnm,snm/"direct.m","def1dr"/

c
c*****  a few comments about the memory allocation. most of the
c*****  large memory is contained in d. the first part of this
c*****  array is taken up with arrays of size nb*nb such as the
c*****  fock matrix. these arrays should end at igarm(12). the rest
c*****  is free memory to evaluate integrals in. the first part of
c*****  this is to store the parameters in.
      ipss = ipars*ipars
      ipps = iparp*ipars
      ippp = iparp*iparp
      ipds = ipard*ipars
      ipdp = ipard*iparp
      ipdd = ipard*ipard
c
c*********************************************
c***** the one electron integral calls********
c*********************************************
c
c
      if (nd.eq.0) ndspl = 0
      nspl = nsspl + npspl
      nspdl = nspl + ndspl
c
c     allocate available memory
c
      k0a = igmem_alloc_all_inf(imense,fnm,snm,'k0a',IGMEM_DEBUG)

c     
      do 30 icart = 1 , 9
         do 20 i = 1 , na
            grad(i,icart) = 0.0d0
 20      continue
 30   continue
c
c***** now the loop over batches for the i position
c
      ic1 = 0
      do 80 ix = 1 , nspdl
         if (ix.le.nsspl) then
            itemp = ix
            i1 = idefs1(itemp)
            i2 = idefs2(itemp)
            ixz1 = 1
            go to 40
         end if
         if (ix.le.nspl) then
            itemp = ix - nsspl
            i1 = idefp1(itemp)
            i2 = idefp2(itemp)
            ixz1 = 2
            go to 40
         end if
         itemp = ix - nspl
         i1 = idefd1(itemp)
         i2 = idefd2(itemp)
         ixz1 = 3
c
c***** now the loop over batches for the j position
c
 40      do 70 jx = 1 , ix
            ic1 = ic1 + 1
            if (jx.le.nsspl) then
               jtemp = jx
               j1 = idefs1(jtemp)
               j2 = idefs2(jtemp)
               jxz1 = 1
               go to 50
            end if
            if (jx.le.nspl) then
               jtemp = jx - nsspl
               j1 = idefp1(jtemp)
               j2 = idefp2(jtemp)
               jxz1 = 2
               go to 50
            end if
            jtemp = jx - nspl
            j1 = idefd1(jtemp)
            j2 = idefd2(jtemp)
            jxz1 = 3
 50         ext1 = .false.
            if (ix.eq.jx) ext1 = .true.
c
            if (ixz1.eq.1) then
               iparty = ipss
c              buff = 'ss'
               go to 60
            end if
            if (ixz1.eq.2) then
               if (jxz1.eq.1) then
                  iparty = ipps
c                 buff = 'ps'
               else
                  iparty = ippp
c                 buff = 'pp'
               end if
               go to 60
            end if
            if (jxz1.eq.1) then
               iparty = ipds
c              buff = 'ds'
            end if
            if (jxz1.eq.2) then
               iparty = ipdp
c              buff = 'dp'
            end if
            if (jxz1.eq.3) then
               iparty = ipdd
c              buff = 'dd'
            end if
c
c***** now call to store ij parameters
c
 60         call memlc2(k00,k01,k02,k03,k04,k05,k06,k07,k0a,iparty)
            ifreem = imense + k06 - k0a
c
c*****  kinetic + overlap terms
c
            call cidone(d(k06),ifreem,ixz1,jxz1,ninty1,kinup,koveru,
     +                  kindow,koverd)
c
            call savone(d(k06),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04)
     +                  ,d(k05),iparty,n11,i1,i2,j1,j2,ext1)
c
            call oneda(d(k06),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04),
     +                 n11,dens,rdens,grad(1,1),grad(1,4),ext1,i1,i2,j1,
     +                 j2,ninty1,kinup,koveru,kindow,koverd)
c
c**** nuclear terms (i|a|j)
c
            ijxz1 = ixz1 + jxz1 - 1
            k0x = k06 + (ijxz1+1)**2
            call cidonb(d(k06),d(k0x),ifreem,ixz1,jxz1,ijxz1,ijxz1)
c
            call savonb(d(k06),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04)
     +                  ,d(k05),iparty,n11,i1,i2,j1,j2,ext1)
c
            call onedb(d(k06),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04),
     +                 n11,dens,grad(1,4),grad(1,7),ext1,i1,i2,j1,
     +                 j2)
 70      continue
 80   continue
c
      if (gdbg(8)) then
         write (iwr,6060)
         write (iwr,6020)
         do 90 i = 1 , na
            write (iwr,6010) i , grad(i,1) , grad(i,2) , grad(i,3)
 90      continue
c
         write (iwr,6050)
         write (iwr,6020)
         do 100 i = 1 , na
            write (iwr,6010) i , grad(i,7) , grad(i,8) , grad(i,9)
 100     continue
      end if
c
      do 120 icart = 7 , 9
         do 110 i = 1 , na
            grad(i,icart-6) = grad(i,icart-6) + grad(i,icart)
 110     continue
 120  continue
c
c********final nuclear repulsion terms
      call nucder(grad,d(k0a),d(k0a+na*na))
c
      if (gdbg(8)) then
         write (iwr,6040)
         write (iwr,6020)
         do 130 i = 1 , na
            write (iwr,6010) i , grad(i,4) , grad(i,5) , grad(i,6)
 130     continue
      end if
c
      do 150 icart = 4 , 6
         do 140 i = 1 , na
            grad(i,icart-3) = grad(i,icart-3) + grad(i,icart)
 140     continue
 150  continue
c
      if (gdbg(8)) then
         write (iwr,6030)
         write (iwr,6020)
         do 160 i = 1 , na
            write (iwr,6010) i , grad(i,1) , grad(i,2) , grad(i,3)
 160     continue
      end if
c     free allocated memory
      call gmem_free_inf(k0a,fnm,snm,'k0a')
      return
c
 6010 format (1x,i3,3(4x,d16.8))
 6020 format (1x,'atom',11x,'x',19x,'y',19x,'z')
 6030 format (//,'  entire 1 electron contribution to gradient'/)
 6040 format (//,'        other 1-electron contributions to gradient'/)
 6050 format (//,'       hellman feynman contribution to gradient'/)
 6060 format (//,'              overlap contribution to gradient'/)
      end
      subroutine defgrd(d,codens,dens,grad,iwr)
c
c***** routine coordinates the calculation of the two electron
c***** contribution to the gradient
c
      implicit REAL  (a-h,o-z)
_IF(parallel)
      logical oipsci
_ENDIF
      character *2 buff
      character *3 ydnam
c
c***** maxat is max number of atoms,maxorb is max num of basis
c***** functions,mxprim is max num of distinct primitives (ie px,py
c***** and pz count only once)
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
      dimension d(*)
      dimension ydnam(3)
      common/igarg/igarm(40)
      common/tempit/zitest,zipass,zjtest,zjpass
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/scra  /iso(mxshel,48),nt
c
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
      dimension grad(na,9),dens(nb,nb),codens(nb,nb)
     &,itesty(3,3,3,3)
      logical glog1,glog2,glog3,glog5,glog6,gdbg,ext4
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10)
c*******zipftd contains the pre-factor test data and
c*******zijkld contains the ijkl test data whilst
c*******zjpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ixz,jxz,kxz,lxz,mxz,ijxz,klxz,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
      logical ext1,ext2,ext3
c     For memory allocation tracking
      character *8 fnm,snm
      data fnm,snm/"direct.m","defgrd"/
c
      data ydnam /'e/x','e/y','e/z'/
c
c*****  a few comments about the memory allocation. most of the
c*****  large memory is contained in agarm. the first part of this
c*****  array is taken up with arrays of size nb*nb such as the
c*****  fock matrix. these arrays should end at igarm(12). the rest
c*****  is free memory to evaluate integrals in. the first part of
c*****  this is to store the parameters in.
c     sum = 0.0d0
c
      ipss = ipars*ipars
      ipps = iparp*ipars
      ippp = iparp*iparp
      ipds = ipard*ipars
      ipdp = ipard*iparp
      ipdd = ipard*ipard
c
c**** indexing an array to tell us which total to add the
c**** integral test data to
c
      icount = 0
      do 50 i = 1 , 3
         do 40 j = 1 , i
            do 30 k = 1 , i
               lmax = k
               if (k.eq.i) lmax = j
               do 20 l = 1 , lmax
                  icount = icount + 1
                  itesty(i,j,k,l) = icount
                  itesty(j,i,k,l) = icount
                  itesty(i,j,l,k) = icount
                  itesty(j,i,l,k) = icount
                  itesty(k,l,i,j) = icount
                  itesty(k,l,j,i) = icount
                  itesty(l,k,i,j) = icount
                  itesty(l,k,j,i) = icount
 20            continue
 30         continue
 40      continue
 50   continue
c
_IF(parallel)
c***   ***node-MPP***
      idum = iipsci()
c***   ***node-MPP***
_ENDIF
      if (nd.eq.0) ndspl = 0
      zitest = 0
      zipass = 0
      zjtest = 0
      zjpass = 0
      nspl = nsspl + npspl
      nspdl = nspl + ndspl
*!!!!!      k0a=igarm(12)
c
c     allocate available memory
c
      k0a = igmem_alloc_all_inf(imense,fnm,snm,'k0a',IGMEM_DEBUG)

*!!!!      if ((nb*nb).lt.(na*9)) k0a=k0a+na*9-nb*nb
      do 70 icart = 4 , 9
         do 60 i = 1 , na
            grad(i,icart) = 0.0d0
 60      continue
 70   continue
c
c***** now the loop over batches for the i position
c
      ic1 = 0
      do 170 ix = 1 , nspdl
         if (ix.le.nsspl) then
            itemp = ix
            i1 = idefs1(itemp)
            i2 = idefs2(itemp)
            ixz1 = 1
            go to 80
         end if
         if (ix.le.nspl) then
            itemp = ix - nsspl
            i1 = idefp1(itemp)
            i2 = idefp2(itemp)
            ixz1 = 2
            go to 80
         end if
         itemp = ix - nspl
         i1 = idefd1(itemp)
         i2 = idefd2(itemp)
         ixz1 = 3
c
c***** now the loop over batches for the j position
c
 80      do 160 jx = 1 , ix
            ic1 = ic1 + 1
            if (jx.le.nsspl) then
               jtemp = jx
               j1 = idefs1(jtemp)
               j2 = idefs2(jtemp)
               jxz1 = 1
               go to 90
            end if
            if (jx.le.nspl) then
               jtemp = jx - nsspl
               j1 = idefp1(jtemp)
               j2 = idefp2(jtemp)
               jxz1 = 2
               go to 90
            end if
            jtemp = jx - nspl
            j1 = idefd1(jtemp)
            j2 = idefd2(jtemp)
            jxz1 = 3
 90         ext1 = .false.
            if (ix.eq.jx) ext1 = .true.
c
            if (ixz1.eq.1) then
               iparty = ipss
               buff = 'ss'
               go to 100
            end if
            if (ixz1.eq.2) then
               if (jxz1.eq.1) then
                  iparty = ipps
                  buff = 'ps'
               else
                  iparty = ippp
                  buff = 'pp'
               end if
               go to 100
            end if
            if (jxz1.eq.1) then
               iparty = ipds
               buff = 'ds'
            end if
            if (jxz1.eq.2) then
               iparty = ipdp
               buff = 'dp'
            end if
            if (jxz1.eq.3) then
               iparty = ipdd
               buff = 'dd'
            end if
c
c***** now call to store ij parameters
c
 100        call memlc1(k00,k01,k02,k03,k04,k05,k06,k0a,iparty)
            length = k06 - k0a
            ifree1 = imense - length
            ext4 = .false.
            call sav6(d(k05),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04),
     +                d(k05),d(k06),iparty,n11,i1,i2,j1,j2,ext1,ext1,
     +                ext4,buff,ifree1)
c
c*****loop over the k batches
c
            jc1 = 0
            do 150 kx = 1 , nspdl
               if (kx.le.nsspl) then
                  ktemp = kx
                  k1 = idefs1(ktemp)
                  k2 = idefs2(ktemp)
                  kxz1 = 1
                  go to 110
               end if
               if (kx.le.nspl) then
                  ktemp = kx - nsspl
                  k1 = idefp1(ktemp)
                  k2 = idefp2(ktemp)
                  kxz1 = 2
                  go to 110
               end if
               ktemp = kx - nspl
               k1 = idefd1(ktemp)
               k2 = idefd2(ktemp)
               kxz1 = 3
c
c***** now the loop over batches for the j position
c
 110           do 140 lx = 1 , kx
                  jc1 = jc1 + 1
                  if (jc1.gt.ic1) go to 140
                  if (lx.le.nsspl) then
                     ltemp = lx
                     l1 = idefs1(ltemp)
                     l2 = idefs2(ltemp)
                     lxz1 = 1
                     go to 120
                  end if
                  if (lx.le.nspl) then
                     ltemp = lx - nsspl
                     l1 = idefp1(ltemp)
                     l2 = idefp2(ltemp)
                     lxz1 = 2
                     go to 120
                  end if
                  ltemp = lx - nspl
                  l1 = idefd1(ltemp)
                  l2 = idefd2(ltemp)
                  lxz1 = 3
_IF(parallel)
c***   ***node-MPP***
 120              continue
                  if( oipsci()) go to 140
                  ext2 = .false.
c***   ***node-MPP***
_ELSE
 120              ext2 = .false.
_ENDIF
                  if (kx.eq.lx) ext2 = .true.
                  ext3 = .false.
                  if (ic1.eq.jc1) ext3 = .true.
c
                  if (kxz1.eq.1) then
                     iparty = ipss
                     buff = 'ss'
                     go to 130
                  end if
                  if (kxz1.eq.2) then
                     if (lxz1.eq.1) then
                        iparty = ipps
                        buff = 'ps'
                     else
                        iparty = ippp
                        buff = 'pp'
                     end if
                     go to 130
                  end if
                  if (lxz1.eq.1) then
                     iparty = ipds
                     buff = 'ds'
                  end if
                  if (lxz1.eq.2) then
                     iparty = ipdp
                     buff = 'dp'
                  end if
                  if (lxz1.eq.3) then
                     iparty = ipdd
                     buff = 'dd'
                  end if
c
c***** now call to store kl parameters
c
 130              k06 = k05 + iparty
                  call memlc2(l00,l01,l02,l03,l04,l05,l06,l07,k06,
     +                        iparty)
                  l08 = l06 + iparty
                  l09 = l08 + iparty
c
                  ifreem = imense - (l09 - k0a)
c
                  ijxz1 = ixz1 + jxz1 - 1
                  klxz1 = kxz1 + lxz1 - 1
                  mxz1 = ijxz1 + klxz1 - 1
                  l0x = l06 + (ijxz+1)*(klxz+1)*(mxz1+1)
                  call cidrv(d(l06),d(l0x),ifreem,ixz1,jxz1,kxz1,lxz1,
     +                       ijxz1,klxz1,mxz1)
c
                  ext4 = .true.
                  call savdrv(d(k05),d(l00),d(l01),d(l02),d(l03),d(l04),
     +                        d(l05),d(l06),d(l07),iparty,n12,k1,k2,l1,
     +                        l2,ext2,ext2,ext4,buff,ifreem,d(k06))
c
                  itempy = itesty(ixz1,jxz1,kxz1,lxz1)
                  zitest = 0.0d0
                  zipass = 0.0d0
                  zjtest = 0.0d0
                  zjpass = 0.0d0
c***** calculation the derivative integrals
                  call bodrv(d(k05),d(l09),d(l09),d(k0a),d(k00),d(k01),
     +                       d(k02),d(k03),d(k04),n11,d(l00),d(l01),
     +                       d(l02),d(l03),d(l04),d(l05),d(l06),n12,
     +                       dens,codens,d(k06),grad(1,4),grad(1,7),
     +                       d(l08),ext1,ext2,ext3,i1,i2,j1,j2,k1,l1,
     +                       l2,prefac(itempy),xinner(itempy))
                  zipftd(itempy,1) = zipftd(itempy,1) + zjtest
                  zipftd(itempy,2) = zipftd(itempy,2) + zjpass
                  zijkld(itempy,1) = zijkld(itempy,1) + zitest
                  zijkld(itempy,2) = zijkld(itempy,2) + zipass
c
 140           continue
 150        continue
 160     continue
 170  continue
_IF(parallel)
c
c***   ***node-MPP***
c...   add 2-e components and send it around
         nbsq = na*3
         call pg_dgop(700,grad(1,4),nbsq,'+')
c***   ***node-MPP***
c
_ENDIF
c**** if symmetry then symmetrise the gradient
c
      if (nt.gt.1) then
         k00 = k0a + 144*3
         call symdrc(grad(1,4),d(k00))
      end if
c
      do 190 icart = 1 , 3
         do 180 i = 1 , na
            grad(i,icart) = grad(i,icart) + grad(i,icart+3)
 180     continue
 190  continue
c
c
      if (gdbg(8)) then
         write (iwr,6030)
         write (iwr,6020)
         do 200 i = 1 , na
            write (iwr,6010) i , grad(i,4) , grad(i,5) , grad(i,6)
 200     continue
      end if
c
      write (iwr,6040)
      mmax = 0
 220  mmin = mmax + 1
      mmax = mmax + 8
      if (mmax.gt.na) mmax = na
      write (iwr,6070)
      write (iwr,6050) (i,i=mmin,mmax)
      write (iwr,6070)
      do 230 n = 1 , 3
         write (iwr,6060) ydnam(n) , (grad(i,n),i=mmin,mmax)
 230  continue
      if (mmax.lt.na) go to 220
c
c     free allocated memory
c
      call gmem_free_inf(k0a,fnm,snm,'k0a')
      return
c
 6010 format (1x,i3,3(4x,d16.8))
 6020 format (1x,'atom',11x,'x',19x,'y',19x,'z')
 6070 format (/)
 6030 format (//,'             2 electron contribution to gradient'/)
 6040 format (/1x,22('=')/1x,'gradient of the energy'/1x,22('='))
 6050 format (1x,'atom',8(6x,i2,7x))
 6060 format (3x,a3,8f15.7)
      end
      subroutine defmp2(d,orbmo,e,nocc,nvir,iwr)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c**** routine coordinates the calculation of the mp2 energy
c
c**** routine has a few flaws- most notably, the calculation of paramete
c**** arrays for the kl primitive pairing takes a significant percentage
c**** of the time since it is repeated for every combination of i and j
c**** shells.
c
      implicit REAL  (a-h,o-z)
      character *2 buffkl,buff
      character *4 char
c
c***** maxat is max number of atoms,maxorb is max num of basis
c***** functions,mxprim is max num of distinct primitives (ie px,py
c***** and pz count only once) ipar2 is the max num of primives on
c***** one index that can be dealt with in an integral batch
c***** ipar3 is max num of basis functions in the batch and ipar
c***** is just ipar2**2
c
c*****  ipar3 must be less than 65 or some shortloop directives
c*****  on cray machines will go haywire
c
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      dimension d(*)
INCLUDE(common/prints)
      common/igarg/igarm(40), ipar,ip2,ip3
      common/tempit/zitest,zipass,zjtest,zjpass
      common/save/info(5,nbrkmx),tmax(nbrkmx)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nnnnt,nt3,
     1itable(8,8),mult(8,8),irr(maxorb),imosb(8,maxorb),irrb(maxorb)
     2,ittt2(8,maxorb),ittt(8,maxorb),ibb(maxorb),ins(8),iss(8)
     3,ivsf(2,8)
      common/scra  /iso(mxshel,48),nt,nauq
*!!!!! nauq read by countd calls
*     &,index,naxi,nauq,nwauq(maxat)
c
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      dimension orbmo(nb,nb),tmaxkl(50)
     &,occmax(maxorb),virmax(maxorb),occim(maxorb),nmemt(maxorb)
     &,itesty(3,3,3,3)
      logical glog1,glog2,glog3,glog5
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3
c*******ipftd contains the pre-factor test data and
c*******ijkld contains the ijkl test data whilst
c*******jpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ixz,jxz,kxz,lxz,mxz,ijxz,klxz,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
      dimension char(21)
      logical ext1,ext2,ext3,ext4
      character*10 charwall
c     For memory allocation tracking
      character *8 fnm,snm
      data fnm,snm/"direct.m","defmp2"/
c
      data char/'ssss','psss','psps','ppss','ppps','pppp'
     &,'dsss','dsps','dspp','dsds','dpss','dpps','dppp'
     &,'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
c
c
      do 30 j = 1 , 2
         do 20 i = 1 , 21
            zipftd(i,j) = 1.0d0
            zijkld(i,j) = 1.0d0
 20      continue
 30   continue
c
c**** indexing an array to tell us which total to add the
c**** integral test data to
c
      icount = 0
      do 70 i = 1 , 3
         do 60 j = 1 , i
            do 50 k = 1 , i
               lmax = k
               if (k.eq.i) lmax = j
               do 40 l = 1 , lmax
                  icount = icount + 1
                  itesty(i,j,k,l) = icount
                  itesty(j,i,k,l) = icount
                  itesty(i,j,l,k) = icount
                  itesty(j,i,l,k) = icount
                  itesty(k,l,i,j) = icount
                  itesty(k,l,j,i) = icount
                  itesty(l,k,i,j) = icount
                  itesty(l,k,j,i) = icount
 40            continue
 50         continue
 60      continue
 70   continue
c
c*****  a few comments about the memory allocation. most of the
c*****  large memory is contained in agarm. the first part of this
c*****  array is taken up with arrays of size nb*nb such as the
c*****  fock matrix. these arrays should end at igarm(2). the rest
c*****  is free memory to evaluate integrals in. the first part of
c*****  this is to store the parameters in.
c
      mp2tar = igmem_alloc_all_inf(imense,fnm,snm,
     +     'mp2tar',IGMEM_DEBUG)

c
c**** add onto m0a the storage af all exchange integrals
      nprims = nsp + npp + ndp
      m0a = mp2tar
      np2 =  nprims*(nprims+1)/2
      m0a = m0a + np2
c
      maxipq = 9
c
c***** the storage requirement for mp2
c
      nsizb = 0
      nsizeb = 0
      do 100 i = 1 , nocc
         do 90 j = 1 , i
            ijsym = mult(irr(i),irr(j))
            do 80 ia = 1 , nvir
               ijasym = mult(ijsym,irr(nocc+ia))
               nsizb = nsizb + ins(ijasym)
 80         continue
 90      continue
         nmemt(i) = nsizb
         nsizeb = nsizeb + nsizb
         nsizb = 0
 100  continue
c
      tim0 = cpulft(1)
      write (iwr,6000)
      write (iwr,6170) tim0 ,charwall()
      write (iwr,6160) nb , nocc , nvir
      if (oprint(55)) then
       write (iwr,6180)
       write (iwr,6070) nsizeb
       write (iwr,6150) nb*nvir*nocc*(nocc+1)/2
      endif
      if (nd.gt.0) maxipq = 36
c
      million = 1000000
      ileft = imense - np2 - nsizeb
      if (ileft.lt.million) then
         ileft = imense - million 
         do 110 i = 1 , nocc
            if (ileft.lt.nmemt(i))
     +          call caserr('not enough memory for mp2 calculation')
 110     continue
         nbatch = nsizeb/ileft + 1
         ileft = nsizeb/(nbatch+1)
         if (ileft.lt.million) ileft = million
      else
         nbatch = 1
         ileft = million
      end if
      write (iwr,6200) nbatch
      write (iwr,6210) ileft
c
c***** now sorting out memory for d function code. memroy needed is only
c***** for the parameter arrays. this is arbitarily assigned to be
c***** a maximum for each array of 1/100 of the available memory.
      xavail = dsqrt((dfloat(ileft)/100.0d0))
      ipard = aint(xavail/6.0d0)
      iparp = aint(xavail/3.0d0)
      ipars = aint(xavail)
      call countd(1,idefs1,idefs2,isatbf,isatp,ipars,nsspl,nauq)
      if(oprint(55)) then
       write (iwr,6020) nsspl
       write (iwr,6060) ipars
       write (iwr,6030)
       write (iwr,6040) (idefs1(j),j=1,nsspl)
       write (iwr,6050)
       write (iwr,6040) (idefs2(j),j=1,nsspl)
       write (iwr,6080)
      endif
      if (np.eq.0) then
         iparp = 0
         npspl = 0
         go to 120
      end if
      call countd(ns+1,idefp1,idefp2,ipatbf,ipatp,iparp,npspl,nauq)
      if(oprint(55)) then
       write (iwr,6130) npspl
       write (iwr,6060) iparp
       write (iwr,6030)
       write (iwr,6040) (idefp1(j),j=1,npspl)
       write (iwr,6050)
       write (iwr,6040) (idefp2(j),j=1,npspl)
       write (iwr,6080)
      endif
 120  if (nd.eq.0) then
         ipard = 0
         ndspl = 0
         go to 130
      end if
      call countd(ns+np+1,idefd1,idefd2,idatbf,idatp,ipard,ndspl,nauq)
      if(oprint(55)) then
       write (iwr,6140) ndspl
       write (iwr,6060) ipard
       write (iwr,6030)
       write (iwr,6040) (idefd1(j),j=1,ndspl)
       write (iwr,6050)
       write (iwr,6040) (idefd2(j),j=1,ndspl)
       write (iwr,6080)
      endif
c
c
c**** sorting out splitting in batches for inner pair
c
 130  ipss = ipars*ipars
      ipps = iparp*ipars
      ippp = iparp*iparp
      ipds = ipard*ipars
      ipdp = ipard*iparp
      ipdd = ipard*ipard
c
c
c**** find the arrays needed for testing
c
      do 140 j = 1 , nb
         virmax(j) = 1.0d-15
         occmax(j) = 1.0d-15
 140  continue
c
      do 160 i = 1 , nocc
         do 150 j = 1 , nb
            occmax(j) = max(occmax(j),dabs(orbmo(j,i)))
 150     continue
 160  continue
      do 180 i = 1 , nvir
         do 170 j = 1 , nb
            virmax(j) = max(virmax(j),dabs(orbmo(j,i+nocc)))
 170     continue
 180  continue
c
      i = ns
      do 190 j = ns + 1 , ns + np
         occmax(j) = max(occmax(i+1),occmax(i+2),occmax(i+3))
         virmax(j) = max(virmax(i+1),virmax(i+2),virmax(i+3))
         i = i + 3
 190  continue
      do 200 j = ns + np + 1 , ns + np + nd
         occmax(j) = max(occmax(i+1),occmax(i+2),occmax(i+3),occmax(i+4)
     +               ,occmax(i+5),occmax(i+6))
         virmax(j) = max(virmax(i+1),virmax(i+2),virmax(i+3),virmax(i+4)
     +               ,virmax(i+5),virmax(i+6))
         i = i + 6
 200  continue
c
      do 210 i = 1 , 21
         prefac(i) = 1.0d-10
         xinner(i) = 1.0d-10
 210  continue
c
c***** multiply all mos associated with dxy,dxz,dyz orbitals by
c***** root(3)
      sqr3 = dsqrt(3.0d0)
      sqr13 = 1.0d0/sqr3
      icmp2 = 0
      do 230 i = ns + np3 + 1 , nb
         icmp2 = icmp2 + 1
         if (icmp2.eq.6) icmp2 = 0
         if (icmp2.ne.1 .and. icmp2.ne.4 .and. icmp2.ne.0) then
            do 220 j = 1 , nocc + nvir
               orbmo(i,j) = orbmo(i,j)*sqr3
 220        continue
         end if
 230  continue
c
c*********************************************
c***** the two electron integral calls********
c*********************************************
c
      zspass = 0
      zstest = 0
      if (nd.eq.0) ndspl = 0
      nspl = nsspl + npspl
      nspdl = nspl + ndspl
      e2 = 0.0d0
c
c***** now the loop over i molecular orbital batches
      nopass = 0
      nxx = 0
      do 480 iocy = 1 , nocc
         nxx = nxx + 1
         iocc = iocy - nxx
c
c*** some messing with the memory allocation
c
         memb = 0
         do 240 nxc = 1 , nxx
            memb = memb + nmemt(iocc+nxc)
 240     continue
         ione = nb*maxipq*nxx + (2+nxx)*nvir + max(ipss,ippp*9,ipdd*36)
         itwo = nvir*nvir
         mem = memb + max(ione,itwo)
c
         ifreem = imense - m0a - million
         if (mem.gt.ifreem) call caserr('error in mp2 memory')
         memp = mem + nb*maxipq + nmemt(iocc+nxx+1)
c
         if (memp.gt.ifreem .or. iocy.ge.nocc) then
c
            nopass = nopass + 1
            tim1 = cpulft(1)
            write (iwr,6080) 
            write (iwr,6090) nopass,tim1 ,charwall()
            write (iwr,6100) mem
            write (iwr,6110) e2
            call timit(3)
            call flushn(iwr)
c***** now work out the starting position in d of all these arrays
c***** b array
            imp2b = m0a
c***** e and c array
            imp2e = m0a + memb
            imp2c = imp2e
c***** d array
            imp2d = imp2c + maxipq*nb*nxx
c***** tp and tq arrays
            imp2tp = imp2d + nxx*nvir
            imp2tq = imp2tp + nvir
c***** temp array
            imp2t = imp2tq + nvir
c
c***** starting position in memory for two electron integrals
c
            k0a = m0a + mem
            ifreem = imense - k0a
c
c**** final array needed for integral testing
c
            do 250 j = 1 , nb
               occim(j) = 1.0d-15
 250        continue
c
            do 270 i = iocc + 1 , iocc + nxx
               do 260 j = 1 , nb
                  occim(j) = max(occim(j),dabs(orbmo(j,i)))
 260           continue
 270        continue
c
            i = ns
            do 280 j = ns + 1 , ns + np
               occim(j) = max(occim(i+1),occim(i+2),occim(i+3))
               i = i + 3
 280        continue
            do 290 j = ns + np + 1 , ns + np + nd
               occim(j) = max(occim(i+1),sqr13*occim(i+2),sqr13*occim(i+
     +                    3),occim(i+4),sqr13*occim(i+5),occim(i+6))
               i = i + 6
 290        continue
c
c*** find the array of exchange integrals for all primitves
c
            icount = igarm(2)
            jc1 = 0
            do 350 kx = 1 , nspdl
               if (kx.le.nsspl) then
                  ktemp = kx
                  k1 = idefs1(ktemp)
                  k2 = idefs2(ktemp)
                  kxz1 = 1
                  go to 300
               end if
               if (kx.le.nspl) then
                  ktemp = kx - nsspl
                  k1 = idefp1(ktemp)
                  k2 = idefp2(ktemp)
                  kxz1 = 2
                  go to 300
               end if
               ktemp = kx - nspl
               k1 = idefd1(ktemp)
               k2 = idefd2(ktemp)
               kxz1 = 3
 300           continue
c              kvl = nvl(kxz1)
c
c***** now the loop over batches for the l position
c
               do 340 lx = 1 , kx
                  jc1 = jc1 + 1
                  if (lx.le.nsspl) then
                     ltemp = lx
                     l1 = idefs1(ltemp)
                     l2 = idefs2(ltemp)
                     lxz1 = 1
                     go to 310
                  end if
                  if (lx.le.nspl) then
                     ltemp = lx - nsspl
                     l1 = idefp1(ltemp)
                     l2 = idefp2(ltemp)
                     lxz1 = 2
                     go to 310
                  end if
                  ltemp = lx - nspl
                  l1 = idefd1(ltemp)
                  l2 = idefd2(ltemp)
                  ispos = ns + np3 + (l1-ns-np-1)*6
                  lxz1 = 3
 310              continue
c                 lvl = nvl(lxz1)
                  ext2 = .false.
                  if (kx.eq.lx) ext2 = .true.
                  ext3 = .false.
c
                  if (kxz1.eq.1) then
                     iparty = ipss
                     buffkl = 'ss'
                     go to 320
                  end if
                  if (kxz1.eq.2) then
                     if (lxz1.eq.1) then
                        iparty = ipps
                        buffkl = 'ps'
                     else
                        iparty = ippp
                        buffkl = 'pp'
                     end if
                     go to 320
                  end if
                  if (lxz1.eq.1) then
                     iparty = ipds
                     buffkl = 'ds'
                  end if
                  if (lxz1.eq.2) then
                     iparty = ipdp
                     buffkl = 'dp'
                  end if
                  if (lxz1.eq.3) then
                     iparty = ipdd
                     buffkl = 'dd'
                  end if
c
c***** now call to store kl parameters
c
 320              call memlc2(l00,l01,l02,l03,l04,l05,l06,l07,k0a,
     +                        iparty)
                  l08 = l06 + iparty
                  ifreem = imense - l08
                  ijxz1 = ixz1 + jxz1 - 1
                  klxz1 = kxz1 + lxz1 - 1
                  mxz1 = ijxz1 + klxz1 - 1
c
                  ext4 = .true.
                  call sav9(d(k0a),d(l00),d(l01),d(icount),d(l03),d(l04)
     +                      ,d(l05),d(l06),d(l07),iparty,n12,k1,k2,l1,
     +                      l2,ext2,buffkl,ifreem,virmax,
     +                      occim)
                  tmaxkl(jc1) = d(icount)
                  do 330 itmax = icount + 1 , icount + n12
                     tmaxkl(jc1) = max(tmaxkl(jc1),d(itmax))
 330              continue
                  if (jc1.gt.50)
     +                call caserr('overwriting in tmaxkl array')
                  icount = icount + n12
c
 340           continue
 350        continue
c
c***** now the loop over basis function batches for the i position
c
            call vclr(d(imp2b),1,memb)
_IF1()c           ibuff = 5
_IF1()c           ibuffo = 1
            do 470 ix = 1 , ns + np + nd
c**** symmetry test on ix shell
               do 360 it = 2 , nt
                  ids = iso(ix,it)
                  if (ids.gt.ix) then
                     zstest = zstest + ix
                     go to 470
                  end if
 360           continue
               if (ix.le.ns) then
                  ippos = ix - 1
                  ixz1 = 1
                  go to 370
               end if
               if (ix.le.(ns+np)) then
                  ippos = ns + (ix-ns-1)*3
                  ixz1 = 2
                  go to 370
               end if
               ippos = ns + np3 + (ix-ns-np-1)*6
               ixz1 = 3
 370           ivl = nvl(ixz1)
_IF1()c      if (ix.gt.ibuff) then
_IF1()c      write(iwr,1711) ibuffo,ibuff
_IF1()c      call timit(1)
_IF1()c      call flushn(iwr)
_IF1()c      ibuffo=ix
_IF1()c      ibuff=ibuff+5
_IF1()c      endif
c
c***** now the loop over batches for the j position
c
               do 460 jx = 1 , ix
                  zstest = zstest + 1
c***** symmetry test on jx shell
                  do 380 it = 2 , nt
                     ids = iso(ix,it)
                     jds = iso(jx,it)
                     if (jds.gt.ix) go to 460
                     if (ids.lt.jds) then
                        nds = ids
                        ids = jds
                        jds = nds
                     end if
                     if (ids.eq.ix .and. jds.gt.jx) go to 460
 380              continue
                  zspass = zspass + 1
c
                  if (jx.le.ns) then
                     jxz1 = 1
                     iqpos = jx - 1
                     go to 390
                  end if
                  if (jx.le.(ns+np)) then
                     iqpos = ns + (jx-ns-1)*3
                     jxz1 = 2
                     go to 390
                  end if
                  iqpos = ns + np3 + (jx-ns-np-1)*6
                  jxz1 = 3
 390              jvl = nvl(jxz1)
                  ext1 = .false.
                  if (ix.eq.jx) ext1 = .true.
c
                  if (ixz1.eq.1) then
                     buff = 'ss'
                     go to 400
                  end if
                  if (ixz1.eq.2) then
                     if (jxz1.eq.1) then
                        buff = 'ps'
                     else
                        buff = 'pp'
                     end if
                     go to 400
                  end if
                  if (jxz1.eq.1) then
                     buff = 'ds'
                  end if
                  if (jxz1.eq.2) then
                     buff = 'dp'
                  end if
                  if (jxz1.eq.3) then
                     buff = 'dd'
                  end if
 400              iparty = 36
c
c***** now call to store ij parameters
c
                  call vclr(d(imp2c),1,ivl*jvl*nb*nxx)
                  call memlc1(k00,k01,k02,k03,k04,k05,k06,k0a,iparty)
                  ifree1 = imense - k06
                  ext4 = .false.
                  call sav9(d(k05),d(k0a),d(k00),d(k01),d(k02),d(k03),
     +                      d(k04),d(k05),d(k06),iparty,n11,ix,ix,jx,jx,
     +                      ext1,buff,ifree1,virmax,occmax)
                  if (n11.gt.iparty) then
                    write(iwr,*) ' this code only seems to cater for ',
     +                           'maximal 6 functions in contraction'
                    write(iwr,*) ' either split your contraction or ',
     +                       'change iparty = 36 to e.g. iparty = 100',
     +                       ' in ** direct.m **'
                    call caserr('internal error in defmp2')
                  end if
c
c*****loop over the k batches
c
                  intdoq = 0
                  icount = igarm(2)
                  jc1 = 0
                  do 450 kx = 1 , nspdl
                     if (kx.le.nsspl) then
                        ktemp = kx
                        k1 = idefs1(ktemp)
                        k2 = idefs2(ktemp)
                        kxz1 = 1
                        irpos = k1 - 1
                        go to 410
                     end if
                     if (kx.le.nspl) then
                        ktemp = kx - nsspl
                        k1 = idefp1(ktemp)
                        k2 = idefp2(ktemp)
                        irpos = ns + (k1-ns-1)*3
                        kxz1 = 2
                        go to 410
                     end if
                     ktemp = kx - nspl
                     k1 = idefd1(ktemp)
                     k2 = idefd2(ktemp)
                     irpos = ns + np3 + (k1-ns-np-1)*6
                     kxz1 = 3
 410                 continue
c                    kvl = nvl(kxz1)
c
c***** now the loop over batches for the j position
c
                     do 440 lx = 1 , kx
                        jc1 = jc1 + 1
                        if (lx.le.nsspl) then
                           ltemp = lx
                           l1 = idefs1(ltemp)
                           l2 = idefs2(ltemp)
                           ispos = l1 - 1
                           lxz1 = 1
                           go to 420
                        end if
                        if (lx.le.nspl) then
                           ltemp = lx - nsspl
                           l1 = idefp1(ltemp)
                           l2 = idefp2(ltemp)
                           ispos = ns + (l1-ns-1)*3
                           lxz1 = 2
                           go to 420
                        end if
                        ltemp = lx - nspl
                        l1 = idefd1(ltemp)
                        l2 = idefd2(ltemp)
                        ispos = ns + np3 + (l1-ns-np-1)*6
                        lxz1 = 3
 420                    continue
c                       lvl = nvl(lxz1)
                        ext2 = .false.
                        if (kx.eq.lx) ext2 = .true.
                        ext3 = .false.
c
                        if (kxz1.eq.1) then
                           iparty = ipss
                           buffkl = 'ss'
                           go to 430
                        end if
                        if (kxz1.eq.2) then
                           if (lxz1.eq.1) then
                              iparty = ipps
                              buffkl = 'ps'
                           else
                              iparty = ippp
                              buffkl = 'pp'
                           end if
                           go to 430
                        end if
                        if (lxz1.eq.1) then
                           iparty = ipds
                           buffkl = 'ds'
                        end if
                        if (lxz1.eq.2) then
                           iparty = ipdp
                           buffkl = 'dp'
                        end if
                        if (lxz1.eq.3) then
                           iparty = ipdd
                           buffkl = 'dd'
                        end if
c
c***** now call to store kl parameters
c
 430                    call memlc2(l00,l01,l02,l03,l04,l05,l06,l07,k05,
     +                              iparty)
                        l08 = l06 + iparty
                        ifreem = imense - l08
                        ijxz1 = ixz1 + jxz1 - 1
                        klxz1 = kxz1 + lxz1 - 1
                        mxz1 = ijxz1 + klxz1 - 1
                        call tizer(d(l06),ifreem,ixz1,jxz1,kxz1,lxz1,
     +                             ijxz1,klxz1,mxz1)
c
                        ext4 = .true.
                        call sav8(d(k05),d(l00),d(l01),d(l02),d(l03),
     +                            d(l04),d(l05),d(l06),iparty,
     +                            n12,k1,k2,l1,l2,ext2,ext4)
c
                        tmax(1) = tmaxkl(jc1)
                        zipass = 0
                        zjpass = 0
                        zitest = 0
                        zjtest = 0
                        itempy = itesty(ixz1,jxz1,kxz1,lxz1)
c*****  call integral routine
                        call bodmp2(d(k05),d(l08),d(l08),d(k0a),d(k00),
     +                              d(k01),d(k02),d(k03),d(k04),n11,
     +                              d(l00),d(l01),d(icount),d(l03),
     +                              d(l04),d(l05),d(l06),n12,orbmo,
     +                              iocc,nxx,irpos,ispos,ext1,ext2,ext3,
     +                              ix,ix,jx,jx,k1,k2,l1,l2,
     +                              prefac(itempy),xinner(itempy),
     +                              d(imp2c),maxipq,d(imp2t),intdoq,iwr)
                        zipftd(itempy,1) = zipftd(itempy,1) + zjtest
                        zipftd(itempy,2) = zipftd(itempy,2) + zjpass
                        zijkld(itempy,1) = zijkld(itempy,1) + zitest
                        zijkld(itempy,2) = zijkld(itempy,2) + zipass
                        icount = icount + n12
c
 440                 continue
 450              continue
                  if (intdoq.ne.0) then
                     imax = 1
                     jmax0 = 1
c*****  at this point the first quarter of 4 index should have been done
c*****  now do second and third part of the transformation
                     call indxbc(ivl,jvl,imax,jmax0,ext1,orbmo,nocc,
     +                           iocc,nxx,ippos,iqpos,d(imp2c),maxipq,
     +                           d(imp2d),nvir,d(imp2b))
                  end if
 460           continue
 470        continue
c*****  do final part of 4 index transformation and then calculate
c*****  the mp2 energy
            call indx4d(orbmo,e,nocc,iocc,nxx,e2,d(imp2b),nvir,d(imp2e))
            nxx = 0
         end if
 480  continue
      if (nt.gt.1.and.oprint(55) ) then
         write (iwr,6080)
         write (iwr,*) ' symmetry test data on outer pair of shells'
         dperc = 100.0d0*zspass/zstest
         write (iwr,6190) dperc
      end if
c
      if(oprint(55)) then
       write (iwr,6080)
       write (iwr,*) ' **************************************'
       write (iwr,*) ' integral test data over all iterations'
       write (iwr,*) ' **************************************'
       write (iwr,6080)
       write (iwr,6120)
       imax = 21
       if (nd.eq.0) imax = 6
       do 490 i = 1 , imax
          dperc1 = zipftd(i,2)/zipftd(i,1)
          dperc2 = zijkld(i,2)/zijkld(i,1)
          dperc3 = dperc1*dperc2
          dperc1 = dperc1*100.0d0
          dperc2 = dperc2*100.0d0
          dperc3 = dperc3*100.0d0
          write (iwr,6010) char(i) , dperc1 , dperc2 , dperc3
 490   continue
      endif
      write (iwr,6005) e2
      call timit(3)
      tim1 = cpulft(1)
      write(iwr,6001) tim1 ,charwall()
c
c     free memory
c
      call gmem_free_inf(mp2tar,fnm,snm,'mp2tar')

c
      return
c
 6000 format (/10x,22('*')/10x,'direct-mp2 calculation'/10x,22('*'))
 6001 format (/1x,'end of direct-mp2 calculation at ',f8.2,
     + ' seconds',a10,' wall'/)
 6005 format (/10x,42('*')
     +        /10x,'mp2 correlation energy     ',f15.8
     +        /10x,42('*'))
 6010 format (2x,a4,10x,f6.2,16x,f6.2,14x,f6.2)
 6020 format (/1x,
     +   'number of batches of s shells                    ',i5)
 6030 format (' starting point(s) of each batch')
 6040 format (1x,12i5)
 6050 format (' ending point(s) of each batch')
 6060 format (1x,
     + 'maximum number of primitive shells in each batch ', i5)
 6070 format (1x,'with symmetry    ',i20)
 6080 format (/)
 6090 format (1x,'commence mp2 pass number',i4,' at ',f8.2,' seconds'
     +        ,a10,' wall')
 6100 format (1x,'memory used in this pass    ',i15)
 6110 format (1x,'mp2 energy before this pass ',f15.8)
_IF1()c1711  format(1x,'p shells',i4,' to',i4,' have been completed in
_IF1()c      transformation (pq|rs) to (ia|jb)')
 6120 format (1x,'routine',5x,'prefactor  pass %',6x,'inner pass %',8x,
     +        'total pass %')
 6130 format (1x,'number of batches of p shells       ',i5)
 6140 format (1x,'number of batches of d shells       ',i5)
 6150 format (1x,'without symmetry ',i20//)
 6160 format (1x,'number of basis functions           ',i6/
     +        1x,'number of active occupied orbitals  ',i6/
     +        1x,'number of active virtual orbitals   ',i6/)
 6170 format (/1x,'commencing direct mp2 calculation at ',
     +        f8.2,' seconds',a10,' wall'/)
 6180 format (1x,
     +   'memory needed by big array in 1 batch, 4 index transformation'
     +   )
 6190 format (1x,'percentage of pairs passing test ',f6.2)
 6200 format (1x,'expected number of batches  ',i5)
 6210 format (1x,'words of memory available for rest of transformation',
     +        i11)
      end
      subroutine defunc(d,fock,dens,codens)
c
c**** routine coordinates the calculation of the two electron integrals
c**** and their addition to the fock matrix
c
c**** fock = fock matrix
c**** dens = density matrix
c**** codens = "contracted" (over shells) density matrix used in the
c****                                                    integral tests
c
      implicit REAL  (a-h,o-z)
      character *2 buff
c
c***** maxat is max number of atoms,maxorb is max num of basis
c***** functions,mxprim is max num of distinct primitives (ie px,py
c***** and pz count only once)
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
      dimension d(*)
      common/igarg/igarm(40)
      common/tempit/zitest,zipass,zjtest,zjpass
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
      dimension fock(nb,nb),dens(nb,nb),codens(nb,nb)
     &,itesty(3,3,3,3)
_IF(parallel)
c***   ***node-MPP***
      common/pardir/next,nbcnt(256)
c***   ***node-MPP***
_ENDIF
      logical glog1,glog2,glog3,glog5
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3
c*******zipftd contains the pre-factor test data and
c*******zijkld contains the ijkl test data whilst
c*******zjpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ixz,jxz,kxz,lxz,mxz,ijxz,klxz,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
      logical ext1,ext2,ext3,ext4
c
c     For tracking memory allocation
      character*8 fnm,snm
      data fnm,snm/"direct.m","defunc"/
c
c
c*****  a few comments about the memory allocation. most of the
c*****  large memory is contained in agarm. the first part of this
c*****  array is taken up with arrays of size nb*nb such as the
c*****  fock matrix. these arrays should end at igarm(12). the rest
c*****  is free memory to evaluate integrals in. the first part of
c*****  this is to store the parameters in.
c     sum = 0.0d0
c
      ipss = ipars*ipars
      ipps = iparp*ipars
      ippp = iparp*iparp
      ipds = ipard*ipars
      ipdp = ipard*iparp
      ipdd = ipard*ipard
c
c*********************************************
c***** the two electron integral calls********
c*********************************************
c
c**** indexing an array to tell us which total to add the
c**** integral test data to
c
      icount = 0
      do 50 i = 1 , 3
         do 40 j = 1 , i
            do 30 k = 1 , i
               lmax = k
               if (k.eq.i) lmax = j
               do 20 l = 1 , lmax
                  icount = icount + 1
                  itesty(i,j,k,l) = icount
                  itesty(j,i,k,l) = icount
                  itesty(i,j,l,k) = icount
                  itesty(j,i,l,k) = icount
                  itesty(k,l,i,j) = icount
                  itesty(k,l,j,i) = icount
                  itesty(l,k,i,j) = icount
                  itesty(l,k,j,i) = icount
 20            continue
 30         continue
 40      continue
 50   continue
c
_IF(parallel)
c***   **MPP**
       call pg_dlbreset
       next = ipg_dlbtask()
c***   **MPP**
_ENDIF(parallel)
      if (nd.eq.0) ndspl = 0
      zitest = 0
      zipass = 0
      zjtest = 0
      zjpass = 0
      nspl = nsspl + npspl
      nspdl = nspl + ndspl
c
c *** explicitly grab all the memory
c
      k0a = igmem_alloc_all_inf(ifree,fnm,snm,'k0a',IGMEM_DEBUG)

c
c***** now the loop over batches for the i position
c
      ic1 = 0
      do 150 ix = 1 , nspdl
         if (ix.le.nsspl) then
            itemp = ix
            i1 = idefs1(itemp)
            i2 = idefs2(itemp)
            ixz1 = 1
            go to 60
         end if
         if (ix.le.nspl) then
            itemp = ix - nsspl
            i1 = idefp1(itemp)
            i2 = idefp2(itemp)
            ixz1 = 2
            go to 60
         end if
         itemp = ix - nspl
         i1 = idefd1(itemp)
         i2 = idefd2(itemp)
         ixz1 = 3
c
c***** now the loop over batches for the j position
c
 60      do 140 jx = 1 , ix
            ic1 = ic1 + 1
            if (jx.le.nsspl) then
               jtemp = jx
               j1 = idefs1(jtemp)
               j2 = idefs2(jtemp)
               jxz1 = 1
               go to 70
            end if
            if (jx.le.nspl) then
               jtemp = jx - nsspl
               j1 = idefp1(jtemp)
               j2 = idefp2(jtemp)
               jxz1 = 2
               go to 70
            end if
            jtemp = jx - nspl
            j1 = idefd1(jtemp)
            j2 = idefd2(jtemp)
            jxz1 = 3
 70         ext1 = .false.
            if (ix.eq.jx) ext1 = .true.
c
            if (ixz1.eq.1) then
               iparty = ipss
               buff = 'ss'
               go to 80
            end if
            if (ixz1.eq.2) then
               if (jxz1.eq.1) then
                  iparty = ipps
                  buff = 'ps'
               else
                  iparty = ippp
                  buff = 'pp'
               end if
               go to 80
            end if
            if (jxz1.eq.1) then
               iparty = ipds
               buff = 'ds'
            end if
            if (jxz1.eq.2) then
               iparty = ipdp
               buff = 'dp'
            end if
            if (jxz1.eq.3) then
               iparty = ipdd
               buff = 'dd'
            end if
c
c***** now call to store ij parameters
c
 80         call memlc1(k00,k01,k02,k03,k04,k05,k06,k0a,iparty)
            ifree1 = ifree - (k06 - k0a)
            ext4 = .false.
            call sav6(d(k05),d(k0a),d(k00),d(k01),d(k02),d(k03),d(k04),
     +                d(k05),d(k06),iparty,n11,i1,i2,j1,j2,ext1,ext1,
     +                ext4,buff,ifree1)
c
c*****loop over the k batches
c
            jc1 = 0
            do 130 kx = 1 , nspdl
               if (kx.le.nsspl) then
                  ktemp = kx
                  k1 = idefs1(ktemp)
                  k2 = idefs2(ktemp)
                  kxz1 = 1
                  go to 90
               end if
               if (kx.le.nspl) then
                  ktemp = kx - nsspl
                  k1 = idefp1(ktemp)
                  k2 = idefp2(ktemp)
                  kxz1 = 2
                  go to 90
               end if
               ktemp = kx - nspl
               k1 = idefd1(ktemp)
               k2 = idefd2(ktemp)
               kxz1 = 3
c
c***** now the loop over batches for the j position
c
 90            do 120 lx = 1 , kx
                  jc1 = jc1 + 1
                  if (jc1.gt.ic1) go to 120
                  if (lx.le.nsspl) then
                     ltemp = lx
                     l1 = idefs1(ltemp)
                     l2 = idefs2(ltemp)
                     lxz1 = 1
                     go to 100
                  end if
                  if (lx.le.nspl) then
                     ltemp = lx - nsspl
                     l1 = idefp1(ltemp)
                     l2 = idefp2(ltemp)
                     lxz1 = 2
                     go to 100
                  end if
                  ltemp = lx - nspl
                  l1 = idefd1(ltemp)
                  l2 = idefd2(ltemp)
                  lxz1 = 3
 100              ext2 = .false.
                  if (kx.eq.lx) ext2 = .true.
                  ext3 = .false.
                  if (ic1.eq.jc1) ext3 = .true.
c
                  if (kxz1.eq.1) then
                     iparty = ipss
                     buff = 'ss'
                     go to 110
                  end if
                  if (kxz1.eq.2) then
                     if (lxz1.eq.1) then
                        iparty = ipps
                        buff = 'ps'
                     else
                        iparty = ippp
                        buff = 'pp'
                     end if
                     go to 110
                  end if
                  if (lxz1.eq.1) then
                     iparty = ipds
                     buff = 'ds'
                  end if
                  if (lxz1.eq.2) then
                     iparty = ipdp
                     buff = 'dp'
                  end if
                  if (lxz1.eq.3) then
                     iparty = ipdd
                     buff = 'dd'
                  end if
c
c***** now call to store kl parameters
c
 110              continue
                  call memlc2(l00,l01,l02,l03,l04,l05,l06,l07,k05,
     +                        iparty)
                  l08 = l06 + iparty
                  ifreem = ifree - (l08 - k0a)
                  ijxz1 = ixz1 + jxz1 - 1
                  klxz1 = kxz1 + lxz1 - 1
                  mxz1 = ijxz1 + klxz1 - 1
                  call tizer(d(l06),ifreem,ixz1,jxz1,kxz1,lxz1,ijxz1,
     +                       klxz1,mxz1)
c
                  ext4 = .true.
                  call sav6(d(k05),d(l00),d(l01),d(l02),d(l03),d(l04),
     +                      d(l05),d(l06),d(l07),iparty,n12,k1,k2,l1,l2,
     +                      ext2,ext2,ext4,buff,ifreem)
c
                  itempy = itesty(ixz1,jxz1,kxz1,lxz1)
                  zitest = 0.0d0
                  zipass = 0.0d0
                  zjtest = 0.0d0
                  zjpass = 0.0d0
c*****  calculates 2 electron integrals and adds them to fock matrix
                  call body(d(k05),d(l08),d(l08),d(k0a),d(k00),d(k01),
     +                      d(k02),d(k03),d(k04),n11,d(l00),d(l01),
     +                      d(l02),d(l03),d(l04),d(l05),d(l06),n12,fock,
     +                      dens,codens,ext1,ext2,ext3,i1,i2,j1,j2,k1,
     +                      l1,l2,prefac(itempy),xinner(itempy))
                  zipftd(itempy,1) = zipftd(itempy,1) + zjtest
                  zipftd(itempy,2) = zipftd(itempy,2) + zjpass
                  zijkld(itempy,1) = zijkld(itempy,1) + zitest
                  zijkld(itempy,2) = zijkld(itempy,2) + zipass
c
 120           continue
 130        continue
 140     continue
 150  continue
_IF(parallel)
c***   ***node-MPP***
       call pg_dlbpush
c...   add fock-matrix and send it around
         nbsq = nb*nb
         call pg_dgop(600,fock,nbsq,'+')
c***   ***node-MPP***
_ENDIF
c
      call gmem_free_inf(k0a,fnm,snm,'k0a')
c
      return
      end
      subroutine dengrd(d,pold)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen
      common/defunk/cobas(maxorb,3),nd
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      dimension d(n,n),pold(n,n)
c
c**** to start with density matrix is in pold
c**** density finishes in pold weighted by root 3
c**** compressed density ends in d
c
      call vclr(d,1,n*n)
      dcmax = 0.0d0
      nsnp = ns + np
      do 50 i = 1 , nsnp + nd
         if (i.le.ns) then
            ivl1 = i
            ivl2 = i
         end if
         if (i.gt.ns .and. i.le.nsnp) then
            ivl1 = ns + (i-ns-1)*3 + 1
            ivl2 = ivl1 + 2
         end if
         if (i.gt.nsnp) then
            ivl1 = ns + np3 + (i-nsnp-1)*6 + 1
            ivl2 = ivl1 + 5
         end if
         do 40 j = 1 , i
            if (j.le.ns) then
               jvl1 = j
               jvl2 = j
            end if
            if (j.gt.ns .and. j.le.nsnp) then
               jvl1 = ns + (j-ns-1)*3 + 1
               jvl2 = jvl1 + 2
            end if
            if (j.gt.nsnp) then
               jvl1 = ns + np3 + (j-nsnp-1)*6 + 1
               jvl2 = jvl1 + 5
            end if
            do 30 ib = ivl1 , ivl2
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
               do 20 jb = jvl1 , jvl2
                  d(i,j) = max(d(i,j),dabs(pold(ib,jb)))
 20            continue
 30         continue
 40      continue
 50   continue
      do 70 i = 1 , nd + nsnp
         do 60 j = 1 , i
            dcmax = max(dcmax,d(i,j))
            d(j,i) = d(i,j)
 60      continue
 70   continue
      dcmax = dcmax*6.0d0
c
      return
      end
      subroutine diisn(f,p,gp,index,nerr,amag)
c
c***** routine performs pulays diis extrapolation
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (nl=10)
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1
      common/worksp/a(nl+1,nl+1),b(nl+1),aa(nl+1,nl+1)
     &,wks1(nl+1),wks2(nl+1),c(nl+1)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n
INCLUDE(common/wrtd)
      common/iofile/iread,iwr
      logical ofull,odiisn
      common/gen3/ofull(100),odiisn
      dimension gp(n,n),f(n,n),p(n,n)
c
c     *****  n is actual number of basis functions ********
c     *****  nl is the maximum number of error vectors ****
c     *****  stored.index is the position in errvec at ****
c     *****  which we will store the vector ***************
c
      call wrt3(f,n*n,idblk(5+index),numdir)
c
      if (.not.odiisn) return
      if (na.lt.3 .and. n.le.20) return
      call rdedx(p,n*n,idblk(4),numdir)
c
      call vclr(gp,1,n*n)
      call mxmb(p,1,n,f,1,n,gp,1,n,n,n,n)
c
      call rdedx(p,n*n,idblk(2),numdir)
c
      call vclr(f,1,n*n)
      call mxmb(p,1,n,gp,1,n,f,1,n,n,n,n)
c
      do 30 i = 1 , n
         do 20 j = 1 , i
            f(i,j) = f(i,j) - f(j,i)
 20      continue
 30   continue
c
      call wrt3(f,n*n,idblk(15+index),numdir)
c
      amag = 0.0d0
      do 50 k = 1 , n
         do 40 l = 1 , k
            amag = amag + f(k,l)*f(k,l)
 40      continue
 50   continue
_IF1()c      write(iwr,747) amag
_IF1()c747   format(1x,'error vector norm ',d12.2)
      if (nerr.lt.2) then
         a(1,1) = 0.0d0
         do 70 k = 1 , n
            do 60 l = 1 , k
               a(1,1) = a(1,1) + f(k,l)*f(k,l)
 60         continue
 70      continue
         call rdedx(f,n*n,idblk(5+index),numdir)
         return
      end if
c
      ndim = nerr + 1
c
c     *****  index is the position of the new parameter & *****
c     *****  error vector in fockm and errvec respectively. ***
c     *****  nerr is how many of the nl available places in ***
c     *****  fockm & errvec contain vectors *******************
c
c     set up diis matrix
c
      do 80 i = 1 , nerr
         a(ndim,i) = -1.0d0
         a(i,ndim) = -1.0d0
         b(i) = 0.0d0
 80   continue
      i = index
      do 110 j = 1 , nerr
         a(i,j) = 0.0d0
         call rdedx(p,n*n,idblk(15+j),numdir)
         do 100 k = 1 , n
            do 90 l = 1 , k
               a(i,j) = a(i,j) + f(k,l)*p(k,l)
 90         continue
 100     continue
         a(j,i) = a(i,j)
 110  continue
      a(ndim,ndim) = 0.0d0
      b(ndim) = -1.0d0
c
c     solve simultaneous equations
c     and use solutions to form new fock matrix
c
      ifail = 0
      call f04atf(a,nl+1,b,ndim,c,aa,nl+1,wks1,wks2,ifail)
      call vclr(f,1,n*n)
c
      do 140 k = 1 , nerr
         call rdedx(p,n*n,idblk(5+k),numdir)
         do 130 i = 1 , n
            do 120 j = 1 , i
               f(i,j) = f(i,j) + c(k)*p(i,j)
 120        continue
 130     continue
 140  continue
c
      do 160 i = 1 , n
         do 150 j = 1 , i - 1
            f(j,i) = f(i,j)
 150     continue
 160  continue
c
      return
      end
      subroutine dnsprm(d,pold,gp,outon)
c
c***** routine calculates the compressed density matrix.
c
c      d is the true density matrix
c
      implicit REAL  (a-h,o-z)
      logical ofull,odiisn,outon
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen
      common/defunk/cobas(maxorb,3),nd
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      common/gen3/ofull(100),odiisn
      common/iofile/iread,iwr
      dimension d(n,n),pold(n,n),gp(n,n)
c
      if (iter.eq.80 .or. iter.eq.90.or.ofull(iter) ) then
c
c***** this corresponds to a complete recalculation of fock matrix
c
         if (outon) write (iwr,200) iter
         do 30 i = 1 , n
            do 20 j = 1 , n
               gp(i,j) = d(i,j)
 20         continue
 30      continue
        dcmax = dabs(gp(1,1)-pold(1,1))
        do 180 i = 1 , n
           do 190 j = 1 , i
              x13 = dabs(gp(i,j)-pold(i,j))
              if (x13.gt.dcmax) dcmax = x13
 190        continue
 180     continue
      else
c
c***** this corresponds to calculating the fock matrix iteratively
c
         do 50 i = 1 , n
            do 40 j = 1 , n
               gp(i,j) = d(i,j) - pold(i,j)
 40         continue
 50      continue
c
        dcmax = dabs(gp(1,1))
        do 70 i = 1 , n
           do 60 j = 1 , i
            x13 = dabs(gp(i,j))
            if (x13.gt.dcmax) dcmax = x13
 60        continue
 70     continue
      end if
c
      do 90 i = 1 , n
         do 80 j = 1 , n
            pold(i,j) = d(i,j)
 80      continue
 90   continue
c
c**** all this code is doing is working out a compressed density
c**** matrix with the absolute maximum element for each shell
c**** combination: result is in d.
c
      call vclr(d,1,n*n)
      nsnp = ns + np
      do 130 i = 1 , nsnp + nd
         if (i.le.ns) then
            ivl1 = i
            ivl2 = i
         end if
         if (i.gt.ns .and. i.le.nsnp) then
            ivl1 = ns + (i-ns-1)*3 + 1
            ivl2 = ivl1 + 2
         end if
         if (i.gt.nsnp) then
            ivl1 = ns + np3 + (i-nsnp-1)*6 + 1
            ivl2 = ivl1 + 5
         end if
         do 120 j = 1 , i
            if (j.le.ns) then
               jvl1 = j
               jvl2 = j
            end if
            if (j.gt.ns .and. j.le.nsnp) then
               jvl1 = ns + (j-ns-1)*3 + 1
               jvl2 = jvl1 + 2
            end if
            if (j.gt.nsnp) then
               jvl1 = ns + np3 + (j-nsnp-1)*6 + 1
               jvl2 = jvl1 + 5
            end if
            do 110 ib = ivl1 , ivl2
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
               do 100 jb = jvl1 , jvl2
                  d(i,j) = max(d(i,j),dabs(gp(ib,jb)))
 100           continue
 110        continue
 120     continue
 130  continue
      do 150 i = 1 , nd + nsnp
         do 140 j = 1 , i
            d(j,i) = d(i,j)
 140     continue
 150  continue
c
c**** code is weighting the dxz dxy dyz element of the density matrix
c**** by root 3
c
      sq3 = dsqrt(3.0d0)
      if (nd.gt.0) then
         ic = 0
         do 170 i = 1 , n
            weigh1 = 1.0d0
            ix = i - ns - np3
            if (ix.gt.0) then
               ic = ic + 1
               if (ic.eq.2 .or. ic.eq.3 .or. ic.eq.5) weigh1 = sq3
               if (ic.eq.6) ic = 0
            end if
            jc = 0
            do 160 j = 1 , n
               weigh2 = 1.0d0
               jx = j - ns - np3
               if (jx.gt.0) then
                  jc = jc + 1
                  if (jc.eq.2 .or. jc.eq.3 .or. jc.eq.5) weigh2 = sq3
                  if (jc.eq.6) jc = 0
               end if
               gp(i,j) = gp(i,j)*weigh1*weigh2
 160        continue
 170     continue
      end if
 200  format(1x,' **** recalculation of fock matrix on iteration ',i3)
      return
      end
      subroutine dp1(xijk,xija,xijp1,xijp2,xijp3,i1,i2
     &,j1,j2,xt1,xt2,wx
     &,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ccc1,ccc2)
c
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (d/p) and (d/h/p) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      dimension xija(n11),xijk(n11),ccc2(n11,3)
     &,s1(n11*3,6),h1(n11*3,6),sx(nb,nb),hx(nb,nb),ccc1(n11,3)
     &,lar(n11),f(8*n11),t(n11),xt1(n14*6,n13*6),xt2(n16*3,n15*3)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
     &,xijp1(n11),xijp2(n11),xijp3(n11)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11*18)
      icount = 0
c
      do 30 i = imin , imax
         do 20 j = jmin , jmax
            icount = icount + 1
            ccc1(icount,1) = xijp1(icount) - coprm(i,1)
            ccc1(icount,2) = xijp2(icount) - coprm(i,2)
            ccc1(icount,3) = xijp3(icount) - coprm(i,3)
            ccc2(icount,1) = xijp1(icount) - coprm(j,1)
            ccc2(icount,2) = xijp2(icount) - coprm(j,2)
            ccc2(icount,3) = xijp3(icount) - coprm(j,3)
 20      continue
 30   continue
      do 60 k = 1 , na
c
         do 40 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 40      continue
c
         call fquik2(f,f,lar,t,4,n11)
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 50 icount = 1 , n11
            index1 = n11 + icount
            index2 = n11 + index1
            index3 = n11 + index2
c           index4 = n11 + index3
c           index5 = n11 + index4
            x7 = xijk(icount)*charge(k)/pi14
            f(icount) = f(icount)*x7
            f(index1) = f(index1)*x7
            f(index2) = f(index2)*x7
            f(index3) = f(index3)*x7
            xij2 = 1.0d0/(xija(icount)+xija(icount))
            x40 = (f(icount)-f(index1))*xij2
            x41 = (f(index1)-f(index2))*xij2
            x1a = xijp1(icount) - coord(1,k)
            x1b = xijp2(icount) - coord(2,k)
            x1c = xijp3(icount) - coord(3,k)
c
c     ***** xx1 has the k dependant factors of (pi/v/s)m=0 ***** 
c     ***** xx2 has the k dependant factors of (pi/v/s)m=1 ***** 
c     ***** xx3 has the k dependant factors of (pi/v/s)m=2 ***** 
c     ***** xx4 contains (1/2@)*( (pi/v/s)m=0 - (pi/v/s)m=1 )
c
            x1x = ccc1(icount,1)*f(icount) - x1a*f(index1)
            x2x = ccc1(icount,1)*f(index1) - x1a*f(index2)
            x3x = ccc1(icount,1)*f(index2) - x1a*f(index3)
            x4x = xij2*(x1x-x2x)
            x1y = ccc1(icount,2)*f(icount) - x1b*f(index1)
            x2y = ccc1(icount,2)*f(index1) - x1b*f(index2)
            x3y = ccc1(icount,2)*f(index2) - x1b*f(index3)
            x4y = xij2*(x1y-x2y)
            x1z = ccc1(icount,3)*f(icount) - x1c*f(index1)
            x2z = ccc1(icount,3)*f(index1) - x1c*f(index2)
            x3z = ccc1(icount,3)*f(index2) - x1c*f(index3)
            x4z = xij2*(x1z-x2z)
c
c***** xd1 and xd2 contain (d|v|s)m=0 and (d|v|s)m=1 respectively
c**** dxx
            xd1a = ccc1(icount,1)*x1x - x1a*x2x + x40
            xd2a = ccc1(icount,1)*x2x - x1a*x3x + x41
c**** dxy
            xd1b = ccc1(icount,2)*x1x - x1b*x2x
            xd2b = ccc1(icount,2)*x2x - x1b*x3x
c**** dxz
            xd1c = ccc1(icount,3)*x1x - x1c*x2x
            xd2c = ccc1(icount,3)*x2x - x1c*x3x
c**** dyy
            xd1d = ccc1(icount,2)*x1y - x1b*x2y + x40
            xd2d = ccc1(icount,2)*x2y - x1b*x3y + x41
c**** dyz
            xd1e = ccc1(icount,3)*x1y - x1c*x2y
            xd2e = ccc1(icount,3)*x2y - x1c*x3y
c**** dzz
            xd1f = ccc1(icount,3)*x1z - x1c*x2z + x40
            xd2f = ccc1(icount,3)*x2z - x1c*x3z + x41
c
c***** (dij|v|pk) goes in h1(i)
c
c***** (dxx|v|px-z)
            h1(icount,1) = h1(icount,1) + ccc2(icount,1)
     +                     *xd1a - x1a*xd2a + x4x + x4x
            h1(index1,1) = h1(index1,1) + ccc2(icount,2)*xd1a - x1b*xd2a
            h1(index2,1) = h1(index2,1) + ccc2(icount,3)*xd1a - x1c*xd2a
c***** (dxy|v|px-z)
            h1(icount,2) = h1(icount,2) + ccc2(icount,1)
     +                     *xd1b - x1a*xd2b + x4y
            h1(index1,2) = h1(index1,2) + ccc2(icount,2)
     +                     *xd1b - x1b*xd2b + x4x
            h1(index2,2) = h1(index2,2) + ccc2(icount,3)*xd1b - x1c*xd2b
c***** (dxz|v|px-z)
            h1(icount,3) = h1(icount,3) + ccc2(icount,1)
     +                     *xd1c - x1a*xd2c + x4z
            h1(index1,3) = h1(index1,3) + ccc2(icount,2)*xd1c - x1b*xd2c
            h1(index2,3) = h1(index2,3) + ccc2(icount,3)
     +                     *xd1c - x1c*xd2c + x4x
c***** (dyy|v|px-z)
            h1(icount,4) = h1(icount,4) + ccc2(icount,1)*xd1d - x1a*xd2d
            h1(index1,4) = h1(index1,4) + ccc2(icount,2)
     +                     *xd1d - x1b*xd2d + x4y + x4y
            h1(index2,4) = h1(index2,4) + ccc2(icount,3)*xd1d - x1c*xd2d
c***** (dyz|v|px-z)
            h1(icount,5) = h1(icount,5) + ccc2(icount,1)*xd1e - x1a*xd2e
            h1(index1,5) = h1(index1,5) + ccc2(icount,2)
     +                     *xd1e - x1b*xd2e + x4z
            h1(index2,5) = h1(index2,5) + ccc2(icount,3)
     +                     *xd1e - x1c*xd2e + x4y
c***** (dzz|v|px-z)
            h1(icount,6) = h1(icount,6) + ccc2(icount,1)*xd1f - x1a*xd2f
            h1(index1,6) = h1(index1,6) + ccc2(icount,2)*xd1f - x1b*xd2f
            h1(index2,6) = h1(index2,6) + ccc2(icount,3)
     +                     *xd1f - x1c*xd2f + x4z + x4z
c
 50      continue
 60   continue
c
      icount = 0
      do 80 i = imin , imax
         do 70 j = jmin , jmax
            icount = icount + 1
            index1 = n11 + icount
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(index1) = pe(i)
            f(icount) = pe(i)*pe(j)
 70      continue
 80   continue
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 90 icount = 1 , n11
         index1 = n11 + icount
         index2 = n11 + index1
c        index3 = n11 + index2
c        index4 = n11 + index3
c        index5 = n11 + index4
c
c     ***** x4 is (s||s), x3 is (s|t|s)
c
         x1 = 1.0d0/xija(icount)
         x2 = f(icount)*x1
         x4 = xijk(icount)*pi14*dsqrt(x1)
         x1 = 0.5d0*x1
         x22 = 2.0d0*x2
         x3 = x2*(3.0d0-x22*t(icount))*x4
         x8 = x2*x4/f(index1)
c
c***** x1x-z is (pi||s) and xtx-z is (pi|v|s)
         x1x = ccc1(icount,1)*x4
         x1y = ccc1(icount,2)*x4
         x1z = ccc1(icount,3)*x4
c
         xtx = ccc1(icount,1)*x3 + x22*x1x
         xty = ccc1(icount,2)*x3 + x22*x1y
         xtz = ccc1(icount,3)*x3 + x22*x1z
c***** xdxx-xdzz is (dij||s) and xtxx-xtzz is (dij||s)
         x6 = x1*x4
         xdxx = ccc1(icount,1)*x1x + x6
         xdxy = ccc1(icount,1)*x1y
         xdxz = ccc1(icount,1)*x1z
         xdyy = ccc1(icount,2)*x1y + x6
         xdyz = ccc1(icount,2)*x1z
         xdzz = ccc1(icount,3)*x1z + x6
c
         xt6 = x1*x3 - x8
         xtxx = ccc1(icount,1)*xtx + xdxx*x22 + xt6
         xtxy = ccc1(icount,1)*xty + xdxy*x22
         xtxz = ccc1(icount,1)*xtz + xdxz*x22
         xtyy = ccc1(icount,2)*xty + xdyy*x22 + xt6
         xtyz = ccc1(icount,2)*xtz + xdyz*x22
         xtzz = ccc1(icount,3)*xtz + xdzz*x22 + xt6
c
c***** now do the overlap integrals in s1 and the hcore integrals
c***** in h1
         x1x = x1x*x1
         x1y = x1y*x1
         x1z = x1z*x1
         xtx = xtx*x1
         xty = xty*x1
         xtz = xtz*x1
c***** (dxx||px-z) and (dxx|h|px-z)
         s1(icount,1) = ccc2(icount,1)*xdxx + x1x + x1x
         s1(index1,1) = ccc2(icount,2)*xdxx
         s1(index2,1) = ccc2(icount,3)*xdxx
         h1(icount,1) = ccc2(icount,1)*xtxx + x22*s1(icount,1)
     +                  - h1(icount,1) + xtx + xtx
         h1(index1,1) = ccc2(icount,2)*xtxx + x22*s1(index1,1)
     +                  - h1(index1,1)
         h1(index2,1) = ccc2(icount,3)*xtxx + x22*s1(index2,1)
     +                  - h1(index2,1)
c***** (dxy||px-z) and (dxy|h|px-z)
         s1(icount,2) = ccc2(icount,1)*xdxy + x1y
         s1(index1,2) = ccc2(icount,2)*xdxy + x1x
         s1(index2,2) = ccc2(icount,3)*xdxy
         h1(icount,2) = ccc2(icount,1)*xtxy + x22*s1(icount,2)
     +                  - h1(icount,2) + xty
         h1(index1,2) = ccc2(icount,2)*xtxy + x22*s1(index1,2)
     +                  - h1(index1,2) + xtx
         h1(index2,2) = ccc2(icount,3)*xtxy + x22*s1(index2,2)
     +                  - h1(index2,2)
c***** (dxz||px-z) and (dxz|h|px-z)
         s1(icount,3) = ccc2(icount,1)*xdxz + x1z
         s1(index1,3) = ccc2(icount,2)*xdxz
         s1(index2,3) = ccc2(icount,3)*xdxz + x1x
         h1(icount,3) = ccc2(icount,1)*xtxz + x22*s1(icount,3)
     +                  - h1(icount,3) + xtz
         h1(index1,3) = ccc2(icount,2)*xtxz + x22*s1(index1,3)
     +                  - h1(index1,3)
         h1(index2,3) = ccc2(icount,3)*xtxz + x22*s1(index2,3)
     +                  - h1(index2,3) + xtx
c***** (dyy||px-z) and (dyy|h|px-z)
         s1(icount,4) = ccc2(icount,1)*xdyy
         s1(index1,4) = ccc2(icount,2)*xdyy + x1y + x1y
         s1(index2,4) = ccc2(icount,3)*xdyy
         h1(icount,4) = ccc2(icount,1)*xtyy + x22*s1(icount,4)
     +                  - h1(icount,4)
         h1(index1,4) = ccc2(icount,2)*xtyy + x22*s1(index1,4)
     +                  - h1(index1,4) + xty + xty
         h1(index2,4) = ccc2(icount,3)*xtyy + x22*s1(index2,4)
     +                  - h1(index2,4)
c***** (dyz||px-z) and (dyz|h|px-z)
         s1(icount,5) = ccc2(icount,1)*xdyz
         s1(index1,5) = ccc2(icount,2)*xdyz + x1z
         s1(index2,5) = ccc2(icount,3)*xdyz + x1y
         h1(icount,5) = ccc2(icount,1)*xtyz + x22*s1(icount,5)
     +                  - h1(icount,5)
         h1(index1,5) = ccc2(icount,2)*xtyz + x22*s1(index1,5)
     +                  - h1(index1,5) + xtz
         h1(index2,5) = ccc2(icount,3)*xtyz + x22*s1(index2,5)
     +                  - h1(index2,5) + xty
c***** (dzz||px-z) and (dzz|h|px-z)
         s1(icount,6) = ccc2(icount,1)*xdzz
         s1(index1,6) = ccc2(icount,2)*xdzz
         s1(index2,6) = ccc2(icount,3)*xdzz + x1z + x1z
         h1(icount,6) = ccc2(icount,1)*xtzz + x22*s1(icount,6)
     +                  - h1(icount,6)
         h1(index1,6) = ccc2(icount,2)*xtzz + x22*s1(index1,6)
     +                  - h1(index1,6)
         h1(index2,6) = ccc2(icount,3)*xtzz + x22*s1(index2,6)
     +                  - h1(index2,6) + xtz + xtz
c
 90   continue
c
c*****  put the result in the matrices to be transformed ***
c
      do 130 ic = 1 , 6
         icount = 0
         do 120 jc = 1 , 3
            do 110 i = 1 , n14
               ip = (i-1)*6
               do 100 j = 1 , n16
                  jp = (j-1)*3
                  icount = icount + 1
                  wx(ip+ic,jp+jc) = s1(icount,ic)
                  wy(ip+ic,jp+jc) = h1(icount,ic)
 100           continue
 110        continue
 120     continue
 130  continue
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16*3,wz,1,ipar2,n14*6,n16*3,n15*3)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wx,ipar2,1,n15*3,n14*6,n13*6)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16*3,wz,1,ipar2,n14*6,n16*3,n15*3)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wy,ipar2,1,n15*3,n14*6,n13*6)
c
      ix1 = 6*i1 - 3*np - ns*5 - 5
      ix2 = ix1 + n13*6 - 1
      jx1 = 3*j1 - ns*2 - 2
      jx2 = jx1 + n15*3 - 1
      ic = 0
      do 150 i = ix1 , ix2
         ic = ic + 1
         jc = 0
         do 140 j = jx1 , jx2
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 140     continue
 150  continue
c
      return
      end
      subroutine drccon(dens,sx,iwksp,iwksp2,labsh,isoc,
     +                  nwauq,istart,nnp)
c
      implicit REAL  (a-h,o-z)
      character *4 char
      logical odiis
INCLUDE(common/sizes)
c
c
c direct initialisation
c
*
* istart is the start of memory for the memory allocation section
*
      dimension dens(nnp,nnp),sx(*)
      dimension iwksp(mxshel,48),iwksp2(mxshel,48),labsh(*),isoc(*)
      dimension nwauq(maxat)
      dimension itrans(6)
*
* these are the gamess common blocks
*
      common/scfopt/maxit,mconv,nconv,npunch,accdi1,accdi2,odiis
     *           ,icoupl(3),dmpcut,
     *            acurcy,en,etot,ehf,ehf0,diff,iterg,icount
     *            ,rshift,exttol,dmptol,vshtol,iextin
     *            ,iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp
      common/transf/psmal,qsmal,rsmal,pnew,qnew,rnew,pp,qp,rp
INCLUDE(common/scra7)
      common/tol/toler,tol2
INCLUDE(common/restar)
INCLUDE(common/prints)
      common/phycon/toang,phyc(83),iunt,nosym
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/infoa)
      logical ozmat
      common/infob/czanr(maxat),czin(maxat),cin(3,maxat),amass(maxat),
     -       c80(maxat3,3),nonsym,map80(maxat),ozmat
INCLUDE(common/nshel)
INCLUDE(common/machin)
INCLUDE(common/runlab)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480)
c
c
c**** the common blocks "data" to "gen" store all the main informayion
c**** used by the direct program.
c
c**** the meanings of the non - obvious arrays in common are:
c  nwa is a poiunter array that indicates to which atom a particular
c  primitive is associated;  isatbf is the number of s basis functions
c  associated with particular atoms ans isatp is the number of primtives
c  shells on the atoms; ipatbf and ipatp are similarly defined for p
c  functions whilst idatbf and idatp are those for d functions.
c  nwshp points to which shell particular primitives belong; nup gives t
c  number (in the primitive list) of the final primitive in a particular
c  contracted basis function (and is heavily used); pe is the primitive
c  exponents and pc the contraction coefficients; whilst the well-used
c  array coprm gives the coordinate position of particular primitives
c  which compares with cobas in common/defunk which gives the coordinate
c  postion of basis funnctions.
c
c  note that generally different arrays are needed for the storage of da
c  to that in gamess because the basis function list that is read in, in
c  gamess format is sorted so that all the s functions occur first, then
c  the p functions, and then the d functions. some of the info on d
c  functions is in the common block defunk.
c  the important non-array data in common blocks data or defunk
c   na.... number of atoms     nb..... number of basis functions
c   ns.... number of s function  nsp number of s primitives
c   np.... number of p function  npp number of p primitves
c   nd and ndp for d funcs
c   (note that where names clash with gamess definitions in this routine
c   the letters dir (for direct) have been added to the names)
*
      common/data/charg(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     - isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim),
     - nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),ndira,ndirb,
     - ndirs,ndirp,nsp,n1,ndirpp,n2,
     - n3,np3,npp3,nprm,n4,nm1,nm2,ndire,repen
c
c*****  gen2 contains the integral test data
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
c*****  gen contains most of the stuff about the options
      logical glog1,glog2,glog5,glog6,glog9,gdbg,gmp2
      logical readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf
      logical gmull
      common/gen/maxitd,iter,concrit,conv,energy,thresh,dcmax,
     -      glog1,glog2,glog5,glog9,glog6,gdbg(10),gmp2,mp2fo,mp2fv,
     -      readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold,
     -      gscf,noptd,npoint,timeg,gmull,iswap,dshift,irhess
INCLUDE(common/wrtd)
      common/gtrans/idtra(maxorb)
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp,ipard
     &,nsspl,npspl,ndspl
c
*!!!!! memory allocator
      common/igarg/imgarg(40),itemp,ileft,iwhole
INCLUDE(common/symtry)
      common/scrtch/ ptremp(3,144),dtremp(6,288)
      common/scra /iso(mxshel,48),ntd,nauq
c
      dimension char(21)
      data itrans/1,4,6,2,3,5/
      data char/'ssss','psss','psps','ppss','ppps','pppp',
     -       'dsss','dsps','dspp','dsds','dpss','dpps','dppp',
     -       'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
c
      thresh = accdi1
c
c*****  fgrid forms the grid form which the incomplete
c*****  gamma function is evaluated via a taylor series
c
      call fgrid
*
* convert data into direct format
*
c
c******** now store the gamess data in a way more **********
c******** suitable for direct **********
c
      do 20 i = 1 , nat
         nwauq(i) = 0
         charg(i) = czan(i)
         coord(1,i) = c(1,i)
         coord(2,i) = c(2,i)
         coord(3,i) = c(3,i)
 20   continue
      do 30 i = 1 , nat
         isatbf(i) = 0
         isatp(i) = 0
         ipatbf(i) = 0
         ipatp(i) = 0
         idatbf(i) = 0
         idatp(i) = 0
 30   continue
*
* read in symmetry
*
      if (gdbg(1)) write (iwr,*) 'setting symmetry labels for shells'
      ntd = nt
      nav = lenwrd()
      call readi(labsh,nw196(5)*nav,ibl196(5),idaf)
*!!!      write(iwr,*)'nshell',nshell
      do 120 i = 1 , nshell
*!!!      write(iwr,*)'shell ',i
         do 110 j = 1 , nt
            iso(i,j) = labsh(i+iliso(j))
 110     continue
*!!!         write(iwr,*)(iso(i,j),j=1,nt)
 120  continue
*
* convert symmetry arrays to direct format
*
      call readi(isoc,nw196(6)*nav,ibl196(6),idaf)
      do 60 i = 1 , nonsym
         do 50 it = 2 , nt
            nnew = isoc(i+ilisoc(it))
*!!!          write(iwr,*)'looking at symmetry operation ',it,' atom ',
*!!!     + nnew
*
* if we find a symmetry equivalent atom set nwauq(i) to 1
*
            if (nwauq(i).eq.0 .and. nnew.ne.i) nwauq(nnew) = 1
 50      continue
c
 60   continue
c
c******* rearrange the data so first s, then p, and then d are considere
c******* loop over the s functions only *****
c
      icount = 0
      jcount = 0
      do 150 i = 1 , nshell
*!!!      write(iwr,*) 'copying symmetry for shell i'
         if (ktype(i).gt.3) call caserr(
     +   'the program can not handle f or g functions')
         if (kmin(i).eq.1) then
            jcount = jcount + 1
            do 130 k = 1 , nt
               iwksp(jcount,k) = iso(i,k)
 130        continue
            labsh(i) = jcount
            kk = katom(i)
            nc(jcount) = kng(i)
            isatbf(kk) = isatbf(kk) + 1
*!!!        write(iwr,*)'copying s shell info'
            do 140 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
               isatp(kk) = isatp(kk) + 1
*!!!         write(iwr,*)icount,ex(k),cs(k)
               pe(icount) = ex(k)
               pc(icount) = cs(k)
 140        continue
         end if
 150  continue
      nsp = icount
      ndirs = jcount
*!!!      write(iwr,*)'storing shell data for s',ndirs
      do 170 it = 1 , nt
         do 160 i = 1 , ndirs
            iwksp2(i,it) = labsh(iwksp(i,it))
 160     continue
 170  continue
c
c******* loop over the p functions only *****
c
      do 200 i = 1 , nshell
         if (ktype(i).eq.2) then
            jcount = jcount + 1
            nc(jcount) = kng(i)
            kk = katom(i)
            ipatbf(kk) = ipatbf(kk) + 1
            do 180 k = 1 , nt
               iwksp(jcount,k) = iso(i,k)
 180        continue
            labsh(i) = jcount
*!!!        write(iwr,*)'copying p shell info'
            do 190 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
*!!!          write(iwr,*)icount,ex(k),cp(k)
               ipatp(kk) = ipatp(kk) + 1
               pe(icount) = ex(k)
               pc(icount) = cp(k)
 190        continue
         end if
 200  continue
      do 220 it = 1 , nt
         do 210 i = ndirs + 1 , jcount
            iwksp2(i,it) = labsh(iwksp(i,it))
 210     continue
 220  continue
      ndirpp = icount - nsp
      ndirp = jcount - ndirs
c
c******* loop over the d functions only *****
c
      do 250 i = 1 , nshell
         if (ktype(i).eq.3) then
            jcount = jcount + 1
            nc(jcount) = kng(i)
            kk = katom(i)
            idatbf(kk) = idatbf(kk) + 1
            do 230 k = 1 , nt
               iwksp(jcount,k) = iso(i,k)
 230        continue
            labsh(i) = jcount
*!!!        write(iwr,*)'copying d shell info'
            do 240 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
               idatp(kk) = idatp(kk) + 1
*!!!          write(iwr,*)icount,ex(k),cd(k)
               pe(icount) = ex(k)
               pc(icount) = cd(k)
 240        continue
         end if
 250  continue
      ndirdp = icount - nsp - ndirpp
      ndp = ndirdp
      ndird = jcount - ndirs - ndirp
      nd = ndird
      do 270 it = 1 , nt
         do 260 i = ndirs + ndirp + 1 , jcount
            iwksp2(i,it) = labsh(iwksp(i,it))
 260     continue
 270  continue
      do 290 it = 1 , nt
         do 280 i = 1 , jcount
            iso(i,it) = iwksp2(i,it)
 280     continue
 290  continue
c
c**** now form transformation matrix that converts gamess
c**** storage of basis fuctions into direct's (or vice versa)
c
      jcount = 0
      nsc = 0
      np3 = 3*ndirp
      npc = ndirs
      ndc = np3 + npc
      do 300 i = 1 , nshell
         if (kmin(i).eq.1) then
            nsc = nsc + 1
            jcount = jcount + 1
            idtra(jcount) = nsc
         end if
         if (ktype(i).eq.2) then
            npc = npc + 1
            jcount = jcount + 1
            idtra(jcount) = npc
            npc = npc + 1
            jcount = jcount + 1
            idtra(jcount) = npc
            npc = npc + 1
            jcount = jcount + 1
            idtra(jcount) = npc
         end if
         if (ktype(i).eq.3) then
            idtra(jcount+1) = ndc + 1
            idtra(jcount+2) = ndc + 4
            idtra(jcount+3) = ndc + 6
            idtra(jcount+4) = ndc + 2
            idtra(jcount+5) = ndc + 3
            idtra(jcount+6) = ndc + 5
            ndc = ndc + 6
            jcount = jcount + 6
         end if
 300  continue
c
      if (gdbg(3)) then
         write (iwr,*) '**************************'
         write (iwr,*) 'switching off symmetry now'
         write (iwr,*) '**************************'
         nt = 1
      end if
c
c*****  rearrange data in ixatbf and ixatp so that only unique atoms
c*****   count
c
      nauq = nat
      if (nt.gt.1) then
         call rdedx(ptr,nw196(1),ibl196(1),idaf)
         call rdedx(dtr,nw196(2),ibl196(2),idaf)
         call rdedx(ftr,nw196(3),ibl196(3),idaf)
c
c .. rearrange dtr array
c
         do 320 it = 1 , nt
            idum = 6*(it-1)
            do 310 i = 1 , 6
               ii = itrans(i)
               dtremp(1,idum+ii) = dtr(1,idum+i)
               dtremp(4,idum+ii) = dtr(2,idum+i)
               dtremp(6,idum+ii) = dtr(3,idum+i)
               dtremp(2,idum+ii) = dtr(4,idum+i)
               dtremp(3,idum+ii) = dtr(5,idum+i)
               dtremp(5,idum+ii) = dtr(6,idum+i)
 310        continue
 320     continue
         call dcopy(1728,dtremp,1,dtr,1)
c
         nauq = 0
         do 330 i = 1 , nat
*!!!        write(iwr,*)'atom nwauq(i)',i,nwauq(i)
            if (nwauq(i).eq.0) then
               nauq = nauq + 1
               isatbf(nauq) = isatbf(i)
               isatp(nauq) = isatp(i)
               ipatbf(nauq) = ipatbf(i)
               ipatp(nauq) = ipatp(i)
               idatbf(nauq) = idatbf(i)
               idatp(nauq) = idatp(i)
            else
               isatbf(nauq) = isatbf(nauq) + isatbf(i)
               isatp(nauq) = isatp(nauq) + isatp(i)
               ipatbf(nauq) = ipatbf(nauq) + ipatbf(i)
               ipatp(nauq) = ipatp(nauq) + ipatp(i)
               idatbf(nauq) = idatbf(nauq) + idatbf(i)
               idatp(nauq) = idatp(nauq) + idatp(i)
            end if
 330     continue
      end if
c
      icount = 0
*!!!      write(iwr,*)'ndirs,ndirp,ndird',ndirs,ndirp,ndird
      do 350 i = 1 , ndirs + ndirp + ndird
         do 340 j = 1 , nc(i)
            icount = icount + 1
            nwshp(icount) = i
 340     continue
 350  continue
      icount = 0
      do 360 i = 1 , ndirs + ndirp + ndird
         icount = icount + nc(i)
         nup(i) = icount
 360  continue
c
      if (mp2fo.eq.0) mp2fo = 1
      if (mp2fv.eq.0) mp2fv = num
      npp3 = ndirpp*3
      nprm = nsp + npp3
      n1 = nsp*(nsp+1)/2
      n2 = nsp*ndirpp
      n3 = ndirpp*(ndirpp+1)/2
      n4 = nprm*(nprm+1)/2
*!!!      write(iwr,*)'ne num nat',ne,num,nat
      ndire = ne
      ndirb = num
      ndira = nat
c
c
      imaxzz = 21
      if (ndird.eq.0) imaxzz = 6
*!!!      write(iwr,*)'imaxzz ',imaxzz
      if (gdbg(1)) write (iwr,6060)
      do 370 i = 1 , imaxzz
         if (gdbg(1)) write (iwr,6010) char(i) , prefac(i) , xinner(i)
 370  continue
      if (gdbg(1)) write (iwr,6020) ndirs , nsp , ndirp*3 , ndirpp*3 ,
     +                             ndird*6 , ndirdp*6 , ndirb , nprm
      if (gdbg(1)) write (iwr,6030)
      if (gdbg(1)) write (iwr,6040) (idtra(i),i=1,num)
      if (gdbg(1)) write (iwr,6050)
      if (gdbg(1)) then
         do 380 i = 1 , ndirs + ndirp + ndird
            write (iwr,*) i , nup(i)
 380     continue
         write (iwr,*) 'ns,nsp,ndirp,ndirpp,nprm,n1,n2,n3,n4 are'
         write (iwr,*) ndirs , nsp , ndirp , ndirpp , nprm , n1 , n2 ,
     +                n3 , n4
         write (iwr,*) 'nd,ndp ' , ndird , ndirdp
         write (iwr,*) 'the contractions and exponents are'
         do 390 j = 1 , ndirpp + nsp + ndirdp
            write (iwr,*) pc(j) , pe(j)
 390     continue
         write (iwr,*) 'the array nwa is'
         write (iwr,*) (nwa(j),j=1,nsp+ndirpp+ndirdp)
         write (iwr,*) 'coordinates of atoms are'
         do 400 k = 1 , ndira
            write (iwr,*) (coord(ix,k),ix=1,3)
 400     continue
         write (iwr,*) 'the number of s functions on each atom'
         write (iwr,*) (isatbf(j),j=1,ndira)
         write (iwr,*) 'the number of s primitives on each atom'
         write (iwr,*) (isatp(j),j=1,ndira)
         write (iwr,*) 'the number of p functions on each atom'
         write (iwr,*) (ipatbf(j),j=1,ndira)
         write (iwr,*) 'the number of p primitives on each atom'
         write (iwr,*) (ipatp(j),j=1,ndira)
         write (iwr,*) 'the number of d functions on each atom'
         write (iwr,*) (idatbf(j),j=1,ndira)
         write (iwr,*) 'the number of d primitives on each atom'
         write (iwr,*) (idatp(j),j=1,ndira)
         write (iwr,*) 'the shell to which each primitive belongs'
         write (iwr,*) (nwshp(j),j=1,nsp+ndirpp+ndirdp)
      end if
      if (gdbg(1) .and. nt.gt.1) then
         write (iwr,*) 'number of elements in group' , nt
         write (iwr,*) 'number of symmetry-unique atoms' , nauq
         write (iwr,*) 'the old iso matrix is'
         do 410 k = 1 , nshell
            write (iwr,*) (iwksp(k,it),it=1,nt)
 410     continue
         write (iwr,*) 'the new iso matrix is'
         do 420 k = 1 , ndirs + ndirp + ndird
            write (iwr,*) (iso(k,it),it=1,nt)
 420     continue
      end if
c
      do 430 i = 1 , nsp + ndirpp + ndirdp
         coprm(i,1) = coord(1,nwa(i))
         coprm(i,2) = coord(2,nwa(i))
         coprm(i,3) = coord(3,nwa(i))
 430  continue
      j = 1
      do 440 i = 1 , ndirs + ndirp + ndird
         cobas(i,1) = coprm(j,1)
         cobas(i,2) = coprm(j,2)
         cobas(i,3) = coprm(j,3)
         j = j + nc(i)
 440  continue
c
      nsq = ndirb*ndirb
      nsq2 = ndirb*(ndirb+1)/2
      imgarg(1) = istart
c***** hessian stored at the beginning of a
      imgarg(2) = imgarg(1) + max(nsq,na*na*9)
      imgarg(3) = imgarg(2) + nsq
      imgarg(4) = imgarg(3) + max(nsq,na*9)
      do 450 inx = 5 , 12
         imgarg(inx) = imgarg(4)
 450  continue
      imense  = igmem_max_memory()
      if (imgarg(12).gt.imense) then
         write (iwr,6005) imense, imgarg(12)
         call caserr('not enough core')
      endif
      ileft = imense - imgarg(12)
c
c***** now sorting out memory for d function code. memroy needed is only
c***** for the parameter arrays. this is arbitarily assigned to be
c***** a maximum for each array of 1/100 of the available memory.
      xavail = dsqrt((dfloat(ileft)/200.0d0))
      ipard = aint(xavail)
      iparp = aint(xavail*dsqrt(2.0d0))
      ipars = aint(xavail*dsqrt(6.0d0))
      call countd(1,idefs1,idefs2,isatbf,isatp,ipars,nsspl,nauq)
      if (gdbg(1)) then
         write (iwr,6070) nsspl
         write (iwr,6110) ipars
         write (iwr,6080)
         write (iwr,6090) (idefs1(j),j=1,nsspl)
         write (iwr,6100)
         write (iwr,6090) (idefs2(j),j=1,nsspl)
         write (iwr,6130)
      end if
      if (ndirp.eq.0) then
         iparp = 0
         npspl = 0
         go to 460
      end if
      call countd(ndirs+1,idefp1,idefp2,ipatbf,ipatp,iparp,npspl,nauq)
      if (gdbg(1)) then
         write (iwr,6140) npspl
         write (iwr,6110) iparp
         write (iwr,6080)
         write (iwr,6090) (idefp1(j),j=1,npspl)
         write (iwr,6100)
         write (iwr,6090) (idefp2(j),j=1,npspl)
         write (iwr,6130)
      end if
 460  if (nd.eq.0) then
         ipard = 0
         ndspl = 0
         go to 470
      end if
      call countd(ndirs+ndirp+1,idefd1,idefd2,idatbf,idatp,ipard,ndspl,
     +            nauq)
      if (gdbg(1)) then
         write (iwr,6150) ndspl
         write (iwr,6110) ipard
         write (iwr,6080)
         write (iwr,6090) (idefd1(j),j=1,ndspl)
         write (iwr,6100)
         write (iwr,6090) (idefd2(j),j=1,ndspl)
         write (iwr,6130)
      end if
c
 470  itemp = max(ipard*6,ipars,iparp*3)
      itemp = itemp*itemp
c***** imgarg(12)=imgarg(11) is lar in main
      imgarg(13) = imgarg(11) + itemp
      imgarg(14) = imgarg(13) + itemp
c***** imgarg(13) and imgarg(14) are wx and wy
      imgarg(15) = imgarg(14) + itemp
      imgarg(16) = imgarg(15) + itemp
c***** imgarg(15) and imgarg(16) are x1 and x2
      imgarg(17) = imgarg(16) + itemp
      call timit(0)
      if (gdbg(1)) then
         write (iwr,6120)
         write (iwr,6130)
         write (iwr,6130)
      end if
c**** ileft is the maximum space for the array d in the 1-electron
c**** part
      ileft = imense - imgarg(17)
      iwhole = imense
*
* transfer 1 -particle arrays into direct......
* the overlap matrix in the a.o. basis ?
*
      if (gdbg(1)) write (iwr,*) 'getting overlap' , ibl7s
      call rdedx(sx,nsq2,ibl7s,num8)
      k = 0
      do 490 i = 1 , ndirb
         id = idtra(i)
         do 480 j = 1 , i
            jd = idtra(j)
            k = k + 1
            dens(id,jd) = sx(k)
            dens(jd,id) = sx(k)
 480     continue
 490  continue
      if (gdbg(1)) write (iwr,*) 'writing overlap '
      call wrt3(dens,nsq,idblk(2),numdir)
c     call prsq(dens,ndirb,ndirb,ndirb)
*
* the core hamiltonian in the a.o.basis
*
      if (gdbg(1)) write (iwr,*) 'getting core matrix' , ibl7f
      call rdedx(sx,nsq2,ibl7f,num8)
      k = 0
      do 510 i = 1 , ndirb
         id = idtra(i)
         do 500 j = 1 , i
            jd = idtra(j)
            k = k + 1
            dens(id,jd) = sx(k)
            dens(jd,id) = sx(k)
 500     continue
 510  continue
      if (gdbg(1)) write (iwr,*) 'writing core matrix'
      call wrt3(dens,nsq,idblk(1),numdir)
c     call prsq(dens,ndirb,ndirb,ndirb)
      if (.not.oentry) then
*
* get the density matrix
*
         if (gdbg(1)) write (iwr,*) 'getting density matrix' , ibl3pa ,
     +                             idaf
         call rdedx(sx,nsq2,ibl3pa,idaf)
         k = 0
         do 530 i = 1 , ndirb
            id = idtra(i)
            do 520 j = 1 , i
               jd = idtra(j)
               k = k + 1
               dens(id,jd) = sx(k)
               dens(jd,id) = sx(k)
 520        continue
 530     continue
         if (gdbg(1)) write (iwr,*) 'writing density matrix'
         call wrt3(dens,nsq,idblk(26),numdir)
         call ctrdir
c      call prsq(dens,ndirb,ndirb,ndirb)
      end if
      oentry = .true.
      return
 6005 format (/' insufficient core for drccon ',i9,' real words ',
     +        ' need ',i9,' real words ')
 6010 format (1x,a8,1x,d10.3,1x,d10.3)
 6020 format (/17x,'basis functions',10x,
     +        'primitives'/'s functions           ',i10,10x,
     +        i10/'p functions           ',i10,10x,
     +        i10/'d functions           ',i10,10x,
     +        i10/'total                 ',i10,10x,i10)
 6030 format (/'direct --> gamess mapping vector')
 6040 format (1x,12i5)
 6050 format (//)
c
 6060 format ('         information on integral factors')
 6070 format (/'number of batches of s shells                    ',i5)
 6080 format ('starting point(s) of each batch')
 6090 format (1x,12i5)
 6100 format ('ending point(s) of each batch')
 6110 format ('maximum number of primitive shellsin each batch       ',
     +        i5)
 6120 format ('finished reading and dealing with the data')
 6130 format (//)
 6140 format (1x,'number of batches of p shells       ',i5)
 6150 format (1x,'number of batches of d shells       ',i5)
      end
      subroutine ctrdir
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/trand/ilifd(maxorb),ntrand(maxorb),itrand(mxorb3),
     + ctrand(mxorb3),otrand
INCLUDE(common/tran)
      common/gtrans/idtra(maxorb)
      common/data/charg(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     + isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     + nwshp(mxprim),
     + nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),ndira,ndirb,
     + ndirs,ndirp,nsp,n1,ndirpp,n2,
     + n3,np3,npp3,nprm,n4,nm1,nm2,ndire,repen
c
      otrand = otran
      if(.not.otrand) then
       do 10 i=1,ndirb
       iii = ilifc(i)
       nt1 = ntran(i)
       ntrand(i) = nt1
       ilifd(i) = iii
       do 20 j=1,nt1
       l = j + iii
       itrand(l) = idtra(itran(l))
       ctrand(l) = ctran(l)
 20    continue
 10    continue
      endif
      return
      end
      block data drcdat
c
c**** the common block index contains critical information used in the
c**** evaluation of the integrals. all the arrays defined here are
c**** used constantly without explanation in the integral code !
c
c**** imagine numbering cartesian basis function components in a shell
c**** in a fixed scheme so that all x parts then all y parts then all
c**** z parts are considered in turn eg.
c******** px   1 : py   2 : pz   3:
c******** dxx  1 : dxy  2 : : : dyz  5 : dzz 6:
c******** fxxx 1 : fxxy 2 : fxxz 3 : fxyy 4 : ......etc
c**************************************************************
c**** this is now the ordering for the components.
c**** the final index (ie. x,y or z equivalent to 1,2 or 3) in the
c**** components of any shell are given by ispdf1(?,1). for instance,
c**** the 6th component in an gshell is xxzz, in an fshell is xzz, in a
c**** dshell is zz; and so ispdf1(6,1) is 3 (z component).
c**** when the final index is removed (ie. the reversal of the applicati
c**** of the vertical recursion relation) a new shell and a new
c**** position results. eg. removing the z component from the 6thg
c**** component in a gshell, xxzz, leaves the fshell component xxz which
c**** from ispdf1(6,2) has the 3rd position in the fshell list.
c
c**** nvl gives the number of cartesian components to particular shells
c**** eg. nvl(5) is a g function with 15 components.
c
c**** ixyz gives information on the number of angular momentum of a
c**** particular type (ie. x,y or z) in a component of a shell
c**** eg. ixyz(2,14,5) asks the question "how many `y' contributions are
c**** there in the 14th component of the h shell?". the 14th component
c**** is hxyzzz. so the answer is 1
c
c**** if the number of components of y angular momentum and the number
c**** of components of z angular momentum are fed into the array, iarray
c**** then the result is the position of that component of any shell in
c**** list of components for that shell. eg. consider contribution with
c**** value of y angular momentum and two of z. this considers
c**** fyzz, gxyzz, and hxxyzz. because of the natural ordering of the sh
c**** components the postion in the list of these components is always 9
c**** and so iarray(1,2)=9
c
      implicit REAL  (a-h,o-z)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      data ispdf1/1,2,3,2,3,3,2,3,3,3,2,3,3,3,3,2,5*3
     &,1,1,1,2,2,3,4,4,5,6,7,7,8,9,10,11,11,12,13,14,15/
      data nvl/1,3,6,10,15,21/
      data ixyz/1,3*0,1,3*0,1,54*0,2,0,0,1,1,0,1,0,1,0,2,0,0,1,1
     &,0,0,2,45*0,3,0,0,2,1,0,2,0,1,1,2,0,1,1,1,1,0,2,0,3,0,0,2,1
     &,0,1,2,0,0,3,33*0,4,0,0,3,1,0,3,0,1,2,2,0,2,1,1,2,0,2,1,3,0
     &,1,2,1,1,1,2,1,0,3,0,4,0,0,3,1,0,2,2,0,1,3,0,0,4,18*0,5,0,0
     &,4,1,0,4,0,1,3,2,0,3,1,1,3,0,2,2,3,0,2,2,1,2,1,2,2,0,3,1,4,0
     &,1,3,1,1,2,2,1,1,3,1,0,4,0,5,0,0,4,1,0,3,2,0,2,3,0,1,4,0,0,5/
      data iarray/1,2,4,7,11,16,3,5,8,12,17,0,6,9,13,18,0,0
     &,10,14,19,0,0,0,15,20,0,0,0,0,21,0,0,0,0,0/
      end
      subroutine drcdab(zmemc,ivl,jvl,kvl,lvl,ic,jc,kmin,kmax
     &,lmin,lmax,l1,l2,ext2
     &,d,ipos,jpos,kpos,lpos,i1,j1,k1)
c
c*****  routine calculates two particle density matrix for gradient
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
c
      dimension zmemc(icsize,ivl*jvl*kvl*lvl),d(nb,nb)
      logical ext2
c
c
      ijkl = 0
      do 80 i = 1 , ivl
         id = (ic-i1)*ivl + i + ipos
         do 70 j = 1 , jvl
            jd = (jc-j1)*jvl + j + jpos
            dij = d(jd,id)*4.0d0
            do 60 k = 1 , kvl
c
               do 50 l = 1 , lvl
                  icount = 0
                  ijkl = ijkl + 1
c
                  do 40 kc = kmin , kmax
                     lmin1 = l1
                     lmax1 = l2
                     if (ext2) lmax1 = kc
                     if (kc.eq.kmax) lmax1 = lmax
                     if (kc.eq.kmin) lmin1 = lmin
                     kd = (kc-k1)*kvl + k + kpos
                     dik = -1.0d0*d(kd,id)
                     djk = -1.0d0*d(kd,jd)
c
                     if (lx.eq.1) then
                        do 20 ld = lmin1 , lmax1
                           icount = icount + 1
                           zmemc(icount,ijkl) = d(ld,kd)*dij + d(ld,jd)
     +                        *dik + d(ld,id)*djk
 20                     continue
                     else
c      ld=l+lpos-lvl
                        ld = (lmin1-l1-1)*lvl + l + lpos
                        do 30 lc = lmin1 , lmax1
                           ld = ld + lvl
                           icount = icount + 1
                           zmemc(icount,ijkl) = d(ld,kd)*dij + d(ld,jd)
     +                        *dik + d(ld,id)*djk
 30                     continue
                     end if
 40               continue
 50            continue
 60         continue
 70      continue
 80   continue
c
      return
      end
      subroutine drcdrv(q,task)
      implicit REAL  (a-h,o-z)
      character*(*) task
INCLUDE(common/sizes)
*
* global driving routine for direct
*
      dimension q(*)
      common/iofile/iread,iwr
INCLUDE(common/funct)
      logical odiis
      common/scfopt/maxxx(4),accdi(2),odiis(4),dmpcut(2),en,etot,ehf
INCLUDE(common/infoa)
      common/igarg/imarg(40),ipar,ileft,iwhole
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     & coprm(mxprim,3),nad,nbd,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ned,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
INCLUDE(common/wrtd)
      common/timez/time1,time2,time3,time4,time5,time6
c
      logical glog1,glog2,glog3,glog5,glog6,gdbg,gmp2
     &,readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf,gmull
c*****  gen contains most of the stuff about the options
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold
     &,gscf,noptd,npoint,timeg,gmull,iswap,dshift
*
*allocate memory
*
      nav = lenwrd()
c     determine memory requirements
      i1 = 0
      i2 = i1 + num*num
      i3 = i2 + num*num
      ishsp = (48*mxshel-1) / nav + 1
      i4 = i3 + ishsp
      i5 = i4 + ishsp
      i6 = i5 + ishsp
      i7 = i6 + (nat*48-1) / nav + 1
      i8 = i7 + maxat
      need = i8
c
c     allocate memory
c
      i1 = igmem_alloc(need)
      istart = i1
      i2 = i1 + num*num
*      space
      i3 = i2 + num*num
*      iwksp1
      i4 = i3 + ishsp
*      iwksp2
      i5 = i4 + ishsp
*      labsh
      i6 = i5 + ishsp
*      isoc
      i7 = i6 + (nat*48-1) / nav + 1
*      nwauq
      need = i7 + maxat
*      fcon
      if (gdbg(1)) write (iwr,*) 'memory in drcdrv' , i1 , i2 , 
     +             i3 , i4 , i5 , i6 , i7 , need
* convert to direct format
      call drccon(q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),
     +            istart,num)
* return memory
      call gmem_free(i1)
c
      need = imarg(17) - imarg(1)
      i1 = igmem_alloc(need)
      if (i1.ne.imarg(1)) then
         write (iwr,*) 'error in direct-scf memory allocation'
         write (iwr,*) 'i1' , i1 , imarg(1)
         call caserr('error in dscf memory allocation')
      end if
      i2 = imarg(2)
*      p
      i3 = imarg(3)
*      gp
      i4 = imarg(4)
*      shalf
      i5 = imarg(5)
*      pold
      i6 = imarg(6)
*      fstore
      i7 = imarg(7)
*      hx
      i8 = imarg(8)
*      sx
*     i9 = imarg(11)
*      ilar
*     i10 = imarg(12)
*      wx
*     i11 = imarg(13)
*      wy
*     i12 = imarg(14)
*      x1
*     i13 = imarg(15)
*      x2
*     i14 = imarg(16)
*      d
*     i15 = imarg(17)
*      the end
      if (gdbg(1)) write (iwr,*) 'drcdrv: ' , task
*
* decide which option has been requested
*
      if (task.eq.'scf') then
* scf
         call drcscf(q(1),q(i1),q(i2),q(i3),q(i4))
         en = repen
         etot = energy
         ehf = energy - repen
*
      else if (task.eq.'gradients') then
*gradients
         call drcgrd(q(1),q(i1),q(i2),q(i3),q(i4),iwr)
c
         do 30 loop = 1 , nat
            iat = i3 - 1 + loop
            do 20 moop = 1 , 3
               egrad((loop-1)*3 + moop) = q(iat)
               iat = iat + nat
 20         continue
 30      continue
c
      else if (task.eq.'mp2') then
*mp2
         call drcmp2(q(1),q(i1),q(i2),q(i3),q(i4),iwr)
      else
         write (iwr,*) 'the task ' , task , ' is not valid'
         call caserr('invalid task specified')
      end if
*
* return memory
*
      call gmem_free(i1)
c
      return
      end
_EXTRACT(drcfun,mips4)
      subroutine drcfun(x,f,npt,dji,madd)
      implicit REAL  (a-h,o-z)
      dimension f(60),madd(20),dji(60)
      dimension xmax(20)
      data xmax/24.0d0,29.0d0,32.0d0,35.0d0,37.0d0,40.0d0,42.0d0,
     &   44.0d0,47.0d0,49.0d0,51.0d0,53.0d0,55.0d0,57.0d0,58.0d0,
     &   59.0d0,59.0d0,59.0d0,59.0d0,59.0d0/
      m = npt + npt - 1
      fact = 4.8d0 + 0.4d0*npt
      if (x.ge.fact) then
         x2 = 0.5d0/x
         xinv = x2 + x2
         f1 = dsqrt(1.5707963267949d0*x2)
         if (x.ge.xmax(m)) then
c...... very high argument
            f(1) = f1
            if (m.gt.0) then
               do 20 i = 1 , m
                  f(i+1) = f(i)*x2
                  x2 = x2 + xinv
 20            continue
               return
            end if
         else
c.....  high argument
            x3 = dexp(-x)
            if (x.gt.21.6d0) then
               f1 = f1 - x3*x2
            else if (x.gt.18.2d0) then
               f1 = f1 - x3*x2*(1.0d0-x2)
            else if (x.gt.(12.0d0+0.1d0*npt)) then
               f1 = ((1.9623264149430d-1*xinv-4.9695241464490d-1)
     +              *xinv-6.0156581186481d-5)*x3 + f1
            else if (x.gt.(9.2d0+0.2d0*npt)) then
               f1 = (((-1.8784686463512d-1*xinv+2.2991849164985d-1)
     +              *xinv-4.9893752514047d-1)*xinv-2.1916512131607d-5)
     +              *x3 + f1
            else
               f1 = f1 +
     +              ((((((4.6897511375022d-1*xinv-6.9955602298985d-1)
     +              *xinv+5.3689283271887d-1)*xinv-3.2883030418398d-1)
     +              *xinv+2.4645596956002d-1)*xinv-4.9984072848436d-1)
     +              *xinv-3.1501078774085d-6)*x3
            end if
            f(1) = f1
            if (m.gt.0) then
               x23 = x2*x3
               do 30 i = 1 , m
                  f(i+1) = f(i)*x2 - x23
                  x2 = x2 + xinv
 30            continue
            end if
         end if
c.....  low argument
      else if (x.lt.1.0d-10) then
         j = m + 1
         do 40 i = 1 , j
            f(i) = dji(i)
 40      continue
         return
      else if (x.lt.1.0d-5) then
         j = m + 1
         do 50 i = 1 , j
            f(i) = dji(i) - x*dji(i+1)
 50      continue
         return
      else
         x2 = x + x
         i = x2 + madd(m)
         j = i
         x3 = dexp(-x)
         f(j+1) = 0.0d0
         do 60 k = 1 , i
            f(j) = (f(j+1)*x2+x3)*dji(j)
            j = j - 1
 60      continue
         return
      end if
      return
      end
_ENDEXTRACT
      subroutine drcgrd(core,fcon,p,gp,shalf,iwr)
      implicit REAL  (a-h,o-z)
      character *4 char
c
c calculates gradients
c
INCLUDE(common/sizes)
INCLUDE(common/statis)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
c
c*** the array d is the temporary array used by 1-electron ints.
c
      dimension shalf(nb,nb),fcon(nb,nb)
      dimension p(nb,nb),gp(nb,nb)
      dimension core(*)
c     dimension a(maxorb)
      common/worksp/adiis(286)
INCLUDE(common/wrtd)
      common/newsav/e(maxorb)
      common/timez/time1,time2,time3,time4,time5,time6
c
      logical glog1,glog2,glog3,glog5,glog6,gdbg,gmp2
     &,readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf,gmull
c*****  gen contains most of the stuff about the options
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold
     &,gscf,noptd,npoint,timeg,gmull,iswap,dshift
c
c*******zipftd contains the pre-factor test data and
c*******zijkld contains the ijkl test data whilst
c*******zjpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
      dimension char(21)
      data char/'ssss','psss','psps','ppss','ppps','pppp'
     &,'dsss','dsps','dspp','dsds','dpss','dpps','dppp'
     &,'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
c
c**** doing a gradient
c
      no = ne/2
c
c***** preparing integral test information
c
      do 30 j = 1 , 2
         do 20 i = 1 , 21
            zipftd(i,j) = 1.0d0
            zijkld(i,j) = 1.0d0
 20      continue
 30   continue
      xfact = 1.0d0
      if (iter.gt.1) xfact = dsqrt(dfloat(iter-1))*xfact
      do 40 i = 1 , 21
         prefac(i) = prefac(i)*xfact
         xinner(i) = xinner(i)*xfact
 40   continue
c
      if (.not.gopt .and. gdbg(8)) then
         imax = 21
         if (nd.eq.0) imax = 6
         write (iwr,6050)
         write (iwr,6020)
         do 50 i = 1 , imax
            write (iwr,6030) char(i) , prefac(i) , xinner(i)
 50      continue
      end if
      if (gdbg(8)) then
         write (iwr,6050)
         write (iwr,6050)
      end if
      iaccur = 4
c
c***** read in old fock matrix to generate new mos
      call rdedx(gp,nb*nb,idblk(5+inddir),numdir)
c
      call scfit(gp,shalf,p,fcon,e)
c
c***** multiply all mos associated with dxy,dxz,dyz orbitals by
c***** root(3)
c
      sqr3 = dsqrt(3.0d0)
      icmp2 = 0
      do 70 i = ns + np3 + 1 , nb
         icmp2 = icmp2 + 1
         if (icmp2.eq.6) icmp2 = 0
         if (icmp2.ne.1 .and. icmp2.ne.4 .and. icmp2.ne.0) then
            do 60 j = 1 , no
               gp(i,j) = gp(i,j)*sqr3
 60         continue
         end if
 70   continue
c
c*****   use e.v.s  in gp to find density - put in fcon
c***** also form matrix r to find orbital contribution to gradient
c
      call vclr(fcon,1,nb*nb)
      call vclr(p,1,nb*nb)
      do 100 i = 1 , nb
         do 90 j = 1 , nb
            do 80 k = 1 , no
               fcon(i,j) = fcon(i,j) + 2.0d0*gp(i,k)*gp(j,k)
               p(i,j) = p(i,j) + 2.0d0*gp(i,k)*gp(j,k)*e(k)
 80         continue
 90      continue
 100  continue
c
c**** now call the first one electron contribution to gradient
c
      if (gdbg(8)) then
         write (iwr,6050)
         write (iwr,6060)
         if (.not.gopt) call timit(0)
         write (iwr,6050)
      end if
      call def1dr(core,p,fcon,gp,iwr)
      call dengrd(p,fcon)
      if (gdbg(8)) then
         write (iwr,6050)
         write (iwr,6070)
      end if
      call timana(6)
      call cpuwal(begin,ebegin)
      if (.not.gopt) call timit(0)
      if (gdbg(8)) write (iwr,6050)
      call defgrd(core,p,fcon,gp,iwr)
      if (gdbg(8)) then
         write (iwr,6050)
         if (.not.gopt) call timit(0)
         write (iwr,6050)
         write (iwr,6040)
         do 110 i = 1 , imax
_IF1()c      write(iwr,*) zipftd(i,1),zipftd(i,2),zijkld(i,1),zijkld(i,2)
            dperc1 = zipftd(i,2)/zipftd(i,1)
            dperc2 = zijkld(i,2)/zijkld(i,1)
            dperc3 = dperc1*dperc2
            dperc1 = dperc1*100.0d0
            dperc2 = dperc2*100.0d0
            dperc3 = dperc3*100.0d0
            write (iwr,6010) char(i) , dperc1 , dperc2 , dperc3
 110     continue
      end if
      call timana(7)
      return
c
 6010 format (2x,a4,10x,f6.2,16x,f6.2,14x,f6.2)
 6020 format (/'      information on integral thresholds'//
     +        ' routine          prefactor              inner'/)
 6030 format (2x,a4,10x,d11.3,10x,d11.3)
 6040 format (1x,'routine',5x,'prefactor  pass %',6x,'inner pass %',8x,
     +        'total pass %')
_IF1()c1703  format(//)
 6050 format (1x)
 6060 format (1x,'beginning gradient calculation'/)
 6070 format (1x,'finished the one electron contribution to gradient'/)
      end
      subroutine drcmp2(core,fcon,p,gp,shalf,iwr)
c
      implicit REAL  (a-h,o-z)
c
c direct mp2 calculation
INCLUDE(common/sizes)
      common/igarg/igxx(40),ipar,ileft,iwhole
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
c
c
      dimension shalf(nb,nb),fcon(nb,nb),p(nb,nb),gp(nb,nb)
      dimension core(*)
c     dimension a(maxorb)
      common/worksp/adiis(286)
INCLUDE(common/wrtd)
      common/newsav/e(maxorb)
      common/timez/time1,time2,time3,time4,time5,time6
c
      logical glog1,glog2,glog3,glog5,glog6,gdbg,gmp2
     &,readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf,gmull
c*****  gen contains most of the stuff about the options
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold
     &,gscf,noptd,npoint,timeg,gmull,iswap,dshift
c
c*******zipftd contains the pre-factor test data and
c*******zijkld contains the ijkl test data whilst
c*******zjpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
c     character *4 char
c     dimension char(21)
c     data char/'ssss','psss','psps','ppss','ppps','pppp'
c    &,'dsss','dsps','dspp','dsds','dpss','dpps','dppp'
c    &,'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
c
c***** doing an mp2 calculation, obtain molecular orbitals
c***** and then call the integral routine defmp2 to do the
c***** calculation
c
      if (.not.(readmo)) then
c        ntri = nb*(nb+1)/2
         call rdedx(fcon,nb*nb,idblk(5+inddir),numdir)
c
         call scfit(fcon,shalf,p,gp,e)
      end if
c
c***** the molecular orbitals are in w
c
c***** find symmetries of the mo's
c
      call rdedx(gp,nb*nb,idblk(2),numdir)
      call trianc(gp,p,nb,nb)
c
      no = ne/2
      nvir = mp2fv - no
      no = no - mp2fo + 1
      iaccur = 4
      call sytype(fcon,e,p,mp2fo,mp2fv,nb,gp,iwr)
      call defmp2(core,fcon(1,mp2fo),e(mp2fo),no,nvir,iwr)
c
c
      return
      end
      subroutine drcscf(core,fcon,p,gp,shalf)
c
      implicit REAL  (a-h,o-z)
      character *4 char
      logical ocvged,outon
c
c direct scf driver
c
INCLUDE(common/sizes)
      logical ofull,odiisn
      parameter (nl=10)
INCLUDE(common/iofile)
      common/gen3/ofull(100),odiisn
INCLUDE(common/restar)
INCLUDE(common/runlab)
      common/igarg/igxx(40),ipar,ileft,iwhole
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
c
c
      dimension shalf(nb,nb),fcon(nb,nb),p(nb,nb),gp(nb,nb)
_IF(parallel)
      logical opg_root
INCLUDE(common/prints)
INCLUDE(common/nshel)
      common/pardir/next,nbcnt(256)
      dimension nhelp(256)
      dimension core(*)
_ELSE
      dimension core(*)
_ENDIF
      dimension a(maxorb)
      common/worksp/adiis(286)
INCLUDE(common/wrtd)
      common/newsav/e(maxorb)
      common/timez/timlim,time2,time3,time4,time5,time6
c
      logical glog1,glog2,glog3,glog5,glog6,gdbg,gmp2
     &,readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf,gmull
c*****  gen contains most of the stuff about the options
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold
     &,gscf,noptd,npoint,timeg,gmull,iswap,dshift
c
c*******zipftd contains the pre-factor test data and
c*******zijkld contains the ijkl test data whilst
c*******zjpftd contains the second  pre-factor test info.
c*******prefac contains the thresholds for the pre-factor tests
c*******xinner contains the thresholds for the inner tests.
c
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
INCLUDE(common/atmol3)
INCLUDE(common/mapper)
c
INCLUDE(common/harmon)
c
      character*10 charwall
      dimension char(21)
      data char/'ssss','psss','psps','ppss','ppps','pppp'
     &,'dsss','dsps','dspp','dsds','dpss','dpps','dppp'
     &,'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/
      data maxf /4/
c
c********** gp contains    shalf (the orthogonaliser)
c********** now gp changed to old HF orbitals (jvl,1997)
c
      outon = nprint.ne. - 5
      skale = 2.0d0
      if (zruntp.eq.'scf') skale = 1.1d0
      if (outon) write (iwr,6100)
      restar = irest.ne.0
      if (restar) write (iwr,6110) irest
      tim0 = cpulft(1)
_IF(parallel)
      l2 = nshell*(nshell+1)/2
      nodet = (l2-1)/ipg_nnodes() + 1
      if (nodet.ge.40) then
       ichunk = nodet / 40
      else
       ichunk = 1
      endif
      if(oprint(57)) then
        print*,'nshell, nnodes, nodet, ichunk = ',
     +          nshell, ipg_nnodes(), nodet, ichunk
      endif
      call pg_dlbchunk(ichunk,.true.)
_ENDIF
*
* read in the overlap matrix
*
      l3 = nb * nb
      l0 = newbas0
      call rdedx(p,l3,idblk(2),numdir)
*
* read in 1-particle core matrix
*
      call rdedx(fcon,l3,idblk(1),numdir)
      do 140 i = 1 , nb
         do 150 j = 1 , i
            p(j,i) = p(i,j)
            fcon(j,i) = fcon(i,j)
 150     continue
 140  continue
      if (gdbg(4)) then
         write (iwr,6050)
         call foutd(fcon,nb)
         write (iwr,6060)
         call foutd(p,nb)
      end if
c...
c...  get the start orbitals as orthogonalisation matrix (level shift)
c...
      nddd = nprint
      nprint = -5
      call getq(gp,gp,gp,nb,nb,0,ieig,ipop,mouta,'newscf')
      nprint = nddd
      dshift = gapa1
c
c     we must still orthogonalise the input vectors, to cure
c     geometry optimisations screwing on the 2nd point!
c
*
* generate overlap in symmetry adapted basis
*
      call tranpd(p,fcon,nb)
      do 120 i = 1 , nb
         do 130 j = 1 , i
            fcon(j,i) = fcon(i,j)
 130     continue
 120  continue
      call trianc(fcon,shalf,nb,nb)
c
      call mult2(gp,fcon,shalf,nb,nb,nb)
      call orfogd(gp,gp,fcon,p,iky,ilifq,l0,nb)
      call wrt3(gp,l3,idblk(3),numdir)
*
* fcon   has the 1-particle core hamiltonian
* p      is the overlap matrix
*
c
      if (restar) then
         call rdedx(adiis,123,idblk(27),numdir)
         iter = nint(adiis(122))
         icount = nint(adiis(123))
         xfact = 1.0d0/dsqrt(dfloat(iter-1))
         do 20 i = 1 , 21
            prefac(i) = prefac(i)*xfact
            xinner(i) = xinner(i)*xfact
 20      continue
      else
         iter = 1
         icount = 0
      end if
c
c read in the density matrix
c
      call rdedx(p,l3,idblk(26),numdir)
c
      repen = 0.0d0
      do 40 i = 2 , na
         do 30 j = 1 , i - 1
            r2 = (coord(1,i)-coord(1,j))**2 + (coord(2,i)-coord(2,j))
     +           **2 + (coord(3,i)-coord(3,j))**2
            repen = repen + charge(i)*charge(j)/dsqrt(r2)
 30      continue
 40   continue
      if (outon) write (iwr,6120) repen
c     write(iwr,1701) repen
      if (gdbg(5)) then
         write (iwr,6070)
         call foutd(p,nb)
      end if
      if (maxit.gt.100) maxit = 100
      if (outon) write (iwr,6130) maxit , concrit
      call timrem(tlefts)
      amag = 0.0d0
      elec0 = 0.0d0
c
c*****the scf iteration return point
c
 50   continue
      call flushn(iwr)
      if (iter.eq.1) then
         call vclr(fcon,1,l3)
      else
         call rdedx(fcon,l3,idblk(4),numdir)
      end if
*
*     use full fock build every maxf iteration
*     stems from convergence problems with delta only
*     approach
*
      if (.not.ofull(iter)) then
         if (mod (iter,maxf).eq.0) then
            ofull(iter) = .true.
         else
            ofull(iter)=  .false.
         endif
      endif
c
      call dnsprm(p,fcon,gp,outon)
c
c***** the maximum change in density matrix elements is in dcmax
c***** (this is the convergence touch paper)
c
      conv = dcmax
      ocvged = conv.lt.concrit
      if (ocvged) then
         irest = 0
         if (outon) write (iwr,6140)
      end if
      if (.not.ocvged .and. iter.eq.maxit) then
         write (iwr,6150)
         energy = 0.0d0
         irest = 3
         adiis(122) = dfloat(iter)
         adiis(123) = dfloat(icount)
         call wrt3(adiis,123,idblk(27),numdir)
         elec = -repen
      end if
      if (ocvged .or. iter.eq.maxit) then
*
* convergence has been achieved
*
         timeg = cpulft(1)
         elec = energy - repen
         write (iwr,6160) iter , timeg , charwall(), 
     1                    elec , repen , energy , conv
         if (gdbg(7)) then
            write (iwr,6030)
            write (iwr,*) '**************************************'
            write (iwr,*) 'integral test data over all iterations'
            write (iwr,*) '**************************************'
            write (iwr,6040)
            write (iwr,6020)
            imax = 21
            if (nd.eq.0) imax = 6
            do 60 i = 1 , imax
               dperc1 = zipfal(i,2)/zipfal(i,1)
               dperc2 = zijkal(i,2)/zijkal(i,1)
               dperc3 = dperc1*dperc2
               dperc1 = dperc1*100.0d0
               dperc2 = dperc2*100.0d0
               dperc3 = dperc3*100.0d0
               write (iwr,6010) char(i) , dperc1 , dperc2 , dperc3
 60         continue
            write (iwr,6030)
         end if
c        no = ne/2
         if (.not.gopt .and. gdbg(7)) then
            write (iwr,*) '****************************************'
            dperc1 = zifall/(zifall+ziffew)
            dperc2 = ziffew/(zifall+ziffew)
            dperc1 = dperc1*100.0d0
            dperc2 = dperc2*100.0d0
            write (iwr,*) 'the info on fquick'
            write (iwr,*) 'percentage of arguments for which the fm(t)'
            write (iwr,*) 'were measured ' , dperc1
            write (iwr,*) 'percentage evaluated by the extrapolation to'
            write (iwr,*) 'infinity formula ' , dperc2
            write (iwr,*) '****************************************'
            write (iwr,6030)
         end if
         call putqd(fcon,shalf,p,gp,e,a,energy,nb,nprint)
         if (outon) then
            cpu = cpulft(1)
            write (iwr,6170) cpu ,charwall()
         end if
         return
      else
         call wrt3(fcon,l3,idblk(4),numdir)
c
c*****  adjust the integral test thresholds for this iteration
c*****  thr(i)=thr(1)/sqrt(i). also adjust accuracy
c ****  initialise iaccur correctly (jk)
c
        iaccur = 4
         if (iter.ne.1) then
            xfact = dfloat(iter)
            xfact = dsqrt(xfact-1.0d0)/dsqrt(xfact)
            do 70 i = 1 , 21
               prefac(i) = prefac(i)*xfact
               xinner(i) = xinner(i)*xfact
 70         continue
*!!!        write(iwr,*)'acvary is ',acvary
            if (acvary) then
               contes = conv*1.0d-5/concrit
               iaccur = 4
               if(.not.ofull(iter)) then
                if (contes.lt.2.0d0) iaccur = 3
                if (contes.lt.4.0d-2) iaccur = 2
                if (contes.lt.6.0d-4) iaccur = 1
               end if
            end if
         end if
_IF1()*        write(iw,*)' iaccur, contes ',iaccur,contes
c
c*********************************************
c***** the two electron integral call********
c*********************************************
c
         call vclr(fcon,1,l3)
_IF(parallel)
         call setsto(ipg_nnodes(),0,nbcnt)
         call defunc(core,fcon,gp,p)
c
c***********two electron integrals done
c
         call pg_igop(333,nbcnt,ipg_nnodes(),'+')
         if (opg_root().and.oprint(57)) 
     +               call outive(nbcnt,ipg_nnodes(),'nbcnt')
_ELSE
         call defunc(core,fcon,gp,p)
c
c***********two electron integrals done
c
_ENDIF
         if (iter.eq.1) then
            call vclr(p,1,l3)
         else
            call rdedx(p,l3,idblk(5),numdir)
         end if
c
c******** form complete fock matrix iteratively and by
c******** symmetrisation if necessary
c
         call fcontr(fcon,p,gp,nb)
c
         call wrt3(p,l3,idblk(5),numdir)
c
c
         call rdedx(p,l3,idblk(4),numdir)
         call rdedx(gp,l3,idblk(1),numdir)
c
c***** establish the energy and convergence criterion
         call fbuild(p,gp,fcon)
c
c
         elec = energy - repen
         de = elec - elec0
         elec0 = elec
         tim1 = cpulft(1)
         delt = tim1 - tim0
         tim0 = tim1
c     write(iwr,90)iter,energy,conv
         if (outon .or. (maxit-iter.lt.10)) write (iwr,6180) iter ,
     +       energy , elec , de , conv , amag , delt , tim1 , dshift
         call timit(0)
c
c...     level shift determining
c
         dshift = gapa1
         if (iter.ge.ibrk) dshift = gapa2
         if (conv.lt.concrit*100.0d0.and.iter.ge.ibrk) dshift = 0.0d0
c
c     write(iwr,1704)
         if (gdbg(6)) then
            write (iwr,6080)
            call foutd(fcon,nb)
            write (iwr,6090)
            call foutd(p,nb)
         end if
         if (gdbg(2)) then
            write (iwr,6020)
            imax = 21
            if (nd.eq.0) imax = 6
            do 80 i = 1 , imax
               dperc1 = zipftd(i,2)/zipftd(i,1)
               dperc2 = zijkld(i,2)/zijkld(i,1)
               dperc3 = dperc1*dperc2
               dperc1 = dperc1*100.0d0
               dperc2 = dperc2*100.0d0
               dperc3 = dperc3*100.0d0
               write (iwr,6010) char(i) , dperc1 , dperc2 , dperc3
 80         continue
            write (iwr,6040)
         end if
         do 100 j = 1 , 2
            do 90 i = 1 , 21
               zipfal(i,j) = zipfal(i,j) + zipftd(i,j)
               zijkal(i,j) = zijkal(i,j) + zijkld(i,j)
               zipftd(i,j) = 1.0d0
               zijkld(i,j) = 1.0d0
 90         continue
 100     continue
c
c***** if required, write out density matrix
c
*      if (gdma) then
*         call denwr(p,gp)
*      endif
c
         iter = iter + 1
c
c         go into diis
c
c        if (conv.lt.thresh .or. iter.gt.3) then
         if (conv.lt.thresh ) then
            thresh = 1.0d+10
            inddir = mod(icount,nl) + 1
            nerr = nl
            if (icount.lt.nl) nerr = icount + 1
            call diisn(fcon,p,gp,inddir,nerr,amag)
            icount = icount + 1
         end if
         call scfit(fcon,shalf,p,gp,e)
c
c***** if not enough time to complete another iteration stop
c
         call timrem(tlefti)
         if (tlefti.lt.skale*(tlefts-tlefti)/iter) then
            irest = 3
            write (iwr,6190)
            call texit(0,irest)
            adiis(122) = dfloat(iter)
            adiis(123) = dfloat(icount)
            call wrt3(adiis,123,idblk(27),numdir)
            call putqd(fcon,shalf,p,gp,e,a,energy,nb,nprint)
            return
         end if
         go to 50
      end if
c  the arrays used to store information from saving have
c  dimensions greater than n1=nsp*(nsp+1)/2, n2=nsp*npp, or
c  n3=npp*(npp+1)/2 for the ss, ps, or pp arrays respectivel
c
_IF1()c90    format(1x,i3,'    energy is ',d20.12,'       conv',d20.12)
 6010 format (2x,a4,10x,f6.2,16x,f6.2,14x,f6.2)
_IF1()c1695  format(1x,'calculation of 1-electron integrals complete')
_IF1()c1701  format(1x,'nuclear repulsion energy',10x,d20.12//)
 6020 format (1x,'routine',5x,'prefactor  pass %',6x,'inner pass %',8x,
     +        'total pass %')
 6030 format (//)
 6040 format (1x)
 6050 format (1x,'the hcore matrix'//)
 6060 format (1x,'the overlap matrix'//)
 6070 format (1x,'the guess density matrix'//)
 6080 format (1x,'the fock matrix'//)
 6090 format (1x,'the density matrix'//)
 6100 format (//16x,39('*')/16x,
     +        'closed-shell rhf direct-scf calculation',' - v2.0'/16x,
     +        39('*')//)
 6110 format (' restart parameter in direct-scf ',i2)
 6120 format (/,' ----- nuclear energy ----- = ',f20.12)
 6130 format (/15x,'convergence data'/15x,16('=')
     +        ///' maximum number of iterations = ',
     +        i10/' convergence criterion        = ',f10.7//1x,112('=')
     +        /3x,'cycle',10x,'total',5x,'electronic',8x,'e conv.',9x,
     +        'tester',11x,'diis',4x,'del(t)',6x,'time',5x,'shift',
     +        /17x,'energy',9x,'energy'/1x,112('='))
 6140 format (/20x,17('-')/20x,'density converged'/20x,17('-'))
 6150 format (/20x,30('-')/20x,'excessive number of iterations'/20x,
     +        30('-'))
 6160 format (//10x,14('-')/10x,'final energies   after',i4,
     +       ' cycles at ',f8.2,' seconds',a10,' wall'/10x,14('-')//10x,
     +        'electronic energy          ',f15.8/10x,
     +        'nuclear energy             ',f15.8/10x,
     +        'total energy               ',f15.8/10x,
     +        'maximum change in densities',f15.8/)
 6170 format (//' end of direct-scf at ',f8.2,' seconds'
     +        ,a10,' wall'//1x,104('-')//)
 6180 format (1x,i5,2x,4f15.8,3x,d12.2,3f10.3)
 6190 format (/10x,26('*')/10x,'*** warning ***'/10x,
     +        'scf has not converged '/10x,
     +        'this job must be restarted'/10x,26('*')/)
      end
      subroutine ds1(xijk,xija,xijp1,xijp2,xijp3,i1,i2
     &,j1,j2,xt1,xt2,wx
     &,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ccc1)
c
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (d/s) and (d/h/s) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      dimension xija(n11),xijk(n11)
     &,s1(n11*6),h1(n11*6),sx(nb,nb),hx(nb,nb),ccc1(n11,3)
     &,lar(n11),f(7*n11),t(n11),xt1(n14*6,n13*6),xt2(n16,n15)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
     &,xijp1(n11),xijp2(n11),xijp3(n11)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11*6)
      icount = 0
c
      do 30 i = imin , imax
         do 20 j = jmin , jmax
            icount = icount + 1
            ccc1(icount,1) = xijp1(icount) - coprm(i,1)
            ccc1(icount,2) = xijp2(icount) - coprm(i,2)
            ccc1(icount,3) = xijp3(icount) - coprm(i,3)
 20      continue
 30   continue
      do 60 k = 1 , na
c
         do 40 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 40      continue
c
         call fquik2(f,f,lar,t,3,n11)
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 50 icount = 1 , n11
            index1 = n11 + icount
            index2 = n11 + index1
            index3 = n11 + index2
            index4 = n11 + index3
            index5 = n11 + index4
            x4 = (f(icount)-f(index1))/(xija(icount)+xija(icount))
            x1a = xijp1(icount) - coord(1,k)
            x1b = xijp2(icount) - coord(2,k)
            x1c = xijp3(icount) - coord(3,k)
c
c     ***** xx1 has the k dependant factors of (pi/v/s)m=0 ***** 
c     ***** xx2 has the k dependant factors of (pi/v/s)m=1 ***** 
c
            x1x = ccc1(icount,1)*f(icount) - x1a*f(index1)
            x2x = ccc1(icount,1)*f(index1) - x1a*f(index2)
            x1y = ccc1(icount,2)*f(icount) - x1b*f(index1)
            x2y = ccc1(icount,2)*f(index1) - x1b*f(index2)
            x1z = ccc1(icount,3)*f(icount) - x1c*f(index1)
            x2z = ccc1(icount,3)*f(index1) - x1c*f(index2)
c
c**** dxx
            h1(icount) = h1(icount) + (ccc1(icount,1)*x1x-x1a*x2x+x4)
     +                   *charge(k)
c**** dxy
            h1(index1) = h1(index1) + (ccc1(icount,2)*x1x-x1b*x2x)
     +                   *charge(k)
c**** dxz
            h1(index2) = h1(index2) + (ccc1(icount,3)*x1x-x1c*x2x)
     +                   *charge(k)
c**** dyy
            h1(index3) = h1(index3) + (ccc1(icount,2)*x1y-x1b*x2y+x4)
     +                   *charge(k)
c**** dyz
            h1(index4) = h1(index4) + (ccc1(icount,3)*x1y-x1c*x2y)
     +                   *charge(k)
c**** dzz
            h1(index5) = h1(index5) + (ccc1(icount,3)*x1z-x1c*x2z+x4)
     +                   *charge(k)
c
 50      continue
 60   continue
c
      icount = 0
      do 80 i = imin , imax
         do 70 j = jmin , jmax
            icount = icount + 1
            index1 = n11 + icount
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(index1) = pe(i)
            f(icount) = pe(i)*pe(j)
 70      continue
 80   continue
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 90 icount = 1 , n11
         index1 = n11 + icount
         index2 = n11 + index1
         index3 = n11 + index2
         index4 = n11 + index3
         index5 = n11 + index4
c
c     ***** x3 contains (s/s) and xx1 has (pi/s) *******
c
         x1 = 1.0d0/xija(icount)
         x2 = f(icount)*x1
         x3 = x2*(7.0d0-2.0d0*x2*t(icount))
         x4 = xijk(icount)*pi14*dsqrt(x1)
         x6 = x4*x1*0.5d0
         x5 = x2*(2.0d0*x6+x4/f(index1))
         x7 = xijk(icount)/pi14
         x1x = ccc1(icount,1)*x4
         x1y = ccc1(icount,2)*x4
         x1z = ccc1(icount,3)*x4
c
c**** dxx
         s1(icount) = ccc1(icount,1)*x1x + x6
         h1(icount) = x3*s1(icount) - x7*h1(icount) - x5
c**** dxy
         s1(index1) = ccc1(icount,1)*x1y
         h1(index1) = x3*s1(index1) - x7*h1(index1)
c**** dxz
         s1(index2) = ccc1(icount,1)*x1z
         h1(index2) = x3*s1(index2) - x7*h1(index2)
c**** dyy
         s1(index3) = ccc1(icount,2)*x1y + x6
         h1(index3) = x3*s1(index3) - x7*h1(index3) - x5
c**** dyz
         s1(index4) = ccc1(icount,2)*x1z
         h1(index4) = x3*s1(index4) - x7*h1(index4)
c**** dzz
         s1(index5) = ccc1(icount,3)*x1z + x6
         h1(index5) = x3*s1(index5) - x7*h1(index5) - x5
 90   continue
c
c*****  put the result in the matrices to be transformed ***
c
      icount = 0
      do 120 ic = 1 , 6
         do 110 i = 1 , n14
            ip = (i-1)*6
            do 100 j = 1 , n16
               icount = icount + 1
               wx(ip+ic,j) = s1(icount)
               wy(ip+ic,j) = h1(icount)
 100        continue
 110     continue
 120  continue
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16,wz,1,ipar2,n14*6,n16,n15)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wx,ipar2,1,n15,n14*6,n13*6)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16,wz,1,ipar2,n14*6,n16,n15)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*6,wy,ipar2,1,n15,n14*6,n13*6)
c
      ix1 = 6*i1 - 3*np - ns*5 - 5
      ix2 = ix1 + n13*6 - 1
      ic = 0
      do 140 i = ix1 , ix2
         ic = ic + 1
         jc = 0
         do 130 j = j1 , j2
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 130     continue
 140  continue
c
      return
      end
      subroutine exchng(zmem,xijk,xija,xijt,xijp1,xijp2
     &,xijp3,nij,ext1,i1,i2,j1,j2)
c
c***** routine calculates exchange integrals (ds|ds), (dp|dp) or (dd|dd)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c**** the fixed parameters in common/indexd/ are explained in the
c**** block data statement (in bdata.f?)
c**** the memory locations and requirements were established in
c**** cidex and are stored in common/cidex/
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto)
     &,xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,ab(3),naind(3),ncind(3),xijt(nij)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cidex/intgrl(8,100),ihex(8)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      common/defunk/cobas(maxorb,3),nd
      integer abind,abindk,abindl
      logical ext1,ext5
c
      sq3 = dsqrt(3.0d0)
      jmax = j2
      icount = 0
      icorig = 0
      ic1 = 0
c
c**** the main contracted loop over i and j shells
c
      do 260 ic = i1 , i2
         if (ext1) jmax = ic
         do 250 jc = j1 , jmax
            ext5 = ext1 .and. ic.eq.jc
            ab(1) = cobas(ic,1) - cobas(jc,1)
            ab(2) = cobas(ic,2) - cobas(jc,2)
            ab(3) = cobas(ic,3) - cobas(jc,3)
c           rab = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
            jpmax = nc(jc)
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
            do 30 iprm = 1 , nc(ic)
               if (ext5) jpmax = iprm
               do 20 jprm = 1 , jpmax
                  icount = icount + 1
                  ic1 = ic1 + 1
                  zmem(ic1,1) = xijp1(icount) - cobas(ic,1)
                  zmem(ic1,2) = xijp2(icount) - cobas(ic,2)
                  zmem(ic1,3) = xijp3(icount) - cobas(ic,3)
                  zmem(ic1,4) = ab(1)
                  zmem(ic1,5) = ab(2)
                  zmem(ic1,6) = ab(3)
 20            continue
 30         continue
c
c***** at this point just counting out the number of integrals
c
            if (jc.lt.jmax) then
               ic2 = ic1 + nc(ic)*nc(jc+1)
            else
               ic2 = ic1 + nc(ic+1)*nc(j1)
            end if
c
            if (.not.(ic2.le.ipsize .and. (ic.ne.i2 .or. jc.ne.j2)))
     +          then
c
c***** now we have a batch of integrals find primitive ints
c
c***** about to form the (ssss)m=0 integrals in zmem(12)
c***** zmem(1,8) contains (1/2)*(1/a+b)
c***** zmem(1,1-3) has pi-ai  zmem(1,4-6) has ai-bi
c***** zmem(1,7) is workspace. zmem(9-11) unused.
c
               icor1 = icorig + 1
               icor2 = icorig + ic1
               mc = 0
               do 40 jcount = icor1 , icor2
                  mc = mc + 1
                  zmem(mc,8) = 0.5d0*(1.0d0/xija(jcount))
                  zmem(mc,7) = xijk(jcount)*xijk(jcount)
     +                         *dsqrt(zmem(mc,8))
 40            continue
c
               do 60 m = 1 , mx
                  mc = 11 + m
                  factor = 1.0d0/(dfloat(2*m-1))
                  do 50 jcount = 1 , ic1
                     zmem(jcount,mc) = zmem(jcount,7)*factor
 50               continue
 60            continue
c
c*** loop over number of distinct integral classes
c*** needed to form the target integrals
c
               memcur = 11 + mx
               do 140 n = 1 , ninty
                  ityp = intgrl(1,n)
                  ktyp = intgrl(2,n)
                  ityp1 = ityp - 1
                  ktyp1 = ktyp - 1
                  ityp2 = ityp - 2
                  ktyp2 = ktyp - 2
                  ivl = nvl(ityp)
                  kvl = nvl(ktyp)
                  if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
                  if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
                  nrecur = intgrl(3,n)
c
                  mem1 = intgrl(4,n)
                  mem2 = intgrl(5,n)
                  mem3 = intgrl(6,n)
                  mem4 = intgrl(7,n)
                  mem5 = intgrl(8,n)
c
c now:
c loop over number of indices associated with i position
                  do 130 i = 1 , ivl
                     nr = nrecur
                     if (nr.lt.5) then
c***** reduction is at the i position
c
                        index = ispdf1(i,1)
                        natred = ispdf1(i,2)
                        itemp = (natred-1)*kvl
                        jmem1 = itemp + mem1
                        jmem5 = mem5 + (natred-1)*kvl1
                        if (nr.ne.1 .and. nr.ne.3) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                           nai = ixyz(index,natred,ityp2)
                           if (nai.eq.0) then
                              nr = nr - 1
                              go to 70
                           end if
                           naind(2) = ixyz(2,natred,ityp2) + 1
                           naind(3) = ixyz(3,natred,ityp2) + 1
                           naind(index) = naind(index) - 1
                           idoubr = iarray(naind(2),naind(3))
                           itemp2 = (idoubr-1)*kvl
                           jmem3 = itemp2 + mem3
                           jmem4 = itemp2 + mem4
                        end if
                     else
c
c**** reduction is at the k position
                        itemp = (i-1)*kvl1
                        jmem1 = itemp + mem1
                        itemp = (i-1)*kvl2
                        jmem3 = itemp + mem3
                        jmem4 = itemp + mem4
                        if (nr.ne.5 .and. nr.ne.6) then
                           naind(1) = ixyz(1,i,ityp1) + 1
                           naind(2) = ixyz(2,i,ityp1) + 1
                           naind(3) = ixyz(3,i,ityp1) + 1
                        end if
                     end if
c
c now:
c loop over number of indices associated with k position
c
 70                  do 120 k = 1 , kvl
                        memcur = memcur + 1
                        if (nr.lt.5) then
                           imem1 = jmem1 + k
                        else
                           index = ispdf1(k,1)
                           natred = ispdf1(k,2)
                           imem1 = jmem1 + natred
                        end if
                        go to (80,90,100,110,80,90,100,110) , nr
 80                     call rex5(zmem(1,memcur),zmem(1,index),
     +                            zmem(1,imem1),ic1)
                        go to 120
 90                     if (nr.lt.5) then
                           imem3 = jmem3 + k
                           imem4 = jmem4 + k
                        else
                           nai = ixyz(index,natred,ktyp2)
* jkendrick
                           if (nai.eq.0) go to 80
                           ncind(2) = ixyz(2,natred,ktyp2) + 1
                           ncind(3) = ixyz(3,natred,ktyp2) + 1
                           ncind(index) = ncind(index) - 1
                           imem3 = jmem3 + iarray(ncind(2),ncind(3))
                           imem4 = jmem4 + iarray(ncind(2),ncind(3))
                        end if
                        call rex6(zmem(1,memcur),zmem(1,index),
     +                            zmem(1,imem1),nai,zmem(1,imem3),
     +                            zmem(1,8),zmem(1,imem4),ic1)
                        go to 120
 100                    if (nr.lt.5) then
                           nci = ixyz(index,k,ktyp1)
                           if (nci.eq.0) go to 80
                           naind(2) = ixyz(2,k,ktyp1) + 1
                           naind(3) = ixyz(3,k,ktyp1) + 1
                           naind(index) = naind(index) - 1
                           imem5 = jmem5 + iarray(naind(2),naind(3))
                           call rex7(zmem(1,memcur),zmem(1,index),
     +                               zmem(1,imem1),zmem(1,8),nci,
     +                               zmem(1,imem5),ic1)
                           go to 120
                        else
                           if (naind(index).eq.1) go to 80
                           naind(index) = naind(index) - 1
                           nci = naind(index)
                           imem5 = (iarray(naind(2),naind(3))-1)
     +                             *kvl1 + natred + mem5
                           naind(index) = naind(index) + 1
                        end if
                        call rex7(zmem(1,memcur),zmem(1,index),
     +                            zmem(1,imem1),zmem(1,8),nci,
     +                            zmem(1,imem5),ic1)
                        go to 120
 110                    if (nr.lt.5) then
                           nci = ixyz(index,k,ktyp1)
                           if (nci.eq.0) go to 90
                           naind(2) = ixyz(2,k,ktyp1) + 1
                           naind(3) = ixyz(3,k,ktyp1) + 1
                           naind(index) = naind(index) - 1
                           imem5 = jmem5 + iarray(naind(2),naind(3))
                           imem3 = jmem3 + k
                           jmem4 = jmem4 + k
                           call rex8(zmem(1,memcur),zmem(1,index),
     +                               zmem(1,imem1),nai,zmem(1,imem3),
     +                               zmem(1,8),zmem(1,imem4),nci,
     +                               zmem(1,imem5),ic1)
                           go to 120
                        else
                           if (naind(index).eq.1) go to 90
                           nai = ixyz(index,natred,ktyp2)
                           if (nai.eq.0) go to 100
                           ncind(2) = ixyz(2,natred,ktyp2) + 1
                           ncind(3) = ixyz(3,natred,ktyp2) + 1
                           ncind(index) = ncind(index) - 1
                           imem3 = jmem3 + iarray(ncind(2),ncind(3))
                           imem4 = jmem4 + iarray(ncind(2),ncind(3))
                           naind(index) = naind(index) - 1
                           nci = naind(index)
                           imem5 = (iarray(naind(2),naind(3))-1)
     +                             *kvl1 + natred + mem5
                           naind(index) = naind(index) + 1
                        end if
                        call rex8(zmem(1,memcur),zmem(1,index),
     +                            zmem(1,imem1),nai,zmem(1,imem3),
     +                            zmem(1,8),zmem(1,imem4),nci,
     +                            zmem(1,imem5),ic1)
 120                 continue
 130              continue
 140           continue
c
c
c***** now we have all the integrals. use the horizontal
c***** recursion relation (if necessary) to find the desired integral
c
c              ihzc = 0
c
               nc12 = 0
               do 150 nc11 = icor1 , icor2
                  nc12 = nc12 + 1
                  zmem(nc12,9) = xijt(nc11)
                  xijt(nc11) = 0.0d0
 150           continue
c
c***** (ds|ds) integral
c
               if (jx.eq.1) then
                  mem1 = ihex(1)
                  do 170 ij = 1 , 6
                     imem1 = mem1 + (ij-1)*6 + ij
                     nc12 = 0
                     weight = 1.0d0
                     if (ij.eq.2 .or. ij.eq.3 .or. ij.eq.5) weight = sq3
                     do 160 nc11 = icor1 , icor2
                        nc12 = nc12 + 1
                        temp = weight*dsqrt(zmem(nc12,imem1))
     +                         *zmem(nc12,9)
                        xijt(nc11) = max(xijt(nc11),temp)
 160                 continue
 170              continue
               end if
               if (jx.eq.2) then
c
c***** (dp|dp) integral
c
                  mem1 = ihex(1)
                  mem2 = ihex(2)
                  mem3 = ihex(3)
                  do 200 ixd = 1 , 6
                     naind(2) = ixyz(2,ixd,2) + 1
                     naind(3) = ixyz(3,ixd,2) + 1
                     weight = 1.0d0
                     if (ixd.eq.2 .or. ixd.eq.3 .or. ixd.eq.5)
     +                   weight = sq3
                     do 190 jxp = 1 , 3
                        naind(jxp) = naind(jxp) + 1
                        ijk = iarray(naind(2),naind(3))
                        imem1 = mem1 + (ijk-1)*10 + ijk
                        imem2 = mem2 + (ijk-1)*6 + ixd
                        imem3 = mem3 + (ixd-1)*6 + ixd
                        naind(jxp) = naind(jxp) - 1
                        nc12 = 0
                        abind = 3 + jxp
                        do 180 nc11 = icor1 , icor2
                           nc12 = nc12 + 1
                           temp = zmem(nc12,imem1)
     +                            + 2.0d0*zmem(nc12,abind)
     +                            *zmem(nc12,imem2) + zmem(nc12,abind)
     +                            *zmem(nc12,abind)*zmem(nc12,imem3)
                           temp = dsqrt(temp)*weight*zmem(nc12,9)
                           xijt(nc11) = max(xijt(nc11),temp)
 180                    continue
 190                 continue
 200              continue
               end if
               if (jx.eq.3) then
c
c***** (dd|dd) integral
c
                  mem1 = ihex(1)
                  mem2 = ihex(2)
                  mem3 = ihex(3)
                  mem4 = ihex(4)
                  mem5 = ihex(5)
                  mem6 = ihex(6)
                  do 240 ij = 1 , 6
                     naind(2) = ixyz(2,ij,2) + 1
                     naind(3) = ixyz(3,ij,2) + 1
                     weight = 1.0d0
                     if (ij.eq.2 .or. ij.eq.3 .or. ij.eq.5) weight = sq3
                     kl = 0
                     do 230 kxp = 1 , 3
                        naind(kxp) = naind(kxp) + 1
                        ijk = iarray(naind(2),naind(3))
                        naind(kxp) = naind(kxp) - 1
                        do 220 lxp = kxp , 3
c
                           kl = kl + 1
                           weigh2 = 1.0d0
                           if (kl.eq.2 .or. kl.eq.3 .or. kl.eq.5)
     +                         weigh2 = sq3
                           naind(lxp) = naind(lxp) + 1
                           ijl = iarray(naind(2),naind(3))
                           naind(kxp) = naind(kxp) + 1
                           ijkl = iarray(naind(2),naind(3))
                           naind(kxp) = naind(kxp) - 1
                           naind(lxp) = naind(lxp) - 1
c
                           imem1 = mem1 + (ijkl-1)*15 + ijkl
                           imem2a = mem2 + (ijkl-1)*10 + ijk
                           imem2b = mem2 + (ijkl-1)*10 + ijl
                           imem3 = mem3 + (ijkl-1)*6 + ij
                           imem4a = mem4 + (ijk-1)*10 + ijk
                           imem4b = mem4 + (ijl-1)*10 + ijl
                           imem4c = mem4 + (ijk-1)*10 + ijl
                           imem5a = mem5 + (ijk-1)*6 + ij
                           imem5b = mem5 + (ijl-1)*6 + ij
                           imem6 = mem6 + (ij-1)*6 + ij
c
                           nc12 = 0
                           abindk = 3 + kxp
                           abindl = 3 + lxp
                           do 210 nc11 = icor1 , icor2
                              nc12 = nc12 + 1
                              abksq = zmem(nc12,abindk)
     +                                *zmem(nc12,abindk)
                              ablsq = zmem(nc12,abindl)
     +                                *zmem(nc12,abindl)
                              abk2 = zmem(nc12,abindk)*2.0d0
                              abkl2 = abk2*zmem(nc12,abindl)
                              abl2 = zmem(nc12,abindl)*2.0d0
c**** imem1  (gijkl | gijkl)     imem2a (gijkl | fijk)
c**** imem2b (gijkl | fijl )     imem3  (gijkl | dij )
c**** imem4a (fijk  | fijk )     imem4b (fijl  | fijl)
c**** imem4c (fijk  | fijl )     imem5a (fijk  | dij )
c**** imem5b (fijl  | dij  )     imem6  (dij   | dij )
c********************************************************
                              temp = zmem(nc12,imem1)
     +                               + abl2*zmem(nc12,imem2a)
     +                               + abk2*zmem(nc12,imem2b)
     +                               + abkl2*zmem(nc12,imem3)
     +                               + ablsq*zmem(nc12,imem4a)
     +                               + abksq*zmem(nc12,imem4b)
     +                               + abkl2*zmem(nc12,imem4c)
     +                               + abk2*ablsq*zmem(nc12,imem5a)
     +                               + abl2*abksq*zmem(nc12,imem5b)
     +                               + abksq*ablsq*zmem(nc12,imem6)
                              temp = dsqrt(temp)
     +                               *weigh2*weight*zmem(nc12,9)
                              xijt(nc11) = max(xijt(nc11),temp)
 210                       continue
 220                    continue
 230                 continue
 240              continue
               end if
               icorig = icor2
               ic1 = 0
               icount = icorig
            end if
 250     continue
 260  continue
c      write(iwr,*) 'the xijt arrray'
c      do 400 i=1,nij
c      temp=xijt(i)**2/2.0d0
c      write(iwr,*) temp
c400   continue
      return
      end
      subroutine fbuild(pold,hx,fcon)
c
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      dimension hx(n,n),pold(n,n)
      dimension fcon(n,n)
c
c
c        form new fock matrix
c
      do 30 i = 1 , n
         do 20 j = 1 , n
            fcon(i,j) = hx(i,j) + fcon(i,j)
 20      continue
 30   continue
c
c          workout energy and convergence property
c
c
      energy = 0.0d0
      do 50 j = 1 , n
         do 40 i = 1 , n
            energy = energy + 0.5d0*pold(i,j)*(hx(i,j)+fcon(i,j))
 40      continue
 50   continue
      energy = energy + repen
      conv = dcmax
c
      return
      end
      subroutine fcontr(fcon,fstore,fock,n)
c
c******* routine forms full fock matrix.
      implicit REAL  (a-h,o-z)
      logical ofull
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3
      common/defunk/cobas(maxorb,3),nd
      common/gen/maxit,iter
INCLUDE(common/mapper)
      common/scra  /iso(mxshel,48),nt
      common/gen3/ofull(100)
      dimension fcon(n,n),fstore(n,n),fock(n*n)
c
c.....now form the lower triangle, incremental, skeletal
c.....                   fock matrix
c
      icount = 0
      do 30 i = 1 , n
         do 20 j = 1 , i
            icount = icount + 1
            fock(icount) = 0.5d0*(fcon(i,j)+fcon(j,i))
 20      continue
 30   continue
c
      if (nd.gt.0) then
c
c....if d functions present weight the contributions from dxy
c.... dxz and dyz by a factor of root 3
c
         sq3 = dsqrt(3.0d0)
         icount = 0
         ic = 0
         do 50 i = 1 , n
            weigh1 = 1.0d0
            ix = i - ns - np3
            if (ix.gt.0) then
               ic = ic + 1
               if (ic.eq.2 .or. ic.eq.3 .or. ic.eq.5) weigh1 = sq3
               if (ic.eq.6) ic = 0
            end if
            jc = 0
            do 40 j = 1 , i
               icount = icount + 1
               weigh2 = 1.0d0
               jx = j - ns - np3
               if (jx.gt.0) then
                  jc = jc + 1
                  if (jc.eq.2 .or. jc.eq.3 .or. jc.eq.5) weigh2 = sq3
                  if (jc.eq.6) jc = 0
               end if
               fock(icount) = fock(icount)*weigh1*weigh2
 40         continue
 50      continue
      end if
c
c**** call the symetrising routine for skeletal fock matrix
c
      if (nt.gt.1) call symhd(fock,fcon,iky)
c
c***** form full square fock matrix from the incremented
c***** lower triangle fock matrix and old fock matrix.
c
      icount = 0
      if (iter.eq.80 .or. iter.eq.90.or.ofull(iter) ) then
c***** recalc of fock matrix ocuured, so fock contains
c***** true matrix
         do 70 i = 1 , n
            do 60 j = 1 , i
               icount = icount + 1
               fstore(i,j) = fock(icount)
               fstore(j,i) = fock(icount)
 60         continue
 70      continue
      else
         do 90 i = 1 , n
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 80 j = 1 , i
               icount = icount + 1
               fstore(i,j) = fstore(i,j) + fock(icount)
               fstore(j,i) = fstore(i,j)
 80         continue
 90      continue
      end if
c
      call dcopy(n*n,fstore,1,fcon,1)
c
      return
      end
      subroutine fgrid
c
c     ****prepare integral grid for fn(t)********
c
c     ***** description of parameters in routine : - ******
c     ***** nn is such that all values of fn(t) between ***
c     ***** t=0.0 and xmax at intervals of xint for *******
c     ***** n = 0 to (nm-1) are included in grid. *********
c
c
      implicit REAL  (a-h,o-z)
_IF1(cfu)      parameter(xint=0.05d0)
_IFN1(cfu)      parameter(xint=0.05d0)
      parameter(nn=1221)
      parameter(nm=15)
      parameter(xmax=61)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      common/ilist/fn(nn,nm),fex(nm),a(nm,4)
      dimension f(60),madd(20),dji(60)
      data madd/17,17,18,18,19,19,20,20,21,21,22,22,
     &   23,23,24,24,25,25,26,26/
      do 19 i=1,15
      do 19 j=1,4
 19   a(i,j) = 0.0d0
c
      a( 1,1) =  9.1d0
      a( 2,1) = 12.1d0
      a( 3,1) = 14.5d0
      a( 4,1) = 16.6d0
      a( 5,1) = 18.6d0
      a( 6,1) = 20.5d0
      a( 7,1) = 22.3d0
      a( 8,1) = 24.1d0
      a( 9,1) = 25.8d0
      a(10,1) = 27.4d0
      a( 1,2) = 13.7d0
      a( 2,2) = 17.1d0
      a( 3,2) = 19.7d0
      a( 4,2) = 22.1d0
      a( 5,2) = 24.3d0
      a( 6,2) = 26.4d0
      a( 7,2) = 28.4d0
      a( 8,2) = 30.3d0
      a( 9,2) = 32.2d0
      a(10,2) = 34.0d0
      a( 1,3) = 18.7d0
      a( 2,3) = 22.5d0
      a( 3,3) = 25.4d0
      a( 4,3) = 28.0d0
      a( 5,3) = 30.4d0
      a( 6,3) = 32.6d0
      a( 7,3) = 34.8d0
      a( 8,3) = 36.9d0
      a( 9,3) = 38.9d0
      a(10,3) = 40.8d0
      a( 1,4) = 24.2d0
      a( 2,4) = 28.1d0
      a( 3,4) = 31.1d0
      a( 4,4) = 33.8d0
      a( 5,4) = 36.3d0
      a( 6,4) = 38.7d0
      a( 7,4) = 41.0d0
      a( 8,4) = 43.2d0
      a( 9,4) = 45.3d0
      a(10,4) = 47.4d0
c
      do 20 i = 1 , 60
         dji(i) = 1.0d0/(i+i-1)
 20   continue
      n = anint(xmax/xint) + 1
      do 40 i = 1 , n
         x = xint*(i-1)
c*** drcfun generates the incomplete gamma function of x for
c**** m=0,npt*2
         npt = (nm+3)/2
         call drcfun(x,f,npt,dji,madd)
         do 30 j = 1 , nm
            fn(i,j) = f(j)
 30      continue
 40   continue
      temp = dsqrt(pi)*0.5d0
      fex(1) = temp
      do 50 i = 1 , nm - 1
         fex(i+1) = fex(i)*(2.0d0*dfloat(i)-1.0d0)*0.5d0
 50   continue
      return
      end
_EXTRACT(fockb,mips4)
      subroutine fockb(g,ivl,jvl,kvl,lvl,ic,jc,kmin,kmax
     &,lmin,lmax,l1,l2,ext2
     &,f,d,temp,memreq,ipos,jpos,kpos,lpos,i1,j1,k1)
c
c
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
c
      dimension g(icsize,memreq),f(nb,nb),d(nb,nb)
     &,temp(icsize*6)
      logical ext2
c
c***** routine to form fock matrix for a given i and j
c***** shell pairing. ivl/jvl are the number of basis functions
c***** in the corresponding shell
c
c***** takes about 25% of direct scf calculation- room for improvement
c
      mem = 0
      morig = 0
      do 90 i = 1 , ivl
         id = (ic-i1)*ivl + i + ipos
         do 80 j = 1 , jvl
            fij = 0.0d0
            jd = (jc-j1)*jvl + j + jpos
            dij = 4.0d0*d(jd,id)
            do 70 k = 1 , kvl
c
               morig = mem
               mem = mem + 1
               jcount = 0
               do 60 kc = kmin , kmax
                  lmin1 = l1
                  lmax1 = l2
                  if (ext2) lmax1 = kc
                  if (kc.eq.kmax) lmax1 = lmax
                  if (kc.eq.kmin) lmin1 = lmin
                  kd = (kc-k1)*kvl + k + kpos
                  djk = -1.0d0*d(kd,jd)
                  dik = -1.0d0*d(kd,id)
c
                  if (lx.eq.1) then
c*********** dealing with an s function **********
                     ld1 = lmin1
                     ld2 = lmax1
                     ld21 = ld2 - ld1 + 1
                     jstart = jcount + 1
                     do 20 ld = ld1 , ld2
                        jcount = jcount + 1
                        ge = g(jcount,mem)
                        f(ld,kd) = f(ld,kd) + ge*dij
                        f(ld,id) = f(ld,id) + ge*djk
                        f(ld,jd) = f(ld,jd) + ge*dik
 20                  continue
      fij = fij + ddot(ld21,d(ld1,kd),1,g(jstart,mem),1)
      f(id,kd) = f(id,kd) - ddot(ld21,d(ld1,jd),1,g(jstart,mem),1)
      f(jd,kd) = f(jd,kd) - ddot(ld21,d(ld1,id),1,g(jstart,mem),1)
                  else
c*****  dealing with a p/d/f function
c*****  the first loop rearranges all the integrals for a given
c*****  i,j,k  basis  function triplet so that they run consecutively
c**** in order that basis functions are in fock and density matrices
                     mem = morig
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                     do 40 l = 1 , lvl
                        jc0 = jcount
                        lcount = l - lvl
                        mem = mem + 1
                        do 30 lc = lmin1 , lmax1
                           lcount = lcount + lvl
                           jc0 = jc0 + 1
                           temp(lcount) = g(jc0,mem)
 30                     continue
 40                  continue
                     jcount = jc0
                     jc0 = 0
                     ld1 = (lmin1-l1)*lvl + 1 + lpos
                     ld2 = (lmax1-l1)*lvl + lvl + lpos
                     ld21 = ld2 - ld1 + 1
                     do 50 ld = ld1 , ld2
                        jc0 = jc0 + 1
                        ge = temp(jc0)
                        f(ld,kd) = f(ld,kd) + ge*dij
                        f(ld,id) = f(ld,id) + ge*djk
                        f(ld,jd) = f(ld,jd) + ge*dik
 50                  continue
           fij = fij + ddot(ld21,d(ld1,kd),1,temp(1),1)
           f(id,kd) = f(id,kd) - ddot(ld21,d(ld1,jd),1,temp(1),1)
           f(jd,kd) = f(jd,kd) - ddot(ld21,d(ld1,id),1,temp(1),1)
                  end if
 60            continue
 70         continue
            f(id,jd) = f(id,jd) + 4.0d0*fij
 80      continue
 90   continue
c
      return
      end
_ENDEXTRACT
      subroutine foutd(f,nb)
      implicit REAL  (a-h,o-z)
      dimension f(nb,nb)
      common/iofile/iread,iwr
      mmax = 8
      imax = 0
 20   imin = imax + 1
      imax = imax + mmax
      if (imax.gt.nb) imax = nb
      write (iwr,6010)
      write (iwr,6020) (i,i=imin,imax)
      write (iwr,6010)
      do 30 j = 1 , nb
         write (iwr,6030) j , (f(i,j),i=imin,imax)
 30   continue
      if (imax.lt.nb) go to 20
      return
 6010 format (/)
 6020 format (5x,8(6x,i3,6x))
 6030 format (i5,8f15.8)
      end
      subroutine fquik2(f,l2,l,t,m,ix)
c
c     *******routine for calculating fn(t)*******
c
c
c     t is the array containing the arguments
c     ix is the number of arguments
c     m is the number of incomplete gammas wanted
c     (i.e.  fn(t)  n = 0 . . . m-1)
c***** f(1) and l2(1) equivalenced
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      parameter (nn=1221)
      parameter (nm=15)
_IF1(cfu)      parameter (xint=0.05d0)
_IFN1(cfu)      parameter (xint=0.05d0)
      common/ilist/fn(nn,nm),fex(nm),aa(nm,4)
      logical glog1,glog2,glog5,glog6,glog9,gdbg,gmp2
     &,readmo,writmo,restar,acvary
      common/gen/maxitd,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog9,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary
      common/gen2/zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     + zipfal(21,2),zijkal(21,2),
     + prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
      dimension f(ix,m+4),t(ix),l(ix),l2(ix)
c
c** routine calculates the incomplete gamma function using a taylor
c** series for low values of t and an extrapolation to infinity for
c** large values. the tabulated values of fm(t') are evaluated by
c** drcfun which is accurate to 1.8d-12 (and relatively accurate to
c** 8.0d-11) for all values of t=0.0-60.0 and m=0-9.
c** iaccur=1 : 5 term taylor : accurate to 7.0d-12   (5.0d-11)
c**        2 : 4 terms       :             2.0d-9    (1.0d-8)
c**        3 : 3 terms       :             4.0d-7    (2.0d-6)
c**        4 : 2 terms       :             6.0d-5    (2.5d-4)
c** the accuracy of the extrapolation formulae is set one magnitude
c** tighter. the accuracy is varied according to closeness to
c** convergence and can be made constant with accuracy directive.
c
c*** addition of another term in the taylor expansion is not justified
c*** without changing drcfun. see obara-saika paper for more discussion
c*** of accuracies
c
      xint1 = 1.0d0/xint
      x12 = 0.5d0
      x13 = 1.0d0/3
      x14 = 0.25d0
      a = aa(m,iaccur)
c
c     *** use taylor series expansion from tabulated **
c     *** values of f, to find fm(t) for low values of t **
c
_IFN1(c)      call dlstlt(ix,t(1),1,a,k,l(1))
_IF1(c)      call whenflt(ix,t(1),1,a,l(1),k)
      if (k.gt.0) then
         call dgthr(k,t(1),f(1,m+2),l(1))
         do 20 i = 1 , k
            l2(i) = anint(f(i,m+2)*xint1)
            f(i,m+3) = dfloat(l2(i))*xint - f(i,m+2)
 20      continue
c
         call dgthr(k,fn(2,m+iaccur),f(1,m+4),l2(1))
         call dgthr(k,fn(2,m+iaccur-1),f(1,m+1),l2(1))
         go to (90,70,50,30) , iaccur
c
 30      do 40 i = 1 , k
            f(i,m+4) = f(i,m+1) + f(i,m+4)*f(i,m+3)*x14
 40      continue
c
         call dgthr(k,fn(2,m+2),f(1,m+1),l2(1))
c
 50      do 60 i = 1 , k
            f(i,m+4) = f(i,m+1) + f(i,m+4)*f(i,m+3)*x13
 60      continue
c
         call dgthr(k,fn(2,m+1),f(1,m+1),l2(1))
c
 70      do 80 i = 1 , k
            f(i,m+4) = f(i,m+1) + f(i,m+4)*f(i,m+3)*x12
 80      continue
c
         call dgthr(k,fn(2,m),f(1,m+1),l2(1))
c
 90      do 100 i = 1 , k
            f(i,m+4) = f(i,m+1) + f(i,m+4)*f(i,m+3)
 100     continue
c
c     ** descramble the f(*,m+4) array containing **
c     ** fm(t) and put the result in f(*,m) **
c
         call dsctr(k,f(1,m+4),l(1),f(1,m))
      end if
c
c
c     *** use extrapolation  ***
c     *** to infinity formula to find fm(t) for large t **
c
      if (k.lt.ix) then
_IFN1(c)         call dlstge(ix,t(1),1,a,kk,l(1))
_IF1(c)         call whenfge(ix,t(1),1,a,l(1),kk)
         call dgthr(kk,t(1),f(1,m+2),l(1))
c*** the next code is avoiding the clearer commented out code
ccccc      nmin1=m-1
ccccc      do 30 i=1,kk
ccccc30    f(i,m+4)=fex(m)/(sqrt(f(i,m+2))*f(i,m+2)**nmin1)
         go to (110,130,150,170,190,210,230,250,270,290) , m
 110     do 120 i = 1 , kk
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2)))
 120     continue
         go to 310
 130     do 140 i = 1 , kk
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*f(i,m+2))
 140     continue
         go to 310
 150     do 160 i = 1 , kk
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*f(i,m+2)*f(i,m+2))
 160     continue
         go to 310
 170     do 180 i = 1 , kk
           f(i,m+4)=fex(m)/(dsqrt(f(i,m+2))*f(i,m+2)*f(i,m+2)*f(i,m+2))
 180     continue
         go to 310
 190     do 200 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x2*x2)
 200     continue
         go to 310
 210     do 220 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x2*x2*f(i,m+2))
 220     continue
         go to 310
 230     do 240 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x2*x2*x2)
 240     continue
         go to 310
 250     do 260 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x2*x2*x2*f(i,m+2))
 260     continue
         go to 310
 270     do 280 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            x4 = x2*x2
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x4*x4)
 280     continue
         go to 310
 290     do 300 i = 1 , kk
            x2 = f(i,m+2)*f(i,m+2)
            x4 = x2*x2
            f(i,m+4) = fex(m)/(dsqrt(f(i,m+2))*x4*x4*f(i,m+2))
 300     continue
c
 310     call dsctr(kk,f(1,m+4),l(1),f(1,m))
         ziffew = ziffew + kk
      end if
      zifall = zifall + ix
c
c     ** use downward recursion to generate full set **
c
      if (m.gt.1) then
         do 320 i = 1 , ix
            f(i,m+1) = 0.5d0*dexp(-t(i))
 320     continue
         do 340 j = 1 , m - 1
            ni = m - j
            cx = 1.0d0/(dfloat(ni)-x12)
            do 330 i = 1 , ix
               f(i,ni) = (f(i,ni+1)*t(i)+f(i,m+1))*cx
 330        continue
 340     continue
      end if
c
      return
      end
      subroutine hrz2(ans,xint1,abi,xint2,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num)
c
      do 20 i = 1 , num
         ans(i) = xint1(i) + abi*xint2(i)
 20   continue
c
      return
      end
      subroutine hrz2c(ans,xint1,cdi,xint2,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),cdi(num)
c
      do 20 i = 1 , num
         ans(i) = xint1(i) + cdi(i)*xint2(i)
 20   continue
c
      return
      end
      subroutine hrz4(ans,xint1,abi,xint2,abj,xint3,xint4,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num)
c
      abij = abi*abj
      do 20 i = 1 , num
         ans(i) = xint1(i) + abi*xint2(i) + abj*xint3(i) + abij*xint4(i)
 20   continue
c
      return
      end
      subroutine hrz4c(ans,xint1,cdi,xint2,cdj,xint3,xint4,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),cdi(num),cdj(num)
c
      do 20 i = 1 , num
         ans(i) = xint1(i) + cdi(i)*xint2(i) + cdj(i)
     +            *(xint3(i)+cdi(i)*xint4(i))
 20   continue
c
      return
      end
      subroutine hrz8(ans,xint1,abi,xint2,abj,xint3,abk
     &,xint4,xint5,xint6,xint7,xint8,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),xint5(num),xint6(num),xint7(num),xint8(num)
c
      abij = abi*abj
      abik = abi*abk
      abjk = abj*abk
      abijk = abij*abk
      do 20 i = 1 , num
         ans(i) = xint1(i) + abi*xint2(i) + abj*xint3(i) + abk*xint4(i)
     +            + abij*xint5(i) + abik*xint6(i) + abjk*xint7(i)
     +            + abijk*xint8(i)
 20   continue
c
      return
      end
      subroutine hrz8c(ans,xint1,cdi,xint2,cdj,xint3,cdk
     &,xint4,xint5,xint6,xint7,xint8,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),xint5(num),xint6(num),xint7(num),xint8(num)
     &,cdi(num),cdj(num),cdk(num)
c
      do 20 i = 1 , num
         ans(i) = xint1(i) + cdk(i)*xint4(i) + cdi(i)
     +            *(xint2(i)+cdk(i)*xint6(i)+cdj(i)
     +            *(xint5(i)+cdk(i)*xint8(i))) + cdj(i)
     +            *(xint3(i)+cdk(i)*xint7(i))
 20   continue
c
      return
      end
      subroutine hrzd0(grad,gamma,xint1,nci,xint2,num,term2)
c
c**** routine calculates the gradient from the quantities related to
c*** the derivcative integrals
c
      implicit REAL  (a-h,o-z)
      logical term2
      dimension grad(num),xint1(num),xint2(num),gamma(num)
c
      if (term2) then
         scalar = dfloat(nci)
         do 20 i = 1 , num
            grad(i) = grad(i) + gamma(i)*(xint1(i)-scalar*xint2(i))
 20      continue
      else
         do 30 i = 1 , num
            grad(i) = grad(i) + gamma(i)*xint1(i)
 30      continue
      end if
c
      return
      end
      subroutine hzcdr(i1,j1,k1,l1,i2,j2,k2,l2,mem,jcount
     &,iw1,iw2,ijxz1,klxz1)
c
c
c**** routine decides on the memory locations and requirements of the
c**** derivative integrals during the hrr phase of the integral eval.
      implicit REAL  (a-h,o-z)
      dimension iw1(ijxz1,klxz1),iw2(ijxz1,klxz1)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,20)
c
c**** the (i1 j1 | k1 l1) term is the larger of the derivatives
c
      if (i1.ne.i2) ijkl2 = i2
      if (j1.ne.j2) write (iwr,*) 'something wrong????'
      if (k1.ne.k2) ijkl2 = k2
      if (l1.ne.l2) ijkl2 = l2
      ij1 = i1 + j1 - 1
      kl1 = k1 + l1 - 1
      ij2 = i2 + j2 - 1
      kl2 = k2 + l2 - 1
c
      jca = jcount
      if (j1.gt.1) then
c**** hrr is applied to the ij position for all components
c**** of the derivative
c**** terms giving (i1 j1 |k1 l1)
         nadd1 = nvl(i1)*nvl(j1)
         do 30 kl = k1 , kl1
            jcount = jcount + 1
            ihz(1,jcount) = mem
            mem = mem + nadd1*nvl(kl)
            do 20 ij = i1 , ij1
               indi = (ij1-ij) + 2
               ihz(indi,jcount) = iw1(ij,kl)
 20         continue
 30      continue
         jcb = jcount
c**** terms giving (i2 j2 | k2 l2)
         if (ijkl2.gt.0) then
            nadd1 = nvl(i2)*nvl(j2)
            do 50 kl = k2 , kl2
               jcount = jcount + 1
               ihz(1,jcount) = mem
               mem = mem + nadd1*nvl(kl)
               do 40 ij = i2 , ij2
                  indi = (ij2-ij) + 2
                  ihz(indi,jcount) = iw2(ij,kl)
 40            continue
 50         continue
         end if
      else
c***** hrr need not be applied to ij position
         do 60 kl = k1 , kl1
            jcount = jcount + 1
            ihz(1,jcount) = iw1(ij1,kl)
 60      continue
         jcb = jcount
         if (ijkl2.gt.0) then
            do 70 kl = k2 , kl2
               jcount = jcount + 1
               ihz(1,jcount) = iw2(ij2,kl)
 70         continue
         end if
      end if
      jcc = jcount
c*** now apply hrr to kl position for (i-1 j | k l)
      if (ijkl2.gt.0) then
         if (l2.eq.1) then
c*** hrr need not be applied
            ihz(1,jcount+1) = ihz(1,jcount)
         else
            jc1 = jcount + 1
            ihz(1,jc1) = mem
            mem = mem + nvl(i2)*nvl(j2)*nvl(k2)*nvl(l2)
            do 80 jc = jcb + 1 , jcc
               indj = jcc - jc + 2
               ihz(indj,jc1) = ihz(1,jc)
 80         continue
         end if
         jcount = jcount + 1
      end if
c*** now apply hrr to kl position for (i+1 j | k l)
      if (l1.eq.1) then
c*** hrr need not be applied
         ihz(1,jcount+1) = ihz(1,jcb)
      else
         jc1 = jcount + 1
         ihz(1,jc1) = mem
         do 90 jc = jca + 1 , jcb
            indj = jcb - jc + 2
            ihz(indj,jc1) = ihz(1,jc)
 90      continue
         mem = mem + nvl(i1)*nvl(j1)*nvl(k1)*nvl(l1)
      end if
      jcount = jcount + 1
_IF1()c      do 200 jc=jca+1,jcount
_IF1()c      write(iwr,210) ihz(1,jc),ihz(2,jc),ihz(3,jc),ihz(4,jc),ihz(5,jc)
_IF1()c210   format(1x,5i10)
_IF1()c200   continue
      return
      end
      subroutine hzcdrd(i1,j1,i2,j2,mem,jcount
     &,iw1,iw2,ijxz1)
c
c**** routine decides on use of hrr for derivative 1-e integrals
      implicit REAL  (a-h,o-z)
      dimension iw1(ijxz1),iw2(ijxz1)
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,20)
c
      if (i1.ne.i2) ijkl2 = i2
      if (j1.ne.j2) ijkl2 = j2
      ij1 = i1 + j1 - 1
      ij2 = i2 + j2 - 1
c
_IF1()      jca = jcount
      jcount = jcount + 1
c*** now apply hrr to kl position for (i-1 j | k l)
      if (ijkl2.gt.0) then
         if (j2.eq.1) then
c*** hrr need not be applied
            ihz(1,jcount) = iw2(ij2)
         else
            ihz(1,jcount) = mem
            mem = mem + nvl(i2)*nvl(j2)
            do 20 ij = i2 , ij2
               indj = ij2 - ij + 2
               ihz(indj,jcount) = iw2(ij)
 20         continue
         end if
      end if
c*** now apply hrr to kl position for (i+1 j | k l)
      jcount = jcount + 1
      if (j1.eq.1) then
c*** hrr need not be applied
         ihz(1,jcount) = iw1(ij1)
      else
         ihz(1,jcount) = mem
         mem = mem + nvl(i1)*nvl(j1)
         do 30 ij = i1 , ij1
            indj = ij1 - ij + 2
            ihz(indj,jcount) = iw1(ij)
 30      continue
      end if
_IF1()      do 40 jc = jca + 1 , jcount
_IF1()c      write(iwr,210) ihz(1,jc),ihz(2,jc),ihz(3,jc),ihz(4,jc),ihz(5,jc)
_IF1()c210   format(1x,5i10)
_IF1() 40   continue
      return
      end
      subroutine hzdrv(i1,j1,k1,l1,i2,j2,k2,l2,ihzc,rab
     &,zmemc,ab,jc1)
c
c**** rotuine coordinates formation of  derivative integrals using the h
c****  and works out their contribution to the gradient.
c
      implicit REAL  (a-h,o-z)
      dimension zmemc(icsize,memcto),ab(3),naind(3),ncind(3)
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,iunx,junx,kunx,lunx,munx,ijunx,klunx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mtgrad
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      logical ider,kder,lder,term2
c
      j1typ1 = j1 - 1
      i1typ1 = i1 - 1
      j2typ1 = j2 - 1
      i2typ1 = i2 - 1
      l1typ1 = l1 - 1
      k1typ1 = k1 - 1
      l2typ1 = l2 - 1
      k2typ1 = k2 - 1
      ivl1 = nvl(i1)
      jvl1 = nvl(j1)
      ivl2 = nvl(i2)
      jvl2 = nvl(j2)
      kvl1 = nvl(k1)
      lvl1 = nvl(l1)
      kvl2 = nvl(k2)
      lvl2 = nvl(l2)
      ider = .false.
      kder = .false.
      lder = .false.
      if (i1.ne.i2) then
         ijkl2 = i2
         ider = .true.
      end if
      if (k1.ne.k2) then
         ijkl2 = k2
         kder = .true.
c      nadd1=lvl1*nvl(k1typ1)
c      if (ijkl2.gt.0) nadd2=lvl1*nvl(k2)
      end if
      if (l1.ne.l2) then
         ijkl2 = l2
         lder = .true.
c      nadd1=kvl1*nvl(ltyp1)
c      if (ijkl2.gt.0) nadd2=kvl1*nvl(l2)
      end if
      if (ijkl2.eq.0) term2 = .false.
c     idouix = 0
c     ij1 = i1 + j1 - 1
      kl1 = k1 + l1 - 1
c     ij2 = i2 + j2 - 1
      kl2 = k2 + l2 - 1
c
      if (j1.gt.1) then
c***** hrr must be applied on the ij position for i1 j1 k1 l1
         do 80 kl = k1 , kl1
            ihzc = ihzc + 1
            mem1 = ihz(1,ihzc)
            mem2 = ihz(2,ihzc)
            mem3 = ihz(3,ihzc)
            mem4 = ihz(4,ihzc)
            mem5 = ihz(5,ihzc)
            kvl = nvl(kl)
            ijc = 0
            do 70 i = 1 , ivl1
               if (i1typ1.eq.0) then
                  naind(2) = 1
                  naind(3) = 1
               else
                  naind(2) = ixyz(2,i,i1typ1) + 1
                  naind(3) = ixyz(3,i,i1typ1) + 1
               end if
               do 60 j = 1 , jvl1
                  ncind(2) = naind(2) + ixyz(2,j,j1typ1)
                  ncind(3) = naind(3) + ixyz(3,j,j1typ1)
                  ij = iarray(ncind(2),ncind(3))
                  imem1 = ijc*kvl + mem1
                  imem2 = (ij-1)*kvl + mem2
                  ijc = ijc + 1
                  if (rab.lt.1.0d-8) then
c********* when (ai-bi) is zero (a+b,0|c+d,0) = (a,b|c+d,0)
                     do 20 k = 1 , kvl
                call dcopy(jc1,zmemc(1,imem2+k),1,zmemc(1,imem1+k),1)
 20                  continue
                     go to 60
                  end if
                  index1 = ispdf1(j,1)
                  ncind(index1) = ncind(index1) - 1
                  imem3a = (iarray(ncind(2),ncind(3))-1)*kvl + mem3
                  if (j1.eq.2) then
*****  a p function is being added to the i/a position.
                     do 30 k = 1 , kvl
                        call hrz2(zmemc(1,imem1+k),zmemc(1,imem2+k),
     +                            ab(index1),zmemc(1,imem3a+k),jc1)
 30                  continue
                     go to 60
                  end if
                  ired1 = ispdf1(j,2)
                  index2 = ispdf1(ired1,1)
                  ncind(index2) = ncind(index2) - 1
                  imem4a = (iarray(ncind(2),ncind(3))-1)*kvl + mem4
                  ncind(index1) = ncind(index1) + 1
                  imem3b = (iarray(ncind(2),ncind(3))-1)*kvl + mem3
                  if (j1.eq.3) then
*****  a d function is being added to the i/a position
                     do 40 k = 1 , kvl
                        call hrz4(zmemc(1,imem1+k),zmemc(1,imem2+k),
     +                            ab(index1),zmemc(1,imem3a+k),
     +                            ab(index2),zmemc(1,imem3b+k),
     +                            zmemc(1,imem4a+k),jc1)
 40                  continue
                     go to 60
                  end if
                  ired2 = ispdf1(ired1,2)
                  index3 = ispdf1(ired2,1)
                  ncind(index3) = ncind(index3) - 1
                  imem4b = (iarray(ncind(2),ncind(3))-1)*kvl + mem4
                  ncind(index2) = ncind(index2) + 1
                  imem3c = (iarray(ncind(2),ncind(3))-1)*kvl + mem3
                  ncind(index1) = ncind(index1) - 1
                  imem4c = (iarray(ncind(2),ncind(3))-1)*kvl + mem4
                  ncind(index2) = ncind(index2) - 1
                  imem5 = (iarray(ncind(2),ncind(3))-1)*kvl + mem5
*****  an f function is being added to the i/a position
                  do 50 k = 1 , kvl
                     call hrz8(zmemc(1,imem1+k),zmemc(1,imem2+k),
     +                         ab(index1),zmemc(1,imem3a+k),ab(index2),
     +                         zmemc(1,imem3b+k),ab(index3),
     +                         zmemc(1,imem3c+k),zmemc(1,imem4a+k),
     +                         zmemc(1,imem4c+k),zmemc(1,imem4b+k),
     +                         zmemc(1,imem5+k),jc1)
 50               continue
 60            continue
 70         continue
 80      continue
         if (ijkl2.gt.0) then
c***** hrr must be applied on the ij position for i1 j1 k1 l1
            do 140 kl = k2 , kl2
               ihzc = ihzc + 1
               mem1 = ihz(1,ihzc)
               mem2 = ihz(2,ihzc)
               mem3 = ihz(3,ihzc)
               mem4 = ihz(4,ihzc)
               mem5 = ihz(5,ihzc)
               kvl = nvl(kl)
               ijc = 0
               do 130 i = 1 , ivl2
                  if (i2typ1.eq.0) then
                     naind(2) = 1
                     naind(3) = 1
                  else
                     naind(2) = ixyz(2,i,i2typ1) + 1
                     naind(3) = ixyz(3,i,i2typ1) + 1
                  end if
                  do 120 j = 1 , jvl1
                     ncind(2) = naind(2) + ixyz(2,j,j2typ1)
                     ncind(3) = naind(3) + ixyz(3,j,j2typ1)
                     ij = iarray(ncind(2),ncind(3))
                     imem1 = ijc*kvl + mem1
                     imem2 = (ij-1)*kvl + mem2
                     ijc = ijc + 1
                     if (rab.lt.1.0d-8) then
c********* when (ai-bi) is zero (a+b,0|c+d,0) = (a,b|c+d,0)
                        do 90 k = 1 , kvl
                 call dcopy(jc1,zmemc(1,imem2+k),1,zmemc(1,imem1+k),1)
 90                     continue
                        go to 120
                     end if
                     index1 = ispdf1(j,1)
                     ncind(index1) = ncind(index1) - 1
                     imem3a = (iarray(ncind(2),ncind(3))-1)*kvl + mem3
                     if (j2.eq.2) then
*****  a p function is being added to the i/a position.
                        do 100 k = 1 , kvl
                           call hrz2(zmemc(1,imem1+k),zmemc(1,imem2+k),
     +                               ab(index1),zmemc(1,imem3a+k),jc1)
 100                    continue
                        go to 120
                     end if
                     ired1 = ispdf1(j,2)
                     index2 = ispdf1(ired1,1)
                     ncind(index2) = ncind(index2) - 1
                     imem4a = (iarray(ncind(2),ncind(3))-1)*kvl + mem4
                     ncind(index1) = ncind(index1) + 1
                     imem3b = (iarray(ncind(2),ncind(3))-1)*kvl + mem3
*****  a d function is being added to the i/a position
                     do 110 k = 1 , kvl
                        call hrz4(zmemc(1,imem1+k),zmemc(1,imem2+k),
     +                            ab(index1),zmemc(1,imem3a+k),
     +                            ab(index2),zmemc(1,imem3b+k),
     +                            zmemc(1,imem4a+k),jc1)
 110                 continue
 120              continue
 130           continue
 140        continue
         end if
c***** at this point the hrr on a and b was necessary and
c***** has been applied.
      else
c***** hrr not needed on a and b position
         ihzc = ihzc + kl1 - k1 + 1
         if (ijkl2.gt.0) ihzc = ihzc + kl2 - k2 + 1
      end if
c
c***** now the hrr is applied to the c and d position
c
      if (ijkl2.gt.0) then
c***** (i2 j2 | k2 l2) exists
c
         ihzc = ihzc + 1
         if (l2.gt.1) then
c***** hrr is needed for c and d in the second derivative term
            klc = 0
            mem1 = ihz(1,ihzc)
            mem2 = ihz(2,ihzc)
            mem3 = ihz(3,ihzc)
            mem4 = ihz(4,ihzc)
            mem5 = ihz(5,ihzc)
            klo1 = kvl2*lvl2
            klo2 = nvl(kl2)
            klo3 = nvl(kl2-1)
            klo4 = nvl(kl2-2)
            klo5 = nvl(kl2-3)
            do 200 k = 1 , kvl2
               if (k2typ1.eq.0) then
                  naind(2) = 1
                  naind(3) = 1
               else
                  naind(2) = ixyz(2,k,k2typ1) + 1
                  naind(3) = ixyz(3,k,k2typ1) + 1
               end if
               do 190 l = 1 , lvl2
                  ncind(2) = ixyz(2,l,l2typ1) + naind(2)
                  ncind(3) = ixyz(3,l,l2typ1) + naind(3)
                  klc = klc + 1
                  imem1 = mem1 + klc
                  kl = iarray(ncind(2),ncind(3))
                  imem2 = mem2 + kl
                  index1 = ispdf1(l,1)
                  indcd1 = mtemp + index1
                  ncind(index1) = ncind(index1) - 1
                  imem3a = mem3 + iarray(ncind(2),ncind(3))
                  ijc = -1
                  if (l2.eq.2) then
*****  a p function is being added to the k/c position
                     do 160 i = 1 , ivl2
                        do 150 j = 1 , jvl2
                           ijc = ijc + 1
                           jmem1 = imem1 + ijc*klo1
                           jmem2 = imem2 + ijc*klo2
                           jmem3 = imem3a + ijc*klo3
                           call hrz2c(zmemc(1,jmem1),zmemc(1,jmem2),
     +                                zmemc(1,indcd1),zmemc(1,jmem3),
     +                                jc1)
 150                    continue
 160                 continue
                     go to 190
                  end if
                  ired1 = ispdf1(l,2)
                  index2 = ispdf1(ired1,1)
                  indcd2 = mtemp + index2
                  ncind(index2) = ncind(index2) - 1
                  imem4a = mem4 + iarray(ncind(2),ncind(3))
                  ncind(index1) = ncind(index1) + 1
                  imem3b = mem3 + iarray(ncind(2),ncind(3))
c**** a d function is being added to the k/c position
                  do 180 i = 1 , ivl2
                     do 170 j = 1 , jvl2
                        ijc = ijc + 1
                        itemp = ijc*klo3
                        jmem1 = imem1 + ijc*klo1
                        jmem2 = imem2 + ijc*klo2
                        jmem3a = imem3a + itemp
                        jmem3b = imem3b + itemp
                        jmem4a = imem4a + ijc*klo4
                        call hrz4c(zmemc(1,jmem1),zmemc(1,jmem2),
     +                             zmemc(1,indcd1),zmemc(1,jmem3a),
     +                             zmemc(1,indcd2),zmemc(1,jmem3b),
     +                             zmemc(1,jmem4a),jc1)
 170                 continue
 180              continue
 190           continue
 200        continue
         end if
      end if
c**** here we have the 2nd term but the hrr must be applied to
c**** first term if necessary
c
c***** (i1 j1 | k1 l1) exists
c
c**** mem00 is the 2nd term in derivative integral
      mem00 = ihz(1,ihzc)
c**** mem0 is position of gamma(i,j,k,l)
      ihzc = ihzc + 1
c**** mem1 is the position of first term derivative integrals
      mem1 = ihz(1,ihzc)
      if (l1.gt.1) then
         klc = 0
         mem2 = ihz(2,ihzc)
         mem3 = ihz(3,ihzc)
         mem4 = ihz(4,ihzc)
         mem5 = ihz(5,ihzc)
         klo1 = kvl1*lvl1
         klo2 = nvl(kl1)
         klo3 = nvl(kl1-1)
         klo4 = nvl(kl1-2)
         klo5 = nvl(kl1-3)
         do 280 k = 1 , kvl1
            if (k1typ1.eq.0) then
               naind(2) = 1
               naind(3) = 1
            else
               naind(2) = ixyz(2,k,k1typ1) + 1
               naind(3) = ixyz(3,k,k1typ1) + 1
            end if
            do 270 l = 1 , lvl1
               ncind(2) = ixyz(2,l,l1typ1) + naind(2)
               ncind(3) = ixyz(3,l,l1typ1) + naind(3)
               klc = klc + 1
               imem1 = mem1 + klc
               kl = iarray(ncind(2),ncind(3))
               imem2 = mem2 + kl
               index1 = ispdf1(l,1)
               indcd1 = mtemp + index1
               ncind(index1) = ncind(index1) - 1
               imem3a = mem3 + iarray(ncind(2),ncind(3))
               ijc = -1
               if (l1.eq.2) then
*****  a p function is being added to the k/c position
                  do 220 i = 1 , ivl1
                     do 210 j = 1 , jvl1
                        ijc = ijc + 1
                        jmem1 = imem1 + ijc*klo1
                        jmem2 = imem2 + ijc*klo2
                        jmem3 = imem3a + ijc*klo3
                        call hrz2c(zmemc(1,jmem1),zmemc(1,jmem2),
     +                             zmemc(1,indcd1),zmemc(1,jmem3),jc1)
 210                 continue
 220              continue
                  go to 270
               end if
               ired1 = ispdf1(l,2)
               index2 = ispdf1(ired1,1)
               indcd2 = mtemp + index2
               ncind(index2) = ncind(index2) - 1
               imem4a = mem4 + iarray(ncind(2),ncind(3))
               ncind(index1) = ncind(index1) + 1
               imem3b = mem3 + iarray(ncind(2),ncind(3))
               if (l1.eq.3) then
c**** a d function is being added to the k/c position
                  do 240 i = 1 , ivl1
                     do 230 j = 1 , jvl1
                        ijc = ijc + 1
                        itemp = ijc*klo3
                        jmem1 = imem1 + ijc*klo1
                        jmem2 = imem2 + ijc*klo2
                        jmem3a = imem3a + itemp
                        jmem3b = imem3b + itemp
                        jmem4a = imem4a + ijc*klo4
                        call hrz4c(zmemc(1,jmem1),zmemc(1,jmem2),
     +                             zmemc(1,indcd1),zmemc(1,jmem3a),
     +                             zmemc(1,indcd2),zmemc(1,jmem3b),
     +                             zmemc(1,jmem4a),jc1)
 230                 continue
 240              continue
                  go to 270
               end if
               ired2 = ispdf1(ired1,2)
               index3 = ispdf1(ired2,1)
               indcd3 = mtemp + index3
               ncind(index3) = ncind(index3) - 1
               imem4b = iarray(ncind(2),ncind(3)) + mem4
               ncind(index2) = ncind(index2) + 1
               imem3c = iarray(ncind(2),ncind(3)) + mem3
               ncind(index1) = ncind(index1) - 1
               imem4c = iarray(ncind(2),ncind(3)) + mem4
               ncind(index2) = ncind(index2) - 1
               imem5 = iarray(ncind(2),ncind(3)) + mem5
               if (l1.eq.4) then
c**** a d function is being added to the k/c position
                  do 260 i = 1 , ivl1
                     do 250 j = 1 , jvl1
                        ijc = ijc + 1
                        itemp = ijc*klo3
                        itemp1 = ijc*klo4
                        jmem1 = imem1 + ijc*klo1
                        jmem2 = imem2 + ijc*klo2
                        jmem3a = imem3a + itemp
                        jmem3b = imem3b + itemp
                        jmem3c = imem3c + itemp
                        jmem4a = imem4a + itemp1
                        jmem4b = imem4b + itemp1
                        jmem4c = imem4c + itemp1
                        jmem5 = imem5 + ijc*klo5
                        call hrz8c(zmemc(1,jmem1),zmemc(1,jmem2),
     +                             zmemc(1,indcd1),zmemc(1,jmem3a),
     +                             zmemc(1,indcd2),zmemc(1,jmem3b),
     +                             zmemc(1,indcd3),zmemc(1,jmem3c),
     +                             zmemc(1,jmem4a),zmemc(1,jmem4c),
     +                             zmemc(1,jmem4b),zmemc(1,jmem5),jc1)
 250                 continue
 260              continue
               end if
 270        continue
 280     continue
      end if
ccccc
c***** hrr not needed and gradient contribution must be found
c***** hrr is needed for c and d in the first derivative term
      ivlx = ivl1
      kvlx = kvl1
      lvlx = lvl1
      if (ider) ivlx = nvl(i1-1)
      if (kder) kvlx = nvl(k1-1)
      if (lder) lvlx = nvl(l1-1)
      do 330 icart = 1 , 3
         mgrad = mtgrad + icart
         if (lder) mgrad = mgrad + 3
         if (kder) mgrad = mgrad + 6
         mem0 = 0
c        ijklc = 0
         ijc = -1
         ijkc = -1
         do 320 i = 1 , ivlx
            if (ider) then
               jklc = 0
               if (i2.eq.0) then
                  naind(2) = 1
                  naind(3) = 1
               else
                  naind(1) = ixyz(1,i,i2) + 1
                  naind(2) = ixyz(2,i,i2) + 1
                  naind(3) = ixyz(3,i,i2) + 1
               end if
               naind(icart) = naind(icart) + 1
               imem1 = mem1 + (iarray(naind(2),naind(3))-1)
     +                 *jvl1*kvl1*lvl1
               if (ijkl2.ne.0) then
                  naind(icart) = naind(icart) - 2
                  nai = naind(icart)
                  term2 = nai.ne.0
                  if (term2) imem00 = mem00 +
     +                                (iarray(naind(2),naind(3))-1)
     +                                *jvl1*kvl1*lvl1
               end if
            end if
            do 310 j = 1 , jvl1
               ijc = ijc + 1
               do 300 k = 1 , kvlx
                  ijkc = ijkc + 1
                  if (kder) then
                     if (k2.eq.0) then
                        naind(2) = 1
                        naind(3) = 1
                     else
                        naind(1) = ixyz(1,k,k2) + 1
                        naind(2) = ixyz(2,k,k2) + 1
                        naind(3) = ixyz(3,k,k2) + 1
                     end if
                     naind(icart) = naind(icart) + 1
                     imem1 = mem1 + ijc*lvl1*kvl1 +
     +                       (iarray(naind(2),naind(3))-1)*lvl1
                     if (ijkl2.ne.0) then
                        naind(icart) = naind(icart) - 2
                        nai = naind(icart)
                        term2 = nai.ne.0
                        if (term2) imem00 = mem00 + ijc*lvl1*kvl2 +
     +                      (iarray(naind(2),naind(3))-1)*lvl1
                     end if
c
                  end if
                  do 290 l = 1 , lvlx
                     jklc = jklc + 1
                     mem0 = mem0 + 1
                     if (lder) then
                        if (l2.eq.0) then
                           naind(2) = 1
                           naind(3) = 1
                        else
                           naind(1) = ixyz(1,l,l2) + 1
                           naind(2) = ixyz(2,l,l2) + 1
                           naind(3) = ixyz(3,l,l2) + 1
                        end if
                        naind(icart) = naind(icart) + 1
                        jmem1 = mem1 + ijkc*lvl1 +
     +                          iarray(naind(2),naind(3))
                        if (ijkl2.ne.0) then
                           naind(icart) = naind(icart) - 2
                           nai = naind(icart)
                           term2 = nai.ne.0
                           if (term2) jmem00 = mem00 + ijkc*lvl2 +
     +                         iarray(naind(2),naind(3))
                        end if
                        call hrzd0(zmemc(1,mgrad),zmemc(1,mem0),
     +                             zmemc(1,jmem1),nai,zmemc(1,jmem00),
     +                             jc1,term2)
                        go to 290
c
                     end if
c*****  a p function is being added to the k/c position
c
                     if (ider) then
                        jmem1 = imem1 + jklc
                        jmem00 = imem00 + jklc
                     else
                        jmem1 = imem1 + l
                        jmem00 = imem00 + l
                     end if
_IF1()c672        do 1001 jcount=1,jc1
_IF1()c           cont=zmemc(jcount,jmem1)
_IF1()c           if (term2) cont=cont-dfloat(nai)*zmemc(jcount,jmem00)
_IF1()c           cont=cont*zmemc(jcount,mem0)
_IF1()c           if (abs(cont).gt.1.0d-6) then
_IF1()c           write(iwr,*) 'contribution is ',cont,mgrad,ider,kder,lder
_IF1()c           write(iwr,*) 'gamma ',zmemc(jcount,mem0),i,j,k,l
_IF1()c           write(iwr,*) zmemc(jcount,jmem1),term2,zmemc(jcount,jmem00)
_IF1()c           endif
_IF1()c1001       continue
                     call hrzd0(zmemc(1,mgrad),zmemc(1,mem0),
     +                          zmemc(1,jmem1),nai,zmemc(1,jmem00),jc1,
     +                          term2)
 290              continue
 300           continue
 310        continue
 320     continue
 330  continue
c
      return
      end
      subroutine hzdrvd(i1,j1,i2,j2,ihzc
     &,zmemc,jc1,mgamma,mtemp,mtgrad)
c**** the fixed parameters in common/indexd/ are explained in the
c**** block data statement (in bdata.f)
c**** the memory locations and requirements were established in
c**** cidonb() and are stored in common/cider/ and cidriv
c
c**** rotuine uses hrr on the derivatives of the nuclear attraction
c**** integrals and finds their contribution to the gradient
      implicit REAL  (a-h,o-z)
      dimension zmemc(ipsize,mempto),naind(3),ncind(3)
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,iunx,junx,kunx,lunx,munx,ijunx,klunx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtempx,ipszt,icszt
     &,ipsize,icsize,memcon
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mgradx
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
      logical ider,jder,term2
c
      j1typ1 = j1 - 1
      i1typ1 = i1 - 1
      j2typ1 = j2 - 1
      i2typ1 = i2 - 1
      ivl1 = nvl(i1)
      jvl1 = nvl(j1)
      ivl2 = nvl(i2)
      jvl2 = nvl(j2)
      ider = .false.
      jder = .false.
      if (i1.ne.i2) then
         ijkl2 = i2
         ider = .true.
      end if
      if (j1.ne.j2) then
         ijkl2 = j2
         jder = .true.
      end if
      if (ijkl2.eq.0) term2 = .false.
c     ij1 = i1 + j1 - 1
c     ij2 = i2 + j2 - 1
c
c***** now the hrr is applied to the c and d position
c
      ihzc = ihzc + 1
      if (ijkl2.gt.0) then
c***** (i2 |a| j2) exists
         if (j2.gt.1) then
c***** hrr is needed for c and d in the second derivative term
            ijc = 0
            mem1 = ihz(1,ihzc)
            mem2 = ihz(2,ihzc)
            mem3 = ihz(3,ihzc)
            mem4 = ihz(4,ihzc)
            mem5 = ihz(5,ihzc)
            do 30 i = 1 , ivl2
               if (i2typ1.eq.0) then
                  naind(2) = 1
                  naind(3) = 1
               else
                  naind(2) = ixyz(2,i,i2typ1) + 1
                  naind(3) = ixyz(3,i,i2typ1) + 1
               end if
               do 20 j = 1 , jvl2
                  ncind(2) = ixyz(2,j,j2typ1) + naind(2)
                  ncind(3) = ixyz(3,j,j2typ1) + naind(3)
                  ijc = ijc + 1
                  imem1 = mem1 + ijc
                  ij = iarray(ncind(2),ncind(3))
                  imem2 = mem2 + ij
                  index1 = ispdf1(j,1)
                  indcd1 = mtemp + index1
                  ncind(index1) = ncind(index1) - 1
                  imem3a = mem3 + iarray(ncind(2),ncind(3))
                  if (j2.eq.2) then
                     call hrz2c(zmemc(1,imem1),zmemc(1,imem2),
     +                          zmemc(1,indcd1),zmemc(1,imem3a),jc1)
                     go to 20
                  end if
                  ired1 = ispdf1(j,2)
                  index2 = ispdf1(ired1,1)
                  indcd2 = mtemp + index2
                  ncind(index2) = ncind(index2) - 1
                  imem4a = mem4 + iarray(ncind(2),ncind(3))
                  ncind(index1) = ncind(index1) + 1
                  imem3b = mem3 + iarray(ncind(2),ncind(3))
c**** a d function is being added to the k/c position
                  call hrz4c(zmemc(1,imem1),zmemc(1,imem2),
     +                       zmemc(1,indcd1),zmemc(1,imem3a),
     +                       zmemc(1,indcd2),zmemc(1,imem3b),
     +                       zmemc(1,imem4a),jc1)
 20            continue
 30         continue
         end if
      end if
c**** here we have the 2nd term but the hrr must be applied to
c**** first term if necessary
c
c**** mem00 is the 2nd term in derivative integral
      mem00 = ihz(1,ihzc)
c**** mem0 is position of gamma(i,j,k,l)
      ihzc = ihzc + 1
c**** mem1 is the position of first term derivative integrals
      mem1 = ihz(1,ihzc)
      if (j1.gt.1) then
c***** (i1 |a| j1) exists
c
         ijc = 0
         mem2 = ihz(2,ihzc)
         mem3 = ihz(3,ihzc)
         mem4 = ihz(4,ihzc)
         mem5 = ihz(5,ihzc)
         do 50 i = 1 , ivl1
            if (i1typ1.eq.0) then
               naind(2) = 1
               naind(3) = 1
            else
               naind(2) = ixyz(2,i,i1typ1) + 1
               naind(3) = ixyz(3,i,i1typ1) + 1
            end if
            do 40 j = 1 , jvl1
               ncind(2) = ixyz(2,j,j1typ1) + naind(2)
               ncind(3) = ixyz(3,j,j1typ1) + naind(3)
               ijc = ijc + 1
               imem1 = mem1 + ijc
               ij = iarray(ncind(2),ncind(3))
               imem2 = mem2 + ij
               index1 = ispdf1(j,1)
               indcd1 = mtemp + index1
               ncind(index1) = ncind(index1) - 1
               imem3a = mem3 + iarray(ncind(2),ncind(3))
               if (j1.eq.2) then
*****  a p function is being added to the k/c position
                  call hrz2c(zmemc(1,imem1),zmemc(1,imem2),
     +                       zmemc(1,indcd1),zmemc(1,imem3a),jc1)
                  go to 40
               end if
               ired1 = ispdf1(j,2)
               index2 = ispdf1(ired1,1)
               indcd2 = mtemp + index2
               ncind(index2) = ncind(index2) - 1
               imem4a = mem4 + iarray(ncind(2),ncind(3))
               ncind(index1) = ncind(index1) + 1
               imem3b = mem3 + iarray(ncind(2),ncind(3))
               if (j1.eq.3) then
c**** a d function is being added to the k/c position
                  call hrz4c(zmemc(1,imem1),zmemc(1,imem2),
     +                       zmemc(1,indcd1),zmemc(1,imem3a),
     +                       zmemc(1,indcd2),zmemc(1,imem3b),
     +                       zmemc(1,imem4a),jc1)
                  go to 40
               end if
               ired2 = ispdf1(ired1,2)
               index3 = ispdf1(ired2,1)
               indcd3 = mtemp + index3
               ncind(index3) = ncind(index3) - 1
               imem4b = iarray(ncind(2),ncind(3)) + mem4
               ncind(index2) = ncind(index2) + 1
               imem3c = iarray(ncind(2),ncind(3)) + mem3
               ncind(index1) = ncind(index1) - 1
               imem4c = iarray(ncind(2),ncind(3)) + mem4
               ncind(index2) = ncind(index2) - 1
               imem5 = iarray(ncind(2),ncind(3)) + mem5
               if (j1.eq.4) then
c**** a d function is being added to the k/c position
                  call hrz8c(zmemc(1,imem1),zmemc(1,imem2),
     +                       zmemc(1,indcd1),zmemc(1,imem3a),
     +                       zmemc(1,indcd2),zmemc(1,imem3b),
     +                       zmemc(1,indcd3),zmemc(1,imem3c),
     +                       zmemc(1,imem4a),zmemc(1,imem4c),
     +                       zmemc(1,imem4b),zmemc(1,imem5),jc1)
               end if
 40         continue
 50      continue
      end if
ccccc
c***** hrr not needed and gradient contribution must be found
c***** hrr is needed for c and d in the first derivative term
      ivlx = ivl1
      jvlx = jvl1
      if (ider) ivlx = nvl(i1-1)
      if (jder) jvlx = nvl(j1-1)
      do 80 icart = 1 , 3
         mgrad = mtgrad + icart
         call vclr(zmemc(1,mgrad),1,jc1)
         mem0 = mgamma
         ijc = 0
         do 70 i = 1 , ivlx
            if (ider) then
               if (i2.eq.0) then
                  naind(2) = 1
                  naind(3) = 1
               else
                  naind(1) = ixyz(1,i,i2) + 1
                  naind(2) = ixyz(2,i,i2) + 1
                  naind(3) = ixyz(3,i,i2) + 1
               end if
               naind(icart) = naind(icart) + 1
               imem1 = mem1 + (iarray(naind(2),naind(3))-1)*jvl1
               if (ijkl2.ne.0) then
                  naind(icart) = naind(icart) - 2
                  nai = naind(icart)
                  term2 = nai.ne.0
                  if (term2) imem00 = mem00 +
     +                                (iarray(naind(2),naind(3))-1)*jvl1
               end if
            end if
            do 60 j = 1 , jvlx
               mem0 = mem0 + 1
               if (jder) then
                  if (j2.eq.0) then
                     naind(2) = 1
                     naind(3) = 1
                  else
                     naind(1) = ixyz(1,j,j2) + 1
                     naind(2) = ixyz(2,j,j2) + 1
                     naind(3) = ixyz(3,j,j2) + 1
                  end if
                  naind(icart) = naind(icart) + 1
                  jmem1 = mem1 + (i-1)*jvl1 + iarray(naind(2),naind(3))
                  if (ijkl2.ne.0) then
                     naind(icart) = naind(icart) - 2
                     nai = naind(icart)
                     term2 = nai.ne.0
                     if (term2) jmem00 = mem00 + (i-1)
     +                   *jvl2 + iarray(naind(2),naind(3))
                  end if
               end if
c*****  a p function is being added to the k/c position
c
               if (ider) then
                  jmem1 = imem1 + j
                  jmem00 = imem00 + j
               end if
_IF1()c           do 1001 jcount=1,jc1
_IF1()c           cont=zmemc(jcount,jmem1)
_IF1()c           if (term2) cont=cont-dfloat(nai)*zmemc(jcount,jmem00)
_IF1()c           cont=cont*zmemc(jcount,mem0)
_IF1()c           if (abs(cont).gt.1.0d-6) then
_IF1()c           write(iwr,*) 'contribution is ',cont,icart,ider,jder
_IF1()c           write(iwr,*) 'gamma ',zmemc(jcount,mem0),i,j
_IF1()c           write(iwr,*) zmemc(jcount,jmem1),term2,zmemc(jcount,jmem00)
_IF1()c           endif
_IF1()c1001       continue
               call hrzd0(zmemc(1,mgrad),zmemc(1,mem0),zmemc(1,jmem1),
     +                    nai,zmemc(1,jmem00),jc1,term2)
 60         continue
 70      continue
 80   continue
c
      return
      end
      subroutine icrs(iw,ix,kx,mx,icc,kc,mc)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix,kx,mx)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up ic
c
      ic = icc - 1
      mc1 = mc + 1
c**** first 2 terms of recursion relation
      if( ic.ge.1) then
       iw(ic,kc,mc) = -1
       iw(ic,kc,mc1) = -1
c**** fifth term of recursion relation
       if (kc.gt.1) iw(ic,kc-1,mc1) = -1
      endif
c**** third and fourth terms
      if (ic.gt.1) then
         ic = ic - 1
         iw(ic,kc,mc) = -1
         iw(ic,kc,mc1) = -1
      end if
      return
      end
      subroutine icrs2(iw,ix,kx,mx,icc,kc,mc,ichild)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix,kx,mx),ichild(6)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up ic
c
      ichild(1) = 4
      ic = icc - 1
      mc1 = mc + 1
      if( ic.ge.1) then
c**** first 2 terms of recursion relation
       ichild(2) = iw(ic,kc,mc)
       ichild(3) = iw(ic,kc,mc1)
c**** fifth term of recursion relation
       if (kc.gt.1) then
          ichild(6) = iw(ic,kc-1,mc1)
       else
          ichild(1) = 2
       end if
      endif
c**** third and fourth terms
      if (ic.gt.1) then
         ic = ic - 1
         ichild(4) = iw(ic,kc,mc)
         ichild(5) = iw(ic,kc,mc1)
      else
         ichild(1) = ichild(1) - 1
      end if
      return
      end
      subroutine icrsd(iw,ix,kx,mx,icc,kc,mc)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix+1,kx+1,mx+1)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up ic
c
      ic = icc - 1
      mc1 = mc + 1
      if( ic.ge.1) then
c**** first 2 terms of recursion relation
       iw(ic,kc,mc) = -1
       iw(ic,kc,mc1) = -1
c**** fifth term of recursion relation
       if (kc.gt.1) iw(ic,kc-1,mc1) = -1
      endif
c**** third and fourth terms
      if (ic.gt.1) then
         ic = ic - 1
         iw(ic,kc,mc) = -1
         iw(ic,kc,mc1) = -1
      end if
      return
      end
      subroutine icrsd2(iw,ix,kx,mx,icc,kc,mc,ichild)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix+1,kx+1,mx+1),ichild(6)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up ic
c
      ichild(1) = 4
      ic = icc - 1
      mc1 = mc + 1
      if( ic.ge.1) then
c**** first 2 terms of recursion relation
       ichild(2) = iw(ic,kc,mc)
       ichild(3) = iw(ic,kc,mc1)
c**** fifth term of recursion relation
       if (kc.gt.1) then
          ichild(6) = iw(ic,kc-1,mc1)
       else
          ichild(1) = 2
       end if
      endif
c**** third and fourth terms
      if (ic.gt.1) then
         ic = ic - 1
         ichild(4) = iw(ic,kc,mc)
         ichild(5) = iw(ic,kc,mc1)
      else
         ichild(1) = ichild(1) - 1
      end if
      return
      end
      subroutine indx4a(g,ivl,jvl,kvl,lvl,kmax,lmax0,ext2
     &,orbmo,iocy,nxx
     &,ipq,irpos,ispos,memreq,c,maxipq,temp)
c
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
c
      dimension g(icsize,memreq),orbmo(nb,nb)
      dimension c(nb,nxx,maxipq),temp(lmax0*lvl,kmax*kvl)
      logical ext2
c
c***** routine to do first part of the 4 index transformation.
c***** transforming (p q |r s) to (p q |iocc s)
c
c***** all rs are transformed for a batch of pq
c
c the labels i j k l are equivalent to p q r s
c
c kmax = number of kshell in batch
c kvl = number of basis functions in one k shell funtion
c irmax is the number of k/r basis functions in batch
c irb is the position of basis function in whole r/k list
c ir is the position of the basis function in the r/k batch
c (so that irb=ir+irpos)
c
c ipq is the counter over the ij basis functions in the entire
c                         ij batch
c answer stored in c(ipq,s) with ipq implying labels p and q
c and the iocc label implicit
      mem = 0
      irmax = kmax*kvl
      ismax = lmax0*lvl
      iocy1 = iocy + 1
      irpos1 = irpos + 1
      ispos1 = ispos + 1
      do 70 i = 1 , ivl
         do 60 j = 1 , jvl
            ipq = ipq + 1
c
            if (ext2) call vclr(temp,1,irmax*ismax)
            do 50 k = 1 , kvl
c
               do 40 l = 1 , lvl
c
                  mem = mem + 1
                  icount = 0
                  ir = k - kvl
                  do 30 kc = 1 , kmax
                     lmax = lmax0
                     if (ext2) lmax = kc
                     ir = ir + kvl
                     is = l - lvl
                     do 20 lc = 1 , lmax
                        is = is + lvl
                        icount = icount + 1
                        temp(is,ir) = g(icount,mem)
 20                  continue
 30               continue
 40            continue
 50         continue
            call mxmb(temp(1,1),ismax,1,orbmo(ispos1,iocy1),1,nb,
     +                c(irpos1,1,ipq),1,nb,irmax,ismax,nxx)
            call mxmb(temp(1,1),1,ismax,orbmo(irpos1,iocy1),1,nb,
     +                c(ispos1,1,ipq),1,nb,ismax,irmax,nxx)

 60      continue
 70   continue
c
      return
      end
      subroutine indx4d(orbmo,e,nocc,iocy,nxx,e2,b,nvir,c)
c
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/scra  /iso(mxshel,48),nt
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nnnnt,nt3,
     1itable(8,8),mult(8,8),irr(maxorb),imosb(8,maxorb),irrb(maxorb)
     2,it2(8,maxorb),it(8,maxorb),ibb(maxorb),ins(8),iss(8),ivsf(2,8)
c
      dimension orbmo(nb,nb),e(nb),temp(maxorb)
      dimension c(nvir,nvir),b(*)
c
c***** routine to do fourth part of the 4 index
c*****  transformation and the mp2 energy evaluation
c
c***** transforming (jocc q |iocc iavir) to (jocc ibvir|iocc iavir)
c
c***** the first array is in b; some of the second,
c***** the intermediate integrals, are in c(nvirtual**2)
c
c***** this routine is called for each iocc orbital
c***** and b contains all jocc; all q;
c***** and all iavir - but only one iocc
c
c**** symmetry adapt the molecular orbitals
c
      do 50 mo = 1 , nocc + nvir
         do 20 i = 1 , nb
            temp(i) = orbmo(i,mo)
            orbmo(i,mo) = 0.0d0
 20      continue
         do 40 i = 1 , nb
            do 30 nsym = 1 , nt
               if (it(nsym,i).eq.0) go to 40
               if (it(nsym,i).gt.0) orbmo(i,mo) = orbmo(i,mo)
     +             + temp(it(nsym,i))
               if (it(nsym,i).lt.0) orbmo(i,mo) = orbmo(i,mo)
     +             - temp(-it(nsym,i))
 30         continue
 40      continue
 50   continue
c
c**** normalise the symmetry adapted basis functions
c
      do 70 i = 1 , nb
         temp(i) = 0.0d0
         do 60 nsym = 1 , nt
            if (it(nsym,i).ne.0) then
               temp(i) = temp(i) + 1.0d0
            end if
 60      continue
         temp(i) = 1.0d0/temp(i)
 70   continue
      do 90 mo = 1 , nocc + nvir
         do 80 i = 1 , nb
            orbmo(i,mo) = orbmo(i,mo)*temp(i)
 80      continue
 90   continue
c
      ijapb = 0
      do 170 nxc = 1 , nxx
         iocc = iocy + nxc
         factor = 2.0d0
         do 160 jocc = 1 , iocc
c
c***** performing the final quarter transformation. care must be taken t
c***** unwind the 3/4 transformed integrals in the right order
c
            ijsym = mult(irr(iocc),irr(jocc))
            if (jocc.eq.iocc) factor = 1.0d0
            call vclr(c,1,nvir*nvir)
c
            do 130 nsymv = 1 , nt
               nsym = mult(ijsym,nsymv)
               do 120 ip = iss(nsym) + 1 , iss(nsym) + ins(nsym)
c
                  do 110 iavir = ivsf(1,nsymv) , ivsf(2,nsymv)
                     ijapb = ijapb + 1
                     do 100 ibvir = ivsf(1,nsym) , ivsf(2,nsym)
                        c(iavir,ibvir) = c(iavir,ibvir) + b(ijapb)
     +                     *orbmo(ip,nocc+ibvir)
 100                 continue
 110              continue
 120           continue
 130        continue
c
            e12 = 0.0d0
            do 150 iavir = 1 , nvir
               do 140 ibvir = 1 , nvir
c........................mp2 energy...............................
                  e12 = e12 + (2.0d0*c(iavir,ibvir)-c(ibvir,iavir))
     +                  *c(iavir,ibvir)
     +                  /(e(iocc)+e(jocc)-e(iavir+nocc)-e(ibvir+nocc))
 140           continue
 150        continue
            e2 = e2 + e12*factor
 160     continue
 170  continue
c
c**** if more than one batch of iocc
c**** and then unsymmetry-adapt the molecular orbitals
c
      if (iocc.lt.nocc) then
         do 210 mo = 1 , nocc + nvir
            do 180 i = 1 , nb
               temp(i) = orbmo(i,mo)
               orbmo(i,mo) = 0.0d0
 180        continue
            do 200 i = 1 , nb
               do 190 nsym = 1 , nt
                  if (it2(nsym,i).ne.0) then
                     if (it2(nsym,i).gt.0) orbmo(i,mo) = orbmo(i,mo)
     +                   + temp(it2(nsym,i))
                     if (it2(nsym,i).lt.0) orbmo(i,mo) = orbmo(i,mo)
     +                   - temp(-it2(nsym,i))
                  end if
 190           continue
 200        continue
 210     continue
      end if
c
      return
      end
      subroutine indxbc(ivl,jvl,imax,jmax0,ext1
     &,orbmo,nocc,iocy,nxx,ippos,iqpos,c,maxipq,d,nvir,b)
c
c
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/scra  /iso(mxshel,48),nt
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nnnnt,nt3,
     1itable(8,8),mult(8,8),irr(maxorb),imosb(8,maxorb),irrb(maxorb)
     2,it2(8,maxorb),it(8,maxorb),ibb(maxorb),ins(8),iss(8),ivsf(2,8)
c
      dimension orbmo(nb,nb)
      dimension c(nb,nxx,maxipq),b(*),d(nvir,nxx)
     &,ivtemp(8),iptemp(8),iqtemp(8)
      logical ext1
c
c***** routine to do second and third part of the 4 index
c*****  transformation.
c***** transforming (p q |iocc s) to (jocc q |iocc iavir)
c***** where q in the second integral is symmetry adapted
c
c***** the first array is in c; the second in b.
c***** the intermediate integrals are in d(nvirtual)
c
c***** this routine is called for each batch of pq until all
c***** batches have been done and b contains all jocc; all q;
c***** and all iavir - but only one iocc
c
c the labels i j k l are equivalent to p q r s
c
c ipq is the counter over the ij basis functions in the entire
c                         ij batch
c
      do 20 nsym = 1 , nt
         ivtemp(nsym) = 1 + ivsf(2,nsym) - ivsf(1,nsym)
 20   continue
c
c**** just a word about some of the arrays in symos that are different
c**** to those in cadpac. the virtual orbitals are sorted according to
c**** symmetry type to help vectorisation (see cwm thesis) and the
c**** starting point and ending point of virtuals of particular
c**** symmetries are given by the array ivsf. ins is the number of symme
c**** adapted atomic orbitals of a given symmetry type. the symmetry
c**** adapted orbitals are also ordered according to symmetry type and
c**** ibb gives the starting position of the first symmetry orbital of a
c**** new type, and this is used for counting purposes.
c
      ipq = 0
      do 150 ic = 1 , imax
         istart = (ic-1)*ivl + ippos
         jmax = jmax0
         if (ext1) jmax = ic
         do 140 jc = 1 , jmax
            jstart = (jc-1)*jvl + iqpos
            do 130 i = 1 , ivl
               ip = istart + i
c********* iptemp array used for counting purposes
               do 30 nsym = 1 , nt
                  iptemp(nsym) = ibb(iabs(it2(nsym,ip))) - 1
 30            continue
               do 120 j = 1 , jvl
                  iq = jstart + j
c********* iqtemp array used for counting purposes
                  do 40 nsym = 1 , nt
                     iqtemp(nsym) = ibb(iabs(it2(nsym,iq))) - 1
 40               continue
c
                  ipq = ipq + 1
                  ijapb = 0
c
c**** now form half-contracted integrals
c
                  call mxmaa(orbmo(1,nocc+1),nb,1,c(1,1,ipq),1,nb,d(1,1)
     +                       ,1,nvir,nvir,nb,nxx)
c
c**** d contains (p q| iavir iocc)
c
                  do 110 nxc = 1 , nxx
                     iocc = iocy + nxc
                     do 100 jocc = 1 , iocc
c
c**** we now form only the fully symmetric 3/4 transformed integrals
c
                        ijsym = mult(irr(iocc),irr(jocc))
c*************ijsym is symmetry of i (cross) j and this defines symmetry
c************* of    a (cross) p/q
                        do 90 nsymv = 1 , nt
c******* we've looped over symmetry types (this really represents the
c******* symmetry of the a virtual orbitals)
                           nsym = mult(ijsym,nsymv)
c******* the symmetry of the to-be-formed s.a. a.o is now fixed as nsym
                           if (it2(nsym,ip).ne.0) then
c******* if here then there are some s.a. a.o. contributed to by p of th
c******* symmetry
c
c**** etc
c
c***** now perform the calculations
c
                              ijapb1 = ijapb + iptemp(nsym)
     +                                 *ivtemp(nsymv)
                              if (it2(nsym,ip).lt.0) then
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                                 do 50 iavir = ivsf(1,nsymv) ,
     +                              ivsf(2,nsymv)
                                    ijapb1 = ijapb1 + 1
                                    b(ijapb1) = b(ijapb1) - d(iavir,nxc)
     +                                 *orbmo(iq,jocc)
 50                              continue
                              else
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                                 do 60 iavir = ivsf(1,nsymv) ,
     +                              ivsf(2,nsymv)
                                    ijapb1 = ijapb1 + 1
                                    b(ijapb1) = b(ijapb1) + d(iavir,nxc)
     +                                 *orbmo(iq,jocc)
 60                              continue
                              end if
                           end if
c
                           if (it2(nsym,iq).ne.0) then
                              ijapb1 = ijapb + iqtemp(nsym)
     +                                 *ivtemp(nsymv)
                              if (it2(nsym,iq).lt.0) then
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                                 do 70 iavir = ivsf(1,nsymv) ,
     +                              ivsf(2,nsymv)
                                    ijapb1 = ijapb1 + 1
                                    b(ijapb1) = b(ijapb1) - d(iavir,nxc)
     +                                 *orbmo(ip,jocc)
 70                              continue
                              else
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                                 do 80 iavir = ivsf(1,nsymv) ,
     +                              ivsf(2,nsymv)
                                    ijapb1 = ijapb1 + 1
                                    b(ijapb1) = b(ijapb1) + d(iavir,nxc)
     +                                 *orbmo(ip,jocc)
 80                              continue
                              end if
                           end if
c
                           ijapb = ijapb + ivtemp(nsymv)*ins(nsym)
c
 90                     continue
 100                 continue
c
c
 110              continue
 120           continue
 130        continue
 140     continue
 150  continue
c
      return
      end
      subroutine kcrs(iw,ix,kx,mx,ic,kcc,mc)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix,kx,mx)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up kc
c
      kc = kcc - 1
      mc1 = mc + 1
c**** first 2 terms of recursion relation
      iw(ic,kc,mc) = -1
      iw(ic,kc,mc1) = -1
c**** fifth term of recursion relation
      if (ic.gt.1) iw(ic-1,kc,mc1) = -1
c**** third and fourth terms
      if (kc.gt.1) then
         kc = kc - 1
         iw(ic,kc,mc) = -1
         iw(ic,kc,mc1) = -1
      end if
      return
      end
      subroutine kcrs2(iw,ix,kx,mx,ic,kcc,mc,ichild)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix,kx,mx),ichild(6)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up kc
c
      kc = kcc - 1
      mc1 = mc + 1
      ichild(1) = 8
c**** first 2 terms of recursion relation
      ichild(2) = iw(ic,kc,mc)
      ichild(3) = iw(ic,kc,mc1)
c**** fifth term of recursion relation
      if (ic.gt.1) then
         ichild(6) = iw(ic-1,kc,mc1)
      else
         ichild(1) = 6
      end if
c**** third and fourth terms
      if (kc.gt.1) then
         kc = kc - 1
         ichild(4) = iw(ic,kc,mc)
         ichild(5) = iw(ic,kc,mc1)
      else
         ichild(1) = ichild(1) - 1
      end if
      return
      end
      subroutine kcrsd(iw,ix,kx,mx,ic,kcc,mc)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix+1,kx+1,mx+1)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up kc
c
      kc = kcc - 1
      mc1 = mc + 1
c**** first 2 terms of recursion relation
      iw(ic,kc,mc) = -1
      iw(ic,kc,mc1) = -1
c**** fifth term of recursion relation
      if (ic.gt.1) iw(ic-1,kc,mc1) = -1
c**** third and fourth terms
      if (kc.gt.1) then
         kc = kc - 1
         iw(ic,kc,mc) = -1
         iw(ic,kc,mc1) = -1
      end if
      return
      end
      subroutine kcrsd2(iw,ix,kx,mx,ic,kcc,mc,ichild)
c
      implicit REAL  (a-h,o-z)
      dimension iw(ix+1,kx+1,mx+1),ichild(6)
c
c**** this routine gives all the child integrals needed
c**** to form integral i(ic,kc,mc) by building up kc
c
      kc = kcc - 1
      mc1 = mc + 1
      ichild(1) = 8
c**** first 2 terms of recursion relation
      ichild(2) = iw(ic,kc,mc)
      ichild(3) = iw(ic,kc,mc1)
c**** fifth term of recursion relation
      if (ic.gt.1) then
         ichild(6) = iw(ic-1,kc,mc1)
      else
         ichild(1) = 6
      end if
c**** third and fourth terms
      if (kc.gt.1) then
         kc = kc - 1
         ichild(4) = iw(ic,kc,mc)
         ichild(5) = iw(ic,kc,mc1)
      else
         ichild(1) = ichild(1) - 1
      end if
      return
      end
_EXTRACT(tst1s,itanium,linux64)
      subroutine kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      logical ext2
c
      if (ext2 .and. lmax1.eq.kc) then
         do 20 lc = lmin1 , lmax1 - 1
            jc1 = jc1 + 1
            jcount = jcount + nc(kc)*nc(lc)
 20      continue
         jc1 = jc1 + 1
         jcount = jcount + nc(kc)*(nc(kc)+1)/2
      else
         do 30 lc = lmin1 , lmax1
            jc1 = jc1 + 1
            jcount = jcount + nc(kc)*nc(lc)
 30      continue
      end if
      return
      end
_ENDEXTRACT
      subroutine lmiss(kpmax,lpmax,ext6,jcount)
      implicit REAL  (a-h,o-z)
      logical ext6
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
      if (ext6) then
         jcount = jcount + kpmax*(kpmax+1)/2
      else
         jcount = jcount + kpmax*lpmax
      end if
c
      return
      end
      subroutine memlc1(k00,k01,k02,k03,k04,k05,k06,kstart,ip1)
      implicit REAL  (a-h,o-z)
c
      k00 = kstart + ip1
      k01 = k00 + ip1
      k02 = k01 + ip1
      k03 = k02 + ip1
      k04 = k03 + ip1
      k05 = k04 + ip1
      k06 = k05 + ip1*3
      return
      end
      subroutine memlc2(k00,k01,k02,k03,k04,k05,k06,k07,kstart,ip1)
      implicit REAL  (a-h,o-z)
c
      k00 = kstart + ip1
      k01 = k00 + ip1
      k02 = k01 + ip1
      k03 = k02 + ip1
      k04 = k03 + ip1
      k05 = k04 + ip1
      k06 = k05 + ip1
      k07 = k06 + ip1*3
      return
      end
      function nposf(n1,nvl,nd,np,ns)
      implicit REAL  (a-h,o-z)
      if (nvl.eq.1) nposf = n1 - 1
      if (nvl.eq.3) nposf = ns + (n1-ns-1)*3
      if (nvl.eq.6) nposf = ns + np*3 + (n1-ns-np-1)*6
      return
      end
      subroutine nucder(grad,drg,dnuc)
c-----------------------------------------------------------
c     nuclear derivatives
c
c-----------------------------------------------------------
      implicit REAL  (a-h,o-z)
c
c      zan( ) - nuclear charge
c      c( ) - cartesian positions
c
INCLUDE(common/sizes)
      common/iofile/iread,iwr
      common/data/zan(maxat),c(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),
     & nat,nb,ns,np,nsp
      logical glog1,glog2,glog3,glog5,glog6,gdbg
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog3,glog6,gdbg(10)
c
c      dnuc - nuclear contribution, grad - total gradient
c
      dimension drg(nat,nat),dnuc(3,nat),grad(nat,3)
      data zzero,one /0.0d0,1.0d0/
c
c     dnuc will hold the answer
c
      do 30 k = 1 , nat
         do 20 kk = 1 , 3
            dnuc(kk,k) = zzero
 20      continue
 30   continue
c
c     calculate r and 1/r**2
c
      drg(1,1) = zzero
      do 60 k = 2 , nat
         drg(k,k) = zzero
         k1 = k - 1
         do 50 l = 1 , k1
            rkl = zzero
            do 40 i = 1 , 3
               rkl = rkl + (c(i,k)-c(i,l))**2
 40         continue
            drg(k,l) = -one/rkl
            drg(l,k) = dsqrt(rkl)
 50      continue
 60   continue
c
c     nuclear contribution to gradient
c
      do 110 kk = 1 , 3
         do 80 k = 2 , nat
            zak = zan(k)
            km1 = k - 1
            do 70 l = 1 , km1
               zal = zan(l)
               pkl = (c(kk,k)-c(kk,l))/drg(l,k)
               dnuc(kk,k) = dnuc(kk,k) + pkl*drg(k,l)*zak*zal
 70         continue
 80      continue
         nat1 = nat - 1
         do 100 k = 1 , nat1
            zak = zan(k)
            kp1 = k + 1
            do 90 l = kp1 , nat
               zal = zan(l)
               pkl = (c(kk,k)-c(kk,l))/drg(k,l)
               dnuc(kk,k) = dnuc(kk,k) + pkl*drg(l,k)*zak*zal
 90         continue
 100     continue
 110  continue
c
      if (gdbg(8)) then
         write (iwr,6030)
         write (iwr,6020)
         do 120 i = 1 , nat
            write (iwr,6010) i , dnuc(1,i) , dnuc(2,i) , dnuc(3,i)
 120     continue
      end if
c
c     add dnuc to total gradient
c
      do 140 k = 1 , nat
         do 130 kk = 1 , 3
            grad(k,kk) = grad(k,kk) + dnuc(kk,k)
 130     continue
 140  continue
      return
 6010 format (1x,i3,3(4x,d16.8))
 6020 format (1x,'atom',11x,'x',19x,'y',19x,'z')
 6030 format (//,' nuclear contribution to gradient'/)
      end
      subroutine oneda(zmem,xija,xijb,xijc,xijp1,xijp2
     &,xijp3,nij,dens,rdens,grad,grtemp
     &,ext1,i1,i2,j1,j2,ninty1,kinup,koveru,kindow,koverd)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c**** the fixed parameters in common/indexd/ are explained in the
c**** block data statement (in bdata.f?)
c**** the memory locations and requirements were established in
c**** cidone and are stored in common/cider/
c
c**** routine finds the contribution to gradient of 1-electron derivativ
c**** integrals (except nuclear attraction - see onedb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto)
     &,xija(nij),xijb(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,xijc(nij),naind(3),ncind(3),grad(na,3),grtemp(na,3)
      dimension rdens(nb,nb),dens(nb,nb)
c********* rdens = energy weighted density matrix
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipszt,icszt
     &,ipsize,icsize,memcon
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/defunk/cobas(maxorb,3),nd
      logical ext1,ext5,term2
c
      ivlx = nvl(ix)
      ivly = nvl(ix-1)
      jvlx = nvl(kx)
      ipos = nposf(i1,ivly,nd,np,ns)
      jpos = nposf(j1,jvlx,nd,np,ns)
      jmax = j2
      icount = 0
c
c**** the main contracted loop over i and j shells
c
c**** note no symmetry or magniotude based tests
c
      do 40 ic = i1 , i2
         if (ext1) jmax = ic
         do 30 jc = j1 , jmax
            ext5 = ext1 .and. ic.eq.jc
            itprim = nc(ic)*nc(jc)
            if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
            do 20 iprim = 1 , itprim
               icount = icount + 1
c
c***** prepare all parameter arrays needed to form
c***** restricted list of intermediate primitive integrals
c
               zmem(icount,3) = xijp1(icount) - cobas(ic,1)
               zmem(icount,4) = xijp2(icount) - cobas(ic,2)
               zmem(icount,5) = xijp3(icount) - cobas(ic,3)
               zmem(icount,6) = xijp1(icount) - cobas(jc,1)
               zmem(icount,7) = xijp2(icount) - cobas(jc,2)
               zmem(icount,8) = xijp3(icount) - cobas(jc,3)
 20         continue
 30      continue
 40   continue
      do 50 num = 1 , nij
         zmem(num,9) = xijb(num)/xijc(num)
         zmem(num,10) = xija(num)*xijb(num)/(1.0d0-xija(num)*xijc(num))
 50   continue
c
c
c***** now do the loops for forming the intermediate
c***** primitive integrals
c
c*** z(3-5) pi-ai ::: z(6-8) pi-bi ::: z(7) c+d :::
c*** xija  1/(2a+2b) ::: xijb  2ab/(a+b) ::: xijc  2a :::
c*** z(9) xijb/2a ::: z(10) xijb/2b ::: z(2) (s|t|s) :::
c*** z(1) (s|s) :::
c
c*** loop over number of distinct integral classes
c*** needed to form the target (s|s) integrals
c
      memcur = 10
      do 210 n = 1 , ninty1
         ityp = intgrl(1,n)
         ktyp = intgrl(2,n)
         ityp1 = ityp - 1
         ktyp1 = ktyp - 1
         ityp2 = ityp - 2
         ktyp2 = ktyp - 2
         ivl = nvl(ityp)
         kvl = nvl(ktyp)
         if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
         if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
         mem1 = intgrl(3,n)
         mem2 = intgrl(4,n)
         mem3 = intgrl(5,n)
         nrecur = 6
         if (ktyp.eq.1 .and. ityp.gt.2) nrecur = 2
         if (ktyp.eq.1 .and. ityp.le.2) nrecur = 1
         if (ktyp.eq.2 .and. ityp.eq.1) nrecur = 3
         if (ktyp.gt.2 .and. ityp.eq.1) nrecur = 4
         if (ktyp.eq.2 .and. ityp.gt.1) nrecur = 5
c
c now:
c loop over number of indices associated with i position
         do 200 i = 1 , ivl
            nr = nrecur
            if (nr.lt.3) then
c***** reduction is at the i position
c
               index = ispdf1(i,1)
               indpa = index + 2
               natred = ispdf1(i,2)
               itemp = (natred-1)*kvl
               jmem1 = itemp + mem1
               if (nr.ne.1) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                  nai = ixyz(index,natred,ityp2)
                  if (nai.eq.0) then
                     nr = nr - 1
                     go to 60
                  end if
                  naind(2) = ixyz(2,natred,ityp2) + 1
                  naind(3) = ixyz(3,natred,ityp2) + 1
                  naind(index) = naind(index) - 1
                  idoubr = iarray(naind(2),naind(3))
                  itemp2 = (idoubr-1)*kvl
                  jmem2 = itemp2 + mem2
               end if
            else
c
c**** reduction is at the k position
               itemp = (i-1)*kvl1
               jmem1 = itemp + mem1
               itemp = (i-1)*kvl2
               jmem2 = itemp + mem2
               if (nr.ne.3 .and. nr.ne.4) then
                  naind(1) = ixyz(1,i,ityp1) + 1
                  naind(2) = ixyz(2,i,ityp1) + 1
                  naind(3) = ixyz(3,i,ityp1) + 1
               end if
            end if
c
c now:
c loop over number of indices associated with k position
c
 60         do 190 k = 1 , kvl
               memcur = memcur + 1
               if (nr.lt.3) then
                  imem1 = jmem1 + k
               else
                  index = ispdf1(k,1)
                  indpb = index + 5
                  natred = ispdf1(k,2)
                  imem1 = jmem1 + natred
               end if
               go to (70,90,110,130,150,170) , nr
 70            do 80 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpa)*zmem(num,imem1)
 80            continue
               go to 190
 90            imem2 = jmem2 + k
               do 100 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpa)*zmem(num,imem1)
     +                               + dfloat(nai)*xija(num)
     +                               *zmem(num,imem2)
 100           continue
               go to 190
 110           do 120 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
 120           continue
               go to 190
 130           nci = ixyz(index,natred,ktyp2)
               if (nci.eq.0) go to 110
               ncind(2) = ixyz(2,natred,ktyp2) + 1
               ncind(3) = ixyz(3,natred,ktyp2) + 1
               ncind(index) = ncind(index) - 1
               imem2 = jmem2 + iarray(ncind(2),ncind(3))
               do 140 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + dfloat(nci)*xija(num)
     +                               *zmem(num,imem2)
 140           continue
               go to 190
 150           if (naind(index).eq.1) go to 110
               naind(index) = naind(index) - 1
               imem3 = (iarray(naind(2),naind(3))-1)*kvl1 + natred +
     +                 mem3
               do 160 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + dfloat(naind(index))*xija(num)
     +                               *zmem(num,imem3)
 160           continue
               naind(index) = naind(index) + 1
               go to 190
 170           if (naind(index).eq.1) go to 130
               nci = ixyz(index,natred,ktyp2)
               if (nci.eq.0) go to 150
               ncind(2) = ixyz(2,natred,ktyp2) + 1
               ncind(3) = ixyz(3,natred,ktyp2) + 1
               ncind(index) = ncind(index) - 1
               imem2 = jmem2 + iarray(ncind(2),ncind(3))
               naind(index) = naind(index) - 1
               imem3 = (iarray(naind(2),naind(3))-1)*kvl1 + natred +
     +                 mem3
               do 180 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + dfloat(nci)*xija(num)
     +                               *zmem(num,imem2)
     +                               + dfloat(naind(index))*xija(num)
     +                               *zmem(num,imem3)
 180           continue
               naind(index) = naind(index) + 1
 190        continue
 200     continue
 210  continue
c
c*** loop over number of distinct integral classes
c*** needed to form the target (s|t|s) integrals
c
      do 370 n = ninty1 + 1 , ninty
         ityp = intgrl(1,n)
         ktyp = intgrl(2,n)
         ityp1 = ityp - 1
         ktyp1 = ktyp - 1
         ityp2 = ityp - 2
         ktyp2 = ktyp - 2
         ivl = nvl(ityp)
         kvl = nvl(ktyp)
         if (ktyp1.gt.0) kvl1 = nvl(ktyp1)
         if (ktyp2.gt.0) kvl2 = nvl(ktyp2)
c
         mem1 = intgrl(3,n)
         mem2 = intgrl(4,n)
         mem3 = intgrl(5,n)
         mem4 = intgrl(6,n)
         mem5 = intgrl(7,n)
         nrecur = 6
         if (ktyp.eq.1 .and. ityp.gt.2) nrecur = 2
         if (ktyp.eq.1 .and. ityp.le.2) nrecur = 1
         if (ktyp.eq.2 .and. ityp.eq.1) nrecur = 3
         if (ktyp.gt.2 .and. ityp.eq.1) nrecur = 4
         if (ktyp.eq.2 .and. ityp.gt.1) nrecur = 5
c
c now:
c loop over number of indices associated with i position
         do 360 i = 1 , ivl
            nr = nrecur
            if (nr.lt.3) then
c***** reduction is at the i position
c
               index = ispdf1(i,1)
               indpa = index + 2
               natred = ispdf1(i,2)
               itemp = (natred-1)*kvl
               jmem1 = itemp + mem1
               if (nr.ne.1) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                  nai = ixyz(index,natred,ityp2)
                  if (nai.eq.0) then
                     nr = nr - 1
                     go to 220
                  end if
                  naind(2) = ixyz(2,natred,ityp2) + 1
                  naind(3) = ixyz(3,natred,ityp2) + 1
                  naind(index) = naind(index) - 1
                  idoubr = iarray(naind(2),naind(3))
                  itemp2 = (idoubr-1)*kvl
                  jmem3 = itemp2 + mem3
                  jmem4 = itemp2 + mem4
               end if
            else
c
c**** reduction is at the k position
               itemp = (i-1)*kvl1
               jmem1 = itemp + mem1
               itemp = (i-1)*kvl2
               jmem3 = itemp + mem3
               jmem4 = itemp + mem4
               if (nr.ne.3 .and. nr.ne.4) then
                  naind(1) = ixyz(1,i,ityp1) + 1
                  naind(2) = ixyz(2,i,ityp1) + 1
                  naind(3) = ixyz(3,i,ityp1) + 1
               end if
            end if
c
c now:
c loop over number of indices associated with k position
c
 220        do 350 k = 1 , kvl
               memcur = memcur + 1
               mem2 = mem2 + 1
               if (nr.lt.3) then
                  imem1 = jmem1 + k
               else
                  index = ispdf1(k,1)
                  indpb = index + 5
                  natred = ispdf1(k,2)
                  imem1 = jmem1 + natred
               end if
               go to (230,250,270,290,310,330) , nr
 230           do 240 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpa)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
 240           continue
               go to 350
 250           imem3 = jmem3 + k
               imem4 = jmem4 + k
               do 260 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpa)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
     +                               + dfloat(nai)
     +                               *(xija(num)*zmem(num,imem3)
     +                               -zmem(num,9)*zmem(num,imem4))
 260           continue
               go to 350
 270           do 280 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
 280           continue
               go to 350
 290           nci = ixyz(index,natred,ktyp2)
               if (nci.eq.0) go to 270
               ncind(2) = ixyz(2,natred,ktyp2) + 1
               ncind(3) = ixyz(3,natred,ktyp2) + 1
               ncind(index) = ncind(index) - 1
               imem3 = jmem3 + iarray(ncind(2),ncind(3))
               imem4 = jmem4 + iarray(ncind(2),ncind(3))
               do 300 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
     +                               + dfloat(nci)
     +                               *(xija(num)*zmem(num,imem3)
     +                               -zmem(num,10)*zmem(num,imem4))
 300           continue
               go to 350
 310           if (naind(index).eq.1) go to 270
               naind(index) = naind(index) - 1
               imem5 = (iarray(naind(2),naind(3))-1)*kvl1 + natred +
     +                 mem5
               do 320 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
     +                               + dfloat(naind(index))*xija(num)
     +                               *zmem(num,imem5)
 320           continue
               naind(index) = naind(index) + 1
               go to 350
 330           if (naind(index).eq.1) go to 290
               nci = ixyz(index,natred,ktyp2)
               if (nci.eq.0) go to 310
               ncind(2) = ixyz(2,natred,ktyp2) + 1
               ncind(3) = ixyz(3,natred,ktyp2) + 1
               ncind(index) = ncind(index) - 1
               imem3 = jmem3 + iarray(ncind(2),ncind(3))
               imem4 = jmem4 + iarray(ncind(2),ncind(3))
               naind(index) = naind(index) - 1
               imem5 = (iarray(naind(2),naind(3))-1)*kvl1 + natred +
     +                 mem5
               do 340 num = 1 , nij
                  zmem(num,memcur) = zmem(num,indpb)*zmem(num,imem1)
     +                               + xijb(num)*zmem(num,mem2)
     +                               + dfloat(nci)
     +                               *(xija(num)*zmem(num,imem3)
     +                               -zmem(num,10)*zmem(num,imem4))
     +                               + dfloat(naind(index))*xija(num)
     +                               *zmem(num,imem5)
 340           continue
               naind(index) = naind(index) + 1
 350        continue
 360     continue
 370  continue
c
c****** form the gamma which will multiply the derivative
c****** overlap integrals
c
      mgamma = memcur
      do 420 i = 1 , ivly
         do 410 j = 1 , jvlx
            icount = 0
            memcur = memcur + 1
            do 400 ic = i1 , i2
               id = (ic-i1)*ivly + i + ipos
               if (ext1) jmax = ic
c
               if (lx.eq.1) then
                  do 380 jd = j1 , jmax
                     icount = icount + 1
                     zmem(icount,memcur) = rdens(jd,id)
 380              continue
               else
                  jd = j + jpos - jvlx
                  do 390 jc = j1 , jmax
                     jd = jd + jvlx
                     icount = icount + 1
                     zmem(icount,memcur) = rdens(jd,id)
 390              continue
               end if
 400        continue
 410     continue
 420  continue
c
      memcur = memcur + 1
      mems = memcur
      moveru = memcur
c
c**** form contracted overup
c
      do 470 m = 1 , jvlx*ivlx
         kin = koveru + m
         memcur = memcur + 1
         do 430 num = 1 , nij
            zmem(num,mems) = xijc(num)*zmem(num,kin)
 430     continue
         icount = 0
         jcount = 0
         do 460 ic = i1 , i2
            if (ext1) jmax = ic
            do 450 jc = j1 , jmax
               ext5 = ext1 .and. ic.eq.jc
               itprim = nc(ic)*nc(jc)
               if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
               jcount = jcount + 1
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
               zmem(jcount,memcur) = 0.0d0
               do 440 iprim = 1 , itprim
                  icount = icount + 1
                  zmem(jcount,memcur) = zmem(jcount,memcur)
     +                                  + zmem(icount,mems)
 440           continue
 450        continue
 460     continue
 470  continue
      if (ix.gt.2) then
c**** form contracted kindown
         moverd = memcur
         do 510 m = 1 , jvlx*ivlx
            kin = koverd + m
            memcur = memcur + 1
            icount = 0
            jcount = 0
            do 500 ic = i1 , i2
               if (ext1) jmax = ic
               do 490 jc = j1 , jmax
                  ext5 = ext1 .and. ic.eq.jc
                  itprim = nc(ic)*nc(jc)
                  if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
                  jcount = jcount + 1
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
                  zmem(jcount,memcur) = 0.0d0
                  do 480 iprim = 1 , itprim
                     icount = icount + 1
                     zmem(jcount,memcur) = zmem(jcount,memcur)
     +                  + zmem(icount,kin)
 480              continue
 490           continue
 500        continue
 510     continue
      end if
c
c**** should now have the required things for gradient contribution
c
      i1typ1 = ix - 2
      mgrad = memcur
      do 540 icart = 1 , 3
         memcur = memcur + 1
         mem0 = mgamma
         call vclr(zmem(1,memcur),1,nij)
         do 530 i = 1 , ivly
            term2 = .false.
            if (i1typ1.eq.0) then
               naind(2) = 1
               naind(3) = 1
            else
               naind(1) = ixyz(1,i,i1typ1) + 1
               naind(2) = ixyz(2,i,i1typ1) + 1
               naind(3) = ixyz(3,i,i1typ1) + 1
            end if
            naind(icart) = naind(icart) + 1
            imem1 = moveru + (iarray(naind(2),naind(3))-1)*jvlx
            if (ix.gt.2) then
               naind(icart) = naind(icart) - 2
               nai = naind(icart)
               term2 = nai.gt.0
               if (term2) imem00 = moverd +
     +                             (iarray(naind(2),naind(3))-1)*jvlx
            end if
            do 520 j = 1 , jvlx
               jmem1 = imem1 + j
               jmem00 = imem00 + j
               mem0 = mem0 + 1
               call hrzd0(zmem(1,memcur),zmem(1,mem0),zmem(1,jmem1),nai,
     +                    zmem(1,jmem00),jcount,term2)
 520        continue
 530     continue
 540  continue
c
c**** now form the gradient
c
      jcount = 0
      do 570 ic = i1 , i2
         if (ext1) jmax = ic
         iatom = nwa(nup(ic))
         do 560 jc = j1 , jmax
            jcount = jcount + 1
            jatom = nwa(nup(jc))
            do 550 icart = 1 , 3
               grad(iatom,icart) = grad(iatom,icart)
     +                             - zmem(jcount,icart+mgrad)
               grad(jatom,icart) = grad(jatom,icart)
     +                             + zmem(jcount,icart+mgrad)
 550        continue
 560     continue
 570  continue
c
c****** form the gamma which will multiply the derivative
c****** kinetic integrals
c
      memcur = mgamma
      do 620 i = 1 , ivly
         do 610 j = 1 , jvlx
            icount = 0
            memcur = memcur + 1
            do 600 ic = i1 , i2
               id = (ic-i1)*ivly + i + ipos
               if (ext1) jmax = ic
c
               if (lx.eq.1) then
                  do 580 jd = j1 , jmax
                     icount = icount + 1
                     zmem(icount,memcur) = dens(jd,id)
 580              continue
               else
                  jd = j + jpos - jvlx
                  do 590 jc = j1 , jmax
                     jd = jd + jvlx
                     icount = icount + 1
                     zmem(icount,memcur) = dens(jd,id)
 590              continue
               end if
 600        continue
 610     continue
 620  continue
c
      memcur = memcur + 1
      mems = memcur
      mkinup = memcur
c
c**** form contracted kinup
c
      do 670 m = 1 , jvlx*ivlx
         kin = kinup + m
         memcur = memcur + 1
         do 630 num = 1 , nij
            zmem(num,mems) = xijc(num)*zmem(num,kin)
 630     continue
         icount = 0
         jcount = 0
         do 660 ic = i1 , i2
            if (ext1) jmax = ic
            do 650 jc = j1 , jmax
               ext5 = ext1 .and. ic.eq.jc
               itprim = nc(ic)*nc(jc)
               if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
               jcount = jcount + 1
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
               zmem(jcount,memcur) = 0.0d0
               do 640 iprim = 1 , itprim
                  icount = icount + 1
                  zmem(jcount,memcur) = zmem(jcount,memcur)
     +                                  + zmem(icount,mems)
 640           continue
 650        continue
 660     continue
 670  continue
      if (ix.gt.2) then
c**** form contracted kindown
         mkind = memcur
         do 710 m = 1 , jvlx*ivlx
            kin = kindow + m
            memcur = memcur + 1
            icount = 0
            jcount = 0
            do 700 ic = i1 , i2
               if (ext1) jmax = ic
               do 690 jc = j1 , jmax
                  ext5 = ext1 .and. ic.eq.jc
                  itprim = nc(ic)*nc(jc)
                  if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
                  jcount = jcount + 1
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
                  zmem(jcount,memcur) = 0.0d0
                  do 680 iprim = 1 , itprim
                     icount = icount + 1
                     zmem(jcount,memcur) = zmem(jcount,memcur)
     +                  + zmem(icount,kin)
 680              continue
 690           continue
 700        continue
 710     continue
      end if
c
c**** should now have the required things for gradient contribution
c
      i1typ1 = ix - 2
      mgrad = memcur
      do 740 icart = 1 , 3
         memcur = memcur + 1
         mem0 = mgamma
         call vclr(zmem(1,memcur),1,nij)
         do 730 i = 1 , ivly
            term2 = .false.
            if (i1typ1.eq.0) then
               naind(2) = 1
               naind(3) = 1
            else
               naind(1) = ixyz(1,i,i1typ1) + 1
               naind(2) = ixyz(2,i,i1typ1) + 1
               naind(3) = ixyz(3,i,i1typ1) + 1
            end if
            naind(icart) = naind(icart) + 1
            imem1 = mkinup + (iarray(naind(2),naind(3))-1)*jvlx
            if (ix.gt.2) then
               naind(icart) = naind(icart) - 2
               nai = naind(icart)
               term2 = nai.gt.0
               if (term2) imem00 = mkind + (iarray(naind(2),naind(3))-1)
     +                             *jvlx
            end if
            do 720 j = 1 , jvlx
               jmem1 = imem1 + j
               jmem00 = imem00 + j
               mem0 = mem0 + 1
               call hrzd0(zmem(1,memcur),zmem(1,mem0),zmem(1,jmem1),nai,
     +                    zmem(1,jmem00),jcount,term2)
 720        continue
 730     continue
 740  continue
c
c**** now form the gradient
c
      jcount = 0
      do 770 ic = i1 , i2
         if (ext1) jmax = ic
         iatom = nwa(nup(ic))
         do 760 jc = j1 , jmax
            jcount = jcount + 1
            jatom = nwa(nup(jc))
            do 750 icart = 1 , 3
               grtemp(iatom,icart) = grtemp(iatom,icart)
     +                               + zmem(jcount,icart+mgrad)
               grtemp(jatom,icart) = grtemp(jatom,icart)
     +                               - zmem(jcount,icart+mgrad)
 750        continue
 760     continue
 770  continue
c
c
      return
      end
      subroutine onedb(zmem,xija,xijk,xijc,xijp1,xijp2
     &,xijp3,nij,dens,grad,grtemp
     &,ext1,i1,i2,j1,j2)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c**** the fixed parameters in common/indexd/ are explained in the
c**** block data statement (in bdata.f?)
c**** the memory locations and requirements were established in
c**** cidonb and are stored in common/cider/ and cidriv
c
c***** routine calculates derivative nuclear attraction integrals
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      dimension zmem(ipsize,mempto)
     &,xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij)
     &,xijc(nij),naind(3),grad(na,3),grtemp(na,3)
     &,dens(nb,nb)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz0(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtempx,ipszt,icszt
     &,ipsize,icsize,memcon
c
      common/cidriv/ncona,icona(5,20),nconb,iconb(5,20),nconc
     &,iconc(5,20),nconx,iconx(5,20),ihz(5,50),mgradx
c
      common/indexd/ispdf1(21,2),nvl(6)
     &,ixyz(3,21,5),iarray(6,6)
      common/iofile/iread,iwr
c
      common/defunk/cobas(maxorb,3),nd
      logical ext1,ext5
c
      ivlx = nvl(ix)
      jvlx = nvl(kx)
      ipos = nposf(i1,ivlx,nd,np,ns)
      jpos = nposf(j1,jvlx,nd,np,ns)
      jmax = j2
c
c     icorig = 0
c
c**** the main contracted loop over i and j shells
c
      do 370 natom = 1 , na
c
         icount = 0
         do 40 ic = i1 , i2
            if (ext1) jmax = ic
            do 30 jc = j1 , jmax
               ext5 = ext1 .and. ic.eq.jc
               itprim = nc(ic)*nc(jc)
               if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
c***** loop over the i and j primitives in this ic,jc basis
c***** function pairing
c
               do 20 iprim = 1 , itprim
                  icount = icount + 1
c
c***** prepare all parameter arrays needed to form
c***** restricted list of intermediate primitive integrals
c
                  zmem(icount,1) = xijp1(icount) - cobas(ic,1)
                  zmem(icount,2) = xijp2(icount) - cobas(ic,2)
                  zmem(icount,3) = xijp3(icount) - cobas(ic,3)
                  zmem(icount,4) = xijp1(icount) - coord(1,natom)
                  zmem(icount,5) = xijp2(icount) - coord(2,natom)
                  zmem(icount,6) = xijp3(icount) - coord(3,natom)
 20            continue
 30         continue
 40      continue
c
         call snterd(zmem(1,4),zmem(1,5),zmem(1,6),xija,xijk,zmem(1,7),
     +               zmem(1,8),zmem(1,9+mx),mx+1,nij,ipsize,
     +               charge(natom))
c
         do 50 num = 1 , nij
            zmem(num,7) = 2.0d0*xija(num)
 50      continue
c
c
c***** now do the loops for forming the intermediate
c***** primitive integrals
c
c*** z(1-3) pi-ai ::: z(4-6) pi-ci ::: z(7) 2a+2b :::
c*** xijc  2a ::: z(8-?) (s|a|s)m
c
c*** loop over number of distinct integral classes
c*** needed to form the target (s|s) integrals
c
         memcur = 8 + mx
         do 110 n = 1 , ninty
            ityp = intgrl(1,n)
c           ityp1 = ityp - 1
            ityp2 = ityp - 2
            ivl = nvl(ityp)
            kvl = 1
c
            mem1 = intgrl(2,n)
            mem2 = intgrl(3,n)
            mem3 = intgrl(4,n)
            mem4 = intgrl(5,n)
            nrecur = 2
            if (ityp.le.2) nrecur = 1
c
c now:
c loop over number of indices associated with i position
            do 100 i = 1 , ivl
               nr = nrecur
c***** reduction is at the i position
c
               index = ispdf1(i,1)
               indpa = index
               indpc = index + 3
               natred = ispdf1(i,2)
               itemp = (natred-1)*kvl
               jmem1 = itemp + mem1
               jmem2 = itemp + mem2
               if (nr.ne.1) then
c***  double reduction at the iposition is a possibility for
c***  this integral class
                  nai = ixyz(index,natred,ityp2)
                  if (nai.eq.0) then
                     nr = nr - 1
                     go to 60
                  end if
                  naind(2) = ixyz(2,natred,ityp2) + 1
                  naind(3) = ixyz(3,natred,ityp2) + 1
                  naind(index) = naind(index) - 1
                  idoubr = iarray(naind(2),naind(3))
                  itemp2 = (idoubr-1)*kvl
                  jmem3 = itemp2 + mem3
                  jmem4 = itemp2 + mem4
               end if
c
c now:
c loop over number of indices associated with k position
c
 60            do 90 k = 1 , kvl
                  memcur = memcur + 1
                  imem1 = jmem1 + k
                  imem2 = jmem2 + k
                  if (nr.eq.1) then
                     do 70 num = 1 , nij
                        zmem(num,memcur) = zmem(num,indpa)
     +                     *zmem(num,imem1) - zmem(num,indpc)
     +                     *zmem(num,imem2)
 70                  continue
                  else
                     imem3 = jmem3 + k
                     imem4 = jmem4 + k
                     do 80 num = 1 , nij
                        zmem(num,memcur) = zmem(num,indpa)
     +                     *zmem(num,imem1) - zmem(num,indpc)
     +                     *zmem(num,imem2) + (dfloat(nai)/zmem(num,7))
     +                     *(zmem(num,imem3)-zmem(num,imem4))
 80                  continue
                  end if
 90            continue
 100        continue
 110     continue
c
c**** form contracted intermediates (ix-1 | jx) and (ix | jx-1)
c
         do 160 n = 1 , nconx
            imem1 = iconx(3,n)
            do 150 i = 1 , nvl(iconx(1,n))
               imem1 = imem1 + 1
               icount = 0
               jcount = 0
               memcur = memcur + 1
c
               do 140 ic = i1 , i2
                  if (ext1) jmax = ic
                  do 130 jc = j1 , jmax
                     ext5 = ext1 .and. ic.eq.jc
                     itprim = nc(ic)*nc(jc)
                     if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
                     jcount = jcount + 1
                     zmem(jcount,memcur) = 0.0d0
c
                     do 120 iprim = 1 , itprim
                        icount = icount + 1
                        zmem(jcount,memcur) = zmem(jcount,memcur)
     +                     + zmem(icount,imem1)
 120                 continue
 130              continue
 140           continue
 150        continue
 160     continue
c
c**** form contracted intermediates (ix-1 | jx) and (ix | jx-1)
c
         do 210 n = 1 , nconb
            imem1 = iconb(3,n)
            do 200 i = 1 , nvl(iconb(1,n))
               imem1 = imem1 + 1
               icount = 0
               jcount = 0
               memcur = memcur + 1
c
               do 190 ic = i1 , i2
                  if (ext1) jmax = ic
                  do 180 jc = j1 , jmax
                     ext5 = ext1 .and. ic.eq.jc
                     itprim = nc(ic)*nc(jc)
                     if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
                     jcount = jcount + 1
                     zmem(jcount,memcur) = 0.0d0
c
                     do 170 iprim = 1 , itprim
                        icount = icount + 1
                        zmem(jcount,memcur) = zmem(jcount,memcur)
     +                     + zmem(icount,imem1)
     +                     *(zmem(icount,7)-xijc(icount))
 170                 continue
 180              continue
 190           continue
 200        continue
 210     continue
         do 260 n = 1 , ncona
            imem1 = icona(3,n)
            do 250 i = 1 , nvl(icona(1,n))
               imem1 = imem1 + 1
               icount = 0
               jcount = 0
               memcur = memcur + 1
c
               do 240 ic = i1 , i2
                  if (ext1) jmax = ic
                  do 230 jc = j1 , jmax
                     ext5 = ext1 .and. ic.eq.jc
                     itprim = nc(ic)*nc(jc)
                     if (ext5) itprim = nc(ic)*(nc(ic)+1)/2
c
                     jcount = jcount + 1
                     zmem(jcount,memcur) = 0.0d0
c
                     do 220 iprim = 1 , itprim
                        icount = icount + 1
                        zmem(jcount,memcur) = zmem(jcount,memcur)
     +                     + zmem(icount,imem1)*xijc(icount)
 220                 continue
 230              continue
 240           continue
 250        continue
 260     continue
c
c****** form the gamma which will multiply the derivative
c****** overlap integrals
c
         mgamma = memcur
         do 310 i = 1 , ivlx
            do 300 j = 1 , jvlx
               icount = 0
               memcur = memcur + 1
               do 290 ic = i1 , i2
                  id = (ic-i1)*ivlx + i + ipos
                  if (ext1) jmax = ic
c
                  if (lx.eq.1) then
                     do 270 jd = j1 , jmax
                        icount = icount + 1
                        zmem(icount,memcur) = dens(jd,id)
 270                 continue
                  else
                     jd = j + jpos - jvlx
                     do 280 jc = j1 , jmax
                        jd = jd + jvlx
                        icount = icount + 1
                        zmem(icount,memcur) = dens(jd,id)
 280                 continue
                  end if
 290           continue
 300        continue
 310     continue
c
c**** form ai-bi in zmemc(1,4-6)
c
         jcount = 0
         do 330 ic = i1 , i2
            if (ext1) jmax = ic
            do 320 jc = j1 , jmax
               jcount = jcount + 1
               zmem(jcount,4) = cobas(ic,1) - cobas(jc,1)
               zmem(jcount,5) = cobas(ic,2) - cobas(jc,2)
               zmem(jcount,6) = cobas(ic,3) - cobas(jc,3)
 320        continue
 330     continue
c**** call the horizontal recursion relations and form gradient
c**** in zmemc(1,mtgrad + 1 to 6)
c
         mtgrad = 0
         mtemp = 3
         ihzc = 0
         call hzdrvd(ix+1,jx,ix-1,jx,ihzc,zmem,jcount,mgamma,mtemp,
     +               mtgrad)
         mtgrad = 3
         call hzdrvd(ix,jx+1,ix,jx-1,ihzc,zmem,jcount,mgamma,mtemp,
     +               mtgrad)
c
c**** now form the gradient
c
         jcount = 0
         do 360 ic = i1 , i2
            if (ext1) jmax = ic
            iatom = nwa(nup(ic))
            do 350 jc = j1 , jmax
               jcount = jcount + 1
               jatom = nwa(nup(jc))
               do 340 icart = 1 , 3
                  grad(iatom,icart) = grad(iatom,icart)
     +                                - zmem(jcount,icart)
                  grad(jatom,icart) = grad(jatom,icart)
     +                                - zmem(jcount,icart+3)
                  grtemp(natom,icart) = grtemp(natom,icart)
     +                                  + zmem(jcount,icart+3)
     +                                  + zmem(jcount,icart)
 340           continue
 350        continue
 360     continue
c
c      do 1 icart=1,3
c      write(iwr,*) 'cartesian component ',icart
c      do 1 i=1,na
c      write(iwr,*) grtemp(i,icart)
c1     continue
c
 370  continue
      return
      end
      subroutine pp1(xijk,xija,xijp1,xijp2,xijp3,i1,i2
     &,j1,j2,xt1,xt2,wx
     &,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ccc1,ccc2,ext1)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (p/p) and (p/h/p) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      dimension xija(n11),xijk(n11),ccc2(n11,3)
     &,s1(n11*9),h1(n11*9),sx(nb,nb),hx(nb,nb),ccc1(n11,3)
     &,lar(n11),f(7*n11),t(n11),xt1(n14*3,n13*3),xt2(n16*3,n15*3)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
     &,xijp1(n11),xijp2(n11),xijp3(n11)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11*9)
      icount = 0
c
      do 30 i = imin , imax
         if (ext1) jmax = i
         do 20 j = jmin , jmax
            icount = icount + 1
            ccc1(icount,1) = xijp1(icount) - coprm(i,1)
            ccc1(icount,2) = xijp2(icount) - coprm(i,2)
            ccc1(icount,3) = xijp3(icount) - coprm(i,3)
c
            ccc2(icount,1) = xijp1(icount) - coprm(j,1)
            ccc2(icount,2) = xijp2(icount) - coprm(j,2)
            ccc2(icount,3) = xijp3(icount) - coprm(j,3)
 20      continue
 30   continue
c
      do 60 k = 1 , na
c
         do 40 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 40      continue
c
         call fquik2(f,f,lar,t,3,n11)
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 50 icount = 1 , n11
            index1 = n11 + icount
            index2 = n11 + index1
            index3 = n11 + index2
            index4 = n11 + index3
            index5 = n11 + index4
            index6 = n11 + index5
            index7 = n11 + index6
            index8 = n11 + index7
            x3 = xija(icount) + xija(icount)
c
c     ***** xx1 has the k dependant factors of (pi/v/s)m=0 ******
c     ***** xx2 has the k dependant factors of (pi/v/s)m=1 ******
c
            x1 = xijp1(icount) - coord(1,k)
            x1x = ccc1(icount,1)*f(icount) - x1*f(index1)
            x2x = ccc1(icount,1)*f(index1) - x1*f(index2)
            x1 = xijp2(icount) - coord(2,k)
            x1y = ccc1(icount,2)*f(icount) - x1*f(index1)
            x2y = ccc1(icount,2)*f(index1) - x1*f(index2)
            x1 = xijp3(icount) - coord(3,k)
            x1z = ccc1(icount,3)*f(icount) - x1*f(index1)
            x2z = ccc1(icount,3)*f(index1) - x1*f(index2)
c
c     **** xxx(3,3) has contributions of pj & cj terms to (pi/v/pj)
c
            x4 = (f(icount)-f(index1))/x3
            x2 = xijp1(icount) - coord(1,k)
            h1(icount) = h1(icount) + (ccc2(icount,1)*x1x-x2*x2x+x4)
     +                   *charge(k)
            h1(index1) = h1(index1) + (ccc2(icount,1)*x1y-x2*x2y)
     +                   *charge(k)
            h1(index2) = h1(index2) + (ccc2(icount,1)*x1z-x2*x2z)
     +                   *charge(k)
            x2 = xijp2(icount) - coord(2,k)
            h1(index3) = h1(index3) + (ccc2(icount,2)*x1x-x2*x2x)
     +                   *charge(k)
            h1(index4) = h1(index4) + (ccc2(icount,2)*x1y-x2*x2y+x4)
     +                   *charge(k)
            h1(index5) = h1(index5) + (ccc2(icount,2)*x1z-x2*x2z)
     +                   *charge(k)
            x2 = xijp3(icount) - coord(3,k)
            h1(index6) = h1(index6) + (ccc2(icount,3)*x1x-x2*x2x)
     +                   *charge(k)
            h1(index7) = h1(index7) + (ccc2(icount,3)*x1y-x2*x2y)
     +                   *charge(k)
            h1(index8) = h1(index8) + (ccc2(icount,3)*x1z-x2*x2z+x4)
     +                   *charge(k)
c
 50      continue
 60   continue
c
      icount = 0
      do 80 i = imin , imax
         if (ext1) jmax = i
         do 70 j = jmin , jmax
            icount = icount + 1
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(icount) = pe(i)*pe(j)
 70      continue
 80   continue
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 90 icount = 1 , n11
         index1 = n11 + icount
         index2 = n11 + index1
         index3 = n11 + index2
         index4 = n11 + index3
         index5 = n11 + index4
         index6 = n11 + index5
         index7 = n11 + index6
         index8 = n11 + index7
c
c     ***** x3 contains (s/s) and xx1 has (pi/s) *******
c
         x1 = 1.0d0/xija(icount)
         x2 = 2.0d0*f(icount)*x1
         x3 = 0.5d0*x2*(7.0d0-x2*t(icount))
         x4 = xijk(icount)*pi14*dsqrt(x1)
         x6 = x4*x1*0.5d0
         x5 = x2*x6
         x7 = xijk(icount)/pi14
         x1x = ccc1(icount,1)*x4
         x1y = ccc1(icount,2)*x4
         x1z = ccc1(icount,3)*x4
c
         s1(icount) = ccc2(icount,1)*x1x + x6
         h1(icount) = x3*s1(icount) - x7*h1(icount) - x5
         s1(index1) = ccc2(icount,1)*x1y
         h1(index1) = x3*s1(index1) - x7*h1(index1)
         s1(index2) = ccc2(icount,1)*x1z
         h1(index2) = x3*s1(index2) - x7*h1(index2)
c
         s1(index3) = ccc2(icount,2)*x1x
         h1(index3) = x3*s1(index3) - x7*h1(index3)
         s1(index4) = ccc2(icount,2)*x1y + x6
         h1(index4) = x3*s1(index4) - x7*h1(index4) - x5
         s1(index5) = ccc2(icount,2)*x1z
         h1(index5) = x3*s1(index5) - x7*h1(index5)
c
         s1(index6) = ccc2(icount,3)*x1x
         h1(index6) = x3*s1(index6) - x7*h1(index6)
         s1(index7) = ccc2(icount,3)*x1y
         h1(index7) = x3*s1(index7) - x7*h1(index7)
         s1(index8) = ccc2(icount,3)*x1z + x6
         h1(index8) = x3*s1(index8) - x7*h1(index8) - x5
 90   continue
c
c****** put the result in the matrices to be transformed ***
c
      icount = 0
      do 130 jc = 1 , 3
         do 120 ic = 1 , 3
            jmax = n16
            do 110 i = 1 , n14
               ip = (i-1)*3
               if (ext1) jmax = i
               do 100 j = 1 , jmax
                  jp = (j-1)*3
                  icount = icount + 1
                  wx(ip+ic,jp+jc) = s1(icount)
                  wy(ip+ic,jp+jc) = h1(icount)
 100           continue
 110        continue
 120     continue
 130  continue
c
      if (ext1) then
         do 150 i = 1 , n14*3
            do 140 j = 1 , i
               wx(j,i) = wx(i,j)
               wy(j,i) = wy(i,j)
 140        continue
 150     continue
      end if
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16*3,wz,1,ipar2,n14*3,n16*3,n15*3)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*3,wx,ipar2,1,n15*3,n14*3,n13*3)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16*3,wz,1,ipar2,n14*3,n16*3,n15*3)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*3,wy,ipar2,1,n15*3,n14*3,n13*3)
c
      ix1 = 3*i1 - ns*2 - 2
      ix2 = ix1 + n13*3 - 1
      jx1 = 3*j1 - ns*2 - 2
      jx2 = jx1 + n15*3 - 1
      ic = 0
      do 170 i = ix1 , ix2
         ic = ic + 1
         jc = 0
         if (ext1) jx2 = i
         do 160 j = jx1 , jx2
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 160     continue
 170  continue
c
      return
      end
      subroutine ppex(xija,xijk,xijp1,xijp2,xijp3,xijt
     &,i1,imax0,j1,jmax0,gx,gy,icmax,ext1)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c     *****************************************************
c     ******* evaluating exchange integrals of type (pp/pp) ********
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      logical ext1
      dimension xija(icmax),xijk(icmax),xijt(icmax)
     &,gx(icmax,3),gy(icmax,3)
     &,xijp1(icmax),xijp2(icmax),xijp3(icmax)
c
c******** n6,n8,n13 and n15 are number of basis functions of *******
c******** k,l,i and j being dealt with in this batch ***************
c******** n7,n9,n14 and n16 are number of primitives of ************
c******** k,l,i and j being dealt with in this batch ***************
c
      jmax = jmax0
      ipos = nup(i1) - nc(i1)
      jpos = nup(j1) - nc(j1)
      jcdash = 0
      do 30 k = 1 , imax0
         kk = k + ipos
         if (ext1) jmax = k
         do 20 l = 1 , jmax
            ll = l + jpos
            jcdash = jcdash + 1
            gy(jcdash,1) = coprm(kk,1)
            gy(jcdash,2) = coprm(kk,2)
            gy(jcdash,3) = coprm(kk,3)
            gx(jcdash,1) = coprm(ll,1)
            gx(jcdash,2) = coprm(ll,2)
            gx(jcdash,3) = coprm(ll,3)
 20      continue
 30   continue
c
      third = 1.0d0/3.0d0
      fifth = 1.0d0/5.0d0
      do 40 icount = 1 , icmax
         z2xi = xijp1(icount) - gy(icount,1)
         z2xj = xijp1(icount) - gx(icount,1)
         z2yi = xijp2(icount) - gy(icount,2)
         z2yj = xijp2(icount) - gx(icount,2)
         z2zi = xijp3(icount) - gy(icount,3)
         z2zj = xijp3(icount) - gx(icount,3)
c
         x1 = 1.0d0/(xija(icount)+xija(icount))
         x2 = xijk(icount)*xijk(icount)*dsqrt(x1)
         x19 = xija(icount)*xija(icount)*x1
         x17 = x19/xija(icount)
         x18 = x19/xija(icount)
         x110 = 0.5d0/xija(icount)
         x11 = 0.5d0/xija(icount)
         x12 = (1.0d0-x17*third)*x11
         x13 = (third-x17*fifth)*x11
         x14 = 0.5d0*x1
         x16 = fifth*x14
c
c     ***** x3 contains (ps/ss)m=0 ; x4 has (ps/ss)m=1; *****
c     ***** x5 has (ps/ss)m=2 ; x6 holds (ps/ss)m=3 *********
c     ***** in x7 is (sp/ss)m=0 ; x8 holds (sp/ss)m=1 *******
c     ***** x9 has (sp/ss)m=2. (not including x2 which is ***
c     ***** common to all terms in recursion relations) *****
c
c        x71 = z2xj
         x81 = third*z2xj
c        x91 = fifth*z2xj
         x31 = z2xi
         x41 = third*z2xi
c        x51 = fifth*z2xi
c
c        x72 = z2yj
         x82 = third*z2yj
c        x92 = fifth*z2yj
         x32 = z2yi
         x42 = third*z2yi
c        x52 = fifth*z2yi
c
c        x73 = z2zj
         x83 = third*z2zj
c        x93 = fifth*z2zj
         x33 = z2zi
         x43 = third*z2zi
c        x53 = fifth*z2zi
c
c
c     *** xx1 contains (pp/ss)m=0 ; xx2 holds (pp/ss)m=1 ***
c     *** xx2 has (pp/ss)m=2 ; xx4 contains (ps/ps)m=1 *****
c     ***********   in xx5 there is (sp/ps)m=1  ************
c
c
         xx111 = x12 + x31*z2xj
         xx121 = x32*z2xj
         xx131 = x33*z2xj
         xx211 = x13 + x41*z2xj
c        xx221 = x42*z2xj
c        xx231 = x43*z2xj
c
         xx112 = x31*z2yj
         xx122 = x12 + x32*z2yj
         xx132 = x33*z2yj
c        xx212 = x41*z2yj
         xx222 = x13 + x42*z2yj
c        xx232 = x43*z2yj
c
         xx113 = x31*z2zj
         xx123 = x32*z2zj
         xx133 = x12 + x33*z2zj
c        xx213 = x41*z2zj
c        xx223 = x42*z2zj
         xx233 = x13 + x43*z2zj
c
         xx411 = x16 + x41*z2xi
c        xx421 = x42*z2xi
c        xx431 = x43*z2xi
         xx511 = x16 + x81*z2xi
c        xx521 = x82*z2xi
c        xx531 = x83*z2xi
c
c        xx412 = x41*z2yi
         xx422 = x16 + x42*z2yi
c        xx432 = x43*z2yi
c        xx512 = x81*z2yi
         xx522 = x16 + x82*z2yi
c        xx532 = x83*z2yi
c
c        xx413 = x41*z2zi
c        xx423 = x42*z2zi
         xx433 = x16 + x43*z2zi
c        xx513 = x81*z2zi
c        xx523 = x82*z2zi
         xx533 = x16 + x83*z2zi
c
c     ********    in xz1 we have (pp/ps)m=0      *********
c
         xz1111 = (x81+x41)*x14 + xx111*z2xi
c        xz1211 = x42*x14 + xx121*z2xi
c        xz1311 = x43*x14 + xx131*z2xi
         xz1121 = x82*x14 + xx112*z2xi
c        xz1221 = xx122*z2xi
c        xz1321 = xx132*z2xi
         xz1131 = x83*x14 + xx113*z2xi
c        xz1231 = xx123*z2xi
c        xz1331 = xx133*z2xi
c
c        xz1112 = xx111*z2yi
         xz1212 = x81*x14 + xx121*z2yi
c        xz1312 = xx131*z2yi
c        xz1122 = x41*x14 + xx112*z2yi
         xz1222 = (x82+x42)*x14 + xx122*z2yi
c        xz1322 = x43*x14 + xx132*z2yi
c        xz1132 = xx113*z2yi
         xz1232 = x83*x14 + xx123*z2yi
c        xz1332 = xx133*z2yi
c
c        xz1113 = xx111*z2zi
c        xz1213 = xx121*z2zi
         xz1313 = x81*x14 + xx131*z2zi
c        xz1123 = xx112*z2zi
c        xz1223 = xx122*z2zi
         xz1323 = x82*x14 + xx132*z2zi
c        xz1133 = x41*x14 + xx113*z2zi
c        xz1233 = x42*x14 + xx123*z2zi
         xz1333 = (x43+x83)*x14 + xx133*z2zi
c
c     ********    find the integrals now     **********
c     ******** kronecker deleta terms occur when **********
c     ******** lp = ip,jp or kp.    *********************
c     ******** ip =lp then factor is ((sp/ps)m=1)/x14 ********
c     ******** jp =lp then factor is ((ps/ps)m=1)/x14 ********
c     ******** kp = lp then  x110(((pp/ss)m=0)-x18((pp/ss)m=1)) *****
c
         gg1 = xz1111*z2xj + (xx511+xx411)*x14 + x110*(xx111-x18*xx211)
         gg2 = xz1212*z2xj + xx422*x14
         gg3 = xz1313*z2xj + xx433*x14
c
         gg4 = xz1121*z2yj + xx411*x14
         gg5 = xz1222*z2yj + (xx522+xx422)*x14 + x110*(xx122-x18*xx222)
         gg6 = xz1323*z2yj + xx433*x14
c
         gg7 = xz1131*z2zj + xx411*x14
         gg8 = xz1232*z2zj + xx422*x14
         gg9 = xz1333*z2zj + (xx433+xx533)*x14 + x110*(xx133-x18*xx233)
c
         x1a = max(gg1,gg2,gg3,gg4,gg5,gg6,gg7,gg8,gg9)
         xijt(icount) = xijt(icount)*dsqrt(x1a*x2)
 40   continue
c
      return
      end
      subroutine ps1(xijk,xija,xijp1,xijp2,xijp3,i1,i2,j1
     &,j2,xt1,xt2,wx,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ccc1)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (p/s) and (p/h/s) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16
      dimension xija(n11),xijk(n11)
     &,s1(n11*3),h1(n11*3),sx(nb,nb),hx(nb,nb),ccc1(n11,3)
     &,lar(n11),f(6*n11),t(n11),xt1(n14*3,n13*3),xt2(n16,n15)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
     &,xijp1(n11),xijp2(n11),xijp3(n11)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
c     ******* (p/s) and (p/h(1)/s) integrals ********
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11*3)
      icount = 0
      do 30 i = imin , imax
         do 20 j = jmin , jmax
            icount = icount + 1
            ccc1(icount,1) = xijp1(icount) - coprm(i,1)
            ccc1(icount,2) = xijp2(icount) - coprm(i,2)
            ccc1(icount,3) = xijp3(icount) - coprm(i,3)
 20      continue
 30   continue
c
      do 60 k = 1 , na
c
         do 40 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 40      continue
         call fquik2(f,f,lar,t,2,n11)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 50 icount = 1 , n11
            index1 = n11 + icount
            index2 = n11 + index1
            h1(icount) = h1(icount)
     +                   + (ccc1(icount,1)*f(icount)-(xijp1(icount)
     +                   -coord(1,k))*f(index1))*charge(k)
            h1(index1) = h1(index1)
     +                   + (ccc1(icount,2)*f(icount)-(xijp2(icount)
     +                   -coord(2,k))*f(index1))*charge(k)
            h1(index2) = h1(index2)
     +                   + (ccc1(icount,3)*f(icount)-(xijp3(icount)
     +                   -coord(3,k))*f(index1))*charge(k)
 50      continue
 60   continue
c
      icount = 0
      do 80 i = imin , imax
         do 70 j = jmin , jmax
            icount = icount + 1
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(icount) = pe(i)*pe(j)
 70      continue
 80   continue
c
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 90 icount = 1 , n11
         index1 = n11 + icount
         index2 = n11 + index1
c
         x1 = 1.0d0/xija(icount)
         x2 = f(icount)*x1
         x3 = xijk(icount)*pi14*dsqrt(x1)
         x4 = x2*(5.0d0-2.0d0*t(icount)*x2)
         x5 = xijk(icount)/pi14
         s1(icount) = ccc1(icount,1)*x3
         h1(icount) = s1(icount)*x4 - h1(icount)*x5
         s1(index1) = ccc1(icount,2)*x3
         h1(index1) = s1(index1)*x4 - h1(index1)*x5
         s1(index2) = ccc1(icount,3)*x3
         h1(index2) = s1(index2)*x4 - h1(index2)*x5
 90   continue
c
c****** put the result in the matrices to be transformed ***
c
      icount = 0
      do 120 ic = 1 , 3
         do 110 i = 1 , n14
            ip = (i-1)*3
            do 100 j = 1 , n16
               icount = icount + 1
               wx(ip+ic,j) = s1(icount)
               wy(ip+ic,j) = h1(icount)
 100        continue
 110     continue
 120  continue
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16,wz,1,ipar2,n14*3,n16,n15)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*3,wx,ipar2,1,n15,n14*3,n13*3)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16,wz,1,ipar2,n14*3,n16,n15)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14*3,wy,ipar2,1,n15,n14*3,n13*3)
c
      ix1 = 3*i1 - ns*2 - 2
      ix2 = ix1 + n13*3 - 1
      ic = 0
      do 140 i = ix1 , ix2
         ic = ic + 1
         jc = 0
         do 130 j = j1 , j2
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 130     continue
 140  continue
c
      return
      end
      subroutine psex(xija,xijk,xijp1,xijp2,xijp3,xijt,gg
     &,icmax,i1,imax0,jmax0)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c     *****************************************************
c     ******* evaluating exchange integrals of type (ps/ps) ********
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      dimension xija(icmax),xijk(icmax),xijt(icmax)
     &,gg(icmax,3),xijp1(icmax),xijp2(icmax),xijp3(icmax)
c
      ipos = nup(i1) - nc(i1)
      jcdash = 0
      do 30 i = 1 , imax0
         kk = i + ipos
         do 20 j = 1 , jmax0
            jcdash = jcdash + 1
            gg(jcdash,1) = coprm(kk,1)
            gg(jcdash,2) = coprm(kk,2)
            gg(jcdash,3) = coprm(kk,3)
 20      continue
 30   continue
c
c************* the integral evaluation ******************
c
      sixth = 1.0d0/6.0d0
      do 40 icount = 1 , icmax
         x1 = 1.0d0/(xija(icount)+xija(icount))
         x12 = sixth*x1
c
         z2xi = xijp1(icount) - gg(icount,1)
         z2xi = z2xi*z2xi
         z2yi = xijp2(icount) - gg(icount,2)
         z2yi = z2yi*z2yi
         z2zi = xijp3(icount) - gg(icount,3)
         z2zi = z2zi*z2zi
         x1a = max(z2xi,z2yi,z2zi) + x12
         xijt(icount) = xijt(icount)*dabs(xijk(icount))
     +                  *dsqrt(x1a*dsqrt(x1))
c
 40   continue
      return
      end
      subroutine putqd (q,shalf,p,scra,eig,pop,energy,nb,nprint)
c
c     store direct-scf vectors, density, eigen values and
c     populations in gamess format, must issue
c     call to tback to generate mos in sabf
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/prints)
INCLUDE(common/atmol3)
INCLUDE(common/wrtd)
INCLUDE(common/tran)
      common/gtrans/idtran(maxorb)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/mapper)
      common/scfopt/maxit(4),accdi(2),icoupl(4),dmpcut(3),etot
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nbas
     &,ns,np,nsp,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,
     & repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
INCLUDE(common/runlab)
      dimension q(nb,nb),eig(*),scra(*),pop(*),p(*),shalf(nb,nb)
      data m1/1/
c
      etot = energy
      nocc = ne/2
c     l3 = nb*nb
      l2 = ikyp(nb)
      call rdedx(q,nb*nb,idblk(5+inddir),numdir)
      call scfit(q,shalf,p,scra,eig)
      call wrt3(p,nb*nb,idblk(26),numdir)
      call vclr(pop,1,nb)
      do 40 loop = 1 , nocc
         pop(loop) = 2.0d0
 40   continue
      call putq(zcom,ztitle,eig,pop,nb,nb,nb,m1,m1,shalf,mouta,iblkv)
      call wrt3(eig,nb,ibl3ea,idaf)
      if (nprint.ne.5 .and. nprint.ne.-5)
     +    call analmo(shalf,eig,pop,ilifq,nb,nb)
      call tdownd(shalf,ilifq,shalf,ilifq,nb)
      do 50 i = 1 , nb
      do 60 j = 1 , nb
       jd = idtran(j)
       q(j,i) = shalf(jd,i)
 60   continue
 50   continue
      call dmtx(p,q,pop,iky,nocc,nb,nb)
      call wrt3(p,l2,ibl3pa,idaf)
      if (nprint.eq.5 .or. nprint.eq.-5) return
      lprnt = nb
      if (.not.oprint(20)) lprnt = min(nocc+5,nb)
      write (iwr,6010)
      call prev(q,eig,lprnt,nb,nb)
      return
 6010 format (/1x,100('-')//50x,12('-')/50x,'eigenvectors'/50x,12('-'))
      end
      subroutine recur1(ans,pai,xint1,wpi,xint2,num)
c
c***** the recur set of routines constitute the inner kernel of the
c***** integral code. they are the vrr of os and the hrr of hgp
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),wpi(num)
c
c***** this routince uses the first 2 terms of recursion relation
c***** building up on i, to generate new integrals from old
c
      do 20 i = 1 , num
         ans(i) = pai*xint1(i) + wpi(i)*xint2(i)
 20   continue
      return
      end
      subroutine recur2(ans,pai,xint1,wpi,xint2,xija
     &,nai,xint3,xkla,rho,xint4,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),wpi(num),xkla(num),rho(num)
c
c***** this routince uses first 4 terms of recursion relation
c***** building up on i, to generate new integrals from old
c
      xmult2 = 0.5d0*dfloat(nai)/xija
c
      do 20 i = 1 , num
         ans(i) = pai*xint1(i) + wpi(i)*xint2(i)
     +            + xmult2*(xint3(i)-rho(i)*xkla(i)*xint4(i))
 20   continue
      return
      end
      subroutine recur3(ans,pai,xint1,wpi,xint2
     &,rho,nci,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num)
     &,xint5(num),wpi(num),rho(num)
c
c***** this routince uses terms 1,2 and 5 of recursion relation
c***** building up on i, to generate new integrals from old
c
      xmult1 = 0.5d0*dfloat(nci)
c
      do 20 i = 1 , num
         ans(i) = pai*xint1(i) + wpi(i)*xint2(i) + xmult1*rho(i)
     +            *xint5(i)
 20   continue
      return
      end
      subroutine recur4(ans,pai,xint1,wpi,xint2,xija
     &,nai,xint3,xkla,rho,xint4,nci,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),xint5(num),wpi(num),xkla(num),rho(num)
c
c***** this routince uses the full recursion relation building
c***** up on i, to generate new integrals from old
c
      xmult1 = 0.5d0*dfloat(nci)
      xmult2 = 0.5d0*dfloat(nai)/xija
c
      do 20 i = 1 , num
         ans(i) = pai*xint1(i) + wpi(i)*xint2(i)
     +            + xmult2*(xint3(i)-rho(i)*xkla(i)*xint4(i))
     +            + xmult1*rho(i)*xint5(i)
 20   continue
      return
      end
      subroutine recur5(ans,qci,xint1,wqi,xint2,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),wqi(num)
     &,qci(num)
c
c***** this routine uses terms 1,2 of recursion relation building
c***** up on k, to generate new integrals from old
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + wqi(i)*xint2(i)
 20   continue
      return
      end
      subroutine recur6(ans,qci,xint1,wqi,xint2,xkla
     &,nci,xint3,xija,rho,xint4,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),wqi(num),xkla(num),rho(num)
     &,qci(num)
c
c***** this routine uses terms 1-4 of recursion relation building
c***** up on k, to generate new integrals from old
c
      xmult1 = 0.5d0*dfloat(nci)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + wqi(i)*xint2(i) + xmult1/xkla(i)
     +            *(xint3(i)-rho(i)*xija*xint4(i))
 20   continue
      return
      end
      subroutine recur7(ans,qci,xint1,wqi,xint2
     &,rho,nai,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num)
     &,xint5(num),wqi(num),rho(num)
     &,qci(num)
c
c***** this routine uses terms 1,2,5 of recursion relation
c***** building up on k, to generate new integrals from old
c
      xmult2 = 0.5d0*dfloat(nai)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + wqi(i)*xint2(i) + xmult2*rho(i)
     +            *xint5(i)
 20   continue
      return
      end
      subroutine recur8(ans,qci,xint1,wqi,xint2,xkla
     &,nci,xint3,xija,rho,xint4,nai,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint2(num),xint3(num)
     &,xint4(num),xint5(num),wqi(num),xkla(num),rho(num)
     &,qci(num)
c
c***** this routine uses the full recursion relation building
c***** up on k, to generate new integrals from old
c
      xmult1 = 0.5d0*dfloat(nci)
      xmult2 = 0.5d0*dfloat(nai)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + wqi(i)*xint2(i) + xmult1/xkla(i)
     +            *(xint3(i)-rho(i)*xija*xint4(i)) + xmult2*rho(i)
     +            *xint5(i)
 20   continue
      return
      end
      subroutine reducd(q,num,nred,v,s,qnew,ichar,iorb,uhf,iwr)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension q(num,nred),v(num,nred),qnew(num,nred),
     &  s(*),ichar(*),vlen(8),irrc(maxorb),iperm(maxorb)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt,nt3,itable(8,8),
     1mult(8,8),irr(maxorb),imosb(8,maxorb),irrb(maxorb)
      common/vectrn/isecv
      common/restri/nfile(63),lds(508),isect(508),ldx(508)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
      logical uhf
      logical easy(maxorb)
c
c     uses projection operators to reduce set of vectors
c     in q(num,nred) to nred vectors of character in
c     ichar
c
c     first deal with all those orbitals in the set whose
c     symmetry is easy to determine
      call vclr(qnew,1,num*nred)
      nproj = 0
      do 20 i = 1 , nred
         irrc(i) = 0
         easy(i) = .false.
 20   continue
      nredd = nred
c
c------------------------easy symmetries-----------------
c
      do 60 ired = 1 , nred
c
c     decompose one orbital into symmetry types
c
         call indsyt(num,q(1,ired),v,s,vlen)
c
         do 50 irep = 1 , nt
            if (dabs(vlen(irep)).gt.0.90d0) then
c     orbital is mostly one symmetry type - seems safe!
               nproj = nproj + 1
               nredd = nredd - 1
               ichar(irep) = ichar(irep) - 1
c     write(iwr,*) 'easy symmetry, nproj, irep', nproj,irep
               if (ichar(irep).lt.0) then
                  write (iwr,*)
     +                         ' something wrong with symmetry (case 1)'
                  call caserr(' give up ')
               end if
               easy(ired) = .true.
               irrc(nproj) = irep
               iperm(nproj) = ired
               do 30 k = 1 , num
                  qnew(k,nproj) = v(k,irep)
 30            continue
               do 40 k = 1 , num
                  q(k,ired) = q(k,ired) - v(k,irep)
 40            continue
               go to 60
            end if
 50      continue
 60   continue
      do 80 i = 2 , nproj
         do 70 j = 1 , i - 1
            call vschmv(num,qnew(1,j),qnew(1,i),s)
 70      continue
 80   continue
      do 90 i = 1 , nproj
         over = vsv(num,qnew(1,i),qnew(1,i),s)
         if (over.le.1.0d-6) then
            write (iwr,*)
     +                   'very short vector produced end of first stage'
            write (iwr,*) 'nred, nredd, nproj, i ' , nred , nredd ,
     +                    nproj , i
            write (iwr,*) 'overlap ' , over
            call caserr('give up')
         end if
         call vrenrm(num,qnew(1,nproj),s)
 90   continue
c------------------------------------------------------------------
c
c     are nredd vectors remaining whose symmetry has not been
c     assigned
 100  j = 0
      do 110 i = 1 , nt
         j = j + ichar(i)
 110  continue
      if (j.ne.nredd) then
         write (iwr,*) ' inconsistency between number of orbitals'
         write (iwr,*) ' and character of representation'
         call caserr('give up')
      end if
      if (nredd.eq.0) then
c
c------------------------------finished ?------------------
c
c
         do 130 k = 2 , nproj
            do 120 l = 1 , k - 1
               call vschmv(num,qnew(1,l),qnew(1,k),s)
 120        continue
 130     continue
         do 140 i = 1 , nproj
            over = vsv(num,qnew(1,i),qnew(1,i),s)
            if (over.le.1.0d-6) then
               write (iwr,*) 'very short vector produced at final stage'
               write (iwr,*) 'nred, nredd, nproj, i ' , nred , nredd ,
     +                       nproj , i
               write (iwr,*) 'overlap ' , over
               call caserr('give up')
            end if
            call vrenrm(num,qnew(1,i),s)
 140     continue
c
c-----------------------------------------------------------
c
         if (uhf .and. isecv.eq.isect(11)) then
            do 170 i = 1 , nred
               irm = iperm(i)
               irrb(iorb+irm-1) = irrc(i)
               do 150 k = 1 , num
                  q(k,irm) = qnew(k,i)
 150           continue
               do 160 k = 1 , nt
                  imosb(k,iorb+irm-1) = itable(irrc(i),k)
 160           continue
 170        continue
         else
            do 200 i = 1 , nred
               irm = iperm(i)
               irr(iorb+irm-1) = irrc(i)
               do 180 k = 1 , num
                  q(k,irm) = qnew(k,i)
 180           continue
               do 190 k = 1 , nt
                  imos(k,iorb+irm-1) = itable(irrc(i),k)
 190           continue
 200        continue
         end if
c
c
c     have completely reduced representation
c     all orbitals have been assigned symmetry type, once and
c     once only, and have been left in an order as close as possible
c     to the original order
c     the orbitals are orthonormal and the projection/orthogonalisation
c     has not produced either a null vector nor the same vector twice
c
         return
      else
c
c   orthogonalise those vectors remaining to those
c   projected out already - this should prevent same
c   vector being projected out twice?
         do 220 i = 1 , nred
            if (.not.easy(i)) then
               do 210 j = 1 , nproj
                  call vschmv(num,qnew(1,j),q(1,i),s)
 210           continue
            end if
 220     continue
         do 240 i = 2 , nred
            if (.not.(easy(i))) then
               do 230 j = 1 , i - 1
                  if (.not.easy(j)) call vschmv(num,q(1,j),q(1,i),s)
 230           continue
            end if
 240     continue
c----------------------------------------------------------------
         do 290 irep = 1 , nt
            if (ichar(irep).gt.0) then
c
               call vclr(v,1,num)
               do 260 j = 1 , nred
                  if (.not.easy(j)) then
                     do 250 i = 1 , num
                        v(i,1) = v(i,1) + q(i,j)
 250                 continue
                  end if
 260           continue
               over = vsv(num,v(1,1),v(1,1),s)
               if (over.le.1.0d-6) then
                  write (iwr,*) 'very short vector prior to projection'
                  write (iwr,*) 'overlap ' , over
                  call caserr('give up')
               end if
c
c     scheme : start with sum of remaining
c              vectors
c            : project out a vector with a symmetry which is
c              known to be present
c            : loop over all symmetry types
c            : repeat if necessary till have enough vectors
c            : overall algorithm amounts to decomposition of
c              original set into vectors along all symmetry
c              directions, then resummation of these components
c              to give new set of vectors
c
               ichar(irep) = ichar(irep) - 1
               nredd = nredd - 1
               nproj = nproj + 1
               irrc(nproj) = irep
               call indprj(v(1,1),num,qnew(1,nproj),s,irep)
c     write(iwr,*)'difficult symmetry, nproj, irep',nproj,irep
c
c     take overlap of projected vector
c     with original orbitals to see which it resembles
c     most
c
               call vrenrm(num,qnew(1,nproj),s)
               do 270 j = 1 , nred
                  if (easy(j)) then
                     vlen(j) = 0.0d0
                  else
                     vlen(j) = vsv(num,q(1,j),qnew(1,nproj),s)
                  end if
 270           continue
               amm = 0.1d0
               irm = 0
               do 280 j = 1 , nred
                  if (dabs(vlen(j)).gt.amm) then
                     amm = dabs(vlen(j))
                     irm = j
                  end if
 280           continue
               if (irm.eq.0) then
c     reason for starting with amm=0.1 is to catch case
c     where projected orbital in no way resembles any of
c     original set
                  write (iwr,*)
     +                      'projected orbtial does not match original?'
                  call caserr(' give up ')
               end if
               call vschmv(num,qnew(1,nproj),q(1,irm),s)
               iperm(nproj) = irm
               easy(irm) = .true.
            end if
 290     continue
         go to 100
      end if
      end
      subroutine rex5(ans,qci,xint1,num)
c
c***** the rex  recursion relations calculate exchange integrals
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),qci(num)
c
c***** this routine uses terms 1 of recursion relation building
c***** up on k or i, to generate new exchange integrals from old
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i)
 20   continue
      return
      end
      subroutine rex6(ans,qci,xint1,nci,xint3,rho,xint4,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint3(num)
     &,xint4(num),rho(num)
     &,qci(num)
c
c***** this routine uses terms 1-4 of recursion relation building
c***** up on k or i, to generate new exchange integrals from old
c
      xmult1 = dfloat(nci)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + xmult1*rho(i)
     +            *(xint3(i)-0.5d0*xint4(i))
 20   continue
      return
      end
      subroutine rex7(ans,qci,xint1,rho,nai,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num)
     &,xint5(num),rho(num)
     &,qci(num)
c
c***** this routine uses terms 1,5 of recursion relation building
c***** up on k or i, to generate new exchange integrals from old
c
      xmult2 = 0.5d0*dfloat(nai)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + xmult2*rho(i)*xint5(i)
 20   continue
      return
      end
      subroutine rex8(ans,qci,xint1
     &,nci,xint3,rho,xint4,nai,xint5,num)
c
      implicit REAL  (a-h,o-z)
      dimension ans(num),xint1(num),xint3(num)
     &,xint4(num),xint5(num),rho(num)
     &,qci(num)
c
c***** this routine uses the full recursion relation building
c***** up on k or i, to generate new exchange integrals from old
c
      xmult1 = dfloat(nci)
      xmult2 = 0.5d0*dfloat(nai)
c
      do 20 i = 1 , num
         ans(i) = qci(i)*xint1(i) + xmult1*rho(i)
     +            *(xint3(i)-0.5d0*xint4(i)) + xmult2*rho(i)*xint5(i)
 20   continue
      return
      end
      subroutine rhrd
      implicit REAL  (a-h,o-z)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480)
      common/hsym/t(20,20),mink,maxk,lkt,minl,maxl,llt,ntr
      dimension u(10)
      equivalence (u(1),u1),(u(2),u2),(u(3),u3),(u(4),u4),
     &    (u(5),u5),(u(6),u6)
c     ----- right multiply  t  by  r,
c           result back in  t
      go to (110,90,70,20) , llt
c     ----- f shell
 20   nf = (ntr-1)*10
      do 60 k = mink , maxk
         do 30 l = 1 , 10
            u(l) = t(k,l+10)
 30      continue
         do 50 l = 1 , 10
            sum = 0.0d0
            do 40 i = 1 , 10
               sum = sum + u(i)*ftr(i,nf+l)
 40         continue
            t(k,l+10) = sum
 50      continue
 60   continue
      go to 110
c     ----- d shell
 70   nd = 6*ntr - 10
      do 80 k = mink , maxk
         u1 = t(k,5)
         u2 = t(k,6)
         u3 = t(k,7)
         u4 = t(k,8)
         u5 = t(k,9)
         u6 = t(k,10)
         t(k,5) = u1*dtr(1,nd+5) + u2*dtr(2,nd+5) + u3*dtr(3,nd+5)
     +            + u4*dtr(4,nd+5) + u5*dtr(5,nd+5) + u6*dtr(6,nd+5)
         t(k,6) = u1*dtr(1,nd+6) + u2*dtr(2,nd+6) + u3*dtr(3,nd+6)
     +            + u4*dtr(4,nd+6) + u5*dtr(5,nd+6) + u6*dtr(6,nd+6)
         t(k,7) = u1*dtr(1,nd+7) + u2*dtr(2,nd+7) + u3*dtr(3,nd+7)
     +            + u4*dtr(4,nd+7) + u5*dtr(5,nd+7) + u6*dtr(6,nd+7)
         t(k,8) = u1*dtr(1,nd+8) + u2*dtr(2,nd+8) + u3*dtr(3,nd+8)
     +            + u4*dtr(4,nd+8) + u5*dtr(5,nd+8) + u6*dtr(6,nd+8)
         t(k,9) = u1*dtr(1,nd+9) + u2*dtr(2,nd+9) + u3*dtr(3,nd+9)
     +            + u4*dtr(4,nd+9) + u5*dtr(5,nd+9) + u6*dtr(6,nd+9)
         t(k,10) = u1*dtr(1,nd+10) + u2*dtr(2,nd+10) + u3*dtr(3,nd+10)
     +             + u4*dtr(4,nd+10) + u5*dtr(5,nd+10) + u6*dtr(6,nd+10)
 80   continue
      go to 110
c     ----- p shell
 90   np = 3*ntr - 4
      do 100 k = mink , maxk
         u1 = t(k,2)
         u2 = t(k,3)
         u3 = t(k,4)
         t(k,2) = u1*ptr(1,np+2) + u2*ptr(2,np+2) + u3*ptr(3,np+2)
         t(k,3) = u1*ptr(1,np+3) + u2*ptr(2,np+3) + u3*ptr(3,np+3)
         t(k,4) = u1*ptr(1,np+4) + u2*ptr(2,np+4) + u3*ptr(3,np+4)
 100  continue
c     ----- left multiply  t  by r
c           result back in  t
 110  go to (210,190,170,120) , lkt
c     ------ f shell
 120  nf = (ntr-1)*10
      do 160 l = minl , maxl
         do 130 k = 1 , 10
            u(k) = t(k+10,l)
 130     continue
         do 150 k = 1 , 10
            sum = 0.0d0
            do 140 i = 1 , 10
               sum = sum + u(i)*ftr(i,nf+k)
 140        continue
            t(k+10,l) = sum
 150     continue
 160  continue
      go to 210
c     ----- d shell
 170  nd = 6*ntr - 10
      do 180 k = minl , maxl
         u1 = t(5,k)
         u2 = t(6,k)
         u3 = t(7,k)
         u4 = t(8,k)
         u5 = t(9,k)
         u6 = t(10,k)
         t(5,k) = u1*dtr(1,nd+5) + u2*dtr(2,nd+5) + u3*dtr(3,nd+5)
     +            + u4*dtr(4,nd+5) + u5*dtr(5,nd+5) + u6*dtr(6,nd+5)
         t(6,k) = u1*dtr(1,nd+6) + u2*dtr(2,nd+6) + u3*dtr(3,nd+6)
     +            + u4*dtr(4,nd+6) + u5*dtr(5,nd+6) + u6*dtr(6,nd+6)
         t(7,k) = u1*dtr(1,nd+7) + u2*dtr(2,nd+7) + u3*dtr(3,nd+7)
     +            + u4*dtr(4,nd+7) + u5*dtr(5,nd+7) + u6*dtr(6,nd+7)
         t(8,k) = u1*dtr(1,nd+8) + u2*dtr(2,nd+8) + u3*dtr(3,nd+8)
     +            + u4*dtr(4,nd+8) + u5*dtr(5,nd+8) + u6*dtr(6,nd+8)
         t(9,k) = u1*dtr(1,nd+9) + u2*dtr(2,nd+9) + u3*dtr(3,nd+9)
     +            + u4*dtr(4,nd+9) + u5*dtr(5,nd+9) + u6*dtr(6,nd+9)
         t(10,k) = u1*dtr(1,nd+10) + u2*dtr(2,nd+10) + u3*dtr(3,nd+10)
     +             + u4*dtr(4,nd+10) + u5*dtr(5,nd+10) + u6*dtr(6,nd+10)
 180  continue
      go to 210
c     ----- p shell
 190  np = 3*ntr - 4
      do 200 k = minl , maxl
         u1 = t(2,k)
         u2 = t(3,k)
         u3 = t(4,k)
         t(2,k) = u1*ptr(1,np+2) + u2*ptr(2,np+2) + u3*ptr(3,np+2)
         t(3,k) = u1*ptr(1,np+3) + u2*ptr(2,np+3) + u3*ptr(3,np+3)
         t(4,k) = u1*ptr(1,np+4) + u2*ptr(2,np+4) + u3*ptr(3,np+4)
 200  continue
 210  return
      end
      subroutine sav6(iscat1,xxk,xxa,xxt,xxp1,xxp2,xxp3,temp1
     &,wksp,ipar,icou,i1,i2,j1,j2,ext1,ext2,ext3,buff,ifreem)
c
c***** parameters stored for use in the two electron integral routines
c***** includes the coordination of the calculation of exchange integral
c
      implicit REAL  (a-h,o-z)
      character *2  ppx,psx,ssx,dsx,dpx,ddx,buff
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext2,gencon,ext4,ext3
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/iofile/iread,iwr
      logical glog,gdbg
      common/gen/maxit(2),concrit(5),glog(5),gdbg(10)
      dimension xxk(ipar),xxa(ipar),xxt(ipar),iscat1(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar)
     &,temp1(ipar*3),wksp(ipar*3)
c
      data ppx,psx,ssx,dsx,dpx,ddx/'pp','ps','ss','ds','dp','dd'/
      pi54 = dsqrt(2.0d0)*pi**1.25d0
      gencon = .false.
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      jmax = nup(j2)
      imax = nup(i2)
      imax0 = imax - imin + 1
      jmax0 = jmax - jmin + 1
      icount = 0
      jcount = 0
c
      jmax = j2
      do 50 i = i1 , i2
         istar1 = nup(i) - nc(i) + 1
         istar2 = nup(i)
         if (ext1) jmax = i
         do 40 j = j1 , jmax
            ext4 = ext1 .and. i.eq.j
            jstar1 = nup(j) - nc(j) + 1
            jstar2 = nup(j)
            jcount = jcount + 1
            r2 = (cobas(j,1)-cobas(i,1))**2 + (cobas(j,2)-cobas(i,2))
     +           **2 + (cobas(j,3)-cobas(i,3))**2
            do 30 ip = istar1 , istar2
               if (ext4) jstar2 = ip
               pic1 = pe(ip)*cobas(i,1)
               pic2 = pe(ip)*cobas(i,2)
               pic3 = pe(ip)*cobas(i,3)
               do 20 jp = jstar1 , jstar2
                  icount = icount + 1
                  temp1(icount) = 1.0d0
                  xxa(icount) = pe(ip) + pe(jp)
                  xxk(icount) = -(r2*pe(ip)*pe(jp))
                  xxt(icount) = pc(ip)*pc(jp)
                  if (ext4 .and. ip.eq.jp) then
                     xxt(icount) = xxt(icount)*0.5d0
                     temp1(icount) = 2.0d0
                  end if
                  xxp1(icount) = pic1 + pe(jp)*cobas(j,1)
                  xxp2(icount) = pic2 + pe(jp)*cobas(j,2)
                  xxp3(icount) = pic3 + pe(jp)*cobas(j,3)
 20            continue
 30         continue
 40      continue
 50   continue
      icou = icount
      if (gencon) then
         do 60 i = 1 , icount
            a = 1.0d0/xxa(i)
            xxp1(i) = xxp1(i)*a
            xxp2(i) = xxp2(i)*a
            xxp3(i) = xxp3(i)*a
            xxk(i) = pi54*dexp(xxk(i)*a)*a
            xxt(i) = dabs(xxt(i))
 60      continue
      else
         do 70 i = 1 , icount
            a = 1.0d0/xxa(i)
            xxp1(i) = xxp1(i)*a
            xxp2(i) = xxp2(i)*a
            xxp3(i) = xxp3(i)*a
            xxk(i) = pi54*xxt(i)*dexp(xxk(i)*a)*a
            xxt(i) = temp1(i)
 70      continue
      end if
c
      if (buff.eq.psx) call psex(xxa,xxk,xxp1,xxp2,xxp3,xxt,temp1,
     +                           icount,i1,imax0,jmax0)
      if (buff.eq.ppx) then
         call ppex(xxa,xxk,xxp1,xxp2,xxp3,xxt,i1,imax0,j1,jmax0,temp1,
     +             wksp,icount,ext1)
      end if
      if (buff.eq.ssx) call ssex(xxa,xxk,xxt,icount)
      if (buff.eq.dpx) then
         call tizex(temp1,ifreem,3,2,3,2,4,4,7)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.ddx) then
         call tizex(temp1,ifreem,3,3,3,3,5,5,9)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.dsx) then
         call tizex(temp1,ifreem,3,1,3,1,3,3,5)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      do 80 i = 1 , icount
         if (xxt(i).lt.1.0d-30) xxt(i) = 1.0d-30
 80   continue
c
c**** some messing needed if a generalised contraction
c
      if (ext2 .and. gencon) then
         jcount = 0
         icount = 0
         do 120 i = i1 , i2
            istar1 = nup(i) - nc(i) + 1
            istar2 = nup(i)
            if (ext1) jmax = i
            do 110 j = j1 , jmax
               ext4 = ext1 .and. i.eq.j
               jstar1 = nup(j) - nc(j) + 1
               jstar2 = nup(j)
               jcount = jcount + 1
               do 100 ip = istar1 , istar2
                  if (ext4) jstar2 = ip
                  do 90 jp = jstar1 , jstar2
                     icount = icount + 1
                     if (ext4) then
                        xxk(jcount) = 2.0d0*xxk(jcount)
                        xxt(jcount) = 2.0d0*xxt(jcount)
                     end if
 90               continue
 100           continue
 110        continue
 120     continue
      end if
c
c****** if this is the kl parameter saving do the necessary
c***** counting to provide the /save/ arrays
c
      if (ext3) then
c
c***** if not enough memory to cope with the batch sizes as they stand
c***** further split the kl batch
c

         ipszb = ipsize
         icszb = icsize
         nbreak = 0
         jmax = j2
         jmin = j1
         imin = i1
         icount = 0
         icold = 0
         icx = 0
         jcount = 0
         tm = 0.0d0
         do 160 i = i1 , i2
            istar1 = nup(i) - nc(i) + 1
            istar2 = nup(i)
            if (ext1) jmax = i
            do 150 j = j1 , jmax
               ext4 = ext1 .and. i.eq.j
               jstar1 = nup(j) - nc(j) + 1
               jstar2 = nup(j)
               jcount = jcount + 1
               do 140 ip = istar1 , istar2
                  if (ext4) jstar2 = ip
                  do 130 jp = jstar1 , jstar2
                     icount = icount + 1
                     icold = icold + 1
                     iscat1(icold) = jcount
                     tm = max(tm,xxt(icount+icx))
 130              continue
 140           continue
               if (j.eq.jmax) then
                  ic = icount + nc(i+1)*nc(j1)
               else
                  ic = icount + nc(i)*nc(j+1)
               end if
               if (jcount.eq.icsize .or. ic.gt.ipsize .or.
     +             (i.eq.i2 .and. j.eq.j2)) then
                  nbreak = nbreak + 1
                  tmax(nbreak) = tm
                  info(1,nbreak) = imin
                  info(2,nbreak) = i
                  info(3,nbreak) = jmin
                  info(4,nbreak) = j
                  info(5,nbreak) = icx
                  icx = icx + icount
                  icst = icount
                  icount = 0
                  jcst = jcount
                  jcount = 0
                  tm = 0.0d0
                  if (j.eq.jmax) then
                     jmin = j1
                     imin = i + 1
                  else
                     imin = i
                     jmin = j + 1
                  end if
               end if
 150        continue
 160     continue
         if (nbreak.gt.nbrkmx) then
            write (iwr,*) 'nbreak is ' , nbreak
            write (iwr,*) 'which is too big for the parameters set'
            call caserr('fatal error detected in integral evaluation')
         end if
         if (nbreak.gt.1 .and. gdbg(1) ) then
            write (iwr,*) 'nbreak is ' , nbreak
         end if
         if (nbreak.eq.1) then
c**** memory taken up by primitives
            ii1 = icst*mempri
c**** memory taken up by contracted (mustnt exceed ii1)
            ji1 = jcst*memcon
            if (ii1.ge.ji1) then
               ipszb = icst
               icszb = jcst
            else
               icszb = jcst
               ipszb = (ji1/mempri) + 1
            end if
         end if 
c 
      end if
      return
      end
      subroutine sav8(iscat1,xxk,xxa,xxt,xxp1,xxp2,xxp3,igath,
     &                ipar,icou,i1,i2,j1,j2,ext1,ext3)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c***** a vectorised version of the parameter saving routine that proved
c***** necessary with the mp2 code. there is still a problem here
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext4,ext3
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
INCLUDE(common/mapper)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/iofile/iread,iwr
      dimension xxk(ipar),xxa(ipar),xxt(ipar),iscat1(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar)
     &,igath(ipar)
c     dimension wksp(ipar*3),array1(nb),array2(nb)
c
c     character *2  ppx,psx,ssx,dsx,dpx,ddx,buff
c     data ppx,psx,ssx,dsx,dpx,ddx/'pp','ps','ss','ds','dp','dd'/
      pi54 = dsqrt(2.0d0)*pi**1.25d0
c
      iminp = nup(i1) - nc(i1) + 1
      jminp = nup(j1) - nc(j1) + 1
      jmaxp = nup(j2)
      imaxp = nup(i2)
      icount = 0
c
      jmax = jmaxp
      do 30 ip = iminp , imaxp
         if (ext1) jmax = ip
         pic1 = pe(ip)*coprm(ip,1)
         pic2 = pe(ip)*coprm(ip,2)
         pic3 = pe(ip)*coprm(ip,3)
         do 20 jp = jminp , jmax
            icount = icount + 1
            xxk(icount) = pe(ip) + pe(jp)
            xxp1(icount) = pc(ip)*pc(jp)
            xxp2(icount) = pic1 + pe(jp)*coprm(jp,1)
            xxp3(icount) = pic2 + pe(jp)*coprm(jp,2)
            xxa(icount) = pic3 + pe(jp)*coprm(jp,3)
            xxt(icount) = -pe(ip)*pe(jp)
     +                    *((coprm(jp,1)-coprm(ip,1))**2+(coprm(jp,2)
     +                    -coprm(ip,2))**2+(coprm(jp,3)-coprm(ip,3))**2)
 20      continue
         if (ext1) xxp1(icount) = xxp1(icount)*0.5d0
 30   continue
      icou = icount
      do 40 i = 1 , icount
         a = 1.0d0/xxk(i)
         xxp2(i) = xxp2(i)*a
         xxp3(i) = xxp3(i)*a
         xxa(i) = xxa(i)*a
         xxt(i) = pi54*xxp1(i)*dexp(xxt(i)*a)*a
 40   continue
c
      jmax = j2
      jcount = 0
      icount = 0
      nupi = nup(i1-1)
      if (i1.eq.1) nupi = 0
      nupj = nup(j1-1)
      if (j1.eq.1) nupj = 0
c
      if (ext1) then
         do 80 i = i1 , i2
            istar2 = nup(i) - nupi
            if (istar2.gt.maxorb) call caserr('oh dear in sav8')
            istar1 = istar2 - nc(i) + 1
            do 70 j = j1 , i
               jcount = jcount + 1
               ext4 = (i.eq.j)
               jstar2 = nup(j) - nupj
               jstar1 = jstar2 - nc(j) + 1
               do 60 ip = istar1 , istar2
                  if (ext4) jstar2 = ip
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                  do 50 jp = jstar1 , jstar2
                     icount = icount + 1
                     igath(icount) = iky(ip) + jp
                     iscat1(icount) = jcount
 50               continue
 60            continue
 70         continue
 80      continue
      else
         jmaxp0 = jmaxp - jminp + 1
         do 120 i = i1 , i2
            istar2 = nup(i) - nupi
            if (istar2.gt.maxorb) call caserr('oh dear in sav8')
            istar1 = istar2 - nc(i) + 1
            do 110 j = j1 , jmax
               jcount = jcount + 1
               ext4 = ext1 .and. i.eq.j
               jstar2 = nup(j) - nupj
               jstar1 = jstar2 - nc(j) + 1
               do 100 ip = istar1 , istar2
                  ipg = (ip-1)*jmaxp0
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                  do 90 jp = jstar1 , jstar2
                     icount = icount + 1
                     igath(icount) = ipg + jp
                     iscat1(icount) = jcount
 90               continue
 100           continue
 110        continue
 120     continue
      end if
c
      call dgthr(icount,xxp2,xxp1,igath)
      call dgthr(icount,xxp3,xxp2,igath)
      call dgthr(icount,xxa,xxp3,igath)
      call dgthr(icount,xxk,xxa,igath)
      call dgthr(icount,xxt,xxk,igath)
c
c****** if this is the kl parameter saving do the necessary
c***** counting to provide the /save/ arrays
c
      if (ext3) then
c******* check that there is enough memory to do this in a single batch.
         ipszb = ipsize
         icszb = icsize
         nbreak = 1
         info(1,nbreak) = i1
         info(2,nbreak) = i2
         info(3,nbreak) = j1
         info(4,nbreak) = j2
         info(5,nbreak) = 0
         if (jcount.eq.icsize .or. icount.gt.ipsize) then
            write (iwr,*)
     +              'the parameter setting for the kl batch are too big'
            write (iwr,*) 'i1,i2,j1,j2 are ' , i1 , i2 , j1 , j2
            call caserr('fatal error detected in integral evaluation')
         end if
c**** memory taken up by primitives
         ii1 = icount*mempri
c**** memory taken up by contracted (mustnt exceed ii1)
         ji1 = jcount*memcon
         if (ii1.ge.ji1) then
            ipszb = icount
            icszb = jcount
         else
            icszb = jcount
            ipszb = (ji1/mempri) + 1
         end if
c
      end if
      return
      end
      subroutine sav9(iscat1,xxk,xxa,xxt,xxp1,xxp2,xxp3,temp1
     &,wksp,ipar,icou,i1,i2,j1,j2,ext1,buff,ifreem
     &,array1,array2)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c**** yet another array that saves parameters. it is used to calculate
c**** exchange integrals prior to an mp2 calculation
c
      implicit REAL  (a-h,o-z)
      character *2  ppx,psx,ssx,dsx,dpx,ddx,buff
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext4
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
      dimension xxk(ipar),xxa(ipar),xxt(ipar),iscat1(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar)
     &,temp1(ipar*3),wksp(ipar*3),array1(nb),array2(nb)
c
      data ppx,psx,ssx,dsx,dpx,ddx/'pp','ps','ss','ds','dp','dd'/
      pi54 = dsqrt(2.0d0)*pi**1.25d0
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      jmax = nup(j2)
      imax = nup(i2)
      imax0 = imax - imin + 1
      jmax0 = jmax - jmin + 1
      icount = 0
c
      jmax = j2
      do 50 i = i1 , i2
         istar1 = nup(i) - nc(i) + 1
         istar2 = nup(i)
         if (ext1) jmax = i
         do 40 j = j1 , jmax
            ext4 = ext1 .and. i.eq.j
            jstar1 = nup(j) - nc(j) + 1
            jstar2 = nup(j)
            r2 = (cobas(j,1)-cobas(i,1))**2 + (cobas(j,2)-cobas(i,2))
     +           **2 + (cobas(j,3)-cobas(i,3))**2
            do 30 ip = istar1 , istar2
               if (ext4) jstar2 = ip
               pic1 = pe(ip)*cobas(i,1)
               pic2 = pe(ip)*cobas(i,2)
               pic3 = pe(ip)*cobas(i,3)
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
               do 20 jp = jstar1 , jstar2
                  icount = icount + 1
                  temp1(icount) = 1.0d0
                  xxa(icount) = pe(ip) + pe(jp)
                  xxk(icount) = -(r2*pe(ip)*pe(jp))
                  xxt(icount) = pc(ip)*pc(jp)
                  if (ext4 .and. ip.eq.jp) then
                     xxt(icount) = xxt(icount)*0.5d0
                     temp1(icount) = 2.0d0
                  end if
                  xxp1(icount) = pic1 + pe(jp)*cobas(j,1)
                  xxp2(icount) = pic2 + pe(jp)*cobas(j,2)
                  xxp3(icount) = pic3 + pe(jp)*cobas(j,3)
 20            continue
 30         continue
 40      continue
 50   continue
      icou = icount
      do 60 i = 1 , icount
         a = 1.0d0/xxa(i)
         xxp1(i) = xxp1(i)*a
         xxp2(i) = xxp2(i)*a
         xxp3(i) = xxp3(i)*a
         xxk(i) = pi54*xxt(i)*dexp(xxk(i)*a)*a
         xxt(i) = temp1(i)
 60   continue
c
      if (buff.eq.psx) call psex(xxa,xxk,xxp1,xxp2,xxp3,xxt,temp1,
     +                           icount,i1,imax0,jmax0)
      if (buff.eq.ppx) then
         call ppex(xxa,xxk,xxp1,xxp2,xxp3,xxt,i1,imax0,j1,jmax0,temp1,
     +             wksp,icount,ext1)
      end if
      if (buff.eq.ssx) call ssex(xxa,xxk,xxt,icount)
      if (buff.eq.dpx) then
         call tizex(temp1,ifreem,3,2,3,2,4,4,7)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.ddx) then
         call tizex(temp1,ifreem,3,3,3,3,5,5,9)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.dsx) then
         call tizex(temp1,ifreem,3,1,3,1,3,3,5)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
c
      jmax = j2
      icount = 0
      jcount = 0
      do 90 i = i1 , i2
         istar1 = nc(i)
         iipmax = istar1*(istar1+1)/2
         if (ext1) jmax = i
         do 80 j = j1 , jmax
            jcount = jcount + 1
            ext4 = ext1 .and. i.eq.j
            jstar1 = nc(j)
            ijpmax = istar1*jstar1
            if (ext4) ijpmax = iipmax
c
            txxx1 = array1(i)*array2(j)
            txxx2 = array1(j)*array2(i)
            txxx = max(txxx1,txxx2)
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
            do 70 ip = 1 , ijpmax
               icount = icount + 1
               iscat1(icount) = jcount
               xxt(icount) = xxt(icount)*txxx
 70         continue
 80      continue
 90   continue
      tm = 0.0d0
      do 100 i = 1 , icount
         tm = max(tm,xxt(i))
         if (xxt(i).lt.1.0d-30) xxt(i) = 1.0d-30
 100  continue
      return
      end
      subroutine savdrv(iscat1,xxk,xxa,xxt,xxp1,xxp2,xxp3,temp1
     &,wksp,ipar,icou,i1,i2,j1,j2,ext1,ext2,ext3,buff,ifreem,xxc)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c******* parameter saving code used for gradient
      implicit REAL  (a-h,o-z)
      character *2  ppx,psx,ssx,dsx,dpx,ddx,buff
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext2,gencon,ext4,ext3
INCLUDE(common/sizes)
      common/save/info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      common/iofile/iread,iwr
      logical glog,gdbg
      common/gen/maxit(2),concrit(5),glog(5),gdbg(10)
      dimension xxk(ipar),xxa(ipar),xxt(ipar),iscat1(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar),xxc(ipar)
     &,temp1(ipar*3),wksp(ipar*3)
c
      data ppx,psx,ssx,dsx,dpx,ddx/'pp','ps','ss','ds','dp','dd'/
      pi54 = dsqrt(2.0d0)*pi**1.25d0
      gencon = .false.
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      jmax = nup(j2)
      imax = nup(i2)
      imax0 = imax - imin + 1
      jmax0 = jmax - jmin + 1
      icount = 0
      jcount = 0
c
      jmax = j2
      do 50 i = i1 , i2
         istar1 = nup(i) - nc(i) + 1
         istar2 = nup(i)
         if (ext1) jmax = i
         do 40 j = j1 , jmax
            ext4 = ext1 .and. i.eq.j
            jstar1 = nup(j) - nc(j) + 1
            jstar2 = nup(j)
            jcount = jcount + 1
            r2 = (cobas(j,1)-cobas(i,1))**2 + (cobas(j,2)-cobas(i,2))
     +           **2 + (cobas(j,3)-cobas(i,3))**2
            do 30 ip = istar1 , istar2
               if (ext4) jstar2 = ip
               pic1 = pe(ip)*cobas(i,1)
               pic2 = pe(ip)*cobas(i,2)
               pic3 = pe(ip)*cobas(i,3)
               do 20 jp = jstar1 , jstar2
                  icount = icount + 1
                  temp1(icount) = 1.0d0
                  xxa(icount) = pe(ip) + pe(jp)
                  xxc(icount) = pe(ip)*2.0d0
                  xxk(icount) = -(r2*pe(ip)*pe(jp))
                  xxt(icount) = pc(ip)*pc(jp)
                  if (ext4 .and. ip.eq.jp) then
                     xxt(icount) = xxt(icount)*0.5d0
                     temp1(icount) = 2.0d0
                  end if
                  xxp1(icount) = pic1 + pe(jp)*cobas(j,1)
                  xxp2(icount) = pic2 + pe(jp)*cobas(j,2)
                  xxp3(icount) = pic3 + pe(jp)*cobas(j,3)
 20            continue
 30         continue
 40      continue
 50   continue
      icou = icount
      if (gencon) then
         do 60 i = 1 , icount
            a = 1.0d0/xxa(i)
            xxp1(i) = xxp1(i)*a
            xxp2(i) = xxp2(i)*a
            xxp3(i) = xxp3(i)*a
            xxk(i) = pi54*dexp(xxk(i)*a)*a
            xxt(i) = dabs(xxt(i))
 60      continue
      else
         do 70 i = 1 , icount
            a = 1.0d0/xxa(i)
            xxp1(i) = xxp1(i)*a
            xxp2(i) = xxp2(i)*a
            xxp3(i) = xxp3(i)*a
            xxk(i) = pi54*xxt(i)*dexp(xxk(i)*a)*a
            xxt(i) = temp1(i)
 70      continue
      end if
c
      if (buff.eq.psx) call psex(xxa,xxk,xxp1,xxp2,xxp3,xxt,temp1,
     +                           icount,i1,imax0,jmax0)
      if (buff.eq.ppx) then
         call ppex(xxa,xxk,xxp1,xxp2,xxp3,xxt,i1,imax0,j1,jmax0,temp1,
     +             wksp,icount,ext1)
      end if
      if (buff.eq.ssx) call ssex(xxa,xxk,xxt,icount)
      if (buff.eq.dpx) then
         call tizex(temp1,ifreem,3,2,3,2,4,4,7)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.ddx) then
         call tizex(temp1,ifreem,3,3,3,3,5,5,9)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      if (buff.eq.dsx) then
         call tizex(temp1,ifreem,3,1,3,1,3,3,5)
         call exchng(temp1,xxk,xxa,xxt,xxp1,xxp2,xxp3,icount,ext1,i1,i2,
     +               j1,j2)
      end if
      do 80 i = 1 , icount
         if (xxt(i).lt.1.0d-30) xxt(i) = 1.0d-30
 80   continue
c
c**** some messing needed if a generalised contraction
c
      if (ext2 .and. gencon) then
         jcount = 0
         icount = 0
         do 120 i = i1 , i2
            istar1 = nup(i) - nc(i) + 1
            istar2 = nup(i)
            if (ext1) jmax = i
            do 110 j = j1 , jmax
               ext4 = ext1 .and. i.eq.j
               jstar1 = nup(j) - nc(j) + 1
               jstar2 = nup(j)
               jcount = jcount + 1
               do 100 ip = istar1 , istar2
                  if (ext4) jstar2 = ip
                  do 90 jp = jstar1 , jstar2
                     icount = icount + 1
                     if (ext4) then
                        xxk(jcount) = 2.0d0*xxk(jcount)
                        xxt(jcount) = 2.0d0*xxt(jcount)
                     end if
 90               continue
 100           continue
 110        continue
 120     continue
      end if
c
c****** if this is the kl parameter saving do the necessary
c***** counting to provide the /save/ arrays
c
      if (ext3) then
         ipszb = ipsize
         icszb = icsize
c        itempo = (jcount+1)/2
         nbreak = 0
         jmax = j2
         jmin = j1
         imin = i1
         icount = 0
         icold = 0
         icx = 0
         jcount = 0
         tm = 0.0d0
         do 160 i = i1 , i2
            istar1 = nup(i) - nc(i) + 1
            istar2 = nup(i)
            if (ext1) jmax = i
            do 150 j = j1 , jmax
               ext4 = ext1 .and. i.eq.j
               jstar1 = nup(j) - nc(j) + 1
               jstar2 = nup(j)
               jcount = jcount + 1
               do 140 ip = istar1 , istar2
                  if (ext4) jstar2 = ip
                  do 130 jp = jstar1 , jstar2
                     icount = icount + 1
                     icold = icold + 1
                     iscat1(icold) = jcount
                     tm = max(tm,xxt(icount+icx))
 130              continue
 140           continue
               if (j.eq.jmax) then
                  ic = icount + nc(i+1)*nc(j1)
               else
                  ic = icount + nc(i)*nc(j+1)
               end if
               if (jcount.eq.icsize .or. ic.gt.ipsize .or.
     +             (i.eq.i2 .and. j.eq.j2)) then
c              if (jcount.eq.itempo.or.ic.gt.ipsize.
c     &             or.(i.eq.i2.and.j.eq.j2)) then
                  nbreak = nbreak + 1
                  tmax(nbreak) = tm
                  info(1,nbreak) = imin
                  info(2,nbreak) = i
                  info(3,nbreak) = jmin
                  info(4,nbreak) = j
                  info(5,nbreak) = icx
                  icx = icx + icount
                  icst = icount
                  icount = 0
                  jcst = jcount
                  jcount = 0
                  tm = 0.0d0
                  if (j.eq.jmax) then
                     jmin = j1
                     imin = i + 1
                  else
                     imin = i
                     jmin = j + 1
                  end if
               end if
 150        continue
 160     continue
         if (nbreak.gt.nbrkmx) then
            write (iwr,*) 'nbreak is ' , nbreak
            write (iwr,*) 'which is too big for the parameters set'
            call caserr('fatal error detected in derivative evaluation')
         end if
         if (nbreak.gt.1 .and. gdbg(1) ) then
            write (iwr,*) 'nbreak is ' , nbreak
         end if
         if (nbreak.eq.1) then
c**** memory taken up by primitives
            ii1 = icst*mempri
c**** memory taken up by contracted (mustnt exceed ii1)
            ji1 = jcst*memcon
            if (ii1.ge.ji1) then
               ipszb = icst
               icszb = jcst
            else
               icszb = jcst
               ipszb = (ji1/mempri) + 1
            end if
         end if
c
      end if
      return
      end
      subroutine savonb(xxb,xxa,xxk,xxc,xxp1,xxp2,xxp3,temp1
     &,ipar,icount,i1,i2,j1,j2,ext1)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c***** save parameters for the calculation of the derivatives of the
c***** nuclear attraction integrals.
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext4
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      dimension xxa(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar),xxc(ipar)
     &,temp1(ipar),xxb(ipar),xxk(ipar)
c
c     pi32 = pi**1.5d0
c
c     imin = nup(i1) - nc(i1) + 1
c     jmin = nup(j1) - nc(j1) + 1
c     jmax = nup(j2)
c     imax = nup(i2)
c     imax0 = imax - imin + 1
c     jmax0 = jmax - jmin + 1
      icount = 0
c
      jmax = j2
      do 50 i = i1 , i2
         istar1 = nup(i) - nc(i) + 1
         istar2 = nup(i)
         if (ext1) jmax = i
         do 40 j = j1 , jmax
            ext4 = ext1 .and. i.eq.j
            jstar1 = nup(j) - nc(j) + 1
            jstar2 = nup(j)
            r2 = -((cobas(j,1)-cobas(i,1))**2+(cobas(j,2)-cobas(i,2))
     +           **2+(cobas(j,3)-cobas(i,3))**2)
            do 30 ip = istar1 , istar2
               if (ext4) jstar2 = ip
               pic1 = pe(ip)*cobas(i,1)
               pic2 = pe(ip)*cobas(i,2)
               pic3 = pe(ip)*cobas(i,3)
               do 20 jp = jstar1 , jstar2
                  icount = icount + 1
                  xxa(icount) = pe(ip) + pe(jp)
                  xxb(icount) = pe(ip)*pe(jp)
                  xxc(icount) = pe(ip)*2.0d0
                  xxk(icount) = r2
                  xxp1(icount) = pic1 + pe(jp)*cobas(j,1)
                  xxp2(icount) = pic2 + pe(jp)*cobas(j,2)
                  xxp3(icount) = pic3 + pe(jp)*cobas(j,3)
                  temp1(icount) = 4.0d0*pc(ip)*pc(jp)
                  if (ext4 .and. ip.eq.jp) then
                     temp1(icount) = 0.5d0*temp1(icount)
                  end if
 20            continue
 30         continue
 40      continue
 50   continue
      do 60 i = 1 , icount
         a = 1.0d0/xxa(i)
         b = xxb(i)*a
         xxp1(i) = xxp1(i)*a
         xxp2(i) = xxp2(i)*a
         xxp3(i) = xxp3(i)*a
         xxk(i) = temp1(i)*pi*a*dexp(xxk(i)*b)
 60   continue
c
      ipszb = icount
      if (ipszb.gt.ipsize) call caserr(
     +            'ran out of space for derivative 1 electron integrals'
     +            )
      return
      end
      subroutine savone(zmem,xxa,xxb,xxc,xxp1,xxp2,xxp3,temp1
     &,ipar,icount,i1,i2,j1,j2,ext1)
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c
c****** routine that saves parameters for the derivative one electron
c****** integrals.
c
      implicit REAL  (a-h,o-z)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      logical ext1,ext4
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4
      common/defunk/cobas(maxorb,3)
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsize,icsize
     &,ipszb,icszb,memcon
      dimension xxa(ipar)
     &,xxp1(ipar),xxp2(ipar),xxp3(ipar),xxc(ipar)
     &,temp1(ipar),xxb(ipar),zmem(ipar*2)
c
      pi32 = pi**1.5d0
c
c     imin = nup(i1) - nc(i1) + 1
c     jmin = nup(j1) - nc(j1) + 1
c     jmax = nup(j2)
c     imax = nup(i2)
c     imax0 = imax - imin + 1
c     jmax0 = jmax - jmin + 1
      icount = 0
c
      jmax = j2
      do 50 i = i1 , i2
         istar1 = nup(i) - nc(i) + 1
         istar2 = nup(i)
         if (ext1) jmax = i
         do 40 j = j1 , jmax
            ext4 = ext1 .and. i.eq.j
            jstar1 = nup(j) - nc(j) + 1
            jstar2 = nup(j)
            r2 = -((cobas(j,1)-cobas(i,1))**2+(cobas(j,2)-cobas(i,2))
     +           **2+(cobas(j,3)-cobas(i,3))**2)
            do 30 ip = istar1 , istar2
               if (ext4) jstar2 = ip
               pic1 = pe(ip)*cobas(i,1)
               pic2 = pe(ip)*cobas(i,2)
               pic3 = pe(ip)*cobas(i,3)
               do 20 jp = jstar1 , jstar2
                  icount = icount + 1
                  xxa(icount) = pe(ip) + pe(jp)
                  xxb(icount) = pe(ip)*pe(jp)
                  xxc(icount) = pe(ip)*2.0d0
                  zmem(icount) = r2
                  xxp1(icount) = pic1 + pe(jp)*cobas(j,1)
                  xxp2(icount) = pic2 + pe(jp)*cobas(j,2)
                  xxp3(icount) = pic3 + pe(jp)*cobas(j,3)
                  temp1(icount) = 2.0d0*pc(ip)*pc(jp)
                  if (ext4 .and. ip.eq.jp) then
                     temp1(icount) = 0.5d0*temp1(icount)
                  end if
 20            continue
 30         continue
 40      continue
 50   continue
      do 60 i = 1 , icount
         a = 1.0d0/xxa(i)
         xxa(i) = a*0.5d0
         b = xxb(i)*a
         xxb(i) = 2.0d0*b
         c = b*(3.0d0+xxb(i)*zmem(i))
         xxp1(i) = xxp1(i)*a
         xxp2(i) = xxp2(i)*a
         xxp3(i) = xxp3(i)*a
         zmem(i) = temp1(i)*(pi32*a**1.5d0)*dexp(zmem(i)*b)
         zmem(i+icount) = zmem(i)*c
 60   continue
c
      ipszb = icount
      if (ipszb.gt.ipsize) call caserr(
     +            'ran out of space for derivative 1 electron integrals'
     +            )
      return
      end
      subroutine scfit(f,shalf,p,gp,e)
c
c
c***** routine just forms the new density matrix from the fock matrix
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen
      logical glog1,glog2,glog5,glog6,glog9,gdbg,gmp2
     &,readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf,gmull
c****** gen contains most of the stuff about the options
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
     &,glog1,glog2,glog5,glog9,glog6,gdbg(10),gmp2,mp2fo,mp2fv
     &,readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,gdold
     &,gscf,noptd,npoint,timeg,gmull,iswap,dshift
INCLUDE(common/wrtd)
      common/iofile/iread,iwr
INCLUDE(common/mapper)
INCLUDE(common/harmon)
c
      dimension p(n,n),shalf(*)
      dimension f(n,n),gp(n,n),e(n)
c
c     *****  ne is number of electrons ; n is actual ***********
c     *****  number of basis functions *************************
c
c          obtain fock matrix in orthogonal basis
c
      call rdedx(gp,n*n,idblk(3),numdir)
      call tdownd(gp,ilifq,gp,ilifq,n)
      no = ne/2
      call vclr(p,1,n*n)
      call mxmb(f,1,n,gp,1,n,p,1,n,n,n,n)
      call vclr(f,1,n*n)
      call mxmb(gp,n,1,p,1,n,f,1,n,n,n,n)
c
      if (oharm) then
c...   we have some 0.0 diagonals that correspond to eliminated 
c...   functions in going from cartesian to harmonic => shift them away
         do i=newbas0+1,newbas1
            if (f(i,i).ne.0.0d0) call caserr(' harmonic problem ')
            f(i,i) = 999999999.0d0
         end do
      end if
c
c...  level shift
c
      ne2 = ne/2
      if (ne2*2.ne.ne) call caserr(' NOT closed shell ')
      do 13 i=ne2+1,n
13    f(i,i) = f(i,i) + dshift
c
c          diagonalise transformed f(now called w)
c          and transform e-vectors back to a.o. basis
c
c     ifail = 0
c     call f02abf(f,n,n,e,p,n,a,ifail)
c     abandoned use of nag routine, due to mxing of
c     degenerate mo coefficients
c
      call trianc(f,shalf,n,n)
c
      diaacc = 1.0d-11
      call jacobi(shalf,iky,n,p,ilifq,n,e,2,2,diaacc)
      call vclr(shalf,1,n*n)
      call rdedx(gp,n*n,idblk(3),numdir)
      call mxmb(gp,1,n,p,1,n,shalf,1,n,n,n,n)
c
c..   write new scf orbitals
c
      call wrt3(shalf,n*n,idblk(3),numdir)
c
c
c          use e.v.s  in gp to find p
c
      call tdownd(f,ilifq,shalf,ilifq,n)
      call vclr(p,1,n*n)
      if (iter.eq.2 .and. iswap.ne.0) then
         do 50 i = 1 , n
            do 40 j = 1 , n
               do 20 k = 1 , no - iswap
                  p(i,j) = p(i,j) + 2.0d0*f(i,k)*f(j,k)
 20            continue
               do 30 k = no + 1 , no + iswap
                  p(i,j) = p(i,j) + 2.0d0*f(i,k)*f(j,k)
 30            continue
 40         continue
 50      continue
         write (iwr,*) iswap , ' orbitals swopped here'
      else
         do 80 i = 1 , n
            do 70 j = 1 , n
               do 60 k = 1 , no
                  p(i,j) = p(i,j) + 2.0d0*f(i,k)*f(j,k)
 60            continue
 70         continue
 80      continue
      end if
c
      return
      end
      subroutine shalft(sx,hx,shalf,e)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,n,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3
INCLUDE(common/mapper)
      dimension sx(n,n),hx(n,n),shalf(n,n),e(n)
c
c symmetrise hx  and sx
c
      do 30 i = 1 , n
         do 20 j = 1 , i
            sx(j,i) = sx(i,j)
            hx(j,i) = hx(i,j)
 20      continue
 30   continue
c
c        find shalf the transformation matrix
c
c     ifail = 0
c     call f02abf(sx,n,n,e,shalf,n,a,ifail)
c     f02abf removed due to mixing of degenerate mos ...
c
      call trianc(sx,shalf,n,n)
      diaacc = 1.0d-11
      call jacobi(shalf,iky,n,sx,ilifq,n,e,2,2,diaacc)
      do 50 i = 1 , n
         do 40 j = 1 , n
            shalf(j,i) = sx(j,i)/(dsqrt(e(i)))
 40      continue
 50   continue
c
      return
      end
      subroutine sinter(xklp1,xklp2,xklp3,wq,xkla,rho,xklk
     &,w,t,temp,f,xijp1,xijp2,xijp3,xija,xijk
     &,mx,jcdash,ips,ff)
c
c****** routine calculates the parameters that depend on i,j k, and l
c****** that are used by the recursion relations to calculate the
c****** integrals. (the (ss|ss)m integrals are calculated here)
c
      implicit REAL  (a-h,o-z)
c
      dimension f(ips,mx),t(ips),rho(ips),temp(ips)
     &,w(ips,3),xklp1(ips),xklp2(ips),xklp3(ips)
     &,xklk(ips),xkla(ips),ff(jcdash,mx),wq(ips,3)
c
c**** note:  w(ips,2) is equivalenced to t
c**** note:  w(ips,3) is equivalenced to temp
c
      do 20 jcount = 1 , jcdash
         rho(jcount) = 1.0d0/(xija+xkla(jcount))
         x1 = (xijp1-xklp1(jcount))**2 + (xijp2-xklp2(jcount))
     +        **2 + (xijp3-xklp3(jcount))**2
c........ t is the asrgument for the incomplete gamma function eval.
         t(jcount) = rho(jcount)*x1*xija*xkla(jcount)
         temp(jcount) = xklk(jcount)*xijk*dsqrt(rho(jcount))
 20   continue
c**** in this next call f is being used as temporary storage
c
      call fquik2(ff,ff,f,t,mx,jcdash)
c
      do 40 m = 1 , mx
         do 30 jcount = 1 , jcdash
c.......(ss|ss) integrals
            f(jcount,m) = temp(jcount)*ff(jcount,m)
 30      continue
 40   continue
c
      psssx = xija*xijp1
      psssy = xija*xijp2
      psssz = xija*xijp3
      do 50 jcount = 1 , jcdash
         wa = rho(jcount)*(psssx+xkla(jcount)*xklp1(jcount))
         wb = rho(jcount)*(psssy+xkla(jcount)*xklp2(jcount))
         wc = rho(jcount)*(psssz+xkla(jcount)*xklp3(jcount))
c........... w-q
         wq(jcount,1) = wa - xklp1(jcount)
         wq(jcount,2) = wb - xklp2(jcount)
         wq(jcount,3) = wc - xklp3(jcount)
c........... w-p
         w(jcount,1) = wa - xijp1
         w(jcount,2) = wb - xijp2
         w(jcount,3) = wc - xijp3
 50   continue
c
      return
      end
      subroutine snterd(xklp1,xklp2,xklp3,xkla,xklk
     &,t,f,ff,mx,jcdash,ips,charge)
c
c**** rotine calculates the (s|a|s)m integrals needed for the calculatio
c**** of the derivative nuclear attraction integrals. some of the
c**** quantities used are defined in sinter.
c
      implicit REAL  (a-h,o-z)
c
      dimension f(ips,mx),t(ips)
     &,xklp1(ips),xklp2(ips),xklp3(ips)
     &,xklk(ips),xkla(ips),ff(jcdash,mx)
c
      do 20 jcount = 1 , jcdash
         t(jcount) = (xklp1(jcount)**2+xklp2(jcount)**2+xklp3(jcount)
     +               **2)*xkla(jcount)
 20   continue
c**** in this next call f is being used as temporary storage
c
      call fquik2(ff,ff,f,t,mx,jcdash)
c
      do 40 m = 1 , mx
         do 30 jcount = 1 , jcdash
            f(jcount,m) = xklk(jcount)*ff(jcount,m)*charge
 30      continue
 40   continue
c
      return
      end
      subroutine ss1(xijk,xija,xijp1,xijp2,xijp3,i1,i2,j1
     &,j2,xt1,xt2,wx,wy,wz,f,t,lar,hx,sx,s1,h1,ipar2,ext1)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c     *****************************************************
c     ******* evaluating 1 electron integrals for  ********
c     **************** (s/s) and (s/h/s) ******************
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
_IFN1(cfu)      parameter (pi=3.1415926535898d0)
_IF1(cfu)      parameter (pi=3.1415926535898d0)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16
      logical ext1
      dimension xija(n11),xijk(n11),xijp1(n11),xijp2(n11)
     &,s1(n11),h1(n11),sx(nb,nb),hx(nb,nb),xijp3(n11)
     &,lar(n11),f(5*n11),t(n11),xt1(n14,n13),xt2(n16,n15)
     &,wx(ipar2,ipar2),wy(ipar2,ipar2),wz(ipar2,ipar2)
      pi14 = pi**0.25d0/dsqrt(2.0d0)
c
c     ******* (s/s) and (s/h(1)/s) integrals ********
c
      imin = nup(i1) - nc(i1) + 1
      jmin = nup(j1) - nc(j1) + 1
      imax = nup(i2)
      jmax = nup(j2)
      call vclr(h1,1,n11)
      do 40 k = 1 , na
c
         do 20 icount = 1 , n11
            t(icount) = xija(icount)
     +                  *((xijp1(icount)-coord(1,k))**2+(xijp2(icount)
     +                  -coord(2,k))**2+(xijp3(icount)-coord(3,k))**2)
 20      continue
         call fquik2(f,f,lar,t,1,n11)
         do 30 icount = 1 , n11
            h1(icount) = h1(icount) + charge(k)*f(icount)
 30      continue
 40   continue
c
      icount = 0
      do 60 i = imin , imax
         if (ext1) jmax = i
         do 50 j = jmin , jmax
            icount = icount + 1
            t(icount) = (coprm(j,1)-coprm(i,1))
     +                  **2 + (coprm(j,2)-coprm(i,2))
     +                  **2 + (coprm(j,3)-coprm(i,3))**2
            f(icount) = pe(i)*pe(j)
 50      continue
 60   continue
c
      do 70 icount = 1 , n11
         x1 = 1.0d0/xija(icount)
         s1(icount) = xijk(icount)*pi14*dsqrt(x1)
         x2 = f(icount)*x1
         h1(icount) = s1(icount)*x2*(3.0d0-2.0d0*t(icount)*x2)
     +                - xijk(icount)*h1(icount)/pi14
 70   continue
c
c****** put the result in the matrices to be transformed ***
c
      if (ext1) then
         icount = 0
         do 90 i = 1 , n14
            do 80 j = 1 , i
               icount = icount + 1
               wx(i,j) = s1(icount)
               wy(i,j) = h1(icount)
               wx(j,i) = s1(icount)
               wy(j,i) = h1(icount)
 80         continue
 90      continue
      else
         icount = 0
         do 110 i = 1 , n14
            do 100 j = 1 , n16
               icount = icount + 1
               wx(i,j) = s1(icount)
               wy(i,j) = h1(icount)
 100        continue
 110     continue
      end if
c
c******** transform the matrices **********
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wx,1,ipar2,xt2,1,n16,wz,1,ipar2,n14,n16,n15)
      call vclr(wx,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14,wx,ipar2,1,n15,n14,n13)
c
      call vclr(wz,1,ipar2*ipar2)
      call mxmb(wy,1,ipar2,xt2,1,n16,wz,1,ipar2,n14,n16,n15)
      call vclr(wy,1,ipar2*ipar2)
      call mxmb(wz,ipar2,1,xt1,1,n14,wy,ipar2,1,n15,n14,n13)
c
      jmax = j2
      ic = 0
      do 130 i = i1 , i2
         ic = ic + 1
         jc = 0
         if (ext1) jmax = i
         do 120 j = j1 , jmax
            jc = jc + 1
            sx(i,j) = wx(ic,jc)
            hx(i,j) = wy(ic,jc)
 120     continue
 130  continue
c
      return
      end
      subroutine ssex(xija,xijk,xijt,icmax)
c
c     *****************************************************
c     ******* evaluating exchange integrals of type (ss/ss) ********
c     *****************************************************
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16
      common/gen/maxit,iter,concrit,conv,energy,thresh,dcmax
      dimension xija(icmax),xijk(icmax),xijt(icmax)
c
c********* the exchange integral evaluation *****
c
      do 20 icount = 1 , icmax
         xijt(icount) = xijt(icount)*dabs(xijk(icount))
     +                  *(xija(icount)+xija(icount))**(-0.25d0)
 20   continue
c
      return
      end
      subroutine symdrc(grad,ict)
c
c========================================
c     symmetrise gradient
c==========================================
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
c
      dimension grad(na,3),ict(na,48)
c
INCLUDE(common/symtry)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480)
c
      nav = lenwrd()
      call readi(ict,nw196(6)*nav,ibl196(6),idaf)
c
c     ----- symmetryze gradient vector -----
c
      zero = 0.0d0
      do 60 ic = 1 , na
         do 20 it = 1 , nt
            if (ict(ic,it).gt.ic) go to 60
 20      continue
         dedx = zero
         dedy = zero
         dedz = zero
         do 30 it = 1 , nt
            icnu = ict(ic,it)
            dedxp = grad(icnu,1)
            dedyp = grad(icnu,2)
            dedzp = grad(icnu,3)
            n = 3*(it-1)
            dedx = dedx + dedxp*ptr(1,n+1) + dedyp*ptr(2,n+1)
     +             + dedzp*ptr(3,n+1)
            dedy = dedy + dedxp*ptr(1,n+2) + dedyp*ptr(2,n+2)
     +             + dedzp*ptr(3,n+2)
            dedz = dedz + dedxp*ptr(1,n+3) + dedyp*ptr(2,n+3)
     +             + dedzp*ptr(3,n+3)
 30      continue
         grad(ic,1) = dedx
         grad(ic,2) = dedy
         grad(ic,3) = dedz
         do 50 it = 1 , nt
            icnu = ict(ic,it)
            if (icnu.ne.ic) then
               if (it.ne.nt) then
                  it1 = it + 1
                  do 40 jt = it1 , nt
                     if (ict(ic,jt).eq.icnu) go to 50
 40               continue
               end if
               jt = invt(it)
               n = 3*(jt-1)
               grad(icnu,1) = grad(ic,1)*ptr(1,n+1) + grad(ic,2)
     +                        *ptr(2,n+1) + grad(ic,3)*ptr(3,n+1)
               grad(icnu,2) = grad(ic,1)*ptr(1,n+2) + grad(ic,2)
     +                        *ptr(2,n+2) + grad(ic,3)*ptr(3,n+2)
               grad(icnu,3) = grad(ic,1)*ptr(1,n+3) + grad(ic,2)
     +                        *ptr(2,n+3) + grad(ic,3)*ptr(3,n+3)
            end if
c
 50      continue
 60   continue
      dum = dfloat(nt)
      do 80 n = 1 , na
         do 70 i = 1 , 3
            grad(n,i) = grad(n,i)/dum
 70      continue
 80   continue
      return
      end
      subroutine symhd(f,h,ia)
      implicit REAL  (a-h,o-z)
      logical iandj
c
c     ----- symmetrize the skeleton fock matrix
c
      dimension f(*),h(*),ia(*)
c
INCLUDE(common/sizes)
INCLUDE(common/symtry)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3
      common/defunk/cobas(maxorb,3),nd
c
      common/scra  /iso(mxshel,48),ntd
      common/hsym/t(20,20),mini,maxi,lit,minj,maxj,ljt,ntr
      dimension indin(48),mi(48)
      data zero,one/0.0d0,1.0d0/
      if (ntd.eq.1) return
      nsnp = ns + np
      nshell = ns + np + nd
      nx = nb*(nb+1)/2
      do 20 i = 1 , nx
         h(i) = zero
 20   continue
c     ----- find a block (i,j)
      do 240 ii = 1 , nshell
c           ****** symmetry check on i
         do 30 itr = 1 , ntd
            ish = iso(ii,itr)
            if (ish.gt.ii) go to 240
            mi(itr) = ish
 30      continue
         if (ii.le.ns) then
            lit = 1
            mini = 1
            maxi = 1
            loci = ii - 1
            go to 40
         end if
         if (ii.gt.nsnp) then
            lit = 3
            mini = 5
            maxi = 10
            loci = (ii-nsnp)*6 + ns + np3 - 10
            go to 40
         end if
         lit = 2
         mini = 2
         maxi = 4
         loci = (ii-ns)*3 + ns - 4
 40      do 230 jj = 1 , ii
c           ****** symmetry check on j
            do 50 itr = 1 , ntd
               indin(itr) = iso(jj,itr)
 50         continue
            do 60 itr = 1 , ntd
               jsh = indin(itr)
               if (jsh.gt.ii) go to 230
               ish = mi(itr)
               if (ish.lt.jsh) then
                  n = ish
                  ish = jsh
                  jsh = n
               end if
               if (ish.eq.ii .and. jsh.gt.jj) go to 230
 60         continue
            if (jj.le.ns) then
               ljt = 1
               minj = 1
               maxj = 1
               locj = jj - 1
               go to 70
            end if
            if (jj.gt.nsnp) then
               ljt = 3
               minj = 5
               maxj = 10
               locj = (jj-nsnp)*6 + ns + np3 - 10
               go to 70
            end if
            ljt = 2
            minj = 2
            maxj = 4
            locj = (jj-ns)*3 + ns - 4
 70         iandj = ii.eq.jj
            jmax = maxj
c     ----- find the equivalent blocks
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
            do 140 itr = 1 , ntd
               ntr = itr
               kk = mi(itr)
               ll = indin(itr)
               if (kk.le.ns) then
                  lock = kk - 1
                  go to 80
               end if
               if (kk.gt.nsnp) then
                  lock = (kk-nsnp)*6 + ns + np3 - 10
                  go to 80
               end if
               lock = (kk-ns)*3 + ns - 4
 80            if (ll.le.ns) then
                  locl = ll - 1
                  go to 90
               end if
               if (ll.gt.nsnp) then
                  locl = (ll-nsnp)*6 + ns + np3 - 10
                  go to 90
               end if
               locl = (ll-ns)*3 + ns - 4
 90            do 110 k = mini , maxi
                  lck = lock + k
                  if (iandj) jmax = k
                  do 100 l = minj , jmax
                     if (ll.gt.kk) then
                        kl = ia(locl+l) + lck
                     else
                        kl = ia(lck) + locl + l
                     end if
                     t(k,l) = f(kl)
                     if (iandj) t(l,k) = f(kl)
 100              continue
 110           continue
               if (lit.gt.1 .or. ljt.gt.1) call rhrd
               do 130 i = mini , maxi
                  lci = ia(loci+i) + locj
                  if (iandj) jmax = i
                  do 120 j = minj , jmax
                     ij = lci + j
                     h(ij) = h(ij) + t(i,j)
 120              continue
 130           continue
 140        continue
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
            do 220 itr = 2 , ntd
               kk = mi(itr)
               ll = indin(itr)
               if (kk.ge.ll) then
                  k = kk
                  l = ll
               else
                  k = ll
                  l = kk
               end if
               if (k.eq.ii .and. l.eq.jj) go to 220
               ntr = itr + 1
               if (ntr.le.ntd) then
                  do 150 it = ntr , ntd
                     i = mi(it)
                     j = indin(it)
                     if (i.lt.j) then
                        ij = i
                        i = j
                        j = ij
                     end if
                     if (i.eq.k .and. j.eq.l) go to 220
 150              continue
               end if
               ntr = invt(itr)
               do 170 i = mini , maxi
                  lci = ia(loci+i) + locj
                  if (iandj) jmax = i
                  do 160 j = minj , jmax
                     t(i,j) = h(lci+j)
                     if (iandj) t(j,i) = h(lci+j)
 160              continue
 170           continue
               if (lit.gt.1 .or. ljt.gt.1) call rhrd
               if (kk.le.ns) then
                  lock = kk - 1
                  go to 180
               end if
               if (kk.gt.nsnp) then
                  lock = (kk-nsnp)*6 + ns + np3 - 10
                  go to 180
               end if
               lock = (kk-ns)*3 + ns - 4
 180           if (ll.le.ns) then
                  locl = ll - 1
                  go to 190
               end if
               if (ll.gt.nsnp) then
                  locl = (ll-nsnp)*6 + ns + np3 - 10
                  go to 190
               end if
               locl = (ll-ns)*3 + ns - 4
 190           do 210 k = mini , maxi
                  lck = lock + k
                  if (iandj) jmax = k
                  do 200 l = minj , jmax
                     if (ll.gt.kk) then
                        kl = ia(locl+l) + lck
                     else
                        kl = ia(lck) + locl + l
                     end if
                     h(kl) = t(k,l)
 200              continue
 210           continue
 220        continue
 230     continue
 240  continue
      dum = one/dfloat(ntd)
      do 250 i = 1 , nx
         f(i) = h(i)*dum
 250  continue
      return
      end
      subroutine sytype(orbmo,eigval,sx,mp2fo,mp2fv,norbit,q,iwr)
c
c     --------------------------------------------
c     constructs symmetry data for 4-index package
c     project m.0.'s to give pure symmetry types
c     only works for simple groups (no degenerate reps.)
c     --------------------------------------------------
c
      implicit REAL  (a-h,o-z)
      logical prout
c
INCLUDE(common/sizes)
c
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nnnnt,nt3,
     1itable(8,8),mult(8,8),irr(maxorb),imosb(8,maxorb),irrb(maxorb)
     2,it2(8,maxorb),it(8,maxorb),ibb(maxorb),ins(8),iss(8),ivsf(2,8)
c
      dimension orbmo(num,num),sx(num*(num+1)/2),eigval(num)
     &,q(num*16)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nzc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),nat,
     & num,ns,np,nsp
     &,n1,nppx,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat)
     &,ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d
     &,idefs1(numspl),idefs2(numspl),idefp1(numspl),idefp2(numspl)
     &,idefd1(numspl),idefd2(numspl),ipars,iparp
     &,ipard,nsspl,npspl,ndspl
c
      character *8 groupy
      common/indsyx/groupy
      character *8 groupx
      common/symtrx/groupx
c
      common/molsym/trg(12),index,naxis
      common/scra  /iso(mxshel,48),nt
      dimension icol(maxorb),ntemp(maxorb)
      common/bufb/ptr(3,144),dtr(6,288)
INCLUDE(common/mapper)
INCLUDE(common/prints)
c
      logical uhf
      character *1 dash
c     character *8 open
      character *8 grp,grpy
      dimension grp(8),grpy(8),ichar(8),
     & char(8)
      data dash/'-'/
      data one/1.0d0/
c     data open/'open'/
      data grpy/'c1','ci','cs','c2','c2v','d2','c2h','d2h'/
      data grp/'c1','ci','cs','cn','cnv','dn','cnh','dnh'/
c
      nactiv = mp2fv - mp2fo + 1
      prout = nt.gt.1 .and. oprint(46)
      uhf = .false.
c     prout = .true.
c
c    -----------------------------------
c    check group and get character table
c    -----------------------------------
c
c     nx = num*(num+1)/2
c
      nshell = ns + np + nd
      do 20 i = 1 , 8
         if (groupx.eq.grp(i)) then
            igrp = i
            go to 40
         end if
 20   continue
 30   write (iwr,6010)
      call caserr(' use non-degenerate group')
 40   go to (60,60,60,50,50,50,50,50) , igrp
 50   if (naxis.ne.2) go to 30
 60   call indtab(igrp,itable,mult)
c     nirr = nt
      nnnnt = nt
      groupy = grpy(igrp)
      if (prout) then
         write (iwr,6020) groupy , (i,i=1,nt)
         ndash = 13 + 4*nt
         write (iwr,6030) (dash,i=1,ndash)
         do 70 i = 1 , nt
            write (iwr,6040) i , (itable(i,j),j=1,nt)
 70      continue
      end if
c
c     --------------------------------------------------
c     construct transformation table for basis functions
c     --------------------------------------------------
c
      do 80 i = 1 , num
         ibasis(1,i) = i
 80   continue
      if (nt.gt.1) then
         do 210 i = 1 , nshell
            kty = 1
c
            if (i.gt.ns) then
               kty = 2
               jstore = ns + (i-ns-1)*3
            end if
c
            if (i.gt.(ns+np)) then
               kty = 3
               jstore = ns + np*3 + (i-ns-np-1)*6
            end if
c
            do 200 iop = 2 , nt
               go to (90,100,150) , kty
c     s functions
 90            ibasis(iop,i) = iso(i,iop)
               go to 200
c     p functions
 100           npp = (iop-1)*3
               kstore = ns + (iso(i,iop)-ns-1)*3
               do 140 j = 1 , 3
                  do 110 k = 1 , 3
                     tr = ptr(k,npp+j)
                     if (dabs(tr-one).lt.1.0d-8) go to 120
                     if (dabs(tr+one).lt.1.0d-8) go to 130
 110              continue
                  call caserr(
     +             'error in sytype - wrong group/orientation')
 120              ibasis(iop,jstore+j) = kstore + k
                  go to 140
 130              ibasis(iop,jstore+j) = -(kstore+k)
 140           continue
               go to 200
c     d functions
 150           npp = (iop-1)*6
               kstore = ns + np*3 + (iso(i,iop)-ns-np-1)*6
               do 190 j = 1 , 6
                  do 160 k = 1 , 6
                     tr = dtr(k,npp+j)
                     if (dabs(tr-one).lt.1.0d-8) go to 170
                     if (dabs(tr+one).lt.1.0d-8) go to 180
 160              continue
                  call caserr(
     +           'error in sytype - wrong group/orientation')
 170              ibasis(iop,jstore+j) = kstore + k
                  go to 190
 180              ibasis(iop,jstore+j) = -(kstore+k)
 190           continue
 200        continue
 210     continue
      end if
      if (prout) then
         write (iwr,6070) (i,i=1,nt)
         write (iwr,6030) (dash,i=1,ndash)
         do 220 i = 1 , num
            write (iwr,6060) i , (ibasis(j,i),j=1,nt)
 220     continue
      end if
c
c
c     character of basis set
      do 240 iop = 1 , nt
         ichar(iop) = 0
         do 230 j = 1 , num
            ib = ibasis(iop,j)
            if (iabs(ib).eq.j) then
               if (ib.lt.0) ichar(iop) = ichar(iop) - 1
               if (ib.gt.0) ichar(iop) = ichar(iop) + 1
            end if
 230     continue
 240  continue
      if (prout) then
       write (iwr,*)' reducible representation of basis functions'
       write (iwr,*) (ichar(i),i=1,nt)
      endif
      issrep = 0
      do 270 irep = 1 , nt
         nirep = 0
         do 250 iop = 1 , nt
            nirep = nirep + ichar(iop)*itable(irep,iop)
 250     continue
         nirep = nirep/nt
         ins(irep) = nirep
         iss(irep) = issrep
         do 260 i = 1 , nirep
            j = i + issrep
            ibb(j) = i
 260     continue
         issrep = issrep + nirep
         if (prout) then
         write (iwr,*) ' representation ' , irep , ' occurs ' , 
     +     nirep ,' times'
         endif
 270  continue
c
c
      nc = 0
      do 290 i = 1 , nt
         do 280 j = 1 , num
            it2(i,j) = 0
            it(i,j) = 0
 280     continue
 290  continue
      do 390 ir = 1 , nt
c     representation ir
         do 380 i = 1 , num
c     loop over basis functions, eliminate non-unique
c     functions
            do 300 iop = 1 , nt
               j = ibasis(iop,i)
               jj = iabs(j)
               if (jj.gt.i) go to 380
 300        continue
            do 310 j = 1 , num
               icol(j) = 0
 310        continue
            do 320 iop = 1 , nt
               jj = ibasis(iop,i)
               j = iabs(jj)
               icol(j) = icol(j) + ksign(jj)*itable(ir,iop)
 320        continue
            do 330 j = 1 , num
               if (icol(j).ne.0) go to 340
 330        continue
            go to 380
 340        do 350 j = 1 , num
               jj = icol(j)
               icol(j) = ksign(jj)
 350        continue
            nc = nc + 1
            irr(nc) = ir
            if (nc.gt.num) call caserr('error2 in sytype')
            do 360 iop = 1 , nt
               imos(iop,nc) = itable(ir,iop)
 360        continue
            id = 0
            do 370 j = 1 , num
               jj = icol(j)
               if (jj.ne.0) then
                  id = id + 1
                  if (id.gt.nt) call caserr('error2 in sytype')
                  it(id,nc) = ksign(jj)*j
               end if
 370        continue
 380     continue
 390  continue
c
c**** rearrange the it array so it shows which symmetry adapted
c**** basis functions a particular b.f. adds to and what the symmetry
c**** is of that basis function
c
      do 410 i = 1 , num
         nirep = irr(i)
         do 400 nsym = 1 , nt
            ipos = it(nsym,i)
            if (ipos.ne.0) then
               ival = i
               if (ipos.lt.0) then
                  ival = -i
                  ipos = -ipos
               end if
               if (it2(nirep,ipos).ne.0)
     +             call caserr('mistake in sytype')
               it2(nirep,ipos) = ival
            end if
 400     continue
 410  continue
c
c     nsq = num*num
c     orbmo=mo's,eigval=eigenvalues,sx=overlap
c     q(i3),q(i4),q(i5)=decompositions of vectors
c
      i3 = 0
      i4 = 0 + num*nt
      nocc = ne/2
      if (nt.eq.1) then
         do 420 i = 1 , nactiv
            irr(i) = 1
            imos(1,i) = 1
 420     continue
         ivsf(2,1) = mp2fv - nocc
         ivsf(1,1) = 1
         go to 580
      end if
c
c     --------------------------------------
c     decide symmetry types of m.o.'s
c     also project out minor symmetry
c     contaminants to ensure exact symmetry
c     for each m.o.
c     --------------------------------------
c
c     ant = one/dfloat(nt)
c
      i = 1
c
c      ii=ilifq(i)
 430  call inddeg(eigval,norbit,i,nred)
      call indcha(orbmo(1,i),num,nred,q(i3+1),sx,char,ichar,iwr)
      call reducd(orbmo(1,i),num,nred,q(i3+1),sx,q(i4+1),ichar,i,uhf,
     +            iwr)
      i = i + nred
c
      if (i.le.norbit) go to 430
c
c***** rearrange the virtual orbitals according to symmetry type
c
      j = nocc
      do 450 nsym = 1 , nt
         ivsf(1,nsym) = j - nocc + 1
         do 440 i = nocc + 1 , mp2fv
            if (irr(i).eq.nsym) then
               j = j + 1
               icol(i) = j
            end if
 440     continue
         ivsf(2,nsym) = j - nocc
 450  continue
      do 460 i = nocc + 1 , mp2fv
         ntemp(i) = irr(i)
 460  continue
      do 470 i = nocc + 1 , mp2fv
         irr(icol(i)) = ntemp(i)
 470  continue
      do 500 k = 1 , nt
         do 480 i = nocc + 1 , mp2fv
            ntemp(i) = imos(k,i)
 480     continue
         do 490 i = nocc + 1 , mp2fv
            imos(k,icol(i)) = ntemp(i)
 490     continue
 500  continue
      do 510 i = nocc + 1 , mp2fv
         q(i) = eigval(i)
 510  continue
      do 520 i = nocc + 1 , mp2fv
         eigval(icol(i)) = q(i)
 520  continue
      do 550 k = 1 , num
         do 530 i = nocc + 1 , mp2fv
            q(i) = orbmo(k,i)
 530     continue
         do 540 i = nocc + 1 , mp2fv
            orbmo(k,icol(i)) = q(i)
 540     continue
 550  continue
c
      do 570 i = 1 , nactiv
         j = mp2fo + i - 1
         irr(i) = irr(j)
         do 560 k = 1 , nt
            imos(k,i) = imos(k,j)
 560     continue
 570  continue
c
c
c
 580  if (prout) then
         write (iwr,6080) (i,i=1,nt)
         ndash = 21 + 4*nt
         write (iwr,6030) (dash,i=1,ndash)
         do 590 i = 1 , nactiv
            write (iwr,6050) i , irr(i) , (imos(j,i),j=1,nt)
 590     continue
      end if
c
c
      return
 6010 format (//10x,'error in symmetry  program'/10x,
     +        'group is not one of - c1,c2,ci,cs,c2v,d2,c2h or d2h')
 6020 format (//10x,'character table -- group  ',a4/10x,30('-')///10x,
     +        'i.r.',4x,':',4x,8i4)
 6030 format (10x,60a1)
 6040 format (18x,':'/18x,':'/10x,i3,5x,':',4x,8i4)
 6050 format (26x,':'/10x,i3,5x,i3,5x,':',4x,8i4)
 6060 format (18x,':'/10x,i3,5x,':',4x,8i4)
 6070 format (//10x,'transformation table for basis functions'/10x,
     +        40('-')//10x,'a.o.',4x,':',4x,8i4)
 6080 format (//10x,'representation table for active',
     +        ' molecular orbitals'/10x,50('-')//10x,'m.o.',4x,'i.r.',
     +        4x,':',4x,8i4)
      end
      subroutine tdownd(qnew,ilifn,q,ilifq,nnn)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension qnew(*),q(*),ilifn(*),ilifq(*)
cfu      common/blkorbs/dum(maxorb)
      dimension dum(maxorb)
      common/trand/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +  ctran(mxorb3),otran
INCLUDE(common/infoa)
c
      if(otran)then
       do 60009 i=1,nnn
       call dcopy(num,q(ilifq(i)+1),1,qnew(ilifn(i)+1),1)
60009  continue
      else
       do 731 i=1,nnn
       m=ilifq(i)
       call vclr(dum,1,num)
       do 733 j=1,num
       n=ntran(j)
       do 733 k=1,n
       l=ilifc(j)    +k
733    dum(itran(l))=ctran(l)*q(m+j)+dum(itran(l))
       call dcopy(num,dum(1),1,qnew(ilifn(i)+1),1)
 731   continue
      endif
      return
      end
      subroutine tranpd(p1,p2,num)
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension p1(num,*),p2(num,*)
      common/blkin/ct1(mxorb3),it1(mxorb3)
      common/trand/ilifc(maxorb),nterm(maxorb),it(mxorb3),ct(mxorb3)
     * ,otran
c
      l3=num*num
      if(otran)goto 100
      do 1 i=1,num
      nt1=nterm(i)
      do 2 k=1,nt1
      l=k+ilifc(i)
      ct1(k)=ct(l)
2     it1(k)=it(l)
      do 1 j=1,i
      p2(i,j)=0.0d0
      nt2=nterm(j)
      do 1 k=1,nt2
      kk=k+ilifc(j)
      l=it(kk)
      do 1 ii=1,nt1
      ll=it1(ii)
      top=p1(max(l,ll),min(ll,l))
1     p2(i,j)=ct1(ii)*ct(kk)*top+p2(i,j)
      return
 100  continue
      call dcopy(l3,p1(1,1),1,p2(1,1),1)
      return
      end
      subroutine tst11(codens,dij,dx,tmax1,test2,xklt,igath,zmem
     &,ic,jc,iscat,kmin,kmax
     &,l1,l2,lmax,lmin,nkl,ipsize,ext2,jcount,jc0,jc1)
c
c**** the data contained in common used by the direct program is
c**** explained in the routine start ***************************
c****
c****** rotuines performs magnitude-based integral test in same
c****** way as ahlrichs & haser.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,4),igath(nkl),xklt(nkl),codens(nb,nb)
     &,iscat(nkl)
c
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext6
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      jc0 = 0
      jc1 = 0
      dx = dij
      do 90 kc = kmin , kmax
         kpmax = nc(kc)
         d1 = max(dij,codens(kc,ic),codens(kc,jc))
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
         if (ext2) lmax1 = lmax1 - 1
         jcold = jc0 + 1
         jc1new = jc1 + 1
         jcold2 = jcount + 1
         do 20 lc = lmin1 , lmax1
            jc1 = jc1 + 1
            jcount = jcount + kpmax*nc(lc)
c........... density matrix weighting to test
            zmem(jc1,4) = max(d1,codens(lc,ic),codens(lc,jc),
     +                    codens(lc,kc))
            dx = max(dx,zmem(jc1,4))
 20      continue
         do 30 lc = jc1new , jc1
            if (zmem(lc,4).lt.1.0d-40) zmem(lc,4) = 1.0d-40
 30      continue
c
         do 40 lc = jc1new , jc1
            zmem(lc,4) = test2/zmem(lc,4)
 40      continue
c
         ngath = jcount - jcold2 + 1
         call dgthr(ngath,zmem(1,4),zmem(jcold2,2),iscat(jcold2))
c
         do 50 icount = jcold2 , jcount
c...................the test itself
            if (xklt(icount).gt.zmem(icount,2)) then
               jc0 = jc0 + 1
               igath(jc0) = icount
            end if
 50      continue
         if (ext2) then
            lc = lmax1 + 1
            jc1 = jc1 + 1
            densm = max(d1,codens(lc,ic),codens(lc,jc),codens(lc,kc))
            dx = max(dx,densm)
            if (densm.lt.1.0d-20) then
               test = tmax1
            else
               test = test2/densm
            end if
            ext6 = ext2 .and. lc.eq.kc
            lpmax = nc(lc)
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
            do 70 kprm = 1 , kpmax
               if (ext6) lpmax = kprm
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
               do 60 lprm = 1 , lpmax
                  jcount = jcount + 1
                  if (xklt(jcount).gt.test) then
                     jc0 = jc0 + 1
                     igath(jc0) = jcount
                  end if
 60            continue
 70         continue
         end if
c
c** need to store the ci coordinates
         do 80 jc01 = jcold , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 80      continue
 90   continue
      return
      end
      subroutine tst1dr(codens,dij,dx,tmax1,test2,xklt,igath,zmem
     &,ic,jc,iscat,kmin,kmax
     &,l1,l2,lmax,lmin,nkl,ipsize,ext2,jcount,jc0,jc1)
c
c
c****** magnitude based integral test for derivative integrals.
c****** see cwm thesis for explanation of test.
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,4),igath(nkl),xklt(nkl),codens(nb,nb)
     &,iscat(nkl)
c
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext6
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      jc0 = 0
      jc1 = 0
      dx = 0.0d0
      do 90 kc = kmin , kmax
         kpmax = nc(kc)
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
         if (ext2) lmax1 = lmax1 - 1
         jcold = jc0 + 1
         jc1new = jc1 + 1
         jcold2 = jcount + 1
         do 20 lc = lmin1 , lmax1
            jc1 = jc1 + 1
            jcount = jcount + kpmax*nc(lc)
            zmem(jc1,4) = dij*codens(lc,kc) + codens(lc,jc)
     +                    *codens(kc,ic) + codens(lc,ic)*codens(kc,jc)
            dx = max(dx,zmem(jc1,4))
 20      continue
         do 30 lc = jc1new , jc1
            if (zmem(lc,4).lt.1.0d-40) zmem(lc,4) = 1.0d-40
 30      continue
c
         do 40 lc = jc1new , jc1
            zmem(lc,4) = test2/zmem(lc,4)
 40      continue
c
         ngath = jcount - jcold2 + 1
         call dgthr(ngath,zmem(1,4),zmem(jcold2,2),iscat(jcold2))
c
         do 50 icount = jcold2 , jcount
            if (xklt(icount).gt.zmem(icount,2)) then
               jc0 = jc0 + 1
               igath(jc0) = icount
            end if
 50      continue
         if (ext2) then
            lc = lmax1 + 1
            jc1 = jc1 + 1
            densm = dij*codens(lc,kc) + codens(lc,jc)*codens(kc,ic)
     +              + codens(lc,ic)*codens(kc,jc)
            dx = max(dx,densm)
            if (densm.lt.1.0d-20) then
               test = tmax1
            else
               test = test2/densm
            end if
            ext6 = ext2 .and. lc.eq.kc
            lpmax = nc(lc)
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
            do 70 kprm = 1 , kpmax
               if (ext6) lpmax = kprm
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
               do 60 lprm = 1 , lpmax
                  jcount = jcount + 1
                  if (xklt(jcount).gt.test) then
                     jc0 = jc0 + 1
                     igath(jc0) = jcount
                  end if
 60            continue
 70         continue
         end if
c
c** need to store the ci coordinates
         do 80 jc01 = jcold , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 80      continue
 90   continue
      return
      end
_EXTRACT(tst1s,itanium,linux64)
      subroutine tst1s(codens,dij,dx,test2,xklt
     &,igath,zmem,ic,jc,iscat,kmin,kmax,m1,m2,m3
     &,l1,l2,lmax,lmin,nkl,ipsize,ext2,ext3,jcount,jc0,jc1)
c
c***** routine misses out small integrals and symmetry equivalent integr
c***** not vectorised. for large molecules can take an enormous amount o
c***** time. the partially vectorised magnitude based test is performed
c***** before the symmetry check on the l shell. if this is not done the
c***** the symmetry checking can take 50% of the calculation  time. for
c***** derivatives and mp2 testing the lshell symmetry check is done pri
c***** to the magnitude based test and should similar time problems be
c***** encountered there, this routine contains the blueprint of how
c***** to deal with it.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,6),igath(nkl),iscat(nkl),codens(nb,nb)
     &,m1(48),m2(48),m3(48),xklt(nkl)
c
      common/scra  /iso(mxshel,48),nt
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext3
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      xnt = dfloat(nt)
      jc0 = 0
      jc1 = 0
      dx = dij
      do 160 kc = kmin , kmax
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
c***** now have the symmetry check on the k loop
c
c        ext8 = .false.
c
         if (ext3) then
            do 20 it = 2 , nt
               kds = iso(kc,it)
               if (kds.gt.ic) then
                  call kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
                  go to 160
               end if
               m3(it) = kds
 20         continue
         else
            do 30 it = 2 , nt
               kds = iso(kc,it)
               m3(it) = kds
               if (m1(it).eq.ic) then
                  if (m2(it).eq.jc) then
                     if (kds.gt.kc) then
                        call kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
                        go to 160
                     end if
                  end if
               end if
 30         continue
         end if
         kpmax1 = nc(kc)*(nc(kc)+1)/2
         kpmax2 = nc(kc)*(nc(kc)-1)/2
         jcold = jc0 + 1
         jc1new = jc1 + 1
         jcxx = jcount
         jcou1 = jcount + 1
         d1 = max(dij,codens(ic,kc),codens(jc,kc))
         do 40 lc = lmin1 , lmax1
            jc1 = jc1 + 1
            jcxx = jcxx + nc(kc)*nc(lc)
            zmem(jc1,5) = max(d1,codens(lc,ic),codens(lc,jc),
     +                    codens(lc,kc))
            dx = max(dx,zmem(jc1,5))
 40      continue
         do 50 lc = jc1new , jc1
            if (zmem(lc,5).lt.1.0d-30) zmem(lc,5) = 1.0d-30
 50      continue
c
         if (lmax1.eq.kc) jcxx = jcxx - kpmax2
         ngath = jcxx - jcount
         call dgthr(ngath,zmem(1,5),zmem(jcou1,6),iscat(jcou1))
c
         do 60 lc = jcou1 , jcxx
            zmem(lc,6) = xklt(lc)*zmem(lc,6)
 60      continue
c
         jcxx = jcount
         do 140 lc = lmin1 , lmax1
            kpmax0 = nc(kc)*nc(lc)
            if (ext2 .and. lc.eq.kc) kpmax0 = kpmax1
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
            jcxx = jcount
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
            do 70 kprm = 1 , kpmax0
               jcxx = jcxx + 1
               if (zmem(jcxx,6).gt.test2) go to 80
 70         continue
            jcount = jcount + kpmax0
            go to 140
c
c****** symetry check on lshell
c
 80         xn4 = 1.0d0
            if (ext3) then
               if (ext2) then
                  do 90 it = 2 , nt
                     lds = iso(lc,it)
                     if (lds.gt.ic) then
                        jcount = jcount + kpmax0
                        go to 140
                     end if
                     kds = m3(it)
                     if (kds.lt.lds) then
                        nds = kds
                        kds = lds
                        lds = nds
                     end if
                     ids = m1(it)
                     if (ids.eq.ic .or. kds.eq.ic) then
                        jds = m2(it)
                        if (ids.le.kds) then
                           if (ids.ne.kds .or. jds.lt.lds) then
                              nds = ids
                              ids = kds
                              kds = nds
                              nds = jds
                              jds = lds
                              lds = nds
                           end if
                        end if
                        if (jds.ge.jc) then
                           if (jds.gt.jc) then
                              jcount = jcount + kpmax0
                              go to 140
                           end if
                           if (kds.ge.kc) then
                              if (kds.gt.kc) then
                                 jcount = jcount + kpmax0
                                 go to 140
                              end if
                              if (lds.ge.lc) then
                                 if (lds.gt.lc) then
                                    jcount = jcount + kpmax0
                                    go to 140
                                 end if
                                 xn4 = xn4 + 1.0d0
                              end if
                           end if
                        end if
                     end if
 90               continue
               else
                  do 100 it = 2 , nt
                     lds = iso(lc,it)
                     ids = m1(it)
                     kds = m3(it)
                     if (ids.eq.ic .or. kds.eq.ic) then
                        jds = m2(it)
                        if (ids.le.kds) then
                           if (ids.ne.kds .or. jds.lt.lds) then
                              nds = ids
                              ids = kds
                              kds = nds
                              nds = jds
                              jds = lds
                              lds = nds
                           end if
                        end if
                        if (jds.ge.jc) then
                           if (jds.gt.jc) then
                              jcount = jcount + kpmax0
                              go to 140
                           end if
                           if (kds.ge.kc) then
                              if (kds.gt.kc) then
                                 jcount = jcount + kpmax0
                                 go to 140
                              end if
                              if (lds.ge.lc) then
                                 if (lds.gt.lc) then
                                    jcount = jcount + kpmax0
                                    go to 140
                                 end if
                                 xn4 = xn4 + 1.0d0
                              end if
                           end if
                        end if
                     end if
 100              continue
               end if
            else if (ext2) then
               do 110 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        lds = iso(lc,it)
                        if (kds.lt.lds) then
                           nds = kds
                           kds = lds
                           lds = nds
                           if (kds.gt.kc) then
                              jcount = jcount + kpmax0
                              go to 140
                           end if
                        end if
                        if (kds.eq.kc) then
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 jcount = jcount + kpmax0
                                 go to 140
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 110           continue
            else
               do 120 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        if (kds.eq.kc) then
                           lds = iso(lc,it)
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 jcount = jcount + kpmax0
                                 go to 140
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 120           continue
            end if
            xmult = xnt/xn4
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
            do 130 kprm = 1 , kpmax0
               jcount = jcount + 1
               jc0 = jc0 + 1
               igath(jc0) = jcount
               zmem(jc0,4) = xmult
 130        continue
c
 140     continue
c** need to store the ci coordinates
         do 150 jc01 = jcold , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 150     continue
 160  continue
      return
      end
_ENDEXTRACT
      subroutine tst1sd(codens,dij,dx,tmax1,test2,xklt
     &,igath,zmem,ic,jc,kmin,kmax,m1,m2,m3,xklk
     &,l1,l2,lmax,lmin,jstart,nkl,ipsize,ext2,ext3,jcount,jc0,jc1)
c
c****** symmetry + magnitude testing for derivative integrals.
c****** if expensive, then the remarks in tst1s should be considered.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,6),igath(nkl),xklt(nkl),codens(nb,nb)
     &,m1(48),m2(48),m3(48),xklk(nkl)
c
      common/scra  /iso(mxshel,48),nt
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext3,ext6
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      xnt = dfloat(nt)
      jc0 = 0
      jc1 = 0
      dx = 0.0d0
      do 110 kc = kmin , kmax
         jcold = jc0
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
c***** now have the symmetry check on the k loop
c
c        ext8 = .false.
c
         if (ext3) then
            do 20 it = 2 , nt
               kds = iso(kc,it)
               if (kds.gt.ic) then
                  call kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
                  go to 110
               end if
               m3(it) = kds
 20         continue
         else
            do 30 it = 2 , nt
               kds = iso(kc,it)
               m3(it) = kds
               if (m1(it).eq.ic) then
                  if (m2(it).eq.jc) then
                     if (kds.gt.kc) then
                        call kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
                        go to 110
                     end if
                  end if
               end if
 30         continue
         end if
         do 90 lc = lmin1 , lmax1
            ext6 = ext2 .and. lc.eq.kc
            kpmax = nc(kc)
            lpmax = nc(lc)
            jc1 = jc1 + 1
c
c****** symetry check on lshell
c
            xn4 = 1.0d0
            if (ext3) then
               if (ext2) then
                  do 40 it = 2 , nt
                     lds = iso(lc,it)
                     if (lds.gt.ic) then
                        call lmiss(kpmax,lpmax,ext6,jcount)
                        go to 90
                     end if
                     kds = m3(it)
                     if (kds.lt.lds) then
                        nds = kds
                        kds = lds
                        lds = nds
                     end if
                     ids = m1(it)
                     if (ids.eq.ic .or. kds.eq.ic) then
                        jds = m2(it)
                        if (ids.le.kds) then
                           if (ids.ne.kds .or. jds.lt.lds) then
                              nds = ids
                              ids = kds
                              kds = nds
                              nds = jds
                              jds = lds
                              lds = nds
                           end if
                        end if
                        if (jds.ge.jc) then
                           if (jds.gt.jc) then
                              call lmiss(kpmax,lpmax,ext6,jcount)
                              go to 90
                           end if
                           if (kds.ge.kc) then
                              if (kds.gt.kc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 90
                              end if
                              if (lds.ge.lc) then
                                 if (lds.gt.lc) then
                                    call lmiss(kpmax,lpmax,ext6,jcount)
                                    go to 90
                                 end if
                                 xn4 = xn4 + 1.0d0
                              end if
                           end if
                        end if
                     end if
 40               continue
               else
                  do 50 it = 2 , nt
                     lds = iso(lc,it)
                     ids = m1(it)
                     kds = m3(it)
                     if (ids.eq.ic .or. kds.eq.ic) then
                        jds = m2(it)
                        if (ids.le.kds) then
                           if (ids.ne.kds .or. jds.lt.lds) then
                              nds = ids
                              ids = kds
                              kds = nds
                              nds = jds
                              jds = lds
                              lds = nds
                           end if
                        end if
                        if (jds.ge.jc) then
                           if (jds.gt.jc) then
                              call lmiss(kpmax,lpmax,ext6,jcount)
                              go to 90
                           end if
                           if (kds.ge.kc) then
                              if (kds.gt.kc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 90
                              end if
                              if (lds.ge.lc) then
                                 if (lds.gt.lc) then
                                    call lmiss(kpmax,lpmax,ext6,jcount)
                                    go to 90
                                 end if
                                 xn4 = xn4 + 1.0d0
                              end if
                           end if
                        end if
                     end if
 50               continue
               end if
            else if (ext2) then
               do 60 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        lds = iso(lc,it)
                        if (kds.lt.lds) then
                           nds = kds
                           kds = lds
                           lds = nds
                           if (kds.gt.kc) then
                              call lmiss(kpmax,lpmax,ext6,jcount)
                              go to 90
                           end if
                        end if
                        if (kds.eq.kc) then
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 90
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 60            continue
            else
               do 70 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        if (kds.eq.kc) then
                           lds = iso(lc,it)
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 90
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 70            continue
            end if
            densm = dij*codens(lc,kc) + codens(lc,jc)*codens(kc,ic)
     +              + codens(lc,ic)*codens(kc,jc)
            dx = max(densm,dx)
            if (densm.lt.1.0d-20) then
               test = tmax1
            else
               test = test2/densm
            end if
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
            klpmax = kpmax*lpmax
            if (ext6) klpmax = kpmax*(kpmax+1)/2
            xmult = xnt/xn4
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
            do 80 kprm = 1 , klpmax
               jcount = jcount + 1
c..............magnitude test
               if (xklt(jcount+jstart).gt.test) then
                  jc0 = jc0 + 1
                  igath(jc0) = jcount
                  zmem(jc0,6) = xklk(jstart+jcount)*xmult
               end if
 80         continue
c
 90      continue
c** need to store the ci coordinates
         do 100 jc01 = jcold + 1 , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 100     continue
 110  continue
      return
      end
      subroutine tst1sm(test2,xklt
     &,igath,zmem,ic,jc,kmin,kmax,m1,m2,m3,xklk
     &,l1,l2,lmax,lmin,jstart,nkl,ipsize,ext2,jcount,jc0,jc1)
c
c***** symmetry test + magnitude test on integrals for an mp2 calculatio
c***** if this routine, proves very expensive in second wave optimisatio
c***** then the remarks at the top of tst1s should be considered. the
c***** ordinary dupuis and king symmetry procedure has been slightly
c***** amended since the ij>kl permutational symmtery is not exploited.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,6),igath(nkl),xklt(nkl)
     &,m1(48),m2(48),m3(48),xklk(nkl)
c
      common/scra  /iso(mxshel,48),nt
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext6
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      xnt = dfloat(nt)
      jc0 = 0
      jc1 = 0
      do 80 kc = kmin , kmax
         jcold = jc0
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
c***** now have the symmetry check on the k loop
c
c        ext8 = .false.
c
         do 20 it = 2 , nt
            kds = iso(kc,it)
            m3(it) = kds
            if (m1(it).eq.ic) then
               if (m2(it).eq.jc) then
                  if (kds.gt.kc) then
                     call kmiss(lmin1,lmax1,kc,ext2,jcount,jc1)
                     go to 80
                  end if
               end if
            end if
 20      continue
c
         do 60 lc = lmin1 , lmax1
            ext6 = ext2 .and. lc.eq.kc
            kpmax = nc(kc)
            lpmax = nc(lc)
            jc1 = jc1 + 1
c
c****** symetry check on lshell
c
            xn4 = 1.0d0
            if (ext2) then
               do 30 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        lds = iso(lc,it)
                        if (kds.lt.lds) then
                           nds = kds
                           kds = lds
                           lds = nds
                           if (kds.gt.kc) then
                              call lmiss(kpmax,lpmax,ext6,jcount)
                              go to 60
                           end if
                        end if
                        if (kds.eq.kc) then
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 60
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 30            continue
            else
               do 40 it = 2 , nt
                  ids = m1(it)
                  if (ids.eq.ic) then
                     jds = m2(it)
                     if (jds.eq.jc) then
                        kds = m3(it)
                        if (kds.eq.kc) then
                           lds = iso(lc,it)
                           if (lds.ge.lc) then
                              if (lds.gt.lc) then
                                 call lmiss(kpmax,lpmax,ext6,jcount)
                                 go to 60
                              end if
                              xn4 = xn4 + 1.0d0
                           end if
                        end if
                     end if
                  end if
 40            continue
            end if
c
c**** the k and l primitive loops for a given kc,lc
c**** pairing. these loops are inner integral test
c
            klpmax = kpmax*lpmax
            if (ext6) klpmax = kpmax*(kpmax+1)/2
            xmult = xnt/xn4
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
            do 50 kprm = 1 , klpmax
               jcount = jcount + 1
c................ magnitude based test
               if (xklt(jcount+jstart).gt.test2) then
                  jc0 = jc0 + 1
                  igath(jc0) = jcount
                  zmem(jc0,6) = xklk(jstart+jcount)*xmult
               end if
 50         continue
c
 60      continue
c** need to store the ci coordinates
         do 70 jc01 = jcold + 1 , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 70      continue
 80   continue
      return
      end
      subroutine tst22(zmem,zmemc,iscatt,jc0,jc1,memcon)
c
c**** this routine fully contracts the half contracted ordinary
c**** or derivative two electron integrals. it also simultaneously
c**** performs a scatter operation. the routine is poorly vectorised
c**** since the vector loop involves paging over a lot of memory
c**** and it takes about 10-15% of the time for the two electron
c**** integrals and about twice that for the derivatives. the basic
c**** problem is that we haven't split our kl batches so that
c**** all pairings in batches have the same degree of contraction, but
c**** have opted for a more general approach allowing extremely long vec
c**** lengths. it is not clear whether this was the right or the wrong
c**** decision.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension zmem(ipsize,mempto),zmemc(icsize,memcto),iscatt(jc0)
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      common/cider/intgrl(8,100),intcon(5,12),ihz(5,5)
     &,ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri
     &,mempto,memcto,mtemp,ipsz,icsz
     &,ipsize,icsize
c
      do 20 m = 1 , memcon
         call vclr(zmemc(1,m),1,jc1)
 20   continue
c
c***** now loop over all the indices associated with the
c***** stored vector loop
c
      do 40 jcount = 1 , jc0
         do 30 m = 1 , memcon
            zmemc(iscatt(jcount),m) = zmemc(iscatt(jcount),m)
     +                                + zmem(jcount,m+mempri)
 30      continue
 40   continue
c
      return
      end
      subroutine tstmp2(test2,xklt,igath,zmem,kmin,kmax
     &,l1,l2,lmax,lmin,jstart,nkl,ipsize,ext2,jcount,jc0,jc1)
c
c****** magnitude based integral test for mp2 calculation. see cwm
c****** thesis for an explanation of how it works.
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/data/charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb)
     &,isatbf(maxat),isatp(maxat),ipatbf(maxat),ipatp(maxat),
     & nwshp(mxprim)
     &,nup(maxorb),pe(mxprim),pc(mxprim),coprm(mxprim,3),na,nb,ns,np,nsp
     &,n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8
     &,n9,n10,n11,n12,n13,n14,n15,n16,n17
c
      dimension zmem(ipsize,3),igath(nkl),xklt(nkl)
c
      common/defunk/cobas(maxorb,3),nd
      logical ext2
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      jc0 = 0
      jc1 = 0
      do 50 kc = kmin , kmax
         kpmax = nc(kc)
c*** a few logicals to determine extents of l loop
         lmin1 = l1
         lmax1 = l2
         if (ext2) lmax1 = kc
         if (kc.eq.kmax) lmax1 = lmax
         if (kc.eq.kmin) lmin1 = lmin
         if (ext2) lmax1 = lmax1 - 1
         jc1 = jc1 + lmax1 - lmin1 + 1
         jcold2 = jcount
         do 20 lc = lmin1 , lmax1
            jcount = jcount + kpmax*nc(lc)
 20      continue
         if (ext2) then
            jcount = jcount + kpmax*(kpmax+1)/2
            jc1 = jc1 + 1
         end if
         jcold = jc0
         do 30 icount = jcold2 + 1 , jcount
            if (xklt(icount+jstart).gt.test2) then
               jc0 = jc0 + 1
               igath(jc0) = icount
            end if
 30      continue
c
c** need to store the ci coordinates
         do 40 jc01 = jcold + 1 , jc0
            zmem(jc01,1) = cobas(kc,1)
            zmem(jc01,2) = cobas(kc,2)
            zmem(jc01,3) = cobas(kc,3)
 40      continue
 50   continue
c..... here!!
      return
      end
      subroutine orfogd(q,qp,b,c,iky,ilifq,newbas,nrow)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),qp(*),b(*),c(*),iky(*),ilifq(*)
      common/blkin/p(1)
c... to orthogonalize the cols of q - result to qp
c... overlap matrix supplied in b - destroyed on exit
c... scratch triangle in c
c... qp can overwrite q
      top=dsqrt(1.0d0/b(1))
      b(1)=top
      do i=2,newbas
        m=iky(i)
        c(m+1)=b(m+1)*top
      enddo
      do k=2,newbas
      m=iky(k)
      top=b(m+k)
      ll=k-1
      do i=1,ll
        bot=c(m+i)
        p(i)=-bot
        top=top-bot*bot
      enddo
      top=dsqrt(1.0d0/top)
      b(m+k)=top
      n=k+1
      if(n.le.newbas) then
        do l=n,newbas
        nsp=iky(l)
        c(nsp+k)=(ddot(ll,p,1,c(nsp+1),1)  +b(nsp+k))*top
        enddo
      endif
      do l=1,ll
       bot=0.0d0
       do  i=l,ll
        bot=p(i)*b(l+iky(i))+bot
       enddo
       b(l+m)=bot*top
      enddo
      enddo
c
      do i=1,nrow
_IF1(c)      call gather(newbas,p,q(i+1),ilifq)
_IFN1(c)      call dgthr(newbas,q(i+1),p,ilifq)
      m=1
      do j=1,newbas
       qp(i+ilifq(j))=ddot(j,b(m),1,p,1)
       m=m+j
      enddo
      enddo
c
      return
      end
      subroutine ver_direct(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/direct.m,v $
     +     "/
      data revision /"$Revision: 6291 $"/
      data date /"$Date: 2013-12-05 14:38:55 +0100 (Thu, 05 Dec 2013) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
