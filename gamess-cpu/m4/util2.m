c 
c  $Author: jvl $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util2.m,v $
c  $State: Exp $
c  
c     deck=util2
c ******************************************************
c ******************************************************
c             =   util2  =
c ******************************************************
c ******************************************************
      subroutine blocki
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/pkfil)
INCLUDE(common/iofile)
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/files)
      common/blkin/gout(510),mword
INCLUDE(common/discc)
INCLUDE(common/restar)
INCLUDE(common/shlt)
INCLUDE(common/atmblk)
INCLUDE(common/maxlen)
      if(omaxb) return 
      mword=icount-1
      if(mword.eq.0)go to 41
      ochek = .false.
c
c...    for incore scf go on to see how much space is needed ...
      if(imaxb_ic.eq.1) go to 50
c
      if(nopk.ne.1) then
      if(.not.opank)then
_IFN1(iv)         call pack(gout(num2ep+1),lab1632,integ,numlabp)
_IF1(iv)         call pak4v(intij,gout(num2ep+1))
      else
_IFN1(iv)         call pack(gout(num2ejk+num2ejk+1),lab1632,
_IFN1(iv)     +             integ,numlabjk)
_IF1(iv)         call pak6v(intij,gout(num2ejk+num2ejk+1))
      endif
c
      else
c
      labss = lab816 + lab816
      numlabi = numlab / 2
_IFN1(iv)      call pack(gout(num2e+1),labss,integ,numlabi)
_IF1(iv)      call pak4v(intij,gout(num2e+1))
      endif
      if(iposun(mainp).ne.iblkmp)callsearch(iblkmp,mainp)
      call put(gout(1),m511,mainp)
c ======================================================
c = debug print for integrals (nprint=4)
c ======================================================
      if(.not.out) go to 240
      if(nopk.eq.1) go to 230
      if(opank) then
      call pr2ejk(mainp,iblkmp)
      else
      call pr2ep (mainp,iblkmp)
      endif
      go to 240
 230  call pr2ei (mainp,iblkmp)
 240  continue
c ============================================================
 50   icount=1
_IFN1(iv)      ic4 = 1
      nrec=nrec+1
      iblkmp=iblkmp+1
      mblp=mblp+1
      if(mblp)41,40,41
 40    m2last(m2file)=iblkmp
       mfilep=mfilep+1
      if(mfilep.le.n2file)go to 1
      if (omem(mainp)) then
         imaxb_ic = 1
c...    memory segment is the only one ??
         mfilep = mfilep - 1
      else
         omaxb = .true.
      end if
c
      if (omaxb) then
_IF(parallel)
      print 3,' ** problem detected on node ',ipg_nodeid(),' **'
3     format(//,a29,i5,a3)
      write(6,2)yed(mainp)
      call flushn(6)
      call caserr('MAINFILE section(s) exceeded')
_ELSE
      write(iwr,2)yed(mainp)
_ENDIF
 2    format(//
     *' *** sorry ***'//
     *' size of mainfile section ',a4,' has been exceeded'//
     *' restart with additional section allocated'//)
      end if
c
      go to 41
  1   mainp=n2tape(mfilep)
      iblkmp=n2blk(mfilep)
      mblp=iblkmp-n2last(mfilep)
      m2file=m2file+1
      m2tape(m2file)=mainp
      m2blk(m2file)=iblkmp
      m2last(m2file)=-1
  41  return
      end
      subroutine qout(g)
c
c ..  qout is now restricted to conventional scf
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/integ(340),intkl(340)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
      common/blkin/gout(510),nword
INCLUDE(common/indez)
      dimension g(*)
c
c
c     ----- pack the 4 indices of two integrals into one word
c     ----- write label + integral on mainfile
c
      ijn = 0
      jmax = maxj
      do 2600 i = mini,maxi
      if (oianj) jmax = i
      i1 = loci + i
      ipack = i4096(i1)
      do 2400 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      i2 = locj + j
      if(i1-i2)10, 20, 20
 10   lab1 = i4096(i2) + i1
      go to 30
 20   lab1 = ipack + i2
 30   lmax = maxl
      kln = 0
      do 2200 k = mink,maxk
      if (okanl) lmax = k
      i3 = lock + k
      kpack = i4096(i3)
      do 2000 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 2400
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 2000
      i4 = locl + l
      if(i3-i4) 40, 50, 50
 40   lab2 = i4096(i4) + i3
      go to 60
 50   lab2 = kpack + i4
 60   gout(icount) = val
_IFN1(iv)      integ(ic4  )=max(lab1,lab2)
_IFN1(iv)      integ(ic4+1)=min(lab1,lab2)
_IFN1(iv)      ic4=ic4+2
_IF1(iv)      integ(icount)=max(lab1,lab2)
_IF1(iv)      intkl(icount)=min(lab1,lab2)
      icount = icount+1
      if (icount .le. nintmx) go to 2000
      call blocki
      if(omaxb)go to 261
 2000 continue
 2200 continue
 2400 continue
 2600 continue
 261  return
      end
_IF(unicos)
_IF(ccpdft)
      subroutine qoutd(fock,dmat,g,
     +    fac1, fac2, facex, ocoul, oexch, odft)
_ELSE
      subroutine qoutd(fock,dmat,g)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
      common/craypk/integ(680),ipad(680)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
      common/blkin/gout(510),nword
INCLUDE(common/indez)
      dimension g(*),fock(*),dmat(*)
_IF(ccpdft)
      if (.not. odft) then
_ENDIF
c
c     ----- pack the 4 indices of two integrals into one word
c     ----- unicos specific to use CAL dbuild
c
      ijn = 0
      jmax = maxj
      ifac = 4096
      do 260 i = mini,maxi
      if (oianj) jmax = i
      i1 = loci + i
      ipack = ifac*i1
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      i2 = locj + j
      if(i1-i2)1, 2, 2
    1 lab1 = ifac*i2 + i1
      go to 3
    2 lab1 = ipack + i2
    3 continue
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      i3 = lock + k
      kpack = ifac*i3
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
         gout(icount) = val
         i4 = locl + l
         if(i3-i4) 4, 5, 5
    4    lab2 = ifac*i4 + i3
         go to 6
    5    lab2 = kpack + i4
    6    integ(ic4  )=max(lab1,lab2)
         integ(ic4+1)=min(lab1,lab2)
         ic4=ic4+2
         icount = icount+1
         if (icount .le. nintmx) go to 200
          ochek=.false.
          nrec=nrec+1
          call dbuild(fock,dmat)
  200 continue
  220 continue
  240 continue
  260 continue
_IF(ccpdft)
      else
c
c     non conventional fock build involved .. cannot use CAL
c
      ijn = 0
      jmax = maxj
      do 360 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 340 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 320 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 300 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 340
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 300
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
           facij=fac1
           fackl=fac2
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
           facij=fac2
           fackl=fac1
        endif
c
c Coulomb term
c
        if(ocoul)then
           itr12=iky(i1)+i2
           itr34=iky(i3)+i4
           val2=val+val
           val4=val2+val2
           fock(itr12) = facij*val4*dmat(itr34) + fock(itr12)
           if(itr12 .ne. itr34)then
              fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
           endif
        endif
c
c exchange term
c
        if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
         itr13=iky(i1)+i3
         itr14=iky(i1)+i4
         itr23=iky(max(i2,i3))+min(i2,i3)
         itr24=iky(max(i2,i4))+min(i2,i4)
         val=val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2
         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
  300 continue
  320 continue
  340 continue
  360 continue
c
      endif
_ENDIF
      return
      end
_ELSE
_IF(ccpdft)
      subroutine dbuild(fock,dmat,g,
     &  fac1, fac2, facex, ocoul, oexch)
_ELSE
      subroutine dbuild(fock,dmat,g)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)

      dimension g(*),fock(*),dmat(*)
c
c     ----- DSCF fock builder
c
      ijn = 0
      jmax = maxj

      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
_IF(ccpdft)
           facij=fac1
           fackl=fac2
_ENDIF
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
_IF(ccpdft)
           facij=fac2
           fackl=fac1
_ENDIF
        endif
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
_IF(ccpdft)
c
c Coulomb term
c
        if(ocoul)then
           val2=val+val
           val4=val2+val2
           fock(itr12) = facij*val4*dmat(itr34) + fock(itr12)
           if(itr12 .ne. itr34)then
              fock(itr34) = fackl*val4*dmat(itr12) + fock(itr34)
           endif
        endif
c
c exchange
c
        if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
         val = val*facex
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2

         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
      endif
_ELSE
c
c this is the HF-only version
c
        val2=val+val
        val4=val2+val2
        val13=val
        val14=val
        if(i1.eq.i3 .or. i2.eq.i4) val13=val2
        if(i2.eq.i3) val14=val2
        f12 = val4*dmat(itr34) + fock(itr12)
        fock(itr34) = val4*dmat(itr12) + fock(itr34)
        fock(itr12) = f12
        f23 = fock(itr23) - val14*dmat(itr14)
        f14 = fock(itr14) - val14*dmat(itr23)
        f13 = fock(itr13) - val13*dmat(itr24)
        fock(itr24) = fock(itr24) - val13*dmat(itr13)
        fock(itr23) = f23
        fock(itr14) = f14
        fock(itr13) = f13
_ENDIF
  200 continue
  220 continue
  240 continue
  260 continue
      end
_ENDIF

_IF(ccpdft)
      subroutine dir_build_uhf(fock,ak,p,q,g,
     &   fac1, fac2, facex, ocoul, oexch)
_ELSE
      subroutine dir_build_uhf(fock,ak,p,q,g)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      dimension g(*),fock(*),ak(*),p(*),q(*)
c
c     ----- Direct open-shell UHF Fock builder 
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
_IF(ccpdft)
           facij=fac1
           fackl=fac2
_ENDIF
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
_IF(ccpdft)
           facij=fac2
           fackl=fac1
_ENDIF
        endif
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
_IF(ccpdft)
c
c Coulomb term
c
        if(ocoul)then

           gik=val
           val2=gik+gik
           val4=val2+val2

           fock(itr12) = facij*val4*p(itr34) + fock(itr12)
           if(itr12 .ne. itr34)then
              fock(itr34) = fackl*val4*p(itr12) + fock(itr34)
           endif

        endif
c
c exchange
c
        if(oexch)then
c
c full term for HF case or weighted term for b3lyp etc
c
           gik=val*facex
           val2=gik+gik
           val4=val2+val2
           gil=gik
           if(i1.eq.i3.or.i2.eq.i4)gik=val2
           if(i2.eq.i3)gil=val2
           if(i2.ge.i3)goto 1
           itr23=ikyk+i2
           if(i2.ge.i4)goto 1
           itr24=iky(i4)+i2
 1         ajk=fock(itr23)-gil*p(itr14)
           bjk=ak(itr23)+gil*q(itr14)
           ail=fock(itr14)-gil*p(itr23)
           bil=ak(itr14)+gil*q(itr23)
           aik=fock(itr13)-gik*p(itr24)
           bik=ak(itr13)+gik*q(itr24)
           fock(itr24)=fock(itr24)-gik*p(itr13)
           ak(itr24)=ak(itr24)+gik*q(itr13)
           fock(itr23)=ajk
           ak(itr23)=bjk
           fock(itr14)=ail
           ak(itr14)=bil
           fock(itr13)=aik
           ak(itr13)=bik

      endif
_ELSE
c
c  non-dft UHF version
c
      gik=val
      val2=gik+gik
      val4=val2+val2
      aij=val4*p(itr34)+fock(itr12)
      fock(itr34)=val4*p(itr12)+fock(itr34)
      fock(itr12)=aij
c... exchange
      gil=gik
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)goto 1
      itr23=ikyk+i2
      if(i2.ge.i4)goto 1
      itr24=iky(i4)+i2
1     ajk=fock(itr23)-gil*p(itr14)
      bjk=ak(itr23)+gil*q(itr14)
      ail=fock(itr14)-gil*p(itr23)
      bil=ak(itr14)+gil*q(itr23)
      aik=fock(itr13)-gik*p(itr24)
      bik=ak(itr13)+gik*q(itr24)
      fock(itr24)=fock(itr24)-gik*p(itr13)
      ak(itr24)=ak(itr24)+gik*q(itr13)
      fock(itr23)=ajk
      ak(itr23)=bjk
      fock(itr14)=ail
      ak(itr14)=bil
      fock(itr13)=aik
      ak(itr13)=bik
_ENDIF
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_IF()
      subroutine dir_build_uhf_orig(fock,ak,p,q,g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      dimension g(*),fock(*),ak(*),p(*),q(*)
c
c     ----- Direct open-shell UHF Fock builder 
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=gik+gik
      val4=val2+val2
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      aij=val4*p(itr34)+fock(itr12)
      fock(itr34)=val4*p(itr12)+fock(itr34)
      fock(itr12)=aij
c... exchange
      gil=gik
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)goto 1
      itr23=ikyk+i2
      if(i2.ge.i4)goto 1
      itr24=iky(i4)+i2
1     ajk=fock(itr23)-gil*p(itr14)
      bjk=ak(itr23)+gil*q(itr14)
      ail=fock(itr14)-gil*p(itr23)
      bil=ak(itr14)+gil*q(itr23)
      aik=fock(itr13)-gik*p(itr24)
      bik=ak(itr13)+gik*q(itr24)
      fock(itr24)=fock(itr24)-gik*p(itr13)
      ak(itr24)=ak(itr24)+gik*q(itr13)
      fock(itr23)=ajk
      ak(itr23)=bjk
      fock(itr14)=ail
      ak(itr14)=bil
      fock(itr13)=aik
      ak(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_ENDIF
      subroutine dir_build_open(coul,exch,dens,g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      dimension g(*),coul(*),exch(*),dens(*)
c
c     ----- Direct open-shell SCF J & K builder (nshell=1)
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        gik=val
        val2=val+val
        ikyi=iky(i1)
        ikyj=iky(i2)
        ikyk=iky(i3)
        itr13=ikyi+i3
        itr14=ikyi+i4
        itr12=ikyi+i2
        itr23=ikyj+i3
        itr24=ikyj+i4
        itr34=ikyk+i4
        gil=val
        if(i1.eq.i3.or.i2.eq.i4)gik=val2
        if(i2.eq.i3)gil=val2
        if(i2.ge.i3)goto 280
        itr23=ikyk+i2
        if(i2.ge.i4)goto 280
        itr24=iky(i4)+i2
  280   bij=val2*dens(itr34)+coul(itr12)
        coul(itr34)=val2*dens(itr12)+coul(itr34)
        coul(itr12)=bij
        bjk=exch(itr23)+gil*dens(itr14)
        bil=exch(itr14)+gil*dens(itr23)
        bik=exch(itr13)+gik*dens(itr24)
        exch(itr24)=exch(itr24)+gik*dens(itr13)
        exch(itr23)=bjk
        exch(itr14)=bil
        exch(itr13)=bik
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
      subroutine dir_build_open2(l2,coul,exch,dens,g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/ijlab)
INCLUDE(common/cslosc)
INCLUDE(common/mapper)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      dimension g(*),coul(*),exch(*),dens(*)
c
c     ----- Direct open-shell SCF J & K builder (nshell > 1)
c
      ijn = 0
      jmax = maxj
      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
      gik=val
      val2=val+val
      ikyi=iky(i1)
      ikyj=iky(i2)
      ikyk=iky(i3)
      itr13=ikyi+i3
      itr14=ikyi+i4
      itr12=ikyi+i2
      itr23=ikyj+i3
      itr24=ikyj+i4
      itr34=ikyk+i4
      gil=val
      if(i1.eq.i3.or.i2.eq.i4)gik=val2
      if(i2.eq.i3)gil=val2
      if(i2.ge.i3)go to 280
      itr23=ikyk+i2
      if(i2.ge.i4)go to 280
      itr24=iky(i4)+i2
  280 continue
_IF1(x)c$dir scalar
      do 300 iiii=1,nsheld
      bij=val2*dens(itr34)+coul(itr12)
      coul(itr34)=val2*dens(itr12)+coul(itr34)
      coul(itr12)=bij
      bjk=exch(itr23)+gil*dens(itr14)
      bil=exch(itr14)+gil*dens(itr23)
      bik=exch(itr13)+gik*dens(itr24)
      exch(itr24)=exch(itr24)+gik*dens(itr13)
      exch(itr23)=bjk
      exch(itr14)=bil
      exch(itr13)=bik
      itr12=itr12+l2
      itr34=itr34+l2
      itr13=itr13+l2
      itr14=itr14+l2
      itr23=itr23+l2
      itr24=itr24+l2
  300 continue
  200 continue
  220 continue
  240 continue
  260 continue
      return
      end
_IF(vector)
      subroutine pkfile(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,q,g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
      dimension q(*),g(*)
INCLUDE(common/sortp)
c ***
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/infoa)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
      common/incrs/kln2(9),lendd,ioff(8),ngm,iadgt
INCLUDE(common/misc)
INCLUDE(common/indez)
      common/blkin/gout(510),nword
INCLUDE(common/pkfil)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      common/flips/ispp(12),ib(mxprms),jb(mxprms),kb(mxprms)
      dimension lb(mxprms)
_IF1(iv)      dimension intpkl(204)
_IF1(iv)      equivalence (intpkl(1),intij(205))
      data dzero,pt5 /0.0d0,0.5d0/
      data lb/0,1,2,3,4,5,6,7,8,9/
c
      istar=lendd*lenint(4)+lenint(nx)+1
c
      do 10 i=1,mxprms
      kb(i)=lb(i)*ngm
      jb(i)=kb(i)*ngm
10    ib(i)=jb(i)*ngm
c
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
c
c     type = 1 for (ii ii ii ii)
c            2     (ii jj jj jj)
c            3     (ii ii kk kk)
c            4     (ii ii ii ll)
c            5     (ii jj kk kk)
c            6     (ii jj jj ll)
c            7     (ii ii kk ll)
c            8     (ii jj kk ll)
c
      if (ii-jj) 100, 100, 240
  100 if (jj-kk) 120, 120, 180
  120 if (kk-ll) 140, 140, 160
  140 ntyp = 1
      go to 380
  160 ntyp = 4
      go to 380
  180 if (kk-ll) 200, 200, 220
  200 ntyp = 3
      go to 380
  220 ntyp = 7
      go to 380
  240 if (jj-kk) 260, 260,320
  260 if (kk-ll) 280,280,300
  280 ntyp = 2
      go to 380
  300 ntyp = 6
      go to 380
  320 if (kk-ll) 340,340,360
  340 ntyp = 5
      go to 380
  360 ntyp = 8
  380 continue
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = istar
      norg2 = norg1+lendd
      norg3 = norg2+lendd
      if (oskpa .and. .not. onpsym) norg2 = istar
      if (oskpc .and. .not. onpsym) norg3 = istar + lendd
      if (oskpb .and. .not. onpsym) norg3 = istar
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      do 1060 i = mini,maxi
      if (oianj) jmax = i
      do 1040 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 1020 k = mink,kmaxx
      if (okanl) lmax = k
      do 1000 l = minl,lmax
      ia = i-mini+1
      ja = j-minj+1
      ka = k-mink+1
      la = l-minl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1 = ib(ia)+jb(ja)+kb(ka)+lb(la)+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2 = ib(ia)+jb(ka)+kb(ja)+lb(la)+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      go to (400,420,400,440,420,420,440,420),ntyp
  400 if (ia .eq. ja) go to 440
  420 n3 = ib(ia)+jb(la)+kb(ja)+lb(ka)+norg3
      go to 460
  440 n3 = ib(ja)+jb(ka)+kb(ia)+lb(la)+norg3
  460 continue
c
c     ----- form first linear combination -----
c
      jump = 1
      i1 = loci+i
      i2 = locj+j
      i3 = lock+k
      i4 = locl+l
      if (i2 .eq. i3) jump = 2
      if ((i2 .eq. i4) .or. (i1 .eq. i3)) jump = 3
      g1 = g(n1)
      g2 = g(n2)
      g3 = g(n3)
      go to (480,500,520,540),ind
  480 valk = g2+g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  500 valk = g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  520 valk = g2
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  540 valk = dzero
      valp = (g1+g1)+(g1+g1)
  560 continue
      g(n1) = valp
c     nn = n1
      go to 820
c
c     ----- form second linear combination -----
c
  580 go to (600,620,640,660),ind
  600 valk = g3+g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  620 valk = g1+g3
      valp = -valk
      go to 680
  640 valk = g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  660 valk = g1
      valp = -valk
  680 continue
      g(n2) = valp
c     nn = n2
      jump = 2
      if ((i1 .eq. i2) .or. (i3 .eq. i4)) jump = 3
      n = i2
      i2 = i3
      i3 = n
      go to 820
c
c     ----- form third linear combination -----
c
  700 go to (720,740,760,780),ind
  720 valk = g1+g2
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  740 valk = g1
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  760 valk = g1+g2
      valp = -valk
      go to 800
  780 valk = g1
      valp = -valk
  800 continue
      g(n3) = valp
c     nn = n3
      i2 = locl+l
      i3 = locj+j
      i4 = lock+k
      jump = 3
  820 continue
c
c     ----- store integral and indices -----
c
      if (opank) go to 880
c
c     ----- -p- supermatrix only. -----
c
      if ( dabs(valp) .lt. cutoff) go to 980
c     valp0 = valp
      if (i1 .eq. i3 .and. i2 .eq. i4) valp = valp*pt5
_IFN1(iv)      integ(ic4  )=iky(i1)+i2
_IFN1(iv)      integ(ic4+1)=iky(i3)+i4
_IFN1(iv)      ic4=ic4+2
_IF1(iv)      intij(icount)=iky(i1)+i2
_IF1(iv)      intkl(icount)=iky(i3)+i4
      gout(icount) = valp
      icount = icount+1
      if (icount .le. nintmx) go to 980
c ***
      if (osortp) then
        call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
      else
        call blocki
      endif
c ***
      if(omaxb)go to 1061
      go to 980
  880 continue
c
c     ----- -p- and -k- supermatrices -----
c
      if (dabs(valp).lt.cutoff .and. dabs(valk).lt.cutoff) go to 980
c     valp0 = valp
c     valk0 = valk
      if (i1 .ne. i3 .or. i2 .ne. i4) go to 900
      valp = valp*pt5
      valk = valk*pt5
  900 continue
_IFN1(iv)      integ(ic4  )=iky(i1)+i2
_IFN1(iv)      integ(ic4+1)=iky(i3)+i4
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1)+i2
_IF1(iv)      intpkl(icount)=iky(i3)+i4
      gout(icount) = valp
      gout(num2ejk+icount) = valk
      icount = icount+1
      if (icount .le. nintmx) go to 980
      call blocki
      if(omaxb)go to 1061
  980 go to (580,700,1000),jump
 1000 continue
 1020 continue
 1040 continue
 1060 continue
 1061 return
      end
_ELSE
      subroutine pkfile(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,
     +                   g,q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
      dimension q(*),g(*)
INCLUDE(common/sortp)
c ***
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
INCLUDE(common/indez)
      common/blkin/gout(510),nword
INCLUDE(common/pkfil)
INCLUDE(common/mapper)
INCLUDE(common/flips)
INCLUDE(common/atmblk)
_IF1(iv)      dimension intpkl(204)
_IF1(iv)      equivalence (intpkl(1),intij(205))
      data dzero,pt5 /0.0d0,0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
c
c     type = 1 for (ii ii ii ii)
c            2     (ii jj jj jj)
c            3     (ii ii kk kk)
c            4     (ii ii ii ll)
c            5     (ii jj kk kk)
c            6     (ii jj jj ll)
c            7     (ii ii kk ll)
c            8     (ii jj kk ll)
c
      if (ii-jj) 100, 100, 240
  100 if (jj-kk) 120, 120, 180
  120 if (kk-ll) 140, 140, 160
  140 ntyp = 1
      go to 380
  160 ntyp = 4
      go to 380
  180 if (kk-ll) 200, 200, 220
  200 ntyp = 3
      go to 380
  220 ntyp = 7
      go to 380
  240 if (jj-kk) 260, 260,320
  260 if (kk-ll) 280,280,300
  280 ntyp = 2
      go to 380
  300 ntyp = 6
      go to 380
  320 if (kk-ll) 340,340,360
  340 ntyp = 5
      go to 380
  360 ntyp = 8
  380 continue
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 50626
      norg3 = 101251
      if (oskpa .and. .not. onpsym) norg2 = 1
      if (oskpc .and. .not. onpsym) norg3 = 50626
      if (oskpb .and. .not. onpsym) norg3 = 1
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      do 1060 i = mini,maxi
      if (oianj) jmax = i
      do 1040 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 1020 k = mink,kmaxx
      if (okanl) lmax = k
      do 1000 l = minl,lmax
      ia = i-mini+1
      ja = j-minj+1
      ka = k-mink+1
      la = l-minl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1 = ibt(ia)+jbt(ja)+kbt(ka)+lbt(la)+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2 = ibt(ia)+jbt(ka)+kbt(ja)+lbt(la)+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      go to (400,420,400,440,420,420,440,420),ntyp
  400 if (ia .eq. ja) go to 440
  420 n3 = ibt(ia)+jbt(la)+kbt(ja)+lbt(ka)+norg3
      go to 460
  440 n3 = ibt(ja)+jbt(ka)+kbt(ia)+lbt(la)+norg3
  460 continue
c
c     ----- form first linear combination -----
c
      jump = 1
      i1 = loci+i
      i2 = locj+j
      i3 = lock+k
      i4 = locl+l
      if (i2 .eq. i3) jump = 2
      if ((i2 .eq. i4) .or. (i1 .eq. i3)) jump = 3
      g1 = g(n1)
      g2 = g(n2)
      g3 = g(n3)
      go to (480,500,520,540),ind
  480 valk = g2+g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  500 valk = g3
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  520 valk = g2
      valp = (g1+g1)+(g1+g1)-valk
      go to 560
  540 valk = dzero
      valp = (g1+g1)+(g1+g1)
  560 continue
      g(n1) = valp
c     nn = n1
      go to 820
c
c     ----- form second linear combination -----
c
  580 go to (600,620,640,660),ind
  600 valk = g3+g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  620 valk = g1+g3
      valp = -valk
      go to 680
  640 valk = g1
      valp = (g2+g2)+(g2+g2)-valk
      go to 680
  660 valk = g1
      valp = -valk
  680 continue
      g(n2) = valp
c     nn = n2
      jump = 2
      if ((i1 .eq. i2) .or. (i3 .eq. i4)) jump = 3
      n = i2
      i2 = i3
      i3 = n
      go to 820
c
c     ----- form third linear combination -----
c
  700 go to (720,740,760,780),ind
  720 valk = g1+g2
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  740 valk = g1
      valp = (g3+g3)+(g3+g3)-valk
      go to 800
  760 valk = g1+g2
      valp = -valk
      go to 800
  780 valk = g1
      valp = -valk
  800 continue
      g(n3) = valp
c     nn = n3
      i2 = locl+l
      i3 = locj+j
      i4 = lock+k
      jump = 3
  820 continue
c
c     ----- store integral and indices -----
c
      if (opank) go to 880
c
c     ----- -p- supermatrix only. -----
c
      if ( dabs(valp) .lt. cutoff) go to 980
c     valp0 = valp
      if (i1 .eq. i3 .and. i2 .eq. i4) valp = valp*pt5
_IFN1(iv)      integ(ic4  )=iky(i1)+i2
_IFN1(iv)      integ(ic4+1)=iky(i3)+i4
_IFN1(iv)      ic4=ic4+2
_IF1(iv)      intij(icount)=iky(i1)+i2
_IF1(iv)      intkl(icount)=iky(i3)+i4
      gout(icount) = valp
      icount = icount+1
      if (icount .le. nintmx) go to 980
c ***
      if (osortp) then
        call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
      else
        call blocki
      endif
c ***
      if(omaxb)go to 1061
      go to 980
  880 continue
c
c     ----- -p- and -k- supermatrices -----
c
      if (dabs(valp).lt.cutoff .and. dabs(valk).lt.cutoff) go to 980
c     valp0 = valp
c     valk0 = valk
      if (i1 .ne. i3 .or. i2 .ne. i4) go to 900
      valp = valp*pt5
      valk = valk*pt5
  900 continue
_IFN1(iv)      integ(ic4  )=iky(i1)+i2
_IFN1(iv)      integ(ic4+1)=iky(i3)+i4
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1)+i2
_IF1(iv)      intpkl(icount)=iky(i3)+i4
      gout(icount) = valp
      gout(num2ejk+icount) = valk
      icount = icount+1
      if (icount .le. nintmx) go to 980
      call blocki
      if(omaxb)go to 1061
  980 go to (580,700,1000),jump
 1000 continue
 1020 continue
 1040 continue
 1060 continue
 1061 return
      end
_ENDIF
      subroutine pfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,
     +                  g,q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
      dimension g(*),q(*)
INCLUDE(common/sortp)
c ***
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
      dimension jump(256),i1(256),i2(256),i3(256),i4(256)
     *,n1(256),n2(256),n3(256)
     *,o12(256),o34(256)
     *,valp1(256),i2out(256),i3out(256)
     *,valp2(256),valp3(256)
      common/blkin/goutx(510),nword
INCLUDE(common/pkfil)
INCLUDE(common/mapper)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 257
      norg3 = 513
      if ( .not. oskpa) go to 100
      if ( .not. onpsym) norg2 = 1
      ib2 = ib1
      jb2 = jb1
      kb2 = kb1
      lb2 = lb1
  100 if ( .not. oskpb) go to 140
      if ( .not. onpsym) norg3 = norg1
      if (ii .ne. kk) go to 120
      if (jj .eq. ll) go to 120
      ib3 = kb1
      jb3 = lb1
      kb3 = jb1
      lb3 = ib1
      go to 180
  120 ib3 = ib1
      jb3 = jb1
      kb3 = lb1
      lb3 = kb1
      go to 180
  140 continue
      if ( .not. oskpc) go to 180
      if ( .not. onpsym) norg3 = norg2
      if (ii .ne. jj) go to 160
      if (kk .eq. ll) go to 160
      ib3 = kb2
      jb3 = lb2
      kb3 = ib2
      lb3 = jb2
      go to 180
  160 ib3 = ib2
      jb3 = jb2
      kb3 = kb2
      lb3 = lb2
  180 continue
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      mmijkl=0
      do 780 i = mini,maxi
      if (oianj) jmax = i
      do 780 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 780 k = mink,kmaxx
      if (okanl) lmax = k
      do 780 l = minl,lmax
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(mmijkl) = ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2(mmijkl) = ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n3(mmijkl) = ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
      do 781 i=1,mmijkl
      jump(i) = 1
      if (i2(i) .eq. i3(i)) jump(i) = 2
      if (i2(i) .eq. i4(i) .or. i1(i) .eq. i3(i)) jump(i) = 3
      o34(i)=i3(i).eq.i4(i)
      o12(i)=i1(i).eq.i2(i)
 781   continue
c
      loop2=0
      loop3=0
c
      go to (200,220,240,260),ind
c
200   do 280 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      g3=g(n3(i))
      valp=4.0d0*g1-g2-g3
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      if(jump(i)-2)2001,2002,280
 2001 valp=4.0d0*g2-g1-g3
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 280
 2002 valp=4.0d0*g3-g1-g2
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 280  continue
      go to 3000
c
220   do 380 i=1,mmijkl
      g1=g(n1(i))
      g3=g(n3(i))
      valp=4.0d0*g1-g3
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      if(jump(i)-2)3001,3002,380
 3001 valp=-g1-g3
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 380
 3002 valp=4.0d0*g3-g1
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 380  continue
      go to 3000
c
240   do 480 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      valp=4.0d0*g1-g2
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      if(jump(i)-2)4001,4002,480
 4001 valp=4.0d0*g2-g1
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 480
 4002 valp=-g1-g2
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 480  continue
      go to 3000
c
260   do 580 i=1,mmijkl
      g1=g(n1(i))
      valp=4.0d0*g1
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      if(jump(i)-2)5001,5002,580
 5001 valp=-g1
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 580
 5002 valp=-g1
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      i3out(loop3)=i
 580  continue
c
c     ----- store integral and indices -----
c
c     ----- -p- supermatrix only. -----
c
3000  do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff) go to 560
      if (i1(i) .eq. i3(i) .and. i2(i) .eq. i4(i))
     * valp1(i)=valp1(i)*pt5
_IFN1(iv)      integ(ic4  )=iky(i1(i))+i2(i)
_IFN1(iv)      integ(ic4+1)=iky(i3(i))+i4(i)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(i))+i2(i)
_IF1(iv)      intkl(icount)=iky(i3(i))+i4(i)
      goutx(icount) = valp1(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 560   continue
c
      do 570 i=1,loop2
      if ( dabs(valp2(i)) .lt. cutoff) go to 570
      iii=i2out(i)
      if (i1(iii) .eq. i2(iii) .and. i3(iii) .eq. i4(iii))
     * valp2(i)=valp2(i)*pt5
_IFN1(iv)      integ(ic4  )=iky(i1(iii))+i3(iii)
_IFN1(iv)      integ(ic4+1)=iky(i2(iii))+i4(iii)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(iii))+i3(iii)
_IF1(iv)      intkl(icount)=iky(i2(iii))+i4(iii)
      goutx(icount) = valp2(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 570   continue
c
      do 590 i=1,loop3
      if ( dabs(valp3(i)) .lt. cutoff) go to 590
      iii=i3out(i)
      if (i1(iii) .eq. i2(iii) .and. i4(iii) .eq. i3(iii))
     * valp3(i)=valp3(i)*pt5
_IFN1(iv)      integ(ic4  )=iky(i1(iii))+i4(iii)
_IFN1(iv)      integ(ic4+1)=iky(i2(iii))+i3(iii)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(iii))+i4(iii)
_IF1(iv)      intkl(icount)=iky(i2(iii))+i3(iii)
      goutx(icount) = valp3(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 590   continue
  782 return
      end
      subroutine pfi70s(ii,jj,kk,ll,g,q)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension g(*),q(*)
INCLUDE(common/sortp)
c ***
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
      dimension i1(256),i2(256),i3(256),i4(256),
     *            n1(3,256),valpg(3,256),
     *            valp1(256),valp2(256),valp3(256)
      common/blkin/goutx(510),nword
INCLUDE(common/pkfil)
INCLUDE(common/mapper)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c     data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
c     ind = 1
      norg1 = 1
      norg2 = 257
      norg3 = 513
      mmijkl=0
      do 780 i = mini,maxi
      do 780 j = minj,maxj
      do 780 k = mink,maxk
      do 780 l = minl,maxl
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(1,mmijkl)=ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n1(2,mmijkl)=ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n1(3,mmijkl)=ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
_IFN1(cx)      call dgthr(mmijkl*3,g,valpg,n1)
_IF1(c)      call gather(mmijkl*3,valpg,g,n1)
c
_IF1(i)cc*vdir: ignore recrdeps
_IFN(convex)
      do 280 i=1,mmijkl
      valp1(i)=4.0d0*valpg(1,i)-valpg(2,i)-valpg(3,i)
      valp2(i)=4.0d0*valpg(2,i)-valpg(1,i)-valpg(3,i)
      valp3(i)=4.0d0*valpg(3,i)-valpg(1,i)-valpg(2,i)
 280  continue
_ENDIF
_IF(convex)
      do 280 i=1,mmijkl
      valp1(i)=4.0d0*g(n1(1,i))-g(n1(2,i))-g(n1(3,i))
      valp2(i)=4.0d0*g(n1(2,i))-g(n1(1,i))-g(n1(3,i))
      valp3(i)=4.0d0*g(n1(3,i))-g(n1(1,i))-g(n1(2,i))
 280  continue
_ENDIF
c
c     ----- store integral and indices -----
c
c     ----- -p- supermatrix only. -----
c
      do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff) go to 560
_IFN1(iv)      integ(ic4  )=iky(i1(i))+i2(i)
_IFN1(iv)      integ(ic4+1)=iky(i3(i))+i4(i)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(i))+i2(i)
_IF1(iv)      intkl(icount)=iky(i3(i))+i4(i)
      goutx(icount) = valp1(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 560   continue
c
      do 570 i=1,mmijkl
      if ( dabs(valp2(i)) .lt. cutoff) go to 570
_IFN1(iv)      integ(ic4  )=iky(i1(i))+i3(i)
_IFN1(iv)      integ(ic4+1)=iky(i2(i))+i4(i)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(i))+i3(i)
_IF1(iv)      intkl(icount)=iky(i2(i))+i4(i)
      goutx(icount) = valp2(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 570   continue
c
      do 590 i=1,mmijkl
      if ( dabs(valp3(i)) .lt. cutoff) go to 590
_IFN1(iv)      integ(ic4  )=iky(i1(i))+i4(i)
_IFN1(iv)      integ(ic4+1)=iky(i2(i))+i3(i)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(i))+i4(i)
_IF1(iv)      intkl(icount)=iky(i2(i))+i3(i)
      goutx(icount) = valp3(i)
      icount = icount+1
      if (icount .gt. nintmx) then
        if (osortp) then
          call zsortp(q(igbf),q(iklbf),q(iiptbf),ijbase)
         else
          call blocki
        endif
        if(omaxb)go to 782
      endif
 590   continue
  782 return
      end
      subroutine pkfi70(ii,jj,kk,ll,oskpa,oskpb,oskpc,onpsym,g)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension g(*)
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/intij(204),intkl(204)
INCLUDE(common/shlt)
INCLUDE(common/restar)
INCLUDE(common/nshel)
INCLUDE(common/shlnos)
INCLUDE(common/misc)
      dimension jump(256),i1(256),i2(256),i3(256),i4(256)
     *, n1(256),n2(256),n3(256)
     *, o12(256),o34(256)
     *, valp1(256),i2out(256),i3out(256)
     *, valp2(256),valp3(256),valk1(256),valk2(256),valk3(256)
      common/blkin/gout(510),nword
INCLUDE(common/pkfil)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data pt5/0.5d0/
      lit = ktype(ii)
      mini = kmin(ii)
      maxi = kmax(ii)
      loci = kloc(ii)-mini
      minj = kmin(jj)
      maxj = kmax(jj)
      locj = kloc(jj)-minj
      lkt = ktype(kk)
      mink = kmin(kk)
      maxk = kmax(kk)
      lock = kloc(kk)-mink
      minl = kmin(ll)
      maxl = kmax(ll)
      locl = kloc(ll)-minl
      oianj = ii .eq. jj
      okanl = kk .eq. ll
      oident = (ii .eq. kk) .and. (jj .eq. ll)
c
c     -----
c
      ind = 1
      if (oskpa .and. onpsym) ind = 2
      if (oskpb .and. onpsym) ind = 3
      if (oskpc .and. onpsym) ind = 3
      if (oskpa .and. oskpb .and. onpsym) ind = 4
      norg1 = 1
      norg2 = 257
      norg3 = 513
      if ( .not. oskpa) go to 100
      if ( .not. onpsym) norg2 = 1
      ib2 = ib1
      jb2 = jb1
      kb2 = kb1
      lb2 = lb1
  100 if ( .not. oskpb) go to 140
      if ( .not. onpsym) norg3 = norg1
      if (ii .ne. kk) go to 120
      if (jj .eq. ll) go to 120
      ib3 = kb1
      jb3 = lb1
      kb3 = jb1
      lb3 = ib1
      go to 180
  120 ib3 = ib1
      jb3 = jb1
      kb3 = lb1
      lb3 = kb1
      go to 180
  140 continue
      if ( .not. oskpc) go to 180
      if ( .not. onpsym) norg3 = norg2
      if (ii .ne. jj) go to 160
      if (kk .eq. ll) go to 160
      ib3 = kb2
      jb3 = lb2
      kb3 = ib2
      lb3 = jb2
      go to 180
  160 ib3 = ib2
      jb3 = jb2
      kb3 = kb2
      lb3 = lb2
  180 continue
      jmax = maxj
      kmaxx = maxk
      lmax = maxl
      mmijkl=0
      do 780 i = mini,maxi
      if (oianj) jmax = i
      do 780 j = minj,jmax
      if (jj .eq. kk) kmaxx = j
      do 780 k = mink,kmaxx
      if (okanl) lmax = k
      do 780 l = minl,lmax
      mmijkl=mmijkl+1
c
c     ----- calculate n1 for (i,j//k,l) -----
c
      n1(mmijkl) = ib(ib1,i )+ib(jb1,j )+ib(kb1,k )+ib(lb1,l )+norg1
c
c     ----- calculate n2 for (i,k//j,l) -----
c
      n2(mmijkl) = ib(ib2,i )+ib(jb2,k )+ib(kb2,j )+ib(lb2,l )+norg2
c
c     ----- calculate n3 for (i,l//j,k) -----
c
      n3(mmijkl) = ib(ib3,i )+ib(jb3,l )+ib(kb3,j )+ib(lb3,k )+norg3
c
c     ----- form first linear combination -----
c
      i1(mmijkl) = loci+i
      i2(mmijkl) = locj+j
      i3(mmijkl) = lock+k
      i4(mmijkl) = locl+l
 780  continue
c
      do 781 i=1,mmijkl
      jump(i) = 1
      if (i2(i) .eq. i3(i)) jump(i) = 2
      if (i2(i) .eq. i4(i) .or. i1(i) .eq. i3(i)) jump(i) = 3
      o34(i)=i3(i).eq.i4(i)
      o12(i)=i1(i).eq.i2(i)
 781   continue
c
      loop2=0
      loop3=0
c
      go to (200,220,240,260),ind
c
200   do 280 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      g3=g(n3(i))
      valk=g2+g3
      valp=4.0d0*g1-valk
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)2001,2002,280
 2001 valk=g1+g3
      valp=4.0d0*g2-valk
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 280
 2002 valk=g1+g2
      valp=4.0d0*g3-valk
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 280  continue
      go to 3000
c
220   do 380 i=1,mmijkl
      g1=g(n1(i))
      g3=g(n3(i))
      valk=g3
      valp=4.0d0*g1-valk
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)3001,3002,380
 3001 valk=g1+g3
      valp=-valk
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 380
 3002 valk=g1
      valp=4.0d0*g3-valk
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 380  continue
      go to 3000
c
240   do 480 i=1,mmijkl
      g1=g(n1(i))
      g2=g(n2(i))
      valk=g2
      valp=4.0d0*g1-valk
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      valk1(i)=valk
      if(jump(i)-2)4001,4002,480
 4001 valk=g1
      valp=4.0d0*g2-valk
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 480
 4002 valk=g1+g2
      valp=-valk
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 480  continue
      go to 3000
c
260   do 580 i=1,mmijkl
      g1=g(n1(i))
      valp=4.0d0*g1
_IF1()      g(n1(i))=valp
      valp1(i)=valp
      valk1(i)=0.0d0
      if(jump(i)-2)5001,5002,580
 5001 valk=g1
      valp=-valk
_IF1()      g(n2(i))=valp
      loop2=loop2+1
      valp2(loop2)=valp
      valk2(loop2)=valk
      i2out(loop2)=i
      if(o12(i).or.o34(i))go to 580
 5002 valk=g1
      valp=-valk
_IF1()      g(n3(i))=valp
      loop3=loop3+1
      valp3(loop3)=valp
      valk3(loop3)=valk
      i3out(loop3)=i
 580  continue
c
c     ----- store integral and indices -----
c
c     ----- -p+k- supermatrix only. -----
c
3000  do 560 i=1,mmijkl
      if ( dabs(valp1(i)) .lt. cutoff.and.  dabs(valk1(i)).lt.cutoff)
     * go to 560
      if (i1(i) .ne. i3(i) .or. i2(i) .ne. i4(i))
     * go to 5601
      valp1(i)=valp1(i)*pt5
      valk1(i)=valk1(i)*pt5
5601  continue
_IFN1(iv)      integ(ic4  )=iky(i1(i))+i2(i)
_IFN1(iv)      integ(ic4+1)=iky(i3(i))+i4(i)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(i))+i2(i)
_IF1(iv)      intkl(icount)=iky(i3(i))+i4(i)
      gout(icount) = valp1(i)
      gout(num2ejk+icount) = valk1(i)
      icount = icount+1
      if (icount .le. nintmx) go to 560
c
      call blocki
      if(omaxb)go to 782
 560   continue
c
      do 570 i=1,loop2
      if ( dabs(valp2(i)) .lt. cutoff. and.  dabs(valk2(i)).lt.cutoff)
     *go to 570
      iii=i2out(i)
      if (i1(iii) .ne. i2(iii) .or. i3(iii) .ne. i4(iii)) go to 5701
       valp2(i)=valp2(i)*pt5
       valk2(i)=valk2(i)*pt5
 5701 continue
_IFN1(iv)      integ(ic4  )=iky(i1(iii))+i3(iii)
_IFN1(iv)      integ(ic4+1)=iky(i2(iii))+i4(iii)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(iii))+i3(iii)
_IF1(iv)      intkl(icount)=iky(i2(iii))+i4(iii)
      gout(icount) = valp2(i)
      gout(num2ejk+icount) = valk2(i)
      icount = icount+1
      if (icount .le. nintmx) go to 570
c
      call blocki
      if(omaxb)go to 782
 570   continue
c
      do 590 i=1,loop3
      if ( dabs(valp3(i)) .lt. cutoff. and . dabs(valk3(i)).lt.cutoff)
     * go to 590
      iii=i3out(i)
      if (i1(iii) .ne. i2(iii) .or. i4(iii) .ne. i3(iii))
     * go to 5901
       valp3(i)=valp3(i)*pt5
       valk3(i)=valk3(i)*pt5
 5901 continue
_IFN1(iv)      integ(ic4  )=iky(i1(iii))+i4(iii)
_IFN1(iv)      integ(ic4+1)=iky(i2(iii))+i3(iii)
_IFN1(iv)      ic4 = ic4 + 2
_IF1(iv)      intij(icount)=iky(i1(iii))+i4(iii)
_IF1(iv)      intkl(icount)=iky(i2(iii))+i3(iii)
      gout(icount) = valp3(i)
      gout(num2ejk+icount) = valk3(i)
      icount = icount+1
      if (icount .le. nintmx) go to 590
c
      call blocki
      if(omaxb)go to 782
 590   continue
  782 return
      end
      subroutine zsortp(gbuf,klbuf,iptbuf,ijbas)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
INCLUDE(common/sortp)
      dimension gbuf(lengbf,*),klbuf(lengbf,*),iptbuf(*)
      dimension gx(340),ijx(340),klx(340)
c ***
INCLUDE(common/shlt)
      common/craypk/intij(340),intkl(340),ipad(680)
INCLUDE(common/restar)
      common/blkin/goutx(510),nword
c ***
c ***
c *** copy goutx, intij, intkl into local work space,
c *** sort contents into buffers writing out as usual
c ***
      mword = icount - 1
      icount = 1
_IFN1(iv)      ic4 = 1
      if(mword.le.0) return
      call dcopy(mword,goutx,1,gx,1)
_IF1(c)      call dcopy(mword,intij(1),2,ijx,1)
_IF1(c)      call dcopy(mword,intij(2),2,klx,1)
_IFN1(civ)      call icopy(mword,intij(1),2,ijx,1)
_IFN1(civ)      call icopy(mword,intij(2),2,klx,1)
_IF1(iv)      call icopy(mword,intij,1,ijx,1)
_IF1(iv)      call icopy(mword,intkl,1,klx,1)
      do 10 i=1,mword
          ibuf = ijx(i) - ijbas
          if(ibuf.le.0 .or. ibuf.gt.ngbf)
     &      call caserr('ibuf value is invalid.')
          iptbuf(ibuf) = iptbuf(ibuf) + 1
          gbuf (iptbuf(ibuf),ibuf) = gx(i)
          klbuf(iptbuf(ibuf),ibuf) = klx(i)
          if(iptbuf(ibuf).eq.lengbf) call clrgbf(ibuf,gbuf,klbuf,iptbuf)
10    continue
      return
      end
      subroutine clrgbf(ibuf,gbuf,klbuf,iptbuf)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c ***
INCLUDE(common/sortp)
      dimension gbuf(lengbf,*),klbuf(lengbf,*),iptbuf(*)
c ***
_IFN1(iv)      common/craypk/intij(680),ipad(680)
_IF1(iv)      common/craypk/intij(340),intkl(340)
INCLUDE(common/restar)
      common/blkin/goutx(340)
INCLUDE(common/shlt)
c ***
c *** empty out buffer ibuf with usual blocki call sequence
c ***
c *** will not work for lengbf other than 340 as /blkin/ and
c *** /crypak/ will be destroyed by zsortp.
c ***
      mword = iptbuf(ibuf)
      if (mword.le.0) return
      ipt = 1
      ijval = ijbase + ibuf
10    nl = nintmx - icount + 1
      nleft = min(nl,mword)
c *** put smallest value at front to minimise precision loss
_IFN1(x)      i = idamin(nleft,gbuf(ipt,ibuf),1)
_IFN1(x)      i=i+ipt-1
_IF1(x)c **** NOTE ... idamin appeared to give spurious results
_IF1(x)c ****          on the convex C1  *** BEWARE ***
_IF1(x)c     lmin = idamin(nleft,gbuf(ipt,ibuf),1)
_IF1(x)      smin = dabs(gbuf(ipt,ibuf))
_IF1(x)      lmin=ipt
_IF1(x)      do 30 loop = lmin+1,lmin+nleft-1
_IF1(x)         if(dabs(gbuf(loop,ibuf)).ge.smin) go to 30
_IF1(x)         lmin = loop
_IF1(x)         smin = dabs(gbuf(loop,ibuf))
_IF1(x)   30 continue
_IF1(x)c ****          on the convex C1  *** BEWARE ***
_IF1(x)      i=lmin
      gx = gbuf(ipt,ibuf)
      gbuf(ipt,ibuf) = gbuf(i,ibuf)
      gbuf(i,ibuf) = gx
      m = klbuf(ipt,ibuf)
      klbuf(ipt,ibuf) = klbuf(i,ibuf)
      klbuf(i,ibuf) = m
_IFN1(iv)      iqq = 2*icount-1
_IFN1(iv)      do 11 i=1,nleft
_IFN1(iv)          intij(iqq)=ijval
_IFN1(iv)          intij(iqq+1) = klbuf(i+(ipt-1),ibuf)
_IFN1(iv)11        iqq = iqq + 2
_IF1(iv)      iqq=icount
_IF1(iv)      do 11 i=1,nleft
_IF1(iv)      intij(iqq)=ijval
_IF1(iv)      intkl(iqq)=klbuf(i+ipt-1,ibuf)
_IF1(iv) 11   iqq=iqq+1
      call dcopy(nleft,gbuf(ipt,ibuf),1,goutx(icount),1)
      call rinsrt(goutx(icount),nleft)
      icount = icount + nleft
_IFN1(iv)      ic4 = 2*icount-1
      if (icount .gt. nintmx) call blocki
      mword = mword - nleft
      ipt = ipt + nleft
      if(mword .gt. 0) goto 10
      iptbuf(ibuf)=0
      return
      end
      subroutine rinsrt(r,m)
      implicit REAL  (a-h,o-z)
c ***
c *** insert last len=9 bits of m into r at position ipos=55.
c ***
_IF(cray)
      r=or(and(compl(511),r),and(511,m))
_ELSEIF(convex)
      integer *8 r,cfiv11,fiv11,m8temp
      dimension m4temp(2)
      equivalence (m8temp,m4temp(1))
      data cfiv11,fiv11/-512,511/
      m4temp(1)=0
      m4temp(2)=m
      r=ior(iand(cfiv11,r),iand(fiv11,m8temp))
_ELSEIF(ksr)
      integer r,cfiv11,fiv11,m8temp
      integer *4  m4temp(2)
      equivalence (m8temp,m4temp(1))
      data cfiv11,fiv11/-512,511/
      m4temp(1)=0
      m4temp(2)=m
      r=or(and(cfiv11,r),and(fiv11,m8temp))
_ELSEIF(i8) 
      integer r,m,mtemp 
      mtemp = IAND64(511,m) 
      r=IAND64(-512,r) 
      r=IOR64(r,mtemp) 
_ELSEIF(titan)
      integer*4 r(2),m
      call mvbits(m,0,9,r(2),0)
_ELSEIF(alliant,dec)
      integer*4 r(2),m
      call mvbits(m,0,9,r(1),0)
_ELSEIF(ibm)
      dimension m8temp(2)
      equivalence (temp,m8temp(1))
      data cfiv11,fiv11/zfffffffffffffe00,
     *                  z1ff/
      data m8temp/0,0/
      m8temp(2)=m
      r=or(and(cfiv11,r),and(fiv11,temp))
_ELSEIF(littleendian)
      dimension m8temp(2)
      equivalence (temp,m8temp(1))
      temp = r
      m8temp(1) = IOR32(IAND32(m8temp(1),not(511)),m)
      r = temp
_ELSE
      integer r
      dimension r(2)
      data mtemp/0/
      r(2)=IAND32(-512,r(2))
      mtemp=IAND32(511,m)
      r(2)=IOR32(r(2),mtemp)
_ENDIF
      return
      end
      subroutine dbutci(ista,jsta,ksta,lsta)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/restar/nprint,itol,icut,normf,normp,nopk,
     *                     nrest,nrec,intloc,ist,jst,kst,lst
     *   ,nspfil(434), n5file,n5tape(20),n5blk(20),n5last(20)
INCLUDE(common/dmisc)
INCLUDE(common/iofile)
INCLUDE(common/filel)
INCLUDE(common/grad2)
INCLUDE(common/dm)
_IF(parallel)
INCLUDE(common/infoa)
_ENDIF
      common/blkin/qqqq(510),mword
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
      write (iwr,9008)
c
c     ----- set starting parameters -----
c
      outd = nprint .eq. -4
      icutd=icut-1
      icutd=max(icutd,8)
      itold=itol
      itold=max(itold,16)
      dcut=done/(ten**icutd)
      tol1=e*(itold+2)
      tol2=e*itold
      tol3=done/ten**(itold+2)
      tol4=done/ten**itold
      onorm=normf.ne.1.or.normp.ne.1
c
c     ----- read in 1e-gradient -----
c
      call rdgrd(de,nrest,ista,jsta,ksta,lsta)
_IF(parallel)
      factor = 1.0d0/(ipg_nnodes()*1.0d0)
      call dscal(nat*3,factor,de,1)
_ENDIF
c
c     ----- position 2-particle density matrix file
c
      lfile=n5file
      do 1 i=1,lfile
      lotape(i)=n5tape(i)
      liblk  (i)=n5blk (i)
 1    llblk  (i)=liblk(i)-n5last(i)
c
      ifiled=1
      mword=0
      iwordd=0
      idmmc=lotape(1)
      iblkmm =liblk(1)
      lblmm  =llblk(1)
      call search(iblkmm,idmmc)
c
      return
 9008 format(/1x,22('-')/1x,'gradient of the energy'/1x,22('-')/)
      end
      subroutine dabci(ii,jj,kk,ll,q4,oeof,abdens)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/dm)
INCLUDE(common/atmblk)
      common/blkin/gin(510),mword
_IFN1(iv)      common/craypk/ijkl(4,340)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
INCLUDE(common/filel)
INCLUDE(common/nshel)
      common/tgrad/dgout(9)
INCLUDE(common/incrd)
      dimension abdens(*)
c
      data eight,four,pt5/8.0d0,4.0d0,0.5d0/
      call vclr(abdens,1,lendd)
      if(oeof)return
      iword=iwordd
      oijdif=ii.ne.jj
      okldif=kk.ne.ll
      oikdif=ii.ne.kk.or.jj.ne.ll
      mini=kloc(ii)
      minj=kloc(jj)
      mink=kloc(kk)
      minl=kloc(ll)
      maxi=kloc(ii+1)-1
      maxj=kloc(jj+1)-1
      maxk=kloc(kk+1)-1
      maxl=kloc(ll+1)-1
1     if(iword.ne.mword)goto10
      if(lblmm.ne.0)go to 5
      ifiled=ifiled+1
      if(ifiled.gt.lfile)go to 888
      idmmc=lotape(ifiled)
      call search(liblk(ifiled),idmmc)
      lblmm=llblk(ifiled)
 5    callfind(idmmc)
      callget(gin,m)
      lblmm=lblmm+1
      if(m)2,888,2
2     iword=0
_IFN1(iv)      call unpack(gin(num2e+1),lab816,ijkl,numlab)
_IF1(iv)      call upak8v(gin(num2e+1),i205)
10    iword=iword+1
_IF(ibm,vax)
      i = i205(iword)
      j = j205(iword)
      k = k205(iword)
      l = l205(iword)
_ELSEIF(littleendian)
      j = ijkl(1,iword)
      i = ijkl(2,iword)
      l = ijkl(3,iword)
      k = ijkl(4,iword)
_ELSE
      i = ijkl(1,iword)
      j = ijkl(2,iword)
      k = ijkl(3,iword)
      l = ijkl(4,iword)
_ENDIF
      if(i.lt.mini)goto1
      if(i.gt.maxi)goto999
      if(j.lt.minj)goto1
      if(j.gt.maxj)goto999
      if(k.lt.mink)goto1
      if(k.gt.maxk)goto999
      if(l.lt.minl)goto1
      if(l.gt.maxl)goto999
      d4=eight
      if(i.eq.j)d4=four
      if(k.eq.l)d4=d4*pt5
      nn=(i-mini)*inc5+(j-minj)*inc4+(k-mink)*inc3+l-minl+inc2
      buff=gin(iword)*q4*d4
      abdens(nn)=buff
      if(oijdif)goto20
      nn=(j-minj)*inc5+(i-mini)*inc4+(k-mink)*inc3+l-minl+inc2
      abdens(nn)=buff
      if(okldif)goto20
      nn=(j-minj)*inc5+(i-mini)*inc4+(l-minl)*inc3+k-mink+inc2
      abdens(nn)=buff
20    if(okldif)goto30
      nn=(i-mini)*inc5+(j-minj)*inc4+(l-minl)*inc3+k-mink+inc2
      abdens(nn)=buff
30    if(oikdif)goto1
      nn=(k-mink)*inc5+(l-minl)*inc4+(i-mini)*inc3+j-minj+inc2
      abdens(nn)=buff
      if(oijdif)goto40
      nn=(k-mink)*inc5+(l-minl)*inc4+(j-minj)*inc3+i-mini+inc2
      abdens(nn)=buff
      if(okldif)goto40
      nn=(l-minl)*inc5+(k-mink)*inc4+(j-minj)*inc3+i-mini+inc2
      abdens(nn)=buff
40    if(okldif)goto1
      nn=(l-minl)*inc5+(k-mink)*inc4+(i-mini)*inc3+j-minj+inc2
      abdens(nn)=buff
      goto1
999   iwordd=iword-1
      return
888   oeof=.true.
      return
      end
      function chrint(i)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c converts an integer to a left justified character string
c without leading zeros, and with minus sign (if negative)
c   works for up to 15 decimal digit integers
      character*16 c,chrint
      character*1 digit(0:9)
      data c/' '/
      data digit/'0','1','2','3','4','5','6','7','8','9'/
      num=iabs(i)
      l=16
      do 1 j=16,2,-1
      new=num/10
      n=num-10*new
      if(n.gt.0) l=j
      c(j:j)=digit(n)
1     num=new
      if(i.lt.0) then
          l=l-1
          c(l:l)='-'
      endif
      chrint=c(l:16)
      return
      end
      function lstchr(a)
      character*(*) a
c routine to return index of last non-blank character in character var
      n=len(a)
      do 10 i=n,1,-1
      if( a(i:i).ne.' ' ) go to 20
10    continue
      i=0
20    lstchr=i
      return
      end
      subroutine wrtc(z,nword,iblk,num3)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension z(*)
INCLUDE(common/junkcc)
_IF(ga)
INCLUDE(common/sizes)
INCLUDE(common/disc)
      If( iwhere( num3 ) .EQ. 6 .Or. iwhere( num3 ) .EQ. 7 ) Then 
         Call wrtc_ga( z, nword, iblk, num3 )
      Else
_ENDIF
      call search(iblk,num3)
      l=511
      j=1
      k=nword
 1    k=k-511
      if(k)10495,2,2
10495 l=k+511
c2    call cmoved(z(j),dchar,l)
 2    do 3 loop=1,l
 3    read(z(j+loop-1),'(a8)') dchar(loop)
      call put(dchar,l,num3)
      j=j+511
      if(k)10500,10500,1
10500 return
_IF(ga)
      End If
_ENDIF
      end
      subroutine rdchr(z,nword,iblk,num3)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/junkcc)
      dimension z(*)
_IF(ga)
INCLUDE(common/sizes)
INCLUDE(common/disc)
      If( iwhere( num3 ) .EQ. 6 .Or. iwhere( num3 ) .EQ. 7 ) Then 
         Call rdchr_ga( z, nword, iblk, num3 )
      Else
_ENDIF
      call search(iblk,num3)
      j=1
10491 call find(num3)
      call get(dchar,l)
c     call dmovec(dchar,z(j),l)
      do 3 loop=1,l
 3    write(z(j+loop-1),'(a8)') dchar(loop)
      j=j+l
      if(j.le.nword)go to 10491
      if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present ', nword, j-1
       call caserr('invalid number of words in rdchr')
      endif
      return
_IF(ga)
      End If
_ENDIF
      end
      subroutine wrtcs(z,nword,num3)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension z(*)
INCLUDE(common/junkcc)
_IF(ga)
INCLUDE(common/sizes)
INCLUDE(common/disc)
      If( iwhere( num3 ) .EQ. 6 .Or. iwhere( num3 ) .EQ. 7 ) Then 
         Call wrtcs_ga( z, nword, num3 )
      Else
_ENDIF
      l=511
      j=1
      k=nword
 1    k=k-511
      if(k)10495,2,2
10495 l=k+511
c2    call cmoved(z(j),dchar,l)
 2    do 3 loop=1,l
 3    read(z(j+loop-1),'(a8)') dchar(loop)
      call put(dchar,l,num3)
      j=j+511
      if(k)10500,10500,1
10500 return
_IF(ga)
      End If
_ENDIF
      end
      subroutine rdchrs(z,nword,num3)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/junkcc)
      dimension z(*)
_IF(ga)
INCLUDE(common/sizes)
INCLUDE(common/disc)
      If( iwhere( num3 ) .EQ. 6 .Or. iwhere( num3 ) .EQ. 7 ) Then 
         Call rdchrs_ga( z, nword, num3 )
      Else
_ENDIF
      j=1
10491 call find(num3)
      call get(dchar,l)
c     call dmovec(dchar,z(j),l)
      do 3 loop=1,l
 3    write(z(j+loop-1),'(a8)') dchar(loop)
      j=j+l
      if(j.le.nword)go to 10491
      if(j-1.ne.nword) then
       write(6,*)
     + 'invalid no of words: requested,present ', nword, j-1
       call caserr('invalid number of words in rdchrs')
      endif
      return
_IF(ga)
      End If
_ENDIF
      end
_IF1(iv)      subroutine setstr(nf,value,array)
_IF1(iv)      implicit REAL  (a-h,o-z),  integer  (i-n)
_IF1(iv)      dimension array(*)
_IF1(iv)      do 1 loop=1,nf
_IF1(iv)    1 array(loop)=value
_IF1(iv)      return
_IF1(iv)      end
      subroutine setstc(nf,ztext,zlab)
      character *(*) ztext,zlab
      dimension zlab(*)
      if(nf.eq.0)go to 2
      do 1 i=1,nf
 1    zlab(i)=ztext
 2    return
      end

_IFN(cray,ibm,convex,ksr)
      function eq(i,j)
c alliant compiler only compares lowest 32 bits of i8 quantities
c so have to use eq if top 32 bits may be significant
      logical eq
      integer i(2),j(2)
      eq = (i(1).eq.j(1)) .and. (i(2).eq.j(2))
      end
_ENDIF

      function locatz(label,nf,itext)
      implicit REAL  (a-h,p-z),integer    (i-n)
_IF(convex)
      integer *8 label,itext
      dimension label(*)
_ELSEIF(cray,ksr)
      dimension label(*)
_ELSEIF(alliant,ibm,vax)
      REAL  label,itext
      dimension label(*)
      logical eq
_ELSE
      integer label(2,*),itext(2)
      logical eq
_ENDIF
      do 1 i=1,nf
_IF(cray,convex,ksr)
      if(label(i).eq.itext)go to 2
_ELSEIF(ibm,vax,alliant)
      if(eq(label(i),itext))go to 2
_ELSE
      if(eq(label(1,i),itext))go to 2
_ENDIF
 1    continue
      locatz=0
      return
 2    locatz=i
      return
      end
      subroutine hstarg(trans,coul,exch,dens,nopk,
     + nclf,nscmf,iscmf,l1,l2,nprint)
c
c     ----- -hstarg- forms multiple fock matrices 
c                    doing only one integral pass -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/ij205(204),kl205(204)
INCLUDE(common/filel)
      common/blkin/gout(510),nint
INCLUDE(common/mapper)
      dimension coul(*),exch(*),dens(*)
INCLUDE(common/infoa)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 imemo
_ENDIF
INCLUDE(common/atmblk)
      common/vinteg/q(1)
      dimension trans(l1,*)
      data pt5,two/0.5d0,2.0d0/
c
      nshell=nclf+nscmf
c     load up the density matrices
      ioff = 1
      if (nclf .gt. 0) then
c      closed shell part
       call dengvb(trans,dens,1,iky,l1)
       ioff=ioff+l2
      endif
c     open shell part
      if (nscmf .gt. 0) then
       ilow = iscmf
       ihi = ilow + nscmf - 1
       do 200 i = ilow,ihi
       call dengvb(trans,dens(ioff),i,iky,l1)
 200   ioff = ioff + l2
      endif
c
c     ----- divide the diagonal density elements by 2
c
      ioff=0
      do 300 i = 1,nshell
      if (nprint.eq.5) then
       write(iwr,*)' density matrix, shell = ',i
       call prtri(dens(ioff+1),l1)
      endif
      do 1 j = 1,l1
      nij = ikyp(j) + ioff
   1  dens(nij)=dens(nij)*pt5
 300  ioff=ioff+l2
c
      if (nopk .ne. 1) go to 380
c
c      integrals in non-supermatrix form
c
      ndzero=nshell*l2
_IF1(civ)      call szero(coul,ndzero)
_IFN1(civ)      call vclr(coul,1,ndzero)
      call vclr(exch,1,ndzero)
c
      do 2000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 2000
c
       if (o255i) then
       call sgmatn(l2,nshell,coul,exch,dens,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       else
       call sgmatn_255(l2,nshell,coul,exch,dens,
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
1050   imemo = imemo + 512
      else
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 2004  jblock=jblock+1
       call get(gout,nw)
       if(nw)2001,2000,2001
 2001  if(jblock.ne.0)call find(main)
c
_IF(ibm,vax)
       call sgmat(l2,nshell,iky,coul,exch,dens)
_ELSE
       if(o255i) then
          call sgmat(l2,nshell,coul,exch,dens)
       else
          call sgmat_255(l2,nshell,coul,exch,dens)
       endif
_ENDIF
c
       if(jblock)2004,2000,2004
      endif
 2000 continue
c
      if(nclf.gt.0) then
       if (nprint.eq.5) then
        write(iwr,*)' closed shell J matrix'
        call prtri(coul,l1)
        write(iwr,*)' closed shell K matrix'
        call prtri(exch,l1)
       endif
c
       call dscal(l2,two,coul,1)
       call vsub(coul,1,exch,1,coul,1,l2)
c
      endif
_IF(parallel)
      go to 750
_ELSE
      return
_ENDIF
c
c     ----- integrals in supermatrix form -----
c
  380 continue
      labss = num2ejk + num2ejk + 1
      ndzero=nshell*l2
      call vclr(coul,1,ndzero)
      call vclr(exch,1,ndzero)
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
c     memory resident integrals
c
       imemo = iqqoff(main) + 1
       do 2020 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 1000
c
       call jksupm(coul,coul(l2+1),exch(l2+1),dens,dens(l2+1),nscmf,
     + l2, q(imemo), q(imemo+num2ejk), q(imemo+num2ejk+num2ejk), 
     +     q(imemo+510) )
2020   imemo = imemo + 512
      else
c
c     file resident integrals
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gout,nw)
       if(nw)1001,1000,1001
1001   if(jblock.ne.0)call find(main)
c
c ----------build the fock matrices--------------
c
_IFN1(iv)       iword=1
       if (nclf.ne.0.and.nscmf.eq.0) then
c
c---------closed shells only------------
c
c     unpack a buffer of labels
c
_IFN1(iv)        call unpack(gout(labss),lab1632,integ,numlabjk)
_IF1(iv)        call upak6v(gout(labss),ij205)
        do 480 m = 1,nint
_IFN1(iv)        ij = integ(iword  )
_IFN1(iv)        kl = integ(iword+1)
_IF1(iv)        ij = ij205(m)
_IF1(iv)        kl = kl205(m)
        goutpm = gout(m)
        coul(ij) = coul(ij) + goutpm*dens(kl)
        coul(kl) = coul(kl) + goutpm*dens(ij)
_IFN1(iv)        iword=iword+2
 480    continue
c      
       else if (nclf.eq.0.and.nscmf.ne.0) then
c
c---------open shells only----------------
c
c     unpack a buffer of labels
c
_IFN1(iv)       call unpack(gout(labss),lab1632,integ,numlabjk)
_IF1(iv)       call upak6v(gout(labss),ij205)
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
        do 500 m = 1,nint
_IFN1(iv)        ij = integ(iword  )
_IFN1(iv)        kl = integ(iword+1)
_IF1(iv)        ij = ij205(m)
_IF1(iv)        kl = kl205(m)
        goutpm = gout(m)
        goutkm = gout(num2ejk+m)
        do 501 i = 1,nscmf
        coul(ij) = coul(ij) + goutpm*dens(kl)
        exch(ij) = exch(ij) + goutkm*dens(kl)
        coul(kl) = coul(kl) + goutpm*dens(ij)
        exch(kl) = exch(kl) + goutkm*dens(ij)
        ij = ij + l2
 501    kl = kl + l2
_IFN1(iv)        iword=iword+2
 500    continue
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
c
       else
c
c---------both closed + open shells------------
c
       call jksupe(coul,coul(l2+1),exch(l2+1),dens,dens(l2+1),nscmf,l2)
      endif
c
c-----------check next buffer-------------------
c
       if(jblock) 1004,1000,1004
      endif
 1000 continue
c
      if(nscmf.gt.0) then
       ioff = 1
       if(nclf.ne.0) ioff = ioff + l2
       do 720 ifo = 1,nscmf
       call vadd(coul(ioff),1,exch(ioff),1,coul(ioff),1,l2)
       call dscal(l2,pt5,coul(ioff),1)
 720   ioff=ioff+l2
      endif
_IF(parallel)
c***   ***node-MPP***
c...  gop results together
 750  ioff = 1
      if (nclf.ne.0) then
         call pg_dgop(300,coul(ioff),l2,'+')
         ioff = ioff + l2
      end if
      do 741 ifo=1,nscmf
         call pg_dgop(400+ioff,coul(ioff),l2,'+')
         call pg_dgop(500+ioff,exch(ioff),l2,'+')
741   ioff = ioff + l2
c***   ***node-MPP***
_ENDIF
      return
      end
_IF1()c ********
_IF1()c ******** hstar - closed shell skeleton fock builder
_IF1()c ********
_IF(ibm,vax)
      subroutine hstar(d,f,exch,nopk)
c
c     ----- -hstar- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order-
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
c     ----- ***********************************  ------
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/filel)
INCLUDE(common/iofile)
      dimension d(*),f(*),exch(*)
INCLUDE(common/infoa)
      common/blkin/gin(510),nint
INCLUDE(common/mapper)
INCLUDE(common/sortp)
INCLUDE(common/atmblk)
c
      data pt5/0.5d0/
c     data two/2.0d0/
c
      call szero(f,nx)
c
      if(nopk.eq.1) go to 2
      do 1 m=1,num
      nij=ikyp(m)
   1  d(nij)=d(nij)*pt5
      go to 280
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
    2 continue
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      call search(liblk(ifile),main)
      call find(main)
      jblock=llblk(ifile)
 1004 jblock=jblock+1
      call get(gin,nw)
      if(nw.eq.0)go to 1005
      if(jblock.ne.0)call find(main)
c
      call gmake(f,d,iky)
c
      if(jblock)1004,1000,1004
 1005 llblk(ifile)=liblk(ifile)-iposun(main)+1
 1000 continue
c
      go to 4
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c
 280  continue
c
      do 2000 ifile=1,lfile
      main=lotape(ifile)
      call search(liblk(ifile),main)
      call find(main)
      jblock=llblk(ifile)
 2004 jblock=jblock+1
      call get(gin,nw)
      if(nw)2001,2005,2001
 2001 if(jblock.ne.0)call find(main)
c
      if(nopk.eq.-1)call gsupa(f,d)
c
      if(nopk.ne.-1) then
        if(osortp) then
          call gsupp(f,d)
        else
          call gsup(f,d)
        endif
      endif
c
      if(jblock)2004,2000,2004
 2005 llblk(ifile)=liblk(ifile)-iposun(main)+1
 2000 continue
c
 360  continue
      do 3   m = 1,num
      nij = ikyp(m)
   3  d(nij) = d(nij)+d(nij)
    4 continue
       call dscal(nx,pt5,f,1)
      return
      end
_ELSEIF(cray,convex,titan)
      subroutine hstar(d,f,exch,nopk,mss)
c
c     ----- -hstar- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order-
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
c     ----- ***********************************  ------
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/filel)
INCLUDE(common/iofile)
      dimension d(*),f(*),exch(*)
_IF(parallel)
c***   ***node-MPP***
      common/scftim/tdiag(3),tfock,tdiis(6)
c***   ***node-MPP***
_ENDIF
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
INCLUDE(common/infoa)
_IF(1s,oldx)
      common/blkin/gin(510),nint
_ENDIF
INCLUDE(common/atmblk)
_IF(unicos)
      common/blkin/gin(340),ijkl(851)
      common/junk2/ii(1020),jj(1020),kk(1020),ll(1020)
      common/bufbb/indf(6,1020),ibias(61),pg(6,1020)
_ENDIF
_IF(convex,titan)
      common/blkin/gin(num2e),gijkl(851)
      common/junk2/ii(1020),jj(1020),kk(1020),ll(1020)
      common/bufbb/indf(6,1020),ibias(62),pg(6,1020)
_ENDIF
INCLUDE(common/mapper)
INCLUDE(common/sortp)
_IF(convex,titan)
      dimension nword(2)
      equivalence (dum,nword(1),m)
_ENDIF
c
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 imemo
_ENDIF
      common/vinteg/q(1)
_IF(cray)
      equivalence (mww,dumm)
_ELSE
      dimension mww(2)
      equivalence (mww(1),dumm)
_ENDIF
c
      data pt5/0.5d0/
c     data two/2.0d0/
_IF(1s,oldx)
      call szero(f,nx)
_ENDIF
_IF(parallel)
      dumtim=dclock()
      call vclr(f,1,nx)
_ENDIF
c
      do 1 m=1,num
      nij=ikyp(m)
   1  d(nij)=d(nij)*pt5
      if (nopk.ne.1) go to 280
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
c    strategy a function of whether modified dft fock matrices are
c    required - if so, abandon the vectorised/cal strategy and
c    fall back to single-block fortran processing
c
_IF(ccpdft)
      odft = CD_active()
      onodft = .not. odft
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ELSE
      odft = .false.
_ENDIF
      ovect = o255i .and. .not. odft
      lenb=nx
      if (ovect) then
_IF(unicos)
         len2=lenb+lenb
         ms=mss
         mw=len2
         msmall=ms*6
         call setsto(4,0,ibias)
         call setsto(2,lenb,ibias(5))
         i=6
4000     j=msmall-i
         if(j.lt.1)goto 3000
         nw=min(i,j)
         call ujau(nw,mw,ibias(i+1),ibias)
         mw=mw+mw
         i=i+nw
         goto 4000
3000     call szero(f,len2*ms)
_ELSE
         len2=lenb+lenb
         ms=mss
         call vclr(f,1,len2*ms)
         call setsto(4080,0,ii)
_ENDIF
      else
_IF(unicos)
         call szero(f,lenb)
_ELSE
         call vclr(f,1,lenb)
_ENDIF
      endif
      mw = 0
c
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if (omem(main)) then
c     memory processing of 2e-integral file
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
       dumm = q(imemo+511)
_IF(cray)
       call fmove(q(imemo),gin(mw+1),511)
       mmm = shiftr(mww,32)
       if(mmm.eq.0) go to 1000
_ELSE
_IF1(xt)       call dcopy(511,q(imemo),1,gin(mw+1),1)
_IF1(x)       if(mww(1).eq.0) go to 1000
_IF1(t)       if(mww(2).eq.0) go to 1000
_ENDIF
c
_IF(unicos)
       m=ijkl(mw+171)
       call upak8z(m,ijkl(mw+1),ii(mw+1),jj(mw+1),
     * kk(mw+1),ll(mw+1))
       mw=mw+m
       if(mw.ge.681)then
        call jkuno(mw,d,gin,ii,jj,kk,ll,indf,pg)
        call jkadd(msmall,mw*6,f,indf,ibias,pg)
        mw=0
       endif
_ELSE
       dum=gijkl(mw+171)
       call upak8z(m,gijkl(mw+1),ii(mw+1),jj(mw+1),
     *  kk(mw+1),ll(mw+1))
       mw=mw+m
       if(mw.ge.681)then
        call jandk2(lenb,1,ms,f,d,ii,jj,kk,ll,gin,pg,indf,mw)
        mw=0
       endif
_ENDIF
c
 1050  imemo = imemo + 512
      else
c      conventional processing of 2e-integral file
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gin(mw+1),nw)
       if(nw.eq.0)go to 1005
       if(jblock.ne.0)call find(main)
c
       if (ovect) then
_IF(unicos)
          m=ijkl(mw+171)
          call upak8z(m,ijkl(mw+1),ii(mw+1),jj(mw+1),
     *    kk(mw+1),ll(mw+1))
          mw=mw+m
          if(mw.ge.681)then
           call jkuno(mw,d,gin,ii,jj,kk,ll,indf,pg)
           call jkadd(msmall,mw*6,f,indf,ibias,pg)
           mw=0
          endif
_ELSE
          dum=gijkl(mw+171)
          call upak8z(m,gijkl(mw+1),ii(mw+1),jj(mw+1),
     *     kk(mw+1),ll(mw+1))
          mw=mw+m
          if(mw.ge.681)then
           call jandk2(lenb,1,ms,f,d,ii,jj,kk,ll,gin,pg,indf,mw)
           mw=0
          endif
_ENDIF
       else if (.not.o255i) then
            call sgmata_255(f,d)
_IF(ccpdft)
       else
            call sgmata_dft(f,d,facex,ocoul,oexch)
_ENDIF
       endif
c
       if(jblock)1004,1000,1004
 1005  llblk(ifile)=liblk(ifile)-iposun(main)+1
      endif
 1000 continue
c
      if (ovect) then
_IF(unicos)
         if(mw.eq.0)goto 6000
         call jkuno(mw,d,gin,ii,jj,kk,ll,indf,pg)
         call jkadd(msmall,mw*6,f,indf,ibias,pg)
6000     if(ms.eq.1)go to 6010
         j=ms/2
         ms=ms-j
         nw=j*len2
         i=ms*len2
         call vadd(f,1,f(i+1),1,f,1,nw)
         goto 6000
6010     call gtrian(lenb,4.0d0,f,f,exch)
_ELSE
         if(mw.eq.0)goto 6000
         call jandk2(lenb,1,ms,f,d,ii,jj,kk,ll,gin,pg,indf,mw)
6000     if(ms.eq.1)go to 6010
         j=ms/2
         ms=ms-j
         nw=j*len2
         i=ms*len2
         call vadd(f,1,f(i+1),1,f,1,nw)
         goto 6000
6010     do 6030 loop=1,lenb
6030     f(loop)=4.0d0*exch(loop)-f(loop)
_ENDIF
      endif
c
      go to 360
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c
 280  continue
_IF(unicos)
      call szero(f,nx)
_ELSE
      call vclr(f,1,nx)
_ENDIF
c
c
       do 2000 ifile=1,lfile
       main=lotape(ifile)
       if(omem(main)) then
        imemo = iqqoff(main) + 1
        if(nopk.ne.-1.and.osortp) then
         do 2010 jblock=1,mxblock
         dumm = q(imemo+511)
_IF(cray)
         mmm = shiftr(mww,32)
         if(mmm.eq.0) go to 2000
_ELSE
_IFN1(x)         if(mww(2).eq.0) go to 2000
_IF1(x)         if(mww(1).eq.0) go to 2000
_ENDIF
         call gsuppm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2010     imemo = imemo + 512
c
        else if(nopk.ne.-1.and..not.osortp) then
c
         do 2020 jblock=1,mxblock
         dumm = q(imemo+511)
_IF(cray)
         mmm = shiftr(mww,32)
         if(mmm.eq.0) go to 2000
_ELSE
_IFN1(x)         if(mww(2).eq.0) go to 2000
_IF1(x)         if(mww(1).eq.0) go to 2000
_ENDIF
         call gsupm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2020     imemo = imemo + 512
c
        else 
c
         do 2030 jblock=1,mxblock
         dumm = q(imemo+511)
_IF(cray)
         mmm = shiftr(mww,32)
         if(mmm.eq.0) go to 2000
_ELSE
_IFN1(x)         if(mww(2).eq.0) go to 2000
_IF1(x)         if(mww(1).eq.0) go to 2000
_ENDIF
         call gsupam(f,d,q(imemo),q(imemo+num2ejk+num2ejk),
     +               q(imemo+510) )
2030     imemo = imemo + 512
        endif
c
       else
c
        call search(liblk(ifile),main)
        call find(main)
        jblock=llblk(ifile)
 2004   jblock=jblock+1
        call get(gin,nw)
        if(nw)2001,2005,2001
 2001   if(jblock.ne.0)call find(main)
c
        if(nopk.eq.-1)call gsupa(f,d)
c
        if(nopk.ne.-1) then
          if(osortp) then
            call gsupp(f,d)
          else
            if(o255i) then
               call gsup(f,d)
            else
               call gsup_255(f,d)
            endif
          endif
        endif
c
        if(jblock)2004,2000,2004
 2005   llblk(ifile)=liblk(ifile)-iposun(main)+1
       endif
 2000  continue
c
c ----
c
 360  continue
      do 3   m = 1,num
      nij = ikyp(m)
   3  d(nij) = d(nij)+d(nij)
      call dscal(nx,pt5,f,1)
_IF(parallel)
c***   ***node-MPP***
c...   add fock-matrix and send it around
       call pg_dgop(600,f,nx,'+')
c***   ***node-MPP***
c
      tfock=tfock+(dclock()-dumtim)
_ENDIF
      return
      end
_ELSE
      subroutine hstar(d,f,exch,nopk)
c
c     ----- -hstar- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order-
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
c     ----- ***********************************  ------
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/filel)
INCLUDE(common/iofile)
      dimension d(*),f(*),exch(*)
_IF(parallel)
c***   ***node-MPP***
      common/scftim/tdiag(3),tfock,tdiis(6)
c***   ***node-MPP***
_ENDIF
INCLUDE(common/infoa)
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
      common/blkin/gin(510),nint
INCLUDE(common/atmblk)
INCLUDE(common/mapper)
INCLUDE(common/sortp)
INCLUDE(common/maxlen)
_IF(64bitpointers)
      integer*8 imemo
_ENDIF
      common/vinteg/q(1)
c
      data pt5/0.5d0/
c
_IF(parallel)
      dumtim=dclock()
      call vclr(f,1,nx)
_ELSE
_IF1(Gapgbdrsk)      call vclr(f,1,nx)
_ENDIF
c
      do 1 m=1,num
      nij=ikyp(m)
   1  d(nij)=d(nij)*pt5
      if (nopk.ne.1) go to 280
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
c    strategy a function of whether modified dft fock matrices are
c    required - if so, abandon the vectorised/cal strategy and
c    fall back to single-block fortran processing
c
_IF(ccpdft)
      odft = CD_active()
      onodft = .not. odft
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ELSE
      odft = .false.
_ENDIF
      ovect = o255i .and. .not. odft
_IF1(Gapgbdrhsk)      call vclr(exch,1,nx)
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
       call upack2(q(imemo+511),ibll,nw)  
       if(nw.eq.0) go to 1000
c
       if (o255i) then
          call sgmtmm(f,d,q(imemo),q(imemo+num2e),q(imemo+510) )
       else
          call sgmtmm_255(f,d,q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
c
 1050  imemo = imemo + 512
      else
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gin,nw)
       if(nw.eq.0)go to 1005
       if(jblock.ne.0)call find(main)
c
       if (ovect) then
          call sgmata(f,d)
       else if(.not.o255i) then
          call sgmata_255(f,d)
_IF(ccpdft)
       else
          call sgmata_dft(f,d,facex,ocoul,oexch)
_ENDIF
       endif
c
       if(jblock)1004,1000,1004
 1005  llblk(ifile)=liblk(ifile)-iposun(main)+1
      endif
 1000 continue
c
      go to 360
c
c     ----- integrals are in supermatrix form (nopk=.false.) -----
c
 280  continue
c
       do 2000 ifile=1,lfile
       main=lotape(ifile)
       if(omem(main)) then
        imemo = iqqoff(main) + 1
        if(nopk.ne.-1.and.osortp) then
         do 2010 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsuppm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2010     imemo = imemo + 512
c
        else if(nopk.ne.-1.and..not.osortp) then
c
         do 2020 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsupm(f,d,q(imemo),q(imemo+num2ep),q(imemo+510) )
2020     imemo = imemo + 512
c
        else 
c
         do 2030 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if(nw.eq.0) go to 2000
c
         call gsupam(f,d,q(imemo),q(imemo+num2ejk+num2ejk),
     +               q(imemo+510) )
2030     imemo = imemo + 512
        endif
c
       else
c
        call search(liblk(ifile),main)
        call find(main)
        jblock=llblk(ifile)
 2004   jblock=jblock+1
        call get(gin,nw)
        if(nw)2001,2005,2001
 2001   if(jblock.ne.0)call find(main)
c
        if(nopk.eq.-1)call gsupa(f,d)
c
        if(nopk.ne.-1) then
          if(osortp) then
            call gsupp(f,d)
          else
            if(o255i) then
               call gsup(f,d)
            else
               call gsup_255(f,d)
            endif
          endif
        endif
c
        if(jblock)2004,2000,2004
 2005   llblk(ifile)=liblk(ifile)-iposun(main)+1
       endif
 2000  continue
c
c ----
c
 360  continue
      do 3   m = 1,num
      nij = ikyp(m)
   3  d(nij) = d(nij)+d(nij)
      call dscal(nx,pt5,f,1)
_IF(parallel)
c***   ***node-MPP***
c...   add fock-matrix and send it around
       call pg_dgop(600,f,nx,'+')
c***   ***node-MPP***
c
      tfock=tfock+(dclock()-dumtim)
_ENDIF
      return
      end
_ENDIF
      subroutine hstaru(da,db,fa,fb,nopk)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
_IFN1(iv)      common/craypk/integ(680),ipad(680)
_IF1(iv)      common/craypk/ij205(204),kl205(204)
INCLUDE(common/mapper)
      common/blkin/gout(510),nint
INCLUDE(common/filel)
INCLUDE(common/maxlen)
INCLUDE(common/atmblk)
_IF(64bitpointers)
      integer*8 imemo
_ENDIF
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
      common/vinteg/q(1)
      dimension da(*),fa(*),db(*),fb(*)
      equivalence (gg1,gg),(ikyj,jl),(ikyk,jk),(ikyi,il)
INCLUDE(common/infoa)
      data pt5/0.5d0/
c     data two/2.0d0/
      do 100 m = 1,nx
      duma=da(m)
      dumb=db(m)
      da(m)=duma+dumb
 100  db(m)=duma-dumb
      call vclr(fa,1,nx)
      call vclr(fb,1,nx)
_IF1(iv)      if(nopk.eq.1) go to 2
      do 1   m=1,num
      nij=ikyp(m)
      da(nij)=da(nij)*pt5
   1  db(nij)=db(nij)*pt5
_IFN1(iv)      if(nopk.ne.1)goto280
_IF1(iv)      go to 280
_IF1(iv)   2  continue
c
c     ----- integrals are not in supermatrix form (nopk=.true.)
c
_IF(ccpdft)
      odft = CD_active()
      onodft = .not. odft
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
      do 1000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
       imemo = iqqoff(main) + 1
       do 1050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 1000
c
       if (o255i) then
       call proc2m(fa,da,fb,db,
_IF(ccpdft)
     +      facex,ocoul,oexch,
_ENDIF
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       else
       call proc2m_255(fa,da,fb,db,
_IF(ccpdft)
     +      facex,ocoul,oexch,
_ENDIF
     +      q(imemo),q(imemo+num2e),q(imemo+510) )
       endif
1050   imemo = imemo + 512
      else
c
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 1004  jblock=jblock+1
       call get(gout,nw)
       if(nw)1001,1000,1001
1001   if(jblock.ne.0)call find(main)
c
_IF(ccpdft)
_IF(ibm)
       if (odft) then
        call proc2f(fa,da,fb,db,iky,facex,ocoul,oexch)
       else
        call proc2(fa,da,fb,db,iky)
       endif
_ELSEIF(vax)
       if (o255i) then
        call proc2f(fa,da,fb,db,iky,facex,ocoul,oexch)
       else
        call proc2_255(fa,da,fb,db,iky,facex,ocoul,oexch)
       endif
_ELSEIF(cray)
       if (o255i) then
         if (odft) then
          call proc2f(fa,da,fb,db,facex,ocoul,oexch)
         else
          call proc2(fa,da,fb,db)
         endif
       else
          call proc2_255(fa,da,fb,db,iky,facex,ocoul,oexch)
       endif
_ELSE
       if (o255i) then
          call proc2f(fa,da,fb,db,facex,ocoul,oexch)
       else
          call proc2_255 (fa,da,fb,db,facex,ocoul,oexch)
       endif
_ENDIF
_ELSE
_IF(ibm,vax)
       call proc2 (fa,da,fb,db,iky)
_ELSE
       if(o255i) then
          call proc2f(fa,da,fb,db)
       else
          call proc2_255 (fa,da,fb,db)
       endif
_ENDIF
_ENDIF
c
       if(jblock)1004,1000,1004
c
      endif
 1000 continue
      do 700 i=1,nx
      duma=fa(i)
      dumb=fb(i)
      fa(i)=duma-dumb
 700  fb(i)=duma+dumb
c
_IFN1(iv)      go to 5
_IF1(iv)      go to 4
c
c     ----- integrals are in a supermatrix form (nopk=.false.)
c
 280  continue
c
      do 2000 ifile=1,lfile
      main=lotape(ifile)
      if(omem(main)) then
c
       imemo = iqqoff(main) + 1
       do 2050 jblock=1,mxblock
c
         call upack2(q(imemo+511),ibll,nw)  
         if (nw.eq.0) go to 2000
c
       call gsupum(fa,fb,da,db,
     +     q(imemo),q(imemo+num2ejk),q(imemo+num2ejk+num2ejk),
     +     q(imemo+510) )
2050   imemo = imemo + 512
      else
       call search(liblk(ifile),main)
       call find(main)
       jblock=llblk(ifile)
 2004  jblock=jblock+1
       call get(gout,nw)
       if(nw)2001,2000,2001
2001   if(jblock.ne.0)call find(main)
_IFN1(iv)       call unpack(gout(num2ejk+num2ejk+1),lab1632,
_IFN1(iv)     +             integ,numlabjk)
_IFN1(iv)       iword=1
_IF1(iv)       call upak6v(gout(num2ejk+num2ejk+1),ij205)
       do 360 m=1,nint
_IFN1(iv)       nij1=integ(iword  )
_IFN1(iv)       nkl1=integ(iword+1)
_IF1(iv)       nij1=ij205(m)
_IF1(iv)       nkl1=kl205(m)
       valp1 = gout(m)
       valk1 = gout(num2ejk+m)
       dump = valp1*da(nkl1)
       dumk = valk1*db(nkl1)
       fa(nij1) = fa(nij1)+dump-dumk
       fb(nij1) = fb(nij1)+dump+dumk
       dump = valp1*da(nij1)
       dumk = valk1*db(nij1)
       fa(nkl1) = fa(nkl1)+dump-dumk
       fb(nkl1) = fb(nkl1)+dump+dumk
_IFN1(iv)       iword=iword+2
 360   continue
c
       if(jblock)2004,2000,2004
      endif
c
 2000 continue
_IFN1(iv)    5 continue
      do   6 m = 1,num
      nij = ikyp(m)
      da(nij) = da(nij)+da(nij)
    6 db(nij) = db(nij)+db(nij)
_IF1(iv)    4 continue
      do 420 m = 1,nx
      duma = da(m)
      dumb = db(m)
      da(m) = (duma+dumb)*pt5
 420  db(m) = (duma-dumb)*pt5
c
      call dscal(nx,pt5,fa,1)
      call dscal(nx,pt5,fb,1)
_IF(parallel)
c***   ***node-MPP***
c...    gather fock-matrices and distribute (watch double size)
      call pg_dgop(700,fa,nx,'+')
      call pg_dgop(701,fb,nx,'+')
c***   ***node-MPP***
_ENDIF
c
      return
      end
_IFN(ibm,vax)
      subroutine gsupam(h,p,g,integ,mword)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IFN(cray,ksr)
      integer *2 integ
      dimension h(*),p(*),g(*),integ(2,*),mword(*)
      do 1 iw=1,mword(1)
      ij=integ(1,iw)
      kl=integ(2,iw)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 continue
_ELSE
      common/craypk/ijkl(1360)
      dimension h(*),p(*),g(*),integ(*)
      call unpack(integ,lab1632,ijkl,numlabjk)
      iword = 1
      do 1 iw=1,mword
      ij=ijkl(iword)
      kl=ijkl(iword+1)
      h(ij)=p(kl)*g(iw)+h(ij)
      h(kl)=p(ij)*g(iw)+h(kl)
    1 iword = iword+2
_ENDIF
      return
      end
      subroutine gsupm(f,d,g,integ,mword)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
_IFN(cray,ksr)
      integer *2 integ
      dimension g(*),integ(2,*),mword(*)
      dimension f(*),d(*)
      do 1 iw=1,mword(1)
      ij=integ(1,iw)
      kl=integ(2,iw)
      f(ij)=d(kl)*g(iw)+f(ij)
      f(kl)=d(ij)*g(iw)+f(kl)
    1 continue
_ELSE
      common/craypk/ijkl(1360)
INCLUDE(common/atmblk)
      dimension g(*),integ(*),f(*),d(*)
      iword = 1
      call unpack(integ,lab1632,ijkl,numlabp)
      do 1 iw=1,mword
      ij=ijkl(iword)
      kl=ijkl(iword+1)
      f(ij)=d(kl)*g(iw)+f(ij)
      f(kl)=d(ij)*g(iw)+f(kl)
    1 iword = iword+2
_ENDIF
      return
      end

_IFN(cray)
_IFN1(at)      subroutine gsuppm(h,p,g,ijkl,mword)
_IF1(at)      subroutine gsuppm(h,p,g,gij,mword)
      implicit REAL  (a-h,o-z)
_IFN1(kx)      integer*4 i8temp(2),fiv11
_IF1(x)      integer *8 fiv11,i8temp
_IF1(k)      integer fiv11,i8temp
_IFN1(at)      integer *2 ijkl
_IFN1(at)      dimension h(*),p(*),g(*),ijkl(2,*),mword(2)
_IF1(at)      dimension h(*),p(*),g(*),gij(*),mword(2)
_IF1(at)      common/craypk/ij205(340),kl205(340),ipad(680)
_IF1(at)      dimension iballs(680)
      equivalence (i8temp,r8temp)
      data fiv11/511/

      if(mword(1).le.0) return
c
_IF1(at)c
_IF1(at)c need kl205 with unit stride to use spdot/spaxpy
_IF1(at)      call unpack(gij,16,iballs,680)
_IF1(at)      iword = 1
_IF1(at)      do 19 i = 1,mword(1)
_IF1(at)          ij205(i) = iballs(iword)
_IF1(at)          kl205(i) = iballs(iword+1)
_IF1(at)19        iword = iword + 2
c
      iwlo = 1
 10   continue
      r8temp=g(iwlo)
c
_IF(convex,ksr)
      nij=IAND64(fiv11,i8temp)       
_ELSEIF(littleendian)
      nij=IAND32(fiv11,i8temp(1))  
_ELSE
      nij=IAND32(fiv11,i8temp(2))
_ENDIF
_IF(alliant,titan)
      ij = ij205(iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword(1)) call caserr('invalid length in gsuppm')
      if(ij.ne.ij205(iwhi)) call caserr('ij mismatch in gsuppm')
_ELSE
      ij = ijkl(1,iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword(1)) call caserr('invalid length in gsuppm')
      if(ij.ne.ijkl(1,iwhi)) call caserr('ij mismatch in gsuppm')
_ENDIF
c
c fortran source
c
       temp=0.0d0
       pij = p(ij)
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
       do 20 iw = iwlo,iwhi
_IF1(at)           kl = kl205(iw)
_IFN1(at)           kl = ijkl(2,iw)

           h(kl) = h(kl) + pij * g(iw)
 20        temp = temp + p(kl) * g(iw)
       h(ij) = h(ij) + temp
c
      iwlo = iwhi + 1
      if(iwlo.le.mword(1)) goto 10
      return
      end
_ENDIF

_IF(cray)
      subroutine gsuppm(h,p,g,iij,mword)
      implicit REAL  (a-h,o-z)
      dimension h(*),p(*),g(*),iij(*)
      common/craypk/ij205(340),kl205(340),ipad(680)
      dimension iballs(680)
      equivalence (i8temp,r8temp)
      data fiv11/511/
      if(mword.le.0) return
c
c
c need kl205 with unit stride to use spdot/spaxpy
      call unpack(iij,16,iballs,680)
      iword = 1
      do 19 i = 1,mword
          ij205(i) = iballs(iword)
          kl205(i) = iballs(iword+1)
19        iword = iword + 2
c
      iwlo = 1
 10   continue
      nij = and(511,g(iwlo))
      ij = ij205(iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword) call caserr('invalid length in gsuppm')
      if(ij.ne.ij205(iwhi)) call caserr('ij mismatch in gsuppm')
c
c cray cal calls ... no faster than cft on xmp. need on 1s.
c
      h(ij) = h(ij) + spdot(nij,p,kl205(iwlo),g(iwlo))
      call spaxpy(nij,p(ij),g(iwlo),h,kl205(iwlo))
c
c
      iwlo = iwhi + 1
      if(iwlo.le.mword) goto 10
      return
      end
_ENDIF

      subroutine gsupum(fa,fb,da,db,goutp,goutk,integ,nint)
      implicit REAL  (a-h,p-w),integer  (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IFN(cray,ksr)
      integer *2 integ
      dimension nint(*)
_ENDIF
_IF(cray,ksr)
      common/craypk/ijkl(1360)
_ENDIF
      dimension goutp(*),goutk(*)
      dimension fa(*),fb(*),da(*),db(*),integ(*)
c
       iword=1
_IFN(cray,ksr)
       do 360 m=1,nint(1)
       nij1=integ(iword  )
       nkl1=integ(iword+1)
_ENDIF
_IF(cray,ksr)
       call unpack(integ,lab1632,ijkl,numlabjk)
       do 360 m=1,nint
       nij1=ijkl(iword  )
       nkl1=ijkl(iword+1)
_ENDIF
       valp1=goutp(m)
       valk1 = goutk(m)
       dump = valp1*da(nkl1)
       dumk = valk1*db(nkl1)
       fa(nij1) = fa(nij1)+dump-dumk
       fb(nij1) = fb(nij1)+dump+dumk
       dump = valp1*da(nij1)
       dumk = valk1*db(nij1)
       fa(nkl1) = fa(nkl1)+dump-dumk
       fb(nkl1) = fb(nkl1)+dump+dumk
       iword=iword+2
 360   continue
       return
       end



      subroutine jksupm(hscm0,hscm1,hscm2,dscm0,dscm1,nscmf,l2,
     + goutp, goutk, integ, nint)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/atmblk)
_IFN(cray,ksr)
      integer *2 integ
_ENDIF
_IF(cray,ksr)
      common/craypk/ijkl(1360)
_ENDIF
      dimension hscm0(*),hscm1(*),hscm2(*),dscm0(*),dscm1(*)
      dimension goutp(*),goutk(*)
_IFN1(ck)      dimension nint(2),integ(2,*)
_IF1(ck)      dimension integ(*)
c
c---------both closed + open shells------------
c
_IF(cray,ksr)
      iword=1
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
      call unpack(integ,lab1632,ijkl,numlabjk)
_ENDIF
      if(nscmf.eq.1) then
c
_IF(cray,ksr)
       do   10 m = 1,nint
       ij = ijkl(iword  )
       kl = ijkl(iword+1)
_ELSE
       do   10 m = 1,nint(1)
       ij = integ(1,m)
       kl = integ(2,m)
_ENDIF
       goutpm = goutp(m)
       goutkm = goutk(m)
       hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
       hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
       hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
       hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
       hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
       hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
_IFN1(ck)   10  continue
_IF1(ck)   10  iword = iword + 2
      else
_IF(cray,ksr)
       do   1 m = 1,nint
       ij = ijkl(iword  )
       kl = ijkl(iword+1)
_ELSE
       do   1 m = 1,nint(1)
       ij = integ(1,m)
       kl = integ(2,m)
_ENDIF
       goutpm = goutp(m)
       goutkm = goutk(m)
       hscm0(ij) = hscm0(ij) + goutpm*dscm0(kl)
       hscm0(kl) = hscm0(kl) + goutpm*dscm0(ij)
_IF1(x)c$dir scalar
       do 2   i = 1,nscmf
       hscm1(ij) = hscm1(ij) + goutpm*dscm1(kl)
       hscm2(ij) = hscm2(ij) + goutkm*dscm1(kl)
       hscm1(kl) = hscm1(kl) + goutpm*dscm1(ij)
       hscm2(kl) = hscm2(kl) + goutkm*dscm1(ij)
       ij = ij + l2
    2  kl = kl + l2
_IFN1(ck)    1  continue
_IF1(ck)    1  iword = iword + 2
      endif
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end

_IF(ccpdft)
_IF(upck-equiv)
      subroutine proc2m(fock,p,ak,q,facex,ocoul,oexch,
     +                  gg,integ,mword)
_ELSE
      subroutine proc2m(fock,p,ak,q,facex,ocoul,oexch,
     +                  gg,gijkl,mword)
_ENDIF
_ELSE
_IF(upck-equiv)
      subroutine proc2m(fock,p,ak,q,gg,integ,mword)
_ELSE
      subroutine proc2m(fock,p,ak,q,gg,gijkl,mword)
_ENDIF
_ENDIF

      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension ak(*),q(*),p(*),fock(*),gg(*)

_IFN(cray,convex,ksr,t3d)
      dimension mword(*)
_ENDIF
_IF(upck-equiv)
      logical *1 integ
      dimension integ(*)
_ELSE
      common/craypk/integ(1360)
INCLUDE(common/atmblk)
      dimension gijkl(*)
_ENDIF
_IF(ccpdft)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
_ENDIF

_IF(upck-equiv)
c  variant uses logical*1 on RHS of unpack
_IF(i8)
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
_ELSE
      logical*1 ilog(4),jlog(4),klog(4),llog(4)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
c
      iword=1
_IFN(upck-equiv)
      call unpack(gijkl,lab816,integ,numlab)
_IF(cray,convex,ksr,t3d)
      do 6000 iw=1,mword
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
       jlog(1) = integ(iword  )
       ilog(1) = integ(iword+1)
       llog(1) = integ(iword+2)
       klog(1) = integ(iword+3)
_ELSEIF(i8)
       ilog(8) = integ(iword  )
       jlog(8) = integ(iword+1)
       klog(8) = integ(iword+2)
       llog(8) = integ(iword+3)
_ELSE
       ilog(4) = integ(iword  )
       jlog(4) = integ(iword+1)
       klog(4) = integ(iword+2)
       llog(4) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
       j = integ(iword  )
       i = integ(iword+1)
       l = integ(iword+2)
       k = integ(iword+3)
_ELSE
       i = integ(iword  )
       j = integ(iword+1)
       k = integ(iword+2)
       l = integ(iword+3)
_ENDIF
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
_ELSE
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
_ENDIF
 6000 iword=iword+4
      return
      end

_IF(ccpdft)
_IF(upck-equiv)
      subroutine proc2m_255(fock,p,ak,q,facex,ocoul,
     +                  oexch,
     +                  gg,integ,mword)
_ELSE
      subroutine proc2m_255(fock,p,ak,q,facex,ocoul,
     +                  oexch,
     +                  gg,gijkl,mword)
_ENDIF
_ELSE
_IF(upck-equiv)
      subroutine proc2m_255(fock,p,ak,q,gg,integ,mword)
_ELSE
      subroutine proc2m_255(fock,p,ak,q,gg,gijkl,mword)
_ENDIF
_ENDIF

      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension ak(*),q(*),p(*),fock(*),gg(*)

_IFN(cray,convex,ksr,t3d)
      dimension mword(*)
_ENDIF
_IF(upck-equiv)
      integer *2 integ
      dimension integ(*)
_ELSE
      common/craypk/integ(1360)
INCLUDE(common/atmblk)
      dimension gijkl(*)
_ENDIF
_IF(ccpdft)
      logical ocoul,oexch
INCLUDE(common/ccpdft.hf77)
_ENDIF

_IF(upck-equiv)
c  variant uses logical*1 on RHS of unpack
_IF(i8)
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
_ELSE
      integer *2 ilog(2),jlog(2),klog(2),llog(2)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
c
      iword=1
_IFN(upck-equiv)
      call unpack(gijkl,lab816,integ,numlab)
_IF(cray,convex,ksr,t3d)
      do 6000 iw=1,mword
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
       jlog(1) = integ(iword  )
       ilog(1) = integ(iword+1)
       llog(1) = integ(iword+2)
       klog(1) = integ(iword+3)
_ELSEIF(i8)
       ilog(4) = integ(iword  )
       jlog(4) = integ(iword+1)
       klog(4) = integ(iword+2)
       llog(4) = integ(iword+3)
_ELSE
       ilog(2) = integ(iword  )
       jlog(2) = integ(iword+1)
       klog(2) = integ(iword+2)
       llog(2) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
       j = integ(iword  )
       i = integ(iword+1)
       l = integ(iword+2)
       k = integ(iword+3)
_ELSE
       i = integ(iword  )
       j = integ(iword+1)
       k = integ(iword+2)
       l = integ(iword+3)
_ENDIF
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
c
c coulomb term
c
      if(ocoul)then
        fock(ij)=g4*p(kl)+fock(ij)
        if(ij.ne.kl) then
         fock(kl)=g4*p(ij)+fock(kl)
        endif
      endif
c
c exchange
c
      if(oexch)then
c
c full term for HF or weighted term for b3lyp etc
c
        gik=gik*facex
        g2 = gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
1       ajk=fock(jk)-gil*p(il)
        bjk=ak(jk)+gil*q(il)
        ail=fock(il)-gil*p(jk)
        bil=ak(il)+gil*q(jk)
        aik=fock(ik)-gik*p(jl)
        bik=ak(ik)+gik*q(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        ak(jl)=ak(jl)+gik*q(ik)
        fock(jk)=ajk
        ak(jk)=bjk
        fock(il)=ail
        ak(il)=bil
        fock(ik)=aik
        ak(ik)=bik
      endif
_ELSE
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
_ENDIF
 6000 iword=iword+4
      return
      end

_IF(upck-equiv)
      subroutine sgmatn(l2,nshell,b,c,p,gg,intin,mword)
_ELSE
      subroutine sgmatn(l2,nshell,b,c,p,gg,gijkl,mword)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension b(*),c(*),p(*),gg(*)
_IFN(cray,convex,ksr,t3d)
      dimension mword(*)
_ENDIF
_IF(upck-equiv)
      logical *1 intin
      dimension intin(*)
_ELSE
INCLUDE(common/atmblk)
      common/craypk/intin(1360)
      dimension gijkl(*)
_ENDIF
_IF(upck-equiv)
_IF(i8)
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
_ELSE
      logical*1 ilog(4),jlog(4),klog(4),llog(4)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
_IFN(upck-equiv)
      call unpack(gijkl,lab816,intin,numlab)
_ENDIF
      int4=1
      if(nshell.gt.1) then
c
_IF(cray,convex,ksr,t3d)
      do 6 iw=1,mword
_ELSE
      do 6 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
       jlog(1)=intin(int4  )
       ilog(1)=intin(int4+1)
       llog(1)=intin(int4+2)
       klog(1)=intin(int4+3)
_ELSEIF(i8)
       ilog(8)=intin(int4  )
       jlog(8)=intin(int4+1)
       klog(8)=intin(int4+2)
       llog(8)=intin(int4+3)
_ELSE
       ilog(4)=intin(int4  )
       jlog(4)=intin(int4+1)
       klog(4)=intin(int4+2)
       llog(4)=intin(int4+3)
_ENDIF
_ELSEIF(littleendian)
       j=intin(int4  )
       i=intin(int4+1)
       l=intin(int4+2)
       k=intin(int4+3)
_ELSE
       i=intin(int4  )
       j=intin(int4+1)
       k=intin(int4+2)
       l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
_IF1(x)c$dir scalar
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
    6 int4=int4+4
      else
_IF(cray,convex,t3d,ksr)
      do 66 iw=1,mword
_ELSE
      do 66 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
      jlog(1)=intin(int4  )
      ilog(1)=intin(int4+1)
      llog(1)=intin(int4+2)
      klog(1)=intin(int4+3)
_ELSEIF(i8)
      ilog(8)=intin(int4  )
      jlog(8)=intin(int4+1)
      klog(8)=intin(int4+2)
      llog(8)=intin(int4+3)
_ELSE
      ilog(4)=intin(int4  )
      jlog(4)=intin(int4+1)
      klog(4)=intin(int4+2)
      llog(4)=intin(int4+3)
_ENDIF
_ELSEIF(littleendian)
      j=intin(int4  )
      i=intin(int4+1)
      l=intin(int4+2)
      k=intin(int4+3)
_ELSE
      i=intin(int4  )
      j=intin(int4+1)
      k=intin(int4+2)
      l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
   66 int4=int4+4
      endif
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end
_IF(upck-equiv)
      subroutine sgmatn_255(l2,nshell,b,c,p,gg,intin,mword)
_ELSE
      subroutine sgmatn_255(l2,nshell,b,c,p,gg,gijkl,mword)
_ENDIF
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
      dimension b(*),c(*),p(*),gg(*)
_IFN(cray,convex,ksr,t3d)
      dimension mword(*)
_ENDIF
_IF(upck-equiv)
      integer *2 intin
      dimension intin(*)
_ELSE
INCLUDE(common/atmblk)
      common/craypk/intin(1360)
      dimension gijkl(*)
_ENDIF
_IF(upck-equiv)
_IF(i8)
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
_ELSE
      integer *2 ilog(2),jlog(2),klog(2),llog(2)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
_IF1(c)cdir$ list
_IF1(c)cdir$ novector
_IFN(upck-equiv)
      call unpack(gijkl,lab816,intin,numlab)
_ENDIF
      int4=1
      if(nshell.gt.1) then
c
_IF(cray,convex,ksr,t3d)
      do 6 iw=1,mword
_ELSE
      do 6 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
       jlog(1)=intin(int4  )
       ilog(1)=intin(int4+1)
       llog(1)=intin(int4+2)
       klog(1)=intin(int4+3)
_ELSEIF(i8)
       ilog(4)=intin(int4  )
       jlog(4)=intin(int4+1)
       klog(4)=intin(int4+2)
       llog(4)=intin(int4+3)
_ELSE
       ilog(2)=intin(int4  )
       jlog(2)=intin(int4+1)
       klog(2)=intin(int4+2)
       llog(2)=intin(int4+3)
_ENDIF
_ELSEIF(littleendian)
       j=intin(int4  )
       i=intin(int4+1)
       l=intin(int4+2)
       k=intin(int4+3)
_ELSE
       i=intin(int4  )
       j=intin(int4+1)
       k=intin(int4+2)
       l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
    1 continue
_IF1(x)c$dir scalar
      do 2 iiii=1,nshell
      bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
      ij=ij+l2
      kl=kl+l2
      ik=ik+l2
      il=il+l2
      jk=jk+l2
    2 jl=jl+l2
    6 int4=int4+4
      else
_IF(cray,convex,t3d,ksr)
      do 66 iw=1,mword
_ELSE
      do 66 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
      jlog(1)=intin(int4  )
      ilog(1)=intin(int4+1)
      llog(1)=intin(int4+2)
      klog(1)=intin(int4+3)
_ELSEIF(i8)
      ilog(4)=intin(int4  )
      jlog(4)=intin(int4+1)
      klog(4)=intin(int4+2)
      llog(4)=intin(int4+3)
_ELSE
      ilog(2)=intin(int4  )
      jlog(2)=intin(int4+1)
      klog(2)=intin(int4+2)
      llog(2)=intin(int4+3)
_ENDIF
_ELSEIF(littleendian)
      j=intin(int4  )
      i=intin(int4+1)
      l=intin(int4+2)
      k=intin(int4+3)
_ELSE
      i=intin(int4  )
      j=intin(int4+1)
      k=intin(int4+2)
      l=intin(int4+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 11
      jk=ikyk+j
      if(j.ge.l)goto 11
      jl=iky(l)+j
 11   bij=g2*p(kl)+b(ij)
      b(kl)=g2*p(ij)+b(kl)
      b(ij)=bij
      bjk=c(jk)+gil*p(il)
      bil=c(il)+gil*p(jk)
      bik=c(ik)+gik*p(jl)
      c(jl)=c(jl)+gik*p(ik)
      c(jk)=bjk
      c(il)=bil
      c(ik)=bik
   66 int4=int4+4
      endif
_IF1(c)cdir$ vector
_IF1(c)cdir$ nolist
      return
      end
_IFN(cray,convex,titan)
_IF(upck-equiv)
      subroutine sgmtmm(fock,p,gg,integ,mword)
_ELSE
      subroutine sgmtmm(fock,p,gg,gijkl,mword)
_ENDIF
      implicit REAL  (a-h,o-z)
      dimension p(*),fock(*),gg(*)
_IFN(t3d)
      dimension mword(*)
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c
_IF(upck-equiv) 
      logical *1 integ 
      dimension integ(*)
_ELSE 
      common/craypk/integ(1360) 
INCLUDE(common/atmblk)
      dimension gijkl(*)
_ENDIF 
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
_IF(upck-equiv)
_IF(i8)
      logical*1 ilog(8),jlog(8),klog(8),llog(8)
_ELSE
      logical*1 ilog(4),jlog(4),klog(4),llog(4)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
c
      iword = 1
_IF(ccpdft)
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
_IFN(upck-equiv)
      call unpack(gijkl,lab816,integ,numlab)
_IF(t3d)
      do 6000 iw=1,mword
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
      jlog(1) = integ(iword  )
      ilog(1) = integ(iword+1)
      llog(1) = integ(iword+2)
      klog(1) = integ(iword+3)
_ELSEIF(i8)
      ilog(8) = integ(iword  )
      jlog(8) = integ(iword+1)
      klog(8) = integ(iword+2)
      llog(8) = integ(iword+3)
_ELSE
      ilog(4) = integ(iword  )
      jlog(4) = integ(iword+1)
      klog(4) = integ(iword+2)
      llog(4) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
_ELSE
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
_ENDIF
 6000 iword=iword+4
      return
      end
_ENDIF
_IF(upck-equiv)
      subroutine sgmtmm_255(fock,p,gg,integ,mword)
_ELSE
      subroutine sgmtmm_255(fock,p,gg,gijkl,mword)
_ENDIF
      implicit REAL  (a-h,o-z)
      dimension p(*),fock(*),gg(*)
_IFN(t3d)
      dimension mword(*)
_ENDIF
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c
_IF(upck-equiv) 
      integer *2 integ 
      dimension integ(*)
_ELSE 
      common/craypk/integ(1360) 
INCLUDE(common/atmblk)
      dimension gijkl(*)
_ENDIF 
_IF(ccpdft)
      logical oexch, ocoul, onodft
INCLUDE(common/ccpdft.hf77)
_ENDIF
c
_IF(upck-equiv)
_IF(i8)
      integer *2 ilog(4),jlog(4),klog(4),llog(4)
_ELSE
      integer *2 ilog(2),jlog(2),klog(2),llog(2)
_ENDIF
      equivalence (i,ilog(1)),(j,jlog(1))
      equivalence (k,klog(1)),(l,llog(1))
      data i,j,k,l/4*0/
_ENDIF
c
      iword = 1
_IF(ccpdft)
      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
_ENDIF
_IFN(upck-equiv)
      call unpack(gijkl,lab816,integ,numlab)
_IF(t3d)
      do 6000 iw=1,mword
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_ELSE
      do 6000 iw=1,mword(1)
_ENDIF
_IF(upck-equiv)
_IF(littleendian)
      jlog(1) = integ(iword  )
      ilog(1) = integ(iword+1)
      llog(1) = integ(iword+2)
      klog(1) = integ(iword+3)
_ELSEIF(i8)
      ilog(4) = integ(iword  )
      jlog(4) = integ(iword+1)
      klog(4) = integ(iword+2)
      llog(4) = integ(iword+3)
_ELSE
      ilog(2) = integ(iword  )
      jlog(2) = integ(iword+1)
      klog(2) = integ(iword+2)
      llog(2) = integ(iword+3)
_ENDIF
_ELSEIF(littleendian)
      j = integ(iword  )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)
_ELSE
      i = integ(iword  )
      j = integ(iword+1)
      k = integ(iword+2)
      l = integ(iword+3)
_ENDIF
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
_IF(ccpdft)
c
c coulomb term
c
      if(ocoul)then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
      if(onodft)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        gik = gik*facex
        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
_ELSE
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      ail=fock(il)-gil*p(jk)
      aik=fock(ik)-gik*p(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      fock(jk)=ajk
      fock(il)=ail
      fock(ik)=aik
_ENDIF
 6000 iword=iword+4
      return
      end
_ENDIF
      subroutine ver_util2(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util2.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
_IF(oldpack)
      subroutine rinsrt(r,m)
      implicit REAL  (a-h,o-z)
c ***
c *** insert last len=9 bits of m into r at position ipos=55.
c ***
_IF1(c)      r=or(and(compl(511),r),and(511,m))
_IF1(x)      integer *8 r,cfiv11,fiv11,m8temp
_IF1(x)      dimension m4temp(2)
_IF1(x)      equivalence (m8temp,m4temp(1))
_IF1(x)      data cfiv11,fiv11/-512,511/
_IF1(x)      m4temp(1)=0
_IF1(x)      m4temp(2)=m
_IF1(x)      r=ior(iand(cfiv11,r),iand(fiv11,m8temp))
_IF1(k)      integer r,cfiv11,fiv11,m8temp
_IF1(k)      integer *4  m4temp(2)
_IF1(k)      equivalence (m8temp,m4temp(1))
_IF1(k)      data cfiv11,fiv11/-512,511/
_IF1(k)      m4temp(1)=0
_IF1(k)      m4temp(2)=m
_IF1(k)      r=or(and(cfiv11,r),and(fiv11,m8temp))
_IF1(adt)      integer*4 r(2),m
_IF1(t)      call mvbits(m,0,9,r(2),0)
_IF1(ad)      call mvbits(m,0,9,r(1),0)
_IF1(i)      dimension m8temp(2)
_IF1(i)      equivalence (temp,m8temp(1))
_IF1(i)      data cfiv11,fiv11/zfffffffffffffe00,
_IF1(i)     *                  z1ff/
_IF1(i)      data m8temp/0,0/
_IF1(i)      m8temp(2)=m
_IF1(i)      r=or(and(cfiv11,r),and(fiv11,temp))
_IF1(pgbrs)      integer r
_IF1(pgbrs)      dimension r(2)
_IF1(pgbrs)      data mtemp/0/
_IF1(pgrs)      r(2)=and(-512,r(2))
_IF1(pgrs)      mtemp=and(511,m)
_IF1(pgrs)      r(2)=or(r(2),mtemp)
_IF1(b)      r(2)=iand(-512,r(2))
_IF1(b)      mtemp=iand(511,m)
_IF1(b)      r(2)=ior(r(2),mtemp)
_IF1(h)      dimension m8temp(2)
_IF1(h)      equivalence (temp,m8temp(1))
_IF1(h)      temp = r
_IF1(h)      m8temp(1) = or(and(m8temp(1),not(511)),m)
_IF1(h)      r = temp
      return
      end
      function locatz(label,nf,itext)
      implicit REAL  (a-h,p-z),integer    (i-n)
_IFN1(cxk)      logical eq
_IF1(x)      integer *8 label,itext
_IF1(tpgbdrhs)      integer label(2,*),itext(2)
_IF1(aiv)      REAL  label,itext
_IFN1(tpgbdrhs)      dimension label(*)
      do 1 i=1,nf
_IF1(cxk)      if(label(i).eq.itext)go to 2
_IF1(iva)      if(eq(label(i),itext))go to 2
_IF1(tpgbdrhs)      if(eq(label(1,i),itext))go to 2
 1    continue
      locatz=0
      return
 2    locatz=i
      return
      end
c
c gsuppm - original version
c
_IFN1(at)      subroutine gsuppm(h,p,g,ijkl,mword)
_IF1(at)      subroutine gsuppm(h,p,g,gij,mword)
      implicit REAL  (a-h,o-z)
_IFN1(kx)      integer i8temp(2),fiv11
_IF1(x)      integer *8 fiv11,i8temp
_IF1(k)      integer fiv11,i8temp
_IFN1(at)      integer *2 ijkl
_IFN1(at)      dimension h(*),p(*),g(*),ijkl(2,*),mword(2)
_IF1(at)      dimension h(*),p(*),g(*),gij(*),mword(2)
_IF1(at)      common/craypk/ij205(340),kl205(340),ipad(680)
_IF1(at)      dimension iballs(680)
      equivalence (i8temp,r8temp)
      data fiv11/511/
      if(mword(1).le.0) return
c
_IF1(at)c
_IF1(at)c need kl205 with unit stride to use spdot/spaxpy
_IF1(at)      call unpack(gij,16,iballs,680)
_IF1(at)      iword = 1
_IF1(at)      do 19 i = 1,mword(1)
_IF1(at)          ij205(i) = iballs(iword)
_IF1(at)          kl205(i) = iballs(iword+1)
_IF1(at)19        iword = iword + 2
c
      iwlo = 1
 10   continue
      r8temp=g(iwlo)
_IF1(x)      nij=iand(fiv11,i8temp)
_IF1(k)      nij=and(fiv11,i8temp)
_IF1(bt)      nij=iand(fiv11,i8temp(2))
_IF1(ad)      nij=iand(fiv11,i8temp(1))
_IF1(pgsr)      nij=and(fiv11,i8temp(2))
_IF1(h)      nij=and(fiv11,i8temp(1))
_IF1(at)      ij = ij205(iwlo)
_IFN1(at)      ij = ijkl(1,iwlo)
      iwhi = iwlo + nij - 1
      if(iwhi.gt.mword(1)) call caserr('invalid length in gsuppm')
_IF1(at)      if(ij.ne.ij205(iwhi)) call caserr('ij mismatch in gsuppm')
_IFN1(at)      if(ij.ne.ijkl(1,iwhi)) call caserr('ij mismatch in gsuppm')
c
c fortran source
c
       temp=0.0d0
       pij = p(ij)
_IF1(x)c$dir no_recurrence
_IF1(a)cvd$  nodepchk
       do 20 iw = iwlo,iwhi
_IF1(at)           kl = kl205(iw)
_IFN1(at)           kl = ijkl(2,iw)
           h(kl) = h(kl) + pij * g(iw)
 20        temp = temp + p(kl) * g(iw)
       h(ij) = h(ij) + temp
c
      iwlo = iwhi + 1
      if(iwlo.le.mword(1)) goto 10
      return
      end

c
cpack proc2m original
c
_IFN1(ck)      subroutine proc2m(fock,p,ak,q,gg,integ,mword)
_IF1(ck)      subroutine proc2m(fock,p,ak,q,gg,ijkl,mword)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
_IFN(cray,ksr)
      logical *1 integ
      dimension ak(*),q(*),p(*),fock(*),gg(*),integ(*),mword(2)
_ENDIF
_IF(cray,ksr)
      common/craypk/integ(1)
      dimension ak(*),q(*),p(*),fock(*),gg(*),ijkl(*)
_ENDIF
INCLUDE(common/mapper)
c
_IFN1(ck)      data i,j,k,l/4*0/
_IF1(r)      logical*1 ilog(4),jlog(4),klog(4),llog(4)
_IF1(r)      equivalence (i,ilog(1)),(j,jlog(1))
_IF1(r)      equivalence (k,klog(1)),(l,llog(1))
c
      iword=1
_IF1(ck)      call unpack(ijkl,8,integ,1360)
_IF1(ck)      do 6000 iw=1,mword
_IFN1(ck)      do 6000 iw=1,mword(1)
_IFN1(adhr)       i = integ(iword  )
_IFN1(adhr)       j = integ(iword+1)
_IFN1(adhr)       k = integ(iword+2)
_IFN1(adhr)       l = integ(iword+3)
_IF1(r)       ilog(4) = integ(iword  )
_IF1(r)       jlog(4) = integ(iword+1)
_IF1(r)       klog(4) = integ(iword+2)
_IF1(r)       llog(4) = integ(iword+3)
_IF1(adh)       j = integ(iword  )
_IF1(adh)       i = integ(iword+1)
_IF1(adh)       l = integ(iword+2)
_IF1(adh)       k = integ(iword+3)
      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
      aij=g4*p(kl)+fock(ij)
      fock(kl)=g4*p(ij)+fock(kl)
      fock(ij)=aij
c... exchange
      gil=gik
      if(i.eq.k.or.j.eq.l)gik=g2
      if(j.eq.k)gil=g2
      if(j.ge.k)goto 1
      jk=ikyk+j
      if(j.ge.l)goto 1
      jl=iky(l)+j
1     ajk=fock(jk)-gil*p(il)
      bjk=ak(jk)+gil*q(il)
      ail=fock(il)-gil*p(jk)
      bil=ak(il)+gil*q(jk)
      aik=fock(ik)-gik*p(jl)
      bik=ak(ik)+gik*q(jl)
      fock(jl)=fock(jl)-gik*p(ik)
      ak(jl)=ak(jl)+gik*q(ik)
      fock(jk)=ajk
      ak(jk)=bjk
      fock(il)=ail
      ak(il)=bil
      fock(ik)=aik
      ak(ik)=bik
 6000 iword=iword+4
      return
      end
_ENDIF
