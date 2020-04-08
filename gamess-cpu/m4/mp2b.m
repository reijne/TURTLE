c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp2b.m,v $
c  $State: Exp $
c  
      subroutine umpmkz(y1,y2,zl,a1,a2,b,nocc,ncoorb,mn,
     1 istrma,istrmb,ibly1,ibly2,iblw)
      implicit REAL  (a-h,o-z)
c
c  construct rest of zl matrix
c
      dimension y1(ncoorb,ncoorb),y2(ncoorb,ncoorb),zl(ncoorb,ncoorb),
     1 a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),b(mn)
c
      ifile1 = 1
      nocc1 = nocc + 1
      nsq = ncoorb*ncoorb
      call rdedx(y1,nsq,ibly1,ifile1)
      call rdedx(y2,nsq,ibly2,ifile1)
      call rdedx(zl,nsq,iblw,ifile1)
      call rewedz(istrma)
      call rewedz(istrmb)
c
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
      tempbl = 0.0d0
      do 90 ki = 1 , ncoorb
         do 80 kj = 1 , ncoorb
            tempbl = tempbl + dabs(zl(ki,kj))
 80      continue
 90   continue
      return
      end
      subroutine cuhf2(eps,ea,eb,amaa,amab,amba,ambb,
     1 nocca,noccb,ncoorb,mna,mn,istmaa,istmab,istmba,istmbb,iblw,ifw)
      implicit REAL  (a-h,o-z)
_IFN1(iv)      common/craypk/labout(1360)
_IF1(iv)      common/craypk/labij(340),labkl(340)
INCLUDE(common/atmblk)
      common/blkin/g(510),nint,nxtr
      dimension amaa(ncoorb,ncoorb),amab(ncoorb,ncoorb),
     1 amba(ncoorb,ncoorb),ambb(ncoorb,ncoorb),ea(ncoorb),
     1 eb(ncoorb),eps(mn)
      logical lstop,skipp
      dimension skipp(100)
      data m0/0/
c
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
      iai = 0
      do 50 ia = noccb + 1 , ncoorb
         do 40 i = 1 , noccb
            iai = iai + 1
            eps(mna+iai) = 1.0d0/(eb(ia)-eb(i))
 40      continue
 50   continue
c
      tempep = 0.0d0
      do 60 ki = 1 , mn
         tempep = tempep + dabs(eps(ki))
 60   continue
      call search(iblw,ifw)
      call rewedz(istmaa)
      call rewedz(istmab)
      call rewedz(istmba)
      call rewedz(istmbb)
      do 140 ip = 1 , ncoorb
         do 130 iq = 1 , ip
            call rdedz(amaa,nsq,istmaa)
            call rdedz(amab,nsq,istmab)
            call rdedz(amba,nsq,istmba)
            call rdedz(ambb,nsq,istmbb)
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
                  iai = (ip-nocca-1)*nocca + iq
                  ibj = 0
                  do 80 ib = nocca + 1 , ncoorb
                     do 70 j = 1 , nocca
                        ibj = ibj + 1
                        if (iai.ge.ibj) then
                           nint = nint + 1
                           g(nint) = amaa(ib,j)
_IFN1(iv)                           nint2 = nint + nint
_IFN1(iv)                           labout(nint2-1) = iai
_IFN1(iv)                           labout(nint2) = ibj
_IF1(iv)                           labij(nint) = iai
_IF1(iv)                           labkl(nint) = ibj
                           if (nint.eq.num2e) then
_IFN1(iv)                           call pack(g(num2ep+1),lab1632,labout,numlabp)
_IF1(iv)                           call pak4v(labij,g(num2ep+1))
                              call put(g,m511,ifw)
                              nint = 0
                           end if
                        end if
 70                  continue
 80               continue
c
                  ibj = 0
                  do 100 ib = noccb + 1 , ncoorb
                     do 90 j = 1 , noccb
                        ibj = ibj + 1
                        nint = nint + 1
                        g(nint) = amab(ib,j)*2.0d0
_IFN1(iv)                        nint2 = nint + nint
_IFN1(iv)                        labout(nint2-1) = iai
_IFN1(iv)                        labout(nint2) = mna + ibj
_IF1(iv)                           labij(nint) = iai
_IF1(iv)                           labkl(nint) = mna + ibj
                        if (nint.eq.num2e) then
_IFN1(iv)                           call pack(g(num2ep+1),lab1632,labout,numlabp)
_IF1(iv)                           call pak4v(labij,g(num2ep+1))
                           call put(g,m511,ifw)
                           nint = 0
                        end if
 90                  continue
 100              continue
               end if
            end if
c
            if (ip.gt.noccb) then
               if (iq.le.noccb) then
                  iai = (ip-noccb-1)*noccb + iq
c
                  ibj = 0
                  do 120 ib = noccb + 1 , ncoorb
                     do 110 j = 1 , noccb
                        ibj = ibj + 1
                        if (iai.ge.ibj) then
                           nint = nint + 1
                           g(nint) = ambb(ib,j)
_IFN1(iv)                           nint2 = nint + nint
_IFN1(iv)                           labout(nint2-1) = mna + iai
_IFN1(iv)                           labout(nint2) = mna + ibj
_IF1(iv)                           labij(nint) = mna + iai
_IF1(iv)                           labkl(nint) = mna + ibj
                           if (nint.eq.num2e) then
_IFN1(iv)                           call pack(g(num2ep+1),lab1632,labout,numlabp)
_IF1(iv)                           call pak4v(labij,g(num2ep+1))
                              call put(g,m511,ifw)
                              nint = 0
                           end if
                        end if
 110                 continue
 120              continue
               end if
            end if
 130     continue
 140  continue
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
      subroutine umpe2(t1,t2,t3,t4,f,g,val3,val4,e2,
     + ifort1,ifort2,ifort3,ifort4,nocca,noccb,ncoorb,
     + ifila)
      implicit REAL  (a-h,o-z)
      dimension t1(ncoorb*ncoorb),t2(ncoorb*ncoorb),f(ncoorb),g(ncoorb)
      dimension t3(ncoorb*ncoorb),t4(ncoorb*ncoorb)
INCLUDE(common/uhfspn)
INCLUDE(common/atmblk)
_IF1()      logical pump2
      common/blkin/yy(510),nint
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
      data quart/0.25d0/
c
      val = 0.0d0
      val1 = 0.0d0
      val2 = 0.0d0
      val3 = 0.0d0
      val4 = 0.0d0
c
      nsq = ncoorb*ncoorb
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      call rewedz(ifort4)
      do 70 ip = 1 , ncoorb
         do 60 iq = 1 , ip
            call rdedz(t1,nsq,ifort1)
            call rdedz(t2,nsq,ifort2)
            call rdedz(t3,nsq,ifort3)
            call rdedz(t4,nsq,ifort4)
c
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
c                 ib = ip
                  j = iq
                  ebj1 = f(ip) - f(iq)
                  do 30 ia = nocca + 1 , ncoorb
                     do 20 i = 1 , nocca
                        iai = (ia-1)*ncoorb + i
                        difen = ebj1 + f(ia) - f(i)
                        ff = 1.0d0/difen
                        val = t1(iai) - t2(iai)
                        val = val*val
                        val1 = val*ff + val1
c
 20                  continue
 30               continue
               end if
            end if
c
            if (ip.gt.noccb) then
               if (iq.le.noccb) then
c                 ib = ip
                  j = iq
                  ebj2 = g(ip) - g(iq)
                  do 50 ia = noccb + 1 , ncoorb
                     do 40 i = 1 , noccb
                        iai = (ia-1)*ncoorb + i
                        difen = ebj2 + g(ia) - g(i)
                        ff = 1.0d0/difen
                        val = t3(iai) - t4(iai)
                        val = val*val
                        val2 = val*ff + val2
 40                  continue
 50               continue
               end if
            end if
 60      continue
 70   continue
c
c     read sequentially down the second list
_IFN1(iv)      call setsto(1360,0,labs)
_IF1(iv)      call setsto(1360,0,i205)
c
      do 100 ipass = 1 , npassm
         iblk1 = ipass*2
         call search(mupblk(iblk1),ifila)
 80      call find(ifila)
         call get(yy,nw)
         if (nw.ne.0) then
c
_IFN1(iv)             call unpack(yy(num2e+1),lab816,labs,numlab)
_IF1(iv)            call upak8v(yy(num2e+1),i205)
c
            do 90 int = 1 , nint
_IFN1(iv)               num = int + int + int + int

_IF(ibm,vax)
               i = i205(int)
               j = j205(int)
               k = k205(int)
               l = l205(int)
_ELSEIF(littleendian)
               i = labs(num-2)
               j = labs(num-3)
               k = labs(num  )
               l = labs(num-1)
_ELSE
               i = labs(num-3)
               j = labs(num-2)
               k = labs(num-1)
               l = labs(num)
_ENDIF
               if (i.gt.noccb) then
                  if (j.le.noccb) then
                     if (k.gt.nocca) then
                        if (l.le.nocca) then
                           val = yy(int)*yy(int)
                           iaa = i
                           ii = j
                           ibb = k
                           jj = l
                           difen = g(iaa) + f(ibb) - g(ii) - f(jj)
                           ff = 1.0d0/difen
                           val3 = val*ff*4 + val3
                        end if
                     end if
                  end if
               end if
c
c
 90         continue
            go to 80
         end if
 100  continue
c
c     read sequentially down third list
c
      do 130 ipass = 1 , npassm
         iblk2 = npassm*2 + ipass*2 - 1
         call search(mupblk(iblk2),ifila)
 110     call find(ifila)
         call get(yy,nw)
         if (nw.ne.0) then
c
_IFN1(iv)             call unpack(yy(num2e+1),lab816,labs,numlab)
_IF1(iv)            call upak8v(yy(num2e+1),i205)
c
            do 120 int = 1 , nint
_IFN1(iv)               num = int + int + int + int
_IF(ibm,vax)
               i = i205(int)
               j = j205(int)
               k = k205(int)
               l = l205(int)
_ELSEIF(littleendian)
               i = labs(num-2)
               j = labs(num-3)
               k = labs(num  )
               l = labs(num-1)
_ELSE
               i = labs(num-3)
               j = labs(num-2)
               k = labs(num-1)
               l = labs(num)
_ENDIF
               if (i.gt.nocca) then
                  if (j.le.nocca) then
                     if (k.gt.noccb) then
                        if (l.le.noccb) then
                           iaa = i
                           ii = j
                           ibb = k
                           jj = l
                           iai = iaa*(iaa-1)/2 + ii
                           ibj = ibb*(ibb-1)/2 + jj
                           if (iai.ne.ibj) then
                              val = yy(int)*yy(int)
c
                              difen = f(iaa) + g(ibb) - f(ii) - g(jj)
                              ff = 1.0d0/difen
                              val4 = val*ff*4 + val4
                           end if
                        end if
                     end if
                  end if
               end if
c
c
 120        continue
            go to 110
         end if
 130  continue
c
c
c
      e2 = -quart*(val1+val2+val3+val4)
      return
      end
      subroutine mpsrt6(q,iq,ncorb,iblki,iblki1,munit,ifort,iblz,
     *iblzz,ifilz,ixoxa,ixoxb)
       implicit REAL  (a-h,o-z)
c
c-------------------------------------------------------------
c     sort out coulomb integrals
c-------------------------------------------------------------
INCLUDE(common/gmempara)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n2,nbuck
INCLUDE(common/uhfspn)
INCLUDE(common/common)
INCLUDE(common/iofile)
      dimension q(*),iq(*)
      parameter (maxbuc=1500)
INCLUDE(common/blksiz)
      common/craypk/labs(1360)
INCLUDE(common/atmblk)
      logical swop1,swop2
c     logical uhf
c     character *8 open
c     data open/'open'/
c
c     ibl5 = no of integrals in block of sortfile
c
c     uhf = scftyp.eq.open
c
      character*6 fnm,snm
      data fnm,snm/'mp2b.m','mpsrt6'/
c
      call setsto(1360,0,labs)
      ibl5 = nsz340
      iilen = nsz340*mach12/2
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nij = ncorb*(ncorb+1)/2
      n2 = ncorb*ncorb
c
c       are going to sort the integrals from file ifili and
c       starting block iblki so that for ij (i.ge.j) all kl
c       integrals are available in a square on stream ifort
c
c       maxt is the number of triangles (squares for exchange ints)
c       which can be held in core (allowing n2 wkspace for reading back)
c       which is the number in each bucket
c
      i1  = igmem_alloc_all_inf(maxq,fnm,snm,'i1',IGMEM_DEBUG)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq-n2)/nij
      nword = (maxq/(1+lenrel(1)))*lenrel(1)
c
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
c
c     nbuck is the number of buckets required
c
      nbuck = (nij/maxt) + 1
      nadd = min(maxt,nij)
      maxa = nbuck*(ibl5+ibl5i) + n2
c
c
      if (nbuck.gt.maxb) then
         write (iwr,6010) maxq , maxa
         call caserr('stop')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      swop1 = .true.
      swop2 = .false.
      do 20 i = 1 , npassm
         call rdsrtj(q(i1),iq(ii1),mupblk(iblki),munit,swop1,swop2,i)
         if (i.lt.npassm) iblki = iblki + 2
 20   continue
c
      ipss = 2
      swop1 = .false.
      swop2 = .true.
      do 30 i = 1 , npassm
         call rdsrtj(q(i1),iq(ii1),mupblk(iblki1),munit,swop1,swop2,
     +               ipss)
         if (i.lt.npassm) iblki1 = iblki1 + 2
 30   continue
c
      if (ixoxa.eq.1) then
         swop1 = .true.
         swop2 = .false.
         call rdsrtj(q(i1),iq(ii1),iblz,ifilz,swop1,swop2,ipss)
      end if
      if (ixoxb.eq.1) then
         swop1 = .false.
         swop2 = .true.
         call rdsrtj(q(i1),iq(ii1),iblzz,ifilz,swop1,swop2,ipss)
      end if
c
c       read through the sort file to give final result
c
      maxqq = maxq - n2
      ipss = 1
      call wtsrtj(q(i1),q(i1+n2),maxqq,ncorb,ifort,ipss)
c
      call closbf(0)
      call gmem_free_inf(i1,fnm,snm,'i1')
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine ump2pd(c,iblw,ifw,ifytr,ifdma,ifdmb,iblkya,iblkyb,
     1 ifint1,ifint2,ifint3,ifint4,ifint5,ifint6,
     1 ifort1,ifort2)
      implicit REAL  (a-h,o-z)
INCLUDE(common/gmempara)
INCLUDE(common/sizes)
      dimension c(*)
c
      logical dpres,gpres
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/nshel)
INCLUDE(common/atmblk)
INCLUDE(common/symtry)      
      logical uhf
      character *8 open
      character *6 fnm,snm
      data fnm,snm/'mp2b.m','ump2pd'/
      data open /'open'/
      uhf = scftyp.eq.open
c
      dpres = .false.
      fpres = .false.
      gpres = .false.
      do 20 i = 1 , nshell
         if (ktype(i).eq.3) dpres = .true.
         if (ktype(i).eq.4) fpres = .true.
         if (ktype(i).eq.5) gpres = .true.
         if (ktype(i).gt.5) call caserr('ump2pd: upto g-functions only')
 20   continue
      nsq = ncoorb*ncoorb
      nova = nocca*nvirta
      novb = noccb*nvirtb
      ntri = ncoorb*(ncoorb+1)/2
      iso = 1
      i1 = iso + nw196(5)
      i2 = i1 + nsq
      i3 = i2 + nsq
      i4 = i3 + nsq
      i5 = i4 + ncoorb
      i6 = i5 + ncoorb
      i7 = i6 + nova
      itop = i7 + novb
      maxq = igmem_max_memory()
      if (itop.gt.maxq) then
         write (iwr,6010) maxq , itop
         call caserr(' not enough core ')
      end if
cjmht - save i4 in isave as required for calculating requirements further on.
c     We cannot use i4 directly as was the case previously, as with the new 
c     memory allocation routines, the offsets returned by igmem_alloc_inf
c     are set to null when the memory is freed.
c
      isave = i4
c
      iso = igmem_alloc_inf(nw196(5),fnm,snm,'iso',IGMEM_DEBUG)
      i1  = igmem_alloc_inf(nsq,fnm,snm,'i1',IGMEM_NORMAL)
      i2  = igmem_alloc_inf(nsq,fnm,snm,'i2',IGMEM_NORMAL)
      i3  = igmem_alloc_inf(nsq,fnm,snm,'i3',IGMEM_NORMAL)
      i4  = igmem_alloc_inf(ncoorb,fnm,snm,'i4',IGMEM_DEBUG)
      i5  = igmem_alloc_inf(ncoorb,fnm,snm,'i5',IGMEM_DEBUG)
      i6  = igmem_alloc_inf(nova,fnm,snm,'i6',IGMEM_DEBUG)
      i7  = igmem_alloc_inf(novb,fnm,snm,'i7',IGMEM_DEBUG)
c
      m9 = 9
      m12 = 12
      call secget(isect(9),m9,isec9)
      call rdedx(c(i4),ncoorb,isec9,ifild)
      call secget(isect(12),m12,isec12)
      call rdedx(c(i5),ncoorb,isec12,ifild)
c
      call umps1(c(i1),c(i2),c(i3),c(i4),c(i5),c(i6),c(i7),
     + nocca,ncoorb,nova,noccb,novb,ifint1,ifint2,ifint4,ifort1,ifort2)
c
      call delfil(ifint1)
      call delfil(ifint2)
      call delfil(ifint4)
      call gmem_free_inf(i7,fnm,snm,'i7')
      call gmem_free_inf(i6,fnm,snm,'i6')
      call gmem_free_inf(i5,fnm,snm,'i5')
      call gmem_free_inf(i4,fnm,snm,'i4')
c
c      i5 = i4 + nsq
      i5 = isave + nsq
      i6 = i5 + nova
      itop = i6 + novb

      if (itop.gt.maxq) then
         write (iwr,6020) maxq , itop
         call caserr(' not enough core ')
      end if

      i4  = igmem_alloc_inf(nsq,fnm,snm,'i4',IGMEM_NORMAL)
      i5  = igmem_alloc_inf(nova,fnm,snm,'i5',IGMEM_DEBUG)
      i6  = igmem_alloc_inf(novb,fnm,snm,'i6',IGMEM_DEBUG)
c
      m8 = 8
      m11 = 11
      isec = isect(8)
      itype = 0
      call secget(isec,itype,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i1),nsq,ibl,ifild)
      isec = isect(11)
      itype = 0
      call secget(isec,m11,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i2),nsq,ibl,ifild)
c
      call umps2(c(i1),c(i2),c(i3),c(i4),c(i5),c(i6),ifort1,ifort2,
     +  ifint1,nocca,nvirta,noccb,nvirtb,ncoorb)

c
      call gmem_free_inf(i6,fnm,snm,'i6')
      call gmem_free_inf(i5,fnm,snm,'i5')
      call gmem_free_inf(i4,fnm,snm,'i4')
c
c      i5 = i4 + ncoorb
      i5 = isave + ncoorb
      i6 = i5 + ncoorb
      i7 = i6 + novb
      itop = i7 + nova
      if (itop.gt.maxq) then
         write (iwr,6030) maxq , itop
         call caserr(' not enough core ')
      end if

      i4  = igmem_alloc_inf(ncoorb,fnm,snm,'i4',IGMEM_DEBUG)
      i5  = igmem_alloc_inf(ncoorb,fnm,snm,'i5',IGMEM_DEBUG)
      i6  = igmem_alloc_inf(novb,fnm,snm,'i6',IGMEM_DEBUG)
      i7  = igmem_alloc_inf(nova,fnm,snm,'i7',IGMEM_DEBUG)
      call secget(isect(12),m12,isec12)
      call rdedx(c(i4),ncoorb,isec12,ifild)
      call secget(isect(9),m9,isec9)
      call rdedx(c(i5),ncoorb,isec9,ifild)
c
      call umps1(c(i1),c(i2),c(i3),c(i4),c(i5),c(i6),c(i7),noccb,ncoorb
     +  ,novb,nocca,nova,ifint5,ifint6,ifint3,ifort1,ifort2)
c
      call delfil(ifint3)
      call delfil(ifint5)
      call delfil(ifint6)
c
      call gmem_free_inf(i7,fnm,snm,'i7')
      call gmem_free_inf(i6,fnm,snm,'i6')
      call gmem_free_inf(i5,fnm,snm,'i5')
      call gmem_free_inf(i4,fnm,snm,'i4')
c
c      i5 = i4 + nsq
      i5 = isave + nsq
      i6 = i5 + novb
      itop = i6 + nova
      if (itop.gt.maxq) then
         write (iwr,6040) maxq , itop
         call caserr(' not enough core ')
      end if
      isave = i6

      i4  = igmem_alloc_inf(nsq,fnm,snm,'i4',IGMEM_NORMAL)
      i5  = igmem_alloc_inf(novb,fnm,snm,'i5',IGMEM_DEBUG)
      i6  = igmem_alloc_inf(nova,fnm,snm,'i6',IGMEM_DEBUG)
c
      isec = isect(11)
      itype = 0
      call secget(isec,m11,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i1),nsq,ibl,ifild)
      isec = isect(8)
      itype = 0
      call secget(isec,m8,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i2),nsq,ibl,ifild)
c
      call umps2(c(i1),c(i2),c(i3),c(i4),c(i5),c(i6),ifort1,ifort2,
     +  ifint2,noccb,nvirtb,nocca,nvirta,ncoorb)
c
      call delfil(ifort1)
      call delfil(ifort2)
      call gmem_free_inf(i6,fnm,snm,'i6')
c     call gmem_free_inf(i5,fnm,snm,'i5')
c
c  transpose matrix
c
      ibase = igmem_alloc_all_inf(maxa,fnm,snm,'ibase',IGMEM_DEBUG)
      call mpsrt0(ifint1,ifort1,nova,nsq,c(ibase),c(ibase),maxa)
      call delfil(ifint1)
      call mpsrt0(ifint2,ifort2,novb,nsq,c(ibase),c(ibase),maxa)
      call delfil(ifint2)
      call gmem_free_inf(ibase,fnm,snm,'ibase')
c
      isec = isect(11)
      itype = 0
      call secget(isec,m11,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i1),nsq,ibl,ifild)
      isec = isect(8)
      itype = 0
      call secget(isec,m8,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i2),nsq,ibl,ifild)
c
      call umps3(c(i2),c(i1),c(i3),c(i4),ifytr,ifdma,ifdmb,nocca,noccb,
     +  ncoorb,iblkya,iblkyb,uhf)
c
c      i7 = i6 + nova
      i7 = isave + nova
      i8 = i7 + nsq
      i9 = i8 + nsq
      i10 = i9 + nsq
      i11 = i10 + nsq
      i12 = i11 + nsq
      itop = i12 + ncoorb*(ncoorb+1)/2*16
      if (dpres) itop = i12 + ncoorb*(ncoorb+1)/2*36
      if (fpres) itop = i12 + ncoorb*(ncoorb+1)/2*100
      if (gpres) itop = i12 + ncoorb*(ncoorb+1)/2*225
      if (itop.gt.maxq) then
         write (iwr,6050) maxq , itop
         call caserr(' not enough core ')
      end if
      i6  = igmem_alloc_inf(nova,fnm,snm,'i6',IGMEM_DEBUG)
      i7  = igmem_alloc_inf(nsq,fnm,snm,'i7',IGMEM_NORMAL)
      i8  = igmem_alloc_inf(nsq,fnm,snm,'i8',IGMEM_NORMAL)
      i9  = igmem_alloc_inf(nsq,fnm,snm,'i9',IGMEM_NORMAL)
      i10 = igmem_alloc_inf(nsq,fnm,snm,'i10',IGMEM_NORMAL)
      i11 = igmem_alloc_inf(nsq,fnm,snm,'i11',IGMEM_NORMAL)
      if (gpres) then
         i12 = igmem_alloc_inf(ncoorb*(ncoorb+1)/2*225,fnm,snm,
     &                         'i12',IGMEM_NORMAL)
      else if (fpres) then
         i12 = igmem_alloc_inf(ncoorb*(ncoorb+1)/2*100,fnm,snm,
     &                         'i12',IGMEM_NORMAL)
      else if (dpres) then
         i12 = igmem_alloc_inf(ncoorb*(ncoorb+1)/2*36,fnm,snm,
     &                         'i12',IGMEM_NORMAL)
      else
         i12 = igmem_alloc_inf(ncoorb*(ncoorb+1)/2*16,fnm,snm,
     &                         'i12',IGMEM_NORMAL)
      endif
c
      isec = isect(8)
      call secget(isec,m8,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i2),nsq,ibl,ifild)
      isec = isect(11)
      call secget(isec,m11,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i1),nsq,ibl,ifild)
c
      call rdedx(c(iso),nw196(5),ibl196(5),ifild)
      cut = 10.0d0*(-icut)
      call umps4(c(iso),nshell,c(i2),c(i1),c(i3),c(i4),c(i5),c(i6),
     +  c(i7),c(i8),c(i9),c(i10),c(i11),c(i12),
     +  nocca,nvirta,noccb,nvirtb,ncoorb,
     +  ntri,ifytr,ifdma,ifdmb,iblw,ifw,ifort1,ifort2,iblkya,iblkyb,
     +  cut)
c
      call delfil(ifort1)
      call delfil(ifort2)
c
      call gmem_free_inf(i12,fnm,snm,'i12')
      call gmem_free_inf(i11,fnm,snm,'i11')
      call gmem_free_inf(i10,fnm,snm,'i10')
      call gmem_free_inf(i9,fnm,snm,'i9')
      call gmem_free_inf(i8,fnm,snm,'i8')
      call gmem_free_inf(i7,fnm,snm,'i7')
      call gmem_free_inf(i6,fnm,snm,'i6')
      call gmem_free_inf(i5,fnm,snm,'i5')
      call gmem_free_inf(i4,fnm,snm,'i4')
      call gmem_free_inf(i3,fnm,snm,'i3')
      call gmem_free_inf(i2,fnm,snm,'i2')
      call gmem_free_inf(i1,fnm,snm,'i1')
      call gmem_free_inf(iso,fnm,snm,'iso')
      return
 6010 format (/' insufficient core for umps1 ',i9,' real words ',
     +        ' need ',i9,' real words ')
 6020 format (/' insufficient core for umps2 ',i9,' real words ',
     +        ' need ',i9,' real words ')
 6030 format (/' insufficient core for umps1b ',i9,' real words ',
     +        ' need ',i9,' real words ')
 6040 format (/' insufficient core for umps2b ',i9,' real words ',
     +        ' need ',i9,' real words ')
 6050 format (/' insufficient core for umps4 ',i9,' real words ',
     +        ' need ',i9,' real words ')
      end
      subroutine umps1(t1,t2,t3,e1,e2,b1,b2,nocc,ncoorb,nov,
     1 noccz,novz,ifint1,ifint2,ifint3,ifort1,ifort2)
      implicit REAL  (a-h,o-z)
      dimension e1(ncoorb),e2(ncoorb),t1(ncoorb,ncoorb),
     1 t2(ncoorb,ncoorb),t3(ncoorb,ncoorb),b1(nov),b2(novz)
c
      nsq = ncoorb*ncoorb
      nocc1 = nocc + 1
      noccz1 = noccz + 1
      call rewedz(ifint1)
      call rewedz(ifint2)
      call rewedz(ifint3)
      call rewedz(ifort1)
      call rewedz(ifort2)
c
      do 70 ip = 1 , ncoorb
         do 60 iq = 1 , ip
            call rdedz(t1,nsq,ifint1)
            call rdedz(t2,nsq,ifint2)
            call rdedz(t3,nsq,ifint3)
            if (ip.gt.nocc) then
               if (iq.le.nocc) then
                  ebj = e1(ip) - e1(iq)
c
                  do 30 ia = nocc1 , ncoorb
                     do 20 i = 1 , nocc
                        iai = (ia-nocc-1)*nocc + i
                        b1(iai) = (t1(i,ia)-t2(i,ia))/(ebj+e1(ia)-e1(i))
 20                  continue
 30               continue
                  call wtedz(b1,nov,ifort1)
c
                  do 50 ia = noccz1 , ncoorb
                     do 40 i = 1 , noccz
                        iai = (ia-noccz-1)*noccz + i
                        b2(iai) = t3(i,ia)/(ebj+e2(ia)-e2(i))
 40                  continue
 50               continue
                  call wtedz(b2,novz,ifort2)
               end if
            end if
 60      continue
 70   continue
      return
      end
      subroutine umps2(vec,vecz,dum,r,b1,b2,ifort1,ifort2,ifort3,
     1 nocc,nvirt,noccz,nvirtz,ncoorb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
      dimension b1(nocc*nvirt),b2(noccz*nvirtz),dum(ncoorb*ncoorb),
     1 vec(ncoorb*ncoorb),vecz(ncoorb*ncoorb),r(ncoorb*ncoorb)
c
      nsq = ncoorb*ncoorb
      nov = nocc*nvirt
      novz = noccz*nvirtz
      call rewedz(ifort1)
      call rewedz(ifort2)
      call rewedz(ifort3)
      do 80 ijkl = 1 , nov
         call rdedz(b1,nov,ifort1)
         call rdedz(b2,novz,ifort2)
         iv1 = ncoorb*nocc + 1
         iv2 = ncoorb*noccz + 1
         call vclr(dum,1,nsq)
c
         call mxmb(b1,1,nocc,vec(iv1),ncoorb,1,dum,1,nocc,nocc,nvirt,
     +             ncoorb)
         call vclr(r,1,nsq)
         call mxmb(vec,1,ncoorb,dum,1,nocc,r,1,ncoorb,ncoorb,nocc,
     +             ncoorb)
         call vclr(dum,1,nsq)
         call mxmb(b2,1,noccz,vecz(iv2),ncoorb,1,dum,1,noccz,noccz,
     +             nvirtz,ncoorb)
         call mxmb(vecz,1,ncoorb,dum,1,noccz,r,1,ncoorb,ncoorb,noccz,
     +             ncoorb)
c
c  symmetrise
c
         do 30 im = 1 , ncoorb
            do 20 in = 1 , ncoorb
               imn = (in-1)*ncoorb + im
               inm = (im-1)*ncoorb + in
               dum(imn) = (r(imn)+r(inm))/2.0d0
 20         continue
 30      continue
c
         icount = 1
         do 70 ii = 1 , nshell
            do 60 jj = 1 , ii
               mini = kmin(ii)
               minj = kmin(jj)
               maxi = kmax(ii)
               maxj = kmax(jj)
               loci = kloc(ii) - mini
               locj = kloc(jj) - minj
c
               do 50 i = mini , maxi
                  i1 = loci + i
                  do 40 j = minj , maxj
                     j1 = locj + j
c
                     ij1 = (j1-1)*ncoorb + i1
                     r(icount) = dum(ij1)
                     icount = icount + 1
 40               continue
 50            continue
 60         continue
 70      continue
         call wtedz(r,nsq,ifort3)
 80   continue
      return
      end
      subroutine umps4(iso,nshels,
     1 veca,vecb,cssa,cssb,d2,d1,dum,ys,b,ya,yb,bs,
     1 nocca,nvirta,noccb,nvirtb,ncoorb,ntri,ifytr,ifdma,ifdmb,
     1 iblw,ifw,ifort1,ifort2,iblkya,iblkyb,cut)
      implicit REAL  (a-h,o-z)
      dimension iso(nshels,*)
      dimension veca(ncoorb*ncoorb),vecb(ncoorb*ncoorb),
     1 cssa(ncoorb*ncoorb),cssb(ncoorb*ncoorb),dum(ncoorb*ncoorb),
     1 ys(ncoorb*ncoorb),b(ncoorb,ncoorb),d1(nocca*nvirta),
     1 d2(noccb*nvirtb),ya(ncoorb*ncoorb),yb(ncoorb*ncoorb),bs(2)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
INCLUDE(common/symtry)
INCLUDE(common/atmblk)
_IFN1(iv)      common/craypk/labout(1360)
_IF1(iv)      common/craypk/labij(340),labkl(340)
      common/blkin/g(510),nint,nxtr
      dimension m0(48)
      logical lab,labc,labcd
      logical ijump,jjump
      data mzero/0/
c
      tempto = 0.0d0
      nsq = ncoorb*ncoorb
      nova = nocca*nvirta
      novb = noccb*nvirtb
      ifile1 = 1
      call rdedx(ys,nsq,ifytr,ifile1)
      call rdedx(ya,nsq,iblkya,ifile1)
      call rdedx(yb,nsq,iblkyb,ifile1)
      call rdedx(cssa,nsq,ifdma,ifile1)
      call rdedx(cssb,nsq,ifdmb,ifile1)
c
      call rewedz(ifort1)
      call rewedz(ifort2)
      call search(iblw,ifw)
      nint = 0
      do 160 ii = 1 , nshell
       iceni = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ijump = id.gt.ii
            m0(it) = id
 30      continue
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
               call rdedz(d1,nova,ifort1)
               call rdedz(d2,novb,ifort2)
               if (.not.(ijump .or. jjump)) then
                  iv1 = ncoorb*nocca + 1
                  iv2 = ncoorb*noccb + 1
                  call vclr(b,1,nsq)
                  call vclr(dum,1,nsq)
c
                  call mxmb(veca(iv1),1,ncoorb,d1,nocca,1,dum,1,ncoorb,
     +                      imax,nvirta,nocca)
                  call mxmb(veca,1,ncoorb,dum,ncoorb,1,b,1,ncoorb,imax,
     +                      nocca,imax)
c
                  call vclr(dum,1,nsq)
                  call mxmb(vecb(iv2),1,ncoorb,d2,noccb,1,dum,1,ncoorb,
     +                      imax,nvirtb,noccb)
                  call mxmb(vecb,1,ncoorb,dum,ncoorb,1,b,1,ncoorb,imax,
     +                      noccb,imax)
c
                  do 70 ms1 = 1 , imax
                     do 60 ms2 = 1 , ms1
                        ms12 = iky(ms1) + ms2
                        bs(ib1+ms12) = (b(ms1,ms2)+b(ms2,ms1))*0.5d0
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
                                    val = -bs(ipos+ikl)
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
                                       g(nint) = val
                                       tempto = tempto + dabs(val)
_IF(ibm,vax)
                                        labij(nint) = j1 + i4096(i1)
                                        labkl(nint) = l1 + i4096(k1)
                                       if (nint.eq.num2e) then
                                         call pak4v(labij,g(num2e+1))
_ELSE
                                       nint4 = nint + nint + nint + nint
                                       labout(nint4-3) = i1
                                       labout(nint4-2) = j1
                                       labout(nint4-1) = k1
                                       labout(nint4) = l1
                                       if (nint.eq.num2e) then
                                         call pack(g(num2e+1),lab816,
     +                                             labout,numlab)
_ENDIF
                                         call put(g,m511,ifw)
                                         nint = 0
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
      subroutine umpmkw(e,y1,y2,w,a1,a2,nocc,ncoorb,
     1 istrma,istrmb,ibly1,ibly2,iblw)
      implicit REAL  (a-h,o-z)
      dimension y1(ncoorb,ncoorb),y2(ncoorb,ncoorb),
     1 a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),w(ncoorb,ncoorb),
     1 e(ncoorb)
c
      nsq = ncoorb*ncoorb
      nocc1 = nocc + 1
      ifile1 = 1
      call rewedz(istrma)
      call rewedz(istrmb)
      call rdedx(y1,nsq,ibly1,ifile1)
      call rdedx(y2,nsq,ibly2,ifile1)
      call rdedx(a1,nsq,iblw,ifile1)
c
      do 30 ir = 1 , nocc
         do 20 is = 1 , nocc
c     w(ir,is)=y1(ir,is)*(e(ir)+e(is))*0.5
c   minus sign as y has minus occ-occ and virt-occ elements
            w(ir,is) = -y1(ir,is)*e(is)
            w(ir,is) = w(ir,is) + a1(ir,is)
 20      continue
 30   continue
c
      do 50 ir = nocc1 , ncoorb
         do 40 is = nocc1 , ncoorb
c     w(ir,is)=-y1(ir,is)*(e(ir)+e(is))*0.5
            w(ir,is) = -y1(ir,is)*e(is)
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
      subroutine umpmky(y,zl,t1,t2,t3,t,tab,e1,e2,nocc,nvirt,ncoorb,
     1 noccz,ifint1,ifint2,ifint3,istrma,
     1 ibly,iblw)
      implicit REAL  (a-h,o-z)
      dimension y(ncoorb,ncoorb),zl(ncoorb,ncoorb),t1(ncoorb,ncoorb),
     1 t3(ncoorb,ncoorb),t(ncoorb,ncoorb),tab(ncoorb,ncoorb),
     1 e1(ncoorb),e2(ncoorb),t2(ncoorb,ncoorb)
c
      tempx = 0.0d0
      tempt = 0.0d0
      tempy = 0.0d0
      tempzl = 0.0d0
      tempw = 0.0d0
      tempv = 0.0d0
      ifile1 = 1
      nocc1 = nocc + 1
      nsq = ncoorb*ncoorb
      call vclr(y,1,nsq)
      call vclr(zl,1,nsq)
      call rewedz(istrma)
      call rewedz(ifint1)
      call rewedz(ifint2)
      call rewedz(ifint3)
c
c  read integrals
c
      do 150 ip = 1 , ncoorb
         do 140 iq = 1 , ip
            call rdedz(t1,nsq,ifint1)
            call rdedz(t2,nsq,ifint2)
            call rdedz(t3,nsq,ifint3)
            do 30 jki = 1 , ncoorb
               do 20 jkl = 1 , ncoorb
                  tempw = tempw + dabs(t2(jki,jkl))
                  tempv = tempv + dabs(t1(jki,jkl))
                  tempx = tempx + dabs(t3(jki,jkl))
 20            continue
 30         continue
c
c  construct a matrix over all orbitals
c
            do 50 i = 1 , ncoorb
               do 40 j = 1 , ncoorb
                  t(i,j) = t1(i,j) + t1(i,j) - t2(i,j) - t2(j,i)
                  tab(i,j) = t3(i,j) + t3(i,j)
 40            continue
 50         continue
            call wtedz(t,nsq,istrma)
            do 70 ki = 1 , ncoorb
               do 60 kj = 1 , ncoorb
                  tempt = tempt + dabs(t(ki,kj))
 60            continue
 70         continue
c
            call vclr(t,1,nsq)
            call vclr(tab,1,nsq)
c
            if (iq.le.nocc) then
               if (ip.gt.nocc) then
                  ebj1 = e1(ip) - e1(iq)
c
c  construct t matrix
c
                  do 90 i = 1 , nocc
                     do 80 ia = nocc1 , ncoorb
                        t(i,ia) = (t1(i,ia)-t2(i,ia))
     +                            /(ebj1+e1(ia)-e1(i))
 80                  continue
 90               continue
c
c  construct part of zl matrix
c
                  call mxmb(t(1,nocc1),1,ncoorb,t1(nocc1,1),1,ncoorb,zl,
     +                      ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t1,1,ncoorb,t(1,nocc1),1,ncoorb,zl(1,nocc1),
     +                      1,ncoorb,ncoorb,nocc,nvirt)
c
c  construct y matrix
c
                  do 110 i = 1 , nocc
                     do 100 ia = nocc1 , ncoorb
                        t1(i,ia) = t1(i,ia)/(ebj1+e1(ia)-e1(i))
 100                 continue
 110              continue
c
                  call mxmb(t(1,nocc1),ncoorb,1,t1(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(t(1,nocc1),1,ncoorb,t1(1,nocc1),ncoorb,1,y,
     +                      1,ncoorb,nocc,nvirt,nocc)
               end if
            end if
c
c
            if (iq.le.noccz) then
               if (ip.gt.noccz) then
                  ebj2 = e2(ip) - e2(iq)
                  do 130 i = 1 , nocc
                     do 120 ia = nocc1 , ncoorb
                        tab(i,ia) = t3(i,ia)/(ebj2+e1(ia)-e1(i))
 120                 continue
 130              continue
c
c   part of zl matrix
c
                  call mxmb(tab(1,nocc1),1,ncoorb,t3(nocc1,1),1,ncoorb,
     +                      zl,ncoorb,1,nocc,nvirt,ncoorb)
c
                  call mxmb(t3,1,ncoorb,tab(1,nocc1),1,ncoorb,
     +                      zl(1,nocc1),1,ncoorb,ncoorb,nocc,nvirt)
c
c   y  matrix
c
                  call mxmb(tab(1,nocc1),ncoorb,1,tab(1,nocc1),1,ncoorb,
     +                      y(nocc1,nocc1),1,ncoorb,nvirt,nocc,nvirt)
c
                  call mxmb(tab(1,nocc1),1,ncoorb,tab(1,nocc1),ncoorb,1,
     +                      y,1,ncoorb,nocc,nvirt,nocc)
               end if
            end if
c
 140     continue
 150  continue
      do 170 i = 1 , nocc
         do 160 j = 1 , nocc
            y(i,j) = -y(i,j)
 160     continue
 170  continue
      call wrt3(y,nsq,ibly,ifile1)
      call wrt3(zl,nsq,iblw,ifile1)
      do 190 ki = 1 , ncoorb
         do 180 kj = 1 , ncoorb
            tempy = tempy + dabs(y(ki,kj))
            tempzl = tempzl + dabs(zl(ki,kj))
 180     continue
 190  continue
      return
      end
      subroutine umpsq(veca,vecb,s,nocca,noccb,num,ntri,ssq)
      implicit REAL  (a-h,o-z)
      dimension veca(num,num),vecb(num,num),s(ntri)
c
      sz = 0.5d0*(nocca-noccb)
      ssq = sz*sz + sz + noccb
      do 50 i = 1 , nocca
         do 40 j = 1 , noccb
            cs = 0
            do 30 m = 1 , num
               do 20 n = 1 , m
                  imn = m*(m-1)/2 + n
                  cs = cs + veca(m,i)*vecb(n,j)*s(imn)
                  if (n.lt.m) then
                     cs = cs + veca(n,i)*vecb(m,j)*s(imn)
                  end if
 20            continue
 30         continue
            css = cs*cs
            ssq = ssq - css
 40      continue
 50   continue
      return
      end
      subroutine uhfmp2(q,iq)
      implicit REAL  (a-h,o-z)
c
c     uhf mp2 routines -
c
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/infoa)
INCLUDE(common/atmblk)
INCLUDE(common/prnprn)
INCLUDE(common/cigrad)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
INCLUDE(common/cndx40)
c
      dimension q(*),iq(*)
      common/scfblk/enucl,etot,ehf
INCLUDE(common/nshel)
INCLUDE(common/uhfspn)
      character*6 fnm
      character*6 snm
c
c     logical swop1,swop2
c     character *8 energy
c
      data m1/1/
c     data m2,m3,m4,l511/2,3,4,511/
c     data energy/'hfscf'/
      data zzero/0.0d0/
      data m13/13/
      data fnm,snm/"mp2b.m","ufhmp2"/
c
c
      if (npassm.gt.1) then
         call caserr(' enlarge the core for this job ')
      end if
      ifort1 = 17
      ifort2 = 18
      ifort3 = 19
      ifort4 = 20
      ifort5 = 21
      ifort6 = 22
      ifort7 = 23
      ifort8 = 24
c
c   to calculate pump2 energies, use ed4
c
      ifort9 = 5
      if (mprest.gt.2) then
c
c   sort the integral lists with mixed spins
c
         ifilz = 1
         iblzz = 1
         iblki = 2
         iblki1 = npassm*2 + 1
         call mpsrt6(q,iq,ncoorb,iblki,iblki1,junitf,ifort5,iblkz,iblzz,
     +  ifilz,ixoxxb,ixoxxa)
         iblki = 2
         iblki1 = npassm*2 + 1
         call mpsrt6(q,iq,ncoorb,iblki1,iblki,junitf,ifort6,iblzz,iblkz,
     +  ifilz,ixoxxa,ixoxxb)
         call delfil(1)
c
         call umpgrd(q,ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,ifort7
     +  ,ifort8)
c
         return
      else
c
c     read in all useful information
c
         call secget(isect(13),m13,isec13)
         call rdedx(enucl,lds(isect(13)),isec13,ifild)
c
c        e0 = ehf
         e1 = ehf
         write (iwr,6010) e1 , enucl , enucl + e1
c
c   sort two lists of integrals into coulomb and exchange terms
c
c     first block contains  alp alp | alp alp
c     second block contains bet bet | alp alp
c     third block contains  alp alp | bet bet
c     fourth block contains bet bet | bet bet
c
c    the above covers all integrals for restrict 1 and restrict 0
c    but for restrict 4  then ed0 contains some mixed spin integrals
c    of the type ( x x | x o ) ........
c
         if (.not.pump2) then
            iblki = 1
c           swop1 = .true.
c           swop2 = .true.
            call mpsrtj(q,iq,ncoorb,mupblk(iblki),junitf,ifort1)
            call mpsrtk(q,iq,ncoorb,mupblk(iblki),junitf,ifort2,m1)
c
            iblki = npassm*2 + 2
            call mpsrtj(q,iq,ncoorb,mupblk(iblki),junitf,ifort3)
            call mpsrtk(q,iq,ncoorb,mupblk(iblki),junitf,ifort4,m1)
         else
c
c    sort to form (ai|jb)-(aj|bi)  for alp alp alp alp and
c      bet bet bet bet
c
            iblki = 1
            nocasq = nocca*nocca
            nocbsq = noccb*noccb
            ntribv = nvirtb*(nvirtb+1)/2
            ntriav = nvirta*(nvirta+1)/2
            call mpsrt1(q,iq,mupblk(iblki),junitf,ifort5,nocasq,ntriav,
     +                  nocca)
            iblki = npassm*2 + 2
            call mpsrt1(q,iq,mupblk(iblki),junitf,ifort8,nocbsq,ntribv,
     +                  noccb)
            iblki = 2
            iblki1 = 3
            nocab = nocca*noccb
            nvab = nvirta*nvirtb
c
c     sort to form  (alp bet alp bet)  term
c
            call mpsrt2(q,iq,mupblk(iblki),mupblk(iblki1),junitf,ifort1,
     +                  nocab,nvab,nocca,noccb,nvirta)
         end if
c
         iq5 = 1
         iq6 = iq5 + num
         iq1 = iq6 + num
         iq2 = iq1 + num*num
         iq3 = iq2 + num*num
         iq4 = iq3 + num*num
         itop = iq4 + num*num
         maxq = igmem_max_memory()
         if (itop.gt.maxq) then
            write (iwr,6020) maxq , itop
            call caserr(' not enough core ')
         end if
c
         iq5 = igmem_alloc_inf(num,fnm,snm,'f',IGMEM_DEBUG)
         iq6 = igmem_alloc_inf(num,fnm,snm,'g',IGMEM_DEBUG)
         iq1 = igmem_alloc_inf(num*num,fnm,snm,'t1',IGMEM_NORMAL)
         iq2 = igmem_alloc_inf(num*num,fnm,snm,'t2',IGMEM_NORMAL)
         iq3 = igmem_alloc_inf(num*num,fnm,snm,'t3',IGMEM_NORMAL)
         iq4 = igmem_alloc_inf(num*num,fnm,snm,'t4',IGMEM_NORMAL)
c
         m9 = 9
         m12 = 12
         call secget(isect(9),m9,isec9)
         call rdedx(q(iq5),lds(isect(9)),isec9,ifild)
         call secget(isect(12),m12,isec12)
         call rdedx(q(iq6),lds(isect(12)),isec12,ifild)
c
c
         e2 = zzero
         if (.not.pump2) then
          call umpe2(q(iq1),q(iq2),q(iq3),q(iq4),q(iq5),q(iq6),val3,
     +               val4,e2,ifort1,ifort2,ifort3,ifort4,nocca,noccb,
     +               ncoorb,junitf)
            etot = enucl + ehf + e2
            m13 = 13
            length = lensec(lds(isect(13)))
            call secput(isect(13),m13,length,isec13)
            call wrt3(enucl,lds(isect(13)),isec13,ifild)
            if(odebug(38)) write (iwr,6030) val3, val4
            write (iwr,6040) e2 , etot
         end if
c
c
         ntri = num*(num+1)/2
         nsq = num*num
         m5 = 5
         m8 = 8
         m11 = 11
         isec = isect(8)
         call secget(isec,m8,ibl)
         ibl = ibl + mvadd
         call vclr(q(iq1),1,nsq)
         call rdedx(q(iq1),num*ncoorb,ibl,ifild)
         isec = isect(11)
         call secget(isec,m11,ibl)
         ibl = ibl + mvadd
         call vclr(q(iq2),1,nsq)
         call rdedx(q(iq2),num*ncoorb,ibl,ifild)
         isec = isect(5)
         call secget(isec,m5,ibl)
         call rdedx(q(iq3),ntri,ibl,ifild)
c
c  s squared expectation value
c
         if (.not.pump2) then
            call umpsq(q(iq1),q(iq2),q(iq3),nocca,noccb,
     +                 num,ntri,ssq)
            call gmem_free_inf(iq4,fnm,snm,'t4')
            call gmem_free_inf(iq3,fnm,snm,'t3')
            call gmem_free_inf(iq2,fnm,snm,'t2')
            call gmem_free_inf(iq1,fnm,snm,'t1')
         else
            call squr(q(iq3),q(iq4),ncoorb)
            call vclr(q(iq3),1,nsq)
            call mxmb(q(iq1),ncoorb,1,q(iq4),1,ncoorb,q(iq3),1,ncoorb,
     +                ncoorb,ncoorb,ncoorb)
c
            call vclr(q(iq1),1,nsq)
            call mxmb(q(iq3),1,ncoorb,q(iq2),1,ncoorb,q(iq1),1,ncoorb,
     +                ncoorb,ncoorb,ncoorb)
c
            call gmem_free_inf(iq4,fnm,snm,'t4')
            call gmem_free_inf(iq3,fnm,snm,'t3')
            call gmem_free_inf(iq2,fnm,snm,'t2')
c
            spin = 0.5d0*(nocca-noccb)
            iq23 = igmem_alloc_all_inf(maxq,fnm,snm,'iq23',IGMEM_DEBUG)
            iq24 = iq23 + nocca*nvirtb
            iq25 = iq24 + nocca*nvirtb
            iq26 = iq25 + nocca*nvirtb
            iq27 = iq26 + nocca*nvirtb
            iq28 = iq27 + nocca*nvirtb
            iq29 = iq28 + nocca*nvirtb
            iq30 = iq29 + nocca*nvirtb
            iq31 = iq30 + nocca*nvirtb
            iq32 = iq31 + nocca*nvirtb
            iq33 = iq32 + nocca*nvirtb
            iq34 = iq33 + nocca*nvirtb
            iq35 = iq34 + nocca*nvirtb
            iq36 = iq35 + ncoorb*ncoorb
            iq37 = iq36 + ncoorb*ncoorb
            iq38 = iq37 + nocca*nvirtb
            iq2 = iq38 + nocca*nvirtb
            nmx = max(nocca,nvirtb)
            iq3 = iq2 + nmx*nmx
            iq4 = iq3 + nmx*nmx
            iq7 = iq4 + nmx*nmx
            iq8 = iq7 + nmx*nmx
            iq9 = iq8 + nmx*nmx
            iq10 = iq9 + nmx*nmx
            iq11 = iq10 + nmx*nmx
            iq12 = iq11 + nmx*nmx
            iq13 = iq12 + nmx*nmx
            iq14 = iq13 + nmx*nmx
            iq15 = iq14 + nmx*nmx
            iq16 = iq15 + nmx*nmx
            iq17 = iq16 + nmx*nmx
            iq18 = iq17 + ncoorb*ncoorb
            iq19 = iq18 + ncoorb*ncoorb
            iq20 = iq19 + ncoorb*nmx
            iq21 = iq20 + ncoorb*nmx
            iq22 = iq21 + ncoorb*nmx
            iq39 = iq22 + ncoorb*nmx
            itop = iq39 + nvirtb*(nvirtb+1)/2
            if(odebug(39))write (iwr,6050) itop
            if (itop.gt.maxq) then
               write (iwr,*) ' insufficient core : require ' , itop ,
     +                     ' words ' , ' have ' , maxq , ' words '
               call caserr(' insufficient core ')
            end if
c
            maxql = maxq - iq2
c the next two lines were there before. I can't see any sense in them
c hope that the third line is an accurate replacement.
chvd        iofst = ncoorb*ncoorb + 2*ncoorb + 1
chvd        call vclr(q(iofst),1,maxq-iofst)
            call vclr(q(iq23),1,maxq-iq23)
            euhf = ehf + enucl
            call umpprj(
     +              q(iq5),q(iq6),q(iq1),euhf,nocca,noccb,ncoorb,spin,
     +                q(iq2),q(iq3),q(iq4),q(iq7),q(iq8),q(iq9),q(iq10),
     +                q(iq11),q(iq12),q(iq13),q(iq14),q(iq15),q(iq16),
     +                q(iq17),q(iq18),q(iq19),q(iq20),q(iq21),q(iq22),
     +                q(iq23),q(iq24),q(iq25),q(iq26),q(iq27),q(iq28),
     +                q(iq29),q(iq30),q(iq31),q(iq32),q(iq33),q(iq34),
     +                q(iq35),q(iq36),q(iq37),q(iq38),nvirta,nvirtb,
     +                ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,ifort7,
     +                ifort8,ifort9,q(iq2),q(iq2),maxql,q(iq39),emp2u)
            call gmem_free_inf(iq23,fnm,snm,'iq23')
            call gmem_free_inf(iq1,fnm,snm,'t1')
            call secget(isect(13),13,iblk13)
            call rdedx(enucl,lds(isect(13)),iblk13,ifild)
            etot = emp2u
            call wrt3(enucl,lds(isect(13)),iblk13,ifild)
c
         end if
c
         call gmem_free_inf(iq6,fnm,snm,'g')
         call gmem_free_inf(iq5,fnm,snm,'f')
c
      end if
 6010 format (/10x,47('*')/
     +        10x,'uhf-mp2 calculation'/10x,47('*')/
     +        10x,'first order electronic energy   ',f15.8/
     +        10x,'nuclear repulsion               ',f15.8/
     +        10x,'first order energy (total)      ',f15.8)
 6030 format (10x,47('*')/
     +        10x,'second order contributions      '/10x,47('*')/
     +        10x,'Val3                            ',f15.8/
     +        10x,'Val4                            ',f15.8)
 6040 format (
     +        10x,'second order perturbation energy',f15.8/10x,47('*')/
     +        10x,'total energy (uhf-mp2)          ',f15.8/
     +        10x,47('*'))
 6020 format (/' insufficient core for e2 ',i9,' real words ',' need ',
     +        i9,' real words ')
 6050  format(/1x,'using ', i6,' words for umpprj (excluding sorts)')
      end
      subroutine umpgrd(q,ifort1,ifort2,ifort3,ifort4,ifort5,ifort6,
     1   ifort7,ifort8)
      implicit REAL  (a-h,o-z)
INCLUDE(common/gmempara)
INCLUDE(common/sizes)
      dimension q(*)
INCLUDE(common/prnprn)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
c
INCLUDE(common/cigrad)
INCLUDE(common/common)
INCLUDE(common/iofile)
INCLUDE(common/cndx41)
c
INCLUDE(common/atmblk)
      logical uhf
      character *8 open
      character *6 fnm,snm
      data fnm,snm/'mp2b.m','umpgrd'/
      data open /'open'/
      uhf = scftyp.eq.open
c     if(mprest.gt.4)goto 109
      mna = nocca*nvirta
      mnb = noccb*nvirtb
      mn = mna + mnb
c     nat3 = nat*3
      ntri = ncoorb*(ncoorb+1)/2
      nsq = ncoorb*ncoorb
      ifile1 = 1
c
c     naba = nvirta*(nvirta+1)/2
c     nabb = nvirtb*(nvirtb+1)/2
c     nija = nocca*(nocca+1)/2
c     nijb = noccb*(noccb+1)/2
c     ndepa = naba + nija
c     ndepb = nabb + nijb
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
c
      call wrt3z(1,ifile1,mpblk(11))
      iyuk = 1 + mn
      i1 = iyuk + ncoorb
      i2 = i1 + ncoorb
      i3 = i2 + nsq
      i4 = i3 + nsq
      i5 = i4 + nsq
      i6 = i5 + nsq
      i7 = i6 + nsq
      i8 = i7 + nsq
      itop = i8 + nsq
      maxq = igmem_max_memory()
      if (itop.gt.maxq) then
         write (iwr,6010) maxq , itop
         call caserr(' not enough core ')
      end if
c
      iyuk = igmem_alloc_inf(ncoorb,fnm,snm,'iyuk',IGMEM_DEBUG)
      i0   = igmem_alloc_inf(mn,fnm,snm,'i0',IGMEM_DEBUG)
      i1   = igmem_alloc_inf(ncoorb,fnm,snm,'i1',IGMEM_DEBUG)
      i2   = igmem_alloc_inf(nsq,fnm,snm,'i2',IGMEM_NORMAL)
      i3   = igmem_alloc_inf(nsq,fnm,snm,'i3',IGMEM_NORMAL)
      i4   = igmem_alloc_inf(nsq,fnm,snm,'i4',IGMEM_NORMAL)
      i5   = igmem_alloc_inf(nsq,fnm,snm,'i5',IGMEM_NORMAL)
      i6   = igmem_alloc_inf(nsq,fnm,snm,'i6',IGMEM_NORMAL)
      i7   = igmem_alloc_inf(nsq,fnm,snm,'i7',IGMEM_NORMAL)
      i8   = igmem_alloc_inf(nsq,fnm,snm,'i8',IGMEM_NORMAL)
c
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(q(iyuk),ncoorb,isec9,ifild)
      m12 = 12
      call secget(isect(12),m12,isec12)
      call rdedx(q(i1),ncoorb,isec12,ifild)
c
      call umpmky(q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),q(i8),q(iyuk),
     +  q(i1),nocca,nvirta,ncoorb,noccb,ifort1,ifort2,ifort5,ifort7,
     +  mpblk(1),mpblk(2))
      call umpmky(q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),q(i8),q(i1),
     +  q(iyuk),noccb,nvirtb,ncoorb,nocca,ifort3,ifort4,ifort6,ifort8,
     +  mpblk(3),mpblk(4))
      call umpmkz(q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),nocca,ncoorb,mna,
     +  ifort7,ifort6,mpblk(1),mpblk(3),mpblk(2))
      call umpmkz(q(i2),q(i3),q(i4),q(i5),q(i6),q(i7+mna),noccb,ncoorb,
     +  mnb,ifort8,ifort5,mpblk(3),mpblk(1),mpblk(4))
      call wrt3(q(i7),mn,iblks,ifils)
c
c  solve cpuhf equations to get z matrix
c
      np = 1
      npstar = 0
      jblk(1) = 1
      call cuhf2(q(i0),q(iyuk),q(i1),q(i2),q(i3),q(i4),q(i5),
     +  nocca,noccb,
     +  ncoorb,mna,mn,ifort7,ifort6,ifort5,ifort8,jblk(1),nofile(1))
c
      call delfil(nofile(1))
      iblkz = iblks + lensec(mn)
      call rdedx(q(i4),mn,iblkz,ifils)
      call umpzy(q(i2),q(i3),q(i4),nocca,noccb,ncoorb,mn,mpblk(1),
     +  mpblk(3))
c
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(q(iyuk),ncoorb,isec9,ifild)
      m12 = 12
      call secget(isect(12),m12,isec12)
      call rdedx(q(i1),ncoorb,isec12,ifild)
c
      call umpmkw(q(iyuk),q(i2),q(i3),q(i4),q(i5),q(i6),nocca,ncoorb,
     +  ifort7,ifort5,mpblk(1),mpblk(3),mpblk(2))
c
      call umpmkw(q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),noccb,ncoorb,
     +  ifort8,ifort6,mpblk(3),mpblk(1),mpblk(4))
c
c ya, yb, wa, wb stored on ed0 in mo basis
c
c
      ityp = 0
      len = lensec(ntri)
      m8 = 8
      call secget(isect(8),m8,iblok)
      iblok = iblok + mvadd
      call rdedx(q(i2),num*ncoorb,iblok,ifild)
      m11 = 11
      call secget(isect(11),m11,iblok)
      iblok = iblok + mvadd
      call rdedx(q(i6),num*ncoorb,iblok,ifild)
      call vclr(q(i5),1,nsq)
      call rdedx(q(i3),nsq,mpblk(1),ifile1)
      call vtamvu(q(i3),q(i2),q(i4),q(i5),ncoorb)
c
c     iblkya = iblkz + lensec(mn)
c     iblkyb = iblkya + lensec(nsq)
      call wrt3(q(i5),nsq,mpblk(5),ifile1)
      call vclr(q(i7),1,nsq)
c
      call rdedx(q(i3),nsq,mpblk(3),ifile1)
      call vtamvu(q(i3),q(i6),q(i4),q(i7),ncoorb)
      call wrt3(q(i7),nsq,mpblk(6),ifile1)
c
      ij = 0
      do 30 i = 1 , ncoorb
         do 20 j = 1 , ncoorb
            q(i5+ij) = q(i5+ij) + q(i7+ij)
            ij = ij + 1
 20      continue
 30   continue
      call wrt3(q(i5),nsq,mpblk(7),ifile1)
      call trsqsq(q(i5),q(i4),ncoorb)
      call secput(isecdd,ityp,len,isdd)
      call wrt3(q(i4),ntri,isdd,ifild)
      if (odebug(38)) then
       write(iwr,*) ' ump2 density matrix contribution '
        call prtris(q(i4),num,iwr)
      endif
c
      call rdedx(q(i3),nsq,mpblk(2),ifile1)
      call vclr(q(i5),1,nsq)
      call vtamvu(q(i3),q(i2),q(i4),q(i5),ncoorb)
      call rdedx(q(i3),nsq,mpblk(4),ifile1)
      call vtamvu(q(i3),q(i6),q(i4),q(i5),ncoorb)
      call wrt3(q(i5),nsq,mpblk(8),ifile1)
      call trsqsq(q(i5),q(i4),ncoorb)
      call secput(isecll,ityp,len,isll)
      call wrt3(q(i4),ntri,isll,ifild)
      if (odebug(38)) then
       write(iwr,*) ' ump2 lagrangian contribution '
       call prtris(q(i4),num,iwr)
      endif
c
      call gmem_free_inf(i8,fnm,snm,'i8')
      call gmem_free_inf(i7,fnm,snm,'i7')
      call gmem_free_inf(i6,fnm,snm,'i6')
      call gmem_free_inf(i5,fnm,snm,'i5')
      call gmem_free_inf(i4,fnm,snm,'i4')
      call gmem_free_inf(i3,fnm,snm,'i3')
      call gmem_free_inf(i2,fnm,snm,'i2')
      call gmem_free_inf(i1,fnm,snm,'i1')
      call gmem_free_inf(i0,fnm,snm,'i0')
      call gmem_free_inf(iyuk,fnm,snm,'iyuk')
c
c   if ump3 then calculate ump2 2-pdm inside yuk3new
c
      if (uhf .and. mp3) return
c
c
      call ump2pd(q,iblk2d,ifil2d,mpblk(7),mpblk(9),mpblk(10),mpblk(5),
     +  mpblk(6),ifort1,ifort2,ifort5,ifort6,ifort3,ifort4,ifort7,
     +  ifort8)
      return
 6010 format (/' insufficient core for umpmky ',i9,' real words ',
     +        ' need ',i9,' real words ')
      end
      subroutine umpzy(ya,yb,z,nocca,noccb,ncoorb,mn,iblya,iblyb)
      implicit REAL  (a-h,o-z)
      dimension ya(ncoorb,ncoorb),yb(ncoorb,ncoorb),z(mn)
c
      nsq = ncoorb*ncoorb
      ifile1 = 1
      call rdedx(ya,nsq,iblya,ifile1)
      call rdedx(yb,nsq,iblyb,ifile1)
      iai = 0
      do 30 ia = nocca + 1 , ncoorb
         do 20 i = 1 , nocca
            iai = iai + 1
            ya(ia,i) = z(iai)
 20      continue
 30   continue
      call wrt3(ya,nsq,iblya,ifile1)
c
      do 50 ia = noccb + 1 , ncoorb
         do 40 i = 1 , noccb
            iai = iai + 1
            yb(ia,i) = z(iai)
 40      continue
 50   continue
      call wrt3(yb,nsq,iblyb,ifile1)
      return
      end
      subroutine ver_mp2b(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mp2b.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
