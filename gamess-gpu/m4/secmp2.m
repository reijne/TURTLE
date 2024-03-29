c 
c  $Author: jmht $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/secmp2.m,v $
c  $State: Exp $
c  
      subroutine ddmpa0(u,a,b,ibu,ifortr,ifortw,nocca,
     +  ncoorb,nsq,n3n,nov,nvirta,ifile1,ifint1)
      implicit REAL  (a-h,o-z)
c
      logical skip
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3)
      dimension u(ncoorb,ncoorb,n3n),a(ncoorb,ncoorb),b(nov*n3n)
      nocc1 = nocca + 1
      nov3 = nov*nunpr
      nat = n3n/3
      call search(1,ifortr)
      call search(1,ifortw)
      call search(1,ifint1)
      call rdedx(u,nsq*n3n,ibu,ifile1)
c
c
      do 50 ixx = n3n - 2 , n3n
         call vclr(u(1,1,ixx),1,ncoorb*ncoorb)
         do 40 k = 1 , nat - 1
            do 30 jj = 1 , ncoorb
               do 20 ii = 1 , ncoorb
                  u(ii,jj,ixx) = u(ii,jj,ixx) - u(ii,jj,ixx-3*k)
 20            continue
 30         continue
 40      continue
 50   continue
c
      do 80 ip = 1 , ncoorb
         do 70 iq = 1 , ip
            call rdedz(a,nsq,ifint1)
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
                  call rdedz(b,nov3,ifortr)
                  do 60 ix = 1 , nunpr
                     ixa = mapnu(ix)
                     iad = (ix-1)*nov + 1
                     call mxmb(a(1,nocc1),ncoorb,1,u(1,1,ixa),1,ncoorb,
     +                         b(iad),nocca,1,nvirta,ncoorb,nocca)
                     call mxmb(u(1,nocc1,ixa),ncoorb,1,a,ncoorb,1,b(iad)
     +                         ,nocca,1,nvirta,ncoorb,nocca)
 60               continue
                  call wtedz(b,nov3,ifortw)
               end if
            end if
 70      continue
 80   continue
      return
      end
      subroutine ddmpa1(u,a,b,ibu,ifortr,ifortw,nocca,
     +    ncoorb,nsq,n3n,nov,nvirta,w,nops,ifint1,ifile1)
      implicit REAL  (a-h,o-z)
c
      logical skip
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3)
      dimension u(ncoorb,ncoorb,n3n),a(ncoorb,ncoorb),b(nov*n3n)
      dimension w(nov*n3n)
      nocc1 = nocca + 1
      nov3 = nov*nunpr
      nat = n3n/3
      call search(1,ifortr)
      call search(1,ifortw)
      call search(1,ifint1)
      call rdedx(u,nsq*n3n,ibu,ifile1)
c
c
      do 50 ixx = n3n - 2 , n3n
         call vclr(u(1,1,ixx),1,ncoorb*ncoorb)
         do 40 k = 1 , nat - 1
            do 30 jj = 1 , ncoorb
               do 20 ii = 1 , ncoorb
                  u(ii,jj,ixx) = u(ii,jj,ixx) - u(ii,jj,ixx-3*k)
 20            continue
 30         continue
 40      continue
 50   continue
c
      do 80 ip = 1 , ncoorb
         do 70 iq = 1 , ip
            call rdedz(a,nsq,ifint1)
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
                  call rdedz(b,nov3,ifortr)
                  call ddmps1(b,w,nov,n3n,nops,ip,iq,nocca,ncoorb)
                  do 60 ix = 1 , nunpr
                     ixa = mapnu(ix)
                     iad = (ix-1)*nov + 1
                     call mxmb(a(1,nocc1),ncoorb,1,u(1,1,ixa),1,ncoorb,
     +                         b(iad),nocca,1,nvirta,ncoorb,nocca)
                     call mxmb(u(1,nocc1,ixa),ncoorb,1,a,ncoorb,1,b(iad)
     +                         ,nocca,1,nvirta,ncoorb,nocca)
 60               continue
                  call wtedz(b,nov3,ifortw)
               end if
            end if
 70      continue
 80   continue
      return
      end
      subroutine ddmpao(a,p,as,u,w,v,nocca,ncoorb,ntri,nsq,
     + n3n,isecdd,ifild,iblu,idintf,num,force,ftp,iblf,ifile1,ldebug,
     + iwrt)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
c     term constructed with a.o derivative integrals
c
      logical ldebug
      dimension a(ntri,n3n),p(ntri),as(nsq),u(ncoorb,ncoorb,n3n),
     1          w(nsq),v(nsq),force(n3n,n3n),ftp(n3n,n3n)
c     parameter (idblk=3120, iidblk=1950)
c     common/dbuf/g(idblk),vall(iidblk),icnt,mxtr
c     common/dlabs/i1(idblk),j1(idblk),k1(idblk),l1(idblk),
c    +            ix1(idblk)
      common/dbuf/icnt,mxtr,g(5118)
      common/dlabs/i1(5,3120)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
      equivalence (gg1,gg),(ikyj,jl),(ikyk,jk),(ikyi,il)
      data m0,zer0/0,0.0d0/
c
      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1
      call setsto(5*idblk,0,i1)
c
      call secget(isecdd,m0,isdd)
      call rdedx(p,ntri,isdd,ifild)
      call vclr(a,1,ntri*n3n)
c
      call search(1,idintf)
      iderb = idblk + iidblk + 2
 20   call reads(dbufs,iderb,idintf)
      if (icnt.eq.0) then
c
c
         call ddaosm(a,u,num,n3n,ntri,nsq)
c
         call rdedx(u,nsq*n3n,iblu,ifile1)
         do 50 ix = 1 , n3n
            call vclr(w,1,nsq)
            call mxmb(v,1,num,u(1,1,ix),1,ncoorb,w,1,num,num,ncoorb,
     +                nocca)
            call vclr(u(1,1,ix),1,nsq)
            call mxmb(w,1,num,v,num,1,u(1,1,ix),1,num,num,nocca,num)
            do 40 il = 1 , num
               do 30 is = 1 , il
                  u(il,is,ix) = u(il,is,ix) + u(is,il,ix)
                  u(is,il,ix) = u(il,is,ix)
 30            continue
 40         continue
 50      continue
c
         call vclr(force,1,n3n*n3n)
         do 60 ix = 1 , n3n
            call square(as,a(1,ix),ncoorb,ncoorb)
            call mxmb(u,nsq,1,as,1,nsq,force(ix,1),n3n,1,n3n,nsq,1)
 60      continue
         nat = n3n/3
         do 80 ix = 1 , n3n
            do 70 iy = 1 , ix
               force(ix,iy) = force(ix,iy) + force(iy,ix)
 70         continue
 80      continue
         call ddmpti(force,nat,nat)
         call rdedx(ftp,n3n*n3n,iblf,ifile1)
         do 100 iy = 1 , n3n
            do 90 ix = 1 , n3n
               ftp(ix,iy) = ftp(ix,iy) + force(ix,iy)
 90         continue
 100     continue
         call wrt3(ftp,n3n*n3n,iblf,ifile1)
c
         if (ldebug) write (iwrt,*) 'in ddmpao'
         if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
         return
      else
         call unpack(g(idblk1),lab816,i1,idblk*5)
         do 120 iw = 1 , icnt
            gg = g(iw)
            gg2 = gg
            gg3 = gg
            gg4 = gg
            gg5 = gg
            gg6 = gg
            gg7 = gg
            gg8 = gg + gg
            gg9 = gg8
            i = i1(1,iw)
            j = i1(2,iw)
            k = i1(3,iw)
            l = i1(4,iw)
            ix = i1(5,iw)
            ikyi = iky(i)
            ij = ikyi + j
            ik = ikyi + k
            il = ikyi + l
            ikyk = iky(k)
            kl = ikyk + l
            ikyj = iky(j)
            if (i.eq.j) then
               gg3 = zer0
               gg4 = zer0
               gg6 = zer0
               gg7 = zer0
               gg9 = gg
            end if
            if (ij.eq.kl) then
               gg5 = zer0
               gg6 = zer0
               gg7 = zer0
               gg9 = zer0
            end if
            if (k.eq.l) then
               gg2 = zer0
               gg4 = zer0
               gg7 = zer0
               gg8 = gg
            end if
            if (i.eq.k) gg1 = gg1 + gg5
            if (j.lt.k) then
               gg3 = gg6
            else if (j.eq.k) then
               gg3 = gg3 + gg6
            else
               jk = ikyj + k
               jl = ikyj + l
               go to 110
            end if
            jk = ikyk + j
            if (j.lt.l) then
               gg4 = gg7
               jl = iky(l) + j
            else if (j.eq.l) then
               gg4 = gg4 + gg7
               jl = ikyj + l
            else
               jl = ikyj + l
            end if
 110        a(ij,ix) = a(ij,ix) + (gg8+gg8)*p(kl)
            a(kl,ix) = a(kl,ix) + (gg9+gg9)*p(ij)
            a(ik,ix) = a(ik,ix) - gg1*p(jl)
            a(jl,ix) = a(jl,ix) - gg4*p(ik)
            a(il,ix) = a(il,ix) - gg2*p(jk)
            a(jk,ix) = a(jk,ix) - gg3*p(il)
c
 120     continue
         go to 20
      end if
      end
      subroutine ddaosm(a,as,num,n3n,ntri,nsq)
      implicit REAL  (a-h,o-z)
c
c     symmetrises a.o. term
c
INCLUDE(common/sizes)
INCLUDE(common/symtry)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      dimension a(ntri,n3n),as(nsq,n3n)
      logical skip
      common/scrtch/
     + ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     + dp(3,3,maxat),qd(6,3,maxat),
     + nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     + iperm(maxat*3),
     + ict(maxat,8),mptr(8,maxat*3),nuniq,nuni1,nuni2,nuni3
INCLUDE(common/mapper)
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
      nat = n3n/3
      if (nt.ne.1) then
c
c     skip symmetrisation if no symmetry !
c
         an = 1.0d0/dfloat(nt)
c
         do 80 n = 1 , nat
            do 70 nc = 1 , 3
               nprt = (n-1)*3 + nc
               do 20 ii = 1 , ntri
                  as(ii,nprt) = a(ii,nprt)
 20            continue
c
               do 60 im = 1 , num
                  do 50 in = 1 , im
                     imn = ind(im,in)
                     do 40 iop = 2 , nt
                        nuat = ict(n,iop)
                        im1 = ibasis(iop,im)
                        in1 = ibasis(iop,in)
                        imm = iabs(im1)
                        inn = iabs(in1)
                        isign1 = im1/imm
                        isign2 = in1/inn
                        isign = isign1*isign2
                        sign = dfloat(isign)
                        immnn = ind(imm,inn)
                        npnc = (iop-1)*3 + nc
                        do 30 k = 1 , 3
                           nuprt = (nuat-1)*3 + k
                           as(imn,nprt) = as(imn,nprt)
     +                        + sign*ptr(k,npnc)*a(immnn,nuprt)
 30                     continue
 40                  continue
 50               continue
 60            continue
 70         continue
 80      continue
c
c
         do 100 ii = 1 , n3n
            do 90 jj = 1 , ntri
               a(jj,ii) = an*as(jj,ii)
 90         continue
 100     continue
      end if
      if (nuniq.eq.0) return
      nuni1 = nuniq*3 - 2
      nuni2 = nuniq*3 - 1
      nuni3 = nuniq*3
      do 110 i = 1 , ntri
         a(i,nuni1) = 0.0d0
         a(i,nuni2) = 0.0d0
         a(i,nuni3) = 0.0d0
 110  continue
      do 130 n = 1 , nat
         if (n.ne.nuniq) then
            do 120 i = 1 , ntri
               a(i,nuni1) = a(i,nuni1) - a(i,n*3-2)
               a(i,nuni2) = a(i,nuni2) - a(i,n*3-1)
               a(i,nuni3) = a(i,nuni3) - a(i,n*3)
 120        continue
         end if
 130  continue
      return
      end
      subroutine ddmpcr(e,cb,ex,force,ftp,eps,z,t,nov,n3n,nocca,
     + nvirta,ncoorb,ible,ibleps,iblf,ifcb,ifex,nops,ifile1,istrma,
     + ldebug,iwrt)
      implicit REAL  (a-h,o-z)
c
      logical ldebug
INCLUDE(common/sizes)
      dimension e(ncoorb),cb(nov,n3n),ex(nov,n3n),force(n3n,n3n)
      dimension ftp(n3n,n3n),eps(ncoorb,ncoorb,n3n),z(nov,n3n),t(nov)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),nunpr
c
      call search(1,ifex)
      call search(1,ifcb)
      call search(1,istrma)
c
      n3n3 = n3n - 3
c     nov3 = nov*n3n3
      nread = nov*nunpr
      call rdedx(eps,ncoorb*ncoorb*n3n,ibleps,ifile1)
      call rdedx(e,ncoorb,ible,ifile1)
      call rdedx(ftp,n3n*n3n,iblf,ifile1)
c
c
      do 40 ix = 1 , n3n3
         do 30 j = 1 , nocca
            do 20 i = 1 , nocca
               eps(i,j,ix) = -4.0d0*eps(i,j,ix)
 20         continue
 30      continue
 40   continue
      do 70 ix = 1 , n3n3
         do 60 ib = nocca + 1 , ncoorb
            do 50 ia = nocca + 1 , ncoorb
               eps(ia,ib,ix) = 4.0d0*eps(ia,ib,ix)
 50         continue
 60      continue
 70   continue
c
c
      call vclr(force,1,n3n*n3n)
      nocc1 = nocca + 1
c
      do 130 ia = nocc1 , ncoorb
         do 120 i = 1 , nocca
            eai = e(ia) - e(i)
            call rdedz(cb,nread,ifcb)
            call rdedz(ex,nread,ifex)
            call ddmps0(cb,z,nov,n3n,nops,ia,i,nocca,ncoorb)
            call ddmps0(ex,z,nov,n3n,nops,ia,i,nocca,ncoorb)
            call rdedz(t,nov,istrma)
            call vclr(z,1,nov*n3n)
c
            do 80 ix = 1 , n3n3
               call mxmb(t,1,nocca,eps(nocc1,nocc1,ix),1,ncoorb,z(1,ix),
     +                   1,nocca,nocca,nvirta,nvirta)
               call mxmb(t,nocca,1,eps(1,1,ix),1,ncoorb,z(1,ix),nocca,1,
     +                   nvirta,nocca,nocca)
 80         continue
c
c
            ibj = 0
            do 110 ib = nocc1 , ncoorb
               do 100 j = 1 , nocca
                  ibj = ibj + 1
                  ddiff = 1.0d0/(eai+e(ib)-e(j))
                  do 90 ix = 1 , n3n3
                     ex(ibj,ix) = (ex(ibj,ix)-cb(ibj,ix)-cb(ibj,ix)+z(
     +                            ibj,ix))*ddiff
 90               continue
 100           continue
 110        continue
c
            call mxmb(cb,nov,1,ex,1,nov,force,1,n3n,n3n3,nov,n3n3)
 120     continue
 130  continue
c
c
      if (ldebug) write (iwrt,*) 'in ddmpcr'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      do 150 ix = 1 , n3n
         do 140 iy = 1 , ix
            ftp(ix,iy) = ftp(ix,iy) + force(ix,iy) + force(iy,ix)
 140     continue
 150  continue
      call wrt3(ftp,n3n*n3n,iblf,ifile1)
c
c
      return
      end
      subroutine ddmpen(istrm9,force,n3n,ifild,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
INCLUDE(common/restri)
      dimension force(n3n,n3n)
      ifile1 = 1
      nat = n3n/3
      call rdedx(force,n3n*n3n,istrm9,ifile1)
c
      do 30 ix = 1 , n3n
         do 20 iy = 1 , ix - 1
            force(iy,ix) = force(ix,iy)
 20      continue
 30   continue
c
      do 60 ix = 1 , n3n - 3
         do 50 iy = n3n - 2 , n3n
            force(iy,ix) = 0.0d0
            do 40 k = 1 , nat - 1
               force(iy,ix) = force(iy,ix) - force(iy-3*k,ix)
 40         continue
 50      continue
 60   continue
c
      call vclr(force(1,n3n-2),1,3*n3n)
      do 90 ix = n3n - 2 , n3n
         do 80 iy = 1 , n3n
            do 70 k = 1 , nat - 1
               force(iy,ix) = force(iy,ix) - force(iy,ix-3*k)
 70         continue
 80      continue
 90   continue
c
      if(ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      isc110 = isect(110)
      m110 = 110
      len = lensec(n3n*n3n)
      call secput(isc110,m110,len,iblok)
      call wrt3(force,n3n*n3n,iblok,ifild)
      call revind
      return
      end
      subroutine ddmp00(igradf,nat3,a)
c---------------------------------------------------------------
c
c     part2
c
c     initialisation and handling of derivative integrals
c
c---------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      dimension a(nat3*nat3)
      call vclr(a,1,nat3*nat3)
      ifile1 = 1
      call wrt3(a,nat3*nat3,igradf,ifile1)
      return
      end
      subroutine ddmptt(iso,nshels,iwr)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
INCLUDE(common/nshel)
INCLUDE(common/symtry)
      logical skip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
INCLUDE(common/infoa)
      dimension iso(nshels,*)
      data one/1.0d0/
c
c     now construct transformation table of perturbations
c
c
      do 40 ii = 1 , nshell
         ic = katom(ii)
         do 30 it = 1 , nt
            id = iso(ii,it)
            ict(ic,it) = katom(id)
 30      continue
 40   continue
c
c
      nat3 = 3*nat
      do 50 iii = 1 , nat3
         mptr(1,iii) = iii
 50   continue
c
      do 110 iat = 1 , nat
         do 100 iop = 2 , nt
            nuat = ict(iat,iop)
            npp = (iop-1)*3
            do 90 j = 1 , 3
               ipr = (iat-1)*3 + j
               do 60 k = 1 , 3
                  tr = ptr(k,npp+j)
                  if (dabs(tr-one).lt.1.0d-8) go to 70
                  if (dabs(tr+one).lt.1.0d-8) go to 80
 60            continue
c
c     not supposed to be here!
c
 70            mptr(iop,ipr) = (nuat-1)*3 + k
               go to 90
 80            mptr(iop,ipr) = -((nuat-1)*3+k)
 90         continue
 100     continue
 110  continue
c
c
      write (iwr,6020)
      do 120 i = 1 , nat3
         write (iwr,6010) (mptr(iop,i),iop=1,nt)
 120  continue
      return
 6010 format (10x,8i5)
 6020 format(/
     + 1x,'transformation table for perturbations'/
     + 1x,'======================================'/)
      end
      subroutine mp2ddrv(q,iq)
      implicit REAL  (a-h,o-z)
c
c     driving routine for mp2 second derivatives
c
      dimension iq(*),q(*)
INCLUDE(common/sizes)
      common/mpases/mpass,ipass,mpst,mpfi,mpadd,iaoxf1,iaoxf2,moderf
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/symtry)
      logical xskip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,xskip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3),nuniq,nuni1,nuni2,nuni3
      common/maxlen/maxq
c
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/iofile)
c
      common/small/eigs(maxorb)
      logical lstop,skipp
      common/lsort/skipp(maxat*3)
INCLUDE(common/atmblk)
INCLUDE(common/prnprn)
      character*10 charwall
      character *8 grad
      data grad/'gradient'/
c
      cutoff = 10.0d0**(-icut)
      l100 = 70 + lenint(60)
      mn = nocca*nvirta
c     nij = nocca*(nocca+1)/2
c     nab = nvirta*(nvirta+1)/2
      nat3 = nat*3
c     nat33 = nat3
      ntri = ncoorb*(ncoorb+1)/2
      nsq = ncoorb*ncoorb
      ifile1 = 1
      ifint1 = 17
      ifint2 = 18
      ifort2 = 19
      istrma = 20
      istrmd = 21
      call search(1,ifint1)
      call search(1,ifint2)
      call search(1,ifort2)
c     ndep = nab + nij
      if (odebug(23). and. runtyp.ne.grad .and. .not.lopti 
     +  .and. .not.lforce) write (iwr,6140) mpflag
      m103 = isect(103)
      ityp = 0
      if (mpflag.eq.6 .and. runtyp.ne.grad .and. .not.lopti .and.
     +    .not.lforce) then
c
         call secget(m103,ityp,iblok)
         call rdedx(cigr,l100,iblok,ifild)
c        symmetry
c     ------------------------------------------------------
c     read in transformation matrices for s,p,d,f functions.
c     ------------------------------------------------------
c
         iso = igmem_alloc(nw196(5))
         call rdedx(ptr,nw196(1),ibl196(1),ifild)
         if (odbas) call rdedx(dtr,nw196(2),ibl196(2),ifild)
         if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),ifild)
         if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),ifild)
         call rdedx(q(iso),nw196(5),ibl196(5),ifild)
c
         call ddmptt(q(iso),nshell,iwr)
         call gmem_free(iso)
c
         i1 = igmem_alloc(nat3*nat3)
         call ddmp00(mpblk(1),nat3,q(i1))
         call gmem_free(i1)
c
         call timit(3)
         iatr   = igmem_alloc(ntri*nat3)
         ipy    = igmem_alloc(ntri)
         iasq   = igmem_alloc(nsq)
         iux    = igmem_alloc(nsq*nat3)
         iwks   = igmem_alloc(nsq)
         ivec   = igmem_alloc(nsq)
         iforce = igmem_alloc(nat3*nat3)
         iftp   = igmem_alloc(nat3*nat3)
c
         call secget(isect(8),8,isec8)
         isec8 = isec8 + mvadd
         call rdedx(q(ivec),num*ncoorb,isec8,ifild)
c
         call ddmpao(q(iatr),q(ipy),q(iasq),q(iux),q(iwks),q(ivec),
     +  nocca,ncoorb,ntri,nsq,nat3,isecdd,ifild,mpblk(6),istrmd,
     +  num,q(iforce),q(iftp),mpblk(1),ifile1,odebug(23),iwr)
c
         call gmem_free(iftp)
         call gmem_free(iforce)
         call gmem_free(ivec)
         call gmem_free(iwks)
         call gmem_free(iux)
         call gmem_free(iasq)
         call gmem_free(ipy)
         call gmem_free(iatr)
c
         write (iwr,6130) cpulft(1) ,charwall()
         call timit(3)
c
c     call of ncmkfc moved to after derivative integral
c     handling to save file space
c
c   perturbation selection
c
         if (odebug(23)) call wrtskp(nat3,iwr)
c
c     multi pass code
c
c >>>>> mpass set so perturbations sorted one at a time
c
         mpass = nat
         iaoxf1 = istrmd
         moderf = ifort2
         if (mpass.eq.1) iaoxf2 = istrmd
         if (mpass.gt.1) then
            iaoxf2 = 22
            call search(1,iaoxf2)
         end if
c
c
         write (iwr,6120)cpulft(1) ,charwall()
         ipass = 0
         do 20 mmmmm = 1 , nat
            mpst = (mmmmm-1)*3+1
            mpfi = mmmmm * 3
            if (.not.(xskip(mpst))) then
               write(iwr,6150) mpst, mpfi
               call setbfa
               ipass = ipass + 1
               call mpsrt7(q,iq,num,iwr)
               call closbf(0)
c
               i1  = igmem_alloc_all(maxa)
               ii1 = lenrel(i1-1)+1
               call mp2dri(ifint2,ifile1,istrma,q(i1),iq(ii1),maxa)
               call gmem_free(i1)
               call delfil(ifint2)
               call timit(3)
               write (iwr,6110) cpulft(1) ,charwall()
            end if
c
 20      continue
c
         call delfil(iaoxf1)
         if (mpass.gt.1) call delfil(iaoxf2)
         call mpsrtj(q,iq,ncoorb,kblk(1),nufile(1),ifint1)
         ipss = 1
         call mpsrtk(q,iq,ncoorb,kblk(1),nufile(1),ifint2,ipss)
         i1 = igmem_alloc_all(maxa)
         call ddmpfc(ncoorb,nat3,nvirt,nocc,mpblk(7),mpblk(6),mpblk(5),
     +  mpblk(2),mpblk(1),q(i1),mpblk(3),mpblk(4),mpblk(8),mpblk(9),
     +  ifile1,ifint1,ifint2,odebug(23),iwr,maxa)
         call gmem_free(i1)
         call delfil(ifint2)
         write (iwr,6090)cpulft(1) ,charwall()
         call timit(3)
c
         iforce = igmem_alloc(nat3*nat3)
         iftp   = igmem_alloc(nat3*nat3)
         iu     = igmem_alloc(nat3*nsq)
         if     = igmem_alloc(nat3*nsq)
         call ddmprs(q(iforce),q(iftp),q(iu),q(if),mpblk(1),mpblk(2),
     +               mpblk(6),nat3,ncoorb,ifile1,odebug(23),iwr)
         call gmem_free(if)
         call gmem_free(iu)
         call gmem_free(iftp)
         call gmem_free(iforce)
c
         nov = nocca*nvirta
         call mpsrtd(q,iq,nocca,nvirta,ifort2,ifort2,nat3,iwr,cutoff)
c
         iu   = igmem_alloc(nsq*nat3)
         iaa  = igmem_alloc(nsq)
         ib   = igmem_alloc(nov*nat3)
         iwks = igmem_alloc(nov*nat3)
         write (iwr,6060) cpulft(1) ,charwall()
         call ddmpa1(q(iu),q(iaa),q(ib),mpblk(6),ifort2,istrmd,nocca,
     +  ncoorb,nsq,nat3,nov,nvirta,q(iwks),nt,ifint1,ifile1)
         call gmem_free(iwks)
         call gmem_free(ib)
         call gmem_free(iaa)
         call gmem_free(iu)
c
         call mpsrt8(q,iq,istrmd,istrmd,nov,iwr,cutoff)
c
         iu   = igmem_alloc(nsq*nat3)
         iaa  = igmem_alloc(nsq)
         ib   = igmem_alloc(nov*nat3)
         iwks = igmem_alloc(nov*nat3)
         call ddmpa0(q(iu),q(iaa),q(ib),mpblk(6),istrmd,ifort2,nocca,
     +  ncoorb,nsq,nat3,nov,nvirta,ifile1,ifint1)
         call gmem_free(iwks)
         call gmem_free(ib)
         call gmem_free(iaa)
         call gmem_free(iu)
c
         call delfil(ifint1)
c
         call mpsrt9(q,iq,ifort2,istrmd,nov,nocca,ncoorb,iwr,cutoff)
c
         iforce = igmem_alloc(nat3*nat3)
         ie     = igmem_alloc(ncoorb)
         icb    = igmem_alloc(nov*nat3)
         iex    = igmem_alloc(nov*nat3)
         iftp   = igmem_alloc(nat3*nat3)
         ieps   = igmem_alloc(nsq*nat3)
         iz     = igmem_alloc(nov*nat3)
         it     = igmem_alloc(nov)
         call ddmpcr(q(ie),q(icb),q(iex),q(iforce),q(iftp),q(ieps),
     +  q(iz),q(it),nov,nat3,nocca,nvirta,ncoorb,mpblk(4),mpblk(5),
     +  mpblk(1),ifort2,istrmd,nt,ifile1,istrma,odebug(23),iwr)
         call gmem_free(it)
         call gmem_free(iz)
         call gmem_free(ieps)
         call gmem_free(iftp)
         call gmem_free(iex)
         call gmem_free(icb)
         call gmem_free(ie)
c
         call timit(3)
         write (iwr,6070) cpulft(1) ,charwall()
         if(odebug(23)) write (iwr,6080)
         call ddmpen(mpblk(1),q(iforce),nat3,ifild,odebug(23),iwr)
         call gmem_free(iforce)
         lds(isect(110)) = lensec(nat*nat*9)
         mpflag = 2
         call revise
         call delfil(istrma)
         call delfil(istrmd)
         call delfil(ifort2)
         call revind
c        call whtps
         return
      else
c
c    set up blocks for ed0 for
c    1) gradient matrix (n3n,n3n)
c    2) workspace for the (rs/bj) derivative integral term
c    3) blank
c    4) e (ncoorb)
c    5) eder(ncoorb,ncoorb,n3n)
c    6) u(ncoorb,ncoorb,n3n)
c    7) y(ncoorb,ncoorb)
c    8) w(ncoorb,ncoorb)
c    9) w2(ncoorb,ncoorb)
c    10) wtrans(nbf,nbf)
c    11) ytrans(nbf,nbf)
c    12) dm(ncoorb,ncoorb)
c
         mpblk(1) = 1
         mpblk(2) = mpblk(1) + lensec(nat3*nat3)
         mpblk(3) = mpblk(2) + lensec(nat3*nsq)
         mpblk(4) = mpblk(3)
         mpblk(5) = mpblk(4) + lensec(ncoorb)
         mpblk(6) = mpblk(5) + lensec(nsq*nat3)
         mpblk(7) = mpblk(6) + lensec(nsq*nat3)
         mpblk(8) = mpblk(7) + lensec(nsq)
         mpblk(9) = mpblk(8) + lensec(nsq)
         mpblk(10) = mpblk(9) + lensec(nsq)
         mpblk(11) = mpblk(10) + lensec(num*num)
         mpblk(12) = mpblk(11) + lensec(num*num)
         mpblk(13) = mpblk(12) + lensec(nsq)
         call revise
         call mpsrtj(q,iq,ncoorb,kblk(1),nufile(1),ifint1)
         ipss = 1
         call mpsrtk(q,iq,ncoorb,kblk(1),nufile(1),ifint2,ipss)
         i1 = igmem_alloc_all(maxq)
         call chfcls(q(i1),maxq)
         call gmem_free(i1)
         write (iwr,6050)cpulft(1) ,charwall()
         call timit(3)
c
c   read in eigenvalues
c
         m9 = 9
         call secget(isect(9),m9,isec9)
         i1 = 1 + ncoorb
         i2 = i1 + nsq
         i3 = i2 + nsq
         i4 = i3 + nsq
         i5 = i4 + nsq
         i6 = i5 + nsq
         i7 = i6 + nsq
         itop = i7 + nocca*nvirta
         if (itop.gt.maxq) then
            write (iwr,6010) maxq , itop
            call caserr(' not enough core')
         end if
         i1 = igmem_alloc(ncoorb)
         i2 = igmem_alloc(nsq)
         i3 = igmem_alloc(nsq)
         i4 = igmem_alloc(nsq)
         i5 = igmem_alloc(nsq)
         i6 = igmem_alloc(nsq)
         i7 = igmem_alloc(nocca*nvirta)
         call rdedx(q(i1),ncoorb,isec9,ifild)
         call mpmaky(q(i2),q(i3),q(i1),ncoorb,q(i4),q(i5),nocca,
     +  mpblk(7),nvirta,q(i6),mpblk(8),iblks,ifils,q(i7),mn,ifile1,
     +  istrma,ifint1,ifint2)
         call gmem_free(i7)
         call gmem_free(i6)
         call gmem_free(i5)
         call gmem_free(i4)
         call gmem_free(i3)
         call gmem_free(i2)
         call gmem_free(i1)
c
         m9 = 9
         call secget(isect(9),m9,isec9)
         call rdedx(eigs,ncoorb,isec9,ifild)
         ieps = igmem_alloc(nocca*nvirta)
         do 40 i = 1 , nocca
            do 30 iaa = nocca + 1 , ncoorb
               iai = (iaa-nocca-1)*nocca + i + ieps-1
               q(iai) = 1.0d0/(eigs(iaa)-eigs(i))
 30         continue
 40      continue
         np = 1
         npstar = 0
         skipp(1) = .false.
         lstop = .false.
c
c     solve for z-matrix
c
         call chfdrv(q(ieps),lstop,skipp)
         call gmem_free(ieps)
c
         call delfil(nofile(1))
         write (iwr,6040) cpulft(1) ,charwall()
         iblz = iblks + lensec(mn)
         iy   = igmem_alloc(nsq)
         iz   = igmem_alloc(nocca*nvirta)
         iww  = igmem_alloc(nsq)
         iww2 = igmem_alloc(nsq)
         ia1  = igmem_alloc(nsq)
         ia2  = igmem_alloc(nsq)
         ie   = igmem_alloc(ncoorb)
         m9 = 9
         call secget(isect(9),m9,isec9)
         call rdedx(q(ie),ncoorb,isec9,ifild)
         call mpmakw(q(iy),q(iz),q(iww),q(iww2),mpblk(7),iblz,mpblk(8),
     +  mpblk(9),q(ia1),q(ia2),nocca,nvirta,ncoorb,ifils,q(ie),ifint1,
     +  ifint2,ifile1)
         call gmem_free(ie)
         call gmem_free(ia2)
         call gmem_free(ia1)
         call gmem_free(iww2)
         call gmem_free(iww)
         call gmem_free(iz)
         call gmem_free(iy)
c
c     y,w,w2 matrices now stored on ed0 in the mo basis
c
         ityp = 0
         len = lensec(ntri)
c
         i1 = igmem_alloc(nsq)
         i2 = igmem_alloc(nsq)
         i3 = igmem_alloc(nsq)
         call secget(isect(8),8,iblok)
         iblok = iblok + mvadd
         call rdedx(q(i2),num*ncoorb,iblok,ifild)
         call rdedx(q(i1),nsq,mpblk(7),ifile1)
         call vtamv(q(i1),q(i2),q(i3),ncoorb)
         call wrt3(q(i1),nsq,mpblk(11),ifile1)
         call trsqsq(q(i1),q(i3),ncoorb)
         call secput(isecdd,ityp,len,isdd)
         call wrt3(q(i3),ntri,isdd,ifild)
         call rdedx(q(i1),nsq,mpblk(8),ifile1)
         call vtamv(q(i1),q(i2),q(i3),ncoorb)
         call wrt3(q(i1),nsq,mpblk(10),ifile1)
c
c   now because for some funny reason we are
c   calculating -w change the sign
c
         do 50 iijj = i1 , nsq+i1-1
            q(iijj) = -q(iijj)
 50      continue
         call trsqsq(q(i1),q(i3),ncoorb)
         call secput(isecll,ityp,len,isll)
         call wrt3(q(i3),ntri,isll,ifild)
         call gmem_free(i3)
         call gmem_free(i2)
         call gmem_free(i1)
c
         write (iwr,6020) cpulft(1) ,charwall()
         i1 = igmem_alloc_all(maxa)
         call ddmptr(q(i1),iblk2d,ifil2d,mpblk(11),mpblk(12),ifile1,
     +               ifort2,ifint1,ifint2,maxa)
         call gmem_free(i1)
         write (iwr,6030) cpulft(1) ,charwall()
         call timit(3)
c
c**** end of first entry
c**** all information available for a gradient at this point
c**** second entry only for second derivatives
         call revind
         call delfil(ifort2)
         call delfil(ifint1)
         call delfil(ifint2)
         if (runtyp.eq.grad .or. lopti .or. lforce) then
            call delfil(1)
         end if
         return
      end if
 6010 format (/' insufficient core for mpmaky : have ',
     +  i8,' real words; need ',i8,' real words')
6020  format(/1x,
     +   '1-particle gradient density matrix complete at ',f8.2,
     +   ' seconds',a10,' wall')
6030  format(/1x,
     +   '2-particle gradient density matrix complete at ',f8.2,
     +   ' seconds',a10,' wall')
6040  format(/1x,'first set of chf equations complete at ',
     +   f8.2,' seconds',a10,' wall')
6050  format(/1x,'sorting complete at ',f8.2,' seconds',a10,' wall')
6060  format(/1x,'commence final assembly of mp2 contribution at ',
     +   f8.2,' seconds',a10,' wall')
6070  format(/1x,'assembly complete at ',f8.2,' seconds',a10,' wall')
6080  format(//1x,'mp2 contribution to force constant matrix'/
     +         1x,'=========================================')
6090  format(/1x,'end of ddmpfc at ',f8.2,' seconds',a10,' wall')
6110  format( 1x,'sorting complete at ',f8.2,' seconds',a10,' wall')
6120  format(/1x,
     +  'commence sorting of perturbations at ',f8.2,' seconds',
     +   a10,' wall')
6130  format(/1x,
     +  'a.o. derivative integral term complete at ',f8.2,' seconds',
     +   a10,' wall')
6140  format(/1x,' restart flag = ',i3)
6150  format(/1x,' sorting of perturbations ', i3,' to ',i3)
      end
      subroutine ddmps2(iso,nshels,e,a1,a2,b,dum,vec,r,mapb,lmap,
     + ifint1,ifint2,ifort2,nocca,nvirta,ncoorb,icount)
c
      implicit REAL  (a-h,o-z)
      logical lmap
INCLUDE(common/sizes)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
      dimension e(ncoorb),a1(ncoorb,ncoorb),a2(ncoorb,ncoorb)
      dimension b(nvirta,nocca),dum(ncoorb*ncoorb)
      dimension vec(ncoorb*ncoorb),r(ncoorb*ncoorb)
      dimension lmap(ncoorb*ncoorb),mapb(ncoorb*ncoorb)
      dimension m0(48),iso(nshels,*)
      logical ijump,jjump
c
      icount = 0
      do 90 ii = 1 , nshell
         ijump = .false.
         do 30 it = 1 , nt
            id = iso(ii,it)
            ijump = ijump .or. id.gt.ii
            m0(it) = id
 30      continue
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
         do 80 jj = 1 , ii
            jjump = .false.
            do 50 it = 1 , nt
               id = m0(it)
               jd = iso(jj,it)
               jjump = jjump .or. jd.gt.ii
               if (id.lt.jd) then
                  nd = id
                  id = jd
                  jd = nd
               end if
               jjump = jjump .or. (id.eq.ii .and. jd.gt.jj)
 50         continue
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
            do 70 i = mini , maxi
               i1 = loci + i
               do 60 j = minj , maxj
                  j1 = locj + j
c
                  ij1 = (j1-1)*ncoorb + i1
                  icount = icount + 1
                  mapb(icount) = ij1
                  lmap(icount) = ijump .or. jjump
 60            continue
 70         continue
c
 80      continue
 90   continue
      call search(1,ifort2)
      call search(1,ifint1)
      call search(1,ifint2)
      nsq = ncoorb*ncoorb
c     nov = nocca*nvirta
c
      do 180 ip = 1 , ncoorb
         do 170 iq = 1 , ip
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            if (ip.gt.nocca) then
               if (iq.le.nocca) then
                  ebj = e(ip) - e(iq)
                  do 110 i = 1 , nocca
                     do 100 ia = nocca + 1 , ncoorb
                        b(ia-nocca,i) = (a1(i,ia)*4.0d0-a2(i,ia)-a2(i,ia
     +                                  ))/(e(ia)+ebj-e(i))
 100                 continue
 110              continue
                  iv1 = ncoorb*nocca + 1
c
                  call vclr(dum,1,ncoorb*ncoorb)
                  call mxmb(vec(iv1),1,ncoorb,b,1,nvirta,dum,nocca,1,
     +                      ncoorb,nvirta,nocca)
                  call vclr(r,1,ncoorb*ncoorb)
                  call mxmb(vec,1,ncoorb,dum,1,nocca,r,1,ncoorb,ncoorb,
     +                      nocca,ncoorb)
c
                  do 120 im = 1 , ncoorb
                     dum((im-1)*ncoorb+im) = r((im-1)*ncoorb+im)
 120              continue
                  do 140 im = 2 , ncoorb
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                     do 130 in = 1 , im - 1
                        imn = (in-1)*ncoorb + im
                        inm = (im-1)*ncoorb + in
                        dum(imn) = (r(imn)+r(inm))*0.5d0
                        dum(inm) = dum(imn)
 130                 continue
 140              continue
c
                  do 150 i = 1 , icount
                     r(i) = dum(mapb(i))
 150              continue
                  do 160 i = 1 , icount
                     if (lmap(i)) r(i) = 0.0d0
 160              continue
c
                  call wtedz(r,icount,ifort2)
               end if
            end if
 170     continue
 180  continue
      return
      end
      subroutine ddmps3(vec,y,ysym,ifytr,ifdm,ncoorb,nocca,ifile1)
c
      implicit REAL  (a-h,o-z)
      dimension y(ncoorb*ncoorb),ysym(ncoorb*ncoorb)
      dimension vec(ncoorb*ncoorb)
c
      call rdedx(y,ncoorb*ncoorb,ifytr,ifile1)
c
      do 30 i = 1 , ncoorb
         do 20 j = 1 , ncoorb
            ij = (j-1)*ncoorb + i
            ji = (i-1)*ncoorb + j
            ysym(ij) = (y(ij)+y(ji))*0.5d0
 20      continue
 30   continue
c
      call vclr(y,1,ncoorb*ncoorb)
      call mxmb(vec,1,ncoorb,vec,ncoorb,1,y,1,ncoorb,ncoorb,nocca,
     +          ncoorb)
c
      call wrt3(ysym,ncoorb*ncoorb,ifytr,ifile1)
      call wrt3(y,ncoorb*ncoorb,ifdm,ifile1)
c
c
      return
      end
      subroutine ddmps4(iso,nshels,vec,dum,d,b,iblw,ifw,ifort,nocca,
     +   nvirta,ncoorb,ys,css,ifytr,ifdm,ntri,wks,ifile1,cut)
c
      implicit REAL  (a-h,o-z)
      dimension iso(nshels,*)
      dimension d(nocca*nvirta),vec(ncoorb*ncoorb),
     +  dum(ncoorb*ncoorb),wks(ncoorb,ncoorb),b(*)
      dimension ys(ncoorb*ncoorb),css(ncoorb*ncoorb)
      logical ijump,jjump
INCLUDE(common/sizes)
INCLUDE(common/symtry)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
INCLUDE(common/atmblk)
_IFN1(iv)      common/craypk/labout(1360)
_IF1(iv)      common/craypk/labij(340),labkl(340)
      common/blkin/g(510),nint,nxtr
      dimension m0(48)
      logical lab,labc,labcd
      data mzero/0/
c
_IFN1(civ)      call izero(1360,labout,1)
_IF1(iv)      call setsto(680,0,labij)
      nov = nocca*nvirta
      call rdedx(ys,ncoorb*ncoorb,ifytr,ifile1)
      call rdedx(css,ncoorb*ncoorb,ifdm,ifile1)
c
      call search(1,ifort)
      call search(iblw,ifw)
      nint = 0
c
      do 160 ii = 1 , nshell
         ijump = .false.
         do 30 it = 1 , nt
            id = iso(ii,it)
            ijump = ijump .or. id.gt.ii
            m0(it) = id
 30      continue
         iceni = katom(ii)
         do 150 jj = 1 , ii
            if (.not.(ijump)) then
               jjump = .false.
               do 50 it = 1 , nt
                  id = m0(it)
                  jd = iso(jj,it)
                  jjump = jjump .or. jd.gt.ii
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
c
            ntimes = (maxi-mini+1)*(maxj-minj+1)
            imax = loci + maxi
            do 80 itimes = 1 , ntimes
               ib1 = (itimes-1)*ntri
c
               call rdedz(d,nov,ifort)
               if (.not.(ijump .or. jjump)) then
                  iv1 = ncoorb*nocca + 1
c
                  call vclr(dum,1,ncoorb*ncoorb)
                  call mxmb(vec(iv1),1,ncoorb,d,nocca,1,dum,1,ncoorb,
     +                      imax,nvirta,nocca)
                  call vclr(wks,1,ncoorb*ncoorb)
                  call mxmb(vec,1,ncoorb,dum,ncoorb,1,wks,1,ncoorb,imax,
     +                      nocca,imax)
c
                  do 70 ms1 = 1 , imax
                     do 60 ms2 = 1 , ms1
                        ms12 = iky(ms1) + ms2
                        b(ib1+ms12) = (wks(ms1,ms2)+wks(ms2,ms1))*0.5d0
 60                  continue
 70               continue
               end if
c
c
 80         continue
            if (.not.(ijump .or. jjump)) then
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
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                                 do 90 l = minl , maxl
                                    l1 = locl + l
c
                                    if (k1.ge.l1) then
                                       ikl = iky(k1) + l1
                                    else
                                       ikl = iky(l1) + k1
                                    end if
                                    ipos = (icount-1)*ntri
                                    val = -b(ipos+ikl)
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
                                    vv1 = ys(imn)*css(ils)
                                    vv2 = ys(ils)*css(imn)
                                    val = val + vv1 + vv2
c
c
                                    vv1 = ys(iml)*css(ins)
                                    vv2 = ys(ims)*css(inl)
                                    vv3 = ys(inl)*css(ims)
                                    vv4 = ys(ins)*css(iml)
                                    val = val - (vv1+vv2+vv3+vv4)*0.25d0
c
c
                                    if (dabs(val).gt.cut) then
                                       nint = nint + 1
                                       g(nint) = val
_IFN1(iv)                                       nint4 = nint + nint + nint + nint
_IFN1(iv)                                       labout(nint4-3) = i1
_IFN1(iv)                                       labout(nint4-2) = j1
_IFN1(iv)                                       labout(nint4-1) = k1
_IFN1(iv)                                       labout(nint4) = l1
_IF1(iv)                                      labij(nint) = j1 + i4096(i1)
_IF1(iv)                                      labkl(nint) = l1 + i4096(k1)
                                       if (nint.eq.num2e) then
_IFN1(iv)                                     call pack(g(num2e+1),lab816,
_IFN1(iv)     +                                  labout,numlab)
_IF1(iv)                                      call pak4v(labij,g(num2e+1))
                                         call put(g,m511,ifw)
                                         nint = 0
_IFN1(civ)                                   call izero(1360,labout,1)
_IF1(iv)                                   call setsto(680,0,labij)
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
c
_IFN1(iv)      call pack(g(num2e+1),lab816,labout,numlab)
_IF1(iv)      call pak4v(labij,g(num2e+1))
      call put(g,m511,ifw)
      call put(g,mzero,ifw)
      call clredx
      return
      end
      subroutine ddmpfc(ncoorb,nat3,nvir,noc,isty,istu,isteps,istt,
     +  istg,a,istt2,isteig,ittw,ittw2,ifile1,ifint1,ifint2,
     +  ldebug,iwrt,maxq)
c-------------------------------------------------------------------
c
c     part 3
c
c     ncmkfc group
c-------------------------------------------------------------------
      implicit REAL  (a-h,o-z)
      logical ldebug
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1   istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension a(*)
c      npass=2
      igradf = istg
      istrmy = isty
      istrmu = istu
      itw = ittw
      itw2 = ittw2
      istrme = isteps
      iste = isteig
      istrmt = istt
      ist2 = istt2
      norbs = ncoorb
      n3n = nat3
      nvirt = nvir
      nocc = noc
      nocc1 = nocc + 1
      nsq = norbs*norbs
      n3nn2 = norbs*norbs*n3n
      iforce = 1
      iforc1 = iforce + n3n*n3n
      iforc2 = iforc1 + n3n*n3n
      iftp = iforc2 + n3n*n3n
      iu = iftp + n3n*n3n
      iepsd = iu + n3nn2
      iwks1 = iepsd + n3nn2
      iwks2 = iwks1 + nsq
      iy = iwks2 + nsq
      ieps = iy + nsq
      iww = ieps + norbs
      iw2 = iww + nsq
      itt = iw2 + nsq
      if3n = itt + nocc*nvirt
      itop = if3n + n3nn2
      call ddmpf1(a(iy),a(iu),a(iepsd),a(iwks1),a(ieps),a(iww),a(iw2),
     +  a(iforce),a(iftp),ifile1,ldebug,iwrt)
      if (itop.gt.maxq) then
         write (iwrt,*) 'ddmpfc can not be done in one pass'
         if (if3n.gt.maxq) call caserr('not possible in two passes')
         call ddmpf3(a(iu),a(iw2),a(iwks1),a(iy),a(iforce),a(iftp),
     +  a(ieps),a(iww),a(iepsd),a(iwks2),a(iforc1),a(itt),a(iforc2),
     +  ifile1,ifint1,ifint2,ldebug,iwrt)
         call ddmpf4(a(iu),a(iw2),a(iwks1),a(iy),a(iforce),a(iftp),
     +  a(ieps),a(iww),a(iepsd),a(iwks2),a(iforc1),a(itt),a(iforc2),
     +  ifile1,ifint1,ifint2,ldebug,iwrt)
      else
         call ddmpf2(a(iu),a(iw2),a(iwks1),a(iy),a(iforce),a(iftp),
     +  a(ieps),a(iepsd),a(iww),a(if3n),a(iwks2),a(iforc1),a(itt),
     +  a(iforc2),ifile1,ifint1,ifint2,ldebug,iwrt)
      end if
      call ddmpfe(a(iftp),ifile1,ldebug,iwrt)
      return
      end
      subroutine ddmptr(c,iblw,ifw,ifytr,ifdm,ifile1,ifort2,ifint1,
     &  ifint2,maxq)
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical dpres,gpres
c
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/iofile)
c
INCLUDE(common/symtry)
INCLUDE(common/atmblk)
INCLUDE(common/nshel)
      dimension c(*)
c
      cutoff=10.0d0**(-icut)
c
      dpres = .false.
      fpres = .false.
      gpres = .false.
      do i = 1 , nshell
         dpres = dpres .or. (ktype(i).eq.3)
      end do
      do i = 1 , nshell
         fpres = fpres .or. (ktype(i).eq.4)
      end do
      do i = 1 , nshell
         gpres = gpres .or. (ktype(i).eq.5)
      end do
      ntrned = 16
      if (dpres) ntrned = 36
      if (fpres) ntrned = 100
      if (gpres) ntrned = 225
      n2 = ncoorb*ncoorb
      nov = nocca*nvirta
      ntri = ncoorb*(ncoorb+1)/2
      ireq = nov + n2*5 + ntrned*ntri + nw196(5)
      if (ireq.gt.maxq) then
         write (iwr,*) ireq , maxq
         call caserr('not enough core')
      end if
c
      i0 = 1
      i1 = i0 + nw196(5)
      call rdedx(c(i0),nw196(5),ibl196(5),ifild)
c
c     i1=1
      i2 = i1 + ncoorb
      i3 = i2 + n2
      i4 = i3 + n2
      i5 = i4 + nov
      i6 = i5 + n2
      i7 = i6 + n2
      i8 = i7 + n2
      i9 = i8 + n2
c
      m9 = 9
      call secget(isect(9),m9,isec9)
      call rdedx(c(i1),ncoorb,isec9,ifild)
      itype = 0
      call secget(isect(8),itype,ibl)
      ibl = ibl + mvadd
      call rdedx(c(i6),n2,ibl,ifild)
      call ddmps2(c(i0),nshell,
     +   c(i1),c(i2),c(i3),c(i4),c(i5),c(i6),c(i7),c(i8),
     +   c(i9),ifint1,ifint2,ifort2,nocca,nvirta,ncoorb,icount)
      call mpsrt0(ifort2,ifort2,nov,icount,c(i1),c(i1),maxq-nw196(5))
      call rdedx(c(i1),n2,ibl,ifild)
      call ddmps3(c(i1),c(n2+i1),c(n2+n2+i1),ifytr,ifdm,ncoorb,nocca,
     +  ifile1)
      i1 = i0 + nw196(5)
      i2 = i1 + n2
      i3 = i2 + nov
      i4 = i3 + n2
      i5 = i4 + n2
      i6 = i5 + n2
      i7 = i6 + n2
c     itop = i7 + n2
      cut2=10.0d0*cutoff
      call ddmps4(c(i0),nshell,
     +  c(i1),c(i3),c(i2),c(i7),iblw,ifw,ifort2,nocca,
     +  nvirta,ncoorb,c(i4),c(i5),ifytr,ifdm,ntri,c(i6),ifile1,cut2)
c
c
      return
      end
      subroutine ddmpfe(ftp,ifile1,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1   istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension ftp(n3n,n3n)
      nat = n3n/3
      call ddmpti(ftp,nat,nat)
      if (ldebug) write (iwrt,*) 'total so far'
      if (ldebug) call prsqm(ftp,n3n,n3n,n3n,iwrt)
      call wrt3(ftp,n3n*n3n,igradf,ifile1)
      return
      end
      subroutine ddmpf1(y,u,epsd,xt,e,w,w2,force,ftp,ifile1,
     +                  ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1 istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension y(norbs,norbs),u(norbs,norbs,n3n),epsd(norbs,norbs,n3n)
      dimension xt(norbs,norbs),e(norbs),w(norbs,norbs),w2(norbs,norbs)
      dimension force(n3n,n3n),ftp(n3n,n3n)
      call vclr(force,1,n3n*n3n)
      nsq = norbs*norbs
      n3n3 = n3n - 3
      n3nn2 = norbs*norbs*n3n
      if (ldebug) write (iwrt,*) 'new ddmpf1 routine'
      call rdedx(ftp,n3n*n3n,igradf,ifile1)
      call rdedx(epsd,n3nn2,istrme,ifile1)
      call rdedx(y,nsq,istrmy,ifile1)
      call rdedx(u,n3nn2,istrmu,ifile1)
c..........(d(x(,d(y),t**2)  term
c.........(y,u,epsd)  term
c
      do 50 ix = 1 , n3n3
         do 40 iy = 1 , ix
            call vclr(xt,1,nsq)
            call mxmb(u(1,1,ix),norbs,1,epsd(1,1,iy),1,norbs,xt,1,norbs,
     +                norbs,norbs,norbs)
            call mxmb(u(1,1,iy),norbs,1,epsd(1,1,ix),1,norbs,xt,1,norbs,
     +                norbs,norbs,norbs)
            do 30 ir = 1 , norbs
               do 20 is = 1 , norbs
                  force(ix,iy) = force(ix,iy) + (y(ir,is)+y(is,ir))
     +                           *xt(ir,is)
 20            continue
 30         continue
 40      continue
 50   continue
      if (ldebug) write (iwrt,*) ' y,u,epsd term'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      do 70 jk = 1 , n3n3
         do 60 ji = jk , n3n3
            ftp(ji,jk) = ftp(ji,jk) + force(ji,jk)
 60      continue
 70   continue
c
      call rdedx(e,norbs,iste,ifile1)
c
      call vclr(force,1,n3n*n3n)
c
      do 120 ix = 1 , n3n3
         do 110 iy = 1 , ix
            call vclr(xt,1,nsq)
            call mxmb(u(1,1,ix),norbs,1,u(1,1,iy),norbs,1,xt,1,norbs,
     +                norbs,norbs,norbs)
            call mxmb(u(1,1,iy),norbs,1,u(1,1,ix),norbs,1,xt,1,norbs,
     +                norbs,norbs,norbs)
            do 100 ip = 1 , norbs
               do 90 ir = 1 , norbs
                  xt(ir,ip) = xt(ir,ip)*e(ip)
                  do 80 is = 1 , norbs
                     xt(ir,ip) = xt(ir,ip) + u(is,ip,ix)*u(is,ir,iy)
     +                           *e(is)
 80               continue
                  force(ix,iy) = force(ix,iy) - (y(ir,ip)+y(ip,ir))
     +                           *xt(ir,ip)
 90            continue
 100        continue
 110     continue
 120  continue
c
c
c.........(y,u,u,a)   term
      if (ldebug) write (iwrt,*) ' y,u,u,eps term'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      do 140 jk = 1 , n3n3
         do 130 ji = jk , n3n3
            ftp(ji,jk) = ftp(ji,jk) + force(ji,jk)
 130     continue
 140  continue
      call vclr(force,1,n3n*n3n)
      call rdedx(w,nsq,itw,ifile1)
c
      do 180 ix = 1 , n3n3
         do 170 iy = 1 , ix
            call vclr(xt,1,nsq)
            call mxmb(u(1,1,ix),1,norbs,u(1,1,iy),1,norbs,xt,1,norbs,
     +                norbs,norbs,norbs)
            call mxmb(u(1,1,iy),1,norbs,u(1,1,ix),1,norbs,xt,1,norbs,
     +                norbs,norbs,norbs)
            call mxmb(u(1,1,ix),norbs,1,u(1,1,iy),1,norbs,xt,1,norbs,
     +                norbs,norbs,norbs)
            do 160 ir = 1 , norbs
               do 150 ip = 1 , norbs
                  force(ix,iy) = force(ix,iy) + (w(ir,ip)+w(ip,ir))
     +                           *xt(ir,ip)
 150           continue
 160        continue
 170     continue
 180  continue
c
c
      if (ldebug) write (iwrt,*) ' w,u,u term'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      do 200 jk = 1 , n3n3
         do 190 ji = jk , n3n3
            ftp(ji,jk) = ftp(ji,jk) + force(ji,jk)
 190     continue
 200  continue
      call vclr(force,1,n3n*n3n)
      call vclr(xt,1,nsq)
      call rdedx(w2,nsq,itw2,ifile1)
      do 240 ix = 1 , n3n3
         call vclr(xt,1,nsq)
         call mxmb(u(1,1,ix),norbs,1,w2(1,1),1,norbs,xt,1,norbs,norbs,
     +             norbs,norbs)
         do 230 j = 1 , nocc
            do 220 it = 1 , norbs
               if (xt(j,it).ne.0.0d0) then
                  do 210 iy = 1 , ix
                     force(ix,iy) = force(ix,iy) - u(it,j,iy)*xt(j,it)
 210              continue
               end if
 220        continue
 230     continue
 240  continue
      if (ldebug) write (iwrt,*) ' w2,u,u term'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      do 260 jk = 1 , n3n3
         do 250 ji = jk , n3n3
            ftp(ji,jk) = ftp(ji,jk) + force(ji,jk)
 250     continue
 260  continue
      call vclr(force,1,n3n*n3n)
      return
      end
      subroutine ddmpf2(u,tt,a1,a2,force,ftp,e,epsd,ff,f,zint,forc1,t,
     1      forc2,ifile1,ifint1,ifint2,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
      dimension ff(norbs,norbs),a1(norbs,norbs),a2(norbs,norbs)
      dimension f(norbs,norbs,n3n),zint(norbs,norbs),t(nocc*nvirt)
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1   istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension tt(norbs,norbs),forc1(n3n,n3n),forc2(n3n,n3n)
      dimension u(norbs,norbs,n3n),force(n3n,n3n),ftp(n3n,n3n)
      dimension e(norbs),epsd(norbs,norbs,n3n)
      if (ldebug) write (iwrt,*)
     +      'new epsilon and u integral terms (ddmpf2)'
c     nvir = nvirt
      nsq = norbs*norbs
      n3nn2 = norbs*norbs*n3n
      n3n3 = n3n - 3
      call rdedx(e,norbs,iste,ifile1)
      call vclr(force,1,n3n*n3n)
      call vclr(forc1,1,n3n*n3n)
      call vclr(forc2,1,n3n*n3n)
      call search(1,ifint1)
      call search(1,ifint2)
      if (ldebug) write (iwrt,*) 'd*d*t2 term'
c.......................d(x)*d(y)*t2 term
      call rdedx(epsd,n3nn2,istrme,ifile1)
      call rdedx(u,n3nn2,istrmu,ifile1)
      do 310 ir = 1 , norbs
         do 300 is = 1 , ir
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            if (ir.le.nocc) then
c
c
               i = ir
               j = is
               do 30 ia = nocc1 , norbs
                  do 20 ib = nocc1 , norbs
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(ia,ib) = 4.0d0*(a2(ia,ib)+a2(ia,ib)-a2(ib,ia))
     +                           *ddiff
                     zint(ia,ib) = a2(ia,ib)*ddiff
 20               continue
 30            continue
               if (i.eq.j) then
                  do 50 iq1 = 1 , norbs
                     do 40 ip1 = 1 , norbs
                        a2(ip1,iq1) = 0.5d0*a2(ip1,iq1)
 40                  continue
 50               continue
               end if
               call vclr(f,1,n3nn2)
               do 80 ix = 1 , n3n3
                  call mxmb(zint(nocc1,nocc1),norbs,1,
     +                      epsd(nocc1,nocc1,ix),norbs,1,
     +                      f(nocc1,nocc1,ix),norbs,1,nvirt,nvirt,nvirt)
                  if (i.ne.j) then
                     call mxmb(zint(nocc1,nocc1),1,norbs,
     +                         epsd(nocc1,nocc1,ix),1,norbs,
     +                         f(nocc1,nocc1,ix),1,norbs,nvirt,nvirt,
     +                         nvirt)
                  end if
                  do 70 ib = nocc1 , norbs
                     do 60 ic = nocc1 , norbs
                        f(ic,ib,ix) = f(ic,ib,ix)
     +                                /(e(ic)+e(ib)-e(i)-e(j))
 60                  continue
 70               continue
 80            continue
               do 90 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt(nocc1,nocc1),1,norbs,epsd(nocc1,nocc1,iy)
     +                      ,1,norbs,ff(nocc1,nocc1),1,norbs,nvirt,
     +                      nvirt,nvirt)
                  call mxmb(tt(nocc1,nocc1),norbs,1,epsd(nocc1,nocc1,iy)
     +                      ,1,norbs,ff(nocc1,nocc1),norbs,1,nvirt,
     +                      nvirt,nvirt)
                  call mxmb(f,nsq,1,ff,1,nsq,force(iy,1),n3n,1,iy,nsq,1)
 90            continue
c
               call vclr(f,1,n3nn2)
               do 100 ix = 1 , n3n3
                  call mxmb(u(1,nocc1,ix),1,norbs,tt(nocc1,nocc1),1,
     +                      norbs,f(nocc1,1,ix),norbs,1,norbs,nvirt,
     +                      nvirt)
 100           continue
               do 110 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(u(1,nocc1,iy),norbs,1,a2,norbs,1,ff(nocc1,1)
     +                      ,1,norbs,nvirt,norbs,norbs)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
c
c
c
 110           continue
            else if (is.gt.nocc) then
               ia = ir
               ib = is
c
c
               do 130 i = 1 , nocc
                  do 120 j = 1 , nocc
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,j) = 4.0d0*(a2(i,j)+a2(i,j)-a2(j,i))*ddiff
                     zint(i,j) = a2(i,j)*ddiff
 120              continue
 130           continue
               if (ia.eq.ib) then
                  do 150 ip1 = 1 , norbs
                     do 140 iq1 = 1 , norbs
                        a2(ip1,iq1) = 0.5d0*a2(ip1,iq1)
 140                 continue
 150              continue
               end if
               call vclr(f,1,n3nn2)
               do 180 ix = 1 , n3n3
                  call mxmb(zint(1,1),norbs,1,epsd(1,1,ix),norbs,1,
     +                      f(1,1,ix),norbs,1,nocc,nocc,nocc)
                  if (ia.ne.ib) then
                     call mxmb(zint(1,1),1,norbs,epsd(1,1,ix),1,norbs,
     +                         f(1,1,ix),1,norbs,nocc,nocc,nocc)
                  end if
                  do 170 j = 1 , nocc
                     do 160 k = 1 , nocc
                        f(k,j,ix) = f(k,j,ix)/(e(ia)+e(ib)-e(k)-e(j))
 160                 continue
 170              continue
 180           continue
               do 190 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt,1,norbs,epsd(1,1,iy),1,norbs,ff,1,norbs,
     +                      nocc,nocc,nocc)
                  call mxmb(tt,norbs,1,epsd(1,1,iy),1,norbs,ff,norbs,1,
     +                      nocc,nocc,nocc)
                  call mxmb(f,nsq,1,ff,1,nsq,force(iy,1),n3n,1,iy,nsq,1)
 190           continue
c
               call vclr(f,1,n3nn2)
               do 200 ix = 1 , n3n3
                  call mxmb(u(1,1,ix),1,norbs,tt,1,norbs,f(1,1,ix),
     +                      norbs,1,norbs,nocc,nocc)
 200           continue
               do 210 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(a2,1,norbs,u(1,1,iy),1,norbs,ff,norbs,1,
     +                      norbs,norbs,nocc)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
 210           continue
            else
               ia = ir
               j = is
               ibi = 0
               do 230 ib = nocc1 , norbs
                  do 220 i = 1 , nocc
                     ibi = ibi + 1
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,ib) = 4.0d0*(a2(i,ib)+a2(i,ib)-a1(i,ib))*ddiff
                     t(ibi) = (a1(i,ib)+a1(i,ib)-a2(i,ib))*ddiff
                     tt(ib,i) = 4.0d0*t(ibi)
                     zint(i,ib) = a2(i,ib)*ddiff
                     zint(ib,i) = a1(ib,i)*ddiff
 220              continue
 230           continue
c
c
               call vclr(f,1,n3nn2)
               do 260 ix = 1 , n3n3
                  call mxmb(zint(1,nocc1),1,norbs,epsd(nocc1,nocc1,ix),
     +                      1,norbs,f(1,nocc1,ix),1,norbs,nocc,nvirt,
     +                      nvirt)
                  call mxmb(zint(nocc1,1),norbs,1,epsd(nocc1,nocc1,ix),
     +                      1,norbs,f(nocc1,1,ix),norbs,1,nocc,nvirt,
     +                      nvirt)
                  do 250 i = 1 , nocc
                     do 240 id = nocc + 1 , norbs
                        ddiff = 1.0d0/(e(ia)+e(id)-e(i)-e(j))
                        f(i,id,ix) = ddiff*f(i,id,ix)
                        f(id,i,ix) = ddiff*f(id,i,ix)
 240                 continue
 250              continue
 260           continue
               do 270 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt(1,nocc1),norbs,1,epsd(1,1,iy),1,norbs,
     +                      ff(1,nocc1),norbs,1,nvirt,nocc,nocc)
                  call mxmb(tt(nocc1,1),1,norbs,epsd(1,1,iy),1,norbs,
     +                      ff(nocc1,1),1,norbs,nvirt,nocc,nocc)
                  call mxmb(f,nsq,1,ff,1,nsq,forc2(iy,1),n3n,1,n3n3,nsq,
     +                      1)
 270           continue
c
c
               call vclr(f,1,n3nn2)
               do 280 ix = 1 , n3n3
                  call mxmb(u(1,nocc1,ix),1,norbs,tt(1,nocc1),norbs,1,
     +                      f(1,1,ix),norbs,1,norbs,nvirt,nocc)
                  call mxmb(u(1,1,ix),1,norbs,tt(nocc1,1),norbs,1,
     +                      f(nocc1,1,ix),norbs,1,norbs,nocc,nvirt)
 280           continue
               do 290 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(u(1,1,iy),norbs,1,a2,1,norbs,ff,1,norbs,
     +                      nocc,norbs,norbs)
                  call mxmb(u(1,nocc1,iy),norbs,1,a1,norbs,1,ff(nocc1,1)
     +                      ,1,norbs,nvirt,norbs,norbs)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
c
c
c
 290           continue
            end if
c
c
c
 300     continue
 310  continue
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc1,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc2,n3n,n3n,n3n,iwrt)
      do 330 jk = 1 , n3n
         do 320 ji = jk , n3n
            ftp(ji,jk) = ftp(ji,jk) - force(ji,jk) - forc1(ji,jk)
     +                   - forc1(jk,ji) + forc2(ji,jk) + forc2(jk,ji)
 320     continue
 330  continue
      return
      end
      subroutine ddmpf3(epsd,tt,a1,a2,force,ftp,e,ff,f,zint,forc1,t,
     &  forc2,ifile1,ifint1,ifint2,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
      dimension ff(norbs,norbs),a1(norbs,norbs),a2(norbs,norbs)
      dimension f(norbs,norbs,n3n),zint(norbs,norbs),t(nocc*nvirt)
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1   istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension tt(norbs,norbs),forc1(n3n,n3n),forc2(n3n,n3n)
      dimension epsd(norbs,norbs,n3n),force(n3n,n3n),ftp(n3n,n3n)
      dimension e(norbs)
      if (ldebug) write (iwrt,*)
     + '**********ddmpf3   pass   1  ***********'
c     nvir = nvirt
      nsq = norbs*norbs
      n3nn2 = norbs*norbs*n3n
      n3n3 = n3n - 3
      call rdedx(e,norbs,iste,ifile1)
      call vclr(force,1,n3n*n3n)
      call vclr(forc1,1,n3n*n3n)
      call vclr(forc2,1,n3n*n3n)
      call search(1,ifint1)
      call search(1,ifint2)
      if (ldebug) write (iwrt,*) 'd*d*t2 term'
c.......................d(x)*d(y)*t2 term
      call rdedx(epsd,n3nn2,istrme,ifile1)
      do 210 ir = 1 , norbs
         do 200 is = 1 , ir
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            if (ir.le.nocc) then
c
c
               i = ir
               j = is
               do 30 ia = nocc1 , norbs
                  do 20 ib = nocc1 , norbs
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(ia,ib) = 4.0d0*(a2(ia,ib)+a2(ia,ib)-a2(ib,ia))
     +                           *ddiff
                     zint(ia,ib) = a2(ia,ib)*ddiff
 20               continue
 30            continue
               call vclr(f,1,n3nn2)
               do 60 ix = 1 , n3n3
                  call mxmb(zint(nocc1,nocc1),norbs,1,
     +                      epsd(nocc1,nocc1,ix),norbs,1,
     +                      f(nocc1,nocc1,ix),norbs,1,nvirt,nvirt,nvirt)
                  if (i.ne.j) then
                     call mxmb(zint(nocc1,nocc1),1,norbs,
     +                         epsd(nocc1,nocc1,ix),1,norbs,
     +                         f(nocc1,nocc1,ix),1,norbs,nvirt,nvirt,
     +                         nvirt)
                  end if
                  do 50 ic = nocc1 , norbs
                     do 40 ib = nocc1 , norbs
                        f(ic,ib,ix) = f(ic,ib,ix)
     +                                /(e(ic)+e(ib)-e(i)-e(j))
 40                  continue
 50               continue
 60            continue
               do 70 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt(nocc1,nocc1),1,norbs,epsd(nocc1,nocc1,iy)
     +                      ,1,norbs,ff(nocc1,nocc1),1,norbs,nvirt,
     +                      nvirt,nvirt)
                  call mxmb(tt(nocc1,nocc1),norbs,1,epsd(nocc1,nocc1,iy)
     +                      ,1,norbs,ff(nocc1,nocc1),norbs,1,nvirt,
     +                      nvirt,nvirt)
                  call mxmb(f,nsq,1,ff,1,nsq,force(iy,1),n3n,1,iy,nsq,1)
c
c
 70            continue
            else if (is.gt.nocc) then
               ia = ir
               ib = is
c
               do 90 i = 1 , nocc
                  do 80 j = 1 , nocc
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,j) = 4.0d0*(a2(i,j)+a2(i,j)-a2(j,i))*ddiff
                     zint(i,j) = a2(i,j)*ddiff
 80               continue
 90            continue
               call vclr(f,1,n3nn2)
               do 120 ix = 1 , n3n3
                  call mxmb(zint(1,1),norbs,1,epsd(1,1,ix),norbs,1,
     +                      f(1,1,ix),norbs,1,nocc,nocc,nocc)
                  if (ia.ne.ib) then
                     call mxmb(zint(1,1),1,norbs,epsd(1,1,ix),1,norbs,
     +                         f(1,1,ix),1,norbs,nocc,nocc,nocc)
                  end if
                  do 110 k = 1 , nocc
                     do 100 j = 1 , nocc
                        f(k,j,ix) = f(k,j,ix)/(e(ia)+e(ib)-e(k)-e(j))
 100                 continue
 110              continue
 120           continue
               do 130 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt,1,norbs,epsd(1,1,iy),1,norbs,ff,1,norbs,
     +                      nocc,nocc,nocc)
                  call mxmb(tt,norbs,1,epsd(1,1,iy),1,norbs,ff,norbs,1,
     +                      nocc,nocc,nocc)
                  call mxmb(f,nsq,1,ff,1,nsq,force(iy,1),n3n,1,iy,nsq,1)
 130           continue
            else
               ia = ir
               j = is
               ibi = 0
               do 150 ib = nocc1 , norbs
                  do 140 i = 1 , nocc
                     ibi = ibi + 1
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,ib) = 4.0d0*(a2(i,ib)+a2(i,ib)-a1(i,ib))*ddiff
                     t(ibi) = (a1(i,ib)+a1(i,ib)-a2(i,ib))*ddiff
                     tt(ib,i) = 4.0d0*t(ibi)
                     zint(i,ib) = a2(i,ib)*ddiff
                     zint(ib,i) = a1(ib,i)*ddiff
 140              continue
 150           continue
c
c
               call vclr(f,1,n3nn2)
               do 180 ix = 1 , n3n3
                  call mxmb(zint(1,nocc1),1,norbs,epsd(nocc1,nocc1,ix),
     +                      1,norbs,f(1,nocc1,ix),1,norbs,nocc,nvirt,
     +                      nvirt)
                  call mxmb(zint(nocc1,1),norbs,1,epsd(nocc1,nocc1,ix),
     +                      1,norbs,f(nocc1,1,ix),norbs,1,nocc,nvirt,
     +                      nvirt)
                  do 170 i = 1 , nocc
                     do 160 id = nocc1 , norbs
                        ddiff = 1.0d0/(e(ia)+e(id)-e(i)-e(j))
                        f(i,id,ix) = ddiff*f(i,id,ix)
                        f(id,i,ix) = ddiff*f(id,i,ix)
 160                 continue
 170              continue
 180           continue
               do 190 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(tt(1,nocc1),norbs,1,epsd(1,1,iy),1,norbs,
     +                      ff(1,nocc1),norbs,1,nvirt,nocc,nocc)
                  call mxmb(tt(nocc1,1),1,norbs,epsd(1,1,iy),1,norbs,
     +                      ff(nocc1,1),1,norbs,nvirt,nocc,nocc)
                  call mxmb(f,nsq,1,ff,1,nsq,forc2(iy,1),n3n,1,n3n3,nsq,
     +                      1)
c
c
 190           continue
            end if
c
c
 200     continue
 210  continue
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc1,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc2,n3n,n3n,n3n,iwrt)
      do 230 ji = 1 , n3n
         do 220 jk = 1 , ji
            ftp(ji,jk) = ftp(ji,jk) - force(ji,jk) - forc1(ji,jk)
     +                   - forc1(jk,ji) + forc2(ji,jk) + forc2(jk,ji)
 220     continue
 230  continue
      return
      end
      subroutine ddmpf4(u,tt,a1,a2,force,ftp,e,ff,f,zint,forc1,t,
     &   forc2,ifile1,ifint1,ifint2,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
      logical ldebug
      dimension ff(norbs,norbs),a1(norbs,norbs),a2(norbs,norbs)
      dimension f(norbs,norbs,n3n),zint(norbs,norbs),t(nocc*nvirt)
      common/parm/norbs,n3n,nvirt,nocc,nab,nij,nocc1,itw,itw2,
     1   istrmy,istrmu,istrme,istrmt,igradf,ist2,iste
      dimension tt(norbs,norbs),forc1(n3n,n3n),forc2(n3n,n3n)
      dimension force(n3n,n3n),ftp(n3n,n3n)
      dimension e(norbs),u(norbs,norbs,n3n)
      if (ldebug) write (iwrt,*)
     + '**********ddmpf4   pass   2  ***********'
c     nvir = nvirt
      nsq = norbs*norbs
      n3nn2 = norbs*norbs*n3n
      n3n3 = n3n - 3
      call rdedx(e,norbs,iste,ifile1)
      call vclr(force,1,n3n*n3n)
      call vclr(forc1,1,n3n*n3n)
      call vclr(forc2,1,n3n*n3n)
      call search(1,ifint1)
      call search(1,ifint2)
      if (ldebug) write (iwrt,*) 'd*d*t2 term'
c.......................d(x)*d(y)*t2 term
      call rdedx(u,n3nn2,istrmu,ifile1)
      do 190 ir = 1 , norbs
         do 180 is = 1 , ir
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            if (ir.le.nocc) then
c
               i = ir
               j = is
               do 30 ia = nocc1 , norbs
                  do 20 ib = nocc1 , norbs
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(ia,ib) = 4.0d0*(a2(ia,ib)+a2(ia,ib)-a2(ib,ia))
     +                           *ddiff
                     zint(ia,ib) = a2(ia,ib)*ddiff
 20               continue
 30            continue
               if (i.eq.j) then
                  do 50 ip1 = 1 , norbs
                     do 40 iq1 = 1 , norbs
                        a2(ip1,iq1) = 0.5d0*a2(ip1,iq1)
 40                  continue
 50               continue
               end if
c
               call vclr(f,1,n3nn2)
               do 60 ix = 1 , n3n3
                  call mxmb(u(1,nocc1,ix),1,norbs,tt(nocc1,nocc1),1,
     +                      norbs,f(nocc1,1,ix),norbs,1,norbs,nvirt,
     +                      nvirt)
 60            continue
               do 70 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(u(1,nocc1,iy),norbs,1,a2,norbs,1,ff(nocc1,1)
     +                      ,1,norbs,nvirt,norbs,norbs)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
c
c
 70            continue
            else if (is.gt.nocc) then
               ia = ir
               ib = is
c
c
               do 90 i = 1 , nocc
                  do 80 j = 1 , nocc
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,j) = 4.0d0*(a2(i,j)+a2(i,j)-a2(j,i))*ddiff
                     zint(i,j) = a2(i,j)*ddiff
 80               continue
 90            continue
               if (ia.eq.ib) then
                  do 110 ip1 = 1 , norbs
                     do 100 iq1 = 1 , norbs
                        a2(ip1,iq1) = 0.5d0*a2(ip1,iq1)
 100                 continue
 110              continue
               end if
c
               call vclr(f,1,n3nn2)
               do 120 ix = 1 , n3n3
                  call mxmb(u(1,1,ix),1,norbs,tt,1,norbs,f(1,1,ix),
     +                      norbs,1,norbs,nocc,nocc)
 120           continue
               do 130 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(u(1,1,iy),norbs,1,a2,norbs,1,ff,1,norbs,
     +                      nocc,norbs,norbs)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
 130           continue
            else
               ia = ir
               j = is
               ibi = 0
               do 150 ib = nocc1 , norbs
                  do 140 i = 1 , nocc
                     ibi = ibi + 1
                     ddiff = 1.0d0/(e(ia)+e(ib)-e(i)-e(j))
                     tt(i,ib) = 4.0d0*(a2(i,ib)+a2(i,ib)-a1(i,ib))*ddiff
                     t(ibi) = (a1(i,ib)+a1(i,ib)-a2(i,ib))*ddiff
                     tt(ib,i) = 4.0d0*t(ibi)
                     zint(i,ib) = a2(i,ib)*ddiff
                     zint(ib,i) = a1(ib,i)*ddiff
 140              continue
 150           continue
c
c
               call vclr(f,1,n3nn2)
               do 160 ix = 1 , n3n3
                  call mxmb(u(1,nocc1,ix),1,norbs,tt(1,nocc1),norbs,1,
     +                      f(1,1,ix),norbs,1,norbs,nvirt,nocc)
                  call mxmb(u(1,1,ix),1,norbs,tt(nocc1,1),norbs,1,
     +                      f(nocc1,1,ix),norbs,1,norbs,nocc,nvirt)
 160           continue
               do 170 iy = 1 , n3n3
                  call vclr(ff,1,nsq)
                  call mxmb(u(1,1,iy),norbs,1,a2,1,norbs,ff,1,norbs,
     +                      nocc,norbs,norbs)
                  call mxmb(u(1,nocc1,iy),norbs,1,a1,norbs,1,ff(nocc1,1)
     +                      ,1,norbs,nvirt,norbs,norbs)
                  call mxmb(f,nsq,1,ff,1,nsq,forc1(1,iy),1,n3n,n3n3,nsq,
     +                      1)
c
 170           continue
            end if
c
 180     continue
 190  continue
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc1,n3n,n3n,n3n,iwrt)
      if (ldebug) call prsqm(forc2,n3n,n3n,n3n,iwrt)
      do 210 ji = 1 , n3n
         do 200 jk = 1 , ji
            ftp(ji,jk) = ftp(ji,jk) - force(ji,jk) - forc1(ji,jk)
     +                   - forc1(jk,ji) + forc2(ji,jk) + forc2(jk,ji)
 200     continue
 210  continue
      return
      end
_EXTRACT(mpmakw,mpmakw,pclinux)
      subroutine mpmakw(y,z,w,w2,ibly,iblz,iblw,iblw2,a1,a2,nocca,
     1   nvirta,ncoorb,ifils,e,ifint1,ifint2,ifile1)
c
      implicit REAL  (a-h,o-z)
      dimension y(ncoorb,ncoorb),w(ncoorb,ncoorb),
     1  z(nocca*nvirta),a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),
     2 e(ncoorb),w2(ncoorb,ncoorb)
c
      nsq = ncoorb*ncoorb
      call vclr(w2,1,nsq)
c
c
      call search(1,ifint1)
      call search(1,ifint2)
      call rdedx(z,nocca*nvirta,iblz,ifils)
      call rdedx(y,nsq,ibly,ifile1)
      do 30 i = 1 , nocca
         do 20 ia = nocca + 1 , ncoorb
            iai = (ia-nocca-1)*nocca + i
            y(ia,i) = -z(iai)
 20      continue
 30   continue
      call wrt3(y,nsq,ibly,ifile1)
c
c
      call rdedx(a1,nsq,iblw,ifile1)
      do 50 is = 1 , ncoorb
         do 40 ir = 1 , ncoorb
            w(ir,is) = y(ir,is)*e(is)
 40      continue
 50   continue
      do 70 j = 1 , nocca
         do 60 i = 1 , nocca
            w(i,j) = w(i,j) - 0.5d0*a1(i,j)
 60      continue
 70   continue
      do 90 j = nocca + 1 , ncoorb
         do 80 i = nocca + 1 , ncoorb
            w(i,j) = w(i,j) - 0.5d0*a1(i,j)
 80      continue
 90   continue
      do 110 ia = nocca + 1 , ncoorb
         do 100 i = 1 , nocca
            w(i,ia) = w(i,ia) - a1(i,ia)
 100     continue
 110  continue
c
c
      do 170 ip = 1 , ncoorb
         do 160 iq = 1 , ip
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            yf = (y(ip,iq)+y(iq,ip))
            if (ip.eq.iq) yf = 0.5d0*yf
            yf1 = 0.5d0*yf
            yf2 = -yf
            if (yf.ne.0.0d0) then
               do 130 k = 1 , nocca
                  do 120 j = 1 , nocca
                     w(j,k) = w(j,k)
     +                        + yf1*(a1(j,k)*4.0d0-a2(j,k)-a2(k,j))
 120              continue
 130           continue
               do 150 k = 1 , ncoorb
                  do 140 j = 1 , ncoorb
                     w2(j,k) = w2(j,k)
     +                         + yf2*(a1(j,k)*4.0d0-a2(j,k)-a2(k,j))
 140              continue
 150           continue
            end if
 160     continue
 170  continue
c
      call wrt3(w,nsq,iblw,ifile1)
      call wrt3(w2,nsq,iblw2,ifile1)
      return
      end
_ENDEXTRACT
      subroutine mpmaky(a1,a2,e,ncoorb,y,zlg,
     1 nocca,istrmy,nvirt,t,
     2 iblw,iblks,ifils,b,mn,ifile1,istrma,ifint1,ifint2)
c
      implicit REAL  (a-h,o-z)
      dimension a1(ncoorb,ncoorb),a2(ncoorb,ncoorb),e(ncoorb),
     + y(ncoorb,ncoorb),t(ncoorb,ncoorb),zlg(ncoorb,ncoorb),b(mn)
c
c
      nov = nocca*nvirt
      nocc = nocca
      norbs = ncoorb
      nsq = norbs*norbs
      nocc1 = nocc + 1
c     nocc2 = nocc + 2
      call vclr(y,1,norbs*norbs)
      call vclr(zlg,1,norbs*norbs)
c
      call search(1,ifint1)
      call search(1,ifint2)
      call search(1,istrma)
      do 70 ip = 1 , norbs
         do 60 iq = 1 , ip
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
c
            if (iq.le.nocca) then
               if (ip.gt.nocca) then
                  call vclr(t,1,norbs*norbs)
                  ebj = e(ip) - e(iq)
c
                  ibi = 0
                  do 30 ia = nocc1 , norbs
                     do 20 i = 1 , nocc
                        ibi = ibi + 1
                        b(ibi) = (a1(i,ia)+a1(i,ia)-a2(i,ia))
     +                           /(ebj+e(ia)-e(i))
                        t(i,ia) = 4.0d0*b(ibi)
 20                  continue
 30               continue
                  call wtedz(b,nov,istrma)
c
                  call mxmb(t(1,nocc1),1,norbs,a1(nocc1,1),1,norbs,zlg,
     +                      norbs,1,nocc,nvirt,norbs)
                  call mxmb(a1,1,norbs,t(1,nocc1),1,norbs,zlg(1,nocc1),
     +                      1,norbs,norbs,nocc,nvirt)
c
                  do 50 ia = nocc1 , norbs
                     do 40 i = 1 , nocc
                        a1(i,ia) = 0.5d0*a1(i,ia)/(ebj+e(ia)-e(i))
 40                  continue
 50               continue
c
                  call mxmb(t(1,nocc1),norbs,1,a1(1,nocc1),1,norbs,
     +                      y(nocc1,nocc1),1,norbs,nvirt,nocc,nvirt)
                  call mxmb(t(1,nocc1),1,norbs,a1(1,nocc1),norbs,1,y,1,
     +                      norbs,nocc,nvirt,nocc)
               end if
            end if
 60      continue
 70   continue
      do 90 ib = nocc1 , norbs
         do 80 ia = nocc1 , norbs
            y(ia,ib) = -y(ia,ib)
 80      continue
 90   continue
c
c.....................integral contribution complete
c..............................to calculate t**2 contributions
c..............................store partial sums j,a,b and i,j,b of t**
c
      do 110 ib = nocc1 , norbs
         do 100 ia = nocc1 , norbs
            zlg(ia,ib) = zlg(ia,ib) + y(ia,ib)*(e(ia)-e(ib))
 100     continue
 110  continue
      do 130 j = 1 , nocca
         do 120 i = 1 , nocca
            zlg(i,j) = zlg(i,j) + y(i,j)*(e(i)-e(j))
 120     continue
 130  continue
c
      call search(1,ifint2)
      call search(1,ifint1)
      do 170 ip = 1 , norbs
         do 160 iq = 1 , ip
            call rdedz(a1,nsq,ifint1)
            call rdedz(a2,nsq,ifint2)
            if (iq.le.nocca) then
               if (ip.gt.nocca) then
                  do 150 j = 1 , norbs
                     do 140 k = 1 , norbs
                        zlg(ip,iq) = zlg(ip,iq) + y(j,k)
     +                               *(a1(j,k)*4.0d0-a2(j,k)-a2(k,j))
 140                 continue
 150              continue
               end if
            end if
 160     continue
 170  continue
c
      do 190 i = nocc1 , norbs
         do 180 j = 1 , nocc
            kt = (i-nocca-1)*nocca + j
            b(kt) = zlg(i,j) - zlg(j,i)
 180     continue
 190  continue
c........................b  is  now  calculated
      do 210 j = 1 , norbs
         do 200 i = 1 , norbs
            y(i,j) = -y(i,j)
 200     continue
 210  continue
      call wrt3(y,ncoorb*ncoorb,istrmy,ifile1)
      call wrt3(zlg,ncoorb*ncoorb,iblw,ifile1)
      call wrt3(b,mn,iblks,ifils)
      return
      end
      subroutine wrtskp(nat3,iwr)
      implicit REAL  (a-h,o-z)
      logical skip
INCLUDE(common/sizes)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3),nuniq
c
      do 20 i = 1 , nat3
         if (skip(i)) write (iwr,6010) i
 20   continue
      return
 6010 format(1x,'perturbation', i4,' will be skipped')
      end
      subroutine rdsrtd(a,ia,ifort,n3n,nocca,nvirta,cut)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical skip
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk1,mkk1,g(5119)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,ncols,nrows,nbuck
      common/maxlen/maxq
      common/symmos/imos(8,maxorb),ibass(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
      dimension ia(*),a(*)
      data lastb/999999/
c
      nov = nocca*nvirta
c
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i) + nov
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      call vclr(g,1,nsz340+nsz170)
c
      call search(1,ifort)
      iblock = 0
c      ninb = no of elements ib bucket(coreload)
      ninb = nadd*ncols
c
      do 50 ipr = 1 , n3n
c
         if (.not.(skip(ipr))) then
c
            do 40 ibj = 1 , nov
               call rdedz(a(1),nov,ifort)
               do 30 iai = 1 , nov
                  val = a(iai)
                  if (dabs(val).ge.cut) then
                     iipr = mapun(ipr)
                     iaddr = (iai-1)*ncols + (iipr-1)*nov + ibj
                     ibuck = (iaddr-1)/ninb
                     iaddr = iaddr - ninb*ibuck
                     ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                     nwb = nwbuck(ibuck) + 1
                     a(ibase(ibuck)+nwb) = val
                     ia(ibasen(ibuck)+nwb) = iaddr
                     nwbuck(ibuck) = nwb
                     if (nwb.eq.ibl5) then
c
c     this block full - empty
c
                        call stopbk
                        mkk1 = mark(ibuck)
                        nkk1 = nwb
_IFN1(iv)                     call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)                     call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)                       call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)                       call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
                        call sttout
                        nwbuck(ibuck) = 0
                        mark(ibuck) = iblock
                        iblock = iblock + nsz
                     end if
                  end if
 30            continue
 40         continue
         end if
 50   continue
c
c      empty anything remaining in buckets
c
      do 60 ibuck = 1 , nbuck
         nwb = nwbuck(ibuck)
         if (nwb.ne.0) then
            call stopbk
            mkk1 = mark(ibuck)
            nkk1 = nwb
_IFN1(iv)            call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)            call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)          call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)          call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
            call sttout
            nwbuck(ibuck) = 0
            mark(ibuck) = iblock
            iblock = iblock + nsz
         end if
 60   continue
c
c
      call stopbk
      return
      end
      subroutine ddmprs(force,ftp,u,f,ibforc,ibf,ibu,n3n,ncoorb,ifile1
     &  ,ldebug,iwrt)
      implicit REAL  (a-h,o-z)
c
      logical ldebug
      dimension force(n3n,n3n),ftp(n3n,n3n),u(ncoorb*ncoorb*n3n)
     1  ,f(ncoorb*ncoorb*n3n)
c
c
c     nat = n3n/3
      n2 = ncoorb*ncoorb
      call vclr(force,1,n3n*n3n)
      call rdedx(f,n3n*n2,ibf,ifile1)
c
      call ddrssy(f,u,ncoorb,n3n)
      call rdedx(u,n3n*n2,ibu,ifile1)
c
      call mxmb(f,n2,1,u,1,n2,force,1,n3n,n3n,n2,n3n)
c
      do 30 ix = 1 , n3n
         do 20 iy = 1 , ix
            force(ix,iy) = -4.0d0*(force(ix,iy)+force(iy,ix))
 20      continue
 30   continue
c
      if (ldebug) write (iwrt,*) 'new rs derint term'
      if (ldebug) call prsqm(force,n3n,n3n,n3n,iwrt)
      call rdedx(ftp,n3n*n3n,ibforc,ifile1)
      do 50 ix = 1 , n3n
         do 40 iy = 1 , ix
            ftp(ix,iy) = ftp(ix,iy) + force(ix,iy)
 40      continue
 50   continue
      call wrt3(ftp,n3n*n3n,ibforc,ifile1)
c
c
      return
      end
      subroutine ddrssy(f,u,ncoorb,n3n)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
c
c     symmetrisation for routine rsterm
c
INCLUDE(common/symtry)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      dimension f(ncoorb,ncoorb,n3n),u(ncoorb,ncoorb,n3n)
      logical skip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3),nuniq,nuni1,nuni2,nuni3
c
      nat = n3n/3
      if (nt.ne.1) then
         an = 1.0d0/dfloat(nt)
c
         do 80 n = 1 , nat
            do 70 nc = 1 , 3
               nprt = (n-1)*3 + nc
               do 30 jj = 1 , ncoorb
                  do 20 ii = 1 , ncoorb
                     u(ii,jj,nprt) = f(ii,jj,nprt)
 20               continue
 30            continue
c
               do 60 iop = 2 , nt
                  nuprt = mptr(iop,nprt)
                  nnuprt = iabs(nuprt)
                  isigny = nuprt/nnuprt
                  signy = dfloat(isigny)
                  do 50 iiq = 1 , ncoorb
                     do 40 iip = 1 , ncoorb
                        isign = imos(iop,iip)*imos(iop,iiq)
                        sign = dfloat(isign)
                        u(iip,iiq,nprt) = u(iip,iiq,nprt)
     +                     + sign*signy*f(iip,iiq,nnuprt)
 40                  continue
 50               continue
 60            continue
 70         continue
 80      continue
c
c
         do 110 kk = 1 , n3n
            do 100 jj = 1 , ncoorb
               do 90 ii = 1 , ncoorb
                  f(ii,jj,kk) = an*u(ii,jj,kk)
 90            continue
 100        continue
 110     continue
      end if
      if (nuniq.eq.0) return
      do 130 jj = 1 , ncoorb
         do 120 ii = 1 , ncoorb
            f(ii,jj,nuni1) = 0.0d0
            f(ii,jj,nuni2) = 0.0d0
            f(ii,jj,nuni3) = 0.0d0
 120     continue
 130  continue
      nat = n3n/3
      do 160 n = 1 , nat
         if (n.ne.nuniq) then
            do 150 jj = 1 , ncoorb
               do 140 ii = 1 , ncoorb
                  f(ii,jj,nuni1) = f(ii,jj,nuni1) - f(ii,jj,n*3-2)
                  f(ii,jj,nuni2) = f(ii,jj,nuni2) - f(ii,jj,n*3-1)
                  f(ii,jj,nuni3) = f(ii,jj,nuni3) - f(ii,jj,n*3)
 140           continue
 150        continue
         end if
 160  continue
      return
      end
      subroutine ddmps0(a,w,nov,n3n,nops,ia,i,nocca,ncoorb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension a(nov,n3n),w(nov,n3n)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      logical skip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3),nuniq,nuni1,nuni2,nuni3
c
      call vclr(w,1,nov*n3n)
      if (nops.ne.1) then
c        an = 1.0d0/dfloat(nops)
c
         do 50 ixu = 1 , nunpr
            ix = mapnu(ixu)
            do 40 iop = 1 , nops
               nupr = mptr(iop,ix)
               nnupr = iabs(nupr)
               if (skip(nnupr)) then
                  isigai = imos(iop,ia)*imos(iop,i)
                  isigpr = nupr/nnupr
                  sign = dfloat(isigai*isigpr)
                  ibj = 0
                  do 30 ib = nocca + 1 , ncoorb
                     signb = sign*dfloat(imos(iop,ib))
                     do 20 j = 1 , nocca
                        ibj = ibj + 1
                        w(ibj,nnupr) = signb*dfloat(imos(iop,j))
     +                                 *a(ibj,ixu)
 20                  continue
 30               continue
               end if
 40         continue
 50      continue
c
c
         do 70 ixu = 1 , nunpr
            ix = mapnu(ixu)
            do 60 ibj = 1 , nov
               w(ibj,ix) = a(ibj,ixu)
 60         continue
 70      continue
         call dcopy(nov*n3n,w,1,a,1)
      end if
      if (nuniq.eq.0) return
      do 80 ii = 1 , nov
         a(ii,nuni1) = 0.0d0
         a(ii,nuni2) = 0.0d0
         a(ii,nuni3) = 0.0d0
 80   continue
      nat = n3n/3
      do 100 n = 1 , nat
         if (n.ne.nuniq) then
            do 90 ii = 1 , nov
               a(ii,nuni1) = a(ii,nuni1) - a(ii,n*3-2)
               a(ii,nuni2) = a(ii,nuni2) - a(ii,n*3-1)
               a(ii,nuni3) = a(ii,nuni3) - a(ii,n*3)
 90         continue
         end if
 100  continue
      return
      end
      subroutine ddmps1(a,w,nov,n3n,nops,ia,i,nocca,ncoorb)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension a(nov,n3n),w(nov,n3n)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      logical skip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
c
      if (nops.eq.1) return
      an = 1.0d0/dfloat(nops)
      call vclr(w,1,nov*n3n)
c
      do 50 ixu = 1 , nunpr
         ix = mapnu(ixu)
         do 40 iop = 1 , nops
            nupr = mptr(iop,ix)
            nnupr = iabs(nupr)
            if (.not.(skip(nnupr))) then
               nnnupr = mapun(nnupr)
               isigai = imos(iop,ia)*imos(iop,i)
               isigpr = nupr/nnupr
               sign = dfloat(isigai*isigpr)
               ibj = 0
               do 30 ib = nocca + 1 , ncoorb
                  signb = sign*dfloat(imos(iop,ib))
                  do 20 j = 1 , nocca
                     ibj = ibj + 1
                     w(ibj,nnnupr) = w(ibj,nnnupr)
     +                               + signb*dfloat(imos(iop,j))
     +                               *a(ibj,ixu)
 20               continue
 30            continue
            end if
 40      continue
 50   continue
c
      do 70 jj = 1 , nunpr
         do 60 ii = 1 , nov
            a(ii,jj) = an*w(ii,jj)
 60      continue
 70   continue
c
      return
      end
      subroutine mpsrtd(q,iq,nocca,nvirta,ifortr,ifortw,n3n,
     +                  iwrt,cutoff)
c
c--------------------------------------------------------------
c
c     part4
c     full derivatives of integrals and final assembly
c
c--------------------------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),nunpr
      common/junke/ibl5,ibl52,maxt,maxb,nadd,ncols,nrows,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
c
      dimension q(*),iq(*)
      common/craypk/labs(1360)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
c
c     ibl5 = no of integrals in block of sortfile
c
      ibl5 = nsz340
      iilen = nsz340*mach12/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
c     nat3 = n3n
c     n3n3 = n3n - 3
c     nat33 = nat3 - 3
      nov = nocca*nvirta
c     nov3 = n3n*nov
c
      ncols = nunpr*nov
      nrows = nov
c
c       are going to sort matrix which has nrows rows
c       and ncols cols to get the transpose
c
c       maxt is the number of columns which can be held in core
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = maxq/ncols
      niqq = maxq - nov
      nword = (niqq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
      nbuck = (nrows/maxt) + 1
      nadd = min(maxt,nrows)
      ires = (nrows/maxt) + 1
      maxa = ires*(ibl5+ibl5i) + nov
c
c     nbuck is the number of buckets actually needed
c
      if (nbuck.gt.maxb) then
         write (iwrt,6010) maxq , maxa
         call caserr('insufficient core')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      call rdsrtd(q(i1),iq(ii1),ifortr,n3n,nocca,nvirta,cutoff)
c
c       read through the sort file to give final result
c
      maxqq = nadd*ncols
      call wtsrtd(q(i1),maxqq,ifortw)
      call closbf(0)
      call gmem_free(i1)
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine ddmpti(fcm,nat,ntr)
      implicit REAL  (a-h,o-z)
      dimension fcm(nat*3,nat*3)
c
c     enforces translational invariance on force constant
c     matrix using atom ntr as 'unique' centre
c
c
c     first make sure symmetric (assumes lower triangle
c     calculated correctly)
c
      nat3 = nat*3
      do 30 i = 2 , nat3
         do 20 j = 1 , i - 1
            fcm(j,i) = fcm(i,j)
 20      continue
 30   continue
      if (ntr.le.0) return
c
c     rows ntr*3-2, ntr*3-1 and ntr*3 are missing, so
c     construct them from other rows
c
      ntr1 = ntr*3 - 2
      ntr2 = ntr*3 - 1
      ntr3 = ntr*3
      do 40 j = 1 , nat3
         fcm(ntr1,j) = 0.0d0
         fcm(ntr2,j) = 0.0d0
         fcm(ntr3,j) = 0.0d0
 40   continue
      do 60 i = 1 , nat
         if (i.ne.ntr) then
            do 50 j = 1 , nat3
               fcm(ntr1,j) = fcm(ntr1,j) - fcm(i*3-2,j)
               fcm(ntr2,j) = fcm(ntr2,j) - fcm(i*3-1,j)
               fcm(ntr3,j) = fcm(ntr3,j) - fcm(i*3,j)
 50         continue
         end if
 60   continue
c
c     now repeat for columns ntr1,ntr2,ntr3
c
      do 70 i = 1 , nat3
         fcm(i,ntr1) = 0.0d0
         fcm(i,ntr2) = 0.0d0
         fcm(i,ntr3) = 0.0d0
 70   continue
      do 90 j = 1 , nat
         if (j.ne.ntr) then
            do 80 i = 1 , nat3
               fcm(i,ntr1) = fcm(i,ntr1) - fcm(i,j*3-2)
               fcm(i,ntr2) = fcm(i,ntr2) - fcm(i,j*3-1)
               fcm(i,ntr3) = fcm(i,ntr3) - fcm(i,j*3)
 80         continue
         end if
 90   continue
      return
      end
      subroutine wtsrtd(q,maxq,ifort)
      implicit REAL  (a-h,o-z)
      dimension q(maxq)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk,mkk,g(5119)
      common/sortpk/labs(1)
INCLUDE(common/three)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,ncols,nrows,nbuck
      common/junk/nwbuck(maxbuc)
_IF1(iv)      dimension lin(2)
_IF1(iv)      equivalence (lin(1),g(1))
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call search(1,ifort)
c     an = 1.0d0/dfloat(nops)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     coulumns min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),ncols,ifort)
               j = j + ncols
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nrows) max = nrows
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
_IFN1(iv)            call unpack(g(nsz341),32,labs,ibl5)
_IFN1(civ)            call dsctr(nkk,g,labs,q)
_IF1(c)            call scatter(nkk,q,labs,g)
_IF1(iv)            ij = ibl5+ibl5+1
_IF1(iv)            do 4000 iword=1,nkk
_IF1(iv)            q(lin(ij)) = g(iword)
_IF1(iv) 4000       ij = ij+1
c
            go to 20
         end if
 40   continue
      return
      end
      subroutine mp2dri(iscrf,ifile1,istrma,q,iq,maxq)
      implicit REAL  (a-h,o-z)
      logical skip,skipit
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3)
c
      dimension q(maxq),iq(*)
      common/mpases/mpass,ipass,mpst,mpfi,mpadd,iaoxf1,iaoxf2,moderf
INCLUDE(common/atmblk)
INCLUDE(common/infoa)
INCLUDE(common/prnprn)
c
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/iofile)
c
      cutoff = 10.0d0**(-icut)
      nij = num*(num+1)/2
      nsq = ncoorb*ncoorb
      nov = nocca*nvirta
      n3n = 3*nat
c
      call search(1,iaoxf2)
      if (ipass.eq.1) call search(1,moderf)
c
c      write(iwr,*)'transformation of derivative integrals pass ',ipass
c      write(iwr,*)'perturbations  ',mpst,'  to ',mpfi
c
      call secget(isect(8),8,isec8)
      isec8 = isec8 + mvadd
      call rdedx(q,num*ncoorb,isec8,ifild)
c
      ivec = 1
      iwt = ivec + num*ncoorb
      iws = iwt + nij
      iff = iws + nsq
      itr = iff + nsq
      if = itr + nsq
      it = if + n3n*ncoorb*ncoorb
      itt = it + nocca*nvirta
c     itop = itt + nocca*nvirta
      if (ipass.eq.1) then
         call vclr(q(if),1,n3n*nsq)
         call wrt3(q(if),n3n*nsq,mpblk(2),ifile1)
      end if
c......
c......
c.............start big loop over perturbations
c......
c......
c
c     structure for batches remains intact but now
c     used one perturbation at a time
c
      do 20 ix = mpst , mpfi
c......
         skipit = .false.
         if (skip(ix)) skipit = .true.
         if (skipit.and.odebug(23)) write (iwr,6010) ix
c
c......
         cut1=10.0d0*cutoff
         call mp2tr0(q,q(iwt),q(iws),q(iff),q(itr),nocca,nvirta,ncoorb,
     +  num,nij,iscrf,iaoxf2,skipit,cut1)
c
         if (.not.(skipit)) then
            call mpsrt0(iscrf,iscrf,nij,nov,q,iq,maxq)
c
            call secget(isect(8),8,isec8)
            isec8 = isec8 + mvadd
            call rdedx(q,num*ncoorb,isec8,ifild)
c
            cut2 = cut1
            call rdedx(q(if),n3n*nsq,mpblk(2),ifile1)
            call mp2tr1(q,q(iwt),q(iws),q(iff),q(itr),q(if),q(it),
     +      q(itt),ix,iscrf,moderf,nocca,nvirta,ncoorb,num,nij,n3n,
     +      istrma,cut2)
c
            call wrt3(q(if),n3n*nsq,mpblk(2),ifile1)
         end if
c
 20   continue
c
      return
 6010 format(1x,'perturbation', i4 , ' omitted')
      end
      subroutine rdsrt8(a,ia,ifort,cut)
      implicit REAL  (a-h,o-z)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk1,mkk1,g(5119)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
      common/maxlen/maxq
      dimension ia(*),a(*)
      data lastb/999999/
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c       the first n2 real words of core are used as input
c       workspace.
c
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i) + nov*n3n
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      call vclr(g,1,nsz340+nsz170)
c
      call search(1,ifort)
      iblock = 0
c     ninb = no of elements in bucket(coreload)
      ninb = nadd*nov*n3n
c
      do 50 i = 1 , nov
c
         call rdedz(a(1),nov*n3n,ifort)
c
         do 40 ix = 1 , n3n
            do 30 j = 1 , nov
               iaddr = (j-1)*nov*n3n + (ix-1)*nov + i
c
c     iaddr is address of (i,j) in transposed matrix
c
               if (dabs(a((ix-1)*nov+j)).ge.cut) then
                  ibuck = (iaddr-1)/ninb
                  iaddr = iaddr - ninb*ibuck
                  ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                  nwb = nwbuck(ibuck) + 1
                  a(ibase(ibuck)+nwb) = a((ix-1)*nov+j)
                  ia(ibasen(ibuck)+nwb) = iaddr
                  nwbuck(ibuck) = nwb
                  if (nwb.eq.ibl5) then
c
c     this block full - empty
c
                     call stopbk
                     mkk1 = mark(ibuck)
                     nkk1 = nwb
_IFN1(iv)                     call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)                     call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)                     call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)                     call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
                     call sttout
                     nwbuck(ibuck) = 0
                     mark(ibuck) = iblock
                     iblock = iblock + nsz
                  end if
               end if
 30         continue
 40      continue
 50   continue
c
c      empty anything remaining in buckets
c
      do 60 ibuck = 1 , nbuck
         nwb = nwbuck(ibuck)
         if (nwb.ne.0) then
            call stopbk
            mkk1 = mark(ibuck)
            nkk1 = nwb
_IFN1(iv)            call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)            call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)            call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)            call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
            call sttout
            nwbuck(ibuck) = 0
            mark(ibuck) = iblock
            iblock = iblock + nsz
         end if
 60   continue
c
      call stopbk
      return
      end
      subroutine rdsrt9(a,ia,ifort,nocca,ncoorb,cut)
      implicit REAL  (a-h,o-z)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk1,mkk1,g(5119)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
      common/maxlen/maxq
      dimension ia(*), a(*)
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
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i) + nov*n3n
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
c
      call vclr(g,1,nsz340+nsz170)
c
      call search(1,ifort)
      iblock = 0
c      ninb = no of elements ib bucket(coreload)
      ninb = nadd*nov*n3n
c
      do 70 ia1 = nocca + 1 , ncoorb
         do 60 i1 = 1 , nocca
            call rdedz(a(1),nov*n3n,ifort)
c
            do 50 ib1 = nocca + 1 , ncoorb
               do 40 j1 = 1 , nocca
                  do 30 ix = 1 , n3n
c
                     jjaa = (ia1-nocca-1)*nocca + j1
                     iibb = (ib1-nocca-1)*nocca + i1
                     jjbb = (ib1-nocca-1)*nocca + j1
                     iaddr = (jjaa-1)*nov*n3n + (ix-1)*nov + iibb
                     val = a((ix-1)*nov+jjbb)
c
c
c     iaddr is address of (i,j) in transposed matrix
c
                     if (dabs(val).ge.cut) then
                        ibuck = (iaddr-1)/ninb
                        iaddr = iaddr - ninb*ibuck
                        ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                        nwb = nwbuck(ibuck) + 1
                        a(ibase(ibuck)+nwb) = val
                        ia(ibasen(ibuck)+nwb) = iaddr
                        nwbuck(ibuck) = nwb
                        if (nwb.eq.ibl5) then
c
c     this block full - empty
c
                           call stopbk
                           mkk1 = mark(ibuck)
                           nkk1 = nwb
_IFN1(iv)                           call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)                           call pack(g(nsz341),32,ia(ibasen(ibuck)+1),
_IFN1(iv)     +                               ibl5)
_IF1(iv)                       call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)                       call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
                           call sttout
                           nwbuck(ibuck) = 0
                           mark(ibuck) = iblock
                           iblock = iblock + nsz
                        end if
                     end if
 30               continue
 40            continue
 50         continue
 60      continue
 70   continue
c
c      empty anything remaining in buckets
c
      do 80 ibuck = 1 , nbuck
         nwb = nwbuck(ibuck)
         if (nwb.ne.0) then
            call stopbk
            mkk1 = mark(ibuck)
            nkk1 = nwb
_IFN1(iv)            call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)            call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)            call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)            call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
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
      subroutine rdsrt7(a,ia,ifortr,mpst,mpfi)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      logical skip
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),
     +    nunpr,skip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
c
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk1,mkk1,g(5119)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n3n,nbuck
      common/maxlen/maxq
c     parameter (idblk=3120, iidblk=1950)
c     common/dbuf/val(idblk),vall(iidblk),icnt,mxtr
c     common/dlabs/i1(idblk),j1(idblk),k1(idblk),l1(idblk),
c    +            ipert(idblk)
INCLUDE(common/atmblk)
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/i1(5,3120)
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
INCLUDE(common/mapper)
      dimension ia(*),a(*)
      data lastb/999999/
_IF1()      ind(i,j) = iky(max(i,j)) + min(i,j)
c
c
c       open the sort file
c       each block consists of ibl5 real and ibl5 integer
c       words.
c       ibase gives offset for start of real part of each
c       block, as elements of a real array
c       ibasen gives offset for start of integer part of
c       each block , as elements of integer array.
c
c
      nij2 = nij*nij
      call search(1,ifortr)
      ibl5i = lenint(ibl5)
      do 20 ibuck = 1 , nbuck
         nwbuck(ibuck) = 0
         mark(ibuck) = lastb
         i = (ibuck-1)*(ibl5+ibl5i)
         ibase(ibuck) = i
         ibasen(ibuck) = lenrel(i+ibl5)
 20   continue
      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1
c
      call vclr(g,1,nsz340+nsz170)
c
      iblock = 0
c      ninb = no of elements in bucket(coreload)
      ninb = nadd*nij
c
      mderb=idblk + iidblk + 2
 30   call reads(dbufs,mderb,ifortr)
      if (icnt.eq.0) then
c
c      empty anything remaining in buckets
c
         do 40 ibuck = 1 , nbuck
            nwb = nwbuck(ibuck)
            if (nwb.ne.0) then
               call stopbk
               mkk1 = mark(ibuck)
               nkk1 = nwb
_IFN1(iv)               call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)               call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)               call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)               call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
               call sttout
               nwbuck(ibuck) = 0
               mark(ibuck) = iblock
               iblock = iblock + nsz
            end if
 40      continue
c
c
c     write(iwrt,*)'position, length of ao-file',iposun(ifortr),
c    &   length(ifortr)
c     write(iwrt,*)'position, length of sortfile',iposun(ifiles),
c    &   length(ifiles)
         call stopbk
         return
      else
         call unpack(val(idblk1),lab816,i1,idblk*5)
         do 50 ic = 1 , icnt
            iprt = i1(5,ic)
c
_IF1()c    removed code - relabeling in routine derwrt
_IF1()c    when integrals written out initially
_IF1()c
_IF1()c      if(.not.skip(iprt))goto 999
_IF1()c      if(nops.eq.1)goto 3
_IF1()c      do 10 iop=1,nops
_IF1()c      nupr=iabs(mptr(iop,iprt))
_IF1()c      if(.not.skip(nupr))then
_IF1()c      nnupr=mptr(iop,iprt)
_IF1()c      iprt=nupr
_IF1()c      i11=ibasis(iop,i1(ic))
_IF1()c      j11=ibasis(iop,j1(ic))
_IF1()c      k11=ibasis(iop,k1(ic))
_IF1()c      l11=ibasis(iop,l1(ic))
_IF1()c      i111=iabs(i11)
_IF1()c      j111=iabs(j11)
_IF1()c      k111=iabs(k11)
_IF1()c      l111=iabs(l11)
_IF1()c      isigi=i11/i111
_IF1()c      isigj=j11/j111
_IF1()c      isigk=k11/k111
_IF1()c      isigl=l11/l111
_IF1()c      isigp=nnupr/nupr
_IF1()c      isign=isigi*isigj*isigk*isigl*isigp
_IF1()c      val(ic)=val(ic)*dfloat(isign)
_IF1()c      ij=ind(i111,j111)
_IF1()c      kl=ind(k111,l111)
_IF1()c      goto 999
_IF1()c      end if
_IF1()c10    continue
_IF1()c      call caserr('symmetry error')
_IF1()c999   continue
_IF1()c
c    structure for sorting in batches , mpst to mpfi , remains
c    but now used one at a time
c
            if (iprt.ge.mpst) then
               if (iprt.le.mpfi) then
                  ij = iky(i1(1,ic)) + i1(2,ic)
                  kl = iky(i1(3,ic)) + i1(4,ic)
                  iprt = (iprt-mpst)*nij2
c
c
                  iaddr = iprt + (ij-1)*nij + kl
c
c     iaddr is address of integral in final sequence
c
                  ibuck = (iaddr-1)/ninb
                  iaddr = iaddr - ninb*ibuck
                  ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                  nwb = nwbuck(ibuck) + 1
                  a(ibase(ibuck)+nwb) = val(ic)
                  ia(ibasen(ibuck)+nwb) = iaddr
                  nwbuck(ibuck) = nwb
                  if (nwb.eq.ibl5) then
c
c     this block full - empty
c
                     call stopbk
                     mkk1 = mark(ibuck)
                     nkk1 = nwb
_IFN1(iv)                     call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)                     call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)                     call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)                     call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
                     call sttout
                     nwbuck(ibuck) = 0
                     mark(ibuck) = iblock
                     iblock = iblock + nsz
                  end if
c
c
c     now check if klij integral is different
c
                  if (ij.ne.kl) then
                     iaddr = iprt + (kl-1)*nij + ij
                     ibuck = (iaddr-1)/ninb
                     iaddr = iaddr - ninb*ibuck
                     ibuck = ibuck + 1
c
c     element goes in bucket ibuck with modified address
c
                     nwb = nwbuck(ibuck) + 1
                     a(ibase(ibuck)+nwb) = val(ic)
                     ia(ibasen(ibuck)+nwb) = iaddr
                     nwbuck(ibuck) = nwb
                     if (nwb.eq.ibl5) then
c
c     this block full - empty
c
                     call stopbk
                     mkk1 = mark(ibuck)
                     nkk1 = nwb
_IFN1(iv)                     call dcopy(ibl5,a(ibase(ibuck)+1),1,g,1)
_IFN1(iv)                     call pack(g(nsz341),32,ia(ibasen(ibuck)+1),ibl5)
_IF1(iv)                       call fmove ( a(ibase(ibuck)+1),g,ibl5)
_IF1(iv)                       call fmove ( ia(ibasen(ibuck)+1),g(nsz341),ibl5i)
                     call sttout
                     nwbuck(ibuck) = 0
                     mark(ibuck) = iblock
                     iblock = iblock + nsz
                     end if
c
c
                  end if
               end if
            end if
c
c
 50      continue
         go to 30
      end if
      end
      subroutine mpsrt7(q,iq,num,iwrt)
      implicit REAL  (a-h,o-z)
      common/mpases/mpass,ipass,mpst,mpfi,mpadd,iaoxf1,iaoxf2,moderf
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n3n,nbuck
c     parameter (idblk = 3120)
c     common/dlabs/labs(5*idblk)
      common/dlabs/labs(5,3120)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
INCLUDE(common/prnprn)
      dimension q(*),iq(*)
c
c     ibl5 = no of integrals in block of sortfile
c
      ibl5 = nsz340
      iilen = nsz340*mach12/2
c
      if(o255i) then
       idblk=3120
      else
       idblk=2264
      endif
      call setsto(idblk*5,0,labs)
c
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      n3n = mpfi - mpst + 1
      nij = num*(num+1)/2
c
c      write(iwrt,*)'sorting of derivative integrals pass ',ipass,
c     + ' of ',  mpass
c       are going to sort the derivative integrals from file ifortr
c       so that for m.ge.n all l.ge.s are available on ifortw
c
c       maxt is the (maximum !!) number of triangles
c       which can be held in core
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = (maxq-lenint(nij))/nij
      nword = (maxq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
c
c     nbuck is the number of buckets required,nadd is the no of
c     triangles actually to be held in core
c
      nbuck = (nij*n3n/maxt) + 1
      nadd = min(maxt,nij*n3n)
c
      maxa = nbuck*(ibl5+ibl5i)
      if (ipass.eq.1) then
         write (iwrt,6020)
         if (odebug(39)) write (iwrt,6030) maxq,maxa
      end if
c
      if (nbuck.gt.maxb) then
         ires = (nij*n3n/maxt) + 1
         maxa = ires*(ibl5+ibl5i)
         write (iwrt,6010) maxq , maxa
         call caserr('insufficient core in mpsrt7')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call rdsrt7(q(i1),iq(ii1),iaoxf1,mpst,mpfi)
c
c       read through the sort file to give final result
c
      maxqq = nadd*nij
      i2 = lenint(nij) 
      call wtsrt7(q(i1+i2),maxqq,iaoxf2)
      call gmem_free(i1)
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
 6020 format(/1x,'sorting derivative integrals')
 6030 format(/1x,'memory available ',i7/
     + 1x,'memory required  ',i7/)
      end
      subroutine mpsrt8(q,iq,ifortr,ifortw,nov1,iwrt,cutoff)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),nunpr
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
      dimension q(*),iq(*)
      common/craypk/labs(1360)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
c
c     ibl5 = no of integrals in block of sortfile
c
      ibl5 = nsz340
      iilen = nsz340*mach12/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nov = nov1
      n3n = nunpr
c
c
c       maxt is the number of columns which can be held in core
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = maxq/(nov*n3n)
      niqq = maxq - nov*n3n
      nword = (niqq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
      nbuck = (nov/maxt) + 1
      nadd = min(maxt,nov)
      ires = (nov/maxt) + 1
      maxa = ires*(ibl5+ibl5i) + nov*n3n
c
c     nbuck is the number of buckets actually needed
c
      if (nbuck.gt.maxb) then
         write (iwrt,6010) maxq , maxa
         call caserr('insufficient core')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      call rdsrt8(q(i1),iq(ii1),ifortr,cutoff)
c
c       read through the sort file to give final result
c
      maxqq = nadd*nov*n3n
      call wtsrt8(q(i1),maxqq,ifortw)
      call closbf(0)
      call gmem_free(i1)
c
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine mpsrt9(q,iq,ifortr,ifortw,nov1,nocca,
     * ncoorb,iwrt,cutoff)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dp(27,maxat),nunpr
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
      common/craypk/labs(1360)
INCLUDE(common/blksiz)
INCLUDE(common/atmblk)
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
c
      dimension q(*),iq(*)
c
c     ibl5 = no of integrals in block of sortfile
c
      ibl5 = nsz340
      iilen = nsz340*mach12/2
      call setsto(1360,0,labs)
      ibl52 = iilen
      ibl5i = lenint(ibl5)
      nov = nov1
      n3n = nunpr
c
c
c       maxt is the number of columns which can be held in core
c       which is the number in each bucket
c
      i1  = igmem_alloc_all(maxq)
      ii1 = lenrel(i1-1)+1
      maxt = maxq/(nov*n3n)
      niqq = maxq - nov*n3n
      nword = (niqq/(1+lenrel(1)))*lenrel(1)
c      maxb is the maximum number of blocks of the sortfile
c      which can be held in core
c      which is the maximum number of buckets
c
      maxb = min(maxbuc,nword/ibl5)
      nbuck = (nov/maxt) + 1
      nadd = min(maxt,nov)
      ires = (nov/maxt) + 1
      maxa = ires*(ibl5+ibl5i) + nov*n3n
c
c     nbuck is the number of buckets actually needed
c
      if (nbuck.gt.maxb) then
         write (iwrt,6010) maxq , maxa
         call caserr('insufficient core')
      end if
c
c       read through original file producing sorted file
c
      call vclr(q(i1),1,maxa)
      call setbfa
      call rdsrt9(q(i1),iq(ii1),ifortr,nocca,ncoorb,cutoff)
c
c       read through the sort file to give final result
c
      maxqq = nadd*nov*n3n
      call wtsrt9(q(i1),maxqq,ifortw)
c
      call closbf(0)
      call gmem_free(i1)
      return
 6010 format (//1x,'insufficient core'/1x,'available',i8,'  required',
     +        i8)
      end
      subroutine mp2tr0(vec,wt,ws,ff,tr,nocca,nvirta,ncoorb,num,nij,
     1                      iscrf,idintf,skipit,cut1)
c
      implicit REAL  (a-h,o-z)
      logical skipit
      dimension vec(num*ncoorb),wt(nij),ws(ncoorb*ncoorb),
     1  ff(ncoorb*ncoorb),tr(ncoorb*ncoorb)
c
      call search(1,iscrf)
      iv1 = num*nocca + 1
      do 30 im = 1 , num
         do 20 in = 1 , im
            call rdedz(wt,nij,idintf)
            if (.not.(skipit)) then
               call square(ws,wt,num,num)
               call vclr(ff,1,ncoorb*ncoorb)
               call mxmb(ws,num,1,vec,1,num,ff,nocca,1,num,num,nocca)
               call vclr(tr,1,nocca*nvirta)
               call mxmb(vec(iv1),num,1,ff,nocca,1,tr,nocca,1,nvirta,
     +                   num,nocca)
               do 40 ii=1,nocca*nvirta
                if (dabs(tr(ii)).le.cut1) tr(ii)=0.0d0
 40            continue
               call wtedz(tr,nocca*nvirta,iscrf)
            end if
 20      continue
 30   continue
      return
      end
      subroutine mp2tr1(vec,wt,ws,ff,tr,f,t,tt,ix,iscrf,moderf,nocca,
     1                nvirta,ncoorb,num,nij,n3n,istrma,cut2)
c
      implicit REAL  (a-h,o-z)
      dimension wt(nij),ws(ncoorb*ncoorb),ff(ncoorb*ncoorb),
     1    tr(ncoorb,ncoorb),f(ncoorb,ncoorb,n3n),vec(num*ncoorb)
     2    ,t(nocca*nvirta),tt(nocca*nvirta)
c
c
      nocc1 = nocca + 1
c     nsq = ncoorb*ncoorb
      call search(1,iscrf)
      call search(1,istrma)
c
c
      do 50 ib = nocc1 , ncoorb
         do 40 j = 1 , nocca
            call rdedz(wt,nij,iscrf)
            call square(ws,wt,num,num)
            call vclr(ff,1,ncoorb*ncoorb)
            call mxmb(vec,num,1,ws,1,num,ff,1,ncoorb,ncoorb,num,num)
            call vclr(tr,1,ncoorb*ncoorb)
            call mxmb(ff,1,ncoorb,vec,1,num,tr,1,ncoorb,ncoorb,num,
     +                ncoorb)
c
c
            call rdedz(tt,nocca*nvirta,istrma)
            call mxmb(tr(1,nocc1),1,ncoorb,tt,nocca,1,f(1,1,ix),1,
     +                ncoorb,ncoorb,nvirta,nocca)
            call mxmb(tt,nocca,1,tr,1,ncoorb,f(1,nocc1,ix),ncoorb,1,
     +                nvirta,nocca,ncoorb)
c
            iai = 0
            do 30 ia = nocc1 , ncoorb
               do 20 i = 1 , nocca
                  iai = iai + 1
                  t(iai) = tr(ia,i)
 20            continue
 30         continue
c
            do 60 ii=1,nocca*nvirta
            if (dabs(t(ii)).le.cut2) t(ii)=0.0d0
 60         continue
c
            call wtedz(t,nocca*nvirta,moderf)
 40      continue
 50   continue
c
      return
      end
      subroutine wtsrt7(q,maxq,ifortw)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      dimension q(maxq)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk,mkk,g(5119)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nij,n3n,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
_IF1(iv)      dimension lin(2)
_IF1(iv)      equivalence (lin(1),g(1))
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call search(1,ifortw)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 50 i = 1 , nbuck
         call vclr(q,1,maxq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     triangles min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),nij,ifortw)
               j = j + nij
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nij*n3n) max = nij*n3n
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
_IFN1(iv)            call unpack(g(nsz341),32,labs,ibl5)
_IFN1(iv)            do 40 iword = 1 , nkk
_IFN1(iv)               q(labs(iword)) = q(labs(iword)) + g(iword)
_IFN1(iv) 40         continue
_IF1(iv)            ij=ibl5+ibl5+1
_IF1(iv)            do 40 iword=1,nkk
_IF1(iv)             q(lin(ij)) = q(lin(ij)) + g(iword)
_IF1(iv) 40         ij = ij+1
            go to 20
         end if
 50   continue
      return
      end
      subroutine wtsrt8(q,maxq,ifort)
      implicit REAL  (a-h,o-z)
      dimension q(maxq)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk,mkk,g(5119)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
_IF1(iv)      dimension lin(2)
_IF1(iv)      equivalence (lin(1),g(1))
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call search(1,ifort)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     coulumns min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),nov*n3n,ifort)
               j = j + nov*n3n
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nov) max = nov
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
_IFN1(iv)            call unpack(g(nsz341),32,labs,ibl5)
_IFN1(civ)            call dsctr(nkk,g,labs,q)
_IF1(c)            call scatter(nkk,q,labs,g)
_IF1(iv)            ij = ibl5+ibl5+1
_IF1(iv)            do 4000 iword=1,nkk
_IF1(iv)            q(lin(ij)) = g(iword)
_IF1(iv) 4000       ij = ij+1
c
            go to 20
         end if
 40   continue
      return
      end
      subroutine wtsrt9(q,maxq,ifort)
      implicit REAL  (a-h,o-z)
      dimension q(maxq)
INCLUDE(common/stak)
INCLUDE(common/blksiz)
      common/bufb/nkk,mkk,g(5119)
      common/sortpk/labs(1)
      common/junke/ibl5,ibl52,maxt,maxb,nadd,nov,n3n,nbuck
INCLUDE(common/three)
      common/junk/nwbuck(maxbuc)
_IF1(iv)      dimension lin(2)
_IF1(iv)      equivalence (lin(1),g(1))
      data lastb/999999/
c
c    read thru the sort file to get core load of elements then
c    write them out on sequential file
c
      call search(1,ifort)
c
      min = 1
      max = nadd
c
c     loop over buckets
c
      do 40 i = 1 , nbuck
         call vclr(q,1,maxq)
         mkk = mark(i)
 20      if (mkk.eq.lastb) then
c
c     coulumns min thru max are in core - clear them out
c
            j = 1
            do 30 n = min , max
               call wtedz(q(j),nov*n3n,ifort)
               j = j + nov*n3n
 30         continue
            min = min + nadd
            max = max + nadd
            if (max.gt.nov) max = nov
         else
c
c     loop over the sortfile blocks comprising this bucket
c
            iblock = mkk
            call rdbak(iblock)
            call stopbk
_IFN1(iv)            call unpack(g(nsz341),32,labs,ibl5)
_IFN1(civ)            call dsctr(nkk,g,labs,q)
_IF1(c)            call scatter(nkk,q,labs,g)
_IF1(iv)            ij = ibl5+ibl5+1
_IF1(iv)            do 4000 iword=1,nkk
_IF1(iv)            q(lin(ij)) = g(iword)
_IF1(iv) 4000       ij = ij+1
c
            go to 20
         end if
 40   continue
      return
      end
      subroutine ver_secmp2(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/secmp2.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
