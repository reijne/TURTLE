c     deck=drv2e
c ******************************************************
c ******************************************************
c             =   drv2e   =
c ******************************************************
c ******************************************************
      subroutine dabab(ii,jj,kk,ll,q4,zscftp,da,db,abdens)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/incrd)
INCLUDE(common/nshel)
      common/restrl/ociopt,ocifor,omp2
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/ghfblk)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
      dimension da(*),db(*),abdens(*)
      data zuhf,zgrhf /'uhf','grhf'/
      data pt5,four /0.5d0,4.0d0/
c
      ouhf = zscftp.eq.zuhf
      ogrhf = zscftp.eq.zgrhf
_IF(ccpdft)
      ohf_coul = (.not.CD_active()).or.CD_HF_coulomb_deriv()
      ohf_exch = (.not.CD_active()).or.CD_HF_exchange()
      hf_wght  = CD_HF_exchange_weight()
_ENDIF
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
c
      if (.not.ogrhf) then
c
         do 60 i = mini , maxi
            nj = ni
            do 50 j = minj , maxj
               nk = nj
               do 40 k = mink , maxk
                  nl = nk
                  do 30 l = minl , maxl
                     nn = nl
                     i1 = loci + i
                     i2 = locj + j
                     i3 = lock + k
                     i4 = locl + l
                     if (i1.lt.i2) then
                        n = i1
                        i1 = i2
                        i2 = n
                     end if
                     if (i3.lt.i4) then
                        n = i3
                        i3 = i4
                        i4 = n
                     end if
                     if (i1.lt.i3) then
                     else if (i1.eq.i3) then
                        if (i2.ge.i4) go to 20
                     else
                        go to 20
                     end if
                     n = i1
                     i1 = i3
                     i3 = n
                     n = i2
                     i2 = i4
                     i4 = n
 20                  mij = iky(i1) + i2
                     mik = iky(i1) + i3
                     mil = iky(i1) + i4
                     mkl = iky(i3) + i4
                     if (i2.lt.i3) then
                        mjk = iky(i3) + i2
                        if (i2.lt.i4) then
                           mjl = iky(i4) + i2
                        else
                           mjl = iky(i2) + i4
                        end if
                     else
                        mjk = iky(i2) + i3
                        mjl = iky(i2) + i4
                     end if
_IF(ccpdft)
                     dfac=0.0d0
                     if(ohf_coul)then
                        dfac = da(mij)*da(mkl)*four 
                     endif
                     if(ohf_exch)then
                        dfac = dfac - hf_wght*(da(mik)*da(mjl) 
     &                                       + da(mil)*da(mjk))
                        if(ouhf)dfac = dfac - hf_wght*(db(mik)*db(mjl)
     &                                               + db(mil)*db(mjk))
                     endif
_ELSE
                     dfac = da(mij)*da(mkl)*four - da(mik)*da(mjl)
     +                      - da(mil)*da(mjk)
                     if (ouhf) dfac = dfac - db(mik)*db(mjl) - db(mil)
     +                                *db(mjk)
_ENDIF
                     if (i1.eq.i2) dfac = dfac*pt5
                     if (i3.eq.i4) dfac = dfac*pt5
                     dfac = dfac*q4
                     if (omp2 .or. mp3) then
                        abdens(nn) = abdens(nn) + dfac
                     else
                        abdens(nn) = dfac
                     end if
                     nl = nl + inc2
 30               continue
                  nk = nk + inc3
 40            continue
               nj = nj + inc4
 50         continue
            ni = ni + inc5
 60      continue
c
      else
c
c     general case
c
         do 120 i = mini , maxi
            nj = ni
            i1 = loci + i
            ii1 = iky(i1)
            do 110 j = minj , maxj
               nk = nj
               j1 = locj + j
               jj1 = iky(j1)
               mij = ii1 + j1
               if (j1.gt.i1) mij = jj1 + i1
               oijeq = i1.eq.j1
               do 100 k = mink , maxk
                  nl = nk
                  k1 = lock + k
                  kk1 = iky(k1)
                  mik = ii1 + k1
                  if (k1.gt.i1) mik = kk1 + i1
                  mjk = jj1 + k1
                  if (k1.gt.j1) mjk = kk1 + j1
                  do 90 l = minl , maxl
                     nn = nl
                     l1 = locl + l
                     ll1 = iky(l1)
                     mkl = kk1 + l1
                     if (l1.gt.k1) mkl = ll1 + k1
                     mjl = jj1 + l1
                     if (l1.gt.j1) mjl = ll1 + j1
                     mil = ii1 + l1
                     if (l1.gt.i1) mil = ll1 + i1
                     okleq = k1.eq.l1
                     dfac = 0.0d0
                     ioff = 0
                     do 80 is = 1 , njk
                        isi = (is-1)*11
                        joff = 0
                        do 70 js = 1 , njk
                           dfac = dfac + 4.0d0*erga(isi+js)
     +                            *(da(ioff+mij)*da(joff+mkl)
     +                            +da(ioff+mkl)*da(joff+mij))
     +                            + 2.0d0*ergb(isi+js)
     +                            *(da(ioff+mik)*da(joff+mjl)
     +                            +da(ioff+mjl)*da(joff+mik)
     +                            +da(ioff+mil)*da(joff+mjk)
     +                            +da(ioff+mjk)*da(joff+mil))
                           joff = joff + nx
 70                     continue
                        ioff = ioff + nx
 80                  continue
                     if (oijeq) dfac = dfac*pt5
                     if (okleq) dfac = dfac*pt5
                     abdens(nn) = dfac*q4
                     nl = nl + inc2
 90               continue
                  nk = nk + inc3
 100           continue
               nj = nj + inc4
 110        continue
            ni = ni + inc5
 120     continue
c
      end if
c
      return
      end
      subroutine dabg(ii,jj,kk,ll,l1,norb,q4,da,v,nconf,
     + onocor,onopen,abdens)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/incrd)
INCLUDE(common/nshel)
INCLUDE(common/scfwfn)
      dimension abdens(*),da(*),v(l1,*),nconf(*)
      data dzero,pt5,eight/0.0d0,0.5d0,8.0d0/
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
      do 100 i = mini , maxi
         nj = ni
         do 90 j = minj , maxj
            nk = nj
            do 80 k = mink , maxk
               nl = nk
               do 70 l = minl , maxl
                  nn = nl
                  i1 = loci + i
                  i2 = locj + j
                  i3 = lock + k
                  i4 = locl + l
                  if (i1.lt.i2) then
                     n = i1
                     i1 = i2
                     i2 = n
                  end if
                  if (i3.lt.i4) then
                     n = i3
                     i3 = i4
                     i4 = n
                  end if
                  if (i1.lt.i3) then
                  else if (i1.eq.i3) then
                     if (i2.ge.i4) go to 20
                  else
                     go to 20
                  end if
                  n = i1
                  i1 = i3
                  i3 = n
                  n = i2
                  i2 = i4
                  i4 = n
 20               mij = iky(i1) + i2
                  mik = iky(i1) + i3
                  mil = iky(i1) + i4
                  mkl = iky(i3) + i4
                  if (i2.lt.i3) then
                     mjk = iky(i3) + i2
                     if (i2.lt.i4) then
                        mjl = iky(i4) + i2
                     else
                        mjl = iky(i2) + i4
                     end if
                  else
                     mjk = iky(i2) + i3
                     mjl = iky(i2) + i4
                  end if
                  dfac = dzero
                  if (.not.(onocor)) then
                     dfac = dfac + alpha(1)*da(mij)*da(mkl)
     +                      + pt5*beta(1)
     +                      *(da(mik)*da(mjl)+da(mil)*da(mjk))
                     if (onopen) go to 60
                     nco1 = nco + 1
                     do 30 io = nco1 , norb
                        iojo = iky(nconf(io)) + 1
                        dfac = dfac + alpha(iojo)
     +                         *(da(mij)*v(i3,io)*v(i4,io)+da(mkl)
     +                         *v(i1,io)*v(i2,io)) + pt5*beta(iojo)
     +                         *(v(i2,io)
     +                         *(da(mik)*v(i4,io)+da(mil)*v(i3,io))
     +                         +v(i1,io)
     +                         *(da(mjl)*v(i3,io)+da(mjk)*v(i4,io)))
 30                  continue
                  end if
                  if (.not.(onopen)) then
                     nco1 = nco + 1
                     do 50 io = nco1 , norb
                        do 40 jo = nco1 , norb
                           iof = nconf(io)
                           jof = nconf(jo)
                           iojo = iky(iof) + jof
                           if (jof.gt.iof) iojo = iky(jof) + iof
                           dfac = dfac + alpha(iojo)*v(i1,io)*v(i2,io)
     +                            *v(i3,jo)*v(i4,jo) + pt5*beta(iojo)
     +                            *v(i1,io)*v(i2,jo)
     +                            *(v(i3,io)*v(i4,jo)+v(i4,io)*v(i3,jo))
 40                     continue
 50                  continue
                  end if
 60               if (i1.eq.i2) dfac = dfac*pt5
                  if (i3.eq.i4) dfac = dfac*pt5
                  dfac = dfac*eight*q4
                  abdens(nn) = dfac
                  nl = nl + inc2
 70            continue
               nk = nk + inc3
 80         continue
            nj = nj + inc4
 90      continue
         ni = ni + inc5
 100  continue
      return
      end
      subroutine dabmc(ii,jj,kk,ll,q4,db,dc,abdens)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/incrd)
INCLUDE(common/nshel)
INCLUDE(common/cigrad)
c
      dimension dc(*),db(*),abdens(*)
      data half,four /0.5d0,4.0d0/
c
      if (ncore.eq.0) return
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
c
c    additions for ci core
c
      ni = 1
      do 50 i = mini , maxi
         nj = ni
         i1 = loci + i
         ii1 = iky(i1)
         do 40 j = minj , maxj
            nk = nj
            j1 = locj + j
            jj1 = iky(j1)
            mij = ii1 + j1
            if (j1.gt.i1) mij = jj1 + i1
            oijeq = i1.eq.j1
            do 30 k = mink , maxk
               nl = nk
               k1 = lock + k
               kk1 = iky(k1)
               mik = ii1 + k1
               if (k1.gt.i1) mik = kk1 + i1
               mjk = jj1 + k1
               if (k1.gt.j1) mjk = kk1 + j1
               do 20 l = minl , maxl
                  nn = nl
                  l1 = locl + l
                  ll1 = iky(l1)
                  mkl = kk1 + l1
                  if (l1.gt.k1) mkl = ll1 + k1
                  mjl = jj1 + l1
                  if (l1.gt.j1) mjl = ll1 + j1
                  mil = ii1 + l1
                  if (l1.gt.i1) mil = ll1 + i1
                  okleq = k1.eq.l1
c
c     db ----- core density
c     dc ----- ci density
c
                  dfaca = db(mij)*db(mkl)*four - db(mik)*db(mjl)
     +                    - db(mjk)*db(mil)
                  dfacb = (db(mij)*dc(mkl)+db(mkl)*dc(mij))
     +                    *four - db(mil)*dc(mjk) - db(mjk)*dc(mil)
     +                    - db(mik)*dc(mjl) - db(mjl)*dc(mik)
                  dfac = dfaca*four + dfacb + dfacb
                  if (oijeq) dfac = dfac*half
                  if (okleq) dfac = dfac*half
                  abdens(nn) = abdens(nn) + dfac*q4
                  nl = nl + inc2
 20            continue
               nk = nk + inc3
 30         continue
            nj = nj + inc4
 40      continue
         ni = ni + inc5
 50   continue
      return
      end
      subroutine dabout(q,odebug,iw)
c -----------------------------------
c    write out results from fokabd
c ------------------------------------
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/specal)
INCLUDE(common/incrd)
INCLUDE(common/infoa)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      dimension q(*)
      data half/0.5d0/
c
c     writes out accumulated derivative fock
c     onto fockfile.
c
      if (odebug) write (iw,6010)
      ioff = ifok
      ifocbl = ibdout
      do 50 n = 1 , nfok
         call rdedx(q(ida+1),nx,ifocbl,iflout)
         ij = 0
         do 30 k = 1 , num
            do 20 l = 1 , k
               ij = ij + 1
               if (k.ne.l) then
                  q(ioff+ij) = half*q(ioff+ij)
               end if
 20         continue
 30      continue
         do 40 k = 1 , nx
            q(ida+k) = q(ida+k) + q(ioff+k)
 40      continue
         call wrt3(q(ida+1),nx,ifocbl,iflout)
         if (odebug) call prtris(q(ida+1),num,iw)
         ioff = ioff + nx
         ifocbl = ifocbl + lenb
 50   continue
c
      return
 6010 format (///1x,'output from dabout'//)
      end
      subroutine ddebut(zscftp,da,db,nconf,dd,fock,
     + l1,l2,l3,ista,jsta,ksta,lsta,q)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/prnprn)
      common/restrl/ociopt(2),omp2
INCLUDE(common/specal)
INCLUDE(common/restar)
INCLUDE(common/cndx41)
INCLUDE(common/dmisc)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/scfwfn)
INCLUDE(common/mapper)
      dimension da(*),db(*),nconf(*),fock(*),dd(*), q(*)
INCLUDE(common/grad2)
INCLUDE(common/ijlab)
INCLUDE(common/ghfblk)
INCLUDE(common/cigrad)
INCLUDE(common/atmblk)
_IF(parallel)
INCLUDE(common/infoa)
_ENDIF
c
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
      data zvb,zmcscf/'vb','mcscf'/
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
c
      if (onocnt) write (iwr,6010)
c
      if (zscftp.eq.zgvb) then
c
c     ----- set up core density matrix, read in  eigenvectors for
c            gvb -----
c
         call rdedx(db,l3,ibl3qa,idaf)
         call tdown(db,ilifq,db,ilifq,l1)
         call dencor(da,db,l1)
c
c     ----- set up mo to fock operator pointers in nconf -----
c
         if (nco.ne.0) call setsto(nco,1,nconf)
         ic = 0
         if (nseto.gt.0) then
            nbase = ncores
            do 30 i = 1 , nseto
               nop = no(i)
               do 20 j = 1 , nop
                  nconf(ic+nco+j) = nbase + 1
 20            continue
               ic = ic + nop
               nbase = nbase + 1
 30         continue
         end if
         if (npair.gt.0) then
            np2 = npair + npair
            do 40 i = 1 , np2
               nconf(i+nco+ic) = ncores + nseto + i
 40         continue
         end if
         norb = nco + npair + npair + ic
         onocor = .false.
         onopen = .false.
         onocor = nco.eq.0
         onopen = nseto.eq.0 .and. npair.eq.0
c
      else if (zscftp.eq.zuhf) then
c
c     ----- read in density matrices (alpha+beta) in uhf -----
c
         call rdedx(da,l2,ibl3pa,idaf)
         call rdedx(db,l2,ibl3pb,idaf)
         do 50 i = 1 , l2
            duma = da(i)
            dumb = db(i)
            da(i) = duma + dumb
            db(i) = duma - dumb
 50      continue
c
      else if (zscftp.eq.zrhf) then
c
c     ----- read in density matrix  in rhf  -----
c
_IF(mp2_parallel)
c gdf:  skip I/O for parallel mp2 gradient
        if (.not.omp2) call rdedx(da,l2,ibl3pa,idaf)
_ELSE
         call rdedx(da,l2,ibl3pa,idaf)
_ENDIF
c
      else if (zscftp.eq.zgrhf) then
c
c     general scf
c
         call rdedx(fock,l3,ibl3qa,idaf)
         call tdown(fock,ilifq,fock,ilifq,l1)
         m = 0
         do 90 is = 1 , njk
            call vclr(da(m+1),1,l2)
            nsi = nbshel(is)
            ils = ilfshl(is)
            do 80 ni = 1 , nsi
               nc = iactiv(ils+ni)
               ncol = (nc-1)*l1
               ij = 0
               do 70 i = 1 , l1
                  dum = fock(i+ncol)
                  do 60 j = 1 , i
                     ij = ij + 1
                     da(ij+m) = da(ij+m) + dum*fock(j+ncol)
 60               continue
 70            continue
 80         continue
            m = m + l2
 90      continue
c
      else if ((zscftp.eq.zmcscf).or.(zscftp.eq.zvb)) then
c
c     mcscf/multi
c
         m = 0
         call secget(isecmo,m,iblok)
         call rdedx(da,l1*ncoorb,iblok+mvadd,idaf)
         m = 0
         call secget(isecdd,m,iblok)
         if(odebug(30)) write (iwr,6040) iblok
         call rdedx(dd,l2,iblok,idaf)
         if(odebug(30)) write (iwr,6050) ifil2d,iblk2d
         call rdedx(fock,nd2mo,iblk2d,ifil2d)
         if (ncore.gt.0) then
            call vclr(db,1,l2)
            ij = 0
            do 120 i = 1 , l1
               do 110 j = 1 , i
                  ij = ij + 1
                  do 100 k = 1 , ncore
                     kk = (k-1)*l1
                     db(ij) = db(ij) + da(i+kk)*da(j+kk)
 100              continue
 110           continue
 120        continue
            write (iwr,6020) ncore
         end if
      else
         call caserr('invalid scftype detected in gradient code')
      end if
c
c     ----- set starting parameters -----
c
      outd= nprint.eq. - 4
      if (ofokab .or. ompir) then
c
c      read in ndenin density matrices from some external file
c      matrices are in mo basis
c
         ioff = 1
         call search(ibden,iflden)
         do 130 loop = 1 , ndenin
            call reads(db(ioff),l2,iflden)
            if (odebug(21)) call prtris(db(ioff),l1,iwr)
            ioff = ioff + l2
 130     continue

         ioffd = igmem_alloc(l2)
         iofft = igmem_alloc(l1)
         ioffv = igmem_alloc(l3)
         m = 0
         call rdedx(q(ioffv),l3,ibl3qa,idaf)
         call tdown(q(ioffv),ilifq,q(ioffv),ilifq,l1)
         ioff = 1
         do 140 loop = 1 , ndenin
            call dcopy(l2,db(ioff),1,q(ioffd),1)
            call demoao(q(ioffd),db(ioff),q(ioffv),q(iofft),l1,
     +  ncoorb,l1)
            ioff = ioff + l2
 140     continue

         call gmem_free(ioffv)
         call gmem_free(iofft)
         call gmem_free(ioffd)

      end if
      if (omp2 .or. ofokab .or. ofock) then
         icutd = icut + 1
      else
         icutd = icut
      end if
      icutd = max(icutd,8)
      if (omp2 .or. ofokab .or. ofock) then
         itold = itol + 2
      else
         itold = itol + 1
      end if
      itold = max(itold,16)
      dcut = done/(ten**icutd)
      tol1 = e*(itold+2)
      tol2 = e*itold
      tol3 = done/ten**(itold+2)
      tol4 = done/ten**itold
      onorm = normf.ne.1 .or. normp.ne.1
c
c     ----- read in 1e-gradient -----
c
_IF(parallel)
c***   ***node-MPP***
      call rdgrd(de,irest,ista,jsta,ksta,lsta)
c***    now divide by number of nodes to get everything ok
      factor = 1.0d0/(ipg_nnodes()*1.0d0)
      call dscal(nat*3,factor,de,1)
c***   ***node-MPP***
_ELSE
      call rdgrd(de,irest,ista,jsta,ksta,lsta)
_ENDIF
c
      if (ofock .or. ofokab) then
         if(odebug(30)) write (iwr,6030)
         call vclr(fock,1,nfok*l2)
      end if
      return
 6010 format (/1x,22('=')/1x,'gradient of the energy'/1x,22('='))
 6020 format (/1x,'number of frozen and core orbitals',i5/)
 6030 format (/1x,'zeroing storage for fock matrices')
 6040 format (/1x,'reading density matrix from block',i5)
 6050 format (/1x,'tpdm from file, block', 2i5)
      end
      subroutine delim(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd,
     * abmax)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/dmisc)
      common/restrl/ociopt(2),omptwo,ohf(7),omp2w
INCLUDE(common/specal)
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      do 20 i = 1 , lendd
         oform(i) = .false.
 20   continue
      abmax = 1.0d0
      if (.not.(ofock .or. ofokab .or. ompir .or. omp2w)) then
         abmax = 0.0d0
         nn = 0
         do 40 i = 1 , ij
            do 30 k = 1 , kl
               nn = nn + 1
               n = ijgt(i) + klgt(k)
               ab = dabs(qa(n))
               if (ab.lt.dcut) oform(nn) = .true.
               if (ab.gt.abmax) abmax = ab
 30         continue
 40      continue
      end if
      call dcopy(lendd,qa,1,qb,1)
      nn = 0
      ijkl = 0
      do 60 i = 1 , ij
         do 50 k = 1 , kl
            nn = nn + 1
            if (.not.(oform(nn))) then
               n = ijgt(i) + klgt(k)
               ijkl = ijkl + 1
               qa(ijkl) = qb(n)
            end if
 50      continue
 60   continue
      return
      end
      subroutine delim2(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      call dcopy(lendd,qa,1,qb,1)
      nn = 0
      ijkl = 0
      do 30 i = 1 , ij
         do 20 k = 1 , kl
            nn = nn + 1
            if (.not.(oform(nn))) then
               n = ijgt(i) + klgt(k)
               ijkl = ijkl + 1
               qa(ijkl) = qb(n)
            end if
 20      continue
 30   continue
      return
      end
      subroutine dencor(da,v,l1)
c     ----- calculate the density matrix for the core orbitals -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/scfwfn)
      dimension da(*),v(l1,*)
c
      if (nco.ne.0) then
         l2 = l1*(l1+1)/2
         call vclr(da,1,l2)
         j = 1
         do 30 i = 1 , l1
            do 20 k = 1 , nco
               fact = f(1)*v(i,k)
               if (fact.ne.0.0d0) then
                  call daxpy(i,fact,v(1,k),1,da(j),1)
               end if
 20         continue
            j = j + i
 30      continue
      end if
      return
      end
      subroutine denfac(dkl,csk,cpk,cdk,cfk,cgk,
     +                      csl,cpl,cdl,cfl,cgl,
     +                     mink,maxk,minl,maxl,okandl,double)
      implicit REAL  (a-h,o-z)
       logical okandl,double
       dimension dkl(*)
      if (.not.double) then
         n = 0
         mmax = maxl
         do 130 k = mink , maxk
            if (okandl) mmax = k
            go to (20,30,60,60,
     +             40,60,60,60,60,60,
     +             50,60,60,60,60,60,60,60,60,60,
     +             55,60,60,60,60,60,60,60,60,60,
     +             60,60,60,60,60) , k
 20         dum1 = csk
            go to 60
 30         dum1 = cpk
            go to 60
 40         dum1 = cdk
            go to 60
 50         dum1 = cfk
            go to 60
 55         dum1 = cgk
 60         do 120 l = minl , mmax
               go to (70,80,110,110,
     +                90,110,110,110,110,110,
     +               100,110,110,110,110,110,110,110,110,110,
     +               105,110,110,110,110,110,110,110,110,110,
     +               110,110,110,110,110) , l
 70            dum2 = dum1*csl
               go to 110
 80            dum2 = dum1*cpl
               go to 110
 90            dum2 = dum1*cdl
               go to 110
 100           dum2 = dum1*cfl
               go to 110
 105           dum2 = dum1*cgl
 110           n = n + 1
               dkl(n) = dum2
 120        continue
 130     continue
      else
         n = 0
         mmax = maxl
         do 250 k = mink , maxk
            if (okandl) mmax = k
            go to (140,150,180,180,
     +             160,180,180,180,180,180,
     +             170,180,180,180,180,180,180,180,180,180,
     +             175,180,180,180,180,180,180,180,180,180,
     +             180,180,180,180,180) , k
 140        dum1 = csk
            go to 180
 150        dum1 = cpk
            go to 180
 160        dum1 = cdk
            go to 180
 170        dum1 = cfk
            go to 180
 175        dum1 = cgk
 180        do 240 l = minl , mmax
               go to (190,200,230,230,
     +                210,230,230,230,230,230,
     +                220,230,230,230,230,230,230,230,230,230,
     +                225,230,230,230,230,230,230,230,230,230,
     +                230,230,230,230,230) , l
 190           dum2 = dum1*csl
               if (k.gt.1) then
                  dum2 = dum2 + csk*cpl
               else
                  dum2 = dum2 + dum2
               end if
               go to 230
 200           dum2 = dum1*(cpl+cpl)
               go to 230
 210           dum2 = dum1*(cdl+cdl)
               go to 230
 220           dum2 = dum1*(cfl+cfl)
               go to 230
 225           dum2 = dum1*(cgl+cgl)
 230           n = n + 1
               dkl(n) = dum2
 240        continue
 250     continue
      end if
      return
      end
_IFN(qmmm)
      subroutine derwrt(qq)
      implicit REAL  (a-h,o-z)
c
INCLUDE(common/sizes)
c
      dimension qq(*)
c
INCLUDE(common/cndx41)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt3,nt4,
     & itable(8,8),mult(8,8),nirr(maxorb)
      logical xskip
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +    dpqd(27,maxat),
     +    nunpr,xskip(maxat*3),mapnu(maxat*3),mapun(maxat*3),
     +    iperm(maxat*3),
     +    ict(maxat,8),mptr(8,maxat*3)
INCLUDE(common/dshlno)
INCLUDE(common/dmisc)
c
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/restar)
INCLUDE(common/specal)
c
INCLUDE(common/atmblk)
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/ipp(5,3120)
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
INCLUDE(common/incrd)
      dimension mpert(12)
c
      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1
c
      cut = 1.0d0 / (10.0d0**(icut-1))
      npp = 0
      do 30 i = 1 , npass + 1
         do 20 j = 1 , 3
            npp = npp + 1
            mpert(npp) = (natomd(i)-1)*3 + j
 20      continue
 30   continue
c
      np3 = (npass+1)*3
      nn = 0
      jmax = maxj
      do 80 i = mini , maxi
         i1 = loci + i
         if (oiandj) jmax = i
         do 70 j = minj , jmax
            i2 = locj + j
            lmax = maxl
            if (osame) maxk = i
            do 60 k = mink , maxk
               i3 = lock + k
               if (okandl) lmax = k
               if (osame .and. i.eq.k) lmax = j
               do 50 l = minl , lmax
                  i4 = locl + l
                  nn = nn + 1
cc
                  if (osame) then
                     ijkl1 = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                       *inc3 + l - minl + 1
                     ijkl2 = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                       *inc3 + j - minj + 1
                  end if
                  fac = 1.0d0
                  if (i1.eq.i2) fac = fac + fac
                  if (i3.eq.i4) fac = fac + fac
cc
                  ioa = 0
                  do 40 npp = 1 , np3
                     derval = qq(nn+ioa)*fac
                     if (osame) derval = qq(ijkl1+ioa) + qq(ijkl2+ioa)
                     if (dabs(derval).ge.cut) then
                        if (iperm(mpert(npp)).ne.0) then
                           icnt = icnt + 1
                           nap = mpert(npp)
                           j1 = i1
                           j2 = i2
                           j3 = i3
                           j4 = i4
                           if (xskip(nap)) then
                              jop = iperm(nap)
                              j1 = ibasis(jop,i1)
                              j2 = ibasis(jop,i2)
                              j3 = ibasis(jop,i3)
                              j4 = ibasis(jop,i4)
                              nupert = mptr(jop,nap)
                              sigh = 
     +                           dsign(1.0d0,dfloat(j1*j2)*
     +                                       dfloat(j3*j4*nupert))
                              derval = derval*sigh
                              jj1 = iabs(j1)
                              jj2 = iabs(j2)
                              jj3 = iabs(j3)
                              jj4 = iabs(j4)
                              j1 = max(jj1,jj2)
                              j2 = min(jj1,jj2)
                              j3 = max(jj3,jj4)
                              j4 = min(jj3,jj4)
                              nap = iabs(nupert)
                           end if
                           val(icnt) = derval
                           if (j1.ge.j3) then
                              ipp(1,icnt) = j1
                              ipp(2,icnt) = j2
                              ipp(3,icnt) = j3
                              ipp(4,icnt) = j4
                           else
                              ipp(1,icnt) = j3
                              ipp(2,icnt) = j4
                              ipp(3,icnt) = j1
                              ipp(4,icnt) = j2
                           end if
                           if (j1.eq.j3 .and. j4.gt.j2) then
                              ipp(4,icnt) = j2
                              ipp(2,icnt) = j4
                           end if
                           ipp(5,icnt) = nap
                           if (icnt.eq.idblk) then
                              call pack(val(idblk1),lab816,ipp,idblk*5)
                              mderb = idblk + iidblk + 2
                              call wrt3s(dbufs,mderb,mpstrm(1))
                              icnt = 0
                           end if
                        end if
                     end if
                     ioa = ioa + lendd
 40               continue
 50            continue
 60         continue
 70      continue
 80   continue
      return
      end
_ELSE
      subroutine derwrt(dum)
         call caserr('not available in qmmm version')
      end
_ENDIF

      subroutine dfinal(q,index,ii)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
_IF(parallel)
c***   ***node-MPP***
c***    writing (and gopping) for index ne 0 is disabled
c***    (or if other integrals still follow)
c***    so no restarts (nb. and no d-gradients!!)
c***    this saves quite a bit of bother
c***    restarts may be reeabled by dividing the partial grads by
c***    the number of nodes after gopping before writing
c***    (the restarting job must have same number then (tricky)
c***   ***node-MPP***
_ENDIF
INCLUDE(common/prnprn)
INCLUDE(common/specal)
INCLUDE(common/cndx41)
INCLUDE(common/restar)
INCLUDE(common/timez)
INCLUDE(common/iofile)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/grad2)
INCLUDE(common/ijlab)
INCLUDE(common/funct)
_IF(charmm)
INCLUDE(common/chmgms)
_ENDIF
_IFN(qmmm)
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
_ENDIF
INCLUDE(common/restri)
INCLUDE(common/prints)
INCLUDE(common/runlab)
INCLUDE(common/runopt)
      dimension f(3,maxat),ydnam(3),q(*)
      dimension ibuffi(8)
      character*10 charwall
      equivalence (egrad(1),f(1,1))
      data ydnam /'e/x','e/y','e/z'/
c
      if (index.eq.1) then
c
c...   this might be elimitated if the gradients are properly bypassed
c
         do n=1,nat
            if (zopti(n).eq.'no') then
                de(1,n) = 0.0d0
                de(2,n) = 0.0d0
                de(3,n) = 0.0d0
            end if
         end do
      end if
c
      if (index.ne.1) then
c
c     ----- get restart data -----
c
         irest = 6
         istnu = 1 + ii
         jstnu = 1
         kstnu = 1
         lstnu = 1
c
c     ----- save gradient + restart data -----
c
_IF(parallel)
c***   ***node-MPP***
c     call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
c***   ***node-MPP***   (clear de to receive further contribs)
_ELSE
      call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
_ENDIF
c
c     ----- check cpu time -----
c
         if (istnu.gt.nshell) return
         call texit(0,irest)
         if (tim.lt.timlim) return
         write (iwr,6010) tim , charwall(), istnu , jstnu , kstnu , lstnu
c
         if (ofokab) call dabout(q,odebug(21),iwr)
         if (ofock) call hsymd(q,iwr)
         write (iwr,6030)
         return
      else if (onocnt) then
c
c     ----- transfer gradient into -egrad- -----
c
         ncoord = 3*nat
_IF(parallel)
c***   ***node-MPP***
c...    now gop everything together
      call pg_dgop(902,de,ncoord,'+')
c***   ***node-MPP***
_ENDIF
         call dcopy(ncoord,de,1,egrad,1)
c
c     ----- save 1e+2e-gradient -----
c
         irest = 0
         ist = 1
         jst = 1
         kst = 1
         lst = 1
         call wrtgrd(egrad,irest,ist,jst,kst,lst)
c
c --- punch out the gradients
c
         call blkgra(nat,egrad)
c
_IF(charmm)
         if(.not. onoatpr)then
_ENDIF
         if (oprint(31)) then
c
c  print out non-bq centres only
c
         maxpr = 0
         do while (maxpr .lt. nat)
            omore = .true.
            nbuff = 0
            do while (omore)
               if (nbuff .eq. 8 .or. maxpr .ge. nat) then
                  omore = .false.
               else
                  maxpr = maxpr + 1
                  if (zaname(maxpr)(1:2) .ne. 'bq')then
                     nbuff = nbuff + 1
                     ibuffi(nbuff) = maxpr
                  endif
               endif
            enddo
            if (nbuff .gt. 0) then
               write (iwr,6020)
               write (iwr,6050) (ibuffi(i),i=1,nbuff)
               write (iwr,6020)
               do n = 1 , 3
                  write (iwr,6060) ydnam(n),(f(n,ibuffi(i)),i=1,nbuff)
               enddo
            endif
         enddo
         else

         maxpr = 0
 20      minpr = maxpr + 1
         maxpr = maxpr + 8
         if (maxpr.gt.nat) maxpr = nat
         write (iwr,6020)
         write (iwr,6050) (i,i=minpr,maxpr)
         write (iwr,6020)
         do 30 n = 1 , 3
            write (iwr,6060) ydnam(n) , (f(n,i),i=minpr,maxpr)
 30      continue
         if (maxpr.lt.nat) go to 20
c
         endif
_IF(charmm)
         endif
_ENDIF

         if (ofokab) then
            call dabout(q,odebug(21),iwr)
            iochf(1) = iochf(1) + nat*ndenin*3*lensec(nx)
         end if
         if (ofock) then
            call hsymd(q,iwr)
            call clredx
         end if
_IFN(qmmm)
         if (ompir) then
            call secget(isect(57),57,iblok)
            do 60 i = 1 , 3
               do 50 j = 1 , 3
                  do 40 k = 1 , nat
                     dipd(i,j,k) = dipd(i,j,k) + dipi(i,j,k)
 40               continue
 50            continue
 60         continue
            call wrt3(dipd,lds(isect(57)),iblok,idaf)
         end if
_ELSE
         if (ompir) then
            call caserr('not available in qmmm version')
         endif
_ENDIF
         call revind
         call clredx
         cpu = cpulft(1)
         if (nprint.ne.-5) write (iwr,6040) cpu ,charwall()
         return
      else
         irest = 6
         istnu = istd
         jstnu = 1
         kstnu = 1
         lstnu = 1
_IF(parallel)
c***   ***node-MPP***
c...    now gop everything together
         ncoord = 3*nat
         call pg_dgop(902,de,ncoord,'+')
         call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
c***   ***node-MPP***
_ELSE
         call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
_ENDIF
         return
      end if
 6010 format (/' insufficient time to complete evaluation of 2-electron'
     +        ,' contribution to gradient'//' job dumped at ',f10.2,
_IF(parallel)
     +        ' seconds',a10,' wall'//' you can forget next batch ',4i5)
_ELSE
     +        ' seconds',a10,' wall'//' next batch ',4i5)
_ENDIF
 6020 format (/)
 6030 format (//10x,27('*')/10x,'*** warning ***'/10x,
     +        'this job must be restarted'/10x,27('*')/)
 6040 format (/' end of calculation of the energy gradient at ',f8.2,
     +        ' seconds',a10,' wall'/)
 6050 format (1x,'atom',8(6x,i3,6x))
 6060 format (3x,a3,8f15.7)
      end
      subroutine dform(x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
      implicit REAL  (a-h,o-z)
      logical unroll
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
INCLUDE(common/root)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/incrd)
INCLUDE(common/dshlno)
      dimension ix(*),iy(*),iz(*)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     1 yd(ncdim,*),zd(ncdim,*),g(*)
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
cvd$r assoc
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5 .and. ijkld.ge.16
      if (.not.unroll .or. ncontr.gt.5) then
         do 30 n = 1 , ijkld
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            do 20 nr = 1 , ncontr
               t1(nr) = xd(nr,mx)*y(nr,my)*z(nr,mz)
               t2(nr) = x(nr,mx)*yd(nr,my)*z(nr,mz)
               t3(nr) = x(nr,mx)*y(nr,my)*zd(nr,mz)
 20         continue
            g(n1+n) = g(n1+n) + dsum(ncontr,t1,1)
            g(n2+n) = g(n2+n) + dsum(ncontr,t2,1)
            g(n3+n) = g(n3+n) + dsum(ncontr,t3,1)
 30      continue
         return
      else
         go to (40,60,80,100,120) , ncontr
      end if
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 40   do 50 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz)
 50   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 60   do 70 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz)
 70   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 80   do 90 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz)
 90   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 100  do 110 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +             *y(4,my)*z(4,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +             *yd(4,my)*z(4,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +             *y(4,my)*zd(4,mz)
 110  continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 120  do 130 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +             *y(4,my)*z(4,mz) + xd(5,mx)*y(5,my)*z(5,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +             *yd(4,my)*z(4,mz) + x(5,mx)*yd(5,my)*z(5,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +             *y(4,my)*zd(4,mz) + x(5,mx)*y(5,my)*zd(5,mz)
 130  continue
      return
      end
      subroutine dforma(spij,spkl,noform,
     *                  x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
      implicit REAL  (a-h,o-z)
      logical noform,spij,spkl,unroll
c
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
      parameter (ncmax=65)
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
INCLUDE(common/root)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/incrd)
INCLUDE(common/dshlno)
      dimension ix(*),iy(*),iz(*)
      dimension noform(*)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     & yd(ncdim,*),zd(ncdim,*),g(*)
c
cvd$r assoc
c
      n = 0
      nn = 0
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5
      if (.not.unroll) then
c
         if (.not.spij) then
            do 40 i = 1 , ijd
               do 30 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     do 20 nr = 1 , ncontr
                        t1(nr) = ddkl(nr,k)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        t2(nr) = ddkl(nr,k)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        t3(nr) = ddkl(nr,k)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 20                  continue
                     g(n+n1) = dsum(ncontr,t1,1) + g(n+n1)
                     g(n+n2) = dsum(ncontr,t2,1) + g(n+n2)
                     g(n+n3) = dsum(ncontr,t3,1) + g(n+n3)
                  end if
 30            continue
 40         continue
         else if (.not.spkl) then
            do 70 i = 1 , ijd
               do 60 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 50 nr = 1 , ncontr
                        s1 = s1 + ddij(nr,i)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 50                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 60            continue
 70         continue
         else
            do 100 i = 1 , ijd
               do 90 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 80 nr = 1 , ncontr
                        s1 = s1 + (ddij(nr,i)*ddkl(nr,k))*xd(nr,mx)
     +                       *y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *y(nr,my)*zd(nr,mz)
 80                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 90            continue
 100        continue
         end if
         return
      else
         go to (110,140,170,200,230) , ncontr
      end if
 110  do 130 i = 1 , ijd
_IF1(x)c$dir scalar
         do 120 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz)
            end if
 120     continue
 130  continue
      return
 140  do 160 i = 1 , ijd
_IF1(x)c$dir scalar
         do 150 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz)
            end if
 150     continue
 160  continue
      return
 170  do 190 i = 1 , ijd
_IF1(x)c$dir scalar
         do 180 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz)
            end if
 180     continue
 190  continue
      return
 200  do 220 i = 1 , ijd
_IF1(x)c$dir scalar
         do 210 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz)
            end if
 210     continue
 220  continue
      return
 230  do 250 i = 1 , ijd
_IF1(x)c$dir scalar
         do 240 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz) + (ddij(5,i)*ddkl(5,k))
     +                   *xd(5,mx)*y(5,my)*z(5,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*yd(5,my)
     +                   *z(5,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*y(5,my)
     +                   *zd(5,mz)
            end if
 240     continue
 250  continue
      return
      end
      subroutine dgenrl(qq,iqq,noform,abmax)
c==========================================
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
      logical trduij,spij,spkl,trdukl
      logical double
      logical sptru,noform
INCLUDE(common/dshlnf)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +            aej(ncmax),aek(ncmax),ael(ncmax),
     +            aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +            dd(4*mxp2),ijden(225),ik(225),
     +            ijx(225),ijy(225),ijz(225),klx(225),kly(225),
     +            klz(225),
     +            dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/dshlno)
INCLUDE(common/root)
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     +   xc00(ncmax),ycp00(ncmax),yc00(ncmax),zcp00(ncmax),
     +   zc00(ncmax),f00(ncmax),
     +   dxij,dyij,dzij,dxkl,dykl,dzkl,
     +   in(12),kn(12),ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2
INCLUDE(common/dmisc)
INCLUDE(common/incrd)
      common/bufb/axak(mxp2),ayak(mxp2),azak(mxp2),axai(mxp2),
     1  ayai(mxp2),azai(mxp2),abv(mxp2),aandbv(mxp2),rhov(mxp2),
     2  xxv(mxp2),c1xv(mxp2),c2xv(mxp2),c3xv(mxp2),c4xv(mxp2),
     3  c1yv(mxp2),c2yv(mxp2),c3yv(mxp2),c4yv(mxp2),
     4  c1zv(mxp2),c2zv(mxp2),c3zv(mxp2),c4zv(mxp2),expev(mxp2)
INCLUDE(common/segm)
c
      dimension iqq(*),qq(*)
      dimension noform(*)
      data one/1.0d0/
c     data zero,pt5/0.0d0,0.5d0/
      data pi252/34.986836655250d0/
c
      if (ijkld.eq.0) return
      if (ijd.eq.1 .and. kld.eq.1) then
         call ssdss(qq)
      else
         ni = lit
         if (oskip(1)) ni = lit - 1
         nj = ljt
         if (oskip(2)) nj = ljt - 1
         nk = lkt
         if (oskip(3)) nk = lkt - 1
         nl = llt
         if (oskip(4)) nl = llt - 1
         kln2 = 1
         kln1 = nl + 1
         ijn2 = kln1*(nk+1)
         ijn1 = ijn2*(nj+1)
         inc1 = ijn1*(ni+1)
c     if(mod(inc1,4).eq.0)inc1=inc1+1
         if (ni.lt.nj) then
            is = ni
            ni = nj
            nj = is
            ij1 = ijn2
            ij2 = ijn1
            xc = xj
            yc = yj
            zc = zj
            dxij = xj - xi
            dyij = yj - yi
            dzij = zj - zi
         else
            ij1 = ijn1
            ij2 = ijn2
            xc = xi
            yc = yi
            zc = zi
            dxij = xi - xj
            dyij = yi - yj
            dzij = zi - zj
         end if
         if (nk.lt.nl) then
            is = nl
            nl = nk
            nk = is
            kl1 = kln2
            kl2 = kln1
            xd = xl
            yd = yl
            zd = zl
            dxkl = xl - xk
            dykl = yl - yk
            dzkl = zl - zk
         else
            xd = xk
            yd = yk
            zd = zk
            dxkl = xk - xl
            dykl = yk - yl
            dzkl = zk - zl
            kl1 = kln1
            kl2 = kln2
         end if
         nmax = ni + nj
         mmax = nk + nl
         max = nmax + 1
         do 20 i = 1 , max
            n = i - 1
            if (n.le.ni) in(i) = ij1*n + 1
            if (n.gt.ni) in(i) = ij1*ni + ij2*(n-ni) + 1
 20      continue
         max = mmax + 1
         do 30 k = 1 , max
            n = k - 1
            if (n.le.nk) kn(k) = kl1*n
            if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 30      continue
c
c     indexing
c
         call indexa(ijx,ijy,ijz,ijd,mini,maxi,minj,maxj,oiandj,ijn1,
     +              ijn2,1)
         call indexa(klx,kly,klz,kld,mink,maxk,minl,maxl,okandl,kln1,
     +              kln2,0)
         nn = 0
         ijkld = 0
         do 50 i = 1 , ijd
            do 40 k = 1 , kld
               nn = nn + 1
               if (.not.(noform(nn))) then
                  ijkld = ijkld + 1
                  iqq(ijkld+ixi-1) = ijx(i) + klx(k)
                  iqq(ijkld+iyi-1) = ijy(i) + kly(k)
                  iqq(ijkld+izi-1) = ijz(i) + klz(k)
               end if
 40         continue
 50      continue
c
         do 60 n = 1 , nij
            axak(n) = aa(n)*(x1(n)-xd)
            ayak(n) = aa(n)*(y1(n)-yd)
            azak(n) = aa(n)*(z1(n)-zd)
            axai(n) = aa(n)*(x1(n)-xc)
            ayai(n) = aa(n)*(y1(n)-yc)
            azai(n) = aa(n)*(z1(n)-zc)
 60      continue
c
c
         trduij = lit.ge.3 .or. ljt.ge.3
         trdukl = lkt.ge.3 .or. llt.ge.3
         spkl = (mink.eq.1 .and. maxk.eq.4) .or.
     +          (minl.eq.1 .and. maxl.eq.4)
         spij = (mini.eq.1 .and. maxi.eq.4) .or.
     +          (minj.eq.1 .and. maxj.eq.4)
         sptru = spkl .or. spij
c
c Now dgenrl uses preallocated memory, max size defined 
c in jkder only pointers ic2 upwards are allocated here
c vect factor of 32 is hardwired
c
c NB there is room above ic7 assuming npass=3, lendd,=lnddm
c
          ncmmm = 32
          ncmmm = (ncmmm/nroots)*nroots

cccccc old memory algorithm cccccccccccccc
c
c     integrals stored at (ic7+1)
c     subsidiary integrals from ic1 onwards
c
cc         ic1 = ic7 + 1 + npass*lendd*3
cc         ncmmm = (nmaxly-ic1-1)/(inc1*6)
cc         ncnnn = ncmax - 1
cc         if (ncmmm.gt.ncnnn) ncmmm = ncnnn
cc         ncmmm = (ncmmm/nroots)*nroots
cc         if (ncmmm.lt.nroots) call caserr('insufficient core in dgenrl')
cc
ccccccc end old memory algorithm cccccccccccc
c
         ic1 = ic7 + 1 + npass*lendd*3
         ic2 = inc1*ncmmm + ic1
         ic3 = inc1*ncmmm + ic2
         ic4 = inc1*ncmmm + ic3
         ic5 = inc1*ncmmm + ic4
         ic6 = inc1*ncmmm + ic5
         inc = npass*lendd*3
         call vclr(qq(ic7+1),1,inc)
         ncontr = 0
c
         maxlg = ngd
         do 170 kg = 1 , ngc
            ak = cgg(kg)
            brrk = ak*rrk
            akxk = ak*xk
            akyk = ak*yk
            akzk = ak*zk
            csk = csc(kg)*pi252
            cpk = cpc(kg)*pi252
            cdk = cdc(kg)*pi252
            cfk = cfc(kg)*pi252
            cgk = cgc(kg)*pi252
c
c     ----- l primitive
c
            if (okandl) maxlg = kg
            do 160 lg = 1 , maxlg
               al = dg(lg)
               b = ak + al
               binv = one/b
               bbrrk = al*brrk*binv
               if ((bbrrk+rsmall).le.tol1) then
                  exkl = dexp(-bbrrk)
                  csl = csd(lg)*binv
                  cpl = cpd(lg)*binv
                  cdl = cdd(lg)*binv
                  cfl = cfd(lg)*binv
                  cgl = cgd(lg)*binv
                  xb = (akxk+al*xl)*binv
                  yb = (akyk+al*yl)*binv
                  zb = (akzk+al*zl)*binv
                  bxbk = b*(xb-xd)
                  bybk = b*(yb-yd)
                  bzbk = b*(zb-zd)
                  bxbi = b*(xb-xc)
                  bybi = b*(yb-yc)
                  bzbi = b*(zb-zc)
c
c     ----- density factor
c
                  double = okandl .and. kg.ne.lg
                  call denfac(dkl,csk,cpk,cdk,cfk,cgk,
     +                        csl,cpl,cdl,cfl,cgl,
     +                        mink,maxk,minl,maxl,okandl,double)
                  dkld = dkl(1)
                  if (sptru) then
                     do 70 k = 1 , kld
                        dkl(k) = dkl(k)/dkld
 70                  continue
                  end if
                  dkld = dkld*exkl
                  if (dabs(dkld*abmax).ge.tol3) then
c
c     ----- pair of i,j primitives
c
                     do 80 n = 1 , nij
                        abv(n) = aa(n)*b
                        aandbv(n) = aa(n) + b
                        expev(n) = exij(n)/dsqrt(aa(n)+b)
                        rhov(n) = abv(n)/aandbv(n)
                        xxv(n) = rhov(n)
     +                           *((x1(n)-xb)**2+(y1(n)-yb)**2+(z1(n)
     +                           -zb)**2)
                        c1xv(n) = bxbk + axak(n)
                        c2xv(n) = bxbk*aa(n)
                        c3xv(n) = bxbi + axai(n)
                        c4xv(n) = b*axai(n)
                        c1yv(n) = bybk + ayak(n)
                        c2yv(n) = bybk*aa(n)
                        c3yv(n) = bybi + ayai(n)
                        c4yv(n) = b*ayai(n)
                        c1zv(n) = bzbk + azak(n)
                        c2zv(n) = bzbk*aa(n)
                        c3zv(n) = bzbi + azai(n)
                        c4zv(n) = b*azai(n)
 80                  continue
c
                     n = 0
                     nn = 0
                     jgmax = ngb
                     do 150 ig = 1 , nga
                        ai = ag(ig)
                        if (oiandj) jgmax = ig
                        do 140 jg = 1 , jgmax
                           n = n + 1
                           if ((bbrrk+r(n)).lt.tol2) then
                              aj = bg(jg)
                              dijd = dd(nn+1)
                              if (sptru) then
                                 dddd = one/dijd
                                 do 90 i = 1 , ijd
                                    dij(i) = dd(ijden(i)+nn)*dddd
 90                              continue
                              end if
                              expe = dkld*dijd*expev(n)
                              if (dabs(expe*abmax).ge.tol4) then
                                 pp = xxv(n)
c
c     ----- roots and weights for quadrature
c
                                 if (nroots.le.3) call rt123
                                 if (nroots.eq.4) call roots4
                                 if (nroots.eq.5) call roots5
                                 if (nroots.gt.5) call rootss
c
c     compute two-electron  integrals for each root
c
                                 nnn0 = ncontr
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                                 do 100 m = 1 , nroots
                                    ncontr = ncontr + 1
                                    u2 = u(m)*rhov(n)
                                    f00(ncontr) = expe*w(m)
                                    dum2 = 0.5d0/(abv(n)+u2*aandbv(n))
                                    dum = dum2 + dum2
                                    bp01(ncontr) = (aa(n)+u2)*dum2
                                    b00(ncontr) = u2*dum2
                                    b10(ncontr) = (b+u2)*dum2
                                    xcp00(ncontr) = (u2*c1xv(n)+c2xv(n))
     +                                 *dum
                                    xc00(ncontr) = (u2*c3xv(n)+c4xv(n))
     +                                 *dum
                                    ycp00(ncontr) = (u2*c1yv(n)+c2yv(n))
     +                                 *dum
                                    yc00(ncontr) = (u2*c3yv(n)+c4yv(n))
     +                                 *dum
                                    zcp00(ncontr) = (u2*c1zv(n)+c2zv(n))
     +                                 *dum
                                    zc00(ncontr) = (u2*c3zv(n)+c4zv(n))
     +                                 *dum
                                    aei(ncontr) = ai
                                    aej(ncontr) = aj
                                    aek(ncontr) = ak
                                    ael(ncontr) = al
 100                             continue
                                 if (sptru) then
                                    ncontr = nnn0
                                    do 130 m = 1 , nroots
                                       ncontr = ncontr + 1
                                       do 110 iii = 1 , ijd
                                         ddij(ncontr,iii) = dij(iii)
 110                                   continue
                                       do 120 iii = 1 , kld
                                         ddkl(ncontr,iii) = dkl(iii)
 120                                   continue
 130                                continue
                                 end if
c
c
c----------------------------------------------
c    defer assembly stage until loop lengths
c    are long enough to vectorise effectively
c----------------------------------------------
c
                                 if (ncontr.ge.ncmmm) then
c     ----- form (i,j//k,l) integrals
c
                                    call dxyz(qq(ic1),qq(ic2),qq(ic3),
     +                                 ncmmm)
                                    ioff = 0
                                    if (.not.(oskip(1))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aei,ljt,lkt,llt,lit,
     +                                    ijn2,kln1,kln2,ijn1,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(2))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aej,lit,lkt,llt,ljt,
     +                                    ijn1,kln1,kln2,ijn2,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(3))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aek,lit,ljt,llt,lkt,
     +                                    ijn1,ijn2,kln2,kln1,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(4))) then
                                       call subsd(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),ael,lit,ljt,lkt,llt,
     +                                    ijn1,ijn2,kln1,kln2,ncmmm)
                                       if (sptru) then
                                         call dforma(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                      qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                         call dform(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                    end if
                                    ncontr = 0
                                 end if
                              end if
                           end if
c
c
c     end of loops over primitives
c
                           nn = nn + 4
 140                    continue
 150                 continue
                  end if
               end if
 160        continue
 170     continue
c
c    tidy up any bits left unassembled
c
c
         if (ncontr.ne.0) then
            call dxyz(qq(ic1),qq(ic2),qq(ic3),ncmmm)
            ioff = 0
            if (.not.(oskip(1))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,
     +                    ijn1,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(2))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,
     +                    ijn2,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(3))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,
     +                    kln1,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(4))) then
               call subsd(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,
     +                    kln2,ncmmm)
               if (sptru) then
                  call dforma(spij,spkl,noform,qq(ic1),qq(ic2),qq(ic3),
     +                        qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi)
     +                        ,iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                       qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
            end if
            ncontr = 0
         end if
c
c ---------------------------------------------
c    fiddle about if first two centres are same
c ---------------------------------------------
c
         if (natomd(1).eq.natomd(2)) then
            inc = lendd*3
            ii = ic7 + inc
            do 180 i = 1 , inc
               qq(ic7+i) = qq(ic7+i) + qq(ii+i)
 180        continue
            iii = ii + inc
            do 190 i = 1 , inc
               qq(ii+i) = qq(iii+i)
 190        continue
            npass = npass - 1
            natomd(2) = natomd(3)
            natomd(3) = natomd(4)
            natomd(4) = 0
         end if
c
c    insert proper normalisation factors in d-functions
c    present (f- and g-functions also)
c
         if (trduij .or. trdukl) then
            if (onorm) then
               call dnorm(qq(ic7+1),noform)
            end if
         end if
      end if
c
c    multiply by density matrix elements
c
      call grdcon(qq,noform)
      return
      end
      subroutine dnorm(qq,noform)
      implicit REAL  (a-h,o-z)
      dimension qq(*)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
INCLUDE(common/dmisc)
INCLUDE(common/dshlno)
      logical noform
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aaa(9*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      dimension noform(*)
INCLUDE(common/incrd)
INCLUDE(common/picon)
      data one/1.0d0/
      n = 0
      max = maxj
      dum1 = one
      do 30 i = mini , maxi
         if (i.eq.8)  dum1 = root3
         if (i.eq.14) dum1 = root5
         if (i.eq.20) dum1 = dum1*root3
         if (i.eq.24) dum1 = root7
         if (i.eq.30) dum1 = dum1*root53
         if (i.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (oiandj) max = i
         do 20 j = minj , max
            if (j.eq.8)  dum2 = dum2*root3
            if (j.eq.14) dum2 = dum2*root5
            if (j.eq.20) dum2 = dum2*root3
            if (j.eq.24) dum2 = dum2*root7
            if (j.eq.30) dum2 = dum2*root53
            if (j.eq.33) dum2 = dum2*root3
            n = n + 1
            dij(n) = dum2
 20      continue
 30   continue
      n = 0
      dum1 = one
      max = maxl
      do 50 k = mink , maxk
         if (k.eq.8)  dum1 = root3
         if (k.eq.14) dum1 = root5
         if (k.eq.20) dum1 = dum1*root3
         if (k.eq.24) dum1 = root7
         if (k.eq.30) dum1 = dum1*root53
         if (k.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (okandl) max = k
         do 40 l = minl , max
            if (l.eq.8)  dum2 = dum2*root3
            if (l.eq.14) dum2 = dum2*root5
            if (l.eq.20) dum2 = dum2*root3
            if (l.eq.24) dum2 = dum2*root7
            if (l.eq.30) dum2 = dum2*root53
            if (l.eq.33) dum2 = dum2*root3
            n = n + 1
            dkl(n) = dum2
 40      continue
 50   continue
      nn = 0
      nnn = 0
      do 80 i = 1 , ijd
         d1 = dij(i)
         do 70 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               nnn = nnn + 1
               ioff = 0
               do 60 npp = 1 , npass*3
                  qq(nnn+ioff) = qq(nnn+ioff)*d1*dkl(k)
                  ioff = ioff + lendd
 60            continue
            end if
 70      continue
 80   continue
      return
      end
      subroutine dpdab1(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit REAL  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/incrd)
INCLUDE(common/atmblk)
      common/blk1/gin(510),mword,nlll1,kwor1,kwor11
_IFN1(iv)      common/sortpk/labs(1360)
_IF1(iv)       common/sortpk/i205(340),j205(340),k205(340),l205(340)
c
      call vclr(abdens,1,lendd)
      ijdiff = ii.ne.jj
      if (kwor1.eq.0) return
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor1)
         iword = 0
_IFN1(iv)         call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)         call upak8v(gin(num2e+1),i205)
         if (kwor1.le.0) go to 30
      end if
      iword = iword + 1
_IFN1(iv)      m = (iword+iword) + (iword+iword)
_IFN1(iv)      i = labs(m-3)
_IFN1(iv)      j = labs(m-2)
_IFN1(iv)      k = labs(m-1)
_IFN1(iv)      l = labs(m)
_IF1(iv)c
_IF1(iv)      i = i205(iword)
_IF1(iv)      j = j205(iword)
_IF1(iv)      k = k205(iword)
_IF1(iv)      l = l205(iword)
_IF1(iv)c
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dpdab2(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit REAL  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/incrd)
INCLUDE(common/atmblk)
_IFN1(iv)      common/three/labs(1360)
_IF1(iv)       common/three/i205(340),j205(340),k205(340),l205(340)
      common/bufc/gin(510),mword,nlll2,kwor2,kwor22
c
      call vclr(abdens,1,lendd)
      if (kwor2.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor2)
         iword = 0
_IF1(iv)         call upak8v(gin(num2e+1),i205)
_IFN1(iv)         call unpack(gin(num2e+1),lab816,labs,numlab)
         if (kwor2.le.0) go to 30
      end if
      iword = iword + 1
_IFN1(iv)      m = (iword+iword) + (iword+iword)
_IFN1(iv)      i = labs(m-3)
_IFN1(iv)      j = labs(m-2)
_IFN1(iv)      k = labs(m-1)
_IFN1(iv)      l = labs(m)
_IF1(iv)      i = i205(iword)
_IF1(iv)      j = j205(iword)
_IF1(iv)      k = k205(iword)
_IF1(iv)      l = l205(iword)
_IF1(iv)c
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dpdab3(abdens,ii,jj,kk,ll,q4,ifil2d,iword)
      implicit REAL  (a-h,o-z)
      dimension abdens(*)
      logical ijdiff,kldiff,ikdiff
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/incrd)
INCLUDE(common/atmblk)
_IFN1(iv)      common/lsort/labs(1360)
_IF1(iv)       common/lsort/i205(340),j205(340),k205(340),l205(340)
      common/bufd/gin(510),mword,nlll3,kwor3,kwor33
c
      call vclr(abdens,1,lendd)
      if (kwor3.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kwor3)
         iword = 0
_IFN1(iv)         call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)         call upak8v(gin(num2e+1),i205)
         if (kwor3.le.0) go to 30
      end if
      iword = iword + 1
_IFN1(iv)      m = (iword+iword) + (iword+iword)
_IFN1(iv)      i = labs(m-3)
_IFN1(iv)      j = labs(m-2)
_IFN1(iv)      k = labs(m-1)
_IFN1(iv)      l = labs(m)
_IF1(iv)c
_IF1(iv)      i = i205(iword)
_IF1(iv)      j = j205(iword)
_IF1(iv)      k = k205(iword)
_IF1(iv)      l = l205(iword)
_IF1(iv)c
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
      subroutine dprim
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/dmisc)
INCLUDE(common/dshlnf)
INCLUDE(common/dshlno)
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +            aej(ncmax),aek(ncmax),ael(ncmax),
     +            a(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),dd(4*mxp2),
     +            ijden(225)
      data one/1.0d0/
      max = maxj
      n = 0
      nn = 0
      do 70 i = mini , maxi
         go to (20,20,30,30,
     +          20,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          30,30,30,30,30) , i
 20      nm = nn
 30      nn = nm
         if (oiandj) max = i
         do 60 j = minj , max
            go to (40,40,50,50,
     +             40,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             50,50,50,50,50) , j
 40         nn = nn + 1
 50         n = n + 1
            ijden(n) = nn
 60      continue
 70   continue
c     ----- i primitive
      nij = 0
      jbmax = ngb
      do 230 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*xi
         ayi = ai*yi
         azi = ai*zi
         csi = csa(ia)
         cpi = cpa(ia)
         cdi = cda(ia)
         cfi = cfa(ia)
         cgi = cga(ia) 
c     ----- j primitive
         if (oiandj) jbmax = ia
         do 220 jb = 1 , jbmax
            aj = bg(jb)
            aa = ai + aj
            ainv = one/aa
            dum = aj*arri*ainv
            csj = csb(jb)*ainv
            cpj = cpb(jb)*ainv
            cdj = cdb(jb)*ainv
            cfj = cfb(jb)*ainv
            cgj = cgb(jb)*ainv
            nm = (nij+nij) + (nij+nij)
            nn = nm
            nij = nij + 1
            r(nij) = dum
            a(nij) = aa
            x1(nij) = (axi+aj*xj)*ainv
            y1(nij) = (ayi+aj*yj)*ainv
            z1(nij) = (azi+aj*zj)*ainv
c     ----- density factor
            do 190 i = mini , maxi
               if (oiandj) max = i
               go to (80,90,190,190,
     +               100,190,190,190,190,190,
     +               110,190,190,190,190,190,190,190,190,190,
     +               115,190,190,190,190,190,190,190,190,190,
     +               190,190,190,190,190) , i
 80            dum1 = csi
               go to 120
 90            dum1 = cpi
               go to 120
 100           dum1 = cdi
               go to 120
 110           dum1 = cfi
               go to 120
 115           dum1 = cgi
 120           do 180 j = minj , max
                  go to (130,140,180,180,
     +                   150,180,180,180,180,180,
     +                   160,180,180,180,180,180,180,180,180,180,
     +                   165,180,180,180,180,180,180,180,180,180,
     +                   180,180,180,180,180) , j
 130              dum2 = dum1*csj
                  go to 170
 140              dum2 = dum1*cpj
                  go to 170
 150              dum2 = dum1*cdj
                  go to 170
 160              dum2 = dum1*cfj
                  go to 170
 165              dum2 = dum1*cgj
 170              nn = nn + 1
                  dd(nn) = dum2
 180           continue
 190        continue
            if (.not.oiandj) go to 220
            if (ia.eq.jb) go to 220
            go to (210,200,210,210,210) , lit
 200        if (mini.ne.2) then
               dd(nm+2) = dd(nm+2) + csi*cpj
               dd(nm+3) = dd(nm+3) + dd(nm+3)
            end if
 210        dd(nm+1) = dd(nm+1) + dd(nm+1)
 220     continue
 230  continue
      if (nij.eq.0) return
      rsmall = r(1)
      do 240 n = 1 , nij
         exij(n) = dexp(-r(n))
 240  continue
      do 250 n = 1 , nij
         if (rsmall.gt.r(n)) rsmall = r(n)
 250  continue
      if (rsmall.ge.tol1) nij = 0
      return
      end
      subroutine dshell(nelec,ish,jsh,ksh,lsh)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/dmisc)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/root)
INCLUDE(common/dshlno)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     +,aek(ncmax),ael(ncmax),
     + aaa(9*mxp2),ijden(225),
     + ik(225),ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/incrd)
INCLUDE(common/dshlnf)
c
      if (nelec.eq.2) then
         okandl = ksh.eq.lsh
         osame = ish.eq.ksh .and. jsh.eq.lsh
         k = katom(ksh)
         xk = c(1,k)
         yk = c(2,k)
         zk = c(3,k)
         k1 = kstart(ksh)
         k2 = k1 + kng(ksh) - 1
         lkt = ktype(ksh)
         mink = kmin(ksh)
         maxk = kmax(ksh)
         lock = kloc(ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex(k)
            csc(ngc) = cs(k)
            cpc(ngc) = cp(k)
            cdc(ngc) = cd(k)
            cfc(ngc) = cf(k)
            cgc(ngc) = cg(k)
 20      continue
         l = katom(lsh)
         xl = c(1,l)
         yl = c(2,l)
         zl = c(3,l)
         l1 = kstart(lsh)
         l2 = l1 + kng(lsh) - 1
         llt = ktype(lsh)
         minl = kmin(lsh)
         maxl = kmax(lsh)
         locl = kloc(lsh) - minl
         ngd = 0
         do 30 l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex(l)
            csd(ngd) = cs(l)
            cpd(ngd) = cp(l)
            cdd(ngd) = cd(l)
            cfd(ngd) = cf(l)
            cgd(ngd) = cg(l)
 30      continue
         nroots = (lit+ljt+lkt+llt-1)/2
         rrk = ((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
c
c     determine various offsets and indexing arrays
c
         inc2 = 1
         inc3 = inc2*(maxl-minl+1)
         inc4 = inc3*(maxk-mink+1)
         inc5 = inc4*(maxj-minj+1)
         lendd = inc5*(maxi-mini+1)
         if (mod(lendd,4).eq.0) lendd = lendd + 1
         ijd = 0
         max = maxj
         do 50 i = mini , maxi
            if (oiandj) max = i
            ittt = inc5*(i-mini) + 1
            do 40 j = minj , max
               ijd = ijd + 1
               ijgt(ijd) = ittt
               ittt = ittt + inc4
 40         continue
 50      continue
         kld = 0
         max = maxl
         do 70 k = mink , maxk
            if (okandl) max = k
            ittt = inc3*(k-mink)
            do 60 l = minl , max
               kld = kld + 1
               klgt(kld) = ittt
               ittt = ittt + inc2
 60         continue
 70      continue
         ijkld = ijd*kld
         ixi = lendd + 1
         iyi = ixi + lendd
         izi = iyi + lendd
         ioff = 0
         return
      else
         oiandj = ish.eq.jsh
         i = katom(ish)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(ish)
         i2 = i1 + kng(ish) - 1
         lit = ktype(ish)
         mini = kmin(ish)
         maxi = kmax(ish)
         loci = kloc(ish) - mini
         nga = 0
         do 80 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex(i)
            csa(nga) = cs(i)
            cpa(nga) = cp(i)
            cda(nga) = cd(i)
            cfa(nga) = cf(i)
            cga(nga) = cg(i)
 80      continue
         j = katom(jsh)
         xj = c(1,j)
         yj = c(2,j)
         zj = c(3,j)
         j1 = kstart(jsh)
         j2 = j1 + kng(jsh) - 1
         ljt = ktype(jsh)
         minj = kmin(jsh)
         maxj = kmax(jsh)
         locj = kloc(jsh) - minj
         ngb = 0
         do 90 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex(j)
            csb(ngb) = cs(j)
            cpb(ngb) = cp(j)
            cdb(ngb) = cd(j)
            cfb(ngb) = cf(j)
            cgb(ngb) = cg(j)
 90      continue
         rri = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
         return
      end if
      end
      subroutine dxyz(x,y,z,ncdim)
      implicit REAL  (a-h,o-z)
      parameter (ncmax=65)
      dimension x(*),y(*),z(*)
      logical n0,n1,m0,m1
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     1   xc00(ncmax),
     2ycp00(ncmax),yc00(ncmax),zcp00(ncmax),zc00(ncmax),f00(ncmax),dxij,
     3dyij,dzij,dxkl,dykl,dzkl,iorg(12),korg(12),
     4nimax,njmax,nkmax,nlmax,nmax,mmax,ij1x,ij2x,kl1x,kl2x
      common/small/ca(ncmax),cb(ncmax)
INCLUDE(common/dshlno)
      dimension i(12),k(12)
c
      data zero,one /0.0d+00,1.0d+00/
c
      do 20 n = 1 , nmax + 1
         i(n) = (iorg(n)-1)*ncdim + 1
 20   continue
      do 30 n = 1 , mmax + 1
         k(n) = korg(n)*ncdim
 30   continue
      ij1 = ij1x*ncdim
      ij2 = ij2x*ncdim
      kl1 = kl1x*ncdim
      kl2 = kl2x*ncdim
      ink = 1
c
      n0 = nmax.eq.0
      n1 = nmax.le.1
      m0 = mmax.eq.0
      m1 = mmax.le.1
      if (n0) then
         i1 = i(1)
         ia = 0
         do 40 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 40      continue
         if (m0) go to 670
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
         do 50 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 50      continue
         if (m1) go to 670
c     ----- i(0,nk) -----
c
         do 60 nc = 1 , ncontr
            ca(nc) = zero
 60      continue
         i3 = i1
         i4 = i1 + k2
         do 90 nk = 2 , mmax
c
            do 70 nc = 1 , ncontr
               ca(nc) = ca(nc) + bp01(nc)
 70         continue
            i5 = i1 + k(nk+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 80 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
               ia = ia + ink
 80         continue
            i3 = i4
            i4 = i5
 90      continue
         if (nlmax.eq.0) go to 670
c     ----- i(0,0,nk,nl) -----
         i5 = k(mmax+1)
         min = nkmax
         go to 610
      else if (m0) then
         i1 = i(1)
         ia = 0
         do 100 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 100     continue
         if (n0) go to 600
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
         do 110 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 110     continue
         if (n1) go to 600
c     ----- i(ni,0) -----
c
         do 120 nc = 1 , ncontr
            ca(nc) = zero
 120     continue
         i3 = i1
         i4 = i2
         do 150 ni = 2 , nmax
            do 130 nc = 1 , ncontr
               ca(nc) = ca(nc) + b10(nc)
 130        continue
            i5 = i(ni+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 140 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
               ia = ia + ink
 140        continue
            i3 = i4
            i4 = i5
 150     continue
         if (njmax.eq.0) go to 600
c     ----- i(ni,nj,0,0) -----
         i5 = i(nmax+1)
         min = nimax
         go to 540
      else
c     ----- i(0,0) -----
         i1 = i(1)
         ia = 0
c
         do 160 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 160     continue
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
c
         do 170 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 170     continue
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
c
         do 180 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 180     continue
c     ----- i(1,1) -----
         i3 = i2 + k2
         ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 190 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)*x(i2+ia) + b00(nc)
            y(i3+ia) = ycp00(nc)*y(i2+ia) + b00(nc)
            z(i3+ia) = zcp00(nc)*z(i2+ia) + b00(nc)*f00(nc)
            ia = ia + ink
 190     continue
         if (.not.(n1)) then
            do 200 nc = 1 , ncontr
               cb(nc) = b00(nc)
               ca(nc) = zero
 200        continue
            i3 = i1
            i4 = i2
            do 250 n = 2 , nmax
c     ----- i(n,0) -----
               i5 = i(n+1)
               do 210 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b10(nc)
 210           continue
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 220 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
                  ia = ia + ink
 220           continue
               do 230 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 230           continue
c     ----- i(n,1) -----
               i3 = i5 + k2
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 240 nc = 1 , ncontr
                  x(i3+ia) = xcp00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = ycp00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zcp00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 240           continue
               i3 = i4
               i4 = i5
 250        continue
         end if
         if (.not.(m1)) then
            do 260 nc = 1 , ncontr
               ca(nc) = zero
               cb(nc) = b00(nc)
 260        continue
            i3 = i1
            i4 = i1 + k2
            do 310 m = 2 , mmax
               do 270 nc = 1 , ncontr
                  ca(nc) = ca(nc) + bp01(nc)
 270           continue
c     ----- i(0,m) -----
               i5 = i1 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 280 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
                  ia = ia + ink
 280           continue
               do 290 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 290           continue
c     ----- i(1,m) -----
               i3 = i2 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 300 nc = 1 , ncontr
                  x(i3+ia) = xc00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = yc00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zc00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 300           continue
               i3 = i4
               i4 = i5
 310        continue
         end if
         if (.not.(n1 .or. m1)) then
c     ----- i(n,m) -----
c
            do 320 nc = 1 , ncontr
               ca(nc) = b00(nc)
 320        continue
            k3 = k2
            do 370 m = 2 , mmax
               k4 = k(m+1)
               do 330 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b00(nc)
                  cb(nc) = b10(nc)
 330           continue
               i3 = i1
               i4 = i2
               do 360 n = 2 , nmax
                  i5 = i(n+1)
                  ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 340 nc = 1 , ncontr
                     x(i5+k4+ia) = cb(nc)*x(i3+k4+ia) + xc00(nc)
     +                             *x(i4+k4+ia) + ca(nc)*x(i4+k3+ia)
                     y(i5+k4+ia) = cb(nc)*y(i3+k4+ia) + yc00(nc)
     +                             *y(i4+k4+ia) + ca(nc)*y(i4+k3+ia)
                     z(i5+k4+ia) = cb(nc)*z(i3+k4+ia) + zc00(nc)
     +                             *z(i4+k4+ia) + ca(nc)*z(i4+k3+ia)
                     ia = ia + ink
 340              continue
                  do 350 nc = 1 , ncontr
                     cb(nc) = cb(nc) + b10(nc)
 350              continue
                  i3 = i4
                  i4 = i5
 360           continue
               k3 = k4
 370        continue
         end if
         if (njmax.eq.0) go to 450
c     ----- i(ni,nj,m) -----
         m = 0
         i5 = i(nmax+1)
      end if
 380  min = nimax
      km = k(m+1)
 390  n = nmax
      i3 = i5 + km
 400  i4 = i(n) + km
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 410 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 410  continue
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 400
      min = min + 1
      if (min.lt.nmax) go to 390
      if (nimax.ne.0) then
         i3 = ij2 + km + i1
         do 440 nj = 1 , njmax
            i4 = i3
            do 430 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 420 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 420           continue
               i4 = i4 + ij1
 430        continue
            i3 = i3 + ij2
 440     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 380
 450  if (nlmax.eq.0) go to 530
c
c     ----- i(ni,nj,nk,nl) -----
c
      i5 = k(mmax+1)
      ia = i1
      ni = 0
 460  nj = 0
      ib = ia
      min = nkmax
 470  m = mmax
      i3 = ib + i5
 480  i4 = ib + k(m)
      ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 490 nc = 1 , ncontr
         x(i3+ic) = x(i3+ic) + dxkl*x(i4+ic)
         y(i3+ic) = y(i3+ic) + dykl*y(i4+ic)
         z(i3+ic) = z(i3+ic) + dzkl*z(i4+ic)
         ic = ic + ink
 490  continue
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 480
      min = min + 1
      if (min.lt.mmax) go to 470
      if (nkmax.ne.0) then
         i3 = ib + kl2
         do 520 nl = 1 , nlmax
            i4 = i3
            do 510 nk = 1 , nkmax
               ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 500 nc = 1 , ncontr
                  x(i4+ic) = x(i4+ic+kl1-kl2) + dxkl*x(i4+ic-kl2)
                  y(i4+ic) = y(i4+ic+kl1-kl2) + dykl*y(i4+ic-kl2)
                  z(i4+ic) = z(i4+ic+kl1-kl2) + dzkl*z(i4+ic-kl2)
                  ic = ic + ink
 500           continue
               i4 = i4 + kl1
 510        continue
            i3 = i3 + kl2
 520     continue
      end if
      nj = nj + 1
      ib = ib + ij2
      if (nj.le.njmax) then
         min = nkmax
         go to 470
      else
         ni = ni + 1
         ia = ia + ij1
         if (ni.le.nimax) go to 460
      end if
 530  return
 540  ni = nmax
      i3 = i5
 550  i4 = i(ni)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 560 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 560  continue
      i3 = i4
      ni = ni - 1
      if (ni.gt.min) go to 550
      min = min + 1
      if (min.lt.nmax) go to 540
      if (nimax.ne.0) then
         i3 = ij2 + i1
         do 590 nj = 1 , njmax
            i4 = i3
            do 580 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 570 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 570           continue
               i4 = i4 + ij1
 580        continue
            i3 = i3 + ij2
 590     continue
      end if
 600  return
 610  nk = mmax
      i3 = i1 + i5
 620  i4 = i1 + k(nk)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 630 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxkl*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dykl*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzkl*z(i4+ia)
         ia = ia + ink
 630  continue
      i3 = i4
      nk = nk - 1
      if (nk.gt.min) go to 620
      min = min + 1
      if (min.lt.mmax) go to 610
      if (nkmax.ne.0) then
         i3 = i1 + kl2
         do 660 nl = 1 , nlmax
            i4 = i3
            do 650 nk = 1 , nkmax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 640 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+kl1-kl2) + dxkl*x(i4+ia-kl2)
                  y(i4+ia) = y(i4+ia+kl1-kl2) + dykl*y(i4+ia-kl2)
                  z(i4+ia) = z(i4+ia+kl1-kl2) + dzkl*z(i4+ia-kl2)
                  ia = ia + ink
 640           continue
               i4 = i4 + kl1
 650        continue
            i3 = i3 + kl2
 660     continue
      end if
 670  return
      end
      subroutine fockd2(qq,noform)
c============================================
c     derivative fock operators, j and k matrices
c     etc. for use in later chf equations
c================================================
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
c     derivative fock operators
c     additions to write integral derivatives
c
      logical noform
INCLUDE(common/dshlno)
INCLUDE(common/dmisc)
INCLUDE(common/specal)
INCLUDE(common/mapper)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  xij(225),xkl(225),ijgt(225),klgt(225)
c
INCLUDE(common/nshel)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/restri)
INCLUDE(common/cndx41)
INCLUDE(common/infoa)
INCLUDE(common/incrd)
_IF(ccpdft)
INCLUDE(common/ccpdft.hf77)
_ENDIF
      dimension qq(*),noform(*)
      data four /4.0d0/
      data small/1.0d-14/
c
c
      if (mp2w) then
         ioa = ic7 + lendd*npass*3
         call vclr(qq(ioa+1),1,3*lendd)
         do 30 i = 1 , npass
            iob = ic7 + (i-1)*lendd*3
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 20 k = 1 , lendd*3
               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
 20         continue
 30      continue
      end if
      idb = ida + nx
      np1 = npass + 1
      nw3 = nx*3
      lenbig = nat*nw3
      len3 = lendd*3
      ntr = ifok + (natomd(np1)-1)*nw3
c
_IF(ccpdft)
      hf_wght = CD_HF_exchange_weight()
_ENDIF
c
      if (ogrhf) then
c
c         general case
c
         nn = 0
         nnn = 0
         jmax = maxj
         do 100 i = mini , maxi
            i1 = loci + i
            ia1 = iky(i1)
            if (oiandj) jmax = i
            do 90 j = minj , jmax
               i2 = locj + j
               ia2 = iky(i2)
               mij = ia1 + i2
               lmax = maxl
               do 80 k = mink , maxk
                  i3 = lock + k
                  ia3 = iky(i3)
                  mik = ia1 + i3
                  if (i3.gt.i1) mik = ia3 + i1
                  mjk = ia2 + i3
                  if (i3.gt.i2) mjk = ia3 + i2
                  if (okandl) lmax = k
                  do 70 l = minl , lmax
                     i4 = locl + l
                     ia4 = iky(i4)
                     mkl = ia3 + i4
                     mil = ia1 + i4
                     if (i4.gt.i1) mil = ia4 + i1
                     mjl = ia2 + i4
                     if (i4.gt.i2) mjl = ia4 + i2
                     nnn = nnn + 1
                     if (.not.(noform(nnn))) then
                        nn = nn + 1
c
                        do 60 npp = 1 , npass
                           ifa = ifok + (natomd(npp)-1)*nw3
                           itrf = ntr
                           ioa = ic7 + (npp-1)*len3
c
                           do 50 icomp = 1 , 3
                              g1 = qq(nn+ioa)
c
                              if (dabs(g1).ge.small) then
                                 g1 = g1 + g1
                                 g14 = g1 + g1
                                 itj = itrf
                                 ifj = ifa
                                 ifd = ida
                                 do 40 ns = 1 , ndens
                                    ifk = ifj + lenbig
                                    itk = itj + lenbig
                                    gdkl = g14*qq(ifd+mkl)
                                    qq(ifj+mij) = qq(ifj+mij) + gdkl
                                    gdij = g14*qq(ifd+mij)
                                    qq(ifj+mkl) = qq(ifj+mkl) + gdij
                                    gdjl = g1*qq(ifd+mjl)
                                    qq(ifk+mik) = qq(ifk+mik) + gdjl
                                    gdjk = g1*qq(ifd+mjk)
                                    qq(ifk+mil) = qq(ifk+mil) + gdjk
                                    gdil = g1*qq(ifd+mil)
                                    qq(ifk+mjk) = qq(ifk+mjk) + gdil
                                    gdik = g1*qq(ifd+mik)
                                    qq(ifk+mjl) = qq(ifk+mjl) + gdik
c     use translational invariance to get missing term
                                    qq(itj+mij) = qq(itj+mij) - gdkl
                                    qq(itj+mkl) = qq(itj+mkl) - gdij
                                    qq(itk+mik) = qq(itk+mik) - gdjl
                                    qq(itk+mil) = qq(itk+mil) - gdjk
                                    qq(itk+mjk) = qq(itk+mjk) - gdil
                                    qq(itk+mjl) = qq(itk+mjl) - gdik
                                    ifj = ifj + lenbig + lenbig
                                    itj = itj + lenbig + lenbig
                                    ifd = ifd + nx
 40                              continue
                              end if
                              ifa = ifa + nx
                              itrf = itrf + nx
                              ioa = ioa + lendd
 50                        continue
 60                     continue
                     end if
 70               continue
 80            continue
 90         continue
 100     continue
      else
         nn = 0
         nnn = 0
         jmax = maxj
         do 160 i = mini , maxi
            i1 = loci + i
            ia1 = iky(i1)
            if (oiandj) jmax = i
            do 150 j = minj , jmax
               i2 = locj + j
               ia2 = iky(i2)
               mij = ia1 + i2
               dij = four*qq(ida+mij)
               lmax = maxl
               do 140 k = mink , maxk
                  i3 = lock + k
                  ia3 = iky(i3)
                  mik = ia1 + i3
                  if (i3.gt.i1) mik = ia3 + i1
                  mjk = ia2 + i3
                  if (i3.gt.i2) mjk = ia3 + i2
_IF(ccpdft)
                  dik = qq(ida+mik)*hf_wght
                  djk = qq(ida+mjk)*hf_wght
_ELSE
                  dik = qq(ida+mik)
                  djk = qq(ida+mjk)
_ENDIF
                  if (okandl) lmax = k
                  do 130 l = minl , lmax
                     i4 = locl + l
                     ia4 = iky(i4)
                     mkl = ia3 + i4
                     mil = ia1 + i4
                     if (i4.gt.i1) mil = ia4 + i1
                     mjl = ia2 + i4
                     if (i4.gt.i2) mjl = ia4 + i2
                     nnn = nnn + 1
                     if (.not.(noform(nnn))) then
                        nn = nn + 1
                        dkl = four*qq(ida+mkl)
_IF(ccpdft)
                        dil = qq(ida+mil)*hf_wght
                        djl = qq(ida+mjl)*hf_wght
_ELSE
                        dil = qq(ida+mil)
                        djl = qq(ida+mjl)
_ENDIF
c
c
                        do 120 npp = 1 , npass
                           ifa = ifok + (natomd(npp)-1)*nw3
                           itrf = ntr
                           ioa = ic7 + (npp-1)*len3
                           do 110 icomp = 1 , 3
                              g1 = qq(nn+ioa)
                              if (dabs(g1).ge.small) then
                                 gdkl = g1*dkl
                                 qq(ifa+mij) = qq(ifa+mij) + gdkl
                                 gdij = g1*dij
                                 qq(ifa+mkl) = qq(ifa+mkl) + gdij
                                 gdjl = g1*djl
                                 qq(ifa+mik) = qq(ifa+mik) - gdjl
                                 gdjk = g1*djk
                                 qq(ifa+mil) = qq(ifa+mil) - gdjk
                                 gdil = g1*dil
                                 qq(ifa+mjk) = qq(ifa+mjk) - gdil
                                 gdik = g1*dik
                                 qq(ifa+mjl) = qq(ifa+mjl) - gdik
c     use translational invariance to get missing term
                                 qq(itrf+mij) = qq(itrf+mij) - gdkl
                                 qq(itrf+mkl) = qq(itrf+mkl) - gdij
                                 qq(itrf+mik) = qq(itrf+mik) + gdjl
                                 qq(itrf+mil) = qq(itrf+mil) + gdjk
                                 qq(itrf+mjk) = qq(itrf+mjk) + gdil
                                 qq(itrf+mjl) = qq(itrf+mjl) + gdik
                                 if (.not.(oclos)) then
                                    ika = ifa + lenbig
                                    itrk = itrf + lenbig
                                    gdjl = g1*qq(idb+mjl)
                                    qq(ika+mik) = qq(ika+mik) + gdjl
                                    gdil = g1*qq(idb+mil)
                                    qq(ika+mjk) = qq(ika+mjk) + gdil
                                    gdjk = g1*qq(idb+mjk)
                                    qq(ika+mil) = qq(ika+mil) + gdjk
                                    gdik = g1*qq(idb+mik)
                                    qq(ika+mjl) = qq(ika+mjl) + gdik
                                    qq(itrk+mik) = qq(itrk+mik) - gdjl
                                    qq(itrk+mjl) = qq(itrk+mjl) - gdik
                                    qq(itrk+mjk) = qq(itrk+mjk) - gdil
                                    qq(itrk+mil) = qq(itrk+mil) - gdjk
                                 end if
                              end if
                              ifa = ifa + nx
                              itrf = itrf + nx
                              ioa = ioa + lendd
 110                       continue
 120                    continue
                     end if
 130              continue
 140           continue
 150        continue
 160     continue
         if (mp2w) call derwrt(qq(ic7+1))
         return
      end if
      end
      subroutine fokabd(qq,noform)
c
c----------------------------------------------------
c    construct fock matrices from derivative integrals
c    and derivative density matrices.
c    used in polarizability derivative calculations
c-----------------------------------------------------
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
      logical noform
INCLUDE(common/dshlno)
INCLUDE(common/dmisc)
INCLUDE(common/mapper)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +   aej(ncmax),aek(ncmax),ael(ncmax),
     +   aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +   dd(4*mxp2),ijden(225),ik(225),
     +   ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +   dij(225),dkl(225),ijgt(225),klgt(225)
c
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/specal)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/restri)
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
INCLUDE(common/cndx41)
INCLUDE(common/incrd)
      dimension qq(*),noform(*)
      data small/1.0d-14/
c
      np1 = npass + 1
      nw3 = nx*3
      lenbig = nat*nw3
      len3 = lendd*3
      ntr = ifok + (natomd(np1)-1)*nw3
c
c       ndens density matrices for each nuclear displacement
c
      nn = 0
      nnn = 0
      jmax = maxj
      do 80 i = mini , maxi
         i1 = loci + i
         ia1 = iky(i1)
         if (oiandj) jmax = i
         do 70 j = minj , jmax
            i2 = locj + j
            ia2 = iky(i2)
            mij = ia1 + i2
            lmax = maxl
            do 60 k = mink , maxk
               i3 = lock + k
               ia3 = iky(i3)
               mik = ia1 + i3
               if (i3.gt.i1) mik = ia3 + i1
               mjk = ia2 + i3
               if (i3.gt.i2) mjk = ia3 + i2
               if (okandl) lmax = k
               do 50 l = minl , lmax
                  i4 = locl + l
                  ia4 = iky(i4)
                  mkl = ia3 + i4
                  mil = ia1 + i4
                  if (i4.gt.i1) mil = ia4 + i1
                  mjl = ia2 + i4
                  if (i4.gt.i2) mjl = ia4 + i2
                  nnn = nnn + 1
                  if (.not.(noform(nnn))) then
                     nn = nn + 1
c
c      loops over nuclei....
c
                     do 40 npp = 1 , npass
                        ifa = ifok + (natomd(npp)-1)*nw3
                        itrf = ntr
                        ioa = ic7 + (npp-1)*len3
                        do 30 icomp = 1 , 3
                           g1 = qq(nn+ioa)
c
                           if (dabs(g1).ge.small) then
                              g14 = (g1+g1) + (g1+g1)
                              itf = itrf
                              iff = ifa
                              ifd = ida + nx
                              do 20 ns = 1 , ndenin
                                 gdkl = g14*qq(ifd+mkl)
                                 qq(iff+mij) = qq(iff+mij) + gdkl
                                 gdij = g14*qq(ifd+mij)
                                 qq(iff+mkl) = qq(iff+mkl) + gdij
                                 gdjl = g1*qq(ifd+mjl)
                                 qq(iff+mik) = qq(iff+mik) - gdjl
                                 gdjk = g1*qq(ifd+mjk)
                                 qq(iff+mil) = qq(iff+mil) - gdjk
                                 gdil = g1*qq(ifd+mil)
                                 qq(iff+mjk) = qq(iff+mjk) - gdil
                                 gdik = g1*qq(ifd+mik)
                                 qq(iff+mjl) = qq(iff+mjl) - gdik
c     use translational invariance to get missing term
                                 qq(itf+mij) = qq(itf+mij) - gdkl
                                 qq(itf+mkl) = qq(itf+mkl) - gdij
                                 qq(itf+mik) = qq(itf+mik) + gdjl
                                 qq(itf+mil) = qq(itf+mil) + gdjk
                                 qq(itf+mjk) = qq(itf+mjk) + gdil
                                 qq(itf+mjl) = qq(itf+mjl) + gdik
                                 iff = iff + lenbig
                                 itf = itf + lenbig
                                 ifd = ifd + nx
 20                           continue
                           end if
                           ifa = ifa + nx
                           itrf = itrf + nx
                           ioa = ioa + lendd
 30                     continue
 40                  continue
                  end if
 50            continue
 60         continue
 70      continue
 80   continue
      return
      end
      subroutine formeg
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/dmisc)
INCLUDE(common/grad2)
      common/tgrad/dgout(3,3)
      data zero /0.0d0/
      dumx = zero
      dumy = zero
      dumz = zero
      do 20 ipass = 1 , npass
         iat = natomd(ipass)
         dum = dgout(1,ipass)
         dumx = dumx + dum
         de(1,iat) = de(1,iat) + dum
         dum = dgout(2,ipass)
         dumy = dumy + dum
         de(2,iat) = de(2,iat) + dum
         dum = dgout(3,ipass)
         dumz = dumz + dum
         de(3,iat) = de(3,iat) + dum
 20   continue
      iat = natomd(npass+1)
      de(1,iat) = de(1,iat) - dumx
      de(2,iat) = de(2,iat) - dumy
      de(3,iat) = de(3,iat) - dumz
      return
      end
      subroutine grdcon(qq,noform)
      implicit REAL  (a-h,o-z)
      logical nofk,noform
INCLUDE(common/sizes)
_IF(rpagrad)
INCLUDE(common/rpadcom)
INCLUDE(common/infoa)
_ENDIF
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/specal)
      common/tgrad/dgout(9)
      parameter(ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     +  aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/dshlno)
INCLUDE(common/incrd)
INCLUDE(common/dmisc)
_IFN(qmmm)
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
_ENDIF
INCLUDE(common/restrl)
INCLUDE(common/iofile)
      dimension noform(*),qq(*)
      data half/0.5d0/
      nofk = .not.(ofock .or. ofokab)
c
c
      if ((mp2w .or. ompir)) then
         ioa = ic7 + lendd*npass*3
         call vclr(qq(ioa+1),1,3*lendd)
         do 30 i = 1 , npass
            iob = ic7 + (i-1)*lendd*3
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 20 k = 1 , lendd*3
               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
 20         continue
 30      continue
      end if
_IFN(qmmm)
      if (ompir) then
         n2 = iabd + lendd + 1
         do 60 i = 1 , 3
            n1 = ic7 + 1
            do 50 k = 1 , npass + 1
               do 40 j = 1 , 3
                  t1 = ddot(ijkld,qq(n2),1,qq(n1),1)
                  dipi(i,j,natomd(k)) = dipi(i,j,natomd(k)) - t1
                  n1 = n1 + lendd
 40            continue
 50         continue
            n2 = n2 + lendd
 60      continue
      end if
_ELSE
      if (ompir) then
         call caserr('not available in qmmm version')
      endif
_ENDIF
c
      if (outd .or. ofock .or. ofokab) then
c
c     used only if extra output required or if need
c     derivatives of fock matrices
c
         jmax = maxj
         nn = 0
         nnn = 0
         do 110 ii = mini , maxi
            i1 = loci + ii
            if (oiandj) jmax = ii
            do 100 jj = minj , jmax
               j1 = locj + jj
c
c
               lmax = maxl
               do 90 kk = mink , maxk
                  k1 = lock + kk
                  if (okandl) lmax = kk
                  do 80 ll = minl , lmax
                     nn = nn + 1
                     if (.not.(noform(nn))) then
                        nnn = nnn + 1
                        l1 = locl + ll
                        densty = qq(iabd+nnn)
                        n0 = 1
                        n1 = ic7
                        do 70 npp = 1 , npass
                           n2 = n1 + lendd
                           n3 = n2 + lendd
                           dumx = qq(n1+nnn)
                           dumy = qq(n2+nnn)
                           dumz = qq(n3+nnn)
                           dgout(n0) = dgout(n0) + densty*dumx
                           dgout(n0+1) = dgout(n0+1) + densty*dumy
                           dgout(n0+2) = dgout(n0+2) + densty*dumz
                           if (outd) write(iwr,6010) npp , i1 , j1 , 
     +                         k1 , l1 , n1 , n2 , n3 , ic7 , lendd , 
     +                         nn , nnn , dumx , dumy , dumz , densty
                           if (.not.(nofk)) then
                              if (i1.eq.j1) then
                                 dumx = dumx*half
                                 dumy = dumy*half
                                 dumz = dumz*half
                              end if
                              if (k1.eq.l1) then
                                 dumx = dumx*half
                                 dumy = dumy*half
                                 dumz = dumz*half
                              end if
                              qq(n1+nnn) = q4*dumx
                              qq(n2+nnn) = q4*dumy
                              qq(n3+nnn) = q4*dumz
                           end if
                           n0 = n0 + 3
                           n1 = n3 + lendd
 70                     continue
                     end if
 80               continue
 90            continue
 100        continue
 110     continue
      else
         n1 = ic7 + 1
         np3 = npass*3
         do 120 n0 = 1 , np3
            dgout(n0) = dgout(n0) + ddot(ijkld,qq(iabd+1),1,qq(n1),1)
            n1 = n1 + lendd
 120     continue
      end if
_IF(rpagrad)
      if (orpagrad) then
c
c        Compute the AO gradient contributions to the excitation 
c        energy gradient.
c
c
         jmax = maxj
         nn = 0
         nnn = 0
         do ii = mini , maxi
            i1 = loci + ii
            ni = (i1-1)*num
            if (oiandj) jmax = ii
            do jj = minj , jmax
               j1  = locj + jj
               nj  = (j1-1)*num
               lij = ni+j1-1
               lji = nj+i1-1
               lmax = maxl
               do kk = mink , maxk
                  k1  = lock + kk
                  nk  = (k1-1)*num
                  lik = ni+k1-1
                  ljk = nj+k1-1
                  lki = nk+i1-1
                  lkj = nk+j1-1
                  if (okandl) lmax = kk
                  do ll = minl , lmax
                     nn = nn + 1
                     if (.not.(noform(nn))) then
                        nnn = nnn + 1
                        l1  = locl + ll
                        nl  = (l1-1)*num
                        lil = ni+l1-1
                        ljl = nj+l1-1
                        lkl = nk+l1-1
                        lli = nl+i1-1
                        llj = nl+j1-1
                        llk = nl+k1-1
                        n1 = ic7
                        do npp = 1 , npass
                           n2 = n1 + lendd
                           n3 = n2 + lendd
                           dumx = qq(n1+nnn)
                           dumy = qq(n2+nnn)
                           dumz = qq(n3+nnn)
                           do ist = 1, nrpastate
                              iqzsu = ipzsu+(ist-1)*num*num
                              iqysu = ipysu+(ist-1)*num*num
c
                              gijkl = 2.0d0*qq(iqzsu+lij)*qq(iqzsu+lkl)
     +                              + 2.0d0*qq(iqysu+lij)*qq(iqysu+lkl)
     +                              - 2.0d0*qq(iqzsu+lij)*qq(iqysu+lkl)
     +                              - 2.0d0*qq(iqysu+lij)*qq(iqzsu+lkl)
     +                              -       qq(iqzsu+lik)*qq(iqzsu+ljl)
     +                              -       qq(iqysu+lik)*qq(iqysu+ljl)
     +                              +       qq(iqzsu+lil)*qq(iqysu+lkj)
     +                              +       qq(iqysu+lil)*qq(iqzsu+lkj)
c
                              gijlk = 2.0d0*qq(iqzsu+lij)*qq(iqzsu+llk)
     +                              + 2.0d0*qq(iqysu+lij)*qq(iqysu+llk)
     +                              - 2.0d0*qq(iqzsu+lij)*qq(iqysu+llk)
     +                              - 2.0d0*qq(iqysu+lij)*qq(iqzsu+llk)
     +                              -       qq(iqzsu+lil)*qq(iqzsu+ljk)
     +                              -       qq(iqysu+lil)*qq(iqysu+ljk)
     +                              +       qq(iqzsu+lik)*qq(iqysu+llj)
     +                              +       qq(iqysu+lik)*qq(iqzsu+llj)
c
                              gjikl = 2.0d0*qq(iqzsu+lji)*qq(iqzsu+lkl)
     +                              + 2.0d0*qq(iqysu+lji)*qq(iqysu+lkl)
     +                              - 2.0d0*qq(iqzsu+lji)*qq(iqysu+lkl)
     +                              - 2.0d0*qq(iqysu+lji)*qq(iqzsu+lkl)
     +                              -       qq(iqzsu+ljk)*qq(iqzsu+lil)
     +                              -       qq(iqysu+ljk)*qq(iqysu+lil)
     +                              +       qq(iqzsu+ljl)*qq(iqysu+lki)
     +                              +       qq(iqysu+ljl)*qq(iqzsu+lki)
c
                              gjilk = 2.0d0*qq(iqzsu+lji)*qq(iqzsu+llk)
     +                              + 2.0d0*qq(iqysu+lji)*qq(iqysu+llk)
     +                              - 2.0d0*qq(iqzsu+lji)*qq(iqysu+llk)
     +                              - 2.0d0*qq(iqysu+lji)*qq(iqzsu+llk)
     +                              -       qq(iqzsu+ljl)*qq(iqzsu+lik)
     +                              -       qq(iqysu+ljl)*qq(iqysu+lik)
     +                              +       qq(iqzsu+ljk)*qq(iqysu+lli)
     +                              +       qq(iqysu+ljk)*qq(iqzsu+lli)
c
                              gklij = 2.0d0*qq(iqzsu+lkl)*qq(iqzsu+lij)
     +                              + 2.0d0*qq(iqysu+lkl)*qq(iqysu+lij)
     +                              - 2.0d0*qq(iqzsu+lkl)*qq(iqysu+lij)
     +                              - 2.0d0*qq(iqysu+lkl)*qq(iqzsu+lij)
     +                              -       qq(iqzsu+lki)*qq(iqzsu+llj)
     +                              -       qq(iqysu+lki)*qq(iqysu+llj)
     +                              +       qq(iqzsu+lkj)*qq(iqysu+lil)
     +                              +       qq(iqysu+lkj)*qq(iqzsu+lil)
c
                              gklji = 2.0d0*qq(iqzsu+lkl)*qq(iqzsu+lji)
     +                              + 2.0d0*qq(iqysu+lkl)*qq(iqysu+lji)
     +                              - 2.0d0*qq(iqzsu+lkl)*qq(iqysu+lji)
     +                              - 2.0d0*qq(iqysu+lkl)*qq(iqzsu+lji)
     +                              -       qq(iqzsu+lkj)*qq(iqzsu+lli)
     +                              -       qq(iqysu+lkj)*qq(iqysu+lli)
     +                              +       qq(iqzsu+lki)*qq(iqysu+ljl)
     +                              +       qq(iqysu+lki)*qq(iqzsu+ljl)
c
                              glkij = 2.0d0*qq(iqzsu+llk)*qq(iqzsu+lij)
     +                              + 2.0d0*qq(iqysu+llk)*qq(iqysu+lij)
     +                              - 2.0d0*qq(iqzsu+llk)*qq(iqysu+lij)
     +                              - 2.0d0*qq(iqysu+llk)*qq(iqzsu+lij)
     +                              -       qq(iqzsu+lli)*qq(iqzsu+lkj)
     +                              -       qq(iqysu+lli)*qq(iqysu+lkj)
     +                              +       qq(iqzsu+llj)*qq(iqysu+lik)
     +                              +       qq(iqysu+llj)*qq(iqzsu+lik)
c
                              glkji = 2.0d0*qq(iqzsu+llk)*qq(iqzsu+lji)
     +                              + 2.0d0*qq(iqysu+llk)*qq(iqysu+lji)
     +                              - 2.0d0*qq(iqzsu+llk)*qq(iqysu+lji)
     +                              - 2.0d0*qq(iqysu+llk)*qq(iqzsu+lji)
     +                              -       qq(iqzsu+llj)*qq(iqzsu+lki)
     +                              -       qq(iqysu+llj)*qq(iqysu+lki)
     +                              +       qq(iqzsu+lli)*qq(iqysu+ljk)
     +                              +       qq(iqysu+lli)*qq(iqzsu+ljk)
c
                              rpade(1,natomd(npp),ist) =
     +                           rpade(1,natomd(npp),ist) + dumx*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                              rpade(2,natomd(npp),ist) =
     +                           rpade(2,natomd(npp),ist) + dumy*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                              rpade(3,natomd(npp),ist) =
     +                           rpade(3,natomd(npp),ist) + dumz*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                              rpade(1,natomd(npass+1),ist) =
     +                           rpade(1,natomd(npass+1),ist) - dumx*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                              rpade(2,natomd(npass+1),ist) =
     +                           rpade(2,natomd(npass+1),ist) - dumy*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                              rpade(3,natomd(npass+1),ist) =
     +                           rpade(3,natomd(npass+1),ist) - dumz*
     +                          (gijkl+gijlk+gjikl+gjilk+
     +                           gklij+gklji+glkij+glkji)
                           enddo
                           n1 = n3 + lendd
                        enddo
                     end if
                  enddo
               enddo
            enddo
         enddo
      endif
_ENDIF
 6010 format (1x,12i6/30x,3e16.8,5x,e16.8)
      end
      subroutine hsymd(q,iwr)
c
c----------------------------------------------
c   write out derivative fock operators
c---------------------------------------------
c
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
INCLUDE(common/restar)
INCLUDE(common/restrl)
INCLUDE(common/restri)
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
c
INCLUDE(common/cndx41)
INCLUDE(common/specal)
INCLUDE(common/incrd)
INCLUDE(common/dmisc)
INCLUDE(common/infoa)
INCLUDE(common/prnprn)
_IF(secd_parallel)
INCLUDE(common/vcore)
_ENDIF
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     + aej(ncmax),aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
      dimension q(*)
      data half/0.5d0/
c
c     writes out accumulated derivative fock and k matrices
c     onto hamfile.
c
      if (odebug(20)) write (iwr,6010)
      ioff = ifok
      ifocbl = iochf(11)
      if (ogrhf) ifocbl = ifocbl + nat3*lenb

_IF(secd_parallel)
csecd
      call pg_dgop(107, q(ioff+1), nx*nfok, '+')
c     write (6,*)'hsymd dgop done',ipg_nodeid()
_ENDIF
      do 50 n = 1 , nfok
_IF(secd_parallel)
csecd
         call fetch(q(ida+1),nx,'dh',n)
_ELSE
         call rdedx(q(ida+1),nx,ifocbl,ifockf)
_ENDIF
         ij = 0
         do 30 k = 1 , num
            do 20 l = 1 , k
               ij = ij + 1
               if (k.ne.l) then
                  q(ioff+ij) = half*q(ioff+ij)
               end if
 20         continue
 30      continue
         do 40 k = 1 , nx
            q(ida+k) = q(ida+k) + q(ioff+k)
 40      continue

_IF(secd_parallel)
csecd
         call stash(q(ida+1),nx,'dh',n)
_ELSE
         call wrt3(q(ida+1),nx,ifocbl,ifockf)
_ENDIF
         if (odebug(20)) call prtris(q(ida+1),num,iwr)
         ioff = ioff + nx
         ifocbl = ifocbl + lenb
c
c
 50   continue
      lfdtrn = .false.
      if(odebug(30)) write (iwr,6020) iochf(11), iochf(12)
      return
 6010 format (/1x,'**** output from hsymd ****'/)
 6020 format (/1x,'hamfile summary'/
     +         1x,'section 11 at block ',i5/
     +         1x,'section 12 at block ',i5/)
      end
      subroutine jkder(zscftp,q,prefac,iso,nshels)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/prnprn)
      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w
INCLUDE(common/cndx41)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
INCLUDE(common/restar)
INCLUDE(common/restri)
c
INCLUDE(common/specal)
      common/blkin/gijkl(510),mword,nlenx,kworx,kworxx
      common/craypk/ijkl(1360)
      common/blk1/gijkl1(510),mwor1,nlen1,kwor1,kwor11
      common/sortpk/ijkl1(1360)
      common/bufc/gijkl2(510),mwor2,nlen2,kwor2,kwor22
      common/three/ijkl2(1360)
      common/bufd/gijkl3(510),mwor3,nlen3,kwor3,kwor33
      common/lsort/ijkl3(1360)
_IFN(qmmm)
      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
_ENDIF
INCLUDE(common/atmblk)
      common/dbuf/icnt,mxtr,val(5118)
      common/dlabs/iao(5,3120)
INCLUDE(common/ghfblk)
INCLUDE(common/cigrad)
      parameter (ncmax=65)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/dshlno)
INCLUDE(common/symtry)
INCLUDE(common/timez)
      common/tgrad/dgout(9)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),
     + aej(ncmax),aek(ncmax),ael(ncmax),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),
     + dd(4*mxp2),ijden(225),ik(225),
     + ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     + dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/ijlab)
INCLUDE(common/dmisc)
INCLUDE(common/segm)
INCLUDE(common/cslosc)
INCLUDE(common/prints)
INCLUDE(common/incrd)
INCLUDE(common/parallel)
_IF(rpagrad)
INCLUDE(common/rpadcom)
INCLUDE(common/sector)
      common/dims  /nvirtx,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      integer nwords(8)
      REAL rwords
      equivalence(rwords,nwords(1))
      integer irpasec , irpatyp
      integer irpasec1, irpatyp1
      integer irpasec2, irpatyp2
      parameter (irpasec =351,irpatyp =120)
      parameter (irpasec1=352,irpatyp1=121)
      parameter (irpasec2=353,irpatyp2=122)
_ENDIF
c
      dimension dbufs(5120)
      equivalence(dbufs(1),icnt)
c
      dimension iso(nshels,*),q(*)
      dimension m0(48),m1(48),m2(48),m3(48)
      dimension prefac(*)
      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zvb,zmcscf/'casscf','vb','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
c
      indq(m,n) = (m-1)*n
c
      ompir = runtyp.eq.zdipd .or. runtyp.eq.zinfra
      ompir = ompir .and. omp2
c
c     ----- check for grhf or gvb cases -----
c
      ouhf = zscftp.eq.zuhf
      orgvb = zscftp.eq.zgvb
      ogrhf = zscftp.eq.zgrhf
      oclos = zscftp.eq.zrhf
      omcscf = zscftp.eq.zmcscf
c
c     ----- check for casscf
c
      ocas = zscftp.eq.zcas.or.zscftp.eq.zvb
      if (ocas .or. mp3 .or. (omp2 .and. .not.ompir)) then
         call setsto(1360,0,ijkl)
      end if
      if (nprint.ne.-5 .and. oprint(57)) then
         write (iwr,6010)
         call writel(prefac,nshels)
      end if
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
_IF(rpagrad)
      if (orpagrad) then
c
c...     Load Z_su and Y_su for the AO contribution to the excitation
c...     energy gradient.
c
         nstate = nrpastate
c        do i = 1, nirr
c           nstate = nstate + rpa_nstates(i)
c        enddo
         ipzsu = igmem_alloc_inf(nstate*num*num,'drv2e.m','jkder',
     &                           'Z_su',IGMEM_NORMAL)
         ipysu = igmem_alloc_inf(nstate*num*num,'drv2e.m','jkder',
     &                           'Y_su',IGMEM_NORMAL)
         call secget(irpasec1,irpatyp1,irpablk1)
         do i = 1, nstate
            call rdedx(nwords,lenint(8),irpablk1,numdu)
            irpablk1 = irpablk1+1+lensec(nwords(1))+lensec(nwords(2))+
     +          lensec(nwords(3))+lensec(nwords(4))+lensec(nwords(5))+
     +          lensec(nwords(6))
            call rdedx(q(ipzsu+(i-1)*num*num),nwords(7),irpablk1,numdu)
            irpablk1 = irpablk1+lensec(nwords(7))
            call rdedx(q(ipysu+(i-1)*num*num),nwords(8),irpablk1,numdu)
            irpablk1 = irpablk1+lensec(nwords(8))
         enddo
      endif
_ENDIF
c
c     ----- set some parameters -----
c
      nav = lenwrd()
      ntpdm = 1
      oeof = .false.
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
_IFN(qmmm)
      if (ompir) then
         ndenin = 3
         ntpdm = 4
         iflden = idaf
         call secget(isect(31),31,ibden)
         iwor1 = 0
         iwor2 = 0
         iwor3 = 0
         mwor1 = 0
         mwor2 = 0
         mwor3 = 0
         ifmp1 = 20
         ifmp2 = 21
         ifmp3 = 22
         ib1 = 1
         call search(ib1,ifmp1)
         call search(ib1,ifmp2)
         call search(ib1,ifmp3)
         call setsto(1360,0,ijkl1)
         call setsto(1360,0,ijkl2)
         call setsto(1360,0,ijkl3)
         call secget(isect(57),57,iblok)
         call rdedx(dipd,lds(isect(57)),iblok,ifild)
      end if
_ELSE
      if (ompir) then
         call caserr('not available in qmmm version')
      endif
_ENDIF

      if (ogrhf) then
         m = 0
         call secget(isect(53),m,iblok)
         call readi(nact,lds(isect(53))*nav,iblok,idaf)
      end if
      nat3 = nat*3
      nbsq = num*num
      lenb = lensec(nx)
      odpres = .false.
      ofpres = .false.
      ogpres = .false.
      do 20 i = 1 , nshels
         if (ktype(i).eq.3) odpres = .true.
         if (ktype(i).eq.4) ofpres = .true.
         if (ktype(i).eq.5) ogpres = .true.
 20   continue
c
      m = ntpdm + 9
      if (omp2w .or. ompir) m = m + 3
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c code for dgenrl preallocation 
c   NOW all version (previously GA version only)
c
c add buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do 21 i = 1 , nshels
         kt_max=max(kt_max,ktype(i))
 21   continue
      inc1_max=(kt_max+1)**4
c
c vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6
c
c     ----- set pointers for partitioning of core -----
c
c
c revised scheme - 
c

      i10 = igmem_alloc_inf(maxorb,'drv2e.m','jkder','mo2fock',
     &                      IGMEM_DEBUG)
cc      i20 = igmem_alloc(l2)
cc      i30 = igmem_alloc(l2)

      ndens = 1
      if (ouhf) then
         ndens = 2
      end if
      if (orgvb) then
         ndens = 3
      end if
      if (ogrhf) then
         ndens = njk
      end if
      if (ocas) then
         ndens = 1
      end if
      if (ofock .or. ompir) then
         ndens = ndens + ndenin
      end if

      if(o255i) then
       idblk=3120
       iidblk=1950
      else
       idblk=2264
       iidblk=2830
      endif
      idblk1=idblk+1

      ida = igmem_alloc_inf(( ndens)*l2,'drv2e.m','jkder','dens-mats',
     &                      IGMEM_NORMAL)
      i20 = ida
      i30 = i20 + l2
      ida = ida - 1

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
         if (ogrhf) nfok = njk*nat*6
      end if
      if (ofokab) nfok = nat3*ndenin
c

      ifok = igmem_alloc_inf(nx*nfok,'drv2e.m','jkder','der-fock',
     &                       IGMEM_NORMAL)
      ifok = ifok - 1

c  core at i00 is for indexing
      i00 = igmem_alloc_inf(lnddm*(3/nav+1),'drv2e.m','jkder',
     &                      'oform-int',IGMEM_DEBUG)

      iabd = igmem_alloc_inf(lnddm*ntpdm,'drv2e.m','jkder',
     &                       'dens-mat-prod',IGMEM_NORMAL)
      iabd = iabd - 1

      ic7 = igmem_alloc_inf(lnddm*9 + inc1_max*ncmmm_max*6,'drv2e.m',
     &                      'jkder','tmp-int',IGMEM_NORMAL)
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

      if(omcscf)then

         i40 = igmem_alloc(l2)

ccc         id3 = i40 + nx
         numcmo = num*ncoorb
         nactp = ncact*(ncact+1)/2
         nd2str = indq(nactp+1,16)
         if (odpres) nd2str = indq(nactp+1,36)
         if (ofpres) nd2str = indq(nactp+1,100)
         if (ogpres) nd2str = indq(nactp+1,225)

         id3 = igmem_alloc(numcmo)
         i30 = i20
ccc         id4 = id3 + numcmo

c     nd2mo=ind(nactp,nactp)
         nd2mo = nactp*(nactp+1)/2
c     write(6,*)' length of tpdm',nd2mo

         id4 = igmem_alloc(nd2str)
ccc         id5 = id4 + nd2str
         id5 = igmem_alloc(nd2mo)
ccc         id6 = id5 + nd2mo
c
         najkl = indq(ncact*4+1,nactp)
         if (odpres) najkl = indq(ncact*6+1,nactp)
         if (ofpres) najkl = indq(ncact*10+1,nactp)
         if (ogpres) najkl = indq(ncact*15+1,nactp)
c
         nabcl = indq(16+1,4*ncact)
         if (odpres) nabcl = indq(36+1,6*ncact)
         if (ofpres) nabcl = indq(100+1,10*ncact)
         if (ogpres) nabcl = indq(225+1,15*ncact)
c
         id6 = igmem_alloc(nabcl)
cc         id7 = id6 + nabcl
         id7 = igmem_alloc(najkl)
cc         id8 = id7 + najkl
         id8 = igmem_alloc(nactp)
cc         id9 = id8 + nactp
c
         nshdm = max(ncact,4)
         if (odpres) nshdm = max(ncact,6)
         if (ofpres) nshdm = max(ncact,10)
         if (ogpres) nshdm = max(ncact,15)
c
         id9 = igmem_alloc(nshdm*nshdm)
cc
cc         iabd = id9 + nshdm*nshdm
c
c   iabd should not be used for mcscf calculations - tpdm in id5
c
cc         i00 = iabd + lnddm*ntpdm
cc         ic7 = i00 + lnddm*(3/nav+1) - 1
      endif


c     tim0 = cpulft(1)
c
      if (ocas) then
         call dbutci(ist,jst,kst,lst)
      else if (omcscf) then
         call ddebut(zscftp,q(id3),q(i30),q(i10),q(i40),q(id5),l1,
     +               l2,l3,ist,jst,kst,lst,q)
      else
         call ddebut(zscftp,q(i20),q(i30),q(i10),q(i10),q(ifok+1),
     +               l1,l2,l3,ist,jst,kst,lst,q)
      end if
      if (mp3 .or. (omp2 .and. .not.ompir)) call search(iblk2d,ifil2d)
      mword = 0
      iword = 0
      kworx = 999
      kwor1 = 999
      kwor2 = 999
      kwor3 = 999
      kloc(nshels+1) = num + 1
c
      icnt = 0
      ib1 = 1
      if (omp2w) then
         call search(ib1,mpstrm(1))
         if (odebug(30)) write (iwr,6020) mpstrm(1)
      end if
_IF(parallel)
      next = ipg_dlbtask()
_ENDIF
      if (ist.le.nshels) then
c
c     ----- ishell -----
c
_IF(parallel)
c         do 140 ii = nshels, ist, -1
         do 140 ii = ist , nshels
_ELSE
         do 140 ii = ist , nshels
_ENDIF
            kadi = kad(ii)
            ijshel = ii*(ii-1)/2
            do 40 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 140
               m0(it) = id
 40         continue
            iceni = katom(ii)
            if (omcscf) call mcajkl(q(id3),q(id5),q(id7),q(id8),q(id9),
     +                              ii,nactp)
c
c     ----- jshell -----
c
            j0 = jst
_IF(parallel)
c            do 130 jj = ii, j0, -1
            do 130 jj = j0 , ii
_ELSE
            do 130 jj = j0 , ii
_ENDIF
               kadij = kadi + kad(jj)
               jst = 1
               tolij = dlntol + prefac(ijshel+jj)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
               if (tolij.gt.-3.401d0) then
                  do 60 it = 1 , nt
                     id = m0(it)
                     jd = iso(jj,it)
                     if (jd.gt.ii) go to 130
                     if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (id.eq.ii .and. jd.gt.jj) go to 130
                     m1(it) = id
                     m2(it) = jd
 60               continue
                  olab = katom(jj).eq.iceni
c
c     store information about the pair (ij)
c
                  call dshell(1,ii,jj,kk,ll)
                  call dprim
                  if (nij.ne.0) then
                     if (omcscf) call mcabkl(q(id3),q(id7),q(id4),ii,jj,
     +                   nactp)
_IF(parallel)
c***   **MPP**
                     icount_dlb = icount_dlb + 1
                     if(icount_dlb . eq. next) then
c***   **MPP**
_ENDIF
c
c     ----- kshell -----
c
                     k0 = kst
                     do 120 kk = k0 , ii
                        kadijk = kadij + kad(kk)
                        kst = 1
                        klshel = kk*(kk-1)/2
                        do 80 it = 1 , nt
                           kd = iso(kk,it)
                           if (kd.gt.ii) go to 120
                           m3(it) = kd
 80                     continue
                        olabc = olab .and. katom(kk).eq.iceni
                        if (omcscf)
     +                      call mcabcl(q(id3),q(id4),q(id6),q(id8),
     +                      q(id9),ii,jj,kk,nactp)
                        l0 = lst
                        maxll = kk
                        if (kk.eq.ii) maxll = jj
                        do 110 ll = l0 , maxll
                           lst = 1
                           if (kadijk+kad(ll).lt.0) then
                              tolijk = tolij + prefac(klshel+ll)
                              if (tolijk.gt.0.0d0) then
                                 olabcd = olabc .and. katom(ll).eq.iceni
                                 if (.not.(olabcd)) then
                                    n4 = 0
                                    do 100 it = 1 , nt
                                       ld = iso(ll,it)
                                       if (ld.gt.ii) go to 110
                                       kd = m3(it)
                                       if (kd.lt.ld) then
                                         nd = kd
                                         kd = ld
                                         ld = nd
                                       end if
                                       id = m1(it)
                                       jd = m2(it)
                                       if (id.eq.ii .or. kd.eq.ii) then
                                         if (kd.ge.id) then
                                         if (kd.ne.id .or. ld.gt.jd)
     +                                      then
                                         nd = id
                                         id = kd
                                         kd = nd
                                         nd = jd
                                         jd = ld
                                         ld = nd
                                         end if
                                         end if
                                         if (jd.ge.jj) then
                                         if (jd.gt.jj) go to 110
                                         if (kd.ge.kk) then
                                         if (kd.gt.kk) go to 110
                                         if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 110
                                         n4 = n4 + 1
                                         end if
                                         end if
                                         end if
                                       end if
 100                                continue
c
c     ----- calculate q4 factor for this group of shells -----
c
                              q4 = dfloat(nt)/dfloat(n4)
c
c     ----- check for redundant combinations -----
c
                              call redund(ii,jj,kk,ll,iwr)
                              if (npass.ne.0) then
c
c     ----- initialize dgout to zero -----
c
                          call vclr(dgout,1,ndgout)
                                 call dshell(2,ii,jj,kk,ll)
c
c     ----- form products of density matrix elements -----
c
                          call vclr(q(iabd+1),1,lendd*ntpdm)
                                 if (omcscf)
     +                              call mcabcd(q(id3),q(id6),
     +                              q(iabd+1),q(id9),ii,jj,kk,ll,
     +                              q4)
                                 if (mp3 .or.
     +                              (omp2 .and. .not.ompir))
     +                              call mcdab(q(iabd+1),ii,jj,kk,
     +                              ll,q4)
                                 if (ompir) then
                                   iabd1 = iabd + lendd + 1
                                   iabd2 = iabd1 + lendd
                                   iabd3 = iabd2 + lendd
                                   call dpdab1(q(iabd1),ii,jj,kk,
     +                                ll,q4,ifmp1,iwor1)
                                   call dpdab2(q(iabd2),ii,jj,kk,
     +                                ll,q4,ifmp2,iwor2)
                                   call dpdab3(q(iabd3),ii,jj,kk,
     +                                ll,q4,ifmp3,iwor3)
                                   call tpdder(q(iabd1),lendd,
     +                                q(i20),q(i30),nx,ndenin,ii,
     +                                jj,kk,ll,q4)
                                 end if
                                 if (.not.orgvb .and.
     +                              .not.ocas .and. .not.omcscf)
     +                              call dabab(ii,jj,kk,ll,q4,
     +                              zscftp,q(i20),q(i30),q(iabd+1)
     +                              )
                                 if (orgvb)
     +                              call dabg(ii,jj,kk,ll,l1,norb,
     +                              q4,q(i20),q(i30),q(i10)
     +                              ,onocor,onopen,q(iabd+1))
                                 if (ocas)
     +                              call dabci(ii,jj,kk,ll,q4,
     +                              oeof,q(iabd+1))
                                 if (omcscf)
     +                              call dabmc(ii,jj,kk,ll,q4, 
     +                              q(i30),q(i40),q(iabd+1))
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
                                 call delim(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)
                                 if (ompir) then
                                   call delim2(q(iabd1),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                   call delim2(q(iabd2),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                   call delim2(q(iabd3),q(ic7+1),
     +                                ijgt,klgt,q(i00),ijd,kld,
     +                                ijkld,lendd)
                                 end if
                                 if (ijkld.ne.0) then
                                   call dgenrl(q(1),q(i00),q(i00),
     +                                abmax)
                                   if (ofock)
     +                                call fockd2(q,q(i00))
                                   if (ofokab) then
                                   call fokabd(q,q(i00))
                                   end if
c
c     ----- generate all 4 partial contributions to the gradient ----
c
                                   call formeg
                                 end if
                              end if
                                 end if
                              end if
                           end if
 110                    continue
 120                 continue
_IF(parallel)
c***   **MPP**
                     next = ipg_dlbtask()
                     endif
c***   **MPP**
_ENDIF
                  end if
               end if
 130        continue
c
c     ----- save gradient and restart data -----
c
            call dfinal(q,0,ii)
            if (tim.ge.timlim) go to 150
 140     continue
_IF(parallel)
            call pg_dlbpush
_ENDIF
      end if
c
c     ----- end of *shell* loops -----
c
      if (omp2w) then
         mderb = idblk + iidblk + 2
         if (icnt.ne.0) then
            call pack(val(idblk1),lab816,iao,idblk*5)
            call wrt3s(dbufs,mderb,mpstrm(1))
         end if
         icnt = 0
         call wrt3s(dbufs,mderb,mpstrm(1))
         call shut1(mpstrm(1))
         if (odebug(30)) write (iwr,6030) mpstrm(1)
      end if
      if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
       isymd = igmem_alloc_inf(nw196(6),'drv2e.m','jkder',
     &                         'trans-atoms-symm',IGMEM_DEBUG)
c
       call symde(q(isymd),nat)
c
c     ----- reset core memory from symde
c
      call gmem_free_inf(isymd,'drv2e.m','jkder','trans-atoms-symm')
c
      endif
      call dfinal(q,1,ii)
      nindmx = 0
 150  call timit(0)
c     dtim = tim - tim0
c
c     ----- reset core memory from jkder -----
c

      if(omcscf)then
         call gmem_free(id9)
         call gmem_free(id8)
         call gmem_free(id7)
         call gmem_free(id6)
         call gmem_free(id5)
         call gmem_free(id4)
         call gmem_free(id3)
         call gmem_free(i40)


      endif
cc
cc revert
cc      call gmem_free(ic1)
      ic7 = ic7 + 1
      call gmem_free_inf(ic7,'drv2e.m','jkder','tmp-int')
      iabd = iabd + 1
      call gmem_free_inf(iabd,'drv2e.m','jkder','dens-mat-prod')
      call gmem_free_inf(i00,'drv2e.m','jkder','oform-int')
      ifok = ifok + 1
      call gmem_free_inf(ifok,'drv2e.m','jkder','der-fock')
      ida = ida + 1
      call gmem_free_inf(ida,'drv2e.m','jkder','dens-mats')
      call gmem_free_inf(i10,'drv2e.m','jkder','mo2fock')
_IF(rpagrad)
      if (orpagrad) then
         call gmem_free_inf(ipysu,'drv2e.m','jkder','Y_su')
         call gmem_free_inf(ipzsu,'drv2e.m','jkder','Z_su')
      endif
_ENDIF
      return
 6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
 6020 format (/1x,'derivative integrals to be output to unit',i3)
 6030 format (/1x,'derivative integrals written to unit',i3)
      end
      subroutine mcabcd(cc,gabcl,gabcd,gsq,ish,jsh,ksh,lsh,q4)
      implicit REAL  (a-h,o-z)
      dimension cc(*),gabcl(*),gabcd(*),gsq(*)
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/cigrad)
INCLUDE(common/incrd)
      data pt5/0.5d0/
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      minc = kloc(ksh)
      mind = kloc(lsh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      maxc = kloc(ksh+1) - 1
      maxd = kloc(lsh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nc = maxc - minc + 1
      nd = maxd - mind + 1
c     nab = na*nb
      ncd = nc*nd
      ncl = nc*ncact
      lmin = ncore*num + mind
      qq4 = 8.0d0*q4
c
      iab = 0
      do 50 ib = minb , maxb
         ib1 = (ib-minb)*inc4
         do 40 ia = mina , maxa
            iaa = (ia-mina)*inc5
            ibb = ib1 + iaa
            iab = iab + 1
c
            dd4 = qq4
            if (ia.eq.ib) dd4 = qq4*pt5
            iabcl = indq(iab,ncl) + 1
            call vclr(gsq(1),1,ncd)
            call mxmb(gabcl(iabcl),1,nc,cc(lmin),num,1,gsq,1,nc,nc,
     +                ncact,nd)
c
            icd = 0
            do 30 id = mind , maxd
               id1 = id - mind + inc2
               do 20 ic = minc , maxc
                  icc = (ic-minc)*inc3 + ibb
                  idd = id1 + icc
c
                  ddd4 = dd4
                  if (ic.eq.id) ddd4 = ddd4*pt5
                  icd = icd + 1
                  gabcd(idd) = gsq(icd)*ddd4
 20            continue
 30         continue
 40      continue
 50   continue
      return
      end
      subroutine mcabcl(cc,gabkl,gabcl,gtr,gsq,ish,jsh,ksh,nkl)
      implicit REAL  (a-h,o-z)
      dimension cc(*),gabkl(*),gabcl(*),gtr(*),gsq(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      minc = kloc(ksh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      maxc = kloc(ksh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nc = maxc - minc + 1
      nab = na*nb
      ncl = nc*ncact
      kkmin = ncore*num + minc
c
      iab = 0
      call vclr(gabcl(1),1,nab*ncl)
      do 40 ib = 1 , nb
         do 30 ia = 1 , na
            iab = iab + 1
            icol2 = indq(iab,ncl) + 1
            ikl = iab
c
            do 20 kl = 1 , nkl
               gtr(kl) = gabkl(ikl)
               ikl = ikl + nab
 20         continue
            call squr(gtr(1),gsq(1),ncact)
            call mxmb(cc(kkmin),1,num,gsq,1,ncact,gabcl(icol2),1,nc,nc,
     +                ncact,ncact)
 30      continue
 40   continue
      return
      end
      subroutine mcabkl(cc,gajkl,gabkl,ish,jsh,nkl)
      implicit REAL  (a-h,o-z)
      dimension cc(*),gajkl(*),gabkl(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
c
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      minb = kloc(jsh)
      maxa = kloc(ish+1) - 1
      maxb = kloc(jsh+1) - 1
      na = maxa - mina + 1
      nb = maxb - minb + 1
      nab = na*nb
      naj = na*ncact
      jmin = ncore*num + minb
c
      call vclr(gabkl(1),1,nab*nkl)
      do 20 kl = 1 , nkl
         icol1 = indq(kl,naj) + 1
         icol2 = indq(kl,nab) + 1
         call mxmb(gajkl(icol1),1,na,cc(jmin),num,1,gabkl(icol2),1,na,
     +             na,ncact,nb)
 20   continue
      return
      end
      subroutine mcajkl(cc,gijkl,gajkl,gtr,gsq,ish,nkl)
      implicit REAL  (a-h,o-z)
      dimension cc(*),gijkl(*),gajkl(*),gtr(*),gsq(*)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/cigrad)
INCLUDE(common/mapper)
c
      ind(i,j) = iky(max(i,j)) + min(i,j)
      indq(m,n) = (m-1)*n
c
      mina = kloc(ish)
      maxa = kloc(ish+1) - 1
      na = maxa - mina + 1
      naj = ncact*na
      imin = ncore*num + mina
c
      call vclr(gajkl(1),1,naj*nkl)
      do 30 kl = 1 , nkl
         icol1 = ind(kl,1)
         icol2 = indq(kl,naj) + 1
         call dcopy(kl,gijkl(icol1),1,gtr(1),1)
c
         icol1 = kl
         do 20 ij = kl + 1 , nkl
            ijkl = ind(ij,kl)
            gtr(icol1+1) = gijkl(ijkl)
            icol1 = icol1 + 1
 20      continue
c
         call squr(gtr(1),gsq(1),ncact)
         call mxmb(cc(imin),1,num,gsq,1,ncact,gajkl(icol2),1,na,na,
     +        ncact,ncact)
 30   continue
      return
      end
      subroutine mcdab(abdens,ii,jj,kk,ll,q4)
c================================================
c    read two particle density matrix from file
c================================================
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
      logical ijdiff,kldiff,ikdiff
INCLUDE(common/cigrad)
INCLUDE(common/incrd)
INCLUDE(common/atmblk)
_IFN1(iv)      common/craypk/labs(1360)
_IF1(iv)      common/craypk/i205(340),j205(340),k205(340),l205(340)
      common/blkin/gin(510),mword,nlenx,kworx,kworxx
c
      dimension abdens(*)
      call vclr(abdens,1,lendd)
      if (kworx.eq.0) return
      ijdiff = ii.ne.jj
      kldiff = kk.ne.ll
      ikdiff = ii.ne.kk .or. jj.ne.ll
      mini = kloc(ii)
      minj = kloc(jj)
      mink = kloc(kk)
      minl = kloc(ll)
      maxi = kloc(ii+1) - 1
      maxj = kloc(jj+1) - 1
      maxk = kloc(kk+1) - 1
      maxl = kloc(ll+1) - 1
 20   if (iword.eq.mword) then
         call find(ifil2d)
         call get(gin,kworx)
         iword = 0
_IFN1(iv)         call unpack(gin(num2e+1),lab816,labs,numlab)
_IF1(iv)         call upak8v(gin(num2e+1),i205)
         if (kworx.le.0) go to 30
      end if
      iword = iword + 1
_IFN1(iv)      m = (iword+iword) + (iword+iword)
_IFN1(iv)      i = labs(m-3)
_IFN1(iv)      j = labs(m-2)
_IFN1(iv)      k = labs(m-1)
_IFN1(iv)      l = labs(m)
_IF1(iv)      i = i205(iword)
_IF1(iv)      j = j205(iword)
_IF1(iv)      k = k205(iword)
_IF1(iv)      l = l205(iword)
      if (i.lt.mini) go to 20
      if (i.le.maxi) then
         if (j.lt.minj) go to 20
         if (j.le.maxj) then
            if (k.lt.mink) go to 20
            if (k.le.maxk) then
               if (l.lt.minl) go to 20
               if (l.le.maxl) then
                  d4 = 8.0d0
                  if (i.eq.j) d4 = 4.0d0
                  if (k.eq.l) d4 = d4*0.5d0
                  nn = (i-mini)*inc5 + (j-minj)*inc4 + (k-mink)
     +                 *inc3 + l - minl + inc2
                  buff = gin(iword)*q4*d4
                  abdens(nn) = buff
                  if (.not.(ijdiff)) then
                     nn = (j-minj)*inc5 + (i-mini)*inc4 + (k-mink)
     +                    *inc3 + l - minl + inc2
                     abdens(nn) = buff
                     if (.not.(kldiff)) then
                        nn = (j-minj)*inc5 + (i-mini)*inc4 + (l-minl)
     +                       *inc3 + k - mink + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  if (.not.(kldiff)) then
                     nn = (i-mini)*inc5 + (j-minj)*inc4 + (l-minl)
     +                    *inc3 + k - mink + inc2
                     abdens(nn) = buff
                  end if
                  if (.not.(ikdiff)) then
                     nn = (k-mink)*inc5 + (l-minl)*inc4 + (i-mini)
     +                    *inc3 + j - minj + inc2
                     abdens(nn) = buff
                     if (.not.(ijdiff)) then
                        nn = (k-mink)*inc5 + (l-minl)*inc4 + (j-minj)
     +                       *inc3 + i - mini + inc2
                        abdens(nn) = buff
                        if (.not.(kldiff)) then
                           nn = (l-minl)*inc5 + (k-mink)*inc4 + (j-minj)
     +                          *inc3 + i - mini + inc2
                           abdens(nn) = buff
                        end if
                     end if
                     if (.not.(kldiff)) then
                        nn = (l-minl)*inc5 + (k-mink)*inc4 + (i-mini)
     +                       *inc3 + j - minj + inc2
                        abdens(nn) = buff
                     end if
                  end if
                  go to 20
               end if
            end if
         end if
      end if
 30   iword = iword - 1
      return
      end
_IF(rpagrad)
      subroutine redund(ii,jj,kk,ll,iw)
      implicit REAL  (a-h,o-z)
      logical ieqj,ieqk,ieql,jeqk,jeql,keql
c
c...  Huub van Dam: 1998-03-19
c
c...  Reimplementation of original subroutine redund. This version 
c...  should be correct and it is hoped that it is easier to read too.
c
c
c...  This routine analyses the 4 functions that define the 
c...  2-electron integral to eliminate redundant derivatives.
c
c...  In case a redundancy is found some coordinates must be chosen 
c...  for which the derivatives with respect to these coordinates will
c...  be skipped. Whenever possible the coordinates should be selected
c...  such that the calculation of the derivative integrals costs the 
c...  least. 
c
c...  In this respect one should remember that the derivative
c...  of a Cartesian Gaussian is a linear combination of Cartesian
c...  Gaussians with polynomial factors of a higher order.
c...  E.g.: (d/dx) x**4 exp(-a*x**2) = (4*x**3 - 2*a*x**5) exp(-a*x**2)
c...  Furthermore, the number of functions in a shell (nf) scales with 
c...  the maximal power [or angular momentum (l)] as
c
c...     nf = (l+1)*(l+2)/2
c
c...  so for given li, lj, lk, ll the number of 2-electron integrals
c...  (nint) will be
c
c...     nint = nfi*nfj*nfk*nfl
c
c...  So if one can choose a center to be skipped in taking the 
c...  derivatives one should the centers such that nint remains 
c...  as small as possible.
c
c...  In this respect one should note that for a given value of the
c...  sum 
c
c...     s = nfi+nfj+nfk+nfl 
c
c...  nint will be the smallest if all nf's except one are equal 1,
c...  whereas nint will be largest if the nf's have approximately
c...  the same value.
c
c...  The redundancies depend on the number of shells being centered
c...  at the same atom. Therefore, we will check which shells are on the
c...  same center. The 4 shell labels result in 6 pairs of centers. 
c...  The centers in a pair may be equal or not equal, but these 
c...  relations are not independent. In the following table all possible
c...  combinations of relations are given. If the centers in the pair
c...  ij are equal this is denoted by 'T' or 't', otherwise there is a
c...  'F' or 'f'. On each line of the table we start a the left chosing
c...  values for relations until our choices + the interdepencies
c...  uniquely define all values. The values we chose are marked 'T' or
c...  'F' whereas the values that follow from the interdepencies are
c...  marked 't' and 'f'.
c
c...     number ij ik il jk jl kl
c...        1   T  T  T  t  t  t
c...        2   T  T  F  t  f  f
c...        3   T  F  T  f  t  f
c...        4   T  F  F  f  f  T
c...        5   T  F  F  f  f  F
c...        6   F  T  T  f  f  t
c...        7   F  T  F  f  T  f
c...        8   F  T  F  f  F  f
c...        9   F  F  T  T  f  f
c...       10   F  F  T  F  f  f
c...       11   F  F  F  T  T  t
c...       12   F  F  F  T  F  f
c...       13   F  F  F  F  T  f
c...       14   F  F  F  F  F  T
c...       15   F  F  F  F  F  F
c
c...  In practice the calling routine jkder enforces the following
c...  relations on the function labels:
c
c...       jj <= ii
c...       kk <= ii
c...       ll <= kk, if kk < ii
c...       ll <= jj, if kk = ii & jj < ii
c...       ll <  jj, otherwise
c
INCLUDE(common/sizes)
INCLUDE(common/dmisc)
INCLUDE(common/nshel)
c
c...  <lll> is used to find the function with lowest angular momentum 
c...  in case all functions are on different centers.
c
      dimension lll(4)
      equivalence (lll(1),lit),(lll(2),ljt),(lll(3),lkt),(lll(4),llt)
c
c...  <lla> is used for easy permutation of the center with lowest 
c...  angular momentum to the last position in case all functions are 
c...  on different centers. This facilitates setting <natomd>
c
      dimension lla(4)
      equivalence (lla(1),iat),(lla(2),jat),(lla(3),kat),(lla(4),lat)
c
c...  See code at label 150 for details on usage of <lll> and <lla>.
c
      oskip(1) = .true.
      oskip(2) = .true.
      oskip(3) = .true.
      oskip(4) = .true.
      npass = 0
      do i = 1 , 4
         natomd(i) = 0
      enddo
      lit = ktype(ii)
      ljt = ktype(jj)
      lkt = ktype(kk)
      llt = ktype(ll)
      iat = katom(ii)
      jat = katom(jj)
      kat = katom(kk)
      lat = katom(ll)
      ieqj = iat.eq.jat
      ieqk = iat.eq.kat
      ieql = iat.eq.lat
      jeqk = jat.eq.kat
      jeql = jat.eq.lat
      keql = kat.eq.lat
      if (ieqj) then
         if (ieqk) then
            if (ieql) then
c
c...           ij=T & ik=T & il=T                               => no  1
c
               go to 10
            else
c
c...           ij=T & ik=T & il=F                               => no  2
c
               go to 20
            endif
         else
            if (ieql) then
c
c...           ij=T & ik=F & il=T                               => no  3
c
               go to 30
            else
               if (keql) then
c
c...              ij=T & ik=F & il=F & kl=T                     => no  4
c
                  go to 40
               else
c
c...              ij=T & ik=F & il=F & kl=F                     => no  5
c
                  go to 50
               endif
            endif
         endif
      else
         if (ieqk) then
            if (ieql) then
c
c...           ij=F & ik=T & il=T                               => no  6
c
               go to 60
            else
               if (jeql) then
c
c...              ij=F & ik=T & il=F & jl=T                     => no  7
c
                  go to 70
               else
c
c...              ij=F & ik=T & il=F & jl=F                     => no  8
c
                  go to 80
               endif
            endif
         else
            if (ieql) then
               if (jeqk) then
c
c...              ij=F & ik=F & il=T & jk=T                     => no  9
c
                  go to 90
               else
c
c...              ij=F & ik=F & il=T & jk=F                     => no 10
c
                  go to 100
               endif
            else
               if (jeqk) then
                  if (jeql) then
c
c...                 ij=F & ik=F & il=F & jk=T & jl=T           => no 11
c
                     go to 110
                  else
c
c...                 ij=F & ik=F & il=F & jk=T & jl=F           => no 12
c
                     go to 120
                  endif
               else
                  if (jeql) then
c
c...                 ij=F & ik=F & il=F & jk=F & jl=T           => no 13
c
                     go to 130
                  else
                     if (keql) then
c
c...                    ij=F & ik=F & il=F & jk=F & jl=F & kl=T => no 14
c
                        go to 140
                     else
c
c...                    ij=F & ik=F & il=F & jk=F & jl=F & kl=F => no 15
c
                        go to 150
                     endif
                  endif
               endif
            endif
         endif
      endif
      call caserr("redund: should never get here!!!")
c
c...  Table entry number 1 
c...  iat = jat = kat = lat
c
 10   continue
c...     All derivatives are redundant so:
c...     skip all, perform no passes.
c
c...     Actually we should never get here because routine jkder
c...     assures that the current routine is never called if 
c...     iat = jat = kat = lat.
         call caserr(" *** REDUND: iat=jat=kat=lat JKDER failed !!!")
      go to 900
c
c...  Table entry number 2 
c...  iat = jat = kat # lat
c
 20   continue
         oskip(4)  = .false.
         natomd(1) = lat
         natomd(2) = iat
         npass     = 1
      go to 900
c
c...  Table entry number 3 
c...  iat = jat = lat # kat
c
 30   continue
         oskip(3)  = .false.
         natomd(1) = kat
         natomd(2) = iat
         npass     = 1
      go to 900
c
c...  Table entry number 4
c...  iat = jat # kat = lat
c
 40   continue
c...     Differentiate with respect to either iat or kat.
c...     Select the most efficient option (n1 and n2 yield an indication
c...     for the number of integrals resulting from differentiating
c...     with respect iat or kat respectively).
         n1 = (lit+1)*(ljt+1)*lkt*llt
         n2 = lit*ljt*(lkt+1)*(llt+1)
         if (n1.lt.n2) then
c...        Differentiate with respect to iat
            oskip(1)  = .false.
            oskip(2)  = .false.
            natomd(1) = iat
            natomd(2) = jat
            natomd(3) = kat
            npass     = 2
         else
c...        Differentiate with respect to kat
            oskip(3)  = .false.
            oskip(4)  = .false.
            natomd(1) = kat
            natomd(2) = lat
            natomd(3) = iat
            npass     = 2
         endif
      go to 900
c
c...  Table entry number  5
c...  iat = jat # kat # lat # iat
c
 50   continue
c...     Differentiate with respect to kat and lat
         oskip(3)  = .false.
         oskip(4)  = .false.
         natomd(1) = kat
         natomd(2) = lat
         natomd(3) = iat
         npass     = 2
      go to 900
c
c...  Table entry number  6
c...  iat = kat = lat # jat
c
 60   continue
c...     Differentiate with respect to jat
         oskip(2)  = .false.
         natomd(1) = jat
         natomd(2) = iat
         npass     = 1
      go to 900
c
c...  Table entry number  7
c...  iat = kat # jat = lat
c
 70   continue
         if (ii.eq.kk.and.jj.eq.ll) then
            if (lit.ge.ljt) then
               oskip(1)  = .false.
               natomd(1) = iat
               natomd(2) = jat
               npass     = 1
            else
               oskip(2)  = .false.
               natomd(1) = jat
               natomd(2) = iat
               npass     = 1
            endif
         else 
c...        Differentiate with respect to either iat or jat.
c...        Select the most efficient option (n1 and n2 yield an 
c...        indication for the number of integrals resulting from 
c...        differentiating with respect iat or jat respectively).
            n1 = (lit+1)*ljt*(lkt+1)*llt
            n2 = lit*(ljt+1)*lkt*(llt+1)
            if (n1.lt.n2) then
c...           Differentiate with respect to iat
               oskip(1)  = .false.
               oskip(3)  = .false.
               natomd(1) = iat
               natomd(2) = kat
               natomd(3) = jat
               npass     = 2
            else
c...           Differentiate with respect to jat
               oskip(2)  = .false.
               oskip(4)  = .false.
               natomd(1) = jat
               natomd(2) = lat
               natomd(3) = iat
               npass     = 2
            endif
         endif
      go to 900
c
c...  Table entry number  8
c...  iat = kat # jat # lat # iat
c
 80   continue
c...     Differentiate with respect to jat and lat
         oskip(2)  = .false.
         oskip(4)  = .false.
         natomd(1) = jat
         natomd(2) = lat
         natomd(3) = iat
         npass     = 2
      go to 900
c
c...  Table entry number 9
c...  iat = lat # jat = kat
c
 90   continue
c
c...     Note: The treatment implemented here differs from the treatment
c...     for Table entry number 7. The reason is that the restrictions
c...     on the function labels imposed by jkder insist that if
c...     ii = ll is to hold then ii = kk must hold too. However, this
c...     implies that iat = kat = lat which contradicts with the path
c...     that led here. So the check (ii.eq.ll.and.jj.eq.kk) can be 
c...     omitted here.
c
c...     Differentiate with respect to either iat or jat.
c...     Select the most efficient option (n1 and n2 yield an indication
c...     for the number of integrals resulting from differentiating
c...     with respect iat or jat respectively).
         n1 = (lit+1)*ljt*lkt*(llt+1)
         n2 = lit*(ljt+1)*(lkt+1)*llt
         if (n1.lt.n2) then
c...        Differentiate with respect to iat
            oskip(1)  = .false.
            oskip(4)  = .false.
            natomd(1) = iat
            natomd(2) = lat
            natomd(3) = jat
            npass     = 2
         else
c...        Differentiate with respect to jat
            oskip(2)  = .false.
            oskip(3)  = .false.
            natomd(1) = jat
            natomd(2) = kat
            natomd(3) = iat
            npass     = 2
         endif
      go to 900
c
c...  Table entry number 10
c...  iat = lat # jat # kat # iat
c
100   continue
c...     Differentiate with respect to jat and kat
         oskip(2)  = .false.
         oskip(3)  = .false.
         natomd(1) = jat
         natomd(2) = kat
         natomd(3) = iat
         npass     = 2
      go to 900
c
c...  Table entry number 11
c...  jat = kat = lat # iat
c
110   continue
c...     Differentiate with respect to iat
         oskip(1)  = .false.
         natomd(1) = iat
         natomd(2) = jat
         npass     = 1
      go to 900
c
c...  Table entry number 12
c...  jat = kat # iat # lat # jat
c
120   continue
c...     Differentiate with respect to iat and lat
         oskip(1)  = .false.
         oskip(4)  = .false.
         natomd(1) = iat
         natomd(2) = lat
         natomd(3) = jat
         npass     = 2
      go to 900
c
c...  Table entry number 13
c...  jat = lat # iat # kat # jat
c
130   continue
c...     Differentiate with respect to iat and kat
         oskip(1)  = .false.
         oskip(3)  = .false.
         natomd(1) = iat
         natomd(2) = kat
         natomd(3) = jat
         npass     = 2
      go to 900
c
c...  Table entry number 14
c...  kat = lat # iat # jat # kat
c
140   continue
c...     Differentiate with respect to iat and jat
         oskip(1)  = .false.
         oskip(2)  = .false.
         natomd(1) = iat
         natomd(2) = jat
         natomd(3) = kat
         npass     = 2
      go to 900
c
c...  Table entry number 15
c...  All centers different
c
150   continue
c...     Differentiate with respect to all centers except one.
         min = lll(1)
         imn = 1
         do i = 2, 4
            if (lll(i).lt.min) then
               min = lll(i)
               imn = i
            endif
         enddo
         mnat     = lla(imn)
         do i = imn+1, 4
           lla(i-1) = lla(i)
         enddo
         lla(4)   = mnat
         do i = 1, 4
           oskip(i)  = .false.
           natomd(i) = lla(i)
         enddo
         oskip(imn) = .true.
         npass      = 3
      go to 900
c
c...  Leave routine
c
900   continue
      if (outd) then
         write (iw,6010)ii, jj, kk, ll, oskip, npass, natomd
         write (*,*)' iat,jat,kat,lat = ',iat,jat,kat,lat
         write (*,*)' lit,ljt,lkt,llt = ',lit,ljt,lkt,llt
      endif
      return
 6010 format (/,' *** REDUND: ii,jj,kk,ll =',4i3,' skip =',4l3,
     +        ' npass =',i2,' centers =',4i5,/)
      end
_ELSE
      subroutine redund(ii,jj,kk,ll,iw)
      implicit REAL  (a-h,o-z)
      logical inej,inek,inel,jnek,jnel,knel
INCLUDE(common/sizes)
INCLUDE(common/dmisc)
INCLUDE(common/nshel)
c
      dimension lll(4)
      equivalence (lll(1),lit)
      oskip(1) = .true.
      oskip(2) = .true.
      oskip(3) = .true.
      oskip(4) = .true.
      npass = 0
      do 20 i = 1 , 4
         natomd(i) = 0
 20   continue
      lit = ktype(ii)
      ljt = ktype(jj)
      lkt = ktype(kk)
      llt = ktype(ll)
      iat = katom(ii)
      jat = katom(jj)
      kat = katom(kk)
      lat = katom(ll)
      inej = iat.ne.jat
      inek = iat.ne.kat
      inel = iat.ne.lat
      jnek = jat.ne.kat
      jnel = jat.ne.lat
      knel = kat.ne.lat
      if (inej) then
         if (.not.(inek)) then
            if (.not.(inel)) go to 40
c      iat=kat    jat=lat
            if (jnel) go to 50
            if (ii.ne.kk .or. jj.ne.ll) then
               n1 = (lit+1)*(lkt+1)*ljt*llt
               n2 = lit*lkt*(ljt+1)*(llt+1)
               if (n1.ge.n2) go to 50
               go to 70
            else
               if (ljt.le.lit) go to 60
               go to 40
            end if
         else if (jnek) then
            if (.not.(jnel)) go to 70
            if (.not.(knel)) go to 80
c     iat # jat # kat # lat  -- omit one centre
            min = lit
            imin = 1
            do 30 iper = 2 , 4
               if (lll(iper).lt.min) then
                  min = lll(iper)
                  imin = iper
               end if
 30         continue
            go to (90,100,110,120) , imin
            go to 90
         else
            if (.not.(jnel)) go to 60
c     ----- jat = kat # iat # lat -----
            oskip(1) = .false.
            oskip(4) = .false.
            natomd(1) = iat
            natomd(2) = lat
            natomd(3) = jat
            npass = 2
            go to 130
         end if
      else if (inek) then
         if (.not.(knel)) then
c     iat=jat  ,  kat=lat   differentiate one pair
            n1 = lit*ljt*(lkt+1)*(llt+1)
            n2 = (lit+1)*(ljt+1)*lkt*llt
            if (n2.lt.n1) go to 80
         end if
      else
         if (inel) then
c     iat=jat=kat   derivative (ij/kl')
            oskip(4) = .false.
            natomd(1) = lat
            natomd(2) = iat
            npass = 1
         end if
         go to 130
      end if
c     iat=jat   derivatives (ij/k'l) and (ij/kl')
      oskip(3) = .false.
      oskip(4) = .false.
      natomd(1) = kat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     iat=kat=lat   derivative (ij'/kl)
 40   oskip(2) = .false.
      natomd(1) = jat
      natomd(2) = iat
      npass = 1
      go to 130
c     iat=kat   derivatives (ij'/kl) and (ij/kl')
 50   oskip(2) = .false.
      oskip(4) = .false.
      natomd(1) = jat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     jat=kat=lat    (i'j/kl)
 60   oskip(1) = .false.
      natomd(1) = iat
      natomd(2) = jat
      npass = 1
      go to 130
c      jat=lat    derivatives (i'j/kl) and (ij/k'l)
 70   oskip(1) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = kat
      natomd(3) = jat
      npass = 2
      go to 130
c     kat=lat   derivatives (i'j/kl) and (ij'/kl)
 80   oskip(1) = .false.
      oskip(2) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      npass = 2
      go to 130
 90   natomd(1) = jat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = iat
      npass = 3
      oskip(2) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 100  natomd(1) = iat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = jat
      npass = 3
      oskip(1) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 110  natomd(1) = iat
      natomd(2) = jat
      natomd(3) = lat
      natomd(4) = kat
      npass = 3
      oskip(1) = .false.
      oskip(2) = .false.
      oskip(4) = .false.
      go to 130
 120  oskip(1) = .false.
      oskip(2) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      natomd(4) = lat
      npass = 3
c     -----
 130  if (.not.outd) return
      write (iw,6010) ii , jj , kk , ll , 
     +    oskip(1) , oskip(2) , oskip(3) , oskip(4) , 
     +    npass , (natomd(i),i=1,4)
      return
 6010 format (/,' ********  ii,jj,kk,ll =',4i3,' skip1,2,3,4 =',4l3,
     +        ' npass =',i2,' centers =',4i5,/)
      end
_ENDIF
      subroutine ssdss(qq)
      implicit REAL  (a-h,o-z)
INCLUDE(common/sizes)
      parameter (mxp2 = mxprms*mxprms)
INCLUDE(common/dshlnf)
      parameter (ncmax=65)
      common/junk/ddij(ncmax,225),ddkl(ncmax,225),aei(ncmax),aej(ncmax)
     + ,aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
INCLUDE(common/incrd)
INCLUDE(common/dmisc)
      dimension qq(*)
      data pt5,pt2/0.5d0,0.2d0/
      data third/0.33333333333333333d0/
      data one/1.0d0/
      data pi252/34.986836655250d0/
      data pie4 /7.85398163397448d-01/
      do 20 i = 1 , 12
         qq(ic7+i) = 0.0d0
 20   continue
c
c     lgmax = ngd
      do 150 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*xk
         akyk = ak*yk
         akzk = ak*zk
         csk = csc(kg)*pi252
c        if (okandl) lgmax = kg
         do 140 lg = 1 , ngd
            al = dg(lg)
            b = ak + al
            b1 = one/b
            bbrrk = al*brrk*b1
            if ((bbrrk+rsmall).le.tol1) then
               exkl = dexp(-bbrrk)*csd(lg)*csk*b1
               xb = (akxk+al*xl)*b1
               yb = (akyk+al*yl)*b1
               zb = (akzk+al*zl)*b1
               nn = 0
               n = 0
               jgmax = ngb
               do 130 ig = 1 , nga
                  ai = ag(ig)
                  if (oiandj) jgmax = ig
                  do 120 jg = 1 , jgmax
                     n = n + 1
                     aj = bg(jg)
                     a = ai + aj
                     if ((bbrrk+r(n)).le.tol2) then
                        ab = a*b
                        aandb = a + b
                        expe = dd(nn+1)*exkl*exij(n)/dsqrt(aandb)
                        if (dabs(expe).ge.tol4) then
                           rho = ab/aandb
                           xa = x1(n)
                           ya = y1(n)
                           za = z1(n)
                           x = rho*((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
                           if (x.gt.5.0d0) then
                              xinv = one/x
                              if (x.le.15.0d0) then
                                 g = dexp(-x)
                                 if (x.gt.10.0d0) then
                                    ww1 = (((-1.8784686463512d-01*xinv+
     +                                 2.2991849164985d-01)
     +                                 *xinv-4.9893752514047d-01)
     +                                 *xinv-2.1916512131607d-05)
     +                                 *g + dsqrt(pie4*xinv)
                                    f1 = (ww1-g)*xinv*pt5
                                    rt1 = f1/(ww1-f1)
                                 else
                                    ww1 =
     +                                 ((((((4.6897511375022d-01*xinv-
     +                                 6.9955602298985d-01)
     +                                 *xinv+5.3689283271887d-01)
     +                                 *xinv-3.2883030418398d-01)
     +                                 *xinv+2.4645596956002d-01)
     +                                 *xinv-4.9984072848436d-01)
     +                                 *xinv-3.1501078774085d-06)
     +                                 *g + dsqrt(pie4*xinv)
                                    f1 = (ww1-g)*xinv*pt5
                                    rt1 = f1/(ww1-f1)
                                 end if
                              else if (x.gt.33.0d0) then
                                 ww1 = dsqrt(pie4*xinv)
                                 rt1 = pt5/(x-pt5)
                              else
                                 g = dexp(-x)
                                 ww1 = ((1.9623264149430d-01*xinv-
     +                                 4.9695241464490d-01)
     +                                 *xinv-6.0156581186481d-05)
     +                                 *g + dsqrt(pie4*xinv)
                                 f1 = (ww1-g)*xinv*pt5
                                 rt1 = f1/(ww1-f1)
                              end if
                           else if (x.gt.1.0d0) then
                              if (x.gt.3.0d0) then
                                 y = x - 4.0d0
                                 f1 = ((((((((((-2.62453564772299d-11*y+
     +                                3.24031041623823d-10)
     +                                *y-3.614965656163d-09)
     +                                *y+3.760256799971d-08)
     +                                *y-3.553558319675d-07)
     +                                *y+3.022556449731d-06)
     +                                *y-2.290098979647d-05)
     +                                *y+1.526537461148d-04)
     +                                *y-8.81947375894379d-04)
     +                                *y+4.33207949514611d-03)
     +                                *y-1.75257821619926d-02)
     +                                *y + 5.28406320615584d-02
                                 ww1 = (x+x)*f1 + dexp(-x)
                                 rt1 = f1/(ww1-f1)
                              else
                                 y = x - 2.0d0
                                 f1 = ((((((((((-1.61702782425558d-10*y+
     +                                1.96215250865776d-09)
     +                                *y-2.14234468198419d-08)
     +                                *y+2.17216556336318d-07)
     +                                *y-1.98850171329371d-06)
     +                                *y+1.62429321438911d-05)
     +                                *y-1.16740298039895d-04)
     +                                *y+7.24888732052332d-04)
     +                                *y-3.79490003707156d-03)
     +                                *y+1.61723488664661d-02)
     +                                *y-5.29428148329736d-02)
     +                                *y + 1.15702180856167d-01
                                 ww1 = (x+x)*f1 + dexp(-x)
                                 rt1 = f1/(ww1-f1)
                              end if
                           else if (x.gt.3.0d-07) then
                              f1 = ((((((((-8.36313918003957d-08*x+
     +                             1.21222603512827d-06)
     +                             *x-1.15662609053481d-05)
     +                             *x+9.25197374512647d-05)
     +                             *x-6.40994113129432d-04)
     +                             *x+3.78787044215009d-03)
     +                             *x-1.85185172458485d-02)
     +                             *x+7.14285713298222d-02)
     +                             *x-1.99999999997023d-01)
     +                             *x + 3.33333333333318d-01
                              ww1 = (x+x)*f1 + dexp(-x)
                              rt1 = f1/(ww1-f1)
                           else
                              rt1 = pt5 - x*pt2
                              ww1 = one - x*third
                           end if
                           u2 = rt1*rho
                           f00 = expe*ww1/(ab+u2*aandb)
                           joff = ic7
                           do 110 iper = 1 , 4
                              if (oskip(iper)) go to 110
                              go to (30,40,50,60) , iper
 30                           alpha = ai + ai
                              xd = xi
                              yd = yi
                              zd = zi
                              go to 70
 40                           alpha = aj + aj
                              xd = xj
                              yd = yj
                              zd = zj
                              go to 70
 50                           alpha = ak + ak
                              xd = xk
                              yd = yk
                              zd = zk
                              go to 70
 60                           alpha = al + al
                              xd = xl
                              yd = yl
                              zd = zl
 70                           bxbd = b*(xb-xd)
                              bybd = b*(yb-yd)
                              bzbd = b*(zb-zd)
                              axad = a*(xa-xd)
                              ayad = a*(ya-yd)
                              azad = a*(za-zd)
                              c1x = bxbd + axad
                              c1y = bybd + ayad
                              c1z = bzbd + azad
                              go to (90,90,80,80) , iper
 80                           c2x = a*bxbd
                              c2y = a*bybd
                              c2z = a*bzbd
                              go to 100
 90                           c2x = b*axad
                              c2y = b*ayad
                              c2z = b*azad
 100                          qq(joff+1) = qq(joff+1) + (u2*c1x+c2x)
     +                           *f00*alpha
                              qq(joff+2) = qq(joff+2) + (u2*c1y+c2y)
     +                           *f00*alpha
                              qq(joff+3) = qq(joff+3) + (u2*c1z+c2z)
     +                           *f00*alpha
                              joff = joff + 3
 110                       continue
                        end if
                     end if
                     nn = nn + 4
 120              continue
 130           continue
            end if
 140     continue
 150  continue
      if (natomd(1).ne.natomd(2)) return
      qq(ic7+1) = qq(ic7+1) + qq(ic7+4)
      qq(ic7+2) = qq(ic7+2) + qq(ic7+5)
      qq(ic7+3) = qq(ic7+3) + qq(ic7+6)
      qq(ic7+4) = qq(ic7+7)
      qq(ic7+5) = qq(ic7+8)
      qq(ic7+6) = qq(ic7+9)
      natomd(2) = natomd(3)
      natomd(3) = natomd(4)
      natomd(4) = 0
      npass = npass - 1
      return
      end
      subroutine subsd(x,y,z,xd,yd,zd,
     *a,m1,m2,m3,m4,i1,i2,k1,k2,ncdim)
      implicit REAL  (a-h,o-z)
INCLUDE(common/incrd)
INCLUDE(common/dshlno)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),
     *         xd(ncdim,*),yd(ncdim,*),zd(ncdim,*)
      parameter (ncmax=65)
      dimension a(ncmax)
c
      do 20 i = 1 , ncontr
         a(i) = a(i) + a(i)
 20   continue
      n1 = 1
c
      do 120 i = 1 , m1
         n2 = n1
         do 110 j = 1 , m2
            n3 = n2
            do 100 k = 1 , m3
               n4 = n3
               do 90 l = 1 , m4
                  go to (30,50,70,70,70,70,70) , l
c
 30               do 40 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2)
 40               continue
                  n4 = n4 + k2
                  go to 90
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 50               do 60 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) - x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) - y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) - z(nr,n4-k2)
 60               continue
                  n4 = n4 + k2
                  go to 90
 70               fac = -dfloat(l-1)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 80 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) + fac*x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) + fac*y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) + fac*z(nr,n4-k2)
 80               continue
c
                  n4 = n4 + k2
 90            continue
               n3 = n3 + k1
 100        continue
            n2 = n2 + i2
 110     continue
         n1 = n1 + i1
 120  continue
      return
      end
      subroutine symde(ict,natoms)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/symtry)
INCLUDE(common/grad2)
_IF(rpagrad)
INCLUDE(common/rpadcom)
_ENDIF
      common/junk/ptr(3,144)
      dimension ict(natoms,*)
      data dzero,done /0.0d0,1.0d0/
c
c     ----- read in tranformation matrices of coordinates. -----
c
c     ----- and
c
c     ----- read in transformation table- atoms versus symmetry operatio
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      nav = lenwrd()
      call readi(ict,nw196(6)*nav,ibl196(6),idaf)
c
c     ----- symmetryze gradient vector -----
c
      do 120 ic = 1 , nat
         do 50 it = 1 , nt
            if (ict(ic,it).gt.ic) go to 120
 50      continue
c
c ps ----start-----
c
c the following code zeros the gradient on atoms with no shells 
c it is suppressed as these gradients are needed for hybrid QM/MM
c calculations
c
c         do 60 ii = 1 , nshell
c            if (katom(ii).eq.ic) go to 80
c 60      continue
c         do 70 it = 1 , nt
c            knu = ict(ic,it)
c            de(1,knu) = 0.0d0
c            de(2,knu) = 0.0d0
c            de(3,knu) = 0.0d0
c 70      continue
c
c  ps ---end----
c
 80      dedx = dzero
         dedy = dzero
         dedz = dzero
         do 90 it = 1 , nt
            icnu = ict(ic,it)
            dedxp = de(1,icnu)
            dedyp = de(2,icnu)
            dedzp = de(3,icnu)
            n = 3*(it-1)
            dedx = dedx + dedxp*ptr(1,n+1) + dedyp*ptr(2,n+1)
     +             + dedzp*ptr(3,n+1)
            dedy = dedy + dedxp*ptr(1,n+2) + dedyp*ptr(2,n+2)
     +             + dedzp*ptr(3,n+2)
            dedz = dedz + dedxp*ptr(1,n+3) + dedyp*ptr(2,n+3)
     +             + dedzp*ptr(3,n+3)
 90      continue
         de(1,ic) = dedx
         de(2,ic) = dedy
         de(3,ic) = dedz
         do 110 it = 1 , nt
            icnu = ict(ic,it)
            if (icnu.ne.ic) then
               if (it.ne.nt) then
                  it1 = it + 1
                  do 100 jt = it1 , nt
                     if (ict(ic,jt).eq.icnu) go to 110
 100              continue
               end if
               jt = invt(it)
               n = 3*(jt-1)
               de(1,icnu) = de(1,ic)*ptr(1,n+1) + de(2,ic)*ptr(2,n+1)
     +                      + de(3,ic)*ptr(3,n+1)
               de(2,icnu) = de(1,ic)*ptr(1,n+2) + de(2,ic)*ptr(2,n+2)
     +                      + de(3,ic)*ptr(3,n+2)
               de(3,icnu) = de(1,ic)*ptr(1,n+3) + de(2,ic)*ptr(2,n+3)
     +                      + de(3,ic)*ptr(3,n+3)
            end if
 110     continue
 120  continue
      dum = done/dfloat(nt)
      do 140 n = 1 , nat
         do 130 i = 1 , 3
            de(i,n) = de(i,n)*dum
 130     continue
 140  continue
_IF(rpagrad)
c
c     ----- symmetryze gradient vector -----
c
      if (orpagrad) then
         do ist = 1, nrpastate
            do 1120 ic = 1 , nat
               do it = 1 , nt
                  if (ict(ic,it).gt.ic) go to 1120
               enddo
               dedx = dzero
               dedy = dzero
               dedz = dzero
               do it = 1 , nt
                  icnu = ict(ic,it)
                  dedxp = rpade(1,icnu,ist)
                  dedyp = rpade(2,icnu,ist)
                  dedzp = rpade(3,icnu,ist)
                  n = 3*(it-1)
                  dedx = dedx + dedxp*ptr(1,n+1) + dedyp*ptr(2,n+1)
     +                 + dedzp*ptr(3,n+1)
                  dedy = dedy + dedxp*ptr(1,n+2) + dedyp*ptr(2,n+2)
     +                 + dedzp*ptr(3,n+2)
                  dedz = dedz + dedxp*ptr(1,n+3) + dedyp*ptr(2,n+3)
     +                 + dedzp*ptr(3,n+3)
               enddo
               rpade(1,ic,ist) = dedx
               rpade(2,ic,ist) = dedy
               rpade(3,ic,ist) = dedz
               do 1110 it = 1 , nt
                  icnu = ict(ic,it)
                  if (icnu.ne.ic) then
                     if (it.ne.nt) then
                        it1 = it + 1
                        do jt = it1 , nt
                           if (ict(ic,jt).eq.icnu) go to 1110
                        enddo
                     end if
                     jt = invt(it)
                     n = 3*(jt-1)
                     rpade(1,icnu,ist) = rpade(1,ic,ist)*ptr(1,n+1) 
     +                                 + rpade(2,ic,ist)*ptr(2,n+1)
     +                                 + rpade(3,ic,ist)*ptr(3,n+1)
                     rpade(2,icnu,ist) = rpade(1,ic,ist)*ptr(1,n+2) 
     +                                 + rpade(2,ic,ist)*ptr(2,n+2)
     +                                 + rpade(3,ic,ist)*ptr(3,n+2)
                     rpade(3,icnu,ist) = rpade(1,ic,ist)*ptr(1,n+3) 
     +                                 + rpade(2,ic,ist)*ptr(2,n+3)
     +                                 + rpade(3,ic,ist)*ptr(3,n+3)
                  end if
 1110          continue
 1120       continue
         enddo
         dum = done/dfloat(nt)
         do ist = 1, nrpastate
            do n = 1 , nat
               do i = 1 , 3
                  rpade(i,n,ist) = rpade(i,n,ist)*dum
               enddo
            enddo
         enddo
      endif
_ENDIF
c
      return
      end
      subroutine tpdder(tpdm,lentpd,da,db,ltri,ndens,ii,jj,kk,ll,q4)
      implicit REAL  (a-h,o-z)
c
c     constructs scf tpdm's from cross products of two
c     sets of one particle density matrices. first , da ,
c     may be thought of as the unperturbed dm, second , in
c     db , may be thought of as a set of derivative density
c     matrices.  intended for use in property derivative
c     calculations
c
c     shells (ii,jj,kk,ll) : q4 =degeneracy (symmetry) factor
c
      logical ijeq,kleq
INCLUDE(common/sizes)
      dimension tpdm(lentpd,ndens),da(ltri),db(ltri,ndens)
INCLUDE(common/nshel)
INCLUDE(common/incrd)
INCLUDE(common/mapper)
c
      mini = kmin(ii)
      minj = kmin(jj)
      mink = kmin(kk)
      minl = kmin(ll)
      maxi = kmax(ii)
      maxj = kmax(jj)
      maxk = kmax(kk)
      maxl = kmax(ll)
      loci = kloc(ii) - mini
      locj = kloc(jj) - minj
      lock = kloc(kk) - mink
      locl = kloc(ll) - minl
      ni = 1
      do 60 i = mini , maxi
         nj = ni
         i1 = loci + i
         ii1 = iky(i1)
         do 50 j = minj , maxj
            nk = nj
            j1 = locj + j
            jj1 = iky(j1)
            if (i1.ge.j1) then
               mij = ii1 + j1
            else
               mij = jj1 + i1
            end if
            ijeq = i1.eq.j1
            do 40 k = mink , maxk
               nl = nk
               k1 = lock + k
               kk1 = iky(k1)
               if (i1.ge.k1) then
                  mik = ii1 + k1
               else
                  mik = kk1 + i1
               end if
               if (j1.ge.k1) then
                  mjk = jj1 + k1
               else
                  mjk = kk1 + j1
               end if
               do 30 l = minl , maxl
                  nn = nl
                  l1 = locl + l
                  ll1 = iky(l1)
                  if (k1.ge.l1) then
                     mkl = kk1 + l1
                  else
                     mkl = ll1 + k1
                  end if
                  if (j1.ge.l1) then
                     mjl = jj1 + l1
                  else
                     mjl = ll1 + j1
                  end if
                  if (i1.ge.l1) then
                     mil = ii1 + l1
                  else
                     mil = ll1 + i1
                  end if
                  kleq = k1.eq.l1
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                  do 20 ndd = 1 , ndens
                     dfac = (da(mij)*db(mkl,ndd)+da(mkl)*db(mij,ndd))
     +                      *4.0d0 - da(mik)*db(mjl,ndd) - da(mjl)
     +                      *db(mik,ndd) - da(mjk)*db(mil,ndd) - da(mil)
     +                      *db(mjk,ndd)
                     if (ijeq) dfac = dfac*0.5d0
                     if (kleq) dfac = dfac*0.5d0
                     tpdm(nn,ndd) = tpdm(nn,ndd) + dfac*q4
 20               continue
                  nl = nl + inc2
 30            continue
               nk = nk + inc3
 40         continue
            nj = nj + inc4
 50      continue
         ni = ni + inc5
 60   continue
      return
      end
      subroutine ver_drv2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/drv2e.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
