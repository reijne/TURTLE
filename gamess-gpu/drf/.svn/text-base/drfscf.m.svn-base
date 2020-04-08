      subroutine drfhamc(f,d,s,dx,dy,dz,iexpc,
     1                   omegas,ijbits,omegao,ijbito)
c------
c      adds reaction field integrals to the fock matrix in the
c      rhf scf procedure
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension f(nchd),d(nchd)
      dimension iexpc(nchd),ijbits(nchd),ijbito(nchd)
      dimension s(nchd),dx(nchd),dy(nchd),dz(nchd)
      dimension omegas(nwtc,nwtc),omegao(nwtc,nwtc)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/opt)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfind)
c
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/shlint)
c
      character*8 direct
c
      dimension sump(4), sumq(4)
c
      data direct /'direct  '/
c
      data zero,pt5,one/0.0d00,0.5d00,1.0d00/
      data tenm1, tenm2 /1.0d-01,1.0d-02/
c
c-----  begin
c
c     --------  expanded external field
c     --------   memory
c
      integr = 0
      ieps = 0
      fexch = gamdrf
      l2 = num*(num+1)/2
      acurc = tenm2*acur
c
c-----  get previous incomplete drf 2-el.hamiltonian
c
c     if (hndtyp .ne. direct) then
cxxx  if ((iodadrf(1,72) .ne. 0) .and.
cxxx 1    (hndtyp .ne. direct)) then
c     if (iodadrf(1,72) .ne. 0) then
c       call daread(idafdrf,iodadrf,f,l2,72)
c     else
c       call clear(f,l2)
c     endif
c     else
      do i = 1, num
c       do j = 1, i-1
c         ij = ia(i) + j
          ij = ia(i) + i
          f(ij) = pt5*f(ij)
c         f(ij) = 2.0*f(ij)
c       enddo
      enddo
c     endif
c
c-----  get difference density
c
c     if (hndtyp .ne. direct) call daread(idafdrf,iodadrf,d,l2,70)
c
c-----  loop over charge distributions ij
c
      do 1700, i = 1, num
c 1-----
        do 1600, j = 1, i
c   2-----
          ij = ia(i) + j
c
          if ((ijbits(ij) .eq. 0) .and. ijbito(ij) .eq. 0) goto 1600
c
          p(1) = dx(ij)
          p(2) = dy(ij)
          p(3) = dz(ij)
          p(4) = s(ij)
          ijexp = iexpc(ij)
c
c    -----  loop over charge distributions kl
c
          do 1500 k = 1, i
c     3-----
            ik = ia(i) + k
            jk = ia(max(j,k)) + min(j,k)
            lmax = k
            if (k .eq. i) lmax = j
            do 1400, l = 1, lmax
c       4-----
              kl = ia(k) + l
              if ((ijbits(kl) .eq. 0) .and. (ijbito(kl) .eq. 0))
     1          goto 1400
              il = ia(i) + l
              jl = ia(max(j,l)) + min(j,l)
c
              q(1) = dx(kl)
              q(2) = dy(kl)
              q(3) = dz(kl)
              q(4) = s(kl)
              klexp = iexpc(kl)
c
c        -----  total dielectric contribution
c
              if ((ijbits(ij) .eq. 0) .or. (ijbits(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegas)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegas)
              call matvec(omgab,q,sumq,4,.false.)
c
              goto 1310
c
 1300         continue
c
c        -----  optic dielectric contribution
c
              if ((ijbito(ij) .eq. 0) .or. (ijbito(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegao)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegao)
              call matvec(omgab,q,sumq,4,.false.)
c
c        -----  reaction field contribution
c
 1310         continue
c
c      --------  p.th. van duijnen, groningen 1990 --------
c
              val = zero
              if (
     1          abs(d(ij)) .gt. acurc .or.
     2          abs(d(kl)) .gt. acurc .or.
     3          abs(d(ik)) .gt. acurc .or.
     4          abs(d(jl)) .gt. acurc .or.
     5          abs(d(il)) .gt. acurc .or.
     6          abs(d(jk)) .gt. acurc ) then
c         5-----
                integr = integr + 1
c
c          -----  note: + sign because source and recipient are electron
c
                do 20, ii = 1, 4
                  val = val + sump(ii)*q(ii)
                  val = val + sumq(ii)*p(ii)
   20           continue
                if(idrfout.eq.5) write (iwr,*) i,j,k,l,val
c         5-----
              endif
c
              if (abs(val) .lt. cutoff) go to 1330
              if (i .eq. j) val = val*pt5
              if (k .eq. l) val = val*pt5
              if ((i .eq. k) .and. (j .eq. l)) val = val*pt5
cxxx          val4 = (val + val) + (val + val)
              val4 = (val + val) 
c
c        -----  skip coulomb contributions for optic dielectric
c
              if ((itwoeps .eq. 1) .and. (ieps .eq. 1)) goto 1320
c
  250         f(ij) = f(ij) + val4*d(kl)
              f(kl) = f(kl) + val4*d(ij)
c
c        -----  skip part of exchange contributions for total dielectric
c
              if ((itwoeps .eq. 1) .and. (ieps .eq. 0)) goto 1330
c
 1320         continue
c
              if (fexch .ne. zero) then
cextra
                val = pt5*val
                f(ik) = f(ik) - val*d(jl)*fexch
                f(il) = f(il) - val*d(jk)*fexch
                f(jk) = f(jk) - val*d(il)*fexch
                f(jl) = f(jl) - val*d(ik)*fexch
              endif
c
 1330         if ((itwoeps .eq. 1) .and. (ieps .eq. 0)
     1             .and. (fexch .ne. zero)) then
                ieps = 1
                goto 1300
              else
                ieps = 0
              endif
c       4-----  next l
 1400       continue
c     3-----  next k
 1500     continue
c   2-----  next j
 1600   continue
c 1-----  next i
 1700 continue
c
      if (idrfout .eq. 5) write(iwr,*) 'drf ints   ',integr
c
c-----  save incomplete fock contribution
c
c     if (hndtyp .ne. 'direct')
c    1 call dawrit(idafdrf,iodadrf,f,l2,72,navdrf)
c
c-----  correct fock matrix
c
c-----  cf hstar, hdir
c
      do 2000, i = 1, num
c       do 2000, j = 1, i-1
c         ij = ia(i) + j
c         f(ij) = pt5*f(ij)
          ij = ia(i) + i
          f(ij) = 2.0d0*f(ij)
 2000 continue
c
      return
      end
      subroutine drfhamu(fa,da,fb,db,s,dx,dy,dz,iexpc,
     1                   omegas,ijbits,omegao,ijbito)
c------
c      adds reaction field integrals to the fock matrix in the
c      uhf scf procedure
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension fa(nchd),da(nchd),fb(nchd),db(nchd)
      dimension s(nchd),dx(nchd),dy(nchd),dz(nchd)
      dimension iexpc(nchd),ijbits(nchd),ijbito(nchd)
      dimension omegas(nwtc,nwtc), omegao(nwtc,nwtc)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/opt)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfind)
c
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/shlint)
c
c
      character*8 direct
c
      dimension sump(4), sumq(4)
c
      data direct /'direct  '/
c
      data zero,pt5,one/0.0d00,0.5d00,1.0d00/
      data tenm1, tenm2 /1.0d-01,1.0d-02/
c
c-----  begin
c     --------  expanded external field
c     --------   memory
c
      fexch = gamdrf
      ieps = 0
      integr = 0
      l2 = num*(num+1)/2
      acurc = tenm2*acur
c
c-----  read previous incomplete fock matrices
c
c     if (hndtyp .ne. 'direct') then
c       if (iodadrf(1,72) .ne. 0) then
cxxx  if ((iodadrf(1,72) .ne. 0) .and.
cxxx 1    (hndtyp .ne. direct)) then
c       call daread(idafdrf,iodadrf,fa,l2,72)
c       call daread(idafdrf,iodadrf,fb,l2,73)
c     else
c       call clear(fa,l2)
c       call clear(fb,l2)
c     endif
c     else
      do i = 1, num
cxxx    do j = 1, i-1
cxxx      ij = ia(i) + j
          ij = ia(i) + i
cxxx      fa(ij) = 2.0*fa(ij)
cxxx      fb(ij) = 2.0*fb(ij)
          fa(ij) = pt5*fa(ij)
          fb(ij) = pt5*fb(ij)
cxxx    enddo
      enddo
c     endif
c
c-----  read difference densities
c       only if not in direct 2-electron integral mode
c
c     if (hndtyp .ne. direct) then
c       call daread(idafdrf,iodadrf,da,l2,70)
c       call daread(idafdrf,iodadrf,db,l2,71)
c     endif
c
      do 1700, i = 1, num
c 1-----
        do 1600, j = 1, i
c   2-----
          ij = ia(i) + j
c
          if ((ijbits(ij) .eq. 0) .and. (ijbito(ij) .eq. 0)) goto 1600
c
          ijexp = iexpc(ij)
c
          p(1) = dx(ij)
          p(2) = dy(ij)
          p(3) = dz(ij)
          p(4) = s(ij)
c
c    -----  loop over charge distributions kl
c
          do 1500 k = 1, i
c     3-----
            ik = ia(i) + k
            jk = ia(max(j,k)) + min(j,k)
            lmax = k
            if (k .eq. i) lmax = j
            do 1400, l = 1, lmax
c       4-----
              kl = ia(k) + l
              if ((ijbits(kl) .eq. 0) .and. (ijbito(kl) .eq. 0))
     1          goto 1400
              il = ia(i) + l
              jl = ia(max(j,l)) + min(j,l)
c
              q(1) = dx(kl)
              q(2) = dy(kl)
              q(3) = dz(kl)
              q(4) = s(kl)
              klexp = iexpc(kl)
c
c        -----  total dielectric contribution
c
              if ((ijbits(ij) .eq. 0) .and. (ijbits(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegas)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegas)
              call matvec(omgab,q,sumq,4,.false.)
c
              goto 1310
c
 1300         continue
c
c        -----  optic dielectric contribution
c
              if ((ijbito(ij) .eq. 0) .and. (ijbito(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegao)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegao)
              call matvec(omgab,q,sumq,4,.false.)
c
c        -----  reaction field contribution
c
 1310         continue
c
c      --------  p.th. van duijnen, groningen 1990 --------
c
              val = zero
              if (
     1          abs(da(ij)) .gt. acurc .or.
     2          abs(da(kl)) .gt. acurc .or.
     3          abs(da(ik)) .gt. acurc .or.
     4          abs(da(jl)) .gt. acurc .or.
     5          abs(da(il)) .gt. acurc .or.
     6          abs(da(jk)) .gt. acurc .or.
     1          abs(db(ij)) .gt. acurc .or.
     2          abs(db(kl)) .gt. acurc .or.
     3          abs(db(ik)) .gt. acurc .or.
     4          abs(db(jl)) .gt. acurc .or.
     5          abs(db(il)) .gt. acurc .or.
     6          abs(db(jk)) .gt. acurc ) then
c         5-----
                integr = integr + 1
c
c          -----  note: + -sign because source and recipient are electro
c
                do 20, ii = 1, 4
                  val = val + sump(ii)*q(ii)
                  val = val + sumq(ii)*p(ii)
   20           continue
                if (idrfout .eq. 5) write (iwr,*) i,j,k,l,val
c         5-----
              endif
c
              if (abs(val) .lt. cutoff) go to 1330
              if (i .eq. j) val = val*pt5
              if (k .eq. l) val = val*pt5
              if ((i .eq. k) .and. (j .eq. l)) val = val*pt5
cxxx          val4 = (val+val)+(val+val)
              val4 = (val+val)
c
c        -----  skip coulomb contributions for optic dielectric
c
              if ((itwoeps .eq. 1) .and. (ieps .eq. 1)) goto 1320
c
              dum = val4*(da(kl)+db(kl))
              fa(ij) = fa(ij) + dum
              fb(ij) = fb(ij) + dum
              dum = val4*(da(ij)+db(ij))
              fa(kl) = fa(kl) + dum
              fb(kl) = fb(kl) + dum
c
c        -----  skip part of exchange contributions for total dielectric
c
 1320         if ((itwoeps .eq. 1) .and. (ieps .eq. 0)) goto 1330
c
              if (fexch .ne. zero) then
c         5-----
cxxx            val2 = val + val
                val2 = val
                fa(ik) = fa(ik) - val2*da(jl)*fexch
                fb(ik) = fb(ik) - val2*db(jl)*fexch
                fa(il) = fa(il) - val2*da(jk)*fexch
                fb(il) = fb(il) - val2*db(jk)*fexch
                fa(jk) = fa(jk) - val2*da(il)*fexch
                fb(jk) = fb(jk) - val2*db(il)*fexch
                fa(jl) = fa(jl) - val2*da(ik)*fexch
                fb(jl) = fb(jl) - val2*db(ik)*fexch
c         5-----
              endif
c
 1330         if ((itwoeps .eq. 1) .and. (ieps .eq. 0)
     1            .and. (fexch .ne. zero)) then
                ieps = 1
                goto 1300
              else
                ieps = 0
              endif
c       4-----  next l
 1400      continue
c     3-----  next k
 1500     continue
c   2-----  next j
 1600   continue
c 1-----  next i
 1700 continue
c
c-----  save incomplete fock matrices
c
c     if (hndtyp .ne. 'direct') then
c     call dawrit(idafdrf,iodadrf,fa,l2,72,navdrf)
c     call dawrit(idafdrf,iodadrf,fb,l2,73,navdrf)
c     endif
c
c-----  correct off-diagonal elements
c
      do 2000, i = 1, num
cxxx    do 2000, j = 1, i-1
cxxx      ij = ia(i) + j
          ij = ia(i) + i
cxxx      fa(ij) = pt5*fa(ij)
cxxx      fb(ij) = pt5*fb(ij)
          fa(ij) = 2.0d0*fa(ij)
          fb(ij) = 2.0d0*fb(ij)
 2000 continue
c
      return
      end
      subroutine drfhamo(f,d,s,dx,dy,dz,iexpc,
     1                   omegas,ijbits,omegao,ijbito,nset)
c------
c      adds reaction field integrals to the fock matrix in the
c      rhfop scf procedure
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension f(nchd),d(nchd)
      dimension iexpc(nchd),ijbits(nchd),ijbito(nchd)
      dimension s(nchd),dx(nchd),dy(nchd),dz(nchd)
      dimension omegas(nwtc,nwtc),omegao(nwtc,nwtc)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfind)
c
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/shlint)
c
      dimension sump(4), sumq(4)
c
      data zero,pt5,one /0.0d00,0.5d00,1.0d00/
      data tenm1, tenm2 /1.0d-01,1.0d-02/
c
c-----  begin
c
      fexch = gamdrf
      integr = 0
      ieps = 0
      l2 = num*(num+1)/2
      acurc = tenm2*acur
c
c-----  read previous incomplete fock matrix
c
      if (iodadrf(1,72) .ne. 0) then
        call daread(idafdrf,iodadrf,f,nset*l2,72)
      else
        call clear(f,nset*l2)
      endif
c
c-----  read difference density
c
      call daread(idafdrf,iodadrf,d,nset*l2,70)
c
c-----  loop over charge distributions ij
c
      do 1700, i = 1, num
c 1-----
        do 1600, j = 1, i
c   2-----
          ij = ia(i) + j
          ijexp = iexpc(ij)
c
          if ((ijbits(ij) .eq. 0) .and. ijbito(ij) .eq. 0) goto 1600
c
          p(1) = dx(ij)
          p(2) = dy(ij)
          p(3) = dz(ij)
          p(4) = s(ij)
c
c    -----  loop over charge distributions kl
c
          do 1500 k = 1, i
c     3-----
            ik = ia(i) + k
            jk = ia(max(j,k)) + min(j,k)
            lmax = k
            if (k .eq. i) lmax = j
            do 1400, l = 1, lmax
c       4-----
              kl = ia(k) + l
              if ((ijbits(kl) .eq. 0) .and. (ijbito(kl) .eq. 0))
     1          goto 1400
              il = ia(i) + l
              jl = ia(max(j,l)) + min(j,l)
c
              q(1) = dx(kl)
              q(2) = dy(kl)
              q(3) = dz(kl)
              q(4) = s(kl)
              klexp = iexpc(kl)
c
c        -----  total dielectric contribution
c
              if ((ijbits(ij) .eq. 0) .and. (ijbits(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegas)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegas)
              call matvec(omgab,q,sumq,4,.false.)
c
              goto 1310
c
 1300         continue
c
c        -----  optic dielectric contribution
c
              if ((ijbito(ij) .eq. 0) .and. (ijbito(kl) .eq. 0))
     1        goto 1330
c
              call drfoab(klexp,ijexp,nwtc,omegao)
              call matvec(omgab,p,sump,4,.false.)
              call drfoab(ijexp,klexp,nwtc,omegao)
              call matvec(omgab,q,sumq,4,.false.)
c
c        -----  reaction field contribution
c
 1310         continue
c
              val = zero
              dmx = max(abs(d(ij)), abs(d(kl)))
              imx = 0
              do 125, iset = 1, nset
                dmx = max(dmx,abs(d(jl+imx)),abs(d(jk+imx)),
     1                        abs(d(il+imx)),abs(d(ik+imx)))
  125         continue
              if (dmx .gt. acurc) then
c         5-----
                integr = integr + 1
                do 150, ii = 1, 4
                   val = val - sump(ii)*q(ii)
                   val = val - sumq(ii)*p(ii)
  150           continue
                if (idrfout .eq. 5) write (iwr,*) i,j,k,l,val
c         5-----
              endif
c
              if (abs(val) .lt. cutoff) go to 1330
              if (i .eq. j) val = val*pt5
              if (k .eq. l) val = val*pt5
              if ((i .eq. k) .and. (j .eq. l)) val = val*pt5
              val4 = (val+val) + (val+val)
c
c        -----  skip coulomb contributions for optic dielectric
c
              if ((itwoeps .eq. 1) .and. (ieps .eq. 1)) goto 1320
c
  250         f(ij) = f(ij) + val4*d(kl)
              f(kl) = f(kl) + val4*d(ij)
c
c        -----  skip part of exchange contributions for total dielectric
c
 1320         if ((itwoeps .eq. 1) .and. (ieps .eq. 0)) goto 1330
c
              if (fexch .ne. zero) then
c         5-----
                imx = 0
                do 300, iset = 1, nset
c           6-----
                  f(ik+imx) = f(ik+imx) - val*d(jl+imx)*fexch
                  f(il+imx) = f(il+imx) - val*d(jk+imx)*fexch
                  f(jk+imx) = f(jk+imx) - val*d(il+imx)*fexch
                  f(jl+imx) = f(jl+imx) - val*d(ik+imx)*fexch
                  imx = imx + l2
c           6-----
  300           continue
c         5-----
              endif
c
 1330         if ((itwoeps .eq. 1) .and. (ieps .eq. 0)
     1            .and. (fexch .ne. zero)) then
                ieps = 1
                goto 1300
              else
                ieps = 0
              endif
c      4-----  next l
 1400       continue
c     3-----  next k
 1500     continue
c   2-----  next j
 1600   continue
c 1-----  next i
 1700 continue
c
c-----  save incomplete fock matrix
c
      call dawrit(idafdrf,iodadrf,f,nset*l2,72,navdrf)
c
c-----  correct off-diagonal elements
c
      do 2000, i = 1, num
        do 2000, j = 1, i-1
          ij = ia(i)+j
          imx = 0
          do 1900, iset = 1, nset
            f(ij+imx) = f(ij+imx)*pt5
            imx = imx + l2
 1900     continue
 2000 continue
c
      return
      end
