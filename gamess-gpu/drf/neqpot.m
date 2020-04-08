      subroutine neqread
c------
c      reads texts defining non-equilibrium reaction fields
c      from input ($neq) and connects to address on da-files
c      containing induced moments or expanded potentials
c      and polarisation energies
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/mollab)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/neqtex)
INCLUDE(comdrf/inddaf)
c
      character*80 blank
      character*80 text
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-neqread'/
      data blank /' '/
c
c-----  begin
c
c-----  get names associated with available induced moments into textneq
c
      if (iodaind(1,1) .eq. 0) then
        do 100 i = 1, maxneq
          textneq(i) = blank
 100    continue
        call neqout(textneq,realneq,maxneq)
        call dawrit(idafind,iodaind,realneq,maxneq*10,1,navind)
      endif
c
      call daread(idafind,iodaind,realneq,maxneq*10,1)
      call neqin(textneq,realneq,maxneq)
c
c
      if ((neqsta .eq. 1) .or. (neqrf .eq. 1)) then
c
c  -----  read $neq for number and names of non-equilibrium rfs
c
            nneqrf = nneqrf + 1
c
c      -----  check if required induced moments or expanded potentials
c             are available on
c             da-file and get index number
c
            call getindx(textrf,textneq,maxneq,ineqix(nneqrf))
c
      endif
c
c-----  check whether there is already a set of induced moments
c       on file with the current title
c
      if (isolsav .eq. 1) then
        write(text,1300) (title(i), i = 1, 10)
 1300   format(10a8)
c
        do 200, i = 1, maxneq
          if (text .eq. textneq(i)) then
            write(iwr,1400) text
 1400       format(/,' set of induced moments by name of:',/,a80,/,
     1  ' is already defined on tape44, choose other title in $basis')
c           call hnderr(3,errmsg)
c           stop
          endif
  200   continue
      endif
c
      return
      end
      subroutine getindx(string,dastr,ntex,index)
c------
c      looks up index of array by name 'string' on da-file,
c      by comparing with text written on da-file
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      character*80 string
      character*80 dastr(ntex)
c
INCLUDE(comdrf/iofil)
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-getindx'/
c
      do 100, i = 1, ntex
        if (string .eq. dastr(i)) goto 200
  100 continue
      write(iwr,1001) string
 1001 format(/,' not found on da-file: string by name:',/,2x,a80)
      call hnderr(3,errmsg)
      stop
  200 continue
c     index = i + 1
      index = 2
c
      return
      end
      subroutine neqsav(idrfout,ndim,ifldin,ngran,idisadd,isolsav,ieps,
     2      imomsav,ineqex,maxneq,nwtr,nwtc,ilvr,nchd,nexpc,neqgrp,
     3      dipind,sdips,vr,ovl,rx,ry,rz,iexpc,vrf,hrf,hrfh)
c------
c      calculates and saves induced moments
c      (or goes on to potential and fields at expansion centra)
c
c      these may be used later to define nonequilibrium rfs
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/mollab)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
INCLUDE(comdrf/neqpar)
c
INCLUDE(comdrf/runpar)
c
      dimension neqgrp(ngran)
      dimension dipind(ndim), sdips((ngran+2)*ndim)
      dimension vr(nwtr,nwtc)
      dimension ovl(nchd),rx(nchd),ry(nchd),rz(nchd)
      dimension hrf(nchd), hrfh(nchd), vrf(nwtc)
      dimension iexpc(nchd)
c
      character*80 blank
c
      character *8 errmsg(3)
c
      data errmsg /'program','stop in','-neqsav-'/
c
      data blank /' '/
c
c-----  begin
c
      if ((isolsav .eq. 1) .and. (.not. mcupdt)) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
c    -----  find the last record containing induced moments
c           to define the index for the present one
c
          do 100, i = 1, maxneq
            if (textneq(i) .eq. blank) goto 200
  100     continue
c
          write(iwr,10001) maxneq
10001     format(/,' already ',i4,'sets of non-equilibrium reaction',
     1             ' fields on neqind etc: clean up')
          call hnderr(3,errmsg)
  200     continue
c
c         index = i + 1
          index = 2
c
          write(textneq(index-1),10002) (title(i), i = 1, 10)
10002     format(10a8)
c
c    -----  write updated register on record -1- of idafind
c
          call neqout(textneq,realneq,maxneq)
          call dawrit(idafind,iodaind,realneq,maxneq*10,1,navind)
c   2-----
        endif
c 1-----
      endif
c
c-----  add nuclear and electronic and external contributions
c
      do 1100, i = 1, ndim
        dipind(i) = sdips(i)
 1100 continue
c
      if (ifldin .eq. 4) then
c 1-----
        do 1200, i = 1, ndim
          dipind(i) = dipind(i) + sdips(ndim+i)
 1200   continue
c
        do 1300, igr = 1, ngran
c   2-----
          if (neqgrp(igr) .eq. 1) then
c     3-----
            do 1250, i = 1, ndim
              dipind(i) = dipind(i) + sdips((igr+1)*ndim+i)
 1250       continue
c     3-----
          endif
c   2-----
 1300   continue
c 1-----
      else
c 1-----
        if (idisadd .eq. 0) then
c   2-----
          do 1400, igr = 1, ngran
c     3-----
            if (neqgrp(igr) .eq. 1) then
c       4-----
              do 1350, i = 1, ndim
                dipind(i) = dipind(i) + sdips(igr*ndim+i)
 1350         continue
c       4-----
            endif
c     3-----
 1400     continue
c   2-----
        endif
c 1-----
      endif
c
      if ((isolsav .eq. 1) .and. (imomsav .eq. 1)) then
c 1-----
c  -----  save induced moments
c
        if (ieps .eq. 0) then
c   2-----
c    -----  write induced moments on record -index- of idafind
c                 (total dielectric contribution)
c
          call dawrit(idafind,iodaind,dipind,ndim,index,navind)
        else
c
c    -----  write induced moments on record -index- of idafino
c                 (optic dielectric contribution)
c
          call dawrit(idafino,iodaino,dipind,ndim,index,navino)
c   2-----
        endif
c 1-----
      endif
c
      if ((ineqex .eq. 1) .or. (imomsav .ne. 1)) then
c 1-----
c  -----  read vr-matrix
c
        call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvr)
c
        if (ineqex .eq. 1) then
c   2-----
c    -----  read overlap and first moment integrals
c
          call daread(idafh,ioda,ovl,nchd,12)
          call daread(idafh,ioda,rx,nchd,53)
          call daread(idafh,ioda,ry,nchd,54)
          call daread(idafh,ioda,rz,nchd,55)
c
c    -----  read assignment vector
c
c         nchd2 = lenint(nchd)
c         call daread(idafdrf,iodadrf,iexpc,nchd2,2)
          call daread(idafdrf,iodadrf,iexpc,nchd,2)
c
c    -----  calculate reaction potential in ao basis
c
          call onerf(vr,ovl,rx,ry,rz,hrf,iexpc,dipind)
c
          if (ieps .eq. 1) then
c     3-----
c      -----  subtract optical component from total component
c
            call daread(idafdrf,iodadrf,hrfh,nchd,91)
c
            do 1700, i = 1, nchd
              hrf(i) = hrfh(i) - hrf(i)
 1700       continue
c     3-----
          endif
c
c    -----  save potential in ao basis on drfdaf, record -91-
c
          call dawrit(idafdrf,iodadrf,hrf,nchd,91,navdrf)
c   2-----
        endif
c
        if (imomsav .ne. 1) then
c   2-----
c    -----  calculate potential and field at expansion centra
c
          call rfpot(idrfout,nwtr,nwtc,ndim,nexpc,ngran,nqdim,
     2               neqgrp,vr,vrf,dipind)
c
c    -----  save reaction potential & field at expansion centra
c
          if (ieps .eq. 0) then
c     3-----
c      -----  write on record -index- of idafind
c                 (total dielectric contribution)
c
            call dawrit(idafind,iodaind,vrf,nqdim,index,navind)
          else
c
c      -----  write on record -index- of idafino
c                 (optic dielectric contribution)
c
            call dawrit(idafino,iodaino,vrf,nqdim,index,navino)
c     3-----
          endif
c   2-----
        endif
c 1-----
      endif
c
      return
      end
      subroutine nqsavs(nwtc,ngran,nzfa,nexpc,nqdim,nqcls,neqgrp,
     2                  uclaseg,uclasdg,uclasrg,extnucg,repmodg,pot)
c------
c      saves total external static potential at expansion centra and
c      classical groups that do not contribute to the saved potential
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension pot(nwtc)
      dimension uclaseg(ngran*(ngran+1)/2),uclasdg(ngran*(ngran+1)/2),
     2          uclasrg(ngran*(ngran+1)/2),
     3          extnucg(ngran),repmodg(ngran)
      dimension neqgrp(ngran)
c
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/auxdrf)
c
cmw      include 'drf/dimpar'
INCLUDE(comdrf/drfzfa)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c
c
      data zero /0.0d00/
c
c-----  begin
c
      call clear (pot,nwtc)
      call clear (auxx,3*ngran)
c
      icnt = 0
c
c-----  add external static potential at expansion centra
c
c-----  collect ad hoc contributions (repulsion, dispersion) as well
c
      do 1100, igr = 1, ngran
c 1-----
        if (neqgrp(igr) .eq. 1) then
c   2-----
c    -----  expanded potential at the expansion centra
c
cnot          call addup(pot,zfa(1,igr),pot,4*nexpc)
c
          icnt = 0
c
          do 1010, kgr = 1, ngran
c     3-----
            if (neqgrp(kgr) .ne. 1) then
c       4-----
c        -----  electrostatic interaction energy between classical group
c
              icnt = icnt + 1
              indxcl = ia(max(igr,kgr)) + min(igr,kgr)
              pot(4*nexpc+icnt) = pot(4*nexpc+icnt) + uclaseg(indxcl)
c
c        -----  dispersion and repulsion
c
              auxx(icnt) = auxx(icnt) + uclasdg(indxcl)
              auxx(icnt+ngran) = auxx(icnt+ngran) + uclasrg(indxcl)
c       4-----
            endif
c     3-----
 1010     continue
c
          if (icnt .eq. 0) then
c     3-----
            icnt = icnt + 1
            pot(4*nexpc+icnt) = zero
            auxx(icnt) = zero
            auxx(icnt + ngran) = zero
c     3-----
          endif
c
c    -----  interaction with nuclei
c
          icnt = icnt + 1
c
c    -----  electrostatic
c
          pot(4*nexpc+icnt) = pot(4*nexpc+icnt) + extnucg(igr)
c
c    -----  repulsion
c
          auxx(icnt+2*ngran) = auxx(icnt+2*ngran) + repmodg(igr)
c   2-----
        endif
c 1-----
 1100 continue
c
      nqdim = 4*nexpc + icnt
      nqcls = icnt - 1
c
c-----  save collected total external potential
c
      call dawrit(idafsta,iodasta,pot(1),nqdim,index,navsta)
c
c-----  save collected dispersion interaction
c
      call dawrit(idafdis,iodadis,auxx,icnt-1,index,navdis)
c
c-----  save collected repulsion interactions
c
      call dawrit(idafrep,iodarep,auxx(ngran+1),icnt-1,index,navrep)
c
      call dawrit(idafrqm,iodarqm,auxx(2*ngran+1),1,index,navrqm)
c
      auxx(1) = dble(nqdim)
      call dawrit(idafinf,iodainf,auxx(1),1,index,navinf)
c
      return
      end
      subroutine nqcost(field,mxgran,ngran,ieps,idisadd,imomsav,
     2           itwoeps,eps1,eps2,kappa1,kappa2,neqgrp,
     3           uclaseg,uclasdg,uclasrg,suclasg,smolxtg,sxtmolg,
     4           snucnuc,stwoel,snucel,selnuc,
     5           upolclg,ucstneq,upoleqs,stabtot)
c------
c      administrates static and reaction field costs and options
c      and saves relevant numbers on da-file
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesahnqcost)
c
      character*8 field
c
      REAL kappa1, kappa2
      dimension neqgrp(ngran)
c
      dimension uclaseg(ngran*(ngran+1)/2),uclasdg(ngran*(ngran+1)/2),
     2 uclasrg(ngran*(ngran+1)/2),smolxtg(ngran),sxtmolg(ngran),
     3 upolclg(ngran)
c
      dimension suclasg(mxgran,ngran)
c
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/auxdrf)
c
INCLUDE(comdrf/neqtex)
INCLUDE(comdrf/inddaf)
c
      data zero, pt5 /0.0d00,5.0d-01/
c
c-----  begin
c
      if (ieps .eq. 0) then
c 1-----
c  -----  static and ad-hoc contributions
c
        elst = zero
        disp = zero
        rep = zero
c
        do 100, igr = 1, ngran
c   2-----
          if (neqgrp(igr) .eq. 1) then
c     3-----
c      -----  self-energy of classical analysis groups:
c
            do 10, jgr = 1, igr
c       4-----
              if (neqgrp(jgr) .eq. 1) then
c         5-----
                indxcl = ia(igr) + jgr
                elst = elst + uclaseg(indxcl)
                disp = disp + uclasdg(indxcl)
                rep = rep + uclasrg(indxcl)
c         5-----
              endif
c       4-----
   10       continue
c     3-----
          endif
c   2-----
  100   continue
c
        ucstneq = elst + disp + rep
c
        auxx(1) = elst
        auxx(2) = disp
        auxx(3) = rep
c
        call dawrit(idafcst,iodacst,auxx,3,index,navcst)
c 1-----
      endif
c
      if (field(5:) .ne. ' ') then
c 1-----
c  -----  reaction field contributions
c
        ucst = zero
c
c  -----  now split contributions
c
        do 1100, igr = 1, ngran
c   2-----
          if (neqgrp(igr) .eq. 1) then
c     3-----
c      -----  screening from other groups defining neqrf
c
            do 1010, jgr = 1, ngran
c       4-----
              if (neqgrp(jgr) .eq. 1)
     2          ucst = ucst - pt5*suclasg(igr,jgr)
c       4-----
 1010       continue
c
c      -----  screening with qm
c
            ucst = ucst - pt5*smolxtg(igr)
            ucst = ucst + pt5*sxtmolg(igr)
c     3-----
          endif
c   2-----
 1100   continue
c
c    -----  quantum mechanical system polarisation energy
c
        ucst = ucst + pt5*snucnuc
        ucst = ucst + pt5*stwoel
        ucst = ucst + pt5*snucel
        ucst = ucst + pt5*selnuc
c
        if (imomsav .eq. 0) then
          polcst = -ucst
        else
          polcst = -pt5*stabtot
        endif
c
        if (ieps .eq. 1) then
c   2-----
c    -----  read equilibrium polarisation energy for total dielectric
c           constant for the whole system
c
          call daread(idafpol,iodapol,auxx,7,index)
c
c    -----  subtract equilibrium polarisation energy for optic dielectri
c           constant for the whole system
c
          upoleqs = auxx(1) - polcst
        else
          upoleqs = polcst
c   2-----
        endif
c
c  -----  write equilibrium polarisation energy and options
c
        auxx(1) = upoleqs
        auxx(2) = eps1
        auxx(3) = kappa1
        auxx(4) = eps2
        auxx(5) = kappa2
        auxx(6) = dble(itwoeps)
        auxx(7) = dble(imomsav)
c
        call dawrit(idafpol,iodapol,auxx,7,index,navpol)
c 1-----
      endif
c
      return
      end
      subroutine neqon2(istat,ieps,vrpot,s,dx,dy,dz,
     2                  vneq,vneqh,iexpc)
c------
c      calculates expanded potential per overlap distribution
c      due to nonequilibrium potential
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/rfene)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
INCLUDE(comdrf/neqpar)
c
INCLUDE(comdrf/runpar)
c
      dimension vrpot(neqdim)
      dimension s(nchd), dx(nchd), dy(nchd), dz(nchd)
      dimension iexpc(nchd)
      dimension vneq(nchd), vneqh(nchd)
c
      dimension f(4)
c
      data zero, pt5 /0.0d00,0.5d00/
c
c-----  begin
c
c-----  read overlap + 1st moment integrals, and assignment vector
c
      call daread(idafh,ioda,s ,nchd,12)
      call daread(idafh,ioda,dx,nchd,53)
      call daread(idafh,ioda,dy,nchd,54)
      call daread(idafh,ioda,dz,nchd,55)
c     nchd2 = lenint(nchd)
c     call daread(idafdrf,iodadrf,iexpc,nchd2,2)
      call daread(idafdrf,iodadrf,iexpc,nchd,2)
c
      if (istat .eq. 1) then
c 1-----
c  -----  external static potential
c
        call daread(idafsta,iodasta,vrpot,neqdim,ineqix(nneq))
c 1-----
      else
c 1-----
c  -----  reaction field potential
c  -----  set record & files and read
c
        if (ieps .eq. 0) then
c   2-----
          call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
c   2-----
        else
c   2-----
          call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
c   2-----
        endif
c 1-----
      endif
c
      call clear(vneq,nchd)
c
c-----  loop over expansion centra
c
      ij = 0
      do 100, i = 1, num
c 1-----
        do 200, j = 1, i
c   2-----
          ij = ij + 1
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
          f(1) = dx(ij)
          f(2) = dy(ij)
          f(3) = dz(ij)
          f(4) = s(ij)
c
c    -----  execute expansion
c
          vneqd = zero
          do 300, k = 1, 4
c     3-----
            vneqd = vneqd + f(k)*vrpot(ip+k)
c     3-----
  300     continue
c
c    -----  note: minus sign because this potential is always coupled ba
c           to electrons
c
          vneq(ij) = vneq(ij) - vneqd
c   2-----
  200   continue
c 1-----
  100 continue
c
      if (idrfout .eq. 3) call hatout(vneq,num,num,3,'vneq')
c
      if (istat .eq. 1) then
c 1-----
c  -----  write electrostatic (nonequilibrium) external potential
c
        call dawrit(idafdrf,iodadrf,vneq,nchd,181,navdrf)
        call neqnuc(nat,neqdim,czan,c,vrpot,unqnuc)
        ustanuc(iactst) = unqnuc
c 1-----
      else
c 1-----
        if (ieps .eq. 1) then
c   2-----
c    -----  read vneqh with total dielectric and subtract the present on
c           (optic dielectric)
c
          call daread(idafdrf,iodadrf,vneqh,nchd,14)
          do 500, i = 1, nchd
c     3-----
            vneq(i) = vneqh(i) - vneq(i)
c     3-----
  500     continue
c
          uneqnt = uneqnuc(iactst)
c   2-----
        endif
c
c  -----  calculate interaction of nuclei with non-equilibrium rf
c
        call neqnuc(nat,nwtc,czan,c,vrpot,unqnuc)
        uneqnuc(iactst) = unqnuc
c
        if (ieps .eq. 1) uneqnuc(iactst) = uneqnt - uneqnuc(iactst)
c
c  -----  write non-equilibrium potential on da31, record -14-
c
        call dawrit(idafdrf,iodadrf,vneq,nchd,14,navdrf)
c 1-----
      endif
c
      return
      end
      subroutine neqone(ieps,momsav,vr,s,dx,dy,dz,vneq,vneqh,
     2                  iexpc,dipind,vrpot)
c------
c      calculates expanded potential due to nonequilibrium rf
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/rfene)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
INCLUDE(comdrf/runpar)
c
      dimension vr(nwtr,nwtc)
      dimension s(nchd), dx(nchd), dy(nchd), dz(nchd)
      dimension iexpc(nchd)
      dimension vneq(nchd), vneqh(nchd), dipind(ndim)
      dimension vrpot(neqdim)
c
      dimension f(4)
      data zero, pt5 /0.0d00,0.5d00/
c
c-----  begin
c
c-----  read overlap + 1st moment integrals, and assignment vector
c
      call daread(idafh,ioda,s ,nchd,12)
      call daread(idafh,ioda,dx,nchd,53)
      call daread(idafh,ioda,dy,nchd,54)
      call daread(idafh,ioda,dz,nchd,55)
c     nchd2 = lenint(nchd)
c     call daread(idafdrf,iodadrf,iexpc,nchd2,2)
      call daread(idafdrf,iodadrf,iexpc,nchd,2)
c
c-----  set record & files and read
c
      if (momsav .eq. 1) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
****
          call daread(idafind,iodaind,dipind,ndim,ineqix(nneq))
          ilvr = 34
c   2-----
        else
c   2-----
          call daread(idafino,iodaino,dipind,ndim,ineqix(nneq))
          ilvr = 36
c   2-----
        endif
c
        if (mcupdt) ilvr = ilvr + 1
c
c  -----  read vr matrix
c
        call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilvr)
c 1-----
      else
c 1-----
        if (ieps .eq. 0) then
c   2-----
          call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
c   2-----
        else
c   2-----
          call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
c   2-----
        endif
c 1-----
      endif
c
      call clear(vneq,nchd)
c
c-----  loop over expansion centra
c
      ij = 0
      do 100, i = 1, num
c 1-----
        do 200, j = 1, i
c   2-----
          ij = ij + 1
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
          f(1) = dx(ij)
          f(2) = dy(ij)
          f(3) = dz(ij)
          f(4) = s(ij)
c
c    -----  execute expansion
c
          vneqd = zero
          do 300, k = 1, 4
c     3-----
            if (momsav .eq. 1) then
c       4-----
              do 400, l = 1, ndim
c         5-----
                vneqd = vneqd + f(k)*vr(l,ip+k)*dipind(l)
c         5-----
  400         continue
c       4-----
            else
c       4-----
              vneqd = vneqd + f(k)*vrpot(ip+k)
c       4-----
            endif
c     3-----
  300     continue
c
c    -----  note: minus sign because this potential is always coupled ba
c           to electrons
c
          vneq(ij) = vneq(ij) - vneqd
c   2-----
  200   continue
c 1-----
  100 continue
c
      if (idrfout .eq. 3) call hatout(vneq,num,num,3,'vneq')
c
      if (ieps .eq. 1) then
c 1-----
c  -----  read vneqh with total dielectric and subtract the present one
c         (optic dielectric)
c
        call daread(idafdrf,iodadrf,vneqh,nchd,14)
        do 500, i = 1, nchd
c   2-----
          vneq(i) = vneqh(i) - vneq(i)
c   2-----
  500   continue
c
        uneqnt = uneqnuc(iactst)
c 1-----
      endif
c
c-----  calculate interaction of nuclei with non-equilibrium rf
c
      if (momsav .eq. 1) then
c       uneqnuc(iactst) = adotb(vr(1,nexp4+ngran+1),dipind,ndim)
        uneqnuc(iactst) = ddot(ndim,vr(1,nexp4+ngran+1),1,dipind,1)
      else
        call neqnuc(nat,neqdim,czan,c,vrpot,uneqnuc(iactst))
      endif
c
      if (ieps .eq. 1) uneqnuc(iactst) = uneqnt - uneqnuc(iactst)
c
c-----  write non-equilibrium potential on da31, record -14-
c
      call dawrit(idafdrf,iodadrf,vneq,nchd,14,navdrf)
c
      return
      end
      subroutine neqnuc(nat,neqdim,czan,c,vrpot,uneqnuc)
c------
c      calculates interaction of nuclei with nonequilibrium
c      potential from saved expansion
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension vrpot(neqdim)
      dimension czan(nat)
      dimension c(3,nat)
c
      dimension q(4)
c
      data zero, one /0.0d00, 1.0d00/
c
c-----  begin
c
      uneqnuc = zero
c
      do 100, n = 1, nat
c 1-----
        q(1) = c(1,n)
        q(2) = c(2,n)
        q(3) = c(3,n)
        q(4) = one
c
        ip = (n-1)*4 + 1
c
c       uneqnuc = uneqnuc + adotb(vrpot(ip),q,4)*czan(n)
        uneqnuc = uneqnuc + ddot(4,vrpot(ip),1,q,1)*czan(n)
c 1-----
  100 continue
c
      return
      end
      subroutine rfpot(idrfout,nwtr,nwtc,ndim,nexpc,ngran,neqdim,
     2                 neqgrp,vr,vrpot,dipind)
c------
c      calculates expanded reaction potential due to induced moments
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
c
      dimension neqgrp(ngran)
c
      dimension vr(nwtr,nwtc)
      dimension vrpot(nwtc), dipind(ndim)
c
      data zero, pt5 /0.0d00,0.5d00/
c
c-----  begin
c
      call clear(vrpot,nwtc)
c
c-----  loop over expansion centra for electrons
c
      do 100, i = 1, nexpc
c 1-----
        ip = (i-1)*4
c
c  -----  execute expansion
c
        do 200, k = 1, 4
c   2-----
          do 300, l = 1, ndim
c     3-----
            vrpot(ip+k) = vrpot(ip+k) + vr(l,ip+k)*dipind(l)
c     3-----
  300     continue
c   2-----
  200   continue
c 1-----
  100 continue
c
c-----  loop over external charge groups
c
      icnt = 0
c
      do 1100, i = 1, ngran
c 1-----
        if (neqgrp(i) .ne. 1) then
c   2-----
          icnt = icnt + 1
          do 1200, l = 1, ndim
            vrpot(nexpc*4+icnt) = vrpot(nexpc*4+icnt)
     2                         + vr(l,(nexpc*4+i))*dipind(l)
 1200     continue
c   2-----
        endif
c 1-----
 1100 continue
c
      if (icnt .eq. 0) then
c 1-----
        icnt = icnt + 1
        vrpot(nexpc*4+icnt) = zero
c 1-----
      endif
c
c-----  nuclei
c
      icnt = icnt + 1
      do 2100, l = 1, ndim
            vrpot(nexpc*4+icnt) = vrpot(nexpc*4+icnt)
     2                         + vr(l,(nexpc*4+ngran+1))*dipind(l)
 2100 continue
c
      neqdim = 4*nexpc + icnt
c
      if (idrfout .eq. 3) call hatout(vrf,neqdim,1,2,'vrpot')
c
      return
      end
      subroutine onerf2(vrf,vrpot,s,dx,dy,dz,iexpc)
c------
c      calculates expanded potential due to non-equilibrium rf
c      for overlap distributions
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
cmw      include 'drf/dimpar'
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/rfene)
c
INCLUDE(comdrf/runpar)
c
      dimension s(nchd), dx(nchd), dy(nchd), dz(nchd)
      dimension iexpc(nchd)
      dimension vrf(nchd), vrpot(neqdim)
c
      dimension f(4)
c
      data zero, pt5 /0.0d00,0.5d00/
c
c-----  begin
c
      call clear(vrf,nchd)
c
c-----  loop over expansion centra
c
      ij = 0
      do 100, i = 1, num
c 1-----
        do 200, j = 1, i
c   2-----
          ij = ij + 1
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
          f(1) = dx(ij)
          f(2) = dy(ij)
          f(3) = dz(ij)
          f(4) = s(ij)
c
c    -----  execute expansion
c
          vrfd = zero
          do 300, k = 1, 4
c     3-----
            vrfd = vrfd + f(k)*vrpot(ip+k)
c     3-----
  300     continue
c
c    -----  note: minus sign because this potential is always coupled ba
c           to electrons
c
          vrf(ij) = vrf(ij) - vrfd
c   2-----
  200   continue
c 1-----
  100 continue
c
      if (idrfout .eq. 3) call hatout(vrf,num,num,3,'vrf')
c
      return
      end
      subroutine neqdip(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           neqdim,neq2eps,itsolv,field,ilvr,vr,vrpot,dip,
     3           dipind,dipino,uneqqm,uneqqmo)
c------
c      drives the calculation of the interaction of
c      dipoles at the expansion centra with the nonequilibrium
c      reaction field
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      logical mcupdt
c
      character*8 field
c
      dimension vr(nwtr,nwtc)
      dimension vrpot(neqdim)
      dimension dip(3*nexp), dipind(npol3+ndimb), dipino(npol3+ndimb)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c
c-----  begin
c
      if (momsav .eq. 1) then
c 1-----
c  -----  the neqrf has been saved in terms of the induced moments
c
        ilvrneq = 34
        if (mcupdt) ilvrneq = ilvrneq + 1
        call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
        call diprf(npol3,ndimb,nexp,vr,dip,dipind,uneqqm)
c
        if (neq2eps .eq. 1) then
c   2-----
          ilvrneq = 36
          if (mcupdt) ilvrneq = ilvrneq + 1
          call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
          call diprf(npol3,ndimb,nexp,vr,dip,dipino,uneqqmo)
          uneqqm = uneqqm - uneqqmo
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  the neqrf has been saved as potential at expansion centra
c
        call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
        call dipneq(nexp,vrpot,dip,uneqqm)
c
        if (neq2eps .eq. 1) then
c   2-----
          call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
          call dipneq(nexp,vrpot,dip,uneqqmo)
          uneqqm = uneqqm - uneqqmo
c   2-----
        endif
c 1-----
      endif
c
      return
      end
      subroutine dipneq(nexp,vrpot,dip,ediprf)
c------
c      this routine calculates the interaction of a
c      nonequilibrium potential with
c      dipoles at the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vrpot(nexp*4),dip(3,nexp)
c
      data zero /0.0d00/
c
c-----  initialize energy contribution
c
      ediprf = zero
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate pointer
c
        ip = (i-1)*4
c
c  -----  loop over x,y,z components
c
        do 200, k = 1, 3
c
c    -----  calculate interaction energy
c
          ediprf = ediprf + dip(k,i)*vrpot(ip+k)
  200   continue
c
  100 continue
      return
      end
      subroutine neqq(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           neqdim,neq2eps,itsolv,field,ilvr,vr,vrpot,q,
     3           dipind,dipino,uneqqm,uneqqmo)
c------
c      drives the calculation of the interaction of
c      charges at the expansion centra with the nonequilibrium
c      reaction field
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      logical mcupdt
c
      character*8 field
c
      dimension vr(nwtr,nwtc)
      dimension vrpot(neqdim)
      dimension q(nexp), dipind(npol3+ndimb), dipino(npol3+ndimb)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c-----  begin
c
      if (momsav .eq. 1) then
c 1-----
c  -----  the neqrf has been saved in terms of the induced moments
c
        ilvrneq = 34
        if (mcupdt) ilvrneq = ilvrneq + 1
        call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
        call qrf(npol3,ndimb,nexp,vr,q,dipind,uneqqm)
c
        if (neq2eps .eq. 1) then
c   2-----
          ilvrneq = 36
          if (mcupdt) ilvrneq = ilvrneq + 1
          call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
          call qrf(npol3,ndimb,nexp,vr,q,dipino,uneqqmo)
          uneqqm = uneqqm - uneqqmo
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  the neqrf has been saved as potential at expansion centra
c
        call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
        call qneq(nexp,vrpot,q,uneqqm)
c
        if (neq2eps .eq. 1) then
c   2-----
          call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
          call qneq(nexp,vrpot,q,uneqqmo)
          uneqqm = uneqqm - uneqqmo
c   2-----
        endif
c 1-----
      endif
c
      return
      end
      subroutine qneq(nexp,vrpot,q,eqrf)
c------
c      this routine calculates the interaction of a
c      nonequilibrium potential with
c      charges at the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vrpot(nexp*4+2)
      dimension q(nexp)
c
      data zero /0.0d00/
c
c-----  initialize energy contrbution
c
      eqrf = zero
c
c-----  loop over expansion centra
c
      do 100, i = 1, nexp
c
c  -----  calculate interaction energy of charges with neq rf
c
        eqrf = eqrf + q(i)*vrpot((i-1)*4+4)
  100 continue
      return
      end
      subroutine neqel(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           ngran,neqdim,neq2eps,itsolv,field,ilvr,
     3           vr,vrpot,dens,ovl,rx,ry,rz,iexpc,dipind,
     5           dipino,uneqel,uneqnuc)
c------
c      drives the calculation of the interaction of
c      electrons and nuclei at the expansion centra with the
c      nonequilibrium reaction field
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      logical mcupdt
c
      character*8 field
c
      dimension vr(nwtr,nwtc)
      dimension vrpot(neqdim)
      dimension dipind(npol3+ndimb), dipino(npol3+ndimb)
c
      dimension dens(num*(num+1)/2), ovl(num*(num+1)/2),
     2  rx(num*(num+1)/2), ry(num*(num+1)/2), rz(num*(num+1)/2)
      dimension iexpc(num*(num+1)/2)
c
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c
c-----  begin
c
c
c-----  calculate interaction of nuclei with non-equilibrium rf
c
      if (momsav .eq. 1) then
c 1-----
c  -----  the neqrf has been saved in terms of the induced moments
c
        ilvrneq = 34
        if (mcupdt) ilvrneq = ilvrneq + 1
        call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
        call elrf(npol3,ndimb,nexp,num,vr,dens,
     2           rx,ry,rz,ovl,iexpc,dipind,uneqel)
c
c       uneqnuc = adotb(vr(1,nexp*4+ngran+1),dipind,npol3+ndimb)
        uneqnuc = ddot(npol3+ndimb,vr(1,nexp*4+ngran+1),1,dipind,1)
c
        if (neq2eps .eq. 1) then
c   2-----
          ilvrneq = 36
          if (mcupdt) ilvrneq = ilvrneq + 1
          call daread(idafdrf,iodadrf,vr,nwtr*nwtc,ilvrneq)
          call elrf(npol3,ndimb,nexp,num,vr,dens,
     2           rx,ry,rz,ovl,iexpc,dipino,uneqelo)
          uneqel = uneqel - uneqelo
c
          uneqnuc = uneqnuc -
     2    ddot(npol3+ndimb,vr(1,nexp*4+ngran+1),1,dipino,1)
c    2    adotb(vr(1,nexp*4+ngran+1),dipino,npol3+ndimb)
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  the neqrf has been saved as potential at expansion centra
c
        call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
        call elneq(nexp,num,vrpot,dens,rx,ry,rz,ovl,
     2                 iexpc,uneqel)
c
        call neqnuc(nat,neqdim,czan,c,vrpot,uneqnuc)
c
        if (neq2eps .eq. 1) then
c   2-----
          call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
          call elneq(nexp,num,vrpot,dens,rx,ry,rz,ovl,
     2                 iexpc,uneqelo)
          uneqel = uneqel - uneqelo
c
          call neqnuc(nat,neqdim,czan,c,vrpot,uneqno)
          uneqnuc = uneqnuc - uneqno
c   2-----
        endif
c 1-----
      endif
c
      return
      end
      subroutine elneq(nexp,num,vrpot,d,rx,ry,rz,s,
     2                 iexpc,eelmol)
c------
c      this routine calculates the interaction energy between
c      the expanded neqrf and the electronic
c      density
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vrpot(nexp*4)
      dimension d(num*(num+1)/2), rx(num*(num+1)/2),ry(num*(num+1)/2),
     1          rz(num*(num+1)/2),s(num*(num+1)/2)
      dimension iexpc(num*(num+1)/2)
c
c----- local arrays
c
      dimension f(4)
c
      data zero, one, two /0.0d00, 1.0d00, 2.0d00/
c
c-----  initialize interaction energy
c
      eelmol = zero
c
c-----  effective loop over charge distributions
c
      ij = 0
      do 100, i = 1, num
        do 200, j = 1, i
          ij = ij + 1
c
          fact = two
          if (j .eq. i) fact = one
c
c    -----  assign expansion centre to charge distribution
c
          iexp = iexpc(ij)
c
c    -----  set pointer
c
          ip = (iexp-1)*4
c
c    -----  premultiply density with dipole and overlap moments
c           of the charge distribution
c
          f(1) = fact*d(ij)*rx(ij)
          f(2) = fact*d(ij)*ry(ij)
          f(3) = fact*d(ij)*rz(ij)
          f(4) = fact*d(ij)*s(ij)
c
c    -----  execute expansion
c
          do 300, k = 1, 4
            eelmol = eelmol - f(k)*vrpot(ip+k)
  300     continue
c
  200   continue
  100 continue
c
      return
      end
      subroutine neqcl(istat,momsav,ngran,ndim,nwtc,nwtr,nexp,nneq,
     2                 neqdim,mcupdt,neq2eps,vr,vrpot,dipind,dipino,
     3                 uneqclg,uneqcl)
c------
c      this routine calculates the interaction energy between
c      the neqrf and the external charges
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      logical mcupdt
c
      dimension vr(nwtr,nwtc)
      dimension vrpot(neqdim), dipind(ndim), dipino(ndim)
      dimension uneqclg(ngran)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c-----  begin
c
      if (istat .eq. 0) then
c 1-----
c  -----  reaction field contributions
c
        if (momsav .eq. 1) then
c   2-----
          do 100, igr = 1, ngran
c     3-----
            ilvrneq = 34
            if (mcupdt) ilvrneq = ilvrneq + 1
            call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilvrneq)
c           uneqclg(igr) = adotb(vr(1,nexp*4+igr),dipind,ndim)
            uneqclg(igr) = ddot(ndim,vr(1,nexp*4+igr),1,dipind,1)
c
            if (neq2eps .eq. 1) then
c       4-----
              ilvrneq = 36
              if (mcupdt) ilvrneq = ilvrneq + 1
              call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilvrneq)
c             uneqcog = adotb(vr(1,nexp*4+igr),dipino,ndim)
              uneqcog = ddot(ndim,vr(1,nexp*4+igr),1,dipino,1)
              uneqclg(igr) = uneqclg(igr) - uneqcog
c       4-----
            endif
c
            uneqcl = uneqcl + uneqclg(igr)
c     3-----
 100      continue
c   2-----
        else
c   2-----
          do 200, igr = 1, ngran
c     3-----
c      -----  the neqrf has been saved directly
c
            call daread(idafind,iodaind,vrpot,neqdim,ineqix(nneq))
            uneqclg(igr) = vrpot(nexp*4+igr)
c
            if (neq2eps .eq. 1) then
c       4-----
              call daread(idafino,iodaino,vrpot,neqdim,ineqix(nneq))
              uneqcog = vrpot(nexp*4+igr)
              uneqclg(igr) = uneqclg(igr) - uneqcog
c       4-----
            endif
c
            uneqcl = uneqcl + uneqclg(igr)
c     3-----
  200     continue
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  interaction with external static field from previous calculati
c
        call daread(idafsta,iodasta,vrpot,neqdim,ineqix(nneq))
c
        do 1100, igr = 1, ngran
c   2-----
          uneqclg(igr) = vrpot(nexp*4+igr)
          uneqcl = uneqcl + uneqclg(igr)
c   2-----
 1100   continue
c 1-----
      endif
c
      return
      end
      subroutine addpot(nqdim,nqcls,itwoeps,field,pot,pot2,hlp)
c------
c      collects total potential at expansion centra
c      during a monte carlo run
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension pot(nqdim), pot2(nqdim), hlp(nqdim)
c
      character*8 field
c
cmw      include 'drf/dimpar'
INCLUDE(comdrf/drfzfa)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
c-----  begin
c
c-----  read collected total external potential at expansion centra
c
      call daread(idafdrf,iodadrf,pot(1),nqdim,151)
      call daread(idafdrf,iodadrf,pot2(1),nqdim,152)
c
c-----  read static external potential for accepted configuration
c
      call daread(idafsta,iodasta,hlp,nqdim,index+1)
c
c-----  add external static potential
c
      call addsp(hlp,nqdim,pot,pot2)
c
c-----  save collected total external potential
c
      call dawrit(idafdrf,iodadrf,pot(1),nqdim,151,navdrf)
      call dawrit(idafdrf,iodadrf,pot2(1),nqdim,152,navdrf)
c
c-----  read collected dispersion interaction
c
      call daread(idafdrf,iodadrf,pot(1),nqcls,161)
      call daread(idafdrf,iodadrf,pot2(1),nqcls,162)
c
      call daread(idafdis,iodadis,hlp,nqcls,index+1)
c
c-----  add dispersion interaction
c
      call addsp(hlp,nqcls,pot,pot2)
c
c-----  save collected dispersion interaction
c
      call dawrit(idafdrf,iodadrf,pot(1),nqcls,161,navdrf)
      call dawrit(idafdrf,iodadrf,pot2(1),nqcls,162,navdrf)
c
c-----  read collected repulsion interaction
c
      call daread(idafdrf,iodadrf,pot(1),nqcls,163)
      call daread(idafdrf,iodadrf,pot2(1),nqcls,164)
c
      call daread(idafrep,iodarep,hlp,nqcls,index+1)
c
c-----  add repulsion interaction
c
      call addsp(hlp,nqcls,pot,pot2)
c
c-----  save collected repulsion interaction
c
      call dawrit(idafdrf,iodadrf,pot(1),nqcls,163,navdrf)
      call dawrit(idafdrf,iodadrf,pot2(1),nqcls,164,navdrf)
c
c-----  read collected qm repulsion interaction
c
      call daread(idafdrf,iodadrf,pot(1),1,165)
      call daread(idafdrf,iodadrf,pot2(1),1,166)
c
      call daread(idafrqm,iodarqm,hlp,1,index+1)
c
c-----  add qm repulsion interaction
c
      call addsp(hlp,1,pot,pot2)
c
c-----  save collected qm repulsion interaction
c
      call dawrit(idafdrf,iodadrf,pot(1),1,165,navdrf)
      call dawrit(idafdrf,iodadrf,pot2(1),1,166,navdrf)
c
c-----  read actual total reaction potential if required
c
      if (field(5:) .ne. ' ') then
c 1-----
        call daread(idafdrf,iodadrf,pot(1),nqdim,153)
        call daread(idafdrf,iodadrf,pot2(1),nqdim,154)
c
c  -----  add total reaction potential
c
        call daread(idafind,iodaind,hlp,nqdim,index+1)
        call addsp(hlp,nqdim,pot,pot2)
c
c  -----  save collected total reaction potential
c
        call dawrit(idafdrf,iodadrf,pot(1),nqdim,153,navdrf)
        call dawrit(idafdrf,iodadrf,pot2(1),nqdim,154,navdrf)
c
        if (itwoeps .eq. 1) then
c   2-----
c    -----  read collected optic potential at expansion centra
c
          call daread(idafdrf,iodadrf,pot(1),nqdim,155)
          call daread(idafdrf,iodadrf,pot2(1),nqdim,156)
c
c    -----  read actual optic reaction potential
c
          call daread(idafino,iodaino,hlp,nqdim,index+1)
c
c    -----  add optic reaction potential
c
          call addsp(hlp,nqdim,pot,pot2)
c
c    -----  save collected optic reaction potential
c
          call dawrit(idafdrf,iodadrf,pot(1),nqdim,155,navdrf)
          call dawrit(idafdrf,iodadrf,pot2(1),nqdim,156,navdrf)
c   2-----
        endif
c
c  -----  polarisation energy cost
c
        call daread(idafdrf,iodadrf,pot(1),1,167)
        call daread(idafdrf,iodadrf,pot2(1),1,168)
        call daread(idafpol,iodapol,hlp,1,index+1)
c
        call addsp(hlp,1,pot,pot2)
c
        call dawrit(idafdrf,iodadrf,pot(1),1,167,navdrf)
        call dawrit(idafdrf,iodadrf,pot2(1),1,168,navdrf)
c 1-----
      endif
c
c-----  classial energy costs
c
      call daread(idafdrf,iodadrf,pot(1),3,169)
      call daread(idafdrf,iodadrf,pot2(1),3,170)
c
      call daread(idafcst,iodacst,hlp,3,index+1)
c
      call addsp(hlp,3,pot,pot2)
c
      call dawrit(idafdrf,iodadrf,pot(1),3,169,navdrf)
      call dawrit(idafdrf,iodadrf,pot2(1),3,170,navdrf)
c
      return
      end
      subroutine addsp(a,n,b,b2)
c------
c      adds a to b and a**2 to b2
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension a(n), b(n), b2(n)
c
      do 100, i = 1, n
c 1-----
        b(i) = b(i) + a(i)
        b2(i) = b2(i) + a(i)**2
c 1-----
  100 continue
c
      return
      end
      subroutine potsav(n,ncl,itwoeps,ntot,maxneq,field,pot,pot2,
     2                  poto,poto2,potf)
c------
c      averages potentials and saves them
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension pot(n),pot2(n),poto(n),poto2(n),potf(n)
c
      character*8 field
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
c
      dimension optneq(7)
c
      character*3 mcwrd1, mcwrd2
c
      data mcwrd1, mcwrd2
     2    /'mc ',  'mcf'/
c
      data thresh, zero /1.0d-20, 0.0d00/
c
c-----  begin
c
      textneq(index) = mcwrd1//textneq(index-1)(4:80)
      textneq(index+1) = mcwrd2//textneq(index-1)(4:80)
c
c-----  read collected total external potential and square
c
      call daread(idafdrf,iodadrf,pot,n,151)
      call daread(idafdrf,iodadrf,pot2,n,152)
c
c-----  average collected total potential and square
c
      do 100, i = 1, n
c 1-----
        pot(i) = pot(i)/(ntot+1)
        pot2(i) = pot2(i)/(ntot+1)
c
c  -----  calculate fluctuation potential
c
        fluct = pot2(i) - pot(i)**2
        if (fluct .lt. thresh) fluct = zero
        potf(i) = sqrt(fluct)
c 1-----
  100 continue
c
c-----  save collected total external potential
c       and fluctuation on da-file "neqsta"
c
      call dawrit(idafsta,iodasta,pot,n,index,navsta)
      call dawrit(idafsta,iodasta,potf,n,index+1,navsta)
c
c-----  dispersion interactions
c
      call daread(idafdrf,iodadrf,pot,ncl,161)
      call daread(idafdrf,iodadrf,pot2,ncl,162)
c
c-----  average collected total potential and square
c
      do 200, i = 1, ncl
c 1-----
        pot(i) = pot(i)/(ntot+1)
        pot2(i) = pot2(i)/(ntot+1)
c
c  -----  calculate fluctuation potential
c
        fluct = pot2(i) - pot(i)**2
        if (fluct .lt. thresh) fluct = zero
        potf(i) = sqrt(fluct)
c 1-----
  200 continue
c
c-----  save collected total dispersion interaction
c       and fluctuation on da-file "neqdis"
c
      call dawrit(idafdis,iodadis,pot,ncl,index,navdis)
      call dawrit(idafdis,iodadis,potf,ncl,index+1,navdis)
c
c-----  repulsion interactions
c
      call daread(idafdrf,iodadrf,pot,ncl,163)
      call daread(idafdrf,iodadrf,pot2,ncl,164)
c
c-----  average collected total potential and square
c
      do 300, i = 1, ncl
c 1-----
        pot(i) = pot(i)/(ntot+1)
        pot2(i) = pot2(i)/(ntot+1)
c
c  -----  calculate fluctuation potential
c
        fluct = pot2(i) - pot(i)**2
        if (fluct .lt. thresh) fluct = zero
        potf(i) = sqrt(fluct)
c 1-----
  300 continue
c
c-----  save collected total repulsion interaction
c       and fluctuation on da-file "neqrep"
c
      call dawrit(idafrep,iodarep,pot,ncl,index,navrep)
      call dawrit(idafrep,iodarep,potf,ncl,index+1,navrep)
c
c-----  qm repulsion interactions
c
      call daread(idafdrf,iodadrf,pot,1,165)
      call daread(idafdrf,iodadrf,pot2,1,166)
c
c-----  average collected total potential and square
c
      do 400, i = 1, 1
c 1-----
        pot(i) = pot(i)/(ntot+1)
        pot2(i) = pot2(i)/(ntot+1)
c
c  -----  calculate fluctuation potential
c
        fluct = pot2(i) - pot(i)**2
        if (fluct .lt. thresh) fluct = zero
        potf(i) = sqrt(fluct)
c 1-----
  400 continue
c
c-----  save collected total dispersion interaction
c       and fluctuation on da-file "neqrqm"
c
      call dawrit(idafrqm,iodarqm,pot,1,index,navrqm)
      call dawrit(idafrqm,iodarqm,potf,1,index+1,navrqm)
c
c----- classical energy costs
c
      call daread(idafdrf,iodadrf,pot,3,169)
      call daread(idafdrf,iodadrf,pot2,3,170)
c
c-----  average collected total potential and square
c
      do 500, i = 1, 3
c 1-----
        pot(i) = pot(i)/(ntot+1)
        pot2(i) = pot2(i)/(ntot+1)
c
c  -----  calculate fluctuation potential
c
        fluct = pot2(i) - pot(i)**2
        if (fluct .lt. thresh) fluct = zero
        potf(i) = sqrt(fluct)
c 1-----
  500 continue
c
c-----  save collected total dispersion interaction
c       and fluctuation on da-file "neqrqm"
c
      call dawrit(idafcst,iodacst,pot,3,index,navcst)
      call dawrit(idafcst,iodacst,potf,3,index+1,navcst)
c
      if (field .ne. ' ') then
c 1-----
c  -----  read collected total reaction potential and square
c
        call daread(idafdrf,iodadrf,pot,n,153)
        call daread(idafdrf,iodadrf,pot2,n,154)
c
        if (itwoeps .eq. 1) then
c   2-----
c    -----  read collected optic reaction potential and square
c
          call daread(idafdrf,iodadrf,poto,n,155)
          call daread(idafdrf,iodadrf,poto2,n,156)
c
          do 1100, i = 1, n
c     3-----
            poto(i) = poto(i)/(ntot+1)
            poto2(i) = poto2(i)/(ntot+1)
c
c      -----  calculate fluctuation potential
c
            fluct = poto2(i) - poto(i)**2
            if (fluct .lt. thresh) fluct = zero
            potf(i) = sqrt(fluct)
c   2-----
 1100   continue
c
c    -----  save collected optic reaction potential and square on
c           da-file "neqino"
c
          call dawrit(idafino,iodaino,poto,n,index,navino)
          call dawrit(idafino,iodaino,potf,n,index+1,navino)
c   2-----
        else
c   2-----
          call clear (poto,n)
          call clear (poto2,n)
c   2-----
        endif
c
c  -----  average collected reaction potential and square
c
        do 1200, i = 1, n
c   2-----
          pot(i) = pot(i)/(ntot+1)
          pot2(i) = pot2(i)/(ntot+1)
c
c    -----  calculate fluctuation potential
c
          fluct = pot2(i) - pot(i)**2
          if (fluct .lt. thresh) fluct = zero
          potf(i) = sqrt(fluct)
c   2-----
 1200   continue
c
c  -----  save collected total external reaction potential
c         and fluctuation on da-file "neqind"
c
        call dawrit(idafind,iodaind,pot,n,index,navind)
        call dawrit(idafind,iodaind,potf,n,index+1,navind)
c
c  -----   polarisation energy cost
c
        call daread(idafdrf,iodadrf,pot,1,167)
        call daread(idafdrf,iodadrf,pot2,1,168)
c
        pot(1) = pot(1)/(ntot+1)
        pot2(1) = pot2(1)/(ntot+1)
        fluct = pot2(1) - pot(1)**2
        if (fluct .lt. thresh) fluct = zero
        potf(1) = sqrt(fluct)
c
        call daread(idafpol,iodapol,optneq,7,index)
c
        optneq(6) = int(zero)
        optneq(1) = pot(1)
c
        call dawrit(idafpol,iodapol,optneq,7,index,navpol)
c
        optneq(1) = potf(1)
        call dawrit(idafpol,iodapol,optneq,7,index+1,navpol)
c 1-----
      endif
c
c-----  copy job name on "neqtex", adding appropriate mc
c       identifiers
c
      call neqout(textneq,realneq,maxneq)
      call dawrit(idafind,iodaind,realneq,maxneq*10,1,navind)
c
      call daread(idafinf,iodainf,optneq,1,index)
      call dawrit(idafinf,iodainf,optneq,1,index+1,navinf)
c
      return
      end
      subroutine neqout(textneq,realneq,maxneq)
      implicit REAL (a-h,o-z)
      character*80 textneq(*)
      character*8 tmpnex
      dimension realneq(*)
c
      ii = 0
      do i = 1, maxneq
       iii = 1
       do loop = 1,10
         tmpnex = textneq(i)(iii:iii+7)
         read (tmpnex,'(a8)') realneq( ii + loop)
         iii = iii + 8
       enddo
       ii = ii + 10
      enddo
c
      return
      end
      subroutine neqin(textneq,realneq,maxneq)
      implicit REAL (a-h,o-z)
      character*80 textneq(*)
      character*8 tmpnex
      dimension realneq(*)
      ii = 0
      do i = 1, maxneq
       iii = 1
       do loop = 1, 10
         write(tmpnex,'(a8)') realneq(ii + loop)
         textneq(i)(iii:iii+7) = tmpnex
         iii = iii + 8
       enddo
      ii = ii + 10
      enddo
c
      return
      end
