      subroutine solvset
c------
c      sets solvent total and optic dielectric constants from solvent da
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/sizescon)
INCLUDE(../m4/common/connolly)
c
INCLUDE(comdrf/iofil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfbem)
INCLUDE(comdrf/drfin)
c
      character*32 solnams
      dimension solnams(50)
      dimension epst(50), epso(50), densty(50)
      REAL molmass(50)
c
      character *8 errmsg(3)
c
      REAL Navog,molmas
c
      data solnams /'water','formamide','ethylene_glycol','methanol',
     1 'n-methylformamide','2-methoxyethanol','n-methylacetamide',
     2 'ethanol','acetic_acid','propanol','isopropanol','nitromethane',
     3 'acetonitrile','dimethyl_sulphoxide','dmso','t-butanol',
     4 'dimethyl_formamide','dimethyl_acetamide','acetone',
     5 'nitrobenzene','dichloromethane','pyridine','acetophenone',
     6 'chloroform','1-2-dimethoxyethane','ethyl_acetate',
     7 'tetrahydrofuran','1-4-dioxane','diethyl_ether','benzene',
     8 'carbon_disulphide','carbon_tetrachloride','cyclohexane',
     9 'n-hexane','hexane','vacuum','x',13*' '/
c
      data epst /78.5d0,109.5d0,29.4d0,32.6d0,182.4d0,15.9d0,175.7d0,
     +           24.3d0,6.2d0,20.1d0,18.3d0,38.6d0,37.5d0,48.9d0,48.9d0,
     +           12.2d0,36.7d0,37.8d0,20.7d0,34.8d0,8.9d0,12.3d0,17.4d0,
     +            4.7d0,7.0d0,6.0d0,7.4d0,2.2d0,4.2d0,2.3d0,2.6d0,2.2d0,
     +            2.0d0,1.9d0,1.9d0,1.0d0,14*0.0d0/
c
      data epso /1.777d0,0.0d0,0.0d0,1.766d0,0.0d0,0.0d0,0.0d0,1.853d0,
     +         0.0d0,1.918d0,0.0d0,1.942d0,1.807d0,2.179d0,2.179d0,
     +         0.0d0,0.0d0,0.0d0,1.846d0,2.422d0,0.0d0,0.0d0,2.360d0,
     +         0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.253d0,0.0d0,
     +         2.150d0,0.0d0,1.890d0,1.890d0,1.0d0,14*0.0d0/
c
      data molmass /18.01528d0,45.04d0,62.07d0,32.04d0,59.07d0,76.1d0,
     +         73.09d0,46.07d0,60.05d0,60.1d0,60.1d0,61.04d0,
     +         41.05d0,78.13d0,78.13d0,74.12d0,73.09d0,87.12d0,
     +         58.08d0,123.11d0,84.93d0,79.1d0,120.15d0,119.38d0,
     +         90.12d0,88.11d0,72.11d0,88.11d0,74.12d0,78.11d0,
     +         76.13d0,153.82d0,84.16d0,86.18d0,
     +         86.18d0,15*0.0d0/
c
      data densty /1.0d0,1.1334d0,1.1088d0,0.7914d0,1.011d0,0.9647d0,
     +        0.9571d0,0.7893d0,1.0492d0,0.8035d0,0.7855d0,
     +        1.1371d0,0.7857d0,1.1014d0,1.1014d0,0.7887d0,0.9487d0,
     +        0.9366d0,0.7899d0,1.2057d0,1.3266d0,0.9819d0,1.0281d0,
     +        1.4832d0,0.8628d0,0.9003d0,0.8892d0,1.0337d0,0.7138d0,
     +        0.8765d0,1.2632d0,1.5940d0,0.7785d0,0.6603d0,0.6603d0,
     +        15*0.0d0/
      data errmsg /'program','stop in','-solvset'/
      data bohr /0.529177249d00/
      data Navog /6.0221367d23/
      data pi /3.1415927d0/
c
c-----  loop over standard solvents
c
      do 100, i = 1, 36
c 1-----
        if (solnam .eq. solnams(i)) then
c   2-----
          if (eps1 .eq. 1.0d0) then
c     3-----
c           No input value given - get from database
c
            if (ioponly .eq. 1) then
c       4-----
c
              if (epso(i) .eq. 0.0d0) then
c         5-----
                write(iwr,1003) solnam
 1003           format(/,' no optic dielectric implemented for ',a32,/,
     1  ' feel free to get it and insert it ',
     2  'in epso in subroutine solvset',/,
     3  ' or give epsopt on input ')
                call hnderr(3,errmsg)
                return
c         5-----
              endif
c
              eps1 = epso(i)
              goto 123
c       4-----
            else
c       4-----
              if (epst(i) .eq. 0.0d0) then
c         5-----
                write(iwr,1002) solnam
 1002         format(/,' no total dielectric implemented for ',a32,/,
     1  ' feel free to get it and insert it ',
     2  'in epst in subroutine solvset',/,
     3  ' or give epsstat on input ')
                call hnderr(3,errmsg)
                return
c         5-----
              endif
c
              eps1 = epst(i)
c       4-----
            endif
c     3-----
          endif
c
          if (itwoeps .eq. 1) then
c     3-----
            if (eps2 .eq. 1.0d0) then
c       4-----
c           No input value given - get from database
c
              if (epso(i) .eq. 0.0d0) then
c         5-----
                write(iwr,1003) solnam
                call hnderr(3,errmsg)
                return
c         5-----
              endif
c
              eps2 = epso(i)
c        4-----
             endif
c     3-----
          endif
c
          goto 123
c
caleko
c
123   continue
      if (solrad .lt. 0.0d0) then
      if ((ibem.eq.5).or.(ibem.eq.3)) then
        if (molmass(i).eq.0.0d0) then
          write(iwr,1004) solnam
1004      format(/,'no molecular mass implemented for ',a32,/,
     1  ' feel free to get it and insert it ',
     2  ' in molmass in subroutine solvset')
          call hnderr(3,errmsg)
          return
        endif
        molmas=molmass(i)
c
        if (densty(i).eq.0.0d0) then
          write(iwr,1005) solnam
1005      format(/,'no density is implemented for ',a32,/,
     1  ' feel free to get it and insert it',
     2  ' in densty in subroutine solvset')
          call hnderr(3,errmsg)
          return
        endif
        density=densty(i)
        vmol=((3.0d0*molmas)/(4.0d0*pi*Navog*density))
        solrad=(vmol**(1.0d0/3.0d0))/(bohr*1.0d-8)
cxxx rprobe is separate parameter
cxxx (for connolly surface, rp is used!)
        if (inpro .eq. 0) rprobe = solrad
        if (inrp .eq. 0) rp = solrad
      endif
      endif
      return
c   2-----
        endif
c
caleko
c
c 1-----  next solvent name
  100 continue
c
      write(iwr,1006) solnam
 1006 format(/,' solvent name ',a32,' not found:',/,
     1 ' give explicit data in directive dielectric for this solvent')
      call caserr('solvent directive input error detected')
c
      return
      end
      subroutine drfgrp(idrfout,nfirst,nlast,mgrp,near,nextpol,npol,
     1      ngrp,igran,ngrnam,namgrp,fgroup,idipcal,afact,ithole,
     2      agrpc,agrpm,agrpe,ieffpol)
c------
c     * * *  routine forms group polarizability for atoms nfirst-nlast
c            and checks the accuracy of this approximation.
c            if successful,or if grouping is forced
c            (fgroup .true.) the corresponding polarizabilities
c            are replaced by a single 3 x 3  point polarizability
c
c     routines called: clear,distab,drfamat,fdagal,dnrm2,ddot,
c                      drfreda,hatout,linv3
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension dip1(3),dip2(3),centpol(3)
      character *10 namgrp
      logical fgroup
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/auxdrf)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/max)
INCLUDE(comdrf/bas)
c
INCLUDE(../m4/common/drfopt)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/drfamb)
INCLUDE(comdrf/grinf)
c
INCLUDE(comdrf/rad)
c
      REAL xscm
      dimension xscm(15*mxnpts**2+8*mxnpts)
c
      logical amb, group
      logical ocalc
      REAL x0
      dimension thispol(6), q(3), x0(3)
c
      dimension polat(mxnpts)
      dimension polatt(6,mxnpts)
c
      data zero,pt5,one/0.0d00,0.5d00,1.0d00/
      data pt3 /0.33333333333333333333333333d00/
c
c-----  begin
c       set memory partitioning
c
      i10 = 1
c
c-----  calculate unit field from atom -near- in selected points
c       center of polarizabilities in -centpol-
c
      totpol = zero
      imp = mpol(nextpol)
      ipgrp = 1
      np = 0
      call clear(centpol,3)
      ipol = nextpol
      nvagrp = 0
*     write(6,*) 'nfirst,nlast = ', nfirst, nlast
      ocalc = .true.
      do 20, i = nfirst, nlast
c 1-----
        if (i .eq. mpol(ipol)) then
c   2-----
c    ----  count  # valence  electrons of this atom
c
          nvagrp = nvagrp + int(vale(nextpol+ipgrp-1)+0.01d0)
          ipgrp = ipgrp + 1
          np = np + 1
          irow = (np-1)*3
          call distab(c(1,near),xpts(1,i),q,dist)
*         write(6,*) ' near, i, dist = ', near,i,dist
*         write(6,*) 'q,one = ', q, one
          if (dist .le. 1.e-10) then
            write(6,*) 'centre of polarizability on qm atom!'
            ocalc = .false.
          endif
          if (ocalc) dmin3 = (one/dist)**3
          do 10, k = 1, 3
            centpol(k) = centpol(k)+polar(ipol)*xpts(k,i)
            if (ocalc) then
              xscm(i10+irow+k-1) = q(k)*dmin3
            else
              xscm(i10+irow+k-1) = one 
            endif
   10     continue
c
          totpol = totpol + polar(ipol)
          imp = imp + 1
          ipol = ipol + 1
c   2-----
        endif
c
        if (ncutpt(i) .ne. 0) amb = .true.
c 1-----
   20 continue
c
cxxx  if (totpol .lt. 1.e-10) totpol = np
c
      i10 = 1
      do 25, k = 1, 3
   25   centpol(k) = centpol(k)/totpol
c
c-----  polarizability matrix of selected points
c
      np3 = np*3
      nnp3 = np3*(np3+1)/2
c
c-----  partition of core
c
c     field in   -i10-
c     a-mat in   -i20-
c     dipoles in -i30-
c
      i20 = i10 + np3
      i30 = i20 + nnp3
      i40 = i30 + np3
c
      call clear(xscm(i20),nnp3)
      call drfamat
     1   (xscm(i20),nextpol,.false.,np,xpts,mpol,polar,idrfout,afact,
     2    ithole)
c
c  -----  invert -a-
c
c     subroutine linv3p can be replaced by standard 'vectorized'
c     matrix inversion subroutine (nag,linpack etc)
c
      call linv3p(xscm(i20),auxx,1,np3,ier)
      if (idrfout .ge. 3)
     1    call hatout(xscm(i20),np3,np3,3,'group-amat')
c
c  -----  calculate induced moments in selected points
c
      call fdagal(xscm(i10),xscm(i20),xscm(i30),ia,np3)
c
c  -----  contract to single dipole moment
c
      call clear(dip1,3)
      do 40, i = 1, np
        do 40, k = 1, 3
   40 dip1(k) = dip1(k) + xscm(i30+(i-1)*3+k-1)
      dip1l = dnrm2(3,dip1,1)
c
c-----  calculate induction energy
c
      e1 = pt5*ddot(np3,xscm(i10),1,xscm(i30),1)
c     e1 = pt5*adotb(xscm(i10),xscm(i30),np3)
c
c-----  make the group polarisability matrix -amat- block diagonal
c       to get effective atom polarisabilities for use in
c       calculation of dispersion energy
c
      call effpol(np,np3,ia,xscm(i20),xscm(i40),
     1            polat,polatt)
c
c-----  contract polarizability matrix to 3x3 point polarizability
c
      call drfreda(xscm(i20),thispol,np)
      if (idrfout .ge. 2) call hatout(thispol,3,3,3,'grouppol')
c
c-----  field of atom -near- in -centpol-
c
      if (ocalc) then
      call distab(c(1,near),centpol,q,dist)
      dmin3 = (one/dist)**3
      do 45, k = 1, 3
   45   q(k) = q(k)*dmin3
c
c-----  calculate induced moment
c
       call fdagal(q,thispol,dip2,ia,3)
       dip2l = dnrm2(3,dip2,1)
       e2 = pt5*ddot(3,q,1,dip2,1)
c      e2 = pt5*adotb(q,dip2,3)
c
c-----  check accuracy
c
      diffe = e1 - e2
      diffm = dip1l - dip2l
c     cos = adotb(dip1,dip2,3)/
      cos = ddot(3,dip1,1,dip2,1)/
     +      (dip1l*dip2l)
      else
        diffe = agrpe
        diffm = agrpm
        diffc = agrpc
      endif
c
      if (field(5:) .eq. ' ') then
c 1-----
        call clear(thispol,6)
        diffe = zero
        diffm = zero
        cos = one
c 1-----
      endif
c
      iwarn = 0
      if (abs(diffe) .ge. agrpe) then
        write(iwr,7011) 
 7011   format(' induction energy criterion ',
     1  'for grouping not met')
        if (fgroup) iwarn=1
      endif
      if (abs(diffm) .ge. agrpm) then
        write(iwr,7021) 
 7021   format(' induced dipole criterion ',
     1  'for grouping not met')
        if (fgroup) iwarn=1
      endif
      if (cos .le. agrpc) then
        write(iwr,7031) 
 7031   format(' cosine induced dipoles ',
     1  'for grouping not met')
        if (fgroup) iwarn=1
      endif
c
      if (iwarn .eq. 1) then
        write(iwr,7041)
 7041 format(' WARNING: grouping forced by user,',
     1 ' but grouping criteria not met')
        if (idrfout .lt. 2) write(iwr,7051)
 7051 format(' grouping information may be obtained by ',
     1 'specifying  the REACT subdirective DRFOUT MORE')
      endif
c
      group = .false.
      if (((abs(diffe) .lt. agrpe) .and.
     1     (abs(diffm) .lt. agrpm) .and.
     2     (cos .gt. agrpc)) .or. fgroup) then
c 1-----
c  -----  group polarisability is to be constructed:
c         shift coordinates down in the list and change -nfirst-
c         as last point read
c
        group = .true.
        ngrp = ngrp + 1
        nlast = nlast + 1
        do 50, k = 1, 3
   50     xpts(k,nlast) = xpts(k,nfirst)
        chrg(nlast) = chrg(nfirst)
        nxcent(nlast) = nxcent(nfirst)
        radext(nlast) = radext(nfirst)
        alfext(nlast) = alfext(nfirst)
        imemgrp(nlast) = ngrp
c
        k = 0
        ipol = nextpol
        do 52, i = nfirst, nlast
c   2-----
          igranl(i) = igran
          if (i .eq. mpol(ipol)) then
c     3-----
            k = k + 1
            valat(i) = vale(nextpol-1+k)
            if (ieffpol .eq. 0) then
              atpol(i) = polar(ipol)
              ij = 0
              do 501, j = 1, 3 
                do 501, l = 1, j
                  ij = ij + 1
                  if (l.eq.j) then
                    atpolt(ij,i) = polar(ipol)/3.0d0
                  else
                    atpolt(ij,i) = zero
                  endif
  501         continue
            else
              atpol(i) = polat(k)
              do 51, j = 1, 6
                atpolt(j,i) = polatt(j,k)
   51         continue
            endif
            ipol = ipol + 1
          else
            valat(i) = zero
            atpol(i) = zero
            call clear(atpolt(1,i),6)
c     3-----
          endif
cxxx
          imemgrp(i) = ngrp
cxxx
c   2-----
   52   continue
c
        valat(nlast) = valat(nfirst)
        atpol(nlast) = atpol(nfirst)
        do 53, k = 1, 6
          atpolt(k,nlast) = atpolt(k,nfirst)
   53   continue
c
c  -----  save -centpol- as new center for group polarizability
c         in position -nfirst-
c
        do 55, k = 1, 3
   55   xpts(k,nfirst) = centpol(k)
        nxcent(nfirst) = 'group '//namgrp
        chrg(nfirst) = zero
        valat(nfirst) = zero
        atpol(nfirst) = zero
        do 56, k = 1, 6
          atpolt(k,nfirst) = zero
   56   continue
c
        imemgrp(nfirst) = ngrp
c
c  -----  discard polarizabilities and valence elns
c         of selected points
c
        npol = nextpol
        mpol(npol) = nfirst
        vale(npol) = nvagrp
        do 57, i = npol+1, maxpol
          mpol(i) = 0
          vale(i) = zero
   57   continue
c
c  -----  save group polarizability
c
        ngrpol = ngrpol + 1
        polav = zero
        do 60, k = 1, 3
          kk = ia(k) + k
          polav = polav + thispol(kk)
          do 60, l = 1, k
          kl = ia(k) + l
   60   grpol(kl,ngrpol) = thispol(kl)
        polar(npol) = pt3*polav
        igrpol(ngrpol) = npol
        nvalgrp(ngrpol) = nvagrp
        radext(nfirst) = polar(npol)**pt3*afact
        alfext(nfirst) = polar(npol)
      else
        write(iwr,*) '- - - criteria for grouping not met,',
     1  ' grouping atoms',nfirst,' to',nlast,' not effected'
      endif
      if (idrfout .ge. 2) then
        write(iwr,8001)
 8001   format(/,'      method ',t35, 'dipole', t60,
     1  ' induction energy')
        write(iwr,8011) 'separate pols:',dip1,e1
 8011   format(1x,a20,4(1x,e12.6))
        write(iwr,8011) 'group polarizability',dip2,e2
cxxx    write(iwr,8021) diffe,diffm,cos
cxxx    write(iwr,*) 'nvagrp ', nvagrp
        if (group) write(iwr,9000) ngrpol,np,nfirst,nlast
      endif
 9000 format(
     1/'- - -  group: ',2i4,' polarizabilites from',i4,' to',i4,
     2 ' contracted'/)
c
c-----  update ambiguous points arrays if necessary
c
      if (amb) then
        if (ncutpt(nfirst) .ne. 0) ncutpt(nlast) = ncutpt(nfirst)
        do 70, i = nfirst + 1, nlast
          if (ncutpt(i) .ne. 0) then
            nambpt(ncutpt(i)) = i
            ncutpt(nfirst) = ncutpt(i)
          endif
   70   continue
      endif
c
c-----  compute group dipole and quadrupole moments if required
c
      if ((group) .and. (idipcal .eq. 1)) then
c 1-----
        do 100, k = 1, 3
          x0(k) = zero
  100   continue
c
        tmass = zero
c
        do 200, i = nfirst + 1, nlast
c   2-----
          amass = amas(nxcent(i)(:2))
          tmass = tmass + amass
          do 300, k = 1, 3
            x0(k) = x0(k) + amass*xpts(k,i)
  300     continue
c   2-----
  200   continue
c
        do 400, k = 1, 3
          x0(k) = x0(k)/tmass
  400   continue
c
        write(iwr,9001) x0
 9001   format(/,'  current group centre of mass:',3f10.5)
c
        call dipcal(nlast-nfirst,x0,xpts(1,nfirst+1),
     1              chrg(nfirst+1),idrfout)
c 1-----
      endif
c
      return
      end
      subroutine effpol(npol,np3,ia,amat,bmat,atpol,atpolt)
c------
c      defines effective atomic polarisabilities &
c      polarisability tensors from the group relay matrix
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/bas)
c
      dimension amat(np3*(np3+1)/2)
      dimension bmat(np3,np3)
      dimension atpol(npol)
      dimension atpolt(6,npol)
      dimension ia(npol*(npol+1)/2)
c
      data pt3 /0.33333333333333333333333333333333d00/
c
c-----  begin
c
      call hexpand(amat,bmat,np3,np3)
c
c-----  loop over atomic polarisabilities - diagonal 3x3 blocks
c
      do 100, i = 1, npol
c 1-----
c  -----  loop over atomic polarisabilities - off-diagonal 3x3 blocks
c
        do 200, j = 1, npol
c   2-----
          if (j .eq. i) goto 200
          do 300, k = 1, 3
c     3-----
            do 400, l = 1, 3
c       4-----
              bmat((i-1)*3+k,(i-1)*3+l) = bmat((i-1)*3+k,(i-1)*3+l)
     1      + bmat((j-1)*3+k,(i-1)*3+l)
c       4-----
  400       continue
c     3-----
  300     continue
c   2-----
  200   continue
c 1-----
  100 continue
c
c-----  extract isotropic and anisotropic effective atomic polarisabilit
c
      call clear (atpol,npol)
      do 1100, i = 1, npol
c 1-----
        do 1200, k = 1, 3
c   2-----
          atpol(i) = atpol(i) + bmat((i-1)*3+k,(i-1)*3+k)
          do 1300, l = 1, k
            kl = ia(k) + l
            atpolt(kl,i) = bmat((i-1)*3+k,(i-1)*3+l)
 1300     continue
c   2-----
 1200   continue
c 1-----
 1100 continue
c
      do 2100, i = 1, npol
        atpol(i) = pt3*atpol(i)
 2100 continue
c
      return
      end
      subroutine dipcal(n,x0,x,q,idrfout)
c------
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension x(3,n)
      dimension x0(3), q(n)
c
INCLUDE(comdrf/iofil)
c
      dimension dip(3)
      dimension th(3,3)
      dimension xr(3,50)
      character *8 errmsg(3)
c
      data zero, pt5, onept5 /0.0d00, 0.5d00, 1.5d00/
      data errmsg /'program', 'stop in','-dipcal-'/
c
      if (n .gt. 50) then
        write(iwr,1000) n
 1000   format(/,' group too large, increase size of array ',
     1  'xr to at least',i4, ' in subroutine dipcal')
        call hnderr(3,errmsg)
        return
      endif
c
      do 100, i = 1, 3
        dip(i) = zero
        do 200, j = 1, n
          xr(i,j) = x(i,j) - x0(i)
  200   continue
        do 300, j = 1, 3
          th(i,j) = zero
  300   continue
  100 continue
      charge = zero
c
      do 400, i = 1, n
        charge = charge + q(i)
        r2 = ddot(3,xr(1,i),1,xr(1,i),1)
c       r2 = adotb(xr(1,i),xr(1,i),3)
c
        do 500, k = 1, 3
          dip(k) = dip(k) + q(i)*xr(k,i)
c
          do 600, l = 1, 3
            th(l,k) = th(l,k) + onept5*q(i)*xr(l,i)*xr(k,i)
  600     continue
c
  500   continue
c
        do 450, k = 1, 3
          th(k,k) = th(k,k) - pt5*q(i)*r2
  450   continue
c
  400 continue
c
      diptot = dnrm2(3,dip,1)*2.54158d0
c
      write(iwr,1001) charge,diptot,dip
 1001 format(/' total charge ',f10.6/' dipole moment',f10.6,
     +       ' debye   (',3f10.6,') au')
      write(iwr,1002) ((th(l,k), k= 1,l), l=1,3)
 1002 format(/' quadrupole moment wrt centre of mass, lower triangle ',/
     1 ' thxx = ', f10.6, ' a.u.',/
     2 ' thxy = ', f10.6, ' a.u.',/
     3 ' thyy = ', f10.6, ' a.u.',/
     4 ' thxz = ', f10.6, ' a.u.',/
     5 ' thyz = ', f10.6, ' a.u.',/
     6 ' thzz = ', f10.6, ' a.u.')
c
      return
      end
      subroutine addamb(elnum,charge,cxyz,name,ibas,iecpx)
c------
c      drives addition of nuclear centres and basis functions/ecp's
c      to qm atom list and changed classical properties
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension cxyz(3)
      character*16 name
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ecp)
c
INCLUDE(comdrf/defpar)
c
      character*2 names(100)
c
      data names/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     *           'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     *           'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     *           'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     *           'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     *           'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     *           'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     *           'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     *           'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     *           'pa','u ','e ','qq','xx',5*'oi'/
c
      data zero /0.0d00/
c
c-----  begin
c
c-----  copy incoming data into common blocks
c-----  coordinates
c
      do 100, i = 1, 3
        cords(i) = cxyz(i)
  100 continue
c
c-----  effective charge; this will be made the nuclear charge of the
c       added atom in the quantum system. the reason for this is that
c       its field will only then be treated exact (whereas expansion
c       would give infinite fields
c
      effchr = charge
c
c-----  chemical symbol & definition of name
c
      namei = name
c
c-----  find nuclear charge
c
      do 200, i = 1, 100
        if(name(:2) .eq. names(i)) goto 220
 200  continue
      write(iwr,*) name,'illegal atom name in readat'
      call caserr('error readat 20')
 220  continue
      chargx = dble(i)
      if (name(:2) .eq. 'e ') chargx = zero
c
      if (elnum .lt. zero) then
        ambel = chargx
      else
        ambel = elnum
      endif
c
      extra = .true.
c
      iecp = iecpx
c
c-----  add to qm list
c
      call atoms(0,ibas,iecp)
c
      return
      end
      subroutine pesinpc
c------
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/mollab)
INCLUDE(comdrf/free)
INCLUDE(comdrf/scm)
c
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/grinf)
c
INCLUDE(comdrf/clas)
c
INCLUDE(comdrf/mcinf)
c
      dimension nstp(4)
      dimension pesvec(3), angstp(3)
c
      character*80 line
      character*16 namei
      character*8 namess
      dimension namess(mxgr1+mxnpts)
c
      dimension cgrps(3,mxgr1+mxnpts)
c
      dimension xfrg(3,mxnpts), xfrgn(3,mxnpts)
      dimension ttrmat(3,3), pmom1(3,3), pmom2(3,3)
      dimension vecn(3), xcent2(3), eula(3)
c
      logical redund, out
c
      dimension vec(3), cxyz(6)
      dimension nfrag(2)
c
      character *8 errmsg(3)
c
      namelist /pesc/ npesgrp, intpesc, npesref, iunit, ipesout,
     2         nstp, vecstp, angstp, itest, pesvec
c
      data errmsg /'program ','stop in ','-pesc  -'/
c      data npesgrp, npesref, iunit /0,0,0/
      data npesref, iunit /0,0/
c      data intpesc,ipesout,itest /0,0,0/
      data itest /0/
      data vecstp /0.0d00/
      data angstp /3*0.0d00/
      data nstp /4*0/
      data pesvec /3*0.0d00/
c
      data zero, one /0.0d+00, 1.0d00/
c
      data bohr /.529177249d00/
c
c      data numpes /0/
c
c     comdrf/clas
c
      npesgrp = 0
      intpesc = 0
      ipesout = 0
      numpes  = 0
c
c--------------------------
c         input description for $pesc
c
c         (compare $pes in pesx)
c
c       npesgrp: fragment to be moved in the pes scan
c                this is the serial number of the group as defined
c                in $external (default=0)
c
c       intpesc:  option for pes scan
c                   0: pes scan of group relative to group defined
c                      by npesref. the relative coordinates are the
c                      distance vector and the euler angles (default)
c                   1: pes scan of group along a vector pesvec,
c                      given on input; nstp(1) steps are made of
c                      length vecstp (pesvec is normalised to unit lengt
c                   2: pes scan according to group coordinates
c                      for npesgrp given on consecutive data groups
c                      $extpes that have the same structure as $external
c
c       npesref: fragment serving as reference with respect to which
c                the group defined by npesgrp will be moved if intpesc=0
c                -1: the quantum fragment serves as reference
c                (default=0)
c
c                conflicts will result in error message and
c                abortion of the program
c
c       ipesout: output flag
c                   0: minimum output (default)
c                >= 5: analysis for each configuration
c
c--------------------------
c
c-----  begin
c
c-----  read namelist $pesc
c
      rewind (ir)
      read (ir,pesc,end=10)
c
      out = ipesout .ge. 3
c
      if ((npesgrp .eq. 0) .or. (npesgrp .gt. ngrpol)) then
c 1-----
        write (iwr, 111) npesgrp
  111   format (/, ' invalid group number: ',i4,' to be moved in',
     2  ' classical pes scan, please check npesgrp')
        call hnderr (3,errmsg)
c 1-----
      endif
c
      if (iunit .eq. 1) vecstp = vecstp/bohr
      lpes = 0
c
      if (intpesc .eq. 0) then
c 1-----
c  -----  pes scan according to relative coordinates
c
        if ((npesref .eq. 0) .or. (npesref .eq. npesgrp)
     2    .or. (npesref .gt. ngrpol)) then
c   2-----
        write (iwr, 1011) npesref
 1011   format (/, ' invalid reference group number: ',i4,' in',
     2  ' classical pes scan, please check npesref')
        call hnderr (3,errmsg)
c   2-----
        endif
c
        if ((vecstp .eq. zero) .and. (angstp(1) .eq. zero)
     1     .and. (angstp(2) .eq. zero) .and. (angstp(3) .eq. zero))
     2     then
c   2-----
          write(iwr,1021)
 1021     format(/, ' all relative coordinates for pes scan are',
     1    ' equal zero; check $pesc input ')
          call hnderr(3,errmsg)
c   2-----
        endif
c
        if ((nstp(1) .eq. 0) .and. (nstp(2) .eq. 0)
     1     .and. (nstp(3) .eq. 0) .and. (nstp(4) .eq. 0))
     2     then
c   2-----
          write(iwr,1031)
 1031     format(/, ' all numbers of steps for pes scan are',
     1    ' equal 0; check $pesc input ')
          call hnderr(3,errmsg)
c   2-----
        endif
c
        if (ngrppt(npesgrp) .gt. mxnpts) then
c   2-----
          write(iwr,1041)
 1041     format(/,
     2  ' group size of group to be moved exceeds maximum: ',/,
     3  ' reduce group size or enlarge mxnpts in ctlnew.f',
cmw     4  ' and drf/dimpar')
     4  ' and rfin/sizesrf')
          call hnderr(3,errmsg)
c   2-----
        endif
c
c  -----  get information on reference group
c
        nats = 0
c
        if (npesref .eq. -1) then
c   2-----
          do 1100, i = 1, nat
c     3-----
            nats = nats + 1
            namess(nats) = anam(i)
c
            do 1010, k = 1, 3
              cgrps(k,nats) = c(k,i)
 1010       continue
c     3-----
 1100     continue
c
          nfrag(1) = nat + 1
c   2-----
        else
c   2-----
          do 1300, i = 1, ngrppt(npesref)
c     3-----
            nats = nats + 1
            namess(nats) = nxcent(igrpst(npesref)+i)(1:8)
c
            do 1210, k = 1, 3
              cgrps(k,nats) = xpts(k,igrpst(npesref)+i)
 1210       continue
c     3-----
 1300     continue
c
          nfrag(1) = ngrppt(npesref) + 1
c   2-----
        endif
c
c  -----  get coordinates of the group to be moved around into
c         cgrps (relative to global origin)
c
        do 1500, i = 1, ngrppt(npesgrp)
c   2-----
          nats = nats + 1
c
          namess(nats) = nxcent(igrpst(npesgrp)+i)(1:8)
c
          do 1410, k = 1, 3
            cgrps(k,nats) = xpts(k,igrpst(npesgrp)+i)
 1410     continue
c   2-----
 1500   continue
c
        nfrag(2) = nfrag(1) + ngrppt(npesgrp) - 1
c
c  -----  calculate centres of mass of both fragments
c         as well as their principle moments of inertia
c
cnot        call fragcal(nats,cgrps,xcent2,nfrag,vec,dist,pmom1,pmom2,
cnot     1               eula,namess,ipesout)
c
c   this routine seems at various points to assume that 
c   fragcal has been called e.g. dist and eula in the print below
c   now the print statement is commented out ...
c       if (out) write(iwr,1311) dist, (eula(k), k = 1, 3)
c1311   format(/,' start geometry inter group distance = ',f8.4,'  au',
c    1         /,' start geometry euler angle theta    = ',f8.4,' deg',
c    2         /,' start geometry euler angle phi      = ',f8.4,' deg',
c    3         /,' start geometry euler angle chi      = ',f8.4,' deg')
c
c  -----  normalise the distance vector between the groups
c
        call normliz (vec,vec,3)
c
c  -----  get coordinates of the group to be moved around into
c         xfrg (relative to centre of mass)
c
        nfrag2 = 0
        do 1700, i = 1, ngrppt(npesgrp)
c   2-----
          nfrag2 = nfrag2 + 1
c
          do 1610, k = 1, 3
c   NOTE: xcent2 appears to be set by call to fragcal above, but this 
c   is commented out. For the moment, do the same below
cnot        xfrg(k,nfrag2) = xpts(k,igrpst(npesgrp)+i) - xcent2(k)
            xfrg(k,nfrag2) = xpts(k,igrpst(npesgrp)+i) 
 1610     continue
c   2-----
 1700   continue
c
c  -----  calculate new coordinates
c
        if (nstp(1) .eq. 0) nstp(1) = 1
        if (nstp(2) .eq. 0) nstp(2) = 1
        if (nstp(3) .eq. 0) nstp(3) = 1
        if (nstp(4) .eq. 0) nstp(4) = 1
c
        do 400, ithet = 1, nstp(2)
c   2-----
          thetan = angstp(1)*(ithet-1)
c
          do 500, iphi = 1, nstp(3)
c     3-----
            phin = angstp(2)*(iphi-1)
c
            redund = .false.
c
            do 600, ichi = 1, nstp(4)
c       4-----
              chin = angstp(3)*(ichi-1)
c
c        -----  calculate euler transformation matrix in the local
c               frame
c        -----  remove redundancies
c
              if (thetan .eq. zero) then
c         5-----
                chin = zero
                if (redund) then
                  goto 600
                endif
                redund = .true.
c         5-----
              endif
c
cnot              call eulmat(thetan,phin,chin,ttrmat)
c
              if (out) call hatout(ttrmat,3,3,2,'trmat1')
c
c        -----  transform euler matrix to cartesian frame
c
cnot              call simtrn(ttrmat,pmom2,ttrmat,3,3,3)
c
              if (out) call hatout(ttrmat,3,3,2,'trmat2')
c
c        -----  perform rotation over euler angles
c
              call mattrns(ttrmat,xfrg,xfrgn,3,3,nfrag2)
c
              if (out) call hatout(xfrgn,3,nfrag2,2,'xfrgn')
c
              do 700, idist = 1, nstp(1)
c         5-----
                vecn(1) = vec(1)*vecstp*(idist-1)
                vecn(2) = vec(2)*vecstp*(idist-1)
                vecn(3) = vec(3)*vecstp*(idist-1)
c
c          -----  new coordinates of group members
c
                do 800, iat = 1, ngrppt(npesgrp)
c           6-----
                  do 900, k = 1, 3
c             7-----
                    lpes = lpes + 1
c   NOTE: xcent2 appears to be set by call to fragcal above, but this 
c   is commented out. For the moment, do the same below
cnot                xscm(lpes) = xfrgn(k,iat) + xcent2(k) + vecn(k)
                    xscm(lpes) = xfrgn(k,iat) + vecn(k)
c             7-----
  900             continue
c           6-----
  800           continue
c
c          -----  count number of pes points
c
                numpes = numpes + 1
c         5-----  next dist
  700         continue
c       4-----  next chi
  600       continue
c     3----- next phi
  500     continue
c    2-----  next theta
  400   continue
c 1-----
      else if (intpesc .eq. 1) then
c 1-----
        if (nstp(1) .eq. 0) then
c   2-----
          write (iwr,2011)
 2011     format (/, ' number of steps along vector in pes scan ',
     2    'equals 0, check nstp(1) in $pesc')
          call hnderr (3,errmsg)
c   2-----
        endif
c
        if (vecstp .eq. zero) then
c   2-----
          write (iwr,2021)
 2021     format (/, ' length of step along vector in pes scan ',
     2    'equals zero, check vecstp in $pesc')
          call hnderr (3,errmsg)
c   2-----
        endif
c
        if ((pesvec(1) .eq. zero) .and.
     2      (pesvec(2) .eq. zero) .and.
     3      (pesvec(3) .eq. zero)) then
c   2-----
          write (iwr,2031)
 2031     format (/, ' all vector components in pes scan ',
     2    'equal zero, check pesvec in $pesc')
          call hnderr (3,errmsg)
c   2-----
        endif
c
        call normliz(pesvec,vec,3)
c
        write(iwr,2041) vec(1), vec(2), vec(3), nstp(1), vecstp
 2041   format(/,
     1  ' the group will be moved along the vector: ',3f10.5,
     2 /,' in ',i3,' steps of length ',f8.4,' a.u.')
c
        numpes = nstp(1)
c
        istrt = igrpst(npesgrp)
c
c  -----  loop over configurations
c
        do 2100, i = 1, numpes
c   2-----
c    -----  calculate new coordinates according to move
c
          do 2010, k = 1, ngrppt(npesgrp) + 1
c     3-----
            do 2001, l = 1, 3
c       4-----
              lpes = lpes + 1
              xscm(lpes) = xpts(l,istrt+k-1) + vec(l)*i*vecstp
c       4-----
 2001       continue
c     3-----
 2010     continue
c   2-----
 2100   continue
c 1-----
      else if (intpesc .eq. 2) then
c 1-----
        maxbl = maxblnk
c
        rewind (ir)
c
 3110   continue
c
        read (ir,3111,err=4100) line
 3111   format (a80)
c
        if (line(1:8) .eq. ' $extpes') then
c   2-----
          ngrp = 0
c
 3210     continue
c
          read (ir,3111) line
c
          if (line(:5) .eq. ' $end') then
c     3-----
            goto 3310
c     3-----
          endif
c
          if (line(:6) .eq. 'blanks') then
c     3-----
            read(line(7:),3221) maxbl
 3221       format (i2)
            write(iwr,*) ' = = = = = maxbl changed to  ',maxbl
            goto 3210
c     3-----
          endif
c
          if (line(:1) .eq. '"') then
c     3-----
            write (iwr,3111) line
            goto 3210
c     3-----
          endif
c
          if (line(:5) .eq. 'group') then
c     3-----
            goto 3210
c     3-----
          endif
c
          ngrp = ngrp + 1
          if (ngrp .gt. ngrppt(npesgrp)) then
c     3-----
            write (iwr,3251) numpes+1
 3251       format (/, ' group information on ',i4,' th instant of',
     2      ' $extpes incompatible with $external, please check')
            call hnderr (3,errmsg)
c     3-----
          endif
c
          call freerd(line,namei,cxyz,6,maxbl)
c
          do 3230, k = 1, 3
            lpes = lpes + 1
            xscm(lpes) = cxyz(k+1)
 3230     continue
c
          goto 3210
c
 3310     numpes = numpes + 1
c   2-----
        endif
c
        goto 3110
c
c  -----  skip to here at end of file
c
 4100   continue
c 1-----
      endif
c
c-----  write configuration information on da-file 31, record 98
c
      call dawrit (idafdrf,iodadrf,xscm,lpes,98,navdrf)
c
c-----  calculate length of each set of new coordinates
c
      lpesc = lpes/numpes
c
      return
c
   10 continue
      write (iwr, 9011)
 9011 format (/, ' no namelist $pesc found on input, please check ')
      call hnderr (3,errmsg)
c
      return
c
      end
      subroutine drfexpc(xscm)
c------
c       in this routine the expansion centra of the qm system
c       are defined.
c       they are the positions of the nuclei plus (optional)
c       the centre of nuclear charge and centra read from input
c       as well as other extra to come out of the assignment of
c       overlap distributions (by assgnij, according to
c       the options set in $expanx)
c
c      namelist $expanx
c
c      iasexp: option for assignment of overlap distributions
c
c          -1: assign all two-centre overlap distributions to the centre
c              of charge
c           0: if no unambiguous assignment can be made, assign
c              distribution to centre of charge
c           1: if no unambiguous assignment can be made, let it go
c              astray, or define new centra (default)
c           2: assign according to least squares fit to potential
c              due to overlap charge and dipole
c
c      ifitx : option for centra at which the potential of overlap
c              distributions is to be calculated for fitting procedure
c              (with iasexp=2 only)
c
c           0: use the same centra as defined by iadexp (default)
c           1: use atomic centra only
c           2: use atomic centra + centre of charge
c           3: use atomic centra + centre of charge + centra read
c              from input ($fitx)
c           4: use atomic centra + centra read from input ($fitx)
c           5: use centra read from input ($fitx) only
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  common blocks
c
      common/expanx/iasexp,ifitx,iradial
INCLUDE(comdrf/assign)
      logical obeen
      common/nottwi/obeen
cINCLUDE(comdrf/iofil)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/ihelp)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/tim)
       common/restar/nprint
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
cahv new call to asgnrf
INCLUDE(../m4/common/nshel)
INCLUDE(../m4/common/symtry)
      common/junk3/ptr(3,144),dtr(6,288),ftr(10,480)
c
      dimension xscm(*)
c
c-----  local variables
c
      dimension xex(3,maxex)
      dimension cm(3)
      logical defnew, each, asscen, oradial
      logical oassign
      character*80 text
      character *8 errmsg(3)
      logical out
      data pt5 /0.5d00/
c
c      namelist /expanx/iasexp,ifitx
c
c      data iasexp, ifitx  /1,0/
      data errmsg /'program', 'stop in', '-drfexpc'/
c
c-----  begin
c
      obeen = .true.
      oassign = .true.
      if ((ibeen .gt. 0) .and. (keepij .ne. 1)) 
     1   oassign = .false.
c
      noprt = 1
c
      out = noprt.eq.1.and.nprint.eq.3
c
c      rewind ir
c      read(ir,expanx,end=999)
c      goto 1
c
c  999 continue
c      write(iwr,*) 'namelist $expanx not found on input: please check'
c      call hnderr(3,errmsg)
c
c-----  check validity of options
c
    1 continue
      if ((icmexp .eq. 0) .and. (iasexp .eq. -1)) then
c       write(iwr,1001) icmexp, iasexp
c1001   format(/,' input parameters icmexp= ',i2,
c    1         ' and iasexp= ',i2,' conflict',/,
c    2         ' please change them on input: $drf and $expanx ')
        write(iwr,1001)
 1001   format(/,' REACT input subdirectives',
     1           ' EXPANDCM and ASSIGN conflict')
        call hnderr(3,errmsg)
      endif
      if (iasexp .eq. 2) then
        if (iasexp .eq. 0 .and. iadexp .ge. 2) then
          write(iwr,1021) iadexp, iasexp
 1021     format(/,' input parameters iadexp= ',i2,
     1         ' and iasexp= ',i2,' conflict',/,
     2         ' please change them on input: $drf and $expanx ')
          call hnderr(3,errmsg)
        endif
      endif
c
      defnew = iadexp .ge. 2
      each = iadexp .eq. 5
      asscen = iadexp .lt. 0
      oradial = iradial .eq. 1
c
c-----  initialize
c
      nexp= nat
      nex = 0
      nfit = 0
c
c-----  set memory partitioning
c
      ixo = igmem_alloc(nx)
      ixdx = igmem_alloc(nx)
      ixdy = igmem_alloc(nx)
      ixdz = igmem_alloc(nx)
      ixis = igmem_alloc(nw196(5))
c     ixis = igmem_alloc(10000)
      ixie = igmem_alloc(nx)
      ixfx = igmem_alloc(nx)
c
c-----  determine number of points where potential is to be
c       calculated
c
c     if (iasexp .eq. 2 .and. ifitx .ne. 0) then
c       ixf = ixfx-1
c
c  -----  atomic centra
c
c       if (ifitx .le. 4) then
c         last = last + 3*nat
c         length = last - ixo
c         need = loc10 + length
c         call setc(need)
c         do 100, i = 1, nat
c           nfit = nfit + 1
c           do 200, k = 1, 3
c             ixf = ixf + 1
c             xscm(ixf) = c(k,i)
c 200       continue
c 100     continue
c       endif
c
c  -----  centre of charge
c
c       if (ifitx .eq. 2 .or. ifitx .eq. 3) then
c         last = last + 3
c         length = last - ixo
c         need = loc10 + length
c         call setc(need)
c         nfit = nfit + 1
c         call drfcm(xscm(ixf+1))
c         ixf = ixf + 3
c       endif
c
c  -----  points read from $fitx
c
c       if (ifitx .ge. 3) then
c         rewind(ir)
c  11     read(ir,1100,end=199) text
c1100     format(a80)
c         if (text(:6) .ne. ' $fitx') goto 11
c  21     read(ir,1100) text
c         if (text(:5) .eq. ' $end') goto 31
c         nfit = nfit + 1
c         last = last + 3
c         length = last - ixo
c         need = loc10 + length
c         call setc(need)
c
c         read(text,1300) duma,(xscm(ixf+k+1), k = 1, 3)
c         ixf = ixf + 3
c1300     format(a10,3f20.10)
c         goto 21
c 199     write(iwr,*)' list $fitx not found on input file, ',
c    1             'please check'
c         call hnderr(3,errmsg)
c  31     continue
c       endif
c     endif
c
c-----  determine the number of expansion centra as far as possible
c
c-----  check if centre of charge is required as expansion centre
c
      if (icmexp .eq. 1) then
        nex = nex + 1
        call drfcm(xex(1,nex))
        icch = nexp + 1
      endif
c
c-----  read expansion centra from input if required
c       (data block beginning with $expx)
c       and store them temporarily in xex
c       currently, a maximum of 50 expansion centra to be read in is all
c
c     if (iadexp .eq. 1 .or. iadexp .eq. 2) then
c       rewind(ir)
c  10   read(ir,1100,end=99) text
c       if (text(:6) .ne. ' $xexpc') goto 10
c  20   read(ir,1100) text
c       if (text(:5) .eq. ' $end') goto 30
c       nex = nex + 1
c
c       if (nex .gt. maxex) then
c         write(iwr,1200) maxex
c1200     format(/' number of extra expansion centra on $xexpc ',
c    1             'greater than ',i4,/,
c    2             ' increase maxex near common blocks /xpanx/')
c         call hnderr(3,errmsg)
c       endif
c
c       read(text,1300) duma,(xex(k,nex),k=1,3)
c       goto 20
c  99   write(iwr,*)' list $xexpc not found on input file, ',
c    1               'please check'
c       call hnderr(3,errmsg)
c  30   continue
cxxx  endif
c
cxxx  nexp = nexp + nex
c
cxxx
      if (iasexp .eq. -3) then
c-----
c      add all midpoints to the list of expansion centres
c
        do i = 1 , nat-1
          do j = i+1, nat
            nex = nex + 1
            if (nex .gt. maxex) then
               write(iwr,1250) maxex
 1250     format(/' number of extra expansion centra ',
     1             'greater than ',i4,/,
     2             ' increase maxex near common blocks /xpanx/')
               call hnderr(3,errmsg)
            endif
            do k = 1, 3
              xex(k,nex) = pt5*(c(k,i) + c(k,j))
            enddo
          enddo
        enddo
        iasexp = 1
      endif
c
      nexp = nexp + nex
cxxx
c     if (iasexp .eq.2) then
c       if (ifitx .eq. 0) nfit = nexp
c       ixel = ixfx + 3*nfit
c       ixdel = ixel + 225*nfit
c       ixc = ixdel + 225*nexp
c     else
c       ixel = ixfx + 1
c       ixdel = ixel + 1
c       ixc = ixdel + 1
c     endif
      ixel = igmem_alloc(1)
      ixdel = igmem_alloc(1)
c     ixc = igmem_alloc(3*nexp)
c     lwor = igmem_max_memory()
      lwor = 3*nexp
      nexpmx = lwor/3
      ixc = igmem_alloc(3*nexpmx)
c     last = ixc + 3*nexp
c     length = last - ixo
c     need = loc10 + length
c     call setc(need)
c
c-----  define expansion centra as far as possible
c
      ixex = ixc - 1
c
c-----  atomic centra
c
      do 300, i = 1, nat
        do 400, k= 1, 3
          ixex = ixex + 1
          xscm(ixex) = c(k,i)
  400   continue
  300 continue
c
c-----  extra expansion centra
c
      do 500, i = 1, nex
        do 600, k = 1, 3
          ixex = ixex + 1
          xscm(ixex) = xex(k,i)
  600   continue
  500 continue
c
c-----  define centra at which potential is to be calculated
c
c     if (iasexp .eq. 2 .and. ifitx .eq. 0) then
c
c  -----  copy expansion centra
c
c       do 700, i = 1, 3*nexp
c         xscm(ixfx+i-1) = xscm(ixc+i-1)
c 700   continue
c     endif
c
c-----  assign expansion centra according to requested criteria
c
c-----  expansion centra relative to center of nuclear charges
c
      call drfcm(cm)
      ixex = ixc - 1
      do 800 n=1,nexp
        do 900 k=1,3
          ixex = ixex + 1
          xscm(ixex)=xscm(ixex)-cm(k)
  900   continue
  800 continue
      if (out)  call hatout(xexp,3,nexp,2,'xexpass')
c
c----- read in transformation matrices for s,p,d,f basis functions.
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      call rdedx(dtr,nw196(2),ibl196(2),idaf)
      call rdedx(ftr,nw196(3),ibl196(3),idaf)
c
c----- read in symmetry array - iso
c
      call rdedx(xscm(ixis),nw196(5),ibl196(5),idaf)
c
      if (oassign) then
        call asgnrf(nfit,nexp,icch,iasexp,defnew,each,asscen,
     1     oradial,nshell,
     2     xscm(ixo),xscm(ixdx),xscm(ixdy),xscm(ixdz),
     1     xscm(ixfx),xscm(ixel),xscm(ixdel),
     2     xscm(ixc),xscm(ixie),xscm(ixis))
c
c
c-----  write expansion centra and assignation to dafile
c
      call dawrit(idafdrf,iodadrf,xscm(ixie),nx,2,navdrf)
        ibeen = ibeen + 1
      else
        if (each) then
          write(iwr,*) 'Cannot use assign .... in a geometry ',
     1 'optimization'
          call hnderr(3,errmsg)
        endif
      endif
      call dawrit(idafdrf,iodadrf,xscm(ixc),3*nexp,1,navdrf)
      if (out)  call imatou(xscm(ixie),num,num,3,'iexp')
c
      call gmem_free(ixc)
      call gmem_free(ixdel)
      call gmem_free(ixel)
      call gmem_free(ixfx)
      call gmem_free(ixie)
      call gmem_free(ixis)
      call gmem_free(ixdz)
      call gmem_free(ixdy)
      call gmem_free(ixdx)
      call gmem_free(ixo)
c
c-----  set some dimensions
c
      nexp4 = 4*nexp
      nwtc = nexp4 + ngran + 1
c
      nzfp = nexp4 + 1
      nzfn = nexp4 + ngran + 1
c
      lomga = nexp4 + ngran + 1
      nomga = lomga**2
c
c-----  reset core memory
c
c     call setc(loadcm)
      return
      end
      subroutine drfzfp(zfp,vrp,eclas,eelst,edisp,erep,
     1                  egrcls,edum,ngrpair)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c          routine forms vectors and scalars describing the inter-    c
c          actions between the external points.                       c
c          since the fields and potentials of the point charges       c
c          will be affected by polarization effects, they are stored  c
c          as the last column of the matrix "wt" which will be        c
c          contracted with 'a' (the polarizability)                   c
c                                                                     c
c       clas is the electrostatic energy between the points           c
c                                                                     c
c                                                                     c
c      interactions between points belonging to the same group        c
c      (i.e. having the same group names) are excluded: this is       c
c      essentially a quantum mechanical interaction, and should be    c
c      treated as such                                                c
c                                                                     c
c      interactions between "overlapping" distributions, as           c
c      defined by their polarizabilities, should be reduced           c
c      consistently with the reduction in routine drftpq              c
c                                                                     c
c     routines called: distab,clear                                   c
c                                                                     c
c                                                                     c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * c
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arguments
c
      dimension zfp(ndim,ngran), vrp(ndim,ngran)
      dimension eelst(ngrpair), edisp(ngrpair), erep(ngrpair),
     1          egrcls(ngrpair), edum(5*ngrpair)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/bas)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/grinf)
INCLUDE(comdrf/drfamb)
c
INCLUDE(comdrf/drfbem)
INCLUDE(comdrf/rad)
c
c-----  local variables
c
      dimension zfpi(3,mxgran)
      dimension pq(3)
      dimension b(3,3), b2(3,3), polten(3,3)
      character*16 namei,namej
      logical poli,polj
c
      data zero,one,two,three,four/0.0d00,1.0d00,2.d0,3.d0,4.d0/
      data pt5,sixth /0.50d00,.1666666666666666666667d00/
      data pt3, pt75, onept5 /
     +     0.333333333333333333333333d0,7.5d-01,1.5d00/
c
      data repcut /1.0d02/
      data small /1.0d-03/
c
c-----  begin
c
      if (ngran.gt.3000) call caserr( 
     + 'drfsub.f - drfzfp: ngran>3000 increase dimension of zfpi >3000')
      do 50, i = 1, ngrpair
c 1-----
        eelst(i) = zero
        edisp(i) = zero
        erep(i) = zero
        egrcls(i) = zero
c 1-----
   50 continue
c
c-----  count polarizable external points ii
c
      imp = 1
c
c-----  loop over external points
c
      do 100, ii = 1, nxtpts
c 1-----
        zi = chrg(ii)
        alfi = alfext(ii)
        if (alfi .eq. zero) alfi = one
        namei = nxcent(ii)
        ipol = mpol(imp)
        poli = ipol .eq. ii
        igrani = igranl(ii)
c
        if (poli) then
c   2-----
          call clear(zfpi,3*ngran)
          irow = (imp-1)*3
c   2-----
        endif
c
c  -----  count polarizable points jj
c
        jmp = 1
c
c  -----  loop over external points
c
        do 90, jj = 1, nxtpts
c   2-----
          jpol = mpol(jmp)
          polj = jpol .eq. jj
          namej = nxcent(jj)
          igranj = igranl(jj)
c
          if (jj .eq. ii) goto 40
c
c    -----  exclude interactions between members of the same group
c
          if (namei(7:(6+ngrnam)) .eq. namej(7:(6+ngrnam))) goto 40
c
c    -----  set index for analysis of interactions
c
          indxen = ia(max(igrani,igranj)) + min(igrani,igranj)
c
          zj = chrg(jj)
          alfj = alfext(jj)
          if (alfj .eq. zero) alfj = one
c
c    -----  distance vector between ii and jj in -pq-
c
          call distab(xpts(1,jj),xpts(1,ii),pq,dist)
          if (dist .lt. small) then
            write(iwr,95) namei, namej, dist
  95        format(/,1x,'WARNING: non-excluded interaction ',
     1   'between ', a16, ' and ', a16, ' at', e15.8,
     2   ' bohr distance: skipped')
            goto 40
          endif
          dmind1 = one/dist
          v = one
          factp = one
          factf = one
          factd = one
c
c    -----  scale interaction between point charges according to thole
c           (optional)
c
          if (modxza .ne. 0) then
c     3-----
            s = (alfi*alfj)**sixth
            v = dist/s
c
            if (ithole .eq. 1) then
c       4-----
              if (v .le. afact) then
                av = v/afact
                factp = av**4 - two*av**3 + two*av
                factf = four*av**3 - three*av**4
                fact5 = 3.0d0*av**4
              factd = (fact5**2 - 2.0d0*fact5*factf + 3.0d0*factf**2) 
     +                       * sixth
              endif
c       4-----
            else
c       4-----
              au = afact*v
              factp = (one - (pt5*au + one)*exp(-au))
              factf = (one - (pt5*au**2 + au + one)*exp(-au))
              term = (pt5*au**3 + onept5*au**2 + three*au + three)
              fact5 = three - term*exp(-au)
              factd = (fact5**2 - 2.0d0*fact5*factf + 3.0d0*factf**2) 
     +                       * sixth
c       4-----
            endif
c     3-----
          endif
c
          if (poli) then
c     3-----
c      -----  add field of charge jj at polarizability imp
c             to total field of charges at polarizability imp
c             sorted w.r.t. the classical energy analysis groups
c
            dmin3 = dmind1**3
c
            do 20, k = 1, 3
              zfpi(k,igranj) = zfpi(k,igranj) + zj*dmin3*pq(k)*factf
  20        continue
c     3-----
          endif
c
          if (iclintd .eq. 1) then
c     3-----
            if (poli .and. polj) then
c       4-----
c        -----  calculate approximate dispersion between
c               classical entities via slater-kirkwood formula
c
              if ((igrppol .eq. 0) .and.
     1           (namei(1:5) .eq. 'group') .and.
     2           (namej(1:5) .eq. 'group')) goto 25
c
              fsk = sqrt(polar(imp)/vale(imp)) +
     1              sqrt(polar(jmp)/vale(jmp))
c
c        -----  isotropic or non-isotropic dispersion may be used
c
              dmin3 = dmind1**3
c
              if (isodis .eq. 1) then
c         5-----
c          -----  isotropic polarisabilities
c
                dmin6 = (dmin3)**two
                discon = onept5*polar(imp)*polar(jmp)*dmin6*factd
c         5-----
              else
c         5-----
c          -----  non-isotropic polarisabilities
c          -----  calculate interaction tensor bq
c
                call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c          -----  multiply polarisability tensors: bq = aj bq ai
c
                call tenmul(b,b,b2,3)
                if (namei(1:5) .eq. 'group') then
                  call hexpand(grpol(1,imp),polten,3,3)
                else
                  call mkdiagm(polar(imp),polten,3)
                endif
                call tenmul(b2,polten,b,3)
                if (namej(1:5) .eq. 'group') then
                  call hexpand(grpol(1,jmp),polten,3,3)
                else
                  call mkdiagm(polar(jmp),polten,3)
                endif
                call tenmul(polten,b,b2,3)
c
c          -----  calculate trace of product tensor bq
c
                call trace(b2,discon,3)
                discon = discon / four
c         5-----
              endif
c
              edisp(indxen) = edisp(indxen) - discon / fsk
  25          continue
c       4-----
            else if (igrppol .eq. 0) then
c       4-----
c           use atom polarizabilities
c
              dmin3 = dmind1**3
c
              if ((atpol(ii) .ne. zero) 
     1      .and. (atpol(jj) .ne. zero)) then
c       4-----
                fsk = sqrt(atpol(ii)/valat(ii)) +
     1                sqrt(atpol(jj)/valat(jj))
c
c        -----  isotropic or non-isotropic dispersion may be used
c
                if (isodis .eq. 1) then
c         5-----
c          -----  isotropic polarisabilities
c
                  dmin6 = (dmin3)**2
                  discon = onept5*atpol(ii)*atpol(jj)*dmin6*factd
c         5-----
                else
c         5-----
c          -----  non-isotropic polarisabilities
c          -----  calculate interaction tensor bq
c
                  call drftpq(pq,dmind1,dmin3,ithole,afact,v,b)
c
c          -----  multiply polarisability tensors: bq = aj bq ai
c
                  call tenmul(b,b,b2,3)
                  call hexpand(atpolt(1,ii),polten,3,3)
                  call tenmul(b2,polten,b,3)
                  call hexpand(atpolt(1,jj),polten,3,3)
                  call tenmul(polten,b,b2,3)
c
c          -----  calculate trace of product tensor bq
c
                  call trace(b2,discon,3)
                  discon = discon / four
c         5-----
                endif
c
                edisp(indxen) = edisp(indxen) - discon / fsk
c       4-----
              endif
c     3-----
            endif
          endif
c
c    -----  calculate electrostatic interaction
c           between external charge groups
c
          eelst(indxen) = eelst(indxen) + zi*zj*dmind1*factp
c
          if ((iclintr .eq. 1) .and. (dist .le. repcut)) then
c     3-----
c      -----  calculate approximate repulsion
c             between external charge groups, with empirical r12-term
c             from charmm
c
            if ((namei(:2) .ne. 'qq') .and.
c    1          (namei(:2) .ne. 'xx') .and.
     2          (namei(:2) .ne. 'e ') .and.
     3          (namei(:2) .ne. ' ') .and.
     4          (namei(:2) .ne. 'gr') .and.
     5          (namej(:2) .ne. 'qq') .and.
c    6          (namej(:2) .ne. 'xx') .and.
     7          (namej(:2) .ne. 'e ') .and.
     8          (namej(:2) .ne. ' ') .and.
     9          (namej(:2) .ne. 'gr')) then
c       4-----
              alfi = alfext(ii)
              call drfnval(namei,nvali,znuc)
              aoverni = sqrt(alfi/nvali)
c
              alfj = alfext(jj)
              call drfnval(namej,nvalj,znuc)
              aovernj = sqrt(alfj/nvalj)
c
              ri = radext(ii)*rfact
              rj = radext(jj)*rfact
c
c          -----  account for h-bonding
c
              if ((ihbond .eq. 1) 
     1            .and. (dist .le. hbondl)) then
c         5-----
                if(namej(:2).eq.'h'.and.
     1            (namei(:2).eq.'n'.or.
     2             namei(:2).eq.'o'.or.
     3             namei(:2).eq.'f')) rj=hbondr
                if(namei(:2).eq.'h'.and.
     1            (namej(:2).eq.'n'.or.
     2             namej(:2).eq.'o'.or.
     3             namej(:2).eq.'f')) ri=hbondr
c         5-----
              endif
c
c        -----  calculate model repulsion
c
              fac1 = aoverni + aovernj
              fac2 = pt75*alfi*alfj
              fac3 = (ri+rj)**6
              fac4 = dmind1**12
              erep(indxen) = erep(indxen) + fac2*fac3*fac4/fac1
c       4-----
            endif
c     3-----
          endif
c
   40     if(polj) jmp=jmp+1
c   2-----
   90   continue
c
c  -----  store zfpi in zfp
c
        if (poli) then
c   2-----
c    -----  field at polarisability imp due to analysis groups
c
          do 60, l = 1, ngran
            do 60 k = 1, 3
              zfp(irow+k,l) = zfp(irow+k,l) + zfpi(k,l)
              vrp(irow+k,l) = vrp(irow+k,l) - zfpi(k,l)
   60     continue
c
          imp=imp+1
c   2-----
        endif
c 1-----
  100 continue
c
c-----  correct energies for double counting and calculate
c       total classical energy
c
      eclas = zero
      do 200, i = 1, ngrpair
c 1-----
        eelst(i) = pt5*eelst(i)
        if (iclinte .eq. 0) eelst(i) = zero
        edisp(i) = pt5*edisp(i)
        erep(i) = pt5*erep(i)
c
        egrcls(i) = eelst(i) + edisp(i) + erep(i)
c
        eclas = eclas + egrcls(i)
c 1-----
  200 continue
c
      return
      end
      subroutine drfzfbe(ieps,ineq,zfbe,vrbe,xsurf,xnorm,area)
c------
c             calculates potentials and fields on boundary from
c             classical point charges
c
c     ndima: dimension of polarization matrix (if present)
c     ndimb: dimension of boundary matrix.
c            =   nbem if kappa.eq.0.0
c            = 2*nbem if kappa.ne.0.0
c     ndim: dimension of complete relay matrix
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension zfbe(ndim,ngran), vrbe(ndim,ngran)
      dimension xsurf(3,nbem), xnorm(3,nbem)
      dimension area(nbem)
c
c-----  common blocks
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
INCLUDE(comdrf/grinf)
INCLUDE(comdrf/neqpar)
c
INCLUDE(comdrf/drfbem)
c
c-----  local arrays
c
      dimension pi(3)
c
c-----  local variables
c
      logical kapnoz
c
      data zero,one,two,four /0.0d00,1.0d00,2.0d00,4.0d00/
c
c-----  begin
c
c-----  set some variables
c
      if (ineq .eq. 0) then
c 1-----
        if (ieps .eq. 0) then
c   2-----
          epsilon = eps1
          kappa = kappa1
        else
          epsilon = eps2
          kappa = kappa2
c   2-----
        endif
c 1-----
      else
c 1-----
        if (ieps .eq. 0) then
c   2-----
          epsilon = epsneq1
          kappa = kapneq1
        else
          epsilon = epsneq2
          kappa = kapneq2
c   2-----
        endif
c 1-----
      endif
c
      pie=four*atan(one)
      epsfact= one/(two*pie*(one+epsilon))
      kapnoz = kappa .ne. zero
      expkd = one
      fkapone = one
c
c-----  loop over boundary elements -ni-
c
      do 2000, ni = 1, nbem
c 1-----
c  -----  loop over external charges -np-
c
        do 1000, np = 1, nxtpts
c   2-----
c    -----  get index of classical analysis group that -np- belongs to
c
          igranp = igranl(np)
c
c    -----  pi:  vector from -np- to -ni- = (i-p)
c           dist:  length of pi
c
          call distab(xpts(1,np),xsurf(1,ni),pi,dist)
          dmind1 = one/dist
          dmin3 = dmind1**2*dmind1
c
          if (kapnoz) then
c     3-----
            expkd = exp(-(kappa*dist))
            fkapone = one + (kappa*dist)
c     3-----
          endif
c
c    -----  fpini:  field of unit (positive) charge
c                   in -np- at -ni-, contracted with
c                   normal vector in -ni- = (i-p).n(i)/dist**3
c                   = f(p;i) . n(i)
c
          fpini = dmin3*ddot(3,pi,1,xnorm(1,ni),1)
c         fpini = dmin3*adotb(pi,xnorm(1,ni),3)
c
c    -----  potential of source charge in -np- at -ni-
c           scaled with 1/(2pi(1+eps)), as input for coupling
c           equations for w(i)
c
          zfbe(ndima+ni,igranp) = zfbe(ndima+ni,igranp)
     1                            + epsfact*chrg(np)*dmind1
c
c    -----  reaction potential energy operator
c           for source charges, containing
c           contribution of unit surface dipoles w(i),
c           already multiplied by source charges
c
c         sum(p) k(i;p)s(i)*q(p)
c
c           the energy is evaluated in subr drfomga,
c           by contracting with the induced dipole density array
c
c    -----  the minus sign is a result of the use of the
c           inverted field: needed is f(i;p), fpini=f(p;i) . n(i)
c
          vrbe(ndima+ni,igranp) = vrbe(ndima+ni,igranp) -
     1         (epsilon*fkapone*expkd - one)*fpini*area(ni)*chrg(np)
c
c    -----  if non-zero ionic strength, a set of equations
c           for z(i) is added
c
          if (kapnoz) then
c     3-----
c      -----  minus field of source charge in -np- at -ni-,
c             contracted with normal vector at -ni-,
c             scaled with eps/(2pi(1+eps)), as input for
c             coupling equations for z(i)
c
            zfbe(ndima+nbem+ni,igranp) = zfbe(ndima+nbem+ni,igranp) -
     1             epsilon*epsfact*chrg(np)*fpini
c
c      -----  reaction potential energy operator
c             for source charges, containing
c             contribution of unit surface charges z(i),
c             already multiplied by source charges
c
c         sum(p) l(i;p)s(i)*q(p)
c
c             the energy is evaluated in subr drfomga,
c             by contacting with the induced charge density array
c
            vrbe(ndima+nbem+ni,igranp) = vrbe(ndima+nbem+ni,igranp) +
     1            (one - expkd)*dmind1*area(ni)*chrg(np)
c     3-----
          endif
c   2-----  next external charge
 1000   continue
c 1-----  next boundary element
 2000 continue
      return
      end
      subroutine drfone(ieps,s,dipx,dipy,dipz,rxx,ryy,rzz,rxy,rxz,ryz,
     1                  omega,xexp,iexpc,velq,vsnue,vseln,
     2                  vself,vext,relay,zfp,aind,vr,ijbit,indx)
c------
c
c       this routine drives the calculation of one-electron
c       integrals of the reaction field.
c
c     vext = static external potential
c     vext = also: screening of electron/external charge
c            interaction due to externally induced dipoles
c
c     velq = screening of electron/external charge
c            interaction due to electronic induced dipoles
c
c     vsnue = screening of  nuclear/electron interaction due to
c             nuclear induced dipoles
c     vseln = screening of  nuclear/electron interaction due to
c             electronic induced dipoles
c
c     vself= interaction of charge distribution with its
c            own reaction field
c
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gvalue)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension s(nchd), dipx(nchd), dipy(nchd), dipz(nchd)
      dimension rxx(nchd), ryy(nchd), rzz(nchd),
     1          rxy(nchd), rxz(nchd), ryz(nchd)
      dimension omega(nwtc,nwtc), relay(ndim,ndim), vr(nwtr,nwtc)
      dimension vext(nchd,ngran),velq(nchd,ngran)
      dimension vsnue(nchd),vseln(nchd),vself(nchd)
      dimension iexpc(nchd), ijbit(nchd), indx(ndim)
      dimension zfp(ndim,ngran)
      dimension aind(ndim)
      dimension xexp(3,nexp)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/sym)
INCLUDE(comdrf/shlint)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(../m4/common/drfopt)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/rfene)
c
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/runpar)
c
c-----  local arrays
c
      dimension zfij(4)
      dimension cdum(3)
c
      logical block,blocks,blocki
      character*8 block1, block2
c
      logical oarf
c
c-----  data statements
c
      data zero,one,pt5,three/0.0d00,1.0d00,0.5d00,3.0d00/
c      data cutoff /1.0d-09/
      data block1 /'blocks  '/
      data block2 /'blocki  '/
c
c     comdrf/shlint
c
      cutoff = 1.0d-09
c
c-----  begin
c
c
      noprt = 0
c
      blocks=blktyp.eq.block1.and.noprt.eq.0
      blocki=blktyp.eq.block2.and.noprt.eq.0
      block =blocks.or.blocki
c
      oarf = iarfcal .eq. 1
c
      ilself = 58
c
      if (iextdip .eq. 1) then
        fact = one
      else
        fact = pt5
      endif
c
      if (ieps .eq. 0) then
c 1-----
        ilomg = 50
        ilelq = 52
        ilsne = 15
        ilexd = 17
        ilsen = 54
        ilij = 56
        ilr = 43
        ilindx = 44
        ilzf = 22
        ilvr = 28
c 1-----
      else
c 1-----
        ilomg = 51
        ilelq = 53
        ilsne = 16
        ilexd = 18
        ilsen = 55
        ilij = 57
        ilr = 47
        ilindx = 48
        ilzf = 24
        ilvr = 32
c 1-----
      endif
c
      if (mcupdt) then
        ilvr = ilvr + 1
      endif
c
c-----  read necessary integrals from da10
c
      call daread(idafh,ioda,s   ,nchd,12)
      call daread(idafh,ioda,dipx,nchd,53)
      call daread(idafh,ioda,dipy,nchd,54)
      call daread(idafh,ioda,dipz,nchd,55)
c
      call daread(idafdrf,iodadrf,xexp ,3*nexp,1)
c     nchd2 = lenint(nchd)
c     call daread(idafdrf,iodadrf,iexpc,nchd2,2)
      call daread(idafdrf,iodadrf,iexpc,nchd,2)
      if (idrfout .eq. 3) call imatou(iexpc,num,num,3,'iexpc')
c
      if (field(5:) .eq. 'scf') then
        call daread(idafdrf,iodadrf,omega,nomga,ilomg)
      endif
c
      if (idrfout .eq. 3) then
         call hatout(s   ,num,num,3,'s'   )
         call hatout(dipx,num,num,3,'dipx')
         call hatout(dipy,num,num,3,'dipy')
         call hatout(dipz,num,num,3,'dipz')
      endif
c
      call clear(vext,ngran*nchd)
c
      if (ieps .eq. 0) then
        if (field(:4) .ne. ' ') then
c 1-----
c  -----  calculate interaction integrals with external charges
c
        if (iexpza .eq. 0) then
c   2-----
          i40 = nchd + 1
          i50 = i40 + 1
c         if(block ) i50 = i40 + l2sym
          i60 = i50 + 1
c         if(blocki) i60 = i50 + l2sym
          call drfexc(vext,cdum(1),cdum(2),cdum(3))
c   2-----
        endif
c
        if (iexpza .ne. 0)
     1    call drfxza(nchd,ngran,vext,s,dipx,dipy,dipz,iexpc)
c
c  -----  write interaction integrals with external charges
c         record -12- of da31
c
        endif
        call dawrit(idafdrf,iodadrf,vext,ngran*nchd,12,navdrf)
c 1-----
      endif
c
      if (idrfout .eq. 3) then
        do 50, igr = 1, ngran
          call hatout(vext(1,igr),num,num,3,'vextc')
   50   continue
      endif
c
      if ((nodiscr .eq. 0) .and. (iextdip .eq. 1)) then
c 1-----
c  -----  field due to dipoles induced by the external charges
c         at the expansion centra, to be added to the
c         static external potential
c
        call clear(vext,ngran*nchd)
c
        do 100, igr = 1, ngran
c   2-----
          call clear(aind,ndim)
c
          if (itsolv .eq. 0) then
c     3-----
            call daread(idafdrf,iodadrf,relay,ndim*ndim,ilr)
c           ndim2 = lenint(ndim)
c           call daread(idafdrf,iodadrf,indx,ndim2,ilindx)
            call daread(idafdrf,iodadrf,indx,ndim,ilindx)
c
c      -----  zfp at vr matrix
c
            call daread(idafdrf,iodadrf,zfp,ngran*ndim,ilzf)
c
c      -----  solve relay-equations for external source field
c
            call luelmf(relay,zfp(1,igr),indx,ndim,ndim,aind)
c
c      -----  zfp is overwritten by vr
c
            call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilvr)
c
c      -----  calculate expanded potential due to externally induced
c             dipoles at the expansion centra, contracted with
c             charge distribution's overlap and first moment integrals
c
            call extdcal(npol3,ndimb,nexp,num,fact,vr,vext(1,igr),
     1                   dipx,dipy,dipz,s,iexpc,aind)
c
            if (idrfout .eq. 3)
     1         call hatout(vext(1,igr),num,num,3,'vextd')
c
            extnucc = fact*ddot(npol3+ndimb,aind,1,vr(1,nwtc),1)
c           extnucc = fact*adotb(aind,vr(1,nwtc),npol3+ndimb)
c
            if (ieps .eq. 0) then
              extnucg(igr,iactst) = extnucg(igr,iactst) + extnucc
              extnuc(iactst) = extnuc(iactst) + extnucc
              sxtnucg(igr,iactst) = zero
            else
              extnucg(igr,iactst) = extnucg(igr,iactst) - extnucc
              extnuc(iactst) = extnuc(iactst) - extnucc
              sextnog(igr,iactst) = zero
            endif
c     3-----
          endif
c
c    -----  add/subtract field due to externally induced dipoles to stat
c
          call daread(idafdrf,iodadrf,velq,ngran*nchd,12)
c
          if (ieps .eq. 0) then
c     3-----
            do 90, i = 1, nchd
              velq(i,igr) = velq(i,igr) + vext(i,igr)
   90       continue
          else
            do 80, i = 1, nchd
              velq(i,igr) = velq(i,igr) - vext(i,igr)
   80       continue
c    3-----
          endif
c
          call dawrit(idafdrf,iodadrf,velq,ngran*nchd,12,navdrf)
c   2-----
  100   continue
c
        sextnuc(iactst) = zero
        sextno(iactst) = zero
c 1-----
      endif
c
c-----  further reaction field scf?
c
      if (field(5:) .ne. 'scf') return
c
      if (oarf) then
c
c  -----  collect fixed contributions to inducing 
c         field for ARF coupling 
c
        call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilwt)
        call clear(zfp,ndim)
        do i = 1, ngran
          do j = 1, ndim
            zfp(j,1) = zfp(j,1) + vr(j,4*nexp+i)
          enddo
        enddo
        do j = 1, ndim
          zfp(j,1) = zfp(j,1) + vr(j,4*nexp+ngran+1)
        enddo
c
        call dawrit(idafdrf,iodadrf,zfp,ndim,101,navdrf)
        call daread(idafdrf,iodadrf,vr,nwtc*nwtr,ilvr)
c
      endif
c
c-----  get nuclear interactions
c
c
      if (ieps .eq. 0) then
c 1-----
        snucnuc(iactst) = omega(nwtc,nwtc)
c
        snucext(iactst) = zero
        sextnuc(iactst) = zero
c
        do 200, igr = 1, ngran
c   2-----
          snucxt = fact*omega(nwtc,(nwtc-ngran-1+igr))
          snucxtg(igr,iactst) =  snucxt
          snucext(iactst) = snucext(iactst) + snucxt
c
          sxtnuc = fact*omega((nwtc-ngran-1+igr),nwtc)
          sxtnucg(igr,iactst) = sxtnuc
          sextnuc(iactst) = sextnuc(iactst) + sxtnuc
c   2-----
  200   continue
c 1-----
      else
c 1-----
        snucno(iactst) = omega(nwtc,nwtc)
c
        snucexo(iactst) = zero
        sextno(iactst) = zero
c
        do 300, igr = 1, ngran
c   2-----
          snucxt = fact*omega(nwtc,(nwtc-ngran-1+igr))
          sncexog(igr,iactst) = snucxt
          snucexo(iactst) = snucexo(iactst) + snucxt
c
          sxtnuc = fact*omega((nwtc-ngran-1+igr),nwtc)
          sextnog(igr,iactst) = sxtnuc
          sextno(iactst) = sextno(iactst) + sxtnuc
c   2-----
  300   continue
c 1-----
      endif
c
c-----  interaction with induced dipoles and
c       induced boundary charges/dipoles due to the external charges
c
      call clear(vext,ngran*nchd)
      call clear(velq,ngran*nchd)
      call clear(vsnue,nchd)
      call clear(vseln,nchd)
c
c-----  set up loop over charge distributions
c
      if (idrfout .eq. 4) write(iwr,*) ' Testing ij/ij RF interractions'
      ij = 0
      do 400, ii = 1, num
c 1-----
        do 390, jj = 1, ii
c   2-----
          ij = ij + 1
          ijbit(ij) = 1
          ijexp = iexpc(ij)
c
          p(1) = dipx(ij)
          p(2) = dipy(ij)
          p(3) = dipz(ij)
          p(4) = s(ij)
c
c    -----  test (ij/ij) value and set ijbit accordingly
c
          call drfoab(ijexp,ijexp,nwtc,omega)
          call matvec(omgab,p,zfij,4,.false.)
          try = abs(ddot(4,zfij,1,p,1))
c         try = abs(adotb(zfij,p,4))
c
c    -----  if integrals of charge distribution -ij-
c           with itself are very small, all other contributions
c           of this distribution will be skipped
c           (also later in drfself and drfdisp2)
c
c         if (try .lt. cutoff) ijbit(ij) = 0
          if (idrfout .eq. 4) write(iwr,*) ii,jj,ij,try,ijbit(ij)
c
c    -----  skip integrals for which (ij/ij) < cutoff
c
          if (ijbit(ij) .eq. 0) goto 390
c
          ip = (ijexp -1)*4
c
c    -----  interaction between charge distribution and
c           dipoles that are induced by the external charges
c
          do 380, igr = 1, ngran
c     3-----
            extpol = zero
            polext = zero
c
            do 370, k = 1, 4
c       4-----
              extpol = extpol + omega(nzfp+igr-1,ip+k)*p(k)
              polext = polext + omega(ip+k,nzfp+igr-1)*p(k)
c       4-----
  370       continue
c
c      -----  note the minus sign, due to the fact that this
c             potential is coupled back to the electrons
c
c      -----  note: factor 0.5 because of double counting
c
            vext(ij,igr) = vext(ij,igr) - fact*extpol
            velq(ij,igr) = velq(ij,igr) - fact*polext
c     3-----
  380     continue
c
c    -----  screening of the nuclear attraction through
c           polarized environment
c
          snue = zero
          seln = zero
c
c    -----  loop over all nuclei
c
          do 360, n = 1, nat
c     3-----
            q(1) = c(1,n)
            q(2) = c(2,n)
            q(3) = c(3,n)
            q(4) = one
c
            call drfoab(ijexp,n,nwtc,omega)
            call matvec(omgab,q,zfij,4,.false.)
            snue = snue + ddot(4,zfij,1,p,1)*czan(n)
c           snue = snue + adotb(zfij,p,4)*czan(n)
c
            call drfoab(n,ijexp,nwtc,omega)
            call matvec(omgab,p,zfij,4,.false.)
            seln = seln + ddot(4,zfij,1,q,1)*czan(n)
c           seln = seln + adotb(zfij,q,4)*czan(n)
c     3-----
  360     continue
c
c    -----  note the minus sign, due to the fact that this
c           potential is coupled back or due to the electrons
c
          vsnue(ij) = vsnue(ij) - snue
          vseln(ij) = vseln(ij) - seln
c   2-----  next jj
  390   continue
c 1-----  next ii
  400 continue
c
      if(idrfout .eq. 3) call imatou(ijbit,num,num,3,'ijbit')
      if(idrfout .eq. 3) then
        do 500, igr = 1, ngran
          call hatout(vext(1,igr),num,num,3,'vextd')
          call hatout(velq(1,igr),num,num,3,'velq')
  500   continue
      endif
      if(idrfout .eq. 3) call hatout(vsnue,num,num,3,'vsnue')
      if(idrfout .eq. 3) call hatout(vseln,num,num,3,'vseln')
c
      call dawrit(idafdrf,iodadrf,vext ,ngran*nchd,ilexd,navdrf)
      call dawrit(idafdrf,iodadrf,velq ,ngran*nchd,ilelq,navdrf)
      call dawrit(idafdrf,iodadrf,vsnue,nchd,ilsne,navdrf)
c     nchd2 = lenint(nchd)
c     call dawrit(idafdrf,iodadrf,ijbit,nchd2,ilij,navdrf)
      call dawrit(idafdrf,iodadrf,ijbit,nchd,ilij,navdrf)
      call dawrit(idafdrf,iodadrf,vseln,nchd,ilsen,navdrf)
c
c-----  calculate the second moment integrals for the expansion
c       of the field
c
      cdum(1) = gx
      cdum(2) = gy
      cdum(3) = gz
      gx = zero
      gy = zero
      gz = zero
c
      call qmints(rxx,ryy,rzz,rxy,rxz,ryz)
      gx = cdum(1) 
      gy = cdum(2) 
      gz = cdum(3) 
c     call qdpint(rxx,ryy,rzz,rxy,rxz,ryz,cdum,num)
c
c-----  save second moment integrals on da31
c
      call dawrit(idafdrf,iodadrf,rxx,nchd,61,navdrf)
      call dawrit(idafdrf,iodadrf,ryy,nchd,62,navdrf)
      call dawrit(idafdrf,iodadrf,rzz,nchd,63,navdrf)
      call dawrit(idafdrf,iodadrf,rxy,nchd,64,navdrf)
      call dawrit(idafdrf,iodadrf,rxz,nchd,65,navdrf)
      call dawrit(idafdrf,iodadrf,ryz,nchd,66,navdrf)
c
c-----  compute the self-energy integrals
c       only if required (dispersion estimate)
c
      call clear(vself,nchd)
c
      if (gamdrf .ne. zero) then
c 1-----
        if ((itwoeps .eq. 1) .and. (ieps .eq. 0)) goto 600
c
        call drfself(s,dipx,dipy,dipz,ijbit,
     1            rxx,ryy,rzz,rxy,rxz,ryz,omega,iexpc,vself)
c
c  -----  save self on record -ilself- of da31
c
        call dawrit(idafdrf,iodadrf,vself,nchd,ilself,navdrf)
c
c  -----  calculate "reverse" dispersion (optional)
c
        call clear(vself,nchd)
        call clear(vext,nchd)
c
        if (irevdis .eq. 1)
     1     call drfdisp2(xexp,iexpc,ijbit,vself,vext)
c
c  -----  save disp2 on record -13- of da31
c
        call dawrit(idafdrf,iodadrf,vself,nchd,13,navdrf)
c
  600   continue
c 1-----
      endif
c
      return
      end
      subroutine arfone(fone,dens,
     1   s,dx,dy,dz,
     2   relay,indx,wt,iexpc,fsf,indmom,sf,enadd)
c------
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfpar)
INCLUDE(../m4/common/infoa)
c
      dimension fone(nchd), dens(nchd)
      dimension s(nchd), dx(nchd), dy(nchd), dz(nchd)
      dimension relay(ndim,ndim)
      dimension indx(ndim)
      dimension wt(nwtr,nwtc)
      dimension iexpc(nchd)
      dimension fsf(ndim)
      dimension sf(ndim)
      REAL indmom
      dimension indmom(ndim)
c
      dimension p(4)
c
      data pt5 /0.5d00/
      data dzero, one, two /0.0d00, 1.0d00, 2.0d00/
c
c-----  BEGIN
c
      enadd = dzero
c
      call daread(idafdrf,iodadrf,wt,nwtc*nwtr,26)
c
      do l = 1, ndim
        sf(l) = fsf(l)
      enddo
c
c-----  effective loop over charge distributions
c
      ij = 0
      do i = 1, num
        do j = 1, i
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
          p(1) = fact*dens(ij)*dx(ij)
          p(2) = fact*dens(ij)*dy(ij)
          p(3) = fact*dens(ij)*dz(ij)
          p(4) = fact*dens(ij)*s(ij)
c
c    -----  execute expansion
c
          do l = 1, ndim
            do k = 1, 4
              sf(l) = sf(l) - p(k)*wt(l,ip+k)
            enddo   
          enddo   
c
        enddo
      enddo
c
      call luelmf(relay,sf,indx,ndim,ndim,indmom)     
      call clear(fone,nchd)
      call daread(idafdrf,iodadrf,wt,nwtc*nwtr,28)
c
c-----  calculate interaction with nuclei and external 
c       charges and add to the total energy
c
      snucen = ddot(ndim,wt(1,nexp*4+ngran+1),1,indmom,1)
c     snucen = adotb(wt(1,nexp*4+ngran+1),indmom,ndim)
      sexten = dzero
      do i = 1, ngran
        sexten = sexten + ddot(ndim,wt(1,nexp*4+i),1,indmom,1)
c       sexten = sexten + adotb(wt(1,nexp*4+i),indmom,ndim)
      enddo
      enadd = snucen + sexten
c
      ij = 0
      do i = 1, num
        do j = 1, i
          ij = ij + 1
c
c         fact = two
c         if (j .eq. i) fact = one
          fact = one
c         if (j .eq. i) fact=two
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
          p(1) = fact*dx(ij)
          p(2) = fact*dy(ij)
          p(3) = fact*dz(ij)
          p(4) = fact*s(ij)
c
c    -----  execute expansion
c
          do l = 1, ndim
            do k = 1, 4
              fone(ij) = fone(ij) 
     1    - p(k)*wt(l,ip+k)*indmom(l)
            enddo   
          enddo
        enddo
      enddo
c
      return
      end
      subroutine drfxza(nchd,ngran,vext,s,dipx,dipy,dipz,iexpc)
c------
c       calculates integrals for expanded external static potential
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension vext(nchd,ngran)
      dimension s(nchd),dipx(nchd),dipy(nchd),dipz(nchd),iexpc(nchd)
c
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
c
INCLUDE(comdrf/drfzfa)
c
      dimension p(4)
c
      data zero, threshs /0.0d00,1.0d-10/
c
      do 100, i = 1, num
c 1-----
        do 90, j = 1, i
c   2-----
          ij = ia(i) + j
          jf = (iexpc(ij)-1)*4 + 1
          p(1) = dipx(ij)
          p(2) = dipy(ij)
          p(3) = dipz(ij)
          p(4) = s(ij)
c
          do 80, igr = 1, ngran
c     3-----
c
c      -----  note the minus sign, due to the fact that this potential
c             is always coupled back to the electrons
c
            vext(ij,igr) = -ddot(4,p,1,zfa(jf,igr),1)
c           vext(ij,igr) = -adotb(p,zfa(jf,igr),4)
c           if (abs(vext(ij,igr)) .lt. threshs) vext(ij,igr) = zero
c     3-----
   80     continue
c   2-----
   90   continue
c 1-----
  100 continue
c
      return
      end
      subroutine extdcal(npol3,ndimb,nexp,num,fact,vr,vext,rx,ry,rz,s,
     1                   iexpc,dipind)
c------
c      this routine calculates the expanded reaction potential due to
c      the induced dipoles and charges at the expansion centra
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension vr(npol3+ndimb,nexp*4+2)
      dimension vext(num*(num+1)/2),s(num*(num+1)/2),
     1          rx(num*(num+1)/2),ry(num*(num+1)/2),rz(num*(num+1)/2)
      dimension iexpc(num*(num+1)/2)
      dimension dipind(npol3+ndimb)
c
c----- local arrays
c
      dimension f(4)
c
      data zero, one, pt5 /0.0d00, 1.0d00, 0.5d00/
c
c-----  effective loop over charge distributions
c
      ij = 0
      do 100, i = 1, num
        do 200, j = 1, i
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
          f(1) = rx(ij)
          f(2) = ry(ij)
          f(3) = rz(ij)
          f(4) = s(ij)
c
c    -----  execute expansion
c
          vextd = zero
          do 300, k = 1, 4
            do 400, l = 1, npol3+ndimb
              vextd = vextd + f(k)*vr(l,ip+k)*dipind(l)
  400       continue
  300     continue
c
          vext(ij) = vext(ij) - fact*vextd
c
  200   continue
  100 continue
      return
      end
      subroutine arfanal(ieps,unucrep,
     +  eoneel,ekin,enua,etwoel,
     +  uqm,uscf,snucnuc,selel,snua,stwoel,
     +  smolnuc,smolel,snucmol,selmol,smolmol,
     +  suqm,upolqm,uneqnuc,uneqel,uneqqm,
     +  ustanuc,ustael,ustaqm,uclase,uclasd,uclasr,uclas,
     +  suclas,upolcl,uneqcl,ustacl,extnuc,extel,sextnuc,
     +  sextel,sextmol,selext,snucext,smolext,
     +  stotnuc,stotel,stotmol,stotext,stabtot,
     +  uelst,suint,uint,udisp,rdisp,
     +  repmod,upoleq,ucstst,ucstpl,uneq,
     +  usta,upolneq,ustaneq,uens,uclasg,uclaseg,
     +  uclasdg,uclasrg,suclasg,upolclg,uneqclg,ustaclg,extnucg,
     +  extelg,sxtnucg,sextelg,sxtmolg,
     +  selextg,snucxtg,smolxtg,stotxtg,
     +  uelstg,suintg,uintg,repmodg,xscm)
c------
c      performs the analysis of the reaction field options
c      the relay-equations are solved for the types of source field
c      required. the external source field may or may not be added
c      to the molecular source field
c      the induced dipoles/surface charges may also be stored for
c      subsequent use in a non-equilibrium rf calculation (isolsav=1)
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ijpair)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/mollab)
c
INCLUDE(comdrf/ene)
INCLUDE(comdrf/opt)
INCLUDE(comdrf/scfopt)
INCLUDE(comdrf/scfpar)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/enrgci)
c
INCLUDE(comdrf/elpinf)
c
INCLUDE(../m4/common/drfopt)
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfzfa)
INCLUDE(comdrf/grinf)
INCLUDE(comdrf/inddaf)
INCLUDE(comdrf/neqtex)
INCLUDE(comdrf/neqpar)
INCLUDE(comdrf/ciopt)
c
INCLUDE(comdrf/drfbem)
c
INCLUDE(comdrf/runpar)
c
      dimension xscm(*)
c
c-----  local variables
c
cmw the ngran in the local dims are all converted to mxgran
      dimension clas((mxgran*(mxgran+1))*5+1)
c
      dimension suclasg(mxgran,ngran)
      dimension upolclg(ngran), uneqclg(ngran), ustaclg(ngran),
     2          extnucg(ngran),
     3          extelg(ngran), sxtnucg(ngran), sextelg(ngran),
     4          sxtmolg(ngran), selextg(ngran), snucxtg(ngran),
     5          smolxtg(ngran), stotxtg(ngran), uelstg(ngran),
     6          suintg(ngran), uintg(ngran), repmodg(ngran)
c
      dimension uclasg(ngran*(ngran+1)/2), uclaseg(ngran*(ngran+1)/2),
     1 uclasdg(ngran*(ngran+1)/2),uclasrg(ngran*(ngran+1)/2)
c
      dimension ixind(mxgran+2)
c
      character*80 blank
      logical uhf,rohf,rgvb,rogvb
      character*8 scfwfn
c
      external qext
c
c-----  data statements
c
      data zero, pt5, one /0.0d00, 0.5d00, 1.0d00/
      data blank /' '/
      data scfwfn /'scf     '/
c
c-----  begin
c
c
c-----  set da-record numbers for reading rf contribution and matrices
c
      ilcl = 21
c
      if (ieps .eq. 0) then
        ilwt = 26
        ilvr = 28
        illur = 43
        ilindx = 44
      else
        ilwt = 30
        ilvr = 32
        illur = 47
        ilindx = 48
      endif
c
      if (mcupdt) then
        ilwt = ilwt + 1
        ilvr = ilvr + 1
        ilcl = 91
      endif
c
c-----  calculate mulliken charges and dipoles and dipole preserving
c       charges if either is required for the analysis
c
      if (((ifldin .le. 2) .or. (ifldout .le. 2))
     1    .and. (ieps .eq. 0)) then
c 1-----
        if (nodpe .eq. 1) then
c   2-----
          ilmc = 121
          ilmd = 122
          ildc = 124
          call dppop(xscm)
        else
          ilmc = 126
          ilmd = 127
          ildc = 129
          call dppope(xscm)
c   2-----
        endif
c 1-----
      endif
c
c-----  set memory partitioning
c
c     call setscm(i10)
c     call cmem(loadcm)
c     loc10 = loccm(xscm(i10))
c ----  allocate all memory
      lwor = igmem_max_memory()
      i10 = igmem_alloc(lwor)
      call arfmem(i10,ixwt,ixvr,ns,itsolv,nwtr,nwtc,ixind,ndim,neqdim,
     1     ifldin,ngran,ixhlp,neqrf,isolsav,ineqex,imomsav,ixhlp2,
     2     neq2eps,ixsf,ixr,idisadd,ixq,ixdip,nexp,ixi,ixd,ixol,nchd,
     3     ixrx,ixry,ixrz,ixone,ixoneh,ixie,last)
      need = last - i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for arfanal'
	call caserr('Insufficient memory for arfanal')
      endif
c     call setc(need)
c
c-----  if the iterative (direct) method for solving the
c       relay equations is not used, read the necessary
c       matrices
c
      if (((field(5:) .ne. ' ') .or. (iextdip .ne. 0))
     1    .and. (itsolv .eq. 0)) then
c 1-----
c
c  -----  wt-matrix
c
        call daread(idafdrf,iodadrf,xscm(ixwt),nwtr*nwtc,ilwt)
        ixzf = ixwt + ndim*nexp4
c
        if (ifldin .le. 2) then
c   2-----
c    -----  lu-decomposed relay and indx
c
          call daread(idafdrf,iodadrf,xscm(ixr),ndim*ndim,illur)
          call daread(idafdrf,iodadrf,xscm(ixi),ndim,ilindx)
c   2-----
        endif
c
        call clear(xscm(ixsf),ndim)
c 1-----
      endif
c
      if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0)) then
c 1-----
c  -----  start procedure for solving response equations
c
        if (ifldin .le. 2) then
c   2-----
c    -----  get distributed charges, located at the expansion centra
c
          if (ifldin .eq. 1) then
c     3-----
c      -----  dipole preserving charges
c
            call daread(idafdrf,iodadrf,xscm(ixq),nexp,ildc)
c     3-----
          else
c     3-----
c      -----  mulliken charges
c
            call daread(idafdrf,iodadrf,xscm(ixq),nexp,ilmc)
c
c      -----  get mulliken dipoles and calculate their field
c             if iterative method is not used
c
            call daread(idafdrf,iodadrf,xscm(ixdip),3*nexp,ilmd)
c
            if (itsolv .eq. 0)
     1      call dipsf(npol3,ndimb,nexp,
     2               xscm(ixwt),xscm(ixdip),xscm(ixsf))
c     3-----
          endif
c
c    -----  calculate field at polarizabilities and
c           boundary elements due to charges (added to that
c           of the dipoles if present)
c           if iterative method is not used
c
          if (itsolv .eq. 0)
     1    call qsf(npol3,ndimb,nexp,xscm(ixwt),xscm(ixq),xscm(ixsf))
c
c    -----  calculate the induced moments
c
          if (idisadd .eq. 0) then
c     3-----
            if (itsolv .eq. 0) then
c       4-----
c        -----  through inversion of relay matrix
c
c        -----  by the molecular field
c
              call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                    ndim,xscm(ixind(1)))
c
c        -----  by the field of the external charges
c
              do 1100, igr = 1, ngran
                call luelmf(xscm(ixr),xscm(ixzf+(igr-1)*ndim),
     1                      xscm(ixi),ndim,ndim,xscm(ixind(igr+1)))
 1100         continue
c       4-----
            else
c       4-----
c        -----  through iterative procedure
c
c        -----  by the molecular field
c
              call solqd(ndim,nexp,xscm(ixq),xscm(ixdip),
     1                   xscm(ixind(1)),xscm(ixhlp))
c
c        -----  by the field of the external charges
c
              do 1200, igr = 1, ngran
                call solext(xscm(ixind(igr+1)),xscm(ixhlp))
 1200         continue
c       4-----
            endif
c     3-----
          else
c     3-----
c      -----  add the field of the external charges and
c             calculate the induced moments
c
            do 1300, i = 1, ndim
              do 1250, igr = 1, ngran
                xscm(ixsf+i-1) = xscm(ixsf+i-1) +
     1                           xscm(ixzf+(igr-1)*ndim+i-1)
 1250         continue
 1300       continue
c
            if (itsolv .eq. 0) then
c       4-----
              call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                    ndim,xscm(ixind(1)))
            else
              call solqde(xscm(ixq),xscm(ixdip),xscm(ixind(1)))
c       4-----
            endif
c     3-----
          endif
c   2-----
        else
c   2-----
c    -----  calculate molecular source field from electron density
c
c    -----  read density
c
          call densrd(xscm(ixd),xscm(ixd+nchd),xscm(ixd+2*nchd))
c
c    -----  read overlap and first moment integrals
c
          call daread(idafh,ioda,xscm(ixol),nchd,12)
          call daread(idafh,ioda,xscm(ixrx),nchd,53)
          call daread(idafh,ioda,xscm(ixry),nchd,54)
          call daread(idafh,ioda,xscm(ixrz),nchd,55)
c
c    -----  read assignment vector
c
          call daread(idafdrf,iodadrf,xscm(ixie),nchd,2)
c
c    -----  sum electronic source field
c
          if (itsolv .eq. 0) then
c     3-----
            call elsf(npol3,ndimb,nexp,num,xscm(ixwt),xscm(ixd),
     1                xscm(ixol),xscm(ixrx),xscm(ixry),xscm(ixrz),
     2                xscm(ixsf),xscm(ixie))
c
c      -----  add nuclear source field if not treated separate
c             from electrons
c
            ixzfn = ixwt + ndim*(nexp4+ngran)
            if (ifldin .eq. 3) then
c       4-----
              do 1400, i = 1, ndim
                xscm(ixsf+i-1) = xscm(ixsf+i-1) + xscm(ixzfn+i-1)
 1400         continue
c       4-----
            endif
c
c      -----  calculate induced moments
c
c      -----  read -relay- and -indx-, overwriting density etc.
c
            call daread(idafdrf,iodadrf,xscm(ixr),ndim*ndim,illur)
            call daread(idafdrf,iodadrf,xscm(ixi),ndim,ilindx)
c
            if (ifldin .eq. 3) then
c       4-----
              if (idisadd .eq. 0) then
c         5-----
c          -----  molecular field only
c
                call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                      ndim,xscm(ixind(1)))
c
c          -----  external source field
c
                do 1500, igr =1 , ngran
                  call luelmf(xscm(ixr),xscm(ixzf+(igr-1)*ndim),
     1            xscm(ixi),ndim,ndim,xscm(ixind(igr+1)))
 1500           continue
c         5-----
              else
c         5-----
                do 1600, i = 1, ndim
                  do 1550, igr = 1, ngran
                    xscm(ixsf+i-1) = xscm(ixsf+i-1) +
     1                               xscm(ixzf+(igr-1)*ndim+i-1)
 1550             continue
 1600           continue
c
c          -----  molecular + classical (external) field
c
                call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                      ndim,xscm(ixind(1)))
c         5-----
              endif
c       4-----
            else
c       4-----
c        -----  electronic induced moments
c
              call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                    ndim,xscm(ixind(1)))
c
c        -----  nuclear induced moments
c
              call luelmf(xscm(ixr),xscm(ixzfn),xscm(ixi),ndim,
     1                    ndim,xscm(ixind(2)))
c
c        -----  externally induced moments
c
              do 1700, igr = 1, ngran
                call luelmf(xscm(ixr),xscm(ixzf+(igr-1)*ndim),
     1          xscm(ixi),ndim,ndim,xscm(ixind(igr+2)))
 1700         continue
c       4-----
            endif
c     3-----
          else
c     3-----
c      -----  iterative method
c
            call solexp
c     3-----
          endif
c   2-----
        endif
c
c  -----  the moments induced by the source fields
c         are in xscm(ixind(k),k=1,...)
c
c  -----  save induced moments for use in non-equilibrium rf calculation
c         if required
c
c       if ((isolsav .eq. 1)
c    2       .or. (ineqex .eq. 1)) then
c   2-----
c         call neqsav(idrfout,ndim,ifldin,ngran,idisadd,isolsav,ieps,
c    2                imomsav,ineqex,maxneq,nwtr,nwtc,ilvr,nchd,nexp,
c    3                neqgrp,xscm(ixhlp),xscm(ixind(1)),xscm(ixvr),
c    4                xscm(ixol),xscm(ixrx),xscm(ixry),xscm(ixrz),
c    5                xscm(ixie),xscm(ixone),xscm(ixone),xscm(ixoneh))
c   2-----
c       endif
c 1-----
      endif
c
c  -----  the moments induced by the source fields
c         are in xscm(ixind(k),k=1,...)
c
c  -----  save induced moments for use in non-equilibrium rf calculation
c         if required
c
        if ((isolsav .eq. 1)
     2       .or. (ineqex .eq. 1)) then
c   2-----
          call neqsav(idrfout,ndim,ifldin,ngran,idisadd,isolsav,ieps,
     2                imomsav,ineqex,maxneq,nwtr,nwtc,ilvr,nchd,nexp,
     3                neqgrp,xscm(ixhlp),xscm(ixind(1)),xscm(ixvr),
     4                xscm(ixol),xscm(ixrx),xscm(ixry),xscm(ixrz),
     5                xscm(ixie),xscm(ixone),xscm(ixone),xscm(ixoneh))
c   2-----
        endif
c
c-----  now calculate reaction field at the expansion centra
c
c-----  reset addresses from ixq onward
c       sf, q, dip, relay/d will be overwritten
c
      ixq = ixsf
c
      if (ifldout .le. 2) then
c 1-----
        ixdip = ixq + nexp
        if (ifldout .eq. 2) then
c   2-----
          last = ixdip + 3*nexp
        else
          last = ixdip + 1
c   2-----
        endif
c 1-----
      else
c 1-----
        ixd = ixq
        ixol = ixd + nchd
        ixone = ixol
        ixrx = ixol + nchd
        ixry = ixrx + nchd
        ixrz = ixry + nchd
        ixie = ixrz + nchd
        last = ixie + nchd
c 1-----
      endif
c
      if ((ifldin .ge. 3) .and. (idisadd .eq. 0)
     1   .and. (irevdis .eq. 1)) then
c 1-----
        if (ifldout .le. 2) then
c   2-----
          ixd = last
          ixaa = ixd + nchd
          ixexp = ixaa + nchd
          ixie = ixexp + 3*nexp
          last = ixie + nchd
c   2-----
        else
c   2-----
          ixaa = ixol
          ixexp = last
          last = ixexp + 3*nexp
c   2-----
        endif
c 1-----
      endif
c
      need = last -i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for arfanal'
	call caserr('Insufficient memory for arfanal')
      endif
c     call setc(need)
c
c-----  read non-equilibrium induced moments if required
c
      if ((neqrf .eq. 1) .and. (momsav .eq. 1)) then
c 1-----
c  -----  read appropriate moments and equilibrium polarisation energy
c
        call daread(idafind,iodaind,xscm(ixhlp),ndim,ineqix(nneq))
c
        if (neq2eps .eq. 1)
     1  call daread(idafino,iodaino,xscm(ixhlp2),ndim,ineqix(nneq))
c 1-----
      endif
c
c-----  initialize some contributions
c
      if (wfntyp .eq. scfwfn) then
        uscf = etot
      else
        uscf = eci(istci)
      endif
c
      if (ieps .eq. 0) uqm = uscf
c
      call initz(smolmol,selmol,snucmol,uclas,uclase,uclasd,uclasr,
     2     suclas,ngran,uintg(1),uelstg(1),sxtmolg(1),suintg(1),
     3     upolclg(1),uneqclg(1),ustaclg(1),suclasg(1,1),
     3     uint,uelst,extel,sextmol,
     4     smolext,suint,udisp,stotmol,stotext,sextel,selext,sextnuc,
     5     snucext,upoleq,uneqqm,ustaqm,uneqcl,ustacl)
c
      if (iextdip .eq. 1) then
        fact = one
      else
        fact = pt5
      endif
c
      if (ifldout .le. 2) then
c 1-----
c  -----  the reaction field is coupled back to distributed
c         charges (and dipoles).
c
c  -----  get distributed charges
c
        if (ifldout .eq. 1) then
c   2-----
c    -----  dipole preserving charges
c
          call daread(idafdrf,iodadrf,xscm(ixq),nexp,ildc)
c   2-----
        else
c   2-----
c    -----  mulliken charges
c
          call daread(idafdrf,iodadrf,xscm(ixq),nexp,ilmc)
c
c    -----  get mulliken dipoles and calculate their
c           contribution to the interaction with the rf
c           and the external charges
c
          call daread(idafdrf,iodadrf,xscm(ixdip),3*nexp,ilmd)
c
          if (ieps .eq. 0) then
c     3-----
            do 3200, igr = 1, ngran
c       4-----
              call dipext(igr,nexp,xscm(ixdip),extdip)
              uelstg(igr) = uelstg(igr) + extdip
              uelst = uelst + extdip
c       4-----
 3200       continue
c
            if (neqsta .eq. 1) then
c       4-----
c        -----  interaction with static external potential from previous
c               calculation
c
              call daread(idafsta,iodasta,xscm(ixvr),neqdim,
     2        ineqix(nneq))
              call dipneq(nexp,xscm(ixvr),xscm(ixdip),edipst)
              ustaqm = ustaqm + edipst
c       4-----
            endif
c     3-----
          endif
c
          if (field(5:) .ne. ' ') then
c     3-----
c      -----  molecular source field
c
            if (ifldin .le. 3) then
c       4-----
              if (idisadd .eq. 0) then
c         5-----
c          -----  separate external source field
c
                call diprf(npol3,ndimb,nexp,
     1               xscm(ixvr),xscm(ixdip),xscm(ixind(1)),smoldip)
                smolmol = smolmol + smoldip
c
                do 3300, igr = 1, ngran
c           6-----
                  call diprf(npol3,ndimb,nexp,
     1            xscm(ixvr),xscm(ixdip),xscm(ixind(igr+1)),sextdip)
c
c            -----  note: factor 0.5 because of double counting
c
c                  sxtmolg(igr) = sxtmolg(igr) + pt5*sextdip
c                  sextmol = sextmol + pt5*sextdip
                  sxtmolg(igr) = sxtmolg(igr) + sextdip
                  sextmol = sextmol + sextdip
c           6-----
 3300           continue
c         5-----
              else
c         5-----
                call diprf(npol3,ndimb,nexp,
     1               xscm(ixvr),xscm(ixdip),xscm(ixind(1)),stotdip)
                stotmol = stotmol + stotdip
c         5-----
              endif
c       4-----
            else
c       4-----
c        -----  separate contributions from electronic,
c               nuclear and external source fields
c
              call diprf(npol3,ndimb,nexp,
     1               xscm(ixvr),xscm(ixdip),xscm(ixind(1)),seldip)
              selmol = selmol + seldip
c
              call diprf(npol3,ndimb,nexp,
     1               xscm(ixvr),xscm(ixdip),xscm(ixind(2)),snucdip)
              snucmol = snucmol + snucdip
c
              do 3400, igr = 1, ngran
c         5-----
                call diprf(npol3,ndimb,nexp,
     1            xscm(ixvr),xscm(ixdip),xscm(ixind(igr+2)),sextdip)
c
c          -----  note: factor 0.5 because of double counting
c
c                sxtmolg(igr) = sxtmolg(igr) + pt5*sextdip
c                sextmol = sextmol + pt5*sextdip
                sxtmolg(igr) = sxtmolg(igr) + sextdip
                sextmol = sextmol + sextdip
c         5-----
 3400         continue
c       4-----
            endif
c     3-----
          endif
c
c    -----  interaction with non-equilibrium reaction field
c
          if ((neqrf .eq. 1) .and. (ieps .eq. 0)) then
c     3-----
            call neqdip(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           neqdim,neq2eps,itsolv,field,ilvr,
     3           xscm(ixvr),xscm(ixvr),xscm(ixdip),xscm(ixhlp),
     4           xscm(ixhlp2),uneqqm,uneqqmo)
c     3-----
          endif
c   2-----
        endif
c
c  -----  calculate interaction of the external charges and
c         rf with the distributed charges
c
        if (ieps .eq. 0) then
c   2-----
          do 3500, igr = 1, ngran
c     3-----
            call qext(igr,nexp,xscm(ixq),extq)
            uelstg(igr) = uelstg(igr) + extq
            uelst = uelst + extq
c     3-----
 3500     continue
c
          if (neqsta .eq. 1) then
c     3-----
c      -----  interaction with static external potential from previous
c             calculation
c
            call daread(idafsta,iodasta,xscm(ixvr),neqdim,ineqix(nneq))
            call qneq(nexp,xscm(ixvr),xscm(ixdip),eqst)
            ustaqm = ustaqm + eqst
c     3-----
          endif
c   2-----
        endif
c
        if (field(5:) .ne. ' ') then
c   2-----
c    -----  read -vr- matrix
c
          call daread(idafdrf,iodadrf,xscm(ixvr),nwtc*nwtr,ilvr)
          if (ifldin .le. 3) then
c     3-----
            if (idisadd .eq. 0) then
c       4-----
              call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(1)),eqrf)
              smolmol = smolmol + eqrf
c
              do 3600, igr = 1, ngran
c         5-----
                call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(igr+1)),eqrf)
c
c          -----  note: factor 0.5 because of double counting
c
c                sxtmolg(igr) = sxtmolg(igr) + pt5*eqrf
c                sextmol = sextmol + pt5*eqrf
                sxtmolg(igr) = sxtmolg(igr) + eqrf
                sextmol = sextmol + eqrf
c
                if (iextdip .ne. 1) then
                  suintg(igr) = suintg(igr) + sextmol
                  suint = suint + sextmol
                else
                  uelstg(igr) = uelstg(igr) + sextmol
                  uelst = uelst + sextmol
                endif
c         5-----
 3600         continue
c       4-----
            else
c       4-----
              call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(1)),eqrf)
              stotmol = stotmol + eqrf
c       4-----
            endif
c     3-----
          else
c     3-----
            call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(1)),eqrf)
            selmol = selmol + eqrf
c
            call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(2)),eqrf)
            snucmol = snucmol + eqrf
c
            smolmol = selmol + snucmol
c
            do 3700, igr = 1, ngran
c       4-----
              call qrf(npol3,ndimb,nexp,xscm(ixvr),xscm(ixq),
     1           xscm(ixind(igr+2)),eqrf)
c
c        -----  note: factor 0.5 because of double counting
c
c              sxtmolg(igr) = sxtmolg(igr) + pt5*eqrf
c              sextmol = sextmol + pt5*eqrf
              sxtmolg(igr) = sxtmolg(igr) + eqrf
              sextmol = sextmol + eqrf
c
              if (iextdip .ne. 1) then
                suintg(igr) = suintg(igr) + sextmol
                suint = suint + sextmol
              else
                uelstg(igr) = uelstg(igr) + sextmol
                uelst = uelst + sextmol
              endif
c       4-----
 3700       continue
c     3-----
          endif
c   2-----
        endif
c
c  -----  interaction with non-equilibrium reaction field
c
        if ((neqrf .eq. 1) .and. (ieps .eq. 0)) then
c   2-----
          call neqq(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           neqdim,neq2eps,itsolv,field,ilvr,
     3           xscm(ixvr),xscm(ixvr),xscm(ixq),xscm(ixhlp),
     4           xscm(ixhlp2),eqrf,eqrfo)
c
          uneqqm = uneqqm + eqrf
c   2-----
        endif
c 1-----
      else
c 1-----
c  -----  calculate interaction energy from electron density
c
c  -----  read density
c
        call densrd(xscm(ixd),xscm(ixd+nchd),xscm(ixd+2*nchd))
c
        if (ieps .eq. 0) then
c   2-----
c    -----  calculate one electron energy from vacuum hamiltonian
c
          call daread(idafdrf,iodadrf,xscm(ixone),nchd,11)
          eoneel = tracep(xscm(ixd),xscm(ixone),num)
c
c    -----  calculate kinetic energy
c
          call daread(idafh,ioda,xscm(ixone),nchd,41)
          ekin = tracep(xscm(ixd),xscm(ixone),num)
c
          enua = eoneel - ekin
          etwoel = ehf - eoneel
c         call daread(idafh,ioda,xscm(ixone),nchd,103)
c         etwoel = pt5*tracep(xscm(ixd),xscm(ixone),num)
c         uqm = unucrep + eoneel + etwoel
c
c    -----  interaction with external charges
c
          call daread(idafdrf,iodadrf,xscm(ixone),ngran*nchd,12)
c
          do 3800, igr = 1, ngran
            extelg(igr) = tracep(xscm(ixd),
     1                    xscm(ixone+(igr-1)*nchd),num)
            extel = extel + extelg(igr)
            uelstg(igr) = extnucg(igr) + extelg(igr)
 3800     continue
c
          uelst = extnuc + extel
c
          if (neqsta .eq. 1) then
c     3-----
c      -----  interaction with static external potential from previous
c             calculation
c
            call daread(idafdrf,iodadrf,xscm(ixone),nchd,181)
            ustael = tracep(xscm(ixd),xscm(ixone),num)
            ustaqm = ustael + ustanuc
c     3-----
          endif
c
          if ((neqrf .eq. 1) .and. (field(:4) .eq. 'scf')) then
c     3-----
c      -----  interaction with static external potential from previous
c             calculation
c
            call daread(idafdrf,iodadrf,xscm(ixone),nchd,14)
            urfel = tracep(xscm(ixd),xscm(ixone),num)
c     3-----
          endif
c   2-----
        endif
c
c  -----  read overlap and first moment integrals
c
        call daread(idafh,ioda,xscm(ixol),nchd,12)
        call daread(idafh,ioda,xscm(ixrx),nchd,53)
        call daread(idafh,ioda,xscm(ixry),nchd,54)
        call daread(idafh,ioda,xscm(ixrz),nchd,55)
c
c  -----  read assignment vector
c
        call daread(idafdrf,iodadrf,xscm(ixie),nchd,2)
c
        ixzfn = ixvr + ndim*(nexp4+ngran)
c
        if (field(5:) .ne. ' ') then
c   2-----
c    -----  read -vr- matrix
c
          call daread(idafdrf,iodadrf,xscm(ixvr),nwtc*nwtr,ilvr)
c
          if (ifldin .le. 3) then
c     3-----
            if (idisadd .eq. 0) then
c       4-----
              if (itsolv .eq. 0) then
c         5-----
c          -----  electronic contributions
c
                call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1             xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2             xscm(ixie),xscm(ixind(1)),smolel)
c
                do 4100, igr = 1, ngran
c           6-----
                  call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1             xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2             xscm(ixie),xscm(ixind(igr+1)),sxtel)
c
c            -----  note: factor 0.5 because of double counting
c
c                  sextelg(igr) = pt5*sxtel
c                  sextel = sextel+ pt5*sxtel
                  sextelg(igr) = sxtel
                  sextel = sextel+ sxtel
c
c          -----  nuclear contributions
c
c          -----  note: factor 0.5 because of double counting
c
             sxtnuc = 
     +       ddot(ndim,xscm(ixzfn),1,xscm(ixind(igr+1)),1)
c    +       adotb(xscm(ixzfn),xscm(ixind(igr+1)),ndim)
c            sxtnucg(igr) = pt5*sxtnuc
c            sextnuc = sextnuc + pt5*sxtnuc
             sxtnucg(igr) = sxtnuc
             sextnuc = sextnuc + sxtnuc
             sxtmolg(igr) = sextelg(igr) + sxtnucg(igr)
c           6-----
 4100           continue
c
                sextmol = sextel + sextnuc
                smolnuc = 
     +          ddot(ndim,xscm(ixzfn),1,xscm(ixind(1)),1)
c    +          adotb(xscm(ixzfn),xscm(ixind(1)),ndim)
c         5-----
              else
c         5-----
                call itrf
c         5-----
              endif
c
              smolmol = smolnuc + smolel
c
              if (iextdip .ne. 1) then
                do 4200, igr = 1, ngran
                  suintg(igr) = suintg(igr)
     1                        + sextelg(igr) + sxtnucg(igr)
 4200           continue
                suint = suint + sextel + sextnuc
              endif
c       4-----
            else
c       4-----
              if (itsolv .eq. 0) then
c         5-----
                call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1             xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2             xscm(ixie),xscm(ixind(1)),stotel)
                stotnuc = ddot(ndim,xscm(ixzfn),1,xscm(ixind(1)),1)
c               stotnuc = adotb(xscm(ixzfn),xscm(ixind(1)),ndim)
                stotmol = stotel + stotnuc
c         5-----
              endif
c       4-----
            endif
c     3-----
          else
c     3-----
            snua = zero
            call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1           xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2           xscm(ixie),xscm(ixind(1)),e2pol)
            stwoel = e2pol
            call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1           xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2           xscm(ixie),xscm(ixind(2)),snucel)
            snua = snua + snucel
c
            do 5100, igr = 1, ngran
c       4-----
              call elrf(npol3,ndimb,nexp,num,xscm(ixvr),xscm(ixd),
     1           xscm(ixrx),xscm(ixry),xscm(ixrz),xscm(ixol),
     2           xscm(ixie),xscm(ixind(igr+2)),sxtel)
c
c        -----  note: factor 0.5 because of double counting
c
c              sextelg(igr) = pt5*sxtel
c              sextel = sextel + pt5*sxtel
              sextelg(igr) = sxtel
              sextel = sextel + sxtel
c
c       -----  note: factor 0.5 because of double counting
c
              sxtnuc = ddot(ndim,xscm(ixzfn),1,xscm(ixind(igr+2)),1)
c             sxtnuc = adotb(xscm(ixzfn),xscm(ixind(igr+2)),ndim)
c              sxtnucg(igr) = pt5*sxtnuc
c              sextnuc = sextnuc + pt5*sxtnuc
              sxtnucg(igr) = sxtnuc
              sextnuc = sextnuc + sxtnuc
c
              sxtmolg(igr) = sxtel + sxtnuc
c       4-----
 5100       continue
c
c      -----  collect contributions
c
            selnuc = ddot(ndim,xscm(ixzfn),1,xscm(ixind(1)),1)
c           selnuc = adotb(xscm(ixzfn),xscm(ixind(1)),ndim)
            snua = snua + selnuc
            snucnuc = ddot(ndim,xscm(ixzfn),1,xscm(ixind(2)),1)
c           snucnuc = adotb(xscm(ixzfn),xscm(ixind(2)),ndim)
c
            sextmol = sextel + sextnuc
            snucmol = snucnuc + snucel
            selmol = stwoel + selnuc
            smolmol = stwoel + snua + snucnuc
c
            if (iextdip .ne. 1) then
              do 5200, igr = 1, ngran
                suintg(igr) = suintg(igr) + sextelg(igr) + sxtnucg(igr)
 5200         continue
              suint = suint + sextel + sextnuc
            endif
c     3-----
          endif
c   2-----
        endif
c
c    -----  interaction with non-equilibrium rf
c
        if ((neqrf .eq. 1) .and. (ieps .eq. 0)) then
c   2-----
          call neqel(momsav,mcupdt,nwtc,nwtr,npol3,ndimb,nexp,nneq,
     2           ngran,neqdim,neq2eps,itsolv,field,ilvr,
     3           xscm(ixvr),xscm(ixvr),xscm(ixd),xscm(ixol),xscm(ixrx),
     4           xscm(ixry),xscm(ixrz),xscm(ixie),xscm(ixhlp),
     5           xscm(ixhlp2),uneqel,uneqnuc)
c
          uneqqm = uneqqm + uneqel + uneqnuc
c
          if ((field(5:) .ne. ' ') .and. (itsolv .ne. 0))
     1        call daread(idafdrf,iodadrf,xscm(ixvr),nwtr*nwtc,ilvr)
c   2-----
        endif
c 1-----
      endif
c
      uclase = zero
      uclasd = zero
      uclasr = zero
c
      if (nodiscr .eq. 0) then
c 1-----
c  -----  contributions at external charges
c
        ngrpair = ngran*(ngran+1)/2
        call daread(idafdrf,iodadrf,clas,(ngrpair*10+1),ilcl)
c
        uclas = clas(1)
        call clear(uclasg,ngrpair)
c
        do 6100, igr = 1, ngran
c   2-----
          do 6050, jgr = 1, igr
c     3-----
            indxcl = ia(igr) + jgr
            if (iclinte .eq. 1) then
              uclaseg(indxcl) = clas(1+indxcl)
              uclasg(indxcl) = uclasg(indxcl) + uclaseg(indxcl)
              uclase = uclase + uclaseg(indxcl)
            endif
c
            if (iclintd .eq. 1) then
              uclasdg(indxcl) = clas(1+ngrpair+indxcl)
              uclasg(indxcl) = uclasg(indxcl) + uclasdg(indxcl)
              uclasd = uclasd + uclasdg(indxcl)
            endif
c
            if (iclintr .eq. 1) then
              uclasrg(indxcl) = clas(1+2*ngrpair+indxcl)
              uclasg(indxcl) = uclasg(indxcl) + uclasrg(indxcl)
              uclasr = uclasr + uclasrg(indxcl)
            endif
c     3-----
 6050     continue
c   2-----
 6100   continue
c
        if (neqsta .eq. 1) then
c   2-----
c    -----  interaction with external classical field from
c           a previous calculation
c
          call neqcl(1,momsav,ngran,ndim,nwtc,nwtr,nexp,nneq,
     2                 neqdim,mcupdt,neq2eps,xscm(ixvr),xscm(ixvr),
     3                 xscm(ixhlp),xscm(ixhlp2),
     4                 ustaclg,ustacl)
c   2-----
        endif
c
        ixzf = ixvr + ndim*(4*nexp)
c
        if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0)) then
c   2-----
c    -----  note: the factors 0.5 are introduced because of
c                 the form of the energy-expression, which
c                 double counts the screening of the interactions
c                 between charges from different sets
c                 (nuclei, electrons, molecule, external)
c
          if (ifldin .le. 3) then
c     3-----
            if (idisadd .eq. 0) then
c       4-----
              if (field(5:) .ne. ' ') then
c         5-----
                do 6200, igr = 1, ngran
c           6-----
                smolxt = 
     +          ddot(ndim,xscm(ixzf+(igr-1)*ndim),1,xscm(ixind(1)),1)
c    +          adotb(xscm(ixzf+(igr-1)*ndim),xscm(ixind(1)),ndim)
c               smolxtg(igr) = pt5*smolxt
c               smolext = smolext + pt5*smolxt
c               suint = suint + pt5*smolxt
                smolxtg(igr) = smolxt
                smolext = smolext + smolxt
                suint = suint + smolxt
c           6-----
 6200           continue
c         5-----
              endif
c
              do 6300, igr = 1, ngran
c         5-----
c          -----  loop over inducing classical groups
c
                do 6250, jgr = 1, ngran
c           6-----
c            -----  loop over receiving classical groups
c
                itmp = ixzf+(jgr-1)*ndim
                sucls = 
     +          ddot(ndim,xscm(itmp),1,xscm(ixind(igr+1)),1)
c    +          adotb(xscm(ixzf+(jgr-1)*ndim),xscm(ixind(igr+1)),ndim)
c
c            -----  factor 0.5 for double counting (cf sextel etc)
c                   for non-diagonal elements only
c
c            -----  for the calculation of the classical polarisation
c                   energy, count diagonal elements only
c
                  if (igr .eq. jgr) then
                    suclasg(igr,jgr) = sucls
                    upolclg(igr) = -pt5*sucls
                  else
c                    suclasg(igr,jgr) = pt5*sucls
                    suclasg(igr,jgr) = sucls
                  endif
c
c            -----  collect classical stabilisation enrgies; in this way
c                   the polarisation energy cost is defined by the total
c                   stabilisation, in stead of the sum of the polarisati
c                   energies of the constituent classical analysis group
c
                  suclas = suclas + sucls
c           6-----
 6250           continue
c         5-----
 6300         continue
c       4-----
            else
c       4-----
              do 6400, igr = 1, ngran
                itmp = ixzf+(igr-1)*ndim
                stotxt = 
     +         ddot(ndim,xscm(itmp),1,xscm(ixind(1)),1)
c    +         adotb(xscm(ixzf+(igr-1)*ndim),xscm(ixind(1)),ndim)
                stotxtg(igr) = stotxt
                stotext = stotext + stotxt
 6400         continue
c       4-----
            endif
c     3-----
          else
c     3-----
            if (field(5:) .ne. ' ') then
c       4-----
              do 6500, igr = 1, ngran
c         5-----
                itmp = ixzf+(igr-1)*ndim
                selxt = 
     +          ddot(ndim,xscm(itmp),1,xscm(ixind(1)),1)
c    +          adotb(xscm(ixzf+(igr-1)*ndim),xscm(ixind(1)),ndim)
                snucxt = 
     +          ddot(ndim,xscm(itmp),1,xscm(ixind(2)),1)
c    +          adotb(xscm(ixzf+(igr-1)*ndim),xscm(ixind(2)),ndim)
c
                selextg(igr) = selxt
                snucxtg(igr) = snucxt
                smolxtg(igr) = selxt + snucxt
                suintg(igr) = suintg(igr) + selxt + snucxt
c
                selext = selext + selxt
                snucext = snucext + snucxt
                smolext = smolext + selxt + snucxt
c         5-----
 6500         continue
c
              suint = suint + selext + snucext
c       4-----
            endif
c
            do 6600, igr = 1, ngran
c       4-----
              do 6700, jgr = 1, ngran
c         5-----
                itmp = ixzf+(jgr-1)*ndim
                sucls = 
     +          ddot(ndim,xscm(itmp),1,xscm(ixind(igr+2)),1)
c    +          adotb(xscm(ixzf+(jgr-1)*ndim),xscm(ixind(igr+2)),ndim)
c
c          -----  correct diagonal elements and gather polarisation ener
c
                if (igr .eq. jgr) then
                  suclasg(igr,jgr) = sucls
                  upolclg(igr) = -pt5*sucls
                else
c                  suclasg(igr,jgr) = pt5*sucls
                  suclasg(igr,jgr) = sucls
                endif
c
c          -----  collect classical stabilisation energies; in this way,
c                 the polarisation energy cost is defined by the total
c                 stabilisation, in stead of the sum of the polarisation
c                 energies of the constituent classical analysis groups
c
                suclas = suclas + sucls
c         5-----
 6700         continue
c       4-----
 6600       continue
c     3-----
          endif
c
          if ((ifldin .ge. 3) .and. (idisadd .eq. 0) .and.
     1        (irevdis .eq. 1) .and. (ieps .eq. 0)) then
c     3-----
c      -----  calculate reverse dispersion
c
            call calrdis(xscm(ixexp),xscm(ixie),xscm(ixd),
     1           xscm(ixaa),rdisp)
            uint = uint + rdisp
c     3-----
          endif
c   2-----
        endif
c
        if ((neqrf .eq. 1) .and. (ieps .eq. 0)) then
c   2-----
c    -----  interaction of external charges with non-equilibrium
c           reaction field
c
          call neqcl(0,momsav,ngran,ndim,nwtc,nwtr,nexp,nneq,neqdim,
     2                 mcupdt,neq2eps,xscm(ixvr),xscm(ixvr),
     3                 xscm(ixhlp),xscm(ixhlp2),uneqclg,uneqcl)
        endif
c 1-----
      endif
c
c-----  dispersion
c
      if ((field(5:) .ne. ' ') .and. (ifldin .ge. 3) .and.
     1    (gamdrf .ne. zero)) then
c 1-----
        if (((itwoeps .eq. 1) .and. (ieps .eq. 1)) .or.
     1       (itwoeps .eq. 0)) then
c   2-----
          uhf = scftyp .eq. 'uhf'
cxxx      rohf = (scftyp.eq. 'rhf') .and. (nopset.gt.0)
          rohf = (scftyp.eq. 'rhf') .and. (nseto.gt.0)
          rgvb = scftyp .eq. 'gvb'
          rogvb = rohf .or. rgvb
c
          ixda = i10
          ixdb = ixda + nchd
          ixv = ixdb + nchd
c
          if (rogvb) then
c     3-----
            ixomg = ixv + num*num
          else
            ixomg = ixv + 1
c     3-----
          endif
c
          ixs = ixomg + nwtc*nwtc
          ixdx = ixs + nchd
          ixdy = ixdx + nchd
          ixdz = ixdy + nchd
          ixrxx = ixdz + nchd
          ixryy = ixrxx + nchd
          ixrzz = ixryy + nchd
          ixrxy = ixrzz + nchd
          ixrxz = ixrxy + nchd
          ixryz = ixrxz + nchd
          ixie = ixryz + nchd
          ixij = ixie + nchd
          last = ixij + nchd
c
          need = last -i10
      if (need .ge. lwor) then
	write(iwr,*) ' Not enough memory for arfanal'
	call caserr('Insufficient memory for arfanal')
      endif
c         call setc(need)
c
          call calself(ieps,xscm(ixda),xscm(ixdb),xscm(ixv),
     1       xscm(ixomg),xscm(ixs),
     2       xscm(ixdx),xscm(ixdy),xscm(ixdz),
     3       xscm(ixrxx),xscm(ixryy),xscm(ixrzz),xscm(ixrxy),
     4       xscm(ixrxz),xscm(ixryz),xscm(ixie),
     5       xscm(ixij),selel,scoul,sexch,
     6       uhf,rohf,rgvb,rogvb)
c
c    -----  factor 0.5 because of inherent polarisation energy cost
c
          udisp = pt5*(selel + sexch)
c
          stwoelo = scoul + sexch
          smolmo = snucnuc + selel + stwoelo + snua
c   2-----
        endif
c 1-----
      endif
c
      if (isolsav .eq. 1) then
c 1-----
c  -----  save static external potential for later use
c
        call nqsavs(nwtc,ngran,nzfa,nexp,nqdim,nqcls,neqgrp,
     2              uclaseg,uclasdg,uclasrg,extnucg,repmodg,
     3              xscm(i10))
c 1-----
      endif
c
c-----  total reaction field energy
c
      stabtot = smolmol + suclas + suint + stotmol + stotext
c
c-----  calculate equilibrium polarisation energy
c       for the calculation of the polarisation energy
c       source and reaction field are to be due to and coupled
c       back to equivalent entities (electrons, nuclei, external
c       charges, molecule, total system)
c
      if ((field(5:) .ne. ' ') .or. (iextdip .eq. 1)) then
c 1-----
c        if (idisadd .eq. 0) then
c   2-----
c    -----  the sources are the molecular (qm treated) system
c           and the external charges
c
c    -----  note: factor 0.5 because of equilibrium condition
c
          if (field(5:) .ne. ' ') then
            upolqm = -pt5*smolmol
          else
            upolqm = zero
          endif
          upolcl = -pt5*suclas
c          upoleq = upolqm + upolcl
c        else
c
c    -----  the source field is the total field of the molecular
c           system and external charges
c
c    -----  note: factor 0.5 because of equilibrium condition
c
          upoleq = -pt5*stabtot
c   2-----
c        endif
c 1-----
      endif
c
c-----  collect contributions
c
      if (ieps .eq. 0) then
c 1-----
        if ((field(5:) .eq. 'scf') .and. (iarfcal .eq. 0)) then
          if (idisadd .eq. 0) then
            uqm = uqm - scffact*smolmol
            uqm = uqm - fact*suint
          else
            uqm = uqm - scffact*stabtot
          endif
        else if (field(5:) .eq. 'scf') then
c         uqm = uqm - smolmol
          uqm = uqm - stabtot
        endif
c
        do 7200, igr = 1, ngran
          uintg(igr) = uintg(igr) + uelstg(igr) + repmodg(igr)
 7200   continue
c
        uint = uint + uelst + repmod
c
        usta = ustaqm + ustacl
c
        uneq = uneqqm + uneqcl
c
        if (field(:4) .eq. 'scf') then
c   2-----
          uqm = uqm - uint
          if ((neqsta .eq. 1) .and. (ieps .eq. 0)) uqm = uqm - ustaqm
          if ((neqrf .eq. 1) .and. (ieps .eq. 0)) uqm = uqm - uneqqm
c   2-----
        endif
c
c-----  configuration total energy
c
        uens = uqm + uclas + uint + stabtot
        if ((field(5:) .ne. ' ') .or. (iextdip .ne. 0))
     1     uens = uens + upoleq
c 1-----
      endif
c
      if ((field(5:) .eq. 'scf') .and. (gamdrf .ne. zero) 
     1  .and. (iarfcal .eq. 0)) then
        uqm = uqm - udisp
        if (scffact .eq. 1.0) then
          uqm = uqm + pt5*selel
          uens = uens + pt5*selel
        endif
      endif
c
c     if ((field(5:) .eq. 'pert') .and. (gamdrf .ne. zero))
      if ((gamdrf .ne. zero) .and.
     1 ((field(5:) .eq. 'pert') .or. (iarfcal .eq. 1)))
     1    uens = uens + udisp
c
      if ((neqsta .eq. 1) .and. (ieps .eq. 0))
     1   uens = uens + usta + ustaneq
c
      if ((neqrf .eq. 1) .and. (ieps .eq. 0))
     1   uens = uens + uneq + upolneq
c
      suqm = smolmol
c
      if (isolsav .eq. 1) then
c 1-----
        call nqcost(field,mxgran,ngran,ieps,idisadd,imomsav,itwoeps,
     2           eps1,eps2,kappa1,kappa2,neqgrp,
     3           uclaseg,uclasdg,uclasrg,suclasg,smolxtg,sxtmolg,
     4           snucnuc,stwoel,snucel,selnuc,
     5           upolclg,ucstst,ucstpl,stabtot)
c 1-----
      endif
c
c-----  reset core memory
c
c     call setc(loadcm)
      call gmem_free(i10)
c
      return
      end
      subroutine arfmem(i10,
     +     ixwt,ixvr,ns,itsolv,nwtr,nwtc,ixind,ndim,neqdim,
     1     ifldin,ngran,ixhlp,neqrf,isolsav,ineqex,imomsav,ixhlp2,
     2     neq2eps,ixsf,ixr,idisadd,ixq,ixdip,nexp,ixi,ixd,ixol,nchd,
     3     ixrx,ixry,ixrz,ixone,ixoneh,ixie,last)
c------
c      sets addresses in work array for reaction field analysis,
c      according to options given by
c      itsolv, ifldin, neqrf, isolsav, ineqex, imomsav, neq2eps and idis
c------
       dimension ixind(ngran+2)
c      include 'rfin/sizesrf'
c      dimension ixind(mxgran+2)
c      dimension ixind(3000)
c
c-----  begin
c
c     if (ngran+2.gt.3000) call caserr(
c    + 'drfsub.f: arfmem: ngran>3000 1 increase dimension of ixind >3000')
      ixwt = i10
      ixvr = i10
      ns = i10
c
c-----  solution via inversion of relay matrix or iterative procedure
c
      if (itsolv .eq. 0) then
c 1-----
        nwtvr = max(nwtr*nwtc,neqdim)
        ixind(1) = ixwt + nwtvr
        ixind(2) = ixind(1) + ndim
c
        if (ifldin .eq. 4) then
c   2-----
          do 100, igr = 1, ngran
            ixind(igr+2) = ixind(igr+1) + ndim
  100     continue
          ixhlp = ixind(ngran+2) + ndim
c
          if ((neqrf .eq. 1) .or. (isolsav .eq. 1)
     1                       .or. (ineqex .eq. 1)) then
c     3-----
            ixhlp2 = ixhlp + ndim
c
            if (neq2eps .eq. 1) then
c       4-----
              ixsf = ixhlp2 + ndim
            else
              ixsf = ixhlp2 + 1
c       4-----
            endif
c     3-----
          else
c     3-----
            ixsf = ixhlp + 1
c     3-----
          endif
c
          ixr = ixsf + ndim
          ns = ns + 2
c   2-----
        else
c   2-----
c    -----  field of external charges (discrete classical system)
c           added to source field or not?
c
          if (idisadd .eq. 0) then
c     3-----
            do 200, igr = 1, ngran
              ixind(igr+1) = ixind(igr) + ndim
  200       continue
            ixhlp = ixind(ngran+1) + ndim
c
            if ((neqrf .eq. 1) .or. (isolsav .eq. 1)
     1                         .or. (ineqex .eq. 1)) then
c       4-----
              ixhlp2 = ixhlp + ndim
c
              if (neq2eps .eq. 1) then
c         5-----
                ixsf = ixhlp2 + ndim
              else
                ixsf = ixhlp2 + 1
c         5-----
              endif
c       4-----
            else
c       4-----
              ixsf = ixhlp + 1
c       4-----
            endif
c
            ns = ns + 1
c     3-----
          else
c     3-----
            ixhlp = ixind(2) + 1
c
            if ((neqrf .eq. 1) .or. (isolsav .eq. 1)
     1                         .or. (ineqex .eq. 1)) then
c       4-----
              ixhlp2 = ixhlp + ndim
c
              if (neq2eps .eq. 1) then
c         5-----
                ixsf = ixhlp2 + ndim
              else
                ixsf = ixhlp2 + 1
c         5-----
              endif
c       4-----
            else
c       4-----
              ixsf = ixhlp + 1
c       4-----
            endif
c     3-----
          endif
c
          ixq = ixsf + ndim
c   2-----
        endif
c
        if (ifldin .le. 2) then
c   2-----
          ixdip = ixq + nexp
c
          if (ifldin .eq. 2) then
c     3-----
            ixr = ixdip + 3*nexp
          else
            ixr = ixdip + 1
c     3-----
          endif
c   2-----
        else if (ifldin .eq. 3) then
c   2-----
          ixr = ixsf + ndim
c   2-----
        endif
c
        ixi = ixr + ndim*ndim
c
c  -----  the density, overlap and first moment integrals
c         are not needed at the same time as the -relay-
c         and -indx- arrays, so they get the same address
c
        ixd = ixr
        lastd = ixd
c
        if ((ifldin .ge. 3) .or. (ineqex .eq. 1)
     2       .or. ((isolsav .eq. 1) .and. (imomsav .ne. 1))) then
c   2-----
          ixol = ixd + nchd
          ixrx = ixol + nchd
          ixry = ixrx + nchd
          ixrz = ixry + nchd
          ixone = ixrz + nchd
c
          if (ineqex .eq. 1) then
            ixoneh = ixone + nchd
          else
            ixoneh = ixone + 1
          endif
c
          if ((isolsav .eq. 1) .and. (imomsav .ne. 1)) then
            ixoneh = ixone + max(nchd,nwtc)
          endif
c
          ixie = ixoneh + nchd
          lastd = ixie + nchd
c   2-----
        endif
c
        last = max(lastd, ixi + ndim)
c 1-----
      else
c 1-----
c  -----  iterative method
c
        ixind(1) = 1
        ixind(2) = ixind(1) + ndim
c
        if (ifldin .eq. 4) then
c   2-----
          do 400, igr = 1, ngran
            ixind(igr+2) = ixind(igr+1) + ndim
  400     continue
          ixhlp = ixind(ngran+2) + ndim
c   2-----
        else
c   2-----
          if (idisadd .eq. 0) then
c     3-----
            do 500, igr = 1, ngran
              ixind(igr+1) = ixind(igr) + ndim
  500       continue
c     3-----
          endif
c
          ixq = ixind(2) + 1
c   2-----
        endif
c
        if (ifldin .le. 2) then
c   2-----
          ixdip = ixq + nexp
c
          if (ifldin .eq. 2) then
c     3-----
            ixhlp = ixdip + 3*nexp
          else
            ixhlp = ixdip + 1
c     3-----
          endif
c   2-----
        else
c   2-----
          ixdip = ixq + 1
          ixhlp = ixdip + 1
c   2-----
        endif
c
        ixd = ixhlp + ndim
c
        if ((ifldin .ge. 3) .or. (ineqex .eq. 1)) then
c   2-----
          ixol = ixd + nchd
          ixrx = ixol + nchd
          ixry = ixrx + nchd
          ixrz = ixry + nchd
          if (ineqex .eq. 1) then
            ixone = ixrz + nchd
            ixoneh = ixone + nchd
          else
            ixone = ixrz + 1
            ixoneh = ixone + 1
          endif
          ixie = ixoneh + nchd
          last = ixie + nchd
c   2-----
        else
c   2-----
          last = ixd
c   2-----
        endif
c 1-----
      endif
c
c     call setc(last)
c
      return
      end
      subroutine densrd(da,db,dx)
c------
c      reads electron density and sums it into da if necessary
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c------
      implicit REAL  (a-h,o-y),integer  (i-n)
      implicit character*8 (z)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension da(nchd), db(nchd), dx(nchd)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/opt)
INCLUDE(comdrf/scfopt)
INCLUDE(comdrf/scfpar)
      character *10 zbnam
      common/runlab/zcom(29),zanam(maxat),zbnam(maxorb),ztag(maxat),
     * zsymm,zgroup,zscftp,zruntp,zguess
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
c
c-----  local variables
c
      logical uhf, rohf, rgvb, rogvb
c
c-----  data statements
c
c-----  begin
c
c-----  determine the wavefunction type
c
      uhf=scftyp.eq.'uhf'
      rohf=scftyp.eq.'rhf'.and.nseto.gt.0
      rgvb=scftyp.eq.'gvb'
      rogvb=rohf.or.rgvb
c
c----- read density
c
      call daread(idafh,ioda,da,nchd,16)
c
      if (zscftp(1:3) .eq. 'uhf') then
c
c  -----  read beta density, and add this to alpha density
c         difference density in beta density
c
        call daread(idafh,ioda,db,nchd,20)
        do 10, i = 1, nchd
          dum=da(i)
          da(i) = dum + db(i)
          db(i) = dum - db(i)
   10   continue
      endif
c
      return
      end
      subroutine initz(smolmol,selmol,snucmol,uclas,uclase,uclasd,
     2  uclasr,suclas,ngran,uintg,uelstg,sxtmolg,suintg,
     3  upolclg,uneqclg,ustaclg,suclasg,uint,uelst,extel,sextmol,
     4  smolext,suint,udisp,stotmol,stotext,sextel,selext,sextnuc,
     5  snucext,upoleq,uneqqm,ustaqm,uneqcl,ustacl)
c------
c      initialises contributions
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension uintg(ngran),uelstg(ngran),sxtmolg(ngran),
     2     suintg(ngran),upolclg(ngran),uneqclg(ngran),ustaclg(ngran)
c
      dimension suclasg(ngran,ngran)
c
      data zero /0.0d00/
c
c-----  begin
c
      smolmol = zero
      selmol = zero
      snucmol = zero
c
      uclas = zero
      uclase = zero
      uclasd = zero
      uclasr = zero
      suclas = zero
c
      do 3100, igr = 1, ngran
c 1-----
        uintg(igr) = zero
        uelstg(igr) = zero
        sxtmolg(igr) = zero
        suintg(igr) = zero
        upolclg(igr) = zero
        uneqclg(igr) = zero
        ustaclg(igr) = zero
c
        do 3050, jgr = 1, ngran
c   2-----
          suclasg(igr,jgr) = zero
c   2-----
 3050   continue
c 1-----
 3100 continue
c
      uint = zero
      uelst = zero
      extel = zero
c
      sextmol = zero
      smolext = zero
      suint = zero
      udisp = zero
c
      stotmol = zero
      stotext = zero
c
      sextel = zero
      selext = zero
c
      sextnuc = zero
      snucext = zero
c
      upoleq = zero
c
      uneqqm = zero
      uneqcl = zero
c
      ustaqm = zero
      ustacl = zero
c
      return
      end
      subroutine dipext(igran,ndip,dip,extdip)
c------
c      calculates the interaction energy between the external charges
c      and the mulliken dipoles, reproducing part of the
c      molecular field
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension dip(3,ndip)
c
c-----  common blocks
c
INCLUDE(comdrf/drfzfa)
c
c-----  data statements
c
      data zero /0.0d00/
c
c-----  begin
c
c-----  initialize energy
c
      extdip = zero
c
c-----  loop over dipoles
c
      do 100, idip = 1, ndip
c
c  -----  calculate pointer in -zfa-
c         -zfa- contains the field and potential due to the
c         external charges at the expansion centra
c
        ip = (idip - 1)*4
c
        do 200, k = 1, 3
          extdip = extdip + dip(k,idip) * zfa(igran,ip + k)
  200   continue
  100 continue
      return
      end
      subroutine qext(igran,nq,q,extq)
c------
c      calculates the interaction energy between the external charges
c      and charges (dp, mulliken, nuclear) at the expansion centra
c      the potential of all external charges at the expansion centra
c      is calculated in drfzfa, and is stored in -zpa-
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension q(nq)
c
c-----  common blocks
c
INCLUDE(comdrf/drfzfa)
c
c-----  data statements
c
      data zero /0.0d00/
c
c-----  begin
c
c-----  initialize energy
c
      extq = zero
c
c-----  loop over charges  (expansion centra)
c
      do 100, iq = 1, nq
        extq = extq + q(iq) * zpa(igran,iq)
  100 continue
      return
      end
      subroutine calrdis(xexp,iexpc,d,aa,rdisp)
c------
c
c       calculates the "reverse" dispersion contribution to
c       the interaction energy between the qm system
c       and the external system. it is put into vdis2.
c
c       this contribution consists of the sum of contributions
c       from the so-called "special atoms", selected from the
c       external system, that are close to the qm system.
c
c       the "reverse" dispersion term is calculated as
c
c                  a(ext)**2/3 sum(ij) <a(ij)>
c     - 0.25  sum -----------------------------
c             ext      r(ext;ij)**6
c
c       in which a is the polarizability and (ij) denotes
c       charge distribution ij, and <> denotes an average
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension xexp(3,nexp)
      dimension iexpc(nchd)
      dimension aa(nchd), d(nchd)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/extinf)
c
INCLUDE(comdrf/rad)
c
c-----  local arrays
c
      dimension pq(3)
      data zero,pt5,one,pt375,two,three,four
     1 /0.0d00,0.5d00,1.0d00,0.375d00,2.0d00,3.0d00,4.0d00/
      data pt6 /.1666666666666666666667d00/
      data pt66/.6666666666666666666666d00/
c
c-----  begin
c
      call densrd(d(1),d(1+nchd),d(1+2*nchd))
c
c-----  get internal reduced average polarizability of
c       distribution ij in -aa-
c
      call daread(idafh,ioda,aa,nchd,213)
      if(idrfout.eq.3) call hatout(aa,num,num,3,'aa/d')
c
      call daread(idafdrf,iodadrf,xexp,3*nexp,1)
c     nchd2 = lenint(nchd)
c     call daread(idafdrf,iodadrf,iexpc,nchd2,2)
      call daread(idafdrf,iodadrf,iexpc,nchd,2)
c
      rdisp = zero
c
c-----  loop over the "special" external atoms
c       these have been selected in subr. drfzfa
c
      do 200, n = 1, nspec
c
c  -----  get pointer to "special" atom
c
        nspe = nspecl(n)
c
c  -----  calculate factor according to polarizability of
c         the "special" atom
c
        fluct = gamdrf*(polars(n)**pt66)
c
c  -----  loop over charge distributions
c
        ij = 0
        do 300, i = 1, num
          do 400, j = 1, i
            ij = ij + 1
            fact = two
            if (j .eq. i) fact = one
            ijexp = iexpc(ij)
            call distab(xexp(1,ijexp),xpts(1,nspe),pq,dist)
c
            factf = one
            v = one
c
c      -----  account for penetration effects (optional)
c
            if (modxza .ne. 0) then
              s = (abs(polars(n)*aa(ij)))**pt6
              if (ithole .eq. 1) then
                if (dist .le. afact*s) then
                  v = dist/(afact*s)
                  factf = (four*v**3 - three*v**4)**2
                endif
              else if (ithole .eq. 2) then
                au = afact*(dist/s)
                factf = (one - (pt5*au**2 + au + one)*exp(-au))**2
              endif
            endif
            dmind1 = one/dist
            dmin6 = factf*(dmind1**6)
            rdisp = rdisp - fluct*dmin6*aa(ij)*fact*d(ij)
  400     continue
  300   continue
  200 continue
c
      return
      end
      subroutine calself(ieps,da,db,v,omega,s,dipx,dipy,dipz,
     1                  rxx,ryy,rzz,rxy,rxz,ryz,iexpc,ijbit,
     2                  selel,scoul,sexch,uhf,rohf,rgvb,rogvb)
c------
c       this routine computes the interaction of a unit
c       charge distribution with the reaction field
c       induced by itself. the contributions are put into
c       selel, scoul and sexch
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gvalue)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension s(nchd), dipx(nchd), dipy(nchd), dipz(nchd)
      dimension rxx(nchd), ryy(nchd), rzz(nchd),
     1          rxy(nchd), rxz(nchd), ryz(nchd)
      dimension v(num,num)
      dimension da(nchd), db(nchd)
      dimension iexpc(nchd), ijbit(nchd)
      dimension omega(nwtc,nwtc)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/ihelp)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/opt)
INCLUDE(comdrf/shlint)
c
INCLUDE(../m4/common/drfopt)
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfind)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfint)
caleko
INCLUDE(comdrf/alpbet)
caleko
c
c-----  local arrays
c
      dimension cdum(3)
      dimension rr(4,4)
      dimension zfij(4)
c
      logical uhf,rohf,rgvb,rogvb
      logical core,open,pair
c
c-----  data statements
c
      data thresh2 /1.0d-06/
c
      data zero,one,two,pt5,pt25/0.0d00,1.0d00,2.0d00,0.5d00,0.25d00/
c
c-----  begin
c
      if (ieps .eq. 0) then
        ilomg = 50
      else
        ilomg = 51
      endif
c
      call densrd(da(1),da(1+nchd),da(1+2*nchd))
c
      core = ncorb.ne.0
      open = nopset.ne.0
      pair = npairs.ne.0
c
c-----  set up mapping of fock matrices
c
      if (rogvb) then
c 1-----
        none = 1
        if (.not. core) none = 0
        if (core) then
c   2-----
          do 10, i = 1, ncorb
            ihlp(i) = 1
   10     continue
c   2-----
        endif
        nop = 0
c
        if (open) then
c   2-----
          do 30, iset = 1, nopset
c     3-----
            iop = nopen(iset)
            do 20, i = 1, iop
              ihlp(ncorb+nop+i) = none + iset
   20       continue
            nop = nop + iop
c     3-----
   30     continue
c   2-----
        endif
c
        if (pair) then
c   2-----
          ngem = 2*npairs
          do 40, igem = 1, ngem
            ihlp(ncorb+nop+igem) = none + nopset + igem
   40     continue
c   2-----
        endif
c
        norb = ncorb + nop + 2*npairs
        nco1 = ncorb + 1
c
c  -----  get vectors and constuct density
c
        call daread(idafh,ioda,v,num*num,15)
        if (core) then
c   2-----
          do 70, i = 1, num
c     3-----
            do 60, j = 1, i
c       4-----
              dum = zero
              do 50, k = 1, ncorb
                dum = dum + v(i,k)*v(j,k)
   50         continue
              ij = ia(i) + j
              db(ij) = dum
c       4-----
   60       continue
c     3-----
   70     continue
c   2-----
        endif
c 1-----
      endif
c
c-----  calculate the second moment integrals for the expansion
c       of the field
c
      cdum(1) = gx
      cdum(2) = gy  
      cdum(3) = gz  
      gx = zero
      gy = zero
      gz = zero
      call qmints(rxx,ryy,rzz,rxy,rxz,ryz)
c     call qdpint(rxx,ryy,rzz,rxy,rxz,ryz,cdum,num)
      gx = cdum(1)
      gy = cdum(2)
      gz = cdum(3)
c
c-----  read necessary integrals from da10
c
      call daread(idafh,ioda,s   ,nchd,12)
      call daread(idafh,ioda,dipx,nchd,53)
      call daread(idafh,ioda,dipy,nchd,54)
      call daread(idafh,ioda,dipz,nchd,55)
c
c-----  read necessary data from da31
c
c     nchd2 = lenint(nchd)
c     call daread(idafdrf,iodadrf,iexpc,nchd2,2)
      call daread(idafdrf,iodadrf,iexpc,nchd,2)
      call daread(idafdrf,iodadrf,omega,nwtc*nwtc,ilomg)
c
c-----  debug and test-output
c
      if(idrfout .eq. 3) then
        call hatout(rxx,num,num,3,'rxx')
        call hatout(ryy,num,num,3,'ryy')
        call hatout(rzz,num,num,3,'rzz')
        call hatout(rxy,num,num,3,'rxy')
        call hatout(rxz,num,num,3,'rxz')
        call hatout(ryz,num,num,3,'ryz')
      endif
c
      selel = zero
      scoul = zero
      sexch = zero
c
c-----  set up loop over charge distributions
c
      do 100, i = 1, num
c 1-----
c  -----  second loop over charge distibutions (lower triangle)
c
        do 200, j = 1, i
c   2-----
          ij = ia(i) + j
c
          ijbit(ij) = 1
c
          factij = two
          if (j .eq. i) factij = one
c
c    -----  get address of charge distribution and assign
c           overlap and first moments
c
          ijexp = iexpc(ij)
          p(1) = dipx(ij)
          p(2) = dipy(ij)
          p(3) = dipz(ij)
          p(4) = s(ij)
c
          call drfoab(ijexp,ijexp,nwtc,omega)
          call matvec(omgab,p,zfij,4,.false.)
          try = abs(ddot(4,zfij,1,p,1))
c         try = abs(adotb(zfij,p,4))
c
c         if (abs(try) .lt. cutoff) ijbit(ij) = 0
c
          if (ijbit(ij) .eq. 0) goto 200
c
c    -----  order overlap, first and second moments in rr
c
          do 210, k = 1, 3
            pk = p(k)
            rr(4,k) = pk
            rr(k,4) = pk
  210     continue
          rr(4,4) = p(4)
          rr(1,1) = rxx(ij)
          rr(2,2) = ryy(ij)
          rr(3,3) = rzz(ij)
          rr(1,2) = rxy(ij)
          rr(2,1) = rr(1,2)
          rr(1,3) = rxz(ij)
          rr(3,1) = rr(1,3)
          rr(2,3) = ryz(ij)
          rr(3,2) = rr(2,3)
c
          self = ddot(16,rr,1,omgab,1)
c         self = adotb(rr,omgab,16)
          self = self*factij*da(ij)
          selel = selel + self
c
c    -----  loop over charge distributions
c
          do 1000, k = 1, i
c     3-----
            lmax = k
            if (i .eq. k) lmax = j
            do 1100, l = 1, lmax
c       4-----
              kl = ia(k) + l
              if (ijbit(kl) .eq. 0) goto 1100
c
c        -----  calculate addresses
c
              ik = ia(i) + k
              il = ia(i) + l
              if (j .lt. k) goto 1150
              jk = ia(j) + k
              jl = ia(j) + l
              goto 1170
 1150         jk = ia(k) + j
              if (j .lt. l) goto 1160
              jl = ia(j) + l
              goto 1170
 1160         jl = ia(l) + j
 1170         continue
c
              ii = i
              jj = j
              kk = k
              ll = l
c
c        -----  calculate density factors
c
              call drfdab(da,db,v,dcoul,dexch,norb,rohf,rgvb,uhf)
c
c
c             if ((abs(dcoul) .gt. thresh2) .or.
c    1            (abs(dexch) .gt. thresh2)) then
c         5-----
                factkl = two
                if (l .eq. k) factkl = one
                fact = factij*factkl
                if ((i .eq. k) .and. (j .eq. l)) fact = pt5*fact
                val = zero
                call drftwo(s,dipx,dipy,dipz,iexpc,omega)
c
c          -----  factor 0.25 for double counting twice
c
                val = val*fact*pt25
                scoul = scoul + dcoul*val
                sexch = sexch + dexch*val
c         5-----
c             endif
c       4-----
 1100       continue
c     3-----
 1000     continue
c   2-----
  200   continue
c 1-----
  100 continue
c
      selel = selel*gamdrf
c
      return
      end
      subroutine drfdab(da,db,v,coul,exch,norb,rohf,rgvb,uhf)
c
c     -----  forms density factors for drf 2-el. integrals
c
c            for closed-shell/uhf: da=total density, db=beta density
c            for rohf/gvb:         db=closed-shell dens., v=vectors
c            coul=coulomb-factor,  exch=exchange-factor
c
c      --------  p.th. van duijnen, groningen, march 1990
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arguments
c
      dimension da(nchd), db(nchd)
      dimension v(num,num)
c
INCLUDE(../m4/common/infoa)
INCLUDE(comdrf/ihelp)
INCLUDE(comdrf/ijpair)
c
INCLUDE(comdrf/opt)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfind)
c
caleko
INCLUDE(comdrf/alpbet)
caleko

      logical core, open, pair
      logical uhf, rohf, rgvb, rogvb
c
      data zero,pt5,pt125 /0.0d+00,0.5d+00,0.125d00/
      data pt25 /0.25d00/
      data one,two,four,eight /1.d00,2.0d00,4.0d00,8.0d00/
c
      core = ncorb.ne.0
      open = nopset.ne.0
      pair = npairs.ne.0
      coul = zero
      exch = zero
      rogvb = rohf.or.rgvb
      nco1 = ncorb + 1
c
      if (rohf) then
        if (core) then
          coul = coul + db(ij)*db(kl)*two
          exch = exch - pt5*(db(ik)*db(jl) + db(il)*db(jk))
          if (open) then
            do 300, io = nco1, norb
              coul = coul + db(ij)*v(kk,io)*v(ll,io)
     1                    + db(kl)*v(ii,io)*v(jj,io)
              exch = exch - pt25*(v(jj,io)*(db(ik)*v(ll,io)
     1                    + db(il)*v(kk,io)) + v(ii,io)*(db(jl)*v(kk,io)
     2                    + db(jk)*v(ll,io)))
  300       continue
          endif
        endif
        if (open) then
          do 400, io = nco1, norb
c           iof = ihlp(io) - 1
            iof = ihlp(io)
            do 450, jo = nco1, norb
c             jof = ihlp(jo) - 1
              jof = ihlp(jo)
              iof0 = max(iof,jof)
              jof0 = min(iof,jof)
              iojo = ia(iof0) + jof0
              coul = coul + pt5*v(ii,io)*v(jj,io)*v(kk,jo)*v(ll,jo)
              exch = exch + pt125*(betaij(iojo) - one)
     1                    *v(ii,io)*v(jj,jo)*
     2                    (v(kk,io)*v(ll,jo) + v(ll,io)*v(kk,jo))
  450       continue
  400     continue
        endif
      else if (rgvb) then
        if (core) then
          coul = coul + alphij(1)*db(ij)*db(kl)
          exch = exch + pt5*betaij(1)*(db(ik)*db(jl) + db(il)*db(jk))
          if (open .or. pair) then
            do 500, io = nco1, norb
              iojo = ia(ihlp(io)) + 1
              coul = coul + alphij(iojo)*
     1        (db(ij)*v(kk,io)*v(ll,io) + db(kl)*v(ii,io)*v(jj,io))
              exch = exch + pt5*betaij(iojo)*
     1               (v(jj,io)*(db(ik)*v(ll,io)+db(il)*v(kk,io)) +
     2                v(ii,io)*(db(jl)*v(kk,io)+db(jk)*v(ll,io)))
  500       continue
          endif
        endif
        if (open .or. pair) then
          do 600, io = nco1, norb
            do 650, jo = nco1, norb
              iof = ihlp(io)
              jof = ihlp(jo)
              iojo = ia(iof) + jof
              if (jof .gt. iof) iojo = ia(jof) + iof
              coul = coul + alphij(iojo)*
     1                      v(ii,io)*v(jj,io)*v(kk,jo)*v(ll,jo)
              exch = exch + pt5*betaij(iojo)*v(ii,io)*v(jj,jo)*
     1                      (v(kk,io)*v(ll,jo) + v(ll,io)*v(kk,jo))
  650       continue
  600     continue
        endif
      else
        coul = coul + da(kl)*da(ij)*four
        exch = exch - (da(il)*da(jk) + da(ik)*da(jl))
        if (uhf) then
          exch = exch - (db(ik)*db(jl) + db(il)*db(jk))
        endif
      endif
      if (rogvb) then
        coul = eight*coul
        exch = eight*exch
      endif
      exch = gamdrf*exch
c
      return
      end
      subroutine drfself(s,dipx,dipy,dipz,ijbit,
     1                  rxx,ryy,rzz,rxy,rxz,ryz,omega,iexpc,vself)
c------
c       this routine computes the interaction of a unit
c       charge distribution with the reaction field
c       induced by itself. the contributions are put into
c       -vself-, later to be contracted with the density
c       matrix to calculate the self-energy.
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension s(nchd), dipx(nchd), dipy(nchd), dipz(nchd)
      dimension rxx(nchd), ryy(nchd), rzz(nchd),
     1          rxy(nchd), rxz(nchd), ryz(nchd)
      dimension omega(nomga)
      dimension vself(nchd)
      dimension iexpc(nchd), ijbit(nchd)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfexp)
c
c-----  local arrays
c
      dimension rr(4,4)
c
c-----  data statements
c
      data zero,one,pt5/0.0d00,1.0d00,0.5d00/
c
c-----  begin
c
c-----  debug and test-output
c
      if(idrfout .eq. 3) then
        call hatout(rxx,num,num,3,'rxx')
        call hatout(ryy,num,num,3,'ryy')
        call hatout(rzz,num,num,3,'rzz')
        call hatout(rxy,num,num,3,'rxy')
        call hatout(rxz,num,num,3,'rxz')
        call hatout(ryz,num,num,3,'ryz')
      endif
c
c-----  set up loop over charge distributions
c
      ij = 0
      do 100, ii = 1, num
c 1-----
c  -----  second loop over charge distibutions (lower triangle)
c
        do 200, jj = 1, ii
c   2-----
          ij = ij + 1
c
c    -----  skip small distributions (as computed in subr. drfone)
c
          if (ijbit(ij) .ne. 1) goto 200
c
c    -----  get address of charge distribution and assign
c           overlap and first moments
c
          ijexp = iexpc(ij)
          p(1) = dipx(ij)
          p(2) = dipy(ij)
          p(3) = dipz(ij)
          p(4) = s(ij)
c
c    -----  get formal self-interaction of distribution
                                                       
          call drfoab(ijexp,ijexp,nwtc,omega)
c
c    -----  order overlap, first and second moments in rr
c
          do 10, k = 1, 3
            pk = p(k)
            rr(4,k) = pk
            rr(k,4) = pk
   10     continue
          rr(4,4) = p(4)
          rr(1,1) = rxx(ij)
          rr(2,2) = ryy(ij)
          rr(3,3) = rzz(ij)
          rr(1,2) = rxy(ij)
          rr(2,1) = rr(1,2)
          rr(1,3) = rxz(ij)
          rr(3,1) = rr(1,3)
          rr(2,3) = ryz(ij)
          rr(3,2) = rr(2,3)
c
c    -----  evaluate self-interaction
c
          self = ddot(16,rr,1,omgab,1)
c         self = adotb(rr,omgab,16)
c
c    -----  collect self-interactions of charge distributions
c    -----  note: factor 0.5 for self-consistency:
c           polarisation is always accounted for!!!
c
          vself(ij) = pt5*self*gamdrf
c   2-----
  200   continue
c 1-----
  100 continue
c
      if(idrfout.eq.3) call hatout(vself,num,num,3,'vself')
      return
      end
      subroutine drftwo(s,dx,dy,dz,iexpc,omega)
c
c      --------  p.th. van duijnen, ibm-kingston 1985 -----
c
c  * * *  routine computes two-electron contribution of reaction field
c
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      dimension s(nchd),dx(nchd),dy(nchd),dz(nchd)
      dimension iexpc(nchd)
      dimension omega(nwtc,nwtc)
c
INCLUDE(comdrf/iofil)
INCLUDE(comdrf/ijpair)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfexp)
INCLUDE(comdrf/drfint)
c
      dimension sump(4),sumq(4)
c
      data zero, pt5 /0.0d00,0.5d00/
c
      ij = ia(max(i,j)) + min(i,j)
      ijexp = iexpc(ij)
      p(1) = dx(ij)
      p(2) = dy(ij)
      p(3) = dz(ij)
      p(4) = s(ij)
c
      kl = ia(max(k,l)) + min(k,l)
      klexp = iexpc(kl)
      q(1) = dx(kl)
      q(2) = dy(kl)
      q(3) = dz(kl)
      q(4) = s(kl)
c
      call drfoab(ijexp,klexp,nwtc,omega)
      call matvec(omgab,q,sumq,4,.false.)
c
      call drfoab(klexp,ijexp,nwtc,omega)
      call matvec(omgab,p,sump,4,.false.)
c
c-----  note: + sign because source and recipient are electrons
c
      val = ddot(4,sump,1,q,1) + 
     +      ddot(4,sumq,1,p,1)
c     val = adotb(sump,q,4) +  adotb(sumq,p,4)
c
c-----  test output
c
      if (idrfout .ge. 4) write(iwr,*) i, j, k, l, val
c
      return
      end
      subroutine drfdisp2(xexp,iexpc,ijbit,vdis2,aa)
c------
c
c       calculates the "reverse" dispersion contribution to
c       the interaction energy between the qm system
c       and the external system. it is put into vdis2.
c
c       this contribution consists of the sum of contributions
c       from the so-called "special atoms", selected from the
c       external system, that are close to the qm system.
c
c       the "reverse" dispersion term is calculated as
c
c                  a(ext)**2/3 sum(ij) <a(ij)>
c     - 0.25  sum -----------------------------
c             ext      r(ext;ij)**6
c
c       in which a is the polarizability and (ij) denotes
c       charge distribution ij, and <> denotes an average
c
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
c-----  dummy arrays
c
      dimension xexp(3,nexp)
      dimension iexpc(nchd), ijbit(nchd)
      dimension vdis2(nchd), aa(nchd)
c
c-----  common blocks
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
c
INCLUDE(comdrf/rad)
c
c-----  local arrays
c
      dimension pq(3)
      data zero,pt5,one,pt375,three,four
     1 /0.0d00,0.5d00,1.0d00,0.375d00,3.0d00,4.0d00/
      data pt66/.6666666666666666666666d00/
      data pt6 /0.166666666666666666667d00/
c
c-----  begin
c
c-----  get internal reduced average polarizability of
c       distribution ij in -aa-
c
      call daread(idafh,ioda,aa,nchd,213)
      if(idrfout.eq.3) call hatout(aa,num,num,3,'aa/d')
c
c-----  loop over the "special" external atoms
c       these have been selected in subr. drfzfa
c
      do 100, n = 1, nspec
c
c  -----  get pointer to "special" atom
c
        nspe=nspecl(n)
c
c  -----  calculate factor according to polarizability of
c         the "special" atom
c
        fluct = gamdrf*(polars(n)**pt66)
c
c  -----  loop over charge distributions
c
        ij = 0
        do 200, i = 1, num
          do 300, j = 1, i
            ij = ij + 1
c
c    -----  skip numerically small charge distributions
c           (according to threshold set in subr. drfone)
c
            if (ijbit(ij) .eq. 0) goto 300
            ijexp = iexpc(ij)
            call distab(xexp(1,ijexp),xpts(1,nspe),pq,dist)
c
            factf = one
            v = one
c
c      -----  account for penetration effects (optional)
c
            if (modxza .ne. 0) then
              s = (abs(polars(n)*aa(ij)))**pt6
              if (ithole .eq. 1) then
                if (dist .le. afact*s) then
                  v = dist/(afact*s)
                  factf = (four*v**3 - three*v**4)**2
                endif
              else if (ithole .eq. 2) then
                au = afact*(dist/s)
                factf = (one - (pt5*au**2 + au + one)*exp(-au))**2
              endif
            endif
            dmind1 = one/dist
            dmin6 = factf*(dmind1**6)
            vdis2(ij) = vdis2(ij) - fluct*dmin6*aa(ij)
  300     continue
  200   continue
  100 continue
      if(idrfout.eq.3) call hatout(vdis2,num,num,3,'disp2')
      return
      end
      subroutine drfexc(gg,sso,gso,eso)
c-----
c       exact electrostatic field *** not implemented
c-----
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
      logical norm,double
      logical block,blocks,blocki
      character*8 block1, block2
      logical some,out
INCLUDE(comdrf/tim)
INCLUDE(comdrf/opt)
      common/iofile/ir,iw,ip,ipadiofil(20)
INCLUDE(comdrf/dafil)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/nshel)
INCLUDE(comdrf/bas)
INCLUDE(comdrf/stv)
INCLUDE(comdrf/ssgg)
INCLUDE(comdrf/ijpair)
INCLUDE(comdrf/sym)
INCLUDE(comdrf/rys)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/extinf)
c
      dimension gg(*),sso(*),gso(*),eso(*)
      dimension dij(225)
      dimension ijx(35),ijy(35),ijz(35)
      dimension xv(5,5,5),yv(5,5,5),zv(5,5,5)
      character *8 errmsg(3)
      data errmsg /'program ','stop in ','-drfexc-'/
      data block1 /'blocks  '/
      data block2 /'blocki  '/
      data maxrys /5/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data pi212  /1.1283791670955d+00/
      data sqrt3  /1.73205080756888d+00/
      data sqrt5  /2.23606797749979d+00/
      data sqrt7  /2.64575131106459d+00/
      data ijx    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      data ijy    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      data ijz    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
c
      tol=rln10*itol
c
      noprt = 0
c
      out =noprt.eq.0.and.nprint.eq.3
      some=noprt.eq.0.and.nprint.ne.-5
      norm=normf.ne.1.or.normp.ne.1
      blocks=blktyp.eq.block1.and.noprt.eq.0
      blocki=blktyp.eq.block2.and.noprt.eq.0
      block =blocks.or.blocki
c
c----- calculate -vext- matrix with exact attraction integrals-----
c
c     ----- ishell -----
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmax(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell -----
c
      do 8000 jj=1,ii
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmax(jj)
      locj=kloc(jj)-minj
c
      ljtmod=ljt+2
c
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt-2)/2+1
      if(nroots.gt.maxrys) then
         write(iw,9997) maxrys,lit,ljt,nroots
         call hnderr(3,errmsg)
      endif
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
       g(ij)=zero
  100 continue
c
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      arri=ai*rr
      axi=ai*xi
      ayi=ai*yi
      azi=ai*zi
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
cmw      cgi=cg(ig)
c
c     ----- j primitive -----
c
      jgmax=j2
      if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
      aj=ex(jg)
      aa=ai+aj
      aa1=one/aa
      dum=aj*arri*aa1
      if(dum.gt.tol) go to 6000
      fac= exp(-dum)
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
cmw      cgj=cg(jg)
      ax=(axi+aj*xj)*aa1
      ay=(ayi+aj*yj)*aa1
      az=(azi+aj*zj)*aa1
c
c     ----- density factor -----
c
      double=iandj.and.ig.ne.jg
      ij=0
      do 360 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi*fac
      go to 220
  120 dum1=cpi*fac
      go to 220
  130 dum1=cdi*fac
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi*fac
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
cmw  180 dum1=cgi*fac
  180 go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      jmax=maxj
      if(iandj) jmax=i
      do 360 j=minj,jmax
      go to (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),j
  230 dum2=dum1*csj
      if(.not.double) go to 350
      if(i.gt.1) go to 240
      dum2=dum2+dum2
      go to 350
  240 dum2=dum2+csi*cpj*fac
      go to 350
  250 dum2=dum1*cpj
      if(double) dum2=dum2+dum2
      go to 350
  260 dum2=dum1*cdj
      if(double) dum2=dum2+dum2
      go to 350
  270 if(norm) dum2=dum2*sqrt3
      go to 350
  280 dum2=dum1*cfj
      if(double) dum2=dum2+dum2
      go to 350
  290 if(norm) dum2=dum2*sqrt5
      go to 350
  300 if(norm) dum2=dum2*sqrt3
      go to 350
cmw  310 dum2=dum1*cgj
  310 if(double) dum2=dum2+dum2
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      ij=ij+1
  360 dij(ij)=dum2
c
c     ----- nuclear attraction -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ic=1,nxtpts
      znuc=-chrg(ic)
      if (znuc .eq. zero) goto 500
      cx=xpts(1,ic)
      cy=xpts(2,ic)
      cz=xpts(3,ic)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      tt=one/(aa+uu)
      t = sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i
      call sxyz
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint*ww
  410 continue
  420 continue
c
      ij=0
      do 450 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 440 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      dum=zero
      do 430 iroot=1,nroots
  430 dum=dum+xv(ix,jx,iroot)*yv(iy,jy,iroot)*zv(iz,jz,iroot)
      dum=dum*(aa1*pi212)
      ij=ij+1
      g(ij)=g(ij)+dum*dij(ij)
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
c     ----- set up overlap and h-core matrices -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=ia(loci+i)+(locj+j)
      gg(nn)= g(ij)
 7500 continue
c
c     ----- salc transformation -----
c
cnot      if(block) call stvsym(blocki,sso,gso,eso)
 8000 continue
 9000 continue
cnot      if(out) call stvprt(blocks,blocki,ss,gg,sso,gso,num)
c
      num2=(num*(num+1))/2
      if(some) write(iw,9998)
      ntim=0
      if(some) ntim=1
      call texit(ntim,ntim)
      return
 9998 format(' ...... end of exact external charge integrals ......')
 9997 format(' in -drfexc- the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i2,/,
     2       ' lit,ljt,nroots= ',3i3)
      end
      subroutine rfupd(totcvg)
c------
c      updates one-electron average rf contribution from
c      reduced one-density matrix and adds it to the
c      1-electron hamiltonian h0 in the ao basis
c
c      also checks convergence of the rf contribution
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
      logical totcvg
c
INCLUDE(comdrf/scm)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/dafil)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/ciopt)
c
INCLUDE(comdrf/runpar)
c
c-----  begin
c
c-----  calculate current 1-electron density matrix
c
cnot      call gugadm
c
c-----  count # of passes through update
c
      irfpas = irfpas + 1
c
      totcvg = .false.
c
      if (mcupdt) then
        ilwt = ilwt + 1
        ilvr = ilvr + 1
      endif
c
      ixd = 1
      ixd2 = ixd + nchd + 1
c
c-----  read density
c
      call densrd(xscm(ixd),xscm(ixd+nchd),xscm(ixd+2*nchd))
c
      call daread(idafdrf,iodadrf,xscm(ixd2),nchd,70)
c
cnot      diff = cvgden(xscm(ixd),xscm(ixd2),nchd)
c
c-----  write actual density on da31 for future reference
c
      call dawrit(idafdrf,iodadrf,xscm(ixd),nchd,70,navdrf)
c
c     if ((diff .lt. rfcvg) .and. (irfpas .gt. 1) .and. (extra)) then
      if (irfpas .gt. 1 .and. extra) then
c 1-----
        totcvg = .true.
c
      else
c 1-----
c       if ((diff .lt. rfcvg) .and. (irfpas .gt. 1)) extra = .true.
        if (irfpas .gt. 1) extra = .true.
c 1-----
      endif
c
c-----  set addresses
c
      ixwt = 1
      ixvr = 1
      ixzf = ixwt + ndim*nexp4
      ixzfn = ixwt + ndim*(nexp4+ngran)
c
      ixind = ixwt + nwtr*nwtc + 1
      ixsf = ixind + ndim + 1
      ixd = ixsf + ndim + 1
      ixol = ixd + nchd + 1
      ixrx = ixol + nchd + 1
      ixry = ixrx + nchd + 1
      ixrz = ixry + nchd + 1
      ixone = ixrz + nchd + 1
      ixr = ixd
      ixi = ixr + ndim*ndim + 1
      ixie = max((ixi + ndim + 1), (ixone + nchd + 1))
      last = ixie + ndim + 1
c
      call setc(last)
c
c-----  initialise source field
c
      call clear(xscm(ixsf),ndim)
c
c-----  read wt-matrix
c
      call daread(idafdrf,iodadrf,xscm(ixwt),nwtr*nwtc,ilwt)
c
c-----  read density
c
      call densrd(xscm(ixd),xscm(ixd+nchd),xscm(ixd+2*nchd))
c
c-----  read overlap and first moment integrals
c
      call daread(idafh,ioda,xscm(ixol),nchd,12)
      call daread(idafh,ioda,xscm(ixrx),nchd,53)
      call daread(idafh,ioda,xscm(ixry),nchd,54)
      call daread(idafh,ioda,xscm(ixrz),nchd,55)
c
c-----  read assignment vector
c
      call daread(idafdrf,iodadrf,xscm(ixie),nchd,2)
c
c-----  sum electronic source field
c
      call elsf(npol3,ndimb,nexp,num,xscm(ixwt),xscm(ixd),
     1                xscm(ixol),xscm(ixrx),xscm(ixry),xscm(ixrz),
     2                xscm(ixsf),xscm(ixie))
c
c-----  calculate induced moments
c
c-----  read -relay- and -indx-, overwriting density etc.
c
      call daread(idafdrf,iodadrf,xscm(ixr),ndim*ndim,illur)
      call daread(idafdrf,iodadrf,xscm(ixi),ndim,ilindx)
c
      call luelmf(xscm(ixr),xscm(ixsf),xscm(ixi),ndim,
     1                      ndim,xscm(ixind))
c
c-----  read vr-matrix
c
      call daread(idafdrf,iodadrf,xscm(ixvr),nwtr*nwtc,ilvr)
c
c-----  read density
c
      call densrd(xscm(ixd),xscm(ixd+nchd),xscm(ixd+2*nchd))
c
c-----  read overlap and first moment integrals
c
      call daread(idafh,ioda,xscm(ixol),nchd,12)
      call daread(idafh,ioda,xscm(ixrx),nchd,53)
      call daread(idafh,ioda,xscm(ixry),nchd,54)
      call daread(idafh,ioda,xscm(ixrz),nchd,55)
c
c-----  read assignment vector
c
      call daread(idafdrf,iodadrf,xscm(ixie),nchd,2)
c
c-----  calculate reaction potential in ao basis
c
      call onerf(xscm(ixvr),xscm(ixol),
     1  xscm(ixrx),xscm(ixry),xscm(ixrz),
     2  xscm(ixone),xscm(ixie),xscm(ixind))
c
c-----  add updated reaction potential to one-electron
c       hamiltonian
c-----  read previous h0 (overwriting overlap matrix)
c
      call daread(idafdrf,iodadrf,xscm(ixol),nchd,19)
c
      do 200, i = 1, nchd
        xscm(ixol+i-1) = xscm(ixol+i-1) + scffact*xscm(ixone+i-1)
  200 continue
c
c-----  store updated one-electron hamiltonian
c
      call dawrit(idafh,ioda,xscm(ixol),nchd,11,nav)
c
      return
      end
      subroutine onerf(vr,s,dx,dy,dz,vrf,iexpc,dipind)
c------
c      calculates expanded potential due to non-equilibrium rf
c------
      implicit REAL  (a-h,o-z),integer  (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(comdrf/iofil)
INCLUDE(../m4/common/infoa)
c
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/rfene)
c
INCLUDE(comdrf/runpar)
c
      dimension vr(nwtr,nwtc)
      dimension s(nchd), dx(nchd), dy(nchd), dz(nchd)
      dimension iexpc(nchd)
      dimension vrf(nchd), dipind(ndim)
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
            do 400, l = 1, ndim
c       4-----
              vrfd = vrfd + f(k)*vr(l,ip+k)*dipind(l)
c       4-----
  400       continue
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
      subroutine mkdiagm(polar,polten,nel)
      implicit none
      REAL polar, polten, zero
      integer nel, i, j
      dimension polten(nel,nel)
      data zero /0.0d00/
      do i = 1, nel
        do j = 1, i
          if (i .eq. j) then
            polten(i,i) = polar
          else
            polten(i,j) = zero
            polten(j,i) = zero
          endif
        enddo
      enddo
      return
      end
      subroutine arfupd(q,enadd,rfield)
c-----
      implicit REAL (a-h, o-z)
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
c
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sector)
INCLUDE(../m4/common/dump3)
INCLUDE(../m4/common/runlab)
INCLUDE(comdrf/darw)
INCLUDE(comdrf/dafil)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfpar)
INCLUDE(comdrf/rfene)
INCLUDE(comdrf/runpar)
      character*8 zrhf, zcas
      character*8 rfield
      dimension zcas(3)
      dimension q(*)
c
      data zrhf /'rhf'/,zcas   /'casscf','mcscf','vb'/
c
      enadd = 0.d00
c
      idens = igmem_alloc(nx)
      ione = igmem_alloc(nx)
      iarf = igmem_alloc(nx)
      iol = igmem_alloc(nx)
      idx = igmem_alloc(nx)
      idy = igmem_alloc(nx)
      idz = igmem_alloc(nx)
      iwt = igmem_alloc(nwtr*nwtc)
      ir = igmem_alloc(ndim*ndim)
      ii = igmem_alloc(ndim)
      ifsf = igmem_alloc(ndim)
      isf = igmem_alloc(ndim)
      iind = igmem_alloc(ndim)
      iexp = igmem_alloc(nx)
c
c
      call daread(idafh,ioda,q(idens),nx,16)
      call dawrit(idafdrf,iodadrf,q(idens),nx,102,navdrf)
c
      if (rfield(5:8) .eq. 'scf ') then
      call daread(idafh,ioda,q(iol),nx,12)
      call daread(idafh,ioda,q(idx),nx,53)
      call daread(idafh,ioda,q(idy),nx,54)
      call daread(idafh,ioda,q(idz),nx,55)
c
      call daread(idafdrf,iodadrf,q(ifsf),ndim,101)
      call daread(idafdrf,iodadrf,q(ii),ndim,44)
      call daread(idafdrf,iodadrf,q(iexp),nx,2)
      call daread(idafdrf,iodadrf,q(ir),ndim*ndim,43)
c
c
c-----  calculate average reaction field
c       from one-electron density
c
      call arfone(
     1 q(iarf),q(idens),
     2 q(iol),q(idx),q(idy),q(idz),
     3 q(ir),q(ii),q(iwt),q(iexp),
     4 q(ifsf),q(iind),q(isf),enadd)
c
c-----  read one-electron hamiltonian
c
      endif
      call daread(idafh,ioda,q(ione),nx,11)
c
c-----  add 1-e ARF contribution and write
c
      if (rfield(5:8) .eq. 'scf ') 
     1 call vadd(q(ione),1,q(iarf),1,q(ione),1,nx)
c
c-----  add 1-e electrostatic contribution and write
c
      if(zscftp.eq.zcas(1).or.zscftp.eq.zcas(2) .or.
     1   zscftp.eq.zcas(3).and.
     2   rfield(1:4) .eq. 'scf') then
          i80  = igmem_alloc(nx*ngran)
          call daread(idafdrf,iodadrf,q(i80),nx*ngran,12)
           do 25, i = 1, ngran
             call addup(q(ione),q(i80+(i-1)*nx),q(ione),nx)
   25      continue
          call gmem_free(i80)
          if (rfield(:4) .eq. 'scf') then
            enadd = enadd - extnuc(iactst) 
            if (neqsta .eq. 1) enadd = enadd - ustanuc(iactst)
            if (neqrf .eq. 1) enadd = enadd - uneqnuc(iactst)
            enadd = enadd - repmod(iactst)
          endif
      endif
      call wrt3(q(ione),nx,ibl3f,numdu)
c
c-----  release memory
c
      call gmem_free(iexp)
      call gmem_free(iind)
      call gmem_free(isf)
      call gmem_free(ifsf)
      call gmem_free(ii)
      call gmem_free(ir)
      call gmem_free(iwt)
      call gmem_free(idz)
      call gmem_free(idy)
      call gmem_free(idx)
      call gmem_free(iol)
      call gmem_free(iarf)
      call gmem_free(ione)
      call gmem_free(idens)
c
      return
      end
      subroutine arfcvg(q,energy,enadd,oarfcvg,narfpts)
      implicit REAL (a-h, o-z)
      dimension q(*)
      logical oarfcvg
c
INCLUDE(../m4/common/sizes)
INCLUDE(comdrf/sizesrf)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/dump3)
INCLUDE(../m4/common/scfopt)
INCLUDE(comdrf/dafil)
INCLUDE(comdrf/darw)
INCLUDE(comdrf/drfdaf)
INCLUDE(comdrf/drfpar)
c
      if (conci_drf.le.0.0d0) conci_drf = acurcy
c
      oarfcvg = .false.
      energy = energy + enadd
c
      i10 = igmem_alloc(nx)
      i20 = igmem_alloc(nx)
      call daread(idafh,ioda,q(i10),nx,16)
c     call rdedx(q(i10),nx,ibl3pa,idaf)
      call daread(idafdrf,iodadrf,q(i20),nx,102)
      diff =  cvgden(q(i10),q(i20),nx)
c
      call gmem_free(i20)
      call gmem_free(i10)
      if (narfpts.lt.1) return
      if ((diff .le. conci_drf).or.
     1    (narfpts+1.gt.maxci_drf)) oarfcvg = .true.
c
      write(iwr,1) energy,diff
1     format(1x,78('='),/,' Average Reaction Field Energy',f20.10,
     1       /,' Convergence on 1-el density',d10.2,/,1x,78('='),/,1x)
      if (narfpts.eq.0) write(iwr,2) maxci_drf,conci_drf
2     format(1x,78('='),/,' maxcycci ',i3,' threshci ',f9.7,/,
     2       1x,78('='))
      if ((diff .lt. conci_drf) .and. oarfcvg)
     1    write(iwr,3)
3     format(1x,78('='),/,' Calculation converged on 1-e density',
     1       /,1x,78('='))
      if ((narfpts+1.gt.maxci_drf) .and. oarfcvg) write(iwr,4)
4     format(1x,78('='),/,' Maximum Average Reaction Field ',
     1      'cycles exceeded',/,1x,78('='))
c
      return
      end
