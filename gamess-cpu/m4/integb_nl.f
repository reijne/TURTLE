c ******************************************************
c ****** ROUTINES for non-LOCAL Pseudopotentials
c ******************************************************
**==xpsnlc.f
      subroutine xpsnlc(hh,f,v,h,s,ict,iso,l2,natom,nbasis,nshels,
     *                  oatomic)
      implicit real*8  (a-h,o-z)
c *** extended to handle g basis functions in the molecular
c *** basis, but only f-terms in the ECP
      logical iandj, oatomic
      character*8 blank,zlib,ztagp
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
      common/junkc/zlib(204),ztagp(10)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
      common/blkin/pxyz(4),ttt(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
     * ,nspacf,chglib(204),maxlib(204),ibllib(204)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/junk/space(450),
     +  xint,yint,zint,ta,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,
     +  apsl(20,4),cpsl(20,4),cpslo(400,4),
     +  libmax(4),libort(4),
     +  xin(225),yin(225),zin(225),dij(225),g(225),m1(225),
     +  ijx(225),ijy(225),ijz(225)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      dimension f(l2),ict(natom,*),s(225,*),iso(nshels,*)
      dimension v(nbasis,*)
      dimension hh(l2),h(l2)
c
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35)
      dimension mmin(5),mmax(5),maz(5)
c
      data m10/10/,m204/204/
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
      data sqrt3/1.73205080756888d0/,pi32/5.5683279968317d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data done/1.0d0/
      data blank/'        '/
      data mmin/1,2,5,11,21/
      data mmax/1,4,10,20,35/
      data maz/1,4, 9,17,29/
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      nav = lenwrd()
      itoli = 15
      tol = 2.30258d0*itoli
c retrieve ict from file
      call readi(ict,nw196(6)*nav,ibl196(6),idaf)
c
c retrieve header blocks from pseudopotential library
c
      mword = 1760 + 8/nav
      m408 = 204 + 408/nav
      call rdchr(zlib,m204,iblkpl,numlib)
      call reads(chglib,m408,numlib)
c
c     loop over atoms
c
      do 580 ic = 1 , natom
c
         if (zpseud(ic).ne.blank) then
c
c     skip if there is an equivalent centre with greater index
c
            do 50 it = 1 , nt
               icd = ict(ic,it)
               if (icd.gt.ic) go to 580
 50         continue
c
            call vclr(f,1,l2)
c
            numps = locatc(zlib,204,zpseud(ic))
            if (numps.eq.0) then
               call caserr2(
     +       'attempt to retrieve unknown pseudopotential from library')
            end if
 60         continue
c           chg = chglib(numps)
            npmax = maxlib(numps)
            call rdchr(ztagp,m10,ibllib(numps)+iblkpl,numlib)
            call reads(apsl,mword,numlib)
            if (npmax.le.0 .or. npmax.gt.4) then
               call caserr2(
     +          'invalid parameters from pseudopotential library')
            else
               do 460 np = 1 , npmax
                  itpmax = libmax(np)
                  itport = libort(np)
c                 ndim = itpmax*itport
                  lkt = np
                  mink = mmin(np)
                  maxk = mmax(np)
                  mazk = maz(np)
                  xi = c(1,ic)
                  yi = c(2,ic)
                  zi = c(3,ic)
                  kdif = mazk - mink + 1
                  kdifnp = kdif*itport*num
                  do 130 ig = 1 , itpmax
                     ax = apsl(ig,np)
                     ax = ax + ax
                     go to (70,80,90,100) , np
 70                  cnorm = pi32/(dsqrt(ax)*ax)
                     go to 110
 80                  cnorm = pi32*0.5d0/(dsqrt(ax)*ax*ax)
                     go to 110
 90                  cnorm = pi32*0.75d0/(dsqrt(ax)*ax**3)
                     go to 110
 100                 cnorm = pi32*1.875d0/(dsqrt(ax)*ax**4)
 110                 cnorm = 1.d0/dsqrt(cnorm)
                     do 120 iip = 1 , itport
                        ipp = (iip-1)*itpmax + ig
                        cpslo(ipp,np) = cpslo(ipp,np)*cnorm
 120                 continue
 130              continue
c
                  do 410 jj=1,nshell
                     if (oatomic.and.katom(jj).ne.ic) go to 410
                     j = katom(jj)
                     xj = c(1,j)
                     yj = c(2,j)
                     zj = c(3,j)
                     j1 = kstart(jj)
                     j2 = j1 + kng(jj) - 1
                     ljt = ktype(jj)
                     minj = kmin(jj)
                     maxj = kmax(jj)
                     locj = kloc(jj) - minj
                     jdif = maxj - minj + 1
                     rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                     ij = 0
                     do 160 i = mink , maxk
                        nnx = ix(i)
                        nny = iy(i)
                        nnz = iz(i)
                        do 150 j = minj , maxj
                           ij = ij + 1
                           ijx(ij) = nnx + jx(j)
                           ijy(ij) = nny + jy(j)
                           ijz(ij) = nnz + jz(j)
                           do 140 ig = 1 , itpmax
                              s(ij,ig) = 0.d0
 140                       continue
 150                    continue
 160                 continue
c
c
                     do 360 ig = 1 , itpmax
                        ai = apsl(ig,np)
                        arri = ai*rr
                        axi = ai*xi
                        ayi = ai*yi
                        azi = ai*zi
c
                        do 350 jg = j1 , j2
                           aj = ex(jg)
                           aa = ai + aj
                           dum = aj*arri/aa
                           if (dum.le.tol) then
                              fac = dexp(-dum)
                              ax = (axi+aj*xj)/aa
                              ay = (ayi+aj*yj)/aa
                              az = (azi+aj*zj)/aa
c
                              nn = 0
                              do 310 i = mink , maxk
                                 go to (170,170,170,170,170,170,170,180,
     +                                  210,210,170,170,170,190,210,210,
     +                                  210,210,210,200) , i
 170                             dum = fac
                                 go to 210
 180                             dum = fac*sqrt3
                                 go to 210
 190                             dum = fac*sqrt5
                                 go to 210
 200                             dum = fac*sqrt5*sqrt3
 210                             do 300 j = minj , maxj
                                 go to (220,230,290,290,
     +                           240,290,290,250,290,290,
     +                           260,290,290,270,290,290,290,290,
     +                           290,280,
     +                           282,290,290,284,290,290,290,290,
     +                           290,286,290,290,288,290,290), j
 220                                  dum1 = cs(jg)
                                      dum2 = dum*dum1
                                      go to 290
 230                                  dum1 = cp(jg)
                                      dum2 = dum*dum1
                                      go to 290
 240                                  dum1 = cd(jg)
                                      dum2 = dum*dum1
                                      go to 290
 250                                  dum2 = dum2*sqrt3
                                      go to 290
 260                                  dum1 = cf(jg)
                                      dum2 = dum*dum1
                                      go to 290
 270                                  dum2 = dum2*sqrt5
                                      go to 290
 280                                  dum2 = dum2*sqrt3
                                      go to 290
 282                                  dum1 = cf(jg)
                                      dum2 = dum*dum1
                                      go to 290
 284                                  dum2 = dum2*sqrt7
                                      go to 290
 286                                  dum2 = dum2*sqrt5/sqrt3
                                      go to 290
 288                                  dum2 = dum2*sqrt3
 290                                  nn = nn + 1
                                      dij(nn) = dum2
 300                             continue
 310                          continue
c
                              aa1 = done/aa
                              ta = dsqrt(aa1)
                              x0 = ax
                              y0 = ay
                              z0 = az
                              in = -5
                              do 330 i = 1 , lkt
                                 in = in + 5
                                 ni = i
                                 do 320 j = 1 , ljt
                                    jn = in + j
                                    nj = j
                                    call stvint
                                    xin(jn) = xint*ta
                                    yin(jn) = yint*ta
                                    zin(jn) = zint*ta
 320                             continue
 330                          continue
                              do 340 i = 1 , ij
                                nnx = ijx(i)
                                nny = ijy(i)
                                nnz = ijz(i)
                                s(i,ig) = s(i,ig) + dij(i)*xin(nnx)
     +                             *yin(nny)*zin(nnz)
 340                          continue
                           end if
 350                    continue
 360                 continue
                     if (lkt.ge.3) call comb(s,itpmax,lkt,jdif)
                     ij = 0
                     do 400 k = mink , mazk
                        do 390 j = minj , maxj
                           ij = ij + 1
                           lj = locj + j
                           do 380 io = 1 , itport
                              lk = (io-1)*kdif + k - mink + 1
                              v(lj,lk) = 0.d0
                              ipp = (io-1)*itpmax
                              do 370 ig = 1 , itpmax
                                 cx = cpslo(ipp+ig,np)
                                 v(lj,lk) = v(lj,lk) + s(ij,ig)*cx
 370                          continue
 380                       continue
 390                    continue
 400                 continue
c
 410              continue
c
                  iijj = 0
                  if (.not.oatomic) then
c
c                   loop over all basis functions
c
                    do 450 ii = 1 , num
                       do 440 jj = 1 , ii
                          iijj = iijj + 1
                          igs = -kdif - mink + 1
                          do 430 ij = 1 , itport
                             igs = igs + kdif
                             do 420 i = mink , mazk
                                i1 = igs + i
                                f(iijj) = f(iijj) + cpsl(ij,np)*v(ii,i1)
     +                                    *v(jj,i1)
 420                         continue
 430                      continue
 440                   continue
 450                continue
                  else
c
c                   one-centre terms only, skip all other terms
c
                    do 650 iii = 1 , nshell
                       ist = kloc(iii)
                       ied = ist+kmax(iii)-kmin(iii)
                       if (katom(iii).ne.ic) then 
                          iijj = ied*(ied+1)/2
                          go to 650
                       endif
                       do ii = ist,ied
                          do 640 jjj = 1 , iii
                             jst = kloc(jjj)
                             jed = jst+kmax(jjj)-kmin(jjj)
                             jed = min(jed,ii)
                             if (katom(jjj).ne.ic) then 
                                iijj = iijj + jed-jst+1
                                go to 640
                             endif
                             do jj = jst,jed
                                iijj = iijj + 1
                                igs = -kdif - mink + 1
                                do 630 ij = 1 , itport
                                   igs = igs + kdif
                                   do 620 i = mink , mazk
                                      i1 = igs + i
                                      f(iijj) = f(iijj) + cpsl(ij,np)
     +                                          *v(ii,i1)*v(jj,i1)
 620                               continue
 630                            continue
                             enddo
 640                      continue
                       enddo
 650                continue
                  endif
c
c
 460           continue
c
c     symmetrisation
c
c --- set up common /junk/ ready for subroutine rhr
               call ptrdtr
c
               do 570 it = 1 , nt
                  if (it.eq.1) then
                     icd = ic
                     call dcopy(l2,f,1,h,1)
                  else
                     icd = ict(ic,it)
                     if (icd.ge.ic) go to 570
                     itm = it - 1
                     do 470 it1 = 1 , itm
                        icd1 = ict(ic,it1)
                        if (icd1.eq.icd) go to 570
 470                 continue
                     ntr = invt(it)
c
                     do 490 ish = 1 , nshell
                        m1(ish) = iso(ish,it)
 490                 continue
c
                     do 550 ii = 1 , nshell
                        if (oatomic.and.katom(ii).ne.ic) go to 550
                        lit = ktype(ii)
                        mini = kmin(ii)
                        maxi = kmax(ii)
                        loci = kloc(ii) - mini
                        id = m1(ii)
                        locid = kloc(id) - mini
                        do 540 jj = 1 , ii
                           if (oatomic.and.katom(jj).ne.ic) go to 540
                           iandj = ii.eq.jj
                           ljt = ktype(jj)
                           minj = kmin(jj)
                           maxj = kmax(jj)
                           locj = kloc(jj) - minj
                           jd = m1(jj)
                           locjd = kloc(jd) - minj
                           maxt = maxj
                           do 510 i = mini , maxi
                              lci = loci + i
                              llc = iky(lci)
                              if (iandj) maxt = i
                              do 500 j = minj , maxt
                                 lcj = locj + j
                                 ttt(i,j) = f(llc+lcj)
                                 if (iandj) ttt(j,i) = ttt(i,j)
 500                          continue
 510                       continue
c
                           call rhr
c
                           do 530 i = mini , maxi
                              lcid = locid + i
                              if (iandj) maxt = i
                              do 520 j = minj , maxt
                                 lcjd = locjd + j
                                 llc = ind(lcid,lcjd)
                                 h(llc) = ttt(i,j)
 520                          continue
 530                       continue
 540                    continue
 550                 continue
                  end if
                  if (nprint.eq.3) then
                     write (iwr,6010) icd
                     do 560 i = 1 , num
                        ji = iky(i) + 1
                        jf = ikyp(i)
                        write (iwr,6020) (h(j),j=ji,jf)
 560                 continue
                  end if
c
                  call vadd(h,1,hh,1,hh,1,l2)
 570           continue
            end if
         end if
 580  continue
      return
 6010 format (//,' pseudopotential integrals /atom',i3)
 6020 format (10f12.8)
c
      end
**==ptrdtr.f
      subroutine ptrdtr
c
c --- restores the arrays ptr and dtr from the dumpfile and
c --- puts them into common /bufb/ where they can be used from
c --- other routines  e.g. rhr
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      logical ciopt, ciforc, mp2, hfgr, bfgs, ump2, lmeth2
      logical ump3, rmp3, ordmo, mp2w, loptor, ladp, lcpf
      logical lopti, lmcscf, lforce, lci, lcart, lmcdat
      logical lfdtrn, unit7, lcontr, lvcd, lgten
      logical ldenom, ignore, ldens, lset, ladapt, lsym, latmol
      logical berny, llibry, limpt, fpres, oss, ldiag, lskip
      logical opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common /restrl/ ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,
     +rmp3,ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
c
      dimension ypname(3),ydname(6)
      data ypname /'x'  ,'y'  ,'z'  /
      data ydname /'xx' ,'yy' ,'zz' ,'xy' ,'xz' ,'yz' /
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
      if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
      if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
      if (nprint.eq.1) then
         write (iwr,6010)
         write (iwr,6030)
         do 40 it = 1 , nt
            write (iwr,6060)
            write (iwr,6070) it
            np = 3*(it-1)
            write (iwr,6040) (ypname(j),j=1,3)
            write (iwr,6080)
            do 20 i = 1 , 3
               write (iwr,6050) ypname(i) , (ptr(i,np+j),j=1,3)
 20         continue
            if (odbas) then
            write (iwr,6020)
            nd = 6*(it-1)
            write (iwr,6040) (ydname(j),j=1,6)
            write (iwr,6080)
            do 30 i = 1 , 6
               write (iwr,6050) ydname(i) , (dtr(i,nd+j),j=1,6)
 30         continue
            write (iwr,6020)
            endif
 40      continue
      end if
c
      return
 6010 format (//'  ptr,dtr,ftr and gtr restored from the dumpfile')
 6020 format (//)
 6030 format (/,' transformation of the basis functions',/)
 6040 format (8x,10(3x,a4,3x))
 6050 format (2x,a4,2x,10f10.6)
 6060 format ('0')
 6070 format (/,21x,'transformation number',i4,/)
 6080 format (/)
      end
**==symd1.f
      subroutine symd1(df,ptr,ic,ict,natom)
      implicit real*8  (a-h,o-z)
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/junk/de(3,maxat)
      dimension ict(natom,*),df(3,*),ptr(3,*)
c
c     zero = 0.0d+00
c     one = 1.0d+00
c
c      write(iwr,9999) ic
c 9999 format(/,' appel a symd1  ic=',i3)
c      write(iwr,9994)
c 9994 format(/' gradients without pseudopotential contributions'/)
c      write(iwr,9998) ((de(i,j),j=1,natom),i=1,3)
c 9998 format(5x,5f12.8)
c      write(iwr,9993)
c 9993 format(/' unsymmetrized pseudopot. contibutions to gradient'/)
c      write(iwr,9998) ((df(i,j),j=1,natom),i=1,3)
c
      do 30 i = 1 , 3
         do 20 k = 1 , natom
            de(i,k) = de(i,k) + df(i,k)
 20      continue
 30   continue
c
      if (nt.eq.1) return
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
c
      do 70 it = 2 , nt
         icnu = ict(ic,it)
         if (icnu.lt.ic) then
            itm = it - 1
            do 40 it1 = 1 , itm
               if (ict(ic,it1).eq.icnu) go to 70
 40         continue
            n = 3*(invt(it)-1)
            do 60 j = 1 , natom
               j1 = ict(j,it)
               do 50 k = 1 , 3
                  de(k,j1) = de(k,j1) + df(1,j)*ptr(1,n+k) + df(2,j)
     +                       *ptr(2,n+k) + df(3,j)*ptr(3,n+k)
 50            continue
 60         continue
         end if
c
c      write(iwr,9996)
c      write(iwr,9998) ((de(i,j),i=1,3),j=1,natom)
c 9996 format(/' total symmetrized gradients'/)
c
 70   continue
      return
      end
**==wder.f
      subroutine wder(dsdw,dsdx,dsdy,dsdz,dd,dumw,dumx,dumy,
     +                dumz,db,ict,scftyp,natom,l2)
      implicit real*8  (a-h,o-z)
      logical out,norm,iskip
      character*8 scftyp,blank
      character*8 zlib,ztagp
      integer maxorb, maxat,  maxvar, maxnz,  mxshel, mxprim
      integer mxgrps, mxgaus, mxgrid, mxcalc, mxplot, mxrest
      integer mxstp,  maxlfn, maxfrt, maxbuf, maxblo, mxprms
      integer nd200,  mxcan1, mxcan2, lenci,  numspl, nbrkmx
      integer mxcsf,  mxnshl, mxroot, mxconf, maxig, mxtrm
      integer mxtda1, mxtda2, mxtda3, mxtda4, mxorb3, maxat3
      integer mxcrec, mxcrc2, mxproc
      integer mcprim, mcfzc
************************************************************************
*   ==========  parameters defining the maximum system size =========
*
*   there are eight    parameters that the programmer need set:
*     maxorb = maximum number of basis functions
*     maxat  = maximum number of atoms (including point charges)
*     maxvar = maximum number of z-matrix variables
*     maxnz  = maximum number of z-matrix cards
*     mxshel = maximum number of shells
*     mxprim = maximum number of shell primitives
*     mxprms = maximum number of primitives in a shell
*
      parameter (maxorb= 4096, maxat=750)
      parameter (maxvar= 2000, maxnz=700)
      parameter (mxshel= 2048, mxprim=8192, mxprms=50)
*
*   following parameters refer to analysis modules
*     mxgaus = maximum number of orbital primitives
*     mxgrps = maximum number of shells
      parameter ( mxgrps = 560, mxgaus = 11600)
*
*   following parameters refer to graphics module
      parameter (mxgrid=10, mxcalc=10, mxplot=10, mxrest=10)
      parameter (mxstp=mxcalc+mxgrid+mxplot+mxrest)
*
*   following parameters refer to I/O system
*   parameters control no. of ed/mt files+ buffers
*
*     maxlfn *  no. of ed/mt streams
*     maxfrt *  no. of fortran data sets
*     maxbuf *  no. of fortran store buffers
*     maxblo *  no. of blocks in 1 buffer
*
      parameter (maxlfn = 40, maxfrt = 60)
      parameter (maxbuf = 9, maxblo=32)
*
*   following parameters refer to direct-CI module
*   max # (external) orbitals
      parameter (nd200 = 255)
*   parameters control canonical set size
*     mxcan1 *  default setting 2508 : high-spin 19606
*     mxcan2 *  default setting 5016 : high-spin 39212
*     parameter (mxcan1  = 2508, mxcan2 = 5016)
      parameter (mxcan1  = 19606, mxcan2 = 39212)
*
*   following parameters refer to full-CI module
*
*     lenci *  default setting 500000
      parameter (lenci = 500000)
*
*   following parameters are needed for DIRECT
*
      parameter (numspl=50)
      parameter (nbrkmx=20)
*
*   following parameters are needed for MRD-CI
*
      parameter (mxcsf=100)
      parameter (mxnshl=30)
      parameter (mxroot=50)
      parameter (mxconf=200000)
      parameter (maxig=400000)
      parameter (mxtrm=600000)
      parameter (mxcrec=2000,mxcrc2=1000)
*
* following parameters are used in the TDA module
*
      parameter (mxtda1=3600)
      parameter (mxtda2=50)
      parameter (mxtda3=20)
      parameter (mxtda4=600)
c
c following parameters are used in the MCSCF
c
      parameter (mcprim=128) ! the max. number of active orbitals
      parameter (mcfzc =512) ! the max. number of frozen core orbitals

************************************************************************
*  
*   for parallel code
*
************************************************************************
      parameter (mxproc=512)

************************************************************************
*
*   the following values should not be altered
*
************************************************************************
      parameter (mxorb3=maxorb*3)
      parameter (maxat3=maxat+3)
************************************************************************
c
c
      character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
      character *10 zbflab
      character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud
      common /runlab/ zcom(19),ztitle(10),zaname(maxat),zbflab(maxorb),
     +   ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,zconf,zstate,
     +   zorb(maxorb),zpseud(maxat)
c
c....  zaname : names of atoms during calculation (after reorder)
c....  ztag   : names of atoms as read in 
c
      common/junkc/zlib(204),ztagp(10)
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      integer nprint, itol, icut, normf, normp, nopk, irest
      integer nrec, ist, jst, kst, lst 
      integer nintmx, nindmx, intg76
      integer mfilep, mainp, mblp, iblkmp
      integer m2file, m2tape, m2blk, m2last
      integer m4file, m4tape, m4blk, m4last
      integer m6file, m6tape, m6blk, m6last
      integer m5file, m5tape, m5blk, m5last
      integer m9file, m9tape, m9blk, m9last
      integer mtfile, mttape, mtblk, mtlast
      integer m1file, m1tape, m1blk, m1last
      integer m11fil, m11tap, m11bl, m11lst
      integer m12fil, m12tap, m12bl, m12lst
      integer m13fil, m13tap, m13bl, m13lst
      integer local, mtask 
      integer itask, irest2, irest3, irest4, irest5, intloc
      integer iblkl, ifill, iblkd, ifild, iblks, ifils, iblkf, ifockf 
      integer nopkr, iofsym, iofrst, idurie, imaxb_ic
      logical omaxb, ognore
      common/restar/nprint,itol,icut,normf,normp,nopk,
     + irest,nrec,omaxb,ist,jst,kst,lst,nintmx,nindmx,intg76,
     + mfilep,mainp,mblp,iblkmp,
     + m2file,m2tape(20),m2blk(20),m2last(20),
     + m4file,m4tape(20),m4blk(20),m4last(20),
     + m6file,m6tape(20),m6blk(20),m6last(20),
     + m5file,m5tape(20),m5blk(20),m5last(20),
     + m9file,m9tape(20),m9blk(20),m9last(20),
     + mtfile,mttape(20),mtblk(20),mtlast(20),
     + m1file,m1tape(20),m1blk(20),m1last(20),
     + m11fil,m11tap(20),m11bl(20),m11lst(20),
     + m12fil,m12tap(20),m12bl(20),m12lst(20),
     + m13fil,m13tap(20),m13bl(20),m13lst(20),
     + local,mtask,itask(50),
     + irest2,irest3,irest4,irest5,intloc,
     + iblkl,ifill,iblkd,ifild,iblks,ifils,iblkf,ifockf,
     + nopkr,iofsym,iofrst,idurie(2),ognore,imaxb_ic
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common /rcombd/ win(36),pxin(36),pyin(36),pzin(36)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c  NB atomic masses now accessed through amass_get
c     array retained here as a placeholder as there
c     are explicit non-included /infoa/ commons in the code
c
      real*8 czan, c, amasold, symz
      integer nat, num, ich, mul, nx, ne, na, nb, imass
      integer nuct, ipseud, lpseud
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,czan(maxat),c(3,maxat)
     +            ,amasold(maxat),
     +             imass(maxat),nuct(maxat),ipseud(maxat),
     +             symz(maxat),lpseud
c
       integer len_infoa
       parameter (len_infoa=8)
c      used: restre(util1),revise(util1),utyp21(server) (first 8)
c
      real*8 ex, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm
      integer isptype
      common /nshel/ ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     +               cf(mxprim),cg(mxprim),
     +               kstart(mxshel),katom(mxshel),ktype(mxshel),
     +               kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
     +               nshell,non,numorb,ndumm,isptype
c
      common/junk/desp(3*maxat),
     +  xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,cx,cy,cz
     + ,chglib(204),maxlib(204),ibllib(204),
     +  apsl(20,4),cpsl(20,4),cpslo(400,4),
     +  libmax(4),libort(4),
     +  s(114),dij(114),ijg(114),ijx(114),ijy(114),ijz(114)
     + ,icent(maxorb),de(3,maxat),ppptr(3,144)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      dimension iskip(20),indx(40)
      dimension ix(20),iy(20),iz(20),jx(10),jy(10),jz(10)
      dimension min(5),max(5),maz(5)
      character *3 dnam
      dimension dnam(3)
c
      dimension dumx(9,*),dumy(9,*),dumz(9,*),dumw(9,*)
      dimension dsdw(*),dsdx(*),dsdy(*),dsdz(*)
      dimension dd(l2),db(l2)
      dimension ict(natom,*)
c
      data sqrt3/1.73205080756888d0/,pi32/5.5683279968317d0/
      data dnam /'e*x','e*y','e*z'/
      data blank/'        '/
      data indx / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     1           -0,-0,-0,-0,-0,-0,-0,-0,-0,-0,
     2           -0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
     3           10,11,12,13,14,15,16,17,18,19/
c     2           -0, 1, 2, 3,-0,-0,-0,-0,-0,-0,
c     3            4, 5, 6, 7, 8, 9,10,11,12,13/
      data jx / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0/
      data ix / 1, 4, 1, 1, 7, 1, 1, 4, 4, 1,
     1         10, 1, 1, 7, 7, 4, 1, 4, 1, 4/
      data jy / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1/
      data iy / 1, 1, 4, 1, 1, 7, 1, 4, 1, 4,
     1          1,10, 1, 4, 1, 7, 7, 1, 4, 4/
      data jz / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1/
      data iz / 1, 1, 1, 4, 1, 1, 7, 1, 4, 4,
     1          1, 1,10, 1, 4, 1, 4, 7, 7, 4/
      data zero,one /0.0d+00,1.0d+00/
      data min /1,2,5,11,21/
      data max /1,4,10,20,35/
      data maz /1,4,9,17,29/
      data m10/10/,m204/204/
c
      call dendd1(scftyp,dd,db,l2)
c
c     ----- read ict from file 
c
      nav = lenwrd()
      call readi(ict,nw196(6)*nav,ibl196(6),idaf)
c
c retrieve header blocks from pseudopotential library
c
      mword = 1760 + 8/nav
      m408 = 204 + 408/nav
      call rdchr(zlib,m204,iblkpl,numlib)
      call reads(chglib,m408,numlib)
c
c     ----- calculate derivatives of the nonlocal potentiel--
c
      tol = 2.30258d+00*itol
      out = nprint.eq. - 3
      norm = normf.ne.1 .or. normp.ne.1
c
c     ----- centres------------
c
      do 620 iat = 1 , natom
c
         if (zpseud(iat).ne.blank) then
            do 60 it = 1 , nt
               icd = ict(iat,it)
               if (icd.gt.iat) go to 620
 60         continue
            do 80 i = 1 , 3
               do 70 k = 1 , natom
                  de(i,k) = zero
 70            continue
 80         continue
c
            xi = c(1,iat)
            yi = c(2,iat)
            zi = c(3,iat)
c
            numps = locatc(zlib,204,zpseud(iat))
            if (numps.eq.0) then
               call caserr2(
     +     'attempt to retrieve unknown pseudopotential from library')
            end if
 90         continue
c           chg = chglib(numps)
            npmax = maxlib(numps)
            call rdchr(ztagp,m10,ibllib(numps)+iblkpl,numlib)
            call reads(apsl,mword,numlib)
            if (npmax.le.0 .or. npmax.gt.3) then
               call caserr2(
     +          'invalid parameters from pseudopotential library')
            else
               do 610 np = 1 , npmax
                  itpmax = libmax(np)
                  itport = libort(np)
c
                  do 150 ig = 1 , itpmax
                     ax = apsl(ig,np)
                     ax = ax + ax
                     go to (100,110,120) , np
 100                 cnorm = pi32/(dsqrt(ax)*ax)
                     go to 130
 110                 cnorm = pi32*0.5d0/(dsqrt(ax)*ax*ax)
                     go to 130
 120                 cnorm = pi32*0.75d0/(dsqrt(ax)*ax**3)
 130                 cnorm = 1.d0/dsqrt(cnorm)
                     do 140 iip = 1 , itport
                        ipp = (iip-1)*itpmax + ig
                        cpslo(ipp,np) = cpslo(ipp,np)*cnorm
 140                 continue
 150              continue
                  lit = np + 1
                  mini = min(np)
                  maxi = max(np)
                  mazi = maz(np)
                  ncomp = mazi - mini + 1
c                 nds = ncomp*num
                  ndsd = 0
                  do 160 i = 1 , 20
                     iskip(i) = .true.
 160              continue
                  do 230 i = mini , maxi
                     go to (170,190,230,230,210,230,230,230,230,230) , i
 170                 do 180 k = 1 , 4
                        iskip(k) = .false.
 180                 continue
                     go to 230
 190                 do 200 k = 1 , 10
                        iskip(k) = .false.
 200                 continue
                     go to 230
 210                 do 220 k = 2 , 20
                        iskip(k) = .false.
 220                 continue
 230              continue
c
c-----------------primitives of the nonlocal potential----------
c
                  do 510 ig = 1 , itpmax
                     ai = apsl(ig,np)
                     axi = ai*xi
                     ayi = ai*yi
                     azi = ai*zi
                     go to (240,250,260,270,270) , np
 240                 csi = one
                     cpi = ai + ai
                     go to 270
 250                 csi = one
                     cpi = one
                     cdi = ai + ai
                     go to 270
 260                 cpi = one
                     cdi = one
                     cfi = ai + ai
c
c     ----- jshell
c
 270                 do 500 jj = 1 , nshell
                        jat = katom(jj)
                        xj = c(1,jat)
                        yj = c(2,jat)
                        zj = c(3,jat)
                        j1 = kstart(jj)
                        j2 = j1 + kng(jj) - 1
                        ljt = ktype(jj)
                        minj = kmin(jj)
                        maxj = kmax(jj)
                        jdif = maxj - minj + 1
                        locj = kloc(jj) - minj
c                       nroots = (lit+ljt)/2
                        rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        arri = ai*rr
c
c     ----- prepare indices for pairs of (i,j) functions
c
                        n0 = 0
                        if (lit.eq.4) n0 = 20
                        ij = 0
                        do 290 i = 1 , 20
                           if (.not.(iskip(i))) then
                              in = indx(i+n0)
                              nnx = ix(i)
                              nny = iy(i)
                              nnz = iz(i)
                              do 280 j = minj , maxj
                                 ij = ij + 1
                                 ijx(ij) = nnx + jx(j)
                                 ijy(ij) = nny + jy(j)
                                 ijz(ij) = nnz + jz(j)
                                 ijg(ij) = in + 19*(j-minj)
 280                          continue
                           end if
 290                    continue
                        do 300 i = 1 , ij
                           n = ijg(i)
                           s(n) = zero
 300                    continue
c
c     ----- j primitive
c
                        do 460 jg = j1 , j2
                           aj = ex(jg)
                           aa = ai + aj
                           dum = aj*arri/aa
                           if (dum.le.tol) then
                              fac = dexp(-dum)
                              csj = cs(jg)
                              cpj = cp(jg)
                              cdj = cd(jg)
                              ax = (axi+aj*xj)/aa
                              ay = (ayi+aj*yj)/aa
                              az = (azi+aj*zj)/aa
c
c     ----- density factor
c
                              nn = 0
                              do 420 i = 1 , 20
                                 if (iskip(i)) go to 420
                                 go to (310,320,350,350,330,350,350,350,
     +                                  350,350,340,350,350,350,350,350,
     +                                  350,350,350,350) , i
 310                             dum1 = csi*fac
                                 go to 350
 320                             dum1 = cpi*fac
                                 go to 350
 330                             dum1 = cdi*fac
                                 go to 350
 340                             dum1 = cfi*fac
 350                             do 410 j = minj , maxj
                                    go to (360,370,400,400,380,400,400,
     +                                 390,400,400) , j
 360                                dum2 = dum1*csj
                                    go to 400
 370                                dum2 = dum1*cpj
                                    go to 400
 380                                dum2 = dum1*cdj
                                    go to 400
 390                                if (norm) dum2 = dum2*sqrt3
 400                                nn = nn + 1
                                    dij(nn) = dum2
 410                             continue
 420                          continue
c
c     ----- overlap
c
                              aa1 = one/aa
                              t = dsqrt(aa1)
                              x0 = ax
                              y0 = ay
                              z0 = az
                              in = -3
                              do 440 i = 1 , lit
                                 in = in + 3
                                 ni = i
                                 do 430 j = 1 , ljt
                                    jn = in + j
                                    nj = j
                                    call vint
                                    pxin(jn) = xint*t
                                    pyin(jn) = yint*t
                                    pzin(jn) = zint*t
 430                             continue
 440                          continue
                              do 450 i = 1 , ij
                               n = ijg(i)
                               nnx = ijx(i)
                               nny = ijy(i)
                               nnz = ijz(i)
                               s(n) = s(n) + dij(i)*pxin(nnx)*pyin(nny)
     +                                *pzin(nnz)
 450                          continue
                           end if
c
c     ----- end of primitive loops -----
c
c debug writes  --- values of s ---
c      write(6,*)' end of primitive loops -- contents of s'
c      write(6,*) s
c
 460                    continue
c
c     ----- form integrals over derivatives -----
c
                        nn = 0
                        n = 1
                        do 470 j = minj , maxj
                           icent(locj+j) = jat
                           if (mini.le.1) then
                              nn = nn + 1
                              win(nn) = s(n)
                              pxin(nn) = s(n+1)
                              pyin(nn) = s(n+2)
                              pzin(nn) = s(n+3)
                              if (maxi.eq.1) then
                                 n = n + 19
                                 go to 470
                              end if
                           end if
                           if (mini.le.2) then
                              nn = nn + 1
                              win(nn) = s(n+1)
                              pxin(nn) = (s(n+4)-s(n))
                              pyin(nn) = s(n+7)
                              pzin(nn) = s(n+8)
                              nn = nn + 1
                              win(nn) = s(n+2)
                              pxin(nn) = s(n+7)
                              pyin(nn) = (s(n+5)-s(n))
                              pzin(nn) = s(n+9)
                              nn = nn + 1
                              win(nn) = s(n+3)
                              pxin(nn) = s(n+8)
                              pyin(nn) = s(n+9)
                              pzin(nn) = (s(n+6)-s(n))
                              if (maxi.eq.4) then
                                 n = n + 19
                                 go to 470
                              end if
                           end if
                           nn = nn + 1
                           win(nn) = s(n+3)
                           pxin(nn) = (s(n+9)-s(n)-s(n))
                           pyin(nn) = s(n+12)
                           pzin(nn) = s(n+13)
                           nn = nn + 1
                           win(nn) = s(n+4)
                           pxin(nn) = s(n+14)
                           pyin(nn) = (s(n+10)-s(n+1)-s(n+1))
                           pzin(nn) = s(n+15)
                           nn = nn + 1
                           win(nn) = s(n+5)
                           pxin(nn) = s(n+16)
                           pyin(nn) = s(n+17)
                           pzin(nn) = (s(n+11)-s(n+2)-s(n+2))
                           nn = nn + 1
                           dum = one
                           if (norm) dum = sqrt3
                           win(nn) = s(n+6)*dum
                           pxin(nn) = dum*(s(n+12)-s(n+1))
                           pyin(nn) = dum*(s(n+14)-s(n))
                           pzin(nn) = dum*s(n+18)
                           nn = nn + 1
                           win(nn) = s(n+7)*dum
                           pxin(nn) = dum*(s(n+13)-s(n+2))
                           pyin(nn) = dum*s(n+18)
                           pzin(nn) = dum*(s(n+16)-s(n))
                           nn = nn + 1
                           win(nn) = s(n+8)*dum
                           pxin(nn) = dum*s(n+18)
                           pyin(nn) = dum*(s(n+15)-s(n+2))
                           pzin(nn) = dum*(s(n+17)-s(n+1))
                           n = n + 19
 470                    continue
c
                        if (np.ge.3) call combd(np,jdif)
c
c
                        n = 0
                        do 490 j = minj , maxj
                           do 480 i = 1 , ncomp
                              n = n + 1
                              ndsd = ndsd + 1
                              dsdw(ndsd) = win(n)
                              dsdx(ndsd) = pxin(n)
                              dsdy(ndsd) = pyin(n)
                              dsdz(ndsd) = pzin(n)
 480                       continue
 490                    continue
 500                 continue
 510              continue
c
c
c  recombine non local potential integrals over primitives-----
c
                  ndsd = 0
c
                  do 600 i = 1 , itport
                     do 530 j = 1 , num
                        do 520 ig = 1 , ncomp
                           dumw(ig,j) = 0.d0
                           dumx(ig,j) = 0.d0
                           dumy(ig,j) = 0.d0
                           dumz(ig,j) = 0.d0
 520                    continue
 530                 continue
                     cps = cpsl(i,np)
c
                     ndsd = 0
                     ipp = (i-1)*itpmax
                     do 560 j = 1 , itpmax
c
                        cx = cpslo(ipp+j,np)
                        do 550 k = 1 , num
                           do 540 ig = 1 , ncomp
                              ndsd = ndsd + 1
                              dumw(ig,k) = dsdw(ndsd)*cx + dumw(ig,k)
                              dumx(ig,k) = dsdx(ndsd)*cx + dumx(ig,k)
                              dumy(ig,k) = dsdy(ndsd)*cx + dumy(ig,k)
                              dumz(ig,k) = dsdz(ndsd)*cx + dumz(ig,k)
 540                       continue
 550                    continue
 560                 continue
c
c--------form the 1.e gradient------------------------
c
                     do 590 j = 1 , num
                        jat = icent(j)
                        do 580 k = 1 , j
                           kat = icent(k)
                           jk = iky(j) + k
                           ddjk = dd(jk)*cps
                           if (j.ne.k) ddjk = ddjk + ddjk
                           dummx = 0.d0
                           dummy = 0.d0
                           dummw = 0.d0
                           dummz = 0.d0
                           dummp = 0.d0
                           dummq = 0.d0
                           dummr = 0.d0
                           do 570 ig = 1 , ncomp
                              dummw = dumw(ig,j)*dumw(ig,k) + dummw
                              dummx = dumx(ig,j)*dumw(ig,k) + dummx
                              dummy = dumy(ig,j)*dumw(ig,k) + dummy
                              dummz = dumz(ig,j)*dumw(ig,k) + dummz
                              dummp = dumw(ig,j)*dumx(ig,k) + dummp
                              dummq = dumw(ig,j)*dumy(ig,k) + dummq
                              dummr = dumw(ig,j)*dumz(ig,k) + dummr
 570                       continue
                           de(1,iat) = de(1,iat) + (dummx+dummp)*ddjk
                           de(2,iat) = de(2,iat) + (dummy+dummq)*ddjk
                           de(3,iat) = de(3,iat) + (dummz+dummr)*ddjk
                           de(1,jat) = de(1,jat) - dummx*ddjk
                           de(2,jat) = de(2,jat) - dummy*ddjk
                           de(3,jat) = de(3,jat) - dummz*ddjk
                           de(1,kat) = de(1,kat) - dummp*ddjk
                           de(2,kat) = de(2,kat) - dummq*ddjk
                           de(3,kat) = de(3,kat) - dummr*ddjk
 580                    continue
 590                 continue
 600              continue
c
c debug write statements -- de at end of each shell ---
c      write(6,*)'de(i,j),i=1,3,j=1,natom'
c      write(6,*)' np = ',np
c      write(6,*)((de(i,j),i=1,3),j=1,natom)
c
c     --- end of pseudopotential shell loop ---
 610           continue
               call symd1(de,ppptr,iat,ict,natom)
            end if
         end if
c     --- end of loop over atoms ---
 620  continue
c
      if (out) then
         write (iwr,6010)
         write (iwr,6020)
         write (iwr,6010)
         write (iwr,6030) (i,i=1,natom)
         do 630 n = 1 , 3
            write (iwr,6040) dnam(n) , (de(n,i),i=1,natom)
 630     continue
      end if
      call timit(0)
      return
 6010 format (/)
 6020 format (20x,'pseudopotential contribution to gradient')
 6030 format (5x,'atom',8(6x,i3,6x))
 6040 format (7x,a3,8e15.7)
      end
**==comb.f
      subroutine comb(s,itpmax,lkt,jd1)
      implicit real*8  (a-h,o-z)
      dimension md1(4),index(10),t(10),s(225,*)
      data rt3ov2,half/0.8660254037844d0,0.5d0/
      data md1/1,3,6,10/
      kd1 = md1(lkt)
c     jk1 = jd1*kd1
      rt5 = dsqrt(5.0d0)
      rt3o40 = dsqrt(3.0d0/40.d0)
      f3o2r5 = 3.0d0/(2.0d0*rt5)
      rt5o8 = dsqrt(5.0d0/8.0d0)
      rt9o8 = dsqrt(9.0d0/8.0d0)
c
      if (lkt.lt.3) then
      else if (lkt.eq.3) then
c                 p   d   f
         do 40 ig = 1 , itpmax
            do 30 j1 = 1 , jd1
               do 20 loop = 1 , kd1
                  index(loop) = j1 + jd1*(loop-1)
                  t(loop) = s(index(loop),ig)
 20            continue
               szz = t(3) - half*(t(1)+t(2))
               sxz = t(5)
               syz = t(6)
               sx2y2 = rt3ov2*(t(1)-t(2))
               sxy = t(4)
               s(index(1),ig) = sx2y2
               s(index(2),ig) = szz
               s(index(3),ig) = sxy
               s(index(4),ig) = sxz
               s(index(5),ig) = syz
 30         continue
 40      continue
      else
         do 70 ig = 1 , itpmax
            do 60 j1 = 1 , jd1
               do 50 loop = 1 , kd1
                  index(loop) = j1 + jd1*(loop-1)
                  t(loop) = s(index(loop),ig)
 50            continue
               szzz = t(3) - f3o2r5*(t(5)+t(7))
               sxzz = rt3o40*(4.0d0*t(8)-rt5*t(1)-t(6))
               syzz = rt3o40*(4.0d0*t(9)-t(4)-rt5*t(2))
               sx2y2z = rt3ov2*(t(5)-t(7))
               sxyz = t(10)
               sxxx = rt5o8*t(1) - rt9o8*t(6)
               syyy = rt5o8*t(2) - rt9o8*t(4)
               s(index(1),ig) = szzz
               s(index(2),ig) = sxzz
               s(index(3),ig) = syzz
               s(index(4),ig) = sx2y2z
               s(index(5),ig) = sxyz
               s(index(6),ig) = sxxx
               s(index(7),ig) = syyy
 60         continue
 70      continue
      end if
      return
      end
**==combd.f
      subroutine combd(lkt,jd1)
      implicit real*8  (a-h,o-z)
      common /rcombd/ win(36),pxin(36),pyin(36),pzin(36)
      dimension md1(4),md2(4)
      dimension t(100)
      dimension s(36,4)
      equivalence(s(1,1),win(1))
      data c1,c2/0.8660254037844d0,0.5d0/
      data md1/1,3,6,10/,md2/1,3,5,7/
      kd1 = md1(lkt)
      kd2 = md2(lkt)
      jk1 = jd1*kd1
c
      do 40 ig = 1 , 4
         call dcopy(jk1,s(1,ig),1,t,1)
         if (lkt.lt.3) go to 40
         if (lkt.ne.3) then
            call caserr2('f-pseudopotential invalid')
         end if
         do 30 j = 1 , jd1
            j1 = (j-1)*kd1 + 1
            j2 = (j-1)*kd2 + 1
            t1 = t(j1)
            t2 = t(j1+1)
            t3 = t(j1+2)
            s(j2,ig) = c1*(t1-t2)
            s(j2+1,ig) = t3 - c2*(t1+t2)
            do 20 l = 3 , 5
               s(j2+l-1,ig) = t(j1+l)
 20         continue
 30      continue
 40   continue
      return
      end
      subroutine ver_integb_nl(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/integb_nl.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
