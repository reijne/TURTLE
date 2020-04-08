
c
c  Two-electron integral driver for CCP1 DFT code
c
c ********
c *** these integral routines have been extracted from GAMESS-UK,
c *** and modified to work as far as possible independently of the 
c *** parent code.
c *** As they now stand, they function (i) only in direct-mode.
c *** (ii) assuming that symmetry has been disabled (i.e must be run 
c *** with nosym), and supermatrices have also been disabled.
c *******
c
c  Subroutines in this file:
c
c
c   subroutine jkint_dft(iso,gout,nshels,
c     &     basi, basj, bask, basl,
c     &     imode,aaa,bbb,ccc,ddd, schwarz_ao, schwarz_cd)
c
c   subroutine jkint_dft_cntgen(iso,gout,nshels,
c     &     basi, basj, bask, basl,
c     &     imode,aaa,bbb,ccc,ddd, schwarz_ao, schwarz_cd,
c     &     te3c_int, nte3c_int, ite3c_stored, nte3c_shl,
c     &     te2c_int, nte2c_int, ite2c_stored, nte2c_shl)
c
c   subroutine jkint_dft_genuse(iso,gout,nshels,
c     &     basi, basj, bask, basl,
c     &     imode,aaa,bbb,ccc,ddd, 
c     &     te3c_int, nte3c_int, ite3c_stored, nte3c_shl,
c     &     te2c_int, nte2c_int, ite2c_stored, nte2c_shl)
c
c   subroutine ijprim_dft(ncentr)
c   subroutine ijss_dft(gout,ncent)
c   subroutine indxa_dft(ijx,ijy,ijz,ij,mini,maxi,
c     &     minj,maxj,iandj,inc1,inc2,inc3)
c   subroutine shells_dft(gout,nelec,
c     &     ish,jsh,ksh,lsh,
c     &     basi, basj, bask, basl, ncentr,
c     &     iexch,nintegral)
c   subroutine coulmb_dft
c   subroutine spchck_dft
c   subroutine spdint_dft(gout)
c   subroutine ssdss_dft(qq)
c   subroutine sskl_dft(gout,ncentr)
c   subroutine ssprim_dft
c   subroutine genral_dft(gout,ncentr)
c   subroutine s0000_dft(g)
c   subroutine xyzint_dft
c   subroutine roots4_dft
c   subroutine roots5_dft
c   subroutine rootss_dft
c   subroutine rt123_dft
c   subroutine ffun_dft(x,f,npt,dji,madd)
c   subroutine final_dft
c   subroutine aux_find(tag)
c   integer function num_3c_ints()
c   subroutine ver_dft_integ2e(s,r,d)
c
c**********************************************************************

c Arguments:
c
c gout - workspace for integrals
c
c     ----- size of gout -
c                         1   if s or k shells
c                        81   if p      shells
c                       256   if      l shells
c                      1296   if d or m shells
c                     10000   if f shells
c                     50625   if g shells
c
c     ----- this version can handle g shells   -----
c
c
c basi  (etc) Basis sets for the 4 indices, -1 for 
c         a redundant index.
c
c imode - how to process integrals
c
c       1   - write them to disk (GAMESS-UK ed2 format) (demolished)
c
c       2   - regular 4-centre fock build
c             contract with density in a,b store in c,d
c
c       3   - 3-centre fock build 
c             contract with fitting coeffients, a,b store in c,d
c
c       4   - 3-centre density expansion
c             contract with density in a, store Tr in c
c
c       5   - 2-centre normalisation compute and store only
c             diagonal integrals (defeat k loop)
c
c       6   - 2-centre integral evaluation (store raw integrals
c             in a)
c
c       61  - 2-centre integral evaluation (store raw integrals
c             in global array)
c
c a,b,c,d   - input/output matrices
c 
c nshels    - leading dimension if iso
c iso       - symmetry array

***************************************************************
* EXTERNAL routines invoked
**************************************************************
*
* Parallel load balancing etc
*
*     iipsci, oipsci
*
*
* Memory allocation routines
*
*     setscm, setc, cmem, loccm
*
* BLAS - linear algebra
*
*      daxpy, vlcr, vadd, setsto
*
* ERROR monitoring
*
*      caserr

***************************************************************
***************************************************************
* common blocks used from main-line code:
*
*   cslosc * files * infoa *  parallel * parcntl 
*   prints * restar * sizes * statis * timez
*
********************************************************************

      subroutine coulmb_dft(clints,gout,ibas)
c
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension clints(*),gout(*)
      dimension ibpop(4,4),ibr(4,4)
c
cINCLUDE(common/flip70)
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
c
      integer maxfitorb
c      parameter(maxfitorb=3000)

      parameter(maxfitorb=maxorb)

      integer iky,ilifq
      common/mapperx/iky(maxfitorb),ilifq(maxfitorb)
c
cINCLUDE(common/restri)
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      integer iwr
      common/dft_iofile/iwr
cINCLUDE(common/nshel)
cINCLUDE(common/shlg70)
c
      real*8 pito52, pidiv4, root3, root5, root53, root7
      common/piconx/pito52,pidiv4,root3,root5,root53,root7
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
cINCLUDE(common/cslosc)
c
      data m25/25/
c     data m22,mword1/22,6000/
c
      data ibpop/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
      data ibr/64,16,4,1,216,36,6,1,1000,100,10,1,
     +         3375,225,15,1 /
      data done/1.0d0/,two/2.0d0/,twopt5/2.5d0/,four/4.0d0/
      data dzero/0.0d0/
c
c     ----- compute all coulomb integrals (ii,jj//ii,jj) -----
c     the largest coulomb integral from each integral block
c     is picked out and saved, so that the schwarz inequality
c     can be used later to avoid entire integral blocks.
c
      call spchck_dft
      outi = nprint.ne.-5
      nclint = (nshell(ibas)+1)*nshell(ibas)/2
c     lensch = lensec(nclint)
c     tim0 = cpulft(1)
c
c === save tolerances set in debut and increase accuracy
c
      tolsv = tol
      tol = 75.0d00
      cutsv = cutoff
      cutoff = 1.0d0/(10.0d0**20)
c
      qq4   = done
      nint  = 0
      call aclear_dp(clints,nclint,0.0d0)
c     if (intg76.ne.0) then
c
c     ----- loop over all shell blocks -----
c     ----- first over rotated axis integrals
c
c        pidiv4 = datan(done)
c        pi = four*pidiv4
c        pito52 = two*pi**twopt5
c        call sinset
c        call filmax
c        do 60 ishell = 1, nshell
c           iiii = iky(ishell)
c           if (kad(ishell).ge.0) then
c              do 50 jshell = 1, ishell
c                 ijij = iiii + jshell
c                 if (kad(jshell).ge.0) then
c
c     use pople code for any pure sp integral blocks,
c     use hondo rys polynomial code for other blocks
c
c                    kshell = ishell
c                    lshell = jshell
c                    call genr70(gout,1,.false.)
c
c     ----- pick out largest coulomb integral for this block -----
c
c                    vmax = dzero
c                    mini = kmin(ishell)
c                    minj = kmin(jshell)
c                    maxi = kmax(ishell)
c                    jmax = kmax(jshell)
c                    oiandj = ishell.eq.jshell
c                    ibb = ib(1,1)
c                    jbb = ib(2,1)
c                    kbb = ib(3,1)
c                    lbb = ib(4,1)
c                    ijn = 0
c                    do 40 i = mini, maxi
c                       if (oiandj) jmax = i
c                       do 30 j = minj, jmax
c                          nn = ibpop(ibb,i) + ibpop(jbb,j)
c    +                        + ibpop(kbb,i) + ibpop(lbb,j) + 1
c                          val = gout(nn)
c                          if (val.gt.dzero) nint = nint+1
c                          if (val.gt.vmax) vmax = val
c30                     continue
c40                  continue
c                    clints(ijij) = vmax
c                 end if
c50            continue
c           end if
c60      continue
c        if (ospbas) go to 110
c     end if
c
      ii = 1
      if (opdbas) then
       ii = 2
      end if
      if (opfbas) then
       ii = 3
      end if
      if (opgbas) then
       ii = 4
      end if
      do 300 loop = 1 , 3
        igt(loop) = ibr(1,ii)
        jgt(loop) = ibr(2,ii)
        kgt(loop) = ibr(3,ii)
        lgt(loop) = ibr(4,ii)
 300  continue
      qq4   = done
c
      do 100 ishell = 1, nshell(ibas)
         iiii = iky(ishell)
         kadi = kad(ibas,ishell)
         do 90 jshell = 1, ishell
            ijij = iiii + jshell
            if ((kadi+kad(ibas,jshell)).lt.0) then
c
c     use hondo rys polynomial code for other blocks
c
               call shells_dft(gout,1,ishell,jshell,ishell,jshell,
     &                         ibas,ibas,ibas,ibas,4,1,idum)
               call ijprim_dft(4)
               if (nij.ne.0) then 
                  call shells_dft(gout,2,ishell,jshell,ishell,jshell,
     &                            ibas,ibas,ibas,ibas,4,1,idum)
                  call genral_dft(gout,4)
c
c     ----- pick out largest coulomb integral for this block -----
c
                  vmax = dzero
                  mini = kmin(ibas,ishell)
                  minj = kmin(ibas,jshell)
                  maxi = kmax(ibas,ishell)
                  jmax = kmax(ibas,jshell)
                  oiandj = ishell.eq.jshell
                  ijn = 0
                  do 80 i = mini, maxi
                     if (oiandj) jmax = i
                     do 70 j = minj, jmax
                        ijn = ijn+1
                        nn = ijgt(ijn) + klgt(ijn)
                        val = gout(nn)
                        if (val.gt.dzero) nint = nint+1
                        if (val.gt.vmax) vmax = val
 70                  continue
 80               continue
                  clints(ijij) = vmax
               end if
            end if
 90      continue
 100  continue
c
 110  if (outv) then
         write(iwr,6020)
         call prtri(clints,nshell(ibas))
      end if
c
      do 120 i = 1, nclint
         clints(i) = dsqrt(clints(i))
 120  continue
c
c     tcoul = cpulft(1)-tim0
c     if (outi) write(iwr,6030) nint, tcoul
c
c     call secput(isect(421),m25,lensch,iblock)
c     call wrt3(clints,nclint,iblock,idaf)
c
c === now reset tolerances
c
      tol = tolsv
      cutoff = cutsv
c
      call spchck_dft
c
      return
 6020 format (/20x,'*****************************'/20x,
     +        'max coulomb integral in shell'/20x,
     +        '*****************************'/)
 6030 format(1x,
     + 'schwarz inequality overhead:',i10,' integrals, t=',
     *       f8.2,' seconds')
      end
c
c note - onorm is hardwired to F in the 2C code
c        (ij density is a fitting function)
c
      subroutine ijprim_dft(ncentr)

      implicit none

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
c normf, normp
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
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c

      integer ncentr

      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
      real*8 cxyz, a, r, p1, q1, r1, dij
      integer ijd
      common/junkx/cxyz(3,5625),
     +     a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     +     ijd(225)
c
      real*8 sqrt3,  sqrt5, sqrt7, done
c
c Local
c
      real*8 ai, arri, axi, ayi, azi, csi, cpi, cdi, cfi, cgi 
      real*8 aj, csj, cpj, cdj, cfj, cgj 
      real*8 aa, aa1, dum, dum1, dum2
      integer i, j, ia, jb, nn, n, nm, max, jbmax
      logical onorm
c
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data done /1.0d0/
c
      if (ij.eq.1) then
c
c     ----- (s,s// -----
c
         call ssprim_dft
         return
      else
         onorm = normf.ne.1 .or. normp.ne.1
         if(ncentr.eq.2)onorm=.false.

         max = maxj
         n = 0
         nn = 0
         do 70 i = mini , maxi
            go to (20,20,30,30,
     +             20,30,30,20,30,30,
     +             20,30,30,20,30,30,30,30,30,20,
     +             20,30,30,20,30,30,30,30,30,20,
     +             30,30,20,30,30) , i
 20         nm = nn
 30         nn = nm
            if (oianj) max = i
            do 60 j = minj , max
               go to (40,40,50,50,
     +                40,50,50,40,50,50,
     +                40,50,50,40,50,50,50,50,50,40,
     +                40,50,50,40,50,50,50,50,50,40,
     +                50,50,40,50,50) , j
 40            nn = nn + 1
 50            n = n + 1
               ijd(n) = nn
 60         continue
 70      continue
c
c     ----- i primitive
c
         nij = 0

         jbmax = ngb
         do 310 ia = 1 , nga
            ai = ag(ia)
            arri = ai*rri
            axi = ai*pi
            ayi = ai*qi
            azi = ai*ri
            csi = csa(ia)
            cpi = cpa(ia)
            cdi = cda(ia)
            cfi = cfa(ia)
            cgi = cga(ia)
c
c     ----- j primitive
c
            if (oianj) jbmax = ia
            do 300 jb = 1 , jbmax
               aj = bg(jb)
               aa = ai + aj
               aa1 = done/aa
               dum = aj*arri*aa1
               if (dum.gt.tol) go to 300
               csj = csb(jb)
               cpj = cpb(jb)
               cdj = cdb(jb)
               cfj = cfb(jb)
               cgj = cgb(jb)
               nm = nij*16
               nn = nm
               nij = nij + 1
               r(nij) = dum
               a(nij) = aa
               p1(nij) = (axi+aj*pj)*aa1
               q1(nij) = (ayi+aj*qj)*aa1
               r1(nij) = (azi+aj*rj)*aa1
c
c     ----- density factor
c
               do 250 i = mini , maxi
                  go to (80,90,250,250,
     +                   100,250,250,110,250,250,
     +                   120,250,250,130,250,250,250,250,250,140,
     +                   142,250,250,144,250,250,250,250,250,146,
     +                   250,250,148,250,250) , i
c
 80               dum1 = csi*aa1
                  go to 150
 90               dum1 = cpi*aa1
                  go to 150
 100              dum1 = cdi*aa1
                  go to 150
 110              if (onorm) dum1 = dum1*sqrt3
                  go to 150
 120              dum1 = cfi*aa1
                  go to 150
 130              if (onorm) dum1 = dum1*sqrt5
                  go to 150
 140              if (onorm) dum1 = dum1*sqrt3
                  go to 150
 142              dum1 = cgi*aa1
                  go to 150
 144              if (onorm) dum1 = dum1*sqrt7
                  go to 150
 146              if (onorm) dum1 = dum1*sqrt5/sqrt3
                  go to 150
 148              if (onorm) dum1 = dum1*sqrt3
 150              if (oianj) max = i
                  do 240 j = minj , max
                     go to (160,170,240,240,
     +                      180,240,240,190,240,240,
     +                      200,240,240,210,240,240,240,240,240,220,
     +                      222,240,240,224,240,240,240,240,240,226,
     +                      240,240,228,240,240),j
 160                 dum2 = dum1*csj
                     go to 230
 170                 dum2 = dum1*cpj
                     go to 230
 180                 dum2 = dum1*cdj
                     go to 230
 190                 if (onorm) dum2 = dum2*sqrt3
                     go to 230
 200                 dum2 = dum1*cfj
                     go to 230
 210                 if (onorm) dum2 = dum2*sqrt5
                     go to 230
 220                 if (onorm) dum2 = dum2*sqrt3
                     go to 230
 222                 dum2 = dum1*cgj
                     go to 230
 224                 if (onorm) dum2 = dum2*sqrt7
                     go to 230
 226                 if (onorm) dum2 = dum2*sqrt5/sqrt3
                     go to 230
 228                 if (onorm) dum2 = dum2*sqrt3
 230                 nn = nn + 1
                     dij(nn) = dum2
 240              continue
 250           continue
               if (.not.oianj) go to 300
               if (ia.eq.jb) go to 300
               go to (290,260,280,270,275) , lit
 260           if (mini.ne.2) then
                  dij(nm+2) = dij(nm+2) + csi*cpj*aa1
                  dij(nm+3) = dij(nm+3) + dij(nm+3)
               end if
               go to 290
 275           dij(nm+10)= dij(nm+10)+ dij(nm+10)
               dij(nm+9) = dij(nm+9) + dij(nm+9)
               dij(nm+8) = dij(nm+8) + dij(nm+8)
               dij(nm+7) = dij(nm+7) + dij(nm+7)
 270           dij(nm+6) = dij(nm+6) + dij(nm+6)
               dij(nm+5) = dij(nm+5) + dij(nm+5)
               dij(nm+4) = dij(nm+4) + dij(nm+4)
 280           dij(nm+2) = dij(nm+2) + dij(nm+2)
               dij(nm+3) = dij(nm+3) + dij(nm+3)
 290           dij(nm+1) = dij(nm+1) + dij(nm+1)
 300        continue
 310     continue
         return
      end if
      end
c
c for i/j pair, therefore 3centre and 4centre codes
c are the same.
c
c GAMESS-UK dependencies
c   sizes
c   restart (okanl??)
c
      subroutine ijss_dft(gout,ncent)

      implicit none

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

      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
      
      real*8 pint, qint, rint, aa, r, p1, q1, r1, dd
      integer ijd
      common/junkx/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)

c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setintx/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                 bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                 f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                 ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denssx /dij(225),dkl(225)
c

      real*8 gout
      dimension gout(*)
      integer ncent

      integer i, n
      integer i1, i2, i3, i4, i5, kg, lg
      integer iii, nnn, jjj, m, mm, min, nn, max
      integer nx, ny, nz, lgmax
      real*8 c10, u2, a, ab, aandb, expe, rho, dum
      real*8 ak, al
      real*8 pa, qa, ra, pb, qb, rb, d2, pc, qc, rc
      real*8 akxk, akyk, akzk
      real*8 axac, ayac, azac
      real*8 bxbc, bybc, bzbc
      real*8 bbrrk, b, b1, brrk
      real*8 c3x, c4x, c3y, c4y, c3z, c4z, csl, csk

      integer ijn1, ijn2
      real*8 pi252, dzero, pt5, done, factor
c
      data ijn1,ijn2 /125,25/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
c
c     ----- select expansion centre for -xyz- integrals -----
c
      if ((ncent.gt.2) .and. (lit.lt.ljt)) then
         ni = ljt - 1
         nj = lit - 1
         ij1 = ijn2
         ij2 = ijn1
         pc = pj
         qc = qj
         rc = rj
         dxij = pj - pi
         dyij = qj - qi
         dzij = rj - ri
      else
         ni = lit - 1
         nj = ljt - 1
         ij1 = ijn1
         ij2 = ijn2
         pc = pi
         qc = qi
         rc = ri
         dxij = pi - pj
         dyij = qi - qj
         dzij = ri - rj
      end if
      nmax = ni + nj
      max = nmax + 1
      do 20 i = 1 , max
         n = i - 1
         if (n.le.ni) in(i) = ij1*n + 1
         if (n.gt.ni) in(i) = ij1*ni + ij2*(n-ni) + 1
 20   continue
c
c     ----- k primitive
c
      lgmax = ngd
      do 200 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*pk
         akyk = ak*qk
         akzk = ak*rk
         csk = csc(kg)*factor
c
c     ----- l primitive
c
         if (okanl) lgmax = kg
         do 190 lg = 1 , lgmax
            al = dg(lg)
            b = ak + al
            b1 = done/b
            bbrrk = al*brrk*b1
            if (bbrrk.le.tol) then
               csl = csd(lg)
               pb = (akxk+al*pl)*b1
               qb = (akyk+al*ql)*b1
               rb = (akzk+al*rl)*b1
               bxbc = b*(pb-pc)
               bybc = b*(qb-qc)
               bzbc = b*(rb-rc)
c
c     ----- density factor
c
               d2 = csk*csl*b1
               if (okanl .and. (kg.gt.lg)) d2 = d2 + d2
c
c     ----- pair of i,j primitives
c
               nn = 0
               do 180 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.gt.tol) then
                     nn = nn + 16
                     go to 180
                  else
                     do 30 i = 1 , ij
                        dij(i) = dd(ijd(i)+nn)
 30                  continue
                     a = aa(n)
                     ab = a*b
                     aandb = a + b
                     expe = d2*dexp(-dum)/dsqrt(aandb)
                     rho = ab/aandb
                     pa = p1(n)
                     qa = q1(n)
                     ra = r1(n)
                     pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                     axac = a*(pa-pc)
                     ayac = a*(qa-qc)
                     azac = a*(ra-rc)
                     c3x = bxbc + axac
                     c4x = b*axac
                     c3y = bybc + ayac
                     c4y = b*ayac
                     c3z = bzbc + azac
                     c4z = b*azac
c
c     ----- roots and weights for quadrature
c
                     if (nroots.le.3) then
                        call rt123_dft
                     else if (nroots.eq.4) then
                        call roots4_dft
                     else
                        call roots5_dft
                     end if
                     mm = 0
c
c     compute two-electron  integrals for each root
c
                     do 90 m = 1 , nroots
                        u2 = u(m)*rho
                        f00 = expe*w(m)
                        dum = done/(ab+u2*aandb)
                        b10 = (b+u2)*pt5*dum
                        pc00 = (u2*c3x+c4x)*dum
                        qc00 = (u2*c3y+c4y)*dum
                        rc00 = (u2*c3z+c4z)*dum
c
c     ----- i(0,0) -----
c
                        i1 = in(1) + mm
                        pint(i1) = done
                        qint(i1) = done
                        rint(i1) = f00
                        if (nmax.eq.0) then
                           mm = mm + 625
                           go to 90
                        else
c
c     ----- i(1,0) -----
c
                           i2 = in(2) + mm
                           pint(i2) = pc00
                           qint(i2) = qc00
                           rint(i2) = rc00*f00
                           if (nmax.le.1) then
                              mm = mm + 625
                              go to 90
                           else
c
c     ----- i(iii,0) -----
c
                              c10 = dzero
                              i3 = i1
                              i4 = i2
                              do 40 iii = 2 , nmax
                                 c10 = c10 + b10
                                 i5 = in(iii+1) + mm
                                 pint(i5) = c10*pint(i3) + pc00*pint(i4)
                                 qint(i5) = c10*qint(i3) + qc00*qint(i4)
                                 rint(i5) = c10*rint(i3) + rc00*rint(i4)
                                 i3 = i4
                                 i4 = i5
 40                           continue
                              if (nj.eq.0) then
                                 mm = mm + 625
                                 go to 90
                              else
c
c     ----- i(iii,jjj,0,0) -----
c
                                 i5 = in(nmax+1) + mm
                                 min = ni
                              end if
                           end if
                        end if
 50                     iii = nmax
                        i3 = i5
 60                     i4 = in(iii) + mm
                        pint(i3) = pint(i3) + dxij*pint(i4)
                        qint(i3) = qint(i3) + dyij*qint(i4)
                        rint(i3) = rint(i3) + dzij*rint(i4)
                        i3 = i4
                        iii = iii - 1
                        if (iii.gt.min) go to 60
                        min = min + 1
                        if (min.lt.nmax) go to 50
                        if (ni.ne.0) then
                           i3 = ij2 + i1
                           do 80 jjj = 1 , nj
                              i4 = i3
                              do 70 iii = 1 , ni
                                 pint(i4) = pint(i4+ij1-ij2)
     +                              + dxij*pint(i4-ij2)
                                 qint(i4) = qint(i4+ij1-ij2)
     +                              + dyij*qint(i4-ij2)
                                 rint(i4) = rint(i4+ij1-ij2)
     +                              + dzij*rint(i4-ij2)
                                 i4 = i4 + ij1
 70                           continue
                              i3 = i3 + ij2
 80                        continue
                        end if
                        mm = mm + 625
 90                  continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                     go to (100,120,140,160,165) , nroots
                  end if
 100              do 110 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz))
     +                           *dij(iii) + gout(nnn)
 110              continue
                  nn = nn + 16
                  go to 180
 120              do 130 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625))
     +                       *dij(iii) + gout(nnn)
 130              continue
                  nn = nn + 16
                  go to 180
 140              do 150 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                      pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                      pint(nx+1250)*qint(ny+1250)*rint(nz+1250))
     +                      *dij(iii) + gout(nnn)
 150              continue
                  nn = nn + 16
                  go to 180
 160              do 170 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                       pint(nx+1250)*qint(ny+1250)*rint(nz+1250)+
     +                       pint(nx+1875)*qint(ny+1875)*rint(nz+1875))
     +                       *dij(iii) + gout(nnn)
 170              continue
                  nn = nn + 16
                  go to 180
 165              do 175 iii = 1 , ij
                     nx = ijx(iii)
                     ny = ijy(iii)
                     nz = ijz(iii)
                     nnn = ijgt(iii)
                     gout(nnn) = (pint(nx)*qint(ny)*rint(nz)+
     +                       pint(nx+625)*qint(ny+625)*rint(nz+625)+
     +                       pint(nx+1250)*qint(ny+1250)*rint(nz+1250)+
     +                       pint(nx+1875)*qint(ny+1875)*rint(nz+1875)+
     +                       pint(nx+2500)*qint(ny+2500)*rint(nz+2500))
     +                       *dij(iii) + gout(nnn)
 175              continue
                  nn = nn + 16
 180           continue
            end if
 190     continue
 200  continue
      return
      end
c
c
c
      subroutine indxa_dft(ijx,ijy,ijz,ij,mini,maxi,
     &     minj,maxj,iandj,inc1,inc2,inc3)
      implicit real*8  (a-h,o-z)
      logical iandj
      dimension ijx(225),ijy(225),ijz(225)
c
      integer ix, iy, iz
      common/dft_inxblk/ix(35),iy(35),iz(35)
c
      dimension jx(35),jy(35),jz(35)
c
      do 20 j = minj , maxj
         jx(j) = ix(j)*inc2
         jy(j) = iy(j)*inc2
         jz(j) = iz(j)*inc2
 20   continue
      ij = 0
      jmax = maxj
      do 40 i = mini , maxi

         nx = ix(i)*inc1 + inc3
         ny = iy(i)*inc1 + inc3
         nz = iz(i)*inc1 + inc3
         if (iandj) jmax = i
         do 30 j = minj , jmax
            ij = ij + 1
            ijx(ij) = nx + jx(j)
            ijy(ij) = ny + jy(j)
            ijz(ij) = nz + jz(j)
 30      continue
 40   continue
      return
      end

      subroutine shells_dft(gout,nelec,
     &     ish,jsh,ksh,lsh,
     &     basi, basj, bask, basl, ncentr,
     &     iexch,nintegral)

      implicit none 

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

      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
      real*8 cxyz, aaa
      integer ijaaa 
      common/junkx/cxyz(3,5625),aaa(21*mxp2),ijaaa(225)

      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c

c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      integer igt, jgt, kgt, lgt
      common /flipsx/ igt(3),jgt(3),kgt(3),lgt(3)
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c

c
c Arguments
c
      real*8 gout
      dimension gout(*)
      integer nelec
      integer iexch, ish, jsh, ksh, lsh
      integer basi, basj, bask, basl, ncentr
      integer nintegral

c
c Local variables
c
      integer i, i1, i2, j, j1, j2, k, k1, k2
      integer l, l1, l2 
      integer iat, jat, kat, lat
      integer nnx, nny, nnz
      integer ngtk, ngtl, max, jmax, lmax, ngij, ngti, ngtj

      integer ix, iy, iz, jx, jy, jz,kx ,ky, kz, lx, ly, lz

      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35),
     +          kx(35),ky(35),kz(35),lx(35),ly(35),lz(35)
      data lx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     +            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     +            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     +            2,  0,  2,  1,  1/
      data kx /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     +           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     +           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     +           10,  0, 10,  5,  5/
      data jx /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     +           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     +          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     +           50,  0, 50, 25, 25/
      data ix /   1,126,  1,  1,251,  1,  1,126,126,  1,
     +          376,  1,  1,251,251,126,  1,126,  1,126,
     +          501,  1,  1,376,376,126,  1,126,  1,251,
     +          251,  1,251,126,126/
      data ly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     +            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     +            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     +            0,  2,  1,  2,  1/
      data ky /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     +            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     +            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     +            0, 10,  5, 10,  5/
      data jy /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     +            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     +            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     +            0, 50, 25, 50, 25/
      data iy /   1,  1,126,  1,  1,251,  1,126,  1,126,
     +            1,376,  1,126,  1,251,251,  1,126,126,
     +            1,501,  1,126,  1,376,376,  1,126,251,
     +            1,251,126,251,126/
      data lz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     +            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     +            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     +            2,  2,  1,  1,  2/
      data kz /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     +            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     +            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     +           10, 10,  5,  5, 10/
      data jz /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     +            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     +            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     +           50, 50, 25, 25, 50/
      data iz /   1,  1,  1,126,  1,  1,251,  1,126,126,
     +            1,  1,376,  1,126,  1,126,251,251,126,
     +            1,  1,501,  1,126,  1,126,376,376,  1,
     +          251,251,126,126,251/

      if (nelec.eq.2) then

         okanl = ksh.eq.lsh .and. (ncentr .eq. 4)
         oident = ish.eq.ksh .and. jsh.eq.lsh .and. (ncentr .eq. 4)

         ngtk = kgt(iexch)
         ngtl = lgt(iexch)
c
c     ----- kshell
c
         kat = katom(bask,ksh)
         pk = c(1,kat)
         qk = c(2,kat)
         rk = c(3,kat)
         k1 = kstart(bask,ksh)
         k2 = k1 + kng(bask,ksh) - 1
         lkt = ktype(bask,ksh)
         mink = kmin(bask,ksh)
         maxk = kmax(bask,ksh)
         lock = kloc(bask,ksh) - mink
         ngc = 0
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex_m(bask,k)
            csc(ngc) = cs(bask,k)
            cpc(ngc) = cp(bask,k)
            cdc(ngc) = cd(bask,k)
            cfc(ngc) = cf(bask,k)
            cgc(ngc) = cg(bask,k)
 20      continue
c
c     ----- lshell
c

         if(basl .lt. 0) then
c
c  Dummy l centr
c
            lat = kat
            pl = c(1,lat)
            ql = c(2,lat)
            rl = c(3,lat)
            l1 = 1
            l2 = 1
            llt = 1
            minl = 1
            maxl = 1
            locl = 0
            ngd = 0
c            if(odbg)write(6,*)'dshell l loop',l,locl
            do l = l1 , l2
               ngd = ngd + 1
               dg(ngd)  = 0.0d0
               csd(ngd) = 1.0d0
               cpd(ngd) = 1.0d0
               cdd(ngd) = 1.0d0
               cfd(ngd) = 1.0d0
               cgd(ngd) = 1.0d0
c               if(odbg)write(6,*)ngc,dg(ngc),csd(ngc),cpd(ngc),cdd(ngc)
            enddo
            rrk = 0.0d0
         else

         lat = katom(basl,lsh)
         pl = c(1,lat)
         ql = c(2,lat)
         rl = c(3,lat)
         l1 = kstart(basl,lsh)
         l2 = l1 + kng(basl,lsh) - 1
         llt = ktype(basl,lsh)
         minl = kmin(basl,lsh)
         maxl = kmax(basl,lsh)
         locl = kloc(basl,lsh) - minl
         ngd = 0
         do l = l1 , l2
            ngd = ngd + 1
            dg(ngd) = ex_m(basl,l)
            csd(ngd) = cs(basl,l)
            cpd(ngd) = cp(basl,l)
            cdd(ngd) = cd(basl,l)
            cfd(ngd) = cf(basl,l)
            cgd(ngd) = cg(basl,l)
         enddo
         rrk = ((pk-pl)**2+(qk-ql)**2+(rk-rl)**2)
        endif
         nroots = (lit+ljt+lkt+llt-4)/2 + 1
c
c     ----- prepare indices for pairs of (k,l) functions
c
         kl = 0
         lmax = maxl
         do 50 k = mink , maxk
            nnx = kx(k)
            nny = ky(k)
            nnz = kz(k)
            if (okanl) lmax = k 
            do 40 l = minl , lmax
               kl = kl + 1
               klx(kl) = nnx + lx(l)
               kly(kl) = nny + ly(l)
               klz(kl) = nnz + lz(l)
               klgt(kl) = ngtk*(k-mink) + ngtl*(l-minl)
 40         continue
 50      continue
         max = kl
         do 60 i = 1 , ij
            if (oident) max = i
            ik(i) = max
 60      continue
         ijkl = ij*kl
         if (oident) ijkl = ij*(ij+1)/2
c
c     zero integral storage
c
         do 65 i = 1,ij
         ngij = ijgt(i)
         nintegral = nintegral + ik(i)
         do 65 k = 1,ik(i)
  65     gout(ngij + klgt(k)) = 0.0d0
         return
      else
c
c@@ note we could triangulate for 3c case
c   but the logic adopted assumes we dont
c
         oianj = ish.eq.jsh .and. ncentr .eq. 4
         ngti = igt(iexch)
         ngtj = jgt(iexch)
c
c     ----- ishell
c
         iat = katom(basi,ish)
         pi = c(1,iat)
         qi = c(2,iat)
         ri = c(3,iat)
         i1 = kstart(basi,ish)
         i2 = i1 + kng(basi,ish) - 1
         lit = ktype(basi,ish)
         mini = kmin(basi,ish)
         maxi = kmax(basi,ish)
         loci = kloc(basi,ish) - mini
         nga = 0
         do 70 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex_m(basi,i)
            csa(nga) = cs(basi,i)
            cpa(nga) = cp(basi,i)
            cda(nga) = cd(basi,i)
            cfa(nga) = cf(basi,i)
            cga(nga) = cg(basi,i)
 70      continue

         if(basj .lt. 0)then
c
c Dummy j centre
c
            jat = iat
            pj = c(1,jat)
            qj = c(2,jat)
            rj = c(3,jat)
            j1 = 1
            j2 = 1
            ljt = 1
            minj = 1
            maxj = 1
            locj = 0
            ngb = 0

c            if(odbg)write(6,*)'dshell j loop',j, locj

            do j = j1 , j2
               ngb = ngb + 1
               bg(ngb)  = 0.0d0
               csb(ngb) = 1.0d0
               cpb(ngb) = 1.0d0
               cdb(ngb) = 1.0d0
               cfb(ngb) = 1.0d0
               cgb(ngb) = 1.0d0

c               if(odbg)write(6,*)ngb,bg(ngb),csb(ngb),cpb(ngb),cdb(ngb)

            enddo
            rri = 0.0d0
         else
c
c     ----- jshell
c
         jat = katom(basj,jsh)
         pj = c(1,jat)
         qj = c(2,jat)
         rj = c(3,jat)
         j1 = kstart(basj,jsh)
         j2 = j1 + kng(basj,jsh) - 1
         ljt = ktype(basj,jsh)
         minj = kmin(basj,jsh)
         maxj = kmax(basj,jsh)
         locj = kloc(basj,jsh) - minj
         ngb = 0
         do 80 j = j1 , j2
            ngb = ngb + 1
            bg(ngb) = ex_m(basj,j)
            csb(ngb) = cs(basj,j)
            cpb(ngb) = cp(basj,j)
            cdb(ngb) = cd(basj,j)
            cfb(ngb) = cf(basj,j)
            cgb(ngb) = cg(basj,j)
 80      continue
         rri = ((pi-pj)**2+(qi-qj)**2+(ri-rj)**2)
         endif
c
c     ----- prepare indices for pairs of (i,j) functions
c
         ij = 0
         jmax = maxj
         do 100 i = mini , maxi
            nnx = ix(i)
            nny = iy(i)
            nnz = iz(i)
            if (oianj) jmax = i
            do 90 j = minj , jmax
               ij = ij + 1
               ijx(ij) = nnx + jx(j)
               ijy(ij) = nny + jy(j)
               ijz(ij) = nnz + jz(j)
               ijgt(ij) = ngti*(i-mini) + ngtj*(j-minj) + 1 
 90         continue
 100     continue
         return
      end if
      end

      subroutine spchck_dft
c
c     ----- check for pure sp basis -----
c
      implicit none

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
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
c for documentation of these functions see the start
c of dft/basis.m
c
      integer BL_create_atomtag
      external BL_create_atomtag

      integer BL_find_atomtag
      external BL_find_atomtag

      integer BL_import_shell
      external BL_import_shell

      integer BL_assign_types_by_z
      external BL_assign_types_by_z

      integer BL_assign_type
      external BL_assign_type

      integer BL_write_basis
      external BL_write_basis

      integer BL_clear_basis_set
      external BL_clear_basis_set

      integer BL_maxang_on_atom
      external BL_maxang_on_atom

      integer BL_basis_size
      external BL_basis_size

      integer BL_max_shell_count
      external BL_max_shell_count

      integer BL_num_sets
      external BL_num_sets

      integer BL_num_types
      external BL_num_types

      integer BL_get_atom_type
      external BL_get_atom_type

      logical BL_atomtyp_exist
      external BL_atomtyp_exist

      integer BL_summarise
      external BL_summarise
c
c
      logical ospbas, onocnt, opdbas, opfbas, opgbas
      integer kad
      common /ijlabx/ kad(4,mxshel),ospbas,onocnt,opdbas,opfbas,opgbas
c
c

      integer lbas, i, ii

      ospbas = .true.
      opdbas = .false.
      opfbas = .false.
      opgbas = .false.
      do lbas=1,BL_num_sets()
        do 20 i = 1 , nshell(lbas)
           kad(lbas,i) = -1
c          ii = i
           if (ktype(lbas,i).gt.2) then
             ospbas = .false.
             if (ktype(lbas,i).eq.3) opdbas = .true.
             if (ktype(lbas,i).eq.4) opfbas = .true.
             if (ktype(lbas,i).eq.5) opgbas = .true.
c            kad(lbas,i) = -1
           end if
 20     continue
c       if (intg76.eq.0) then
c         do 30 i = 1 , nshell(lbas)
c            kad(lbas,i) = -1
c30       continue
c       end if
      enddo
      return
      end
      subroutine spdint_dft(gout)
c
c     ----- form integrals over functions -----
c
      implicit none
c
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
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
      real*8 pin, qin, rin, aaa
      integer ijaaa
      common/junkx/pin(5625),qin(5625),rin(5625),aaa(21*mxp2),
     &     ijaaa(225)

c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 dij, dkl
      common /denssx /dij(225),dkl(225)
c

      real*8 gout
      dimension gout(*)
c
c Local
c
      integer i, k, n, n1, nx, ny, nz, mx, my, mz, max
      real*8 d1

      go to (20,50,80,110,140,170,200,230,260) , nroots
 20   do 40 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 30 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz))*d1*dkl(k) + gout(n)
 30      continue
 40   continue
      return
 50   do 70 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 60 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+pin(mx+625)*qin(my+625)
     +                *rin(mz+625))*d1*dkl(k) + gout(n)
 60      continue
 70   continue
      return
 80   do 100 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 90 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250))
     +                *d1*dkl(k) + gout(n)
 90      continue
 100  continue
      return
 110  do 130 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 120 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875))
     +                *d1*dkl(k) + gout(n)
 120     continue
 130  continue
      return
 140  do 160 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 150 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500))
     +                *d1*dkl(k) + gout(n)
 150     continue
 160  continue
      return
 170  do 190 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 180 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125))
     +                *d1*dkl(k) + gout(n)
 180     continue
 190  continue
      return
 200  do 220 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 210 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750))
     +                *d1*dkl(k) + gout(n)
 210     continue
 220  continue
      return
 230  do 250 i = 1 , ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 240 k = 1 , max
            mx = nx + klx(k)
            my = ny + kly(k)
            mz = nz + klz(k)
            n = n1 + klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750)+
     +                 pin(mx+4375)*qin(my+4375)*rin(mz+4375))
     +                *d1*dkl(k)+gout(n)
 240     continue
 250  continue
      return
c
 260  do 280 i = 1,ij
         d1 = dij(i)
         nx = ijx(i)
         ny = ijy(i)
         nz = ijz(i)
         n1 = ijgt(i)
         max = ik(i)
         do 270 k = 1,max
            mx = nx+klx(k)
            my = ny+kly(k)
            mz = nz+klz(k)
            n = n1+klgt(k)
            gout(n) = (pin(mx)*qin(my)*rin(mz)+
     +                 pin(mx+625)*qin(my+625)*rin(mz+625)+
     +                 pin(mx+1250)*qin(my+1250)*rin(mz+1250)+
     +                 pin(mx+1875)*qin(my+1875)*rin(mz+1875)+
     +                 pin(mx+2500)*qin(my+2500)*rin(mz+2500)+
     +                 pin(mx+3125)*qin(my+3125)*rin(mz+3125)+
     +                 pin(mx+3750)*qin(my+3750)*rin(mz+3750)+
     +                 pin(mx+4375)*qin(my+4375)*rin(mz+4375)+
     +                 pin(mx+5000)*qin(my+5000)*rin(mz+5000))
     +                *d1*dkl(k)+gout(n)
 270     continue
 280  continue
      return
      end
      subroutine ssdss_dft(qq)
      implicit none

      real*8 qq(*)

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
      real*8 ag,csa,cpa,cda,cfa,cga,bg,csb,cpb,cdb,cfb,cgb
      real*8 cgg,csc,cpc,cdc,cfc,cgc,dg,csd,cpd,cdd,cfd,cgd
      real*8 xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk
      real*8 exij,rsmall
      integer nga,ngb,ngc,ngd
      common/dshlnf/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +              cfa(mxprms),cga(mxprms),
     +              bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +              cfb(mxprms),cgb(mxprms),
     +             cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +              cfc(mxprms),cgc(mxprms),
     +              dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +              cfd(mxprms),cgd(mxprms),
     +              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     +              nga,ngb,ngc,ngd,exij(mxprms*mxprms),rsmall
c
      integer mxp2
      parameter (mxp2 = mxprms*mxprms)

      integer ncmax
      parameter(ncmax=65)

      real*8 ddij,ddkl,aei,aej,aek,ael,aa,r,x1,y1,z1,dd
      integer ijden,ik,ijx,ijy,ijz,klx,kly,klz
      real*8 dij,dkl
      integer ijgt,klgt

      common/dft_d2escr/ddij(ncmax,225),ddkl(ncmax,225),
     +  aei(ncmax),aej(ncmax),aek(ncmax),ael(ncmax),
     +  aa(mxp2),r(mxp2),x1(mxp2),y1(mxp2),z1(mxp2),
     +  dd(4*mxp2),ijden(225),ik(225),
     +  ijx(225),ijy(225),ijz(225),klx(225),kly(225),klz(225),
     +  dij(225),dkl(225),ijgt(225),klgt(225)
      integer kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      common/dft_incrd/kln2,kln1,ijn2,ijn1,inc1,inc2,inc3,inc4,inc5,
     + lendd,ioff,ixi,iyi,izi,iabd,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ida,
     + ifok, itop
      logical outd,oskip,oiandj,okandl,osame,onorm,oclos,ogrhf
      logical onocor,onopen
      real*8 tol1,tol2,tol3,tol4,dcut
      integer natomd,npass, norb
      common/dft_dmisc/tol1,tol2,tol3,tol4,dcut,natomd(4),npass,outd,
     + oskip(4),oiandj,okandl,osame,onorm,oclos,ogrhf,
     + norb,onocor,onopen

      real*8 f1, rt1, g, ww1, f00, y, u2, ai, aj, ak, al
      real*8 xa, ya, za, xb, yb, zb, xd, yd, zd
      real*8 rho, xinv, x, a, b, b1
      real*8 c1x, c1y, c1z, c2x, c2y, c2z
      real*8 axad, ayad, azad
      real*8 akxk, akyk, akzk
      real*8 bxbd, bybd, bzbd
      real*8 brrk, bbrrk, exkl, expe, csk
      real*8 aandb, ab, alpha
      integer joff, n, nn, ig, jg, kg, lg, iper, i, jgmax

      real*8 pt5, pt2, third, one, pi252, pie4
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
c
c  original code - for ncentr == 2, the apply k/l switch
c  has been deleted (only the else clause remains)
c

      subroutine sskl_dft(gout,ncentr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
      common/junkx/pint(5625),qint(5625),rint(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indezx/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +                ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setintx/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                 bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                 f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                 ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denssx /dij(225),dkl(225)
c
      dimension gout(*)
      data kln1,kln2 /5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d0/
      data sqrt7 /2.64575131106459d0/
      data pi252 /34.986836655250d0/
      data dzero,pt5,done /0.0d0,0.5d0,1.0d0/
      factor = pi252*qq4
      onorm = normf.ne.1 .or. normp.ne.1
      if(ncentr .ne. 4) onorm = .false.
      if ((lkt.lt.llt)) then

c shouldn't get here
         if(ncentr .ne. 4)call caserr('sskl')

         nk = llt - 1
         nl = lkt - 1
         kl1 = kln2
         kl2 = kln1
         pd = pl
         qd = ql
         rd = rl
         dxkl = pl - pk
         dykl = ql - qk
         dzkl = rl - rk
      else
c
c     ----- select expansion centre for -xyz- integrals -----
c
         nk = lkt - 1
         nl = llt - 1
         kl1 = kln1
         kl2 = kln2
         pd = pk
         qd = qk
         rd = rk
         dxkl = pk - pl
         dykl = qk - ql
         dzkl = rk - rl
      end if
      mmax = nk + nl
      max = mmax + 1
      do 20 k = 1 , max
         n = k - 1
         if (n.le.nk) kn(k) = kl1*n
         if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 20   continue
      lgmax = ngd
c
c     ----- k primitive
c
      do 370 kg = 1 , ngc
         ak = cgg(kg)
         brrk = ak*rrk
         akxk = ak*pk
         akyk = ak*qk
         akzk = ak*rk
         csk = csc(kg)*factor
         cpk = cpc(kg)*factor
         cdk = cdc(kg)*factor
         cfk = cfc(kg)*factor
         cgk = cgc(kg)*factor
         if (okanl) lgmax = kg
c
c     ----- l primitive
c
         do 360 lg = 1 , lgmax
            al = dg(lg)
            b = ak + al
            b1 = done/b
            bbrrk = al*brrk*b1
            if (bbrrk.le.tol) then
               csl = csd(lg)
               cpl = cpd(lg)
               cdl = cdd(lg)
               cfl = cfd(lg)
               cgl = cgd(lg)
               pb = (akxk+al*pl)*b1
               qb = (akyk+al*ql)*b1
               rb = (akzk+al*rl)*b1
               bxbd = b*(pb-pd)
               bybd = b*(qb-qd)
               bzbd = b*(rb-rd)
               odoub = okanl .and. kg.gt.lg
c
c     ----- density factor
c
               n = 0
               max = maxl
               do 200 k = mink , maxk
                  go to (30,40,100,100,
     +                   50,100,100,60,100,100,
     +                   70,100,100,80,100,100,100,100,100,90,
     +                   92,100,100,94,100,100,100,100,100,96,
     +                  100,100,98,100,100), k
c
 30               dum1 = csk*b1
                  go to 100
 40               dum1 = cpk*b1
                  go to 100
 50               dum1 = cdk*b1
                  go to 100
 60               if (onorm) dum1 = dum1*sqrt3
                  go to 100
 70               dum1 = cfk*b1
                  go to 100
 80               if (onorm) dum1 = dum1*sqrt5
                  go to 100
 90               if (onorm) dum1 = dum1*sqrt3
                  go to 100
 92               dum1 = cgk*b1
                  go to 100
 94               if (onorm) dum1 = dum1*sqrt7
                  go to 100
 96               if (onorm) dum1 = dum1*sqrt5/sqrt3
                  go to 100
 98               if (onorm) dum1 = dum1*sqrt3
 100              if (okanl) max = k
                  do 190 l = minl , max
                     go to (110,120,180,180,
     +                      130,180,180,140,180,180,
     +                      150,180,180,160,180,180,180,180,180,170,
     +                      182,180,180,184,180,180,180,180,180,186,
     +                      180,180,188,180,180), l
 110                 dum2 = dum1*csl
                     if (odoub) then
                        if (k.gt.1) then
                           dum2 = dum2 + csk*cpl*b1
                        else
                           dum2 = dum2 + dum2
                        end if
                     end if
                     go to 180
 120                 dum2 = dum1*cpl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 130                 dum2 = dum1*cdl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 140                 if (onorm) dum2 = dum2*sqrt3
                     go to 180
 150                 dum2 = dum1*cfl
                     if (odoub) dum2 = dum2 + dum2
                     go to 180
 160                 if (onorm) dum2 = dum2*sqrt5
                     go to 180
 170                 if (onorm) dum2 = dum2*sqrt3
                     go to 180
 182                 dum2 = dum1*cgl
                     if (odoub) dum2 = dum2+dum2
                     go to 180
 184                 if (onorm) dum2 = dum2 * sqrt7
                     go to 180
 186                 if (onorm) dum2 = dum2 * sqrt5/sqrt3
                     go to 180
 188                 if (onorm) dum2 = dum2 * sqrt3
 180                 n = n + 1
                     dkl(n) = dum2
 190              continue
 200           continue
               nn = 0
c
c     ----- pair of i,j primitives
c
               do 350 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.gt.tol) then
                     nn = nn + 16
                     go to 350
                  else
                     a = aa(n)
                     ab = a*b
                     aandb = a + b
                     expe = dd(1+nn)*dexp(-dum)/dsqrt(aandb)
                     rho = ab/aandb
                     pa = p1(n)
                     qa = q1(n)
                     ra = r1(n)
                     pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                     axad = a*(pa-pd)
                     ayad = a*(qa-qd)
                     azad = a*(ra-rd)
                     c1x = bxbd + axad
                     c2x = a*bxbd
                     c1y = bybd + ayad
                     c2y = a*bybd
                     c1z = bzbd + azad
                     c2z = a*bzbd
                     if (nroots.le.3) then
                        call rt123_dft
                     else if (nroots.eq.4) then
                        call roots4_dft
                     else
                        call roots5_dft
                     end if
c
c     ----- roots and weights for quadrature
c
                     mm = 0
                     do 260 m = 1 , nroots
c
c     compute two-electron  integrals for each root
c
                        u2 = u(m)*rho
                        f00 = expe*w(m)
                        dum = done/(ab+u2*aandb)
                        bp01 = (a+u2)*pt5*dum
                        pcp00 = (u2*c1x+c2x)*dum
                        qcp00 = (u2*c1y+c2y)*dum
                        rcp00 = (u2*c1z+c2z)*dum
c
c     ----- i(0,0) -----
c
                        i1 = 1 + mm
                        pint(i1) = done
                        qint(i1) = done
                        rint(i1) = f00
                        if (mmax.eq.0) then
                           mm = mm + 625
                           go to 260
                        else
c
c     ----- i(0,1) -----
c
                           k2 = kn(2)
                           i3 = i1 + k2
                           pint(i3) = pcp00
                           qint(i3) = qcp00
                           rint(i3) = rcp00*f00
                           if (mmax.le.1) then
                              mm = mm + 625
                              go to 260
                           else
c
c     ----- i(0,kkk) -----
c
                              cp01 = dzero
                              i3 = i1
                              i4 = i1 + k2
                              do 210 kkk = 2 , mmax
                                 cp01 = cp01 + bp01
                                 i5 = i1 + kn(kkk+1)
                                 pint(i5) = cp01*pint(i3)
     +                              + pcp00*pint(i4)
                                 qint(i5) = cp01*qint(i3)
     +                              + qcp00*qint(i4)
                                 rint(i5) = cp01*rint(i3)
     +                              + rcp00*rint(i4)
                                 i3 = i4
                                 i4 = i5
 210                          continue
                              if (nl.eq.0) then
                                 mm = mm + 625
                                 go to 260
                              else
c
c     ----- i(0,0,kkk,lll) -----
c
                                 i5 = kn(mmax+1)
                                 min = nk
                              end if
                           end if
                        end if
 220                    kkk = mmax
                        i3 = i1 + i5
 230                    i4 = i1 + kn(kkk)
                        pint(i3) = pint(i3) + dxkl*pint(i4)
                        qint(i3) = qint(i3) + dykl*qint(i4)
                        rint(i3) = rint(i3) + dzkl*rint(i4)
                        i3 = i4
                        kkk = kkk - 1
                        if (kkk.gt.min) go to 230
                        min = min + 1
                        if (min.lt.mmax) go to 220
                        if (nk.ne.0) then
                           i3 = i1 + kl2
                           do 250 lll = 1 , nl
                              i4 = i3
                              do 240 kkk = 1 , nk
                                 pint(i4) = pint(i4+kl1-kl2)
     +                              + dxkl*pint(i4-kl2)
                                 qint(i4) = qint(i4+kl1-kl2)
     +                              + dykl*qint(i4-kl2)
                                 rint(i4) = rint(i4+kl1-kl2)
     +                              + dzkl*rint(i4-kl2)
                                 i4 = i4 + kl1
 240                          continue
                              i3 = i3 + kl2
 250                       continue
                        end if
                        mm = mm + 625
 260                 continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                     go to (270,290,310,330,335) , nroots
                  end if
 270              do 280 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz))
     +                         *dkl(kkk) + gout(nnn)
 280              continue
                  nn = nn + 16
                  go to 350
 290              do 300 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625))
     +                       *dkl(kkk) + gout(nnn)
 300              continue
                  nn = nn + 16
                  go to 350
 310              do 320 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250))
     +                       *dkl(kkk) + gout(nnn)
 320              continue
                  nn = nn + 16
                  go to 350
 330              do 340 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250)+
     +                       pint(mx+1875)*qint(my+1875)*rint(mz+1875))
     +                         *dkl(kkk) + gout(nnn)
 340              continue
                  nn = nn + 16
                  go to 350
 335              do 345 kkk = 1 , kl
                     mx = 1 + klx(kkk)
                     my = 1 + kly(kkk)
                     mz = 1 + klz(kkk)
                     nnn = 1 + klgt(kkk)
                     gout(nnn) = (pint(mx)*qint(my)*rint(mz)+
     +                       pint(mx+625)*qint(my+625)*rint(mz+625)+
     +                       pint(mx+1250)*qint(my+1250)*rint(mz+1250)+
     +                       pint(mx+1875)*qint(my+1875)*rint(mz+1875)+
     +                       pint(mx+2500)*qint(my+2500)*rint(mz+2500))
     +                         *dkl(kkk) + gout(nnn)
 345              continue
                  nn = nn + 16
 350           continue
            end if
 360     continue
 370  continue
      return
      end
      subroutine ssprim_dft

      implicit none

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
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c

      real*8 cxyz,a,r,p1,q1,r1,dij
      integer ijd
      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
      common/junkx/cxyz(3,5625),
     +   a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     +   ijd(225)
      real*8 done
      real*8 ai, arri, axi, ayi, azi
      real*8 aj, aa, aa1, dum, csi
      integer nn, ia, jb, jbmax
c
      data done /1.0d0/

c
c     ----- i primitive
c
      nij = 0
      jbmax = ngb

cc      write(6,*)'ijprim',nga, ngb

      do 30 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*pi
         ayi = ai*qi
         azi = ai*ri
         csi = csa(ia)
         if (oianj) jbmax = ia
c
c     ----- j primitive
c
         do 20 jb = 1 , jbmax
            aj = bg(jb)
            aa = ai + aj
            aa1 = done/aa
            dum = aj*arri*aa1

cc            write(6,*)'ijprim ag, bg', ag(ia), bg(jb), dum, tol

            if (dum.le.tol) then
               nn = 1 + nij*16
               nij = nij + 1
               r(nij) = dum
               a(nij) = aa
               p1(nij) = (axi+aj*pj)*aa1
               q1(nij) = (ayi+aj*qj)*aa1
               r1(nij) = (azi+aj*rj)*aa1
               dij(nn) = csi*csb(jb)*aa1
               if (oianj .and. (ia.gt.jb)) dij(nn) = dij(nn) + dij(nn)
            end if
 20      continue
 30   continue
      return
      end
c
c ncent - class of integral
c
c onorm - only applies to 4 centre case
c
      subroutine genral_dft(gout,ncentr)
c
      implicit none
c
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

ccccINCLUDE(common/cslosc)

c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c

      integer mxp2
      parameter (mxp2 = mxprms * mxprms)

      real*8 pin, qin, rin,aa,r,p1,q1,r1,dd
      integer ijd
      common/junkx/pin(5625), qin(5625), rin(5625),
     + aa(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dd(16*mxp2),
     + ijd(225)

c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
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
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setintx/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                 bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                 f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                 ij1,ij2,kl1,kl2
c
c
      real*8 dij, dkl
      common /denssx /dij(225),dkl(225)
c

      real*8 gout
      dimension gout(*)
      integer ncentr
c
      real*8 axac, ayac, azac
      real*8 axad, ayad, azad, bxbd, bybd, bzbd
      real*8 bxbc, bybc, bzbc
      real*8 csk, cpk, cdk, cfk, cgk
      real*8 csl, cpl, cdl, cfl, cgl
      real*8 c1x, c2x, c3x, c4x
      real*8 c1y, c2y, c3y, c4y
      real*8 c1z, c2z, c3z, c4z
      real*8 brrk, bbrrk, akxk, akyk, akzk
      real*8 dum, dum1, dum2, a, b, ab, aandb
      real*8 ai, aj, ak, al
      real*8 pb, qb, rb
      real*8 pa, qa, ra
      real*8 pc, qc, rc
      real*8 pd, qd, rd
      real*8 a1, b1, rho, expe, factor, u2
      logical odoub, onorm
      integer i, l, k, kg, lg, lgmax, mm, nn, m, n, max
      integer in1(9)

      integer ijn1,ijn2,kln1,kln2

      real*8 sqrt3,  sqrt5, sqrt7, done, pi252, pt5

      data ijn1,ijn2,kln1,kln2 /125,25,5,1/
      data sqrt3 /1.73205080756888d0/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data pi252 /34.986836655250d0/
      data pt5,done /0.5d0,1.0d0/
c
      if (ijkl.eq.1) then
c
c     ----- (s,s//s,s) -----
c
         call s0000_dft(gout)
      else if (ij.eq.1) then
c
c     ----- (s,s//k,l) -----
c
         call sskl_dft(gout,ncentr)
      else if (kl.eq.1) then
         call ijss_dft(gout,ncentr)
      else
         factor = pi252*qq4
         onorm = normf.ne.1 .or. normp.ne.1
         if(ncentr .lt. 4) onorm = .false.
c
c     ----- select expansion centre for -xyz- integrals -----
c
         if ((ncentr .gt. 2) .and. (lit.lt.ljt)) then
            ni = ljt - 1
            nj = lit - 1
            ij1 = ijn2
            ij2 = ijn1
            pc = pj
            qc = qj
            rc = rj
            dxij = pj - pi
            dyij = qj - qi
            dzij = rj - ri
         else
            ni = lit - 1
            nj = ljt - 1
            ij1 = ijn1
            ij2 = ijn2
            pc = pi
            qc = qi
            rc = ri
            dxij = pi - pj
            dyij = qi - qj
            dzij = ri - rj
         end if
         if (lkt.lt.llt) then

            if(ncentr .ne. 4)call caserr('genral')

            nk = llt - 1
            nl = lkt - 1
            kl1 = kln2
            kl2 = kln1
            pd = pl
            qd = ql
            rd = rl
            dxkl = pl - pk
            dykl = ql - qk
            dzkl = rl - rk
         else
            nk = lkt - 1
            nl = llt - 1
            kl1 = kln1
            kl2 = kln2
            pd = pk
            qd = qk
            rd = rk
            dxkl = pk - pl
            dykl = qk - ql
            dzkl = rk - rl
         end if
         nmax = ni + nj
         mmax = nk + nl
         max = nmax + 1
         do 20 i = 1 , max
            n = i - 1
            if (n.le.ni) in1(i) = ij1*n + 1
            if (n.gt.ni) in1(i) = ij1*ni + ij2*(n-ni) + 1
 20      continue
         max = mmax + 1
         do 30 k = 1 , max
            n = k - 1
            if (n.le.nk) kn(k) = kl1*n
            if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 30      continue
c
c     ----- k primitive
c
         lgmax = ngd
         do 270 kg = 1 , ngc
            ak = cgg(kg)
            brrk = ak*rrk
            akxk = ak*pk
            akyk = ak*qk
            akzk = ak*rk
            csk = csc(kg)*factor
            cpk = cpc(kg)*factor
            cdk = cdc(kg)*factor
            cfk = cfc(kg)*factor
            cgk = cgc(kg)*factor
c
c     ----- l primitive
c
            if (okanl) lgmax = kg
            do 260 lg = 1 , lgmax
               al = dg(lg)
               b = ak + al
               b1 = done/b
               bbrrk = al*brrk*b1
               if (bbrrk.le.tol) then
                  csl = csd(lg)
                  cpl = cpd(lg)
                  cdl = cdd(lg)
                  cfl = cfd(lg)
                  cgl = cgd(lg)
                  pb = (akxk+al*pl)*b1
                  qb = (akyk+al*ql)*b1
                  rb = (akzk+al*rl)*b1
                  bxbd = b*(pb-pd)
                  bybd = b*(qb-qd)
                  bzbd = b*(rb-rd)
                  bxbc = b*(pb-pc)
                  bybc = b*(qb-qc)
                  bzbc = b*(rb-rc)
c
c     ----- density factor
c
                  odoub = okanl .and. kg.gt.lg
                  n = 0
                  max = maxl
                  do 210 k = mink , maxk
                     go to (40,50,110,110,
     +                 60,110,110,70,110,110,
     +                 80,110,110,90,110,110,110,110,110,100,
     +                 92,110,110,94,110,110,110,110,110,96,
     +                110,110,98,110,110), k
 40                  dum1 = csk*b1
                     go to 110
 50                  dum1 = cpk*b1
                     go to 110
 60                  dum1 = cdk*b1
                     go to 110
 70                  if (onorm) dum1 = dum1*sqrt3
                     go to 110
 80                  dum1 = cfk*b1
                     go to 110
 90                  if (onorm) dum1 = dum1*sqrt5
                     go to 110
 100                 if (onorm) dum1 = dum1*sqrt3
                     go to 110
 92                  dum1 = cgk*b1
                     go to 110
 94                  if (onorm) dum1 = dum1*sqrt7
                     go to 110
 96                  if (onorm) dum1 = dum1*sqrt5/sqrt3
                     go to 110
 98                  if (onorm) dum1 = dum1*sqrt3
 110                 if (okanl) max = k
                     do 200 l = minl , max
                        go to (120,130,190,190,
     +                    140,190,190,150,190,190,
     +                    160,190,190,170,190,190,190,190,190,180,
     +                    182,190,190,184,190,190,190,190,190,186,
     +                    190,190,188,190,190), l
c
 120                    dum2 = dum1*csl
                        if (odoub) then
                           if (k.gt.1) then
                              dum2 = dum2 + csk*cpl*b1
                           else
                              dum2 = dum2 + dum2
                           end if
                        end if
                        go to 190
 130                    dum2 = dum1*cpl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 140                    dum2 = dum1*cdl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 150                    if (onorm) dum2 = dum2 * sqrt3
                        go to 190
 160                    dum2 = dum1*cfl
                        if (odoub) dum2 = dum2 + dum2
                        go to 190
 170                    if (onorm) dum2 = dum2 * sqrt5
                        go to 190
 180                    if (onorm) dum2 = dum2 * sqrt3
                        go to 190
 182                    dum2 = dum1*cgl
                        if (odoub) dum2 = dum2+dum2
                        go to 190
 184                   if (onorm) dum2 = dum2 * sqrt7
                        go to 190
 186                   if (onorm) dum2 = dum2 * sqrt5/sqrt3
                        go to 190
 188                   if (onorm) dum2 = dum2 * sqrt3
 190                    n = n + 1
                        dkl(n) = dum2
 200                 continue
 210              continue
c
c     ----- pair of i,j primitives
c
                  nn = 0
                  do 250 n = 1 , nij
                     dum = bbrrk + r(n)
                     if (dum.le.tol) then
                        do 220 i = 1 , ij
                           dij(i) = dd(ijd(i)+nn)
 220                    continue
                        a = aa(n)
                        ab = a*b
                        aandb = a + b
                        expe = dexp(-dum)/dsqrt(aandb)
                        rho = ab/aandb
                        pa = p1(n)
                        qa = q1(n)
                        ra = r1(n)
                        pp = rho*((pa-pb)**2+(qa-qb)**2+(ra-rb)**2)
                        axad = a*(pa-pd)
                        ayad = a*(qa-qd)
                        azad = a*(ra-rd)
                        axac = a*(pa-pc)
                        ayac = a*(qa-qc)
                        azac = a*(ra-rc)
                        c1x = bxbd + axad
                        c2x = a*bxbd
                        c3x = bxbc + axac
                        c4x = b*axac
                        c1y = bybd + ayad
                        c2y = a*bybd
                        c3y = bybc + ayac
                        c4y = b*ayac
                        c1z = bzbd + azad
                        c2z = a*bzbd
                        c3z = bzbc + azac
                        c4z = b*azac
c
c     ----- roots and weights for quadrature
c
                        if (nroots.le.3) call rt123_dft
                        if (nroots.eq.4) call roots4_dft
                        if (nroots.eq.5) call roots5_dft
                        if (nroots.ge.6) call rootss_dft
                        mm = 0
                        max = nmax + 1
                        do 240 m = 1 , nroots
                           u2 = u(m)*rho
                           f00 = expe*w(m)
                           do 230 i = 1 , max
                              in(i) = in1(i) + mm
 230                       continue
                           dum = done/(ab+u2*aandb)
                           dum2 = pt5*dum
                           bp01 = (a+u2)*dum2
                           b00 = u2*dum2
                           b10 = (b+u2)*dum2
                           pcp00 = (u2*c1x+c2x)*dum
                           pc00 = (u2*c3x+c4x)*dum
                           qcp00 = (u2*c1y+c2y)*dum
                           qc00 = (u2*c3y+c4y)*dum
                           rcp00 = (u2*c1z+c2z)*dum
                           rc00 = (u2*c3z+c4z)*dum
                           call xyzint_dft
                           mm = mm + 625
 240                    continue
c
c     ----- form (i,j//k,l) integrals over functions
c
                        call spdint_dft(gout)
                     end if
                     nn = nn + 16
 250              continue
               end if
 260        continue
 270     continue
      end if
      return
      end
      subroutine s0000_dft(g)

      implicit none

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
      integer mxp2
      parameter (mxp2 = mxprms * mxprms)
c
      real*8 ag, csa, cpa, cda, cfa, cga
      real*8 bg, csb, cpb, cdb, cfb, cgb
      real*8 cgg, csc, cpc, cdc, cfc, cgc
      real*8 dg, csd, cpd ,cdd, cfd, cgd
      real*8 pi, qi, ri, pj, qj, rj, rri, pk, qk, rk, pl, ql, rl, rrk
      integer nga, ngb, ngc, ngd
      common /shlinfx/ag(mxprms),csa(mxprms),cpa(mxprms),cda(mxprms),
     +                cfa(mxprms),cga(mxprms),
     +                bg(mxprms),csb(mxprms),cpb(mxprms),cdb(mxprms),
     +                cfb(mxprms),cgb(mxprms),
     +                cgg(mxprms),csc(mxprms),cpc(mxprms),cdc(mxprms),
     +                cfc(mxprms),cgc(mxprms),
     +                dg(mxprms),csd(mxprms),cpd(mxprms),cdd(mxprms),
     +                cfd(mxprms),cgd(mxprms),
     +                pi,qi,ri,pj,qj,rj,rri,pk,qk,rk,pl,ql,rl,rrk,
     +                nga,ngb,ngc,ngd
c
c
      real*8 tol, cutoff
      integer icount, ic4, isti, jsti, ksti, lsti, lastb, lastu
      integer len4, lennx
      logical out, outv
      common /shltx/ tol,cutoff,icount,ic4,out,
     +               isti, jsti, ksti, lsti, lastb, lastu, outv,
     +               len4, lennx
c

      real*8 cxyz, a, r, p1, q1, r1, dij
      integer ijd
      common/junkx/cxyz(3,5625),
     + a(mxp2),r(mxp2),p1(mxp2),q1(mxp2),r1(mxp2),dij(16*mxp2),
     + ijd(225)

c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnosx/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /miscx/ oianj,okanl,oident,omisc,
     +               oham,opdipd,omp2,ipos1,ipos2
c

      real*8 g
      dimension g(*)

      real*8 e, f1, pp, p, q, aa, ab, ww1, expe
      real*8 csk, bxk, byk, bzk, bk, brrk
      real*8 bbx, bby, bbz, sum, bbrrk, dum, bb, bl
      real*8 bb1, d2, gout

      integer kg, lg, nn, n, lgmax

      real*8 pi252, pie4, dzero, done
      data pi252 /34.986836655250d0/
      data pie4 /7.85398163397448d-01/
      data dzero,done /0.0d0,1.0d0/

      gout = dzero
      lgmax = ngd
      do 40 kg = 1 , ngc
         bk = cgg(kg)
         brrk = bk*rrk
         bxk = bk*pk
         byk = bk*qk
         bzk = bk*rk
         csk = csc(kg)
         if (okanl) lgmax = kg
         do 30 lg = 1 , lgmax
            bl = dg(lg)
            bb = bk + bl
            bb1 = done/bb
            dum = bl*brrk*bb1
            if (dum.le.tol) then
               bbrrk = dum
               d2 = csd(lg)*csk*bb1
               if (okanl .and. lg.ne.kg) d2 = d2 + d2
               bbx = (bxk+bl*pl)*bb1
               bby = (byk+bl*ql)*bb1
               bbz = (bzk+bl*rl)*bb1
               sum = dzero
               nn = 1
               do 20 n = 1 , nij
                  dum = bbrrk + r(n)
                  if (dum.le.tol) then
                     expe = dexp(-dum)
                     aa = a(n)
                     ab = done/(aa+bb)
                     dum = p1(n) - bbx
                     pp = dum*dum
                     dum = q1(n) - bby
                     pp = dum*dum + pp
                     dum = r1(n) - bbz
                     pp = dum*dum + pp
                     p = pp*aa*bb*ab
                     if (p.gt.5.0d0) then
                        if (p.le.15.0d0) then
                           e = dexp(-p)
                           if (p.gt.10.0d0) then
                              ww1 = (((-1.8784686463512d-01/p+
     +                              2.2991849164985d-01)
     +                              /p-4.9893752514047d-01)
     +                              /p-2.1916512131607d-05)
     +                              *e + dsqrt(pie4/p)
                           else
                              ww1 = ((((((4.6897511375022d-01/p-
     +                              6.9955602298985d-01)
     +                              /p+5.3689283271887d-01)
     +                              /p-3.2883030418398d-01)
     +                              /p+2.4645596956002d-01)
     +                              /p-4.9984072848436d-01)
     +                              /p-3.1501078774085d-06)
     +                              *e + dsqrt(pie4/p)
                           end if
                        else if (p.gt.33.0d0) then
                           ww1 = dsqrt(pie4/p)
                        else
                           e = dexp(-p)
                           ww1 = ((1.9623264149430d-01/p-
     +                           4.9695241464490d-01)
     +                           /p-6.0156581186481d-05)
     +                           *e + dsqrt(pie4/p)
                        end if
                     else if (p.gt.1.0d0) then
                        if (p.gt.3.0d0) then
                           q = p - 4.0d0
                           f1 = ((((((((((-2.62453564772299d-11*q+
     +                          3.24031041623823d-10)
     +                          *q-3.614965656163d-09)
     +                          *q+3.760256799971d-08)
     +                          *q-3.553558319675d-07)
     +                          *q+3.022556449731d-06)
     +                          *q-2.290098979647d-05)
     +                          *q+1.526537461148d-04)
     +                          *q-8.81947375894379d-04)
     +                          *q+4.33207949514611d-03)
     +                          *q-1.75257821619926d-02)
     +                          *q + 5.28406320615584d-02
                           ww1 = (p+p)*f1 + dexp(-p)
                        else
                           q = p - 2.0d0
                           f1 = ((((((((((-1.61702782425558d-10*q+
     +                          1.96215250865776d-09)
     +                          *q-2.14234468198419d-08)
     +                          *q+2.17216556336318d-07)
     +                          *q-1.98850171329371d-06)
     +                          *q+1.62429321438911d-05)
     +                          *q-1.16740298039895d-04)
     +                          *q+7.24888732052332d-04)
     +                          *q-3.79490003707156d-03)
     +                          *q+1.61723488664661d-02)
     +                          *q-5.29428148329736d-02)
     +                          *q + 1.15702180856167d-01
                           ww1 = (p+p)*f1 + dexp(-p)
                        end if
                     else if (p.gt.3.0d-07) then
                        f1 = ((((((((-8.36313918003957d-08*p+
     +                       1.21222603512827d-06)
     +                       *p-1.15662609053481d-05)
     +                       *p+9.25197374512647d-05)
     +                       *p-6.40994113129432d-04)
     +                       *p+3.78787044215009d-03)
     +                       *p-1.85185172458485d-02)
     +                       *p+7.14285713298222d-02)
     +                       *p-1.99999999997023d-01)
     +                       *p + 3.33333333333318d-01
                        ww1 = (p+p)*f1 + dexp(-p)
                     else
                        ww1 = 1.0d0 - p/3.0d0
                     end if
                     sum = sum + dij(nn)*ww1*expe*dsqrt(ab)
                  end if
                  nn = nn + 16
 20            continue
               gout = gout + d2*sum
            end if
 30      continue
 40   continue
      g(1) = gout*pi252*qq4
      return
      end
      subroutine xyzint_dft

      implicit none
c
c
      integer in, kn, ni, nj, nk, nl, nmax, mmax
      integer ij1, ij2, kl1, kl2
      real*8 bp01, b00, b10, pcp00, pc00, qcp00, qc00, rcp00, rc00
      real*8 f00, dxij, dyij, dzij, dxkl, dykl, dzkl
      common /setintx/ in(9),kn(9),ni,nj,nk,nl,nmax,mmax,
     +                 bp01,b00,b10,pcp00,pc00,qcp00,qc00,rcp00,rc00,
     +                 f00,dxij,dyij,dzij,dxkl,dykl,dzkl,
     +                 ij1,ij2,kl1,kl2
c
      real*8 pint, qint, rint
      common/junkx/pint(5625),qint(5625),rint(5625)
c
      integer m, n, i1, i2, i3, i4, i5, ia, ib
      integer k1, k2, k3, k4, km
      integer mi, mj, mk, ml, min

      real*8 c01, c10, cp10, cp01
      logical on0, on1, om0, om1

      real*8 dzero, done
      data dzero,done /0.0d0,1.0d0/
c
      on0 = nmax.eq.0
      on1 = nmax.le.1
      om0 = mmax.eq.0
      om1 = mmax.le.1
c
c     ----- in(0,0) -----
c
      i1 = in(1)
      pint(i1) = done
      qint(i1) = done
      rint(i1) = f00
      if (on0 .and. om0) return
      if (.not.(on0)) then
c
c     ----- in(1,0) -----
c
         i2 = in(2)
         pint(i2) = pc00
         qint(i2) = qc00
         rint(i2) = rc00*f00
      end if
      if (.not.(om0)) then
c
c     ----- in(0,1) -----
c
         k2 = kn(2)
         i3 = i1 + k2
         pint(i3) = pcp00
         qint(i3) = qcp00
         rint(i3) = rcp00*f00
         if (.not.(on0)) then
c
c     ----- in(1,1) -----
c
            i3 = i2 + k2
            cp10 = b00
            pint(i3) = pcp00*pint(i2) + cp10
            qint(i3) = qcp00*qint(i2) + cp10
            rint(i3) = rcp00*rint(i2) + cp10*f00
         end if
      end if
      if (.not.(on1)) then
         c10 = dzero
         i3 = i1
         i4 = i2
         do 20 n = 2 , nmax
            c10 = c10 + b10
c
c     ----- in(n,0) -----
c
            i5 = in(n+1)
            pint(i5) = c10*pint(i3) + pc00*pint(i4)
            qint(i5) = c10*qint(i3) + qc00*qint(i4)
            rint(i5) = c10*rint(i3) + rc00*rint(i4)
            if (.not.(om0)) then
               cp10 = cp10 + b00
c
c     ----- in(n,1) -----
c
               i3 = i5 + k2
               pint(i3) = pcp00*pint(i5) + cp10*pint(i4)
               qint(i3) = qcp00*qint(i5) + cp10*qint(i4)
               rint(i3) = rcp00*rint(i5) + cp10*rint(i4)
            end if
            i3 = i4
            i4 = i5
 20      continue
      end if
      if (.not.(om1)) then
         cp01 = dzero
         c01 = b00
         i3 = i1
         i4 = i1 + k2
         do 30 m = 2 , mmax
            cp01 = cp01 + bp01
c
c     ----- in(0,m) -----
c
            i5 = i1 + kn(m+1)
            pint(i5) = cp01*pint(i3) + pcp00*pint(i4)
            qint(i5) = cp01*qint(i3) + qcp00*qint(i4)
            rint(i5) = cp01*rint(i3) + rcp00*rint(i4)
            if (.not.(on0)) then
               c01 = c01 + b00
c
c     ----- in(1,m) -----
c
               i3 = i2 + kn(m+1)
               pint(i3) = pc00*pint(i5) + c01*pint(i4)
               qint(i3) = qc00*qint(i5) + c01*qint(i4)
               rint(i3) = rc00*rint(i5) + c01*rint(i4)
            end if
            i3 = i4
            i4 = i5
 30      continue
      end if
      if (.not.(on1 .or. om1)) then
c
c     ----- in(n,m) -----
c
         c01 = b00
         k3 = k2
         do 50 m = 2 , mmax
            k4 = kn(m+1)
            c01 = c01 + b00
            i3 = i1
            i4 = i2
            c10 = b10
            do 40 n = 2 , nmax
               i5 = in(n+1)
               pint(i5+k4) = c10*pint(i3+k4) + pc00*pint(i4+k4)
     +                       + c01*pint(i4+k3)
               qint(i5+k4) = c10*qint(i3+k4) + qc00*qint(i4+k4)
     +                       + c01*qint(i4+k3)
               rint(i5+k4) = c10*rint(i3+k4) + rc00*rint(i4+k4)
     +                       + c01*rint(i4+k3)
               c10 = c10 + b10
               i3 = i4
               i4 = i5
 40         continue
            k3 = k4
 50      continue
      end if
      if (nj.eq.0) go to 110
c
c     ----- in(mi,mj,m) -----
c
      m = 0
      i5 = in(nmax+1)
 60   min = ni
      km = kn(m+1)
 70   n = nmax
      i3 = i5 + km
 80   i4 = in(n) + km
      pint(i3) = pint(i3) + dxij*pint(i4)
      qint(i3) = qint(i3) + dyij*qint(i4)
      rint(i3) = rint(i3) + dzij*rint(i4)
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 80
      min = min + 1
      if (min.lt.nmax) go to 70
      if (ni.ne.0) then
         i3 = ij2 + i1 + km
         do 100 mj = 1 , nj
            i4 = i3
            do 90 mi = 1 , ni
               pint(i4) = pint(i4+ij1-ij2) + dxij*pint(i4-ij2)
               qint(i4) = qint(i4+ij1-ij2) + dyij*qint(i4-ij2)
               rint(i4) = rint(i4+ij1-ij2) + dzij*rint(i4-ij2)
               i4 = i4 + ij1
 90         continue
            i3 = i3 + ij2
 100     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 60
 110  if (nl.eq.0) go to 170
c
c     ----- in(mi,mj,mk,ml) -----
c
      i5 = kn(mmax+1)
      ia = i1
      mi = 0
 120  mj = 0
      ib = ia
      min = nk
 130  m = mmax
      i3 = ib + i5
 140  i4 = ib + kn(m)
      pint(i3) = pint(i3) + dxkl*pint(i4)
      qint(i3) = qint(i3) + dykl*qint(i4)
      rint(i3) = rint(i3) + dzkl*rint(i4)
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 140
      min = min + 1
      if (min.lt.mmax) go to 130
      if (nk.ne.0) then
         i3 = ib + kl2
         do 160 ml = 1 , nl
            i4 = i3
            do 150 mk = 1 , nk
               pint(i4) = pint(i4+kl1-kl2) + dxkl*pint(i4-kl2)
               qint(i4) = qint(i4+kl1-kl2) + dykl*qint(i4-kl2)
               rint(i4) = rint(i4+kl1-kl2) + dzkl*rint(i4-kl2)
               i4 = i4 + kl1
 150        continue
            i3 = i3 + kl2
 160     continue
      end if
      mj = mj + 1
      ib = ib + ij2
      if (mj.le.nj) then
         min = nk
         go to 130
      else
         mi = mi + 1
         ia = ia + ij1
         if (mi.le.ni) go to 120
      end if
 170  return
      end
      subroutine roots4_dft
c          *****   version february 16,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r14,pie4/1.45303521503316d-01, 7.85398163397448d-01/
      data r24,w24/ 1.33909728812636d+00, 2.34479815323517d-01/
      data r34,w34/ 3.92696350135829d+00, 1.92704402415764d-02/
      data r44,w44/ 8.58863568901199d+00, 2.25229076750736d-04/
      if (pp.gt.15.0d+00) then
         ww1 = dsqrt(pie4/pp)
         if (pp.gt.35.0d+00) then
            if (pp.gt.53.0d+00) then
               rt1 = r14/(pp-r14)
               rt2 = r24/(pp-r24)
               rt3 = r34/(pp-r34)
               rt4 = r44/(pp-r44)
               ww4 = w44*ww1
               ww3 = w34*ww1
               ww2 = w24*ww1
               ww1 = ww1 - ww2 - ww3 - ww4
               return
            else
               g = dexp(-pp)*(pp*pp)**2
c     x=35.0 to 53.0                             nroots = 4
               rt4 = ((-2.19135070169653d-03*pp-1.19108256987623d-01)
     +               *pp-7.50238795695573d-01)*g + r44/(pp-r44)
               rt3 = ((-9.65842534508637d-04*pp-4.49822013469279d-02)
     +               *pp+6.08784033347757d-01)*g + r34/(pp-r34)
               rt2 = ((-3.62569791162153d-04*pp-9.09231717268466d-03)
     +               *pp+1.84336760556262d-01)*g + r24/(pp-r24)
               rt1 = ((-4.07557525914600d-05*pp-6.88846864931685d-04)
     +               *pp+1.74725309199384d-02)*g + r14/(pp-r14)
               ww4 = ((5.76631982000990d-06*pp-7.89187283804890d-05)
     +               *pp+3.28297971853126d-04)*g + w44*ww1
               ww3 = ((2.08294969857230d-04*pp-3.77489954837361d-03)
     +               *pp+2.09857151617436d-02)*g + w34*ww1
               ww2 = ((6.16374517326469d-04*pp-1.26711744680092d-02)
     +               *pp+8.14504890732155d-02)*g + w24*ww1
               ww1 = ww1 - ww2 - ww3 - ww4
               return
            end if
         else if (pp.gt.20.0d+00) then
c     x=20.0 to 35.0                             nroots = 4
            g = dexp(-pp)
          rt1 = ((((((-4.45711399441838d-05*pp+1.27267770241379d-03)*pp-
     +            2.36954961381262d-01)*pp+1.54330657903756d+01)
     +            *pp-5.22799159267808d+02)*pp+1.05951216669313d+04)
     +            *pp+(-2.51177235556236d+06/pp+8.72975373557709d+05)
     +            /pp-1.29194382386499d+05)*g + r14/(pp-r14)
          rt2 = (((((-7.85617372254488d-02*pp+6.35653573484868d+00)*pp-
     +            3.38296938763990d+02)*pp+1.25120495802096d+04)
     +            *pp-3.16847570511637d+05)
     +            *pp+((-1.02427466127427d+09/pp+3.70104713293016d+08)
     +            /pp-5.87119005093822d+07)/pp+5.38614211391604d+06)
     +            *g + r24/(pp-r24)
          rt3 = (((((-2.37900485051067d-01*pp+1.84122184400896d+01)*pp-
     +            1.00200731304146d+03)*pp+3.75151841595736d+04)
     +            *pp-9.50626663390130d+05)
     +            *pp+((-2.88139014651985d+09/pp+1.06625915044526d+09)
     +            /pp-1.72465289687396d+08)/pp+1.60419390230055d+07)
     +            *g + r34/(pp-r34)
          rt4 = ((((((-6.00691586407385d-04*pp-3.64479545338439d-01)*pp+
     +            1.57496131755179d+01)*pp-6.54944248734901d+02)
     +            *pp+1.70830039597097d+04)*pp-2.90517939780207d+05)
     +            *pp+(+3.49059698304732d+07/pp-1.64944522586065d+07)
     +            /pp+2.96817940164703d+06)*g + r44/(pp-r44)
            if (pp.le.25.0d+00) ww4 = (((((((2.33766206773151d-07*pp-
     +                               3.81542906607063d-05)
     +                               *pp+3.51416601267000d-03)
     +                               *pp-1.66538571864728d-01)
     +                               *pp+4.80006136831847d+00)
     +                               *pp-8.73165934223603d+01)
     +                               *pp+9.77683627474638d+02)
     +                               *pp+1.66000945117640d+04/pp-
     +                               6.14479071209961d+03)*g + w44*ww1
            if (pp.gt.25.0d+00) ww4 = ((((((5.74245945342286d-06*pp-
     +                               7.58735928102351d-05)
     +                               *pp+2.35072857922892d-04)
     +                               *pp-3.78812134013125d-03)
     +                               *pp+3.09871652785805d-01)
     +                               *pp-7.11108633061306d+00)
     +                               *pp+5.55297573149528d+01)
     +                               *g + w44*ww1
           ww3 = ((((((2.36392855180768d-04*pp-9.16785337967013d-03)*pp+
     +            4.62186525041313d-01)*pp-1.96943786006540d+01)
     +            *pp+4.99169195295559d+02)*pp-6.21419845845090d+03)
     +            *pp+((+5.21445053212414d+07/pp-1.34113464389309d+07)
     +            /pp+1.13673298305631d+06)/pp-2.81501182042707d+03)
     +            *g + w34*ww1
           ww2 = ((((((7.29841848989391d-04*pp-3.53899555749875d-02)*pp+
     +            2.07797425718513d+00)*pp-1.00464709786287d+02)
     +            *pp+3.15206108877819d+03)*pp-6.27054715090012d+04)
     +            *pp+(+1.54721246264919d+07/pp-5.26074391316381d+06)
     +            /pp+7.67135400969617d+05)*g + w24*ww1
           ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*g + ww1 - ww2 - ww3 - ww4
            return
         else
c     x=15.0 to 20.0                             nroots = 4
            q = pp - 17.5d+00
            rt1 = (((((((((((4.36701759531398d-17*q-1.12860600219889d-16
     +            )*q-6.149849164164d-15)*q+5.820231579541d-14)
     +            *q+4.396602872143d-13)*q-1.24330365320172d-11)
     +            *q+6.71083474044549d-11)*q+2.43865205376067d-10)
     +            *q+1.67559587099969d-08)*q-9.32738632357572d-07)
     +            *q+2.39030487004977d-05)*q-4.68648206591515d-04)
     +            *q + 8.34977776583956d-03
            rt2 = (((((((((((4.98913142288158d-16*q-2.60732537093612d-16
     +            )*q-7.775156445127d-14)*q+5.766105220086d-13)
     +            *q+6.432696729600d-12)*q-1.39571683725792d-10)
     +            *q+5.95451479522191d-10)*q+2.42471442836205d-09)
     +            *q+2.47485710143120d-07)*q-1.14710398652091d-05)
     +            *q+2.71252453754519d-04)*q-4.96812745851408d-03)
     +            *q + 8.26020602026780d-02
            rt3 = (((((((((((1.91498302509009d-15*q+1.48840394311115d-14
     +            )*q-4.316925145767d-13)*q+1.186495793471d-12)
     +            *q+4.615806713055d-11)*q-5.54336148667141d-10)
     +            *q+3.48789978951367d-10)*q-2.79188977451042d-09)
     +            *q+2.09563208958551d-06)*q-6.76512715080324d-05)
     +            *q+1.32129867629062d-03)*q-2.05062147771513d-02)
     +            *q + 2.88068671894324d-01
            rt4 = (((((((((((-5.43697691672942d-15*q-
     +            1.12483395714468d-13)*q+2.826607936174d-12)
     +            *q-1.266734493280d-11)*q-4.258722866437d-10)
     +            *q+9.45486578503261d-09)*q-5.86635622821309d-08)
     +            *q-1.28835028104639d-06)*q+4.41413815691885d-05)
     +            *q-7.61738385590776d-04)*q+9.66090902985550d-03)
     +            *q-1.01410568057649d-01)*q + 9.54714798156712d-01
            ww4 = ((((((((((((-7.56882223582704d-19*q+
     +            7.53541779268175d-18)*q-1.157318032236d-16)
     +            *q+2.411195002314d-15)*q-3.601794386996d-14)
     +            *q+4.082150659615d-13)*q-4.289542980767d-12)
     +            *q+5.086829642731d-11)*q-6.35435561050807d-10)
     +            *q+6.82309323251123d-09)*q-5.63374555753167d-08)
     +            *q+3.57005361100431d-07)*q-2.40050045173721d-06)
     +            *q + 4.94171300536397d-05
            ww3 = (((((((((((-5.54451040921657d-17*q+
     +            2.68748367250999d-16)*q+1.349020069254d-14)
     +            *q-2.507452792892d-13)*q+1.944339743818d-12)
     +            *q-1.29816917658823d-11)*q+3.49977768819641d-10)
     +            *q-8.67270669346398d-09)*q+1.31381116840118d-07)
     +            *q-1.36790720600822d-06)*q+1.19210697673160d-05)
     +            *q-1.42181943986587d-04)*q + 4.12615396191829d-03
            ww2 = (((((((((((-1.86506057729700d-16*q+
     +            1.16661114435809d-15)*q+2.563712856363d-14)
     +            *q-4.498350984631d-13)*q+1.765194089338d-12)
     +            *q+9.04483676345625d-12)*q+4.98930345609785d-10)
     +            *q-2.11964170928181d-08)*q+3.98295476005614d-07)
     +            *q-5.49390160829409d-06)*q+7.74065155353262d-05)
     +            *q-1.48201933009105d-03)*q + 4.97836392625268d-02
            ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*dexp(-pp) + ww1 - ww2 - ww3 -
     +            ww4
            return
         end if
      else if (pp.gt.5.0d+00) then
         if (pp.gt.10.0d+00) then
c     x=10.0 to 15.0                             nroots = 4
            q = pp - 12.5d+00
            rt1 = (((((((((((4.94869622744119d-17*q+8.03568805739160d-16
     +            )*q-5.599125915431d-15)*q-1.378685560217d-13)
     +            *q+7.006511663249d-13)*q+1.30391406991118d-11)
     +            *q+8.06987313467541d-11)*q-5.20644072732933d-09)
     +            *q+7.72794187755457d-08)*q-1.61512612564194d-06)
     +            *q+4.15083811185831d-05)*q-7.87855975560199d-04)
     +            *q + 1.14189319050009d-02
            rt2 = (((((((((((4.89224285522336d-16*q+1.06390248099712d-14
     +            )*q-5.446260182933d-14)*q-1.613630106295d-12)
     +            *q+3.910179118937d-12)*q+1.90712434258806d-10)
     +            *q+8.78470199094761d-10)*q-5.97332993206797d-08)
     +            *q+9.25750831481589d-07)*q-2.02362185197088d-05)
     +            *q+4.92341968336776d-04)*q-8.68438439874703d-03)
     +            *q + 1.15825965127958d-01
            rt3 = ((((((((((6.12419396208408d-14*q+1.12328861406073d-13)
     +            *q-9.051094103059d-12)*q-4.781797525341d-11)
     +            *q+1.660828868694d-09)*q+4.499058798868d-10)
     +            *q-2.519549641933d-07)*q+4.977444040180d-06)
     +            *q-1.25858350034589d-04)*q+2.70279176970044d-03)
     +            *q-3.99327850801083d-02)*q + 4.33467200855434d-01
            rt4 = (((((((((((4.63414725924048d-14*q-4.72757262693062d-14
     +            )*q-1.001926833832d-11)*q+6.074107718414d-11)
     +            *q+1.576976911942d-09)*q-2.01186401974027d-08)
     +            *q-1.84530195217118d-07)*q+5.02333087806827d-06)
     +            *q+9.66961790843006d-06)*q-1.58522208889528d-03)
     +            *q+2.80539673938339d-02)*q-2.78953904330072d-01)
     +            *q + 1.82835655238235d+00
            ww4 = (((((((((((((2.90401781000996d-18*q-
     +            4.63389683098251d-17)*q+6.274018198326d-16)
     +            *q-8.936002188168d-15)*q+1.194719074934d-13)
     +            *q-1.45501321259466d-12)*q+1.64090830181013d-11)
     +            *q-1.71987745310181d-10)*q+1.63738403295718d-09)
     +            *q-1.39237504892842d-08)*q+1.06527318142151d-07)
     +            *q-7.27634957230524d-07)*q+4.12159381310339d-06)
     +            *q-1.74648169719173d-05)*q + 8.50290130067818d-05
            ww3 = ((((((((((((-4.19569145459480d-17*q+
     +            5.94344180261644d-16)*q-1.148797566469d-14)
     +            *q+1.881303962576d-13)*q-2.413554618391d-12)
     +            *q+3.372127423047d-11)*q-4.933988617784d-10)
     +            *q+6.116545396281d-09)*q-6.69965691739299d-08)
     +            *q+7.52380085447161d-07)*q-8.08708393262321d-06)
     +            *q+6.88603417296672d-05)*q-4.67067112993427d-04)
     +            *q + 5.42313365864597d-03
            ww2 = ((((((((((-6.22272689880615d-15*q+1.04126809657554d-13
     +            )*q-6.842418230913d-13)*q+1.576841731919d-11)
     +            *q-4.203948834175d-10)*q+6.287255934781d-09)
     +            *q-8.307159819228d-08)*q+1.356478091922d-06)
     +            *q-2.08065576105639d-05)*q+2.52396730332340d-04)
     +            *q-2.94484050194539d-03)*q + 6.01396183129168d-02
            ww1 = (((-1.8784686463512d-01/pp+2.2991849164985d-01)
     +            /pp-4.9893752514047d-01)/pp-2.1916512131607d-05)
     +            *dexp(-pp) + dsqrt(pie4/pp) - ww4 - ww3 - ww2
            return
         else
            q = pp - 7.5d+00
c     x=5.0 to 10.0                              nroots = 4
            rt1 = (((((((((4.64217329776215d-15*q-6.27892383644164d-15)*
     +            q+3.462236347446d-13)*q-2.927229355350d-11)
     +            *q+5.090355371676d-10)*q-9.97272656345253d-09)
     +            *q+2.37835295639281d-07)*q-4.60301761310921d-06)
     +            *q+8.42824204233222d-05)*q-1.37983082233081d-03)
     +            *q + 1.66630865869375d-02
            rt2 = (((((((((2.93981127919047d-14*q+8.47635639065744d-13)*
     +            q-1.446314544774d-11)*q-6.149155555753d-12)
     +            *q+8.484275604612d-10)*q-6.10898827887652d-08)
     +            *q+2.39156093611106d-06)*q-5.35837089462592d-05)
     +            *q+1.00967602595557d-03)*q-1.57769317127372d-02)
     +            *q + 1.74853819464285d-01
            rt3 = ((((((((((2.93523563363000d-14*q-6.40041776667020d-14)
     +            *q-2.695740446312d-12)*q+1.027082960169d-10)
     +            *q-5.822038656780d-10)*q-3.159991002539d-08)
     +            *q+4.327249251331d-07)*q+4.856768455119d-06)
     +            *q-2.54617989427762d-04)*q+5.54843378106589d-03)
     +            *q-7.95013029486684d-02)*q + 7.20206142703162d-01
            rt4 = (((((((((((-1.62212382394553d-14*q+
     +            7.68943641360593d-13)*q+5.764015756615d-12)
     +            *q-1.380635298784d-10)*q-1.476849808675d-09)
     +            *q+1.84347052385605d-08)*q+3.34382940759405d-07)
     +            *q-1.39428366421645d-06)*q-7.50249313713996d-05)
     +            *q-6.26495899187507d-04)*q+4.69716410901162d-02)
     +            *q-6.66871297428209d-01)*q + 4.11207530217806d+00
            ww1 = ((((((((((-1.65995045235997d-15*q+6.91838935879598d-14
     +            )*q-9.131223418888d-13)*q+1.403341829454d-11)
     +            *q-3.672235069444d-10)*q+6.366962546990d-09)
     +            *q-1.039220021671d-07)*q+1.959098751715d-06)
     +            *q-3.33474893152939d-05)*q+5.72164211151013d-04)
     +            *q-1.05583210553392d-02)*q + 2.26696066029591d-01
            ww2 = ((((((((((((-3.57248951192047d-16*q+
     +            6.25708409149331d-15)*q-9.657033089714d-14)
     +            *q+1.507864898748d-12)*q-2.332522256110d-11)
     +            *q+3.428545616603d-10)*q-4.698730937661d-09)
     +            *q+6.219977635130d-08)*q-7.83008889613661d-07)
     +            *q+9.08621687041567d-06)*q-9.86368311253873d-05)
     +            *q+9.69632496710088d-04)*q-8.14594214284187d-03)
     +            *q + 8.50218447733457d-02
            ww3 = (((((((((((((1.64742458534277d-16*q-
     +            2.68512265928410d-15)*q+3.788890667676d-14)
     +            *q-5.508918529823d-13)*q+7.555896810069d-12)
     +            *q-9.69039768312637d-11)*q+1.16034263529672d-09)
     +            *q-1.28771698573873d-08)*q+1.31949431805798d-07)
     +            *q-1.23673915616005d-06)*q+1.04189803544936d-05)
     +            *q-7.79566003744742d-05)*q+5.03162624754434d-04)
     +            *q-2.55138844587555d-03)*q + 1.13250730954014d-02
            ww4 = ((((((((((((((-1.55714130075679d-17*q+
     +            2.57193722698891d-16)*q-3.626606654097d-15)
     +            *q+5.234734676175d-14)*q-7.067105402134d-13)
     +            *q+8.793512664890d-12)*q-1.006088923498d-10)
     +            *q+1.050565098393d-09)*q-9.91517881772662d-09)
     +            *q+8.35835975882941d-08)*q-6.19785782240693d-07)
     +            *q+3.95841149373135d-06)*q-2.11366761402403d-05)
     +            *q+9.00474771229507d-05)*q-2.78777909813289d-04)
     +            *q + 5.26543779837487d-04
            return
         end if
      else if (pp.gt.1.0d+00) then
c     x= 1.0 to 5.0                              nroots = 4
         q = pp - 3.0d+00
         rt1 = (((((((((-1.48570633747284d-15*q-1.33273068108777d-13)*q+
     +         4.068543696670d-12)*q-9.163164161821d-11)
     +         *q+2.046819017845d-09)*q-4.03076426299031d-08)
     +         *q+7.29407420660149d-07)*q-1.23118059980833d-05)
     +         *q+1.88796581246938d-04)*q-2.53262912046853d-03)
     +         *q + 2.51198234505021d-02
         rt2 = (((((((((1.35830583483312d-13*q-2.29772605964836d-12)*q-
     +         3.821500128045d-12)*q+6.844424214735d-10)
     +         *q-1.048063352259d-08)*q+1.50083186233363d-08)
     +         *q+3.48848942324454d-06)*q-1.08694174399193d-04)
     +         *q+2.08048885251999d-03)*q-2.91205805373793d-02)
     +         *q + 2.72276489515713d-01
         rt3 = (((((((((5.02799392850289d-13*q+1.07461812944084d-11)*q-
     +         1.482277886411d-10)*q-2.153585661215d-09)
     +         *q+3.654087802817d-08)*q+5.15929575830120d-07)
     +         *q-9.52388379435709d-06)*q-2.16552440036426d-04)
     +         *q+9.03551469568320d-03)*q-1.45505469175613d-01)
     +         *q + 1.21449092319186d+00
         rt4 = (((((((((-1.08510370291979d-12*q+6.41492397277798d-11)*q+
     +         7.542387436125d-10)*q-2.213111836647d-09)
     +         *q-1.448228963549d-07)*q-1.95670833237101d-06)
     +         *q-1.07481314670844d-05)*q+1.49335941252765d-04)
     +         *q+4.87791531990593d-02)*q-1.10559909038653d+00)
     +         *q + 8.09502028611780d+00
         ww1 = ((((((((((-4.65801912689961d-14*q+7.58669507106800d-13)*q
     +         -1.186387548048d-11)*q+1.862334710665d-10)
     +         *q-2.799399389539d-09)*q+4.148972684255d-08)
     +         *q-5.933568079600d-07)*q+8.168349266115d-06)
     +         *q-1.08989176177409d-04)*q+1.41357961729531d-03)
     +         *q-1.87588361833659d-02)*q + 2.89898651436026d-01
         ww2 = ((((((((((((-1.46345073267549d-14*q+2.25644205432182d-13)
     +         *q-3.116258693847d-12)*q+4.321908756610d-11)
     +         *q-5.673270062669d-10)*q+7.006295962960d-09)
     +         *q-8.120186517000d-08)*q+8.775294645770d-07)
     +         *q-8.77829235749024d-06)*q+8.04372147732379d-05)
     +         *q-6.64149238804153d-04)*q+4.81181506827225d-03)
     +         *q-2.88982669486183d-02)*q + 1.56247249979288d-01
         ww3 = (((((((((((((9.06812118895365d-15*q-1.40541322766087d-13)
     +         *q+1.919270015269d-12)*q-2.605135739010d-11)
     +         *q+3.299685839012d-10)*q-3.86354139348735d-09)
     +         *q+4.16265847927498d-08)*q-4.09462835471470d-07)
     +         *q+3.64018881086111d-06)*q-2.88665153269386d-05)
     +         *q+2.00515819789028d-04)*q-1.18791896897934d-03)
     +         *q+5.75223633388589d-03)*q-2.09400418772687d-02)
     +         *q + 4.85368861938873d-02
         ww4 = ((((((((((((((-9.74835552342257d-16*q+
     +         1.57857099317175d-14)*q-2.249993780112d-13)
     +         *q+3.173422008953d-12)*q-4.161159459680d-11)
     +         *q+5.021343560166d-10)*q-5.545047534808d-09)
     +         *q+5.554146993491d-08)*q-4.99048696190133d-07)
     +         *q+3.96650392371311d-06)*q-2.73816413291214d-05)
     +         *q+1.60106988333186d-04)*q-7.64560567879592d-04)
     +         *q+2.81330044426892d-03)*q-7.16227030134947d-03)
     +         *q + 9.66077262223353d-03
         return
      else if (pp.gt.3.0d-07) then
c     x=0.0 to 1.0                               nroots = 4
        rt1 = ((((((-1.95309614628539d-10*pp+5.19765728707592d-09)*pp-
     +         1.01756452250573d-07)*pp+1.72365935872131d-06)
     +         *pp-2.61203523522184d-05)*pp+3.52921308769880d-04)
     +         *pp-4.09645850658433d-03)*pp + 3.48198973061469d-02
        rt2 = (((((-1.89554881382342d-08*pp+3.07583114342365d-07)*pp+
     +         1.270981734393d-06)*pp-1.417298563884d-04)
     +         *pp+3.226979163176d-03)*pp-4.48902570678178d-02)
     +         *pp + 3.81567185080039d-01
        rt3 = ((((((1.77280535300416d-09*pp+3.36524958870615d-08)*pp-
     +         2.58341529013893d-07)*pp-1.13644895662320d-05)
     +         *pp-7.91549618884063d-05)*pp+1.03825827346828d-02)
     +         *pp-2.04389090525137d-01)*pp + 1.73730726945889d+00
        rt4 = (((((-5.61188882415248d-08*pp-2.49480733072460d-07)*pp+
     +         3.428685057114d-06)*pp+1.679007454539d-04)
     +         *pp+4.722855585715d-02)*pp-1.39368301737828d+00)
     +         *pp + 1.18463056481543d+01
        ww1 = ((((((-1.14649303201279d-08*pp+1.88015570196787d-07)*pp-
     +         2.33305875372323d-06)*pp+2.68880044371597d-05)
     +         *pp-2.94268428977387d-04)*pp+3.06548909776613d-03)
     +         *pp-3.13844305680096d-02)*pp + 3.62683783378335d-01
        ww2 = ((((((((-4.11720483772634d-09*pp+6.54963481852134d-08)*pp-
     +         7.20045285129626d-07)*pp+6.93779646721723d-06)
     +         *pp-6.05367572016373d-05)*pp+4.74241566251899d-04)
     +         *pp-3.26956188125316d-03)*pp+1.91883866626681d-02)
     +         *pp-8.98046242565811d-02)*pp + 3.13706645877886d-01
        ww3 = ((((((((-3.41688436990215d-08*pp+5.07238960340773d-07)*pp-
     +         5.01675628408220d-06)*pp+4.20363420922845d-05)
     +         *pp-3.08040221166823d-04)*pp+1.94431864731239d-03)
     +         *pp-1.02477820460278d-02)*pp+4.28670143840073d-02)
     +         *pp-1.29314370962569d-01)*pp + 2.22381034453369d-01
        ww4 = (((((((((4.99660550769508d-09*pp-7.94585963310120d-08)*pp+
     +         8.359072409485d-07)*pp-7.422369210610d-06)
     +         *pp+5.763374308160d-05)*pp-3.86645606718233d-04)
     +         *pp+2.18417516259781d-03)*pp-9.99791027771119d-03)
     +         *pp+3.48791097377370d-02)*pp-8.28299075413889d-02)
     +         *pp + 1.01228536290376d-01
         return
      else
c     x is approximately zero.                   nroots = 4
         rt1 = 3.48198973061471d-02 - 4.09645850660395d-03*pp
         rt2 = 3.81567185080042d-01 - 4.48902570656719d-02*pp
         rt3 = 1.73730726945891d+00 - 2.04389090547327d-01*pp
         rt4 = 1.18463056481549d+01 - 1.39368301742312d+00*pp
         ww1 = 3.62683783378362d-01 - 3.13844305713928d-02*pp
         ww2 = 3.13706645877886d-01 - 8.98046242557724d-02*pp
         ww3 = 2.22381034453372d-01 - 1.29314370958973d-01*pp
         ww4 = 1.01228536290376d-01 - 8.28299075414321d-02*pp
         return
      end if
      end
      subroutine roots5_dft
c          *****   version  february 27,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r15,pie4/1.17581320211778d-01, 7.85398163397448d-01/
      data r25,w25/ 1.07456201243690d+00, 2.70967405960535d-01/
      data r35,w35/ 3.08593744371754d+00, 3.82231610015404d-02/
      data r45,w45/ 6.41472973366203d+00, 1.51614186862443d-03/
      data r55,w55/ 1.18071894899717d+01, 8.62130526143657d-06/
      if (pp.gt.15.0d+00) then
         if (pp.gt.25.0d+00) then
            ww1 = dsqrt(pie4/pp)
            if (pp.le.40.0d+00) then
c     x=25.0 to 40.0                             nroots = 5
               g = dexp(-pp)
              rt1 = ((((((((-1.73363958895356d-06*pp+
     +               1.19921331441483d-04)*pp-1.59437614121125d-02)
     +               *pp+1.13467897349442d+00)*pp-4.47216460864586d+01)
     +               *pp+1.06251216612604d+03)*pp-1.52073917378512d+04)
     +               *pp+1.20662887111273d+05)*pp-4.07186366852475d+05)
     +               *g + r15/(pp-r15)
              rt2 = ((((((((-1.60102542621710d-05*pp+
     +               1.10331262112395d-03)*pp-1.50043662589017d-01)
     +               *pp+1.05563640866077d+01)*pp-4.10468817024806d+02)
     +               *pp+9.62604416506819d+03)*pp-1.35888069838270d+05)
     +               *pp+1.06107577038340d+06)*pp-3.51190792816119d+06)
     +               *g + r25/(pp-r25)
              rt3 = ((((((((-4.48880032128422d-05*pp+
     +              2.69025112122177d-03)*pp-4.01048115525954d-01)
     +               *pp+2.78360021977405d+01)*pp-1.04891729356965d+03)
     +               *pp+2.36985942687423d+04)*pp-3.19504627257548d+05)
     +               *pp+2.34879693563358d+06)*pp-7.16341568174085d+06)
     +               *g + r35/(pp-r35)
              rt4 = ((((((((-6.38526371092582d-05*pp-
     +               2.29263585792626d-03)*pp-7.65735935499627d-02)
     +               *pp+9.12692349152792d+00)*pp-2.32077034386717d+02)
     +               *pp+2.81839578728845d+02)*pp+9.59529683876419d+04)
     +               *pp-1.77638956809518d+06)*pp+1.02489759645410d+07)
     +               *g + r45/(pp-r45)
              rt5 = ((((((((-3.59049364231569d-05*pp-
     +               2.25963977930044d-02)*pp+1.12594870794668d+00)
     +               *pp-4.56752462103909d+01)*pp+1.05804526830637d+03)
     +               *pp-1.16003199605875d+04)*pp-4.07297627297272d+04)
     +               *pp+2.22215528319857d+06)*pp-1.61196455032613d+07)
     +               *g + r55/(pp-r55)
              ww5 = (((((((((-4.61100906133970d-10*pp+
     +               1.43069932644286d-07)*pp-1.63960915431080d-05)
     +               *pp+1.15791154612838d-03)*pp-5.30573476742071d-02)
     +               *pp+1.61156533367153d+00)*pp-3.23248143316007d+01)
     +               *pp+4.12007318109157d+02)*pp-3.02260070158372d+03)
     +               *pp+9.71575094154768d+03)*g + w55*ww1
              ww4 = (((((((((-2.40799435809950d-08*pp+
     +               8.12621667601546d-06)*pp-9.04491430884113d-04)
     +               *pp+6.37686375770059d-02)*pp-2.96135703135647d+00)
     +               *pp+9.15142356996330d+01)*pp-1.86971865249111d+03)
     +               *pp+2.42945528916947d+04)*pp-1.81852473229081d+05)
     +               *pp+5.96854758661427d+05)*g + w45*ww1
              ww3 = ((((((((1.83574464457207d-05*pp-1.54837969489927d-03
     +               )*pp+1.18520453711586d-01)*pp-6.69649981309161d+00)
     +               *pp+2.44789386487321d+02)*pp-5.68832664556359d+03)
     +               *pp+8.14507604229357d+04)*pp-6.55181056671474d+05)
     +               *pp+2.26410896607237d+06)*g + w35*ww1
              ww2 = ((((((((2.77778345870650d-05*pp-2.22835017655890d-03
     +               )*pp+1.61077633475573d-01)*pp-8.96743743396132d+00)
     +               *pp+3.28062687293374d+02)*pp-7.65722701219557d+03)
     +               *pp+1.10255055017664d+05)*pp-8.92528122219324d+05)
     +               *pp+3.10638627744347d+06)*g + w25*ww1
              ww1 = ww1 - 0.01962d+00*g - ww2 - ww3 - ww4 - ww5
               return
            else if (pp.gt.59.0d+00) then
c     x=59.0 to infinity                         nroots = 5
               rt1 = r15/(pp-r15)
               rt2 = r25/(pp-r25)
               rt3 = r35/(pp-r35)
               rt4 = r45/(pp-r45)
               rt5 = r55/(pp-r55)
               ww2 = w25*ww1
               ww3 = w35*ww1
               ww4 = w45*ww1
               ww5 = w55*ww1
               ww1 = ww1 - ww2 - ww3 - ww4 - ww5
               return
            else
c     x=40.0 to 59.0                             nroots = 5
               ppp = pp**3
               g = ppp*dexp(-pp)
               rt1 = (((-2.43758528330205d-02*pp+2.07301567989771d+00)
     +               *pp-6.45964225381113d+01)*pp+7.14160088655470d+02)
     +               *g + r15/(pp-r15)
               rt2 = (((-2.28861955413636d-01*pp+1.93190784733691d+01)
     +               *pp-5.99774730340912d+02)*pp+6.61844165304871d+03)
     +               *g + r25/(pp-r25)
               rt3 = (((-6.95053039285586d-01*pp+5.76874090316016d+01)
     +               *pp-1.77704143225520d+03)*pp+1.95366082947811d+04)
     +               *g + r35/(pp-r35)
               rt4 = (((-1.58072809087018d+00*pp+1.27050801091948d+02)
     +               *pp-3.86687350914280d+03)*pp+4.23024828121420d+04)
     +               *g + r45/(pp-r45)
               rt5 = (((-3.33963830405396d+00*pp+2.51830424600204d+02)
     +               *pp-7.57728527654961d+03)*pp+8.21966816595690d+04)
     +               *g + r55/(pp-r55)
               g = ppp*g
               ww5 = ((1.35482430510942d-08*pp-3.27722199212781d-07)
     +               *pp+2.41522703684296d-06)*g + w55*ww1
               ww4 = ((1.23464092261605d-06*pp-3.55224564275590d-05)
     +               *pp+3.03274662192286d-04)*g + w45*ww1
               ww3 = ((1.34547929260279d-05*pp-4.19389884772726d-04)
     +               *pp+3.87706687610809d-03)*g + w35*ww1
               ww2 = ((2.09539509123135d-05*pp-6.87646614786982d-04)
     +               *pp+6.68743788585688d-03)*g + w25*ww1
               ww1 = ww1 - ww2 - ww3 - ww4 - ww5
               return
            end if
         else if (pp.gt.20.0d+00) then
c     x=20.0 to 25.0                             nroots = 5
            q = pp - 22.5d+00
            rt1 = (((((((((-1.13927848238726d-15*q+7.39404133595713d-15)
     +            *q+1.445982921243d-13)*q-2.676703245252d-12)
     +            *q+5.823521627177d-12)*q+2.17264723874381d-10)
     +            *q+3.56242145897468d-09)*q-3.03763737404491d-07)
     +            *q+9.46859114120901d-06)*q-2.30896753853196d-04)
     +            *q + 5.24663913001114d-03
            rt2 = ((((((((((2.89872355524581d-16*q-1.22296292045864d-14)
     +            *q+6.184065097200d-14)*q+1.649846591230d-12)
     +            *q-2.729713905266d-11)*q+3.709913790650d-11)
     +            *q+2.216486288382d-09)*q+4.616160236414d-08)
     +            *q-3.32380270861364d-06)*q+9.84635072633776d-05)
     +            *q-2.30092118015697d-03)*q + 5.00845183695073d-02
            rt3 = ((((((((((1.97068646590923d-15*q-4.89419270626800d-14)
     +            *q+1.136466605916d-13)*q+7.546203883874d-12)
     +            *q-9.635646767455d-11)*q-8.295965491209d-11)
     +            *q+7.534109114453d-09)*q+2.699970652707d-07)
     +            *q-1.42982334217081d-05)*q+3.78290946669264d-04)
     +            *q-8.03133015084373d-03)*q + 1.58689469640791d-01
            rt4 = ((((((((((1.33642069941389d-14*q-1.55850612605745d-13)
     +            *q-7.522712577474d-13)*q+3.209520801187d-11)
     +            *q-2.075594313618d-10)*q-2.070575894402d-09)
     +            *q+7.323046997451d-09)*q+1.851491550417d-06)
     +            *q-6.37524802411383d-05)*q+1.36795464918785d-03)
     +            *q-2.42051126993146d-02)*q + 3.97847167557815d-01
            rt5 = ((((((((((-6.07053986130526d-14*q+1.04447493138843d-12
     +            )*q-4.286617818951d-13)*q-2.632066100073d-10)
     +            *q+4.804518986559d-09)*q-1.835675889421d-08)
     +            *q-1.068175391334d-06)*q+3.292234974141d-05)
     +            *q-5.94805357558251d-04)*q+8.29382168612791d-03)
     +            *q-9.93122509049447d-02)*q + 1.09857804755042d+00
            ww1 = (((((((((-9.10338640266542d-15*q+1.00438927627833d-13)
     +            *q+7.817349237071d-13)*q-2.547619474232d-11)
     +            *q+1.479321506529d-10)*q+1.52314028857627d-09)
     +            *q+9.20072040917242d-09)*q-2.19427111221848d-06)
     +            *q+8.65797782880311d-05)*q-2.82718629312875d-03)
     +            *q + 1.28718310443295d-01
            ww2 = (((((((((5.52380927618760d-15*q-6.43424400204124d-14)*
     +            q-2.358734508092d-13)*q+8.261326648131d-12)
     +            *q+9.229645304956d-11)*q-5.68108973828949d-09)
     +            *q+1.22477891136278d-07)*q-2.11919643127927d-06)
     +            *q+4.23605032368922d-05)*q-1.14423444576221d-03)
     +            *q + 5.06607252890186d-02
            ww3 = (((((((((3.99457454087556d-15*q-5.11826702824182d-14)*
     +            q-4.157593182747d-14)*q+4.214670817758d-12)
     +            *q+6.705582751532d-11)*q-3.36086411698418d-09)
     +            *q+6.07453633298986d-08)*q-7.40736211041247d-07)
     +            *q+8.84176371665149d-06)*q-1.72559275066834d-04)
     +            *q + 7.16639814253567d-03
            ww4 = (((((((((((-2.14649508112234d-18*q-
     +            2.45525846412281d-18)*q+6.126212599772d-16)
     +            *q-8.526651626939d-15)*q+4.826636065733d-14)
     +            *q-3.39554163649740d-13)*q+1.67070784862985d-11)
     +            *q-4.42671979311163d-10)*q+6.77368055908400d-09)
     +            *q-7.03520999708859d-08)*q+6.04993294708874d-07)
     +            *q-7.80555094280483d-06)*q + 2.85954806605017d-04
            ww5 = ((((((((((((-5.63938733073804d-21*q+
     +            6.92182516324628d-20)*q-1.586937691507d-18)
     +            *q+3.357639744582d-17)*q-4.810285046442d-16)
     +            *q+5.386312669975d-15)*q-6.117895297439d-14)
     +            *q+8.441808227634d-13)*q-1.18527596836592d-11)
     +            *q+1.36296870441445d-10)*q-1.17842611094141d-09)
     +            *q+7.80430641995926d-09)*q-5.97767417400540d-08)
     +            *q + 1.65186146094969d-06
            return
         else
c     x=15.0 to 20.0                             nroots = 5
            q = pp - 17.5d+00
            rt1 = ((((((((((1.91875764545740d-16*q+7.8357401095707d-16)*
     +            q-3.260875931644d-14)*q-1.186752035569d-13)
     +            *q+4.275180095653d-12)*q+3.357056136731d-11)
     +            *q-1.123776903884d-09)*q+1.231203269887d-08)
     +            *q-3.99851421361031d-07)*q+1.45418822817771d-05)
     +            *q-3.49912254976317d-04)*q + 6.67768703938812d-03
            rt2 = ((((((((((2.02778478673555d-15*q+1.01640716785099d-14)
     +            *q-3.385363492036d-13)*q-1.615655871159d-12)
     +            *q+4.527419140333d-11)*q+3.853670706486d-10)
     +            *q-1.184607130107d-08)*q+1.347873288827d-07)
     +            *q-4.47788241748377d-06)*q+1.54942754358273d-04)
     +            *q-3.55524254280266d-03)*q + 6.44912219301603d-02
            rt3 = ((((((((((7.79850771456444d-15*q+6.00464406395001d-14)
     +            *q-1.249779730869d-12)*q-1.020720636353d-11)
     +            *q+1.814709816693d-10)*q+1.766397336977d-09)
     +            *q-4.603559449010d-08)*q+5.863956443581d-07)
     +            *q-2.03797212506691d-05)*q+6.31405161185185d-04)
     +            *q-1.30102750145071d-02)*q + 2.10244289044705d-01
            rt4 = (((((((((((-2.92397030777912d-15*q+
     +            1.94152129078465d-14)*q+4.859447665850d-13)
     +            *q-3.217227223463d-12)*q-7.484522135512d-11)
     +            *q+7.19101516047753d-10)*q+6.88409355245582d-09)
     +            *q-1.44374545515769d-07)*q+2.74941013315834d-06)
     +            *q-1.02790452049013d-04)*q+2.59924221372643d-03)
     +            *q-4.35712368303551d-02)*q + 5.62170709585029d-01
            rt5 = (((((((((((1.17976126840060d-14*q+1.24156229350669d-13
     +            )*q-3.892741622280d-12)*q-7.755793199043d-12)
     +            *q+9.492190032313d-10)*q-4.98680128123353d-09)
     +            *q-1.81502268782664d-07)*q+2.69463269394888d-06)
     +            *q+2.50032154421640d-05)*q-1.33684303917681d-03)
     +            *q+2.29121951862538d-02)*q-2.45653725061323d-01)
     +            *q + 1.89999883453047d+00
            ww1 = ((((((((((1.74841995087592d-15*q-6.95671892641256d-16)
     +            *q-3.000659497257d-13)*q+2.021279817961d-13)
     +            *q+3.853596935400d-11)*q+1.461418533652d-10)
     +            *q-1.014517563435d-08)*q+1.132736008979d-07)
     +            *q-2.86605475073259d-06)*q+1.21958354908768d-04)
     +            *q-3.86293751153466d-03)*q + 1.45298342081522d-01
            ww2 = ((((((((((-1.11199320525573d-15*q+1.85007587796671d-15
     +            )*q+1.220613939709d-13)*q+1.275068098526d-12)
     +            *q-5.341838883262d-11)*q+6.161037256669d-10)
     +            *q-1.009147879750d-08)*q+2.907862965346d-07)
     +            *q-6.12300038720919d-06)*q+1.00104454489518d-04)
     +            *q-1.80677298502757d-03)*q + 5.78009914536630d-02
            ww3 = ((((((((((-9.49816486853687d-16*q+6.67922080354234d-15
     +            )*q+2.606163540537d-15)*q+1.983799950150d-12)
     +            *q-5.400548574357d-11)*q+6.638043374114d-10)
     +            *q-8.799518866802d-09)*q+1.791418482685d-07)
     +            *q-2.96075397351101d-06)*q+3.38028206156144d-05)
     +            *q-3.58426847857878d-04)*q + 8.39213709428516d-03
            ww4 = (((((((((((1.33829971060180d-17*q-3.44841877844140d-16
     +            )*q+4.745009557656d-15)*q-6.033814209875d-14)
     +            *q+1.049256040808d-12)*q-1.70859789556117d-11)
     +            *q+2.15219425727959d-10)*q-2.52746574206884d-09)
     +            *q+3.27761714422960d-08)*q-3.90387662925193d-07)
     +            *q+3.46340204593870d-06)*q-2.43236345136782d-05)
     +            *q + 3.54846978585226d-04
            ww5 = (((((((((((((2.69412277020887d-20*q-
     +            4.24837886165685d-19)*q+6.030500065438d-18)
     +            *q-9.069722758289d-17)*q+1.246599177672d-15)
     +            *q-1.56872999797549d-14)*q+1.87305099552692d-13)
     +            *q-2.09498886675861d-12)*q+2.11630022068394d-11)
     +            *q-1.92566242323525d-10)*q+1.62012436344069d-09)
     +            *q-1.23621614171556d-08)*q+7.72165684563049d-08)
     +            *q-3.59858901591047d-07)*q + 2.43682618601000d-06
            return
         end if
      else if (pp.gt.5.0d+00) then
         if (pp.gt.10.0d+00) then
c     x=10.0 to 15.0                             nroots = 5
            q = pp - 12.5d+00
            rt1 = ((((((((((-4.16387977337393d-17*q+7.20872997373860d-16
     +            )*q+1.395993802064d-14)*q+3.660484641252d-14)
     +            *q-4.154857548139d-12)*q+2.301379846544d-11)
     +            *q-1.033307012866d-09)*q+3.997777641049d-08)
     +            *q-9.35118186333939d-07)*q+2.38589932752937d-05)
     +            *q-5.35185183652937d-04)*q + 8.85218988709735d-03
            rt2 = ((((((((((-4.56279214732217d-16*q+6.24941647247927d-15
     +            )*q+1.737896339191d-13)*q+8.964205979517d-14)
     +            *q-3.538906780633d-11)*q+9.561341254948d-11)
     +            *q-9.772831891310d-09)*q+4.240340194620d-07)
     +            *q-1.02384302866534d-05)*q+2.57987709704822d-04)
     +            *q-5.54735977651677d-03)*q + 8.68245143991948d-02
            rt3 = ((((((((((-2.52879337929239d-15*q+2.13925810087833d-14
     +            )*q+7.884307667104d-13)*q-9.023398159510d-13)
     +            *q-5.814101544957d-11)*q-1.333480437968d-09)
     +            *q-2.217064940373d-08)*q+1.643290788086d-06)
     +            *q-4.39602147345028d-05)*q+1.08648982748911d-03)
     +            *q-2.13014521653498d-02)*q + 2.94150684465425d-01
            rt4 = ((((((((((-6.42391438038888d-15*q+5.37848223438815d-15
     +            )*q+8.960828117859d-13)*q+5.214153461337d-11)
     +            *q-1.106601744067d-10)*q-2.007890743962d-08)
     +            *q+1.543764346501d-07)*q+4.520749076914d-06)
     +            *q-1.88893338587047d-04)*q+4.73264487389288d-03)
     +            *q-7.91197893350253d-02)*q + 8.60057928514554d-01
            rt5 = (((((((((((-2.24366166957225d-14*q+
     +            4.87224967526081d-14)*q+5.587369053655d-12)
     +            *q-3.045253104617d-12)*q-1.223983883080d-09)
     +            *q-2.05603889396319d-09)*q+2.58604071603561d-07)
     +            *q+1.34240904266268d-06)*q-5.72877569731162d-05)
     +            *q-9.56275105032191d-04)*q+4.23367010370921d-02)
     +            *q-5.76800927133412d-01)*q + 3.87328263873381d+00
            ww1 = (((((((((8.98007931950169d-15*q+7.25673623859497d-14)*
     +            q+5.851494250405d-14)*q-4.234204823846d-11)
     +            *q+3.911507312679d-10)*q-9.65094802088511d-09)
     +            *q+3.42197444235714d-07)*q-7.51821178144509d-06)
     +            *q+1.94218051498662d-04)*q-5.38533819142287d-03)
     +            *q + 1.68122596736809d-01
            ww2 = ((((((((((-1.05490525395105d-15*q+1.96855386549388d-14
     +            )*q-5.500330153548d-13)*q+1.003849567976d-11)
     +            *q-1.720997242621d-10)*q+3.533277061402d-09)
     +            *q-6.389171736029d-08)*q+1.046236652393d-06)
     +            *q-1.73148206795827d-05)*q+2.57820531617185d-04)
     +            *q-3.46188265338350d-03)*q + 7.03302497508176d-02
            ww3 = (((((((((((3.60020423754545d-16*q-6.24245825017148d-15
     +            )*q+9.945311467434d-14)*q-1.749051512721d-12)
     +            *q+2.768503957853d-11)*q-4.08688551136506d-10)
     +            *q+6.04189063303610d-09)*q-8.23540111024147d-08)
     +            *q+1.01503783870262d-06)*q-1.20490761741576d-05)
     +            *q+1.26928442448148d-04)*q-1.05539461930597d-03)
     +            *q + 1.15543698537013d-02
            ww4 = (((((((((((((2.51163533058925d-18*q-
     +            4.31723745510697d-17)*q+6.557620865832d-16)
     +            *q-1.016528519495d-14)*q+1.491302084832d-13)
     +            *q-2.06638666222265d-12)*q+2.67958697789258d-11)
     +            *q-3.23322654638336d-10)*q+3.63722952167779d-09)
     +            *q-3.75484943783021d-08)*q+3.49164261987184d-07)
     +            *q-2.92658670674908d-06)*q+2.12937256719543d-05)
     +            *q-1.19434130620929d-04)*q + 6.45524336158384d-04
            ww5 = ((((((((((((((-1.29043630202811d-19*q+
     +            2.16234952241296d-18)*q-3.107631557965d-17)
     +            *q+4.570804313173d-16)*q-6.301348858104d-15)
     +            *q+8.031304476153d-14)*q-9.446196472547d-13)
     +            *q+1.018245804339d-11)*q-9.96995451348129d-11)
     +            *q+8.77489010276305d-10)*q-6.84655877575364d-09)
     +            *q+4.64460857084983d-08)*q-2.66924538268397d-07)
     +            *q+1.24621276265907d-06)*q-4.30868944351523d-06)
     +            *q + 9.94307982432868d-06
            return
         else
c     x=5.0 to 10.0                              nroots = 5
            q = pp - 7.5d+00
            rt1 = ((((((((-1.13825201010775d-14*q+1.89737681670375d-13)*
     +            q-4.81561201185876d-12)*q+1.56666512163407d-10)
     +            *q-3.73782213255083d-09)*q+9.15858355075147d-08)
     +            *q-2.13775073585629d-06)*q+4.56547356365536d-05)
     +            *q-8.68003909323740d-04)*q + 1.22703754069176d-02
            rt2 = (((((((((-3.67160504428358d-15*q+1.27876280158297d-14)
     +            *q-1.296476623788d-12)*q+1.477175434354d-11)
     +            *q+5.464102147892d-10)*q-2.42538340602723d-08)
     +            *q+8.20460740637617d-07)*q-2.20379304598661d-05)
     +            *q+4.90295372978785d-04)*q-9.14294111576119d-03)
     +            *q + 1.22590403403690d-01
            rt3 = (((((((((1.39017367502123d-14*q-6.96391385426890d-13)*
     +            q+1.176946020731d-12)*q+1.725627235645d-10)
     +            *q-3.686383856300d-09)*q+2.87495324207095d-08)
     +            *q+1.71307311000282d-06)*q-7.94273603184629d-05)
     +            *q+2.00938064965897d-03)*q-3.63329491677178d-02)
     +            *q + 4.34393683888443d-01
            rt4 = ((((((((((-1.27815158195209d-14*q+1.99910415869821d-14
     +            )*q+3.753542914426d-12)*q-2.708018219579d-11)
     +            *q-1.190574776587d-09)*q+1.106696436509d-08)
     +            *q+3.954955671326d-07)*q-4.398596059588d-06)
     +            *q-2.01087998907735d-04)*q+7.89092425542937d-03)
     +            *q-1.42056749162695d-01)*q + 1.39964149420683d+00
            rt5 = ((((((((((-1.19442341030461d-13*q-2.34074833275956d-12
     +            )*q+6.861649627426d-12)*q+6.082671496226d-10)
     +            *q+5.381160105420d-09)*q-6.253297138700d-08)
     +            *q-2.135966835050d-06)*q-2.373394341886d-05)
     +            *q+2.88711171412814d-06)*q+4.85221195290753d-02)
     +            *q-1.04346091985269d+00)*q + 7.89901551676692d+00
            ww1 = (((((((((7.95526040108997d-15*q-2.48593096128045d-13)*
     +            q+4.761246208720d-12)*q-9.535763686605d-11)
     +            *q+2.225273630974d-09)*q-4.49796778054865d-08)
     +            *q+9.17812870287386d-07)*q-1.86764236490502d-05)
     +            *q+3.76807779068053d-04)*q-8.10456360143408d-03)
     +            *q + 2.01097936411496d-01
            ww2 = (((((((((((1.25678686624734d-15*q-2.34266248891173d-14
     +            )*q+3.973252415832d-13)*q-6.830539401049d-12)
     +            *q+1.140771033372d-10)*q-1.82546185762009d-09)
     +            *q+2.77209637550134d-08)*q-4.01726946190383d-07)
     +            *q+5.48227244014763d-06)*q-6.95676245982121d-05)
     +            *q+8.05193921815776d-04)*q-8.15528438784469d-03)
     +            *q + 9.71769901268114d-02
            ww3 = ((((((((((((-8.20929494859896d-16*q+
     +            1.37356038393016d-14)*q-2.022863065220d-13)
     +            *q+3.058055403795d-12)*q-4.387890955243d-11)
     +            *q+5.923946274445d-10)*q-7.503659964159d-09)
     +            *q+8.851599803902d-08)*q-9.65561998415038d-07)
     +            *q+9.60884622778092d-06)*q-8.56551787594404d-05)
     +            *q+6.66057194311179d-04)*q-4.17753183902198d-03)
     +            *q + 2.25443826852447d-02
            ww4 = ((((((((((((((-1.08764612488790d-17*q+
     +            1.85299909689937d-16)*q-2.730195628655d-15)
     +            *q+4.127368817265d-14)*q-5.881379088074d-13)
     +            *q+7.805245193391d-12)*q-9.632707991704d-11)
     +            *q+1.099047050624d-09)*q-1.15042731790748d-08)
     +            *q+1.09415155268932d-07)*q-9.33687124875935d-07)
     +            *q+7.02338477986218d-06)*q-4.53759748787756d-05)
     +            *q+2.41722511389146d-04)*q-9.75935943447037d-04)
     +            *q + 2.57520532789644d-03
            ww5 = (((((((((((((((7.28996979748849d-19*q-
     +            1.26518146195173d-17)*q+1.886145834486d-16)
     +            *q-2.876728287383d-15)*q+4.114588668138d-14)
     +            *q-5.44436631413933d-13)*q+6.64976446790959d-12)
     +            *q-7.44560069974940d-11)*q+7.57553198166848d-10)
     +            *q-6.92956101109829d-09)*q+5.62222859033624d-08)
     +            *q-3.97500114084351d-07)*q+2.39039126138140d-06)
     +            *q-1.18023950002105d-05)*q+4.52254031046244d-05)
     +            *q-1.21113782150370d-04)*q + 1.75013126731224d-04
            return
         end if
      else if (pp.gt.1.0d+00) then
c     x=1.0 to 5.0                               nroots = 5
         q = pp - 3.0d+00
         rt1 = ((((((((-2.58163897135138d-14*q+8.14127461488273d-13)*q-
     +         2.11414838976129d-11)*q+5.09822003260014d-10)
     +         *q-1.16002134438663d-08)*q+2.46810694414540d-07)
     +         *q-4.92556826124502d-06)*q+9.02580687971053d-05)
     +         *q-1.45190025120726d-03)*q + 1.73416786387475d-02
         rt2 = (((((((((1.04525287289788d-14*q+5.44611782010773d-14)*q-
     +         4.831059411392d-12)*q+1.136643908832d-10)
     +         *q-1.104373076913d-09)*q-2.35346740649916d-08)
     +         *q+1.43772622028764d-06)*q-4.23405023015273d-05)
     +         *q+9.12034574793379d-04)*q-1.52479441718739d-02)
     +         *q + 1.76055265928744d-01
         rt3 = (((((((((-6.89693150857911d-14*q+5.92064260918861d-13)*q+
     +         1.847170956043d-11)*q-3.390752744265d-10)
     +         *q-2.995532064116d-09)*q+1.57456141058535d-07)
     +         *q-3.95859409711346d-07)*q-9.58924580919747d-05)
     +         *q+3.23551502557785d-03)*q-5.97587007636479d-02)
     +         *q + 6.46432853383057d-01
         rt4 = ((((((((-3.61293809667763d-12*q-2.70803518291085d-11)*q+
     +         8.83758848468769d-10)*q+1.59166632851267d-08)
     +         *q-1.32581997983422d-07)*q-7.60223407443995d-06)
     +         *q-7.41019244900952d-05)*q+9.81432631743423d-03)
     +         *q-2.23055570487771d-01)*q + 2.21460798080643d+00
         rt5 = (((((((((7.12332088345321d-13*q+3.16578501501894d-12)*q-
     +         8.776668218053d-11)*q-2.342817613343d-09)
     +         *q-3.496962018025d-08)*q-3.03172870136802d-07)
     +         *q+1.50511293969805d-06)*q+1.37704919387696d-04)
     +         *q+4.70723869619745d-02)*q-1.47486623003693d+00)
     +         *q + 1.35704792175847d+01
         ww1 = (((((((((1.04348658616398d-13*q-1.94147461891055d-12)*q+
     +         3.485512360993d-11)*q-6.277497362235d-10)
     +         *q+1.100758247388d-08)*q-1.88329804969573d-07)
     +         *q+3.12338120839468d-06)*q-5.04404167403568d-05)
     +         *q+8.00338056610995d-04)*q-1.30892406559521d-02)
     +         *q + 2.47383140241103d-01
         ww2 = (((((((((((3.23496149760478d-14*q-5.24314473469311d-13)*q
     +         +7.743219385056d-12)*q-1.146022750992d-10)
     +         *q+1.615238462197d-09)*q-2.15479017572233d-08)
     +         *q+2.70933462557631d-07)*q-3.18750295288531d-06)
     +         *q+3.47425221210099d-05)*q-3.45558237388223d-04)
     +         *q+3.05779768191621d-03)*q-2.29118251223003d-02)
     +         *q + 1.59834227924213d-01
         ww3 = ((((((((((((-3.42790561802876d-14*q+5.26475736681542d-13)
     +         *q-7.184330797139d-12)*q+9.763932908544d-11)
     +         *q-1.244014559219d-09)*q+1.472744068942d-08)
     +         *q-1.611749975234d-07)*q+1.616487851917d-06)
     +         *q-1.46852359124154d-05)*q+1.18900349101069d-04)
     +         *q-8.37562373221756d-04)*q+4.93752683045845d-03)
     +         *q-2.25514728915673d-02)*q + 6.95211812453929d-02
         ww4 = (((((((((((((1.04072340345039d-14*q-1.60808044529211d-13)
     +         *q+2.183534866798d-12)*q-2.939403008391d-11)
     +         *q+3.679254029085d-10)*q-4.23775673047899d-09)
     +         *q+4.46559231067006d-08)*q-4.26488836563267d-07)
     +         *q+3.64721335274973d-06)*q-2.74868382777722d-05)
     +         *q+1.78586118867488d-04)*q-9.68428981886534d-04)
     +         *q+4.16002324339929d-03)*q-1.28290192663141d-02)
     +         *q + 2.22353727685016d-02
         ww5 = ((((((((((((((-8.16770412525963d-16*q+
     +         1.31376515047977d-14)*q-1.856950818865d-13)
     +         *q+2.596836515749d-12)*q-3.372639523006d-11)
     +         *q+4.025371849467d-10)*q-4.389453269417d-09)
     +         *q+4.332753856271d-08)*q-3.82673275931962d-07)
     +         *q+2.98006900751543d-06)*q-2.00718990300052d-05)
     +         *q+1.13876001386361d-04)*q-5.23627942443563d-04)
     +         *q+1.83524565118203d-03)*q-4.37785737450783d-03)
     +         *q + 5.36963805223095d-03
         return
      else if (pp.gt.3.0d-07) then
c     x=0.0 to 1.0                               nroots = 5
        rt1 = ((((((-4.46679165328413d-11*pp+1.21879111988031d-09)*pp-
     +         2.62975022612104d-08)*pp+5.15106194905897d-07)
     +         *pp-9.27933625824749d-06)*pp+1.51794097682482d-04)
     +         *pp-2.15865967920301d-03)*pp + 2.26659266316985d-02
        rt2 = ((((((1.93117331714174d-10*pp-4.57267589660699d-09)*pp+
     +         2.48339908218932d-08)*pp+1.50716729438474d-06)
     +         *pp-6.07268757707381d-05)*pp+1.37506939145643d-03)
     +         *pp-2.20258754419939d-02)*pp + 2.31271692140905d-01
        rt3 = (((((4.84989776180094d-09*pp+1.31538893944284d-07)*pp-
     +         2.766753852879d-06)*pp-7.651163510626d-05)
     +         *pp+4.033058545972d-03)*pp-8.16520022916145d-02)
     +         *pp + 8.57346024118779d-01
        rt4 = ((((-2.48581772214623d-07*pp-4.34482635782585d-06)*pp-
     +         7.46018257987630d-07)*pp+1.01210776517279d-02)
     +         *pp-2.83193369640005d-01)*pp + 2.97353038120345d+00
        rt5 = (((((-8.92432153868554d-09*pp+1.77288899268988d-08)*pp+
     +         3.040754680666d-06)*pp+1.058229325071d-04)
     +         *pp+4.596379534985d-02)*pp-1.75382723579114d+00)
     +         *pp + 1.84151859759049d+01
        ww1 = ((((((-2.03822632771791d-09*pp+3.89110229133810d-08)*pp-
     +         5.84914787904823d-07)*pp+8.30316168666696d-06)
     +         *pp-1.13218402310546d-04)*pp+1.49128888586790d-03)
     +         *pp-1.96867576904816d-02)*pp + 2.95524224714749d-01
        ww2 = (((((((8.62848118397570d-09*pp-1.38975551148989d-07)*pp+
     +         1.602894068228d-06)*pp-1.646364300836d-05)
     +         *pp+1.538445806778d-04)*pp-1.28848868034502d-03)
     +         *pp+9.38866933338584d-03)*pp-5.61737590178812d-02)
     +         *pp + 2.69266719309991d-01
        ww3 = ((((((((-9.41953204205665d-09*pp+1.47452251067755d-07)*pp-
     +         1.57456991199322d-06)*pp+1.45098401798393d-05)
     +         *pp-1.18858834181513d-04)*pp+8.53697675984210d-04)
     +         *pp-5.22877807397165d-03)*pp+2.60854524809786d-02)
     +         *pp-9.71152726809059d-02)*pp + 2.19086362515979d-01
        ww4 = ((((((((-3.84961617022042d-08*pp+5.66595396544470d-07)*pp-
     +         5.52351805403748d-06)*pp+4.53160377546073d-05)
     +         *pp-3.22542784865557d-04)*pp+1.95682017370967d-03)
     +         *pp-9.77232537679229d-03)*pp+3.79455945268632d-02)
     +         *pp-1.02979262192227d-01)*pp + 1.49451349150573d-01
        ww5 = (((((((((4.09594812521430d-09*pp-6.47097874264417d-08)*pp+
     +         6.743541482689d-07)*pp-5.917993920224d-06)
     +         *pp+4.531969237381d-05)*pp-2.99102856679638d-04)
     +         *pp+1.65695765202643d-03)*pp-7.40671222520653d-03)
     +         *pp+2.50889946832192d-02)*pp-5.73782817487958d-02)
     +         *pp + 6.66713443086877d-02
         return
      else
c     x is approximately zero.                   nroots = 5
         rt1 = 2.26659266316985d-02 - 2.15865967920897d-03*pp
         rt2 = 2.31271692140903d-01 - 2.20258754389745d-02*pp
         rt3 = 8.57346024118836d-01 - 8.16520023025515d-02*pp
         rt4 = 2.97353038120346d+00 - 2.83193369647137d-01*pp
         rt5 = 1.84151859759051d+01 - 1.75382723579439d+00*pp
         ww1 = 2.95524224714752d-01 - 1.96867576909777d-02*pp
         ww2 = 2.69266719309995d-01 - 5.61737590184721d-02*pp
         ww3 = 2.19086362515981d-01 - 9.71152726793658d-02*pp
         ww4 = 1.49451349150580d-01 - 1.02979262193565d-01*pp
         ww5 = 6.66713443086877d-02 - 5.73782817488315d-02*pp
         return
      end if
      end
      subroutine rootss_dft
      implicit real*8  (a-h,o-z)
      logical jump
      dimension f(60)
c
      real*8 whi, wlow, rhi, rlow, rfac, amps
      integer ipoint, madd
      common /dft_rtdata/ whi(45),wlow(45),rhi(45),rlow(45),rfac(45),
     +                amps(9),ipoint(9),madd(20)
c
      dimension  beta(12),gamma(12),v(12),d(12)
c
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
c
c      finds weights and points for rys integration
c      using legendre or hermite points as starting
c      estimates.
c
      data m1/1/
      data ms2,m24/-2,24/
      data eps/1.0d-12/
      data d(1)/1.0d0/
c
c
      nroot2 = nroots - 2
      mprime = nroots + nroots - 1
      mprim2 = mprime - 1
      iflag = ipoint(nroots)
      xmax = amps(nroots)
      if (pp.le.1.0d-11) then
         do 20 ipt = 1 , nroots
            u(ipt) = rlow(iflag+ipt)
            w(ipt) = wlow(iflag+ipt)
 20      continue
      else if (pp.ge.xmax) then
         top = 1.0d0/pp
         scale = dsqrt(top)
         do 30 ipt = 1 , nroots
            u(ipt) = rhi(iflag+ipt)*top
            w(ipt) = whi(iflag+ipt)*scale
 30      continue
      else
         call ffun_dft(pp,f,nroots,dji,madd)
         a = f(2)/f(1)
         if (nroots.gt.m1) then
            a1 = a
            do 40 mm = 1 , mprim2
               f(mm+24) = f(mm+2) - a*f(mm+1)
 40         continue
            ibas = 0
            jbas = m24
            max = mprim2
            do 60 k = 2 , nroots
               cfac = f(jbas+1)
               ax = f(jbas+2)/cfac
               bet = a - ax
               beta(k-1) = bet
               a = ax
               gam = -cfac/f(ibas+1)
               gamma(k-1) = gam
               if (k.eq.nroots) go to 70
               max = max + ms2
               do 50 mm = 1 , max
                  f(ibas+mm) = bet*f(jbas+mm+1) + gam*f(ibas+mm+2)
     +                         + f(jbas+mm+2)
 50            continue
               ibas = jbas
               jbas = m24 - jbas
 60         continue
 70         gam = 0.0d0
            do 80 k = 1 , nroots
               r = rlow(iflag+k)/(rfac(iflag+k)*pp+1.0d0)
               if (pp.gt.17.0d0) r = rhi(iflag+k)/pp
               u(k) = r
               gam = gam + r
 80         continue
            epsi = eps*ax
            do 120 mm = 1 , nroots
               rrr = u(mm)
               r = rrr*ax/gam
               jump = .false.
               do 100 icount = 1 , 50
                  v(1) = r - a1
                  top = r + beta(1)
                  dd = v(1) + top
                  d(2) = dd
                  vv = v(1)*top + gamma(1)
                  if (nroot2.gt.0) then
                     v(2) = vv
                     do 90 k = 3 , nroots
                        top = r + beta(k-1)
                        dd = top*dd + gamma(k-1)*d(k-2) + vv
                        d(k) = dd
                        vv = top*vv + gamma(k-1)*v(k-2)
                        v(k) = vv
 90                  continue
                  end if
                  scale = vv/dd
                  r = r - scale
                  if (jump) go to 110
                  if (dabs(scale).lt.epsi) jump = .true.
 100           continue
               write (6,6010) pp , nroots
 110           u(mm) = r
               w(mm) = cfac/(dd*v(nroots-1))
               ax = ax - r
               gam = gam - rrr
 120        continue
         else
            u(1) = a
            w(1) = f(1)
         end if
      end if
      do 130 i = 1 , nroots
         rr = u(i)
         u(i) = rr/(1.0d0-rr)
 130  continue
      return
 6010 format (1x,'no convergence',f15.6,i6)
      end
      subroutine rt123_dft
c             *****   version february 13,1975   *****
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      real*8 pp, u, w, dji
      integer nroots, nroot2
      common /dft_root/ pp,u(12),w(12),nroots,nroot2,dji(60)
c
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r12,pie4/2.75255128608411d-01, 7.85398163397448d-01/
      data r22,w22/ 2.72474487139158d+00, 9.17517095361369d-02/
      data r13/     1.90163509193487d-01/
      data r23,w23/ 1.78449274854325d+00, 1.77231492083829d-01/
      data r33,w33/ 5.52534374226326d+00, 5.11156880411248d-03/
      if (pp.gt.5.0d+00) then
         if (pp.le.15.0d+00) then
            g = dexp(-pp)
            if (pp.gt.10.0d+00) then
c     x = 10.0 to 15.0                 nroots=1,2, or 3
               ww1 = (((-1.8784686463512d-01/pp+2.2991849164985d-01)
     +               /pp-4.9893752514047d-01)/pp-2.1916512131607d-05)
     +               *g + dsqrt(pie4/pp)
               f1 = (ww1-g)/(pp+pp)
               if (nroots.lt.2) then
                  rt1 = f1/(ww1-f1)
                  return
               else if (nroots.eq.2) then
                 rt1 = ((((-1.01041157064226d-05*pp+1.19483054115173d-03
     +                  )*pp-6.73760231824074d-02)
     +                  *pp+1.25705571069895d+00)
     +                  *pp+(((-8.57609422987199d+03/pp+
     +                  5.91005939591842d+03)/pp-1.70807677109425d+03)
     +                /pp+2.64536689959503d+02)/pp-2.38570496490846d+01)
     +                  *g + r12/(pp-r12)
                 rt2 = (((3.39024225137123d-04*pp-9.34976436343509d-02)
     +                  *pp-4.22216483306320d+00)
     +                  *pp+(((-2.08457050986847d+03/pp-
     +                  1.04999071905664d+03)/pp+3.39891508992661d+02)
     +                /pp-1.56184800325063d+02)/pp+8.00839033297501d+00)
     +                  *g + r22/(pp-r22)
                 ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  f2 = (f1+f1+f1-g)/(pp+pp)
                  q = pp - 12.5d+00
                  rt1 = (((((((((((4.42133001283090d-16*q-
     +                  2.77189767070441d-15)*q-4.084026087887d-14)
     +                  *q+5.379885121517d-13)*q+1.882093066702d-12)
     +                  *q-8.67286219861085d-11)*q+7.11372337079797d-10)
     +                  *q-3.55578027040563d-09)*q+1.29454702851936d-07)
     +                  *q-4.14222202791434d-06)*q+8.04427643593792d-05)
     +                  *q-1.18587782909876d-03)
     +                  *q + 1.53435577063174d-02
                  rt2 = (((((((((((6.85146742119357d-15*q-
     +                  1.08257654410279d-14)*q-8.579165965128d-13)
     +                  *q+6.642452485783d-12)*q+4.798806828724d-11)
     +                  *q-1.13413908163831d-09)*q+7.08558457182751d-09)
     +                  *q-5.59678576054633d-08)*q+2.51020389884249d-06)
     +                  *q-6.63678914608681d-05)*q+1.11888323089714d-03)
     +                  *q-1.45361636398178d-02)
     +                  *q + 1.65077877454402d-01
                  rt3 = ((((((((((((3.20622388697743d-15*q-
     +                  2.73458804864628d-14)*q-3.157134329361d-13)
     +                  *q+8.654129268056d-12)*q-5.625235879301d-11)
     +                  *q-7.718080513708d-10)*q+2.064664199164d-08)
     +                  *q-1.567725007761d-07)*q-1.57938204115055d-06)
     +                  *q+6.27436306915967d-05)*q-1.01308723606946d-03)
     +                  *q+1.13901881430697d-02)*q-1.01449652899450d-01)
     +                  *q + 7.77203937334739d-01
                  go to 20
               end if
            else
c     x = 5.0 to 10.0                  nroots =1,2, or 3
             ww1 = ((((((4.6897511375022d-01/pp-6.9955602298985d-01)/pp+
     +               5.3689283271887d-01)/pp-3.2883030418398d-01)
     +               /pp+2.4645596956002d-01)/pp-4.9984072848436d-01)
     +               /pp-3.1501078774085d-06)*g + dsqrt(pie4/pp)
               f1 = (ww1-g)/(pp+pp)
               if (nroots.lt.2) then
                  rt1 = f1/(ww1-f1)
                  return
               else if (nroots.eq.2) then
                  q = pp - 7.5d+00
                  rt1 = (((((((((((((-1.43632730148572d-16*q+
     +                  2.38198922570405d-16)*q+1.358319618800d-14)
     +                  *q-7.064522786879d-14)*q-7.719300212748d-13)
     +                  *q+7.802544789997d-12)*q+6.628721099436d-11)
     +                  *q-1.775564159743d-09)*q+1.713828823990d-08)
     +                  *q-1.497500187053d-07)*q+2.283485114279d-06)
     +                  *q-3.76953869614706d-05)*q+4.74791204651451d-04)
     +                  *q-4.60448960876139d-03)
     +                  *q + 3.72458587837249d-02
                  rt2 = ((((((((((((2.48791622798900d-14*q-
     +                  1.36113510175724d-13)*q-2.224334349799d-12)
     +                  *q+4.190559455515d-11)*q-2.222722579924d-10)
     +                  *q-2.624183464275d-09)*q+6.128153450169d-08)
     +                  *q-4.383376014528d-07)*q-2.49952200232910d-06)
     +                  *q+1.03236647888320d-04)*q-1.44614664924989d-03)
     +                  *q+1.35094294917224d-02)*q-9.53478510453887d-02)
     +                  *q + 5.44765245686790d-01
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  f2 = (f1+f1+f1-g)/(pp+pp)
                  q = pp - 7.5d+00
                  rt1 = (((((((((((5.74429401360115d-16*q+
     +                  7.11884203790984d-16)*q-6.736701449826d-14)
     +                  *q-6.264613873998d-13)*q+1.315418927040d-11)
     +                  *q-4.23879635610964d-11)*q+1.39032379769474d-09)
     +                  *q-4.65449552856856d-08)*q+7.34609900170759d-07)
     +                  *q-1.08656008854077d-05)*q+1.77930381549953d-04)
     +                  *q-2.39864911618015d-03)
     +                  *q + 2.39112249488821d-02
                  rt2 = (((((((((((1.13464096209120d-14*q+
     +                  6.99375313934242d-15)*q-8.595618132088d-13)
     +                  *q-5.293620408757d-12)*q-2.492175211635d-11)
     +                  *q+2.73681574882729d-09)*q-1.06656985608482d-08)
     +                  *q-4.40252529648056d-07)*q+9.68100917793911d-06)
     +                  *q-1.68211091755327d-04)*q+2.69443611274173d-03)
     +                  *q-3.23845035189063d-02)
     +                  *q + 2.75969447451882d-01
                  rt3 = ((((((((((((6.66339416996191d-15*q+
     +                  1.84955640200794d-13)*q-1.985141104444d-12)
     +                  *q-2.309293727603d-11)*q+3.917984522103d-10)
     +                  *q+1.663165279876d-09)*q-6.205591993923d-08)
     +                  *q+8.769581622041d-09)*q+8.97224398620038d-06)
     +                  *q-3.14232666170796d-05)*q-1.83917335649633d-03)
     +                  *q+3.51246831672571d-02)*q-3.22335051270860d-01)
     +                  *q + 1.73582831755430d+00
                  go to 20
               end if
            end if
         else if (pp.gt.33.0d+00) then
c     x = 33.0 to infinity             nroots=1,2, or 3
            ww1 = dsqrt(pie4/pp)
            if (nroots.lt.2) then
               rt1 = 0.5d+00/(pp-0.5d+00)
               return
            else if (nroots.eq.2) then
               if (pp.gt.40.0d+00) then
                  rt1 = r12/(pp-r12)
                  rt2 = r22/(pp-r22)
                  ww2 = w22*ww1
                  ww1 = ww1 - ww2
                  return
               else
                  g = dexp(-pp)
                  rt1 = (-8.78947307498880d-01*pp+1.09243702330261d+01)
     +                  *g + r12/(pp-r12)
                  rt2 = (-9.28903924275977d+00*pp+8.10642367843811d+01)
     +                  *g + r22/(pp-r22)
                  ww2 = (4.46857389308400d+00*pp-7.79250653461045d+01)
     +                  *g + w22*ww1
                  ww1 = ww1 - ww2
                  return
               end if
            else if (pp.gt.47.0d+00) then
               rt1 = r13/(pp-r13)
               rt2 = r23/(pp-r23)
               rt3 = r33/(pp-r33)
               ww2 = w23*ww1
               ww3 = w33*ww1
               ww1 = ww1 - ww2 - ww3
               return
            else
               g = dexp(-pp)
               rt1 = ((-7.39058467995275d+00*pp+3.21318352526305d+02)
     +               *pp-3.99433696473658d+03)*g + r13/(pp-r13)
               rt2 = ((-7.38726243906513d+01*pp+3.13569966333873d+03)
     +               *pp-3.86862867311321d+04)*g + r23/(pp-r23)
               rt3 = ((-2.63750565461336d+02*pp+1.04412168692352d+04)
     +               *pp-1.28094577915394d+05)*g + r33/(pp-r33)
               ww3 = (((1.52258947224714d-01*pp-8.30661900042651d+00)
     +               *pp+1.92977367967984d+02)*pp-1.67787926005344d+03)
     +               *g + w33*ww1
               ww2 = ((6.15072615497811d+01*pp-2.91980647450269d+03)
     +               *pp+3.80794303087338d+04)*g + w23*ww1
               ww1 = ww1 - ww2 - ww3
               return
            end if
         else
c     x = 15.0 to 33.0                 nroots=1,2, or 3
            g = dexp(-pp)
            ww1 = ((1.9623264149430d-01/pp-4.9695241464490d-01)
     +            /pp-6.0156581186481d-05)*g + dsqrt(pie4/pp)
            f1 = (ww1-g)/(pp+pp)
            if (nroots.lt.2) then
               rt1 = f1/(ww1-f1)
               return
            else if (nroots.eq.2) then
             rt1 = ((((-1.14906395546354d-06*pp+1.76003409708332d-04)*pp
     +               -1.71984023644904d-02)*pp-1.37292644149838d-01)
     +               *pp+(-4.75742064274859d+01/pp+9.21005186542857d+00)
     +               /pp-2.31080873898939d-02)*g + r12/(pp-r12)
             rt2 = (((3.64921633404158d-04*pp-9.71850973831558d-02)
     +               *pp-4.02886174850252d+00)
     +               *pp+(-1.35831002139173d+02/pp-8.66891724287962d+01)
     +               /pp+2.98011277766958d+00)*g + r22/(pp-r22)
               ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
               ww1 = ww1 - ww2
               return
            else
               f2 = (f1+f1+f1-g)/(pp+pp)
               if (pp.gt.20.0d+00) then
                 rt1 = ((((-4.97561537069643d-04*pp-5.00929599665316d-02
     +                  )*pp+1.31099142238996d+00)
     +                  *pp-1.88336409225481d+01)
     +                 *pp-6.60344754467191d+02/pp+1.64931462413877d+02)
     +                  *g + r13/(pp-r13)
                 rt2 = ((((-4.48218898474906d-03*pp-5.17373211334924d-01
     +                  )*pp+1.13691058739678d+01)
     +                  *pp-1.65426392885291d+02)
     +                 *pp-6.30909125686731d+03/pp+1.52231757709236d+03)
     +                  *g + r23/(pp-r23)
                 rt3 = ((((-1.38368602394293d-02*pp-1.77293428863008d+00
     +                  )*pp+1.73639054044562d+01)
     +                  *pp-3.57615122086961d+02)
     +                 *pp-1.45734701095912d+04/pp+2.69831813951849d+03)
     +                  *g + r33/(pp-r33)
               else
                 rt1 = ((((((-2.43270989903742d-06*pp+
     +                  3.57901398988359d-04)*pp-2.34112415981143d-02)
     +                *pp+7.81425144913975d-01)*pp-1.73209218219175d+01)
     +                  *pp+2.43517435690398d+02)
     +                  *pp+(-1.97611541576986d+04/pp+
     +                  9.82441363463929d+03)/pp-2.07970687843258d+03)
     +                  *g + r13/(pp-r13)
                 rt2 = (((((-2.62627010965435d-04*pp+
     +                  3.49187925428138d-02)*pp-3.09337618731880d+00)
     +                *pp+1.07037141010778d+02)*pp-2.36659637247087d+03)
     +                  *pp+
     +                  ((-2.91669113681020d+06/pp+1.41129505262758d+06)
     +                /pp-2.91532335433779d+05)/pp+3.35202872835409d+04)
     +                  *g + r23/(pp-r23)
                 rt3 = (((((9.31856404738601d-05*pp-2.87029400759565d-02
     +                  )*pp-7.83503697918455d-01)
     +                *pp-1.84338896480695d+01)*pp+4.04996712650414d+02)
     +                  *pp+
     +                  (-1.89829509315154d+05/pp+5.11498390849158d+04)
     +                  /pp-6.88145821789955d+03)*g + r33/(pp-r33)
               end if
               go to 20
            end if
         end if
      else if (pp.gt.1.0d+00) then
         if (pp.gt.3.0d+00) then
c     x = 3.0 to 5.0                   nroots =1,2, or 3
            q = pp - 4.0d+00
            if (nroots.eq.3) then
               rt1 = (((((((1.44265709189601d-11*q-4.66622033006074d-10)
     +               *q+7.649155832025d-09)*q-1.229940017368d-07)
     +               *q+2.026002142457d-06)*q-2.87048671521677d-05)
     +               *q+3.70326938096287d-04)*q-4.21006346373634d-03)
     +               *q + 3.50898470729044d-02
               rt2 = ((((((((-2.65526039155651d-11*q+
     +               1.97549041402552d-10)*q+2.15971131403034d-09)
     +               *q-7.95045680685193d-08)*q+5.15021914287057d-07)
     +               *q+1.11788717230514d-05)*q-3.33739312603632d-04)
     +               *q+5.30601428208358d-03)*q-5.93483267268959d-02)
     +               *q + 4.31180523260239d-01
               rt3 = ((((((((-3.92833750584041d-10*q-
     +               4.16423229782280d-09)*q+4.42413039572867d-08)
     +               *q+6.40574545989551d-07)*q-3.05512456576552d-06)
     +               *q-1.05296443527943d-04)*q-6.14120969315617d-04)
     +               *q+4.89665802767005d-02)*q-6.24498381002855d-01)
     +               *q + 3.36412312243724d+00
               f2 = ((((((((((-2.36788772599074d-11*q+
     +              2.89147476459092d-10)*q-3.18111322308846d-09)
     +              *q+3.25336816562485d-08)*q-3.00873821471489d-07)
     +              *q+2.48749160874431d-06)*q-1.81353179793672d-05)
     +              *q+1.14504948737066d-04)*q-6.10614987696677d-04)
     +              *q+2.64584212770942d-03)*q-8.66415899015349d-03)
     +              *q + 1.75257821619922d-02
            else
               f1 = ((((((((((-2.62453564772299d-11*q+
     +              3.24031041623823d-10)*q-3.614965656163d-09)
     +              *q+3.760256799971d-08)*q-3.553558319675d-07)
     +              *q+3.022556449731d-06)*q-2.290098979647d-05)
     +              *q+1.526537461148d-04)*q-8.81947375894379d-04)
     +              *q+4.33207949514611d-03)*q-1.75257821619926d-02)
     +              *q + 5.28406320615584d-02
               ww1 = (pp+pp)*f1 + dexp(-pp)
               if (nroots.eq.2) then
                  rt1 = ((((((((-4.11560117487296d-12*q+
     +                  7.10910223886747d-11)*q-1.73508862390291d-09)
     +                  *q+5.93066856324744d-08)*q-9.76085576741771d-07)
     +                  *q+1.08484384385679d-05)*q-1.12608004981982d-04)
     +                  *q+1.16210907653515d-03)*q-9.89572595720351d-03)
     +                  *q + 6.12589701086408d-02
                  rt2 = (((((((((-1.80555625241001d-10*q+
     +                  5.44072475994123d-10)*q+1.603498045240d-08)
     +                  *q-1.497986283037d-07)*q-7.017002532106d-07)
     +                  *q+1.85882653064034d-05)*q-2.04685420150802d-05)
     +                  *q-2.49327728643089d-03)*q+3.56550690684281d-02)
     +                  *q-2.60417417692375d-01)
     +                  *q + 1.12155283108289d+00
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  rt1 = f1/(ww1-f1)
                  return
               end if
            end if
         else
c     x = 1.0 to 3.0                   nroots=1,2, or 3
            q = pp - 2.0d+00
            if (nroots.eq.3) then
               rt1 = ((((((((1.44687969563318d-12*q+4.85300143926755d-12
     +               )*q-6.55098264095516d-10)*q+1.56592951656828d-08)
     +               *q-2.60122498274734d-07)*q+3.86118485517386d-06)
     +               *q-5.13430986707889d-05)*q+6.03194524398109d-04)
     +               *q-6.11219349825090d-03)*q + 4.52578254679079d-02
               rt2 = (((((((6.95964248788138d-10*q-5.35281831445517d-09)
     +               *q-6.745205954533d-08)*q+1.502366784525d-06)
     +               *q+9.923326947376d-07)*q-3.89147469249594d-04)
     +               *q+7.51549330892401d-03)*q-8.48778120363400d-02)
     +               *q + 5.73928229597613d-01
               rt3 = ((((((((-2.81496588401439d-10*q+
     +               3.61058041895031d-09)*q+4.53631789436255d-08)
     +               *q-1.40971837780847d-07)*q-6.05865557561067d-06)
     +               *q-5.15964042227127d-05)*q+3.34761560498171d-05)
     +               *q+5.04871005319119d-02)*q-8.24708946991557d-01)
     +               *q + 4.81234667357205d+00
               f2 = ((((((((((-1.48044231072140d-10*q+
     +              1.78157031325097d-09)*q-1.92514145088973d-08)
     +              *q+1.92804632038796d-07)*q-1.73806555021045d-06)
     +              *q+1.39195169625425d-05)*q-9.74574633246452d-05)
     +              *q+5.83701488646511d-04)*q-2.89955494844975d-03)
     +              *q+1.13847001113810d-02)*q-3.23446977320647d-02)
     +              *q + 5.29428148329709d-02
            else
               f1 = ((((((((((-1.61702782425558d-10*q+
     +              1.96215250865776d-09)*q-2.14234468198419d-08)
     +              *q+2.17216556336318d-07)*q-1.98850171329371d-06)
     +              *q+1.62429321438911d-05)*q-1.16740298039895d-04)
     +              *q+7.24888732052332d-04)*q-3.79490003707156d-03)
     +              *q+1.61723488664661d-02)*q-5.29428148329736d-02)
     +              *q + 1.15702180856167d-01
               ww1 = (pp+pp)*f1 + dexp(-pp)
               if (nroots.eq.2) then
                  rt1 = (((((((((-6.36859636616415d-12*q+
     +                  8.47417064776270d-11)*q-5.152207846962d-10)
     +                  *q-3.846389873308d-10)*q+8.472253388380d-08)
     +                  *q-1.85306035634293d-06)*q+2.47191693238413d-05)
     +                  *q-2.49018321709815d-04)*q+2.19173220020161d-03)
     +                  *q-1.63329339286794d-02)
     +                  *q + 8.68085688285261d-02
                  rt2 = (((((((((1.45331350488343d-10*q+
     +                  2.07111465297976d-09)*q-1.878920917404d-08)
     +                  *q-1.725838516261d-07)*q+2.247389642339d-06)
     +                  *q+9.76783813082564d-06)*q-1.93160765581969d-04)
     +                  *q-1.58064140671893d-03)*q+4.85928174507904d-02)
     +                  *q-4.30761584997596d-01)
     +                  *q + 1.80400974537950d+00
                  ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
                  ww1 = ww1 - ww2
                  return
               else
                  rt1 = f1/(ww1-f1)
                  return
               end if
            end if
         end if
      else if (pp.gt.3.0d-07) then
c     x = 0.0 to 1.0                   nroots=1,2, or 3
         if (nroots.eq.3) then
          rt1 = ((((((-5.10186691538870d-10*pp+2.40134415703450d-08)*pp-
     +            5.01081057744427d-07)*pp+7.58291285499256d-06)
     +            *pp-9.55085533670919d-05)*pp+1.02893039315878d-03)
     +            *pp-9.28875764374337d-03)*pp + 6.03769246832810d-02
          rt2 = ((((((-1.29646524960555d-08*pp+7.74602292865683d-08)*pp+
     +            1.56022811158727d-06)*pp-1.58051990661661d-05)
     +            *pp-3.30447806384059d-04)*pp+9.74266885190267d-03)
     +            *pp-1.19511285526388d-01)*pp + 7.76823355931033d-01
          rt3 = ((((((-9.28536484109606d-09*pp-3.02786290067014d-07)*pp-
     +            2.50734477064200d-06)*pp-7.32728109752881d-06)
     +            *pp+2.44217481700129d-04)*pp+4.94758452357327d-02)
     +            *pp-1.02504611065774d+00)*pp + 6.66279971938553d+00
          f2 = ((((((((-7.60911486098850d-08*pp+1.09552870123182d-06)*pp
     +           -1.03463270693454d-05)*pp+8.16324851790106d-05)
     +           *pp-5.55526624875562d-04)*pp+3.20512054753924d-03)
     +           *pp-1.51515139838540d-02)*pp+5.55555554649585d-02)
     +           *pp-1.42857142854412d-01)*pp + 1.99999999999986d-01
         else
          f1 = ((((((((-8.36313918003957d-08*pp+1.21222603512827d-06)*pp
     +           -1.15662609053481d-05)*pp+9.25197374512647d-05)
     +           *pp-6.40994113129432d-04)*pp+3.78787044215009d-03)
     +           *pp-1.85185172458485d-02)*pp+7.14285713298222d-02)
     +           *pp-1.99999999997023d-01)*pp + 3.33333333333318d-01
            ww1 = (pp+pp)*f1 + dexp(-pp)
            if (nroots.eq.2) then
              rt1 = (((((((-2.35234358048491d-09*pp+2.49173650389842d-08
     +               )*pp-4.558315364581d-08)*pp-2.447252174587d-06)
     +               *pp+4.743292959463d-05)*pp-5.33184749432408d-04)
     +               *pp+4.44654947116579d-03)*pp-2.90430236084697d-02)
     +               *pp + 1.30693606237085d-01
              rt2 = (((((((-2.47404902329170d-08*pp+2.36809910635906d-07
     +               )*pp+1.835367736310d-06)*pp-2.066168802076d-05)
     +               *pp-1.345693393936d-04)*pp-5.88154362858038d-05)
     +               *pp+5.32735082098139d-02)*pp-6.37623643056745d-01)
     +               *pp + 2.86930639376289d+00
               ww2 = ((f1-ww1)*rt1+f1)*(1.0d+00+rt2)/(rt2-rt1)
               ww1 = ww1 - ww2
               return
            else
               rt1 = f1/(ww1-f1)
               return
            end if
         end if
c     x is approximately zero.         nroots=1,2, or 3
      else if (nroots.lt.2) then
         rt1 = 0.5d+00 - pp/5.0d+00
         ww1 = 1.0d+00 - pp/3.0d+00
         return
      else if (nroots.eq.2) then
         rt1 = 1.30693606237085d-01 - 2.90430236082028d-02*pp
         rt2 = 2.86930639376291d+00 - 6.37623643058102d-01*pp
         ww1 = 6.52145154862545d-01 - 1.22713621927067d-01*pp
         ww2 = 3.47854845137453d-01 - 2.10619711404725d-01*pp
         return
      else
         rt1 = 6.03769246832797d-02 - 9.28875764357368d-03*pp
         rt2 = 7.76823355931043d-01 - 1.19511285527878d-01*pp
         rt3 = 6.66279971938567d+00 - 1.02504611068957d+00*pp
         ww1 = 4.67913934572691d-01 - 5.64876917232519d-02*pp
         ww2 = 3.60761573048137d-01 - 1.49077186455208d-01*pp
         ww3 = 1.71324492379169d-01 - 1.27768455150979d-01*pp
         return
      end if
      g = dexp(-pp)
      f1 = ((pp+pp)*f2+g)/3.0d+00
      ww1 = (pp+pp)*f1 + g
 20   t1 = rt1/(rt1+1.0d+00)
      t2 = rt2/(rt2+1.0d+00)
      t3 = rt3/(rt3+1.0d+00)
      a2 = f2 - t1*f1
      a1 = f1 - t1*ww1
      ww3 = (a2-t2*a1)/((t3-t2)*(t3-t1))
      ww2 = (t3*a1-a2)/((t3-t2)*(t2-t1))
      ww1 = ww1 - ww2 - ww3
      return
      end
      subroutine ffun_dft(x,f,npt,dji,madd)
      implicit real*8  (a-h,o-z)
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
      subroutine final_dft

      implicit none
c
c No termination phase yet
c
      return
      end

      subroutine aux_find(tag)
      implicit none
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
      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
      real*8 auxg
      common/dft_auxg/auxg(2000)
      integer tag
      real*8 exad_find
      integer lshl,lprm
      integer pfirst,qfirst
      real*8 ex1,ex2,ex_min
      do lshl=1,nshell(tag)
        pfirst=kstart(tag,lshl)
        ex1=ex_m(tag,pfirst)
        qfirst=pfirst
        do lprm=1,kng(tag,lshl)
          ex2=ex_m(tag,qfirst)
          ex_min=min(ex1,ex2)
          qfirst=qfirst+1
        enddo
        auxg(lshl)=ex_min
      enddo
      return
      end

      integer function num_3c_ints()
      implicit none
c
c     This function computes the maximum number of 3-centre integrals.
c     The assumption is that we compute all the integrals for 
c
c         AO basis                        Fitting basis
c      ishell   jshell <= ishell          kshell
c
c
c      i.e. in AO basis we have a jagged triangle (with full diagnoal
c      blocks) instead of the usual smooth triangle.
c     
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
C *Parameter module for CCP1/DFT            
C *Memory information		
C *----------------
C * max_block 	- maximum number of blocks which can be allocate
      integer max_block
      parameter(max_block=20)
C *							
C *Basis sets information                              
C *----------------------				
C *max_tag 	- maximum number of basis sets which can be alloc
C *max_atype	- maximum number of centre types	
C *max_gtype	- maximum number of grid types (at least every element 
C                 has different grid type)
C *max_grids    - maximum number of grids (different terms may have
C                 different grids, i.e. CPKS equations use a different
C                 grid than the one used for the KS-matrix)
C *max_shel	- maximum number of shells on centre
C *max_prm	- maximum number of primitives for any given centre
C *maxL		- maximum angular momentum allowed	
C *max_func	- maximum number of basis functions for any centre
      integer max_tag,max_atype,max_shel,max_prm,max_ang
      integer max_gtype,max_grids
      parameter(max_tag=3,max_atype=30,max_shel=500,max_prm=5000)
      parameter(max_ang=5,max_gtype=10,max_grids=2)
C *								
C *Geometry information					
C *--------------------					
C *max_atom	- maximum number of atoms in system
      integer max_atom
      parameter(max_atom=750)
C *
C *Accuracy information				
C ---------------------			
C *global_accuracy 	- global accuracy	
      real*8  global_accuracy
      parameter(global_accuracy=1.0d-14)

C *Grid information
c *----------------
c *maxradzn - The maximum number of radial zones. 
      integer maxgpt,maxrad,maxfpt,maxang,maxradzn,maxtablerows
      parameter(maxgpt=2900,maxrad=50,maxang=302,maxfpt=100)
      parameter(maxradzn=35,maxtablerows=7)
c
C *Basis sets include file
c
      real*8  alpha
      real*8  cont_coeff
      integer num_bset,bset_tags
      common/basis_sets/alpha(max_tag,max_atype,max_prm),
     &                  cont_coeff(max_tag,max_atype,max_prm,max_ang),
     &                  num_bset,bset_tags(max_tag)

      integer Ashl, Aprm, Abfn, totshl, totprm, totbfn, size_shlA
      integer size_basA, size_primA, maxi_shlA, maxi_basA, maxi_primA
      integer num_types, atom_tag , num_shl, atm_typ, nprim, angmom
      integer hybrid, pstart

      common/basis_size_info/Ashl(max_tag,max_atype),
     &                       Aprm(max_tag,max_atype),
     &                       Abfn(max_tag,max_atype),
     &                       totshl(max_tag),
     &                       totprm(max_tag),
     &                       totbfn(max_tag),
     &                       size_shlA,maxi_shlA,
     &                       size_basA,maxi_basA,
     &                       size_primA,maxi_primA
C
C Descriptions
C 
C  num_types		-	number of types of atoms for a given basis
C  atom_tag		-	list of atom type id numbers
C
      common/basis_cent_info/num_types(max_tag),
     &                       atom_tag(max_tag,max_atom)
C
C Descriptions
C
C  num_shl		-	number of shells for given atom type
C  atm_typ		-	atomic number of atom type
C  nprim		-  	number of primitives 
C  angmom		-	angular momentum
C  hybrid		-	level of hybrid shell. Is same as angmom if
C				shell is not a hybrid. For hybrid shells is
C				always less than angmom e.g. for an sp shell
C				angmom=2, hybrid=1
C  pstart		-	start of exponents and contraction coeffs 
C				contained in include basis.hf77, for a given
C				shell
C
      common/basis_cont_info/num_shl(max_tag,max_atype),
     &                       atm_typ(max_tag,max_atype),
     &                       nprim(max_tag,max_atype,max_shel),
     &                       angmom(max_tag,max_atype,max_shel),
     &                       hybrid(max_tag,max_atype,max_shel),
     &                       pstart(max_tag,max_atype,max_shel)

      integer mxbas,maxiprm,maxishl
      parameter(mxbas=3,maxiprm=8192,maxishl=3*2048)
c
      real*8 ex_m, cs, cp, cd, cf, cg
      integer kstart, katom, ktype, kng, kloc, kmin, kmax
      integer nshell, non, numorb, ndumm, nbasfn
      common /mbasis/ex_m(mxbas,maxiprm),cs(mxbas,maxiprm),
     +               cp(mxbas,maxiprm),cd(mxbas,maxiprm),
     +               cf(mxbas,maxiprm),cg(mxbas,maxiprm),
     +               kstart(mxbas,maxishl),katom(mxbas,maxishl),
     +               ktype(mxbas,maxishl),kng(mxbas,maxishl),
     +               kloc(mxbas,maxishl),kmin(mxbas,maxishl),
     +               kmax(mxbas,maxishl),
     +               nshell(mxbas),non(mxbas),numorb(mxbas),
     +               ndumm(mxbas),nbasfn(mxbas)
c
c
c     Local:
c
      integer nints   ! the number of integrals
      integer nfunc   ! the number of functions in a shell
      integer ao_tag  ! the tag for the AO basis set
      integer cd_tag  ! the tag for the Coulomb fitting basis set
      integer ishl    ! shell counter
c
c     Code:
c
      ao_tag = bset_tags(1)
      cd_tag = bset_tags(2)
c
c     Rough estimate of AO function pairs
c
      nints = (totbfn(ao_tag)+1)*totbfn(ao_tag)/2
c
c     Correct the estimate by adding the missing diagonal blocks
c
      do ishl = 1, totshl(ao_tag)
        nfunc = kmax(ao_tag,ishl)-kmin(ao_tag,ishl)+1
        nints = nints + (nfunc-1)*nfunc/2
      enddo
c
c     Now multiply by the number of fitting functions
c
      nints = nints * totbfn(cd_tag)
c
      num_3c_ints = nints
      return
      end

      subroutine ver_dft_integ2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/integ2e.m,v $
     +     "/
      data revision /
     +     "$Revision: 6317 $"
     +      /
      data date /
     +     "$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
