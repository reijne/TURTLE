      subroutine dirdef (q,d1,d2,u1,u2,utemp,codens)
c
c----------------------------------------------------------------------
c   Coordinates the calculation of the two electron integrals
c   and their contribution to the matrix u. Adapted from GAMESS
c   routine defunc.
c   (c) Carsten Fuchs 1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character *2 buff
      logical glog1,glog2,glog3,glog5
      logical ext1,ext2,ext3,ext4
      character *8 fnm
      character *6 snm
      data fnm/"dirrpa.m"/
      data snm/"dirdef"/
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
      common/cider /intgrl(8,100),intcon(5,12),ihz(5,5),
     +              ixz,jxz,kxz,lxz,mxz,ijxz,klxz,ninty,ncon,memans,
     +              mempri,mempto,memcto,mtemp,ipsize,icsize,
     +              ipszb,icszb,memcon
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),
     +              nc(maxorb),isatbf(maxat),isatp(maxat),
     +              ipatbf(maxat),ipatp(maxat),
     +              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     +              coprm(mxprim,3),na,nb,ns,np,nsp,n1,npp,n2,n3,np3,
     +              npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8,
     +              n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat),
     +              ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d,
     +              idefs1(numspl),idefs2(numspl),idefp1(numspl),
     +              idefp2(numspl),
     +              idefd1(numspl),idefd2(numspl),ipars,iparp,
     +              ipard,nsspl,npspl,ndspl
      common/gen   /maxit,iter,concrit,conv,energy,thresh,dcmax,
     +              glog1,glog2,glog5,glog3
      common/gen2  /zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     +              zipfal(21,2),zijkal(21,2),
     +              prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      common/tempit/zitest,zipass,zjtest,zjpass

      dimension q(*),d1(*),d2(*)
      dimension u1(nb,nb,maxvec),u2(nb,nb,maxvec),utemp(*),codens(*)
      dimension itesty(3,3,3,3)

      if (odebug(60)) print*,'entering dirdef'
c     sum = 0.0d0
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
c
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
      nbsq = nb*nb
      ifree = igmem_max_memory()
     +      - memreq_pg_dgop(nbsq,'+')
     +      - 2*igmem_overhead()
      k0a = igmem_alloc_inf(ifree,fnm,snm,"k0a",IGMEM_NORMAL)
      if (odebug(60)) print*,'dirdef: new pointer is on ',k0a
c
      if (odebug(60)) print*,'dirdef: ifree is ',ifree
      if (odebug(60)) print*,'dirdef: occupied core is now ',k0a+ifree-1
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
            if (odebug(60)) 
     +          print*,'k00,k01,k02,k03,k04,k05,k06,k0a,iparty = '
            if (odebug(60)) 
     +          print*,k00,k01,k02,k03,k04,k05,k06,k0a,iparty
            ifree1 = ifree - (k06 - k0a)
            if (odebug(60)) print*,'ifree1 = ',ifree1
            ext4 = .false.
            if (odebug(60)) print*,'first call of sav6'
            call sav6(q(k05),q(k0a),q(k00),q(k01),q(k02),q(k03),q(k04),
     +                q(k05),q(k06),iparty,n11,i1,i2,j1,j2,ext1,ext1,
     +                ext4,buff,ifree1)
            if (odebug(60)) print*,'... done'
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
 110              call memlc2(l00,l01,l02,l03,l04,l05,l06,l07,k05,
     +                        iparty)
                  if (odebug(60)) 
     +          print*,'l00,l01,l02,l03,l04,l05,l06,l07,k05,iparty = '
                  if (odebug(60)) 
     +          print*,l00,l01,l02,l03,l04,l05,l06,l07,k05,iparty
                  l08 = l06 + iparty
                  ifreem = ifree - (l08 - k0a)
                  if (odebug(60)) print*,'ifreem = ',ifreem
                  ijxz1 = ixz1 + jxz1 - 1
                  klxz1 = kxz1 + lxz1 - 1
                  mxz1 = ijxz1 + klxz1 - 1
                  time = cpulft(1)
                  if (odebug(60)) print*,'calling tizer'
                  call tizer(q(l06),ifreem,ixz1,jxz1,kxz1,lxz1,ijxz1,
     +                       klxz1,mxz1)
                  if (odebug(60)) print*,'... done'
c
                  ext4 = .true.
                  if (odebug(60)) print*,'second call of sav6'
                  call sav6(q(k05),q(l00),q(l01),q(l02),q(l03),q(l04),
     +                      q(l05),q(l06),q(l07),iparty,n12,k1,k2,l1,l2,
     +                      ext2,ext2,ext4,buff,ifreem)
                  if (odebug(60)) print*,'... done'
                  timing(10) = timing(10) + cpulft(1) - time
c
                  itempy = itesty(ixz1,jxz1,kxz1,lxz1)
                  zitest = 0.0d0
                  zipass = 0.0d0
                  zjtest = 0.0d0
                  zjpass = 0.0d0
c****** calculates 2 electron integrals and adds them to matrix u
                  if (odebug(60)) print*,'calling rpbody'
                  call rpbody
     +			   (q(k05),q(l08),q(l08),q(k0a),q(k00),q(k01),
     +                      q(k02),q(k03),q(k04),n11,q(l00),q(l01),
     +                      q(l02),q(l03),q(l04),q(l05),q(l06),n12,
     +                      d1,d2,u1,u2,utemp,codens,
     +                      ext1,ext2,ext3,i1,i2,j1,j2,k1,
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
c
c     free memory
c
      call gmem_free_inf(k0a,fnm,snm,"k0a")
c
      return
      end
      subroutine dirrpm (ipo,nbasis,dip,c,d1,d2,u1,u2,utemp,scr,q,iq,
     +                   qanal)
c
c----------------------------------------------------------------------
c   Performs Direct-RPA calculation
c   (c) Carsten Fuchs 1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character*80 dmpfil,rstfil
      character*2 char2i
      character*4 name
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      common/blkin /gin(511)
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      common/infoa /nat,ich,mul,nba,nx,ne,ncore
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
c
      real*8 freq, pola
      integer nfreq, ipola, mxfreq
      parameter (mxfreq=100)
      common /polars/ nfreq,ipola(3,3),freq(mxfreq),pola(6,mxfreq)
c
      common/rpacom/modus,nevlo(8),nevhi(8),nev(8),
     +              maxr,maxit,maxrr,maxitt,nbegin,
     +              eps,epss,epsana,epstab,
     +              sytol,tolmy,fetda,ferpa,tottda,totrpa,ipref,iblo13
c
      integer itdspc, irpspc, itdtex, irptex, itable, inorms, ivect
      integer itmfil, ldump, idump, lrest, irest, mrest
      common /rpfile/ itdspc,irpspc,itdtex,irptex,itable,inorms,ivect,
     +                itmfil,ldump,idump,lrest,irest,mrest
c
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      common/table /dmpfil,rstfil,rlines(2)
c
      dimension dip(3*nvc),c(nba,nba)
      dimension d1(*),d2(*),u1(*),u2(*),utemp(*),scr(*),q(*),iq(*)
      dimension qanal(*)
      dimension ipo(nbasis,*)
      character *8 fnm
      character *6 snm
      data fnm/'dirrpa.m'/
      data snm/'dirrpm'/
c
      iblo7 = 1
      name(1:2) = 'tm'
      nav = lenwrd()
c
c   perform symmetry blocked random phase calculation
c
      do 1000 isym = 1 , nirr
         isize = icount(isym)
         neig = nev(isym)
         if (neig .gt. 0) write(iwr,120) isym,isize
120      format (///44('-')/'dimension of rpa matrix for symmetry ',
     ?           i1,': ',i4/44('-')/)
         call flushn(iwr)
         if (isize .eq. 0) go to 1000
         otda = .true.
         if (modus .eq. 3) otda = .false.
c
c   tm file
c
         itmfil = 40 + isym
         name(3:4) = char2i(itmfil)
         if (opg_root()) then
          open (itmfil,file=name,form='unformatted',status='unknown')
         endif
         call tmfile (d1,scr)
c
         if (ototal) then
            inozer = 0
            k = iofsym(isym) + 1
            do 130 icoord = 1 , 3
               if (dnrm2(isize,dip(k),1) .gt. tolmy) then
                  inozer = inozer + 1
                  nozer(inozer) = icoord
               end if
               k = k + nvc
130         continue
            if (odebug(75)) call outive (nozer,inozer,'nozer')
         end if
c
         nstart = nevhi(isym)
         nevlo1 = nevlo(isym) - 1
         isize2 = isize + isize
         isize3 = isize + isize2
c
c        determine core requirements
c
         maxrn = (maxr-1)/nav+1
         neign = (neig-1)/nav+1
         isizen = (isize-1)/nav+1
c
         need =  maxrn +           !   iblb
     +           maxrn +           !   ible
     +           maxrn +           !   iblbb
     +           maxrn +           !   iblee
     +           neign +           !   ibleig
     +           neign             !   iconv
c
         iused  = igmem_alloc_inf(need,fnm,snm,'b+e+bb+ee+eig+conv',
     +                            IGMEM_DEBUG)
c
         need2 = isizen +          !   indx
     +           isize2 +          !   idiag
     +           isize2            !   ibrun

         iused2 = igmem_alloc_inf(need2,fnm,snm,'indx+diag+brun',
     +                            IGMEM_DEBUG)
c
         iblb   = (iused-1)*nav + 1
         ible   = iblb + maxr
         iblbb  = ible + maxr
         iblee  = iblbb + maxr
         ibleig = iblee + maxr
         iconv  = ibleig + neig

         iblock = iblb
         do 307 j = 1 , maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(isize2)
            iblock = iblock + 1
307      continue
         if (odebug(33)) call outive (iq(iblb),maxr,'b blocks')

         iblock = ible
         do 308 j = 1 , maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(isize2)
            iblock = iblock + 1
308      continue
         if (odebug(33)) call outive (iq(ible),maxr,'e blocks')

         iblock = iblbb
         do 309 j = 1 , maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(isize2)
            iblock = iblock + 1
309      continue
         if (odebug(33)) call outive (iq(iblbb),maxr,'bb blocks')

         iblock = iblee
         do 310 j = 1 , maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(isize2)
            iblock = iblock + 1
310      continue
         if (odebug(33)) call outive (iq(iblee),maxr,'ee blocks')

         iblock = ibleig
         do 311 j = 1 , neig
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(isize2)
            iblock = iblock + 1
311      continue
         if (odebug(33)) call outive (iq(ibleig),neig,'eig blocks')

c   get diagonal and search for the NSTART smallest elements

630      continue
c        indx   = iconv + neig
         indx   = (iused2-1)*nav + 1
         indx1  = indx - 1
c        iqneed = indx + isize - iblb
c        idiag  = iused + (iqneed-1)/nav + 2
         idiag = iused2 + isizen
         idiag1 = idiag - 1
         call getdia (.true.,ipo,nbasis,dummy,iq(indx),q(idiag))
         if (neig .eq. 0) go to 800

c   unit starting vectors or restart starting vectors

         if (otda .and. lrest.gt.0 .and. mrest.eq.2) go to 899
         if (iblo13 .gt. 0) go to 790
         iblock = iblb
         if (lrest .gt. 0) then
          if (opg_root()) then
            if (odebug(10)) print*,'opening restore file',rstfil
            open (irest,file=rstfil(1:lrest),form='unformatted',
     +      status='unknown')
            rewind irest
          end if
         end if
         do 1300 i = 1 , nstart
            if (lrest .gt. 0) then
               if (opg_root()) then
               if (odebug(10)) print*,'reading restart vector ',i
               read (irest) (q(idiag1+ii),ii=1,isize2)
               end if
            else
               call vclr(q(idiag),1,isize2)
               q(idiag1+iq(indx1+i)) = 1.0d0
               if (odebug(70)) write (iwr,1299) i,iq(indx1+i)
1299           format ('trial vector ',i3,
     +         ' corresponds to diagonal entry ',i4)
            end if
            call wrt3 (q(idiag),isize2,iq(iblock),num8)
            iblock = iblock + 1
1300     continue
         if (lrest .gt. 0) then
            if (opg_root()) then
             close (irest)
            end if
            lrest = 0
         end if
         go to 800

c   special starting vectors

790      iunit = 14
         jblo13 = 0
1200     jblo13 = jblo13 + 1
         call rdedx (gin,511,jblo13,iunit)
         inew = 0
         do 1210 m = 1 , 511
            if (gin(m) .eq. 0.0d0) then
               if (inew.gt.0 .and. jsym.eq.isym) then
                  call wrt3 (q(idiag),isize2,iq(iblb-1+jroot),num8)
               end if
               inew = 1
               jsym = 0
               go to 1210
            end if
            if (inew .eq. 1) then
               jsym = dnint(gin(m))
               inew = 2
               go to 1210
            end if
            if (inew .eq. 2) then
               jroot = dnint(gin(m))
               inew = 3
               if (jsym .eq. isym) 
     +         call vclr(q(idiag),1,isize2)
               go to 1210
            end if
            if (inew .eq. 3) then
               inew = 4
               if (jsym .ne. isym) go to 1210
               temp = gin(m)
               go to 1210
            end if
            if (inew .eq. 4) then
               inew = 5
               if (jsym .ne. isym) go to 1210
               iitemp = dnint(gin(m))
               go to 1210
            end if
            if (inew .eq. 5) then
               inew = 3
               if (jsym .ne. isym) go to 1210
               iatemp = dnint(gin(m))
               k = ipo(iitemp,iatemp)
               q(idiag1+k) = temp
               go to 1210
            end if
1210     continue
         if (jblo13 .lt. iblo13) go to 1200
         if (odebug(74)) then
            do 1250 j = 1 , nstart
               write (iwr,1220) isym,j
1220           format ('symmetry ',i1,
     +         ', special starting vector ',i3,':')
               call rdedx (q(idiag),isize2,iq(iblb-1+j),num8)
               call outvec (q(idiag),isize,' ')
1250        continue
         end if

c   orthonormalize special starting vectors

         ibnew = idiag
         ibrun = idiag + isize2
         do 1280 j = 1 , nstart
            call rdedx (q(ibnew),isize2,iq(iblb-1+j),num8)
            do 1270 k = 1 , j-1
               call rdedx (q(ibrun),isize2,iq(iblb-1+k),num8)
               t = - ddot(isize,q(ibnew),1,q(ibrun),1)
               call daxpy(isize,t,q(ibrun),1,q(ibnew),1)
1270        continue
            t = dnrm2(isize,q(ibnew),1)
            if (1.0d0+t .eq. 1.0d0) then
               call caserr ('special trial vector must be discarded.')
            end if
            t = 1.0d0 / t
            call dscal(isize,t,q(ibnew),1)
            call wrt3 (q(ibnew),isize2,iq(iblb-1+j),num8)
1280     continue

c   pointers for partitioning of core

800      continue

c        iqneed = indx - iblb
c        ib     = iused + (iqneed-1)/nav + 2
c
         call gmem_free_inf(iused2,fnm,snm,'indx+diag+brun')
c
c        determine core requirements
c
         ib = 0
         isred1 = ib + isize2*maxvec
         iered1 = isred1 + maxrs
         ibred  = iered1 + maxrs
         isred2 = ibred + maxr
         iered2 = isred2 + maxr2s
         ietare = iered2 + maxr2s
         ietaim = ietare + maxr2
         isigma = ietaim + maxr2
         ixred  = isigma + maxr2
         ie     = ixred + maxr2s
         iv     = ie + isize2*maxvec
         ires   = iv + isize2
         need3  = ires + neig
c
c        allocate and assign pointers
c
         ib = igmem_alloc_inf(need3,fnm,snm,'td/rp-algd',IGMEM_DEBUG)
c
         isred1 = ib + isize2*maxvec
         iered1 = isred1 + maxrs
         ibred  = iered1 + maxrs
         isred2 = ibred + maxr
         iered2 = isred2 + maxr2s
         ietare = iered2 + maxr2s
         ietaim = ietare + maxr2
         isigma = ietaim + maxr2
         ixred  = isigma + maxr2
         ie     = ixred + maxr2s
         iv     = ie + isize2*maxvec
         ires   = iv + isize2
c

         idone = 0
         if (neig .eq. 0) then
            if (oall) then
               go to 660
            else
               go to 900
            end if
         end if

c   perform TDA / RPA calculation

         call setsto(neig,0,iq(iconv))
         niter = 0
         call vclr(tim,1,4)
         tim(1) = cpulft(1)
         call vclr(timing,1,10)
         ierr = 0
         if (opg_root()) then
          rewind itable
          rewind inorms
         endif
         if (otda) then
            write (iwr,131) 'TDA'
         else
            write (iwr,131) 'RPA'
         end if
150      if (otda) then
            call tdalgd (ipo,nbasis,dip,c,d1,d2,u1,u2,utemp,scr,
     ?                   q(ib),q(ie),q(iv),
     ?                   iq(iblb),iq(ible),iq(iblbb),iq(iblee),
     ?                   iq(ibleig),iq(iconv),
     ?                   q(iered1),q(ietare),q(ixred),q(ires),q)
         else
            call rpalgd (ipo,nbasis,dip,c,d1,d2,u1,u2,utemp,scr,
     ?                   q(ib),q(ie),q(iv),
     ?                   iq(iblb),iq(ible),iq(iblbb),iq(iblee),
     ?                   iq(ibleig),iq(iconv),
     ?                   q(isred1),q(iered1),q(ibred),q(isred2),
     +                   q(iered2),
     ?                   q(ietare),q(ietaim),q(isigma),q(ixred),
     ?                   q(ires),q)
         end if
         nconv = isum(neig,iq(iconv),1)
         if (nconv.lt.neig .and. ierr.eq.0 .and. niter.lt.maxit) then
            call restrt (otda,q(ib),q(ie),q(iv),
     ?                   iq(iblb),iq(ible),iq(iblbb),iq(iblee))
            go to 150
         end if

c   output of convergence behaviour

         tim(4) = cpulft(1) - tim(1)
         write (iwr,132)
132      format ('iteration  reduced    eigenvalues'/
     ?           'number     dimension'/89('-'))
         call outeig (itable,iwr,dummy,neig,niter,ntot,nevlo1,2)
         if (nconv .lt. neig) then
            write (iwr,1310)
1310        format 
     +      (/'maximum number of iterations reached, no convergence.'/)
         else if (niter .eq. 1) then
            write (iwr,100)
100         format (//89('*')/'*',10x,'all eigenvalues have converged ',
     ?              'after one iteration.',26x,'*'/89('*')//)
         else
            write (iwr,201) niter
201         format (//89('*')/'*',10x,
     +                   'all eigenvalues have converged after ',
     ?              i2,' iterations.',26x,'*'/89('*')//)
         end if 
         write (iwr,133)
133      format (/'iteration             norms of residue vectors'
     ?           /'number'/80('-'))
         call outres (iwr,inorms,dummy,nevlo1,neig,niter,2)

c   write table of results

         if (nconv .gt. 0) then
            ienerg = ietare + nevlo1
            oincor = .false.
            call rpanal (ipo,nbasis,q(ienerg),q(iv),dummy,q(ietaim),
     +                   q(isigma),
     +                   dip(iofsym(isym)+1),iq(ibleig),iq(iconv),
     +                   otda,oincor,2,qanal)
         end if

c   timing analysis

         write (iwr,205) tim(2),tim(3),tim(4)-tim(2)-tim(3),tim(4)
205      format (/'timing analysis:'/
     ?  'multiplication of E with trial vectors: ',f14.4,' seconds'/
     ?  'solution of reduced eigenvalue problems:',f14.4,' seconds'/
     ?  'remaining tasks:                        ',f14.4,' seconds'/
     ?  'total iterative procedure:              ',f14.4,' seconds')
         timing(7) = timing(3) - timing(6)
         timing(8) = tim(2)
         do 206 itime = 1 , 5
            timing(8) = timing(8) - timing(itime)
206      continue
         fact = 100.0d0 / tim(2)
         write (iwr,207) timing(1) * fact,
     ?                  timing(2) * fact,
     ?                  timing(7) * fact,
     ?                  timing(6) * fact,
     ?                  timing(4) * fact,
     ?                  timing(5) * fact,
     ?                  timing(8) * fact
207      format (/'analysis of multiplication routine:'/
     ?            'construction of d matrices:             ',f4.1,' %'/
     ?            'weighting of d and compressing:         ',f4.1,' %'/
     ?            'integral evaluation:                    ',f4.1,' %'/
     ?            'integral contribution to u matrices:    ',f4.1,' %'/
     ?            'weighting and symmetrization of u:      ',f4.1,' %'/
     ?            'construction of e,f:                    ',f4.1,' %'/
     ?            'remaining tasks:                        ',f4.1,' %'/)
         if(odebug(60)) write (iwr,208) 100.0d0*timing(10)/timing(7)
208      format (/'tizer/sav6 took ',f4.1,' % of integral evaluation'/)

c   restart ?

899      if (modus .eq. 1) go to 660
         if (modus .eq .2 .and. otda) then
            otda = .false.
            call gmem_free_inf(ib,fnm,snm,'td/rp-algd')
            iused2 = igmem_alloc_inf(need2,fnm,snm,'indx+diag+brun',
     &                               IGMEM_DEBUG)
            go to 630
         end if

c   compute total oscillator strength (see C. Fuchs, PhD thesis)

660      if (ototal .and. idone.lt.inozer) then
            if (inozer-idone .le. maxvec) then
               ibb = ib
               do 930 icoord = idone+1 , inozer
                  imy = iofsym(isym) + (nozer(icoord)-1)*nvc + 1
                  call dcopy(isize,dip(imy),1,q(ibb),1)
                  if (modus .ge. 2) 
     +     call dcopy(isize,dip(imy),1,q(ibb+isize),1)
                  ibb = ibb + isize2
930            continue
               mvec = inozer - idone
               if (modus .le. 2) then
                  call rpmuld (.true.,ipo,nbasis,c,d1,d2,u1,u2,
     +                         utemp,scr,q(ib),q(ie),q)
                  ibb = ib
                  iee = ie
                  do 940 icoord = idone+1 , inozer
                     tottda = tottda + ddot(isize,q(ibb),1,q(iee),1)
                     ibb = ibb + isize2
                     iee = iee + isize2
940               continue
               end if
               if (modus .ge. 2) then
                  call rpmuld (.false.,ipo,nbasis,c,d1,d2,u1,u2,
     +                         utemp,scr,q(ib),q(ie),q)
                  ibb = ib
                  iee = ie
                  do 950 icoord = idone+1 , inozer
                     if (odebug(75)) 
     +               call outvec (q(ibb),isize,'vector my')
                     totrpa = totrpa + ddot(isize2,q(ibb),1,q(iee),1)
                     ibb = ibb + isize2
                     iee = iee + isize2
950               continue
               end if
            else
               write (iwr,960) inozer-idone,maxvec
960            format(//'WARNING: cannot complete calculation of total'/
     ?                   'oscillator strength because there are more'/
     ?                   'remaining coordinates (',i1,') than trial'/
     ?                   'vectors allowed in one batch (',i1,').'//)
            end if
         end if
900      continue
         if (nfreq .gt. 0) then

c   polarizability calculation

           print*,'polarizability not yet implemented -- sorry'
           write (iwr,910)
910        format (//122('=')//)
         endif

c   end loop over symmetries

         if (opg_root()) then
          close (itmfil)
         endif
         call gmem_free_inf(ib,fnm,snm,'td/rp-algd')
         call gmem_free_inf(iused,fnm,snm,'b+e+bb+ee+eig+conv')
1000  continue

      return
131   format (///'iterative ',a3,' procedure:'/24('*')//)
      end
c
      subroutine rpalgd (ipo,nbasis,dip,c,d1,d2,u1,u2,utemp,scr,x,e,v,
     +                         iblb,ible,iblbb,iblee,ibleig,iconv,
     +                         sred1,ered1,bred,sred2,ered2,
     +                             etare,etaim,sigma,xred,
     +                   resarr,q)
c
c----------------------------------------------------------------------
c   Performs the iterative solution of the generalized RPA eigenvalue
c   equation
c                 
c       (         |  1   0  |   |  A   B  | ) | y |   | 0 |
c       ( omega * |         | - |         | ) |   | = |   |
c       (         |  0  -1  |   |  B   A  | ) | z |   | 0 |
c
c   using the iterative algorithm by Olsen & Jorgensen 
c   (J. Comp. Phys. 74, 265 (1988))
c   Adapted from routine rpalgo for the direct implementation.
c   (c) Carsten Fuchs 1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character*80 dmpfil,rstfil
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
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
c     A common block to pass the actual time periods around in the
c     response theory code. This is needed because in subroutine 
c     such as rpmuld and rpbody the time periods that are used 
c     should be different for RPA and TDA calculations. Hence the
c     usual parameters cannot be used directly.
c
      integer TP_RESP_MO2AO
      integer TP_RESP_INT
      integer TP_RESP_CNTRCT
      integer TP_RESP_AO2MO
      common/response_timep/TP_RESP_MO2AO,TP_RESP_INT,TP_RESP_CNTRCT,
     +                      TP_RESP_AO2MO
c
c
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,newbas,ncol,
     ?              ivalue,ioccup,ispa
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     ?              lcore(8),lvirt(8),iofset(8,8),npair(8),
     ?              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
      common/infoa /nat,ich,mul,nba,nx,ne,ncore
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/rpacom/modus,nevlo(8),nevhi(8),nev(8),
     ?              maxr,maxit,maxrr,maxitt,nbegin,
     +              eps,epss,epsana,epstab,
     ?              sytol,tolmy,fetda,ferpa,tottda,totrpa,ipref,iblo13
c
      integer itdspc, irpspc, itdtex, irptex, itable, inorms, ivect
      integer itmfil, ldump, idump, lrest, irest, mrest
      common /rpfile/ itdspc,irpspc,itdtex,irptex,itable,inorms,ivect,
     +                itmfil,ldump,idump,lrest,irest,mrest
c
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      common/table /dmpfil,rstfil,rlines(2)

c
      dimension dip(3*nvc),c(nba,nba),ipo(nbasis,*)
      dimension d1(*),d2(*),u1(*),u2(*),utemp(*),scr(*)
      dimension x(isize2,maxvec),e(isize2,maxvec),v(isize2)
      dimension iblb(maxr),ible(maxr),iblbb(maxr),iblee(maxr),
     ?          ibleig(neig),iconv(neig)
      dimension sred1(maxr,maxr),ered1(maxr,maxr),bred(maxr)
      dimension sred2(maxr2,maxr2),ered2(maxr2,maxr2)
      dimension etare(maxr2),etaim(maxr2),sigma(maxr2),xred(maxr2,maxr2)
      dimension resarr(*),q(*)
c
      data small/1.0d-9/
      data ieigen,iresid / 20,21 /
c
      call start_time_period(TP_RPA)
      TP_RESP_MO2AO  = TP_RPA_MO2AO
      TP_RESP_INT    = TP_RPA_INT
      TP_RESP_CNTRCT = TP_RPA_CNTRCT
      TP_RESP_AO2MO  = TP_RPA_AO2MO
c
      ntot1 = 1
      ntot = nstart
      n1 = 1
      n2 = isize + 1
10    niter = niter + 1
      orest = oskip .and. ntot.eq.nstart .and. niter.gt.1
      if (niter .gt. maxit) go to 1000
      if (ldump .gt. 0) then
        if (opg_root()) then
         open (idump,file=dmpfil(1:ldump),form='unformatted',
     +         status = 'unknown')
         rewind idump
        end if
      end if

c   calculate the vectors e * b(i) and fill up matrices sred1,ered1
c
c      m0    ... index of first trial vector in the current batch
c      m1    ... index of last trial vector in the current batch
c      mvec  ... number of trial vectors in the current batch
c      m,n   ... true indices of trial vectors
c      i,j   ... temporary indices of trial vectors
c                 (with respect to the current batch)

      m0 = ntot1
20    m1 = min(m0+maxvec-1,ntot)
      mvec = m1 - m0 + 1
      m = m0
      do 30 i = 1 , mvec
         if (odebug(43)) write (iwr,1100) 30,'trial',m,iblb(m)
         call rdedx (x(1,i),isize2,iblb(m),num8)
         if (odebug(30)) then
            print*,'trial vector ',m,':'
            call outvec (x(1,i),isize2,'x')
         end if
         if (orest) then
            call rdedx (e(1,i),isize2,ible(m),num8)
         end if
         m = m + 1
30    continue
      if (orest) go to 35

c   the next piece of code adds the dipole moment vector(s) to the
c   trial vector(s) if the total oscillator strength is to be calculated

      mvec0 = mvec
      if (ototal .and. idone.lt.inozer .and. mvec0.lt.maxvec) then
         if (odebug(77)) write (iwr,301) isym,m0,m1,mvec0,maxvec
301      format ('additional space in symmetry ',i1/
     ?           'when calculating trial vectors ',i2,' to ',i2/
     ?           'mvec is ',i2,', and maxvec is ',i4)
         do 32 mv = 1 , maxvec-mvec0
            mvec = mvec + 1
            idone = idone + 1
            if (odebug(77)) write (iwr,302) nozer(idone),mvec
302         format ('adding my vector ',i1,' as trial vector ',i2)
            ibegin = iofsym(isym) + (nozer(idone)-1)*nvc + 1
            call dcopy(isize,dip(ibegin),1,x(n1,mvec),1)
            call dcopy(isize,dip(ibegin),1,x(n2,mvec),1)
            if (idone .eq. inozer) go to 33
32       continue
33       continue
      end if

c   compute product vectors

      time = cpulft(1)
      call rpmuld (.false.,ipo,nbasis,c,d1,d2,u1,u2,utemp,scr,x,e,q)
      tim(2) = tim(2) + cpulft(1) - time

c   again worry about total oscillator strength

      if (mvec .gt. mvec0) then
         do 34 mv = mvec0+1 , mvec
            totrpa = totrpa + ddot(isize2,x(1,mv),1,e(1,mv),1)
34       continue
      end if
      mvec = mvec0

c   fill up the matrices sred1 & ered1

35    m = m0
      do 60 i = 1 , mvec
         sred1(m,m) = ddot(isize,x(n1,i),1,x(n1,i),1) - 
     ?                ddot(isize,x(n2,i),1,x(n2,i),1)
         ered1(m,m) = ddot(isize,x(n1,i),1,e(n1,i),1) + 
     ?                ddot(isize,x(n2,i),1,e(n2,i),1)
         bred(m)    = ddot(isize,x(n1,i),1,e(n2,i),1) + 
     ?                ddot(isize,x(n2,i),1,e(n1,i),1)
         if (odebug(17)) then
            print*,'sred1 & ered1 for i = ',i
            call outsqr (sred1,maxr,ntot,ntot,'sred1')
            call outsqr (ered1,maxr,ntot,ntot,'ered1')
         end if
         if (odebug(30)) then
            print*,'ex for trial vector ',m,':'
            call outvec (e(1,i),isize2,'ex')
         end if
         if (.not. orest) call wrt3 (e(1,i),isize2,ible(m),num8)
         do 40 n = 1, m0-1
            if (odebug(43)) write (iwr,1100) 40,'trial',n,iblb(n)
            call rdedx (v,isize2,iblb(n),num8)
            ered1(n,m) = ddot(isize,v(n1),1,e(n1,i),1) + 
     ?                   ddot(isize,v(n2),1,e(n2,i),1)
            ered1(m,n) = ddot(isize,v(n1),1,e(n2,i),1) + 
     ?                   ddot(isize,v(n2),1,e(n1,i),1)
            sred1(n,m) = ddot(isize,v(n1),1,x(n1,i),1) - 
     ?                   ddot(isize,v(n2),1,x(n2,i),1)
            sred1(m,n) = ddot(isize,v(n2),1,x(n1,i),1) - 
     ?                         ddot(isize,v(n1),1,x(n2,i),1)
40       continue
	 n = m0
         do 50 j = 1 , i-1
            ered1(n,m) = ddot(isize,x(n1,j),1,e(n1,i),1) + 
     ?			 ddot(isize,x(n2,j),1,e(n2,i),1)
            ered1(m,n) = ddot(isize,x(n1,j),1,e(n2,i),1) + 
     ?			 ddot(isize,x(n2,j),1,e(n1,i),1)
            sred1(n,m) = ddot(isize,x(n1,j),1,x(n1,i),1) - 
     ?			 ddot(isize,x(n2,j),1,x(n2,i),1)
            sred1(m,n) = ddot(isize,x(n2,j),1,x(n1,i),1) - 
     ?			 ddot(isize,x(n1,j),1,x(n2,i),1)
            n = n + 1
50       continue
         if (odebug(17)) then
            print*,'sred1 & ered1 after loop 50'
	    call outsqr (sred1,maxr,ntot,ntot,'sred1')
	    call outsqr (ered1,maxr,ntot,ntot,'ered1')
         end if
         m = m + 1
60    continue
      if (m1 .eq. ntot) go to 70
      m0 = m0 + maxvec
      go to 20

c   construct ered2,sred2 for library programs

70    if (odebug(17)) then
         print*,'sred1 & ered1, iteration step ',niter
	 call outsqr (sred1,maxr,ntot,ntot,'sred1')
	 call outsqr (ered1,maxr,ntot,ntot,'ered1')
      end if
      call redeig (maxr,ntot,sred1,ered1,bred,sred2,ered2)
      ntot2 = 2 * ntot
      if (odebug(17)) then
         print*,'sred2 & ered2, iteration step ',niter
	 call outsqr (sred2,maxr2,ntot2,ntot2,'sred2')
	 call outsqr (ered2,maxr2,ntot2,ntot2,'ered2')
      end if

c   solve reduced eigenvalue problem

      if (oreduc) go to 120

c   qz algorithm

      time = cpulft(1)
      call qzhes (maxr2,ntot2,ered2,sred2,.true.,xred)
      call qzit (maxr2,ntot2,ered2,sred2,0.0d0,.true.,xred,ierr)
      call qzval (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,.true.,xred)
      if (ierr .ne. 0) go to 400
      call qzvec (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,xred)
      tim(3) = tim(3) + cpulft(1) - time

c   collect real positive eigenvalues

      nrealp = 0
      do 80 i = 1 , ntot2
         tt = etare(i)
         uu = etaim(i)
         vv = sigma(i)
         if (1.0d0+uu .ne. 1.0d0 .or. 1.0d0+vv .eq. 1.0d0) go to 80
         tt = tt / vv
         if (tt .lt. 0.0d0) go to 80
         nrealp = nrealp + 1
         etare(nrealp) = tt
         if (i.ne.nrealp) call dcopy(ntot2,xred(1,i),1,xred(1,nrealp),1)
80    continue
      if (nrealp .eq. 0) go to 500
      if (nrealp .lt. nstart) then
         write (iwr,90) niter,nrealp
90       format ('warning: in iteration step ',i2,' there are only ',
     ?           i2,' positive real eigenvalues.'/
     ?           'trying to continue ...')
         call vclr(etare(nrealp+1),1,nstart-nrealp)
      end if
      nev2 = min(nrealp,nstart)

c   sort eigenvalues

      do 110 i = 1 , nev2
         k = i
         if (nrealp .gt. 1) then
	    emin = etare(i)
	    do 100 j = i+1 , nrealp
	       if (etare(j) .lt. emin) then
	          k = j
	          emin = etare(j)
               end if
100         continue
         end if
         t = etare(i)
         etare(i) = etare(k)
         etare(k) = t
         call dswap(ntot2,xred(1,i),1,xred(1,k),1)
110   continue
      go to 140

c   reduction to ordinary eigenvalue problem

120   time = cpulft(1)
      call reduc (maxr2,ntot2,sred2,ered2,etaim,ierr)
      if (ierr .ne. 0) go to 600
      call tred2(maxr2,ntot2,sred2,etare,sigma,xred)
      call tql2(maxr2,ntot2,etare,sigma,xred,ierr)
      if (ierr .ne. 0) go to 700
      call rebak (maxr2,ntot2,ered2,etaim,ntot2,xred)
      tim(3) = tim(3) + cpulft(1) - time
      do 130 i = 1 , nstart
         etare(i) = 1.0d0 / etare(ntot2+1-i)
         call dcopy(ntot2,xred(1,ntot2+1-i),1,xred(1,i),1)
130   continue
      nev2 = nstart

140   call outeig (itable,iwr,etare(nevlo1+1),neig,niter,ntot,nevlo1,1)
      if (odebug(42)) then
         write (iwr,141) niter
141      format ('iteration ',i4)
         call outvec (etare,neig,'eigenvalues')
      end if
      if (opg_root()) then
       open (ieigen,file='roots',status='unknown')
      endif
      call outeig (itable,ieigen,dummy,neig,niter,ntot,nevlo1,2)
      if (opg_root()) close (ieigen)

c   new trial vectors to be expected

      ntot0 = ntot
      ntot1 = ntot + 1
      call vfill(-2.0d0,resarr,1,neig)
      do 230 j = 1 , nev2
         if (j .lt. nevlo(isym)) then
            if (ldump.eq.0 .or. orest) then
               go to 230
            else
               go to 145
            end if
         end if
         if (iconv(j-nevlo1) .eq. 1) then
	    resarr(j-nevlo1) = -1.0d0
	    if (ldump.eq.0 .or. orest) go to 230
         end if

c   compute residue vector

145      if (orest) then
            call rdedx (v,isize2,iblbb(j),num8)
            call rdedx (e(1,1),isize2,iblee(j),num8)
         else
            call vclr(e(1,1),1,isize2)
            call vclr(v,1,isize2)
            do 170 k = 1 , ntot0
               x1kj = xred(k      ,j)
               x2kj = xred(k+ntot0,j)
               if (odebug(43)) write (iwr,1100) 170,'trial',k,iblb(k)
               call rdedx (x(1,1),isize2,iblb(k),num8)
               do 150 i1 = 1 , isize
                  i2 = i1 + isize
                  v(i1) = v(i1) + x(i1,1) * x1kj + x(i2,1) * x2kj
                  v(i2) = v(i2) + x(i2,1) * x1kj + x(i1,1) * x2kj
150            continue
               if (j .lt. nevlo(isym)) go to 170
               if (odebug(43)) write (iwr,1100) 170,'result',k,ible(k)
               call rdedx (x(1,1),isize2,ible(k),num8)
               do 160 i1 = 1 , isize
                  i2 = i1 + isize
                  e(i1,1) = e(i1,1) + x(i1,1) * x1kj + x(i2,1) * x2kj
                  e(i2,1) = e(i2,1) + x(i2,1) * x1kj + x(i1,1) * x2kj
160            continue
170         continue
            if (ldump .gt. 0) then
              if (opg_root()) then
               write (idump) (v(i),i=1,isize2)
              end if
            end if
            if (j.lt.nevlo(isym) .or. iconv(j-nevlo1).eq.1) go to 230
         end if

c   r = omega * s - e

         freq = etare(j)
         call vsmsb (v(n1),1, freq,e(n1,1),1,e(n1,1),1,isize)
         call vsmsb (v(n2),1,-freq,e(n2,1),1,e(n2,1),1,isize)
         if (odebug(47)) then
            write (iwr,61) j
61          format ('residue vector no.',i3,':')
            call outvec (e(1,1),isize2,' ')
         end if
         t = dnrm2(isize2,e(1,1),1)
         resarr(j-nevlo1) = t

c   threshold ?

         if (t .le. eps) then
	    iconv(j-nevlo1) = 1
	    call wrt3 (v,isize2,ibleig(j-nevlo1),num8)
	    go to 230
         end if
         if (ntot.ge.maxr .or. ntot.ge.isize) go to 230
         ntot = ntot + 1

c   update residue vector 

         do 190 ii = 1 , ncore
            do 180 ia = nc1 , nba
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
	          k1 = ipo(ii,ia)
	          k2 = k1 + isize
	          diff = evalue(ia) - evalue(ii)
	          denom1 =  freq - diff
	          denom2 = -freq - diff
	          if (dabs(denom1) .gt. small) then
	             e(k1,1) = e(k1,1) / denom1
                  else
	             e(k1,1) = 1.0d0
                  end if
	          if (dabs(denom2) .gt. small) then
	             e(k2,1) = e(k2,1) / denom2
                  else
	             e(k2,1) = 1.0d0
                  end if
               end if
180         continue
190      continue

c   gram-schmidt orthogonalisation

         do 200 k = 1 , ntot-1
            if (odebug(43)) write (iwr,1100) 200,'trial',k,iblb(k)
            call rdedx (x(1,1),isize2,iblb(k),num8)
            t = - ddot(isize,e(n1,1),1,x(n1,1),1) 
     ?	        - ddot(isize,e(n2,1),1,x(n2,1),1)
            call daxpy(isize,t,x(n1,1),1,e(n1,1),1)
            call daxpy(isize,t,x(n2,1),1,e(n2,1),1)
            t = - ddot(isize,e(n1,1),1,x(n2,1),1) 
     ?	        - ddot(isize,e(n2,1),1,x(n1,1),1)
            call daxpy(isize,t,x(n2,1),1,e(n1,1),1)
            call daxpy(isize,t,x(n1,1),1,e(n2,1),1)
200      continue
         t = dnrm2(isize2,e(1,1),1)
         if (1.0d0+t .eq. 1.0d0) then
            write (iwr,210) j,niter
210         format ('discard new trial vector no. ',i2,' in iteration',
     ?              ' step ',i2/'because of linear dependence from',
     ?              ' the previous trial vectors.')
            ntot = ntot - 1
         else
            call symort (isize,e(n1,1),e(n2,1),sytol,oldep)
            if (oldep) then
               write (iwr,220) j,niter
220            format ('discard new trial vector no. ',i2,' in',
     ?                 ' iteration step ',i2/'because b and b#',
     ?                 ' are (nearly) linearly dependent.')
               ntot = ntot - 1
	    else
	       call wrt3 (e(1,1),isize2,iblb(ntot),num8)
            end if
         end if
230   continue
      if (ldump .gt. 0) then
        if (opg_root()) then
         close (idump)
        end if
      end if
      call outres (iwr,inorms,resarr,nevlo1,neig,niter,1)
      if (opg_root()) then
       open (iresid,file='residues',status='unknown')
      endif
      call outres (iresid,inorms,dummy,nevlo1,neig,niter,2)
      if (opg_root()) close (iresid)
      nconv = isum(neig,iconv,1)
      if (nconv .eq. neig) go to 1000

c   if at least one new trial vector has been added,
c   start next iteration step

      if (ntot .ge. ntot1) go to 10
      if (ntot .ge. maxr) then
         do 260 j = 1 , nev2
            call vclr(v,1,isize2)
            do 250 k = 1 , ntot0
               x1kj = xred(k      ,j)
               x2kj = xred(k+ntot0,j)
               if (odebug(43)) write (iwr,1100) 250,'trial',k,iblb(k)
               call rdedx (x(1,1),isize2,iblb(k),num8)
               do 240 i1 = 1 , isize
                  i2 = i1 + isize
                  v(i1) = v(i1) + x(i1,1) * x1kj + x(i2,1) * x2kj
                  v(i2) = v(i2) + x(i2,1) * x1kj + x(i1,1) * x2kj
240            continue
250         continue
            if (odebug(43)) write (iwr,117) j,iblbb(j)
117         format (
     +            'diralg, loop 260: writing restart trial vector ',i2,
     ?            ' to block ',i12)
            call wrt3 (v,isize2,iblbb(j),num8)
            if (odebug(44)) 
     +      call outvec (v,isize2,'restart trial vector')
            call vclr(v,1,isize2)
            do 255 k = 1 , ntot0
               x1kj = xred(k      ,j)
               x2kj = xred(k+ntot0,j)
               if (odebug(43)) write (iwr,1100) 255,'product',k,ible(k)
               call rdedx (e(1,1),isize2,ible(k),num8)
               do 245 i1 = 1 , isize
                  i2 = i1 + isize
                  v(i1) = v(i1) + e(i1,1) * x1kj + e(i2,1) * x2kj
                  v(i2) = v(i2) + e(i2,1) * x1kj + e(i1,1) * x2kj
245            continue
255         continue
            if (odebug(43)) write (iwr,119) j,iblee(j)
119         format (
     +       'diralg, loop 260: writing restart product vector ',i2,
     ?       ' to block ',i12)
            call wrt3 (v,isize2,iblee(j),num8)
            if (odebug(44)) 
     +      call outvec (v,isize2,'restart product vector')
260      continue
         call end_time_period(TP_RPA)
         return
      end if
      write (iwr,270) niter
270   format('no new trial vector has been added in iteration step ',i2,
     ? ':'/'unable to finish iterative calculation of rpa eigenvalues.')
      ierr = 1
      go to 1000

c   error messages

400   write (iwr,410) niter,ierr
410   format (/'convergence problems with QZIT in iteration step',
     ?        i3,': ierr = ',i3/)
      go to 1000

500   write (iwr,510) niter
510   format (/
     + 'in iteration step ',i3,' there are no real eigenvalues.'/)
      ierr = 1
      go to 1000

600   write (iwr,610) niter,ierr
610   format (/
     +' problem in iteration step ',i2,': matrix ERED2 is not positive',
     ?	       ' definite.'/'ierr = ',i4)
      go to 1000

700   write (iwr,710) niter,ierr
710   format (/'convergence problems with TQL2 in iteration step',
     ?        i2,': ierr = ',i2/)

1000  niter = min(niter,maxit)

      call end_time_period(TP_RPA)
      return
1100  format ('loop ',i5,': trying to read ',a8,' vector ',i2,
     ?        ' from block ',i12)
      end
      subroutine rpbody (iscat2,zmem,zmemc,xijk,xija,xijt,
     ?		 xijp1,xijp2,xijp3,nij,xklk,xkla,xklt,
     ?           xklp1,xklp2,xklp3,igath,nkl,d1,d2,u1,u2,utemp,codens,
     ?           ext1,ext2,ext3,i1,i2,j1,j2,k1,l1,l2,prefac,xinner)

c---------------------------------------------------------------------
c   Computes the two-electron integrals. Virtually identical to
c   GAMESS routine body.
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      logical ext1,ext2,ext3,ext4,ext5,ext7

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
c     A common block to pass the actual time periods around in the
c     response theory code. This is needed because in subroutine 
c     such as rpmuld and rpbody the time periods that are used 
c     should be different for RPA and TDA calculations. Hence the
c     usual parameters cannot be used directly.
c
      integer TP_RESP_MO2AO
      integer TP_RESP_INT
      integer TP_RESP_CNTRCT
      integer TP_RESP_AO2MO
      common/response_timep/TP_RESP_MO2AO,TP_RESP_INT,TP_RESP_CNTRCT,
     +                      TP_RESP_AO2MO
c

      common/cider /intgrl(8,100),intcon(5,12),ihz(5,5),
     ?              ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri,
     ?              mempto,memcto,mtemp,ipszt,icszt,
     ?              ipsize,icsize,memcon
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),
     +              nc(maxorb),isatbf(maxat),isatp(maxat),
     +              ipatbf(maxat),ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,n1,npp,n2,n3,np3,
     ?              npp3,nprm,n4,nm1,nm2,ne,repen,n5,n6,n7,n8,
     ?              n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd
      common/gen   /maxit,iter,concrit,conv,energy,thresh,dmax
      common/help  /gg(65536)
      common/indexd/ispdf1(21,2),nvl(6),ixyz(3,21,5),iarray(6,6)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/save  /info(5,nbrkmx),tmax(nbrkmx),nbreak
      common/scra  /iso(maxorb,48)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/tempit/zitest,zipass,zjtest,zjpass

      dimension zmem(ipsize,mempto),zmemc(icsize,memcto),
     ?     xija(nij),xijk(nij),xijp1(nij),xijp2(nij),xijp3(nij),
     ?     xkla(nkl),xklk(nkl),xklp1(nkl),xklp2(nkl),xklp3(nkl),
     ?     ab(3),naind(3),ncind(3),pa(3),igath(nkl),xijt(nij),xklt(nkl),
     ?     m0(48),m1(48),m2(48),m3(48)
      dimension d1(*),d2(*),u1(*),u2(*),utemp(*),codens(nb,nb),
     ?		iscat2(nkl)

      if (odebug(37)) print*,'entering rpbody'
      if (odebug(37)) print*,'dmax = ',dmax
      if (odebug(37)) print*,'nbreak = ',nbreak
      if (odebug(37)) print*,'nb = ',nb
      if (odebug(37)) call outvec (tmax,nbrkmx,'tmax')
      if (odebug(37)) print*,'ipsize = ',ipsize
      if (odebug(37)) print*,'icsize = ',icsize
      if (odebug(37)) print*,'mempto = ',mempto
      if (odebug(37)) print*,'memcto = ',memcto

      call start_time_period(TP_RESP_INT)

      ivlx = nvl(ix)
      jvlx = nvl(jx)
      kvlx = nvl(kx)
      lvlx = nvl(lx)
      ipos = nposf(i1,ivlx,nd,np,ns)
      jpos = nposf(j1,jvlx,nd,np,ns)
      kpos = nposf(k1,kvlx,nd,np,ns)
      lpos = nposf(l1,lvlx,nd,np,ns)
      if (odebug(37)) then
	 print*,'ivlx,jvlx,kvlx,lvlx = ',ivlx,jvlx,kvlx,lvlx
	 print*,'ipos,jpos,kpos,lpos = ',ipos,jpos,kpos,lpos
      end if
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
            dij = 2.0d0 * (codens(ic,jc) + codens(jc,ic))
            ab(1) = cobas(ic,1) - cobas(jc,1)
            ab(2) = cobas(ic,2) - cobas(jc,2)
            ab(3) = cobas(ic,3) - cobas(jc,3)
            rab = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
c
c           jpmax = nc(jc)
c
c****** now the loop over the number of different kl
c****** mini-batches in this particular maxi-batch
c
c            if (odebug(83)) print*,'ic,jc = ',ic,jc
            do 440 nbrk = 1 , nbreak
               if (odebug(84)) 
     +                print*,'ix,jx,kx,lx,ic,jc,nbrk,itprim = ',
     ?                        ix,jx,kx,lx,ic,jc,nbrk,itprim
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
	       if (odebug(37)) print*,'test1, tout = ',test1,tout
               if (test1.lt.tout) then
                  test2 = xinner/xijtm
                  if (odebug(37)) print*,'test2 = ',test2
                  tmax1 = tmax(nbrk)*10.0d0
                  jstart = info(5,nbrk)
c********* call to integral tests
                  if (nt.eq.1) then
		     if (odebug(37)) print*,'calling test11'
                     call test11(codens,dij,dmaxij,tmax1,test2,
     +                      xklt(jstart+1),igath,zmem(1,4),ic,jc,
     +                      iscat2(jstart+1),kmin,kmax,l1,l2,
     +                      lmax,lmin,nkl,ipsize,ext2,
     +                      jcount,jc0,jc1)
                  else
		     if (odebug(37)) print*,'calling test1s'
                     call test1s(codens,dij,dmaxij,test2,
     +                      xklt(jstart+1),igath,zmem(1,4),ic,jc,
     +                      iscat2(jstart+1),kmin,kmax,m1,m2,
     +                      m3,l1,l2,lmax,lmin,nkl,ipsize,
     +                      ext2,ext3,jcount,jc0,jc1)
                  end if
		  if (odebug(37)) print*,'jc0 = ',jc0
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
c...carsten0
		  if (odebug(79)) then
                   write(iwr,9000)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9000               format(3i3,2x,'after test1s : zmem(1,.) = ',10f12.7)
                  end if
c...carsten1
	if (odebug(38)) then
	   call outvec (zmem(1,1),jc0,'zmem(1,1) in rpabod')
	   call outvec (zmem(1,2),jc0,'zmem(1,2) in rpabod')
	   call outvec (zmem(1,3),jc0,'zmem(1,3) in rpabod')
        endif
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
c...carsten0
		  if (odebug(79)) then
                   write(iwr,9001)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9001               format(3i3,2x,'after loop 80: zmem(1,.) = ',10f12.7)
                  end if
c...carsten1
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
c...carsten0
		  if (odebug(79)) then
                   write(iwr,9002)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9002               format(3i3,2x,'after sinter : zmem(1,.) = ',10f12.7)
                  end if
c...carsten1
c
                        if (ext7) zmem(jc0,9) = zmem(jc0,9)*2.0d0
c
c*** loop over number of distinct integral classes
c*** needed to form the target integrals
c
c       if (odebug(83)) 
c    + call outsqr(zmem,ipsize,ipsize,mempto,'zmem before 220')
        if (odebug(83)) 
     + write (iwr,219) ix,jx,kx,lx,ic,jc,nbrk,itprim,(zmem(1,iz),iz=1,8)
219          format(8i3,8f9.4)
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
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9010)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9010                format(3i3,2x,'after recur1: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
                                 go to 200
 130                             imem3 = jmem3 + k
                                 imem4 = jmem4 + k
                                 call recur2(zmem(1,memcur),pa(index),
     +                              zmem(1,imem1),zmem(1,indwp),
     +                              zmem(1,imem2),xija(icount),nai,
     +                              zmem(1,imem3),zmem(1,7),zmem(1,8),
     +                              zmem(1,imem4),jc0)
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9020)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9020                format(3i3,2x,'after recur2: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur2')
c...carsten1
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
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9030)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9030                format(3i3,2x,'after recur3: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur3')
c...carsten1
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
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9040)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9040                format(3i3,2x,'after recur4: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur4')
c...carsten1
                                 go to 200
 160                             call recur5(zmem(1,memcur),
     +                              zmem(1,indqc),zmem(1,imem1),
     +                              zmem(1,indwq),zmem(1,imem2),jc0)
c...carsten0
		   if (odebug(79)) then
                      write(iwr,9050)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9050                 format(3i3,2x,'after recur5: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur5')
c...carsten1
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
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9060)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9060                format(3i3,2x,'after recur6: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur6')
c...carsten1
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
c...carsten0
		   if (odebug(79)) then
                   write(iwr,9070)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9070               format(3i3,2x,'after recur7: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur7')
c...carsten1
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
c...carsten0
		   if (odebug(79)) then
                    write(iwr,9080)ic,jc,nbrk,(zmem(1,iz),iz=1,10)
9080                format(3i3,2x,'after recur8: zmem(1,.) = ',10f12.7)
                   end if
c...carsten1
c...carsten0
c		if (odebug(81)) call outvec(zmem(1,memcur),jc0,'recur8')
c...carsten1
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
                        if (ncr .eq. 1) then
                        if (ext4) then
                        call vclr(zmem(jc0,imem2m)
     +                  ,1,jc0+itprim)
                        end if
                        call dcopy(jc0
     +                  ,zmem(1,imem1m),1
     +                  ,zmem(1,imem2m),1)
                        else
                         call vadd(zmem(1,imem2m),1
     +                            ,zmem(1,imem1m),1
     +                            ,zmem(1,imem2m),1
     +                            ,jc0)
                        if (odebug(80)) 
     +                  call outvec (zmem(1,imem2m),jc0,
     +                  'zmem(1,imem2m) in loop 230')
                        end if
 230                          continue
 240                       continue
 250                    continue
c
 260                 continue
        if (odebug(82)) 
     +      call outsqr(zmem,ipsize,ipsize,mempto,'zmem after 260')
 265                 continue
c
c****** now that the loops over the i and j primitives for the one
c****** particular ic,jc contracted pairing have been closed, we
c****** can contract the kl part of the integrals.
c
                     if (jc0.gt.0.and.ncr.ne.0) then
          call igthr(jc0,iscat2(jstart),igath(1),igath(1))
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
****** a d function is being added to the i/a position
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
c***** send the integrals to the u build!
                        call end_time_period(TP_RESP_INT)
                        call start_time_period(TP_RESP_CNTRCT)
                        if (onew) then
                           kdim = (kmax-kmin+1)*kvlx
                           ldim = (l2-l1+1)*lvlx
                           call ubuil2(zmemc(1,memans+1),ivlx,jvlx,kvlx,
     +                          lvlx,ic,jc,kmin,kmax,lmin,lmax,l1,l2,
     +                          ext2,u1,u2,utemp,d1,d2,gg,memreq,
     +                          ipos,jpos,kpos,lpos,i1,j1,k1,kdim,ldim)
                        else
                           call ubuild(zmemc(1,memans+1),ivlx,jvlx,kvlx,
     +                          lvlx,ic,jc,kmin,kmax,lmin,lmax,l1,l2,
     +                          ext2,u1,d1,d2,zmemc(1,mtemp),memreq,
     +                          ipos,jpos,kpos,lpos,i1,j1,k1)
                        end if
                        call end_time_period(TP_RESP_CNTRCT)
                        call start_time_period(TP_RESP_INT)
                     end if
                  end if
               end if
c
 440        continue
 450        icorig = icorig + itprim
            zjtest = zjtest + itprim
 460     continue
 470  continue
      call end_time_period(TP_RESP_INT)
      return
      end
      block data dirrpa_data
      implicit real*8 (a-h,p-y), integer (i-n), logical (o)
      common/rpacom/modus,nevlo(8),nevhi(8),nev(8),
     +              maxr,maxit,maxrr,maxitt,nbegin,
     +              eps,epss,epsana,epstab,
     +              sytol,tolmy,fetda,ferpa,tottda,totrpa,ipref,iblo13
      data tolmy/1.0d-9/
      end
      subroutine rpdriv (q,ipo)
c
c----------------------------------------------------------------------
c   Driving routine for RPA
c   (c) Carsten Fuchs 1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-y), integer (i-n), logical (o)
      implicit character*8 (z)
      character *8 obnam
      character*80 dmpfil,rstfil
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
      common/blkin /potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/block /iky(maxorb+1),ilifq(maxorb+1)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypl/obnam(29)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall,ioneel,jasym
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      common/gen2  /sipftd(21,2),sijkld(21,2),sjpftd(21,2),
     +              sipfal(21,2),sijkal(21,2),
     +              prefac(21),xinner(21),sifall,siffew,sicfb1,sicfb2
      common/gtrans/idtra(maxorb)
      common/infoa /nat,ich,mul,nba,nx,ne,ncore,nee,czan(maxat),
     +              c(3,maxat),amass(maxat),iannn(2,maxat),
     +              ipseud(maxat),symz(maxat),lpseud
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr
      common/junkb /csi(mxprim),cpi(mxprim),cdi(mxprim),cfi(mxprim),
     +              cznuc(maxat)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/molsym/tr(3,3),trx,try,trz,index,jaxis,igroup,ispmol,
     +              prmoms(3),aprev(maxat3,3)
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
      common/phycon/toang,phyc(83),iunt,nosym
c
      real*8 freq, pola
      integer nfreq, ipola, mxfreq
      parameter (mxfreq=100)
      common /polars/ nfreq,ipola(3,3),freq(mxfreq),pola(6,mxfreq)
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/restar/nprint,itol,icut,normf,normp,nopk,iresta,
     +              nrec,intloc,
     +		    ist,jst,kst,lst,jspace(458),icoord
      common/rpacom/modus,nevlo(8),nevhi(8),nev(8),
     +              maxr,maxit,maxrr,maxitt,nbegin,
     +              eps,epss,epsana,epstab,
     +              sytol,tolmy,fetda,ferpa,tottda,totrpa,ipref,iblo13
c
      integer itdspc, irpspc, itdtex, irptex, itable, inorms, ivect
      integer itmfil, ldump, idump, lrest, irest, mrest
      common /rpfile/ itdspc,irpspc,itdtex,irptex,itable,inorms,ivect,
     +                itmfil,ldump,idump,lrest,irest,mrest
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
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/table /dmpfil,rstfil,rlines(2)
      common/tran  /ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +              ctran(mxorb3),otran,ift(2)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
c
      integer ltask, lcoord, iplot, iprop, imull, idma, inbo
      logical opass1, opass2, opass3, orege, opass4, offsym
      logical opass5, opass6, opass8, opass9, odisc, oirc
      logical orgall, orestj, opass10, opass11
      logical omrdci, ofulci, omech, omopac
      logical orpa, omclr, odirpa, occsd, occsdt, oqcisd, oqcisdt
      logical omrdcid,opark,oclunp,ofill,odiesel
c
      common /restrj/ ltask,lcoord,opass1,opass2,opass3,orege,iplot,
     + iprop,imull,opass4,offsym,opass5,opass6,opass10,idma,inbo,
     + omrdci,ofulci,omech,omopac,opass8,opass9,odisc,
     + oirc,orpa,omclr,odirpa,occsd,occsdt,oqcisd,oqcisdt,orgall,
     + omrdcid,orestj(7),opark,oclunp,ofill(6),odiesel, opass11
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
      common/atmol3/minv(2),mouta,moutb
      common/fsave/ofiles(2)
      character*10 charwall
      character *6 snm
      character *8 fnm
      data fnm/'dirrpa.m'/
      data snm/'rpdriv'/
c
      dimension q(*),ipo(*)
c
      data m51/51/
c
      time = cpulft(1)
      call cpuwal(begin,ebegin)
      odirct = odirpa
      call start_time_period(TP_RESPONSE)
c
c   header
c
      write (iwr,520)
520   format (///19x,83('*')/19x,'*',81x,'*')
      iform = modus
      if (odirct) iform = iform + 3
      if (iform .eq. 1) write (iwr,521)
      if (iform .eq. 2) write (iwr,522)
      if (iform .eq. 3) write (iwr,523)
      if (iform .eq. 4) write (iwr,524)
      if (iform .eq. 5) write (iwr,525)
      if (iform .eq. 6) write (iwr,526)
521   format (19x,'*    TDA CALCULATION             ',49x,'*')
522   format (19x,'*    TDA  & RPA CALCULATION      ',49x,'*')
523   format (19x,'*    RPA CALCULATION             ',49x,'*')
524   format (19x,'*    DIRECT TDA CALCULATION      ',49x,'*')
525   format (19x,'*    DIRECT TDA & RPA CALCULATION',49x,'*')
526   format (19x,'*    DIRECT RPA CALCULATION      ',49x,'*')
      write (iwr,540)
540   format (19x,'*',81x,'*'/19x,83('*')//122('=')//)
c
c   direct RPA: check point group
c
      if (odirct) then
         if (igroup.eq.1                   .or.                    ! C1
     ?       igroup.eq.2                   .or.                    ! Cs
     ?      (igroup.eq.6 .and. jaxis.eq.2) .or.                    ! C2h
     ?      (igroup.eq.7 .and. jaxis.eq.2) .or.                    ! C2v
     ?      (igroup.eq.8 .and. jaxis.eq.2) .or.                    ! D2
     ?      (igroup.eq.9 .and. jaxis.eq.2)     ) go to 3           ! D2h
         write (iwr,2000)
2000     format(//'---> You cannot use this point group in the DIRECT',
     +            ' RPA program.    <---'
     +           /'---> Please use one of the following groups: C1,Cs',
     +            'C2h,C2v,D2,D2h. <---'//)
         call caserr ('invalid group for direct rpa program')
3        continue 
c
c   prefactor tolerances for integrals
c
         if (ipref .eq. 0) then
            do 5 i = 1 , 21
               prefac(i) = 0.0d0
               xinner(i) = 0.0d0
5           continue
         else
            do 6 i = 1 , 21
               prefac(i) = 10.0d0 ** (-ipref)
               xinner(i) = 10.0d0 ** (-ipref)
6           continue
         end if
      end if
c
c     determine maximum available memory
c
      nmaxly  = igmem_max_memory()
c
      nat = non
      nba = numorb
      nbatri = (nba * (nba+1)) / 2
      nbasq = nba * nba
      nx = nbatri
c
c   get symmetry information
c
      call secget(isect(490),m51,iblks)
      nav = lenwrd()
      call readi(mmmm,mach(13)*nav,iblks,idaf)
c
c   read in orbital energies and information about symmetry adaption
c
      call secget (mouta,3,iblko)
      call rdchr (obnam,29,iblko,idaf)
      call reads (evalue,mach(8),idaf)
      call readis (ilifc,mach(9)*nav,idaf)
      a = 0.5d0 * dsum(nba,hocc,1) + 0.1d0
      ncore = a
      nc1 = ncore + 1
      nvirt = nba - ncore
      nvc = nvirt * ncore
      nvc2 = nvc + nvc
c
      lenipo = lenint(nbasq)
      istart = igmem_alloc_inf(lenipo,fnm,snm,'ipo',IGMEM_DEBUG)
      iipo = (istart-1)*nav
      idip = igmem_alloc_inf(3*nvc,fnm,snm,'idip',IGMEM_DEBUG)
      if(ofiles(1)) then
c       allocate extra memory for subroutines rpanal and outtdm
c       if RPA Transition Density Matrices are output.
        itdm = igmem_alloc_inf(3*nbasq+2*nbatri,fnm,snm,'itdm',
     &                         IGMEM_DEBUG)
      else
        itdm = idip
      endif
      icoef = igmem_alloc_inf(nbasq,fnm,snm,'icoef',IGMEM_DEBUG)


c
c   read in scf mos
c
      iblkq = iblko + lensec(mach(8)) + lensec(mach(9)) + 1
      call rdedx (q(icoef),nbasq,iblkq,idaf)
      if (odebug(5))
     +   call outsqr (q(icoef),nba,nba,nba,'original mo coefficients')
c
      do 10 i = 1, nba+1
         ilifq(i) = (i-1)*nba
10    continue
c
c   determine mo symmetries
c
      do 15 i = 1 , nba
         ii = icoef + ilifq(i)
         ibig = idamax(nba,q(ii),1)
         isymmo(i) = isymao(ibig)
15    continue
c
c   now isymao is not needed any longer, so isymao(i) will be the
c   number of mo i in its irrep
c
      call setsto(8,0,npair)
      do 17 i = nbegin , nba
         k = isymmo(i)
         l = npair(k) + 1
         isymao(i) = l
         npair(k) = l
17    continue
      call setsto(8,0,npair)
c
c   number of irreducible representations
c
      nir = imaxf(nba,isymmo,1)
      if (nir .gt. 0) nirr = nir
      if (nir .gt. 2) nirr = 4
      if (nir .gt. 4) nirr = 8
      nirs = nirr - 2
      if (nirr .eq. 1) nirs = 0
c
c   lookup arrays
c
      k = 0
      do 50 isym = 1 , nirr
	 icnt = 0
	 do 40 i = 1 , nirr
	    j = iperm(isym,i)
	    iofset(isym,i) = icnt + 1
	    do 30 ii = 1 , ncore
	       if (isymmo(ii) .eq. j) then
	          do 20 ia = nc1 , nba
	             if (isymmo(ia) .eq. i) then
		        icnt = icnt + 1
                        ipo(ii+ilifq(ia)+iipo) = icnt
                        ipo(ia+ilifq(ii)+iipo) = icnt
                     end if
20                continue
	       end if
30          continue
40       continue
         icount(isym) = icnt
         iofsym(isym) = k
         k = k + icnt
50    continue

c
c   number of eigenvalues wanted
c
      do 60 i = 1 , nirr
         nevhi(i) = min(icount(i),nevhi(i))
         nev(i) = nevhi(i) - nevlo(i) + 1
         if (nev(i).lt.0 .or. nevhi(i).eq.0) nev(i) = 0
60    continue
c
      if (.not. odirct) go to 200
c
c   memory handling
c
c       isspac(1) ... mo coefficients & dipole integrals & one (nb,nb)
c                     array for compressed density matrix
c       isspac(2) ... iterative rpa procedure (without b,e)
c       isspac(3) ... memory for integrals
c       isspac(4) ... memory for ipar arrays (part of isspac(3))
c       isspac(5) ... d,dd & u matrices and b,e from diralg
c       isspac(6) ... isspac(1) + isspac(2)
c       isspac(7) ... isspac(1) + isspac(2) + isspac(3)

      isspac(1) = 2*nbasq + 3*nvc + lenipo
      maxr2 = maxr + maxr
      maxrs = maxr * maxr
      maxr2s = 4 * maxrs
      isspac(2) = (4*maxr)/nav + 2*maxrs + maxr + 3*maxr2s + 3*maxr2
      k = 0
      do 70 i = 1 , nirr
	 kk = 0
	 if (nev(i) .gt. 0) then
	    kk = (2*nev(i))/nav + 2*icount(i) + nev(i)
         end if
         k = max(k,kk)
70    continue
      isspac(2) = isspac(2) + k
      write (iwr,601) nmaxly,isspac(1),isspac(2)
601   format (/72('*')/'*',70x,'*'/
     ?'*  total space available:                         ',i14,7x,'*'/
     ?'*  space for mo coefficients and dipole integrals:',i14,7x,'*'/
     ?'*  space for iterative rpa procedure:             ',i14,7x,'*')
      isspac(6) = isspac(1) + isspac(2)
      if (isspac(6) .gt. nmaxly) then
	 write (iwr,602)
602	 format ('too much space needed for direct rpa -- ',
     ?           'increase memory or decrease maxr !')
         call caserr ('error in direct rpa program')
      end if
      write (iwr,6025) isspac(3)
6025  format (
     +'*  space for two-electron integrals etc.:         ',i14,7x,'*')
      isspac(7) = isspac(6) + isspac(3)
      if (isspac(7) .gt. nmaxly) then
	 write (iwr,603)
603	 format
     +     ('too much space needed for integrals in direct rpa -- '/
     ?      'increase memory or decrease maxr or decrease isspac(3) !')
         call caserr ('error in direct rpa program')
      end if
      isspac(5) = nmaxly - isspac(7)
      k = 0
      do 80 i = 1 , nirr
	 if (nev(i) .gt. 0) k = max(k,icount(i))
80    continue
      if (onew) then
         maxvec = isspac(5) / (4*nbasq + 4*k + 1)
      else
         maxvec = isspac(5) / (3*nbasq + 4*k)
      end if
      write (iwr,6035) isspac(5),maxvec
6035  format (
     + '*  space for trial vectors & matrices d,dd,u:     ',i14,7x,'*'/
     + '*  maximal no of trial vectors in one batch:      ',i14,7x,'*')
      if (maxvec .lt. 1) then
	 write (iwr,604)
604	 format ('not enough space for matrices d,dd,u associated to ',
     +      'trial vectors -- '/
     +      'increase memory or decrease maxr or decrease isspac(3) !')
         call caserr ('error in direct rpa program')
      end if
      if (mymaxv .ne. 0) maxvec = mymaxv
      write (iwr,605)
605   format ('*',70x,'*'/72('*')/)

c   transform scf mos to non-symmetry-adapted ao basis

200   otran = .false.
      i1 = igmem_alloc_inf(nbasq,fnm,snm,'tdown',IGMEM_DEBUG)
      call tdown (q(i1),ilifq,q(icoef),ilifq,nba)
      call dcopy(nbasq,q(i1),1,q(icoef),1)
      if (odebug(5))
     ?	 call outsqr (q(icoef),nba,nba,nba,
     +                'original mo coefficients, non-symmetry-adapted')
      call gmem_free_inf(i1,fnm,snm,'tdown')
c
c   read in dipole integrals
c
      ix = igmem_alloc_inf(3*nbatri,fnm,snm,'rpdipo1',IGMEM_DEBUG)
      i1 = igmem_alloc_inf(nbasq,fnm,snm,'rpdipo2',IGMEM_DEBUG)
      i2 = igmem_alloc_inf(nbasq,fnm,snm,'rpdipo3',IGMEM_DEBUG)
      call rpdipo (ipo(iipo+1),nba,q(ix),q(icoef),q(i1),q(i2),q(idip))
      call gmem_free_inf(i2,fnm,snm,'rpdipo3')
      call gmem_free_inf(i1,fnm,snm,'rpdipo2')
      call gmem_free_inf(ix,fnm,snm,'rpdipo1')
c
      if (.not. odirct) go to 365
      m48 = 48 * maxorb
      i1 = igmem_alloc_inf(m48,fnm,snm,'dircon1',IGMEM_DEBUG)
      i2 = igmem_alloc_inf(m48,fnm,snm,'dircon2',IGMEM_DEBUG)
      i3 = igmem_alloc_inf(m48,fnm,snm,'dircon3',IGMEM_DEBUG)
      i4 = igmem_alloc_inf(maxat,fnm,snm,'dircon4',IGMEM_DEBUG)
      call dircon (q(i1),q(i2),q(i3),q(i4))
      call gmem_free_inf(i4,fnm,snm,'dircon4')
      call gmem_free_inf(i3,fnm,snm,'dircon3')
      call gmem_free_inf(i2,fnm,snm,'dircon2')
      call gmem_free_inf(i1,fnm,snm,'dircon1')
c
c   permute the columns of c so that the mo's are grouped 
c   according to irreps
c
      i1 = igmem_alloc_inf(nbasq,fnm,snm,'permute',IGMEM_DEBUG)
      k1 = i1
      do 350 i = 1 , nirr
         k0 = icoef
         do 340 j = 1 , nba
            if (isymmo(j) .eq. i) then
               call dcopy(nba,q(k0),1,q(k1),1)
               k1 = k1 + nba
            end if
            k0 = k0 + nba
340      continue
350   continue
      if (odebug(6)) 
     +	call outsqr (q(i1),nba,nba,nba,'semipermuted mo coefficients')
c
c   permute the rows of c in order to convert to direct format
c
      iold = i1 - 1
      inew = icoef - 1
      do 360 i = 1 , nba
         call dcopy(nba,q(iold+i),nba,q(inew+idtra(i)),nba)
360   continue
c
      if (odebug(6)) 
     ?	 call outsqr (q(icoef),nba,nba,nba,'permuted mo coefficients')
c
      call gmem_free_inf(i1,fnm,snm,'permute')
c
      i1 = igmem_alloc_inf(nbasq*maxvec,fnm,snm,'d1',IGMEM_DEBUG)
      i2 = igmem_alloc_inf(nbasq*maxvec,fnm,snm,'d2',IGMEM_DEBUG)
      i3 = igmem_alloc_inf(nbasq*maxvec,fnm,snm,'u1',IGMEM_DEBUG)
c
      if (onew) then
         i4 = igmem_alloc_inf(nbasq*maxvec,fnm,snm,'u2',IGMEM_DEBUG)
         i5 = igmem_alloc_inf(maxvec,fnm,snm,'utemp',IGMEM_DEBUG)
      else
         i4 = igmem_alloc_inf(1,fnm,snm,'u2',IGMEM_DEBUG)
         i5 = i4
      end if
      i6 = igmem_alloc_inf(nbasq,fnm,snm,'scr',IGMEM_DEBUG)
c
c   offsets
c
365   continue
      do 370 i = 1 , nirr
	 kcore(i) = 0
	 kvirt(i) = 0
370   continue
      do 380 i = 1 , ncore
         k = isymmo(i)
         kcore(k) = kcore(k) + 1
380   continue
      do 390 i = nc1 , nba
         k = isymmo(i)
         kvirt(k) = kvirt(k) + 1
390   continue
      do 400 i = 1 , nirr
	 ksymm(i) = kcore(i) + kvirt(i)
400   continue
      n = 1
      do 410 i = 1 , nirr
         lcore(i) = n
         n = n + kcore(i)
	 lvirt(i) = n
	 n = n + kvirt(i)
410   continue
c
c   print out information
c
      write (iwr,550) nba,ncore,maxr,maxit,eps,epsana
550   format ('number of basis functions:           ',7x,i4/
     ?        'number of core orbitals:             ',7x,i4/
     ?        'maximal size of reduced matrices:    ',7x,i4/
     ?        'maximum number of iterations:        ',7x,i4/
     ?        'residue tolerance:                   ',5x,f6.4/
     ?        'analysis threshold:                  ',5x,f6.4//)
      write (iwr,710) (i,i=1,nirr)
710   format ('core orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
720   format (56a1)
      write (iwr,730) (kcore(i),i=1,nirr)
730   format ('kcore  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,740) (i,i=1,nirr)
740   format (//'virtual orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,750) (kvirt(i),i=1,nirr)
750   format ('kvirt  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,760) (i,i=1,nirr)
760   format (//'all orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,770) (ksymm(i),i=1,nirr)
770   format ('ksymm  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,780) (i,i=1,nirr)
780   format (//'dimension of the matrices A,B:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,790) (icount(i),i=1,nirr)
790   format ('icount |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,800)
800   format (//'the following roots are to be computed:'//
     ?        'irrep        roots')
      write (iwr,810)
810   format ('---------------------')
      do 830 i = 1 , nirr
         if (nev(i) .gt. 0) write (iwr,820) i,nevlo(i),nevhi(i)
820      format (i3,8x,i3,' to ',i3)
830   continue
      write (iwr,810)

      write (iwr,560)
560   format (//23x,18('=')/23x,'Specified options:'/23x,18('=')/)
      if (oanal) write (iwr,570)
570   format (23x,
     +        'analyse          ---  detailed analysis of excitations')
      if (nbegin .gt. 1) write (iwr,580) nbegin,nbegin
580   format (23x,'begin',i4,8x,
     ?            '---  symmetry counting of orbitals starts with',
     ?            ' orbital',i4)
      if (modus .ge. 2) then
         if (oreduc) then
            write (iwr,590)
590         format (23x,'diagonalization  ---  reduce')
         else
            write (iwr,600)
600         format (23x,'diagonalization  ---  qz')
         end if
      end if
      if (lrest .gt. 0) write (iwr,610) (rstfil(i:i),i=1,lrest)
610   format (23x,'restore          ---  use initial vectors from file'/
     ?        45x,80a1)
      if (ldump .gt. 0) write (iwr,620) (dmpfil(i:i),i=1,ldump)
620   format (23x,
     +     'dump             ---  dump temporary eigenvectors to file'/
     ?        45x,80a1)
      if (modus .eq. 2) write (iwr,630)
630   format (23x,'tda              ---  do also TDA calculation')
      if (modus .eq. 1) write (iwr,640)
640   format (23x,'tda only         ---  do only TDA calculation')
      if (ouse) write (iwr,58)
58    format (23x,'use tda vectors  ---  start RPA calculation',
     ?            ' with TDA eigenvectors as starting vectors')

      sytol = 10.0d0**(-12)
      fetda = 0.0d0
      ferpa = 0.0d0
      tottda = 0.0d0
      totrpa = 0.0d0
c
c  open fortran files
c
c     now conducted at code outset
c
      if (opg_root()) then
       open (itable,form='unformatted',status='unknown')
       open (inorms,form='unformatted',status='unknown')
       if (modus .le. 2) then
          open (itdspc,file='tda_spectrum',form='formatted',
     +          status='unknown')
          open (itdtex,file='tda_table.tex',form='formatted',
     +          status='unknown')
       end if
       if (modus .ge. 2) then
          open (irpspc,file='rpa_spectrum',form='formatted',
     +          status='unknown')
          open (irptex,file='rpa_table.tex',form='formatted',
     +          status='unknown')
       end if
       if (okeep) open (ivect,form='unformatted',status='unknown')
       if (odebug(20)) then
	ioneel = 20
  	open (ioneel,form='unformatted',status='unknown')
       end if
      endif
c
      if (nfreq .gt. 0) call vclr(pola(1,1),1,6*nfreq)
c
c   initialize files
c
      call rpanal (ipo(iipo+1),nba,dummy,dummy,dummy,dummy,dummy,
     +             dummy,dummy,dummy,dummy,dummy,1,q(itdm))
c
c   perform rpa calculation
c
      if (odirct) then
         call dirrpm (ipo(iipo+1),nba,q(idip),q(icoef),
     +                q(i1),q(i2),q(i3),q(i4),q(i5),
     +                q(i6),q,q,q(itdm))
         call gmem_free_inf(i6,fnm,snm,'scr')
         if(onew) call gmem_free_inf(i5,fnm,snm,'utemp')
         call gmem_free_inf(i4,fnm,snm,'u2')
         call gmem_free_inf(i3,fnm,snm,'u1')
         call gmem_free_inf(i2,fnm,snm,'d2')
         call gmem_free_inf(i1,fnm,snm,'d1')
      else
*        iscr = igmem_alloc(nvc)
         iscr = igmem_alloc_all(nscra)
         call rpmain (ipo(iipo+1),nba,q(idip),q(icoef),q(iscr),q(itdm),
     +                nscra)
         call gmem_free(iscr)
      end if
c
c   output of total oscillator strengths
c
      call outosc
c
c   output of polarizability tensors
c
      if (nfreq .gt. 0) call outpol (nfreq,freq,pola,ipola,gin)
c
c   close spectrum and tex files
c
      call rpanal (ipo(iipo+1),nba,dummy,dummy,dummy,dummy,dummy,
     ?             dummy,dummy,dummy,dummy,dummy,3,q(itdm))
c
      write (iwr,900) cpulft(1)-time ,charwall()
900   format(/'total time for tda/rpa calculation: ',f14.4,' seconds'
     ?        ,a10,' wall'///
     ?	      19x,83('=')/19x,'=',81x,'='/
     ?	      19x,'=    END OF TDA / RPA CALCULATION',49x,'='/
     ?	      19x,'=',81x,'='/19x,83('='))
c
c     free memory
c
       call gmem_free_inf(icoef,fnm,snm,'icoef')
       if( ofiles(1)) call gmem_free_inf(itdm,fnm,snm,'itdm')
       call gmem_free_inf(idip,fnm,snm,'idip')
       call gmem_free_inf(istart,fnm,snm,'ipo')
c
c     close and delete itable and inorms
c
      if(opg_root()) then
        close(unit=itable,status='delete')
        close(unit=inorms,status='delete')
      endif
c
      call end_time_period(TP_RESPONSE)
      call timana(30)
c
      return
      end
      subroutine rpmuld (otda,ipo,nbasis,c,d1,d2,u1,u2,utemp,dcom,x,e,q)

c---------------------------------------------------------------------
c                                       | e |     | A  B | | y |
c   Computes the matrix-vector product  |   |  =  |      | |   |  using the
c                                       | f |     | B  A | | z |
c
c   two-electron integrals in the AO basis.
c   Algorithm: (C. Fuchs, PhD thesis)
c
c                           T        T T
c   1. Compute      d = C yC  - (C zC )
c                        v  c     v  c
c
c      where C is the MO coefficient matrix, divided into core (c) and
c      virtual (v) part, and T denotes transposition.
c
c                                         AO        AO
c   2. Act on d with the integral matrix A  :  u = A  d
c
c             AO
c      where A      = 2 (mn|pq) - (mp|nq)
c             mn,pq
c
c      This is done in routine ubuild.
c
c   3. Obtain e and f as
c
c           T               T T
c      e = C uC       f = -C u C
c           v  c            v   c
c
c   In the TDA procedure (otda = .true.), only  e = Ay  needs to be 
c   calculated.
c
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)

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
c     A common block to pass the actual time periods around in the
c     response theory code. This is needed because in subroutine 
c     such as rpmuld and rpbody the time periods that are used 
c     should be different for RPA and TDA calculations. Hence the
c     usual parameters cannot be used directly.
c
      integer TP_RESP_MO2AO
      integer TP_RESP_INT
      integer TP_RESP_CNTRCT
      integer TP_RESP_AO2MO
      common/response_timep/TP_RESP_MO2AO,TP_RESP_INT,TP_RESP_CNTRCT,
     +                      TP_RESP_AO2MO
c

      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     ?              lcore(8),lvirt(8),iofset(8,8),npair(8),
     ?              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,n1,npp,n2,n3,
     +              np3,npp3
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
      common/infoa /nat,ich,mul,nba,nx,ne,ncore
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension c(nba,nba)
      dimension d1(nba,nba,maxvec),d2(nba,nba,maxvec),
     ?          u1(nba,nba,maxvec),u2(nba,nba,maxvec),
     ?		utemp(maxvec),dcom(*)
      dimension x(isize2,maxvec),e(isize2,maxvec)
      dimension q(*),ipo(nbasis,*)

c
c     calculate d
c
      time = cpulft(1)
      call start_time_period(TP_RESP_MO2AO)
      if (odebug(60)) write (iwr,10)
10    format ('entering rpmuld')
      call vclr(d1,1,nbasq*mvec)
      do 120 m = 1 , mvec
         do 100 i = 1 , nirr
            j = iperm(i,isym)
c...cv*y
            call mxmaa (c(1,lvirt(i)),1,nba,
     ?	                x(iofset(isym,i),m),1,kvirt(i),
     ?		        u1(1,1,m),1,nba,
     ?                  nba,kvirt(i),kcore(j))
c...(cv*y)*cc(t)
            call mxmb  (u1(1,1,m),1,nba,
     ?	                c(1,lcore(j)),nba,1,
     ?		        d1(1,1,m),1,nba,
     ?                  nba,kcore(j),nba)
            if (otda) go to 100
c...cv*z
            call mxmaa (c(1,lvirt(i)),1,nba,
     ?	                x(isize+iofset(isym,i),m),1,kvirt(i),
     ?	   	        u1(1,1,m),1,nba,
     ?                  nba,kvirt(i),kcore(j))
c...cc*(cv*z)(t)
             if(kcore(j).ne.0) then
              call mxmbn (c(1,lcore(j)),1,nba,
     +                    u1(1,1,m),nba,1,
     +                    d1(1,1,m),1,nba,
     +                    nba,kcore(j),nba)
             endif
100      continue
120   continue
      if (odebug(60)) write (iwr,125)
125   format ('finished calculation of d matrix')
      time2 = cpulft(1)
      timing(1) = timing(1) + time2 - time
      if (odebug(31)) then
	 do m=1,mvec
	 write (iwr,126) m
126      format ('trial vector',i3)
         call outsqr (d1(1,1,m),nba,nba,nba,'d after construction')
	 enddo
      end if

      call weight (d1)
      if (onew) then
         do 130 m = 1 , mvec
            call trnsps (d1(1,1,m),d2(1,1,m),nba,nba)
130      continue
      end if
      if (odebug(60)) write (iwr,131)
131   format ('finished weighting and transposing d matrix')
      call dircod (d1,d2,dcom)
      time = cpulft(1)
      timing(2) = timing(2) + time - time2

c   calculate the matrix u

      call vclr(u1,1,nbasq*mvec)
      if (onew) call vclr(u2,1,nbasq*mvec)
      call end_time_period(TP_RESP_MO2AO)
      call dirdef (q,d1,d2,u1,u2,utemp,dcom)
      call start_time_period(TP_RESP_AO2MO)
      time2 = cpulft(1)
      timing(3) = timing(3) + time2 - time
      if (onew) then
         do 140 m = 1 , mvec
            call trnadd (u1(1,1,m),u2(1,1,m),nba)
140      continue
      end if
      call weight (u1)
      if (odebug(32)) then
	 do m=1,mvec
	 write (iwr,141) m
141      format ('u matrix ',i3)
	 call outsqr (u1(1,1,m),nba,nba,nba,'non symmetrized u')
	 enddo
      endif

c   symmetrize u

      do 150 m = 1 , mvec
         call symtru (u1(1,1,m),d1(1,1,m))
         if (odebug(32)) then
	    write (iwr,151) m
151         format ('u matrix ',i3)
	    call outsqr (u1(1,1,m),nba,nba,nba,'symmetrized u')
         end if
150   continue
      time = cpulft(1)
      timing(4) = timing(4) + time - time2

c   calculate (e,f)

      call vclr(e,1,isize2*mvec)
      do 220 m = 1 , mvec
         do 200 i = 1 , nirr
            j = iperm(i,isym)
c...u*cc
            call mxmaa (u1(1,1,m),1,nba,
     ?	                c(1,lcore(j)),1,nba,
     ?	   	        d1(1,1,m),1,nba,
     ?		        nba,nba,kcore(j))
c...cv(t)*(u*cc)
            call mxmb  (c(1,lvirt(i)),nba,1,
     ?	                d1(1,1,m),1,nba,
     ?		        e(iofset(isym,i),m),1,kvirt(i),
     ?                  kvirt(i),nba,kcore(j))
            if (otda) go to 200
c...u(t)*cc
            call mxmaa (u1(1,1,m),nba,1,
     ?	                c(1,lcore(j)),1,nba,
     ?		        d1(1,1,m),1,nba,
     ?		        nba,nba,kcore(j))
c...cv(t)*(u(t)*cc)
            if (kvirt(i).ne.0.and.kcore(j).ne.0) then
             call mxmbn (c(1,lvirt(i)),nba,1,
     +                   d1(1,1,m),1,nba,
     +                   e(isize+iofset(isym,i),m),1,kvirt(i),
     +                   kvirt(i),nba,kcore(j))
            endif
200      continue
         if (odebug(29)) then
            write (iwr,201) m
201         format ('trial vector ',i3)
            length = isize2
            if (otda) length = isize
            call outvec (e(1,m),length,'... before fock contributions')
         endif
220   continue
      timing(5) = timing(5) + cpulft(1) - time

c   fock eigenvalue contributions

      do 500 m = 1 , mvec
         do 400 ii = 1 , ncore
            do 300 ia = nc1 , nba
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
	          k1 = ipo(ii,ia)
                  diff = evalue(ia) - evalue(ii)
                  e(k1,m) = e(k1,m) + diff * x(k1,m)
                  if (otda) go to 250
                  k2 = k1 + isize
                  e(k2,m) = e(k2,m) + diff * x(k2,m)
250               continue
               end if
300         continue
400      continue
500   continue

      call end_time_period(TP_RESP_AO2MO)
      return
      end
      subroutine tdalgd (ipo,nbasis,dip,c,d1,d2,u1,u2,utemp,scr,x,e,v,
     +			 iblb,ible,iblbb,iblee,ibleig,iconv,
     +			 ared,eta,xred,resarr,q)
c---------------------------------------------------------------------
c   Performs the iterative calculation of the lowest eigenvalues
c   plus eigenvectors of the TDA matrix A using the Davidson
c   algorithm.
c   Adapted from routine tdalgo for the direct implementation.
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character*80 dmpfil,rstfil
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
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
c
c     A common block to pass the actual time periods around in the
c     response theory code. This is needed because in subroutine 
c     such as rpmuld and rpbody the time periods that are used 
c     should be different for RPA and TDA calculations. Hence the
c     usual parameters cannot be used directly.
c
      integer TP_RESP_MO2AO
      integer TP_RESP_INT
      integer TP_RESP_CNTRCT
      integer TP_RESP_AO2MO
      common/response_timep/TP_RESP_MO2AO,TP_RESP_INT,TP_RESP_CNTRCT,
     +                      TP_RESP_AO2MO
c
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      common/infoa /nat,ich,mul,nba,nx,ne,ncore
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/rpacom/modus,nevlo(8),nevhi(8),nev(8),
     +              maxr,maxit,maxrr,maxitt,nbegin,
     +              eps,epss,epsana,epstab,
     +              sytol,tolmy,fetda,ferpa,tottda,totrpa,ipref,iblo13
c
      integer itdspc, irpspc, itdtex, irptex, itable, inorms, ivect
      integer itmfil, ldump, idump, lrest, irest, mrest
      common /rpfile/ itdspc,irpspc,itdtex,irptex,itable,inorms,ivect,
     +                itmfil,ldump,idump,lrest,irest,mrest
c
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      common/table /dmpfil,rstfil,rlines(2)
c
c
      dimension dip(3*nvc),c(nbasis,nbasis),ipo(nbasis,*)
      dimension d1(*),d2(*),u1(*),u2(*),utemp(*),scr(*)
      dimension x(isize2,maxvec),e(isize2,maxvec),v(isize2)
      dimension iblb(maxr),ible(maxr),iblbb(maxr),iblee(maxr),
     ?          ibleig(neig),iconv(neig)
      dimension ared(maxr,maxr),eta(maxr),xred(maxr,maxr)
c
      dimension resarr(*),q(*)
c
      data small/1.0d-9/
      data ieigen,iresid / 20,21 /
c
      call start_time_period(TP_TDA)
      TP_RESP_MO2AO  = TP_TDA_MO2AO
      TP_RESP_INT    = TP_TDA_INT
      TP_RESP_CNTRCT = TP_TDA_CNTRCT
      TP_RESP_AO2MO  = TP_TDA_AO2MO
c
      ntot1 = 1
      ntot = nstart
10    niter = niter + 1
      orest = oskip .and. ntot.eq.nstart .and. niter.gt.1
      if (niter .gt. maxit) go to 1000
      if (ldump .gt. 0) then
        if (opg_root()) then
         open (idump,file=dmpfil(1:ldump),form='unformatted',
     +   status = 'unknown')
         rewind idump
        end if
      end if

c   calculate the vectors a * b(i) and fill up matrix ared
c
c      m0    ... index of first trial vector in the current batch
c      m1    ... index of last trial vector in the current batch
c      mvec  ... number of trial vectors in the current batch
c      m,n   ... true indices of trial vectors
c      i,j   ... temporary indices of trial vectors
c		 (with respect to the current batch)

      m0 = ntot1
20    m1 = min(m0+maxvec-1,ntot)
      mvec = m1 - m0 + 1
      m = m0
      do 30 i = 1 , mvec
         if (odebug(43)) write (iwr,1100) 30,'trial',m,iblb(m)
         call rdedx_less (x(1,i),isize,iblb(m),num8)
         if (odebug(30)) then
            print*,'trial vector ',m,':'
            call outvec (x(1,i),isize,'x')
         end if
         if (orest) then
            call rdedx_less (e(1,i),isize,ible(m),num8)
         end if
	 m = m + 1
30    continue
      if (orest) go to 35

c   the next piece of code adds the dipole moment vector(s) to the
c   trial vector(s) if the total oscillator strength is to be calculated

      mvec0 = mvec
      if (ototal .and. idone.lt.inozer .and. mvec0.lt.maxvec) then
         if (odebug(77)) write (iwr,301) isym,m0,m1,mvec0,maxvec
301      format ('additional space in symmetry ',i1/
     ?           'when calculating trial vectors ',i2,' to ',i2/
     ?           'mvec is ',i2,', and maxvec is ',i4)
         do 32 mv = 1 , maxvec-mvec0
            mvec = mvec + 1
            idone = idone + 1
            if (odebug(77)) write (iwr,302) nozer(idone),mvec
302         format ('adding my vector ',i1,' as trial vector ',i2)
            ibegin = iofsym(isym) + (nozer(idone)-1)*nvc + 1
            call dcopy(isize,dip(ibegin),1,x(1,mvec),1)
            if (idone .eq. inozer) go to 33
32       continue
33       continue
      end if

c   compute product vectors

      time = cpulft(1)
      call rpmuld (.true.,ipo,nbasis,c,d1,d2,u1,u2,utemp,scr,x,e,q)
      tim(2) = tim(2) + cpulft(1) - time

c   again worry about total oscillator strength

      if (mvec .gt. mvec0) then
         do 34 mv = mvec0+1 , mvec
            tottda = tottda + ddot(isize,x(1,mv),1,e(1,mv),1)
34       continue
      end if
      mvec = mvec0

35    m = m0
      do 60 i = 1 , mvec
         if (.not. orest) call wrt3 (e(1,i),isize,ible(m),num8)
         do 40 n = 1, m0-1
            if (odebug(43)) write (iwr,1100) 40,'trial',n,iblb(n)
            call rdedx_less(v,isize,iblb(n),num8)
            ared(n,m) = ddot(isize,v,1,e(1,i),1)
            ared(m,n) = ared(n,m)
40       continue
	 n = m0
         do 50 j = 1 , i
            ared(n,m) = ddot(isize,x(1,j),1,e(1,i),1)
            ared(m,n) = ared(n,m)
            n = n + 1
50       continue
         if (odebug(17)) 
     +       call outsqr (ared,maxr,ntot,ntot,'ared after loop 50')
         m = m + 1
60    continue
      if (m1 .eq. ntot) go to 70
      m0 = m0 + maxvec
      go to 20

c   solve reduced eigenvalue problem

70    time = cpulft(1)
      call tred2(maxr,ntot,ared,eta,v,xred)
      call tql2(maxr,ntot,eta,v,xred,ierr)
      if (ierr .ne. 0) go to 400
      tim(3) = tim(3) + cpulft(1) - time

c   output of eigenvalues

      call outeig (itable,iwr,eta(nevlo1+1),neig,niter,ntot,nevlo1,1)
      if (odebug(42)) then
         print*,'niter = ',niter
         call outvec (eta,neig,'eigenvalues')
      end if
      open (ieigen,file='roots', status = 'unknown')
      call outeig (itable,ieigen,dummy,neig,niter,ntot,nevlo1,2)
      close (ieigen)

c   new trial vectors to be expected

      ntot0 = ntot
      ntot1 = ntot + 1
      call vfill(-2.0d0,resarr,1,neig)
      do 230 j = 1 , nevhi(isym)
         if (j .lt. nevlo(isym)) then
            if (ldump.eq.0 .or. orest) then
               go to 230
            else
               go to 145
            end if
         end if
         if (iconv(j-nevlo1) .eq. 1) then
	    resarr(j-nevlo1) = -1.0d0
	    if (ldump.eq.0 .or. orest) go to 230
         end if

c   compute residue vector

145      if (orest) then
            call rdedx (x(1,1),isize,iblbb(j),num8)
            call rdedx (e(1,1),isize,iblee(j),num8)
         else
            call vclr(x(1,1),1,isize)
            call vclr(e(1,1),1,isize)
            do 170 k = 1 , ntot0
               if (odebug(43)) write (iwr,1100) 170,'trial',k,iblb(k)
               call rdedx_less (v,isize,iblb(k),num8)
               call daxpy(isize,xred(k,j),v,1,x(1,1),1)
               if (j .lt. nevlo(isym)) go to 170
               if (odebug(43)) write (iwr,1100) 170,'result',k,ible(k)
               call rdedx_less (v,isize,ible(k),num8)
               call daxpy(isize,xred(k,j),v,1,e(1,1),1)
170         continue
            if (ldump .gt. 0) then
              if (opg_root()) then
               write (idump) (x(i,1),i=1,isize)
              end if
            end if
            if (j.lt.nevlo(isym) .or. iconv(j-nevlo1).eq.1) go to 230
         end if

c   r = omega * x - e

         call daxmyz (isize,eta(j),x(1,1),e(1,1),e(1,1))
         if (odebug(47)) then
            write (iwr,61) j
61          format ('residue vector no.',i3,':')
            call outvec (e(1,1),isize,' ')
         end if
         t = dnrm2(isize,e(1,1),1)
         resarr(j-nevlo1) = t

c   threshold ?

         if (t .le. eps) then
	    iconv(j-nevlo1) = 1
	    call wrt3 (x(1,1),isize,ibleig(j-nevlo1),num8)
	    go to 230
         end if
         if (ntot.ge.maxr .or. ntot.ge.isize) go to 230
         ntot = ntot + 1

c   update residue vector 

         do 190 ii = 1 , ncore
            do 180 ia = nc1 , nba
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
	          k1 = ipo(ii,ia)
	          diff = evalue(ia) - evalue(ii)
	          denom1 =  eta(j) - diff
	          if (dabs(denom1) .gt. small) then
	             e(k1,1) = e(k1,1) / denom1
                  else
	             e(k1,1) = 1.0d0
                  end if
               end if
180         continue
190      continue

c   gram-schmidt orthogonalisation

         do 200 k = 1 , ntot-1
            if (odebug(43)) write (iwr,1100) 200,'trial',k,iblb(k)
            call rdedx_less(x(1,1),isize,iblb(k),num8)
            call orth (isize,e(1,1),x(1,1))
200      continue
         t = dnrm2(isize,e(1,1),1)
         if (1.0d0+t .eq. 1.0d0) then
            write (iwr,210) j,niter
210         format ('discard new trial vector no. ',i2,' in iteration',
     ?              ' step ',i2/'because of linear dependence from',
     ?              ' the previous trial vectors.')
            ntot = ntot - 1
         else
            call dscal(isize,1.0d0/t,e(1,1),1)
	    call wrt3 (e(1,1),isize,iblb(ntot),num8)
         end if
230   continue
      if (ldump .gt. 0) then
        if (opg_root()) then
         close (idump)
        end if
      end if
      call outres (iwr,inorms,resarr,nevlo1,neig,niter,1)
      if (opg_root()) then
       open (iresid,file='residues',status='unknown')
      endif
      call outres (iresid,inorms,dummy,nevlo1,neig,niter,2)
      if (opg_root()) close (iresid)
      nconv = isum(neig,iconv,1)
      if (nconv .eq. neig) go to 1000

c   if at least one new trial vector has been added,
c   start next iteration step

      if (ntot .ge. ntot1) go to 10
      if (ntot .ge. maxr) then
         do 260 j = 1 , nevhi(isym)
            call vclr(x(1,1),1,isize)
            do 250 k = 1 , ntot0
               if (odebug(43)) write (iwr,1100) 250,'trial',k,iblb(k)
               call rdedx_less (v,isize,iblb(k),num8)
               call daxpy(isize,xred(k,j),v,1,x(1,1),1)
250         continue
            if (odebug(43)) write (iwr,117) j,iblbb(j)
117         format('tdalgd, loop 260: writing restart trial vector ',i2,
     ?              ' to block ',i12)
            call wrt3 (x(1,1),isize,iblbb(j),num8)
            if (odebug(44)) 
     +          call outvec (x(1,1),isize,'restart trial vector')
            call vclr(e(1,1),1,isize)
            do 255 k = 1 , ntot0
               if (odebug(43)) write (iwr,1100) 255,'product',k,ible(k)
               call rdedx_less (v,isize,ible(k),num8)
               call daxpy(isize,xred(k,j),v,1,e(1,1),1)
255         continue
            if (odebug(43)) write (iwr,119) j,iblee(j)
119         format (
     +       'diralg, loop 260: writing restart product vector ',i2,
     ?       ' to block ',i12)
            call wrt3 (e(1,1),isize,iblee(j),num8)
            if (odebug(44)) 
     +          call outvec (e(1,1),isize,'restart product vector')
260      continue
         call end_time_period(TP_TDA)
         return
      end if

c   error messages

      write (iwr,270) niter
270   format('no new trial vector has been added in iteration step ',i3,
     ? ':'/'unable to finish iterative calculation of tda eigenvalues.')
      ierr = 1
      go to 1000

400   write (iwr,410) niter,ierr
410   format (/'convergence problems with TQL2 in iteration step ',
     ?        i3,': ierr = ',i3/)

1000  niter = min(niter,maxit)

      call end_time_period(TP_TDA)
      return
1100  format ('loop ',i5,': trying to read ',a8,' vector ',i2,
     ?        ' from block ',i12)
      end
      subroutine dircod (d,dd,dcom)

c---------------------------------------------------------------------
c   Calculates the compressed density matrix.
c   Adapted from GAMESS routine dnsprm.
c   d is the true density matrix.
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8  (a-h,p-z), integer (i-n), logical (o)
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
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,n,ns,np,nsp,
     ?              n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd
      common/gen   /maxit,iter,concrit,conv,energy,thresh,dcmax
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension d(n,n,maxvec),dd(n,n,maxvec),dcom(n,n)

      dcmax = 0.0d0
      if (odebug(35)) 
     +    call outsqr (d(1,1,1),n,n,n,'first d matrix in dircod')
      do 30 m = 1 , mvec
         do 20 i = 1 , n
            do 10 j = 1 , i
	       t = 2.0d0 * (dabs(d(i,j,m)) + dabs(d(j,i,m)))
               if (onew) go to 5
	       dd(i,j,m) = t
	       dd(j,i,m) = t
 5	       dcmax = dmax1(dcmax,t)
10          continue
20       continue
30    continue
      if (odebug(40)) print*,'dcmax = ',dcmax

      call vclr(dcom,1,n*n)
      nsnp = ns + np
c      do 70 m=1,mvec
c      do 60 i=1,n
c	 do 50 j=1,n
c	    dcom(i,j,m) = 1.0d0
c50       continue
c60    continue
c70    continue
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
         do 120 j = 1 , nsnp + nd
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
               do 100 jb = jvl1 , jvl2
		  do 90 m = 1 , mvec
                     dcom(i,j) = max(dcom(i,j),dabs(d(ib,jb,m)))
  90              continue
 100           continue
 110        continue
 120     continue
 130  continue

      if (odebug(97)) call outsqr(dcom,n,n,n,'codens')
      return
      end
      subroutine dircon (iwksp,iwksp2,labsh,nwauq)

c---------------------------------------------------------------------
c   Converts some arrays to another format which is used in the
c   direct implementation. Adapted from GAMESS routine drccon.
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character *4 char
      logical glog1,glog2,glog5,glog6,glog9,gdbg,gmp2
      logical readmo,writmo,restar,acvary,grdt,gopt,gdma,gdold,gscf
      logical gmull

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

      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/data  /charg(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),ndira,ndirb,ndirs,ndirp,nsp,n1,
     +              ndirpp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ndire,repen
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd,ndp,idatbf(maxat),idatp(maxat),
     ?              ijdd1(numspl),ijdd2(numspl),n18,ipar2d,ipar3d,
     ?              idefs1(numspl),idefs2(numspl),idefp1(numspl),
     +              idefp2(numspl),
     ?              idefd1(numspl),idefd2(numspl),ipars,iparp,ipard,
     ?              nsspl,npspl,ndspl
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
      common/gen   /maxitd,iter,concrit,conv,energy,thresh,dcmax,
     ?              glog1,glog2,glog5,glog9,glog6,gdbg(10),
     +              gmp2,mp2fo,mp2fv,
     ?              readmo,writmo,restar,iaccur,acvary,grdt,gopt,gdma,
     +              gdold,
     ?              gscf,noptd,npoint,timeg,gmull,iswap,irhess
      common/gen2  /zipftd(21,2),zijkld(21,2),zjpftd(21,2),
     ?              zipfal(21,2),zijkal(21,2),
     ?              prefac(21),xinner(21),zifall,ziffew,zicfb1,zicfb2
      common/gtrans/idtra(maxorb)
      common/infoa /nat,ich,mul,num,nx,ne,ncore,nb,czan(maxat),
     +              c(3,maxat),amas(maxat),imass(maxat),nuct(maxat)
      common/infob /czanr(maxat),czin(maxat),cin(3,maxat),amass(maxat),
     ?              c80(maxat3,3),nonsym,map80(maxat),ozmat
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      common/phycon/toang,phyc(83),iunt,nosym
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      common/restar/nprint,itol,icut,normf,normp,nopk,
     ?              irest,nrec,intloc,ist,jst,kst,lst,jspace(682)
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
c
      common/scfopt/maxit,mconv,nconv,npunch,accdi1,accdi2,odiis,
     ?              icoupl(3),dmpcut,acurcy,en,etot,ehf,ehf0,diff,iterg,
     ?              icount,rshift,exttol,dmptol,vshtol,iextin,
     ?              iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp
      common/scra  /iso(maxorb,48),ntd,nauq,ict(maxat,48)
      common/scra7 /iblks,iblkt,iblkf,iblkst,iblktt,iblkft,iblkqs,
     ?              isp(18),kblkla
      common/scrtch/ptremp(3,144),dtremp(6,288)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
      common/tol   /toler,tol2
      common/transf/psmal,qsmal,rsmal,pnew,qnew,rnew,pp,qp,rp
      common/work  /jrec,jump
      common/wrtd  /iblk(30),num88,index
 
      dimension char(21)
      dimension iwksp(maxorb,48),iwksp2(maxorb,48),labsh(*)
      dimension nwauq(maxat)
      dimension itrans(6)

      data itrans/1,4,6,2,3,5/
      data char/'ssss','psss','psps','ppss','ppps','pppp',
     ?          'dsss','dsps','dspp','dsds','dpss','dpps','dppp',
     ?          'dpds','dpdp','ddss','ddps','ddpp','ddds','dddp','dddd'/

c
c****** fgrid forms the grid form which the incomplete
c****** gamma function is evaluated via a taylor series
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
      if (odebug(1)) write (iwr,*) 'setting symmetry labels for shells'
      ntd = nt
      nav = lenwrd()
      call readi(labsh,nw196(6)*nav,ibl196(6),idaf)
      do 90 i = 1 , nat
         if (odebug(2)) write(iwr,*)'atom ',i
         do 80 j = 1 , nt
            ict(i,j) = labsh(i+ilisoc(j))
 80      continue
            if (odebug(2)) write(iwr,*)(ict(i,j),j=1,nt)
 90   continue
      call readi(labsh,nw196(5)*nav,ibl196(5),idaf)
      if (odebug(2)) write(iwr,*)'nshell',nshell
      do 120 i = 1 , nshell
      if (odebug(2)) write(iwr,*)'shell ',i
         do 110 j = 1 , nt
            iso(i,j) = labsh(i+iliso(j))
 110     continue
         if (odebug(2)) write(iwr,*)(iso(i,j),j=1,nt)
 120  continue
*
* convert symmetry arrays to direct format
*
      do 60 i = 1 , nonsym
         do 50 it = 2 , nt
            nnew = ict(i,it)
            if (nwauq(i).eq.0 .and. nnew.ne.i) nwauq(nnew) = 1
 50      continue
 60   continue
c
*****
c
c******* rearrange the data so first s, then p, and then d are 
c******* considered
c******* loop over the s functions only *****
c
      icount = 0
      jcount = 0
      do 150 i = 1 , nshell
      if (odebug(2)) write(iwr,*) 'copying symmetry for shell i'
         if (ktype(i).gt.3) call caserr(
     +                     'the program cannot handle f or g functions'
     +                     )
         if (kmin(i).eq.1) then
            jcount = jcount + 1
            do 130 k = 1 , nt
               iwksp(jcount,k) = iso(i,k)
 130        continue
            labsh(i) = jcount
            kk = katom(i)
            nc(jcount) = kng(i)
            isatbf(kk) = isatbf(kk) + 1
            if (odebug(2)) write(iwr,*)'copying s shell info'
            do 140 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
               isatp(kk) = isatp(kk) + 1
               if (odebug(2)) write(iwr,*)icount,ex(k),cs(k)
               pe(icount) = ex(k)
               pc(icount) = cs(k)
 140        continue
         end if
 150  continue
      nsp = icount
      ndirs = jcount
      if (odebug(2)) write(iwr,*)'storing shell data for s',ndirs
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
            if (odebug(2)) write(iwr,*)'copying p shell info'
            do 190 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
              if (odebug(2)) write(iwr,*)icount,ex(k),cp(k)
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
            if (odebug(2)) write(iwr,*)'copying d shell info'
            do 240 k = kstart(i) , kstart(i) + kng(i) - 1
               icount = icount + 1
               nwa(icount) = kk
               idatp(kk) = idatp(kk) + 1
              if (odebug(2)) write(iwr,*)icount,ex(k),cd(k)
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
c      if (num .le. 20) then
c	 do nfunc=1,num
c	 invrs(idtra(nfunc)) = nfunc
c	 end do
c	 call outive(invrs,num,'invrs')
c      end if
      if (odebug(2)) then
         write (iwr,6030)
	 write (iwr,6040) (idtra(i),i=1,num)
      end if
      if (odebug(3)) then
         write (iwr,*) '**************************'
         write (iwr,*) 'switching off symmetry now'
         write (iwr,*) '**************************'
         nt = 1
	 nirr = 1
	 call setsto(num,1,isymmo)
      end if
c
c****** rearrange data in ixatbf and ixatp so that only unique atoms
c******* count
c
      nauq = nat
      if (nt.gt.1) then
         call rdedx(ptr,nw196(1),ibl196(1),idaf)
         call rdedx(dtr,nw196(2),ibl196(2),idaf)
         call rdedx(ftr,nw196(3),ibl196(3),idaf)
         call rdedx(gtr,nw196(4),ibl196(4),idaf)
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
      if (odebug(2)) write(iwr,*)'ndirs,ndirp,ndird',ndirs,ndirp,ndird
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
          if (odebug(2)) write(iwr,*)'ne num nat',ne,num,nat
      ndire = ne
      ndirb = num
      ndira = nat
c
c
      imaxzz = 21
      if (ndird.eq.0) imaxzz = 6
          if (odebug(2)) write(iwr,*)'imaxzz ',imaxzz
      if (odebug(1)) write (iwr,6060)
      do 370 i = 1 , imaxzz
         if (odebug(1)) write (iwr,6010) char(i) , prefac(i) , xinner(i)
 370  continue
      if (odebug(1)) write (iwr,6020) ndirs , nsp , ndirp*3 , ndirpp*3 ,
     +                             ndird*6 , ndirdp*6 , ndirb , nprm
      if (odebug(1)) write (iwr,6030)
      if (odebug(1)) write (iwr,6040) (idtra(i),i=1,num)
      if (odebug(1)) write (iwr,6050)
      if (odebug(1)) then
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
      if (odebug(1) .and. nt.gt.1) then
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
c     nsq = ndirb*ndirb
c     nsq2 = ndirb*(ndirb+1)/2
c
      ileft = isspac(4)
      xavail = dsqrt(dble(ileft))
      ipard = aint(xavail)
      iparp = aint(xavail)
      ipars = aint(xavail)
      call countd(1,idefs1,idefs2,isatbf,isatp,ipars,nsspl,nauq)
      if (odebug(1)) then
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
      call countd(ndirs+1,idefp1,idefp2,ipatbf,ipatp,iparp,
     ?             npspl,nauq)
      if (odebug(1)) then
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
      call countd(ndirs+ndirp+1,idefd1,idefd2,idatbf,
     ?             idatp,ipard,ndspl,nauq)
      if (odebug(1)) then
         write (iwr,6150) ndspl
         write (iwr,6110) ipard
         write (iwr,6080)
         write (iwr,6090) (idefd1(j),j=1,ndspl)
         write (iwr,6100)
         write (iwr,6090) (idefd2(j),j=1,ndspl)
         write (iwr,6130)
      end if

 470  call timit(0)
      if (odebug(1)) then
         write (iwr,6120)
         write (iwr,6130)
         write (iwr,6130)
      end if

      return
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
 6110 format ('maximum number of primitive shells in each batch       ',
     +        i5)
 6120 format ('finished reading and dealing with the data')
 6130 format (//)
 6140 format (1x,'number of batches of p shells       ',i5)
 6150 format (1x,'number of batches of d shells       ',i5)
      end
      subroutine restrt (otda,b,e,v,iblb,ible,iblbb,iblee)

c---------------------------------------------------------------------
c   Performs restart of TDA or RPA calculation when the maximal
c   dimension of the reduced matrices has been reached
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
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
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/orthoc/sytol
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      dimension b(*),e(*),v(*),iblb(*),ible(*),iblbb(*),iblee(*)

      if (otda) then
         if (oskip) then
            do 200 j = 1 , nstart
               if (odebug(43)) write (iwr,152) j,iblbb(j)
               call rdedx (b,isize,iblbb(j),num8)
               if (odebug(44)) 
     +             call outvec (b,isize,'restart trial vector')
               if (odebug(43)) write (iwr,153) j,iblee(j)
               call rdedx (e,isize,iblee(j),num8)
               if (odebug(44)) 
     +             call outvec (e,isize,'restart product vector')
               do 100 k = 1 , j-1
                  call rdedx (v,isize,iblb(k),num8)
                  t = - ddot(isize,b,1,v,1) 
                  call daxpy(isize,t,v,1,b,1)
                  call rdedx (v,isize,ible(k),num8)
                  call daxpy(isize,t,v,1,e,1)
100            continue
               s = dnrm2(isize,b,1)
               if (1.0d0+s .eq. 1.0d0) 
     +         call caserr ('singularity in restart')
               s = 1.0d0 / s
               call dscal(isize,s,b,1)
               call dscal(isize,s,e,1)
               call wrt3 (b,isize,iblb(j),num8)
               call wrt3 (e,isize,ible(j),num8)
               if (odebug(45)) call outvec (b,isize,'restart vector')
200         continue
         else
            do 400 j = 1 , nstart
               call rdedx (b,isize,iblbb(j),num8)
               do 300 k = 1 , j-1
                  call rdedx (v,isize,iblb(k),num8)
                  call orth (isize,b,v)
300            continue
               s = dnrm2(isize,b,1)
               if (1.0d0+s .eq. 1.0d0) 
     +         call caserr ('singularity in restart')
               s = 1.0d0 / s
               call dscal(isize,s,b,1)
               call wrt3 (b,isize,iblb(j),num8)
               if (odebug(45)) call outvec (b,isize,'restart vector')
400         continue
         end if
      else
         i1 = 1
         i2 = 1 + isize
         if (oskip) then
            do 700 j = 1 , nstart
               if (odebug(43)) write (iwr,152) j,iblbb(j)
               call rdedx (b,isize2,iblbb(j),num8)
               if (odebug(44)) 
     +         call outvec (b,isize2,'restart trial vector')
               if (odebug(43)) write (iwr,153) j,iblee(j)
               call rdedx (e,isize2,iblee(j),num8)
               if (odebug(44)) 
     +         call outvec (e,isize2,'restart product vector')
               do 500 k = 1 , j-1
                  call rdedx (v,isize2,iblb(k),num8)
                  t1 = - ddot(isize,b(i1),1,v(i1),1) 
     ?		       - ddot(isize,b(i2),1,v(i2),1)
                  call daxpy(isize,t1,v(i1),1,b(i1),1)
                  call daxpy(isize,t1,v(i2),1,b(i2),1)
                  t2 = - ddot(isize,b(i1),1,v(i2),1) 
     ?	   	       - ddot(isize,b(i2),1,v(i1),1)
                  call daxpy(isize,t2,v(i2),1,b(i1),1)
                  call daxpy(isize,t2,v(i1),1,b(i2),1)
                  call rdedx (v,isize2,ible(k),num8)
                  call daxpy(isize,t1,v(i1),1,e(i1),1)
                  call daxpy(isize,t1,v(i2),1,e(i2),1)
                  call daxpy(isize,t2,v(i2),1,e(i1),1)
                  call daxpy(isize,t2,v(i1),1,e(i2),1)
500            continue
               s = ddot(isize,b(i1),1,b(i1),1) 
     +           + ddot(isize,b(i2),1,b(i2),1)
               if (1.0d0+s .eq. 1.0d0) 
     +         call caserr ('singularity in restart')
               t = 2.0d0 * ddot(isize,b(i1),1,b(i2),1)
               if (s-dabs(t) .lt. sytol) then
                  call caserr ('singularity in restart symort call')
               end if
               alpha = dsqrt( s + t )
               beta  = dsqrt( s - t )
               do 600 k1 = 1 , isize
                  k2 = k1 + isize
                  y = b(k1)
                  z = b(k2)
                  b(k1) = y + z
                  b(k2) = y - z
                  y = b(k1) / alpha
                  z = b(k2) / beta
                  b(k1) = 0.5d0 * ( y + z )
                  b(k2) = 0.5d0 * ( y - z )
                  y = e(k1)
                  z = e(k2)
                  e(k1) = y + z
                  e(k2) = y - z
                  y = e(k1) / alpha
                  z = e(k2) / beta
                  e(k1) = 0.5d0 * ( y + z )
                  e(k2) = 0.5d0 * ( y - z )
600            continue
               call wrt3 (b,isize2,iblb(j),num8)
               call wrt3 (e,isize2,ible(j),num8)
               if (odebug(45)) call outvec (b,isize2,'restart vector')
700         continue
         else
            do 900 j = 1 , nstart
               call rdedx (b,isize2,iblbb(j),num8)
               do 800 k = 1 , j-1
                  call rdedx (v,isize2,iblb(k),num8)
                  t = - ddot(isize,b(i1),1,v(i1),1) 
     ?		      - ddot(isize,b(i2),1,v(i2),1)
                  call daxpy(isize,t,v(i1),1,b(i1),1)
                  call daxpy(isize,t,v(i2),1,b(i2),1)
                  t = - ddot(isize,b(i1),1,v(i2),1) 
     ?		      - ddot(isize,b(i2),1,v(i1),1)
                  call daxpy(isize,t,v(i2),1,b(i1),1)
                  call daxpy(isize,t,v(i1),1,b(i2),1)
800            continue
               t = dnrm2(isize2,b,1)
               if (1.0d0+t .eq. 1.0d0) 
     +         call caserr ('singularity in restart')
               call symort (isize,b(i1),b(i2),sytol,oldep)
               if (oldep) 
     +         call caserr ('singularity in restart symort call')
               call wrt3 (b,isize2,iblb(j),num8)
               if (odebug(45)) call outvec (b,isize2,'restart vector')
900         continue
         end if
      end if
      return
152   format (
     + 'trying to read restart trial vector   ',i2,' from block ',i7)
153   format (
     + 'trying to read restart product vector ',i2,' from block ',i7)
      end
      subroutine rpdipo (ipo,nbasis,dip,coef,q1,q2,vmy)

c---------------------------------------------------------------------
c   Reads in the dipole moment integrals, transforms them to the MO basis
c   and constructs the vectors my(x), my(y), my(z) needed for the
c   calculation of the oscillator strengths
c
c   dip ..... the dipole moment integrals
c   coef .... the MO coefficients
c   q1,q2 ... scratch vectors
c   vmy ..... the vectors my(x),my(y),my(z)
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-y), integer (i-n), logical (o)
      implicit character*8 (z)

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

      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/block /iky(maxorb+1),ilifq(maxorb+1)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     ?              lcore(8),lvirt(8),iofset(8,8),npair(8),
     ?              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
      common/infoa /nat,ich,mul,nba,nx,ne,ncore,nee,czan(maxat),
     +              c(3,maxat),amass(maxat),iannn(2,maxat),
     +              ipseud(maxat),symz(maxat),lpseud
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension o1e(6)
      dimension dip(*),coef(*),q1(*),q2(*),vmy(*),ipo(nbasis,*)

      data o1e/.false.,.false.,.false.,.true.,.true.,.true./

c   read in dipole integrals

      call getmat(dip(1),dip(1),dip(1),dip(1),dip(nbatri+1),
     +            dip(nbatri*2+1),potnuc,nba,o1e,isect(492))

      icoor = 1
      imy = 0
      do 2500 icoord = 1 , 3

c   transform dipole integrals to mo basis

         call square (q1,dip(icoor),nba,nba)
         call mxmaa (q1,1,nba, coef,1,nba, q2,1,nba, nba,nba,nba)
         call mxmaa (coef,nba,1, q2,1,nba, q1,1,nba, nba,nba,nba)
         if (odebug(85)) then
            print*,'dipole integrals, coordinate ',icoord
            call outsqr (q1,nba,nba,nba,' ')
         end if

c   cut out the dipole integrals <a|my|i> needed for rpa

         do 100 ii = 1 , ncore
            do 90 ia = nc1 , nba
               k = iperm(isymmo(ii),isymmo(ia))
               vmy(imy+iofsym(k)+ipo(ii,ia)) = q1(ilifq(ii)+ia)
90          continue
100      continue

         icoor = icoor + nbatri
         imy = imy + nvc
2500  continue

      return
      end
      subroutine symtru (u,usym)

c--------------------------------------------------------------------------
c   Symmetrizes the skeleton "Fock" matrix u, with essentially the
c   technique described by M. Dupuis and H. F. King,
c   Int. J. Quantum Chem. 11, 613 (1977).
c   Adapted from GAMESS routine symhd.
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)

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

      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,n1,npp,n2,n3,
     +              np3,npp3
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
      common/hsym  /t(20,20),mini,maxi,lit,minj,maxj,ljt,ntr
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common/molsym/tr(3,3),trx,try,trz,index,jaxis,igroup,ispmol,
     ?              prmoms(3),aprev(maxat3,3)
      common/scra  /iso(maxorb,48)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      dimension u(nb,nb),usym(nb,nb)
      dimension mi(48)
      dimension ichar(8,22) ! 1st index symm.op., second index irrep
      dimension ioff(19)

      data ioff / 0, 0, 0, 0, 0, 2, 6, 10, 14, 0, 0, 
     +            0, 0, 0, 0, 0, 0, 0, 0 /

      data ichar/ 1, 1, 0, 0, 0, 0, 0, 0,      ! A'        Cs
     ?            1,-1, 0, 0, 0, 0, 0, 0,      ! A''
c--------------------------------------------------------------------------
     ?            1, 1, 1, 1, 0, 0, 0, 0,      ! Ag
     ?            1, 1,-1,-1, 0, 0, 0, 0,      ! Bg        C2h
     ?            1,-1, 1,-1, 0, 0, 0, 0,      ! Bu 
     ?            1,-1,-1, 1, 0, 0, 0, 0,      ! Au
c--------------------------------------------------------------------------
     ?            1, 1, 1, 1, 0, 0, 0, 0,      ! A1
     ?            1,-1,-1, 1, 0, 0, 0, 0,      ! B1        C2v
     ?            1,-1, 1,-1, 0, 0, 0, 0,      ! B2
     ?            1, 1,-1,-1, 0, 0, 0, 0,      ! A2
c--------------------------------------------------------------------------
     ?            1, 1, 1, 1, 0, 0, 0, 0,      ! A1
     ?            1,-1, 1,-1, 0, 0, 0, 0,      ! B2        D2
     ?            1, 1,-1,-1, 0, 0, 0, 0,      ! B1
     ?            1,-1,-1, 1, 0, 0, 0, 0,      ! B3
c--------------------------------------------------------------------------
     ?		  1, 1, 1, 1, 1, 1, 1, 1,      ! Ag
     ?            1,-1, 1,-1,-1, 1,-1, 1,      ! B1u
     ?            1,-1, 1,-1, 1,-1, 1,-1,      ! B2u
     ?            1, 1, 1, 1,-1,-1,-1,-1,      ! B3g       D2h
     ?            1, 1,-1,-1,-1,-1, 1, 1,      ! B3u
     ?            1,-1,-1, 1, 1,-1,-1, 1,      ! B2g
     ?            1,-1,-1, 1,-1, 1, 1,-1,      ! B1g
     ?            1, 1,-1,-1, 1, 1,-1,-1/      ! Au

      if (nt .eq. 1) return
      nsnp = ns + np
      nshell = ns + np + nd
      call vclr(usym,1,nb*nb)
c      print*,'nb = ',nb
c      call outsqr (ptr,3,3,144,'ptr')
c
c   find a block (i,j)

      do 240 ii = 1 , nshell
         do 30 itr = 1 , nt
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
 40      do 230 jj = 1 , nshell
            if (odebug(4)) print*,'examining block ',ii,jj
            do 60 itr = 1 , nt
               ish = mi(itr)
               jsh = iso(jj,itr)
               if (odebug(4)) print*,'(ish,jsh) = ',ish,jsh
               if (ish.eq.ii .and. jsh.gt.jj) go to 230
 60         continue
            if (odebug(4)) print*,'representative block',ii,jj
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
 70         continue
c     ----- find the equivalent blocks
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
            do 140 itr = 1 , nt
               ntr = itr
               kk = mi(itr)
               ll = iso(jj,itr)
               if (odebug(4)) 
     +         print*,'   symm. op. ',itr,' -- equivalent block ',kk,ll
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
                  do 100 l = minj , maxj
                     lcl = locl + l
                     t(k,l) = u(lck,lcl)
 100              continue
 110           continue
               if (odebug(4)) then
        	  print*,'block from u:'
                  do ik=mini,maxi
                  write (iwr,1101) (t(ik,il),il=minj,maxj)
1101              format (20f14.9)
                  end do
               end if
               if (lit.gt.1 .or. ljt.gt.1) call rhrd
               if (odebug(4)) then
                  print*,'rotated block:'
                  do ik=mini,maxi
                  write (iwr,1101) (t(ik,il),il=minj,maxj)
                  end do
               end if
               do 130 i = mini , maxi
                  lci = loci + i
                  if (ichar(itr,ioff(igroup)+isym) .gt. 0) then
                     do 120 j = minj , maxj
                        lcj = locj + j
                        usym(lci,lcj) = usym(lci,lcj) + t(i,j)
 120                 continue
                  else
                     do 125 j = minj , maxj
                        lcj = locj + j
                        usym(lci,lcj) = usym(lci,lcj) - t(i,j)
 125                 continue
                  end if
 130           continue
 140        continue
            if (odebug(4)) then
                print*,'final block:'
                do ik=mini,maxi
                   write (iwr,1101) (usym(loci+ik,locj+il),il=minj,maxj)
                end do
            end if
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
            do 220 itr = 2 , nt
               kk = mi(itr)
               ll = iso(jj,itr)
               if (kk.eq.ii .and. ll.eq.jj) go to 220
               ntr = itr + 1
               if (ntr.le.nt) then
                  do 150 it = ntr , nt
                     i = mi(it)
                     j = iso(jj,it)
                     if (i.eq.kk .and. j.eq.ll) go to 220
 150              continue
               end if
               ntr = invt(itr)
               do 170 i = mini , maxi
        	  lci = loci + i
                  do 160 j = minj , maxj
		     lcj = locj + j
                     t(i,j) = usym(lci,lcj)
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
		  if (ichar(ntr,ioff(igroup)+isym) .gt. 0) then
                     do 200 l = minj , maxj
		        lcl = locl + l
		        usym(lck,lcl) = t(k,l)
 200                 continue
		  else
                     do 205 l = minj , maxj
		        lcl = locl + l
		        usym(lck,lcl) = - t(k,l)
 205                 continue
		  end if
 210           continue
 220        continue
 230     continue
 240  continue

      dum = 1.0d0/dble(nt)
      do 260 i=1,nb
	 do 250 j=1,nb
	    u(i,j) = dum * usym(i,j)
250      continue
260   continue

      return
      end
      subroutine test11(codens,dij,dx,tmax1,test2,xklt,igath,zmem
     &,ic,jc,iscat,kmin,kmax
     &,l1,l2,lmax,lmin,nkl,ipsize,ext2,jcount,jc0,jc1)

c--------------------------------------------------------------------------
c   Performs magnitude-based integral test.
c   Virtually identical to GAMESS routine tst11 except that the
c   lack of symmetry of the "density" matrices is taken into account.
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------------

      implicit real*8  (a-h,o-z), integer (i-n)
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
	 d1 = max(dij,
     ?		  dabs(codens(ic,kc)),dabs(codens(kc,ic)),
     ?		  dabs(codens(jc,kc)),dabs(codens(kc,jc)))
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
            zmem(jc1,4) = max(d1,codens(ic,lc),codens(lc,ic),
     ?                           codens(jc,lc),codens(lc,jc),
     ?                           2.0d0*(codens(kc,lc)+codens(lc,kc)))
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
      subroutine test1s(codens,dij,dx,test2,xklt
     &,igath,zmem,ic,jc,iscat,kmin,kmax,m1,m2,m3
     &,l1,l2,lmax,lmin,nkl,ipsize,ext2,ext3,jcount,jc0,jc1)

c--------------------------------------------------------------------------
c   Sorts out small and symmetry-equivalent integrals.
c   Virtually identical to GAMESS routine tst1s except that the
c   lack of symmetry of the "density" matrices is taken into account.
c   You should also read the notes about the efficiency of tst1s
c   in the corresponding source code file.
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------------

      implicit real*8  (a-h,o-z), integer (i-n)
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
      common/scra  /iso(maxorb,48),nt
      common/defunk/cobas(maxorb,3),nd
      logical ext2,ext3
c     logical ext8
c
c***** the k and l contracted loops for inner integral
c***** test and counting purposes only
c
c..... here!
      jcount = 0
      xnt = dble(nt)
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
	 d1 = max(dij,
     ?		  dabs(codens(ic,kc)),dabs(codens(kc,ic)),
     ?		  dabs(codens(jc,kc)),dabs(codens(kc,jc)))
         do 40 lc = lmin1 , lmax1
            jc1 = jc1 + 1
            jcxx = jcxx + nc(kc)*nc(lc)
            zmem(jc1,5) = max(d1,codens(ic,lc),codens(lc,ic),
     ?                           codens(jc,lc),codens(lc,jc),
     ?                           2.0d0*(codens(kc,lc)+codens(lc,kc)))
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
c$dir scalar
            do 70 kprm = 1 , kpmax0
               jcxx = jcxx + 1
               if (dabs(zmem(jcxx,6)).gt.test2) go to 80
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
      subroutine trnadd (a,b,n)
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension a(n,n),b(n,n)
      do 20 i = 1 , n
         do 10 j = 1 , n
            a(i,j) = a(i,j) + b(j,i)
10       continue
20    continue
      return
      end
      subroutine ubuil2 (g,ivl,jvl,kvl,lvl,ic,jc,kmin,kmax,
     ?       lmin,lmax,l1,l2,ext2,
     ?       u1,u2,utemp,d1,d2,gg,memreq,ipos,jpos,kpos,lpos,i1,j1,k1,
     ?       kdim,ldim)

c---------------------------------------------------------------------
c   Forms u matrix for a given i and j shell pairing. 
c   This is the "modern" version which treats several density matrices
c   in one step by means of matrix-matrix multiplications.
c   Adapted from GAMESS routine fockb.
c   ivl/jvl are the number of basis functions in the corresponding shell.
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      logical ext2

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

      common/cider /intgrl(8,100),intcon(5,12),ihz(5,5),
     ?              ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri,
     ?              mempto,memcto,mtemp,ipszt,icszt,
     ?              ipsize,icsize,memcon
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     ?              lcore(8),lvirt(8),iofset(8,8),npair(8),
     ?              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,
     ?              n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,
     +              n6,n7,n8,
     ?              n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension u1(nb,nb,maxvec),u2(nb,nb,maxvec),utemp(maxvec)
      dimension d1(nb,nb,maxvec),d2(nb,nb,maxvec)
      dimension g(icsize,memreq)
      dimension gg(kdim,ldim)

      time = cpulft(1)
c-----------------------------------------------------------------------
c   arrange integrals in matrix gg(k,l)                                |
c-----------------------------------------------------------------------
      morig = 0
      mem = kdim * ldim
      if (mem.gt.65536) then
       write(iwr,*) ' **** increase dimension of gg in rpbody ****'
       write(iwr,*) ' **** size required = ', mem
       call caserr('dimensioning problem in rpbody')
      endif
      mem = 0
      call vclr(gg(1,1),1,kdim*ldim)
      if (odebug(91)) then
         print*,'ic,jc = ',ic,jc
         print*,'ivl,jvl,kvl,lvl = ',ivl,jvl,kvl,lvl 
         print*,'kmin,kmax,lmin,lmax = ',kmin,kmax,lmin,lmax
         print*,'ipos,jpos,kpos,lpos = ',ipos,jpos,kpos,lpos
         print*,'l1,l2,kdim,ldim = ',l1,l2,kdim,ldim
	 print*,'mvec = ',mvec
      end if
      do 90 i = 1 , ivl
         id = (ic-i1)*ivl + i + ipos
         if (odebug(91)) print*,'   i,id = ',i,id
         do 80 j = 1 , jvl
            jd = (jc-j1)*jvl + j + jpos
            if (odebug(91)) print*,'      j,jd = ',j,jd
            do 70 k = 1 , kvl
               morig = mem
               mem = mem + 1
               jcount = 0
               do 60 kc = kmin , kmax
                  lmin1 = l1
                  lmax1 = l2
                  if (ext2) lmax1 = kc
                  if (kc .eq. kmax) lmax1 = lmax
                  if (kc .eq. kmin) lmin1 = lmin
                  kd = (kc-kmin)*kvl + k
                  if (odebug(91)) print*,'         kc,kd = ',kc,kd
                  if (odebug(91)) 
     +            print*,'         lmin1,lmax1 = ',lmin1,lmax1
                  if (lx .eq. 1) then
		     ld21 = lmax1 - lmin1 + 1
                     call dcopy(ld21
     +               ,g(jcount+1,mem),1,gg(kd,lmin1-l1+1),kdim)
		     jcount = jcount + ld21
                  else
                     mem = morig
                     l0 = (lmin1-l1-1)*lvl
                     do 40 l = 1 , lvl
                        jc0 = jcount
                        lcount = l + l0
                        mem = mem + 1
                        do 30 lc = lmin1 , lmax1
                           lcount = lcount + lvl
                           jc0 = jc0 + 1
                           gg(kd,lcount) = g(jc0,mem)
c                           if (odebug(91)) print*,'gg(',kd,lcount,
c                               ') = g(',jc0,mem,') = ',g(jc0,mem)
30                      continue
40                   continue
                     jcount = jc0
c                     ld21 = lmax1 - lmin1 + 1
c                     joff = jcount + 1
c                     loff = (lmin1-l1) * lvl
c                     lstrid = kdim * lvl
c                     do 40 l = 1 , lvl
c                        call dcopy(ld21,g(joff,morig+l),1,
c                                         gg(kd,loff+l),lstrid)
c40                   continue
c                     jcount = jcount + ld21
c                     mem = morig + lvl
                  end if
60             continue
70          continue
            if (odebug(91)) then
               call outsqr (gg,kdim,kdim,ldim,'gg in ubuil2')
            end if

            kstart = (kmin-k1)*kvl + 1 + kpos
	    if (lx .eq. 1) then
	       ld1 = l1
	       ld2 = l2
            else
	       ld1 = 1 + lpos
	       ld2 = (l2-l1)*lvl + lvl + lpos
            end if
            call vclr(utemp,1,mvec)
c-----------------------------------------------------------------------
c   utemp(m) = sum_k,l  [ d1(k,l,m) + d2(k,l,m) ] * g(k,l)             |
c-----------------------------------------------------------------------
            do 120 ld = ld1 , ld2
               call mxmb (d1(kstart,ld,1), nbasq, 1,
     ?                    gg(1,ld-ld1+1), 1, kdim,
     ?                    utemp, 1, mvec,
     ?                    mvec, kdim, 1)                       
               call mxmb (d2(kstart,ld,1), nbasq, 1,
     ?                    gg(1,ld-ld1+1), 1, kdim,
     ?                    utemp, 1, mvec,
     ?                    mvec, kdim, 1)                       
120         continue
c-----------------------------------------------------------------------
c   u(i,j,m) <-- u(i,j,m) + 2 * utemp(m)                               |
c   u(j,i,m) <-- u(j,i,m) + 2 * utemp(m)                               |
c-----------------------------------------------------------------------
            call daxpy(mvec, 2.0d0, utemp, 1, u1(id,jd,1), nbasq)
            call daxpy(mvec, 2.0d0, utemp, 1, u2(id,jd,1), nbasq)
c-----------------------------------------------------------------------
c   utemp(m) = 2 * [ d(i,j,m) + d(j,i,m) ]                             |
c-----------------------------------------------------------------------
            call vadd(d1(id,jd,1)
     +                ,nbasq, d2(id,jd,1), nbasq, utemp, 1, mvec)
            call dscal(mvec, 2.0d0, utemp, 1)
c----------------------------------------------------------------------
c   u(k,l,m) <-- u(k,l,m) + g(k,l) * utemp(m)                          |
c   u(l,k,m) <-- u(l,k,m) + g(k,l) * utemp(m)                          |
c-----------------------------------------------------------------------
            do 130 ld = ld1 , ld2
               call mxmb (gg(1,ld-ld1+1), 1, kdim,
     ?                    utemp, mvec, 1,
     ?                    u1(kstart,ld,1), 1, nbasq,
     ?                    kdim, 1, mvec)
               call mxmb (gg(1,ld-ld1+1), 1, kdim,
     ?                    utemp, mvec, 1,
     ?                    u2(kstart,ld,1), 1, nbasq,
     ?                    kdim, 1, mvec)
130         continue
c            kd = kstart
c            do 150 k = 1 , kdim
c               do 140 m = 1 , mvec
c                  utemp = utemp(m)
c                  ld = ld1
c                  do 130 l = 1 , ldim
c                     prod = gg(k,l) * utemp
c                     u1(kd,ld,m) = u1(kd,ld,m) + prod
c                     u2(kd,ld,m) = u2(kd,ld,m) + prod
c                     ld = ld + 1
c130               continue
c140            continue
c               kd = kd + 1
c150         continue
c----------------------------------------------------------------------
c   u(i,k,m) <-- u(i,k,m) - g(k,l) * d(j,l,m)                         |
c----------------------------------------------------------------------
            call mxmbn (gg, 1, kdim,
     ?                  d2(ld1,jd,1), 1, nbasq,
     ?                  u2(kstart,id,1), 1, nbasq,
     ?                  kdim, ldim, mvec)
c----------------------------------------------------------------------
c   u(j,k,m) <-- u(j,k,m) - g(k,l) * d(i,l,m)                         |
c----------------------------------------------------------------------
            call mxmbn (gg, 1, kdim,
     ?                  d2(ld1,id,1), 1, nbasq,
     ?                  u2(kstart,jd,1), 1, nbasq,
     ?                  kdim, ldim, mvec)
c-----------------------------------------------------------------------
c   u(i,l,m) <-- u(i,l,m) - g(k,l) * d(j,k,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, kdim, 1,
     ?                  d2(kstart,jd,1), 1, nbasq,
     ?                  u2(ld1,id,1), 1, nbasq,
     ?                  ldim, kdim, mvec)
c-----------------------------------------------------------------------
c   u(j,l,m) <-- u(j,l,m) - g(k,l) * d(i,k,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, kdim, 1,
     ?                  d2(kstart,id,1), 1, nbasq,
     ?                  u2(ld1,jd,1), 1, nbasq,
     ?                  ldim, kdim, mvec)
c-----------------------------------------------------------------------
c   u(k,i,m) <-- u(k,i,m) - g(k,l) * d(l,j,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, 1, kdim,
     ?                  d1(ld1,jd,1), 1, nbasq,
     ?                  u1(kstart,id,1), 1, nbasq,
     ?                  kdim, ldim, mvec)
c-----------------------------------------------------------------------
c   u(k,j,m) <-- u(k,j,m) - g(k,l) * d(l,i,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, 1, kdim,
     ?                  d1(ld1,id,1), 1, nbasq,
     ?                  u1(kstart,jd,1), 1, nbasq,
     ?                  kdim, ldim, mvec)
c-----------------------------------------------------------------------
c   u(l,i,m) <-- u(l,i,m) - g(k,l) * d(k,j,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, kdim, 1,
     ?                  d1(kstart,jd,1), 1, nbasq,
     ?                  u1(ld1,id,1), 1, nbasq,
     ?                  ldim, kdim, mvec)
c-----------------------------------------------------------------------
c   u(l,j,m) <-- u(l,j,m) - g(k,l) * d(k,i,m)                          |
c-----------------------------------------------------------------------
            call mxmbn (gg, kdim, 1,
     ?                  d1(kstart,id,1), 1, nbasq,
     ?                  u1(ld1,jd,1), 1, nbasq,
     ?                  ldim, kdim, mvec)
80       continue
90    continue
      timing(6) = timing(6) + cpulft(1) - time
      return
      end
      subroutine ubuild (g,ivl,jvl,kvl,lvl,ic,jc,kmin,kmax,
     ?           lmin,lmax,l1,l2,ext2,
     ?           u,d,dd,temp,memreq,ipos,jpos,kpos,lpos,i1,j1,k1)

c--------------------------------------------------------------------------
c   Forms u matrix for a given i and j shell pairing.
c   This is the old version which computes each u matrix separately.
c   Adapted from GAMESS routine fockb.
c   ivl/jvl are the number of basis functions in the corresponding shell.
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      logical ext2

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

      common/cider /intgrl(8,100),intcon(5,12),ihz(5,5),
     ?              ix,jx,kx,lx,mx,ijx,klx,ninty,ncon,memans,mempri,
     ?              mempto,memcto,mtemp,ipszt,icszt,
     ?              ipsize,icsize,memcon
      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),nc(maxorb),
     ?              isatbf(maxat),isatp(maxat),ipatbf(maxat),
     +              ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,
     ?              n1,npp,n2,n3,np3,npp3,nprm,n4,nm1,nm2,ne,repen,n5,
     +              n6,n7,n8,
     ?              n9,n10,n11,n12,n13,n14,n15,n16,n17
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr,
     +              timing(10)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension u(nb,nb,maxvec),d(nb,nb,maxvec),dd(nb,nb,maxvec)
      dimension g(icsize,memreq),temp(icsize*6)

      time = cpulft(1)
      if (odebug(39)) then
         print*,'ic,jc = ',ic,jc
         print*,'ivl,jvl,kvl,lvl = ',ivl,jvl,kvl,lvl 
         print*,'kmin,kmax,lmin,lmax,l1,l2 = ',kmin,kmax,lmin,lmax,l1,l2
         print*,'ipos,jpos,kpos,lpos = ',ipos,jpos,kpos,lpos
	 print*,'mvec = ',mvec
      end if
      if (odebug(34)) then
         jcomax=0
         do kc=kmin,kmax
            lmin1=l1
            lmax1=l2
            if (ext2) lmax1=kc
            if (kc.eq.kmax) lmax1=lmax
            if (kc.eq.kmin) lmin1=lmin
            jcomax=jcomax+lmax1-lmin1+1
         enddo
         memmax=ivl*jvl*kvl*lvl
         m=1
         n=8
510      if (memmax.lt.m) go to 560
         n=min(n,memmax)
         do 520 j=1,jcomax
         write(iwr,550) ic,jc,ivl,jvl,kvl,lvl,j,m,n,(g(j,i),i=m,n)
520      continue
         m=m+8
         n=n+8
         goto 510
550      format(2i3,4i2,3i4,8f14.7)
560      continue
c	 call outsqr (g,icsize,icsize,memreq,'g in ubuild')
      end if
      mem = 0
      morig = 0
      do 90 i = 1 , ivl
         id = (ic-i1)*ivl + i + ipos
         do 80 j = 1 , jvl
            jd = (jc-j1)*jvl + j + jpos
            do 70 k = 1 , kvl
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
                  if (lx.eq.1) then
                     ld1 = lmin1
                     ld2 = lmax1
                     ld21 = ld2 - ld1 + 1
                     jstart = jcount + 1
		     do 210 m = 1 , mvec
                        jcount = jstart - 1
c		        t = d(id,jd,m) + d(jd,id,m)
c		        t = t + t
                        t = dd(id,jd,m)
		        do 120 ld = ld1 , ld2
                           jcount = jcount + 1
                           ge = g(jcount,mem)
			   u(ld,kd,m) = u(ld,kd,m) + ge * t
			   u(ld,id,m) = u(ld,id,m) - ge * d(kd,jd,m)
			   u(ld,jd,m) = u(ld,jd,m) - ge * d(kd,id,m)
			   u(kd,ld,m) = u(kd,ld,m) + ge * t
			   u(id,ld,m) = u(id,ld,m) - ge * d(jd,kd,m)
			   u(jd,ld,m) = u(jd,ld,m) - ge * d(id,kd,m)
120                     continue                        
		t = ddot(ld21,dd(ld1,kd,m),1,g(jstart,mem),1)
		u(id,jd,m) = u(id,jd,m) + t
		u(jd,id,m) = u(jd,id,m) + t
	        u(id,kd,m) = u(id,kd,m) - 
     ?  	ddot(ld21,d(jd,ld1,m),nb,g(jstart,mem),1)
		u(kd,id,m) = u(kd,id,m) -
     ?		ddot(ld21,d(ld1,jd,m),1,g(jstart,mem),1)
		u(jd,kd,m) = u(jd,kd,m) -
     ?		ddot(ld21,d(id,ld1,m),nb,g(jstart,mem),1)
		u(kd,jd,m) = u(kd,jd,m) -
     ?		ddot(ld21,d(ld1,id,m),1,g(jstart,mem),1)
210                  continue
                  else
                     mem = morig
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
		     do 310 m = 1 , mvec
                        jc0 = 0
                        t = dd(id,jd,m)
		        do 150 ld = ld1 , ld2
                           jc0 = jc0 + 1
                           ge = temp(jc0)
			   u(ld,kd,m) = u(ld,kd,m) + ge * t
			   u(ld,id,m) = u(ld,id,m) - ge * d(kd,jd,m)
			   u(ld,jd,m) = u(ld,jd,m) - ge * d(kd,id,m)
			   u(kd,ld,m) = u(kd,ld,m) + ge * t
			   u(id,ld,m) = u(id,ld,m) - ge * d(jd,kd,m)
			   u(jd,ld,m) = u(jd,ld,m) - ge * d(id,kd,m)
150                     continue                        
		        t = ddot(ld21,dd(ld1,kd,m),1,temp(1),1)
		        u(id,jd,m) = u(id,jd,m) + t
		        u(jd,id,m) = u(jd,id,m) + t
		        u(id,kd,m) = u(id,kd,m) - 
     ?			ddot(ld21,d(jd,ld1,m),nb,temp(1),1)
		        u(kd,id,m) = u(kd,id,m) -
     ?			ddot(ld21,d(ld1,jd,m),1,temp(1),1)
		        u(jd,kd,m) = u(jd,kd,m) -
     ?			ddot(ld21,d(id,ld1,m),nb,temp(1),1)
		        u(kd,jd,m) = u(kd,jd,m) -
     ?			ddot(ld21,d(ld1,id,m),1,temp(1),1)
310                  continue
                  end if
60             continue
70          continue
80       continue
90    continue
      timing(6) = timing(6) + cpulft(1) - time
      return
      end
      subroutine weight (d)

c--------------------------------------------------------------------------
c   Code is weighting the dxz dxy dyz elements of the density matrix
c   by root 3 (see GAMESS routine dnsprm).
c   (c) Carsten Fuchs 1993
c--------------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)

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

      common/data  /charge(maxat),coord(3,maxat),nwa(mxprim),
     +              nc(maxorb),isatbf(maxat),isatp(maxat),
     +              ipatbf(maxat),ipatp(maxat),
     ?              nwshp(mxprim),nup(maxorb),pe(mxprim),pc(mxprim),
     ?              coprm(mxprim,3),na,nb,ns,np,nsp,n1,npp,n2,n3,np3,
     +              npp3
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/defunk/cobas(maxorb,3),nd
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      dimension d(nb,nb,maxvec)


      if (nd .eq. 0) return
      sq3 = dsqrt(3.0d0)
      ic = 0
      do 30 i = 1 , nb
         weigh1 = 1.0d0
         ix = i - ns - np3
         if (ix.gt.0) then
            ic = ic + 1
            if (ic.eq.2 .or. ic.eq.3 .or. ic.eq.5) weigh1 = sq3
            if (ic.eq.6) ic = 0
         end if
         jc = 0
         do 20 j = 1 , nb
            weigh2 = 1.0d0
            jx = j - ns - np3
            if (jx.gt.0) then
               jc = jc + 1
               if (jc.eq.2 .or. jc.eq.3 .or. jc.eq.5) weigh2 = sq3
               if (jc.eq.6) jc = 0
            end if
            do 10 m = 1 , mvec
               d(i,j,m) = d(i,j,m) * weigh1 * weigh2
10          continue
20       continue
30    continue

      return
      end
      function char2i (i)
c--------------------------------------------------------------------------
c   Converts a number with at most two digits to a
c   character string
c   (c) Carsten Fuchs 1991-1993
c--------------------------------------------------------------------------
      character*2 c,char2i
      character*1 digit(0:9)
      data digit/'0','1','2','3','4','5','6','7','8','9'/
      num = iabs(i)
      char2i = '  '
      if (num .ge. 100) return
      l = 2
      do 1 j=2,1,-1
      new = num / 10
      n = num - 10 * new
      if (n .gt. 0) l = j
      c(j:j) = digit(n)
1     num = new
      goto (2,3),l
2     char2i(1:2) = c(1:2)
      go to 5
3     char2i(1:1) = c(2:2)
      char2i(2:2) = ' '
5     return
      end
      subroutine ver_dirrpa(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/dirrpa.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
