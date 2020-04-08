      subroutine addone(p,d,ioff,isymrs)
c---------------------------------------------------------------------
c   Adds the one-electron part to the two-particle density matrix
c   (code obsolete, you are invited to improve it!)
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z),logical(o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension p(*),d(*),ioff(*)
c...for vector computers: no recurrence in this loop
      do 10 itux=1,nact3
      it=(itux-1)/nact2+1
      iux=itux-(it-1)*nact2
      iu=(iux-1)/nact+1
      ix=iux-ilifa(iu)
      if(mult(itypea(it),itypea(ix)).eq.isymrs)then
      k=ipos(ilifa(it)+iu)+ioff(iux)
      p(k)=p(k)-d(ipos(ilifa(it)+ix))
      endif
10    continue
      return
      end
      subroutine adiagc (diaga,diagt,zint,zb,z,fi)
c
c---------------------------------------------------------------------
c   Constructs the configuration diagonal of the matrix A
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension diaga(*),diagt(*)
      dimension zint(*),zb(maxbb,nact),z(maxbb,nact),fi(*)
      dimension icga(31),icgb(31),za(31),ioff(8)

      iccc = 0
      if (odebug(43)) write (iwr,1) ibzint
1     format ('adiagc: trying to read zint from block ',i6)
      nints = n0int2
      if (ozint) nints = nint2
c
c     this appears to be reading just 2-electron integrals
c     yet both 1e and 2e-integrals have been written here
c     in lrzint
      call rdedx_less(zint,nints,ibzint,ned7)
      call dscal(nints,2.0d0,zint,1)
      k = 0
      do 10 isyma=1,nirr
      ioff(isyma) = k
10    k = k + nsymm(isyma)**2

c   fill diagonal with constant entries

      call vfill(const,diagt,1,nci1)
      call vfill(const,diaga,1,nci1)
      excapp = excmax
      if (oexca) excapp = excave
      if (odebug(7)) print*,'excapp = ',excapp

      do 170 isymb=1,nirr
      isyma = mult(isymb,isym)
      nsa = nstra(isyma)
      nsb = nstrb(isymb)
      if (nsa*nsb .le. 0) go to 170
      call vclr(zb,1,nact*maxbb)

c   store beta string occupancies

      mtb = 0
      do 40 ib=1,nsb
20    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      if (jmtb .ne. isymb) go to 20
      do 30 ii=1,nb
30    zb(ib,icgb(ii)) = 1.0d0
40    continue

c   loop over alpha strings

      mta = 0
      do 160 ia=1,nsa
50    call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      if (jmta .ne. isyma) go to 50
      call vclr(za,1,nact)
      do 60 ii=1,na
60    za(icga(ii)) = 1.0d0

c   occupation numbers of determinants

      do 70 i=1,nact
      do 70 ib=1,nsb
70    z(ib,i) = zb(ib,i) + za(i)

c   contributions from inactive fock matrix

      do 80 i=nst,nprim
         k = itype(i)
         fitt = fi(ioff(k)+(npoint(i)+mcore(k)-1)*nsymm(k)
     +          +npoint(i)+mcore(k))
         do 80 ib=1,nsb
            diaga(iccc+ib) = diaga(iccc+ib) + z(ib,i-ncore) * fitt
80          diagt(iccc+ib) = diagt(iccc+ib) + z(ib,i-ncore) * fitt

c   coulomb and exchange contribution

      do 110 i=1,nact
      do 110 j=1,i
      ij = ilifa(i) + j
      if(ozint)ij=icf(i)+j
      if(ozint)then
      coul = zint(ic2e((i*(i+1))/2)+ic3e((j*(j+1))/2))
      exch = zint(ic2e(ij)+ic3e(ij))
      else
      coul = zint(ipos(ilifa(i)+i)+ioff0(ilifa(j)+j))
      exch = zint(ipos(ij)+ioff0(ij))
      endif
      zz = coul - exch
      if(j.lt.i) then
         do 90 ib=1,nsb
            val1=z(ib,i)*z(ib,j)*zz
            val2=za(i)*zb(ib,j)+za(j)*zb(ib,i)
c           diaga(iccc+ib)=diaga(iccc+ib)+val1
c           if (odebug(7)) print*,'excapp,exch = ',excapp,exch
            diaga(iccc+ib)=diaga(iccc+ib)+val1+val2*excapp
90          diagt(iccc+ib)=diagt(iccc+ib)+val1+val2*exch
      else
         do 100 ib=1,nsb
         diaga(iccc+ib)=diaga(iccc+ib)+za(i)*zb(ib,i)*exch
100      diagt(iccc+ib)=diagt(iccc+ib)+za(i)*zb(ib,i)*exch
      endif
110   continue
160   iccc = iccc + nsb
170   continue
      return
      end
      subroutine blowup (ntot0,v,b,xred,msplit,iblb)
c-----------------------------------------------------------------------
c   Expands an eigenvector of the reduced MCLR equation to the
c   full space
c   (c) Carsten Fuchs 1991-1993
c-----------------------------------------------------------------------
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension v(*),b(*),xred(maxr2),msplit(maxr),iblb(maxr)
      call vclr(v,1,npar2)
      do 200 k = 1 , ntot0
         n = leng(msplit(k))
         n1 = noff(msplit(k)) - 1
         n2 = n1 + npar
         call rdedx (b,2*n,iblb(k),ned7)
         x1kj = xred(k)
         x2kj = xred(k+ntot0)
c...for vector computers: no recurrence in this loop
         do 100 i = 1 , n
            b1ik = b(i)
            b2ik = b(n+i)
            v(n1+i) = v(n1+i) + b1ik * x1kj + b2ik * x2kj
            v(n2+i) = v(n2+i) + b2ik * x1kj + b1ik * x2kj
100      continue
200   continue
      return
      end
      subroutine canoni (f,d,q)
c
c-----------------------------------------------------------------------
c   Computes the pseudocanonical and natural orbitals
c   (c) Carsten Fuchs 1991-1993
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension f(*),d(*),q(*)
c
      ibase = icorr(0)
      ntrans = 0
      do 10 i = 1 , nirr
         ntrans = ntrans + mcore(i)**2 + mact(i)**2 + mvirt(i)**2
10    continue
      iuu = icorr(ntrans)
c...next line test
      call vclr(q(iuu),1,ntrans)
      kk = 0

c   pseudocanonical core orbitals

      ifok = 1
      do 540 i = 1 , nirr
         m = mcore(i)
         n = nsymm(i)
         if (m .eq. 0) go to 535
         ieig = icorr(n)
         iscr = icorr(n)
         ivec = icorr(n*n)
         call tred2(n,m,f(ifok),q(ieig),q(iscr),q(ivec))
         call tql2(n,m,q(ieig),q(iscr),q(ivec),ierr)
         if (ierr .ne. 0) call caserr ('ierr .ne. 0 in canoni')
         if (odebug(69)) then
            write (iwr,220) i
         do 525 j=1,m
525         write (iwr,260) q(ieig-1+j),(q(ivec-1+(j-1)*n+k),k=1,m)
         endif
         iofcor(i) = kk
         do 530 k = 1 , m
            call fmove (q(ivec+(k-1)*n),q(iuu+iofcor(i)+(k-1)*m),m)
530      continue
         kk = kk + m*m
         call corlsr (ieig)
535      ifok = ifok + n*n
540   continue
c
c   natural orbitals
c
      call dscal(n0int1,-1.0d0,d,1)
      iden = 1
      do 570 i = 1 , nirr
         m = mact(i)
         if (m .eq. 0) go to 570
         iocc = icorr(m)
         iscr = icorr(m)
         ivec = icorr(m*m)
         call tred2(m,m,d(iden),q(iocc),q(iscr),q(ivec))
         call tql2(m,m,q(iocc),q(iscr),q(ivec),ierr)
         if (ierr .ne. 0) call caserr ('ierr .ne. 0 in canoni')
         if (odebug(69)) then
            write (iwr,320) i
         do 550 j=1,m
550         write (iwr,260) q(iocc-1+j),(q(ivec-1+(j-1)*m+k),k=1,m)
         end if
         iofact(i) = kk
         call dcopy(m*m,q(ivec),1,q(iuu+iofact(i)),1)
         kk = kk + m*m
         ioc = iocc
         do 590 it = 1 , nact
            if (itypea(it) .ne. i) go to 590
            ddiag(it) = -q(ioc)
            ioc = ioc + 1
590      continue
         call corlsr (iocc)
         iden = iden + m*m
570   continue
      call dscal(n0int1,-1.0d0,d,1)

c   pseudocanonical virtual orbitals

      ifok = 1
      do 160 i = 1 , nirr
         m = mvirt(i)
         n = nsymm(i)
         if (m .eq. 0) go to 150
         ieig = icorr(n)
         iscr = icorr(n)
         ivec = icorr(n*n)
         ifck = ifok + nprm(i)*n + nprm(i)
         call tred2(n,m,f(ifck),q(ieig),q(iscr),q(ivec))
         call tql2(n,m,q(ieig),q(iscr),q(ivec),ierr)
         if (ierr .ne. 0) call caserr ('ierr .ne. 0 in canoni')
         if (odebug(69)) then
            write (iwr,420) i
            do 130 j=1,m
130         write (iwr,260) q(ieig-1+j),(q(ivec-1+(j-1)*n+k),k=1,m)
         end if
         iofvir(i) = kk
         do 140 k = 1 , m
            call fmove (q(ivec+(k-1)*n),q(iuu+iofvir(i)+(k-1)*m),m)
140      continue
         kk = kk + m*m
         call corlsr (ieig)
150      ifok = ifok + n*n
160   continue
c
      if (kk .ne. ntrans) call caserr ('something wrong in canoni')
c
c   write transformation matrix to block IBCAN
c
      ibcan = iblo7
      call wrt3 (q(iuu),ntrans,ibcan,ned7)
      if (odebug(90)) print*,'output canonical orbitals on blocks',
     ?                ibcan,' -',ibcan+lensec(ntrans)-1
      iblo7 = iblo7 + lensec(ntrans)
      call corlsr (ibase)
c
220   format (/'pseudocanonical core orbitals, symmetry ',i1,':'/)
260   format (f9.4,6x,20f9.4)
320   format (/'natural orbitals, symmetry ',i1,':'/)
420   format (/'pseudocanonical virtual orbitals, symmetry ',i1,':'/)
      return
      end
      subroutine ciadrs (mode,isymm,iofa,iofb,inter,ndim)
c----------------------------------------------------------------------
c   Constructs indexing arrays for the addressing of the
c   CI vectors
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension iofa(*),iofb(*),inter(ndim,nact,2)
      dimension icg(31)
      if (mode .eq. 1) go to 50

c   addresses of intermediate states

      nab = na
      do 40 kk=1,2
      if (kk .eq. 2) nab = nb
      do 20 k=1,nab
      inter(k,k,kk) = 0
      do 10 l=k,nact-1+k-nab
10    inter(k,l+1,kk) = inter(k,l,kk) + ibinom(nact-l,nab-k)
20    continue
      do 30 i=1,nab-1
      do 30 k=1,nact-1
30    inter(i,k,kk) = inter(i,k,kk) - inter(i+1,k+1,kk)
40    continue

c   addressing of ci vector

c   loop over beta string symmetries

50    nci = 0
      do 130 iz=1,nirr

c   first do beta string address offsets (iofb)

      kt = 0
      mt = 0
90    call strong (mt,nact,nb,icg,itypea,jmt,iu,mult)
      if (mt .eq. 0) go to 100
      if (jmt .ne. iz) go to 90
      kt = kt + 1
      iofb(mt) = kt
      go to 90
100   kstrb(iz) = kt

c   now loop over valid alpha strings

      iaz = mult(isymm,iz)
      mt = 0
      kstra(iaz) = 0
110   call strong (mt,nact,na,icg,itypea,jmt,iu,mult)
      if (mt .eq. 0) go to 120
      if (jmt .ne. iaz) go to 110
      iofa(mt) = nci
      kstra(iaz) = kstra(iaz) + 1
      nci = nci + kt
      go to 110
120   continue
130   continue
      if (odebug(81)) 
     +      call outive (iofa,nstraa,'alpha offsets from ciadrs:')
      return
      end
      subroutine coco (yo,yc,zo,zc,vo,vc,wo,wc,a,b)

c---------------------------------------------------------------------
c   Computes the coupling parts (o-c & c-o) of the vector E*x
c   in the case where the matrices A(oc),B(oc) are stored on disk
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension yo(*),yc(*),zo(*),zc(*),vo(*),vc(*),wo(*)
      dimension wc(*),a(*),b(*)

      npiece = (norb-1)/lenrs + 1
      lrs = lenrs
      nrs = 0
      iblka = ibco
      do 300 ipiece=1,npiece
      if (ipiece .eq. npiece) lrs = norb - nrs
      lenco = lrs * nci1
      iblkb = iblka + lensec(lenco)
      call rdedx (a,lenco,iblka,ned7)
      call rdedx (b,lenco,iblkb,ned7)
      do 200 n=1,lrs
      do 100 i=1,nci1
      vo(nrs+n) = vo(nrs+n) + a((n-1)*nci1+i) * yc(i) + 
     +                        b((n-1)*nci1+i) * zc(i)
      wo(nrs+n) = wo(nrs+n) + b((n-1)*nci1+i) * yc(i) +
     +                        a((n-1)*nci1+i) * zc(i)
      vc(i) = vc(i) + a((n-1)*nci1+i) * yo(nrs+n) + 
     +                b((n-1)*nci1+i) * zo(nrs+n)
100   wc(i) = wc(i) + b((n-1)*nci1+i) * yo(nrs+n) + 
     +                a((n-1)*nci1+i) * zo(nrs+n)
200   continue
      iblka = iblka + 2 * lensec(lenco)
300   nrs = nrs + lrs

      return
      end
      subroutine comout (kc,ki,comm,io,jo,oend)
c---------------------------------------------------------------------
c   Writes the entries of A(oo),B(oo) on disk, later they
c   are reordered in routine cosort.
c   The indices pq,rs of A(oo)_pq,rs are stored in packed form.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension comm(*),io(*),jo(*)
      nc = kc - 1
      ni = nc / 2
      if (oend) goto 20
      nblk = nc / 408
      if (nblk .eq. 0) return
      jc = 1
      ji = 1
      do 10 iblk=1,nblk
      call fmove (comm(jc),gin,408)
      call pack22 (204,gin(409),io(ji),jo(ji))
      gin(511) = dble(204)
      call wrt3 (gin,511,iblo4,ned4)
      iblo4 = iblo4 + 1
      jc = jc + 408
10    ji = ji + 204
      ncrest = nc - nblk * 408
      nirest = ni - nblk * 204
      call fmove (comm(jc),comm(1),ncrest)
      call icopy (nirest,io(ji),1,io(1),1)
      call icopy (nirest,jo(ji),1,jo(1),1)
      kc = ncrest + 1
      ki = nirest + 1
      return
20    call vclr(gin,1,511)
      call fmove (comm,gin,nc)
      call pack22 (ni,gin(409),io(1),jo(1))
      gin(511) = dble(ni)
      call wrt3 (gin,511,iblo4,ned4)
      iblo4 = iblo4 + 1
      return
      end
      subroutine copart (isplit,c,d,e1,e2,fv,fw,gv,gw,yc,zc,vc,wc,
     ?                   dc,pc,qc,rc,ioa0,ioa1,iob0,iob1,inter,
     +                   iwa0,iwa1,iwb,
     ?                   nwa0,nwa1,nwb,nablk,nbblk)
c
c----------------------------------------------------------------------
c   Constructs the c-o part of the vector E*x and the transition density
c   matrices D(c),P(c) needed for the o-c part. Usually the most time
c   consuming routine !
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical(o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension c(1),d(1),e1(1),e2(1),fv(1),fw(1),gv(1),gw(1)
      dimension dc(1),pc(1),qc(1),rc(1)
      dimension yc(1),zc(1),vc(1),wc(1)
      dimension ioa0(1),ioa1(1),iob0(1),iob1(1),inter(na,nact,2)
      dimension iwa0(1),iwa1(1),iwb(1),nwa0(1),nwa1(1),nwb(1)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
c
      iuaa = 0
      iubb = 0
      do 600 isyma=1,nirr
      if (nstra(isyma) .eq. 0) go to 600
      nalock = (nstra(isyma)-1)/nablk+1
      do 590 isymb=1,nirr
      if (nstrb(isymb) .eq. 0) go to 590
      nblock = (nstrb(isymb)-1)/nbblk+1
c
c   symmetry of intermediate state
c
      isymk = mult(isyma,isymb)
c
c   symmetry of orbital excitations
c
      isymr0 = mult(isss,isymk)
c
c   symmetry of excited beta strings
c
      isymb0 = mult (isymb,isymr0)
      nstrb0 = nstrb(isymb0)
      mtaa = 0
      do 580 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)
c
c   generate alpha replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isss
c
      time = cpulft(1)
      mta = mtaa
      do 20 i=1,na
20    icga(i) = icgaa(i)
      iua = iuaa
      iofiwa = 1
      do 40 ia=1,naa
30    call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      if (jmta .ne. isyma) go to 30
      if (ia .eq. 1) ia0 = ioa0(mta)
      call rplace (icga,isymr0,iwa0(iofiwa),nwa0(ia),na,ioa0,inter)
40    iofiwa = iofiwa + maxrp3
c
c   generate alpha replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isym
c
      isymr1 = mult(isym,isymk)
      isymb1 = mult(isymb,isymr1)
      nstrb1 = nstrb(isymb1)
      if (isym .eq. isss) then
         call icopy (nablk,nwa0,1,nwa1,1)
         call icopy (maxrp3*nablk,iwa0,1,iwa1,1)
         ia1 = ia0
      else
         mta = mtaa
         do 50 i=1,na
50       icga(i) = icgaa(i)
         iua = iuaa
         iofiwa = 1
         do 70 ia=1,naa
60       call string (mta,nact,na,icga,itypea,jmta,iua,mult)
         if (jmta .ne. isyma) go to 60
         if (ia .eq. 1) ia1 = ioa1(mta)
         call rplace (icga,isymr1,iwa1(iofiwa),nwa1(ia),na,ioa1,inter)
70       iofiwa = iofiwa + maxrp3
         if (odebug(57)) then
            write (iwr,73)isyma,isymb,ialock
73          format ('alpha replacements for symmetries (',i1,',',i1,
     +              ') -',' block ',i2)
            iofiw1 = 0
            do 75 ia=1,naa
            iofiw2 = iofiw1 + maxrpl
            iofiw3 = iofiw2 + maxrpl
            do 72 iw=1,nwa1(ia)
            itu = iwa1(iofiw3+iw)
            it = (itu-1)/nact + 1
            iu = itu - ilifa(it)
            write (iwr,71) ia,it,iu,iwa1(iofiw1+iw),iwa1(iofiw2+iw)
71          format ('<',i2,'/ea_',2i2,'/',i2,'> = ',i3)
72          continue
75          iofiw1 = iofiw1 + maxrp3
         end if
         continue
      end if
      cotim(1) = cotim(1)+cpulft(1)-time
      mtbb = 0
      do 560 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa * nbb
      nabp0 = naabb * npair(isymr0)
      nabp1 = naabb * npair(isymr1)
c
c   generate beta replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isss
c
      time = cpulft(1)
      mtb = mtbb  
      iub = iubb
      do 80 i=1,nb
80    icgb(i)=icgbb(i)
      iofiwb = 1
      do 100 ib=1,nbb
90    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      if (jmtb .ne. isymb) go to 90
      if (ib .eq. 1) ib0 = iob0(mtb) - 1
      call rplace (icgb,isymr0,iwb(iofiwb),nwb(ib),nb,iob0,inter(1,1,2))
100   iofiwb = iofiwb + maxrp3
      cotim(1) = cotim(1) + cpulft(1) - time
c
c   construct d matrix
c
      time = cpulft(1)
      call vclr(d,1,nabp0)
      iofiw1 = 0
      iofca = ib0
      iofd = - naabb
      do 140 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 130 iw=1,nwa0(ia)
      ioffc = iofca + iwa0(iofiw1+iw)
      ioffd = iofd + ipos(iwa0(iofiw3+iw))*naabb
      if (iwa0(iofiw2+iw) .gt. 0) then
         do 110 ib=1,nbb
110      d(ioffd+ib) = d(ioffd+ib) + c(ioffc+ib)
      else
         do 120 ib=1,nbb
120      d(ioffd+ib) = d(ioffd+ib) - c(ioffc+ib)
      end if
130   continue
      iofiw1 = iofiw1 + maxrp3
140   iofd = iofd + nbb
      iofd = 1 - nbb - naabb
      iofcb = ia0 - nstrb0
      iofiw1 = 0
      do 180 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 170 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + ipos(iwb(iofiw3+iw))*naabb
      if (iwb(iofiw2+iw) .gt. 0) then
         do 150 ia=1,naa
150      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c(ioffc+ia*nstrb0)
      else
         do 160 ia=1,naa
160      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c(ioffc+ia*nstrb0)
      end if
170   continue
      iofiw1 = iofiw1 + maxrp3
180   iofd = iofd + 1
      cotim(2) = cotim(2) + cpulft(1) - time
      if (odebug(33)) then
         write (iwr,951) isyma,isymb,ialock,iblock
951      format ('d matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 955 iv=1,nact
         do 955 ix=1,nact
         if (mult(itypea(iv),itypea(ix)) .eq. isymr0) then
            write (iwr,952) iv,ix
952         format ('pair (',i2,',',i2,'):')
            k = (ipos(ilifa(iv)+ix)-1)*naabb
            do 953 ib=1,nbb
953         write (iwr,954) (d(k+(ia-1)*nbb+ib),ia=1,naa)
954         format (10f9.4)
         end if
955      continue
      end if
c
c   generate beta replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isym
c
      time = cpulft(1)
      if (isym .eq. isss) then
         ib1 = ib0
      else
         mtb = mtbb
         iub = iubb
         do 190 i=1,nb
190      icgb(i) = icgbb(i)
         iofiwb = 1
         do 210 ib=1,nbb
200      call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
         if (jmtb .ne. isymb) go to 200
         if (ib .eq. 1) ib1 = iob1(mtb) - 1
         call rplace (icgb,isymr1,iwb(iofiwb),nwb(ib),nb,iob1,
     +                inter(1,1,2))
210      iofiwb = iofiwb + maxrp3
         if (odebug(57)) then
            write (iwr,83)isyma,isymb,iblock
83          format ('beta replacements for symmetries (',i1,',',i1,
     +              ') -',' block ',i2)
            iofiw1 = 0
            do 85 ib=1,nbb
            iofiw2 = iofiw1 + maxrpl
            iofiw3 = iofiw2 + maxrpl
            do 82 iw=1,nwb(ib)
            itu = iwb(iofiw3+iw)
            it = (itu-1)/nact + 1
            iu = itu - ilifa(it)
            write (iwr,81) ib,it,iu,iwb(iofiw1+iw),iwb(iofiw2+iw)
81          format ('<',i2,'/eb_',2i2,'/',i2,'> = ',i3)
82          continue
85          iofiw1 = iofiw1 + maxrp3
         end if
         continue
      end if
      cotim(1) = cotim(1)+cpulft(1)-time
      if (isymk .ne. isym) go to 332
      if (isplit .eq. 2) go to 300
c
c   one electron contributions to vectors vco and wco
c
c      iofd = 0
      time = cpulft(1)
      do 290 it=1,nact
      do 290 iu=1,nact
      if (mult(itypea(it),itypea(iu)) .eq. isymr) then
         iofvw = ia1 + ib1
         itu = ipos(ilifa(it)+iu)
         iofd = (itu-1) * naabb
         fvtu = fv(itu)
         fwtu = fw(itu)
         do 280 ia=1,naa
         do 270 ib=1,nbb
         vc(iofvw+ib) = vc(iofvw+ib) + fvtu * d(iofd+ib)
270      wc(iofvw+ib) = wc(iofvw+ib) + fwtu * d(iofd+ib)
         iofd = iofd + nbb
280      iofvw = iofvw + nstrb(isymb)
      end if
290   continue
      cotim(5) = cotim(5) + cpulft(1) - time
      if (isplit .eq. 1) go to 335
c
c   transition density matrix dc
c
300   continue
      time = cpulft(1)
      do 330 it=1,nact
      do 330 iu=1,nact
      if (mult(itypea(it),itypea(iu)) .eq. isymr) then
         iofyz = ia1 + ib1
         itu = ipos(ilifa(it)+iu)
         iut = ipos(ilifa(iu)+it)
         iofd1 = (itu-1) * naabb
         iofd2 = (iut-1) * naabb
         do 320 ia=1,naa
         do 310 ib=1,nbb
310      dc(itu) = dc(itu) + d(iofd1+ib) * yc(iofyz+ib)
     ?                     - d(iofd2+ib) * zc(iofyz+ib)
         iofd1 = iofd1 + nbb
         iofd2 = iofd2 + nbb
320      iofyz = iofyz + nstrb(isymb)
      end if
330   continue
      cotim(4) = cotim(4) + cpulft(1) - time
332   if (isplit .eq. 2) go to 420
c
c   construct matrices e1 and e2
c
335   continue
      time = cpulft(1)
      call vclr(e1,1,nabp1)
      call vclr(e2,1,nabp1)
      if (odebug(55)) then
      do 7100 it=1,nact
      do 7100 iu=1,nact
      if (mult(itypea(it),itypea(iu)) .ne. isymr1) go to 7100
      do 7090 k=1,naabb
      do 7080 iv=1,nact
      do 7080 ix=1,nact
      if (mult(itypea(iv),itypea(ix)) .ne. isymr0) go to 7080
      e1((ipos(ilifa(it)+iu)-1)*naabb+k) = 
     +                  e1((ipos(ilifa(it)+iu)-1)*naabb+k)
     ?   + d((ipos(ilifa(iv)+ix)-1)*naabb+k) 
     ?   * gv(ipos(ilifa(it)+iu)+ioff1(ilifa(iv)+ix))
      e2((ipos(ilifa(it)+iu)-1)*naabb+k) = 
     +                  e2((ipos(ilifa(it)+iu)-1)*naabb+k)
     ?   + d((ipos(ilifa(iv)+ix)-1)*naabb+k) 
     ?   * gw(ipos(ilifa(it)+iu)+ioff1(ilifa(iv)+ix))
7080  continue
7090  continue
7100  continue
      else
      call mxmaa (d,1,naabb,gv(intof1(isymr0)+1),npair(isymr1),1,
     +            e1,1,naabb,
     ?            naabb,npair(isymr0),npair(isymr1))
      call mxmaa (d,1,naabb,gw(intof1(isymr0)+1),npair(isymr1),1,
     +            e2,1,naabb,
     ?            naabb,npair(isymr0),npair(isymr1))
      end if
      cotim(3) = cotim(3) + cpulft(1) - time
      if (odebug(42)) then
         write (iwr,751) isyma,isymb,ialock,iblock
751      format ('e1 matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 755 iv=1,nact
         do 755 ix=1,nact
         if (mult(itypea(iv),itypea(ix)) .eq. isymr1) then
            write (iwr,752) iv,ix
752         format ('pair (',i2,',',i2,'):')
            k = (ipos(ilifa(iv)+ix)-1)*naabb
            do 753 ib=1,nbb
753         write (iwr,754) (e1(k+(ia-1)*nbb+ib),ia=1,naa)
754         format (10f9.4)
         end if
755      continue
      end if
c
c   contributions to v(co) and w(co)
c
      time = cpulft(1)      
      iofvwa = ib1
      iofe = -naabb
      iofiw1 = 0
      do 370 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 360 iw=1,nwa1(ia)
      iut = iwa1(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      ioffvw = iofvwa + iwa1(iofiw1+iw)
      ioffe = iofe + ipos(ilifa(it)+iu)*naabb
      if (iwa1(iofiw2+iw) .gt. 0) then
         do 340 ib=1,nbb
         vc(ioffvw+ib) = vc(ioffvw+ib) + e1(ioffe+ib)
340      wc(ioffvw+ib) = wc(ioffvw+ib) + e2(ioffe+ib)
      else
         do 350 ib=1,nbb
         vc(ioffvw+ib) = vc(ioffvw+ib) - e1(ioffe+ib)
350      wc(ioffvw+ib) = wc(ioffvw+ib) - e2(ioffe+ib)
      end if
360   continue
      iofiw1 = iofiw1 + maxrp3
370   iofe = iofe + nbb
      iofvwb = ia1 - nstrb1
      iofe = 1 - nbb - naabb
      iofiw1 = 0
      do 410 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 400 iw=1,nwb(ib)
      iut = iwb(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      ioffvw = iofvwb + iwb(iofiw1+iw)
      ioffe = iofe + ipos(ilifa(it)+iu)*naabb
      if (iwb(iofiw2+iw) .gt. 0) then
         do 380 ia=1,naa
         vc(ioffvw+ia*nstrb1) = vc(ioffvw+ia*nstrb1) + e1(ioffe+ia*nbb)
380      wc(ioffvw+ia*nstrb1) = wc(ioffvw+ia*nstrb1) + e2(ioffe+ia*nbb)
      else
         do 390 ia=1,naa
         vc(ioffvw+ia*nstrb1) = vc(ioffvw+ia*nstrb1) - e1(ioffe+ia*nbb)
390      wc(ioffvw+ia*nstrb1) = wc(ioffvw+ia*nstrb1) - e2(ioffe+ia*nbb)
      end if
400   continue
      iofiw1 = iofiw1 + maxrp3
410   iofe = iofe + 1
      cotim(5) = cotim(5) + cpulft(1) - time
      if (isplit .eq. 1) go to 540
c
c   construct matrices e(y),e(z)
c
420   continue
      time = cpulft(1)
      call vclr(e1,1,nabp1)
      call vclr(e2,1,nabp1)
      iofiw1 = 0
      iofyza = ib1
      iofe = -naabb
      do 460 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 450 iw=1,nwa1(ia)
      ioffyz = iofyza + iwa1(iofiw1+iw)
      iut = iwa1(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      ioffe = iofe + ipos(ilifa(it)+iu)*naabb
      if (iwa1(iofiw2+iw) .gt. 0) then
         do 430 ib=1,nbb
         e1(ioffe+ib) = e1(ioffe+ib) + yc(ioffyz+ib)
430      e2(ioffe+ib) = e2(ioffe+ib) + zc(ioffyz+ib)
      else
         do 440 ib=1,nbb
         e1(ioffe+ib) = e1(ioffe+ib) - yc(ioffyz+ib)
440      e2(ioffe+ib) = e2(ioffe+ib) - zc(ioffyz+ib)
      end if
450   continue
      iofiw1 = iofiw1 + maxrp3
460   iofe = iofe + nbb
      iofyzb = ia1 - nstrb1
      iofe = 1 - nbb - naabb
      iofiw1 = 0
      do 500 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 490 iw=1,nwb(ib)
      ioffyz = iofyzb + iwb(iofiw1+iw)
      iut = iwb(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      ioffe = iofe + ipos(ilifa(it)+iu)*naabb
      if (iwb(iofiw2+iw) .gt. 0) then
         do 470 ia=1,naa
         e1(ioffe+ia*nbb) = e1(ioffe+ia*nbb) + yc(ioffyz+ia*nstrb1)
470      e2(ioffe+ia*nbb) = e2(ioffe+ia*nbb) + zc(ioffyz+ia*nstrb1)
      else
         do 480 ia=1,naa
         e1(ioffe+ia*nbb) = e1(ioffe+ia*nbb) - yc(ioffyz+ia*nstrb1)
480      e2(ioffe+ia*nbb) = e2(ioffe+ia*nbb) - zc(ioffyz+ia*nstrb1)
      end if
490   continue
      iofiw1 = iofiw1 + maxrp3
500   iofe = iofe + 1
      cotim(6) = cotim(6)+cpulft(1)-time
      if (odebug(41)) then
         write (iwr,551) isyma,isymb,ialock,iblock
551      format ('matrix e(y) for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 555 it=1,nact
         do 555 iu=1,nact
         if (mult(itypea(it),itypea(iu)) .eq. isymr1) then
            k = (ipos(ilifa(it)+iu)-1)*naabb
            do 553 ib=1,nbb
            do 553 ia=1,naa
            t = e1(k+(ia-1)*nbb+ib)
            if (1.0d0+t .ne. 1.0d0) write (iwr,554) ia,ib,it,iu,t
554         format ('<',2i2,'|E_',2i2,'|n> = ',f14.9)
553         continue
         end if
555      continue
      end if
c
c   contributions to pc and qc
c
      time = cpulft(1)
      call mxmb (e1,naabb,1,d,1,naabb,pc(intof1(isymr0)+1),1,
     +  npair(isymr1),npair(isymr1),naabb,npair(isymr0))
      call mxmb (e2,naabb,1,d,1,naabb,qc(intof1(isymr0)+1),1,
     +  npair(isymr1),npair(isymr1),naabb,npair(isymr0))
      cotim(4) = cotim(4) + cpulft(1) - time
540   continue
c
c   end of beta block
c
      do 550 i=1,nb
550   icgbb(i) = icgb(i)
      iubb = iub
560   mtbb = mtb
c
c   end of alpha block
c
      do 570 i=1,na
570   icgaa(i) = icga(i)
      iuaa = iua
580   mtaa = mta
c
c   end of beta symmetry
c
590   continue
c
c   end of alpha symmetry
c
600   continue
      time = cpulft(1)
      if (isplit .ge. 2) then
         call invers (qc,rc)
         call daxpy(n1int2,-1.0d0,rc,1,pc,1)
      end if
      shtim(7) = shtim(7) + cpulft(1) - time
c
c   correction term for isym=isss
c
      if (isym.eq.isss .and. isplit.ne.2) then
         if (odebug(95)) print*,'fvisum = ',fvisum
         call daxpy(nci0, fvisum,c,1,vc,1) 
         call daxpy(nci0,-fvisum,c,1,wc,1) 
      end if
      if (isplit.ge.2 .and. odebug(36) .and. nact2.le.511) then
         write (iwr,603)
603      format (/'transition density matrix:'/)
         call vclr(gin,1,nact2)
         do 604 i=1,nact
         do 604 j=1,nact
         if (mult(itypea(i),itypea(j)) .eq. isymr) then
            gin(ilifa(i)+j) = dc(ipos(ilifa(i)+j))
         end if
604      continue
         do 605 i=1,nact
605      write (iwr,606)(gin(ilifa(i)+j),j=1,nact)
606      format (10f9.4)
      end if
      if (isplit .ge. 2) call addone (pc,dc,ioff1,isymr)
      if (isplit .ge. 2 .and. odebug(39) .and. niter .le. 2) then
         write (iwr,611)
611      format (/'two particle transition density matrix after ',
     +            ' addition of one electron part:'/)
         call outden (2,isym,ioff1,dc,pc)
      end if
      if (isplit .ge. 2) then
         do 710 iv = 2 , nact
           do 700 ix = 1 , iv-1
           isymvx = mult(itypea(iv),itypea(ix))
           isymtu = mult(isymvx,isymr)
           iaddv = ioff1(ilifa(ix)+iv)+1
           iaddx = ioff1(ilifa(iv)+ix)+1
           call daxpy(npair(isymtu)
     +                ,1.0d0,pc(iaddv),1,pc(iaddx),1)
700        continue
710      continue
      end if
c
      return
      end
      subroutine cosort (q)

c----------------------------------------------------------------------
c   Sorts the entries of the orbital-orbital parts A(oo),B(oo).
c   The unsorted values <0|[E_pq,[H_0,E_rs]]|0> which were constructed
c   in routine OOPART are read from ED4. The sorting procedure is
c   similar to that employed in integral transformations (see
c   G.H.F. Diercksen, Theoret. Chim. Acta 33, 1 (1974)). The ordered
c   entries are columns in the matrices A(oo),B(oo) and are stored
c   on ED7.
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------

      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      logical btest
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      equivalence (dum,m)
      integer*2 itemp(4)
      equivalence (pack,itemp(1))
      data itemp / 4*0 /
c
      dimension q(*),iad(980),ifill(980)
      dimension ijij(2,340)
c
      time = cpulft(1)
      mmax = icorrm() - 1
      ntot = (norb-1) / lp + 1
      if (ntot .gt. 980) call caserr ('too many boxes in COSORT !')
      isp = ntot * 511
      if (odebug(99)) then
         write (iwr,10) norb,lp,ntot,isp,mmax
10       format (/'cosort parameters:'/
     ?            '------------------'/
     ?            'no. of orbital pairs:      ',i8/
     ?            'no. of pairs in one batch: ',i8/
     ?            'no. of batches (boxes):    ',i8/
     ?            'space needed:              ',i8/
     ?            'space available:           ',i8/)
      end if
      if (isp .gt. mmax) 
     ?   call caserr ('insufficient space for commutator reordering')
      call setsto(ntot,0,iad)
      call setsto(ntot,0,ifill)
      call setsto(680,0,ijij)
      call vclr(q,1,isp)
c
c   read in commutators
c
      iblin = 1
      iblout = iblo4
      nabtri = (norb*(norb+1))/2
      ntimes = (nabtri-1)/204 + 1
      ncomm = 204
      do 200 itimes=1,ntimes
      if (itimes .eq. ntimes) ncomm = nabtri-(ntimes-1)*204
      ncomm2 = ncomm + ncomm
      call rdedx (gin,511,iblin,ned4)
      call upak4t (ncomm2,gin(409),ijij)
      kk = 1
      do 100 k=1,ncomm
      ind1 = ijij(1,k)
      ind2 = ijij(2,k)
      irun = 1
50    ibox = (ind1-1)/lp + 1
      if (odebug(99)) print*,'ind1,ind2,ibox = ',ind1,ind2,ibox
      mlower = (ibox-1)*511
      j = ifill(ibox) + 1
      q(mlower+j) = gin(kk)
      q(mlower+204+j) = gin(kk+1)
      kpos = mlower+408+j/2
      if (btest(j,0)) then
         itemp(1) = ind1
         itemp(2) = ind2
         q(kpos+1) = pack
      else
         pack = q(kpos)
         itemp(3) = ind1
         itemp(4) = ind2
         q(kpos) = pack
      end if
      if (j .ge. 204) then
         q(mlower+511) = dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblout,ned4)
         iad(ibox) = iblout*1000 + 204
         iblout = iblout + 1
         call vclr(q(mlower+1),1,511)
         ifill(ibox) = 0
      else
         ifill(ibox) = j
      end if
      if (irun .eq. 2) go to 100
      ind1 = ijij(2,k)
      ind2 = ijij(1,k)
      irun = 2
      go to 50
100   kk = kk+2
200   iblin = iblin + 1

c   finally write partially filled blocks onto external file

      do 120 ibox=1,ntot
      j = ifill(ibox)
      if (j .eq. 0) go to 120
      mlower = (ibox-1)*511
      if (btest(j,0)) then
         kpos = mlower + 340 + j/2 + 1
         pack = q(kpos)
         itemp(3) = 0
         itemp(4) = 0
         q(kpos) = pack
      end if
      q(mlower+511) = dble(iad(ibox))
      call wrt3 (q(mlower+1),511,iblout,ned4)
      iad(ibox) = iblout*1000 + j
      iblout = iblout + 1
120   continue

c   reordering step 2

      if (lp*norb2 .gt. mmax) 
     ?   call caserr ('insufficient space for reordering step 2')
      ibloab = iblo7
      ipair = 0
      np = lp
      do 650 ibox=1,ntot
      if (ibox .eq. ntot) np = norb - ipair
      npp = np * norb
      npp2 = npp + npp
      call vclr(q,1,npp2)
      ii = iad(ibox)
600   iblo = ii/1000
      if (iblo .eq. 0) go to 620
      m = ii - 1000*iblo
      call rdedx (gin,511,iblo,ned4)
      call upak4t (408,gin(409),ijij)
c...for vector computers: no recurrence in this loop
      do 610 k = 1 , m
         kpos = (ijij(1,k)-ipair-1) * norb + ijij(2,k)
         q(    kpos) = gin(    k)
         q(npp+kpos) = gin(204+k)
610   continue
      ii = idint(gin(511))
      go to 600
620   call wrt3 (q(1),npp2,iblo7,ned7)
      if (odebug(97)) then
         write (iwr,623) ipair+1,ipair+np
623      format (/'reordered commutators, pairs ',i3,' to ',i3,' :'//
     ?            'A commutators'/)
         do 624 i=1,np
624      call outvec (q((i-1)*norb+1),norb,' ')
         write (iwr,625)
625      format (/'B commutators'/)
         do 626 i=1,np
626      call outvec (q(npp+(i-1)*norb+1),norb,' ')
      end if
      ipair = ipair + np
650   iblo7 = iblo7 + lensec(npp2)
c
      write (iwr,300) cpulft(1)-time
300   format ('reordering of commutators ',28('.'),f9.4,' seconds')
c
      return
      end
      subroutine couplx (a,b,c,d,e,f,g,q,w,
     ?                   ioa0,ioa1,iob0,iob1,inter,iwa0,iwa1,iwb,
     ?                   nwa0,nwa1,nwb,nablk,nbblk)

c---------------------------------------------------------------------
c   Constructs the entries of the matrices A(co),B(co) and writes them
c   to disk. This is not the default because of extensive I/O.
c   a --- pieces of A(co)
c   b --- pieces of B(co)
c   c --- reference CI vector
c   d --- matrix D, with D(k,vx) = <k|E_vx|0>
c   e --- matrix E, with E(k,rs) = sum_vx  (rs|vx) * D(k,vx)
c   f --- inactive Fock matrix
c   g --- modified inactive Fock matrix
c   q --- scratch space for integrals
c   w --- W matrix, with W(pi,tu) = 2 (pi|tu) - (pt|ui) (p active/virt.)
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),q(*),w(*)
      dimension ioa0(*),ioa1(*),iob0(*),iob1(*),inter(na,nact,2)
      dimension iwa0(*),iwa1(*),iwb(*),nwa0(*),nwa1(*),nwb(*)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
      dimension mk1(8),mk2(8),mk3(8),mk4(8),mk5(8)

      ind(i,j) = ((i-1)*(i-2))/2 + j

      time = cpulft(1)
      npiece = (norb-1)/lenrs + 1
      write (iwr,1) lenrs,npiece
1     format (//
     +     'maximal length of (rs) pieces in routine COUPLE: ',i5/
     +     'number of pieces:                                ',i5//)
      i1 = icorr(nbasq)
      i2 = icorr(nbasq)
      ibco = iblo7
      write (iwr,2) ibco
2     format ('starting storage of A(co),B(co) at block no. ',i7//)
      k = 0
      do 10 isyma=1,nirr
      mk3(isyma) = k
10    k = k + nsymm(isyma)**2

c   read in inactive fock matrix
 
      call rdedx (f,nfock,iblofi,ned7)

c   read in w matrix

      irow = norbti + norbai
      icol = npair(isymr)
      isize = irow * icol
      if (isize .gt. 0) call rdedx (w,isize,iblkw,ned7)

c   read in integrals (ti|vx),(au|vx) and construct modified fock matrices

      call fmove (f,g,nfock)
      intpos = 0
      intfil = ned0
      iblf = iblkgv
      intmod = 0
      do 100 iv=1,nact
      isymv = itypea(iv)
      do 100 ix=1,iv
      isymx = itypea(ix)
      isymvx = mult(isymv,isymx)
      j1 = i1
      do 20 isym1=1,nirr
      isym2 = mult(isym1,isymvx)
      if (isym2 .gt. isym1) go to 20
      mc1 = mcore(isym1)
      ma1 = mact (isym1)
      mv1 = mvirt(isym1)
      mc2 = mcore(isym2)
      ma2 = mact (isym2)
      mv2 = mvirt(isym2)
      m = ma1 * (mc2+mv2) + ma2 * (mc1+mv1)
      if (m .eq. 0) go to 20
      call intin (q(j1),m)
      mk1(isym1) = j1
      mk2(isym1) = j1 + ma1*mc2 + ma2*mc1
      j1 = j1 + m
      if (isym2 .eq. isym1) go to 20
      mk1(isym2) = mk1(isym1) + ma1*mc2
      mk2(isym2) = mk2(isym1) + ma2*mv1
      call trnsps (q(mk1(isym2)),q(i2),mc1,ma2)
      call fmove (q(i2),q(mk1(isym2)),mc1*ma2)
      call trnsps (q(mk2(isym2)),q(i2),ma1,mv2)
      call fmove (q(i2),q(mk2(isym2)),ma1*mv2)
20    continue

c   FI(ix) - (vi|vx)  for all i

      ioffx = mk3(isymx) + (npoint(ix+ncore)+mcore(isymx)-1)
     +                   *nsymm(isymx)
      iofqv = mk1(isymv)-1 + npoint(iv+ncore) - mact(isymv)
      do 30 ii=1,mcore(isymx)
30    g(ioffx+ii) = g(ioffx+ii) - q(iofqv+ii*mact(isymv))

c   FI(ax) - (av|vx)  for all a

      ioffx = ioffx + nprm(isymx)
      iofqv = mk2(isymx)-1 + (npoint(iv+ncore)-1)*mvirt(isymx)
      do 40 ia=1,mvirt(isymx)
40    g(ioffx+ia) = g(ioffx+ia) - q(iofqv+ia)

      if (ix .eq. iv) go to 70

c   FI(iv) - (xi|vx)  for all i, if  iv.ne.ix

      ioffv = mk3(isymv) + (npoint(iv+ncore)+mcore(isymv)-1)
     +                   *nsymm(isymv)
      iofqx = mk1(isymx)-1 + npoint(ix+ncore) - mact(isymx)
      do 50 ii=1,mcore(isymv)
50    g(ioffv+ii) = g(ioffv+ii) - q(iofqx+ii*mact(isymx))

c   FI(av) - (ax|vx)  for all a, if  iv.ne.ix

      ioffv = ioffv + nprm(isymv)
      iofqx = mk2(isymv)-1 + (npoint(ix+ncore)-1)*mvirt(isymv)
      do 60 ia=1,mvirt(isymv)
60    g(ioffv+ia) = g(ioffv+ia) - q(iofqx+ia)

70    continue

100   continue

      iuaa = 0
      iubb = 0
      icount = 0
      do 1600 isyma=1,nirr
      if (nstra(isyma) .eq. 0) go to 1600
      nalock = (nstra(isyma)-1)/nablk + 1
      do 1590 isymb=1,nirr
      if (nstrb(isymb) .eq. 0) go to 1590
      icount = icount + 1
      nblock = (nstrb(isymb)-1)/nbblk + 1

c   symmetry of intermediate states

      isymk = mult(isyma,isymb)

c   symmetry of orbital excitations

      isymr0 = mult(isss,isymk)

c   symmetry of excited beta strings

      isymb0 = mult(isymb,isymr0)
      nstrb0 = nstrb(isymb0)
      mtaa = 0

c   loop over alpha strings in block

      do 1580 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)

c   generate alpha replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isss

      mta = mtaa
      do 120 i=1,na
120   icga(i) = icgaa(i)
      iua = iuaa
      iofiwa = 1
      do 140 ia=1,naa
130   call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      if (jmta .ne. isyma) go to 130
      if (ia .eq. 1) ia0 = ioa0(mta)
      call rplace (icga,isymr0,iwa0(iofiwa),nwa0(ia),na,ioa0,inter)
140   iofiwa = iofiwa + maxrp3

c   generate alpha replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isym

      isymr1 = mult(isym,isymk)
      isymb1 = mult(isymb,isymr1)
      nstrb1 = nstrb(isymb1)
      if (isym .eq. isss) then
         call icopy (nablk,nwa0,1,nwa1,1)
         call icopy (maxrp3*nablk,iwa0,1,iwa1,1)
         ia1 = ia0
      else
         mta = mtaa
         do 150 i=1,na
150      icga(i) = icgaa(i)
         iua = iuaa
         iofiwa = 1
         do 170 ia=1,naa
160      call string (mta,nact,na,icga,itypea,jmta,iua,mult)
         if (jmta .ne. isyma) go to 160
         if (ia .eq. 1) ia1 = ioa1(mta)
         call rplace (icga,isymr1,iwa1(iofiwa),nwa1(ia),na,ioa1,inter)
170      iofiwa = iofiwa + maxrp3
      end if

      mtbb = 0

c   loop over beta strings in block

      do 1560 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa * nbb
      nabp0 = naabb * npair(isymr0)

c   generate beta replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isss

      mtb = mtbb  
      iub = iubb
      do 180 i=1,nb
180   icgb(i) = icgbb(i)
      iofiwb = 1
      do 200 ib=1,nbb
190   call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      if (jmtb .ne. isymb) go to 190
      if (ib .eq. 1) ib0 = iob0(mtb) - 1
      call rplace (icgb,isymr0,iwb(iofiwb),nwb(ib),nb,iob0,inter(1,1,2))
200   iofiwb = iofiwb + maxrp3

c   construct d matrix

      call vclr(d,1,nabp0)
      iofiw1 = 0
      iofca = ib0
      iofd = - naabb
      do 240 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 230 iw=1,nwa0(ia)
      ioffc = iofca + iwa0(iofiw1+iw)
      ioffd = iofd + ipos(iwa0(iofiw3+iw)) * naabb
      if (iwa0(iofiw2+iw) .gt. 0) then
         do 210 ib=1,nbb
210      d(ioffd+ib) = d(ioffd+ib) + c(ioffc+ib)
      else
         do 220 ib=1,nbb
220      d(ioffd+ib) = d(ioffd+ib) - c(ioffc+ib)
      end if
230   continue
      iofiw1 = iofiw1 + maxrp3
240   iofd = iofd + nbb
      iofd = 1 - nbb - naabb
      iofcb = ia0 - nstrb0
      iofiw1 = 0
      do 280 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 270 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + ipos(iwb(iofiw3+iw)) * naabb
      if (iwb(iofiw2+iw) .gt. 0) then
         do 250 ia=1,naa
250      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c(ioffc+ia*nstrb0)
      else
         do 260 ia=1,naa
260      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c(ioffc+ia*nstrb0)
      end if
270   continue
      iofiw1 = iofiw1 + maxrp3
280   iofd = iofd + 1
      if (odebug(33)) then
         write (iwr,951) isyma,isymb,ialock,iblock
951      format ('d matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 955 iv=1,nact
         do 955 ix=1,nact
         if (mult(itypea(iv),itypea(ix)) .eq. isymr0) then
            write (iwr,952) iv,ix
952         format ('pair (',i2,',',i2,'):')
            k = (ipos(ilifa(iv)+ix)-1)*naabb
            do 953 ib=1,nbb
953         write (iwr,954) (d(k+(ia-1)*nbb+ib),ia=1,naa)
954         format (10f9.4)
         end if
955      continue
      end if

c   generate beta replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isym

      if (isym .eq. isss) then
         ib1 = ib0
      else
         mtb = mtbb
         iub = iubb
         do 290 i=1,nb
290      icgb(i) = icgbb(i)
         iofiwb = 1
         do 310 ib=1,nbb
300      call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
         if (jmtb .ne. isymb) go to 300
         if (ib .eq. 1) ib1 = iob1(mtb) - 1
         call rplace (icgb,isymr1,iwb(iofiwb),nwb(ib),nb,iob1,
     +                inter(1,1,2))
310      iofiwb = iofiwb + maxrp3
      end if

c   loop over (rs) pieces of matrices A(co),B(co)

      nrs = 0
      lrs = lenrs
      iblka = ibco
      do 990 ipiece=1,npiece
      if (ipiece .eq. npiece) lrs = norb - nrs
      lenco = lrs * nci1
      iblkb = iblka + lensec(lenco)

c   if first block, initiate A(co) and B(co) parts to zero,
c   else read A(co) and B(co) parts from external file

      if (icount.eq.1 .and. ialock.eq.1 .and. iblock.eq.1) then
         call vclr(a,1,lenco)
         call vclr(b,1,lenco)
c
c   add c contributions
c
         if (isym .eq. isss) then
            do 330 ip=nst,nba
            isymp = itype(ip)
            do 320 ii=1,ncore
            if (itype(ii) .eq. isymp) then
               k = ipo(ind(ip,ii))
               if (k.gt.nrs .and. k.le.nrs+lrs) then
                  if (ip .le. nprim) then
                     iadd = mcore(isymp)
                  else
                     iadd = nprm(isymp)
                  end if
                  fock = f(mk3(isymp)+(npoint(ip)+iadd-1)*nsymm(isymp)
     ?                    +npoint(ii))
                  fock = fock + fock
                  call daxpy(nci0, fock,c,1,a((k-nrs-1)*nci0+1),1)
                  call daxpy(nci0,-fock,c,1,b((k-nrs-1)*nci0+1),1)
               end if
            end if
320         continue
330         continue
         end if
      else
         call rdedx (a,lenco,iblka,ned7)
         call rdedx (b,lenco,iblkb,ned7)
      end if

c   d contributions that do not depend on integrals

      if (isymk .ne. isym) go to 460

c   fock matrix contributions

      do 380 it=nst,nprim
      isymt = itype(it)
      isymi = mult(isymt,isymr)
      do 370 ii=1,ncore
      if (itype(ii) .eq. isymi) then
         k = ipo(ind(it,ii))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 360 ix=nst,nprim
            if (itype(ix) .ne. isymi) go to 360
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
            iofda = (ipos(ilifa(it-ncore)+ix-ncore)-1)*naabb
            iofdb = (ipos(ilifa(ix-ncore)+it-ncore)-1)*naabb
            if = mk3(isymi) + (npoint(ix)+mcore(isymi)-1)*nsymm(isymi) 
     ?                      + npoint(ii)
            focka = g(if)
            fockb = f(if)
            do 350 ia=1,naa
            do 340 ib=1,nbb
            a(iofab+ib) = a(iofab+ib) - focka * d(iofda+ib)
340         b(iofab+ib) = b(iofab+ib) + fockb * d(iofdb+ib)
            iofda = iofda + nbb
            iofdb = iofdb + nbb
350         iofab = iofab + nstrb(isymb)
360         continue
         end if
      end if
370   continue
380   continue

      do 430 ic=nprimp,nba
      isymc = itype(ic)
      isymt = mult(isymc,isymr)
      do 420 it=nst,nprim
      if (itype(it) .eq. isymt) then
         k = ipo(ind(ic,it))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 410 ix=nst,nprim
            if (itype(ix) .ne. isymc) go to 410
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
c           iofhe = ia1 + ib1
            iofda = (ipos(ilifa(ix-ncore)+it-ncore)-1)*naabb
            iofdb = (ipos(ilifa(it-ncore)+ix-ncore)-1)*naabb
            if = mk3(isymc) + (npoint(ix)+mcore(isymc)-1)*nsymm(isymc) 
     ?                      + npoint(ic)+nprm(isymc)
            focka = f(if)
            fockb = g(if)
            do 400 ia=1,naa
            do 390 ib=1,nbb
            a(iofab+ib) = a(iofab+ib) + focka * d(iofda+ib)
390         b(iofab+ib) = b(iofab+ib) - fockb * d(iofdb+ib)
            iofda = iofda + nbb
            iofdb = iofdb + nbb
400         iofab = iofab + nstrb(isymb)
410         continue
         end if
      end if
420   continue
430   continue

c   w matrix contributions

      if (isize .eq. 0) go to 460
      lenw = min(lrs,irow-nrs)
      if (lenw .le. 0) go to 460
      do 450 iv=1,nact
      do 450 ix=1,nact
      if (mult(itypea(iv),itypea(ix)) .eq. isymr) then
         ivx = ipos(ilifa(iv)+ix)
         ixv = ipos(ilifa(ix)+iv)
         iofda0 = (ivx-1) * naabb
         iofdb0 = (ixv-1) * naabb
         iofab0 = ia1 + ib1
         ioffw = (ivx-1)*irow + nrs
         do 440 irs=1,lenw
         iofab = iofab0
         iofda = iofda0
         iofdb = iofdb0
         ww = w(ioffw+irs)
         do 435 ia=1,naa
         do 433 ib=1,nbb
         a(iofab+ib) = a(iofab+ib) + ww * d(iofda+ib)
433      b(iofab+ib) = b(iofab+ib) - ww * d(iofdb+ib)
         iofab = iofab + nstrb(isymb)
         iofda = iofda + nbb
435      iofdb = iofdb + nbb
440      iofab0 = iofab0 + nci1
      end if
450   continue

c   read in integrals (ti|vx),(at|vx)

460   np = 0
      do 470 isym1=1,nirr
      isym2 = mult(isym1,isymr0)
      mk4(isym1) = np
      np = np + mact(isym1) * mcore(isym2)
      mk5(isym1) = np
470   np = np + mvirt(isym1) * mact(isym2)
      nabp1 = naabb * np
      call vclr(e,1,nabp1)

      intpos = 0
      intfil = ned0
      iblf = iblkgv
      intmod = 0
      do 800 iv=1,nact
      isymv = itypea(iv)
      do 800 ix=1,iv
      isymx = itypea(ix)
      isymvx = mult(isymv,isymx)
      j1 = i1
      do 480 isym1=1,nirr
      isym2 = mult(isym1,isymvx)
      if (isym2 .gt. isym1) go to 480
      mc1 = mcore(isym1)
      ma1 = mact (isym1)
      mv1 = mvirt(isym1)
      mc2 = mcore(isym2)
      ma2 = mact (isym2)
      mv2 = mvirt(isym2)
      m = ma1 * (mc2+mv2) + ma2 * (mc1+mv1)
      if (m .eq. 0) go to 480
      call intin (q(j1),m)
      mk1(isym1) = j1
      mk2(isym1) = j1 + ma1*mc2 + ma2*mc1
      j1 = j1 + m
      if (isym2 .eq. isym1) go to 480
      mk1(isym2) = mk1(isym1) + ma1*mc2
      mk2(isym2) = mk2(isym1) + ma2*mv1
      call trnsps (q(mk1(isym2)),q(i2),mc1,ma2)
      call fmove (q(i2),q(mk1(isym2)),mc1*ma2)
      call trnsps (q(mk2(isym2)),q(i2),ma1,mv2)
      call fmove (q(i2),q(mk2(isym2)),ma1*mv2)
480   continue
      if (isymvx .ne. isymr0) go to 500

c   contribution to e matrix

      ivx = (ipos(ilifa(iv)+ix)-1) * naabb
      ixv = (ipos(ilifa(ix)+iv)-1) * naabb
      do 490 isym1=1,nirr
      isym2 = mult(isym1,isymr0)
      ma1mc2 = mact(isym1) * mcore(isym2)
      mv1ma2 = mvirt(isym1) * mact(isym2)
      if (odebug(42)) then
         write (iwr,481) isyma,isymb,ialock,iblock,
     +                   iv+ncore,ix+ncore,isym1,isym2
481      format ('symk = ',i1,' x ',i1,' - blocks ',i2,',',i2,
     ?           ' - orbital pair ',2i2,' - vectors of symm ',2i2)
         call outvec (q(mk1(isym1)),ma1mc2,'ti integrals')
         call outvec (q(mk2(isym1)),mv1ma2,'at integrals')
         call outsqr (d(ivx),naa,naa,nbb,'d(vx) matrix')
      end if
      iofeti = mk4(isym1) * naabb
      iofeat = mk5(isym1) * naabb
      iofzti = mk1(isym1) - 1
      iofzat = mk2(isym1) - 1
      do 483 irs=1,ma1mc2
      do 482 i=1,naabb
482   e(iofeti+i) = e(iofeti+i) + d(ivx+i) * q(iofzti+irs)
483   iofeti = iofeti + naabb
      do 485 irs=1,mv1ma2
      do 484 i=1,naabb
484   e(iofeat+i) = e(iofeat+i) + d(ivx+i) * q(iofzat+irs)
485   iofeat = iofeat + naabb
c         call mxmb (d(ivx),1,0, q(mk1(isym1)),0,1,
c     ?              e(mk4(isym1)*naabb+1),1,naabb, naabb,1,ma1mc2)
c         call mxmb (d(ivx),1,0, q(mk2(isym1)),0,1,
c     ?              e(mk5(isym1)*naabb+1),1,naabb, naabb,1,mv1ma2)
      if (ix .eq. iv) go to 490
      if (odebug(42)) call outsqr (d(ixv),naa,naa,nbb,'d(xv) matrix')
      iofeti = mk4(isym1) * naabb
      iofeat = mk5(isym1) * naabb
      do 487 irs=1,ma1mc2
      do 486 i=1,naabb
486   e(iofeti+i) = e(iofeti+i) + d(ixv+i) * q(iofzti+irs)
487   iofeti = iofeti + naabb
      do 489 irs=1,mv1ma2
      do 488 i=1,naabb
488   e(iofeat+i) = e(iofeat+i) + d(ixv+i) * q(iofzat+irs)
489   iofeat = iofeat + naabb
c         call mxmb (d(ixv),1,0, q(mk1(isym1)),0,1,
c     ?              e(mk4(isym1)*naabb+1),1,naabb, naabb,1,ma1mc2)
c         call mxmb (d(ixv),1,0, q(mk2(isym1)),0,1,
c     ?              e(mk5(isym1)*naabb+1),1,naabb, naabb,1,mv1ma2)
490   continue

c   d contributions that depend on integrals

500   continue
      if (isymk .ne. isym) go to 800
      isymav = mult(isymr,isymv)
      isymax = mult(isymr,isymx)

      do 600 ic=nprimp,nba

c   A(co)(m,av) <-- A(co)(m,av) - (at|vx) * D(m,tx)

      if (itype(ic) .eq. isymav) then
         k = ipo(ind(ic,iv+ncore))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 530 it=1,nact
            isymt = itypea(it)
            if (mult(isymt,isymx) .ne. isymr0) go to 530
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
            iofda = (ipos(ilifa(it)+ix)-1)*naabb
            z = q(mk2(isymav)-1+(npoint(it+ncore)-1)*mvirt(isymav)
     +                         +npoint(ic))
            do 520 ia=1,naa
            do 510 ib=1,nbb
510         a(iofab+ib) = a(iofab+ib) - z * d(iofda+ib)
            iofda = iofda + nbb
520         iofab = iofab + nstrb(isymb)
530         continue
         end if
      end if
      if (iv .eq. ix) go to 600

c   A(co)(m,ax) <-- A(co)(m,ax) - (at|vx) * D(m,tv)   (if ix.ne.iv)

      if (itype(ic) .eq. isymax) then
         k = ipo(ind(ic,ix+ncore))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 560 it=1,nact
            isymt = itypea(it)
            if (mult(isymt,isymv) .ne. isymr0) go to 560
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
            iofda = (ipos(ilifa(it)+iv)-1)*naabb
            z = q(mk2(isymax)-1+(npoint(it+ncore)-1)*
     +                           mvirt(isymax)+npoint(ic))
            do 550 ia=1,naa
            do 540 ib=1,nbb
540         a(iofab+ib) = a(iofab+ib) - z * d(iofda+ib)
            iofda = iofda + nbb
550         iofab = iofab + nstrb(isymb)
560         continue
         end if
      end if
600   continue

      isymiv = mult(isymr,isymv)
      isymix = mult(isymr,isymx)

      do 700 ii=1,ncore

c   B(co)(m,vi) <-- B(co)(m,vi) - (ti|vx) * D(m,tx)

      if (itype(ii) .eq. isymiv) then
         k = ipo(ind(iv+ncore,ii))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 630 it=1,nact
            isymt = itypea(it)
            if (mult(isymt,isymx) .ne. isymr0) go to 630
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
            iofdb = (ipos(ilifa(it)+ix)-1)*naabb
            z = q(mk1(isymt)-1+(npoint(ii)-1)*mact(isymt)+
     +                          npoint(it+ncore))
            do 620 ia=1,naa
            do 610 ib=1,nbb
610         b(iofab+ib) = b(iofab+ib) - z * d(iofdb+ib)
            iofdb = iofdb + nbb
620         iofab = iofab + nstrb(isymb)
630         continue
         end if
      end if
      if (iv .eq. ix) go to 700

c   B(co)(m,xi) <-- B(co)(m,xi) - (ti|vx) * D(m,tv)   (if ix.ne.iv)

      if (itype(ii) .eq. isymix) then
         k = ipo(ind(ix+ncore,ii))
         if (k.gt.nrs .and. k.le.nrs+lrs) then
            do 660 it=1,nact
            isymt = itypea(it)
            if (mult(isymt,isymv) .ne. isymr0) go to 660
            iofab = (k-nrs-1)*nci1 + ia1 + ib1
            iofdb = (ipos(ilifa(it)+iv)-1)*naabb
            z = q(mk1(isymt)-1+(npoint(ii)-1)*mact(isymt)+
     +                          npoint(it+ncore))
            do 650 ia=1,naa
            do 640 ib=1,nbb
640         b(iofab+ib) = b(iofab+ib) - z * d(iofdb+ib)
            iofdb = iofdb + nbb
650         iofab = iofab + nstrb(isymb)
660         continue
         end if
      end if
700   continue
800   continue
      if (odebug(42)) then
         write (iwr,751) isyma,isymb,ialock,iblock
751      format ('e matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 755 ir=nst,nba
         isym1 = itype(ir)
         if (ir .le. nprim) then
            is0 = 1
            is1 = ncore
            k = mk4(isym1)
            m = mact(isym1)
         else
            is0 = nst
            is1 = nprim
            k = mk5(isym1)
            m = mvirt(isym1)
         end if
         do 755 is=is0,is1
         isym2 = itype(is)
         if (mult(isym1,isym2) .eq. isymr0) then
            iofe = (k+(npoint(is)-1)*m+npoint(ir)-1) * naabb
            write (iwr,752) ir,is
752         format ('pair (',i2,',',i2,'):')
            do 753 ib=1,nbb
753         write (iwr,754) (e(iofe+(ia-1)*nbb+ib),ia=1,naa)
754         format (100f9.4)
         end if
755      continue
      end if

c   e matrix contributions to a and b
      
      iofaba = ib1
      iofiw1 = 0
      iofe = 0
      do 880 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 870 iw=1,nwa1(ia)
      iut = iwa1(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      isymt = itypea(it)
      isymu = itypea(iu)
      iofmab = iofaba + iwa1(iofiw1+iw)
c      if (odebug(57)) write (iwr,801) iofmab,it,iu,ia,iwa1(iofiw2+iw)
c801   format ('alpha replacement: <',i2,'|E(',i2,',',i2,')|',i2,'> = ',i3)

      isymi = mult(isymt,isymr)
      do 830 ii=1,ncore
      if (itype(ii) .ne. isymi) go to 830
      k = ipo(ind(it+ncore,ii))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofa = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk4(isymu)+(npoint(ii)-1)*mact(isymu)+
     +                        npoint(iu+ncore)-1)* naabb + iofe

c   A(co)(m,ti) <-- A(co)(m,ti) - E(k,ui) * <m|E(alpha)_tu|k>

         if (iwa1(iofiw2+iw) .gt. 0) then
            do 810 ib=1,nbb
c            if (odebug(57)) write (iwr,802) it+ncore,ii,e(ioffe+ib)
c802         format ('a(',i2,',',i2,'): subtracting ',f9.4)
810         a(iofa+ib) = a(iofa+ib) - e(ioffe+ib)
         else
            do 820 ib=1,nbb
c            if (odebug(57)) write (iwr,803) it+ncore,ii,e(ioffe+ib)
c803         format ('a(',i2,',',i2,'): adding      ',f9.4)
820         a(iofa+ib) = a(iofa+ib) + e(ioffe+ib)
         end if
      end if
830   continue

      isymi = mult(isymu,isymr)
      do 860 ii=1,ncore
      if (itype(ii) .ne. isymi) go to 860
      k = ipo(ind(iu+ncore,ii))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofb = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk4(isymt)+(npoint(ii)-1)*mact(isymt)+
     +                        npoint(it+ncore)-1) * naabb + iofe

c   B(co)(m,ui) <-- B(co)(m,ui) + E(k,ti) * <m|E(alpha)_tu|k>

         if (iwa1(iofiw2+iw) .gt. 0) then
            do 840 ib=1,nbb
c            if (odebug(57)) write (iwr,1041) iu+ncore,ii,e(ioffe+ib)
840         b(iofb+ib) = b(iofb+ib) + e(ioffe+ib)
         else
            do 850 ib=1,nbb
c            if (odebug(57)) write (iwr,1031) iu+ncore,ii,e(ioffe+ib)
850         b(iofb+ib) = b(iofb+ib) - e(ioffe+ib)
         end if
      end if
860   continue

      isymc = mult(isymu,isymr)
      do 1030 ic=nprimp,nba
      if (itype(ic) .ne. isymc) go to 1030
      k = ipo(ind(ic,iu+ncore))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofa = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk5(isymc)+(npoint(it+ncore)-1)*mvirt(isymc)+
     +                        npoint(ic)-1) * naabb + iofe

c   A(co)(m,au) <-- A(co)(m,au) + E(k,at) * <m|E(alpha)_tu|k>

         if (iwa1(iofiw2+iw) .gt. 0) then
            do 1010 ib=1,nbb
c            if (odebug(57)) write (iwr,803) ic,iu+ncore,e(ioffe+ib)
1010        a(iofa+ib) = a(iofa+ib) + e(ioffe+ib)
         else
            do 1020 ib=1,nbb
c            if (odebug(57)) write (iwr,802) ic,iu+ncore,e(ioffe+ib)
1020        a(iofa+ib) = a(iofa+ib) - e(ioffe+ib)
         end if
      end if
1030  continue

      isymc = mult(isymt,isymr)
      do 1060 ic=nprimp,nba
      if (itype(ic) .ne. isymc) go to 1060
      k = ipo(ind(ic,it+ncore))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofb = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk5(isymc)+(npoint(iu+ncore)-1)*mvirt(isymc)+
     +                        npoint(ic)-1) * naabb + iofe

c   B(co)(m,at) <-- B(co)(m,at) - E(k,au) * <m|E(alpha)_tu|k>

         if (iwa1(iofiw2+iw) .gt. 0) then
            do 1040 ib=1,nbb
c            if (odebug(57)) write (iwr,1031) ic,it+ncore,e(ioffe+ib)
c1031        format ('b(',i2,',',i2,'): subtracting ',f9.4)
1040        b(iofb+ib) = b(iofb+ib) - e(ioffe+ib)
         else
            do 1050 ib=1,nbb
c            if (odebug(57)) write (iwr,1041) ic,it+ncore,e(ioffe+ib)
c1041        format ('b(',i2,',',i2,'): adding      ',f9.4)
1050        b(iofb+ib) = b(iofb+ib) + e(ioffe+ib)
         end if
      end if
1060  continue
870   continue
      iofiw1 = iofiw1 + maxrp3
880   iofe = iofe + nbb

      iofabb = ia1 - nstrb1
      if (odebug(57)) print*,'ia1 = ',ia1,'    nstrb1 = ',nstrb1
      iofiw1 = 0
      iofe = 1 - nbb
      do 980 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 970 iw=1,nwb(ib)
      iut = iwb(iofiw3+iw)
      iu = (iut-1)/nact + 1
      it = iut - ilifa(iu)
      isymt = itypea(it)
      isymu = itypea(iu)
      iofmab = iofabb + iwb(iofiw1+iw)
c      if (odebug(57)) write (iwr,901) iwb(iofiw1+iw),it,iu,ib,iwb(iofiw2+iw)
c901   format ('beta  replacement: <',i2,'|E(',i2,',',i2,')|',i2,'> = ',i3)

      isymi = mult(isymt,isymr)
      do 930 ii=1,ncore
      if (itype(ii) .ne. isymi) go to 930
      k = ipo(ind(it+ncore,ii))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofa = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk4(isymu)+(npoint(ii)-1)*mact(isymu)+
     +                        npoint(iu+ncore)-1) * naabb + iofe

c   A(co)(m,ti) <-- A(co)(m,ti) - E(k,ui) * <m|E(beta)_tu|k>

         if (iwb(iofiw2+iw) .gt. 0) then
            do 910 ia=1,naa
c            if (odebug(57)) write (iwr,802) it+ncore,ii,e(ioffe+ia*nbb)
910         a(iofa+ia*nstrb1) = a(iofa+ia*nstrb1) - e(ioffe+ia*nbb)
         else
            do 920 ia=1,naa
c            if (odebug(57)) write (iwr,803) it+ncore,ii,e(ioffe+ia*nbb)
920         a(iofa+ia*nstrb1) = a(iofa+ia*nstrb1) + e(ioffe+ia*nbb)
         end if
      end if
930   continue

      isymi = mult(isymu,isymr)
      do 960 ii=1,ncore
      if (itype(ii) .ne. isymi) go to 960
      k = ipo(ind(iu+ncore,ii))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofb = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk4(isymt)+(npoint(ii)-1)*mact(isymt)+
     +                        npoint(it+ncore)-1) * naabb + iofe

c   B(co)(m,ui) <-- B(co)(m,ui) + E(k,ti) * <m|E(beta)_tu|k>

         if (iwb(iofiw2+iw) .gt. 0) then
            do 940 ia=1,naa
c            if (odebug(57)) write (iwr,1041) iu+ncore,ii,e(ioffe+ia*nbb)
940         b(iofb+ia*nstrb1) = b(iofb+ia*nstrb1) + e(ioffe+ia*nbb)
         else
            do 950 ia=1,naa
c            if (odebug(57)) write (iwr,1031) iu+ncore,ii,e(ioffe+ia*nbb)
950         b(iofb+ia*nstrb1) = b(iofb+ia*nstrb1) - e(ioffe+ia*nbb)
         end if
      end if
960   continue

      isymc = mult(isymu,isymr)
      do 1130 ic=nprimp,nba
      if (itype(ic) .ne. isymc) go to 1130
      k = ipo(ind(ic,iu+ncore))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofa = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk5(isymc)+(npoint(it+ncore)-1)*mvirt(isymc)+
     +                        npoint(ic)-1) * naabb + iofe

c   A(co)(m,au) <-- A(co)(m,au) + E(k,at) * <m|E(beta)_tu|k>

         if (iwb(iofiw2+iw) .gt. 0) then
            do 1110 ia=1,naa
c            if (odebug(57)) write (iwr,803) ic,iu+ncore,e(ioffe+ia*nbb)
1110        a(iofa+ia*nstrb1) = a(iofa+ia*nstrb1) + e(ioffe+ia*nbb)
         else
            do 1120 ia=1,naa
c            if (odebug(57)) write (iwr,802) ic,iu+ncore,e(ioffe+ia*nbb)
1120        a(iofa+ia*nstrb1) = a(iofa+ia*nstrb1) - e(ioffe+ia*nbb)
         end if
      end if
1130  continue

      isymc = mult(isymt,isymr)
      do 1160 ic=nprimp,nba
      if (itype(ic) .ne. isymc) go to 1160
      k = ipo(ind(ic,it+ncore))
      if (k.gt.nrs .and. k.le.nrs+lrs) then
         iofb = (k-nrs-1)*nci1 + iofmab
         ioffe = (mk5(isymc)+(npoint(iu+ncore)-1)*mvirt(isymc)+
     +                        npoint(ic)-1) * naabb + iofe

c   B(co)(m,at) <-- B(co)(m,at) - E(k,au) * <m|E(beta)_tu|k>

         if (iwb(iofiw2+iw) .gt. 0) then
            do 1140 ia=1,naa
c            if (odebug(57)) write (iwr,1031) ic,it+ncore,e(ioffe+ia*nbb)
1140        b(iofb+ia*nstrb1) = b(iofb+ia*nstrb1) - e(ioffe+ia*nbb)
         else
            do 1150 ia=1,naa
c            if (odebug(57)) write (iwr,1041) ic,it+ncore,e(ioffe+ia*nbb)
1150        b(iofb+ia*nstrb1) = b(iofb+ia*nstrb1) + e(ioffe+ia*nbb)
         end if
      end if
1160  continue
970   continue
      iofiw1 = iofiw1 + maxrp3
980   iofe = iofe + 1

c   write updated (rs) piece to external file

      call wrt3 (a,lenco,iblka,ned7)
      call wrt3 (b,lenco,iblkb,ned7)
      iblka = iblka + 2 * lensec(lenco)

c   end of loop over (rs) pieces

990   nrs = nrs + lrs

c   end of beta block

      do 1550 i=1,nb
1550  icgbb(i) = icgb(i)
      iubb = iub
1560  mtbb = mtb

c   end of alpha block

      do 1570 i=1,na
1570  icgaa(i) = icga(i)
      iuaa = iua
1580  mtaa = mta

c   end of beta symmetry

1590  continue

c   end of alpha symmetry

1600  continue

      iblo7 = iblka

      if (odebug(34)) then
         nrs = 0
         lrs = lenrs
         iblka = ibco
         do 999 ipiece=1,npiece
         if (ipiece .eq. npiece) lrs = norb - nrs
         lenco = lrs * nci1
         iblkb = iblka + lensec(lenco)
         call rdedx (a,lenco,iblka,ned7)
         call rdedx (b,lenco,iblkb,ned7)
         do irs=1,lrs
         write (iwr,996) nrs+irs
996      format ('A(co) for orbital pair ',i4,':')
         call outvec (a((irs-1)*nci1+1),nci1,' ')
         write (iwr,997) nrs+irs
997      format ('B(co) for orbital pair ',i4,':')
         call outvec (b((irs-1)*nci1+1),nci1,' ')
         end do
         iblka = iblka + 2 * lensec(lenco)
999      nrs = nrs + lrs
      end if

      call corlsr (i1)
      write (iwr,1700) cpulft(1)-time
1700  format (//'construction of co part: ',f16.9,' seconds'//)

      return
      end
      subroutine denmat (c,g,p2,p3,ioa,iob,inter,iwa,iwb,d,
     +                   nablk,nbblk,nwa,nwb)
c
c---------------------------------------------------------------------
c   Constructs the one- and two-particle density matrices
c   and two permuted two-particle matrices p2, p3 which are needed
c   during the construction of the orbital-orbital part
c   of the matrices A and B
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      dimension c(*),g(*),d(*),nwa(*),nwb(*)
      dimension p2(*),p3(*)
      dimension ioa(*),iob(*),inter(na,nact,2),iwa(*),iwb(*)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
c
      call vclr(g,1,n0int)
      call vclr(p2,1,n0int2)
      call vclr(p3,1,n0int2)
c
      iuaa = 0
      iubb = 0
c
      do 220 isyma=1,nirr
      if (nstra(isyma).eq.0) goto 220
      nalock = (nstra(isyma)-1)/nablk+1
      do 210 isymb=1,nirr
      if (nstrb(isymb).eq.0) goto 210
      nblock = (nstrb(isymb)-1)/nbblk+1
c...symmetry of intermediate state
      isymk = mult(isyma,isymb)
c...symmetry of orbital excitations
      isymrs = mult(isss,isymk)
c...symmetry of excited beta strings
      isymbe = mult(isymb,isymrs)
      nstrbe = nstrb(isymbe)
      mtaa = 0
      do 200 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)
c...generate alpha replacements of required symmetry
      mta = mtaa
      do 10 i=1,na
10    icga(i) = icgaa(i)
      iua = iuaa
      iofiwa = 1
      do 30 ia=1,naa
20    call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      if (jmta.ne.isyma) goto 20
      if (ia.eq.1) ia1 = ioa(mta)
      call rplace(icga,isymrs,iwa(iofiwa),nwa(ia),na,ioa,inter)
30    iofiwa = iofiwa + maxrp3
c
      mtbb = 0
      do 180 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa * nbb
      naabbp = naabb * npair(isymrs)
c...loop over beta strings in block storing single replacements
      mtb = mtbb
      iub = iubb
      do 40 i=1,nb
40    icgb(i) = icgbb(i)
      iofiwb = 1
      do 60 ib=1,nbb
50    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      if (jmtb.ne.isymb) goto 50
      if (ib.eq.1) ib1 = iob(mtb) - 1
      call rplace(icgb,isymrs,iwb(iofiwb),nwb(ib),nb,iob,inter(1,1,2))
60    iofiwb = iofiwb + maxrp3
c
c...d matrix
      call vclr(d,1,naabbp)
      iofiw1 = 0
      iofca = ib1
      iofd = - naabb
      do 850 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 840 iw=1,nwa(ia)
      ioffc = iofca + iwa(iofiw1+iw)
      ioffd = iofd + ipos(iwa(iofiw3+iw))*naabb
      if (iwa(iofiw2+iw).gt.0) then
         do 810 ib=1,nbb
810      d(ioffd+ib) = d(ioffd+ib) + c(ioffc+ib)
      else
         do 830 ib=1,nbb
830      d(ioffd+ib) = d(ioffd+ib) - c(ioffc+ib)
      endif
840   continue
      iofiw1 = iofiw1 + maxrp3
850   iofd = iofd + nbb
      iofd = 1 - nbb - naabb
      iofcb = ia1 - nstrbe
      iofiw1 = 0
      do 950 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 940 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + ipos(iwb(iofiw3+iw))*naabb
      if (iwb(iofiw2+iw).gt.0) then
         do 910 ia=1,naa
910      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c(ioffc+ia*nstrbe)
      else
         do 930 ia=1,naa
930      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c(ioffc+ia*nstrbe)
      endif
940   continue
      iofiw1 = iofiw1 + maxrp3
950   iofd = iofd + 1
      if (odebug(22)) then
         write (iwr,951) isyma,isymb,ialock,iblock
951      format ('d matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 955 iv=1,nact
         ix1=nact
         do 955 ix=1,ix1
         if (mult(itypea(iv),itypea(ix)).eq.isymrs) then
            write (iwr,952) iv,ix
952         format ('pair (',i2,',',i2,'):')
            k = (ipos(ilifa(iv)+ix)-1)*naabb
            do 953 ib=1,nbb
953         write (iwr,954) (d(k+(ia-1)*nbb+ib),ia=1,naa)
954         format (10f9.4)
         endif
955      continue
      endif
c...one and two particle density matrix
      do 530 it=1,nact
      do 530 iu=1,nact
      if (mult(itypea(it),itypea(iu)).eq.isymrs) then
         iofut=(ipos(ilifa(iu)+it)-1)*naabb
         do 520 iv=1,nact
         do 520 ix=1,nact
         if(mult(itypea(iv),itypea(ix)).eq.isymrs) then
            ituvx=ipos(ilifa(it)+iu)+ioff0(ilifa(iv)+ix)
            iuvtx=ipos(ilifa(iu)+iv)+ioff0(ilifa(it)+ix)
            iuxtv=ipos(ilifa(iu)+ix)+ioff0(ilifa(it)+iv)
            iofvx=(ipos(ilifa(iv)+ix)-1)*naabb
            t=0.0d0
            do 510 i=1,naabb
510         t=t+d(iofut+i)*d(iofvx+i)
            g(ituvx)=g(ituvx)+t 
            p2(iuvtx)=p2(iuvtx)+t
            p3(iuxtv)=p3(iuxtv)+t
         endif
520      continue
      endif
530   continue
      if (isymrs.ne.1) goto 150
      iof1d = 0
      do 90 ip=n0int2+1,n0int
      iof1c = ia1 + ib1
      do 80 ia=1,naa
      do 70 ib=1,nbb
70    g(ip) = g(ip) + d(iof1d+ib) * c(iof1c+ib)
      iof1d = iof1d + nbb
80    iof1c = iof1c + nstrbe
90    continue
150   continue
c...end of beta block
      do 170 i=1,nb
170   icgbb(i) = icgb(i)
      iubb = iub
180   mtbb = mtb
c...end of alpha block
      do 190 i=1,na
190   icgaa(i) = icga(i)
      iuaa = iua
200   mtaa = mta
c...end of beta symmetry
210   continue
c...end of alpha symmetry
220   continue
c      call addone(g,g(n0int2+1),ioff0,1)
c...for vector computers: no recurrence in this loop
      do 310 itux=1,nact3
      it=(itux-1)/nact2+1
      iux=itux-(it-1)*nact2
      iu=(iux-1)/nact+1
      ix=iux-ilifa(iu)
      if(itypea(it).eq.itypea(ix))then
         ituux=ipos(ilifa(it)+iu)+ioff0(iux)
         iuutx=ipos(ilifa(iu)+iu)+ioff0(ilifa(it)+ix)
         iutxu=ipos(ilifa(iu)+it)+ioff0(ilifa(ix)+iu)
         dtx=g(n0int2+ipos(ilifa(it)+ix))
         g(ituux)=g(ituux)-dtx
         p2(iuutx)=p2(iuutx)-dtx
         p3(iutxu)=p3(iutxu)-dtx
      endif
310   continue
      return
      end
      subroutine dia (a1,a2,a,u,diag,q)
c---------------------------------------------------------------------
c   Transforms the diagonal of A(oo) to the pseudocanonical/natural
c   orbital basis. Ought to be rewritten but is fast enough.
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension a1(*),a2(*),a(*),u(*),diag(*),q(*)
      time = cpulft(1)
      call rdedx (u,ntrans,ibcan,ned7)
      if (odebug(94)) call outvec (u,ntrans,'u in dia:')
      call vclr(diag,1,norb)
      iblo = ibloab
      isym1 = 1
c mrow = # rows of actual symmetry isym1
      mrow = 0                  
c          that have been read in previous runs
      iofa = 0                  
      ipair = 0
      icase = 1
      np = lp
      ntot = (norb-1)/lp + 1
      do 100 ibox=1,ntot
      if (ibox .eq. ntot) np = norb - ipair
      mp = np          
c   mp = # rows available that have not been treated yet
      kpair = 0        
c kpair = # rows of symmetries < isym1
      npp = np * norb  
c           that have been treated in the current run
      if (odebug(94)) print*,'npp = ',npp
      call rdedx_less(q,npp,iblo,ned7)
10    isym2 = mult(isym1,isymr)
      m1 = mvirt(isym1)
      if (icase .eq. 1) m1 = mact(isym1)
      m2 = mcore(isym2)
      if (icase .eq. 3) m2 = mact(isym2)
      if (icase .eq. 1) iofset = iorbti(isym1)
      if (icase .eq. 2) iofset = iorbai(isym1)
      if (icase .eq. 3) iofset = iorbat(isym1)
      if (odebug(94)) print*,'iofset = ',iofset
      m = m1 * m2
      if (m .eq. 0) go to 30
      mrest = m - mrow
      nrow = min(mp,mrest)
      if (odebug(94)) print*,'m*nrow = ',m*nrow
      do 20 ipqrs=1,m*nrow
      ipqm1 = (ipqrs-1)/m
      irs = ipqrs - ipqm1 * m
      if (odebug(94)) print*,'kpair+ipqm1 = ',kpair+ipqm1
20    a(iofa+ipqrs) = q((kpair+ipqm1)*norb+iofset+irs)
      iofa = iofa + m * nrow
      mrow = mrow + nrow
      if (mp .lt. mrest) go to 90
      mp = mp - nrow
      kpair = kpair + nrow
30    isym1 = isym1 + 1
      if (isym1 .gt. nirr) go to 40
      mrow = 0
      if (mp .gt. 0) go to 10
      go to 90
40    if (icase .eq. 3)go to 90
      icase = icase + 1
      isym1 = 1
      mrow = 0
      if (mp .gt. 0) go to 10
90    ipair = ipair + np
100   iblo = iblo + lensec(2*npp)
      if (odebug(94)) call outvec (a,iofa,'vector a:')

c   transformation

c   t-i

      iofa=0
      do 190 isymac=1,nirr
      isymco = mult(isymac,isymr)
      ma = mact(isymac)
      mc = mcore(isymco)
      m = ma * mc
      if (m .eq. 0) go to 190
      call vclr(a2,1,m*m)
      do 110 ik=1,mc
      do 110 ijuvi=1,m*m
      iju = (ijuvi-1)/m+1
      ivi = ijuvi - (iju-1)*m
      iv = (ivi-1)/mc + 1
      ii = ivi - (iv-1)*mc
c      print*,'a2(',iju,',',ivi,'): adding ',a(iofa+ijuvi-ii+ik),'*',u(iofcor(isymco)+(ii-1)*mc+ik)
      a2(ijuvi) = a2(ijuvi) + a(iofa+ijuvi-ii+ik)*
     +                        u(iofcor(isymco)+(ii-1)*mc+ik)
110   continue
      if (odebug(94)) call outvec (a2,m*m,
     ?    'vector a2 (1st step of ti transformation):')
      call vclr(a1,1,m*ma)
      do 120 ij=1,mc
      do 120 iuvi=1,m*ma
      iu = (iuvi-1)/m + 1
      ivi = iuvi - (iu-1)*m
      iv = (ivi-1)/mc + 1
      ii = ivi-(iv-1)*mc
      a1(iuvi) = a1(iuvi) + a2(((iu-1)*mc+ij-1)*m+ivi)
     ?                    * u(iofcor(isymco)+(ii-1)*mc+ij)
120   continue
      if (odebug(94)) call outvec (a1,m*ma,
     ?    'vector a1 (2nd step of ti transformation):')
      call vclr(a2,1,m*ma)
      do 130 iv=1,ma
      do 130 iuti=1,m*ma
      iu=(iuti-1)/m+1
      iti=iuti-(iu-1)*m
      it=(iti-1)/mc+1
      ii=iti-(it-1)*mc
c      write(iwr,119)ii,iu,it,ii,a1(iuti+(iv-it)*mc),u(iofact(isymac)+(it-1)*ma+iv)
c119   format('a2(',4i2,'): adding ',f9.4,' * ',f9.4)
      a2(iuti)=a2(iuti)+a1(iuti+(iv-it)*mc)*
     +                   u(iofact(isymac)+(it-1)*ma+iv)
130   continue
      if (odebug(94)) call outvec (a2,m*ma,
     ?    'vector a2 (3rd step of ti transformation):')
      iofd=iorbti(isymac)
      if (odebug(94)) print*,'iofd = ',iofd
      do 140 iu=1,ma
      do 140 iti=1,m
      itmin1=(iti-1)/mc
      diag(iofd+iti)=diag(iofd+iti)+a2((iu-1)*m+iti) 
     ?                             *u(iofact(isymac)+itmin1*ma+iu)
140   continue
190   iofa=iofa+m*m
c...a-i
      do 290 isymvi=1,nirr
      isymco=mult(isymvi,isymr)
      mv=mvirt(isymvi)
      mc=mcore(isymco)
      m=mv*mc
      if(m.eq.0) go to 290
      call vclr(a2,1,m*m)
      do 210 ik=1,mc
      do 210 ijbci=1,m*m
      ijb=(ijbci-1)/m+1
      ici=ijbci-(ijb-1)*m
      ic=(ici-1)/mc+1
      ii=ici-(ic-1)*mc
      a2(ijbci)=a2(ijbci)+a(iofa+ijbci-ii+ik)*
     +                    u(iofcor(isymco)+(ii-1)*mc+ik)
210   continue
      if (odebug(94)) call outvec (a2,m*m,
     ?    'vector a2 (1st step of ai transformation):')
      call vclr(a1,1,m*mv)
      do 220 ij=1,mc
      do 220 ibci=1,m*mv
      ib=(ibci-1)/m+1
      ici=ibci-(ib-1)*m
      ic=(ici-1)/mc+1
      ii=ici-(ic-1)*mc
      a1(ibci)=a1(ibci)+a2(((ib-1)*mc+ij-1)*m+ici)
     ?                 *u(iofcor(isymco)+(ii-1)*mc+ij)
220   continue
      if (odebug(94)) call outvec (a1,m*mv,
     ?    'vector a1 (2nd step of ai transformation):')
      call vclr(a2,1,m*mv)
      do 230 ic=1,mv
      do 230 ibai=1,m*mv
      ib=(ibai-1)/m+1
      iai=ibai-(ib-1)*m
      ia=(iai-1)/mc+1
      ii=iai-(ia-1)*mc
c      write(iwr,119)ii,iu,it,ii,a1(iuti+(iv-it)*mc),u(iofact(isymac)+(it-1)*ma+iv)
c119   format('a2(',4i2,'): adding ',f9.4,' * ',f9.4)
      a2(ibai)=a2(ibai)+a1(ibai+(ic-ia)*mc)*
     +                  u(iofvir(isymvi)+(ia-1)*mv+ic)
230   continue
      if (odebug(94)) call outvec (a2,m*ma,
     ?    'vector a2 (3rd step of ai transformation):')
      iofd=iorbai(isymvi)
      if (odebug(94)) print*,'iofd = ',iofd
      do 240 ib=1,mv
      do 240 iai=1,m
      iamin1=(iai-1)/mc
      diag(iofd+iai)=diag(iofd+iai)+a2((ib-1)*m+iai) 
     ?                             *u(iofvir(isymvi)+iamin1*mv+ib)
240   continue
290   iofa=iofa+m*m
c...a-t
      do 390 isymvi=1,nirr
      isymac = mult(isymvi,isymr)
      mv = mvirt(isymvi)
      ma = mact(isymac)
      m = mv * ma
      if (m .eq. 0) go to 390
      call vclr(a2,1,m*m)
      do 310 iv=1,ma
      do 310 iubct=1,m*m
      iub=(iubct-1)/m+1
      ict=iubct-(iub-1)*m
      ic=(ict-1)/ma+1
      it=ict-(ic-1)*ma
      a2(iubct)=a2(iubct)+a(iofa+iubct-it+iv)*
     +                    u(iofact(isymac)+(it-1)*ma+iv)
310   continue
      if (odebug(94)) call outvec (a2,m*m,
     ?    'vector a2 (1st step of at transformation):')
      call vclr(a1,1,m*mv)
      do 320 iu=1,ma
      do 320 ibct=1,m*mv
      ib=(ibct-1)/m+1
      ict=ibct-(ib-1)*m
      ic=(ict-1)/ma+1
      it=ict-(ic-1)*ma
      a1(ibct)=a1(ibct)+a2(((ib-1)*ma+iu-1)*m+ict)
     ?                 *u(iofact(isymac)+(it-1)*ma+iu)
320   continue
      if (odebug(94)) call outvec (a1,m*mv,
     ?    'vector a1 (2nd step of at transformation):')
      call vclr(a2,1,m*mv)
      do 330 ic=1,mv
      do 330 ibat=1,m*mv
      ib=(ibat-1)/m+1
      iat=ibat-(ib-1)*m
      ia=(iat-1)/ma+1
      ii=iat-(ia-1)*ma
c      write(iwr,119)ii,iu,it,ii,a1(iuti+(iv-it)*mc),u(iofact(isymac)+(it-1)*ma+iv)
c119   format('a2(',4i2,'): adding ',f9.4,' * ',f9.4)
      a2(ibat)=a2(ibat)+a1(ibat+(ic-ia)*ma)*
     +                  u(iofvir(isymvi)+(ia-1)*mv+ic)
330   continue
      if (odebug(94)) call outvec (a2,m*ma,
     ?     'vector a2 (3rd step of at transformation):')
      iofd=iorbat(isymvi)
      if (odebug(94)) print*,'iofd = ',iofd
      do 340 ib=1,mv
      do 340 iat=1,m
      iamin1 = (iat-1)/ma
      diag(iofd+iat) = diag(iofd+iat) + a2((ib-1)*m+iat) 
     ?                                * u(iofvir(isymvi)+iamin1*mv+ib)
340   continue
390   iofa = iofa + m * m

      write (iwr,500) cpulft(1)-time
500   format ('transformation of diagonal ',27('.'),f9.4,' seconds'/)

      return
      end
      subroutine dipcon (c,d,dipx,dipy,dipz,vmyx,vmyy,vmyz,
     ?           ioa,iob,inter,iwa,iwb,nwa,nwb,nablk,nbblk)
c
c---------------------------------------------------------------------
c   Constructs the elements <0|my(p)|n>, p=x,y,z, which enter the
c   vector <0|[Gamma_j,my(p)]|0> needed for the calculation of 
c   oscillator strengths
c   (c) Carsten Fuchs 1991-1992
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      dimension c(*),d(*)
      dimension dipx(*),dipy(*),dipz(*),vmyx(*),vmyy(*),vmyz(*)
      dimension ioa(*),iob(*),inter(na,nact,2)
      dimension iwa(*),iwb(*),nwa(*),nwb(*)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
c
      iuaa = 0
      iubb = 0
      ia1  = 0
      call vclr(vmyx,1,nci1)
      call vclr(vmyy,1,nci1)
      call vclr(vmyz,1,nci1)
c
      do 220 isymb=1,nirr
      if (nstrb(isymb) .eq. 0) go to 220
      nblock = (nstrb(isymb)-1)/nbblk + 1
      isyma = mult(isymb,isym)
      if (nstra(isyma) .eq. 0) go to 220
      nalock = (nstra(isyma)-1)/nablk + 1
c
c   symmetry of excited beta strings
c
      isymbe = mult(isymb,isymr)
      nstrbe = nstrb(isymbe)
      mtaa = 0
      do 210 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)
c
c   generate alpha replacements <k|E_rs|n> with
c   symm(k_a)=isyma, symm(k_b)=isymb, symm(n)=isss
c
      mta = mtaa
      do 10 i=1,na
10    icga(i) = icgaa(i)
      iua = iuaa
      iofiwa = 1
      do 30 ia=1,naa
20    call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      if (jmta .ne. isyma) go to 20
      if (ia .eq. 1) ia0 = ioa(mta)
      call rplace (icga,isymr,iwa(iofiwa),nwa(ia),na,ioa,inter)
30    iofiwa = iofiwa + maxrp3
      ib1 = 0
      mtbb = 0
      do 190 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa * nbb
      naabbp = naabb * npair(isymr)
c
c   generate beta replacements <k|E_rs|n> with
c   symm(k_a)=isyma,symm(k_b)=isymb,symm(n)=isss
c
      mtb = mtbb
      iub = iubb
      do 40 i=1,nb
40    icgb(i) = icgbb(i)
      iofiwb = 1
      do 60 ib=1,nbb
50    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      if (jmtb .ne. isymb) go to 50
      if (ib .eq. 1) ib0 = iob(mtb) - 1
      call rplace (icgb,isymr,iwb(iofiwb),nwb(ib),nb,iob,inter(1,1,2))
60    iofiwb = iofiwb + maxrp3
c
c   construct d matrix
c
      call vclr(d,1,naabbp)
      iofiw1 = 0
      iofca = ib0
      iofd = - naabb
      do 100 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 90 iw=1,nwa(ia)
      ioffc = iofca + iwa(iofiw1+iw)
      ioffd = iofd + ipos(iwa(iofiw3+iw))*naabb
      if (iwa(iofiw2+iw) .gt. 0) then
         do 70 ib=1,nbb
70       d(ioffd+ib) = d(ioffd+ib) + c(ioffc+ib)
      else
         do 80 ib=1,nbb
80       d(ioffd+ib) = d(ioffd+ib) - c(ioffc+ib)
      end if
90    continue
      iofiw1 = iofiw1 + maxrp3
100   iofd = iofd + nbb
      iofd = 1 - nbb - naabb
      iofcb = ia0 - nstrbe
      iofiw1 = 0
      do 140 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 130 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + ipos(iwb(iofiw3+iw))*naabb
      if (iwb(iofiw2+iw) .gt. 0) then
         do 110 ia=1,naa
110      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c(ioffc+ia*nstrbe)
      else
         do 120 ia=1,naa
120      d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c(ioffc+ia*nstrbe)
      end if
130   continue
      iofiw1 = iofiw1 + maxrp3
140   iofd = iofd + 1
      if (odebug(68)) then
         print*,'d matrix from dipcon:'
         write (iwr,951) isyma,isymb,ialock,iblock
951      format ('d matrix for symmetries (',i1,',',i1,') -',
     ?           ' block (',i2,',',i2,'):')
         do 955 iv=1,nact
         do 955 ix=1,nact
         if (mult(itypea(iv),itypea(ix)) .eq. isymr) then
            write (iwr,952) iv,ix
952         format ('pair (',i2,',',i2,'):')
            k = (ipos(ilifa(iv)+ix)-1)*naabb
            do 953 ib=1,nbb
953         write (iwr,954) (d(k+(ia-1)*nbb+ib),ia=1,naa)
954         format (10f9.4)
         end if
955      continue
      end if
c
c   contributions to vmyx,vmyy,vmyz
c
      do 170 it=1,nact
      do 170 iu=1,it
      if (mult(itypea(it),itypea(iu)) .eq. isymr) then
         iofd1 = (ipos(ilifa(it)+iu)-1)*naabb
         iofd2 = (ipos(ilifa(iu)+it)-1)*naabb
         iofmy = ia1 + ib1
         k = (it+ncore+nfreez-1)*nbasis + iu+ncore+nfreez
         do 160 ia=1,naa
         do 150 ib=1,nbb
         vmyx(iofmy+ib) = vmyx(iofmy+ib) + 
     +                    dipx(k) * (d(iofd1+ib) + d(iofd2+ib))
         vmyy(iofmy+ib) = vmyy(iofmy+ib) + 
     +                    dipy(k) * (d(iofd1+ib) + d(iofd2+ib))
150      vmyz(iofmy+ib) = vmyz(iofmy+ib) + 
     +                    dipz(k) * (d(iofd1+ib) + d(iofd2+ib))
         iofd1 = iofd1 + nbb
         iofd2 = iofd2 + nbb
160      iofmy = iofmy + nstrb(isymb)
      end if
170   continue 
c
c   end of beta block
c
      do 180 i=1,nb
180   icgbb(i) = icgb(i)
      ib1 = ib1 + nbb
      iubb = iub
190   mtbb = mtb
c
c   end of alpha block
c
      do 200 i=1,na
200   icgaa(i) = icga(i)
      ia1 = ia1 + naa * nstrb(isymb)
      iuaa = iua
210   mtaa = mta
c
c   end of beta symmetry
c
220   continue
c
      return
      end
      subroutine diporb (dipx,dipy,dipz,vmyx,vmyy,vmyz,d)

c---------------------------------------------------------------------
c   Constructs the vector with entries  <0|[A,E_pq]|0>
c   for A = mu_x, mu_y, mu_z (components of the dipole moment operator).
c   This is the orbital part of the vector which enters in the
c   calculation of the oscillator strength.
c   (c) Carsten Fuchs 1991-1992
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension dipx(*),dipy(*),dipz(*),vmyx(*),vmyy(*),vmyz(*)
      dimension d(nact,nact)

      ind(i,j) = ((i-1)*(i-2))/2 + j

      do 210 ipi=1,ncore*(nact+nvirt)
      ip = (ipi-1)/ncore + 1
      ii = ipi - (ip-1)*ncore
      ip = ip + ncore
      if (mult(itype(ip),itype(ii)) .eq. isymr) then
         k = ipo(ind(ip,ii))
         l = (ip+nfreez-1)*nbasis + ii+nfreez
         vmyx(k) = 2.0d0 * dipx(l)
         vmyy(k) = 2.0d0 * dipy(l)
         vmyz(k) = 2.0d0 * dipz(l)
      endif
210   continue
      do 230 it=1,nact
      do 230 ii=1,ncore
      if (mult(itypea(it),itype(ii)) .eq. isymr) then
         k = ipo(ind(it+ncore,ii))
         do 220 iu=1,nact
         if (itypea(iu) .eq. itypea(it)) then
            l = (iu+ncore+nfreez-1)*nbasis + ii+nfreez
            vmyx(k) = vmyx(k) - dipx(l) * d(it,iu)
            vmyy(k) = vmyy(k) - dipy(l) * d(it,iu)
            vmyz(k) = vmyz(k) - dipz(l) * d(it,iu)
         endif
220      continue
      endif
230   continue
      k = norbti + norbai + 1
      call vclr(vmyx(k),1,norbat)
      call vclr(vmyy(k),1,norbat)
      call vclr(vmyz(k),1,norbat)
      do 250 ia=nprimp,nba
      do 250 it=1,nact
      if (mult(itype(ia),itypea(it)) .eq. isymr) then
         k = ipo(ind(ia,it+ncore))
         do 240 ix=1,nact
         if (itypea(ix) .eq. itypea(it)) then
            l = (ia+nfreez-1)*nbasis + ix+ncore+nfreez
            vmyx(k) = vmyx(k) + dipx(l) * d(it,ix)
            vmyy(k) = vmyy(k) + dipy(l) * d(it,ix)
            vmyz(k) = vmyz(k) + dipz(l) * d(it,ix)
         endif
240      continue
      endif
250   continue

      return
      end
      subroutine getdip (q)
c
c----------------------------------------------------------------------
c  Reads the dipole integrals in the ao basis and transforms 
c  them to the mcscf mo basis --- mcscf mos must first be 
c  transformed back to the non-symmetry-adapted basis
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer numdu, iblkdu, maxb, iblkla
      logical orevis
      real*8 apos
      common/sector/numdu,iblkdu,orevis(2),
     + apos(508),maxb,iblkla
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical odebug
      common /debug/ odebug(100)
c
      character*8 obnam
      common /craypl/ obnam(29)
c
      integer mxorb1
      integer iky, ilifq, mapie, ilifm, i4096
      parameter (mxorb1=maxorb+1)
      common /blocko/ iky(mxorb1),ilifq(mxorb1),mapie(mxorb1),
     +               ilifm(mxorb1),i4096(maxorb)
c
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      dimension o1e(6)
      dimension q(*)
      data o1e/.false.,.false.,.true.,.true.,.true.,.true./

c   read in information about symmetry adaption

      call secget (isecao,3,iblk)
      call rdchr (obnam,29,iblk,ned3)
      call reads (evalue,mach(8),ned3)
      nav = lenwrd()
      call readis (ilifc,mach(9)*nav,ned3)
      nb1 = nbasis * nbasis
      nb2 = (nbasis*(nbasis+1))/2
      i1 = 1
      i2 = i1 + nb1
      ih = i2 + nb1
      ix = ih + nb2
      iy = ix + nb2
      iz = iy + nb2
      i3 = iz + nb2
      i4 = i3 + nb1
      do 10 i=1,nbasis+1
10    ilifq(i) = (i-1)*nbasis
      if (odebug(84)) then
         call outive (ilifc,nbasis,'ilifc')
         call outive (ntran,nbasis,'ntran')
         call outive (itran,nbasis,'itran')
         call outvec (ctran,nbasis,'ctran')
      end if

c   read in mcscf mos

      call secget (isecb,3,iblko)
      iblkq = iblko + lensec(mach(8)) +lensec(mach(9)) + 1
      call rdedx (q(i1),nb1,iblkq,ned3)
      num = nbasis
      otran = .false.

c   transform mcscf mos to non-symmetry-adapted ao basis

      call tdown (q(i2),ilifq,q(i1),ilifq,nbasis)
      if (odebug(84)) 
     ?   call outsqr (q(i2),nbasis,nbasis,nbasis,'scattered mcscf mos')

c   read in t+v plus dipole integrals

      call getmat(q(ih),q(ih),q(ih),q(ix),q(iy),q(iz),
     +            potnuc,nbasis,o1e,isect(492))
c
      core = potnuc
c
      isec134 = isect(134)
      m80 = 80
      call secput (isec134,m80,4*lensec(nb1),idiblo)

c   transform to mo basis

      iq = ih
      do 2500 icoord=1,4
      call square (q(i3),q(iq),nbasis,nbasis)
      call mxmaa (q(i3),1,nbasis, q(i2),1,nbasis, q(i4),1,nbasis,
     ?            nbasis,nbasis,nbasis)
      call mxmaa (q(i2),nbasis,1, q(i4),1,nbasis, q(i3),1,nbasis,
     ?            nbasis,nbasis,nbasis)
      if (odebug(84)) then
         print*,'dipole integrals, coordinate ',icoord
         call outsqr (q(i3),nbasis,nbasis,nbasis,' ')
      end if
      call wrt3 (q(i3),nb1,idiblo,ned3)
      iq = iq + nb2
      idiblo = idiblo + lensec(nb1)
2500  continue

      return
      end
      subroutine intord (q,d)
c
c---------------------------------------------------------------------
c   Reorders the two-electron integrals
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
      integer seclen
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      dimension q(*),d(nact,nact),fac(2)
      dimension intwm(8),jwm(8),kwm(8),nwm(8),jblkwm(8)
      dimension mk1(8),mk2(8),mk3(8),mk4(8),mk5(8)
      dimension ibx1of(8,8),ibx2of(8,8),iof(8)
      dimension iad(980),ifill(980)
      dimension ihlp(3,8,50)
      data fac/2.0d0,-1.0d0/
c
      if(odebug(50).or.odebug(78)) then
         do i=1,ncore
         ihlp(1,itype(i),npoint(i))=i
         end do
         do i=nst,nprim
         ihlp(2,itype(i),npoint(i))=i
         end do
         do i=nprimp,nba
         ihlp(3,itype(i),npoint(i))=i
         end do
      end if
      ibase = icorr(0)
      time = cpulft(1)
      iblo4 = 1
      k = 0
      do 2050 isyma=1,nirr
         iof(isyma) = k
         k = k + nsymm(isyma)**2
2050  continue
c
c   space
c
      i1 = icorr(nbasq)
      i2 = icorr(nbasq)
      ifi = icorr(nfock)
      ift = icorr(nfock)
      nbasr = nbasq + 511
      joo = icorr(nbasr)
      jgv = icorr(nbasr)
      jaa = icorr(nbasr)
      do 90 i=1,nirr
         jwm(i) = icorr(nbasr)
90    continue
      mmax = icorrm()-1000
      maxhlf = mmax/2
      ir = icorr(maxhlf)
      if (odebug(78)) print*,'ir = ',ir
      iw = icorr(maxhlf)
      ir1 = ir-1
c     iw1 = iw-1
      n = 0
      do 5 i=1,nirr
         k = 0
         do 4 isyma=1,nirr
            isymb = mult(i,isyma)
            k = k + 2*mact(isyma)*mact(isymb) 
     ?            + 3*mact(isyma)*mcore(isymb)
4        continue
         n = max(n,k)
5     continue
      nprbox = maxhlf/n
c...define box offsets
      do 42 isym = 1 , nirr
         do 41 isyma = 1 , nirr
            ibx1of(isym,isyma) = 0
            ibx2of(isym,isyma) = 0
41       continue
42    continue
      iboff = 0
      do 22 isym = 1 , nirr
         do 2 isyma = 1 , nirr
            isymb = mult(isyma,isym)
            ibx1of(isym,isyma) = iboff
            m = mvirt(isyma) * mcore(isymb)
            if (m .gt. 0) iboff = iboff + (m-1)/nprbox + 1
2        continue
22    continue
      do 33 isym = 1 , nirr
         do 3 isyma = 1 , nirr
            isymb = mult(isyma,isym)
            if (isymb .gt. isyma) go to 3
            ibx2of(isym,isyma) = iboff
            m = mvirt(isyma) * mvirt(isymb)
            if (isyma .eq. isymb) m = (mvirt(isyma)*(mvirt(isyma)+1))/2
            if (m .gt. 0) iboff = iboff + (m-1)/nprbox + 1
3        continue
33    continue
      ntot = iboff
      write (iwr,1) mmax,n,nprbox,ntot
1     format (/'INTINT reordering parameters:'//
     ?         'size of scratch space = ',i8/
     ?         'max no of (tu) pairs  = ',i8/
     ?         '# pairs per box       = ',i8/
     ?         '# boxes needed        = ',i8/)
      print*,'offsets of boxes (ai):'
      do 11 i=1,nirr
11    write(iwr,12)(ibx1of(i,isyma),isyma=1,nirr)
12    format(8i7)
      print*,'offsets of boxes (ab):'
      do 13 i=1,nirr
13    write(iwr,12)(ibx2of(i,isyma),isyma=1,nirr)
      if (ntot*511 .gt. mmax)
     ?   call caserr ('insufficient space for integral reordering')
      if (ntot .gt. 980)
     ?   call caserr ('more than 980 boxes in INTINT !')
      call setsto(ntot,0,iad)
      call setsto(ntot,0,ifill)
c...count integrals
c
c...integrals for double commutators
c...blocks 1,2,4
      int124 = 0
      do 10 i=1,ncore
      do 10 j=1,i
      isymij = mult(itype(i),itype(j))
      do 10 isyma=1,nirr
      isymb = mult(isyma,isymij)
      m = (mact(isyma)+mvirt(isyma))*(mact(isymb)+mvirt(isymb))
      if (isymb .le. isyma) then
         int124 = int124 + 2*m
      else
         int124 = int124 + m
      end if
10    continue
c...blocks 3,5,6
      int356 = 0
      do 20 i=nst,nprim
      do 20 j=1,ncore
      isymij = mult(itype(i),itype(j))
      do 20 isyma=1,nirr
      isymb = mult(isyma,isymij)
      int356 = int356 + 2*(mvirt(isyma)*mact(isymb)
     ?                + (mact(isyma)+mvirt(isyma))*mvirt(isymb))
20    continue
      do 30 i=nst,nprim
      do 30 j=nst,i
      isymij = mult(itype(i),itype(j))
      do 30 isyma=1,nirr
      isymb = mult(isyma,isymij)
      int356 = int356 + 2*mvirt(isyma)*mvirt(isymb)
30    continue
c...integrals for gv,gw
      intgv = 0
      do 40 i=nst,nprim
      do 40 j=nst,i
      isymij = mult(itype(i),itype(j))
      do 40 isyma=1,nirr
      isymb = mult(isyma,isymij)
      if (isymb .gt. isyma) go to 40
      intgv = intgv + mact(isyma)*(mcore(isymb)+mvirt(isymb))
     ?              + mact(isymb)*(mcore(isyma)+mvirt(isyma))
40    continue
c...all-active integrals
      intaa = 0
      do 55 i=nst,nprim
      do 55 j=nst,i
      isymij = mult(itype(i),itype(j))
      do 50 isyma=1,nirr
      isymb = mult(isyma,isymij)
      if (isymb .gt. isyma) go to 50
      intaa = intaa + mact(isyma)*mact(isymb)
50    continue
      if (isymij .ne. 1) go to 55
      do 52 isyma=1,nirr
52    intaa = intaa + mact(isyma)**2
55    continue
c...integrals for matrix W
c      do 70 i=1,nirr
c      irow=0
c      icol=0
c      do 60 isyma=1,nirr
c      isymb=mult(isyma,i)
c      irow=irow+(mact(isyma)+mvirt(isyma))*mcore(isymb)
c      icol=icol+mact(isyma)*mact(isymb)
c60    continue
c70    intwm(i)=irow*icol
      call setsto(nirr,0,intwm)
      do 70 i=nst,nprim
      do 70 j=nst,i
      isymij = mult(itype(i),itype(j))
      do 60 isyma=1,nirr
      isymb = mult(isyma,isymij)
      m = (mact(isyma)+mvirt(isyma))*mcore(isymb)+
     ?    (mact(isymb)+mvirt(isymb))*mcore(isyma)
      if (isymb .le. isyma) then
         intwm(isymij) = intwm(isymij) + 2*m
      else
         intwm(isymij) = intwm(isymij) + m
      end if
60    continue
70    continue
c...define block addresses on ed0
      ibl124 = 1
      ibl356 = ibl124 + seclen(int124)
      iblkgv = ibl356 + seclen(int356)
      iblkaa = iblkgv + seclen(intgv)
      iblk   = iblkaa + seclen(intaa)
      do 80 i=1,nirr
      iblkwm(i) = iblk
      iblk = iblk + seclen(intwm(i))
80    continue
c...dummy block after last W matrix so that INTIN causes no trouble
      zero = 0.0d0
      call wrt3 (zero,1,iblk,ned0)
      if (odebug(50) .or. odebug(78) .or. odebug(83)) then
         print*,'number of integrals:'
         print*,'int124 = ',int124
         print*,'int356 = ',int356
         print*,'intgv = ',intgv
         print*,'intaa = ',intaa
         do i=1,nirr
         print*,'intwm(',i,') = ',intwm(i)
         end do
         print*,'needing ',iblk-1,' blocks on ed0'
         print*,'starting blocks:'
         print*,'ibl124 = ',ibl124
         print*,'ibl356 = ',ibl356
         print*,'iblkgv = ',iblkgv
         print*,'iblkaa = ',iblkaa
         do i=1,nirr
         print*,'iblkwm(',i,') = ',iblkwm(i)
         end do
      end if
      if (iblk .gt. 100000) 
     ?   call caserr ('insufficient space for integrals')
c
      call vclr(q(ift),1,nfock)
c...read in integrals from ed6
      intpos = 0
      intfil = ned6
      iblf = 1
      intmod = 0
      koo = joo
      noo = 0
      jblkoo = ibl124
c...integrals (ij|..),(i.|.j)
*****************************
      do 1000 i=1,ncore
      do 1000 j=1,i
      isymij = mult(itype(i),itype(j))
      iexch = 0
      isyma = 0
110   isyma = isyma+1
      if (isyma .gt. nirr) go to 190
      isymb = mult(isyma,isymij)
      n = nsymm(isyma)*nsymm(isymb)
      if (n .eq. 0) go to 110
      if (iexch .eq. 1) go to 120
      if (isymb-isyma) 120,130,110
120   call intin (q(i1),n)
      go to 140
130   call intin (q(i2),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call square (q(i1),q(i2),nsymm(isyma),nsymm(isyma))
c...integrals (ij|tu),(iu|tj)
140   ip = i1+mcore(isymb)*nsymm(isyma)+mcore(isyma)
      do 150 k=1,mact(isymb)
      call fmove (q(ip),q(koo),mact(isyma))
      if (odebug(82)) then
         print*,'moved (ii|tt) integrals:'
         call outvec (q(koo),mact(isyma),' ')
      end if
      ip = ip+nsymm(isyma)
      koo = koo+mact(isyma)
150   noo = noo+mact(isyma)
c...integrals (ij|at),(it|aj)
      ip = i1+mcore(isymb)*nsymm(isyma)+nprm(isyma)
      do 160 k=1,mact(isymb)
      call fmove (q(ip),q(koo),mvirt(isyma))
      if(odebug(82))then
         print*,'moved (ii|at) integrals:'
         call outvec (q(koo),mvirt(isyma),' ')
      end if
      ip = ip+nsymm(isyma)
      koo = koo+mvirt(isyma)
160   noo = noo+mvirt(isyma)
c...integrals (ij|ta),(ia|tj)
      ip = i1+nprm(isymb)*nsymm(isyma)+mcore(isyma)
      do 170 k=1,mvirt(isymb)
      call fmove (q(ip),q(koo),mact(isyma))
      if(odebug(82))then
         print*,'moved (ii|ta) integrals:'
         call outvec (q(koo),mact(isyma),' ')
      end if
      ip = ip+nsymm(isyma)
      koo = koo+mact(isyma)
170   noo = noo+mact(isyma)
c...integrals (ij|ab),(ib|aj)
      ip = i1+nprm(isymb)*nsymm(isyma)+nprm(isyma)
      do 180 k=1,mvirt(isymb)
      call fmove (q(ip),q(koo),mvirt(isyma))
      if(odebug(82))then
         print*,'moved (ii|aa) integrals:'
         call outvec (q(koo),mvirt(isyma),' ')
      end if
      ip = ip+nsymm(isyma)
      koo = koo+mvirt(isyma)
180   noo = noo+mvirt(isyma)
      if(odebug(11)) print*,'calling outint with oo integrals'
      call outint (q(joo),noo,jblkoo,ned0)
      koo = joo+noo
c...contribution to core
      if (i .eq. j) 
     ?   core = core - 
     ?   fac(iexch+1)*dsum(mcore(isyma),q(i1),nsymm(isyma)+1)
      go to 110
190   continue
      if (iexch .eq. 1) go to 1000
      iexch = 1
      isyma = 0
      go to 110
1000  continue
      if (noo .gt. 0) call wrt3 (q(joo),noo,jblkoo,ned0)
      if (odebug(83)) print*,'writing ',noo,
     ?                       ' 124-integrals on block ',jblkoo
c...integrals (ti|..),(t.|.i),(tu|..),(t.|.u)
*********************************************
      kgv = jgv
      kaa = jaa
      do 92 i=1,nirr
92    kwm(i) = jwm(i)
      ngv = 0
      naa = 0
      do 94 i=1,nirr
94    nwm(i) = 0
      jblkgv = iblkgv
      jblkaa = iblkaa
      do 96 i=1,nirr
96    jblkwm(i) = iblkwm(i)
      do 2000 i=nst,nprim
      do 2000 j=1,i
      isymij = mult(itype(i),itype(j))
      if(odebug(50))print*,'entering pair ',i,j,' --- symmetry ',isymij
      iexch = 0
200   j1 = i1
c...read in integrals (all symmetries)
      do 230 isyma=1,nirr
      isymb = mult(isyma,isymij)
      n = nsymm(isyma)*nsymm(isymb)
      if (n .eq. 0) go to 230
      if (iexch .eq. 1) go to 220
      if (isymb-isyma) 220,210,230
210   call intin (q(i2),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call square (q(j1),q(i2),nsymm(isyma),nsymm(isyma))
      go to 225
220   call intin (q(j1),n)
225   mk1(isyma) = j1
      j1 = j1+n
230   continue
c...contributions to active fock matrix
      if (j.le.ncore .or. isymij.ne.1) go to 235
      factor = 2.0d0
      if (i .eq. j) factor = 1.0d0
      if (iexch .eq. 1) factor = -0.5d0
      factor = factor*d(i-ncore,j-ncore)
      do 232 isyma=1,nirr
      call daxpy(nsymm(isyma)**2
     +          ,factor,q(mk1(isyma)),1
     +          ,q(ift+iof(isyma)),1)
      if (iexch.eq.0 .or. i.eq.j) go to 232
      call trnsps (q(mk1(isyma)),q(i2),nsymm(isyma),nsymm(isyma))
      call daxpy(nsymm(isyma)**2
     +           ,factor,q(i2),1
     +           ,q(ift+iof(isyma)),1)
232   continue
235   if (iexch .eq. 1) go to 360
c========================================================
c   distribute coulomb integrals (ti|ab) resp. (tu|ab)
c========================================================
      do 260 isyma=1,nirr
      isymb = mult(isyma,isymij)
      if (isymb .gt. isyma) go to 260
      ip = mk1(isyma)+nprm(isymb)*nsymm(isyma)+nprm(isyma)-1
      ib1 = mvirt(isymb)
      do 250 ia=1,mvirt(isyma)
      if (isymb .eq. isyma) ib1=ia
      do 240 ib=1,ib1
      ibox = ibx2of(isymij,isyma) + 
     ?       ((ib-1)*mvirt(isyma)+ia-1)/nprbox + 1
      if (isymb .eq. isyma) 
     ?   ibox = ibx2of(isymij,isyma) + ((ia*(ia-1))/2+ib-1)/nprbox + 1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+(ib-1)*nsymm(isyma)+ia)
      if (odebug(50)) write(iwr,237) i,j,ihlp(3,isyma,ia),
     ?                               ihlp(3,isymb,ib),q(mlower+kk),ibox
237   format('putting (',2i3,'|',2i3,') = ',f14.7,' in box ',i3)
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
240   continue
250   continue
260   continue
      if(j.le.ncore) go to 1999
c==========================================
c   distribute coulomb integrals (tu|ai)
c==========================================
      do 350 isyma=1,nirr
      isymb=mult(isyma,isymij)
      if(isymb.gt.isyma) go to 350
      ip=mk1(isyma)-1+nprm(isyma)
c...distribute integrals (tu|ai) with symm(a)=isyma, symm(i)=isymb
      do 320 ii=1,mcore(isymb)
      do 310 ia=1,mvirt(isyma)
      ibox=ibx1of(isymij,isyma)+((ii-1)*mvirt(isyma)+ia-1)/nprbox+1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+ia)
      if (odebug(50)) write (iwr,237) i,j,ihlp(3,isyma,ia),
     ?                                ihlp(1,isymb,ii),q(ip+ia),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
310   continue
320   ip=ip+nsymm(isyma)
      if(isyma.eq.isymb) go to 350
c...distribute integrals (tu|ia) with symm(i)=isyma, symm(a)=isymb
      ip=mk1(isyma)-1+nprm(isymb)*nsymm(isyma)
      do 340 ii=1,mcore(isyma)
      do 330 ia=1,mvirt(isymb)
      ibox=ibx1of(isymij,isymb)+((ii-1)*mvirt(isymb)+ia-1)/nprbox+1
      kk=ifill(ibox)+1
      mlower=(ibox-1)*511
      q(mlower+kk)=q(ip+(ia-1)*nsymm(isyma)+ii)
      if (odebug(50)) write (iwr,237) i,j,ihlp(1,isyma,ii),
     ?                                ihlp(3,isymb,ia),q(mlower+kk),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
330   continue
340   continue
350   continue
      go to 390
360   continue
      if(j.le.ncore) go to 369
c=====================================================
c   distribute exchange integrals (ta|iu) = (at|iu)
c                                  -   -     -  -
c=====================================================
      do 365 isyma=1,nirr
      isymb=mult(isyma,isymij)
      ip=mk1(isymb)-1+nprm(isyma)*nsymm(isymb)
      do 364 ii=1,mcore(isymb)
      do 363 ia=1,mvirt(isyma)
      ibox=ibx1of(isymij,isyma)+((ii-1)*mvirt(isyma)+ia-1)/nprbox+1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+(ia-1)*nsymm(isymb)+ii)
      if (odebug(50)) write (iwr,237) ihlp(3,isyma,ia),i,
     ?                                ihlp(1,isymb,ii),j,
     ?                                q(mlower+kk),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
363   continue
364   continue
365   continue
c=====================================================
c   distribute exchange integrals (ti|au) = (au|it)
c                                  -   -     -  -
c=====================================================
      do 368 isyma=1,nirr
      isymb=mult(isyma,isymij)
      ip=mk1(isyma)-1+nprm(isyma)
      do 367 ii=1,mcore(isymb)
      do 366 ia=1,mvirt(isyma)
      ibox=ibx1of(isymij,isyma)+((ii-1)*mvirt(isyma)+ia-1)/nprbox+1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+ia)
      if (odebug(50)) write (iwr,237) ihlp(3,isyma,ia),j,
     ?                                ihlp(1,isymb,ii),i,q(ip+ia),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
366   continue
367   ip=ip+nsymm(isyma)
368   continue
369   continue
c=======================================================================
c   distribute exchange integrals (ta|bi) = (at|bi) , (ta|bu) = (at|bu)
c                                  -   -     -  -      -   -     -  -
c=======================================================================
      do 380 isyma=1,nirr
      isymb=mult(isyma,isymij)
      if(isymb.gt.isyma) go to 380
      ip=mk1(isymb)-1+nprm(isyma)*nsymm(isymb)+nprm(isymb)
      ib1=mvirt(isymb)
      do 375 ia=1,mvirt(isyma)
      if(isymb.eq.isyma) ib1=ia
      do 370 ib=1,ib1
      ibox=ibx2of(isymij,isyma)+((ib-1)*mvirt(isyma)+ia-1)/nprbox+1
      if (isymb .eq. isyma) 
     ?   ibox = ibx2of(isymij,isyma) + ((ia*(ia-1))/2+ib-1)/nprbox + 1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+ib)
      if (odebug(50)) write (iwr,237) ihlp(3,isyma,ia),i,
     ?                                ihlp(3,isymb,ib),j,q(ip+ib),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
370   continue
375   ip=ip+nsymm(isymb)
380   continue
c=======================================================================
c   distribute exchange integrals (tb|ai) = (ai|bt) , (tb|au) = (au|bt)
c                                  -   -     -  -      -   -     -  -
c=======================================================================
      do 388 isyma=1,nirr
      isymb=mult(isyma,isymij)
      if(isymb.gt.isyma) go to 388
      ip=mk1(isyma)-1+nprm(isymb)*nsymm(isyma)+nprm(isyma)
      ib1=mvirt(isymb)
      do 385 ia=1,mvirt(isyma)
      if(isymb.eq.isyma) ib1=ia
      do 383 ib=1,ib1
      ibox=ibx2of(isymij,isyma)+((ib-1)*mvirt(isyma)+ia-1)/nprbox+1
      if (isymb .eq. isyma) 
     ?   ibox = ibx2of(isymij,isyma) + ((ia*(ia-1))/2+ib-1)/nprbox + 1
      kk=ifill(ibox)+1
      mlower=ir1+(ibox-1)*511
      q(mlower+kk)=q(ip+(ib-1)*nsymm(isyma)+ia)
      if (odebug(50)) write (iwr,237) ihlp(3,isyma,ia),j,
     ?                                ihlp(3,isymb,ib),i,
     ?                                q(mlower+kk),ibox
      if(kk.ge.510) then
         q(mlower+511)=dble(iad(ibox))
         call wrt3 (q(mlower+1),511,iblo4,ned4)
         iad(ibox)=iblo4*1000+510
         iblo4=iblo4+1
         call vclr(q(mlower+1),1,511)
         ifill(ibox)=0
      else
         ifill(ibox)=kk
      end if
383   continue
385   continue
388   continue
390   continue
      if(j.le.ncore) go to 1999
c...integrals (tt|..) for remaining tasks
      do 580 isyma=1,nirr
      isymb=mult(isyma,isymij)
      if(iexch.eq.0.and.isymb.gt.isyma) go to 580
      ma=nsymm(isyma)
      mca=mcore(isyma)
      mcb=mcore(isymb)
      maa=mact(isyma)
      mab=mact(isymb)
      mva=mvirt(isyma)
      mvb=mvirt(isymb)
c     mava=maa+mva
      if(iexch.eq.1) go to 475
c...integrals (tt|ti) for gv
      ip=mk1(isyma)+mca
      do 400 k=1,mcb
      call fmove (q(ip),q(kgv),maa)
      ip=ip+ma
      kgv=kgv+maa
400   ngv=ngv+maa
c...integrals (tt|it) for gv
      ip=mk1(isyma)+mcb*ma
      do 430 k=1,mab
      call fmove (q(ip),q(kgv),mca)
      ip=ip+ma
      kgv=kgv+mca
430   ngv=ngv+mca
c...integrals (tt|at) for gv
      ip=mk1(isyma)+mcb*ma+mca+maa
      do 450 k=1,mab
      call fmove (q(ip),q(kgv),mva)
      ip=ip+ma
      kgv=kgv+mva
450   ngv=ngv+mva
c...integrals (tt|ta) for gv
      ip=mk1(isyma)+(mcb+mab)*ma+mca
      do 470 k=1,mvb
      call fmove (q(ip),q(kgv),maa)
      ip=ip+ma
      kgv=kgv+maa
470   ngv=ngv+maa
c...all-active integrals
475   if(iexch.eq.1.and.isymij.ne.1) go to 445
      ip=mk1(isyma)+mcb*ma+mca
      do 440 k=1,mab
      call fmove (q(ip),q(kaa),maa)
      ip=ip+ma
      kaa=kaa+maa
440   naa=naa+maa
c...integrals (tt|ti) for W matrix
445   ip=mk1(isyma)+mca
      do 410 k=1,maa
      call dcopy(mcb,q(ip),ma,q(kwm(isymij)),1)
      ip=ip+1
      kwm(isymij)=kwm(isymij)+mcb
410   nwm(isymij)=nwm(isymij)+mcb
c...integrals (tt|it) for W matrix
      ip=mk1(isyma)+mcb*ma
      do 437 k=1,mab
      call fmove (q(ip),q(kwm(isymij)),mca)
      ip=ip+ma
      kwm(isymij)=kwm(isymij)+mca
437   nwm(isymij)=nwm(isymij)+mca
c...integrals (tt|ai) for W matrix
      ip=mk1(isyma)+mca+maa
      do 420 k=1,mva
      call dcopy(mcb,q(ip),ma,q(kwm(isymij)),1)
      ip=ip+1
      kwm(isymij)=kwm(isymij)+mcb
420   nwm(isymij)=nwm(isymij)+mcb
c...integrals (tt|ia) for W matrix
      ip=mk1(isyma)+(mcb+mab)*ma
      do 460 k=1,mvb
      call fmove (q(ip),q(kwm(isymij)),mca)
      ip=ip+ma
      kwm(isymij)=kwm(isymij)+mca
460   nwm(isymij)=nwm(isymij)+mca
580   continue
      if(odebug(11)) print*,'calling outint with gv integrals'
      call outint (q(jgv),ngv,jblkgv,ned0)
      kgv=jgv+ngv
      if(odebug(11)) print*,'calling outint with aa integrals'
      call outint (q(jaa),naa,jblkaa,ned0)
      kaa=jaa+naa
      if (odebug(11)) 
     ?   print*,'calling outint with wm integrals, symmetry',isymij
      call outint (q(jwm(isymij)),nwm(isymij),jblkwm(isymij),ned0)
      kwm(isymij)=jwm(isymij)+nwm(isymij)
1999  continue
      if(iexch.eq.1) go to 2000
      iexch=1
      go to 200
2000  continue
      if (ngv .gt. 0) call wrt3 (q(jgv),ngv,jblkgv,ned0)
      if (odebug(83)) print*,'writing ',ngv,
     ?                       '  gv-integrals on block ',jblkgv
      if (naa .gt. 0) call wrt3 (q(jaa),naa,jblkaa,ned0)
      if (odebug(83)) print*,'writing ',naa,
     ?                       '  aa-integrals on block ',jblkaa
      do 510 i=1,nirr
      if (odebug(83)) print*,'writing ',nwm(i),
     ?                       '  wm-integrals of symmetry ',i,
     ?                       ' on block ',jblkwm(i)
510   if (nwm(i) .gt. 0) call wrt3 (q(jwm(i)),nwm(i),jblkwm(i),ned0)
c...finally write partially filled blocks onto external file
      do 520 ibox=1,ntot
      kk = ifill(ibox)
      if (kk .eq. 0) go to 520
      mlower = ir1+(ibox-1)*511
      q(mlower+511) = dble(iad(ibox))
      call wrt3 (q(mlower+1),511,iblo4,ned4)
      iad(ibox) = iblo4*1000+kk
      iblo4 = iblo4+1
520   continue
c...read in inactive fock matrix
      ifock = ifi
      do 530 isyma=1,nirr
      if (nsymm(isyma) .eq. 0) go to 530
      call intin (q(i1),(nsymm(isyma)*(nsymm(isyma)+1))/2)
      call square (q(ifock),q(i1),nsymm(isyma),nsymm(isyma))
      core = core + 2.0d0 * dsum(mcore(isyma),q(ifock),nsymm(isyma)+1)
      ifock = ifock+nsymm(isyma)**2
530   continue
      if (odebug(12)) print*,'core = ',core
c-----------------------
c   reordering step 2
c-----------------------
      iblk=ibl356
      iws=iw
      nw=0
c...first (ai) boxes
      do 730 isym=1,nirr
c...determine how many integrals (ai|tu),(at|iu) 
c...there will be for a pair (ai) of symmetry ISYM, and set marks
      listai=0
      do 610 isyma=1,nirr
      isymb=mult(isyma,isym)
      if(isymb.gt.isyma) go to 610
      mk1(isyma)=listai
      if(isymb.eq.isyma) listai=listai+(mact(isyma)*(mact(isyma)+1))/2
      if(isymb.lt.isyma) listai=listai+mact(isyma)*mact(isymb)
610   continue
      do 620 isyma=1,nirr
      isymb=mult(isyma,isym)
      mk2(isyma)=listai
620   listai=listai+mact(isyma)*mact(isymb)
      inlist=listai
      if(isym.eq.1) inlist=inlist+nact
      if(odebug(78))print*,'inlist = ',inlist
      if (inlist .eq. 0) go to 730
c...read in partially ordered integrals
      do 720 isyma=1,nirr
      isymb = mult(isyma,isym)
      ipair = 0
      np = nprbox
      npr = mvirt(isyma)*mcore(isymb)
      if (npr .eq. 0) go to 720
      mbox = (npr-1)/nprbox+1
      do 710 ibox=ibx1of(isym,isyma)+1,ibx1of(isym,isyma)+mbox
      if(odebug(78))print*,'treating box ',ibox
      if (ibox .eq. ibx1of(isym,isyma)+mbox) np = npr-ipair
      nints = np*inlist
      if(odebug(78))print*,'np = ',np,'   nints = ',nints
      temp = 0.0d0
      j1 = ir+510*((nints-1)/510)
      kk = iad(ibox)
630   iblo = kk/1000
      if(odebug(78))print*,'iblo = ',iblo
      if(odebug(78))print*,'j1 = ',j1
      if (iblo .eq. 0) go to 640
      m = ii-iblo*1000
      call rdedx (q(j1),511,iblo,ned4)
      kk=idint(q(j1+510))
      q(j1+510)=temp
      temp=q(j1)
      j1=j1-510
      go to 630
c...distribute integrals
640   if(j1+510.ne.ir) call caserr ('something wrong in intint ...')
      kpos1=ir1
      do 700 it=nst,nprim
      do 700 iu=nst,it
      if(mult(itype(it),itype(iu)).ne.isym) go to 700
      kpos2=kpos1+np
      kpos3=kpos2+np
      isymt=itype(it)
      isymu=itype(iu)
      if(isymt-isymu) 670,660,650
650   l1=mk1(isymt)+(npoint(iu)-1)*mact(isymt)+npoint(it)
      go to 680
660   l1=mk1(isymt)+(npoint(it)*(npoint(it)-1))/2+npoint(iu)
      go to 680
670   l1 = mk1(isymu) + (npoint(it)-1) * mact(isymu) + npoint(iu)
680   l2 = mk2(isymt) + (npoint(iu)-1) * mact(isymt) + npoint(it)
      l3 = mk2(isymu) + (npoint(it)-1) * mact(isymu) + npoint(iu)
c...for vector computers: no recurrence in this loop
      do 690 n = 1 , np
         q(iws-1+(n-1)*listai+l1) = q(kpos1+n)
         q(iws-1+(n-1)*listai+l2) = q(kpos2+n)
         q(iws-1+(n-1)*listai+l3) = q(kpos3+n)
690   continue
      kpos1 = kpos1 + 3 * np
700   continue
c...output reordered integrals
      nw=nw+np*listai
      call outint (q(iw),nw,iblk,ned0)
      iws=iw+nw
710   ipair=ipair+np
720   continue
730   continue
c...now (ab) boxes
      do 930 isym=1,nirr
c...determine how many integrals (ab|ti),(at|bi),(ai|bt),(ab|tu),(at|bu)
c...there will be for a pair (ab) of symmetry ISYM, and set marks
      listab=0
      do 750 isyma=1,nirr
      isymb=mult(isyma,isym)
      mk1(isyma)=listab
750   listab=listab+mact(isyma)*mcore(isymb)
      do 760 isyma=1,nirr
      isymb=mult(isyma,isym)
      mk2(isyma)=listab
760   listab=listab+mact(isyma)*mcore(isymb)
      do 770 isyma=1,nirr
      isymb=mult(isyma,isym)
      mk3(isyma)=listab
770   listab=listab+mact(isyma)*mcore(isymb)
      do 780 isyma=1,nirr
      isymb=mult(isyma,isym)
      if(isymb.gt.isyma) go to 780
      mk4(isyma)=listab
      if(isymb.eq.isyma) listab=listab+(mact(isyma)*(mact(isyma)+1))/2
      if(isymb.lt.isyma) listab=listab+mact(isyma)*mact(isymb)
780   continue
      do 790 isyma=1,nirr
      isymb=mult(isyma,isym)
      mk5(isyma)=listab
790   listab=listab+mact(isyma)*mact(isymb)
      inlist=listab
      if(isym.eq.1) inlist=inlist+nact
      if(odebug(78))print*,'inlist = ',inlist
      if (inlist .eq. 0) go to 930
c...read in partially ordered integrals
      do 920 isyma=1,nirr
      isymb = mult(isyma,isym)
      if (isymb .gt. isyma) go to 920
      ipair = 0
      np = nprbox
      npr = mvirt(isyma)*mvirt(isymb)
      if (isyma .eq. isymb) npr = (mvirt(isyma)*(mvirt(isyma)+1))/2
      if (npr .eq. 0) go to 920
      mbox = (npr-1)/nprbox+1
      do 910 ibox=ibx2of(isym,isyma)+1,ibx2of(isym,isyma)+mbox
      if (odebug(78)) print*,'treating box ',ibox
      if (ibox .eq. ibx2of(isym,isyma)+mbox) np = npr-ipair
      nints = np*inlist
      if (odebug(78)) print*,'np = ',np,'   nints = ',nints
      temp = 0.0d0
      j1 = ir + 510*((nints-1)/510)
      kk = iad(ibox)
810   iblo = kk/1000
      if (odebug(78)) print*,'iblo = ',iblo
      if (odebug(78)) print*,'j1 = ',j1
      if (iblo .eq. 0) go to 820
      m = ii-iblo*1000
      call rdedx (q(j1),511,iblo,ned4)
      kk = idint(q(j1+510))
      q(j1+510) = temp
      temp = q(j1)
      j1 = j1 - 510
      go to 810
c...distribute integrals
820   if(j1+510.ne.ir) call caserr ('something wrong in intint !')
      kpos1=ir1
c      if(odebug(78))then
c  print*,'semiordered integrals:'
c         do l=1,inlist
c  call outvec (q(ir+(l-1)*3*np),3*np,' ')
c  end do
c      end if
      do 900 it=nst,nprim
      do 900 iu=1,it
      if (mult(itype(it),itype(iu)) .ne. isym) go to 900
      kpos2 = kpos1 + np
      kpos3 = kpos2 + np
      isymt = itype(it)
      isymu = itype(iu)
      if (iu .gt. ncore) go to 830
      m = (npoint(iu)-1)*mact(isymt)+npoint(it)
      l1 = mk1(isymt)+m
      l2 = mk2(isymt)+m
      l3 = mk3(isymt)+m
      go to 880
830   if (isymt-isymu) 860,850,840
840   l1 = mk4(isymt) + (npoint(iu)-1) * mact(isymt) + npoint(it)
      go to 870
850   l1 = mk4(isymt) + (npoint(it)*(npoint(it)-1))/2 + npoint(iu)
      go to 870
860   l1 = mk4(isymu) + (npoint(it)-1) * mact(isymu) + npoint(iu)
870   l2 = mk5(isymt) + (npoint(iu)-1) * mact(isymt) + npoint(it)
      l3 = mk5(isymu) + (npoint(it)-1) * mact(isymu) + npoint(iu)
880   continue
c      if(odebug(78))print*,'pair ',it,iu,': l1,l2,l3 = ',l1,l2,l3
c...for vector computers: no recurrence in this loop
      do 890 n = 1 , np
         q(iws-1+(n-1)*listab+l1) = q(kpos1+n)
         q(iws-1+(n-1)*listab+l2) = q(kpos2+n)
         q(iws-1+(n-1)*listab+l3) = q(kpos3+n)
890   continue
      kpos1 = kpos1 + 3 * np
900   continue
c...output reordered integrals
      nw = nw + np*listab
      call outint (q(iw),nw,iblk,ned0)
      iws = iw + nw
910   ipair = ipair + np
920   continue
930   continue
c...write last block
      if (nw .gt. 0) call wrt3 (q(iw),nw,iblk,ned0)
      if (nw .gt. 0) print*,'last block: iblk = ',iblk,'     nw = ',nw
c...write dummy block
      if (iblk+1 .ne. iblkgv) call wrt3 (zero,1,iblk+1,ned0)
      if (odebug(5)) then
         k=ift
         do isyma=1,nirr
         print*,'active fock matrix, symmetry',isyma
         call outsqr (q(k),nsymm(isyma),nsymm(isyma),nsymm(isyma),' ')
         k=k+nsymm(isyma)**2
         end do
      end if
c...total fock matrix
      call daxpy(nfock,1.0d0,q(ifi),1,q(ift),1)
      iblofi = iblo7
      call wrt3 (q(ifi),nfock,iblofi,ned7)
      iblo7 = iblo7 + lensec(nfock)
      iblofa = iblo7
      call wrt3 (q(ift),nfock,iblofa,ned7)
      iblo7 = iblo7 + lensec(nfock)
      if (odebug(4)) then
         k = ifi
         do isyma=1,nirr
         print*,'inactive fock matrix, symmetry',isyma
         call outsqr (q(k),nsymm(isyma),nsymm(isyma),nsymm(isyma),' ')
         k = k + nsymm(isyma)**2
         end do
      end if
      if (odebug(9)) then
         k = ift
         do isyma=1,nirr
         print*,'total fock matrix, symmetry',isyma
         call outsqr (q(k),nsymm(isyma),nsymm(isyma),nsymm(isyma),' ')
         k = k + nsymm(isyma)**2
         end do
      end if
c...
      call corlsr (ibase)
      write (iwr,80000) cpulft(1)-time
80000 format ('integral reordering ',34('.'),f9.4,' seconds')
      return
      end
      subroutine invers (f,g)
c
c---------------------------------------------------------------------
c   Constructs the matrix  g(tuvx) = f(xvut)
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical(o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension f(*),g(*)

c   f(xvut) --> g(utxv)

      do 100 isyma = 1 , nirr
         isymb = mult(isyma,isymr)
         if (npair(isyma)*npair(isymb) .eq. 0) go to 100
         call trnsps (f(intof1(isyma)+1),g(intof1(isymb)+1),
     ?          npair(isymb),npair(isyma))
100   continue

c   g(utxv) --> f(utvx)

      do 250 iv = 1 , nact
         do 200 ix = 1 , nact
            isymvx = mult(itypea(iv),itypea(ix))
            isymtu = mult(isymvx,isymr)
            if (npair(isymtu) .eq. 0) go to 200
            call fmove (g(ioff1(ilifa(ix)+iv)+1),
     ?                  f(ioff1(ilifa(iv)+ix)+1),npair(isymtu))
200      continue
250   continue

c   f(utvx) --> g(tuvx)

      do 450 iv = 1 , nact
         do 400 ix = 1 , nact
            isymvx = mult(itypea(iv),itypea(ix))
            isymtu = mult(isymvx,isymr)
            if (npair(isymtu) .eq. 0) go to 400
            do 300 isyma = 1 , nirr
               isymb = mult(isyma,isymtu)
               if (mact(isyma)*mact(isymb) .eq. 0) go to 300
               call trnsps (
     +              f(ioff1(ilifa(iv)+ix)+1+iofsd(isymtu,isyma)),
     ?              g(ioff1(ilifa(iv)+ix)+1+iofsd(isymtu,isymb)),
     ?                      mact(isymb),mact(isyma))
300         continue
400      continue
450   continue

      return
      end
      subroutine lookup (symbol)
c
c----------------------------------------------------------------------
c   Constructs indexing arrays that are needed throughout the program
c   -- for example, the symmetries of the active orbitals, the indexing
c   of active orbital pairs (t,u) or the indexing of the all-active
c   two-electron integrals (tu|vx)
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
      character *1 digit(8)
      character *3 codes(10)
      character *3 symbol(14)
      character *8 word
      character*24 type(7)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      character *8 orbtag
      common /mcrun /orbtag(maxorb)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
c
      dimension icga(31),icgb(31)
c
      data digit/'1','2','3','4','5','6','7','8'/
      data codes/'fzc','fzv','cor','vir','doc','uoc','alp','bet',
     ?           'spe','end'/
      data type /'  frozen doubly occupied',
     ?           'inactive doubly occupied',
     ?           '  active doubly occupied',
     ?           '  active unoccupied     ',
     ?           '   alpha occupied       ',
     ?           '    beta occupied       ',
     ?           'inactive unoccupied     '/
      data m51 / 51 /
c
      multip = 1
c
c   read in symmetry information
c
      call secget(isect(490),m51,ibl490)
      nav = lenwrd()
      call readi(nirri,mach(13)*nav,ibl490,ned3)
      nirr = nirri
      nirs = max(0,nirr-2)
      do 1 i = 1 , nirr
         norbs(i)  = 0
         ifreez(i) = 0
         mcore(i)  = 0
         mact(i)   = 0
         nsymm(i)  = 0
1     continue
      do 2 i=1,maxorb
      if (isymai(i) .lt. 1) go to 3
2     nsymm(isymai(i)) = nsymm(isymai(i)) + 1
3     nbasis = i - 1
      nfreez = 0
      ncore = 0
      nact = 0

c   orbitals

      write (iwr,10)
10    format (//15x,28('*'),' orbital information ',28('*')/
     ?          15x,'*',75x,'*'/
     ?          15x,'*',5x,'orbitals',11x,'type of orbital',8x,
     ?                     'symmetry',5x,'label',10x,'*'/
     ?          15x,'*',5x,'--------',11x,'---------------',8x,
     ?                     '--------',5x,'-----',10x,'*')
      write (iwr,101)
101   format (15x,'*',75x,'*')
      itag = 0
4     itag = itag + 1      
      word = orbtag(itag)
      do 5 i=1,10
5     if (word(1:3) .eq. codes(i)) go to 6
      call caserr ('unrecognized orbital specification')
6     if (i .eq. 10) go to 20
      do 1105 k=1,nirr
1105  if (word(4:4) .eq. digit(k)) go to 1106
      call caserr ('unrecognized symmetry specification')
1106  if (itag .ge. nbegin) then
         norbs(k) = norbs(k) + 1
      end if
      go to (11,12,13,14,15,15,15,15,19,20), i

c   frozen core

11    nfreez = nfreez + 1
      ifreez(k) = ifreez(k) + 1
      ntype = 1
      write (iwr,102) itag,type(ntype),k
102   format (15x,'*',7x,i3,10x,a24,7x,i1,23x,'*')
      go to 4

12    call caserr ('fzv not possible in mclr')
      return

c   core

13    ncore = ncore + 1
      if (nfreez.gt.0 .and. ncore.eq.1) write (iwr,101)
      mcore(k) = mcore(k) + 1
      itype(itag-nfreez) = k
      npoint(itag-nfreez) = mcore(k)
      isymmo(itag-nfreez) = norbs(k)
      rocc(itag-nfreez) = 2.0d0
      ntype = 2
      write (iwr,103) itag,type(ntype),k,isymmo(itag-nfreez),
     ?                symbol(nirs+k)
103   format (15x,'*',7x,i3,10x,a24,7x,i1,7x,i2,a3,11x,'*')
      go to 4

14    call caserr ('vir not possible in mclr')
      return

c   active

15    nact = nact + 1
      if (ncore.gt.0 .and. nact.eq.1) write (iwr,101)
      mact(k) = mact(k) + 1
      itype(itag-nfreez) = k
      npoint(itag-nfreez) = mact(k)
      isymmo(itag-nfreez) = norbs(k)
      rocc(itag-nfreez) = 1.0d0
      if (i .eq. 5) rocc(itag-nfreez) = 2.0d0
      if (i .eq. 6) rocc(itag-nfreez) = 0.0d0
      ntype = i - 2
      ibfcod(nact) = i
      if (i .eq. 7) multip = multip + 1
      if (i .eq. 8) multip = multip - 1
      write (iwr,103) itag,type(ntype),k,isymmo(itag-nfreez),
     ?                symbol(nirs+k)
      go to 4

19    call caserr ('spe not possible in mclr')
      return

20    nba = nbasis - nfreez
      nprim = ncore + nact
      nvirt = nba - nprim
      do 21 j=1,nirr
      nsymm(j) = nsymm(j) - ifreez(j)
      nprm(j) = mcore(j) + mact(j)
21    mvirt(j) = nsymm(j) - nprm(j)
      write (iwr,101)
      k = nprim
      do 23 j=1,nirr
      if (mvirt(j) .eq. 0) go to 23
      if (mvirt(j) .eq. 1) then
         write (iwr,2111) nfreez+k+1,type(7),j,norbs(j)+1,symbol(nirs+j)
2111     format (15x,'*',7x,i3,10x,a24,7x,i1,7x,i2,a3,11x,'*')
      else
         write (iwr,2112) nfreez+k+1,nfreez+k+mvirt(j),
     ?                    type(7),j,norbs(j)+1,symbol(nirs+j),
     ?                    norbs(j)+mvirt(j),symbol(nirs+j)
2112     format (15x,'*',4x,i3,' to',i3,7x,a24,7x,i1,7x,i2,a3,
     ?           ' to ',i2,a3,'  *')
      end if
      do 22 iorb=1,mvirt(j)
      k = k + 1
      itype(k) = j
      npoint(k) = iorb
22    isymmo(k) = norbs(j) + iorb
23    continue
      write (iwr,101)
      write (iwr,231)
231   format (15x,77('*')//)

      do 25 i=1,nact
25    itypea(i) = itype(ncore+i)
      do 30 i=1,nba
30    ilifa(i) = (i-1)*nact

c   obtain numbers of electrons and space symmetry

      na = 0
      nb = 0
      isss = 1
      do 40 i=1,nact
      ii = ibfcod(i)
      if (ii.eq.5 .or. ii.eq.7) na = na + 1
      if (ii.eq.5 .or. ii.eq.8) nb = nb + 1
      if (ii.eq.7 .or .ii.eq.8) isss = mult(isss,itypea(i))
40    continue

c   define additional constants

      nprimp = nprim + 1
      nst    = ncore + 1
      nact2  = nact * nact
      nact3  = nact2 * nact
      nact4  = nact3 * nact
      nconac = ncore * nact
      nconco = ncore * ncore
      nvinco = nvirt * ncore
      nvinac = nvirt * nact
      nbasq  = nba * nba
      nbatri = (nba*(nba+1))/2
      intlen = na * nact * 2
      nstraa = ibinom(nact,na)
      nstrbb = ibinom(nact,nb)
      call vclr(rocc(nprimp),1,nvirt)

      write (iwr,200)
200   format ('important parameters:'/21('-')/)
      write (iwr,201) na,nb,multip,nirr,isss
201   format ('number of alpha electrons ............. ',i3/
     ?        'number of beta electrons .............. ',i3/
     ?        'spin multiplicity (2S+1) .............. ',i3//
     ?        'number of irreducible representations . ',i3/
     ?        'space symmetry of the ground state .... ',i3/)
      write (iwr,202) nba,nfreez,ncore,nact,nvirt
202   format ('number of basis functions ............. ',i3/
     ?        'number of frozen core orbitals ........ ',i3/
     ?        'number of core orbitals ............... ',i3/
     ?        'number of active orbitals ............. ',i3/
     ?        'number of virtual orbitals ............ ',i3/)
      write (iwr,710) (i,i=1,nirr)
710   format (//'core orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
720   format (56a1)
      write (iwr,730) (mcore(i),i=1,nirr)
730   format ('mcore  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,735) (i,i=1,nirr)
735   format (//'active orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,737) (mact(i),i=1,nirr)
737   format ('mact   |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,740) (i,i=1,nirr)
740   format (//'virtual orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,750) (mvirt(i),i=1,nirr)
750   format ('mvirt  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,760) (i,i=1,nirr)
760   format (//'all orbitals according to symmetries:'//
     ?        'irrep  |',8(i4,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
      write (iwr,770) (nsymm(i),i=1,nirr)
770   format ('nsymm  |',8(1x,i3,2x))
      write (iwr,720) ('-',i=1,8+6*nirr)
c      write (iwr,780) (i,i=1,nirr)
c780   format (//'dimension of the matrices A,B:'//
c     ?        'irrep  |',8(i4,2x))
c      write (iwr,720) ('-',i=1,8+6*nirr)
c      write (iwr,790) (icount(i),i=1,nirr)
c790   format ('icount |',8(i4,2x))
c      write (iwr,720) ('-',i=1,8+6*nirr)
c      write (iwr,800)
c800   format (//'the following roots are to be computed:'//
c     ?        'irrep        roots')
c      write (iwr,810)
c810   format ('---------------------')
c      do 830 i = 1 , nirr
c         if (nev(i) .gt. 0) write (iwr,820) i,nevlo(i),nevhi(i)
c820      format (i3,8x,i3,' to ',i3)
c830   continue
c      write (iwr,810)

c   determine length of the ci vectors of different symmetries

      call setsto(16,0,nstra)
      call setsto(16,0,lengci)
      mt = 0
110   call string (mt,nact,nb,icgb,itypea,jmt,iu,mult)
      if (mt .eq. 0) go to 120
      nstrb(jmt) = nstrb(jmt) + 1
      go to 110
120   continue
      mt = 0
130   call string (mt,nact,na,icga,itypea,jmt,iu,mult)
      if (mt .eq. 0) go to 140
      nstra(jmt) = nstra(jmt) + 1
      go to 130
140   continue
      do 150 isyma=1,nirr
      do 150 isymb=1,nirr
      isym = mult(isyma,isymb)
150   lengci(isym) = lengci(isym) + nstra(isyma)*nstrb(isymb)
      write (iwr,203) nstraa,nstrbb
203   format (//'number of intermediate alpha strings .. ',i6/
     ?        'number of intermediate beta strings ... ',i6//
     ?        'distributions of strings into symmetry:'/)
      write (iwr,205) (i,i=1,nirr)
205   format ('irrep   |',8(i6,2x))
      write (iwr,207) ('-',i=1,8+8*nirr)
207   format (72a1)
      write (iwr,209) 'nstra ',(nstra(i),i=1,nirr)
      write (iwr,209) 'nstrb ',(nstrb(i),i=1,nirr)
209   format (a6,2x,'|',8(i6,2x))
      write (iwr,207) ('-',i=1,8+8*nirr)
      write (iwr,211)
211   format (//'length of the ci vectors'/)
      write (iwr,205) (i,i=1,nirr)
      write (iwr,207) ('-',i=1,8+8*nirr)
      write (iwr,209) 'lengci',(lengci(i),i=1,nirr)
      write (iwr,207) ('-',i=1,8+8*nirr)
      write (iwr,213)
213   format (//)

      nci0 = lengci(isss)
      maxaa = 0
      maxbb = 0
      do 160 i = 1 , nirr
         maxaa = max(maxaa,nstra(i))
         maxbb = max(maxbb,nstrb(i))
160   continue

c   max no of singles from a string

      maxrpl = max(na*(nact+1-na),nb*(nact+1-nb))
      maxrp3 = maxrpl*3

c   addressing of active pairs (t,u)

      do 180 isym = 1,nirr
         k = 0
         do 169 isymt = 1 , nirr
            iofsd(isym,isymt) = k
            isymu = mult(isym,isymt)
            do 168 it = 1 , nact
               if (itypea(it) .ne. isymt) go to 168
               do 167 iu = 1 , nact
                  if (itypea(iu) .ne. isymu) go to 167
                  k = k + 1
                  ipos(ilifa(it)+iu) = k
167            continue
168         continue
169      continue
         npair(isym) = k
180   continue
      n0int1 = npair(1)
      if (odebug(38)) then
         write (iwr,1605)
1605     format (/'array ipos:'/)
         do 161 i=1,nact
161      write (iwr,162) (ipos(ilifa(i)+j),j=1,nact)
162      format (31i6)
      end if
      if (odebug(10)) write (iwr,163) (npair(i),i=1,8)
163   format (/'npair  = ',8i4)

c   set up addressing array ioff0 for symmetry isss

      call orbadr (isss,ioff0,n0int2,intof0,npair,nact,
     ?             itypea,mult,nirr)
      n0int = n0int1 + n0int2
      if (odebug(38)) then
          write (iwr,164) 
164      format (/'array ioff0:'/)
         do 165 i=1,nact
165      write (iwr,162) (ioff0(ilifa(i)+j),j=1,nact)
      end if

c   max no of orbital pairs

      maxpar = 0
      do 190 i = 1 , nirr
         maxpar = max(maxpar,npair(i))
190   continue

c   dimension of fock matrix

      nfock = 0
      do 310 i = 1 , nirr
         nfock = nfock + nsymm(i)**2
310   continue

      do 210 i = 1 , nact+1
         icf(i) = (i*(i-1))/2
210   continue

c   offsets and addresses for all-active integrals ZINT

      nint2 = 0
      do 240 iz = 1 , nirr
         intoff(iz) = nint2
         kt = 0
         ij = 0
         do 225 i = 1 , nact
            do 220 j = 1 , i
               ij = ij + 1
               if (iz .ne. mult(itypea(i),itypea(j))) go to 220
               kt = kt + 1
               ic3e(ij) = kt
220         continue
225      continue
         ij = 0
         do 235 i = 1 , nact
            do 230 j = 1 , i
               ij = ij + 1
               if (iz .ne. mult(itypea(i),itypea(j))) go to 230
               ic2e(ij) = nint2
               nint2 = nint2 + kt
230         continue
235      continue
         npairs(iz) = kt
240   continue
      nint1 = npairs(1)

      return
      end
      subroutine lralgo (nstrt,civec,d,x,s,e,sred1,ered1,bred,
     +      sred2,ered2,
     +      etare,etaim,sigma,xred,res,reso,resc,denom,eigen,msplit,
     +      iblb,ibls,ible,ibleig,iconv,q)
c
c----------------------------------------------------------------------
c   Solves the generalized linear response eigenvalue problem
c   (omega * S - E) x = 0  with the iterative MCLR algorithm
c   by Olsen et al. (J. Comp. Phys. 74, 265 (1988))
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
c
      integer ispec, itex, itable, inorms, ivect
      common /lrfile/ ispec,itex,itable,inorms,ivect
c
      common/lost  /dlost(100),ocvgd(50)
c
      dimension civec(*),d(*),x(*),s(*),e(*)
      dimension sred1(maxr,maxr),ered1(maxr,maxr),bred(maxr)
      dimension sred2(maxr2,maxr2),ered2(maxr2,maxr2)
      dimension etare(maxr2),etaim(maxr2),sigma(maxr2)
      dimension xred(maxr2,maxr2)
c      dimension resarr(neig),denom(neig),eigen(neig),q(*)
      dimension res(neig),reso(neig),resc(neig),denom(neig)
      dimension eigen(neig),q(*)
      dimension msplit(maxr),iblb(maxr),ibls(maxr),ible(maxr)
      dimension ibleig(neig),iconv(neig)
c
      data sigtol/0.00001d0/
      data ieigen,iresid / 20,21 /
c
      ntot1 = 1
      ntot = nstrt
10    niter = niter + 1
      if (niter.gt.maxit) return
c
c   fill up the reduced matrices S~,E~ (the compressed versions)
c
      do 40 i=ntot1,ntot
      do 40 j=i,1,-1
      n = leng(msplit(j))
      ix1 = 1
      ix2 = n + 1
      ih1 = noff(msplit(j))
      ih2 = ih1 + npar
      if (odebug(43)) write (iwr,111) j,iblb(j)
111   format ('mclr  : trying to read trial vector ',i2,
     +        ' from block ',i12)
      if (msplit(j)-2) 11,12,13
11    call rdedx (x,2*n,iblb(j),ned7)
      call bmove (x(ix2),x(io2),norb)
      call vclr(x(ic1),1,nci1)
      call vclr(x(ic2),1,nci1)
      go to 20
12    call rdedx (x,2*n,iblb(j),ned7)
      call bmove (x(ix2),x(ic2),nci1)
      call bmove (x(ix1),x(ic1),nci1)
      call vclr(x(io1),1,norb)
      call vclr(x(io2),1,norb)
      go to 20
13    call rdedx (x,2*n,iblb(j),ned7)
20    if (j .lt. i) go to 30
      time = cpulft(1)
      call lrmult (msplit(i),.false.,civec,d,x,s,e,q,q)
      tim(2) = tim(2) + cpulft(1) - time
      sred1(i,i) = ddot(n,x(ih1),1,s(ih1),1) + 
     +             ddot(n,x(ih2),1,s(ih2),1)
      ered1(i,i) = ddot(n,x(ih1),1,e(ih1),1) + 
     +             ddot(n,x(ih2),1,e(ih2),1)
      bred(i)    = ddot(n,x(ih1),1,e(ih2),1) + 
     +             ddot(n,x(ih2),1,e(ih1),1)
      if (odebug(30)) then
         print*,'sx and ex for trial vector ',i,':'
         call outvec (s,npar2,'sx')
         call outvec (e,npar2,'ex')
      end if
      if (odebug(71)) go to 40
      time = cpulft(1)
      if (msplit(i) .ne. 2) call trnsfn (s,2,2,npar,.true.,q)
      call trnsfn (e,2,2,npar,.true.,q)
      if (odebug(96)) then
         print*,'trial vector ',i,':'
         call outvec (x,npar2,'x  (with transformed orbital part)')
         call outvec (s,npar2,'sx (with transformed orbital part)')
         call outvec (e,npar2,'ex (with transformed orbital part)')
      end if
      shtim(1) = shtim(1) + cpulft(1) - time
      if (ostore) then
         call wrt3 (s,npar2,ibls(i),ned7)
         call wrt3 (e,npar2,ible(i),ned7)
      end if
      go to 40
30    ered1(j,i) = ddot(n,x(ih1),1,e(ih1),1) + 
     +             ddot(n,x(ih2),1,e(ih2),1)
      ered1(i,j) = ddot(n,x(ih1),1,e(ih2),1) + 
     +             ddot(n,x(ih2),1,e(ih1),1)
c      if (msplit(j) .eq. msplit(i)) then
         sred1(j,i) = ddot(n,x(ih1),1,s(ih1),1) + 
     +                ddot(n,x(ih2),1,s(ih2),1)
         sred1(i,j) = ddot(n,x(ih1),1,s(ih2),1) + 
     +                ddot(n,x(ih2),1,s(ih1),1)
40    continue
c
c   construct the reduced matrices S~,E~ (uncompressed versions)
c
      call redeig (maxr,ntot,sred1,ered1,bred,sred2,ered2)
      ntot2 = 2 * ntot
      if (odebug(17)) then
         print*,'sred2 & ered2, iteration step ',niter
         call outsqr (sred2,maxr2,ntot2,ntot2,'sred2')
         call outsqr (ered2,maxr2,ntot2,ntot2,'ered2')
      end if
c
c   solve reduced eigenvalue problem
c
      if (oreduc) go to 80
c
c   QZ algorithm
c
      time = cpulft(1)
      call qzhes (maxr2,ntot2,ered2,sred2,.true.,xred)
      call qzit (maxr2,ntot2,ered2,sred2,0.0d0,.true.,xred,ierr)
      call qzval (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,.true.,xred)
      if (ierr .ne. 0) go to 400
      call qzvec (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,xred)
      tim(3) = tim(3) + cpulft(1) - time
      if (odebug(18)) then
         do i=1,ntot2
         write (iwr,3666) i,niter
3666     format ('eigenvector ',i2,
     +           ' from QZ algorithm, iteration step',i2)
         call outvec (xred(1,i),ntot2,' ')
         end do
      end if
c
c   collect real positive eigenvalues
c
      nrealp = 0
      do 50 i=1,ntot2
      if (1.0d0+etaim(i).ne.1.0d0.or.dabs(sigma(i)).lt.sigtol) go to 50
      etare(i) = etare(i) / sigma(i)
      if (etare(i).lt.0.0d0) go to 50
      nrealp = nrealp + 1
      etare(nrealp) = etare(i)
      if (i.ne.nrealp) call dcopy(ntot2,xred(1,i),1,xred(1,nrealp),1)
50    continue
      if (nrealp .eq. 0) go to 500
      if (nrealp .lt. neig) then
         write(iwr,45) niter,nrealp
45       format ('warning: in iteration step ',i2,' there are only ',i2,
     ?           ' positive real eigenvalues.'/'trying to continue ...')
         call vfill(0.0d0,etare(nrealp+1),1,neig-nrealp)
      end if
      neig1 = min(nrealp,neig)
c
c   sort eigenvalues and eigenvectors
c
      do 70 i=1,neig
      k = i
      if (nrealp .gt. 1) then
         emin = etare(i)
         do 60 j=i+1,nrealp
         if (etare(j) .lt. emin) then
            k = j
            emin = etare(j)
         end if
60       continue
      end if
      t = etare(i)
      etare(i) = etare(k)
      etare(k) = t
      call dswap(ntot2,xred(1,i),1,xred(1,k),1)
70    continue
      go to 100
c
c   reduction to standard form
c
80    time = cpulft(1)
      call reduc (maxr2,ntot2,sred2,ered2,etaim,ierr)
      if (ierr .ne. 0) go to 600
      call tred2(maxr2,ntot2,sred2,etare,sigma,xred)
      call tql2(maxr2,ntot2,etare,sigma,xred,ierr)
      if (ierr .ne. 0) go to 700
      call rebak (maxr2,ntot2,ered2,etaim,ntot2,xred)
      tim(3) = tim(3) + cpulft(1) - time
      do 90 i=1,neig
      etare(i) = 1.0d0 / etare(ntot2+1-i)
90    call dcopy(ntot2,xred(1,ntot2+1-i),1,xred(1,i),1)
      neig1 = neig
c
c   output of eigenvalues
c
100   continue
c
c   new trial vectors to be expected
c
      ntot0 = ntot
      ntot1 = ntot + 1
      do 200 j=1,neig1
      if (odebug(54)) then
         print*,'x vector no ',j
         call outvec (xred(1,j),ntot2,' ')
      end if
      if (odebug(51)) print*,'iconv(',j,') = ',iconv(j)
      if (iconv(j) .eq. 1) then
c         resarr(j) = -1.0d0
         res(j)  = -1.0d0
         reso(j) = -1.0d0
         resc(j) = -1.0d0
         if (j .le. 50) ocvgd(j) = .false.
         go to 200
      end if
c
c   compute residue vector
c
      if (ostore) go to 110
      time = cpulft(1)
      call lrmult (3,.false.,civec,d,x,s,e,q,q)
      tim(2) = tim(2) + cpulft(1) - time
      if (odebug(47) .and. niter.le.2) then
         call outvec (s,npar2,'vector s*v')
         call outvec (e,npar2,'vector h*v')
      end if
      go to 150
c
c   store version
c
110   call vclr(s,1,npar2)
      call vclr(e,1,npar2)
      do 140 k=1,ntot0
      x1kj = xred(k      ,j)
      x2kj = xred(k+ntot0,j)
      call rdedx (x,npar2,ibls(k),ned7)
c...for vector computers: no recurrence in this loop
      do 120 i1 = 1 , npar
         i2 = i1 + npar
         s(i1) = s(i1) + x(i1) * x1kj - x(i2) * x2kj
         s(i2) = s(i2) + x(i2) * x1kj - x(i1) * x2kj
120   continue
      call rdedx (x,npar2,ible(k),ned7)
c...for vector computers: no recurrence in this loop
      do 130 i1 = 1 , npar
         i2 = i1 + npar
         e(i1) = e(i1) + x(i1) * x1kj + x(i2) * x2kj
         e(i2) = e(i2) + x(i2) * x1kj + x(i1) * x2kj
130   continue
140   continue
c
c   r = omega * s - e
c
150   call vsmsb(s,1,etare(j)-cshift,e,1,e,1,npar2)
c
c   residue vector now in e
c
      if (odebug(47)) then
         write (iwr,61) j
61       format ('residue vector no.',i3,':')
         call outvec (e,npar2,' ')
      end if
c
c   compute euclidean norm of residue vector
c
      ronorm = dnrsq(norb,e(io1),1) + dnrsq(norb,e(io2),1)
      rcnorm = dnrsq(nci1,e(ic1),1) + dnrsq(nci1,e(ic2),1)
      rnorm  = dsqrt(ronorm+rcnorm)
      ronorm = dsqrt(ronorm)
      rcnorm = dsqrt(rcnorm)
c      resarr(j) = rnorm
      res(j)  = rnorm
      reso(j) = ronorm
      resc(j) = rcnorm
      if (odebug(63)) write (iwr,753) j,ronorm,rcnorm,rnorm
753   format ('residue vector ',i2,':'/'ronorm = ',f14.9/'rcnorm = ',
     ?        f14.9/'rnorm  = ',f14.9)
c
c   threshold ?
c
      if (rnorm .gt. eps) go to 160
      iconv(j) = 1
      if (j .le. 50) ocvgd(j) = .true.
      eigen(j) = etare(j)
c
c   compute eigenvector v(j) (residue vector in e no longer needed)
c
      call blowup (ntot0,x,e,xred(1,j),msplit,iblb)
c
c   write vector v(j) onto external file
c
      call wrt3 (x,npar2,ibleig(j),ned7)
      if (odebug(90)) 
     +        print*,'output vector v(',j,') on blocks',ibleig(j),
     ?                       ' -',ibleig(j)+lensec(npar2)-1
      if (odebug(46)) then
         write (iwr,606) j,niter
606      format ('v vector no ',i2,', iteration step',i3)
         call outvec (x,npar2,' ')
      end if
      denom(j) = ddot(npar2,x,1,s,1)
      denomy(j) = ddot(npar,x,1,s,1)
      go to 200
160   if (j .le. 50) ocvgd(j) = .false.
      if (ntot.ge.maxr .or. ntot.ge.npar) go to 200
      ntot = ntot + 1
      if (ostore .or. odebug(71)) go to 170
      time = cpulft(1)
      call trnsfn (e,2,2,npar,.true.,q)
      shtim(1) = shtim(1) + cpulft(1) - time
170   if (onospl .or. niter.lt.nsplit) then
         msplit(ntot) = 3
      else if (ronorm .ge. rcnorm) then
         msplit(ntot) = 1
         if (odebug(63)) print*,'use orbital part'
      else
         msplit(ntot) = 2
         if (odebug(63)) print*,'use configuration part'
      end if
c      if (odebug(92)) then
         msave = msplit(ntot)
         msplit(ntot) = 3
c      end if
      n = leng(msplit(ntot))
      n2 = n + n
      ix1 = 1
      ix2 = n + 1
      ih1 = noff(msplit(ntot))
      ih2 = ih1 + npar
c
c   update residue vector (updated vector in e)
c
      if (ootva .and. msave.eq.1 .and. nsplit.eq.0 .and. ostore) then
         if (rcnorm .gt. eps) then
            thresh = tfac * rcnorm
         else
            thresh = dmax1(dsqrt(eps*eps-rcnorm*rcnorm),0.1d0*eps)
         end if
c         thresh = dmax1(eps,tfac*rcnorm)
         write (iwr,320) niter,j,rcnorm,thresh
320      format (/'starting optimal orbital trial vector algorithm in ',
     +      'iteration ',i2/'for root no ',i2,' --- rcnorm = ',f9.5,
     +      ' --- thresh = ',f9.5//' microiteration       residue'/
     +                   '---------------------------------')
         if = icorr(norb2)
         call optorb (etare(j),thresh,civec,d,x,s,e,q(if),q,iq,
     ?                msplit,ntot0,xred(1,j),ible)
         if (odebug(16)) then
c read v(j)
            call rdedx (x,npar2,ibleig(j),ned7)         
c replace orbital parts
            call dcopy(norb,q(if),1,x(io1),1)          
            call dcopy(norb,q(if+norb),1,x(io2),1)
c compute
            call lrmult (3,.false.,civec,d,x,s,e,q,q)   
c (omega * S - E) v~(j)
      call vsmsb(s,1,etare(j),e,1,e,1,npar2)
            t = dsqrt(dnrsq(norb,e(io1),1)+dnrsq(norb,e(io2),1))
            write (iwr,322) niter,j,t,thresh
322         format (/'iteration step ',i2,
     +               ': norm of orbital part of modified ',
     +              'residue vector ',i2,' is'/f14.9,
     +              ' (should be < ',f14.9,' )'/)
         end if
         call dcopy(norb,q(if),1,e(io1),1)
         call dcopy(norb,q(if+norb),1,e(io2),1)
         call vclr(e(ic1),1,nci1)
         call vclr(e(ic2),1,nci1)
         call corlsr (if)
      else
         call lrupda (3,etare(j),e(io1),e(io2),
     +                  x(io1),x(io2),e(io1),e(io2))
      end if
      if (odebug(48)) call outvec (e,npar2,'updated residue vector:')
c
c   orthogonalize against reference ci vector if necessary

      if (isym.eq.isss .and. msave.ge.2) then
         call orth (nci0,e(ic1),civec)
         call orth (nci0,e(ic2),civec)
      end if

c   gram-schmidt orthogonalisation

      do 180 k=1,ntot-1
      if (msplit(k)*msplit(ntot) .eq. 2) go to 180
      if (odebug(43)) write (iwr,111) k,iblb(k)
      call rdedx (x,2*leng(msplit(k)),iblb(k),ned7)
      m = min(msplit(k),msplit(ntot))
      n = leng(m)
      ix1 = 1
      if (msplit(ntot).eq.2 .and. msplit(k).eq.3) ix1 = ic1
      ix2 = ix1 + leng(msplit(k))
      ih1 = noff(m)
      ih2 = ih1 + npar
      t = - ddot(n,e(ih1),1,x(ix1),1) - 
     +      ddot(n,e(ih2),1,x(ix2),1)
      call daxpy(n,t,x(ix1),1,e(ih1),1)
      call daxpy(n,t,x(ix2),1,e(ih2),1)
      t = - ddot(n,e(ih1),1,x(ix2),1) - 
     +      ddot(n,e(ih2),1,x(ix1),1)
      call daxpy(n,t,x(ix2),1,e(ih1),1)
      call daxpy(n,t,x(ix1),1,e(ih2),1)
180   continue
      msplit(ntot) = msave
      if (odebug(92) .and. msave.lt.3) then
         if (msave .eq. 1) then
            call vclr(e(ic1),1,nci1)
            call vclr(e(ic2),1,nci1)
         else
            call vclr(e(io1),1,norb)
            call vclr(e(io2),1,norb)
         end if
         msplit(ntot) = 3
      end if
      n = leng(msplit(ntot))
      n2 = n + n
      ih1 = noff(msplit(ntot))
      ih2 = ih1 + npar
      t = dsqrt(dnrsq(n,e(ih1),1) + dnrsq(n,e(ih2),1) )
c      t = dnrm2(npar2,e,1)
      if (1.0d0+t .eq. 1.0d0) then
         write(iwr,81) j,niter
81       format ('discard new trial vector no. ',i2,' in iteration',
     ?           ' step ',i2/'because of linear dependence from',
     ?           ' the previous trial vectors.')
         ntot = ntot - 1
      else ! perform symmetric orthonormalisation
        if (odebug(49)) then
        call outvec (e(ih1),n,'first  part of vector to be "symorted":')
        call outvec (e(ih2),n,'second part of vector to be "symorted":')
        end if
         call symort (n,e(ih1),e(ih2),sytol,oldep)
         if (oldep) then
            write (iwr,82) j,niter
82          format ('discard new trial vector no. ',i2,' in',
     ?              ' iteration step ',i2/'because b and b#',
     ?              ' are (nearly) linearly dependent.')
            ntot = ntot - 1
         else ! write new trial vector onto external file
c           if (n .ne. npar) call dcopy(n,e(ih2),1,e(ih1+n),1)
c           The call above is fatally flawed as e(ih2) and e(ih1+n)
c           are partially aliased whereas compiling dcopy it is assumed
c           they are not.
            if (n .ne. npar) then
              do i1 = 0, n-1
                e(ih1+n+i1) = e(ih2+i1)
              enddo
            endif
c           iblost(j,niter+1) = iblo7
            call wrt3 (e(ih1),n2,iblb(ntot),ned7)
            if (odebug(90)) print*,'MCLR: output trial vector',ntot,
     ?                             ' on blocks',iblb(ntot),' - ',
     ?                             iblb(ntot)+lensec(n2)-1
            if (odebug(75)) then
               write (iwr,8200) j,niter
8200           format ('new trial vector ',i2,
     +                 ', iteration step ',i3,':')
               call outvec (e(ih1),n2,' ')
            end if
         end if
      end if
200   continue
      call lreign (itable,iwr,etare,ocvgd,neig,niter,ntot0,1)
      call residu (iwr,inorms,res,reso,resc,neig,niter,1)
      open (ieigen,file='roots',status='unknown')
      call lreign (itable,ieigen,dummy,dummy,neig,niter,ntot0,2)
      close (ieigen)
      open (iresid,file='residues',status='unknown')
      call residu (iresid,inorms,dummy,dummy,dummy,neig,niter,2)
      close (iresid)
      nconv = isum(neig,iconv,1)
      if (nconv .eq. neig) then

c   write last eigenvectors to external file in case some roots got lost

         do 204 j = 1 , neig1
            call vclr(s,1,npar2) ! to contain Sv
            call vclr(e,1,npar2) ! to contain v
            do 240 k = 1 , ntot0
               x1kj = xred(k      ,j)
               x2kj = xred(k+ntot0,j)
               call rdedx (x,npar2,ibls(k),ned7)
c...for vector computers: no recurrence in this loop
               do 220 i1 = 1 , npar
                  i2 = i1 + npar
                  s(i1) = s(i1) + x(i1) * x1kj - x(i2) * x2kj
                  s(i2) = s(i2) + x(i2) * x1kj - x(i1) * x2kj
220            continue
               n = leng(msplit(k))
               n1 = noff(msplit(k)) - 1
               n2 = n1 + npar
               call rdedx (x,2*n,iblb(k),ned7)
c...for vector computers: no recurrence in this loop
               do 230 i1 = 1 , n
                  i2 = i1 + n
                  e(n1+i1) = e(n1+i1) + x(i1) * x1kj + x(i2) * x2kj
                  e(n2+i1) = e(n2+i1) + x(i2) * x1kj + x(i1) * x2kj
230            continue
240         continue
            dlost(j) = ddot(npar2,e,1,s,1)
            call wrt3 (e,npar2,ible(j),ned7)
204      continue
         go to 1000
      end if

c   if at least one new trial vector has been added, start next iteration step

      if (ntot .ge. ntot1) go to 10
      if (ntot .ge. maxr) then

c   use approximate eigenvectors v(j) as new starting vectors

         do 205 j=1,neig1
            call blowup (ntot0,x,e,xred(1,j),msplit,iblb)
            call wrt3 (x,npar2,ible(j),ned7)
205      continue
         go to 1000
      end if

c   error messages

      write (iwr,210) niter
210   format (
     +  'no new trial vector has been added in iteration step ',i2,':'/
     +        'unable to finish iterative mclr calculation.')
      ierr = 1
      go to 1000

400   write (iwr,410) niter,ierr
410   format(/'convergence problems with qzit in iteration step ',i2,
     +        ': ierr = ',i2/)
      go to 1000

500   write (iwr,510) niter
510   format (/'in iteration step ',i2,
     +         ' there are no real eigenvalues.'/)
      ierr = 1
      go to 1000

600   write (iwr,610) niter,ierr
610   format (/'problem in iteration step ',i2,
     +         ': ered2 is not positive ',
     ?         'definite.'/'ierr = ',i4)
      go to 1000

700   write (iwr,710) niter,ierr
710   format (/'convergence problems with tql2g in iteration step',
     ?        i2,': ierr = ',i2/)
1000  continue

      return
      end
      subroutine lranal (civec,x,q,iq,energy,fe,iblock,iconv,
     ?                   symbol,zymbol,mode)
c
c----------------------------------------------------------------------
c   Analysis of the MCLR eigenvectors and output of tables to
c   stdout, TeX file and spectrum file
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), integer (i-n), logical (o)
      character*1   sgn,digit(10)
      character*2   cnum,char2i
      character*3   symbol,zymbol
      character*6   symtex(14),zymtex(14)
      character*80  outlin
      parameter(maxout=80)
      character*132 texlin
      parameter(maxtex=132)
c     maxtex is length of strings texlin
c     maxout is length of strings outlin
      logical btest
      real*8 occold,occnew
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
c
      integer ispec, itex, itable, inorms, ivect
      common /lrfile/ ispec,itex,itable,inorms,ivect
c
c
      dimension civec(*),x(*),q(*),iq(*),energy(*),
     ?          fe(*),iblock(*),iconv(*)
      dimension symbol(14),zymbol(14)
      dimension iocc(62),icga(31),icgb(31)
      dimension comp(4,20)
      dimension ihole(2),ipart(2)
      dimension kstraa(8),kstrba(8),kstrab(8),kstrbb(8)
      dimension ioff(8),ioffb(8)
c
      equivalence (sout(1),ihole(1)),(evalue(1),ipart(1))

      data symtex/'a''','a''''','a_1','b_1','b_2',
     ?            'a_2','a_g','b_{1u}','b_{2u}','b_{3g}',
     ?            'b_{3u}','b_{2g}','b_{1g}','a_u'/
      data zymtex/'A''','A''''','A_1','B_1','B_2',
     ?            'A_2','A_g','B_{1u}','B_{2u}','B_{3g}',
     ?            'B_{3u}','B_{2g}','B_{1g}','A_u'/
      data digit/'0','1','2','3','4','5','6','7','8','9'/
      data ev/27.21165d0/
      data pageth/45.0d0/
c
      ind(i,j) = ((i-1)*(i-2))/2 + j
c
      call vclr(comp,1,4*20)
      go to (1000,2000,3000), mode
c
1000  continue
c
c   initialize tex file
c
      write (itex,1010)
      write (itex,1020) 'MCLR'
      write (itex,1030)
      rlines = 6.3d0 
c
c   initialize spectrum file
c
      write (ispec,1040) multip,zymbol(nirs+isss),0.0d0,0.0d0
      if (nev(isss) .eq. 0) write (ispec,*)
      return
c
2000  continue
c
      time = cpulft(1)
      ncorem = (na + nb + 1)/2
c     thrcon = 0.0005d0
      epscau = 0.25d0 * epstab
      k = 0
      do 2010 isymb=1,nirr
      isyma = mult(isymb,isss)
      ioff(isymb) = k
2010  k = k + nstra(isyma) * nstrb(isymb)

c   write table of results

      write (itex,2020)
      rlines = rlines + 0.3d0

      write (iwr,2030) isym,zymbol(nirs+isym),epstab

      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t4 = 0.0d0
      t5 = 0.0d0
      t6 = 0.0d0
c     t7 = 0.0d0
c     t8 = 0.0d0
      do 240 j=1,neig
      if (iconv(j) .ne. 1) then
         write (iwr,2035) j
         write (iwr,3120)
         go to 240
      end if
      j1 = j
      if (isym .eq. isss) j1 = j1 + 1
      call rdedx (x,npar2,iblock(j),ned7)

c   determine y/z and o/c weights

      xnorsq = dnrsq(npar2,x,1)
      if (xnorsq .le. 1.0d-7) go to 2040
      if (j .gt. 20) go to 2040
      comp(1,j) = dnrsq(norb,x(io1),1) / xnorsq
      comp(2,j) = dnrsq(nci1,x(ic1),1) / xnorsq
      comp(3,j) = dnrsq(norb,x(io2),1) / xnorsq
      comp(4,j) = dnrsq(nci1,x(ic2),1) / xnorsq

2040  call trnsfn (x,1,1,npar,.true.,q)
      if (odebug(73)) write (iwr,12) j,j,denomy(j)
      if (dabs(denomy(j)) .lt. dsmall) then
         write (iwr,2050) j
         go to 240
      end if
      if (denomy(j) .lt. 0.0d0) then
         write (iwr,2060) j
         go to 240
      end if
      t = 1.0d0 / dsqrt(denomy(j))
      call dscal(npar,t,x,1)

      n = 0
      total = 0.0d0

c   E_ai|0>

      t0 = cpulft(1)
      do 2120 ia=nprimp,nba
      do 2120 ii=1,ncore
      if (mult(itype(ia),itype(ii)) .ne. isymr) go to 2120
      yai = x(ipo(ind(ia,ii)))
      if (dabs(yai) .lt. epstab) go to 2120

c   ATTENTION: next loop contains no overflow warning
c   which would prohibit vectorization

      do 2110 i=1,nci0
      t = yai * civec(i)
      if (dabs(t) .ge. epscau) then
         n = n + 1
         iii(n) = (1-ia) * nba - ii
         jjj(n) = i
         kkk(n) = 0
         gin(n) = t
         total = total + t * t
      end if
2110  continue
      if (n .gt. 340) go to 340
2120  continue
c     nai = n
      total = total + total
      t1 = t1 + cpulft(1) - t0

      ioa = icori(nstraa)
      iob = icori(nstrbb)
      ndim = max(na,nb)
      inter = icori(na*nact*2)
      call setsto(na*nact*2,0,iq(inter))
      call ciadrs (2,isss,iq(ioa),iq(iob),iq(inter),na)

c   E_ti|0>

      t0 = cpulft(1)
      nti = 0
      do 2390 isymi=1,nirr
      if (odebug(56)) print*,'entering i symmetry',isymi
      if (mcore(isymi) .eq. 0) go to 2390
      isymt = mult(isymr,isymi)
      if (mact (isymt) .eq. 0) go to 2390
      isymex = mult(isss,isymt)

c   lookup arrays for configurations with one additional alpha electron

      na = na + 1
      nstran = ibinom(nact,na)
      nstrbn = ibinom(nact,nb)
      iofaa = icori(nstran)
      iofba = icori(nstrbn)
      ndim = max(na,nb)
      intera = icori(ndim*nact*2)
      call setsto(ndim*nact*2,0,iq(intera))
      call ciadrs (2,isymex,iq(iofaa),iq(iofba),iq(intera),ndim)
      ncia = nci
      newa = icorr(ncia)
      newa1 = newa - 1
      do 2210 i=1,nirr
      kstraa(i) = kstra(i)
2210  kstrba(i) = kstrb(i)
      k = 0
      do 2220 i=1,nirr
2220  k = k + kstrba(i) * kstraa(mult(i,isymex))
      na = na - 1

c   lookup arrays for configurations with one additional beta electron

      nb = nb + 1
      nstran = ibinom(nact,na)
      nstrbn = ibinom(nact,nb)
      iofab = icori(nstran)
      iofbb = icori(nstrbn)
      ndim = max(na,nb)
      interb = icori(ndim*nact*2)
      call setsto(ndim*nact*2,0,iq(interb))
      call ciadrs (2,isymex,iq(iofab),iq(iofbb),iq(interb),ndim)
      ncib = nci
      newb = icorr(ncib)
      newb1 = newb - 1
      do 2230 i=1,nirr
      kstrab(i) = kstra(i)
2230  kstrbb(i) = kstrb(i)
      k = 0
      do 2240 i=1,nirr
      ioffb(i) = k
2240  k = k + kstrbb(i) * kstrab(mult(i,isymex))
      nb = nb - 1

      do 2320 ii=1,ncore
      if (itype(ii) .ne. isymi) go to 2320
      call vclr(q(newa),1,ncia)
      call vclr(q(newb),1,ncib)
      do 2310 it=1,nact
      if (odebug(56)) print*,'entering active orbitals'
      if (itypea(it) .ne. isymt) go to 2310
      xti = x(ipo(ind(it+ncore,ii)))

c   alpha replacements

      mta = 0
      iua = 0
      do 2270 ia=1,nstraa
      call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      do 2250 i=1,na
      if (icga(i) .eq. it) go to 2270
2250  continue
      icga(na+1) = it
      ipar = 1
      naold = na
      na = max(naold+1,nb)
      mtaa = istrad(ipar,iq(intera),icga,naold+1)
      na = naold
      isymb = mult(jmta,isss)
      iadold = iq(ioa-1+mta)
      iadnew = newa1 + iq(iofaa-1+mtaa)
      if (ipar .gt. 0) then
         do 2260 ib=1,nstrb(isymb)
2260     q(iadnew+ib) = q(iadnew+ib) + civec(iadold+ib) * xti
      else
         do 2265 ib=1,nstrb(isymb)
2265     q(iadnew+ib) = q(iadnew+ib) - civec(iadold+ib) * xti
      end if
      icga(na+1) = 0
2270  continue

c   beta replacements

      mtb = 0
      iub = 0
      do 2300 ib=1,nstrbb
      call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      do 2280 i=1,nb
      if (icgb(i) .eq. it) go to 2300
2280  continue
      icgb(nb+1) = it
      ipar = 1
      naold = na
      na = max(na,nb+1)
      mtbb = istrad(ipar,iq(interb+na*nact),icgb,nb+1)
      na = naold
      isyma = mult(jmtb,isss)
      isymb = mult(jmtb,isymt)
      nstrb0 = nstrb(jmtb)
      nstrb1 = kstrbb(isymb)
      iadold = ioff(jmtb) - nstrb0 + iq(iob-1+mtb)
      iadnew = newb1 + ioffb(isymb) - nstrb1 + iq(iofbb-1+mtbb)
      if (ipar .gt. 0) then
         do 2290 ia=1,nstra(isyma)
2290     q(iadnew+ia*nstrb1) = q(iadnew+ia*nstrb1) 
     ?                         + civec(iadold+ia*nstrb0) * xti
      else
         do 2295 ia=1,nstra(isyma)
2295     q(iadnew+ia*nstrb1) = q(iadnew+ia*nstrb1) 
     ?                         - civec(iadold+ia*nstrb0) * xti
      end if
2300  continue

c   end loop over active orbitals it

2310  continue

c   gather components .ge. threshold

      do 2350 k=1,ncia
      if (dabs(q(newa1+k)) .le. epscau) go to 2350
      n = n + 1
      nti = nti + 1
      if (n .gt. 340) go to 340
      call occdet (k,icga,icgb,na+1,nb,kstraa,kstrba,isymex)
c      if (odebug(56)) write (iwr,134) 
c     ?  ii,q(newa1+k),(icga(iaa),iaa=1,na+1),(icgb(ibb),ibb=1,nb)
c134   format ('orbital ',i2,':',f10.4,'(',<na+1>i2,' - ',<nb>i2,')')
      inta = 0
      intb = 0
      do 2330 iaa=1,na+1
2330  inta = ibset (inta,icga(iaa))
      do 2340 ibb=1,nb
2340  intb = ibset (intb,icgb(ibb))
      iii(n) = ii
      jjj(n) = inta
      kkk(n) = intb
      gin(n) = q(newa1+k)
      total = total + gin(n)**2
2350  continue
      do 2380 k=1,ncib
      if (dabs(q(newb1+k)) .le. epscau) go to 2380
      n = n + 1
      nti = nti + 1
      if (n .gt. 340) go to 340
      call occdet (k,icga,icgb,na,nb+1,kstrab,kstrbb,isymex)
c      if (odebug(56)) write (iwr,135) 
c     ?  ii,q(newb1+k),(icga(iaa),iaa=1,na),(icgb(ibb),ibb=1,nb+1)
c135   format ('orbital ',i2,':',f10.4,'(',<na>i2,' - ',<nb+1>i2,')')
      inta = 0
      intb = 0
      do 2360 iaa=1,na
2360  inta = ibset (inta,icga(iaa))
      do 2370 ibb=1,nb+1
2370  intb = ibset (intb,icgb(ibb))
      iii(n) = ii
      jjj(n) = inta
      kkk(n) = intb
      gin(n) = q(newb1+k)
      total = total + gin(n)**2
2380  continue

c   end loop over inactive orbitals ii

2320  continue

c   end loop over symmetries of inactive orbitals

      call corlsi (iofaa)
2390  continue

      t2 = t2 + cpulft(1) - t0

c   E_at|0>

      t0 = cpulft(1)
      nat = 0
      do 2710 isymc=1,nirr
      if (odebug(56)) print*,'entering a symmetry',isymc
      if (mvirt(isymc) .eq. 0) go to 2710
      isymt = mult(isymr,isymc)
      if (mact (isymt) .eq. 0) go to 2710
c      kvirt = 0
      isymex = mult(isss,isymt)

c   lookup arrays for configurations with one missing alpha electron

      na = na - 1
      nstran = ibinom(nact,na)
      nstrbn = ibinom(nact,nb)
      iofaa = icori(nstran)
      iofba = icori(nstrbn)
      ndim = max(na,nb)
      intera = icori(ndim*nact*2)
      call setsto(ndim*nact*2,0,iq(intera))
      call ciadrs (2,isymex,iq(iofaa),iq(iofba),iq(intera),ndim)
      ncia = nci
      if (odebug(56)) write (iwr,308) isymc,ncia
      newa = icorr(ncia)
      newa1 = newa - 1
      do 2490 i=1,nirr
      kstraa(i) = kstra(i)
2490  kstrba(i) = kstrb(i)
      k = 0
      do 2500 i=1,nirr
2500  k = k + kstrba(i) * kstraa(mult(i,isymex))
      na = na + 1

c   lookup arrays for configurations with one missing beta electron

      if (nb .eq. 0) go to 2525
      nb = nb - 1
      nstran = ibinom(nact,na)
      nstrbn = ibinom(nact,nb)
      iofab = icori(nstran)
      iofbb = icori(nstrbn)
      ndim = max(na,nb)
      interb = icori(ndim*nact*2)
      call setsto(ndim*nact*2,0,iq(interb))
      call ciadrs (2,isymex,iq(iofab),iq(iofbb),iq(interb),ndim)
      ncib = nci
      if (odebug(56)) write (iwr,309) isymc,ncib
      newb = icorr(ncib)
      newb1 = newb - 1
      do 2510 i=1,nirr
      kstrab(i) = kstra(i)
2510  kstrbb(i) = kstrb(i)
      k = 0
      do 2520 i=1,nirr
      ioffb(i) = k
2520  k = k + kstrbb(i) * kstrab(mult(i,isymex))
      nb = nb + 1
2525  continue

      do 2640 ic=nprimp,nba
      if (itype(ic) .ne. isymc) go to 2640
      if (odebug(56)) print*,'entering ic = ',ic
      call vclr(q(newa),1,ncia)
      call vclr(q(newb),1,ncib)
      do 2630 it=1,nact
      if (itypea(it) .ne. isymt) go to 2630
      if (odebug(56)) print*,'   entering it = ',it+ncore
      xai = x(ipo(ind(ic,it+ncore)))

c   alpha replacements

      mta = 0
      iua = 0
      do 2570 ia=1,nstraa
      call string (mta,nact,na,icga,itypea,jmta,iua,mult)
      do 2530 i=1,na
      if (icga(i) .eq. it) go to 2540
2530  continue
      go to 2570
2540  icga(i) = icga(na)
      ipar = 1
      naold = na
      na = max(naold-1,nb)
      mtaa = istrad(ipar,iq(intera),icga,naold-1)
      na = naold
      isymb = mult(jmta,isss)
      iadold = iq(ioa-1+mta)
      iadnew = newa1 + iq(iofaa-1+mtaa)
      if (ipar .gt. 0) then
         do 2550 ib=1,nstrb(isymb)
2550     q(iadnew+ib) = q(iadnew+ib) + civec(iadold+ib) * xai
      else
         do 2560 ib=1,nstrb(isymb)
2560     q(iadnew+ib) = q(iadnew+ib) - civec(iadold+ib) * xai
      end if
      icga(i) = it
2570  continue

c   beta replacements

      if (nb .eq. 0) go to 2625
      mtb = 0
      iub = 0
      do 2620 ib=1,nstrbb
      call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mult)
      do 2580 i=1,nb
      if (icgb(i) .eq. it) go to 2590
2580  continue
      go to 2620
2590  icgb(i) = icgb(nb)
      ipar = 1
      naold = na
      na = max(na,nb-1)
      mtbb = istrad(ipar,iq(interb+na*nact),icgb,nb-1)
      na = naold
      isyma = mult(jmtb,isss)
      isymb = mult(jmtb,isymt)
      nstrb0 = nstrb(jmtb)
      nstrb1 = kstrbb(isymb)
      iadold = ioff(jmtb) - nstrb0 + iq(iob-1+mtb)
      iadnew = newb1 + ioffb(isymb) - nstrb1 + iq(iofbb-1+mtbb)
      if (ipar .gt. 0) then
         do 2600 ia=1,nstra(isyma)
2600     q(iadnew+ia*nstrb1) = q(iadnew+ia*nstrb1)
     ?                         + civec(iadold+ia*nstrb0) * xai
      else
         do 2610 ia=1,nstra(isyma)
2610     q(iadnew+ia*nstrb1) = q(iadnew+ia*nstrb1)
     ?                         - civec(iadold+ia*nstrb0) * xai
      end if
      icgb(i) = it
2620  continue
2625  continue

c   end loop over active orbitals it

2630  continue

c   gather components .ge. threshold

      do 2670 k=1,ncia
      if (dabs(q(newa1+k)) .le. epscau) go to 2670
      n = n + 1
      nat = nat + 1
      if (n .gt. 340) go to 340
      call occdet (k,icga,icgb,na-1,nb,kstraa,kstrba,isymex)
c      if (odebug(56)) write (iwr,234) 
c     ?  ic,q(newa1+k),(icga(iaa),iaa=1,na-1),(icgb(ibb),ibb=1,nb)
c234   format ('orbital ',i2,':',f10.4,'(',<na-1>i2,' - ',<nb>i2,')')
      inta = 0
      intb = 0
      do 2650 iaa=1,na-1
2650  inta = ibset (inta,icga(iaa))
      do 2660 ibb=1,nb
2660  intb = ibset (intb,icgb(ibb))
      iii(n) = ic
      jjj(n) = inta
      kkk(n) = intb
      gin(n) = q(newa1+k)
      total = total + gin(n)**2
2670  continue
      do 2700 k=1,ncib
      if (dabs(q(newb1+k)) .le. epscau) go to 2700
      n = n + 1
      nat = nat + 1
      if (n .gt. 340) go to 340
      call occdet (k,icga,icgb,na,nb-1,kstrab,kstrbb,isymex)
c      if (odebug(56)) write (iwr,235) 
c     ?  ic,q(newb1+k),(icga(iaa),iaa=1,na),(icgb(ibb),ibb=1,nb-1)
c235   format ('orbital ',i2,':',f10.4,'(',<na>i2,' - ',<nb-1>i2,')')
      inta = 0
      intb = 0
      do 2680 iaa=1,na
2680  inta = ibset (inta,icga(iaa))
      do 2690 ibb=1,nb-1
2690  intb = ibset (intb,icgb(ibb))
      iii(n) = ic
      jjj(n) = inta
      kkk(n) = intb
      gin(n) = q(newb1+k)
      total = total + gin(n)**2
2700  continue

c   end loop over virtual orbitals ic

2640  continue

c   end loop over symmetries of virtual orbitals

      call corlsi (iofaa)
2710  continue

      t3 = t3 + cpulft(1) - t0

c   configuration components
 
c   ATTENTION: next loop contains no overflow warning
c   which would prohibit vectorization

      t0 = cpulft(1)
      do 2800 i=1,nci1
      if (dabs(x(norb+i)) .ge. epscau) then
         n = n + 1
         iii(n) = 0
         jjj(n) = i
         kkk(n) = 0
         gin(n) = x(norb+i)
         total = total + x(norb+i)**2
      end if
2800  continue
      t4 = t4 + cpulft(1) - t0
c      write (iwr,231) total
c231   format ('total = ',f14.9)

      go to 350
340   call caserr (
     +     'too many large components in an mclr eigenvector!')
350   continue

c   sort the non-spin-adapted components

      t0 = cpulft(1)
      do 500 i=1,n
      k = i
      if (n .gt. 1) then
         cmax = dabs(gin(i))
         do 400 i1=i+1,n
         if (dabs(gin(i1)) .gt. cmax) then
            k = i1
            cmax = dabs(gin(i1))
         end if
400      continue
      end if
      ghelp = gin(i)
      gin(i) = gin(k)
      gin(k) = ghelp
      ihelp = iii(i)
      iii(i) = iii(k)
      iii(k) = ihelp
      jhelp = jjj(i)
      jjj(i) = jjj(k)
      jjj(k) = jhelp
      khelp = kkk(i)
      kkk(i) = kkk(k)
      kkk(k) = khelp
500   continue
      if (gin(1) .lt. 0.0d0) call dscal(n,-1.0d0,gin,1)
      if (onoad) then
         nn = n
         go to 4501
      end if

c     tot = 0.0d0
      k = 1
      l = 1
4400  if (k .gt. n) go to 4405
      nhole = 0
      npart = 0
      nopsh = 0

c   evaluate excitation operator contribution if present

      if (iii(k) .ne. 0) then
         if (iii(k) .lt. 0) then
            ireal = - iii(k)
            ia = (ireal-1)/nba + 1
            ii = ireal - (ia-1)*nba
            nhole = nhole + 1
            ihole(nhole) = ii
         else
            if (iii(k) .le. ncore) then
               nhole = nhole + 1
               ihole(nhole) = iii(k)
            else
               ia = iii(k)
            end if
         end if
      end if

c   evaluate electron distribution in the active space

      call setsto(nact,0,iocc)
      if (iii(k) .le. 0) then
         if (iii(k) .lt. 0) then
            isym1 = isym
            isym = isss
         end if
         call detocc (jjj(k),icga,icgb)
         if (iii(k) .lt. 0) then
            isym = isym1
         end if
         do 280 iaa=1,na
280      iocc(icga(iaa)) = iocc(icga(iaa)) + 1
         do 290 ibb=1,nb
290      iocc(icgb(ibb)) = iocc(icgb(ibb)) + 1
      else
         do 270 it=1,nact
         if (btest(jjj(k),it)) iocc(it) = iocc(it) + 1
         if (btest(kkk(k),it)) iocc(it) = iocc(it) + 1
270      continue
      end if
      do 300 iorb=1,ncorem
      if (iocc(iorb) .lt. 2) then
         nhole = nhole + 1
         ihole(nhole) = ncore + iorb
         if (iocc(iorb) .lt. 1) ihole(nhole) = - ihole(nhole)
      end if
300   continue
      do 310 iorb=ncorem+1,nact
      if (iocc(iorb) .gt. 0) then
         npart = npart + 1
         ipart(npart) = ncore + iorb
         if (iocc(iorb) .gt. 1) ipart(npart) = - ipart(npart)
      end if
310   continue
      if (iii(k).lt.0 .or. iii(k).gt.ncore) then
         npart = npart + 1
         ipart(npart) = ia
      end if
      do 311 khole=1,nhole
311   if (ihole(khole) .gt. 0) nopsh = nopsh + 1
      do 312 kpart=1,npart
312   if (ipart(kpart) .gt. 0) nopsh = nopsh + 1
      if (nopsh .eq. 0) then
         factor = 1.0d0
      else
         factor = dsqrt(dble(nopsh))
      end if
      gin(l) = factor * gin(k)
      iii(l) = iii(k)
      jjj(l) = jjj(k)
      kkk(l) = kkk(k)
      nopshf = nopsh / 2
      if (nopsh.gt.0 .and. iii(k).lt.0) nopshf = nopshf - 1
      ipower = 2 ** nopshf
      k = k + ipower
      l = l + 1
      go to 4400

4405  continue
      t5 = t5 + cpulft(1) - t0
      t0 = cpulft(1)
      nn = l - 1

c   sort the spin-adapted components

      do 4500 i=1,nn
      k = i
      if (nn .gt. 1) then
         cmax = dabs(gin(i))
         do 4000 i1=i+1,nn
         if (dabs(gin(i1)) .gt. cmax) then
            k = i1
            cmax = dabs(gin(i1))
         end if
4000     continue
      end if
      ghelp = gin(i)
      gin(i) = gin(k)
      gin(k) = ghelp
      ihelp = iii(i)
      iii(i) = iii(k)
      iii(k) = ihelp
      jhelp = jjj(i)
      jjj(i) = jjj(k)
      jjj(k) = jhelp
      khelp = kkk(i)
      kkk(i) = kkk(k)
      kkk(k) = khelp
4500  continue
4501  continue

      do 4600 k=1,nn
      if (dabs(gin(k)) .lt. epstab) go to 4600
      nhole = 0
      npart = 0

c   evaluate excitation operator contribution if present

      if (iii(k) .ne. 0) then
         if (iii(k) .lt. 0) then
            ireal = - iii(k)
            ia = (ireal-1)/nba + 1
            ii = ireal - (ia-1)*nba
            nhole = nhole + 1
            ihole(nhole) = ii
         else
            if (iii(k) .le. ncore) then
               nhole = nhole + 1
               ihole(nhole) = iii(k)
            else
               ia = iii(k)
            end if
         end if
      end if

c   evaluate electron distribution in the active space

      call setsto(nact,0,iocc)
      if (iii(k) .le. 0) then
         if (iii(k) .lt. 0) then
            isym1 = isym
            isym = isss
         end if
         call detocc (jjj(k),icga,icgb)
         if (iii(k) .lt. 0) then
            isym = isym1
         end if
         do 380 iaa=1,na
380      iocc(icga(iaa)) = iocc(icga(iaa)) + 1
         do 390 ibb=1,nb
390      iocc(icgb(ibb)) = iocc(icgb(ibb)) + 1
      else
         do 370 it=1,nact
         if (btest(jjj(k),it)) iocc(it) = iocc(it) + 1
         if (btest(kkk(k),it)) iocc(it) = iocc(it) + 1
370      continue
      end if
      do 410 iorb=1,nact
      jorb = ncore + iorb
      occold = rocc(jorb)
      occnew = dble(iocc(iorb))
      if (occnew .lt. occold) then
         nhole = nhole + 1
         ihole(nhole) = jorb
         if (occnew .lt. occold-1.0d0) ihole(nhole) = - ihole(nhole)
      else if (occnew .gt. occold) then
         npart = npart + 1
         ipart(npart) = jorb
         if (occnew .gt. occold+1.0d0) ipart(npart) = - ipart(npart)
      end if
410   continue
      if (iii(k).lt.0 .or. iii(k).gt.ncore) then
         npart = npart + 1
         ipart(npart) = ia
      end if

c   construct line

      ico = 1
      ict = 1

c   Hartree-Fock configuration

      if (nhole*npart .eq. 0) then
         do 520 iorb=1,ncore+ncorem
         if (iorb .gt. 1) then
            texlin(ict:ict+1) = '\\:'
            ict = ict + 2
         end if
         cnum = char2i(isymmo(iorb))
         outlin(ico:ico+1) = cnum
         texlin(ict:ict+1) = cnum
         ico = ico + 2
         ict = ict + 2

c   protect a',a''

         if (nirs+itype(iorb) .le. 2) then
           texlin(ict:ict) = '{'
           ict = ict + 1
         end if
         outlin(ico:ico+2) = symbol(nirs+itype(iorb))
         texlin(ict:ict+5) = symtex(nirs+itype(iorb))
         ico = ico + 3
         ict = ict + 6
         if (nirs+itype(iorb) .le. 2) then
           texlin(ict:ict) = '}'
           ict = ict + 1
         end if
         outlin(ico:ico) = '('
         texlin(ict:ict) = '^'
         ico = ico + 1
         ict = ict + 1
         kocc = nint(rocc(iorb)) + 1
         outlin(ico:ico) = digit(kocc)
         texlin(ict:ict) = digit(kocc)
         ico = ico + 1
         ict = ict + 1
         outlin(ico:ico+1) = ') '
         ico = ico + 2
520      continue
         go to 565
      end if
      do 550 khole=1,nhole
      if (khole .gt. 1) then
         outlin(ico:ico+1) = ', '
         ico = ico + 2
         texlin(ict:ict+2) = ',\\,'
         ict = ict + 3
      end if
      cnum = char2i(isymmo(iabs(ihole(khole))))
      outlin(ico:ico+1) = cnum
      texlin(ict:ict+1) = cnum
      ico = ico + 2
      ict = ict + 2
      if (nirs+itype(iabs(ihole(khole))) .le. 2) then
         texlin(ict:ict) = '{'
         ict = ict + 1
      end if
      outlin(ico:ico+2) = symbol(nirs+itype(iabs(ihole(khole))))
      texlin(ict:ict+5) = symtex(nirs+itype(iabs(ihole(khole))))
      ico = ico + 3
      ict = ict + 6
      if (nirs+itype(iabs(ihole(khole))) .le. 2) then
         texlin(ict:ict) = '}'
         ict = ict + 1
      end if
      if (ihole(khole) .lt. 0) then
         outlin(ico:ico+2) = '(2)'
         ico = ico + 3
         texlin(ict:ict+1) = '^2'
         ict = ict + 2
      end if
550   continue
      outlin(ico:ico+4) = ' --> '
      ico = ico + 5
      texlin(ict:ict+11) = '\\rightarrow '
      ict = ict + 12
      do 560 kpart=1,npart
      if (kpart .gt. 1) then
         outlin(ico:ico+1) = ', '
         ico = ico + 2
         texlin(ict:ict+2) = ',\\,'
         ict = ict + 3
      end if
      cnum = char2i(isymmo(iabs(ipart(kpart))))
      outlin(ico:ico+1) = cnum
      texlin(ict:ict+1) = cnum
      ico = ico + 2
      ict = ict + 2
      if (nirs+itype(iabs(ipart(kpart))) .le. 2) then
         texlin(ict:ict) = '{'
         ict = ict + 1
      end if
      outlin(ico:ico+2) = symbol(nirs+itype(iabs(ipart(kpart))))
      texlin(ict:ict+5) = symtex(nirs+itype(iabs(ipart(kpart))))
      ico = ico + 3
      ict = ict + 6
      if (nirs+itype(iabs(ipart(kpart))) .le. 2) then
         texlin(ict:ict) = '}'
         ict = ict + 1
      end if
      if (ipart(kpart) .lt. 0) then
         outlin(ico:ico+2) = '(2)'
         ico = ico + 3
         texlin(ict:ict+1) = '^2'
         ict = ict + 2
      end if
560   continue
565   lenout = ico - 1
      lentex = ict - 1
      if (lentex.gt.maxtex) then
         write(iwr,*)'Maximum length of texlin is ',maxtex,' characters'
         write(iwr,*)'Current length of texlin is ',lentex,' characters'
         call caserr('lranal: lentex exceeds maxtex')
      endif
      if (lenout.gt.maxout) then
         write(iwr,*)'Maximum length of outlin is ',maxout,' characters'
         write(iwr,*)'Current length of outlin is ',lenout,' characters'
         call caserr('lranal: lenout exceeds maxout')
      endif
      if (k .eq. 1) then
         sgn = ' '
         write (iwr,750) j1,multip,zymbol(nirs+isym),energy(j)*ev,
     +   fe(j), sgn,dabs(gin(k)),(outlin(l:l),l=1,lenout)
         write (ispec,155) multip,zymbol(nirs+isym),energy(j),
     +   fe(j)
         write (itex,3110) j1,multip,zymtex(nirs+isym),
     +                     energy(j)*ev,fe(j),
     +                     dabs(gin(k)),(texlin(l:l),l=1,lentex)
         write (itex,3111)
         rlines = rlines + 1.0d0
      else
         sgn = '+'
         if (gin(k) .lt. 0.0d0) sgn = '-'
         write (iwr,760) sgn,dabs(gin(k)),(outlin(l:l),l=1,lenout)
         write (itex,4550) sgn,dabs(gin(k)),(texlin(l:l),l=1,lentex)
         write (itex,3111)
         rlines = rlines + 1.0d0
      end if
4600  continue
      if (nn .eq. 0) then
         write (iwr,751) j1,multip,zymbol(nirs+isym),energy(j)*ev,
     +                   fe(j)
         write (ispec,155) multip,zymbol(nirs+isym),energy(j),fe(j)
         write (itex,3112) j1,multip,zymtex(nirs+isym),
     ?                     energy(j)*ev,fe(j)
         rlines = rlines + 1.0d0
      end if
      write (iwr,3120)
      write (itex,2020)
      rlines = rlines + 0.3d0
      if (rlines .ge. pageth) then
         write (itex,4700)
         rlines = 3.3d0
      end if
      t6 = t6 + cpulft(1) - t0
240   continue
      write (ispec,*)
      write (iwr,359) t1,t2,t3,t4,t5,t6
      write (iwr,401)
      do 440 j=1,min(neig,20)
440   write (iwr,450) j,comp(1,j)+comp(2,j),comp(3,j)+comp(4,j),
     ?                  comp(1,j)+comp(3,j),comp(2,j)+comp(4,j),
     ?                  comp(1,j),comp(2,j),comp(3,j),comp(4,j)
      write (iwr,360) cpulft(1)-time
      return

c   closing statements

3000  continue
      write (itex,2020)
      rlines = rlines + 0.3d0
      write (itex,3010)
      return

2050  format ('eigenvector ',i3,' cannot be scaled.')
12    format ('<y(',i2,'),Sigma y(',i2,')> = ',f14.9)
2035  format (/'eigenvector ',i2,' has not converged and',
     ?            ' cannot be analysed.'/)
2030  format (//35('=')/' MCLR results for symmetry',i2,' (',a3,')'
     ?         /35('=')//' state       energy  oscillator',6x,
     ?         'leading configurations of first order perturbed',
     ?         ' wavefunction'/'              (eV)    strength',7x,
     ?         '(|c| > ',f5.2,')'/110('-'))
1040  format(4x,i1,a3,2f11.6)
2060  format ('scaling factor for eigenvector ',i3,' is negative !')
308   format ('isymc = ',i1,': ',i5,' new (na-1,nb) configurations')
309   format ('isymc = ',i1,': ',i5,' new (na,nb-1) configurations')
750   format (i3,'(',i1,')',a3,2f10.4,6x,a1,f5.2,6x,80a1)
760   format (35x,a1,f5.2,6x,80a1)
155   format (4x,i1,a3,2f11.6)
751   format (i3,'(',i1,')',a3,2f10.4)
3120  format (110('-'))
359   format ('a-i contributions: ',f9.4/
     ?        't-i contributions: ',f9.4/
     ?        'a-t contributions: ',f9.4/
     ?        'configurations:    ',f9.4//
     ?        '1st sorting:       ',f9.4/
     ?        '2nd sorting:       ',f9.4)
360   format (//'analysis of results took ',f16.9,' seconds'//)
450   format (i3,8f9.4)
401   format (//'analysis of orbital/configuration and y/z weights:'
     ?        //9x,'y',8x,'z',8x,'o',8x,'c',7x,
     ?            'yo',7x,'yc',7x,'zo',7x,'zc'/)
3112  format ('$ ',i3,'\\,^',i1,'\\!',a6,' $ & $ ',f5.2,
     ?        ' $ & $ ',f6.3,' $ & \\\\ ')
4550  format (' & & & $ ',a1,f4.2,' \\;(\\:',132a1)
4700  format ('\\end{tabular}'/'\\newpage'/
     ? '\\begin{tabular}{|c|c|c|l|} \\hline \\hline & & & \\\\'/
     ? 'State & $T_e$ [eV] & $f_e$ & $\\Psi\\approx$ \\\\'/
     ? ' & & & \\\\ \\hline \\hline')
1010  format ('\\documentstyle[11pt]{article}'/
     ?        '\\pagestyle{empty}'/
     ?        '\\topmargin=-1.8cm'/
     ?        '\\oddsidemargin=-0.5cm'/
     ?        '\\evensidemargin=-0.5cm'/
     ?        '\\textwidth=17cm'/
     ?        '\\textheight=26cm'/
     ?        '\\parindent=0cm'/
     ?        '\\begin{document}')
1020  format ('{\\Large Table of ',a4,' results for $\\ldots$}'/)
1030  format ('\\begin{tabular}{|c|c|c|l|} \\hline \\hline & & & \\\\'/
     ?        'State & $T_e$ [eV] & $f_e$ & $\\Psi\\approx$ \\\\'/
     ?        ' & & & \\\\ \\hline')
2020  format ('\\hline')
3110  format ('$ ',i3,'\\,^',i1,'\\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +        f6.3,' $ & $\\hspace{0.8em} ',f4.2,' \\;(\\:',132a1)
3111  format ('\\:) $ \\\\')
3010  format ('\\end{tabular}'/'\\end{document}')
      end
      subroutine lreign (itable,iw,eigen,ocvgd,neig,niter,ntot,mode)
c
c---------------------------------------------------------------------
c   MODE = 1 :  Write eigenvalues of each iteration to fortran 
c               file ITABLE
c   MODE = 2 :  Read eigenvalues from ITABLE and write them to
c               standard output in the correct order
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
      character*1 now(5)
      dimension eigen(*),ocvgd(*),buf(5),obuf(5)
      m = 1
      n = 5
      kk = 1
      go to (10,40), mode

c   write mode

10    nn = min(n,neig)
      write (itable) niter,ntot,(eigen(i),ocvgd(i),i=m,nn)
      if (neig .le. n) go to 30
      m = m + 5
      n = n + 5
      go to 10
30    return

c   read mode

40    nn = min(n,neig)
      write (iw,50) (i,i=m,nn)
50    format (15x,5i14)
      rewind itable
      do 80 i = 1 , niter
         do 60 j = 1 , kk-1
            read (itable)
60       continue
         read (itable) ibuf,jbuf,(buf(k),obuf(k),k=1,nn-m+1)
         do 65 k = 1 , nn-m+1
            if (obuf(k)) then
               now(k) = '*'
            else
               now(k) = ' '
            end if
65       continue
         write (iw,66) ibuf,jbuf,(buf(k),now(k),k=1,nn-m+1)
66       format (i5,i11,3x,f14.9,a1,4(f13.9,a1))
         do 70 j = kk+1 , (neig-1)/5+1
            read (itable)
70       continue
80    continue
      if (neig .le. n) return
      write (iw,90)
90    format (/)
      m = m + 5
      n = n + 5
      kk = kk + 1
      go to 40

      end
      subroutine lrinit 
c
c---------------------------------------------------------------------
c   Sets defaults for MCLR calculation and
c   reads in directives
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-w), logical (o)
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
c     mxtoken - the maximum number of tokens on an input line
c               for space separated tokens: 132/2 plus 1 to allow
c               for look ahead (132 is the line-length defined in
c               common/workc)
c
      integer mxtoken
      parameter(mxtoken=67)
      integer jrec, jump, istrt, inumb, iwidth, nend, nstart
      integer nline, noline, jwidth, nerr
      logical oswit, oterm, oflush
      common/work/jrec,jump,istrt(mxtoken),inumb(mxtoken),iwidth,
     + nend(mxtoken),nstart(mxtoken), nline,noline,jwidth,nerr,
     + oswit,oterm,oflush
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      character *8 orbtag
      common /mcrun /orbtag(maxorb)
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
      logical odebug
      common /debug/ odebug(100)
c
c
      real*8 freq, pola
      integer nfreq, ipola, mxfreq
      parameter (mxfreq=100)
      common /polars/ nfreq,ipola(3,3),freq(mxfreq),pola(6,mxfreq)
c
      dimension ymclr(50)
      data ymclr/'dirc','symm','maxr','maxi','eps ',   ! 01 - 05
     ?           'epsa','vecm','end ','debu','noca',   ! 06 - 10
     ?           'begi','nosp','pola','anal','redu',   ! 11 - 15
     ?           'true','exca','----','nozi','ipos',   ! 16 - 20
     ?           'reve','supp','inve','aos ','----',   ! 21 - 25
     ?           'sect','scf ','mcsc','cano','cive',   ! 26 - 30
     ?           'rpa ','orbi','file','----','safe',   ! 31 - 35
     ?           'thre','eige','eqsy','tabl','----',   ! 36 - 40
     ?           'olda','noan','coco','spli','ootv',   ! 41 - 45
     ?           'noad','dsma','shif','----','----'/   ! 46 - 50

c   set up defaults

      maxr = 50
      maxit = 30
      maxro = 20
      maxito = 30
      eps = 0.001d0
      epsana = 0.01d0
      epstab = 0.2d0
      tfac = 0.3d0
      nbegin = 1
      ocanon = .true.
      onospl = .true.
      oanal = .false.
      oreduc = .false.
      oexca = .false.
      ozint = .true.
      orevrs = .false.
      ostore = .true.
      osect = .false.
      orpa = .false.
      onewa = .true.
      onoana = .false.
      ococo = .false.
      ootva = .false.
      onoad = .false.
      dsmall = 1.0d-7
      oshift = .false.
      cshift = 0.0d0
      call setsto(8,0,nev)
      do 110 i = 1 , 100
         odebug(i) = .false.
110   continue

c   read in directives

100   call input
      call inpa4 (ytext)
      i = locatc (ymclr,50,ytext)
      if (i .eq. 0) go to 50
      go to ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     ?       11,12,13,14,15,16,17,18,19,20,
     ?       21,22,23,24,25,26,27,28,29,30,
     ?       31,32,33,34,35,36,37,38,39,40,
     ?       41,42,43,44,45,46,47,48      ),i
1     write (iwr,1001)
1001  format ('keyword dirc no longer necessary')
      go to 100
2     call inpi (nrsym)
      call inpi (nev(nrsym))
      go to 100
3     call inpi (maxr)
      go to 100
4     call inpi (maxit)
      go to 100
5     call inpf (eps)
      go to 100
6     call inpf (epsana)
      go to 100
7     call inpi (iseca)
      call inpi (isecb)
      go to 100
9     call inpi (idebug)
      odebug(idebug) = .true.
      go to 100
10    ocanon = .false.
      go to 100
11    call inpi (nbegin)
      go to 100
12    write (iwr,1205)
1205  format (/'nosplit keyword no longer necessary'/)
      go to 100
c
c   polarizabilities
c
13    continue
      if (nfreq .gt. mxfreq) go to 100
      do 13020 i = 2 , jump
         nfreq = nfreq + 1
         if (nfreq .le. mxfreq) then
            call inpf (freq(nfreq))
         else
            write (iwr,13010) mxfreq
13010       format (/
     +  'WARNING: Dynamic polarizabilities will be evaluated only for',
     +  ' the first ',i6,' specified frequencies.')
            go to 100
         end if
13020 continue
      go to 100
c
c   detailed analysis of eigenvectors

14    oanal = .true.
      go to 100
c
c   reduce small eigenvalue problems to standard form
c
15    oreduc = .true.
      go to 100
16    go to 100
17    oexca = .true.
      go to 100
18    go to 100
19    ozint = .false.
      go to 100
20    go to 100
21    orevrs = .true.
      go to 100
c
c   suppress storing of vectors Sb,Eb
c
22    ostore = .false.
      go to 100
c
23    go to 100
24    call inpi (isecao)
      go to 100
25    go to 100
26    osect = .true.
      go to 100
27    if (osect) then
         call inpi (isecao)
      else
         write (iwr,60)
         go to 50
      end if
      go to 100
28    if (osect) then
         call inpi (isecb)
      else
         write (iwr,60)
         go to 50
      end if
      go to 100
29    if (osect) then
         call inpi (iseccn)
      else
         write (iwr,60)
         go to 50
      end if
      go to 100
30    if (osect) then 
         call inpi (iseca)
      else
         write (iwr,60)
         go to 50
      end if
      go to 100
31    orpa = .true.
      go to 100
c
c   orbitals
c
32    iorbit = 1
      oslash = .false.
292   call input
      do 291 i = 1 , jump
         call inpa (orbtag(iorbit))
         if (index(orbtag(iorbit),'/') .ne. 0) oslash = .true.
         if (orbtag(iorbit) .eq. 'end') then
            if (oslash) then
               oslash = .false.
               go to 291
            else
               if (iorbit .eq. 1) 
     ?            call caserr (
     +            'You forgot to specify the active orbitals !')
               go to 100
            end if
         end if
         iorbit = iorbit + 1
291   continue
      go to 292
33    m6file = 1
      call inpi (m6tape(1))
      call inpi (m6blk(1))
      call inpi (m6last(1))
      go to 100
34    go to 100
35    go to 100
c
c   thresholds
c
36    othrsh = .true.
      go to 100
c
c   ... for iterative eigenvalue search
c
37    if (othrsh) then
         call inpf (eps)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
c   ... for polarizability calculation
c
38    if (othrsh) then
         call inpf (epss)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
c   ... for table of results
c
39    if (othrsh) then
         call inpf (epstab)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
40    go to 100
c
41    onewa = .false.
      go to 100
c
42    onoana = .true.
      go to 100
c
43    ococo = .true.
      go to 100
c
44    onospl = .false.
      call inpi (nsplit)
      go to 100 
c
45    ootva = .true.
      if (jump .ge. 2) call inpi (maxro)
      if (jump .ge. 3) call inpf (tfac)
      go to 100
c
46    onoad = .true.
      go to 100
c
47    call inpf (dsmall)
      go to 100
c
48    oshift = .true.
      call inpf (cshift)
      go to 100
c
8     return
50    call caserr (
     +         'unrecognised directive or faulty directive ordering')
60    format (/'no section directive was presented'/)
      return
      end


      block data lr_data
      implicit real*8 (a-h,p-z), logical (o)
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
      integer ispec, itex, itable, inorms, ivect
      common /lrfile/ ispec,itex,itable,inorms,ivect
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      data mult/1,2,3,4,5,6,7,8,
     1          2,1,4,3,6,5,8,7,
     2          3,4,1,2,7,8,5,6,
     3          4,3,2,1,8,7,6,5,
     4          5,6,7,8,1,2,3,4,
     5          6,5,8,7,2,1,4,3,
     6          7,8,5,6,3,4,1,2,
     7          8,7,6,5,4,3,2,1/
      data ispec,itex,itable,inorms,ivect/11,12,14,15,16/
      data ned0,ned3,ned4,ned6,ned7,sytol/1,4,5,7,8,1.0d-12/
      end

      subroutine lrmain (q,iq,nmaxly)
c
c---------------------------------------------------------------------
c   Driving routine for MCLR calculation
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
      character *3 symbol(14)
      character *3 zymbol(14)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      real*8 freq, pola
      integer nfreq, ipola, mxfreq
      parameter (mxfreq=100)
      common /polars/ nfreq,ipola(3,3),freq(mxfreq),pola(6,mxfreq)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical odebug
      common /debug/ odebug(100)
c
c
      integer ispec, itex, itable, inorms, ivect
      common /lrfile/ ispec,itex,itable,inorms,ivect
c
c
      integer jjfile, notape, iblk, lblk
      integer nnfile, nofile, jblk, mblk
      integer mmfile, nufile, kblk, nblk
      integer ione, lone, lds, isect, ldsect, iacsct
      common/restri/jjfile,notape(4),iblk(4),lblk(4),
     +              nnfile,nofile(4),jblk(4),mblk(4),
     +              mmfile,nufile(4),kblk(4),nblk(4),
     +              ione(12),lone(12),
     +              lds(508),isect(508),ldsect(508),iacsct(508)
c
      common/lost  /dlost(100),ocvgd(50)
      common/statis/begin,ebegin
c
      dimension q(*),iq(*)
c
      data m511/511/
      data tolmy/1.0d-9/
      data symbol/'a''','a"','a1','b1','b2','a2',
     ?            'ag','b1u','b2u','b3g','b3u','b2g','b1g','au'/
      data zymbol/'A''','A"','A1','B1','B2','A2',
     ?            'Ag','B1u','B2u','B3g','B3u','B2g','B1g','Au'/

      ind(i,j) = ((i-1)*(i-2)) / 2 + j
c
      write (iwr,4100)
4100  format(///19x,83('-')/19x,'-',81x,'-'/19x,
     ? '-    MULTICONFIGURATIONAL LINEAR RESPONSE CALCULATION',29x,'-'/
     ? 19x,'-',81x,'-'/19x,83('-'))
      ipola(1,1) = 1
      ipola(2,1) = 2
      ipola(3,1) = 3
      ipola(2,1) = 2
      ipola(2,2) = 4
      ipola(2,3) = 5
      ipola(3,1) = 3
      ipola(3,2) = 5
      ipola(3,3) = 6
c
      call cpuwal(begin,ebegin)
c
      i = inicor(nmaxly)
      write (iwr,101) i,intrel
101   format (/'space available:',i12,' reals'/
     ?         'intrel = ',i2/)
      time = cpulft(1)
c
c   some initializations
c
      iblo7 = 1
      fesum = 0.0d0
      write (iwr,910)
      if (nfreq .gt. 0) call vclr(pola(1,1),1,6*nfreq)
c
c   open fortran files
c
      open (ispec,file='spectrum',status='unknown')
      open (itex ,file='table.tex',status='unknown')
      open (itable,form='unformatted',status='unknown')
      open (inorms,form='unformatted',status='unknown')
c
c   get dipole integrals
c
      call getdip (q)
c
c   get lookups
c
      call lookup (symbol)
      if (orpa .or. odebug(87)) write (iwr,87)
87    format (/'sorry -- getuni not available on HP because ',
     ?         ' routines DGEFA and DGESL are missing'/)
c
c   initialize tex and spectrum file
c
      call lranal (dummy,dummy,dummy,dummy,
     ?             dummy,dummy,dummy,dummy,
     ?             symbol,zymbol,1)
c
      maxr2 = maxr + maxr
      maxrs = maxr * maxr
      maxr2s = maxr2 * maxr2
      c23 = 2.0d0 / 3.0d0
c
c   determine bounds of scratch space
c
      call mclr_space
c
c   read in total energy
c
      call secget (iseca,172,iblkd)
      call rdedx (gin,m511,iblkd,ned3)
      e0 = gin(14)
      if (odebug(12)) print*,'e0 = ',e0
c
c   read in ci vector
c
      idens = icorr(nact2)
      idens1 = icorr(n0int1)
      icivec = icorr(nci0)
      iblkc = iblkd + lensec(mach(18))
      call rdedx (q(icivec),nci0,iblkc,ned3)
      if (odebug(2)) call outvec (q(icivec),nci0,'ci vector:')
c
c   construct density matrices
c
      t = cpulft(1)
      ig     = icorr(n0int)
      ip2    = icorr(n0int2)
      ip3    = icorr(n0int2)
      iofa   = icori(nstraa)
      iofb   = icori(nstrbb)
      inter  = icori(intlen)
      call setsto(intlen,0,iq(inter))
      call ciadrs (2,isss,iq(iofa),iq(iofb),iq(inter),na)
      nablk = (maxaa/2)*2 + 1
      nbblk = (maxbb/2)*2 + 1
      mmax = icorrm()-1
240   istore = (maxrp3+1)*(nablk+nbblk) + maxpar*nablk*nbblk
      if (istore .le. mmax) go to 260
      if (nablk .gt. nbblk) go to 250
      nbblk = nbblk - 2
      go to 240
250   nablk = nablk - 2
      go to 240
260   nwa = icori(nablk)
      nwb = icori(nbblk)
      iwa = icori(maxrp3*nablk)
      iwb = icori(maxrp3*nbblk)
      id = icorr(maxpar*nablk*nbblk)
      call denmat (q(icivec),q(ig),q(ip2),q(ip3),
     +             iq(iofa),iq(iofb),iq(inter),
     ?             iq(iwa),iq(iwb),q(id),nablk,nbblk,iq(nwa),iq(nwb))
      write (iwr,267) cpulft(1)-t
267   format ('construction of density matrices ',21('.'),
     +         f9.4,' seconds')
      ibdp = iblo7
      call wrt3 (q(ig),n0int,ibdp,ned7)
      if (odebug(90)) print*,'output density matrices on blocks',
     ?                ibdp,' -',ibdp+lensec(n0int)-1
      iblo7 = iblo7 + lensec(n0int)
      iblkpp = iblo7
      call wrt3 (q(ip2),2*n0int2,iblkpp,ned7)
      if (odebug(90)) print*,'output matrices p2 & p3 on blocks',
     ?                iblkpp,' -',iblkpp+lensec(2*n0int2)-1
      iblo7 = iblo7 + lensec(2*n0int2)
      if (odebug(6)) call outden (1,1,ioff0,q(ig+n0int2),q(ig))
      call vclr(q(idens),1,nact2)
      call fmove (q(ig+n0int2),q(idens1),n0int1)
      igoff = ig - 1 + n0int2
      do 2000 itu = 1 , nact2
         it = (itu-1)/nact + 1
         iu = itu-ilifa(it)
         if (itypea(it).eq.itypea(iu)) 
     +                           q(idens-1+itu) = q(igoff+ipos(itu))
2000  continue
      do 2010 it = 1 , nact
         ddiag(it) = q(idens-1+ilifa(it)+it)
2010  continue
      call corlsr (ig)

c   perform integral reordering

      call intord (q,q(idens))
      call lrzint (q)
      call rmatrx (q)
      ifock = icorr(nfock)
      call rdedx (q(ifock),nfock,iblofa,ned7)
      if (ocanon) call canoni (q(ifock),q(idens1),q)
      call corlsr (ifock)
      if (nfreez .eq. 0) then
         const = core - e0
      else                         !  compute <c,Hc>
         nlen = nci0
         isym = isss
         izint = icorr(nint1+nint2)
         isol = icorr(nci0)
         call rdedx (q(izint),nint1+nint2,ibzint,ned7)
         iofa = icori(nstraa)
         iofb = icori(nstrbb)
         inter = icori(intlen)
         call setsto(intlen,0,iq(inter))
         call ciadrs (2,isss,iq(iofa),iq(iofb),iq(inter),na)
         nablk = (maxaa/2)*2+1
         nbblk = (maxbb/2)*2+1
         mmax = icorrm()-1
940      istore = (maxrp3+1)*(nablk+nbblk) + 2*maxpar*nablk*nbblk
         if (istore .le. mmax) go to 960
         if (nablk .gt. nbblk) go to 950
         nbblk = nbblk - 2
         go to 940
950      nablk = nablk - 2
         go to 940
960      nwa = icori(nablk)
         nwb = icori(nbblk)
         iwa = icori(maxrp3*nablk)
         iwb = icori(maxrp3*nbblk)
         id = icorr(maxpar*nablk*nbblk)
         ie = icorr(maxpar*nablk*nbblk)
         call vclr(q(isol),1,nci0)
         call cislow (q(icivec),dummy,q(isol),q(izint),dummy,
     +                1,1,0,0,ic3e,iq(iofa),iq(iofb),iq(inter),
     +                iq(iwa),iq(iwb),q(id),q(ie),
     ?                nablk,nbblk,iq(nwa),iq(nwb))
         const = - ddot(nci0,q(icivec),1,q(isol),1)
         call corlsr (izint)
      end if
      if (odebug(12)) print*,'const = ',const
      iblo7s = iblo7
      if (orpa) rewind ivect

c   perform symmetry blocked linear response calculation

      do 1000 isym = 1 , nirr
         lrstar = icorr(0)
         if (nev(isym).eq.0 .and. nfreq.eq.0) go to 1000
         lrspce = icorrm()-1
         write (iwr,910)
         call flushn(iwr)
         iblo7 = iblo7s
         isymr = mult(isym,isss)
         n1int1 = npair(isymr)

c   addresses of orbital parameters

         k = 0
         do 51 isymt = 1 , nirr
            isymi = mult(isymt,isymr)
            iorbti(isymt) = k
            do 50 it = nst , nprim
               if (itype(it) .eq. isymt) then
                  do 49 ii = 1 , ncore
                     if (itype(ii) .eq. isymi) then
                        k = k + 1
                        ipo(ind(it,ii)) = k
                     end if
49             continue
               end if
50          continue
51       continue
         norbti = k
         do 61 isyma = 1 , nirr
            isymi = mult(isyma,isymr)
            iorbai(isyma) = k
            do 60 ia = nprimp , nba
               if (itype(ia) .eq. isyma) then
                  do 59 ii = 1 , ncore
                     if (itype(ii) .eq. isymi) then
                        k = k + 1
                        ipo(ind(ia,ii)) = k
                     end if
59                continue
               end if
60          continue
61       continue
         norbai = k - norbti
         do 71 isyma = 1 , nirr
            isymt = mult(isyma,isymr)
            iorbat(isyma) = k
            do 70 ia = nprimp , nba
               if (itype(ia) .eq. isyma) then
                  do 69 it = nst , nprim
                     if (itype(it) .eq. isymt) then
                        k = k + 1
                        ipo(ind(ia,it)) = k
                     end if
69                continue
               end if
70          continue
71       continue
         norbat = k - norbti - norbai
         norb = k

c   number of parameters

         norb2 = norb + norb
         nci1 = lengci(isym)
         npar = norb + nci1
         npar2 = npar + npar
         nconf = numcsf(nact,na,nb,isym,itypea,mult)
         if (isym .eq. isss) nconf = nconf - 1
         nparam = norb + nconf
         write (iwr,310) isym,zymbol(nirs+isym),norb,nci1,
     +                   nconf,npar2,lrstar
310      format (/81('*')/16x,
     +           'LINEAR RESPONSE CALCULATION FOR SYMMETRY',
     ?           i3,' (',a3,')'/81('*')//
     ?           'number of excitation operators: ',i7/
     ?           'length of ci vector:            ',i7/
     ?           'number of CSFs:                 ',i7/
     ?           'total number of parameters:     ',i7//
     ?           'pointer of scratch space is on  ',i7/)
         leng(1) = norb
         leng(2) = nci1
         leng(3) = npar
         noff(1) = 1
         noff(2) = norb + 1
         noff(3) = 1
         io1 = 1
         io2 = npar + 1
         ic1 = norb + 1
         ic2 = io2 + norb
         neig = min(nev(isym),nparam)
         if (odebug(60)) then
            neig = npar
            if (isym .eq. isss) neig = npar - 1
         end if
         if (odebug(65)) neig = nparam

c   set up addressing array ioff1 for symmetry isymr

         call orbadr (isymr,ioff1,n1int2,intof1,npair,nact,
     +                itypea,mult,nirr)
         if (odebug(40)) then
            write (iwr,264) isym
264         format (/'array ioff1 for symmetry ',i1,':'/)
            do 265 i=1,nact
265         write (iwr,162) (ioff1(ilifa(i)+j),j=1,nact)
162         format (31i6)
         end if
         n1int = n1int1 + n1int2
         iconv = icori(neig)
         iblb = icori(maxr)
         ibls = icori(maxr)
         ible = icori(maxr)
c         iblbb = icori(maxr)
c         iblss = icori(maxr)
c         iblee = icori(maxr)
         imspli = icori(maxr)

c   space for commutators

         ibegin = (neig-1)/intrel+1 + 4 * ((maxr-1)/intrel+1)
         mclrsp = ibegin + 3 * maxr2 + (neig-1)/intrel+1
     +                   + 3 * npar2 + 2 * maxrs + maxr 
     +                   + 3 * maxr2s + 5 * neig
         is1 = 0
         is2 = 0
         do 79 isym1=1,nirr
            isym2 = mult(isym1,isymr)
            m1 = (mact (isym1) * mcore(isym2)) ** 2
            m2 = (mvirt(isym1) * mcore(isym2)) ** 2
            m3 = (mvirt(isym1) * mact (isym2)) ** 2
            is1 = max(is1,m1,m2,m3)
            is2 = is2 + m1 + m2 + m3
79       continue
         kspace = lrspce - max(mclrsp,ibegin+nfock+norb+
     +                          2*is1+is2+ntrans)
         if (kspace .lt. norb2) then
            write (iwr,7901) lrspce,kspace,norb2
7901        format('memory problems:'/
     ?             'space available for MCLR ......... ',i9/
     ?             'space available for o-o .......... ',i9/
     ?             'norb2 ............................ ',i9)
            call caserr ('insufficient space for o-o handling')
         else
            lp = 1
            if (norb .gt. 0) lp = min(norb,kspace/norb2)
         end if

c   construction of A(oo) and B(oo): 
c   the double commutators  <0|[E_pq,[H_0,E_rs]]|0>

         ifi = icorr(nfock)
         idiag = icorr(norb)
         if (ncore .gt. 0) call wmatrx (q)
         ip1 = icorr(n0int)
         ip2 = icorr(n0int2)
         ip3 = icorr(n0int2)
         ift = icorr(nfock)
         icomm = icorr(2*(nbasq+204))
         iio = icori(nbasq+204)
         ijo = icori(nbasq+204)
         call rdedx (q(ip1),n0int,ibdp,ned7)
         call rdedx (q(ip2),2*n0int2,iblkpp,ned7)
         call oopart (q(idens1),q(ip1),q(ip2),q(ip3),q(ifi),q(ift),
     ?        q(icomm),iq(iio),iq(ijo),q)
         call corlsr (ip1)

c   reordering of the double commutators

         is = icorr(0)
         call cosort (q(is))

c   transform diagonal of A(oo)

         if (ocanon) then
            ia1 = icorr(is1)
            ia2 = icorr(is1)
            ia3 = icorr(is2)
            iuu = icorr(ntrans)
            iqq = icorr(0)
            if (odebug(94)) print*,'pointer is on ',iqq,
     +                    ' before routine DIA'
            call dia (q(ia1),q(ia2),q(ia3),q(iuu),q(idiag),q(iqq))
         end if

c   construct diagonal of matrix A

         call corlsr (idiag)
         idiaga = icorr(npar)
         idiagt = icorr(nci1)
         izint = icorr(n0int)
         izb   = icorr(maxbb*nact)
         iz    = icorr(maxbb*nact)
         call adiagc (q(idiaga+norb),q(idiagt),q(izint),q(izb),
     +                q(iz),q(ifi))
         ibdiag = iblo7
         call wrt3 (q(idiaga),npar,ibdiag,ned7)
         if (odebug(90)) print*,'output diagonal on blocks',ibdiag,'- ',
     ?                   ibdiag+lensec(npar)-1
         iblo7 = iblo7 + lensec(npar)
         call corlsr (izint)
         if (odebug(19)) then
            call outvec (q(idiaga),norb,'orbital diagonal:')
            call outvec (q(idiaga+norb),nci1,
     +                   'approximate configuration diagonal:')
            call outvec (q(idiagt),nci1,'true configuration diagonal:')
         end if
         call dcopy(nci1,q(idiagt),1,q(idiaga+norb),1)

c   blocks for trial vectors b, vectors Sb and Eb

         iblock = iblb
         do 307 j=1,maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(npar2)
            iblock = iblock + 1
307      continue

         iblock = ibls
         do 308 j=1,maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(npar2)
            iblock = iblock + 1
308      continue

         iblock = ible
         do 309 j=1,maxr
            iq(iblock) = iblo7
            iblo7 = iblo7 + lensec(npar2)
            iblock = iblock + 1
309      continue

c      iblock = iblbb
c      do 317 j=1,maxr
c      iq(iblock) = iblo7
c      iblo7 = iblo7 + lensec(npar2)
c317   iblock = iblock + 1
c
c      iblock = iblss
c      do 318 j=1,maxr
c      iq(iblock) = iblo7
c      iblo7 = iblo7 + lensec(npar2)
c318   iblock = iblock + 1
c
c      iblock = iblee
c      do 319 j=1,maxr
c      iq(iblock) = iblo7
c      iblo7 = iblo7 + lensec(npar2)
c319   iblock = iblock + 1

c   choose starting vectors

         if (nev(isym) .eq. 0) go to 3095
         call vectrs (q(idiaga),q(icivec),iq(iblb),iq(imspli),
     +                symbol,ivect,q,iq)
3095     call corlsr (ifi)

      if (ococo) then

c   first worry about space

         jcush = 10000
         ntotal = lrspce - jcush
         if (mclrsp .ge. ntotal) 
     +      call caserr ('not enough store for mclr')
         lenrs1 = (ntotal - mclrsp) / (2 * nci1)
         mxelen = 0
         do i=1,nirr
            np = 0
            do isym1=1,nirr
               isym2 = mult(i,isym1)
               np = np + mact(isym1) * mcore(isym2) + 
     +                   mvirt(isym1) * mact(isym2)
            enddo
            mxelen = max(mxelen,np)
         enddo
         iwarn = 0
         nblkmn = 20
286      mcplsp = ibegin + (intlen-1)/intrel + 1
     ?                   + 2 * ((nstraa-1)/intrel+1 
     +                   + (nstrbb-1)/intrel+1)
     ?                   + 3 * ((nblkmn-1)/intrel + 1)
     ?                   + 3 * ((maxrp3*nblkmn-1)/intrel + 1)
     ?                   + (maxpar+mxelen) * nblkmn * nblkmn
     ?                   + 2 * nfock
     ?                   + (norbti+norbai) * npair(isymr)
     ?                   + 2 * nbasq
         if (mcplsp .lt. ntotal) go to 586
         write (iwr,386) 
386      format ('warning: too much space needed for',
     +           ' construction of co part'/
     ?           'trying to decrease NBLKMN'//)
         iwarn = 1
         nblkmn = nblkmn - 1
         if (nblkmn .ge. 1) go to 286
         write (iwr,486)
486      format ('impossible to perform construction of co part'/
     ?           'increase core !')
         call caserr ('not enough store')
586      if (iwarn .eq. 1) write (iwr,686) nblkmn
686      format (' ... OK ! NBLKMN has been decreased to ',i2/)
         lenrs2 = (ntotal - mcplsp) / (2 * nci1)
         lenrs = min(lenrs1,lenrs2)

         inter = icori(intlen)
         iofa0 = icori(nstraa)
         iofb0 = icori(nstrbb)
         call setsto(intlen,0,iq(inter))
         call ciadrs (2,isss,iq(iofa0),iq(iofb0),iq(inter),na)
         if (isym .eq. isss) then
            iofa1 = iofa0
            iofb1 = iofb0
         else
            iofa1 = icori(nstraa)
            iofb1 = icori(nstrbb)
            call ciadrs (1,isym,iq(iofa1),iq(iofb1),iq(inter),na)
         end if
         nablk = (maxaa/2)*2 + 1
         nbblk = (maxbb/2)*2 + 1
         mmax   = icorrm() - 1
140      istore = 2*((nablk-1)/intrel+1) + (nbblk-1)/intrel+1
     ?          + 2*((maxrp3*nablk-1)/intrel+1) 
     +          + (maxrp3*nbblk-1)/intrel+1
     ?          + 2*lenrs*nci1
     ?          + (maxpar+mxelen)*nablk*nbblk + 2*(nfock+nbasq)
     ?          + (norbti+norbai)*npair(isymr)
         if (istore .le. mmax) go to 160
         if (nablk .gt. nbblk) go to 150
         nbblk = nbblk - 2
         go to 140
150      nablk = nablk - 2
         go to 140
160      nwa0 = icori(nablk)
         nwa1 = icori(nablk)
         nwb = icori(nbblk)
         iwa0 = icori(maxrp3*nablk)
         iwa1 = icori(maxrp3*nablk)
         iwb = icori(maxrp3*nbblk)
         ia = icorr(lenrs*nci1)
         ib = icorr(lenrs*nci1)
         id = icorr(maxpar*nablk*nbblk)
         ie = icorr(mxelen*nablk*nbblk)
         if = icorr(nfock)
         ig = icorr(nfock)
         iw = icorr((norbti+norbai)*npair(isymr))
         call couplx (q(ia),q(ib),q(icivec),q(id),q(ie),q(if),
     +                q(ig),q,q(iw),
     ?                iq(iofa0),iq(iofa1),iq(iofb0),iq(iofb1),iq(inter),
     ?                iq(iwa0),iq(iwa1),iq(iwb),iq(nwa0),iq(nwa1),
     +                iq(nwb), nablk,nbblk)
         call corlsi (inter)
      end if
      if (nev(isym) .eq. 0) go to 320

c   pointers for partitioning of q

      idenom = icorr(neig)
      ieigen = icorr(neig)
      ibleig = icori(neig)
c     placed here because of root catching (see below)
      ietare = icorr(maxr2) 
      ix     = icorr(npar2)
      isx    = icorr(npar2)
      ihx    = icorr(npar2)
      isred1 = icorr(maxrs)
      ihred1 = icorr(maxrs)
      ibred  = icorr(maxr)
      isred2 = icorr(maxr2s)
      ihred2 = icorr(maxr2s)
      ietaim = icorr(maxr2)
      isigma = icorr(maxr2)
      ixred  = icorr(maxr2s)
      ires   = icorr(neig)
      ireso  = icorr(neig)
      iresc  = icorr(neig)

      niter = 0
      call setsto(neig,0,iq(iconv))

      iblock = ibleig
      do 311 j=1,neig
      iq(iblock) = iblo7
      iblo7 = iblo7 + lensec(npar2)
311   iblock = iblock + 1

      tim(1) = cpulft(1)
      tim(2) = 0.0d0
      tim(3) = 0.0d0
      call vclr(shtim,1,10)
      call vclr(cotim,1,10)
      do 2671 i=1,neig
         ocvgd(i) = .false.
2671  continue
      ierr = 0
      rewind itable
      rewind inorms
      igam = idens1
      nstrt = neig
      write (iwr,131)
131   format (///89('*')/'*',10x,'MCLR iterations',62x,'*'/89('*')/)
1300  call lralgo (nstrt,q(icivec),q(igam),q(ix),q(isx),q(ihx),
     ?             q(isred1),q(ihred1),q(ibred),q(isred2),q(ihred2),
     ?             q(ietare),q(ietaim),q(isigma),q(ixred),
     ?             q(ires),q(ireso),q(iresc),q(idenom),q(ieigen),
     ?             iq(imspli),iq(iblb),iq(ibls),iq(ible),
     +             iq(ibleig),iq(iconv),q)
      if (ierr .ne. 0) go to 200
      nconv = isum (neig,iq(iconv),1)
      if (nconv.eq.neig .or. niter.gt.maxit) go to 200

c   restart iteration process with lowest eigenvectors of last iteration

      ix1 = ix
      ix2 = ix + npar
      iy1 = isx
      iy2 = isx + npar
      do 760 j=1,neig
      if (odebug(43)) print*,
     +               'trying to read preliminary e-vector from block',
     ?               iq(ibleig-1+j)
c      call rdedx (q(ix1),npar2,iq(ibleig-1+j),ned7)
      call rdedx (q(ix1),npar2,iq(ible-1+j),ned7)
      do 750 k=1,j-1
      call rdedx (q(iy1),npar2,iq(iblb-1+k),ned7)
      t = - ddot(npar,q(ix1),1,q(iy1),1) - 
     +      ddot(npar,q(ix2),1,q(iy2),1)
      call daxpy(npar,t,q(iy1),1,q(ix1),1)
      call daxpy(npar,t,q(iy2),1,q(ix2),1)
      t = - ddot(npar,q(ix1),1,q(iy2),1) - 
     +      ddot(npar,q(ix2),1,q(iy1),1)
      call daxpy(npar,t,q(iy2),1,q(ix1),1)
      call daxpy(npar,t,q(iy1),1,q(ix2),1)
750   continue
      t = dnrm2(npar2,q(ix1),1)
      if (1.0d0+t .eq. 1.0d0) 
     +      call caserr ('singularity in mclr restart')
      call symort (npar,q(ix1),q(ix2),sytol,oldep)
      if (oldep) call caserr ('singularity in restart symort call')
      if (odebug(90)) print*,
     +                'output restart trial vector',j,' on blocks',
     ?                iq(iblb-1+j),' -',iq(iblb-1+j)+lensec(npar2)-1
      call wrt3 (q(ix1),npar2,iq(iblb-1+j),ned7)
760   continue
      if (.not.onospl .and. nsplit.eq.0) then
         ocrazy = .false.
         do 770 j=neig,1,-1
         iq(imspli-1+2*j-1) = 1
         iq(imspli-1+2*j) = 2
         call rdedx (q(ix),npar2,iq(iblb-1+j),ned7)
         call dcopy(norb,q(ix),1,q(isx),1)
         call dcopy(norb,q(ix+npar),1,q(isx+norb),1)
         t = ddot(norb2,q(isx),1,q(isx),1)+1.0d0
         if (t.eq.1.0d0) then
            ocrazy = .true.
            iq(imspli-1+2*j-1) = 0
         else
            call wrt3 (q(isx),norb2,iq(iblb-1+2*j-1),ned7)
         endif
         call dcopy(nci1,q(ix+norb),1,q(isx),1)
         call dcopy(nci1,q(ix+npar+norb),1,q(isx+nci1),1)
         t = ddot(nci1,q(isx),1,q(isx),1)+1.0d0
         if (t.eq.1.0d0) then
            ocrazy = .true.
            iq(imspli-1+2*j) = 0
         else
            call wrt3 (q(isx),2*nci1,iq(iblb-1+2*j),ned7)
         endif
770      continue
         nstrt = 2 * neig
         if (ocrazy) then
c
c           this is insane, lets get rid of the zero vectors!!!
c
            i=0
            do j=1,nstrt
               if (iq(imspli-1+j).eq.1.and.j.gt.i+1) then
                  i=i+1
                  iq(imspli-1+i)=1
                  call rdedx (q(isx),norb2,iq(iblb-1+j),ned7)
                  call wrt3  (q(isx),norb2,iq(iblb-1+i),ned7)
               else if (iq(imspli-1+j).eq.2.and.j.gt.i+1) then
                  i=i+1
                  iq(imspli-1+i)=2
                  call rdedx (q(isx),2*nci1,iq(iblb-1+j),ned7)
                  call wrt3  (q(isx),2*nci1,iq(iblb-1+i),ned7)
               else
c                 skip zero vector
               endif
            enddo
            nstrt = i
         endif
      else
         call setsto(neig,3,iq(imspli))
      end if
      go to 1300

c   output of convergence behaviour

200   tim(4) = cpulft(1) - tim(1)
      niters = min(niter,maxit)
      write (iwr,312)
312   format ('iteration  reduced    eigenvalues'/
     +        'number     dimension'/ 89('-'))
      call lreign (itable,iwr,dummy,dummy,neig,niters,ntot,2)
      if (niter .gt. maxit) then
         write (iwr,1310)
1310     format (/
     +   'maximum number of iterations reached, no convergence.'/)
      else if (ierr .eq. 0) then
         write (iwr,99)
99       format (//89('*'))
         if (niter .eq. 1) write (iwr,100)
         if (niter .gt. 1) write (iwr,201) niter
100      format ('*',10x,
     +   'all eigenvalues have converged after one iteration.'
     +   ,26x,'*'/89('*')//)
201      format ('*',10x,'all eigenvalues have converged after ',i2, 
     ?           ' iterations.',26x,'*'/89('*')//)
      end if
c      write (iwr,313)
c313   format ('iteration             residues'/'number'/89('-'))
c      call outres (iwr,inorms,dummy,0,neig,niters,2)
      call residu (iwr,inorms,dummy,dummy,dummy,neig,niters,2)
      write (iwr,205) tim(2),tim(3),tim(4)-tim(2)-tim(3),tim(4)
205   format (/'timing analysis:'/
     +    'multiplication with trial vectors:      ',f14.4,' seconds'/
     +    'solution of reduced eigenvalue problems:',f14.4,' seconds'/
     +    'remaining tasks:                        ',f14.4,' seconds'/
     +    'total iterative procedure:              ',f14.4,' seconds')
         write(iwr,207)(shtim(i),i=1,7)
207   format (/10x,'time for the specific routines:'/
     ?        /10x,'trf ........',f14.4,' seconds'
     ?        /10x,'rpmult .....',f14.4,' seconds'
     ?        /10x,'modify .....',f14.4,' seconds'
     ?        /10x,'co .........',f14.4,' seconds'
     ?        /10x,'ocpart .....',f14.4,' seconds'
     ?        /10x,'cislow .....',f14.4,' seconds'
     ?        /10x,'invers .....',f14.4,' seconds')
      write (iwr,208) cotim(1),cotim(2),cotim(3),
     +                cotim(5),cotim(6),cotim(4)
208   format (20x,'analysis of co routine:'//
     +   20x,'generation of replacements ........ ',f14.4,' seconds'/
     +   20x,'formation of D matrix ............. ',f14.4,' seconds'/
     +   20x,'formation of matrices E1,E2 ....... ',f14.4,' seconds'/
     +   20x,'contribution to vco,wco ........... ',f14.4,' seconds'/
     +   20x,'formation of matrices Ex,Ey ....... ',f14.4,' seconds'/
     +   20x,'contributions to Dc,Pc ............ ',f14.4,' seconds'/)
      if (ierr .ne. 0) go to 470
      call corlsr (ix)

c   preparations for the computation of oscillator strengths

320   ix = icorr(npar2)
c      write (iwr,102) ix
c102   format (/'pointer is on ',i7,' before calculating oscillator strength'/)
      ife = icorr(neig)
c...next pointer for oscillator strengths of the very last eigenvectors
c...in case some roots got lost
      ife1 = icorr(neig)
      imyox = icorr(norb)
      imycx = icorr(nci1)
      imyoy = icorr(norb)
      imycy = icorr(nci1)
      imyoz = icorr(norb)
      imycz = icorr(nci1)
      nb1 = nbasis * nbasis
      lendip = nb1
      idipx = icorr(lendip)
      idipy = icorr(lendip)
      idipz = icorr(lendip)
      if (odebug(72)) print*,'myo addresses:',imyox,imyoy,imyoz
      if (odebug(72)) print*,'myc addresses:',imycx,imycy,imycz
      if (odebug(72)) print*,'dip addresses:',idipx,idipy,idipz
c      imyox1 = icorr(norb)
c      imyoy1 = icorr(norb)
c      imyoz1 = icorr(norb)
       isec134 = isect(134)
       m80 = 80
       call secget (isec134,m80,idiblo)
      idip = idipx
      do 400 icoord=1,3
      idiblo = idiblo + lensec(lendip)
      if (odebug(84)) print*,'dipole integrals, coordinate ',icoord
      call rdedx (q(idip),lendip,idiblo,ned3)
      if (odebug(84)) call outsqr(q(idip),nbasis,nbasis,nbasis,' ')
400   idip = idip + lendip
      call diporb (q(idipx),q(idipy),q(idipz),
     ?             q(imyox),q(imyoy),q(imyoz),q(idens))
c      if (odebug(72)) then
c  call outvec (q(imyox),norb,'myo for coordinate x:')
c  call outvec (q(imyoy),norb,'myo for coordinate y:')
c  call outvec (q(imyoz),norb,'myo for coordinate z:')
c      end if
      if (isym .ne. isss) go to 406
      if (odebug(72)) print*,'orthogonal correction:'
      idipx1 = idipx - 1
      idipy1 = idipy - 1
      idipz1 = idipz - 1
      dxcore = 0.0d0
      dycore = 0.0d0
      dzcore = 0.0d0
      do 405 i=1,ncore
      ii = (i-1)*nbasis + i
      dxcore = dxcore + q(idipx1+ii)
      dycore = dycore + q(idipy1+ii)
405   dzcore = dzcore + q(idipz1+ii)
      dxcore = dxcore + dxcore
      dycore = dycore + dycore
      dzcore = dzcore + dzcore
      if (odebug(72)) print*,'dxcore = ',dxcore
      if (odebug(72)) print*,'dycore = ',dycore
      if (odebug(72)) print*,'dzcore = ',dzcore
      if (odebug(72)) then
         dmx = dxcore
         dmy = dycore
         dmz = dzcore
         do it=1,nact
         do iu=1,nact
         dtu = q(idens-1+(it-1)*nact+iu)
         xtu = q(idipx-1+(it+ncore+nfreez-1)*nbasis+iu+ncore+nfreez)
         ytu = q(idipy-1+(it+ncore+nfreez-1)*nbasis+iu+ncore+nfreez)
         ztu = q(idipz-1+(it+ncore+nfreez-1)*nbasis+iu+ncore+nfreez)
         write(iwr,2037) it,iu,dtu,xtu,ytu,ztu
2037     format(2i3,4f14.9)
         dmx = dmx + dtu * xtu
         dmy = dmy + dtu * ytu
         dmz = dmz + dtu * ztu
         enddo
         enddo
         print*,'dmxyz with density matrices:'
         print*,dmx,dmy,dmz
      end if
406   idip = idipx - 1
      do 408 icoord=1,3
      do 407 it=nst+nfreez,nprim+nfreez
      itt = idip + (it-1)*nbasis + it
407   q(itt) = 0.5d0 * q(itt)
408   idip = idip + lendip
      inter  = icori(intlen)
      iofa   = icori(nstraa)
      iofb   = icori(nstrbb)
      call setsto(intlen,0,iq(inter))
      call ciadrs (2,isss,iq(iofa),iq(iofb),iq(inter),na)
      nablk = (maxaa/2)*2 + 1
      nbblk = (maxbb/2)*2 + 1
      mmax = icorrm()-1
410   istore = (nablk-1)/2 + 1 + (nbblk-1)/2 + 1
     ?         + (maxrp3*nablk-1)/2 + 1 + (maxrp3*nbblk-1)/2 + 1
     ?         + maxpar*nablk*nbblk
      if (istore.le.mmax) go to 430
      if (nablk.gt.nbblk) go to 420
      nbblk = nbblk - 2
      go to 410
420   nablk = nablk - 2
      go to 410
430   nwa = icori(nablk)
      nwb = icori(nbblk)
      iwa = icori(maxrp3*nablk)
      iwb = icori(maxrp3*nbblk)
      id = icorr(maxpar*nablk*nbblk)
      if (odebug(72)) print*,'id address: ',id
      call dipcon (q(icivec),q(id),q(idipx),q(idipy),q(idipz),
     +    q(imycx),q(imycy),q(imycz),iq(iofa),iq(iofb),
     +   iq(inter),iq(iwa),iq(iwb),iq(nwa),iq(nwb),nablk,nbblk)
      call corlsr (idipx)
      if (isym .eq. isss) then
         call daxpy(nci0,dxcore,q(icivec),1,q(imycx),1)
         call daxpy(nci0,dycore,q(icivec),1,q(imycy),1)
         call daxpy(nci0,dzcore,q(icivec),1,q(imycz),1)
      end if
      if (odebug(72)) then
         call outvec (q(imycx),nci1,'myc for coordinate x:')
         call outvec (q(imycy),nci1,'myc for coordinate y:')
         call outvec (q(imycz),nci1,'myc for coordinate z:')
      end if
c      if (niter .gt. maxit) go to 470
      if (nev(isym) .eq. 0) go to 900

c   compute oscillator strengths

      call vclr(q(ife),1,neig)
      ioo1 = ix
      icc1 = ix + norb
      ioo2 = ix + npar
      icc2 = icc1 + npar
      ib = ibleig
      ic = iconv
      id = idenom
      ie = ieigen
      if = ife
      if(odebug(73))call outvec (q(idenom),neig,'denominators <x,Sx>:')
      if(odebug(73))
     +            call outvec (denomy,neig,'denominators <y,Sigma y>:')
      do 460 j=1,neig
         if (iq(ic) .ne. 1) go to 450
         if (dabs(q(id)) .lt. dsmall) then
            write (iwr,3204) j
3204        format ('eigenvector ',i3,' cannot be scaled.')
            go to 450
         end if
         if (odebug(43).or.odebug(72)) write (iwr,438) iq(ib)
438      format ('lrmain: trying to read eigenvector from block ',i6)
         call rdedx (q(ix),npar2,iq(ib),ned7)
         call trnsfn (q(ix),2,1,npar,.true.,q)
         imyoq = imyox
         imycq = imycx
         do 440 icoord=1,3
            t =    ddot(norb,q(imyoq),1,q(ioo1),1) 
     ?           + ddot(nci1,q(imycq),1,q(icc1),1)
     ?           - ddot(norb,q(imyoq),1,q(ioo2),1) 
     ?           - ddot(nci1,q(imycq),1,q(icc2),1)
            q(if) = q(if) + t * t
            imyoq = imyoq + npar
            imycq = imycq + npar
440      continue
         q(if) = c23 * q(ie) * q(if) / q(id)
         fesum = fesum + q(if)
450      ib = ib + 1
         ic = ic + 1
         id = id + 1
         ie = ie + 1
         if = if + 1
460   continue

c   repeat the thing for the very last eigenvectors 
c   in case some roots got lost

      call vclr(q(ife1),1,neig)
      ioo1 = ix
      icc1 = ix + norb
      ioo2 = ix + npar
      icc2 = icc1 + npar
      ib = ible
      ie = ietare
      if = ife1
      print*,'check of root catching:'
      call outvec (dlost,neig,'lost denominators <x,Sx>:')
      call outive (iq(ib),neig,'block addresses')
      call outvec (q(ie),neig,'last eigenvalues:')
      do 560 j=1,neig
         if (dabs(dlost(j)) .lt. dsmall) then
            write (iwr,5204) j
5204        format ('very last eigenvector ',i3,' cannot be scaled.')
            go to 550
         end if
         if (odebug(43).or.odebug(72)) write (iwr,538) iq(ib)
538      format (
     +   'lrmain: trying to read very last eigenvector from block ',i6)
         call rdedx (q(ix),npar2,iq(ib),ned7)
         call trnsfn (q(ix),2,1,npar,.true.,q)
         imyoq = imyox
         imycq = imycx
         do 540 icoord=1,3
            t =    ddot(norb,q(imyoq),1,q(ioo1),1) 
     ?           + ddot(nci1,q(imycq),1,q(icc1),1)
     ?           - ddot(norb,q(imyoq),1,q(ioo2),1) 
     ?           - ddot(nci1,q(imycq),1,q(icc2),1)
            q(if) = q(if) + t * t
            imyoq = imyoq + npar
            imycq = imycq + npar
540      continue
         q(if) = c23 * q(ie) * q(if) / dlost(j)
550      ib = ib + 1
         ie = ie + 1
         if = if + 1
560   continue
      call outvec 
     ?   (q(ife1),neig,'oscillator strengths of very last eigenvectors')

c   output and analysis of results

c      ibas = icorr(0)
c      write (iwr,103) ibas
c103   format (/'pointer is on ',i7,' before lranal'/)
      call lranal (q(icivec),q(ix), q,        iq,
     ?             q(ieigen),q(ife),iq(ibleig),iq(iconv),
     ?             symbol,zymbol,2)

c   polarizability calculation

470   if (nfreq .eq. 0) go to 1000

      write (iwr,910)
910   format (//122('=')//)

900   continue
      isx    = icorr(npar2)
      ihx    = icorr(npar2)
      isred  = icorr(maxrs)
      ihred  = icorr(maxrs)
      ibred  = icorr(maxr)
      iared  = icorr(maxr2s)
      ifred1 = icorr(maxr)
      ifred2 = icorr(maxr2)
      iwork  = icori(maxr2)
      if (isym .eq. isss) then
         dmx = ddot(nci0,q(imycx),1,q(icivec),1)
         dmy = ddot(nci0,q(imycy),1,q(icivec),1)
         dmz = ddot(nci0,q(imycz),1,q(icivec),1)
         if (odebug(72)) print*,'dmxyz = ',dmx,dmy,dmz
         call daxpy(nci0,-dmx,q(icivec),1,q(imycx),1)
         call daxpy(nci0,-dmy,q(icivec),1,q(imycy),1)
         call daxpy(nci0,-dmz,q(icivec),1,q(imycz),1)
      end if
      thresh = eps
      igam = idens1
      imy1 = imyox
      do 850 icoor1=1,3
      if (dnrm2(npar,q(imy1),1) .le. tolmy) go to 850
      if (odebug(89)) then
         call outvec (q(imy1),npar,'vector my before trf:')
         if (isym .eq. isss) then
            temp = ddot(nci0,q(icivec),1,q(imy1+norb),1)
            print*,'overlap of my with c: ',temp
         end if
      end if
      call trnsfn (q(imy1),1,2,npar,.true.,q)
      do 840 i=1,nfreq

c   starting vector

      if(odebug(89))
     +       call outvec (q(imy1),npar,'vector my before lrupda:')
c      call vclr(q(isx),1,npar)
      call fmove (q(imy1),q(isx),npar)
      call lrupda (3,freq(i),q(imy1),q(isx),q(ihx),q(ihx+npar),
     +                       q(ix),q(ix+npar))
      if (freq(i) .eq. 0.0d0) then
c  if (isym .eq. isss) call orth (nci0,q(ix),q(icivec))
         t = dnrm2(npar,q(ix),1)
         if (1.0d0+t .eq. 1.0d0) call caserr ('vector f equals zero !')
         call dscal(npar,1.0d0/t,q(ix),1)
         call vclr(q(ix+npar),1,npar)
      else
c  if (isym .eq. isss) then
c     call orth (nci0,q(ix),q(icivec))
c     call orth (nci0,q(ix+npar),q(icivec))
c  end if
         call dscal(npar,-1.0d0,q(ix+npar),1)
         call symort (npar,q(ix),q(ix+npar),sytol,oldep)
         if (oldep) call caserr ('oldep = .true. !')
      end if
c      iq(iblb) = iblo7
      iq(imspli) = 3
      call wrt3 (q(ix),npar2,iq(iblb),ned7)
      if (odebug(90)) print*,
     +            'output initial trial vector for pola on blocks',
     ?                       iq(iblb),' -',iq(iblb)+lensec(npar2)-1
c      iblo7 = iblo7 + lensec(npar2)

c   solve linear response equation system

      write (iwr,823) isym,icoor1,freq(i)
823   format (/'Polarizability calculation for symmetry',i2,
     +  ', coordinate',i2,
     +  ', frequency',f9.4,' a.u.'/80('-')//6x,
     +  'iteration       orb residue   conf residue  total residue'/
     +  6x,'---------------------------------------------------------')
      niter = 0
c     icgvd = 0
      if (odebug(89)) 
     +       print*,'block addresses:',iq(iblb),iq(ibls),iq(ible)
824   call lrsolv (q(icivec),q(igam),q(ix),q(isx),q(ihx),freq(i),thresh,
     +        q(isred),q(ihred),q(ibred),q(iared),q(imy1),
     +        q(ifred1),q(ifred2),iq(iwork),iq(iblb),iq(ibls),iq(ible),
     +        iq(imspli),q)
      if (ierr .ne. 0) go to 840
      if (icvgd .eq. 1) go to 825
      if (niter .gt. maxit) go to 8245
      if (freq(i) .eq. 0.0d0) then
         t = dnrm2(npar,q(ix),1)
         if (1.0d0+t.eq.1.0d0) call caserr
     +                  ('restart trial vector equals zero !')
         call dscal(npar,1.0d0/t,q(ix),1)
         call vclr(q(ix+npar),1,npar)
      else
         call symort (npar,q(ix),q(ix+npar),sytol,oldep)
         if (oldep) call caserr 
     +                  ('symort problems in restart run of lrsolv')
      end if
      go to 824
8245  write (iwr,1310)
      go to 840
825   call trnsfn (q(ix),2,2,npar,.true.,q)
      imy2 = imyox
      do 830 icoor2=1,icoor1
      k = ipola(icoor1,icoor2)
      pola(k,i) = pola(k,i) - ddot(npar,q(ix     ),1,q(imy2),1)
     ?                      + ddot(npar,q(ix+npar),1,q(imy2),1)
830   imy2 = imy2 + npar
840   continue
850   imy1 = imy1 + npar

c   End loop over symmetries ISYM

1000  call corlsr (lrstar)

c   closing statements for TeX file

      call lranal (dummy,dummy,dummy,dummy,
     ?             dummy,dummy,dummy,dummy,
     ?             dummy,dummy,3)

      write (iwr,1010) fesum
1010  format (/'total oscillator strength:',5x,f9.4/)

      if (nfreq .gt. 0) call outpol (nfreq,freq,pola,ipola,gin)

      write (iwr,1015) cpulft(1)-time
1015  format (/
     +'total time for linear response calculation: ',f14.4,' seconds'//
     + 75('*')/'*',15x,'END OF LINEAR RESPONSE CALCULATION',24x,'*'/
     + 75('*'))
c
      call timana(30)
c
      return
      end
      subroutine lrmult (isplit,onlyoo,civec,d,x,s,e,q,iq)
c
c---------------------------------------------------------------------
c   Control routine for multiplication of S and E by a trial vector
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer iblbo, iblso, ibleo
      common /optima/ iblbo(50),iblso(50),ibleo(50)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
      dimension civec(*),d(*),x(*),s(*),e(*),q(*),iq(*)
c
      ibase = icorr(0)
      call vclr(e,1,npar2)
      call vclr(s,1,npar2)
      if (onlyoo) then
         nstrid = norb
         ix1 = 1
         ix2 = norb + 1
      else
         nstrid = npar
         ix1 = io1
         ix2 = io2
      end if
c
c   transform orbital parts of trial vector
c
      if (odebug(71) .or. isplit.eq.2 .or. norb.eq.0) go to 10
      time = cpulft(1)
      if (odebug(43) .and. onlyoo) then
         call outive (iblbo,maxro,'lrmult: b block addresses at start')
         call outive (iblso,maxro,'lrmult: s block addresses at start')
         call outive (ibleo,maxro,'lrmult: e block addresses at start')
      end if
      call trnsfn (x,2,1,nstrid,.true.,q)
      if (odebug(43) .and. onlyoo) then
      call outive (iblbo,maxro,'lrmult: b block addresses after transf')
      call outive (iblso,maxro,'lrmult: s block addresses after transf')
      call outive (ibleo,maxro,'lrmult: e block addresses after transf')
      end if
      shtim(1) = shtim(1) + cpulft(1) - time
10    continue
c
c   compute vector Sx
c
c   ... orbital part
c
      if (isplit .ne. 2) then
         call sorbit (x(ix1),x(ix2),s(ix1),s(ix2),d)
         if (odebug(43) .and. onlyoo) then
         call outive (iblbo,maxro,'lrmult: b block addresses after soo')
         call outive (iblso,maxro,'lrmult: s block addresses after soo')
         call outive (ibleo,maxro,'lrmult: e block addresses after soo')
         end if
      end if
      if (odebug(28)) call outvec (s,npar2,'sx after soo:')
c
c   ... configuration part
c
      if (isplit .ne. 1) then
         call dcopy(nci1,x(ic1),1,s(ic1),1)
         call vneg (x(ic2),1,s(ic2),1,nci1)
      end if
c
c   compute orbital-orbital parts v(oo),w(oo)
c
      is = icorr(0)
      mmax = icorrm() - 1
      if (lp*norb2 .gt. mmax) call caserr 
     +                       ('insufficient space for oo routine')
      time = cpulft(1)
      if (odebug(43)) then
         write (iwr,433) is,ibloab,lp,norb
433      format ('lrmult: is = ',i9,'  ibloab = ',i7,
     +           '   lp = ',i7,',  norb = ',i8)
      end if
      call rpmult (.true.,ibloab,lp,norb,ned7,x(ix1),x(ix2),
     +                                        e(ix1),e(ix2),q(is))
      if (odebug(43) .and. onlyoo) then
      call outive (iblbo,maxro,'lrmult: b block addresses after rpmult')
      call outive (iblso,maxro,'lrmult: s block addresses after rpmult')
      call outive (ibleo,maxro,'lrmult: e block addresses after rpmult')
      end if
      shtim(2) = shtim(2) + cpulft(1) - time
      if (odebug(26)) then
         call outvec (e(ix1),norb,'vo after oo:')
         call outvec (e(ix2),norb,'wo after oo:')
      end if
      call corlsr (ibase)
      if (onlyoo) return

c   compute coefficients of modified hamiltonian [H,X]

      if (ococo) then
         ia = icorr(lenrs*nci1) 
         ib = icorr(lenrs*nci1) 
         call coco (x(io1),x(ic1),x(io2),x(ic2),
     ?              e(io1),e(ic1),e(io2),e(ic2),q(ia),q(ib))
         call corlsr (ia)
         go to 300
      end if
      ifv  = icorr(n1int1)
      ifw  = icorr(n1int1)
      igv  = icorr(n1int2)
      igw  = icorr(n1int2)
      ifv1 = icorr(n1int1)
      ifw1 = icorr(n1int1)
      igv1 = icorr(n1int2)
      igw1 = icorr(n1int2)
      time = cpulft(1)
      call modify (x(io1),x(io2),q(ifv1),q(ifw1),q(ifv),q(ifw),
     ?                              q(igv1),q(igw1),q(igv),q(igw),q)
      shtim(3) = shtim(3)+cpulft(1)-time
      call corlsr (ifv1)

c   compute configuration-orbital part v(co),w(co) 
c   and transition density matrices Dc,Pc

      inter = icori(intlen)
      iofa0 = icori(nstraa)
      iofb0 = icori(nstrbb)
      call setsto(intlen,0,iq(inter))
      call ciadrs (2,isss,iq(iofa0),iq(iofb0),iq(inter),na)
      if (isym .eq. isss) then
         iofa1 = iofa0
         iofb1 = iofb0
      else
         iofa1 = icori(nstraa)
         iofb1 = icori(nstrbb)
         call ciadrs (1,isym,iq(iofa1),iq(iofb1),iq(inter),na)
      end if
      idc   = icorr(n1int1)
      ipc   = icorr(n1int2)
      iqc   = icorr(n1int2)
      irc   = icorr(n1int2)
      nablk = (maxaa/2)*2 + 1
      nbblk = (maxbb/2)*2 + 1
      mmax   = icorrm() - 1
140   istore = 2*((nablk-1)/2+1) + (nbblk-1)/2+1
     ?       + 2*((maxrp3*nablk-1)/2+1) + (maxrp3*nbblk-1)/2+1
     ?       + 3*maxpar*nablk*nbblk
      if (istore.le.mmax) go to 160
      if (nablk.gt.nbblk) go to 150
      nbblk = nbblk - 2
      go to 140
150   nablk = nablk - 2
      go to 140
160   nwa0   = icori(nablk)
      nwa1   = icori(nablk)
      nwb    = icori(nbblk)
      iwa0   = icori(maxrp3*nablk)
      iwa1   = icori(maxrp3*nablk)
      iwb    = icori(maxrp3*nbblk)
      id     = icorr(maxpar*nablk*nbblk)
      ie1    = icorr(maxpar*nablk*nbblk)
      ie2    = icorr(maxpar*nablk*nbblk)
      if (isplit .ge. 2) then
         call vclr(q(idc),1,n1int)
         call vclr(q(iqc),1,n1int2)
      end if
      time = cpulft(1)
      call copart (isplit,civec,q(id),q(ie1),q(ie2),
     ?             q(ifv),q(ifw),q(igv),q(igw),
     ?             x(ic1),x(ic2),e(ic1),e(ic2),
     ?             q(idc),q(ipc),q(iqc),q(irc),
     ?             iq(iofa0),iq(iofa1),iq(iofb0),iq(iofb1),iq(inter),
     ?             iq(iwa0),iq(iwa1),iq(iwb),
     ?             iq(nwa0),iq(nwa1),iq(nwb),nablk,nbblk)
      shtim(4) = shtim(4) + cpulft(1) - time
      if (isplit .eq. 1) go to 200
      call corlsi (nwa0)
      if (isym .eq. isss) then
         val = ddot(nci0,x(ic1),1,civec,1) - 
     +         ddot(nci0,x(ic2),1,civec,1)
         val = val + val
      end if
      if (odebug(32)) then
         call outvec (q(idc),n1int1,'dc after co:')
         call outvec (q(ipc),n1int2,'pc after co:')
      end if
      i1 = icorr(n1int1)
      i2 = icorr(n1int2)
      call dcopy(n1int1,q(idc),1,q(i1),1)
      call dcopy(n1int2,q(ipc),1,q(i2),1)
      time = cpulft(1)
      call ocpart (e(io1),e(io2),q(i1),q(i2),q)
      shtim(5) = shtim(5)+cpulft(1)-time
      if (odebug(26)) then
         call outvec (e(io1),norb,'vo after oc:')
         call outvec (e(io2),norb,'wo after oc:')
      end if
300   call corlsr (ibase)
      if (odebug(52)) then
         call outvec (e(ic1),nci1,'vc before cc part:')
         call outvec (e(ic2),nci1,'wc before cc part:')
      end if

c   compute configuration-configuration parts v(cc),w(cc)

      time = cpulft(1)
c...next line for common/multic/ so that CISLOW works
      nlen = npar
      izint = icorr(nint1+nint2)
      call rdedx (q(izint),nint1+nint2,ibzint,ned7)
      iofa = icori(nstraa)
      iofb = icori(nstrbb)
      inter = icori(intlen)
      call setsto(intlen,0,iq(inter))
      call ciadrs (2,isym,iq(iofa),iq(iofb),iq(inter),na)
      nablk = (maxaa/2)*2+1
      nbblk = (maxbb/2)*2+1
      mmax = icorrm()-1
340   istore = (maxrp3+1)*(nablk+nbblk) + 2*maxpar*nablk*nbblk
      if (istore .le. mmax) go to 360
      if (nablk .gt. nbblk) go to 350
      nbblk = nbblk - 2
      go to 340
350   nablk = nablk - 2
      go to 340
360   nwa = icori(nablk)
      nwb = icori(nbblk)
      iwa = icori(maxrp3*nablk)
      iwb = icori(maxrp3*nbblk)
      id = icorr(maxpar*nablk*nbblk)
      ie = icorr(maxpar*nablk*nbblk)
      call cislow (x(ic1),dummy,e(ic1),q(izint),dummy,
     +        2,2,0,0,ic3e,
     +        iq(iofa),iq(iofb),iq(inter),iq(iwa),iq(iwb),q(id),q(ie),
     +        nablk,nbblk,iq(nwa),iq(nwb))
      shtim(6) = shtim(6) + cpulft(1) - time
      if (odebug(98)) then
         print*,'hamiltonian on x before addition of core part:'
         call outvec (e(ic1),nci1,'1st vector:')
         call outvec (e(ic2),nci1,'2nd vector:')
      end if
      call daxpy(nci1,const,x(ic1),1,e(ic1),1)
      call daxpy(nci1,const,x(ic2),1,e(ic2),1)
      call corlsr (izint)
200   continue
      if (isym.eq.isss .and. odebug(35)) then
         print*,'orthogonality test:'
         temp = ddot(nci0,civec,1,e(ic1),1)
         print*,'<c,vc> = ',temp
      endif
      if (oshift) then
         call daxpy(npar2,cshift,s,1,e,1)
      end if
      call corlsr (ibase)
      return
      end
      subroutine lrsolv (civec,d,x,s,e,frqncy,thresh,sred,ered,bred,
     +                   ared,
     +                   f,fred1,fred2,iwork,iblb,ibls,ible,msplit,q)
c
c---------------------------------------------------------------------
c   Solves the multiconfigurational linear response equation system 
c
c                    (omega * S - E) c = f
c
c   using the iterative algorithm by Olsen & Joergensen.
c
c   On entry, x contains trial solution, on return, x contains converged 
c   solution.
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension civec(*),d(*),x(*),s(*),e(*),q(*)
      dimension sred(maxr,maxr),ered(maxr,maxr),bred(maxr)
      dimension ared(maxr2,maxr2),f(*),fred1(*),fred2(*)
      dimension iwork(maxr2),iblb(maxr),ibls(maxr),ible(maxr)
      dimension msplit(maxr)

      time = cpulft(1)
      ntot = 1
10    niter = niter + 1
      if (niter .gt. maxit) go to 1000

c   fill up the matrices sred,ered

      do 20 j=ntot,1,-1
      m = msplit(j)
      n = leng(m)
      ix1 = 1
      ix2 = n + 1
      ie1 = noff(m)
      ie2 = ie1 + npar
      if (j .lt. ntot) go to 14
      fred1(ntot) = ddot(n,x(ix1),1,f(ie1),1) - 
     +              ddot(n,x(ix2),1,f(ie1),1)
      go to (11,12,13),m
11    call bmove (x(ix2),x(io2),norb)
      call vclr(x(ic1),1,nci1)
      call vclr(x(ic2),1,nci1)
      go to 13
12    call bmove (x(ix2),x(ic2),nci1)
      call bmove (x(ix1),x(ic1),nci1)
      call vclr(x(io1),1,norb)
      call vclr(x(io2),1,norb)
13    call lrmult (m,.false.,civec,d,x,s,e,q,q)
      if (odebug(30)) call outvec (e,npar2,'e:')
      bred(ntot)      = ddot(n,x(ix1),1,e(ie2),1) + 
     +                  ddot(n,x(ix2),1,e(ie1),1)
      sred(ntot,ntot) = ddot(n,x(ix1),1,s(ie1),1) + 
     +                  ddot(n,x(ix2),1,s(ie2),1)
      ered(ntot,ntot) = ddot(n,x(ix1),1,e(ie1),1) + 
     +                  ddot(n,x(ix2),1,e(ie2),1)
      if (odebug(30)) print*,'h(i,i) = ',ered(ntot,ntot)
      if (m .ne. 2) call trnsfn (s,2,2,npar,.true.,q)
      call trnsfn (e,2,2,npar,.true.,q)
      if (ostore) then
         call wrt3 (s,npar2,ibls(ntot),ned7)
         call wrt3 (e,npar2,ible(ntot),ned7)
      end if
      go to 20
14    call rdedx (x,2*n,iblb(j),ned7)
      sred(j,ntot) = ddot(n,x(ix1),1,s(ie1),1) + 
     +               ddot(n,x(ix2),1,s(ie2),1)
      ered(j,ntot) = ddot(n,x(ix1),1,e(ie1),1) + 
     +               ddot(n,x(ix2),1,e(ie2),1)
      sred(ntot,j) = ddot(n,x(ix1),1,s(ie2),1) + 
     +               ddot(n,x(ix2),1,s(ie1),1)
      ered(ntot,j) = ddot(n,x(ix1),1,e(ie2),1) + 
     +               ddot(n,x(ix2),1,e(ie1),1)
20    continue

c   solve reduced equation system

c...for vector computers: no recurrence in this loop
      do 200 i = 1 , ntot
         fred2(i)      =  fred1(i)
         fred2(ntot+i) = -fred1(i)
200   continue
      call redeqs (maxr,ntot,frqncy,sred,ered,bred,ared,
     +                              fred2,iwork,ierr)
      if (odebug(91)) call outvec (fred2,2*ntot,'reduced solution:')
      if (ierr .ne. 0) go to 400

c   compute residue vector

      if (ostore) go to 110
      call blowup (ntot,x,e,fred2,msplit,iblb)
      call lrmult (3,.false.,civec,d,x,s,e,q,q)
      call trnsfn (s,2,2,npar,.true.,q)
      call trnsfn (e,2,2,npar,.true.,q)
      if (odebug(47)) then
         call outvec (s,npar2,'vector s*v')
         call outvec (e,npar2,'vector e*v')
      end if
      go to 150

c   store version

110   call vclr(s,1,npar2)
      call vclr(e,1,npar2)
      do 140 k=1,ntot
      call rdedx (x,npar2,ibls(k),ned7)
c...for vector computers: no recurrence in this loop
      do 120 i1=1,npar
      i2 = i1 + npar
      s(i1) = s(i1) + x(i1) * fred2(k) - x(i2) * fred2(ntot+k)
      s(i2) = s(i2) + x(i2) * fred2(k) - x(i1) * fred2(ntot+k)
120   continue
      call rdedx (x,npar2,ible(k),ned7)
c...for vector computers: no recurrence in this loop
      do 130 i1=1,npar
      i2 = i1 + npar
      e(i1) = e(i1) + x(i1) * fred2(k) + x(i2) * fred2(ntot+k)
      e(i2) = e(i2) + x(i2) * fred2(k) + x(i1) * fred2(ntot+k)
130   continue
140   continue

c   r = (omega * s - e) - f

150   continue
c...for vector computers: no recurrence in this loop
      do 50 i = 1 , npar
         e(     i) = frqncy * s(     i) - e(     i) - f(i)
         e(npar+i) = frqncy * s(npar+i) - e(npar+i) + f(i)
50    continue

c   residue vector now in e

      if (odebug(47)) then
         write (iwr,61) niter
61       format ('LRSOLV: residue vector no.',i3,':')
         call outvec (e,npar2,' ')
      end if

c   compute euclidean norm of residue vector

      ronorm = dnrsq(norb,e(io1),1) + dnrsq(norb,e(io2),1)
      rcnorm = dnrsq(nci1,e(ic1),1) + dnrsq(nci1,e(ic2),1)
      rnorm  = dsqrt(ronorm+rcnorm)
      ronorm = dsqrt(ronorm)
      rcnorm = dsqrt(rcnorm)
      write (iwr,60) niter,ronorm,rcnorm,rnorm
60    format (8x,i3,8x,3f14.9)

c   threshold ?

      if (rnorm .le. thresh) go to 230
      if (ntot.ge.maxr .or. ntot.ge.npar) go to 1000
      ntot = ntot + 1
      if (onospl .or. niter.lt.nsplit) then
         msplit(ntot) = 3
      else if (ronorm .ge. rcnorm) then
         msplit(ntot) = 1
      else
         msplit(ntot) = 2
      end if
      msave = msplit(ntot)
      msplit(ntot) = 3

c   update residue vector (updated vector in e)

      call lrupda (3,frqncy,e(io1),e(io2),x(io1),x(io2),e(io1),e(io2))
      if (odebug(48)) call outvec (e,npar2,'updated residue vector:')

c   gram-schmidt orthogonalisation

      if (frqncy .eq. 0.0d0) go to 100
      do 70 k = 1 , ntot-1
         if (msplit(k)*msplit(ntot) .eq. 2) go to 70
         call rdedx (x,2*leng(msplit(k)),iblb(k),ned7)
         m = min(msplit(k),msplit(ntot))
         n = leng(m)
         ix1 = 1
         if (msplit(ntot).eq.2 .and. msplit(k).eq.3) ix1 = ic1
         ix2 = ix1 + leng(msplit(k))
         ie1 = noff(m)
         ie2 = ie1 + npar
         t = - ddot(n,e(ie1),1,x(ix1),1) - 
     +         ddot(n,e(ie2),1,x(ix2),1)
         call daxpy(n,t,x(ix1),1,e(ie1),1)
         call daxpy(n,t,x(ix2),1,e(ie2),1)
         t = - ddot(n,e(ie1),1,x(ix2),1) - 
     +         ddot(n,e(ie2),1,x(ix1),1)
         call daxpy(n,t,x(ix2),1,e(ie1),1)
         call daxpy(n,t,x(ix1),1,e(ie2),1)
70    continue
      if (odebug(90)) 
     +    call outvec (e,npar2,'orthogonalized residue vector:')
      n = leng(msplit(ntot))
      n2 = 2 * n
      i1 = noff(msplit(ntot))
      i2 = i1 + npar
      t = dsqrt(dnrsq(n,e(i1),1)+dnrsq(n,e(i2),1))
      if (odebug(90)) print*,'norm: ',t
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,80) niter
80       format (
     +   'unable to finish iterative solution after iteration step ',i2/
     +   'since new trial vector is linear dependent.')
         ierr = 1
         go to 1000
      else ! perform symmetric orthonormalisation
         if (odebug(49)) then
         call outvec (e(i1),n,'first  part of vector to be "symorted":')
         call outvec (e(i2),n,'second part of vector to be "symorted":')
         end if
         call symort (n,e(i1),e(i2),sytol,oldep)
         if (oldep) then
            write (iwr,90) niter
90          format (
     +      'unable to finish iterative solution after iteration step ',
     +       i2/'because b and b# are (nearly) linearly dependent.')
            ierr = 1
            go to 1000
         end if
      end if
      go to 220
100   continue
      do 210 k = 1 , ntot-1
         if (msplit(k)*msplit(ntot) .eq. 2) go to 210
         call rdedx (x,2*leng(msplit(k)),iblb(k),ned7)
         m = min(msplit(k),msplit(ntot))
         n = leng(m)
         ix1 = 1
         if (msplit(ntot).eq.2 .and. msplit(k).eq.3) ix1 = ic1
         ie1 = noff(m)
         call orth (n,e(ie1),x(ix1))
210   continue
      if (odebug(90)) 
     +    call outvec (e,npar,'orthogonalized residue vector:')
      n = leng(msplit(ntot))
      n2 = 2 * n
      i1 = noff(msplit(ntot))
      i2 = i1 + npar
      t = dnrm2(n,e(i1),1)
      if (odebug(90)) print*,'norm: ',t
      if (1.0d0+t.eq.1.0d0) then
         write (iwr,80) niter
         ierr = 1
         go to 1000
      else
         call dscal(n,1.0d0/t,e(i1),1)
         call vclr(e(i2),1,n)
      end if

c   write new trial vector onto external file

220   msplit(ntot) = msave
      if (odebug(92) .and. msave.lt.3) then
         if (msave .eq. 1) then
            call vclr(e(ic1),1,nci1)
            call vclr(e(ic2),1,nci1)
         else
            call vclr(e(io1),1,norb)
            call vclr(e(io2),1,norb)
         end if
         msplit(ntot) = 3
      end if
      n = leng(msplit(ntot))
      n2 = 2 * n
      go to (121,122,123),msplit(ntot)
121   call dcopy(norb,e(io1),1,x(io1),1)
      call dcopy(norb,e(io2),1,x(ic1),1)
      go to 124
122   call dcopy(nci1,e(ic1),1,x(io1),1)
      call dcopy(nci1,e(ic2),1,x(io1+nci1),1)
      go to 124
123   call dcopy(npar2,e,1,x,1)
124   call wrt3 (x,n2,iblb(ntot),ned7)
      if (odebug(90)) print*,
     +    'LRSOLV: output trial vector',ntot,' on blocks',
     ?                iblb(ntot),iblb(ntot)+lensec(n2)-1
      go to 10

c   convergence reached

230   icvgd = 1
      if (ostore) then
         call blowup (ntot,x,e,fred2,msplit,iblb)
         call trnsfn (x,2,1,npar,.true.,q)
      end if
      if (niter .eq. 1) write (iwr,240) cpulft(1)-time
      if (niter .gt. 1) write (iwr,250) niter,cpulft(1)-time
240   format (/
     + 'convergence reached after one iteration in ',f9.4,' seconds.'/)
250   format (/'convergence reached after',i3,
     ?         ' iterations in',f9.4,' seconds.'/)
      go to 1000

c   error message

400   write (iwr,410) niter,ierr
410   format (/
     + 'iteration step ',i2,': condition number too small: ierr = ',i3/)
1000  continue
      return
      end
      subroutine lrupda (isplit,frequ,res1,res2,diag1,diag2,upda1,upda2)
c
c---------------------------------------------------------------------
c   Updates residual vector (res1,res2) by componentwise division
c   by the diagonal of  frequ * S - E
c
c   diag1,diag2 is space for diagonal, result will be in upda1,upda2
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension res1(*),res2(*),diag1(*),diag2(*),upda1(*),upda2(*)
      ind(i,j) = ((i-1)*(i-2))/2 + j
c
c   read in diagonal of matrix a
c
      if (odebug(43)) write (iwr,1) ibdiag
1     format ('lrupda: trying to read diagonal from block ',i6)
      call rdedx (diag1,npar,ibdiag,ned7)
      if (isplit .eq. 2) go to 50
c...for vector computers: no recurrence in this loop
      do 10 iit=1,nconac
      ii = (iit-1)/nact + 1
      it = iit - ilifa(ii)
      if (mult(itypea(it),itype(ii)) .eq. isymr) then
         k = ipo(ind(it+ncore,ii))
         sdiag = 2.0d0 - ddiag(it)
         u = diag1(k)
         diag1(k) =  (frequ - cshift) * sdiag - u
         diag2(k) = (-frequ - cshift) * sdiag - u
      end if
10    continue
c     t = frequ + frequ
      do 20 i = norbti+1 , norbti+norbai
         u = diag1(i)
         diag1(i) =  (frequ - cshift) * 2.0d0 - u
         diag2(i) = (-frequ - cshift) * 2.0d0 - u
20    continue
c...for vector computers: no recurrence in this loop
      do 30 iat=1,nvinac
      ia = (iat-1)/nact + 1
      it = iat - ilifa(ia)
      if (mult(itype(ia+nprim),itypea(it)) .eq. isymr) then
         k = ipo(ind(ia+nprim,it+ncore))
         sdiag = ddiag(it) 
         u = diag1(k)
         diag1(k) =  (frequ - cshift) * sdiag - u
         diag2(k) = (-frequ - cshift) * sdiag - u
      end if
30    continue
      if (isplit .eq. 1) go to 70
50    continue
      do 60 i = norb+1 , npar
         u = diag1(i)
         diag1(i) = ( frequ - cshift) - u
         diag2(i) = (-frequ - cshift) - u
60    continue
70    n = noff(isplit)
      if (odebug(89)) call outvec (diag1,npar,'diag1 in lrupda:')
      if (odebug(89)) print*,'first  call of VDIVZ'
      call vdivz (res1(n),1,diag1(n),1,1.0d0,upda1(n),1,leng(isplit))
      if (odebug(89)) print*,'second call of VDIVZ'
      call vdivz (res2(n),1,diag2(n),1,1.0d0,upda2(n),1,leng(isplit))

      return
      end
      subroutine lrzint (q)
c
c---------------------------------------------------------------------
c   Collects the active-active one-electron integrals
c
c       h(tu) = FI(tu) - 0.5 * sum_x (tx|xu)
c
c   and the all-active two-electron integrals
c
c       0.5 * (tu|vx)
c
c   and puts them to scratch file ED7, beginning at block IBZINT
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      dimension q(*)
      dimension mk1(8),mk2(8)
      ibase = icorr(0)
      time = cpulft(1)
      excmax = 0.0d0
      ifi = icorr(nfock)
      call rdedx (q(ifi),nfock,iblofi,ned7)
      ifock = ifi - 1
      do 10 i=1,nirr
      mk2(i) = ifock
10    ifock = ifock + nsymm(i)**2
      izint = icorr(nint1+nint2)
      izint1 = izint - 1
      i1 = icorr(nact2)
      ij = 0
      ij2 = 0
      intpos = 0
      intfil = ned0
      iblf = iblkaa
      intmod = 0
      do 90 i=1,nact
      do 90 j=1,i
      j1 = i1
      isymij = mult(itypea(i),itypea(j))
      ij = ij + 1
      ijkl = ic2e(ij)
      do 20 isyma=1,nirr
      isymb = mult(isyma,isymij)
      if (isymb .gt. isyma) go to 20
      n = mact(isyma) * mact(isymb)
      if (n .eq. 0) go to 20
      call intin (q(j1),n)
      mk1(isyma) = j1-1
      j1 = j1 + n
20    continue
      do 31 k=nst,nprim
      isymk = itype(k)
      do 30 l=nst,k
      isyml = itype(l)
      if (mult(isymk,isyml) .ne. isymij) go to 30
      ijkl = ijkl + 1
      if (isymk .ge. isyml) then
         ip = mk1(isymk) + (npoint(l)-1)*mact(isymk) + npoint(k)
      else
         ip = mk1(isyml) + (npoint(k)-1)*mact(isyml) + npoint(l)
      end if
      q(izint1+ijkl) = 0.5d0 * q(ip)
30    continue
31    continue
      if (isymij .ne. 1) go to 90
      ij2 = ij2 + 1
      isymi = itypea(i)
      mci = mcore(isymi)
      g = q(mk2(isymi)+
     +    (npoint(j+ncore)+mci-1)*nsymm(isymi)+npoint(i+ncore)+mci)
      if (odebug(76)) print*,'fock matrix element',ij2,': ',g
      do 80 isyma=1,nirr
      if (mact(isyma) .eq. 0) go to 80
      call intin (q(i1),mact(isyma)**2)
      if (odebug(76)) then
         print*,'exchange integrals, symmetry ',isyma
         call outsqr (q(i1),mact(isyma),mact(isyma),mact(isyma),' ')
      end if
      do 40 k=1,mact(isyma)
40    g = g - 0.5d0 * q(i1-1+(k-1)*mact(isyma)+k)
      if (i .ne. j) go to 80
      if (isyma .ne. isymi) go to 60
      do 50 k=1,mact(isyma)
      if (k .ne. npoint(i+ncore)) then
         if (odebug(7)) 
     +   write (iwr,51) excmax,q(i1-1+(k-1)*mact(isyma)+k)
51       format ('old excmax: ',f9.4,
     +           '  current exchange integral: ',f9.4)
         excmax = dmax1(excmax,q(i1-1+(k-1)*mact(isyma)+k))
      end if
50    continue
      go to 80
60    continue
      do 70 k=1,mact(isyma)
      if (odebug(7)) write (iwr,51) excmax,q(i1-1+(k-1)*mact(isyma)+k)
70    excmax = dmax1(excmax,q(i1-1+(k-1)*mact(isyma)+k))
80    continue
      q(izint1+nint2+ij2) = g
90    continue
      ibzint = iblo7
      call wrt3 (q(izint),nint1+nint2,ibzint,ned7)
      if (odebug(90)) print*,'output zint to blocks ',iblo7,
     ?                       ' -',iblo7+lensec(nint1+nint2)-1
      iblo7 = iblo7 + lensec(nint1+nint2)
      if (odebug(76)) call outvec (q(izint),nint1+nint2,'zint:')
      write (iwr,100) cpulft(1)-time
100   format ('making zint ',42('.'),f9.4,' seconds')
      call corlsr (ibase)
      return
      end
      subroutine modify (yo,zo,fv1,fw1,fv2,fw2,gv1,gw1,gv2,gw2,q)
c
c---------------------------------------------------------------------
c   Computes the coefficients of the modified Hamiltonians
c   [H_0,X] and [H_0,X+]. The final coefficients are contained
c   in arrays FV2, FW2 (one-electron) and GV2, GW2 (two-electron).
c   FV1, FW1, GV1, GW1 are auxiliary arrays.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
c
      dimension yo(*),zo(*),q(*)
      dimension fv1(*),fw1(*),fv2(*),fw2(*),gv1(*),gw1(*),gv2(*),gw2(*)
      dimension mk1(8),mk2(8),ifvw(8)
      ibase = icorr(0)
c
c   read in inactive fock matrix
c
      ifock = icorr(nfock)
      call rdedx (q(ifock),nfock,iblofi,ned7)
      if (odebug(44)) then
         call outvec (yo,norb,'yo in MODIFY')
         call outvec (zo,norb,'zo in MODIFY')
      end if
      call vclr(fv1,1,n1int1)
      call vclr(fw1,1,n1int1)

c   construct fock matrix contributions to fv,fw

      ifi = ifock
      ifv = 1
      do 10 isyma = 1 , nirr
         isymb = mult(isyma,isymr)
         iit = ifi + nsymm(isyma) * mcore(isyma)
         iat = iit + nprm(isyma)
c---------------------------------------------------------------------
c     fv(tu)  <--  sum_i  zo(ui) * fi(it)                            |
c---------------------------------------------------------------------
         call mxmaa (zo(iorbti(isymb)+1),mcore(isyma),1,
     +               q(iit),1,nsymm(isyma),fv1(ifv),1,mact(isymb),
     +               mact(isymb),mcore(isyma),mact(isyma))
c---------------------------------------------------------------------
c     fw(tu)  <--  sum_i  yo(ui) * fi(it)                            |
c---------------------------------------------------------------------
         call mxmaa (yo(iorbti(isymb)+1),mcore(isyma),1,
     +               q(iit),1,nsymm(isyma),fw1(ifv),1,mact(isymb),
     +               mact(isymb),mcore(isyma),mact(isyma))
c----------------------------------------------------------------------
c     fv(tu)  <--  sum_a  yo(au) * fi(at)                            |
c----------------------------------------------------------------------
         call mxmb  (yo(iorbat(isyma)+1),1,mact(isymb),
     +               q(iat),1,nsymm(isyma),fv1(ifv),1,mact(isymb),
     +               mact(isymb),mvirt(isyma),mact(isyma))
c----------------------------------------------------------------------
c     fw(tu)  <--  sum_a  zo(au) * fi(at)                            |
c----------------------------------------------------------------------
         call mxmb  (zo(iorbat(isyma)+1),1,mact(isymb),
     +               q(iat),1,nsymm(isyma),fw1(ifv),1,mact(isymb),
     +               mact(isymb),mvirt(isyma),mact(isyma))
         ifi = ifi + nsymm(isyma)**2
         ifv = ifv + mact(isyma) * mact(isymb)
10    continue
      if (odebug(44).and.nact2.le.511) then
         call vclr(gin,1,nact2)
         do i=1,nact
         do j=1,nact
         if(mult(itypea(i),itypea(j)).eq.isymr) 
     ?           gin(ilifa(i)+j)=fv1(ipos(ilifa(i)+j))
         enddo
         enddo
         write (iwr,281)
281      format(/'array fv before integral contributions:'/)
         do 282 i=1,nact
282      write (iwr,183)(gin(ilifa(i)+j),j=1,nact)
      end if

c   correction term for c-o part

      if (isym .ne. isss) go to 70
      fvisum = 0.0d0
      if1 = ifock - 1
      do 60 isyma = 1 , nirr
         if2 = if1 + mcore(isyma) * nsymm(isyma)
         iyz = iorbti(isyma)
         do 30 it = 1 , mact(isyma)
            do 20 ii=1,mcore(isyma)
               fvisum = fvisum + q(if2+ii) * (yo(iyz+ii)-zo(iyz+ii))
20          continue
            if2 = if2 + nsymm(isyma)
            iyz = iyz + mcore(isyma)
30       continue
         if2 = if1 + nprm(isyma) * nsymm(isyma)
         iyz = iorbai(isyma)
         do 50 ia = 1 , mvirt(isyma)
            do 40 ii = 1 , mcore(isyma)
               fvisum = fvisum + q(if2+ii) * (yo(iyz+ii)-zo(iyz+ii))
40          continue
            if2 = if2 + nsymm(isyma)
            iyz = iyz + mcore(isyma)
50       continue
         if1 = if1 + nsymm(isyma)**2
60    continue
      fvisum = fvisum + fvisum
70    call corlsr (ifock)

c   read in W matrix and construct W(T) x yo, W(T) x zo contributions 

      irow = norbti + norbai
      icol = npair(isymr)
      isize = irow * icol
      if (isize .eq. 0) go to 80
      iw = icorr(isize)
      if (odebug(43)) print*,'trying to read W matrix from block',iblkw,
     ?                       ' on ED7'
      call rdedx (q(iw),isize,iblkw,ned7)
      call mxmb (q(iw),irow,1, yo,1,irow, fv1,1,icol, icol,irow,1)
      call mxmb (q(iw),irow,1, zo,1,irow, fw1,1,icol, icol,irow,1)
      call corlsr (iw)

c   construct gv1,gw1

80    call vclr(gv1,1,n1int2)
      call vclr(gw1,1,n1int2)
      intpos = 0
      intfil = ned0
      iblf = iblkgv
      intmod = 0
      i1 = icorr(nbasq)
      do 2000 iv = 1 , nact
         do 1000 ix = 1 , iv
            j1 = i1
            isymvx = mult(itypea(iv),itypea(ix))
            isymtu = mult(isymvx,isymr)
            do 200 isyma = 1 , nirr
               isymb = mult(isyma,isymvx)
               if (isymb .gt. isyma) go to 200
               mca = mcore(isyma)
               maa = mact(isyma)
               mva = mvirt(isyma)
               mcb = mcore(isymb)
               mab = mact(isymb)
               mvb = mvirt(isymb)
               m = maa * (mcb+mvb) + mab * (mca+mva)
               if (m .eq. 0) go to 200
               call intin (q(j1),m)
               mk1(isyma) = j1
               mk2(isyma) = j1 + maa*mcb + mab*mca
               j1 = j1 + m
               if (isymb .eq. isyma) go to 200
               mk1(isymb) = mk1(isyma) + maa * mcb
               mk2(isymb) = mk2(isyma) + mab * mva
200         continue
            ioffg0 = ioff1(ilifa(iv)+ix) + 1
            ioffg = ioffg0
            do 400 isyma = 1 , nirr
               isymb = mult(isyma,isymvx)
               isymu = mult(isymb,isymr)
               mca = mcore(isyma)
               maa = mact(isyma)
               mva = mvirt(isyma)
               mcb = mcore(isymb)
               mab = mact(isymb)
               mvb = mvirt(isymb)
               mau = mact(isymu)
c
c     gv(tuvx) <--- sum_i (ti|vx) * zo(ui)
c     gw(tuvx) <--- sum_i (ti|vx) * yo(ui)
c
c         t            i            t
c       -----        -----        -----
c     u | g |   =  u | x |  x   i | I |
c       -----        -----        -----
c     (mau,maa)    (mau,mcb)    (mcb,maa)
c
c                               isyma .ge. isymb : 't'
c                               isyma .lt. isymb : 'n'
c
              if (isyma .ge. isymb) then
               call mxmb (zo(iorbti(isymu)+1),mcb,1,q(mk1(isyma)),maa,1,
     ?              gv1(ioffg),1,mau,mau,mcb,maa)
               call mxmb (yo(iorbti(isymu)+1),mcb,1,q(mk1(isyma)),maa,1,
     ?              gw1(ioffg),1,mau,mau,mcb,maa)
              else
               call mxmb (zo(iorbti(isymu)+1),mcb,1,q(mk1(isyma)),1,mcb,
     ?              gv1(ioffg),1,mau,mau,mcb,maa)
               call mxmb (yo(iorbti(isymu)+1),mcb,1,q(mk1(isyma)),1,mcb,
     ?              gw1(ioffg),1,mau,mau,mcb,maa)
               end if
c
c     gv(tuvx) <--- sum_i (at|vx) * yo(au)
c     gw(tuvx) <--- sum_i (at|vx) * zo(au)
c
c         t            a            t
c       -----        -----        -----
c     u | g |   =  u | x |  x   a | I |
c       -----        -----        -----
c     (mau,maa)    (mau,mvb)    (mvb,maa)
c
c                               isymb .ge. isyma : 'n'
c                               isymb .lt. isyma : 't'
c
              if (isymb .ge. isyma) then
               call mxmb (yo(iorbat(isymb)+1),1,mau,q(mk2(isymb)),1,mvb,
     ?                    gv1(ioffg),1,mau,mau,mvb,maa)
               call mxmb (zo(iorbat(isymb)+1),1,mau,q(mk2(isymb)),1,mvb,
     ?                    gw1(ioffg),1,mau,mau,mvb,maa)
              else
               call mxmb (yo(iorbat(isymb)+1),1,mau,q(mk2(isymb)),maa,1,
     ?                    gv1(ioffg),1,mau,mau,mvb,maa)
               call mxmb (zo(iorbat(isymb)+1),1,mau,q(mk2(isymb)),maa,1,
     ?                    gw1(ioffg),1,mau,mau,mvb,maa)
               end if
               ioffg = ioffg + maa * mau
400         continue

c   copy (vx) to (xv)

            if (ix .eq. iv) go to 1000
            call fmove (gv1(ioffg0),gv1(ioff1(ilifa(ix)+iv)+1),
     +                  npair(isymtu))
            call fmove (gw1(ioffg0),gw1(ioff1(ilifa(ix)+iv)+1),
     +                  npair(isymtu))
1000     continue
2000  continue

c   move one-electron contribution of gv to fv

      do 180 ix = 1 , nact
c...for vector computers: no recurrence in this loop
         do 170 itu = 1 , nact2
            it = (itu-1)/nact + 1
            iu = itu - ilifa(it)
            if (mult(itypea(it),itypea(iu)) .eq. isymr) then
               fv1(ipos(itu)) = fv1(ipos(itu)) - 
     ?                  gv1(ipos(ilifa(it)+ix)+ioff1(ilifa(ix)+iu))
               fw1(ipos(itu)) = fw1(ipos(itu)) - 
     ?                  gw1(ipos(ilifa(it)+ix)+ioff1(ilifa(ix)+iu))
            end if
170      continue
180   continue

c   construct final fv,fw

      k = 1
      do 510 isyma = 1 , nirr
         ifvw(isyma) = k
         isymb = mult(isyma,isymr)
         k = k + mact(isyma) * mact(isymb)
510   continue
      do 520 isyma = 1 , nirr
         isymb = mult(isyma,isymr)
         m = mact(isyma) * mact(isymb)
         if (m .eq. 0) go to 520
         isymc = ifvw(isymb)
         call trnsps (fw1(ifvw(isyma)),fv2(ifvw(isymb)),mact(isymb),
     +                mact(isyma))
         call daxpy(m
     +             ,-1.0d0,fv1(isymc),1,fv2(isymc),1)
         call trnsps (fv2(ifvw(isymb)),fw2(ifvw(isyma)),mact(isyma),
     +                mact(isymb))
520   continue
      call dscal(npair(isymr),-1.0d0,fv2,1)
      if (odebug(44).and.nact2.le.511) then
         call vclr(gin,1,nact2)
         do i=1,nact
         do j=1,nact
         if(mult(itypea(i),itypea(j)).eq.isymr) 
     ?      gin(ilifa(i)+j)=fv2(ipos(ilifa(i)+j))
         enddo
         enddo
         write (iwr,181)
181      format (/'array fv:'/)
         do 182 i=1,nact
182      write (iwr,183)(gin(ilifa(i)+j),j=1,nact)
183      format (10f9.4)
      end if

c   construct final gv,gw

      call invers (gw1,gv2)
      call daxpy(n1int2,-1.0d0,gv1,1,gv2,1)
      call fmove (gv2,gv1,n1int2)
      call invers (gv1,gw2)
      call dscal(n1int2,-1.0d0,gv2,1)
      if (odebug(45)) call outvec (gv2,n1int2,'array gv from modify:')
      call corlsr (ibase)

      return
      end
      function numcsf (nact,na,nb,isymm,itypea,mult)

c---------------------------------------------------------------------
c   Computes the number of spin-adapted configurations.
c   Modified version of the GAMESS routine ncsf. To my opinion,
c   there was a mistake with open shells. Anyway, the function
c   value is never used, it's just for information.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      integer icga(31),itypea(nact),mult(8,8),ib(33)
      numb = 0
      do 60 nopsh=na-nb,na+nb,2
      ndoc = (na+nb-nopsh)/2
      nuoc = nact-nopsh-ndoc
      if (nuoc.lt.0) goto 70
      if (nopsh.eq.0) then
         if (isymm.eq.1) nsing = 1
         if (isymm.ne.1) nsing = 0
      else
         nsing = 0
         mt = 0
10       call string (mt,nact,nopsh,icga,itypea,jmta,iua,mult)
         if (mt.eq.0) goto 20
         if (jmta.eq.isymm) nsing = nsing + 1
         goto 10
20       continue
      endif
c...scan branching diagram for number of spin eigenfunctions
      do 30 j=2,nopsh+2
30    ib(j) = 0
      ib(1) = 1
      do 50 i=1,nopsh
      do 40 j=i+1,2,-2
40    ib(j) = ib(j-1) + ib(j+1)
50    ib(1) = ib(2)
      nbrnch = ib(na-nb+1)
      numb = numb + nsing * nbrnch * ibinom(nact-nopsh,ndoc)
60    continue
70    numcsf = numb
      return
      end
      subroutine occdet (icc,icga,icgb,nalpha,nbeta,nstrga,nstrgb,
     +                   isymex)

c---------------------------------------------------------------------
c   Generates the occupation pattern for determinant icc.
c   Identical to GAMESS routine detocc except that a modified version
c   of the string generator (namely, "strong") is called.
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension icga(nalpha),icgb(nbeta),nstrga(8),nstrgb(8)
      ii = 0
      do 10 ibs=1,nirr
      ias = mult(ibs,isymex)
      jj = ii + nstrga(ias) * nstrgb(ibs)
      if (ii.lt.icc .and. jj.ge.icc)  go to 20
10    ii = jj
20    ibs = mult(ias,isymex)
      ia = (icc-ii-1) / nstrgb(ibs) + 1
      ib = icc - ii - (ia-1) * nstrgb(ibs)
      mta = 0
      do 40 ii=1,ia
30    call strong (mta,nact,nalpha,icga,itypea,jmta,iua,mult)
      if (jmta .ne. ias) go to 30
40    continue
      mtb = 0
      do 60 ii=1,ib
50    call strong (mtb,nact,nbeta,icgb,itypea,jmtb,iub,mult)
      if (jmtb .ne. ibs) go to 50
60    continue
      return
      end
      subroutine ocpart (vo,wo,dc,pc,q)
c
c---------------------------------------------------------------------
c   Constructs the o-c part of the vector E*x
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
c
      dimension vo(1),wo(1),dc(1),pc(1),q(1)
      dimension mk1(8),mk2(8),ioff(8)
      ibase = icorr(0)
c
c   read in inactive fock matrix
c
      ifi = icorr(nfock)
      call rdedx (q(ifi),nfock,iblofi,ned7)
      k = 0
      do 10 isyma=1,nirr
      ioff(isyma) = k
10    k = k + nsymm(isyma)**2
c
c   non-orthogonality contribution
c
      if (isym .ne. isss) go to 19
      do 17 isyma=1,nirr
      if = ifi + ioff(isyma) + mcore(isyma) * nsymm(isyma)
      iv = iorbti(isyma) + 1
      do 12 it=1,mact(isyma)
      call daxpy(mcore(isyma), val,q(if),1,vo(iv),1)
      call daxpy(mcore(isyma),-val,q(if),1,wo(iv),1)
      if = if + nsymm(isyma)
12    iv = iv + mcore(isyma)
      iv = iorbai(isyma) + 1
      do 13 ia=1,mvirt(isyma)
      call daxpy(mcore(isyma), val,q(if),1,vo(iv),1)
      call daxpy(mcore(isyma),-val,q(if),1,wo(iv),1)
      if = if + nsymm(isyma)
13    iv = iv + mcore(isyma)
17    continue
c
c   FI contributions to vo,wo
c
19    do 20 isyma=1,nirr
      isymb = mult(isyma,isymr)
      ifiix = ifi + ioff(isymb) + mcore(isymb) * nsymm(isymb)
      ifixa = ifi + ioff(isyma) + nprm(isyma) * nsymm(isyma) 
     +            + mcore(isyma)
c---------------------------------------------------------------------
c     vo(ti)  <--  - sum_x  fi(ix) * dc(tx)                           |
c---------------------------------------------------------------------
      if(mcore(isymb).ne.0 .and. mact(isymb).ne.0
     +   .and. mact(isyma).ne.0) then
        call mxmbn (q(ifiix),1,nsymm(isymb),
     +              dc(iofsd(isymr,isyma)+1),1,mact(isymb),
     +              vo(iorbti(isyma)+1),1,mcore(isymb),
     +              mcore(isymb),mact(isymb),mact(isyma))
      endif
c---------------------------------------------------------------------
c     vo(at)  <--    sum_x  dc(xt) * fi(xa)                           |
c----------------------------------------------------------------------
      if (odebug(87)) then
         print*,'dc for symmetry ',isyma
         call outsqr (dc(iofsd(isymr,isyma)+1),
     ?                mact(isymb),mact(isymb),mact(isyma),' ')
         print*,'fi(x,a) for symmetry ',isyma
         call outsqr (q(ifixa),nsymm(isyma),mact(isyma),
     +                mvirt(isyma),' ')
      end if
      call mxmb (dc(iofsd(isymr,isyma)+1),1,mact(isymb),
     ?           q(ifixa),1,nsymm(isyma),
     ?           vo(iorbat(isyma)+1),1,mact(isymb),
     ?           mact(isymb),mact(isyma),mvirt(isyma))
c----------------------------------------------------------------------
c     wo(ti)  <--    sum_x  fi(ix) * dc(xt)                           |
c----------------------------------------------------------------------
      call mxmb (q(ifiix),1,nsymm(isymb),
     ?           dc(iofsd(isymr,isymb)+1),mact(isyma),1,
     ?           wo(iorbti(isyma)+1),1,mcore(isymb),
     ?           mcore(isymb),mact(isymb),mact(isyma))
c----------------------------------------------------------------------
c     wo(at)  <--  - sum_x  dc(tx) * fi(xa)                           |
c----------------------------------------------------------------------
      if (mvirt(isyma).ne.0 .and. mact(isymb).ne.0
     +   .and. mact(isyma).ne.0) then
        call mxmbn (dc(iofsd(isymr,isymb)+1),mact(isyma),1,
     +              q(ifixa),1,nsymm(isyma),
     +              wo(iorbat(isyma)+1),1,mact(isymb),
     +              mact(isymb),mact(isyma),mvirt(isyma))
      endif
20    continue
      if (odebug(58)) call outvec (vo,norb,'vo in ocpart after fi sums')
      if (odebug(58)) call outvec (wo,norb,'wo in ocpart after fi sums')
      call corlsr (ifi)

c   W matrix contributions

      irow = norbti + norbai
      icol = npair(isymr)
      isize = irow * icol
      if (isize .eq. 0) go to 40

c   construct Dc(T)

      idct = icorr(npair(isymr))
      do 30 isyma=1,nirr
      isymb = mult(isyma,isymr)
30    call trnsps (dc(iofsd(isymr,isyma)+1),q(idct+iofsd(isymr,isymb)),
     ?             mact(isymb),mact(isyma))

c   read in W matrix and construct contributions W x Dc, W x Dc(T) of vo,wo

      iw = icorr(isize)
      if (odebug(43)) print*,'trying to read W matrix from block',iblkw,
     ?                       ' on ED7'
      call rdedx (q(iw),isize,iblkw,ned7)
      if (irow.ne.0.and.icol.ne.0) then
      call mxmb  (q(iw),1,irow, dc     ,1,icol, vo,1,irow, irow,icol,1)
      call mxmbn (q(iw),1,irow, q(idct),1,icol, wo,1,irow, irow,icol,1)
      if (odebug(58)) call outvec(vo,norb,'vo in ocpart after W matrix')
      if (odebug(58)) call outvec(wo,norb,'wo in ocpart after W matrix')
      endif
      call corlsr (idct)

c   read in integrals (ui|vx),(au|vx) for the I x Pc contributions of vo,wo

40    intpos = 0
      intfil = ned0
      iblf = iblkgv
      intmod = 0
      i1 = icorr(nbasq)
      do 1000 iv=1,nact
      do 1000 ix=1,iv
      j1 = i1
      isymvx = mult(itypea(iv),itypea(ix))
      isymtu = mult(isymvx,isymr)
      do 200 isyma=1,nirr
      isymb = mult(isyma,isymvx)
      if (isymb .gt. isyma) go to 200
      mca = mcore(isyma)
      maa = mact(isyma)
      mva = mvirt(isyma)
      mcb = mcore(isymb)
      mab = mact(isymb)
      mvb = mvirt(isymb)
      m = maa * (mcb+mvb) + mab * (mca+mva)
      if (m .eq. 0) go to 200
      call intin (q(j1),m)
      mk1(isyma) = j1
      mk2(isyma) = j1 + maa*mcb + mab*mca
      j1 = j1 + m
      if (isymb .eq. isyma) go to 200
      mk1(isymb) = mk1(isyma) + maa*mcb
      mk2(isymb) = mk2(isyma) + mab*mva
200   continue
      ioffpc = ioff1(ilifa(iv)+ix) + 1
      do 400 isymt=1,nirr
      isymu = mult(isymt,isymtu)
      isyma = mult(isymt,isymr)
      isymb = mult(isymu,isymr)
      mca = mcore(isyma)
      mcb = mcore(isymb)
      mva = mvirt(isyma)
      mvb = mvirt(isymb)
      mat = mact(isymt)
      mau = mact(isymu)
c                         t            u            t
c                       -----        -----        -----
c     voc(ti)         i | v |   =  i | I |  x   u | P |
c                       -----        -----        -----
c                     (mca,mat)    (mca,mau)    (mau,mat)
c
c                                  isymu .ge. isyma : 't'
c                                  isymu .lt. isyma : 'n'
c
      if (mca.ne.0. and. mau.ne.0. and. mat.ne.0) then
       if (isymu .ge. isyma) then
        call mxmbn (q(mk1(isymu)),mau,1, pc(ioffpc),1,mau,
     +              vo(iorbti(isymt)+1),1,mca, mca,mau,mat)
       else
        call mxmbn (q(mk1(isymu)),1,mca, pc(ioffpc),1,mau,
     +              vo(iorbti(isymt)+1),1,mca, mca,mau,mat)
       end if
      end if
c                         u            t            u
c                       -----        -----        -----
c     woc(ui)         i | w |   =  i | I |  x   t | P |
c                       -----        -----        -----
c                     (mcb,mau)    (mcb,mat)    (mat,mau)
c
c                                  isymt .ge. isymb : 't'
c                                  isymt .lt. isymb : 'n'
c
      if (isymt .ge. isymb) then
         call mxmb (q(mk1(isymt)),mat,1, pc(ioffpc),mau,1,
     ?              wo(iorbti(isymu)+1),1,mcb, mcb,mat,mau)
      else
         call mxmb (q(mk1(isymt)),1,mcb, pc(ioffpc),mau,1,
     ?              wo(iorbti(isymu)+1),1,mcb, mcb,mat,mau)
      end if
c                         a            t            a
c                       -----        -----        -----
c     voc(au)         u | v |   =  u | P |  x   t | I |
c                       -----        -----        -----
c                     (mau,mvb)    (mau,mat)    (mat,mvb)
c
c                                               isymb .ge. isymt : 't'
c                                               isymb .lt. isymt : 'n'
c
      if (isymb .ge. isymt) then
         call mxmb (pc(ioffpc),1,mau, q(mk2(isymb)),mvb,1,
     ?              vo(iorbat(isymb)+1),1,mau, mau,mat,mvb)
      else
         call mxmb (pc(ioffpc),1,mau, q(mk2(isymb)),1,mat,
     ?              vo(iorbat(isymb)+1),1,mau, mau,mat,mvb)
      end if
c                         a            u            a
c                       -----        -----        -----
c     woc(at)         t | w |   =  t | P |  x   u | I |
c                       -----        -----        -----
c                     (mat,mva)    (mat,mau)    (mau,mva)
c
c                                               isyma .ge. isymu : 't'
c                                               isyma .lt. isymu : 'n'
c
      if (mat.ne.0 .and. mau.ne.0 .and. mva.ne.0) then
       if (isyma .ge. isymu) then
        call mxmbn (pc(ioffpc),mau,1, q(mk2(isyma)),mva,1,
     +              wo(iorbat(isyma)+1),1,mat, mat,mau,mva)
       else
        call mxmbn (pc(ioffpc),mau,1, q(mk2(isyma)),1,mau,
     +              wo(iorbat(isyma)+1),1,mat, mat,mau,mva)
       end if
      endif
400   ioffpc = ioffpc + mau * mat
1000  continue
      call corlsr (ibase)
      return
      end
      subroutine oopart (d,p1,p2,p3,fi,ft,comm,io,jo,q)

c---------------------------------------------------------------------
c   Constructs the orbital-orbital parts of the matrices A and B
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
      integer*4 tabv,tavb,vabi,vaib
      character *1 xn,xt
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      dimension d(*),p1(*),p2(*),p3(*),fi(*),ft(*),comm(*),io(*),
     +          jo(*),q(*)
      dimension iofd(8),ioff(8),mk1(8),mk2(8),mk3(8),mk4(8),mk5(8),
     +          mk6(8),mk7(8)
      dimension ihlp(3,8,50)
c
      data xn, xt / 'n', 't'/
      if (odebug(27).or.odebug(77).or.odebug(86)) then
         do i=1,ncore
         ihlp(1,itype(i),npoint(i))=i
         end do
         do i=nst,nprim
         ihlp(2,itype(i),npoint(i))=i
         end do
         do i=nprimp,nba
         ihlp(3,itype(i),npoint(i))=i
         end do
      end if
c      if (odebug(27))then
c  print*,'p3 matrix'
c  call outden(2,isymr,ioff1,dum,p3)
c      end if
      ibase = icorr(0)
      time = cpulft(1)
      iblo4 = 1
      kc = 1
      ki = 1
      ndens = 0
      nfock = 0
      do 10 isyma=1,nirr
      iofd(isyma) = ndens
      ioff(isyma) = nfock
      ndens = ndens + mact(isyma)**2
10    nfock = nfock + nsymm(isyma)**2
      if (odebug(86)) call outive (iofd,nirr,'iofd')
      if (odebug(86)) call outive (ioff,nirr,'ioff')

c   read in inactive fock matrix

      call rdedx (fi,nfock,iblofi,ned7)
      if (odebug(74)) then
         ioffi = 1
         do 51 i1=1,nirr
         print*,'inactive fock matrix for symmetry ',i1
         call outsqr (fi(ioffi),nsymm(i1),nsymm(i1),nsymm(i1),' ')
51       ioffi = ioffi + nsymm(i1)**2
      end if

c   read in total fock matrix

      call rdedx (ft,nfock,iblofa,ned7)
      if (odebug(74))then
         ioffi=1
         do 52 i1=1,nirr
         print*,'total fock matrix for symmetry ',i1
         call outsqr(ft(ioffi),nsymm(i1),nsymm(i1),nsymm(i1),
     +               'fock matrix ')
52       ioffi=ioffi+nsymm(i1)**2
      end if

c   construct matrices D' and D''

      id1 = icorr(ndens)
      id2 = icorr(ndens)
      if (odebug(86)) print*,'id1,id2 = ',id1,id2
      call fmove (d,q(id1),ndens)
      call fmove (d,q(id2),ndens)
      do 20 isyma=1,nirr
      iyy1 = id1 + iofd(isyma)
      iyy2 = id2 + iofd(isyma)
       do loop = 1, mact(isyma)
        q(iyy1) = q(iyy1) - 1.0d0
        q(iyy2) = q(iyy2) - 2.0d0
        iyy1 = iyy1 + mact(isyma) + 1
        iyy2 = iyy2 + mact(isyma) + 1
       enddo
20    continue
c      if (odebug(27))then
c  do isyma=1,nirr
c  print*,'matrix D for symmetry',isyma
c  call outsqr(d(iofd(isyma) + 1),mact(isyma),mact(isyma),mact(isyma),' ')
c  print*,'matrix D1 for symmetry',isyma
c  call outsqr(q(id1+iofd(isyma)),mact(isyma),mact(isyma),mact(isyma),' ')
c  print*,'matrix D2 for symmetry',isyma
c  call outsqr(q(id2+iofd(isyma)),mact(isyma),mact(isyma),mact(isyma),' ')
c         end do
c      end if

c   read in matrix r

      ir = icorr(ndens)
      call rdedx (q(ir),ndens,iblkr,ned7)

c   reserve space

      i1 = icorr(2*nbasq)     !  (input integrals)
      i2 = icorr(2*nbasq)     !  (output commutators)
      i3 = icorr(nbasq)       !  (scratch space)
      i4 = icorr(nbasq)       !  (scratch space)
      if (odebug(86)) print*,'i1,i2,i3,i4 = ',i1,i2,i3,i4

c   read in integrals

c================================
c   integrals (ij|..), (i.|j.)
c================================

      intpos = 0
      intfil = ned0
      iblf = ibl124
      intmod = 0
      do 1000 i=1,ncore
      isymi = itype(i)
      isym1 = mult(isymi,isymr)
      ma1 = mact(isym1)
      mv1 = mvirt(isym1)
      do 1000 j=1,i
      isymj = itype(j)
      isym2 = mult(isymj,isymr)
      ma2 = mact(isym2)
      mv2 = mvirt(isym2)
      isymij = mult(isymi,isymj)
      j1 = i1

c   coulomb

      do 150 isyma=1,nirr
      isymb = mult(isyma,isymij)
      if (isymb .gt. isyma) go to 150
      maa = mact(isyma)
      mva = mvirt(isyma)
      mab = mact(isymb)
      mvb = mvirt(isymb)
      m = (maa+mva)*(mab+mvb)
      if (m .eq. 0) go to 150
      call intin (q(j1),m)
      mk1(isyma) = j1
      mk2(isyma) = j1+maa*mab
      if (isymb .ne. isyma) mk2(isymb) = mk2(isyma) + mva*mab
      mk3(isyma) = mk2(isyma) + mva*mab + maa*mvb
      j1 = j1 + m
150   continue

c   exchange

      do 200 isyma=1,nirr
      isymb = mult(isyma,isymij)
      maa = mact(isyma)
      mva = mvirt(isyma)
      mab = mact(isymb)
      mvb = mvirt(isymb)
      m = (maa+mva)*(mab+mvb)
      if (m .eq. 0) go to 200
      call intin (q(j1),m)
      mk4(isyma) = j1
      mk5(isyma) = mk4(isyma) + maa*mab
      mk6(isyma) = mk5(isyma) + mva*mab
      mk7(isyma) = mk6(isyma) + maa*mvb
      j1 = j1 + m
200   continue
      itvj = i2
      ibtj = itvj + ma1*ma2
      itbj = ibtj + mv1*ma2
      iabj = itbj + ma1*mv2
      itjv = iabj + mv1*mv2
      ibjt = itjv + ma1*ma2
      itjb = ibjt + mv1*ma2
      iajb = itjb + ma1*mv2
      ma12 = ma1*ma2
      if (isymij .ne. 1) go to 250
      ij = ioff(isymi) + (npoint(j)-1)*nsymm(isymi)+npoint(i)
      fiij = fi(ij)
      ftij = ft(ij)

***********************************************************************
*   commutators it-vj, it-jv (block 1)                                *
***********************************************************************
c----------------------------------------------------------------------
c       sum_x  D'(tx) * (ix|vj)  +  sum_x  (it|xj) * D'(xv)  in buffer|
c----------------------------------------------------------------------
250   call mxmaa (q(id1+iofd(isym1)),1,ma1,q(mk4(isym2)),ma2,1,
     +            q(i4),1,ma1, ma1,ma1,ma2)
      call mxmb  (q(mk4(isym2)),ma2,1,q(id1+iofd(isym2)),1,ma2,
     +            q(i4),1,ma1, ma1,ma2,ma2)
c---------------------------------------------------------------------
c                   sum_x  D'(tx) * (xv|ij)  for  it-vj               |
c---------------------------------------------------------------------
      if (isym1 .ge. isym2) then
         call mxmaa (q(id1+iofd(isym1)),1,ma1,q(mk1(isym1)),1,ma1,
     +               q(itvj),1,ma1, ma1,ma1,ma2)
      else
         call mxmaa (q(id1+iofd(isym1)),1,ma1,q(mk1(isym2)),ma2,1,
     +               q(itvj),1,ma1, ma1,ma1,ma2)
      end if
c---------------------------------------------------------------------
c                   sum_x  (tx|ij) * D'(xv)  for  it-vj              |
c---------------------------------------------------------------------
      if (isym1 .ge. isym2) then
         call mxmb (q(mk1(isym1)),1,ma1,q(id1+iofd(isym2)),1,ma2,
     +              q(itvj),1,ma1, ma1,ma2,ma2)
      else
         call mxmb (q(mk1(isym2)),ma2,1,q(id1+iofd(isym2)),1,ma2,
     +              q(itvj),1,ma1, ma1,ma2,ma2)
      end if
c---------------------------------------------------------------------
c                 - sum_x  (ix|tj) * D'(xv)  for  it-jv              |
c---------------------------------------------------------------------
      if (ma12 .gt. 0) then
      call dgemm(xn,xn,ma1,ma2,ma2
     +     ,-1.0d0,q(mk4(isym1)),ma1
     +     ,q(id1+iofd(isym2)),ma2,0.0d0,q(itjv),ma1)
      end if
c---------------------------------------------------------------------
c                 - sum_x  D'(tx) * (iv|xj)  for  it-jv              |
c---------------------------------------------------------------------
      if (ma1.ne.0 .and. ma2.ne.0) 
     + call mxmbn (q(id1+iofd(isym1)),1,ma1,q(mk4(isym1)),1,ma1,
     +            q(itjv),1,ma1, ma1,ma1,ma2)
c--------------------------------------------------------------------
c                   add/subtract buffer sums                         |
c--------------------------------------------------------------------
      call daxpy(ma12,-2.0d0,q(i4),1,q(itvj),1)
      call daxpy(ma12, 2.0d0,q(i4),1,q(itjv),1)
c--------------------------------------------------------------------
c                   sum_xy  p1(tv,xy) * (ij|xy)  for  it-vj          |
c--------------------------------------------------------------------
      ip = intof0(isymij) + iofsd(isymij,isym1) + 1
      do 280 isyma=1,nirr
      isymb = mult(isyma,isymij)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 280
      if (isyma .le. isymb) then
         call mxmaa (p1(ip),1,npair(isymij), q(mk1(isymb)),1,m, 
     ?               q(i4),1,ma12, ma12,m,1)
      else
         call trnsps (q(mk1(isyma)),q(i3),mact(isyma),mact(isymb))
         call mxmaa (p1(ip),1,npair(isymij), q(i3),1,m, 
     +               q(i4),1,ma12, ma12,m,1)
      end if
      call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(ma12,1.0d0,q(i3),1,q(itvj),1)
      ip = ip + npair(isymij) * m
280   continue
c---------------------------------------------------------------------
c                   sum_xy  p2(tv,xy) * (ix|jy)  for  it-vj          |
c---------------------------------------------------------------------
      ip = intof0(isymij) + iofsd(isymij,isym1) + 1
      do 290 isyma=1,nirr
      isymb = mult(isyma,isymij)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 290
      call mxmaa (p2(ip),1,npair(isymij), q(mk4(isymb)),1,m, 
     ?            q(i4),1,ma12, ma12,m,1)
      call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(ma12,1.0d0,q(i3),1,q(itvj),1)
      ip = ip + npair(isymij) * m
290   continue
c---------------------------------------------------------------------
c                 - sum_xy  p3(tv,xy) * (ix|jy)  for  it-jv          |
c---------------------------------------------------------------------
      ip = intof0(isymij) + iofsd(isymij,isym1) + 1
      do 300 isyma=1,nirr
      isymb = mult(isyma,isymij)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 300
      call mxmaa (p3(ip),1,npair(isymij), q(mk4(isymb)),1,m, 
     ?            q(i4),1,ma12, ma12,m,1)
      call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(ma12,-1.0d0,q(i3),1,q(itjv),1)
      ip = ip + npair(isymij) * m
300   continue

      if (isymij .ne. 1) go to 350
c---------------------------------------------------------------------
c                   fi(ij) * d(tv)                                   |
c---------------------------------------------------------------------
      call daxpy(ma12,fiij,d(iofd(isym1)+1),1,q(itvj),1)
c---------------------------------------------------------------------
c            -  2 * delta(tv) * f(ij)                                |
c---------------------------------------------------------------------
      call daxpy(ma1,-2.0d0,ftij,0,q(itvj),ma1+1)
c---------------------------------------------------------------------
c                 - delta(ij) * r(tv)                                |
c---------------------------------------------------------------------
      if (i .ne. j) go to 350
      call daxpy(ma12,-1.0d0,q(ir+iofd(isym1)),1,q(itvj),1)
c---------------------------------------------------------------------
c             + 2 * delta(ij) * f(tv)                                |
c---------------------------------------------------------------------
      iftv = ioff(isym1) + mcore(isym1)*nsymm(isym1) + mcore(isym1) + 1
      icom = itvj
      do 320 iv=1,ma1
      call daxpy(ma1,2.0d0,ft(iftv),1,q(icom),1)
      iftv = iftv + nsymm(isym1)
320   icom = icom + ma1
350   continue
      if (odebug(27)) then
         write(iwr,351)i,j
351      format(/'commutators it-vj & it-jv for i = ',
     +            i2,', j = ',i2,' :'/)
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write(iwr,2002)ihlp(2,isym1,iz),
     +                  (q(itvj-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write(iwr,2002)ihlp(2,isym1,iz),
     +                  (q(itjv-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
      end if

**********************************************************************
*   commutators it-bj, it-jb, ib-tj, ib-jt (block 2)                 *
**********************************************************************

      if (odebug(77))then
         print*,'integrals (ij|at) for i,j = ',i,j
         call outsqr(q(mk2(isym1)),mv1,mv1,ma2,' ')
         print*,'integrals (it|aj) for i,j = ',i,j
         call outsqr(q(mk5(isym1)),mv1,mv1,ma2,' ')
         print*,'integrals (ia|tj) for i,j = ',i,j
         call outsqr(q(mk6(isym1)),ma1,ma1,mv2,' ')
      end if
c---------------------------------------------------------------------
c               2 * sum_x  D''(tx) * (ix|jb)  for  it-jb             |
c---------------------------------------------------------------------
      if (ma1*mv2 .gt. 0) then
         call dgemm(xn,xt,ma1,mv2,ma1,2.0d0
     +            ,q(id2+iofd(isym1)),ma1
     +            ,q(mk5(isym2)),mv2,0.0d0,q(itjb),ma1)
      end if
c---------------------------------------------------------------------
c                   sum_x  D''(tx) * (ij|xb)  for  it-bj             |
c---------------------------------------------------------------------
      if (isym1.gt.isym2) then
         call mxmaa (q(id2+iofd(isym1)),1,ma1,q(mk2(isym2)),1,ma1,
     +               q(itbj),1,ma1, ma1,ma1,mv2)
      else
         call mxmaa (q(id2+iofd(isym1)),1,ma1,q(mk2(isym2)),mv2,1,
     +               q(itbj),1,ma1, ma1,ma1,mv2)
      end if
c---------------------------------------------------------------------
c             - 2 * sum_x  D''(tx) * (ix|bj)  for  it-bj             |
c---------------------------------------------------------------------
      call daxpy(ma1*mv2,-1.0d0,q(itjb),1,q(itbj),1)
c---------------------------------------------------------------------
c                 - sum_x  D''(tx) * (ib|jx)  for  it-jb             |
c---------------------------------------------------------------------
      if(ma1.ne.0 .and. mv2.ne.0) 
     +  call mxmbn (q(id2+iofd(isym1)),1,ma1,q(mk6(isym1)),1,ma1,
     +              q(itjb),1,ma1, ma1,ma1,mv2)
c---------------------------------------------------------------------
c               2 * sum_x  (ib|xj) * D''(xt)  for  ib-jt             |
c---------------------------------------------------------------------
      if (mv1*ma2 .gt. 0) then
        call dgemm(xt,xn,mv1,ma2,ma2,2.0d0
     +       ,q(mk6(isym2)),ma2
     +       ,q(id2+iofd(isym2)),ma2,0.0d0,q(ibjt),mv1)
      end if
c---------------------------------------------------------------------
c                   sum_x  (ij|bx) * D''(xt)  for  ib-tj             |
c---------------------------------------------------------------------
      if (isym1 .ge. isym2) then
         call mxmaa (q(mk2(isym1)),1,mv1,q(id2+iofd(isym2)),1,ma2,
     +               q(ibtj),1,mv1, mv1,ma2,ma2)
      else
         call mxmaa (q(mk2(isym1)),ma2,1,q(id2+iofd(isym2)),1,ma2,
     +               q(ibtj),1,mv1, mv1,ma2,ma2)
      end if
c---------------------------------------------------------------------
c             - 2 * sum_x  (ib|xj) * D''(xt)  for  ib-tj             |
c---------------------------------------------------------------------
      call daxpy(mv1*ma2,-1.0d0,q(ibjt),1,q(ibtj),1)
c---------------------------------------------------------------------
c                 - sum_x  (ix|bj) * D''(xt)  for  ib-jt             |
c---------------------------------------------------------------------
      if (mv1.ne.0 .and. ma2.ne.0)
     + call mxmbn (q(mk5(isym2)),1,mv1,q(id2+iofd(isym2)),1,ma2,
     +             q(ibjt),1,mv1, mv1,ma2,ma2)
c---------------------------------------------------------------------
c               2 * delta(ij) * f(bt)                                |
c---------------------------------------------------------------------
      if (i .ne. j) go to 450
      ifbt = ioff(isym1) + mcore(isym1)*nsymm(isym1) + nprm(isym1) + 1
      icom = ibtj
      do 400 it=1,ma1
      call daxpy(mv1,2.0d0,ft(ifbt),1,q(icom),1)
      ifbt = ifbt + nsymm(isym1)
400   icom = icom + mv1
450   continue
      if (odebug(27)) then
         if (i .eq. j)go to 452
         write(iwr,451)i,j
451      format(/'commutators it-bj & it-jb for i = ',i2,
     +                 ', j = ',i2,' :'/)
         write(iwr,1003)(ihlp(3,isym2,jz),jz=1,mv2)
         do iz=1,ma1
         write(iwr,2003)ihlp(2,isym1,iz),
     +                  (q(itbj-1+(jz-1)*ma1+iz),jz=1,mv2)
         end do
         write(iwr,1003)(ihlp(3,isym2,jz),jz=1,mv2)
         do iz=1,ma1
         write(iwr,2003)ihlp(2,isym1,iz),
     +                  (q(itjb-1+(jz-1)*ma1+iz),jz=1,mv2)
         end do
452      write(iwr,454)i,j
454      format(/'commutators ib-tj & ib-jt for i = ',i2,
     +                 ', j = ',i2,' :'/)
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,mv1
         write(iwr,2002)ihlp(3,isym1,iz),
     +                  (q(ibtj-1+(jz-1)*mv1+iz),jz=1,ma2)
         end do
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,mv1
         write(iwr,2002)ihlp(3,isym1,iz),
     +                  (q(ibjt-1+(jz-1)*mv1+iz),jz=1,ma2)
         end do
      end if

**********************************************************************
*   commutators ia-bj, ia-jb (block 4)                               *
**********************************************************************
c---------------------------------------------------------------------
c                   move (ab|ij) to ia-bj                            |
c---------------------------------------------------------------------
      if (isym1 .ge. isym2) then
         call fmove  (q(mk3(isym1)),q(iabj),mv1*mv2)
      else
         call trnsps (q(mk3(isym2)),q(iabj),mv2,mv1)
      end if
c---------------------------------------------------------------------
c                   transpose (ai|bj)=(ia|bj) to ia-jb               |
c---------------------------------------------------------------------
      call trnsps (q(mk7(isym2)),q(iajb),mv2,mv1)
c---------------------------------------------------------------------
c                   add (ai|bj) to ia-bj                             |
c---------------------------------------------------------------------
      call daxpy(mv1*mv2,-2.0d0,q(iajb),1,q(iabj),1)
c---------------------------------------------------------------------
c                   add (aj|bi)=(ib|aj) to ia-jb                     |
c---------------------------------------------------------------------
      call daxpy(mv1*mv2,-0.5d0,q(mk7(isym1)),1,q(iajb),1)
      call dscal(mv1*mv2,-2.0d0,q(iabj),1)
      call dscal(mv1*mv2,-4.0d0,q(iajb),1)
c---------------------------------------------------------------------
c             - 2 * delta(ab) * f(ij)                                |
c---------------------------------------------------------------------
      if (isymij .ne. 1) go to 550
      call daxpy(mv1,-2.0d0,ftij,0,q(iabj),mv1+1)
c---------------------------------------------------------------------
c             + 2 * delta(ij) * f(ab)                                |
c---------------------------------------------------------------------
      if (i .ne. j) go to 550
      ifab = ioff(isym1) + nprm(isym1)*nsymm(isym1) + nprm(isym1) + 1
      icom = iabj
      do 500 ib=1,mv1
      call daxpy(mv1,2.0d0,ft(ifab),1,q(icom),1)
      ifab = ifab + nsymm(isym1)
500   icom = icom + mv1
550   continue
      if (odebug(27)) then
         write(iwr,551)i,j
551      format(/'commutators ia-bj & ia-jb for i = ',i2,
     +           ', j = ',i2,' :'/)
         write(iwr,1003)(ihlp(3,isym2,jz),jz=1,mv2)
         do iz=1,mv1
         write(iwr,2003)ihlp(3,isym1,iz),
     +                  (q(iabj-1+(jz-1)*mv1+iz),jz=1,mv2)
         end do
         write(iwr,1003)(ihlp(3,isym2,jz),jz=1,mv2)
         do iz=1,mv1
         write(iwr,2003)ihlp(3,isym1,iz),
     +                  (q(iajb-1+(jz-1)*mv1+iz),jz=1,mv2)
         end do
      end if

c   output of commutators

      if (i .eq. j) go to 700
      call dcopy(ma12,q(itvj),1,comm(kc  ),2)
      call dcopy(ma12,q(itjv),1,comm(kc+1),2)
      kc = kc + 2*ma12
      call dcopy(ma1*mv2,q(itbj),1,comm(kc  ),2)
      call dcopy(ma1*mv2,q(itjb),1,comm(kc+1),2)
      kc = kc + 2*ma1*mv2
      call dcopy(mv1*ma2,q(ibtj),1,comm(kc  ),2)
      call dcopy(mv1*ma2,q(ibjt),1,comm(kc+1),2)
      kc = kc + 2*mv1*ma2
      call dcopy(mv1*mv2,q(iabj),1,comm(kc  ),2)
      call dcopy(mv1*mv2,q(iajb),1,comm(kc+1),2)
      kc = kc + 2*mv1*mv2
      if (kc.gt.2*(nbasq+204)) call caserr ('comm exceeds array bounds')

c   indices it-vj

      iad2 = iorbti(isym2) + npoint(j)
      do 610 iv=1,ma2
      iad1 = iorbti(isym1) + npoint(i)
      do 600 it=1,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
600   iad1 = iad1 + mcore(isymi)
610   iad2 = iad2 + mcore(isymj)

c   indices it-bj

      iad2 = iorbai(isym2) + npoint(j)
      do 630 ib=1,mv2
      iad1 = iorbti(isym1) + npoint(i)
      do 620 it=1,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
620   iad1 = iad1 + mcore(isymi)
630   iad2 = iad2 + mcore(isymj)

c   indices ib-tj

      iad2 = iorbti(isym2) + npoint(j)
      do 650 it=1,ma2
      iad1 = iorbai(isym1) + npoint(i)
      do 640 ib=1,mv1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
640   iad1 = iad1 + mcore(isymi)
650   iad2 = iad2 + mcore(isymj)

c   indices ia-bj

      iad2 = iorbai(isym2) + npoint(j)
      do 670 ib=1,mv2
      iad1 = iorbai(isym1) + npoint(i)
      do 660 ia=1,mv1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
660   iad1 = iad1 + mcore(isymi)
670   iad2 = iad2 + mcore(isymj)
      go to 800
700   continue
      do 710 iv=1,ma1
      call dcopy(ma1+1-iv,q(itvj-1+(iv-1)*ma1+iv),1,comm(kc  ),2)
      call dcopy(ma1+1-iv,q(itjv-1+(iv-1)*ma1+iv),1,comm(kc+1),2)
710   kc = kc + 2*(ma1+1-iv)
      call dcopy(mv1*ma1,q(ibtj),1,comm(kc  ),2)
      call dcopy(mv1*ma1,q(ibjt),1,comm(kc+1),2)
      kc = kc + 2*mv1*ma1
      do 720 ib=1,mv1
      call dcopy(mv1+1-ib,q(iabj-1+(ib-1)*mv1+ib),1,comm(kc  ),2)
      call dcopy(mv1+1-ib,q(iajb-1+(ib-1)*mv1+ib),1,comm(kc+1),2)
720   kc = kc + 2*(mv1+1-ib)
      if (kc.gt.2*(nbasq+204)) call caserr ('comm exceeds array bounds')

c   indices it-vi

      iad2 = iorbti(isym1) + npoint(i)
      do 740 iv=1,ma1
      iad1 = iorbti(isym1) + (iv-1)*mcore(isymi) + npoint(i)
      do 730 it=iv,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
730   iad1 = iad1 + mcore(isymi)
740   iad2 = iad2 + mcore(isymi)

c   indices ib-ti

      iad2 = iorbti(isym1) + npoint(i)
      do 760 it=1,ma1
      iad1 = iorbai(isym1) + npoint(i)
      do 750 ib=1,mv1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
750   iad1 = iad1 + mcore(isymi)
760   iad2 = iad2 + mcore(isymi)

c   indices ia-bi

      iad2 = iorbai(isym1) + npoint(i)
      do 780 ib=1,mv1
      iad1 = iorbai(isym1) + (ib-1)*mcore(isymi) + npoint(i)
      do 770 ia=ib,mv1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
770   iad1 = iad1 + mcore(isymi)
780   iad2 = iad2 + mcore(isymi)

c   output

800   oend = .false.
      call comout (kc,ki,comm,io,jo,oend)
1000  continue

c================================
c   integrals (ai|..), (a.|i.)
c================================

      intpos = 0
      intfil = ned0
      iblf = ibl356
      intmod = 0
      do 2010 isymci=1,nirr
      do 2010 isymc=1,nirr
      isymi = mult(isymci,isymc)
      isym1 = mult(isymi,isymr)
      isym2 = mult(isymc,isymr)
      ma1 = mact(isym1)
      ma2 = mact(isym2)
      ma12 = ma1*ma2
c      if (odebug(86)) write (iwr,862) ma1,ma2,isym1,isym2
c862   format ('ma1 = ',i3,3x,'ma2 = ',i3,3x,'isym1 = ',i1,3x,'isym2 = ',i1)
      do 2010 ii=1,mcore(isymi)
      do 2000 ia=1,mvirt(isymc)
      j1 = i1

c   coulomb

      do 1200 isyma=1,nirr
      isymb = mult(isyma,isymci)
      if (isymb .gt. isyma) go to 1200
      maa = mact(isyma)
      mab = mact(isymb)
      m = maa * mab
      if (m .eq. 0) go to 1200
      mk1(isyma) = j1
      if (isyma .gt. isymb) then
         call intin (q(j1),m)
      else
         call intin (q(i3),(maa*(maa+1))/2)
         call square (q(j1),q(i3),maa,maa)
      end if
      if (odebug(77)) then
         print*,'integrals (ai|tu) for a,i = ',
     ?          ihlp(3,isymc,ia),ihlp(1,isymi,ii),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mact(isymb),' ')
      end if
      j1 = j1 + m
1200  continue

c   exchange

      do 1300 isyma=1,nirr
      isymb = mult(isyma,isymci)
      maa = mact(isyma)
      mab = mact(isymb)
      m = maa * mab
      if (m .eq. 0) go to 1300
      mk2(isyma) = j1
      call intin (q(j1),m)
      if (odebug(77)) then
         print*,'integrals (at|iu) for a,i = ',
     ?          ihlp(3,isymc,ia),ihlp(1,isymi,ii),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mact(isymb),' ')
      end if
      j1 = j1 + m
1300  continue

      if (ma12 .eq. 0) go to 2000
      itbv = i2
      itvb = itbv + ma12

**********************************************************************
*   commutators it-bv, it-vb (block 3)                               *
**********************************************************************
c---------------------------------------------------------------------
c               2 * sum_x  (bx|it) * D(xv)  for  it-bv               |
c---------------------------------------------------------------------
      call dgemm(xt,xn,ma1,ma2,ma2,2.0d0
     +        ,q(mk2(isym2)),ma2
     +        ,d(iofd(isym2)+1),ma2,0.0d0,q(itbv),ma1)
c      call mxmaa (q(mk2(isym2)),ma2,1, d(iofd(isym2)+1),1,ma2,
c     ?            q(itbv),1,ma1, ma1,ma2,ma2)
c      call dscal(ma12,2.0d0,q(itbv),1)
c---------------------------------------------------------------------
c                   sum_x  (bi|tx) * D(xv)  for  it-vb               |
c---------------------------------------------------------------------
      if (isym1 .ge. isym2) then
         call mxmaa (q(mk1(isym1)),1,ma1, d(iofd(isym2)+1),1,ma2, 
     +               q(itvb),1,ma1, ma1,ma2,ma2)
      else
         call mxmaa (q(mk1(isym2)),ma2,1, d(iofd(isym2)+1),1,ma2, 
     +               q(itvb),1,ma1, ma1,ma2,ma2)
      end if
c---------------------------------------------------------------------
c             - 2 * sum_x  (bx|it) * D(xv)  for  it-vb               |
c---------------------------------------------------------------------
      call daxpy(ma12,-1.0d0,q(itbv),1,q(itvb),1)
c---------------------------------------------------------------------
c                 - sum_x  (bt|ix) * D(xv)  for  it-bv               |
c---------------------------------------------------------------------
      if (ma1.ne.0 .and. ma2.ne.0) 
     + call mxmbn (q(mk2(isym1)),1,ma1, d(iofd(isym2)+1),1,ma2, 
     +             q(itbv),1,ma1, ma1,ma2,ma2)
c---------------------------------------------------------------------
c                   sum_xy  p1(tv,xy) * (bi|xy)  for  it-vb          |
c---------------------------------------------------------------------
      ip = intof0(isymci) + iofsd(isymci,isym1) + 1
      do 1450 isyma=1,nirr
      isymb = mult(isyma,isymci)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 1450
      if (isyma .le. isymb) go to 1350
      call trnsps (q(mk1(isyma)),q(i3),mact(isyma),mact(isymb))
      call mxmaa (p1(ip),1,npair(isymci), q(i3),1,m, 
     +            q(i4),1,ma12, ma12,m,1)
      go to 1400
1350  call mxmaa (p1(ip),1,npair(isymci), q(mk1(isymb)),1,m, 
     +            q(i4),1,ma12, ma12,m,1)
1400  call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(ma12,1.0d0,q(i3),1,q(itvb),1)
      ip = ip + npair(isymci) * m
1450  continue
c      if (odebug(27))then
c  print*,'p1 sums for it-vb, i = ',ii,' b = ',ia
c  call outsqr(q(i3),ma1,ma1,ma2,' ')
c      end if
c---------------------------------------------------------------------
c                   sum_xy  p2(vt,xy) * (bx|iy)  for  it-vb          |
c---------------------------------------------------------------------
c      if (odebug(27))call vclr(gin,1,ma12)
      ip = intof0(isymci) + iofsd(isymci,isym2) + 1
      do 1500 isyma=1,nirr
      isymb = mult(isyma,isymci)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 1500
      call trnsps (q(mk2(isyma)),q(i3),mact(isyma),mact(isymb))
c      if (odebug(27))then
c  print*,'p2 sums for it-vb, symmetry ',isyma,' --- p2 matrix:'
c  call outsqr(p2(ip),npair(isymci),ma12,m,' ')
c  print*,'p2 sums for it-vb, symmetry ',isyma,' --- integrals:'
c  call outsqr(q(i3),mact(isymb),mact(isymb),mact(isyma),' ')
c      end if
      call mxmb (p2(ip),1,npair(isymci), q(i3),1,m, 
     +           q(itvb),1,ma12, ma12,m,1)
c      if (odebug(27))call mxmb (p2(ip),1,npair(isymci),q(i3),1,m,gin,1,ma12,ma12,m,1)
      ip = ip + npair(isymci) * m
1500  continue
c      if (odebug(27))then
c  print*,'p2 sums for it-vb, i = ',ii,' b = ',ia
c  call outsqr(gin,ma1,ma1,ma2,' ')
c      end if
c---------------------------------------------------------------------
c                 - sum_xy  p3(vt,xy) * (bx|iy)  for  it-bv          |
c---------------------------------------------------------------------
      ip = intof0(isymci) + iofsd(isymci,isym2) + 1
c      if (odebug(27))call vclr(gin,1,ma12)
      do 1550 isyma=1,nirr
      isymb = mult(isyma,isymci)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 1550
      call trnsps (q(mk2(isyma)),q(i3),mact(isyma),mact(isymb))
c      if (odebug(27))then
c  print*,'p3 sums for it-vb, symmetry ',isyma,' --- p3 matrix:'
c  call outsqr(p3(ip),npair(isymci),ma12,m,' ')
c  print*,'p3 sums for it-vb, symmetry ',isyma,' --- integrals:'
c  call outsqr(q(i3),mact(isymb),mact(isymb),mact(isyma),' ')
c      end if
      if (ma12.ne.0) 
     + call mxmbn (p3(ip),1,npair(isymci), q(i3),1,m, 
     +            q(itbv),1,ma12, ma12,m,1)
c      if (odebug(27))
c    +call mxmb (p3(ip),1,npair(isymci),q(i3),1,m,gin,1,ma12, ma12,m,1)
      ip = ip + npair(isymci) * m
1550  continue
c      if (odebug(27))then
c  print*,'p3 sums for it-vb, i = ',ii,' b = ',ia
c  call outsqr(gin,ma1,ma1,ma2,' ')
c      end if
c---------------------------------------------------------------------
c                   fi(bi) * d(tv)                                   |
c---------------------------------------------------------------------
      if (isymci .ne. 1) go to 1600
      fibi = fi(ioff(isymc)+(ii-1)*nsymm(isymc)+nprm(isymc)+ia)
      call daxpy(ma12,fibi,d(iofd(isym1)+1),1,q(itvb),1)
1600  continue
      if (odebug(27)) then
         write (iwr,1601) ihlp(1,isymi,ii),ihlp(3,isymc,ia)
1601     format (/'commutators it-bv & it-vb for i = ',i2,
     +            ', b = ',i2,' :'/)
         write (iwr,1002) (ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write (iwr,2002) ihlp(2,isym1,iz),
     +                    (q(itbv-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
         write(iwr,1002) (ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write(iwr,2002) ihlp(2,isym1,iz),
     +                   (q(itvb-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
      end if

c   output of commutators

      call dcopy(ma12,q(itbv),1,comm(kc  ),2)
      call dcopy(ma12,q(itvb),1,comm(kc+1),2)
      kc = kc + 2 * ma12
      if (kc.gt.2*(nbasq+204)) call caserr ('comm exceeds array bounds')

c   indices it-bv

      iad2 = iorbat(isymc) + (ia-1)*ma2 + 1
      do 1700 iv=1,ma2
      iad1 = iorbti(isym1) + npoint(i)
      do 1650 it=1,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
1650  iad1 = iad1 + mcore(isymi)
1700  iad2 = iad2 + 1
      oend = .false.
      call comout (kc,ki,comm,io,jo,oend)
2000  continue
2010  continue

c================================
c   integrals (ab|..), (a.|b.)
c================================

      do 4000 isymcd=1,nirr
      do 4000 isymc=1,nirr
      isymd = mult(isymcd,isymc)
      if (isymd .gt. isymc) go to 4000
      isym1 = mult(isymc,isymr)
      isym2 = mult(isymd,isymr)
      mc1 = mcore(isym1)
      mc2 = mcore(isym2)
      ma1 = mact(isym1)
      ma2 = mact(isym2)
      maa = ma1*ma2
      mca = mc1*ma2
      mac = ma1*mc2
      ib1 = mvirt(isymd)
      do 3000 ia=1,mvirt(isymc)
      if (isymcd .eq. 1) ib1 = ia
      do 3000 ib=1,ib1
      j1 = i1
c...integrals (ab|ti)
      do 2050 isyma=1,nirr
      isymb = mult(isymcd,isyma)
      m = mact(isyma) * mcore(isymb)
      if (m .eq. 0) go to 2050
      mk1(isyma) = j1
      call intin (q(j1),m)
      if (odebug(77)) then
         print*,'integrals (ab|ti) for a,b = ',
     ?           ihlp(3,isymc,ia),ihlp(3,isymd,ib),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mcore(isymb),' ')
      end if
      j1 = j1 + m
2050  continue

c   integrals (at|bi)

      do 2100 isyma=1,nirr
      isymb = mult(isymcd,isyma)
      m = mact(isyma) * mcore(isymb)
      if (m .eq. 0) go to 2100
      mk2(isyma) = j1
      call intin (q(j1),m)
      if (odebug(77)) then
         print*,'integrals (at|bi) for a,b = ',
     ?          ihlp(3,isymc,ia),ihlp(3,isymd,ib),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mcore(isymb),' ')
      end if
      j1 = j1 + m
2100  continue

c   integrals (ai|bt)

      do 2150 isyma=1,nirr
      isymb = mult(isymcd,isyma)
      m = mact(isyma) * mcore(isymb)
      if (m .eq. 0) go to 2150
      mk3(isyma) = j1
      call intin (q(j1),m)
      if (odebug(77)) then
         print*,'integrals (ai|bt) for a,b = ',
     ?          ihlp(3,isymc,ia),ihlp(3,isymd,ib),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mcore(isymb),' ')
      end if
      j1 = j1 + m
2150  continue

c   integrals (ab|tu)

      do 2300 isyma=1,nirr
      isymb = mult(isyma,isymcd)
      if (isymb .gt. isyma) go to 2300
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 2300
      mk4(isyma) = j1
      if (isyma .gt. isymb) then
         call intin (q(j1),m)      
      else
         call intin (q(i3),(mact(isyma)*(mact(isyma)+1))/2)
         call square (q(j1),q(i3),mact(isyma),mact(isyma))
      end if
      if (odebug(77)) then
         print*,'integrals (ab|tu) for a,b = ',
     ?          ihlp(3,isymc,ia),ihlp(3,isymd,ib),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mact(isymb),' ')
      end if
      j1 = j1 + m
2300  continue

c   integrals (at|bu)

      do 2350 isyma=1,nirr
      isymb = mult(isymcd,isyma)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 2350
      mk5(isyma) = j1
      call intin (q(j1),m)
      if (odebug(77)) then
         print*,'integrals (at|bu) for a,b = ',
     ?          ihlp(3,isymc,ia),ihlp(3,isymd,ib),'   symmetry ',isyma
         call outsqr (q(j1),mact(isyma),mact(isyma),mact(isymb),' ')
      end if
      j1 = j1 + m
2350  continue

      iabv = i2
      iavb = iabv + mca
      vabi = iavb + mca
      vaib = vabi + mac
      tabv = vaib + mac
      tavb = tabv + maa

**********************************************************************
*   commutators ia-bv, ia-vb, va-bi, va-ib (block 5)                 *
**********************************************************************
c---------------------------------------------------------------------
c                   sum_x  (ax|bi) * D(xv)  for  ia-vb               |
c---------------------------------------------------------------------
      call mxmaa (q(mk2(isym2)),ma2,1,d(iofd(isym2)+1),1,ma2,
     +            q(iavb),1,mc1, mc1,ma2,ma2)
c---------------------------------------------------------------------
c               2 * sum_x  (ai|bx) * D(xv)  for  ia-bv               |
c---------------------------------------------------------------------
      if (mca .gt. 0) then
         call dgemm(xt,xn,mc1,ma2,ma2,2.0d0
     +        ,q(mk3(isym2)),ma2
     +        ,d(iofd(isym2)+1),ma2,0.0d0,q(iabv),mc1)
      end if
c---------------------------------------------------------------------
c             - 2 * sum_x  (ai|bx) * D(xv)  for  ia-vb               |
c---------------------------------------------------------------------
      call daxpy(mca,-1.0d0,q(iabv),1,q(iavb),1)
c---------------------------------------------------------------------
c                 - sum_x  (ab|xi) * D(xv)  for  ia-bv               |
c---------------------------------------------------------------------
      if (mc1.ne.0 .and. ma2.ne.0)
     + call mxmbn (q(mk1(isym2)),ma2,1,d(iofd(isym2)+1),1,ma2,
     +            q(iabv),1,mc1, mc1,ma2,ma2)
c---------------------------------------------------------------------
c                   sum_x  D(vx) * (ai|bx)  for  va-ib               |
c---------------------------------------------------------------------
      call mxmaa (d(iofd(isym1)+1),1,ma1,q(mk3(isym1)),1,ma1,
     +            q(vaib),1,ma1, ma1,ma1,mc2)
c---------------------------------------------------------------------
c               2 * sum_x  D(vx) * (ax|bi)  for  va-bi               |
c---------------------------------------------------------------------
      if (mac .gt. 0) then
         call dgemm(xn,xn,ma1,mc2,ma1,2.0d0
     +         ,d(iofd(isym1)+1),ma1
     +         ,q(mk2(isym1)),ma1,0.0d0,q(vabi),ma1)
      end if
c---------------------------------------------------------------------
c             - 2 * sum_x  D(vx) * (ax|bi)  for  va-ib               |
c---------------------------------------------------------------------
      call daxpy(mac,-1.0d0,q(vabi),1,q(vaib),1)
c---------------------------------------------------------------------
c                 - sum_x  D(vx) * (ab|xi)  for  va-bi               |
c---------------------------------------------------------------------
      if (mc2.ne.0 .and. ma1.ne.0)
     + call mxmbn (d(iofd(isym1) + 1),1,ma1,q(mk1(isym1)),1,ma1,
     +             q(vabi),1,ma1, ma1,ma1,mc2)
c---------------------------------------------------------------------
c             - 2 * delta(ab) * f(vi)  (only to ia-bv !)             |
c---------------------------------------------------------------------
      if (isymcd .ne. 1.or.ia .ne. ib) go to 2450
      ifiv = ioff(isym1)+mc1*nsymm(isym1) + 1
      icom = iabv
      do 2400 iv=1,ma1
      call daxpy(mc1,-2.0d0,ft(ifiv),1,q(icom),1)
      ifiv = ifiv+nsymm(isym1)
2400  icom = icom+mc1
2450  continue

**********************************************************************
*   commutators ta-bv, ta-vb (block 6)                               *
**********************************************************************
      call vclr(q(tabv),1,maa)
      call vclr(q(tavb),1,maa)
c---------------------------------------------------------------------
c                   sum_xy  p1(tv,xy) * (ab|xy)  for  ta-bv          |
c---------------------------------------------------------------------
      ip = intof0(isymcd) + iofsd(isymcd,isym1) + 1
      do 2600 isyma=1,nirr
      isymb = mult(isyma,isymcd)
      m = mact(isyma)*mact(isymb)
      if (m .eq. 0) go to 2600
      if (isyma.le.isymb) go to 2500
      call trnsps (q(mk4(isyma)),q(i3),mact(isyma),mact(isymb))
      call mxmaa (p1(ip),1,npair(isymcd),q(i3),1,m,q(i4),1,maa,maa,m,1)
      go to 2550
2500  call mxmaa (p1(ip),1,npair(isymcd),q(mk4(isymb)),1,m,
     +            q(i4),1,maa,maa,m,1)
2550  call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(maa,1.0d0,q(i3),1,q(tabv),1)
      ip = ip + npair(isymcd) * m
2600  continue
c---------------------------------------------------------------------
c                   sum_xy  p2(tv,xy) * (ax|by)  for  ta-bv          |
c---------------------------------------------------------------------
      ip = intof0(isymcd) + iofsd(isymcd,isym1) + 1
      do 2650 isyma=1,nirr
      isymb = mult(isyma,isymcd)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 2650
      call trnsps (q(mk5(isyma)),q(i3),mact(isyma),mact(isymb))
      call mxmaa (p2(ip),1,npair(isymcd),q(i3),1,m,q(i4),1,maa,maa,m,1)
      call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(maa,1.0d0,q(i3),1,q(tabv),1)
      ip = ip + npair(isymcd) * m
2650  continue
c---------------------------------------------------------------------
c                 - sum_xy  p3(tv,xy) * (ax|by)  for  ta-vb          |
c---------------------------------------------------------------------
      ip = intof0(isymcd) + iofsd(isymcd,isym1) + 1
      do 2700 isyma=1,nirr
      isymb = mult(isyma,isymcd)
      m = mact(isyma) * mact(isymb)
      if (m .eq. 0) go to 2700
      call trnsps (q(mk5(isyma)),q(i3),mact(isyma),mact(isymb))
c      if (odebug(27))then
c  print*,'p3 sums for ta-vb, symmetry ',isyma,' --- p3 matrix:'
c  call outsqr(p3(ip),npair(isymcd),maa,m,' ')
c  print*,'p3 sums for ta-vb, symmetry ',isyma,' --- integrals:'
c  call outsqr(q(i3),mact(isymb),mact(isymb),mact(isyma),' ')
c      end if
      call mxmaa (p3(ip),1,npair(isymcd),q(i3),1,m,q(i4),1,maa,maa,m,1)
      call trnsps (q(i4),q(i3),ma2,ma1)
      call daxpy(maa,-1.0d0,q(i3),1,q(tavb),1)
      ip = ip + npair(isymcd) * m
2700  continue
c---------------------------------------------------------------------
c                   fi(ab) * d(tv)                                   |
c---------------------------------------------------------------------
      if (isymcd .ne. 1) go to 2705
      fiab = fi(ioff(isymc)+(nprm(isymc)+ia-1)*nsymm(isymc)
     +                     +nprm(isymc)+ib)
      call daxpy(maa,fiab,d(iofd(isym1)+1),1,q(tabv),1)
c---------------------------------------------------------------------
c                 - delta(ab) * r(tv)                                |
c---------------------------------------------------------------------
      if (ia .eq. ib) call daxpy(maa
     +    ,-1.0d0,q(ir+iofd(isym1)),1
     +    ,q(tabv),1)
2705  continue
      if (odebug(27)) then
         write(iwr,2706)ihlp(3,isymc,ia),ihlp(3,isymd,ib)
2706     format(/'commutators ia-bv & ia-vb for a = ',i2,
     +           ', b = ',i2,' :'/)
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,mc1
         write(iwr,2002)ihlp(1,isym1,iz),
     +                  (q(iabv-1+(jz-1)*mc1+iz),jz=1,ma2)
         end do
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,mc1
         write(iwr,2002)ihlp(1,isym1,iz),
     +                  (q(iavb-1+(jz-1)*mc1+iz),jz=1,ma2)
         end do
         if (isymcd .eq. 1.and.ia .eq. ib)go to 27075
         write(iwr,2707)ihlp(3,isymc,ia),ihlp(3,isymd,ib)
2707     format(/'commutators va-bi & va-ib for a = ',i2,
     +           ', b = ',i2,' :'/)
         write(iwr,1001)(ihlp(1,isym2,jz),jz=1,mc2)
         do iz=1,ma1
         write(iwr,2001)ihlp(2,isym1,iz),
     +                  (q(vabi-1+(jz-1)*ma1+iz),jz=1,mc2)
         end do
         write(iwr,1001)(ihlp(1,isym2,jz),jz=1,mc2)
         do iz=1,ma1
         write(iwr,2001)ihlp(2,isym1,iz),
     +                  (q(vaib-1+(jz-1)*ma1+iz),jz=1,mc2)
         end do
27075    write(iwr,2708)ihlp(3,isymc,ia),ihlp(3,isymd,ib)
2708     format(/'commutators ta-bv & ta-vb for a = ',i2,
     +           ', b = ',i2,' :'/)
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write(iwr,2002)ihlp(2,isym1,iz),
     +                  (q(tabv-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
         write(iwr,1002)(ihlp(2,isym2,jz),jz=1,ma2)
         do iz=1,ma1
         write(iwr,2002)ihlp(2,isym1,iz),
     +                  (q(tavb-1+(jz-1)*ma1+iz),jz=1,ma2)
         end do
      end if

c   output of commutators

      if (isymcd.eq.1 .and. ia.eq.ib) go to 2800
      call dcopy(mca,q(iabv),1,comm(kc  ),2)
      call dcopy(mca,q(iavb),1,comm(kc+1),2)
      kc = kc + 2 * mca
      call dcopy(mac,q(vabi),1,comm(kc  ),2)
      call dcopy(mac,q(vaib),1,comm(kc+1),2)
      kc = kc + 2 * mac
      call dcopy(maa,q(tabv),1,comm(kc  ),2)
      call dcopy(maa,q(tavb),1,comm(kc+1),2)
      kc = kc + 2 * maa
      if (kc.gt.2*(nbasq+204)) call caserr ('comm exceeds array bounds')

c   indices ia-bv

      iad2 = iorbat(isymd) + (ib-1)*ma2 + 1
      do 2720 iv=1,ma2
      iad1 = iorbai(isymc) + (ia-1)*mc1 + 1
      do 2710 ii=1,mc1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
2710  iad1 = iad1 + 1
2720  iad2 = iad2 + 1

c   indices va-bi

      iad2 = iorbai(isymd) + (ib-1)*mc2 + 1
      do 2740 ii=1,mc2
      iad1 = iorbat(isymc) + (ia-1)*ma1 + 1
      do 2730 iv=1,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
2730  iad1 = iad1 + 1
2740  iad2 = iad2 + 1

c   indices ta-bv

      iad2 = iorbat(isymd) + (ib-1)*ma2 + 1
      do 2760 iv=1,ma2
      iad1 = iorbat(isymc) + (ia-1)*ma1 + 1
      do 2750 it=1,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
2750  iad1 = iad1 + 1
2760  iad2 = iad2 + 1
      go to 2900
2800  continue
      call dcopy(mca,q(iabv),1,comm(kc  ),2)
      call dcopy(mca,q(iavb),1,comm(kc+1),2)
      kc = kc + 2 * mca
      do 2810 iv=1,ma1
      call dcopy(ma1+1-iv,q(tabv-1+(iv-1)*ma1+iv),1,comm(kc  ),2)
      call dcopy(ma1+1-iv,q(tavb-1+(iv-1)*ma1+iv),1,comm(kc+1),2)
2810  kc = kc + 2*(ma1+1-iv)
      if (kc.gt.2*(nbasq+204)) call caserr ('comm exceeds array bounds')

c   indices ia-av

      iad2 = iorbat(isymc) + (ia-1)*ma1 + 1
      do 2830 iv=1,ma1
      iad1 = iorbai(isymc) + (ia-1)*mc1 + 1
      do 2820 ii=1,mc1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
2820  iad1 = iad1 + 1
2830  iad2 = iad2 + 1

c   indices ta-av

      iad2 = iorbat(isymc) + (ia-1)*ma1 + 1
      do 2850 iv=1,ma1
      iad1 = iorbat(isymc) + (ia-1)*ma1 + iv
      do 2840 it=iv,ma1
      io(ki) = iad1
      jo(ki) = iad2
      ki = ki + 1
2840  iad1 = iad1 + 1
2850  iad2 = iad2 + 1

c   output

2900  oend = .false.
      call comout (kc,ki,comm,io,jo,oend)
3000  continue
4000  continue

      oend = .true.
      call comout (kc,ki,comm,io,jo,oend)

      write (iwr,9999) cpulft(1)-time
9999  format ('construction of commutators ',26('.'),f9.4,' seconds')

      call corlsr (ibase)

1001  format(4x,100i14)
1002  format(4x,100i14)
1003  format(4x,100i14)
2001  format(i2,2x,100f14.9)
2002  format(i2,2x,100f14.9)
2003  format(i2,2x,100f14.9)

      return
      end
      subroutine oosolv (civec,d,x,s,e,frqncy,thresh,sred,ered,
     +                   bred,ared,
     ?                   f,fred1,fred2,iwork,q)
c
c---------------------------------------------------------------------
c   Solves the equation system 
c
c   (          | Sigma(oo)      0     |   | A(oo)  B(oo) | ) ( y(oo) )   ( v )
c   ( frqncy * |                      | - |              | ) (       ) = (   )
c   (          |     0     -Sigma(oo) |   | B(oo)  A(oo) | ) ( z(oo) )   ( w )
c
c   for the optimal orbital trial vector algorithm.
c
c   On entry, x contains trial solution, on return, x contains converged 
c   solution.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      logical odebug
      common /debug/ odebug(100)
c
c
      integer iblbo, iblso, ibleo
      common /optima/ iblbo(50),iblso(50),ibleo(50)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension civec(*),d(*),x(*),s(*),e(*),q(*)
      dimension sred(maxro,maxro),ered(maxro,maxro),bred(maxro)
      dimension ared(2*maxro,2*maxro),f(*),fred1(*),fred2(*)
      dimension iwork(2*maxro)
c
      time = cpulft(1)
      ntot = 1
10    nitero = nitero + 1
      if (nitero .gt. maxito) go to 1000
c
c   fill up the matrices sred,ered
c
      if (odebug(43)) then
         call outive (iblbo,maxro,'oosolv: b block addresses at start')
         call outive (iblso,maxro,'oosolv: s block addresses at start')
         call outive (ibleo,maxro,'oosolv: e block addresses at start')
      end if
      do 20 j=ntot,1,-1
      ix1 = 1
      ix2 = norb + 1
      if (j .lt. ntot) go to 14
      fred1(ntot)       = ddot(norb,x(ix1),1,f(ix1),1) + 
     ?                    ddot(norb,x(ix2),1,f(ix2),1)
      fred1(maxro+ntot) = ddot(norb,x(ix1),1,f(ix2),1) + 
     ?                    ddot(norb,x(ix2),1,f(ix1),1)
      call lrmult (1,.true.,civec,d,x,s,e,q,q)
      if (odebug(30)) call outvec (e,norb2,'e:')
      bred(ntot)      = ddot(norb,x(ix1),1,e(ix2),1) + 
     ?                  ddot(norb,x(ix2),1,e(ix1),1)
      sred(ntot,ntot) = ddot(norb,x(ix1),1,s(ix1),1) + 
     ?                  ddot(norb,x(ix2),1,s(ix2),1)
      ered(ntot,ntot) = ddot(norb,x(ix1),1,e(ix1),1) + 
     ?                  ddot(norb,x(ix2),1,e(ix2),1)
      if (odebug(30)) print*,'h(i,i) = ',ered(ntot,ntot)
      if (odebug(43)) then
      call outive(iblbo,maxro,'oosolv: b block addresses before transf')
      call outive(iblso,maxro,'oosolv: s block addresses before transf')
      call outive(ibleo,maxro,'oosolv: e block addresses before transf')
      end if
      call trnsfn (s,2,2,norb,.true.,q)
      call trnsfn (e,2,2,norb,.true.,q)
      call wrt3 (s,norb2,iblso(ntot),ned7)
      call wrt3 (e,norb2,ibleo(ntot),ned7)
      go to 20
14    if (odebug(43)) write (iwr,111) j,iblbo(j)
111   format ('oosolv: trying to read trial vector ',i2,
     +        ' from block ',i12)
      call rdedx (x,norb2,iblbo(j),ned7)
      sred(j,ntot) = ddot(norb,x(ix1),1,s(ix1),1) + 
     +               ddot(norb,x(ix2),1,s(ix2),1)
      ered(j,ntot) = ddot(norb,x(ix1),1,e(ix1),1) + 
     +               ddot(norb,x(ix2),1,e(ix2),1)
      sred(ntot,j) = ddot(norb,x(ix1),1,s(ix2),1) + 
     +               ddot(norb,x(ix2),1,s(ix1),1)
      ered(ntot,j) = ddot(norb,x(ix1),1,e(ix2),1) + 
     +               ddot(norb,x(ix2),1,e(ix1),1)
20    continue

c   solve reduced equation system
      if (odebug(43)) then
      call outive(iblbo,maxro,'oosolv: b block addresses before redeqs')
      call outive(iblso,maxro,'oosolv: s block addresses before redeqs')
      call outive(ibleo,maxro,'oosolv: e block addresses before redeqs')
      end if

      call dcopy(ntot,fred1(      1),1,fred2(     1),1)
      call dcopy(ntot,fred1(maxro+1),1,fred2(ntot+1),1)
      if (odebug(91)) then
         print*,'OOSOLV: iteration step',nitero
         call outvec (fred2,2*ntot,'right hand side')
      end if
      call redeqs (maxro,ntot,frqncy,sred,ered,bred,ared,
     +             fred2,iwork,ierr)
c      if (odebug(91)) call outvec (fred2,2*ntot,'reduced solution:')
      if (ierr .ne. 0) go to 400

c   compute residue vector

      if (odebug(43)) then
      call outive(iblbo,maxro,'oosolv: b block addresses after redeqs')
      call outive(iblso,maxro,'oosolv: s block addresses after redeqs')
      call outive(ibleo,maxro,'oosolv: e block addresses after redeqs')
      end if

      call vclr(s,1,norb2)
      call vclr(e,1,norb2)
      do 140 k=1,ntot
      if (odebug(43)) write (iwr,112) k,iblso(k)
112   format ('oosolv: trying to vector sx no. ',i2,' from block ',i12)
      call rdedx (x,norb2,iblso(k),ned7)
c...for vector computers: no recurrence in this loop
      do 120 i1=1,norb
      i2 = i1 + norb
      s(i1) = s(i1) + x(i1) * fred2(k) - x(i2) * fred2(ntot+k)
      s(i2) = s(i2) + x(i2) * fred2(k) - x(i1) * fred2(ntot+k)
120   continue
      if (odebug(43)) write (iwr,113) k,ibleo(k)
113   format ('oosolv: trying to vector ex no. ',i2,' from block ',i12)
      call rdedx (x,norb2,ibleo(k),ned7)
c...for vector computers: no recurrence in this loop
      do 130 i1=1,norb
      i2 = i1 + norb
      e(i1) = e(i1) + x(i1) * fred2(k) + x(i2) * fred2(ntot+k)
      e(i2) = e(i2) + x(i2) * fred2(k) + x(i1) * fred2(ntot+k)
130   continue
140   continue

c   r = (omega * s - e) - f

      do 50 i=1,norb2
50    e(i) = frqncy * s(i) - e(i) - f(i)

c   residue vector now in e

      if (odebug(47)) then
         write (iwr,61) nitero
61       format ('OPTORB: residue vector no.',i3,':')
         call outvec (e,norb2,' ')
      end if

c   compute euclidean norm of residue vector

      rnorm = dnrm2(norb2,e,1)
      if (ntot.eq.1 .and. nitero.gt.1) then
         write (iwr,60) nitero,rnorm
60       format (6x,i3,8x,f14.9,8x,'RESTART')
      else
         write (iwr,65) nitero,rnorm
65       format (6x,i3,8x,f14.9)
      end if

c   threshold ?

      if (rnorm .le. thresh) go to 230
      if (ntot.ge.maxro .or. ntot.ge.norb) go to 900
      ntot = ntot + 1

c   update residue vector (updated vector in e)

      ix1 = 1
      ix2 = norb + 1
      call lrupda (1,frqncy,e(ix1),e(ix2),x(io1),x(io2),e(ix1),e(ix2))
      if (odebug(48)) call outvec (e,norb2,'updated residue vector:')

c   gram-schmidt orthogonalisation

      do 70 k=1,ntot-1
      if (odebug(43)) write (iwr,114) k,iblbo(k)
114   format ('oosolv (gso): trying to trial vector no. ',i2,
     +        ' from block ',i12)
      call rdedx (x,norb2,iblbo(k),ned7)
      t = - ddot(norb,e(ix1),1,x(ix1),1) - 
     +      ddot(norb,e(ix2),1,x(ix2),1)
      call daxpy(norb,t,x(ix1),1,e(ix1),1)
      call daxpy(norb,t,x(ix2),1,e(ix2),1)
      t = - ddot(norb,e(ix1),1,x(ix2),1) - 
     +      ddot(norb,e(ix2),1,x(ix1),1)
      call daxpy(norb,t,x(ix2),1,e(ix1),1)
      call daxpy(norb,t,x(ix1),1,e(ix2),1)
70    continue
      if (odebug(90)) 
     +     call outvec (e,norb2,'orthogonalized residue vector:')
      t = dnrm2(norb2,e,1)
      if (odebug(90)) print*,'norm: ',t
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,80) nitero
80       format('unable to finish OPTORB after microiteration step ',i2/
     ?          'since new trial vector is linear dependent.')
         ierr = 1
         go to 1000
      else ! perform symmetric orthonormalisation
         if (odebug(49)) then
            call outvec (e(ix1),norb,
     +                   'first  part of vector to be "symorted":')
            call outvec (e(ix2),norb,
     +                   'second part of vector to be "symorted":')
         end if
         call symort (norb,e(ix1),e(ix2),sytol,oldep)
         if (oldep) then
           write (iwr,90) nitero
90         format ('unable to finish OPTORB after microiteration step ',
     ?           i2/'because b and b# are (nearly) linearly dependent.')
           ierr = 1
           go to 1000
         end if
      end if

c   write new trial vector onto external file

      call dcopy(norb2,e,1,x,1)
      call wrt3 (x,norb2,iblbo(ntot),ned7)
      if (odebug(90)) 
     +        print*,'LRSOLV: output trial vector',ntot,' on blocks',
     ?                iblbo(ntot),iblbo(ntot)+lensec(norb2)-1
      go to 10

c   convergence reached

230   icvgd = 1
900   call vclr(x,1,norb2)
      do 110 k = 1 , ntot
         call rdedx (e,norb2,iblbo(k),ned7)
         x1kj = fred2(k)
         x2kj = fred2(k+ntot)
c...for vector computers: no recurrence in this loop
         do 100 i = 1 , norb
            b1ik = e(i)
            b2ik = e(norb+i)
            x(i     ) = x(i     ) + b1ik * x1kj + b2ik * x2kj
            x(norb+i) = x(norb+i) + b2ik * x1kj + b1ik * x2kj
100      continue
110   continue
c      if (odebug(93)) call trnsfn (x,2,1,norb,.true.,q)
      if (icvgd .eq. 0) go to 1000
      if (nitero .eq. 1) write (iwr,240) cpulft(1)-time
      if (nitero .gt. 1) write (iwr,250) nitero,cpulft(1)-time
240   format (/'convergence reached after one microiteration in ',
     ?        f9.4,' seconds.'/)
250   format (/'convergence reached after',i3,
     ?         ' microiterations in',f9.4,' seconds.'/)
      go to 1000

c   error message

400   write (iwr,410) nitero,ierr
410   format (/'microiteration step ',i2,
     ?         ': condition number too small: ierr = ',i3/)
1000  continue

      return
      end
      subroutine optorb (frqncy,thresh,civec,d,x,s,e,f,q,iq,
     ?                   msplit,ntot0,xred,ible)
c
c---------------------------------------------------------------------
c   Control routine for optimal orbital trial vector algorithm
c   (see Joergensen, Jensen, Olsen, J. Chem. Phys. 89, 3654 (1988))
c   (Has not been extensively tested !)
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      integer iblbo, iblso, ibleo
      common /optima/ iblbo(50),iblso(50),ibleo(50)
c
c
      integer io1, io2, ic1, ic2
      common /lradrs/ io1,io2,ic1,ic2
c
      dimension civec(*),d(*),x(*),s(*),e(*),f(*),q(*),iq(*)
      dimension msplit(maxr),xred(maxr2),ible(maxr)
c
      ibase = icorr(0)
c
c   construct vector [v(co),w(co)]
c
      call vclr(f,1,norb2)
      do 140 k = 1 , ntot0
         if (odebug(1)) print*,'OPTORB: msplit(',k,') = ',msplit(k)
         if (msplit(k) .ne. 2) go to 140
         x1kj = xred(k      )
         x2kj = xred(k+ntot0)
         call rdedx (e,npar2,ible(k),ned7)
         if (odebug(1)) call outvec (e,npar2,'vector Ex')
c...for vector computers: no recurrence in this loop
         do 130 i1 = 1 , norb
         f(i1     ) = f(i1     ) + e(i1     ) * x1kj + e(i1+npar) * x2kj
         f(i1+norb) = f(i1+norb) + e(i1+npar) * x1kj + e(i1     ) * x2kj
130      continue
140   continue

c   workspace

      maxro2 = 2 * maxro
      maxros = maxro * maxro
      maxroq = maxro2 * maxro2
      isred  = icorr(maxros)
      iered  = icorr(maxros)
      ibred  = icorr(maxro)
      iared  = icorr(maxroq)
      ifred1 = icorr(maxro2)
      ifred2 = icorr(maxro2)
      iwork  = icori(maxro2)
      iblo7s = iblo7

      do 307 j=1,maxro
      iblbo(j) = iblo7
307   iblo7 = iblo7 + lensec(norb2)

      do 308 j=1,maxro
      iblso(j) = iblo7
308   iblo7 = iblo7 + lensec(norb2)

      do 309 j=1,maxro
      ibleo(j) = iblo7
309   iblo7 = iblo7 + lensec(norb2)
 
      if (odebug(43)) then
         call outive (iblbo,maxro,'optorb: b block addresses')
         call outive (iblso,maxro,'optorb: s block addresses')
         call outive (ibleo,maxro,'optorb: e block addresses')
      end if

c   trial solution

      ix1 = 1
      ix2 = norb + 1
c      call lrupda (1,frqncy,q(if),q(if+norb),s(io1),s(io2),x(ix1),x(ix2))
      call lrupda (1,frqncy,f(ix1),f(ix2),s(io1),s(io2),x(ix1),x(ix2))
c      call symort (norb,x(ix1),x(ix2),sytol,oldep)  ! not necessary
      call wrt3 (x,norb2,iblbo(1),ned7)

      nitero = 0
      icvgd = 0
500   call oosolv (civec,d,x,s,e,frqncy,thresh,q(isred),q(iered),
     ?             q(ibred),q(iared),f,q(ifred1),q(ifred2),iq(iwork),
     ?             q)
      if (ierr .ne. 0) call caserr ('ierr .ne. 0 in OPTORB')
      if (icvgd .eq. 1) go to 520
      if (nitero .gt. maxito)
     ?  call caserr 
     +    ('no convergence of optimal orbital trial vector algorithm')
c      call symort (norb,x(ix1),x(ix2),sytol,oldep)   ! not necessary
c      if (oldep) call caserr 
c    +    ('symort problems in restart run of oosolv')
      call wrt3 (x,norb2,iblbo(1),ned7)
      go to 500
520   continue

c   solution in x --- copy to e

c      call dcopy(norb,x(ix1),1,e(io1),1)
c      call dcopy(norb,x(ix2),1,e(io2),1)
c      call vclr(e(ic1),1,nci1)
c      call vclr(e(ic2),1,nci1)

      call dcopy(norb2,x,1,f,1)

      iblo7 = iblo7s
      call corlsr (ibase)

      return
      end
      subroutine orbadr (isymm,ioff,n,intoff,npair,nact,itypea,
     +                   mult,nirr)
c---------------------------------------------------------------------
c   Constructs the addresses of active orbital pairs
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      implicit integer (i-n)
      dimension ioff(*),intoff(8),npair(8),itypea(nact),mult(8,8)
      n = 0
      do 30 isymvx = 1 , nirr
         isymtu = mult(isymm,isymvx)
         intoff(isymvx) = n
         do 17 isymv = 1 , nirr
            isymx = mult(isymv,isymvx)
            do 16 iv = 1 , nact
               if (itypea(iv) .ne. isymv) go to 16
               do 15 ix = 1 , nact
                  if (itypea(ix) .ne. isymx) go to 15
                  ioff((iv-1)*nact+ix) = n
                  n = n + npair(isymtu)
15             continue
16          continue
17       continue
30    continue
      return
      end
      subroutine outden (mode,isymm,ioff,gam1,gam2)

c---------------------------------------------------------------------
c   Utility routine for output of the density matrices
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension ioff(1),gam1(1),gam2(1)
      if (mode.eq.2) goto 1000
      do 165 isymd=1,nirr
      write(iwr,162)isymd
162   format(/'one particle density matrix for symmetry ',i1,':'/)
      do 164 i=1,nact
      if(itypea(i).ne.isymd) goto 164
      k=0
      do 163 j=1,nact
      if(itypea(j).ne.isymd) goto 163
      k=k+1
      gin(k)=gam1(ipos(ilifa(i)+j))
163   continue
      write(iwr,1635)(gin(kk),kk=1,k)
1635  format(20f14.9)
164   continue
165   continue
1000  continue
      do 299 iv=1,nact
      do 299 ix=1,nact
      call vclr(gin,1,nact2)
      isymvx=mult(itypea(iv),itypea(ix))
      isymtu=mult(isymm,isymvx)
      do 295 it=1,nact
      do 295 iu=1,nact
      if(mult(itypea(it),itypea(iu)).ne.isymtu) goto 295
      gin(ilifa(it)+iu)=gam2(ipos(ilifa(it)+iu)+ioff(ilifa(iv)+ix))
295   continue
      write(iwr,296)iv,ix
296   format(/'two particle density matrix P(..,..,',i2,',',i2,'):'/)
      do 297 i=1,nact
297   write(iwr,298)(gin(ii),ii=ilifa(i)+1,i*nact)
298   format(20f14.9)
299   continue
      return
      end
      subroutine outint (q,n,iblk,ned)
c---------------------------------------------------------------------
c   Writes two-electron integrals to disk.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z), logical (o)
c
      logical odebug
      common /debug/ odebug(100)
c
      dimension q(*)
      nblk = n / 511
      if(odebug(11))
     +     print*,'outint called with n = ',n,'  --> nblk = ',nblk
      if (nblk .eq. 0) return
      nout = nblk * 511
      call wrt3 (q,nout,iblk,ned)
      if(odebug(11))print*,'outint: writing on block ',iblk
      n = n - nout
      call fmove (q(nout+1),q(1),n)
      iblk = iblk + nblk
      return
      end
      subroutine pack22 (nword,iii,i,j)
c---------------------------------------------------------------------
c   Packs four integer*2 numbers into one word.
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
      integer*2 iii(4,*),i(4,*),j(4,*)
      nw = (nword-1)/2 + 1
      kk = 1
      do 10 k=1,nw
      iii(1,k) = i(1,kk)
      iii(2,k) = j(1,kk)
      iii(3,k) = i(1,kk+1)
      iii(4,k) = j(1,kk+1)
10    kk = kk + 2
      return
      end
      subroutine residu (iw,inorms,res,reso,resc,neig,niter,mode)

c---------------------------------------------------------------------
c  MODE = 1: Write norms of residue vectors of each iteration to fortran
c            file INORMS
c  MODE = 2: Read norms of residue vectors from INORMS and write them to
c            standard output in the correct order, representing
c            converged roots as blanks
c  (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z)
      character*130 line
      dimension res(*),reso(*),resc(*),buf(60)
      character*1 digit(10)
      data digit/'0','1','2','3','4','5','6','7','8','9'/

c   for MCLR: 4 roots per line, format f9.5

      data nroots,iform1,iform2 / 4,9,5 /

      m = 1
      n = nroots
      kk = 1
      go to (10,30), mode

c   write mode

10    nn = min(n,neig)
      write (inorms) niter,(reso(i),resc(i),res(i),i=m,nn)
      if (neig .le. n) then
         return
      end if
      m = m + nroots
      n = n + nroots
      go to 10

c   read mode

30    if (nroots .gt. 20) return
      nn = min(n,neig)
      write (iw,40) (i,i=m,nn)
c...next format should be changed to ('iteration',2x,<nroots>(i14,13x))
c...if the value of nroots is altered
40    format ('iteration',2x,4(i14,13x))
      write (iw,41) ('-------------------------',i=1,nn-m+1)
c...in the next two formats, 4(...) should be changed to <nroots>(...)
c...if the value of nroots is altered
41    format ('number',5x,4(a25,2x))
      write (iw,42) ('orb      conf    total',i=1,nn-m+1)
42    format (13x,4(a22,5x))
      rewind inorms
      do 140 i = 1 , niter
         do 50 j = 1 , kk-1
            read (inorms)
50       continue
         read (inorms) idum,(buf(k),k=1,3*(nn-m+1))
         do 110 k = 1 , 3*(nn-m+1)
            temp = buf(k) + 0.5d0 * (10.0d0**(-iform2))
            if (temp .lt. 0.0d0) go to 90
            if (temp .ge. 10.0d0**(iform1-iform2-2)) go to 70
            t = temp
            ipos = k * iform1 - iform2
            line(ipos:ipos) = '.'
            do 60 ii = iform1-iform2-2 , -iform2 , -1
               power = 10.0d0**ii
               u = t / power
               iu = u
               ipos = k * iform1 - iform2 - 1 - ii + min(1,max(0,-ii))
               line(ipos:ipos) = digit(iu+1)
               if (iu.eq.0 .and. ii.ge.1) line(ipos:ipos) = ' '
               t = t - dble(iu) * power
60          continue
            if ((k/3)*3 .eq. k) then
               if (buf(k-2) .ge. buf(k-1)) then
                  ipos = (k-2) * iform1 + 1
               else
                  ipos = (k-1) * iform1 + 1
               end if
               line(ipos:ipos) = '*'
            end if
            go to 110
70          do 80 ipos = (k-1)*iform1+1 , k*iform1
               line(ipos:ipos) = '*'
80          continue
            go to 110
90          do 100 ipos=(k-1)*iform1+1,k*iform1
               line(ipos:ipos) = ' '
100         continue
110      continue
         write (iw,120) i,(line(ipos:ipos),ipos=1,3*(nn-m+1)*iform1)
120      format (i5,4x,130a1)
         do 130 j = kk+1 , (neig-1)/nroots+1
            read (inorms)
130      continue
140   continue
      if (neig .le. n) return
      write (iw,150)
150   format (/)
      m = m + nroots
      n = n + nroots
      kk = kk + 1
      go to 30

      end
      subroutine rmatrx(q)

c---------------------------------------------------------------------
c   Constructs matrix R defined by
c   R(tv)  =  sum_x  D(tx) * FI(xv)  +  sum_uxy  P(utxy) * (xy|vu)
c   (needed in routine oopart)
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      dimension q(*)
      dimension mk1(8)
      ibase = icorr(0)
      time = cpulft(1)
c
c   prepare space
c
      ir = icorr(n0int1)
c
c   read in density matrices
c
      idens = icorr(n0int)
      call rdedx (q(idens),n0int,ibdp,ned7)
      if (odebug(24)) 
     ?   call outvec (q(idens+n0int2),n0int1,
     ?                'one particle density matrix from rmatrx:')
      ip = idens
      id = idens + n0int2
c
c   scratch space
c
      i1 = icorr(max(2*nact2,nfock))
c
c   read in inactive fock matrix
c
      call rdedx (q(i1),nfock,iblofi,ned7)
c
c   FI contributions
c
      ifock = i1
      idens = id
      irmat = ir
      do 50 isyma=1,nirr
         ifi = ifock + mcore(isyma)*nsymm(isyma) + mcore(isyma)
         ma = mact(isyma)
         call mxmaa (q(idens),1,ma,
     ?               q(ifi),1,nsymm(isyma),
     ?               q(irmat),1,ma,
     ?               ma,ma,ma)
         ifock = ifock + nsymm(isyma)**2
         idens = idens + ma*ma
         irmat = irmat + ma*ma
50    continue
      if (odebug(24)) then
         print*,'r matrix, first part'
         irmat = ir
         do i=1,nirr
            print*,'symmetry',i
            call outsqr (q(irmat),mact(i),mact(i),mact(i),' ')
            irmat = irmat + mact(i)**2
         end do
      end if

c   read in all-active integrals

      intpos = 0
      intfil = ned0
      iblf = iblkaa
      intmod = 0
      do 1500 ix=1,nact
         do 1000 iy=1,ix
            j1 = i1
            isymxy = mult(itypea(ix),itypea(iy))
            if (ix .ne. iy) then
               call daxpy(npair(isymxy)
     +         ,1.0d0,q(ip+ioff0(ilifa(iy)+ix)),1 
     +         ,q(ip+ioff0(ilifa(ix)+iy)),1)
            end if

c   coulomb

            do 200 isyma=1,nirr
               isymb = mult(isyma,isymxy)
               if (isymb .gt. isyma) go to 200
               m = mact(isyma) * mact(isymb)
               if (m .eq. 0) go to 200
               call intin (q(j1),m)
               mk1(isyma) = j1
               j1 = j1 + m
200         continue

c   exchange (dummy !)

            if (isymxy .ne. 1) go to 350
            do 300 isyma=1,nirr
               if (mact(isyma) .ne. 0) then
                  call intin (q(j1),mact(isyma)**2)
               end if
300         continue
350         irmat = ir
            do 400 isyma=1,nirr
               isymb = mult(isyma,isymxy)
               if (isyma .ge. isymb) then
                  call mxmb (q(ip+ioff0(ilifa(ix)+iy)
     ?                           +iofsd(isymxy,isymb)),1,mact(isyma),
     ?                       q(mk1(isyma)),mact(isyma),1,
     ?                       q(irmat),1,mact(isyma),
     ?                       mact(isyma),mact(isymb),mact(isyma))
               else
                  call mxmb (q(ip+ioff0(ilifa(ix)+iy)
     ?                           +iofsd(isymxy,isymb)),1,mact(isyma),
     ?                       q(mk1(isymb)),1,mact(isymb),
     ?                       q(irmat),1,mact(isyma),
     ?                       mact(isyma),mact(isymb),mact(isyma))
               end if
               irmat = irmat + mact(isyma)**2
400         continue
1000     continue
1500  continue
      iblkr = iblo7
      call wrt3 (q(ir),n0int1,iblkr,ned7)
      if (odebug(90)) print*,'output r matrix on blocks',iblkr,
     ?                       ' -',iblkr+lensec(n0int1)-1
      iblo7 = iblo7 + lensec(n0int1)
      if (odebug(24)) then
         print*,'r matrix'
         irmat = ir
         do i=1,nirr
            print*,'symmetry',i
            call outsqr (q(irmat),mact(i),mact(i),mact(i),' ')
            irmat = irmat + mact(i)**2
         end do
      end if
      write (iwr,2000) cpulft(1)-time
2000  format (/'construction of r matrix: ',f10.4,' seconds')
      call corlsr (ibase)

      return
      end
      subroutine rplace (icg,isymrs,iw,nw,n,iad,inter)

c---------------------------------------------------------------------
c   Generates one-electron replacements e(pq)|string>.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8  (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension icg(n),iw(maxrpl,3),iad(1),inter(na,nact)
      dimension ia(32)
      do 10 i=1,n
10    ia(i) = icg(i)
      ia1 = ia(1)
      ia(n+1) = 0
      iz = 1
      nw = 0
      do 40 i=1,n
      isymj = mult(isymrs,itypea(ia1))
      kk = 1
      iz2 = iz
      jj = 1
      do 20 j=1,n
20    jj = jj + inter(j,ia(j))
      do 30 j1=1,nact
      ii = jj - inter(kk,ia(kk)) + inter(kk,j1)
      if (j1.eq.ia(kk+1)) then
         kk = kk + 1
         jj = ii
         iz2 = -iz2
      else if (itypea(j1).eq.isymj) then
         nw = nw + 1
         iw(nw,1) = iad(ii)
         iw(nw,2) = iz2
         iw(nw,3) = ilifa(ia1) + j1
      endif
30    continue
      iz = -iz
      ia1 = ia(i+1)
      ia(i+1) = ia(1)
40    ia(1) = ia1
      return
      end
      function seclen (nword)
c---------------------------------------------------------------------
c   Determines how many blocks a section comprises if an array
c   of nword real*8 numbers is to be stored.
c   Identical to GAMESS function lensec except that nword=0 is
c   treated correctly.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-w), integer (i-n), logical (o)
      integer seclen
      if (nword .eq. 0) then
         seclen = 0
      else
         seclen = (nword-1)/511 + 1
      end if
      return
      end
      subroutine sorbit (yo,zo,so,to,d)

c---------------------------------------------------------------------
c   Constructs the orbital parts of the vector Sx
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension yo(*),zo(*),so(*),to(*),d(*)

      do 10 i=1,norbti+norbai
      so(i) =  2.0d0 * yo(i)
10    to(i) = -2.0d0 * zo(i)

      do 20 isymt=1,nirr
      isymi = mult(isymt,isymr)
      i1 = iorbti(isymt) + 1
      i2 = iofsd(1,isymt) + 1
      if (mcore(isymi).ne.0 .and. mact(isymt).ne.0) then
       call mxmbn (yo(i1),1,mcore(isymi), d(i2),1,mact(isymt),
     +             so(i1),1,mcore(isymi), 
     +             mcore(isymi),mact(isymt),mact(isymt))
       call mxmb  (zo(i1),1,mcore(isymi), d(i2),1,mact(isymt),
     +             to(i1),1,mcore(isymi), 
     +             mcore(isymi),mact(isymt),mact(isymt))
      endif
20    continue

      do 30 isyma=1,nirr
      isymt = mult(isyma,isymr)
      i1 = iorbat(isyma) + 1
      i2 = iofsd(1,isymt) + 1
      if( mact(isymt).ne.0 .and. mvirt(isyma).ne.0 ) then
       call mxmb  (d(i2),1,mact(isymt), yo(i1),1,mact(isymt),
     +             so(i1),1,mact(isymt), 
     +             mact(isymt),mact(isymt),mvirt(isyma))
       call mxmbn (d(i2),1,mact(isymt), zo(i1),1,mact(isymt),
     +             to(i1),1,mact(isymt), 
     +             mact(isymt),mact(isymt),mvirt(isyma))
      endif
30    continue

      return
      end
      subroutine mclr_space
c---------------------------------------------------------------------
c   Computes scratch space requirements.
c   Should be obsolete by now.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      data icush/10000/
      istore = icorrm() - icush
      nevmax = imaxf(8,nev,1)
      norbmx = 0
      do 5 i=1,nirr
      k = 0
      do 4 isyma=1,nirr
      isymb = mult(i,isyma)
4     k = k + mcore(isyma)*(mact(isymb)+mvirt(isymb)) 
     +      + mact(isyma)*mvirt(isymb)
5     norbmx = max(norbmx,k)
      ncimax = imaxf(8,lengci,1)
      npamx2 = 2 * (norbmx + ncimax)
      ibegin = (nevmax-1)/intrel + 1         !  iconv = icori(neig)
     ?         + 4 * ((maxr-1)/intrel + 1)   !  iblb,ibls,ible,imspli
      nintmx = 0
      do 20 i=1,nirr
      n = 0
      do 10 isymtu=1,nirr
10    n = n + npair(isymtu) * npair(mult(i,isymtu))
20    nintmx = max(n,nintmx)
      mclrsp = ibegin 
     ?         + maxr2                     !  ietare = icorr(maxr2)
     ?         + nevmax                    !  idenom = icorr(neig)
     ?         + (nevmax-1)/intrel + 1     !  ibleig = icori(neig)
     ?         + npamx2                    !  ix     = icorr(npar2)
     ?         + npamx2                    !  isx    = icorr(npar2)
     ?         + npamx2                    !  ihx    = icorr(npar2)
     ?         + maxrs                     !  isred1 = icorr(maxrs)
     ?         + maxrs                     !  ihred1 = icorr(maxrs)
     ?         + maxr                      !  ibred  = icorr(maxr)
     ?         + maxr2s                    !  isred2 = icorr(maxr2s)
     ?         + maxr2s                    !  ihred2 = icorr(maxr2s)
     ?         + maxr2                     !  ietaim = icorr(maxr2)
     ?         + maxr2                     !  isigma = icorr(maxr2)
     ?         + maxr2s                    !  ixred  = icorr(maxr2s)
     ?         + maxit*neig                !  ires   = icorr(maxit*neig)
      intsp1 = istore - 7 * ((nba-1)/intrel + 1) - nact2 - 2 * nbatri
      intsp2 = istore - mclrsp - 2 * nact2 - 2 * nintmx 
     +                         - 2 * nact4 - nbasq
     ?         - ((max(ncore,nact)-1)/intrel + 1)
     ?         - ((max(ncore,nact,nvirt)-1)/intrel + 1)
     ?         - ((nvirt-1)/intrel + 1)
      intspc = min(intsp1,intsp2)
      komsp1 = (istore - mclrsp) / 7
      komsp2 = (istore - intspc - ibegin - norbmx
     ?                 - 6 * ((nba-1)/intrel + 1)
     ?                 - nact2 - nact4 - 2 * nbatri) / 7
      k0h    = min(komsp1,komsp2)
      k0     = k0h + k0h
      k02    = k0 + k0
      k025   = k02 + k0h
      kspace = istore - nact2 - n0int1 - nci0 - mclrsp
      if (odebug(3)) write (iwr,100) mclrsp,intspc,k0,kspace
100   format (/'space statistics: '/
     ?         'mclrsp (space for mclr algorithm) :',i8,' words'/
     ?         'intspc (space for integrals)      :',i8,' words'/
     ?         'k0     (space for commutators)    :',i8,' words'/
     ?         'kspace (new space for commutators):',i8,' words')
      return
      end
      subroutine startv (q,n,k,index)
c
c---------------------------------------------------------------------
c   Searches the k lowest entries of the diagonal of A
c   in order to construct starting vectors.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension q(*),index(*)
c
      osort = .false.
      call quick (q,n,k,index,osort)
      ilen = k
      if (odebug(15)) call outvec (q,ilen,'q in STARTV')
c
      do 30 i=1,ilen
      m = i
      qmin = q(i)
      do 20 j=i+1,ilen
      if (q(j) .lt. qmin) then
         m = j
         qmin = q(j)
      end if
20    continue
      t = q(i)
      q(i) = q(m)
      q(m) = t
      l = index(i)
      index(i) = index(m)
30    index(m) = l
c
      if (odebug(15)) call outvec (q,ilen,'q in STARTV after ordering')
c
      return
      end
      subroutine strong (mt,m,n,icg,itype,jmt,iu,mult)

c---------------------------------------------------------------------
c   Generates a new alpha or beta string.
c   Identical to GAMESS routine string except that the case
c   of 0 electrons is handled correctly.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8  (a-h,o-z)
      dimension icg(n),itype(m),mult(8,8)
c
c...carsten0
      if (n.eq.0) then
         if (mt .gt. 0) go to 60
         mt = 1
         jmt = 1
         return
      end if
c...carsten1
      if (mt.gt.0) goto 20
      do 10 i=1,n
10    icg(i)=i
      iu=n
      goto 40
c
20    icg(iu) = icg(iu) + 1
      if (icg(iu).le.m+iu-n) goto 30
      iu = iu-1
      if (iu.eq.0) goto 60
      goto 20
30    if (iu.eq.n) goto 40
      iu = iu+1
      icg(iu) = icg(iu-1)+1
      goto 30
c
40    mt = mt+1
      jmt = 1
      do 50 j=1,n
50    jmt = mult(jmt,itype(icg(j)))
      return
60    mt = 0
      return
      end
      subroutine trnsfn (x,nvec,mode,nstrid,onat,q)

c---------------------------------------------------------------------
c   Transforms NVEC orbital coefficient vectors, supplied in 
c
c     X((I-1)*NSTRID+1),...,X((I-1)*NSTRID+NORB), I=1,...,NVEC,
c
c   to the canonical orbital basis or back
c
c   MODE = 1 : form  U*X*U(T)
c   MODE = 2 : form  U(T)*X*U
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit real*8 (a-h,p-z), logical (o)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      dimension x(*),q(*),ia(3)

      if (.not.ocanon .or. norb.eq.0) return

      is = icorr(nvec*norb)
      iu = icorr(ntrans)
      call rdedx (q(iu),ntrans,ibcan,ned7)
      if (onat) go to 5
      done = 1.0d0
      k = iu
      do 1 i=1,nirr
1     k = k + mcore(i)**2
      do 2 i=1,nirr
      if (mact(i) .eq. 0) go to 2
      call vclr(q(k),1,mact(i)**2)
      call dcopy(mact(i),done,0,q(k),mact(i)+1)
2     continue
5     ia(1) = 1
      ia(3) = 1

c   first part of transformation

c   t-i

      do 10 isymac=1,nirr
      isymco = mult(isymac,isymr)
      ia(2) = mcore(isymco)
      do 10 i=1,nvec
10    call mxmaa (q(iu+iofcor(isymco)),ia(mode),ia(mode+1),
     ?            x((i-1)*nstrid+iorbti(isymac)+1),1,mcore(isymco),
     ?            q(is+(i-1)*norb+iorbti(isymac)),1,mcore(isymco),
     ?            mcore(isymco),mcore(isymco),mact(isymac))

c   a-i

      do 20 isymvi=1,nirr
      isymco = mult(isymvi,isymr)
      ia(2) = mcore(isymco)
      do 20 i=1,nvec
20    call mxmaa (q(iu+iofcor(isymco)),ia(mode),ia(mode+1),
     ?            x((i-1)*nstrid+iorbai(isymvi)+1),1,mcore(isymco),
     ?            q(is+(i-1)*norb+iorbai(isymvi)),1,mcore(isymco),
     ?            mcore(isymco),mcore(isymco),mvirt(isymvi))

c   a-t

      do 30 isymvi=1,nirr
      isymac = mult(isymvi,isymr)
      ia(2) = mact(isymac)
      do 30 i=1,nvec
30    call mxmaa (q(iu+iofact(isymac)),ia(mode),ia(mode+1),
     ?            x((i-1)*nstrid+iorbat(isymvi)+1),1,mact(isymac),
     ?            q(is+(i-1)*norb+iorbat(isymvi)),1,mact(isymac),
     ?            mact(isymac),mact(isymac),mvirt(isymvi))

c   second part of transformation

c   t-i

      do 40 isymac=1,nirr
      isymco = mult(isymac,isymr)
      ia(2) = mact(isymac)
      do 40 i=1,nvec
40    call mxmaa (q(is+(i-1)*norb+iorbti(isymac)),1,mcore(isymco),
     ?            q(iu+iofact(isymac)),ia(3-mode),ia(4-mode),
     ?            x((i-1)*nstrid+iorbti(isymac)+1),1,mcore(isymco),
     ?            mcore(isymco),mact(isymac),mact(isymac))

c   a-i

      do 50 isymvi=1,nirr
      isymco = mult(isymvi,isymr)
      ia(2) = mvirt(isymvi)
      do 50 i=1,nvec
50    call mxmaa (q(is+(i-1)*norb+iorbai(isymvi)),1,mcore(isymco),
     ?            q(iu+iofvir(isymvi)),ia(3-mode),ia(4-mode),
     ?            x((i-1)*nstrid+iorbai(isymvi)+1),1,mcore(isymco),
     ?            mcore(isymco),mvirt(isymvi),mvirt(isymvi))

c   a-t

      do 60 isymvi=1,nirr
      isymac = mult(isymvi,isymr)
      ia(2) = mvirt(isymvi)
      do 60 i=1,nvec
60    call mxmaa (q(is+(i-1)*norb+iorbat(isymvi)),1,mact(isymac),
     ?            q(iu+iofvir(isymvi)),ia(3-mode),ia(4-mode),
     ?            x((i-1)*nstrid+iorbat(isymvi)+1),1,mact(isymac),
     ?            mact(isymac),mvirt(isymvi),mvirt(isymvi))

      call corlsr (is)
      return
      end
      subroutine upak4t (nword,iii,ij)
c---------------------------------------------------------------------
c   Unpacks the indices of the entries of A(oo),B(oo)
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
      integer*2 iii(*),ij(4,*)
      do 10 loop=1,nword
         ij(1,loop) = iii(loop)
10    continue
      return
      end
      subroutine vectrs (diag,civec,iblock,msplit,symbol,ifile,q,iq)
c
c---------------------------------------------------------------------
c   Chooses the starting vectors for the iterative MCLR procedure.
c   Usually, the unit vectors corresponding to the lowest values
c   of the vector diag(A)/diag(Sigma) are chosen as initial guesses.
c   Configuration trial vectors are spin-adapted.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
      character *3 symbol(14)
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
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
      logical odebug
      common /debug/ odebug(100)
c
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      dimension diag(*),civec(*),iblock(*),msplit(*),q(*),iq(*)
      dimension icga(31),icgb(31),iocc(62),icount(8),lcore(8)
      dimension lvirt(8),mk1(8)
      data m51,small/51,1.0d-7/
      ind(i,j) = ((i-1)*(i-2))/2 + j
c
      ibase = icorr(0)
      ono = onospl .or. nsplit.gt.1
      ix = icorr(npar2)
      iy = ix
      iz = ix + npar
      iy1 = iy - 1
      iz1 = iz - 1
      ixx = icorr(npar2)
      iofa = icori(nstraa)
      iofb = icori(nstrbb)
      inter = icori(intlen)
      call setsto(intlen,0,iq(inter))
      call ciadrs (2,isym,iq(iofa),iq(iofb),iq(inter),na)
c
c   next line for common/multic/ so that SPINAD & ISTRAD work
c
      nlen = nci1
      if (orpa) go to 500
c
c   search for the NEIG smallest elements of diagonal
c
      do 15 ir = nst , nba
         if (ir .le. nprim) then
            is1 = ncore
         else
            is1 = nprim
         end if
         xoccr = 0.0d0
         if (ir .le. nprim) xoccr = ddiag(ir-ncore)
         do 10 is = 1 , is1
            if (mult(itype(ir),itype(is)) .eq. isymr) then
               k = ipo(ind(ir,is))
               xoccs = 2.0d0
               if (is .gt. ncore) xoccs = ddiag(is-ncore)
               xoccd = xoccs - xoccr
               if (dabs(xoccd) .ge. 1.0d-7) then
                  diag(k) = diag(k) / xoccd
               else
                  diag(k) = 1.0d0
               end if
               if (odebug(14)) 
     +         print*,'xocc(',is,') - xocc(',ir,') = ',xoccd
            end if
10       continue
15    continue
      if (odebug(14)) call outvec (diag,npar,'modified diagonal:')
      call vabs (diag,1,diag,1,npar)
      iquick = icori(npar)
      if (odebug(60)) then
         do 20 i = 1 , npar
            iq(iquick-1+i) = i
20       continue
      else
         ngsafe = min(5*neig,npar)
         call startv (diag,npar,ngsafe,iq(iquick))
      end if

c   choose starting vectors

      write (iwr,30)
30    format (/65('*')/'*',20x,
     +         'Choice of starting vectors',17x,'*'/65('*')/)
      istart = 0
      do 120 ii = 1 , ngsafe
         index = iq(iquick-1+ii)
         call vclr(q(ix),1,npar2)
         if (index .le. norb) then
            istart = istart + 1
            do 45 ir = nst , nba
               if (ir .le. nprim) then
                  is1 = ncore
               else
                  is1 = nprim
               end if
               do 40 is = 1 , is1
                  if (mult(itype(ir),itype(is)) .eq. isymr) then
                     k = ipo(ind(ir,is))
                     if (k. eq. index) go to 50
                  end if
40             continue
45          continue
50          xoccr = 0.0d0
            if (ir .le. nprim) xoccr = ddiag(ir-ncore)
            xoccs = 2.0d0
            if (is .gt. ncore) xoccs = ddiag(is-ncore)
            xoccd = xoccs - xoccr
            write (iwr,60) istart,index,
     ?                     isymmo(is),symbol(nirs+itype(is)),
     ?                     isymmo(ir),symbol(nirs+itype(ir)),
     ?                     diag(ii),diag(ii)*xoccd
60          format ('trial vector ',i2,
     +              ' corresponds to orbital rotation ',i4/
     ?              '     orbitals:                      ',i2,a3,
     +              ' -->',i3,a3/
     ?              '     value of diagonal element:      ',f14.9,
     ?              '  (diag(A)/diag(Sigma))'/
     ?              '                                     ',f14.9,
     ?              '  (diag(A))'/)
            msplit(istart) = 1
            if (ono) msplit(istart) = 3
            q(ix-1+index) = 1.0d0
         else
            ix2 = ix
            if (ono) ix2 = ix + norb
            ixx2 = ixx
            if (ono) ixx2 = ixx + norb
            iconf = index - norb
            if (odebug(60)) then
               q(ix2-1+iconf) = 1.0d0
               go to 80
            end if
            call spinad (q(ix2),iconf,iq(iofa),iq(iofb),iq(inter))
            if (odebug(13)) 
     +      call outvec (q(ix2),nci1,'spin adapted vector')
            if (isym .eq. isss) call orth (nci0,q(ix2),civec)
            if (odebug(13)) 
     +      call outvec (q(ix2),nci1,'orthogonalized vector')
            do 70 jstart = 1 , istart
               if (msplit(jstart) .eq. 1) go to 70
               call rdedx (q(ixx),2*leng(msplit(jstart)),
     +                            iblock(jstart),ned7)
               call orth (nci1,q(ix2),q(ixx2))
70          continue
            t = dnrm2(nci1,q(ix2),1)
            if (t .le. small) then
               write (iwr,71) iconf
71             format ('discard configuration trial vector ',i7)
               go to 120
            end if
            call dscal(nci1,1.0d0/t,q(ix2),1)
80          istart = istart + 1
            msplit(istart) = 2
            if (ono) msplit(istart) = 3
            call setsto(2*nact,0,iocc)
            call detocc (iconf,icga,icgb)
            do 90 iaa = 1 , na
               iocc(2*icga(iaa)-1) = 1
90          continue
            do 100 ibb = 1 , nb
               iocc(2*icgb(ibb)) = 1
100         continue
            write (iwr,110) istart,iconf,(iocc(k),k=1,2*nact)
110         format ('trial vector ',i2,
     +              ' corresponds to configuration ',i6/
     ?              '     occupation pattern:             ',31(2i1,1x))
            write (iwr,115) diag(ii)
115         format ('     value of diagonal element:      ',f14.9/)
         end if
         n2 = 2*leng(msplit(istart))
         call wrt3 (q(ix),n2,iblock(istart),ned7)
         if (odebug(90)) print*,'VECTRS: output starting vector',istart,
     ?                  ' on blocks',iblock(istart),' -',iblock(istart)
     ?                                                +lensec(n2)-1
         if (istart .ge. neig) go to 200
120   continue
      write (iwr,130)
130   format (
     +'unable to construct a sufficient number of starting vectors ...'/
     ?'either there is something wrong or you have to increase the',
     ?' value of NGSAFE in the program !')
      call caserr ('error in mclr program')
200   continue
      go to 1000
500   continue
c...choose rpa eigenvectors as starting vectors
      ncrpa = ncore + na
      nvrpa = nba - ncrpa
      if (odebug(88)) print*,'ncrpa = ',ncrpa,'      nvrpa = ',nvrpa
      fac = 1.0d0 / dsqrt(2.0d0)
c...read in symmetries of scf mos
      nav = lenwrd()
      jsymao = icori(65+2*maxorb+nav)
      call secget (isect(490),m51,iblks)
      call readi(iq(jsymao),mach(13)*nav,iblks,ned3)
      jsymmo = jsymao + 64 + maxorb
      if (odebug(88)) 
     +    call outive(iq(jsymmo+1),nbasis,'scf mo symmetries')
c...how many core and virtual orbitals in the individual symmetries ?
      do 510 i=1,nirr
      lcore(i) = 0
510   lvirt(i) = 0
      do 520 ii=1,ncrpa
      k = iq(jsymmo+ii)
520   lcore(k) = lcore(k)+1
      do 530 ia=ncrpa+1,nba
      k = iq(jsymmo+ia)
530   lvirt(k) = lvirt(k)+1
      iu = icorr(nfock)
      call rdedx (q(iu),nfock,iblku,ned7)
      k = 0
      do 540 isyma=1,nirr
      mk1(isyma) = k
540   k = k+nsymm(isyma)**2
      read (ifile) isymrp,neigrp
      if (isymrp .ne. isym .or. neigrp .ne. neig) 
     ?   call caserr ('file with rpa vectors has incompatible format')
      ilook1 = icori(ncrpa*nvrpa)
      ipo1 = ilook1 - 1
c...rpa lookup table
      call setsto(8,0,icount)
      do 560 ii=1,ncrpa
      do 560 ia=1,nvrpa
      k = mult(iq(jsymmo+ii),iq(jsymmo+ia+ncrpa))
      l = icount(k)+1
      iq(ipo1+(ii-1)*nvrpa+ia) = l
560   icount(k) = l
c...read vectors from fortran file
      isize = icount(isym)
      if(odebug(88))print*,'isize = ',isize
      iyvec = icorr(isize)
      izvec = icorr(isize)
      iyvec1 = iyvec - 1
      izvec1 = izvec - 1
      do 800 j=1,neig
      iblock(j) = iblo7
      call vclr(q(ix),1,npar2)
      read (ifile) (q(iyvec1+k),k=1,isize)
      read (ifile) (q(izvec1+k),k=1,isize)
      if (odebug(88)) then
         print*,'rpa eigenvector, symmetry',isym,' no. ',j
         call outvec(q(iyvec),isize,' ')
         call outvec(q(izvec),isize,' ')
      end if
c...transform rpa eigenvector
      do 630 isyma=1,nirr
      isymb = mult(isyma,isymr)
      if (lvirt(isyma)*lcore(isymb).eq.0) go to 630
      iya = icorr(nsymm(isyma)*nsymm(isymb))
      iyb = icorr(nsymm(isyma)*nsymm(isymb))
      iza = icorr(nsymm(isyma)*nsymm(isymb))
      izb = icorr(nsymm(isyma)*nsymm(isymb))
      call vclr(q(iya),1,nsymm(isyma)*nsymm(isymb))
      call vclr(q(iza),1,nsymm(isyma)*nsymm(isymb))
      icol = 0
      do 625 ii=1,ncrpa
      if (iq(jsymmo+ii).ne.isymb) go to 625
      icol = icol+1
      irow = 0
      do 620 ia=ncrpa+1,nba
      if (iq(jsymmo+ia).ne.isyma) go to 620
      irow = irow+1
      iold = iq(ipo1+(ii-1)*nvrpa+ia-ncrpa)
      inew = (icol-1)*nsymm(isyma)+lcore(isyma)+irow-1
      q(iya+inew) = q(iyvec1+iold)
      q(iza+inew) = q(izvec1+iold)
620   continue
625   continue
      call mxms (q(iya),1,nsymm(isyma),q(iu+mk1(isymb)),1,nsymm(isymb),
     +           q(iyb),1,nsymm(isyma),
     +           nsymm(isyma),nsymm(isymb),nsymm(isymb))
      call mxms (q(iu+mk1(isyma)),nsymm(isyma),1,q(iyb),1,nsymm(isyma),
     +           q(iya),1,nsymm(isyma),
     +           nsymm(isyma),nsymm(isyma),nsymm(isymb))
      call mxms (q(iza),1,nsymm(isyma),q(iu+mk1(isymb)),1,nsymm(isymb),
     +           q(izb),1,nsymm(isyma),
     +           nsymm(isyma),nsymm(isymb),nsymm(isymb))
      call mxms (q(iu+mk1(isyma)),nsymm(isyma),1,q(izb),1,nsymm(isyma),
     +           q(iza),1,nsymm(isyma),
     +           nsymm(isyma),nsymm(isyma),nsymm(isymb))
      if (odebug(88)) then
         print*,'VECTRS: symmetry',isym,', y eigenvector ',j,
     +          ', part',isyma
         call outsqr (q(iya),nsymm(isyma),nsymm(isyma),nsymm(isymb),' ')
         print*,'VECTRS: symmetry',isym,', z eigenvector ',j,
     +          ', part',isyma
         call outsqr (q(iza),nsymm(isyma),nsymm(isyma),nsymm(isymb),' ')
      end if
      icol = 0
      do 665 is=1,nprim
      if (itype(is) .ne. isymb) go to 665
      icol = icol+1
      irow = mcore(isyma)
      do 660 ir=nst,nba
      if (itype(ir) .ne. isyma) go to 660
      irow = irow+1
      if (ir .le. is) go to 660
      k = (icol-1)*nsymm(isyma)+irow-1
      if (is .le. ncore .or. ir .gt. nprim) then
         q(iy1+ipo(ind(ir,is))) = q(iya+k)
         q(iz1+ipo(ind(ir,is))) = q(iza+k)
      else if (is .le. ncrpa .and. ir .gt. ncrpa) then
         it = is-ncore
         do 640 iorb=1,na-1
         icga(iorb) = iorb
640      if (iorb .ge. it) icga(iorb) = iorb+1
         icga(na) = ir-ncore
         do 650 iorb=1,na
650      icgb(iorb) = iorb
         ipar = 1
         mta = istrad (ipar,iq(inter        ),icga,na)
         mtb = istrad (ipar,iq(inter+na*nact),icgb,nb)
         if (ipar .ne. 1) print*,'WARNING: ipar .ne. 1 in VECTRS !'
         ff = dble(ipar) * fac
         q(iy1+norb+iq(iofa-1+mta)+iq(iofb-1+mtb)) = ff * q(iya+k)
         q(iy1+norb+iq(iofa-1+mtb)+iq(iofb-1+mta)) = ff * q(iya+k)
         q(iz1+norb+iq(iofa-1+mta)+iq(iofb-1+mtb)) = ff * q(iza+k)
         q(iz1+norb+iq(iofa-1+mtb)+iq(iofb-1+mta)) = ff * q(iza+k)
      end if
660   continue
665   continue
      call corlsr (iya)
630   continue
      if (isym .eq. isss) then
         call orth (nci0,q(iy+norb),civec)
         call orth (nci0,q(iz+norb),civec)
      end if
      do 790 k=1,j-1
      call rdedx (q(ixx),npar2,iblock(k),ned7)
790   call orth (npar2,q(ix),q(ixx))
      t = dnrm2(npar2,q(ix),1)
      if (1.0d0+t .eq. 1.0d0) call caserr ('singularity ...')
      call dscal(npar2,1.0d0/t,q(ix),1)
      call wrt3 (q(ix),npar2,iblock(j),ned7)
800   iblo7 = iblo7 + lensec(npar2)
      do 810 j=1,neig
810   msplit(j) = 3
1000  call corlsr (ibase)
      return
      end
      subroutine wmatrx (q)
c
c---------------------------------------------------------------------
c   Constructs the matrix   W(pi,tu) = 2 (pi|tu) - (pt|ui)   (p active/virtual)
c   needed for the modified one-electron coefficients fv(tu) (--> MODIFY)
c   and the configuration-orbital part (--> COPART)
c   (pi) and (tu) are both of symmetry ISYMR 
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit real*8 (a-h,p-z), logical (o)
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
      logical odebug
      common /debug/ odebug(100)
c
      character*80 dmpfil,rstfil
      common/blkin /acorev(512),potnuc,dx,dy,dz,sout(510)
      common/blkin1/evalue(maxorb),rocc(maxorb),etot,nbasis,newbas,
     +              ncol,ivalue,ioccup,ispa
      common/cntrol/iseca,isecb,maxr,maxr2,maxrs,maxr2s,maxit,ispc8,
     +              eps,epsana,sytol,
     +              iblkc,ned3,ned4,ned6,ned7,itech,ocanon,nbegin,
     +              onospl,ono4,
     +              epstab,oanal,oreduc,otrue,oexca,iwast(2),multip,
     +              ointin,
     +              ozint,oipos,orevrs,ostore,cshift,oinvrs,isecao,
     +              nitst,iseccn,orpa,odirc,olsen,osafe,onewa,onoana,
     +              ococo,nsplit,ootva,ootvsp,
     +              tfac,maxro,maxito,onoad,onoad4,dsmall,oshift
      common/craypk/mmmm(65),ibfcod(maxorb),itype(maxorb)
      common/detcic/na,nb,isym,icf(32),npairs(8),nstra(8),nstrb(8),
     +              intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb,
     +              maxrpl,maxpar
      common/gjs   /nirri,multi(8,8),isymai(maxorb),isymmo(maxorb)
      common/infoa /nat,ich,mul,num
      common/insort/ned0,ibl124,ibl356,iblkgv,iblkaa,iblkwm(8)
      common/intbuf/intpos,intfil,intmod
      common/iter  /niter,neig,tim(4),leng(3),noff(3),shtim(10),
     +              ierr,icvgd,
     +              cotim(10),nitero
      common/junk  /gin(511),ddiag(31)
      common/junk2 /iii(340),jjj(340),kkk(340),lll(340)
      common/lineax/nacnac,nconac,nconco,nvinco,nvinac,nbasq,nbatri,
     +              nact2,nact3,nact4
      common/mcblok/iblo7,ibdp,iblkr,ibzint,iblofi,iblofa,ibdiag,
     +              ibloab,ibcan,
     +              iblkpp,iblkw,iblku,ibco
      common/mccore/intrel,lword,ltop,lmax,lmin
      common/mcints/list(7),lpair(7),nbox(7),nbox0(7),nbox1(7),nbl(7)
      common/multic/waste(40),iwaste(24),nba,ncore,nact,nlen,norbr,
     +              nprim,
     +              nvirt,nst,nfreez,ifreez(8),nprimp,nirr
      common/posit /ipo((maxorb*(maxorb-1))/2),ipos(961),
     +              ioff0(961),ioff1(961)
      common/scrach/intspc,k0h,k0,k02,k025,kspace,lp,iblo4
      common/syminf/mult(8,8),npoint(maxorb),kstart(8),nfin(8),nsymm(8),
     +              nprm(8),mcore(8),mact(8),mvirt(8),ic1d,ne,
     +              itypea(31),
     +              ilifa(maxorb),jwaste(558)
      common/symlok/ic2e(496),ic3e(496)
      common/symsca/isss,isymr,maxrp3,nci1,n0int1,n0int2,n0int,
     +              n1int1,n1int2,
     +              n1int,intlen,norb,norb2,nci0,npar,npar2,
     +              norbti,norbai,norbat,ntrans,nirs,nfock,nparam,nci,
     +              lenrs
      common/symvek/nev(8),lengci(8),npair(8),
     +              norbs(8),intof0(8),intof1(8),
     +              iorbti(8),iorbai(8),iorbat(8),iofcor(8),
     +              iofact(8),iofvir(8),
     +              ivioff(8),iofsd(8,8),kstra(8),kstrb(8)
      common/table /dmpfil,rstfil,rlines
      common/values/excmax,core,e0,val,fvisum,excave,const,trycor,
     +              denomy(50)
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      dimension q(1)
      dimension mk1(8),mk2(8),mk3(8),mk4(8),mk5(8),mk6(8)
      time = cpulft(1)
      ibase = icorr(0)
      intpos = 0
      intfil = ned0
      iblf = iblkwm(isymr)
      intmod = 0
      i1 = icorr(2*nbasq)
      irow = norbti+norbai
      icol = npair(isymr)
      isize = irow*icol
      if (isize .eq. 0) go to 9999
      write (iwr,10) irow,icol
10    format (/'size of W matrix:',i4,' x ',i4/)
      i2 = icorr(isize)
      do 1000 it=1,nact
      do 1000 iu=1,it
      if (mult(itypea(it),itypea(iu)) .ne. isymr) go to 1000
      j1 = i1
c
c   coulomb integrals (tu|..)
c
      do 200 isyma=1,nirr
      isymb = mult(isyma,isymr)
      if (isymb .gt. isyma) go to 200
      mca = mcore(isyma)
      maa = mact(isyma)
      mva = mvirt(isyma)
      mcb = mcore(isymb)
      mab = mact(isymb)
      mvb = mvirt(isymb)
      m = (maa+mva)*mcb+(mab+mvb)*mca
      if (m .eq. 0) go to 200
      call intin (q(j1),m)
      if (odebug(80)) then
         write (iwr,11) it+ncore,iu+ncore,isyma,isymb
11       format(/'wmatrx: active-active pair ',2i3,
     ?           ' --- coulomb integrals for symmetries ',2i2,':'/)
         call outsqr (q(j1),mcb,mcb,maa,'(tu|vi) - active/core:')
         call outsqr (q(j1+maa*mcb),mca,mca,mab,
     +                                  '(tu|iv) - core/active:')
         call outsqr (q(j1+maa*mcb+mca*mab),mcb,mcb,mva,
     ?                '(tu|ai) - virtual/core:')
         call outsqr (q(j1+maa*mcb+mva*mcb+mca*mab),mca,mca,mvb,
     ?                '(tu|ia) - core/virtual:')
      end if
      mk1(isyma) = j1
      if (isymb .ne. isyma) mk1(isymb) = mk1(isyma)+maa*mcb
      mk2(isyma) = j1+maa*mcb+mca*mab
      if (isymb .ne. isyma) mk2(isymb) = mk2(isyma)+mva*mcb
      j1 = j1+m
200   continue

c   exchange integrals (t.|.u)

      do 300 isyma=1,nirr
      isymb = mult(isyma,isymr)
      mca = mcore(isyma)
      maa = mact(isyma)
      mva = mvirt(isyma)
      mcb = mcore(isymb)
      mab = mact(isymb)
      mvb = mvirt(isymb)
      m = (maa+mva)*mcb+(mab+mvb)*mca
      if (m .eq. 0) go to 300
      call intin (q(j1),m)
      if (odebug(80)) then
         write (iwr,21) it+ncore,iu+ncore,isyma,isymb
21       format (/'wmatrx: active-active pair ',2i3,
     ?            ' --- exchange integrals for symmetries ',2i2,':'/)
         call outsqr (q(j1),mcb,mcb,maa,'(ti|vu) - active/core:')
         call outsqr (q(j1+maa*mcb),mca,mca,mab,
     +                                  '(tv|iu) - core/active:')
         call outsqr (q(j1+maa*mcb+mca*mab),mcb,mcb,mva,
     ?                '(ti|au) - virtual/core:')
         call outsqr (q(j1+maa*mcb+mva*mcb+mca*mab),mca,mca,mvb,
     ?                '(ta|iu) - core/virtual:')
      end if
      mk3(isyma) = j1
      mk5(isymb) = mk3(isyma)+maa*mcb
      mk4(isyma) = mk5(isymb)+mca*mab
      mk6(isymb) = mk4(isyma)+mva*mcb
      j1 = j1+m
300   continue

c   compute columns (tu),(ut)

      iofftu = i2+(ipos(ilifa(it)+iu)-1)*irow
      ioffut = i2+(ipos(ilifa(iu)+it)-1)*irow
      do 400 isyma=1,nirr
      isymb = mult(isyma,isymr)
      maa = mact(isyma)
      mcb = mcore(isymb)
      mti = maa*mcb
      if (mti .eq. 0) go to 370
      itu = iofftu+iorbti(isyma)
      iut = ioffut+iorbti(isyma)
c---------------------------------------------------------------------
c     W(vi,tu) = 2 (vi|tu) - (vt|ui) = 2 (tu|vi) - (tv|iu)           |
c---------------------------------------------------------------------
      call daxmyz (mti,2.0d0,q(mk1(isyma)),q(mk5(isyma)),q(itu))
      if (it .eq. iu) go to 370
c---------------------------------------------------------------------
c     W(vi,ut) = 2 (vi|ut) - (vu|ti) = 2 (tu|vi) - (ti|vu)           |
c---------------------------------------------------------------------
      call daxmyz (mti,2.0d0,q(mk1(isyma)),q(mk3(isyma)),q(iut))
370   mva = mvirt(isyma)
      mai = mva*mcb
      if (mai .eq. 0) go to 400
      itu = iofftu+iorbai(isyma)
      iut = ioffut+iorbai(isyma)
c---------------------------------------------------------------------
c     W(ai,tu) = 2 (ai|tu) - (at|ui) = 2 (tu|ai) - (ta|iu)           |
c---------------------------------------------------------------------
      call daxmyz (mai,2.0d0,q(mk2(isyma)),q(mk6(isyma)),q(itu))
      if (it .eq. iu) go to 400
c---------------------------------------------------------------------
c     W(ai,ut) = 2 (ai|ut) - (au|ti) = 2 (tu|ai) - (ti|au)           |
c---------------------------------------------------------------------
      call daxmyz (mai,2.0d0,q(mk2(isyma)),q(mk4(isyma)),q(iut))
400   continue
1000  continue
      iblkw = iblo7
      call wrt3 (q(i2),isize,iblkw,ned7)
      if (odebug(90)) print*,'output w matrix on blocks',iblkw,
     ?                       ' -',iblkw+lensec(isize)-1
      iblo7 = iblo7+lensec(isize)
9999  call corlsr (ibase)
      if (odebug(80)) call outsqr (q(i2),irow,irow,icol,'w matrix')
      write (iwr,10000) cpulft(1)-time
10000 format ('construction of W matrix ',29('.'),f9.4,' seconds')
      return
      end
      subroutine ver_mclr(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mclr.m,v $
     +     "/
      data revision /"$Revision: 6317 $"/
      data date /"$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
