      subroutine g_rd2el(ptocc,orbsym,buf,bkt,ibkt,ni,nj,nk,nl,
     &                    ibktsp,lnbuf,idu,bkof1,bkof2,tfile,nbf,
     &                    lenbf,maxval)
      implicit real*8  (a-h,o-z)
      integer ptocc(nbf),orbsym(nbf)
      real*8 buf(lnbuf),bkt(ibktsp)
      integer ibkt(ibktsp),ni(lnbuf),nj(lnbuf),nk(lnbuf),nl(lnbuf)
      integer idu(nbf),bkof1(16),bkof2(16),tfile(16)
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
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
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common /cor/ core
      common/blkin/gin(510),numi
      common/craypk/i205(1360)
c     
c... this is the interface to the one and two electron integrals
c... reads them from transformed integral file
c
      call setsto(1360,0,i205)
      do 8 ifile=1,lfile
      iunit=lotape(ifile)
      call search(liblk(ifile),iunit)
      call find(iunit)
      lbl=llblk(ifile)
    6 lbl=lbl+1
      call get(gin(1),mw)
      if(mw.eq.0) go to 8
      if(lbl.ne.0) call find(iunit)
      icnt=0
      int4=1
      call unpack(gin(num2e+1),lab816,i205,numlab)
      do 30 num=1,numi
      j=i205(int4  )
      i=i205(int4+1)
      l=i205(int4+2)
      k=i205(int4+3)
      newi=ptocc(i)
      isnew=orbsym(newi)
      newj=ptocc(j)
      jsnew=orbsym(newj)
      ijsnw=ieor(isnew,jsnew)
      newk=ptocc(k)
      ksnew=orbsym(newk)
      newl=ptocc(l)
      lsnew=orbsym(newl)
      klsnw=ieor(ksnew,lsnew)
      if (ijsnw.eq.klsnw)then
       icnt=icnt+1
       newij=or(min(newi,newj),ishft(max(newi,newj),8))
       ijov=or(idu(min(newi,newj)),ishft(idu(max(newi,newj)),1))
       ni(icnt) = newij
       nj(icnt) = ijov
       nk(icnt) = max(newk,newl)
       nl(icnt) = min(newk,newl)
       buf(icnt)=gin(num)
      endif
      int4=int4+4
   30 continue
      call prcss(buf,bkt,ibkt,ni,nj,nk,nl,icnt,ibktsp,lnbuf,
     &             idu,bkof1,bkof2,tfile,nbf,lenbf,maxval)
      if(lbl.ne.0) go to 6
    8 continue
c
      return
      end
      subroutine g_rd1el(ptocc,orbsym,honeel,nbf)
      implicit real*8  (a-h,o-z)
      integer ptocc(*),orbsym(*),nbf
      real*8 honeel(*)
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
      common /cor/ core
      common/blkin/gin(510),numi
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
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
      integer nact, nna, nnb, isss, nci, ischem, itype
      integer maxit, iacc1, iacc2, ivec1, ivec2, nguess, nroots, nsecor
      integer ivec3, ivec4
      integer ncor
      common /savem/ nact,nna,nnb,isss,nci,ischem,itype(maxorb),
     +              maxit,iacc1,iacc2,ivec1,ivec2,nguess,nroots,nsecor,
     +              ivec3,ivec4,ncor
c
      common/junk/occ(maxorb),potn,corree,ncolo(3),ncore,
     *mapcie(maxorb),map2(maxorb),nactt,mapaie(maxorb),mapaei(maxorb)
     *,iqsec
      external fget
c...
c... restore 1-elec integrals
c...
      call secget(nsecor,1004,jblkk)
      call rdedx(occ,mach(15),jblkk,idaf)
      write(iwr,1100)nsecor,ibl3d,yed(idaf)
1100  format(
     */' transformed 1-electron integrals restored from section',i4,/,
     *' of dumpfile starting at block',i8,' of ',a4)
      core=corree
      nin=0
      jblkk=jblkk+lensec(mach(15))
      call search(jblkk,idaf)
      call fget(gin,kword,idaf)
      do i=1,nbf
      inew=ptocc(i)
      isnew=orbsym(inew)
      do j=1,i
      jnew=ptocc(j)
      jsnew=orbsym(jnew)
      nin=nin+1
      if(nin.gt.kword)then
       call fget(gin,kword,idaf)
       nin=1
      endif
      ij = ((max(inew,jnew)-1)*max(inew,jnew))/2 + min(inew,jnew)
c      print *,'i,j,inew,jnew,isnew,jsnew,nin ',
c     &         i,j,inew,jnew,isnew,jsnew,nin,gin(nin)
      if(isnew.eq.jsnew)then
        honeel(ij)=gin(nin)
      else
        honeel(ij)=0.0d00
      endif
      enddo
      enddo
      return
      end
      subroutine titandrv(q,energy)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/t_adata/a0,a1,a2,a3,a4,a5,a6,a9,ap5
c
      character*10 charwall
      dimension q(*)
      call cpuwal(begin,ebegin)
      write(iwr,1100)
 1100 format(/1x,104('=')/
     *40x,22('*')/
     *40x,'coupled cluster module'/
     *40x,22('*')//)
c
      ibase = igmem_alloc_all(lword)
c
      if(lword.lt.1)call caserr(
     +        'insufficient memory for coupled cluster module')
c
c set t_adata
c
      a0 = 0.0d0
      a1 = 1.0d0
      a2 = 2.0d0
      a3 = 3.0d0
      a4 = 2.0d0
      a5 = 3.0d0
      a6 = 6.0d0
      a9 = 9.0d0
      ap5= 0.5d0
c
      call titandrv2(q(ibase),lword,energy)
c
      at=cpulft(1)
      write(iwr,1300) at ,charwall()
 1300 format(/1x,'end of coupled cluster module at',f10.2,' seconds'
     *,a10,' wall'/
     */1x,104('=')/)
      call clredx
      call timana(29)
      call gmem_free(ibase)
      return
      end
      subroutine titandrv2(q,llword,energy)
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
      dimension q(llword)
      common /mccore/ intrel,lword,ltop,lmax,lmin
      common /ftimes / t(10)
c
      integer nact, nna, nnb, isss, nci, ischem, itype
      integer maxit, iacc1, iacc2, ivec1, ivec2, nguess, nroots, nsecor
      integer ivec3, ivec4
      integer ncor
      common /savem/ nact,nna,nnb,isss,nci,ischem,itype(maxorb),
     +              maxit,iacc1,iacc2,ivec1,ivec2,nguess,nroots,nsecor,
     +              ivec3,ivec4,ncor
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common /cor/ core
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
      integer nsymx,noccx(8),nvirx(8)
      real*8 potnucx
c
      integer nccit,iccty,iccth
      common/tit_inp/nccit,iccty,iccth
c
c
      call gettim(startc,starte)
      lword = llword
      call initcc(q,lword,potnucx,ncorx)
      if (nact.gt.maxorb .or.nna.ne.nnb.or.nna.le.0)then
        write(iwr,*) nact,nna,nnb
        call caserr('invalid parameters in ccsd module')
      endif
      write (iwr,30) nact,nna-ncorx,isss,(itype(i),i=1,nact)
30    format(/
     &' number of orbitals          =',i3/
     &' number of occupied orbtials =',i3/
     &' space symmetry              =',i3//
     &' orbital symmetries           '/3(40i2/1x)/)
      write (iwr,40) lword
40    format(' memory available =',i10,' reals'/)
c
      call icopy(8,0,0,noccx,1)
      call icopy(8,0,0,nvirx,1)
      nsymx=0
      do i=1,nna-ncorx
       is=itype(i)
       noccx(is)=noccx(is)+1
      enddo
      do i=nna-ncorx+1,nact
       is=itype(i)
       nvirx(is)=nvirx(is)+1
      enddo
      is=8
      do while(noccx(is).eq.0.and.nvirx(is).eq.0)
       is=is-1
      enddo
      nsymx=1
      do while(nsymx.lt.is)
       nsymx=nsymx*2
      enddo
      write(iwr,50)(noccx(i),i=1,nsymx)
 50   format(' occupied per symmetry ',8i4)
      write(iwr,60)(nvirx(i),i=1,nsymx)
 60   format(' virtual per symmetry  ',8i4)
c
      call flushn(iwr)
      call tsort(q,q,lword,nsymx,noccx,nvirx,potnucx,itype)
      call flushn(iwr)
      iconvi=iccth
      maxit=nccit+1
      if (iccty.eq.1)iopt=0
      if (iccty.eq.2)iopt=3
      if (iccty.eq.3)iopt=5
      if (iccty.eq.4)iopt=6
      call vccsd(q,q,lword,iconvi,maxit,iopt,energy)
      call flushn(iwr)
      return
      end
      subroutine initcc(q,lword,potnucx,ncorx)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
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
      real*8 potnucx
      parameter (mxorb1=maxorb+1)
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/junk/pop(maxorb),potn,core,ncolo(3),ncore,
     *mapcie(maxorb),map2(maxorb),nactt,mapaie(maxorb),mapaei(maxorb)
     *,iqsec,nacta(maxorb),nactb(maxorb),nactc(5),isecor,
     *evalue(maxorb),eocc(mxorb1),nbas,newb,ncol,ieval,ipop,ispp
     *,nirr,mult(8,8),isymao(maxorb),isymmo(maxorb),nsp
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
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
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      common/junkc/zjob,zdate,ztime,zprog,ztype,zspace(14),ztext(10)
c
      integer nact, nna, nnb, isss, nci, ischem, itype
      integer maxit, iacc1, iacc2, ivec1, ivec2, nguess, nroots, nsecor
      integer ivec3, ivec4
      integer ncor
      common /savem/ nact,nna,nnb,isss,nci,ischem,itype(maxorb),
     +              maxit,iacc1,iacc2,ivec1,ivec2,nguess,nroots,nsecor,
     +              ivec3,ivec4,ncor
c
      dimension q(*)
      data thresh/1.0d-5/
      data m0/0/
      data m51,m29/51,29/
      data wwww/0.0d0/
c
      nav=lenwrd()
      i=inicor(lword)
      call secget(isect(470),1005,iblka)
      call readi(nacta,mach(14)*nav,iblka,idaf)
      call timer(wwww)
      nsecor=isecor
      call secget(isecor,1004,iblka)
      call rdedx(pop,mach(15),iblka,idaf)
      write(iwr,5004) yed(idaf),ibl3d,isecor
 5004 format(/1x,'dumpfile resides on ',a4,' at block ',i6/
     * 1x,'core hamiltonian to be retrieved from section ',i4)
      write(iwr,2)nact,ncore,potn,core
    2 format(1x,'header block information :'/
     *       1x,'=========================='/
     *       1x,'number of active orbitals ',i4/
     *       1x,'number of core   orbitals ',i4/
     *       1x,'nuclear repulsion energy ',e21.12/
     *       1x,'core energy              ',e21.12)
      potnucx=core
      ncorx=ncore
      if(nact.gt.1) go to 10005
      call caserr('invalid number of active orbitals')
80010 call caserr('parameter error in ccsd preprocessor')
10005 do 82003 i=1,nact
      j=mapaie(i)
82003 mapaei(j)=i
      lfile=m6file
      do 60001 i=1,lfile
      lotape(i)=m6tape(i)
      liblk  (i)=m6blk(i)
60001 llblk  (i)=liblk(i)-m6last(i)
      nact=nactt
      call secget(isect(490),m51,iblka)
      call readi(nirr,mach(13)*nav,iblka,idaf)
      call setsto(nact,m0,itype)
      call secget(iqsec,3,iblka)
      call rdchr(zjob,m29,iblka,idaf)
      call reads(evalue,mach(8),idaf)
      write(iwr,89005) iqsec,ztype,ztime,zdate,zjob,ztext,nbas,ncol
89005 format(/1x,'scf mo specifications restored from section ',i4,
     *          ' of dumpfile'/
     *       1x,'header block information :'/
     *       1x,'=========================='/
     *       1x,a7,' vectors created at ',
     *          a8,' on ',a8,' in the job ',a8/
     *       1x,'with the title :',1x,10a8/
     *       1x,'no. of gtos        ',i4/
     *       1x,'no. of scf vectors ',i4)
      call symvec(q,isymao,isymmo,nbas,ncol,iblka)
      do 44 i=1,nact
      itype(i)=isymmo(mapaie(i))
 44   continue
      nocc=0
      nunocc=0
      do 89006 i=1,nact
      jscf=mapaie(i)
      qio=eocc(jscf)
      if(qio.le.thresh) nunocc=nunocc+1
      if(qio.gt.thresh) nocc=nocc+1
89006 continue
      if(nocc+nunocc.ne.nact) go to 80010
      if(nunocc.le.0) go to 80010
      write(iwr,1)
 1    format(/1x,'transformed integral files'/1x,26('*'))
      call filprn(m6file,m6blk,m6last,m6tape)
      return
      end
      subroutine ccsdin
      implicit real*8  (a-h,o-z),integer (i-n)
      character *8 ztext
      character *4 itext,ifd
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
c     **** input routine for ccsd ****
c
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/infoa/nat(2),mul,num,nxx,nelec,naa,nbb
      common/data1/vlist(400),newpro(206),louta(2,maxorb),norbt,
     * norbta(maxorb),norbtp,norbc
      common/scfwfn/cicoef(699),nop(10),nco,nseto,npair
     *,ibms(4),nope
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
      integer nact, nna, nnb, isss, nci, ischem, itype
      integer maxit, iacc1, iacc2, ivec1, ivec2, nguess, nroots, nsecor
      integer ivec3, ivec4
      integer ncor
      common /savem/ nact,nna,nnb,isss,nci,ischem,itype(maxorb),
     +              maxit,iacc1,iacc2,ivec1,ivec2,nguess,nroots,nsecor,
     +              ivec3,ivec4,ncor
c
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
      common/timez/timlim,ti(3),safety(6),isecss
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer nccit,iccty,iccth
      common/tit_inp/nccit,iccty,iccth
      dimension ifd(6)
      data ifd/
     * 'ccit','ccty','ccth','****','****','****'/
c
c
      top=cpulft(1)
      write(iwr,444)top
 444  format(//1x,104('=')//
     *' **** coupled cluster input processor called at',f9.2,' secs')
      if(jump.eq.1) then
c       mspin=mul
        nint=nelec-nope-npair-npair
        nint=nint/2-norbc
        if(nseto.ne.0)then
          do 8002 i=1,nseto
8002      nint=nint+nop(i)
        endif
        nint=nint+npair+npair
        next=norbt-nint
        nact = nint+next
        nna = naa -norbc
        nnb = nbb - norbc
      else
        call inpi(nact)
        call inpi(nna)
        call inpi(nnb)
      endif
      write(iwr,446)nact,nna,nnb
446   format(/' # orbitals       =',i6,
     *//      ' # alpha electrons=',i6,
     *//      ' # beta electrons =',i6)
      if(nact.gt.maxorb)call caserr(
     * 'invalid number of active orbitals')
      if(nnb.ne.nna.or.nnb.lt.1)call caserr(
     * 'invalid number of electrons (CCSD closed shell only)')
c
      nna = nna + norbc
      nnb = nnb + norbc
      isss = 1
      nccit=20
      iccth=10
      jrec = 0
      call inpa(ztext)
      if(ztext.eq.'ccsd(t)' ) then
       iccty = 2
      else if (ztext.eq.'qcisd') then
       iccty = 3
      else if (ztext.eq.'qcisd(t)') then
       iccty = 4
      else
       iccty = 1
      endif
93    call input
c.... see what password is on it
      call inpa4(itext)
      ii=locatc(ifd,3,itext)
      if (ii.gt.0) go to 99
      jrec=jrec-1
      ii = locatc(ydd(101),limit2,itext)
      if(ii)9995,9996,9995
 9996 call caserr(
     *'unrecognised directive or invalid directive sequence')
c.... go to proper place
99    go to (3,4,5),ii
c.... iterations
 3    call inpi(nccit)
      if(nccit.le.0)nccit=20
      go to 93
c.... cc type
 4    call inpi(iccty)
      if(iccty.le.0.or.iccty.ge.5)iccty=1
      go to 93
c.... threshold
 5    call inpi(iccth)
      if(iccth.le.0)iccth=10
      go to 93
c
9995  continue
c
      return
      end
      subroutine ver_ccsd(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/ccsd.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
