      subroutine hfprop3(zscft,core)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      character *8 zscft
      logical opg_root
      logical     some
      logical     out 
      character*8 scf,mp2
      character*8 prop
      character*4 iblk 
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
c
c
      logical lspac2,mp3,pump2,lcanon
      integer ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint
      integer iconvv,np,mn,imolp,iorder,ipropc,nruns,nmol,natre,nch
      integer nmul,nbas,nsh,nelect,iopp,norder,nsys,itoli,icuti
      integer n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf
      integer npstar,npfin,minvec,icflag,mpflag,mpstrm
      integer mpblk,ispare,irest6
      integer irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv
      integer irblok,irunit,irfile,irbl,mcrest,ngpts
      integer nvirtb,nsb,mpfill,mprest
      integer len_cndx41
c     
      common/cndx41/
     + ncoorb,nocc,nocca,noccb,nvirt,nvirta,ipol,iprint,iconvv,np,
     + mn,imolp,iorder,ipropc,nruns,nmol,natre(2),nch(2),
     + nmul(2),nbas(2),nsh(2),nelect(2),iopp(15),norder,nsys,itoli,icuti
     +,n1st,nrec0,intlo0,iresti,iposf1,iposf2,iposm,iochf(100),
     + npstar,npfin,minvec,icflag,mpflag,mpstrm(20),mpblk(20),
     + mp3,pump2,lcanon,lspac2(16),ispare(52),irest6,
     + irest1,irestp,mppas,mcicfl,nps1,nps2,nrefs,ionsv,irblok,irunit,
     + irfile,irbl,mcrest,ngpts,nvirtb,nsb,mpfill,mprest
       parameter(len_cndx41=289)
c      used: paczer(master),restre(util1),revise(util1),utyp21(server)
      logical ociopt, omp2
      common/restrl/ ociopt(2),omp2
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc 
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
      common/junk2/xyzmap(3,100),zznuc(maxat),nuc(maxat)
c
      dimension core(*)
c
      data scf     /'scf     '/
      data mp2     /'mp2     '/
      data prop    /'propty  '/
      data iblk    /'   '/
      data nodip,noqdp,nootp,nopop,nospin,noloc,nodpl,nodhp,nod2hp,
     1     nofpl,nofhp,nof2hp,nosos,noelp,noelf,noelfg,noden,nonmr,
     2     noston
     3    /  0,    0,    0,    0,    1,    1,    0,    1,    1,
     4       1,    1,    1,    1,    0,    0,    0,    0,    1,
     5       1/
c
      out   =.false.
      some  =.true.        
      some  =some.or.out
      some  =some.and.opg_root()
c
c     set up interface values
c
      nbfnd = num
c
      idena = 0
      noscf = .false.
c
      nefc  =0
      iefc  =iblk
      do iat=1,nat
         nuc(iat)=nint(czan(iat))
      enddo
      if (omp2.or.mp3) then
       wfntyp = mp2
      else
       wfntyp = scf  
      endif
c
      if(nodhp.eq.0.or.nod2hp.eq.0) nodpl=0
      if(nofhp.eq.0.or.nof2hp.eq.0) nofpl=0
      if(out) write(iwr,9999) nodip,noqdp,nootp,nopop,nospin,noloc,
     1                        nodpl,nodhp,nod2hp,nofpl,nofhp,nof2hp,
     2                        nosos,noelp,noelf,noelfg,noden,nonmr,
     3                        noston
      if(some) write(iwr,9998) zruntp,zscft
      if(nodip.eq.0.or.noloc.eq.0.or.nodpl.eq.0)write(iwr,12340)
      if(noqdp.eq.0) write(iwr,12341)
      if(nootp.eq.0) write(iwr,12342)
      if(noelp.eq.0) write(iwr,12343)
      if(noelf.eq.0) write(iwr,12344)
      if(noelfg.eq.0) write(iwr,12345)
      if(noden.eq.0) write(iwr,12346)
      if(nopop.eq.0) write(iwr,12347)
c
c
c     ----- dipole moment -----
c
      if(nodip.eq.0.or.noloc.eq.0.or.nodpl.eq.0) then
        write(iwr,12348)
        call dipol2(core)
      endif
c
c     ----- quadrupole moment -----
c
      if(zruntp.ne.prop.and.noqdp.eq.0) then
       write(iwr,12348)
       call qdpole(core)
      endif

c
c     ----- octupole moment -----
c
      if(zruntp.ne.prop.and.nootp.eq.0) then
       write(iwr,12348)
       call octpol(core)
      endif
c
c
c     ----- mulliken population analysis -----
c
      if(nopop.eq.0) then
       write(iwr,12348)
       call mulken2(core)
      endif
c
c
c     ----- electrostatic potential -----
c
      if(zruntp.ne.prop.and.noelp.eq.0) then
       write(iwr,12348)
       call elpmap(core)
      endif
c
c     ----- electric field -----
c
      if(zruntp.ne.prop.and.noelf.eq.0) then
       write(iwr,12348)
       call elfmap(core)
      endif
c
c     ----- electric field gradient -----
c
      if(zruntp.ne.prop.and.noelfg.eq.0) then
       write(iwr,12348)
       call efgmap(core)
      endif
c
c     ----- electron and spin densities -----
c
      if(noden.eq.0) then
       write(iwr,12348)
       call eldmap(core,zscftp)
       call spind2(core)
      endif
c
      return
 9999 format(/' ----- property selection -----'/
     1 ' nodip =',i2,' noqdp =',i2,' nootp =',i2,
     2 ' nopop =',i2,' nospin=',i2,
     3 ' noloc =',i2,' nodpl =',i2/' nodhp =',i2,' nod2hp=',i2,
     4 ' nofpl =',i2,' nofhp =',i2,' nof2hp=',i2,' nosos =',i2,
     5 ' noelp =',i2/' noelf =',i2,' noelfg=',i2,' noden =',i2,
     6 ' nonmr =',i2,' noston=',i2)
 9998 format(/1x,'runtype = ',a8/
     1        1x,'scftype = ',a8//5x,
     1'compute the following with respect to atomic centres'/
     2 5x,
     3'++++++++++++++++++++++++++++++++++++++++++++++++++++'/)
12348 format(/1x,80('=')/)
12340 format(10x,'dipole moment')
12341 format(10x,'second and quadrupole moments')
12342 format(10x,'third and octupole moments')
12343 format(10x,'electrostatic potential')
12344 format(10x,'electric field')
12345 format(10x,'electric field gradient')
12346 format(10x,'electron and spin densities')
12347 format(10x,' .. plus a mulliken population analysis')

      end
      subroutine dipol2(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical opg_root, ono, oroot
      character*8 scf, uhf, mp2, prop
      character*4 keyefc
      logical     some,out, ofield
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      common/three/am(maxat),
     +             dmxe(3),dmq(3),dmqefc(3),dxt(3),dxtefc(3),
     +             axs(3,3),thdum(176)
      dimension x(*)
      character *1 xx
      dimension xx(3)
      dimension cm(3)
      data xx /'x','y','z'/
      data scf,uhf /'scf     ','uhf     '/
      data mp2     /'mp2     '/
      data debye   /2.54176568d+00/
      data zero    /0.0d+00/
      data tol     /1.0d-06/
      data prop    /'propty  '/
      data keyefc  /' efc'/
      data small/1.0d-16/
      data m3/3/
c
      ofield = iefc.eq.keyefc
      out =.false.
      ono = .false.
      some=.true.
      oroot = opg_root()
      some=(some.and.oroot).or.out
      if(zruntp.eq.prop) some=.false.
      if (oroot) write(iwr,9993)
      if (some) write(iwr,9994)
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
c
c     ----- allocate memory block -----
c
      need = l2+l2+l2+l2+l2+l3+l3
      init = igmem_alloc(need)
c
      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of core -----
c
      i10=init
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
c     ----- memory block for da -----
      i50=i40+l2
c     ----- memory block for db -----
      i60=i50+l2
c     ----- memory block for -a- vectors -----
      i70=i60+l3
c     ----- memory block for -b- vectors -----
      last=i70+l3
c
      ncall=1
   10 continue
      call prpams(am,ncall,ncode,some)
      if(ncode.ne.0) go to 100
c
      do 20 i=1,3
   20 cm(i)=zero
      sum=zero
      do 30 iat=1,nat
      sum=sum+am(iat)
      do 30 i=1,3
   30 cm(i)=cm(i)+am(iat)*c(i,iat)
      if( abs(sum).lt.tol) go to 45
      do 40 i=1,3
   40 cm(i)=cm(i)/sum
   45 continue
      if(some) write(iwr,9997) cm(1),cm(2),cm(3)
c
c     use centre of mass, or origin nominated under gauge?
c
      if (orgin) then
       if (some) write(iwr,1234) (ggx(i),i=1,3)
       call dcopy(3,ggx(1),1,cm(1),1)
      else
       if (some) write(iwr,1235)
      endif
c
c     ----- calculate dipole moment integrals -----
c
      call dipint2(x(i10),x(i20),x(i30),cm,l1)
c
      ibl7x = iposun(num8)
      call wrt3(x(i10),l2,ibl7x,num8)
      ibl7y = iposun(num8)
      call wrt3(x(i20),l2,ibl7y,num8)
      ibl7z = iposun(num8)
      call wrt3(x(i30),l2,ibl7z,num8)
      ibl7la = iposun(num8)
c
      if(out) then
         call writel(x(i10),l1)
         call writel(x(i20),l1)
         call writel(x(i30),l1)
      endif
c
      if(zruntp.eq.prop) go to 110
c
c     ----- electronic contribution to dipole moment -----
c
c     ----- density matrix for x(i40)=da,x(i50)=db -----
c     ----- calculate density matrix -----

      call prpdens(x(i40),x(i50),x(i60),x(i70),zscftp,ono,oroot)
      if (ono) return
c
      if(out) then
         call writel(x(i40),l1)
         call writel(x(i50),l1)
      endif
c
      dmxe(1)=-tracep(x(i40),x(i10),l1)
      dmxe(2)=-tracep(x(i40),x(i20),l1)
      dmxe(3)=-tracep(x(i40),x(i30),l1)
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 50
       if(iuno.ne.0.and.iunosp.ne.0) go to 50
        dmxe(1)=-tracep(x(i50),x(i10),l1)+dmxe(1)
        dmxe(2)=-tracep(x(i50),x(i20),l1)+dmxe(2)
        dmxe(3)=-tracep(x(i50),x(i30),l1)+dmxe(3)
   50 continue
c
      do  i=1,3
         dmq(i)=zero
         dmqefc(i)=zero
      enddo
c
c     ----- nuclear contribution
c
      do 60 i=1,nat
      dmq(1)=dmq(1)+czan(i)*(c(1,i)-cm(1))
      dmq(2)=dmq(2)+czan(i)*(c(2,i)-cm(2))
      dmq(3)=dmq(3)+czan(i)*(c(3,i)-cm(3))
   60 continue
c
c     ----- -efc- contribution -----
c
      if(ofield) then
        do i=1,nefc
         dmqefc(1)=dmqefc(1)+efcz(i)*(efcc(1,i)-cm(1))
         dmqefc(2)=dmqefc(2)+efcz(i)*(efcc(2,i)-cm(2))
         dmqefc(3)=dmqefc(3)+efcz(i)*(efcc(3,i)-cm(3))
       enddo
      endif
c
      do i = 1 , 3
       dxt(i) = dmq(i) + dmxe(i)
       dxtefc(i) = dxt(i) + dmqefc(i)
      enddo
      dipol = ddot(m3,dxt,1,dxt,1)
      if (dipol.gt.small) dipol = dsqrt(dipol)
      dipefc = ddot(m3,dmqefc,1,dmqefc,1)
      if (dipefc.gt.small) dipefc = dsqrt(dipefc)
      diptot = ddot(m3,dxtefc,1,dxtefc,1)
      if (diptot.gt.small) diptot = dsqrt(diptot)
c
      if(some) then
       write (iwr,6000)
       if(ofield) then
         write(iwr,6005)
         do i = 1 , 3
           write (iwr,6010)xx(i),dmq(i),dmqefc(i),dmxe(i),dxtefc(i)
         enddo
         write (iwr,6020) dipefc, diptot
       else
         write (iwr,6030)
         do i = 1 , 3
           write (iwr,6010)xx(i),dmq(i) ,dmxe(i) ,dxt(i)
         enddo
         write (iwr,6040) diptot
       endif
      endif
      do i = 1, 3
       dmq(i) = dmq(i) * debye
       dmqefc(i) = dmqefc(i) * debye
       dmxe(i) = dmxe(i) * debye
       dxtefc(i) = dxtefc(i) * debye
       dxt(i) = dxt(i) * debye
      enddo
      dipol=dipol*debye
      dipefc=dipefc*debye
      diptot=diptot*debye
      if(some) then
        write(iwr,6015)
       if(ofield) then
         write (iwr,6005)
         do i = 1 , 3
           write (iwr,6010)xx(i),dmq(i),dmqefc(i),dmxe(i),dxtefc(i)
         enddo
         write (iwr,6060) dipefc, diptot
       else
         write (iwr,6030)
         do i = 1 , 3
           write (iwr,6010)xx(i),dmq(i) ,dmxe(i) ,dxt(i)
         enddo
         write (iwr,6070) diptot
       endif
      endif
c
      ncall=ncall+1
      ncode=1 
c     go to 10
c
  100 continue
c
c     ----- project dipole components onto inertial axes -----
c
      if(some.and.ncode.eq.0) 
     &   call inrtia(am,cm,axs,1)
      if(some.and.ncode.eq.0)
     &   write (iwr,9999)
      if(some.and.ncode.eq.0)
     &   call dipaxs(axs,'a','b','c',dxt(1),dxt(2),dxt(3))
c
c
  110 continue
c
c     ----- release memory block -----
c
      call gmem_free(i10)
c
      return
c
 9993 format(/10x,13('-')
     1       /10x,'dipole moment',
     2       /10x,13('-')/)
 9994 format(' 1 a.u. = 2.541766 debyes ')
 9999 format(/2x,'dipole moments along inertial axes')
 9998 format(' not enough core in -dipole- . stop ')
 9997 format(/2x,' centre of mass',
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 6000 format (/2x,'dipole moment (a.u.) (w.r.t. above origin)')
 6005 format(/11x,'nuclear',11x,'-efc-',6x,'electronic',11x,'total'/)
 6010 format (1x,a1,4f16.7)
 6020 format (/1x,'-efc- dipole moment = ',f16.7,' (a.u.)'/
     +         1x,'total dipole moment = ',f16.7,' (a.u.)')
 6030 format (/11x,'nuclear',6x,'electronic',11x,'total'/)
 6040 format (/1x,'total dipole moment = ',f16.7,' (a.u.)')
 6015 format (/2x,'dipole moment (debye) (w.r.t. to above origin)')
 6060 format (/1x,'-efc- dipole moment = ',f16.7,' (debye)'/
     +         1x,'total dipole moment = ',f16.7,' (debye)')
 6070 format (/1x,'total dipole moment = ',f16.7,' (debye)')
 1234 format (/1x,'origin for property evaluation:', 
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 1235 format (/1x,'origin for property evaluation is centre of mass')
      end
      subroutine prpams(bmass,ncall,ncode,out)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical out
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
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      dimension bmass(*)
      data zero /0.0d+00/
c
      if(ncall .gt. mass_numvec())then
        ncode = 1
        return
      endif
      do i = 1,nat
        bmass(i) = amass_get(ncall,i)
c
        if (bmass(i).lt.zero) then
           write (iwr,6010)
           ncode = 1
           return
        end if
      enddo
c
      if(out) then

      write (iwr,6020)
       do  iat = 1 , nat
          write (iwr,6030) iat , zaname(iat) , bmass(iat)
       enddo
      end if
c
      ncode = 0
      return
c
 6010 format (' *** zero isotope mass encountered  ***')
 6020 format(/10x,21('=')/10x,'atomic weights (a.u.)'
     1       /10x,21('=')/)
 6030 format(i5,5x,a8,2x,f15.5)
      end
      subroutine dipint2(xs,ys,zs,cm,l1)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical out
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
      common/xyzdip/xint,yint,zint,xintx,yinty,zintz,t,xc,yc,zc,
     1              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension xs(*),ys(*),zs(*)
      common/junk3/sx(225),sy(225),sz(225),
     + dij(225),
     + xin(5,5), yin(5,5), zin(5,5),
     + xxin(5,5),yyin(5,5),zzin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension cm(3)
      data zero  /0.0d+00/
      data one   /1.0e+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      norm=normf.ne.1.or.normp.ne.1
c
      xc=cm(1)
      yc=cm(2)
      zc=cm(3)
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
      sx(ij)=zero
      sy(ij)=zero
      sz(ij)=zero
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
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
c     ----- dipole moment integrals -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljt
      nj=j
      do 370 i=1,lit
      ni=i
      call dipxyz2
       xin(i,j)=xint*t1
       yin(i,j)=yint*t1
       zin(i,j)=zint*t1
      xxin(i,j)=xintx*t1
      yyin(i,j)=yinty*t1
      zzin(i,j)=zintz*t1
  370 continue
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 390 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      ij=ij+1
      sx(ij)=sx(ij)+dij(ij)*(xxin(ix,jx)* yin(iy,jy)* zin(iz,jz))
      sy(ij)=sy(ij)+dij(ij)*( xin(ix,jx)*yyin(iy,jy)* zin(iz,jz))
      sz(ij)=sz(ij)+dij(ij)*( xin(ix,jx)* yin(iy,jy)*zzin(iz,jz))
  390 continue
c
 6000 continue
 7000 continue
c
c     ----- set up dipole moment matrices -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=iky(loci+i)+(locj+j)
      xs(nn)=sx(ij)
      ys(nn)=sy(ij)
      zs(nn)=sz(ij)
 7500 continue
c
 8000 continue
 9000 continue
c
      if(out) then
         write(iwr,9999)
         call writel(xs,l1)
         write(iwr,9998)
         call writel(ys,l1)
         write(iwr,9997)
         call writel(zs,l1)
      endif
c
      return
 9999 format(/10x,25('-')/10x,'x-dipole moment integrals'/
     1         10x,25('-'))
 9998 format(/10x,25('-')/10x,'y-dipole moment integrals'/
     1         10x,25('-'))
 9997 format(/10x,25('-')/10x,'z-dipole moment integrals'/
     1         10x,25('-'))
      end
      subroutine dipxyz2
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzdip/xint,yint,zint,xintx,yinty,zintz,t,xc,yc,zc,
     1 x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min /1,2,4, 7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      xintx=zero
      yinty=zero
      zintz=zero
      npts=(ni+nj-2+1)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
      xintx=xintx+px*(ptx-xc)
      yinty=yinty+py*(pty-yc)
      zintz=zintz+pz*(ptz-zc)
   16 continue
      return
      end
      subroutine prpdens(densa,densb,veca,vecb,zscft,ono,oroot)
c
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *8 anos, bnos
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/bufb/eiga(maxorb),frocca(maxorb),eigb(maxorb),
     +            froccb(maxorb)
      common/junkc/zcom1(29),zcomm(19),zbuff(10),zatom(maxat)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
c
      common/multic/radius(40),irad(57+mcfzc),itype(maxorb),isecn
      character*8 wfntyp, mp2
      common/prop3c/wfntyp, iefc
      dimension densa(*),densb(*),veca(*),vecb(*),ytype(2)
      data mp2     /'mp2     '/
      data anos,bnos/' nos',' snos'/
      data ytype/'-a-','-b-'/
      data itvec/3/
      data m990/990/
c
      oexist = .false.
      odens = .false.
      ono = .false.
c
      if(zruntp.ne.'ci'.and.zscftp.eq.'casscf') then
       odens = .true.
       idena = isecda
         call secloc(isecda,oexist,iblko)
         if (.not.oexist) then
            if(oroot) write (iwr,6035)
            ono = .true.
         end if
      endif
      if (idena.gt.0 ) odens = .true.
      newbas = num
c
c     decide on section for retrieving appropriate orbital set
c**** need to resolve this when in analysis mode
c
      minva = mouta
      minvb = moutb
c
c     check for CI calculation - then load spinfree and spin NOS
c     if direct-ci
c
      if (zruntp.eq.'ci') then
      if (omrdci.or.opark) then
        if(nwroot.gt.0) then
c        note this must be modified for multiple roots
         minva = nosec(1)
c        no spin NOs for table-driven MRDCI 
         minvb = 0
        else
         if(oroot) write(iwr,6025)
         ono = .true.
         return
        endif
      else if (ofulci) then
        if(oroot) write(iwr,6020)
        ono = .true.
        return
      else if (occsd.or.oqcisd) then
        if(oroot) write(iwr,6000)
        ono = .true.
        return
      else
c     direct-ci
        minva = ispacg
       if (ispacg .le. 0 ) then
        ono = .true.
        if(oroot) write(iwr,6025)
        return
       endif
       if (isping.gt.0 ) then
        minvb = isping
       endif
      endif
      endif
c
c     check for MP2 or MP3 calculation
c
      if (wfntyp.eq.mp2) then
       if(iuno.gt.0) then
        minva = iuno
       else
        if(oroot) write(iwr,6025)
        ono = .true.
        return
       endif
       if(iunosp.gt.0) then
        minvb = iunosp
       endif
      endif
c
c     check for MCSCF calculation
c
      if (zruntp.ne.'ci'.and.zscftp.eq.'mcscf') then
       if(isecn.eq.0) then
        if(oroot) write(iwr,6015)
        ono = .true.
        return
       endif
       minva = isecn
      endif
c
      lenb = num*(num+1)/2
      l3 = num*num
      if (.not.odens) then
         if (minva.le.0) then
            call caserr('invalid section specified for input vectors')
         end if
         if(wfntyp.ne.mp2.and.zscft.eq.'uhf'.and.
     +      iuno.ne.0.and.iunosp.ne.0) then
          call dmano(veca,densa,iuno,iunopp,anos)
          if (zscft.eq.'uhf') call dmano(vecb,densb,iunosp,iunspp,bnos)
          old = otran
          go to 20
         endif
         call secget(minva,itvec,iblkv)
         call getqp(zcomm,zbuff,eiga,frocca,nbas,newbas,ncola,ieiga,
     +              idiffa,maxorb,iblkv)
         old = otran
         if(oroot) write (iwr,6030) ytype(1) , minva , ibl3d , 
     +                   (zcomm(7-i),i=1,6) , zbuff
         call rdedx(veca,l3,iblkv,idaf)
         if (.not.otran) call tdown(veca,ilifq,veca,ilifq,ncola)
         if ((zscft.eq.'gvb' .or. zscft.eq.'grhf').and.
     +        (zruntp.ne.'ci'))  then
            call reden2(veca,densa,densb,frocca,newbas,lenb)
            call vadd(densa,1,densb,1,densa,1,lenb)
         else
            call dmtx(densa,veca,frocca,iky,newbas,newbas,newbas)
         end if
         if (zscft.eq.'uhf'. or. 
     +       (zruntp.eq.'ci'. and. isping.gt.0) ) then
            if (minvb.le.0) then
             call caserr('invalid section specified for input vectors')
            else
               call secget(minvb,itvec,iblkv)
               call getqp(zcomm,zbuff,eigb,froccb,nbas,newbas,ncolb,
     +                    ieigb,idiffb,maxorb,iblkv)
               if(oroot) write (iwr,6030) ytype(2) , minvb , ibl3d , 
     +                         (zcomm(7-i),i=1,6) , zbuff
               call rdedx(vecb,l3,iblkv,idaf)
               if (.not.otran) call tdown(vecb,ilifq,vecb,ilifq,
     +             ncolb)
               call dmtx(densb,vecb,froccb,iky,newbas,newbas,newbas)
c              call vadd(densa,1,densb,1,densa,1,lenb)
            end if
         end if
      else
         if(oroot) write (iwr,6040) idena
         call secget(idena,m990,iblkv)
         call rdedx(densa,lenb,iblkv,idaf)
      end if
20    otran = .true.
      ibl7la= iposun(num8)
      otran = old
      return
6025  format(/10x,'spinfree NATORB section not specified')
6035  format(/10x,'CASSCF 1-particle density matrix not available')
6015  format(/10x,'MCSCF natural orbitals not available')
6000  format(/10x,'NOS not available for CCSD calculation'/)
6020  format(/10x,'NOS not available for FULL-CI calculation'/)
6030  format (//1x,a3,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6/
     +        ' header block information : '/
     +        ' vectors created under account ',a8/1x,a7,
     +        'vectors created by ',a8,' program at ',a8,' on ',a8,
     +        ' in the job ',a8/' with the title : ',10a8)
6040  format (/1x,'1-particle density matrix restored from section ',i3)
      end
      subroutine inrtia(am,cm,paxs,nocofm)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
c     ----- program to get inertial axes. constructs inertia -----
c           tensor, then calls diagonalization routine to get
c           direction cosines for principal inertia axes and
c           inertia moments. 
c
      parameter (zero = 0.0d+00)
      logical some
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      dimension am(*),cm(3),paxs(3,3)
      dimension a(6),pmom(3)
      character*1 dxyz(3),dabc(3)
      save nrtflg
      data nrtflg /0/
      data dxyz   /'x','y','z'/
      data dabc   /'a','b','c'/
c
      some=nprint.ne.-5
      if(nrtflg.gt.0) some=.false.
      do 10 i=1,6
         a(i)=zero
   10 continue
c
c     ----- calculate inertia tensor (wrt to centre of mass) -----
c
      do 20 iat=1,nat
         wt = am(iat)
         x = c(1,iat) - cm(1)
         y = c(2,iat) - cm(2)
         z = c(3,iat) - cm(3)
         a(1) = a(1) + wt*( y*y + z*z )
         a(2) = a(2) - wt*x*y
         a(3) = a(3) + wt*( x*x + z*z )
         a(4) = a(4) - wt*x*z
         a(5) = a(5) - wt*y*z
         a(6) = a(6) + wt*( x*x + y*y )
   20 continue
c
c     ----- diagonalize to find moments and principle axes -----
c
      no=3
      call diaaxsh(a,paxs,pmom,iky,no,no,no)
      if(some.and.nocofm.eq.0) write(iwr,9998) cm(1),cm(2),cm(3)
      if(some) write(iwr,9999)
      if(some) call praxsh(paxs,pmom,dxyz,dabc,no,0)
      nrtflg = 1
      return
 9999 format(/8x,' principle moments of inertia and principle axes',
     1       /8x,'       (transpose of axis rotation matrix)')
 9998 format(/2x,' centre of mass',
     1 ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
      end
      subroutine diaaxsh(a,vec,eig,ia,nvec,n,ndim)
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- routine to substitute diagiv for diagonalization -----
c           of symmetric 3x3 matrix a  in triangular form
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      dimension a(6),h(3,3),vec(3,3),eig(3),ia(*),big(6)
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data conv   /1.0d-10/
c
      if(n.eq.3.and.ndim.eq.3) go to 10
         write(iwr,9999)
         call caserr('program stop in -diaaxs-')
   10 continue
c
      do 30 i=1,3
         do 20 j=1,3
            vec(j,i)=zero
   20    continue
         vec(i,i)=one
   30 continue
      ij=0
      do 40 i=1,3
      do 40 j=1,i
         ij=ij+1
         h(i,j)=a(ij)
   40    h(j,i)=a(ij)
      call diajac(h,vec,eig,n,n,big)
c
c     ----- check for right handedness, correct if not -----
c
      test =   vec(1,3)*( vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2) )
     1       + vec(2,3)*( vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2) )
     2       + vec(3,3)*( vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2) )
      if(test.gt.zero) return
      if( abs(eig(1)-eig(2)).gt.conv) go to 60
         t = eig(1)
         eig(1) = eig(2)
         eig(2) = t
         do 50 i=1,3
            t = vec(i,1)
            vec(i,1) = vec(i,2)
            vec(i,2) = t
   50    continue
         return
   60 if( abs(eig(2)-eig(3)).gt.conv) go to 80
         t = eig(2)
         eig(2) = eig(3)
         eig(3) = t
         do 70 i=1,3
            t = vec(i,2)
            vec(i,2) = vec(i,3)
            vec(i,3) = t
   70    continue
         return
   80 do 90 i=1,3
         vec(i,3) = - vec(i,3)
   90 continue
      return
 9999 format(/' -diaaxs- diagonalization only set up for 3x3 matrix')
      end
      subroutine praxsh(v,e,dir1,dir2,ndim,list)
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- print out axis rotation matrix -----
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension v(ndim,*),e(*)
      character*1 dir1(3), dir2(3)
      write (iwr,9008)
      if(list.eq.0) write (iwr,9068) (e(i),i = 1,3)
      if(list.eq.1) write (iwr,9168) (e(i),i = 1,3)
      if(list.eq.2) write (iwr,9268) (e(i),i = 1,3)
      write (iwr,9008)
      if(list.eq.0) write (iwr,9028) (dir2(i),i = 1,3)
      if(list.eq.1) write (iwr,9128) (dir2(i),i = 1,3)
      if(list.eq.2) write (iwr,9228) (dir2(i),i = 1,3)
      write (iwr,9008)
      do 120 j = 1,3
      if(list.eq.0) write (iwr,9048) j,dir1(j),(v(j,i),i = 1,3)
      if(list.eq.1) write (iwr,9148) j,dir1(j),(v(j,i),i = 1,3)
      if(list.eq.2) write (iwr,9248) j,dir1(j),(v(j,i),i = 1,3)
  120 continue
      return
 9008 format(1x)
 9028 format(15x,10(7x,a1,3x))
 9048 format(i5,5x,a1,4x,10f11.5)
 9068 format(15x,10f11.3)
 9128 format(15x,7(9x,a1,5x))
 9148 format(i5,5x,a1,4x,7f15.10)
 9168 format(15x,7f15.10)
 9228 format(15x,7(9x,a1,5x))
 9248 format(i5,5x,a1,4x,7e15.8)
 9268 format(15x,7e15.8)
      end
      subroutine diajac(a,v,d,n,ndim,t)
      implicit real*8 (a-h,o-z), integer   (i-n)
c
      dimension a(ndim,ndim),v(ndim,ndim),d(ndim),t(ndim,2)
 
      call jacobih(a,v,d,n,ndim,nrot,t(1,1),t(1,2))
      call jacsrth(v,d,n,ndim)
 
      return
      end
      subroutine jacobih(a,v,d,n,np,nrot,b,z)
      implicit real*8 (a-h,o-z), integer   (i-n)
      dimension a(np,np),v(np,np),d(np),b(np),z(np)
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pt2    /0.2d+00/
      data pt5    /0.5d+00/
      data hundrd /1.0d+02/
c
      do 12 ip=1,n
         do 11 iq=1,n
            v(iq,ip)=zero
   11    continue
         v(ip,ip)=one
   12 continue
      do 13 ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=zero
   13 continue
      nrot=0
      if(n.eq.1) return
      do 24 i=1,50
         sm=zero
         do 15 ip=1,n-1
            do 14 iq=ip+1,n
               sm=sm+abs(a(ip,iq))
   14       continue
   15    continue
         if(sm.eq.zero) return
         if(i.lt.4) then
            thres=pt2*sm/n**2
         else
            thres=zero
         endif
         do 22 ip=1,n-1
            do 21 iq=ip+1,n
               g=hundrd*abs(a(ip,iq))
               if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.
     1                         (abs(d(iq))+g.eq.abs(d(iq)))) then
                  a(ip,iq)=zero
               else if(abs(a(ip,iq)).gt.thres) then
                  h=d(iq)-d(ip)
                  if(abs(h)+g.eq.abs(h)) then
                     t=a(ip,iq)/h
                  else
                     theta=pt5*h/a(ip,iq)
                     t=one/(abs(theta)+sqrt(one+theta**2))
                     if(theta.lt.zero) t=-t
                  endif
                  c=one/sqrt(one+t**2)
                  s=t*c
                  tau=s/(one+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=zero
                  do 16 j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
   16             continue
                  do 17 j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
   17             continue
                  do 18 j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
   18             continue
                  do 19 j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
   19             continue
                  nrot=nrot+1
               endif
   21       continue
   22    continue
         do 23 ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=zero
   23    continue
   24 continue
      call caserr(' iteration count in jacobih exceeded')
      return
      end
      subroutine jacsrth(v,d,n,np)
      implicit real*8 (a-h,o-z), integer   (i-n)
      dimension v(np,np),d(np)
      if(n.eq.1) return
      do 13 i=1,n-1
         k=i
         p=d(i)
         do 11 j=i+1,n
            if(d(j).ge.p) then
               k=j
               p=d(j)
            endif
   11    continue
         if(k.ne.i) then
            d(k)=d(i)
            d(i)=p
            do 12 j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
   12       continue
         endif
   13 continue
      return
      end
      subroutine dipaxs(a,d1,d2,d3,dmx,dmy,dmz)
c
c     ----- program to rotate dipole tensor coordinates -----
c           matrix a is transponse of rotation matrix
c           of direction cosines. 
c
      implicit real*8 (a-h,o-z), integer   (i-n)
      parameter (debye = 2.54176568d+00)
      parameter (zero = 0.0d+00)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*1 d1,d2,d3,d(3)
      dimension a(3,3),dip(3)
      d(1)=d1
      d(2)=d2
      d(3)=d3
      dip(1)=dmx/debye
      dip(2)=dmy/debye
      dip(3)=dmz/debye
      write(iwr,9000)
      do 20 na=1,3
         zz=zero
         do 10 i=1,3
            zz = zz + a(i,na)*dip(i)
   10    continue
         write(iwr,9001) d(na),zz,zz*debye
   20 continue
      write(iwr,9002)
      return
 9000 format(/10x,'component',9x,'bohr*e',15x,'debye')
 9001 format(13x,'mu ',a1,2f20.10)
 9002 format(1x)
      end
      subroutine qdpole(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, ono, oroot
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
      character*8 scf, uhf, mp2
      character*4 keyefc
      logical out, ofield, some
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      common/three/am(maxat),
     + qmn(6),qme(6),qmt(6),qmnefc(6),qmtefc(6),
     + qmnb(6),qmeb(6),qmtb(6),qmnefcb(6),qmtefcb(6),
     + thn(6),the(6),tht(6),thnefc(6),thtefc(6),
     + thnb(6),theb(6),thtb(6),thnefcb(6),thtefcb(6),
     + thdum(80)
      dimension x(*)
      character *2 xx
      dimension xx(7)
      dimension cm(3)
      dimension a(6),vec(3,3),eig(3)
      data scf,uhf /'scf     ','uhf     '/
      data mp2     /'mp2     '/
      data buck    /1.344911d+00/
      data tol     /1.0d-06/
      data zero    /0.0d+00/
      data two     /2.0d+00/
      data three   /3.0d+00/
      data keyefc  /' efc'/
      data xx / 'xx', 'yy', 'zz', 'xy', 'xz', 'yz', 'rr'/
c
      out=.false.      
      ono = .false.
      ofield = iefc.eq.keyefc
      oroot = opg_root()
      some = .true.
      some = some .and. oroot
c
c
      if (oroot) then
       write(iwr,9991)
       write(iwr,9993)
      endif
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
c
c     ----- allocate memory block -----
c
      need = l2+l2+l2+l2+l2+l2+l2+l2+l3+l3
      init = igmem_alloc(need)

      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of core -----
c
      i10=  init
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
      i50=i40+l2
      i60=i50+l2
      i70=i60+l2
c     ----- memory block for density matrix -----
      i80=i70+l2
      i90=i80+l2
c     ----- memory block for -a- vectors  -----
      i91=i90+l3
c     ----- memory block for -b- vectors -----
      last=i91+l3
c
      ncall=1
   10 continue
      call prpams(am,ncall,ncode,some)
      if(ncode.ne.0) go to 100
c
      sum=zero
      do 20 i=1,3
   20 cm(i)=zero
      do 30 iat=1,nat
      sum=sum+am(iat)
      do 30 i=1,3
   30 cm(i)=cm(i)+am(iat)*c(i,iat)
      if( abs(sum).lt.tol) go to 45
      do 40 i=1,3
   40 cm(i)=cm(i)/sum
   45 continue
      if(some) write(iwr,9998) cm(1),cm(2),cm(3)
c
c     use centre of mass, or origin nominated under gauge?
c
      if (orgin) then
       if (some) write(iwr,1234) (ggx(i),i=1,3)
       call dcopy(3,ggx(1),1,cm(1),1)
      else
       if (some) write(iwr,1235)
      endif
c
c     ----- calculate second moment integrals -----
c
      call qdpint(x(i10),x(i20),x(i30),
     1                x(i40),x(i50),x(i60),cm,l1)
c
c     ----- electronic contribution to second moments -----
c     
c     ----- density matrix for x(i70)=da,x(i80)=db -----

c     ----- calculate density matrix -----

      call prpdens(x(i70),x(i80),x(i90),x(i91),zscftp,ono,oroot)
      if (ono) return
c
c     xx=1, yy=2, zz=3, xy=4, xz=5, yz=6
c
      qme(1)=-tracep(x(i70),x(i10),l1)
      qme(2)=-tracep(x(i70),x(i20),l1)
      qme(3)=-tracep(x(i70),x(i30),l1)
      qme(4)=-tracep(x(i70),x(i40),l1)
      qme(5)=-tracep(x(i70),x(i50),l1)
      qme(6)=-tracep(x(i70),x(i60),l1)
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 50
       if(iuno.ne.0.and.iunosp.ne.0) go to 50
        qme(1)=-tracep(x(i80),x(i10),l1)+qme(1)
        qme(2)=-tracep(x(i80),x(i20),l1)+qme(2)
        qme(3)=-tracep(x(i80),x(i30),l1)+qme(3)
        qme(4)=-tracep(x(i80),x(i40),l1)+qme(4)
        qme(5)=-tracep(x(i80),x(i50),l1)+qme(5)
        qme(6)=-tracep(x(i80),x(i60),l1)+qme(6)
   50 continue
c
c     ----- < r**2 > = diamagnetic susceptibility -----
c
      rsquar=-(qme(1)+qme(2)+qme(3))
      if (oroot) write(iwr,9992) rsquar
c
c     ----- nuclear contribution
c
      do i = 1, 6
       qmn(i) = zero
       qmnefc(i) = zero
      enddo
c
      do 60 iat=1,nat
      qmn(1)=qmn(1)+czan(iat)*(c(1,iat)-cm(1))*(c(1,iat)-cm(1))
      qmn(2)=qmn(2)+czan(iat)*(c(2,iat)-cm(2))*(c(2,iat)-cm(2))
      qmn(3)=qmn(3)+czan(iat)*(c(3,iat)-cm(3))*(c(3,iat)-cm(3))
      qmn(4)=qmn(4)+czan(iat)*(c(1,iat)-cm(1))*(c(2,iat)-cm(2))
      qmn(5)=qmn(5)+czan(iat)*(c(1,iat)-cm(1))*(c(3,iat)-cm(3))
      qmn(6)=qmn(6)+czan(iat)*(c(2,iat)-cm(2))*(c(3,iat)-cm(3))
   60 continue
c
c     ----- -efc- contribution -----
c
      if(ofield) then
      do iat=1,nefc
       qmnefc(1)=qmnefc(1)+efcz(iat)*(efcc(1,iat)-cm(1))*
     +                               (efcc(1,iat)-cm(1))
       qmnefc(2)=qmnefc(2)+efcz(iat)*(efcc(2,iat)-cm(2))*
     +                               (efcc(2,iat)-cm(2))
       qmnefc(3)=qmnefc(3)+efcz(iat)*(efcc(3,iat)-cm(3))*
     +                               (efcc(3,iat)-cm(3))
       qmnefc(4)=qmnefc(4)+efcz(iat)*(efcc(1,iat)-cm(1))*
     +                               (efcc(2,iat)-cm(2))
       qmnefc(5)=qmnefc(5)+efcz(iat)*(efcc(1,iat)-cm(1))*
     +                               (efcc(3,iat)-cm(3))
       qmnefc(6)=qmnefc(6)+efcz(iat)*(efcc(2,iat)-cm(2))*
     +                               (efcc(3,iat)-cm(3))
      enddo
      endif
c
      do i = 1,6
       qmt(i) = qme(i) + qmn(i)
       qmtefc(i) = qmt(i) + qmnefc(i)
c       convert to buckinghams
       qmeb(i) = qme(i) * buck
       qmnb(i) = qmn(i) * buck
       qmtb(i) = qmt(i) * buck
       qmnefcb(i) = qmnefc(i) * buck
       qmtefcb(i) = qmtefc(i) * buck
      enddo
c
      if (oroot) then
        write(iwr,6000)
        if (ofield) then
         write(iwr,6005)
         do i = 1 , 6
          write (iwr,6010)xx(i),qmn(i),qmnefc(i),qme(i),qmtefc(i)
         enddo
         write(iwr,6090)
         do i = 1 , 6
          write (iwr,6010)xx(i),qmnb(i),qmnefcb(i),qmeb(i),qmtefcb(i)
         enddo
c
        else
         write (iwr,6040)
         write (iwr,6030)
         do i = 1 , 6
           write (iwr,6010)xx(i),qmn(i) ,qme(i) ,qmt(i)
         enddo
         rsquare=qme(1)+qme(2)+qme(3)
         rsquarn=qmn(1)+qmn(2)+qmn(3)
         rsquart=qmt(1)+qmt(2)+qmt(3)
         write (iwr,6011) xx(7), rsquarn, rsquare, rsquart
c     principal axis transformation
         call paxis2(qmn,qme,qmt,iwr,1)
c
         write(iwr,6090)
         write(iwr,6040)
         write(iwr,6030)
         do i = 1 , 6
           write (iwr,6010)xx(i),qmnb(i) ,qmeb(i) ,qmtb(i)
         enddo
c     principal axis transformation
         call paxis2(qmnb,qmeb,qmtb,iwr,1)
        endif
      endif
c
      dum = qmn(1) + qmn(2) + qmn(3)
      thn(1) = (three * qmn(1) - dum) / two
      thn(2) = (three * qmn(2) - dum) / two
      thn(3) = (three * qmn(3) - dum) / two
      thn(4) =  three * qmn(4) / two
      thn(5) =  three * qmn(5) / two
      thn(6) =  three * qmn(6) / two
c
      dum = qme(1) + qme(2) + qme(3)
      the(1) = (three * qme(1) - dum) / two
      the(2) = (three * qme(2) - dum) / two
      the(3) = (three * qme(3) - dum) / two
      the(4) =  three * qme(4) / two
      the(5) =  three * qme(5) / two
      the(6) =  three * qme(6) / two
c
      dum = qmnefc(1) + qmnefc(2) + qmnefc(3)
      thnefc(1) = (three * qmnefc(1) - dum) / two
      thnefc(2) = (three * qmnefc(2) - dum) / two
      thnefc(3) = (three * qmnefc(3) - dum) / two
      thnefc(4) =  three * qmnefc(4) / two
      thnefc(5) =  three * qmnefc(5) / two
      thnefc(6) =  three * qmnefc(6) / two
c
      do i = 1,6
       tht(i) = the(i) + thn(i)
       thtefc(i) = tht(i) + thnefc(i)
c       convert to buckinghams
       theb(i) = the(i) * buck
       thnb(i) = thn(i) * buck
       thtb(i) = tht(i) * buck
       thnefcb(i) = thnefc(i) * buck
       thtefcb(i) = thtefc(i) * buck
      enddo
      if (oroot) then
       write(iwr,1235) 
       write(iwr,7000)
        if (ofield) then
         write(iwr,6005)
         do i = 1 , 6
          write (iwr,6010)xx(i),thn(i),thnefc(i),the(i),thtefc(i)
         enddo
         write(iwr,7090)
         do i = 1 , 6
          write (iwr,6010)xx(i),thnb(i),thnefcb(i),theb(i),thtefcb(i)
         enddo
c
        else
         write (iwr,6040)
         write (iwr,6030)
         do i = 1 , 6
           write (iwr,6010)xx(i),thn(i) ,the(i) ,tht(i)
         enddo
c     principal axis transformation - quadrupole
         call paxis2(thn,the,tht,iwr,1)
c
         write (iwr,7090)
         write (iwr,6040)
         write (iwr,6030)
         do i = 1 , 6
           write (iwr,6010)xx(i),thnb(i) ,theb(i) ,thtb(i)
         enddo
c     principal axis transformation - quadrupole
         call paxis2(thnb,theb,thtb,iwr,1)
        endif
      endif
c
      ncall=ncall+1
c
c
  100 continue
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
c
 9991 format(/10x,29('-')
     +       /10x,'second and quadrupole moments',
     +       /10x,29('-')/)
 9993 format(/' 1 a.u. = 1.344911 buckinghams ',
     +        '= 1.344911 10**(-26) esu*cm**2 ')
 9998 format(/2x,' centre of mass',
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 9996 format(' not enough core in -qdpole- . stop ')
 6000 format(/2x,
     + 'second moments in atomic units (w.r.t. above origin)'/2x,
     + '----------------------------------------------------')
 6040 format (4x,'**** expectation values --- ',
     +             ' global axis coordinates')
 6090 format(/2x,
     + 'second moments in Buckinghams (w.r.t. above origin)'/2x,
     + '--------------------------------------------------'/)
 6005 format(/12x,'nuclear',11x,'-efc-',6x,'electronic',11x,'total'/)
 6030 format(/12x,'nuclear',6x,'electronic',11x,'total'/)
 6010 format(1x,a2,4f16.7)
 6011 format(/1x,a2,4f16.7)
 7000 format(/2x,
     + 'quadrupole moments in atomic units (w.r.t. above origin)'/2x,
     + '-------------------------------------------------------')
 7090 format(/2x,
     + 'quadrupole moments in Buckinghams (w.r.t. above origin)'/2x,
     + '-------------------------------------------------------'/)
 9992 format(/' < r**2 > = ',f10.6,' a.u. ',
     + ' ( 1 a.u. = 0.280023 10**(-16) cm**2 ) '/
     + ' ( also called diamagnetic susceptibility ) ')
 1234 format (/1x,'origin for property evaluation:',
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 1235 format (/1x,'origin for property evaluation is centre of mass')
      end
      subroutine paxis2(pnucin,pelein,ptotin,iwr,mm)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension map(6),pnucin(*),pelein(*),ptotin(*)
      common/miscop/p(3,3),prnct(3,3),prnc(3,3),q(3,3),buf(6),
     + bufn(6),bufe(6),pnuc(6),pele(6),ptot(6),u(3,3),iq(3,3)
c     xx=1, yy=2, zz=3, xy=4, xz=5, yz=6
      character *2 xx
      dimension xx(6)
      data xx / 'xx', 'yy', 'zz', 'xy', 'xz', 'yz'/
      data map/1,4,5,2,6,3/
      data m3,nr/3,3/
c
      write(iwr,6000)
      do i = 1 , 6
       pnuc(i) = pnucin(i)
       pele(i) = pelein(i)
       ptot(i) = ptotin(i)
      enddo
      ij = 0
      do i = 1 , 3
         do j = i , 3
            ij = ij + 1
            prnct(i,j) = ptot(map(ij))
         enddo
      enddo
      call hdiag(prnct,u,p,iq,3,3,0,nr)
      ij = 0
      do i = 1 , 3
         do j = i , 3
            ij = ij + 1
            ptot(map(ij)) = prnct(i,j)
         enddo
      enddo
      call mwt(u,m3,m3,m3)
c     nuclear (n=1) then electronic (n=2)
      do n = 1 , 2
         ij = 0
         do i = 1 , 3
            do j = i , 3
               ij = ij + 1
               if(n.eq.1) then
               prnc(i,j) = pnuc(map(ij))
               else
               prnc(i,j) = pele(map(ij))
               endif
            enddo
         enddo
         call mabat(u,prnc,m3,m3,q,m3)
         ij = 0
         do i = 1 , 3
            do j = i , 3
               ij = ij + 1
               if(n.eq.1) then
                pnuc(map(ij)) = prnc(i,j)
               else
                pele(map(ij)) = prnc(i,j)
               endif
            enddo
         enddo
         if (n.eq.1) then
            do m = 1 , 6
               bufn(m) = pnuc(m)
            enddo
         else
            do m = 1 , 6
               bufe(m) = pele(m)
            enddo
         end if
      enddo

      do m = 1 , 6
         buf(m) = ptot(m)
      enddo
      write (iwr,6030)
      do ij = 1 , 6
        write (iwr,6080)xx(ij),bufn(ij),bufe(ij),buf(ij)
      enddo

      if (mm.eq.1) then
         pp = dabs(buf(1))
         qq = dabs(buf(2))
         rr = dabs(buf(3))
         if (pp.lt.qq) then
            temp = qq
            qq = pp
            pp = temp
         end if
         if (qq.lt.rr) then
            temp = rr
            rr = qq
            qq = temp
            if (pp.lt.qq) then
               temp = qq
               qq = pp
               pp = temp
            end if
         end if
         if (pp.ne.0.0d0) then
            asym = (qq-rr)/pp
            write (iwr,6010) asym
         end if
      end if

      write (iwr,6060)
      do n = 1 , 3
         pp = buf(n)
         write (iwr,6050) pp , (u(n,j),j=1,3)
      enddo
c
      return
 6000 format (/4x,'**** expectation values --- ',
     +             ' principal axis coordinates'/)
 6010 format (/1x,'asymmetry parameter =',f10.3)
 6030 format(12x,'nuclear',6x,'electronic',11x,'total'/)
c    +        f15.7/)
 6050 format (f12.6,5x,3(f12.6,3x))
 6060 format (/1x,'rotation matrix'//3x,'eigenvalue',23x,
     +        'eigenvectors'/)
 6080 format(1x,a2,4f16.7)
      end
      subroutine qdpint(xxs,yys,zzs,xys,xzs,yzs,cm,l1)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical out
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
      common/xyzqdp/xint,yint,zint,xintx,yinty,zintz,
     1              xintxx,yintyy,zintzz,t,xc,yc,zc,
     2              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension xxs(*),yys(*),zzs(*),xys(*),xzs(*),yzs(*)
      common/junk3/sxx(225),syy(225),szz(225),
     +             sxy(225),sxz(225),syz(225),
     +             dij(225),
     +             xin(5,5),  yin(5,5),  zin(5,5),
     +             xxin(5,5), yyin(5,5), zzin(5,5),
     +             xxxin(5,5),yyyin(5,5),zzzin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension cm(3)
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      norm=normf.ne.1.or.normp.ne.1
c
      xc=cm(1)
      yc=cm(2)
      zc=cm(3)
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
      sxx(ij)=zero
      syy(ij)=zero
      szz(ij)=zero
      sxy(ij)=zero
      sxz(ij)=zero
      syz(ij)=zero
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
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
c     ----- quadrupole moment integrals -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljt
      nj=j
      do 370 i=1,lit
      ni=i
      call qdpxyz
        xin(i,j)=xint*t1
        yin(i,j)=yint*t1
        zin(i,j)=zint*t1
       xxin(i,j)=xintx*t1
       yyin(i,j)=yinty*t1
       zzin(i,j)=zintz*t1
      xxxin(i,j)=xintxx*t1
      yyyin(i,j)=yintyy*t1
      zzzin(i,j)=zintzz*t1
  370 continue
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 390 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      ij=ij+1
      sxx(ij)=sxx(ij)+dij(ij)*(xxxin(ix,jx)*  yin(iy,jy)*  zin(iz,jz))
      syy(ij)=syy(ij)+dij(ij)*(  xin(ix,jx)*yyyin(iy,jy)*  zin(iz,jz))
      szz(ij)=szz(ij)+dij(ij)*(  xin(ix,jx)*  yin(iy,jy)*zzzin(iz,jz))
      sxy(ij)=sxy(ij)+dij(ij)*( xxin(ix,jx)* yyin(iy,jy)*  zin(iz,jz))
      sxz(ij)=sxz(ij)+dij(ij)*( xxin(ix,jx)*  yin(iy,jy)* zzin(iz,jz))
      syz(ij)=syz(ij)+dij(ij)*(  xin(ix,jx)* yyin(iy,jy)* zzin(iz,jz))
  390 continue
 6000 continue
 7000 continue
c
c     ----- set up quadrupole moment matrices -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=iky(loci+i)+(locj+j)
      xxs(nn)=sxx(ij)
      yys(nn)=syy(ij)
      zzs(nn)=szz(ij)
      xys(nn)=sxy(ij)
      xzs(nn)=sxz(ij)
      yzs(nn)=syz(ij)
 7500 continue
c
 8000 continue
 9000 continue
c
      if(out) then
         write(iwr,9999)
         call writel(xxs,l1)
         write(iwr,9998)
         call writel(yys,l1)
         write(iwr,9997)
         call writel(zzs,l1)
         write(iwr,9996)
         call writel(xys,l1)
         write(iwr,9995)
         call writel(xzs,l1)
         write(iwr,9994)
         call writel(yzs,l1)
      endif
c
      return
 9999 format(/10x,26('-')/10x,'xx-second moment integrals'/
     +         10x,26('-'))
 9998 format(/10x,26('-')/10x,'yy-second moment integrals'/
     +         10x,26('-'))
 9997 format(/10x,26('-')/10x,'zz-second moment integrals'/
     +         10x,26('-'))
 9996 format(/10x,26('-')/10x,'xy-second moment integrals'/
     +         10x,26('-'))
 9995 format(/10x,26('-')/10x,'xz-second moment integrals'/
     +         10x,26('-'))
 9994 format(/10x,26('-')/10x,'yz-second moment integrals'/
     +         10x,26('-'))
      end
      subroutine qdpxyz
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzqdp/xint,yint,zint,xintx,yinty,zintz,
     1              xintxx,yintyy,zintzz,t,xc,yc,zc,
     2              x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min /1,2,4, 7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      xintx=zero
      yinty=zero
      zintz=zero
      xintxx=zero
      yintyy=zero
      zintzz=zero
      npts=(ni+nj-2+2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
      xintx=xintx+px*(ptx-xc)
      yinty=yinty+py*(pty-yc)
      zintz=zintz+pz*(ptz-zc)
      xintxx=xintxx+px*(ptx-xc)*(ptx-xc)
      yintyy=yintyy+py*(pty-yc)*(pty-yc)
      zintzz=zintzz+pz*(ptz-zc)*(ptz-zc)
   16 continue
      return
      end
      subroutine efgmap(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, oroot
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
      logical some,out, ono
      character*8 scf, uhf, mp2
      character*4 keyefc
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
      dimension x(*)
      dimension a(6),vec(3,3),eig(3)
      data keyefc /' efc'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data three  /3.0d+00/
      data tenm3  /1.0d-03/
      data angle  /180.0d+00/
      data uhf    /'uhf     '/
      data mp2    /'mp2     '/
      data scf    /'scf     '/
c
c     ----- calculate electric field gradient -----
c
      oroot = opg_root()
c
      if (oroot) write(iwr,9999)
      if (oroot) write(iwr,9994)
      some=.true.
      out =.false.
      ono = .false.
c
      pi  = acos(-one)
      deg = angle/pi
c
c     ----- set block for partitioning of memory -----
c
      nptmax = nat
      l1     = num
      l2     = (l1*(l1+1))/2
      l3     = l1*l1
      i00    =          1
      i10    = i00    + nptmax*3
      i20    = i10    + nptmax*6
      i30    = i20    + (num*(num+1))/2
      i40    = i30    + (num*(num+1))/2
      i50    = i40    + l3
      i60    = i50    + l3
      last   = i60
      need   = last   - 1
c
c     ----- allocate memory block -----
c
      init = igmem_alloc(need)

      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of memory -----
c
      i00    =          init
      i10    = i00    + nptmax*3
      i20    = i10    + nptmax*6
      i30    = i20    + (num*(num+1))/2
      i40    = i30    + (num*(num+1))/2
      i50    = i40    + l3
      i60    = i50    + l3
c
c     ----- define points for calculation -----
c           1. nuclei
c
      nptmax=0
      nptmax=nptmax+nat
      do 30 iat=1,nat
         x(init  +3*(iat-1))=c(1,iat)
         x(init+1+3*(iat-1))=c(2,iat)
         x(init+2+3*(iat-1))=c(3,iat)
   30 continue
c
c     ----- density matrix for x(i20)=da,x(i30)=db(0?) -----
c

c     ----- calculate density matrix -----

      call prpdens(x(i20),x(i30),x(i40),x(i50),zscftp,ono,oroot)
      if (ono) return

      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
       if(iuno.eq.0.or.iunosp.eq.0) then
         do 120 i=1,l2
            x(i+i20-1)=x(i+i20-1)+x(i+i30-1)
  120    continue
       endif
      endif
c
c     ----- calculate electric field gradient at all points -----
c
      call elfgrd(nptmax,x(i00),x(i10),x(i20))
c
c     ----- output results of all points -----
c
      if (oroot) write(iwr,9997)
c
      do 230  ipt=1,nptmax
         xp = x(1+3*(ipt-1)+i00-1)
         yp = x(2+3*(ipt-1)+i00-1)
         zp = x(3+3*(ipt-1)+i00-1)
c
c     ----- form nuclear contribution -----
c
         efgrxx = zero
         efgryy = zero
         efgrzz = zero
         efgrxy = zero
         efgrxz = zero
         efgryz = zero
         do 210 i = 1,nat
            xn = c(1,i) - xp
            yn = c(2,i) - yp
            zn = c(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.tenm3) go to 210
               rr5=rr*rr*rr*rr*rr
               efgrxx = efgrxx - czan(i)*xn*xn/rr5
               efgryy = efgryy - czan(i)*yn*yn/rr5
               efgrzz = efgrzz - czan(i)*zn*zn/rr5
               efgrxy = efgrxy - czan(i)*xn*yn/rr5
               efgrxz = efgrxz - czan(i)*xn*zn/rr5
               efgryz = efgryz - czan(i)*yn*zn/rr5
  210    continue
c
c     ----- form -efc- contribution -----
c
      if(iefc.eq.keyefc) then
         do 220 i = 1,nefc
            xn = efcc(1,i) - xp
            yn = efcc(2,i) - yp
            zn = efcc(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.tenm3.and.oroot) 
     &        write(iwr,9993) xp,yp,zp,i
            if(rr.lt.tenm3) go to 220
               rr5=rr*rr*rr*rr*rr
               efgrxx = efgrxx - efcz(i)*xn*xn/rr5
               efgryy = efgryy - efcz(i)*yn*yn/rr5
               efgrzz = efgrzz - efcz(i)*zn*zn/rr5
               efgrxy = efgrxy - efcz(i)*xn*yn/rr5
               efgrxz = efgrxz - efcz(i)*xn*zn/rr5
               efgryz = efgryz - efcz(i)*yn*zn/rr5
  220    continue
      endif
c
       efgrxx  = x(1+6*(ipt-1)+i10-1)/three + efgrxx
       efgryy  = x(2+6*(ipt-1)+i10-1)/three + efgryy
       efgrzz  = x(3+6*(ipt-1)+i10-1)/three + efgrzz
       efgrxy  = x(4+6*(ipt-1)+i10-1)/three + efgrxy
       efgrxz  = x(5+6*(ipt-1)+i10-1)/three + efgrxz
       efgryz  = x(6+6*(ipt-1)+i10-1)/three + efgryz
c
       efgxx  =   two*efgrxx - efgryy - efgrzz
       efgyy  =   two*efgryy - efgrxx - efgrzz
       efgzz  =   two*efgrzz - efgrxx - efgryy
       efgxy  = three*efgrxy
       efgxz  = three*efgrxz
       efgyz  = three*efgryz
c
       if (oroot) write(iwr,9995)
     + ipt,zaname(ipt),efgxx,efgyy,efgzz,efgxy,efgxz,efgyz
c
       x(1+6*(ipt-1)+i10-1) = efgxx
       x(2+6*(ipt-1)+i10-1) = efgyy
       x(3+6*(ipt-1)+i10-1) = efgzz
       x(4+6*(ipt-1)+i10-1) = efgxy
       x(5+6*(ipt-1)+i10-1) = efgxz
       x(6+6*(ipt-1)+i10-1) = efgyz
c
  230  continue
c
      if (oroot) write(iwr,9992)
c
      do 240 ipt=1,nptmax
      a(1) = x(1+6*(ipt-1)+i10-1)
      a(2) = x(4+6*(ipt-1)+i10-1)
      a(3) = x(2+6*(ipt-1)+i10-1)
      a(4) = x(5+6*(ipt-1)+i10-1)
      a(5) = x(6+6*(ipt-1)+i10-1)
      a(6) = x(3+6*(ipt-1)+i10-1)
      call diaefg(a,vec,eig,3)
      eta  = abs( (eig(1)-eig(2)) / eig(3) )
      if (oroot) write(iwr,9991)
     &   ipt,zaname(ipt),eig(1),eig(2),eig(3),eta
  240 continue
c
      if (oroot) write(iwr,9990)
c
      do ipt=1,nptmax
      a(1) = x(1+6*(ipt-1)+i10-1)
      a(2) = x(4+6*(ipt-1)+i10-1)
      a(3) = x(2+6*(ipt-1)+i10-1)
      a(4) = x(5+6*(ipt-1)+i10-1)
      a(5) = x(6+6*(ipt-1)+i10-1)
      a(6) = x(3+6*(ipt-1)+i10-1)
      call diaefg(a,vec,eig,3)
c     do j=1,3
c      do i=1,3
c       vec(i,j)= acos(vec(i,j))*deg
c      enddo
c     enddo
      if (oroot) then
       write(iwr,9989) ipt,zaname(ipt)
       write(iwr,9987)
       do j = 1,3
        write(iwr,9988) eig(j),(vec(i,j),i=1,3)
       enddo
      endif
      enddo
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
 9999 format(/10x,23('-')/10x,'electric field gradient'
     1       /10x,23('-')/)
 9998 format(' not enough core in -efgmap-')
 9997 format(/2x,104('-')/
     + 1x,' atom',5x,8x,1x,
     + 19x,'electric field gradient in molecular frame (a.u.)'/
     + 17x,8x,'xx',13x,'yy',13x,'zz',13x,'xy',13x,'xz',13x,'yz'/
     + 2x,104('-'))
 9996 format(' --- warning - electric field gradient at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(1x,i3,1x,a8,3x,6f15.6)
 9994 format(' 1 a.u. = 0.324123 10**(16) esu/cm**3 ',
     1       ' ( or statvolts/cm**2 )',' = 0.97174 10**(22) v/m**2 ')
 9993 format(' --- warning - electric field gradient at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
 9992 format(/1x,' atom',11x,'principal components of the -efg- ',
     1       'tensor (a.u.)',9x,'asymmetry parameter eta'/
     2       2x,94('-'))
 9991 format(1x,i3,1x,a8,1x,3f15.6,14x,f15.6)
 9990 format(/' rotation of the principal axis',
     1        ' of -efg- tensor w.r.t absolute (molecular) frame')
 9989 format(/' atom = ',i5, ' (',1x,a8,')'/' ------ ')
 9987 format(/1x,'rotation matrix'/3x,'eigenvalue',23x,
     +        'eigenvectors'/)
 9988 format(f12.6,5x,3(f12.6,3x))
      end
      subroutine diaefg(a,vec,eig,ndim)
      implicit real*8 (a-h,o-z), integer   (i-n)
      dimension a(6),h(3,3),vec(3,3),eig(3),big(6)
      data zero   /0.0d+00/
      data one    /1.0d+00/
c
      do 20 i=1,3
         do 10 j=1,3
   10       vec(j,i)=zero
   20    vec(i,i)=one
      ij=0
      do 25 i=1,3
      do 25 j=1,i
         ij=ij+1
         h(i,j)=a(ij)
   25    h(j,i)=a(ij)
      call diajac(h,vec,eig,ndim,ndim,big)
c
c     ----- put the principal components of the tensor and -----
c           the corresponding eigenvectors in increasing
c           order according to the absolute value of the
c           eigenvalue ( the principal component ) .
c
      do 50 i=1,3
         jj=i
         do 30 j=i,3
            if( abs(eig(j)).lt. abs(eig(jj))) jj=j
   30       continue
         if(jj.eq.i) go to 50
c
         xx=eig(jj)
         eig(jj)=eig(i)
         eig(i)=xx
         do 40 j=1,3
            xx=vec(j,jj)
            vec(j,jj)=vec(j,i)
            vec(j,i)=xx
   40       continue
c
   50    continue
      return
      end
      subroutine elfgrd(npt,xyzpt,elfgr,dab)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical some,out
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
      common/root/xx,u(12),w(12),nroots
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension w2(6),w4(6)
      dimension xyzpt(3,*),elfgr(6,*),dab(*)
      common /junk3/ dij(225),
     +               xv(5,5,6),  yv(5,5,6),  zv(5,5,6),
     +              dxv(5,5,6), dyv(5,5,6), dzv(5,5,6),
     +             ddxv(5,5,6),ddyv(5,5,6),ddzv(5,5,6)
      dimension ijx(35),ijy(35),ijz(35)
      data maxrys /6/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data four   /4.0d+00/
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
      out =nprint.eq.6
      some=nprint.ne.-5
      norm=normf.ne.1.or.normp.ne.1
c
c
      do 10 ipt=1,npt
      do 10 i=1,6
   10 elfgr(i,ipt)=zero
c
      nder=2
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(iwr,9997) maxrys,lit,ljt,nroots
         call caserr('program stop in -elfgrd-')
      endif
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      nn=iky(loci+i)+(locj+j)
      den=dab(nn)
      if(.not.iandj.or.j.ne.i) den=den+den
      ij=ij+1
  360 dij(ij)=dum2*den
c
c     ----- electric field gradient term -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ipt=1,npt
      znuc=one
      cx=xyzpt(1,ipt)
      cy=xyzpt(2,ipt)
      cz=xyzpt(3,ipt)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      if (nroots.ge.6) call rootss
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      u2=uu
      u4=uu*uu
      ww=w(iroot)*znuc
      w2(iroot)=ww*u2*two
      w4(iroot)=ww*u4*four
      tt=one/(aa+uu)
      t = sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i
      call dsxyzh
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint
      call dvxyzh
      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint
      call ddvxyzh
      ddxv(i,j,iroot)=xint
      ddyv(i,j,iroot)=yint
      ddzv(i,j,iroot)=zint
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
      dumx=zero
      dumy=zero
      dumz=zero
      dumxx=zero
      dumyy=zero
      dumzz=zero
      dumxy=zero
      dumxz=zero
      dumyz=zero
      do 430 iroot=1,nroots
      dumx =dumx+
     1    dxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumy =dumy+
     1     xv(ix,jx,iroot)* dyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumz =dumz+
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)* dzv(iz,jz,iroot)*w2(iroot)
      dum  =
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w2(iroot)
      dumxx=dumxx-dum+
     1   ddxv(ix,jx,iroot)*  yv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumyy=dumyy-dum+
     1     xv(ix,jx,iroot)*ddyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumzz=dumzz-dum+
     1     xv(ix,jx,iroot)*  yv(iy,jy,iroot)*ddzv(iz,jz,iroot)*w4(iroot)
      dumxy=dumxy+
     1    dxv(ix,jx,iroot)* dyv(iy,jy,iroot)*  zv(iz,jz,iroot)*w4(iroot)
      dumxz=dumxz+
     1    dxv(ix,jx,iroot)*  yv(iy,jy,iroot)* dzv(iz,jz,iroot)*w4(iroot)
      dumyz=dumyz+
     1     xv(ix,jx,iroot)* dyv(iy,jy,iroot)* dzv(iz,jz,iroot)*w4(iroot)
  430 continue
      ij=ij+1
      dum=pi212*aa1*dij(ij)
      elfgr(1,ipt)=elfgr(1,ipt)+dumxx*dum
      elfgr(2,ipt)=elfgr(2,ipt)+dumyy*dum
      elfgr(3,ipt)=elfgr(3,ipt)+dumzz*dum
      elfgr(4,ipt)=elfgr(4,ipt)+dumxy*dum
      elfgr(5,ipt)=elfgr(5,ipt)+dumxz*dum
      elfgr(6,ipt)=elfgr(6,ipt)+dumyz*dum
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
 8000 continue
 9000 continue
c
c
      return
 9997 format(' in -elfgrd- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3/
     2       ' lit,ljt,nroots = ',3i3)
      end
      subroutine dsxyzh
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj
     +                 ,ni,nj,cx,cy,cz
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min  /1,2,4, 7,11,16,22/
      data max  /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj-2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)*t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
   16 continue
      return
      end
c
      subroutine dvxyzh
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj
     +                 ,ni,nj,cx,cy,cz
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min  /1,2,4, 7,11,16,22/
      data max  /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj+1-2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 13 i=imin,imax
      dum=h(i)*t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      px=ptx-cx
      py=pty-cy
      pz=ptz-cz
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 go to (12,11,10,9,8,7),nj
    7 px=px*bx
      py=py*by
      pz=pz*bz
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 dum=w(i)
      xint=xint+dum*px
      yint=yint+dum*py
      zint=zint+dum*pz
   13 continue
      return
      end
c
      subroutine ddvxyzh
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj
     +                 ,ni,nj,cx,cy,cz
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min  /1,2,4, 7,11,16,22/
      data max  /1,3,6,10,15,21,28/
      data zero /0.0d+00/
      xint = zero
      yint = zero
      zint = zero
      npts = (ni+nj+2-2)/2+1
      imin = min(npts)
      imax = max(npts)
      do 11 i = imin,imax
      dum = h(i)*t
      ptx = dum+x0
      pty = dum+y0
      ptz = dum+z0
      px = ptx-cx
      py = pty-cy
      pz = ptz-cz
      px = px*px
      py = py*py
      pz = pz*pz
      ax = ptx-xi
      ay = pty-yi
      az = ptz-zi
      bx = ptx-xj
      by = pty-yj
      bz = ptz-zj
      go to (5,4,3,2,1),ni
    1 px = px*ax
      py = py*ay
      pz = pz*az
    2 px = px*ax
      py = py*ay
      pz = pz*az
    3 px = px*ax
      py = py*ay
      pz = pz*az
    4 px = px*ax
      py = py*ay
      pz = pz*az
    5 go to (10,9,8,7,6),nj
    6 px = px*bx
      py = py*by
      pz = pz*bz
    7 px = px*bx
      py = py*by
      pz = pz*bz
    8 px = px*bx
      py = py*by
      pz = pz*bz
    9 px = px*bx
      py = py*by
      pz = pz*bz
   10 dum = w(i)
      xint = xint+dum*px
      yint = yint+dum*py
      zint = zint+dum*pz
   11 continue
      return
      end
      subroutine elfmap(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root,oroot
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
      logical some,out,ono
      character*8 scf, uhf, mp2
      character*4 keyefc
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
      dimension x(*)
      data keyefc /' efc'/
      data zero   /0.0d+00/
      data tenm3  /1.0d-03/
      data uhf    /'uhf     '/
      data mp2    /'mp2     '/
      data scf    /'scf     '/
c
c     ----- calculate electric field -----
c
      oroot = opg_root()
      if (oroot) write(iwr,9999)
      if (oroot) write(iwr,9994)
      some=.true.
      out =.false.       
      ono =.false.
c
c     ----- set block for partitioning of memory -----
c
      nptmax = nat
      l1     = num
      l2     = (l1*(l1+1))/2
      l3     = l1*l1
      i00    =          1
      i10    = i00    + nptmax*3
      i20    = i10    + nptmax*3
      i30    = i20    + (num*(num+1))/2
      i40    = i30    + (num*(num+1))/2
      i50    = i40    + l3
      i60    = i50    + l3
      last   = i60
      need   = last   - 1
c
c     ----- allocate memory block -----
c
      init = igmem_alloc(need)

      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of memory -----
c
      i00    =          init
      i10    = i00    + nptmax*3
      i20    = i10    + nptmax*3
      i30    = i20    + (num*(num+1))/2
      i40    = i30    + (num*(num+1))/2
      i50    = i40    + l3
      i60    = i50    + l3
c
c     ----- define points for calculation -----
c           1. nuclei
c
      nptmax=0
c
      nptmax=nptmax+nat
      do 30 iat=1,nat
         x(init  +3*(iat-1))=c(1,iat)
         x(init+1+3*(iat-1))=c(2,iat)
         x(init+2+3*(iat-1))=c(3,iat)
   30 continue
c
c     ----- density matrix for x(i20)=da,x(i30)=db -----
c

c     ----- calculate density matrix -----

      call prpdens(x(i20),x(i30),x(i40),x(i50),zscftp,ono,oroot)
      if (ono) return

c
c     ----- get total density matrix -----
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
       if(iuno.eq.0.or.iunosp.eq.0) then
         do 120 i=1,l2
            x(i+i20-1)=x(i+i20-1)+x(i+i30-1)
  120    continue
       endif
      endif
c
c     ----- calculate electric field at all points -----
c
      call elfldc(nptmax,x(i00),x(i10),x(i20))
c
c     ----- output results of all points -----
c
      if (oroot) write(iwr,9997)
c
      do 230  ipt=1,nptmax
         xp = x(1+3*(ipt-1)+i00-1)
         yp = x(2+3*(ipt-1)+i00-1)
         zp = x(3+3*(ipt-1)+i00-1)
c
c     ----- form nuclear contribution -----
c
         elfldx = zero
         elfldy = zero
         elfldz = zero
         do 210 i = 1,nat
            xn = c(1,i) - xp
            yn = c(2,i) - yp
            zn = c(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.tenm3) go to 210
               rr3=rr*rr*rr
               elfldx = elfldx - czan(i)*xn/rr3
               elfldy = elfldy - czan(i)*yn/rr3
               elfldz = elfldz - czan(i)*zn/rr3
  210    continue
c
c     ----- form -efc- contribution -----
c
      if(iefc.eq.keyefc) then
         do 220 i = 1,nefc
            xn = efcc(1,i) - xp
            yn = efcc(2,i) - yp
            zn = efcc(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.tenm3.and.oroot)
     &         write(iwr,9993) xp,yp,zp,i
            if(rr.lt.tenm3) go to 220
               rr3=rr*rr*rr
               elfldx = elfldx - efcz(i)*xn/rr3
               elfldy = elfldy - efcz(i)*yn/rr3
               elfldz = elfldz - efcz(i)*zn/rr3
  220    continue
      endif
c
         elfldx  = x(1+3*(ipt-1)+i10-1) + elfldx
         elfldy  = x(2+3*(ipt-1)+i10-1) + elfldy
         elfldz  = x(3+3*(ipt-1)+i10-1) + elfldz
         elfld   = sqrt(elfldx*elfldx + elfldy*elfldy + elfldz*elfldz)
         if (oroot) write(iwr,9995) 
     &      ipt,zaname(ipt),elfldx,elfldy,elfldz,elfld
  230    continue

c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
 9999 format(/10x,14('-')/10x,'electric field',
     1       /10x,14('-')/)
 9998 format(' not enough core in -elfmap-')
 9997 format(/2x,73('-')/
     + 3x,' atom',3x,8x,1x,
     + 15x,'electric field (a.u.)'/
     + 16x,8x,'Ex',13x,'Ey',13x,'Ez',11x,'Field'/
     + 2x,73('-'))
 9996 format(' --- warning - electric field at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(1x,i3,1x,a8,2x,4f15.7)
 9994 format(' 1 a.u. = 0.171524 10**(-8) dyn/esu ')
 9993 format(' --- warning - electric field at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
      end
      subroutine elfldc(npt,xyzpt,elfld,dab)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical some,out
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
      common/root/xx,u(12),w(12),nroots
      common/xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     +                  ni,nj,cx,cy,cz
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension elfld(3,1),xyzpt(3,1),dab(1)
      common /junk3/ dij(225),
     + xv(5,5,5), yv(5,5,5), zv(5,5,5),
     + dxv(5,5,5),dyv(5,5,5),dzv(5,5,5)
      dimension ijx(35),ijy(35),ijz(35)
      data maxrys /6/
      data rln10  /2.30258d+00/
      data zero   /0.0d+00/
      data one    /1.0d+00/
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
      out =nprint.eq.6
      some=nprint.ne.-5
      norm=normf.ne.1.or.normp.ne.1
c
      do 10 ipt=1,npt
      do 10 i=1,3
   10 elfld(i,ipt)=zero
c
      nder=1
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt+nder-2)/2+1
      if(nroots.gt.maxrys) then
         write(iwr,9997) maxrys,lit,ljt,nroots
         call caserr('program stop in -elfldc-')
      endif
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
      go to 350
  320 if(norm) dum2=dum2*sqrt7
      go to 350
  330 if(norm) dum2=dum2*sqrt5/sqrt3
      go to 350
  340 if(norm) dum2=dum2*sqrt3
  350 continue
c
      nn=iky(loci+i)+(locj+j)
      den=dab(nn)
      if(.not.iandj.or.j.ne.i) den=den+den
      ij=ij+1
  360 dij(ij)=dum2*den
c
c     ----- electric field term -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ipt=1,npt
      znuc=one
      cx=xyzpt(1,ipt)
      cy=xyzpt(2,ipt)
      cz=xyzpt(3,ipt)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      if (nroots.ge.6) call rootss
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      ww=ww*(uu+uu)
      tt=one/(aa+uu)
      t = sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i
      call dsxyzh
      xv(i,j,iroot)=xint
      yv(i,j,iroot)=yint
      zv(i,j,iroot)=zint*ww
      call dvxyzh
      dxv(i,j,iroot)=xint
      dyv(i,j,iroot)=yint
      dzv(i,j,iroot)=zint*ww
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
      dumx=zero
      dumy=zero
      dumz=zero
      do 430 iroot=1,nroots
      dumx=dumx+dxv(ix,jx,iroot)* yv(iy,jy,iroot)* zv(iz,jz,iroot)
      dumy=dumy+ xv(ix,jx,iroot)*dyv(iy,jy,iroot)* zv(iz,jz,iroot)
      dumz=dumz+ xv(ix,jx,iroot)* yv(iy,jy,iroot)*dzv(iz,jz,iroot)
  430 continue
      ij=ij+1
      dum=pi212*aa1*dij(ij)
      elfld(1,ipt)=elfld(1,ipt)+dumx*dum
      elfld(2,ipt)=elfld(2,ipt)+dumy*dum
      elfld(3,ipt)=elfld(3,ipt)+dumz*dum
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
 8000 continue
 9000 continue
c
c
      return
 9997 format(' in -elfldc- , the rys quadrature is not implemented',
     1       ' beyond -nroots- = ',i3/
     2       ' lit,ljt,nroots= ',3i3)
      end
      subroutine elpmap(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, oroot
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
      character*8 scf, uhf, mp2
      character*4 keyefc
      logical some,out,ono
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
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      dimension x(*)
      common/three/am(maxat),thdum(200)
      common/junk2/xyzmap(3,100),zznuc(maxat),nuc(maxat)
      dimension cm(3)
      dimension mo(100),moa(100)
      dimension xrange(3),yrange(3),zrange(3),xvec(3),yvec(3),zvec(3)
      dimension origin(3)
      dimension iatom(4)
      equivalence (mo(1),moa(1))
      data keyefc /' efc'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data tenm3  /1.0d-03/
      data ten12  /1.0d+12/
      data uhf    /'uhf     '/
      data mp2    /'mp2     '/
      data scf    /'scf     '/
      data moa    /100*0/
      data nptxyz /100/
      data iunit  /-1/
      data bohr   /0.529177249d+00/
      data xrange /0.0d+00,0.0d+00,1.0d+00/
      data yrange /0.0d+00,0.0d+00,1.0d+00/
      data zrange /0.0d+00,0.0d+00,1.0d+00/
      data xvec   /1.0d+00,0.0d+00,0.0d+00/
      data yvec   /0.0d+00,1.0d+00,0.0d+00/
      data zvec   /0.0d+00,0.0d+00,1.0d+00/
      data origin /0.0d+00,0.0d+00,0.0d+00/
      data iatom  /0,0,0,0/
c
c
c     ----- calculate electrostatic potential -----
c
      oroot = opg_root()
      if (oroot) write(iwr,9999)
      if (oroot) write(iwr,9994)
      some=.true.
      out =.false.     
      out = out .and. oroot
      ono =.false.
      call vfill(ten12,xyzmap,1,300)
c
c     ----- set pointers for partitioning of memory -----
c
      l1 = num
      l2 = (l1*(l1+1))/2
      l3 = l1*l1
c note: nptmax calculated based on zero grid points 
      nptmax = nat +1
c
      i00 =          1
      i10 = i00    + nptmax*3
      i20 = i10    + nptmax
      i30 = i20    + (num*(num+1))/2
      i40 = i30    + (num*(num+1))/2
c     ----- memory block for -a- vectors  -----
      i50 = i40+l3
c     ----- memory block for -b- vectors  -----
      i60 = i50+l3
      last= i60
c     need= last   - 1
      need= last
c
c     ----- allocate memory block -----
c
      init = igmem_alloc(need)
      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of memory -----
c
      i00 =          init
      i10 = i00    + nptmax*3
      i20 = i10    + nptmax
      i30 = i20    + (num*(num+1))/2
      i40 = i30    + (num*(num+1))/2
c     ----- memory block for -a- vectors  -----
      i50 = i40    +  l3
c     ----- memory block for -b- vectors  -----
      i60 = i50    +  l3
c
c     ----- density matrix for x(i20)=da,x(i30)=db -----
c

c     ----- calculate density matrix -----

      call prpdens(x(i20),x(i30),x(i40),x(i50),zscftp,ono,oroot)
      if (ono) return
c
c     ----- define points for calculation -----
c           1. grid points
c           2. specific points -xyzmap-
c           3. nuclei
c           4. centre of mass
c
      nptmax=0
      nptgrd=0
      nptmax=nptmax+nptgrd
c
c
c     if(nptgrd.gt.0.and.oroot) then
c        write(iwr,9985) nxgrd,nygrd,nzgrd,
c    1    iatom,origin,xvec,xrange(1),xrange(2),xrange(3),
c    2                 yvec,yrange(1),yrange(2),yrange(3),
c    3                 zvec,zrange(1),zrange(2),zrange(3)
c     else
c        if (out) write(iwr,9986)
c     endif
c
      fac=one
      if(iunit.gt.0) fac=bohr
      npt=0
      do 20 ipt=1,nptxyz
      if(xyzmap(1,ipt).eq.ten12.or.xyzmap(2,ipt).eq.ten12.or.
     1   xyzmap(3,ipt).eq.ten12) go to 20
      npt=npt+1
      x(init  +3*(npt-1+nptmax))=xyzmap(1,ipt)/fac
      x(init+1+3*(npt-1+nptmax))=xyzmap(2,ipt)/fac
      x(init+2+3*(npt-1+nptmax))=xyzmap(3,ipt)/fac

   20 continue

      nptmax=nptmax+npt

c
      do 30 iat=1,nat
      x(init  +3*(iat-1+nptmax))=c(1,iat)
      x(init+1+3*(iat-1+nptmax))=c(2,iat)
      x(init+2+3*(iat-1+nptmax))=c(3,iat)

   30 continue
      nptmax=nptmax+nat
c
      ncall=1
      call prpams(am,ncall,ncode,out)
      xcm=zero
      ycm=zero
      zcm=zero
      tcm=zero
      do 40 iat=1,nat
      xcm=xcm+am(iat)*c(1,iat)
      ycm=ycm+am(iat)*c(2,iat)
      zcm=zcm+am(iat)*c(3,iat)
   40 tcm=tcm+am(iat)
      cm(1)=xcm/tcm
      cm(2)=ycm/tcm
      cm(3)=zcm/tcm
c
c     use centre of mass, or origin nominated under gauge?
c
      if (orgin) then
       call dcopy(3,ggx(1),1,cm(1),1)
      endif
      x(init  +3*(1-1+nptmax))=cm(1)
      x(init+1+3*(1-1+nptmax))=cm(2)
      x(init+2+3*(1-1+nptmax))=cm(3)
      nptmax=nptmax+1
c
c
c     ----- get total density matrix -----
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
       if(iuno.eq.0.or.iunosp.eq.0) then
         do 120 i=1,l2
            x(i+i20-1)=x(i+i20-1)+x(i+i30-1)
  120    continue
        endif
      endif

c
c     ----- calculate electronic contribution at all points -----
c
      call elpotc(nptmax,x(i00),x(i10),x(i20))
c
c     ----- output results of all points -----
c
      if (oroot) write(iwr,9997)
c
c
      do 230  ipt=1,nptmax
         xp = x(1+3*(ipt-1)+i00-1)
         yp = x(2+3*(ipt-1)+i00-1)
         zp = x(3+3*(ipt-1)+i00-1)
c
c     ----- form nuclear contribution -----
c
         elpotn = zero
         do 210 i = 1,nat
            xn = c(1,i) - xp
            yn = c(2,i) - yp
            zn = c(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(oroot.and.
     1         (rr.lt.tenm3).and.(ipt.lt.(nptmax-1-nat).or.
     1                          ipt.gt.(nptmax-1    )    ) )
     2      write(iwr,9996) xp,yp,zp,i
            if(rr.lt.tenm3) go to 210
               elpotn = elpotn + czan(i)/rr
  210    continue
c
c     ----- form -efc- contribution -----
c
      if(iefc.eq.keyefc) then
         do 220 i = 1,nefc
            xn = efcc(1,i) - xp
            yn = efcc(2,i) - yp
            zn = efcc(3,i) - zp
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.tenm3.and.oroot)
     &         write(iwr,9993) xp,yp,zp,i
            if(rr.lt.tenm3) go to 220
               elpotn = elpotn + efcz(i)/rr
  220    continue
      endif
c
         elpote = x(ipt+i10-1)
         elpott = - elpote + elpotn
         x(ipt+i10-1)=elpott
         if(oroot.and.(ipt.gt.nptgrd.or.out)) then
           if (ipt.lt.nptmax) then
           write(iwr,9995) ipt,zaname(ipt),xp,yp,zp,elpott,elpote
           else
           write(iwr,9984) ipt,xp,yp,zp,elpott,elpote
           endif
         endif
  230    continue
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
 9999 format(/10x,45('-')
     1       /10x,'electrostatic potential/diamagnetic shielding',
     2       /10x,45('-')/)
 9997 format(/3x,89('-')/
     +       3x,'point',13x,'x',9x,'y',9x,'z',5x,'potential(a.u.)',
     +       3x,'diamagnetic shielding(a.u.)'/
     +       3x,89('-'))
 9996 format(' --- warning - electrostatic potential at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9995 format(1x,i3,1x,a8,2x,3f10.5,f15.6,6x,f15.6)
 9984 format(1x,i3,11x,3f10.5,f15.6,6x,f15.6)
 9994 format(' 1 a.u. = 9.07618 esu/cm ( or statvolts ) ')
 9993 format(' --- warning - electrostatic potential at ',
     1 3f10.5,' . contribution from  -efc-  ',i3,' ignored')
 9990 format(3f13.6,i5)
 9989 format(6f13.6,2x)
c9986 format(/' no grid defined. ')
c9985 format(/1x,'( ',i4,' by ',i4,' by ',i4,' ) grid defined.',
c    1 /' iatom  = ',4i4,
c    2 /' origin = ',3f10.5,
c    3 /' xvec   = ',3f10.5,' xmin=',f8.2,' xmax=',f8.2,' xstep=',f8.2,
c    4 /' yvec   = ',3f10.5,' ymin=',f8.2,' ymax=',f8.2,' ystep=',f8.2,
c    5 /' zvec   = ',3f10.5,' zmin=',f8.2,' zmax=',f8.2,' zstep=',f8.2,
c    6 /)
      end
      subroutine elpotc(npt,xyzpt,elpot,dab)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical out
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
      common/root/xx,u(12),w(12),nroots
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     &                  ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension elpot(*),xyzpt(3,*),dab(*)
      common /junk3/ dij(225),
     +               xv(5,5,5),yv(5,5,5),zv(5,5,5)
      dimension ijx(35),ijy(35),ijz(35)
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data pi212 /1.1283791670955d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
c
      do 10 ipt=1,npt
   10 elpot(ipt)=zero
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      nroots=(lit+ljt-2)/2+1
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
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
c     ----- electrostatic potential -----
c
      aax=aa*ax
      aay=aa*ay
      aaz=aa*az
      do 500 ipt=1,npt
      znuc=one
      cx=xyzpt(1,ipt)
      cy=xyzpt(2,ipt)
      cz=xyzpt(3,ipt)
      xx=aa*((ax-cx)**2+(ay-cy)**2+(az-cz)**2)
      if(nroots.le.3) call rt123
      if(nroots.eq.4) call roots4
      if(nroots.eq.5) call roots5
      do 420 iroot=1,nroots
      uu=u(iroot)*aa
      ww=w(iroot)*znuc
      tt=one/(aa+uu)
      t= sqrt(tt)
      x0=(aax+uu*cx)*tt
      y0=(aay+uu*cy)*tt
      z0=(aaz+uu*cz)*tt
      do 410 j=1,ljt
      nj=j
      do 410 i=1,lit
      ni=i
      call sxyzh
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
      dum=dum*dij(ij)
      den=dab(iky(loci+i)+(locj+j))
      if(.not.iandj.or.j.ne.i) den=den+den
      elpot(ipt)=elpot(ipt)+dum*den
  440 continue
  450 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
 8000 continue
 9000 continue
c
      return
      end
      subroutine sxyzh
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     &                  ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min  /1,2,4, 7,11,16,22/
      data max  /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj-2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)*t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
   16 continue
      return
      end
      subroutine eldmap(x,zscft)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, oroot
      character *8 zscft
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
      character*8 rhf, uhf, mp2, scf, prop
c     character*8 gvb
      logical some,out,dbug,ono
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
      dimension x(*)
      common/three/am(maxat),thdum(200)
      common/junk2/xyzmap(3,100),zznuc(maxat),nuc(maxat)
      dimension cm(3)
      dimension mo(100),moa(100),mob(100)
      dimension xrange(3),yrange(3),zrange(3),xvec(3),yvec(3),zvec(3)
      dimension origin(3)
      dimension iatom(4)
      equivalence (mo(1),moa(1))
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      data ten12  /1.0d+12/
      data prop   /'propty  '/
      data rhf    /'rhf     '/
      data uhf    /'uhf     '/
      data mp2    /'mp2     '/
c     data gvb    /'gvb     '/
      data scf    /'scf     '/
      data moa    /100*0/
      data mob    /100*0/
      data maxmo  /100/
      data nptxyz /100/
      data iunit  /-1/
      data bohr   /0.529177249d+00/
      data xrange /0.0d+00,0.0d+00,1.0d+00/
      data yrange /0.0d+00,0.0d+00,1.0d+00/
      data zrange /0.0d+00,0.0d+00,1.0d+00/
      data xvec   /1.0d+00,0.0d+00,0.0d+00/
      data yvec   /0.0d+00,1.0d+00,0.0d+00/
      data zvec   /0.0d+00,0.0d+00,1.0d+00/
      data origin /0.0d+00,0.0d+00,0.0d+00/
      data iatom  /0,0,0,0/
c
c     ----- calculate electron or spin density -----
c
      oroot = opg_root()
      if (oroot) write(iwr,9999)
      if (oroot) write(iwr,9996)
      dbug=.false.
      out =.false.       
      out =(out.or.dbug).and.oroot
      ono=.false.
      some=.true.
      some=some.or.out
      if(zruntp.eq.prop) then
         nco=na
      endif
      call vfill(ten12,xyzmap,1,300)
c
c     ----- define points for calculation -----
c           1. grid points
c           2. specific points -xyzmap-
c           3. nuclei
c           4. centre of mass
c
      nptmax=0
      nptgrd=0
      nptmax=nptmax+nptgrd
c
c     if(nptgrd.gt.0.and.oroot) then
c        write(iwr,9985) nxgrd,nygrd,nzgrd,
c    1    iatom,origin,xvec,xrange(1),xrange(2),xrange(3),
c    2                 yvec,yrange(1),yrange(2),yrange(3),
c    3                 zvec,zrange(1),zrange(2),zrange(3)
c     else
c        if (out) write(iwr,9986)
c     endif
c
      fac=one
      if(iunit.gt.0) fac=bohr
      npt=0
      do 20 ipt=1,nptxyz
      if(xyzmap(1,ipt).eq.ten12.or.xyzmap(2,ipt).eq.ten12.or.
     1   xyzmap(3,ipt).eq.ten12) go to 20
      npt=npt+1
   20 continue
      nptmax=nptmax+npt
c
      nptmax=nptmax+nat
c
      ncall=1
      call prpams(am,ncall,ncode,out)
      xcm=zero
      ycm=zero
      zcm=zero
      tcm=zero
      do 40 iat=1,nat
      xcm=xcm+am(iat)*c(1,iat)
      ycm=ycm+am(iat)*c(2,iat)
      zcm=zcm+am(iat)*c(3,iat)
   40 tcm=tcm+am(iat)
      cm(1)=xcm/tcm
      cm(2)=ycm/tcm
      cm(3)=zcm/tcm
c
c     use centre of mass, or origin nominated under gauge?
c
      if (orgin) then
       call dcopy(3,ggx(1),1,cm(1),1)
      endif
c
      nptmax=nptmax+1
c
c     ----- we are going to use the general case -----
c
      l1= num
      l2=(num*(num+1))/2
      l3= num* num
      need= l2 + l2 + l3 + l3 + 4*nptmax
c
c     ----- allocate memory block -----
c
      init = igmem_alloc(need)
c
      call vclr(x(init),1,need)
c
      nptmax=0           
      do 30 iat=1,nat
      x(init  +3*(iat-1+nptmax))=c(1,iat)
      x(init+1+3*(iat-1+nptmax))=c(2,iat)
      x(init+2+3*(iat-1+nptmax))=c(3,iat)
   30 continue
      nptmax=nptmax+nat
      x(init  +3*(1-1+nptmax))=cm(1)
      x(init+1+3*(1-1+nptmax))=cm(2)
      x(init+2+3*(1-1+nptmax))=cm(3)
      nptmax=nptmax+1
c
      if(dbug) then
         write(6,*) 'need   = ',need
         write(6,*) 'runtyp = ',zruntp
         write(6,*) 'wfntyp = ',wfntyp
         write(6,*) 'scftyp = ',zscft
      endif
c
c--   if(wfntyp.eq.scf.and.zscftp.ne.gvb) go to 200
c
c     ----- general case = one particle reduced density matrix -----
c           is not diagonal.
c
      if(zruntp.eq.prop) then
         if (oroot) write(iwr,9983)
         call caserr('program stop in -eldmap-')
      endif
      if(dbug) then
         write(6,*) 'general case ...'
      endif
c
c     ----- set pointers for partitioning of memory -----
c
      l1= num
      l2=(num*(num+1))/2
      l3= num* num
      i00 =          init
      i10 = i00    + nptmax*3
      i20 = i10    + nptmax
      i30 = i20    + (num*(num+1))/2
      i40 = i30    + (num*(num+1))/2
c     ----- memory block for -a- vectors  -----
      i50 = i40    + l3
c     ----- memory block for -b- vectors  -----
      i60 = i50    + l3
c
c     ----- get total density matrix -----
c     ----- calculate density matrix -----

      call prpdens(x(i20),x(i30),x(i40),x(i50),zscft,ono,oroot)
      if (ono) return

      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
         if (iuno.gt.0.and.iunosp.gt.0) then
c     x(i20) and x(i30) contain the spinfree and spin density matrix
         else
c     x(i20) and x(i30) contain the -a- and -b- density matrix
         do 120 i=1,l2
            duma      =x(i+i20-1)
            dumb      =x(i+i30-1)
            x(i+i20-1)=duma+dumb             
            x(i+i30-1)=duma-dumb             
  120    continue
         endif
      else if (zruntp.eq.'ci'.and.isping.gt.0) then
          if(oroot) write(iwr,9979)
      else
      endif
c
c     ----- calculate electronic contribution at all points -----
c
      call eldenc(nptmax,x(i00),x(i10),x(i20))
c
c     ----- output results for all points -----
c
      if (oroot) then
         write(iwr,9993)
         write(iwr,9997)
c
      endif
c
      do 130  ipt=1,nptmax
         xp = x(1+3*(ipt-1)+i00-1)
         yp = x(2+3*(ipt-1)+i00-1)
         zp = x(3+3*(ipt-1)+i00-1)
         eldene = x(ipt+i10-1)
         if(oroot.and.(ipt.gt.nptgrd.or.out)) then
            if (ipt.eq.nptmax) then
            write(iwr,9975) ipt,xp,yp,zp,eldene
            else
            write(iwr,9995) ipt,zaname(ipt),xp,yp,zp,eldene
            endif
         endif
  130    continue
c
c     ----- spin density now -----
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf).or.
     2   (zruntp.eq.'ci'.and.isping.gt.0)) then
         call eldenc(nptmax,x(i00),x(i10),x(i30))
         if (oroot) then
            write(iwr,9992)
            write(iwr,9997)
         endif
         do 140  ipt=1,nptmax
            xp = x(1+3*(ipt-1)+i00-1)
            yp = x(2+3*(ipt-1)+i00-1)
            zp = x(3+3*(ipt-1)+i00-1)
            eldene = x(ipt+i10-1)
            if(oroot.and.(ipt.gt.nptgrd.or.out)) then
               if (ipt.eq.nptmax) then
               write(iwr,9975) ipt,xp,yp,zp,eldene
               else
               write(iwr,9995) ipt,zaname(ipt),xp,yp,zp,eldene
               endif
            endif
  140    continue
      endif
      go to 500
  200 continue
c
c     ----- special case = one particle reduced density matrix -----
c           is diagonal.
c
      l1= num
      l2=(num*(num+1))/2
      l3= num* num
      i00 =    init
      i10 =i00+nptmax*3
      i11 =i10+nptmax
      i12 =i11+nptmax
      i20 =i12+nptmax
      i30 =i20+num*num
      i40  = i30    + num*num
c     ----- memory block for vectors  -----
      i50  = i40    +  l3
c     ----- memory block for eigen values -----
      last  = i50    +  l1
      need=last-1
c
c     ----- read molecular orbitals -----
c
c     ----- get mo vectors from file -----
c     ----- read eig-vec x(i40); eng-val in x(i50) for tmp. -----

       call prpvect(x(i20),x(i30),x(i50),x(i50),
     +                     zscft,oroot)
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
      endif
c
c     ----- individual (alpha) orbital density -----
c
      do 230 imo=1,maxmo
      if(mo(imo).eq.0) go to 230
      if (oroot) then
         if(mo(imo).gt.0) write(iwr,9994) abs(mo(imo))
         if(mo(imo).gt.0) write(iwr,9997)
         if(mo(imo).lt.0) write(iwr,9982) abs(mo(imo))
         if(mo(imo).lt.0) write(iwr,9980)
c
      endif
c
      call elwfnc(nptmax,x(i00),x(i12),x(i20+
     1                l1*(abs(mo(imo))-1)))
c
      do 220 ipt=1,nptmax
      if(mo(imo).gt.0) x(ipt+i12-1)=x(ipt+i12-1)**2
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      eldene=x(ipt+i12-1)
      if(oroot.and.ipt.gt.nptgrd.or.out)
     &   write(iwr,9975) ipt,xp,yp,zp,eldene
  220 continue
  230 continue
c
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 300
c
c     ----- individual (beta) orbital density -----
c
      do 250 imo=1,maxmo
      if (oroot) then
         if(mob(imo).eq.0) go to 250
         if(mob(imo).gt.0) write(iwr,9984) abs(mob(imo))
         if(mob(imo).gt.0) write(iwr,9997)
         if(mob(imo).lt.0) write(iwr,9981) abs(mob(imo))
         if(mob(imo).lt.0) write(iwr,9980)
c
      endif
c
      call elwfnc(nptmax,x(i00),x(i12),x(i30+
     1                l1*(abs(mob(imo))-1)))
c
      do 240 ipt=1,nptmax
      if(mob(imo).gt.0) x(ipt+i12-1)=x(ipt+i12-1)**2
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      eldene=x(ipt+i12-1)
      if(oroot.and.ipt.gt.nptgrd.or.out) 
     &  write(iwr,9975) ipt,xp,yp,zp,eldene
  240 continue
  250 continue
c
  300 continue
c
c     ----- total electron density and spin density -----
c
      do 310 ipt=1,nptmax
      x(ipt+i10-1)=zero
  310 x(ipt+i11-1)=zero
c
      do 330 imo=1,na
      call elwfnc(nptmax,x(i00),x(i12),x(i20+l1*(imo-1)))
      if((wfntyp.eq.scf.and.zscftp.eq.rhf).and.imo.le.nco) occ=two
      if((wfntyp.eq.scf.and.zscftp.eq.rhf).and.imo.gt.nco) occ=one
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)               ) occ=one
      do 320 ipt=1,nptmax
      x(ipt+i10-1)=x(ipt+i10-1)+occ*x(ipt+i12-1)**2
  320 x(ipt+i11-1)=x(ipt+i11-1)+occ*x(ipt+i12-1)**2
  330 continue
c
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 360
      do 350 imo=1,nb
      call elwfnc(nptmax,x(i00),x(i12),x(i30+l1*(imo-1)))
      occ=one
      do 340 ipt=1,nptmax
      x(ipt+i10-1)=x(ipt+i10-1)+occ*x(ipt+i12-1)**2
  340 x(ipt+i11-1)=x(ipt+i11-1)-occ*x(ipt+i12-1)**2
  350 continue
  360 continue
c
c     ----- print total electron and spin density -----
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
       if (oroot) then
         write(iwr,9992)
         write(iwr,9997)
c
       endif
c
      do 370 ipt=1,nptmax
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      spnden=x(ipt+i11-1)
      if(oroot.and.(ipt.gt.nptgrd.or.out)) then
        if (ipt.eq.nptmax) then
         write(iwr,9975) ipt,xp,yp,zp,spnden
        else
         write(iwr,9995) ipt,zaname(ipt),xp,yp,zp,spnden
        endif
      endif
  370 continue
      endif
c
      if (oroot) then
         write(iwr,9993)
         write(iwr,9997)
c
      endif
c
      do 380 ipt=1,nptmax
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      eldene=x(ipt+i10-1)
      if(oroot.and.(ipt.gt.nptgrd.or.out)) then
        if (ipt.eq.nptmax) then
         write(iwr,9975) ipt,xp,yp,zp,eldene
        else
         write(iwr,9995) ipt,zaname(ipt),xp,yp,zp,eldene
        endif
      endif
  380 continue
c
  500 continue
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
 9999 format(/10x,21('-')/10x,'electron/spin density'
     1       /10x,21('-')/)
 9998 format(' not enough core in -eldmap-')
 9997 format(/1x,59('-')/
     +        3x,'point',13x,'x',9x,'y',9x,'z',6x,'density (a.u.)'/
     +       1x,59('-'))
 9996 format(' 1 a.u. = 4.80286 10**(-10) esu ')
 9975 format(1x,i3,11x,3f10.5,f15.6)
 9995 format(1x,i3,1x,a8,2x,3f10.5,f15.6)
 9994 format(/' density for (alpha) orbital = ',i5/
     1        ' ----------------------------------- ')
 9993 format(/' total electron density '/' ---------------------- ')
 9992 format(/' total spin density     '/' ------------------     ')
 9990 format(3f13.6,i5)
 9989 format(6f13.6,2x)
c9986 format(/' no grid defined. ')
c9985 format(/1x,'( ',i4,' by ',i4,' by ',i4,' ) grid defined.',
c    1 /' iatom  = ',4i4,
c    2 /' origin = ',3f10.5,
c    3 /' xvec   = ',3f10.5,' xmin=',f8.2,' xmax=',f8.2,' xstep=',f8.2,
c    4 /' yvec   = ',3f10.5,' ymin=',f8.2,' ymax=',f8.2,' ystep=',f8.2,
c    5 /' zvec   = ',3f10.5,' zmin=',f8.2,' zmax=',f8.2,' zstep=',f8.2,
c    6 /)
 9984 format(/' density for ( beta) orbital = ',i5/
     1        ' ----------------------------------- ')
 9983 format(/' for -runtyp- .eq. -propty  - there is no total',
     1         ' density. stop')
 9982 format(/' wavefunction for (alpha) orbital = ',i5/
     1         ' ---------------------------------------- ')
 9981 format(/' wavefunction for ( beta) orbital = ',i5/
     1         ' ---------------------------------------- ')
 9980 format(3x,'point',6x,'x',9x,'y',9x,'z',6x,'  wfn   (a.u.)')
 9979 format(1x,'evaluate spin density based on spin NOs')
      end
      subroutine elwfnc(npt,xyzpt,elden,v)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical norm
      logical out
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
      common/xyzwfn/xwfn,ywfn,zwfn,x0,y0,z0,xi,yi,zi,ni
      dimension elden(*),xyzpt(3,*),v(*)
      common /junk3/ di(15),xg(5),yg(5),zg(5)
      dimension ijx(35),ijy(35),ijz(35)
      data zero  /0.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
      do 10 ipt=1,npt
   10 elden(ipt)=zero
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
c     ----- i primitive -----
c
      do 7000 ig=i1,i2
      ai=ex(ig)
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c     ----- density factor -----
c
      inum=0
      do 230 i=mini,maxi
      go to (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),i
  110 dum1=csi
      go to 220
  120 dum1=cpi
      go to 220
  130 dum1=cdi
      go to 220
  140 if(norm) dum1=dum1*sqrt3
      go to 220
  150 dum1=cfi
      go to 220
  160 if(norm) dum1=dum1*sqrt5
      go to 220
  170 if(norm) dum1=dum1*sqrt3
      go to 220
  180 dum1=cgi
      go to 220
  190 if(norm) dum1=dum1*sqrt7
      go to 220
  200 if(norm) dum1=dum1*sqrt5/sqrt3
      go to 220
  210 if(norm) dum1=dum1*sqrt3
  220 continue
c
      inum=inum+1
  230 di(inum)=dum1*v(loci+i)
c
c     ----- loop over points to be evaluated -----
c
      do 500 ipt=1,npt
      x0 = xyzpt(1,ipt)
      y0 = xyzpt(2,ipt)
      z0 = xyzpt(3,ipt)
c
      dum = ai*((x0-xi)**2+(y0-yi)**2+(z0-zi)**2)
      if(dum.gt.tol) go to 500
      fac = exp(-dum)
c
c     ----- wavefunction values -----
c
      do 370 i=1,lit
      ni=i
      call wfnxyz
      xg(i)=xwfn
      yg(i)=ywfn
      zg(i)=zwfn
  370 continue
c
      inum=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      inum=inum+1
      wfn=fac*di(inum)*xg(ix)*yg(iy)*zg(iz)
      elden(ipt)=elden(ipt) + wfn
  390 continue
c
  500 continue
c
 7000 continue
c
 9000 continue
c
      return
      end
      subroutine wfnxyz
      implicit real*8 (a-h,o-z), integer   (i-n)
      common/xyzwfn/xwfn,ywfn,zwfn,x0,y0,z0,xi,yi,zi,ni
      data one /1.0d+00/
c
      xwfn=one
      ywfn=one
      zwfn=one
      ptxi=x0-xi
      ptyi=y0-yi
      ptzi=z0-zi
      go to (7,6,5,4,3,2,1),ni
    1 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    2 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    3 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    4 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    5 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    6 xwfn=xwfn*ptxi
      ywfn=ywfn*ptyi
      zwfn=zwfn*ptzi
    7 continue
      return
      end
      subroutine eldenc(npt,xyzpt,elden,dab)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical out
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
      common/xyzden/xint,yint,zint,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension elden(*),xyzpt(3,npt),dab(*)
      common /junk3/ dij(225),
     +               xin(5,5),yin(5,5),zin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
c
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      tol=rln10*itol
      norm=normf.ne.1.or.normp.ne.1
c
c
      do 10 ipt=1,npt
   10 elden(ipt)=zero
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
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
c     ----- electron density -----
c
      do 500 ipt=1,npt
      x0 = xyzpt(1,ipt)
      y0 = xyzpt(2,ipt)
      z0 = xyzpt(3,ipt)
c
      dum = aa*((x0-ax)**2+(y0-ay)**2+(z0-az)**2)
      if(dum.gt.tol) go to 500
      fac = exp(-dum)
c
c     ----- density values -----
c
      do 370 j=1,ljt
      nj=j
      do 370 i=1,lit
      ni=i
      call denint2
      xin(i,j)=xint
      yin(i,j)=yint
      zin(i,j)=zint
  370 continue
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 390 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      ij=ij+1
      nn=iky(loci+i)+(locj+j)
      eld=fac*dij(ij)*xin(ix,jx)*yin(iy,jy)*zin(iz,jz)
      den=dab(iky(loci+i)+(locj+j))
      if(.not.iandj.or.j.ne.i) den=den+den
      elden(ipt)=elden(ipt)+eld*den
  390 continue
c
  500 continue
c
 6000 continue
 7000 continue
c
 8000 continue
 9000 continue
c
c
      return
      end
      subroutine denint2
      implicit real*8 (a-h,o-z), integer   (i-n)
      common/xyzden/xint,yint,zint,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj
      data one /1.0d+00/
c
      xint=one
      yint=one
      zint=one
      ptxi=x0-xi
      ptyi=y0-yi
      ptzi=z0-zi
      ptxj=x0-xj
      ptyj=y0-yj
      ptzj=z0-zj
      go to (7,6,5,4,3,2,1),ni
    1 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    2 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    3 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    4 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    5 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    6 xint=xint*ptxi
      yint=yint*ptyi
      zint=zint*ptzi
    7 go to (14,13,12,11,10,9,8),nj
    8 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
    9 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   10 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   11 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   12 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   13 xint=xint*ptxj
      yint=yint*ptyj
      zint=zint*ptzj
   14 continue
      return
      end
      subroutine prpvect(veca,vecb,eiga,eigb,zscft,oroot)
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
      integer ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs
      integer ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb
      integer ibl7x, ibl7y, ibl7z
      integer ibl7so,ibl7o, ibl7la,ibl7ec,ibl7ha
      common /scra7/ ibl7s, ibl7t, ibl7f, ibl7st,ibl7tt,ibl7ft,ibl3qs,
     +               ibl7qa,ibl7pa,ibl7ea,ibl7qb,ibl7pb,ibl7eb,
     +               ibl7x, ibl7y, ibl7z,
     +               ibl7so,ibl7o(8), ibl7la, ibl7ec, ibl7ha
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
c
      logical oming, oextg, odirec, ovcfre, osemi, ostopm, olevd
      logical osym_op,osym_dens
      integer mina, minb, mouta, moutb, lock, ibrk
      integer numg, isecg, iblk3g, mextra, nsymo, iorbsy, isunor
      real*8  gapa1, gapa2, gapb1, gapb2, scaleg, symtag
      common/atmol3/mina,minb,mouta,moutb,lock,ibrk,
     + numg(2),isecg(2),iblk3g(2),mextra,nsymo,iorbsy(100),
     + oming, oextg,
     + gapa1,gapa2,gapb1,gapb2,scaleg,symtag(9),odirec(50),
     + isunor,ovcfre,osemi,ostopm,olevd,osym_op,osym_dens
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/bufb/frocca(maxorb),froccb(maxorb)
      common/junkc/zcom1(29),zcomm(19),zbuff(10),zatom(maxat)
c
      real*8 ctran
      integer ilifc, ntran, itran
      logical otran,otri
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     +            ctran(mxorb3),otran,otri
c
c
      integer mp, lwordp, lvec, ied3, newbas, lenb
      integer nbfnd, lenbad, lentrd, len3, ngrpm
      integer idena, nwroot, nosec
      logical noscf
      common /limy/ mp(10),lwordp,lvec,ied3,newbas,lenb,
     +              nbfnd,lenbad,lentrd,len3,ngrpm,
     +              idena,noscf,nwroot,nosec(mxroot)
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      dimension veca(*),vecb(*),eiga(*),eigb(*),ytype(2)
      data ytype/'-a-','-b-'/
      data itvec/3/
c
      newbas = num
c
c**** need to resolve this when in analysis mode
c
      minva = mouta
      minvb = moutb
c
      lenb = num*(num+1)/2
      l3 = num*num
      if (minva.le.0) then
         call caserr('invalid section specified for input vectors')
      end if
      call secget(minva,itvec,iblkv)
      call getqp(zcomm,zbuff,eiga,frocca,nbas,newbas,ncola,ieiga,
     +           idiffa,maxorb,iblkv)
      old = otran
      if(oroot) write (iwr,6030) ytype(1) , minva , ibl3d , yed(idaf)
      call rdedx(veca,l3,iblkv,idaf)
      if (.not.otran) call tdown(veca,ilifq,veca,ilifq,ncola)
      if (zscft.eq.'uhf') then
        if (minvb.le.0) then
         call caserr('invalid section specified for input vectors')
        else
           call secget(minvb,itvec,iblkv)
           call getqp(zcomm,zbuff,eigb,froccb,nbas,newbas,ncolb,
     +                ieigb,idiffb,maxorb,iblkv)
           if(oroot) write (iwr,6030) ytype(2) , minvb , ibl3d , 
     +                                yed(idaf)
           call rdedx(vecb,l3,iblkv,idaf)
           if (.not.otran) call tdown(vecb,ilifq,vecb,ilifq,
     +         ncolb)
        end if
      end if
      otran = .true.
      ibl7la= iposun(num8)
      otran = old
      return
 6030 format (/1x,a3,'vectors restored from section',i4,
     +        ' of dumpfile starting at block',i6,' of ',a4)
      end
      subroutine octpol(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, oroot
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
      character*8 scf, uhf, mp2
      character*4 keyefc
      logical out,otpout,ofield,ono
      common/efcpar/efcc(3,10),efcz(10),efclab(10),nefc
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      real*8 gx, gy, gz, ggx
      logical orgin
      common/gvalue/gx,gy,gz, ggx(3), orgin
      common/three/am(maxat),
     + thn(10),the(10),tht(10),thnefc(10),thtefc(10),
     + thnb(10),theb(10),thtb(10),thnefcb(10),thtefcb(10),
     + octn(10),octe(10),octt(10),octnefc(10),octtefc(10),
     + octnb(10),octeb(10),octtb(10),octnefcb(10),octtefcb(10)
      dimension x(*)
      dimension cm(3)
      character *3 xxx
      dimension xxx(10)
      data xxx/
     +   'xxx', 'yyy', 'zzz', 'xxy', 'xxz', 'yyx', 'yyz', 
     +   'zzx', 'zzy', 'xyz' /
      data scf,uhf /'scf     ','uhf     '/
      data mp2     /'mp2     '/
      data octo    /0.711688d+00/
      data tol     /1.0d-06/
      data zero    /0.0d+00/
      data two     /2.0d+00/
      data three   /3.0d+00/
      data four    /4.0d+00/
      data five    /5.0d+00/
      data keyefc  /' efc'/
c
      otpout=.true. 
      if(.not.otpout) return
      oroot = opg_root()
      otpout = otpout.and.oroot
      out=.false.       
      ono=.false.
      ofield = iefc.eq.keyefc
c
      if (oroot) then
       write(iwr,9991)
       write(iwr,9993)
      endif
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
c
c     ----- allocate memory block -----
c
      need = l2+l2+l2+l2+l2+l2+l2+l2+l2+l2+l2+l2+l3+l3
c
      init = igmem_alloc(need)
      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of core -----
c
      i10=  init
      i20=i10+l2
      i30=i20+l2
      i40=i30+l2
      i50=i40+l2
      i60=i50+l2
      i70=i60+l2
      i80=i70+l2
      i90=i80+l2
      i100=i90+l2
      i110=i100+l2
c     ----- memory block for density matrix -----
      i120=i110+l2
      i130=i120+l2
c     ----- memory block for -a- vectors  -----
      i140=i130+l3
c     ----- memory block for -b- vectors  -----
      last=i140+l3
c
      ncall=1
   10 continue
      call prpams(am,ncall,ncode,otpout)
      if(ncode.ne.0) go to 100
c
      sum=zero
      do 20 i=1,3
   20 cm(i)=zero
      do 30 iat=1,nat
      sum=sum+am(iat)
      do 30 i=1,3
   30 cm(i)=cm(i)+am(iat)*c(i,iat)
      if( abs(sum).lt.tol) go to 45
      do 40 i=1,3
   40 cm(i)=cm(i)/sum
   45 continue
      if(otpout) write(iwr,9998) cm(1),cm(2),cm(3)
c
c     use centre of mass, or origin nominated under gauge?
c
      if (orgin) then
       if (otpout) write(iwr,1234) (ggx(i),i=1,3)
       call dcopy(3,ggx(1),1,cm(1),1)
      else
       if (otpout) write(iwr,1235)
      endif
c
c     ----- calculate third moment integrals -----
c
      call octpint(x(i10),x(i20),x(i30),x(i40),x(i50),x(i60),
     1                x(i70),x(i80),x(i90),x(i100),
     2                cm,l1)
c
c     ----- electronic contribution to third moments -----
c
c     ----- density matrix for x(i110)=da,x(i120)=db(0?) -----

c     ----- calculate density matrix -----

      call prpdens(x(i110),x(i120),x(i130),x(i140),zscftp,ono,oroot)
      if (ono) return
c
c   xxx=1, yyy=2, zzz=3, xxy=4, xxz=5, yyx=6, yyz=7, 
c   zzx=8, zzy=9, xyz=10
c
      the(1) =-tracep(x(i110),x(i10),l1)
      the(2) =-tracep(x(i110),x(i20),l1)
      the(3) =-tracep(x(i110),x(i30),l1)
      the(4) =-tracep(x(i110),x(i40),l1)
      the(5) =-tracep(x(i110),x(i50),l1)
      the(6) =-tracep(x(i110),x(i60),l1)
      the(7) =-tracep(x(i110),x(i70),l1)
      the(8) =-tracep(x(i110),x(i80),l1)
      the(9) =-tracep(x(i110),x(i90),l1)
      the(10)=-tracep(x(i110),x(i100),l1)
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 50
       if(iuno.ne.0.and.iunosp.ne.0) go to 50
        the(1) =-tracep(x(i120),x(i10),l1)+the(1)
        the(2) =-tracep(x(i120),x(i20),l1)+the(2)
        the(3) =-tracep(x(i120),x(i30),l1)+the(3)
        the(4) =-tracep(x(i120),x(i40),l1)+the(4)
        the(5) =-tracep(x(i120),x(i50),l1)+the(5)
        the(6) =-tracep(x(i120),x(i60),l1)+the(6)
        the(7) =-tracep(x(i120),x(i70),l1)+the(7)
        the(8) =-tracep(x(i120),x(i80),l1)+the(8)
        the(9) =-tracep(x(i120),x(i90),l1)+the(9)
        the(10)=-tracep(x(i120),x(i100),l1)+the(10)
   50 continue
c
c     ----- nuclear contribution -----
c
      do i = 1, 10
       thn(i) = zero
       thnefc(i) = zero
      enddo
c
      do 60 iat=1,nat
      thn(1)=thn(1)+czan(iat)*
     +       (c(1,iat)-cm(1))*(c(1,iat)-cm(1))
     +                       *(c(1,iat)-cm(1))
      thn(2)=thn(2)+czan(iat)*
     +       (c(2,iat)-cm(2))*(c(2,iat)-cm(2))
     +                       *(c(2,iat)-cm(2))
      thn(3)=thn(3)+czan(iat)*
     +       (c(3,iat)-cm(3))*(c(3,iat)-cm(3))
     +                       *(c(3,iat)-cm(3))
      thn(4)=thn(4)+czan(iat)*
     +       (c(1,iat)-cm(1))*(c(1,iat)-cm(1))
     +                       *(c(2,iat)-cm(2))
      thn(5)=thn(5)+czan(iat)*
     +       (c(1,iat)-cm(1))*(c(1,iat)-cm(1))
     +                       *(c(3,iat)-cm(3))
      thn(6)=thn(6)+czan(iat)*
     +       (c(2,iat)-cm(2))*(c(2,iat)-cm(2))
     +                       *(c(1,iat)-cm(1))
      thn(7)=thn(7)+czan(iat)*
     +       (c(2,iat)-cm(2))*(c(2,iat)-cm(2))
     +                       *(c(3,iat)-cm(3))
      thn(8)=thn(8)+czan(iat)*
     +       (c(3,iat)-cm(3))*(c(3,iat)-cm(3))
     +                       *(c(1,iat)-cm(1))
      thn(9)=thn(9)+czan(iat)*
     +       (c(3,iat)-cm(3))*(c(3,iat)-cm(3))
     +                       *(c(2,iat)-cm(2))
      thn(10)=thn(10)+czan(iat)*
     +       (c(1,iat)-cm(1))*(c(2,iat)-cm(2))
     +                       *(c(3,iat)-cm(3))
   60 continue
c
c     ----- -efc- contribution -----
c
      if (ofield) then
      do iat=1,nefc
      thnefc(1)=thnefc(1)+efcz(iat)*
     +                   (efcc(1,iat)-cm(1))*(efcc(1,iat)-cm(1))
     +                  *(efcc(1,iat)-cm(1))
      thnefc(2)=thnefc(2)+efcz(iat)*
     +                   (efcc(2,iat)-cm(2))*(efcc(2,iat)-cm(2))
     +                  *(efcc(2,iat)-cm(2))
      thnefc(3)=thnefc(3)+efcz(iat)*
     +                   (efcc(3,iat)-cm(3))*(efcc(3,iat)-cm(3))
     +                  *(efcc(3,iat)-cm(3))
      thnefc(4)=thnefc(4)+efcz(iat)*
     +                   (efcc(1,iat)-cm(1))*(efcc(1,iat)-cm(1))
     +                  *(efcc(2,iat)-cm(2))
      thnefc(5)=thnefc(5)+efcz(iat)*
     +                   (efcc(1,iat)-cm(1))*(efcc(1,iat)-cm(1))
     +                  *(efcc(3,iat)-cm(3))
      thnefc(6)=thnefc(6)+efcz(iat)*
     +                   (efcc(2,iat)-cm(2))*(efcc(2,iat)-cm(2))
     +                  *(efcc(1,iat)-cm(1))
      thnefc(7)=thnefc(7)+efcz(iat)*
     +                   (efcc(2,iat)-cm(2))*(efcc(2,iat)-cm(2))
     +                  *(efcc(3,iat)-cm(3))
      thnefc(8)=thnefc(8)+efcz(iat)*
     +                   (efcc(3,iat)-cm(3))*(efcc(3,iat)-cm(3))
     +                  *(efcc(1,iat)-cm(1))
      thnefc(9)=thnefc(9)+efcz(iat)*
     +                   (efcc(3,iat)-cm(3))*(efcc(3,iat)-cm(3))
     +                  *(efcc(2,iat)-cm(2))
      thnefc(10)=thnefc(10)+efcz(iat)*
     +                   (efcc(1,iat)-cm(1))*(efcc(2,iat)-cm(2))
     +                  *(efcc(3,iat)-cm(3))
      enddo
      endif
c
      do i = 1,10
       tht(i) = the(i) + thn(i)
       thtefc(i) = tht(i) + thnefc(i)
c       convert to 10**(-34) esu*cm**3
       theb(i) = the(i) * octo
       thnb(i) = thn(i) * octo
       thtb(i) = tht(i) * octo
       thnefcb(i) = thnefc(i) * octo
       thtefcb(i) = thtefc(i) * octo
      enddo
c
      if (oroot) then
        write(iwr,1235) 
        write(iwr,6000)
        if (ofield) then
         write(iwr,6005)
         do i = 1 , 10
          write (iwr,6010)xxx(i),thn(i),thnefc(i),the(i),
     +                           thtefc(i)
         enddo
         write(iwr,6090)
         do i = 1 , 10
          write (iwr,6010)xxx(i),thnb(i),thnefcb(i),theb(i),
     +                          thtefcb(i)
         enddo
c
        else
         write (iwr,6030)
         do i = 1 , 10
           write (iwr,6010)xxx(i),thn(i) ,the(i) ,tht(i)
         enddo
         write(iwr,6090)
         do i = 1 , 10
           write (iwr,6010)xxx(i),thnb(i) ,theb(i) ,thtb(i)
         enddo
        endif
      endif
c
      octn(1) =      thn(1)-three*(thn(6)+thn(8))/two
      octn(2) =      thn(2)-three*(thn(4)+thn(9))/two
      octn(3) =      thn(3)-three*(thn(5)+thn(7))/two
      octn(4) =(four*thn(4)-thn(2)-thn(9))/two
      octn(5) =(four*thn(5)-thn(7)-thn(3))/two
      octn(6) =(four*thn(6)-thn(1)-thn(8))/two
      octn(7) =(four*thn(7)-thn(5)-thn(3))/two
      octn(8) =(four*thn(8)-thn(1)-thn(6))/two
      octn(9) =(four*thn(9)-thn(4)-thn(2))/two
      octn(10)= five*thn(10)/two
c
      octe(1) =      the(1)-three*(the(6)+the(8))/two
      octe(2) =      the(2)-three*(the(4)+the(9))/two
      octe(3) =      the(3)-three*(the(5)+the(7))/two
      octe(4) =(four*the(4)-the(2)-the(9))/two
      octe(5) =(four*the(5)-the(7)-the(3))/two
      octe(6) =(four*the(6)-the(1)-the(8))/two
      octe(7) =(four*the(7)-the(5)-the(3))/two
      octe(8) =(four*the(8)-the(1)-the(6))/two
      octe(9) =(four*the(9)-the(4)-the(2))/two
      octe(10)= five*the(10)/two
c
      octnefc(1) =      thnefc(1)-three*(thnefc(6)+thnefc(8))/two
      octnefc(2) =      thnefc(2)-three*(thnefc(4)+thnefc(9))/two
      octnefc(3) =      thnefc(3)-three*(thnefc(5)+thnefc(7))/two
      octnefc(4) =(four*thnefc(4)-thnefc(2)-thnefc(9))/two
      octnefc(5) =(four*thnefc(5)-thnefc(7)-thnefc(3))/two
      octnefc(6) =(four*thnefc(6)-thnefc(1)-thnefc(8))/two
      octnefc(7) =(four*thnefc(7)-thnefc(5)-thnefc(3))/two
      octnefc(8) =(four*thnefc(8)-thnefc(1)-thnefc(6))/two
      octnefc(9) =(four*thnefc(9)-thnefc(4)-thnefc(2))/two
      octnefc(10)= five*thnefc(10)/two
c
      do i = 1,10
       octt(i) = octe(i) + octn(i)
       octtefc(i) = octt(i) + octnefc(i)
c       convert to 10**(-34) esu*cm**3
       octeb(i) = octe(i) * octo
       octnb(i) = octn(i) * octo
       octtb(i) = octt(i) * octo
       octnefcb(i) = octnefc(i) * octo
       octtefcb(i) = octtefc(i) * octo
      enddo
c
      if (oroot) then
       write(iwr,1235) 
       write(iwr,7000)
        if (ofield) then
         write(iwr,6005)
         do i = 1 , 10
          write (iwr,6010)xxx(i),octn(i),octnefc(i),octe(i),
     +                           octtefc(i)
         enddo
         write(iwr,7090)
         do i = 1 , 10
          write (iwr,6010)xxx(i),octnb(i),octnefcb(i),octeb(i),
     +                           octtefcb(i)
         enddo
c
        else
         write (iwr,6030)
         do i = 1 , 10
           write (iwr,6010)xxx(i),octn(i) ,octe(i) ,octt(i)
         enddo
         write(iwr,7090)
         do i = 1 , 10
           write (iwr,6010)xxx(i),octnb(i) ,octeb(i) ,octtb(i)
         enddo
        endif
      endif
c
      ncall=ncall+1
c     go to 10
c
  100 continue
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
c
 9991 format(/10x,26('-')
     +       /10x,'third and octupole moments',
     +       /10x,26('-')/)
 9993 format(/' 1 a.u. = 0.711688 10**(-34) esu*cm**3 ')
 6000 format(/2x,
     + 'third moments in atomic units (w.r.t. above origin)'/2x,
     + '---------------------------------------------------')
 6090 format(/2x,
     +'third moments in 10**(-34) esu*cm**3 (w.r.t. above origin)'/
     + 2x,
     +'----------------------------------------------------------'/)
 6005 format(/13x,'nuclear',11x,'-efc-',6x,'electronic',11x,'total'/)
 6030 format(/13x,'nuclear',6x,'electronic',11x,'total'/)
 6010 format(1x,a3,4f16.7)
 7000 format(/2x,
     + 'octupole moments in atomic units (w.r.t. above origin)'/2x,
     + '------------------------------------------------------')
 7090 format(/2x,'octupole moments in 10**(-34) esu*cm**3 ',
     + '(w.r.t. above origin)'/
     +        2x,'----------------------------------------',
     + '---------------------'/)
 9998 format(/2x,' centre of mass',
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 1234 format (/1x,'origin for property evaluation:',
     + ' x = ',f15.7,' y = ',f15.7,' z = ',f15.7)
 1235 format (/1x,'origin for property evaluation is centre of mass')
 9996 format(' not enough core in -otpole- . stop ')
      end
      subroutine octpint(xxxs,yyys,zzzs,xxys,xxzs,yyxs,yyzs,
     1                  zzxs,zzys,xyzs,cm,l1)
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      logical iandj
      logical norm,double
      logical out
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
      common/xyzoctp/xint,yint,zint,xintx,yinty,zintz,
     1              xintxx,yintyy,zintzz,xinxxx,yinyyy,zinzzz,
     2              t,xc,yc,zc,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      dimension xxxs(*),yyys(*),zzzs(*),xxys(*),xxzs(*),yyxs(*),yyzs(*)
      dimension zzxs(*),zzys(*),xyzs(*)
      common /junk3/ sxxx(225),syyy(225),szzz(225),
     +               sxxy(225),sxxz(225),syyx(225),syyz(225),
     +               szzx(225),szzy(225),sxyz(225),
     +               dij(225),
     +               xin(5,5),   yin(5,5),   zin(5,5),
     +               xxin(5,5),  yyin(5,5),  zzin(5,5),
     +               xxxin(5,5), yyyin(5,5), zzzin(5,5),
     +               xxxxin(5,5),yyyyin(5,5),zzzzin(5,5)
      dimension ijx(35),ijy(35),ijz(35)
      dimension cm(3)
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data rln10 /2.30258d+00/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
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
      out=nprint.eq.6
      norm=normf.ne.1.or.normp.ne.1
c
c
      xc=cm(1)
      yc=cm(2)
      zc=cm(3)
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
      iandj=ii.eq.jj
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
c
c     ----- prepare indices for pairs of (i,j) functions -----
c
      ij=0
      do 100 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=ij+1
      sxxx(ij)=zero
      syyy(ij)=zero
      szzz(ij)=zero
      sxxy(ij)=zero
      sxxz(ij)=zero
      syyx(ij)=zero
      syyz(ij)=zero
      szzx(ij)=zero
      szzy(ij)=zero
      sxyz(ij)=zero
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
      cgi=cg(ig)
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
      cgj=cg(jg)
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
  180 dum1=cgi*fac
      go to 220
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
  310 dum2=dum1*cgj
      if(double) dum2=dum2+dum2
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
c     ----- octupole moment integrals -----
c
      t = sqrt(aa)
      t1=one/t
      x0=ax
      y0=ay
      z0=az
      do 370 j=1,ljt
      nj=j
      do 370 i=1,lit
      ni=i
      call octpxyz
         xin(i,j)=xint*t1
         yin(i,j)=yint*t1
         zin(i,j)=zint*t1
        xxin(i,j)=xintx*t1
        yyin(i,j)=yinty*t1
        zzin(i,j)=zintz*t1
       xxxin(i,j)=xintxx*t1
       yyyin(i,j)=yintyy*t1
       zzzin(i,j)=zintzz*t1
      xxxxin(i,j)=xinxxx*t1
      yyyyin(i,j)=yinyyy*t1
      zzzzin(i,j)=zinzzz*t1
  370 continue
c
      ij=0
      do 390 i=mini,maxi
      ix=ijx(i)
      iy=ijy(i)
      iz=ijz(i)
      jmax=maxj
      if(iandj) jmax=i
      do 390 j=minj,jmax
      jx=ijx(j)
      jy=ijy(j)
      jz=ijz(j)
      ij=ij+1
      sxxx(ij)=sxxx(ij)+
     1         dij(ij)*(xxxxin(ix,jx)*   yin(iy,jy)*   zin(iz,jz))
      syyy(ij)=syyy(ij)+
     1         dij(ij)*(   xin(ix,jx)*yyyyin(iy,jy)*   zin(iz,jz))
      szzz(ij)=szzz(ij)+
     1         dij(ij)*(   xin(ix,jx)*   yin(iy,jy)*zzzzin(iz,jz))
      sxxy(ij)=sxxy(ij)+
     1         dij(ij)*( xxxin(ix,jx)*  yyin(iy,jy)*   zin(iz,jz))
      sxxz(ij)=sxxz(ij)+
     1         dij(ij)*( xxxin(ix,jx)*   yin(iy,jy)*  zzin(iz,jz))
      syyx(ij)=syyx(ij)+
     1         dij(ij)*(  xxin(ix,jx)* yyyin(iy,jy)*   zin(iz,jz))
      syyz(ij)=syyz(ij)+
     1         dij(ij)*(   xin(ix,jx)* yyyin(iy,jy)*  zzin(iz,jz))
      szzx(ij)=szzx(ij)+
     1         dij(ij)*(  xxin(ix,jx)*   yin(iy,jy)* zzzin(iz,jz))
      szzy(ij)=szzy(ij)+
     1         dij(ij)*(   xin(ix,jx)*  yyin(iy,jy)* zzzin(iz,jz))
      sxyz(ij)=sxyz(ij)+
     1         dij(ij)*(  xxin(ix,jx)*  yyin(iy,jy)*  zzin(iz,jz))
  390 continue
 6000 continue
 7000 continue
c
c     ----- set up octupole moment matrices -----
c
      ij=0
      do 7500 i=mini,maxi
      jmax=maxj
      if(iandj) jmax=i
      do 7500 j=minj,jmax
      ij=ij+1
      nn=iky(loci+i)+(locj+j)
      xxxs(nn)=sxxx(ij)
      yyys(nn)=syyy(ij)
      zzzs(nn)=szzz(ij)
      xxys(nn)=sxxy(ij)
      xxzs(nn)=sxxz(ij)
      yyxs(nn)=syyx(ij)
      yyzs(nn)=syyz(ij)
      zzxs(nn)=szzx(ij)
      zzys(nn)=szzy(ij)
      xyzs(nn)=sxyz(ij)
 7500 continue
c
 8000 continue
 9000 continue
c
c
      if(out) then
         write(iwr,9999)
         call writel(xxxs,l1)
         write(iwr,9998)
         call writel(yyys,l1)
         write(iwr,9997)
         call writel(zzzs,l1)
         write(iwr,9996)
         call writel(xxys,l1)
         write(iwr,9995)
         call writel(xxzs,l1)
         write(iwr,9994)
         call writel(yyxs,l1)
         write(iwr,9993)
         call writel(yyzs,l1)
         write(iwr,9992)
         call writel(zzxs,l1)
         write(iwr,9991)
         call writel(zzys,l1)
         write(iwr,9990)
         call writel(xyzs,l1)
      endif
c
      return
 9999 format(/10x,26('*')/10x,'xxx-third moment integrals'/
     +         10x,26('*'))
 9998 format(/10x,26('*')/10x,'yyy-third moment integrals'/
     +         10x,26('*'))
 9997 format(/10x,26('*')/10x,'zzz-third moment integrals'/
     +         10x,26('*'))
 9996 format(/10x,26('*')/10x,'xxy-third moment integrals'/
     +         10x,26('*'))
 9995 format(/10x,26('*')/10x,'xxz-third moment integrals'/
     +         10x,26('*'))
 9994 format(/10x,26('*')/10x,'yyx-third moment integrals'/
     +         10x,26('*'))
 9993 format(/10x,26('*')/10x,'yyz-third moment integrals'/
     +         10x,26('*'))
 9992 format(/10x,26('*')/10x,'zzx-third moment integrals'/
     +         10x,26('*'))
 9991 format(/10x,26('*')/10x,'zzy-third moment integrals'/
     +         10x,26('*'))
 9990 format(/10x,26('*')/10x,'xyz-third moment integrals'/
     +         10x,26('*'))
      end
      subroutine octpxyz
      implicit real*8 (a-h,o-z), integer   (i-n)
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
      common/xyzoctp/xint,yint,zint,xintx,yinty,zintz,
     1              xintxx,yintyy,zintzz,xinxxx,yinyyy,zinzzz,
     2              t,xc,yc,zc,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(7),max(7)
      data min /1,2,4, 7,11,16,22/
      data max /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      xintx=zero
      yinty=zero
      zintz=zero
      xintxx=zero
      yintyy=zero
      zintzz=zero
      xinxxx=zero
      yinyyy=zero
      zinzzz=zero
      npts=(ni+nj-2+3)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
      xintx=xintx+px*(ptx-xc)
      yinty=yinty+py*(pty-yc)
      zintz=zintz+pz*(ptz-zc)
      xintxx=xintxx+px*(ptx-xc)*(ptx-xc)
      yintyy=yintyy+py*(pty-yc)*(pty-yc)
      zintzz=zintzz+pz*(ptz-zc)*(ptz-zc)
      xinxxx=xinxxx+px*(ptx-xc)*(ptx-xc)*(ptx-xc)
      yinyyy=yinyyy+py*(pty-yc)*(pty-yc)*(pty-yc)
      zinzzz=zinzzz+pz*(ptz-zc)*(ptz-zc)*(ptz-zc)
   16 continue
      return
      end
      subroutine mulken2(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, oroot
      logical o1e, uno, ono
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
      logical   opass
      character*8 scf, uhf, mp2, prop
      character*8 elec, alpha, beta, all, spin
      logical last
      logical out
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
      common/junk2/xyzmap(3,100),zznuc(maxat),nuc(maxat)
      common/junk/limlow(maxat),limsup(maxat)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z
      integer ibl3rs, ibl3qa, ibl3pa, ibl3ea, ibl3qb, ibl3pb
      integer ibl3eb, ibl3g, ibl3hs, ionsec, ibl3op, ions2
      integer isecqm
      common /dump3/ ibl3s, ibl3t, ibl3f, ibl3x, ibl3y, ibl3z,
     +               ibl3rs,ibl3qa,ibl3pa,ibl3ea,ibl3qb,ibl3pb,
     +               ibl3eb,ibl3g,ibl3hs,ionsec,ibl3op,ions2,
     +               isecqm
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common/blkcore/corev(512),array(10)
      dimension x(*)
      dimension o1e(6)
      data alpha,beta,all /'   alpha','    beta','     all'/
      data spin           /'    spin' /
      data scf,uhf        /'scf     ','uhf     '/
      data mp2            /'mp2     '/
      data prop           /'propty  '/
c
      oroot = opg_root()
      out=.false.
      ono=.false.
      if (oroot) then
       write(iwr,9991)
      endif
      ispher = -1
c
      l1= num
      l2=(num*(num+1))/2
      l3= num*num
c
c     ----- allocate memory block -----
c
      need = l3+l1+l2  + 
     &       (nat*(nat+1))/2+l2+l1+nat  +
     &       l2+l1+nat +
     &       l3+l1     +
     &       l3+nat*nat+nat+nat+l3+nat+nat  +
     &       l3+l3

      init = igmem_alloc(need)
c
      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of core -----
c
      i00=  init
      i01=i00+l3
      i10=i01+l1
      i20=i10+l2
c
      i30=i20
      i40=i30+(nat*(nat+1))/2
      i50=i40+l2
      i60=i50+l1
      i70=i60+nat
c
      i41=i70
      i51=i41+l2
      i61=i51+l1
      i71=i61+nat
c
      i42=i71+l3
      i52=i42+l1
c
      i53=i52+l3
      i63=i53+nat*nat
      i73=i63+nat
      i83=i73+nat
      i93=i83+l3
      i94=i93+nat
      i95=i94+nat
c     ----- memory block for -a- vectors -----
      i96=i95+l3
c     ----- memory block for -b- vectors -----
      lastw=i96+l3
c
      do iat=1,nat
         zznuc(iat)=nuc(iat)
      enddo
c
c     ----- determine number of basis functions per atom -----
c
      call atmbfn
c
c     -----         read in overlap matrix              -----
c
c     -  s- at x(i10)
c
      o1e(1) = .true.
      do loop = 2,6
       o1e(loop) = .false.
      enddo
      call getmat(x(i10),x(i10),x(i10),x(i10),x(i10),x(i10),
     +            array,l1,o1e,ionsec)
      if(out.and.oroot) then
         write(iwr,9989)
         call prtri(x(i10),l1)
      endif
c
c     are spinfree and spin nos available for uhf calculation
c
      uno = iuno.ne.0.and.iunosp.ne.0
c
      if(zruntp.eq.prop) go to 1000
c
c     -----            read density matrices           -----
c
c               -da- at x(i40),  -db- at x(i41)
c
c     -----  or for uhf calculations with natorb,
c           spinfree at x(i40), spin nos at x(i41)
c
c     ----- calculate density matrix -----

      call prpdens(x(i40),x(i41),x(i95),x(i96),zscftp,ono,oroot)
      if (ono) return
      if(zscftp.eq.uhf) then
         wfntyp=scf
      endif
      if(out) then
         write(iwr,*) 'wfntyp = ',wfntyp
         write(iwr,*) 'scftyp = ',zscftp
      endif
c
      ipass=1
      iden=i40
      igoc=i50
      igac=i60
      last=.false.
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)) then
       elec = all
      else if(uno) then
       elec = all
      else
       elec = alpha
      endif
  100 continue
c
      opass = .not.(ipass.eq.2.and.uno)
c
c     ----- do a mulliken population analysis -----
c              calculate overlap population
c
      call vmul(x(iden),1,x(i10),1,x(iden),1,l2)
  200 if (oroot.and.opass) write(iwr,9999) elec
      if(out.and.oroot) then
         write(iwr,9988)
         call prtri (x(iden),l1)
      endif
c
c     ----- calculate total gross population in ao's ----
c
      call grossc(x(iden),x(igoc),l1)
      if (oroot.and.opass) then
         write(iwr,9998)
         write(iwr,9997) (i, zbflab(i),x(i-1+igoc),i=1,l1)
      endif
c
      if (oroot.and.opass) write(iwr,9996)
c
c     ----- compress from -ao- to shells -----
c
      call shlpop(x(igoc),ispher)
c
c     ----- compress from -ao- to atoms -----
c
      call atpop(x(iden),x(i30),nat)
      if(out) then
         call prtri(x(i30),nat)
      endif
c
c     ----- calculate total gross population on atoms -----
c
      call grossc(x(i30),x(igac),nat)
      if (oroot.and.opass) then
         write(iwr,9995)
         write(iwr,9992) (i,zaname(i),zznuc(i),czan(i),
     1                x(i-1+igac),i=1,nat)
      endif
c
      if((wfntyp.ne.scf.or.zscftp.ne.uhf).and.
     1   (wfntyp.ne.mp2.or.zscftp.ne.uhf)     ) go to 400
      if(last) go to 400
      if(ipass.eq.2) go to 300
      ipass=2
      iden=i41
      igoc=i51
      igac=i61
      elec=beta
      if(uno) elec = spin
      go to 100
  300 continue
c
c     ----- calculate orbital and atomic spin densities -----
c
      if (uno) then
        itmp1 = i51
        itmp2 = i61
      else
        itmp1 = i50
        itmp2 = i60
        do i=1,l1
           x(i-1+i50)=x(i-1+i50)-x(i-1+i51)
        enddo
        do i=1,nat
           x(i-1+i60)=x(i-1+i60)-x(i-1+i61)
        enddo
      endif
      if (oroot) then
         write(iwr,9994)
         write(iwr,9997) (i, zbflab(i),x(i-1+itmp1),i=1,num)
         write(iwr,9993)
         write(iwr,9992) (i,zaname(i),zznuc(i),czan(i),
     1                x(i-1+itmp2),i=1,nat)
      endif
c
c     ----- do all electrons -----
c
      if (uno) go to 400
      ipass=1
      last=.true.
      elec=all
      do i=1,l2
         x(i-1+i40)=x(i-1+i40)+x(i-1+i41)
      enddo
      iden=i40
      igoc=i50
      igac=i60
      go to 200
  400 continue
c
c     ----- bond index and valency analysis -----
c
      call vclr(x(i41),1,l2)
c
c     ----- calculate density matrix -----

      call prpdens(x(i40),x(i41),x(i95),x(i96),zscftp,ono,oroot)
      if(ono)return
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
       if(iuno.eq.0.or.iunosp.eq.0) then
         do i=1,l2
            duma=x(i+i40-1)
            dumb=x(i+i41-1)
            x(i+i40-1)=duma+dumb
            x(i+i41-1)=duma-dumb
         enddo
       endif
      endif
      call apspop(x(i10),x(i40),x(i41),x(i51),x(i83),
     1            x(i53),x(i63),x(i73),x(i93),x(i94),iky,nat,num)
c
 1000 continue
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return
 9991 format(/10x,28('-')
     1       /10x,'mulliken population analysis',
     2       /10x,28('-')/)
 9999 format(//10x,28('-')/10x,'mulliken population analysis',
     1 10x,a8,' electrons'/10x,28('-'))
 9998 format(/10x,'----- total gross population in ao -----'/)
 9997 format(4(i5,1x,a10,f10.4))
 9996 format(/10x,'----- condensed to atoms -----'/)
 9995 format(/10x,'----- total gross population on atoms ----'/)
 9994 format(/10x,'----- spin population in ao -----'/)
 9993 format(/10x,'----- atomic spin population -----'/)
 9992 format(10x,i5,2x,a8,2x,f6.1,f6.1,f12.5)
 9989 format(/10x,'----- overlap matrix in -mulken- -----')
 9988 format(/10x,'----- overlap population -----')
      end
c
      subroutine atmbfn
      implicit real*8 (a-h,o-z), integer   (i-n)
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
      common/junk/limlow(maxat),limsup(maxat)
c
      do i=1,nat
         limlow(i)=num+1
         limsup(i)=0
         do j=1,nshell
            if(katom(j).eq.i) then
               limlow(i)=min(limlow(i),kloc(j)                )
               limsup(i)=max(limsup(i),kloc(j)+kmax(j)-kmin(j))
            endif
         enddo
      enddo
c
      return
      end
      subroutine apspop(s,d,ds,p,ps,b,f,v,sb,fv,ia,nat,num)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root,oroot
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
      common/junk/limlow(maxat),limsup(maxat)
      dimension s(*),d(*),p(num,*),b(nat,1),f(*),v(*),ia(*)
      dimension ds(*),ps(num,1),sb(*),fv(*)
      dimension ii(4),jj(4)
      data zero /0.0d+00/
      data two  /2.0d+00/
      data btol /0.1d+00/
c
      oroot = opg_root()
c
c     ----- closed and open shell systems population analysis -----
c
c     original reference for ndo systems:
c     d.r. armstrong, p.g. perkins and j.j.p. stewart
c     j.c.s. dalton trans. 838 (1973)
c     extension  first principles: i. mayer, chem. phys. lett.
c     97,270(1983).
c     number of free electrons: o.g.stradella, h.o. villar, e.a.castro
c     theoret. chim. acta(berl.) 70, 67(1986)
c
c       d * s = p (non - hermitean)
c
      do i=1,num
         do j=1,num
            dums=zero
            dum =zero
            do k=1,num
               ik=ia(max(i,k))+min(i,k)
               kj=ia(max(k,j))+min(k,j)
               dums=dums+ds(ik)*s(kj)
               dum =dum + d(ik)*s(kj)
            enddo
            ps(i,j)=dums
             p(i,j)=dum
         enddo
      enddo
c
c     ----- -b-, valency , and free electron -----
c
      do iat=1,nat
         lli=limlow(iat)
         lui=limsup(iat)
         do jat=1,nat
            llj=limlow(jat)
            luj=limsup(jat)
            dum=zero
            if(lui.ne.0.and.luj.ne.0) then
               do l = llj,luj
                  do k = lli,lui
                     dum=dum+p(l,k)*p(k,l)+ps(l,k)*ps(k,l)
                  enddo
               enddo
            endif
            b(iat,jat)=dum
         enddo
         dums=zero
         dum =zero
         do l = lli,lui
            dum=dum+two*p(l,l)
            do m = lli,lui
               dums=dums-ps(l,m)*ps(m,l)
               dum =dum - p(l,m)* p(m,l)
            enddo
         enddo
         v(iat)=dum
         dum=zero
         do jat=1,nat
            if(iat.ne.jat) then
               dum=dum-b(iat,jat)
            endif
         enddo
         f(iat)=dum+(v(iat)+dums+b(iat,iat))/two
         b(iat,iat)=zero
      enddo
c
      if (oroot) write(iwr,9999)
c
c     ----- print all bond indices -----
c
      mmax=0
   10 mmin=mmax+1
      mmax=mmax+8
      if(mmax.gt.nat) mmax=nat
      if (oroot) write(iwr,9995)
      do i =1,nat
       if (oroot) write(iwr,9998) (i,j,b(i,j),j=mmin,mmax)
      enddo
      if(mmax.lt.nat) go to 10
c
c     ----- now print only the largest one -----
c
      if (oroot) write(iwr,9993)
      kk=0
      do i=1,nat
         do j=i+1,nat
            if(b(i,j).gt.btol) then
               kk = kk+1
               ii(kk) = i
               jj(kk) = j
               if(kk.eq.4) then

                 if (oroot)  write(iwr,9994)
     1              (ii(k),zaname(ii(k)),jj(k),zaname(jj(k)),
     1                            b(ii(k),jj(k)),k=1,4)

                  kk = 0
               endif
            endif
         enddo
         if(kk.ne.0) then

            if (oroot) write(iwr,9994)
     1          (ii(k),zaname(ii(k)),jj(k),zaname(jj(k)),
     1                      b(ii(k),jj(k)),k=1,kk)

         endif
         kk=0
      enddo
c
c     ----- sum bond indices for a given atom -----
c
      do jat=1,nat
         dum=zero
         do iat=1,nat
            dum=dum+b(iat,jat)
         enddo
         b(jat,jat)=dum
      enddo
c
c     ----- v(a) = n(a) + sum of b(a,b) is the mulliken charge -----
c           spin(a) = sum of b(a,b)-v(a)
c
      do iat=1,nat
         sb(iat)=b(iat,iat)+f(iat)
         fv(iat)=v(iat)-b(iat,iat)
      enddo
c
c     ----- print analysis -----
c
      if (oroot) write(iwr,9997)
      do i =1,nat
         if (oroot)
     & write(iwr,9996)i,zaname(i),v(i),f(i),b(i,i),sb(i),fv(i)
      enddo
c
      return
 9999 format(/20x,'----- Bond Indices -----'//
     +        10x,'O.G.Stradella, H.O. Villar, and E.A. Castro'/
     +        10x,'Theoret. Chim. Acta (Berl.) 70, 67 (1986)'/)
 9998 format(8(1x,i2,'-',i2,f10.5))
 9997 format(/70x,'  Free electrons',5x,'        Valency'/
     1         36x,'  Number of   ',3x,'   Sum of    ',
     2          6x,'+ bond indices ',5x,'   - bond indices    '/
     3         21x,'   Valency',5x,'Free electrons',4x,'Bond indices',
     4          5x,'=Mulliken charge',5x,'= Net spin population')
 9996 format(1x,i5,2x,a8,2x,3x,f10.5,5x,f10.5,7x,f10.5,8x,f10.5,
     1      14x,f10.5)
 9995 format(/)
 9994 format(5(i3,1x,a4,'-',i3,1x,a4,f10.5))
 9993 format(/1x,'large bond indices'/1x,18('-'))
      end
c
      subroutine shlpop(aopop,ispher)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root,oroot
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
      parameter (zero=0.0d+00)
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
      dimension aopop(*)
      dimension shl(5)
c
      oroot = opg_root()
      if (oroot) write(iwr,9999)
      if(ispher.lt.0) then
c
c     ----- in -ao- basis -----
c
         do iat=1,nat
            do ishl=1,5
               shl(ishl)=zero
            enddo
            do ish=1,nshell
               if(katom(ish).eq.iat) then
                  lit=ktype(ish)
                  min=kmin(ish)
                  max=kmax(ish)
                  loc=kloc(ish)-min
                  if(lit.eq.2) then
                     if(min.eq.1) then
                        shl(1)=shl(1)+aopop( 1+loc)
                        min=min+1
                     endif
                  endif
                  dum=zero
                  do i=min,max
                     dum=dum+aopop(i+loc)
                  enddo
                  shl(lit)=shl(lit)+dum
               endif
            enddo
            sum=shl(1)+shl(2)+shl(3)+shl(4)+shl(5)

            if (oroot) write(iwr,9998)
     &         iat,zaname(iat),(shl(i),i=1,5),sum

         enddo
      else
c
c     ----- in -sph.harm.- basis -----
c
         do iat=1,nat
            do ishl=1,5
               shl(ishl)=zero
            enddo
            do ish=1,nshell
               if(katom(ish).eq.iat) then
                  lit=ktype(ish)
                  min=kmin(ish)
                  max=kmax(ish)
                  loc=kloc(ish)-min
                  if(lit.eq.2) then
                     if(min.eq.1) then
                        shl(1)=shl(1)+aopop( 1+loc)
                        min=min+1
                     endif
                  elseif(lit.eq.3) then
                     shl(1)=shl(1)+aopop(10+loc)
                     max=max-1
                  elseif(lit.eq.4) then
                     shl(2)=shl(2)+aopop(18+loc)
     1                            +aopop(19+loc)+aopop(20+loc)
                     max=max-3
                  elseif(lit.eq.5) then
                     shl(1)=shl(1)+aopop(30+loc)
                     shl(3)=shl(3)+aopop(31+loc)+aopop(32+loc)
     1                            +aopop(33+loc)+aopop(34+loc)
     2                            +aopop(35+loc)
                     max=max-6
                  endif
                  dum=zero
                  do i=min,max
                     dum=dum+aopop(i+loc)
                  enddo
                  shl(lit)=shl(lit)+dum
               endif
            enddo
            sum=shl(1)+shl(2)+shl(3)+shl(4)+shl(5)

            if (oroot) write(iwr,9998) 
     &         iat,zaname(iat),(shl(i),i=1,5),sum

         enddo
      endif
c
      return
 9999 format(/1x,'total s,p,d,... shell population'/1x,32('-')
     1       /4x,'     atom      ','    s     ','    p     ',
     2    '    d     ','    f     ','    g     ','  total   ',
     3       /1x,76('-'))
 9998 format(1x,i4,1x,a8,2x,6f10.5)
      end
      subroutine spind2(x)
      implicit real*8 (a-h,o-z), integer   (i-n)
      logical opg_root, ono, oroot
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
      character*8 scf, uhf, mp2
      logical some,out
      logical nospin
      character*4 iefc
      character*8 wfntyp
      common/prop3c/wfntyp, iefc
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
      real*8 cicoef, f, alpha, beta
      integer no, nco, nseto, npair, ncores, ibm, nset, nopen
      integer  nope, noe
      logical old
      common/scfwfn/cicoef(2,12),f(25),alpha(325),beta(325),
     + no(10),nco,nseto,npair,ncores,ibm,old,nset,nopen,nope,noe(10)
c
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
c
c...   for UHF natorbs
c
      integer iuno, iunopp, iunosp, iunspp
      logical  oanil
      common/unocas/iuno,iunopp,iunosp,iunspp,oanil
c
      common/blkorbs/ deig(maxorb),docc(maxorb+1),ipad(6)
      common/bufb/eiga(maxorb),frocca(maxorb),eigb(maxorb),
     +            froccb(maxorb)
c
      real*8 cconv, atwt, cspin
      integer lin, kin, nin, ngmx, nfmx, nomx, ncmx, ntmx, ntcol, npmx
      integer icon, ncom, maxnt,itps, maxcon
c
      common /datprp/ cconv(11),atwt(30),cspin(30),
     + lin(84),kin(84),nin(84),
     + ngmx,nfmx,nomx,ncmx,ntmx,ntcol,npmx,
     + icon(24),ncom(22),maxnt,itps,maxcon
c
      dimension x(*)
      dimension a(6),vec(3,3),eig(3)
c     dimension isotop(18),gn(18)
      data one    /1.0d+00/
      data two    /2.0d+00/
      data three  /3.0d+00/
      data four   /4.0d+00/
      data uhf    /'uhf'/
      data mp2    /'mp2'/
      data scf    /'scf'/
      data angle  /180.d+00/ 
      data ge     /2.002319304386d+00/
      data emf    /1836.152701d+00/
      data vl     /137.0359895d+00/
      data auev   /27.2113961d+00/
      data evmhz  /2.41798836d+08/
      data gmhz   /2.8025d+00/
c     data gn /5.58556d+00  , -4.255d+00 ,2.1709333d+00,-0.7800666d+00, 
c    1         1.7923333d+00,  1.4048d+00,0.4036d+00   ,-0.75748d+00  , 
c    2         5.2576d+00   , -0.4412d+00,1.4783333d+00,-0.34212d+00  , 
c    3         1.45656d+00  ,-1.11106d+00,2.2634d+00   , 0.4288666d+00, 
c    4         0.5478866d+00,-0.3714285d+00/
c     data isotop / 1, 3, 7, 9,11,13,14,17,19,
c    1             21,23,25,27,29,31,33,35,39/
c
      nospin=nprint.eq.-5
      if(nospin) return
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     +   (wfntyp.eq.mp2.and.zscftp.eq.uhf).or.
     +   (zruntp.eq.'ci'.and.isping.gt.0)) then
       nospin = .false.
      else
       nospin = .true.
      endif
      oroot = opg_root()
      if(nospin.and.oroot) then
       if(na.ne.nb) write(iwr,9993)
      endif
      if(nospin) return
c
      some=.true.
      out =nprint.eq.6
      ono = .false.
      if (oroot) then
        write(iwr,9999)
      endif
c
c     --- calculate constants and conversion terms --- 
c
      pi   = acos(-one)
      fac  =(four*pi/three)
      betae=one/(two*vl)
      betan=betae/emf
      convf=auev*evmhz
      con  =ge*betae*betan*convf
      deg  =angle/pi
c     ----- calculate spin density at the nuclei    -----
c           this is for the isotropic interaction
c           ( fermi contact term ) of the hyperfine
c           coupling tensor.
c
c     ----- special case = one particle reduced density matrix -----
c           is diagonal.
c
      l1= num
      l2=(num*(num+1))/2
      l3= num* num
      nptmax = nat
c
c     ----- allocate memory block -----
c
      need = l3 + l3 +  l2 + l2 + 10 * nptmax + 6 * nat
      init = igmem_alloc(need)
 
      call vclr(x(init),1,need)
c
c     ----- set pointers for partitioning of core -----
c
      i00 =init
      i10 =i00+nptmax*3
      i11 =i10+nptmax
      i20 =i11+nptmax*6
      i30 =i20 + l3
      i40 =i30 + l3
      i50 =i40 + l2
      i80 =i50 + l2
      last = i80 + 6 * nat
      need=last-1
c
c     ----- define points for calculation -----
c 
      do iat=1,nat
        x(1+3*(iat-1)+i00-1)=c(1,iat)
        x(2+3*(iat-1)+i00-1)=c(2,iat)
        x(3+3*(iat-1)+i00-1)=c(3,iat)
      enddo
c
c     ----- read molecular orbitals -----
c     ----- density matrix for x(i40)=da,x(i50)=db -----
c     ----- calculate density matrix -----
 
      call prpdens(x(i40),x(i50),x(i20),x(i30),zscftp,ono,oroot)
      if (ono) return
c
      if((wfntyp.eq.scf.and.zscftp.eq.uhf).or.
     1   (wfntyp.eq.mp2.and.zscftp.eq.uhf)    ) then
         if (iuno.gt.0.and.iunosp.gt.0) then
c     x(i40) and x(i50) contain the spinfree and spin density matrix
         else
c     x(i40) and x(i50) contain the -a- and -b- density matrix
         do i=1,l2
            duma      =x(i+i40-1)
            dumb      =x(i+i50-1)
            x(i+i40-1)=duma+dumb
            x(i+i50-1)=duma-dumb
         enddo
         endif
      else if (zruntp.eq.'ci'.and.isping.gt.0) then
          write(iwr,9979)
      else
      endif
c
      call vclr(x(i10),1,nptmax)
c
      if (wfntyp.eq.scf.and.zscftp.eq.uhf) then
        if (iuno.gt.0.and.iunosp.gt.0) then
c
c     use spin nos at x(i30)
c
        do imo=1,l1
        call elwfnc(nptmax,x(i00),x(i11),x(i30+l1*(imo-1)))
        dum = docc(imo)
         do ipt=1,nptmax
           x(ipt+i10-1)=x(ipt+i10-1)+dum*x(ipt+i11-1)**2
         enddo
        enddo
c
        else
c
        do imo=1,na 
        call elwfnc(nptmax,x(i00),x(i11),x(i20+l1*(imo-1)))
         do ipt=1,nptmax
           x(ipt+i10-1)=x(ipt+i10-1)+x(ipt+i11-1)**2
         enddo
        enddo
c
        do imo=1,nb
        call elwfnc(nptmax,x(i00),x(i11),x(i30+l1*(imo-1)))
          do ipt=1,nptmax
           x(ipt+i10-1)=x(ipt+i10-1)-x(ipt+i11-1)**2
          enddo
        enddo
c
        endif
c
      else if (isping.gt.0.and.zruntp.eq.'ci') then
c
c       use spin nos at x(i30)
c
        do imo=1,l1
        call elwfnc(nptmax,x(i00),x(i11),x(i30+l1*(imo-1)))
        dum = froccb(imo)
         do ipt=1,nptmax
           x(ipt+i10-1)=x(ipt+i10-1)+dum*x(ipt+i11-1)**2
         enddo
        enddo
c
      else
c
      endif
c
c     ----- print spin density and -aiso- at the nuclei ----- 
c
      if (oroot) then
        write(iwr,9994)
        write(iwr,9997)
      endif
c
      do ipt=1,nptmax
      iaz=nint(czan(ipt))
        if(iaz.le.30) then
        gnu=cspin(iaz) + cspin(iaz)
      else
        write(iwr,9985)
        gnu=one
      endif
c
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      spnden=x(ipt+i10-1)
      aiso  =two*fac*con*gnu*spnden
      if (oroot) 
     + write(iwr,9995) ipt,zaname(ipt),spnden,aiso,aiso/gmhz
      enddo
      if (oroot) write(iwr,9984)
c
c     ----- calculate spin-dipolar term at the nuclei ----- 
c           this is for the anisotropic interaction
c           of the hyperfine coupling tensor. 
c     ----- use spind density matrix at x(i50)
c
c     ----- calculate spin-dipolar term at all points ----- 
c
      call elfgrd(nptmax,x(i00),x(i11),x(i50))
c
      if (oroot) then
        write(iwr,9996)
        write(iwr,9992)
      endif
c 
      do ipt=1,nptmax
      iaz=nint(czan(ipt))
      xp=x(1+3*(ipt-1)+i00-1)
      yp=x(2+3*(ipt-1)+i00-1)
      zp=x(3+3*(ipt-1)+i00-1)
      spndxx  = x(1+6*(ipt-1)+i11-1) + fac*x(ipt+i10-1)
      spndyy  = x(2+6*(ipt-1)+i11-1) + fac*x(ipt+i10-1)
      spndzz  = x(3+6*(ipt-1)+i11-1) + fac*x(ipt+i10-1)
      spndxy  = x(4+6*(ipt-1)+i11-1)
      spndxz  = x(5+6*(ipt-1)+i11-1)
      spndyz  = x(6+6*(ipt-1)+i11-1)
      if (oroot)
     + write(iwr,9991) ipt,zaname(ipt),spndxx,spndyy,spndzz,
     +                                 spndxy,spndxz,spndyz
      x(1+6*(ipt-1)+i11-1) = spndxx
      x(2+6*(ipt-1)+i11-1) = spndyy
      x(3+6*(ipt-1)+i11-1) = spndzz
      x(4+6*(ipt-1)+i11-1) = spndxy
      x(5+6*(ipt-1)+i11-1) = spndxz 
      x(6+6*(ipt-1)+i11-1) = spndyz
      enddo
c
      if (oroot) then
        write(iwr,9974)
        write(iwr,9990)
      endif
c
      do ipt=1,nptmax
      iaz=nint(czan(ipt))
      if(iaz.le.30) then
        gnu=cspin(iaz) + cspin(iaz)
      else
        write(iwr,9985)
        gnu=one
      endif
      a(1) = x(1+6*(ipt-1)+i11-1)
      a(2) = x(4+6*(ipt-1)+i11-1)
      a(3) = x(2+6*(ipt-1)+i11-1)
      a(4) = x(5+6*(ipt-1)+i11-1)
      a(5) = x(6+6*(ipt-1)+i11-1)
      a(6) = x(3+6*(ipt-1)+i11-1)
      call diaspn(a,vec,eig,3)
      fac1 = con*gnu
      fac2 = con*gnu/gmhz
      x(1 + 6*(ipt-1) + i80 -1 ) = eig(1) * fac1
      x(2 + 6*(ipt-1) + i80 -1 ) = eig(2) * fac1
      x(3 + 6*(ipt-1) + i80 -1 ) = eig(3) * fac1
      x(4 + 6*(ipt-1) + i80 -1 ) = eig(1) * fac2
      x(5 + 6*(ipt-1) + i80 -1 ) = eig(2) * fac2
      x(6 + 6*(ipt-1) + i80 -1 ) = eig(3) * fac2
      if (oroot)
     + write(iwr,9989) ipt,zaname(ipt),eig(1),eig(2),eig(3)
      enddo
c
      if (oroot) then
       write(iwr,9970)
       do ipt = 1, nptmax
        ij = 6*(ipt-1) + i80 -1
        write(iwr,9989) ipt,zaname(ipt),(x(loop+ij),loop=1,3)
       enddo
       write(iwr,9960)
       do ipt = 1, nptmax
        ij = 6*(ipt-1) + i80 -1
        write(iwr,9989) ipt,zaname(ipt),(x(loop+ij),loop=4,6)
       enddo
      endif
c
      if (oroot) write(iwr,9988)
c
      do ipt=1,nptmax
      a(1) = x(1+6*(ipt-1)+i11-1)
      a(2) = x(4+6*(ipt-1)+i11-1)
      a(3) = x(2+6*(ipt-1)+i11-1)
      a(4) = x(5+6*(ipt-1)+i11-1)
      a(5) = x(6+6*(ipt-1)+i11-1)
      a(6) = x(3+6*(ipt-1)+i11-1)
      call diaspn(a,vec,eig,3)
      do j=1,3
       do i=1,3
        vec(i,j)= acos(vec(i,j))*deg
       enddo
      enddo
      if (oroot) then
       write(iwr,9987) ipt,zaname(ipt)
       write(iwr,9986) ((vec(i,j),j=1,3),i=1,3)
      endif
      enddo
c
c     ----- release memory block -----
c
      call gmem_free(init)
c
      return 
 9999 format(/10x,22('=')/10x,'hyperfine interactions'
     +       /10x,22('=')/)
 9997 format(3x,'atom ',11x,
     +          'density (a.u.)',4x,'aiso(mhz)',5x,'aiso(gauss)'/
     +          1x,60('-'))
 9996 format(/1x,'-------------------------------------------'
     +       /1x,'anisotropic interaction (spin-dipolar term)'
     +       /1x,'-------------------------------------------')
 9995 format(1x,i4,1x,a8,1x,f13.6,5x,f13.6,2x,f13.6)
 9984 format(1x,60('-'))
 9994 format(/1x,'---------------------------------------'
     +       /1x,'total spin density (fermi contact term)'
     +       /1x,'--------------------------------------')
 9993 format(/1x,
     +  'hyperfine interactions available for -uhf- and -ci- only.')
 9992 format(3x,'atom ',6x,'spin-dipolar term (a.u.)'/
     + 24x,'xx',13x,'yy',13x,'zz',13x,'xy',13x,'xz',13x,'yz'/
     + 1x,104('-'))
 9991 format(1x,i4,1x,a8,1x,6f15.6)
 9974 format(1x,104('-'))
 9990 format(/10x,'principal comp. of -hf- tensor (au)'/
     +        10x,'-----------------------------------')
 9970 format(/10x,'anisotropic -hfcc- (mhz)'/
     +        10x,'------------------------')
 9960 format(/10x,'anisotropic -hfcc- (gauss)'/
     +        10x,'--------------------------')
 9989 format(1x,i4,1x,a8,1x,3f13.6)
 9988 format(/' orientation (in degrees) of the principal axis',
     +        ' of hyperfine tensor w.r.t absolute (molecular) frame')
 9987 format(/' atom = ',i5, ' (',1x,a8,')'/' ------ ')
 9986 format(1x,3f15.6)
 9985 format(' ..... for atoms with z.gt.30 , neither the -aiso- ',
     +       ' nor the anisotropic -hfcc- have been multiplied',
     +       ' by -g(n)- ')
 9979 format(1x,'evaluate spin density based on spin NOs')
      end
      subroutine diaspn(a,vec,eig,ndim)
      implicit real*8 (a-h,o-z), integer   (i-n)
      dimension a(*),h(3,3),vec(3,3),eig(3),big(6)
      data zero   /0.0d+00/
      data one    /1.0d+00/
c          
      do i=1,3
         do j=1,3
            vec(j,i)=zero
         enddo
       vec(i,i)=one
      enddo
      ij=0
      do i=1,3
       do j=1,i
        ij=ij+1
        h(i,j)=a(ij)
        h(j,i)=a(ij)
       enddo
      enddo
      call diajac(h,vec,eig,ndim,ndim,big)
      return
      end
      subroutine ver_analg(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/analg.m,v $
     +     "/
      data revision /"$Revision: 6027 $"/
      data date /"$Date: 2009-08-03 12:26:49 +0200 (Mon, 03 Aug 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
