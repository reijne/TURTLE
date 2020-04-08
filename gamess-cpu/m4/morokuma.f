c     
c     Morokuma Energy Decomposition Analysis 
c
      subroutine morokscf(q)
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
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
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
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
      real*8 corev, array
      common/blkcore/corev(512),array(10)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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

c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      character *8 title,guess
      common/restrz/title(12),guess
c
      real*8 ptr, dtr, ftr, gtr, ppp, sta, cca, ra, scalea
      real*8 stb, ccb, rb, scaleb, ekk
      real*8 derror
      integer iposa, nstora, mpa, nspaca, iposb, nstorb
      integer mpb, nspacb, nsss, igvbo, igvb1, igsp
      integer intci, iso
      logical ondiis, odiisb
c
      common/junk/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     &     ppp(maxorb),
     &     sta(210),cca(20),ra(19),derror,scalea(20),
     &     iposa(20),nstora,mpa,ondiis,
     &     nspaca,stb(210),ccb(20),rb(20),scaleb(20),iposb(20),nstorb,
     &     mpb,odiisb,nspacb,nsss(30),
     &     igvbo(maxorb),igvb1(maxorb),igsp(4),
     &     ekk(63),intci(150), iso(1)

      integer ijkl
      common/craypk/ijkl(1360)
c

      real*8 q(*)

      character*30 filename
      logical o1e
      dimension o1e(6)
      character*8 zrhf
      character*4 yavr

      integer m1, m10, m16, igs
      integer l0, l1, l2, l3, l4, len2, lscdf, len
      integer i10,i20, i30, i21, ioc, ien
c     integer i40, i31
      integer i50, i60, iblk16

      integer ifq, ifo, ife, ifov
      integer ih0x, ihx

      integer isiz
      integer nw1d, nw2d, nex
      real*8 diaacc
      integer istat
      logical outon, odamph, oshift
      logical opr_skel, opr_symh, opr_vecs, opr_dens, opr_fock
      logical opr_conv, opr_parm
      integer ii, ij, j, k, icount, ix, ifr
      integer is, iff, lockt, lprnt, loop

      real*8 cpulft
      external cpulft

      real*8 chg_fa, coor_fa(3)
      real*8 ehf_frag(2), efrag
      integer nat_f1, num_f1
      integer nocmx
      logical oerr


      integer jblkfx
c, jblkqa
c  jblkpa, jblkea
      integer jblkh, jblkh0, jblkd
c  jblkqs
c jblkst, jblks, jblkf
      integer iblkh0, iblkqq
      integer iblkh

      integer iov, iovt, ike, ih0, ip
      integer iv, ivsave, ivt
      integer itmp

      real*8 timeit
      integer i,ibase, lwor,last
      integer nss, lword, nav, ndaf

c     integer mpunch
      integer lensec

      integer igmem_alloc, igmem_max_memory
      external igmem_alloc, igmem_max_memory
      
      integer lenwrd
      external lenwrd
      
      logical ochkfa
      external ochkfa

      real*8 eel0, eel, ees, eesx, exprime, epl, eex, eint, ect
      real*8 emix
c     real*8 eexpl, eplx

      real*8 dum

      real*8 enucf
      external enucf

      data zrhf /'rhf'/
c     data zcas   /'casscf','mcscf'/
c     data zuhf /'uhf'/

      data m1,m10,m16/1,10,16/

c     data zscf/'scf'/
c     data yblnk,yav/' ',' av'/
      data igs/5/

      call cpuwal(begin,ebegin)
      if(nprint.ne.-5)write(iwr,180)
 180  format(//1x,104('-'))

      nav = lenwrd()
c
c - most parms now set in iterate
c   maybe these should be to
c
      odiisb=.false.
      nstorb= 0
c     
c     ----- read in transformation matrices for s,p,d,f,g basis functions.
c     
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
      if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
      if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
      call readi(iso,nw196(5)*nav,ibl196(5),idaf)

      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = l2+l2
      len2=l2+l2

      iov = igmem_alloc(l2)     ! S
      ih0 = igmem_alloc(l2)     ! core hamiltonian
      ih0x = igmem_alloc(l2)    ! modified core hamiltonian
      ihx = igmem_alloc(l2)     ! modified 2e hamiltonian
      iovt = igmem_alloc(l2)    ! S (symmetry adapted basis)
      ike = igmem_alloc(l2)
c
      i10 = igmem_alloc(l2)

      do i =1,3
       o1e(i) = .true.
       o1e(i+3) = .false.
      enddo
      call getmat(q(iov),q(ike),q(ih0),
     +     q(i10),q(i10),q(i10),
     +     array,num,o1e,ionsec)

      call gmem_free(i10)
c
c     ----- transform the s-matrix
c     
      call tranp(q(iov),q(iovt))
c
c     ----- establish buffers for the two electron integral file(s)
c     
      if(nopk .ne. 1)
     &     call caserr('internal integral storage error')
      if (zscftp .ne. zrhf) call caserr('bad scf type')
      if (na .ne. nb) call caserr('must be closed shell')
c
      if(odscf) call caserr('morokuma cannot be direct')

      lfile=m2file
      do 141 i=1,lfile
         lotape(i)=m2tape(i)
         liblk(i)  =m2blk(i)
         llblk(i)  =liblk(i)-m2last(i)
 141  continue
c     
c     clear buffer for integral output
c     
      call icopy(1360,0,0,ijkl,1)
c     
c     atomic density startup for open and closed shells
c     modified for correct restart
c     
      if (guess.eq.'atoms') then
c     
c     do one cycle scf with the density matrix from denat
c     
         call denscf(q,zscftp)
c         if (irest.ne.0) go to 160
c     
c     set zguess to 'anything' so denscf will be called but once
c     
         zguess = 'anything'
         guess  = 'anything'
      end if
c     
c     ----- closed shell restricted hartree fock -----
c     
c     determine available memory
      lwor = igmem_max_memory()
c
c  we can can this when all memory is 
c  dynamic - for now it is over-generous
c
      nex = 0
      lscdf = 0
      if(lwor.gt.((35+nex)*l2+num+lscdf)) then
         len=(31+nex)*l2+num+lscdf
      else
         call caserr('not enough memory')
      endif
      nss=2
c some cockup here..
      lword=nss*len2+len+1000
      ibase = igmem_alloc(lword)
c
c  start of rhf driver
c
c     out = nprint .eq. 5
      outon = nprint .ne. -5
      timeit = cpulft(1)
c     
c     ----- set pointers for partitioning of core -----
c     
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav

      i50 = ibase
      i60 = i50+l2
c      i80 = i60+l2
c      i90 = i60+l1
      i10 = i60+l2
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
c     i31 = i20+l3
c     i40 = i30+l2
c     extra triangle for crystal energy correction
      jblkfx = i10+nss*l4
      isiz=0
      if(ocryst)isiz=l2
      jblkh0= jblkfx+isiz
      jblkh = jblkh0+l2
      jblkd = jblkh+l2   !symass
      ii    = jblkd+3*l2

      last = ii-ibase+1
      if(last.gt.lword)then
         write(iwr,9309)lword,last
         call caserr('insufficient memory available')
      endif

      lprnt = l1  ! print all virtuals

c
c  dynamic core
c
      ioc = igmem_alloc(l1)
      ien = igmem_alloc(l1)
      ip  = igmem_alloc(l2) 
      iv  = igmem_alloc(l3) 
      ivt  = igmem_alloc(l3) 
      ivsave  = igmem_alloc(l3) 
c
c     ----- occupation numbers -----
c     
c     -o- at q(ioc)

c     
c     ----- initialize variables -----
c     
      if (maxcyc .le. 0) maxcyc = 30
c     mpunch = 1
c     if (npunch .ne. 1) mpunch = 0
c     if (nprint .eq. 7) mpunch = 1
      i=lensec(l2)
      iblkh0 = ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)
      if(numdis.eq.0) numdis=num8

c
c  some of these are now reset within iterate
c  (shambles)
c
      if(irest.ne.3.or.numdis.eq.num8) then
         ehf = 0.0d0
         ehf0 = 0.0d0
         iter = 0
         kcount = 0
         iterv = 0
         damp = 0.0d0
         damp0 = 0.0d0
         rshift = 0.0d0
         diff = 0.0d0
         diffd = 0.0d0
         diffp = 0.0d0
         diffpp = 0.0d0
         de = 0.0d0
c        dep = 0.0d0
         deavg = 0.0d0
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(sta,nw2d,numdis)
      endif
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0

c     oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
c     
      if (odamph .or. oshift) damp = 1.0d0
c     

      if(outon)write (iwr,9348)
 9348 format(/
     +     40x,40('*') /
     +     40x,'    closed-shell rhf scf calculation  '/
     +     40x,' morokuma energy decomposition version'/
     +     40x,40('*'))

c
c  allocate additional pointers
c
c     ifq (fragment vectors)
c     ifo (fragment occupations)
c     ife (fragment orbital energies)
c     ifs (fragment overlap matrix - used for integrity test only)
c
      ifq = igmem_alloc(l3)
      ifo = igmem_alloc(l1)
      ife = igmem_alloc(l1)
      ifov = igmem_alloc(l2)
c
c
      imoro = ALL_BLOCKS

c
c print type of run,and load info from previous runs
c
      if(imode.eq.1)then
         write(iwr,*)'first fragment calculation fragment name=',
     &        fragname(1)
      else if(imode.eq.2)then
         write(iwr,*)'second fragment calculation fragment name=',
     &        fragname(2)
      else if(imode.eq.3)then
         write(iwr,*)'perform interaction energy analysis'
c
c  check fragment geometries and import vectors
c
         icount = 1
         ix = 0

         call dcopy(l3, 0.0d0, 0, q(ifq), 1)
         call dcopy(l2, 0.0d0, 0, q(ifov), 1)

         do ifr = 1,2

            write(iwr,*)'fragment ',ifr,' = ',fragname(ifr)
            
            call charlim(fragname(ifr),is,iff)
            filename=fragname(ifr)(is:iff)
            open(unit=57,file=filename,form='formatted',
     +       status='unknown')

            read(57,101,end=998)nat_f1
            do i=1, nat_f1
               read(57,103,end=998)chg_fa,coor_fa
               if(.not.ochkfa(  czan(icount), c(1,icount),
     &              chg_fa, coor_fa))then

                  write(iwr,*)'frag',ifr
                  write(iwr,*)'frag atom',i
                  write(iwr,*)'molec atom',icount
                  call caserr('atom mismatch')
               endif
               icount = icount + 1
            enddo

c
c  load fragment wavefunction
c
            read(57,101,end=998)num_f1
            do i=1,num_f1
               do j=1,num_f1
                  read(57,102,end=998)
     &                 q(ifq+(ix+i-1)*num+(ix+j-1))
               enddo
            enddo


c orbital energies
            read(57,101,end=998)num_f1
            do j=1,num_f1
               read(57,102,end=998)q(ife + (ix+j-1))
            enddo

c occupations..
            read(57,101,end=998)num_f1
            do j=1,num_f1
               read(57,102,end=998)dum
               q(ifo + (ix+j-1)) = dum

               if(dabs(dum) .lt. 1.0d-6)then
                  ifocc(ix+j) = 0
               else if(dabs(dum - 2.0d0) .lt. 1.0d-6)then
                  ifocc(ix+j) = 1
               else
                  call caserr('fragments must be closed shell')
               endif

            enddo
c
c total energy
            read(57,101,end=998)num_f1
            read(57,102,end=998)ehf_frag(ifr)
c
c fragment overlap matrix
            read(57,101,end=998)num_f1
            ij=(ix*(ix+1))/2 
            do i=1,num_f1
               ij=ij+ix
               do j=1,i
                  read(57,102,end=998)q(ifov +ij)
                  ij=ij+1
               enddo
            enddo

c set ifrag elements
            do j=1,num_f1
               ifrag(ix+j) = ifr
            enddo

            ix = ix + num_f1
            close(57)

         enddo


         write(iwr,*)'ifrag array',(ifrag(j),j=1,num)


         if(icount .ne. nat +1)call caserr('atom count error')
         if(ix .ne. num)call caserr('basis fn count error')

         write(iwr,*)'fragment atom lists were ok'

         write (iwr,*)'Fragment MO matrix'
         call prsq(q(ifq), l1, l1, l1)

         len = (l1-1)/lenwrd() + 1


c         write (iwr,*)'Fragment Overlap matrix'
c         call prtril(q(ifov),l1)

         itmp = igmem_alloc(l2)

         call dcopy(l2,q(iov),1,q(itmp),1)
         call ifzero(q(itmp),l1)

c         write (iwr,*)'Fragment Overlap matrix - from molecular '
c         call prtril(q(itmp),l1)

         oerr = .false.
         ij = 0
         do i = 1, l1
            do j = 1, i
               if(dabs(q(itmp + ij) 
     &              - q(ifov + ij)) .gt. 1.0d-8)then
                  write(iwr,*)i, j, q(itmp + ij), 
     &                 q(ifov + ij)
                  oerr = .true.
               endif
               ij = ij + 1
            enddo
         enddo

         if(oerr)then
            write(iwr,*)
     &           'Intrafragment overlap does not match that from'
            write(iwr,*)
     &           'fragments. Check that basis functions are in the'
            write(iwr,*)
     &           'same order in both fragments and complex.'
            call caserr('mismatched fragment/molecular overlaps')
         endif

         call gmem_free(itmp)

      else
         write(iwr,*)'test: mode flag=',imode
      endif

c     
c     ----- nuclear energy
c     
c  ocryst: en = crnuc(nat,czan,c,ecrn)
      en = enucf(nat,czan,c)
      write (iwr,9008) en
 9008 format(/,30h ----- nuclear energy ----- = ,f20.12)


c
c  some useful matrices for all cases
c
c
c  retain all virtuals
c
      l0=num

      if(imode.le.2)then
c
c   perform an SCF calculation on the fragment and save the results

c     
c restore vectors from guess (SAB)
c
         call rdedx(q(ivt),l3,ibl3qa,idaf)
c
c read orbital energies -> q(ien) and compute occupancies 
c
         call rdedx(q(ien),l1,ibl3ea,idaf)
         call llvmo(q(ien),q(ioc),na,nocmx,l1)
         call dscal(l1,2.0d0,q(ioc),1)
c
c      yavr=yblnk
c      if(nocmx.gt.na)yavr=yav
       lprnt = l1
c         if(.not.oprint(20)) lprnt=min(nocmx+5,l1)

c
c  full 2e hamiltonian
         oifc = .true.
         oife = .true.
         oifh = .true.
c
c no printing
         opr_skel  = .false.
         opr_symh  = .false.
         opr_fock  = .false.
         opr_vecs  = .false.
         opr_dens  = .false.
         opr_conv  = .true.
         opr_parm  = .true.
         lockt=lock

         call iterate(q,
     &        q(ih0), q(i10), q(ip), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstora, derror, iposa, diaacc, lprnt,
     &        iwr, num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .false.,
     &        istat)

         lock = lockt

c
c save molecular solution on dumpfile

         if(numdis.eq.num8)numdis=0
         ndaf = mouta
         call scfsav(q(ivt),q(ip),
     *        q(ien),q(ioc),ndaf,l1,
     *        l2,ibl3pa,ibl3ea)

         call dcopy(l3,q(ivt),1,q(i10),1)
         call dcopy(l1,q(ien),1,q(i21),1)

         call analmo(q(i10),q(i21),q(ioc),ilifq,l0,l1)
         call tdown(q(i10),ilifq,q(i10),ilifq,l0)

         if(otsym) 
     +        call symass(q(i10),q(i21),q(ioc),q)

         if(.not. oprint(25))then
            write (iwr,9148)
 9148       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors'/
     +           50x,12('-'))
            call prev(q(i10),q(i21),lprnt,l1,l1)
c     
            if(.not.outon)then
               write (iwr,9168)
 9168       format(/20x,14('*')/20x,'density matrix'/20x,14('*'))
               call prtril(q(i30),l1)
            endif
         endif

         if(istat .ne. 0)call caserr('fatal SCF convergence problem')
c
c  write results out for use in subsequent interaction analysis run
c
         call charlim(fragname(imode),is,iff)
         filename=fragname(imode)(is:iff)
         write(iwr,*)'*** saving vectors to ',filename
         open(unit=57,file=filename,form='formatted',status='unknown')
         write(57,101)nat
         do i=1,nat
            write(57,103)czan(i),(c(j,i),j=1,3)
         enddo

         write(57,101)num
c
c !!! should be in AO basis here
c
         write(57,102)(q(ivt + k - 1),k=1,l3)
         write(57,101)num
         write(57,102)(q(ien + k - 1),k=1,l1)
         write(57,101)num
         write(57,102)(q(ioc + k - 1),k=1,l1)
         write(57,101)num
         write(57,102)etot
         write(57,101)num
         write(57,102)(q(iov + k -1 ),k=1,l2)

 101     format(1x,i3)
 102     format(1x,e20.14)
 103     format(1x,4e20.14)

         close(57)

         array(1) = en
         array(2) = ehf
         array(3) = etot
         do loop = 4,10
          array(loop) = 0.0d0
         enddo
         call secput(isect(494),m16,m1,iblk16)
         call wrt3(array,m10,iblk16,idaf)

      else if(imode.eq.3) then

         write(iwr,*)'sum of fragment energies:'
         write(iwr,*)ehf_frag(1) + ehf_frag(2)

         write(iwr,*)'use fragment MO to startup'

c    try othog in different order
c         call dcopy(l3,q(ifq),1,q(ivt),1)
c         write(6,*)'Schmidt 2'
c ????? need to convert to SAB here
c         call dcopy(l3,q(ifq),1,q(ivt),1)
c hardwire for this case
c         ilo = 14
c         ihi = 26
c         jlo = 1
c         jhi = 5
c         call schm(q,q(ivt),q(iovt),ilo,ihi,jlo,jhi,l1)
c         write(6,*)'check transformed s'
c         call tran1eo(q, q(iovt), q(i10), q(ivt), l1) 
c         call prtril(q(i10),l1)
c
c !! need to convert to SAB here
c
         call dcopy(l3,q(ifq),1,q(ivt),1)
c
c  use occ and energy from fragment calculations
c
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num
c
c  use full 2e hamiltonian
         oifc = .true.  
         oife = .true.
         oifh = .true.
c
c no printing
         opr_skel  = .false.
         opr_symh  = .false.
         opr_fock  = .false.
         opr_vecs  = .false.
         opr_dens  = .false.
         opr_conv  = .true.
         opr_parm  = .true.
         lockt=lock

         imoro = ALL_BLOCKS
         
         call iterate(q, 
     &        q(ih0), q(i10), q(ip), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstora, derror, iposa, diaacc, lprnt,
     &        iwr, num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat)

         emorok(ALL_BLOCKS,YES_K) = eel + en

         emorok(ESX_EX,YES_K) = eel0 + en

c
c save molecular solution on dumpfile for possible further
c analysis
c
         if(numdis.eq.num8)numdis=0
         ndaf = mouta
         call scfsav(q(ivt),q(ip),
     *        q(ien),q(ioc),ndaf,l1,
     *        l2,ibl3pa,ibl3ea)

         array(1) = en
         array(2) = ehf
         array(3) = etot
         do loop = 4,10
          array(loop) = 0.0d0
         enddo
         call secput(isect(494),m16,m1,iblk16)
         call wrt3(array,m10,iblk16,idaf)
c

c ????? need to convert to SAB here

c
c  load fragment vectors, occupancies and energies
c  into the molecular arrays
c
         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oifc = .true.  
         oife = .true.         
         oifh = .true.
c flag retention of ESX (diagonal) blocks only
         imoro = ESX            

         write(iwr,*)'ESX'
         lockt=lock

         call iterate(q,
     &        q(ih0), q(i10), q(ip), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstora, derror, iposa, diaacc, lprnt,
     &        iwr, num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat)
         lock=lockt

         emorok(ESX,YES_K) = eel + en

         write(iwr,*)'ESX + PLX without exchange'

         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oife = .false. 
         maxcyc=30

         imoro = ESX_PLX

         call iterate(q,
     &        q(ih0), q(i10), q(ip), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstora, derror, iposa, diaacc, lprnt,
     &        iwr, num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat)

         emorok(ESX_PLX,NO_K) = eel + en

         emorok(ESX,NO_K)     = eel0 + en

         write(iwr,*)'ESX + CT'

         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oife = .true. 
         maxcyc=30

         imoro = ESX_CTT
         call iterate(q,
     &        q(ih0), q(i10), q(ip), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstora, derror, iposa, diaacc, lprnt,
     &        iwr, num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat)

         emorok(ESX_CTT,YES_K) = eel + en

         write(iwr,*)
         write(iwr,*)'       Morokuma  Analysis'
         write(iwr,*)

         efrag = ehf_frag(1) + ehf_frag(2)

         write(iwr,*)'Results of calculations: '

         write(iwr,104)
         call prtrow('Fragment 1          ',ehf_frag(1))
         call prtrow('Fragment 2          ',ehf_frag(2))

         call prtrow(' complex            ',emorok(0,1) - efrag)

         call prtrow(' exs                ',emorok(1,1) - efrag)
         call prtrow(' exs + plx          ',emorok(2,1) - efrag)
         call prtrow(' exs + ct           ',emorok(3,1) - efrag)
         call prtrow(' exs + ex^          ',emorok(4,1) - efrag)

         call prtrow(' exs no k           ',emorok(1,2) - efrag)
         call prtrow(' exs + plx no k     ',emorok(2,2) - efrag)
         call prtrow(' exs + ct no k      ',emorok(3,2) - efrag)
         call prtrow(' exs + ex^ no k     ',emorok(4,2) - efrag)

         write(iwr,104)

         write(iwr,*)

         write(iwr,100)

         ees = emorok(ESX,NO_K) - efrag
         call prtrow('Electrostatic  (Ees)',ees)

         eesx = emorok(ESX,YES_K) - efrag

cc         call prtrow('              (Eexs)',eesx)
cc         call prtrow('          (Eesx-Ees)',eesx-ees)

         exprime =  emorok(ESX_EX,YES_K) - emorok(ESX,YES_K)
cc         call prtrow('              (Eex^)',exprime)
cc         call prtrow('       No K   (Eex^)',
cc     &        emorok(ESX,NO_K) -  emorok(ESX_EX,NO_K))

         eex = exprime + eesx - ees
         call prtrow('Exchange       (Eex)',eex)

         epl = emorok(ESX_PLX,NO_K) - efrag - ees
         call prtrow('Polarisation   (Epl)',epl)

         ect = emorok(ESX_CTT,YES_K) -  emorok(ESX,YES_K)
         call prtrow('Charge Transfer(Ect)',ect)

cc         call prtrow('Chg. Tran. no K(Ect)',
cc     &        emorok(ESX_CTT,NO_K) -  emorok(ESX,NO_K))

c         eplx =  emorok(ESX_PLX,YES_K) - emorok(ESX,YES_K) 
c         call prtrow('              (Eplx)',eplx)

c eexpl currently missed out
c         eexpl = eplx - epl
c         call prtrow('Exchange-Pol (Eexpl)',eexpl)

         eint = emorok(ALL_BLOCKS,YES_K) - efrag

c eexpl currently missed out
c         emix = eint - ( ees + epl + eex + ect + eexpl)

         emix = eint - ( ees + epl + eex + ect)
         call prtrow('Coupling Term (Emix)',emix)

         write(iwr,104)
         call prtrow('Total  Interaction  ',eint)
         write(iwr,104)

 100     format(//1x,
     &'      Analysis following Morokuma - IJQC v. X p. 325 (1976)',/,
     &'      Term                E(H)      E(Kj/mol)  E(kcal/mol)',
     &        /,1x,60('-'))
 104     format(1x,60('-'))

      endif

c
c   free dynamic core
c      
      call gmem_free(ifov)
      call gmem_free(ife)
      call gmem_free(ifo)
      call gmem_free(ifq)
      
      call gmem_free(ivsave)
      call gmem_free(ivt)
      call gmem_free(iv)
      call gmem_free(ip)
      call gmem_free(ien)
      call gmem_free(ioc)

      call gmem_free(ibase)
      
      call gmem_free(ike)
      call gmem_free(iovt)
      call gmem_free(ihx)
      call gmem_free(ih0x)
      call gmem_free(ih0)
      call gmem_free(iov)
      
      return

 998  write(6,*)'unexpected end of stream on ',filename
      call caserr('fragment file error')

c999  write(6,*)'problem with ',filename
c     call caserr('error opening fragment file')

 9309 format(//
     *     1x,'words available ',i8/
     *     1x,'words requested ',i8//)
c9177 format(10x,a36,f8.2)
c9178 format(//
c    *     10x,'********************************************'/
c    *     10x,'Scf timing statistics              (seconds)'/
c    *     10x,'********************************************'/
c    *     10x,'fock formation                      ',f8.2/
c    *     10x,'fock diagonalisation                ',f8.2/
c    *     10x,'diis                                ',f8.2/
c    *     10x,'fock transformation                 ',f8.2/
c    *     10x,'orthogonalisation                   ',f8.2/
c    *     10x,'symmetrisation, damping, extrpn.    ',f8.2/
c    *     10x,'density matrix, etc.                ',f8.2/
c    *     10x,'initialisation                      ',f8.2/
c    *     10x,'dft                                 ',f8.2/
c    *     10x,'********************************************'/
c    *     10x,'Total                               ',f8.2/
c    *     10x,'********************************************'/)
      end

      subroutine charlim(s,iss,iff)
      character s*(*)
      iss=0
      iff=0
      do ii = 1,len(s)
         if(s(ii:ii) .ne. ' ')then
            if (iss.eq.0) iss = ii
         else if (iff.eq.0 .and. iss.ne.0)then
            iff=ii-1
         endif
      enddo
      return
      end

      logical function ochkfa(cz,c,czchk,cchk)
      implicit none
      real*8 cz,c(3),czchk,cchk(3),tol
      tol=1.0d-6
      if( (dabs(cz-czchk) .gt. tol) .or.
     &     (dabs(c(1)-cchk(1)).gt. tol) .or.
     &     (dabs(c(2)-cchk(2)).gt. tol) .or.
     &     (dabs(c(3)-cchk(3)).gt. tol) )then
         ochkfa = .false.
      else
         ochkfa = .true.
      endif
      return
      end
c

      subroutine iterate(q,
     &     h0,
     &     h  ,
     &     p  ,
     &     e  ,
     &     st  ,
     &     v  ,
     &     occ  ,
     &     hsav  ,
     &     vsav  ,  
     &     h00  ,
     &     gapa1, gapa2,
     &     ibrk,
     &     igs,
     &     pr_skel  ,
     &     pr_symh  ,
     &     pr_fock  ,
     &     pr_vecs  ,
     &     pr_dens  ,
     &     pr_parm  ,
     &     pr_conv  ,
     &     lock  ,
     &     ondiis  ,
     &     nstore  ,
     &     derror  ,
     &     iposit, 
     &     diaacc,
     &     lprnt  ,
     &     iwr  ,
     &     num  ,
     &     iky,
     &     ilifq,
     &     nocmx  ,
     &     na  ,
     &     yavr,
     &     eel0,
     &     eel,
     &     fq,
     &     lock_occupations,
     &     istat)
c     
      implicit none
c     
      character*10 charwall
c     arguments

      real*8 q(*)                 ! base memory
      real*8 h0(*)                ! core hamiltonian
      real*8 h(*)                 ! full hamiltonian
      real*8 p(*)                 ! density
      real*8 e(*)                 ! orbital energies
      real*8 st(*)                ! overlap q(jblkst)
      real*8 v(*)                 ! vectors  ????
      real*8 occ(*)               ! occupancies
      real*8 hsav(*)              ! saved hamiltonian
      real*8 vsav(*)              ! vectors q(jblkqa)
      real*8 h00(*)               ! q(jblkh0) - used in diis
      real*8 gapa1, gapa2         ! level shifter
      integer ibrk              ! when to switch level shifters
      integer igs               ! how often to reorthogonalise
      logical pr_skel           ! print skeleton fock
      logical pr_symh           ! print symmetried fock
      logical pr_fock           ! print total fock
      logical pr_vecs           ! print eigenvectors
      logical pr_dens           ! print eigenvectors
      logical pr_parm           ! print scf parameters
      logical pr_conv           ! print convergence data
      integer lock              ! 
      logical ondiis            ! diis common /junk/
      integer nstore            ! number of diis vectors in use
      real*8 derror               ! diis error
      integer iposit(*)         ! addresses of diis fock matrices
      real*8 diaacc

      integer lprnt             ! number of orbitals to print
      integer iwr               ! output stream
      integer num               ! size of basis
      integer iky(*)
      integer ilifq(*)          ! transformation table
      integer nocmx             ! number of occ orbs
      integer na                ! number of alpha electrons
c to flag partial occ of degenerate homos
      character*4 yavr         
c electronic energy of input orbitals
      real*8 eel0                 
c electronic energy of converged solution
      real*8 eel                 

      real*8 fq(*)                ! fragment orbitals
c force FMO occupation pattern on molecule
      logical lock_occupations  
      integer istat             ! return code
c     
c     local variables
c
      logical finished
      integer iscr              ! scratch memory pointer
      integer iscr1, iscr2      ! scratch memory pointer
      integer iscr1x, iscr2x, iscr3x ! scratch memory pointer
      integer ifsave            ! scratch memory pointer
      integer iqdiis            ! diis array
      integer iht, i30

      integer ifqi, ihx, ih0x, istx, iovq

      integer l0, l1, l2, l3
c     integer l4
      integer i, j, ii, ij, loop
      integer ig, m2, len

      real*8 diffdp,diffdpp
      real*8 dep
      real*8 dmptlc, dmplim
      real*8 ehf1, ehf2
      real*8 tim1
      real*8 skale, tlefti, tlefts
      real*8 deter

      logical oshift
      logical odamph
      logical oextra
      logical ocvged

      integer lockt

      logical odbg
      logical omemck

      real*8 dum

      real*8 qmax
      integer iforb, ifpvi, ivsav
      logical assign
      logical first_cycle
c     
c dynamic core
c     
c     external functions
c
      external igmem_alloc
      integer igmem_alloc

      integer lenwrd
      external lenwrd

      external tracep
      real*8 tracep

      external idamax
      integer idamax

      real*8 cpulft
      external cpulft

c
c convergence control parameters
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c

      integer select(maxorb), nsel
      integer ioccsav(maxorb)

      character*4 yblnk, yav

      data yblnk,yav/' ',' av'/

      diffdpp = 1.0d4
      diffdp = 1.0d4

      l0=num                    ! number of canonical orbitals kept
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
c     l4 = l2+l2

      odbg = .false.
      omemck = .false.

      dmptlc =1.0d-2

      lockt=lock

      istat=0

      call timrem(tlefts)

c
c     ----- decode convergence method -----
c           mconv is in /scfopt/
c
      oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
c
c    ---- some tolerancies
c
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
c
c  old fock matrices for extrap
      ifsave = igmem_alloc(3*l2 + 1000)

c
c flag for FMO-based locking

      first_cycle = .true.

      if(odiis)then
c
c old fock matrices for diis 

         iqdiis = igmem_alloc(16*l2 + 1000)
         ii=1
         do 180 loop=1,8
            iposit(loop)=ii
            ii=ii+l2+l2
 180     continue
      endif
c
c  inverse eigenvector matrix for backtransformations
c
      ifqi = igmem_alloc(l3)

c
c old vectors for locking
c
      ivsav = igmem_alloc(l3)

      len = (l1-1)/lenwrd() + 1

      iscr1 = igmem_alloc(len)
      iscr2 = igmem_alloc(len)

      call dcopy(l1*l1,fq,1,q(ifqi),1)

      if(imode.gt.2)then
         call minvrt(q(ifqi),l1,deter,
     &        q(iscr1),q(iscr2))
      endif

      write(iwr,*)'deter',deter

      call gmem_free(iscr2)
      call gmem_free(iscr1)
c
c  these probably should be l2, but for mult2 in diis
c
      iscr1x = igmem_alloc(l3) 
      iscr2x = igmem_alloc(l3) 
      iht = igmem_alloc(l3)     ! transformed H
      i30 = igmem_alloc(l3)     ! scratch

      write(iwr,*)'copy ie'
c
c  local copies of S, and H0
c
      ih0x = igmem_alloc(l2) 
      call dcopy(l2,h0,1,q(ih0x),1)

      ihx = igmem_alloc(l2)
cccc      call dcopy(l2,h,1,q(ihx),1)

      istx = igmem_alloc(l2)
      call dcopy(l2,st,1,q(istx),1)

c      write(iwr,*)'mod 1e'
c
c  modify S, H0 according to Morokuma scheme
c

      call modmoro( q, q(istx), fq, q(ifqi), l1)
      call modmoro( q, q(ih0x), fq, q(ifqi), l1)

cdbg
c      write(iwr,*)'modified overlap'
c      call prtril(q(istx),l1)

cdbg
c      write(iwr,*)'modified h0'
c      call prtril(q(ih0x),l1)
c
c  qq(ivoff+) eigenvectors
c
      write(iwr,*)'diag q'
      iovq = igmem_alloc(l3)

      call dcopy(l2,st,1,q(iscr1x),1)

      iscr3x = igmem_alloc(max(l3,l1+1))
      call qmat(q,q(iscr1x),q(iovq),q(iscr2x),
     &     q(iscr3x),iky,l0,l1,l3,l1,.false.)
      call gmem_free(iscr3x)

      write(iwr,*)'symmetrisation'
c
c  orthogonalise input orbitals according to modified S
c  the symmetric scheme
c

      nsel = 0
      do i = 1,l1

         dum = occ(i)
         if(dabs(dum) .lt. 1.0d-6)then
            select(i) = 0
         else if(dabs(dum - 2.0d0) .lt. 1.0d-6)then
            select(i) = 1
            nsel = nsel + 1
         else
            call caserr('systems must be closed shell')
         endif

      enddo

c      write(iwr,*)'initial check s'
c      call tran1eo(q(istx), q(iscr1x), v, l1) 
c      call prtril(q(iscr1x),l1)

      call dcopy(l3,v,1,q(iscr1x),1)
      call symmorth2(q,q(iscr1x), v, 
     &     select,q(istx),l1,nsel,iwr)

c      write(iwr,*)'check transformed s'
c      call tran1eo(q(istx), q(iscr1x), v, l1) 
c      call prtril(q(iscr1x),l1)

c      write(iwr,*)'input orbs'
c      call prsq(v,l1,l1,l1)
c
c  - initial density matrix  !!! symmetry adaption
c
c
c transform vectors (SA basis) ->  (AO basis)
c
c         call tdown(q(iv),ilifq,q(ivt),ilifq,l1)

      call dmtx(p,v,occ,iky,nocmx,l1,l1)

cc      write(iwr,*)'initial density'
cc      call prtril(p,l1)

      finished = .false. 

      if(pr_parm)then
         write(iwr,*)'SCF parameters'
         write(iwr,*)'maxcyc =',maxcyc
         write(iwr,*)'oifc   =',oifc
         write(iwr,*)'oife   =',oife
         write(iwr,*)'oifh   =',oifh
         write(iwr,*)'odiis  =',odiis
         write(iwr,*)'accdi1 =',accdi1
         write(iwr,*)'acurcy =',acurcy
         write(iwr,*)'dmpcut =',dmpcut
         write(iwr,*)'oshift =',oshift
         write(iwr,*)'odamph =',odamph
         write(iwr,*)'oextra =',oextra
         write(iwr,*)'shifters',gapa1,gapa2,ibrk
      endif

      lockt=lock
c
      if(pr_conv)then
         write (iwr,9028)
 9028    format(/1x,100('=')/
     *        3x,'cycle',10x,'total',5x,'electronic',8x,'e conv.',
     *        9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     *        17x,'energy',9x,'energy',35x,'shift'/1x,100('='))
      endif
c
c   some iteration control parameters (see also calling program)
c
c     
c     ----- set convergence criteria -----
c     - see also subroutine iterate
c

c diis
      nstore= 0
      accdi2=acurcy*acurcy
      ondiis=.false.

      ehf = 0.0d0
      ehf0 = 0.0d0
      iter = 0
      kcount = 0
      iterv = 0
      damp = 0.0d0
      damp0 = 0.0d0
      rshift = 0.0d0
      diff = 0.0d0
      diffd = 0.0d0
      diffp = 0.0d0
      diffpp = 0.0d0
      de = 0.0d0
      dep = 0.0d0
      deavg = 0.0d0
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0
      if (odamph .or. oshift) damp = 1.0d0
c
      call dcopy(l3, v, 1, vsav, 1)
c
      do while (.not. finished) 

         if(omemck)call ma_summarize_allocated_blocks
c
c   construct skeletonised 2e part of fock matrix
c
         call mhstar(p, h)

         if (pr_skel) then
            write (iwr,9068)
 9068       format(/
     +           20x,20('*')/
     +           20x,'skeleton fock matrix'/
     +           20x,20('*'))
            call prtril(h,l1)
         endif
c     
c     ----- symmetrize skeleton fock matrix -----
c     
         call symh(h,q(iscr1x),iky,0,0)
c
c  remove inter-fragment terms from h in AO basis
c
         if(.not. oifh)call ifzero(h,l1)
c
c   scheme 2 more comprehensive tinkering
c
         call modmoro(q,  h, fq, q(ifqi), l1)

         if (pr_symh) then
            write (iwr,9088)
 9088       format(/
     +           20x,23('*')/
     +           20x,'symmetrized fock matrix'/
     +           20x,23('-'))
            call prtril(h,l1)
         endif
c     
c     ----- read in core hamiltonian matrix
c     and calculate hf energy -----
c     
c     for xtal field calculations, the core hamiltonian
c     includes the crystal potential, but for the energy 
c     all terms modelling molecule-molecule interactions 
c     must be divided by two 
c     

         ehf0 = ehf
         call vadd(h,1,q(ih0x),1,h,1,l2)
         ehf1 = tracep(p,q(ih0x),l1)
         ehf2 = tracep(p,h,l1)
         ehf = (ehf1+ehf2) * 0.5d0

         if(iter.eq.0)eel0 = ehf
c     
c     if(ocryst)then
c     call secget(isecfx,m53,iblk)
c     call rdedx(q(jblkfx),l2,iblk,idaf)
c     ehf3 = tracep(p,q(jblkfx),l1)
c     write(6,*)'ehf3',ehf3
c     esave = ehf
c     ehf = (ehf1+ehf2-ehf3)*pt5
c     ecre = esave - ehf
c     write(6,*)'new ehf, ecre',ehf,ecre
c     endif
c     
c     ----- save fock matrix h -> hsav
c     
         call dcopy(l2,h,1,hsav,1)

         if (pr_fock) then
            write (iwr,9048)
 9048       format(
     +           20x,11('*'),/
     +           20x,'fock matrix'/
     +           20x,11('*'))
            call prtril(h,l1)
         endif
c     
         if(maxcyc .eq. 0)then
            finished = .true.
            goto 400
         endif

         iter = iter+1
         etot = ehf+en
         dep = de
         de = ehf-ehf0
         if (iter .eq. 1) then
            deavg = 0.0d0
         else if (iter .eq. 2) then
            deavg =  dabs(de)
         else if (iter .ge. 3) then
            deavg = (  dabs(de)+  dabs(dep)+0.2d0*deavg)/2.2d0
         endif
c     
c     ----- damp and extrapolate hamiltonian matrix -----
c     
c     -h - at x(i10)     hamiltonian matrix (n th iteration)
c     -ho- at x(i20)     old h matrix (n-1 th iteration)
c     -ha- at x(i30)     ancient h matrix (n-2 th iteration)
c     -hp- at x(i40)     prehistoric h matrix (n-3 th iteration)
c     
         if (iter .gt. 2) call dampd(de,dep,deavg,damp,
     +        acurcy,diff,diffp,dmptlc)
c     
         if (damp .lt. dmpcut) damp = dmpcut

c         write(iwr,*)'iter, iterv',iter, iterv
c         write(iwr,*)'extrpm',de, damp, damp0

         if(omemck)write(iwr,*)'2'
         if(omemck)call ma_summarize_allocated_blocks

      iscr3x = igmem_alloc(l3)
         call extrpmm(de,damp,damp0,h,
     &        q(iscr1x),q(iscr2x),q(iscr3x),
     &        l1, l2, iterv,1,1,q(ifsave))
      call gmem_free(iscr3x)

         diffpp=diffp
         diffp=diff
         diffdpp=diffdp
         diffdp=diffd

         if(omemck)write(iwr,*)'2.5'
         if(omemck)call ma_summarize_allocated_blocks

c         write(iwr,*)'h after extrap'
c         call prtril(h,l1)
c         write(iwr,*)'odiis',odiis
c
c    takes over scratch arrays from extrpmm
c  - should replace these when extrpmm loses
c    its l2*3 old matrices
c
c    h is copied to n**2 space to allow for
c    mult2 workspace
c
         if(odiis) then
c
c   st should be in non-SAB basis here!!!

            call dcopy(l2,h,1,q(iscr1x),1)

      iscr3x = igmem_alloc(l3)
            call diiscmm(q(iqdiis),q(iscr1x),
     &           q(iscr2x), q(iscr1x),q(iscr3x),
     &           h00,vsav,q(istx),diff)
c            write(iwr,*)'diis res',ondiis,diff

      call gmem_free(iscr3x)

            call dcopy(l2,q(iscr1x),1,h,1)

         endif
         diffd=diff

         if(odbg)write(iwr,*)'extrap fock'
         if(odbg)call prtril(h,l1)

         if(omemck)write(iwr,*)'3'
         if(omemck)call ma_summarize_allocated_blocks
c     
c     ------ take precautions if tester is increasing again
c     
         if (iter.ne.1 .and. ondiis .and. 
     1       diffd.ge.diffdp .and. diffdp.ge.diffdpp) then
            nstore=0             
            if(pr_conv)write(iwr,*)'tester up!'
            call dcopy(l2,hsav,1,h,1)
            ondiis=.false.
         endif
c     
c     ----- read in fock tranformation matrix -----
c     transform hamiltonian matrix
c     
c     - q- at x(i30) orthonormalizing transformation
c     - h- at x(i10)
c     -h'- at x(i50) transformed h matrix
c     
c     ----- if vshift is true, use the previous set of
c     molecular orbitals to transform the hamiltonian matrix
c     
c     tdown   q(i30) - result
c     q(jblkqa,s) vsav, qs - input
c     
         if(oshift) then
            call tdown(q(i30),ilifq,vsav,ilifq,l0)
         else
            call tdown(q(i30),ilifq,q(iovq),ilifq,l0)
         endif
c     
c     q(iht) = Qt H qq(ivoff+),  H= h, qq(ivoff+)=q(i30)
c     
         call tran1eo(q, h, q(iht), q(i30), l1) 

         if(omemck)write(iwr,*)'4'
         if(omemck)call ma_summarize_allocated_blocks

         if(odbg)then
            write(iwr,*)'transformed fock'
            call prtril(q(iht),l1)
         endif

         diff=0.0d0
         ii=nocmx+1
         if(ii.le.l0)then
            do 250 i=ii,l0
               ij=iky(i)+iht
               loop=idamax(nocmx,q(ij),1)
               if(loop.gt.0) diff=dmax1(diff,dabs(q(ij+loop-1)))
 250        continue
         endif
         if(ondiis)diff=diffd
         dmplim = dmax1(dmpcut,2.0d0)

c         write(iwr,*)'after check fock',nocmx,diff

         ocvged = (damp.lt.dmplim) 
     &        .and. (diff.lt.acurcy) 
     &        .and. (iter.gt.1)

         if (ocvged) go to 321
c     
c     shift the diagonal of the transformed
c     h matrix by rshift for the virtual part.
c     
         rshift=0.0d0
         if(.not.ondiis)
     *        call shiftq(q(iht),nocmx,0,l0,
     *        de,dep,iterv,1,gapa1,ibrk,gapa2)

c         write(iwr,*)'after shiftq',nocmx,rshift
c     
c     ----- diagonalize new hamiltonian matrix -----
c     
c     -h- at x(i50)
c     -v- at x(i30)
c     
         m2=2
         diaacc = diff*5.0d-3
         if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
         if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c     
         call jacobi(q(iht),iky,l0,v,ilifq,l1,
     &        e,m2,lock,diaacc)
c     
c     ----- back-transform the eigenvectors -----
c     
c     -v- at x(i30)
c     -d- at x(i41)
c     -q- at x(i10)
c     scratch area at x(i50)
c     
c     transform vt with q(jblka,s)
c     
         if(omemck)call ma_summarize_allocated_blocks

         iscr = igmem_alloc(l3)
         if(oshift) then
            call tfsqc(v,vsav,q(iscr),l0,l1,l1)
         else
            call tfsqc(v,q(iovq),q(iscr),l0,l1,l1)
         endif
         call gmem_free(iscr)
c     
c     ----- if vshift is true, reorthogonalize the vectors and
c     thereby the transformation matrix every igs th iteration.
c     
         ig = mod(iter,igs)
         if (ig .eq. 0) then

            if(pr_conv)write(iwr,*)'re-orthogonalize..'

            iscr1 = igmem_alloc(l3)
            iscr2 = igmem_alloc(l3)
c
c  tranform S
c
            call tran1eo(q, q(istx), q(iscr1), v, l1)
c     
c     orthogonalise v  according to S in q(iscr1)
c 
            call orfog(v,v,q(iscr1),q(iscr2),
     &           iky,ilifq,l1,l1,1)

            call gmem_free(iscr2)
            call gmem_free(iscr1)

            if(omemck)call ma_summarize_allocated_blocks

         endif
c     
c     copy orthogonalised vectors -> vsav
c     
         call dcopy(l3,v,1,vsav,1)
c     
c     transform v to AO basis
c     
         call tdown(v,ilifq,v,ilifq,l0)
c     
c     unshift virtual energies
c     
         ig=nocmx+1
         if(ig.le.l0)call vsadd(e(nocmx+1),1,-rshift,
     &        e(nocmx+1),1,l0-nocmx)

         if (pr_vecs) then
c     
c     eigenvectors
c     
            write (iwr,9148) 
 9148       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors'/
     +           50x,12('-'))
            call prev(v,e,lprnt,l1,l1)
         endif
c     
c     ----- form density matrix -----
c     
c     -v- at x(i30)
c     -o- at x(i80)
c     
c     load occs -> occ
c
      if(lock_occupations)then

c
c disable the llvmo call below
         assign = .false.
c
c locking of occupations according to fragment occupations
c
         if(first_cycle)then

c
c determine MO occupancies from FMO parentage
c
c transform to fragment orbital basis
            iscr = igmem_alloc(l3)
            call tfsqc(v,q(ifqi),q(iscr),l0,l1,l1)

c            write (iwr,91481) 
c91481       format(//1x,100('-')//
c     +           50x,12('-')/
c     +           50x,'eigenvectors in frag MO basis'/
c     +           50x,12('-'))
c            call prev(v,e,lprnt,l1,l1)
            yavr=yblnk
            do i = 1, num
               qmax = 0.0d0
               do j = 1, num
                  if( abs(v((i-1)*num + j)) .gt. abs(qmax))then
                     qmax = v((i-1)*num + j)
                     iforb = j
                  endif
               enddo
               occ(i) = 2.0d0 * dble( ifocc(iforb) )
               if(odbg)write(6,*)'MO ',i,' FMO ',iforb, qmax, occ(i)
            enddo
c
c check counts, abort if number of electrons changed
            j = 0
            do i = 1, num
               if(occ(i) .gt. 1.0d-10)then
c                  write(6,*)'occ',i,occ(i)
                  nocmx = i
                  j = j + 1
               endif
            enddo
            if(j .ne. na)then
               write(6,*)'counts',j,na
               write(6,*)'problem setting occupancies'
               assign = .true.
            endif
c back transform the vectors 
            call tfsqc(v,fq,q(iscr),l0,l1,l1)

            call gmem_free(iscr)
c
c store current orbs and occupancies for next time
c
            call dcopy(l3,v,1,q(ivsav),1)
            do i = 1, num
               ioccsav(i)=0
               if(occ(i) .gt. 1.0d-10)ioccsav(i) = 1
            enddo

            first_cycle = .false.

         else
c
c try and set occs based on orbitals of previous cycle
c

            iscr = igmem_alloc(l3)
            ifpvi = igmem_alloc(l3)
            iscr1 = igmem_alloc(len)
            iscr2 = igmem_alloc(len)

            call dcopy(l1*l1,q(ivsav),1,q(ifpvi),1)
            call minvrt(q(ifpvi),l1,deter,
     &           q(iscr1),q(iscr2))

            if(odbg)write(iwr,*)'deter pv',deter

c     use orbitals of previous cycle
            call tfsqc(v,q(ifpvi),q(iscr),l0,l1,l1)

            call gmem_free(iscr2)
            call gmem_free(iscr1)
            call gmem_free(ifpvi)

            if(odbg)then
            write (iwr,91482) 
91482       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors in old MO basis'/
     +           50x,12('-'))
            call prev(v,e,lprnt,l1,l1)
            endif
c
            yavr=yblnk

            do i = 1, num
               qmax = 0.0d0
               do j = 1, num
                  if( abs(v((i-1)*num + j)) .gt. abs(qmax))then
                     qmax = v((i-1)*num + j)
                     iforb = j
                  endif
               enddo
               occ(i) = 2.0d0 * dble( ioccsav(iforb) )
               if(odbg)
     &          write(iwr,*)'MO ',i,' old MO ',iforb, qmax, occ(i)
            enddo

            j = 0
            do i = 1, num
               if(occ(i) .gt. 1.0d-10)then
                  nocmx = i
                  j = j + 1
               endif
            enddo
            if(j .ne. na)then
               write(iwr,*)'problem setting occupancies'
               write(iwr,*)'counts',j,na
               assign = .true.
            endif

c previous cycle
            call tfsqc(v,q(ivsav),q(iscr),l0,l1,l1)
            call gmem_free(iscr)

c
c now copy the current orbitals and occupations
            call dcopy(l3,v,1,q(ivsav),1)
            do i = 1, num
               ioccsav(i)=0
               if(occ(i) .gt. 1.0d-10)ioccsav(i) = 1
            enddo

         endif
      else
c no locking, use llvmo
         assign = .true.
      endif
c
c make assignment on basis of energy
c

      if(assign)then
         call llvmo(e,occ,na,nocmx,l1)
         call dscal(l1,2.0d0,occ,1)
         yavr=yblnk
         if(nocmx.gt.na)yavr=yav
      endif

         call dmtx(p,v,occ,iky,nocmx,l1,l1)
         if (pr_dens) then
            write (iwr,9168)
 9168       format(/20x,14('*')/20x,'density matrix'/20x,14('*'))
            call prtril(p,l1)
         endif

 321      continue

         tim1=cpulft(1)
c     
c     ----- check and print convergence behavior -----
c     
         if(pr_conv)then
            write (iwr,9188) iter,kcount,etot,
     1           ehf,de,diff,rshift,damp,derror,yavr
 9188       format(1x,i3,1x,i3,4f15.8,f10.3,f8.3,f15.9,a3)
         endif

         dmplim = dmax1(dmpcut,2.0d0)

         ocvged = (damp.lt.dmplim) .and. 
     &        (diff.lt.acurcy) .and. 
     &        (iter.gt.1)

         if (ocvged) go to 400
c     
c     ----- exit in case of time limit -----
c     ----- punch out restart data -----
c     
         call timit(0)
         call timrem(tlefti)

         skale=1.1d0
         if (tlefti .lt. skale*(tlefts-tlefti)/iter) then
            write(iwr,*)'scf timed out'
            istat = 1
            finished = .true.
         else if (iter .gt. maxcyc) then
c     
c     too many iterations
c     
            write (iwr,9288)
 9288       format(/
     +           20x,30('-')/
     +           20x,'excessive number of iterations'/
     +           20x,30('-'))
            etot = 0.0d0
            ehf0 = ehf
            ehf = -en
            istat = 2
            finished = .true.
         endif

 400     continue

         if(ocvged)finished = .true.

      enddo
c
      if(ocvged)eel = ehf

      if(pr_conv)then

         if (ocvged)then
            write (iwr,9208)
 9208       format(/20x,16('-')/20x,'energy converged'/,20x,16('-'))
         else
            write (iwr,9209)
 9209    format(/20x,22('!')/
     +           20x,'*** no convergence ***'/,20x,22('!'))
         endif

         write (iwr,9368) iter,tim1,charwall(),ehf,en,etot
         if(yavr.ne.yblnk)write(iwr,9369)
      endif

 9368 format(  /20x,14('-')/
     * 20x,'final energies   after',i4,' cycles at ',f8.2,' seconds'
     *     ,a10,' wall'/
     *     20x,14('-')//
     *     20x,'electronic energy',f18.10/
     *     20x,'nuclear energy   ',f18.10/
     *     20x,'total energy     ',f18.10)
 9369 format(//20x,90('*')/
     *     20x,'*',88x,'*'/20x,'*',88x,'*'/
     *     20x,'*',31x,'warning  state is averaged',31x,'*'/
     =     20x,'*',88x,'*'/20x,90('*')//)

c
c  return SA vectors in v as well as vsav
c
      call dcopy(l3, vsav, 1, v, 1)
c
c  release dynamic memory
c
      call gmem_free(iovq)
      call gmem_free(istx)
      call gmem_free(ihx)
      call gmem_free(ih0x)
      call gmem_free(i30)
      call gmem_free(iht)
ccc      call gmem_free(iscr3x)
      call gmem_free(iscr2x)
      call gmem_free(iscr1x)
      call gmem_free(ivsav)
      call gmem_free(ifqi)
      if(odiis)call gmem_free(iqdiis)
      call gmem_free(ifsave)
c
c  this helps catch memory corruption 'till we
c  have diis/extrap sorted
c
c      write(iwr,*)'memory check:'
c      
c      call ma_summarize_allocated_blocks

      return

      end
c
c  slightly modified diis, to reflect change in
c  memory addressing - h00, v, and s passed directly by
c  address rather than offset
c     
c  until mult2 situation is clear, h1, h2 should probably
c  be given a num*num square of memory for mult2 work.
c
c ===
c ===  6-triangle variant of diis - h0 overlays h2
c ===
c
      subroutine diiscmm(q,h0,h1,h2,h3,h00,v,s,diff0)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension q(*),h0(*),h1(*),h2(*),h3(*),s(*),v(*),h00(*)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
      common/junk/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     * pp(maxorb),
     * st(210),ct(20),r(19),derror,scale(20),iposit(20),nstore
     *,mp,ondiis
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
c
      common/scftim/tdiag(4),tdiis,tmult(5)
c
c
c  diiscm is specifically for the closed shell  cases
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written out to memory
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored in memory, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c jdiis  logical is set to true only if the diis solution
c        is to be used in the scf program
c h0     current array (not yet stored)
c h1     initially the array from the previous cycle
c h2+h3  scratch arrays --- taken over from extrap
c ndaf   the starting location for the vectors .. need space for 20
c accdi1  diis comes in only when diff is less than this
c accdi2  the diis solution is only used when the residuals
c          are less than thia
c
c num    the dimension of the arrays
c nx     the size of the arrays
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      derror=0.0d0
c     l3=num*num
      nmin=3
      iblk=1
      nmax=8
      len2=nx+nx
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
c      write(iwr,*)'diis test',nmax, nmin, ispace
      if(nmin.lt.1) goto 300
      call dcopy(nx,h0,1,h00,1)

c      write(iwr,*)'diis input fock'
c      call prtril(h0,num)
c
c ----- sort out the vectors
c
      call tdown(h3,ilifq,v,ilifq,num)
c
c ----- calculate the error vector
c
      call mult2(h3,h1,h2,num,num,num)

c      write(iwr,*)'diis tr fock'
c      call prtril(h1,num)

      diff0=0.0d0
      ij=0
      do 5 i=1,num
      do 4 j=1,i
      ij=ij+1
      if(i.le.na.and.j.le.na) h1(ij)=0.0d0
      if(i.gt.na.and.j.gt.na) h1(ij)=0.0d0
      if( dabs(h1(ij)).gt.diff0) diff0= dabs(h1(ij))
 4    continue
 5    h1(ij)=0.0d0
      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) then
c         write(iwr,*)'fail diis test',diff0, accdi1
         goto 310
      else
c         write(iwr,*)'pass diis test',diff0, accdi1
      endif
c
c ------ transform back to the ao basis
c
      do 10 i=1,num
      do 10 j=1,i
      dum=h3(ilifq(i)+j)
      h3(ilifq(i)+j)=h3(ilifq(j)+i)
 10   h3(ilifq(j)+i)=dum
      call mult2(h3,h2,h1,num,num,num)
      call square(h3,s,num,num)
      call mult2(h3,h1,h2,num,num,num)
      call dcopy(nx,h00,1,h0,1)
      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1
      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.lt.1) goto 100
      do 50 i=1,ipos1
      ibl=iposit(i)+nx
      mp=mp+1
      st(mp)=tracep(h1,q(ibl),num)
 50   continue
 100  mp=mp+1
      st(mp)=tracep(h1,h1,num)
      itot=nstore
      if(nstore.gt.nmax)itot=nmax
      call dcopy(nx,h0,1,q(iposit(ipos)),1)
      call dcopy(nx,h1,1,q(iposit(ipos)+nx),1)
      ipos1=ipos+1
      if(ipos1.gt.itot) goto 115
      do 110 i=ipos1,itot
      ibl=iposit(i)+nx
      mp=mp+i-1
      st(mp)=tracep(h1,q(ibl),num)
 110  continue
 115  continue
c
c --- now solve the diis equations
c
      if(nstore.le.nmin)then
c         write(6,*)'diis tfail est nstore',nstore,nmin
         goto 400
      endif
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h3,h3(ndim+1),ifail)
c     call fixnag(ifail,'f04atf-diiscm')
      if(ifail.ne.0.and.oprint(47))write(6,129)ifail
 129   format(//1x,'diis failure in f04atf .... ifail= ',i2)
      if(ifail.ne.0) goto 210
      do 114 i=1,ndim
  114  ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c --- now read through the diis file again
c
c     del=0.0d0
      call vclr(h0,1,nx)
      do 150 k=1,itot
      ibl=iposit(k)
      call daxpy(nx,ct(k),q(ibl),1,h0,1)
 150  continue
      if(.not.ondiis) kcount=0
c      write(6,*)'diis on'
      ondiis=.true.
      go to 9999
 210  ibl=iposit(ipos)
      call dcopy(nx,q(ibl),1,h0,1)
      call dcopy(len2,q(ibl),1,q(iblk),1)
      nstore=1
      itot=1
      mp=iky(ipos)+ipos
      st(1)=st(mp)
      go to 410
 310  continue
      call dcopy(nx,h00,1,h0,1)
 300  nstore=0
      itot=0
      mp=0
      go to 410
 400  call dcopy(nx,h00,1,h0,1)
 410  ondiis=.false.
c      write(6,*)'diis is off'
9999  continue
      return
      end

c         call extrpmm(de,damp,damp0,h,
c     &        qq(ivoff+iscr1x),qq(ivoff+iscr2x),qq(ivoff+iscr3x),
c     &        l1,l2,
c     +        1,                ! index into core 
c     +        iterv,1,1,qq(ivoff+ifsave) )

c
c  modified so as not to assume that h1, h2, and h3
c  are contiguous in memory, and remove any dependence
c  on memory addressing
c
c  old matrices do assume this, and are passed as a 3*l2 
c  array (hold) to reflect this
c
      subroutine extrpmm(dee,dampp,damp00,h0,h1,h2,h3,l1,l2,
     +     iterl,ncall,ityp,hold)
c
c     -----  use either davidson's damping and extrapolation
c            or pople's extrapolation -----
c
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
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
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
      common/extsav/odampr,oextpr
      dimension hold(l2,*),h0(*),h1(*),h2(*),h3(*)
      data dzero,done,two,four/0.0d0,1.0d0,2.0d0,4.0d0/
      data pt5/0.5d0/
      data tol,tol1,tol2/1.0d-07,1.9d0,0.99d0/
      data shrnkf,dmpmin,rsmin/100.0d0,0.01d0,0.8d0/
      data zgvb,zgrhf /'gvb','grhf'/
c
c      ndaf2 = ndaf1+l2
c      ndaf3 = ndaf2+l2

      len3 = l2*3
c
c     ----- decode convergence method -----
c
      oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3

      if (ncall .eq. 1) kcount = kcount+1

      if (iter .eq. 1) then
c
c     ----- iter = 1 -----
c
c
c  - load stored matrices
c
         call dcopy(l2,h0,1,hold(1,1),1)
         call dcopy(l2,h0,1,hold(1,2),1)
         call dcopy(l2,h0,1,hold(1,3),1)
         return
      else
c
c     ----- current fock matrix is in -h0-
c           get previous fock matrices (h1,h2,h3) -----
c
         call dcopy(l2,hold(1,1),1,h1,1)
         call dcopy(l2,hold(1,2),1,h2,1)
         call dcopy(l2,hold(1,3),1,h3,1)

         if (iter .eq. 2) then
c     
c     ----- iter = 2 -----
c     
            if (odamph) go to 320
            if ( .not. oshift) go to 420
            if ((zscftp.eq.zgvb.or.zscftp.eq.zgrhf)
     &           .and. odiis) go to 420
            dampp = done
            go to 320
            
         endif
c     
c     ----- iter > 2 -----
c     
         if (oshift) go to 220
         if (odampr) go to 200
         if ( .not. odamph) go to 420
         if ( dabs(dee) .gt. exttol) go to 180
         if (dee .gt. dzero .and. kcount .gt. iextin) go to 180
         go to 420
 180     dmptol = dmptol/shrnkf
         exttol = exttol/shrnkf
         go to 300
 200     if ( dabs(dee) .gt. dmptol) go to 300
         if (dee .gt. dzero) go to 300
         if (dampp .gt. dmpmin) go to 300
         go to 420
         
 220     if ( .not. oextpr) go to 260
         if ( dabs(dee) .gt. exttol) go to 240
         if (dee .gt. dzero .and. kcount .gt. iextin) go to 240
         go to 420
 240     exttol = exttol/shrnkf
         dmptol = dmptol/shrnkf
         rshift = rsmin
         iterl = 0
         if (odamph) go to 300
         go to 620
 260     if ( dabs(dee) .gt. dmptol) go to 280
         if (dee .gt. dzero) go to 280
         if (odamph .and. dampp .gt. dmpmin) go to 280
         if ( .not. oextra) go to 280
         if (rshift .ge. vshtol) go to 280
         if (iterl .eq. 0) go to 280
         go to 420
 280     if (odamph) go to 300
         if (iterl .lt. 2) go to 620
         if ( .not. oextra) go to 620
         go to 440
c     
c     ----- davidson's damping -----
c     
 300     if (kcount .lt. iextin) go to 320
c     
         if ( .not. oshift .or. iterl .ge. 2) go to 360
 320     temp=done/(done+dampp)
         do 340 i = 1,l2
            h0(i) = (h0(i)+dampp*h1(i))*temp
 340     continue
         go to 400
c     
c     ----- iter > 2 , damping -----
c     
 360     cutoff = pt5-damp00
         temp=done/(done+dampp)
         do 380 i = 1,l2
            h0(i) = (h0(i)+dampp*h1(i))*temp
 380     continue
         if ( .not. oextra .or. cutoff .lt. dzero) go to 400
         odampr = .true.
         oextpr = .true.
c     
         go to 460
 400     if (ncall .ne. 1) go to 660
         odampr = .true.
         oextpr = .false.
         go to 660
 420     if ( .not. oextra) go to 620
c     
c     ----- pople's extrapolation procedure
c     
         if (ncall .ne. 1) go to 460
         odampr = .false.
         oextpr = .true.
 440     dampp = dzero
c     
c     ----- skip to end if first cycle or after extrapolation -----
c     
 460     call dcopy(len3,h0,1,hold(1,1),1)
         call dcopy(len3,h1,1,hold(1,2),1)
         call dcopy(len3,h2,1,hold(1,3),1)
         call vsub(h2,1,h3,1,h3,1,l2)
         call vsub(h1,1,h2,1,h2,1,l2)
         call vsub(h0,1,h1,1,h1,1,l2)
         if (kcount .lt. iextin .or. iter .lt. 4) go to 680
c     
c     ----- find displacement dp1,dp2,dp3 -----
c     
         if (ityp .eq. 2) then
            sp11 = ddot(l2,h1,1,h1,1)
            sp12 = ddot(l2,h2,1,h1,1)
            sp13 = ddot(l2,h3,1,h1,1)
            sp22 = ddot(l2,h2,1,h2,1)
            sp23 = ddot(l2,h3,1,h2,1)
            sp33 = ddot(l2,h3,1,h3,1)
         else
            sp11 = tracep(h1,h1,l1)
            sp12 = tracep(h2,h1,l1)
            sp13 = tracep(h3,h1,l1)
            sp22 = tracep(h2,h2,l1)
            sp23 = tracep(h3,h2,l1)
            sp33 = tracep(h3,h3,l1)
         endif
         dp1 = dsqrt(sp11)
         dp2 = dsqrt(sp22)
c     
c     ----- find cosine of angle between successive displacements -----
c     
         dp3 = dsqrt(sp33)
c     
c     ----- find cosine of angle between -dp(3)- and
c     plane of =dp(1)- and -dp(2)-.
c     
         cosphi = sp12/(dp1*dp2)
         r = sp11*sp22-sp12*sp12
         p = (sp13*sp22-sp12*sp23)/r
         q = (sp23*sp11-sp12*sp13)/r
c     
c     ----- do not extrapolate unless -4- consecutive points are
c     nearly coplanar -----
c     
         cospsi = dsqrt(p*p*sp11+q*q*sp22+two*p*q*sp12)/dp3
         if (cospsi .le. tol) go to 680
c     
c     ----- express -dp(1)- as x*dp(3)(projected)+y*dp(2) -----
c     
         if (dampp .gt. dmpmin) go to 680
         q = -q/p
c     
c     ----- test if 2*2 matrix has real eigenvalues
c     between -tol/2 and +tol/2 -----
c     
         p = done/p
         pq = q*q+four*p
         if (pq .lt. dzero) go to 680
         pq =  dabs(q)+dsqrt(pq)
c     
c     ----- if -4- point extrapolation is not possible,
c     try -3- point
c     
         if (pq .le. tol1) go to 560
         if ( dabs(cosphi) .le. tol2) go to 680
         p = dp1/(dp2*cosphi-dp1)
         call daxpy(l2,p,h1,1,h0,1)
         go to 600
 560     ppp = p/(done-p-q)
         qqq = (p+q)/(done-p-q)
         do 580 i = 1,l2
            h0(i) = h0(i)+ppp*h2(i)+qqq*h1(i)
 580     continue
c     
 600     kcount = 0
         call dcopy(l2,h0,1,hold(1,1),1)
         go to 680
         
c     
c     ----- no damping or extrapolation -----
c     
 620     if (ncall .ne. 1) go to 640
         odampr = .false.
         oextpr = .false.
c     
c     ----- save new (modified) fock matrix -----
c     
 640     dampp = dzero
c     
 660     call dcopy(l2,h0,1,hold(1,1),1)
         call dcopy(l2,h1,1,hold(1,2),1)
         call dcopy(l2,h2,1,hold(1,3),1)

      endif

  680 return
      end


      subroutine mhstar(d, f)
c
c
c     ----- -hstar- forms a skeleton matrix -----
c                   f=( h* + h )/2
c
c              f(i,j)=(h**(i,j) + h**(j,i))/2
c
c     indices in labels are in standard order-
c      i.ge.j , k.ge.l , (ij).ge.(kl)
c
c     all contributions are made into lower half of
c     skeleton matrix.
c     only off-diagonal elements need be divided by two,
c     to obtain the correct f matrix.
c
c  mhstar - modified and simplified for morokuma energy 
c           analysis.
c
c  assumes that supermatrix option is off
c
c  oife controls inclusion of inter-fragment exchange integrals
c  oifc controls inclusion of inter-fragment coulomb  integrals
c
c  basis functions are classified 
c
c     ----- ***********************************  ------
c
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
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      dimension d(*),f(*)
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
      common/blkin/gin(510),nint
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical osortp,oschw
      common/sortpdat/
     + igbf,iklbf,iiptbf,ijbase,ngbf,lengbf,osortp,oschw
c
      logical omem, ofirt
      integer maxq, maxqq, iqqoff, mxblock
c
      common/maxlen/maxq,maxqq,omem(maxlfn),iqqoff(maxlfn),
     +              ofirt,mxblock
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      call vclr(f,1,nx)
c
      ifile=1
      main=lotape(ifile)

      do 1 m=1,num
         nij=ikyp(m)
         d(nij)=d(nij)*0.5d0
 1    continue
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
      do 1000 ifile=1,lfile

         main=lotape(ifile)
         if(omem(main))
     &      call caserr('integral storage error - no incore option!')
         call search(liblk(ifile),main)
         call find(main)
         jblock=llblk(ifile)
 1004    jblock=jblock+1
         call get(gin,nw)
         if(nw.eq.0)go to 1005
         if(jblock.ne.0)call find(main)

         if(o255i) then
            call sgmata_morok(f,d)
         else
c
c Apply some checks
c Assume we are building with the DFT code, but DFT is 
c not switched on
c
          call sgmata_255_morok(f,d)
       endif
c
         if(jblock)1004,1000,1004
 1005    llblk(ifile)=liblk(ifile)-iposun(main)+1

 1000 continue
c
      do 3   m = 1,num
         nij = ikyp(m)
         d(nij) = d(nij)+d(nij)
 3    continue
      call dscal(nx,0.5d0,f,1)
c
c  parallel dgop here
c
      return
      end

      subroutine sgmata_morok(fock, p)
      implicit real*8  (a-h,o-z)
      dimension p(*),fock(*)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      logical otestif, oif
      external otestif
c
      common/blkin/gg(510),mword
      common/craypk/integ(1)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
      iword = 1
      do 6000 iw=1,mword

         j = integ(iword )
         i = integ(iword+1)
         l = integ(iword+2)
         k = integ(iword+3)

         oif = otestif(i,j,k,l)

c... coulomb
         gik=gg(iw)
         g2=gik+gik
         g4=g2+g2
         ikyi=iky(i)
         ikyj=iky(j)
         ikyk=iky(k)
         ik=ikyi+k
         il=ikyi+l
         ij=ikyi+j
         jk=ikyj+k
         jl=ikyj+l
         kl=ikyk+l
         if((.not.  oif) .or. oifc)then
            aij=g4*p(kl)+fock(ij)
            fock(kl)=g4*p(ij)+fock(kl)
            fock(ij)=aij
         endif

c... exchange

         if((.not.  oif) .or. oife)then
            gil=gik
            if(i.eq.k.or.j.eq.l)gik=g2
            if(j.eq.k)gil=g2
            if(j.ge.k)goto 1
            jk=ikyk+j
            if(j.ge.l)goto 1
            jl=iky(l)+j
 1          ajk=fock(jk)-gil*p(il)
            ail=fock(il)-gil*p(jk)
            aik=fock(ik)-gik*p(jl)
            fock(jl)=fock(jl)-gik*p(ik)
            fock(jk)=ajk
            fock(il)=ail
            fock(ik)=aik
         endif

         iword=iword+4

 6000    continue
         return
      end

      subroutine sgmata_255_morok(fock,p)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

      logical oexch, ocoul, onodft
c
c flag specifying if we are doing a dft calc
c
      logical occpdft
      common/ccpdft/occpdft
c
c  CCPDFT API declarations
c
      logical CD_2e
      integer CD_4c2eon
      integer CD_abort
      integer CD_accuracy
      logical CD_active
      integer CD_assign_grid
      integer CD_ang_npoints_row
      integer CD_auto_ang_prune
      logical CD_check_print
      integer CD_chf_dksm_mo
      integer CD_chf_lhs_ao
      integer CD_chf_lhs_mo
      integer CD_chf_rhs_ao
      integer CD_chf_rhs_mo
      integer CD_clone_grid
      integer CD_conv_prune_on
      integer CD_create_grid
      integer CD_debug
      integer CD_defaults
      integer CD_defaults_old
      integer CD_dksm_exp_ao
      integer CD_dksm_exp_mo
      integer CD_euleron
      integer CD_gausslon
      integer CD_generation
      logical CD_gradcorr
      integer CD_gradquad
      integer CD_gridatomradius
      integer CD_gridscale
      integer CD_hess_ao
      integer CD_hess_mo
      logical CD_HF_coulomb
      logical CD_HF_coulomb_deriv
      logical CD_has_HF_exchange
      real*8    CD_has_HF_exchange_weight
      logical CD_HF_exchange
      real*8    CD_HF_exchange_weight
      logical CD_ignore_accuracy
      integer CD_init
      integer CD_inttol
      logical CD_is_rks
      logical CD_is_jfiton
      logical CD_is_jfitmem
      integer CD_jfit_clean1
      integer CD_jfit_clean2
      logical CD_jfit_incore
      integer CD_jfit_init1
      integer CD_jfit_init2
      integer CD_jfitoff
      integer CD_jfiton
      integer CD_jfitgon
      integer CD_jmulton
      integer CD_lebedevon
      integer CD_logon
      integer CD_lypon
      integer CD_memreq_chf_dksm_ao
      integer CD_memreq_chf_dksm_mo
      integer CD_memreq_chf_lhs_mo
      integer CD_memreq_chf_rhs_ao
      integer CD_memreq_chf_rhs_mo
      integer CD_memreq_energy
      integer CD_memreq_energy_ao
      integer CD_memreq_energy_mo
      integer CD_MHL_ang_prune
      integer CD_over
      integer CD_pener
      integer CD_pole
      integer CD_pruneatomradius
      integer CD_psitol
      integer CD_rad_npoints_row
      integer CD_radscale_scheme
      integer CD_request
      logical CD_request_multstate
      integer CD_reset_2e
      integer CD_rks
      integer CD_schwarz
      integer CD_screen
      integer CD_screenatomradius
      integer CD_set_2e
      integer CD_set_functional
      integer CD_set_ignore_accuracy
      integer CD_set_print_level
      integer CD_set_weight
      integer CD_sortpoints
      integer CD_energy
      integer CD_energy_ao
      integer CD_energy_mo
      integer CD_forces_ao
      integer CD_uks
      integer CD_xcfiton
      integer CD_import_geom
      integer CD_update_geom
      integer CD_warn
      integer CD_weightatomradius
c
      integer gden_init
      integer gden_energy
      integer gden_forces
c
c     declare API routines as external
c
      external CD_2e
      external CD_4c2eon
      external CD_abort
      external CD_accuracy
      external CD_active
      external CD_assign_grid
      external CD_auto_ang_prune
      external CD_check_print
      external CD_chf_dksm_ao
      external CD_chf_dksm_mo
      external CD_chf_lhs_ao
      external CD_chf_lhs_mo
      external CD_chf_rhs_ao
      external CD_chf_rhs_mo
      external CD_clone_grid
      external CD_conv_prune_on
      external CD_create_grid
      external CD_debug
      external CD_defaults
      external CD_defaults_old
      external CD_dksm_exp_ao
      external CD_dksm_exp_mo
      external CD_euleron
      external CD_gausslon
      external CD_generation
      external CD_gradcorr
      external CD_gradquad
      external CD_gridatomradius
      external CD_gridscale
      external CD_hess_ao
      external CD_hess_mo
      external CD_HF_coulomb
      external CD_HF_coulomb_deriv
      external CD_has_HF_exchange
      external CD_has_HF_exchange_weight
      external CD_HF_exchange
      external CD_HF_exchange_weight
      external CD_ignore_accuracy
      external CD_init
      external CD_inttol
      external CD_is_rks
      external CD_is_jfiton
      external CD_is_jfitmem
      external CD_jfit_clean1
      external CD_jfit_clean2
      external CD_jfit_incore
      external CD_jfit_init1
      external CD_jfit_init2
      external CD_jfitoff
      external CD_jfiton
      external CD_jfitgon
      external CD_jmulton
      external CD_lebedevon
      external CD_logon
      external CD_memreq_chf_dksm_mo
      external CD_memreq_chf_lhs_mo
      external CD_memreq_chf_rhs_mo
      external CD_memreq_energy
      external CD_memreq_energy_ao
      external CD_memreq_energy_mo
      external CD_MHL_ang_prune
      external CD_over
      external CD_pener
      external CD_pole
      external CD_psitol
      external CD_request
      external CD_request_multstate
      external CD_reset_2e
      external CD_rks
      external CD_schwarz
      external CD_screen
      external CD_screenatomradius
      external CD_set_2e
      external CD_set_functional
      external CD_set_ignore_accuracy
      external CD_set_print_level
      external CD_set_weight
      external CD_sortpoints
      external CD_energy
      external CD_energy_ao
      external CD_energy_mo
      external CD_forces_ao
      external CD_uks
      external CD_xcfiton
      external CD_import_geom
      external CD_update_geom
      external CD_warn
      external CD_weightatomradius
c
      external gden_init
      external gden_energy
      external gden_forces
c
c print control
c
      integer PRINT_NONE
      parameter(PRINT_NONE=0)

      integer PRINT_LOW
      parameter(PRINT_LOW=2)

      integer PRINT_DEFAULT
      parameter(PRINT_DEFAULT=5)

      integer PRINT_HIGH
      parameter(PRINT_HIGH=7)

      integer PRINT_ALL
      parameter(PRINT_ALL=10)
      dimension p(*),fock(*)
      integer *2 integ
      common/blkin/gg(255),integ(1020),mword
c
      logical oif, otestif
      external otestif


      onodft = .not. CD_active()
      ocoul = onodft .or. CD_HF_coulomb()
      oexch = CD_HF_exchange()
      if(oexch) then
        facex = CD_HF_exchange_weight()
      else
        facex = 1.0d0
      endif
c
      iword = 1
      do 6000 iw=1,mword
      j = integ(iword )
      i = integ(iword+1)
      l = integ(iword+2)
      k = integ(iword+3)

      oif = otestif(i,j,k,l)

c... coulomb

      gik=gg(iw)
      g2=gik+gik
      g4=g2+g2
      ikyi=iky(i)
      ikyj=iky(j)
      ikyk=iky(k)
      ik=ikyi+k
      il=ikyi+l
      ij=ikyi+j
      jk=ikyj+k
      jl=ikyj+l
      kl=ikyk+l
c
c coulomb term
c
      if(ocoul .and. ((.not.  oif) .or. oifc))then
        aij=g4*p(kl)+fock(ij)
        fock(kl)=g4*p(ij)+fock(kl)
        fock(ij)=aij
      endif
c
c exchange
c
        if(onodft)then

        if((.not.  oif) .or. oife)then
c
c full term for HF
c
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1010
        jk=ikyk+j
        if(j.ge.l)goto 1010
        jl=iky(l)+j
 1010   ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik

        endif

      else if(oexch)then
c
c DFT exchange with scaling factor
c
        call caserr('morokuma not working with dft')

        gik = gik*facex

        g2=gik+gik
        gil=gik
        if(i.eq.k.or.j.eq.l)gik=g2
        if(j.eq.k)gil=g2
        if(j.ge.k)goto 1
        jk=ikyk+j
        if(j.ge.l)goto 1
        jl=iky(l)+j
 1      ajk=fock(jk)-gil*p(il)
        ail=fock(il)-gil*p(jk)
        aik=fock(ik)-gik*p(jl)
        fock(jl)=fock(jl)-gik*p(ik)
        fock(jk)=ajk
        fock(il)=ail
        fock(ik)=aik
      endif
 6000 iword=iword+4
      return
      end


      logical function otestif(i,j,k,l)

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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)


      ifi = ifrag(i) 
      ifj = ifrag(j) 

      if(ifi .ne. ifj) goto 100

      ifk = ifrag(k) 
      ifl = ifrag(l) 

      if(ifk .ne. ifl) goto 100

      if(ifi .ne. ifk) goto 100

      otestif = .false.
      return

 100  otestif = .true.
      return

      end
c
c   zeros the interfragment part of
c   a one-electron operator in MO basis
c   (assumes triangular form)
c
      subroutine ifzero(m,n)
      implicit none
      real*8 m(*)
      integer n

      integer i, j, ij

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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)


      ij = 1
      do i = 1, n
         do j= 1, i
            if(ifrag(i) .ne. ifrag(j))m(ij) = 0.0d0
            ij = ij + 1
         enddo
      enddo

      end
c
c  uses a set of MO vectors to transform a one-electron
c  opeator to the MO basis
c
      subroutine tran1eo(q,in,out,v,n)
      implicit none

      real*8 in(*), out(*)
      integer n
      real*8 v(n,n),q(*)

      integer igmem_alloc
      external igmem_alloc

      integer l2, l3, itmp1, itmp2

      l2 = (n*(n+1))/2
      l3 = n*n

      itmp1 = igmem_alloc(l3)
      itmp2 = igmem_alloc(l3)

      call dcopy(l2,in,1,q(itmp1),1)

      call mult2(v,q(itmp2),q(itmp1),n,n,n)

      call dcopy(l2,q(itmp2),1,out,1)

      call gmem_free(itmp2)
      call gmem_free(itmp1)

      return
      end

c
c************************************************************************
c
c   subroutine symmorth2
c
c  symmetric orthogonalisation using s**1/2 - 
c  this version designed to handle only a subset 
c  of orbitals and uses schmidt orthogonalisation
c  to orthogonalise the remainder
c
c  s is provided in the *AO* basis 
c
c  set select(i)=1 for the symm-orth set only
c
c************************************************************************
      subroutine symmorth2(q,vin,vout,select,s,n,nsel,iwr)

      implicit none

      integer n, nsel, select(n),iwr
      real*8 vin(n,*), vout(n,*), s(*), q(*)

      integer igmem_alloc
      external igmem_alloc

      real*8 ddot
      external ddot


      integer itr, iscr1, iscr2, iscr3, iscr4, iscr5, iscr6
      integer i, j, ij, k

      real*8 dum, fact

c??   itr = igmem_alloc(10*n*n)
      itr = igmem_alloc(n*n)
      write(iwr,*)'itr',itr

      iscr6 = igmem_alloc((n*(n+1))/2)
      write(iwr,*)'transform S to MO basis'
c
c transform S to MO basis
c
      call tran1eo(q, s, q(iscr6), vin, n)
c
c expand to square
c
      write(iwr,*)'expand to square'
      call square(q(itr),q(iscr6),n,n)

      call gmem_free(iscr6)


      iscr1 = igmem_alloc(nsel*nsel)
      iscr2 = igmem_alloc(nsel*nsel)
      iscr3 = igmem_alloc(nsel)
      iscr4 = igmem_alloc(nsel)
      iscr5 = igmem_alloc(nsel*nsel)
c
c  copy required elements to matrix nsel*nsel 
c
      ij = 0
      do i=1,n
         if(select(i).eq.1)then
            do j=1,n
               if(select(j).eq.1)then
                  q(iscr5 + ij) 
     &          = q(itr + (i-1)*n + (j-1))
                  ij = ij + 1
               endif
            enddo
         endif
      enddo
c
c  find s**-1/2
c
      write(iwr,*)'find s**-1/2'
      call rootmtm(q(iscr5),q(iscr1),
     &     q(iscr2),q(iscr3),
     &     q(iscr4),nsel,nsel,1,0,iwr)
c
c  construct a transformation matrix from nonzero 
c  s**-.5 blocks and unit diagonal elsewhere
c
      write(iwr,*)'return from rootmtm'
      call dcopy(n*n,0.0d0,0,q(itr),1)
      ij = 0
      do i = 1, n
         if(select(i).eq.1)then
            do j = 1, n
               if(select(j).eq.1)then
                  q(itr + (i-1)*n + (j-1)) = 
     &                 q(iscr5 + ij) 
                  ij = ij + 1
               endif
            enddo
         else
            q(itr + (i-1)*n + (i-1)) =  1.0d0
         endif
      enddo

      call gmem_free(iscr5)
      call gmem_free(iscr4)
      call gmem_free(iscr3)
      call gmem_free(iscr2)
      call gmem_free(iscr1)
c
c  transform vectors
c
         call mxmg(vin,1,n,q(itr),1,n,vout,1,n,n,n,n)

c
c now schmidt orthogonalise all unselected orbitals
c
c  square S in AO basis
c
      call square(q(itr),s,n,n)
      iscr1 = igmem_alloc(n*n)
******
      do i = 1,n
       do j = 1,n
        q(iscr1+i+(j-1)*n-1) = ddot(n,q(itr+(i-1)*n),1,vout(1,j),1)
       enddo
      enddo

      do i = 1,n
         if(select(i).eq.0)then
c
c  schmidt
c
            do j = 1,n
               if( (select(j) .ne. 0) .or. (j .lt. i)) then
                  dum = 0.0d0
                  do  k = 1 , n
                   fact = vout(k,i)
                   if (fact.ne.0.0d0) dum = dum +
     +               q(iscr1+k+(j-1)*n-1)* fact
                  enddo 
                 call daxpy(n
     +                 ,-dum,vout(1,j),1,vout(1,i),1)
               endif
            enddo
c
c  normalise 
c
            dum = 0.0d0
            do j = 1 , n
               if (vout(j,i).ne.0.0d0) dum = dum +
     +              ddot(n
     +              ,q(itr+(j-1)*n),1
     +              ,vout(1,i),1)*vout(j,i)
            enddo
            dum = 1.0d0/dsqrt(dum)
c            write(iwr,*)'shmidt norm',i,dum
            call dscal(n,dum,vout(1,i),1)
         endif
      enddo

      write(iwr,*)'free',iscr1
      call gmem_free(iscr1)
      write(iwr,*)'free',itr
      call gmem_free(itr)

      end
**==rootmtm.f
      subroutine rootmtm(a,b,res,aa,bb,mdim,nbas,inv,nignore,iwr)
c
c     this routine takes the matrix currently in a and returns
c     either its square root or its inverse square root in a.  the
c     eigenvectors of the matrix initially in a are returned in b.
c     inv=0 means the square root of a is to be formed, and inv=1 means
c     that the inverse square root is to be formed.
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(mdim,mdim),b(mdim,mdim),aa(mdim),bb(mdim)
      dimension res(mdim,mdim)
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
      dimension ilifq(maxorb),iky(maxorb)
      data done/1.0d0/
c
c     call f02abf(a,mdim,nbas,aa,b,mdim,bb,ifail)
      l3 = mdim*mdim
      do i=1,mdim
         ilifq(i) = (i-1)*mdim
         iky(i) = i*(i-1)/2
      end do
      diaacc = 1.0d-11
      call trianc(a,res,mdim,mdim)
      write(iwr,*)'call jacobi'
      call jacobi(res,iky,mdim,b,ilifq,mdim,aa,2,2,diaacc)
      write(iwr,*)'return from jacobi'
c
      call vclr(res,1,l3)
      ignore = 0
      do i = 1 , nbas
         if (aa(i).le.1.0d-12) then
             aa(i) = 1.0d0
             ignore = ignore + 1
         end if
         term = dsqrt(aa(i))
         if (inv.ne.0) term = done/term
c
         call vsmul(b(1,i),1,term,res(1,i),1,nbas)
      enddo
c
      if (nignore.lt.ignore)
     +  call caserr('old negative eigen value in projection matrix.')
c
      call mxmg(res,1,mdim,b,mdim,1,a,1,mdim,nbas,nbas,nbas)
      return
      end
c
c  print a row of the energy table
c
      subroutine prtrow(s,e)
      implicit none
      real*8 e,fac1,fac2
      character s*(20)
      real*8 toang, ams
      integer ifau, nosymm, iseczz
      common/phycon/toang(30),ams(54),ifau,nosymm,iseczz
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c

      fac1 = toang(5)*toang(8)/1000.0d0
      fac2 = fac1 / toang(6)

      if(dabs(e).lt.1.0)then
         write(iwr,101)s,e,e*fac1,e*fac2
      else
         write(iwr,101)s,e
      endif
 101     format(1x,a20,f12.6,2f12.3)
      end
c
c  Morokuma scheme 2 modifications
c  - requires operator and fragment MO in AO basis
c
      subroutine modmoro(q, h, fq, fqi, l1)

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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)


      integer igmem_alloc
      external igmem_alloc

      real*8 h(*), fq(*), fqi(*), q(*)
      integer i, j, ij, ix, l2, l1

cccc      write(6,*)'mod key',imoro
      if(imoro .ne. ALL_BLOCKS)then

         l2 = (l1*(l1+1))/2
         ix = igmem_alloc(l2)

         call tran1eo(q, h, q(ix), fq, l1) 

cdbg
c         write(6,*)'before mod'
c         call prtril(q(ix),l1)

         ij=0
         do i=1,l1
            do j=1,i
               if (imoro .eq. ESX_CTT) then 

                  if( (ifrag(i) .eq. ifrag(j)) .and.
     &                 (ifocc(i) .eq. ifocc(j)))then
c     retain ESX block
                  else if ((ifrag(i) .ne. ifrag(j)) .and.
     &                    (ifocc(i) .ne .ifocc(j)))then
c     retain CT block
                  else
c     delete
                     q(ix+ij) = 0.0d0
                  endif
               else if(imoro.eq.ESX_PLX)then 
                  if(ifrag(i).ne.ifrag(j))then
                     q(ix+ij) = 0.0d0
                  endif
               else if(imoro.eq.ESX)then 
                  if((ifrag(i).ne.ifrag(j)) .or.
     &                 (ifocc(i) .ne .ifocc(j)) )then
                     q(ix+ij) = 0.0d0
                  endif
               else if(imoro.eq.ESX_EX)then 
                  if(ifocc(i) .ne .ifocc(j))then
                     q(ix+ij) = 0.0d0
                  endif
               else
                  call caserr('flag error in modmoro')
               endif
               ij = ij + 1
            enddo
         enddo
c
c  back transform 

         call tran1eo(q, q(ix), h, fqi, l1) 

         call gmem_free(ix)

      endif

      end
c     
c     Morokuma Energy Decomposition Analysis - Direct SCF implementation
c
      subroutine morokdscf(q)
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
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
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
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
      real*8 enrgy, egrad
      common/gms_funct/enrgy,egrad(3*maxat)
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
c
      integer lfile, lotape, liblk, llblk
      common /filel/ lfile,lotape(20),liblk(20),llblk(20)
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
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
      real*8 corev, array
      common/blkcore/corev(512),array(10)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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

c crystal field
      integer isecfx 
      logical ocryst
      common/xfield/isecfx,ocryst
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      character *8 title,guess
      common/restrz/title(12),guess
c
      integer iso
      common/scra/iso(mxshel,48)
c
      real*8 ptr, dtr, ftr, gtr
      common/junk/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
c
      real*8 st, cdiis, rdiis, derror, sdiis
      integer iposit,nstore,mp,junkj
      logical ondiis
      common/diisd/st(210),cdiis(20),rdiis(19),derror,sdiis(20),
     +             iposit(20),nstore,mp,ondiis,junkj

      integer ijkl
      common/craypk/ijkl(1360)
c

      real*8 q(*)

      character*30 filename
      logical o1e
      dimension o1e(6)
      character*8 zrhf
      character*4 yavr

      integer m1, m10, m16, igs
      integer l0, l1, l2, l3, l4, len2, lscdf, len
      integer i10,i20, i30, i21, ioc, ien
c     integer i40, i31
      integer i50, i60, iblk16

      integer ifq, ifo, ife, ifov
      integer ih0x, ihx

      integer nw1d, nw2d
      real*8 diaacc
      integer istat
      logical outon, odamph, oshift
      logical opr_skel, opr_symh, opr_vecs, opr_dens, opr_fock
      logical opr_conv, opr_parm
      integer ii, ij, j, k, icount, ix, ifr
      integer is, iff, lockt, lprnt, loop
      integer nshblk, nshtri, len171, ln171

      real*8 cpulft
      external cpulft

      real*8 chg_fa, coor_fa(3)
      real*8 ehf_frag(2), efrag
      integer nat_f1, num_f1
      integer nocmx
      logical oerr


c, jblkqa
c  jblkpa, jblkea
      integer jblkh, jblkh0, jblkd
c  jblkqs
c jblkst, jblks, jblkf
      integer iblkh0, iblkqq
      integer iblkh

      integer iov, iovt, ike, ih0, ip
      integer iv, ivsave, ivt
      integer itmp

      real*8 timeit
      integer i,ibase, lwor,last
      integer lword, nav, ndaf

      integer irdmat,iprefa,ifock,idmat

c     integer mpunch
      integer lensec

      integer igmem_alloc, igmem_max_memory
      external igmem_alloc, igmem_max_memory
      
      integer lenwrd
      external lenwrd
      
      logical ochkfa
      external ochkfa

      real*8 eel0, eel, ees, eesx, exprime, epl, eex, eint, ect
      real*8 emix
c     real*8 eexpl, eplx

      real*8 dum

      real*8 enucf
      external enucf

      data zrhf /'rhf'/

      data m1,m10,m16/1,10,16/

      data igs/5/

      call cpuwal(begin,ebegin)
      outon = nprint .ne. -5
      if(outon)write(iwr,180)
 180  format(//1x,104('-'))

      nav = lenwrd()


c
c     ----- read in transformation matrices for s,p,d,f,g basis functions.
c     
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
      if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
      if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
      call readi(iso,nw196(5)*nav,ibl196(5),idaf)

      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = l2+l2
      len2=l2+l2

      iov = igmem_alloc(l2)     ! S
      ih0 = igmem_alloc(l2)     ! core hamiltonian
      ih0x = igmem_alloc(l2)    ! modified core hamiltonian
      ihx = igmem_alloc(l2)     ! modified 2e hamiltonian
      iovt = igmem_alloc(l2)    ! S (symmetry adapted basis)
      ike = igmem_alloc(l2)
c
      do i =1,3
       o1e(i) = .true.
       o1e(i+3) = .false.
      enddo
      call getmat(q(iov),q(ike),q(ih0),
     +     q(iov),q(iov),q(iov),
     +     array,num,o1e,ionsec)
       call wrt3(q(iov),l2,ibl7s,num8)
       call wrt3(q(ike),l2,ibl7t,num8)
       call wrt3(q(ih0),l2,ibl7f,num8)
c
c     ----- transform the s-matrix
c     
      call tranp(q(iov),q(iovt))
      call wrt3(q(iovt),l2,ibl7st,num8)
c
c     ----- establish buffers for the two electron integral file(s)
c     
      if(nopk .ne. 1)
     &     call caserr('internal integral storage error')
      if (zscftp .ne. zrhf) call caserr('bad scf type')
      if (na .ne. nb) call caserr('must be closed shell')
c
      if(.not.odscf) then
      lfile=m2file
      do 141 i=1,lfile
         lotape(i)=m2tape(i)
         liblk(i)  =m2blk(i)
         llblk(i)  =liblk(i)-m2last(i)
 141  continue
      endif
c     
c     clear buffer for integral output
c     
c     call icopy(1360,0,0,ijkl,1)
c     
c     atomic density startup for open and closed shells
c     modified for correct restart
c     
      if (guess.eq.'atoms') then
c     
c     do one cycle scf with the density matrix from denat
c     
         call denscf(q,zscftp)
c         if (irest.ne.0) go to 160
c     
c     set zguess to 'anything' so denscf will be called but once
c     
         zguess = 'anything'
         guess  = 'anything'
      end if
c     
c     ----- closed shell restricted hartree fock -----
c     
c     determine available memory
      lwor = igmem_max_memory()
c
c  attempt to miminimise original over-generous estimate
c
      
      lword = 10*l2+ 2 * ikyp(nshell)
      if (lwor.lt.lword) then
         call caserr('not enough memory')
      endif
      ibase = igmem_alloc(lword)
c
c  start of rhf driver
c
c     out = nprint .eq. 5
      timeit = cpulft(1)
c     
c     ----- set pointers for partitioning of core -----
c     
      nw1d = 16+4/nav
      nw2d = 270 + 24/nav
c     allow for reduced density matrices etc.
c
      irdmat = ibase
      iprefa=irdmat+ikyp(nshell)
      i50 = iprefa+ikyp(nshell)
      i60 = i50+l2
      i10 = i60+l2
      ifock = i10
      i20 = i10+l2
      i21 = i10+l3
      i30 = i20+l2
      jblkh0= i30+l2
      jblkh = jblkh0+l2
      jblkd = jblkh+l2   !symass + extrpmm
      ii    = jblkd+3*l2

      last = ii-ibase
      if(last.gt.lword)then
         write(iwr,9309)lword,last
         call caserr('insufficient memory available')
      endif

      lprnt = l1  ! print all virtuals

      dlnmxd=0.0d0
      dlntol=tolitr(3)
      call start_time_period(TP_RDMAT)
c
c     specification of section isect(471)
c     only output delta-difference data if is this is to be invoked
c     i.e. only rdmat and prefac traffic involved in default
c *** section isect(471) now holds:
c *** prefac mat., reduced dmat, delta f, delta d, last f, current d
c
      m171t=171
      len171=  lensec(l2)
      nshtri=ikyp(nshell)
      nshblk=lensec(nshtri)
      iof171(1)=0
      iof171(2)= iof171(1) + nshblk
      iof171(3)= iof171(2) + nshblk
      do i=4,6
      iof171(i)=iof171(i-1)+len171
      enddo
      call rdmake(q(iprefa))
      ln171=nshblk*2
      if(odelta) then
       ln171=4*len171+ln171
      endif
      call secput(isect(471),m171t,ln171,ibl171)
      if(outon) then
      write(iwr,2980) ibl171,ln171
      endif
      call wrt3(q(iprefa),nshtri,ibl171,idaf)
      if(odelta) then
        call zer171(q(i10),l2,4,ibl171+iof171(3),len171,idaf)
      endif
      call clredx
c
      call end_time_period(TP_RDMAT)
c
c  dynamic core
c
      ioc = igmem_alloc(l1)
      ien = igmem_alloc(l1)
      ip  = igmem_alloc(l2) 
      idmat = ip
      iv  = igmem_alloc(l3) 
      ivt  = igmem_alloc(l3) 
      ivsave  = igmem_alloc(l3) 
c
c     ----- occupation numbers -----
c     
c     -o- at q(ioc)

c     
c     ----- initialize variables -----
c     
      if (maxcyc .le. 0) maxcyc = 30
c     mpunch = 1
c     if (npunch .ne. 1) mpunch = 0
c     if (nprint .eq. 7) mpunch = 1
      i=lensec(l2)

      if(numdis.eq.0) numdis=num8
c
c  some of these are now reset within iterate
c  (shambles)
c
      if(irest.ne.3.or.numdis.eq.num8) then
         ehf = 0.0d0
         ehf0 = 0.0d0
         iter = 0
         kcount = 0
         iterv = 0
         damp = 0.0d0
         damp0 = 0.0d0
         rshift = 0.0d0
         diff = 0.0d0
         diffd = 0.0d0
         diffp = 0.0d0
         diffpp = 0.0d0
         de = 0.0d0
c        dep = 0.0d0
         deavg = 0.0d0
      else
         call rdedx(en,nw1d,ibldis,numdis)
         call reads(st,nw2d,numdis)
      endif
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0

c     oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
c     
      if (odamph .or. oshift) damp = 1.0d0
c     

      if(outon)write (iwr,9348)
 9348 format(/
     +     40x,43('*') /
     +     40x,'* closed-shell rhf direct-scf calculation *'/
     +     40x,'*  morokuma energy decomposition version  *'/
     +     40x,43('*'))

c
c  allocate additional pointers
c
c     ifq (fragment vectors)
c     ifo (fragment occupations)
c     ife (fragment orbital energies)
c     ifs (fragment overlap matrix - used for integrity test only)
c
      ifq = igmem_alloc(l3)
      ifo = igmem_alloc(l1)
      ife = igmem_alloc(l1)
      ifov = igmem_alloc(l2)
c
c
      imoro = ALL_BLOCKS

c
c print type of run,and load info from previous runs
c
      if(imode.eq.1)then
       write(iwr,*)'first fragment calculation fragment name=',
     &        fragname(1)
      else if(imode.eq.2)then
       write(iwr,*)'second fragment calculation fragment name=',
     &        fragname(2)
      else if(imode.eq.3)then
       write(iwr,*)'perform interaction energy analysis'
c
c  check fragment geometries and import vectors
c
         icount = 1
         ix = 0

         call dcopy(l3, 0.0d0, 0, q(ifq), 1)
         call dcopy(l2, 0.0d0, 0, q(ifov), 1)

         do ifr = 1,2

            write(iwr,*)'fragment ',ifr,' = ',fragname(ifr)
            
            call charlim(fragname(ifr),is,iff)
            filename=fragname(ifr)(is:iff)
            open(unit=57,file=filename,form='formatted',
     +       status='unknown')

            read(57,101,end=998)nat_f1
            do i=1, nat_f1
               read(57,103,end=998)chg_fa,coor_fa
               if(.not.ochkfa(  czan(icount), c(1,icount),
     &              chg_fa, coor_fa))then

                  write(iwr,*)'frag',ifr
                  write(iwr,*)'frag atom',i
                  write(iwr,*)'molec atom',icount
                  call caserr('atom mismatch')
               endif
               icount = icount + 1
            enddo

c
c  load fragment wavefunction
c
            read(57,101,end=998)num_f1
            do i=1,num_f1
               do j=1,num_f1
                  read(57,102,end=998)
     &                 q(ifq+(ix+i-1)*num+(ix+j-1))
               enddo
            enddo


c orbital energies
            read(57,101,end=998)num_f1
            do j=1,num_f1
               read(57,102,end=998)q(ife + (ix+j-1))
            enddo

c occupations..
            read(57,101,end=998)num_f1
            do j=1,num_f1
               read(57,102,end=998)dum
               q(ifo + (ix+j-1)) = dum

               if(dabs(dum) .lt. 1.0d-6)then
                  ifocc(ix+j) = 0
               else if(dabs(dum - 2.0d0) .lt. 1.0d-6)then
                  ifocc(ix+j) = 1
               else
                  call caserr('fragments must be closed shell')
               endif

            enddo
c
c total energy
            read(57,101,end=998)num_f1
            read(57,102,end=998)ehf_frag(ifr)
c
c fragment overlap matrix
            read(57,101,end=998)num_f1
            ij=(ix*(ix+1))/2 
            do i=1,num_f1
               ij=ij+ix
               do j=1,i
                  read(57,102,end=998)q(ifov +ij)
                  ij=ij+1
               enddo
            enddo

c set ifrag elements
            do j=1,num_f1
               ifrag(ix+j) = ifr
            enddo

            ix = ix + num_f1
            close(57)

         enddo


         write(iwr,*)'ifrag array',(ifrag(j),j=1,num)


         if(icount .ne. nat +1)call caserr('atom count error')
         if(ix .ne. num)call caserr('basis fn count error')

         write(iwr,*)'fragment atom lists were ok'

         write (iwr,*)'Fragment MO matrix'
         call prsq(q(ifq), l1, l1, l1)

         len = (l1-1)/lenwrd() + 1


c         write (iwr,*)'Fragment Overlap matrix'
c         call prtril(q(ifov),l1)

         itmp = igmem_alloc(l2)

         call dcopy(l2,q(iov),1,q(itmp),1)
         call ifzero(q(itmp),l1)

c         write (iwr,*)'Fragment Overlap matrix - from molecular '
c         call prtril(q(itmp),l1)

         oerr = .false.
         ij = 0
         do i = 1, l1
            do j = 1, i
               if(dabs(q(itmp + ij) 
     &              - q(ifov + ij)) .gt. 1.0d-8)then
                  write(iwr,*)i, j, q(itmp + ij), 
     &                 q(ifov + ij)
                  oerr = .true.
               endif
               ij = ij + 1
            enddo
         enddo

         if(oerr)then
            write(iwr,*)
     &           'Intrafragment overlap does not match that from'
            write(iwr,*)
     &           'fragments. Check that basis functions are in the'
            write(iwr,*)
     &           'same order in both fragments and complex.'
            call caserr('mismatched fragment/molecular overlaps')
         endif

         call gmem_free(itmp)

      else
         write(iwr,*)'test: mode flag=',imode
      endif

c     
c     ----- nuclear energy
c     
c  ocryst: en = crnuc(nat,czan,c,ecrn)
      en = enucf(nat,czan,c)
      write (iwr,9008) en
 9008 format(/,30h ----- nuclear energy ----- = ,f20.12)


c
c  some useful matrices for all cases
c
c
c  retain all virtuals
c
      l0=num

      if(imode.le.2)then
c
c   perform an SCF calculation on the fragment and save the results

c     
c restore vectors from guess (SAB)
c
         call rdedx(q(ivt),l3,ibl3qa,idaf)
c
c read orbital energies -> q(ien) and compute occupancies 
c
         call rdedx(q(ien),l1,ibl3ea,idaf)
         call llvmo(q(ien),q(ioc),na,nocmx,l1)
         call dscal(l1,2.0d0,q(ioc),1)
c
c      yavr=yblnk
c      if(nocmx.gt.na)yavr=yav
       lprnt = l1
c         if(.not.oprint(20)) lprnt=min(nocmx+5,l1)

c
c  full 2e hamiltonian
         oifc = .true.
         oife = .true.
         oifh = .true.
c
c no printing
         opr_skel  = .false.
         opr_symh  = .false.
         opr_fock  = .false.
         opr_vecs  = .false.
         opr_dens  = .false.
         opr_conv  = .true.
         opr_parm  = .true.
         lockt=lock

         call diterate(q,
     &        q(ih0), q(i10), q(ip), q(irdmat), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        q(jblkd),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstore, derror, iposit, diaacc, lprnt,
     &        num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .false.,
     &        istat,
     &        nshtri,
     &        zscftp)

         lock = lockt

c
c save molecular solution on dumpfile

         if(numdis.eq.num8)numdis=0
         ndaf = mouta
         call scfsav(q(ivt),q(ip),
     *        q(ien),q(ioc),ndaf,l1,
     *        l2,ibl3pa,ibl3ea)

         call dcopy(l3,q(ivt),1,q(i10),1)
         call dcopy(l1,q(ien),1,q(i21),1)

         call analmo(q(i10),q(i21),q(ioc),ilifq,l0,l1)
         call tdown(q(i10),ilifq,q(i10),ilifq,l0)

         if(otsym)call symass(q(i10),q(i21),q(ioc),q)

         if(.not. oprint(25))then
            write (iwr,9148)
 9148       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors'/
     +           50x,12('-'))
            call prev(q(i10),q(i21),lprnt,l1,l1)
c     
            if(.not.outon)then
               write (iwr,9168)
 9168       format(/20x,14('*')/20x,'density matrix'/20x,14('*'))
               call prtril(q(i30),l1)
            endif
         endif

c        if(istat .ne. 0)call caserr('fatal SCF convergence problem')
c        above brain dead approach will not do..

         array(1) = en
         array(2) = ehf
         array(3) = etot
         do loop = 4,10
          array(loop) = 0.0d0
         enddo
         call secput(isect(494),m16,m1,iblk16)
         call wrt3(array,m10,iblk16,idaf)
c
         if(istat .ne. 0) go to 8000
c
c  write results out for use in subsequent interaction analysis run
c
         call charlim(fragname(imode),is,iff)
         filename=fragname(imode)(is:iff)
         write(iwr,*)'*** saving vectors to ',filename
         open(unit=57,file=filename,form='formatted',status='unknown')
         write(57,101)nat
         do i=1,nat
            write(57,103)czan(i),(c(j,i),j=1,3)
         enddo

         write(57,101)num
c
c !!! should be in AO basis here
c
         write(57,102)(q(ivt + k - 1),k=1,l3)
         write(57,101)num
         write(57,102)(q(ien + k - 1),k=1,l1)
         write(57,101)num
         write(57,102)(q(ioc + k - 1),k=1,l1)
         write(57,101)num
         write(57,102)etot
         write(57,101)num
         write(57,102)(q(iov + k -1 ),k=1,l2)

 101     format(1x,i3)
 102     format(1x,e20.14)
 103     format(1x,4e20.14)

         close(57)

      else if(imode.eq.3) then

         write(iwr,*)'sum of fragment energies:'
         write(iwr,*)ehf_frag(1) + ehf_frag(2)

         write(iwr,*)'use fragment MO to startup'

c    try othog in different order
c         call dcopy(l3,q(ifq),1,q(ivt),1)
c         write(iwr,*)'Schmidt 2'
c ????? need to convert to SAB here
c         call dcopy(l3,q(ifq),1,q(ivt),1)
c hardwire for this case
c         ilo = 14
c         ihi = 26
c         jlo = 1
c         jhi = 5
c         call schm(q,q(ivt),q(iovt),ilo,ihi,jlo,jhi,l1)
c         write(iwr,*)'check transformed s'
c         call tran1eo(q, q(iovt), q(i10), q(ivt), l1) 
c         call prtril(q(i10),l1)
c
c !! need to convert to SAB here
c
         call dcopy(l3,q(ifq),1,q(ivt),1)
c
c  use occ and energy from fragment calculations
c
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num
c
c  use full 2e hamiltonian
         oifc = .true.  
         oife = .true.
         oifh = .true.
c
c no printing
         opr_skel  = .false.
         opr_symh  = .false.
         opr_fock  = .false.
         opr_vecs  = .false.
         opr_dens  = .false.
         opr_conv  = .true.
         opr_parm  = .true.
         lockt=lock

         imoro = ALL_BLOCKS
         
         call diterate(q,
     &        q(ih0), q(i10), q(ip), q(irdmat), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        q(jblkd),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstore, derror, iposit, diaacc, lprnt,
     &        num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat,
     &        nshtri,
     &        zscftp)

         emorok(ALL_BLOCKS,YES_K) = eel + en

         emorok(ESX_EX,YES_K) = eel0 + en

c
c save molecular solution on dumpfile for possible further
c analysis
c
         if(numdis.eq.num8)numdis=0
         ndaf = mouta
         call scfsav(q(ivt),q(ip),
     *        q(ien),q(ioc),ndaf,l1,
     *        l2,ibl3pa,ibl3ea)

         array(1) = en
         array(2) = ehf
         array(3) = etot
         do loop = 4,10
          array(loop) = 0.0d0
         enddo
         call secput(isect(494),m16,m1,iblk16)
         call wrt3(array,m10,iblk16,idaf)
c
         if(istat .ne. 0) go to 8000
c

c ????? need to convert to SAB here

c
c  load fragment vectors, occupancies and energies
c  into the molecular arrays
c
         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oifc = .true.  
         oife = .true.         
         oifh = .true.
c flag retention of ESX (diagonal) blocks only
         imoro = ESX            

         write(iwr,*)'ESX'
         lockt=lock

         call diterate(q,
     &        q(ih0), q(i10), q(ip), q(irdmat), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        q(jblkd),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstore, derror, iposit, diaacc, lprnt,
     &        num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat,
     &        nshtri,
     &        zscftp)

         lock=lockt

         emorok(ESX,YES_K) = eel + en


         write(iwr,*)'ESX + PLX without exchange'

         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oife = .false. 
         maxcyc=30

c
         if(istat .ne. 0) go to 8000
c
         imoro = ESX_PLX

         call diterate(q,
     &        q(ih0), q(i10), q(ip), q(irdmat), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        q(jblkd),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstore, derror, iposit, diaacc, lprnt,
     &        num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat,
     &        nshtri,
     &        zscftp)

         emorok(ESX_PLX,NO_K) = eel + en

         emorok(ESX,NO_K)     = eel0 + en

         write(iwr,*)'ESX + CT'

         call dcopy(l3,q(ifq),1,q(ivt),1)
         call dcopy(l1,q(ifo),1,q(ioc),1)
         call dcopy(l1,q(ife),1,q(ien),1)
         nocmx=num

         oife = .true. 
         maxcyc=30

         imoro = ESX_CTT
         call diterate(q,
     &        q(ih0), q(i10), q(ip), q(irdmat), 
     &        q(ien), q(iovt),
     &        q(ivt)  , q(ioc),
     &        q(jblkh), q(ivsave), q(jblkh0),
     &        q(jblkd),
     &        gapa1, gapa2, ibrk, igs,
     &        opr_skel, opr_symh, opr_fock, opr_vecs, opr_dens,
     &        opr_parm, opr_conv,
     &        lock, ondiis, nstore, derror, iposit, diaacc, lprnt,
     &        num, iky, ilifq, nocmx, na, yavr, 
     &        eel0, eel,
     &        q(ifq),
     &        .true.,
     &        istat,
     &        nshtri,
     &        zscftp)
c
         if(istat .ne. 0) go to 8000
c

         emorok(ESX_CTT,YES_K) = eel + en

         write(iwr,*)
         write(iwr,*)'       Morokuma  Analysis'
         write(iwr,*)

         efrag = ehf_frag(1) + ehf_frag(2)

         write(iwr,*)'Results of calculations: '

         write(iwr,104)
         call prtrow('Fragment 1          ',ehf_frag(1))
         call prtrow('Fragment 2          ',ehf_frag(2))

         call prtrow(' complex            ',emorok(0,1) - efrag)

         call prtrow(' exs                ',emorok(1,1) - efrag)
         call prtrow(' exs + plx          ',emorok(2,1) - efrag)
         call prtrow(' exs + ct           ',emorok(3,1) - efrag)
         call prtrow(' exs + ex^          ',emorok(4,1) - efrag)

         call prtrow(' exs no k           ',emorok(1,2) - efrag)
         call prtrow(' exs + plx no k     ',emorok(2,2) - efrag)
         call prtrow(' exs + ct no k      ',emorok(3,2) - efrag)
         call prtrow(' exs + ex^ no k     ',emorok(4,2) - efrag)

         write(iwr,104)

         write(iwr,*)

         write(iwr,100)

         ees = emorok(ESX,NO_K) - efrag
         call prtrow('Electrostatic  (Ees)',ees)

         eesx = emorok(ESX,YES_K) - efrag

         exprime =  emorok(ESX_EX,YES_K) - emorok(ESX,YES_K)

         eex = exprime + eesx - ees
         call prtrow('Exchange       (Eex)',eex)

         epl = emorok(ESX_PLX,NO_K) - efrag - ees
         call prtrow('Polarisation   (Epl)',epl)

         ect = emorok(ESX_CTT,YES_K) -  emorok(ESX,YES_K)
         call prtrow('Charge Transfer(Ect)',ect)

c eexpl currently missed out
c         eexpl = eplx - epl
c         call prtrow('Exchange-Pol (Eexpl)',eexpl)

         eint = emorok(ALL_BLOCKS,YES_K) - efrag

c eexpl currently missed out
c         emix = eint - ( ees + epl + eex + ect + eexpl)

         emix = eint - ( ees + epl + eex + ect)
         call prtrow('Coupling Term (Emix)',emix)

         write(iwr,104)
         call prtrow('Total  Interaction  ',eint)
         write(iwr,104)

 100     format(//1x,
     &'      Analysis following Morokuma - IJQC v. X p. 325 (1976)',/,
     &'      Term                E(H)      E(Kj/mol)  E(kcal/mol)',
     &        /,1x,60('-'))
 104     format(1x,60('-'))

      endif

8000  continue
c
c   free dynamic core
c      
      call gmem_free(ifov)
      call gmem_free(ife)
      call gmem_free(ifo)
      call gmem_free(ifq)
      
      call gmem_free(ivsave)
      call gmem_free(ivt)
      call gmem_free(iv)
      call gmem_free(ip)
      call gmem_free(ien)
      call gmem_free(ioc)

      call gmem_free(ibase)
      
      call gmem_free(ike)
      call gmem_free(iovt)
      call gmem_free(ihx)
      call gmem_free(ih0x)
      call gmem_free(ih0)
      call gmem_free(iov)
      
      return

 998  write(iwr,*)'unexpected end of stream on ',filename
      call caserr('fragment file error')

c999  write(iwr,*)'problem with ',filename
c     call caserr('error opening fragment file')

 9309 format(//
     *     1x,'words available ',i8/
     *     1x,'words requested ',i8//)
c9177 format(10x,a36,f8.2)
c9178 format(//
c    *     10x,'********************************************'/
c    *     10x,'Scf timing statistics              (seconds)'/
c    *     10x,'********************************************'/
c    *     10x,'fock formation                      ',f8.2/
c    *     10x,'fock diagonalisation                ',f8.2/
c    *     10x,'diis                                ',f8.2/
c    *     10x,'fock transformation                 ',f8.2/
c    *     10x,'orthogonalisation                   ',f8.2/
c    *     10x,'symmetrisation, damping, extrpn.    ',f8.2/
c    *     10x,'density matrix, etc.                ',f8.2/
c    *     10x,'initialisation                      ',f8.2/
c    *     10x,'dft                                 ',f8.2/
c    *     10x,'********************************************'/
c    *     10x,'Total                               ',f8.2/
c    *     10x,'********************************************'/)
 2980 format(/1x,
     &  'output section 171 to block ',i5/1x,
     &  'section length              ',i5/)
      end


      subroutine diterate(q,
     &     h0,
     &     h  ,
     &     p  ,
     &     rdmat ,
     &     e  ,
     &     st  ,
     &     v  ,
     &     occ  ,
     &     hsav  ,
     &     vsav  ,  
     &     h00  ,
     &     extrap ,
     &     gapa1, gapa2,
     &     ibrk,
     &     igs,
     &     pr_skel  ,
     &     pr_symh  ,
     &     pr_fock  ,
     &     pr_vecs  ,
     &     pr_dens  ,
     &     pr_parm  ,
     &     pr_conv  ,
     &     lock  ,
     &     ondiis  ,
     &     nstore  ,
     &     derror  ,
     &     iposit, 
     &     diaacc,
     &     lprnt  ,
     &     num  ,
     &     iky,
     &     ilifq,
     &     nocmx  ,
     &     na  ,
     &     yavr,
     &     eel0,
     &     eel,
     &     fq,
     &     lock_occupations,
     &     istat,
     &     nshtri,
     &     zscftp)
c     
      implicit none
c     
c     arguments

      real*8 q(*)                 ! base memory
      real*8 h0(*)                ! core hamiltonian
c     logical oifc, oife, oifh  ! interfragment interaction flags
      real*8 h(*)                 ! full hamiltonian
      real*8 p(*)                 ! density
      real*8 rdmat(*)             ! reduced density
      real*8 e(*)                 ! orbital energies
      real*8 st(*)                ! overlap q(jblkst)
      real*8 v(*)                 ! vectors  ????
      real*8 occ(*)               ! occupancies
      real*8 hsav(*)              ! saved hamiltonian
      real*8 vsav(*)              ! vectors q(jblkqa)
      real*8 h00(*)               ! q(jblkh0) - used in diis
      real*8 extrap(*)            ! q(jblkd) - used in extrpmm
      real*8 gapa1, gapa2         ! level shifter
      integer ibrk              ! when to switch level shifters
      integer igs               ! how often to reorthogonalise
      logical pr_skel           ! print skeleton fock
      logical pr_symh           ! print symmetried fock
      logical pr_fock           ! print total fock
      logical pr_vecs           ! print eigenvectors
      logical pr_dens           ! print eigenvectors
      logical pr_parm           ! print scf parameters
      logical pr_conv           ! print convergence data
      integer lock              ! 
      logical ondiis            ! diis common /junk/
      integer nstore            ! number of diis vectors in use
      real*8 derror               ! diis error
      integer iposit(*)         ! addresses of diis fock matrices
      real*8 diaacc

      integer lprnt             ! number of orbitals to print
      integer nshtri            ! triangle size over shells
      integer num               ! size of basis
      integer iky(*)
      integer ilifq(*)          ! transformation table
      integer nocmx             ! number of occ orbs
      integer na                ! number of alpha electrons
c to flag partial occ of degenerate homos
      character*4 yavr         
c electronic energy of input orbitals
      real*8 eel0                 
c electronic energy of converged solution
      real*8 eel                 

      real*8 fq(*)                ! fragment orbitals
c force FMO occupation pattern on molecule
      logical lock_occupations  
      integer istat             ! return code
      character *8 zscftp       ! scftype
c     
c     local variables
c
      logical finished
      integer iscr              ! scratch memory pointer
      integer iscr1, iscr2      ! scratch memory pointer
      integer iscr1x, iscr2x, iscr3x ! scratch memory pointer
c     integer ifsave            ! scratch memory pointer
      integer iqdiis            ! diis array
      integer iht, i30

      integer ifqi, ihx, ih0x, istx, iovq

      integer l0, l1, l2, l3
      integer l4
      integer i, j, ii, ij, loop
      integer ig, m2, len

      real*8 diffdp,diffdpp
      real*8 dep
      real*8 dmptlc, dmplim
      real*8 ehf1, ehf2
      real*8 tim1
      real*8 skale, tlefti, tlefts
      real*8 deter

      logical oshift
      logical odamph
      logical oextra
      logical ocvged
      logical outon

      integer lockt

      logical odbg
      logical omemck

      real*8 dum

      real*8 qmax
      integer iforb, ifpvi, ivsav
      logical assign
      logical first_cycle
      logical force_serial
 
      integer iblkh0, iblkqq, iblkh, ndafd, ndafdi
      integer domo,vmo
c     
c dynamic core
c     
c     external functions
c
      external igmem_alloc
      integer igmem_alloc

      integer lensec
      external lensec

      external tracep
      real*8 tracep

      external idamax
      integer idamax
      external lenwrd
      integer lenwrd
      external ipg_nodeid
      integer ipg_nodeid

      real*8 cpulft
      external cpulft

c
c convergence control parameters
c
c
      logical odiis, optester, odynamic, omaxcyc,onoor,otestd
      integer maxcyc,mconv,nconv,npunch,icoupl,ifolow,irotb
      integer iter,kcount,iextin,iterv,idiisf
      real*8 accdi1,accdi2,dmpcut,acurcy,en,etot
      real*8 ehf,ehf0,diff,rshift,exttol,dmptol,vshtol
      real*8 damp,damp0,diffd,diffp,de,deavg,diffsp
      real*8 ek, vir,diffpp
      common/scfopt/maxcyc,mconv,nconv,npunch,accdi1,accdi2,odiis,
     +      icoupl,ifolow,irotb,dmpcut,acurcy,en,etot,ehf,ehf0,diff,
     +      iter,kcount,rshift,exttol,dmptol,vshtol,iextin,
     +      iterv,damp,damp0,diffd,diffp,diffpp,de,deavg,diffsp,
     +      ek,vir,idiisf,optester,odynamic,omaxcyc,onoor,otestd
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
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
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
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      integer select(maxorb), nsel
      integer ioccsav(maxorb)

      character*4 yblnk, yav

      data yblnk,yav/' ',' av'/

      diffdpp = 1.0d4
      diffdp = 1.0d4

      l0=num                    ! number of canonical orbitals kept
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      l4 = l2+l2

      outon = nprint .ne. -5
      odbg = .false.
      omemck = .false.

      dmptlc =1.0d-2

      lockt=lock

      istat=0


      call timrem(tlefts)

c
c     ----- decode convergence method -----
c           mconv is in /scfopt/
c
      oextra = mod(mconv,2) .eq. 0
      odamph = mod(mconv,4) .gt. 1
      oshift = mod(mconv,8) .gt. 3
c
c    ---- some tolerancies
c
      iextin = 4
      exttol = 1.0d-3
      dmptol = 1.0d-4
      vshtol = 0.4d0
      if(nconv.le.0)nconv=5
      acurcy=10.0d0**(-nconv)
c
c     is delta density ever to be invoked? if not, we can
c     reduce i/o activity around section 471 ..
c
c     odelta = deltol .gt. dlog(acurcy)
c
      write(iwr,*) 'odelta = ', odelta
c
c  old fock matrices for extrap
c     ifsave = igmem_alloc(3*l2 + 1000)

c
c flag for FMO-based locking

      first_cycle = .true.

c     if(odiis)then
      i=lensec(l2)
      iblkh0 = ibl7la
      iblkqq = iblkh0 + i
      iblkh=iblkqq+lensec(l3)
      if(numdis.eq.0.or.numdis.eq.num8)then
         numdis=num8
         ndafd=iblkh+i
      else
         ndafd=ibldis+2
      endif
      ndafdi=ndafd+3*i

c     if(odiis)then
c
c old fock matrices for diis 
c
c        iqdiis = igmem_alloc(16*l2 + 1000)
c        ii=1
c        do 180 loop=1,8
c           iposit(loop)=ii
c           ii=ii+l2+l2
c180     continue
c     endif
c
c  inverse eigenvector matrix for backtransformations
c
      ifqi = igmem_alloc(l3)

c
c old vectors for locking
c
      ivsav = igmem_alloc(l3)

      len = (l1-1)/lenwrd() + 1

      iscr1 = igmem_alloc(len)
      iscr2 = igmem_alloc(len)

      call dcopy(l1*l1,fq,1,q(ifqi),1)

      if(imode.gt.2)then
         call minvrt(q(ifqi),l1,deter,
     &        q(iscr1),q(iscr2))
      endif

      write(iwr,*)'deter',deter

      call gmem_free(iscr2)
      call gmem_free(iscr1)
c
c  these probably should be l2, but for mult2 in diis
c
      iscr1x = igmem_alloc(l4) 
      iscr2x = igmem_alloc(l4) 
      iht = igmem_alloc(l3)     ! transformed H
      i30 = igmem_alloc(l3)     ! scratch

      write(iwr,*)'copy ie'
c
c  local copies of S, and H0
c
      ih0x = igmem_alloc(l2) 
      call dcopy(l2,h0,1,q(ih0x),1)

      ihx = igmem_alloc(l2)
cccc      call dcopy(l2,h,1,q(ihx),1)

      istx = igmem_alloc(l2)
      call dcopy(l2,st,1,q(istx),1)

c      write(iwr,*)'mod 1e'
c
c  modify S, H0 according to Morokuma scheme
c

      call modmoro( q, q(istx), fq, q(ifqi), l1)
      call modmoro( q, q(ih0x), fq, q(ifqi), l1)

cdbg
c      write(iwr,*)'modified overlap'
c      call prtril(q(istx),l1)
c      onto num8 for diiscd
       call wrt3(q(istx),l2,ibl7s,num8)
cdbg
c      write(iwr,*)'modified h0'
c      call prtril(q(ih0x),l1)
c
c  qq(ivoff+) eigenvectors
c
      write(iwr,*)'diag q'
      iovq = igmem_alloc(l3)

      call dcopy(l2,st,1,q(iscr1x),1)

      iscr3x = igmem_alloc(max(l3,l1+1))
      call qmat(q,q(iscr1x),q(iovq),q(iscr2x),
     &     q(iscr3x),iky,l0,l1,l3,l1,.false.)
      call gmem_free(iscr3x)

      write(iwr,*)'symmetrisation'
c
c  orthogonalise input orbitals according to modified S
c  the symmetric scheme
c

      nsel = 0
      do i = 1,l1

         dum = occ(i)
         if(dabs(dum) .lt. 1.0d-6)then
            select(i) = 0
         else if(dabs(dum - 2.0d0) .lt. 1.0d-6)then
            select(i) = 1
            nsel = nsel + 1
         else
            call caserr('systems must be closed shell')
         endif

      enddo

c      write(iwr,*)'initial check s'
c      call tran1eo(q(istx), q(iscr1x), v, l1) 
c      call prtril(q(iscr1x),l1)

      call dcopy(l3,v,1,q(iscr1x),1)
      call symmorth2(q,q(iscr1x), v, 
     &     select,q(istx),l1,nsel,iwr)

c      write(iwr,*)'check transformed s'
c      call tran1eo(q(istx), q(iscr1x), v, l1) 
c      call prtril(q(iscr1x),l1)

c      write(iwr,*)'input orbs'
c      call prsq(v,l1,l1,l1)
c
c  - initial density matrix  !!! symmetry adaption
c
c
c transform vectors (SA basis) ->  (AO basis)
c
c         call tdown(q(iv),ilifq,q(ivt),ilifq,l1)

      call dmtx(p,v,occ,iky,nocmx,l1,l1)
c     write(iwr,*) 'sumup: dmtxp'
c     call sumup(p,l2,iwr)
cc      write(iwr,*)'initial density'
cc      call prtril(p,l1)

      finished = .false. 

      if(pr_parm)then
         write(iwr,*)'SCF parameters'
         write(iwr,*)'maxcyc =',maxcyc
         write(iwr,*)'oifc   =',oifc
         write(iwr,*)'oife   =',oife
         write(iwr,*)'oifh   =',oifh
         write(iwr,*)'odiis  =',odiis
         write(iwr,*)'accdi1 =',accdi1
         write(iwr,*)'acurcy =',acurcy
         write(iwr,*)'dmpcut =',dmpcut
         write(iwr,*)'oshift =',oshift
         write(iwr,*)'odamph =',odamph
         write(iwr,*)'oextra =',oextra
         write(iwr,*)'shifters',gapa1,gapa2,ibrk
      endif

      lockt=lock
c
c     if(pr_conv)then
c        write (iwr,9028)
c9028    format(/1x,100('=')/
c    *        3x,'cycle',10x,'total',5x,'electronic',8x,'e conv.',
c    *        9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
c    *        17x,'energy',9x,'energy',35x,'shift'/1x,100('='))
c     endif
c
c   some iteration control parameters (see also calling program)
c
c     
c     ----- set convergence criteria -----
c     - see also subroutine iterate
c

c diis
      nstore= 0
      accdi2=acurcy*acurcy
      ondiis=.false.

      ehf = 0.0d0
      ehf0 = 0.0d0
      iter = 0
      kcount = 0
      iterv = 0
      damp = 0.0d0
      damp0 = 0.0d0
      rshift = 0.0d0
      diff = 0.0d0
      diffd = 0.0d0
      diffp = 0.0d0
      diffpp = 0.0d0
      de = 0.0d0
      dep = 0.0d0
      deavg = 0.0d0
      if (dmpcut .le. 0.0d0) dmpcut = 0.0d0
      if (odamph .or. oshift) damp = 1.0d0
c
      call dcopy(l3, v, 1, vsav, 1)
      call wrt3(vsav,l3,ibl3qa,idaf)
c
      do while (.not. finished) 
c

         if(omemck)call ma_summarize_allocated_blocks
c ***
c *** now form increment to density matrix
c *** in section 171 is kept:
c *** delta f, delta d, last f, current d, reduced dmat, prefac mat.
c *** note that delta f and delta d are NOT accessed if
c *** delta-SCF has not been requested : odelta = .false.)
c
      call start_time_period(TP_RDMAT)
      if(outon) write(iwr,3030)
       if (odelta) then
        call wrt3(p,l2,ibl171+iof171(6),idaf)
       endif
       if( (dlnmxd.gt.deltol) ) then
       if(outon) write(iwr,3000)
        call mkrdmt(zscftp,rdmat,p,l2,nprint)
        call clredx
       endif
       call wrt3(rdmat,nshtri,ibl171+iof171(2),idaf)
c
      if(iter.lt.itrtol(1)) then
         dlntol=tolitr(1)
      else if(iter.lt.itrtol(2)) then
         dlntol=tolitr(2)
      else
         dlntol=tolitr(3)
      endif
      call end_time_period(TP_RDMAT)
      dlntol=dlntol - dmin1(dmax1(dlnmxd,delfac),0.0d0)
      if(outon) then
         write(iwr,3020) dlntol
      endif
c
c   construct skeletonised 2e part of fock matrix
c
         call dmhstar(q, p, h, rdmat, istat)
c

         if(istat.ne.0) go to 460
         if(odelta) then
          call wrt3(h,l2,ibl171+iof171(5),idaf)
         endif

         if (pr_skel) then
            write (iwr,9068)
 9068       format(/
     +           20x,20('*')/
     +           20x,'skeleton fock matrix'/
     +           20x,20('*'))
            call prtril(h,l1)
         endif
c     
c     ----- symmetrize skeleton fock matrix -----
c     
         call symh(h,q(iscr1x),iky,0,0)
c
c  remove inter-fragment terms from h in AO basis
c
         if(.not. oifh)call ifzero(h,l1)
c
c   scheme 2 more comprehensive tinkering
c
         call modmoro(q,  h, fq, q(ifqi), l1)

         if (pr_symh) then
            write (iwr,9088)
 9088       format(/
     +           20x,23('*')/
     +           20x,'symmetrized fock matrix'/
     +           20x,23('-'))
            call prtril(h,l1)
         endif
c     
c     ----- read in core hamiltonian matrix
c     and calculate hf energy -----
c     
c     for xtal field calculations, the core hamiltonian
c     includes the crystal potential, but for the energy 
c     all terms modelling molecule-molecule interactions 
c     must be divided by two 
c     

         ehf0 = ehf
         call vadd(h,1,q(ih0x),1,h,1,l2)
c        call adonee(h,q(ih0x),q(iprefa))
         ehf1 = tracep(p,q(ih0x),l1)
         ehf2 = tracep(p,h,l1)
         ehf = (ehf1+ehf2) * 0.5d0

         if(iter.eq.0)eel0 = ehf
c     
c     if(ocryst)then
c     call secget(isecfx,m53,iblk)
c     call rdedx(q(jblkfx),l2,iblk,idaf)
c     ehf3 = tracep(p,q(jblkfx),l1)
c     write(iwr,*)'ehf3',ehf3
c     esave = ehf
c     ehf = (ehf1+ehf2-ehf3)*pt5
c     ecre = esave - ehf
c     write(iwr,*)'new ehf, ecre',ehf,ecre
c     endif
c     
c     ----- save fock matrix h -> hsav
c     
         call dcopy(l2,h,1,hsav,1)

         if (pr_fock) then
            write (iwr,9048)
 9048       format(
     +           20x,11('*'),/
     +           20x,'fock matrix'/
     +           20x,11('*'))
            call prtril(h,l1)
         endif
c     
         if(maxcyc .eq. 0)then
            finished = .true.
            goto 400
         endif

         iter = iter+1
         etot = ehf+en
         dep = de
         de = ehf-ehf0
         if (iter .eq. 1) then
            deavg = 0.0d0
         else if (iter .eq. 2) then
            deavg =  dabs(de)
         else if (iter .ge. 3) then
            deavg = (  dabs(de)+  dabs(dep)+0.2d0*deavg)/2.2d0
         endif
c     
c     ----- damp and extrapolate hamiltonian matrix -----
c     
c     -h - at x(i10)     hamiltonian matrix (n th iteration)
c     -ho- at x(i20)     old h matrix (n-1 th iteration)
c     -ha- at x(i30)     ancient h matrix (n-2 th iteration)
c     -hp- at x(i40)     prehistoric h matrix (n-3 th iteration)
c     
         if (iter .gt. 2) call dampd(de,dep,deavg,damp,
     +        acurcy,diff,diffp,dmptlc)
c     
         if (damp .lt. dmpcut) damp = dmpcut

c         write(iwr,*)'iter, iterv',iter, iterv
c         write(iwr,*)'extrpm',de, damp, damp0

         if(omemck)write(iwr,*)'2'
         if(omemck)call ma_summarize_allocated_blocks

      iscr3x = igmem_alloc(l3)
         call extrpmm(de,damp,damp0,h,
     &        q(iscr1x),q(iscr2x),q(iscr3x),
     &        l1, l2, iterv,1,1,extrap)
      call gmem_free(iscr3x)

         diffpp=diffp
         diffp=diff
         diffdpp=diffdp
         diffdp=diffd

         if(omemck)write(iwr,*)'2.5'
         if(omemck)call ma_summarize_allocated_blocks

c         write(iwr,*)'h after extrap'
c         call prtril(h,l1)
c         write(iwr,*)'odiis',odiis
c
c    takes over scratch arrays from extrpmm
c  - should replace these when extrpmm loses
c    its l2*3 old matrices
c
c    h is copied to n**2 space to allow for
c    mult2 workspace
c
      if(odiis) then
      call start_time_period(TP_DIIS)
c
c   st should be in non-SAB basis here!!!

c           call dcopy(l2,h,1,q(iscr1x),1)
c     iscr3x = igmem_alloc(l3)
c           call diiscmm(q(iqdiis),q(iscr1x),
c    &           q(iscr2x), q(iscr1x),q(iscr3x),
c    &           h00,vsav,q(istx),diff)
c            write(iwr,*)'diis res',ondiis,diff
c     call gmem_free(iscr3x)

      call dcopy(l2,h,1,q(iscr1x),1)
      iscr3x = igmem_alloc(l4)
      call diiscd(q(iscr1x),q(iscr2x),q(iscr1x),q(iscr3x),
     +        iblkh0,ibl3qa,
     +        ndafdi,numdis,diff,lockt,domo,vmo)
      call gmem_free(iscr3x)
c        write(iwr,*) 'sumup: h after diis'
c        call sumup(q(iscr1x),l2,iwr)
         call dcopy(l2,q(iscr1x),1,h,1)
         call end_time_period(TP_DIIS)
         endif
         diffd=diff

         if(odbg)write(iwr,*)'extrap fock'
         if(odbg)call prtril(h,l1)

         if(omemck)write(iwr,*)'3'
         if(omemck)call ma_summarize_allocated_blocks
c     
c     ------ take precautions if tester is increasing again
c     
         if(iter.ne.1 .and. ondiis .and. 
     1      diffd.ge.diffdp .and. diffdp.ge.diffdpp) then
            nstore=0             
            if(pr_conv)write(iwr,*)'tester up!'
            call dcopy(l2,hsav,1,h,1)
            ondiis=.false.
         endif
c     
c     ----- read in fock tranformation matrix -----
c     transform hamiltonian matrix
c     
c     - q- at x(i30) orthonormalizing transformation
c     - h- at x(i10)
c     -h'- at x(i50) transformed h matrix
c     
c     ----- if vshift is true, use the previous set of
c     molecular orbitals to transform the hamiltonian matrix
c     
c     tdown   q(i30) - result
c     q(jblkqa,s) vsav, qs - input
c     
         if(oshift) then
            call tdown(q(i30),ilifq,vsav,ilifq,l0)
         else
            call tdown(q(i30),ilifq,q(iovq),ilifq,l0)
         endif
c     
c     q(iht) = Qt H qq(ivoff+),  H= h, qq(ivoff+)=q(i30)
c     
         call tran1eo(q, h, q(iht), q(i30), l1) 
c     write(iwr,*) 'sumup: tran1eo-1'
c     call sumup(q(iht),l2,iwr)
         if(omemck)write(iwr,*)'4'
         if(omemck)call ma_summarize_allocated_blocks

         if(odbg)then
            write(iwr,*)'transformed fock'
            call prtril(q(iht),l1)
         endif

         diff=0.0d0
         ii=nocmx+1
         if(ii.le.l0)then
            do 250 i=ii,l0
               ij=iky(i)+iht
               loop=idamax(nocmx,q(ij),1)
               if(loop.gt.0) diff=dmax1(diff,dabs(q(ij+loop-1)))
 250        continue
         endif
         if(ondiis)diff=diffd
         dmplim = dmax1(dmpcut,2.0d0)

c         write(iwr,*)'after check fock',nocmx,diff

         ocvged = (damp.lt.dmplim) 
     &        .and. (diff.lt.acurcy) 
     &        .and. (iter.gt.1)

         if (ocvged) go to 321
c     
c     shift the diagonal of the transformed
c     h matrix by rshift for the virtual part.
c     
         rshift=0.0d0
         if(.not.ondiis)
     *        call shiftq(q(iht),nocmx,0,l0,
     *        de,dep,iterv,1,gapa1,ibrk,gapa2)

c         write(iwr,*)'after shiftq',nocmx,rshift
c
      force_serial = .true.
c     
c     ----- diagonalize new hamiltonian matrix -----
c     
c     -h- at x(i50)
c     -v- at x(i30)
c     
         m2=2
         diaacc = diff*5.0d-3
         if(diaacc.gt.5.0d-5) diaacc = 5.0d-5
         if(diaacc.lt.1.0d-11) diaacc = 1.0d-11
c     
         call start_time_period(TP_DIAG)
c        write(iwr,*) 'sumup: into jacobi'
c        call sumup(q(iht),l2,iwr)
         call jacobi(q(iht),iky,l0,v,ilifq,l1,
     &        e,m2,lock,diaacc)
c        write(iwr,*) 'sumup: back from jacobi'
c        call sumup(v,l3,iwr)
         call end_time_period(TP_DIAG)
c     
c     ----- back-transform the eigenvectors -----
c     
c     -v- at x(i30)
c     -d- at x(i41)
c     -q- at x(i10)
c     scratch area at x(i50)
c     
c     transform vt with q(jblka,s)
c     
         if(omemck)call ma_summarize_allocated_blocks

         iscr = igmem_alloc(l3)
         if(oshift) then
           call tfsqc(v,vsav,q(iscr),l0,l1,l1)
         else
           call tfsqc(v,q(iovq),q(iscr),l0,l1,l1)
         endif
         call gmem_free(iscr)
c     write(iwr,*) 'sumup: tfsqc-1'
c     call sumup(v,l3,iwr)
c
c     
c     ----- if vshift is true, reorthogonalize the vectors and
c     thereby the transformation matrix every igs th iteration.
c     
         ig = mod(iter,igs)
         if (ig .eq. 0) then
         call start_time_period(TP_ORFOG)

            if(pr_conv)write(iwr,*)'re-orthogonalize..'

            iscr1 = igmem_alloc(l3)
            iscr2 = igmem_alloc(l3)
c
c  tranform S
c
            call tran1eo(q, q(istx), q(iscr1), v, l1)
c     write(iwr,*) 'sumup: transform-s'
c     call sumup(q(iscr1),l2,iwr)
c     
c     orthogonalise v  according to S in q(iscr1)
c 
            call orfog(v,v,q(iscr1),q(iscr2),
     &           iky,ilifq,l1,l1,1)

            call gmem_free(iscr2)
            call gmem_free(iscr1)

            if(omemck)call ma_summarize_allocated_blocks

         call end_time_period(TP_ORFOG)
         endif
c     
c     copy orthogonalised vectors -> vsav
c     
         call dcopy(l3,v,1,vsav,1)
         call wrt3(vsav,l3,ibl3qa,idaf)
c     
c     transform v to AO basis
c     
         call tdown(v,ilifq,v,ilifq,l0)
c     
c     unshift virtual energies
c     
         ig=nocmx+1
         if(ig.le.l0)call vsadd(e(nocmx+1),1,-rshift,
     &        e(nocmx+1),1,l0-nocmx)

         if (pr_vecs) then
c     
c     eigenvectors
c     
            write (iwr,9148) 
 9148       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors'/
     +           50x,12('-'))
            call prev(v,e,lprnt,l1,l1)
         endif
c     
c     ----- form density matrix -----
c     
c     -v- at x(i30)
c     -o- at x(i80)
c     
c     load occs -> occ
c
      if(lock_occupations)then

c
c disable the llvmo call below
         assign = .false.
c
c locking of occupations according to fragment occupations
c
         if(first_cycle)then

c
c determine MO occupancies from FMO parentage
c
c transform to fragment orbital basis
            iscr = igmem_alloc(l3)
            call tfsqc(v,q(ifqi),q(iscr),l0,l1,l1)

c            write (iwr,91481) 
c91481       format(//1x,100('-')//
c     +           50x,12('-')/
c     +           50x,'eigenvectors in frag MO basis'/
c     +           50x,12('-'))
c            call prev(v,e,lprnt,l1,l1)
            yavr=yblnk
            do i = 1, num
               qmax = 0.0d0
               do j = 1, num
                  if( abs(v((i-1)*num + j)) .gt. abs(qmax))then
                     qmax = v((i-1)*num + j)
                     iforb = j
                  endif
               enddo
               occ(i) = 2.0d0 * dble( ifocc(iforb) )
               if(odbg)write(iwr,*)'MO ',i,' FMO ',iforb, qmax, occ(i)
            enddo
c
c check counts, abort if number of electrons changed
            j = 0
            do i = 1, num
               if(occ(i) .gt. 1.0d-10)then
c                  write(iwr,*)'occ',i,occ(i)
                  nocmx = i
                  j = j + 1
               endif
            enddo
            if(j .ne. na)then
               write(iwr,*)'counts',j,na
               write(iwr,*)'problem setting occupancies'
               assign = .true.
            endif
c back transform the vectors 
            call tfsqc(v,fq,q(iscr),l0,l1,l1)

            call gmem_free(iscr)
c
c store current orbs and occupancies for next time
c
            call dcopy(l3,v,1,q(ivsav),1)
            do i = 1, num
               ioccsav(i)=0
               if(occ(i) .gt. 1.0d-10)ioccsav(i) = 1
            enddo

            first_cycle = .false.

         else
c
c try and set occs based on orbitals of previous cycle
c

            iscr = igmem_alloc(l3)
            ifpvi = igmem_alloc(l3)
            iscr1 = igmem_alloc(len)
            iscr2 = igmem_alloc(len)

            call dcopy(l1*l1,q(ivsav),1,q(ifpvi),1)
            call minvrt(q(ifpvi),l1,deter,
     &           q(iscr1),q(iscr2))

            if(odbg)write(iwr,*)'deter pv',deter

c     use orbitals of previous cycle
            call tfsqc(v,q(ifpvi),q(iscr),l0,l1,l1)

            call gmem_free(iscr2)
            call gmem_free(iscr1)
            call gmem_free(ifpvi)

            if(odbg)then
            write (iwr,91482) 
91482       format(//1x,100('-')//
     +           50x,12('-')/
     +           50x,'eigenvectors in old MO basis'/
     +           50x,12('-'))
            call prev(v,e,lprnt,l1,l1)
            endif
c
            yavr=yblnk

            do i = 1, num
               qmax = 0.0d0
               do j = 1, num
                  if( abs(v((i-1)*num + j)) .gt. abs(qmax))then
                     qmax = v((i-1)*num + j)
                     iforb = j
                  endif
               enddo
               occ(i) = 2.0d0 * dble( ioccsav(iforb) )
               if(odbg)
     &              write(iwr,*)'MO ',i,' old MO ',iforb, qmax, occ(i)
            enddo

            j = 0
            do i = 1, num
               if(occ(i) .gt. 1.0d-10)then
                  nocmx = i
                  j = j + 1
               endif
            enddo
            if(j .ne. na)then
               write(iwr,*)'problem setting occupancies'
               write(iwr,*)'counts',j,na
               assign = .true.
            endif

c previous cycle
            call tfsqc(v,q(ivsav),q(iscr),l0,l1,l1)
c     write(iwr,*) 'sumup: previous-cycle - v'
c     call sumup(v,l3,iwr)
            call gmem_free(iscr)

c
c now copy the current orbitals and occupations
            call dcopy(l3,v,1,q(ivsav),1)
            do i = 1, num
               ioccsav(i)=0
               if(occ(i) .gt. 1.0d-10)ioccsav(i) = 1
            enddo

         endif
      else
c no locking, use llvmo
         assign = .true.
      endif
c
c make assignment on basis of energy
c

      if(assign)then
         call llvmo(e,occ,na,nocmx,l1)
         call dscal(l1,2.0d0,occ,1)
         yavr=yblnk
         if(nocmx.gt.na)yavr=yav
      endif
         call dmtx(p,v,occ,iky,nocmx,l1,l1)
c     write(iwr,*) 'sumup: p - density'
c     call sumup(p,l2,iwr)
         if (pr_dens) then
            write (iwr,9168)
 9168       format(/20x,14('*')/20x,'density matrix'/20x,14('*'))
            call prtril(p,l1)
         endif

 321      continue

         tim1=cpulft(1)
c     
c     ----- check and print convergence behavior -----
c     
         if(pr_conv)then
c           write (iwr,9188) iter,kcount,etot,
c    1      ehf,de,diff,rshift,damp,derror,delt,tim1,yavr
c9188       format(1x,i3,1x,i3,3f16.8,f15.8,f10.3,f8.3,f15.9,
c    +             2f10.3,a3)
            write(iwr,90282)
90282       format(/1x,103('=')/
     *     3x,'cycle',11x,'total',6x,'electronic',9x,'e conv.',
     *     9x,'tester',3x,'virtual',1x,'damping',11x,'diis'/
     *     18x,'energy',10x,'energy',36x,'shift'/1x,103('='))
            write (iwr,9189) iter,kcount,etot,ehf,de,diff,
     *      rshift,damp,derror
 9189       format(1x,i3,1x,i3,3f16.8,f15.8,f10.3,f8.3,f15.9)
         endif

         dmplim = dmax1(dmpcut,2.0d0)
         ocvged = (damp.lt.dmplim) .and. 
     &        (diff.lt.acurcy) .and. 
     &        (iter.gt.1)

         if (ocvged) go to 400
c     
c     ----- exit in case of time limit -----
c     ----- punch out restart data -----
c     
         call timit(0)
         call timrem(tlefti)
         skale=1.1d0
         if (tlefti .lt. skale*(tlefts-tlefti)/iter) then
           write(iwr,*)'scf timed out'
           istat = 1
           finished = .true.
           irest=3
           nindmx = 0
           write(iwr,9367)
           call texit(0,irest)
         else if (iter .gt. maxcyc) then
c     too many iterations
           write (iwr,9288)
           etot = 0.0d0
           irest=3
           nindmx = 0
           ehf0 = ehf
           ehf = -en
           istat = 2
           finished = .true.
         endif
 400     continue

         if(ocvged)finished = .true.

      enddo
c
 460  continue
      if(ocvged) then
       eel = ehf
       istat = 0
      endif
      if(pr_conv)then
         if (ocvged)then
            write (iwr,9208)
 9208       format(/20x,16('-')/20x,'energy converged'/,20x,16('-'))
         else
            write (iwr,9209)
 9209    format(/20x,22('!')/
     +           20x,'*** no convergence ***'/,20x,22('!'))
         endif
         write (iwr,9368) iter,tim1,ehf,en,etot
         if(yavr.ne.yblnk)write(iwr,9369)
      endif

 9368 format(  /20x,14('-')/
     * 20x,'final energies   after',i4,' cycles at ',f8.2,' seconds'
     *     ,a10,' wall'/
     *     20x,14('-')//
     *     20x,'electronic energy',f18.10/
     *     20x,'nuclear energy   ',f18.10/
     *     20x,'total energy     ',f18.10)
 9369 format(//20x,90('*')/
     *     20x,'*',88x,'*'/20x,'*',88x,'*'/
     *     20x,'*',31x,'warning  state is averaged',31x,'*'/
     =     20x,'*',88x,'*'/20x,90('*')//)
 9367 format(/10x,26('*')/
     +     10x,'*** warning ***'/
     +     10x,'scf has not converged '/
     +     10x,'this job must be restarted'/
     +     10x,26('*')/)
 9288 format(/
     +     20x,30('-')/
     +     20x,'excessive number of iterations'/
     +     20x,30('-'))
 3030 format(/5x,25('*'))
 3000 format(5x,'using full density matrix')
 3020 format(5x,'dlntol =  ',f15.10/5x,25('*'))

c
c  return SA vectors in v as well as vsav
c
      call dcopy(l3, vsav, 1, v, 1)
c
c  release dynamic memory
c
      call gmem_free(iovq)
      call gmem_free(istx)
      call gmem_free(ihx)
      call gmem_free(ih0x)
      call gmem_free(i30)
      call gmem_free(iht)
ccc      call gmem_free(iscr3x)
      call gmem_free(iscr2x)
      call gmem_free(iscr1x)
      call gmem_free(ivsav)
      call gmem_free(ifqi)
c     if(odiis)call gmem_free(iqdiis)
c     call gmem_free(ifsave)
c
c  this helps catch memory corruption 'till we
c  have diis/extrap sorted
c
c      write(iwr,*)'memory check:'
c      
c      call ma_summarize_allocated_blocks

      return

      end

      subroutine dmhstar(q, d, f, rdmat, istat)
c
c
c  oife controls inclusion of inter-fragment exchange integrals
c  oifc controls inclusion of inter-fragment coulomb  integrals
c
c  basis functions are classified 
c
c     ----- ***********************************  ------
c
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
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c ...
c ... timing statistics
c ...
      real*8 begin, ebegin, timsec, walsec, tstart, estart
      common/statis/begin,ebegin,timsec(50),walsec(50),
     + tstart,estart
c
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
      dimension d(*),f(*),q(*),rdmat(*)
      character*10 fnm
      character*7 snm
      data fnm,snm/"morokuma.m","dmhstar"/
c
      inull = igmem_null()
      if(istat.ne.1) call vclr(f,1,nx)
c
      do m=1,num
         nij=ikyp(m)
         d(nij)=d(nij)*0.5d0
      enddo
c
      call start_time_period(TP_DHSTAR)
      call timana(5)
      iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,"prefac",
     +                         IGMEM_DEBUG)
      call rdmake(q(iprefa))
      call jandk(zscftp,q,f,q(inull),q(inull),d,q(inull),q(iprefa),
     +           rdmat)
      call gmem_free_inf(iprefa,fnm,snm,"prefac")
      call end_time_period(TP_DHSTAR)
      call cpuwal(begin,ebegin)
c
      do 3   m = 1,num
         nij = ikyp(m)
         d(nij) = d(nij)+d(nij)
 3    continue
      call dscal(nx,0.5d0,f,1)
c
      return
      end

      subroutine dbuild_morok(fock,dmat,g)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension g(*),fock(*),dmat(*)
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
c     dlnmxd: the maximum value of the logarithm of the reduced
c             density matrix.
c     dlnmxs: the maximum value of the logarithm of the Schwarz 
c             integrals.
c     dlntol: minus the logarithm of the direct SCF integral cutoff.
c
      real*8 dlnmxd, dlnmxs, dlntol, tolitr, deltol, delfac
      integer intcut, intmag
      integer ibl171, lentri, itrtol, iexch, iof171, m171t
      integer nhamd, nsheld
      logical odscf, oimag, odnew, opref, odelta
      common/cslosc/dlnmxd,dlnmxs,dlntol,tolitr(3),deltol,delfac,
     + odscf,
     + intcut(10),oimag,intmag(1060),
     + ibl171,lentri,itrtol(3),iexch,odnew,opref,iof171(6),m171t,
     + nhamd, nsheld,odelta
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
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
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
      integer ik, ijx, ijy, ijz, klx, kly, klz
      integer ijgt, klgt
      common /indez/ ijgt(225),ijx(225),ijy(225),ijz(225),
     +       ik(225),klgt(225),klx(225),kly(225),klz(225)
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      logical otestif, oif
      external otestif

c
c     ----- DSCF fock builder for morokuma energy decomposition
c
      ijn = 0
      jmax = maxj

      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      n1 = ijgt(ijn)
      int2 = locj + j
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+klgt(kln)
      val = g(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        oif = otestif(i1,i2,i3,i4)
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
c
c Coulomb term
c
        if((.not.  oif) .or. oifc)then
           val2=val+val
           val4=val2+val2
           f12 = val4*dmat(itr34) + fock(itr12)
           fock(itr34) = val4*dmat(itr12) + fock(itr34)
           fock(itr12) = f12
        endif
c
c exchange
c
        if((.not.  oif) .or. oife)then
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2

         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
        endif
  200 continue
  220 continue
  240 continue
  260 continue
      end

      subroutine dbuild70_morok(fock,dmat,gout)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension gout(*),fock(*),dmat(*)
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
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      logical out, outv
      common/shlt/tol,cutoff,icount,ic4,out,
     +           isti, jsti, ksti, lsti, lastb, lastu, outv,
     +           len4, lennx
c
      real*8 qq4
      integer lit, ljt, lkt, llt, loci, locj, lock, locl
      integer mini, minj, mink, minl, maxi, maxj, maxk, maxl
      integer nij, ij, kl, ijkl, ncontr
      common /shlnos/ qq4,lit,ljt,lkt,llt,loci,locj,lock,locl,
     + mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     + nij,ij,kl,ijkl,ncontr
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
      common/flip70/ib1,jb1,kb1,lb1,ib2,jb2,kb2,lb2,ib3,jb3,kb3,lb3
c
      integer ishell, jshell, kshell, lshell
      integer inew, jnew, knew, lnew
      common/shlg70/ishell,jshell,kshell,lshell,inew,jnew,knew,lnew
c
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
c  control parameters for morokuma energy decomposition
c  analysis
c
      logical omorok, oifc, oife, oifh
      integer imode, ifrag, ifocc, imoro
      real*8 emorok
      common/morok1/emorok(0:4,2),imode,ifrag(maxorb),
     &     ifocc(maxorb),imoro, omorok, 
     &     oifc, oife, oifh
      character*20 fragname
      common/morokc/fragname(3)

c
c  settings
c
      integer ALL_BLOCKS
      parameter (ALL_BLOCKS=0)
      
      integer ESX
      parameter (ESX=1)

      integer ESX_PLX
      parameter (ESX_PLX=2)

      integer ESX_CTT
      parameter (ESX_CTT=3)

      integer ESX_EX
      parameter (ESX_EX=4)

      integer YES_K
      parameter (YES_K=1)

      integer NO_K
      parameter (NO_K=2)

c
      dimension ib(4,4)
      data ib/0,0,0,0,64,16,4,1,128,32,8,2,192,48,12,3/
c
      logical otestif, oif
      external otestif

c
c     ----- DSCF fock builder for morokuma energy decomposition
c
      ijn = 0
      oident = ishell .eq. kshell .and. jshell .eq. lshell
      mini = kmin(ishell)
      minj = kmin(jshell)
      mink = kmin(kshell)
      minl = kmin(lshell)
      maxi = kmax(ishell)
      maxj = kmax(jshell)
      maxk = kmax(kshell)
      maxl = kmax(lshell)
      loci = kloc(ishell)-mini
      locj = kloc(jshell)-minj
      lock = kloc(kshell)-mink
      locl = kloc(lshell)-minl
      oianj = ishell .eq. jshell
      okanl = kshell .eq. lshell
      jmax = maxj

      do 260 i = mini,maxi
      if (oianj) jmax = i
      int1 = loci + i
      do 240 j = minj,jmax
      ijn = ijn+1
      int2 = locj + j
      n1 = ib(ib1,i)+ib(jb1,j)+1
      ii1 = max(int1,int2)
      ii2 = min(int1,int2)
      lmax = maxl
      kln = 0
      do 220 k = mink,maxk
      if (okanl) lmax = k
      int3 = lock + k
      do 200 l = minl,lmax
      kln = kln+1
      if (oident .and. kln .gt. ijn) go to 240
      nn = n1+ib(kb1,k)+ib(lb1,l)
      val = gout(nn)
      if ( dabs(val) .lt. cutoff) go to 200
        int4 = locl + l
        ii3 = max(int3,int4)
        ii4 = min(int3,int4)
        if(ii1.ge.ii3) then
           i1 = ii1
           i2 = ii2
           i3 = ii3
           i4 = ii4
        else
           i1 = ii3
           i2 = ii4
           i3 = ii1
           i4 = ii2
        endif
        oif = otestif(i1,i2,i3,i4)
        itr12=iky(i1)+i2
        itr13=iky(i1)+i3
        itr14=iky(i1)+i4
        itr34=iky(i3)+i4
        itr23=iky(max(i2,i3))+min(i2,i3)
        itr24=iky(max(i2,i4))+min(i2,i4)
c
c Coulomb term
c
        if((.not.  oif) .or. oifc)then
           val2=val+val
           val4=val2+val2
           f12 = val4*dmat(itr34) + fock(itr12)
           fock(itr34) = val4*dmat(itr12) + fock(itr34)
           fock(itr12) = f12
        endif
c
c exchange
c
        if((.not.  oif) .or. oife)then
         val2=val+val
         val13=val
         val14=val
         if(i1.eq.i3 .or. i2.eq.i4) val13=val2
         if(i2.eq.i3) val14=val2

         f23 = fock(itr23) - val14*dmat(itr14)
         f14 = fock(itr14) - val14*dmat(itr23)
         f13 = fock(itr13) - val13*dmat(itr24)
         fock(itr24) = fock(itr24) - val13*dmat(itr13)
         fock(itr23) = f23
         fock(itr14) = f14
         fock(itr13) = f13
        endif
  200 continue
  220 continue
  240 continue
  260 continue
      end
      subroutine ver_morokuma(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/morokuma.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
