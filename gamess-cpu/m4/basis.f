c 
c  $Author: jmht $
c  $Date: 2010-05-31 13:43:04 +0200 (Mon, 31 May 2010) $
c  $Locker:  $
c  $Revision: 6135 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/basis.m,v $
c  $State: Exp $
c  
c     deck=basis
      subroutine atoms(q,csinp,cpinp,cdinp,cfinp,cginp,
     + ztita,czana,zatagg,
     + csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     + ict,newsh,iso,isoc,natmax,nshmax,ngsmax,obasin,onel,
     + oshell)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character *10 zbasis, ztype
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
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
      common/junk/trmat(432),ptr(3,144),dtr(6,288),
     + ftr(10,480),gtr(15,720),
     + intyp(mxshel),ns(maxat),ks(maxat),
     + intypi(mxshel),nsi(maxat),ksi(maxat),
     + exi(mxprim),csi(mxprim),cpi(mxprim),cdi(mxprim),
     +             cfi(mxprim),cgi(mxprim),
     + kstari(mxshel),katomi(mxshel),ktypi(mxshel),kngi(mxshel),
     + kloci(mxshel),kmini(mxshel),kmaxi(mxshel),nshi
c
      real*8 stos1, stos2, stos3, stos4, stos5, stos6, stos7
      real*8 rleesf
      real*8 c11al, c11be, c31al, c31be, c33al, c33be
      integer nbfs, minf, maxf, nangm
c
      common /datbas/ stos1(54),stos2(54),stos3(54),stos4(54),
     + stos5(54),stos6(54),stos7(54),
     + rleesf(36,2),
     + nbfs(24),minf(24),maxf(24),nangm(24),
     + c11al(12),c11be(12),c31al(12),c31be(12),c33al(12),c33be(12)
c
c
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
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
      common/junkc/ ylabel(26)
      dimension ict(natmax,*),newsh(nshmax,*),iso(*),isoc(*)
      dimension ztita(*),czana(*),zatagg(*)
      dimension exx(mxprms),css(mxprms),cpp(mxprms),cdd(mxprms)
      dimension scc(8)
      dimension csinpi(ngsmax),cpinpi(ngsmax),cdinpi(ngsmax)
      dimension cfinpi(ngsmax),cginpi(ngsmax)
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax)
      dimension cfinp(ngsmax),cginp(ngsmax)
c
      parameter(numbas = 32)
      dimension ztagbs(numbas),zshell(5),zecp(8)
c
      dimension q(*)
c
      data ztagbs/
     +   'sv',  'tzv',  'tzvp',  'mini1','midi1',
     +'mini4','midi4','ecptzv','ecptzvp',  'svp',
     +   'dz',  'dzp','ecpmin',  'ecpdz',  'sto',
     + 'exte',  'svr',  'tzvr',  'tzvpr', 'svpr',
     +  'dzr', 'dzpr',
     +  'cc-pvdz','cc-pvtz','cc-pvqz','cc-pv5z','dft',
     +  'ecp', 'a1-dgaus', 'a2-dgaus', 'demon', 'f-ahlric'  /
      data zshell /'sv','dz','tz','qz','5z' /
      data zecp / 'sbkjc','cep','lanl','lanl2','crenbl',
     +            'crenbs','strlc','strsc' /
c
      data ybas1/'sto'/
      data pi32 /5.56832799683170d0/
      data toll/1.0d-10/
      data pt187,pt6562 /1.875d+00,6.5625d+00/
      data pt5,pt75 /0.5d0,0.75d0/
      data zblank /' '       /
      data yend,yfini,ystop/'end','fini','stop'/
      data z321/'3-21g'/
      data ntype/25/
c
c     read in unique centers and atomic basis sets grouped
c     in shells
c
      write(iwr,9068)
      igauss = 0
      ng5 = ngsmax * 5
      call vclr(exi,1,ngsmax)
      call vclr(csi,1,ng5)
      call vclr(csinpi,1,ng5)
c
      call setsto(nonsym,0,nsi)
      nshi   = 0
      loci   = 0
      ngai   = 0
      ocep=.false.
      olanl=.false.
      opol=.false.
      opolx=.false.
      opolh=.false.
      orydx=.false.
      orydh=.false.
      oryd =.false.
      odunn=.true.
      oecpmin=.true.
      oecpdz =.false.
      oignore = .false.
      zbasis=zblank
      osplit=.false.
      oshell=.false.
      ztype =zblank
      sc=0.0d0
      oinvrt = .false.
      onwchem = .false.
      oend = .false.
      do 3000 i=1,8
 3000 scc(i)=0.0d0
c
      if(obasin) go to 3050
c ----- default basis is 3-21g
      zbasis=z321
      ztype=ztagbs(1)
      ityp=25
      go to 3051
 3050 if(jump.eq.1)go to 160
c
c     check for ignore sub-directive. I.e. whether it is okay to have
c     no basis set on (some of) the atoms
c
      call inpan(ztype)
      if(ztype.eq.'ignore') then
        oignore = .true.
      else
        jrec = jrec-1
      endif
c
c     check for requested inversion of exps and coefficients
c     or non-standard basis specs e.g. nwchem or g94
c
      call inpan(ztype)
      if(ztype.eq.'reverse'.or.ztype.eq.'invert') then
        oinvrt = .true.
        go to 160
      else if (ztype.eq.'nwchem'.or.ztype.eq.'g94') then 
        onwchem = .true.
        go to 160
      else if (ztype.eq.' ') then
        go to 160
      else
        jrec = jrec -1
      endif

c
c ----- shorthand basis specification
c       e.g. basis dunning, basis sto3g etc
c       or alternatives e.g. tz/tzv, dzv/dz etc
c       decide type
c       default scaling only
c
      call inpan(ztype)
c
c ... sv/tzv/tzvp/
c
      ytype=ytrunc(ztype)
      jumpj=locatc(ztagbs,numbas,ztype)
      if(jumpj.eq.0)  then
c check for sto-ng
         if(ztype(1:3).eq.'sto') then
            call atosto(ztype,ybasis,igauss)
            mfgbas = 11
            jumpj=15
c check for extended
         else if(ztype(1:4).eq.'exte') then
            mfgbas = 12
            jumpj = 16
c check for tz rather than tzv
         else if(ztype.eq.'tz') then
            jumpj = 2
c check for tzp rather than tzvp
         else if(ztype.eq.'tzp') then
            jumpj = 3
c check for dzv rather than dz
         else if(ztype.eq.'dzv') then
            jumpj = 11
c check for dzvp rather than dzp
         else if(ztype.eq.'dzvp') then
            jumpj = 12
c check for dunning
         else if(ztype(1:4).eq.'dunn')then
            mfgbas = 3
            jumpj=1
            odunn = .true.
c check for ahlrichs
         else if(ztype(1:4).eq.'ahlr')then
            mfgbas = 3
            jumpj=1
            odunn = .false.
c check for 6-311g
         else if(ztype(1:6).eq.'6-311g')then
            mfgbas = 6
            if(ztype(7:7).eq.'*') opolx=.true.
            if(ztype(8:8).eq.'*') opolh=.true.
            jumpj=2
            if(opolx.and.opolh) jumpj=3
         else if(ztype(5:5).eq.'g'.or.ztype(6:6).eq.'g'.or. 
     +           ztype(7:7).eq.'g') then
c check for n-m1g
            call ato321(ztype,ibasis,igauss,opolh,opolx,
     +      orydh,orydx,odunn,oerr1,oerr2,oerr3)
            if(oerr1) go to 1530
            if(oerr2) go to 1540
            if(oerr3) then
               zbasis = ztype
               go to 1550
            endif
            mfgbas = ibasis+1
            jumpj=1
         else
            call caserr2('unrecognised basis set option specified')
         endif
          zbasis=ztype
          ztype=ztagbs(jumpj)
          jrec=jrec-1
      endif
      ityp=jumpj+ntype
      go to (10, 50, 50, 70, 80, 90,100,110,120,
     +       10, 30, 30,130,140,141,142, 10, 50,
     +       50, 10, 30, 30, 35, 35, 35, 35, 36, 
     +       37, 35, 35, 35, 35 ), jumpj
c ... sv/svp/svpr/svr
 10    call bastyp(ztype,opol,oryd)
       call inpan(zbasis)
3051   call ato321(zbasis,ibasis,igauss,opolh,opolx,
     +             orydh,orydx,odunn,oerr1,oerr2,oerr3)
      if(oerr1) go to 1530
      if(oerr2) go to 1540
      if(oerr3) go to 1550
      if(ibasis.eq.2) then
        if(odunn) then
          zbasis='dunning'
        else
          zbasis='ahlrichs'
        endif
      endif
       mfgbas=ibasis+1
       go to 15000
c ... dz/dzp/dzpr/dzr
 30    call bastyp(ztype,opol,oryd)
       mfgbas=4
       call inpan(zbasis)
       if(zbasis.eq.'ahlrichs') then
        odunn = .false.
       else
        odunn = .true.
       endif
       go to 15000
c ... ccpv-dz/tz/qz/5z or fitting basis
 35    mfgbas = jumpj - 10
       go to 15000
c ... dft basis (dzvp, dzvp2 or tzvp)
 36    mfgbas = jumpj - 10
       call inpa(ztemp)
       zbasis = ztemp
       go to 15000
c ... ecp basis (sbkjc/cep,lanl/lanl2,crenbl,crenbs,strlc,strsc)
 37    mfgbas = jumpj - 10
       call inpa(ztemp)
       zbasis = ztemp
       go to 15000
c ... tzv/tzvp/tzvr/tzvpr
 50    call bastyp(ztype,opol,oryd)
       mfgbas=5
       call inpan(zbasis)
       if(zbasis(1:6).eq.'6-311g') then
         if(zbasis(7:7).eq.'*') opolx=.true.
         if(zbasis(8:8).eq.'*') opolh=.true.
         mfgbas=mfgbas+1
       else if(zbasis.eq.'ahlrichs') then
         odunn = .false.
       else
         odunn = .true.
       endif
       go to 15000
c ... midi1
 80    osplit=.true.
c ... mini1
 70    mfgbas=7
       go to 15000
c ... midi4
 100   osplit=.true.
c ... mini4
 90    mfgbas=8
       go to 15000
c ... ecptzvp
 120   opol=.true.
c ... ecptzv
 110   mfgbas=10
       go to 15000
c ... ecpdz
 140   oecpmin=.false.
       oecpdz =.true.
       go to 1401
c ... ecpmin
 130   oecpmin=.true.
       oecpdz =.false.
 1401  call inpan(zbasis)
       if(zbasis.eq.'lanl2') then
          olanl=.true.
       else if(zbasis.eq.'cep') then
          ocep=.true.
       else
          zbasis=' '
       endif
       mfgbas=9
       go to 15000
c ... sto
 141   mfgbas=11
       go to 15000
c ... extended
 142   zbasis='ahlrichs'
c
c -----
c
15000  continue

      write(iwr,9569)ztype,zbasis
9569  format(1x,40('*')/
     * 1x,'* basis selected is ',a7,1x,a10,' *'/
     *1x,40('*')/)
c
      if(mfgbas.eq.1.or.mfgbas.eq.2) oextg = .true.
      if(mfgbas.eq.11) oming = .true.
c
c ... drive through nonsym centres
c ... eliminating those without basis functions
c
      if(opol) then
        opolh=.true.
        opolx=.true.
      endif
      if(oryd) then
        orydh=.true.
        orydx=.true.
      endif
      do 3005 nati=1,nonsym
      zatagi=ztag(nati)
      j=nati-1
      if(j.eq.0)go to 3040
      natj=locatc(ztag,j,zatagi)
      if(natj)3005,3040,3005
3040  ierr1 = 0
      ierr2 = 0
      igau  = igauss
      nucz = jsubst(zatagi)
      if(nucz.le.0)go to 3005
       if(nucz.eq.1) then
       opol=opolh
       oryd=orydh
       else
       opol=opolx
       oryd=orydx
       endif
c
      go to (13000, 13010, 13020,   13030,   13040, 13050,
c            ------SV/SVP-------   -DZ/DZP-  --TZV/TZVP--
     *       13060, 13070, 13080, 13090, 13100, 13200,
c            -MINI- -MIDI- -----ECP----  -STO-  -EXTE-
     +       13300, 13400, 13500, 13600, 13640, 13645, 
c          CC-PVDZ  -PVTZ  -PVQZ  -PV5Z   -DFT- -ECP-
     +       13650, 13660, 13670, 13680, 13700), mfgbas
c         a1-dgaus a2-dgaus demon ahlrichs
c
c ---EXTENDED-
c
13200  call ahext(csinpi,cpinpi,cdinpi,sc,scc,nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,ngsmax,
     + ierr1,ierr2,nati)
       go to 3010
c
c ---CC-PVDZ/TZ/QZ/5Z
c
13300  ztemp = zshell(2)
       go to 13700
13400  ztemp = zshell(3)
       go to 13700
13500  ztemp = zshell(4)
       go to 13700
13600  ztemp = zshell(5)
13700  call duncc(ztemp,csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     +  nucz,intypi,nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +  ierr1,ierr2,nati)
       go to 3010
c
c ---DFT---
c
13640 call dftorbs(ztemp,csinpi,cpinpi,cdinpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c ---ECP---
c
13645 call ecpbas(ztemp,csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c    test a1-dgauss fitting basis sets 
c
13650  call dgauss_a1_fit(csinpi,cpinpi,cdinpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c    test a2-dgauss fitting basis sets 
c
13660  call dgauss_a2_fit(csinpi,cpinpi,cdinpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c    test demon fitting basis sets 
c
13670  call demon_fit(csinpi,cpinpi,cdinpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c    test ahrichs fitting basis sets 
c
13680  call ahlrichs_fit(csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 3010
c
c ---SV/SVP---
c
13000  call n21g(csinpi,cpinpi,cdinpi,opol,oryd,
     +  sc,scc,igau,nucz,intypi,
     +  nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +  nshmax,ngsmax,ierr1,ierr2,nati,iwr)
       oshell = .true.
       go to 3010
13010  call n31g(csinpi,cpinpi,cdinpi,cfinpi,opol,oryd,
     +  sc,scc,rleesf,igau,nucz,intypi,nangm,
     +  nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,ierr1,ierr2,
     +  nati,iwr)
       oshell = .true.
       go to 3010
13020  if (odunn) then
        call dunbas(zshell(1),csinpi,cpinpi,cdinpi,opol,oryd,
     +  sc,scc,nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       else
        call ahbas(zshell(1),csinpi,cpinpi,cdinpi,opol,oryd,
     +  sc,scc,nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ngsmax,ierr1,ierr2,nati)
       endif
       go to 3010
c
c ---DZ/DZP---
c
13030 if(odunn) then
       call dunbas(zshell(2),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
      else
       call ahbas(zshell(2),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ngsmax,ierr1,ierr2,nati)
      endif
      go to 3010
c
c ---TZV/TZVP
c
13040 if(odunn) then
       call dunbas(zshell(3),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
      else
       call ahbas(zshell(3),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,ngsmax,
     + ierr1,ierr2,nati)
      endif
      go to 3010
13050 call poptzv(csinpi,cpinpi,cdinpi,opol,oryd,
     + nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     + ierr1,ierr2,nati,iwr)
      go to 3010
c
c ---MINI---
c
13060 call mini1(csinpi,cpinpi,cdinpi,cfinpi,osplit,sc,scc,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ngsmax,nati)
      go to 3030
c
c ---MIDI---
c
13070 call mini4(csinpi,cpinpi,cdinpi,osplit,sc,nucz,intypi,nangm,nbfs,
     +minf,maxf,loci,ngai,nsi,nshmax,ngsmax,ierr1,ierr2,nati,iwr)
      go to 3010
c
c ---ECP---
c
13080 if (ocep) then
       call ecpsbkjc(csinpi,cpinpi,cdinpi,cfinpi,oecpmin,oecpdz,
     + sc,scc,nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     + ierr1,ierr2,nati,stos1,iwr)
      else
       call ecphw(olanl,csinpi,cpinpi,cdinpi,cfinpi,oecpmin,sc,scc,
     + nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     + ierr1,ierr2,nati,stos1,iwr)
      endif
       oshell = .true.
      go to 3010
13090  call potzv(csinpi,cpinpi,cdinpi,opol,sc,scc,nucz,intypi,
     * nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     * ierr1,ierr2,nati,iwr)
       oshell = .true.
      go to 3010
c
c ---STO---
c
13100  continue
      call stong(csinpi,cpinpi,cdinpi,igauss,
     * sc,scc,nucz,intypi,
     * nangm,nbfs,minf,maxf,
     * stos1,stos2,stos3,stos4,stos5,stos6,stos7,
     *loci,ngai,nsi,nshmax,ngsmax,
     * ierr1,ierr2,nati,iwr)
       oshell = .true.
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
3010  if (ierr1 .ne. 0) go to 1500
      if (ierr2 .ne. 0) go to 1520
3030  sc=0.0d0
      do 33000 i=1,8
33000 scc(i)=0.0d0
 3005  continue



       go to 780
c
c     ----- read in shells -----
c
 160  call input
      oshell=.false.
      oecpmin=.true.
      oecpdz =.false.
      ocep=.false.
      olanl=.false.
      opol=.false.
      oryd = .false.
      opolh=.false.
      opolx=.false.
      orydh=.false.
      orydx=.false.
      osplit=.false.
      odunn = .true.
      ybasis = '    '
      if (onwchem) then
       call inpa (ztagc)
       ytmp=ytrunc(ztagc)
       if(ytmp.eq.yend.or.ytmp.eq.yfini.or.ytmp.eq.ystop)
     +  go to 780
       call inpan(ztype)
       ytype=ytrunc(ztype)
      else
       call inpan(ztype)
       ytype=ytrunc(ztype)
       if(ytype.eq.yend.or.ytype.eq.yfini.or.ytype.eq.ystop)
     +  go to 780
      endif
8013  continue
      if(.not.onwchem) then
       call inpa(ztagc)
      endif
c
c ----- check on valid type
c
      ityp=locatc(ylabel,ntype,ytype)
      if (ityp.eq.25) then
         ityp = 22
         ytype = 'l'
         ztype = 'l'
      endif
      if(ityp.le.0) then
        j=locatc(ztagbs,numbas,ztype)
        if(j.eq.0) then
         if(ztype(1:3).eq.'sto') then
          call atosto(ztype,ybasis,igauss)
          j=15
         else if(ztype.eq.'tz') then
          j=2
         else if(ztype.eq.'tzp') then
          j=3
         else if(ztype.eq.'dzv') then
          j=11
         else if(ztype.eq.'dzvp') then
          j=12
         else if(ztype(1:4).eq.'exte') then
          j=16
         else
          call caserr2(
     *    'unknown orbital type specified in basis directive')
         endif
        endif
        ityp=j+ntype
      endif
c
c ----- check on valid centre
c
      nati=locatc(ztag,nonsym,ztagc)
      if(nati.le.0) then
        call caserr2(
     *'attempting to site basis functions on unrecognised centre')
      endif
c
c ----- now extract igauss value for stong and sv basis sets
c ----- input scaling factors
c
      if(ityp.le.15.or.ityp.eq.24)then
c
c ----- sto ng
        call inpan(zbasis)
        call atosto(zbasis,ybasis,igauss)
        go to 7006
      endif
c
      if(ityp.le.ntype) go to 7005
      jumper=ityp-ntype
c             sv   tzv   tzvp  mini1 midi1
      go to (7004, 5100, 5100, 5300, 5200,
c            mini4 midi4 ecptzv ecptzvp
     +       5400, 5500, 5600, 5700,
c             svp    dz    dzp  ecpmin ecpdz   sto  extended
     +       7004,16000,16000, 17000, 17500, 18000, 18100,
c             svr  tzvr  tzvpr  svpr   dzr    dzpr
     +       7004, 5100, 5100,  7004, 16000, 16000, 
c         cc-pvdz  cc-pvtz   cc-pvqz  cc-pv5z  dft    ecp
     +      16100,   16200,  16300,   16400, 16350, 16340 ),
     +                          jumper
c
c
c ----- SV / SVP
c ----- n-31g/n-21g/dunnings sv/ahlrichs sv
c
c     ----- n-31g basis set -----
c           j.chem.phys. 54,  724 (1971)
c                        56, 2257 (1972)
c     4-31g available for hydrogen, boron to fluorine,
c           and phosphorus to chlorine.
c     5-31g and 6-31g available for hydrogen, and
c          and boron to fluorine.
c
 7004 call bastyp(ztype,opol,oryd)
      call inpan(zbasis)
      call ato321(zbasis,ibasis,igauss,opolh,opolx,orydh,orydx,
     +            odunn,oerr1,oerr2,oerr3)
      if (oerr1) go to 1530
      if (oerr2) go to 1540
      if (oerr3) go to 1550
      ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
      if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
      endif
c
c  Allow for the following variations for specifying diffuse and
c  polarisation functions e.g.
c  pol:     (i)   svp atom 3-21g   (including h's)
c           (ii)  sv  atom 3-21g*  (excluding h's)
c           (iii) sv  atom 3-21g** (including h's)
c  diffuse: (i)   svr atom 3-21g   (including h's)
c           (ii)  sv  atom 3-21+g  (excluding h's)
c           (iii) sv  atom 3-21++g (including h's)
c  The above I think is what users are most likely to expect. At least
c  it is consistent with the documentation "PART 3. DATA INPUT - 
c  Pre-directives and CLASS 1 Directives", subsection "Polarisation
c  Basis Sets"
c
      if (.not.opol) then
         if (nucz.eq.1) then
            opol = opolh
         else
            opol = opolx
         endif
      endif
      if (.not.oryd) then
         if (nucz.eq.1) then
            oryd = orydh
         else
            oryd = orydx
         endif
      endif
c
      call inpf(sc)
      do 7007 i=1,8
 7007 call inpf(scc(i))
c
      if (ibasis-1) 7008,7009,7003
7009  call n31g(csinpi,cpinpi,cdinpi,cfinpi,opol,oryd,
     +          sc,scc,rleesf,igauss,nucz,intypi,nangm,
     +          nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     +          ierr1,ierr2,nati,iwr)
      go to 7012
c
 7008 call n21g(csinpi,cpinpi,cdinpi,opol,oryd,
     +          sc,scc,igauss,nucz,intypi,
     +          nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +          nshmax,ngsmax,ierr1,ierr2,nati,iwr)
7012  if (ierr1 .ne. 0) go to 1500
      if (ierr2 .ne. 0) go to 1520
      go to 160
c
c SV and SVP basis sets (Dunning or Ahlrichs)
c
 7003 if(odunn) then
c
c ----- dunnings sv basis
c ----- modern theoretical chem, vol 3, p1 (plenum 1977)
c
c    (2s) hydrogen
c    (3s2p) for li to ne
c    (6s4p) for al,si,p,s,cl
c
       call dunbas(zshell(1),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
      else
c
c ----- ahlrichs sv basis
c       Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c
       call ahbas(zshell(1),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ngsmax,ierr1,ierr2,nati)
      endif
c
      go to 7012
c
c DZ and DZP basis sets (Dunning or Ahlrichs)
c
16000 call bastyp(ztype,opol,oryd)
      ierr1 = 0
      ierr2 = 0
      call inpa(ztemp)
      if(ztemp(1:4).eq.'dunn') then
       odunn = .true.
      else if(ztemp(1:4).eq.'ahlr') then
       odunn = .false.
      else
       odunn = .true.
       jrec=jrec-1
      endif
      atnum = czann(nati)
      nucz =  nint(atnum)
        if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
        endif
c
      call inpf(sc)
c
      do 17007 i=1,6
17007 call inpf(scc(i))
c
      if(odunn) then
c
c ----- dunnings dz and dzp  basis
c     (2s) hydrogen
c     (4s2p) for li to ne
c     (6s4p) for al,si,p,s,cl
c
      call dunbas(zshell(2),csinpi,cpinpi,cdinpi,opol,oryd,
     + sc,scc,nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
c
      else
c
c ----- ahlrichs sv basis
c       Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c
       call ahbas(zshell(2),csinpi,cpinpi,cdinpi,opol,oryd,
     * sc,scc,nucz,intypi,nangm,nbfs,
     * minf,maxf,loci,ngai,nsi,ngsmax,ierr1,ierr2,nati)
c
      endif
c
      go to 7012
c
c ----- dunnings cc basis sets + wachters +f for tm 1st row
c ----- CC-PVDZ/TZ/QZ/5Z
c
16100 ztemp = zshell(2)
       go to 16500
16200  ztemp = zshell(3)
       go to 16500
16300  ztemp = zshell(4)
       go to 16500
16400  ztemp = zshell(5)
16500 ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
      if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
      endif
      call duncc(ztemp,csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     +  nucz,intypi,nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +  ierr1,ierr2,nati)
      go to 7012
c
c ---ECP---
c
16340 call inpa(ztemp)
      ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
      if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
      endif
      call ecpbas(ztemp,csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     + nucz,intypi,nangm,nbfs,
     + minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
       go to 7012
c
c ---- DGauss and DeMon polarized DFT orbitals
c
16350 call inpa(ztemp)
      ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
      if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
      endif
      call dftorbs(ztemp,csinpi,cpinpi,cdinpi,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ierr1,ierr2,nati)
      go to 7012
c
c ----- sto-ng basis sets
c
18000 ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
        if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
        endif
c
      call inpf(sc)
c
      do 17003 i=1,8
17003 call inpf(scc(i))
c
       call stong(csinpi,cpinpi,cdinpi,igauss,
     * sc,scc,nucz,intypi,
     * nangm,nbfs,minf,maxf,
     * stos1,stos2,stos3,stos4,stos5,stos6,stos7,
     *loci,ngai,nsi,nshmax,ngsmax,
     * ierr1,ierr2,nati,iwr)
c
      go to 7012
c
c ----- ecp basis sets / minimal or DZ
c
17500 oecpmin=.false.
      oecpdz =.true.
      go to 17501
17000 oecpmin = .true.
      oecpdz = .false.
17501 ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
      call inpa(ztemp)
        if(ztemp.eq.'cep') then
        ocep=.true.
        else if (ztemp.eq.'lanl2') then
        olanl=.true.
        else
        jrec=jrec-1
        endif
        if(nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
        endif
c
      call inpf(sc)
c
      do 17001 i=1,6
17001 call inpf(scc(i))
      if (ocep) then
       call ecpsbkjc(csinpi,cpinpi,cdinpi,cfinpi,oecpmin,oecpdz,
     + sc,scc,nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     + ierr1,ierr2,nati,stos1,iwr)
      else
       call ecphw(olanl,csinpi,cpinpi,cdinpi,cfinpi,oecpmin,sc,scc,
     + nucz,intypi,
     + nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     + ierr1,ierr2,nati,stos1,iwr)
      endif
c
      go to 7012
c
c ----  default
c ----- TZV or TZVP basis sets, 6-311G, Dunning or Ahlrichs  -----
c
5100   call bastyp(ztype,opol,oryd)
       ierr1 = 0
       ierr2 = 0
       call inpa(ztemp)
         if(ztemp(1:6).eq.'6-311g') then
          opople = .true.
         else if(ztemp(1:4).eq.'dunn') then
          odunn = .true.
          opople = .false.
         else if(ztemp(1:4).eq.'ahlr') then
          odunn = .false.
          opople = .false.
         else
          opople=.false.
          odunn = .true.
          jrec=jrec-1
         endif
       atnum = czann(nati)
       nucz =   nint(atnum)
         if(nucz.le.0)then
         call inpa(ztagbq)
         nucz=jsubst(ztagbq)
         endif
c
       call inpf(sc)
       do 10002 i=1,6
10002  call inpf(scc(i))
c
       if(opople) then
         call poptzv(csinpi,cpinpi,cdinpi,opol,oryd,
     +   nucz,intypi,
     +   nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     +   ierr1,ierr2,nati,iwr)
       else if(odunn) then
c
c ----- tzv -----
c ----- dunnings (10s6p/5s3p) for li to ne
c ----- mcleans (12s9p/6s5p) for na to ar
c ----- dunnings (5s/3s) for h
c ----- wachters (14s9p5d/10s8p3d)
c
c ----- tzvp ---
c       polarisation functions from ahlrichs/taylor (1st+2nd row)
c       j.chim.phys. 1981, 78.
c
         call dunbas(zshell(3),csinpi,cpinpi,cdinpi,opol,oryd,
     +   sc,scc,nucz,intypi,nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +   ierr1,ierr2,nati)
       else
c
c     Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c     tzv basis sets 
c
         call ahbas(zshell(3),csinpi,cpinpi,cdinpi,opol,oryd,
     +   sc,scc,nucz,intypi,nangm,nbfs,minf,maxf,loci,ngai,nsi,
     +   ngsmax,ierr1,ierr2,nati)
       endif
c
       go to 7012
c
c +++ EXTENDED
c     Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c
18100 ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =   nint(atnum)
      if(nucz.le.0)then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
      endif
      call inpf(sc)
      do 10004 i=1,6
10004 call inpf(scc(i))
      call ahext(csinpi,cpinpi,cdinpi,sc,scc,nucz,intypi,
     +nangm,nbfs,minf,maxf,loci,ngai,nsi,ngsmax,
     +ierr1,ierr2,nati)
      go to 7012
c
c +++ mini-1 and midi-1
c
5200  osplit=.true.
5300  atnum = czann(nati)
      nucz =  nint(atnum)
        if (nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
        endif
      call inpf(sc)
      do 13004 i=1,2
13004 call inpf(scc(i))
      call mini1(csinpi,cpinpi,cdinpi,cfinpi,osplit,sc,scc,
     +  nucz,intypi,nangm,nbfs,
     +  minf,maxf,loci,ngai,nsi,ngsmax,nati)
      go to 160
c
c +++ mini-4 and midi-4
c
5400  osplit=.true.
5500  ierr1 = 0
      ierr2 = 0
      atnum = czann(nati)
      nucz =  nint(atnum)
        if (nucz.le.0) then
        call inpa(ztagbq)
        nucz=jsubst(ztagbq)
        endif
      call inpf(sc)
      call mini4(csinpi,cpinpi,cdinpi,osplit,sc,nucz,intypi,nangm,nbfs,
     +minf,maxf,loci,ngai,nsi,nshmax,ngsmax,ierr1,ierr2,nati,iwr)
c
       go to 7012
c
c ----
c ----- ecptzv and ecptzvp -----
c
5700   opol=.true.
5600   ierr1 = 0
       ierr2 = 0
       atnum = czann(nati)
       nucz =   nint(atnum)
         if(nucz.le.0)then
         call inpa(ztagbq)
         nucz=jsubst(ztagbq)
         endif
       call inpf(sc)
       do 10003 i=1,2
10003  call inpf(scc(i))
c
       call potzv(csinpi,cpinpi,cdinpi,opol,sc,scc,nucz,intypi,
     * nangm,nbfs,minf,maxf,loci,ngai,nsi,nshmax,ngsmax,
     * ierr1,ierr2,nati,iwr)
c
       go to 7012
c
7005  ibasis=0
7006  nshi = nshi+1
      if (nshi .gt. nshmax) go to 1500
      kmini(nshi) = minf(ityp)
      kmaxi(nshi) = maxf(ityp)
      kstari(nshi) = ngai+1
      katomi(nshi) = nati
      ktypi(nshi) = nangm(ityp)
      intypi(nshi) = ityp
      kloci(nshi) = loci+1
      loci = loci+nbfs(ityp)
      nsi(nati) = nsi(nati)+1
      if (ybasis .eq. ybas1) go to 260
c
c     ----- general basis set -----
c
      k=kstari(nshi)
      igauss=0
7010  call input
      ybasis = '    '
      call inpa(ztest)
      ytest=ytrunc(ztest)
      if(ytest.ne.yend. and . ytest.ne.yfini. and .ytest.ne.ystop)
     * go to 8000
      oend=.true.
      go to 8010
 8000 if(onwchem) then
       ztagc = ztest
       call inpa(ztest)
       ytest=ytrunc(ztest)
      endif

      i=locatc(ylabel,ntype,ytest)
      if(i)8010,9650,8010
c
c ----- input primitive
c
 9650 i=locatc(ztagbs,numbas,ztest)
      if (i.ne.0.or.ztest(1:3).eq.'sto' )go to 8010
      jrec=0
      if (oinvrt.or.onwchem) then
       call inpf(exi(k))
       call inpf(c1)
       call inpf(c2)
      else
       call inpf(c1)
       call inpf(exi(k))
       call inpf(c2)
      endif
      igauss=igauss+1
      if (ytype .eq. ylabel(16)) csinpi(k) = c1
      if (ytype .eq. ylabel(17)) cpinpi(k) = c1
      if (ytype .eq. ylabel(18)) cdinpi(k) = c1
      if (ytype .eq. ylabel(19)) cfinpi(k) = c1
      if (ytype .eq. ylabel(20)) cginpi(k) = c1
      if (ytype .eq. ylabel(21)) csinpi(k) = c1
      if (ytype .eq. ylabel(22)) csinpi(k) = c1
      if (ytype .eq. ylabel(22)) cpinpi(k) = c2
      if (ytype .eq. ylabel(23)) cdinpi(k) = c1
      csi(k) = csinpi(k)
      cpi(k) = cpinpi(k)
      cdi(k) = cdinpi(k)
      cfi(k) = cfinpi(k)
      cgi(k) = cginpi(k)
      k=k+1
      go to 7010
c
 8010 if(igauss.le.0.or.igauss.gt.mxprms)call caserr2(
     *'invalid number of primitives specified in shell')
      kngi(nshi)=igauss
      ngai=ngai+igauss
      if(ngai.gt.ngsmax)go to 1520
      ytype=ytest
      ztype=ztest
      k1=kstari(nshi)
      k2=k1+kngi(nshi)-1
      go to 660
c
c     ----- sto-ng basis set -----
c           j.chem.phys. 38, 2686 (1963)
c                        51, 2657 (1969)
c                        52, 2769 (1970)
c     note - call error routine if atom heavier than zenon
c
  260 kngi(nshi)=igauss
      ngai=ngai+igauss
      k1=kstari(nshi)
      k2=k1+igauss-1
      atnum = czann(nati)
      nucz =  nint(atnum)
      if(nucz.le.0) then
      call inpa(ztagbq)
      nucz=jsubst(ztagbq)
      endif
      call inpf(sc)
      if (nucz .gt. 54) call caserr2(
     *'attempting to site sto function on invalid centre')
      call vclr(exx,1,igauss)
      call vclr(css,1,igauss)
      call vclr(cpp,1,igauss)
      call vclr(cdd,1,igauss)
      go to (300,320,340,360,380,400,420,440,460,480,500,520,540,560,
     +     580,600,600,600,600,600,600,600,600,590),ityp
  300 call s1s(exx,css,igauss)
       stosc = stos1(nucz)
      go to 620
  320 call s2s(exx,css,igauss)
       stosc = stos2(nucz)
      go to 620
  340 call s2p(exx,cpp,igauss)
       stosc = stos2(nucz)
      go to 620
  360 call s2sp(exx,css,cpp,igauss)
       stosc = stos2(nucz)
      go to 620
  380 call s3s(exx,css,igauss)
       stosc = stos3(nucz)
      go to 620
  400 call s3p(exx,cpp,igauss)
       stosc = stos3(nucz)
      go to 620
  420 call s3sp(exx,css,cpp,igauss)
       stosc = stos3(nucz)
      go to 620
  440 call s3d(exx,cdd,igauss)
       stosc = stos4(nucz)
      go to 620
  460 call s4s(exx,css,igauss)
       stosc = stos5(nucz)
      go to 620
  480 call s4p(exx,cpp,igauss)
       stosc = stos5(nucz)
      go to 620
  500 call s4sp(exx,css,cpp,igauss)
       stosc = stos5(nucz)
      go to 620
  520 call s4d(exx,cdd,igauss)
       stosc = stos6(nucz)
      go to 620
  540 call s5s(exx,css,igauss)
       stosc = stos7(nucz)
      go to 620
  560 call s5p(exx,cpp,igauss)
       stosc = stos7(nucz)
      go to 620
  580 call s5d(exx,cdd,igauss)
      go to 620
  590 call s5sp(exx,css,cpp,igauss)
       stosc = stos7(nucz)
      go to 620
  600 call caserr2(
     *'invalid sto function specified')
  620 continue
      if (sc .eq. 0.0d0) sc = stosc
      if (sc .eq. 0.0d0) call caserr2(
     *'scale factor must be supplied for sto function')
      do 640 i = 1,igauss
      k = k1+i-1
      exi(k) = exx(i)*sc**2
      csinpi(k) = css(i)
      cpinpi(k) = cpp(i)
      cdinpi(k) = cdd(i)
      csi(k) = csinpi(k)
      cpi(k) = cpinpi(k)
      cdi(k) = cdinpi(k)
  640 continue
c
c     if(normp.ne.1) ... unnormalization of the primitive functions.
c     if contraction coefficients are given in terms of normalized
c     primitive functions, change them to go with unnormalized
c     primitives .
c     for d shells, the input coefficients cd must be the coefficients
c     corresponding to the normalized primitive x**2 *exp(-a*r**2).
c     for f shells, the input coefficients cf must be the coefficients
c     corresponding to the normalized primitive x**3 *exp(-a*r**2).
c     for g shells, the input coefficients cg must be the coefficients
c     corresponding to the normalized primitive x**4 *exp(-a*r**2).
c
  660 continue
      if (normp .eq. 1) go to 700
      do 680 ig = k1,k2
      ee = exi(ig)+exi(ig)
      facs = pi32/(ee*dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      facf = pt187*facs/(ee**3)
      facg = pt6562*facs/(ee**4)
      csi(ig) = csi(ig)/dsqrt(facs)
      cpi(ig) = cpi(ig)/dsqrt(facp)
      cdi(ig) = cdi(ig)/dsqrt(facd)
      cfi(ig) = cfi(ig)/dsqrt(facf)
      cgi(ig) = cgi(ig)/dsqrt(facg)
  680 continue
  700 continue
c
c     if(normf.ne.1) normalize the contracted basis functions.
c
      if (normf .eq. 1) go to 8100
      facs = 0.0d0
      facp = 0.0d0
      facd = 0.0d0
      facf = 0.0d0
      facg = 0.0d0
      do 740 ig = k1,k2
      do 740 jg = k1,ig
      ee = exi(ig)+exi(jg)
      fac = ee*dsqrt(ee)
      dums = csi(ig)*csi(jg)/fac
      dump = pt5*cpi(ig)*cpi(jg)/(ee*fac)
      dumd = pt75*cdi(ig)*cdi(jg)/(ee**2*fac)
      dumf = pt187*cfi(ig)*cfi(jg)/(ee**3*fac)
      dumg = pt6562*cgi(ig)*cgi(jg)/(ee**4*fac)
      if (ig .eq. jg) go to 720
      dums = dums+dums
      dump = dump+dump
      dumd = dumd+dumd
      dumf = dumf+dumf
      dumg = dumg+dumg
  720 facs = facs+dums
      facp = facp+dump
      facd = facd+dumd
      facf = facf+dumf
      facg = facg+dumg
  740 continue
      do 760 ig = k1,k2
      if (facs .gt. toll) csi(ig) = csi(ig)/dsqrt(facs*pi32)
      if (facp .gt. toll) cpi(ig) = cpi(ig)/dsqrt(facp*pi32)
      if (facd .gt. toll) cdi(ig) = cdi(ig)/dsqrt(facd*pi32)
      if (facf .gt. toll) cfi(ig) = cfi(ig)/dsqrt(facf*pi32)
      if (facg .gt. toll) cgi(ig) = cgi(ig)/dsqrt(facg*pi32)
  760 continue
c
 8100  continue
       if(oend)go to 780
       if(ybasis.eq.ybas1)go to 160
       go to 8013
c
c     generate equivalent centers
c
c ----- first generate equivalent centres
c ----- reordering input shells en route
c
 780  continue


      call atoms2(csinp,cpinp,cdinp,cfinp,cginp,
     * ztita,czana,zatagg,
     *csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     *ict,newsh,iso,isoc,natmax,nshmax,ngsmax,onel)
c
c...  generate extended internal basis and symmetry info for ZORA
c
      if (ozora) call basis_zora(q,csinp,cpinp,cdinp,cfinp,cginp,
     * ztita,czana,zatagg,
     *csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     *ict,newsh,iso,isoc,natmax,nshmax,ngsmax,onel)
c
cjvl
c   store contents of czan in czanr
c   needed since czanr is used as real nuclear charge later
      call dcopy(nat,czan,1,czanr,1)
cjvl
      call chkbas(oignore)
      return
 1500 call prsize(iwr)
      call caserr2(
     *'excessive number of shells')
 1520 call prsize(iwr)
      call caserr2(
     *'excessive number of contracted primitive functions')
 1530 call caserr2('invalid n specified in n-m1g keyword')
 1540 call caserr2('invalid m specified in n-m1g keyword')
 1550 write(iwr,1560)zbasis
 1560 format(1x,'expected keyword in one of the formats: '/,1x,
     +       'n-m1g, n-1m+g, n-1m++g, n-1mg*, n-1m+g*, n-1m++g*, ',
     +       'n-1mg**, n-1m+g**, or n-1m++g**'/,
     +       1x,'but got: ',a8)
      call caserr2('invalid n-m1g keyword')
 9068  format(/1x,104('-')//
     *40x,18('*')/
     *40x,'molecular geometry'/
     *40x,18('*')/)
      return
      end
      subroutine chkbas(oignore)
      implicit none
      logical oignore
c
c     This subroutine checks whether each physical atom has got a
c     basis set assigned to it. If an atom doesn't have a basis set
c     then that is considered to be a fatal error.
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
c
c     local variables
c
      integer iat, atmno, ishl
      logical ofound, oallfound
c
c     functions
c
      integer isubst
      external isubst
c
      oallfound = .true.
      do iat = 1, nat
         atmno = isubst(zaname(iat))
         if (atmno.gt.0) then
c
c           The current "atom" is truly a physical atom.
c
            ofound = .false.
            ishl   = 0
 10         if (ishl.lt.nshell.and..not.ofound) then
               ishl = ishl+1
               ofound = ofound.or.(iat.eq.katom(ishl))
               goto 10
            endif
            oallfound = oallfound.and.ofound
            if (.not.ofound) then
               if (oignore) then
                  write(iwr,60)"WARNING:",iat,zaname(iat)
               else
                  write(iwr,60)"ERROR:  ",iat,zaname(iat)
               endif
            endif
         endif
      enddo
      if (.not.oignore.and..not.oallfound) then
         write(iwr,70)
         call caserr2('No basis set was defined for some atoms')
      endif
 60   format(a8,'No basis set on atom',i6,' called ',a8)
 70   format('If for some reason you want to run the calculation',/,
     +       'anyway then use the "ignore" sub-directive')
      end
      subroutine atosto(ztext,ybas,igauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      integer j
      character*(*) ztext
      dimension xng(6)
      data xblnk,xng/' ',
     *'1','2','3','4','5','6'/
c
      j = 4
      if (ztext(j:j).eq.'-') j = j + 1 
      xmg(1:1)=ztext(j:j)
      do 1 igauss=1,6
      if(xmg.eq.xng(igauss))go to 2
 1    continue
      call caserr2('invalid n specified in stong keyword')
 2    ybas(1:3)=ztext(1:3)
      ybas(4:4)=xblnk
c
      return
      end
      subroutine ato321(ztext,ibasis,igauss,opolh,opolx,
     +                  orydh,orydx,odunn,oerr1,oerr2,oerr3)
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) ztext
      dimension xng(6)
      data xng/
     *'1','2','3','4','5','6'/
      data yblnk,ydunn,yahlr/' ','dunn','ahlr' /
c
      opolh=.false.
      opolx=.false.
      orydh = .false.
      orydx = .false.
      oerr1=.false.
      oerr2=.false.
      oerr3=.false.
      odunn=.true.
      ytext=ytrunc(ztext)
      if(ytext.eq.ydunn.or.ytext.eq.yahlr.or.ytext.eq.yblnk)
     +     go to 5
      xmg(1:1)=ztext(1:1)
      do 1 igauss=1,6
         if(xmg.eq.xng(igauss))go to 2
 1    continue
      oerr1=.true.
      return
 2    if (ztext(2:2).ne.'-') then
         oerr3 = .true.
      endif
      xmg(1:1)=ztext(3:3)
      do 3 ibasis=2,3
         if(xmg.eq.xng(ibasis))go to 4
 3    continue
      oerr2=.true.
      return
 4    ibasis=ibasis-2
c
c     now check for any of the following syntaxes 
c     N-M1+G, N-M1++G, N-M1+G*, N-M1++G*, N-M1++G**, N-M1++G**
c     N-M1G*, N-M1G**
c
c     first diffuse functions (+ or ++)
c
      ipol = 6
      if(ztext(5:5).eq.'+') then
        orydx=.true.
        ipol = ipol + 1
        if(ztext(6:6).eq.'+') then
          orydh=.true.
          ipol = ipol + 1
          if (ztext(7:7).ne.'g') then
            oerr3 = .true.
          endif
        else if (ztext(6:6).ne.'g') then
          oerr3 = .true.
        endif
      else if (ztext(5:5).ne.'g') then
        oerr3 = .true.
      endif
      if(ztext(ipol:ipol).eq.'*') then
         opolx=.true.
         if(ztext(ipol+1:ipol+1).eq.'*') then
           opolh=.true.
         else if (ztext(ipol+1:ipol+1).ne.' ') then
            oerr3 = .true.
         endif
      else if (ztext(ipol:ipol).ne.' ') then
         oerr3 = .true.
      endif
      return
 5    if(ytext.eq.yahlr) odunn = .false.
      ibasis=2
c
      return
      end
      subroutine atoms2(csinp,cpinp,cdinp,cfinp,cginp,
     + ztita,czana,zatagg,
     + csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     + ict,newsh,iso,isoc,natmax,nshmax,ngsmax,onel)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c ktype is nangm(ityp)
c nbfs(ityp) is number of functions in the shell
c ylabel(ityp) see in block data mains/all.f
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
      parameter (maxgrorder=48)
      parameter(ishthr=100)
c
      real*8 degecr
      integer iscsym
      logical otsym,oingam,omydbg,ostart
      common/fsymas/degecr,otsym,oingam,omydbg,ostart,iscsym
c
c
      real*8 toler, toler2
      integer isymtl
      common /tol/ toler,toler2,isymtl
c
      common/junkc/ ylabel(26)
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
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508)
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
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
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      common/transf/psmal,qsmal,rsmal,pnew,qnew,rnew,pp,qp,rp
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
c     ZORA common
c
      logical ozora,oscalz,onlyp,onlyd,onlyf,onlyg,opzora,
     1        op2zora,onlys,ocontr,oint_zora,o1e_zora,osmall_zora,
     2        opre_zora,oso,onoext_z,oioraz,oatint_z,oatscf_z,
     3        oscalatz
      real*8 cspeed,critat_z
      integer nbas_zora,ibls_z,iblv_z,nwv_z,ibiso_z,nwiso_z,ibcor_z,
     1        ibsx_z,icoul_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     2        ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,ibscale_z,ibcoul_z,
     3        niter_z,int_zora,num_ext,igauge_z,nwcor_z,isexch_z,
     4        nat_z,ibdens_z,nwdens_z,ibscalao_z,ibshat_z,is_z,
     5        irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,ibl7ew_z
      common/zorac/ cspeed,critat_z,ozora,oscalz,
     1              onlyp,onlyd,onlyf,onlyg,opzora,
     1              op2zora,onlys,ocontr,nbas_zora,
     2              oint_zora,o1e_zora,icoul_z,
     3              ibls_z,iblv_z,nwv_z,
     3              ibiso_z,nwiso_z,ibcor_z,nwcor_z,
     4              ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z,
     5              ibtrin_z,ibtrout_z,ibscalf_z,nwtrat_z,
     6              ibscale_z,ibcoul_z,ibdens_z,nwdens_z,osmall_zora,
     7              niter_z,opre_zora,int_zora,num_ext,isexch_z,
     8              oso,onoext_z,is_z,igauge_z,oioraz,
     9              nat_z,ibscalao_z,oatint_z,oatscf_z,ibshat_z,
     1              oscalatz,irest_z,numdu_z,ibldu_z,iblock_z,lshift_z,
     2              ibl7ew_z
c....    icoul_z : 1 : Full coulomb
c...               2 : atomic coulomb (default)
c...               3 : no coulomb (just 1-electron ints)
c...               4 : small basis coulomb
c...     isexch_z  0 : no exchange in zora correction (default)
c...               1 : exchange in zora correction (like in dft)
c...     niter_z : # iterations between coulomb matrix update
c...     int_zora :internal basis gen. : 0 none, 1, automatic, 2 by hand
c
c...     oscalz  : scaling on or off (true or false)
c...     oioraz  : change metric to perform inf. order scaling 
c...               (Dyall/EvLenthe (not quite))
c...     oscalatz: atoms only scaled zora (together with get atom)
c...     only'x' : internal basis only op to 'x'-type functions
c...     op(2)zora: print flags (op2zora+opzora: prints *everything*)
c...     ocontr  : adds l+1 block to internal basis in contracted way 
c...     nbas_zora: size of internal basis
c...     oint_zora: in internal or external mode???
c...     o1e_zora: one electron integrals have been calculated in *this*
c...               internal basis
c...     osmall_zora: small swap for projected coulomb
c...     opre_zora: pre_zora has been called
c...     oso     : spin orbit calculation 
c...     onoext_z: no external functions copied to internal basis
c...     is_z : section for orbitals used for  coulomb matrix
c...     igauge_z : add so many s-marices to h-matrix
c...                the energy should change by igauge
c...     nat_z  : indicates (if ne 0) that starting density is used
c...              to calc. coulomb; gives # it. to get atomic density
c...     critat_z : accuracy of atomic energy for atomic zora
c...     oatint_z : indicates that integrals to be calculated are atomic
c...     oatscf_z : indicates we are doing a set of atomic SCFs
c...     ibshat_z : block on num8 to store blocked s/h
c...     for disk-adress and lengths (ib.., nw.. ) see subroutine atoms2
c...     irest_z  : .ne.0 indicates atomic zora corrections from dumpfile; 1 not on restart
c...     numdu_z,ibldu_z : if ne 0 other dumpfile to read zora corrections from
c...     iblock_z is the block number on ed3, where the zora corrections (section 425) reside
c...     ibl7ew is where the energy weighted density matrix resides
c
c...  counterpoise common
c
      character *8 ztag00
      integer maxghost,nghost
      parameter (maxghost=100)
      common/ghost/ ztag00(maxghost)
      common/ghostn/ nghost
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
      common/junk/trmat(432),
     +   ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +   intyp(mxshel),ns(maxat),ks(maxat),
     +   intypi(mxshel),nsi(maxat),ksi(maxat),
     +   exi(mxprim),
     +   csi(mxprim),cpi(mxprim),cdi(mxprim),cfi(mxprim),cgi(mxprim),
     +   kstari(mxshel),katomi(mxshel),ktypi(mxshel),kngi(mxshel),
     +   kloci(mxshel),kmini(mxshel),kmaxi(mxshel),nshi
c
      real*8 stos1, stos2, stos3, stos4, stos5, stos6, stos7
      real*8 rleesf
      real*8 c11al, c11be, c31al, c31be, c33al, c33be
      integer nbfs, minf, maxf, nangm
c
      common /datbas/ stos1(54),stos2(54),stos3(54),stos4(54),
     + stos5(54),stos6(54),stos7(54),
     + rleesf(36,2),
     + nbfs(24),minf(24),maxf(24),nangm(24),
     + c11al(12),c11be(12),c31al(12),c31be(12),c33al(12),c33be(12)
c

      dimension kmytyp(mxshel)
      dimension ict(natmax,*),newsh(nshmax,*),iso(*),isoc(*)
      dimension ztita(*),czana(*),zatagg(*)
      dimension ytemp(mxshel)
      dimension csinpi(ngsmax),cpinpi(ngsmax),cdinpi(ngsmax)
      dimension cfinpi(ngsmax),cginpi(ngsmax)
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax)
      dimension cfinp(ngsmax),cginp(ngsmax)
      
      data toll/1.0d-10/
      data zline,zlinev /' * * * *','       *'/
      data tm4/1.0d-04/
      data m1,m18/1,18/
c
      m1420=5*mxprim+maxat
      m110 =10+maxat
      rne = 0.0d0
      if (nghost.gt.0) call ghosts
c
c     generate equivalent centers
c     set table ( centers versus transformations )
c     set table ( shells versus transformations )


c
c ----- first generate equivalent centres
c ----- reordering input shells en route
c
      loc=0
      ngauss=0
      nshell=0
c...  initialise isptype to an illegal value so we can trap failures
c...  to initialise it properly... (hvd)
      isptype=-1
      nat=0
c
      do 8190 i=1,nonsym
 8190 map80(i)=0
      nav = lenwrd()
c
      do 9000 i=1,nonsym
      if(map80(i).eq.1)then
      go to 9000
      endif
      nat=nat+1
      p=cat(1,i)
      q=cat(2,i)
      r=cat(3,i)
      test=dsqrt(p*p+q*q+r*r)
      if(test.ge.toler)go to 7777
      p=0.0d0
      q=0.0d0
      r=0.0d0
7777  rne=rne+czann(i)
      c(1,nat)=p
      c(2,nat)=q
      c(3,nat)=r
      czan(nat)=czann(i)
      nsnat=0
      zaname(nat)=ztag(i)
c
c-weights now see function amass_get for nuclear masses
c
c      amas(nat)=amass(i)

c
c ----- abstract and order shells on unique centre i
c
      ks(nat)=nshell+1
      do 9001 j=1,nshi

      zagat=ztag(katomi(j))
      if(zagat.ne.ztag(i))go to 9001

      nshell=nshell+1
      nsnat=nsnat+1
      ityp=intypi(j)
      intyp(nshell)=ityp
      kmin(nshell)=kmini(j)
      kmax(nshell)=kmaxi(j)
      katom(nshell)=nat
      ktype(nshell)=ktypi(j)
      igauss=kngi(j)
      kng(nshell)=igauss
      kstart(nshell)=ngauss+1
      kloc(nshell)=loc+1
      loc=loc+nbfs(ityp)
c
      l=kstari(j)
c
      do 9002 k=1,igauss
      kk=k+ngauss
      ex(kk)=exi(l)
      cs(kk)=csi(l)
      cp(kk)=cpi(l)
      cd(kk)=cdi(l)
      cf(kk)=cfi(l)
      cg(kk)=cgi(l)
      csinp(kk)=csinpi(l)
      cpinp(kk)=cpinpi(l)
      cdinp(kk)=cdinpi(l)
      cfinp(kk)=cfinpi(l)
      cginp(kk)=cginpi(l)
 9002 l=l+1
c
      ngauss=ngauss+igauss
      if(ngauss.gt.ngsmax)go to 1520
 9001 continue
c
      ns(nat)=nsnat
c
      call lframe(p,q,r,ps,qs,rs)
      psmal = ps
      qsmal = qs
      rsmal = rs
      nat0 = nat
      do 840 it = 1,nt
c  skip identity
      if (it .eq. 1) go to 840
      nn = 9*(it-1)
      call trans1(nn)
      call rot
c first test only in the group of already processed atoms -
c in array c - nat and not nonsym
      do 800 iat = 1,nat
      test = (pp-c(1,iat))**2+(qp-c(2,iat))**2+(rp-c(3,iat))**2
      if (test .le. toler) then
      go to 840
      endif
  800 continue
c
c ----- check for atom in cat
c
      do 810 k=1,nonsym
      if(map80(k).eq.1)go to 810
      test = (pp-cat(1,k))**2
     *      +(qp-cat(2,k))**2
     *      +(rp-cat(3,k))**2
      if(test.gt.toler)go to 810
c flag that this one was already incorporated
c in the list i array c
      map80(k)=1
      go to 812
 810  continue
      write(iwr,811)i,it,toler,pp,qp,rp
c
c i think it should not be added even when using xyz - perhaps
c it could help sometimes, but usually then nat will .ne. nonsym
c
c     if(ozmat)go to 840
      goto 840
c 
c
c add the newly found center
 812  nuat = nat
      nat = nat+1
      if (nat .gt. natmax) go to 1480
      c(1,nat) = c(1,nuat)
      c(2,nat) = c(2,nuat)
      c(3,nat) = c(3,nuat)
      c(1,nuat) = pp
      c(2,nuat) = qp
      c(3,nuat) = rp
c copy all properties of the atom from the origin
      ns(nat) = ns(nat0)
      ks(nat) = ks(nat-1)+ns(nat0)
      zaname(nat) = zaname(nat0)
      czan(nat) = czan(nat0)
c
c-weights
c
c  now see amass_get for nuclear masses
c
c      amas(nat)=amas(nat0)
c
      rne = rne+czan(nat)
      nshell = nshell+ns(nat)
      if (nshell .gt. nshmax) go to 1500
      ns1 = ns(nat)
      if(ns1.eq.0)go to 840
      do 820 k = 1,ns1
      j = ks(nat0)+k-1
      jj = ks(nat)+k-1
      kmin(jj) = kmin(j)
      kmax(jj) = kmax(j)
      kstart(jj) = kstart(j)
      ktype(jj) = ktype(j)
      kng(jj) = kng(j)
      katom(jj) = nat
      kloc(jj) = loc+1
      intyp(jj) = intyp(j)
  820 loc = loc+nbfs(intyp(jj))
  840 continue
c
 9000 continue
c
c     if(normp.ne.1) the contracted basis functions have been
c     expressed in terms of unnormalized primitive functions.
c
      if (normp .ne. 1 . and . oprint(48)) write (iwr,9408)
c
c     if(normf.ne.1) the contracted basis functions have been
c     normalized to unity
c
      if (normf .ne. 1 . and . oprint(48)) write (iwr,9508)
c
c     form transformation tables for atoms and shells.
c
c over all atoms
c
      if(omydbg)then
      write(iwr,*)' atoms2: calculating transformation table of ',
     +      nat,'   atoms'
      write(iwr,*)' new array of coordinates'
      do 945 jj=1,nat
      write(iwr,998)(c(ii,jj),ii=1,3)
945   continue
      endif
      do 960 iat = 1,nat
      p = c(1,iat)
      q = c(2,iat)
      r = c(3,iat)
      ns1 = ks(iat)-1
      ns2 = ns(iat)
      call lframe(p,q,r,ps,qs,rs)
      psmal = ps
      qsmal = qs
      rsmal = rs
c over all group elements
      do 940 it = 1,nt
      nn = 9*(it-1)
      call trans1(nn)
      call rot
c over all atoms, first assign 0 to see if no position was found
c and to detect later an error
      ictr=0
      do 880 i = 1,nat
      test = (pp-c(1,i))**2+(qp-c(2,i))**2+(rp-c(3,i))**2
      if (test .gt. toler) go to 880
      ictr = i
      go to 900
  880 continue
  900 ict(iat,it) = ictr
      if (omydbg) then
      write(iwr,999)iat,it,pp,qp,rp,ictr
      endif
c permutation of atoms generated
      ns3 = ks(ictr)-1
      if (ns2 .eq. 0) go to 940
      do 920 ish = 1,ns2
  920 newsh(ns1+ish,it) = ns3+ish
  940 continue
c permutation of shells generated
  960 continue
c
c write transf of atoms and shells
c
      if(otsym) then
        ocant=.false.
        do 5003,i=1,nshell
        kmytyp(i)=ktype(i)
        inti=intyp(i)
c see block data mains/all.f for these magic numbers
        if(inti.eq.4.or.inti.eq.7.or.inti.eq.11.or.inti.eq.22.or.
     +     inti.eq.24) kmytyp(i)=ishthr+1
        if(inti.eq.23.or.inti.ge.25) then
        ocant=.true.
        kmytyp(i)=ishthr+2
        endif
5003    continue
        if(ocant) then
        write(iwr,5004)
        otsym=.false.
        endif
        if (.not.oint_zora) call putsym(iscsym,nat,nshell,loc,nt,
     +   kloc,katom,kmytyp,ict,newsh,natmax,nshmax)
      endif
c
      ns1=0
      iat=0
      do 2040 i=1,48
      iliso(i)=ns1
      ilisoc(i)=iat
      ns1=ns1+nshell
 2040 iat=iat+nat
c
      do 1040 i = 1,nshell
      ns1 = 0
        do 1000 it = 1,nt
           iso(i+ns1) = newsh(i,it)
           ns1 = ns1 + nshell
 1000   continue
 1040 continue
c
      do 1120 i = 1,nat
      ns1 = 0
         do 1080 it = 1,nt
         isoc(i+ns1) = ict(i,it)
         ns1 = ns1 + nat
 1080    continue
 1120 continue
c
c ----- if zmat invoked, deduce g80->games mapping vector
c
      if(ozmat)then
         do 6400 i=1,nat
            map80(i)=0
            p=cnew(i,1)
            q=cnew(i,2)
            r=cnew(i,3)
            do 6401 j=1,nat
               rij=(p-c(1,j))**2 + (q-c(2,j))**2 + (r-c(3,j))**2
               if( dabs(rij).gt.toler)go to 6401
c              check we have not used this atom before
               do k=1,i-1
                  if (map80(k).eq.j) goto 6401
               enddo
               map80(i)=j
               go to 6400
 6401       continue
 6400    continue
c
c apply atom permutation to tables of isotopic substitution
c patterns
c
         call mass_permute(map80)
c
      else
c
c apply atom permutation to tables of isotopic substitution
c patterns based on coordinates in cat
c
         do 6406 i=1,nat
            map80(i)=0
            p=cat(1,i)
            q=cat(2,i)
            r=cat(3,i)
            do 6408 j=1,nat
               rij=(p-c(1,j))**2 + (q-c(2,j))**2 + (r-c(3,j))**2
               if( dabs(rij).le.toler) then
c                 check we have not used this atom before
                  do k=1,i-1
                     if (map80(k).eq.j) goto 6408
                  enddo
                  map80(i)=j
                  goto 6406
               endif
 6408       continue
 6406    continue
         call mass_permute(map80)
      endif
c
c     ----- evaluate length of section 496
c
      if (.not.oint_zora) then
         mword1=(nshell*nt-1)/nav+1
         ibl1=lensec(mword1)
         mword2=(nat*nt-1)/nav+1
         ibl2=lensec(mword2)
         ibl3=lensec(1728)
         ibl4=lensec(4800)
         ibl5=lensec(10800)
c
         nw196(1)=432
         nw196(2)=1728
         nw196(3)=4800
         nw196(4)=10800
         nw196(5)=mword1
         nw196(6)=mword2
c
         len196=2+ibl1+ibl2+ibl3+ibl4+ibl5
c
         call secput(isect(496),m18,len196,ibl496)
c
         ibl196(1)=ibl496+1
         ibl196(2)=ibl196(1)+1
         ibl196(3)=ibl196(2)+ibl3
         ibl196(4)=ibl196(3)+ibl4
         ibl196(5)=ibl196(4)+ibl5
         ibl196(6)=ibl196(5)+ibl1
c
c     ----- if restart job, read in actual coordinates -----
c
         if (lcoord .lt. 1) go to 1160
         call rdrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
         if (lcoord.eq.2) go to 1160
         if (nprint .ne. 7) go to 2000
         do 1140 i = 1,nat
 1140    write (ipu,9028) zaname(i),czan(i),c(1,i),c(2,i),c(3,i)
         go to 2000
c
c  ----- write -isoc-  -iso- -ptr- -dtr- -ftr- and -gtr- to section 496
c
1160     nwnt = 134
         call wrt3i(invt,nwnt,ibl496,idaf)
         call wrt3(ptr,nw196(1),ibl196(1),idaf)
         call wrt3(dtr,nw196(2),ibl196(2),idaf)
         call wrt3(ftr,nw196(3),ibl196(3),idaf)
         call wrt3(gtr,nw196(4),ibl196(4),idaf)
         call wrt3i(iso,nw196(5)*nav,ibl196(5),idaf)
         call wrt3i(isoc,nw196(6)*nav,ibl196(6),idaf)
      endif
c
c     nprint=7   ... punch basis set input data + mo*s
c     nprint=1   ... extra printing for basis set + symmetry.
c     nprint=0   ... normal printing
c
2000  continue
      if (nprint .eq. 6. or.(oint_zora.and..not.opzora)) go to 1540
      write (iwr,9088)
      nbqnop = 0
      do 1180 iat = 1,nat
      if (.not. oprint(31) .or. (zaname(iat)(1:2) .ne. 'bq'))then
      write (iwr,9108)zaname(iat),czan(iat),c(1,iat),c(2,iat),
     * c(3,iat),ns(iat)
      ns2=ns(iat)
      if(ns2.eq.0)go to 1180
      i=(ns2-1)/3+1
      i3=3*i
      do 7401 j=1,i3
7401  ytemp(j)=' '
      i=0
      ns1 = ks(iat)
      ns2 = ns1+ns2-1
      do 7402 ish=ns1,ns2
      i=i+1
 7402 ytemp(i)=ylabel(intyp(ish))
      write (iwr,9268) (ytemp(ish),ish=1,i3)
      else
        nbqnop = nbqnop + 1
      endif
 1180 continue
      if (nbqnop .gt. 0)then
         write (iwr,9109) nbqnop
      endif

      write(iwr,9129)
c noprint basis
      if(oprint(22))go to 2100
      if (oint_zora) then
         write(iwr,9130)
      else
         write (iwr,9128)
      endif
      do 1360 iat = 1,nat
      do 1200 it = 1,nt
      if (ict(iat,it) .gt. iat) go to 1360
 1200 continue
c this makes skip over all atoms which are sym equiv to some later atom
c so that basis only over nonequivalent atoms is printed
      ns2=ns(iat)
      if(ns2.eq.0)go to 1360
      write (iwr,9548) zaname(iat)
      ns1 = ks(iat)
      ns2 = ns1+ns2-1
c
      if (ns2.gt.mxshel) then
         call prsize(iwr)
         call caserr2('maximum number of shells exceeded')
      end if
c
      do 1340 ish = ns1,ns2
      write (iwr,9468)
      i1 = kstart(ish)
      i2 = i1+kng(ish)-1
      ityp = intyp(ish)
      do 1320 ig = i1,i2
      go to (1220,1220,1240,1280,1220,1240,1280,1260,1220,1240,1280,
     +     1260,1220,1240,1260,1220,1240,1260,1250,1270,
     +     1220,1280,1260,1280),ityp
 1220 c1 = cs(ig)
      c2 = csinp(ig)
      go to 1300
 1240 c1 = cp(ig)
      c2 = cpinp(ig)
      go to 1300
 1260 c1 = cd(ig)
      c2 = cdinp(ig)
      go to 1300
 1250 c1 = cf(ig)
      c2 = cfinp(ig)
      go to 1300
 1270 c1 = cg(ig)
      c2 = cginp(ig)
      go to 1300
 1280 c1 = cs(ig)
      c2 = csinp(ig)
      c3 = cp(ig)
      c4 = cpinp(ig)
c this is the main output of basis - composite shells
c check for zero values (probably programming error!)
      write (iwr,9528) ish,ylabel(ityp),ig,ex(ig),c1,c2,c3,c4
c
      if (dabs(ex(ig)).lt.toll .or. dabs(c1).lt.toll .or.
     +    dabs(c2).lt.toll .or. dabs(c3).lt.toll .or.
     +    dabs(c4).lt.toll) then
       call caserr2('zero exponent or coefficient detected')
      endif
      go to 1320
c
c this is the main output of basis - simple shells
c
 1300 write (iwr,9148) ish,ylabel(ityp),ig,ex(ig),c1,c2
      if (dabs(ex(ig)).lt.toll .or. dabs(c1).lt.toll .or.
     +    dabs(c2).lt.toll ) then
       call caserr2('zero exponent or coefficient detected')
      endif
 1320 continue
 1340 continue
 1360 continue
      write(iwr,9700)
2100  if (.not.oprint(48)) go to 1540
      write (iwr,9428)
      write (iwr,9288)
      imax = 0
 1380 imin = imax+1
      imax = imax+15
      if (imax .gt. nt) imax = nt
      imax1 = imax+1
      write (iwr,9308)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9348) zlinev,(i,i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
      do 1400 iat = 1,nat
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9368) iat,(ict(iat,i),i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
 1400 continue
      if (imax .lt. nt) go to 1380
      write (iwr,9428)
      write (iwr,9388)
      imax = 0
 1420 imin = imax+1
      imax = imax+15
      if (imax .gt. nt) imax = nt
      imax1 = imax+1
      write (iwr,9308)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9348) zlinev,(i,i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
      do 1440 ish = 1,nshell
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9368) ish,(newsh(ish,i),i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
 1440 continue
      if (imax .lt. nt) go to 1420
      write (iwr,9428)
      write (iwr,9488)
      imax = 0
 1460 imin = imax+1
      imax = imax+15
      if (imax .gt. nt) imax = nt
      imax1 = imax+1
      write (iwr,9308)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9348) zlinev,(i,i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9348) zlinev,(invt(i),i = imin,imax)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zlinev,i = imin,imax1)
      write (iwr,9328) (zline,i = imin,imax1)
      if (imax .lt. nt) go to 1460
      go to 1540
 1480 call caserr2(
     *'excessive number of atoms')
 1500 call caserr2(
     *'excessive number of shells')
 1520 call caserr2(
     *'excessive number of contracted primitive functions')
 1540 continue
c
c here is a good place to test the transformation table of atoms
c
      do 955 ii=1,nat
      do 965 jj=1,nt
       if(ict(ii,jj).eq.0) then
         write(iwr,9569)ii,jj,ict(ii,jj)
         write(iwr,9570) 
         call caserr2('error in symmetry package')
         endif
 965  continue
 955  continue
c
c check also number of atoms inputted and generated by symmetry
c
      if(nonsym.ne.nat)then
      write(iwr,*)' atoms2: nonsym,nat:',nonsym,nat
      call caserr2(
     +  'symm.- no. of generated atoms differ from no. input')
      endif
      num = loc

      if (onel) then
c
c      make sure the charge is set to a value consistent with the
c      number of electrons
c
       ich = nint(rne) - ne
      else
       ne = nint(rne)
       testne = rne - ne
       if ( dabs(testne) .gt. tm4) then
        write (iwr,9568) rne,1.0d0*ne,testne
        call caserr2('non-integral charge detected')
       endif
      endif

      if(.not.ozmat)go to 6404
c
c ----- if zmat invoked,print g80->gamess atom mapping vector
c
      if(oprint(48))write(iwr,6402)
      if(oprint(48))write(iwr,6403)(map80(i),i=1,nat)
c
c ----- classify nucleii into types ...
c
 6404 call build_nuct(nat,ns,zaname,czan,nuct)
c
c ----- output basis specification
c
c     if (.not.oint_zora) then
      if(irest.ge.1)go to 850
      non=nat
      numorb=num
      call dcopy(nat,czan,1,czana,1)
      do 6409 i=1,nat
 6409 zatagg(i)=zaname(i)
      do 6410 i=1,10
 6410 ztita(i)=ztitle(i)
c
      if (.not.oint_zora) then
         i=lensec(mxprim)+lensec(m110)+lensec(m1420)+lensec(mach(2))
         call secput(isect(491),m1,i,ibl491)
         call wrt3(ex,ngsmax,ibl491,idaf)
         call wrt3s(csinp,m1420,idaf)
         call wrtcs(ztita,m110,idaf)
         call wrt3is(kstart,mach(2)*nav,idaf)
      endif
c
 850  if(.not.onel)ne = ne-ich
c
c     Work out the number of alpha (na) and beta (nb) electrons from
c     the total number (ne) of electrons and the spin multiplicity. The
c     rules are:
c     1. na >= nb
c     2. na+nb = ne
c     The following few lines guarantee that these rules are met in
c     all cases.
c
      nb = (ne-mul+1)/2
      na = ne - nb
c
      nx = num*(num+1)/2
c
c...  now check for overflow of specifications in sizes
c     and multiplicity specification
c
      if (nshell.gt.mxshel.or.num.gt.maxorb.or.nat.gt.maxat) 
     +   call prsize(iwr)
      if (nshell.gt.mxshel) call caserr2('Too many shells cf. sizes')
      if (num.gt.maxorb)
     +   call caserr2('Too many basis functions cf. sizes')
      if (nat.gt.maxat) call caserr2('Too many atoms cf. sizes')
c
      if (.not. oint_zora) call preharm
c
      if ((.not.oharm.and..not.oint_zora).or.(oint_zora.and.opzora))
     +   write (iwr,9168) nshell,num,ne,ich,mul,na,nb,nat
      if (oharm.and..not.oint_zora) 
     1    write(iwr,9169) nshell,num,newbash,ne,ich,mul,na,nb,nat
      if (oint_zora.and..not.opzora)
     +   write (iwr,9170) nshell,num
 9700 format(1x,113('='))
c
      if (oint_zora) then
c
c...  extra zora info is dumped to ed7 (num8)
c...  make disk layout
c...  ibls_z,iblv_z,nwv_z : blocks and # elements for s,v (intern basis)
c...  ibiso_z,nwiso_z : block / # words for iso ( internal basis)
c...  ibcor_z,nwcor_z : block / # words for H-correction (extern)
c...  ibsx_z,ibsy_z,ibsz_z,ibsmin_z,ibsplus_z,nwmat_z: blocks + length
c...       for x,y,z,s**-0.5, s+**0.5, in internal basis
c...  ibtrin_z, ibtrout_z, nwtrat_z : transformation internal to extern
c...  ibcoul_z : coulomb matrix in internal basis (nwv_z)
c...  ibdens_z : density matrix (external basis) used in zora (nwdens_z)
c...  ibscalf_z: scale factors (orbital expectations values)  (num)
c...  ibscalao_z: above scale factors in ao-basis
c...  all lengths in real*8 words
c
c     ----- evaluate lengths  (in 8-byte words)
c
         nav = lenwrd()
         nwiso_z = (nshell*nt-1)/nav+1
         nwcor_z = num_ext*(num_ext+1)/2
         nwv_z = num*(num+1)/2
         nwmat_z = num*num
         nwtrat_z = num*num_ext
         nwdens_z = num_ext*(num_ext+1)/2
c
         ibiso_z = ibl7la
         ibcor_z = ibiso_z + lensec(nwiso_z)
         ibls_z = ibcor_z + lensec(nwcor_z)
         iblv_z = ibls_z + lensec(nwv_z)
c
         ibsx_z = iblv_z + lensec(nwv_z)
         ibsy_z = ibsx_z + lensec(nwmat_z)
         ibsz_z = ibsy_z + lensec(nwmat_z)
         ibsmin_z = ibsz_z + lensec(nwmat_z)
         ibsplus_z = ibsmin_z + lensec(nwmat_z)
         ibtrin_z = ibsplus_z + lensec(nwmat_z)
         ibtrout_z = ibtrin_z + lensec(nwtrat_z)
c
         ibscale_z = ibtrout_z + lensec(nwtrat_z)
         ibscalf_z = ibscale_z + lensec(nwcor_z)
         ibscalao_z = ibscalf_z + lensec(num_ext)
         ibcoul_z =  ibscalao_z + lensec(nwcor_z)
         ibdens_z = ibcoul_z + lensec(nwv_z) 
         ibl7la = ibdens_z + lensec(nwdens_z)
c
         if (oscalatz) then
            ibshat_z = ibl7la
            ibl7la = ibshat_z + lensec(2*nwcor_z)
         end if
c
c
c
c  ----- write -isoc-  -iso-  to ed7
c   (-ptr- -dtr- -ftr- and -gtr-  are written to section 496 in atoms)
c
         call wrt3i(iso,nwiso_z*nav,ibiso_z,num8)
c
         o1e_zora = .false.
      endif
c
      return
 9028 format(a8,2x,f5.0,3f20.12)
 9088 format(9x,79('*')/9x,'*',77x,'*'/
     *9x,'*',5x,'atom',3x,'atomic',16x,'coordinates',17x,'number of',
     *6x,'*'/
     *9x,'*',12x,'charge',7x,'x',13x,'y',14x,'z',7x,'shells',9x,'*'/
     *9x,'*',77x,'*'/9x,79('*'))
 9108 format(9x,'*',77x,'*'/9x,'*',77x,'*'/
     *        9x,'*',4x,a8,f5.1,3(f12.7,3x),i5,10x,'*')
 9109 format (9x,'*   Output of ',i5,' BQ centres suppressed',36x,'*')
 9129 format(9x,'*',77x,'*'/9x,79('*')/)
 9128 format(/40x,19('*')/
     *40x,'molecular basis set'/
     *40x,19('*')///40x,30('=')/
     +40x,'contracted primitive functions'/40x,30('=')//
     +1x,'atom',8x,'shell',3x,'type',2x,'prim',7x,'exponents',
     +12x,'contraction coefficients'/1x,113('='))
 9130 format(/40x,23('*')/
     *40x,'ZORA internal basis set'/
     *40x,23('*')///40x,30('=')/
     +40x,'contracted primitive functions'/40x,30('=')//
     +1x,'atom',8x,'shell',3x,'type',2x,'prim',7x,'exponents',
     +12x,'contraction coefficients'/1x,113('='))
 9148 format(13x,i5,3x,a4,1x,i5,1x,2f15.6,3h  (,f12.6,3h  ))
 9168 format(/
     *' total number of shells',15x,i5/
     +' total number of basis functions',6x,i5/
     +  ' number of electrons',18x,i5/
     + ' charge of molecule',19x,i5/
     + ' state multiplicity',19x,i5/
     + ' number of occupied orbitals (alpha)',2x,i5/
     +    ' number of occupied orbitals (beta )',2x,i5/
     +    ' total number of atoms',16x,i5/)
 9169 format(/
     *' total number of shells',22x,i5/
     +' total number of cartesian basis functions',3x,i5/
     +' total number of harmonic basis functions',4x,i5/
     +' number of electrons',25x,i5/
     +' charge of molecule',26x,i5/
     +' state multiplicity',26x,i5/
     +' number of occupied orbitals (alpha)',9x,i5/
     +' number of occupied orbitals (beta )',9x,i5/
     +' total number of atoms',23x,i5/)
 9170 format(' **** ZORA Internal Basis ****',/,
     *' total number of shells',15x,i5/
     +' total number of basis functions',6x,i5)
 9268 format(9x,'*',64x,3a4,1x,'*')
 9288 format(/,1x,5h*****,30h transformation table of atoms, 5h*****,/,
     +     30x,15h rows are atoms,/,30x,
     +     32h columns are symmetry operations)
 9308 format(//)
 9328 format(1x,16a8)
 9348 format(1x,a8,15(2x,i2,3x,1h*))
 9368 format(1x,16(2x,i2,3x,1h*))
 9388 format(/1x,'***** transformation table of shells *****'/
     +       30x,'rows are shells'/
     +       30x,'columns are symmetry operations')
 9408 format(//' contracted primitives have been unnormalised')
 9428 format(1h0)
 9468 format(/)
 9488 format(/1x,'***** inverse transformations *****'/)
 9508 format(//' contracted basis functions are normalized to unity'/)
 9528 format(13x,i5,3x,a4,1x,i5,1x,2f15.6,3h  (,f12.6,3h  ), f15.6,
     +     3h  (,f12.6,3h  ))
 9548 format(/1x,a8)
 9568 format(1x,'non-integral charge found',1x,3e15.8)
 9569 format(1x,'atoms2: atom,operation,picture:',3i7)
 9570 format(1x,
     +'There is an inconsistency in symmetry package - symmetry '/1x,
     +'molecule was not properly recognized or the orientation'/1x,
     +'is wrong or the group operations are incorrectly generated'//1x,
     +'Maybe the molecule has accidentally equal two main moments'/1x,
     +'of inertia (or charge) - try to change second parameter after'/1x
     +,'symm directive or perturb your geometry a little.'//1x,
     +'If this does not help, rerun with iprint symm syma and'/1x,
     +'check the transformation table of atoms - 0 means error'/1x,
     +'- atom does not have symmetry counterpart.'//1x,
     +'This could be caused by an bug in generator of group'/1x,
     +'operations - check transformation matrices and try to'/1x,
     +'fix the bug or simply decrease artificially your symmetry'/1x,
     +'to conduct calculation in a different group'//)
811   format(/1x,53('*'),' warning ',55('*')/
     *1x,'no symmetry related atom for picture of atom ',
     *i3,' in new list under operation',i3/
     *1x,'tolerance = ',f15.13,' x y z coords =',3f12.8/
     *1x,'this can happen if the epsilon for test of positions'/
     +1x,'is too small'/
     *1x,'but it can mean also a bug in the symmetry package'/
     *1x,117('*')/)
998   format(3f15.8)
999   format('atom',i3,'  oper',i3,'  xyz',3f15.8,'   found',i4)
5004  format(/
     +  ' WARNING: symmetry assignment switched off because of shell'/
     +  ' type which is not recognized by the symass'/)
6402  format(//' gaussian80-gamess atom mapping vector'/)
6403  format(/10x,20i4)
      end
      subroutine prsize(iwr)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c...  print sizes
c...  could be printed as special debug 
c...  makes version identification easier
c
      write(iwr,600) maxorb,maxat,maxvar,maxnz,mxshel,mxprim,
     1               mxprms,mxgaus,mxgrps,
     2               mxgrid, mxcalc,mxplot,mxrest,
     3               maxlfn,maxfrt,maxbuf,maxblo,
     4               nd200,mxcan1,mxcan2,
     5               lenci
600   format(//1x,
     + ' =======  parameters defining the maximum system size ======',/,
     + ' maxorb = maximum number of basis functions           ',i6,/,
     + ' maxat  = maximum number of atoms (incl. point charges',i6,/,
     + ' maxvar = maximum number of z-matrix variables        ',i6,/,
     + ' maxnz  = maximum number of z-matrix cards            ',i6,/,
     + ' mxshel = maximum number of shells                    ',i6,/,
     + ' mxprim = maximum number of shell primitives          ',i6,/,
     + ' mxprms = maximum number of primitives in a shell     ',i6,/,
     + ' ',/,
     + ' following parameters refer to analysis modules ',/,
     + ' mxgaus = maximum number of orbital primitives        ',i6,/,
     + ' mxgrps = maximum number of shells                    ',i6,/,
     + ' ',/,
     + ' following parameters refer to graphics module ',/,
     + ' mxgrid, mxcalc,mxplot,mxrest ',4i6,/,
     + ' ',/,
     + ' following parameters refer to I/O system ',/,
     + ' parameters control no. of ed/mt files+ buffers ',/,
     + ' maxlfn *  no. of ed/mt streams         ',i6,/,
     + ' maxfrt *  no. of fortran data sets     ',i6,/,
     + ' maxbuf *  no. of fortran store buffers ',i6,/,
     + ' maxblo *  no. of blocks in 1 buffer    ',i6,/,
     + ' ',/,
     + ' following parameters refer to direct-CI module ',/,
     + ' nd200 = maximum number of orbitals in CI            ',i6,/,
     + ' parameters control canonical set size ',/,
     + ' mxcan1 *  default setting 2508 : high-spin 19606    ',i6,/,
     + ' mxcan2 *  default setting 5016 : high-spin 39212    ',i6,/,
     + ' ',/,
     + ' following parameters refer to full-CI module ',/,
     + ' lenci *  default setting 500000     ',i8,/,
     + ' *************************************************************',
     + /,1x)
c
      return
      end
      subroutine stong(csinp,cpinp,cdinp,nsto,sc,scc,
     +nucz,intyp,nangm,nbfs,minf,maxf,
     + stos1,stos2,stos3,stos4,stos5,stos6,stos7,
     +loc,ngauss,ns,nshmax,ngsmax,ierr1,ierr2,nat,iwr)
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
      common/blkin/eex(10),ccs(10),ccp(10),ccd(10)
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      dimension csinp(*),cpinp(*),cdinp(*)
      dimension scc(*)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      dimension stos1(*), stos2(*), stos3(*), stos4(*)
      dimension stos5(*), stos6(*), stos7(*)
c
      data pt75/0.75d0/
      data pt5,pi32,tol /0.5d+00,
     + 5.56832799683170d+00,1.0d-10/
c
      call vclr(eex,1,40)
c
       if(sc.le.0.0d0) sc=stos1(nucz)
       if(scc(1).le.0.0d0) scc(1)=stos2(nucz)
       if(scc(2).le.0.0d0) scc(2)=stos3(nucz)
       if(scc(3).le.0.0d0) scc(3)=stos4(nucz)
       if(scc(4).le.0.0d0) scc(4)=stos5(nucz)
       if(scc(5).le.0.0d0) scc(5)=stos6(nucz)
       if(scc(6).le.0.0d0) scc(6)=stos7(nucz)
      if (nucz .le. 2) then
c
c     ----- hydrogen and helium -----
c
      mxpass=1
c
      else if (nucz .le. 10) then
c
c     ----- lithium to neon -----
c
      mxpass=2
c
      else if (nucz .le. 18) then
c
c     ----- sodium to argon ------------
c
      mxpass=3
c
      else if (nucz .le. 20) then
c
c     ----- potassium and calcium------------
c
      mxpass=4
c
      else if (nucz .le. 36) then
c
c     ----- scandium to krypton
c
      mxpass=5
c
      else if (nucz .le. 38) then
c
c     ----- rubidium and strontium
c
      mxpass=6
c
      else if (nucz .le. 54) then
c
c     ----- yttrium to xenon
c
      mxpass=7
c
      else
      if (opg_root()) then
         write(iwr,*)'*** nuclear charge = ',nucz
      endif
      call caserr2('requested basis set not available')
      endif
c
c     ----- loop over shells -----
c
      igauss=nsto
      ipass = 0
  180 ipass = ipass+1
      if (nucz .gt. 2) go to 240
c
c     ----- h and he-----
c
      ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  240 if (nucz .gt. 10) go to 360
c
c     ----- li - ne ---
c
      go to (260,280),ipass
  260 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  280 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
c
  360 if (nucz .gt. 18) go to 440
c
c     ----- na - to - ar -----
c
      go to (380,400,420),ipass
  380 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  400 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
  420 ityp = 7
      scale=scc(2)
      call s3sp(eex,ccs,ccp,igauss)
      go to 700
c
  440 if (nucz .gt. 20) go to 500
c
c     ----- k - to - ca -----
c
      go to (450,460,470,480),ipass
  450 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  460 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
  470 ityp = 7
      scale=scc(2)
      call s3sp(eex,ccs,ccp,igauss)
      go to 700
  480 ityp = 11
      scale=scc(4)
      call s4sp(eex,ccs,ccp,igauss)
      go to 700
  500 if (nucz .gt. 36) go to 600
c
c     ----- sc to kr ----
c
      go to (510,520,530,540,550),ipass
  510 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  520 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
  530 ityp = 7
      scale=scc(2)
      call s3sp(eex,ccs,ccp,igauss)
      go to 700
  540 ityp = 8
      scale=scc(3)
      call s3d(eex,ccd,igauss)
      go to 700
  550 ityp = 11
      scale=scc(4)
      call s4sp(eex,ccs,ccp,igauss)
      go to 700
c
  600 if (nucz .gt. 38) go to 720
c
c     ----- rb to sr ----
c
      go to (610,620,630,640,650,660),ipass
  610 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  620 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
  630 ityp = 7
      scale=scc(2)
      call s3sp(eex,ccs,ccp,igauss)
      go to 700
  640 ityp = 8
      scale=scc(3)
      call s3d(eex,ccd,igauss)
      go to 700
  650 ityp = 11
      scale=scc(4)
      call s4sp(eex,ccs,ccp,igauss)
      go to 700
  660 ityp = 11
      scale=scc(6)
      call s5sp(eex,ccs,ccp,igauss)
      go to 700
c
c
c     ----- y to xe ----
c
  720 go to (810,820,830,840,850,860,870),ipass
  810 ityp = 1
      scale=sc
      call s1s(eex,ccs,igauss)
      go to 700
  820 ityp = 4
      scale=scc(1)
      call s2sp(eex,ccs,ccp,igauss)
      go to 700
  830 ityp = 7
      scale=scc(2)
      call s3sp(eex,ccs,ccp,igauss)
      go to 700
  840 ityp = 8
      scale=scc(3)
      call s3d(eex,ccd,igauss)
      go to 700
  850 ityp = 11
      scale=scc(4)
      call s4sp(eex,ccs,ccp,igauss)
      go to 700
  860 ityp = 12
      scale=scc(5)
      call s4d(eex,ccd,igauss)
      go to 700
  870 ityp = 11
      scale=scc(6)
      call s5sp(eex,ccs,ccp,igauss)
  700 continue
      nshell = nshell+1
      if (nshell .gt. nshmax) ierr1 = 1
      if (ierr1 .ne. 0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss .gt. ngsmax) ierr2 = 1
      if (ierr2 .ne. 0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      sc2 = scale*scale
      do 730 i = 1,igauss
      k = k1+i-1
      ex(k) = eex(i)*sc2
      csinp(k) = ccs(i)
      cpinp(k) = ccp(i)
      cdinp(k) = ccd(i)
      cs(k) = csinp(k)
      cp(k) = cpinp(k)
  730 cd(k) = cdinp(k)
c
c     ----- always unnormalize primitives... -----
c
      do 740 k = k1,k2
      ee = ex(k)+ex(k)
      facs = pi32/(ee*  dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      cs(k) = cs(k)/  dsqrt(facs)
      cp(k) = cp(k)/  dsqrt(facp)
  740 cd(k) = cd(k)/  dsqrt(facd)
c
      if (normf .eq. 1) go to 3034
c
      facs = 0.0d0
      facp = 0.0d0
      facd = 0.0d0
      do 3031 ig = k1,k2
      do 3031 jg = k1,ig
      ee = ex(ig)+ex(jg)
      fac = ee*dsqrt(ee)
      dums = cs(ig)*cs(jg)/fac
      dump = pt5*cp(ig)*cp(jg)/(ee*fac)
      dumd = pt75*cd(ig)*cd(jg)/(ee**2*fac)
      if (ig .eq. jg) go to 3032
      dums = dums+dums
      dump = dump+dump
      dumd = dumd+dumd
 3032 facs = facs+dums
      facp = facp+dump
      facd = facd+dumd
 3031  continue
      do 3033 ig = k1,k2
      if (facs .gt. tol) cs(ig) = cs(ig)/dsqrt(facs*pi32)
      if (facp .gt. tol) cp(ig) = cp(ig)/dsqrt(facp*pi32)
      if (facd .gt. tol) cd(ig) = cd(ig)/dsqrt(facd*pi32)
 3033 continue
 3034 if (ipass .lt. mxpass) go to 180
      return
      end
c
      subroutine s1s(exx,cs,ngauss)
c
c        *****  least squares fits to slater 1s functions       *****
c        *****  by gaussian 1s functions                        *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(1s---1(1s))                                 *****
c        *****  er1 4.3191e-02                                  *****
c
100   exx(1) = 2.709498091d-01
      cs(1) = 1.000000000d+00
      return
c
c        *****  sto(1s---2(1s))                                 *****
c        *****  er1 3.1606e-03                                  *****
c
 120  exx(1) = 8.518186635d-01
      cs(1) = 4.301284983d-01
      exx(2) = 1.516232927d-01
      cs(2) = 6.789135305d-01
      return
c
c        *****  sto(1s---3(1s))                                 *****
c        *****  er1 3.3053e-04                                  *****
c
 140  exx(1) = 2.227660584d+00
      cs(1) = 1.543289673d-01
      exx(2) = 4.057711562d-01
      cs(2) = 5.353281423d-01
      exx(3) = 1.098175104d-01
      cs(3) = 4.446345422d-01
      return
c
c        *****  sto(1s---4(1s))                                 *****
c        *****  er1 4.3763e-05                                  *****
c
 160  exx(1) = 5.216844534d+00
      cs(1) = 5.675242080d-02
      exx(2) = 9.546182760d-01
      cs(2) = 2.601413550d-01
      exx(3) = 2.652034102d-01
      cs(3) = 5.328461143d-01
      exx(4) = 8.801862774d-02
      cs(4) = 2.916254405d-01
      return
c
c        *****  sto(1s---5(1s))                                 *****
c        *****  er1 6.8840e-06                                  *****
c
 180  exx(1) = 1.130563696d+01
      cs(1) = 2.214055312d-02
      exx(2) = 2.071728178d+00
      cs(2) = 1.135411520d-01
      exx(3) = 5.786484833d-01
      cs(3) = 3.318161484d-01
      exx(4) = 1.975724573d-01
      cs(4) = 4.825700713d-01
      exx(5) = 7.445271746d-02
      cs(5) = 1.935721966d-01
      return
c
c        *****  sto(1s---6(1s))                                 *****
c        *****  er1 1.2372e-06                                  *****
c
 200  exx(1) = 2.310303149d+01
      cs(1) = 9.163596280d-03
      exx(2) = 4.235915534d+00
      cs(2) = 4.936149294d-02
      exx(3) = 1.185056519d+00
      cs(3) = 1.685383049d-01
      exx(4) = 4.070988982d-01
      cs(4) = 3.705627997d-01
      exx(5) = 1.580884151d-01
      cs(5) = 4.164915298d-01
      exx(6) = 6.510953954d-02
      cs(6) = 1.303340841d-01
c
      return
c
      end
c
      subroutine s2s(exx,cs,ngauss)
c
c        *****  least squares fits to slater 2s functions       *****
c        *****  by guassian 1s functions                        *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(2s---1(1s))                                 *****
c        *****  er1 6.2923e-03                                  *****
c
 100  exx(1) = 1.012151084d-01
      cs(1) = 1.000000000d+00
      return
c
c        *****  sto(2s---2(1s))                                 *****
c        *****  er1 2.6728e-03                                  *****
c
 120  exx(1) = 1.292278611d-01
      cs(1) = 7.470867124d-01
      exx(2) = 4.908584205d-02
      cs(2) = 2.855980556d-01
      return
c
c        *****  sto(2s---3(1s))                                 *****
c        *****  er1 6.8661e-05                                  *****
c
 140  exx(1) = 2.581578398d+00
      cs(1) = -5.994474934d-02
      exx(2) = 1.567622104d-01
      cs(2) = 5.960385398d-01
      exx(3) = 6.018332272d-02
      cs(3) = 4.581786291d-01
      return
c
c        *****  sto(2s---4(1s))                                 *****
c        *****  er1 2.6974e-05                                  *****
c
 160  exx(1) = 1.161525551d+01
      cs(1) = -1.198411747d-02
      exx(2) = 2.000243111d+00
      cs(2) = -5.472052539d-02
      exx(3) = 1.607280687d-01
      cs(3) = 5.805587176d-01
      exx(4) = 6.125744532d-02
      cs(4) = 4.770079976d-01
      return
c
c        *****  sto(2s---5(1s))                                 *****
c        *****  er1 2.7009e-06                                  *****
c
 180  exx(1) = 8.984956862d+00
      cs(1) = -1.596349096d-02
      exx(2) = 1.673710636d+00
      cs(2) = -5.685884883d-02
      exx(3) = 1.944726668d-01
      cs(3) = 3.698265599d-01
      exx(4) = 8.806345634d-02
      cs(4) = 5.480512593d-01
      exx(5) = 4.249068522d-02
      cs(5) = 1.472634893d-01
      return
c
c        *****  sto(2s---6(1s))                                 *****
c        *****  er1 3.2726e-07                                  *****
c
 200  exx(1) = 2.768496241d+01
      cs(1) = -4.151277819d-03
      exx(2) = 5.077140627d+00
      cs(2) = -2.067024148d-02
      exx(3) = 1.426786050d+00
      cs(3) = -5.150303337d-02
      exx(4) = 2.040335729d-01
      cs(4) = 3.346271174d-01
      exx(5) = 9.260298399d-02
      cs(5) = 5.621061301d-01
      exx(6) = 4.416183978d-02
      cs(6) = 1.712994697d-01
c
      return
c
      end
c
      subroutine s2p(exx,cp,ngauss)
c
c        *****  least squares fits to slater 2p functions       *****
c        *****  by gaussian 2p functions                        *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(2p---1(2p))                                 *****
c        *****  er1 4.8232e-02                                  *****
c
  100 exx(1) = 1.759666885d-01
      cp(1) = 1.000000000d+00
      return
c
c        *****  sto(2p---2(2p))                                 *****
c        *****  er1 3.0947e-03                                  *****
c
  120 exx(1) = 4.323908358d-01
      cp(1) = 4.522627513d-01
      exx(2) = 1.069439065d-01
      cp(2) = 6.713122642d-01
      return
c
c        *****  sto(2p---3(2p))                                 *****
c        *****  er1 2.6858e-04                                  *****
c
  140 exx(1) = 9.192379002d-01
      cp(1) = 1.623948553d-01
      exx(2) = 2.359194503d-01
      cp(2) = 5.661708862d-01
      exx(3) = 8.009805746d-02
      cp(3) = 4.223071752d-01
      return
c
c        *****  sto(2p---4(2p))                                 *****
c        *****  er1 2.9037e-05                                  *****
c
 160  exx(1) = 1.798260992d+00
      cp(1) = 5.713170255d-02
      exx(2) = 4.662622228d-01
      cp(2) = 2.857455515d-01
      exx(3) = 1.643718620d-01
      cp(3) = 5.517873105d-01
      exx(4) = 6.543927065d-02
      cp(4) = 2.632314924d-01
      return
c
c        *****  sto(2p---5(2p))                                 *****
c        *****  er1 3.7171e-06                                  *****
c
 180  exx(1) = 3.320386533d+00
      cp(1) = 2.079051117d-02
      exx(2) = 8.643257633d-01
      cp(2) = 1.235472099d-01
      exx(3) = 3.079819284d-01
      cp(3) = 3.667738886d-01
      exx(4) = 1.273309895d-01
      cp(4) = 4.834930290d-01
      exx(5) = 5.606243164d-02
      cp(5) = 1.653444074d-01
      return
c
c        *****  sto(2p---6(2p))                                 *****
c        *****  er1 5.4444e-07                                  *****
c
  200 exx(1) = 5.868285913d+00
      cp(1) = 7.924233646d-03
      exx(2) = 1.530329631d+00
      cp(2) = 5.144104825d-02
      exx(3) = 5.475665231d-01
      cp(3) = 1.898400060d-01
      exx(4) = 2.288932733d-01
      cp(4) = 4.049863191d-01
      exx(5) = 1.046655969d-01
      cp(5) = 4.012362861d-01
      exx(6) = 4.948220127d-02
      cp(6) = 1.051855189d-01
c
      return
c
      end
c
      subroutine s2sp(exx,cs,cp,ngauss)
c
c        *****  simultaneous least squares fits to slater 2s    *****
c        *****  and 2p functions by gaussian 1s + 2p functions  *****
c        *****  with common gaussian exponents                  *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
  100 call caserr2('requested basis set unavailable')
c        *****  sto(2s---2(1s),2p---2(2p))                      *****
c
c        *****  er1 5.5955e-03      er2 3.4793e-03              *****
c
  120 exx(1) = 3.842442531d-01
      cs(1) = 4.947176920d-02
      cp(1) = 5.115407076d-01
      exx(2) = 9.745448900d-02
      cs(2) = 9.637824081d-01
      cp(2) = 6.128198961d-01
      return
c
c        *****  sto(2s---3(1s),2p---3(2p))                      *****
c        *****  er1 6.4229e-04      er2 3.6036e-04              *****
c
 140  exx(1) = 9.942027296d-01
      cs(1) = -9.996722919d-02
      cp(1) = 1.559162750d-01
      exx(2) = 2.310313333d-01
      cs(2) = 3.995128261d-01
      cp(2) = 6.076837186d-01
      exx(3) = 7.513856000d-02
      cs(3) = 7.001154689d-01
      cp(3) = 3.919573931d-01
      return
c
c        *****  sto(2s---4(1s),2p---4(2p))                      *****
c        *****  er1 7.9474e-05      er2 4.8191e-05              *****
c
  160 exx(1) = 2.323503675d+00
      cs(1) = -6.220714565d-02
      cp(1) = 4.368434884d-02
      exx(2) = 5.029886906d-01
      cs(2) = 2.976804596d-05
      cp(2) = 2.863793984d-01
      exx(3) = 1.635406719d-01
      cs(3) = 5.588549221d-01
      cp(3) = 5.835753141d-01
      exx(4) = 6.281044213d-02
      cs(4) = 4.977673218d-01
      cp(4) = 2.463134378d-01
      return
c
c        *****  sto(2s---5(1s),2p---5(2p))                      *****
c        *****  er1 1.1708e-05      er2 7.2753e-06              *****
c
  180 exx(1) = 5.036294248d+00
      cs(1) = -2.940855088d-02
      cp(1) = 1.255609128d-02
      exx(2) = 1.032503477d+00
      cs(2) = -6.532746883d-02
      cp(2) = 1.075576962d-01
      exx(3) = 3.290598322d-01
      cs(3) = 1.289973181d-01
      cp(3) = 3.735975367d-01
      exx(4) = 1.279200125d-01
      cs(4) = 6.122899938d-01
      cp(4) = 5.102395637d-01
      exx(5) = 5.449486448d-02
      cs(5) = 3.461205655d-01
      cp(5) = 1.568281801d-01
      return
c
c        *****  sto(2s---6(1s),2p---6(2p))                      *****
c        *****  er1 2.0132e-06      er2 1.2173e-06              *****
c
  200 exx(1) = 1.030869372d+01
      cs(1) = -1.325278809d-02
      cp(1) = 3.759696623d-03
      exx(2) = 2.040359519d+00
      cs(2) = -4.699171014d-02
      cp(2) = 3.767936984d-02
      exx(3) = 6.341422177d-01
      cs(3) = -3.378537151d-02
      cp(3) = 1.738967435d-01
      exx(4) = 2.439773685d-01
      cs(4) = 2.502417861d-01
      cp(4) = 4.180364347d-01
      exx(5) = 1.059595374d-01
      cs(5) = 5.951172526d-01
      cp(5) = 4.258595477d-01
      exx(6) = 4.856900860d-02
      cs(6) = 2.407061763d-01
      cp(6) = 1.017082955d-01
      return
c
      end
c
      subroutine s3s(exx,cs,ngauss)
c
c        *****  least squares fits to slater 3s functions by    *****
c        *****  gaussian 1s functions                           *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(3s---1(1s))                                 *****
c        *****  er1 1.6764e-02                                  *****
c
 100  exx(1) = 5.296881757d-02
      cs(1) = 1.000000000d+00
      return
c
c        *****  sto(3s---2(1s))                                 *****
c        *****  er1 7.6424e-04                                  *****
c
 120  exx(1) = 6.694095822d-01
      cs(1) = -1.529645716d-01
      exx(2) = 5.837135094d-02
      cs(2) = 1.051370110d+00
      return
c
c        *****  sto(3s---3(1s))                                 *****
c        *****  er1 8.0718e-05                                  *****
c
  140 exx(1) = 5.641487709d-01
      cs(1) = -1.782577972d-01
      exx(2) = 6.924421391d-02
      cs(2) = 8.612761663d-01
      exx(3) = 3.269529097d-02
      cs(3) = 2.261841969d-01
      return
c
c        *****  sto(3s---4(1s))                                 *****
c        *****  er1 1.7254e-06                                  *****
c
 160  exx(1) = 1.513265591d+00
      cs(1) = -3.295496352d-02
      exx(2) = 4.262497508d-01
      cs(2) = -1.724516959d-01
      exx(3) = 7.643320863d-02
      cs(3) = 7.518511194d-01
      exx(4) = 3.760545063d-02
      cs(4) = 3.589627317d-01
      return
c
c        *****  sto(3s---5(1s))                                 *****
c        *****  er1 7.9816e-07                                  *****
c
  180 exx(1) = 4.275877914d+00
      cs(1) = -3.920358850d-03
      exx(2) = 1.132409433d+00
      cs(2) = -4.168430506d-02
      exx(3) = 4.016256968d-01
      cs(3) = -1.637440990d-01
      exx(4) = 7.732370620d-02
      cs(4) = 7.419373723d-01
      exx(5) = 3.800708627d-02
      cs(5) = 3.724364929d-01
      return
c
c        *****  sto(3s---6(1s))                                 *****
c        *****  er1 4.0662e-08                                  *****
c
 200  exx(1) = 3.273031938d+00
      cs(1) = -6.775596947d-03
      exx(2) = 9.200611311d-01
      cs(2) = -5.639325779d-02
      exx(3) = 3.593349765d-01
      cs(3) = -1.587856086d-01
      exx(4) = 8.636686991d-02
      cs(4) = 5.534527651d-01
      exx(5) = 4.797373812d-02
      cs(5) = 5.015351020d-01
      exx(6) = 2.724741144d-02
      cs(6) = 7.223633674d-02
      return
c
      end
c
      subroutine s3p(exx,cp,ngauss)
c
c        *****  least squares fits to slater 3p functions       *****
c        *****  with gaussian 2p functions                      *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(3p---1(2p))                                 *****
c        *****  er1 1.2745e-02                                  *****
c
 100  exx(1) = 9.113614253d-02
      cp(1) = 1.000000000d+00
      return
c
c        *****  sto(3p---2(2p))                                 *****
c        *****  er1 2.3771e-04                                  *****
c
 120  exx(1) = 1.458620964d-01
      cp(1) = 5.349653144d-01
      exx(2) = 5.664210742d-02
      cp(2) = 5.299607212d-01
      return
c
c        *****  sto(3p---3(2p))                                 *****
c        *****  er1 1.3487e-04                                  *****
c
 140  exx(1) = 2.692880368d+00
      cp(1) = -1.061945788d-02
      exx(2) = 1.489359592d-01
      cp(2) = 5.218564264d-01
      exx(3) = 5.739585040d-02
      cp(3) = 5.450015143d-01
      return
c
c        *****  sto(3p---4(2p))                                 *****
c        *****  er1 2.9785e-06                                  *****
c
  160 exx(1) = 1.853180239d+00
      cp(1) = -1.434249391d-02
      exx(2) = 1.915075719d-01
      cp(2) = 2.755177589d-01
      exx(3) = 8.655487938d-02
      cp(3) = 5.846750879d-01
      exx(4) = 4.184253862d-02
      cp(4) = 2.144986514d-01
      return
c
c        *****  sto(3p---5(2p))                                 *****
c        *****  er1 1.3387e-06                                  *****
c
  180 exx(1) = 6.466803859d+00
      cp(1) = -2.329023747d-03
      exx(2) = 1.555914802d+00
      cp(2) = -1.357395221d-02
      exx(3) = 1.955925255d-01
      cp(3) = 2.632185383d-01
      exx(4) = 8.809647701d-02
      cp(4) = 5.880427024d-01
      exx(5) = 4.234835707d-02
      cp(5) = 2.242794445d-01
      return
c
c        *****  sto(3p---6(2p))                                 *****
c        *****  er1 7.9285e-08                                  *****
c
  200 exx(1) = 5.077973607d+00
      cp(1) = -3.329929840d-03
      exx(2) = 1.340786940d+00
      cp(2) = -1.419488340d-02
      exx(3) = 2.248434849d-01
      cp(3) = 1.639395770d-01
      exx(4) = 1.131741848d-01
      cp(4) = 4.485358256d-01
      exx(5) = 6.076408893d-02
      cp(5) = 3.908813050d-01
      exx(6) = 3.315424265d-02
      cp(6) = 7.411456232d-02
      return
c
      end
c
      subroutine s3sp(exx,cs,cp,ngauss)
c
c        *****  simultaneous least squares fits to slater       *****
c        *****  3s and 3p functions with gaussian 1s and 2p     *****
c        *****  functions using common gaussian exponents       *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
c
c        *****  sto(3s---1(1s),3p---1(2p))                      *****
c        *****  er1 5.2955e-02      er2 4.1793e-02              *****
c
  100 exx(1) = 7.163507065d-02
      cs(1) = 1.000000000d+00
      cp(1) = 1.000000000d+00
      return
c
c        *****  sto(3s---2(1s),3p---2(2p))                      *****
c        *****  er1 7.3952e-03      er2 1.0928e-03              *****
c
  120 exx(1) = 1.939790354d-01
      cs(1) = -2.983986044d-01
      cp(1) = 3.480471912d-01
      exx(2) = 6.655227113d-02
      cs(2) = 1.227982887d+00
      cp(2) = 7.222523221d-01
      return
c
c        *****  sto(3s---3(1s),3p---3(2p))                      *****
c        *****  er1 3.1253e-04      er2 2.7910e-04              *****
c
  140 exx(1) = 4.828540806d-01
      cs(1) = -2.196203690d-01
      cp(1) = 1.058760429d-02
      exx(2) = 1.347150629d-01
      cs(2) = 2.255954336d-01
      cp(2) = 5.951670053d-01
      exx(3) = 5.272656258d-02
      cs(3) = 9.003984260d-01
      cp(3) = 4.620010120d-01
      return
c
c        *****  sto(3s---4(1s),3p---4(2p))                      *****
c        *****  er1 2.2403e-05      er2 2.8610e-05              *****
c
  160 exx(1) = 9.217072439d-01
      cs(1) = -8.529019644d-02
      cp(1) = -2.504945181d-02
      exx(2) = 2.534730123d-01
      cs(2) = -2.132074034d-01
      cp(2) = 1.686604461d-01
      exx(3) = 9.976476472d-02
      cs(3) = 5.920843928d-01
      cp(3) = 6.409553151d-01
      exx(4) = 4.439990833d-02
      cs(4) = 6.115584746d-01
      cp(4) = 2.779508957d-01
      return
c
c        *****  sto(3s---5(1s),3p---5(2p))                      *****
c        *****  er1 1.9857e-06      er2 3.5565e-06              *****
c
  180 exx(1) = 1.709113496d+00
      cs(1) = -2.662203391d-02
      cp(1) = -1.566883448d-02
      exx(2) = 4.654245585d-01
      cs(2) = -1.603484072d-01
      cp(2) = 7.214920506d-03
      exx(3) = 1.785129144d-01
      cs(3) = -4.779473307d-02
      cp(3) = 3.170854762d-01
      exx(4) = 8.067420191d-02
      cs(4) = 7.275158441d-01
      cp(4) = 5.818821382d-01
      exx(5) = 3.914669014d-02
      cs(5) = 4.123846408d-01
      cp(5) = 1.701799824d-01
      return
c
c        *****  sto(3s---6(1s),3p---6(2p))                      *****
c        *****  er1 2.2088e-07      er2 4.9043e-07              *****
c
  200 exx(1) = 3.080165240d+00
      cs(1) = -7.943126362d-03
      cp(1) = -7.139358907d-03
      exx(2) = 8.248959202d-01
      cs(2) = -7.100264172d-02
      cp(2) = -1.829277070d-02
      exx(3) = 3.093447349d-01
      cs(3) = -1.785026925d-01
      cp(3) = 7.621621428d-02
      exx(4) = 1.384683897d-01
      cs(4) = 1.510635058d-01
      cp(4) = 4.145098597d-01
      exx(5) = 6.852094951d-02
      cs(5) = 7.354914767d-01
      cp(5) = 4.889621471d-01
      exx(6) = 3.531333690d-02
      cs(6) = 2.760593123d-01
      cp(6) = 1.058816521d-01
      return
c
      end
c
      subroutine s3d(exx,cs,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(3d---1(3d))
c*********           er1 5.0728e-02
c
  100 exx(1) = 1.302270363d-01
      cs(1) = 1.000000000d+00
      return
c
c*****************      sto(3d---2(3d))
c*********           er1 2.9834e-03
c
  120 exx(1) = 2.777427345d-01
      cs(1) = 4.666137923d-01
      exx(2) = 8.336507714d-02
      cs(2) = 6.644706516d-01
      return
c
c*****************      sto(3d---3(3d))
c*********           er1 2.2687e-04
c
  140 exx(1) = 5.229112225d-01
      cs(1) = 1.686596060d-01
      exx(2) = 1.639595876d-01
      cs(2) = 5.847984817d-01
      exx(3) = 6.386630021d-02
      cs(3) = 4.056779523d-01
      return
c
c*****************      sto(3d---4(3d))
c*********           er1 2.1116e-05
  160 exx(1) = 9.185846715d-01
      cs(1) = 5.799057705d-02
      exx(2) = 2.920461109d-01
      cs(2) = 3.045581349d-01
      exx(3) = 1.187568890d-01
      cs(3) = 5.601358038d-01
      exx(4) = 5.286755896d-02
      cs(4) = 2.432423313d-01
      return
c
c*****************      sto(3d---5(3d))
c*********           er1 2.3115e-06
c
  180 exx(1) = 1.539033958d+00
      cs(1) = 2.020869128d-02
      exx(2) = 4.922090297d-01
      cs(2) = 1.321157923d-01
      exx(3) = 2.029756928d-01
      cs(3) = 3.911240346d-01
      exx(4) = 9.424112917d-02
      cs(4) = 4.779609701d-01
      exx(5) = 4.569058269d-02
      cs(5) = 1.463662294d-01
      return
c
c*****************      sto(3d---6(3d))
c*********           er1 2.8899e-07
c
  200 exx(1) = 2.488296923d+00
      cs(1) = 7.283828112d-03
      exx(2) = 7.981487853d-01
      cs(2) = 5.386799363d-02
      exx(3) = 3.311327490d-01
      cs(3) = 2.072139149d-01
      exx(4) = 1.559114463d-01
      cs(4) = 4.266269092d-01
      exx(5) = 7.877734732d-02
      cs(5) = 3.843100204d-01
      exx(6) = 4.058484363d-02
      cs(6) = 8.902827546d-02
      return
c
      end
c
      subroutine s4s(exx,cs,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(4s---1(1s))
c*********           er1 4.1760e-02
c
  100 exx(1) = 3.264600274d-02
      cs(1) = 1.000000000d+00
      return
c
c*****************      sto(4s---2(1s))
c*********           er1 1.7480e-04
c
  120 exx(1) = 2.441785453d-01
      cs(1) = -3.046656896d-01
      exx(2) = 4.051097664d-02
      cs(2) = 1.146877294d+00
      return
c
c*****************      sto(4s---3(1s))
c*********           er1 1.6867e-05
c
  140 exx(1) = 2.267938753d-01
      cs(1) = -3.349048323d-01
      exx(2) = 4.448178019d-02
      cs(2) = 1.056744667d+00
      exx(3) = 2.195294664d-02
      cs(3) = 1.256661680d-01
      return
c
c*****************      sto(4s---4(1s))
c*********           er1 8.4105e-07
c
  160 exx(1) = 3.242212833d-01
      cs(1) = -1.120682822d-01
      exx(2) = 1.663217177d-01
      cs(2) = -2.845426863d-01
      exx(3) = 5.081097451d-02
      cs(3) = 8.909873788d-01
      exx(4) = 2.829066600d-02
      cs(4) = 3.517811205d-01
      return
c
c*****************      sto(3d---6(3d))
c*********           er1 2.8899e-07
c
  180 exx(1) = 2.980263783d+00
      cs(1) = 1.513948997d-03
      exx(2) = 3.792228833d-01
      cs(2) = -7.316801548d-02
      exx(3) = 1.789717224d-01
      cs(3) = -3.143703799d-01
      exx(4) = 5.002110360d-02
      cs(4) = 9.032615169d-01
      exx(5) = 2.789361681d-02
      cs(5) = 3.294210848d-01
      return
c
c*****************      sto(4s---6(1s))
c*********           er1 5.9700e-09
c
  200 exx(1) = 3.232838646d+00
      cs(1) = 1.374817488d-03
      exx(2) = 3.605788802d-01
      cs(2) = -8.666390043d-02
      exx(3) = 1.717905487d-01
      cs(3) = -3.130627309d-01
      exx(4) = 5.277666487d-02
      cs(4) = 7.812787397d-01
      exx(5) = 3.163400284d-02
      cs(5) = 4.389247988d-01
      exx(6) = 1.874093091d-02
      cs(6) = 2.487178756d-02
      return
c
      end
c
      subroutine s4p(exx,cp,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(4p---1(2p))
c*********           er1 3.9516e-03
c
  100 exx(1) = 5.578350235d-02
      cp(1) = 1.000000000d+00
      return
c
c*****************      sto(4p---2(2p))
c*********           er1 2.1579e-03
c
  120 exx(1) = 6.190052680d-02
      cp(1) = 8.743116767d-01
      exx(2) = 2.648418407d-02
      cp(2) = 1.513640107d-01
      return
c
c*****************      sto(4p---3(2p))
c*********           er1 1.2621e-05
c
  140 exx(1) = 4.859692220d-01
      cp(1) = -6.147823411d-02
      exx(2) = 7.430216918d-02
      cp(2) = 6.604172234d-01
      exx(3) = 3.653340923d-02
      cp(3) = 3.932639495d-01
      return
c
c*****************      sto(4p---4(2p))
c*********           er1 5.4740e-06
c
  160 exx(1) = 1.492607880d+00
      cp(1) = -6.035216774d-03
      exx(2) = 4.327619272d-01
      cp(2) = -6.013310874d-02
      exx(3) = 7.553156064d-02
      cp(3) = 6.451518200d-01
      exx(4) = 3.706272183d-02
      cp(4) = 4.117923820d-01
      return
c
c*****************      sto(4p---5(2p))
c*********           er1 1.3595e-07
c
  180 exx(1) = 1.091977298d+00
      cp(1) = -1.143929558d-02
      exx(2) = 3.719985051d-01
      cp(2) = -6.322651538d-02
      exx(3) = 8.590019352d-02
      cp(3) = 4.398907721d-01
      exx(4) = 4.786503860d-02
      cp(4) = 5.245859166d-01
      exx(5) = 2.730479990d-02
      cp(5) = 1.017072253d-01
      return
c
c*****************      sto(4p---6(2p))
c*********           er1 1.3897e-08
c
  200 exx(1) = 2.389722618d+00
      cp(1) = -1.665913575d-03
      exx(2) = 7.960947826d-01
      cp(2) = -1.657464971d-02
      exx(3) = 3.415541380d-01
      cp(3) = -5.958513378d-02
      exx(4) = 8.847434525d-02
      cp(4) = 4.053115554d-01
      exx(5) = 4.958248334d-02
      cp(5) = 5.433958189d-01
      exx(6) = 2.816929784d-02
      cp(6) = 1.204970491d-01
      return
c
      end
c
      subroutine s4sp(exx,cs,cp,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
  100 call caserr2('requested basis set unavailable')
c
c*****************      sto(4s---2(1s),4p---2(2p))
c*********           er1 5.1660e-03      er2 4.1650e-03
c
  120 exx(1) = 1.242349951d-01
      cs(1) = -6.615603018d-01
      cp(1) = 1.309890515d-01
      exx(2) = 4.999746691d-02
      cs(2) = 1.510754134d+00
      cp(2) = 8.946431268d-01
      return
c
c*****************      sto(4s---3(1s),4p---3(2p))
c*********           er1 1.6710e-04      er2 2.9851e-04
c
  140 exx(1) = 2.464581400d-01
      cs(1) = -3.088441215d-01
      cp(1) = -1.215468600d-01
      exx(2) = 9.095855374d-02
      cs(2) = 1.960641166d-02
      cp(2) = 5.715227604d-01
      exx(3) = 4.016825636d-02
      cs(3) = 1.131034442d+00
      cp(3) = 5.498949471d-01
      return
c
c*****************      sto(4s---4(1s),4p---4(2p))
c*********           er1 1.5073e-05      er2 1.6717e-05
c
  160 exx(1) = 4.389527791d-01
      cs(1) = -5.893559813d-02
      cp(1) = -7.180522657d-02
      exx(2) = 1.624449952d-01
      cs(2) = -4.289230261d-01
      cp(2) = 3.726863425d-02
      exx(3) = 6.983976374d-02
      cs(3) = 5.472998468d-01
      cp(3) = 6.965842249d-01
      exx(4) = 3.469204610d-02
      cs(4) = 7.853692030d-01
      cp(4) = 3.310236850d-01
      return
c
c*****************      sto(4s---5(1s),4p---5(2p))
c*********           er1 1.6981e-06      er2 1.1140e-06
c
  180 exx(1) = 7.857764224d-01
      cs(1) = -1.628249797d-03
      cp(1) = -2.399519595d-02
      exx(2) = 2.685947984d-01
      cs(2) = -2.138931336d-01
      cp(2) = -8.118321551d-02
      exx(3) = 1.155680976d-01
      cs(3) = -3.039583456d-01
      cp(3) = 2.238628944d-01
      exx(4) = 5.839870334d-02
      cs(4) = 8.491175450d-01
      cp(4) = 6.726900408d-01
      exx(5) = 3.082651925d-02
      cs(5) = 5.164147260d-01
      cp(5) = 2.003802746d-01
      return
c
c*****************      sto(3s---6(1s),3p---6(2p))
c    ------no er1 and er2 values given (?)   renate------
c
  200 exx(1) = 3.080165240d+00
      exx(2) = 8.248959202d-01
      exx(3) = 3.093447349d-01
      exx(4) = 1.384683897d-01
      exx(5) = 6.852094951d-02
      exx(6) = 3.531333690d-02
      cs(1) = -7.943126362d-03
      cs(2) = -7.100264172d-02
      cs(3) = -1.785026925d-01
      cs(4) = 1.510635058d-01
      cs(5) = 7.354914767d-01
      cs(6) = 2.760593123d-01
      cp(1) = -7.139358907d-03
      cp(2) = -1.829277070d-02
      cp(3) = 7.621621428d-02
      cp(4) = 4.145098597d-01
      cp(5) = 4.889621471d-01
      cp(6) = 1.058816521d-01
      return
c
      end
c
      subroutine s4d(exx,cd,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cd(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(4d---1(3d))
c*********           er1 2.0360e-02
c
  100 exx(1) = 7.941656339d-02
      cd(1) = 1.000000000d+00
      return
c
c*****************      sto(4d---2(3d))
c*********           er1 3.7488e-04
c
  120 exx(1) = 1.330958892d-01
      cd(1) = 4.932764167d-01
      exx(2) = 5.272119659d-02
      cd(2) = 5.918727866d-01
      return
c
c*****************      sto(4d---3(3d))
c*********           er1 1.2729e-05
c
  140 exx(1) = 1.777717219d-01
      cd(1) = 2.308552718d-01
      exx(2) = 8.040647350d-02
      cd(2) = 6.042409177d-01
      exx(3) = 3.949855551d-02
      cd(3) = 2.595768926d-01
      return
c
c*****************      sto(4d---4(3d))
c*********           er1 5.7655e-06
c
 160  exx(1) = 1.995825422d+00
      cd(1) = -2.816702620d-03
      exx(2) = 1.823461280d-01
      cd(2) = 2.177095871d-01
      exx(3) = 8.197240896d-02
      cd(3) = 6.058047348d-01
      exx(4) = 4.000634951d-02
      cd(4) = 2.717811257d-01
      return
c
c*****************      sto(4d---5(3d))
c*********           er1 1.5122e-07
c
  180 exx(1) = 1.522122079d+00
      cd(1) = -3.673711876d-03
      exx(2) = 2.173041823d-01
      cd(2) = 1.167122499d-01
      exx(3) = 1.084876577d-01
      cd(3) = 4.216476416d-01
      exx(4) = 5.836797641d-02
      cd(4) = 4.547673415d-01
      exx(5) = 3.206682246d-02
      cd(5) = 1.037803318d-01
      return
c
c*****************      sto(4d---6(3d))
c*********           er1 8.0888e-08
c
  200 exx(1) = 4.634239420d+00
      cd(1) = -4.749842876d-04
      exx(2) = 1.341648295d+00
      cd(2) = -3.566777891d-03
      exx(3) = 2.209593028d-01
      cd(3) = 1.108670481d-01
      exx(4) = 1.101467943d-01
      cd(4) = 4.159646930d-01
      exx(5) = 5.904190370d-02
      cd(5) = 4.621672517d-01
      exx(6) = 3.232628887d-02
      cd(6) = 1.081250196d-01
      return
c
      end
c
      subroutine s5s(exx,cs,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(5s---1(1s))
c*********           er1 7.1317e-02
c
  100 exx(1) = 2.216912938d-02
      cs(1) = 1.000000000d+00
      return
c
c*****************      sto(5s---2(1s))
c*********           er1 1.4289e-04
c
  120 exx(1) = 1.213425654d-01
      cs(1) = -5.114756049d-01
      exx(2) = 3.133152144d-02
      cs(2) = 1.307377277d+00
      return
c
c*****************      sto(5s---3(1s))
c*********           er1 4.1124e-05
c
  140 exx(1) = 1.080198458d-01
      cs(1) = -6.617401158d-01
      exx(2) = 4.408119382d-02
      cs(2) = 7.467595004d-01
      exx(3) = 2.610811810d-02
      cs(3) = 7.146490945d-01
      return
c
c*****************      sto(5s---4(1s))
c*********           er1 5.4159e-07
c
  160 exx(1) = 8.602284252d-01
      cs(1) = 1.103657561d-02
      exx(2) = 1.189050200d-01
      cs(2) = -5.606519023d-01
      exx(3) = 3.446076176d-02
      cs(3) = 1.179429987d+00
      exx(4) = 1.974798796d-02
      cs(4) = 1.734974376d-01
      return
c
c*****************      sto(5s---5(1s))
c*********           er1 7.0816e-08
c
  180 exx(1) = 7.403763257d-01
      cs(1) = 1.375523371d-02
      exx(2) = 1.367990863d-01
      cs(2) = -3.097344179d-01
      exx(3) = 9.135301779d-02
      cs(3) = -3.199192259d-01
      exx(4) = 3.726907315d-02
      cs(4) = 1.084547038d+00
      exx(5) = 2.241490836d-02
      cs(5) = 3.345288361d-01
      return
c
c*****************      sto(5s---6(1s))
c*********           er1 7.9988e-09
c
  200 exx(1) = 1.410128298d+00
      cs(1) = 2.695439582d-03
      exx(2) = 5.077878915d-01
      cs(2) = 1.850157487d-02
      exx(3) = 1.847926858d-01
      cs(3) = -9.588628125d-02
      exx(4) = 1.061070594d-01
      cs(4) = -5.200673560d-01
      exx(5) = 3.669584901d-02
      cs(5) = 1.087619490d+00
      exx(6) = 2.213558430d-02
      cs(6) = 3.103964343d-01
      return
c
      end
c
      subroutine s5p(exx,cp,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cp(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(5p---1(2p))
c*********           er1 8.0406e-03
c
  100 exx(1) = 3.769845216d-02
      cp(1) = 1.000000000d+00
      return
c
c*****************      sto(5p---2(2p))
c*********           er1 6.8912e-04
c
  120 exx(1) = 2.691294191d-01
      cp(1) = -1.034227010d-01
      exx(2) = 3.980805011d-02
      cp(2) = 1.033376378d+00
      return
c
c*****************      sto(5p---3(2p))
c*********           er1 6.5766e-06
c
  140 exx(1) = 2.127482317d-01
      cp(1) = -1.389529695d-01
      exx(2) = 4.729648620d-02
      cp(2) = 8.076691064d-01
      exx(3) = 2.604865324d-02
      cp(3) = 2.726029342d-01
      return
c
c*****************      sto(5p---4(2p))
c*********           er1 4.4493e-07
c
  160 exx(1) = 3.962838833d-01
      cp(1) = -1.801459207d-02
      exx(2) = 1.838858552d-01
      cp(2) = -1.360777372d-01
      exx(3) = 4.943555157d-02
      cp(3) = 7.533973719d-01
      exx(4) = 2.750222273d-02
      cp(4) = 3.409304859d-01
      return
c
c*****************      sto(5p---5(2p))
c
c*********           er1 1.3483e-08
c
  180 exx(1) = 3.422168934d-01
c
      cp(1) = -3.113958289d-02
      exx(2) = 1.665099900d-01
      cp(2) = -1.374007017d-01
      exx(3) = 5.443732013d-02
      cp(3) = 5.573881018d-01
      exx(4) = 3.367775277d-02
      cp(4) = 4.855428100d-01
      exx(5) = 2.091949042d-02
      cp(5) = 6.605423564d-02
      return
c
c*****************      sto(5p---6(2p))
c*********           er1 3.9301e-09
c
  200 exx(1) = 3.778623374d+00
      cp(1) = 1.163246387d-04
      exx(2) = 3.499121109d-01
      cp(2) = -2.920771322d-02
      exx(3) = 1.683175469d-01
      cp(3) = -1.381051233d-01
      exx(4) = 5.404070736d-02
      cp(4) = 5.706134877d-01
      exx(5) = 3.328911801d-02
      cp(5) = 4.768808140d-01
      exx(6) = 2.063815019d-02
      cp(6) = 6.021665516d-02
      return
c
      end
      subroutine s5sp(exx,cs,cp,ngauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cs(*),cp(*)
c
c     from pietro, et.al. inorg. chem. 20(11), 3650 (1981)
      go to (1,2,3,4,5,6),ngauss
c     *****************************************************************
    1 call caserr2('requested basis set unavailable')
c     *****************************************************************
    2 continue
c*****************      sto(5s---2(1s),5p---2(2p))
c*********           er1 5.1660e-03      er2 4.1650d-03
      exx(1)= 8.124323517d-02
      cs(1) =-1.093689792d+00
      cp(1) =-1.047403535d-01
      exx(2)= 3.930355482d-02
      cs(2) = 1.881374656d+00
      cp(2) = 1.087626793d+00
      return
c     *****************************************************************
    3 continue
c*****************      sto(5s---3(1s),5p---3(2p))
c*********           er1 1.6710e-04      er2 2.9851d-04
      exx(1)= 1.349013828d-01
      cs(1) =-3.842642607d-01
      cp(1) =-3.481691526d-01
      exx(2)= 7.263605443d-02
      cs(2) =-1.972567438d-01
      cp(2) = 6.290323690d-01
      exx(3)= 3.208462257d-02
      cs(3) = 1.375495512d+00
      cp(3) = 6.662832743d-01
      return
c     *****************************************************************
    4 continue
c*****************      sto(5s---4(1s),5p---4(2p))
      exx(1)= 2.577707269d-01
      cs(1) = 4.045111204d-02
      cp(1) =-8.586149652d-02
      exx(2)= 1.189639508d-01
      cs(2) =-6.576691017d-01
      cp(2) =-1.090154130d-01
      exx(3)= 5.270776581d-02
      cs(3) = 3.792524265d-01
      cp(3) = 7.234042095d-01
      exx(4)= 2.870356027d-02
      cs(4) = 1.038589734d+00
      cp(4) = 4.117428940d-01
      return
c     *****************************************************************
    5 continue
c*****************      sto(5s---5(1s),5p---5(2p))
      exx(1)= 4.822618266d-01
      cs(1) = 3.037662332d-02
      cp(1) =-1.257672570d-02
      exx(2)= 1.813516259d-01
      cs(2) =-1.825060513d-01
      cp(2) =-1.563935424d-01
      exx(3)= 8.588912730d-02
      cs(3) =-6.591492319d-01
      cp(3) = 6.500746946d-02
      exx(4)= 4.617413271d-02
      cs(4) = 9.186932397d-01
      cp(4) = 7.859820497d-01
      exx(5)= 2.582552665d-02
      cs(5) = 6.952739594d-01
      cp(5) = 2.582283808d-01
      return
c     *****************************************************************
    6 continue
c******************    sto(5s---6(1s),5p---6(2p))
      exx(1)= 7.701420258d-01
      cs(1) = 1.267447151d-02
      cp(1) =-1.105673292d-03
      exx(2)= 2.756268915d-01
      cs(2) = 3.266734789d-03
      cp(2) =-6.243132446d-02
      exx(3)= 1.301847480d-01
      cs(3) =-4.307553999d-01
      cp(3) =-1.628476766d-01
      exx(4)= 6.953441940d-02
      cs(4) =-3.231998963d-01
      cp(4) = 3.210328714d-01
      exx(5)= 4.002545502d-02
      cs(5) = 1.104322879d+00
      cp(5) = 6.964579592d-01
      exx(6)= 2.348388309d-02
      cs(6) = 4.368498703d-01
      cp(6) = 1.493146125d-01
      return
      end
      subroutine s5d(exx,cd,ngauss)
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension exx(*),cd(*)
      go to (100,120,140,160,180,200),ngauss
c
c*****************      sto(5d---1(3d))
c*********           er1 6.7308e-03
c
  100 exx(1) = 5.352200793d-02
      cd(1) = 1.000000000d+00
      return
c
c*****************      sto(5d---2(3d))
c*********           er1 3.1257e-04
c
  120 exx(1) = 6.906014388d-02
      cd(1) = 6.539405185d-01
      exx(2) = 3.399457777d-02
      cd(2) = 3.948945302d-01
      return
c
c*****************      sto(5d---3(3d))
c*********           er1 2.4664e-05
c
  140 exx(1) = 4.913352950d-01
      cd(1) = -2.010175008d-02
      exx(2) = 7.329090601d-02
      cd(2) = 5.899370608d-01
      exx(3) = 3.594209290d-02
      cd(3) = 4.658445960d-01
      return
c
c*****************      sto(5d---4(3d))
c*********           er1 1.2737e-06
c
  160 exx(1) = 4.230617826d-01
      cd(1) = -2.421626009d-02
      exx(2) = 8.293863702d-02
      cd(2) = 3.937644956d-01
      exx(3) = 4.590326388d-02
      cd(3) = 5.489520286d-01
      exx(4) = 2.628744797d-02
      cd(4) = 1.190436963d-01
      return
c
c*****************      sto(5d---5(3d))
c*********           er1 8.6929e-08
c
  180 exx(1) = 9.702946470d-01
      cd(1) = -3.231527611d-03
      exx(2) = 3.603270196d-01
      cd(2) = -2.434931372d-02
      exx(3) = 8.668717752d-02
      cd(3) = 3.440817054d-01
      exx(4) = 4.833708379d-02
      cd(4) = 5.693674376d-01
      exx(5) = 2.751899341d-02
      cd(5) = 1.511340183d-01
      return
c
c*****************      sto(5d---6(3d))
c*********           er1 1.4086e-08
c
  200 exx(1) = 8.820520428d-01
      cd(1) = -4.097377019d-03
      exx(2) = 3.410838409d-01
      cd(2) = -2.508271857d-02
      exx(3) = 9.204308840d-02
      cd(3) = 2.648458555d-01
      exx(4) = 5.472831774d-02
      cd(4) = 5.097437054d-01
      exx(5) = 3.391202830d-02
      cd(5) = 2.654483467d-01
      exx(6) = 2.108227374d-02
      cd(6) = 2.623132212d-02
      return
c
      end
c
      subroutine n31g(csinp,cpinp,cdinp,cfinp,opol,oryd,
     +  sc,scc,rleesf,igauss,nucz,intyp,nangm,
     +  nbfs,minf,maxf,loc,ngauss,ns,nshmax,ngsmax,ierr1,
     +  ierr2,nat,iwr)
c
c     ----- set up the n-31g basis set -----
c
c  j.chem.phys. 
c    54,  724- 728 (1971)    -31g for h
c                           4-31g for c,n,o,f
c    56, 2257-2261 (1972) 5,6-31g for c,n,o,f
c    56, 4233-4234 (1972)   4-31g for b
c    62, 2921-2923 (1975) 4,6-31g for li
c                           6-31g for b
c    66,  879- 880 (1977) 4,6-31g for be
c    66,  879- 880 (1977)   6-31g for be
c    56, 5255-5257 (1972)   4-31g for p,s,cl
c    77, 3654-3665 (1982)   6-31g for na,mg,al,p-ar
c  chem.phys.lett. 
c    76,  163- 168 (1980)   6-31g for si
c  -31g for he and 4,6-31g for ne are apparently unpublished.
c  these basis sets were taken directly from gaussian 82.
c    jcp 109 (4) 1223-1229, 1998.  6-31g* for k, ca, sc-zn
c    V. Rassolov, J.A. Pople, M. Ratner and T.L. Windus,
c
c     4-31g available for h-he, li-ne, p-cl
c           note that 4-31g for li and be is actually 5-21g
c     5-31g available for h-he, c-f.
c     6-31g available for h-he, li-ne, na-ar, k-zn
c
c
c  Both polarisation and diffuse functions may be added
c  to these and the N-21G sets (see dn21g, and n21g). 
c  The following polarisation functions are applied to the
c  Pople Basis sets:
c
c                POPLE    POPN311   DUNNING   HUZINAGA    HONDO7
c                ------   -------   -------   --------    ------
c           H    1.1(p)    0.75(p)   1.0(p)     1.0(p)    1.0(p)
c           He   1.1(p)    0.75(p)   1.0(p)     1.0(p)    1.0(p)
c
c           Li   0.2       0.200                0.076(p)
c           Be   0.4       0.255                0.164(p)  0.32
c           B    0.6       0.401     0.70       0.388     0.50
c           C    0.8       0.626     0.75       0.600     0.72
c           N    0.8       0.913     0.80       0.864     0.98
c           O    0.8       1.292     0.85       1.154     1.28
c           F    0.8       1.750     0.90       1.496     1.62
c           Ne   0.8       2.304     1.00       1.888     2.00
c
c           Na   0.175                          0.061(p)  0.157
c           Mg   0.175                          0.101(p)  0.234
c           Al   0.325                          0.198     0.311
c           Si   0.395                          0.262     0.388
c           P    0.55                           0.340     0.465
c           S    0.65                           0.421     0.542
c           Cl   0.75                           0.514     0.619
c           Ar   0.85                           0.617     0.696
c
c           K    0.2 (d)                        0.039(p)
c           Ca   0.2 (d)                        0.059(p)
c        Sc-Zn   0.8 (f)
c           Ga   0.207               0.141
c           Ge   0.246               0.202
c           As   0.293               0.273
c           Se   0.338               0.315
c           Br   0.389               0.338
c           Kr   0.443               0.318
c
c           Rb   0.11                           0.034(p)
c           Sr   0.11                           0.048(p)
c
c  The following  diffuse functions are applied:
c
c             H
c          0.0360
c            Li      Be       B       C       N       O       F
c          0.0074  0.0207  0.0315  0.0438  0.0639  0.0845  0.1076
c            Na      Mg      Al      Si       P       S      Cl
c          0.0076  0.0146  0.0318  0.0331  0.0348  0.0405  0.0483
c                            Ga      Ge      As      Se      Br
c                          0.0205  0.0222  0.0287  0.0318  0.0376
c                            In      Sn      Sb      Te       I
c                          0.0223  0.0231  0.0259  0.0306  0.0368
c                            Tl      Pb      Bi      Po      At
c                          0.0170  0.0171  0.0215  0.0230  0.0294
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
      common/blkin/e123(48),cs123(48),cp23(48),cd23(48),cf23(48)
      common/junk /ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      dimension csinp(*),cpinp(*),cdinp(*),cfinp(*)
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
      dimension  e1(6), e2(6), e3(6), e4(6), e5(6), e6(6), e7(6), e8(6)
      dimension cs1(6),cs2(6),cs3(6),cs4(6),cs5(6),cs6(6),cs7(6),cs8(6)
      dimension cp1(6),cp2(6),cp3(6),cp4(6),cp5(6),cp6(6),cp7(6),cp8(6)
      dimension cd1(6),cd2(6),cd3(6),cd4(6),cd5(6),cd6(6),cd7(6),cd8(6)
      dimension cf1(6),cf2(6),cf3(6),cf4(6),cf5(6),cf6(6),cf7(6),cf8(6)
      dimension scalf(8),scc(*),rleesf(36,2)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      equivalence (e123(1),e1(1)), (cs123(1),cs1(1)), (cp23(1),cp1(1)),
     +            (e123(7),e2(1)), (cs123(7),cs2(1)), (cp23(7),cp2(1)),
     +           (e123(13),e3(1)),(cs123(13),cs3(1)),(cp23(13),cp3(1)),
     +           (e123(19),e4(1)),(cs123(19),cs4(1)),(cp23(19),cp4(1)),
     +           (e123(25),e5(1)),(cs123(25),cs5(1)),(cp23(25),cp5(1)),
     +           (e123(31),e6(1)),(cs123(31),cs6(1)),(cp23(31),cp6(1)),
     +           (e123(37),e7(1)),(cs123(37),cs7(1)),(cp23(37),cp7(1)),
     +           (e123(43),e8(1)),(cs123(43),cs8(1)),(cp23(43),cp8(1))
      equivalence (cd23(1) ,cd1(1)), (cf23(1),  cf1(1)),
     +            (cd23(7) ,cd2(1)), (cf23(7),  cf2(1)),
     +           (cd23(13) ,cd3(1)), (cf23(13), cf3(1)),
     +           (cd23(19) ,cd4(1)), (cf23(19), cf4(1)),
     +           (cd23(25) ,cd5(1)), (cf23(25), cf5(1)),
     +           (cd23(31) ,cd6(1)), (cf23(31), cf6(1)),
     +           (cd23(37) ,cd7(1)), (cf23(37), cf7(1)),
     +           (cd23(43) ,cd8(1)), (cf23(43), cf8(1))
      data pt75/0.75d0/
      data pt1875/1.875d+00/
      data done,pt5,pi32,tol /1.0d+00,0.5d+00,
     + 5.56832799683170d+00,1.0d-10/
c
c     ----- get exponents, contraction coefficients, and
c          scale factors. -----
c
c
      call vclr(e123,1,240)
      igdum = igauss
      mxpass=7
      if(opol)mxpass=8
      if (nucz .le. 2) then
c
c     ----- hydrogen,helium -----
c
      call ezero(e1,e2,cs1,cp2,nucz)
      scalf(5) = sc
      scalf(6) = scc(1)
      scalf(7) = scc(2)
      scalf(8) = scc(3)
      if (scalf(5) .eq. 0.0d0) scalf(5) = rleesf(nucz,1)
      if (scalf(6) .eq. 0.0d0) scalf(6) = rleesf(nucz,2)
      if (scalf(7) .eq. 0.0d0) scalf(7) = done
c
      else if (nucz .le. 10)  then
c
c     ----- lithium to neon -----
c
      call eone(e1,e2,e3,cs1,cs2,cp2,cd3,igauss,nucz)
      scalf(4) = sc
      scalf(5) = scc(1)
      scalf(6) = scc(2)
      scalf(7) = scc(3)
      scalf(8) = scc(4)
      if (scalf(4) .eq. 0.0d0) scalf(4) = done
      if (scalf(5) .eq. 0.0d0) scalf(5) = rleesf(nucz,1)
      if (scalf(6) .eq. 0.0d0) scalf(6) = rleesf(nucz,2)
      if (scalf(7) .eq. 0.0d0) scalf(7) = done
c
      else if (nucz .le. 18) then
c
c     ----- sodium to argon -----
c
      call etwo(e1,e2,e3,e4,cs1,cs2,cs3,cp2,cp3,cd4,igauss,nucz,
     +          iwr)
      scalf(3) = sc
      scalf(4) = scc(1)
      scalf(5) = scc(2)
      scalf(6) = scc(3)
      scalf(7) = scc(4)
      scalf(8) = scc(5)
      if (scalf(3) .eq. 0.0d0) scalf(3) = done
      if (scalf(4) .eq. 0.0d0) scalf(4) = done
      if (scalf(5) .eq. 0.0d0) scalf(5) = rleesf(nucz,1)
      if (scalf(6) .eq. 0.0d0) scalf(6) = rleesf(nucz,2)
      if (scalf(7) .eq. 0.0d0) scalf(7) = done
c
      else if (nucz .le. 30) then
c
c     ----- potassium to zinc  -----
c
      if (oryd) then
        if (opg_root()) then
           write(iwr,*)'*** nuclear charge = ',nucz
        endif
        call caserr2('requested basis set not available')
      endif
      call ethree(e1,e2,e3,e4,e5,e6,e7,e8,
     +            cs1,cs2,cs3,cs4,cs5,cp2,cp3,
     +            cp4,cp5,cd6,cd7,cf8,igauss,nucz)
       if (nucz .le. 20) then
c
c     ----- potassium, calcium -----
c
        scalf(3) = sc
        scalf(4) = scc(1)
        scalf(5) = scc(2)
        scalf(6) = scc(3)
        scalf(7) = scc(4)
        scalf(8) = scc(5)
        if (scalf(3) .eq. 0.0d0) scalf(3) = done
        if (scalf(4) .eq. 0.0d0) scalf(4) = done
        if (scalf(5) .eq. 0.0d0) scalf(5) = done
        if (scalf(6) .eq. 0.0d0) scalf(6) = rleesf(nucz,1)
        if (scalf(7) .eq. 0.0d0) scalf(7) = rleesf(nucz,2)
       else
c
c     ----- scandium to zinc ----
c
        scalf(1) = sc
        scalf(2) = scc(1)
        scalf(3) = scc(2)
        scalf(4) = scc(3)
        scalf(5) = scc(4)
        scalf(6) = scc(5)
        scalf(7) = scc(6)
        scalf(8) = scc(7)
        if (scalf(1) .eq. 0.0d0) scalf(1) = done
        if (scalf(2) .eq. 0.0d0) scalf(2) = done
        if (scalf(3) .eq. 0.0d0) scalf(3) = done
        if (scalf(4) .eq. 0.0d0) scalf(4) = rleesf(nucz,1)
        if (scalf(5) .eq. 0.0d0) scalf(5) = rleesf(nucz,2)
        if (scalf(6) .eq. 0.0d0) scalf(6) = done
        if (scalf(7) .eq. 0.0d0) scalf(7) = done
       endif
      else
        if (opg_root()) then
           write(iwr,*)'*** nuclear charge = ',nucz
        endif
        call caserr2('requested basis set not available')
      endif
      if (scalf(8) .eq. 0.0d0) scalf(8) = done
c
c     ----- loop over 3,4,5 or 7 shells. ----- (n-31+g)
c     ----- loop over 4,5,6,or 8 shells. ----- (n-31+g*)
c
      ipass = 0
  180 ipass = ipass+1
      if (nucz .gt. 2) go to 240
c
c     ----- hydrogen, helium -----
c
      go to (180,180,180,180,200,220,230,235),ipass
  200 ityp = 1
      igauss = 3
      ig = 0
      go to 420
  220 ityp = 1
      igauss = 1
      ig = 3
      go to 420
  230 ityp = 1
      igauss = 1
      ig = 4
      if(oryd) then
       go to 420
      else
       go to 640
      endif
  235 ityp = 17
      igauss = 1
      ig = 6
      go to 420
  240 if (nucz .gt. 4) go to 620
c
c     ----- lithium to beryllium -----
c     note that 4-31g for li and be is actually a 5-21g set
c
      go to (180,180,180,560,580,600,610,615),ipass
  560 ityp = 1
      ig = 0
      if(igdum.eq.4) igauss=5
      go to 420
  580 ityp = 4
      igauss = 3
      if(igdum.eq.4.or.igdum.eq.5)  igauss=2
      ig = 6
      go to 420
  600 ityp = 4
      igauss = 1
      ig = 9
      go to 420
  610 ityp = 22
      igauss = 1
      ig = 10
      if(oryd) then
       go to 420
      else
       go to 640
      endif
  615 ityp = 18
      igauss = 1
      ig = 12
      go to 420
c
  620 if (nucz .gt. 10) go to 320
c
c     ----- boron to fluorine -----
c
      go to (180,180,180,260,280,300,310,315),ipass
  260 ityp = 1
      ig = 0
      go to 420
  280 ityp = 4
      igauss = 3
      ig = 6
      go to 420
  300 ityp = 4
      igauss = 1
      ig = 9
      go to 420
  310 ityp = 22
      igauss = 1
      ig = 10
      if(oryd) then
       go to 420
      else
       go to 640
      endif
  315 ityp = 18
      igauss = 1
      ig = 12
      go to 420
  320 if (nucz .gt. 18) go to 710
c
c     ----- phosphorus to chlorine -----
c
      go to (180,180,340,360,380,400,410,415),ipass
  340 ityp = 1
      ig = 0
      go to 420
  360 ityp = 4
      ig = 6
      go to 420
  380 ityp = 7
      igauss = 3
      ig = 12
      go to 420
  400 ityp = 7
      igauss = 1
      ig = 15
      go to 420
  410 ityp = 22
      igauss = 1
      ig = 16
      if(oryd) then
       go to 420
      else
       go to 640
      endif
  415 ityp = 18
      igauss = 1
      ig = 18
      go to 420
c
  710 if (nucz .gt. 20) go to 810
c
c     ----- potassium, calcium -----
c
      go to (180,180,730,740,750,760,770,780),ipass
  730 ityp = 1
      ig = 0
      go to 420
  740 ityp = 4
      ig = 6
      go to 420
  750 ityp = 7
      ig = 12
      go to 420
  760 ityp = 11
      igauss = 3
      ig = 18
      go to 420
  770 ityp = 11
      igauss = 1
      ig = 24
      go to 420
  780 ityp = 18
      igauss = 1
      ig = 30
      go to 420
c
  810 if (nucz .gt. 30) then
        if (opg_root()) then
           write(iwr,*)'*** nuclear charge = ',nucz
        endif
        call caserr2('requested basis set not available')
      endif
c
c     ----- scandium to zinc
c
      go to (830,840,850,860,870,880,890,900),ipass
  830 ityp = 1
      ig = 0
      go to 420
  840 ityp = 4
      ig = 6
      go to 420
  850 ityp = 7
      ig = 12
      go to 420
  860 ityp = 11
      igauss = 3
      ig = 18
      go to 420
  870 ityp = 11
      igauss = 1
      ig = 24
      go to 420
  880 ityp = 18
      igauss = 3
      ig = 30
      go to 420
  890 ityp = 18
      igauss = 1
      ig = 36
      go to 420
  900 ityp = 19
      igauss = 1
      ig = 42
c
  420 continue
      nshell = nshell+1
      if (nshell .gt. nshmax) ierr1 = 1
      if (ierr1 .ne. 0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss .gt. ngsmax) ierr2 = 1
      if (ierr2 .ne. 0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do 440 i = 1,igauss
      k = k1+i-1
      ex(k) = e123(ig+i)*scalf(ipass)**2
      csinp(k) = cs123(ig+i)
      cpinp(k) = cp23(ig+i)
      cdinp(k) = cd23(ig+i)
      cfinp(k) = cf23(ig+i)
      cs(k) = csinp(k)
      cp(k) = cpinp(k)
      cd(k) = cdinp(k)
  440 cf(k) = cfinp(k)
c
c     ----- always unnormalize primitives... -----
c
      do 460 k = k1,k2
      ee = ex(k)+ex(k)
      facs = pi32/(ee*dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      facf = pt1875*facs/(ee**3)
      cs(k) = cs(k)/dsqrt(facs)
      cp(k) = cp(k)/dsqrt(facp)
      if(ityp.eq.18)cd(k) = cd(k)/dsqrt(facd)
      if(ityp.eq.19)cf(k) = cf(k)/dsqrt(facf)
  460 continue
c
c     ----- if(normf.eq.0) normalize basis functions. -----
c
      if (normf .eq. 1) go to 540
      facs = 0.0d0
      facp = 0.0d0
      facd = 0.0d0
      facf = 0.0d0
      do 500 ig = k1,k2
      do 500 jg = k1,ig
      ee = ex(ig)+ex(jg)
      fac = ee*dsqrt(ee)
      dums = cs(ig)*cs(jg)/fac
      dump = pt5*cp(ig)*cp(jg)/(ee*fac)
      dumd = pt75*cd(ig)*cd(jg)/(ee**2*fac)
      dumf = pt1875*cf(ig)*cf(jg)/(ee*ee*ee*fac)
      if (ig .eq. jg) go to 480
      dums = dums+dums
      dump = dump+dump
      dumd = dumd+dumd
      dumf = dumf+dumf
  480 facs = facs+dums
      facp = facp+dump
      facd = facd+dumd
      facf = facf+dumf
  500 continue
c
      do 520 ig = k1,k2
      if (facs .gt. tol) cs(ig) = cs(ig)/dsqrt(facs*pi32)
      if (facp .gt. tol) cp(ig) = cp(ig)/dsqrt(facp*pi32)
      if (facd .gt. tol) cd(ig) = cd(ig)/dsqrt(facd*pi32)
      if (facf .gt. tol) cf(ig) = cf(ig)/dsqrt(facf*pi32)
  520 continue
  540 continue
c
  640 if (ipass .lt. mxpass) go to 180
      return
c
      end
c
      subroutine ezero(e1,e2,cs1,cp2,ia)
c
c    *****  energy optimized 31g basis for hydrogen & helium ****
c    *****  compatible with 4-31g 5-31g and 6-31g series    *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e1(*),e2(*),cs1(*),cp2(*)
c
      go to (100,200),ia
c
c        *****  hydrogen 
c  
 100  e1(1) = 1.300773400d+01
      cs1(1) = 3.349460434d-02
      e1(2) = 1.962079420d+00
      cs1(2) = 2.347269535d-01
      e1(3) = 4.445289530d-01
      cs1(3) = 8.137573262d-01
      e1(4) = 1.219491560d-01
      cs1(4) = 1.000000000d+00
      e1(5) = 0.0360d+00
      cs1(5) = 1.000000000d+00
      e2(1) = 1.10d+00
      cp2(1) = 1.000000000d+00
      return
c
c        ***** helium
c
 200  e1(1) = 38.421634d+00
      cs1(1) = 0.023766d+00
      e1(2) = 5.77803d+00
      cs1(2) = 0.154679d+00
      e1(3) = 1.241774d+00
      cs1(3) = 0.469630d+00
      e1(4) = 0.297964d+00
      cs1(4) = 1.0d+00
      e1(5) = 0.0360d+00
      cs1(5) = 1.000000000d+00
      e2(1) = 1.10d+00
      cp2(1) = 1.000000000d+00
      return
c
      end
c
      subroutine eone(e1,e2,e3,cs1,cs2,cp2,cd3,ngauss,ia)
c
c        *****  energy optimized 4-31g 5-31g and 6-31g bases    *****
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e1(*),e2(*),cs1(*),cs2(*),cp2(*)
      dimension e3(*),cd3(*)
c
      go to (800,800,810,910, 100,240,380,520,660,670),ia
c
c        *****  lithium
c
  810 e3(1) = 0.2d+00
      cd3(1)= 1.0d+00
      go to (815,815,815,816,816,820),ngauss
  815 call caserr2('requested basis set for lithium not available')
      return
  816 continue
c     there is no 4-31g basis set.  4-31g molecular calculations
c     employ the 5-21g basis for lithium and beryllium.
c     5-21g
      e1(1)=.275394444d03
      cs1(1)=.612184691d-02
      e1(2)=.414351754d02
      cs1(2)=.451129615d-01
      e1(3)=.936699378d01
      cs1(3)=.192694150d00
      e1(4)=.253772533d01
      cs1(4)=.468544208d00
      e1(5)=.746636540d00
      cs1(5)=.440607515d00
      e2(1)=.692397267d00
      cs2(1)=-.252536797d00
      cp2(1)=.143591732d00
      e2(2)=.821924442d-01
      cs2(2)=.109734080d01
      cp2(2)=.947803050d00
      e2(4)=.322374501d-01
      go to 817
  820 continue
c        -----  6-31g                                           -----
      e1(1)=6.42418915d+02
      cs1(1)=2.14260781d-03
      e1(2)=9.67985153d+01
      cs1(2)=1.62088715d-02
      e1(3)=2.20911212d+01
      cs1(3)=7.73155725d-02
      e1(4)=6.20107025d+00
      cs1(4)=2.45786052d-01
      e1(5)=1.93511768d+00
      cs1(5)=4.70189004d-01
      e1(6)=6.36735789d-01
      cs1(6)=3.45470845d-01
      e2(1)=2.19145858d+00
      cs2(1)=-3.50917459d-02
      cp2(1)=8.94150804d-03
      e2(2)=5.96126266d-01
      cs2(2)=-1.91232844d-01
      cp2(2)=1.41009464d-01
      e2(3)=7.45154442d-02
      cs2(3)=1.08398780d+00
      cp2(3)=9.45363695d-01
      e2(4)=2.86686637d-02
  817 cs2(4)=1.00000000d+00
      cp2(4)=1.00000000d+00
      e2(5)=0.0074d+00
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  beryllium
c
  910 e3(1) = 0.4d+00
      cd3(1)= 1.0d+00
      go to (915,915,915,916,916,920),ngauss
  915 call caserr
     *('requested basis set for beryllium not available')
      return
  916 continue
c
c     4-31g  // use 5-21g basis
c
      e1(1)=554.010d0
      e1(2)=83.2631d0
      e1(3)=18.8635d0
      e1(4)=5.17782d0
      e1(5)=1.55602d0
      cs1(1)=.00540997d0
      cs1(2)=.0402515d0
      cs1(3)=.176858d0
      cs1(4)=.452559d0
      cs1(5)=.470293d0
      e2(1)=1.35899d0
      e2(2)=.284533d0
      cs2(1)=-.477429d0
      cs2(2)=1.24745d0
      cp2(1)=.201142d0
      cp2(2)=.884483d0
      e2(4)=.0804858d0
      go to 917
  920 continue
c        -----  6-31g                                           -----
      e1(1)=1.26706795d+03
      cs1(1)=1.93997980d-03
      e1(2)=1.90356049d+02
      cs1(2)=1.47864428d-02
      e1(3)=4.32959013d+01
      cs1(3)=7.17950984d-02
      e1(4)=1.21440986d+01
      cs1(4)=2.36347608d-01
      e1(5)=3.80922922d+00
      cs1(5)=4.71763178d-01
      e1(6)=1.26846854d+00
      cs1(6)=3.55183381d-01
      e2(1)=5.69387692d+00
      cs2(1)=-2.88764407d-02
      cp2(1)=4.83496858d-03
      e2(2)=1.55562520d+00
      cs2(2)=-1.77564983d-01
      cp2(2)=1.44044900d-01
      e2(3)=1.71854950d-01
      cs2(3)=1.07163052d+00
      cp2(3)=9.49691875d-01
      e2(4)=5.71804680d-02
  917 cs2(4)=1.00000000d+00
      cp2(4)=1.00000000d+00
      e2(5)= 0.0207d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  boron
c
  100 e3(1) = 0.6d+00
      cd3(1)= 1.0d+00
      go to (120,120,120,180,120,121),ngauss
  120 call caserr2('requested basis set for boron not available')
c
c        *****  4-31g                                           *****
c
  180 e1(1) = 3.307528520d+02
      cs1(1) = 1.799417960d-02
      e1(2) = 4.984386500d+01
      cs1(2) = 1.246937000d-01
      e1(3) = 1.111705350d+01
      cs1(3) = 4.343353750d-01
      e1(4) = 2.922724310d+00
      cs1(4) = 5.609793740d-01
      e2(1) = 5.355136790d+00
      cs2(1) = -1.303870780d-01
      cp2(1) = 6.374292250d-02
      e2(2) = 1.370915820d+00
      cs2(2) = -2.514343900d-01
      cp2(2) = 2.761330530d-01
      e2(3) = 4.037878930d-01
      cs2(3) = 1.205129200d+00
      cp2(3) = 7.773865960d-01
      e2(4) = 1.149706420d-01
      go to 122
  121 continue
c        -----  6-31g                                           -----
      e1(1)=2.06888d+03
      cs1(1)=1.86627d-03
      e1(2)=3.10650d+02
      cs1(2)=1.42515d-02
      e1(3)=7.06830d+01
      cs1(3)=6.95516d-02
      e1(4)=1.98611d+01
      cs1(4)=2.32573d-01
      e1(5)=6.29930d+00
      cs1(5)=4.67079d-01
      e1(6)=2.12703d+00
      cs1(6)=3.63431d-01
      e2(1)=4.45657d+00
      cs2(1)=-1.30394d-01
      cp2(1)=7.45976d-02
      e2(2)=1.12201d+00
      cs2(2)=-1.30789d-01
      cp2(2)=3.07847d-01
      e2(3)=3.38780d-01
      cs2(3)=1.13094d+00
      cp2(3)=7.43457d-01
      e2(4)=1.01045d-01
  122 cs2(4)=1.00000d+00
      cp2(4)=1.00000d+00
      e2(5)= 0.0315d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  carbon
c
  240 e3(1) = 0.8d+00
      cd3(1)= 1.0d+00
      go to (260,260,260,320,340,360),ngauss
  260 call caserr2('requested basis set for carbon not available')
c
c        *****  4-31g
c
  320 e1(1) = 4.869669280d+02
      cs1(1) = 1.772582290d-02
      e1(2) = 7.337109420d+01
      cs1(2) = 1.234778670d-01
      e1(3) = 1.641345790d+01
      cs1(3) = 4.338754000d-01
      e1(4) = 4.344983560d+00
      cs1(4) = 5.615041970d-01
      e2(1) = 8.673525310d+00
      cs2(1) = -1.213837490d-01
      cp2(1) = 6.354538410d-02
      e2(2) = 2.096619260d+00
      cs2(2) = -2.273384980d-01
      cp2(2) = 2.982677570d-01
      e2(3) = 6.046513290d-01
      cs2(3) = 1.185173920d+00
      cp2(3) = 7.621032280d-01
      e2(4) = 1.697095320d-01
      go to 265
c
c        *****  5-31g                                           *****
c
  340 e1(1) = 1.264250200d+03
      cs1(1) = 5.473495530d-03
      e1(2) = 1.901442980d+02
      cs1(2) = 4.079115360d-02
      e1(3) = 4.312858670d+01
      cs1(3) = 1.812203490d-01
      e1(4) = 1.194438200d+01
      cs1(4) = 4.634824780d-01
      e1(5) = 3.651484710d+00
      cs1(5) = 4.524711990d-01
      e2(1) = 7.942730520d+00
      cs2(1) = -1.207731490d-01
      cp2(1) = 6.867749980d-02
      e2(2) = 1.907237930d+00
      cs2(2) = -1.697932090d-01
      cp2(2) = 3.141028910d-01
      e2(3) = 5.535773850d-01
      cs2(3) = 1.149811710d+00
      cp2(3) = 7.459685180d-01
      e2(4) = 1.585119750d-01
      go to 265
c
c        *****  6-31g                                           *****
c
  360 e1(1) = 3.047524880d+03
      cs1(1) = 1.834737130d-03
      e1(2) = 4.573695180d+02
      cs1(2) = 1.403732280d-02
      e1(3) = 1.039486850d+02
      cs1(3) = 6.884262220d-02
      e1(4) = 2.921015530d+01
      cs1(4) = 2.321844430d-01
      e1(5) = 9.286662960d+00
      cs1(5) = 4.679413480d-01
      e1(6) = 3.163926960d+00
      cs1(6) = 3.623119850d-01
      e2(1) = 7.868272350d+00
      cs2(1) = -1.193324200d-01
      cp2(1) = 6.899906660d-02
      e2(2) = 1.881288540d+00
      cs2(2) = -1.608541520d-01
      cp2(2) = 3.164239610d-01
      e2(3) = 5.442492580d-01
      cs2(3) = 1.143456440d+00
      cp2(3) = 7.443082910d-01
      e2(4) = 1.559860190d-01
  265 cs2(4) = 1.000000000d+00
      cp2(4) = 1.000000000d+00
      e2(5)= 0.0438d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  nitrogen
c
  380 e3(1) = 0.8d+00
      cd3(1)= 1.0d+00
      go to (400,400,400,460,480,500),ngauss
  400 call caserr2('requested basis set for nitrogen not available')
c
c        *****  4-31g
c
  460 e1(1) = 6.712795030d+02
      cs1(1) = 1.759825110d-02
      e1(2) = 1.012016620d+02
      cs1(2) = 1.228462410d-01
      e1(3) = 2.269996590d+01
      cs1(3) = 4.337821410d-01
      e1(4) = 6.040609000d+00
      cs1(4) = 5.614182170d-01
      e2(1) = 1.264524000d+01
      cs2(1) = -1.174892990d-01
      cp2(1) = 6.402034430d-02
      e2(2) = 2.981719040d+00
      cs2(2) = -2.139940160d-01
      cp2(2) = 3.112025550d-01
      e2(3) = 8.494317690d-01
      cs2(3) = 1.174502110d+00
      cp2(3) = 7.527482390d-01
      e2(4) = 2.352813130d-01
      go to 490
c
c        *****  5-31g                                           *****
c
  480 e1(1) = 1.745284990d+03
      cs1(1) = 5.421222340d-03
      e1(2) = 2.625413270d+02
      cs1(2) = 4.043458220d-02
      e1(3) = 5.958251390d+01
      cs1(3) = 1.804448930d-01
      e1(4) = 1.654434660d+01
      cs1(4) = 4.634396240d-01
      e1(5) = 5.084876110d+00
      cs1(5) = 4.526306780d-01
      e2(1) = 1.194496010d+01
      cs2(1) = -1.161469430d-01
      cp2(1) = 6.732112530d-02
      e2(2) = 2.800090140d+00
      cs2(2) = -1.757313530d-01
      cp2(2) = 3.221237760d-01
      e2(3) = 7.981580700d-01
      cs2(3) = 1.150529060d+00
      cp2(3) = 7.421561640d-01
      e2(4) = 2.235978380d-01
      go to 490
c
c        *****  6-31g
c
  500 e1(1) = 4.173511460d+03
      cs1(1) = 1.834772160d-03
      e1(2) = 6.274579110d+02
      cs1(2) = 1.399462700d-02
      e1(3) = 1.429020930d+02
      cs1(3) = 6.858655180d-02
      e1(4) = 4.023432930d+01
      cs1(4) = 2.322408730d-01
      e1(5) = 1.282021290d+01
      cs1(5) = 4.690699480d-01
      e1(6) = 4.390437010d+00
      cs1(6) = 3.604551990d-01
      e2(1) = 1.186242410d+01
      cs2(1) = -1.149611820d-01
      cp2(1) = 6.757974390d-02
      e2(2) = 2.771431290d+00
      cs2(2) = -1.691174790d-01
      cp2(2) = 3.239072960d-01
      e2(3) = 7.878975580d-01
      cs2(3) = 1.145851950d+00
      cp2(3) = 7.408951400d-01
      e2(4) = 2.207741540d-01
  490 cs2(4) = 1.000000000d+00
      cp2(4) = 1.000000000d+00
      e2(5)= 0.0639d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  oxygen
c
  520 e3(1) = 0.8d+00
      cd3(1)= 1.0d+00
      go to (540,540,540,600,620,640),ngauss
 540  call caserr2('requested basis set for oxygen not available')
c
c        *****  4-31g                                           *****
c
 600  e1(1) = 8.832728600d+02
      cs1(1) = 1.755062800d-02
      e1(2) = 1.331292800d+02
      cs1(2) = 1.228292230d-01
      e1(3) = 2.990640790d+01
      cs1(3) = 4.348835840d-01
      e1(4) = 7.978677160d+00
      cs1(4) = 5.600108040d-01
      e2(1) = 1.652325950d+01
      cs2(1) = -1.134010030d-01
      cp2(1) = 6.854527470d-02
      e2(2) = 3.856837080d+00
      cs2(2) = -1.772864660d-01
      cp2(2) = 3.312254350d-01
      e2(3) = 1.092728880d+00
      cs2(3) = 1.150407930d+00
      cp2(3) = 7.346078780d-01
      e2(4) = 2.955850070d-01
      go to 630
c
c        *****  5-31g                                           *****
c
  620 e1(1) = 2.296705350d+03
      cs1(1) = 5.402590830d-03
      e1(2) = 3.454369370d+02
      cs1(2) = 4.033911540d-02
      e1(3) = 7.840108180d+01
      cs1(3) = 1.805909170d-01
      e1(4) = 2.181041640d+01
      cs1(4) = 4.643773710d-01
      e1(5) = 6.721723480d+00
      cs1(5) = 4.511585830d-01
      e2(1) = 1.591091840d+01
      cs2(1) = -1.118074560d-01
      cp2(1) = 7.079811160d-02
      e2(2) = 3.695812810d+00
      cs2(2) = -1.519256420d-01
      cp2(2) = 3.386958530d-01
      e2(3) = 1.043638440d+00
      cs2(3) = 1.133714340d+00
      cp2(3) = 7.277198890d-01
      e2(4) = 2.838735860d-01
      go to 630
c
c        *****  6-31g                                           *****
c
  640 e1(1) = 5.484671660d+03
      cs1(1) = 1.831074430d-03
      e1(2) = 8.252349460d+02
      cs1(2) = 1.395017220d-02
      e1(3) = 1.880469580d+02
      cs1(3) = 6.844507810d-02
      e1(4) = 5.296450000d+01
      cs1(4) = 2.327143360d-01
      e1(5) = 1.689757040d+01
      cs1(5) = 4.701928980d-01
      e1(6) = 5.799635340d+00
      cs1(6) = 3.585208530d-01
      e2(1) = 1.585513340d+01
      cs2(1) = -1.107775490d-01
      cp2(1) = 7.087426820d-02
      e2(2) = 3.673026820d+00
      cs2(2) = -1.480262620d-01
      cp2(2) = 3.397528390d-01
      e2(3) = 1.034345220d+00
      cs2(3) = 1.130767010d+00
      cp2(3) = 7.271585770d-01
      e2(4) = 2.811389240d-01
  630 cs2(4) = 1.000000000d+00
      cp2(4) = 1.000000000d+00
      e2(5)= 0.0845d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
      return
c
c        *****  fluorine
c
  660 e3(1) = 0.8d+00
      cd3(1)= 1.0d+00
      go to (680,680,680,740,760,780),ngauss
  680 call caserr2('requested basis set for fluorine not available')
c
c     4-31g
c
  740 e1(1) = 1.126162690d+03
      cs1(1) = 1.747576090d-02
      e1(2) = 1.697431570d+02
      cs1(2) = 1.225230890d-01
      e1(3) = 3.818151120d+01
      cs1(3) = 4.349985020d-01
      e1(4) = 1.021203590d+01
      cs1(4) = 5.598121670d-01
      e2(1) = 2.149536670d+01
      cs2(1) = -1.110570790d-01
      cp2(1) = 6.988875080d-02
      e2(2) = 4.989777570d+00
      cs2(2) = -1.683221010d-01
      cp2(2) = 3.393875100d-01
      e2(3) = 1.403573860d+00
      cs2(3) = 1.143625550d+00
      cp2(3) = 7.279589810d-01
      e2(4) = 3.730318350d-01
      go to 770
c
c        *****  5-31g                                           *****
c
  760 e1(1) = 2.927122920d+03
      cs1(1) = 5.380259070d-03
      e1(2) = 4.402360610d+02
      cs1(2) = 4.020133150d-02
      e1(3) = 9.993134980d+01
      cs1(3) = 1.804364600d-01
      e1(4) = 2.783870100d+01
      cs1(4) = 4.647456520d-01
      e1(5) = 8.602577920d+00
      cs1(5) = 4.506029000d-01
      e2(1) = 2.090380370d+01
      cs2(1) = -1.094195750d-01
      cp2(1) = 7.156794460d-02
      e2(2) = 4.831641890d+00
      cs2(2) = -1.494113590d-01
      cp2(2) = 3.450879300d-01
      e2(3) = 1.353751870d+00
      cs2(3) = 1.130961860d+00
      cp2(3) = 7.228825230d-01
      e2(4) = 3.610197530d-01
      go to 770
c
c     6-31g
c
  780 e1(1) = 7.001713090d+03
      cs1(1) = 1.819616900d-03
      e1(2) = 1.051366090d+03
      cs1(2) = 1.391607960d-02
      e1(3) = 2.392856900d+02
      cs1(3) = 6.840532450d-02
      e1(4) = 6.739744530d+01
      cs1(4) = 2.331857600d-01
      e1(5) = 2.151995730d+01
      cs1(5) = 4.712674390d-01
      e1(6) = 7.403101300d+00
      cs1(6) = 3.566185460d-01
      e2(1) = 2.084795280d+01
      cs2(1) = -1.085069750d-01
      cp2(1) = 7.162872430d-02
      e2(2) = 4.808308340d+00
      cs2(2) = -1.464516580d-01
      cp2(2) = 3.459121030d-01
      e2(3) = 1.344069860d+00
      cs2(3) = 1.128688580d+00
      cp2(3) = 7.224699570d-01
      e2(4) = 3.581513930d-01
  770 cs2(4) = 1.000000000d+00
      cp2(4) = 1.000000000d+00
      e2(5)=0.1076d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
c
      return
c
c        *****  neon
c
  670 e3(1) = 0.8d+00
      cd3(1)= 1.0d+00
      go to (675,675,675,880,675,890),ngauss
  675 call caserr2('requested basis set for neon not available')
  880 continue
c     ----- 4-31g -----
      e1(1)=0.139793208d+04
      cs1(1)=0.174238054d-01
      e1(2)=0.210769781d+03
      cs1(2)=0.122272745d+00
      e1(3)=0.474672569d+02
      cs1(3)=0.435014232d+00
      e1(4)=0.127226263d+02
      cs1(4)=0.559714642d+00
      e2(1)=0.272130332d+02
      cs2(1)=-0.109609439d+00
      cp2(1)=0.704403067d-01
      e2(2)=0.629413435d+01
      cs2(2)=-0.164124890d+00
      cp2(2)=0.343993047d+00
      e2(3)=0.176005125d+01
      cs2(3)=0.114015159d+01
      cp2(3)=0.724514960d+00
      e2(4)=0.461866992d+00
      go to 685
  890 continue
c     ----- 6-31g -----
      e1(1)=8425.85153d+00
      cs1(1)=0.00188434805d+00
      e1(2)=1268.5194d+00
      cs1(2)=0.0143368994d+00
      e1(3)=289.621414d+00
      cs1(3)=0.0701096233d+00
      e1(4)=81.8590040d+00
      cs1(4)=0.237373266d+00
      e1(5)=26.2515079d+00
      cs1(5)=0.473007126d+00
      e1 (6)=9.09472051d+00
      cs1(6)=0.348401241d+00
      e2(1)=26.532131d+00
      cs2(1)=-0.107118287d+00
      cp2(1)=0.0719095885d+00
      e2(2)=6.10175501d+00
      cs2(2)=-0.146163821d+00
      cp2(2)=0.349513372d+00
      e2(3)=1.69627153d+00
      cs2(3)=1.12777350d+00
      cp2(3)=0.719940512d+00
      e2(4)=0.445818700d+00
  685 cs2(4)=1.0d+00
      cp2(4)=1.0d+00
      e2(5)= 0.11d0
      cs2(5)=1.00000000d+00
      cp2(5)=1.00000000d+00
 800  return
c
      end
c
      subroutine etwo(e1,e2,e3,e4,cs1,cs2,cs3,cp2,cp3,cd4,
     *     ngauss,ia,iwr)
c
c        *****  energy optimized 4-31g basis for second
c        *****  row atoms
c        ***** polarisation functions from
c        ***** francl,pietro,hehre,binkley,gordon,defrees and pople
c        ***** j. Chem. Phys 77 (1982) 3654
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e1(*),e2(*),e3(*),e4(*)
      dimension cs1(*),cs2(*),cs3(*),cp2(*),cp3(*)
      dimension cd4(*)
c
      go to (560,560,560,560,560,560,560,560,560,560,10 , 50,100,120,
     +     140,280,420,550),ia
c
c        *****  sodium
c
   10 e4(1) = 0.175d0
      cd4(1)= 1.000d0
      go to (15,15,15,15,15,20),ngauss
   15 call caserr2('requested basis set for sodium not available')
   20 continue
c        -----  66-31g                                          -----
      e1(1)=9.99320d+03
      cs1(1)=1.93766d-03
      e1(2)=1.49989d+03
      cs1(2)=1.48070d-02
      e1(3)=3.41951d+02
      cs1(3)=7.27055d-02
      e1(4)=9.46796d+01
      cs1(4)=2.52629d-01
      e1(5)=2.97345d+01
      cs1(5)=4.93242d-01
      e1(6)=1.00063d+01
      cs1(6)=3.13169d-01
      e2(1)=1.50963d+02
      cs2(1)=-3.54208d-03
      cp2(1)=5.00166d-03
      e2(2)=3.55878d+01
      cs2(2)=-4.39588d-02
      cp2(2)=3.55109d-02
      e2(3)=1.11683d+01
      cs2(3)=-1.09752d-01
      cp2(3)=1.42825d-01
      e2(4)=3.90201d+00
      cs2(4)=1.87398d-01
      cp2(4)=3.38620d-01
      e2(5)=1.38177d+00
      cs2(5)=6.46699d-01
      cp2(5)=4.51579d-01
      e2(6)=4.66382d-01
      cs2(6)=3.06058d-01
      cp2(6)=2.73271d-01
      e3(1)=4.97966d-01
      cs3(1)=-2.48503d-01
      cp3(1)=-2.30225d-02
      e3(2)=8.43529d-02
      cs3(2)=-1.31704d-01
      cp3(2)=9.50359d-01
      e3(3)=6.66350d-02
      cs3(3)=1.23352d+00
      cp3(3)=5.98579d-02
      e3(4)=2.59544d-02
      cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)= 0.0076d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  magnesium                                       *****
c
   50 e4(1) = 0.175d0
      cd4(1)= 1.000d0
      go to (55,55,55,55,55,60),ngauss
   55 call caserr2('requested basis set for magnesium not available')
   60 continue
c        -----  66-31g                                          -----
      e1(1)=1.17228d+04
      cs1(1)=1.97783d-03
      e1(2)=1.75993d+03
      cs1(2)=1.51140d-02
      e1(3)=4.00846d+02
      cs1(3)=7.39108d-02
      e1(4)=1.12807d+02
      cs1(4)=2.49191d-01
      e1(5)=3.59997d+01
      cs1(5)=4.87928d-01
      e1(6)=1.21828d+01
      cs1(6)=3.19662d-01
      e2(1)=1.89180d+02
      cs2(1)=-3.23717d-03
      cp2(1)=4.92813d-03
      e2(2)=4.52119d+01
      cs2(2)=-4.10079d-02
      cp2(2)=3.49888d-02
      e2(3)=1.43563d+01
      cs2(3)=-1.12600d-01
      cp2(3)=1.40725d-01
      e2(4)=5.13886d+00
      cs2(4)=1.48633d-01
      cp2(4)=3.33642d-01
      e2(5)=1.90652d+00
      cs2(5)=6.16497d-01
      cp2(5)=4.44940d-01
      e2(6)=7.05887d-01
      cs2(6)=3.64829d-01
      cp2(6)=2.69254d-01
      e3(1)=9.29340d-01
      cs3(1)=-2.12290d-01
      cp3(1)=-2.24192d-02
      e3(2)=2.69035d-01
      cs3(2)=-1.07985d-01
      cp3(2)=1.92270d-01
      e3(3)=1.17379d-01
      cs3(3)=1.17584d+00
      cp3(3)=8.46181d-01
      e3(4)=4.21061d-02
      cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0146d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  aluminium                                       *****
c
  100 e4(1) = 0.325d0
      cd4(1)= 1.000d0
      go to(110,110,110,110,110,115),ngauss
  110 call caserr2('requested basis set for aluminium not available')
  115 continue
c        -----  66-31g                                          -----
      e1(1)=1.39831d+04
      cs1(1)=1.94267d-03
      e1(2)=2.09875d+03
      cs1(2)=1.48599d-02
      e1(3)=4.77705d+02
      cs1(3)=7.28494d-02
      e1(4)=1.34360d+02
      cs1(4)=2.46830d-01
      e1(5)=4.28709d+01
      cs1(5)=4.87258d-01
      e1(6)=1.45189d+01
      cs1(6)=3.23496d-01
      e2(1)=2.39668d+02
      cs2(1)=-2.92619d-03
      cp2(1)=4.60285d-03
      e2(2)=5.74419d+01
      cs2(2)=-3.74083d-02
      cp2(2)=3.31990d-02
      e2(3)=1.82859d+01
      cs2(3)=-1.14487d-01
      cp2(3)=1.36282d-01
      e2(4)=6.59914d+00
      cs2(4)=1.15635d-01
      cp2(4)=3.30476d-01
      e2(5)=2.49049d+00
      cs2(5)=6.12595d-01
      cp2(5)=4.49146d-01
      e2(6)=9.44545d-01
      cs2(6)=3.93799d-01
      cp2(6)=2.65704d-01
      e3(1)=1.27790d+00
      cs3(1)=-2.27606d-01
      cp3(1)=-1.75126d-02
      e3(2)=3.97590d-01
      cs3(2)=1.44583d-03
      cp3(2)=2.44533d-01
      e3(3)=1.60095d-01
      cs3(3)=1.09279d+00
      cp3(3)=8.04934d-01
      e3(4)=5.56577d-02
      cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0318d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  silicon                                         *****
c
  120 e4(1) = 0.450d0
      cd4(1)= 1.000d0
      go to(125,125,125,125,125,130),ngauss
  125 call caserr2('requested basis set for silicon not available')
  130 continue
c        ----- gordon's 66-31g                                  -----
      e1(1)=1.61921d+04
      cs1(1)=1.94924d-03
      e1(2)=2.43609d+03
      cs1(2)=1.48559d-02
      e1(3)=5.56001d+02
      cs1(3)=7.25689d-02
      e1(4)=1.56813d+02
      cs1(4)=2.45655d-01
      e1(5)=5.01692d+01
      cs1(5)=4.86060d-01
      e1(6)=1.70300d+01
      cs1(6)=3.25720d-01
      e2(1)=2.93350d+02
      cs2(1)=-2.82991d-03
      cp2(1)=4.43334d-03
      e2(2)=7.01173d+01
      cs2(2)=-3.60737d-02
      cp2(2)=3.24402d-02
      e2(3)=2.24301d+01
      cs2(3)=-1.16808d-01
      cp2(3)=1.33719d-01
      e2(4)=8.19425d+00
      cs2(4)=9.35768d-02
      cp2(4)=3.26780d-01
      e2(5)=3.14768d+00
      cs2(5)=6.01705d-01
      cp2(5)=4.51139d-01
      e2(6)=1.21515d+00
      cs2(6)=4.22207d-01
      cp2(6)=2.64105d-01
      e3(1)=1.65370d+00
      cs3(1)=-2.40600d-01
      cp3(1)=-1.51774d-02
      e3(2)=5.40760d-01
      cs3(2)=7.37953d-02
      cp3(2)=2.75139d-01
      e3(3)=2.04406d-01
      cs3(3)=1.04094d+00
      cp3(3)=7.83008d-01
      e3(4)=7.23837d-02
      cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0331d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  phosphorous                                     *****
c
  140 e4(1) = 0.550d0
      cd4(1)= 1.000d0
      go to (160,160,160,220,160,225),ngauss
  160 call caserr2('requested basis set for phosphorus not available')
c
c     44-31g
c
 220  e1(1) = 3.018671780d+03
      cs1(1) = 1.852131370d-02
      e1(2) = 4.551271210d+02
      cs1(2) = 1.299048640d-01
      e1(3) = 1.023147300d+02
      cs1(3) = 4.551002880d-01
      e1(4) = 2.761784730d+01
      cs1(4) = 5.331318610d-01
      e2(1) = 1.144294010d+02
      cs2(1) = -2.475029610d-02
      cp2(1) = 2.741400250d-02
      e2(2) = 2.658229590d+01
      cs2(2) = -1.350924600d-01
      cp2(2) = 1.690791420d-01
      e2(3) = 7.871888900d+00
      cs2(3) = 2.277360800d-01
      cp2(3) = 4.691020890d-01
      e2(4) = 2.487857250d+00
      cs2(4) = 8.755931160d-01
      cp2(4) = 5.181530590d-01
      e3(1) = 5.075061900d+01
      cs3(1) = -4.511922300d-02
      cp3(1) = 3.779071180d-03
      e3(2) = 1.672862420d+00
      cs3(2) = -8.504729900d-01
      cp3(2) = -4.634384050d-02
      e3(3) = 6.210974120d-01
      cs3(3) = 1.596285850d+00
      cp3(3) = 1.033944290d+00
      e3(4) = 1.670160070d-01
      go to 227
  225 continue
c        -----  66-31g                                          -----
      e1(1)=1.94133d+04
      cs1(1)=1.85160d-03
      e1(2)=2.90942d+03
      cs1(2)=1.42062d-02
      e1(3)=6.61364d+02
      cs1(3)=6.99995d-02
      e1(4)=1.85759d+02
      cs1(4)=2.40079d-01
      e1(5)=5.91943d+01
      cs1(5)=4.84762d-01
      e1(6)=2.00310d+01
      cs1(6)=3.35200d-01
      e2(1)=3.39478d+02
      cs2(1)=-2.78217d-03
      cp2(1)=4.56462d-03
      e2(2)=8.10101d+01
      cs2(2)=-3.60499d-02
      cp2(2)=3.36936d-02
      e2(3)=2.58780d+01
      cs2(3)=-1.16631d-01
      cp2(3)=1.39755d-01
      e2(4)=9.45221d+00
      cs2(4)=9.68328d-02
      cp2(4)=3.39362d-01
      e2(5)=3.66566d+00
      cs2(5)=6.14418d-01
      cp2(5)=4.50921d-01
      e2(6)=1.46746d+00
      cs2(6)=4.03798d-01
      cp2(6)=2.38586d-01
      e3(1)=2.15623d+00
      cs3(1)=-2.52923d-01
      cp3(1)=-1.77653d-02
      e3(2)=7.48997d-01
      cs3(2)=3.28517d-02
      cp3(2)=2.74058d-01
      e3(3)=2.83145d-01
      cs3(3)=1.08125d+00
      cp3(3)=7.85421d-01
      e3(4)=9.98317d-02
 227  cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0348d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  sulfur                                          *****
c
  280 e4(1) = 0.650d0
      cd4(1)= 1.000d0
      go to (300,300,300,360,300,310),ngauss
  300 call caserr2('requested basis set for sulphur not available')
c
c        *****  44-31g                                          *****
c
  360 e1(1) = 3.442124410d+03
      cs1(1) = 1.849212360d-02
      e1(2) = 5.189131000d+02
      cs1(2) = 1.298220220d-01
      e1(3) = 1.166909030d+02
      cs1(3) = 4.550417870d-01
      e1(4) = 3.157164720d+01
      cs1(4) = 5.330083560d-01
      e2(1) = 1.274405760d+02
      cs2(1) = -2.726461060d-02
      cp2(1) = 2.915199950d-02
      e2(2) = 2.974766730d+01
      cs2(2) = -1.424834150d-01
      cp2(2) = 1.779596760d-01
      e2(3) = 8.834664280d+00
      cs2(3) = 2.597043520d-01
      cp2(3) = 4.836237120d-01
      e2(4) = 2.817389820d+00
      cs2(4) = 8.525472950d-01
      cp2(4) = 4.942553020d-01
      e3(1) = 3.729185370d+00
      cs3(1) = -2.775315230d-01
      cp3(1) = -3.375092630d-02
      e3(2) = 1.406770170d+00
      cs3(2) = -4.576434550d-01
      cp3(2) = 1.457110450d-01
      e3(3) = 5.481099690d-01
      cs3(3) = 1.431684270d+00
      cp3(3) = 8.982887430d-01
      e3(4) = 1.703809050d-01
      go to 315
  310 continue
c        -----  66-31g                                          -----
      e1(1)=2.19171d+04
      cs1(1)=1.86924d-03
      e1(2)=3.30149d+03
      cs1(2)=1.42303d-02
      e1(3)=7.54146d+02
      cs1(3)=6.96962d-02
      e1(4)=2.12711d+02
      cs1(4)=2.38487d-01
      e1(5)=6.79896d+01
      cs1(5)=4.83307d-01
      e1(6)=2.30515d+01
      cs1(6)=3.38074d-01
      e2(1)=4.23735d+02
      cs2(1)=-2.37677d-03
      cp2(1)=4.06101d-03
      e2(2)=1.00710d+02
      cs2(2)=-3.16930d-02
      cp2(2)=3.06813d-02
      e2(3)=3.21599d+01
      cs2(3)=-1.13317d-01
      cp2(3)=1.30452d-01
      e2(4)=1.18079d+01
      cs2(4)=5.60900d-02
      cp2(4)=3.27205d-01
      e2(5)=4.63110d+00
      cs2(5)=5.92255d-01
      cp2(5)=4.52851d-01
      e2(6)=1.87025d+00
      cs2(6)=4.55006d-01
      cp2(6)=2.56042d-01
      e3(1)=2.61584d+00
      cs3(1)=-2.50374d-01
      cp3(1)=-1.45105d-02
      e3(2)=9.22167d-01
      cs3(2)=6.69570d-02
      cp3(2)=3.10263d-01
      e3(3)=3.41287d-01
      cs3(3)=1.05451d+00
      cp3(3)=7.54483d-01
      e3(4)=1.17167d-01
 315  cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0405d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  chlorine                                        *****
c
  420 e4(1) = 0.75d0
      cd4(1)= 1.000d0
      go to (440,440,440,500,440,445),ngauss
  440 call caserr2('requested basis set for chlorine not available')
c
c        *****  44-31g                                          *****
c
  500 e1(1) = 3.910302690d+03
      cs1(1) = 1.837943110d-02
      e1(2) = 5.895518070d+02
      cs1(2) = 1.291401230d-01
      e1(3) = 1.325939240d+02
      cs1(3) = 4.540448900d-01
      e1(4) = 3.590354250d+01
      cs1(4) = 5.344394360d-01
      e2(1) = 1.477653530d+02
      cs2(1) = -2.674332300d-02
      cp2(1) = 2.886446880d-02
      e2(2) = 3.450607530d+01
      cs2(2) = -1.446911820d-01
      cp2(2) = 1.779646700d-01
      e2(3) = 1.028647150d+01
      cs2(3) = 2.517035690d-01
      cp2(3) = 4.869998070d-01
      e2(4) = 3.311147380d+00
      cs2(4) = 8.598203810d-01
      cp2(4) = 4.890184500d-01
      e3(1) = 4.280284910d+00
      cs3(1) = -2.703962750d-01
      cp3(1) = -3.670288510d-02
      e3(2) = 1.641016670d+00
      cs3(2) = -3.416297190d-01
      cp3(2) = 1.918492420d-01
      e3(3) = 6.144785030d-01
      cs3(3) = 1.350024480d+00
      cp3(3) = 8.643376810d-01
      e3(4) = 1.956594110d-01
      go to 447
  445 continue
c        -----  66-31g                                          -----
      e1(1)=2.51801d+04
      cs1(1)=1.83296d-03
      e1(2)=3.78035d+03
      cs1(2)=1.40342d-02
      e1(3)=8.60474d+02
      cs1(3)=6.90974d-02
      e1(4)=2.42145d+02
      cs1(4)=2.37452d-01
      e1(5)=7.73349d+01
      cs1(5)=4.83034d-01
      e1(6)=2.62470d+01
      cs1(6)=3.39856d-01
      e2(1)=4.91765d+02
      cs2(1)=-2.29739d-03
      cp2(1)=3.98940d-03
      e2(2)=1.16984d+02
      cs2(2)=-3.07137d-02
      cp2(2)=3.03177d-02
      e2(3)=3.74153d+01
      cs2(3)=-1.12528d-01
      cp2(3)=1.29880d-01
      e2(4)=1.37834d+01
      cs2(4)=4.50163d-02
      cp2(4)=3.27951d-01
      e2(5)=5.45215d+00
      cs2(5)=5.89353d-01
      cp2(5)=4.53527d-01
      e2(6)=2.22588d+00
      cs2(6)=4.65206d-01
      cp2(6)=2.52154d-01
      e3(1)=3.18649d+00
      cs3(1)=-2.51830d-01
      cp3(1)=-1.42993d-02
      e3(2)=1.14427d+00
      cs3(2)=6.15890d-02
      cp3(2)=3.23572d-01
      e3(3)=4.20377d-01
      cs3(3)=1.06018d+00
      cp3(3)=7.43507d-01
      e3(4)=1.42657d-01
  447 cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.0483d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
c        *****  argon                                           *****
c
  550 e4(1) = 0.85d0
      cd4(1)= 1.000d0
      go to (551,551,551,551,551,552),ngauss
  551 call caserr2('requested basis set for argon not available')
  552 continue
c        -----  66-31g                                          -----
      e1(1)=2.83483d+04
      cs1(1)=1.82526d-03
      e1(2)=4.25762d+03
      cs1(2)=1.39686d-02
      e1(3)=9.69857d+02
      cs1(3)=6.87073d-02
      e1(4)=2.73263d+02
      cs1(4)=2.36204d-01
      e1(5)=8.73695d+01
      cs1(5)=4.82214d-01
      e1(6)=2.96867d+01
      cs1(6)=3.42043d-01
      e2(1)=5.75891d+02
      cs2(1)=-2.15972d-03
      cp2(1)=3.80665d-03
      e2(2)=1.36816d+02
      cs2(2)=-2.90775d-02
      cp2(2)=2.92305d-02
      e2(3)=4.38098d+01
      cs2(3)=-1.10827d-01
      cp2(3)=1.26467d-01
      e2(4)=1.62094d+01
      cs2(4)=2.76999d-02
      cp2(4)=3.23510d-01
      e2(5)=6.46084d+00
      cs2(5)=5.77613d-01
      cp2(5)=4.54896d-01
      e2(6)=2.65114d+00
      cs2(6)=4.88688d-01
      cp2(6)=2.56630d-01
      e3(1)=3.86028d+00
      cs3(1)=-2.55592d-01
      cp3(1)=-1.59197d-02
      e3(2)=1.41373d+00
      cs3(2)=3.78066d-02
      cp3(2)=3.24646d-01
      e3(3)=5.16646d-01
      cs3(3)=1.08056d+00
      cp3(3)=7.43990d-01
      e3(4)=1.73888d-01
      cs3(4)=1.00000d+00
      cp3(4)=1.00000d+00
      e3(5)=0.05d0
      cs3(5)=1.00000d+00
      cp3(5)=1.00000d+00
      return
c
  560 continue
      if (opg_root()) then
         write(iwr,*)'*** nuclear charge = ',ia
      endif
      call caserr2('requested basis set not available')
      return
c
      end
      subroutine ethree(e1,e2,e3,e4,e5,e6,e7,e8,
     +         cs1,cs2,cs3,cs4,cs5,
     +             cp2,cp3,cp4,cp5,
     +             cd6,cd7,cf8,ngauss,ia)
c
c        ***** energy optimized 6-31g* basis for 
c        ***** k, ca and first row transition metals sc - zn
c        ***** V. Rassolov, J.A. Pople, M. Ratner and T.L. Windus, 
c        ***** J. Chem. Phys.  109 (4) 1223-1229, 1998. 
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e1(*),e2(*),e3(*),e4(*)
      dimension e5(*),e6(*),e7(*),e8(*)
      dimension cs1(*),cs2(*),cs3(*),cs4(*),cs5(*)
      dimension        cp2(*),cp3(*),cp4(*),cp5(*)
      dimension        cd6(*),cd7(*),cf8(*)
c
      ian = ia - 18
      go to (10,30,50,70,90,110,130,150,170,190,210,230),ian
c
c        *****  potassium (22s,16p,1d) -> [5s,4p,1d]
c
   10 go to (15,15,15,15,15,20),ngauss
   15 call caserr2('requested basis set for potassium not available')
   20 continue
c        -----  6-31g  ------
      e1(1)= 31594.42000000d0
      e1(2)=  4744.33000000d0
      e1(3)=  1080.41900000d0
      e1(4)=   304.23380000d0
      e1(5)=    97.24586000d0
      e1(6)=    33.02495000d0
      cs1(1)=    0.00182801d0
      cs1(2)=    0.01399403d0
      cs1(3)=    0.06887129d0
      cs1(4)=    0.23697600d0
      cs1(5)=    0.48290400d0
      cs1(6)=    0.34047950d0
c
      e2(1)= 622.76250000d0
      e2(2)= 147.88390000d0
      e2(3)=  47.32735000d0
      e2(4)=  17.51495000d0
      e2(5)=   6.92272200d0
      e2(6)=   2.76827700d0
      cs2(1)= -0.00250298d0
      cs2(2)= -0.03315550d0
      cs2(3)= -0.12263870d0
      cs2(4)=  0.05353643d0
      cs2(5)=  0.61938600d0
      cs2(6)=  0.43458780d0
      cp2(1)=  0.00409464d0
      cp2(2)=  0.03145199d0
      cp2(3)=  0.13515580d0
      cp2(4)=  0.33905000d0
      cp2(5)=  0.46294550d0
      cp2(6)=  0.22426380d0
c
      e3(1)=  11.84802000d0
      e3(2)=   4.07921100d0
      e3(3)=   1.76348100d0
      e3(4)=   0.78892700d0
      e3(5)=   0.35038700d0
      e3(6)=   0.14634400d0
      cs3(1)=  0.01277689d0
      cs3(2)=  0.20987670d0
      cs3(3)= -0.00309527d0
      cs3(4)= -0.55938840d0
      cs3(5)= -0.51347600d0
      cs3(6)= -0.06598035d0
      cp3(1)= -0.01221377d0
      cp3(2)= -0.00690054d0
      cp3(3)=  0.20074660d0
      cp3(4)=  0.42813320d0
      cp3(5)=  0.39701560d0
      cp3(6)=  0.11047180d0
c
      e4(1)=   0.71680100d0
      e4(2)=   0.23374100d0
      e4(3)=   0.03867500d0
      cs4(1)= -0.05237772d0
      cs4(2)= -0.27985030d0
      cs4(3)=  1.14154700d0
      cp4(1)=  0.03164300d0
      cp4(2)= -0.04046160d0
      cp4(3)=  1.01202900d0
c
      e5(1)= 0.01652100d0
      cs5(1)=1.00000d0
      cp5(1)=1.00000d0
c
      e6(1) = 0.200d0
      cd6(1)= 1.000d0
      return
c
c        *****  calcium (22s,16p1d) -> [5s,4p1d]
c
   30 go to (35,35,35,35,35,40),ngauss
   35 call caserr2('requested basis set for calcium not available')
   40 continue
c        -----  6-31g*  -----
      e1(1)= 35264.86000000d0
      e1(2)=  5295.50300000d0
      e1(3)=  1206.02000000d0
      e1(4)=   339.68390000d0
      e1(5)=   108.62640000d0
      e1(6)=    36.92103000d0
      cs1(1)=    0.00181350d0
      cs1(2)=    0.01388493d0
      cs1(3)=    0.06836162d0
      cs1(4)=    0.23561880d0
      cs1(5)=    0.48206390d0
      cs1(6)=    0.34298190d0
c
      e2(1)= 706.30960000d0
      e2(2)= 167.81870000d0
      e2(3)=  53.82558000d0
      e2(4)=  20.01638000d0
      e2(5)=   7.97027900d0
      e2(6)=   3.21205900d0
      cs2(1)=  0.00244823d0
      cs2(2)=  0.03241504d0
      cs2(3)=  0.12262190d0
      cs2(4)= -0.04316965d0
      cs2(5)= -0.61269950d0
      cs2(6)= -0.44875400d0
      cp2(1)=  0.00402037d0
      cp2(2)=  0.03100601d0
      cp2(3)=  0.13372790d0
      cp2(4)=  0.33679830d0
      cp2(5)=  0.46312810d0
      cp2(6)=  0.22575320d0
c
      e3(1)=  14.19518000d0
      e3(2)=   4.88082800d0
      e3(3)=   2.16039000d0
      e3(4)=   0.98789900d0
      e3(5)=   0.44951700d0
      e3(6)=   0.18738700d0
      cs3(1)=  0.01084500d0
      cs3(2)=  0.20883330d0
      cs3(3)=  0.03150338d0
      cs3(4)= -0.55265180d0
      cs3(5)= -0.54379970d0
      cs3(6)= -0.06669342d0
      cp3(1)= -0.01289621d0
      cp3(2)= -0.01025198d0
      cp3(3)=  0.19597810d0
      cp3(4)=  0.43579330d0
      cp3(5)=  0.39964520d0
      cp3(6)=  0.09713636d0
c
      e4(1)=   1.03227100d0
      e4(2)=   0.38117100d0
      e4(3)=   0.06513100d0
      cs4(1)= -0.04439720d0
      cs4(2)= -0.32845630d0
      cs4(3)=  1.16301000d0
      cp4(1)= -0.42986210d0
      cp4(2)=  0.00693583d0
      cp4(3)=  0.97059330d0
c
      e5(1)= 0.02601000d0
      cs5(1)=1.00000d0
      cp5(1)=1.00000d0
c
      e6(1) = 0.200d0
      cd6(1)= 1.000d0
      return
c
c        *****  scandium (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
   50 go to (55,55,55,55,55,60),ngauss
   55 call caserr2('requested basis set for scandium not available')
   60 continue
      e1(1)= 39088.98000000d0
      e1(2)=  5869.79200000d0
      e1(3)=  1336.91000000d0
      e1(4)=   376.60310000d0
      e1(5)=   120.46790000d0
      e1(6)=    40.98032000d0
      cs1(1)=    0.00180326d0
      cs1(2)=    0.01380769d0
      cs1(3)=    0.06800396d0
      cs1(4)=    0.23470990d0
      cs1(5)=    0.48156900d0
      cs1(6)=    0.34456520d0
c
      e2(1) =  786.28520000d0
      e2(2) =  186.88700000d0
      e2(3) =   60.00935000d0
      e2(4) =   22.25883000d0
      e2(5) =    8.88514900d0
      e2(6) =    3.60921100d0
      cs2(1)=    0.00245186d0
      cs2(2)=    0.03259579d0
      cs2(3)=    0.12382420d0
      cs2(4)=   -0.04359890d0
      cs2(5)=   -0.61771810d0
      cs2(6)=   -0.44328230d0
      cp2(1) =   0.00403953d0
      cp2(2) =   0.03122570d0
      cp2(3) =   0.13498330d0
      cp2(4) =   0.34247930d0
      cp2(5) =   0.46231130d0
      cp2(6) =   0.21775240d0
c
      e3(1) =   29.84355000d0
      e3(2) =    9.54238300d0
      e3(3) =    4.05679000d0
      e3(4) =    1.70470300d0
      e3(5) =    0.70623400d0
      e3(6) =    0.27953600d0
      cs3(1) =  -0.00258630d0
      cs3(2) =   0.07188424d0
      cs3(3) =   0.25032600d0
      cs3(4) =  -0.29910030d0
      cs3(5) =  -0.74468180d0
      cs3(6) =  -0.17997760d0
      cp3(1) =  -0.00609665d0
      cp3(2) =  -0.02628884d0
      cp3(3) =   0.05091001d0
      cp3(4) =   0.37980970d0
      cp3(5) =   0.51708830d0
      cp3(6) =   0.18297720d0
c
      e4(1) =  1.06560900d0
      e4(2) =  0.42593300d0
      e4(3) =  0.07632000d0
      cs4(1) = 0.06482978d0
      cs4(2) = 0.32537560d0
      cs4(3) =-1.17080600d0
      cp4(1) =-0.29384400d0
      cp4(2) = 0.09235323d0
      cp4(3) = 0.98479300d0
c
      e5(1)=    0.02959400d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  11.14701000d0
      e6(2) =   2.82104300d0
      e6(3) =   0.81962000d0
      cd6(1)=   0.08747672d0
      cd6(2)=   0.37956350d0
      cd6(3)=   0.71803930d0
      e7(1) =   0.22146800d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  titanium (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
   70 go to (75,75,75,75,75,80),ngauss
   75 call caserr2('requested basis set for titanium not available')
   80 continue
      e1(1)= 43152.95000000d0
      e1(2)=  6479.57100000d0
      e1(3)=  1475.67500000d0
      e1(4)=   415.69910000d0
      e1(5)=   133.00060000d0
      e1(6)=    45.27222000d0
      cs1(1)=    0.00179187d0
      cs1(2)=    0.01372392d0
      cs1(3)=    0.06762830d0
      cs1(4)=    0.23376420d0
      cs1(5)=    0.48106960d0
      cs1(6)=    0.34622800d0
c
      e2(1) =  874.68260000d0
      e2(2) =  207.97850000d0
      e2(3) =   66.87918000d0
      e2(4) =   24.87347000d0
      e2(5) =    9.96844100d0
      e2(6) =    4.06382600d0
      cs2(1)=    0.00243101d0
      cs2(2)=    0.03233027d0
      cs2(3)=    0.12425200d0
      cs2(4)=   -0.03903905d0
      cs2(5)=   -0.61717890d0
      cs2(6)=   -0.44730970d0
      cp2(1) =   0.00401768d0
      cp2(2) =   0.03113966d0
      cp2(3) =   0.13490770d0
      cp2(4) =   0.34316720d0
      cp2(5) =   0.46257600d0
      cp2(6) =   0.21546030d0
c
      e3(1) =   33.64363000d0
      e3(2) =   10.87565000d0
      e3(3) =    4.62822500d0
      e3(4) =    1.95012600d0
      e3(5) =    0.80945200d0
      e3(6) =    0.32047400d0
      cs3(1) =  -0.00294036d0
      cs3(2) =   0.07163103d0
      cs3(3) =   0.25289150d0
      cs3(4) =  -0.29664010d0
      cs3(5) =  -0.74322150d0
      cs3(6) =  -0.18535200d0
      cp3(1) =  -0.00631162d0
      cp3(2) =  -0.02697638d0
      cp3(3) =   0.05316847d0
      cp3(4) =   0.38455490d0
      cp3(5) =   0.51276620d0
      cp3(6) =   0.18111350d0
c
      e4(1) =  1.22414800d0
      e4(2) =  0.48426300d0
      e4(3) =  0.08409600d0
      cs4(1) = 0.06351465d0
      cs4(2) = 0.31514040d0
      cs4(3) =-1.16259500d0
      cp4(1) =-0.21120700d0
      cp4(2) = 0.07771998d0
      cp4(3) = 0.98982140d0
c
      e5(1)=    0.03203600d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  13.69085000d0
      e6(2) =   3.51315400d0
      e6(3) =   1.04043400d0
      cd6(1)=   0.08589418d0
      cd6(2)=   0.37846710d0
      cd6(3)=   0.71612390d0
      e7(1) =   0.28696200d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  vanadium (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
   90 go to (95,95,95,95,95,100),ngauss
   95 call caserr2('requested basis set for vanadium not available')
  100 continue
      e1(1)= 47354.33000000d0
      e1(2)=  7110.78700000d0
      e1(3)=  1619.59100000d0
      e1(4)=   456.33790000d0
      e1(5)=   146.06060000d0
      e1(6)=    49.75791000d0
      cs1(1)=    0.00178451d0
      cs1(2)=    0.01366754d0
      cs1(3)=    0.06736122d0
      cs1(4)=    0.23305520d0
      cs1(5)=    0.48063160d0
      cs1(6)=    0.34748020d0
c
      e2(1) =  968.14840000d0
      e2(2) =  230.28210000d0
      e2(3) =   74.14591000d0
      e2(4) =   27.64107000d0
      e2(5) =   11.11475000d0
      e2(6) =    4.54311300d0
      cs2(1)=    0.00241060d0
      cs2(2)=    0.03207243d0
      cs2(3)=    0.12459420d0
      cs2(4)=   -0.03482177d0
      cs2(5)=   -0.61673740d0
      cs2(6)=   -0.45098440d0
      cp2(1) =   0.00399501d0
      cp2(2) =   0.03104061d0
      cp2(3) =   0.13477470d0
      cp2(4) =   0.34372790d0
      cp2(5) =   0.46287590d0
      cp2(6) =   0.21355470d0
c
      e3(1) =   37.64050000d0
      e3(2) =   12.28238000d0
      e3(3) =    5.23336600d0
      e3(4) =    2.20895000d0
      e3(5) =    0.91788000d0
      e3(6) =    0.36341200d0
      cs3(1) =  -0.00323320d0
      cs3(2) =   0.07130744d0
      cs3(3) =   0.25438200d0
      cs3(4) =  -0.29338870d0
      cs3(5) =  -0.74156950d0
      cs3(6) =  -0.19094100d0
      cp3(1) =  -0.00649406d0
      cp3(2) =  -0.02753453d0
      cp3(3) =   0.05516284d0
      cp3(4) =   0.38796720d0
      cp3(5) =   0.50902580d0
      cp3(6) =   0.18038400d0
c
      e4(1) =  1.39278100d0
      e4(2) =  0.54391300d0
      e4(3) =  0.09147600d0
      cs4(1) = 0.06139703d0
      cs4(2) = 0.30611300d0
      cs4(3) =-1.15489000d0
      cp4(1) =-0.18912650d0
      cp4(2) = 0.08005453d0
      cp4(3) = 0.98773990d0
c
      e5(1)=    0.03431200d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  16.05025000d0
      e6(2) =   4.16006300d0
      e6(3) =   1.24326500d0
      cd6(1)=   0.08599899d0
      cd6(2)=   0.38029960d0
      cd6(3)=   0.71276590d0
      e7(1) =   0.34427700d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  chromium (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  110 go to (115,115,115,115,115,120),ngauss
  115 call caserr2('requested basis set for chromium not available')
  120 continue
      e1(1)= 51789.81000000d0
      e1(2)=  7776.84900000d0
      e1(3)=  1771.38500000d0
      e1(4)=   499.15880000d0
      e1(5)=   159.79820000d0
      e1(6)=    54.47021000d0
      cs1(1)=    0.00177618d0
      cs1(2)=    0.01360476d0
      cs1(3)=    0.06706925d0
      cs1(4)=    0.23231040d0
      cs1(5)=    0.48024100d0
      cs1(6)=    0.34876530d0
c
      e2(1) = 1064.32800000d0
      e2(2) =  253.21380000d0
      e2(3) =   81.60924000d0
      e2(4) =   30.48193000d0
      e2(5) =   12.29439000d0
      e2(6) =    5.03772200d0
      cs2(1)=    0.00239967d0
      cs2(2)=    0.03194886d0
      cs2(3)=    0.12508680d0
      cs2(4)=   -0.03221866d0
      cs2(5)=   -0.61722840d0
      cs2(6)=   -0.45259360d0
      cp2(1) =   0.00398700d0
      cp2(2) =   0.03104662d0
      cp2(3) =   0.13505180d0
      cp2(4) =   0.34488650d0
      cp2(5) =   0.46285710d0
      cp2(6) =   0.21104260d0
c
      e3(1) =   41.56291000d0
      e3(2) =   13.67627000d0
      e3(3) =    5.84439000d0
      e3(4) =    2.47160900d0
      e3(5) =    1.02830800d0
      e3(6) =    0.40725000d0
      cs3(1) =  -0.00345422d0
      cs3(2) =   0.07218428d0
      cs3(3) =   0.25448200d0
      cs3(4) =  -0.29345340d0
      cs3(5) =  -0.73854550d0
      cs3(6) =  -0.19471570d0
      cp3(1) =  -0.00672250d0
      cp3(2) =  -0.02806471d0
      cp3(3) =   0.05820028d0
      cp3(4) =   0.39169880d0
      cp3(5) =   0.50478230d0
      cp3(6) =   0.17902900d0
c
      e4(1) =  1.57146400d0
      e4(2) =  0.60558000d0
      e4(3) =  0.09856100d0
      cs4(1) = 0.05892219d0
      cs4(2) = 0.29760550d0
      cs4(3) =-1.14750600d0
      cp4(1) =-0.19301000d0
      cp4(2) = 0.09605620d0
      cp4(3) = 0.98176090d0
c
      e5(1)=    0.03645900d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  18.41930000d0
      e6(2) =   4.81266100d0
      e6(3) =   1.44644700d0
      cd6(1)=   0.08650816d0
      cd6(2)=   0.38266990d0
      cd6(3)=   0.70937720d0
      e7(1) =   0.40041300d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  manganese (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  130 go to (135,135,135,135,135,140),ngauss
  135 call caserr2('requested basis set for manganese not available')
  140 continue
      e1(1)= 56347.14000000d0
      e1(2)=  8460.94300000d0
      e1(3)=  1927.32500000d0
      e1(4)=   543.23430000d0
      e1(5)=   173.99050000d0
      e1(6)=    59.36005000d0
      cs1(1)=    0.00177158d0
      cs1(2)=    0.01357081d0
      cs1(3)=    0.06690605d0
      cs1(4)=    0.23185410d0
      cs1(5)=    0.47990460d0
      cs1(6)=    0.34957370d0
c
      e2(1) = 1165.41200000d0
      e2(2) =  277.32760000d0
      e2(3) =   89.47278000d0
      e2(4) =   33.48256000d0
      e2(5) =   13.54037000d0
      e2(6) =    5.55797200d0
      cs2(1)=    0.00238875d0
      cs2(2)=    0.03181708d0
      cs2(3)=    0.12546700d0
      cs2(4)=   -0.02955431d0
      cs2(5)=   -0.61751600d0
      cs2(6)=   -0.45444580d0
      cp2(1) =   0.00397732d0
      cp2(2) =   0.03103112d0
      cp2(3) =   0.13518940d0
      cp2(4) =   0.34573870d0
      cp2(5) =   0.46292050d0
      cp2(6) =   0.20905920d0
c
      e3(1) =   45.83532000d0
      e3(2) =   15.18777000d0
      e3(3) =    6.50071000d0
      e3(4) =    2.75158300d0
      e3(5) =    1.14540400d0
      e3(6) =    0.45368700d0
      cs3(1) =  -0.00366586d0
      cs3(2) =   0.07231971d0
      cs3(3) =   0.25444860d0
      cs3(4) =  -0.29103800d0
      cs3(5) =  -0.73598600d0
      cs3(6) =  -0.19976170d0
      cp3(1) =  -0.00688758d0
      cp3(2) =  -0.02846816d0
      cp3(3) =   0.06031832d0
      cp3(4) =   0.39389610d0
      cp3(5) =   0.50137690d0
      cp3(6) =   0.17922640d0
c
      e4(1) =  1.75799900d0
      e4(2) =  0.66702200d0
      e4(3) =  0.10512900d0
      cs4(1) = 0.05628572d0
      cs4(2) = 0.28974910d0
      cs4(3) =-1.14065300d0
      cp4(1) =-0.50350240d0
      cp4(2) = 0.23450110d0
      cp4(3) = 0.91412570d0
c
      e5(1)=    0.03841800d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  20.94355000d0
      e6(2) =   5.51048600d0
      e6(3) =   1.66503800d0
      cd6(1)=   0.08672702d0
      cd6(2)=   0.38418830d0
      cd6(3)=   0.70690710d0
      e7(1) =   0.46173300d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  iron (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  150 go to (155,155,155,155,155,160),ngauss
  155 call caserr2('requested basis set for iron not available')
  160 continue
      e1(1)= 61132.62000000d0
      e1(2)=  9179.34200000d0
      e1(3)=  2090.85700000d0
      e1(4)=   589.24790000d0
      e1(5)=   188.75430000d0
      e1(6)=    64.44629000d0
      cs1(1)=    0.00176611d0
      cs1(2)=    0.01353038d0
      cs1(3)=    0.06673128d0
      cs1(4)=    0.23148230d0
      cs1(5)=    0.47970580d0
      cs1(6)=    0.35019760d0
c
      e2(1) = 1259.98000000d0
      e2(2) =  299.87610000d0
      e2(3) =   96.84917000d0
      e2(4) =   36.31020000d0
      e2(5) =   14.72996000d0
      e2(6) =    6.06607500d0
      cs2(1)=    0.00243801d0
      cs2(2)=    0.03224048d0
      cs2(3)=    0.12657240d0
      cs2(4)=   -0.03139902d0
      cs2(5)=   -0.62075930d0
      cs2(6)=   -0.45029140d0
      cp2(1) =   0.00402802d0
      cp2(2) =   0.03144647d0
      cp2(3) =   0.13683170d0
      cp2(4) =   0.34872360d0
      cp2(5) =   0.46179310d0
      cp2(6) =   0.20430580d0
c
      e3(1) =   50.43485000d0
      e3(2) =   16.83929000d0
      e3(3) =    7.19208600d0
      e3(4) =    3.05342000d0
      e3(5) =    1.27364300d0
      e3(6) =    0.50409100d0
      cs3(1) =  -0.00387326d0
      cs3(2) =   0.07196598d0
      cs3(3) =   0.25565910d0
      cs3(4) =  -0.28828370d0
      cs3(5) =  -0.73428220d0
      cs3(6) =  -0.20493530d0
      cp3(1) =  -0.00701713d0
      cp3(2) =  -0.02877660d0
      cp3(3) =   0.06181383d0
      cp3(4) =   0.39549460d0
      cp3(5) =   0.49890590d0
      cp3(6) =   0.17912510d0
c
      e4(1) =  1.95031600d0
      e4(2) =  0.73672100d0
      e4(3) =  0.11417700d0
      cs4(1) = 0.05694869d0
      cs4(2) = 0.28829150d0
      cs4(3) =-1.13815900d0
      cp4(1) =-0.45937960d0
      cp4(2) = 0.28521390d0
      cp4(3) = 0.90764850d0
c
      e5(1)=    0.04114800d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  23.14994000d0
      e6(2) =   6.12236800d0
      e6(3) =   1.84660100d0
      cd6(1)=   0.08876935d0
      cd6(2)=   0.38963190d0
      cd6(3)=   0.70148160d0
      e7(1) =   0.50436100d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  cobalt (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  170 go to (175,175,175,175,175,180),ngauss
  175 call caserr2('requested basis set for cobalt not available')
  180 continue
c        -----  6-31g*  -----
      e1(1)= 66148.99000000d0
      e1(2)=  9933.07700000d0
      e1(3)=  2262.81600000d0
      e1(4)=   637.91540000d0
      e1(5)=   204.41220000d0
      e1(6)=    69.82538000d0
      cs1(1)=    0.00175979d0
      cs1(2)=    0.01348162d0
      cs1(3)=    0.06649342d0
      cs1(4)=    0.23079390d0
      cs1(5)=    0.47929190d0
      cs1(6)=    0.35140970d0
c
      e2(1)= 1378.84100000d0
      e2(2)=  328.26940000d0
      e2(3)=  106.09460000d0
      e2(4)=   39.83275000d0
      e2(5)=   16.18622000d0
      e2(6)=    6.66778800d0
      cs2(1)=   0.00237628d0
      cs2(2)=   0.03167450d0
      cs2(3)=   0.12628880d0
      cs2(4)=  -0.02584552d0
      cs2(5)=  -0.61834910d0
      cs2(6)=  -0.45670080d0
      cp2(1)=   0.00397149d0
      cp2(2)=   0.03108174d0
      cp2(3)=   0.13574390d0
      cp2(4)=   0.34768270d0
      cp2(5)=   0.46263400d0
      cp2(6)=   0.20516320d0
c
      e3(1)=   54.52355000d0
      e3(2)=   18.29783000d0
      e3(3)=    7.86734800d0
      e3(4)=    3.34053400d0
      e3(5)=    1.39375600d0
      e3(6)=    0.55132600d0
      cs3(1)=  -0.00399300d0
      cs3(2)=   0.07409663d0
      cs3(3)=   0.25420000d0
      cs3(4)=  -0.29216570d0
      cs3(5)=  -0.73187030d0
      cs3(6)=  -0.20407840d0
      cp3(1)=  -0.00729077d0
      cp3(2)=  -0.02926027d0
      cp3(3)=   0.06564150d0
      cp3(4)=   0.40006520d0
      cp3(5)=   0.49502360d0
      cp3(6)=   0.17582400d0
c
      e4(1)=    2.15194700d0
      e4(2)=    0.81106300d0
      e4(3)=    0.12101700d0
      cs4(1)=   0.05379843d0
      cs4(2)=   0.27599710d0
      cs4(3)=  -1.12969200d0
      cp4(1)=  -0.21654960d0
      cp4(2)=   0.12404880d0
      cp4(3)=   0.97240640d0
c
      e5(1)=    0.04303700d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  25.59306000d0
      e6(2) =   6.80099000d0
      e6(3) =   2.05164700d0
      cd6(1)=   0.09004748d0
      cd6(2)=   0.39317030d0
      cd6(3)=   0.69768440d0
      e7(1) =   0.55567100d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  nickel (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  190 go to (195,195,195,195,195,200),ngauss
  195 call caserr2('requested basis set for nickel not available')
  200 continue
c        -----  6-31g*  -----
      e1(1)= 71396.35000000d0
      e1(2)= 10720.84000000d0
      e1(3)=  2442.12900000d0
      e1(4)=   688.42650000d0
      e1(5)=   220.61530000d0
      e1(6)=    75.39373000d0
      cs1(1)=    0.00175300d0
      cs1(2)=    0.01343122d0
      cs1(3)=    0.06627041d0
      cs1(4)=    0.23025080d0
      cs1(5)=    0.47901860d0
      cs1(6)=    0.35234440d0
c
      e2(1)= 1492.53200000d0
      e2(2)=  355.40130000d0
      e2(3)=  114.95340000d0
      e2(4)=   43.22043000d0
      e2(5)=   17.59710000d0
      e2(6)=    7.25776500d0
      cs2(1)=   0.00237071d0
      cs2(2)=   0.03160566d0
      cs2(3)=   0.12663350d0
      cs2(4)=  -0.02417037d0
      cs2(5)=  -0.61877750d0
      cs2(6)=  -0.45767700d0
      cp2(1)=   0.00396755d0
      cp2(2)=   0.03109479d0
      cp2(3)=   0.13595170d0
      cp2(4)=   0.34851360d0
      cp2(5)=   0.46254980d0
      cp2(6)=   0.20351860d0
c
      e3(1)=   59.35261000d0
      e3(2)=   20.02181000d0
      e3(3)=    8.61456100d0
      e3(4)=    3.66053100d0
      e3(5)=    1.52811100d0
      e3(6)=    0.60405700d0
      cs3(1)=  -0.00416200d0
      cs3(2)=   0.07425111d0
      cs3(3)=   0.25413600d0
      cs3(4)=  -0.29034770d0
      cs3(5)=  -0.73021210d0
      cs3(6)=  -0.20760570d0
      cp3(1)=  -0.00742145d0
      cp3(2)=  -0.02953410d0
      cp3(3)=   0.06731852d0
      cp3(4)=   0.40166600d0
      cp3(5)=   0.49266230d0
      cp3(6)=   0.17568930d0
c
      e4(1)=    2.37927600d0
      e4(2)=    0.88583900d0
      e4(3)=    0.12852900d0
      cs4(1)=   0.05157888d0
      cs4(2)=   0.27076110d0
      cs4(3)=  -1.12477000d0
      cp4(1)=  -0.18876630d0
      cp4(2)=   0.10151990d0
      cp4(3)=   0.97909060d0
c
      e5(1)=    0.04519500d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  28.19147000d0
      e6(2) =   7.52358400d0
      e6(3) =   2.27122800d0
      cd6(1)=   0.09098881d0
      cd6(2)=   0.39582080d0
      cd6(3)=   0.69471540d0
      e7(1) =   0.61160300d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
      return
c
c        *****  copper (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  210 go to (215,215,215,215,215,220),ngauss
  215 call caserr2('requested basis set for copper not available')
  220 continue
c        -----  6-31g*  -----
      e1(1)= 76794.38000000d0
      e1(2)= 11530.70000000d0
      e1(3)=  2626.57500000d0
      e1(4)=   740.49030000d0
      e1(5)=   237.35280000d0
      e1(6)=    81.15818000d0
      cs1(1)=    0.00174816d0
      cs1(2)=    0.01339602d0
      cs1(3)=    0.06610885d0
      cs1(4)=    0.22982650d0
      cs1(5)=    0.47876750d0
      cs1(6)=    0.35307390d0
c
      e2(1)= 1610.81400000d0
      e2(2)=  383.63670000d0
      e2(3)=  124.17330000d0
      e2(4)=   46.74678000d0
      e2(5)=   19.06569000d0
      e2(6)=    7.87156700d0
      cs2(1)=   0.00236406d0
      cs2(2)=   0.03153635d0
      cs2(3)=   0.12694520d0
      cs2(4)=  -0.02262840d0
      cs2(5)=  -0.61920800d0
      cs2(6)=  -0.45853930d0
      cp2(1)=   0.00396331d0
      cp2(2)=   0.03110223d0
      cp2(3)=   0.13613500d0
      cp2(4)=   0.34929140d0
      cp2(5)=   0.46247800d0
      cp2(6)=   0.20201020d0
c
      e3(1)=   64.45732000d0
      e3(2)=   21.85212000d0
      e3(3)=    9.40534300d0
      e3(4)=    3.99916800d0
      e3(5)=    1.67029700d0
      e3(6)=    0.65962700d0
      cs3(1)=  -0.00433108d0
      cs3(2)=   0.07412307d0
      cs3(3)=   0.25421080d0
      cs3(4)=  -0.28748430d0
      cs3(5)=  -0.72914360d0
      cs3(6)=  -0.21139510d0
      cp3(1)=  -0.00752373d0
      cp3(2)=  -0.02975687d0
      cp3(3)=   0.06849654d0
      cp3(4)=   0.40271410d0
      cp3(5)=   0.49084900d0
      cp3(6)=   0.17592680d0
c
      e4(1)=    2.60008800d0
      e4(2)=    0.96309400d0
      e4(3)=    0.13616100d0
      cs4(1)=   0.05027577d0
      cs4(2)=   0.26500400d0
      cs4(3)=  -1.12015500d0
      cp4(1)=  -0.17029110d0
      cp4(2)=   0.09310133d0
      cp4(3)=   0.98143360d0
c
      e5(1)=    0.04733200d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  30.85341000d0
      e6(2) =   8.26498500d0
      e6(3) =   2.49533200d0
      cd6(1)=   0.09199905d0
      cd6(2)=   0.39850210d0
      cd6(3)=   0.69178970d0
      e7(1) =   0.66765800d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
c
      return
c
c        *****  zinc (22s,16p,4d,1f) -> [5s,4p,2d,1f]
c
  230 go to (235,235,235,235,235,240),ngauss
  235 call caserr2('requested basis set for zinc not available')
  240 continue
c        -----  6-31g*  -----
      e1(1)= 82400.94000000d0
      e1(2)= 12372.55000000d0
      e1(3)=  2818.35100000d0
      e1(4)=   794.57170000d0
      e1(5)=   254.72320000d0
      e1(6)=    87.13880000d0
      cs1(1)=    0.00174333d0
      cs1(2)=    0.01335966d0
      cs1(3)=    0.06594365d0
      cs1(4)=    0.22941510d0
      cs1(5)=    0.47854530d0
      cs1(6)=    0.35377530d0
c
      e2(1)= 1732.56900000d0
      e2(2)=  412.71490000d0
      e2(3)=  133.67800000d0
      e2(4)=   50.38585000d0
      e2(5)=   20.58358000d0
      e2(6)=    8.50594000d0
      cs2(1)=   0.00236146d0
      cs2(2)=   0.03150177d0
      cs2(3)=   0.12727740d0
      cs2(4)=  -0.02145928d0
      cs2(5)=  -0.61976520d0
      cs2(6)=  -0.45901800d0
      cp2(1)=   0.00396313d0
      cp2(2)=   0.03113411d0
      cp2(3)=   0.13639310d0
      cp2(4)=   0.35012660d0
      cp2(5)=   0.46231790d0
      cp2(6)=   0.20049950d0
c
      e3(1)=   69.36492000d0
      e3(2)=   23.62082000d0
      e3(3)=   10.18471000d0
      e3(4)=    4.33408200d0
      e3(5)=    1.81091800d0
      e3(6)=    0.71484100d0
      cs3(1)=  -0.00444010d0
      cs3(2)=   0.07505253d0
      cs3(3)=   0.25331110d0
      cs3(4)=  -0.28818970d0
      cs3(5)=  -0.72670520d0
      cs3(6)=  -0.21334390d0
      cp3(1)=  -0.00768926d0
      cp3(2)=  -0.02997982d0
      cp3(3)=   0.07082411d0
      cp3(4)=   0.40461410d0
      cp3(5)=   0.48823250d0
      cp3(6)=   0.17519700d0
c
      e4(1)=    2.82384200d0
      e4(2)=    1.03954300d0
      e4(3)=    0.14326400d0
      cs4(1)=   0.04898543d0
      cs4(2)=   0.25927930d0
      cs4(3)=  -1.11571100d0
      cp4(1)=  -0.15867630d0
      cp4(2)=   0.08379327d0
      cp4(3)=   0.98405470d0
c
      e5(1)=    0.04929600d0
      cs5(1)=   1.00000d0
      cp5(1)=   1.00000d0
c
      e6(1) =  33.70764000d0
      e6(2) =   9.06110600d0
      e6(3) =   2.73838300d0
      cd6(1)=   0.09262648d0
      cd6(2)=   0.40029800d0
      cd6(3)=   0.68966080d0
      e7(1) =   0.73029400d0
      cd7(1)=   1.00000000d0
c
      e8(1) =   0.80000000d0
      cf8(1)=   1.00000000d0
c
      return
      end
      subroutine n21g(csinp,cpinp,cdinp,opol,oryd,sc,scc,
     + igauss,nucz,intyp,nangm,
     + nbfs,minf,maxf,loc,ngauss,ns,nshmax,ngsmax,ierr1,ierr2,
     + nat,iw)
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
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
      common/blkin/e123(21),cs123(21),cp23(21),cd3(21)
      dimension csinp(*),cpinp(*),cdinp(*)
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
      dimension scalf(9),scc(*)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      data pt75/0.75d+00/
      data pt5,done/0.5d+00,1.0d+00/
      data pi32/5.56832799683170d+00/
      data hscal,hescal/1.1d+00,1.0d+00/
      tol = 1.0d-10
c
c     ----- set up the n-21g basis set -----
c
c     h-ne   n=3,6  j.s.binkley, j.a.pople, w.j.hehre
c                   j. chem. phys. 1980, 102, 939-947.
c     na-ar  n=3,6  m.s.gordon, j.s.binkley, j.a.pople, w.j.pietro,
c                   w.j.hehre   j.am.chem.soc. 1982, 104, 2797-2803.
c     k,ca,ga-kr,rb,sr,in-xe   n=3    k.d.dobbs, w.j.hehre
c                   j.comput.chem. 7, 359-378(1986)
c     sc-zn  n=3    k.d.dobbs, w.j.hehre
c                   j.comput.chem. 8, 861-879(1987)
c     y-cd   n=3    k.d.dobbs, w.j.hehre
c                   j.comput.chem. 8, 880-893(1987)
c
c     the bases are in the order 1s,2sp,3sp,3d,4sp,4d,5sp
c
c     ----- get exponents, contraction coefficients, and
c          scale factors. -----
c
c
      call vclr(e123,1,84)
      if (nucz .le. 2) then
c
c     ----- hydrogen and helium -----
c
      call dn21g(e123,cs123,cp23(igauss+2),cd3(igauss+2),
     +           nucz,igauss)
      mxpass = 3
      if (opol) mxpass = mxpass + 1
c
      else if (nucz .le. 10) then
c
c     ----- lithium to neon -----
c
      call dn21g(e123,cs123,cp23(igauss+1),cd3(igauss+5),
     +           nucz,igauss)
      mxpass = 4
      if (opol) mxpass = mxpass + 1
c
      else if (nucz .le. 18) then
c
c     ----- sodium to argon  -----
c
      call dn21g(e123,cs123,cp23(igauss+1),cd3(2*igauss+5),
     +           nucz,igauss)
      mxpass = 5
      if (opol) mxpass = mxpass + 1
c
      else if (nucz .le. 30) then
c
c     ----- potassium to zinc
c
      call t1n21g(e123,cs123,cp23(igauss+1),cd3(3*igauss+4),
     * nucz,igauss,iw)
      mxpass = 5
      if(nucz.gt.20) mxpass=mxpass+2
c
      else if (nucz .le. 36) then
c
c     ----- gallium to krypton
c
      call en21g(e123,cs123,cp23(igauss+1),cd3(3*igauss+4),
     * nucz,igauss)
      mxpass = 6
      else if (nucz .le. 48) then
c
c     ----- rubidium to cadmium
c
      call t2n21g(e123,cs123,cp23(igauss+1),cd3(4*igauss+4),
     * nucz,igauss)
      mxpass = 7
      if(nucz.gt.38) mxpass=mxpass+2
c
      else if (nucz .le. 54) then
c
c     ----- indium to xenon
c
      call fn21g(e123,cs123,cp23(igauss+1),cd3(4*igauss+4),
     * nucz,igauss)
      mxpass = 8
      else
         if (opg_root()) then
            write(iw,*)'*** nuclear charge = ',nucz
         endif
         call caserr2('requested basis set not available')
      endif
      if(nucz.gt.18.and.opol) then
         if (opg_root()) then
            write(iw,*)'*** nuclear charge = ',nucz
         endif
         call caserr2('standard polarisation set not available')
      endif
      if(nucz.gt.18.and.oryd) then
         if (opg_root()) then
            write(iw,*)'*** nuclear charge = ',nucz
         endif
         call caserr2('standard diffuse set not available')
      endif
      if(nucz.eq.1) then
        scalf(1) = hscal
        scalf(2) = hscal
        scalf(3) = 1.0d0
        scalf(4) = 1.0d0
      else if(nucz.eq.2) then
        scalf(1) = hescal
        scalf(2) = hescal
        scalf(3) = 1.0d0
        scalf(4) = 1.0d0
      else
        scalf(1)=sc
        do 172 loop=2,mxpass
  172   scalf(loop)=scc(loop-1)
        do 167 loop=1,mxpass
  167   if(scalf(loop).eq.0.0d0) scalf(loop)=done
      endif
c
c     ----- loop over 2,3,4,7 or 9 shells. -----
c
      ipass = 0
  180 ipass = ipass+1
      if (nucz .gt. 2) go to 240
c
c     ----- hydrogen and helium-----
c
      go to (200,220,230,235),ipass
  200 ityp = 1
      igau = 2
      ig = 0
      go to 420
  220 ityp = 1
      igau = 1
      ig = 2
      go to 420
  230 ityp = 1
      igau = 1
      ig =  3
      if (oryd) then
       go to 420
      else
       go to 560
      endif
  235 ityp = 17
      igau = 1
      ig = 4
      go to 420
  240 if (nucz .gt. 10) go to 320
c
c     ----- lithium to neon -----
c
      go to (260,280,300,310,315),ipass
  260 ityp = 1
      igau = igauss
      ig = 0
      go to 420
  280 ityp = 4
      igau = 2
      ig = igauss
      go to 420
  300 ityp = 4
      igau = 1
      ig = igauss + 2
      go to 420
  310 ityp = 22
      igau = 1
      ig = igauss + 3
      if (oryd) then
       go to 420
      else
       go to 560
      endif
  315 ityp = 18
      igau = 1
      ig = igauss + 4
      go to 420
c
  320 if (nucz .gt. 18) go to 590
c
c     ----- sodium to argon -----
c
      go to (340,360,380,400,410,415),ipass
  340 ityp = 1
      igau = igauss
      ig = 0
      go to 420
  360 ityp = 4
      ig = igauss
      go to 420
  380 ityp = 7
      igau = 2
      ig = igauss + igauss
      go to 420
  400 ityp = 7
      igau = 1
      ig = igauss + igauss + 2
      go to 420
  410 ityp = 22
      igau = 1
      ig = igauss + igauss+  3
      if (oryd) then
       go to 420
      else
       go to 560
      endif
  415 ityp = 18
      igau = 1
      ig = igauss + igauss + 4
      go to 420
c
  590 if(nucz.gt.36) go to 275
c
c     ----- potassium to krypton -----
c
      go to (600,610,620,630,640,650,660),ipass
 600  ityp=1
      igau=igauss
      ig=0
      go to 420
 610  ityp=4
      ig=igauss
      go to 420
 620  ityp=7
      ig=igauss+igauss
      go to 420
 630  ityp=11
      igau=2
      ig=igauss*3
      go to 420
 640  ityp=11
      igau=1
      ig=igauss*3+2
      go to 420
 650  ityp=8
      ig=igauss*3+3
        if(nucz.gt.30) then
        igau=3
        else
        igau=2
        endif
      go to 420
 660  ityp=8
      igau=1
      ig=igauss*3+5
      go to 420
 275  if(nucz.gt.54) then
         if (opg_root()) then
            write(iw,*)'*** nuclear charge = ',nucz
         endif
         call caserr2('requested basis set not available')
      endif
      go to (800,810,820,830,840,850,860,870,880),ipass
 800  ityp=1
      igau=igauss
      ig=0
      go to 420
 810  ityp=4
      ig=igauss
      go to 420
 820  ityp=7
      ig=igauss+igauss
      go to 420
 830  ityp=11
      ig=igauss*3
      go to 420
 840  ityp=24
      igau=2
      ig=igauss*4
      go to 420
 850  ityp=24
      igau=1
      ig=igauss*4+2
      go to 420
 860  ityp=8
      igau=3
      ig=igauss*4+3
      go to 420
 870  ityp=12
      ig=igauss*5+3
       if(nucz.gt.48) then
       igau=3
       else
       igau=2
       endif
      go to 420
 880  ityp=12
      igau=1
      ig=igauss*5+5
  420 continue
      nshell = nshell+1
      if(nshell.gt.nshmax) ierr1=1
      if(ierr1.ne.0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igau
      kloc(nshell) = loc+1
      ngauss = ngauss+igau
      if(ngauss.gt.ngsmax) ierr2=1
      if(ierr2.ne.0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do 440 i = 1,igau
      k = k1+i-1
      ex(k) = e123(ig+i)*scalf(ipass)**2
      csinp(k) = cs123(ig+i)
      cpinp(k) = cp23(ig+i)
      cdinp(k) = cd3(ig+i)
      cs(k) = csinp(k)
      cp(k) = cpinp(k)
  440 cd(k) = cdinp(k)
c
c     ----- always unnormalize primitives... -----
c
      do 460 k = k1,k2
      ee = ex(k)+ex(k)
      facs = pi32/(ee*dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      cs(k) = cs(k)/dsqrt(facs)
      cp(k) = cp(k)/dsqrt(facp)
  460 cd(k) = cd(k)/dsqrt(facd)
c
c     ----- if(normf.eq.0) normalize basis functions. -----
c
      if (normf .eq. 1) go to 540
      facs = 0.0d0
      facp = 0.0d0
      facd = 0.0d0
      do 500 ig = k1,k2
      do 500 jg = k1,ig
      ee = ex(ig)+ex(jg)
      fac = ee*dsqrt(ee)
      dums = cs(ig)*cs(jg)/fac
      dump = pt5*cp(ig)*cp(jg)/(ee*fac)
      dumd = pt75*cd(ig)*cd(jg)/(ee*ee*fac)
      if (ig .eq. jg) go to 480
      dums = dums+dums
      dump = dump+dump
      dumd = dumd+dumd
  480 facs = facs+dums
      facp = facp+dump
      facd = facd+dumd
  500 continue
c
      do 520 ig = k1,k2
      if (facs .gt. tol) cs(ig) = cs(ig)/dsqrt(facs*pi32)
      if (facp .gt. tol) cp(ig) = cp(ig)/dsqrt(facp*pi32)
  520 if (facd .gt. tol) cd(ig) = cd(ig)/dsqrt(facd*pi32)
  540 continue
  560 if (ipass .lt. mxpass) go to 180
c
      return
      end
      subroutine t1n21g(e,cs,cp,cd,ia,ngauss,iwr)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),cs(*),cp(*),cd(*)
      iaa=ia-18
      go to (10,20,100,200,300,400,500,600,700,800,900,1000),iaa
  10  go to ( 11, 11, 12, 11, 11, 11),ngauss
  11  call caserr2('requested basis set not available for potassium')
  12  e( 1)=1721.1755d+00
      e( 2)=260.01633d+00
      e( 3)=56.624554d+00
      e( 4)=71.55720d+00
      e( 5)=15.43894d+00
      e( 6)=4.474551d+00
      e( 7)=4.121275d+00
      e( 8)=1.188621d+00
      e( 9)=0.3756738d+00
      e(10)=2.445766d-01
      e(11)=3.897175d-02
      e(12)=1.606255d-02
      cs( 1)=0.0648747d+00
      cs( 2)=0.3808593d+00
      cs( 3)=0.6773681d+00
      cs( 4)=-0.1093429d+00
      cs( 5)=0.1130640d+00
      cs( 6)=0.9462575d+00
      cs( 7)=-0.2699730d+00
      cs( 8)=0.3646323d+00
      cs( 9)=0.8107533d+00
      cs(10)=-2.688250d-01
      cs(11)=1.128983d+00
      cs(12)=1.0d+00
      cp( 1)=0.1339654d+00
      cp( 2)=0.5302673d+00
      cp( 3)=0.5117992d+00
      cp( 4)=0.01994922d+00
      cp( 5)=0.4340213d+00
      cp( 6)=0.6453226d+00
      cp( 7)=3.081035d-04
      cp( 8)=9.998787d-01
      cp( 9)=1.0d+00
      go to 2380
  20  go to ( 21, 21, 22, 21, 21, 21),ngauss
  21  call caserr2('requested basis set not available for calcium')
  22  e( 1)=1915.4348d+00
      e( 2)=289.53324d+00
      e( 3)=63.106352d+00
      e( 4)=80.39744d+00
      e( 5)=17.33075d+00
      e( 6)=5.083624d+00
      e( 7)=4.782229d+00
      e( 8)=1.462558d+00
      e( 9)=0.4792230d+00
      e(10)=4.396824d-01
      e(11)=5.913040d-02
      e(12)=2.389701d-02
      cs( 1)=0.0646237d+00
      cs( 2)=0.3798376d+00
      cs( 3)=0.6783294d+00
      cs( 4)=-0.1093028d+00
      cs( 5)=0.1088996d+00
      cs( 6)=0.9492768d+00
      cs( 7)=-0.2816074d+00
      cs( 8)=0.3410510d+00
      cs( 9)=0.8381044d+00
      cs(10)=-2.697049d-01
      cs(11)=1.113293d+00
      cs(12)=1.0d+00
      cp( 1)=0.1354332d+00
      cp( 2)=0.5372215d+00
      cp( 3)=0.5018044d+00
      cp( 4)=0.01900928d+00
      cp( 5)=0.4360377d+00
      cp( 6)=0.6386709d+00
      cp( 7)=3.081105d-04
      cp( 8)=9.998964d-01
      cp( 9)=1.0d+00
      go to 2380
  100 go to ( 110, 110, 120, 110, 110, 110),ngauss
 110  call caserr2('requested basis set not available for scandium')
 120  e( 1) = 2.119887d+03
      e( 2) = 3.204299d+02
      e( 3) = 6.989893d+01
      e( 4) = 8.976450d+01
      e( 5) = 1.938510d+01
      e( 6) = 5.731423d+00
      e( 7) = 5.491938d+00
      e( 8) = 1.743742d+00
      e( 9) = 5.662273d-01
      e(10) = 5.168015d-01
      e(11) = 6.721404d-02
      e(12) = 2.598452d-02
      e(13) = 5.722215d+00
      e(14) = 1.360849d+00
      e(15) = 3.226516d-01
      cs( 1) = 6.442080d-02
      cs( 2) = 3.791603d-01
      cs( 3) = 6.789629d-01
      cs( 4) = -1.093837d-01
      cs( 5) = 1.050699d-01
      cs( 6) = 9.522045d-01
      cs( 7) = -2.852107d-01
      cs( 8) = 3.241555d-01
      cs( 9) = 8.565921d-01
      cs(10) = -2.626780d-01
      cs(11) = 1.108079d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.363278d-01
      cp( 2) = 5.418598d-01
      cp( 3) = 4.950551d-01
      cp( 4) = 1.761356d-02
      cp( 5) = 4.336448d-01
      cp( 6) = 6.425507d-01
      cp( 7) = 3.270567d-04
      cp( 8) = 9.998935d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.652364d-01
      cd( 2) = 8.558605d-01
      cd( 3) = 0.100000d+01
      go to 2380
  200 go to ( 210, 210, 220, 210, 210, 210),ngauss
 210  call caserr2('requested basis set not available for titanium')
 220  e( 1) = 2.335020d+03
      e( 2) = 3.530441d+02
      e( 3) = 7.705845d+01
      e( 4) = 9.957387d+01
      e( 5) = 2.154671d+01
      e( 6) = 6.413965d+00
      e( 7) = 6.238279d+00
      e( 8) = 1.996108d+00
      e( 9) = 6.464899d-01
      e(10) = 5.732849d-01
      e(11) = 7.311942d-02
      e(12) = 2.653794d-02
      e(13) = 7.083666d+00
      e(14) = 1.709634d+00
      e(15) = 4.141225d-01
      cs( 1) = 6.421660d-02
      cs( 2) = 3.784120d-01
      cs( 3) = 6.796813d-01
      cs( 4) = -1.094719d-01
      cs( 5) = 1.019427d-01
      cs( 6) = 9.546224d-01
      cs( 7) = -2.861372d-01
      cs( 8) = 3.218728d-01
      cs( 9) = 8.595511d-01
      cs(10) = -2.424500d-01
      cs(11) = 1.100075d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.372966d-01
      cp( 2) = 5.458753d-01
      cp( 3) = 4.890681d-01
      cp( 4) = 1.923665d-02
      cp( 5) = 4.404422d-01
      cp( 6) = 6.356195d-01
      cp( 7) = 2.920158d-04
      cp( 8) = 9.999067d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.629210d-01
      cd( 2) = 8.557721d-01
      cd( 3) = 0.100000d+01
      go to 2380
  300 go to ( 310, 310, 320, 310, 310, 310),ngauss
  310 call caserr2('requested basis set not available for vanadium')
 320  e( 1) = 2.563877d+03
      e( 2) = 3.875340d+02
      e( 3) = 8.459823d+01
      e( 4) = 1.097938d+02
      e( 5) = 2.376921d+01
      e( 6) = 7.122961d+00
      e( 7) = 6.981204d+00
      e( 8) = 2.219839d+00
      e( 9) = 7.198030d-01
      e(10) = 6.312620d-01
      e(11) = 8.006166d-02
      e(12) = 2.886489d-02
      e(13) = 8.342917d+00
      e(14) = 2.032944d+00
      e(15) = 4.957115d-01
      cs( 1) = 6.394750d-02
      cs( 2) = 3.775940d-01
      cs( 3) = 6.805421d-01
      cs( 4) = -1.098355d-01
      cs( 5) = 1.007070d-01
      cs( 6) = 9.556327d-01
      cs( 7) = -2.884588d-01
      cs( 8) = 3.364357d-01
      cs( 9) = 8.481903d-01
      cs(10) = -2.364899d-01
      cs(11) = 1.097721d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.384210d-01
      cp( 2) = 5.504894d-01
      cp( 3) = 4.824165d-01
      cp( 4) = 2.182075d-02
      cp( 5) = 4.567616d-01
      cp( 6) = 6.186750d-01
      cp( 7) = 1.899536d-04
      cp( 8) = 9.999396d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.640620d-01
      cd( 2) = 8.539665d-01
      cd( 3) = 0.100000d+01
      go to 2380
  400 go to ( 410, 410, 420, 410, 410, 410),ngauss
  410 call caserr2('requested basis set not available for chromium')
 420  e( 1) = 2.798294d+03
      e( 2) = 4.231370d+02
      e( 3) = 9.243886d+01
      e( 4) = 1.202806d+02
      e( 5) = 2.603727d+01
      e( 6) = 7.844172d+00
      e( 7) = 7.793276d+00
      e( 8) = 2.497196d+00
      e( 9) = 8.051419d-01
      e(10) = 7.039206d-01
      e(11) = 8.616195d-02
      e(12) = 3.219882d-02
      e(13) = 9.625339d+00
      e(14) = 2.362264d+00
      e(15) = 5.770944d-01
      cs( 1) = 6.382380d-02
      cs( 2) = 3.770840d-01
      cs( 3) = 6.809889d-01
      cs( 4) = -1.177790d-01
      cs( 5) = 1.014311d-01
      cs( 6) = 9.571981d-01
      cs( 7) = -2.888567d-01
      cs( 8) = 3.351147d-01
      cs( 9) = 8.502481d-01
      cs(10) = -2.322508d-01
      cs(11) = 1.093671d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.398782d-01
      cp( 2) = 5.559834d-01
      cp( 3) = 4.748183d-01
      cp( 4) = 2.218471d-02
      cp( 5) = 4.616250d-01
      cp( 6) = 6.145386d-01
      cp( 7) = 1.799645d-04
      cp( 8) = 9.999448d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.655959d-01
      cd( 2) = 8.521557d-01
      cd( 3) = 0.100000d+01
      go to 2380
  500 go to ( 510, 510, 520, 510, 510, 510),ngauss
  510 call caserr2('requested basis set not available for manganese')
 520  e( 1) = 3.041686d+03
      e( 2) = 4.600901d+02
      e( 3) = 1.005958d+02
      e( 4) = 1.317673d+02
      e( 5) = 2.856915d+01
      e( 6) = 8.660501d+00
      e( 7) = 8.569081d+00
      e( 8) = 2.768178d+00
      e( 9) = 8.872882d-01
      e(10) = 7.674426d-01
      e(11) = 9.202527d-02
      e(12) = 3.326490d-02
      e(13) = 1.106884d+01
      e(14) = 2.730707d+00
      e(15) = 6.685095d-01
      cs( 1) = 6.374490d-02
      cs( 2) = 3.767490d-01
      cs( 3) = 6.812474d-01
      cs( 4) = -1.102964d-01
      cs( 5) = 9.818963d-02
      cs( 6) = 9.576595d-01
      cs( 7) = -2.917135d-01
      cs( 8) = 3.439630d-01
      cs( 9) = 8.451975d-01
      cs(10) = -2.300039d-01
      cs(11) = 1.091450d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.404540d-01
      cp( 2) = 5.578022d-01
      cp( 3) = 4.715006d-01
      cp( 4) = 2.422379d-02
      cp( 5) = 4.686598d-01
      cp( 6) = 6.074211d-01
      cp( 7) = 3.078886d-04
      cp( 8) = 9.999074d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.652718d-01
      cd( 2) = 8.517945d-01
      cd( 3) = 0.100000d+01
      go to 2380
  600 go to ( 610, 610, 620, 610, 610, 610),ngauss
 610  call caserr2('requested basis set not available for iron')
 620  e( 1) = 3.299184d+03
      e( 2) = 4.990886d+02
      e( 3) = 1.091614d+02
      e( 4) = 1.434652d+02
      e( 5) = 3.116858d+01
      e( 6) = 9.483612d+00
      e( 7) = 9.464565d+00
      e( 8) = 3.100373d+00
      e( 9) = 9.864930d-01
      e(10) = 8.534123d-01
      e(11) = 9.881222d-02
      e(12) = 3.644214d-02
      e(13) = 1.235449d+01
      e(14) = 3.055605d+00
      e(15) = 7.385909d-01
      cs( 1) = 6.358590d-02
      cs( 2) = 3.762016d-01
      cs( 3) = 6.817845d-01
      cs( 4) = -1.105517d-01
      cs( 5) = 9.684681d-02
      cs( 6) = 9.587974d-01
      cs( 7) = -2.920555d-01
      cs( 8) = 3.375236d-01
      cs( 9) = 8.519416d-01
      cs(10) = -2.279441d-01
      cs(11) = 1.088287d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.411006d-01
      cp( 2) = 5.603874d-01
      cp( 3) = 4.676444d-01
      cp( 4) = 2.376201d-02
      cp( 5) = 4.689113d-01
      cp( 6) = 6.083113d-01
      cp( 7) =-4.262652d-04
      cp( 8) = 1.000124d+00
      cp( 9) = 0.100000d+01
      cd( 1) = 2.686110d-01
      cd( 2) = 8.492717d-01
      cd( 3) = 0.100000d+01
      go to 2380
  700 go to ( 710 , 710 , 720 , 710 , 710 , 710 ),ngauss
  710 call caserr2('requested basis set not available for cobalt')
 720  e( 1) = 3.564762d+03
      e( 2) = 5.393908d+02
      e( 3) = 1.180449d+02
      e( 4) = 1.554382d+02
      e( 5) = 3.381561d+01
      e( 6) = 1.033323d+01
      e( 7) = 1.038152d+01
      e( 8) = 3.382714d+00
      e( 9) = 1.076954d+00
      e(10) = 9.090155d-01
c original     e(11) = 1.050406d-02
c AdM
      e(11) = 1.050406d-01
      e(12) = 3.725658d-02
      e(13) = 1.374070d+01
      e(14) = 3.408983d+00
      e(15) = 8.186409d-01
      cs( 1) = 6.348660d-02
      cs( 2) = 3.758181d-01
      cs( 3) = 6.821217d-01
      cs( 4) = -1.109867d-01
      cs( 5) = 9.676742d-02
      cs( 6) = 9.589921d-01
      cs( 7) = -2.922622d-01
      cs( 8) = 3.435207d-01
      cs( 9) = 8.469634d-01
      cs(10) = -2.174599d-01
      cs(11) = 1.084998d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.420642d-01
      cp( 2) = 5.634439d-01
      cp( 3) = 4.630244d-01
      cp( 4) = 2.631326d-02
      cp( 5) = 4.769170d-01
      cp( 6) = 5.991543d-01
      cp( 7) = 2.284428d-04
      cp( 8) = 9.999337d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.709550d-01
      cd( 2) = 8.473421d-01
      cd( 3) = 0.100000d+01
      go to 2380
 800  go to (810,810, 820, 810,810, 810),ngauss
 810  call caserr2('requested basis set not available for nickel')
 820  e( 1) = 3.848005d+03
      e( 2) = 5.820307d+02
      e( 3) = 1.273674d+02
      e( 4) = 1.682896d+02
      e( 5) = 3.665633d+01
      e( 6) = 1.123212d+01
      e( 7) = 1.135877d+01
      e( 8) = 3.738846d+00
      e( 9) = 1.182463d+00
      e(10) = 9.889038d-01
      e(11) = 1.110250d-01
      e(12) = 3.925822d-02
      e(13) = 1.522069d+01
      e(14) = 3.786020d+00
      e(15) = 9.045574d-01
      cs( 1) = 6.326610d-02
      cs( 2) = 3.751710d-01
c cs(3) is given in Hehre (J.Comput.Chem 10)
      cs( 3) = 6.828238d-01
      cs( 4) = -1.111151d-01
      cs( 5) = 9.532380d-02
      cs( 6) = 9.601613d-01
      cs( 7) = -2.920604d-01
      cs( 8) = 3.375407d-01
      cs( 9) = 8.525330d-01
      cs(10) = -2.136872d-01
      cs(11) = 1.081933d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.424905d-01
      cp( 2) = 5.655470d-01
      cp( 3) = 4.599926d-01
c original      cp( 4) = 1.613762d-02
c correction by Hehre (J.Comput.Chem 10)
      cp( 4) = 2.613762d-02
      cp( 5) = 4.765980d-01
      cp( 6) = 6.003798d-01
      cp( 7) = 2.943514d-04
      cp( 8) = 9.999170d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.726060d-01
      cd( 2) = 8.459279d-01
      cd( 3) = 0.100000d+01
      go to 2380
 900  go to (910,910,920,910,910,910),ngauss
 910  call caserr2('requested basis set not available for cupper')
 920  write (iwr,'(A)') ' Cu basis refit AdM/JvL'
      e( 1) = 4.134302d+03
      e( 2) = 6.254912d+02
      e( 3) = 1.369556d+02
      e( 4) = 1.814960d+02
      e( 5) = 3.957431d+01
      e( 6) = 1.216246d+01
      e( 7) = 1.235111d+01
      e( 8) = 4.049651d+00
      e( 9) = 1.279225d+00
*     e(10) = 9.943270d-01
*     e(11) = 1.169329d-01
*     e(12) = 7.506963d-04
c begin AdM
      e(10) = 0.106860d+01
      e(11) = 0.117067d+00
      e(12) = 0.417710d-01
c end AdM
      e(13) = 1.675938d+01
      e(14) = 4.178977d+00
      e(15) = 9.943270d-01
      cs( 1) = 6.318780d-02
      cs( 2) = 3.748448d-01
      cs( 3) = 6.831002d-01
      cs( 4) = -1.113198d-01
      cs( 5) = 9.448679d-02
      cs( 6) = 9.608790d-01
      cs( 7) = -2.922231d-01
      cs( 8) = 3.429909d-01
      cs( 9) = 8.479463d-01
*     cs(10) = -2.065079d-01
*     cs(11) = 1.079273d+00
c begin AdM
      cs(10) =  -0.206583d+00
      cs(11) = 0.107825d+01
c end AdM
      cs(12) = 0.100000d+01
      cp( 1) = 1.430844d-01
      cp( 2) = 5.677561d-01
      cp( 3) = 4.567141d-01
      cp( 4) = 2.772714d-02
      cp( 5) = 4.835244d-01
      cp( 6) = 5.929779d-01
*     cp( 7) = 1.385032d-04
*     cp( 8) = 9.999613d-01
c begin AdM
      cp( 7) = -0.120234d-03
      cp( 8) = 0.100003d+01
c end AdM
      cp( 9) = 0.100000d+01
      cd( 1) = 2.741125d-01
      cd( 2) = 8.446245d-01
      cd( 3) = 0.100000d+01
      go to 2380
 1000 go to (1010,1010,1020,1010,1010,1010),ngauss
 1010 call caserr2('requested basis set not available for zinc')
 1020 e( 1) = 4.432288d+03
      e( 2) = 6.706601d+02
      e( 3) = 1.469024d+02
      e( 4) = 1.950042d+02
      e( 5) = 4.256889d+01
      e( 6) = 1.312143d+01
      e( 7) = 1.340231d+01
      e( 8) = 4.399906d+00
      e( 9) = 1.385148d+00
      e(10) = 1.121558d+00
      e(11) = 1.229436d-01
      e(12) = 4.219327d-02
      e(13) = 1.836820d+01
      e(14) = 4.591304d+00
      e(15) = 1.090203d+00
      cs( 1) = 6.309280d-02
      cs( 2) = 3.745038d-01
      cs( 3) = 6.834160d-01
      cs( 4) = -1.116283d-01
      cs( 5) = 9.433553d-02
      cs( 6) = 9.611002d-01
      cs( 7) = -2.917811d-01
      cs( 8) = 3.426145d-01
      cs( 9) = 8.482840d-01
      cs(10) = -2.023706d-01
      cs(11) = 1.077035d+00
      cs(12) = 0.100000d+01
      cp( 1) = 1.438055d-01
      cp( 2) = 5.700019d-01
      cp( 3) = 4.533119d-01
      cp( 4) = 2.870528d-02
      cp( 5) = 4.862515d-01
      cp( 6) = 5.902353d-01
      cp( 7) = 3.440941d-04
      cp( 8) = 9.999053d-01
      cp( 9) = 0.100000d+01
      cd( 1) = 2.753856d-01
      cd( 2) = 8.434773d-01
      cd( 3) = 0.100000d+01
 2380 continue
      return
      end
      subroutine en21g(e,cs,cp,cd,ia,ngauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),cs(*),cp(*),cd(*)
      iaa=ia-30
      go to (300,400,500,600,700,800),iaa
  300 go to ( 310, 310, 320, 310, 310, 310),ngauss
 310  call caserr2('requested basis set not available for gallium')
 320  e( 1)=4751.8979d+00
      e( 2)=718.92054d+00
      e( 3)=157.44592d+00
      e( 4)=209.5834d+00
      e( 5)=45.69171d+00
      e( 6)=14.13297d+00
      e( 7)=14.59954d+00
      e( 8)=4.860842d+00
      e( 9)=1.549111d+00
      e(10)=1.267943d+00
      e(11)=0.1883995d+00
      e(12)=0.05723676d+00
      e(13)=21.292530d+00
      e(14)=5.3931662d+00
      e(15)=1.3338828d+00
      cs( 1)=0.0628396d+00
      cs( 2)=0.3736112d+00
      cs( 3)=0.6843626d+00
      cs( 4)=-0.1115162d+00
      cs( 5)=0.09269636d+00
      cs( 6)=0.9622870d+00
      cs( 7)=0.2910292d+00
      cs( 8)=-0.3231876d+00
      cs( 9)=-0.8643910d+00
      cs(10)=-0.2851306d+00
      cs(11)=1.128022d+00
      cs(12)=1.0d+00
      cp( 1)=0.1442658d+00
      cp( 2)=0.5731775d+00
      cp( 3)=0.4490858d+00
      cp( 4)=0.02656186d+00
      cp( 5)=0.4833137d+00
      cp( 6)=0.5924304d+00
      cp( 7)=0.03018346d+00
      cp( 8)=0.9884658d+00
      cp( 9)=1.0d+00
      cd( 1)=0.1619895d+00
      cd( 2)=0.5116739d+00
      cd( 3)=0.5898732d+00
      go to 2380
  400 go to ( 410, 410, 420, 410, 410, 410),ngauss
 410  call caserr2('requested basis set not available for germanium')
 420  e( 1)=5073.7499d+00
      e( 2)=767.72417d+00
      e( 3)=168.18881d+00
      e( 4)=224.4360d+00
      e( 5)=48.95543d+00
      e( 6)=15.18370d+00
      e( 7)=15.91257d+00
      e( 8)=5.441437d+00
      e( 9)=1.742603d+00
      e(10)=1.466538d+00
      e(11)=0.2630934d+00
      e(12)=0.08482072d+00
      e(13)=24.321421d+00
      e(14)=6.2238135d+00
      e(15)=1.5887375d+00
      cs( 1)=0.0627249d+00
      cs( 2)=0.3731671d+00
      cs( 3)=0.6847867d+00
      cs( 4)=-0.1115150d+00
      cs( 5)=0.09120021d+00
      cs( 6)=0.9634491d+00
      cs( 7)=-0.2895652d+00
      cs( 8)=0.2938828d+00
      cs( 9)=0.8891993d+00
      cs(10)=-0.3967339d+00
      cs(11)=1.190670d+00
      cs(12)=1.0d+00
      cp( 1)=0.1446395d+00
      cp( 2)=0.5753796d+00
      cp( 3)=0.4459949d+00
      cp( 4)=0.02297302d+00
      cp( 5)=0.4732446d+00
      cp( 6)=0.6032779d+00
      cp( 7)=0.02789294d+00
      cp( 8)=0.9874901d+00
      cp( 9)=1.0d+00
      cd( 1)=0.1577985d+00
      cd( 2)=0.5114922d+00
      cd( 3)=0.5857703d+00
      go to 2380
  500 go to ( 510, 510, 520, 510, 510, 510),ngauss
  510 call caserr2('requested basis set not available for arsenic')
  520 e( 1)=5407.6138d+00
      e( 2)=818.17436d+00
      e( 3)=179.26569d+00
      e( 4)=237.7783d+00
      e( 5)=54.25662d+00
      e( 6)=16.32803d+00
      e( 7)=17.10185d+00
      e( 8)=5.805144d+00
      e( 9)=1.902084d+00
      e(10)=1.675404d+00
      e(11)=0.3416557d+00
      e(12)=0.1136303d+00
      e(13)=27.437209d+00
      e(14)=7.0840440d+00
      e(15)=1.8558226d+00
      cs( 1)=0.0626011d+00
      cs( 2)=0.3727790d+00
      cs( 3)=0.6851842d+00
      cs( 4)=-0.1128384d+00
      cs( 5)=0.08722744d+00
      cs( 6)=0.9681883d+00
      cs( 7)=-0.2914537d+00
      cs( 8)=0.2969619d+00
      cs( 9)=0.8865791d+00
      cs(10)=-0.5057610d+00
      cs(11)=1.251764d+00
      cs(12)=1.0d+00
      cp( 1)=0.1496798d+00
      cp( 2)=0.5623223d+00
      cp( 3)=0.4593235d+00
      cp( 4)=0.02568559d+00
      cp( 5)=0.4833968d+00
      cp( 6)=0.5887854d+00
      cp( 7)=0.02528246d+00
      cp( 8)=0.9874328d+00
      cp( 9)=1.0d+00
      cd( 1)=0.1544952d+00
      cd( 2)=0.5114318d+00
      cd( 3)=0.5821935d+00
      go to 2380
  600 go to ( 610, 610, 620, 610, 610, 610),ngauss
  610 call caserr2('requested basis set not available for selenium')
  620 e( 1)=5751.3215d+00
      e( 2)=870.25721d+00
      e( 3)=190.72949d+00
      e( 4)=255.0164d+00
      e( 5)=55.57654d+00
      e( 6)=17.35661d+00
      e( 7)=18.44568d+00
      e( 8)=6.328759d+00
      e( 9)=2.096758d+00
      e(10)=1.872633d+00
      e(11)=0.4174736d+00
      e(12)=0.1370907d+00
      e(13)=30.627464d+00
      e(14)=7.9712764d+00
      e(15)=2.1348097d+00
      cs( 1)=0.0624934d+00
      cs( 2)=0.3723683d+00
      cs( 3)=0.6855799d+00
      cs( 4)=-0.1119076d+00
      cs( 5)=0.09099936d+00
      cs( 6)=0.9636682d+00
      cs( 7)=-0.2917925d+00
      cs( 8)=0.2846212d+00
      cs( 9)=0.8973052d+00
      cs(10)=-0.5677639d+00
      cs(11)=1.294127d+00
      cs(12)=1.0d+00
      cp( 1)=0.1461488d+00
      cp( 2)=0.5813714d+00
      cp( 3)=0.4374597d+00
      cp( 4)=0.02442141d+00
      cp( 5)=0.4833648d+00
      cp( 6)=0.5879038d+00
      cp( 7)=0.02825548d+00
      cp( 8)=0.9849060d+00
      cp( 9)=1.0d+00
      cd( 1)=0.1519858d+00
      cd( 2)=0.5116403d+00
      cd( 3)=0.5786936d+00
      go to 2380
  700 go to ( 710, 710, 720, 710, 710, 710),ngauss
  710 call caserr2('requested basis set not available for bromine')
  720 e( 1)=6103.2899d+00
      e( 2)=923.69743d+00
      e( 3)=202.52031d+00
      e( 4)=270.6015d+00
      e( 5)=58.25357d+00
      e( 6)=18.46933d+00
      e( 7)=19.76142d+00
      e( 8)=6.821752d+00
      e( 9)=2.291629d+00
      e(10)=2.131206d+00
      e(11)=0.4993537d+00
      e(12)=0.1647637d+00
      e(13)=33.965097d+00
      e(14)=8.9008312d+00
      e(15)=2.4284360d+00
      cs( 1)=0.0624175d+00
      cs( 2)=0.3720414d+00
      cs( 3)=0.6858728d+00
      cs( 4)=-0.1121487d+00
      cs( 5)=0.09314451d+00
      cs( 6)=0.9616794d+00
      cs( 7)=-0.2938704d+00
      cs( 8)=0.2802663d+00
      cs( 9)=0.9020357d+00
      cs(10)=-0.6518031d+00
      cs(11)=1.336012d+00
      cs(12)=1.000000d+00
      cp( 1)=0.1477514d+00
      cp( 2)=0.6010557d+00
      cp( 3)=0.4128704d+00
      cp( 4)=0.02500708d+00
      cp( 5)=0.4866098d+00
      cp( 6)=0.5824234d+00
      cp( 7)=0.02870833d+00
      cp( 8)=0.9840695d+00
      cp( 9)=1.0000000d+00
      cd( 1)=0.1496666d+00
      cd( 2)=0.5117475d+00
      cd( 3)=0.5759148d+00
      go to 2380
  800 go to ( 810, 810, 820, 810, 810, 810),ngauss
  810 call caserr2('requested basis set not available for krypton')
  820 e( 1)=6446.6307d+00
      e( 2)=976.87570d+00
      e( 3)=214.47955d+00
      e( 4)=287.6446d+00
      e( 5)=62.62009d+00
      e( 6)=19.69174d+00
      e( 7)=21.12321d+00
      e( 8)=7.303286d+00
      e( 9)=2.488850d+00
      e(10)=2.361374d+00
      e(11)=0.5860160d+00
      e(12)=0.1944473d+00
      e(13)=37.368103d+00
      e(14)=9.8543131d+00
      e(15)=2.7327955d+00
      cs( 1)=0.0625398d+00
      cs( 2)=0.3721075d+00
      cs( 3)=0.6856107d+00
      cs( 4)=-0.1120607d+00
      cs( 5)=0.09013913d+00
      cs( 6)=0.9643301d+00
      cs( 7)=-0.2958173d+00
      cs( 8)=0.2792167d+00
      cs( 9)=0.9037303d+00
      cs(10)=-0.7202454d+00
      cs(11)=1.376846d+00
      cs(12)=1.0d+00
      cp( 1)=0.1475279d+00
      cp( 2)=0.5868918d+00
      cp( 3)=0.4295068d+00
      cp( 4)=0.02606955d+00
      cp( 5)=0.4922498d+00
      cp( 6)=0.5742737d+00
      cp( 7)=0.02877518d+00
      cp( 8)=0.9833391d+00
      cp( 9)=1.0d+00
      cd( 1)=0.1479466d+00
      cd( 2)=0.5121719d+00
      cd( 3)=0.5729498d+00
 2380 continue
      return
      end
      subroutine fn21g(e,cs,cp,cd,ia,ngauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),cs(*),cp(*),cd(*)
      iaa=ia-48
      go to (300,400,500,600,700,800),iaa
  300 go to ( 310, 310, 320, 310, 310, 310),ngauss
  310 call caserr2('requested basis set not available for indium')
  320 e( 1)=12214.547d+00
      e( 2)=1848.9136d+00
      e( 3)=406.36833d+00
      e( 4)=550.4423d+00
      e( 5)=119.7744d+00
      e( 6)=38.66927d+00
      e( 7)=47.02931d+00
      e( 8)=22.49642d+00
      e( 9)=6.697117d+00
      e(10)=6.572360d+00
      e(11)=2.502158d+00
      e(12)=0.9420246d+00
      e(13)=1.001221d+00
      e(14)=0.1659704d+00
      e(15)=0.05433974d+00
      e(16)=102.17356d+00
      e(17)=28.394632d+00
      e(18)=8.9248045d+00
      e(19)=4.5353637d+00
      e(20)=1.5371481d+00
      e(21)=0.49949226d+00
      cs( 1)=0.0612476d+00
      cs( 2)=0.3676754d+00
      cs( 3)=0.6901359d+00
      cs( 4)=-0.1127094d+00
      cs( 5)=0.08344350d+00
      cs( 6)=0.9696880d+00
      cs( 7)=-0.2758954d+00
      cs( 8)=0.05977348d+00
      cs( 9)=1.082148d+00
      cs(10)=-0.4284831d+00
      cs(11)=0.4633644d+00
      cs(12)=0.8219679d+00
      cs(13)=-0.4364172d+00
      cs(14)=1.189893d+00
      cs(15)=1.0d+00
      cp( 1)=0.1523703d+00
      cp( 2)=0.6096508d+00
      cp( 3)=0.3970249d+00
      cp( 4)=-0.1408485d+00
      cp( 5)=0.5290867d+00
      cp( 6)=0.6620681d+00
      cp( 7)=0.01091305d+00
      cp( 8)=0.5036759d+00
      cp( 9)=0.5581809d+00
      cp(10)=0.02316334d+00
      cp(11)=0.9903309d+00
      cp(12)=1.0d+00
      cd( 1)=0.1205559d+00
      cd( 2)=0.4884976d+00
      cd( 3)=0.5850190d+00
      cd( 4)=0.2508574d+00
      cd( 5)=0.5693113d+00
      cd( 6)=0.3840635d+00
      go to 2380
  400 go to ( 410, 410, 420, 410, 410, 410),ngauss
  410 call caserr2('requested basis set not available for tin')
  420 e( 1)=12741.674d+00
      e( 2)=1928.4692d+00
      e( 3)=423.80797d+00
      e( 4)=574.2875d+00
      e( 5)=124.9537d+00
      e( 6)=40.39576d+00
      e( 7)=48.80662d+00
      e( 8)=23.83588d+00
      e( 9)=7.048296d+00
      e(10)=6.973378d+00
      e(11)=2.693040d+00
      e(12)=1.025958d+00
      e(13)=1.131463d+00
      e(14)=0.2034092d+00
      e(15)=0.07056383d+00
      e(16)=108.05630d+00
      e(17)=30.131576d+00
      e(18)=9.5300359d+00
      e(19)=4.9626098d+00
      e(20)=1.7120829d+00
      e(21)=0.57719451d+00
      cs( 1)=0.0611353d+00
      cs( 2)=0.3672929d+00
      cs( 3)=0.6905447d+00
      cs( 4)=-0.1127462d+00
      cs( 5)=0.08286347d+00
      cs( 6)=0.9701505d+00
      cs( 7)=-0.2824534d+00
      cs( 8)=0.06605595d+00
      cs( 9)=1.081987d+00
      cs(10)=-0.4340356d+00
      cs(11)=0.4610286d+00
      cs(12)=0.8285579d+00
      cs(13)=-0.5252085d+00
      cs(14)=1.229226d+00
      cs(15)=1.0d+00
      cp( 1)=0.1525798d+00
      cp( 2)=0.6110106d+00
      cp( 3)=0.3951549d+00
      cp( 4)=-0.1509628d+00
      cp( 5)=0.5399677d+00
      cp( 6)=0.6604823d+00
      cp( 7)=0.01195130d+00
      cp( 8)=0.5067195d+00
      cp( 9)=0.5529106d+00
      cp(10)=0.02107053d+00
      cp(11)=0.9905914d+00
      cp(12)=1.0d+00
      cd( 1)=0.1198237d+00
      cd( 2)=0.4875910d+00
      cd( 3)=0.5849874d+00
      cd( 4)=0.2529487d+00
      cd( 5)=0.5727612d+00
      cd( 6)=0.3690387d+00
      go to 2380
  500 go to ( 510, 510, 520, 510, 510, 510),ngauss
  510 call caserr2('requested basis set not available for antimony')
  520 e(1)=13289.383d+00
      e(2)=2010.5218d+00
      e(3)=441.69815d+00
      e(4)=598.8890d+00
      e(5)=130.0385d+00
      e(6)=42.13286d+00
      e(7)=51.51333d+00
      e(8)=24.43595d+00
      e(9)=7.420931d+00
      e(10)=7.314235d+00
      e(11)=2.844053d+00
      e(12)=1.105855d+00
      e(13)=1.278637d+00
      e(14)=0.2412321d+00
      e(15)=0.08662967d+00
      e(16)=115.80955d+00
      e(17)=32.305835d+00
      e(18)=10.250328d+00
      e(19)=5.4862102d+00
      e(20)=1.9216196d+00
      e(21)=0.66606265d+00
      cs(1)=0.0609843d+00
      cs(2)=0.3668487d+00
      cs(3)=0.6910501d+00
      cs(4)=-0.1127201d+00
      cs(5)=0.08264433d+00
      cs(6)=0.9702579d+00
      cs(7)=-0.2770433d+00
      cs(8)=0.05750323d+00
      cs(9)=1.084703d+00
      cs(10)=-0.4403812d+00
      cs(11)=0.4737341d+00
      cs(12)=0.822135d+00
      cs(13)=-0.6016951d+00
      cs(14)=1.258692d+00
      cs(15)=1.0d+00
      cp(1)=0.1530671d+00
      cp(2)=0.6135972d+00
      cp(3)=0.3916990d+00
      cp(4)=-0.1378699d+00
      cp(5)=0.536355d+00
      cp(6)=0.6508676d+00
      cp(7)=0.01530518d+00
      cp(8)=0.5160832d+00
      cp(9)=0.538757d+00
      cp(10)=0.02225276d+00
      cp(11)=0.9896434d+00
      cp(12)=1.0d+00
      cd( 1)=0.1166279d+00
      cd( 2)=0.4834363d+00
      cd( 3)=0.5901395d+00
      cd( 4)=0.2483656d+00
      cd( 5)=0.5743154d+00
      cd( 6)=0.3643044d+00
      go to 2380
  600 go to ( 610, 610, 620, 610, 610, 610),ngauss
 610  call caserr2('requested basis set not available for tellurium')
 620  e(1)=13796.56d+00
      e(2)=2088.8798d+00
      e(3)=459.39319d+00
      e(4)=623.2631d+00
      e(5)=135.36d+00
      e(6)=44.00048d+00
      e(7)=54.19078d+00
      e(8)=25.82039d+00
      e(9)=7.809583d+00
      e(10)=7.764217d+00
      e(11)=3.043932d+00
      e(12)=1.199253d+00
      e(13)=1.340364d+00
      e(14)=0.2780884d+00
      e(15)=0.09672607d+00
      e(16)=121.4083d+00
      e(17)=34.015217d+00
      e(18)=10.869138d+00
      e(19)=5.8031113d+00
      e(20)=2.0580658d+00
      e(21)=0.73283021d+00
      cs(1)=0.0610862d+00
      cs(2)=0.3669629d+00
      cs(3)=0.6907944d+00
      cs(4)=-0.11282d+00
      cs(5)=0.08225243d+00
      cs(6)=0.9706007d+00
      cs(7)=-0.274412d+00
      cs(8)=0.05182968d+00
      cs(9)=1.087622d+00
      cs(10)=-0.4467306d+00
      cs(11)=0.4694704d+00
      cs(12)=0.8298105d+00
      cs(13)=-0.5904699d+00
      cs(14)=1.281968d+00
      cs(15)=1.0d+00
      cp(1)=0.1534195d+00
      cp(2)=0.6148996d+00
      cp(3)=0.3895162d+00
      cp(4)=-0.1433312d+00
      cp(5)=0.5391881d+00
      cp(6)=0.6522851d+00
      cp(7)=0.0127462d+00
      cp(8)=0.5221231d+00
      cp(9)=0.5326614d+00
      cp(10)=0.02558834d+00
      cp(11)=0.9871016d+00
      cp(12)=1.0d+00
      cd( 1)=0.1169136d+00
      cd( 2)=0.4835555d+00
      cd( 3)=0.5883864d+00
      cd( 4)=0.2601939d+00
      cd( 5)=0.5797758d+00
      cd( 6)=0.3405580d+00
      go to 2380
  700 go to ( 710 , 710 , 720 , 710 , 710 , 710 ),ngauss
  710 call caserr2('requested basis set not available for iodine')
  720 e( 1)=14351.186d+00
      e( 2)=2173.0741d+00
      e( 3)=477.87205d+00
      e( 4)=648.1887d+00
      e( 5)=140.3064d+00
      e( 6)=45.69880d+00
      e( 7)=56.69469d+00
      e( 8)=27.48875d+00
      e( 9)=8.209096d+00
      e(10)=8.191679d+00
      e(11)=3.244596d+00
      e(12)=1.300489d+00
      e(13)=1.451380d+00
      e(14)=0.3281033d+00
      e(15)=0.1150759d+00
      e(16)=128.09026d+00
      e(17)=35.982378d+00
      e(18)=11.551116d+00
      e(19)=6.1461523d+00
      e(20)=2.2209370d+00
      e(21)=0.80991202d+00
      cs( 1)=0.0610028d+00
      cs( 2)=0.3666398d+00
      cs( 3)=0.6911306d+00
      cs( 4)=-0.1128507d+00
      cs( 5)=0.08322834d+00
      cs( 6)=0.9697516d+00
      cs( 7)=-0.2736965d+00
      cs( 8)=0.04649968d+00
      cs( 9)=1.091576d+00
      cs(10)=-0.4508277d+00
      cs(11)=0.4632094d+00
      cs(12)=0.8386360d+00
      cs(13)=-0.6658515d+00
      cs(14)=1.328584d+00
      cs(15)=1.0d+00
      cp( 1)=0.1541145d+00
      cp( 2)=0.6194618d+00
      cp( 3)=0.3837583d+00
      cp( 4)=-0.1523217d+00
      cp( 5)=0.5437686d+00
      cp( 6)=0.6561679d+00
      cp( 7)=0.01186988d+00
      cp( 8)=0.5265246d+00
      cp( 9)=0.5266076d+00
      cp(10)=0.02754131d+00
      cp(11)=0.9851368d+00
      cp(12)=1.0d+00
      cd( 1)=0.1158636d+00
      cd( 2)=0.4820494d+00
      cd( 3)=0.5894448d+00
      cd( 4)=0.2681817d+00
      cd( 5)=0.5800614d+00
      cd( 6)=0.3262263d+00
      go to 2380
 800  go to (810,810, 820, 810,810, 810),ngauss
 810  call caserr2('requested basis set not available for xenon')
 820  e( 1)=14902.236d+00
      e( 2)=2256.5383d+00
      e( 3)=496.37317d+00
      e( 4)=673.6611d+00
      e( 5)=145.8491d+00
      e( 6)=47.57708d+00
      e( 7)=59.16752d+00
      e( 8)=28.61159d+00
      e( 9)=8.596596d+00
      e(10)=8.638676d+00
      e(11)=3.462818d+00
      e(12)=1.401040d+00
      e(13)=1.578474d+00
      e(14)=0.3750814d+00
      e(15)=0.1331790d+00
      e(16)=134.91331d+00
      e(17)=37.956387d+00
      e(18)=12.227475d+00
      e(19)=6.6004928d+00
      e(20)=2.3980513d+00
      e(21)=0.88648239d+00
      cs( 1)=0.0609969d+00
      cs( 2)=0.3666290d+00
      cs( 3)=0.6911155d+00
      cs( 4)=-0.1129128d+00
      cs( 5)=0.08290529d+00
      cs( 6)=0.9700289d+00
      cs( 7)=-0.2739774d+00
      cs( 8)=0.04553006d+00
      cs( 9)=1.092528d+00
      cs(10)=-0.4558408d+00
      cs(11)=0.4617355d+00
      cs(12)=0.8442883d+00
      cs(13)=-0.7277718d+00
      cs(14)=1.362797d+00
      cs(15)=1.0d+00
      cp( 1)=0.1544275d+00
      cp( 2)=0.6206174d+00
      cp( 3)=0.3820041d+00
      cp( 4)=-0.1518572d+00
      cp( 5)=0.5471507d+00
      cp( 6)=0.6519540d+00
      cp( 7)=0.009585056d+00
      cp( 8)=0.5298193d+00
      cp( 9)=0.5235985d+00
      cp(10)=0.02807244d+00
      cp(11)=0.9842642d+00
      cp(12)=1.0d+00
      cd( 1)=0.1150105d+00
      cd( 2)=0.4815952d+00
      cd( 3)=0.5896133d+00
      cd( 4)=0.2718844d+00
      cd( 5)=0.5855569d+00
      cd( 6)=0.3127456d+00
 2380 continue
      return
      end
      subroutine t2n21g(e,cs,cp,cd,ia,ngauss)
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),cs(*),cp(*),cd(*)
      iaa=ia-36
      go to (10,20,100,200,300,400,500,600,700,800,900,1000),iaa
  10  go to ( 11, 11, 12, 11, 11, 11 ),ngauss
  11  call caserr2('requested basis set not available for rubidium')
  12  e( 1)=6816.7225d+00
      e( 2)=1033.0007d+00
      e( 3)=226.90861d+00
      e( 4)=304.1283d+00
      e( 5)=66.26058d+00
      e( 6)=20.91945d+00
      e( 7)=22.46533d+00
      e( 8)=7.877468d+00
      e( 9)=2.705271d+00
      e(10)=2.692116d+00
      e(11)=0.7230563d+00
      e(12)=0.2598383d+00
      e(13)=1.897140d-01
      e(14)=3.399726d-02
      e(15)=1.471231d-02
      e(16)=40.866031d+00
      e(17)=10.840885d+00
      e(18)=3.0508341d+00
      cs( 1)=0.0624962d+00
      cs( 2)=0.3719500d+00
      cs( 3)=0.6857293d+00
      cs( 4)=-0.1123296d+00
      cs( 5)=0.09075080d+00
      cs( 6)=0.9639410d+00
      cs( 7)=-0.3004850d+00
      cs( 8)=0.2783566d+00
      cs( 9)=0.9076099d+00
      cs(10)=-0.3311623d+00
      cs(11)=0.5096991d+00
      cs(12)=0.6982461d+00
      cs(13)=-2.711927d-01
      cs(14)=1.141550d+00
      cs(15)=1.0d+00
      cp( 1)=0.1484409d+00
      cp( 2)=0.5891247d+00
      cp( 3)=0.4258251d+00
      cp( 4)=0.02445410d+00
      cp( 5)=0.4944538d+00
      cp( 6)=0.5718567d+00
      cp( 7)=0.01190051d+00
      cp( 8)=0.4951731d+00
      cp( 9)=0.5737244d+00
      cp(10)=3.081009d-04
      cp(11)=9.998654d-01
      cp(12)=1.0d+00
      cd( 1)=0.1466037d+00
      cd( 2)=0.5127252d+00
      cd( 3)=0.5699804d+00
      go to 2380
  20  go to ( 21, 21, 22, 21, 21, 21 ),ngauss
  21  call caserr2('requested basis set not available for strontium')
  22  e( 1)=7215.4735d+00
      e( 2)=1092.8519d+00
      e( 3)=239.98182d+00
      e( 4)=322.1246d+00
      e( 5)=70.09046d+00
      e( 6)=22.17641d+00
      e( 7)=23.92763d+00
      e( 8)=8.475114d+00
      e( 9)=2.942934d+00
      e(10)=2.940966d+00
      e(11)=0.8523559d+00
      e(12)=0.3215375d+00
      e(13)=3.480419d-01
      e(14)=4.817715d-02
      e(15)=2.180335d-02
      e(16)=44.566115d+00
      e(17)=11.881489d+00
      e(18)=3.3875579d+00
      cs( 1)=0.0622818d+00
      cs( 2)=0.3713101d+00
      cs( 3)=0.6864439d+00
      cs( 4)=-0.1122353d+00
      cs( 5)=0.08954360d+00
      cs( 6)=0.9648135d+00
      cs( 7)=-0.3024718d+00
      cs( 8)=0.2700841d+00
      cs( 9)=0.9159200d+00
      cs(10)=-0.3519846d+00
      cs(11)=0.4972546d+00
      cs(12)=0.7238598d+00
      cs(13)=-2.851468d-01
      cs(14)=1.120939d+00
      cs(15)=1.0d+00
      cp( 1)=0.1488369d+00
      cp( 2)=0.5919466d+00
      cp( 3)=0.4221715d+00
      cp( 4)=0.02483697d+00
      cp( 5)=0.4934775d+00
      cp( 6)=0.5709829d+00
      cp( 7)=0.009723336d+00
      cp( 8)=0.4983220d+00
      cp( 9)=0.5650561d+00
      cp(10)=3.081096d-04
      cp(11)=9.998935d-01
      cp(12)=1.0d+00
      cd( 1)=0.1451271d+00
      cd( 2)=0.5130677d+00
      cd( 3)=0.5676640d+00
      go to 2380
  100 go to ( 110, 110, 120, 110, 110, 110),ngauss
 110  call caserr2('requested basis set not available for yttrium')
 120  e( 1) = 7.646421d+03
      e( 2) = 1.156863d+03
      e( 3) = 2.537152d+02
      e( 4) = 3.418540d+02
      e( 5) = 7.420986d+01
      e( 6) = 2.351352d+01
      e( 7) = 1.886260d+01
      e( 8) = 1.645405d+01
      e( 9) = 3.484500d+00
      e(10) = 3.221733d+00
      e(11) = 1.050705d+00
      e(12) = 3.925923d-01
      e(13) = 4.327637d-01
      e(14) = 5.701219d-02
      e(15) = 2.375370d-02
      e(16) = 5.035375d+01
      e(17) = 1.353078d+01
      e(18) = 3.944996d+00
      e(19) = 1.530137d+00
      e(20) = 6.300673d-01
      e(21) = 2.165884d-01
      cs( 1) = 6.189050d-02
      cs( 2) = 3.702068d-01
      cs( 3) = 6.877560d-01
      cs( 4) = -1.119001d-01
      cs( 5) = 8.680524d-02
      cs( 6) = 9.667847d-01
      cs( 7) = -1.477873d+00
      cs( 8) =  1.347259d+00
      cs( 9) =  1.006231d+00
      cs(10) = -3.699580d-01
      cs(11) = 4.308639d-01
      cs(12) = 8.020874d-01
      cs(13) = -3.464582d-01
      cs(14) =  1.132777d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.485717d-01
      cp( 2) = 5.943067d-01
      cp( 3) = 4.196040d-01
      cp( 4) = -7.041407d-01
      cp( 5) = 1.057862d+00
      cp( 6) = 7.393821d-01
      cp( 7) = 2.494435d-03
      cp( 8) = 4.537623d-01
      cp( 9) = 6.130679d-01
      cp(10) =-1.336559d-03
      cp(11) = 1.000440d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.367956d-01
      cd( 2) = 5.019062d-01
      cd( 3) = 5.788598d-01
      cd( 4) = 3.384033d-01
      cd( 5) = 7.293289d-01
      cd( 6) = 1.000000d+00
      go to 2380
  200 go to ( 210, 210, 220, 210, 210, 210),ngauss
 210  call caserr2('requested basis set not available for zirconium')
 220  e( 1) = 8.084592d+03
      e( 2) = 1.221668d+03
      e( 3) = 2.676917d+02
      e( 4) = 3.610212d+02
      e( 5) = 7.830495d+01
      e( 6) = 2.484523d+01
      e( 7) = 2.000628d+01
      e( 8) = 1.757415d+01
      e( 9) = 3.742985d+00
      e(10) = 3.554788d+00
      e(11) = 1.178992d+00
      e(12) = 4.446966d-01
      e(13) = 5.050488d-01
      e(14) = 6.211612d-02
      e(15) = 2.557955d-02
      e(16) = 5.472323d+01
      e(17) = 1.477416d+01
      e(18) = 4.358961d+00
      e(19) = 1.862842d+00
      e(20) = 6.433135d-01
      e(21) = 1.993954d-01
      cs( 1) = 6.157760d-02
      cs( 2) = 3.693989d-01
      cs( 3) = 6.877280d-01
      cs( 4) = -1.119067d-01
      cs( 5) = 8.583993d-02
      cs( 6) = 9.675130d-01
      cs( 7) = -1.544349d+00
      cs( 8) =  1.409596d+00
      cs( 9) =  1.009349d+00
      cs(10) = -3.793873d-01
      cs(11) = 4.232847d-01
      cs(12) = 8.140676d-01
      cs(13) = -3.295118d-01
      cs(14) =  1.120709d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.487616d-01
      cp( 2) = 5.965690d-01
      cp( 3) = 4.167849d-01
      cp( 4) = -7.568480d-01
      cp( 5) = 1.112090d+00
      cp( 6) = 7.368024d-01
      cp( 7) = 2.599455d-03
      cp( 8) = 4.599758d-01
      cp( 9) = 6.058485d-01
      cp(10) =-1.248930d-03
      cp(11) = 1.000384d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.348240d-01
      cd( 2) = 5.005544d-01
      cd( 3) = 5.787823d-01
      cd( 4) = 2.850320d-01
      cd( 5) = 7.792074d-01
      cd( 6) = 1.000000d+00
      go to 2380
  300 go to ( 310, 310, 320, 310, 310, 310),ngauss
  310 call caserr2('requested basis set not available for niobium')
 320  e( 1) = 8.466517d+03
      e( 2) = 1.281261d+03
      e( 3) = 2.812311d+02
      e( 4) = 3.794729d+02
      e( 5) = 8.233589d+01
      e( 6) = 2.622248d+01
      e( 7) = 2.116294d+01
      e( 8) = 1.858978d+01
      e( 9) = 4.009981d+00
      e(10) = 3.836375d+00
      e(11) = 1.303325d+00
      e(12) = 4.934306d-01
      e(13) = 5.723734d-01
      e(14) = 6.820320d-02
      e(15) = 2.715715d-02
      e(16) = 5.901219d+01
      e(17) = 1.601279d+01
      e(18) = 4.777185d+00
      e(19) = 1.970443d+00
      e(20) = 6.619347d-01
      e(21) = 2.059972d-01
      cs( 1) = 6.180380d-02
      cs( 2) = 3.698049d-01
      cs( 3) = 6.880795d-01
      cs( 4) = -1.121063d-01
      cs( 5) = 8.650279d-02
      cs( 6) = 9.670574d-01
      cs( 7) = -1.555131d+00
      cs( 8) =  1.417939d+00
      cs( 9) =  1.010620d+00
      cs(10) = -3.891037d-01
      cs(11) = 4.349696d-01
      cs(12) = 8.115899d-01
      cs(13) = -3.156094d-01
      cs(14) =  1.114047d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.496674d-01
      cp( 2) = 5.987182d-01
      cp( 3) = 4.132382d-01
      cp( 4) = -7.554622d-01
      cp( 5) = 1.113966d+00
      cp( 6) = 7.327672d-01
      cp( 7) = 3.290347d-03
      cp( 8) = 4.716451d-01
      cp( 9) = 5.936990d-01
      cp(10) =-1.133018d-03
      cp(11) = 1.000338d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.337104d-01
      cd( 2) = 5.000390d-01
      cd( 3) = 5.775384d-01
      cd( 4) = 3.106809d-01
      cd( 5) = 7.800691d-01
      cd( 6) = 1.000000d+00
      go to 2380
  400 go to ( 410, 410, 420, 410, 410, 410),ngauss
  410 call caserr2('requested basis set not available for molybdenum')
 420  e( 1) = 8.899491d+03
      e( 2) = 1.346764d+03
      e( 3) = 2.956352d+02
      e( 4) = 3.993139d+02
      e( 5) = 8.659356d+01
      e( 6) = 2.763904d+01
      e( 7) = 2.250292d+01
      e( 8) = 1.949171d+01
      e( 9) = 4.278180d+00
      e(10) = 4.163021d+00
      e(11) = 1.435305d+00
      e(12) = 5.437821d-01
      e(13) = 6.318014d-01
      e(14) = 7.325791d-02
      e(15) = 2.802515d-02
      e(16) = 6.378045d+01
      e(17) = 1.737358d+01
      e(18) = 5.230784d+00
      e(19) = 2.270937d+00
      e(20) = 7.546530d-01
      e(21) = 2.351422d-01
      cs( 1) = 6.170640d-02
      cs( 2) = 3.694536d-01
      cs( 3) = 6.884343d-01
      cs( 4) = -1.121440d-01
      cs( 5) = 8.601148d-02
      cs( 6) = 9.674335d-01
      cs( 7) = -1.422306d+00
      cs( 8) =  1.284135d+00
      cs( 9) =  1.010866d+00
      cs(10) = -3.964233d-01
      cs(11) = 4.370792d-01
      cs(12) = 8.148462d-01
      cs(13) = -3.033617d-01
      cs(14) =  1.108413d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.500665d-01
      cp( 2) = 6.007695d-01
      cp( 3) = 4.103865d-01
      cp( 4) = -6.680661d-01
      cp( 5) = 1.030346d+00
      cp( 6) = 7.283480d-01
      cp( 7) = 2.962629d-03
      cp( 8) = 4.791471d-01
      cp( 9) = 5.864869d-01
      cp(10) =-1.079133d-03
      cp(11) = 1.000313d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.317388d-01
      cd( 2) = 4.985316d-01
      cd( 3) = 5.781775d-01
      cd( 4) = 3.112644d-01
      cd( 5) = 7.810342d-01
      cd( 6) = 1.000000d+00
      go to 2380
  500 go to ( 510, 510, 520, 510, 510, 510),ngauss
  510 call caserr2('requested basis set not available for technecium')
 520  e( 1) = 9.329482d+03
      e( 2) = 1.412506d+03
      e( 3) = 3.102643d+02
      e( 4) = 4.188175d+02
      e( 5) = 9.125078d+01
      e( 6) = 2.911212d+01
      e( 7) = 2.591064d+01
      e( 8) = 2.326770d+01
      e( 9) = 4.707083d+00
      e(10) = 4.441138d+00
      e(11) = 1.595639d+00
      e(12) = 5.955598d-01
      e(13) = 6.738812d-01
      e(14) = 7.724070d-02
      e(15) = 2.869556d-02
      e(16) = 6.878375d+01
      e(17) = 1.880389d+01
      e(18) = 5.705228d+00
      e(19) = 2.599164d+00
      e(20) = 8.622757d-01
      e(21) = 2.706073d-01
      cs( 1) = 6.171190d-02
      cs( 2) = 3.693370d-01
      cs( 3) = 6.884724d-01
      cs( 4) = -1.124025d-01
      cs( 5) = 8.531816d-02
      cs( 6) = 9.681774d-01
      cs( 7) = -1.380446d+00
      cs( 8) =  1.197900d+00
      cs( 9) =  1.052649d+00
      cs(10) = -4.041144d-01
      cs(11) = 4.398379d-01
      cs(12) = 8.219360d-01
      cs(13) = -2.700028d-01
      cs(14) =  1.099150d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.500719d-01
      cp( 2) = 6.000567d-01
      cp( 3) = 4.109856d-01
      cp( 4) = -1.655271d+00
      cp( 5) = 1.986020d+00
      cp( 6) = 7.290339d-01
      cp( 7) = 1.229062d-02
      cp( 8) = 4.632067d-01
      cp( 9) = 5.983826d-01
      cp(10) =-9.197676d-04
      cp(11) = 1.000264d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.296930d-01
      cd( 2) = 4.966191d-01
      cd( 3) = 5.795466d-01
      cd( 4) = 3.092195d-01
      cd( 5) = 7.829056d-01
      cd( 6) = 1.000000d+00
      go to 2380
  600 go to ( 610, 610, 620, 610, 610, 610),ngauss
 610  call caserr2('requested basis set not available for ruthenium')
 620  e( 1) = 9.786161d+03
      e( 2) = 1.481477d+03
      e( 3) = 3.254122d+02
      e( 4) = 4.398665d+02
      e( 5) = 9.576273d+01
      e( 6) = 3.060566d+01
      e( 7) = 2.727737d+01
      e( 8) = 2.451082d+01
      e( 9) = 5.008946d+00
      e(10) = 4.765812d+00
      e(11) = 1.734531d+00
      e(12) = 6.466355d-01
      e(13) = 7.406620d-01
      e(14) = 8.217096d-02
      e(15) = 3.009659d-02
      e(16) = 7.398330d+01
      e(17) = 2.028149d+01
      e(18) = 6.194298d+00
      e(19) = 2.889108d+00
      e(20) = 9.539610d-01
      e(21) = 2.958807d-01
      cs( 1) = 6.160520d-02
      cs( 2) = 3.689816d-01
      cs( 3) = 6.888451d-01
      cs( 4) = -1.123912d-01
      cs( 5) = 8.469449d-02
      cs( 6) = 9.686379d-01
      cs( 7) = -1.395553d+00
      cs( 8) =  1.210852d+00
      cs( 9) =  1.054045d+00
      cs(10) = -4.103627d-01
      cs(11) = 4.480025d-01
      cs(12) = 8.198083d-01
      cs(13) = -2.639655d-01
      cs(14) =  1.094857d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.503791d-01
      cp( 2) = 6.019294d-01
      cp( 3) = 4.084639d-01
      cp( 4) = -1.668618d+00
      cp( 5) = 2.002799d+00
      cp( 6) = 7.251435d-01
      cp( 7) = 1.127417d-02
      cp( 8) = 4.727032d-01
      cp( 9) = 5.898430d-01
      cp(10) =-7.620441d-04
      cp(11) = 1.000212d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.277598d-01
      cd( 2) = 4.951474d-01
      cd( 3) = 5.806552d-01
      cd( 4) = 3.159938d-01
      cd( 5) = 7.780656d-01
      cd( 6) = 1.000000d+00
      go to 2380
  700 go to ( 710 , 710 , 720 , 710 , 710 , 710 ),ngauss
  710 call caserr2('requested basis set not available for rhodium')
 720  e( 1) = 1.021771d+04
      e( 2) = 1.548412d+03
      e( 3) = 3.404990d+02
      e( 4) = 4.607593d+02
      e( 5) = 1.003289d+02
      e( 6) = 3.213971d+01
      e( 7) = 2.879329d+01
      e( 8) = 2.591768d+01
      e( 9) = 5.320640d+00
      e(10) = 5.109748d+00
      e(11) = 1.875414d+00
      e(12) = 6.995578d-01
      e(13) = 8.005711d-01
      e(14) = 8.732134d-02
      e(15) = 3.140693d-02
      e(16) = 7.925597d+01
      e(17) = 2.178945d+01
      e(18) = 6.697518d+00
      e(19) = 3.190908d+00
      e(20) = 1.054575d+00
      e(21) = 3.260791d-01
      cs( 1) = 6.173240d-02
      cs( 2) = 3.691533d-01
      cs( 3) = 6.885138d-01
      cs( 4) = -1.124461d-01
      cs( 5) = 8.438113d-02
      cs( 6) = 9.689016d-01
      cs( 7) = -1.404091d+00
      cs( 8) =  1.216169d+00
      cs( 9) =  1.056555d+00
      cs(10) = -4.126469d-01
      cs(11) = 4.518853d-01
      cs(12) = 8.188989d-01
      cs(13) = -2.553480d-01
      cs(14) =  1.091308d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.508582d-01
      cp( 2) = 6.035139d-01
      cp( 3) = 4.060249d-01
      cp( 4) = -1.712218d+00
      cp( 5) = 2.047603d+00
      cp( 6) = 7.229831d-01
      cp( 7) = 9.374064d-03
      cp( 8) = 4.815544d-01
      cp( 9) = 5.822312d-01
      cp(10) =-7.759711d-04
      cp(11) = 1.000212d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.261896d-01
      cd( 2) = 4.939541d-01
      cd( 3) = 5.813296d-01
      cd( 4) = 3.210399d-01
      cd( 5) = 7.738518d-01
      cd( 6) = 1.000000d+00
      go to 2380
 800  go to (810,810, 820, 810,810, 810),ngauss
 810  call caserr2('requested basis set not available for palladium')
 820  e( 1) = 1.072874d+04
      e( 2) = 1.624074d+03
      e( 3) = 3.567937d+02
      e( 4) = 4.824783d+02
      e( 5) = 1.050590d+02
      e( 6) = 3.368145d+01
      e( 7) = 3.018654d+01
      e( 8) = 2.716642d+01
      e( 9) = 5.635934d+00
      e(10) = 5.475374d+00
      e(11) = 1.997604d+00
      e(12) = 7.439302d-01
      e(13) = 8.901632d-01
      e(14) = 9.282090d-02
      e(15) = 3.377394d-02
      e(16) = 8.423691d+01
      e(17) = 2.324919d+01
      e(18) = 7.196760d+00
      e(19) = 3.473077d+00
      e(20) = 1.148050d+00
      e(21) = 3.548106d-01
      cs( 1) = 6.142950d-02
      cs( 2) = 3.683282d-01
      cs( 3) = 6.895025d-01
      cs( 4) = -1.126789d-01
      cs( 5) = 8.461197d-02
      cs( 6) = 9.687868d-01
      cs( 7) = -1.418547d+00
      cs( 8) =  1.229444d+00
      cs( 9) =  1.057083d+00
      cs(10) = -4.172610d-01
      cs(11) = 4.705878d-01
      cs(12) = 8.046363d-01
      cs(13) = -2.784324d-01
      cs(14) =  1.093028d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.510793d-01
      cp( 2) = 6.050916d-01
      cp( 3) = 4.039804d-01
      cp( 4) = -1.709817d+00
      cp( 5) = 2.049308d+00
      cp( 6) = 7.186296d-01
      cp( 7) = 1.158391d-02
      cp( 8) = 4.974548d-01
      cp( 9) = 5.655272d-01
      cp(10) =-1.271580d-03
      cp(11) = 1.000332d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.256429d-01
      cd( 2) = 4.937201d-01
      cd( 3) = 5.803431d-01
      cd( 4) = 3.281543d-01
      cd( 5) = 7.680267d-01
      cd( 6) = 1.000000d+00
      go to 2380
 900  go to (910,910,920,910,910,910),ngauss
 910  call caserr2('requested basis set not available for silver')
 920  e( 1) = 1.119078d+04
      e( 2) = 1.695077d+03
      e( 3) = 3.726752d+02
      e( 4) = 5.046162d+02
      e( 5) = 1.098718d+02
      e( 6) = 3.529513d+01
      e( 7) = 3.156877d+01
      e( 8) = 2.834397d+01
      e( 9) = 5.945127d+00
      e(10) = 5.800256d+00
      e(11) = 2.127256d+00
      e(12) = 7.935512d-01
      e(13) = 9.285445d-01
      e(14) = 9.725467d-02
      e(15) = 3.493292d-02
      e(16) = 8.993335d+01
      e(17) = 2.487496d+01
      e(18) = 7.738191d+00
      e(19) = 3.796557d+00
      e(20) = 1.256644d+00
      e(21) = 3.881333d-01
      cs( 1) = 6.149480d-02
      cs( 2) = 3.684053d-01
      cs( 3) = 6.893247d-01
      cs( 4) = -1.126577d-01
      cs( 5) = 8.402784d-02
      cs( 6) = 9.692344d-01
      cs( 7) = -1.422028d+00
      cs( 8) =  1.234098d+00
      cs( 9) =  1.055683d+00
      cs(10) = -4.196171d-01
      cs(11) = 4.843501d-01
      cs(12) = 7.952035d-01
      cs(13) = -2.523005d-01
      cs(14) =  1.087392d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.514798d-01
      cp( 2) = 6.065142d-01
      cp( 3) = 4.018302d-01
      cp( 4) = -1.673366d+00
      cp( 5) = 2.018976d+00
      cp( 6) = 7.126889d-01
      cp( 7) = 1.430410d-02
      cp( 8) = 5.071943d-01
      cp( 9) = 5.539736d-01
      cp(10) =-1.480711d-03
      cp(11) = 1.000388d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.240159d-01
      cd( 2) = 4.923831d-01
      cd( 3) = 5.814968d-01
      cd( 4) = 3.314259d-01
      cd( 5) = 7.651634d-01
      cd( 6) = 1.000000d+00
      go to 2380
 1000 go to (1010,1010,1020,1010,1010,1010),ngauss
 1010 call caserr2('requested basis set not available for cadmium')
 1020 e( 1) = 1.168609d+04
      e( 2) = 1.770111d+03
      e( 3) = 3.892090d+02
      e( 4) = 5.276004d+02
      e( 5) = 1.148329d+02
      e( 6) = 3.695829d+01
      e( 7) = 3.301548d+01
      e( 8) = 2.954543d+01
      e( 9) = 6.278508d+00
      e(10) = 6.150596d+00
      e(11) = 2.259746d+00
      e(12) = 8.414261d-01
      e(13) = 9.490686d-01
      e(14) = 1.014878d-01
      e(15) = 3.598726d-02
      e(16) = 9.547274d+01
      e(17) = 2.648196d+01
      e(18) = 8.282886d+00
      e(19) = 4.082141d+00
      e(20) = 1.357279d+00
      e(21) = 4.208308d-01
      cs( 1) = 6.142650d-02
      cs( 2) = 3.681567d-01
      cs( 3) = 6.895722d-01
      cs( 4) = -1.125925d-01
      cs( 5) = 8.326963d-02
      cs( 6) = 9.697978d-01
      cs( 7) = -1.406471d+00
      cs( 8) =  1.218156d+00
      cs( 9) =  1.055520d+00
      cs(10) = -4.229209d-01
      cs(11) = 4.987714d-01
      cs(12) = 7.850775d-01
      cs(13) = -2.215547d-01
      cs(14) =  1.080944d+00
      cs(15) =  1.000000d+00
      cp( 1) = 1.518051d-01
      cp( 2) = 6.077598d-01
      cp( 3) = 3.999632d-01
      cp( 4) = -1.609024d+00
      cp( 5) = 1.959568d+00
      cp( 6) = 7.080271d-01
      cp( 7) = 1.448229d-02
      cp( 8) = 5.186611d-01
      cp( 9) = 5.426658d-01
      cp(10) =-1.540266d-03
      cp(11) = 1.000412d+00
      cp(12) = 1.000000d+00
      cd( 1) = 1.230828d-01
      cd( 2) = 4.916768d-01
      cd( 3) = 5.815408d-01
      cd( 4) = 3.379410d-01
      cd( 5) = 7.591679d-01
      cd( 6) = 1.000000d+00
 2380 continue
      return
      end
      subroutine dn21g(e,cs,cp,cd,ia,ngauss)
c
c
c***********************************************************************
c     routine to contain the data for the n-21g and nm-21g bases.
c     consult routine n21g for calling details.
c     n-21g may be augmented with both polarisation and
c     diffuse functions (see routine n31).
c***********************************************************************
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c     branch on input atomic number.
c
      dimension e(*),cs(*),cp(*),cd(*)
      go to (100,120,140,280,420,560,700,840,980,1120,1260,1400,1540,
     +     1680,1820,1960,2100,2240),ia
c
c***********************************************************************
c          hydrogen
c***********************************************************************
c
  100 continue
      e( 1) = 0.450180d+01
      e( 2) = 0.681444d+00
      e( 3) = 0.151398d+00
      e( 4) = 0.0360d+00
      e( 5) = 1.1d0
      cs( 1) = 0.156285d+00
      cs( 2) = 0.904691d+00
      cs( 3) = 0.100000d+01
      cs( 4) = 0.100000d+01
      cp( 1) = 0.100000d+01
      go to 2380
  120 continue
c
c***********************************************************************
c          helium
c***********************************************************************
c
      e( 1) = 0.136267d+02
      e( 2) = 0.199935d+01
      e( 3) = 0.382993d+00
      e( 4) = 0.0360d+00
      e( 5) = 1.1d0
      cs( 1) = 0.175230d+00
      cs( 2) = 0.893483d+00
      cs( 3) = 0.100000d+01
      cs( 4) = 0.100000d+01
      cp( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          lithium
c**********************************************************************
c
  140 go to ( 160, 180, 220, 240, 200, 260),ngauss
  160 continue
  180 continue
  200 continue
      call caserr2('requested basis set not available for lithium')
  220 continue
      e( 1) = 0.368382d+02
      e( 2) = 0.548172d+01
      e( 3) = 0.111327d+01
      e( 4) = 0.540205d+00
      e( 5) = 0.102255d+00
      e( 6) = 0.285645d-01
      e( 7) = 0.0074d+00
      e( 8) = 0.20d0
      cs( 1) = 0.696686d-01
      cs( 2) = 0.381346d+00
      cs( 3) = 0.681702d+00
      cs( 4) = -0.263127d+00
      cs( 5) = 0.114339d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 270
  240 continue
      e( 1) = 0.109353d+03
      e( 2) = 0.164228d+02
      e( 3) = 0.359415d+01
      e( 4) = 0.905297d+00
      e( 5) = 0.540205d+00
      e( 6) = 0.102255d+00
      e( 7) = 0.285645d-01
      e( 8) = 0.0074d+00
      e( 9) = 0.20d0
      cs( 1) = 0.190277d-01
      cs( 2) = 0.130276d+00
      cs( 3) = 0.439082d+00
      cs( 4) = 0.557314d+00
      cs( 5) = -0.263127d+00
      cs( 6) = 0.114339d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 270
  260 continue
      e( 1) = 0.642418d+03
      e( 2) = 0.965164d+02
      e( 3) = 0.220174d+02
      e( 4) = 0.617645d+01
      e( 5) = 0.193511d+01
      e( 6) = 0.639577d+00
      e( 7) = 0.540205d+00
      e( 8) = 0.102255d+00
      e( 9) = 0.285645d-01
      e(10) = 0.0074d+00
      e(11) = 0.20d0
      cs( 1) = 0.215096d-02
      cs( 2) = 0.162677d-01
      cs( 3) = 0.776383d-01
      cs( 4) = 0.246495d+00
      cs( 5) = 0.467506d+00
      cs( 6) = 0.346915d+00
      cs( 7) = -0.263127d+00
      cs( 8) = 0.114339d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  270 continue
      cp( 1) = 0.161546d+00
      cp( 2) = 0.915663d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          beryllium
c**********************************************************************
c
  280 go to ( 300, 320, 360, 380, 340, 400),ngauss
  300 continue
  320 continue
  340 continue
      call caserr2('requested basis set not available for beryllium')
  360 continue
      e( 1) = 0.718876d+02
      e( 2) = 0.107289d+02
      e( 3) = 0.222205d+01
      e( 4) = 0.129548d+01
      e( 5) = 0.268881d+00
      e( 6) = 0.773501d-01
      e( 7) = 0.0207d0
      e( 8) = 0.40d0
      cs( 1) = 0.644263d-01
      cs( 2) = 0.366096d+00
      cs( 3) = 0.695934d+00
      cs( 4) = -0.421064d+00
      cs( 5) = 0.122407d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 410
  380 continue
      e( 1) = 0.207980d+03
      e( 2) = 0.314316d+02
      e( 3) = 0.699419d+01
      e( 4) = 0.181295d+01
      e( 5) = 0.129548d+01
      e( 6) = 0.268881d+00
      e( 7) = 0.773501d-01
      e( 8) = 0.0207d0
      e( 9) = 0.40d0
      cs( 1) = 0.180520d-01
      cs( 2) = 0.123804d+00
      cs( 3) = 0.427580d+00
      cs( 4) = 0.569901d+00
      cs( 5) = -0.421064d+00
      cs( 6) = 0.122407d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 410
  400 continue
      e( 1) = 0.126450d+04
      e( 2) = 0.189930d+03
      e( 3) = 0.431275d+02
      e( 4) = 0.120889d+02
      e( 5) = 0.380790d+01
      e( 6) = 0.128266d+01
      e( 7) = 0.129548d+01
      e( 8) = 0.268881d+00
      e( 9) = 0.773501d-01
      e(10) = 0.0207d0
      e(11) = 0.40d0
      cs( 1) = 0.194336d-02
      cs( 2) = 0.148251d-01
      cs( 3) = 0.720662d-01
      cs( 4) = 0.237022d+00
      cs( 5) = 0.468789d+00
      cs( 6) = 0.356382d+00
      cs( 7) = -0.421064d+00
      cs( 8) = 0.122407d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  410 continue
      cp( 1) = 0.205132d+00
      cp( 2) = 0.882528d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          boron
c**********************************************************************
c
  420 go to ( 440, 460, 500, 520, 480, 540),ngauss
  440 continue
  460 continue
  480 continue
      call caserr2('requested basis set not available for boron')
  500 continue
      e( 1) = 0.116434d+03
      e( 2) = 0.174314d+02
      e( 3) = 0.368016d+01
      e( 4) = 0.228187d+01
      e( 5) = 0.465248d+00
      e( 6) = 0.124328d+00
      e( 7) = 0.0315d0
      e( 8) = 0.60d0
      cs( 1) = 0.629605d-01
      cs( 2) = 0.363304d+00
      cs( 3) = 0.697255d+00
      cs( 4) = -0.368662d+00
      cs( 5) = 0.119944d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 550
  520 continue
      e( 1) = 0.345445d+03
      e( 2) = 0.519156d+02
      e( 3) = 0.115337d+02
      e( 4) = 0.301981d+01
      e( 5) = 0.228187d+01
      e( 6) = 0.465248d+00
      e( 7) = 0.124328d+00
      e( 8) = 0.0315d0
      e( 9) = 0.60d0
      cs( 1) = 0.170348d-01
      cs( 2) = 0.119266d+00
      cs( 3) = 0.424702d+00
      cs( 4) = 0.575037d+00
      cs( 5) = -0.368662d+00
      cs( 6) = 0.119944d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 550
  540 continue
      e( 1) = 0.208212d+04
      e( 2) = 0.312310d+03
      e( 3) = 0.708874d+02
      e( 4) = 0.198525d+02
      e( 5) = 0.629161d+01
      e( 6) = 0.212862d+01
      e( 7) = 0.228187d+01
      e( 8) = 0.465248d+00
      e( 9) = 0.124328d+00
      e(10) = 0.0315d0
      e(11) = 0.60d0
      cs( 1) = 0.184986d-02
      cs( 2) = 0.141277d-01
      cs( 3) = 0.692697d-01
      cs( 4) = 0.232393d+00
      cs( 5) = 0.470154d+00
      cs( 6) = 0.360288d+00
      cs( 7) = -0.368662d+00
      cs( 8) = 0.119944d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  550 continue
      cp( 1) = 0.231152d+00
      cp( 2) = 0.866764d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          carbon
c**********************************************************************
c
  560 go to ( 580, 600, 640, 660, 620, 680),ngauss
  580 continue
  600 continue
  620 continue
      call caserr2('requested basis set not available for carbon')
  640 continue
      e( 1) = 0.172256d+03
      e( 2) = 0.259109d+02
      e( 3) = 0.553335d+01
      e( 4) = 0.366498d+01
      e( 5) = 0.770545d+00
      e( 6) = 0.195857d+00
      e( 7) = 0.0438d0
      e( 8) = 0.80d0
      cs( 1) = 0.617669d-01
      cs( 2) = 0.358794d+00
      cs( 3) = 0.700713d+00
      cs( 4) = -0.395897d+00
      cs( 5) = 0.121584d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 690
  660 continue
      e( 1) = 0.514836d+03
      e( 2) = 0.773470d+02
      e( 3) = 0.172534d+02
      e( 4) = 0.455754d+01
      e( 5) = 0.366498d+01
      e( 6) = 0.770545d+00
      e( 7) = 0.195857d+00
      e( 8) = 0.0438d0
      e( 9) = 0.80d0
      cs( 1) = 0.165399d-01
      cs( 2) = 0.116447d+00
      cs( 3) = 0.419945d+00
      cs( 4) = 0.580709d+00
      cs( 5) = -0.395897d+00
      cs( 6) = 0.121584d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 690
  680 continue
      e( 1) = 0.304752d+04
      e( 2) = 0.456424d+03
      e( 3) = 0.103653d+03
      e( 4) = 0.292258d+02
      e( 5) = 0.934863d+01
      e( 6) = 0.318904d+01
      e( 7) = 0.366498d+01
      e( 8) = 0.770545d+00
      e( 9) = 0.195857d+00
      e(10) = 0.0438d0
      e(11) = 0.80d0
      cs( 1) = 0.182588d-02
      cs( 2) = 0.140566d-01
      cs( 3) = 0.687570d-01
      cs( 4) = 0.230422d+00
      cs( 5) = 0.468463d+00
      cs( 6) = 0.362780d+00
      cs( 7) = -0.395897d+00
      cs( 8) = 0.121584d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  690 continue
      cp( 1) = 0.236460d+00
      cp( 2) = 0.860619d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          nitrogen
c**********************************************************************
c
  700 go to ( 720, 740, 780, 800, 760, 820),ngauss
  720 continue
  740 continue
  760 continue
      call caserr2('requested basis set not available for nitrogen')
  780 continue
      e( 1) = 0.242766d+03
      e( 2) = 0.364851d+02
      e( 3) = 0.781449d+01
      e( 4) = 0.542522d+01
      e( 5) = 0.114915d+01
      e( 6) = 0.283205d+00
      e( 7) = 0.0639d0
      e( 8) = 0.80d0
      cs( 1) = 0.598657d-01
      cs( 2) = 0.352955d+00
      cs( 3) = 0.706513d+00
      cs( 4) = -0.413301d+00
      cs( 5) = 0.122442d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 790
  800 continue
      e( 1) = 0.715345d+03
      e( 2) = 0.107490d+03
      e( 3) = 0.240414d+02
      e( 4) = 0.639437d+01
      e( 5) = 0.542522d+01
      e( 6) = 0.114915d+01
      e( 7) = 0.283205d+00
      e( 8) = 0.0639d0
      e( 9) = 0.80d0
      cs( 1) = 0.162587d-01
      cs( 2) = 0.114910d+00
      cs( 3) = 0.417387d+00
      cs( 4) = 0.583539d+00
      cs( 5) = -0.413301d+00
      cs( 6) = 0.122442d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 790
  820 continue
      e( 1) = 0.415011d+04
      e( 2) = 0.620084d+03
      e( 3) = 0.141688d+03
      e( 4) = 0.403367d+02
      e( 5) = 0.130267d+02
      e( 6) = 0.447003d+01
      e( 7) = 0.542522d+01
      e( 8) = 0.114915d+01
      e( 9) = 0.283205d+00
      e(10) = 0.0639d0
      e(11) = 0.80d0
      cs( 1) = 0.184541d-02
      cs( 2) = 0.141645d-01
      cs( 3) = 0.686325d-01
      cs( 4) = 0.228574d+00
      cs( 5) = 0.466162d+00
      cs( 6) = 0.365672d+00
      cs( 7) = -0.413301d+00
      cs( 8) = 0.122442d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  790 continue
      cp( 1) = 0.237972d+00
      cp( 2) = 0.858953d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          oxygen
c**********************************************************************
c
  840 go to ( 860, 880, 920, 940, 900, 960),ngauss
  860 continue
  880 continue
  900 continue
      call caserr2('requested basis set not available for oxygen')
  920 continue
      e( 1) = 0.322037d+03
      e( 2) = 0.484308d+02
      e( 3) = 0.104206d+02
      e( 4) = 0.740294d+01
      e( 5) = 0.157620d+01
      e( 6) = 0.373684d+00
      e( 7) = 0.0845d0
      e( 8) = 0.80d0
      cs( 1) = 0.592394d-01
      cs( 2) = 0.351500d+00
      cs( 3) = 0.707658d+00
      cs( 4) = -0.404453d+00
      cs( 5) = 0.122156d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 970
  940 continue
      e( 1) = 0.938318d+03
      e( 2) = 0.141662d+03
      e( 3) = 0.318308d+02
      e( 4) = 0.851101d+01
      e( 5) = 0.740294d+01
      e( 6) = 0.157620d+01
      e( 7) = 0.373684d+00
      e( 8) = 0.0845d0
      e( 9) = 0.80d0
      cs( 1) = 0.162714d-01
      cs( 2) = 0.114340d+00
      cs( 3) = 0.416787d+00
      cs( 4) = 0.583808d+00
      cs( 5) = -0.404453d+00
      cs( 6) = 0.122156d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 970
  960 continue
      e( 1) = 0.547227d+04
      e( 2) = 0.817806d+03
      e( 3) = 0.186446d+03
      e( 4) = 0.530230d+02
      e( 5) = 0.171800d+02
      e( 6) = 0.591196d+01
      e( 7) = 0.740294d+01
      e( 8) = 0.157620d+01
      e( 9) = 0.373684d+00
      e(10) = 0.0845d0
      e(11) = 0.80d0
      cs( 1) = 0.183217d-02
      cs( 2) = 0.141047d-01
      cs( 3) = 0.686262d-01
      cs( 4) = 0.229376d+00
      cs( 5) = 0.466399d+00
      cs( 6) = 0.364173d+00
      cs( 7) = -0.404453d+00
      cs( 8) = 0.122156d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
  970 continue
      cp( 1) = 0.244586d+00
      cp( 2) = 0.853955d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          fluorine
c**********************************************************************
c
  980 go to ( 1000, 1020, 1060, 1080, 1040, 1100),ngauss
 1000 continue
 1020 continue
 1040 continue
      call caserr2('requested basis set not available for fluorine')
 1060 continue
      e( 1) = 0.413801d+03
      e( 2) = 0.622446d+02
      e( 3) = 0.134340d+02
      e( 4) = 0.977759d+01
      e( 5) = 0.208617d+01
      e( 6) = 0.482383d+00
      e (7) = 0.1076d0
      e( 8) = 0.80d0
      cs( 1) = 0.585483d-01
      cs( 2) = 0.349308d+00
      cs( 3) = 0.709632d+00
      cs( 4) = -0.407327d+00
      cs( 5) = 0.122314d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 1090
 1080 continue
      e( 1) = 0.120491d+04
      e( 2) = 0.181792d+03
      e( 3) = 0.408879d+02
      e( 4) = 0.109630d+02
      e( 5) = 0.977759d+01
      e( 6) = 0.208617d+01
      e( 7) = 0.482383d+00
      e( 8) = 0.1076d0
      e( 9) = 0.80d0
      cs( 1) = 0.160678d-01
      cs( 2) = 0.113306d+00
      cs( 3) = 0.415239d+00
      cs( 4) = 0.585754d+00
      cs( 5) = -0.407327d+00
      cs( 6) = 0.122314d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 1090
 1100 continue
      e( 1) = 0.678319d+04
      e( 2) = 0.104244d+04
      e( 3) = 0.242398d+03
      e( 4) = 0.696320d+02
      e( 5) = 0.226894d+02
      e( 6) = 0.779636d+01
      e( 7) = 0.977759d+01
      e( 8) = 0.208617d+01
      e( 9) = 0.482383d+00
      e(10) = 0.1076d0
      e(11) = 0.80d0
      cs( 1) = 0.188463d-02
      cs( 2) = 0.138121d-01
      cs( 3) = 0.662493d-01
      cs( 4) = 0.221875d+00
      cs( 5) = 0.460842d+00
      cs( 6) = 0.378453d+00
      cs( 7) = -0.407327d+00
      cs( 8) = 0.122314d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
 1090 continue
      cp( 1) = 0.246680d+00
      cp( 2) = 0.852321d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          neon
c**********************************************************************
c
 1120 go to (1140,1160,1200,1220,1180,1240),ngauss
 1140 continue
 1160 continue
 1180 continue
      call caserr2('requested basis set not available for neon')
 1200 continue
      e( 1) = 0.515724d+03
      e( 2) = 0.776538d+02
      e( 3) = 0.168136d+02
      e( 4) = 0.124830d+02
      e( 5) = 0.266451d+01
      e( 6) = 0.606250d+00
      e( 7) = 0.11d0
      e( 8) = 0.80d0
      cs( 1) = 0.581430d-01
      cs( 2) = 0.347951d+00
      cs( 3) = 0.710714d+00
      cs( 4) = -0.409922d+00
      cs( 5) = 0.122431d+01
      cs( 6) = 0.100000d+01
      cs( 7) = 0.100000d+01
      go to 1230
 1220 continue
      e( 1) = 0.151275d+04
      e( 2) = 0.227368d+03
      e( 3) = 0.510739d+02
      e( 4) = 0.137213d+02
      e( 5) = 0.124830d+02
      e( 6) = 0.266451d+01
      e( 7) = 0.606250d+00
      e( 8) = 0.11d0
      e( 9) = 0.80d0
      cs( 1) = 0.158096d-01
      cs( 2) = 0.112430d+00
      cs( 3) = 0.414266d+00
      cs( 4) = 0.587193d+00
      cs( 5) = -0.409922d+00
      cs( 6) = 0.122431d+01
      cs( 7) = 0.100000d+01
      cs( 8) = 0.100000d+01
      go to 1230
 1240 continue
      e( 1) = 0.878583d+04
      e( 2) = 0.132390d+04
      e( 3) = 0.300795d+03
      e( 4) = 0.851891d+02
      e( 5) = 0.276534d+02
      e( 6) = 0.953039d+01
      e( 7) = 0.124830d+02
      e( 8) = 0.266451d+01
      e( 9) = 0.606250d+00
      e(10) = 0.11d0
      e(11) = 0.80d0
      cs( 1) = 0.178077d-02
      cs( 2) = 0.135790d-01
      cs( 3) = 0.670847d-01
      cs( 4) = 0.226825d+00
      cs( 5) = 0.465053d+00
      cs( 6) = 0.368995d+00
      cs( 7) = -0.409922d+00
      cs( 8) = 0.122431d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
 1230 continue
      cp( 1) = 0.247460d+00
      cp( 2) = 0.851743d+00
      cp( 3) = 0.100000d+01
      cp( 4) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          sodium
c***********************************************************************
c
 1260 go to (1280,1300,1360,1320,1340,1380),ngauss
 1280 continue
 1300 continue
 1320 continue
 1340 continue
      call caserr2('requested basis set not available for sodium')
 1360 e( 1) = 0.547613d+03
      e( 2) = 0.820678d+02
      e( 3) = 0.176917d+02
      e( 4) = 0.175407d+02
      e( 5) = 0.379398d+01
      e( 6) = 0.906441d+00
      e( 7) = 0.501824d+00
      e( 8) = 0.609458d-01
      e( 9) = 0.244349d-01
      e(10) = 0.0076d0
      e(11) = 0.175d0
      cs( 1) = 0.674911d-01
      cs( 2) = 0.393505d+00
      cs( 3) = 0.665605d+00
      cs( 4) = -0.111937d+00
      cs( 5) = 0.254654d+00
      cs( 6) = 0.844417d+00
      cs( 7) = -0.219660d+00
      cs( 8) = 0.108912d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 1) = 0.128233d+00
      cp( 2) = 0.471533d+00
      cp( 3) = 0.604273d+00
      cp( 4) = 0.906649d-02
      cp( 5) = 0.997202d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 1380 continue
      e( 1) = 0.999320d+04
      e( 2) = 0.149989d+04
      e( 3) = 0.341951d+03
      e( 4) = 0.946796d+02
      e( 5) = 0.297345d+02
      e( 6) = 0.100063d+02
      e( 7) = 0.150963d+03
      e( 8) = 0.355878d+02
      e( 9) = 0.111683d+02
      e(10) = 0.390201d+01
      e(11) = 0.138177d+01
      e(12) = 0.466382d+00
      e(13) = 0.501824d+00
      e(14) = 0.609458d-01
      e(15) = 0.244349d-01
      e(16) = 0.0076d0
      e(17) = 0.175d0
      cs( 1) = 0.193766d-02
      cs( 2) = 0.148070d-01
      cs( 3) = 0.727055d-01
      cs( 4) = 0.252629d+00
      cs( 5) = 0.493242d+00
      cs( 6) = 0.313169d+00
      cs( 7) = -0.354208d-02
      cs( 8) = -0.439588d-01
      cs( 9) = -0.109752d+00
      cs(10) = 0.187398d+00
      cs(11) = 0.646699d+00
      cs(12) = 0.306058d+00
      cs(13) = -0.219660d+00
      cs(14) = 0.108912d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.500166d-02
      cp( 2) = 0.355109d-01
      cp( 3) = 0.142825d+00
      cp( 4) = 0.338620d+00
      cp( 5) = 0.451579d+00
      cp( 6) = 0.273271d+00
      cp( 7) = 0.906649d-02
      cp( 8) = 0.997202d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          magnesium
c***********************************************************************
c
 1400 go to (1420,1440,1500,1460,1480,1520),ngauss
 1420 continue
 1440 continue
 1460 continue
 1480 continue
      call caserr2('requested basis set not available for magnesium')
 1500 continue
      e( 1) = 0.652841d+03
      e( 2) = 0.983805d+02
      e( 3) = 0.212996d+02
      e( 4) = 0.233727d+02
      e( 5) = 0.519953d+01
      e( 6) = 0.131508d+01
      e( 7) = 0.611349d+00
      e( 8) = 0.141841d+00
      e( 9) = 0.464011d-01
      e(10) = 0.0146d0
      e(11) = 0.175d0
      cs( 1) = 0.675982d-01
      cs( 2) = 0.391778d+00
      cs( 3) = 0.666661d+00
      cs( 4) = -0.110246d+00
      cs( 5) = 0.184119d+00
      cs( 6) = 0.896399d+00
      cp( 1) = 0.121014d+00
      cp( 2) = 0.462810d+00
      cp( 3) = 0.606907d+00
      cs( 7) = -0.361101d+00
      cs( 8) = 0.121505d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 4) = 0.242633d-01
      cp( 5) = 0.986673d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 1520 e( 1) = 0.117228d+05
      e( 2) = 0.175993d+04
      e( 3) = 0.400846d+03
      e( 4) = 0.112807d+03
      e( 5) = 0.359997d+02
      e( 6) = 0.121828d+02
      e( 7) = 0.189180d+03
      e( 8) = 0.452119d+02
      e( 9) = 0.143563d+02
      e(10) = 0.513886d+01
      e(11) = 0.190652d+01
      e(12) = 0.705887d+00
      e(13) = 0.611349d+00
      e(14) = 0.141841d+00
      e(15) = 0.464011d-01
      e(16) = 0.0146d0
      e(17) = 0.175d0
      cs( 1) = 0.197783d-02
      cs( 2) = 0.151140d-01
      cs( 3) = 0.739108d-01
      cs( 4) = 0.249191d+00
      cs( 5) = 0.487928d+00
      cs( 6) = 0.319662d+00
      cs( 7) = -0.323717d-02
      cs( 8) = -0.410079d-01
      cs( 9) = -0.112600d+00
      cs(10) = 0.148633d+00
      cs(11) = 0.616497d+00
      cs(12) = 0.364829d+00
      cs(13) = -0.361101d+00
      cs(14) = 0.121505d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.492813d-02
      cp( 2) = 0.349888d-01
      cp( 3) = 0.140725d+00
      cp( 4) = 0.333642d+00
      cp( 5) = 0.444940d+00
      cp( 6) = 0.269254d+00
      cp( 7) = 0.242633d-01
      cp( 8) = 0.986673d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c***********************************************************************
c          aluminium
c***********************************************************************
c
 1540 go to (1560,1580,1640,1600,1620,1660),ngauss
 1560 continue
 1580 continue
 1600 continue
 1620 continue
      call caserr2('requested basis set not available for aluminium')
 1640 e( 1) = 0.775737d+03
      e( 2) = 0.116952d+03
      e( 3) = 0.253326d+02
      e( 4) = 0.294796d+02
      e( 5) = 0.663314d+01
      e( 6) = 0.172675d+01
      e( 7) = 0.946160d+00
      e( 8) = 0.202506d+00
      e( 9) = 0.639088d-01
      e(10) = 0.0318d0
      e(11) = 0.325d0
      cs( 1) = 0.668347d-01
      cs( 2) = 0.389061d+00
      cs( 3) = 0.669468d+00
      cs( 4) = -0.107902d+00
      cs( 5) = 0.146245d+00
      cs( 6) = 0.923730d+00
      cp( 1) = 0.117574d+00
      cp( 2) = 0.461174d+00
      cp( 3) = 0.605535d+00
      cs( 7) = -0.320327d+00
      cs( 8) = 0.118412d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 4) = 0.519383d-01
      cp( 5) = 0.972660d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 1660 continue
      e( 1) = 0.139831d+05
      e( 2) = 0.209875d+04
      e( 3) = 0.477705d+03
      e( 4) = 0.134360d+03
      e( 5) = 0.428709d+02
      e( 6) = 0.145189d+02
      e( 7) = 0.239668d+03
      e( 8) = 0.574419d+02
      e( 9) = 0.182859d+02
      e(10) = 0.659914d+01
      e(11) = 0.249049d+01
      e(12) = 0.944545d+00
      e(13) = 0.946160d+00
      e(14) = 0.202506d+00
      e(15) = 0.639088d-01
      e(16) = 0.0318d0
      e(17) = 0.325d0
      cs( 1) = 0.194267d-02
      cs( 2) = 0.148599d-01
      cs( 3) = 0.728494d-01
      cs( 4) = 0.246830d+00
      cs( 5) = 0.487258d+00
      cs( 6) = 0.323496d+00
      cs( 7) = -0.292619d-02
      cs( 8) = -0.374083d-01
      cs( 9) = -0.114487d+00
      cs(10) = 0.115635d+00
      cs(11) = 0.612595d+00
      cs(12) = 0.393799d+00
      cs(13) = -0.320327d+00
      cs(14) = 0.118412d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.460285d-02
      cp( 2) = 0.331990d-01
      cp( 3) = 0.136282d+00
      cp( 4) = 0.330476d+00
      cp( 5) = 0.449146d+00
      cp( 6) = 0.265704d+00
      cp( 7) = 0.519383d-01
      cp( 8) = 0.972660d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c********************************************************************
c          silicon
c********************************************************************
c
 1680 go to (1700,1720,1780,1740,1760,1800),ngauss
 1700 continue
 1720 continue
 1740 continue
 1760 continue
      call caserr2('requested basis set not available for silicon')
 1780 e( 1) = 0.910655d+03
      e( 2) = 0.137336d+03
      e( 3) = 0.297601d+02
      e( 4) = 0.366716d+02
      e( 5) = 0.831729d+01
      e( 6) = 0.221645d+01
      e( 7) = 0.107913d+01
      e( 8) = 0.302422d+00
      e( 9) = 0.933392d-01
      e(10) = 0.0331d0
      e(11) = 0.395d0
      cs( 1) = 0.660823d-01
      cs( 2) = 0.386229d+00
      cs( 3) = 0.672380d+00
      cs( 4) = -0.104511d+00
      cs( 5) = 0.107410d+00
      cs( 6) = 0.951446d+00
      cp( 1) = 0.113355d+00
      cp( 2) = 0.457578d+00
      cp( 3) = 0.607427d+00
      cs( 7) = -0.376108d+00
      cs( 8) = 0.125165d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 4) = 0.671030d-01
      cp( 5) = 0.956883d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 1800 continue
      e( 1) = 0.161159d+05
      e( 2) = 0.242558d+04
      e( 3) = 0.553867d+03
      e( 4) = 0.156340d+03
      e( 5) = 0.500683d+02
      e( 6) = 0.170178d+02
      e( 7) = 0.292718d+03
      e( 8) = 0.698731d+02
      e( 9) = 0.223363d+02
      e(10) = 0.815039d+01
      e(11) = 0.313458d+01
      e(12) = 0.122543d+01
      e(13) = 0.107913d+01
      e(14) = 0.302422d+00
      e(15) = 0.933392d-01
      e(16) = 0.0331d0
      e(17) = 0.395d0
      cs( 1) = 0.195948d-02
      cs( 2) = 0.149288d-01
      cs( 3) = 0.728478d-01
      cs( 4) = 0.246130d+00
      cs( 5) = 0.485914d+00
      cs( 6) = 0.325002d+00
      cs( 7) = -0.278094d-02
      cs( 8) = -0.357146d-01
      cs( 9) = -0.114985d+00
      cs(10) = 0.935634d-01
      cs(11) = 0.603017d+00
      cs(12) = 0.418959d+00
      cs(13) = -0.376108d+00
      cs(14) = 0.125165d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.443826d-02
      cp( 2) = 0.326679d-01
      cp( 3) = 0.134721d+00
      cp( 4) = 0.328678d+00
      cp( 5) = 0.449640d+00
      cp( 6) = 0.261372d+00
      cp( 7) = 0.671030d-01
      cp( 8) = 0.956883d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c********************************************************************
c          phosphorus
c***********************************************************************
c
 1820 go to (1840,1860,1920,1880,1900,1940),ngauss
 1840 continue
 1860 continue
 1880 continue
 1900 continue
      call caserr2('requested basis set not available for phosphorus')
 1920 continue
      e( 1) = 0.105490d+04
      e( 2) = 0.159195d+03
      e( 3) = 0.345304d+02
      e( 4) = 0.442866d+02
      e( 5) = 0.101019d+02
      e( 6) = 0.273997d+01
      e( 7) = 0.121865d+01
      e( 8) = 0.395546d+00
      e( 9) = 0.122811d+00
      e(10) = 0.0348d0
      e(11) = 0.55d0
      cs( 1) = 0.655407d-01
      cs( 2) = 0.384036d+00
      cs( 3) = 0.674541d+00
      cs( 4) = -0.102130d+00
      cs( 5) = 0.815922d-01
      cs( 6) = 0.969788d+00
      cs( 7) = -0.371495d+00
      cs( 8) = 0.127099d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 1) = 0.110851d+00
      cp( 2) = 0.456495d+00
      cp( 3) = 0.606936d+00
      cp( 4) = 0.915823d-01
      cp( 5) = 0.934924d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 1940 continue
      e( 1) = 0.194133d+05
      e( 2) = 0.290942d+04
      e( 3) = 0.661364d+03
      e( 4) = 0.185759d+03
      e( 5) = 0.591943d+02
      e( 6) = 0.200310d+02
      e( 7) = 0.339478d+03
      e( 8) = 0.810101d+02
      e( 9) = 0.258780d+02
      e(10) = 0.945221d+01
      e(11) = 0.366566d+01
      e(12) = 0.146746d+01
      e(13) = 0.121865d+01
      e(14) = 0.395546d+00
      e(15) = 0.122811d+00
      e(16) = 0.0348d0
      e(17) = 0.55d0
      cs( 1) = 0.185160d-02
      cs( 2) = 0.142062d-01
      cs( 3) = 0.699995d-01
      cs( 4) = 0.240079d+00
      cs( 5) = 0.484762d+00
      cs( 6) = 0.335200d+00
      cs( 7) = -0.278217d-02
      cs( 8) = -0.360499d-01
      cs( 9) = -0.116631d+00
      cs(10) = 0.968328d-01
      cs(11) = 0.614418d+00
      cs(12) = 0.403798d+00
      cs(13) = -0.371495d+00
      cs(14) = 0.127099d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.456462d-02
      cp( 2) = 0.336936d-01
      cp( 3) = 0.139755d+00
      cp( 4) = 0.339362d+00
      cp( 5) = 0.450921d+00
      cp( 6) = 0.238586d+00
      cp( 7) = 0.915823d-01
      cp( 8) = 0.934924d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c********************************************************************
c          sulphur
c*******************************************************************
c
 1960 go to (1980,2000,2060,2020,2040,2080),ngauss
 1980 continue
 2000 continue
 2020 continue
 2040 continue
      call caserr2('requested basis set not available for sulphur')
 2060 continue
      e( 1) = 0.121062d+04
      e( 2) = 0.182747d+03
      e( 3) = 0.396673d+02
      e( 4) = 0.522236d+02
      e( 5) = 0.119629d+02
      e( 6) = 0.328911d+01
      e( 7) = 0.122384d+01
      e( 8) = 0.457303d+00
      e( 9) = 0.142269d+00
      e(10) = 0.0405d0
      e(11) = 0.65d0
      cs( 1) = 0.650071d-01
      cs( 2) = 0.382040d+00
      cs( 3) = 0.676545d+00
      cs( 4) = -0.100310d+00
      cs( 5) = 0.650877d-01
      cs( 6) = 0.981455d+00
      cs( 7) = -0.286089d+00
      cs( 8) = 0.122806d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 1) = 0.109646d+00
      cp( 2) = 0.457649d+00
      cp( 3) = 0.604261d+00
      cp( 4) = 0.164777d+00
      cp( 5) = 0.870855d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 2080 continue
      e( 1) = 0.219171d+05
      e( 2) = 0.330149d+04
      e( 3) = 0.754146d+03
      e( 4) = 0.212711d+03
      e( 5) = 0.679896d+02
      e( 6) = 0.230515d+02
      e( 7) = 0.423735d+03
      e( 8) = 0.100710d+03
      e( 9) = 0.321599d+02
      e(10) = 0.118079d+02
      e(11) = 0.463110d+01
      e(12) = 0.187025d+01
      e(13) = 0.122384d+01
      e(14) = 0.457303d+00
      e(15) = 0.142269d+00
      e(16) = 0.0405d0
      e(17) = 0.65d0
      cs( 1) = 0.186924d-02
      cs( 2) = 0.142303d-01
      cs( 3) = 0.696962d-01
      cs( 4) = 0.238487d+00
      cs( 5) = 0.483307d+00
      cs( 6) = 0.338074d+00
      cs( 7) = -0.237677d-02
      cs( 8) = -0.316930d-01
      cs( 9) = -0.113317d+00
      cs(10) = 0.560900d-01
      cs(11) = 0.592255d+00
      cs(12) = 0.455006d+00
      cs(13) = -0.286089d+00
      cs(14) = 0.122806d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.406101d-02
      cp( 2) = 0.306813d-01
      cp( 3) = 0.130452d+00
      cp( 4) = 0.327205d+00
      cp( 5) = 0.452851d+00
      cp( 6) = 0.256042d+00
      cp( 7) = 0.164777d+00
      cp( 8) = 0.870855d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c*****************************************************************
c          chlorine
c*****************************************************************
c
 2100 go to (2120,2140,2200,2160,2180,2220),ngauss
 2120 continue
 2140 continue
 2160 continue
 2180 continue
      call caserr2('requested basis set not available for chlorine')
 2200 continue
      e( 1) = 0.137640d+04
      e( 2) = 0.207857d+03
      e( 3) = 0.451554d+02
      e( 4) = 0.608014d+02
      e( 5) = 0.139765d+02
      e( 6) = 0.388710d+01
      e( 7) = 0.135299d+01
      e( 8) = 0.526955d+00
      e( 9) = 0.166714d+00
      e(10) = 0.0483d0
      e(11) = 0.75d0
      cs( 1) = 0.645827d-01
      cs( 2) = 0.380363d+00
      cs( 3) = 0.678190d+00
      cs( 4) = -0.987639d-01
      cs( 5) = 0.511338d-01
      cs( 6) = 0.991337d+00
      cs( 7) = -0.222401d+00
      cs( 8) = 0.118252d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 1) = 0.108598d+00
      cp( 2) = 0.458682d+00
      cp( 3) = 0.601962d+00
      cp( 4) = 0.219216d+00
      cp( 5) = 0.822321d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 2220 continue
      e( 1) = 0.251801d+05
      e( 2) = 0.378035d+04
      e( 3) = 0.860474d+03
      e( 4) = 0.242145d+03
      e( 5) = 0.773349d+02
      e( 6) = 0.262470d+02
      e( 7) = 0.491765d+03
      e( 8) = 0.116984d+03
      e( 9) = 0.374153d+02
      e(10) = 0.137834d+02
      e(11) = 0.545215d+01
      e(12) = 0.222588d+01
      e(13) = 0.135299d+01
      e(14) = 0.526955d+00
      e(15) = 0.166714d+00
      e(16) = 0.0483d0
      e(17) = 0.75d0
      cs( 1) = 0.183296d-02
      cs( 2) = 0.140342d-01
      cs( 3) = 0.690974d-01
      cs( 4) = 0.237452d+00
      cs( 5) = 0.483034d+00
      cs( 6) = 0.339856d+00
      cs( 7) = -0.229739d-02
      cs( 8) = -0.307137d-01
      cs( 9) = -0.112528d+00
      cs(10) = 0.450163d-01
      cs(11) = 0.589353d+00
      cs(12) = 0.465206d+00
      cs(13) = -0.222401d+00
      cs(14) = 0.118252d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.398940d-02
      cp( 2) = 0.303177d-01
      cp( 3) = 0.129880d+00
      cp( 4) = 0.327951d+00
      cp( 5) = 0.453527d+00
      cp( 6) = 0.252154d+00
      cp( 7) = 0.219216d+00
      cp( 8) = 0.822321d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
c
c****************************************************************
c          argon
c******************************************************************
c
 2240 go to (2260,2280,2340,2300,2320,2360),ngauss
 2260 continue
 2280 continue
 2300 continue
 2320 continue
      call caserr2('requested basis set not available for argon')
 2340 continue
      e( 1) = 0.155371d+04
      e( 2) = 0.234678d+03
      e( 3) = 0.510121d+02
      e( 4) = 0.700453d+02
      e( 5) = 0.161473d+02
      e( 6) = 0.453492d+01
      e( 7) = 0.154209d+01
      e( 8) = 0.607267d+00
      e( 9) = 0.195373d+00
      e(10) = 0.05d0
      e(11) = 0.85d0
      cs( 1) = 0.641707d-01
      cs( 2) = 0.378797d+00
      cs( 3) = 0.679752d+00
      cs( 4) = -0.974661d-01
      cs( 5) = 0.390569d-01
      cs( 6) = 0.999916d+00
      cs( 7) = -0.176866d+00
      cs( 8) = 0.114690d+01
      cs( 9) = 0.100000d+01
      cs(10) = 0.100000d+01
      cp( 1) = 0.107619d+00
      cp( 2) = 0.459576d+00
      cp( 3) = 0.600041d+00
      cp( 4) = 0.255687d+00
      cp( 5) = 0.789842d+00
      cp( 6) = 0.100000d+01
      cp( 7) = 0.100000d+01
      cd( 1) = 0.100000d+01
      go to 2380
 2360 continue
      e( 1) = 0.283483d+05
      e( 2) = 0.425762d+04
      e( 3) = 0.969857d+03
      e( 4) = 0.273263d+03
      e( 5) = 0.873695d+02
      e( 6) = 0.296867d+02
      e( 7) = 0.575891d+03
      e( 8) = 0.136816d+03
      e( 9) = 0.438098d+02
      e(10) = 0.162094d+02
      e(11) = 0.646084d+01
      e(12) = 0.265114d+01
      e(13) = 0.154209d+01
      e(14) = 0.607267d+00
      e(15) = 0.195373d+00
      e(16) = 0.05d0
      e(17) = 0.85d0
      cs( 1) = 0.182526d-02
      cs( 2) = 0.139686d-01
      cs( 3) = 0.687073d-01
      cs( 4) = 0.236204d+00
      cs( 5) = 0.482214d+00
      cs( 6) = 0.342043d+00
      cs( 7) = -0.215972d-02
      cs( 8) = -0.290775d-01
      cs( 9) = -0.110827d+00
      cs(10) = 0.276999d-01
      cs(11) = 0.577613d+00
      cs(12) = 0.488688d+00
      cs(13) = -0.176866d+00
      cs(14) = 0.114690d+01
      cs(15) = 0.100000d+01
      cs(16) = 0.100000d+01
      cp( 1) = 0.380665d-02
      cp( 2) = 0.292305d-01
      cp( 3) = 0.126467d+00
      cp( 4) = 0.323510d+00
      cp( 5) = 0.454896d+00
      cp( 6) = 0.256630d+00
      cp( 7) = 0.255687d+00
      cp( 8) = 0.789842d+00
      cp( 9) = 0.100000d+01
      cp(10) = 0.100000d+01
      cd( 1) = 0.100000d+01
 2380 continue
      return
c
      end
      subroutine dunbas(ztype,csinp,cpinp,cdinp,opol,oryd,
     +      sc,scc,nucz,intyp,nangm,
     +      nbfs,minf,maxf,loc,ngauss,ns,ierr1,ierr2,nat)
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
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
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      dimension csinp(*),cpinp(*),cdinp(*)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      dimension scc(*)
      data pt75 /0.75d+00/
      data done/1.0d+00/
      data pt5,pi32,tol /0.5d+00,
     + 5.56832799683170d+00,1.0d-10/
c
c     ----- get exponents, contraction coefficients, and
c          scale factors. -----
c
      ng = -2**20
      igauss = ng
      ityp = ng
      ierr1=0
      ierr2=0
      odone = .false.
      call vclr(eex,1,200)
c
c     decide whether sv (i.e.dzv), dz or tzv
c
      osv = .false.
      odz = .false.
      otz = .false.
      if(ztype.eq.'sv') then
        osv = .true.
      else if(ztype.eq.'dz') then
        odz = .true.
      else if(ztype.eq.'tz') then
        otz = .true.
      else
        call caserr2('unrecognised dunning basis set')
      endif
c
      if (nucz .le. 2) then
c
c     ----- hydrogen -----
c
       if(osv.or.odz) then
        scalh1 = 1.20d0
        scalh2 = 1.15d0
        scalh3 = 1.00d0
        call ddzero(eex,ccs,ccp,nucz)
       else if (otz) then
        scalh1 = 1.00d0
        scalh2 = 1.00d0
        scalh3 = 1.00d0
        call fzero(eex,ccs,ccp,nucz)
       endif
       if (sc .le. 0.0d0) sc = scalh1
       if (scc(1) .le. 0.0d0) scc(1) = scalh2
       if (scc(2) .le. 0.0d0) scc(2) = scalh3
       eex(1) = eex(1)*sc**2
       eex(2) = eex(2)*sc**2
       eex(3) = eex(3)*sc**2
       eex(4) = eex(4)*scc(1)**2
       eex(5) = eex(5)*scc(2)**2
c
      else if (nucz .le. 10) then
c
c     ----- lithium to neon -----
c
       if(osv) then
         call ddone(eex,ccs,ccp,ccd,nucz)
       else if (odz) then
         call dzone(eex,ccs,ccp,ccd,nucz)
       else if (otz) then
         call fone(eex,ccs,ccp,ccd,nucz)
       endif
c
      else if (nucz .le. 18) then
c
c     ----- phosphorus to argon -----
c
       if(osv.or.odz) then
         call ddtwo(eex,ccs,ccp,ccd,nucz)
       else if (otz) then
         call ftwo(eex,ccs,ccp,ccd,nucz)
       endif
c
      else if (nucz .le. 30) then
c
c     ----- potassium to zinc -----
c
       if(oryd)call caserr2(
     +   'requested rydberg basis set not available')
       if(osv.or.odz) then
         call ddthree(eex,ccs,ccp,ccd,nucz)
       else if (otz) then
         call fthree(eex,ccs,ccp,ccd,nucz)
       endif
c
      else if (nucz .le. 36) then
c
c     ----- gallium to krypton -----
c
       if(oryd)call caserr2(
     +   'requested rydberg basis set not available')
       if(osv.or.odz) then
         call ddfour(eex,ccs,ccp,ccd,nucz)
       else if (otz) then
         call ffour(eex,ccs,ccp,ccd,nucz)
       endif
c
      else
       call caserr2('requested dunning  basis set not available')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  180 ipass = ipass+1
      if (osv) then
        call svsh(nucz,ipass,opol,oryd,ityp,igauss,ng,odone)
      else if(odz) then
        call dzsh(nucz,ipass,opol,oryd,ityp,igauss,ng,odone)
      else if(otz) then
        call tzvsh(nucz,ipass,opol,oryd,ityp,igauss,ng,odone)
      else
        call caserr2('unrecognised dunning basis set')
      endif
c
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      if (nshell.gt.mxshel) ierr1 = 1
      if (ierr1.ne.0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss.gt.mxprim) ierr2 = 1
      if (ierr2.ne.0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do 440 i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i) 
         if(ityp.eq.16) csinp(k) = ccs(ng+i)
         if(ityp.eq.17) cpinp(k) = ccp(ng+i)
         if(ityp.eq.18) cdinp(k) = ccd(ng+i)
         cs(k) = 0.0d0
         cp(k) = 0.0d0
         cd(k) = 0.0d0
         cf(k) = 0.0d0
         if(ityp.eq.16) cs(k) = csinp(k)
         if(ityp.eq.17) cp(k) = cpinp(k)
         if(ityp.eq.18) cd(k) = cdinp(k)
  440 continue
c
c     ----- always unnormalize primitives -----
c
      do 460 k = k1,k2
         ee = ex(k)+ex(k)
         facs = pi32/(ee*sqrt(ee))
         facp = pt5*facs/ee
         facd = pt75*facs/(ee*ee)
         if(ityp.eq.16) cs(k) = cs(k)/sqrt(facs)
         if(ityp.eq.17) cp(k) = cp(k)/sqrt(facp)
         if(ityp.eq.18) cd(k) = cd(k)/sqrt(facd)
  460 continue
c
c     ----- if(normf.eq.0) normalize basis functions -----
c
      if (normf .eq. 1) go to 180
      facs = 0.0d0
      facp = 0.0d0
      facd = 0.0d0
      do 510 ig = k1,k2
         do 500 jg = k1,ig
            ee = ex(ig)+ex(jg)
            fac = ee*sqrt(ee)
            dums = cs(ig)*cs(jg)/fac
            dump = pt5*cp(ig)*cp(jg)/(ee*fac)
            dumd = pt75*cd(ig)*cd(jg)/(ee*ee*fac)
            if (ig .eq. jg) go to 480
               dums = dums+dums
               dump = dump+dump
               dumd = dumd+dumd
  480       continue
            facs = facs+dums
            facp = facp+dump
            facd = facd+dumd
  500    continue
  510 continue
c
      fac=0.0d0
      if(ityp.eq.16.and. facs.gt.tol) fac=done/sqrt(facs*pi32)
      if(ityp.eq.17.and. facp.gt.tol) fac=done/sqrt(facp*pi32)
      if(ityp.eq.18.and. facd.gt.tol) fac=done/sqrt(facd*pi32)
c
      do 550 ig = k1,k2
         if(ityp.eq.16) cs(ig) = fac*cs(ig)
         if(ityp.eq.17) cp(ig) = fac*cp(ig)
         if(ityp.eq.18) cd(ig) = fac*cd(ig)
  550 continue
      go to 180
c
  220 continue
      return
      end
c
c
      subroutine ddone(e,s,p,d,n)
c
c     ----- dunning's contraction of huzinaga's (9s,5p) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li -----
c
  100 continue
      e(1) = 9.213d+02
      s(1) = 0.001367d+00
      e(2) = 1.387d+02
      s(2) = 0.010425d+00
      e(3) = 3.194d+01
      s(3) = 0.049859d+00
      e(4) = 9.353d+00
      s(4) = 0.160701d+00
      e(5) = 3.158d+00
      s(5) = 0.344604d+00
      e(6) = 1.157d+00
      s(6) = 0.425197d+00
      e(7) = 4.446d-01
      s(7) = 0.169468d+00
      e(8) = 4.446d-01
      s(8) = -0.222311d+00
      e(9) = 7.666d-02
      s(9) = 1.116477d+00
      e(10) = 2.864d-02
      s(10) = 1.000000d+00
      e(11) = 1.488d+00
      p(11) = 0.038770d+00
      e(12) = 2.667d-01
      p(12) = 0.236257d+00
      e(13) = 7.201d-02
      p(13) = 0.830448d+00
      e(14) = 2.370d-02
      p(14) = 1.000000d+00
c
      e(15) = 2.000d-01
      d(15) = 1.000000d+00
c
      e(16) = 8.0d-03
      s(16) = 1.000000d+00
      e(17) = 6.0d-03
      p(17) = 1.000000d+00
      e(18) = 6.0d-03
      d(18) = 1.000000d+00
c
      return
c
c     ----- be -----
c
  120 continue
      e(1) = 1.741d+03
      s(1) = 0.001305d+00
      e(2) = 2.621d+02
      s(2) = 0.009955d+00
      e(3) = 6.033d+01
      s(3) = 0.048031d+00
      e(4) = 1.762d+01
      s(4) = 0.158577d+00
      e(5) = 5.933d+00
      s(5) = 0.351325d+00
      e(6) = 2.185d+00
      s(6) = 0.427006d+00
      e(7) = 8.590d-01
      s(7) = 0.160490d+00
      e(8) = 2.185d+00
      s(8) = -0.185294d+00
      e(9) = 1.806d-01
      s(9) = 1.057014d+00
      e(10) = 5.835d-02
      s(10) = 1.000000d+00
      e(11) = 6.710d+00
      p(11) = 0.016378d+00
      e(12) = 1.442d+00
      p(12) = 0.091553d+00
      e(13) = 4.103d-01
      p(13) = 0.341469d+00
      e(14) = 1.397d-01
      p(14) = 0.685428d+00
      e(15) = 4.922d-02
      p(15) = 1.000000d+00
c
      e(16) = 2.550d-01
      d(16) = 1.000000d+00
c
      e(17) = 1.6d-02
      s(17) = 1.000000d+00
      e(18) = 1.2d-02
      p(18) = 1.000000d+00
      e(19) = 1.2d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- b  -----
c
  140 continue
      e(1) = 2.788d+03
      s(1) = 0.001288d+00
      e(2) = 4.190d+02
      s(2) = 0.009835d+00
      e(3) = 9.647d+01
      s(3) = 0.047648d+00
      e(4) = 2.807d+01
      s(4) = 0.160069d+00
      e(5) = 9.376d+00
      s(5) = 0.362894d+00
      e(6) = 3.406d+00
      s(6) = 0.433582d+00
      e(7) = 1.306d+00
      s(7) = 0.140082d+00
      e(8) = 3.406d+00
      s(8) = -0.179330d+00
      e(9) = 3.245d-01
      s(9) = 1.062594d+00
      e(10) = 1.022d-01
      s(10) = 1.000000d+00
      e(11) = 1.134d+01
      p(11) = 0.017988d+00
      e(12) = 2.436d+00
      p(12) = 0.110343d+00
      e(13) = 6.836d-01
      p(13) = 0.383072d+00
      e(14) = 2.134d-01
      p(14) = 0.647895d+00
      e(15) = 7.011d-02
      p(15) = 1.000000d+00
c
      e(16) = 7.000d-01
      d(16) = 1.000000d+00
c
      e(17) = 1.9d-02
      s(17) = 1.000000d+00
      e(18) = 1.5d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- c  -----
c
  160 continue
      e(1) = 4.233d+03
      s(1) = 0.001220d+00
      e(2) = 6.349d+02
      s(2) = 0.009342d+00
      e(3) = 1.461d+02
      s(3) = 0.045452d+00
      e(4) = 4.250d+01
      s(4) = 0.154657d+00
      e(5) = 1.419d+01
      s(5) = 0.358866d+00
      e(6) = 5.148d+00
      s(6) = 0.438632d+00
      e(7) = 1.967d+00
      s(7) = 0.145918d+00
      e(8) = 5.148d+00
      s(8) = -0.168367d+00
      e(9) = 4.962d-01
      s(9) = 1.060091d+00
      e(10) = 1.533d-01
      s(10) = 1.000000d+00
      e(11) = 1.816d+01
      p(11) = 0.018539d+00
      e(12) = 3.986d+00
      p(12) = 0.115436d+00
      e(13) = 1.143d+00
      p(13) = 0.386188d+00
      e(14) = 3.594d-01
      p(14) = 0.640114d+00
      e(15) = 1.146d-01
      p(15) = 1.000000d+00
c
      e(16) = 7.500d-01
      d(16) = 1.000000d+00
c
      e(17) = 2.3d-02
      s(17) = 1.000000d+00
      e(18) = 2.1d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- n -----
c
  180 continue
      e(1) = 5.909d+03
      s(1) = 0.001190d+00
      e(2) = 8.875d+02
      s(2) = 0.009099d+00
      e(3) = 2.047d+02
      s(3) = 0.044145d+00
      e(4) = 5.984d+01
      s(4) = 0.150464d+00
      e(5) = 2.000d+01
      s(5) = 0.356741d+00
      e(6) = 7.193d+00
      s(6) = 0.446533d+00
      e(7) = 2.686d+00
      s(7) = 0.145603d+00
      e(8) = 7.193d+00
      s(8) = -0.160405d+00
      e(9) = 7.000d-01
      s(9) = 1.058215d+00
      e(10) = 2.133d-01
      s(10) = 1.000000d+00
      e(11) = 2.679d+01
      p(11) = 0.018254d+00
      e(12) = 5.956d+00
      p(12) = 0.116461d+00
      e(13) = 1.707d+00
      p(13) = 0.390178d+00
      e(14) = 5.314d-01
      p(14) = 0.637102d+00
      e(15) = 1.654d-01
      p(15) = 1.000000d+00
c
      e(16) = 8.000d-01
      d(16) = 1.000000d+00
c
      e(17) = 2.8d-02
      s(17) = 1.000000d+00
      e(18) = 2.5d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- o  ------
c
  200 continue
      e(1) = 7.817d+03
      s(1) = 0.001176d+00
      e(2) = 1.176d+03
      s(2) = 0.008968d+00
      e(3) = 2.732d+02
      s(3) = 0.042868d+00
      e(4) = 8.117d+01
      s(4) = 0.143930d+00
      e(5) = 2.718d+01
      s(5) = 0.355630d+00
      e(6) = 9.532d+00
      s(6) = 0.461248d+00
      e(7) = 3.414d+00
      s(7) = 0.140206d+00
      e(8) = 9.532d+00
      s(8) = -0.154153d+00
      e(9) = 9.398d-01
      s(9) = 1.056914d+00
      e(10) = 2.846d-01
      s(10) = 1.000000d+00
      e(11) = 3.518d+01
      p(11) = 0.019580d+00
      e(12) = 7.904d+00
      p(12) = 0.124200d+00
      e(13) = 2.305d+00
      p(13) = 0.394714d+00
      e(14) = 7.171d-01
      p(14) = 0.627376d+00
      e(15) = 2.137d-01
      p(15) = 1.000000d+00
c
      e(16) = 8.500d-01
      d(16) = 1.000000d+00
c
      e(17) = 3.2d-02
      s(17) = 1.000000d+00
      e(18) = 2.8d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- f  -----
c
  220 continue
      e(1) = 9.995d+03
      s(1) = 0.001166d+00
      e(2) = 1.506d+03
      s(2) = 0.008876d+00
      e(3) = 3.503d+02
      s(3) = 0.042380d+00
      e(4) = 1.041d+02
      s(4) = 0.142929d+00
      e(5) = 3.484d+01
      s(5) = 0.355372d+00
      e(6) = 1.222d+01
      s(6) = 0.462085d+00
      e(7) = 4.369d+00
      s(7) = 0.140848d+00
      e(8) = 1.222d+01
      s(8) = -0.148452d+00
      e(9) = 1.208d+00
      s(9) = 1.055270d+00
      e(10) = 3.634d-01
      s(10) = 1.000000d+00
      e(11) = 4.436d+01
      p(11) = 0.020876d+00
      e(12) = 1.008d+01
      p(12) = 0.130107d+00
      e(13) = 2.996d+00
      p(13) = 0.396166d+00
      e(14) = 9.383d-01
      p(14) = 0.620404d+00
      e(15) = 2.733d-01
      p(15) = 1.000000d+00
c
      e(16) = 9.000d-01
      d(16) = 1.000000d+00
c
      e(17) = 3.6d-02
      s(17) = 1.000000d+00
      e(18) = 2.9d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- ne -----
c
  240 continue
      e(1) = 1.210d+04
      s(1) = 0.001200d+00
      e(2) = 1.821d+03
      s(2) = 0.009092d+00
      e(3) = 4.328d+02
      s(3) = 0.041305d+00
      e(4) = 1.325d+02
      s(4) = 0.137867d+00
      e(5) = 4.377d+01
      s(5) = 0.362433d+00
      e(6) = 1.491d+01
      s(6) = 0.472247d+00
      e(7) = 5.127d+00
      s(7) = 0.130035d+00
      e(8) = 1.491d+01
      s(8) = -0.140810d+00
      e(9) = 1.491d+00
      s(9) = 1.053327d+00
      e(10) = 4.468d-01
      s(10) = 1.000000d+00
      e(11) = 5.645d+01
      p(11) = 0.020875d+00
      e(12) = 1.292d+01
      p(12) = 0.130032d+00
      e(13) = 3.865d+00
      p(13) = 0.395679d+00
      e(14) = 1.203d+00
      p(14) = 0.621450d+00
      e(15) = 3.444d-01
      p(15) = 1.000000d+00
c
      e(16) = 1.000d-00
      d(16) = 1.000000d+00
c
      e(17) = 4.0d-02
      s(17) = 1.000000d+00
      e(18) = 3.0d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
      end
c
      subroutine ddtwo(e,s,p,d,n)
c
c     ----- dunning's contraction of huzinaga's (11s,7p) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
  100 call caserr2('requested basis set not available for sodium')
  120 call caserr2('requested basis set not available for magnesium')
c
c     ----- al -----
c
  140 continue
      e(1) = 2.349d+04
      s(1) = 0.002509d+00
      e(2) = 3.548d+03
      s(2) = 0.018986d+00
      e(3) = 8.235d+02
      s(3) = 0.092914d+00
      e(4) = 2.377d+02
      s(4) = 0.335935d+00
      e(5) = 7.860d+01
      s(5) = 0.647391d+00
      e(6) = 7.860d+01
      s(6) = 0.111937d+00
      e(7) = 2.905d+01
      s(7) = 0.655976d+00
      e(8) = 1.162d+01
      s(8) = 0.283349d+00
      e(9) = 3.465d+00
      s(9) = 1.000000d+00
      e(10) = 1.233d+00
      s(10) = 1.000000d+00
      e(11) = 2.018d-01
      s(11) = 1.000000d+00
      e(12) = 7.805d-02
      s(12) = 1.000000d+00
      e(13) = 1.415d+02
      p(13) = 0.017882d+00
      e(14) = 3.322d+01
      p(14) = 0.120375d+00
      e(15) = 1.039d+01
      p(15) = 0.411580d+00
      e(16) = 3.593d+00
      p(16) = 0.595353d+00
      e(17) = 3.593d+00
      p(17) = 0.211758d+00
      e(18) = 1.242d+00
      p(18) = 0.837795d+00
      e(19) = 3.040d-01
      p(19) = 1.000000d+00
      e(20) = 7.629d-02
      p(20) = 1.000000d+00
c
      e(21) = 2.5d-01
      d(21) = 1.000000d+00
c
      e(22) = 1.4d-02
      s(22) = 1.000000d+00
      e(23) = 1.1d-02
      p(23) = 1.000000d+00
      e(24) = 1.5d-02
      d(24) = 1.000000d+00
c
      return
c
c     ----- si -----
c
  160 continue
      e(1) = 2.674d+04
      s(1) = 0.002583d+00
      e(2) = 4.076d+03
      s(2) = 0.019237d+00
      e(3) = 9.533d+02
      s(3) = 0.093843d+00
      e(4) = 2.746d+02
      s(4) = 0.341235d+00
      e(5) = 9.068d+01
      s(5) = 0.641675d+00
      e(6) = 9.068d+01
      s(6) = 0.121439d+00
      e(7) = 3.353d+01
      s(7) = 0.653143d+00
      e(8) = 1.346d+01
      s(8) = 0.277624d+00
      e(9) = 4.051d+00
      s(9) = 1.000000d+00
      e(10) = 1.484d+00
      s(10) = 1.000000d+00
      e(11) = 2.704d-01
      s(11) = 1.000000d+00
      e(12) = 9.932d-02
      s(12) = 1.000000d+00
      e(13) = 1.637d+02
      p(13) = 0.011498d+00
      e(14) = 3.835d+01
      p(14) = 0.077726d+00
      e(15) = 1.202d+01
      p(15) = 0.263595d+00
      e(16) = 4.185d+00
      p(16) = 0.758269d+00
      e(17) = 4.185d+00
      p(17) = -1.173045d+00
      e(18) = 1.483d+00
      p(18) = 1.438335d+00
      e(19) = 3.350d-01
      p(19) = 1.000000d+00
      e(20) = 9.699d-02
      p(20) = 1.000000d+00
c
      e(21) = 3.25d-01
      d(21) = 1.000000d+00
c
      e(22) = 1.7d-02
      s(22) = 1.000000d+00
      e(23) = 1.4d-02
      p(23) = 1.000000d+00
      e(24) = 1.5d-02
      d(24) = 1.000000d+00
c
      return
c
c     ----- p  -----
c
  180 continue
      e(1) = 3.063d+04
      s(1) = 0.002619d+00
      e(2) = 4.684d+03
      s(2) = 0.019479d+00
      e(3) = 1.094d+03
      s(3) = 0.095207d+00
      e(4) = 3.153d+02
      s(4) = 0.345742d+00
      e(5) = 1.041d+02
      s(5) = 0.636288d+00
      e(6) = 1.041d+02
      s(6) = 0.130706d+00
      e(7) = 3.842d+01
      s(7) = 0.650274d+00
      e(8) = 1.545d+01
      s(8) = 0.272308d+00
      e(9) = 4.656d+00
      s(9) = 1.000000d+00
      e(10) = 1.759d+00
      s(10) = 1.000000d+00
      e(11) = 3.409d-01
      s(11) = 1.000000d+00
      e(12) = 1.238d-01
      s(12) = 1.000000d+00
      e(13) = 1.877d+02
      p(13) = 0.013158d+00
      e(14) = 4.363d+01
      p(14) = 0.090494d+00
      e(15) = 1.360d+01
      p(15) = 0.305054d+00
      e(16) = 4.766d+00
      p(16) = 0.713579d+00
      e(17) = 4.766d+00
      p(17) = -0.792573d+00
      e(18) = 1.743d+00
      p(18) = 1.429987d+00
      e(19) = 4.192d-01
      p(19) = 1.000000d+00
      e(20) = 1.245d-01
      p(20) = 1.000000d+00
c
      e(21) = 3.70d-01
      d(21) = 1.000000d+00
c
      e(22) = 2.0d-02
      s(22) = 1.000000d+00
      e(23) = 1.7d-02
      p(23) = 1.000000d+00
      e(24) = 1.5d-02
      d(24) = 1.000000d+00
c
      return
c
c     ----- s  -----
c
  200 continue
      e(1) = 3.571d+04
      s(1) = 0.002565d+00
      e(2) = 5.397d+03
      s(2) = 0.019405d+00
      e(3) = 1.250d+03
      s(3) = 0.095595d+00
      e(4) = 3.599d+02
      s(4) = 0.345793d+00
      e(5) = 1.192d+02
      s(5) = 0.635794d+00
      e(6) = 1.192d+02
      s(6) = 0.130096d+00
      e(7) = 4.398d+01
      s(7) = 0.651301d+00
      e(8) = 1.763d+01
      s(8) = 0.271955d+00
      e(9) = 5.420d+00
      s(9) = 1.000000d+00
      e(10) = 2.074d+00
      s(10) = 1.000000d+00
      e(11) = 4.246d-01
      s(11) = 1.000000d+00
      e(12) = 1.519d-01
      s(12) = 1.000000d+00
      e(13) = 2.129d+02
      p(13) = 0.014091d+00
      e(14) = 4.960d+01
      p(14) = 0.096685d+00
      e(15) = 1.552d+01
      p(15) = 0.323874d+00
      e(16) = 5.476d+00
      p(16) = 0.691756d+00
      e(17) = 5.476d+00
      p(17) = -0.626737d+00
      e(18) = 2.044d+00
      p(18) = 1.377051d+00
      e(19) = 5.218d-01
      p(19) = 1.000000d+00
      e(20) = 1.506d-01
      p(20) = 1.000000d+00
c
      e(21) = 5.32d-01
      d(21) = 1.000000d+00
c
      e(22) = 2.3d-02
      s(22) = 1.000000d+00
      e(23) = 2.0d-02
      p(23) = 1.000000d+00
      e(24) = 1.5d-02
      d(24) = 1.000000d+00
c
      return
c
c     ----- cl -----
c
  220 continue
      e(1) = 4.085d+04
      s(1) = 0.002532d+00
      e(2) = 6.179d+03
      s(2) = 0.019207d+00
      e(3) = 1.425d+03
      s(3) = 0.095257d+00
      e(4) = 4.092d+02
      s(4) = 0.345589d+00
      e(5) = 1.355d+02
      s(5) = 0.636401d+00
      e(6) = 1.355d+02
      s(6) = 0.120956d+00
      e(7) = 5.013d+01
      s(7) = 0.648511d+00
      e(8) = 2.021d+01
      s(8) = 0.275487d+00
      e(9) = 6.283d+00
      s(9) = 1.000000d+00
      e(10) = 2.460d+00
      s(10) = 1.000000d+00
      e(11) = 5.271d-01
      s(11) = 1.000000d+00
      e(12) = 1.884d-01
      s(12) = 1.000000d+00
      e(13) = 2.408d+02
      p(13) = 0.014595d+00
      e(14) = 5.656d+01
      p(14) = 0.099047d+00
      e(15) = 1.785d+01
      p(15) = 0.330562d+00
      e(16) = 6.350d+00
      p(16) = 0.682874d+00
      e(17) = 6.350d+00
      p(17) = -0.561785d+00
      e(18) = 2.403d+00
      p(18) = 1.351901d+00
      e(19) = 6.410d-01
      p(19) = 1.000000d+00
      e(20) = 1.838d-01
      p(20) = 1.000000d+00
      e(21) = 6.00d-01
      d(21) = 1.000000d+00
c
      e(22) = 2.5d-02
      s(22) = 1.000000d+00
      e(23) = 2.0d-02
      p(23) = 1.000000d+00
      e(24) = 1.5d-02
      d(24) = 1.000000d+00
c
      return
  240 call caserr2('requested basis set not available for argon')
      return
c
      end
c
      subroutine ddthree(e,s,p,d,n)
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-18
      go to (100,120,140,160,180,200,220,240,260,280,300,320),nn
  100 call caserr2('requested basis set not available for potassium')
  120 call caserr2('requested basis set not available for calcium')
  140 continue
c
c ----- scandium
c
      e(1) = 4.58226d+04
      s(1) = 0.00147d+00
      e(2) = 6.82496d+03
      s(2) = 0.01150d+00
      e(3) = 1.53960d+03
      s(3) = 0.05763d+00
      e(4) = 4.33931d+02
      s(4) = 0.20300d+00
      e(5) = 1.40439d+02
      s(5) = 1.00000d+00
      e(6) = 4.92303d+01
      s(6) = 1.00000d+00
      e(7) = 1.20090d+01
      s(7) = 1.00000d+00
      e(8) = 4.92315d+00
      s(8) = 1.00000d+00
      e(9) = 1.21920d+00
      s(9) = 1.00000d+00
      e(10) = 5.15736d-01
      s(10) = 1.00000d+00
      e(11) = 7.64080d-02
      s(11) = 1.00000d+00
      e(12) = 2.59687d+02
      p(12) = 0.02513d+00
      e(13) = 6.09968d+01
      p(13) = 0.15943d+00
      e(14) = 1.88155d+01
      p(14) = 0.44943d+00
      e(15) = 6.40111d+00
      p(15) = 1.00000d+00
      e(16) = 1.65914d+00
      p(16) = 1.00000d+00
      e(17) = 5.27407d-01
      p(17) = 1.00000d+00
      e(18) = 2.068d-01
      p(18) = 1.00000d+00
      e(19) = 8.11d-02
      p(19) = 1.00000d+00
      e(20) = 1.519d+01
      d(20) = 2.89190492d-02
      e(21) = 3.710d+00
      d(21) = 1.49174065d-01
      e(22) = 1.084d+00
      d(22) = 3.34660015d-01
      e(23) = 2.974d-01
      d(23) = 4.60697426d-01
      e(24) = 6.936d-02
      d(24) = 1.00000d+00
      return
  160 continue
c
c ----- titanium
c
      e(1) = 4.86800d+04
      s(1) = 0.00156d+00
      e(2) = 7.14404d+03
      s(2) = 0.01212d+00
      e(3) = 1.63802d+03
      s(3) = 0.05931d+00
      e(4) = 4.65527d+02
      s(4) = 0.20638d+00
      e(5) = 1.52394d+02
      s(5) = 1.00000d+00
      e(6) = 5.40349d+01
      s(6) = 1.00000d+00
      e(7) = 1.32356d+01
      s(7) = 1.00000d+00
      e(8) = 5.47417d+00
      s(8) = 1.00000d+00
      e(9) = 1.38870d+00
      s(9) = 1.00000d+00
      e(10) = 5.83370d-01
      s(10) = 1.00000d+00
      e(11) = 8.0663d-02
      s(11) = 1.00000d+00
      e(12) = 2.84565d+02
      p(12) = 0.02568d+00
      e(13) = 6.68999d+01
      p(13) = 0.16259d+00
      e(14) = 2.07458d+01
      p(14) = 0.45281d+00
      e(15) = 7.13330d+00
      p(15) = 1.00000d+00
      e(16) = 1.87140d+00
      p(16) = 1.00000d+00
      e(17) = 5.96298d-01
      p(17) = 1.00000d+00
      e(18) = 2.338d-01
      p(18) = 1.00000d+00
      e(19) = 9.17d-02
      p(19) = 1.00000d+00
      e(20) = 2.073d+01
      d(20) = 2.66517493d-02
      e(21) = 5.090d+00
      d(21) = 1.46864218d-01
      e(22) = 1.545d+00
      d(22) = 3.40392761d-01
      e(23) = 4.585d-01
      d(23) = 4.56804981d-01
      e(24) = 1.102d-01
      d(24) = 1.00000d+00
      return
  180 continue
c
c ----- vanadium
c
      e(1) = 5.02813d+04
      s(1) = 0.00167d+00
      e(2) = 7.42464d+03
      s(2) = 0.01284d+00
      e(3) = 1.71655d+03
      s(3) = 0.06201d+00
      e(4) = 4.92456d+02
      s(4) = 0.21113d+00
      e(5) = 1.63130d+02
      s(5) = 1.00000d+00
      e(6) = 5.82952d+01
      s(6) = 1.00000d+00
      e(7) = 1.48239d+01
      s(7) = 1.00000d+00
      e(8) = 6.09769d+00
      s(8) = 1.00000d+00
      e(9) = 1.57364d+00
      s(9) = 1.00000d+00
      e(10) = 6.50938d-01
      s(10) = 1.00000d+00
      e(11) = 8.8861d-02
      s(11) = 1.00000d+00
      e(12) = 3.11142d+02
      p(12) = 0.02609d+00
      e(13) = 7.32585d+01
      p(13) = 0.16498d+00
      e(14) = 2.28156d+01
      p(14) = 0.45557d+00
      e(15) = 7.91458d+00
      p(15) = 1.00000d+00
      e(16) = 2.10690d+00
      p(16) = 1.00000d+00
      e(17) = 6.74850d-01
      p(17) = 1.00000d+00
      e(18) = 2.646d-01
      p(18) = 1.00000d+00
      e(19) = 1.037d-01
      p(19) = 1.00000d+00
      e(20) = 2.118d+01
      d(20) = 3.35595342d-02
      e(21) = 5.566d+00
      d(21) = 1.64567490d-01
      e(22) = 1.753d+00
      d(22) = 3.65791138d-01
      e(23) = 5.256d-01
      d(23) = 4.58772599d-01
      e(24) = 1.336d-01
      d(24) = 1.00000d+00
      return
  200 continue
c
c ----- chromium
c
      e(1) = 5.56434d+04
      s(1) = 0.00162d+00
      e(2) = 8.24033d+03
      s(2) = 0.01268d+00
      e(3) = 1.87603d+03
      s(3) = 0.06209d+00
      e(4) = 5.37177d+02
      s(4) = 0.21149d+00
      e(5) = 1.77239d+02
      s(5) = 1.00000d+00
      e(6) = 6.33965d+01
      s(6) = 1.00000d+00
      e(7) = 1.62072d+01
      s(7) = 1.00000d+00
      e(8) = 6.73669d+00
      s(8) = 1.00000d+00
      e(9) = 1.74645d+00
      s(9) = 1.00000d+00
      e(10) = 7.17674d-01
      s(10) = 1.00000d+00
      e(11) = 9.4109d-02
      s(11) = 1.00000d+00
      e(12) = 3.57615d+02
      p(12) = 0.02440d+00
      e(13) = 8.36031d+01
      p(13) = 0.15873d+00
      e(14) = 2.59210d+01
      p(14) = 0.45090d+00
      e(15) = 8.97595d+00
      p(15) = 1.00000d+00
      e(16) = 2.41609d+00
      p(16) = 1.00000d+00
      e(17) = 7.66685d-01
      p(17) = 1.00000d+00
      e(18) = 3.006d-01
      p(18) = 1.00000d+00
      e(19) = 1.179d-01
      p(19) = 1.00000d+00
      e(20) = 2.800d+01
      d(20) = 2.82278486d-02
      e(21) = 7.213d+00
      d(21) = 1.53903211d-01
      e(22) = 2.241d+00
      d(22) = 3.68451130d-01
      e(23) = 6.612d-01
      d(23) = 4.68309376d-01
      e(24) = 1.620d-01
      d(24) = 1.00000d+00
      return
  220 continue
c
c ----- manganese
c
      e(1) = 6.03705d+04
      s(1) = 0.00162d+00
      e(2) = 8.91051d+03
      s(2) = 0.01284d+00
      e(3) = 2.00893d+03
      s(3) = 0.06309d+00
      e(4) = 5.79251d+02
      s(4) = 0.21138d+00
      e(5) = 1.92434d+02
      s(5) = 1.00000d+00
      e(6) = 6.90690d+01
      s(6) = 1.00000d+00
      e(7) = 1.77053d+01
      s(7) = 1.00000d+00
      e(8) = 7.39915d+00
      s(8) = 1.00000d+00
      e(9) = 1.94442d+00
      s(9) = 1.00000d+00
      e(10) = 7.91057d-01
      s(10) = 1.00000d+00
      e(11) = 1.01306d-01
      s(11) = 1.00000d+00
      e(12) = 3.83732d+02
      p(12) = 0.02506d+00
      e(13) = 9.05004d+01
      p(13) = 0.16096d+00
      e(14) = 2.82205d+01
      p(14) = 0.45387d+00
      e(15) = 9.82707d+00
      p(15) = 1.00000d+00
      e(16) = 2.64911d+00
      p(16) = 1.00000d+00
      e(17) = 8.39238d-01
      p(17) = 1.00000d+00
      e(18) = 3.290d-01
      p(18) = 1.00000d+00
      e(19) = 1.290d-01
      p(19) = 1.00000d+00
      e(20) = 2.815d+01
      d(20) = 3.42869181d-02
      e(21) = 7.564d+00
      d(21) = 1.70996105d-01
      e(22) = 2.410d+00
      d(22) = 3.79828354d-01
      e(23) = 7.282d-01
      d(23) = 4.58334035d-01
      e(24) = 1.799d-01
      d(24) = 1.00000d+00
      return
  240 continue
c
c ----- iron
c
      e(1) = 6.35098d+04
      s(1) = 0.00169d+00
      e(2) = 9.34459d+03
      s(2) = 0.01320d+00
      e(3) = 2.15666d+03
      s(3) = 0.06284d+00
      e(4) = 6.27320d+02
      s(4) = 0.21100d+00
      e(5) = 2.08380d+02
      s(5) = 1.00000d+00
      e(6) = 7.48588d+01
      s(6) = 1.00000d+00
      e(7) = 1.92876d+01
      s(7) = 1.00000d+00
      e(8) = 8.09403d+00
      s(8) = 1.00000d+00
      e(9) = 2.14813d+00
      s(9) = 1.00000d+00
      e(10) = 8.68997d-01
      s(10) = 1.00000d+00
      e(11) = 1.10224d-01
      s(11) = 1.00000d+00
      e(12) = 4.24039d+02
      p(12) = 0.02472d+00
      e(13) = 9.89427d+01
      p(13) = 0.16155d+00
      e(14) = 3.08593d+01
      p(14) = 0.45455d+00
      e(15) = 1.07975d+01
      p(15) = 1.00000d+00
      e(16) = 2.91785d+00
      p(16) = 1.00000d+00
      e(17) = 9.24889d-01
      p(17) = 1.00000d+00
      e(18) = 3.626d-01
      p(18) = 1.00000d+00
      e(19) = 1.422d-01
      p(19) = 1.00000d+00
      e(20) = 3.178d+01
      d(20) = 3.54066628d-02
      e(21) = 8.426d+00
      d(21) = 1.78115212d-01
      e(22) = 2.719d+00
      d(22) = 3.8343646d-01
      e(23) = 8.335d-01
      d(23) = 4.51535378d-01
      e(24) = 2.113d-01
      d(24) = 1.00000d+00
      return
  260 continue
c
c ----- cobalt
c
      e(1) = 7.12839d+04
      s(1) = 0.00163d+00
      e(2) = 1.03287d+04
      s(2) = 0.01286d+00
      e(3) = 2.36330d+03
      s(3) = 0.06242d+00
      e(4) = 6.77873d+02
      s(4) = 0.21218d+00
      e(5) = 2.24950d+02
      s(5) = 1.00000d+00
      e(6) = 8.10468d+01
      s(6) = 1.00000d+00
      e(7) = 2.10939d+01
      s(7) = 1.00000d+00
      e(8) = 8.87011d+00
      s(8) = 1.00000d+00
      e(9) = 2.36802d+00
      s(9) = 1.00000d+00
      e(10) = 9.53468d-01
      s(10) = 1.00000d+00
      e(11) = 1.18259d-01
      s(11) = 1.00000d+00
      e(12) = 4.59697d+02
      p(12) = 0.02447d+00
      e(13) = 1.09051d+02
      p(13) = 0.15854d+00
      e(14) = 3.39790d+01
      p(14) = 0.45458d+00
      e(15) = 1.18744d+01
      p(15) = 1.00000d+00
      e(16) = 3.24082d+00
      p(16) = 1.00000d+00
      e(17) = 1.02446d+00
      p(17) = 1.00000d+00
      e(18) = 4.017d-01
      p(18) = 1.00000d+00
      e(19) = 1.575d-01
      p(19) = 1.00000d+00
      e(20) = 3.497d+01
      d(20) = 3.62293812d-02
      e(21) = 9.463d+00
      d(21) = 1.81096242d-01
      e(22) = 3.051d+00
      d(22) = 3.90134959d-01
      e(23) = 9.365d-01
      d(23) = 4.50939364d-01
      e(24) = 2.350d-01
      d(24) = 1.00000d+00
      return
  280 continue
c
c ----- nickel
c
      e(1) = 7.38504d+04
      s(1) = 0.00168d+00
      e(2) = 1.09392d+04
      s(2) = 0.01307d+00
      e(3) = 2.50461d+03
      s(3) = 0.06348d+00
      e(4) = 7.20706d+02
      s(4) = 0.21402d+00
      e(5) = 2.40114d+02
      s(5) = 1.00000d+00
      e(6) = 8.68503d+01
      s(6) = 1.00000d+00
      e(7) = 2.24686d+01
      s(7) = 1.00000d+00
      e(8) = 9.50508d+00
      s(8) = 1.00000d+00
      e(9) = 2.57910d+00
      s(9) = 1.00000d+00
      e(10) = 1.02958d+00
      s(10) = 1.00000d+00
      e(11) = 1.25776d-01
      s(11) = 1.00000d+00
      e(12) = 4.93637d+02
      p(12) = 0.02496d+00
      e(13) = 1.16498d+02
      p(13) = 0.16147d+00
      e(14) = 3.65669d+01
      p(14) = 0.45551d+00
      e(15) = 1.28762d+01
      p(15) = 1.00000d+00
      e(16) = 3.54585d+00
      p(16) = 1.00000d+00
      e(17) = 1.12160d+00
      p(17) = 1.00000d+00
      e(18) = 4.398d-01
      p(18) = 1.00000d+00
      e(19) = 1.724d-01
      p(19) = 1.00000d+00
      e(20) = 3.949d+01
      d(20) = 3.55460797d-02
      e(21) = 1.075d+01
      d(21) = 1.79996302d-01
      e(22) = 3.475d+00
      d(22) = 3.92560290d-01
      e(23) = 1.065d+00
      d(23) = 4.54418307d-01
      e(24) = 2.641d-01
      d(24) = 1.00000d+00
      return
  300 continue
c
c ----- cupper
c
      e(1) = 7.65691d+04
      s(1) = 0.00170d+00
      e(2) = 1.19447d+04
      s(2) = 0.01270d+00
      e(3) = 2.72644d+03
      s(3) = 0.06257d+00
      e(4) = 7.85697d+02
      s(4) = 0.21049d+00
      e(5) = 2.61625d+02
      s(5) = 1.00000d+00
      e(6) = 9.43546d+01
      s(6) = 1.00000d+00
      e(7) = 2.45815d+01
      s(7) = 1.00000d+00
      e(8) = 1.04216d+01
      s(8) = 1.00000d+00
      e(9) = 2.80706d+00
      s(9) = 1.00000d+00
      e(10) = 1.12010d+00
      s(10) = 1.00000d+00
      e(11) = 1.31029d-01
      s(11) = 1.00000d+00
      e(12) = 5.37269d+02
      p(12) = 0.02448d+00
      e(13) = 1.27621d+02
      p(13) = 0.15879d+00
      e(14) = 3.99997d+01
      p(14) = 0.45462d+00
      e(15) = 1.40748d+01
      p(15) = 1.00000d+00
      e(16) = 3.88005d+00
      p(16) = 1.00000d+00
      e(17) = 1.22287d+00
      p(17) = 1.00000d+00
      e(18) = 4.795d-01
      p(18) = 1.00000d+00
      e(19) = 1.880d-01
      p(19) = 1.00000d+00
      e(20) = 4.366d+01
      d(20) = 3.56934638d-02
      e(21) = 1.197d+01
      d(21) = 1.80263235d-01
      e(22) = 3.916d+00
      d(22) = 3.89610112d-01
      e(23) = 1.222d+00
      d(23) = 4.50798523d-01
      e(24) = 3.066d-01
      d(24) = 1.00000d+00
      return
  320 call caserr2('requested basis set not available for zinc')
      return
      end
      subroutine ahbas(ztype,csinp,cpinp,cdinp,opol,oryd,
     + sc,scc,nucz,intyp,nangm,nbfs,minf,maxf,loc,ngauss,
     + ns,ngsmax,ierr1,ierr2,nat)
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
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
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax),
     +  intyp(*),nangm(*),nbfs(*),minf(*),maxf(*),ns(*),scc(*)
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      data pt5,pt75/
     + 0.5d+00,0.75d+00/
      data pi32,tm10/5.56832799683170d+00,1.0d-10/
c     data tm6/1.0d-06/
      data scal1,scal2,scal3/1.00d0,1.00d0,1.00d0/
c
c     Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c     SV/DZ/TZ basis sets 
c
      osv = .false.
      odz = .false.
      otz = .false.
      if(ztype.eq.'sv') then
        osv = .true.
      else if(ztype.eq.'dz') then
        odz = .true.
      else if(ztype.eq.'tz') then
        otz = .true.
      else
        call caserr2('unrecognised ahlrichs basis set')
      endif
      if(oryd) call caserr2('no diffuse functions available')
      scale = 0.0d+00
      ng = -2**20
      igauss = ng
      ityp = ng
      ierr1=0
      ierr2=0
      odone = .false.
      call vclr(eex,1,200)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
        if (osv.or.odz) then
          call ahsv0(eex,ccs,ccp,nucz)
        else if(otz) then
          call ahtz0(eex,ccs,ccp,nucz)
        endif
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
        if (osv) then
          call ahsv1(eex,ccs,ccp,ccd,nucz)
        else if(odz) then
          call ahdz1(eex,ccs,ccp,ccd,nucz)
        else if(otz) then
          call ahtz1(eex,ccs,ccp,ccd,nucz)
        endif
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
        if (osv) then
          call ahsv2(eex,ccs,ccp,ccd,nucz)
        else if(odz) then
          call ahdz2(eex,ccs,ccp,ccd,nucz)
        else if(otz) then
          call ahtz2(eex,ccs,ccp,ccd,nucz)
        endif
c
c     ----- potassium to zinc -----
c
      else if(nucz.le.30) then
        if (osv) then
          call ahsv3(eex,ccs,ccp,ccd,nucz)
        else if(odz) then
          call ahdz3(eex,ccs,ccp,ccd,nucz)
        else if(otz) then
          call caserr2('ahlrichs TZ function not available')
        endif
c
c     ----- gallium to krypton -----
c
      else if(nucz.le.36) then
        if (osv) then
          call ahsv4(eex,ccs,ccp,ccd,nucz)
        else if(odz) then
          call ahdz4(eex,ccs,ccp,ccd,nucz)
        else if(otz) then
          call caserr2('ahlrichs TZ function not available')
        endif
c
c     ----- past krypton does not exist
c
      else
        call caserr2(
     +   'attempting to site ahlrichs function on invalid centre')
      endif
c
      if (sc .le. 0.0d0) sc = scal1
      if (scc(1) .le. 0.0d0) scc(1) = scal2
      if (scc(2) .le. 0.0d0) scc(2) = scal3
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
      scs = sc
      scp = scc(1)
      scd = scc(2)
      if(osv) then
        call ahsvsh(nucz,ipass,opol,ityp,igauss,ng,scale,
     +            scs,scp,scd,odone)
      else if(odz) then
        call ahdzsh(nucz,ipass,opol,ityp,igauss,ng,scale,
     +            scs,scp,scd,odone)
      else if(otz) then
        call ahtzsh(nucz,ipass,opol,ityp,igauss,ng,scale,
     +            scs,scp,scd,odone)
      else
        call caserr2('unrecognised ahlrichs basis set')
      endif
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      if (nshell.gt.mxshel) ierr1 = 1
      if (ierr1.ne.0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss.gt.mxprim) ierr2 = 1
      if (ierr2.ne.0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      scale = scale*scale
      do 440 i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i) * scale
         if(ityp.eq.16) csinp(k) = ccs(ng+i)
         if(ityp.eq.17) cpinp(k) = ccp(ng+i)
         if(ityp.eq.18) cdinp(k) = ccd(ng+i)
         cs(k) = 0.0d+00
         cp(k) = 0.0d+00
         cd(k) = 0.0d+00
         cf(k) = 0.0d+00
         if(ityp.eq.16) cs(k) = csinp(k)
         if(ityp.eq.17) cp(k) = cpinp(k)
         if(ityp.eq.18) cd(k) = cdinp(k)
  440 continue
c
c     ----- always unnormalize primitives -----
c
      do 460 k = k1,k2
         ee = ex(k)+ex(k)
         facs = pi32/(ee*sqrt(ee))
         facp = pt5*facs/ee
         facd = pt75*facs/(ee*ee)
         if(ityp.eq.16) cs(k) = cs(k)/sqrt(facs)
         if(ityp.eq.17) cp(k) = cp(k)/sqrt(facp)
         if(ityp.eq.18) cd(k) = cd(k)/sqrt(facd)
  460 continue
c
c     ----- if(normf.eq.0) normalize basis functions -----
c
      if (normf .eq. 1) go to 210
      facs = 0.0d+00
      facp = 0.0d+00
      facd = 0.0d+00
      do 510 ig = k1,k2
         do 500 jg = k1,ig
            ee = ex(ig)+ex(jg)
            fac = ee*sqrt(ee)
            dums = cs(ig)*cs(jg)/fac
            dump = pt5*cp(ig)*cp(jg)/(ee*fac)
            dumd = pt75*cd(ig)*cd(jg)/(ee*ee*fac)
            if (ig .eq. jg) go to 480
               dums = dums+dums
               dump = dump+dump
               dumd = dumd+dumd
  480       continue
            facs = facs+dums
            facp = facp+dump
            facd = facd+dumd
  500    continue
  510 continue
c
      fac=0.0d+00
      if(ityp.eq.16.and. facs.gt.tm10) fac=1.0d+00/sqrt(facs*pi32)
      if(ityp.eq.17.and. facp.gt.tm10) fac=1.0d+00/sqrt(facp*pi32)
      if(ityp.eq.18.and. facd.gt.tm10) fac=1.0d+00/sqrt(facd*pi32)
c
      do 550 ig = k1,k2
         if(ityp.eq.16) cs(ig) = fac*cs(ig)
         if(ityp.eq.17) cp(ig) = fac*cp(ig)
         if(ityp.eq.18) cd(ig) = fac*cd(ig)
  550 continue
      go to 210
c
  220 continue
      return
      end
c
      subroutine ddzero(e,s,p,n)
c
c     ----- dunning's contraction of huzinaga's (4s) -----
c
c     ----- hydrogen (4s)/(2s) -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      go to (100,120),n
c
c     ----- h -----
c
  100 continue
      e(1) = 1.336d+01
      s(1) = 0.032828d+00
      e(2) = 2.013d+00
      s(2) = 0.231204d+00
      e(3) = 4.538d-01
      s(3) = 0.817226d+00
      e(4) = 1.233d-01
      s(4) = 1.000000d+00
      e(5) = 1.000000d+00
      p(5) = 1.000000d+00
      return
c
c     ----- he -----
c
  120 call caserr2('requested basis set not available for helium')
      return
c
      end
c
      subroutine dzone(e,s,p,d,n)
c
c     ----- dunning's DZ contraction of huzinaga's (9s,5p) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li -----
c
  100 continue
      e(1) = 9.213d+02
      s(1) = 0.001367d+00
      e(2) = 1.387d+02
      s(2) = 0.010425d+00
      e(3) = 3.194d+01
      s(3) = 0.049859d+00
      e(4) = 9.353d+00
      s(4) = 0.160701d+00
      e(5) = 3.158d+00
      s(5) = 0.344604d+00
      e(6) = 1.157d+00
      s(6) = 0.425197d+00
      e(7) = 4.446d-01
      s(7) = 0.169468d+00
      e(8) = 4.446d-01
      s(8) = -0.222311d+00
      e(9) = 7.666d-02
      s(9) = 1.116477d+00
      e(10) = 2.864d-02
      s(10) = 1.000000d+00
      e(11) = 1.488d+00
      p(11) = 0.038770d+00
      e(12) = 2.667d-01
      p(12) = 0.236257d+00
      e(13) = 7.201d-02
      p(13) = 0.830448d+00
      e(14) = 2.370d-02
      p(14) = 1.000000d+00
      e(15) = 2.0d-01
      d(15) = 1.000000d+00
c
      e(16) = 8.0d-03
      s(16) = 1.000000d+00
      e(17) = 6.0d-03
      p(17) = 1.000000d+00
      e(18) = 6.0d-03
      d(18) = 1.000000d+00
c
      return
c
c     ----- be -----
c
  120 continue
      e(1) = 1.741d+03
      s(1) = 0.001305d+00
      e(2) = 2.621d+02
      s(2) = 0.009955d+00
      e(3) = 6.033d+01
      s(3) = 0.048031d+00
      e(4) = 1.762d+01
      s(4) = 0.158577d+00
      e(5) = 5.933d+00
      s(5) = 0.351325d+00
      e(6) = 2.185d+00
      s(6) = 0.427006d+00
      e(7) = 8.590d-01
      s(7) = 0.160490d+00
      e(8) = 2.185d+00
      s(8) = -0.185294d+00
      e(9) = 1.806d-01
      s(9) = 1.057014d+00
      e(10) = 5.835d-02
      s(10) = 1.000000d+00
      e(11) = 6.710d+00
      p(11) = 0.016378d+00
      e(12) = 1.442d+00
      p(12) = 0.091553d+00
      e(13) = 4.103d-01
      p(13) = 0.341469d+00
      e(14) = 1.397d-01
      p(14) = 0.685428d+00
      e(15) = 4.922d-02
      p(15) = 1.000000d+00
      e(16) = 2.55d-01
      d(16) = 1.000000d+00
c
      e(17) = 1.6d-02
      s(17) = 1.000000d+00
      e(18) = 1.2d-02
      p(18) = 1.000000d+00
      e(19) = 1.2d-02
      d(19) = 1.000000d+00
c
      return
c
c     ----- b  -----
c
  140 continue
      e(1) = 2.78841d+03
      s(1) = 0.002122d+00
      e(2) = 4.19039d+02
      s(2) = 0.016171d+00
      e(3) = 9.64683d+01
      s(3) = 0.078356d+00
      e(4) = 2.80694d+01
      s(4) = 0.263250d+00
      e(5) = 9.3760d+00
      s(5) = 0.596729d+00
      e(6) = 1.3057d+00
      s(6) = 0.230397d+00
      e(7) = 3.4062d+00
      s(7) = 1.0d+00
      e(8) = 0.3245d+00
      s(8) = 1.d+00
      e(9) = 0.1022d+00
      s(9) = 1.d+00
      e(10) = 1.13413d+01
      p(10) = 0.017987d+00
      e(11) = 2.436d+00
      p(11) = 0.110339d+00
      e(12) = 6.836d-01
      p(12) = 0.383111d+00
      e(13) = 2.134d-01
      p(13) = 0.647860d+00
      e(14) = 7.01d-02
      p(14) = 1.000000d+00
      e(15) = 7.00d-01
      d(15) = 1.000000d+00
c
      e(16) = 1.9d-02
      s(16) = 1.000000d+00
      e(17) = 1.5d-02
      p(17) = 1.000000d+00
      e(18) = 1.5d-02
      d(18) = 1.000000d+00
c
      return
c
c     ----- c  -----
c
  160 continue
      e(1) = 4.23261d+03
      s(1) = 0.002029d+00
      e(2) = 6.34882d+02
      s(2) = 0.015535d+00
      e(3) = 1.46097d+02
      s(3) = 0.075411d0
      e(4) = 4.24974d+01
      s(4) = 0.257121d+00
      e(5) = 1.41892d+01
      s(5) = 0.596555d+00
      e(6) = 1.9666d+00
      s(6) = 0.242517d0
      e(7) = 5.1477d+00
      s(7) = 1.0d+00
      e(8) = 4.962d-01
      s(8) = 1.d+00
      e(9) = 1.533d-01
      s(9) = 1.000000d+00
      e(10) = 1.81557d+01
      p(10) = 0.018534d+00
      e(11) = 3.9864d+00
      p(11) = 0.115442d+00
      e(12) = 1.1429d+00
      p(12) = 0.386206d+00
      e(13) = 3.594d-01
      p(13) = 0.640089d+00
      e(14) = 1.146d-01
      p(14) = 1.000000d+00
      e(15) = 7.50d-01
      d(15) = 1.000000d+00
c
      e(16) = 2.3d-02
      s(16) = 1.000000d+00
      e(17) = 2.1d-02
      p(17) = 1.000000d+00
      e(18) = 1.5d-02
      d(18) = 1.000000d+00
c
      return
c
c     ----- n -----
c
  180 continue
      e(1) = 5.90944d+03
      s(1) = 0.002004d+00
      e(2) = 8.874510d+02
      s(2) = 0.015310d+00
      e(3) = 2.0479d+02
      s(3) = 0.074293d+00
      e(4) = 5.98376d+01
      s(4) = 0.253364d+00
      e(5) = 19.9981d+00
      s(5) = 0.600576d+00
      e(6) = 2.686d+00
      s(6) = 0.245111d+00
      e(7) = 7.1927d+00
      s(7) = 1.d+00
      e(8) = 7.000d-01
      s(8) = 1.d+00
      e(9) = 2.133d-01
      s(9) = 1.000000d+00
      e(10) = 2.6786d+01
      p(10) = 0.018257d+00
      e(11) = 5.9564d+00
      p(11) = 0.116407d+00
      e(12) = 1.7074d+00
      p(12) = 0.390111d+00
      e(13) = 5.314d-01
      p(13) = 0.637221d+00
      e(14) = 1.654d-01
      p(14) = 1.000000d+00
      e(15) = 8.00d-01
      d(15) = 1.000000d+00
c
      e(16) = 2.8d-02
      s(16) = 1.000000d+00
      e(17) = 2.5d-02
      p(17) = 1.000000d+00
      e(18) = 1.5d-02
      d(18) = 1.000000d+00
c
      return
c
c     ----- o  ------
c
  200 continue
      e(1) = 7.81654d+03
      s(1) = 0.002031d+00
      e(2) = 1.17582d+03
      s(2) = 0.015436d+00
      e(3) = 2.73188d+02
      s(3) = 0.073771d+00
      e(4) = 8.11696d+01
      s(4) = 0.247606d+00
      e(5) = 2.71836d+01
      s(5) = 0.611832d+00
      e(6) = 3.4136d+00
      s(6) = 0.241205d+00
      e(7) = 9.5322d+00
      s(7) = 1.d+00
      e(8) = 9.398d-01
      s(8) = 1.d+00
      e(9) = 2.846d-01
      s(9) = 1.000000d+00
      e(10) = 3.51832d+01
      p(10) = 0.019580d+00
      e(11) = 7.904d+00
      p(11) = 0.124189d+00
      e(12) = 2.3051d+00
      p(12) = 0.394727d+00
      e(13) = 7.171d-01
      p(13) = 0.627375d+00
      e(14) = 2.137d-01
      p(14) = 1.000000d+00
c
      e(15) = 8.50d-01
      d(15) = 1.000000d+00
c
      e(16) = 3.2d-02
      s(16) = 1.000000d+00
      e(17) = 2.8d-02
      p(17) = 1.000000d+00
      e(18) = 1.5d-02
      d(18) = 1.000000d+00
c
      return
c
c     ----- f  -----
c
  220 continue
      e(1) = 9.99479d+03
      s(1) = 0.002017d+00
      e(2) = 1.50603d+03
      s(2) = 0.015295d+00
      e(3) = 350.269d+00
      s(3) = 0.073110d+00
      e(4) = 1.04053d+02
      s(4) = 0.246420d+00
      e(5) = 3.48432d+01
      s(5) = 0.612593d+00
      e(6) = 4.3688d+00
      s(6) = 0.242489d+00
      e(7) = 1.22164d+01
      s(7) = 1.d+00
      e(8) = 1.2078d+00
      s(8) = 1.d+00
      e(9) = 3.634d-01
      s(9) = 1.000000d+00
      e(10) = 4.43555d+01
      p(10) = 0.020868d+00
      e(11) = 1.0082d+01
      p(11) = 0.130092d+00
      e(12) = 2.9959d+00
      p(12) = 0.396219d+00
      e(13) = 9.383d-01
      p(13) = 0.620368d+00
      e(14) = 2.733d-01
      p(14) = 1.000000d+00
      e(15) = 9.00d-01
      d(15) = 1.000000d+00
c
      e(16) = 3.6d-02
      s(16) = 1.000000d+00
      e(17) = 2.9d-02
      p(17) = 1.000000d+00
      e(18) = 1.5d-02
      d(18) = 1.000000d+00
c
      return
c
c     ----- ne -----
c
  240 continue
      e(1) = 1.210d+04
      s(1) = 0.001200d+00
      e(2) = 1.821d+03
      s(2) = 0.009092d+00
      e(3) = 4.328d+02
      s(3) = 0.041305d+00
      e(4) = 1.325d+02
      s(4) = 0.137867d+00
      e(5) = 4.377d+01
      s(5) = 0.362433d+00
      e(6) = 1.491d+01
      s(6) = 0.472247d+00
      e(7) = 5.127d+00
      s(7) = 0.130035d+00
      e(8) = 1.491d+01
      s(8) = -0.140810d+00
      e(9) = 1.491d+00
      s(9) = 1.053327d+00
      e(10) = 4.468d-01
      s(10) = 1.000000d+00
      e(11) = 5.645d+01
      p(11) = 0.020875d+00
      e(12) = 1.292d+01
      p(12) = 0.130032d+00
      e(13) = 3.865d+00
      p(13) = 0.395679d+00
      e(14) = 1.203d+00
      p(14) = 0.621450d+00
      e(15) = 3.444d-01
      p(15) = 1.000000d+00
      e(16) = 1.00d+00
      d(16) = 1.000000d+00
c
      e(17) = 4.0d-02
      s(17) = 1.000000d+00
      e(18) = 3.0d-02
      p(18) = 1.000000d+00
      e(19) = 1.5d-02
      d(19) = 1.000000d+00
c
      return
      end
c
      subroutine fzero(e,s,p,n)
c
c     ----- dunning's contraction of huzinaga's (5s) -----
c
c     ----- hydrogen (5s)/(3s) -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      go to (100,120),n
c
c     ----- h -----
c
  100 continue
      e(1) = 33.64d0
      s(1) = 0.025374d0
      e(2) = 5.058d0
      s(2) = 0.189684d0
      e(3) = 1.147d0
      s(3) = 0.852933d0
      e(4) = 0.3211d0
      s(4) = 1.0d0
      e(5) = 0.1013d0
      s(5) = 1.0d0
      e(6) = 1.0d0
      p(6) = 1.0d0
      return
c
c     ----- he -----
c
  120 call caserr2('requested basis set not available for helium')
      return
c
      end
c
      subroutine fone(e,s,p,d,n)
c
c     ----- dunning's contraction of huzinaga's (10s,6p) ----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li -----
c
 100  call caserr2('requested basis set not available for lithium')
      return
c
c     ----- be -----
c
  120 continue
      e(1) = 3630.0d0
      s(1) = 0.000839d0
      e(2) = 532.3d0
      s(2) = 0.006735d0
      e(3) = 117.8d0
      s(3) = 0.035726d0
      e(4) = 32.66d0
      s(4) = 0.138635d0
      e(5) = 10.48d0
      s(5) = 0.385399d0
      e(6) = 3.668d0
      s(6) = 0.547688d0
      e(7) = 3.668d0
      s(7) = 0.213406d0
      e(8) = 1.354d0
      s(8) = 0.814692d0
      e(9) = 0.3890d0
      s(9) = 1.0d0
      e(10) = 0.1502d0
      s(10) = 1.0d0
      e(11) = 0.05241d0
      s(11) = 1.0d0
      e(12) = 3.202d0
      p(12) = 0.052912d0
      e(13) = 0.6923d0
      p(13) = 0.267659d0
      e(14) = 0.2016d0
      p(14) = 0.792085d0
      e(15) = 0.1183d0
      p(15) = 1.0d0
      e(16) = 0.0694d0
      p(16) = 1.0d0
c
      e(17) = 0.32d0
      d(17) = 1.0d0
c
      e(18) = 1.6d-02
      s(18) = 1.000000d+00
      e(19) = 1.2d-02
      p(19) = 1.000000d+00
      e(20) = 1.2d-02
      d(20) = 1.000000d+00
c
      return
c
c     ----- b  -----
c
  140 continue
      e(1) = 6250.0d0
      s(1) = 0.000798d0
      e(2) = 916.1d0
      s(2) = 0.006410d0
      e(3) = 202.2d0
      s(3) = 0.034299d0
      e(4) = 55.83d0
      s(4) = 0.135487d0
      e(5) = 17.86d0
      s(5) = 0.388532d0
      e(6) = 6.253d0
      s(6) = 0.547758d0
      e(7) = 6.253d0
      s(7) = 0.232643d0
      e(8) = 2.312d0
      s(8) = 0.797219d0
      e(9) = 0.6824d0
      s(9) = 1.0d0
      e(10) = 0.2604d0
      s(10) = 1.0d0
      e(11) = 0.0894d0
      s(11) = 1.0d0
      e(12) = 15.46d0
      p(12) = 0.016822d0
      e(13) = 3.483d0
      p(13) = 0.100878d0
      e(14) = 1.066d0
      p(14) = 0.336895d0
      e(15) = 0.3928d0
      p(15) = 0.672317d0
      e(16) = 0.1503d0
      p(16) = 1.0d0
      e(17) = 0.05722d0
      p(17) = 1.0d0
c
      e(18) = 0.50d0
      d(18) = 1.0d0
c
      e(19) = 1.9d-02
      s(19) = 1.000000d+00
      e(20) = 1.5d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
c
c     ----- c  -----
c
  160 continue
      e(1) = 9471.0d0
      s(1) = 0.000776d0
      e(2) = 1398.0d0
      s(2) = 0.006218d0
      e(3) = 307.5d0
      s(3) = 0.033575d0
      e(4) = 84.54d0
      s(4) = 0.134278d0
      e(5) = 26.91d0
      s(5) = 0.393668d0
      e(6) = 9.409d0
      s(6) = 0.544169d0
      e(7) = 9.409d0
      s(7) = 0.248075d0
      e(8) = 3.500d0
      s(8) = 0.782844d0
      e(9) = 1.068d0
      s(9) = 1.0d0
      e(10) = 0.4002d0
      s(10) = 1.0d0
      e(11) = 0.1351d0
      s(11) = 1.0d0
      e(12) = 25.37d0
      p(12) = 0.016295d0
      e(13) = 5.776d0
      p(13) = 0.102098d0
      e(14) = 1.787d0
      p(14) = 0.340228d0
      e(15) = 0.6577d0
      p(15) = 0.668269d0
      e(16) = 0.2480d0
      p(16) = 1.0d0
      e(17) = 0.09106d0
      p(17) = 1.0d0
c
      e(18) = 0.72d0
      d(18) = 1.0d0
c
      e(19) = 2.3d-02
      s(19) = 1.000000d+00
      e(20) = 2.1d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
c
c     ----- n -----
c
  180 continue
      e(1) = 13520.0d0
      s(1) = 0.000760d0
      e(2) = 1999.0d0
      s(2) = 0.006076d0
      e(3) = 440.0d0
      s(3) = 0.032847d0
      e(4) = 120.9d0
      s(4) = 0.132396d0
      e(5) = 38.470d0
      s(5) = 0.393261d0
      e(6) = 13.46d0
      s(6) = 0.546339d0
      e(7) = 13.46d0
      s(7) = 0.252036d0
      e(8) = 4.993d0
      s(8) = 0.779385d0
      e(9) = 1.569d0
      s(9) = 1.0d0
      e(10) = 0.5800d0
      s(10) = 1.0d0
      e(11) = 0.1923d0
      s(11) = 1.0d0
      e(12) = 35.91d0
      p(12) = 0.016916d0
      e(13) = 8.480d0
      p(13) = 0.102200d0
      e(14) = 2.706d0
      p(14) = 0.338134d0
      e(15) = 0.9921d0
      p(15) = 0.669281d0
      e(16) = 0.3727d0
      p(16) = 1.0d0
      e(17) = 0.1346d0
      p(17) = 1.0d0
c
      e(18) = 0.98d0
      d(18) = 1.0d0
c
      e(19) = 2.8d-02
      s(19) = 1.000000d+00
      e(20) = 2.5d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
c
c     ----- o  ------
c
  200 continue
      e(1) = 18050.0d0
      s(1) = 0.000757d0
      e(2) = 2660.0d0
      s(2) = 0.006066d0
      e(3) = 585.7d0
      s(3) = 0.032782d0
      e(4) = 160.9d0
      s(4) = 0.132609d0
      e(5) = 51.160d0
      s(5) = 0.396839d0
      e(6) = 17.90d0
      s(6) = 0.542572d0
      e(7) = 17.90d0
      s(7) = 0.262490d0
      e(8) = 6.639d0
      s(8) = 0.769828d0
      e(9) = 2.077d0
      s(9) = 1.0d0
      e(10) = 0.7736d0
      s(10) = 1.0d0
      e(11) = 0.2558d0
      s(11) = 1.0d0
      e(12) = 49.83d0
      p(12) = 0.016358d0
      e(13) = 11.49d0
      p(13) = 0.106453d0
      e(14) = 3.609d0
      p(14) = 0.349302d0
      e(15) = 1.321d0
      p(15) = 0.657183d0
      e(16) = 0.4821d0
      p(16) = 1.0d0
      e(17) = 0.1651d0
      p(17) = 1.0d0
c
      e(18) = 1.28d0
      d(18) = 1.0d0
c
      e(19) = 3.2d-02
      s(19) = 1.000000d+00
      e(20) = 2.8d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
c
c     ----- f  -----
c
  220 continue
      e(1) = 23340.0d0
      s(1) = 0.000757d0
      e(2) = 3431.0d0
      s(2) = 0.006081d0
      e(3) = 757.7d0
      s(3) = 0.032636d0
      e(4) = 209.2d0
      s(4) = 0.131704d0
      e(5) = 66.73d0
      s(5) = 0.396240d0
      e(6) = 23.37d0
      s(6) = 0.543672d0
      e(7) = 23.37d0
      s(7) = 0.264893d0
      e(8) = 8.624d0
      s(8) = 0.767925d0
      e(9) = 2.692d0
      s(9) = 1.0d0
      e(10) = 1.009d0
      s(10) = 1.0d0
      e(11) = 0.3312d0
      s(11) = 1.0d0
      e(12) = 65.66d0
      p(12) = 0.016037d0
      e(13) = 15.22d0
      p(13) = 0.105697d0
      e(14) = 4.788d0
      p(14) = 0.350227d0
      e(15) = 1.732d0
      p(15) = 0.658195d0
      e(16) = 0.6206d0
      p(16) = 1.0d0
      e(17) = 0.2070d0
      p(17) = 1.0d0
c
      e(18) = 1.62d0
      d(18) = 1.0d0
c
      e(19) = 3.6d-02
      s(19) = 1.000000d+00
      e(20) = 2.9d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
c
c     ----- ne -----
c
  240 continue
      e(1) = 28660.0d0
      s(1) = 0.000767d0
      e(2) = 4263.0d0
      s(2) = 0.006068d0
      e(3) = 946.8d0
      s(3) = 0.032474d0
      e(4) = 261.5d0
      s(4) = 0.131468d0
      e(5) = 83.34d0
      s(5) = 0.397723d0
      e(6) = 29.17d0
      s(6) = 0.542491d0
      e(7) = 29.17d0
      s(7) = 0.269065d0
      e(8) = 10.76d0
      s(8) = 0.764121d0
      e(9) = 3.343d0
      s(9) = 1.0d0
      e(10) = 1.241d0
      s(10) = 1.0d0
      e(11) = 0.4063d0
      s(11) = 1.0d0
      e(12) = 84.84d0
      p(12) = 0.015550d0
      e(13) = 19.71d0
      p(13) = 0.103011d0
      e(14) = 6.219d0
      p(14) = 0.349215d0
      e(15) = 2.211d0
      p(15) = 0.662839d0
      e(16) = 0.7853d0
      p(16) = 1.0d0
      e(17) = 0.2566d0
      p(17) = 1.0d0
c
      e(18) = 2.0d0
      d(18) = 1.0d0
c
      e(19) = 4.0d-02
      s(19) = 1.000000d+00
      e(20) = 3.0d-02
      p(20) = 1.000000d+00
      e(21) = 1.5d-02
      d(21) = 1.000000d+00
c
      return
      end
      subroutine ftwo(e,s,p,d,n)
c
c     ----- mclean's  contraction of huzinaga's (12s,9p) -----
c           polarisation + diffuse functions (dunning)
c
      implicit real*8 (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- na ----
c
100   continue
      e(1) = 36166.4d0
      s(1) = 0.001032d0
      e(2) = 5372.58d0
      s(2) = 0.008071d0
      e(3) = 1213.21d0
      s(3) = 0.042129d0
      e(4) = 339.62d0
      s(4) = 0.169789d0
      e(5) = 109.55d0
      s(5) = 0.514621d0
      e(6) = 38.777d0
      s(6) = 0.379817d0
      e(7) = 38.777d0
      s(7) = 0.374762d0
      e(8) = 14.576d0
      s(8) = 0.575769d0
      e(9) = 5.2699d0
      s(9) = 0.112933d0
      e(10) = 1.8278d0
      s(10) = 1.0d0
      e(11) = 0.6199d0
      s(11) = 1.0d0
      e(12) = 0.0572d0
      s(12) = 1.0d0
      e(13) = 0.0240d0
      s(13) = 1.0d0
      e(14) = 144.645d0
      p(14) = 0.011485d0
      e(15) = 33.907d0
      p(15) = 0.082383d0
      e(16) = 10.629d0
      p(16) = 0.319658d0
      e(17) = 3.8239d0
      p(17) = 0.701295d0
      e(18) = 1.4443d0
      p(18) = 0.638506d0
      e(19) = 0.5526d0
      p(19) = 0.425365d0
      e(20) = 0.1887d0
      p(20) = 1.0d0
      e(21) = 0.0465d0
      p(21) = 1.0d0
      e(22) = 0.0163d0
      p(22) = 1.0d0
c
      e(23) = 0.157d0
      d(23) = 1.0d0
c
      e(24) = 8.0d-03
      s(24) = 1.000000d+00
      e(25) = 6.0d-03
      p(25) = 1.000000d+00
      e(26) = 6.0d-03
      d(26) = 1.000000d+00
      return
c
c     ----- mg ----
c
 120  continue
      e(1) = 43866.5d0
      s(1) = 0.000918d0
      e(2) = 6605.37d0
      s(2) = 0.007047d0
      e(3) = 1513.26d0
      s(3) = 0.035941d0
      e(4) = 432.317d0
      s(4) = 0.141461d0
      e(5) = 142.149d0
      s(5) = 0.426764d0
      e(6) = 51.398d0
      s(6) = 0.497975d0
      e(7) = 51.398d0
      s(7) = 0.251355d0
      e(8) = 19.920d0
      s(8) = 0.618671d0
      e(9) = 8.0247d0
      s(9) = 0.188417d0
      e(10) = 2.5082d0
      s(10) = 1.0d0
      e(11) = 0.8715d0
      s(11) = 1.0d0
      e(12) = 0.1082d0
      s(12) = 1.0d0
      e(13) = 0.0401d0
      s(13) = 1.0d0
      e(14) = 193.854d0
      p(14) = 0.010188d0
      e(15) = 45.442d0
      p(15) = 0.075360d0
      e(16) = 14.186d0
      p(16) = 0.307419d0
      e(17) = 5.0575d0
      p(17) = 0.717575d0
      e(18) = 1.8886d0
      p(18) = 0.667339d0
      e(19) = 0.7227d0
      p(19) = 0.394649d0
      e(20) = 0.2364d0
      p(20) = 1.0d0
      e(21) = 0.0934d0
      p(21) = 1.0d0
      e(22) = 0.0348d0
      p(22) = 1.0d0
c
      e(23) = 0.234d0
      d(23) = 1.0d0
c
      e(24) = 1.2d-02
      s(24) = 1.000000d+00
      e(25) = 9.0d-03
      p(25) = 1.000000d+00
      e(26) = 9.0d-03
      d(26) = 1.000000d+00
      return
c
c     ----- al -----
c
  140 continue
      e(1) = 54866.49d0
      s(1) = 0.000839d0
      e(2) = 8211.77d0
      s(2) = 0.006527d0
      e(3) = 1866.18d0
      s(3) = 0.033666d0
      e(4) = 531.129d0
      s(4) = 0.132902d0
      e(5) = 175.118d0
      s(5) = 0.401266d0
      e(6) = 64.0055d0
      s(6) = 0.531338d0
      e(7) = 64.0055d0
      s(7) = 0.202305d0
      e(8) = 25.2925d0
      s(8) = 0.624790d0
      e(9) = 10.5349d0
      s(9) = 0.227439d0
      e(10) = 3.2067d0
      s(10) = 1.0d0
      e(11) = 1.1526d0
      s(11) = 1.0d0
      e(12) = 0.1767d0
      s(12) = 1.0d0
      e(13) = 0.0652d0
      s(13) = 1.0d0
      e(14) = 259.284d0
      p(14) = 0.009448d0
      e(15) = 61.0769d0
      p(15) = 0.070974d0
      e(16) = 19.3032d0
      p(16) = 0.295636d0
      e(17) = 7.0109d0
      p(17) = 0.728219d0
      e(18) = 2.6739d0
      p(18) = 0.644467d0
      e(19) = 1.0366d0
      p(19) = 0.417413d0
      e(20) = 0.3168d0
      p(20) = 1.0d0
      e(21) = 0.1143d0
      p(21) = 1.0d0
      e(22) = 0.0414d0
      p(22) = 1.0d0
c
      e(23) = 0.311d0
      d(23) = 1.0d0
c
      e(24) = 1.4d-02
      s(24) = 1.000000d+00
      e(25) = 1.1d-02
      p(25) = 1.000000d+00
      e(26) = 1.5d-02
      d(26) = 1.000000d+00
      return
c
c     ----- si -----
c
  160 continue
      e(1) = 69379.23d0
      s(1) = 0.000757d0
      e(2) = 10354.94d0
      s(2) = 0.005932d0
      e(3) = 2333.88d0
      s(3) = 0.031088d0
      e(4) = 657.14d0
      s(4) = 0.124967d0
      e(5) = 214.30d0
      s(5) = 0.386897d0
      e(6) = 77.629d0
      s(6) = 0.554888d0
      e(7) = 77.629d0
      s(7) = 0.177881d0
      e(8) = 30.631d0
      s(8) = 0.627765d0
      e(9) = 12.801d0
      s(9) = 0.247623d0
      e(10) = 3.9269d0
      s(10) = 1.0d0
      e(11) = 1.4523d0
      s(11) = 1.0d0
      e(12) = 0.2562d0
      s(12) = 1.0d0
      e(13) = 0.0943d0
      s(13) = 1.0d0
      e(14) = 335.48d0
      p(14) = 0.008866d0
      e(15) = 78.90d0
      p(15) = 0.068299d0
      e(16) = 24.988d0
      p(16) = 0.290958d0
      e(17) = 9.2197d0
      p(17) = 0.732116d0
      e(18) = 3.6211d0
      p(18) = 0.619879d0
      e(19) = 1.4513d0
      p(19) = 0.439148d0
      e(20) = 0.5050d0
      p(20) = 1.0d0
      e(21) = 0.1863d0
      p(21) = 1.0d0
      e(22) = 0.0654d0
      p(22) = 1.0d0
c
      e(23) = 0.388d0
      d(23) = 1.0d0
c
      e(24) = 1.7d-02
      s(24) = 1.000000d+00
      e(25) = 1.4d-02
      p(25) = 1.000000d+00
      e(26) = 1.5d-02
      d(26) = 1.000000d+00
      return
c
c     ----- p  -----
c
  180 continue
      e(1) = 77492.43d0
      s(1) = 0.000787d0
      e(2) = 11605.79d0
      s(2) = 0.006107d0
      e(3) = 2645.96d0
      s(3) = 0.031373d0
      e(4) = 754.98d0
      s(4) = 0.124239d0
      e(5) = 248.75d0
      s(5) = 0.380893d0
      e(6) = 91.157d0
      s(6) = 0.559812d0
      e(7) = 91.157d0
      s(7) = 0.163986d0
      e(8) = 36.226d0
      s(8) = 0.625950d0
      e(9) = 15.211d0
      s(9) = 0.262196d0
      e(10) = 4.7138d0
      s(10) = 1.0d0
      e(11) = 1.7827d0
      s(11) = 1.0d0
      e(12) = 0.3425d0
      s(12) = 1.0d0
      e(13) = 0.1246d0
      s(13) = 1.0d0
      e(14) = 384.84d0
      p(14) = 0.003240d0
      e(15) = 90.552d0
      p(15) = 0.024925d0
      e(16) = 28.806d0
      p(16) = 0.105697d0
      e(17) = 10.688d0
      p(17) = 0.263229d0
      e(18) = 4.2521d0
      p(18) = 0.719053d0
      e(19) = 4.2521d0
      p(19) = -1.612739d0
      e(20) = 1.7405d0
      p(20) = 1.205083d0
      e(21) = 0.5979d0
      p(21) = 1.0d0
      e(22) = 0.2292d0
      p(22) = 1.0d0
      e(23) = 0.0838d0
      p(23) = 1.0d0
c
      e(24) = 0.465d0
      d(24) = 1.0d0
c
      e(25) = 2.0d-02
      s(25) = 1.000000d+00
      e(26) = 1.7d-02
      p(26) = 1.000000d+00
      e(27) = 1.5d-02
      d(27) = 1.000000d+00
      return
c
c     ----- s  -----
c
  200 continue
      e(1) = 93413.4d0
      s(1) = 0.000742d0
      e(2) = 13961.7d0
      s(2) = 0.005790d0
      e(3) = 3169.9d0
      s(3) = 0.029945d0
      e(4) = 902.46d0
      s(4) = 0.118971d0
      e(5) = 297.16d0
      s(5) = 0.368273d0
      e(6) = 108.702d0
      s(6) = 0.577507d0
      e(7) = 108.702d0
      s(7) = 0.142943d0
      e(8) = 43.155d0
      s(8) = 0.624606d0
      e(9) = 18.108d0
      s(9) = 0.283438d0
      e(10) = 5.5705d0
      s(10) = 1.0d0
      e(11) = 2.1427d0
      s(11) = 1.0d0
      e(12) = 0.4340d0
      s(12) = 1.0d0
      e(13) = 0.1570d0
      s(13) = 1.0d0
      e(14) = 495.04d0
      p(14) = 0.002720d0
      e(15) = 117.22d0
      p(15) = 0.021090d0
      e(16) = 37.507d0
      p(16) = 0.092369d0
      e(17) = 13.910d0
      p(17) = 0.246781d0
      e(18) = 5.5045d0
      p(18) = 0.743840d0
      e(19) = 5.5045d0
      p(19) = -1.608719d0
      e(20) = 2.2433d0
      p(20) = 1.223255d0
      e(21) = 0.7762d0
      p(21) = 1.0d0
      e(22) = 0.2919d0
      p(22) = 1.0d0
      e(23) = 0.1029d0
      p(23) = 1.0d0
c
      e(24) = 0.542d0
      d(24) = 1.0d0
c
      e(25) = 2.3d-02
      s(25) = 1.000000d+00
      e(26) = 2.0d-02
      p(26) = 1.000000d+00
      e(27) = 1.5d-02
      d(27) = 1.000000d+00
      return
c
c     ----- cl ----
c
 220  continue
      e(1) = 105818.8d0
      s(1) = 0.000743d0
      e(2) = 15872.0d0
      s(2) = 0.005753d0
      e(3) = 3619.7d0
      s(3) = 0.029676d0
      e(4) = 1030.8d0
      s(4) = 0.118010d0
      e(5) = 339.91d0
      s(5) = 0.365230d0
      e(6) = 124.538d0
      s(6) = 0.581221d0
      e(7) = 124.538d0
      s(7) = 0.137548d0
      e(8) = 49.514d0
      s(8) = 0.622881d0
      e(9) = 20.806d0
      s(9) = 0.290143d0
      e(10) = 6.4648d0
      s(10) = 1.0d0
      e(11) = 2.5254d0
      s(11) = 1.0d0
      e(12) = 0.5378d0
      s(12) = 1.0d0
      e(13) = 0.1935d0
      s(13) = 1.0d0
      e(14) = 589.78d0
      p(14) = 0.002760d0
      e(15) = 139.85d0
      p(15) = 0.021536d0
      e(16) = 44.795d0
      p(16) = 0.095916d0
      e(17) = 16.612d0
      p(17) = 0.262315d0
      e(18) = 6.5995d0
      p(18) = 0.726811d0
      e(19) = 6.5995d0
      p(19) = -1.564294d0
      e(20) = 2.7141d0
      p(20) = 1.495778d0
      e(21) = 0.9528d0
      p(21) = 1.0d0
      e(22) = 0.3580d0
      p(22) = 1.0d0
      e(23) = 0.1250d0
      p(23) = 1.0d0
c
      e(24) = 0.619d0
      d(24) = 1.0d0
c
      e(25) = 2.5d-02
      s(25) = 1.000000d+00
      e(26) = 2.0d-02
      p(26) = 1.000000d+00
      e(27) = 1.5d-02
      d(27) = 1.000000d+00
      return
c
c     ----- ar ----
c
 240  continue
      e(1) = 118022.4d0
      s(1) = 0.000747d0
      e(2) = 17683.5d0
      s(2) = 0.005790d0
      e(3) = 4027.8d0
      s(3) = 0.029919d0
      e(4) = 1145.40d0
      s(4) = 0.119206d0
      e(5) = 377.16d0
      s(5) = 0.369028d0
      e(6) = 138.160d0
      s(6) = 0.576459d0
      e(7) = 138.160d0
      s(7) = 0.143927d0
      e(8) = 54.989d0
      s(8) = 0.622938d0
      e(9) = 23.171d0
      s(9) = 0.283964d0
      e(10) = 7.3779d0
      s(10) = 1.0d0
      e(11) = 2.9237d0
      s(11) = 1.0d0
      e(12) = 0.6504d0
      s(12) = 1.0d0
      e(13) = 0.2328d0
      s(13) = 1.0d0
      e(14) = 663.06d0
      p(14) = 0.003082d0
      e(15) = 157.09d0
      p(15) = 0.024165d0
      e(16) = 50.231d0
      p(16) = 0.108223d0
      e(17) = 18.635d0
      p(17) = 0.294192d0
      e(18) = 7.4465d0
      p(18) = 0.687862d0
      e(19) = 7.4465d0
      p(19) = -1.214482d0
      e(20) = 3.0957d0
      p(20) = 1.632370d0
      e(21) = 1.1065d0
      p(21) = 1.0d0
      e(22) = 0.4156d0
      p(22) = 1.0d0
      e(23) = 0.1454d0
      p(23) = 1.0d0
c
      e(24) = 0.696d0
      d(24) = 1.0d0
c
      e(25) = 2.7d-02
      s(25) = 1.000000d+00
      e(26) = 2.0d-02
      p(26) = 1.000000d+00
      e(27) = 1.5d-02
      d(27) = 1.000000d+00
c
      return
      end
c
      subroutine fthree(e,s,p,d,n)
c
c     ----- <10s8p3d> contraction for k-zn .....
c
c     wachters 14s,9p,5d modified as follows ...
c         1. most diffuse s removed ..
c         2. additional s spanning the 3s-4s region
c         3. 2 additional p's to describe 4p ...
c         4. (5d/4d,1d) + hay diffuse d function
c
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-18
      go to (100,120,130,140,150,160,170,180,190,
     *       200,210,220), nn
c
c     ----- k  ----
c
100   call caserr2('requested basis set not available for potassium')
      return
c
c     ----- ca ----
c
 120  call caserr2('requested basis set not available for calcium')
      return
c
c ----- sc
c
 130  e(1)=188961.0d0
      s(1)=0.00025d0
      e(2)=28491.4d0
      s(2)=0.00193d0
      e(3)=6572.43d0
      s(3)=0.00980d0
      e(4)=1881.03d0
      s(4)=0.03940d0
      e(5)=617.979d0
      s(5)=0.12532d0
      e(6)=225.127d0
      s(6)=1.0d0
      e(7)=88.5664d0
      s(7)=1.0d0
      e(8)=36.5819d0
      s(8)=1.0d0
      e(9)=11.0358d0
      s(9)=1.0d0
      e(10)=4.49063d0
      s(10)=1.0d0
      e(11)=1.12935d0
      s(11)=1.0d0
      e(12)=0.454613d0
      s(12)=1.0d0
      e(13)=0.188d0
      s(13)=1.0d0
      e(14)=0.077533d0
      s(14)=1.0d0
      e(15)=1113.82d0
      p(15)=0.00223d0
      e(16)=266.244d0
      p(16)=0.01764d0
      e(17)=86.5763d0
      p(17)=0.08123d0
      e(18)=32.5934d0
      p(18)=0.24040d0
      e(19)=13.2190d0
      p(19)=1.0d0
      e(20)=5.58357d0
      p(20)=1.0d0
      e(21)=2.18594d0
      p(21)=1.0d0
      e(22)=0.895011d0
      p(22)=1.0d0
      e(23)=0.351975d0
      p(23)=1.0d0
      e(24)=0.137d0
      p(24)=1.0d0
      e(25)=0.053d0
      p(25)=1.0d0
      e(26)=22.4283d0
      d(26)=0.02140d0
      e(27)=6.03479d0
      d(27)=0.11304d0
      e(28)=1.97905d0
      d(28)=0.31393d0
      e(29)=0.666750d0
      d(29)=0.45818d0
      e(30)=0.218231d0
      d(30)=1.0d0
      e(31)=0.0588d0
      d(31)=1.0d0
      return
c
c ----- ti
c
 140  e(1)=206082.0d0
      s(1)=0.00025d0
      e(2)=31226.8d0
      s(2)=0.00193d0
      e(3)=7199.32d0
      s(3)=0.00988d0
      e(4)=2048.75d0
      s(4)=0.03994d0
      e(5)=670.790d0
      s(5)=0.12737d0
      e(6)=243.650d0
      s(6)=1.0d0
      e(7)=95.9250d0
      s(7)=1.0d0
      e(8)=39.8101d0
      s(8)=1.0d0
      e(9)=12.2205d0
      s(9)=1.0d0
      e(10)=5.00882d0
      s(10)=1.0d0
      e(11)=1.28569d0
      s(11)=1.0d0
      e(12)=0.512806d0
      s(12)=1.0d0
      e(13)=0.209d0
      s(13)=1.0d0
      e(14)=0.085576d0
      s(14)=1.0d0
      e(15)=1264.70d0
      p(15)=0.00214d0
      e(16)=301.230d0
      p(16)=0.01725d0
      e(17)=96.9777d0
      p(17)=0.08117d0
      e(18)=36.3727d0
      p(18)=0.24143d0
      e(19)=14.7814d0
      p(19)=1.0d0
      e(20)=6.27465d0
      p(20)=1.0d0
      e(21)=2.47878d0
      p(21)=1.0d0
      e(22)=1.01618d0
      p(22)=1.0d0
      e(23)=0.398162d0
      p(23)=1.0d0
      e(24)=0.156d0
      p(24)=1.0d0
      e(25)=0.0611d0
      p(25)=1.0d0
      e(26)=25.9924d0
      d(26)=0.02253d0
      e(27)=7.08634d0
      d(27)=0.11961d0
      e(28)=2.34871d0
      d(28)=0.32700d0
      e(29)=0.800198d0
      d(29)=0.46626d0
      e(30)=0.262049d0
      d(30)=1.0d0
      e(31)=0.0720d0
      d(31)=1.0d0
      return
c
c ----- v
c
 150  e(1)=226878.0d0
      s(1)=0.00025d0
      e(2)=33899.6d0
      s(2)=0.00197d0
      e(3)=7720.24d0
      s(3)=0.01017d0
      e(4)=2191.53d0
      s(4)=0.04098d0
      e(5)=719.169d0
      s(5)=0.12980d0
      e(6)=262.086d0
      s(6)=1.0d0
      e(7)=103.653d0
      s(7)=1.0d0
      e(8)=43.2548d0
      s(8)=1.0d0
      e(9)=13.5088d0
      s(9)=1.0d0
      e(10)=5.56715d0
      s(10)=1.0d0
      e(11)=1.45174d0
      s(11)=1.0d0
      e(12)=0.574445d0
      s(12)=1.0d0
      e(13)=0.231d0
      s(13)=1.0d0
      e(14)=0.092820d0
      s(14)=1.0d0
      e(15)=1398.43d0
      p(15)=0.00214d0
      e(16)=331.571d0
      p(16)=0.01735d0
      e(17)=107.002d0
      p(17)=0.08128d0
      e(18)=40.3183d0
      p(18)=0.24130d0
      e(19)=16.4635d0
      p(19)=1.0d0
      e(20)=7.02352d0
      p(20)=1.0d0
      e(21)=2.79025d0
      p(21)=1.0d0
      e(22)=1.14609d0
      p(22)=1.0d0
      e(23)=0.447272d0
      p(23)=1.0d0
      e(24)=0.175d0
      p(24)=1.0d0
      e(25)=0.068d0
      p(25)=1.0d0
      e(26)=30.2212d0
      d(26)=0.02270d0
      e(27)=8.27218d0
      d(27)=0.12248d0
      e(28)=2.75823d0
      d(28)=0.33286d0
      e(29)=0.942137d0
      d(29)=0.46989d0
      e(30)=0.305298d0
      d(30)=1.0d0
      e(31)=0.0820d0
      d(31)=1.0d0
      return
c
c ----- cr
c
 160  e(1)=236658.0d0
      s(1)=0.00027d0
      e(2)=35364.0d0
      s(2)=0.00208d0
      e(3)=8058.31d0
      s(3)=0.01071d0
      e(4)=2294.23d0
      s(4)=0.04284d0
      e(5)=756.118d0
      s(5)=0.13420d0
      e(6)=277.000d0
      s(6)=1.0d0
      e(7)=110.179d0
      s(7)=1.0d0
      e(8)=46.3710d0
      s(8)=1.0d0
      e(9)=14.8215d0
      s(9)=1.0d0
      e(10)=6.13262d0
      s(10)=1.0d0
      e(11)=1.62959d0
      s(11)=1.0d0
      e(12)=0.641177d0
      s(12)=1.0d0
      e(13)=0.253d0
      s(13)=1.0d0
      e(14)=0.099511d0
      s(14)=1.0d0
      e(15)=1478.77d0
      p(15)=0.00228d0
      e(16)=351.490d0
      p(16)=0.01836d0
      e(17)=113.826d0
      p(17)=0.08503d0
      e(18)=43.1567d0
      p(18)=0.24761d0
      e(19)=17.7775d0
      p(19)=1.0d0
      e(20)=7.66128d0
      p(20)=1.0d0
      e(21)=3.07765d0
      p(21)=1.0d0
      e(22)=1.26619d0
      p(22)=1.0d0
      e(23)=0.493534d0
      p(23)=1.0d0
      e(24)=0.192d0
      p(24)=1.0d0
      e(25)=0.075d0
      p(25)=1.0d0
      e(26)=34.0221d0
      d(26)=0.02328d0
      e(27)=9.43161d0
      d(27)=0.12564d0
      e(28)=3.15918d0
      d(28)=0.33897d0
      e(29)=1.07995d0
      d(29)=0.47139d0
      e(30)=0.346582d0
      d(30)=1.0d0
      e(31)=0.0912d0
      d(31)=1.0d0
      return
c
c ----- mn
c
 170  e(1)=243694.0d0
      s(1)=0.00029d0
      e(2)=35995.0d0
      s(2)=0.00225d0
      e(3)=8223.56d0
      s(3)=0.01152d0
      e(4)=2353.12d0
      s(4)=0.04559d0
      e(5)=780.965d0
      s(5)=0.14039d0
      e(6)=288.519d0
      s(6)=1.0d0
      e(7)=115.701d0
      s(7)=1.0d0
      e(8)=49.1175d0
      s(8)=1.0d0
      e(9)=16.0885d0
      s(9)=1.0d0
      e(10)=6.70430d0
      s(10)=1.0d0
      e(11)=1.80517d0
      s(11)=1.0d0
      e(12)=0.703011d0
      s(12)=1.0d0
      e(13)=0.273d0
      s(13)=1.0d0
      e(14)=0.106385d0
      s(14)=1.0d0
      e(15)=1500.39d0
      p(15)=0.00258d0
      e(16)=358.800d0
      p(16)=0.02048d0
      e(17)=116.699d0
      p(17)=0.09293d0
      e(18)=44.6132d0
      p(18)=0.26074d0
      e(19)=18.5985d0
      p(19)=1.0d0
      e(20)=8.13778d0
      p(20)=1.0d0
      e(21)=3.33734d0
      p(21)=1.0d0
      e(22)=1.37895d0
      p(22)=1.0d0
      e(23)=0.538639d0
      p(23)=1.0d0
      e(24)=0.210d0
      p(24)=1.0d0
      e(25)=0.082d0
      p(25)=1.0d0
      e(26)=37.8977d0
      d(26)=0.02409d0
      e(27)=10.5201d0
      d(27)=0.13073d0
      e(28)=3.53764d0
      d(28)=0.34648d0
      e(29)=1.21217d0
      d(29)=0.47163d0
      e(30)=0.387912d0
      d(30)=1.0d0
      e(31)=0.1054d0
      d(31)=1.0d0
      return
c
c ----- fe
c
 180  e(1)=257539.0d0
      s(1)=0.00029d0
      e(2)=38636.9d0
      s(2)=0.00226d0
      e(3)=8891.44d0
      s(3)=0.01152d0
      e(4)=2544.01d0
      s(4)=0.04566d0
      e(5)=844.777d0
      s(5)=0.14035d0
      e(6)=312.527d0
      s(6)=1.0d0
      e(7)=125.593d0
      s(7)=1.0d0
      e(8)=53.4987d0
      s(8)=1.0d0
      e(9)=17.7151d0
      s(9)=1.0d0
      e(10)=7.37677d0
      s(10)=1.0d0
      e(11)=2.01847d0
      s(11)=1.0d0
      e(12)=0.779935d0
      s(12)=1.0d0
      e(13)=0.298d0
      s(13)=1.0d0
      e(14)=0.114220d0
      s(14)=1.0d0
      e(15)=1678.40d0
      p(15)=0.00249d0
      e(16)=396.392d0
      p(16)=0.02015d0
      e(17)=128.588d0
      p(17)=0.09199d0
      e(18)=49.1158d0
      p(18)=0.25991d0
      e(19)=20.5035d0
      p(19)=1.0d0
      e(20)=8.98712d0
      p(20)=1.0d0
      e(21)=3.68249d0
      p(21)=1.0d0
      e(22)=1.52175d0
      p(22)=1.0d0
      e(23)=0.592684d0
      p(23)=1.0d0
      e(24)=0.231d0
      p(24)=1.0d0
      e(25)=0.0899d0
      p(25)=1.0d0
      e(26)=41.4526d0
      d(26)=0.02511d0
      e(27)=11.5403d0
      d(27)=0.13626d0
      e(28)=3.88543d0
      d(28)=0.35323d0
      e(29)=1.32380d0
      d(29)=0.46867d0
      e(30)=0.416680d0
      d(30)=1.0d0
      e(31)=0.1133d0
      d(31)=1.0d0
      return
c
c ----- co
c
 190  e(1)=270991.0d0
      s(1)=0.00031d0
      e(2)=39734.8d0
      s(2)=0.00242d0
      e(3)=9057.46d0
      s(3)=0.01238d0
      e(4)=2598.21d0
      s(4)=0.04849d0
      e(5)=868.200d0
      s(5)=0.14672d0
      e(6)=323.431d0
      s(6)=1.0d0
      e(7)=130.860d0
      s(7)=1.0d0
      e(8)=56.1219d0
      s(8)=1.0d0
      e(9)=18.9219d0
      s(9)=1.0d0
      e(10)=7.95238d0
      s(10)=1.0d0
      e(11)=2.19754d0
      s(11)=1.0d0
      e(12)=0.846713d0
      s(12)=1.0d0
      e(13)=0.321d0
      s(13)=1.0d0
      e(14)=0.122266d0
      s(14)=1.0d0
      e(15)=1636.21d0
      p(15)=0.00296d0
      e(16)=390.903d0
      p(16)=0.02336d0
      e(17)=127.884d0
      p(17)=0.10343d0
      e(18)=49.2413d0
      p(18)=0.27954d0
      e(19)=20.7512d0
      p(19)=1.0d0
      e(20)=9.20368d0
      p(20)=1.0d0
      e(21)=3.81779d0
      p(21)=1.0d0
      e(22)=1.58762d0
      p(22)=1.0d0
      e(23)=0.624660d0
      p(23)=1.0d0
      e(24)=0.2458d0
      p(24)=1.0d0
      e(25)=0.0967d0
      p(25)=1.0d0
      e(26)=44.9774d0
      d(26)=0.02627d0
      e(27)=12.5690d0
      d(27)=0.14182d0
      e(28)=4.24422d0
      d(28)=0.35942d0
      e(29)=1.44330d0
      d(29)=0.46613d0
      e(30)=0.449965d0
      d(30)=1.0d0
      e(31)=0.1219d0
      d(31)=1.0d0
      return
c
c ----- ni
c
 200  e(1)=284878.0d0
      s(1)=0.00032d0
      e(2)=41997.9d0
      s(2)=0.00246d0
      e(3)=9627.67d0
      s(3)=0.01254d0
      e(4)=2761.96d0
      s(4)=0.04926d0
      e(5)=920.488d0
      s(5)=0.14950d0
      e(6)=341.805d0
      s(6)=1.0d0
      e(7)=138.023d0
      s(7)=1.0d0
      e(8)=59.2587d0
      s(8)=1.0d0
      e(9)=20.3712d0
      s(9)=1.0d0
      e(10)=8.5940d0
      s(10)=1.0d0
      e(11)=2.39417d0
      s(11)=1.0d0
      e(12)=0.918169d0
      s(12)=1.0d0
      e(13)=0.346d0
      s(13)=1.0d0
      e(14)=0.130176d0
      s(14)=1.0d0
      e(15)=1774.18d0
      p(15)=0.00295d0
      e(16)=423.403d0
      p(16)=0.02337d0
      e(17)=138.311d0
      p(17)=0.10406d0
      e(18)=53.1703d0
      p(18)=0.28226d0
      e(19)=22.3874d0
      p(19)=1.0d0
      e(20)=9.92848d0
      p(20)=1.0d0
      e(21)=4.11625d0
      p(21)=1.0d0
      e(22)=1.71031d0
      p(22)=1.0d0
      e(23)=0.672528d0
      p(23)=1.0d0
      e(24)=0.264d0
      p(24)=1.0d0
      e(25)=0.104d0
      p(25)=1.0d0
      e(26)=48.9403d0
      d(26)=0.02706d0
      e(27)=13.7169d0
      d(27)=0.14598d0
      e(28)=4.63951d0
      d(28)=0.36418d0
      e(29)=1.57433d0
      d(29)=0.46348d0
      e(30)=0.486409d0
      d(30)=1.0d0
      e(31)=0.1316d0
      d(31)=1.0d0
      return
c
c ----- cu (2d)
c
 210  e(1)=307637.0d0
      s(1)=0.00031d0
      e(2)=46592.9d0
      s(2)=0.00236d0
      e(3)=10651.1d0
      s(3)=0.01213d0
      e(4)=3043.31d0
      s(4)=0.04791d0
      e(5)=1010.62d0
      s(5)=0.14652d0
      e(6)=374.252d0
      s(6)=1.0d0
      e(7)=150.796d0
      s(7)=1.0d0
      e(8)=64.6268d0
      s(8)=1.0d0
      e(9)=22.1381d0
      s(9)=1.0d0
      e(10)=9.34746d0
      s(10)=1.0d0
      e(11)=2.60863d0
      s(11)=1.0d0
      e(12)=0.997212d0
      s(12)=1.0d0
      e(13)=0.374d0
      s(13)=1.0d0
      e(14)=0.140120d0
      s(14)=1.0d0
      e(15)=2026.05d0
      p(15)=0.00267d0
      e(16)=484.207d0
      p(16)=0.02133d0
      e(17)=157.803d0
      p(17)=0.09690d0
      e(18)=60.4844d0
      p(18)=0.27103d0
      e(19)=25.3593d0
      p(19)=1.0d0
      e(20)=11.1685d0
      p(20)=1.0d0
      e(21)=4.56411d0
      p(21)=1.0d0
      e(22)=1.88440d0
      p(22)=1.0d0
      e(23)=0.734735d0
      p(23)=1.0d0
      e(24)=0.286d0
      p(24)=1.0d0
      e(25)=0.111d0
      p(25)=1.0d0
      e(26)=53.6478d0
      d(26)=0.02727d0
      e(27)=15.0747d0
      d(27)=0.14780d0
      e(28)=5.10392d0
      d(28)=0.36659d0
      e(29)=1.72743d0
      d(29)=0.46395d0
      e(30)=0.528322d0
      d(30)=1.0d0
      e(31)=0.1491d0
      d(31)=1.0d0
      return
c
c ----- zn (1s)
c
 220  e(1)=316336.0d0
      s(1)=0.00032d0
      e(2)=48561.0d0
      s(2)=0.00242d0
      e(3)=11157.4d0
      s(3)=0.01241d0
      e(4)=3205.01d0
      s(4)=0.04864d0
      e(5)=1068.58d0
      s(5)=0.14807d0
      e(6)=396.394d0
      s(6)=1.0d0
      e(7)=159.806d0
      s(7)=1.0d0
      e(8)=68.5890d0
      s(8)=1.0d0
      e(9)=23.7081d0
      s(9)=1.0d0
      e(10)=10.0372d0
      s(10)=1.0d0
      e(11)=2.81043d0
      s(11)=1.0d0
      e(12)=1.06964d0
      s(12)=1.0d0
      e(13)=0.396d0
      s(13)=1.0d0
      e(14)=0.146951d0
      s(14)=1.0d0
      e(15)=2213.18d0
      p(15)=0.00262d0
      e(16)=527.050d0
      p(16)=0.02091d0
      e(17)=172.293d0
      p(17)=0.09501d0
      e(18)=66.0814d0
      p(18)=0.26855d0
      e(19)=27.6863d0
      p(19)=1.0d0
      e(20)=12.1841d0
      p(20)=1.0d0
      e(21)=4.98796d0
      p(21)=1.0d0
      e(22)=2.05791d0
      p(22)=1.0d0
      e(23)=0.798609d0
      p(23)=1.0d0
      e(24)=0.311d0
      p(24)=1.0d0
      e(25)=0.121d0
      p(25)=1.0d0
      e(26)=58.4084d0
      d(26)=0.02759d0
      e(27)=16.4492d0
      d(27)=0.14994d0
      e(28)=5.57570d0
      d(28)=0.36929d0
      e(29)=1.88441d0
      d(29)=0.46343d0
      e(30)=0.572305d0
      d(30)=1.0d0
      e(31)=0.155d0
      d(31)=1.0d0
      return
c
      end
      subroutine poptzv(csinp,cpinp,cdinp,opol,oryd,
     +    nucz,intyp,nangm,nbfs,minf,maxf,loc,ngauss,ns,
     +    nshmax,ngsmax,ierr1,ierr2,nat,iwr)
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
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     *ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      dimension csinp(*),cpinp(*),cdinp(*)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      data pt75/0.75d0/
      data pt5,pi32/0.5d+00,5.56832799683170d+00/
c     data tol/1.0d-10/
c
      if(oryd) call caserr2('no diffuse functions available')
      call vclr(eex,1,200)
c
      if (nucz .le. 2) then
c
c     ----- hydrogen and helium
c
      mxpass=3
      if(opol) mxpass=mxpass+1
c
      else if (nucz .le. 10) then
c
c     ----- lithium to neon -----
c
      mxpass=4
      if(opol) mxpass=mxpass+1
c
      else if (nucz .le. 18) then
c
c     ----- sodium to argon ------------
c
      mxpass=11
      if(opol) mxpass=mxpass+1
c
       else
          if (opg_root()) then
             write(iwr,*)'*** nuclear charge = ',nucz
          endif
          call caserr ('requested basis set not available')
       endif
c
      call f6311(eex,ccs,ccp,ccd,nucz)
c
c     ----- loop over shells -----
c
      ipass = 0
  180 ipass = ipass+1
      if (nucz .gt. 2) go to 240
c
c     ----- h and he -----
c
      go to (200,220,230,235),ipass
  200 ityp = 1
      igauss = 3
      ig = 0
      go to 700
  220 ityp = 1
      igauss = 1
      ig = 3
      go to 700
  230 ityp = 1
      igauss = 1
      ig = 4
      go to 700
  235 ityp = 17
      igauss = 1
      ig = 5
      go to 700
  240 if (nucz .gt. 10) go to 480
c
c     ----- li  - to - ne -----
c           4(sp) or 5(spd) shells
c
      go to (260,280,300,320,340),ipass
  260 ityp = 1
      igauss = 6
      ig = 0
      go to 700
  280 ityp = 4
      igauss = 3
      ig = 6
      go to 700
  300 ityp = 4
      igauss = 1
      ig = 9
      go to 700
  320 ityp = 4
      igauss = 1
      ig = 10
      go to 700
  340 ityp = 18
      igauss = 1
      ig = 11
      go to 700
  480 if (nucz .ge. 17) go to 900
c
c     ----- na  - to - s -----
c           11(sp) or 12(spd) shells
c
      go to (
     *500,520,540,560,580,600,610,620,630,640,650,655),ipass
 500  ityp = 16
      igauss = 6
      ig = 0
      go to 700
 520  ityp = 16
      igauss = 3
      ig = 6
      go to 700
 540  ityp = 16
      igauss = 1
      ig = 9
      go to 700
 560  ityp = 16
      igauss = 1
      ig = 10
      go to 700
 580  ityp = 16
      igauss = 1
      ig = 11
      go to 700
 600  ityp = 16
      igauss = 1
      ig = 12
      go to 700
 610  ityp = 17
      igauss = 4
      ig = 13
      go to 700
 620  ityp = 17
      igauss = 2
      ig = 17
      go to 700
 630  ityp = 17
      igauss = 1
      ig = 19
      go to 700
 640  ityp = 17
      igauss = 1
      ig = 20
      go to 700
 650  ityp = 17
      igauss = 1
      ig = 21
      go to 700
 655  ityp = 18
      igauss = 1
      ig = 22
      go to 700
c
c     ----- cl and ar ----
c
 900  go to (
     *500,520,540,560,580,600,910,920,930,940,950,955),ipass
 910  ityp = 17
      igauss = 5
      ig = 13
      go to 700
 920  ityp = 17
      igauss = 2
      ig = 18
      go to 700
 930  ityp = 17
      igauss = 1
      ig = 20
      go to 700
 940  ityp = 17
      igauss = 1
      ig = 21
      go to 700
 950  ityp = 17
      igauss = 1
      ig = 22
      go to 700
 955  ityp = 18
      igauss = 1
      ig = 23
c
  700 continue
      nshell = nshell+1
      if (nshell .gt. nshmax) ierr1 = 1
      if (ierr1 .ne. 0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss .gt. ngsmax) ierr2 = 1
      if (ierr2 .ne. 0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do 720 i = 1,igauss
      k = k1+i-1
      ex(k) = eex(ig+i)
      csinp(k) = ccs(ig+i)
      cpinp(k) = ccp(ig+i)
      cdinp(k) = ccd(ig+i)
      cs(k) = csinp(k)
      cp(k) = cpinp(k)
  720 cd(k) = cdinp(k)
c
c     ----- always unnormalize primitives... -----
c
      do 740 k = k1,k2
      ee = ex(k)+ex(k)
      facs = pi32/(ee*  dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      cs(k) = cs(k)/  dsqrt(facs)
      cp(k) = cp(k)/  dsqrt(facp)
  740 cd(k) = cd(k)/  dsqrt(facd)
c
      if (ipass .lt. mxpass) go to 180
      return
      end
      subroutine f6311(e,s,p,d,n)
c
c     ----- 6-311g/McLean basis (TZV)
c     R.Krishnan, J.S.Binkley, R.Seeger and J.A.Pople
c     J. Chem. Phys. 72 (1980) 650
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
c
      go to (1 ,2 ,
     *   10, 20, 40, 60, 80, 85, 90, 95,
     *  100,120,140,160,180,200,220,240),n
c
c     ----- h ----
c
1     continue
      e( 1)= 0.338650d+02
      e( 2)= 0.509479d+01
      e( 3)= 0.115879d+01
      e( 4)= 0.325840d00
      e( 5)= 0.102741d00
      e( 6)= 0.75d00
      s( 1)= 0.254938d-01
      s( 2)= 0.190373d00
      s( 3)= 0.852161d00
      s( 4)=1.0d0
      s( 5)=1.0d0
      p( 6)=1.0d0
      return
c
c     ----- he ----
c
 2    continue
      e( 1)= 0.981243d+02
      e( 2)= 0.147689d+02
      e( 3)= 0.331883d+01
      e( 4)= 0.874047d00
      e( 5)= 0.244564d00
      e( 6)= 0.75d00
      s( 1)= 0.287452d-01
      s( 2)= 0.208061d00
      s( 3)= 0.837635d00
      s( 4)=+1.0d0
      s( 5)=+1.0d0
      p( 6)=1.0d0
      return
c
c     ----- li ----
c
 10   continue
      e( 1)= 0.900460d+03
      e( 2)= 0.134433d+03
      e( 3)= 0.304365d+02
      e( 4)= 0.862639d+01
      e( 5)= 0.248332d+01
      e( 6)= 0.303179d00
      e( 7)= 0.486890d+01
      e( 8)= 0.856924d00
      e( 9)= 0.243227d00
      e(10)= 0.635070d-01
      e(11)= 0.243683d-01
      e(12)= 0.200000d00
      s( 1)= 0.228704d-02
      s( 2)= 0.176350d-01
      s( 3)= 0.873434d-01
      s( 4)= 0.280977d00
      s( 5)= 0.658741d00
      s( 6)= 0.118712d00
      s( 7)= 0.933293d-01
      s( 8)= 0.943045d00
      s( 9)=-0.279827d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.327661d-01
      p( 8)= 0.159792d00
      p( 9)= 0.885667d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- be ----
c
 20   continue
      e( 1)= 0.168280d+04
      e( 2)= 0.251715d+03
      e( 3)= 0.574116d+02
      e( 4)= 0.165171d+02
      e( 5)= 0.485364d+01
      e( 6)= 0.626863d00
      e( 7)= 0.830938d+01
      e( 8)= 0.174075d+01
      e( 9)= 0.485816d00
      e(10)= 0.163613d00
      e(11)= 0.567285d-01
      e(12)= 0.250000d00
      s( 1)= 0.228574d-02
      s( 2)= 0.175938d-01
      s( 3)= 0.863315d-01
      s( 4)= 0.281835d00
      s( 5)= 0.640594d00
      s( 6)= 0.144467d00
      s( 7)= 0.108621d00
      s( 8)= 0.927301d00
      s( 9)=-0.297169d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.361344d-01
      p( 8)= 0.216958d00
      p( 9)= 0.841839d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- b  -----
c
 40   continue
      e( 1)= 0.285889d+04
      e( 2)= 0.428140d+03
      e( 3)= 0.975282d+02
      e( 4)= 0.279693d+02
      e( 5)= 0.821577d+01
      e( 6)= 0.111278d+01
      e( 7)= 0.132415d+02
      e( 8)= 0.300166d+01
      e( 9)= 0.912856d00
      e(10)= 0.315454d00
      e(11)= 0.988563d-01
      e(12)= 0.401000d00
      s( 1)= 0.215375d-02
      s( 2)= 0.165823d-01
      s( 3)= 0.821870d-01
      s( 4)= 0.276618d00
      s( 5)= 0.629316d00
      s( 6)= 0.173770d00
      s( 7)= 0.117443d00
      s( 8)= 0.918002d00
      s( 9)=-0.265105d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.418100d-01
      p( 8)= 0.236575d00
      p( 9)= 0.816214d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- c -----
c
 60   continue
      e( 1)= 0.456324d+04
      e( 2)= 0.682024d+03
      e( 3)= 0.154973d+03
      e( 4)= 0.444553d+02
      e( 5)= 0.130290d+02
      e( 6)= 0.182773d+01
      e( 7)= 0.209642d+02
      e( 8)= 0.480331d+01
      e( 9)= 0.145933d+01
      e(10)= 0.483456d00
      e(11)= 0.145585d00
      e(12)= 0.626000d00
      s( 1)= 0.196665d-02
      s( 2)= 0.152306d-01
      s( 3)= 0.761269d-01
      s( 4)= 0.260801d00
      s( 5)= 0.616462d00
      s( 6)= 0.221006d00
      s( 7)= 0.114660d00
      s( 8)= 0.919999d00
      s( 9)=-0.303068d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.402487d-01
      p( 8)= 0.237594d00
      p( 9)= 0.815854d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- n  -----
c
 80   continue
      e( 1)= 0.629348d+04
      e( 2)= 0.949044d+03
      e( 3)= 0.218776d+03
      e( 4)= 0.636916d+02
      e( 5)= 0.188282d+02
      e( 6)= 0.272023d+01
      e( 7)= 0.306331d+02
      e( 8)= 0.702614d+01
      e( 9)= 0.211205d+01
      e(10)= 0.684009d00
      e(11)= 0.200878d00
      e(12)= 0.913000d00
      s( 1)= 0.196979d-02
      s( 2)= 0.149613d-01
      s( 3)= 0.735006d-01
      s( 4)= 0.248937d00
      s( 5)= 0.602460d00
      s( 6)= 0.256202d00
      s( 7)= 0.111906d00
      s( 8)= 0.921666d00
      s( 9)=-0.256919d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.383119d-01
      p( 8)= 0.237403d00
      p( 9)= 0.817592d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- o  -----
c
 85   continue
      e( 1)= 0.858850d+04
      e( 2)= 0.129723d+04
      e( 3)= 0.299296d+03
      e( 4)= 0.873771d+02
      e( 5)= 0.256789d+02
      e( 6)= 0.374004d+01
      e( 7)= 0.421175d+02
      e( 8)= 0.962837d+01
      e( 9)= 0.285332d+01
      e(10)= 0.905661d00
      e(11)= 0.255611d00
      e(12)= 0.129200d01
      s( 1)= 0.189515d-02
      s( 2)= 0.143859d-01
      s( 3)= 0.707320d-01
      s( 4)= 0.240001d00
      s( 5)= 0.594797d00
      s( 6)= 0.280802d00
      s( 7)= 0.113889d00
      s( 8)= 0.920811d00
      s( 9)=-0.327447d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.365114d-01
      p( 8)= 0.237153d00
      p( 9)= 0.819702d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- f  ----
c
 90   continue
      e( 1)= 0.114271d+05
      e( 2)= 0.172235d+04
      e( 3)= 0.395746d+03
      e( 4)= 0.115139d+03
      e( 5)= 0.336026d+02
      e( 6)= 0.491901d+01
      e( 7)= 0.554441d+02
      e( 8)= 0.126323d+02
      e( 9)= 0.371756d+01
      e(10)= 0.116545d+01
      e(11)= 0.321892d00
      e(12)= 0.175000d01
      s( 1)= 0.180093d-02
      s( 2)= 0.137419d-01
      s( 3)= 0.681334d-01
      s( 4)= 0.233325d00
      s( 5)= 0.589086d00
      s( 6)= 0.299505d00
      s( 7)= 0.114536d00
      s( 8)= 0.920512d00
      s( 9)=-0.337804d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.354609d-01
      p( 8)= 0.237451d00
      p( 9)= 0.820458d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- ne ----
c
 95   continue
      e( 1)= 0.139957d+05
      e( 2)= 0.211710d+04
      e( 3)= 0.490425d+03
      e( 4)= 0.143833d+03
      e( 5)= 0.419265d+02
      e( 6)= 0.615684d+01
      e( 7)= 0.691211d+02
      e( 8)= 0.158350d+02
      e( 9)= 0.467326d+01
      e(10)= 0.145756d+01
      e(11)= 0.397057d00
      e(12)= 0.230400d01
      s( 1)= 0.183276d-02
      s( 2)= 0.138827d-01
      s( 3)= 0.680687d-01
      s( 4)= 0.231328d00
      s( 5)= 0.585890d00
      s( 6)= 0.305883d00
      s( 7)= 0.119149d00
      s( 8)= 0.917375d00
      s( 9)=-0.405839d-02
      s(10)=1.0d0
      s(11)=1.0d0
      p( 7)= 0.356574d-01
      p( 8)= 0.239477d00
      p( 9)= 0.818461d00
      p(10)=1.0d0
      p(11)=1.0d0
      d(12)=1.0d0
      return
c
c     ----- na ----
c
100   continue
      e( 1) = 0.3616640000d+05
      s( 1) = 0.1032000000d-02
      e( 2) = 0.5372580000d+04
      s( 2) = 0.8071000000d-02
      e( 3) = 0.1213210000d+04
      s( 3) = 0.4212900000d-01
      e( 4) = 0.3396230000d+03
      s( 4) = 0.1697890000d+00
      e( 5) = 0.1095530000d+03
      s( 5) = 0.5146210000d+00
      e( 6) = 0.3877730000d+02
      s( 6) = 0.3798170000d+00
      e( 7) = 0.3877730000d+02
      s( 7) = 0.3747620000d+00
      e( 8) = 0.1457590000d+02
      s( 8) = 0.5757690000d+00
      e( 9) = 0.5269930000d+01
      s( 9) = 0.1129330000d+00
      e(10) = 0.1827770000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.6199480000d+00
      s(11) = 0.1000000000d+01
      e(12) = 0.5724000000d-01
      s(12) = 0.1000000000d+01
      e(13) = 0.2404800000d-01
      s(13) = 0.1000000000d+01
      e(14) = 0.1446450000d+03
      p(14) = 0.1148500000d-01
      e(15) = 0.3390740000d+02
      p(15) = 0.8238300000d-01
      e(16) = 0.1062850000d+02
      p(16) = 0.3196580000d+00
      e(17) = 0.3823890000d+01
      p(17) = 0.7012950000d+00
      e(18) = 0.1444290000d+01
      p(18) = 0.6385060000d+00
      e(19) = 0.5526210000d+00
      p(19) = 0.4253650000d+00
      e(20) = 0.1887200000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.4650100000d-01
      p(21) = 0.1000000000d+01
      e(22) = 0.1628500000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.157d0
      d(23) = 1.0d0
      return
c
c     ----- mg ----
c
 120  continue
      e( 1) = 0.4386650000d+05
      s( 1) = 0.9180000000d-03
      e( 2) = 0.6605370000d+04
      s( 2) = 0.7047000000d-02
      e( 3) = 0.1513260000d+04
      s( 3) = 0.3594100000d-01
      e( 4) = 0.4323170000d+03
      s( 4) = 0.1414610000d+00
      e( 5) = 0.1421490000d+03
      s( 5) = 0.4267640000d+00
      e( 6) = 0.5139830000d+02
      s( 6) = 0.4979750000d+00
      e( 7) = 0.5139830000d+02
      s( 7) = 0.2513550000d+00
      e( 8) = 0.1991960000d+02
      s( 8) = 0.6186710000d+00
      e( 9) = 0.8024740000d+01
      s( 9) = 0.1884170000d+00
      e(10) = 0.2508170000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.8715310000d+00
      s(11) = 0.1000000000d+01
      e(12) = 0.1081880000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.4013000000d-01
      s(13) = 0.1000000000d+01
      e(14) = 0.1938540000d+03
      p(14) = 0.1018800000d-01
      e(15) = 0.4544200000d+02
      p(15) = 0.7536000000d-01
      e(16) = 0.1418640000d+02
      p(16) = 0.3074190000d+00
      e(17) = 0.5057510000d+01
      p(17) = 0.7175750000d+00
      e(18) = 0.1888610000d+01
      p(18) = 0.6673390000d+00
      e(19) = 0.7226520000d+00
      p(19) = 0.3946490000d+00
      e(20) = 0.2364170000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.9335800000d-01
      p(21) = 0.1000000000d+01
      e(22) = 0.3480900000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.234d0
      d(23) = 1.0d0
      return
c
c     ----- al -----
c
  140 continue
      e( 1) = 0.5486648900d+05
      s( 1) = 0.8390000000d-03
      e( 2) = 0.8211766500d+04
      s( 2) = 0.6527000000d-02
      e( 3) = 0.1866176100d+04
      s( 3) = 0.3366600000d-01
      e( 4) = 0.5311293400d+03
      s( 4) = 0.1329020000d+00
      e( 5) = 0.1751179700d+03
      s( 5) = 0.4012660000d+00
      e( 6) = 0.6400550000d+02
      s( 6) = 0.5313380000d+00
      e( 7) = 0.6400550000d+02
      s( 7) = 0.2023050000d+00
      e( 8) = 0.2529250700d+02
      s( 8) = 0.6247900000d+00
      e( 9) = 0.1053491000d+02
      s( 9) = 0.2274390000d+00
      e(10) = 0.3206711000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.1152555000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.1766780000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.6523700000d-01
      s(13) = 0.1000000000d+01
      e(14) = 0.2592836200d+03
      p(14) = 0.9448000000d-02
      e(15) = 0.6107687000d+02
      p(15) = 0.7097400000d-01
      e(16) = 0.1930323700d+02
      p(16) = 0.2956360000d+00
      e(17) = 0.7010882000d+01
      p(17) = 0.7282190000d+00
      e(18) = 0.2673865000d+01
      p(18) = 0.6444670000d+00
      e(19) = 0.1036596000d+01
      p(19) = 0.4174130000d+00
      e(20) = 0.3168190000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.1142570000d+00
      p(21) = 0.1000000000d+01
      e(22) = 0.4139700000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.311d0
      d(23) = 1.0d0
      return
c
c     ----- si -----
c
  160 continue
      e( 1) = 0.6937923000d+05
      s( 1) = 0.7570000000d-03
      e( 2) = 0.1035494000d+05
      s( 2) = 0.5932000000d-02
      e( 3) = 0.2333879600d+04
      s( 3) = 0.3108800000d-01
      e( 4) = 0.6571429500d+03
      s( 4) = 0.1249670000d+00
      e( 5) = 0.2143011300d+03
      s( 5) = 0.3868970000d+00
      e( 6) = 0.7762916800d+02
      s( 6) = 0.5548880000d+00
      e( 7) = 0.7762916800d+02
      s( 7) = 0.1778810000d+00
      e( 8) = 0.3063080700d+02
      s( 8) = 0.6277650000d+00
      e( 9) = 0.1280129500d+02
      s( 9) = 0.2476230000d+00
      e(10) = 0.3926866000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.1452343000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.2562340000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.9427900000d-01
      s(13) = 0.1000000000d+01
      e(14) = 0.3354831900d+03
      p(14) = 0.8866000000d-02
      e(15) = 0.7890036600d+02
      p(15) = 0.6829900000d-01
      e(16) = 0.2498815000d+02
      p(16) = 0.2909580000d+00
      e(17) = 0.9219711000d+01
      p(17) = 0.7321170000d+00
      e(18) = 0.3621140000d+01
      p(18) = 0.6198790000d+00
      e(19) = 0.1451310000d+01
      p(19) = 0.4391480000d+00
      e(20) = 0.5049770000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.1863170000d+00
      p(21) = 0.1000000000d+01
      e(22) = 0.6543200000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.388d0
      d(23) = 1.0d0
      return
c
c     ----- p  -----
c
  180 continue
      e( 1) = 0.7749240000d+05
      s( 1) = 0.7810000000d-03
      e( 2) = 0.1160580000d+05
      s( 2) = 0.6068000000d-02
      e( 3) = 0.2645960000d+04
      s( 3) = 0.3116000000d-01
      e( 4) = 0.7549760000d+03
      s( 4) = 0.1234310000d+00
      e( 5) = 0.2487550000d+03
      s( 5) = 0.3782090000d+00
      e( 6) = 0.9115650000d+02
      s( 6) = 0.5632620000d+00
      e( 7) = 0.9115650000d+02
      s( 7) = 0.1602550000d+00
      e( 8) = 0.3622570000d+02
      s( 8) = 0.6276470000d+00
      e( 9) = 0.1521130000d+02
      s( 9) = 0.2638490000d+00
      e(10) = 0.4794170000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.1807930000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.3568160000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.1147830000d+00
      s(13) = 0.1000000000d+01
      e(14) = 0.3848430000d+03
      p(14) = 0.9206000000d-02
      e(15) = 0.9055210000d+02
      p(15) = 0.6987400000d-01
      e(16) = 0.2913390000d+02
      p(16) = 0.2924700000d+00
      e(17) = 0.1088620000d+02
      p(17) = 0.7281030000d+00
      e(18) = 0.4352590000d+01
      p(18) = 0.6283490000d+00
      e(19) = 0.1777060000d+01
      p(19) = 0.4280440000d+00
      e(20) = 0.6970050000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.2535320000d+00
      p(21) = 0.1000000000d+01
      e(22) = 0.6849300000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.465d0
      d(23) = 1.0d0
      return
c
c     ----- s  -----
c
  200 continue
      e( 1) = 0.9341340000d+05
      s( 1) = 0.7430000000d-03
      e( 2) = 0.1396170000d+05
      s( 2) = 0.5793000000d-02
      e( 3) = 0.3169910000d+04
      s( 3) = 0.2995400000d-01
      e( 4) = 0.9024560000d+03
      s( 4) = 0.1190280000d+00
      e( 5) = 0.2971580000d+03
      s( 5) = 0.3684320000d+00
      e( 6) = 0.1087020000d+03
      s( 6) = 0.5772990000d+00
      e( 7) = 0.1087020000d+03
      s( 7) = 0.1431860000d+00
      e( 8) = 0.4315530000d+02
      s( 8) = 0.6244650000d+00
      e( 9) = 0.1810790000d+02
      s( 9) = 0.2833660000d+00
      e(10) = 0.5560090000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.2131830000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.4204030000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.1360450000d+00
      s(13) = 0.1000000000d+01
      e(14) = 0.4950400000d+03
      p(14) = 0.8309000000d-02
      e(15) = 0.1172210000d+03
      p(15) = 0.6402400000d-01
      e(16) = 0.3777490000d+02
      p(16) = 0.2776140000d+00
      e(17) = 0.1405840000d+02
      p(17) = 0.7450760000d+00
      e(18) = 0.5565740000d+01
      p(18) = 0.6137120000d+00
      e(19) = 0.2262970000d+01
      p(19) = 0.4438180000d+00
      e(20) = 0.8079940000d+00
      p(20) = 0.1000000000d+01
      e(21) = 0.2774600000d+00
      p(21) = 0.1000000000d+01
      e(22) = 0.7714100000d-01
      p(22) = 0.1000000000d+01
      e(23) = 0.542d0
      d(23) = 1.0d0
      return
c
c     ----- cl ----
c
 220  continue
      e( 1) = 0.1058190000d+06
      s( 1) = 0.7380000000d-03
      e( 2) = 0.1587200000d+05
      s( 2) = 0.5718000000d-02
      e( 3) = 0.3619650000d+04
      s( 3) = 0.2949500000d-01
      e( 4) = 0.1030800000d+04
      s( 4) = 0.1172860000d+00
      e( 5) = 0.3399080000d+03
      s( 5) = 0.3629490000d+00
      e( 6) = 0.1245380000d+03
      s( 6) = 0.5841490000d+00
      e( 7) = 0.1245380000d+03
      s( 7) = 0.1341770000d+00
      e( 8) = 0.4951350000d+02
      s( 8) = 0.6242500000d+00
      e( 9) = 0.2080560000d+02
      s( 9) = 0.2917560000d+00
      e(10) = 0.6583460000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.2564680000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.5597630000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.1832730000d+00
      s(13) = 0.1000000000d+01
      e(14) = 0.5897760000d+03
      p(14) = 0.2391000000d-02
      e(15) = 0.1398490000d+03
      p(15) = 0.1850400000d-01
      e(16) = 0.4514130000d+02
      p(16) = 0.8137700000d-01
      e(17) = 0.1687330000d+02
      p(17) = 0.2215520000d+00
      e(18) = 0.6741100000d+01
      p(18) = 0.7725690000d+00
      e(19) = 0.6741100000d+01
      p(19) =-0.1572244000d+01
      e(20) = 0.2771520000d+01
      p(20) = 0.9923890000d+00
      e(21) = 0.1023870000d+01
      p(21) = 0.1000000000d+01
      e(22) = 0.3813680000d+00
      p(22) = 0.1000000000d+01
      e(23) = 0.1094370000d+00
      p(23) = 0.1000000000d+01
      e(24) = 0.619d0
      d(24) = 1.0d0
      return
c
c     ----- ar ----
c
 240  continue
      e( 1) = 0.1180223800d+06
      s( 1) = 0.7470000000d-03
      e( 2) = 0.1768354100d+05
      s( 2) = 0.5790000000d-02
      e( 3) = 0.4027765700d+04
      s( 3) = 0.2991900000d-01
      e( 4) = 0.1145397700d+04
      s( 4) = 0.1192060000d+00
      e( 5) = 0.3771637500d+03
      s( 5) = 0.3690280000d+00
      e( 6) = 0.1381596900d+03
      s( 6) = 0.5764590000d+00
      e( 7) = 0.1381596900d+03
      s( 7) = 0.1439270000d+00
      e( 8) = 0.5498911700d+02
      s( 8) = 0.6229380000d+00
      e( 9) = 0.2317066700d+02
      s( 9) = 0.2839640000d+00
      e(10) = 0.7377860000d+01
      s(10) = 0.1000000000d+01
      e(11) = 0.2923688000d+01
      s(11) = 0.1000000000d+01
      e(12) = 0.6504050000d+00
      s(12) = 0.1000000000d+01
      e(13) = 0.2328250000d+00
      s(13) = 0.1000000000d+01
      e(14) = 0.6630620100d+03
      p(14) = 0.3082000000d-02
      e(15) = 0.1570928100d+03
      p(15) = 0.2416500000d-01
      e(16) = 0.5023110000d+02
      p(16) = 0.1082230000d+00
      e(17) = 0.1863534500d+02
      p(17) = 0.2941920000d+00
      e(18) = 0.7446537000d+01
      p(18) = 0.6878620000d+00
      e(19) = 0.7446537000d+01
      p(19) =-0.1214482000d+01
      e(20) = 0.3095698000d+01
      p(20) = 0.1632370000d+01
      e(21) = 0.1106463000d+01
      p(21) = 0.1000000000d+01
      e(22) = 0.4156010000d+00
      p(22) = 0.1000000000d+01
      e(23) = 0.1454490000d+00
      p(23) = 0.1000000000d+01
      e(24) = 0.696d0
      d(24) = 1.0d0
      return
      end
c
      subroutine mini1(csinp,cpinp,cdinp,cfinp,osplit,
     +                 scfacs,scfac,nucz,intyp,nangm,nbfs,minf,maxf,
     +                 loc,ngauss,ns,ngsmax,nat)
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
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax),cfinp(ngsmax),
     +          intyp(*),nangm(*),nbfs(*),minf(*),maxf(*),ns(*),scfac(*)
c
      common/blkin/eex(45),coef(45)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     * ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      data dzero,done,pt5,pt75,pt1875/
     + 0.0d+00,1.0d+00,0.5d+00,0.75d+00,1.875d+00/
      data pi32,tm6,tm10/5.56832799683170d+00,1.0d-06,1.0d-10/
c
c     ----- get the huzinaga mini or midi basis set -----
c
c     these basis sets are taken from "the green book",
c     "gaussian basis sets for molecular calculations"
c     edited by s. huzinaga, with coauthors j. andzelm,
c     m. klobukowski, e.radzio-andzelm, y.sakai, h.tatewaki.
c     elsevier, amsterdam, 1984.
c
c     the mini basis set is a minimal basis set, 3 gaussians
c     per atomic orbital (the same as sto-3g). however the exponents
c     are optimal for each atom, and s and p exponents are different.
c     midi basis sets are derived from the mini by floating the
c     outermost primitive of the valence orbitals, and then
c     renormalizing the remaining 2 gaussians.
c
      scale = dzero
      ng = -2**20
      igauss = ng
      ityp = ng
      ierr3=0
      odone = .false.
      call vclr(eex,1,45)
      call vclr(coef,1,45)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
        call hone(eex,coef,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
        call htwo(eex,coef,nucz)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
        call hthree(eex,coef,nucz)
c
c     ----- potassium to krypton -----
c
      else if(nucz.le.36) then
        call hfour(eex,coef,nucz)
c
c     ----- rubidium to xenon -----
c
      else if(nucz.le.54) then
        call hfive(eex,coef,nucz)
c
c     ----- cesium to radon -----
c
      else if(nucz.le.86) then
        call hsix(eex,coef,nucz)
c
      else
c
c     ----- past radon does not exist -----
c
      call caserr2('attempting to site function on invalid centre')
      endif
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
      scs = scfacs
      scp = scfac(1)
      scd = scfac(2)
      call hshell(nucz,ipass,osplit,ityp,igauss,ng,scale,
     *            scs,scp,scd,odone)
      if(odone. and. ierr3.ne.0) 
     + call caserr2('error detected in specifying huzinaga basis')
      if(odone) return
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      if(nshell.gt.mxshel) then
       call caserr2('excessive number of shells')
      end if
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if(ngauss.gt.mxprim) then
       call caserr2(
     * 'excessive number of contracted primitive functions')
      end if
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      scale = scale*scale
      do 440 i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i) * scale
         if(ityp.eq.16) csinp(k) = coef(ng+i)
         if(ityp.eq.17) cpinp(k) = coef(ng+i)
         if(ityp.eq.18) cdinp(k) = coef(ng+i)
         if(ityp.eq.19) cfinp(k) = coef(ng+i)
         cs(k) = dzero
         cp(k) = dzero
         cd(k) = dzero
         cf(k) = dzero
         if(ityp.eq.16) cs(k) = csinp(k)
         if(ityp.eq.17) cp(k) = cpinp(k)
         if(ityp.eq.18) cd(k) = cdinp(k)
         if(ityp.eq.19) cf(k) = cfinp(k)
  440 continue
c
c     ----- always unnormalize primitives -----
c
      do 460 k = k1,k2
         ee = ex(k)+ex(k)
         facs = pi32/(ee*sqrt(ee))
         facp = pt5*facs/ee
         facd = pt75*facs/(ee*ee)
         facf = pt1875*facs/(ee**3)
         if(ityp.eq.16) cs(k) = cs(k)/sqrt(facs)
         if(ityp.eq.17) cp(k) = cp(k)/sqrt(facp)
         if(ityp.eq.18) cd(k) = cd(k)/sqrt(facd)
         if(ityp.eq.19) cf(k) = cf(k)/sqrt(facf)
  460 continue
c
c     ----- if(normf.eq.0) normalize basis functions -----
c
      if (normf .eq. 1) go to 210
      facs = dzero
      facp = dzero
      facd = dzero
      facf = dzero
      do 510 ig = k1,k2
         do 500 jg = k1,ig
            ee = ex(ig)+ex(jg)
            fac = ee*sqrt(ee)
            dums = cs(ig)*cs(jg)/fac
            dump = pt5*cp(ig)*cp(jg)/(ee*fac)
            dumd = pt75*cd(ig)*cd(jg)/(ee*ee*fac)
            dumf = pt1875*cf(ig)*cf(jg)/(ee*ee*ee*fac)
            if (ig .eq. jg) go to 480
               dums = dums+dums
               dump = dump+dump
               dumd = dumd+dumd
               dumf = dumf+dumf
  480       continue
            facs = facs+dums
            facp = facp+dump
            facd = facd+dumd
            facf = facf+dumf
  500    continue
  510 continue
c
      fac=dzero
      if(ityp.eq.16.and. facs.gt.tm10) fac=done/sqrt(facs*pi32)
      if(ityp.eq.17.and. facp.gt.tm10) fac=done/sqrt(facp*pi32)
      if(ityp.eq.18.and. facd.gt.tm10) fac=done/sqrt(facd*pi32)
      if(ityp.eq.19.and. facf.gt.tm10) fac=done/sqrt(facf*pi32)
c                                 verify normalization
      tnorm = abs(done-fac)
      if(tnorm.lt.tm6) go to 520
      if(osplit  .and.  igauss.lt.3) go to 520
      if(osplit  .and.  ityp.eq.19 .and.  igauss.eq.3) go to 520
c
       write(iwr,9000) fac,nucz,ipass,ityp
       write(iwr,9010) (eex(ng+k),k=1,igauss)
       write(iwr,9010) (coef(ng+k),k=1,igauss)
      ierr3=ierr3+1
c
  520 continue
      do 550 ig = k1,k2
         if(ityp.eq.16) cs(ig) = fac*cs(ig)
         if(ityp.eq.17) cp(ig) = fac*cp(ig)
         if(ityp.eq.18) cd(ig) = fac*cd(ig)
         if(ityp.eq.19) cf(ig) = fac*cf(ig)
  550 continue
      go to 210
c
 9000 format(/1x,'error!!! normalization factor=',e16.8/
     *     1x,'for atom z=',i3,' shell no.',i3,' ityp=',i4/
     *     1x,'(ityp of 16,17,18,19 means an s,p,d,f shell)'/
     *     1x,'check built in huzinaga exponents and cont. coefs')
 9010 format(3f16.10)
      end
      subroutine hone(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
c     first row basis sets taken from f. b. van duijneveldt,
c     ibm technical report rj-945.
c
      if(nucz.eq.2) go to 200
c
c           hydrogen
c
      e(1)=4.501800d+00
      e(2)=0.681444d+00
      e(3)=0.151398d+00
      c(1)=0.070452d+00
      c(2)=0.407826d+00
      c(3)=0.647752d+00
      return
c
c           helium
c
  200 continue
      e(1)=13.626736d+00
      e(2)=1.999349d+00
      e(3)=0.382993d+00
      c(1)=0.080241d+00
      c(2)=0.409143d+00
      c(3)=0.657278d+00
      return
      end
      subroutine htwo(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      ibr=nucz-2
      go to (100,200,300,400,500,600,700,800), ibr
c
c           lithium
c
  100 continue
      e(1)=35.04615d+00
      e(2)=5.20169d+00
      e(3)=1.05624d+00
      c(1)=0.0737603d+00
      c(2)=0.3974715d+00
      c(3)=0.6650924d+00
      e(4)=0.71517d+00
      e(5)=0.07053d+00
      e(6)=0.02735d+00
      c(4)=-0.09397d+00
      c(5)=0.5701001d+00
      c(6)=0.4997501d+00
      return
c
c           beryllium
c
  200 continue
      e(1)=66.953540d+00
      e(2)=9.9392900d+00
      e(3)=2.05713d+00
      c(1)=0.0702000d+00
      c(2)=0.3919101d+00
      c(3)=0.6699701d+00
      e(4)=1.9122100d+00
      e(5)=0.16132d+00
      e(6)=0.05524d+00
      c(4)=-0.0828204d+00
      c(5)=0.5575528d+00
      c(6)=0.5160426d+00
      return
c
c           boron
c
  300 continue
      e(1)=108.43704d+00
      e(2)=16.120560d+00
      e(3)=3.3734300d+00
      c(1)=0.0686505d+00
      c(2)=0.3899329d+00
      c(3)=0.6713951d+00
      e(4)=3.5665d+00
      e(5)=0.29547d+00
      e(6)=0.09805d+00
      c(4)=-0.0824192d+00
      c(5)=0.5590644d+00
      c(6)=0.5167948d+00
      e(7)=2.57207d+00
      e(8)=0.51694d+00
      e(9)=0.12314d+00
      c(7)=0.1059001d+00
      c(8)=0.4571805d+00
      c(9)=0.6318607d+00
      return
c
c           carbon
c
  400 continue
      e(1)=153.17226d+00
      e(2)=23.073030d+00
      e(3)=4.9232900d+00
      c(1)=0.0707401d+00
      c(2)=0.3953804d+00
      c(3)=0.6633107d+00
      e(4)=5.7255700d+00
      e(5)=0.45504000d+00
      e(6)=0.14707d+00
      c(4)=-0.0813805d+00
      c(5)=0.5748532d+00
      c(6)=0.5024128d+00
      e(7)=4.25131d+00
      e(8)=0.86327d+00
      e(9)=0.20135d+00
      c(7)=0.1099306d+00
      c(8)=0.4627127d+00
      c(9)=0.6275137d+00
      return
c
c           nitrogen
c
  500 continue
      e(1)=218.36449d+00
      e(2)=32.59889d+00
      e(3)=6.91739d+00
      c(1)=0.0678703d+00
      c(2)=0.3902018d+00
      c(3)=0.6700832d+00
      e(4)=8.32638d+00
      e(5)=0.65919d+00
      e(6)=0.21009d+00
      c(4)=-0.0808903d+00
      c(5)=0.5672020d+00
      c(6)=0.5110918d+00
      e(7)=6.12035d+00
      e(8)=1.25938d+00
      e(9)=0.29145d+00
      c(7)=0.1159195d+00
      c(8)=0.4699582d+00
      c(9)=0.6184476d+00
      return
c
c           oxygen
c
  600 continue
      e(1)=281.86658d+00
      e(2)=42.416d+00
      e(3)=9.09562d+00
      c(1)=0.0690599d+00
      c(2)=0.3931595d+00
      c(3)=0.6656691d+00
      e(4)=11.46603d+00
      e(5)=0.88786d+00
      e(6)=0.2788d+00
      c(4)=-0.0808199d+00
      c(5)=0.5820895d+00
      c(6)=0.4971596d+00
      e(7)=8.04724d+00
      e(8)=1.66842d+00
      e(9)=0.37251d+00
      c(7)=0.1242709d+00
      c(8)=0.4765935d+00
      c(9)=0.6130445d+00
      return
c
c           fluorine
c
  700 continue
      e(1)=368.37112d+00
      e(2)=55.06106d+00
      e(3)=11.74767d+00
      c(1)=0.0670398d+00
      c(2)=0.3892491d+00
      c(3)=0.6707884d+00
      e(4)=15.15184d+00
      e(5)=1.15137d+00
      e(6)=0.35811d+00
      c(4)=-0.0805499d+00
      c(5)=0.587729d+00
      c(6)=0.4919792d+00
      e(7)=10.57707d+00
      e(8)=2.19498d+00
      e(9)=0.47937d+00
      c(7)=0.1262695d+00
      c(8)=0.4779482d+00
      c(9)=0.6140077d+00
      return
c
c           neon
c
  800 continue
      e(1)=456.95285d+00
      e(2)=68.365430d+00
      e(3)=14.619760d+00
      c(1)=0.0669098d+00
      c(2)=0.3893491d+00
      c(3)=0.6705184d+00
      e(4)=19.327190d+00
      e(5)=1.4418200d+00
      e(6)=0.44408000d+00
      c(4)=-0.0802497d+00
      c(5)=0.5952977d+00
      c(6)=0.4848681d+00
      e(7)=13.352520d+00
      e(8)=2.77947d+00
      e(9)=0.60097d+00
      c(7)=0.1288403d+00
      c(8)=0.4804412d+00
      c(9)=0.6116716d+00
      return
      end
      subroutine hthree(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      ibr=nucz-10
      go to (100,200,300,400,500,600,700,800), ibr
c
c           sodium
c
  100 continue
      e(1)=542.76053d+00
      e(2)=81.959470d+00
      e(3)=17.723770d+00
      c(1)=0.0684104d+00
      c(2)=0.3920923d+00
      c(3)=0.6660838d+00
      e(4)=23.280420d+00
      e(5)=1.8683400d+00
      e(6)=0.6232500d+00
      c(4)=-0.0838006d+00
      c(5)=0.5827943d+00
      c(6)=0.4924737d+00
      e(7)=0.5067300d+00
      e(8)=0.0535100d+00
      e(9)=0.0208000d+00
      c(7)=-0.1157622d+00
      c(8)=0.6958634d+00
      c(9)=0.3810473d+00
      e(10)=17.836360d+00
      e(11)=3.7956900d+00
      e(12)=0.8775100d+00
      c(10)=0.1257103d+00
      c(11)=0.4804611d+00
      c(12)=0.6022814d+00
      return
c
c           magnesium
c
  200 continue
      e(1)=650.643367d+00
      e(2)=98.370780d+00
      e(3)=21.322490d+00
      c(1)=0.0680297d+00
      c(2)=0.3907384d+00
      c(3)=0.6672673d+00
      e(4)=27.977380d+00
      e(5)=2.3265200d+00
      e(6)=0.8180800d+00
      c(4)=-0.0867195d+00
      c(5)=0.5856969d+00
      c(6)=0.4864974d+00
      e(7)=0.9113400d+00
      e(8)=0.0994000d+00
      e(9)=0.0362300d+00
      c(7)=-0.1276505d+00
      c(8)=0.6507728d+00
      c(9)=0.4362719d+00
      e(10)=23.216620d+00
      e(11)=5.0022200d+00
      e(12)=1.2046500d+00
      c(10)=0.1214603d+00
      c(11)=0.4792914d+00
      c(12)=0.5989417d+00
      return
c
c           aluminum
c
  300 continue
      e(1)=777.44334d+00
      e(2)=117.23153d+00
      e(3)=25.376297d+00
      c(1)=0.0668865d+00
      c(2)=0.3877681d+00
      c(3)=0.6707027d+00
      e(4)=33.356253d+00
      e(5)=2.8013153d+00
      e(6)=1.0227334d+00
      c(4)=-0.0889564d+00
      c(5)=0.6010610d+00
      c(6)=0.4687856d+00
      e(7)=1.3291880d+00
      e(8)=0.1616630d+00
      e(9)=0.05898876d+00
      c(7)=-0.1513892d+00
      c(8)=0.6593862d+00
      c(9)=0.4386928d+00
      e(10)=30.569580d+00
      e(11)=6.6447007d+00
      e(12)=1.6543954d+00
      c(10)=0.1123529d+00
      c(11)=0.4674668d+00
      c(12)=0.6097823d+00
      e(13)=0.0446518d+00
      e(14)=0.29965103d+00
      e(15)=0.11584816d+00
      c(13)=0.3964698d+00
      c(14)=0.2264326d+00
      c(15)=0.5005859d+00
      return
c
c           silicon
c
  400 continue
      e(1)=909.23487d+00
      e(2)=137.12456d+00
      e(3)=29.714810d+00
      c(1)=0.0664053d+00
      c(2)=0.3862219d+00
      c(3)=0.6722400d+00
      e(4)=39.129423d+00
      e(5)=3.3359807d+00
      e(6)=1.2512594d+00
      c(4)=-0.0909994d+00
      c(5)=0.6116149d+00
      c(6)=0.4568602d+00
      e(7)=1.8031012d+00
      e(8)=0.22638923d+00
      e(9)=0.08239593d+00
      c(7)=-0.1687328d+00
      c(8)=0.6754527d+00
      c(9)=0.4294192d+00
      e(10)=37.881761d+00
      e(11)=8.3045980d+00
      e(12)=2.1207924d+00
      c(10)=0.1087525d+00
      c(11)=0.4635152d+00
      c(12)=0.6113343d+00
      e(13)=0.44780246d+00
      e(14)=0.06236175d+00
      e(15)=0.17083752d+00
      c(13)=0.2389130d+00
      c(14)=0.3454528d+00
      c(15)=0.5422949d+00
      return
c
c           phosphorus
c
  500 continue
      e(1)=1053.2658d+00
      e(2)=158.79044d+00
      e(3)=34.424407d+00
      c(1)=0.0658649d+00
      c(2)=0.3845776d+00
      c(3)=0.6739634d+00
      e(4)=45.450377d+00
      e(5)=3.8999265d+00
      e(6)=1.4885071d+00
      c(4)=-0.0926548d+00
      c(5)=0.6265125d+00
      c(6)=0.4410392d+00
      e(7)=2.3008408d+00
      e(8)=0.29895928d+00
      e(9)=0.10885358d+00
      c(7)=-0.1805490d+00
      c(8)=0.6809516d+00
      c(9)=0.4291418d+00
      e(10)=46.100019d+00
      e(11)=10.165057d+00
      e(12)=2.6447940d+00
      c(10)=0.1053880d+00
      c(11)=0.4597120d+00
      c(12)=0.6137136d+00
      e(13)=0.63268608d+00
      e(14)=0.24021895d+00
      e(15)=0.08644653d+00
      c(13)=0.2358855d+00
      c(14)=0.5541603d+00
      c(15)=0.3365297d+00
      return
c
c           sulfur
c
  600 continue
      e(1)=1201.4584d+00
      e(2)=181.39212d+00
      e(3)=39.404795d+00
      c(1)=0.0657647d+00
      c(2)=0.3839476d+00
      c(3)=0.6743724d+00
      e(4)=52.139030d+00
      e(5)=4.5287990d+00
      e(6)=1.7549380d+00
      c(4)=-0.0942319d+00
      c(5)=0.6354680d+00
      c(6)=0.4315059d+00
      e(7)=2.8404372d+00
      e(8)=0.38143256d+00
      e(9)=0.13878583d+00
      c(7)=0.1900423d+00
      c(8)=-0.6855267d+00
      c(9)=-0.4292724d+00
      e(10)=54.644071d+00
      e(11)=12.122902d+00
      e(12)=3.2065037d+00
      c(10)=0.1036728d+00
      c(11)=0.4581901d+00
      c(12)=0.6134001d+00
      e(13)=0.86327394d+00
      e(14)=0.10867870d+00
      e(15)=0.31812973d+00
      c(13)=0.2294364d+00
      c(14)=0.3537001d+00
      c(15)=0.5529600d+00
      return
c
c           chlorine
c
  700 continue
      e(1)=1362.0220d+00
      e(2)=205.81110d+00
      e(3)=44.772167d+00
      c(1)=0.0655436d+00
      c(2)=0.3829867d+00
      c(3)=0.6752103d+00
      e(4)=59.225732d+00
      e(5)=5.2139024d+00
      e(6)=2.0473460d+00
      c(4)=-0.0956202d+00
      c(5)=0.6414256d+00
      c(6)=0.4251534d+00
      e(7)=3.4197116d+00
      e(8)=0.47001739d+00
      e(9)=0.16995836d+00
      c(7)=0.1964008d+00
      c(8)=-0.6923599d+00
      c(9)=-0.4261933d+00
      e(10)=64.099958d+00
      e(11)=14.287139d+00
      e(12)=3.8281352d+00
      c(10)=0.1017893d+00
      c(11)=0.4561069d+00
      c(12)=0.6142825d+00
      e(13)=1.0951257d+00
      e(14)=0.13217647d+00
      e(15)=0.39600382d+00
      c(13)=0.2359029d+00
      c(14)=0.3465996d+00
      c(15)=0.5580658d+00
      return
c
c           argon
c
  800 continue
      e(1)=1536.9325d+00
      e(2)=232.17765d+00
      e(3)=50.521685d+00
      c(1)=0.0651589d+00
      c(2)=0.3818067d+00
      c(3)=0.6764462d+00
      e(4)=66.933949d+00
      e(5)=5.9185516d+00
      e(6)=2.3393420d+00
      c(4)=-0.0967405d+00
      c(5)=0.6527492d+00
      c(6)=0.4135728d+00
      e(7)=4.0453072d+00
      e(8)=0.56570104d+00
      e(9)=0.20406479d+00
      c(7)=0.2007356d+00
      c(8)=-0.6962680d+00
      c(9)=-0.4248435d+00
      e(10)=74.352915d+00
      e(11)=16.631346d+00
      e(12)=4.5039275d+00
      c(10)=0.1000791d+00
      c(11)=0.4542258d+00
      c(12)=0.6152592d+00
      e(13)=1.3570912d+00
      e(14)=0.48811284d+00
      e(15)=0.16212589d+00
      c(13)=0.2372757d+00
      c(14)=0.5583599d+00
      c(15)=0.3461648d+00
      return
      end
      subroutine hfour(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      if(nucz.eq.19) then
c
c           potassium
c
      e(1)=1721.1755d+00
      e(2)=260.01633d+00
      e(3)=56.624554d+00
      c(1)=0.0648747d+00
      c(2)=0.3808593d+00
      c(3)=0.6773681d+00
      e(4)=75.055600d+00
      e(5)=6.6911626d+00
      e(6)=2.6671665d+00
      c(4)=-0.0978731d+00
      c(5)=0.6595598d+00
      c(6)=0.4065295d+00
      e(7)=4.6689367d+00
      e(8)=0.70001332d+00
      e(9)=0.27533364d+00
      c(7)=-0.2129302d+00
      c(8)=0.6892491d+00
      c(9)=0.4339254d+00
      e(10)=0.25231153d+00
      e(11)=0.03763114d+00
      e(12)=0.01621835d+00
      c(10)=-0.1529800d+00
      c(11)=0.6817284d+00
      c(12)=0.4144829d+00
      e(13)=85.789846d+00
      e(14)=19.254794d+00
      e(15)=5.2686239d+00
      c(13)=0.0980354d+00
      c(14)=0.4510408d+00
      c(15)=0.6177076d+00
      e(16)=1.6831452d+00
      e(17)=0.62580945d+00
      e(18)=0.22389802d+00
      c(16)=0.2386360d+00
      c(17)=0.5710568d+00
      c(18)=0.3165357d+00
c
      else if(nucz.eq.20) then
c
c           calcium
c
      e(1)=1915.4348d+00
      e(2)=289.53324d+00
      e(3)=63.106352d+00
      c(1)=0.0646237d+00
      c(2)=0.3798376d+00
      c(3)=0.6783294d+00
      e(4)=83.633281d+00
      e(5)=7.5118396d+00
      e(6)=3.0145846d+00
      c(4)=-0.0988822d+00
      c(5)=0.6660273d+00
      c(6)=0.3999119d+00
      e(7)=5.3707540d+00
      e(8)=0.83837962d+00
      e(9)=0.34622613d+00
      c(7)=-0.2235309d+00
      c(8)=0.7028863d+00
      c(9)=0.4235525d+00
      e(10)=0.42974055d+00
      e(11)=0.06217342d+00
      e(12)=0.02494803d+00
      c(10)=-0.1793098d+00
      c(11)=0.6728792d+00
      c(12)=0.4356664d+00
      e(13)=97.974592d+00
      e(14)=22.067384d+00
      e(15)=6.0938756d+00
      c(13)=0.0963163d+00
      c(14)=0.4481076d+00
      c(15)=0.6199212d+00
      e(16)=2.0178863d+00
      e(17)=0.76665059d+00
      e(18)=0.28431923d+00
      c(16)=0.2450293d+00
      c(17)=0.5843849d+00
      c(18)=0.2876662d+00
c
      else if(nucz.le.30) then
        call htmone(e,c,nucz)
      else
      ibr=nucz-30
      go to (300,400,500,600,700,800), ibr
c
c           gallium
c
  300 continue
      e(1)=4751.8979d+00
      e(2)=718.92054d+00
      e(3)=157.44592d+00
      c(1)=0.0628396d+00
      c(2)=0.3736112d+00
      c(3)=0.6843626d+00
      e(4)=210.53833d+00
      e(5)=19.794805d+00
      e(6)=8.0497965d+00
      c(4)=-0.1054027d+00
      c(5)=0.7070778d+00
      c(6)=0.3599497d+00
      e(7)=16.424933d+00
      e(8)=2.6576527d+00
      e(9)=1.0609873d+00
      c(7)=0.2270732d+00
      c(8)=-0.7111408d+00
      c(9)=-0.4213859d+00
      e(10)=1.3533154d+00
      e(11)=0.17615808d+00
      e(12)=0.06371947d+00
      c(10)=-0.1507084d+00
      c(11)=0.6413441d+00
      c(12)=0.4600643d+00
      e(13)=273.64584d+00
      e(14)=62.815151d+00
      e(15)=18.134221d+00
      c(13)=0.0908956d+00
      c(14)=0.4456493d+00
      c(15)=0.6179026d+00
      e(16)=6.8387139d+00
      e(17)=2.5958668d+00
      e(18)=0.98389370d+00
      c(16)=0.2739015d+00
      c(17)=0.5689512d+00
      c(18)=0.2756075d+00
      e(19)=0.31568846d+00
      e(20)=0.10961938d+00
      e(21)=0.03994387d+00
      c(19)=0.2062964d+00
      c(20)=0.5563565d+00
      c(21)=0.3670988d+00
      e(22)=21.292530d+00
      e(23)=5.3931662d+00
      e(24)=1.3338828d+00
      c(22)=0.1619895d+00
      c(23)=0.5116739d+00
      c(24)=0.5898732d+00
      return
c
c           germanium
c
  400 continue
      e(1)=5073.7499d+00
      e(2)=767.72417d+00
      e(3)=168.18881d+00
      c(1)=0.0627249d+00
      c(2)=0.3731671d+00
      c(3)=0.6847867d+00
      e(4)=224.98366d+00
      e(5)=21.187659d+00
      e(6)=8.6256370d+00
      c(4)=-0.1057096d+00
      c(5)=0.7094592d+00
      c(6)=0.3575586d+00
      e(7)=17.764305d+00
      e(8)=2.8551503d+00
      e(9)=1.1585827d+00
      c(7)=-0.2280577d+00
      c(8)=0.7250544d+00
      c(9)=0.4064506d+00
      e(10)=1.5859493d+00
      e(11)=0.22692409d+00
      e(12)=0.08372494d+00
      c(10)=-0.1769108d+00
      c(11)=0.6604420d+00
      c(12)=0.4525469d+00
      e(13)=293.73174d+00
      e(14)=67.552546d+00
      e(15)=19.568631d+00
      c(13)=0.0906975d+00
      c(14)=0.4451940d+00
      c(15)=0.6178940d+00
      e(16)=7.4326168d+00
      e(17)=2.8426842d+00
      e(18)=1.1041168d+00
      c(16)=0.2751638d+00
      c(17)=0.5703353d+00
      c(18)=0.2690321d+00
      e(19)=0.38224145d+00
      e(20)=0.14147296d+00
      e(21)=0.05262765d+00
      c(19)=0.2627754d+00
      c(20)=0.5701236d+00
      c(21)=0.2901737d+00
      e(22)=24.321421d+00
      e(23)=6.2238135d+00
      e(24)=1.5887375d+00
      c(22)=0.1577985d+00
      c(23)=0.5114922d+00
      c(24)=0.5857703d+00
      return
c
c           arsenic
c
  500 continue
      e(1)=5407.6138d+00
      e(2)=818.17436d+00
      e(3)=179.26569d+00
      c(1)=0.0626011d+00
      c(2)=0.3727790d+00
      c(3)=0.6851842d+00
      e(4)=239.85464d+00
      e(5)=22.64906d+00
      e(6)=9.2280039d+00
      c(4)=-0.1060438d+00
      c(5)=0.7107219d+00
      c(6)=0.3563753d+00
      e(7)=19.058432d+00
      e(8)=3.0854133d+00
      e(9)=1.2727348d+00
      c(7)=-0.2308017d+00
      c(8)=0.7325742d+00
      c(9)=0.3994769d+00
      e(10)=1.8198631d+00
      e(11)=0.27883046d+00
      e(12)=0.10372930d+00
      c(10)=-0.1983584d+00
      c(11)=0.6759512d+00
      c(12)=0.4476004d+00
      e(13)=315.02252d+00
      e(14)=72.557901d+00
      e(15)=21.079855d+00
      c(13)=0.0903237d+00
      c(14)=0.4443841d+00
      c(15)=0.6183953d+00
      e(16)=8.1044166d+00
      e(17)=3.1282025d+00
      e(18)=1.2391193d+00
      c(16)=0.2728475d+00
      c(17)=0.5717222d+00
      c(18)=0.2664488d+00
      e(19)=0.45260113d+00
      e(20)=0.17204402d+00
      e(21)=0.06573731d+00
      c(19)=0.3086721d+00
      c(20)=0.5671752d+00
      c(21)=0.2410732d+00
      e(22)=27.437209d+00
      e(23)=7.0840440d+00
      e(24)=1.8558226d+00
      c(22)=0.1544952d+00
      c(23)=0.5114318d+00
      c(24)=0.5821935d+00
      return
c
c           selenium
c
  600 continue
      e(1)=5751.3215d+00
      e(2)=870.25721d+00
      e(3)=190.72949d+00
      c(1)=0.0624934d+00
      c(2)=0.3723683d+00
      c(3)=0.6855799d+00
      e(4)=255.18802d+00
      e(5)=24.155395d+00
      e(6)=9.8525390d+00
      c(4)=-0.1063403d+00
      c(5)=0.7120577d+00
      c(6)=0.3550792d+00
      e(7)=20.463481d+00
      e(8)=3.3125921d+00
      e(9)=1.3791545d+00
      c(7)=-0.2326525d+00
      c(8)=0.7465053d+00
      c(9)=0.3857492d+00
      e(10)=2.0492750d+00
      e(11)=0.33641155d+00
      e(12)=0.12612688d+00
      c(10)=-0.2169517d+00
      c(11)=0.6852937d+00
      c(12)=0.4486601d+00
      e(13)=337.08753d+00
      e(14)=77.748281d+00
      e(15)=22.650400d+00
      c(13)=0.0899627d+00
      c(14)=0.4435785d+00
      c(15)=0.6189101d+00
      e(16)=8.6914000d+00
      e(17)=3.3668794d+00
      e(18)=1.3546532d+00
      c(16)=0.2789803d+00
      c(17)=0.5760239d+00
      c(18)=0.2529785d+00
      e(19)=0.54258450d+00
      e(20)=0.20662960d+00
      e(21)=0.07760208d+00
      c(19)=0.3286253d+00
      c(20)=0.5489677d+00
      c(21)=0.2444426d+00
      e(22)=30.627464d+00
      e(23)=7.9712764d+00
      e(24)=2.1348097d+00
      c(22)=0.1519858d+00
      c(23)=0.5116403d+00
      c(24)=0.5786936d+00
      return
c
c           bromine
c
  700 continue
      e(1)=6103.2899d+00
      e(2)=923.69743d+00
      e(3)=202.52031d+00
      c(1)=0.0624175d+00
      c(2)=0.3720414d+00
      c(3)=0.6858728d+00
      e(4)=270.99649d+00
      e(5)=25.706395d+00
      e(6)=10.501175d+00
      c(4)=-0.1066146d+00
      c(5)=0.7133792d+00
      c(6)=0.3537553d+00
      e(7)=21.859513d+00
      e(8)=3.5681587d+00
      e(9)=1.5090061d+00
      c(7)=-0.2354458d+00
      c(8)=0.7528340d+00
      c(9)=0.3801608d+00
      e(10)=2.3404284d+00
      e(11)=0.38825413d+00
      e(12)=0.14598986d+00
      c(10)=-0.2363212d+00
      c(11)=0.7085531d+00
      c(12)=0.4329347d+00
      e(13)=360.49973d+00
      e(14)=83.244591d+00
      e(15)=24.311838d+00
      c(13)=0.0894281d+00
      c(14)=0.4423065d+00
      c(15)=0.6200053d+00
      e(16)=9.4221830d+00
      e(17)=3.6771841d+00
      e(18)=1.4935778d+00
      c(16)=0.2767522d+00
      c(17)=0.5795347d+00
      c(18)=0.2491162d+00
      e(19)=0.63454669d+00
      e(20)=0.24252803d+00
      e(21)=0.09053005d+00
      c(19)=0.3466285d+00
      c(20)=0.5413298d+00
      c(21)=0.2347667d+00
      e(22)=33.965097d+00
      e(23)=8.9008312d+00
      e(24)=2.4284360d+00
      c(22)=0.1496666d+00
      c(23)=0.5117475d+00
      c(24)=0.5759148d+00
      return
c
c           krypton
c
  800 continue
      e(1)=6446.6307d+00
      e(2)=976.87570d+00
      e(3)=214.47955d+00
      c(1)=0.0625398d+00
      c(2)=0.3721075d+00
      c(3)=0.6856107d+00
      e(4)=287.01158d+00
      e(5)=27.328672d+00
      e(6)=11.174440d+00
      c(4)=-0.1068933d+00
      c(5)=0.7139493d+00
      c(6)=0.3532914d+00
      e(7)=23.331340d+00
      e(8)=3.8133080d+00
      e(9)=1.6134940d+00
      c(7)=-0.2376841d+00
      c(8)=0.7698488d+00
      c(9)=0.3641119d+00
      e(10)=2.5964926d+00
      e(11)=0.4490900d+00
      e(12)=0.17012362d+00
      c(10)=-0.2498594d+00
      c(11)=0.7128783d+00
      c(12)=0.4361011d+00
      e(13)=384.37950d+00
      e(14)=88.868156d+00
      e(15)=26.019696d+00
      c(13)=0.0890288d+00
      c(14)=0.4413316d+00
      c(15)=0.6207416d+00
      e(16)=10.184792d+00
      e(17)=4.0094356d+00
      e(18)=1.6459659d+00
      c(16)=0.2745449d+00
      c(17)=0.5813334d+00
      c(18)=0.2470820d+00
      e(19)=0.72962570d+00
      e(20)=0.28060560d+00
      e(21)=0.10521386d+00
      c(19)=0.3620055d+00
      c(20)=0.5352020d+00
      c(21)=0.2242007d+00
      e(22)=37.368103d+00
      e(23)=9.8543131d+00
      e(24)=2.7327955d+00
      c(22)=0.1479466d+00
      c(23)=0.5121719d+00
      c(24)=0.5729498d+00
      endif
      return
      end
      subroutine hfive(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      if(nucz.eq.37) then 
c
c           rubidium
c
      e(1)=6816.7225d+00
      e(2)=1033.0007d+00
      e(3)=226.90861d+00
      c(1)=0.0624962d+00
      c(2)=0.3719500d+00
      c(3)=0.6857293d+00
      e(4)=302.83628d+00
      e(5)=29.201467d+00
      e(6)=12.082845d+00
      c(4)=-0.1074004d+00
      c(5)=0.7049410d+00
      c(6)=0.3622692d+00
      e(7)=24.692611d+00
      e(8)=4.2016508d+00
      e(9)=1.8522603d+00
      c(7)=-0.2430518d+00
      c(8)=0.7375963d+00
      c(9)=0.3988455d+00
      e(10)=2.9359876d+00
      e(11)=0.52542640d+00
      e(12)=0.21916888d+00
      c(10)=-0.2806358d+00
      c(11)=0.7159840d+00
      c(12)=0.4402934d+00
      e(13)=0.19166785d+00
      e(14)=0.03343655d+00
      e(15)=0.01473352d+00
      c(13)=0.1580968d+00
      c(14)=-0.6855144d+00
      c(15)=-0.4191493d+00
      e(16)=407.86841d+00
      e(17)=94.455492d+00
      e(18)=27.745510d+00
      c(16)=0.0890053d+00
      c(17)=0.4411370d+00
      c(18)=0.6204330d+00
      e(19)=10.715052d+00
      e(20)=4.2311626d+00
      e(21)=1.7434701d+00
      c(19)=0.2867814d+00
      c(20)=0.5863660d+00
      c(21)=0.2277478d+00
      e(22)=0.84991897d+00
      e(23)=0.34353368d+00
      e(24)=0.13941053d+00
      c(22)=0.3796538d+00
      c(23)=0.5370593d+00
      c(24)=0.1869077d+00
      e(25)=40.866031d+00
      e(26)=10.840885d+00
      e(27)=3.0508341d+00
      c(25)=0.1466037d+00
      c(26)=0.5127252d+00
      c(27)=0.5699804d+00
c
      else if(nucz.eq.38) then
c
c           strontium
c
      e(1)=7215.4735d+00
      e(2)=1092.8519d+00
      e(3)=239.98182d+00
      c(1)=0.0622818d+00
      c(2)=0.3713101d+00
      c(3)=0.6864439d+00
      e(4)=320.46089d+00
      e(5)=30.742778d+00
      e(6)=12.612079d+00
      c(4)=-0.1074604d+00
      c(5)=0.7140333d+00
      c(6)=0.3533614d+00
      e(7)=26.234117d+00
      e(8)=4.4599505d+00
      e(9)=1.9839250d+00
      c(7)=-0.2455863d+00
      c(8)=0.7556468d+00
      c(9)=0.3814041d+00
      e(10)=3.1962046d+00
      e(11)=0.60836651d+00
      e(12)=0.26573519d+00
      c(10)=-0.2974230d+00
      c(11)=0.7364992d+00
      c(12)=0.4298103d+00
      e(13)=0.33420177d+00
      e(14)=0.05169948d+00
      e(15)=0.02174977d+00
      c(13)=0.1961392d+00
      c(14)=-0.6888851d+00
      c(15)=-0.4263861d+00
      e(16)=434.56915d+00
      e(17)=100.70360d+00
      e(18)=29.625034d+00
      c(16)=0.0882772d+00
      c(17)=0.4395651d+00
      c(18)=0.6220621d+00
      e(19)=11.572090d+00
      e(20)=4.5937454d+00
      e(21)=1.9060808d+00
      c(19)=0.2829866d+00
      c(20)=0.5896899d+00
      c(21)=0.2263017d+00
      e(22)=0.95506174d+00
      e(23)=0.39453919d+00
      e(24)=0.15759333d+00
      c(22)=0.4130687d+00
      c(23)=0.5407330d+00
      c(24)=0.1403073d+00
      e(25)=44.566115d+00
      e(26)=11.881489d+00
      e(27)=3.3875579d+00
      c(25)=0.1451271d+00
      c(26)=0.5130677d+00
      c(27)=0.5676640d+00
c
      else if(nucz.le.48) then
       call htmtwo(e,c,nucz)
      else 
      ibr=nucz-48
      go to (300,400,500,600,700,800), ibr
c
c           indium
c
  300 continue
      e(1)=12214.547d+00
      e(2)=1848.9136d+00
      e(3)=406.36833d+00
      c(1)=0.0612476d+00
      c(2)=0.3676754d+00
      c(3)=0.6901359d+00
      e(4)=541.34239d+00
      e(5)=53.018379d+00
      e(6)=22.190283d+00
      c(4)=-0.1096575d+00
      c(5)=0.7156786d+00
      c(6)=0.3516058d+00
      e(7)=46.414454d+00
      e(8)=8.3531661d+00
      e(9)=3.8649088d+00
      c(7)=-0.2690027d+00
      c(8)=0.8333355d+00
      c(9)=0.3162081d+00
      e(10)=7.1775944d+00
      e(11)=1.621446d+00
      e(12)=0.71996916d+00
      c(10)=-0.3439307d+00
      c(11)=0.7769627d+00
      c(12)=0.4273271d+00
      e(13)=1.0654870d+00
      e(14)=0.14431502d+00
      e(15)=0.05458735d+00
      c(13)=-0.2008092d+00
      c(14)=0.6599925d+00
      c(15)=0.4553545d+00
      e(16)=781.86921d+00
      e(17)=182.79818d+00
      e(18)=54.211142d+00
      c(16)=0.0821886d+00
      c(17)=0.4267795d+00
      c(18)=0.6357629d+00
      e(19)=114.03634d+00
      e(20)=18.713126d+00
      e(21)=7.1660478d+00
      c(19)=-0.0311407d+00
      c(20)=0.4133453d+00
      c(21)=0.6578615d+00
      e(22)=2.79637815d+00
      e(23)=1.1841215d+00
      e(24)=0.48773463d+00
      c(22)=0.4068683d+00
      c(23)=0.5479905d+00
      c(24)=0.1323252d+00
      e(25)=0.22993404d+00
      e(26)=0.08355539d+00
      e(27)=0.03249889d+00
      c(25)=0.2630070d+00
      c(26)=0.5733118d+00
      c(27)=0.2829178d+00
      e(28)=102.17356d+00
      e(29)=28.394632d+00
      e(30)=8.9248045d+00
      c(28)=0.1205559d+00
      c(29)=0.4884976d+00
      c(30)=0.5850190d+00
      e(31)=4.5353637d+00
      e(32)=1.5371481d+00
      e(33)=0.49949226d+00
      c(31)=0.2508574d+00
      c(32)=0.5693113d+00
      c(33)=0.3840635d+00
      return
c
c           tin
c
  400 continue
      e(1)=12741.674d+00
      e(2)=1928.4692d+00
      e(3)=423.80797d+00
      c(1)=0.0611353d+00
      c(2)=0.3672929d+00
      c(3)=0.6905447d+00
      e(4)=564.68220d+00
      e(5)=55.241831d+00
      e(6)=23.121421d+00
      c(4)=-0.1097488d+00
      c(5)=0.7181120d+00
      c(6)=0.3490973d+00
      e(7)=48.493098d+00
      e(8)=8.7703153d+00
      e(9)=4.0378104d+00
      c(7)=-0.2708984d+00
      c(8)=0.8409064d+00
      c(9)=0.3101030d+00
      e(10)=7.6154955d+00
      e(11)=1.7298929d+00
      e(12)=0.77188815d+00
      c(10)=-0.3475656d+00
      c(11)=0.7937583d+00
      c(12)=0.4132006d+00
      e(13)=1.1976547d+00
      e(14)=0.17303676d+00
      e(15)=0.06695953d+00
      c(13)=-0.2385589d+00
      c(14)=0.7059861d+00
      c(15)=0.4236124d+00
      e(16)=816.77323d+00
      e(17)=191.12772d+00
      e(18)=56.771595d+00
      c(16)=0.0820729d+00
      c(17)=0.4264331d+00
      c(18)=0.6359016d+00
      e(19)=121.24487d+00
      e(20)=19.577904d+00
      e(21)=7.5717131d+00
      c(19)=-0.0315643d+00
      c(20)=0.4151667d+00
      c(21)=0.6549931d+00
      e(22)=3.0474031d+00
      e(23)=1.3022184d+00
      e(24)=0.54462907d+00
      c(22)=0.3982670d+00
      c(23)=0.5557792d+00
      c(24)=0.1301395d+00
      e(25)=0.26752856d+00
      e(26)=0.10418254d+00
      e(27)=0.04162604d+00
      c(25)=0.3210906d+00
      c(26)=0.5659043d+00
      c(27)=0.2218115d+00
      e(28)=108.05630d+00
      e(29)=30.131576d+00
      e(30)=9.5300359d+00
      c(28)=0.1198237d+00
      c(29)=0.4875910d+00
      c(30)=0.5849874d+00
      e(31)=4.9626098d+00
      e(32)=1.7120829d+00
      e(33)=0.57719451d+00
      c(31)=0.2529487d+00
      c(32)=0.5727612d+00
      c(33)=0.3690387d+00
      return
c
c           antimony
c
  500 continue
      e(1)=13289.383d+00
      e(2)=2010.5218d+00
      e(3)=441.69815d+00
      c(1)=0.0609843d+00
      c(2)=0.3668487d+00
      c(3)=0.6910501d+00
      e(4)=588.38790d+00
      e(5)=57.523192d+00
      e(6)=24.036693d+00
      c(4)=-0.1098365d+00
      c(5)=0.7203891d+00
      c(6)=0.3468674d+00
      e(7)=50.787743d+00
      e(8)=9.1473306d+00
      e(9)=4.1723695d+00
      c(7)=-0.2714342d+00
      c(8)=0.8547971d+00
      c(9)=0.2964847d+00
      e(10)=7.9995999d+00
      e(11)=1.8381194d+00
      e(12)=0.83126332d+00
      c(10)=-0.3512204d+00
      c(11)=0.8121005d+00
      c(12)=0.3980049d+00
      e(13)=1.3445260d+00
      e(14)=0.21146841d+00
      e(15)=0.08689142d+00
      c(13)=-0.2714728d+00
      c(14)=0.6682413d+00
      c(15)=0.4725305d+00
      e(16)=860.08072d+00
      e(17)=201.01421d+00
      e(18)=59.657596d+00
      c(16)=0.0809758d+00
      c(17)=0.4241487d+00
      c(18)=0.6387257d+00
      e(19)=122.27109d+00
      e(20)=20.929665d+00
      e(21)=8.0980350d+00
      c(19)=-0.0318112d+00
      c(20)=0.4045588d+00
      c(21)=0.6653749d+00
      e(22)=3.3003687d+00
      e(23)=1.4287737d+00
      e(24)=0.59505134d+00
      c(22)=0.3889519d+00
      c(23)=0.5640495d+00
      c(24)=0.1289287d+00
      e(25)=0.32801846d+00
      e(26)=0.13661293d+00
      e(27)=0.05411368d+00
      c(25)=0.3110910d+00
      c(26)=0.5652580d+00
      c(27)=0.2257358d+00
      e(28)=115.80955d+00
      e(29)=32.305835d+00
      e(30)=10.250328d+00
      c(28)=0.1166279d+00
      c(29)=0.4834363d+00
      c(30)=0.5901395d+00
      e(31)=5.4862102d+00
      e(32)=1.9216196d+00
      e(33)=0.66606265d+00
      c(31)=0.2483656d+00
      c(32)=0.5743154d+00
      c(33)=0.3643044d+00
      return
c
c           tellurium
c
  600 continue
      e(1)=13796.560d+00
      e(2)=2088.8798d+00
      e(3)=459.39319d+00
      c(1)=0.0610862d+00
      c(2)=0.3669629d+00
      c(3)=0.6907944d+00
      e(4)=611.58438d+00
      e(5)=60.221194d+00
      e(6)=25.491369d+00
      c(4)=-0.1101441d+00
      c(5)=0.7134762d+00
      c(6)=0.3534463d+00
      e(7)=53.042328d+00
      e(8)=9.6632025d+00
      e(9)=4.4188576d+00
      c(7)=-0.2728363d+00
      c(8)=0.8475640d+00
      c(9)=0.3052385d+00
      e(10)=8.4582805d+00
      e(11)=1.9753929d+00
      e(12)=0.90818331d+00
      c(10)=-0.3587862d+00
      c(11)=0.8112858d+00
      c(12)=0.4033116d+00
      e(13)=1.4270256d+00
      e(14)=0.23880222d+00
      e(15)=0.09319910d+00
      c(13)=-0.2806370d+00
      c(14)=0.7509597d+00
      c(15)=0.4038957d+00
      e(16)=895.07497d+00
      e(17)=209.57862d+00
      e(18)=62.405210d+00
      c(16)=0.0811099d+00
      c(17)=0.4236824d+00
      c(18)=0.6386030d+00
      e(19)=131.44799d+00
      e(20)=21.663180d+00
      e(21)=8.4667304d+00
      c(19)=-0.0320210d+00
      c(20)=0.4119187d+00
      c(21)=0.6571079d+00
      e(22)=3.4281d+00
      e(23)=1.4943404d+00
      e(24)=0.63848516d+00
      c(22)=0.4159312d+00
      c(23)=0.5489605d+00
      c(24)=0.1133174d+00
      e(25)=0.36614463d+00
      e(26)=0.14887871d+00
      e(27)=0.05976329d+00
      c(25)=0.3690254d+00
      c(26)=0.5355049d+00
      c(27)=0.2011813d+00
      e(28)=121.40830d+00
      e(29)=34.015217d+00
      e(30)=10.869138d+00
      c(28)=0.1169136d+00
      c(29)=0.4835555d+00
      c(30)=0.5883864d+00
      e(31)=5.8031113d+00
      e(32)=2.0580658d+00
      e(33)=0.73283021d+00
      c(31)=0.2601939d+00
      c(32)=0.5797758d+00
      c(33)=0.3405580d+00
      return
c
c           iodine
c
  700 continue
      e(1)=14351.186d+00
      e(2)=2173.0741d+00
      e(3)=477.87205d+00
      c(1)=0.0610028d+00
      c(2)=0.3666398d+00
      c(3)=0.6911306d+00
      e(4)=636.28295d+00
      e(5)=62.368944d+00
      e(6)=26.035562d+00
      c(4)=-0.1100819d+00
      c(5)=0.7219464d+00
      c(6)=0.3454775d+00
      e(7)=55.135577d+00
      e(8)=10.126452d+00
      e(9)=4.5928632d+00
      c(7)=-0.2754113d+00
      c(8)=0.8543504d+00
      c(9)=0.3006473d+00
      e(10)=8.9049791d+00
      e(11)=2.0925689d+00
      e(12)=0.9735899d+00
      c(10)=-0.3630477d+00
      c(11)=0.8279763d+00
      c(12)=0.3894209d+00
      e(13)=1.5539994d+00
      e(14)=0.27488087d+00
      e(15)=0.10826672d+00
      c(13)=0.2967653d+00
      c(14)=-0.7583335d+00
      c(15)=-0.4067087d+00
      e(16)=933.48281d+00
      e(17)=218.50559d+00
      e(18)=65.325776d+00
      c(16)=0.0810928d+00
      c(17)=0.4231510d+00
      c(18)=0.6385674d+00
      e(19)=158.48285d+00
      e(20)=21.215135d+00
      e(21)=8.4623149d+00
      c(19)=-0.0283465d+00
      c(20)=0.4622270d+00
      c(21)=0.6046503d+00
      e(22)=3.6295595d+00
      e(23)=1.6133740d+00
      e(24)=0.70733879d+00
      c(22)=0.4210620d+00
      c(23)=0.5427157d+00
      c(24)=0.1108784d+00
      e(25)=0.42307578d+00
      e(26)=0.17128391d+00
      e(27)=0.06847909d+00
      c(25)=0.3814674d+00
      c(26)=0.5321926d+00
      c(27)=0.1922032d+00
      e(28)=128.09026d+00
      e(29)=35.982378d+00
      e(30)=11.551116d+00
      c(28)=0.1158636d+00
      c(29)=0.4820494d+00
      c(30)=0.5894448d+00
      e(31)=6.1461523d+00
      e(32)=2.2209370d+00
      e(33)=0.80991202d+00
      c(31)=0.2681817d+00
      c(32)=0.5800614d+00
      c(33)=0.3262263d+00
      return
c
c           xenon
c
  800 continue
      e(1)=14902.236d+00
      e(2)=2256.5383d+00
      e(3)=496.37317d+00
      c(1)=0.0609969d+00
      c(2)=0.3666290d+00
      c(3)=0.6911155d+00
      e(4)=660.94947d+00
      e(5)=64.863039d+00
      e(6)=27.111403d+00
      c(4)=-0.1102247d+00
      c(5)=0.7225832d+00
      c(6)=0.3448171d+00
      e(7)=57.515261d+00
      e(8)=10.536085d+00
      e(9)=4.7042710d+00
      c(7)=-0.2761463d+00
      c(8)=0.8678121d+00
      c(9)=0.2879474d+00
      e(10)=9.3786082d+00
      e(11)=2.2401181d+00
      e(12)=1.0615968d+00
      c(10)=-0.3673742d+00
      c(11)=0.8242752d+00
      c(12)=0.3963011d+00
      e(13)=1.6958079d+00
      e(14)=0.30651180d+00
      e(15)=0.12049354d+00
      c(13)=-0.3127449d+00
      c(14)=0.7841594d+00
      c(15)=0.3890762d+00
      e(16)=971.98851d+00
      e(17)=227.72243d+00
      e(18)=68.162563d+00
      c(16)=0.0809562d+00
      c(17)=0.4227290d+00
      c(18)=0.6388428d+00
      e(19)=164.53820d+00
      e(20)=22.210478d+00
      e(21)=8.90000556d+00
      c(19)=-0.0285155d+00
      c(20)=0.4626758d+00
      c(21)=0.6037699d+00
      e(22)=3.7649228d+00
      e(23)=1.6522375d+00
      e(24)=0.67649625d+00
      c(22)=0.4503453d+00
      c(23)=0.5400257d+00
      c(24)=0.0828350d+00
      e(25)=0.47842852d+00
      e(26)=0.19781294d+00
      e(27)=0.07915387d+00
      c(25)=-0.3900001d+00
      c(26)=-0.5233464d+00
      c(27)=-0.1906046d+00
      e(28)=134.91331d+00
      e(29)=37.956387d+00
      e(30)=12.227475d+00
      c(28)=0.1150105d+00
      c(29)=0.4815952d+00
      c(30)=0.5896133d+00
      e(31)=6.6004928d+00
      e(32)=2.3980513d+00
      e(33)=0.88648239d+00
      c(31)=0.2718844d+00
      c(32)=0.5855569d+00
      c(33)=0.3127456d+00
      endif
      return
      end
      subroutine hsix(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
c     huzinaga sixth row main group elements
c
      if(nucz.eq.55) then
c
c     cesium
c
      e(1)= 1.5525866d+04
      e(2)= 2.3490055d+03
      e(3)= 5.1623555d+02
      c(1)= 6.07242d-02
      c(2)= 3.658416d-01
      c(3)= 6.92047d-01
      e(4)= 6.8750505d+02
      e(5)= 6.7191157d+01
      e(6)= 2.7992003d+01
      c(4)=-1.102127d-01
      c(5)= 7.272034d-01
      c(6)= 3.401429d-01
      e(7)= 1.0242833d+03
      e(8)= 2.3964889d+02
      e(9)= 7.1255548d+01
      c(7)= 7.91852d-02
      c(8)= 4.199333d-01
      c(9)= 6.432514d-01
      e(10)= 6.0168188d+01
      e(11)= 1.0868694d+01
      e(12)= 4.7407047d+00
      c(10)=-2.753151d-01
      c(11)= 8.901281d-01
      c(12)= 2.648633d-01
      e(13)= 1.3945864d+02
      e(14)= 2.581081d+01
      e(15)= 1.0139452d+01
      c(13)=-3.31228d-02
      c(14)= 3.891083d-01
      c(15)= 6.790325d-01
      e(16)= 1.4405139d+02
      e(17)= 4.0505146d+01
      e(18)= 1.3086244d+01
      c(16)= 1.118302d-01
      c(17)= 4.768349d-01
      c(18)= 5.954174d-01
      e(19)= 9.8871626d+00
      e(20)= 2.3693425d+00
      e(21)= 1.1416761d+00
      c(19)=-3.688268d-01
      c(20)= 8.325824d-01
      c(21)= 3.887891d-01
      e(22)= 4.3242445d+00
      e(23)= 1.9364152d+00
      e(24)= 8.145053d-01
      c(22)= 3.816688d-01
      c(23)= 5.769763d-01
      c(24)= 1.154047d-01
      e(25)= 7.3715777d+00
      e(26)= 2.7248237d+00
      e(27)= 1.0266781d+00
      c(25)= 2.537255d-01
      c(26)= 5.78886d-01
      c(27)= 3.325981d-01
      e(28)= 1.7690784d+00
      e(29)= 3.8362031d-01
      e(30)= 1.7130634d-01
      c(28)=-3.409418d-01
      c(29)= 7.000971d-01
      c(30)= 4.934524d-01
      e(31)= 5.9458425d-01
      e(32)= 2.7499182d-01
      e(33)= 1.115529d-01
      c(31)= 3.199069d-01
      c(32)= 5.531867d-01
      c(33)= 2.160885d-01
      e(34)= 2.1184878d-01
      e(35)= 4.032421d-02
      e(36)= 1.705814d-02
      c(34)=-2.750217d-01
      c(35)= 1.93604d-01
      c(36)= 9.141935d-01
c
      else if(nucz.eq.56) then
c
c     barium
c
      e(1)= 1.6098886d+04
      e(2)= 2.436296d+03
      e(3)= 5.3562078d+02
      c(1)= 6.07199d-02
      c(2)= 3.657475d-01
      c(3)= 6.920963d-01
      e(4)= 7.1309716d+02
      e(5)= 6.9806392d+01
      e(6)= 2.9031078d+01
      c(4)=-1.103332d-01
      c(5)= 7.275923d-01
      c(6)= 3.399361d-01
      e(7)= 1.0653684d+03
      e(8)= 2.4949087d+02
      e(9)= 7.4311477d+01
      c(7)= 7.90271d-02
      c(8)= 4.191952d-01
      c(9)= 6.43758d-01
      e(10)= 6.2480417d+01
      e(11)= 1.1412912d+01
      e(12)= 4.9805876d+00
      c(10)=-2.773775d-01
      c(11)= 8.867966d-01
      c(12)= 2.701325d-01
      e(13)= 1.4706367d+02
      e(14)= 2.679279d+01
      e(15)= 1.0604816d+01
      c(13)=-3.33375d-02
      c(14)= 3.922373d-01
      c(15)= 6.75227d-01
      e(16)= 1.5169168d+02
      e(17)= 4.2753274d+01
      e(18)= 1.3850108d+01
      c(16)= 1.105491d-01
      c(17)= 4.750537d-01
      c(18)= 5.972132d-01
      e(19)= 1.0265539d+01
      e(20)= 2.5292321d+00
      e(21)= 1.2273505d+00
      c(19)=-3.77872d-01
      c(20)= 8.405296d-01
      c(21)= 3.887798d-01
      e(22)= 4.433544d+00
      e(23)= 1.9758409d+00
      e(24)= 8.1191629d-01
      c(22)= 4.153937d-01
      c(23)= 5.642885d-01
      c(24)= 9.2372d-02
      e(25)= 7.7435537d+00
      e(26)= 2.8709772d+00
      e(27)= 1.0976557d+00
      c(25)= 2.643998d-01
      c(26)= 5.858296d-01
      c(27)= 3.112276d-01
      e(28)= 1.9672783d+00
      e(29)= 4.1058204d-01
      e(30)= 1.8712746d-01
      c(28)=-3.632895d-01
      c(29)= 8.040338d-01
      c(30)= 3.967839d-01
      e(31)= 6.1380042d-01
      e(32)= 2.7282586d-01
      e(33)= 1.1006559d-01
      c(31)= 4.303242d-01
      c(32)= 5.355877d-01
      c(33)= 1.136744d-01
      e(34)= 2.1270737d-01
      e(35)= 4.180388d-02
      e(36)= 1.768408d-02
      c(34)=-2.013774d-01
      c(35)= 7.154536d-01
      c(36)= 4.176434d-01
c
      else if(nucz.le.70) then
c
         call hlanth(e,c,nucz)
c
      else if(nucz.le.80) then
c
         call htmthr(e,c,nucz)
c
      else
c
      ihop=nucz-80
      go to (300,400,500,600,700,800), ihop
c
c     thallium
c
  300 continue
      e(1)= 3.4131365d+04
      e(2)= 5.1723898d+03
      e(3)= 1.1401778d+03
      c(1)= 6.0001d-02
      c(2)= 3.626063d-01
      c(3)= 6.949969d-01
      e(4)= 1.5117882d+03
      e(5)= 1.5097546d+02
      e(6)= 6.4320792d+01
      c(4)=-1.127989d-01
      c(5)= 7.369186d-01
      c(6)= 3.301527d-01
      e(7)= 2.3775987d+03
      e(8)= 5.6119485d+02
      e(9)= 1.697312d+02
      c(7)= 7.53926d-02
      c(8)= 4.099688d-01
      c(9)= 6.518046d-01
      e(10)= 1.4010346d+02
      e(11)= 2.7595196d+01
      e(12)= 1.1845526d+01
      c(10)=-2.971921d-01
      c(11)= 9.065841d-01
      c(12)= 2.669625d-01
      e(13)= 3.1631847d+02
      e(14)= 6.4888544d+01
      e(15)= 2.673353d+01
      c(13)=-3.74019d-02
      c(14)= 3.911708d-01
      c(15)= 6.743423d-01
      e(16)= 3.8360621d+02
      e(17)= 1.110916d+02
      e(18)= 3.7869799d+01
      c(16)= 1.021961d-01
      c(17)= 4.68311d-01
      c(18)= 5.976181d-01
      e(19)= 2.6902991d+01
      e(20)= 6.6871496d+00
      e(21)= 3.3176591d+00
      c(19)=-3.793081d-01
      c(20)= 8.769824d-01
      c(21)= 3.550289d-01
      e(22)= 1.2328211d+01
      e(23)= 5.5114323d+00
      e(24)= 2.2169643d+00
      c(22)= 4.360675d-01
      c(23)= 5.574551d-01
      c(24)= 7.5431d-02
      e(25)= 2.2556467d+01
      e(26)= 9.0361081d+00
      e(27)= 3.8880752d+00
      c(25)= 2.872732d-01
      c(26)= 5.606543d-01
      c(27)= 2.917691d-01
      e(28)= 5.6851804d+01
      e(29)= 1.7181822d+01
      e(30)= 5.3079932d+00
      c(28)= 1.704893d-01
      c(29)= 5.344428d-01
      c(30)= 5.513373d-01
      e(31)= 5.6090779d+00
      e(32)= 1.2974407d+00
      e(33)= 5.5988532d-01
      c(31)=-3.948161d-01
      c(32)= 8.823134d-01
      c(33)= 3.507943d-01
      e(34)= 1.993276d+00
      e(35)= 8.9758441d-01
      e(36)= 3.7280215d-01
      c(34)= 4.847093d-01
      c(35)= 4.987045d-01
      c(36)= 8.95761d-02
      e(37)= 2.5003711d+00
      e(38)= 9.4037699d-01
      e(39)= 3.2314233d-01
      c(37)= 3.399404d-01
      c(38)= 5.676636d-01
      c(39)= 2.723492d-01
      e(40)= 9.2812701d-01
      e(41)= 1.1022919d-01
      e(42)= 4.05724d-02
      c(40)=-1.732925d-01
      c(41)= 8.697721d-01
      c(42)= 2.258045d-01
      e(43)= 2.2282221d-01
      e(44)= 9.313363d-02
      e(45)= 3.725345d-02
      c(43)= 2.145554d-01
      c(44)= 5.125545d-01
      c(45)= 3.817007d-01
      return
c
c     lead
c
  400 continue
      e(1)= 3.507524d+04
      e(2)= 5.3112719d+03
      e(3)= 1.1700875d+03
      c(1)= 5.98351d-02
      c(2)= 3.621904d-01
      c(3)= 6.955115d-01
      e(4)= 1.5518174d+03
      e(5)= 1.5431441d+02
      e(6)= 6.548108d+01
      c(4)=-1.12727d-01
      c(5)= 7.418872d-01
      c(6)= 3.251188d-01
      e(7)= 2.4494365d+03
      e(8)= 5.7781442d+02
      e(9)= 1.7467073d+02
      c(7)= 7.49361d-02
      c(8)= 4.089706d-01
      c(9)= 6.530269d-01
      e(10)= 1.4395513d+02
      e(11)= 2.8198551d+01
      e(12)= 1.2139813d+01
      c(10)=-2.977604d-01
      c(11)= 9.161297d-01
      c(12)= 2.568361d-01
      e(13)= 3.2439605d+02
      e(14)= 6.7056357d+01
      e(15)= 2.7668571d+01
      c(13)=-3.77804d-02
      c(14)= 3.872369d-01
      c(15)= 6.780358d-01
      e(16)= 3.9602788d+02
      e(17)= 1.1478406d+02
      e(18)= 3.9165515d+01
      c(16)= 1.015518d-01
      c(17)= 4.672091d-01
      c(18)= 5.988161d-01
      e(19)= 2.7546889d+01
      e(20)= 6.8739521d+00
      e(21)= 3.4374831d+00
      c(19)=-3.839922d-01
      c(20)= 8.909138d-01
      c(21)= 3.440137d-01
      e(22)= 1.3013204d+01
      e(23)= 5.8680654d+00
      e(24)= 2.2720421d+00
      c(22)= 4.135703d-01
      c(23)= 5.716137d-01
      c(24)= 8.58261d-02
      e(25)= 2.3480294d+01
      e(26)= 9.4592198d+00
      e(27)= 4.1021059d+00
      c(25)= 2.843855d-01
      c(26)= 5.58818d-01
      c(27)= 2.94868d-01
      e(28)= 5.9277049d+01
      e(29)= 1.7991863d+01
      e(30)= 5.6059964d+00
      c(28)= 1.703522d-01
      c(29)= 5.352027d-01
      c(30)= 5.483777d-01
      e(31)= 6.0560593d+00
      e(32)= 1.2999385d+00
      e(33)= 5.8010673d-01
      c(31)=-3.675232d-01
      c(32)= 9.261441d-01
      c(33)= 2.848743d-01
      e(34)= 2.132532d+00
      e(35)= 9.8157824d-01
      e(36)= 3.9467207d-01
      c(34)= 4.708301d-01
      c(35)= 5.131364d-01
      c(36)= 8.62959d-02
      e(37)= 2.6854228d+00
      e(38)= 1.0321544d+00
      e(39)= 3.6406946d-01
      c(37)= 3.437064d-01
      c(38)= 5.720556d-01
      c(39)= 2.539847d-01
      e(40)= 1.0041909d+00
      e(41)= 1.5157282d-01
      e(42)= 4.158037d-02
      c(40)=-2.345531d-01
      c(41)= 8.702589d-01
      c(42)= 2.826179d-01
      e(43)= 2.3589377d-01
      e(44)= 9.702808d-02
      e(45)= 3.817897d-02
      c(43)= 3.212519d-01
      c(44)= 5.604186d-01
      c(45)= 2.230238d-01
      return
c
c     bisthmus
c
  500 continue
      e(1)= 3.596366d+04
      e(2)= 5.4454286d+03
      e(3)= 1.1996058d+03
      c(1)= 5.97899d-02
      c(2)= 3.620306d-01
      c(3)= 6.956797d-01
      e(4)= 1.5912502d+03
      e(5)= 1.5816746d+02
      e(6)= 6.6985958d+01
      c(4)=-1.127545d-01
      c(5)= 7.427995d-01
      c(6)= 3.2429d-01
      e(7)= 2.5172545d+03
      e(8)= 5.9388843d+02
      e(9)= 1.7956047d+02
      c(7)= 7.46957d-02
      c(8)= 4.083137d-01
      c(9)= 6.537242d-01
      e(10)= 1.4755203d+02
      e(11)= 2.9134674d+01
      e(12)= 1.2498034d+01
      c(10)=-2.984864d-01
      c(11)= 9.096834d-01
      c(12)= 2.646867d-01
      e(13)= 3.3093361d+02
      e(14)= 6.907541d+01
      e(15)= 2.8522477d+01
      c(13)=-3.79046d-02
      c(14)= 3.862204d-01
      c(15)= 6.790555d-01
      e(16)= 4.0774836d+02
      e(17)= 1.1836555d+02
      e(18)= 4.0481465d+01
      c(16)= 1.012945d-01
      c(17)= 4.662523d-01
      c(18)= 5.99408d-01
      e(19)= 2.8216255d+01
      e(20)= 7.1495803d+00
      e(21)= 3.5218354d+00
      c(19)=-3.886382d-01
      c(20)= 8.948072d-01
      c(21)= 3.451839d-01
      e(22)= 1.3329817d+01
      e(23)= 6.0302971d+00
      e(24)= 2.4456258d+00
      c(22)= 4.230622d-01
      c(23)= 5.634345d-01
      c(24)= 8.20663d-02
      e(25)= 2.4322803d+01
      e(26)= 9.8293386d+00
      e(27)= 4.2764203d+00
      c(25)= 2.836262d-01
      c(26)= 5.599295d-01
      c(27)= 2.933151d-01
      e(28)= 6.1976718d+01
      e(29)= 1.8871879d+01
      e(30)= 5.9233061d+00
      c(28)= 1.693912d-01
      c(29)= 5.355429d-01
      c(30)= 5.466911d-01
      e(31)= 6.2636009d+00
      e(32)= 1.4271684d+00
      e(33)= 6.1662978d-01
      c(31)=-3.81391d-01
      c(32)= 9.092481d-01
      c(33)= 3.171021d-01
      e(34)= 2.2105424d+00
      e(35)= 1.0056903d+00
      e(36)= 4.0840048d-01
      c(34)= 5.015939d-01
      c(35)= 4.977996d-01
      c(36)= 6.7302d-02
      e(37)= 2.7618843d+00
      e(38)= 1.0712666d+00
      e(39)= 3.8408917d-01
      c(37)= 3.726898d-01
      c(38)= 5.718582d-01
      c(39)= 2.160448d-01
      e(40)= 1.0835409d+00
      e(41)= 1.6361768d-01
      e(42)= 4.260071d-02
      c(40)=-2.66303d-01
      c(41)= 9.547685d-01
      c(42)= 2.038896d-01
      e(43)= 2.4893283d-01
      e(44)= 9.886529d-02
      e(45)= 3.911584d-02
      c(43)= 4.353007d-01
      c(44)= 5.414271d-01
      c(45)= 1.184326d-01
      return
c
c     polonium
c
  600 continue
      e(1)= 3.6867892d+04
      e(2)= 5.5815605d+03
      e(3)= 1.2295467d+03
      c(1)= 5.97396d-02
      c(2)= 3.618671d-01
      c(3)= 6.95855d-01
      e(4)= 1.6313133d+03
      e(5)= 1.618207d+02
      e(6)= 6.7853196d+01
      c(4)=-1.126817d-01
      c(5)= 7.46108d-01
      c(6)= 3.213951d-01
      e(7)= 2.5859828d+03
      e(8)= 6.1013029d+02
      e(9)= 1.8448204d+02
      c(7)= 7.44629d-02
      c(8)= 4.077334d-01
      c(9)= 6.543629d-01
      e(10)= 1.5133201d+02
      e(11)= 2.9857813d+01
      e(12)= 1.2666596d+01
      c(10)=-2.989084d-01
      c(11)= 9.149115d-01
      c(12)= 2.600228d-01
      e(13)= 3.4055674d+02
      e(14)= 7.1276662d+01
      e(15)= 2.9491315d+01
      c(13)=-3.82561d-02
      c(14)= 3.823956d-01
      c(15)= 6.825233d-01
      e(16)= 4.217038d+02
      e(17)= 1.2238581d+02
      e(18)= 4.1880233d+01
      c(16)= 1.003116d-01
      c(17)= 4.647457d-01
      c(18)= 6.01276d-01
      e(19)= 2.9409549d+01
      e(20)= 7.2021546d+00
      e(21)= 3.5479184d+00
      c(19)=-3.820933d-01
      c(20)= 9.275939d-01
      c(21)= 3.053887d-01
      e(22)= 1.4029785d+01
      e(23)= 6.4245656d+00
      e(24)= 2.6492696d+00
      c(22)= 4.029453d-01
      c(23)= 5.715225d-01
      c(24)= 9.42772d-02
      e(25)= 2.5183871d+01
      e(26)= 1.0214589d+01
      e(27)= 4.4950523d+00
      c(25)= 2.83553d-01
      c(26)= 5.580633d-01
      c(27)= 2.934985d-01
      e(28)= 6.4795841d+01
      e(29)= 1.9791217d+01
      e(30)= 6.2541857d+00
      c(28)= 1.681614d-01
      c(29)= 5.355027d-01
      c(30)= 5.456886d-01
      e(31)= 6.4693577d+00
      e(32)= 1.4613166d+00
      e(33)= 6.449106d-01
      c(31)=-3.830339d-01
      c(32)= 9.531587d-01
      c(33)= 2.725793d-01
      e(34)= 2.3710859d+00
      e(35)= 1.1111472d+00
      e(36)= 4.1830074d-01
      c(34)= 4.781188d-01
      c(35)= 5.168473d-01
      c(36)= 7.11346d-02
      e(37)= 3.0692218d+00
      e(38)= 1.219541d+00
      e(39)= 4.4471254d-01
      c(37)= 3.500265d-01
      c(38)= 5.809209d-01
      c(39)= 2.231577d-01
      e(40)= 1.2428201d+00
      e(41)= 1.8277504d-01
      e(42)= 4.363342d-02
      c(40)=-2.907553d-01
      c(41)= 9.91198d-01
      c(42)= 1.728911d-01
      e(43)= 2.7680626d-01
      e(44)= 1.0126194d-01
      e(45)= 4.006407d-02
      c(43)= 4.983952d-01
      c(44)= 5.216655d-01
      c(45)= 7.70469d-02
      return
c
c     astatine
c
  700 continue
      e(1)= 3.7776571d+04
      e(2)= 5.7185732d+03
      e(3)= 1.2597282d+03
      c(1)= 5.97016d-02
      c(2)= 3.617435d-01
      c(3)= 6.959848d-01
      e(4)= 1.6696097d+03
      e(5)= 1.658396d+02
      e(6)= 6.9717694d+01
      c(4)=-1.127567d-01
      c(5)= 7.462998d-01
      c(6)= 3.211019d-01
      e(7)= 2.6553837d+03
      e(8)= 6.2633825d+02
      e(9)= 1.8939113d+02
      c(7)= 7.42646d-02
      c(8)= 4.073604d-01
      c(9)= 6.548009d-01
      e(10)= 1.5565678d+02
      e(11)= 3.0520449d+01
      e(12)= 1.2699037d+01
      c(10)=-2.973679d-01
      c(11)= 9.201904d-01
      c(12)= 2.542612d-01
      e(13)= 3.4649889d+02
      e(14)= 7.3495813d+01
      e(15)= 3.0402772d+01
      c(13)=-3.8432d-02
      c(14)= 3.80631d-01
      c(15)= 6.843942d-01
      e(16)= 4.3427898d+02
      e(17)= 1.261555d+02
      e(18)= 4.3240302d+01
      c(16)= 9.99131d-02
      c(17)= 4.638397d-01
      c(18)= 6.020334d-01
      e(19)= 3.0256338d+01
      e(20)= 7.4299444d+00
      e(21)= 3.6548245d+00
      c(19)=-3.835086d-01
      c(20)= 9.334748d-01
      c(21)= 3.007817d-01
      e(22)= 1.4334294d+01
      e(23)= 6.5437078d+00
      e(24)= 2.6726334d+00
      c(22)= 4.138476d-01
      c(23)= 5.676362d-01
      c(24)= 8.65671d-02
      e(25)= 2.6081763d+01
      e(26)= 1.0651015d+01
      e(27)= 4.7470935d+00
      c(25)= 2.823875d-01
      c(26)= 5.537894d-01
      c(27)= 2.970252d-01
      e(28)= 6.7670974d+01
      e(29)= 2.0720618d+01
      e(30)= 6.5872235d+00
      c(28)= 1.670327d-01
      c(29)= 5.357835d-01
      c(30)= 5.445086d-01
      e(31)= 6.799067d+00
      e(32)= 1.5270767d+00
      e(33)= 6.6035701d-01
      c(31)=-3.767089d-01
      c(32)= 9.732817d-01
      c(33)= 2.500856d-01
      e(34)= 2.4797553d+00
      e(35)= 1.1577118d+00
      e(36)= 4.2831956d-01
      c(34)= 4.940629d-01
      c(35)= 5.099957d-01
      c(36)= 5.96195d-02
      e(37)= 3.1664963d+00
      e(38)= 1.2549894d+00
      e(39)= 4.7297912d-01
      c(37)= 3.785868d-01
      c(38)= 5.763374d-01
      c(39)= 1.915263d-01
      e(40)= 1.3540229d+00
      e(41)= 2.0700628d-01
      e(42)= 4.467849d-02
      c(40)=-2.865345d-01
      c(41)= 1.010233d+00
      c(42)= 1.589027d-01
      e(43)= 2.9831651d-01
      e(44)= 1.0695344d-01
      e(45)= 4.102365d-02
      c(43)= 5.541128d-01
      c(44)= 4.849933d-01
      c(45)= 5.59909d-02
      return
c
c     radon
c
  800 continue
      e(1)= 3.8634492d+04
      e(2)= 5.852336d+03
      e(3)= 1.2897187d+03
      c(1)= 5.97563d-02
      c(2)= 3.617782d-01
      c(3)= 6.958837d-01
      e(4)= 1.710508d+03
      e(5)= 1.6955891d+02
      e(6)= 7.0548371d+01
      c(4)=-1.126837d-01
      c(5)= 7.496104d-01
      c(6)= 3.182236d-01
      e(7)= 2.7332176d+03
      e(8)= 6.4418252d+02
      e(9)= 1.9494224d+02
      c(7)= 7.3863d-02
      c(8)= 4.061085d-01
      c(9)= 6.56063d-01
      e(10)= 1.5892829d+02
      e(11)= 3.1748256d+01
      e(12)= 1.4197685d+01
      c(10)=-3.022314d-01
      c(11)= 9.042464d-01
      c(12)= 2.711234d-01
      e(13)= 3.6084753d+02
      e(14)= 7.7080714d+01
      e(15)= 3.1547775d+01
      c(13)=-3.86689d-02
      c(14)= 3.67417d-01
      c(15)= 6.977983d-01
      e(16)= 4.5168521d+02
      e(17)= 1.311073d+02
      e(18)= 4.4898191d+01
      c(16)= 9.79008d-02
      c(17)= 4.603345d-01
      c(18)= 6.066244d-01
      e(19)= 3.0871286d+01
      e(20)= 7.9177164d+00
      e(21)= 3.9643042d+00
      c(19)=-3.929887d-01
      c(20)= 8.867217d-01
      c(21)= 3.562264d-01
      e(22)= 5.5835907d+01
      e(23)= 1.2265749d+01
      e(24)= 5.6389792d+00
      c(22)=-1.14409d-01
      c(23)= 5.688023d-01
      c(24)= 5.191901d-01
      e(25)= 2.6989804d+01
      e(26)= 1.0979775d+01
      e(27)= 4.8990999d+00
      c(25)= 2.829638d-01
      c(26)= 5.586436d-01
      c(27)= 2.90951d-01
      e(28)= 7.1255348d+01
      e(29)= 2.1844942d+01
      e(30)= 6.9772755d+00
      c(28)= 1.637813d-01
      c(29)= 5.341357d-01
      c(30)= 5.469714d-01
      e(31)= 6.6237042d+00
      e(32)= 1.7190412d+00
      e(33)= 7.9915287d-01
      c(31)=-4.3546d-01
      c(32)= 9.268309d-01
      c(33)= 3.407903d-01
      e(34)= 2.50277d+00
      e(35)= 1.1043323d+00
      e(36)= 4.1119d-01
      c(34)= 5.662413d-01
      c(35)= 4.638213d-01
      c(36)= 3.02746d-02
      e(37)= 3.6329003d+00
      e(38)= 1.5131756d+00
      e(39)= 5.8490621d-01
      c(37)= 3.218022d-01
      c(38)= 5.8162d-01
      c(39)= 2.380017d-01
      e(40)= 1.342131d+00
      e(41)= 2.6767369d-01
      e(42)= 1.0544732d-01
      c(40)=-3.405358d-01
      c(41)= 8.102372d-01
      c(42)= 3.844879d-01
      e(43)= 3.81996d-01
      e(44)= 1.60392d-01
      e(45)= 6.61157d-02
      c(43)= 4.321029d-01
      c(44)= 5.050096d-01
      c(45)= 1.595361d-01
c
      endif
c
      return
      end
      subroutine htmone(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      ibr=nucz-20
      go to (100,200,300,400,500,600,700,800,900,1000), ibr
c
c     scandium (4-f)
c
  100 continue
      e(1)=2121.1153d+00
      e(2)=320.52014d+00
      e(3)=69.895734d+00
      c(1)=0.0643811d+00
      c(2)=0.3791356d+00
      c(3)=0.6790384d+00
      e(4)=92.827788d+00
      e(5)=8.3661954d+00
      e(6)=3.3619693d+00
      c(4)=-0.0997668d+00
      c(5)=0.6751087d+00
      c(6)=0.3908713d+00
      e(7)=6.1732210d+00
      e(8)=0.95574656d+00
      e(9)=0.37934186d+00
      c(7)=-0.2224250d+00
      c(8)=0.71394d+00
      c(9)=0.4142398d+00
      e(10)=0.45254574d+00
      e(11)=0.05968487d+00
      e(12)=0.02352408d+00
      c(10)=-0.1502541d+00
      c(11)=0.651097d+00
      c(12)=0.4444372d+00
      e(13)=110.16022d+00
      e(14)=24.886185d+00
      e(15)=6.9221188d+00
      c(13)=0.0957608d+00
      c(14)=0.4482512d+00
      c(15)=0.6189104d+00
      e(16)=2.3903679d+00
      e(17)=0.90248233d+00
      e(18)=0.32475569d+00
      c(16)=0.2403846d+00
      c(17)=0.5797456d+00
      c(18)=0.3020693d+00
      e(19)=4.526837d+00
      e(20)=1.0394154d+00
      e(21)=0.21850767d+00
      c(19)=0.1644737d+00
      c(20)=0.4787843d+00
      c(21)=0.6556389d+00
      return
c
c     titanium (5-f)
c
  200 continue
      e( 1)=2336.7153d+00
      e( 2)=353.15331d+00
      e( 3)=77.059812d+00
      c( 1)=0.0641708d+00
      c( 2)=0.3783819d+00
      c( 3)=0.6797627d+00
      e( 4)=102.45710d+00
      e( 5)=9.2769361d+00
      e( 6)=3.7368558d+00
      c( 4)=-0.1005637d+00
      c( 5)=0.6815896d+00
      c( 6)=0.3844102d+00
      e( 7)=7.0008762d+00
      e( 8)=1.0875122d+00
      e( 9)=0.42697745d+00
      c( 7)=-0.2231174d+00
      c( 8)=0.7199595d+00
      c( 9)=0.4094673d+00
      e(10)=0.54148021d+00
      e(11)=0.06648139d+00
      e(12)=0.02576476d+00
      c(10)=-0.1434300d+00
      c(11)=0.6030041d+00
      c(12)=0.4881755d+00
      e(13)=123.33273d+00
      e(14)=27.926218d+00
      e(15)=7.8132163d+00
      c(13)=0.0949569d+00
      c(14)=0.4477326d+00
      c(15)=0.6188790d+00
      e(16)=2.75228040d+00
      e(17)=1.0394584d+00
      e(18)=0.37298279d+00
      c(16)=0.2438232d+00
      c(17)=0.5781427d+00
      c(18)=0.3010109d+00
      e(19)=5.9418882d+00
      e(20)=1.3963895d+00
      e(21)=0.30602356d+00
      c(19)=0.1638319d+00
      c(20)=0.4905779d+00
      c(21)=0.6362293d+00
      return
c
c     vanadium (6-d)
c
  300 continue
      e( 1)=2563.4389d+00
      e( 2)=387.46176d+00
      e( 3)=84.577619d+00
      c( 1)=0.0639578d+00
      c( 2)=0.3776589d+00
      c( 3)=0.6804836d+00
      e( 4)=112.57688d+00
      e( 5)=10.242910d+00
      e( 6)=4.1346585d+00
      c( 4)=-0.1013124d+00
      c( 5)=0.6865452d+00
      c( 6)=0.3795191d+00
      e( 7)=7.8770045d+00
      e( 8)=1.2266238d+00
      e( 9)=0.47757702d+00
      c( 7)=-0.2233367d+00
      c( 8)=0.7226660d+00
      c( 9)=0.4075932d+00
      e(10)=0.61583188d+00
      e(11)=0.07247550d+00
      e(12)=0.02753023d+00
      c(10)=-0.1329102d+00
      c(11)=0.5794076d+00
      c(12)=0.5081949d+00
      e(13)=136.91997d+00
      e(14)=31.081206d+00
      e(15)=8.7448393d+00
      c(13)=0.0945211d+00
      c(14)=0.4477583d+00
      c(15)=0.6181758d+00
      e(16)=3.1310711d+00
      e(17)=1.1822484d+00
      e(18)=0.42352383d+00
      c(16)=0.2473272d+00
      c(17)=0.5762748d+00
      c(18)=0.3001113d+00
      e(19)=7.2698474d+00
      e(20)=1.7314426d+00
      e(21)=0.38675247d+00
      c(19)=0.1642244d+00
      c(20)=0.4967508d+00
      c(21)=0.6258979d+00
      return
c
c     chromium (7-s)
c
  400 continue
      e( 1)=2799.7104d+00
      e( 2)=423.21821d+00
      e( 3)=92.434349d+00
      c( 1)=0.0637916d+00
      c( 2)=0.3770772d+00
      c( 3)=0.6810364d+00
      e( 4)=123.22634d+00
      e( 5)=11.246630d+00
      e( 6)=4.5372546d+00
      c( 4)=-0.1019598d+00
      c( 5)=0.6925012d+00
      c( 6)=0.3737091d+00
      e( 7)=8.7804852d+00
      e( 8)=1.3793596d+00
      e( 9)=0.53502318d+00
      c( 7)=-0.2235841d+00
      c( 8)=0.7175520d+00
      c( 9)=0.4136269d+00
      e(10)=0.72655885d+00
      e(11)=0.07639020d+00
      e(12)=0.02860476d+00
      c(10)=-0.1294804d+00
      c(11)=0.5711774d+00
      c(12)=0.5126118d+00
      e(13)=151.28977d+00
      e(14)=34.410311d+00
      e(15)=9.7258583d+00
      c(13)=0.0940601d+00
      c(14)=0.4477828d+00
      c(15)=0.6176532d+00
      e(16)=3.4776582d+00
      e(17)=1.3079405d+00
      e(18)=0.46753229d+00
      c(16)=0.2581111d+00
      c(17)=0.5757324d+00
      c(18)=0.2908906d+00
      e(19)=8.6262099d+00
      e(20)=2.0754716d+00
      e(21)=0.47058738d+00
      c(19)=0.1647178d+00
      c(20)=0.5011371d+00
      c(21)=0.6180587d+00
      return
c
c     manganese (6-d)
c
  500 continue
      e( 1)=3044.3379d+00
      e( 2)=460.32835d+00
      e( 3)=100.60952d+00
      c( 1)=0.0636857d+00
      c( 2)=0.3766602d+00
      c( 3)=0.6813959d+00
      e( 4)=134.04812d+00
      e( 5)=12.355650d+00
      e( 6)=5.0118328d+00
      c( 4)=-0.1026805d+00
      c( 5)=0.6915740d+00
      c( 6)=0.3747452d+00
      e( 7)=9.7428339d+00
      e( 8)=1.5375154d+00
      e( 9)=0.59408091d+00
      c( 7)=-0.2239616d+00
      c( 8)=0.7171970d+00
      c( 9)=0.4146924d+00
      e(10)=0.76261459d+00
      e(11)=0.07836474d+00
      e(12)=0.02930169d+00
      c(10)=-0.1159459d+00
      c(11)=0.5638312d+00
      c(12)=0.5157105d+00
      e(13)=166.73180d+00
      e(14)=37.977416d+00
      e(15)=10.774333d+00
      c(13)=0.0933798d+00
      c(14)=0.4472796d+00
      c(15)=0.6178799d+00
      e(16)=3.9731007d+00
      e(17)=1.5014460d+00
      e(18)=0.53680586d+00
      c(16)=0.2500751d+00
      c(17)=0.5725861d+00
      c(18)=0.3020932d+00
      e(19)=9.8656162d+00
      e(20)=2.3866048d+00
      e(21)=0.53817952d+00
      c(19)=0.1674847d+00
      c(20)=0.5041187d+00
      c(21)=0.6143492d+00
      return
c
c     iron (5-f)
c
  600 continue
      e( 1)=3302.0297d+00
      e( 2)=499.32535d+00
      e( 3)=109.17229d+00
      c( 1)=0.0635283d+00
      c( 2)=0.3761252d+00
      c( 3)=0.6819206d+00
      e( 4)=145.61613d+00
      e( 5)=13.468079d+00
      e( 6)=5.4624842d+00
      c( 4)=-0.1032386d+00
      c( 5)=0.6955826d+00
      c( 6)=0.3708869d+00
      e( 7)=10.736268d+00
      e( 8)=1.7009576d+00
      e( 9)=0.65364599d+00
      c( 7)=-0.2241245d+00
      c( 8)=0.7177415d+00
      c( 9)=0.4148646d+00
      e(10)=0.81194929d+00
      e(11)=0.07967335d+00
      e(12)=0.02982454d+00
      c(10)=-0.1083327d+00
      c(11)=0.5566047d+00
      c(12)=0.5194924d+00
      e(13)=182.55775d+00
      e(14)=41.647822d+00
      e(15)=11.857806d+00
      c(13)=0.0929956d+00
      c(14)=0.4473347d+00
      c(15)=0.6174492d+00
      e(16)=4.4051684d+00
      e(17)=1.6609361d+00
      e(18)=0.59330383d+00
      c(16)=0.2543695d+00
      c(17)=0.5722388d+00
      c(18)=0.2986741d+00
      e(19)=11.198379d+00
      e(20)=2.7235552d+00
      e(21)=0.61395523d+00
      c(19)=0.1694177d+00
      c(20)=0.5066518d+00
      c(21)=0.6105528d+00
      return
c
c     cobalt (4-f)
c
  700 continue
      e( 1)=3569.1023d+00
      e( 2)=539.75057d+00
      e( 3)=118.06735d+00
      c( 1)=0.0634072d+00
      c( 2)=0.3757053d+00
      c( 3)=0.6823127d+00
      e( 4)=157.60903d+00
      e( 5)=14.633143d+00
      e( 6)=5.9354094d+00
      c( 4)=-0.1037551d+00
      c( 5)=0.6988574d+00
      c( 6)=0.3677667d+00
      e( 7)=11.790403d+00
      e( 8)=1.8690733d+00
      e( 9)=0.71451064d+00
      c( 7)=-0.2238506d+00
      c( 8)=0.7186408d+00
      c( 9)=0.4142824d+00
      e(10)=0.83085499d+00
      e(11)=0.08093542d+00
      e(12)=0.03021553d+00
      c(10)=-0.0999744d+00
      c(11)=0.5545757d+00
      c(12)=0.5193358d+00
      e(13)=198.85110d+00
      e(14)=45.443908d+00
      e(15)=12.984531d+00
      c(13)=0.0927924d+00
      c(14)=0.4476071d+00
      c(15)=0.6167211d+00
      e(16)=4.8835867d+00
      e(17)=1.8407375d+00
      e(18)=0.65687534d+00
      c(16)=0.2550456d+00
      c(17)=0.5708333d+00
      c(18)=0.2998891d+00
      e(19)=12.634196d+00
      e(20)=3.0847767d+00
      e(21)=0.69473162d+00
      c(19)=0.1704471d+00
      c(20)=0.5085159d+00
      c(21)=0.6080700d+00
      return
c
c     nickel (3-d)
c
  800 continue
      e( 1)=3849.2417d+00
      e( 2)=582.08183d+00
      e( 3)=127.35676d+00
      c( 1)=0.0632462d+00
      c( 2)=0.3751737d+00
      c( 3)=0.6828494d+00
      e( 4)=170.11363d+00
      e( 5)=15.852695d+00
      e( 6)=6.4342950d+00
      c( 4)=-0.1042343d+00
      c( 5)=0.7012328d+00
      c( 6)=0.3655163d+00
      e( 7)=12.889647d+00
      e( 8)=2.0461387d+00
      e( 9)=0.77923389d+00
      c( 7)=-0.2236128d+00
      c( 8)=0.7180189d+00
      c( 9)=0.4151786d+00
      e(10)=0.84339749d+00
      e(11)=0.08247463d+00
      e(12)=0.03073497d+00
      c(10)=-0.0938590d+00
      c(11)=0.5512875d+00
      c(12)=0.5211700d+00
      e(13)=216.01321d+00
      e(14)=49.433819d+00
      e(15)=14.166185d+00
      c(13)=0.0925038d+00
      c(14)=0.4477120d+00
      c(15)=0.6162927d+00
      e(16)=5.3861988d+00
      e(17)=2.0323051d+00
      e(18)=0.72556764d+00
      c(16)=0.2552417d+00
      c(17)=0.5684974d+00
      c(18)=0.3023861d+00
      e(19)=14.087823d+00
      e(20)=3.4535186d+00
      e(21)=0.77782386d+00
      c(19)=0.1719093d+00
      c(20)=0.5102939d+00
      c(21)=0.6052430d+00
      return
c
c     copper (2-s)
c
  900 continue
      e(1)=4137.5632d+00
      e(2)=625.73294d+00
      e(3)=136.95969d+00
      c(1)=0.0631365d+00
      c(2)=0.3747929d+00
      c(3)=0.6832095d+00
      e(4)=183.1051d+00
      e(5)=17.122418d+00
      e(6)=6.9491752d+00
      c(4)=-0.1046861d+00
      c(5)=0.7034840d+00
      c(6)=0.3634281d+00
      e(7)=13.999937d+00
      e(8)=2.2431225d+00
      e(9)=0.85310594d+00
      c(7)=-0.2239634d+00
      c(8)=0.7113033d+00
      c(9)=0.4227015d+00
      e(10)=0.83040542d+00
      e(11)=0.08540467d+00
      e(12)=0.03181982d+00
      c(10)=-0.0890494d+00
      c(11)=0.5338648d+00
      c(12)=0.5380713d+00
      e(13)=233.62353d+00
      e(14)=53.545569d+00
      e(15)=15.388041d+00
      c(13)=0.0923643d+00
      c(14)=0.4480449d+00
      c(15)=0.6155840d+00
      e(16)=5.8732783d+00
      e(17)=2.2084129d+00
      e(18)=0.78648338d+00
      c(16)=0.2593135d+00
      c(17)=0.5687950d+00
      c(18)=0.2987106d+00
      e(19)=15.586486d+00
      e(20)=3.8355780d+00
      e(21)=0.86504248d+00
      c(19)=0.1734535d+00
      c(20)=0.5121846d+00
      c(21)=0.6020146d+00
      return
c
c     zinc (1-s)
c
 1000 continue
      e( 1)=4432.2885d+00
      e( 2)=670.66012d+00
      e( 3)=146.90245d+00
      c( 1)=0.0630928d+00
      c( 2)=0.3745038d+00
      c( 3)=0.6834160d+00
      e( 4)=196.48911d+00
      e( 5)=18.441663d+00
      e( 6)=7.4929688d+00
      c( 4)=-0.1050863d+00
      c( 5)=0.7049994d+00
      c( 6)=0.3620077d+00
      e( 7)=15.180864d+00
      e( 8)=2.4467928d+00
      e( 9)=0.95293505d+00
      c( 7)=-0.2254217d+00
      c( 8)=0.7096746d+00
      c( 9)=0.4236169d+00
      e(10)=1.1210370d+00
      e(11)=0.1251838d+00
      e(12)=0.04406728d+00
      c(10)=-0.1166693d+00
      c(11)=0.5950219d+00
      c(12)=0.4923302d+00
      e(13)=252.47779d+00
      e(14)=57.934825d+00
      e(15)=16.689667d+00
      c(13)=0.0919491d+00
      c(14)=0.4476288d+00
      c(15)=0.6158034d+00
      e(16)=6.3410797d+00
      e(17)=2.3880503d+00
      e(18)=0.87307830d+00
      c(16)=0.2675452d+00
      c(17)=0.5701029d+00
      c(18)=0.2856214d+00
      e(19)=18.368202d+00
      e(20)=4.5913041d+00
      e(21)=1.0902026d+00
      c(19)=0.1672532d+00
      c(20)=0.5122790d+00
      c(21)=0.5948177d+00
      return
      end
      subroutine htmtwo(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
      ibr=nucz-38
      go to (100,200,300,400,500,600,700,800,900,1000), ibr
c
c     yttrium (4-f)
c
  100 continue
      e( 1)=7647.6932d+00
      e( 2)=1157.0020d+00
      e( 3)=253.75700d+00
      c( 1)=0.0618818d+00
      c( 2)=0.3701570d+00
      c( 3)=0.6878017d+00
      e( 4)=339.10364d+00
      e( 5)=32.403264d+00
      e( 6)=13.281887d+00
      c( 4)=-0.1075786d+00
      c( 5)=0.7181663d+00
      c( 6)=0.3490842d+00
      e( 7)=27.996257d+00
      e( 8)=4.6829212d+00
      e( 9)=2.0566163d+00
      c( 7)=-0.2457049d+00
      c( 8)=0.7878323d+00
      c( 9)=0.3488907d+00
      e(10)=3.4864304d+00
      e(11)=0.66744514d+00
      e(12)=0.27959697d+00
      c(10)=-0.3000605d+00
      c(11)=0.7858715d+00
      c(12)=0.3857544d+00
      e(13)=0.46867231d+00
      e(14)=0.04978055d+00
      e(15)=0.02040820d+00
      c(13)=-0.2438980d+00
      c(14)=0.7576396d+00
      c(15)=0.3475036d+00
      e(16)=464.62389d+00
      e(17)=107.54806d+00
      e(18)=31.650408d+00
      c(16)=0.0870002d+00
      c(17)=0.4370123d+00
      c(18)=0.6250957d+00
      e(19)=12.443905d+00
      e(20)=4.9934558d+00
      e(21)=2.1033868d+00
      c(19)=0.2788493d+00
      c(20)=0.5885743d+00
      c(21)=0.2292208d+00
      e(22)=1.0947920d+00
      e(23)=0.46072626d+00
      e(24)=0.18428969d+00
      c(22)=0.4060140d+00
      c(23)=0.5320897d+00
      c(24)=0.1575146d+00
      e(25)=50.348619d+00
      e(26)=13.497089d+00
      e(27)=3.9166237d+00
      c(25)=0.1367063d+00
      c(26)=0.5029592d+00
      c(27)=0.5788425d+00
      e(28)=1.4046174d+00
      e(29)=0.45202594d+00
      e(30)=0.12497283d+00
      c(28)=0.1778019d+00
      c(29)=0.4652198d+00
      c(30)=0.5916440d+00
      return
c
c     zirconium (5-f)
c
  200 continue
      e( 1)=8074.5440d+00
      e( 2)=1220.5300d+00
      e( 3)=267.54751d+00
      c( 1)=0.0616606d+00
      c( 2)=0.3696140d+00
      c( 3)=0.6884549d+00
      e( 4)=357.68109d+00
      e( 5)=34.142568d+00
      e( 6)=13.947566d+00
      c( 4)=-0.1077339d+00
      c( 5)=0.7216851d+00
      c( 6)=0.3456921d+00
      e( 7)=29.610522d+00
      e( 8)=4.9809391d+00
      e( 9)=2.1850378d+00
      c( 7)=-0.2483862d+00
      c( 8)=0.8005982d+00
      c( 9)=0.3376774d+00
      e(10)=3.8399904d+00
      e(11)=0.73978700d+00
      e(12)=0.30961916d+00
      c(10)=-0.3104305d+00
      c(11)=0.8022937d+00
      c(12)=0.3741656d+00
      e(13)=0.55502132d+00
      e(14)=0.05532887d+00
      e(15)=0.02137166d+00
      c(13)=-0.2288857d+00
      c(14)=0.7817074d+00
      c(15)=0.3195193d+00
      e(16)=491.62198d+00
      e(17)=113.93040d+00
      e(18)=33.605515d+00
      c(16)=0.0867409d+00
      c(17)=0.4362728d+00
      c(18)=0.6255724d+00
      e(19)=13.247292d+00
      e(20)=5.3389914d+00
      e(21)=2.2525019d+00
      c(19)=0.2803844d+00
      c(20)=0.5922702d+00
      c(21)=0.2225147d+00
      e(22)=1.2314752d+00
      e(23)=0.52978846d+00
      e(24)=0.21408758d+00
      c(22)=0.4058518d+00
      c(23)=0.5246687d+00
      c(24)=0.1633953d+00
      e(25)=54.716068d+00
      e(26)=14.743888d+00
      e(27)=4.3333204d+00
      c(25)=0.1347642d+00
      c(26)=0.5014995d+00
      c(27)=0.5786696d+00
      e(28)=1.6859057d+00
      e(29)=0.53743517d+00
      e(30)=0.15022104d+00
      c(28)=0.1925337d+00
      c(29)=0.5052917d+00
      c(30)=0.5385693d+00
      return
c
c     niobium (6-s)
c
  300 continue
      e(1)=8465.0874d+00
      e(2)=1281.1163d+00
      e(3)=281.21159d+00
      c(1)=0.0618144d+00
      c(2)=0.3698298d+00
      c(3)=0.6880474d+00
      e(4)=375.85025d+00
      e(5)=36.053370d+00
      e(6)=14.816642d+00
      c(4)=-0.1080131d+00
      c(5)=0.7199938d+00
      c(6)=0.3472657d+00
      e(7)=31.365510d+00
      e(8)=5.2718932d+00
      e(9)=2.3070864d+00
      c(7)=-0.2499774d+00
      c(8)=0.8156328d+00
      c(9)=0.3232937d+00
      e(10)=4.1677904d+00
      e(11)=0.8215584d+00
      e(12)=0.34324945d+00
      c(10)=-0.3163261d+00
      c(11)=0.8073422d+00
      c(12)=0.3740729d+00
      e(13)=0.58716698d+00
      e(14)=0.06616704d+00
      e(15)=0.02593466d+00
      c(13)=-0.2224673d+00
      c(14)=0.7048521d+00
      c(15)=0.4026700d+00
      e(16)=519.28752d+00
      e(17)=120.45899d+00
      e(18)=35.608105d+00
      c(16)=0.0865238d+00
      c(17)=0.4357148d+00
      c(18)=0.6258786d+00
      e(19)=13.976575d+00
      e(20)=5.6388612d+00
      e(21)=2.3712274d+00
      c(19)=0.2867436d+00
      c(20)=0.5975733d+00
      c(21)=0.2096306d+00
      e(22)=1.3455223d+00
      e(23)=0.56846749d+00
      e(24)=0.22696230d+00
      c(22)=0.4308749d+00
      c(23)=0.5215640d+00
      c(24)=0.1404843d+00
      e(25)=59.035804d+00
      e(26)=15.997120d+00
      e(27)=4.7574796d+00
      c(25)=0.1335596d+00
      c(26)=0.5006254d+00
      c(27)=0.5776791d+00
      e(28)=1.9792676d+00
      e(29)=0.62926194d+00
      e(30)=0.17752293d+00
      c(28)=0.2015902d+00
      c(29)=0.5271681d+00
      c(30)=0.5065745d+00
      return
c
c     molybdenum (7-s)
c
  400 continue
      e(1)=8898.9752d+00
      e(2)=1346.6930d+00
      e(3)=295.62876d+00
      c(1)=0.0617106d+00
      c(2)=0.3694639d+00
      c(3)=0.6884201d+00
      e(4)=395.08162d+00
      e(5)=37.964578d+00
      e(6)=15.608898d+00
      c(4)=-0.1082198d+00
      c(5)=0.7205080d+00
      c(6)=0.3468176d+00
      e(7)=33.085228d+00
      e(8)=5.6180351d+00
      e(9)=2.4825318d+00
      c(7)=-0.2526233d+00
      c(8)=0.8177164d+00
      c(9)=0.3225675d+00
      e(10)=4.5000086d+00
      e(11)=0.90593940d+00
      e(12)=0.37777344d+00
      c(10)=-0.3206422d+00
      c(11)=0.8092539d+00
      c(12)=0.3763588d+00
      e(13)=0.65419456d+00
      e(14)=0.07258102d+00
      e(15)=0.02797177d+00
      c(13)=-0.2116815d+00
      c(14)=0.6864946d+00
      c(15)=0.4187549d+00
      e(16)=548.86388d+00
      e(17)=127.39268d+00
      e(18)=37.712842d+00
      c(16)=0.0860665d+00
      c(17)=0.4346698d+00
      c(18)=0.6268737d+00
      e(19)=14.744588d+00
      e(20)=5.9600351d+00
      e(21)=2.4937534d+00
      c(19)=0.2917026d+00
      c(20)=0.6020321d+00
      c(21)=0.1991546d+00
      e(22)=1.4621554d+00
      e(23)=0.60827575d+00
      e(24)=0.23912471d+00
      c(22)=0.4518315d+00
      c(23)=0.5172613d+00
      c(24)=0.1232182d+00
      e(25)=63.792046d+00
      e(26)=17.356419d+00
      e(27)=5.2109516d+00
      c(25)=0.1316317d+00
      c(26)=0.4990747d+00
      c(27)=0.5782636d+00
      e(28)=2.2861428d+00
      e(29)=0.72931951d+00
      e(30)=0.2094156d+00
      c(28)=0.2072519d+00
      c(29)=0.5396099d+00
      c(30)=0.4846239d+00
      return
c
c     technetium (6-d)
c
  500 continue
      e(1)=9331.8660d+00
      e(2)=1412.7444d+00
      e(3)=310.28294d+00
      c(1)=0.0616947d+00
      c(2)=0.3693039d+00
      c(3)=0.6885214d+00
      e(4)=414.56635d+00
      e(5)=39.944899d+00
      e(6)=16.451697d+00
      c(4)=-0.1084438d+00
      c(5)=0.7206125d+00
      c(6)=0.3467395d+00
      e(7)=34.887843d+00
      e(8)=5.9470507d+00
      e(9)=2.6240416d+00
      c(7)=-0.2545651d+00
      c(8)=0.8276066d+00
      c(9)=0.3138042d+00
      e(10)=4.8541288d+00
      e(11)=0.98792524d+00
      e(12)=0.41093096d+00
      c(10)=-0.3241571d+00
      c(11)=0.8151655d+00
      c(12)=0.3733329d+00
      e(13)=0.71023291d+00
      e(14)=0.07580316d+00
      e(15)=0.02868756d+00
      c(13)=-0.1967391d+00
      c(14)=0.6579450d+00
      c(15)=0.4429358d+00
      e(16)=580.94360d+00
      e(17)=135.36036d+00
      e(18)=39.889938d+00
      c(16)=0.0845326d+00
      c(17)=0.4321914d+00
      c(18)=0.6304900d+00
      e(19)=87.347195d+00
      e(20)=13.287886d+00
      e(21)=4.9090562d+00
      c(19)=-0.0295523d+00
      c(20)=0.4201419d+00
      c(21)=0.6547726d+00
      e(22)=1.9258524d+00
      e(23)=0.84385557d+00
      e(24)=0.33618541d+00
      c(22)=0.3011680d+00
      c(23)=0.5688683d+00
      c(24)=0.2252629d+00
      e(25)=69.016818d+00
      e(26)=18.830107d+00
      e(27)=5.6925288d+00
      c(25)=0.1290297d+00
      c(26)=0.4966878d+00
      c(27)=0.5805887d+00
      e(28)=2.6360931d+00
      e(29)=0.84019534d+00
      e(30)=0.23965135d+00
      c(28)=0.2078727d+00
      c(29)=0.5439701d+00
      c(30)=0.4804963d+00
      return
c
c     ruthenium (5-f)
c
  600 continue
      e(1)=9786.1429d+00
      e(2)=1481.4714d+00
      e(3)=325.40260d+00
      c(1)=0.0616049d+00
      c(2)=0.3689840d+00
      c(3)=0.6888455d+00
      e(4)=434.79530d+00
      e(5)=41.918783d+00
      e(6)=17.257795d+00
      c(4)=-0.1086056d+00
      c(5)=0.7222292d+00
      c(6)=0.3451630d+00
      e(7)=36.735183d+00
      e(8)=6.2780444d+00
      e(9)=2.7519759d+00
      c(7)=-0.2564513d+00
      c(8)=0.8399185d+00
      c(9)=0.3027171d+00
      e(10)=5.2173878d+00
      e(11)=1.0759512d+00
      e(12)=0.44710855d+00
      c(10)=-0.3265667d+00
      c(11)=0.8168014d+00
      c(12)=0.3742870d+00
      e(13)=0.79343625d+00
      e(14)=0.07742209d+00
      e(15)=0.02934842d+00
      c(13)=-0.1864361d+00
      c(14)=0.6339363d+00
      c(15)=0.4599697d+00
      e(16)=613.23511d+00
      e(17)=142.89492d+00
      e(18)=42.136712d+00
      c(16)=0.0839306d+00
      c(17)=0.4310037d+00
      c(18)=0.6318141d+00
      e(19)=91.243693d+00
      e(20)=14.138100d+00
      e(21)=5.2568481d+00
      c(19)=-0.0298490d+00
      c(20)=0.4184301d+00
      c(21)=0.6558081d+00
      e(22)=2.0542379d+00
      e(23)=0.88489893d+00
      e(24)=0.35005918d+00
      c(22)=0.3259420d+00
      c(23)=0.5693077d+00
      c(24)=0.2001975d+00
      e(25)=74.088526d+00
      e(26)=20.283721d+00
      e(27)=6.1778559d+00
      c(25)=0.1274352d+00
      c(26)=0.4953767d+00
      c(27)=0.5811668d+00
      e(28)=2.8983401d+00
      e(29)=0.92735408d+00
      e(30)=0.26618991d+00
      c(28)=0.2189550d+00
      c(29)=0.5502884d+00
      c(30)=0.4631699d+00
      return
c
c     rhodium (4-f)
c
  700 continue
      e(1)=10225.249d+00
      e(2)=1549.1390d+00
      e(3)=340.56851d+00
      c(1)=0.0616846d+00
      c(2)=0.3690559d+00
      c(3)=0.6886519d+00
      e(4)=455.17966d+00
      e(5)=43.986950d+00
      e(6)=18.132961d+00
      c(4)=-0.1088006d+00
      c(5)=0.7224507d+00
      c(6)=0.3449743d+00
      e(7)=38.751877d+00
      e(8)=6.6037791d+00
      e(9)=2.8952705d+00
      c(7)=-0.2573715d+00
      c(8)=0.8517085d+00
      c(9)=0.2909376d+00
      e(10)=5.6152907d+00
      e(11)=1.1613908d+00
      e(12)=0.48020132d+00
      c(10)=-0.3266525d+00
      c(11)=0.8214561d+00
      c(12)=0.3704174d+00
      e(13)=0.85451427d+00
      e(14)=0.07947300d+00
      e(15)=0.03009295d+00
      c(13)=-0.1758253d+00
      c(14)=0.6050953d+00
      c(15)=0.4842742d+00
      e(16)=644.23213d+00
      e(17)=150.28013d+00
      e(18)=44.384507d+00
      c(16)=0.0837350d+00
      c(17)=0.4304508d+00
      c(18)=0.6321828d+00
      e(19)=95.833303d+00
      e(20)=14.963302d+00
      e(21)=5.6013156d+00
      c(19)=-0.0300917d+00
      c(20)=0.4189458d+00
      c(21)=0.6546474d+00
      e(22)=2.1513076d+00
      e(23)=0.90932095d+00
      e(24)=0.35848243d+00
      c(22)=0.3612858d+00
      c(23)=0.5622628d+00
      c(24)=0.1716521d+00
      e(25)=79.284881d+00
      e(26)=21.777632d+00
      e(27)=6.6786285d+00
      c(25)=0.1260477d+00
      c(26)=0.4942743d+00
      c(27)=0.5815451d+00
      e(28)=3.2225210d+00
      e(29)=1.0349052d+00
      e(30)=0.29821448d+00
      c(28)=0.2230469d+00
      c(29)=0.5533923d+00
      c(30)=0.4551932d+00
      return
c
c     palladium (3-d)
c
  800 continue
      e(1)=10732.268d+00
      e(2)=1624.2557d+00
      e(3)=356.77443d+00
      c(1)=0.0614091d+00
      c(2)=0.3683372d+00
      c(3)=0.6895216d+00
      e(4)=476.45981d+00
      e(5)=46.089079d+00
      e(6)=19.011457d+00
      c(4)=-0.1089646d+00
      c(5)=0.7229184d+00
      c(6)=0.3445285d+00
      e(7)=40.501816d+00
      e(8)=7.0317180d+00
      e(9)=3.1102881d+00
      c(7)=-0.2610873d+00
      c(8)=0.8473020d+00
      c(9)=0.2980309d+00
      e(10)=5.9535538d+00
      e(11)=1.2737777d+00
      e(12)=0.52932056d+00
      c(10)=-0.3315604d+00
      c(11)=0.8065730d+00
      c(12)=0.3906202d+00
      e(13)=0.91680596d+00
      e(14)=0.08014658d+00
      e(15)=0.03033468d+00
      c(13)=-0.1665563d+00
      c(14)=0.5960441d+00
      c(15)=0.4889716d+00
      e(16)=676.16554d+00
      e(17)=157.85238d+00
      e(18)=46.699134d+00
      c(16)=0.0835399d+00
      c(17)=0.4299341d+00
      c(18)=0.6325113d+00
      e(19)=101.25606d+00
      e(20)=15.760330d+00
      e(21)=5.9416372d+00
      c(19)=-0.0303593d+00
      c(20)=0.4214283d+00
      c(21)=0.6515363d+00
      e(22)=2.3082437d+00
      e(23)=0.97075110d+00
      e(24)=0.38294812d+00
      c(22)=0.3727721d+00
      c(23)=0.5581097d+00
      c(24)=0.1642278d+00
      e(25)=84.456904d+00
      e(26)=23.277977d+00
      e(27)=7.1869656d+00
      c(25)=0.1251206d+00
      c(26)=0.4936143d+00
      c(27)=0.5812637d+00
      e(28)=3.4886548d+00
      e(29)=1.1266213d+00
      e(30)=0.32670438d+00
      c(28)=0.2326503d+00
      c(29)=0.5566155d+00
      c(30)=0.4419046d+00
      return
c
c     silver (2-s)
c
  900 continue
      e(1)=11195.393d+00
      e(2)=1695.6114d+00
      e(3)=372.71564d+00
      c(1)=0.0614652d+00
      c(2)=0.3683407d+00
      c(3)=0.6894170d+00
      e(4)=498.08789d+00
      e(5)=48.189273d+00
      e(6)=19.879233d+00
      c(4)=-0.1090984d+00
      c(5)=0.7248356d+00
      c(6)=0.3425905d+00
      e(7)=42.617304d+00
      e(8)=7.3549843d+00
      e(9)=3.2193663d+00
      c(7)=-0.2616127d+00
      c(8)=0.8638930d+00
      c(9)=0.2814949d+00
      e(10)=6.4471843d+00
      e(11)=1.3406340d+00
      e(12)=0.55082491d+00
      c(10)=-0.3264484d+00
      c(11)=0.8244194d+00
      c(12)=0.3684570d+00
      e(13)=0.896024d+00
      e(14)=0.07873945d+00
      e(15)=0.02942569d+00
      c(13)=-0.1464338d+00
      c(14)=0.6078123d+00
      c(15)=0.474578d+00
      e(16)=709.79348d+00
      e(17)=165.75065d+00
      e(18)=49.079374d+00
      c(16)=0.0832016d+00
      c(17)=0.4292951d+00
      c(18)=0.6331456d+00
      e(19)=105.90189d+00
      e(20)=16.746072d+00
      e(21)=6.3433500d+00
      c(19)=-0.0305549d+00
      c(20)=0.4170359d+00
      c(21)=0.6552444d+00
      e(22)=2.511517d+00
      e(23)=1.0560803d+00
      e(24)=0.41595731d+00
      c(22)=0.3678284d+00
      c(23)=0.5598336d+00
      c(24)=0.1679555d+00
      e(25)=90.022978d+00
      e(26)=24.875706d+00
      e(27)=7.7224910d+00
      c(25)=0.1237838d+00
      c(26)=0.4925467d+00
      c(27)=0.5818767d+00
      e(28)=3.7854315d+00
      e(29)=1.2278562d+00
      e(30)=0.35902458d+00
      c(28)=0.2395459d+00
      c(29)=0.5591501d+00
      c(30)=0.4311417d+00
      return
c
c     cadmium (1-s)
c
 1000 continue
      e(1)=11686.086d+00
      e(2)=1770.1114d+00
      e(3)=389.20899d+00
      c(1)=0.0614265d+00
      c(2)=0.3681567d+00
      c(3)=0.6895722d+00
      e(4)=519.91106d+00
      e(5)=50.397993d+00
      e(6)=20.834269d+00
      c(4)=-0.1092676d+00
      c(5)=0.7250389d+00
      c(6)=0.3423508d+00
      e(7)=44.652532d+00
      e(8)=7.7506717d+00
      e(9)=3.3818076d+00
      c(7)=-0.2634441d+00
      c(8)=0.8694815d+00
      c(9)=0.2772414d+00
      e(10)=6.7897003d+00
      e(11)=1.4664135d+00
      e(12)=0.61806063d+00
      c(10)=-0.3343146d+00
      c(11)=0.8177517d+00
      c(12)=0.3812846d+00
      e(13)=0.91773482d+00
      e(14)=0.11202343d+00
      e(15)=0.04082857d+00
      c(13)=-0.1509803d+00
      c(14)=0.6086784d+00
      c(15)=0.4891503d+00
      e(16)=745.10479d+00
      e(17)=174.04696d+00
      e(18)=51.552438d+00
      c(16)=0.0827347d+00
      c(17)=0.4283685d+00
      c(18)=0.6341921d+00
      e(19)=106.27841d+00
      e(20)=17.924992d+00
      e(21)=6.7904099d+00
      c(19)=-0.0305639d+00
      c(20)=0.4096720d+00
      c(21)=0.6626091d+00
      e(22)=2.5792466d+00
      e(23)=1.0573897d+00
      e(24)=0.40668464d+00
      c(22)=0.415516d+00
      c(23)=0.5559935d+00
      c(24)=0.1202681d+00
      e(25)=95.472743d+00
      e(26)=26.481959d+00
      e(27)=8.2828858d+00
      c(25)=0.1230828d+00
      c(26)=0.4916768d+00
      c(27)=0.5815408d+00
      e(28)=4.0821413d+00
      e(29)=1.3572792d+00
      e(30)=0.42083076d+00
      c(28)=0.2516043d+00
      c(29)=0.5652168d+00
      c(30)=0.3995299d+00
      return
      end
      subroutine htmthr(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
c     huzinaga third row transition series...5c(n),6c(2) states.
c
      ihop = nucz-70
      go to (100,200,300,400,500,600,700,800,900,1000), ihop
c
c     lutetium 2-d
c
  100 continue
      e(1)= 2.6213072d+04
      e(2)= 3.9645756d+03
      e(3)= 8.7202321d+02
      c(1)= 5.99916d-02
      c(2)= 3.631627d-01
      c(3)= 6.947236d-01
      e(4)= 1.1596313d+03
      e(5)= 1.1482931d+02
      e(6)= 4.8944774d+01
      c(4)=-1.120702d-01
      c(5)= 7.324106d-01
      c(6)= 3.342924d-01
      e(7)= 1.7955607d+03
      e(8)= 4.2201312d+02
      e(9)= 1.2679229d+02
      c(7)= 7.62034d-02
      c(8)= 4.127937d-01
      c(9)= 6.497893d-01
      e(10)= 1.0517731d+02
      e(11)= 2.0488594d+01
      e(12)= 9.4407118d+00
      c(10)=-2.944235d-01
      c(11)= 8.862065d-01
      c(12)= 2.818423d-01
      e(13)= 2.4023172d+02
      e(14)= 4.7779063d+01
      e(15)= 1.9398857d+01
      c(13)=-3.63245d-02
      c(14)= 3.912159d-01
      c(15)= 6.751724d-01
      e(16)= 2.7718839d+02
      e(17)= 7.9463212d+01
      e(18)= 2.6610039d+01
      c(16)= 1.055301d-01
      c(17)= 4.738686d-01
      c(18)= 5.941071d-01
      e(19)= 2.029241d+01
      e(20)= 4.6894803d+00
      e(21)= 2.5387204d+00
      c(19)=-3.484181d-01
      c(20)= 7.636206d-01
      c(21)= 4.37517d-01
      e(22)= 9.5967927d+00
      e(23)= 4.4367903d+00
      e(24)= 2.0246813d+00
      c(22)= 3.360259d-01
      c(23)= 5.591485d-01
      c(24)= 1.807147d-01
      e(25)= 1.5600726d+01
      e(26)= 5.9803128d+00
      e(27)= 2.3344255d+00
      c(25)= 2.805518d-01
      c(26)= 5.66057d-01
      c(27)= 3.127337d-01
      e(28)= 3.1797464d+01
      e(29)= 9.14725d+00
      e(30)= 2.4942441d+00
      c(28)= 1.889543d-01
      c(29)= 5.299306d-01
      c(30)= 5.769141d-01
      e(31)= 5.0301304d+00
      e(32)= 8.8480792d-01
      e(33)= 4.4784708d-01
      c(31)=-4.28184d-02
      c(32)= 2.562323d-01
      c(33)= 7.786912d-01
      e(34)= 2.1297889d+00
      e(35)= 8.9173211d-01
      e(36)= 3.2974496d-01
      c(34)= 4.09732d-02
      c(35)= 5.605793d-01
      c(36)= 4.812806d-01
      e(37)= 2.0469529d+00
      e(38)= 6.8141564d-01
      e(39)= 2.0136947d-01
      c(37)= 8.70835d-02
      c(38)= 3.990061d-01
      c(39)= 6.900484d-01
      e(40)= 1.0345315d+00
      e(41)= 8.244872d-02
      e(42)= 3.059794d-02
      c(40)= 8.1952d-02
      c(41)= 3.276736d-01
      c(42)= 6.843155d-01
      return
c
c     hafnium 3-f
c
  200 continue
      e(1)= 2.6969327d+04
      e(2)= 4.0796473d+03
      e(3)= 8.9744182d+02
      c(1)= 5.99635d-02
      c(2)= 3.630082d-01
      c(3)= 6.948626d-01
      e(4)= 1.1911456d+03
      e(5)= 1.178739d+02
      e(6)= 4.9501542d+01
      c(4)=-1.119877d-01
      c(5)= 7.381808d-01
      c(6)= 3.292139d-01
      e(7)= 1.8523764d+03
      e(8)= 4.3542673d+02
      e(9)= 1.3087367d+02
      c(7)= 7.59814d-02
      c(8)= 4.122427d-01
      c(9)= 6.503501d-01
      e(10)= 1.0779325d+02
      e(11)= 2.1253125d+01
      e(12)= 9.660219d+00
      c(10)=-2.97236d-01
      c(11)= 8.867886d-01
      c(12)= 2.843174d-01
      e(13)= 2.4845071d+02
      e(14)= 4.9377736d+01
      e(15)= 2.0077755d+01
      c(13)=-3.65325d-02
      c(14)= 3.90915d-01
      c(15)= 6.753456d-01
      e(16)= 2.8685428d+02
      e(17)= 8.2345093d+01
      e(18)= 2.7632269d+01
      c(16)= 1.052524d-01
      c(17)= 4.734033d-01
      c(18)= 5.942747d-01
      e(19)= 2.0900824d+01
      e(20)= 4.8004147d+00
      e(21)= 2.5851469d+00
      c(19)=-3.510341d-01
      c(20)= 7.977843d-01
      c(21)= 4.052346d-01
      e(22)= 9.8265412d+00
      e(23)= 4.5507968d+00
      e(24)= 2.1136634d+00
      c(22)= 3.475344d-01
      c(23)= 5.527717d-01
      c(24)= 1.742289d-01
      e(25)= 1.6251487d+01
      e(26)= 6.2625999d+00
      e(27)= 2.4873801d+00
      c(25)= 2.809288d-01
      c(26)= 5.630773d-01
      c(27)= 3.126199d-01
      e(28)= 3.4279107d+01
      e(29)= 9.9276463d+00
      e(30)= 2.7558598d+00
      c(28)= 1.855064d-01
      c(29)= 5.305358d-01
      c(30)= 5.738195d-01
      e(31)= 4.8165887d+00
      e(32)= 9.190986d-01
      e(33)= 4.6756457d-01
      c(31)=-1.1025d-01
      c(32)= 3.161871d-01
      c(33)= 7.519785d-01
      e(34)= 2.0644787d+00
      e(35)= 9.3120381d-01
      e(36)= 3.426192d-01
      c(34)= 5.02644d-02
      c(35)= 5.912507d-01
      c(36)= 4.410394d-01
      e(37)= 2.0540301d+00
      e(38)= 6.9039011d-01
      e(39)= 2.1235875d-01
      c(37)= 1.125023d-01
      c(38)= 4.581519d-01
      c(39)= 6.142745d-01
      e(40)= 1.0028075d+00
      e(41)= 8.478757d-02
      e(42)= 3.146592d-02
      c(40)= 9.64616d-02
      c(41)= 4.044906d-01
      c(42)= 6.034115d-01
      return
c
c     tantalum 4-f
c
  300 continue
      e(1)= 2.7707063d+04
      e(2)= 4.1922778d+03
      e(3)= 9.2272356d+02
      c(1)= 6.0013d-02
      c(2)= 3.630896d-01
      c(3)= 6.947137d-01
      e(4)= 1.2259011d+03
      e(5)= 1.2142987d+02
      e(6)= 5.0889951d+01
      c(4)=-1.120743d-01
      c(5)= 7.373676d-01
      c(6)= 3.302371d-01
      e(7)= 1.9081438d+03
      e(8)= 4.4862494d+02
      e(9)= 1.3493024d+02
      c(7)= 7.58846d-02
      c(8)= 4.120017d-01
      c(9)= 6.505247d-01
      e(10)= 1.1155105d+02
      e(11)= 2.174652d+01
      e(12)= 1.0030729d+01
      c(10)=-2.968041d-01
      c(11)= 8.976579d-01
      c(12)= 2.714252d-01
      e(13)= 2.5483949d+02
      e(14)= 5.1093731d+01
      e(15)= 2.0793577d+01
      c(13)=-3.66507d-02
      c(14)= 3.898915d-01
      c(15)= 6.763348d-01
      e(16)= 2.9767187d+02
      e(17)= 8.5493895d+01
      e(18)= 2.872252d+01
      c(16)= 1.044578d-01
      c(17)= 4.722716d-01
      c(18)= 5.955515d-01
      e(19)= 2.1501135d+01
      e(20)= 4.989021d+00
      e(21)= 2.6049632d+00
      c(19)=-3.545572d-01
      c(20)= 8.174897d-01
      c(21)= 3.900778d-01
      e(22)= 1.0057908d+01
      e(23)= 4.5799849d+00
      e(24)= 2.1300786d+00
      c(22)= 3.645967d-01
      c(23)= 5.571217d-01
      c(24)= 1.518756d-01
      e(25)= 1.6876653d+01
      e(26)= 6.514203d+00
      e(27)= 2.6052822d+00
      c(25)= 2.825487d-01
      c(26)= 5.651022d-01
      c(27)= 3.072128d-01
      e(28)= 3.669239d+01
      e(29)= 1.0694725d+01
      e(30)= 3.0231809d+00
      c(28)= 1.829896d-01
      c(29)= 5.308347d-01
      c(30)= 5.704735d-01
      e(31)= 4.7611306d+00
      e(32)= 9.7852208d-01
      e(33)= 4.8069074d-01
      c(31)=-1.740521d-01
      c(32)= 3.866417d-01
      c(33)= 7.160257d-01
      e(34)= 1.9903475d+00
      e(35)= 8.8952373d-01
      e(36)= 3.4184436d-01
      c(34)= 1.092912d-01
      c(35)= 6.171303d-01
      c(36)= 3.579635d-01
      e(37)= 2.0602332d+00
      e(38)= 7.0272676d-01
      e(39)= 2.1834769d-01
      c(37)= 1.391422d-01
      c(38)= 5.065767d-01
      c(39)= 5.467383d-01
      e(40)= 9.7168133d-01
      e(41)= 8.715914d-02
      e(42)= 3.234604d-02
      c(40)= 1.133863d-01
      c(41)= 4.626334d-01
      c(42)= 5.371199d-01
      return
c
c     tungsten 5-d
c
  400 continue
      e(1)= 2.848347d+04
      e(2)= 4.3106522d+03
      e(3)= 9.4880189d+02
      c(1)= 5.99847d-02
      c(2)= 3.629424d-01
      c(3)= 6.948544d-01
      e(4)= 1.2601431d+03
      e(5)= 1.2499189d+02
      e(6)= 5.3264926d+01
      c(4)=-1.122722d-01
      c(5)= 7.356067d-01
      c(6)= 3.311238d-01
      e(7)= 1.9659337d+03
      e(8)= 4.6228149d+02
      e(9)= 1.3907004d+02
      c(7)= 7.57123d-02
      c(8)= 4.116367d-01
      c(9)= 6.509095d-01
      e(10)= 1.1485307d+02
      e(11)= 2.2490909d+01
      e(12)= 1.0205372d+01
      c(10)=-2.968151d-01
      c(11)= 8.959454d-01
      c(12)= 2.743324d-01
      e(13)= 2.603359d+02
      e(14)= 5.2992087d+01
      e(15)= 2.1564622d+01
      c(13)=-3.68773d-02
      c(14)= 3.872042d-01
      c(15)= 6.790686d-01
      e(16)= 3.0802405d+02
      e(17)= 8.8529624d+01
      e(18)= 2.9802451d+01
      c(16)= 1.040791d-01
      c(17)= 4.717159d-01
      c(18)= 5.958943d-01
      e(19)= 2.1986566d+01
      e(20)= 5.1545586d+00
      e(21)= 2.6238159d+00
      c(19)=-3.6002d-01
      c(20)= 8.471585d-01
      c(21)= 3.660872d-01
      e(22)= 1.0448376d+01
      e(23)= 4.7961135d+00
      e(24)= 2.2000558d+00
      c(22)= 3.599774d-01
      c(23)= 5.60685d-01
      c(24)= 1.526327d-01
      e(25)= 1.7714492d+01
      e(26)= 6.8804872d+00
      e(27)= 2.8091248d+00
      c(25)= 2.790386d-01
      c(26)= 5.599718d-01
      c(27)= 3.12865d-01
      e(28)= 3.9064396d+01
      e(29)= 1.1456039d+01
      e(30)= 3.2885779d+00
      c(28)= 1.810851d-01
      c(29)= 5.314735d-01
      c(30)= 5.670319d-01
      e(31)= 4.5567679d+00
      e(32)= 9.9070497d-01
      e(33)= 4.9395055d-01
      c(31)=-3.082678d-01
      c(32)= 4.976985d-01
      c(33)= 6.665086d-01
      e(34)= 1.9472206d+00
      e(35)= 8.6995199d-01
      e(36)= 3.2961937d-01
      c(34)= 1.570695d-01
      c(35)= 6.436159d-01
      c(36)= 2.852549d-01
      e(37)= 2.076216d+00
      e(38)= 7.3677302d-01
      e(39)= 2.2781597d-01
      c(37)= 1.61212d-01
      c(38)= 5.329476d-01
      c(39)= 5.016985d-01
      e(40)= 9.4116795d-01
      e(41)= 8.956342d-02
      e(42)= 3.323831d-02
      c(40)= 1.050188d-01
      c(41)= 4.969138d-01
      c(42)= 5.046213d-01
      return
c
c     rhenium 6-s
c
  500 continue
      e(1)= 2.9263742d+04
      e(2)= 4.4298015d+03
      e(3)= 9.7511876d+02
      c(1)= 5.99726d-02
      c(2)= 3.628568d-01
      c(3)= 6.949255d-01
      e(4)= 1.2932196d+03
      e(5)= 1.2860003d+02
      e(6)= 5.4807887d+01
      c(4)=-1.123781d-01
      c(5)= 7.355024d-01
      c(6)= 3.31334d-01
      e(7)= 2.0244843d+03
      e(8)= 4.7614657d+02
      e(9)= 1.4333076d+02
      c(7)= 7.55597d-02
      c(8)= 4.112017d-01
      c(9)= 6.512966d-01
      e(10)= 1.1845628d+02
      e(11)= 2.3141641d+01
      e(12)= 1.0379264d+01
      c(10)=-2.962528d-01
      c(11)= 9.005201d-01
      c(12)= 2.697741d-01
      e(13)= 2.7029118d+02
      e(14)= 5.4486952d+01
      e(15)= 2.2248123d+01
      c(13)=-3.70188d-02
      c(14)= 3.885717d-01
      c(15)= 6.773764d-01
      e(16)= 3.1855059d+02
      e(17)= 9.1648758d+01
      e(18)= 3.0906581d+01
      c(16)= 1.036622d-01
      c(17)= 4.710062d-01
      c(18)= 5.964425d-01
      e(19)= 2.2971336d+01
      e(20)= 5.3257089d+00
      e(21)= 2.7498632d+00
      c(19)=-3.550665d-01
      c(20)= 8.4405d-01
      c(21)= 3.649955d-01
      e(22)= 1.0747778d+01
      e(23)= 4.8975554d+00
      e(24)= 2.2377669d+00
      c(22)= 3.693555d-01
      c(23)= 5.63724d-01
      c(24)= 1.393835d-01
      e(25)= 1.8298578d+01
      e(26)= 7.1624476d+00
      e(27)= 2.9489065d+00
      c(25)= 2.810588d-01
      c(26)= 5.588154d-01
      c(27)= 3.099727d-01
      e(28)= 4.1648317d+01
      e(29)= 1.2273287d+01
      e(30)= 3.5704082d+00
      c(28)= 1.783847d-01
      c(29)= 5.314959d-01
      c(30)= 5.652308d-01
      e(31)= 4.7595836d+00
      e(32)= 1.0178146d+00
      e(33)= 5.0739077d-01
      c(31)=-3.10233d-01
      c(32)= 5.625827d-01
      c(33)= 6.054301d-01
      e(34)= 1.8856196d+00
      e(35)= 8.8924388d-01
      e(36)= 3.4366701d-01
      c(34)= 1.999751d-01
      c(35)= 6.225831d-01
      c(36)= 2.615926d-01
      e(37)= 2.0604954d+00
      e(38)= 7.4967705d-01
      e(39)= 2.424829d-01
      c(37)= 1.961924d-01
      c(38)= 5.459967d-01
      c(39)= 4.501411d-01
      e(40)= 9.1127985d-01
      e(41)= 9.200041d-02
      e(42)= 3.414271d-02
      c(40)= 1.036582d-01
      c(41)= 5.125143d-01
      c(42)= 4.878821d-01
      return
c
c     osmium 5-d
c
  600 continue
      e(1)= 3.0105978d+04
      e(2)= 4.5546731d+03
      e(3)= 1.0023444d+03
      c(1)= 5.98649d-02
      c(2)= 3.625769d-01
      c(3)= 6.952575d-01
      e(4)= 1.3299061d+03
      e(5)= 1.3189507d+02
      e(6)= 5.5352028d+01
      c(4)=-1.122784d-01
      c(5)= 7.397582d-01
      c(6)= 3.277809d-01
      e(7)= 2.0883111d+03
      e(8)= 4.9108568d+02
      e(9)= 1.4776879d+02
      c(7)= 7.51719d-02
      c(8)= 4.103016d-01
      c(9)= 6.523726d-01
      e(10)= 1.2184431d+02
      e(11)= 2.3821993d+01
      e(12)= 1.0499604d+01
      c(10)=-2.966599d-01
      c(11)= 9.060556d-01
      c(12)= 2.652479d-01
      e(13)= 2.7266734d+02
      e(14)= 5.6742347d+01
      e(15)= 2.3097886d+01
      c(13)=-3.72868d-02
      c(14)= 3.835103d-01
      c(15)= 6.82828d-01
      e(16)= 3.300562d+02
      e(17)= 9.5019652d+01
      e(18)= 3.208058d+01
      c(16)= 1.028721d-01
      c(17)= 4.696901d-01
      c(18)= 5.978924d-01
      e(19)= 2.3220256d+01
      e(20)= 5.6648299d+00
      e(21)= 2.8808691d+00
      c(19)=-3.683901d-01
      c(20)= 8.288351d-01
      c(21)= 3.928901d-01
      e(22)= 1.0934214d+01
      e(23)= 4.965104d+00
      e(24)= 2.1875843d+00
      c(22)= 3.834679d-01
      c(23)= 5.658388d-01
      c(24)= 1.230106d-01
      e(25)= 1.9043432d+01
      e(26)= 7.4801362d+00
      e(27)= 3.1008998d+00
      c(25)= 2.804308d-01
      c(26)= 5.593308d-01
      c(27)= 3.083838d-01
      e(28)= 4.4179317d+01
      e(29)= 1.3087026d+01
      e(30)= 3.8529209d+00
      c(28)= 1.763139d-01
      c(29)= 5.316512d-01
      c(30)= 5.632001d-01
      e(31)= 4.6824502d+00
      e(32)= 1.077075d+00
      e(33)= 5.2101139d-01
      c(31)=-3.721096d-01
      c(32)= 6.39539d-01
      c(33)= 5.653353d-01
      e(34)= 1.8723995d+00
      e(35)= 8.9080145d-01
      e(36)= 3.5473666d-01
      c(34)= 2.522785d-01
      c(35)= 6.021502d-01
      c(36)= 2.279055d-01
      e(37)= 2.0861668d+00
      e(38)= 7.6997534d-01
      e(39)= 2.5281533d-01
      c(37)= 2.275213d-01
      c(38)= 5.470991d-01
      c(39)= 4.189813d-01
      e(40)= 8.6522667d-01
      e(41)= 9.542435d-02
      e(42)= 3.505925d-02
      c(40)=-8.1871d-03
      c(41)= 5.441173d-01
      c(42)= 5.030192d-01
      return
c
c     iridium 4-f
c
  700 continue
      e(1)= 3.0907827d+04
      e(2)= 4.6769981d+03
      e(3)= 1.0294994d+03
      c(1)= 5.98596d-02
      c(2)= 3.624691d-01
      c(3)= 6.953337d-01
      e(4)= 1.3673566d+03
      e(5)= 1.3548412d+02
      e(6)= 5.7389281d+01
      c(4)=-1.123837d-01
      c(5)= 7.396342d-01
      c(6)= 3.273403d-01
      e(7)= 2.1451999d+03
      e(8)= 5.0478331d+02
      e(9)= 1.5203985d+02
      c(7)= 7.51929d-02
      c(8)= 4.102334d-01
      c(9)= 6.522788d-01
      e(10)= 1.2557265d+02
      e(11)= 2.4514298d+01
      e(12)= 1.0831887d+01
      c(10)=-2.964114d-01
      c(11)= 9.069686d-01
      c(12)= 2.639532d-01
      e(13)= 2.8416628d+02
      e(14)= 5.8149227d+01
      e(15)= 2.3772619d+01
      c(13)=-3.73807d-02
      c(14)= 3.862631d-01
      c(15)= 6.796648d-01
      e(16)= 3.4065395d+02
      e(17)= 9.8218233d+01
      e(18)= 3.3213551d+01
      c(16)= 1.025775d-01
      c(17)= 4.69076d-01
      c(18)= 5.982939d-01
      e(19)= 2.4082749d+01
      e(20)= 5.786534d+00
      e(21)= 2.9424679d+00
      c(19)=-3.668321d-01
      c(20)= 8.518878d-01
      c(21)= 3.679458d-01
      e(22)= 1.1374593d+01
      e(23)= 5.1149536d+00
      e(24)= 2.223298d+00
      c(22)= 3.847966d-01
      c(23)= 5.73833d-01
      c(24)= 1.132512d-01
      e(25)= 1.9801434d+01
      e(26)= 7.8062497d+00
      e(27)= 3.242046d+00
      c(25)= 2.795121d-01
      c(26)= 5.612916d-01
      c(27)= 3.060363d-01
      e(28)= 4.6700526d+01
      e(29)= 1.3898596d+01
      e(30)= 4.1392565d+00
      c(28)= 1.746611d-01
      c(29)= 5.319165d-01
      c(30)= 5.609526d-01
      e(31)= 4.9603751d+00
      e(32)= 1.100215d+00
      e(33)= 5.0934516d-01
      c(31)=-3.599398d-01
      c(32)= 7.202836d-01
      c(33)= 4.825747d-01
      e(34)= 1.8664367d+00
      e(35)= 9.1439782d-01
      e(36)= 3.7323659d-01
      c(34)= 2.954146d-01
      c(35)= 5.703553d-01
      c(36)= 2.147306d-01
      e(37)= 2.1418014d+00
      e(38)= 8.0244678d-01
      e(39)= 2.5966359d-01
      c(37)= 2.501547d-01
      c(38)= 5.534561d-01
      c(39)= 3.922531d-01
      e(40)= 8.6227728d-01
      e(41)= 9.795204d-02
      e(42)= 3.598793d-02
      c(40)=-1.661545d-01
      c(41)= 5.908108d-01
      c(42)= 5.070257d-01
      return
c
c     platinum 3-f
c
  800 continue
      e(1)= 3.1715151d+04
      e(2)= 4.800447d+03
      e(3)= 1.0568133d+03
      c(1)= 5.98597d-02
      c(2)= 3.624256d-01
      c(3)= 6.953573d-01
      e(4)= 1.402488d+03
      e(5)= 1.391957d+02
      e(6)= 5.8595143d+01
      c(4)=-1.124334d-01
      c(5)= 7.402836d-01
      c(6)= 3.271179d-01
      e(7)= 2.2068799d+03
      e(8)= 5.1935114d+02
      e(9)= 1.565111d+02
      c(7)= 7.50316d-02
      c(8)= 4.097815d-01
      c(9)= 6.527002d-01
      e(10)= 1.2903195d+02
      e(11)= 2.5401783d+01
      e(12)= 1.1170637d+01
      c(10)=-2.968866d-01
      c(11)= 8.998716d-01
      c(12)= 2.724105d-01
      e(13)= 2.9310449d+02
      e(14)= 5.9918114d+01
      e(15)= 2.4556392d+01
      c(13)=-3.75407d-02
      c(14)= 3.853595d-01
      c(15)= 6.802771d-01
      e(16)= 3.5185991d+02
      e(17)= 1.0148732d+02
      e(18)= 3.4372939d+01
      c(16)= 1.021866d-01
      c(17)= 4.686069d-01
      c(18)= 5.986499d-01
      e(19)= 2.4521199d+01
      e(20)= 6.1520604d+00
      e(21)= 3.0498908d+00
      c(19)=-3.76516d-01
      c(20)= 8.367552d-01
      c(21)= 3.937522d-01
      e(22)= 1.1628939d+01
      e(23)= 5.2412093d+00
      e(24)= 2.2477054d+00
      c(22)= 3.943199d-01
      c(23)= 5.710945d-01
      c(24)= 1.057293d-01
      e(25)= 2.0563355d+01
      e(26)= 8.1266279d+00
      e(27)= 3.4000049d+00
      c(25)= 2.796901d-01
      c(26)= 5.621668d-01
      c(27)= 3.03321d-01
      e(28)= 4.9223375d+01
      e(29)= 1.4713828d+01
      e(30)= 4.4259372d+00
      c(28)= 1.732586d-01
      c(29)= 5.323738d-01
      c(30)= 5.587198d-01
      e(31)= 5.0181161d+00
      e(32)= 1.1523707d+00
      e(33)= 5.0991299d-01
      c(31)=-3.948925d-01
      c(32)= 7.877904d-01
      c(33)= 4.384118d-01
      e(34)= 1.8605783d+00
      e(35)= 8.9302366d-01
      e(36)= 3.7742695d-01
      c(34)= 3.623975d-01
      c(35)= 5.386728d-01
      c(36)= 1.779594d-01
      e(37)= 2.242815d+00
      e(38)= 8.3045525d-01
      e(39)= 2.6917025d-01
      c(37)= 2.71483d-01
      c(38)= 5.621581d-01
      c(39)= 3.627509d-01
      e(40)= 8.5904808d-01
      e(41)= 9.951759d-02
      e(42)= 3.692875d-02
      c(40)=-2.480406d-01
      c(41)= 6.300501d-01
      c(42)= 4.869011d-01
      return
c
c     gold 2-d
c
  900 continue
      e(1)= 3.2538784d+04
      e(2)= 4.9255753d+03
      e(3)= 1.0846021d+03
      c(1)= 5.9856d-02
      c(2)= 3.623664d-01
      c(3)= 6.953909d-01
      e(4)= 1.4392918d+03
      e(5)= 1.4287377d+02
      e(6)= 6.0413806d+01
      c(4)=-1.125149d-01
      c(5)= 7.407681d-01
      c(6)= 3.263717d-01
      e(7)= 2.2683556d+03
      e(8)= 5.3384887d+02
      e(9)= 1.6096667d+02
      c(7)= 7.49274d-02
      c(8)= 4.095379d-01
      c(9)= 6.529009d-01
      e(10)= 1.3290783d+02
      e(11)= 2.5950372d+01
      e(12)= 1.1561528d+01
      c(10)=-2.972133d-01
      c(11)= 9.102609d-01
      c(12)= 2.606326d-01
      e(13)= 3.0198242d+02
      e(14)= 6.1409206d+01
      e(15)= 2.5215656d+01
      c(13)=-3.75287d-02
      c(14)= 3.888529d-01
      c(15)= 6.767782d-01
      e(16)= 3.6341101d+02
      e(17)= 1.0488558d+02
      e(18)= 3.5572d+01
      c(16)= 1.017057d-01
      c(17)= 4.678669d-01
      c(18)= 5.993448d-01
      e(19)= 2.5227768d+01
      e(20)= 6.3776449d+00
      e(21)= 3.1285944d+00
      c(19)=-3.796176d-01
      c(20)= 8.433751d-01
      c(21)= 3.902769d-01
      e(22)= 1.1983824d+01
      e(23)= 5.3684081d+00
      e(24)= 2.2381172d+00
      c(22)= 4.010601d-01
      c(23)= 5.743307d-01
      c(24)= 9.54635d-02
      e(25)= 2.1143756d+01
      e(26)= 8.3920687d+00
      e(27)= 3.5488617d+00
      c(25)= 2.840649d-01
      c(26)= 5.617205d-01
      c(27)= 2.974496d-01
      e(28)= 5.1866259d+01
      e(29)= 1.5554932d+01
      e(30)= 4.7228707d+00
      c(28)= 1.71727d-01
      c(29)= 5.328346d-01
      c(30)= 5.567749d-01
      e(31)= 4.9230485d+00
      e(32)= 1.1881975d+00
      e(33)= 5.2307149d-01
      c(31)=-3.845678d-01
      c(32)= 8.484104d-01
      c(33)= 3.848926d-01
      e(34)= 1.8782203d+00
      e(35)= 9.0272768d-01
      e(36)= 3.7578948d-01
      c(34)= 4.042684d-01
      c(35)= 5.182477d-01
      c(36)= 1.55614d-01
      e(37)= 2.3358427d+00
      e(38)= 8.691822d-01
      e(39)= 2.8462575d-01
      c(37)= 2.897666d-01
      c(38)= 5.607057d-01
      c(39)= 3.443139d-01
      e(40)= 8.8981019d-01
      e(41)= 1.0617371d-01
      e(42)= 3.788171d-02
      c(40)=-2.042568d-01
      c(41)= 6.180411d-01
      c(42)= 4.942361d-01
      return
c
c     mercury 1-s
c
 1000 continue
      e(1)= 3.3404144d+04
      e(2)= 5.054976d+03
      e(3)= 1.1129314d+03
      c(1)= 5.97939d-02
      c(2)= 3.622279d-01
      c(3)= 6.95561d-01
      e(4)= 1.4756246d+03
      e(5)= 1.4640265d+02
      e(6)= 6.1389131d+01
      c(4)=-1.124854d-01
      c(5)= 7.439363d-01
      c(6)= 3.235869d-01
      e(7)= 2.3328646d+03
      e(8)= 5.491452d+02
      e(9)= 1.656402d+02
      c(7)= 7.47191d-02
      c(8)= 4.088872d-01
      c(9)= 6.53549d-01
      e(10)= 1.3621257d+02
      e(11)= 2.6696685d+01
      e(12)= 1.1554849d+01
      c(10)=-2.979513d-01
      c(11)= 9.145658d-01
      c(12)= 2.583104d-01
      e(13)= 3.0619849d+02
      e(14)= 6.3373202d+01
      e(15)= 2.6012143d+01
      c(13)=-3.75406d-02
      c(14)= 3.880999d-01
      c(15)= 6.777274d-01
      e(16)= 3.7583702d+02
      e(17)= 1.0848287d+02
      e(18)= 3.6833664d+01
      c(16)= 1.009775d-01
      c(17)= 4.667307d-01
      c(18)= 6.006099d-01
      e(19)= 2.5982997d+01
      e(20)= 6.5233802d+00
      e(21)= 3.2362472d+00
      c(19)=-3.806873d-01
      c(20)= 8.614964d-01
      c(21)= 3.719927d-01
      e(22)= 1.1965861d+01
      e(23)= 5.3228262d+00
      e(24)= 2.1625624d+00
      c(22)= 4.339413d-01
      c(23)= 5.593747d-01
      c(24)= 7.59121d-02
      e(25)= 2.1948085d+01
      e(26)= 8.7311154d+00
      e(27)= 3.6993476d+00
      c(25)= 2.83437d-01
      c(26)= 5.638837d-01
      c(27)= 2.947653d-01
      e(28)= 5.420784d+01
      e(29)= 1.6330801d+01
      e(30)= 5.0053841d+00
      c(28)= 1.715987d-01
      c(29)= 5.338727d-01
      c(30)= 5.533129d-01
      e(31)= 5.3106531d+00
      e(32)= 1.2593063d+00
      e(33)= 5.4614633d-01
      c(31)=-3.956055d-01
      c(32)= 8.415402d-01
      c(33)= 3.937369d-01
      e(34)= 1.9253031d+00
      e(35)= 8.7120263d-01
      e(36)= 3.6365398d-01
      c(34)= 4.572396d-01
      c(35)= 5.06407d-01
      c(36)= 1.138995d-01
      e(37)= 2.3911943d+00
      e(38)= 8.9025734d-01
      e(39)= 2.9559019d-01
      c(37)= 3.162181d-01
      c(38)= 5.584083d-01
      c(39)= 3.179613d-01
      e(40)= 8.9197217d-01
      e(41)= 1.0752428d-01
      e(42)= 3.93799d-02
      c(40)=-1.998916d-01
      c(41)= 6.354779d-01
      c(42)= 4.747992d-01
      return
      end
      subroutine hlanth(e,c,nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),c(*)
c
c     huzinaga lanthanides...all are 4c(n),5c(0),6c(2) states.
c     note also that the f expansions hold 4 gaussians.
c
      ihop = nucz-56
      go to (100,200,300,400,500,600,700,800,900,1000,
     *       1100,1200,1300,1400), ihop
c
c     lanthanum 2-f
c
  100 continue
      e(1)= 1.662544d+04
      e(2)= 2.5192937d+03
      e(3)= 5.546769d+02
      c(1)= 6.09359d-02
      c(2)= 3.661543d-01
      c(3)= 6.914798d-01
      e(4)= 7.3806161d+02
      e(5)= 7.2542668d+01
      e(6)= 3.0332746d+01
      c(4)=-1.105265d-01
      c(5)= 7.265141d-01
      c(6)= 3.408799d-01
      e(7)= 1.100734d+03
      e(8)= 2.581062d+02
      e(9)= 7.7064873d+01
      c(7)= 7.95299d-02
      c(8)= 4.204338d-01
      c(9)= 6.419901d-01
      e(10)= 6.5022157d+01
      e(11)= 1.188051d+01
      e(12)= 5.0799848d+00
      c(10)=-2.778738d-01
      c(11)= 8.951231d-01
      c(12)= 2.629053d-01
      e(13)= 1.5729247d+02
      e(14)= 2.7442045d+01
      e(15)= 1.0932262d+01
      c(13)=-3.33548d-02
      c(14)= 4.057066d-01
      c(15)= 6.616068d-01
      e(16)= 1.5717832d+02
      e(17)= 4.439378d+01
      e(18)= 1.4460783d+01
      c(16)= 1.120244d-01
      c(17)= 4.782146d-01
      c(18)= 5.922607d-01
      e(19)= 1.0833873d+01
      e(20)= 2.6683331d+00
      e(21)= 1.279477d+00
      c(19)=-3.775915d-01
      c(20)= 8.430458d-01
      c(21)= 3.86647d-01
      e(22)= 4.597872d+00
      e(23)= 2.0136263d+00
      e(24)= 7.7692443d-01
      c(22)= 4.39394d-01
      c(23)= 5.57974d-01
      c(24)= 7.52893d-02
      e(25)= 7.932049d+00
      e(26)= 2.9378166d+00
      e(27)= 1.123846d+00
      c(25)= 2.829726d-01
      c(26)= 5.906859d-01
      c(27)= 2.872377d-01
      e(28)= 2.0541754d+00
      e(29)= 4.2145344d-01
      e(30)= 1.8969742d-01
      c(28)=-3.560968d-01
      c(29)= 8.367243d-01
      c(30)= 3.610625d-01
      e(31)= 6.3938885d-01
      e(32)= 2.7659869d-01
      e(33)= 1.0974625d-01
      c(31)= 4.469577d-01
      c(32)= 5.338706d-01
      c(33)= 9.96811d-02
      e(34)= 2.3231033d+01
      e(35)= 6.8221675d+00
      e(36)= 2.1843543d+00
      e(37)= 6.2500149d-01
      c(34)= 8.48026d-02
      c(35)= 3.185478d-01
      c(36)= 5.215605d-01
      c(37)= 4.458382d-01
      e(38)= 2.1333194d-01
      e(39)= 4.33102d-02
      e(40)= 1.832129d-02
      c(38)=-1.919227d-01
      c(39)= 7.0088d-01
      c(40)= 4.300542d-01
      return
c
c     cerium 3-h
c
  200 continue
      e(1)= 1.7261118d+04
      e(2)= 2.6142239d+03
      e(3)= 5.7525356d+02
      c(1)= 6.07638d-02
      c(2)= 3.656823d-01
      c(3)= 6.920432d-01
      e(4)= 7.6572152d+02
      e(5)= 7.5195756d+01
      e(6)= 3.1381854d+01
      c(4)=-1.106185d-01
      c(5)= 7.27871d-01
      c(6)= 3.39605d-01
      e(7)= 1.1410638d+03
      e(8)= 2.6781001d+02
      e(9)= 8.0147542d+01
      c(7)= 7.96317d-02
      c(8)= 4.202625d-01
      c(9)= 6.417642d-01
      e(10)= 6.7475838d+01
      e(11)= 1.2438016d+01
      e(12)= 5.3575335d+00
      c(10)=-2.798715d-01
      c(11)= 8.931133d-01
      c(12)= 2.662826d-01
      e(13)= 1.5936103d+02
      e(14)= 2.9022085d+01
      e(15)= 1.156846d+01
      c(13)=-3.32807d-02
      c(14)= 3.970831d-01
      c(15)= 6.699493d-01
      e(16)= 1.641335d+02
      e(17)= 4.6508119d+01
      e(18)= 1.5198705d+01
      c(16)= 1.118343d-01
      c(17)= 4.779431d-01
      c(18)= 5.918922d-01
      e(19)= 1.1413568d+01
      e(20)= 2.7973046d+00
      e(21)= 1.3446362d+00
      c(19)=-3.759205d-01
      c(20)= 8.430916d-01
      c(21)= 3.85076d-01
      e(22)= 4.8235492d+00
      e(23)= 2.1015449d+00
      e(24)= 7.8594489d-01
      c(22)= 4.438528d-01
      c(23)= 5.57609d-01
      c(24)= 7.19006d-02
      e(25)= 8.4807186d+00
      e(26)= 3.1599453d+00
      e(27)= 1.2103814d+00
      c(25)= 2.785615d-01
      c(26)= 5.871626d-01
      c(27)= 2.949838d-01
      e(28)= 2.2097198d+00
      e(29)= 4.3697325d-01
      e(30)= 1.9857828d-01
      c(28)=-3.554107d-01
      c(29)= 8.258265d-01
      c(30)= 3.659862d-01
      e(31)= 6.6808886d-01
      e(32)= 2.889284d-01
      e(33)= 1.1363078d-01
      c(31)= 4.45396d-01
      c(32)= 5.33501d-01
      c(33)= 1.027921d-01
      e(34)= 2.6139561d+01
      e(35)= 7.7491838d+00
      e(36)= 2.5050383d+00
      e(37)= 7.3333938d-01
      c(34)= 8.09736d-02
      c(35)= 3.139003d-01
      c(36)= 5.232838d-01
      c(37)= 4.426836d-01
      e(38)= 2.245275d-01
      e(39)= 4.484319d-02
      e(40)= 1.896978d-02
      c(38)=-2.004514d-01
      c(39)= 6.844032d-01
      c(40)= 4.488729d-01
      return
c
c     praseodymium 4-i
c
  300 continue
      e(1)= 1.7886843d+04
      e(2)= 2.7083926d+03
      e(3)= 5.9592268d+02
      c(1)= 6.0681d-02
      c(2)= 3.654327d-01
      c(3)= 6.923172d-01
      e(4)= 7.9303026d+02
      e(5)= 7.7889637d+01
      e(6)= 3.2526404d+01
      c(4)=-1.10719d-01
      c(5)= 7.292466d-01
      c(6)= 3.381817d-01
      e(7)= 1.1917574d+03
      e(8)= 2.7946559d+02
      e(9)= 8.3496655d+01
      c(7)= 7.87516d-02
      c(8)= 4.18651d-01
      c(9)= 6.439389d-01
      e(10)= 7.0100725d+01
      e(11)= 1.2959578d+01
      e(12)= 5.5771203d+00
      c(10)=-2.809137d-01
      c(11)= 8.959656d-01
      c(12)= 2.641582d-01
      e(13)= 1.6474089d+02
      e(14)= 3.0359081d+01
      e(15)= 1.2121698d+01
      c(13)=-3.387d-02
      c(14)= 3.95047d-01
      c(15)= 6.719295d-01
      e(16)= 1.7177896d+02
      e(17)= 4.8703353d+01
      e(18)= 1.5963562d+01
      c(16)= 1.113158d-01
      c(17)= 4.777975d-01
      c(18)= 5.91725d-01
      e(19)= 1.1952788d+01
      e(20)= 2.9175114d+00
      e(21)= 1.3900121d+00
      c(19)=-3.748941d-01
      c(20)= 8.546478d-01
      c(21)= 3.733135d-01
      e(22)= 5.1385128d+00
      e(23)= 2.2538269d+00
      e(24)= 8.7313773d-01
      c(22)= 4.321463d-01
      c(23)= 5.611356d-01
      c(24)= 8.01383d-02
      e(25)= 8.9281848d+00
      e(26)= 3.331993d+00
      e(27)= 1.2748513d+00
      c(25)= 2.809806d-01
      c(26)= 5.869343d-01
      c(27)= 2.928538d-01
      e(28)= 2.3097738d+00
      e(29)= 4.5748947d-01
      e(30)= 2.0271971d-01
      c(28)=-3.516056d-01
      c(29)= 8.231239d-01
      c(30)= 3.686808d-01
      e(31)= 7.0148489d-01
      e(32)= 3.0047369d-01
      e(33)= 1.2131011d-01
      c(31)= 4.422667d-01
      c(32)= 5.336878d-01
      c(33)= 1.066328d-01
      e(34)= 2.7475121d+01
      e(35)= 8.1971451d+00
      e(36)= 2.6760726d+00
      e(37)= 7.9088718d-01
      c(34)= 8.42526d-02
      c(35)= 3.215852d-01
      c(36)= 5.239577d-01
      c(37)= 4.293941d-01
      e(38)= 2.3370796d-01
      e(39)= 4.605482d-02
      e(40)= 1.977678d-02
      c(38)=-2.047518d-01
      c(39)= 6.691118d-01
      c(40)= 4.637875d-01
      return
c
c     neodymium 5-i
c
  400 continue
      e(1)= 1.8526642d+04
      e(2)= 2.8048423d+03
      e(3)= 6.1704358d+02
      c(1)= 6.05876d-02
      c(2)= 3.651271d-01
      c(3)= 6.926526d-01
      e(4)= 8.2124079d+02
      e(5)= 8.0575481d+01
      e(6)= 3.3535519d+01
      c(4)=-1.107753d-01
      c(5)= 7.31765d-01
      c(6)= 3.357808d-01
      e(7)= 1.2345803d+03
      e(8)= 2.8969646d+02
      e(9)= 8.6648623d+01
      c(7)= 7.87454d-02
      c(8)= 4.18602d-01
      c(9)= 6.438205d-01
      e(10)= 7.2745461d+01
      e(11)= 1.3509141d+01
      e(12)= 5.8627242d+00
      c(10)=-2.823472d-01
      c(11)= 8.966473d-01
      c(12)= 2.641148d-01
      e(13)= 1.6858307d+02
      e(14)= 3.1866276d+01
      e(15)= 1.2721128d+01
      c(13)=-3.40075d-02
      c(14)= 3.909916d-01
      c(15)= 6.759948d-01
      e(16)= 1.796566d+02
      e(17)= 5.0987916d+01
      e(18)= 1.6747778d+01
      c(16)= 1.106871d-01
      c(17)= 4.774459d-01
      c(18)= 5.919556d-01
      e(19)= 1.2536668d+01
      e(20)= 3.0529657d+00
      e(21)= 1.4459091d+00
      c(19)=-3.73329d-01
      c(20)= 8.556638d-01
      c(21)= 3.715675d-01
      e(22)= 5.5014686d+00
      e(23)= 2.441901d+00
      e(24)= 1.0013475d+00
      c(22)= 4.135184d-01
      c(23)= 5.644878d-01
      c(24)= 9.51912d-02
      e(25)= 9.4309935d+00
      e(26)= 3.5273865d+00
      e(27)= 1.3506401d+00
      c(25)= 2.808573d-01
      c(26)= 5.852852d-01
      c(27)= 2.94646d-01
      e(28)= 2.3675661d+00
      e(29)= 4.8032361d-01
      e(30)= 2.1321496d-01
      c(28)=-3.452861d-01
      c(29)= 8.094339d-01
      c(30)= 3.827644d-01
      e(31)= 7.4285559d-01
      e(32)= 3.2091091d-01
      e(33)= 1.2917825d-01
      c(31)= 4.228918d-01
      c(32)= 5.407912d-01
      c(33)= 1.205496d-01
      e(34)= 2.9610554d+01
      e(35)= 8.8796987d+00
      e(36)= 2.911521d+00
      e(37)= 8.6225689d-01
      c(34)= 8.36033d-02
      c(35)= 3.215007d-01
      c(36)= 5.245325d-01
      c(37)= 4.276097d-01
      e(38)= 2.4169741d-01
      e(39)= 4.856361d-02
      e(40)= 2.045286d-02
      c(38)=-2.123201d-01
      c(39)= 6.410192d-01
      c(40)= 4.966707d-01
      return
c
c     promethium 6-h
c
  500 continue
      e(1)= 1.9159059d+04
      e(2)= 2.9008327d+03
      e(3)= 6.3827198d+02
      c(1)= 6.05619d-02
      c(2)= 3.650056d-01
      c(3)= 6.927569d-01
      e(4)= 8.4958277d+02
      e(5)= 8.3341451d+01
      e(6)= 3.4625265d+01
      c(4)=-1.108554d-01
      c(5)= 7.336528d-01
      c(6)= 3.339623d-01
      e(7)= 1.2809289d+03
      e(8)= 3.0067027d+02
      e(9)= 8.9959488d+01
      c(7)= 7.84782d-02
      c(8)= 4.180326d-01
      c(9)= 6.444278d-01
      e(10)= 7.5520257d+01
      e(11)= 1.4029012d+01
      e(12)= 6.0294403d+00
      c(10)=-2.828345d-01
      c(11)= 9.021283d-01
      c(12)= 2.592573d-01
      e(13)= 1.7515156d+02
      e(14)= 3.3032232d+01
      e(15)= 1.3226298d+01
      c(13)=-3.42461d-02
      c(14)= 3.940065d-01
      c(15)= 6.729218d-01
      e(16)= 1.8722437d+02
      e(17)= 5.32178d+01
      e(18)= 1.7525139d+01
      c(16)= 1.104956d-01
      c(17)= 4.776191d-01
      c(18)= 5.913436d-01
      e(19)= 1.3179882d+01
      e(20)= 3.1664005d+00
      e(21)= 1.5049793d+00
      c(19)=-3.692896d-01
      c(20)= 8.617535d-01
      c(21)= 3.618931d-01
      e(22)= 5.7636903d+00
      e(23)= 2.5431625d+00
      e(24)= 1.0364114d+00
      c(22)= 4.174654d-01
      c(23)= 5.6393d-01
      c(24)= 9.21161d-02
      e(25)= 9.8543016d+00
      e(26)= 3.6872548d+00
      e(27)= 1.4117181d+00
      c(25)= 2.853567d-01
      c(26)= 5.852009d-01
      c(27)= 2.902922d-01
      e(28)= 2.5039828d+00
      e(29)= 5.0626814d-01
      e(30)= 2.2606303d-01
      c(28)=-3.467544d-01
      c(29)= 7.774164d-01
      c(30)= 4.133012d-01
      e(31)= 7.7440027d-01
      e(32)= 3.3432315d-01
      e(33)= 1.3759787d-01
      c(31)= 4.197609d-01
      c(32)= 5.371351d-01
      c(33)= 1.280269d-01
      e(34)= 3.0839766d+01
      e(35)= 9.3015739d+00
      e(36)= 3.0627573d+00
      e(37)= 9.081981d-01
      c(34)= 8.66158d-02
      c(35)= 3.277359d-01
      c(36)= 5.24085d-01
      c(37)= 4.196829d-01
      e(38)= 2.4982112d-01
      e(39)= 4.97224d-02
      e(40)= 2.11403d-02
      c(38)=-2.139898d-01
      c(39)= 6.289105d-01
      c(40)= 5.077477d-01
      return
c
c     samarium 7-f
c
  600 continue
      e(1)= 1.9842708d+04
      e(2)= 3.0021144d+03
      e(3)= 6.60288d+02
      c(1)= 6.04129d-02
      c(2)= 3.646398d-01
      c(3)= 6.932013d-01
      e(4)= 8.7904586d+02
      e(5)= 8.6285511d+01
      e(6)= 3.6009643d+01
      c(4)=-1.109973d-01
      c(5)= 7.328357d-01
      c(6)= 3.345842d-01
      e(7)= 1.329827d+03
      e(8)= 3.1219809d+02
      e(9)= 9.345362d+01
      c(7)= 7.81022d-02
      c(8)= 4.170192d-01
      c(9)= 6.45485d-01
      e(10)= 7.8329115d+01
      e(11)= 1.4573156d+01
      e(12)= 6.2287472d+00
      c(10)=-2.83507d-01
      c(11)= 9.054729d-01
      c(12)= 2.565598d-01
      e(13)= 1.8254042d+02
      e(14)= 3.4346164d+01
      e(15)= 1.3784402d+01
      c(13)=-3.45341d-02
      c(14)= 3.937909d-01
      c(15)= 6.729263d-01
      e(16)= 1.9529689d+02
      e(17)= 5.5593987d+01
      e(18)= 1.834806d+01
      c(16)= 1.09991d-01
      c(17)= 4.771221d-01
      c(18)= 5.916223d-01
      e(19)= 1.3752721d+01
      e(20)= 3.3563672d+00
      e(21)= 1.5718294d+00
      c(19)=-3.70342d-01
      c(20)= 8.468583d-01
      c(21)= 3.795108d-01
      e(22)= 6.1634636d+00
      e(23)= 2.729236d+00
      e(24)= 1.0706707d+00
      c(22)= 3.975439d-01
      c(23)= 5.769811d-01
      c(24)= 1.012375d-01
      e(25)= 1.0389082d+01
      e(26)= 3.8703504d+00
      e(27)= 1.4811196d+00
      c(25)= 2.859941d-01
      c(26)= 5.863899d-01
      c(27)= 2.888684d-01
      e(28)= 2.6439751d+00
      e(29)= 5.2293987d-01
      e(30)= 2.3610458d-01
      c(28)=-3.54785d-01
      c(29)= 7.671895d-01
      c(30)= 4.223481d-01
      e(31)= 8.0855734d-01
      e(32)= 3.4537441d-01
      e(33)= 1.4370986d-01
      c(31)= 4.152969d-01
      c(32)= 5.410361d-01
      c(33)= 1.292299d-01
      e(34)= 3.3212664d+01
      e(35)= 1.004506d+01
      e(36)= 3.3169803d+00
      e(37)= 9.7947735d-01
      c(34)= 8.51107d-02
      c(35)= 3.257303d-01
      c(36)= 5.249144d-01
      c(37)= 4.214426d-01
      e(38)= 2.5807912d-01
      e(39)= 5.136601d-02
      e(40)= 2.183911d-02
      c(38)=-2.15899d-01
      c(39)= 6.094654d-01
      c(40)= 5.271336d-01
      return
c
c     europium 8-s
c
  700 continue
      e(1)= 2.0507003d+04
      e(2)= 3.1028009d+03
      e(3)= 6.8235369d+02
      c(1)= 6.03526d-02
      c(2)= 3.644263d-01
      c(3)= 6.934313d-01
      e(4)= 9.0788537d+02
      e(5)= 8.9668374d+01
      e(6)= 3.8112692d+01
      c(4)=-1.113376d-01
      c(5)= 7.253497d-01
      c(6)= 3.414365d-01
      e(7)= 1.3845829d+03
      e(8)= 3.2459262d+02
      e(9)= 9.7074434d+01
      c(7)= 7.7376d-02
      c(8)= 4.157355d-01
      c(9)= 6.47184d-01
      e(10)= 8.0647184d+01
      e(11)= 1.5647728d+01
      e(12)= 7.4065456d+00
      c(10)=-2.893875d-01
      c(11)= 8.501431d-01
      c(12)= 3.140286d-01
      e(13)= 1.8869784d+02
      e(14)= 3.5858186d+01
      e(15)= 1.4439022d+01
      c(13)=-3.50543d-02
      c(14)= 3.895085d-01
      c(15)= 6.767807d-01
      e(16)= 2.059212d+02
      e(17)= 5.8526401d+01
      e(18)= 1.9299158d+01
      c(16)= 1.076601d-01
      c(17)= 4.744727d-01
      c(18)= 5.956181d-01
      e(19)= 1.4082371d+01
      e(20)= 3.9851714d+00
      e(21)= 2.1380518d+00
      c(19)=-3.914035d-01
      c(20)= 5.756313d-01
      c(21)= 6.601493d-01
      e(22)= 7.0626748d+00
      e(23)= 3.4820688d+00
      e(24)= 1.6742596d+00
      c(22)= 3.001791d-01
      c(23)= 5.330359d-01
      c(24)= 2.383721d-01
      e(25)= 1.1011484d+01
      e(26)= 4.1362487d+00
      e(27)= 1.5896167d+00
      c(25)= 2.799785d-01
      c(26)= 5.796519d-01
      c(27)= 3.014955d-01
      e(28)= 4.54898d+00
      e(29)= 7.0187288d-01
      e(30)= 3.5287489d-01
      c(28)=-4.1773d-03
      c(29)= 2.346292d-01
      c(30)= 7.822717d-01
      e(31)= 1.9741929d+00
      e(32)= 7.7456119d-01
      e(33)= 2.7333192d-01
      c(31)=-1.23731d-02
      c(32)= 5.332838d-01
      c(33)= 5.520593d-01
      e(34)= 3.7125837d+01
      e(35)= 1.1248423d+01
      e(36)= 3.7191113d+00
      e(37)= 1.0960099d+00
      c(34)= 7.85951d-02
      c(35)= 3.140544d-01
      c(36)= 5.258441d-01
      c(37)= 4.342082d-01
      e(38)= 8.3915121d-01
      e(39)= 6.491549d-02
      e(40)= 2.469337d-02
      c(38)=-8.19894d-02
      c(39)=-3.165455d-01
      c(40)=-6.935623d-01
      return
c
c     gadolinium 7-f
c
  800 continue
      e(1)= 2.116154d+04
      e(2)= 3.2022178d+03
      e(3)= 7.0445174d+02
      c(1)= 6.03678d-02
      c(2)= 3.644456d-01
      c(3)= 6.933749d-01
      e(4)= 9.3554786d+02
      e(5)= 9.2934461d+01
      e(6)= 3.9928439d+01
      c(4)=-1.115875d-01
      c(5)= 7.215131d-01
      c(6)= 3.449564d-01
      e(7)= 1.4301271d+03
      e(8)= 3.3566767d+02
      e(9)= 1.0050588d+02
      c(7)= 7.74131d-02
      c(8)= 4.155803d-01
      c(9)= 6.471221d-01
      e(10)= 8.3470797d+01
      e(11)= 1.6165723d+01
      e(12)= 7.452228d+00
      c(10)=-2.897088d-01
      c(11)= 8.635201d-01
      c(12)= 3.018628d-01
      e(13)= 1.9218592d+02
      e(14)= 3.7398973d+01
      e(15)= 1.5032371d+01
      c(13)=-3.50449d-02
      c(14)= 3.886697d-01
      c(15)= 6.779961d-01
      e(16)= 2.1333543d+02
      e(17)= 6.0761796d+01
      e(18)= 2.0097751d+01
      c(16)= 1.080439d-01
      c(17)= 4.753318d-01
      c(18)= 5.939344d-01
      e(19)= 1.4303089d+01
      e(20)= 4.0502079d+00
      e(21)= 2.1005199d+00
      c(19)=-3.962516d-01
      c(20)= 6.657319d-01
      c(21)= 5.798678d-01
      e(22)= 7.2245576d+00
      e(23)= 3.5923145d+00
      e(24)= 1.7459549d+00
      c(22)= 3.154012d-01
      c(23)= 5.187296d-01
      c(24)= 2.372441d-01
      e(25)= 1.1465387d+01
      e(26)= 4.3293004d+00
      e(27)= 1.6672499d+00
      c(25)= 2.828563d-01
      c(26)= 5.769317d-01
      c(27)= 3.010715d-01
      e(28)= 4.6945382d+00
      e(29)= 7.2433139d-01
      e(30)= 3.7987327d-01
      c(28)=-1.15343d-02
      c(29)= 1.844031d-01
      c(30)= 8.315934d-01
      e(31)= 2.0373631d+00
      e(32)= 7.735558d-01
      e(33)= 2.7298819d-01
      c(31)=-1.0986d-03
      c(32)= 5.498008d-01
      c(33)= 5.283374d-01
      e(34)= 3.8819868d+01
      e(35)= 1.1714379d+01
      e(36)= 3.8626271d+00
      e(37)= 1.1414087d+00
      c(34)= 8.02236d-02
      c(35)= 3.202909d-01
      c(36)= 5.238402d-01
      c(37)= 4.302543d-01
      e(38)= 8.6600235d-01
      e(39)= 6.699265d-02
      e(40)= 2.548351d-02
      c(38)= 7.47141d-02
      c(39)= 2.936293d-01
      c(40)= 7.177832d-01
      return
c
c     terbium 6-h
c
  900 continue
      e(1)= 2.1835019d+04
      e(2)= 3.3056159d+03
      e(3)= 7.2731432d+02
      c(1)= 6.03439d-02
      c(2)= 3.642347d-01
      c(3)= 6.935561d-01
      e(4)= 9.6508624d+02
      e(5)= 9.568196d+01
      e(6)= 4.0815799d+01
      c(4)=-1.115621d-01
      c(5)= 7.271652d-01
      c(6)= 3.395253d-01
      e(7)= 1.4802159d+03
      e(8)= 3.4751808d+02
      e(9)= 1.0411187d+02
      c(7)= 7.71847d-02
      c(8)= 4.150012d-01
      c(9)= 6.476882d-01
      e(10)= 8.6613515d+01
      e(11)= 1.6578152d+01
      e(12)= 7.8003711d+00
      c(10)=-2.901632d-01
      c(11)= 8.787749d-01
      c(12)= 2.847502d-01
      e(13)= 2.0260666d+02
      e(14)= 3.8524834d+01
      e(15)= 1.5570199d+01
      c(13)=-3.5362d-02
      c(14)= 3.920956d-01
      c(15)= 6.740818d-01
      e(16)= 2.2188313d+02
      e(17)= 6.3252113d+01
      e(18)= 2.0957007d+01
      c(16)= 1.076858d-01
      c(17)= 4.752801d-01
      c(18)= 5.938155d-01
      e(19)= 1.577566d+01
      e(20)= 4.0878441d+00
      e(21)= 2.1592946d+00
      c(19)=-3.703324d-01
      c(20)= 6.57369d-01
      c(21)= 5.630353d-01
      e(22)= 7.7988359d+00
      e(23)= 3.7575148d+00
      e(24)= 1.763071d+00
      c(22)= 2.983771d-01
      c(23)= 5.507163d-01
      c(24)= 2.242048d-01
      e(25)= 1.2013693d+01
      e(26)= 4.5266838d+00
      e(27)= 1.7459841d+00
      c(25)= 2.839041d-01
      c(26)= 5.768815d-01
      c(27)= 3.001895d-01
      e(28)= 4.8423886d+00
      e(29)= 7.4714359d-01
      e(30)= 3.8423134d-01
      c(28)= 3.6219d-03
      c(29)= 2.157088d-01
      c(30)= 7.961193d-01
      e(31)= 2.1015281d+00
      e(32)= 7.8244623d-01
      e(33)= 2.7860506d-01
      c(31)= 8.2196d-03
      c(32)= 5.531276d-01
      c(33)= 5.178321d-01
      e(34)= 4.0680481d+01
      e(35)= 1.2337528d+01
      e(36)= 4.0635453d+00
      e(37)= 1.1838629d+00
      c(34)= 8.0844d-02
      c(35)= 3.214826d-01
      c(36)= 5.255148d-01
      c(37)= 4.297361d-01
      e(38)= 8.9327635d-01
      e(39)= 6.910253d-02
      e(40)= 2.628609d-02
      c(38)= 7.66729d-02
      c(39)= 2.782983d-01
      c(40)= 7.316191d-01
      return
c
c     dysprosium 5-i
c
 1000 continue
      e(1)= 2.2512038d+04
      e(2)= 3.4087516d+03
      e(3)= 7.5021737d+02
      c(1)= 6.03497d-02
      c(2)= 3.642008d-01
      c(3)= 6.935559d-01
      e(4)= 9.9666912d+02
      e(5)= 9.8727277d+01
      e(6)= 4.2083165d+01
      c(4)=-1.116446d-01
      c(5)= 7.284609d-01
      c(6)= 3.382378d-01
      e(7)= 1.5272025d+03
      e(8)= 3.5881603d+02
      e(9)= 1.0760264d+02
      c(7)= 7.72444d-02
      c(8)= 4.151054d-01
      c(9)= 6.474102d-01
      e(10)= 8.9147986d+01
      e(11)= 1.7526082d+01
      e(12)= 8.1030013d+00
      c(10)=-2.922379d-01
      c(11)= 8.565066d-01
      c(12)= 3.113017d-01
      e(13)= 2.068709d+02
      e(14)= 4.0194535d+01
      e(15)= 1.6221668d+01
      c(13)=-3.54586d-02
      c(14)= 3.894833d-01
      c(15)= 6.769116d-01
      e(16)= 2.3167736d+02
      e(17)= 6.6102671d+01
      e(18)= 2.1916854d+01
      c(16)= 1.065078d-01
      c(17)= 4.73606d-01
      c(18)= 5.958856d-01
      e(19)= 1.5711507d+01
      e(20)= 4.2363528d+00
      e(21)= 2.2096177d+00
      c(19)=-3.883642d-01
      c(20)= 7.041276d-01
      c(21)= 5.332091d-01
      e(22)= 7.8825557d+00
      e(23)= 3.7978777d+00
      e(24)= 1.8177367d+00
      c(22)= 3.234649d-01
      c(23)= 5.35763d-01
      c(24)= 2.136097d-01
      e(25)= 1.2545087d+01
      e(26)= 4.7409523d+00
      e(27)= 1.8318939d+00
      c(25)= 2.845373d-01
      c(26)= 5.74427d-01
      c(27)= 3.019803d-01
      e(28)= 4.8826712d+00
      e(29)= 7.7030946d-01
      e(30)= 3.9223241d-01
      c(28)=-1.00639d-02
      c(29)= 2.164299d-01
      c(30)= 8.016135d-01
      e(31)= 2.1190102d+00
      e(32)= 8.0074647d-01
      e(33)= 2.8941946d-01
      c(31)= 1.5719d-02
      c(32)= 5.448934d-01
      c(33)= 5.19235d-01
      e(34)= 4.2262676d+01
      e(35)= 1.2869755d+01
      e(36)= 4.2445315d+00
      e(37)= 1.2421218d+00
      c(34)= 8.24456d-02
      c(35)= 3.249725d-01
      c(36)= 5.236025d-01
      c(37)= 4.268708d-01
      e(38)= 9.2097321d-01
      e(39)= 7.124511d-02
      e(40)= 2.710112d-02
      c(38)= 7.84492d-02
      c(39)= 2.655523d-01
      c(40)= 7.430275d-01
      return
c
c     holmium 4-i
c
 1100 continue
      e(1)= 2.3180533d+04
      e(2)= 3.5116702d+03
      e(3)= 7.7332413d+02
      c(1)= 6.04086d-02
      c(2)= 3.642593d-01
      c(3)= 6.934164d-01
      e(4)= 1.0292619d+03
      e(5)= 1.0170128d+02
      e(6)= 4.3328224d+01
      c(4)=-1.116859d-01
      c(5)= 7.311447d-01
      c(6)= 3.354456d-01
      e(7)= 1.5737315d+03
      e(8)= 3.699692d+02
      e(9)= 1.1106992d+02
      c(7)= 7.73985d-02
      c(8)= 4.155441d-01
      c(9)= 6.467612d-01
      e(10)= 9.2657169d+01
      e(11)= 1.7937269d+01
      e(12)= 8.477776d+00
      c(10)=-2.915066d-01
      c(11)= 8.703011d-01
      c(12)= 2.947321d-01
      e(13)= 2.0941284d+02
      e(14)= 4.1838408d+01
      e(15)= 1.6853075d+01
      c(13)=-3.54127d-02
      c(14)= 3.889577d-01
      c(15)= 6.778885d-01
      e(16)= 2.372149d+02
      e(17)= 6.7889543d+01
      e(18)= 2.2617285d+01
      c(16)= 1.085255d-01
      c(17)= 4.770762d-01
      c(18)= 5.904551d-01
      e(19)= 1.6692436d+01
      e(20)= 4.335891d+00
      e(21)= 2.2322713d+00
      c(19)=-3.779939d-01
      c(20)= 7.18871d-01
      c(21)= 5.097373d-01
      e(22)= 8.1634159d+00
      e(23)= 3.8567519d+00
      e(24)= 1.8073015d+00
      c(22)= 3.329288d-01
      c(23)= 5.475342d-01
      c(24)= 1.934593d-01
      e(25)= 1.3093652d+01
      e(26)= 4.9213964d+00
      e(27)= 1.8661708d+00
      c(25)= 2.858557d-01
      c(26)= 5.819554d-01
      c(27)= 2.947871d-01
      e(28)= 4.885196d+00
      e(29)= 7.9382901d-01
      e(30)= 3.9985838d-01
      c(28)=-6.7471d-03
      c(29)= 2.293095d-01
      c(30)= 7.884315d-01
      e(31)= 2.1201059d+00
      e(32)= 8.2864216d-01
      e(33)= 2.9950199d-01
      c(31)= 1.66251d-02
      c(32)= 5.383539d-01
      c(33)= 5.250268d-01
      e(34)= 4.4012473d+01
      e(35)= 1.341627d+01
      e(36)= 4.4270768d+00
      e(37)= 1.2896098d+00
      c(34)= 8.37922d-02
      c(35)= 3.283186d-01
      c(36)= 5.233255d-01
      c(37)= 4.242354d-01
      e(38)= 9.4909291d-01
      e(39)= 7.342041d-02
      e(40)= 2.792858d-02
      c(38)= 7.98273d-02
      c(39)= 2.501738d-01
      c(40)= 7.569843d-01
      return
c
c     erbium 3-h
c
 1200 continue
      e(1)= 2.3932986d+04
      e(2)= 3.6231807d+03
      e(3)= 7.9744649d+02
      c(1)= 6.02683d-02
      c(2)= 3.639567d-01
      c(3)= 6.938113d-01
      e(4)= 1.05807d+03
      e(5)= 1.0560404d+02
      e(6)= 4.5715154d+01
      c(4)=-1.120662d-01
      c(5)= 7.216806d-01
      c(6)= 3.446059d-01
      e(7)= 1.631852d+03
      e(8)= 3.8338971d+02
      e(9)= 1.1506026d+02
      c(7)= 7.68215d-02
      c(8)= 4.142441d-01
      c(9)= 6.483076d-01
      e(10)= 9.5445255d+01
      e(11)= 1.8624316d+01
      e(12)= 8.6540627d+00
      c(10)=-2.931741d-01
      c(11)= 8.716745d-01
      c(12)= 2.956139d-01
      e(13)= 2.1745963d+02
      e(14)= 4.3178334d+01
      e(15)= 1.7445255d+01
      c(13)=-3.56639d-02
      c(14)= 3.907248d-01
      c(15)= 6.759383d-01
      e(16)= 2.4776101d+02
      e(17)= 7.0854934d+01
      e(18)= 2.3602151d+01
      c(16)= 1.071608d-01
      c(17)= 4.757213d-01
      c(18)= 5.925282d-01
      e(19)= 1.7427739d+01
      e(20)= 4.4976911d+00
      e(21)= 2.3319556d+00
      c(19)=-3.754235d-01
      c(20)= 7.118474d-01
      c(21)= 5.141121d-01
      e(22)= 8.4116514d+00
      e(23)= 4.0197899d+00
      e(24)= 1.9188907d+00
      c(22)= 3.387219d-01
      c(23)= 5.328756d-01
      c(24)= 2.018587d-01
      e(25)= 1.3592062d+01
      e(26)= 5.1225912d+00
      e(27)= 1.9684932d+00
      c(25)= 2.891028d-01
      c(26)= 5.765507d-01
      c(27)= 2.96061d-01
      e(28)= 4.8804329d+00
      e(29)= 8.1161512d-01
      e(30)= 4.1312407d-01
      c(28)=-1.68679d-02
      c(29)= 2.22692d-01
      c(30)= 7.987453d-01
      e(31)= 2.1201431d+00
      e(32)= 8.4677323d-01
      e(33)= 2.9517372d-01
      c(31)= 1.38675d-02
      c(32)= 5.58206d-01
      c(33)= 5.114888d-01
      e(34)= 4.5343394d+01
      e(35)= 1.3864352d+01
      e(36)= 4.5763128d+00
      e(37)= 1.327826d+00
      c(34)= 8.62039d-02
      c(35)= 3.331165d-01
      c(36)= 5.228664d-01
      c(37)= 4.197267d-01
      e(38)= 9.923d-01
      e(39)= 7.562842d-02
      e(40)= 2.876849d-02
      c(38)= 7.67945d-02
      c(39)= 2.395752d-01
      c(40)= 7.680281d-01
      return
c
c     thulium 2-f
c
 1300 continue
      e(1)= 2.465641d+04
      e(2)= 3.7325062d+03
      e(3)= 8.2163815d+02
      c(1)= 6.02408d-02
      c(2)= 3.638374d-01
      c(3)= 6.939176d-01
      e(4)= 1.0914725d+03
      e(5)= 1.0854804d+02
      e(6)= 4.6629354d+01
      c(4)=-1.120302d-01
      c(5)= 7.267395d-01
      c(6)= 3.397656d-01
      e(7)= 1.6854211d+03
      e(8)= 3.9599853d+02
      e(9)= 1.1891754d+02
      c(7)= 7.66255d-02
      c(8)= 4.137828d-01
      c(9)= 6.487444d-01
      e(10)= 9.8840728d+01
      e(11)= 1.927499d+01
      e(12)= 8.9158303d+00
      c(10)=-2.924475d-01
      c(11)= 8.725193d-01
      c(12)= 2.946143d-01
      e(13)= 2.2782181d+02
      e(14)= 4.4528319d+01
      e(15)= 1.8046201d+01
      c(13)=-3.59085d-02
      c(14)= 3.923215d-01
      c(15)= 6.740286d-01
      e(16)= 2.5688245d+02
      e(17)= 7.3519651d+01
      e(18)= 2.4520879d+01
      c(16)= 1.068666d-01
      c(17)= 4.757448d-01
      c(18)= 5.923866d-01
      e(19)= 1.7990337d+01
      e(20)= 4.7415251d+00
      e(21)= 2.41551d+00
      c(19)=-3.78638d-01
      c(20)= 7.019566d-01
      c(21)= 5.28208d-01
      e(22)= 8.9691876d+00
      e(23)= 4.2464848d+00
      e(24)= 1.9503073d+00
      c(22)= 3.202135d-01
      c(23)= 5.533706d-01
      c(24)= 2.01616d-01
      e(25)= 1.4183386d+01
      e(26)= 5.334892d+00
      e(27)= 2.0382571d+00
      c(25)= 2.893569d-01
      c(26)= 5.783934d-01
      c(27)= 2.946384d-01
      e(28)= 4.9925787d+00
      e(29)= 8.3566168d-01
      e(30)= 4.1594385d-01
      c(28)=-1.2713d-02
      c(29)= 2.45873d-01
      c(30)= 7.758076d-01
      e(31)= 2.134923d+00
      e(32)= 8.6763438d-01
      e(33)= 3.1114728d-01
      c(31)= 2.19159d-02
      c(32)= 5.421133d-01
      c(33)= 5.18372d-01
      e(34)= 4.6858648d+01
      e(35)= 1.4376483d+01
      e(36)= 4.7433937d+00
      e(37)= 1.3719529d+00
      c(34)= 8.80212d-02
      c(35)= 3.368047d-01
      c(36)= 5.227706d-01
      c(37)= 4.159762d-01
      e(38)= 1.0216999d+00
      e(39)= 7.786914d-02
      e(40)= 2.962085d-02
      c(38)= 7.96995d-02
      c(39)= 2.235851d-01
      c(40)= 7.819576d-01
      return
c
c     ytterbium 1-s
c
 1400 continue
      e(1)= 2.5385694d+04
      e(2)= 3.8436032d+03
      e(3)= 8.4615684d+02
      c(1)= 6.02178d-02
      c(2)= 3.637205d-01
      c(3)= 6.940245d-01
      e(4)= 1.1234055d+03
      e(5)= 1.1169297d+02
      e(6)= 4.7233256d+01
      c(4)=-1.120119d-01
      c(5)= 7.308335d-01
      c(6)= 3.365089d-01
      e(7)= 1.7377854d+03
      e(8)= 4.0850276d+02
      e(9)= 1.2273231d+02
      c(7)= 7.65445d-02
      c(8)= 4.135681d-01
      c(9)= 6.489045d-01
      e(10)= 1.0178796d+02
      e(11)= 1.9813955d+01
      e(12)= 9.1947339d+00
      c(10)=-2.946531d-01
      c(11)= 8.851695d-01
      c(12)= 2.825866d-01
      e(13)= 2.3449301d+02
      e(14)= 4.6044797d+01
      e(15)= 1.8679254d+01
      c(13)=-3.60684d-02
      c(14)= 3.93155d-01
      c(15)= 6.732472d-01
      e(16)= 2.6549751d+02
      e(17)= 7.6107209d+01
      e(18)= 2.5437233d+01
      c(16)= 1.069945d-01
      c(17)= 4.761187d-01
      c(18)= 5.915046d-01
      e(19)= 1.8855724d+01
      e(20)= 4.8006027d+00
      e(21)= 2.4675257d+00
      c(19)=-3.703201d-01
      c(20)= 7.219045d-01
      c(21)= 5.008053d-01
      e(22)= 9.1334438d+00
      e(23)= 4.3362256d+00
      e(24)= 2.0173973d+00
      c(22)= 3.354143d-01
      c(23)= 5.418041d-01
      c(24)= 1.976589d-01
      e(25)= 1.4703007d+01
      e(26)= 5.5153232d+00
      e(27)= 2.0915122d+00
      c(25)= 2.922065d-01
      c(26)= 5.811741d-01
      c(27)= 2.89788d-01
      e(28)= 5.1383398d+00
      e(29)= 8.6005928d-01
      e(30)= 4.3522801d-01
      c(28)=-1.98422d-02
      c(29)= 2.196768d-01
      c(30)= 8.031232d-01
      e(31)= 2.1972533d+00
      e(32)= 8.7979164d-01
      e(33)= 3.1464744d-01
      c(31)= 2.66514d-02
      c(32)= 5.465593d-01
      c(33)= 5.11104d-01
      e(34)= 4.8314683d+01
      e(35)= 1.4868502d+01
      e(36)= 4.9219249d+00
      e(37)= 1.4244136d+00
      c(34)= 9.01533d-02
      c(35)= 3.403088d-01
      c(36)= 5.210745d-01
      c(37)= 4.128335d-01
      e(38)= 1.051529d+00
      e(39)= 8.014258d-02
      e(40)= 3.048564d-02
      c(38)= 7.47088d-02
      c(39)= 2.034147d-01
      c(40)= 8.021211d-01
      return
      end
      function hscale(nucz)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension row1(2),row2(8),row3(8),row4(8),row5(8)
      data row1/1.25d+00,1.00d+00/
      data row2/1.091d+00,1.105d+00,1.118d+00,1.075d+00,
     *          1.035d+00,1.014d+00,1.007d+00,1.000d+00/
      data row3/1.104d+00,1.091d+00,1.126d+00,1.104d+00,
     *          1.036d+00,1.014d+00,1.004d+00,1.000d+00/
      data row4/1.131d+00,1.086d+00,1.158d+00,1.116d+00,
     *          1.040d+00,1.017d+00,1.006d+00,1.000d+00/
      data row5/1.134d+00,1.078d+00,1.200d+00,1.131d+00,
     *          1.053d+00,1.024d+00,1.011d+00,1.000d+00/
      data done/1.0d+00/
c
c     provide default scale factors for huzinaga mini basis sets.
c     the main group scale factors are derived from rhf
c     geometry optimizations on sample molecules, such as
c     hydrides, carried out at north dakota state university
c     in the summer of 1985 by john deisz.
c     s and p scale factors are the same, which does not seem
c     to be a severe constaint. for a given element, the s scale
c     is sometimes smaller than the p in one molecule, larger in
c     another, and about equal in a third.
c     other elements (transition metals, z>54) are unscaled.
c
      if(nucz.le.2) then
         hscale=row1(nucz)
      else if(nucz.le.10) then
         nuc=nucz-2
         hscale=row2(nuc)
      else if(nucz.le.18) then
         nuc=nucz-10
         hscale=row3(nuc)
      else if(nucz.le.20) then
         nuc=nucz-18
         hscale=row4(nuc)
      else if(nucz.le.30) then
         hscale=done
      else if(nucz.le.36) then
         nuc=nucz-28
         hscale=row4(nuc)
      else if(nucz.le.38) then
         nuc=nucz-36
         hscale=row5(nuc)
      else if(nucz.le.48) then
         hscale=done
      else if(nucz.le.54) then
         nuc=nucz-46
         hscale=row5(nuc)
      else
       hscale=done
      endif
      return
      end
      subroutine hshell(nucz,ipass,omidi,itype,igauss,ng,
     *                  scale,scs,scp,scd,odone)
      implicit real*8 (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension ioff(15,15),joff(17,15),ktyp(15,3)
c
      data ktyp/1, 1,2, 1,2,3, 1,2,3,   1,2,3,  1,2*0,
     *          1, 1,2, 1,2,3, 1,2,3,   1,2, 4, 1,2*0,
     *          1, 1,2, 1,2,3, 1,2,3,4, 1,2,3,  1,2/
c                         pointers to the mini basis
      data ioff/0,  14*0,
     2          0, 3,  13*0,
     *          0, 3, 6,  12*0,
     3          0, 3, 9, 6,  11*0,
     *          0, 3, 9, 6,12,  10*0,
     4          0, 3,12, 6,15, 9,  9*0,
     *          0, 3,12, 6,15,18, 9,  8*0,
     *          0, 3,12, 6,15,21, 9,18,  7*0,
     5          0, 3,15, 6,18,24, 9,21,12,  6*0,
     *          0, 3,15, 6,18,24, 9,21,27,12,  5*0,
     *          0, 3,15, 6,18,27, 9,21,30,12,24,  4*0,
     6          0, 3, 6, 9,12,15,18,21,24,27,30,33,  3*0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,37, 2*0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39, 0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,42/
c                         pointers to the midi basis
      data joff/0, 2, 15*0,
     2          0, 3, 5,  14*0,
     *          0, 3, 5, 6, 8,  12*0,
     3          0, 3, 9, 6, 8,  12*0,
     *          0, 3, 9, 6, 8,12,14,  10*0,
     4          0, 3,12, 6,15, 9,11,  10*0,
     *          0, 3,12, 6,15,18,20, 9,11,  8*0,
     *          0, 3,12, 6,15,21, 9,11,18,20,  7*0,
     5          0, 3,15, 6,18,24, 9,21,12,14,  7*0,
     *          0, 3,15, 6,18,24, 9,21,27,29,12,14,  5*0,
     *          0, 3,15, 6,18,27, 9,21,30,12,14,24,26,  4*0,
     6          0, 3, 6, 9,12,15,18,21,24,27,30,33,35,  4*0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,36,37,39,  2*0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,36,38,39,41, 0,
     *          0, 3, 6, 9,12,15,18,21,24,27,30,33,36,39,41,42,44/
      data dzero,done/0.0d+00,1.0d+00/
c
c     set values for the current huzinaga mini/midi shell
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3,4 for s,p,d,f shell
c            ng     = offset in eex,coef arrays for current shell
c            scale  = scaling factor for current shell
c
      kind=1
      if(nucz.gt. 2) kind=2
      if(nucz.gt. 4) kind=3
      if(nucz.gt.10) kind=4
      if(nucz.gt.12) kind=5
      if(nucz.gt.18) kind=6
      if(nucz.gt.20) kind=7
      if(nucz.gt.30) kind=8
      if(nucz.gt.36) kind=9
      if(nucz.gt.38) kind=10
      if(nucz.gt.48) kind=11
      if(nucz.gt.54) kind=12
      if(nucz.gt.56) kind=13
      if(nucz.gt.70) kind=14
      if(nucz.gt.80) kind=15
      if(nucz.gt.86) then
         call caserr2('huzinaga basis sets only extend to radon')
      end if
c
      kind2 = 1
      if(nucz.gt.56) kind2=2
      if(nucz.gt.70) kind2=3
c
      mxpass=kind
c                 all atoms have two valence shells, except alkalis
      nsplit=2
      if(kind.le.2  .or.  kind.eq.4  .or.  kind.eq.6  .or.
     *   kind.eq.9  .or.  kind.eq.12) nsplit=1
      if(omidi) mxpass=mxpass+nsplit
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c              ******
c
      igauss = 3
      if(omidi) then
         mi = mxpass-ipass
         if(mi.lt.2*nsplit) then
            if(mod(mi,2).eq.1) igauss=2
            if(mod(mi,2).eq.0) igauss=1
            if(mi.ge.2) then
               ityp=ktyp(kind-1,kind2)
            else
               ityp=ktyp(kind,kind2)
            end if
         else
            ityp=ktyp(ipass,kind2)
         end if
         ng  =joff(ipass,kind)
      else
         ityp=ktyp(ipass,kind2)
         ng  =ioff(ipass,kind)
      end if
c                lanthanides have 4 (or 3,1) gaussians in valence f
      if(kind2.eq.2  .and.  ityp.eq.4  .and.  igauss.ne.1)
     *   igauss=igauss+1
c                alkali's don't yet fill the previous d level
      if(kind.eq.6  .and.  ipass.ge.kind) ityp=1
      if(kind.eq.9  .and.  ipass.ge.kind) ityp=1
      if(kind.eq.12 .and.  ipass.ge.kind) ityp=1
c
c          do not scale midi, or mini inner shells
c
      scale=done
      if(omidi) go to 100
      if(ipass.le.(kind-nsplit)) go to 100
      if(ityp.eq.1) scale=scs
      if(ityp.eq.2) scale=scp
      if(ityp.eq.3) scale=scd
      if(scale.eq.dzero) scale=hscale(nucz)
100   itype = ityp + 15
      return
      end
      subroutine mini4(csinp,cpinp,cdinp,osplit,sc,nucz,intyp,nangm,
     +      nbfs,minf,maxf,loc,ngauss,ns,nshmax,ngsmax,ierr1,ierr2,
     +      nat,iwr)
c
c  +++ mini-4 basisset huzinaga
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
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     * ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     *cf(mxprim),cg(mxprim),
     +kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     +kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      common/blkin/eex(27),ccs(27),ccp(27),ccd(27)
      dimension csinp(*),cpinp(*),cdinp(*)
      dimension ns(*),intyp(*),nangm(*),nbfs(*),minf(*),maxf(*)
      data pt75/0.75d0/
      data dzero,pt5,pi32,tol /0.0d0,0.5d0,
     + 5.56832799683170d0,1.0d-10/
      data scalh1/1.00d0/
      call vclr(eex,1,108)
      mpass = 11
      if (nucz.gt.2) go to 120
      call gzero(eex,ccs,nucz)
      if (sc.le.dzero) sc = scalh1
      eex(1) = eex(1)*sc**2
      eex(2) = eex(2)*sc**2
      eex(3) = eex(3)*sc**2
      eex(4) = eex(4)*sc**2
      go to 170
120   if (nucz.gt.10) go to 130
      call gone(eex,ccs,ccp,nucz)
      go to 170
130   if (nucz.gt.18) go to 160
      call gtwo(eex,ccs,ccp,nucz)
      go to 170
160   if (nucz.gt.36) then
         if (opg_root()) then
            write(iwr,*)'*** nuclear charge = ',nucz
         endif
         call caserr2('requested basis set not available')
      endif
      call gthree(eex,ccs,ccp,ccd,nucz)
170   continue
      ipass = 0
180   ipass = ipass+1
      if (nucz.gt.2) go to 240
      if (osplit) go to 210
      go to (180,180,180,180,180,180,180,180,180,180,200),ipass
200   ityp = 16
      igauss = 4
      ig = 0
      go to 700
210   continue
      go to (180,180,180,180,180,180,180,180,180,220,230),ipass
220   ityp = 16
      igauss = 3
      ig = 0
      go to 700
230   ityp = 16
      igauss = 1
      ig = 3
      go to 700
240   if (nucz.gt.4) go to 310
      if (osplit) go to 270
      go to (180,180,180,180,180,180,180,180,180,200,260),ipass
260   ityp = 16
      igauss = 3
      ig = 4
      go to 700
270   continue
      go to (180,180,180,180,180,180,180,180,200,280,300),ipass
280   ityp = 16
      igauss = 2
      ig = 4
      go to 700
300   ityp = 16
      igauss = 1
      ig = 6
      go to 700
310   if (nucz.gt.10) go to 360
      if (osplit) go to 330
      go to (180,180,180,180,180,180,180,180,200,260,320),ipass
320   ityp = 17
      igauss = 4
      ig = 7
      go to 700
330   continue
      go to (180,180,180,180,180,180,200,280,300,340,350),ipass
340   ityp = 17
      igauss = 3
      ig = 7
      go to 700
350   ityp = 17
      igauss = 1
      ig = 10
      go to 700
360   if (nucz.gt.12) go to 440
      if (osplit) go to 390
      go to (180,180,180,180,180,180,180,200,260,370,380),ipass
370   ityp = 16
      igauss = 3
      ig = 7
      go to 700
380   ityp = 17
      igauss = 4
      ig = 10
      go to 700
390   continue
      go to (180,180,180,180,180,200,260,400,410,420,430),ipass
400   ityp = 16
      igauss = 2
      ig = 7
      go to 700
410   ityp = 16
      igauss = 1
      ig = 9
      go to 700
420   ityp = 17
      igauss = 3
      ig = 10
      go to 700
430   ityp = 17
      igauss = 1
      ig = 13
      go to 700
440   if (nucz.gt.18) go to 490
      if (osplit) go to 460
      go to (180,180,180,180,180,180,200,260,370,380,450),ipass
450   ityp = 17
      igauss = 3
      ig = 14
      go to 700
460   continue
      go to (180,180,180,180,200,260,400,410,380,470,480),ipass
470   ityp = 17
      igauss = 2
      ig = 14
      go to 700
480   ityp = 17
      igauss = 1
      ig = 16
      go to 700
490   if (nucz.gt.20) go to 590
      if (osplit) go to 530
      go to (180,180,180,180,180,200,260,370,500,510,520),ipass
500   ityp = 16
      igauss = 3
      ig = 10
      go to 700
510   ityp = 17
      igauss = 4
      ig = 13
      go to 700
520   ityp = 17
      igauss = 3
      ig = 17
      go to 700
530   continue
      go to (180,180,180,200,260,370,540,550,510,570,580),ipass
540   ityp = 16
      igauss = 2
      ig = 10
      go to 700
550   ityp = 16
      igauss = 1
      ig = 12
      go to 700
570   ityp = 17
      igauss = 2
      ig = 17
      go to 700
580   ityp = 17
      igauss = 1
      ig = 19
      go to 700
590   if (nucz.gt.30) go to 640
      if (osplit) go to 610
      go to (180,180,180,180,200,260,370,500,510,520,600),ipass
600   ityp = 18
      igauss = 4
      ig = 20
      go to 700
610   continue
      go to (180,200,260,370,540,550,510,570,580,620,630),ipass
620   ityp = 18
      igauss = 3
      ig = 20
      go to 700
630   ityp = 18
      igauss = 1
      ig = 23
      go to 700
640   if (osplit) go to 670
      go to (180,180,180,200,260,370,500,510,520,650,660),ipass
650   ityp = 17
      igauss = 3
      ig = 20
      go to 700
660   ityp = 18
      igauss = 4
      ig = 23
      go to 700
670   continue
      go to (200,260,370,540,550,510,520,680,685,690,695),ipass
680   ityp = 17
      igauss = 2
      ig = 20
      go to 700
685   ityp = 17
      igauss = 1
      ig = 22
      go to 700
690   ityp = 18
      igauss = 3
      ig = 23
      go to 700
695   ityp = 18
      igauss = 1
      ig = 26
700   continue
      nshell = nshell+1
      if (nshell.gt.nshmax) ierr1 = 1
      if (ierr1.ne.0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if (ngauss.gt.ngsmax) ierr2 = 1
      if (ierr2.ne.0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      do 720 i = 1,igauss
      k = k1+i-1
      ex(k) = eex(ig+i)
      csinp(k) = ccs(ig+i)
      cpinp(k) = ccp(ig+i)
      cdinp(k) = ccd(ig+i)
      cs(k) = csinp(k)
      cp(k) = cpinp(k)
720   cd(k) = cdinp(k)
      do 740 k = k1,k2
      ee = ex(k)+ex(k)
      facs = pi32/(ee*dsqrt(ee))
      facp = pt5*facs/ee
      facd = pt75*facs/(ee*ee)
      cs(k) = cs(k)/dsqrt(facs)
      cp(k) = cp(k)/dsqrt(facp)
740   cd(k) = cd(k)/dsqrt(facd)
      if (normf.eq.1) go to 820
      facs = dzero
      facp = dzero
      facd = dzero
      do 780 ig = k1,k2
      do 780 jg = k1,ig
      ee = ex(ig)+ex(jg)
      fac = ee*dsqrt(ee)
      dums = cs(ig)*cs(jg)/fac
      dump = pt5*cp(ig)*cp(jg)/(ee*fac)
      dumd = pt75*cd(ig)*cd(jg)/(ee**2*fac)
      if (ig.eq.jg) go to 760
      dums = dums+dums
      dump = dump+dump
      dumd = dumd+dumd
760   facs = facs+dums
      facp = facp+dump
      facd = facd+dumd
780   continue
      do 800 ig = k1,k2
      if (facs.gt.tol) cs(ig) = cs(ig)/dsqrt(facs*pi32)
      if (facp.gt.tol) cp(ig) = cp(ig)/dsqrt(facp*pi32)
800   if (facd.gt.tol) cd(ig) = cd(ig)/dsqrt(facd*pi32)
820   continue
      if (ipass.lt.mpass) go to 180
      return
      end
      subroutine gzero(e,s,n)
c
c  +++  mini basisset huzinaga
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*)
      go to (100,120),n
100   continue
      e(1) = 13.013372d0
      s(1) = 0.019678d0
      e(2) = 1.962496d0
      s(2) = 0.137952d0
      e(3) = 0.444569d0
      s(3) = 0.478313d0
      e(4) = 0.121953d0
      s(4) = 0.501131d0
      return
120   continue
      call caserr2('requested basis set not available for helium')
      return
      end
      subroutine gone(e,s,p,n)
c
c  +++ mini basisset huzinaga
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
100   continue
      e(1) = 104.72613d0
      s(1) = 0.02011d0
      e(2) = 15.66225d0
      s(2) = 0.13752d0
      e(3) = 3.39955d0
      s(3) = 0.45714d0
      e(4) = 0.85039d0
      s(4) = 0.53469d0
      e(5) = 0.75753d0
      s(5) = -0.09166d0
      e(6) = 0.06590d0
      s(6) = 0.64701d0
      e(7) = 0.02529d0
      s(7) = 0.41987d0
      return
120   continue
      e(1) = 193.43735d0
      s(1) = 0.01977d0
      e(2) = 29.03932d0
      s(2) = 0.13558d0
      e(3) = 6.38974d0
      s(3) = 0.45787d0
      e(4) = 1.64101d0
      s(4) = 0.53205d0
      e(5) = 1.90163d0
      s(5) = -0.08406d0
      e(6) = 0.15889d0
      s(6) = 0.58668d0
      e(7) = 0.05357d0
      s(7) = 0.48847d0
      return
140   continue
      e(1) = 308.88024d0
      s(1) = 0.01959d0
      e(2) = 46.44338d0
      s(2) = 0.13472d0
      e(3) = 10.30027d0
      s(3) = 0.45857d0
      e(4) = 2.68081d0
      s(4) = 0.53027d0
      e(5) = 3.44091d0
      s(5) = -0.08412d0
      e(6) = 0.29355d0
      s(6) = 0.58328d0
      e(7) = 0.09551d0
      s(7) = 0.49536d0
      e(8) = 6.03322d0
      p(8) = 0.03520d0
      e(9) = 1.24581d0
      p(9) = 0.19728d0
      e(10) = 0.33699d0
      p(10) = 0.50551d0
      e(11) = 0.09542d0
      p(11) = 0.48022d0
      return
160   continue
      e(1) = 451.90190d0
      s(1) = 0.01943d0
      e(2) = 67.99577d0
      s(2) = 0.13400d0
      e(3) = 15.14423d0
      s(3) = 0.45910d0
      e(4) = 3.96843d0
      s(4) = 0.52937d0
      e(5) = 5.40077d0
      s(5) = -0.08349d0
      e(6) = 0.46395d0
      s(6) = 0.57839d0
      e(7) = 0.14736d0
      s(7) = 0.50233d0
      e(8) = 9.47460d0
      p(8) = 0.03790d0
      e(9) = 2.00714d0
      p(9) = 0.20895d0
      e(10) = 0.54621d0
      p(10) = 0.50882d0
      e(11) = 0.15171d0
      p(11) = 0.46929d0
      return
180   continue
      e(1) = 621.71806d0
      s(1) = 0.01933d0
      e(2) = 93.58168d0
      s(2) = 0.13358d0
      e(3) = 20.90232d0
      s(3) = 0.45967d0
      e(4) = 5.50117d0
      s(4) = 0.52854d0
      e(5) = 7.78145d0
      s(5) = -0.08288d0
      e(6) = 0.67090d0
      s(6) = 0.57420d0
      e(7) = 0.20942d0
      s(7) = 0.50797d0
      e(8) = 13.58580d0
      p(8) = 0.03974d0
      e(9) = 2.92376d0
      p(9) = 0.21671d0
      e(10) = 0.79806d0
      p(10) = 0.51054d0
      e(11) = 0.21858d0
      p(11) = 0.46259d0
      return
200   continue
      e(1) = 821.83934d0
      s(1) = 0.01917d0
      e(2) = 123.68182d0
      s(2) = 0.13281d0
      e(3) = 27.66617d0
      s(3) = 0.45946d0
      e(4) = 7.29957d0
      s(4) = 0.52901d0
      e(5) = 10.60696d0
      s(5) = -0.08293d0
      e(6) = 0.91764d0
      s(6) = 0.57659d0
      e(7) = 0.28317d0
      s(7) = 0.50670d0
      e(8) = 17.75037d0
      p(8) = 0.04281d0
      e(9) = 3.86468d0
      p(9) = 0.22833d0
      e(10) = 1.04772d0
      p(10) = 0.50884d0
      e(11) = 0.27560d0
      p(11) = 0.46111d0
      return
220   continue
      e(1) = 1044.02406d0
      s(1) = 0.01916d0
      e(2) = 157.21281d0
      s(2) = 0.13282d0
      e(3) = 35.24628d0
      s(3) = 0.45997d0
      e(4) = 9.32771d0
      s(4) = 0.52807d0
      e(5) = 13.83395d0
      s(5) = -0.08284d0
      e(6) = 1.20268d0
      s(6) = 0.57583d0
      e(7) = 0.36776d0
      s(7) = 0.50835d0
      e(8) = 22.74336d0
      p(8) = 0.04466d0
      e(9) = 4.99043d0
      p(9) = 0.23521d0
      e(10) = 1.34908d0
      p(10) = 0.50864d0
      e(11) = 0.34749d0
      p(11) = 0.45864d0
      return
240   continue
      e(1) = 1295.58161d0
      s(1) = 0.01911d0
      e(2) = 195.10741d0
      s(2) = 0.13262d0
      e(3) = 43.79454d0
      s(3) = 0.46020d0
      e(4) = 11.60971d0
      s(4) = 0.52776d0
      e(5) = 17.49415d0
      s(5) = -0.08265d0
      e(6) = 1.52322d0
      s(6) = 0.57512d0
      e(7) = 0.46227d0
      s(7) = 0.50972d0
      e(8) = 28.44174d0
      p(8) = 0.04596d0
      e(9) = 6.27937d0
      p(9) = 0.23988d0
      e(10) = 1.69636d0
      p(10) = 0.50886d0
      e(11) = 0.43189d0
      p(11) = 0.45595d0
      return
      end
      subroutine gtwo(e,s,p,n)
c
c  +++  mini basisset huzinaga
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
100   continue
      e(1) = 1584.62613d0
      s(1) = 0.01896d0
      e(2) = 237.93226d0
      s(2) = 0.13204d0
      e(3) = 53.48784d0
      s(3) = 0.45899d0
      e(4) = 14.22088d0
      s(4) = 0.52910d0
      e(5) = 21.42052d0
      s(5) = -0.08586d0
      e(6) = 1.94477d0
      s(6) = 0.57234d0
      e(7) = 0.63758d0
      s(7) = 0.50724d0
      e(8) = 0.52440d0
      s(8) = -0.11339d0
      e(9) = 0.05407d0
      s(9) = 0.70084d0
      e(10) = 0.02106d0
      s(10) = 0.37434d0
      e(11) = 38.22544d0
      p(11) = 0.04370d0
      e(12) = 8.52713d0
      p(12) = 0.23443d0
      e(13) = 2.36602d0
      p(13) = 0.51289d0
      e(14) = 0.64731d0
      p(14) = 0.44280d0
      return
120   continue
      e(1) = 1878.63824d0
      s(1) = 0.01905d0
      e(2) = 284.11346d0
      s(2) = 0.13153d0
      e(3) = 64.18707d0
      s(3) = 0.45735d0
      e(4) = 17.11882d0
      s(4) = 0.53059d0
      e(5) = 25.86957d0
      s(5) = -0.08867d0
      e(6) = 2.40558d0
      s(6) = 0.58078d0
      e(7) = 0.82946d0
      s(7) = 0.49576d0
      e(8) = 0.90918d0
      s(8) = -0.12634d0
      e(9) = 0.10146d0
      s(9) = 0.64488d0
      e(10) = 0.03698d0
      s(10) = 0.44259d0
      e(11) = 49.11481d0
      p(11) = 0.04215d0
      e(12) = 11.09001d0
      p(12) = 0.22979d0
      e(13) = 3.13902d0
      p(13) = 0.51616d0
      e(14) = 0.89691d0
      p(14) = 0.43473d0
      return
140   continue
      e(1) = 2220.96184d0
      s(1) = 0.01894d0
      e(2) = 335.40618d0
      s(2) = 0.13107d0
      e(3) = 75.86390d0
      s(3) = 0.45638d0
      e(4) = 20.28413d0
      s(4) = 0.53163d0
      e(5) = 30.64994d0
      s(5) = -0.09129d0
      e(6) = 2.93625d0
      s(6) = 0.58299d0
      e(7) = 1.05426d0
      s(7) = 0.49157d0
      e(8) = 1.31578d0
      s(8) = -0.14841d0
      e(9) = 0.16269d0
      s(9) = 0.65887d0
      e(10) = 0.05946d0
      s(10) = 0.43887d0
      e(11) = 67.36412d0
      p(11) = 0.03572d0
      e(12) = 15.29472d0
      p(12) = 0.20630d0
      e(13) = 4.38339d0
      p(13) = 0.50534d0
      e(14) = 1.28709d0
      p(14) = 0.46254d0
      e(15) = 4.27052d0
      p(15) = -0.00855d0
      e(16) = 0.22965d0
      p(16) = 0.42253d0
      e(17) = 0.06459d0
      p(17) = 0.68105d0
      return
160   continue
      e(1) = 2581.38826d0
      s(1) = 0.01893d0
      e(2) = 390.50243d0
      s(2) = 0.13071d0
      e(3) = 88.51192d0
      s(3) = 0.45541d0
      e(4) = 23.72339d0
      s(4) = 0.53255d0
      e(5) = 35.78064d0
      s(5) = -0.09358d0
      e(6) = 3.51980d0
      s(6) =  0.58666d0
      e(7) = 1.30161d0
      s(7) = 0.48679d0
      e(8) = 1.75621d0
      s(8) = -0.16340d0
      e(9) = 0.22983d0
      s(9) = 0.66613d0
      e(10) = 0.08409d0
      s(10) = 0.43872d0
      e(11) = 82.63764d0
      p(11) = 0.03467d0
      e(12) = 18.80634d0
      p(12) = 0.20378d0
      e(13) = 5.45008d0
      p(13) = 0.50592d0
      e(14) = 1.64722d0
      p(14) = 0.45887d0
      e(15) = 5.64039d0
      p(15) = -0.01147d0
      e(16) = 0.33106d0
      p(16) = 0.46349d0
      e(17) = 0.09475d0
      p(17) = 0.64160d0
      return
180   continue
      e(1) = 2952.46717d0
      s(1) = 0.01905d0
      e(2) = 447.22377d0
      s(2) = 0.13114d0
      e(3) = 101.74040d0
      s(3) = 0.45512d0
      e(4) = 27.37811d0
      s(4) = 0.53198d0
      e(5) = 41.44401d0
      s(5) = -0.09549d0
      e(6) = 4.14795d0
      s(6) = 0.59214d0
      e(7) = 1.56898d0
      s(7) = 0.48049d0
      e(8) = 2.26237d0
      s(8) = -0.17484d0
      e(9) = 0.29976d0
      s(9) = 0.68195d0
      e(10) = 0.10950d0
      s(10) = 0.42719d0
      e(11) = 97.51908d0
      p(11) = 0.03471d0
      e(12) = 22.35582d0
      p(12) = 0.20420d0
      e(13) = 6.56167d0
      p(13) = 0.50704d0
      e(14) = 2.03529d0
      p(14) = 0.45230d0
      e(15) = 7.37347d0
      p(15) = -0.01337d0
      e(16) = 0.45273d0
      p(16) = 0.47940d0
      e(17) = 0.13089d0
      p(17) = 0.62553d0
      return
200   continue
      e(1) = 3423.41131d0
      s(1) = 0.01867d0
      e(2) = 516.87666d0
      s(2) = 0.12937d0
      e(3) = 117.31748d0
      s(3) = 0.45280d0
      e(4) = 31.52309d0
      s(4) = 0.53591d0
      e(5) = 47.63257d0
      s(5) = -0.09656d0
      e(6) = 4.79276d0
      s(6) = 0.59627d0
      e(7) = 1.83831d0
      s(7) = 0.47574d0
      e(8) = 2.77567d0
      s(8) = -0.18823d0
      e(9) = 0.38467d0
      s(9) = 0.67301d0
      e(10) = 0.14030d0
      s(10) = 0.44277d0
      e(11) = 115.12771d0
      p(11) = 0.03410d0
      e(12) = 26.46023d0
      p(12) = 0.20292d0
      e(13) = 7.81884d0
      p(13) = 0.50851d0
      e(14) = 2.46658d0
      p(14) = 0.44875d0
      e(15) = 9.30690d0
      p(15) = -0.01522d0
      e(16) = 0.59482d0
      p(16) = 0.48462d0
      e(17) = 0.16512d0
      p(17) = 0.62662d0
      return
220   continue
      e(1) = 3855.56170d0
      s(1) = 0.01873d0
      e(2) = 584.23357d0
      s(2) = 0.12920d0
      e(3) = 132.93186d0
      s(3) = 0.45217d0
      e(4) = 35.79090d0
      s(4) = 0.53636d0
      e(5) = 54.15434d0
      s(5) = -0.09844d0
      e(6) = 5.52527d0
      s(6) = 0.60937d0
      e(7) = 2.15040d0
      s(7) = 0.46242d0
      e(8) = 3.34449d0
      s(8) = -0.19079d0
      e(9) = 0.47494d0
      s(9) = 0.68441d0
      e(10) = 0.17308d0
      s(10) = 0.43372d0
      e(11) = 131.59905d0
      p(11) = 0.03449d0
      e(12) = 30.40899d0
      p(12) = 0.20499d0
      e(13) = 9.05879d0
      p(13) = 0.51141d0
      e(14) = 2.90681d0
      p(14) = 0.44052d0
      e(15) = 11.28714d0
      p(15) = -0.01652d0
      e(16) = 0.75025d0
      p(16) = 0.48988d0
      e(17) = 0.20510d0
      p(17) = 0.62406d0
      return
240   continue
      e(1) = 4379.83516d0
      s(1) = 0.01846d0
      e(2) = 663.38769d0
      s(2) = 0.12767d0
      e(3) = 150.74592d0
      s(3) = 0.45013d0
      e(4) = 40.52929d0
      s(4) = 0.53976d0
      e(5) = 61.25439d0
      s(5) = -0.09961d0
      e(6) = 6.28070d0
      s(6) = 0.61793d0
      e(7) = 2.46782d0
      s(7) = 0.45359d0
      e(8) = 3.95647d0
      s(8) = -0.19594d0
      e(9) = 0.57301d0
      s(9) = 0.68535d0
      e(10) = 0.20874d0
      s(10) = 0.43561d0
      e(11) = 152.61138d0
      p(11) = 0.03368d0
      e(12) = 35.31198d0
      p(12) = 0.20282d0
      e(13) = 10.55209d0
      p(13) = 0.51225d0
      e(14) = 3.41986d0
      p(14) = 0.43980d0
      e(15) = 13.40059d0
      p(15) = -0.01774d0
      e(16) = 0.91883d0
      p(16) = 0.49558d0
      e(17) = 0.25042d0
      p(17) = 0.61931d0
      return
      end
      subroutine gthree(e,s,p,d,n)
c
c  +++ mini basisset huzinaga
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-18
      go to (100,120,130,140,150,160,170,180,190,200,
     +210,220,230,240,250,260,270,280),nn
100   continue
      e(1) = 5012.49720d0
      s(1) = 0.01792d0
      e(2) = 755.53034d0
      s(2) = 0.12517d0
      e(3) = 170.85707d0
      s(3) = 0.44723d0
      e(4) = 45.75320d0
      s(4) = 0.54515d0
      e(5) = 68.48405d0
      s(5) = -0.10116d0
      e(6) = 7.19166d0
      s(6) = 0.60909d0
      e(7) = 2.87354d0
      s(7) = 0.46229d0
      e(8) = 4.56739d0
      s(8) = -0.20912d0
      e(9) = 0.71123d0
      s(9) = 0.67448d0
      e(10) = 0.28239d0
      s(10) = 0.44910d0
      e(11) = 0.23828d0
      s(11) = 0.14941d0
      e(12) = 0.03646d0
      s(12) = -0.69630d0
      e(13) = 0.01576d0
      s(13) = -0.39957d0
      e(14) = 176.05906d0
      p(14) = 0.03275d0
      e(15) = 40.87042d0
      p(15) = 0.19887d0
      e(16) = 12.28523d0
      p(16) = 0.51003d0
      e(17) = 4.01859d0
      p(17) = 0.44382d0
      e(18) = 15.14810d0
      p(18) = -0.01784d0
      e(19) = 1.16433d0
      p(19) = 0.49478d0
      e(20) = 0.34300d0
      p(20) = 0.60921d0
      return
120   continue
      e(1) = 5622.34128d0
      s(1) = 0.01768d0
      e(2) = 846.96003d0
      s(2) = 0.12380d0
      e(3) = 191.31099d0
      s(3) = 0.44524d0
      e(4) = 51.17270d0
      s(4) = 0.54832d0
      e(5) = 76.52458d0
      s(5) = -0.10206d0
      e(6) = 8.04094d0
      s(6) = 0.61867d0
      e(7) = 3.23797d0
      s(7) = 0.45242d0
      e(8) = 5.26751d0
      s(8) = -0.21989d0
      e(9) = 0.84626d0
      s(9) = 0.69444d0
      e(10) = 0.35254d0
      s(10) = 0.43216d0
      e(11) = 0.41520d0
      s(11) = -0.17025d0
      e(12) = 0.06157d0
      s(12) = 0.66893d0
      e(13) = 0.02483d0
      s(13) = 0.43748d0
      e(14) = 202.98255d0
      p(14) = 0.03153d0
      e(15) = 47.10648d0
      p(15) = 0.19460d0
      e(16) = 14.18286d0
      p(16) = 0.50873d0
      e(17) = 4.66623d0
      p(17) = 0.44796d0
      e(18) = 17.08207d0
      p(18) = -0.01837d0
      e(19) = 1.42216d0
      p(19) = 0.49898d0
      e(20) = 0.44121d0
      p(20) = 0.59816d0
      return
130   continue
      e(1) = 6301.57893d0
      s(1) = 0.01737d0
      e(2) = 946.20216d0
      s(2) = 0.12239d0
      e(3) = 213.08479d0
      s(3) = 0.44367d0
      e(4) = 56.88553d0
      s(4) = 0.55126d0
      e(5) = 84.60882d0
      s(5) = -0.10340d0
      e(6) = 9.05720d0
      s(6) = 0.61508d0
      e(7) = 3.69065d0
      s(7) = 0.45607d0
      e(8) = 6.02937d0
      s(8) = -0.22308d0
      e(9) = 0.98035d0
      s(9) = 0.70087d0
      e(10) = 0.40764d0
      s(10) = 0.42791d0
      e(11) = 0.49191d0
      s(11) = -0.16090d0
      e(12) = 0.06983d0
      s(12) = 0.66617d0
      e(13) = 0.02747d0
      s(13) = 0.43651d0
      e(14) = 239.43133d0
      p(14) = 0.02905d0
      e(15) = 55.29443d0
      p(15) = 0.18529d0
      e(16) = 16.57390d0
      p(16) = 0.50471d0
      e(17) = 5.43466d0
      p(17) = 0.46109d0
      e(18) = 19.52532d0
      p(18) = -0.01890d0
      e(19) = 1.66790d0
      p(19) = 0.50154d0
      e(20) = 0.51969d0
      p(20) = 0.59524d0
      e(21) = 11.56957d0
      d(21) = 0.05926d0
      e(22) = 2.91387d0
      d(22) = 0.26363d0
      e(23) = 0.84235d0
      d(23) = 0.51110d0
      e(24) = 0.22647d0
      d(24) = 0.48303d0
      return
140   continue
      e(1) = 6878.30916d0
      s(1) = 0.01749d0
      e(2) = 1035.14430d0
      s(2) = 0.12278d0
      e(3) = 233.75774d0
      s(3) = 0.44363d0
      e(4) = 62.57798d0
      s(4) = 0.55062d0
      e(5) = 93.42842d0
      s(5) = -0.10405d0
      e(6) = 10.01838d0
      s(6) = 0.62160d0
      e(7) = 4.09094d0
      s(7) = 0.44963d0
      e(8) = 6.83736d0
      s(8) = -0.22536d0
      e(9) = 1.11487d0
      s(9) = 0.70431d0
      e(10) = 0.45976d0
      s(10) = 0.42603d0
      e(11) = 0.55462d0
      s(11) = -0.15456d0
      e(12) = 0.07749d0
      s(12) = 0.64855d0
      e(13) = 0.02996d0
      s(13) = 0.45234d0
      e(14) = 261.68259d0
      p(14) = 0.02975d0
      e(15) = 60.68987d0
      p(15) = 0.18847d0
      e(16) = 18.30848d0
      p(16) = 0.50773d0
      e(17) = 6.06318d0
      p(17) = 0.45318d0
      e(18) = 22.23974d0
      p(18) = -0.01978d0
      e(19) = 1.91392d0
      p(19) = 0.50772d0
      e(20) = 0.59525d0
      p(20) = 0.58968d0
      e(21) = 13.69844d0
      d(21) = 0.06184d0
      e(22) = 3.51187d0
      d(22) = 0.27274d0
      e(23) = 1.03911d0
      d(23) = 0.51596d0
      e(24) = 0.28709d0
      d(24) = 0.46068d0
      return
150   continue
      e(1) = 7538.50148d0
      s(1) = 0.01745d0
      e(2) = 1134.73776d0
      s(2) = 0.12248d0
      e(3) = 256.39486d0
      s(3) = 0.44299d0
      e(4) = 68.69291d0
      s(4) = 0.55138d0
      e(5) = 102.70394d0
      s(5) = -0.10485d0
      e(6) = 11.05970d0
      s(6) = 0.62650d0
      e(7) = 4.53443d0
      s(7) = 0.44479d0
      e(8) = 7.71280d0
      s(8) = -0.22567d0
      e(9) = 1.24857d0
      s(9) = 0.71482d0
      e(10) = 0.50982d0
      s(10) = 0.41588d0
      e(11) = 0.61771d0
      s(11) = -0.14698d0
      e(12) = 0.08441d0
      s(12) = 0.64017d0
      e(13) = 0.03210d0
      s(13) = 0.45833d0
      e(14) = 293.05988d0
      p(14) = 0.02913d0
      e(15) = 67.96779d0
      p(15) = 0.18637d0
      e(16) = 20.52353d0
      p(16) = 0.50768d0
      e(17) = 6.81398d0
      p(17) = 0.45457d0
      e(18) = 25.00415d0
      p(18) = -0.02030d0
      e(19) = 2.18110d0
      p(19) = 0.50948d0
      e(20) = 0.67669d0
      p(20) = 0.58843d0
      e(21) = 16.00851d0
      d(21) = 0.06267d0
      e(22) = 4.14854d0
      d(22) = 0.27667d0
      e(23) = 1.23972d0
      d(23) = 0.51736d0
      e(24) = 0.34410d0
      d(24) = 0.45302d0
      return
160   continue
      e(1) = 8259.66169d0
      s(1) = 0.01734d0
      e(2) = 1242.65581d0
      s(2) = 0.12186d0
      e(3) = 280.62986d0
      s(3) = 0.44213d0
      e(4) = 75.16762d0
      s(4) = 0.55274d0
      e(5) = 112.31259d0
      s(5) = -0.10552d0
      e(6) = 12.16673d0
      s(6) = 0.62906d0
      e(7) = 4.99791d0
      s(7) = 0.44245d0
      e(8) = 8.61352d0
      s(8) = -0.22380d0
      e(9) = 1.39148d0
      s(9) = 0.70870d0
      e(10) = 0.54903d0
      s(10) = 0.42334d0
      e(11) = 0.67041d0
      s(11) = -0.12211d0
      e(12) = 0.07860d0
      s(12) = 0.61508d0
      e(13) = 0.02944d0
      s(13) = 0.47068d0
      e(14) = 326.34485d0
      p(14) = 0.02853d0
      e(15) = 75.71993d0
      p(15) = 0.18440d0
      e(16) = 22.86154d0
      p(16) = 0.50804d0
      e(17) = 7.59831d0
      p(17) = 0.45568d0
      e(18) = 28.77843d0
      p(18) = -0.02103d0
      e(19) = 2.42253d0
      p(19) = 0.51611d0
      e(20) = 0.73195d0
      p(20) = 0.58574d0
      e(21) = 16.96163d0
      d(21) = 0.06604d0
      e(22) = 4.36311d0
      d(22) = 0.28517d0
      e(23) = 1.26148d0
      d(23) = 0.51091d0
      e(24) = 0.31945d0
      d(24) = 0.47400d0
      return
170   continue
      e(1) = 8952.97049d0
      s(1) = 0.01736d0
      e(2) = 1349.07130d0
      s(2) = 0.12179d0
      e(3) = 305.10198d0
      s(3) = 0.44165d0
      e(4) = 81.82873d0
      s(4) = 0.55308d0
      e(5) = 122.42058d0
      s(5) = -0.10616d0
      e(6) = 13.32666d0
      s(6) = 0.63120d0
      e(7) = 5.49330d0
      s(7) = 0.44040d0
      e(8) = 9.50157d0
      s(8) = -0.22770d0
      e(9) = 1.58022d0
      s(9) = 0.69368d0
      e(10) = 0.64232d0
      s(10) = 0.43963d0
      e(11) = 0.76522d0
      s(11) = -0.13801d0
      e(12) = 0.09568d0
      s(12) = 0.63778d0
      e(13) = 0.03543d0
      s(13) = 0.45618d0
      e(14) = 367.14467d0
      p(14) = 0.02732d0
      e(15) = 85.22744d0
      p(15) = 0.17889d0
      e(16) = 25.72239d0
      p(16) = 0.50478d0
      e(17) = 8.53416d0
      p(17) = 0.46401d0
      e(18) = 31.03165d0
      p(18) = -0.02106d0
      e(19) = 2.76359d0
      p(19) = 0.51018d0
      e(20) = 0.85153d0
      p(20) = 0.58890d0
      e(21) = 20.91329d0
      d(21) = 0.06381d0
      e(22) = 5.50492d0
      d(22) = 0.28215d0
      e(23) = 1.66445d0
      d(23) = 0.51856d0
      e(24) = 0.46258d0
      d(24) = 0.44430d0
      return
180   continue
      e(1) = 9678.75072d0
      s(1) = 0.01739d0
      e(2) = 1458.96956d0
      s(2) = 0.12185d0
      e(3) = 330.32731d0
      s(3) = 0.44148d0
      e(4) = 88.71446d0
      s(4) = 0.55303d0
      e(5) = 132.78442d0
      s(5) = -0.10679d0
      e(6) = 14.54392d0
      s(6) = 0.63381d0
      e(7) = 6.00098d0
      s(7) = 0.43808d0
      e(8) = 10.46384d0
      s(8) = 0.22856d0
      e(9) = 1.75242d0
      s(9) = -0.69246d0
      e(10) = 0.70942d0
      s(10) = -0.44189d0
      e(11) = 0.85008d0
      s(11) = -0.13378d0
      e(12) = 0.10263d0
      s(12) = 0.63617d0
      e(13) = 0.03741d0
      s(13) = 0.45635d0
      e(14) = 395.28902d0
      p(14) = 0.02785d0
      e(15) = 91.98856d0
      p(15) = 0.18152d0
      e(16) = 27.88387d0
      p(16) = 0.50738d0
      e(17) = 9.31097d0
      p(17) = 0.45777d0
      e(18) = 34.42680d0
      p(18) = -0.02151d0
      e(19) = 3.07012d0
      p(19) = 0.51350d0
      e(20) = 0.94469d0
      p(20) = 0.58595d0
      e(21) = 23.15678d0
      d(21) = 0.06553d0
      e(22) = 6.12648d0
      d(22) = 0.28743d0
      e(23) = 1.84895d0
      d(23) = 0.51754d0
      e(24) = 0.50594d0
      d(24) = 0.44332d0
      return
190   continue
      e(1) = 10425.41800d0
      s(1) = 0.01742d0
      e(2) = 1572.52255d0
      s(2) = 0.12195d0
      e(3) = 356.47734d0
      s(3) = 0.44133d0
      e(4) = 95.87448d0
      s(4) = 0.55292d0
      e(5) = 143.78973d0
      s(5) = -0.10728d0
      e(6) = 15.78481d0
      s(6) = 0.63800d0
      e(7) = 6.51752d0
      s(7) = 0.43403d0
      e(8) = 11.50024d0
      s(8) = -0.22855d0
      e(9) = 1.92358d0
      s(9) = 0.69452d0
      e(10) = 0.77475d0
      s(10) = 0.44010d0
      e(11) = 0.93055d0
      s(11) = -0.12976d0
      e(12) = 0.10971d0
      s(12) = 0.62967d0
      e(13) = 0.03949d0
      s(13) = 0.46174d0
      e(14) = 424.82412d0
      p(14) = 0.02832d0
      e(15) = 99.03732d0
      p(15) = 0.18394d0
      e(16) = 30.13287d0
      p(16) = 0.50983d0
      e(17) = 10.11710d0
      p(17) = 0.45202d0
      e(18) = 37.92911d0
      p(18) = -0.02185d0
      e(19) = 3.39318d0
      p(19) = 0.51561d0
      e(20) = 1.04180d0
      p(20) = 0.58427d0
      e(21) = 25.54420d0
      d(21) = 0.06700d0
      e(22) = 6.79072d0
      d(22) = 0.29192d0
      e(23) = 2.04992d0
      d(23) = 0.51711d0
      e(24) = 0.55639d0
      d(24) = 0.44067d0
      return
200   continue
      e(1) = 11199.40300d0
      s(1) = 0.01745d0
      e(2) = 1690.16398d0
      s(2) = 0.12206d0
      e(3) = 383.57266d0
      s(3) = 0.44123d0
      e(4) = 103.30348d0
      s(4) = 0.55276d0
      e(5) = 155.22941d0
      s(5) = -0.10773d0
      e(6) = 17.07780d0
      s(6) = 0.64183d0
      e(7) = 7.05392d0
      s(7) = 0.43035d0
      e(8) = 12.57705d0
      s(8) = -0.22849d0
      e(9) = 2.10314d0
      s(9) = 0.69550d0
      e(10) = 0.84301d0
      s(10) = 0.43942d0
      e(11) = 1.01244d0
      s(11) = -0.12606d0
      e(12) = 0.11682d0
      s(12) = 0.62268d0
      e(13) = 0.04156d0
      s(13) = 0.46779d0
      e(14) = 455.74189d0
      p(14) = 0.02872d0
      e(15) = 106.48141d0
      p(15) = 0.18588d0
      e(16) = 32.50784d0
      p(16) = 0.51183d0
      e(17) = 10.96437d0
      p(17) = 0.44727d0
      e(18) = 41.64391d0
      p(18) = -0.02221d0
      e(19) = 3.72798d0
      p(19) = 0.51779d0
      e(20) = 1.14212d0
      p(20) = 0.58251d0
      e(21) = 28.03443d0
      d(21) = 0.06831d0
      e(22) = 7.48866d0
      d(22) = 0.29567d0
      e(23) = 2.26217d0
      d(23) = 0.51672d0
      e(24) = 0.61010d0
      d(24) = 0.43825d0
      return
210   continue
      e(1) = 11901.91700d0
      s(1) = 0.01764d0
      e(2) = 1802.84331d0
      s(2) = 0.12261d0
      e(3) = 410.72130d0
      s(3) = 0.44119d0
      e(4) = 110.91181d0
      s(4) = 0.55197d0
      e(5) = 167.83708d0
      s(5) = -0.10778d0
      e(6) = 18.24447d0
      s(6) = 0.65712d0
      e(7) = 7.48671d0
      s(7) = 0.41505d0
      e(8) = 13.81524d0
      s(8) = -0.22432d0
      e(9) = 2.22960d0
      s(9) = 0.71868d0
      e(10) = 0.85996d0
      s(10) = 0.41511d0
      e(11) = 0.91788d0
      s(11) = -0.09319d0
      e(12) = 0.09217d0
      s(12) = 0.57582d0
      e(13) = 0.03316d0
      s(13) = 0.50010d0
      e(14) = 485.34916d0
      p(14) = 0.02930d0
      e(15) = 113.70531d0
      p(15) = 0.18869d0
      e(16) = 34.82864d0
      p(16) = 0.51465d0
      e(17) = 11.79948d0
      p(17) = 0.44088d0
      e(18) = 45.94963d0
      p(18) = -0.02301d0
      e(19) = 4.02529d0
      p(19) = 0.52484d0
      e(20) = 1.20961d0
      p(20) = 0.57843d0
      e(21) = 28.72223d0
      d(21) = 0.07337d0
      e(22) = 7.63462d0
      d(22) = 0.30693d0
      e(23) = 2.24888d0
      d(23) = 0.51187d0
      e(24) = 0.56497d0
      d(24) = 0.44890d0
      return
220   continue
      e(1) = 12838.19400d0
      s(1) = 0.01750d0
      e(2) = 1939.18903d0
      s(2) = 0.12219d0
      e(3) = 440.92411d0
      s(3) = 0.44096d0
      e(4) = 119.02647d0
      s(4) = 0.55263d0
      e(5) = 179.46791d0
      s(5) = -0.10851d0
      e(6) = 19.81290d0
      s(6) = 0.64868d0
      e(7) = 8.18921d0
      s(7) = 0.42375d0
      e(8) = 14.93354d0
      s(8) = 0.22711d0
      e(9) = 2.45812d0
      s(9) = -0.70658d0
      e(10) = 0.97304d0
      s(10) = -0.42771d0
      e(11) = 1.18291d0
      s(11) = -0.11915d0
      e(12) = 0.12953d0
      s(12) = 0.61710d0
      e(13) = 0.04504d0
      s(13) = 0.47147d0
      e(14) = 519.11459d0
      p(14) = 0.02959d0
      e(15) = 121.83170d0
      p(15) = 0.18998d0
      e(16) = 37.43978d0
      p(16) = 0.51577d0
      e(17) = 12.73124d0
      p(17) = 0.43775d0
      e(18) = 49.41241d0
      p(18) = -0.02278d0
      e(19) = 4.44378d0
      p(19) = 0.52079d0
      e(20) = 1.35486d0
      p(20) = 0.58038d0
      e(21) = 33.68740d0
      d(21) = 0.06949d0
      e(22) = 9.06147d0
      d(22) = 0.29985d0
      e(23) = 2.74053d0
      d(23) = 0.51635d0
      e(24) = 0.73211d0
      d(24) = 0.43585d0
      return
230   continue
      e(1) = 13931.937d0
      s(1) = 0.0171814d0
      e(2) = 2100.4471d0
      s(2) = 0.1205193d0
      e(3) = 476.12690d0
      s(3) = 0.4388564d0
      e(4) = 128.12429d0
      s(4) = 0.5563888d0
      e(5) = 192.22228d0
      s(5) = -0.1089698d0
      e(6) = 21.347910d0
      s(6) = 0.6448008d0
      e(7) = 8.8805409d0
      s(7) = 0.4275402d0
      e(8) = 16.184874d0
      s(8) = 0.2286846d0
      e(9) = 2.6512096d0
      s(9) = -0.7181421d0
      e(10) = 1.0674061d0
      s(10) = -0.4154128d0
      e(11) = 1.3554783d0
      s(11) = -0.1482103d0
      e(12) = 0.18472212d0
      s(12) = 0.6387621d0
      e(13) = 0.06639980d0
      s(13) = 0.4642895d0
      e(14) = 579.86288d0
      p(14) = 0.0279630d0
      e(15) = 135.76337d0
      p(15) = 0.1813498d0
      e(16) = 41.889217d0
      p(16) = 0.5067964d0
      e(17) = 14.174081d0
      p(17) = 0.4552253d0
      e(18) = 6.3268159d0
      p(18) = 0.3221193d0
      e(19) = 2.4320543d0
      p(19) = 0.5553634d0
      e(20) = 0.93392650d0
      p(20) = 0.2397043d0
      e(21) = 0.30045657d0
      p(21) = 0.2463414d0
      e(22) = 0.10605282d0
      p(22) = 0.5505357d0
      e(23) = 0.03842867d0
      p(23) = 0.3376249d0
      e(24) = 38.593817d0
      d(24) = 0.0671457d0
      e(25) = 10.487668d0
      d(25) = 0.2959598d0
      e(26) = 3.2393388d0
      d(26) = 0.5205164d0
      e(27) = 0.90878475d0
      d(27) = 0.4227369d0
      return
240   continue
      e(1) = 14881.232d0
      s(1) = 0.0171408d0
      e(2) = 2243.5833d0
      s(2) = 0.1202577d0
      e(3) = 508.65980d0
      s(3) = 0.4383460d0
      e(4) = 136.91849d0
      s(4) = 0.5570547d0
      e(5) = 205.38629d0
      s(5) = -0.1093452d0
      e(6) = 22.881048d0
      s(6) = 0.6454654d0
      e(7) = 9.5419406d0
      s(7) = 0.4269113d0
      e(8) = 17.442959d0
      s(8) = -0.2310171d0
      e(9) = 2.8747550d0
      s(9) = 0.7219613d0
      e(10) = 1.1811841d0
      s(10) = 0.4116321d0
      e(11) = 1.5871377d0
      s(11) = -0.1721057d0
      e(12) = 0.23235290d0
      s(12) = 0.6680776d0
      e(13) = 0.08493708d0
      s(13) = 0.4452923d0
      e(14) = 623.84392d0
      p(14) = 0.0277654d0
      e(15) = 146.15819d0
      p(15) = 0.1805044d0
      e(16) = 45.169851d0
      p(16) = 0.5061559d0
      e(17) = 15.320233d0
      p(17) = 0.4561635d0
      e(18) = 6.8683736d0
      p(18) = 0.3244056d0
      e(19) = 2.6634716d0
      p(19) = 0.5556493d0
      e(20) = 1.0480038d0
      p(20) = 0.2332192d0
      e(21) = 0.05234136d0
      p(21) = 0.2775625d0
      e(22) = 0.37181803d0
      p(22) = 0.2941707d0
      e(23) = 0.13927348d0
      p(23) = 0.5528443d0
      e(24) = 44.106656d0
      d(24) = 0.0642230d0
      e(25) = 12.077536d0
      d(25) = 0.2903450d0
      e(26) = 3.7899474d0
      d(26) = 0.5237000d0
      e(27) = 1.1019081d0
      d(27) = 0.4159562d0
      return
250   continue
      e(1) = 15860.364d0
      s(1) = 0.0171048d0
      e(2) = 2391.2232d0
      s(2) = 0.1200249d0
      e(3) = 542.23120d0
      s(3) = 0.4378747d0
      e(4) = 146.00269d0
      s(4) = 0.5576566d0
      e(5) = 219.08007d0
      s(5) = -0.1096535d0
      e(6) = 24.437631d0
      s(6) = 0.6475717d0
      e(7) = 10.206970d0
      s(7) = 0.4248194d0
      e(8) = 18.744848d0
      s(8) = -0.2334702d0
      e(9) = 3.1047892d0
      s(9) = 0.7294911d0
      e(10) = 1.2970598d0
      s(10) = 0.4044278d0
      e(11) = 1.8301844d0
      s(11) = -0.1924620d0
      e(12) = 0.28147770d0
      s(12) = 0.6866084d0
      e(13) = 0.10437810d0
      s(13) = 0.4355541d0
      e(14) = 670.45052d0
      p(14) = 0.0275156d0
      e(15) = 157.17065d0
      p(15) = 0.1794001d0
      e(16) = 48.639537d0
      p(16) = 0.5052466d0
      e(17) = 16.532308d0
      p(17) = 0.4576592d0
      e(18) = 7.4556466d0
      p(18) = 0.3248901d0
      e(19) = 2.9185276d0
      p(19) = 0.5568808d0
      e(20) = 1.1701325d0
      p(20) = 0.2280292d0
      e(21) = 0.17189969d0
      p(21) = 0.5495546d0
      e(22) = 0.44648872d0
      p(22) = 0.3297273d0
      e(23) = 0.06621832d0
      p(23) = 0.2383460d0
      e(24) = 49.802186d0
      d(24) = 0.0618851d0
      e(25) = 13.722924d0
      d(25) = 0.2858504d0
      e(26) = 4.3622725d0
      d(26) = 0.5264330d0
      e(27) = 1.3044304d0
      d(27) = 0.4103582d0
      return
260   continue
      e(1) = 16870.587d0
      s(1) = 0.0170711d0
      e(2) = 2543.6306d0
      s(2) = 0.1197991d0
      e(3) = 576.90712d0
      s(3) = 0.4374221d0
      e(4) = 155.38134d0
      s(4) = 0.5582388d0
      e(5) = 232.94274d0
      s(5) = -0.1100492d0
      e(6) = 26.119010d0
      s(6) = 0.6460109d0
      e(7) = 10.950779d0
      s(7) = 0.4264111d0
      e(8) = 20.093365d0
      s(8) = -0.2361244d0
      e(9) = 3.3498536d0
      s(9) = 0.7360716d0
      e(10) = 1.4214673d0
      s(10) = 0.3984420d0
      e(11) = 2.0742664d0
      s(11) = -0.2101056d0
      e(12) = 0.33403502d0
      s(12) = 0.7060172d0
      e(13) = 0.12457287d0
      s(13) = 0.4249450d0
      e(14) = 717.79796d0
      p(14) = 0.0273385d0
      e(15) = 168.39415d0
      p(15) = 0.1785904d0
      e(16) = 52.192242d0
      p(16) = 0.5045380d0
      e(17) = 17.779509d0
      p(17) = 0.4586622d0
      e(18) = 8.0417034d0
      p(18) = 0.3281021d0
      e(19) = 3.1683176d0
      p(19) = 0.5593898d0
      e(20) = 1.2871289d0
      p(20) = 0.2193119d0
      e(21) = 0.53499306d0
      p(21) = 0.3487730d0
      e(22) = 0.20451022d0
      p(22) = 0.5359349d0
      e(23) = 0.07749696d0
      p(23) = 0.2376894d0
      e(24) = 55.836292d0
      d(24) = 0.0597210d0
      e(25) = 15.468666d0
      d(25) = 0.2815170d0
      e(26) = 4.9700412d0
      d(26) = 0.5287895d0
      e(27) = 1.5203984d0
      d(27) = 0.4063501d0
      return
270   continue
      e(1) = 17914.012d0
      s(1) = 0.0170373d0
      e(2) = 2700.8820d0
      s(2) = 0.1195879d0
      e(3) = 612.63764d0
      s(3) = 0.4370146d0
      e(4) = 165.04835d0
      s(4) = 0.5587696d0
      e(5) = 247.33791d0
      s(5) = -0.1103867d0
      e(6) = 27.825882d0
      s(6) = 0.6460753d0
      e(7) = 11.692350d0
      s(7) = 0.4264025d0
      e(8) = 21.472689d0
      s(8) = -0.2389743d0
      e(9) = 3.6072844d0
      s(9) = 0.7433529d0
      e(10) = 1.5510670d0
      s(10) = 0.3921021d0
      e(11) = 2.3397075d0
      s(11) = -0.2266599d0
      e(12) = 0.38832539d0
      s(12) = 0.7202304d0
      e(13) = 0.14562498d0
      s(13) = 0.4183848d0
      e(14) = 769.88290d0
      p(14) = 0.0270082d0
      e(15) = 180.60286d0
      p(15) = 0.1771929d0
      e(16) = 56.012187d0
      p(16) = 0.5034818d0
      e(17) = 19.108991d0
      p(17) = 0.4607090d0
      e(18) = 8.6662920d0
      p(18) = 0.3301375d0
      e(19) = 3.4372753d0
      p(19) = 0.5621199d0
      e(20) = 1.4102279d0
      p(20) = 0.2118896d0
      e(21) = 0.62509542d0
      p(21) = 0.3662229d0
      e(22) = 0.23852551d0
      p(22) = 0.5323300d0
      e(23) = 0.08924869d0
      p(23) = 0.2245335d0
      e(24) = 61.699454d0
      d(24) = 0.0584980d0
      e(25) = 17.176489d0
      d(25) = 0.2795209d0
      e(26) = 5.5706174d0
      d(26) = 0.5311305d0
      e(27) = 1.7383942d0
      d(27) = 0.4006201d0
      return
280   continue
      e(1) = 18990.321d0
      s(1) = 0.0170036d0
      e(2) = 2863.2884d0
      s(2) = 0.1193618d0
      e(3) = 649.58505d0
      s(3) = 0.4365613d0
      e(4) = 175.03786d0
      s(4) = 0.5593600d0
      e(5) = 262.24460d0
      s(5) = -0.1106834d0
      e(6) = 29.573679d0
      s(6) = 0.6464273d0
      e(7) = 12.459722d0
      s(7) = 0.4260356d0
      e(8) = 22.925626d0
      s(8) = -0.2414622d0
      e(9) = 3.8683898d0
      s(9) = 0.7529089d0
      e(10) = 1.6809618d0
      s(10) = 0.3833081d0
      e(11) = 2.6079105d0
      s(11) = -0.2392303d0
      e(12) = 0.44479912d0
      s(12) = 0.7301479d0
      e(13) = 0.16799905d0
      s(13) = 0.4145776d0
      e(14) = 820.11422d0
      p(14) = 0.0268826d0
      e(15) = 192.54933d0
      p(15) = 0.1765619d0
      e(16) = 59.813230d0
      p(16) = 0.5028124d0
      e(17) = 20.451806d0
      p(17) = 0.4615069d0
      e(18) = 9.3219035d0
      p(18) = 0.3311094d0
      e(19) = 3.7297257d0
      p(19) = 0.5638281d0
      e(20) = 1.5422929d0
      p(20) = 0.2068282d0
      e(21) = 0.72202413d0
      p(21) = 0.3764629d0
      e(22) = 0.27801808d0
      p(22) = 0.5270024d0
      e(23) = 0.10446752d0
      p(23) = 0.2180729d0
      e(24) = 68.889745d0
      d(24) = 0.0559334d0
      e(25) = 19.235851d0
      d(25) = 0.2738328d0
      e(26) = 6.2805882d0
      d(26) = 0.5328788d0
      e(27) = 1.9877968d0
      d(27) = 0.4009652d0
      return
      end
      subroutine ahdz1(e,s,p,d,n)
c
c  ----- alhrich's DZ  (8s4p) / [4s2p] {5111/31} ------
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li DZ -----
c
c li    (8s) / [4s]     {5111}
c scf energy is       -7.4309987271 a.u.
c
  100 continue
      e( 1) =    1484.2786000000d+00
      s( 1) =  0.7625174800d-03
      e( 2) =     222.6999200000d+00
      s( 2) =  0.5872458200d-02
      e( 3) =      50.6845460000d+00
      s( 3) =  0.2968575000d-01
      e( 4) =      14.3166660000d+00
      s( 4) =  0.1098372400d+00
      e( 5) =       4.5943260000d+00
      s( 5) =  0.2853470200d+00
      e( 6) =       1.5829911000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       0.5622588600d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.0478710540d+00
      s( 8) =  0.1000000000d+01
c
c additional p-GTO
c Ref.: see Ahlrichs and Taylor, JCP, 78 (1981), 315 
c geometric mean of DP-exponents.
c
      e( 9) =  0.17d+00
      p( 9) =  1.00d+00
c
      return
c
c     ----- be DZ -----
c
c be    (8s) / [4s]     {5111}
c scf energy is      -14.5703354622 a.u.
c
  120 continue
      e( 1) =    1429.9977041000d+00
      s( 1) =  0.1822530144d-02
      e( 2) =     214.7270646600d+00
      s( 2) =  0.1392938830d-01
      e( 3) =      48.8393147500d+00
      s( 3) =  0.6829046151d-01
      e( 4) =      13.6975120920d+00
      s( 4) =  0.2318872924d+00
      e( 5) =       4.2820865621d+00
      s( 5) =  0.4999316813d+00
      e( 6) =       1.3937553677d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       0.1904399024d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.0618134048d+00
      s( 8) =  0.1000000000d+01
c
c p set taken from DZ.3P basis
c
      e( 9) =   3.6316917145d+00
      p( 9) =    -.29033998305d-01
      e(10) =   .71695694366d+00
      p(10) =    -.16877854032d+00
      e(11) =   .19541932860d+00
      p(11) =    -.51403419628d+00
      e(12) =   .60515465890d-01  
      p(12) =  1.0000000000d+00
c
      return
c
c     ----- b DZ  -----
c
c B     (8s4p) / [4s2p]     {5111/31}
c scf energy is      -24.5228300102 a.u.
c
  140 continue
      e( 1) =    2410.2061000000d+00
      s( 1) =  0.1754991500d-02
      e( 2) =     361.8699200000d+00
      s( 2) =  0.1343640500d-01
      e( 3) =      82.3085930000d+00
      s( 3) =  0.6636902500d-01
      e( 4) =      23.1099420000d+00
      s( 4) =  0.2302529700d+00
      e( 5) =       7.2517259000d+00
      s( 5) =  0.5149821200d+00
      e( 6) =       2.3666469000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       0.3598707800d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.1116288900d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       6.0017073000d+00
      p( 9) = -0.3554532000d-01
      e(10) =       1.2401097000d+00
      p(10) = -0.1983450500d+00
      e(11) =       0.3366802400d+00
      p(11) = -0.5047872900d+00
      e(12) =       0.0955908620d+00
      p(12) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(13) =   .50d+00
      d(13) =    1.00d+00
c
      return
c
c     ----- c DZ  -----
c
c C     (8s4p) / [4s2p]     {5111/31}
c scf energy is      -37.6763463602 a.u.
c
  160 continue
      e( 1) =    3623.8613000000d+00
      s( 1) =  0.1633919100d-02
      e( 2) =     544.0462100000d+00
      s( 2) =  0.1252170100d-01
      e( 3) =     123.7433800000d+00
      s( 3) =  0.6211391400d-01
      e( 4) =      34.7632090000d+00
      s( 4) =  0.2181772900d+00
      e( 5) =      10.9333330000d+00
      s( 5) =  0.4980043100d+00
      e( 6) =       3.5744765000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       0.5748324500d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.1730364000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       9.4432819000d+00
      p( 9) =  0.3789545100d-01
      e(10) =       2.0017986000d+00
      p(10) =  0.2081817700d+00
      e(11) =       0.5462971800d+00
      p(11) =  0.5047416600d+00
      e(12) =       0.1520268400d+00
      p(12) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(13) =   .80d+00
      d(13) =   1.00d+00
c
      return
c
c     ----- n DZ -----
c
c N     (8s4p) / [4s2p]     {5111/31}
c scf energy is      -54.3794821043 a.u.
c
  180 continue
      e( 1) =    5071.9892000000d+00
      s( 1) =  0.1706723800d-02
      e( 2) =     761.4179100000d+00
      s( 2) =  0.1308769600d-01
      e( 3) =     173.1841800000d+00
      s( 3) =  0.6509825200d-01
      e( 4) =      48.6703900000d+00
      s( 4) =  0.2305160300d+00
      e( 5) =      15.3314480000d+00
      s( 5) =  0.5330047400d+00
      e( 6) =       5.0186710000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       0.8355477200d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.2462681400d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      13.5507220000d+00
      p( 9) =  0.4063084000d-01
      e(10) =       2.9178682000d+00
      p(10) =  0.2208579900d+00
      e(11) =       0.7983125200d+00
      p(11) =  0.5186014200d+00
      e(12) =       0.2190012500d+00
      p(12) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(13) =   1.00d+00
      d(13) =   1.00d+00
c
      return
c
c     ----- o  DZ ------
c
c O     (8s4p) / [4s2p]     {5111/31}
c scf energy is      -74.7737321990 a.u.
c
  200 continue
      e( 1) =    6773.3747000000d+00
      s( 1) =  0.1726613000d-02
      e( 2) =    1016.7970000000d+00
      s( 2) =  0.1324648400d-01
      e( 3) =     231.2673800000d+00
      s( 3) =  0.6602615700d-01
      e( 4) =      65.0084540000d+00
      s( 4) =  0.2352836100d+00
      e( 5) =      20.4998080000d+00
      s( 5) =  0.5497623600d+00
      e( 6) =       6.7160664000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       1.1471572000d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.3342325100d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      17.6942640000d+00
      p( 9) =  0.4359022400d-01
      e(10) =       3.8536027000d+00
      p(10) =  0.2317808700d+00
      e(11) =       1.0467465000d+00
      p(11) =  0.5145596900d+00
      e(12) =       0.2758604300d+00
      p(12) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(13) =   1.20d+00
      d(13) =   1.00d+00
c
      return
c
c     ----- f  DZ -----
c
c f     (8s4p) / [4s2p]     {5111/31}
c scf energy is      -99.3545187656 a.u.
c
  220 continue
      e( 1) =    8709.5462000000d+00
      s( 1) =  0.1393573900d-02
      e( 2) =    1307.4134000000d+00
      s( 2) =  0.1069507300d-01
      e( 3) =     297.3644100000d+00
      s( 3) =  0.5338930800d-01
      e( 4) =      83.6021540000d+00
      s( 4) =  0.1911247900d+00
      e( 5) =      26.3851210000d+00
      s( 5) =  0.4498753700d+00
      e( 6) =       8.6499174000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       1.5043291000d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.4340840500d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      22.6659940000d+00
      p( 9) =  0.3402887200d-01
      e(10) =       4.9760341000d+00
      p(10) =  0.1786531500d+00
      e(11) =       1.3476543000d+00
      p(11) =  0.3849732300d+00
      e(12) =       0.3477516800d+00
      p(12) =  0.1000000000d+01
c
c     ***** f DZP
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c depending on atomic charge, values down to 0.4 result in molecular
c optimization, e.g. for Li2F.
c
      d(13)=     1.40d0
      e(13)=     1.00d0
c
      return
c
c     ----- ne DZ -----
c
c ne    (8s4p) / [4s2p]     {5111/31}
c scf energy is     -128.4675279880 a.u.
c
  240 continue
      e( 1) =   10880.9393060000d+00
      s( 1) = -0.1701620401d-02
      e( 2) =    1633.3374648000d+00
      s( 2) = -0.1306253110d-01
      e( 3) =     371.4906217800d+00
      s( 3) = -0.6528233414d-01
      e( 4) =     104.4554596400d+00
      s( 4) = -0.2345130799d+00
      e( 5) =      32.9885758720d+00
      s( 5) = -0.5550008290d+00
      e( 6) =      10.8206201570d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       1.9072526872d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.5458668198d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      28.3915075350d+00
      p( 9) =  0.4644649139d-01
      e(10) =       6.2702151539d+00
      p(10) =  0.2419362441d+00
      e(11) =       1.6963316642d+00
      p(11) =  0.5120699660d+00
      e(12) =       0.4324828638d+00
      p(12) =  0.1000000000d+01
c 
c additional d-GTO 
c Ref.: PSD 16 (Huzinaga), p. 23
c
      e(13) =   1.888d+00
      d(13) = 1.000d+00
c
      return
      end
      subroutine ahdz2(e,s,p,d,n)
c
c  --- ahlrich's DZ (11s7p)/[6s4p] {521111/4111} contraction ---
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- na DZ -----
c
c na    (11s5p) / [6s2p]     {521111/41}
c scf energy is     -161.8381829223 a.u.
c
  100 continue
      e( 1) =   28048.8847970000d+00
      s( 1) =  0.4675237569d-03
      e( 2) =    4207.0284406000d+00
      s( 2) =  0.3613331347d-02
      e( 3) =     957.8847776600d+00
      s( 3) =  0.1852918001d-01
      e( 4) =     272.1529647500d+00
      s( 4) =  0.7129397506d-01
      e( 5) =      91.0913121190d+00
      s( 5) =  0.1904644741d+00
      e( 6) =      38.0786796090d+00
      s( 6) =  0.4182358213d+00
      e( 7) =      20.3069998720d+00
      s( 7) =  0.4544171396d+00
      e( 8) =       9.2665418351d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       2.2160241629d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.6920501915d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.0405244972d+00
      s(11) =  0.1000000000d+01
      e(12) =      75.4305991610d+00
      p(12) =  0.1542495086d-01
      e(13) =      17.2779118600d+00
      p(13) =  0.9975477048d-01
      e(14) =       5.1829377222d+00
      p(14) =  0.3121674127d+00
      e(15) =       1.6602128738d+00
      p(15) =  0.4926848330d+00
      e(16) =       0.5127119427d+00
      p(16) =  0.1000000000d+01
c
c additional p-function
c Ref.: R. Ahlrichs 12.08.92, as steepr s-GTO of DZ.
c which is optimal for NanClm
c
      e(17) = .052d+00
      p(17) =  1.00d+00
c
      return
c
c     ----- mg DZ -----
c
c mg    (11s5p) / [6s2p]     {521111/41}
c scf energy is      -199.5915911791 a.u.
c 
  120 continue
      e( 1) =   27472.3341500000d+00
      s( 1) =  0.6606035799d-03
      e( 2) =    4121.2426211000d+00
      s( 2) =  0.5099531712d-02
      e( 3) =     938.0872553900d+00
      s( 3) =  0.2603752216d-01
      e( 4) =     265.7469489400d+00
      s( 4) =  0.9932009521d-01
      e( 5) =      86.8908153910d+00
      s( 5) =  0.2705621898d+00
      e( 6) =      31.5180850080d+00
      s( 6) =  0.2007775183d+00
      e( 7) =      12.1050112980d+00
      s( 7) =  0.1283739710d+00
      e( 8) =       2.6398546036d+00
      s( 8) =  1.0000000000d+00
      e( 9) =       0.9014319486d+00
      s( 9) =  1.0000000000d+00
      e(10) =       0.1069914310d+00
      s(10) =  1.0000000000d+00
      e(11) =       0.0400575895d+00
      s(11) =  1.0000000000d+00
      e(12) =      98.0605066090d+00
      p(12) =  0.1448240121d-01
      e(13) =      22.5847139410d+00
      p(13) =  0.9554979322d-01
      e(14) =       6.8377400326d+00
      p(14) =  0.3079644881d+00
      e(15) =       2.2333606339d+00
      p(15) =  0.4994545351d+00
      e(16) =       0.7162575708d+00
      p(16) =  0.3147289558d+00
c
c two additional p-functions
c taken from DZ.3p basis set
c
      e(17) = .18914796195d+00
      p(17) =  1.00d+00
      e(18) = .53768755187d-01
      p(18) =  1.00d+00
      return
c
c     ----- al DZ -----
c
c al    (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -241.8551662242 a.u.
c 
  140 continue
      e( 1) =   32386.3810000000d+00
      s( 1) =  0.6266288800d-03
      e( 2) =    4858.4056000000d+00
      s( 2) =  0.4837597300d-02
      e( 3) =    1105.9177000000d+00
      s( 3) =  0.2470575000d-01
      e( 4) =     313.3923900000d+00
      s( 4) =  0.9426714800d-01
      e( 5) =     102.6549500000d+00
      s( 5) =  0.2565490100d+00
      e( 6) =      37.4098730000d+00
      s( 6) =  0.4546557900d+00
      e( 7) =      14.4578780000d+00
      s( 7) =  0.2958714400d+00
      e( 8) =       3.2405266000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       1.1620812000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.1769117000d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.0654139930d+00
      s(11) =  0.1000000000d+01
      e(12) =     145.6757800000d+00
      p(12) =  0.1057249000d-01
      e(13) =      33.8475080000d+00
      p(13) =  0.7309954900d-01
      e(14) =      10.4135470000d+00
      p(14) =  0.2579905100d+00
      e(15) =       3.5311030000d+00
      p(15) =  0.4753817700d+00
      e(16) =       1.2051750000d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.2718885000d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.0714646610d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(19) = .30d+00
      d(19) =  1.00d+00
c
      return
c
c     ----- si DZ -----
c
c si    (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -288.8288401073 a.u.
c
  160 continue
      e( 1) =   37778.3210000000d+00
      s( 1) =  0.6138218400d-03
      e( 2) =    5667.2439000000d+00
      s( 2) =  0.4739080100d-02
      e( 3) =    1290.0636000000d+00
      s( 3) =  0.2420910000d-01
      e( 4) =     365.6640500000d+00
      s( 4) =  0.9242286600d-01
      e( 5) =     119.9406500000d+00
      s( 5) =  0.2515488100d+00
      e( 6) =      43.8630830000d+00
      s( 6) =  0.4530613100d+00
      e( 7) =      17.0325370000d+00
      s( 7) =  0.2995574900d+00
      e( 8) =       3.8907716000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       1.4469558000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.2553205900d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.0933278120d+00
      s(11) =  0.1000000000d+01
      e(12) =     178.3497700000d+00
      p(12) =  0.1025016100d-01
      e(13) =      41.5494270000d+00
      p(13) =  0.7174730800d-01
      e(14) =      12.8429440000d+00
      p(14) =  0.2574644700d+00
      e(15) =       4.3972841000d+00
      p(15) =  0.4821427500d+00
      e(16) =       1.5311894000d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.3444351600d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.0978095200d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(19) = .35d+00
      d(19) =  1.00d+00
c
      return
c
c     ----- p DZ -----
c
c p     (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -340.6881043463 a.u.
c
  180 continue
      e( 1) =   43629.1463840000d+00
      s( 1) =  0.6029069432d-03
      e( 2) =    6544.9101761000d+00
      s( 2) =  0.4655146804d-02
      e( 3) =    1489.8784631000d+00
      s( 3) =  0.2378670338d-01
      e( 4) =     422.3830019300d+00
      s( 4) =  0.9086664055d-01
      e( 5) =     138.6974796600d+00
      s( 5) =  0.2474379471d+00
      e( 6) =      50.8653866990d+00
      s( 6) =  0.4585474744d+00
      e( 7) =      19.8269368190d+00
      s( 7) =  0.3075056301d+00
      e( 8) =       4.5933029207d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       1.7570187795d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.3435880478d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1245145268d+00
      s(11) =  0.1000000000d+01
      e(12) =     215.9706794900d+00
      p(12) =  0.9970358783d-02
      e(13) =      50.4289869420d+00
      p(13) =  0.7061426604d-01
      e(14) =      15.6536040090d+00
      p(14) =  0.2577344302d+00
      e(15) =       5.4085976775d+00
      p(15) =  0.4915306763d+00
      e(16) =       1.9189363936d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.4542082208d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.1319093708d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(19) = .45d+00
      d(19) = 1.00d+00
      return
c
c     ----- s DZ  -----
c
c s     (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -397.4662309332 a.u. 
c
  200 continue
      e( 1) =   49931.3060000000d+00
      s( 1) =  0.6327450100d-03
      e( 2) =    7490.2711000000d+00
      s( 2) =  0.4885864500d-02
      e( 3) =    1705.1055000000d+00
      s( 3) =  0.2497198200d-01
      e( 4) =     483.4788400000d+00
      s( 4) =  0.9545408500d-01
      e( 5) =     158.9058900000d+00
      s( 5) =  0.2601099800d+00
      e( 6) =      58.4131860000d+00
      s( 6) =  0.4609598200d+00
      e( 7) =      22.8417670000d+00
      s( 7) =  0.3130622900d+00
      e( 8) =       5.3498518000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       2.0927789000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.4422500000d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1587434200d+00
      s(11) =  0.1000000000d+01
      e(12) =     253.9644200000d+00
      p(12) =  0.9839016000d-02
      e(13) =      59.3945760000d+00
      p(13) =  0.7024181400d-01
      e(14) =      18.4938710000d+00
      p(14) =  0.2588755700d+00
      e(15) =       6.4332816000d+00
      p(15) =  0.4971731700d+00
      e(16) =       2.3192912000d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.5741082300d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.1615119300d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(19) = .55d+00
      d(19) = 1.00d+00
      return
c
c     ----- cl DZ -----
c
c cl    (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -459.4338276452 a.u. 
c
  220 continue
      e( 1) =   56684.8380000000d+00
      s( 1) =  0.6891293000d-03
      e( 2) =    8503.3349000000d+00
      s( 2) =  0.5321568400d-02
      e( 3) =    1935.7459000000d+00
      s( 3) =  0.2720543000d-01
      e( 4) =     548.9509300000d+00
      s( 4) =  0.1040553000d+00
      e( 5) =     180.5637200000d+00
      s( 5) =  0.2837805700d+00
      e( 6) =      66.5033430000d+00
      s( 6) =  0.4622066900d+00
      e( 7) =      26.0741590000d+00
      s( 7) =  0.3174546700d+00
      e( 8) =       6.1597909000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       2.4530894000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.5502628500d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1960931800d+00
      s(11) =  0.1000000000d+01
      e(12) =     293.0510100000d+00
      p(12) =  0.1064151800d-01
      e(13) =      68.6181320000d+00
      p(13) =  0.7642011100d-01
      e(14) =      21.4183990000d+00
      p(14) =  0.2833538200d+00
      e(15) =       7.4912517000d+00
      p(15) =  0.5450146900d+00
      e(16) =       2.7386125000d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.7084048600d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.1975506000d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, 
c depending on atomic charge, values down to 0.2 result in molecular
c optimization, e.g. for NanClm.
c
      e(19) =   .65d+00
      d(19) =  1.00d+00
c
      return
c
c     ----- ar DZ -----
c
c ar    (11s7p) / [6s4p]     {521111/4111}
c scf energy is     -526.7581276401 a.u.
c
  240 continue
      e( 1) =   63887.4728280000d+00
      s( 1) =  0.6067939377d-03
      e( 2) =    9583.7593589000d+00
      s( 2) =  0.4686048050d-02
      e( 3) =    2181.7219111000d+00
      s( 3) =  0.2396183814d-01
      e( 4) =     618.7769560400d+00
      s( 4) =  0.9170348506d-01
      e( 5) =     203.6640541300d+00
      s( 5) =  0.2503134202d+00
      e( 6) =      75.1338844930d+00
      s( 6) =  0.4626430101d+00
      e( 7) =      29.5236413720d+00
      s( 7) =  0.3209477320d+00
      e( 8) =       7.0236865832d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       2.8378923141d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.6675325064d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.2365339870d+00
      s(11) =  0.1000000000d+01
      e(12) =     332.0659199700d+00
      p(12) =  0.1010187759d-01
      e(13) =      77.8221111520d+00
      p(13) =  0.7282545709d-01
      e(14) =      24.3367346900d+00
      p(14) =  0.2706449036d+00
      e(15) =       8.5477344907d+00
      p(15) =  0.5182478798d+00
      e(16) =       3.1621383978d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.8540548107d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.2383338259d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c taken from tzp basis set in hondo basis set library
c
      e(19) =  .696d+00
      d(19) =  1.00d+00
c
      return
      end
      subroutine ahdz3(e,s,p,d,n)
c
c -- ahlrich's DZ (14s9p5d)/[8s5p3d] {62111111/33111/311} contr. --
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-18
      go to (100,120,140,160,180,200,220,240,260,280,300,320),nn
c
c k     (14s9p) / [8s5p]     {62111111/33111}
c optimized for atomic ground state k(2s)
c scf energy is     -599.1483848866 a.u.
c
  100 continue
      e( 1) =  182594.2732400000d+00
      s( 1) = -0.2277473593d-03
      e( 2) =   27369.0049860000d+00
      s( 2) = -0.1766398713d-02
      e( 3) =    6229.1737971000d+00
      s( 3) = -0.9194971209d-02
      e( 4) =    1764.5823486000d+00
      s( 4) = -0.3745508434d-01
      e( 5) =     577.0512058800d+00
      s( 5) = -0.1220445444d+00
      e( 6) =     210.2493336400d+00
      s( 6) = -0.2989897208d+00
      e( 7) =      82.6178306410d+00
      s( 7) = -0.4051469703d+00
      e( 8) =      33.2331845220d+00
      s( 8) = -0.2925315143d+00
      e( 9) =       8.1064902102d+00
      s( 9) =  0.1000000000d+01
      e(10) =       3.3340258818d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.8455442226d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.3282162653d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0364035285d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0176462865d+00
      s(14) =  0.1000000000d+01
      e(15) =     891.0544687200d+00
      p(15) =  0.2184291424d-02
      e(16) =     211.0158199600d+00
      p(16) =  0.1758911082d-01
      e(17) =      67.6714123780d+00
      p(17) =  0.8177747312d-01
      e(18) =      25.2714881000d+00
      p(18) =  0.2456563188d+00
      e(19) =      10.1389589690d+00
      p(19) =  0.4339837689d+00
      e(20) =       4.2018627167d+00
      p(20) =  0.3623766494d+00
      e(21) =       1.6029467214d+00
      p(21) =  0.1000000000d+01
      e(22) =       0.6146125305d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.2219930079d+00
      p(23) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Wachters, JCP, 52 (1970), 1033
c
      e(24) =   .041737d+00
      p(24) =  1.000000d+00
c
      return
c
c ca    (14s9p) / [8s5p]     {62111111/33111}
c optimized for atomic ground state ca(1s)
c scf energy is     -676.7402559879 a.u.
c 
  120 continue
      e( 1) =  202698.9551900000d+00
      s( 1) =  0.2229637833d-03
      e( 2) =   30382.4659150000d+00
      s( 2) =  0.1729318657d-02
      e( 3) =    6915.0833150000d+00
      s( 3) =  0.9002263178d-02
      e( 4) =    1959.0210984000d+00
      s( 4) =  0.3666988842d-01
      e( 5) =     640.9360207500d+00
      s( 5) =  0.1194098721d+00
      e( 6) =     233.9769683600d+00
      s( 6) =  0.2918252607d+00
      e( 7) =      92.2891983120d+00
      s( 7) =  0.4044152799d+00
      e( 8) =      37.2544799060d+00
      s( 8) =  0.2963128886d+00
      e( 9) =       9.1319843289d+00
      s( 9) =  0.1000000000d+01
      e(10) =       3.8177857927d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.0493475069d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.4286599152d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0628226390d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0260162379d+00
      s(14) =  0.1000000000d+01
      e(15) =    1019.7607638000d+00
      p(15) = -0.2059857528d-02
      e(16) =     241.5963075900d+00
      p(16) = -0.1665007982d-01
      e(17) =      77.6370120840d+00
      p(17) = -0.7776460231d-01
      e(18) =      29.1154194250d+00
      p(18) = -0.2418056124d+00
      e(19) =      11.7625935910d+00
      p(19) = -0.4325779868d+00
      e(20) =       4.9228924219d+00
      p(20) = -0.3673248154d+00
      e(21) =       1.9223257245d+00
      p(21) =  0.1000000000d+01
      e(22) =       0.7579619694d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.2843019910d+00
      p(23) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(24) =   .074979d+00
      p(24) =  1.000000d+00
c
      return
c
c sc    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c optimized for atomic ground state sc(2d)
c scf energy is     -759.7135729955 a.u. 
c 
  140 continue
      e( 1) =  224356.9429700000d+00
      s( 1) =  0.2199595627d-03
      e( 2) =   33628.7060600000d+00
      s( 2) =  0.1706041889d-02
      e( 3) =    7654.0171483000d+00
      s( 3) =  0.8881529467d-02
      e( 4) =    2168.4396566000d+00
      s( 4) =  0.3618370500d-01
      e( 5) =     709.6266385700d+00
      s( 5) =  0.1178392675d+00
      e( 6) =     259.3173132700d+00
      s( 6) =  0.2879164514d+00
      e( 7) =     102.4780375800d+00
      s( 7) =  0.3670513872d+00
      e( 8) =      41.4431396910d+00
      s( 8) =  0.2716199591d+00
      e( 9) =      10.2018556240d+00
      s( 9) =  0.1000000000d+01
      e(10) =       4.3044982698d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.2219457424d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.4972518144d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0727280775d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0293938293d+00
      s(14) =  0.1000000000d+01
      e(15) =    1154.9909592000d+00
      p(15) =  0.2010272272d-02
      e(16) =     273.7049661200d+00
      p(16) =  0.1630506985d-01
      e(17) =      88.0693273110d+00
      p(17) =  0.7652788902d-01
      e(18) =      33.1107600260d+00
      p(18) =  0.2395090404d+00
      e(19) =      13.4364407510d+00
      p(19) =  0.4323538023d+00
      e(20) =       5.6580984638d+00
      p(20) =  0.3696942939d+00
      e(21) =       2.2436704314d+00
      p(21) =  0.1000000000d+01
      e(22) =       0.8879800428d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.3328240973d+00
      p(23) =  0.1000000000d+01
      e(24) =      19.2385240380d+00
      d(24) =  0.2705466472d-01
      e(25) =       5.1173305876d+00
      d(25) =  0.1381324997d+00
      e(26) =       1.6549883906d+00
      d(26) =  0.3489530302d+00
      e(27) =       0.5401760459d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.1621062530d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .089748d+00
      p(29) = 1.000000d+00
c
      return
c
c ti    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is     -848.3790671511 a.u.
c optimized for atomic ground state ti(3f)
c calc. performed in oh-symmetry with 3d-occupation eg**2;
c 
  160 continue
      e( 1) =  247254.2047200000d+00
      s( 1) = -0.2239632908d-03
      e( 2) =   37060.7475770000d+00
      s( 2) = -0.1737119899d-02
      e( 3) =    8435.1700312000d+00
      s( 3) = -0.9043901944d-02
      e( 4) =    2389.8299020000d+00
      s( 4) = -0.3685081726d-01
      e( 5) =     782.2311987000d+00
      s( 5) = -0.1200420772d+00
      e( 6) =     286.0677640600d+00
      s( 6) = -0.2933689401d+00
      e( 7) =     113.2064433000d+00
      s( 7) = -0.3981781479d+00
      e( 8) =      45.8457203900d+00
      s( 8) = -0.2971303420d+00
      e( 9) =      11.3285248130d+00
      s( 9) =  0.1000000000d+01
      e(10) =       4.8125219381d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.3957658690d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.5632789519d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0805129593d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0319264578d+00
      s(14) =  0.1000000000d+01
      e(15) =    1298.4162107000d+00
      p(15) =  0.2024735219d-02
      e(16) =     307.7562578100d+00
      p(16) =  0.1647258601d-01
      e(17) =      99.1312339270d+00
      p(17) =  0.7766728462d-01
      e(18) =      37.3453055620d+00
      p(18) =  0.2432881999d+00
      e(19) =      15.2107126020d+00
      p(19) =  0.4425935919d+00
      e(20) =       6.4368981418d+00
      p(20) =  0.3804963554d+00
      e(21) =       2.5843928911d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.0229993308d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.3823395902d+00
      p(23) =  0.1000000000d+01
      e(24) =      23.4657510010d+00
      d(24) =  0.2654283272d-01
      e(25) =       6.3332300693d+00
      d(25) =  0.1380228459d+00
      e(26) =       2.0764079893d+00
      d(26) =  0.3532784293d+00
      e(27) =       0.6904116039d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.2109251522d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .101561d+00
      p(29) = 1.000000d+00
c
      return
c
c v    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is     -942.8517631586 a.u.
c optimized for atomic ground state v(4f)
c calculation performed in oh-symmetry with 3d-occupation t2g**3;
c 
  180 continue
      e( 1) =  271348.1594400000d+00
      s( 1) =  0.2196709488d-03
      e( 2) =   40672.0050780000d+00
      s( 2) =  0.1703860646d-02
      e( 3) =    9257.0761385000d+00
      s( 3) =  0.8871372088d-02
      e( 4) =    2622.7328773000d+00
      s( 4) =  0.3615432452d-01
      e( 5) =     858.5948018500d+00
      s( 5) =  0.1178089076d+00
      e( 6) =     314.1941865000d+00
      s( 6) =  0.2880253616d+00
      e( 7) =     124.4786212700d+00
      s( 7) =  0.4275791727d+00
      e( 8) =      50.4693993090d+00
      s( 8) =  0.3214423133d+00
      e( 9) =      12.5138469110d+00
      s( 9) =  0.1000000000d+01
      e(10) =       5.3447758082d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.5755631295d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.6306903307d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0879544416d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0342969197d+00
      s(14) =  0.1000000000d+01
      e(15) =    1450.0319327000d+00
      p(15) = -0.1951739710d-02
      e(16) =     343.7504988400d+00
      p(16) = -0.1592201424d-01
      e(17) =     110.8252170800d+00
      p(17) = -0.7538526025d-01
      e(18) =      41.8224645460d+00
      p(18) = -0.2366564132d+00
      e(19) =      17.0880737100d+00
      p(19) = -0.4334691987d+00
      e(20) =       7.2613218939d+00
      p(20) = -0.3743300135d+00
      e(21) =       2.9458145262d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.1650209680d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.4339937839d+00
      p(23) =  0.1000000000d+01
      e(24) =      27.3603358640d+00
      d(24) =  0.2664427003d-01
      e(25) =       7.4542765527d+00
      d(25) =  0.1399943647d+00
      e(26) =       2.4632270217d+00
      d(26) =  0.3576118644d+00
      e(27) =       0.8249954244d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.2526297057d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .111248d+00
      p(29) = 1.000000d+00
c
      return
c
c cr    (14s9p5d) / [8s5p3d]       {62111111/33111/311}
c scf energy is    -1043.3098073415 a.u. 
c optimized for atomic ground state cr(7s)
c
  200 continue
      e( 1) =  296910.3817500000d+00
      s( 1) =  0.2163182810d-03
      e( 2) =   44503.2979810000d+00
      s( 2) =  0.1677889232d-02
      e( 3) =   10129.0837310000d+00
      s( 3) =  0.8736768207d-02
      e( 4) =    2869.8125395000d+00
      s( 4) =  0.3561401119d-01
      e( 5) =     939.5285182800d+00
      s( 5) =  0.1161208156d+00
      e( 6) =     343.8766756200d+00
      s( 6) =  0.2842936594d+00
      e( 7) =     136.2724208100d+00
      s( 7) =  0.3961880476d+00
      e( 8) =      55.2749359100d+00
      s( 8) =  0.2993414345d+00
      e( 9) =      13.7537150650d+00
      s( 9) =  0.1000000000d+01
      e(10) =       5.8893064458d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.7366781917d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.6761901817d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.0921876160d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0361201071d+00
      s(14) =  0.1000000000d+01
      e(15) =    1608.8743670000d+00
      p(15) =  0.1936446367d-02
      e(16) =     381.4298700600d+00
      p(16) =  0.1583900730d-01
      e(17) =     123.0192684100d+00
      p(17) =  0.7535015415d-01
      e(18) =      46.4467821280d+00
      p(18) =  0.2394945140d+00
      e(19) =      19.0034947600d+00
      p(19) =  0.4409940312d+00
      e(20) =       8.0889965397d+00
      p(20) =  0.3801997556d+00
      e(21) =       3.3013703884d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.2877624138d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.4654684426d+00
      p(23) =  0.1000000000d+01
      e(24) =      27.5467337310d+00
      d(24) =  0.3063326043d-01
      e(25) =       7.4647940726d+00
      d(25) =  0.1560406567d+00
      e(26) =       2.4327524296d+00
      d(26) =  0.3700305778d+00
      e(27) =       0.7818781773d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.2197737763d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .120675d+00
      p(29) = 1.000000d+00
c
      return
c
c mn    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1149.8197211745 a.u. 
c 
  220 continue
      e( 1) =  323081.7447500000d+00
      s( 1) = -0.2185907002d-03
      e( 2) =   48425.9871470000d+00
      s( 2) = -0.1695536359d-02
      e( 3) =   11021.9284210000d+00
      s( 3) = -0.8829048037d-02
      e( 4) =    3122.8651683000d+00
      s( 4) = -0.3599386059d-01
      e( 5) =    1022.5702211000d+00
      s( 5) = -0.1173644143d+00
      e( 6) =     374.5666773200d+00
      s( 6) = -0.2872594276d+00
      e( 7) =     148.6501395300d+00
      s( 7) = -0.3802598779d+00
      e( 8) =      60.3790567850d+00
      s( 8) = -0.2894541087d+00
      e( 9) =      15.0612890370d+00
      s( 9) =  0.1000000000d+01
      e(10) =       6.4827687206d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.9535131494d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.7701785025d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1019202536d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0386197484d+00
      s(14) =  0.1000000000d+01
      e(15) =    1777.8786824000d+00
      p(15) = -0.1933161753d-02
      e(16) =     421.5775418000d+00
      p(16) = -0.1584400717d-01
      e(17) =     136.1133399900d+00
      p(17) = -0.7556898117d-01
      e(18) =      51.5067786340d+00
      p(18) = -0.2359743511d+00
      e(19) =      21.1532254410d+00
      p(19) = -0.4371879157d+00
      e(20) =       9.0477144185d+00
      p(20) = -0.3802419250d+00
      e(21) =       3.7308563614d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.4705331536d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.5439184098d+00
      p(23) =  0.1000000000d+01
      e(24) =      35.4248166460d+00
      d(24) =  0.2827075194d-01
      e(25) =       9.7813282131d+00
      d(25) =  0.1507166546d+00
      e(26) =       3.2670712909d+00
      d(26) =  0.3815855202d+00
      e(27) =       1.1028271910d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.3374808657d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .127650d+00
      p(29) = 1.000000d+00
c
      return
c
c fe    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1262.3875745795 a.u. 
c optimized for atomic ground state fe(5d)
c 
  240 continue
      e( 1) =  350644.1260700000d+00
      s( 1) =  0.2193042561d-03
      e( 2) =   52557.1044770000d+00
      s( 2) =  0.1701093534d-02
      e( 3) =   11962.1576030000d+00
      s( 3) =  0.8858517366d-02
      e( 4) =    3389.2987902000d+00
      s( 4) =  0.3611936931d-01
      e( 5) =    1109.9409814000d+00
      s( 5) =  0.1178052803d+00
      e( 6) =     406.7589212400d+00
      s( 6) =  0.2884662617d+00
      e( 7) =     161.5541919600d+00
      s( 7) =  0.3816201933d+00
      e( 8) =      65.6752111800d+00
      s( 8) =  0.2921257515d+00
      e( 9) =      16.4257132240d+00
      s( 9) =  0.1000000000d+01
      e(10) =       7.0930281780d+00
      s(10) =  0.1000000000d+01
      e(11) =       2.1594285091d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.8477134907d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1107949621d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0413168043d+00
      s(14) =  0.1000000000d+01
      e(15) =    1953.6045319000d+00
      p(15) =  0.1938161054d-02
      e(16) =     463.2931406600d+00
      p(16) =  0.1591619666d-01
      e(17) =     149.6749333100d+00
      p(17) =  0.7614915048d-01
      e(18) =      56.7071120310d+00
      p(18) =  0.2362443046d+00
      e(19) =      23.3410288660d+00
      p(19) =  0.4397888463d+00
      e(20) =      10.0113931640d+00
      p(20) =  0.3836258556d+00
      e(21) =       4.1564019580d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.6369582567d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.6037629747d+00
      p(23) =  0.1000000000d+01
      e(24) =      38.9664423920d+00
      d(24) =  0.2788406771d-01
      e(25) =      10.7989708770d+00
      d(25) =  0.1486313493d+00
      e(26) =       3.6129748398d+00
      d(26) =  0.3691295870d+00
      e(27) =       1.2130285825d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.3652397846d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .134915d+00
      p(29) = 1.000000d+00
c
      return
c
c co    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1381.3478575032 a.u. 
c optimized for atomic ground state co(4f)
c calculation performed in Oh-symm with 3d-occn eg**4 t2g**3;
c 
  260 continue
      e( 1) =  379411.7746100002d+00
      s( 1) = -0.2154791018d-03
      e( 2) =   56868.9343960000d+00
      s( 2) = -0.1671443777d-02
      e( 3) =   12943.6084620000d+00
      s( 3) = -0.8704447298d-02
      e( 4) =    3667.4532601000d+00
      s( 4) = -0.3549617452d-01
      e( 5) =    1201.1404591000d+00
      s( 5) = -0.1158102324d+00
      e( 6) =     440.3391947800d+00
      s( 6) = -0.2837460804d+00
      e( 7) =     175.0003529700d+00
      s( 7) = -0.3840417252d+00
      e( 8) =      71.1909815490d+00
      s( 8) = -0.2954541634d+00
      e( 9) =      17.8490079240d+00
      s( 9) =  0.1000000000d+01
      e(10) =       7.7276070252d+00
      s(10) =  0.1000000000d+01
      e(11) =       2.3708060829d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.9262827161d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1190862053d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0437823922d+00
      s(14) =  0.1000000000d+01
      e(15) =    2137.6007768000d+00
      p(15) = -0.1905793976d-02
      e(16) =     506.9729217400d+00
      p(16) = -0.1567838788d-01
      e(17) =     163.8752488100d+00
      p(17) = -0.7522962999d-01
      e(18) =      62.1518989910d+00
      p(18) = -0.2336923138d+00
      e(19) =      25.6323314350d+00
      p(19) = -0.4369310321d+00
      e(20) =      11.0207390270d+00
      p(20) = -0.3821086308d+00
      e(21) =       4.6024341479d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.8102832176d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.6656922371d+00
      p(23) =  0.1000000000d+01
      e(24) =      42.9243727770d+00
      d(24) =  0.2975487132d-01
      e(25) =      11.9408096010d+00
      d(25) =  0.1588571448d+00
      e(26) =       4.0037210467d+00
      d(26) =  0.3897061776d+00
      e(27) =       1.3412566852d+00
      d(27) =  0.4965261055d+00
      e(28) =       0.4001124188d+00
      d(28) =  0.3271993275d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =  .141308d+00
      p(29) = 1.000000d+00
c
      return
c
c ni    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1506.7922576508 a.u. 
c optimized for atomic ground state ni(3f)
c calculation performed in Oh-symm with 3d-occn t2g**6 eg**2;
c 
  280 continue
      e( 1) =  409347.8689900000d+00
      s( 1) =  0.2125177567d-03
      e( 2) =   61355.8392550000d+00
      s( 2) =  0.1648492989d-02
      e( 3) =   13964.8581280000d+00
      s( 3) =  0.8585328712d-02
      e( 4) =    3956.8598699000d+00
      s( 4) =  0.3501533398d-01
      e( 5) =    1296.0290474000d+00
      s( 5) =  0.1142767272d+00
      e( 6) =     475.2787651600d+00
      s( 6) =  0.2801524952d+00
      e( 7) =     188.9883877500d+00
      s( 7) =  0.3924628424d+00
      e( 8) =      76.9287369760d+00
      s( 8) =  0.3033130218d+00
      e( 9) =      19.3318287300d+00
      s( 9) =  0.1000000000d+01
      e(10) =       8.3876399732d+00
      s(10) =  0.1000000000d+01
      e(11) =       2.5894342897d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.0071374699d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1273257312d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0461898063d+00
      s(14) =  0.1000000000d+01
      e(15) =    2329.8384386000d+00
      p(15) = -0.1841418999d-02
      e(16) =     552.6069301500d+00
      p(16) = -0.1517368907d-01
      e(17) =     178.7119358400d+00
      p(17) = -0.7300466965d-01
      e(18) =      67.8417508130d+00
      p(18) = -0.2322484485d+00
      e(19) =      28.0280335330d+00
      p(19) = -0.4359604704d+00
      e(20) =      12.0764691990d+00
      p(20) = -0.3821378573d+00
      e(21) =       5.0694186684d+00
      p(21) =  0.1000000000d+01
      e(22) =       1.9911895236d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.7300688782d+00
      p(23) =  0.1000000000d+01
      e(24) =      47.0883036450d+00
      d(24) =  0.3152258030d-01
      e(25) =      13.1440672260d+00
      d(25) =  0.1685532569d+00
      e(26) =       4.4158589541d+00
      d(26) =  0.4093000465d+00
      e(27) =       1.4769996640d+00
      d(27) =  0.5150152778d+00
      e(28) =       0.4372896476d+00
      d(28) =  0.3396267123d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =   .146588d+00
      p(29) =  1.000000d+00
c
      return
c
c cu    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1638.8537855703 a.u. 
c optimized for atomic state cu(2s)
c
  300 continue
      e( 1) =  441087.2507000000d+00
      s( 1) = -0.2181417310d-03
      e( 2) =   66112.0211870000d+00
      s( 2) = -0.1692193588d-02
      e( 3) =   15047.0114250000d+00
      s( 3) = -0.8813713163d-02
      e( 4) =    4263.4273084000d+00
      s( 4) = -0.3595465952d-01
      e( 5) =    1396.3815797000d+00
      s( 5) = -0.1174297047d+00
      e( 6) =     511.9605578800d+00
      s( 6) = -0.2884426743d+00
      e( 7) =     203.4542694800d+00
      s( 7) = -0.4267889893d+00
      e( 8) =      82.7923370270d+00
      s( 8) = -0.3304412855d+00
      e( 9) =      20.8542856340d+00
      s( 9) =  0.1000000000d+01
      e(10) =       9.0410679584d+00
      s(10) =  0.1000000000d+01
      e(11) =       2.7518135173d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.0434856515d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1117229244d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0410410203d+00
      s(14) =  0.1000000000d+01
      e(15) =    2530.0965671000d+00
      p(15) =  0.1913794126d-02
      e(16) =     600.0979295400d+00
      p(16) =  0.1579767096d-01
      e(17) =     194.0820447900d+00
      p(17) =  0.7627125962d-01
      e(18) =      73.6718213770d+00
      p(18) =  0.2388145231d+00
      e(19) =      30.4473696900d+00
      p(19) =  0.4498001590d+00
      e(20) =      13.1227148750d+00
      p(20) =  0.3933768239d+00
      e(21) =       5.5214839972d+00
      p(21) =  0.1000000000d+01
      e(22) =       2.1457922130d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.7679748870d+00
      p(23) =  0.1000000000d+01
      e(24) =      47.3137437030d+00
      d(24) =  0.3239976043d-01
      e(25) =      13.1546884490d+00
      d(25) =  0.1682270653d+00
      e(26) =       4.3662885749d+00
      d(26) =  0.3849442960d+00
      e(27) =       1.4122065936d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.3884071300d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033  (2P)
c
      e(29) =  .155065d+00
      p(29) = 1.000000d+00
c
      return
c
c zn    (14s9p5d) / [8s5p3d]     {62111111/33111/311}
c scf energy is    -1777.7414723727 a.u. 
c 
  320 continue
      e( 1) =  472745.9416400000d+00
      s( 1) = -0.2170099403d-03
      e( 2) =   70858.0994840000d+00
      s( 2) = -0.1683375662d-02
      e( 3) =   16127.6397330000d+00
      s( 3) = -0.8767749579d-02
      e( 4) =    4569.7614320000d+00
      s( 4) = -0.3576876080d-01
      e( 5) =    1496.9765379000d+00
      s( 5) = -0.1168069854d+00
      e( 6) =     549.2564688800d+00
      s( 6) = -0.2867065307d+00
      e( 7) =     218.5901101200d+00
      s( 7) = -0.4268431703d+00
      e( 8) =      89.0688233280d+00
      s( 8) = -0.3325198992d+00
      e( 9) =      22.4757516930d+00
      s( 9) =  0.1000000000d+01
      e(10) =       9.7836990852d+00
      s(10) =  0.1000000000d+01
      e(11) =       3.0476097903d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.1751310173d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1433673810d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0507615043d+00
      s(14) =  0.1000000000d+01
      e(15) =    2738.9909849000d+00
      p(15) = -0.1870911026d-02
      e(16) =     649.7299519800d+00
      p(16) = -0.1546153148d-01
      e(17) =     210.2912907400d+00
      p(17) = -0.7475262879d-01
      e(18) =      79.9543603690d+00
      p(18) = -0.2328695546d+00
      e(19) =      33.1312904480d+00
      p(19) = -0.4401935958d+00
      e(20) =      14.3263958370d+00
      p(20) = -0.3873639145d+00
      e(21) =       6.0658646628d+00
      p(21) =  0.1000000000d+01
      e(22) =       2.3753719617d+00
      p(22) =  0.1000000000d+01
      e(23) =       0.8659631726d+00
      p(23) =  0.1000000000d+01
      e(24) =      56.0796648030d+00
      d(24) =  0.2959823379d-01
      e(25) =      15.7482590980d+00
      d(25) =  0.1587918488d+00
      e(26) =       5.3098927689d+00
      d(26) =  0.3798417808d+00
      e(27) =       1.7734652510d+00
      d(27) =  0.1000000000d+01
      e(28) =       0.5196323393d+00
      d(28) =  0.1000000000d+01
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(29) =   .162455d+00
      p(29) =  1.000000d+00
c
      return
      end
      subroutine ahdz4(e,s,p,d,n)
c
c -- ahlrich's DZ (14s10p5d)/[5s4p2d] {63311/5311/41} contractions --
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-30
      go to (100,120,140,160,180,200),nn
c
c ga    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state ga(2p)
c scf energy is    -1923.1570192810 a.u. 
c
  100 continue
      e( 1) =  505335.6329600000d+00
      s( 1) =  0.2097373573d-03
      e( 2) =   75742.7838440000d+00
      s( 2) =  0.1626971736d-02
      e( 3) =   17239.4663050000d+00
      s( 3) =  0.8474121088d-02
      e( 4) =    4884.9379362000d+00
      s( 4) =  0.3457114825d-01
      e( 5) =    1600.5426736000d+00
      s( 5) =  0.1128666588d+00
      e( 6) =     587.7569228800d+00
      s( 6) =  0.2767363154d+00
      e( 7) =     234.3010979700d+00
      s( 7) =  0.3851844447d+00
      e( 8) =      95.6116180320d+00
      s( 8) =  0.3019956501d+00
      e( 9) =      24.1697502180d+00
      s( 9) =  0.1000000000d+01
      e(10) =      10.5706034070d+00
      s(10) =  0.1000000000d+01
      e(11) =       3.3873246666d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.3329947879d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1924928200d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0705519234d+00
      s(14) =  0.1000000000d+01
      e(15) =    2971.6078970000d+00
      p(15) =  0.1599173568d-02
      e(16) =     705.1990378700d+00
      p(16) =  0.1322247523d-01
      e(17) =     228.7652045100d+00
      p(17) =  0.6379157288d-01
      e(18) =      87.4563910770d+00
      p(18) =  0.2177508273d+00
      e(19) =      36.5402319660d+00
      p(19) =  0.4139491890d+00
      e(20) =      15.9745797230d+00
      p(20) =  0.3753823569d+00
      e(21) =       6.8870611525d+00
      p(21) =  0.1000000000d+01
      e(22) =       2.8318589081d+00
      p(22) =  0.5855251845d+00
      e(23) =       1.1078009985d+00
      p(23) =  0.3367154269d+00
      e(24) =       0.2345085866d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.0636268812d+00
      p(25) =  0.1000000000d+01
      e(26) =      65.3446418070d+00
      d(26) =  0.2737834475d-01
      e(27) =      18.5005775170d+00
      d(27) =  0.1510566877d+00
      e(28) =       6.3162058406d+00
      d(28) =  0.3749184331d+00
      e(29) =       2.1638832514d+00
      d(29) =  0.4751026537d+00
      e(30) =       0.6668837843d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .207d+00
      d(31) =   1.000d+00
c
      return
c
c ge    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state ge(3p)
c scf energy is    -2075.2541874462 a.u.
c 
  120 continue
      e( 1) =  538981.2085900000d+00
      s( 1) = -0.2120086041d-03
      e( 2) =   80785.7689540000d+00
      s( 2) = -0.1644598990d-02
      e( 3) =   18387.3215390000d+00
      s( 3) = -0.8566070205d-02
      e( 4) =    5210.3369514000d+00
      s( 4) = -0.3494637061d-01
      e( 5) =    1707.4890772000d+00
      s( 5) = -0.1140615208d+00
      e( 6) =     627.5485407300d+00
      s( 6) = -0.2793632796d+00
      e( 7) =     250.5669613300d+00
      s( 7) = -0.3910094336d+00
      e( 8) =     102.3933310100d+00
      s( 8) = -0.3084509205d+00
      e( 9) =      25.9239233440d+00
      s( 9) =  0.1000000000d+01
      e(10) =      11.3890811910d+00
      s(10) =  0.1000000000d+01
      e(11) =       3.7472808901d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.5027733307d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.2444103389d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0912053935d+00
      s(14) =  0.1000000000d+01
      e(15) =    3196.8571212000d+00
      p(15) =  0.1509296599d-02
      e(16) =     758.7857470700d+00
      p(16) =  0.1248854362d-01
      e(17) =     246.4018616200d+00
      p(17) =  0.6025328674d-01
      e(18) =      94.4290894810d+00
      p(18) =  0.2263660966d+00
      e(19) =      39.6023332660d+00
      p(19) =  0.4315237187d+00
      e(20) =      17.3940232190d+00
      p(20) =  0.3949642301d+00
      e(21) =       7.5468872383d+00
      p(21) =  0.1000000000d+01
      e(22) =       3.1497947773d+00
      p(22) =  0.5977365286d+00
      e(23) =       1.2651485767d+00
      p(23) =  0.3494677277d+00
      e(24) =       0.3002403173d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.0873191872d+00
      p(25) =  0.1000000000d+01
      e(26) =      74.7764479030d+00
      d(26) =  0.2576022196d-01
      e(27) =      21.3077627270d+00
      d(27) =  0.1454099178d+00
      e(28) =       7.3452966023d+00
      d(28) =  0.3713636770d+00
      e(29) =       2.5656881541d+00
      d(29) =  0.4800315047d+00
      e(30) =       0.8198858734d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .246d+00
      d(31) =   1.000d+00
c
      return
c
c as    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state as(4s)
c scf energy is    -2234.1300405905 a.u.
c 
  140 continue
      e( 1) =  573690.5114000000d+00
      s( 1) = -0.2142347047d-03
      e( 2) =   85988.1834550000d+00
      s( 2) = -0.1661874968d-02
      e( 3) =   19571.4742700000d+00
      s( 3) = -0.8656182068d-02
      e( 4) =    5546.0309804000d+00
      s( 4) = -0.3531404437d-01
      e( 5) =    1817.8361947000d+00
      s( 5) = -0.1152313940d+00
      e( 6) =     668.6367549200d+00
      s( 6) = -0.2819264713d+00
      e( 7) =     267.3879594800d+00
      s( 7) = -0.3913353032d+00
      e( 8) =     109.4133525300d+00
      s( 8) = -0.3105265226d+00
      e( 9) =      27.7377668180d+00
      s( 9) =  0.1000000000d+01
      e(10) =      12.2385286750d+00
      s(10) =  0.1000000000d+01
      e(11) =       4.1265920360d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.6835020306d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.2995400986d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1131886664d+00
      s(14) =  0.1000000000d+01
      e(15) =    3425.0975358000d+00
      p(15) = -0.1473737753d-02
      e(16) =     813.0879709800d+00
      p(16) = -0.1220146593d-01
      e(17) =     264.2841476900d+00
      p(17) = -0.5885733672d-01
      e(18) =     101.5126370200d+00
      p(18) = -0.2057314885d+00
      e(19) =      42.7223506770d+00
      p(19) = -0.3928634381d+00
      e(20) =      18.8451139770d+00
      p(20) = -0.3619639761d+00
      e(21) =       8.2254455309d+00
      p(21) =  0.1000000000d+01
      e(22) =       3.4766363492d+00
      p(22) = -0.5964455342d+00
      e(23) =       1.4285918301d+00
      p(23) = -0.3502359236d+00
      e(24) =       0.3714979509d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.1124476071d+00
      p(25) =  0.1000000000d+01
      e(26) =      84.4426216440d+00
      d(26) =  0.2452046892d-01
      e(27) =      24.1880517260d+00
      d(27) =  0.1411022432d+00
      e(28) =       8.4037833603d+00
      d(28) =  0.3687438523d+00
      e(29) =       2.9811949434d+00
      d(29) =  0.4841029006d+00
      e(30) =       0.9792312379d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .293d+00
      d(31) =   1.000d+00
c
      return
c
c se    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state 3p
c scf energy is    -2399.7539361822 a.u.
c
  160 continue
      e( 1) =  609461.5616500000d+00
      s( 1) =  0.2066633429d-03
      e( 2) =   91349.7284690001d+00
      s( 2) =  0.1603149593d-02
      e( 3) =   20791.8449100000d+00
      s( 3) =  0.8350413760d-02
      e( 4) =    5892.0009402000d+00
      s( 4) =  0.3406658739d-01
      e( 5) =    1931.5795523000d+00
      s( 5) =  0.1111319992d+00
      e( 6) =     711.0217329900d+00
      s( 6) =  0.2716087745d+00
      e( 7) =     284.7656300000d+00
      s( 7) =  0.4035335848d+00
      e( 8) =     116.6728033600d+00
      s( 8) =  0.3220153758d+00
      e( 9) =      29.6113798470d+00
      s( 9) =  0.1000000000d+01
      e(10) =      13.1189215960d+00
      s(10) =  0.1000000000d+01
      e(11) =       4.5252274573d+00
      s(11) =  0.1000000000d+01
      e(12) =       1.8747678530d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.3585773767d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1361774385d+00
      s(14) =  0.1000000000d+01
      e(15) =    3656.2658173000d+00
      p(15) =  0.1471412900d-02
      e(16) =     868.0821167300d+00
      p(16) =  0.1218823195d-01
      e(17) =     282.3887909000d+00
      p(17) =  0.5878046890d-01
      e(18) =     108.6824951800d+00
      p(18) =  0.2238512571d+00
      e(19) =      45.8804283730d+00
      p(19) =  0.4278403799d+00
      e(20) =      20.3136332100d+00
      p(20) =  0.3956986798d+00
      e(21) =       8.9138677790d+00
      p(21) =  0.1000000000d+01
      e(22) =       3.8029430087d+00
      p(22) =  0.5889501703d+00
      e(23) =       1.5911663461d+00
      p(23) =  0.3414700424d+00
      e(24) =       0.4475557901d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.1323631907d+00
      p(25) =  0.1000000000d+01
      e(26) =      94.4937939930d+00
      d(26) =  0.2349045169d-01
      e(27) =      27.1865089330d+00
      d(27) =  0.1374936562d+00
      e(28) =       9.5088741989d+00
      d(28) =  0.3664676278d+00
      e(29) =       3.4175441801d+00
      d(29) =  0.4875229874d+00
      e(30) =       1.1481458665d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .338d+00
      d(31) =   1.000d+00
c
      return
c
c br    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state br(2p)
c scf energy is    -2572.3216875618 a.u. 
c 
  180 continue
      e( 1) =  646312.4103900000d+00
      s( 1) =  0.2084417341d-03
      e( 2) =   96873.0389000001d+00
      s( 2) =  0.1616955676d-02
      e( 3) =   22048.9675740000d+00
      s( 3) =  0.8422493476d-02
      e( 4) =    6248.3557482000d+00
      s( 4) =  0.3436087680d-01
      e( 5) =    2048.7467977000d+00
      s( 5) =  0.1120636106d+00
      e( 6) =     754.7082696300d+00
      s( 6) =  0.2736043108d+00
      e( 7) =     302.6968910300d+00
      s( 7) =  0.4056127980d+00
      e( 8) =     124.1686689800d+00
      s( 8) =  0.3254150105d+00
      e( 9) =      31.5442810700d+00
      s( 9) =  0.1000000000d+01
      e(10) =      14.0294359810d+00
      s(10) =  0.1000000000d+01
      e(11) =       4.9416679228d+00
      s(11) =  0.1000000000d+01
      e(12) =       2.0756618020d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.4210147444d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1605811437d+00
      s(14) =  0.1000000000d+01
      e(15) =    3890.7181798000d+00
      p(15) =  0.1511978712d-02
      e(16) =     923.8586949300d+00
      p(16) =  0.1252930300d-01
      e(17) =     300.7553708500d+00
      p(17) =  0.6040558694d-01
      e(18) =     115.9637691400d+00
      p(18) =  0.2110228896d+00
      e(19) =      49.0932710780d+00
      p(19) =  0.4034128832d+00
      e(20) =      21.8103920520d+00
      p(20) =  0.3738648104d+00
      e(21) =       9.6190135611d+00
      p(21) =  0.1000000000d+01
      e(22) =       4.1341549028d+00
      p(22) =  0.5993162501d+00
      e(23) =       1.7549598851d+00
      p(23) =  0.3397940021d+00
      e(24) =       0.5280153232d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.1560483796d+00
      p(25) =  0.1000000000d+01
      e(26) =     104.8573297000d+00
      d(26) =  0.2264931304d-01
      e(27) =      30.2800805400d+00
      d(27) =  0.1345624953d+00
      e(28) =      10.6514971270d+00
      d(28) =  0.3646967466d+00
      e(29) =       3.8705893430d+00
      d(29) =  0.4904647689d+00
      e(30) =       1.3242996938d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .389d+00
      d(31) =   1.000d+00
c
      return
c
c kr    (14s11p5d) / [8s6p2d]     {62111111/331211/41}
c optimized for atomic ground state kr(1s)
c scf energy is    -2751.9286711260
c
  200 continue
      e( 1) =  684221.9970000000d+00
      s( 1) =  0.2093888105d-03
      e( 2) =  102555.2004300001d+00
      s( 2) =  0.1624305928d-02
      e( 3) =   23342.3937750000d+00
      s( 3) =  0.8460820558d-02
      e( 4) =    6615.0835507000d+00
      s( 4) =  0.3451693173d-01
      e( 5) =    2169.3438677000d+00
      s( 5) =  0.1125458118d+00
      e( 6) =     799.6960795700d+00
      s( 6) =  0.2745085888d+00
      e( 7) =     321.1821210500d+00
      s( 7) =  0.3392875532d+00
      e( 8) =     131.9011167000d+00
      s( 8) =  0.2735993221d+00
      e( 9) =      33.5364730120d+00
      s( 9) =  0.1000000000d+01
      e(10) =      14.9698885520d+00
      s(10) =  0.1000000000d+01
      e(11) =       5.3755395083d+00
      s(11) =  0.1000000000d+01
      e(12) =       2.2858479626d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.4869023713d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1864101171d+00
      s(14) =  0.1000000000d+01
      e(15) =    4128.6993199000d+00
      p(15) = -0.1563604463d-02
      e(16) =     980.4753440500d+00
      p(16) = -0.1296144307d-01
      e(17) =     319.4013582000d+00
      p(17) = -0.6246445396d-01
      e(18) =     123.3618104000d+00
      p(18) = -0.2124567682d+00
      e(19) =      52.3623187440d+00
      p(19) = -0.4060431581d+00
      e(20) =      23.3354250930d+00
      p(20) = -0.3765220810d+00
      e(21) =      10.3407450510d+00
      p(21) =  0.1000000000d+01
      e(22) =       4.4694078019d+00
      p(22) =  0.5567193067d+00
      e(23) =       1.9177854904d+00
      p(23) =  0.3060445329d+00
      e(24) =       0.6125906362d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.1823640502d+00
      p(25) =  0.1000000000d+01
      e(26) =     115.5623666700d+00
      d(26) =  0.2194405819d-01
      e(27) =      33.4772471490d+00
      d(27) =  0.1321146024d+00
      e(28) =      11.8347997360d+00
      d(28) =  0.3632765030d+00
      e(29) =       4.3415492653d+00
      d(29) =  0.4930286009d+00
      e(30) =       1.5081356984d+00
      d(30) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(31) =   .443d+00
      d(31) =   1.000d+00
c
      return
      end
      subroutine ahdzsh(nucz,ipass,opol,itype,igauss,ng,
     *                  scale,scs,scp,scd,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(16,10), itypes(16,10), ioff(16,10)
      dimension igaup(2,10), itypep(2,10), ioffp(2,10)
      dimension kind(10)
      data kind/2,4,4,6,8,8,10,13,16,16/
      data igau /
     + 3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,3,1,0,0,0,0,0,0,0,0,0,0,
     + 5,2,1,1,1,1,4,1,0,0,0,0,0,0,0,0,
     + 5,2,1,1,1,1,4,1,0,0,0,0,0,0,0,0,
     + 5,2,1,1,1,1,4,1,1,1,0,0,0,0,0,0,
     + 6,2,1,1,1,1,1,1,3,3,1,1,1,0,0,0,
     + 6,2,1,1,1,1,1,1,3,3,1,1,1,3,1,1,
     + 6,2,1,1,1,1,1,1,3,3,1,2,1,1,4,1/
      data igaup /
     + 1,0,
     + 1,0,
     + 3,1,
     + 1,0,
     + 1,0,
     + 1,1,
     + 1,0,
     + 1,0,
     + 1,0,
     + 1,0/
      data itypes /
     + 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3 /
      data itypep/
     + 2,0,
     + 2,0,
     + 2,2,
     + 3,0,
     + 2,0,
     + 2,2,
     + 3,0,
     + 2,0,
     + 2,0,
     + 3,0 /
c
c    pointers to the dz basis
c
      data ioff/0, 3, 14*0,
     +          0, 5, 6, 7, 12*0,
     +          0, 5, 6, 7, 12*0,
     +          0, 5, 6, 7, 8,11,10*0,
     +          0, 5, 7, 8, 9,10,11,15,8*0,
     +          0, 5, 7, 8, 9,10,11,15,8*0,
     +          0, 5, 7, 8, 9,10,11,15,16,17,6*0,
     +          0, 6, 8, 9,10,11,12,13,14,17,20,21,22,3*0,
     +          0, 6, 8, 9,10,11,12,13,14,17,20,21,22,23,26,27,
     +          0, 6, 8, 9,10,11,12,13,14,17,20,21,23,24,25,29/
c
c    pointers to the dz+p basis
c
      data ioffp/4, 0,
     +           8, 0,
     +           8, 11,
     +          12, 0,
     +          16, 0,
     +          16, 17,
     +          18, 0,
     +          23, 0,
     +          28, 0,
     +          30, 0 /
c
      data done/1.0d+00/
c     data dzero/0.0d+00/
c
c     set values for the current ahlrichs dz  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.gt. 4) ind=4
      if(nucz.gt.10) ind=5
      if(nucz.gt.11) ind=6
      if(nucz.gt.12) ind=7
      if(nucz.gt.18) ind=8
      if(nucz.gt.20) ind=9
      if(nucz.gt.30) ind=10
      if(nucz.gt.36) then
         call caserr2(
     +          'ahlrichs dz basis sets only extend to krypton')
      end if
c
c
      mxpass=kind(ind)
c                 all atoms have one polarisation function excpt be
      npol = 1
      if(nucz.eq.4.or.nucz.eq.12 ) npol = 2
      if(opol) mxpass=mxpass+npol
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c              ******
c
      kkkk = ipass-kind(ind)
      if (kkkk.le.0) then
       igauss = igau(ipass,ind)
       ng = ioff(ipass,ind) 
       ityp = itypes(ipass,ind)
      else
       igauss = igaup(kkkk,ind)
       ng = ioffp(kkkk,ind) 
       ityp = itypep(kkkk,ind)
      endif
c
c          do not scale midi, or mini inner shells
c
      scale=done
      if(ityp.eq.1) scale=scs
      if(ityp.eq.2) scale=scp
      if(ityp.eq.3) scale=scd
100   itype = ityp + 15
      return
      end
      subroutine ahexsh(nucz,ipass,itype,igauss,ng,
     *                  scale,scs,scp,scd,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(18,8), itypes(18,8), ioff(18,8)
      dimension kind(8)
      data kind/4,9,9,11,14,14,11,18/
      data igau /
     + 4,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,1,1,4,1,1,0,0,0,0,0,0,0,0,0,
     + 6,1,1,1,1,1,4,1,1,1,1,0,0,0,0,0,0,0,
     + 6,1,1,1,1,1,1,1,4,1,1,1,1,1,0,0,0,0,
     + 6,1,1,1,1,1,1,1,5,1,1,1,1,1,0,0,0,0,
     + 5,1,2,1,1,1,1,5,1,1,1,0,0,0,0,0,0,0,
     + 6,1,1,1,1,1,1,1,1,1,6,1,1,1,1,1,1,1 /
      data itypes /
     + 1,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,3,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,3,0,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,2,0,0,0,0,
     + 1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3 /
c
c    pointers to the extended basis
c
      data ioff/0, 4, 5, 6, 14*0,
     +          0, 5, 6, 7, 8, 9,10,11,12,9*0,
     +          0, 5, 6, 7, 8, 9,10,14,15,9*0,
     +          0, 6, 7, 8, 9,10,11,15,16,17,18,7*0,
     +          0, 6, 7, 8, 9,10,11,12,13,17,18,19,20,21,4*0,
     +          0, 6, 7, 8, 9,10,11,12,13,18,19,20,21,22,4*0,
     +          0, 5, 6, 8, 9,10,11,12,17,18,19,7*0,
     +          0, 6, 7, 8, 9,10,11,12,13,14,15,21,22,23,24,25,26,27
     +            /
c
      data done/1.0d+00/
c     data dzero/0.0d+00/
c
c     set values for the current ahlrichs extended  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.eq.5.or.nucz.eq.6.or.nucz.eq.10) ind=4
      if(nucz.eq.7.or.nucz.eq.8.or.nucz.eq.9) ind=5
      if(nucz.eq.11) ind = 6
      if(nucz.eq.12) ind = 7
      if(nucz.ge.13) ind = 8
      if(nucz.gt.18) then
         call caserr2(
     +          'ahlrichs extended basis sets only extend to argon')
      end if
c
c
      mxpass=kind(ind)
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c
      igauss = igau(ipass,ind)
      ng = ioff(ipass,ind) 
      ityp = itypes(ipass,ind)
c
c          do not scale midi, or mini inner shells
c
      scale=done
      if(ityp.eq.1) scale=scs
      if(ityp.eq.2) scale=scp
      if(ityp.eq.3) scale=scd
100   itype = ityp + 15
      return
      end
      subroutine ahext(csinp,cpinp,cdinp,sc,scc,
     + nucz,intyp,nangm,nbfs,minf,maxf,loc,ngauss,ns,
     + ngsmax,ierr1,ierr2,nat)
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax),
     +  intyp(*),nangm(*),nbfs(*),minf(*),maxf(*),ns(*),scc(*)
c
      common/blkin/eex(50),ccs(50),ccp(50),ccd(50)
      common/junk/ptr(18192),iptr(4,maxat),iptrs(2,mxshel),
     + ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
     + cf(mxprim),cg(mxprim),
     + kstart(mxshel),katom(mxshel),ktype(mxshel),kng(mxshel),
     + kloc(mxshel),kmin(mxshel),kmax(mxshel),nshell
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
      data dzero,done,pt5,pt75/
     + 0.0d+00,1.0d+00,0.5d+00,0.75d+00/
      data pi32,tm10/5.56832799683170d+00,1.0d-10/
c     data tm6/1.0d-06/
      data scal1,scal2,scal3/1.00d0,1.00d0,1.00d0/
c
c     Ahlrichs, Scaefer and Horn, JCP 92 (1992) 2571
c     extended basis sets 
c
      scale = dzero
      ng = -2**20
      igauss = ng
      ityp = ng
      ierr1=0
      ierr2=0
      odone = .false.
      call vclr(eex,1,200)
c
c     ----- hydrogen to helium -----
c
      if (nucz .le. 2) then
        call ahext0(eex,ccs,ccp,nucz)
c
c     ----- lithium to neon -----
c
      else if (nucz .le. 10) then
        call ahext1(eex,ccs,ccp,ccd,nucz)
c
c     ----- sodium to argon -----
c
      else if (nucz .le. 18) then
        call ahext2(eex,ccs,ccp,ccd,nucz)
c
c     ----- past argon does not exist  -----
c
      else 
        call caserr2(
     +  'attempting to site ahlrichs ext. function on invalid centre')
      endif
c
      if (sc .le. 0.0d0) sc = scal1
      if (scc(1) .le. 0.0d0) scc(1) = scal2
      if (scc(2) .le. 0.0d0) scc(2) = scal3
c
c     ----- loop over each shell -----
c
      ipass = 0
  210 ipass = ipass+1
      scs = sc
      scp = scc(1)
      scd = scc(2)
      call ahexsh(nucz,ipass,ityp,igauss,ng,scale,
     *            scs,scp,scd,odone)
      if(odone) go to 220
c
c     ----- define the current shell -----
c
      nshell = nshell+1
      if(nshell.gt.mxshel) ierr1 = 1
      if(ierr1.ne.0) return
      ns(nat) = ns(nat)+1
      kmin(nshell) = minf(ityp)
      kmax(nshell) = maxf(ityp)
      kstart(nshell) = ngauss+1
      katom(nshell) = nat
      ktype(nshell) = nangm(ityp)
      intyp(nshell) = ityp
      kng(nshell) = igauss
      kloc(nshell) = loc+1
      ngauss = ngauss+igauss
      if(ngauss.gt.mxprim) ierr2 = 1
      if(ierr2.ne.0) return
      loc = loc+nbfs(ityp)
      k1 = kstart(nshell)
      k2 = k1+kng(nshell)-1
      scale = scale*scale
      do 440 i = 1,igauss
         k = k1+i-1
         ex(k) = eex(ng+i) * scale
         if(ityp.eq.16) csinp(k) = ccs(ng+i)
         if(ityp.eq.17) cpinp(k) = ccp(ng+i)
         if(ityp.eq.18) cdinp(k) = ccd(ng+i)
         cs(k) = dzero
         cp(k) = dzero
         cd(k) = dzero
         cf(k) = dzero
         if(ityp.eq.16) cs(k) = csinp(k)
         if(ityp.eq.17) cp(k) = cpinp(k)
         if(ityp.eq.18) cd(k) = cdinp(k)
  440 continue
c
c     ----- always unnormalize primitives -----
c
      do 460 k = k1,k2
         ee = ex(k)+ex(k)
         facs = pi32/(ee*sqrt(ee))
         facp = pt5*facs/ee
         facd = pt75*facs/(ee*ee)
         if(ityp.eq.16) cs(k) = cs(k)/sqrt(facs)
         if(ityp.eq.17) cp(k) = cp(k)/sqrt(facp)
         if(ityp.eq.18) cd(k) = cd(k)/sqrt(facd)
  460 continue
c
c     ----- if(normf.eq.0) normalize basis functions -----
c
      if (normf .eq. 1) go to 210
      facs = dzero
      facp = dzero
      facd = dzero
      do 510 ig = k1,k2
         do 500 jg = k1,ig
            ee = ex(ig)+ex(jg)
            fac = ee*sqrt(ee)
            dums = cs(ig)*cs(jg)/fac
            dump = pt5*cp(ig)*cp(jg)/(ee*fac)
            dumd = pt75*cd(ig)*cd(jg)/(ee*ee*fac)
            if (ig .eq. jg) go to 480
               dums = dums+dums
               dump = dump+dump
               dumd = dumd+dumd
  480       continue
            facs = facs+dums
            facp = facp+dump
            facd = facd+dumd
  500    continue
  510 continue
c
      fac=dzero
      if(ityp.eq.16.and. facs.gt.tm10) fac=done/sqrt(facs*pi32)
      if(ityp.eq.17.and. facp.gt.tm10) fac=done/sqrt(facp*pi32)
      if(ityp.eq.18.and. facd.gt.tm10) fac=done/sqrt(facd*pi32)
c
      do 550 ig = k1,k2
         if(ityp.eq.16) cs(ig) = fac*cs(ig)
         if(ityp.eq.17) cp(ig) = fac*cp(ig)
         if(ityp.eq.18) cd(ig) = fac*cd(ig)
  550 continue
      go to 210
c
  220 continue
      return
      end
      subroutine ahext0(e,s,p,n)
c
c     ----- ahlrich's contraction of (5s) -----
c     ----- hydrogen + helium (5s)/(3s) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      go to (100,120),n
c
c     ----- h extended -----
c
c h     (6s) / [3s]     {411}
c scf energy is        -.4999455704 a.u.
c
  100 continue
      e( 1) =      82.9215680000d+00
      s( 1) =  0.1997641000d-02
      e( 2) =      12.4523920000d+00
      s( 2) =  0.1528135600d-01
      e( 3) =       2.8330483000d+00
      s( 3) =  0.7530533200d-01
      e( 4) =       0.8000088400d+00
      s( 4) =  0.2563240600d+00
      e( 5) =       0.2585943900d+00
      s( 5) =  0.4976255500d+00
      e( 6) =       0.0899689190d+00
      s( 6) =  0.2967690300d+00
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e( 7) =  .80d+00
      p( 7) = 1.00d+00
c
      return
c
c     ----- he extended -----
c
c he    (6s) / [3s]     {411}
c scf energy is       -2.8611533453 a.u.
c
  120 continue
      e( 1) =     234.0637236100d+00
      s( 1) =  0.2586302956d-02
      e( 2) =      35.1740491020d+00
      s( 2) =  0.1952388427d-01
      e( 3) =       7.9911142747d+00
      s( 3) =  0.9098180782d-01
      e( 4) =       2.2124231186d+00
      s( 4) =  0.2720055756d+00
      e( 5) =       0.6670698619d+00
      s( 5) =  0.1000000000d+01
      e( 6) =       0.2089475538d+00
      s( 6) =  0.1000000000d+01
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e( 7) =   1.000d+00
      p( 7) =   1.000d+00
      return
c
      end
      subroutine ahext1(e,s,p,d,n)
c
c  ----- alhrich's extended basis sets
c
c li    (11s) / [7s]     {5111111}
c be TZ (10s) / [6s]     {511111}
c b     (11s7p) / [6s4p]     {611111/4111}
c c     (11s7p) / [6s4p]     {611111/4111}
c n     (13s8p) / [8s5p]     {61111111/41111}
c o     (13s8p) / [8s5p]     {61111111/41111}
c f     (13s8p) / [8s5p]     {61111111/41111}
c ne    (11s7p) / [6s4p]     {611111/4111}
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li extended -----
c
c atomic SCF calculation of Li(2S) 11s (uncontracted)
c SCF energy is       -7.4326805262 a.u. 
c li    (11s) / [7s]     {5111111}
c scf energy is    
c
  100 continue
      e( 1) =     7099.1720356d+00
      s( 1) =      .10728854441d-03
      e( 2) =     1063.8120419d+00
      s( 2) =      .83262807537d-03
      e( 3) =     242.11050680d+00
      s( 3) =      .43469069312d-02
      e( 4) =     68.556971437d+00
      s( 4) =      .17882632778d-01
      e( 5) =     22.340295421d+00
      s( 5) =      .60194574581d-01
      e( 6) =     8.0134749986d+00
      s( 6) =     1.00000000000d+00
      e( 7) =     3.0536958806d+00
      s( 7) =     1.00000000000d+00
      e( 8) =     1.2075774490d+00
      s( 8) =     1.00000000000d+00
      e( 9) =     .48202250929d+00
      s( 9) =     1.00000000000d+00
      e(10) =     .75096259666d-01
      s(10) =     1.00000000000d+00
      e(11) =     .28318829792d-01
      s(11) =     1.00000000000d+00
c
c two additional p-GTOs
c Ref.: C. Ochsenfeld - optimised for LinFm
c
      e(12) =  0.40d+00
      p(12) =  1.00d+00
      e(13) =  0.06d+00
      p(13) =  1.00d+00
c
      return
c
c     ----- be TZ -----
c
c Be    (10s) / [6s]     {511111}
c scf energy is      -14.5726557698 a.u.
c
  120 continue
      e( 1) =    6526.1477796000d+00
      s( 1) =  0.4373173971d-03
      e( 2) =     978.3549641700d+00
      s( 2) =  0.3388207195d-02
      e( 3) =     222.6586575600d+00
      s( 3) =  0.1756308508d-01
      e( 4) =      63.0130282810d+00
      s( 4) =  0.7058048133d-01
      e( 5) =      20.4687201820d+00
      s( 5) =  0.2226383623d+00
      e( 6) =       7.2692487387d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       2.7393701115d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       1.0615791118d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.1815844019d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.0588268200d+00
      s(10) =  0.1000000000d+01
c
c p set taken from TZ.3P basis
c be    (10s6p) / [6s3p]     {511111/411}
c optimized for 3P state, i.e. 1s(2)2s(1)2p(1)
c
      e(11) =  14.018558828d+00
      p(11) =      -.67974011741d-02
      e(12) =  3.1647356439d+00
      p(12) =      -.42491136104d-01
      e(13) =  .90157670087d+00
      p(13) =      -.17210557921d+00
      e(14) =  .30279153007d+00
      p(14) =      -.51791035465d+00
      e(15) =  .11279491554d+00
      p(15) =      1.0000000000d+00
      e(16) =  .42973443575d-01  
      p(16) =      1.0000000000d+00
c
      return
c
c     ----- b extended  -----
c
c B     (11s7p) / [6s4p]     {611111/4111}
c scf energy is      -24.5286838372 a.u. 
c
  140 continue
      e( 1) =   13036.4895210000d+00
      s( 1) =  0.2128998705d-03
      e( 2) =    1953.9818502000d+00
      s( 2) =  0.1650992395d-02
      e( 3) =     444.6779894300d+00
      s( 3) =  0.8590282445d-02
      e( 4) =     125.8369861000d+00
      s( 4) =  0.3498550959d-01
      e( 5) =      40.8628348630d+00
      s( 5) =  0.1147942475d+00
      e( 6) =      14.4630143870d+00
      s( 6) =  0.2912866912d+00
      e( 7) =       5.3553844892d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       2.0214609873d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.4695345401d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.2008719504d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.0776474044d+00
      s(11) =  0.1000000000d+01
      e(12) =      40.0062556110d+00
      p(12) =  0.1967896987d-02
      e(13) =       9.3434421395d+00
      p(13) =  0.1401390880d-01
      e(14) =       2.8365144955d+00
      p(14) =  0.5910436760d-01
      e(15) =       0.9913530774d+00
      p(15) =  0.1815970136d+00
      e(16) =       0.3759004240d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.1457794124d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.0561127384d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(19) =   .50d+00
      d(19) =    1.00d+00
c
      return
c
c     ----- c TZ  -----
c
c C     (11s7p) / [6s4p]     {611111/4111}
c scf energy is      -37.6880139943 a.u.
c
  160 continue
      e( 1) =   19841.3209220000d+00
      s( 1) =  0.2046894143d-03
      e( 2) =    2973.8160397000d+00
      s( 2) =  0.1587643858d-02
      e( 3) =     676.7627345800d+00
      s( 3) =  0.8267607387d-02
      e( 4) =     191.5065286000d+00
      s( 4) =  0.3377581584d-01
      e( 5) =      62.1831346210d+00
      s( 5) =  0.1118539369d+00
      e( 6) =      22.0160682120d+00
      s( 6) =  0.2892470480d+00
      e( 7) =       8.1594761071d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       3.0814906390d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.7549393891d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.3154547453d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1184049458d+00
      s(11) =  0.1000000000d+01
      e(12) =      61.6429463680d+00
      p(12) =  0.2079720580d-02
      e(13) =      14.4467635950d+00
      p(13) =  0.1520821294d-01
      e(14) =       4.4418116527d+00
      p(14) =  0.6521038115d-01
      e(15) =       1.5716359069d+00
      p(15) =  0.1944082183d+00
      e(16) =       0.5978320874d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.2301391977d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.0866295242d+00
      p(18) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(19) =   .80d+00
      d(19) =   1.00d+00
c
      return
c
c     ----- n extended -----
c
c N     (13s8p) / [8s5p]     {61111111/41111}
c scf energy is      -54.4008073560 a.u.
c
  180 continue
      e( 1) =   95630.8094580000d+00
      s( 1) =  0.3912904112d-04
      e( 2) =   14346.1749930000d+00
      s( 2) =  0.3034870465d-03
      e( 3) =    3266.3047300000d+00
      s( 3) =  0.1592520711d-02
      e( 4) =     924.8837440400d+00
      s( 4) =  0.6672144219d-02
      e( 5) =     301.5424456600d+00
      s( 5) =  0.2368543127d-01
      e( 6) =     108.6740554800d+00
      s( 6) =  0.7200608100d-01
      e( 7) =      42.1516475130d+00
      s( 7) =  0.1000000000d+01
      e( 8) =      17.2711312470d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       7.3970846393d+00
      s( 9) =  0.1000000000d+01
      e(10) =       3.2275502735d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.0161595792d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.3970171732d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.1514180421d+00
      s(13) =  0.1000000000d+01
      e(14) =     150.0821550900d+00
      p(14) = -0.8612151398d-03
      e(15) =      35.5068025700d+00
      p(15) = -0.6845460287d-02
      e(16) =      11.2482900530d+00
      p(16) = -0.3177164382d-01
      e(17) =       4.0898739977d+00
      p(17) = -0.1052638111d+00
      e(18) =       1.6217042186d+00
      p(18) =  0.1000000000d+01
      e(19) =       0.6639588641d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.2707180163d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1067810877d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(22) =   1.00d+00
      d(22) =   1.00d+00
c
      return
c
c     ----- o  extended ------
c
c O     (13s8p) / [8s5p]     {61111111/41111}
c scf energy is      -74.8091962776 a.u.
c
  200 continue
      e( 1) =  125915.6144200000d+00
      s( 1) =  0.3896449855d-04
      e( 2) =   18852.6248400000d+00
      s( 2) =  0.3029682311d-03
      e( 3) =    4289.0658708000d+00
      s( 3) =  0.1590406506d-02
      e( 4) =    1214.2628862000d+00
      s( 4) =  0.6665219319d-02
      e( 5) =     395.7828657500d+00
      s( 5) =  0.2368547120d-01
      e( 6) =     142.5727042600d+00
      s( 6) =  0.7216520458d-01
      e( 7) =      55.2758925880d+00
      s( 7) =  0.1000000000d+01
      e( 8) =      22.6545067020d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       9.7161217805d+00
      s( 9) =  0.1000000000d+01
      e(10) =       4.2513578184d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.3698362649d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.5331785722d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.2013026678d+00
      s(13) =  0.1000000000d+01
      e(14) =     191.0999422200d+00
      p(14) =  0.8567159544d-03
      e(15) =      45.2137232180d+00
      p(15) =  0.6835663454d-02
      e(16) =      14.3498879750d+00
      p(16) =  0.3191352301d-01
      e(17) =       5.2421108385d+00
      p(17) =  0.1043635562d+00
      e(18) =       2.0791956302d+00
      p(18) =  0.1000000000d+01
      e(19) =       0.8426971542d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.3361151273d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1286321834d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(22) =   1.20d+00
      d(22) =   1.00d+00
c
      return
c
c     ----- f  extended -----
c
c F     (13s8p) / [8s5p]     {61111111/41111}
c scf energy is      -99.4090463076 a.u.
c
  220 continue
      e( 1) =  159499.794970000d+00
      s( 1) =  0.3910257108d-04
      e( 2) =   23891.6072210000d+00
      s( 2) =  0.3038379742d-03
      e( 3) =    5440.0598370000d+00
      s( 3) =  0.1592711356d-02
      e( 4) =    1541.5172577000d+00
      s( 4) =  0.6668547704d-02
      e( 5) =     502.7160162400d+00
      s( 5) =  0.2370041137d-01
      e( 6) =     181.0928901700d+00
      s( 6) =  0.7232915975d-01
      e( 7) =      70.1904245180d+00
      s( 7) =  0.1000000000d+01
      e( 8) =      28.7684434500d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      12.3465130650d+00
      s( 9) =  0.1000000000d+01
      e(10) =       5.4103890967d+00
      s(10) =  0.1000000000d+01
      e(11) =       1.7707869891d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.6868535961d+00
      s(12) =  0.1000000000d+01
      e(13) =       0.2573325891d+00
      s(13) =  0.1000000000d+01
      e(14) =     241.1554214100d+00
      p(14) =  0.8499352567d-03
      e(15) =      57.0624883720d+00
      p(15) =  0.6813261121d-02
      e(16) =      18.1418257640d+00
      p(16) =  0.3200347052d-01
      e(17) =       6.6521340126d+00
      p(17) =  0.1036770503d+00
      e(18) =       2.6399926877d+00
      p(18) =  0.1000000000d+01
      e(19) =       1.0645959566d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.4195655711d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1575643121d+00
      p(21) =  0.1000000000d+01
c
c     ***** f TZP
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c depending on atomic charge, values down to 0.4 result in molecular
c optimization, e.g. for Li2F.
c
      d(22)=     1.40d0
      e(22)=     1.00d0
c
      return
c
c     ----- ne extended -----
c
c Ne    (11s7p) / [6s4p]     {611111/4111}
c scf energy is     -128.5445531091 a.u.
c
  240 continue
      e( 1) =   61048.0001200000d+00
      s( 1) = -0.1968202913d-03
      e( 2) =    9149.3065762000d+00
      s( 2) = -0.1527132572d-02
      e( 3) =    2082.0949756000d+00
      s( 3) = -0.7964229209d-02
      e( 4) =     589.1124426800d+00
      s( 4) = -0.3271437995d-01
      e( 5) =     191.2157297300d+00
      s( 5) = -0.1101660576d+00
      e( 6) =      67.6928695980d+00
      s( 6) = -0.2951088499d+00
      e( 7) =      25.1079239330d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       9.4812724488d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       2.5044729033d+00
      s( 9) =  0.1000000000d+01
      e(10) =       1.0074778194d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.3589118918d+00
      s(11) =  0.1000000000d+01
      e(12) =     174.9945784200d+00
      p(12) = -0.2595569832d-02
      e(13) =      41.1833058170d+00
      p(13) = -0.1966692327d-01
      e(14) =      12.9057153280d+00
      p(14) = -0.8522300900d-01
      e(15) =       4.6431715537d+00
      p(15) = -0.2309168175d+00
      e(16) =       1.7476880161d+00
      p(16) =  0.1000000000d+01
      e(17) =       0.6432317589d+00
      p(17) =  0.1000000000d+01
      e(18) =       0.2229674064d+00
      p(18) =  0.1000000000d+01
c 
c additional d-GTO 
c Ref.: PSD 16 (Huzinaga), p. 23
c
      e(19) =   1.888d+00
      d(19) = 1.000d+00
c
      return
      end
      subroutine ahext2(e,s,p,d,n)
c
c  ----- alhrich's extended basis sets
c
c Na    (14s9p) / [8s5p]     {61111111/51111}
c Mg TZ (12s6p) / [7s2p]     {5121111/51}
c Al    (15s12p) / [10s7p]     {6111111111/6111111}
c Si    (15s12p) / [10s7p]     {6111111111/6111111}
c P     (15s12p) / [10s7p]     {6111111111/6111111}
c S     (15s12p) / [10s7p]     {6111111111/6111111}
c Cl    (15s12p) / [10s7p]     {6111111111/6111111}
c Ar    (15s12p) / [10s7p]     {6111111111/6111111}
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- na extended -----
c
c Na    (14s9p) / [8s5p]     {61111111/51111}
c atomic SCF calculation of --->  Na  14s9p (uncontracted) <---
c SCF energy is     -161.8584370818 a.u. 
c
  100 continue
      e( 1) =    128983.96170d+00
      s( 1) =     0.79638839013d-04
      e( 2) =    19325.397745d+00
      s( 2) =     0.61835158793d-03
      e( 3) =    4398.1782277d+00
      s( 3) =     0.32377123408d-02
      e( 4) =    1245.5002338d+00
      s( 4) =     0.13441868300d-01
      e( 5) =    406.12829947d+00
      s( 5) =     0.46584429784d-01
      e( 6) =    146.38046993d+00
      s( 6) =     0.13280079867d+00
      e( 7) =    56.770624892d+00
      s( 7) =     0.10000000000d+01
      e( 8) =    23.238128393d+00
      s( 8) =     0.10000000000d+01
      e( 9) =    9.7383862833d+00
      s( 9) =     0.10000000000d+01
      e(10) =    3.1455878231d+00
      s(10) =     0.10000000000d+01
      e(11) =    1.2721344880d+00
      s(11) =     0.10000000000d+01
      e(12) =    .49609983518d+00
      s(12) =     0.10000000000d+01
      e(13) =    .58570908644d-01
      s(13) =     0.10000000000d+01
      e(14) =    .23344006904d-01
      s(14) =     0.10000000000d+01
      e(15) =    690.02808061d+00
      p(15) =    0.37562228078d-03
      e(16) =    163.56938429d+00
      p(16) =    0.31866976044d-02
      e(17) =    52.783476894d+00
      p(17) =    0.16379003297d-01
      e(18) =    19.777245383d+00
      p(18) =    0.59888023090d-01
      e(19) =    8.1199829651d+00
      p(19) =    0.15895776291d+00
      e(20) =    3.4924788410d+00
      p(20) =    0.10000000000d+01
      e(21) =    1.5092323157d+00
      p(21) =    0.10000000000d+01
      e(22) =    .64335696718d+00
      p(22) =    0.10000000000d+01
      e(23) =    .26090569846d+00
      p(23) =    0.10000000000d+01
c
c additional p-function
c Ref.: R. Ahlrichs 12.08.92, as steepr s-GTO of TZ.
c which is optimal for NanClm
c
      e(24) = .052d+00
      p(24) =  1.00d+00
c
      return
c
c     ----- mg TZ -----
c
c mg    (12s6p) / [7s2p]     {5121111/51}
c scf energy is      -199.6074427183 a.u.
c 
  120 continue
      e( 1) =   58088.3497780000d+00
      s( 1) =  0.2921796983d-03
      e( 2) =    8708.1180998000d+00
      s( 2) =  0.2264423666d-02
      e( 3) =    1981.6224182000d+00
      s( 3) =  0.1175603987d-01
      e( 4) =     560.4398599700d+00
      s( 4) =  0.4756994318d-01
      e( 5) =     181.7297610800d+00
      s( 5) =  0.1530804790d+00
      e( 6) =      64.6033357290d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      24.7421446700d+00
      s( 7) =  0.3988356060d+00
      e( 8) =      10.0246647160d+00
      s( 8) =  0.1818881521d+00
      e( 9) =       2.6947441296d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.9114525192d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1068967301d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.0400916171d+00
      s(12) =  0.1000000000d+01
      e(13) =     179.8337272400d+00
      p(13) =  0.5381976916d-02
      e(14) =      42.1130966640d+00
      p(14) =  0.3932575177d-01
      e(15) =      13.1173537750d+00
      p(15) =  0.1574904406d+00
      e(16) =       4.6229435410d+00
      p(16) =  0.3594056082d+00
      e(17) =       1.6685785801d+00
      p(17) =  0.4551622999d+00
      e(18) =       0.5854779723d+00
      p(18) =  0.2197049159d+00
c
c two additional p-functions
c taken from TZ.3p basis set
c
      e(19) = .18914796195d+00
      p(19) =  1.00d+00
      e(20) = .53768755187d-01
      p(20) =  1.00d+00
      return
c
c     ----- al extended -----
c
c Al    (15s12p) / [10s7p]     {6111111111/6111111}
c atomic SCF calculation of  Al(2P) 15s12p 
c SCF energy is     -241.8763584439 a.u. 
c 
  140 continue
c
      e( 1) =    390911.62280d+00
      s( 1) =    0.30390020231d-04
      e( 2) =    58552.131470d+00
      s( 2) =    0.23623499579d-03
      e( 3) =    13325.231083d+00
      s( 3) =    0.12398092467d-02
      e( 4) =    3773.9107696d+00
      s( 4) =    0.52022904942d-02
      e( 5) =    1230.9252794d+00
      s( 5) =    0.18542121223d-01
      e( 6) =    444.18207327d+00
      s( 6) =    0.57057843844d-01
      e( 7) =    172.98075234d+00
      s( 7) =    0.10000000000d+01
      e( 8) =    71.404651440d+00
      s( 8) =    0.10000000000d+01
      e( 9) =    30.815061920d+00
      s( 9) =    0.10000000000d+01
      e(10) =    13.549854113d+00
      s(10) =    0.10000000000d+01
      e(11) =    4.9931204147d+00
      s(11) =    0.10000000000d+01
      e(12) =    2.1702954031d+00
      s(12) =    0.10000000000d+01
      e(13) =    .92151393471d+00
      s(13) =    0.10000000000d+01
      e(14) =    .18192946185d+00
      s(14) =    0.10000000000d+01
      e(15) =    .65013418045d-01
      s(15) =    0.10000000000d+01
      e(16) =    1487.9288001d+00
      p(16) =    0.20212806816d-03
      e(17) =    352.61610178d+00
      p(17) =    0.17538339255d-02
      e(18) =    114.27847825d+00
      p(18) =    0.94578586690d-02
      e(19) =    43.263402008d+00
      p(19) =    0.36921543721d-01
      e(20) =    18.006948808d+00
      p(20) =    0.10906663749d+00
      e(21) =    7.9577668647d+00
      p(21) =    0.23289734653d+00
      e(22) =    3.6041351087d+00
      p(22) =    0.10000000000d+01
      e(23) =    1.6432322430d+00
      p(23) =    0.10000000000d+01
      e(24) =    .73789494707d+00
      p(24) =    0.10000000000d+01
      e(25) =    .25685663476d+00
      p(25) =    0.10000000000d+01
      e(26) =    .97317359939d-01
      p(26) =    0.10000000000d+01
      e(27) =    .36874004690d-01
      p(27) =    0.10000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(28) = .30d+00
      d(28) =  1.00d+00
c
      return
c
c     ----- si extended -----
c
c Si    (15s12p) / [10s7p]     {6111111111/6111111}
c atomic SCF calculation of --->  Si(3P) 15s12p  <---
c SCF energy is     -288.8539191577 a.u. (virial theorem =  2.000000000)
c
  160 continue
      e( 1) =    457646.72830d+00
      s( 1) =    0.30089847252d-04
      e( 2) =    68547.741767d+00
      s( 2) =    0.23390575839d-03
      e( 3) =    15600.029060d+00
      s( 3) =    0.12275993843d-02
      e( 4) =    4418.1985722d+00
      s( 4) =    0.51516011270d-02
      e( 5) =    1441.0974848d+00
      s( 5) =    0.18365987345d-01
      e( 6) =    520.05534417d+00
      s( 6) =    0.56559194116d-01
      e( 7) =    202.57427107d+00
      s( 7) =    0.10000000000d+01
      e( 8) =    83.675573696d+00
      s( 8) =    0.10000000000d+01
      e( 9) =    36.159765066d+00
      s( 9) =    0.10000000000d+01
      e(10) =    15.940674610d+00
      s(10) =    0.10000000000d+01
      e(11) =    5.9175528299d+00
      s(11) =    0.10000000000d+01
      e(12) =    2.6251161612d+00
      s(12) =    0.10000000000d+01
      e(13) =    1.1485428954d+00
      s(13) =    0.10000000000d+01
      e(14) =    .25688300480d+00
      s(14) =    0.10000000000d+01
      e(15) =    .91470687927d-01
      s(15) =    0.10000000000d+01
      e(16) =    1775.4637147d+00
      p(16) =    0.20211958447d-03
      e(17) =    420.73841187d+00
      p(17) =    0.17565536424d-02
      e(18) =    136.39078264d+00
      p(18) =    0.95147880868d-02
      e(19) =    51.688850321d+00
      p(19) =    0.37364206934d-01
      e(20) =    21.553982620d+00
      p(20) =    0.11095802480d+00
      e(21) =    9.5526930168d+00
      p(21) =    0.23774324891d+00
      e(22) =    4.3511289127d+00
      p(22) =    0.1000000000d+01
      e(23) =    2.0085914577d+00
      p(23) =    0.1000000000d+01
      e(24) =    .92200555058d+00
      p(24) =    0.1000000000d+01
      e(25) =    .34965464699d+00
      p(25) =    0.1000000000d+01
      e(26) =    .13782307603d+00
      p(26) =    0.1000000000d+01
      e(27) =    .53423784518d-01
      p(27) =    0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(28) = .35d+00
      d(28) =  1.00d+00
c
      return
c
c     ----- p extended -----
c
c P     (15s12p) / [10s7p]     {6111111111/6111111}
c scf energy is     -340.7181726894 a.u.
c
  180 continue
c
      e( 1) =  561557.2068700000d+00
      s( 1) =  0.3068244822d-04
      e( 2) =   84013.6540740000d+00
      s( 2) =  0.2390039674d-03
      e( 3) =   19095.5908820000d+00
      s( 3) =  0.1256574853d-02
      e( 4) =    5405.0469480000d+00
      s( 4) =  0.5277307456d-02
      e( 5) =    1762.2518102000d+00
      s( 5) =  0.1888350540d-01
      e( 6) =     634.5854892200d+00
      s( 6) =  0.5879938567d-01
      e( 7) =     245.7227608600d+00
      s( 7) =  0.1000000000d+01
      e( 8) =     100.6294056300d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      43.0942761610d+00
      s( 9) =  0.1000000000d+01
      e(10) =      18.8336982820d+00
      s(10) =  0.1000000000d+01
      e(11) =       6.8319094658d+00
      s(11) =  0.1000000000d+01
      e(12) =       3.0821107780d+00
      s(12) =  0.1000000000d+01
      e(13) =       1.3805386253d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.3383059852d+00
      s(14) =  0.1000000000d+01
      e(15) =       0.1202471108d+00
      s(15) =  0.1000000000d+01
      e(16) =    2227.1309170000d+00
      p(16) =  0.2055988048d-03
      e(17) =     526.4780937500d+00
      p(17) =  0.1800649076d-02
      e(18) =     170.6930632900d+00
      p(18) =  0.9826462411d-02
      e(19) =      64.8084343650d+00
      p(19) =  0.3910532967d-01
      e(20) =      27.0708330770d+00
      p(20) =  0.1186517742d+00
      e(21) =      12.0059023760d+00
      p(21) =  0.2635493413d+00
      e(22) =       5.4598894501d+00
      p(22) =  0.1000000000d+01
      e(23) =       2.5123390482d+00
      p(23) =  0.1000000000d+01
      e(24) =       1.1518272435d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.4700650854d+00
      p(25) =  0.1000000000d+01
      e(26) =       0.1884759898d+00
      p(26) =  0.1000000000d+01
      e(27) =       0.0732865965d+00
      p(27) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(28) = .45d+00
      d(28) = 1.00d+00
      return
c
c     ----- s extended  -----
c
c S    (15s12p) / [10s7p]     {6111111111/6111111}
c atomic SCF calculation S(3P) 15s12p  (uncontracted)
c SCF energy is     -397.5042520829
c
  200 continue
      e( 1) =    562582.94217d+00
      s( 1) =    0.32555313078d-04
      e( 2) =    84267.095293d+00
      s( 2) =    0.25307878126d-03
      e( 3) =    19177.573939d+00
      s( 3) =    0.13278889689d-02
      e( 4) =    5431.4887347d+00
      s( 4) =    0.55703647034d-02
      e( 5) =    1771.6941577d+00
      s( 5) =    0.19825647440d-01
      e( 6) =    639.46818436d+00
      s( 6) =    0.60824663880d-01
      e( 7) =    249.22540228d+00
      s( 7) =    1.00000000000d+00
      e( 8) =    103.12848700d+00
      s( 8) =    1.00000000000d+00
      e( 9) =    44.772288583d+00
      s( 9) =    1.00000000000d+00
      e(10) =    19.916090786d+00
      s(10) =    1.00000000000d+00
      e(11) =    7.2809413999d+00
      s(11) =    1.00000000000d+00
      e(12) =    3.2736230089d+00
      s(12) =    1.00000000000d+00
      e(13) =    1.5107815628d+00
      s(13) =    1.00000000000d+00
      e(14) =    .42276050040d+00
      s(14) =    1.00000000000d+00
      e(15) =    .14997087473d+00
      s(15) =    1.00000000000d+00
      e(16) =    2204.6101917d+00
      p(16) =    0.23823840683d-03
      e(17) =    522.43278113d+00
      p(17) =    0.20696556363d-02
      e(18) =    169.32765074d+00
      p(18) =    0.11203451408d-01
      e(19) =    64.176605416d+00
      p(19) =    0.43929126070d-01
      e(20) =    26.777252629d+00
      p(20) =    0.12885459179d+00
      e(21) =    11.854802616d+00
      p(21) =    0.26861316013d+00
      e(22) =     5.3915268953d+00
      p(22) =    0.10000000000d+01
      e(23) =     2.4887142591d+00
      p(23) =    0.10000000000d+01
      e(24) =     1.1191956887d+00
      p(24) =    0.10000000000d+01
      e(25) =     .48549015863d+00
      p(25) =    0.10000000000d+01
      e(26) =     .20042374307d+00
      p(26) =    0.10000000000d+01
      e(27) =     .79533765190d-01
      p(27) =    0.10000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(28) = .55d+00
      d(28) = 1.00d+00
      return
c
c     ----- cl extended -----
c
c Cl    (15s12p) / [10s7p]     {6111111111/6111111}
c scf energy is     -459.4811379727 a.u.
c
  220 continue
c 
      e( 1) =  676289.2509400004d+00
      s( 1) =  0.4733489843d-04
      e( 2) =  101324.8846700001d+00
      s( 2) =  0.3677275430d-03
      e( 3) =   23068.6874070000d+00
      s( 3) =  0.1929378774d-02
      e( 4) =    6536.5232367000d+00
      s( 4) =  0.8092053477d-02
      e( 5) =    2132.1154941000d+00
      s( 5) =  0.2891652107d-01
      e( 6) =     767.8623491000d+00
      s( 6) =  0.8977290647d-01
      e( 7) =     297.3065362500d+00
      s( 7) =  0.1000000000d+01
      e( 8) =     121.8337343900d+00
      s( 8) =  0.1000000000d+01
      e( 9) =      52.3488890430d+00
      s( 9) =  0.1000000000d+01
      e(10) =      23.0473669030d+00
      s(10) =  0.1000000000d+01
      e(11) =       8.1455363849d+00
      s(11) =  0.1000000000d+01
      e(12) =       3.6536757267d+00
      s(12) =  0.1000000000d+01
      e(13) =       1.6859777627d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.5124638621d+00
      s(14) =  0.1000000000d+01
      e(15) =       0.1819687006d+00
      s(15) =  0.1000000000d+01
      e(16) =    2951.9051373000d+00
      p(16) =  0.2290263555d-03
      e(17) =     699.1669237900d+00
      p(17) =  0.2001703065d-02
      e(18) =     226.9702838300d+00
      p(18) =  0.1097258076d-01
      e(19) =      86.2870692030d+00
      p(19) =  0.4406097980d-01
      e(20) =      36.0682529270d+00
      p(20) =  0.1354338698d+00
      e(21) =      15.9508255140d+00
      p(21) =  0.3064715799d+00
      e(22) =       3.2630665314d+00
      p(22) =  0.1000000000d+01
      e(23) =       7.1919747202d+00
      p(23) =  0.1000000000d+01
      e(24) =       1.4294058648d+00
      p(24) =  0.1000000000d+01
      e(25) =       0.6342979079d+00
      p(25) =  0.1000000000d+01
      e(26) =       0.2583807966d+00
      p(26) =  0.1000000000d+01
      e(27) =       0.0997627521d+00
      p(27) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, 
c depending on atomic charge, values down to 0.2 result in molecular
c optimization, e.g. for NanClm.
c
      e(28) =   .65d+00
      d(28) =  1.00d+00
c
      return
c
c     ----- ar extended -----
c
c Ar    (15s12p) / [10s7p]     {6111111111/6111111}
c atomic SCF calculation of  Ar(1S) 15s12p
c SCF energy is     -526.8167163760 a.u.
c
  240 continue
      e( 1) =    551271.36124d+00
      s( 1) =    0.44927502366d-04
      e( 2) =    82579.979599d+00
      s( 2) =    0.34918318169d-03
      e( 3) =    18793.603762d+00
      s( 3) =    0.18309112499d-02
      e( 4) =    5322.3829641d+00
      s( 4) =    0.76618405772d-02
      e( 5) =    1735.7943827d+00
      s( 5) =    0.27081134666d-01
      e( 6) =    626.20006837d+00
      s( 6) =    0.81565846399d-01
      e( 7) =    243.72642332d+00
      s( 7) =    0.10000000000d+01
      e( 8) =    100.60836087d+00
      s( 8) =    0.10000000000d+01
      e( 9) =    43.512854225d+00
      s( 9) =    0.10000000000d+01
      e(10) =    19.336809634d+00
      s(10) =    0.10000000000d+01
      e(11) =    7.5452518952d+00
      s(11) =    0.10000000000d+01
      e(12) =    3.2346169053d+00
      s(12) =    0.10000000000d+01
      e(13) =    1.2267045870d+00
      s(13) =    0.10000000000d+01
      e(14) =    .52677375475d+00
      s(14) =    0.10000000000d+01
      e(15) =    .19678029796d+00
      s(15) =    0.10000000000d+01
      e(16) =    2872.8216904d+00
      p(16) =    0.23972548748d-03
      e(17) =    680.74215068d+00
      p(17) =    0.20873049741d-02
      e(18) =    220.69703811d+00
      p(18) =    0.11353165704d-01
      e(19) =    83.744147790d+00
      p(19) =    0.44846404581d-01
      e(20) =    35.015348182d+00
      p(20) =    0.13223285188d+00
      e(21) =    15.550076569d+00
      p(21) =    0.27680246911d+00
      e(22) =    7.1123829575d+00
      p(22) =    0.10000000000d+01
      e(23) =    3.3084364096d+00
      p(23) =    0.10000000000d+01
      e(24) =    1.4741831979d+00
      p(24) =    0.10000000000d+01
      e(25) =    .65340322850d+00
      p(25) =    0.10000000000d+01
      e(26) =    .27445125440d+00
      p(26) =    0.10000000000d+01
      e(27) =    .10971473038d+00
      p(27) =    0.10000000000d+01
c
c additional d-GTO
c taken from tzp basis set in hondo basis set library
c
      e(28) =  .696d+00
      d(28) =  1.00d+00
c
      return
      end
      subroutine ahsv0(e,s,p,n)
c
c     ----- ahlrich's contraction of huzinaga's (4s) -----
c     ----- hydrogen (4s)/(2s) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      go to (100,120),n
c
c     ----- h SV,DZ -----
c
c h     (4s) / [2s]     {31}
c SCF energy is        -.4992784057 a.u. 
c
  100 continue
      e( 1) =     13.0107010000d+00
      s( 1) = 0.1968215800d-01
      e( 2) =      1.9622572000d+00
      s( 2) = 0.1379652400d+00
      e( 3) =      0.4445379600d+00
      s( 3) = 0.4783193500d+00
      e( 4) =      0.1219496200d+00
      s( 4) = 0.1000000000d+01
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e(5) =  .80d+00
      p(5) = 1.00d+00
c
      return
c
c     ----- he SV,DZ -----
c
c he    (4s) / [2s]     {31}
c SCF energy is       -2.8551604793 a.u. 
c
  120 continue
      e( 1) =     38.3549367370d+00
      s( 1) = 0.2381428891d-01
      e( 2) =      5.7689081479d+00
      s( 2) = 0.1549090678d+00
      e( 3) =      1.2399407035d+00
      s( 3) = 0.4699809663d+00
      e( 4) =      0.2975781595d+00
      s( 4) = 0.1000000000d+01
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e(5) =   1.000d+00
      p(5) =   1.000d+00
      return
c
      end
      subroutine ahsv1(e,s,p,d,n)
c
c  ----- alhrich's SV (511/31) contraction of primitive (7s,4p) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li SV -----
c
c li    (7s) / [3s]     {511}
c SCF energy is       -7.4250640446 a.u. 
c
  100 continue
      e( 1) =    266.2778551600d+00
      s( 1) = 0.6492015033d-02
      e( 2) =     40.0697834470d+00
      s( 2) = 0.4774786322d-01
      e( 3) =      9.0559944389d+00
      s( 3) = 0.2026879611d+00
      e( 4) =      2.4503009051d+00
      s( 4) = 0.4860657482d+00
      e( 5) =      0.7220957186d+00
      s( 5) = 0.4362697796d+00
      e( 6) =      0.0528108847d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.0209609488d+00
      s( 7) = 0.1000000000d+01
c
c additional p-GTO
c Ref.: see Ahlrichs and Taylor, JCP, 78 (1981), 315 
c geometric mean of DP-exponents.
c
      e( 8) =  0.17d+00
      p( 8) =  1.00d+00
c
      return
c
c     ----- be SV -----
c
c be    (7s) / [3s]     {511}
c SCF energy is      -14.5537268707 a.u. 
c
  120 continue
      e( 1) =    515.1861612500d+00
      s( 1) = 0.5561530798d-02
      e( 2) =     77.5110375950d+00
      s( 2) = 0.4119006806d-01
      e( 3) =     17.5524816930d+00
      s( 3) = 0.1791337811d+00
      e( 4) =      4.8028940596d+00
      s( 4) = 0.4473671646d+00
      e( 5) =      1.4516214316d+00
      s( 5) = 0.4200958192d+00
      e( 6) =      0.1328163362d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.0458373722d+00
      s( 7) = 0.1000000000d+01
c
c p set taken from SV.3P basis
c
      e( 8) =   3.6316917145d+00
      p( 8) =    -.29033998305d-01
      e( 9) =   .71695694366d+00
      p( 9) =    -.16877854032d+00
      e(10) =   .19541932860d+00
      p(10) =    -.51403419628d+00
      e(11) =   .60515465890d-01  
      p(11) =  1.0000000000d+00
c
      return
c
c     ----- b SV  -----
c
c B     (7s4p) / [3s2p]     {511/31}
c SCF energy is     -24.49755868320 a.u.  
c
  140 continue
      e( 1) =    839.3183008600d+00
      s( 1) = -.5592920107d-02
      e( 2) =    126.2646484300d+00
      s( 2) = -.4156552077d-01
      e( 3) =     28.6206007630d+00
      s( 3) = -.1829981698d+00
      e( 4) =      7.8793722710d+00
      s( 4) = -.4654039187d+00
      e( 5) =      2.4088857172d+00
      s( 5) = -.4417388479d+00
      e( 6) =      0.2510510904d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.0836488661d+00
      s( 7) = 0.1000000000d+01
      e( 8) =      6.0332223619d+00
      p( 8) = -.3560367246d-01
      e( 9) =      1.2499157866d+00
      p( 9) = -.1989577577d+00
      e(10) =      0.3387167635d+00
      p(10) = -.5085020262d+00
      e(11) =      0.0964156324d+00
      p(11) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(12) =   .50d+00
      d(12) =    1.00d+00
c
      return
c
c     ----- c SV  -----
c
c C     (7s4p) / [3s2p]     {511/31}
c SCF energy is      -37.6411433792 a.u.  
c
  160 continue
      e( 1) =   1238.4016938000d+00
      s( 1) = 0.5456883208d-02
      e( 2) =    186.2900499200d+00
      s( 2) = 0.4063840921d-01
      e( 3) =     42.2511763460d+00
      s( 3) = 0.1802559389d+00
      e( 4) =     11.6765579320d+00
      s( 4) = 0.4631512176d+00
      e( 5) =      3.5930506482d+00
      s( 5) = 0.4408717331d+00
      e( 6) =      0.4024514736d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.1309018267d+00
      s( 7) = 0.1000000000d+01
      e( 8) =      9.4680970621d+00
      p( 8) = 0.3838787173d-01
      e( 9) =      2.0103545142d+00
      p( 9) = 0.2111702511d+00
      e(10) =      0.5477100471d+00
      p(10) = 0.5132817211d+00
      e(11) =      0.1526861380d+00
      p(11) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(12) =   .80d+00
      d(12) =   1.00d+00
c
      return
c
c     ----- n SV -----
c
c N     (7s4p) / [3s2p]     {511/31}
c SCF energy is      -54.3329250250 a.u. 
c
  180 continue
      e( 1) =   1712.8415853000d+00
      s( 1) = -.5393412531d-02
      e( 2) =    257.6481267700d+00
      s( 2) = -.4022158112d-01
      e( 3) =     58.4582458530d+00
      s( 3) = -.1793114499d+00
      e( 4) =     16.1983679050d+00
      s( 4) = -.4637631782d+00
      e( 5) =      5.0052600809d+00
      s( 5) = -.4417142266d+00
      e( 6) =      0.5873185657d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.1876459225d+00
      s( 7) = 0.1000000000d+01
      e( 8) =     13.5714702330d+00
      p( 8) = -.4007239885d-01
      e( 9) =      2.9257372874d+00
      p( 9) = -.2180704503d+00
      e(10) =      0.7992775075d+00
      p(10) = -.5129446605d+00
      e(11) =      0.2195434803d+00
      p(11) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(12) =   1.00d+00
      d(12) =   1.00d+00
c
      return
c
c     ----- o  SV ------
c
c O     (7s4p) / [3s2p]     {511/31}
c SCF energy is     -74.71374601625 a.u.  
c
  200 continue
      e( 1) =   2266.1767785000d+00
      s( 1) = -.5343180993d-02
      e( 2) =    340.8701019100d+00
      s( 2) = -.3989003923d-01
      e( 3) =     77.3631351670d+00
      s( 3) = -.1785391199d+00
      e( 4) =     21.4796449400d+00
      s( 4) = -.4642768496d+00
      e( 5) =      6.6589433124d+00
      s( 5) = -.4430974517d+00
      e( 6) =      0.8097597567d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.2553077223d+00
      s( 7) = 0.1000000000d+01
      e( 8) =     17.7215043170d+00
      p( 8) = 0.4339457319d-01
      e( 9) =      3.8635505440d+00
      p( 9) = 0.2309412077d+00
      e(10) =      1.0480920883d+00
      p(10) = 0.5137531106d+00
      e(11) =      0.2764154441d+00
      p(11) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(12) =   1.20d+00
      d(12) =   1.00d+00
c
      return
c
c     ----- f  SV -----
c
c f     (7s4p) / [3s2p]     {511/31}
c  SCF energy is      -99.2796895705 a.u.  
c
  220 continue
      e( 1) =   2894.8325990000d+00
      s( 1) = -.5340825552d-02
      e( 2) =    435.4193912000d+00
      s( 2) = -.3990425887d-01
      e( 3) =     98.8433288660d+00
      s( 3) = -.1791276804d+00
      e( 4) =     27.4851980010d+00
      s( 4) = -.4675809083d+00
      e( 5) =      8.5405498171d+00
      s( 5) = -.4465313102d+00
      e( 6) =      1.0654578038d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.3324734675d+00
      s( 7) = 0.1000000000d+01
      e( 8) =     22.6966339240d+00
      p( 8) = -.4521287444d-01
      e( 9) =      4.9872339257d+00
      p( 9) = -.2375431707d+00
      e(10) =      1.3491613954d+00
      p(10) = -.5128735359d+00
      e(11) =      0.3482988198d+00
      p(11) = 0.1000000000d+01
c
c     ***** f SVP
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c depending on atomic charge, values down to 0.4 result in molecular
c optimization, e.g. for Li2F.
c
      d(12)=     1.40d0
      e(12)=     1.00d0
c
      return
c
c     ----- ne SV -----
c
c ne    (7s4p) / [3s2p]     {511/31}
c SCF energy is     -128.37640681003 a.u. 
c
  240 continue
      e( 1) =   3598.9736625000d+00
      s( 1) = -.5325929700d-02
      e( 2) =    541.3207311200d+00
      s( 2) = -.3981741797d-01
      e( 3) =    122.9045006200d+00
      s( 3) = -.1791435819d+00
      e( 4) =     34.2166170220d+00
      s( 4) = -.4689358298d+00
      e( 5) =     10.6505841240d+00
      s( 5) = -.4478253758d+00
      e( 6) =      1.3545953960d+00
      s( 6) = 0.1000000000d+01
      e( 7) =      0.4191936264d+00
      s( 7) = 0.1000000000d+01
      e( 8) =     28.4240537850d+00
      p( 8) = -.4603194480d-01
      e( 9) =      6.2822510953d+00
      p( 9) = -.2399318304d+00
      e(10) =      1.6978715079d+00
      p(10) = -.5087172496d+00
      e(11) =      0.4330070017d+00
      p(11) = 0.1000000000d+01
c 
c additional d-GTO 
c Ref.: PSD 16 (Huzinaga), p. 23
c
      e(12) =   1.888d+00
      d(12) = 1.000d+00
c
      return
      end
      subroutine ahsv2(e,s,p,d,n)
c
c     ----- ahlrich's SV (10s7p/4s3p) contractions  -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- na SV -----
c
c na    (10s5p) / [4s1p]     {5311/5}
c only a factor for splitting the exponent of the 2s-shell in the sz9s5p
c was optimized. therefor this is no minimum on the energy surface!
c scf energy is     -161.7917031467 a.u.
c
  100 continue
      e( 1) =   4098.2003908000d+00
      s( 1) = -.5853591188d-02
      e( 2) =    616.4937403100d+00
      s( 2) = -.4364716187d-01
      e( 3) =    139.9664400100d+00
      s( 3) = -.1943146588d+00
      e( 4) =     39.0734410510d+00
      s( 4) = -.4868506573d+00
      e( 5) =     11.9298472050d+00
      s( 5) = -.4188170514d+00
      e( 6) =     20.6599660300d+00
      s( 6) = 0.8594968985d-01
      e( 7) =      1.9838860978d+00
      s( 7) = -.5635914404d+00
      e( 8) =      0.6483632394d+00
      s( 8) = -.5195400905d+00
      e( 9) =      0.0524439674d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.0280481607d+00
      s(10) = 0.1000000000d+01
      e(11) =     75.4018620170d+00
      p(11) = 0.1543536253d-01
      e(12) =     17.2748189780d+00
      p(12) = 0.9973829318d-01
      e(13) =      5.1842347425d+00
      p(13) = 0.3120959397d+00
      e(14) =      1.6601211973d+00
      p(14) = 0.4929567481d+00
      e(15) =      0.5123252896d+00
      p(15) = 0.3242039832d+00
c
c additional p-function
c Ref.: R. Ahlrichs 12.08.92, as steepr s-GTO of SV.
c which is optimal for NanClm
c
      e(16) = .052d+00
      p(16) =  1.00d+00
c
      return
c
c     ----- mg SV -----
c
c mg    (10s5p) / [4s1p]     {5311/5}
c scf energy is     -199.53359311668 a.u.
c 
  120 continue
      e( 1) =   4953.8339196000d+00
      s( 1) = -.5777896750d-02
      e( 2) =    745.1804415400d+00
      s( 2) = -.4312476108d-01
      e( 3) =    169.2160497200d+00
      s( 3) = -.1926821699d+00
      e( 4) =     47.3006720190d+00
      s( 4) = -.4864143912d+00
      e( 5) =     14.4613369730d+00
      s( 5) = -.4255089408d+00
      e( 6) =     24.7681747890d+00
      s( 6) = 0.8795696998d-01
      e( 7) =      2.4940945349d+00
      s( 7) = -.5516505813d+00
      e( 8) =      0.8780758453d+00
      s( 8) = -.5344329483d+00
      e( 9) =      0.0872127825d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.0335992938d+00
      s(10) = 0.1000000000d+01
      e(11) =     98.0530104940d+00
      p(11) = -.1448056460d-01
      e(12) =     22.5869322770d+00
      p(12) = -.9549575079d-01
      e(13) =      6.8391509842d+00
      p(13) = -.3078767265d+00
      e(14) =      2.2332843818d+00
      p(14) = -.4993629289d+00
      e(15) =      0.7160659939d+00
      p(15) = -.3150347621d+00
c
c additional p-function
c taken from SV.3p basis set, the steeper p-GTO
c
      e(16) = .18914796195d+00
      p(16) =  1.00d+00
c
      return
c
c     ----- al SV -----
c
c al    (10s7p) / [4s3p]     {5311/511}
c scf energy is     -241.78686952827 a.u. 
c 
  140 continue
      e( 1) =   5887.5727030000d+00
      s( 1) = 0.1348334799d-02
      e( 2) =    885.6122599600d+00
      s( 2) = 0.1007157681d-01
      e( 3) =    201.1360489900d+00
      s( 3) = 0.4513245406d-01
      e( 4) =     56.2849746740d+00
      s( 4) = 0.1146126804d+00
      e( 5) =     17.2295512430d+00
      s( 5) = 0.1015960894d+00
      e( 6) =     29.3402499220d+00
      s( 6) = 0.6934745421d-01
      e( 7) =      3.0439630420d+00
      s( 7) = -.4252811768d+00
      e( 8) =      1.1285539518d+00
      s( 8) = -.4144983221d+00
      e( 9) =      0.1423417516d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.0544001923d+00
      s(10) = 0.1000000000d+01
      e(11) =    145.1191880900d+00
      p(11) = 0.6396337313d-02
      e(12) =     33.7178948330d+00
      p(12) = 0.4418935997d-01
      e(13) =     10.3698630830d+00
      p(13) = 0.1558157599d+00
      e(14) =      3.5135616036d+00
      p(14) = 0.2863528695d+00
      e(15) =      1.1980050273d+00
      p(15) = 0.2292142325d+00
      e(16) =      0.2658300591d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.0710033620d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(18) = .30d+00
      d(18) =  1.00d+00
c
      return
c
c     ----- si SV -----
c
c si    (10s7p) / [4s3p]     {5311/511}
c scf energy is     -288.74978456230 a.u.
c
  160 continue
      e( 1) =   6903.7118686000d+00
      s( 1) = 0.1337396300d-02
      e( 2) =   1038.4346419000d+00
      s( 2) = 0.9996654624d-02
      e( 3) =    235.8758148000d+00
      s( 3) = 0.4491016510d-01
      e( 4) =     66.0693851690d+00
      s( 4) = 0.1146363854d+00
      e( 5) =     20.2479457610d+00
      s( 5) = 0.1028006386d+00
      e( 6) =     34.3534817300d+00
      s( 6) = 0.7083728501d-01
      e( 7) =      3.6370788192d+00
      s( 7) = -.4302883625d+00
      e( 8) =      1.4002048599d+00
      s( 8) = -.4138277497d+00
      e( 9) =      0.2048441481d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.0779940955d+00
      s(10) = 0.1000000000d+01
      e(11) =    179.8390737300d+00
      p(11) = 0.6191665646d-02
      e(12) =     41.9072588460d+00
      p(12) = 0.4339943198d-01
      e(13) =     12.9552943670d+00
      p(13) = 0.1563201935d+00
      e(14) =      4.4383267393d+00
      p(14) = 0.2941999698d+00
      e(15) =      1.5462247904d+00
      p(15) = 0.2353682381d+00
      e(16) =      0.3560761230d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.1000851376d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(18) = .35d+00
      d(18) =  1.00d+00
c
      return
c
c     ----- p SV -----
c
c p     (10s7p) / [4s3p]     {5311/511}
c scf energy is     -340.59736844515 a.u. 
c
  180 continue
      e( 1) =   8002.4795106000d+00
      s( 1) = 0.5750348930d-02
      e( 2) =   1203.6813590000d+00
      s( 2) = 0.4300762896d-01
      e( 3) =    273.4422703100d+00
      s( 3) = 0.1936398625d+00
      e( 4) =     76.6555415170d+00
      s( 4) = 0.4965169340d+00
      e( 5) =     23.5169274350d+00
      s( 5) = 0.4498326248d+00
      e( 6) =     39.7916834390d+00
      s( 6) = 0.9518812979d-01
      e( 7) =      4.2770343323d+00
      s( 7) = -.5764984037d+00
      e( 8) =      1.6940256888d+00
      s( 8) = -.5423958387d+00
      e( 9) =      0.2756767464d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.1049559010d+00
      s(10) = 0.1000000000d+01
      e(11) =    219.5075582300d+00
      p(11) = 0.9210056526d-02
      e(12) =     51.2741550300d+00
      p(12) = 0.6540976575d-01
      e(13) =     15.9215958920d+00
      p(13) = 0.2403373028d+00
      e(14) =      5.5069913481d+00
      p(14) = 0.4631832179d+00
      e(15) =      1.9537719426d+00
      p(15) = 0.3739256338d+00
      e(16) =      0.4780339793d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.1365795262d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(18) = .45d+00
      d(18) = 1.00d+00
      return
c
c     ----- s SV  -----
c
c s     (10s7p) / [4s3p]     {5311/511}
c scf energy is     -397.36223459768 a.u. 
c
  200 continue
      e( 1) =   9184.9303010000d+00
      s( 1) = -.2229438776d-02
      e( 2) =   1381.5105503000d+00
      s( 2) = -.1668302994d-01
      e( 3) =    313.8714758000d+00
      s( 3) = -.7526243612d-01
      e( 4) =     88.0538706230d+00
      s( 4) = -.1937682704d+00
      e( 5) =     27.0399149050d+00
      s( 5) = -.1771802080d+00
      e( 6) =     45.6487313030d+00
      s( 6) = -.1073606257d+00
      e( 7) =      4.9664522326d+00
      s( 7) = 0.6506629302d+00
      e( 8) =      2.0116242047d+00
      s( 8) = 0.5971215535d+00
      e( 9) =      0.3566107701d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.1350722148d+00
      s(10) = 0.1000000000d+01
      e(11) =    261.9823343900d+00
      p(11) = -.9272992982d-02
      e(12) =     61.3068947360d+00
      p(12) = -.6654766924d-01
      e(13) =     19.1037298870d+00
      p(13) = -.2482859590d+00
      e(14) =      6.6567720378d+00
      p(14) = -.4870384740d+00
      e(15) =      2.3959635161d+00
      p(15) = -.3933785031d+00
      e(16) =      0.6177616168d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.1699337687d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(18) = .55d+00
      d(18) = 1.00d+00
      return
c
c     ----- cl SV -----
c
c cl    (10s7p) / [4s3p]     {5311/511}
c SCF energy is     -459.31519371363 a.u.  
c
  220 continue
      e( 1) =  10449.8275660000d+00
      s( 1) = 0.1970836248d-02
      e( 2) =   1571.7365221000d+00
      s( 2) = 0.1475472798d-01
      e( 3) =    357.1206552300d+00
      s( 3) = 0.6667911288d-01
      e( 4) =    100.2518593500d+00
      s( 4) = 0.1722892408d+00
      e( 5) =     30.8127275540d+00
      s( 5) = 0.1588378610d+00
      e( 6) =     51.9237894340d+00
      s( 6) = -.1000929891d+00
      e( 7) =      5.7045760975d+00
      s( 7) = 0.6084175275d+00
      e( 8) =      2.3508376809d+00
      s( 8) = 0.5435215336d+00
      e( 9) =      0.4460512467d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.1684885619d+00
      s(10) = 0.1000000000d+01
      e(11) =    307.6679056900d+00
      p(11) = -.8780148412d-02
      e(12) =     72.1020155150d+00
      p(12) = -.6356335547d-01
      e(13) =     22.5326802620d+00
      p(13) = -.2401642828d+00
      e(14) =      7.8991765444d+00
      p(14) = -.4779886656d+00
      e(15) =      2.8767268321d+00
      p(15) = -.3851585001d+00
      e(16) =      0.7745936396d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.2103769970d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, 
c depending on atomic charge, values down to 0.2 result in molecular
c optimization, e.g. for NanClm.
c
      e(18) =   .65d+00
      d(18) =  1.00d+00
c
      return
c
c     ----- ar SV -----
c
c ar    (10s7p) / [4s3p]     {5311/511}
c scf energy is     -526.6233849100 a.u. 
c
  240 continue
      e( 1) =  11797.1197650000d+00
      s( 1) = 0.2021447997d-02
      e( 2) =   1774.3522753000d+00
      s( 2) = 0.1513985269d-01
      e( 3) =    403.1887573300d+00
      s( 3) = 0.6852540087d-01
      e( 4) =    113.2493399900d+00
      s( 4) = 0.1776292921d+00
      e( 5) =     34.8352981930d+00
      s( 5) = 0.1649649500d+00
      e( 6) =     58.6147750420d+00
      s( 6) = -.1034339400d+00
      e( 7) =      6.4922045078d+00
      s( 7) = 0.6313336590d+00
      e( 8) =      2.7117014403d+00
      s( 8) = 0.5488757230d+00
      e( 9) =      0.5441297454d+00
      s( 9) = 0.1000000000d+01
      e(10) =      0.2051741114d+00
      s(10) = 0.1000000000d+01
      e(11) =    356.2870725600d+00
      p(11) = -.8732178327d-02
      e(12) =     83.5931328660d+00
      p(12) = -.6368031744d-01
      e(13) =     26.1867040290d+00
      p(13) = -.2431190689d+00
      e(14) =      9.2257275925d+00
      p(14) = -.4895606965d+00
      e(15) =      3.3922754954d+00
      p(15) = -.3922919004d+00
      e(16) =      0.9474053401d+00
      p(16) = 0.1000000000d+01
      e(17) =      0.2565913506d+00
      p(17) = 0.1000000000d+01
c
c additional d-GTO
c taken from tzp basis set in hondo basis set library
c
      e(18) =  .696d+00
      d(18) =  1.00d+00
      return
c
      end
      subroutine ahsv3(e,s,p,d,n)
c
c ----- ahlrich's SV (14s8p5d/5s2p2d) {63311/53/41}contractions  ---
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-18
      go to (100,120,140,160,180,200,220,240,260,280,300,320),nn
c
c k     (14s8p) / [5s2p]     {63311/53}
c optimized for atomic ground state k(2s)
c scf energy is     -599.0769622789 a.u.
c
  100 continue
      e( 1) =  31478.7467640000d+00
      s( 1) = 0.3983865399d-02
      e( 2) =   4726.8876066000d+00
      s( 2) = 0.3050175976d-01
      e( 3) =   1075.4345353000d+00
      s( 3) = 0.1507375262d+00
      e( 4) =    303.3981102300d+00
      s( 4) = 0.5191293980d+00
      e( 5) =     98.3271128310d+00
      s( 5) = 0.1036695701d+01
      e( 6) =     33.6362221770d+00
      s( 6) = 0.7639896320d+00
      e( 7) =     65.6392099620d+00
      s( 7) = -.2824261711d+00
      e( 8) =      7.3162592218d+00
      s( 8) = 0.1691493586d+01
      e( 9) =      2.8902580135d+00
      s( 9) = 0.1296533195d+01
      e(10) =      4.5459748965d+00
      s(10) = -.7634355527d-02
      e(11) =      0.7040412406d+00
      s(11) = 0.2563571896d-01
      e(12) =      0.2826688896d+00
      s(12) = 0.1660685921d-01
      e(13) =      0.0290581640d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0121116382d+00
      s(14) = 0.1000000000d+01
      e(15) =    361.2249215400d+00
      p(15) = 0.2090647982d-01
      e(16) =     84.6702221660d+00
      p(16) = 0.1504364174d+00
      e(17) =     26.4690882360d+00
      p(17) = 0.5544006108d+00
      e(18) =      9.2658077615d+00
      p(18) = 0.1040900999d+01
      e(19) =      3.3423388293d+00
      p(19) = 0.6782534119d+00
      e(20) =      1.5100876104d+00
      p(20) = 0.7524819115d+00
      e(21) =      0.5656837516d+00
      p(21) = 0.1370858503d+01
      e(22) =      0.2081700850d+00
      p(22) = 0.6604763308d+00
c
c additional p-GTO, steeper one of DP
c Wachters, JCP, 52 (1970), 1033
c
      e(23) =   .041737d+00
      p(23) =  1.000000d+00
c
      return
c
c ca    (14s8p) / [5s2p]     {63311/53}
c optimized for atomic ground state ca(1s)
c scf energy is     -676.6600696944 a.u.
c 
  120 continue
      e( 1) =  35138.7139290000d+00
      s( 1) = 0.3948252074d-02
      e( 2) =   5276.4111348000d+00
      s( 2) = 0.3023424355d-01
      e( 3) =   1200.4692589000d+00
      s( 3) = 0.1495201968d+00
      e( 4) =    338.7181054200d+00
      s( 4) = 0.5159734571d+00
      e( 5) =    109.8538592200d+00
      s( 5) = 0.1033951030d+01
      e( 6) =     37.6088802990d+00
      s( 6) = 0.7693793353d+00
      e( 7) =     73.1079775550d+00
      s( 7) = -.2826852501d+00
      e( 8) =      8.2407705688d+00
      s( 8) = 0.1679609214d+01
      e( 9) =      3.2959812993d+00
      s( 9) = 0.1280376602d+01
      e(10) =      5.2341800914d+00
      s(10) = -.7686860456d-02
      e(11) =      0.8418722052d+00
      s(11) = 0.2538237598d-01
      e(12) =      0.3651029403d+00
      s(12) = 0.1651217151d-01
      e(13) =      0.0512224029d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0198251114d+00
      s(14) = 0.1000000000d+01
      e(15) =    413.1131389300d+00
      p(15) = 0.2032713535d-01
      e(16) =     96.9357862240d+00
      p(16) = 0.1473027636d+00
      e(17) =     30.3721546590d+00
      p(17) = 0.5488716732d+00
      e(18) =     10.6847768300d+00
      p(18) = 0.1044065982d+01
      e(19) =      3.8821258350d+00
      p(19) = 0.6865349068d+00
      e(20) =      1.7993016295d+00
      p(20) = 0.7541024687d+00
      e(21) =      0.6918905653d+00
      p(21) = 0.1340929660d+01
      e(22) =      0.2636402410d+00
      p(22) = 0.5639198944d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(23) =   .074979d+00
      p(23) =  1.000000d+00
c
      return
c
c sc    (14s8p5d) / [5s2p2d]     {63311/53/41}
c optimized for atomic ground state sc(2d)
c scf energy is     -759.6252351669 a.u.
c 
  140 continue
      e( 1) =  38956.0818040000d+00
      s( 1) = 0.3929320586d-02
      e( 2) =   5849.5733637000d+00
      s( 2) = 0.3009322123d-01
      e( 3) =   1330.8813154000d+00
      s( 3) = 0.1489037095d+00
      e( 4) =    375.5553416500d+00
      s( 4) = 0.5146428291d+00
      e( 5) =    121.8726137000d+00
      s( 5) = 0.1033770807d+01
      e( 6) =     41.7602437290d+00
      s( 6) = 0.7743685163d+00
      e( 7) =     81.0606339530d+00
      s( 7) = -.2831855232d+00
      e( 8) =      9.2059823972d+00
      s( 8) = 0.1677080698d+01
      e( 9) =      3.7063215732d+00
      s( 9) = 0.1259473368d+01
      e(10) =      5.9888909988d+00
      s(10) = -.7782136749d-02
      e(11) =      0.9736343238d+00
      s(11) = 0.2549969275d-01
      e(12) =      0.4204101922d+00
      s(12) = 0.1619156082d-01
      e(13) =      0.0594405519d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0228978063d+00
      s(14) = 0.1000000000d+01
      e(15) =    466.3148126200d+00
      p(15) = 0.1998430039d-01
      e(16) =    109.5121709700d+00
      p(16) = 0.1456104307d+00
      e(17) =     34.3759218270d+00
      p(17) = 0.5468746622d+00
      e(18) =     12.1420969550d+00
      p(18) = 0.1047900601d+01
      e(19) =      4.4336767669d+00
      p(19) = 0.6889489033d+00
      e(20) =      2.0971291866d+00
      p(20) = 0.7561921472d+00
      e(21) =      0.8097760696d+00
      p(21) = 0.1317821224d+01
      e(22) =      0.3083404659d+00
      p(22) = 0.5431226817d+00
      e(23) =     19.2403349280d+00
      d(23) = 0.2703908214d-01
      e(24) =      5.1178995899d+00
      d(24) = 0.1380368474d+00
      e(25) =      1.6554278827d+00
      d(25) = 0.3486908640d+00
      e(26) =      0.5401635561d+00
      d(26) = 0.4859418572d+00
      e(27) =      0.1621121452d+00
      d(27) = 0.3437444969d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .089748d+00
      p(28) = 1.000000d+00
c
      return
c
c ti    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is     -848.2820957739 a.u.
c optimized for atomic ground state ti(3f)
c calc. performed in oh-symmetry with 3d-occupation eg**2;
c 
  160 continue
      e( 1) =  42961.5121850000d+00
      s( 1) = 0.3912763536d-02
      e( 2) =   6450.9759169000d+00
      s( 2) = 0.2996982049d-01
      e( 3) =   1467.7210915000d+00
      s( 3) = 0.1483635271d+00
      e( 4) =    414.2099735500d+00
      s( 4) = 0.5134728532d+00
      e( 5) =    134.4871584000d+00
      s( 5) = 0.1033536548d+01
      e( 6) =     46.1222097960d+00
      s( 6) = 0.7785423393d+00
      e( 7) =     89.4477625430d+00
      s( 7) = -.2838540126d+00
      e( 8) =     10.2233460600d+00
      s( 8) = 0.1677278533d+01
      e( 9) =      4.1353774271d+00
      s( 9) = 0.1241192846d+01
      e(10) =      6.7896181452d+00
      s(10) = -.7839999452d-02
      e(11) =      1.1106730691d+00
      s(11) = 0.2549549302d-01
      e(12) =      0.4756597558d+00
      s(12) = 0.1606117289d-01
      e(13) =      0.0659869569d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0252103423d+00
      s(14) = 0.1000000000d+01
      e(15) =    522.0368478200d+00
      p(15) = 0.1975417964d-01
      e(16) =    122.6864948900d+00
      p(16) = 0.1446067762d+00
      e(17) =     38.5729036110d+00
      p(17) = 0.5466900417d+00
      e(18) =     13.6721693190d+00
      p(18) = 0.1053164754d+01
      e(19) =      5.0118529359d+00
      p(19) = 0.6911121336d+00
      e(20) =      2.4131928282d+00
      p(20) = 0.7580343714d+00
      e(21) =      0.9325227005d+00
      p(21) = 0.1303624140d+01
      e(22) =      0.3542905839d+00
      p(22) = 0.5363865330d+00
      e(23) =     23.4651259570d+00
      d(23) = 0.2653638012d-01
      e(24) =      6.3332593832d+00
      d(24) = 0.1379645396d+00
      e(25) =      2.0766489946d+00
      d(25) = 0.3531264423d+00
      e(26) =      0.6902736195d+00
      d(26) = 0.4864712417d+00
      e(27) =      0.2108873855d+00
      d(27) = 0.3302631426d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .101561d+00
      p(28) = 1.000000d+00
c
      return
c
c v    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is     -942.7457186340 a.u. 
c optimized for atomic ground state v(4f)
c calculation performed in oh-symmetry with 3d-occupation t2g**3;
c 
  180 continue
      e( 1) =  47160.3760600000d+00
      s( 1) = 0.1449868891d-02
      e( 2) =   7081.4110871000d+00
      s( 2) = 0.1110643525d-01
      e( 3) =   1611.1621223000d+00
      s( 3) = 0.5500542359d-01
      e( 4) =    454.7294055100d+00
      s( 4) = 0.1906025259d+00
      e( 5) =    147.7132120800d+00
      s( 5) = 0.3843502296d+00
      e( 6) =     50.6995389500d+00
      s( 6) = 0.2909554679d+00
      e( 7) =     98.2624926690d+00
      s( 7) = -.1094233786d+00
      e( 8) =     11.2942930990d+00
      s( 8) = 0.6453949040d+00
      e( 9) =      4.5853360105d+00
      s( 9) = 0.4711788078d+00
      e(10) =      7.6359689588d+00
      s(10) = -.2245494905d+00
      e(11) =      1.2539836689d+00
      s(11) = 0.7259485276d+00
      e(12) =      0.5327193539d+00
      s(12) = 0.4556058275d+00
      e(13) =      0.0722462396d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0273580874d+00
      s(14) = 0.1000000000d+01
      e(15) =    580.5504498800d+00
      p(15) = 0.9731511092d-02
      e(16) =    136.5234112700d+00
      p(16) = 0.7153124114d-01
      e(17) =     42.9839588200d+00
      p(17) = 0.2719768841d+00
      e(18) =     15.2827987630d+00
      p(18) = 0.5261898889d+00
      e(19) =      5.6202495154d+00
      p(19) = 0.3445253350d+00
      e(20) =      2.7485386415d+00
      p(20) = 0.3404039650d+00
      e(21) =      1.0618550073d+00
      p(21) = 0.5798399612d+00
      e(22) =      0.4023551865d+00
      p(22) = 0.2391164308d+00
      e(23) =     27.3584340170d+00
      d(23) = 0.2664192705d-01
      e(24) =      7.4540604253d+00
      d(24) = 0.1399531173d+00
      e(25) =      2.4633917847d+00
      d(25) = 0.3575106664d+00
      e(26) =      0.8248092528d+00
      d(26) = 0.4848835415d+00
      e(27) =      0.2525790474d+00
      d(27) = 0.3233284500d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .111248d+00
      p(28) = 1.000000d+00
c
      return
c
c cr    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1043.1960697534 a.u. 
c optimized for atomic ground state cr(7s)
c
  200 continue
      e( 1) =  51528.0863490000d+00
      s( 1) = 0.1440582311d-02
      e( 2) =   7737.2103487000d+00
      s( 2) = 0.1103620229d-01
      e( 3) =   1760.3748470000d+00
      s( 3) = 0.5467665181d-01
      e( 4) =    496.8770654400d+00
      s( 4) = 0.1896503810d+00
      e( 5) =    161.4652059800d+00
      s( 5) = 0.3829541285d+00
      e( 6) =     55.4663522680d+00
      s( 6) = 0.2909005067d+00
      e( 7) =    107.5473299900d+00
      s( 7) = -.1093228110d+00
      e( 8) =     12.4086718970d+00
      s( 8) = 0.6447259947d+00
      e( 9) =      5.0423628826d+00
      s( 9) = 0.4626271256d+00
      e(10) =      8.5461640165d+00
      s(10) = -.2271101329d+00
      e(11) =      1.3900441221d+00
      s(11) = 0.7330152759d+00
      e(12) =      0.5606660288d+00
      s(12) = 0.4422556543d+00
      e(13) =      0.0714837060d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0282506876d+00
      s(14) = 0.1000000000d+01
      e(15) =    640.4853609600d+00
      p(15) = 0.9612671520d-02
      e(16) =    150.6971119400d+00
      p(16) = 0.7088983466d-01
      e(17) =     47.5037552960d+00
      p(17) = 0.2706525899d+00
      e(18) =     16.9341201650d+00
      p(18) = 0.5243734341d+00
      e(19) =      6.2409680590d+00
      p(19) = 0.3410799471d+00
      e(20) =      3.0885463206d+00
      p(20) = 0.3397398690d+00
      e(21) =      1.1791047769d+00
      p(21) = 0.5727206293d+00
      e(22) =      0.4336977443d+00
      p(22) = 0.2458272821d+00
      e(23) =     27.5594794260d+00
      d(23) = 0.3061248804d-01
      e(24) =      7.4687020327d+00
      d(24) = 0.1559327094d+00
      e(25) =      2.4345903574d+00
      d(25) = 0.3698442128d+00
      e(26) =      0.7824475481d+00
      d(26) = 0.4707111808d+00
      e(27) =      0.2199577431d+00
      d(27) = 0.3394164989d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .120675d+00
      p(28) = 1.000000d+00
c
      return
c
c mn    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1149.6942409653 a.u.
c optimized for atomic ground state mn(6s)
c derived from  mn 13s8p5d by adding  steep s functn to '5s' bl
c 
  220 continue
      e( 1) =  56137.0090370000d+00
      s( 1) = 0.1432130470d-02
      e( 2) =   8429.2063943000d+00
      s( 2) = 0.1097250916d-01
      e( 3) =   1917.8277233000d+00
      s( 3) = 0.5438246871d-01
      e( 4) =    541.3623019800d+00
      s( 4) = 0.1888433513d+00
      e( 5) =    176.0006914200d+00
      s( 5) = 0.3819802505d+00
      e( 6) =     60.5004770100d+00
      s( 6) = 0.2915677260d+00
      e( 7) =    117.1728288200d+00
      s( 7) = -.1093366133d+00
      e( 8) =     13.5969733680d+00
      s( 8) = 0.6430503943d+00
      e( 9) =      5.5483996341d+00
      s( 9) = 0.4584897058d+00
      e(10) =      9.4662853530d+00
      s(10) = -.2253897726d+00
      e(11) =      1.5595006070d+00
      s(11) = 0.7230775866d+00
      e(12) =      0.6523020587d+00
      s(12) = 0.4530072154d+00
      e(13) =      0.0840037345d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0312560986d+00
      s(14) = 0.1000000000d+01
      e(15) =    706.0049753500d+00
      p(15) = 0.9505551817d-02
      e(16) =    166.1972882000d+00
      p(16) = 0.7035627114d-01
      e(17) =     52.4520619060d+00
      p(17) = 0.2700555698d+00
      e(18) =     18.7469328620d+00
      p(18) = 0.5257434460d+00
      e(19) =      6.9282991622d+00
      p(19) = 0.3425403322d+00
      e(20) =      3.4772204938d+00
      p(20) = 0.3399407374d+00
      e(21) =      1.3406906449d+00
      p(21) = 0.5720383625d+00
      e(22) =      0.5049880304d+00
      p(22) = 0.2384760583d+00
      e(23) =     35.4232649350d+00
      d(23) = 0.2698530411d-01
      e(24) =      9.7814221451d+00
      d(24) = 0.1438345865d+00
      e(25) =      3.2673488767d+00
      d(25) = 0.3641895838d+00
      e(26) =      1.1026472189d+00
      d(26) = 0.4815267066d+00
      e(27) =      0.3374320593d+00
      d(27) = 0.3145875436d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .127650d+00
      p(28) = 1.000000d+00
c
      return
c
c fe    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1262.2517936909 a.u. 
c optimized for atomic ground state fe(5d)
c derived from huzinaga fe 14s8p5d 5333/53/5 by replacing
c the last '3s'-contn by '1s1s' taken from 11s6p3d 4322/42/3
c and full optimization; '5d' now contracted to 41 
c 
  240 continue
      e( 1) =  60923.6406430000d+00
      s( 1) = 0.1430225447d-02
      e( 2) =   9147.8893982000d+00
      s( 2) = 0.1095879004d-01
      e( 3) =   2081.3505927000d+00
      s( 3) = 0.5433255425d-01
      e( 4) =    587.5597706700d+00
      s( 4) = 0.1888499501d+00
      e( 5) =    191.0904399000d+00
      s( 5) = 0.3825306995d+00
      e( 6) =     65.7327301120d+00
      s( 6) = 0.2930833598d+00
      e( 7) =    127.2589192800d+00
      s( 7) = -.1096456493d+00
      e( 8) =     14.8309130100d+00
      s( 8) = 0.6438763133d+00
      e( 9) =      6.0653307408d+00
      s( 9) = 0.4547234732d+00
      e(10) =     10.4499437100d+00
      s(10) = -.2253963995d+00
      e(11) =      1.7245228003d+00
      s(11) = 0.7216439816d+00
      e(12) =      0.7177217733d+00
      s(12) = 0.4498549292d+00
      e(13) =      0.0914498283d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0337066910d+00
      s(14) = 0.1000000000d+01
      e(15) =    773.4375099500d+00
      p(15) = 0.9432573514d-02
      e(16) =    182.1514971400d+00
      p(16) = 0.7002962058d-01
      e(17) =     57.5472727580d+00
      p(17) = 0.2699365200d+00
      e(18) =     20.6149889350d+00
      p(18) = 0.5270001105d+00
      e(19) =      7.6348557890d+00
      p(19) = 0.3428414803d+00
      e(20) =      3.8719327990d+00
      p(20) = 0.3397440299d+00
      e(21) =      1.4924724132d+00
      p(21) = 0.5684259401d+00
      e(22) =      0.5606128496d+00
      p(22) = 0.2364936584d+00
      e(23) =     38.9681334190d+00
      d(23) = 0.2787966438d-01
      e(24) =     10.8000670780d+00
      d(24) = 0.1485831998d+00
      e(25) =      3.6136457999d+00
      d(25) = 0.3690547950d+00
      e(26) =      1.2129967888d+00
      d(26) = 0.4774510088d+00
      e(27) =      0.3652439317d+00
      d(27) = 0.3141814230d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .134915d+00
      p(28) = 1.000000d+00
c
      return
c
c co    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1381.2013463555 a.u. 
c optimized for atomic ground state co(4f)
c calculation performed in Oh-symm with 3d-occn eg**4 t2g**3;
c 
  260 continue
      e( 1) =  65902.2082570000d+00
      s( 1) = 0.1428461494d-02
      e( 2) =   9895.3896027000d+00
      s( 2) = 0.1094607278d-01
      e( 3) =   2251.4305789000d+00
      s( 3) = 0.5428595389d-01
      e( 4) =    635.6109708400d+00
      s( 4) = 0.1888517908d+00
      e( 5) =    206.7882068100d+00
      s( 5) = 0.3830163499d+00
      e( 6) =     71.1792429710d+00
      s( 6) = 0.2944355127d+00
      e( 7) =    137.7726804000d+00
      s( 7) = -.1099022174d+00
      e( 8) =     16.1180792430d+00
      s( 8) = 0.6445553740d+00
      e( 9) =      6.6030327710d+00
      s( 9) = 0.4511678792d+00
      e(10) =     11.4799157880d+00
      s(10) = -.2259384691d+00
      e(11) =      1.8956426324d+00
      s(11) = 0.7223140901d+00
      e(12) =      0.7846623207d+00
      s(12) = 0.4490381230d+00
      e(13) =      0.0984257744d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0359457419d+00
      s(14) = 0.1000000000d+01
      e(15) =    843.6435857500d+00
      p(15) = 0.9386609725d-02
      e(16) =    198.7638699400d+00
      p(16) = 0.6988020872d-01
      e(17) =     62.8549630980d+00
      p(17) = 0.2703707035d+00
      e(18) =     22.5628422800d+00
      p(18) = 0.5290478688d+00
      e(19) =      8.3713209127d+00
      p(19) = 0.3435702958d+00
      e(20) =      4.2858719800d+00
      p(20) = 0.3402799904d+00
      e(21) =      1.6508041817d+00
      p(21) = 0.5669339238d+00
      e(22) =      0.6183423110d+00
      p(22) = 0.2361797978d+00
      e(23) =     42.9278676120d+00
      d(23) = 0.2848778837d-01
      e(24) =     11.9425330530d+00
      d(24) = 0.1520695128d+00
      e(25) =      4.0046495664d+00
      d(25) = 0.3731091400d+00
      e(26) =      1.3413193804d+00
      d(26) = 0.4754983768d+00
      e(27) =      0.4001500974d+00
      d(27) = 0.3134683142d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =  .141308d+00
      p(28) = 1.000000d+00
c
      return
c
c ni    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1506.6346013010 a.u.
c optimized for atomic ground state ni(3f)
c calculation performed in Oh-symm with 3d-occn t2g**6 eg**2;
c 
  280 continue
      e( 1) =  71074.8032110000d+00
      s( 1) = 0.1426038673d-02
      e( 2) =  10672.0209410000d+00
      s( 2) = 0.1092823699d-01
      e( 3) =   2428.1389007000d+00
      s( 3) = 0.5421262694d-01
      e( 4) =    685.5359514800d+00
      s( 4) = 0.1887476890d+00
      e( 5) =    223.1007286300d+00
      s( 5) = 0.3832461699d+00
      e( 6) =     76.8420140420d+00
      s( 6) = 0.2955063714d+00
      e( 7) =    148.7112201600d+00
      s( 7) = -.1101444306d+00
      e( 8) =     17.4591549870d+00
      s( 8) = 0.6452142699d+00
      e( 9) =      7.1625280665d+00
      s( 9) = 0.4479783810d+00
      e(10) =     12.5561371250d+00
      s(10) = -.2264540322d+00
      e(11) =      2.0735740488d+00
      s(11) = 0.7232095929d+00
      e(12) =      0.8538264060d+00
      s(12) = 0.4486802648d+00
      e(13) =      0.1053676627d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0381340877d+00
      s(14) = 0.1000000000d+01
      e(15) =    916.7360866200d+00
      p(15) = 0.9343963561d-02
      e(16) =    216.0613991300d+00
      p(16) = 0.6973737490d-01
      e(17) =     68.3839148170d+00
      p(17) = 0.2707349501d+00
      e(18) =     24.5938439520d+00
      p(18) = 0.5307830155d+00
      e(19) =      9.1392960204d+00
      p(19) = 0.3441022944d+00
      e(20) =      4.7193371746d+00
      p(20) = 0.3407608202d+00
      e(21) =      1.8161849234d+00
      p(21) = 0.5658016961d+00
      e(22) =      0.6784075072d+00
      p(22) = 0.2361671736d+00
      e(23) =     47.0938321080d+00
      d(23) = 0.2898231695d-01
      e(24) =     13.1464639750d+00
      d(24) = 0.1549499595d+00
      e(25) =      4.4170548925d+00
      d(25) = 0.3763311511d+00
      e(26) =      1.4771565078d+00
      d(26) = 0.4736509601d+00
      e(27) =      0.4373592179d+00
      d(27) = 0.3124783783d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =   .146588d+00
      p(28) =  1.000000d+00
c
      return
c
c cu    (14s8p5d) / [5s2p2d]     {63311/53/41}
c scf energy is    -1638.6871046889 a.u. 
c optimized for atomic state cu(2s)
c
  300 continue
      e( 1) =  76381.3480560000d+00
      s( 1) = 0.1433607990d-02
      e( 2) =  11468.7774990000d+00
      s( 2) = 0.1098674987d-01
      e( 3) =   2609.4246495000d+00
      s( 3) = 0.5451365247d-01
      e( 4) =    736.7503309800d+00
      s( 4) = 0.1899012826d+00
      e( 5) =    239.8241995800d+00
      s( 5) = 0.3858195921d+00
      e( 6) =     82.6568292520d+00
      s( 6) = 0.2979060750d+00
      e( 7) =    160.1354419600d+00
      s( 7) = -.1114677857d+00
      e( 8) =     18.8341776950d+00
      s( 8) = 0.6534930103d+00
      e( 9) =      7.7176595741d+00
      s( 9) = 0.4477053442d+00
      e(10) =     13.7108467170d+00
      s(10) = -.2287091112d+00
      e(11) =      2.2349895670d+00
      s(11) = 0.7346442303d+00
      e(12) =      0.8781836007d+00
      s(12) = 0.4327307087d+00
      e(13) =      0.0871874581d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0329691147d+00
      s(14) = 0.1000000000d+01
      e(15) =    991.2407578200d+00
      p(15) = 0.9387849880d-02
      e(16) =    233.6937611600d+00
      p(16) = 0.7020828246d-01
      e(17) =     74.0209309270d+00
      p(17) = 0.2732352222d+00
      e(18) =     26.6649674470d+00
      p(18) = 0.5358079273d+00
      e(19) =      9.9192087478d+00
      p(19) = 0.3457579491d+00
      e(20) =      5.1519553926d+00
      p(20) = 0.3422910808d+00
      e(21) =      1.9638205828d+00
      p(21) = 0.5645659248d+00
      e(22) =      0.7156009704d+00
      p(22) = 0.2407858432d+00
      e(23) =     47.3350495900d+00
      d(23) = 0.3237554776d-01
      e(24) =     13.1616660770d+00
      d(24) = 0.1681021868d+00
      e(25) =      4.3693777244d+00
      d(25) = 0.3847770798d+00
      e(26) =      1.4132925109d+00
      d(26) = 0.4614788018d+00
      e(27) =      0.3887800145d+00
      d(27) = 0.3238887326d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033  (2P)
c
      e(28) =  .155065d+00
      p(28) = 1.000000d+00
c
      return
c
c zn    (13s8p5d) / [5s2p2d]     {53311/53/41}
c scf energy is    -1777.5602821596 a.u. 
c optimized for atomic ground state zn(1s) 
c 
  320 continue
      e( 1) =  82000.7116290000d+00
      s( 1) = 0.1421076400d-02
      e( 2) =  12312.4717770000d+00
      s( 2) = 0.1089149949d-01
      e( 3) =   2801.3944193000d+00
      s( 3) = 0.5405718806d-01
      e( 4) =    790.9942430200d+00
      s( 4) = 0.1884746390d+00
      e( 5) =    257.5655107900d+00
      s( 5) = 0.3834654935d+00
      e( 6) =     88.8149334000d+00
      s( 6) = 0.2972379420d+00
      e( 7) =    171.8635371600d+00
      s( 7) = -.1105184952d+00
      e( 8) =     20.3025347850d+00
      s( 8) = 0.6460771698d+00
      e( 9) =      8.3464123068d+00
      s( 9) = 0.4422011732d+00
      e(10) =     14.8475369400d+00
      s(10) = -.2270530928d+00
      e(11) =      2.4495029507d+00
      s(11) = 0.7243321794d+00
      e(12) =      0.9984582182d+00
      s(12) = 0.4483649559d+00
      e(13) =      0.1189130794d+00
      s(13) = 0.1000000000d+01
      e(14) =      0.0422974288d+00
      s(14) = 0.1000000000d+01
      e(15) =   1071.5185372000d+00
      p(15) = 0.9276779724d-02
      e(16) =    252.6971215200d+00
      p(16) = 0.6954114943d-01
      e(17) =     80.1008291260d+00
      p(17) = 0.2715677256d+00
      e(18) =     28.9033931720d+00
      p(18) = 0.5340135557d+00
      e(19) =     10.7688998790d+00
      p(19) = 0.3450132345d+00
      e(20) =      5.6446212530d+00
      p(20) = 0.3412960016d+00
      e(21) =      2.1678291347d+00
      p(21) = 0.5639052197d+00
      e(22) =      0.8054089834d+00
      p(22) = 0.2367610974d+00
      e(23) =     56.0889391910d+00
      d(23) = 0.2958886914d-01
      e(24) =     15.7519089170d+00
      d(24) = 0.1587257140d+00
      e(25) =      5.3115812367d+00
      d(25) = 0.3797622916d+00
      e(26) =      1.7737904917d+00
      d(26) = 0.4689895917d+00
      e(27) =      0.5197558367d+00
      d(27) = 0.3090714908d+00
c
c additional p-GTO, steeper one of DP
c Ref.: Wachters, JCP, 52 (1970), 1033
c
      e(28) =   .162455d+00
      p(28) =  1.000000d+00
c
      return
      end
      subroutine ahsv4(e,s,p,d,n)
c
c -- ahlrich's SV (14s10p5d)/[5s4p2d] {63311/5311/41} contractions --
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-30
      go to (100,120,140,160,180,200),nn
c
c ga    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c optimized for atomic ground state ga(2p)
c scf energy is    -1922.9626417170 a.u.
c
  100 continue
      e( 1) =   87842.1262960000d+00
      s( 1) =  0.1427158880d-02
      e( 2) =   13189.4964020000d+00
      s( 2) =  0.1093889410d-01
      e( 3) =    3000.9482141000d+00
      s( 3) =  0.5430820013d-01
      e( 4) =     847.3842596600d+00
      s( 4) =  0.1895105598d+00
      e( 5) =     276.0098064200d+00
      s( 5) =  0.3861518592d+00
      e( 6) =      95.2166720710d+00
      s( 6) =  0.3005108272d+00
      e( 7) =     184.0070396800d+00
      s( 7) = -0.1112446196d+00
      e( 8) =      21.8270514330d+00
      s( 8) =  0.6483132292d+00
      e( 9) =       9.0026963761d+00
      s( 9) =  0.4435838359d+00
      e(10) =      16.0331769380d+00
      s(10) = -0.2301429956d+00
      e(11) =       2.6707724878d+00
      s(11) =  0.7294686152d+00
      e(12) =       1.1252834170d+00
      s(12) =  0.4621497695d+00
      e(13) =       0.1596760883d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0576430448d+00
      s(14) =  0.1000000000d+01
      e(15) =    1167.2665844000d+00
      p(15) =  0.9097402131d-02
      e(16) =     275.3806278900d+00
      p(16) =  0.6845627334d-01
      e(17) =      87.3750730700d+00
      p(17) =  0.2692229389d+00
      e(18) =      31.5972546700d+00
      p(18) =  0.5350789788d+00
      e(19) =      11.8241227980d+00
      p(19) =  0.3527906276d+00
      e(20) =       6.2881845133d+00
      p(20) =  0.3393824453d+00
      e(21) =       2.5199853240d+00
      p(21) =  0.5687194761d+00
      e(22) =       1.0169726306d+00
      p(22) =  0.2771769494d+00
      e(23) =       0.2474837360d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.0666852942d+00
      p(24) =  0.1000000000d+01
      e(25) =      65.3542376740d+00
      d(25) =  0.2737014638d-01
      e(26) =      18.5046567470d+00
      d(26) =  0.1509946398d+00
      e(27) =       6.3179620632d+00
      d(27) =  0.3748532855d+00
      e(28) =       2.1641389426d+00
      d(28) =  0.4751060685d+00
      e(29) =       0.6669309275d+00
      d(29) =  0.2984328740d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .207d+00
      d(30) =   1.000d+00
c
      return
c
c ge    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c optimized for atomic ground state ge(3p)
c scf energy is    -2075.0461658896 a.u. 
c 
  120 continue
      e( 1) =   93889.8366420000d+00
      s( 1) =  0.1423397606d-02
      e( 2) =   14097.4975280000d+00
      s( 2) =  0.1091079565d-01
      e( 3) =    3207.5477309000d+00
      s( 3) =  0.5418370594d-01
      e( 4) =     905.7672726900d+00
      s( 4) =  0.1892282035d+00
      e( 5) =     295.1101469300d+00
      s( 5) =  0.3861284700d+00
      e( 6) =     101.8471314100d+00
      s( 6) =  0.3016405074d+00
      e( 7) =     196.5671966200d+00
      s( 7) = -0.1111877094d+00
      e( 8) =      23.4052925220d+00
      s( 8) =  0.6461600737d+00
      e( 9) =       9.6839116702d+00
      s( 9) =  0.4418890457d+00
      e(10) =      17.2697365440d+00
      s(10) = -0.2302742138d+00
      e(11) =       2.8964622160d+00
      s(11) =  0.7301716940d+00
      e(12) =       1.2553621412d+00
      s(12) =  0.4619722226d+00
      e(13) =       0.2021308149d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0738679110d+00
      s(14) =  0.1000000000d+01
      e(15) =    1259.2085995000d+00
      p(15) =  0.9011546425d-02
      e(16) =     297.1562638200d+00
      p(16) =  0.6798684169d-01
      e(17) =      94.3533875220d+00
      p(17) =  0.2685385649d+00
      e(18) =      34.1763296770d+00
      p(18) =  0.5365964922d+00
      e(19) =      12.8161396150d+00
      p(19) =  0.3563351496d+00
      e(20) =       6.8471029784d+00
      p(20) =  0.3390069312d+00
      e(21) =       2.7717363939d+00
      p(21) =  0.5680936526d+00
      e(22) =       1.1458418175d+00
      p(22) =  0.2724653988d+00
      e(23) =       0.3067963154d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.0892836441d+00
      p(24) =  0.1000000000d+01
      e(25) =      74.7821681770d+00
      d(25) =  0.2575586021d-01
      e(26) =      21.3108497590d+00
      d(26) =  0.1453681613d+00
      e(27) =       7.3464792363d+00
      d(27) =  0.3713420986d+00
      e(28) =       2.5656271395d+00
      d(28) =  0.4800299844d+00
      e(29) =       0.8198177307d+00
      d(29) =  0.2897879074d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .246d+00
      d(30) =   1.000d+00
c
      return
c
c as    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c optimized for atomic ground state as(4s)
c scf energy is    -2233.9080800483 a.u.
c 
  140 continue
      e( 1) =  100146.5255400000d+00
      s( 1) =  0.1425834962d-02
      e( 2) =   15036.8617110000d+00
      s( 2) =  0.1093017696d-01
      e( 3) =    3421.2902833000d+00
      s( 3) =  0.5429417461d-01
      e( 4) =     966.1696571700d+00
      s( 4) =  0.1897607815d+00
      e( 5) =     314.8739402600d+00
      s( 5) =  0.3877519545d+00
      e( 6) =     108.7082379000d+00
      s( 6) =  0.3040281204d+00
      e( 7) =     209.5423895000d+00
      s( 7) = -0.1116209420d+00
      e( 8) =      25.0382211390d+00
      s( 8) =  0.6469760776d+00
      e( 9) =      10.3909643430d+00
      s( 9) =  0.4422360867d+00
      e(10) =      18.5550900930d+00
      s(10) = -0.2299419057d+00
      e(11) =       3.1281217449d+00
      s(11) =  0.7331910761d+00
      e(12) =       1.3884885073d+00
      s(12) =  0.4553365394d+00
      e(13) =       0.2471436214d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.0914294287d+00
      s(14) =  0.1000000000d+01
      e(15) =    1355.6443507000d+00
      p(15) =  0.8918250790d-02
      e(16) =     319.9992927000d+00
      p(16) =  0.6745475072d-01
      e(17) =     101.6773409200d+00
      p(17) =  0.2675977211d+00
      e(18) =      36.8863238450d+00
      p(18) =  0.5377684452d+00
      e(19) =      13.8611159090d+00
      p(19) =  0.3599257024d+00
      e(20) =       7.4260666912d+00
      p(20) =  0.3403684964d+00
      e(21) =       3.0316247187d+00
      p(21) =  0.5703014933d+00
      e(22) =       1.2783078340d+00
      p(22) =  0.2660617024d+00
      e(23) =       0.3756850336d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.1139480545d+00
      p(24) =  0.1000000000d+01
      e(25) =      84.4455145390d+00
      d(25) =  0.2451840272d-01
      e(26) =      24.1904161020d+00
      d(26) =  0.1410745468d+00
      e(27) =       8.4045015119d+00
      d(27) =  0.3687522892d+00
      e(28) =       2.9808970748d+00
      d(28) =  0.4840956136d+00
      e(29) =       0.9790924336d+00
      d(29) =  0.2825026878d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .293d+00
      d(30) =   1.000d+00
c
      return
c
c se    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c optimized for atomic ground state 3p
c scf energy is    -2399.5175767938 a.u.
c
  160 continue
      e( 1) =  106612.2002700000d+00
      s( 1) =  0.1427488911d-02
      e( 2) =   16007.6047010000d+00
      s( 2) =  0.1094352511d-01
      e( 3) =    3642.1699707000d+00
      s( 3) =  0.5437417160d-01
      e( 4) =    1028.5912993000d+00
      s( 4) =  0.1901809295d+00
      e( 5) =     335.3029888800d+00
      s( 5) =  0.3891302170d+00
      e( 6) =     115.8012915400d+00
      s( 6) =  0.3062020709d+00
      e( 7) =     222.9332502000d+00
      s( 7) = -0.1119880801d+00
      e( 8) =      26.7262579340d+00
      s( 8) =  0.6475212421d+00
      e( 9) =      11.1245019230d+00
      s( 9) =  0.4424197665d+00
      e(10) =      19.8885200610d+00
      s(10) = -0.2285722776d+00
      e(11) =       3.3668473803d+00
      s(11) =  0.7359135995d+00
      e(12) =       1.5249277839d+00
      s(12) =  0.4433019958d+00
      e(13) =       0.2963003791d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1100928897d+00
      s(14) =  0.1000000000d+01
      e(15) =    1455.9068120000d+00
      p(15) =  0.8820359704d-02
      e(16) =     343.7510183100d+00
      p(16) =  0.6687585197d-01
      e(17) =     109.2955496400d+00
      p(17) =  0.2664057808d+00
      e(18) =      39.7077110220d+00
      p(18) =  0.5383492842d+00
      e(19) =      14.9501852320d+00
      p(19) =  0.3630328199d+00
      e(20) =       8.0208962094d+00
      p(20) =  0.3415380703d+00
      e(21) =       3.2934649756d+00
      p(21) =  0.5725790658d+00
      e(22) =       1.4058602438d+00
      p(22) =  0.2554981322d+00
      e(23) =       0.4507612323d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.1335341333d+00
      p(24) =  0.1000000000d+01
      e(25) =      94.4940240440d+00
      d(25) =  0.2349010110d-01
      e(26) =      27.1881852600d+00
      d(26) =  0.1374773577d+00
      e(27) =       9.5091567352d+00
      d(27) =  0.3664992912d+00
      e(28) =       3.4170516853d+00
      d(28) =  0.4875098988d+00
      e(29) =       1.1479590083d+00
      d(29) =  0.2765794338d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .338d+00
      d(30) =   1.000d+00
c
      return
c
c br    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c optimized for atomic ground state br(2p)
c scf energy is    -2572.0704481746 a.u. 
c 
  180 continue
      e( 1) =  113286.3877600000d+00
      s( 1) =  0.1428303778d-02
      e( 2) =   17009.6263030000d+00
      s( 2) =  0.1095041750d-01
      e( 3) =    3870.1842567000d+00
      s( 3) =  0.5442100660d-01
      e( 4) =    1093.0357227000d+00
      s( 4) =  0.1904790770d+00
      e( 5) =     356.3972179700d+00
      s( 5) =  0.3902464274d+00
      e( 6) =     123.1253964300d+00
      s( 6) =  0.3081443251d+00
      e( 7) =     236.7408400700d+00
      s( 7) = -0.1122806567d+00
      e( 8) =      28.4686610700d+00
      s( 8) =  0.6477596231d+00
      e( 9) =      11.8834437220d+00
      s( 9) =  0.4423557599d+00
      e(10) =      21.2696333120d+00
      s(10) = -0.2264257632d+00
      e(11) =       3.6129226841d+00
      s(11) =  0.7382371201d+00
      e(12) =       1.6626648969d+00
      s(12) =  0.4268386869d+00
      e(13) =       0.3482379323d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1301903139d+00
      s(14) =  0.1000000000d+01
      e(15) =    1560.2801881000d+00
      p(15) =  0.8716666907d-02
      e(16) =     368.4785920500d+00
      p(16) =  0.6624363742d-01
      e(17) =     117.2297884900d+00
      p(17) =  0.2649561039d+00
      e(18) =      42.6489092480d+00
      p(18) =  0.5383916059d+00
      e(19) =      16.0872250960d+00
      p(19) =  0.3657938789d+00
      e(20) =       8.6352810058d+00
      p(20) =  0.3424878737d+00
      e(21) =       3.5613665502d+00
      p(21) =  0.5750067821d+00
      e(22) =       1.5292626609d+00
      p(22) =  0.2433039417d+00
      e(23) =       0.5306429485d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.1570275897d+00
      p(24) =  0.1000000000d+01
      e(25) =     104.8551864200d+00
      d(25) =  0.2265014758d-01
      e(26) =      30.2811436880d+00
      d(26) =  0.1345548323d+00
      e(27) =      10.6513942670d+00
      d(27) =  0.3647445454d+00
      e(28) =       3.8699456233d+00
      d(28) =  0.4904458706d+00
      e(29) =       1.3240876762d+00
      d(29) =  0.2713728904d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .389d+00
      d(30) =   1.000d+00
c
      return
c
c kr    (14s10p5d) / [5s4p2d]     {63311/5311/41}
c scf energy is    -2751.6619907725 a.u. 
c
  200 continue
      e( 1) =  120165.6487500000d+00
      s( 1) =  0.1429423450d-02
      e( 2) =   18042.5001690000d+00
      s( 2) =  0.1095959702d-01
      e( 3) =    4105.1800388000d+00
      s( 3) =  0.5447964747d-01
      e( 4) =    1159.4447248000d+00
      s( 4) =  0.1908123912d+00
      e( 5) =     378.1381034600d+00
      s( 5) =  0.3914044595d+00
      e( 6) =     130.6761004500d+00
      s( 6) =  0.3100734669d+00
      e( 7) =     250.9637335300d+00
      s( 7) = -0.1125860401d+00
      e( 8) =      30.2659303800d+00
      s( 8) =  0.6481632252d+00
      e( 9) =      12.6681107570d+00
      s( 9) =  0.4424023899d+00
      e(10) =      22.6979401010d+00
      s(10) = -0.2239196491d+00
      e(11) =       3.8673017398d+00
      s(11) =  0.7406095303d+00
      e(12) =       1.8013706671d+00
      s(12) =  0.4080365845d+00
      e(13) =       0.4030863914d+00
      s(13) =  0.1000000000d+01
      e(14) =       0.1517052775d+00
      s(14) =  0.1000000000d+01
      e(15) =    1668.5736623000d+00
      p(15) =  0.8621500980d-02
      e(16) =     394.1386279800d+00
      p(16) =  0.6566550688d-01
      e(17) =     125.4664485400d+00
      p(17) =  0.2636657670d+00
      e(18) =      45.7049927290d+00
      p(18) =  0.5386765857d+00
      e(19) =      17.2701332290d+00
      p(19) =  0.3686516750d+00
      e(20) =       9.2699487770d+00
      p(20) =  0.3434618902d+00
      e(21) =       3.8364617719d+00
      p(21) =  0.5778477620d+00
      e(22) =       1.6467163312d+00
      p(22) =  0.2307561520d+00
      e(23) =       0.6147968649d+00
      p(23) =  0.1000000000d+01
      e(24) =       0.1831946741d+00
      p(24) =  0.1000000000d+01
      e(25) =     115.5586860300d+00
      d(25) =  0.2194548533d-01
      e(26) =      33.4779858640d+00
      d(26) =  0.1321121288d+00
      e(27) =      11.8344502070d+00
      d(27) =  0.3633330696d+00
      e(28) =       4.3408254894d+00
      d(28) =  0.4930048445d+00
      e(29) =       1.5079273157d+00
      d(29) =  0.2667685891d+00
c
c additional d-GTO
c Ref.: PSD 16 (Huzinaga), p. 23 
c
      e(30) =   .443d+00
      d(30) =   1.000d+00
c
      return
      end
      subroutine ahsvsh(nucz,ipass,opol,itype,igauss,ng,
     *                  scale,scs,scp,scd,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(11,9), itypes(11,9), ioff(11,9)
      dimension igaup(2,9), itypep(2,9), ioffp(2,9)
      dimension kind(9)
      data kind/2,3,3,5,5,7,7,9,11/
      data igau /
     + 3,1,0,0,0,0,0,0,0,0,0,
     + 5,1,1,0,0,0,0,0,0,0,0,
     + 5,1,1,0,0,0,0,0,0,0,0,
     + 5,1,1,3,1,0,0,0,0,0,0,
     + 5,3,1,1,5,0,0,0,0,0,0,
     + 5,3,1,1,5,1,1,0,0,0,0,
     + 6,3,3,1,1,5,3,0,0,0,0,
     + 6,3,3,1,1,5,3,4,1,0,0,
     + 6,3,3,1,1,5,3,1,1,4,1/
      data igaup /
     + 1,0,
     + 1,0,
     + 3,1,
     + 1,0,
     + 1,0,
     + 1,0,
     + 1,0,
     + 1,0,
     + 1,0/
      data itypes /
     + 1,1,0,0,0,0,0,0,0,0,0,
     + 1,1,1,0,0,0,0,0,0,0,0,
     + 1,1,1,0,0,0,0,0,0,0,0,
     + 1,1,1,2,2,0,0,0,0,0,0,
     + 1,1,1,1,2,0,0,0,0,0,0,
     + 1,1,1,1,2,2,2,0,0,0,0,
     + 1,1,1,1,1,2,2,0,0,0,0,
     + 1,1,1,1,1,2,2,3,3,0,0,
     + 1,1,1,1,1,2,2,2,2,3,3 /
      data itypep/
     + 2,0,
     + 2,0,
     + 2,2,
     + 3,0,
     + 2,0,
     + 3,0,
     + 2,0,
     + 2,0,
     + 3,0 /
c
c    pointers to the sv basis
c
      data ioff/0, 3, 9*0,
     +          0, 5, 6, 8*0,
     +          0, 5, 6, 8*0,
     +          0, 5, 6, 7,10,6*0,
     +          0, 5, 8, 9,10,6*0,
     +          0, 5, 8, 9,10,15,16,4*0,
     +          0, 6, 9,12,13,14,19,4*0,
     +          0, 6, 9,12,13,14,19,22,26,2*0,
     +          0, 6, 9,12,13,14,19,22,23,24,28/ 
c
c    pointers to the sv+p basis
c
      data ioffp/4, 0,
     +           7, 0,
     +           7, 10,
     +          11, 0,
     +          15, 0,
     +          17, 0,
     +          22, 0,
     +          27, 0,
     +          29, 0 /
c
      data done/1.0d+00/
c     data dzero/0.0d+00/
c
c     set values for the current ahlrichs sv  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.gt. 4) ind=4
      if(nucz.gt.10) ind=5
      if(nucz.gt.12) ind=6
      if(nucz.gt.18) ind=7
      if(nucz.gt.20) ind=8
      if(nucz.gt.30) ind=9
      if(nucz.gt.36) then
         call caserr2('ahlrichs basis sets only extend to krypton')
      end if
c
c
      mxpass=kind(ind)
c                 all atoms have one polarisation function excpt be
      npol = 1
      if(nucz.eq.4 ) npol = 2
      if(opol) mxpass=mxpass+npol
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c              ******
c
      kkkk = ipass-kind(ind)
      if (kkkk.le.0) then
       igauss = igau(ipass,ind)
       ng = ioff(ipass,ind) 
       ityp = itypes(ipass,ind)
      else
       igauss = igaup(kkkk,ind)
       ng = ioffp(kkkk,ind) 
       ityp = itypep(kkkk,ind)
      endif
c
c          do not scale midi, or mini inner shells
c
      scale=done
      if(ityp.eq.1) scale=scs
      if(ityp.eq.2) scale=scp
      if(ityp.eq.3) scale=scd
100   itype = ityp + 15
      return
      end
      subroutine ahtz0(e,s,p,n)
c
c     ----- ahlrich's contraction of (5s) -----
c     ----- hydrogen + helium (5s)/(3s) -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*)
      go to (100,120),n
c
c     ----- h TZ -----
c
c h     (5s) / [3s]     {311}
c scf energy is        -.4998098322 a.u.
c
  100 continue
      e( 1) =      34.0613410000d+00
      s( 1) =  0.6025197800d-02
      e( 2) =       5.1235746000d+00
      s( 2) =  0.4502109400d-01
      e( 3) =       1.1646626000d+00
      s( 3) =  0.2018972600d+00
      e( 4) =       0.3272304100d+00
      s( 4) =  0.1000000000d+01
      e( 5) =       0.1030724100d+00
      s( 5) =  0.1000000000d+01
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e( 6) =  .80d+00
      p( 6) = 1.00d+00
c
      return
c
c     ----- he TZ -----
c
c he    (5s) / [3s]     {311}
c scf energy is       -2.8598954257 a.u.
c
  120 continue
      e( 1) =      98.0783216160d+00
      s( 1) =  0.7580306497d-02
      e( 2) =      14.7644042470d+00
      s( 2) =  0.5484862094d-01
      e( 3) =       3.3185831473d+00
      s( 3) =  0.2207438219d+00
      e( 4) =       0.8741386955d+00
      s( 4) =  0.1000000000d+01
      e( 5) =       0.2445989721d+00
      s( 5) =  0.1000000000d+01
c
c additional p-GTO
c Ref.: R. Ahlrichs 13.08.92
c
      e( 6) =   1.000d+00
      p( 6) =   1.000d+00
      return
c
      end
      subroutine ahtz1(e,s,p,d,n)
c
c  ----- alhrich's TZ (10s6p)/[6s3p]{511111/411} -----
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-2
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- li TZ -----
c
c li    (10s) / [6s]     {511111}
c scf energy is       -7.4326172485 a.u.
c
  100 continue
      e( 1) =    3341.4028812000d+00
      s( 1) =  0.6882573495d-03
      e( 2) =     500.9515931100d+00
      s( 2) =  0.5329642996d-02
      e( 3) =     114.0110310600d+00
      s( 3) =  0.2756437364d-01
      e( 4) =      32.2654730600d+00
      s( 4) =  0.1098824256d+00
      e( 5) =      10.4720585230d+00
      s( 5) =  0.3394304925d+00
      e( 6) =       3.6978284154d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       1.3712085293d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.5190845988d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.0745328871d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.0282336001d+00
      s(10) =  0.1000000000d+01
c
c two additional p-GTOs
c Ref.: C. Ochsenfeld - optimised for LinFm
c
      e(11) =  0.40d+00
      p(11) =  1.00d+00
      e(12) =  0.06d+00
      p(12) =  1.00d+00
c
      return
c
c     ----- be TZ -----
c
c Be    (10s) / [6s]     {511111}
c scf energy is      -14.5726557698 a.u.
c
  120 continue
      e( 1) =    6526.1477796000d+00
      s( 1) =  0.4373173971d-03
      e( 2) =     978.3549641700d+00
      s( 2) =  0.3388207195d-02
      e( 3) =     222.6586575600d+00
      s( 3) =  0.1756308508d-01
      e( 4) =      63.0130282810d+00
      s( 4) =  0.7058048133d-01
      e( 5) =      20.4687201820d+00
      s( 5) =  0.2226383623d+00
      e( 6) =       7.2692487387d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       2.7393701115d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       1.0615791118d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.1815844019d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.0588268200d+00
      s(10) =  0.1000000000d+01
c
c p set taken from TZ.3P basis
c be    (10s6p) / [6s3p]     {511111/411}
c optimized for 3P state, i.e. 1s(2)2s(1)2p(1)
c
      e(11) =  14.018558828d+00
      p(11) =      -.67974011741d-02
      e(12) =  3.1647356439d+00
      p(12) =      -.42491136104d-01
      e(13) =  .90157670087d+00
      p(13) =      -.17210557921d+00
      e(14) =  .30279153007d+00
      p(14) =      -.51791035465d+00
      e(15) =  .11279491554d+00
      p(15) =      1.0000000000d+00
      e(16) =  .42973443575d-01
      p(16) =      1.0000000000d+00
c
      return
c
c     ----- b TZ  -----
c
c B     (10s6p) / [6s3p]     {511111/411}
c scf energy is      -24.5282379803 a.u.
c
  140 continue
      e( 1) =    5658.4992000000d+00
      s( 1) =  0.5520832700d-03
      e( 2) =     848.6970300000d+00
      s( 2) =  0.4265142000d-02
      e( 3) =     193.1365600000d+00
      s( 3) =  0.2185316500d-01
      e( 4) =      54.5863170000d+00
      s( 4) =  0.8457067800d-01
      e( 5) =      17.6001260000d+00
      s( 5) =  0.2430476400d+00
      e( 6) =       6.1269618000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       2.2245270000d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.5565035200d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.2387622900d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.0857939550d+00
      s(10) =  0.1000000000d+01
      e(11) =      22.4558020000d+00
      p(11) = -0.5026062100d-02
      e(12) =       5.1030675000d+00
      p(12) = -0.3282697600d-01
      e(13) =       1.4984096000d+00
      p(13) = -0.1314094100d+00
      e(14) =       0.5098907300d+00
      p(14) = -0.3312526400d+00
      e(15) =       0.1819135000d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.0648271830d+00
      p(16) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(17) =   .50d+00
      d(17) =    1.00d+00
c
      return
c
c     ----- c TZ  -----
c
c C     (10s6p) / [6s3p]     {511111/411}
c scf energy is      -37.6872906161 a.u.
c
  160 continue
      e( 1) =    8506.0384000000d+00
      s( 1) =  0.5337366400d-03
      e( 2) =    1275.7329000000d+00
      s( 2) =  0.4125023200d-02
      e( 3) =     290.3118700000d+00
      s( 3) =  0.2117133700d-01
      e( 4) =      82.0562000000d+00
      s( 4) =  0.8241786000d-01
      e( 5) =      26.4796410000d+00
      s( 5) =  0.2401285800d+00
      e( 6) =       9.2414585000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       3.3643530000d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       0.8717416400d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.3635235200d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.1287313500d+00
      s(10) =  0.1000000000d+01
      e(11) =      34.7094960000d+00
      p(11) =  0.5330097400d-02
      e(12) =       7.9590883000d+00
      p(12) =  0.3586581400d-01
      e(13) =       2.3786972000d+00
      p(13) =  0.1420029900d+00
      e(14) =       0.8154006500d+00
      p(14) =  0.3420310500d+00
      e(15) =       0.2895378500d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.1008475400d+00
      p(16) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(17) =   .80d+00
      d(17) =   1.00d+00
c
      return
c
c     ----- n TZ -----
c
c N     (10s6p) / [6s3p]     {511111/411}
c scf energy is      -54.3988903986 a.u.
c
  180 continue
      e( 1) =   11913.4167560000d+00
      s( 1) = -0.5229701737d-03
      e( 2) =    1786.7213834000d+00
      s( 2) = -0.4042840404d-02
      e( 3) =     406.5901283400d+00
      s( 3) = -0.2077271513d-01
      e( 4) =     114.9252506500d+00
      s( 4) = -0.8118313785d-01
      e( 5) =      37.1058834220d+00
      s( 5) = -0.2387149740d+00
      e( 6) =      12.9716761980d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       4.7302291164d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       1.2525184258d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.5126007123d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.1793971400d+00
      s(10) =  0.1000000000d+01
      e(11) =      49.2187583920d+00
      p(11) =  0.5552695348d-02
      e(12) =      11.3489353040d+00
      p(12) =  0.3805461796d-01
      e(13) =       3.4285088246d+00
      p(13) =  0.1494141224d+00
      e(14) =       1.1799512559d+00
      p(14) =  0.3489818737d+00
      e(15) =       0.4172612258d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.1429513128d+00
      p(16) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(17) =   1.00d+00
      d(17) =   1.00d+00
c
      return
c
c     ----- o  TZ ------
c
c O     (10s6p) / [6s3p]     {511111/311}
c scf energy is      -74.8062891051 a.u.
c
  200 continue
      e( 1) =   15902.6474590000d+00
      s( 1) =  0.5149980370d-03
      e( 2) =    2384.9537829000d+00
      s( 2) =  0.3981976443d-02
      e( 3) =     542.7195718200d+00
      s( 3) =  0.2047697192d-01
      e( 4) =     153.4040787400d+00
      s( 4) =  0.8026236792d-01
      e( 5) =      49.5457161400d+00
      s( 5) =  0.2376683995d+00
      e( 6) =      17.3396498970d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       6.3303355272d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       1.6995882201d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.6895449127d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.2393602818d+00
      s(10) =  0.1000000000d+01
      e(11) =      63.2705240110d+00
      p(11) =  0.6070920596d-02
      e(12) =      14.6233122950d+00
      p(12) =  0.4194768872d-01
      e(13) =       4.4489518003d+00
      p(13) =  0.1615688399d+00
      e(14) =       1.5281513180d+00
      p(14) =  0.3568277929d+00
      e(15) =       0.5299731587d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.1750944600d+00
      p(16) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c
      e(17) =   1.20d+00
      d(17) =   1.00d+00
c
      return
c
c     ----- f  TZ -----
c
c f     (10s6p) / [6s3p]     {511111/411}
c scf energy is      -99.4048383121 a.u.
c
  220 continue
      e( 1) =   20450.4890000000d+00
      s( 1) =  0.5110349500d-03
      e( 2) =    3066.9547000000d+00
      s( 2) =  0.3951882000d-02
      e( 3) =     697.9100300000d+00
      s( 3) =  0.2033455300d-01
      e( 4) =     197.2702000000d+00
      s( 4) =  0.7987648000d-01
      e( 5) =      63.7283430000d+00
      s( 5) =  0.2377560100d+00
      e( 6) =      22.3218090000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =       8.1557609000d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       2.2114295000d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       0.8903856700d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.3069660400d+00
      s(10) =  0.1000000000d+01
      e(11) =      80.2180200000d+00
      p(11) =  0.6374446400d-02
      e(12) =      18.5872810000d+00
      p(12) =  0.4436019100d-01
      e(13) =       5.6844581000d+00
      p(13) =  0.1688003800d+00
      e(14) =       1.9512781000d+00
      p(14) =  0.3616297900d+00
      e(15) =       0.6702411400d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.2168225200d+00
      p(16) =  0.1000000000d+01
c
c     ***** f TZP
c additional d-GTO
c Ref.: R. Ahlrichs 12.08.92
c depending on atomic charge, values down to 0.4 result in molecular
c optimization, e.g. for Li2F.
c
      d(17)=     1.40d0
      e(17)=     1.00d0
c
      return
c
c     ----- ne TZ -----
c
c ne    (10s6p) / [6s3p]     {511111/411}
c scf energy is     -128.5408109343 a.u.
c
  240 continue
      e( 1) =   25558.3983200000d+00
      s( 1) = -0.5057269197d-03
      e( 2) =    3832.9404804000d+00
      s( 2) = -0.3911248585d-02
      e( 3) =     872.2107377700d+00
      s( 3) = -0.2013482353d-01
      e( 4) =     246.5375541400d+00
      s( 4) = -0.7922342728d-01
      e( 5) =      79.6583948640d+00
      s( 5) = -0.2367579848d+00
      e( 6) =      27.9197996460d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      10.2071316640d+00
      s( 7) =  0.1000000000d+01
      e( 8) =       2.7883186671d+00
      s( 8) =  0.1000000000d+01
      e( 9) =       1.1153994748d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.3823401054d+00
      s(10) =  0.1000000000d+01
      e(11) =      99.7609664310d+00
      p(11) =  0.6564866162d-02
      e(12) =      23.1676889530d+00
      p(12) =  0.4595815153d-01
      e(13) =       7.1133898669d+00
      p(13) =  0.1734998103d+00
      e(14) =       2.4418574365d+00
      p(14) =  0.3648607515d+00
      e(15) =       0.8344206317d+00
      p(15) =  0.1000000000d+01
      e(16) =       0.2663988109d+00
      p(16) =  0.1000000000d+01
c 
c additional d-GTO 
c Ref.: PSD 16 (Huzinaga), p. 23
c
      e(17) =   1.888d+00
      d(17) = 1.000d+00
c
      return
      end
      subroutine ahtz2(e,s,p,d,n)
c
c  --- ahlrich's TZ (12s9p)/[7s5p] {5121111/51111} contraction ---
c
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
      nn = n-10
      go to (100,120,140,160,180,200,220,240),nn
c
c     ----- na TZ -----
c
c na    (12s6p) / [7s2p]     {5121111/51}
c scf energy is     -161.8519976833 a.u.
c
  100 continue
      e( 1) =   48876.1981470000d+00
      s( 1) =  0.3049570984d-03
      e( 2) =    7327.2863708000d+00
      s( 2) =  0.2363283260d-02
      e( 3) =    1667.4194652000d+00
      s( 3) =  0.1226845080d-01
      e( 4) =     471.5810119600d+00
      s( 4) =  0.4962517661d-01
      e( 5) =     152.9366981100d+00
      s( 5) =  0.1594959178d+00
      e( 6) =      54.3558554660d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      20.7652228110d+00
      s( 7) =  0.4007204643d+00
      e( 8) =       8.3599495733d+00
      s( 8) =  0.1833737271d+00
      e( 9) =       2.1480945935d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.6881941424d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.0604771041d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.0271340943d+00
      s(12) =  0.1000000000d+01
      e(13) =     137.9625476100d+00
      p(13) =  0.5926462784d-02
      e(14) =      32.2091053610d+00
      p(14) =  0.4254399271d-01
      e(15) =       9.9746379096d+00
      p(15) =  0.1664453483d+00
      e(16) =       3.4780930022d+00
      p(16) =  0.3681491453d+00
      e(17) =       1.2284371548d+00
      p(17) =  0.4578103525d+00
      e(18) =       0.4173193641d+00
      p(18) =  0.1000000000d+01
c
c additional p-function
c Ref.: R. Ahlrichs 12.08.92, as steepr s-GTO of TZ.
c which is optimal for NanClm
c
      e(19) = .052d+00
      p(19) =  1.00d+00
c
      return
c
c     ----- mg TZ -----
c
c mg    (12s6p) / [7s2p]     {5121111/51}
c scf energy is      -199.6074427183 a.u.
c 
  120 continue
      e( 1) =   58088.3497780000d+00
      s( 1) =  0.2921796983d-03
      e( 2) =    8708.1180998000d+00
      s( 2) =  0.2264423666d-02
      e( 3) =    1981.6224182000d+00
      s( 3) =  0.1175603987d-01
      e( 4) =     560.4398599700d+00
      s( 4) =  0.4756994318d-01
      e( 5) =     181.7297610800d+00
      s( 5) =  0.1530804790d+00
      e( 6) =      64.6033357290d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      24.7421446700d+00
      s( 7) =  0.3988356060d+00
      e( 8) =      10.0246647160d+00
      s( 8) =  0.1818881521d+00
      e( 9) =       2.6947441296d+00
      s( 9) =  0.1000000000d+01
      e(10) =       0.9114525192d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1068967301d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.0400916171d+00
      s(12) =  0.1000000000d+01
      e(13) =     179.8337272400d+00
      p(13) =  0.5381976916d-02
      e(14) =      42.1130966640d+00
      p(14) =  0.3932575177d-01
      e(15) =      13.1173537750d+00
      p(15) =  0.1574904406d+00
      e(16) =       4.6229435410d+00
      p(16) =  0.3594056082d+00
      e(17) =       1.6685785801d+00
      p(17) =  0.4551622999d+00
      e(18) =       0.5854779723d+00
      p(18) =  0.2197049159d+00
c
c two additional p-functions
c taken from TZ.3p basis set
c
      e(19) = .18914796195d+00
      p(19) =  1.00d+00
      e(20) = .53768755187d-01
      p(20) =  1.00d+00
      return
c
c     ----- al TZ -----
c
c al    (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -241.8718992955 a.u.
c 
  140 continue
      e( 1) =   68090.4116770000d+00
      s( 1) =  0.2841290542d-03
      e( 2) =   10207.5408220000d+00
      s( 2) =  0.2202069974d-02
      e( 3) =    2322.8045854000d+00
      s( 3) =  0.1143344385d-01
      e( 4) =     656.8865383800d+00
      s( 4) =  0.4628458562d-01
      e( 5) =     212.9669410400d+00
      s( 5) =  0.1491056579d+00
      e( 6) =      75.7247165850d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      29.0763200240d+00
      s( 7) =  0.4166351811d+00
      e( 8) =      11.8596239340d+00
      s( 8) =  0.1897837136d+00
      e( 9) =       3.3162153181d+00
      s( 9) =  0.1000000000d+01
      e(10) =       1.1754788791d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.1752052855d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.0647587465d+00
      s(12) =  0.1000000000d+01
      e(13) =     442.7918258700d+00
      p(13) =  0.1639111843d-02
      e(14) =     104.7598883700d+00
      p(14) =  0.1314351546d-01
      e(15) =      33.3759794470d+00
      p(15) =  0.6149770806d-01
      e(16) =      12.3053879570d+00
      p(16) =  0.1882590885d+00
      e(17) =       4.8642475706d+00
      p(17) =  0.3606261120d+00
      e(18) =       1.9604027100d+00
      p(18) =  0.1000000000d+01
      e(19) =       0.7843673382d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.1901416508d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.0558743730d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(22) = .30d+00
      d(22) =  1.00d+00
c
      return
c
c     ----- si TZ -----
c
c si    (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -288.8485540288 a.u.
c
  160 continue
      e( 1) =   79079.4340000000d+00
      s( 1) =  0.2643138600d-03
      e( 2) =   11855.0100000000d+00
      s( 2) =  0.2048514300d-02
      e( 3) =    2697.7051000000d+00
      s( 3) =  0.1063724100d-01
      e( 4) =     762.8722700000d+00
      s( 4) =  0.4308247700d-01
      e( 5) =     247.2845500000d+00
      s( 5) =  0.1389827900d+00
      e( 6) =      87.9312400000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      33.8232840000d+00
      s( 7) =  0.4407154300d+00
      e( 8) =      13.8681080000d+00
      s( 8) =  0.2009116500d+00
      e( 9) =       3.9920017000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       1.4659925000d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.2527108600d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.0924916730d+00
      s(12) =  0.1000000000d+01
      e(13) =     483.0235200000d+00
      p(13) =  0.1916154700d-02
      e(14) =     114.2508100000d+00
      p(14) =  0.1530976500d-01
      e(15) =      36.3877860000d+00
      p(15) =  0.7109435800d-01
      e(16) =      13.4117040000d+00
      p(16) =  0.2124324400d+00
      e(17) =       5.2884033000d+00
      p(17) =  0.3897630200d+00
      e(18) =       2.1374219000d+00
      p(18) =  0.1000000000d+01
      e(19) =       0.8646846300d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.2548985500d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.0793970310d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(22) = .35d+00
      d(22) =  1.00d+00
c
      return
c
c     ----- p TZ -----
c
c p     (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -340.7119214717 a.u.
c
  180 continue
      e( 1) =   91049.0495250000d+00
      s( 1) =  0.2538017465d-03
      e( 2) =   13649.2488640000d+00
      s( 2) =  0.1967123901d-02
      e( 3) =    3105.9432074000d+00
      s( 3) =  0.1021606304d-01
      e( 4) =     878.2598263000d+00
      s( 4) =  0.4139916202d-01
      e( 5) =     284.6339290700d+00
      s( 5) =  0.1337541498d+00
      e( 6) =     101.2061000500d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      38.9778414240d+00
      s( 7) =  0.4541611301d+00
      e( 8) =      16.0484940730d+00
      s( 8) =  0.2074219800d+00
      e( 9) =       4.7225339088d+00
      s( 9) =  0.1000000000d+01
      e(10) =       1.7825221852d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.3396296788d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.1233385985d+00
      s(12) =  0.1000000000d+01
      e(13) =     528.4438932800d+00
      p(13) = -0.2281519015d-02
      e(14) =     124.9714995100d+00
      p(14) = -0.1818364204d-01
      e(15) =      39.7957092480d+00
      p(15) = -0.8394070053d-01
      e(16) =      14.6575449730d+00
      p(16) = -0.2460471387d+00
      e(17) =       5.7578892754d+00
      p(17) = -0.4365566753d+00
      e(18) =       2.3123200955d+00
      p(18) =  0.1000000000d+01
      e(19) =       0.8899378100d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.3046714991d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1007815396d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(22) = .45d+00
      d(22) = 1.00d+00
      return
c
c     ----- s TZ  -----
c
c s     (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -397.4966169684 a.u.
c
  200 continue
      e( 1) =  103953.9500000000d+00
      s( 1) =  0.2474271700d-03
      e( 2) =   15583.7860000000d+00
      s( 2) =  0.1917768300d-02
      e( 3) =    3546.1293000000d+00
      s( 3) =  0.9960999500d-02
      e( 4) =    1002.6808000000d+00
      s( 4) =  0.4038642900d-01
      e( 5) =     324.9028800000d+00
      s( 5) =  0.1306752300d+00
      e( 6) =     115.5122500000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      44.5282100000d+00
      s( 7) =  0.5040355100d+00
      e( 8) =      18.3978920000d+00
      s( 8) =  0.2306805100d+00
      e( 9) =       5.5100683000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       2.1259867000d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.4369189300d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.1573085100d+00
      s(12) =  0.1000000000d+01
      e(13) =     606.6936700000d+00
      p(13) =  0.2322964500d-02
      e(14) =     143.5069600000d+00
      p(14) =  0.1857053200d-01
      e(15) =      45.7461600000d+00
      p(15) =  0.8603274000d-01
      e(16) =      16.8729120000d+00
      p(16) =  0.2524843100d+00
      e(17) =       6.6399196000d+00
      p(17) =  0.4463274200d+00
      e(18) =       2.6727135000d+00
      p(18) =  0.1000000000d+01
      e(19) =       1.0000089000d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.3543894200d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1167130500d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, old library.
c
      e(22) = .55d+00
      d(22) = 1.00d+00
      return
c
c     ----- cl TZ -----
c
c cl    (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -459.4730301720 a.u.
c
  220 continue
      e( 1) =  117805.7900000000d+00
      s( 1) =  0.2418075400d-03
      e( 2) =   17660.2650000000d+00
      s( 2) =  0.1874268000d-02
      e( 3) =    4018.5973000000d+00
      s( 3) =  0.9736293500d-02
      e( 4) =    1136.2231000000d+00
      s( 4) =  0.3949474600d-01
      e( 5) =     368.1206200000d+00
      s( 5) =  0.1279716600d+00
      e( 6) =     130.8615400000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      50.4790110000d+00
      s( 7) =  0.4287409400d+00
      e( 8) =      20.9168050000d+00
      s( 8) =  0.1966849300d+00
      e( 9) =       6.3531388000d+00
      s( 9) =  0.1000000000d+01
      e(10) =       2.4948014000d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.5433594800d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.1943436600d+00
      s(12) =  0.1000000000d+01
      e(13) =     681.0687900000d+00
      p(13) =  0.2365893000d-02
      e(14) =     161.1135900000d+00
      p(14) =  0.1894071700d-01
      e(15) =      51.3866370000d+00
      p(15) =  0.8784513400d-01
      e(16) =      18.9585120000d+00
      p(16) =  0.2570743500d+00
      e(17) =       3.0035158000d+00
      p(17) =  0.3715236100d+00
      e(18) =       7.4565293000d+00
      p(18) =  0.1000000000d+01
      e(19) =       1.0609361000d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.3945202300d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1332327900d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c Ref.: R. Ahlrichs 13.08.92, 
c depending on atomic charge, values down to 0.2 result in molecular
c optimization, e.g. for NanClm.
c
      e(22) =   .65d+00
      d(22) =  1.00d+00
c
      return
c
c     ----- ar TZ -----
c
c ar    (12s9p) / [7s5p]     {5121111/51111}
c scf energy is     -526.8056097367 a.u.
c
  240 continue
      e( 1) =  132594.7425500000d+00
      s( 1) =  0.2548411516d-03
      e( 2) =   19877.2100510000d+00
      s( 2) =  0.1975334879d-02
      e( 3) =    4523.0525989000d+00
      s( 3) =  0.1026245695d-01
      e( 4) =    1278.8031045000d+00
      s( 4) =  0.4164867475d-01
      e( 5) =     414.2561894000d+00
      s( 5) =  0.1351337705d+00
      e( 6) =     147.2419411000d+00
      s( 6) =  0.1000000000d+01
      e( 7) =      56.8261553810d+00
      s( 7) =  0.4514163288d+00
      e( 8) =      23.6040859760d+00
      s( 8) =  0.2075883242d+00
      e( 9) =       7.2524952042d+00
      s( 9) =  0.1000000000d+01
      e(10) =       2.8892979871d+00
      s(10) =  0.1000000000d+01
      e(11) =       0.6588798270d+00
      s(11) =  0.1000000000d+01
      e(12) =       0.2343956108d+00
      s(12) =  0.1000000000d+01
      e(13) =     828.3396980400d+00
      p(13) =  0.2177303934d-02
      e(14) =     196.0651075800d+00
      p(14) =  0.1759744324d-01
      e(15) =      62.6812275850d+00
      p(15) =  0.8290592139d-01
      e(16) =      23.2250778440d+00
      p(16) =  0.2493410292d+00
      e(17) =       9.2167815494d+00
      p(17) =  0.4532728064d+00
      e(18) =       3.7595532955d+00
      p(18) =  0.1000000000d+01
      e(19) =       1.4316785909d+00
      p(19) =  0.1000000000d+01
      e(20) =       0.5129826610d+00
      p(20) =  0.1000000000d+01
      e(21) =       0.1682148695d+00
      p(21) =  0.1000000000d+01
c
c additional d-GTO
c taken from tzp basis set in hondo basis set library
c
      e(22) =  .696d+00
      d(22) =  1.00d+00
c
      return
      end
      subroutine ahtzsh(nucz,ipass,opol,itype,igauss,ng,
     *                  scale,scs,scp,scd,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(12,7), itypes(12,7), ioff(12,7)
      dimension igaup(3,7), itypep(3,7), ioffp(3,7)
      dimension kind(7)
      data kind/3,6,6,9,9,9,12/
      data igau /
     + 3,1,1,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,1,1,0,0,0,0,0,0,
     + 5,1,1,1,1,1,0,0,0,0,0,0,
     + 5,1,1,1,1,1,4,1,1,0,0,0,
     + 5,1,2,1,1,1,1,5,1,0,0,0,
     + 5,1,2,1,1,1,1,5,1,0,0,0,
     + 5,1,2,1,1,1,1,5,1,1,1,1 /
      data igaup /
     + 1,0,0,
     + 1,1,0,
     + 4,1,1,
     + 1,0,0,
     + 1,0,0,
     + 1,1,0,
     + 1,0,0 /
      data itypes /
     + 1,1,1,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,0,0,0,0,0,0,
     + 1,1,1,1,1,1,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,0,0,0,
     + 1,1,1,1,1,1,1,2,2,0,0,0,
     + 1,1,1,1,1,1,1,2,2,0,0,0,
     + 1,1,1,1,1,1,1,2,2,2,2,2 /
      data itypep/
     + 2,0,0,
     + 2,2,0,
     + 2,2,2,
     + 3,0,0,
     + 2,0,0,
     + 2,2,0,
     + 3,0,0/
c
c    pointers to the tz basis
c
      data ioff/0, 3, 4, 9*0,
     +          0, 5, 6, 7, 8, 9, 6*0,
     +          0, 5, 6, 7, 8, 9, 6*0,
     +          0, 5, 6, 7, 8, 9,10,14,15,3*0,
     +          0, 5, 6, 8, 9,10,11,12,17,3*0,
     +          0, 5, 6, 8, 9,10,11,12,17,3*0,
     +          0, 5, 6, 8, 9,10,11,12,17,18,19,20/
c
c    pointers to the tz+p basis
c
      data ioffp/5, 0, 0,
     +          10,11, 0,
     +          10,14,15,
     +          16, 0, 0,
     +          18, 0, 0,
     +          18,19, 0,
     +          21, 0, 0 /
c
      data done/1.0d+00/
c     data dzero/0.0d+00/
c
c     set values for the current ahlrichs tz  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.gt. 4) ind=4
      if(nucz.gt.10) ind=5
      if(nucz.gt.11) ind=6
      if(nucz.gt.12) ind=7
      if(nucz.gt.18) then
         call caserr2(
     +          'ahlrichs tz basis sets only extend to argon')
      end if
c
c
      mxpass=kind(ind)
c                 all atoms have one polarisation function 
c                 except be and mg
      npol = 1
      if(nucz.eq.4 ) npol = 3
      if(nucz.eq.12 ) npol = 2
      if(opol) mxpass=mxpass+npol
c
      if(ipass.gt.mxpass) odone=.true.
      if(odone) go to 100
c
      kkkk = ipass-kind(ind)
      if (kkkk.le.0) then
       igauss = igau(ipass,ind)
       ng = ioff(ipass,ind) 
       ityp = itypes(ipass,ind)
      else
       igauss = igaup(kkkk,ind)
       ng = ioffp(kkkk,ind) 
       ityp = itypep(kkkk,ind)
      endif
c
c          do not scale midi, or mini inner shells
c
      scale=done
      if(ityp.eq.1) scale=scs
      if(ityp.eq.2) scale=scp
      if(ityp.eq.3) scale=scd
100   itype = ityp + 15
      return
      end
      subroutine bastyp(ztext,opol,oryd)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) ztext
c
      opol=.false.
      oryd=.false.
      if(index(ztext,'p').gt.0) opol = .true.
      if(index(ztext,'r').gt.0) oryd = .true.
c
      return
      end
      subroutine ddfour(e,s,p,d,n)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
c
c     ----- binning/curtiss [6s,4p,1d] contraction -----
c     ----- of dunning's (14s,11p,5d) primitives -----
c           j.comp.chem. 11, 1206-1216( 1990)
c
c     there are four typos in the binning/curtiss paper,
c     ge:p(18), as:p(21), br:s(3),  kr:e(11).
c     the values here are correct.
c
      nn = n-30
      s(11) = 1.0d+00
      s(12) = 1.0d+00
      s(13) = 1.0d+00
      s(14) = 1.0d+00
      p(24) = 1.0d+00
      p(25) = 1.0d+00
      go to (300,400,500,600,700,800),nn
c
c        gallium
c
  300 continue
      e(1) = 457600.0d+00
      e(2) = 68470.0d+00
      e(3) = 15590.0d+00
      e(4) = 4450.0d+00
      e(5) = 1472.0d+00
      e(6) = 541.3d+00
      e(7) = 214.8d+00
      e(8) = 88.81d+00
      s(1) = 0.00022088d+00
      s(2) = 0.0017224d+00
      s(3) = 0.0088661d+00
      s(4) = 0.035919d+00
      s(5) = 0.11446d+00
      s(6) = 0.28163d+00
      s(7) = 0.44810d+00
      s(8) = 0.29552d+00
      e(9) = 27.18d+00
      e(10)= 11.54d+00
      s(9) = 0.48025d+00
      s(10)= 1.14780d+00
      e(11)= 3.303d+00
      e(12)= 1.334d+00
      e(13)= 0.1947d+00
      e(14)= 0.07158d+00
      e(15)= 3274.0d+00
      e(16)= 765.4d+00
      e(17)= 241.6d+00
      e(18)= 89.39d+00
      e(19)= 36.36d+00
      e(20)= 15.60d+00
      p(15)= 0.0014743d+00
      p(16)= 0.013270d+00
      p(17)= 0.065384d+00
      p(18)= 0.22961d+00
      p(19)= 0.39929d+00
      p(20)= 0.43593d+00
      e(21)= 6.472d+00
      e(22)= 2.748d+00
      e(23)= 1.090d+00
      p(21)= 0.26105d+00
      p(22)= 0.48347d+00
      p(23)= 0.28721d+00
      e(24)= 0.2202d+00
      e(25)= 0.06130d+00
      e(26)= 59.66d+00
      e(27)= 17.10d+00
      e(28)= 6.030d+00
      e(29)= 2.171d+00
      e(30)= 0.6844d+00
      d(26)= 0.031960d+00
      d(27)= 0.16358d+00
      d(28)= 0.36720d+00
      d(29)= 0.45704d+00
      d(30)= 0.30477d+00
c
      e(31) = 0.141d0
      d(31) = 1.000d0
      return
c
c        germanium
c
  400 continue
      e(1) = 489000.0d+00
      e(2) = 73100.0d+00
      e(3) = 16640.0d+00
      e(4) = 4742.0d+00
      e(5) = 1569.0d+00
      e(6) = 577.0d+00
      e(7) = 229.0d+00
      e(8) = 94.81d+00
      s(1) = 0.000222d+00
      s(2) = 0.001728d+00
      s(3) = 0.008950d+00
      s(4) = 0.035905d+00
      s(5) = 0.11284d+00
      s(6) = 0.28748d+00
      s(7) = 0.43449d+00
      s(8) = 0.30683d+00
      e(9) = 29.22d+00
      e(10)= 12.45d+00
      s(9) = 0.47759d+00
      s(10)= 1.11445d+00
      e(11)= 3.642d+00
      e(12)= 1.502d+00
      e(13)= 0.2462d+00
      e(14)= 0.09209d+00
      e(15)= 3596.0d+00
      e(16)= 843.7d+00
      e(17)= 266.2d+00
      e(18)= 98.28d+00
      e(19)= 39.93d+00
      e(20)= 17.14d+00
      p(15)= 0.001442d+00
      p(16)= 0.012442d+00
      p(17)= 0.064694d+00
      p(18)= 0.218665d+00
      p(19)= 0.39335d+00
      p(20)= 0.43453d+00
      e(21)= 7.157d+00
      e(22)= 3.068d+00
      e(23)= 1.246d+00
      p(21)= 0.249975d+00
      p(22)= 0.46938d+00
      p(23)= 0.28896d+00
      e(24)= 0.2795d+00
      e(25)= 0.08340d+00
      e(26)= 70.18d+00
      e(27)= 20.07d+00
      e(28)= 7.059d+00
      e(29)= 2.553d+00
      e(30)= 0.8301d+00
      d(26)= 0.02877d+00
      d(27)= 0.15525d+00
      d(28)= 0.36814d+00
      d(29)= 0.46578d+00
      d(30)= 0.28648d+00
c
      e(31) = 0.202d0
      d(31) = 1.000d0
      return
c
c        arsenic
c
  500 continue
      e(1) = 526200.0d+00
      e(2) = 78310.0d+00
      e(3) = 17800.0d+00
      e(4) = 5070.0d+00
      e(5) = 1677.0d+00
      e(6) = 617.0d+00
      e(7) = 245.0d+00
      e(8) = 101.5d+00
      s(1) = 0.000220d+00
      s(2) = 0.00172d+00
      s(3) = 0.00889d+00
      s(4) = 0.03535d+00
      s(5) = 0.11226d+00
      s(6) = 0.28625d+00
      s(7) = 0.43406d+00
      s(8) = 0.30880d+00
      e(9) = 31.39d+00
      e(10)= 13.43d+00
      s(9) = 0.47299d+00
      s(10)= 1.12375d+00
      e(11)= 4.0d+00
      e(12)= 1.683d+00
      e(13)= 0.3003d+00
      e(14)= 0.1137d+00
      e(15)= 3909.0d+00
      e(16)= 924.0d+00
      e(17)= 291.8d+00
      e(18)= 107.5d+00
      e(19)= 43.62d+00
      e(20)= 18.73d+00
      p(15)= 0.00133d+00
      p(16)= 0.01140d+00
      p(17)= 0.05981d+00
      p(18)= 0.20511d+00
      p(19)= 0.37432d+00
      p(20)= 0.42141d+00
      e(21)= 7.841d+00
      e(22)= 3.391d+00
      e(23)= 1.405d+00
      p(21)= 0.21901d+00
      p(22)= 0.42955d+00
      p(23)= 0.26400d+00
      e(24)= 0.3441d+00
      e(25)= 0.1068d+00
      e(26)= 81.59d+00
      e(27)= 23.26d+00
      e(28)= 8.148d+00
      e(29)= 2.954d+00
      e(30)= 0.9827d+00
      d(26)= 0.02558d+00
      d(27)= 0.14492d+00
      d(28)= 0.36044d+00
      d(29)= 0.46255d+00
      d(30)= 0.26813d+00
c
      e(31) = 0.273d0
      d(31) = 1.000d0
      return
c
c        selenium
c
  600 continue
      e(1) = 560600.0d+00
      e(2) = 84010.0d+00
      e(3) = 19030.0d+00
      e(4) = 5419.0d+00
      e(5) = 1792.0d+00
      e(6) = 659.5d+00
      e(7) = 262.0d+00
      e(8) = 108.6d+00
      s(1) = 0.00022000d+00
      s(2) = 0.0017000d+00
      s(3) = 0.0088300d+00
      s(4) = 0.035080d+00
      s(5) = 0.11266d+00
      s(6) = 0.28471d+00
      s(7) = 0.43348d+00
      s(8) = 0.31140d+00
      e(9) = 33.66d+00
      e(10)= 14.45d+00
      s(9) = 0.46823d+00
      s(10)= 1.09261d+00
      e(11)= 4.378d+00
      e(12)= 1.876d+00
      e(13)= 0.3592d+00
      e(14)= 0.13670d+00
      e(15)= 4114.0d+00
      e(16)= 980.4d+00
      e(17)= 310.7d+00
      e(18)= 114.20d+00
      e(19)= 46.26d+00
      e(20)= 19.87d+00
      p(15)= 0.0011700d+00
      p(16)= 0.009850d+00
      p(17)= 0.051770d+00
      p(18)= 0.17869d+00
      p(19)= 0.32504d+00
      p(20)= 0.360890d+00
      e(21)= 8.309d+00
      e(22)= 3.598d+00
      e(23)= 1.522d+00
      p(21)= 0.21488d+00
      p(22)= 0.39629d+00
      p(23)= 0.22508d+00
      e(24)= 0.4032d+00
      e(25)= 0.12310d+00
      e(26)= 94.03d+00
      e(27)= 26.79d+00
      e(28)= 9.336d+00
      e(29)= 3.383d+00
      e(30)= 1.1450d+00
      d(26)= 0.020630d+00
      d(27)= 0.12227d+00
      d(28)= 0.31972d+00
      d(29)= 0.41738d+00
      d(30)= 0.23184d+00
c
      e(31) = 0.315d0
      d(31) = 1.000d0
      return
c
c        bromine
c
  700 continue
      e(1) = 574300.0d+00
      e(2) = 89070.0d+00
      e(3) = 20210.0d+00
      e(4) = 5736.0d+00
      e(5) = 1899.0d+00
      e(6) = 698.7d+00
      e(7) = 277.8d+00
      e(8) = 115.2d+00
      s(1) = 0.00022040d+00
      s(2) = 0.0016871d+00
      s(3) = 0.0087457d+00
      s(4) = 0.035597d+00
      s(5) = 0.11352d+00
      s(6) = 0.27973d+00
      s(7) = 0.44938d+00
      s(8) = 0.30023d+00
      e(9) = 35.97d+00
      e(10)= 15.50d+00
      s(9) = 0.45726d+00
      s(10)= 1.22978d+00
      e(11)= 4.771d+00
      e(12)= 2.077d+00
      e(13)= 0.4211d+00
      e(14)= 0.1610d+00
      e(15)= 4406.0d+00
      e(16)= 1042.0d+00
      e(17)= 332.1d+00
      e(18)= 121.9d+00
      e(19)= 49.24d+00
      e(20)= 21.16d+00
      p(15)= 0.0013766d+00
      p(16)= 0.012207d+00
      p(17)= 0.060190d+00
      p(18)= 0.22370d+00
      p(19)= 0.40047d+00
      p(20)= 0.44458d+00
      e(21)= 8.836d+00
      e(22)= 3.829d+00
      e(23)= 1.643d+00
      p(21)= 0.27245d+00
      p(22)= 0.49469d+00
      p(23)= 0.25684d+00
      e(24)= 0.465d+00
      e(25)= 0.1427d+00
      e(26)= 108.4d+00
      e(27)= 30.71d+00
      e(28)= 10.66d+00
      e(29)= 3.851d+00
      e(30)= 1.317d+00
      d(26)= 0.021521d+00
      d(27)= 0.13376d+00
      d(28)= 0.36673d+00
      d(29)= 0.49037d+00
      d(30)= 0.26749d+00
c
      e(31) = 0.338d0
      d(31) = 1.000d0
      return
c
c        krypton
c
  800 continue
      e(1) = 605700.0d+00
      e(2) = 90300.0d+00
      e(3) = 20920.0d+00
      e(4) = 5889.0d+00
      e(5) = 1950.0d+00
      e(6) = 718.2d+00
      e(7) = 285.4d+00
      e(8) = 118.6d+00
      s(1) = 0.00022997d+00
      s(2) = 0.0017432d+00
      s(3) = 0.0090459d+00
      s(4) = 0.036901d+00
      s(5) = 0.11743d+00
      s(6) = 0.28617d+00
      s(7) = 0.45408d+00
      s(8) = 0.28501d+00
      e(9) = 38.16d+00
      e(10)= 16.45d+00
      s(9) = 0.46465d+00
      s(10)= 1.24448d+00
      e(11)= 5.211d+00
      e(12)= 2.291d+00
      e(13)= 0.4837d+00
      e(14)= 0.1855d+00
      e(15)= 4678.0d+00
      e(16)= 1120.0d+00
      e(17)= 357.1d+00
      e(18)= 131.4d+00
      e(19)= 52.86d+00
      e(20)= 22.70d+00
      p(15)= 0.0013437d+00
      p(16)= 0.011821d+00
      p(17)= 0.058277d+00
      p(18)= 0.21854d+00
      p(19)= 0.39757d+00
      p(20)= 0.44119d+00
      e(21)= 9.547d+00
      e(22)= 4.167d+00
      e(23)= 1.811d+00
      p(21)= 0.26887d+00
      p(22)= 0.49381d+00
      p(23)= 0.25069d+00
      e(24)= 0.5337d+00
      e(25)= 0.1654d+00
      e(26)= 125.6d+00
      e(27)= 35.31d+00
      e(28)= 12.15d+00
      e(29)= 4.35d+00
      e(30)= 1.494d+00
      d(26)= 0.015695d+00
      d(27)= 0.10301d+00
      d(28)= 0.29968d+00
      d(29)= 0.41171d+00
      d(30)= 0.21548d+00
c
      e(31) = 0.318d0
      d(31) = 1.000d0
      return
      end
      subroutine dzsh(nucz,ipass,opol,oryd,itype,igauss,ng,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(16,7), itypes(16,7), ioff(16,7)
      dimension igaup(7), itypep(7), ioffp(7)
      dimension igaur(3,7), ityper(3,7), ioffr(3,7)
      dimension kind(7),kindr(7),kindp(7)
      data kind/2,5,5,6,10,16,11/
      data kindp/1,1,1,1,1,0,1/
      data kindr/0,3,3,3,3,0,0/
      data igau /
     + 3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 7,2,1,3,1,0,0,0,0,0,0,0,0,0,0,0,
     + 7,2,1,4,1,0,0,0,0,0,0,0,0,0,0,0,
     + 6,1,1,1,4,1,0,0,0,0,0,0,0,0,0,0,
     + 5,3,1,1,1,1,4,2,1,1,0,0,0,0,0,0,
     + 4,1,1,1,1,1,1,1,3,1,1,1,1,1,4,1,
     + 8,2,1,1,1,1,6,3,1,1,5,0,0,0,0,0 /
      data igaup /1,1,1,1,1,0,1/
      data igaur /
     + 0,0,0,
     + 1,1,1,
     + 1,1,1,
     + 1,1,1,
     + 1,1,1,
     + 0,0,0,
     + 0,0,0 /
      data itypes /
     + 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,
     + 1,1,1,1,1,1,2,2,2,2,3,0,0,0,0,0 /
      data itypep/2,3,3,3,3,0,3/
      data ityper /
     + 0,0,0,
     + 1,2,3,
     + 1,2,3,
     + 1,2,3,
     + 1,2,3,
     + 0,0,0,
     + 0,0,0 /
c
c    pointers to the dz basis
c
      data ioff/0, 3, 14*0,
     +          0, 7, 9,10,13,11*0,
     +          0, 7, 9,10,14,11*0,
     +          0, 6, 7, 8, 9,13,10*0,
     +          0, 5, 8, 9,10,11,12,16,18,19,6*0,
     +          0, 4, 5, 6, 7, 8, 9,10,11,14,15,16,17,18,19,23,
     +          0, 8,10,11,12,13,14,20,23,24,25, 5*0/
c
c    pointers to the dz+p basis
c
      data ioffp/4,14,15,14,20, 0,30/
c
c    pointers to the dz+r basis
c
      data ioffr/
     +  0, 0, 0,
     + 15,16,17,
     + 16,17,18,
     + 15,16,17,
     + 21,22,23,
     +  0, 0, 0,
     +  0, 0, 0 /
c
c     data dzero,done/0.0d+00,1.0d+00/
c
c     set values for the current dz  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=4
      if(nucz.eq. 4.or.nucz.eq.10) ind=3
      if(nucz.gt.10) ind=5
      if(nucz.gt.18) ind=6
      if(nucz.gt.30) ind=7
      if(nucz.gt.36) then
         call caserr2(
     +          'dunning dz basis sets only extend to krypton')
      end if
c
c
      npol = 0
      nryd = 0
      mxpass=kind(ind)
      if (opol) then
       npol = kindp(ind)
      endif
      if (oryd) then
       nryd = kindr(ind)
      endif
      mxpass=mxpass+npol+nryd
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not.odone) then
c
        kkkk = ipass-kind(ind)
        if (kkkk.le.0) then
         igauss = igau(ipass,ind)
         ng = ioff(ipass,ind) 
         ityp = itypes(ipass,ind)
        else
         kkkk = kkkk - npol
           if(kkkk.le.0) then
            igauss = igaup(ind)
            ng = ioffp(ind) 
            ityp = itypep(ind)
           else
            igauss = igaur(kkkk,ind)
            ng = ioffr(kkkk,ind) 
            ityp = ityper(kkkk,ind)
           endif
        endif
      endif
      itype = ityp + 15
      return
      end
      subroutine ffour(e,s,p,d,n)
      implicit real*8  (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension e(*),s(*),p(*),d(*)
c
c     ----- binning/curtiss [9s,6p,2d] contraction -----
c     ----- of dunning's (14s,11p,5d) primitives -----
c           j.comp.chem. 11, 1206-1216( 1990)
c
c     there are four typos in the binning/curtiss paper,
c     ge:p(18), as:p(21), br:s(3),  kr:e(11).
c     the values here are correct.
c
      nn = n-30
      s(7)  = 1.0d+00
      s(8)  = 1.0d+00
      s(9)  = 1.0d+00
      s(10)  = 1.0d+00
      s(11) = 1.0d+00
      s(12) = 1.0d+00
      s(13) = 1.0d+00
      s(14) = 1.0d+00
      p(21) = 1.0d+00
      p(22) = 1.0d+00
      p(23) = 1.0d+00
      p(24) = 1.0d+00
      p(25) = 1.0d+00
      d(30) = 1.0d+00
      d(31) = 1.0d+00
      go to (300,400,500,600,700,800),nn
c
c        gallium
c
  300 continue
      e(1) = 457600.0d+00
      e(2) = 68470.0d+00
      e(3) = 15590.0d+00
      e(4) = 4450.0d+00
      e(5) = 1472.0d+00
      e(6) = 541.3d+00
      s(1) = 0.00022088d+00
      s(2) = 0.0017224d+00
      s(3) = 0.0088661d+00
      s(4) = 0.035919d+00
      s(5) = 0.11446d+00
      s(6) = 0.28163d+00
      e(7) = 214.8d+00
      e(8) = 88.81d+00
      e(9) = 27.18d+00
      e(10)= 11.54d+00
      e(11)= 3.303d+00
      e(12)= 1.334d+00
      e(13)= 0.1947d+00
      e(14)= 0.07158d+00
      e(15)= 3274.0d+00
      e(16)= 765.4d+00
      e(17)= 241.6d+00
      e(18)= 89.39d+00
      e(19)= 36.36d+00
      e(20)= 15.60d+00
      p(15)= 0.0014743d+00
      p(16)= 0.013270d+00
      p(17)= 0.065384d+00
      p(18)= 0.22961d+00
      p(19)= 0.39929d+00
      p(20)= 0.43593d+00
      e(21)= 6.472d+00
      e(22)= 2.748d+00
      e(23)= 1.090d+00
      e(24)= 0.2202d+00
      e(25)= 0.06130d+00
      e(26)= 59.66d+00
      e(27)= 17.10d+00
      e(28)= 6.030d+00
      e(29)= 2.171d+00
      e(30)= 0.6844d+00
      d(26)= 0.031960d+00
      d(27)= 0.16358d+00
      d(28)= 0.36720d+00
      d(29)= 0.45704d+00
c
      e(31) = 0.141d0
      return
c
c        germanium
c
  400 continue
      e(1) = 489000.0d+00
      e(2) = 73100.0d+00
      e(3) = 16640.0d+00
      e(4) = 4742.0d+00
      e(5) = 1569.0d+00
      e(6) = 577.0d+00
      e(7) = 229.0d+00
      e(8) = 94.81d+00
      s(1) = 0.000222d+00
      s(2) = 0.001728d+00
      s(3) = 0.008950d+00
      s(4) = 0.035905d+00
      s(5) = 0.11284d+00
      s(6) = 0.28748d+00
      e(9) = 29.22d+00
      e(10)= 12.45d+00
      e(11)= 3.642d+00
      e(12)= 1.502d+00
      e(13)= 0.2462d+00
      e(14)= 0.09209d+00
      e(15)= 3596.0d+00
      e(16)= 843.7d+00
      e(17)= 266.2d+00
      e(18)= 98.28d+00
      e(19)= 39.93d+00
      e(20)= 17.14d+00
      p(15)= 0.001442d+00
      p(16)= 0.012442d+00
      p(17)= 0.064694d+00
      p(18)= 0.218665d+00
      p(19)= 0.39335d+00
      p(20)= 0.43453d+00
      e(21)= 7.157d+00
      e(22)= 3.068d+00
      e(23)= 1.246d+00
      e(24)= 0.2795d+00
      e(25)= 0.08340d+00
      e(26)= 70.18d+00
      e(27)= 20.07d+00
      e(28)= 7.059d+00
      e(29)= 2.553d+00
      e(30)= 0.8301d+00
      d(26)= 0.02877d+00
      d(27)= 0.15525d+00
      d(28)= 0.36814d+00
      d(29)= 0.46578d+00
c
      e(31) = 0.202d0
      return
c
c        arsenic
c
  500 continue
      e(1) = 526200.0d+00
      e(2) = 78310.0d+00
      e(3) = 17800.0d+00
      e(4) = 5070.0d+00
      e(5) = 1677.0d+00
      e(6) = 617.0d+00
      e(7) = 245.0d+00
      e(8) = 101.5d+00
      s(1) = 0.000220d+00
      s(2) = 0.00172d+00
      s(3) = 0.00889d+00
      s(4) = 0.03535d+00
      s(5) = 0.11226d+00
      s(6) = 0.28625d+00
      e(9) = 31.39d+00
      e(10)= 13.43d+00
      e(11)= 4.0d+00
      e(12)= 1.683d+00
      e(13)= 0.3003d+00
      e(14)= 0.1137d+00
      e(15)= 3909.0d+00
      e(16)= 924.0d+00
      e(17)= 291.8d+00
      e(18)= 107.5d+00
      e(19)= 43.62d+00
      e(20)= 18.73d+00
      p(15)= 0.00133d+00
      p(16)= 0.01140d+00
      p(17)= 0.05981d+00
      p(18)= 0.20511d+00
      p(19)= 0.37432d+00
      p(20)= 0.42141d+00
      e(21)= 7.841d+00
      e(22)= 3.391d+00
      e(23)= 1.405d+00
      e(24)= 0.3441d+00
      e(25)= 0.1068d+00
      e(26)= 81.59d+00
      e(27)= 23.26d+00
      e(28)= 8.148d+00
      e(29)= 2.954d+00
      e(30)= 0.9827d+00
      d(26)= 0.02558d+00
      d(27)= 0.14492d+00
      d(28)= 0.36044d+00
      d(29)= 0.46255d+00
c
      e(31) = 0.273d0
      return
c
c        selenium
c
  600 continue
      e(1) = 560600.0d+00
      e(2) = 84010.0d+00
      e(3) = 19030.0d+00
      e(4) = 5419.0d+00
      e(5) = 1792.0d+00
      e(6) = 659.5d+00
      e(7) = 262.0d+00
      e(8) = 108.6d+00
      s(1) = 0.00022000d+00
      s(2) = 0.0017000d+00
      s(3) = 0.0088300d+00
      s(4) = 0.035080d+00
      s(5) = 0.11266d+00
      s(6) = 0.28471d+00
      e(9) = 33.66d+00
      e(10)= 14.45d+00
      e(11)= 4.378d+00
      e(12)= 1.876d+00
      e(13)= 0.3592d+00
      e(14)= 0.13670d+00
      e(15)= 4114.0d+00
      e(16)= 980.4d+00
      e(17)= 310.7d+00
      e(18)= 114.20d+00
      e(19)= 46.26d+00
      e(20)= 19.87d+00
      p(15)= 0.0011700d+00
      p(16)= 0.009850d+00
      p(17)= 0.051770d+00
      p(18)= 0.17869d+00
      p(19)= 0.32504d+00
      p(20)= 0.360890d+00
      e(21)= 8.309d+00
      e(22)= 3.598d+00
      e(23)= 1.522d+00
      e(24)= 0.4032d+00
      e(25)= 0.12310d+00
      e(26)= 94.03d+00
      e(27)= 26.79d+00
      e(28)= 9.336d+00
      e(29)= 3.383d+00
      e(30)= 1.1450d+00
      d(26)= 0.020630d+00
      d(27)= 0.12227d+00
      d(28)= 0.31972d+00
      d(29)= 0.41738d+00
c
      e(31) = 0.315d0
      return
c
c        bromine
c
  700 continue
      e(1) = 574300.0d+00
      e(2) = 89070.0d+00
      e(3) = 20210.0d+00
      e(4) = 5736.0d+00
      e(5) = 1899.0d+00
      e(6) = 698.7d+00
      e(7) = 277.8d+00
      e(8) = 115.2d+00
      s(1) = 0.00022040d+00
      s(2) = 0.0016871d+00
      s(3) = 0.0087457d+00
      s(4) = 0.035597d+00
      s(5) = 0.11352d+00
      s(6) = 0.27973d+00
      e(9) = 35.97d+00
      e(10)= 15.50d+00
      e(11)= 4.771d+00
      e(12)= 2.077d+00
      e(13)= 0.4211d+00
      e(14)= 0.1610d+00
      e(15)= 4406.0d+00
      e(16)= 1042.0d+00
      e(17)= 332.1d+00
      e(18)= 121.9d+00
      e(19)= 49.24d+00
      e(20)= 21.16d+00
      p(15)= 0.0013766d+00
      p(16)= 0.012207d+00
      p(17)= 0.060190d+00
      p(18)= 0.22370d+00
      p(19)= 0.40047d+00
      p(20)= 0.44458d+00
      e(21)= 8.836d+00
      e(22)= 3.829d+00
      e(23)= 1.643d+00
      e(24)= 0.465d+00
      e(25)= 0.1427d+00
      e(26)= 108.4d+00
      e(27)= 30.71d+00
      e(28)= 10.66d+00
      e(29)= 3.851d+00
      e(30)= 1.317d+00
      d(26)= 0.021521d+00
      d(27)= 0.13376d+00
      d(28)= 0.36673d+00
      d(29)= 0.49037d+00
c
      e(31) = 0.338d0
      return
c
c        krypton
c
  800 continue
      e(1) = 605700.0d+00
      e(2) = 90300.0d+00
      e(3) = 20920.0d+00
      e(4) = 5889.0d+00
      e(5) = 1950.0d+00
      e(6) = 718.2d+00
      e(7) = 285.4d+00
      e(8) = 118.6d+00
      s(1) = 0.00022997d+00
      s(2) = 0.0017432d+00
      s(3) = 0.0090459d+00
      s(4) = 0.036901d+00
      s(5) = 0.11743d+00
      s(6) = 0.28617d+00
      e(9) = 38.16d+00
      e(10)= 16.45d+00
      e(11)= 5.211d+00
      e(12)= 2.291d+00
      e(13)= 0.4837d+00
      e(14)= 0.1855d+00
      e(15)= 4678.0d+00
      e(16)= 1120.0d+00
      e(17)= 357.1d+00
      e(18)= 131.4d+00
      e(19)= 52.86d+00
      e(20)= 22.70d+00
      p(15)= 0.0013437d+00
      p(16)= 0.011821d+00
      p(17)= 0.058277d+00
      p(18)= 0.21854d+00
      p(19)= 0.39757d+00
      p(20)= 0.44119d+00
      e(21)= 9.547d+00
      e(22)= 4.167d+00
      e(23)= 1.811d+00
      e(24)= 0.5337d+00
      e(25)= 0.1654d+00
      e(26)= 125.6d+00
      e(27)= 35.31d+00
      e(28)= 12.15d+00
      e(29)= 4.35d+00
      e(30)= 1.494d+00
      d(26)= 0.015695d+00
      d(27)= 0.10301d+00
      d(28)= 0.29968d+00
      d(29)= 0.41171d+00
c
      e(31) = 0.318d0
      return
      end
      subroutine svsh(nucz,ipass,opol,oryd,itype,igauss,ng,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(16,6), itypes(16,6), ioff(16,6)
      dimension igaup(6), itypep(6), ioffp(6)
      dimension igaur(3,6), ityper(3,6), ioffr(3,6)
      dimension kind(6),kindr(6),kindp(6)
      data kind/2,5,5,10,16,11/
      data kindp/1,1,1,1,0,1/
      data kindr/0,3,3,3,0,0/
      data igau /
     + 3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 7,2,1,3,1,0,0,0,0,0,0,0,0,0,0,0,
     + 7,2,1,4,1,0,0,0,0,0,0,0,0,0,0,0,
     + 5,3,1,1,1,1,4,2,1,1,0,0,0,0,0,0,
     + 4,1,1,1,1,1,1,1,3,1,1,1,1,1,4,1,
     + 8,2,1,1,1,1,6,3,1,1,5,0,0,0,0,0 /
      data igaup /1,1,1,1,0,1/
      data igaur /
     + 0,0,0,
     + 1,1,1,
     + 1,1,1,
     + 1,1,1,
     + 0,0,0,
     + 0,0,0 /
      data itypes /
     + 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,
     + 1,1,1,1,1,1,2,2,2,2,3,0,0,0,0,0 /
      data itypep/2,3,3,3,0,3/
      data ityper /
     + 0,0,0,
     + 1,2,3,
     + 1,2,3,
     + 1,2,3,
     + 0,0,0,
     + 0,0,0 /
c
c    pointers to the sv basis
c
      data ioff/0, 3, 14*0,
     +          0, 7, 9,10,13,11*0,
     +          0, 7, 9,10,14,11*0,
     +          0, 5, 8, 9,10,11,12,16,18,19,6*0,
     +          0, 4, 5, 6, 7, 8, 9,10,11,14,15,16,17,18,19,23,
     +          0, 8,10,11,12,13,14,20,23,24,25, 5*0/
c
c    pointers to the sv+p basis
c
      data ioffp/4,14,15,20, 0,30/
c
c    pointers to the sv+r basis
c
      data ioffr/
     +  0, 0, 0,
     + 15,16,17,
     + 16,17,18,
     + 21,22,23,
     +  0, 0, 0,
     +  0, 0, 0 /
c
c     data dzero,done/0.0d+00,1.0d+00/
c
c     set values for the current sv  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 3) ind=3
      if(nucz.gt.10) ind=4
      if(nucz.gt.18) ind=5
      if(nucz.gt.30) ind=6
      if(nucz.gt.36) then
         call caserr2(
     +          'dunning sv basis sets only extend to krypton')
      end if
c
c
      npol = 0
      nryd = 0
      mxpass=kind(ind)
      if (opol) then
       npol = kindp(ind)
      endif
      if (oryd) then
       nryd = kindr(ind)
      endif
      mxpass=mxpass+npol+nryd
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not.odone) then
c
        kkkk = ipass-kind(ind)
        if (kkkk.le.0) then
         igauss = igau(ipass,ind)
         ng = ioff(ipass,ind) 
         ityp = itypes(ipass,ind)
        else
         kkkk = kkkk - npol
           if(kkkk.le.0) then
            igauss = igaup(ind)
            ng = ioffp(ind) 
            ityp = itypep(ind)
           else
            igauss = igaur(kkkk,ind)
            ng = ioffr(kkkk,ind) 
            ityp = ityper(kkkk,ind)
           endif
        endif
      endif
      itype = ityp + 15
      return
      end
      subroutine tzvsh(nucz,ipass,opol,oryd,itype,igauss,ng,odone)
      implicit real*8 (a-h,p-w),integer (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension igau(20,7), itypes(20,7), ioff(20,7)
      dimension igaup(7), itypep(7), ioffp(7)
      dimension igaur(3,7), ityper(3,7), ioffr(3,7)
      dimension kind(7),kindr(7),kindp(7)
      data kind/3,8,8,11,11,20,17 /
      data kindp/1,1,1,1,1,1,1/
      data kindr/0,3,3,3,3,0,0/
      data igau /
     + 3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 6,2,1,1,1,3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 6,2,1,1,1,4,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
     + 6,3,1,1,1,1,4,2,1,1,1,0,0,0,0,0,0,0,0,0,
     + 6,3,1,1,1,1,5,2,1,1,1,0,0,0,0,0,0,0,0,0,
     + 5,1,1,1,1,1,1,1,1,1,4,1,1,1,1,1,1,1,4,1,
     + 6,1,1,1,1,1,1,1,1,6,1,1,1,1,1,4,1,0,0,0  /
      data igaup /1,1,1,1,1,1,1/
      data igaur /
     + 0,0,0,
     + 1,1,1,
     + 1,1,1,
     + 1,1,1,
     + 1,1,1,
     + 0,0,0,
     + 0,0,0 /
      data itypes /
     + 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,0,0,0,0,
     + 1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,
     + 1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,0,0,0  /
      data itypep/2,3,3,3,3,3,3/
      data ityper /
     + 0,0,0,
     + 1,2,3,
     + 1,2,3,
     + 1,2,3,
     + 1,2,3,
     + 0,0,0,
     + 0,0,0 /
c
c    pointers to the tzv basis
c
      data ioff/
     + 0, 3, 4,17*0,
     + 0, 6, 8, 9,10,11,14,15,12*0,
     + 0, 6, 8, 9,10,11,15,16,12*0,
     + 0, 6, 9,10,11,12,13,17,19,20,21,9*0,
     + 0, 6, 9,10,11,12,13,18,20,21,22,9*0,
     + 0, 5, 6, 7, 8, 9,10,11,12,13,14,18,19,20,21,22,23,24,25,29,
     + 0, 6, 7, 8, 9,10,11,12,13,14,20,21,22,23,24,25,29,3*0 /
c
c    pointers to the tzv+p basis
c
      data ioffp/5,16,17,22,23,30,30/
c
c
c    pointers to the tz+r basis
c
      data ioffr/
     +  0, 0, 0,
     + 17,18,19,
     + 18,19,20,
     + 23,24,25,
     + 24,25,26,
     +  0, 0, 0,
     +  0, 0, 0 /
c     data dzero,done/0.0d+00,1.0d+00/
c
c     set values for the current tzv  shell
c
c     return igauss = number of gaussians in current shell
c            ityp   = 1,2,3 for s,p,d shell
c            ng     = offset in e,cs,cp,cd arrays for current shell
c            scale  = scaling factor for current shell
c
c
      ind=1
      if(nucz.gt. 2) ind=2
      if(nucz.gt. 4) ind=3
      if(nucz.gt.10) ind=4
      if(nucz.gt.14) ind=5
      if(nucz.gt.18) ind=6
      if(nucz.gt.30) ind=7
      if(nucz.gt.36) then
         call caserr2(
     +          'dunning tzv basis set only extends to krypton')
      end if
c
c
      npol = 0
      nryd = 0
      mxpass=kind(ind)
      if (opol) then
       npol = kindp(ind)
      endif
      if (oryd) then
       nryd = kindr(ind)
      endif
      mxpass=mxpass+npol+nryd
c
      if(ipass.gt.mxpass) odone=.true.
      if(.not.odone) then
c
        kkkk = ipass-kind(ind)
        if (kkkk.le.0) then
         igauss = igau(ipass,ind)
         ng = ioff(ipass,ind) 
         ityp = itypes(ipass,ind)
        else
         kkkk = kkkk - npol
           if(kkkk.le.0) then
            igauss = igaup(ind)
            ng = ioffp(ind) 
            ityp = itypep(ind)
           else
            igauss = igaur(kkkk,ind)
            ng = ioffr(kkkk,ind) 
            ityp = ityper(kkkk,ind)
           endif
        endif
      endif
      itype = ityp + 15
      return
      end
      subroutine ghosts
c
c...  eliminate ghost atoms according to specification in /ghost/
c     ghost (nuct = 3)
c        all centers without a charge but with basis functions
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
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
c from March 97, amass is used only for user-input
c mass values
c czanr contains a copy of czan made in basis.m after reorganisation and
c before pseudo potential; so real charges in real order 
c beware, infob is explicit in direct.m and  dirrpa.m
c
      real*8 czann, czanr, cat, amass, cnew
      integer nonsym, map80
      logical ozmat
      common/infob/czanr(maxat),czann(maxat),cat(3,maxat),amass(maxat),
     +             cnew(maxat3,3),nonsym,map80(maxat),ozmat

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
c...  counterpoise common
c
      character *8 ztag00
      integer maxghost,nghost
      parameter (maxghost=100)
      common/ghost/ ztag00(maxghost)
      common/ghostn/ nghost
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      write(iwr,600) (ztag00(i),i=1,nghost)
600   format(1x,'** Ghost centers :',(t20,10(1x,a8)))
      write(iwr,'(1x)')
c
      do i=1,nghost
         do j=1,nonsym
            omatch = .true.
            k = index(ztag00(i)//' ',' ') 
            l = index(ztag(j)//' ',' ') 
            maxkl = max(k,l)-1
            do m=maxkl,1,-1
               if (ztag00(i)(k:k).ne.'*'.or.k.ne.1) k = k - 1
               l = l - 1
               if (ztag(j)(l:l).ne.ztag00(i)(k:k)
     1             .and.ztag00(i)(k:k).ne.'*') then
                    omatch = .false.
               end if
            end do
            if (l.ne.1.or.k.ne.1) omatch = .false.
            if (omatch) then
               nuct(j) = 3
               ipseud(j) = 0
               czann(j) = 0
            end if
         end do
      end do
c
      return
      end
      subroutine ver_basis(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/basis.m,v $
     +     "/
      data revision /"$Revision: 6135 $"/
      data date /
     +  "$Date: 2010-05-31 13:43:04 +0200 (Mon, 31 May 2010) $"
     +  /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
