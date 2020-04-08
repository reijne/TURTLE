c 
c  $Author: hvd $
c  $Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
c  $Locker:  $
c  $Revision: 6090 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mainci.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   blockcas   =
c ******************************************************
c ******************************************************
      block data mcscf
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
      common/mlngth/nana(6)
      common/stak  /btri,mmmtri(4)
      common/stak2/ ctri,nnntri(4)
      common/bufc/cccc(510),icccc(4)
      common/testop/oto(20)
c
      integer nwt, ispavs
      real*8 weight
      common /avstat/ nwt,ispavs(3),weight(5)
c
c
      integer n20cas
      integer ifsimu, isimul, iam, noci, ipople
      parameter (n20cas = 100)
      common /simul/ ifsimu,isimul(n20cas),iam(n20cas),noci(n20cas),
     +               ipople
c
c
      integer igrad,ifout1,ifout2,ifmola,if2,if4,iflagr
      integer isecdm,isecda,iseclm,isecla
      integer idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
      common /dm/ igrad,ifout1,ifout2,ifmola,if2,if4,iflagr,
     +            isecdm,isecda,iseclm,isecla,
     +            idmmc,iblkmm,lblmm,ifiled,iwordd,iseclg
c
c
      real*8 fudgit, fmax, swnr, swsimu
      integer nhesnr
      logical ojust, osuped
      common /ciconv/ fudgit(10),fmax,ojust,osuped,swnr,swsimu,nhesnr
c
c
      real*8 cccnv
      logical oconv, oconvf
      integer icanon
c
      common /gms_finish/ cccnv,oconv,icanon,oconvf
c
      integer mode, nrever, ifhes, isort, iaugm, isp, ibuffn
      common /ctrl/ mode,nrever,ifhes,isort(n20cas),iaugm,
     +              isp,ibuffn
c
      common/degen/ifdeg,idegen(100),ophase
c
      real*8 ssz
      logical ofdav, odoubl
      integer ibuffc, nbuffc, nbasec
      integer lencon,lenint
      common /caspar/ ssz,
     +                ofdav,ibuffc,nbuffc,nbasec,
     +                lencon,lenint,
     +                odoubl
c
      common/cone/m,na,nb,mf,nm1a,nm1b,mta1,mtb1,n82,iperm(8,8),ims(8)
      common/exc   /km,kkk,ql(5),idom(5),iocca(5,15),ioccb(5,15),irestu
     * , ifstr
      common/qpar/
     1ma,n2e,ntype,ncore,nact,nprim,nsec,nprimp,naa,nab,nst,n1e,isym,
     2mults,macro,maxc,itypci,itydrt,iduqpar
      common/cipar/clevc,accci,timeci,nitsci,nfudge,ncfl,ncfu,nprin
      common/blbpar/accblb,clev,timblb,nitblb,nprinb
      common/dims /ippd(12)
      common/popnos/ppo(300)
      common/add/ippa(27)
c
c... ic2e, ic3e and ic4e have room for 70 orbitals. Previously 15
c... orbitals could be accommodated for casscf. Now more are needed
c... for vb.
c
      integer ic1e, ic2e, ic3e, ic4e,maxcas
      parameter (maxcas=150)
      common /qice / ic1e(maxorb,2),ic2e(maxcas*(maxcas+1)/2+1),
     #               ic3e(maxcas*(maxcas+1)/2),ic4e(maxcas*(maxcas+1)/2)
c
      common/actlen/lena(108),ilifp(maxorb),lena2(24)
      common/potn/ppp(511)
      common/tab3/itab(4)
      common /drtlnk/ idrtln(6),labdrt(maxorb)
      common /loops / rloop(5),loops (15)
      common /form  / form  (8)
      common /valen / valen (2)
      common /casscg/ icas(4),icas2(maxorb)
      common /junkg / junkg (18)
      common/junke/maxt,ires,ipass,
     1nteff,mpass1,mpass2,lentri,
     2nbuck,mloww,mhi,ntri,jacc,ijunke(8)
      data oto/20*.false./
      data ipople/0/
      data itypci/5/
      data odoubl/.true./
      data nwt/1/
      data ifsimu,isimul,iam/0,n20cas*0,n20cas*0/
      data noci/n20cas*0/
      data iaugm/0/
      data fudgit,ojust,fmax/10*1.0d-5,.true.,.15d0/
      data osuped,swnr,swsimu,nhesnr/.true.,0.05d0,0.0d0,1/
      data ophase/.false./
      data macro/1/
cjvl  total # must be n20cas (=100)
      data isort/4*0,2,3*1,2,3*1,2,3*1,2,3*1,
     +           2,3*1,2,3*1,2,3*1,2,3*1,2,3*1,
     +           2,3*1,2,3*1,2,3*1,2,3*1,2,3*1,
     +           2,3*1,2,3*1,2,3*1,2,3*1,2,3*1,
     +           2,3*1,2,3*1,2,3*1,2,3*1,2,3*1/
      data nrever/0/
      data mode/0/
      data ifdeg,idegen/101*0/
      data ssz/1.3d0/
      data iperm/1,2,3,4,5,6,7,8,
     1           2,1,4,3,6,5,8,7,
     2           3,4,1,2,7,8,5,6,
     3           4,3,2,1,8,7,6,5,
     4           5,6,7,8,1,2,3,4,
     5           6,5,8,7,2,1,4,3,
     6           7,8,5,6,3,4,1,2,
     7           8,7,6,5,4,3,2,1/
      data idom(1),km,iocca,ioccb,irestu/1,1,151*0/
      data maxc/20/
      data icanon/3/
      data oconv/.false./
      data cccnv/1.0d-4/
      data timblb,nitblb,accblb,clev,nprinb/20.d0,40,
     *  1.0d-05,2.0d0,1000/
      data timeci,nitsci,accci,nfudge,ncfl,ncfu,nprin,clevc/
     *  30.d0,20,1.d-4,30,0,0,32000,0.0d0/
      end
c ******************************************************
c ******************************************************
c             =   blockci   =
c ******************************************************
c ******************************************************
      block data  mrdci
      implicit real*8  (a-h,p-z),integer   (i-n),logical     (o)
      integer  mms,ms,nd32,n32m,mmms
      integer*4 mms_i,ms_i,mmms_i,nd32_i,n32m_i
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
      integer n24, magic, n32, n64, n8, n16, n48
      common /cdcryi/ n24,magic,n32,n64,n8,n16,n48
c
      common /cdcryz/ mms,ms,nd32,n32m,mmms
c
c    following common block definitions are to allow
c    gamess blocks to replace those in direct-ci
c    as follows
c    locol->junk; lhsrhs->three; symbol->scra; loco->lsort
c
c    to minimize core in scf gamess revert to definitions
c    in member blockall of this pds
c
      common/lsort/jlhs(mxcan1),jrhs(mxcan1),jr(mxcan1),
     *xj(mxcan1),
     *yj(mxcan1),nt9(16),jl9(16)
      common /scra/ real2(3400),ireal2(6,3400)
      common /three/ int3(4500)
c     common /junk / real4(4000),int4(9224)
c
      logical latorb, mspin1, ispind
      integer ispacg, isping, iprwh, iprcon
      common /natorb/ ispacg,isping,latorb,mspin1,ispind,
     +                iprwh,iprcon
c
      common/symchk/crtsym(2),nrep,osymd,osadap
      common/pager/ipage(20)
c
      integer khgot2,lump1,lump2,lump3,lump4,lump5,lump6,lump7,
     + lump8,lump9,lump10,lump11,lump12,lump13,lump14,lump15,lump16,
     + lump17,lump18,lump19,lump20,lump21,lump22,lump23,lump24,lump25,
     + lump26,lump27,lump28,
     + mgot0,mgot1,mgot2,lump29,lump30,lump31,lump32,
     + lump33,lump38,lump39,lump40,lump41,lump47,lump48,lump49,lump50,
     + lump51,lump52,lump53,lump54,lump55,lump56,lump57,lump58,lump99,
     + kumpc,kumpe,lump59
      common/outctl/khgot2,lump1,lump2,lump3,lump4,lump5,lump6,lump7,
     + lump8,lump9,lump10,lump11,lump12,lump13,lump14,lump15,lump16,
     + lump17,lump18,lump19,lump20,lump21,lump22,lump23,lump24,lump25,
     + lump26,lump27,lump28,
     + mgot0,mgot1,mgot2,lump29,lump30,lump31,lump32,
     + lump33,lump38,lump39,lump40,lump41,lump47,lump48,lump49,lump50,
     + lump51,lump52,lump53,lump54,lump55,lump56,lump57,lump58,lump99,
     + kumpc,kumpe,lump59
c
c
      integer mestri
      integer minstr, maxstr, nesvac, nesuno, nesdue
      integer nesref, mxstri
      logical screan
      common /hold/ mestri(2,10,4),minstr(10,4),maxstr(10,4),
     +              nesvac,nesuno,nesdue,nesref,mxstri,screan
c
c
      real*8 tvec
      integer ivec, nvec, mxvec, iselvc, nselvc, istvc
      integer iprvc
      common /trial/ tvec(20),ivec(20),nvec,mxvec,
     +               iselvc,nselvc,istvc,iprvc
c
      common/erg/ptnuc(nd200+22),ipotnc(5)
c
      real*8 timbeg,  tmlast,  timsec
      real*8 wtimbeg, wtmlast, wtimsec
      real*8 cpwa
      integer nodel
      common/timanb/timbeg,tmlast,timsec(14),
     *wtimbeg,wtmlast,wtimsec(14),cpwa,nodel
c
c
      character *8 rname, sname, text11, text22, anumt, bnumt
      common/cic/rname(8),sname,text11(14),text22(14),anumt(8),bnumt(8)
c
      common/dopey/nshif(127)
      common/xynumb/xnumb(6)
      common/expanc/con12(10),npopbb(2240)
      common/spew/rmodel(6)
c
      integer iro, iky, niky, m1e, m2e, m1ia, m2coul
      integer m2exch, m2fock, m2psup, m2qsup
      integer mbase, nbase, npairi, norbi, mbasi, mnbasi
      integer iap, mbfock, mbexch, mbcoul, mbpsup, mbqsup
      integer norbas, m2exab, m2em1e
      common /helpr/ iro(nd200),iky(nd200),niky(8),m1e,m2e,m1ia,m2coul,
     +   m2exch,m2fock,m2psup,m2qsup,
     +   mbase(288),nbase(8),npairi(36),norbi(8),mbasi(8),mnbasi(8),
     +   iap(62),mbfock(62),mbexch(36),mbcoul(36),mbpsup(36),mbqsup(36),
     +   norbas(nd200),m2exab,m2em1e
c
      common/symcon/nornor(27)
c
      integer nuse, ngot
      common /corctl/ nuse, ngot
c
      common /ccntl/ nnref(15)
      common/stoctl/iiiiii(32)
      common/moco/imoco(64),rmoco(64)
      common/cntrl/nval(14)
c
      integer lcall,  mrstor
      logical modev, moded, modest, instor, loijka, msort, loijab
      logical loiajb, lvar
      common/diactl/modev,moded,modest,instor,
     +              loijka,msort,loijab,loiajb,lcall,mrstor,lvar
c
c
      integer norbe, norbsi, norbtr, norbsq
      integer ibassi, ibastr, ibassq, mubsi, nubsi
      integer mbassi, mbassq, ibass3
      integer n127, n127sq, nincr, maxsq, mubsq, madsi, madsq, nadsq
      integer nubtr, nubsq, nvmax
c
      common /symci/ norbe(8),norbsi(8),norbtr(8),norbsq(8),
     + ibassi(8,8),ibastr(8,8),ibassq(8,8),mubsi,nubsi,
     + mbassi(8,8),mbassq(8,8),ibass3(8,8),
     + n127,n127sq,nincr,maxsq,mubsq,madsi(8,8),madsq(8,8),nadsq(8,8),
     + nubtr,nubsq,nvmax
c
c
      integer nconst, ncond, maxpa, maxpat, maxsin, maxepa
      integer mxdvpa, mxddpa, nspips, nspipt, mxsspa, mxttpa
      integer mxstpa, nsitot, ntrtot
      integer mx12pa, mxvvpa, mxv2pa, maxpsq, max12
      common /auxh/ nconst(8),ncond(8),maxpa(4),maxpat,maxsin(3),
     +              maxepa(4),mxdvpa,mxddpa,nspips(8),nspipt(8),
     +              mxsspa,mxttpa,mxstpa,nsitot(8),ntrtot(8),
     +              mx12pa,mxvvpa,mxv2pa,maxpsq,max12
c
      common /table/ nirrrr(68)
      common/disktl/nxblk(204)
c
      real*8 valp
      integer link, nwbuf, lnumb, iblkp, indxi, indxj, indxk
      integer indxl, istakp
      common /presrt/ link(78),nwbuf(78),lnumb(78),ibasen(78),
     +                iblkp,indxi,indxj,indxk,indxl,istakp,valp
c
       common/mapp/mapei(nd200*2+1)
c
      real*8 critcep, crycep
      integer ifcepa, icepa , lcepaul, lsinsh, itcepa, lprcep, indcep
      integer nrfcsf, mcycep, ibrst,icorr_mr,mrdoc
      common/ccepa/critcep,crycep,ifcepa,icepa,lcepaul,lsinsh,itcepa,
     1                lprcep,indcep,nrfcsf,mcycep,ibrst,icorr_mr,mrdoc
c
c
      real*8 cradd, weighb
      integer iex1,iex2
      integer iblkcc, iblkcr, iblkcv, ncv, nbuenk
      integer mbuenk
      logical ifbuen
      common /comrjb/ cradd(3),weighb,iex1,iex2,ifbuen,
     +                iblkcc,iblkcr,iblkcv,ncv,nbuenk,mbuenk
c
c
      common/count/ iwrit,irea,iblkw,iblkr,iou2,lenou2,inn2,lenin2
      common/intbu2/ iblff(6),gggg(511)
      common/intbuf /mintb(4),gintb(511)
      common/couple /coupl(511),icoupl
      common/entryp /rorf(maxorb+1),irorbf
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      common /mcopt / var,varc,thzr,done,dtwo
      common /mccore/ intrel,lword,ltop,lmax,lmin,mreal
      integer mults, symaos, istart, mfin, nsymm ,nprm, ncor
      integer nactt, nsecc, ic1d, ne, itypea, ilifa, lentr
      integer lentrc, lentri, lenprm, lenrec, lensq, maxprm, maxbas
      integer lensqr, iorbsm, nrot, nrottu, nrotit, nrotia, nrotta
      integer irottu, lentca, lentra, nsymao, nbasao, ltri, ltrimo
      common /syminf/ mults(8,8),symaos(maxorb),istart(8),mfin(8)
     +         ,nsymm(8),nprm(8),ncor(8),nactt(8),nsecc(8)
     +         ,ic1d,ne,itypea(31),ilifa(31)
     +         ,lentr(8),lentrc(8),lentri
     +         ,lenprm,lenrec,lensq,maxprm,maxbas
     +         ,lensqr(8),iorbsm(maxorb)
     +         ,nrot,nrottu,nrotit,nrotia,nrotta,irottu(465)
     +         ,lentca(8),lentra(8),nsymao(8),nbasao,ltri,ltrimo
      real*8  radius,trust1,tfac1,trust2,tfac2,sparse,conv
      real*8  econv,sconv,glast,glast2,elast,elast2,enext,slast
      real*8  weight,auto1,auto2,auto3,gfak1,gfak2,gfak3
      real*8  drmax,varmin,disvar,varmax,copvar,select,augvar
      real*8  cishft,drdamp,ciacc,thrdiv,ciderr,sparec
      integer itmaxr,igvec,ntexp,ipri,maxdis,idstrt
      integer idstep,maxaug,idsci,maxci,icstep,icimax,irdamp
      integer maxito,maxitc,maxrep,nitrep,iuprod,igwgt,icstrt
      integer iroot1,icinat,icimx1,icimx2
      integer nbasis,ncoremc,nact,nci,norbr,nprim,nsec,nst
      integer nfreez,ifreez,nprimp,nirrr,lenbas,nblkq,nstate
      integer itype,ifzsym,numa,num2,num4,num6,num3,iblk3,isec,iblkq
      integer isecd,iblkdmc,iblkc,isecn,iblkn,isecnc,isecp,numft,nblkd1
      integer nblko1,nblkd2,nblko2,numscr,iblka,iblk2,iblk4,iblk6
      integer iblk8,iblkft,iffmtp,ifsort,ifwvfn,ifanal,itrsfm
      integer maxcyc,iter,itinfo,ianal,iprint,lprint
      integer idumpo,icani,icant,icana,icang,n1elec,i1elec
      integer iexc,iexcv,nref,bfkey,irestr,ispare,ibfcod
      common/multic/radius,trust1,tfac1,trust2,tfac2,sparse,conv
     +             ,econv,sconv,glast,glast2,elast,elast2,enext,slast
     +             ,weight(5),auto1,auto2,auto3,gfak1,gfak2,gfak3
     +             ,drmax,varmin,disvar,varmax,copvar,select,augvar
     +             ,cishft,drdamp,ciacc,thrdiv,ciderr,sparec(2)
     +             ,itmaxr,igvec,ntexp,ipri,maxdis,idstrt
     +             ,idstep,maxaug,idsci,maxci,icstep,icimax,irdamp
     +             ,maxito,maxitc,maxrep,nitrep,iuprod,igwgt,icstrt
     +             ,iroot1,icinat,icimx1,icimx2
     +             ,nbasis,ncoremc,nact,nci,norbr,nprim,nsec,nst
     +             ,nfreez,ifreez(8),nprimp
     +             ,nirrr,lenbas,nblkq,nstate,itype(maxorb)
     +             ,ifzsym(mcfzc)
     +             ,numa,num2,num4,num6,num3,iblk3,isec,iblkq,isecd
     +             ,iblkdmc,iblkc,isecn,iblkn,isecnc,isecp,numft,nblkd1
     +             ,nblko1,nblkd2,nblko2,numscr,iblka,iblk2,iblk4,iblk6
     +             ,iblk8,iblkft,iffmtp,ifsort,ifwvfn,ifanal,itrsfm
     +             ,maxcyc,iter,itinfo(40),ianal,iprint,lprint
     +             ,idumpo,icani,icant,icana,icang,n1elec,i1elec(20)
     +             ,iexc,iexcv,nref,bfkey(31),irestr(100)
     +             ,ispare(5),ibfcod(1356)
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      common/rotat/akkkk(60),nswap
      common/linkmc/lnkmc(12)
      common/mcaddr/iaddr(7)
      character*3 codes
      common /drtcoc/ codes(9)
      integer dela,delb,delele,virtul,occupd,valocc,rescor,resvir
     >       ,frozen,valvir,opensh,multi,speshl,multrf,valenc
     >       ,fzc,fzv,cor,vir,doc,uoc,alp,bet,spe
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     1,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     2,               valvir,opensh,multi,speshl,multrf,valenc
     3,   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
c *** mrd/ci common blocks
      common/aplus/ra1(900),ra2(751),ra3(255),ira1(300),
     *   ira2(6,256),ira3(2040),ira4(12,8),ira5(2,9),ira6(3,4)
      common/linkmr/linkm(514)
      common/jany/jerk(10),jbnk(10),jbun(10)
      common/aaaa/icom(21)
      common/stak3/ dtri,llltri(4)
      common/bufd/dddd(510),idddd(4)
      common/scrtch/cisd(185000)
      dimension mms_i(2),ms_i(2),mmms_i(2),nd32_i(2),n32m_i(2)
      equivalence (mms,mms_i),(ms,ms_i),
     * (mmms,mmms_i),(nd32,nd32_i),(n32m,n32m_i)
c
      data icom /1,2,1,3,3,1,4,6,4,1,5,10,10,5,1,6,15,20,15,6,1/
      data jerk/0,30,118,158,248,274,321,0,0,0/
      data jbnk/4,17,7,15,3,5,7,0,0,0/
      data jbun/ 3,3,2,2,1,1,1,0,0,0/
c
      data n24,magic,n32,n64,n8,n16,n48/24,4,32,64,8,16,48/
c     /cic/
      data anumt/'one','two','three','four',
     *           'five','six','seven','eight'/
      data bnumt/'1','2','3','4',
     *           '5','6','7','8'/
      data mms_i/z'ffffff',0/,ms_i/z'ffff',0/
      data mmms_i/z'ff',0/
      data nd32_i/z'00000000',z'aaaa0000'/
      data n32m_i/z'aaaaaaaa',z'aaaaaaaa'/
      end
      block data saveci
c
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c...  block data to initialise variables that could be better
c...  handled by save statements in subroutines
c...  unfortunately some compilers have trouble with those
c...  (variables have been checked not to cause interfeence
c...   in any of the other routines mentioned below)
c...   this is for ci-part
c
c...  subroutines affected :
c     initex    (casb) : common/cassav/mfrm
c
      common/cassav/mfrm
c
      end
      subroutine ver_mainci(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mainci.m,v $
     +     "/
      data revision /"$Revision: 6090 $"/
      data date /"$Date: 2009-11-04 16:57:24 +0100 (Wed, 04 Nov 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
