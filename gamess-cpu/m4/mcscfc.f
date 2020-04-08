c 
c  $Author: jmht $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfc.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   mcfmtp     =
c ******************************************************
c ******************************************************
      subroutine editft(iq,iwrite)
c...   this routine deletes interacting singles from off diagonal
c     parts of formula tape
c..   the code is inefficient and must not be used on large cases
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
      integer g,packft
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     >               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      common /scra / out(511)
      parameter (maxsng=511)
      dimension iq(1),isingl(maxsng),ioccn(31)
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
c...  scan on-diagonal part of formula tape in order to determine
c...  reference configuration
      nocc = nact*nci
      iocc = icori(nocc)
      do 270 i=1,nact
      ioc = 0
      if (ibfcod(i).eq.5) ioc=-2
      if (ibfcod(i).eq.7 .or. ibfcod(i).eq.8) ioc=-1
      ioccn(i)=-ioc
      do 270 j=1,nci
270   iq(iocc-1+j+(i-1)*nci)=ioc
      call rdedx(surd,511,iblkft,numft)
      do 280 i=1,nblkd1
      call find(numft)
      call get (g,nw)
      call unpkft(nw)
      do 280 iw=1,nw
      io = iorb(iw)
      iii = ii(iw)
280    iq(iocc-1+(io-1)*nci+iii) = iq(iocc-1+(io-1)*nci+iii)
     >         + nint(surd(ival(iw)))
c      call outisq (iq(iocc),nci,nci,nact,'iocc')
c...  now the ref config should have zeroes in iocc
      do 190 iref=1,nci
      io=0
      do 191 j=1,nact
191   io = io + iabs(iq(iocc-1+(j-1)*nci+iref))
      if (io.eq.0) goto 192
190   continue
      call caserr('help')
192   continue
c
c...  now form list of interacting singles
      nsing=0
      do 20 i=1,nblko1
      call find(numft)
      call get (g,nw)
      call unpkft(nw)
      do 21 iw=1,nw
      if (ii(iw).ne.iref.and.jj(iw).ne.iref) goto 21
      iii = ii(iw)
      if (iii.eq.iref) iii=jj(iw)
      jmin=0
      do 22 j=1,nsing
      if (isingl(j)-iii) 23,21,22
23    jmin=j
22    continue
      nsing=nsing+1
      if (nsing.gt.maxsng) call caserr('overflow in editft')
      do 25 j=nsing,jmin+2,-1
25    isingl(j) = isingl(j-1)
      isingl(jmin+1)=iii
21    continue
20    continue
      write(iwrite,643)iref,(ioccn(k),k=1,nact)
643   format(' the following configurations interact through single'
     <,' excitations'/' with the reference configuration number',i7
     >,' =',31i2/)
      do 644 j=1,nsing
644   write(iwrite,645)isingl(j),
     > (iq(iocc-1+(k-1)*nci+isingl(j))+ioccn(k),k=1,nact)
645   format(i8,5x,31i3)
      call corlsi(iocc)
c
c...  now we can edit the formula tape
c....  first, copy the whole tape to num4
      call search (iblk4,num4)
      call search (iblkft,numft)
      do 143 i=1,1+nblkd1+nblko1+1+nblkd2+nblko2+1
      call find (numft)
      call get (g,nw)
143   call put (g,nw,num4)
      do 89 izero=1,511
      if ( dabs(surd(izero)).lt.1.0d-9) goto 91
89    continue
      call caserr('nasty in editft')
91     continue
      call search (1+iblkft,numft)
      call search (1+iblk4,num4)
      do 145 i=1,nblkd1+nblko1+1+nblkd2+nblko2+1
      call find(num4)
      call get (g,nw)
      call unpkft(nw)
      do 146 iw=1,nw
      do 147 j=1,nsing
      if (ii(iw).eq.isingl(j) .or. jj(iw).eq.isingl(j)) goto 148
147   continue
      goto 146
148   g(iw) = packft(ii(iw),jj(iw),izero,iorb(iw),jorb(iw),
     +               korb(iw),lorb(iw) )
146   continue
145   call put (g,nw,numft)
      write(iwrite,853)
853   format(' the formula tape has been edited to remove interacting si
     >ngle excitation configurations')
      return
      end
      subroutine fmtp(iq,iwrite)
      implicit integer (a-z)
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
      external mcindx
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
      integer orbfrm,symorb,bfnum,bfsym,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      character*3 codes
      common /drtcoc/ codes(9)
      integer dela,delb,delele,virtul,occupd,valocc,rescor,resvir
     >       ,frozen,valvir,opensh,multi,speshl,multrf,valenc
     >       ,fzc,fzv,cor,vir,doc,uoc,alp,bet,spe
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     1,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     2,               valvir,opensh,multi,speshl,multrf,valenc
     3,   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      dimension iq(*)
c
      call accnt('formtape',1)
      nsym=nirrr
      nrefs=nref
      nbf=nact
      norbs=nact
      optio(2)=iexc
      optio(3)=iexcv
      optio(4)=1
      optio(5)=0
      optio(6)=6
      optio(7)=0
      if (optio(3).eq.0.and.optio(2).ne.32767)write(iwrite,10) optio(2)
10    format(' excitation level =',i3)
      if (optio(3).ne.0) write(iwrite,20) optio(2),optio(3)
20    format(' general excitation level =',i3
     1      /' valence excitation level =',i3)
      if (optio(2).eq.32767) write(iwrite,30)
30    format(' full excitation in active space')
      maxb=16
      mmin = icori(0)
      mmax = icorim()
      do 40 i=mmin,mmax
40    iq(i)=0
      nrefs = max(nrefs,1)
c
      ntype = icori(ntypes)
      junk = ntypes*nbf
      bfnum = icori(junk)
      junk = ntypes*nsym
      numsym = icori(junk)
      bfsym = icori(nbf)
      orbtbf = icori(norbs)
      iout   = icori(nbf)
      orbsym = icori(norbs)
      call mcreor(iq(ntype),iq(bfnum),iq(numsym),iq(bfsym)
     1,           ibfcod,iq(orbtbf),iq(iout),iq(orbsym))
      spec=0
      nspc=0
      if (spec.gt.0) nspc=4**spec
      nlevs=norbs+1
c      levfrm = 0
      nlevoc=nlevs-levfrm+1
      levpt  = icori(nlevs)
      levnr  = icori(nlevs)
      nelecs = icori(nrefs)
      spc    = icori(nspc)
      nrowsp = icorim()/(8+nrefs) - 1
      nrws4p = nrowsp*4
      a      = icori(-nrowsp)
      bp     = icori(-nrowsp)
      sp     = icori(-nrowsp)
      arcp   = icori(-nrws4p)
      nlwksp = icori(-nrowsp)
      junk   = nrefs*nrowsp
      x      = icori(-junk)
      if(spec.gt.0) 
     +              call caserr('special codes not implemented')
      call mcpald(iq(bfsym),bfkey,ibfcod
     1,           iq(a),iq(bp),iq(sp),iq(levpt),iq(levnr),iq(arcp)
     2,           iq(nlwksp),iq(x),iq(nelecs),iq(spc),iwrite)
      nci = nwks
      call corlsi (a)
      a      = icori(nrows)
      b      = icori(nrows)
      s      = icori(nrows)
      arc    = icori(nrows4)
      nlwks  = icori(nrows)
      nuwks  = icori(nrows)
      puwk   = icori(nrows)
      indx   = icori(nwks)
      nwght  = icori(nrows4)
      call mcindx(iq(a),iq(bp),iq(sp),iq(arcp),iq(nlwksp),
     +            iq(b),iq(s),iq(arc),iq(nlwks),iq(nuwks),iq(puwk),
     +            iq(indx),
     +            nrowsp,iq(nwght),nrws4p,iq(levnr),iq(levpt),iq(iout),
     +            iq(orbsym),ibfcod,iwrite)
      symorb=nsym*norbs
      call mcloop(iq(b),iq(nlwks),iq(nuwks),iq(puwk),iq(nwght),
     +iq(a),iq(indx),
     +iq(arc),iq(levnr),iq(levpt),iq(orbsym),iwrite)
      call corlsi(ntype)
      iffmtp=1
      call accnt(' ',1)
      if (lto(4)) call testft(iq(1),iwrite)
      if (lto(2)) call editft(iq(1),iwrite)
      return
      end
      subroutine mcindx (nabca,nabcbo,nabcso,iarco,nlwkso,
     +  nabcb,nabcs,iarc,nlwks,nuwks,puwk,indx,nrowso,nwght,
     +  nrowsf,levnr,levpt,iout,isym,bfcode,iwrite)
      implicit real*8  (a-h,o-z)
      logical btest
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
      parameter (maxnr=1022)
c   at end of subroutine, nelec is returned in nabca
      character*4 jarcs
      integer bfcode
      common/entryp/icode(1)
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
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
      character*3 codes
      common /drtcoc/ codes(9)
      integer dela,delb,delele,virtul,occupd,valocc,rescor,resvir
     >       ,frozen,valvir,opensh,multi,speshl,multrf,valenc
     >       ,fzc,fzv,cor,vir,doc,uoc,alp,bet,spe
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     1,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     2,               valvir,opensh,multi,speshl,multrf,valenc
     3,   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
      dimension bfcode(nbf,nrefs)
      dimension levpt(nlevs),levnr(nlevs),iout(nbf),isym(nbf)
      integer puwk
      dimension iarco(nrows4), nlwkso(nrowso), iarc(4,nrows),nlwks(nrows
     1), nuwks(nrows), puwk(nrows), indx(nwks)
      dimension nabca(nrows), nabcb(nrows), nabcs(nrows),
     1nabcbo(nrowso), nabcso(nrowso)
      common /lsort / iwght(maxnr),itest(maxnr)
      dimension lwght(32), levir(32)
      dimension ishift(4), jarcs(4)
c     dimension ishfto(4)
      dimension nwght(nrowsf)
c      data ismt/0,1,8,9,64,65,72,73/
c  initialization and reduction of core space
      call accnt('mcindx',2)
      nindx = nwks
      do 10 i=1,nrows
10    nabcb(i)=nabcbo(i)
      do 20 i=1,nrows
20    nabcs(i)=nabcso(i)
      do 30 k=1,4
      ishift(k)=(k-1)*nrows
c     ishfto(k)=(k-1)*nrowso
30    continue
      do 40 i=1,nrows
      do 40 k=1,4
      iarc(k,i)=iarco((i-1)*4+k)
40    continue
      do 50 i=1,nrows
50    nlwks(i)=nlwkso(i)
      do 60 i=1,nindx
60    indx(i)=0
      do 70 i=1,nrows
      nuwks(i)=0
      puwk(i)=0
70    continue
c  generate number of upper walks array (nuwks)
c  generate primary upper walk array  (puwk)
      levm=nlevs
      nuwks(1)=1
      puwk(1)=1
80    lev=levm
      levm=lev-1
      if (levm.eq.0) go to 140
      nr=levnr(lev)
      nrm=levnr(levm)
      nptmx=levpt(levm)
      if (nr.gt.maxnr) call caserr('too many rows in one level')
      do 90 i=1,nr
90    iwght(i)=0
      do 130 k=1,4
      nptm=levpt(levm)
      do 100 i=1,nrm
      nptm=nptm+1
      itest(i)=0
      if (puwk(nptm).gt.0) itest(i)=-1
100   continue
      npt=levpt(lev)
      do 120 i=1,nr
      npt=npt+1
      jarc=iarc(k,npt)
      if (jarc.eq.0) go to 120
      ipt=jarc+nptmx
      nuwks(ipt)=nuwks(ipt)+nuwks(npt)
      if (itest(jarc).lt.0) go to 110
      itest(jarc)=1
      puwk(ipt)=puwk(npt)+iwght(i)
110   iwght(i)=iwght(i)+nlwks(ipt)
c      print*,'i,npt,ipt,jarc,itest(jarc),iwght(i),nlwks(ipt),puwk(ipt)'
c     >      , i,npt,ipt,jarc,itest(jarc),iwght(i),nlwks(ipt),puwk(ipt)
120   continue
130   continue
      go to 80
140   continue
c  write out arrays
c     nrows8=nrows*8
c  print out arrays if requested
      if (btest(iprint,1)) then
      optio(1)=0
      write(iwrite,330)
      lev=nlevs
150   continue
      nr=levnr(lev)
      npt=levpt(lev)
      levm=lev-1
      write(iwrite,340)
      write(iwrite,350)
      do 170 i=1,nr
      npt=npt+1
      do 160 k=1,4
      write (jarcs(k),'(i4)') iarc(k,npt)
      if (iarc(k,npt).eq.0) jarcs(k) = ' '
160   continue
      ia=nabca(npt)
      ib=nabcb(npt)
      ic=levm-ia-ib
      ism=nabcs(npt)
      write(iwrite,360) levm,i,ia,ib,ic,ism,jarcs,nlwks(npt),nuwks(npt),
     1puwk(npt)
170   continue
      write(iwrite,340)
      lev=lev-1
      if (lev.eq.0) go to 190
      icd=bfcode(icode(lev),1)
      do 180 i=1,nbf
      if (iout(i).eq.lev) iorb=i
180   continue
      ism=isym(lev)
c      ism=ismt(ism)
      write(iwrite,370) lev,iorb,ism,codes(icd)
      go to 150
190   continue
      end if
c  generate the indexing array  (indx)
c      print*,'generate indx'
      lev=1
      iwk=0
      levm=1
      ir=1
      lwght(1)=0
      iw=1
200   continue
c      print*,'200; lev,levm,lwght(levm),iw,ir ',lev,levm,lwght(levm),iw
c     >,ir
      lwght(lev)=lwght(levm)+iw
      if (lev.eq.nlevs) go to 240
      levir(lev)=ir
      levm=lev
      lev=levm+1
      levir(lev)=levnr(lev)+1
210   continue
      ir=levir(lev)
      nptx=levpt(lev)
      nptmx=levpt(levm)
      irm=levir(levm)
220   continue
      ir=ir-1
      if (ir.eq.0) go to 250
      npt=ir+nptx
      iw=0
      do 230 k=1,4
      jarc=iarc(k,npt)
      if (jarc.eq.0) go to 230
      if (irm.eq.jarc) go to 200
      ipt=jarc+nptmx
      iw=iw+nlwks(ipt)
230   continue
      go to 220
240   continue
      iwk=iwk+1
      iw=lwght(nlevs)
c      print*,'setting indx(',iw,') = ',iwk
      indx(iw)=iwk
250   continue
      lev=levm
      levm=lev-1
      if (levm.gt.0) go to 210
      do 280 k=1,4
      do 280 i=1,nrows
      iarpt=i+ishift(k)
      nwght(iarpt)=0
280   continue
      do 320 lev=2,nlevs
      levm=lev-1
      nr=levnr(lev)
      npt=levpt(lev)
      nptm=levpt(levm)
      do 310 i=1,nr
      npt=npt+1
      do 300 k=1,4
      iarpt=npt+ishift(k)
      j=iarc(k,npt)
      iarptx=iarpt+nrows
      if (j.gt.0) go to 290
      if (k.eq.4) go to 300
      nwght(iarptx)=nwght(iarpt)
      go to 300
290   continue
      iptm=j+nptm
      iarc(k,npt)=iptm
      if (k.eq.4) go to 300
      nwght(iarptx)=nwght(iarpt)+nlwks(iptm)
300   continue
310   continue
320   continue
      call accnt(' ',2)
      return
c      call dump
c
330   format(//25x,22('*')//25x,'the distinct row table'//25x,22('*')/)
340   format (1x)
350   format (8x,'i   j      a  b  c   sym      k0  k1  k2  k3      nlwk
     1s    nuwks     puwk')
360   format (5x,2i4,4x,3i3,3x,i2,5x,4a4,2x,3i9)
370   format (6h level,i3,9h  orbital,i3,6h  sym ,i2,4x,a4)
      end
      subroutine mcincs
      implicit real*8  (a-h,o-z)
      logical diagon
      integer orbtbf
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
      common/entryp/orbtbf(maxorb),nloop
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk
     1,juwk,itrak,isegt,diagon,nw1,nw2,nform1,nform2,ijklty,iflg11,iflg5
     2,kseg
      nloop = 0
      call mcinmx
      return
      end
      subroutine mcputo(index)
      implicit real*8  (a-h,o-z)
      logical diagon
      integer orbtbf
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
      common/entryp/orbtbf(maxorb),nloop
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk
     1,juwk,itrak,isegt,diagon,nw1,nw2,nform1,nform2,ijklty,iflg11,iflg5
     2,kseg
      dimension index (nwks),itrac(2,18)
      data itrac /1,0, 2,0, 3,0, 9,9, 9,9, 9,9, 3,1, 0,0, 1,2, 2,1,
     1            1,3, 3,2, 3,2, 10*9/
      if ((diagon.and.iuwk.ne.juwk).or.((.not.diagon).and.iuwk.eq.juwk))
     1 return
c      call accnt('mcputo  ',2)
      val1 = acf
      val2 = d
      nloop=nloop+1
      itr1=itrac(1,itrak)
      itr2=itrac(2,itrak)
      if (itr2.eq.0) val2=0.0d0
      if (itrak.eq.13 .or. itrak.eq.7) val2=1.0d0
      if (itr2+itr1.gt.5) call caserr('illegal track value')
      itype=18
      if (itrak.eq.8) itype=19
      val2=val1*val2
      do 400 i=1,nlwk
      ii=index(iuwk)
      jj=index(juwk)
      iuwk=iuwk+1
      juwk=juwk+1
      do 390 j=1,nuwk
      goto (10,120,150,210,260,350,370),ijklty
c
c...  ikjl  ijkl  iljk
c
10    if (itr2.ne.3) goto 20
      if (itr1.ne.1) call caserr('ijkl error type 1,3')
c.... type 3a,3b
      call outmx2 (val1,iorb,korb,jorb,lorb,ii,jj)
      call outmx2 (val2,jorb,korb,iorb,lorb,ii,jj)
      goto 380
20    goto (30,40,60),itr1
c.... type 1a
30    call outmx2 (val1,korb,iorb,jorb,lorb,jj,ii)
      goto 80
40    if (isegt.le.20) goto 50
c..  raising op closure
c.... type 2b
      call outmx2 (val1,iorb,jorb,korb,lorb,ii,jj)
      goto 80
c..  lowring op closure
c.... type 1b
50    call outmx2 (val1,jorb,iorb,korb,lorb,jj,ii)
      goto 80
60    if (iflg11.ne.1) goto 70
c.... type 3b
      call outmx2 (val1,jorb,korb,iorb,lorb,ii,jj)
      goto 80
c.... type 2a
70    call outmx2 (val1,korb,jorb,iorb,lorb,ii,jj)
80    if (itr2.eq.0) goto 380
      if (itr2.ne.2) call caserr('error ijkl tr2 #2')
      goto (90,100,110),itr1
c.... type 1b
90    call outmx2 (val2,jorb,iorb,korb,lorb,jj,ii)
      goto 380
100   continue
c.... type 2b
110   call outmx2 (val2,iorb,jorb,korb,lorb,ii,jj)
      goto 380
c
c....  ijjl  iljj
120   continue
      goto (130,140),itr1
c.... type 4
130   call outmx2 (val1,jorb,iorb,jorb,lorb,jj,ii)
      goto 380
c.... type 5b
140   continue
c      write (iwrite,*) iorb,jorb,lorb,ii,jj,itr1,itr2,val1,val2,kseg
c     * ,isegt
      call outmx2 (val1,jorb,jorb,iorb,lorb,ii,jj)
c.... type 5a
c      if (itr2.eq.0) write (iwrite,*) val1,val2,itr1,itr2
      if (itr2.ne.0) call outmx2 (val2,jorb,lorb,iorb,jorb,ii,jj)
      goto 380
c
c....  ikil  iikl    loops 6 7a 7b
150   if (itr2.ne.0) goto 200
      goto (160,190),itr1
160   if (isegt.ge.21) goto 170
      call outmx2 (val1,iorb,korb,lorb,iorb,ii,jj)
      goto 380
170   if (iflg5.ne.0) goto 180
      call outmx2 (val1,iorb,lorb,korb,iorb,ii,jj)
      goto 380
180   call outmx2 (val1,iorb,korb,iorb,lorb,ii,jj)
      goto 380
190   call outmx2 (val1,iorb,iorb,korb,lorb,ii,jj)
      goto 380
200   call outmx2 (val1,iorb,lorb,korb,iorb,ii,jj)
      call outmx2 (val2,iorb,iorb,korb,lorb,ii,jj)
      goto 380
c
c...   iljl  ijll
210   continue
      if (kseg.ne.138) goto 220
      call outmx2 (val1,iorb,lorb,jorb,lorb,ii,jj)
      goto 380
220   if (itr2.ne.0) goto 250
      goto (230,240),itr1
c.... type 9
230   call outmx2 (val1,lorb,iorb,jorb,lorb,jj,ii)
      goto 380
240   call outmx2 (val1,lorb,lorb,iorb,jorb,ii,jj)
      goto 380
250   call outmx2 (val1,lorb,jorb,iorb,lorb,ii,jj)
      call outmx2 (val2,lorb,lorb,iorb,jorb,ii,jj)
      goto 380
c
c....   iiil  illl  il
260   if (itype.eq.19) goto 340
      goto (270,280,290),itr1
270   call caserr('improper iiil element')
280   call caserr('improper illl element')
290   call outmx1 (val1,iorb,lorb,ii,jj)
c      if (ii.eq.492.or.jj.eq.492)
c     >print*,'outmx1 53 kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb'
c     >,kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb
      if (itr2.eq.0) goto 380
      goto (310,320,330),itr2
310   call outmx2 (val1,iorb,lorb,iorb,iorb,ii,jj)
      goto 380
320   call outmx2 (val1,lorb,lorb,iorb,lorb,ii,jj)
      goto 380
330   call caserr('improper <i/h/l> element')
340    call outmx2 (val1,iorb,lorb,iorb,iorb,ii,jj)
      call outmx2 (val1,lorb,lorb,iorb,lorb,ii,jj)
      call outmx1 (val1,iorb,lorb,ii,jj)
c      if (ii.eq.492.or.jj.eq.492)
c     >print*,'outmx1 58 kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb'
c     >,kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb
      goto 380
c
c...   ilil  iill
350   continue
      if (kseg.ne.129) goto 360
c.... type 12
      call outmx2 (val1,iorb,lorb,iorb,lorb,ii,jj)
      goto 380
c.... type 13a/13b
360   call outmx2 (val1,lorb,iorb,iorb,lorb,ii,jj)
      call outmx2 (val2,iorb,iorb,lorb,lorb,ii,jj)
      goto 380
c
c...  iiii  ii
370     call outmx1 (val1,iorb,iorb,ii,jj)
c      if (ii.eq.492.or.jj.eq.492)
c     >print*,'outmx1 70 kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb'
c     >,kseg,val1,val2,itrak,ii,jj,iorb,jorb,korb,lorb
      if (itr2.ne.0) call outmx2 (val2,iorb,jorb,korb,lorb,ii,jj)
       goto 380
c
380   ii=ii+1
390   jj=jj+1
400    continue
c      call accnt('mcloop  ',2)
      return
      end
      subroutine mcfino(iwrite)
      implicit real*8  (a-h,o-z)
      logical diagon
      integer orbtbf
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
      common/entryp/orbtbf(maxorb),nloop
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk
     1,juwk,itrak,isegt,diagon,nw1,nw2,nform1,nform2,ijklty,iflg11,iflg5
     2,kseg
      character *3 iword,ion,ioff
      data ion,ioff/' on','off'/
      iword=ioff
      if (diagon) iword=ion
      write (iwrite,410) nloop,iword
410   format(i8,a5,' diagonal loops generated')
      nloop=0
      call mcfinx(iwrite)
      return
      end
      subroutine mcinmx
      implicit real*8  (a-h,o-z)
      logical diagon
      integer orbtbf
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
      parameter (mxorb1=maxorb+1)
      common/entryp/orbtbf(mxorb1),i0
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
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common /couple/ surd(511)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk,juwk
     1 ,itrak,isegt,diagon,nw1,nw2,nform1,nform2
      i0=mcsrd1(1.0d0,-1.0d0)
      i0=mcsurd(0.0d0)
      nform1=0
      nform2=0
      nw1=0
      nw2=0
      call search (iblk6,num6)
      call search (iblk4,num4)
      return
      end
      subroutine outmx1 (val,i,j,ii,jj)
      implicit real*8  (a-h,o-z)
      logical diagon
      integer out1,out2,packft
      integer orbtbf
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
      parameter (mxorb1=maxorb+1)
      common/entryp/orbtbf(mxorb1),i0
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
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common /couple/ surd(511)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk,juwk
     1 ,itrak,isegt,diagon,nw1,nw2,nform1,nform2
      common /lsort / space(511,3),out1(511),out2(511)
      ival=mcsurd(val)
      if (ival.eq.i0) return
      nw1=nw1+1
      ival=mcsurd(val)
      out1(nw1)=packft(ii,jj,ival,orbtbf(i),orbtbf(j),0,0)
      if (nw1.lt.511) return
      call put (out1,nw1,num6)
      nw1=0
      nform1=nform1+511
      return
      end
      subroutine outmx2 (val,i,j,k,l,ii,jj)
      implicit real*8  (a-h,o-z)
      logical diagon
      integer out1,out2,packft
      integer orbtbf
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
      parameter (mxorb1=maxorb+1)
      common/entryp/orbtbf(mxorb1),i0
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
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common /couple/ surd(511)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk,juwk
     1 ,itrak,isegt,diagon,nw1,nw2,nform1,nform2
      common /lsort / space(511,3),out1(511),out2(511)
      ival=mcsurd(val)
      if (ival.eq.i0) return
      nw2=nw2+1
      ival=mcsurd(val)
      out2(nw2)=
     +  packft(ii,jj,ival,orbtbf(i),orbtbf(j),orbtbf(k),orbtbf(l))
      if (nw2.lt.511) return
      call put (out2,nw2,num4)
      nw2=0
      nform2=nform2+511
      return
      end
      subroutine mcfinx(iwrite)
      implicit real*8  (a-h,o-z)
      logical diagon
      integer  out1,out2
      integer orbtbf
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
      common/entryp/orbtbf(1)
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
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
      common /couple/ surd(511)
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk,juwk
     1 ,itrak,isegt,diagon,nw1,nw2,nform1,nform2
      common /lsort / space(511,3),out1(511),out2(511)
      character*3 iword
      nform1=nform1+nw1
      nform2=nform2+nw2
      nblk1=lensec(nform1)
      nblk2=lensec(nform2)
      if (nform1.le.0) nblk1=0
      if (nform2.le.0) nblk2=0
      if (nw1.gt.0) call put (out1,nw1,num6)
      if (nw2.gt.0) call put (out2,nw2,num4)
      nw1=0
      nw2=0
      if (diagon) goto 10
      nblko1=nblk1
      nblko2=nblk2
      nblk1=nblk1+nblkd1
      nblk2=nblk2+nblkd2
      call search (iblk6,num6)
      call search (iblk4,num4)
      call wrt3 (surd,511,iblkft,numft)
      do 11 iblk=1,nblk1
      call find (num6)
      call get (out1,nw)
11    call put (out1,nw,numft)
      call put (out1,0,numft)
      do 12 iblk=1,nblk2
      call find (num4)
      call get (out2,nw)
12    call put (out2,nw,numft)
      call put (out2,0,numft)
      call clredx
      iword='off'
      goto 20
10    nblkd1=nblk1
      nblkd2=nblk2
      iword='on'
20    write (iwrite,30) nform1,iword,'one'
      write (iwrite,30) nform2,iword,'two'
30    format(i8,a5,' diagonal',a4,' particle formulae')
      nform1=0
      nform2=0
      nblk=iposun(numft)-iblkft
      if (.not.diagon) write (iwrite,55) nblk
55    format(' length of formula tape =',i7,' blocks')
      return
      end
      function mcsrd1(val1st,val2nd)
      implicit real*8  (a-h,o-z)
      parameter (small=1.0d-6)
      common/couple/surd(511),nsurd
      common/lsort /link(511,2),ibal(511)
      val1=val1st
      val2=val2nd
      if(val1.lt.val2)goto 80
      val1=val2nd
      val2=val1st
80    link(1,1)=0
      link(1,2)=2
      link(2,1)=0
      link(2,2)=0
      ibal(1)=2
      ibal(2)=0
      link(511,1)=2
      link(511,2)=1
      surd(1)=val1
      surd(2)=val2
      nsurd=2
      mcsrd1=0
      return
      end
      function mcsurd(val)
      implicit real*8  (a-h,o-z)
      parameter (small=1.0d-6)
      common/couple/surd(511),nsurd
      common/lsort /link(511,2),ibal(511)
c     locate the index of coupling coefficient val by searching a balanc
c     binary tree. if unsuccessful, tree is updated, and re-balanced if
c     necessary. (see knuth, vol3, p451).                 pjk/15-iii-82
c
      linki=link(511,2)
      is=linki
      it=511
10    mcsurd=linki
      test=val-surd(mcsurd)
      if( dabs(test).lt.small)return
      lr=1
      if(test.gt.0.0d0)lr=2
      linki=link(mcsurd,lr)
      if(linki.eq.0)goto 20
      if(ibal(linki).eq.0)goto 10
      it=mcsurd
      is=linki
      goto 10
20    nsurd=nsurd+1
      if(nsurd.gt.511)call caserr('table of coefficients too big.')
      link(mcsurd,lr)=nsurd
      surd(nsurd)=val
      link(nsurd,1)=0
      link(nsurd,2)=0
      ibal(nsurd)=0
c.....adjust balance factors
      lr=1
      if(val.gt.surd(is))lr=2
      ir=link(is,lr)
      mcsurd=ir
30    if(mcsurd.eq.nsurd)goto 40
      lr=1
      if(val.gt.surd(mcsurd))lr=2
      ibal(mcsurd)=lr
      mcsurd=link(mcsurd,lr)
      goto 30
40    lr=1
      if(val.gt.surd(is))lr=2
      ibals=ibal(is)
      if(ibals.eq.lr)goto 50
      ibal(is)=0
      if(ibals.ne.0)return
      ibal(is)=lr
      link(511,1)=link(511,1)+1
      return
c.....rebalance the tree
50    if(ibal(ir).ne.lr)goto 60
c.... single rotation
      mcsurd=ir
      link(is,lr)=link(ir,3-lr)
      link(ir,3-lr)=is
      ibal(is)=0
      ibal(ir)=0
      goto 70
c.....double rotation
60    mcsurd=link(ir,3-lr)
      link(ir,3-lr)=link(mcsurd,lr)
      link(mcsurd,lr)=ir
      link(is,lr)=link(mcsurd,3-lr)
      link(mcsurd,3-lr)=is
      ibalss=ibal(mcsurd)
      ibal(is)=0
      ibal(ir)=0
      if(ibalss.eq.lr)ibal(is)=3-lr
      if(ibalss.eq.3-lr)ibal(ir)=lr
      ibal(mcsurd)=0
c.....tidy up
70    lr=1
      if(is.eq.link(it,2))lr=2
      link(it,lr)=mcsurd
      mcsurd=nsurd
      return
      end
      subroutine mcloop(nabc,nlwks,nuwks,puwk,iwght,nelec,
     +                  index,iarc,levnr,levpt,isym,iwrite)
      implicit real*8  (a-h,o-z)
      integer orbtbf
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
      common/entryp/orbtbf(1)
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      logical diagon,intrac
      common/multil/acf,ccf,d,iorb,jorb,korb,lorb,nuwk,nlwk,iuwk,juwk
     > ,itrak,isegt,diagon,nw1,nw2,nform1,nform2,ijklty,iflg11,iflg5
     > ,kseg
      dimension isegm(32),jsegm(32),imain(32),isub(32)
      dimension iuwkmn(32),iuwksb(32),itrack(32),ishift(4),isym(32)
      dimension acoef(32),bcoef(32),levnr(32),levpt(32)
      dimension iarc (4,nrows), index(nwks)
c      dimension ismoff(8)
      dimension lmin(32)
      dimension nabc(nrows),nlwks(nrows),nuwks(nrows),
     *puwk(nrows),iwght(nrows4)
      dimension nelec(nrows),nlcsmn(22)
       dimension lcond(8)
       dimension jsegnr(22),jsegpt(22),iarcmn(228),iarcsb(228)
       dimension itrk(228),jcond(228),kcond(228),nxtseg(228)
       dimension coeffs(20,21),cfs(420)
       dimension jsegpx(3)
      dimension mults(8),lkupsm(64)
      equivalence (coeffs(1,1),cfs(1))
       integer puwk
       data mults/0,8,16,24,32,40,48,56/
       data lkupsm/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,
     a 2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,
     b 8,7,6,5,4,3,2,1/
       data jsegnr/16,34,52,63,75,92,102,118,128,137,148,155,162,172,
     a 179,186,193,200,207,214,221,228/
       data jcond/12*1,4*-1,13*1,5*-1,13*1,5*-1,176*0/
       data kcond/7*1,9*0,5*1,13*0,5*1,13*0,7*1,4*0,7*1,5*0,11*1,6*0,
     a 5*1,5*0,10*1,6*0,5*1,5*0,5*1,4*0,5*1,6*0,1,1,5*0,1,1,5*0,4*1,
     b 62*0/
       data itrk/10,2,2,3,7,3,7,1,3*9,5*1,0,2*10,2,10,4*11,2*9,3,3,6*0,
     a  10,10,2,10,4*11,2*9,2*3,27*0,1,3*0,4*12,8*0,2*3,1,2*0,2*1,5*0,
     b  2*1,4*0,4*1,6*0,1,0,0,1,1,14*0,1,35*0,8,6*0,8,6*0,13,6*0,13,33*0
     c /
       data nxtseg/3*0,17,15,18,16,10,5,5,4,7,2,2,3,3,20,3*21,22,11,11,
     a 12,12,6,6,7,7,4*2,3,19,3*22,21,11,11,13,13,6,6,9,9,4*3,2,3*0,
     b 21,21,22,22,4*4,3*0,21,21,22,22,4*5,7,3*0,21,21,22,22,19,19,20,
     c 20,4*6,7,9,0,21,21,20,20,4*7,8,0,0,21,21,22,22,19,19,20,20,4*8,
     d 7,9,0,22,22,19,19,4*9,8,0,21,21,22,22,4*10,0,21,21,22,22,4*11,
     e 12,13,21,21,3*12,12,14,22,22,4*13,14,21,21,22,22,4*14,12,13,
     f 0,0,4*15,16,0,0,4*16,15,0,0,4*17     ,18,0,0,4*18,17,0,0,4*19,
     g 20,0,0,20,3*20,19,0,0,4*21,22,0,0,4*22,21/
      data iarcmn/4,3,2,3,4,2,4,4,2,3,4,3,3,4,2,4,1,3,4,3*2,4,3,4,1,2,1,
     a 3,1,2,3,4,2,1,2,4,3*3,4,2,4,1,3,1,2,1,2,3,4,3,2,3,4,3,4,2,4,
     b 1,2,3,4,2,3,4,3,4,2,4,1,2,3,4,3,2,3,4,3,4,2,4,1,2,1,3,1,2,3,4,
     c 3,2,2,2,4,1,2,1,2,3,4,2,2,3,3,4,2,4,1,2,1,3,1,2,3,4,3,2,3,3,4,
     d 1,3,1,2,3,4,3,1,1,3,1,2,1,2,3,4,1,1,3,1,2,1,2,3,4,3,2,1,2,1,2,
     e 3,4,2,1,3,1,2,3,4,3,1,3,1,2,1,2,3,4,3,2,1,2,1,2,3,4,2,1,3,1,2,
     f 3,4,3,1,2,1,2,3,4,2,1,3,1,2,3,4,3,3,4,1,2,3,4,3,2,4,1,2,3,4,2,
     g 1,2,1,2,3,4,2,1,3,1,2,3,4,3/
      data iarcsb/4,3,2,1,2,1,3,1,2,3,4,2,1,2,1,3,4,3,4,2,3,1,3,1,2,3,4,
     a 2,4,1,2,3,4,3,4,2,4,3,2,1,2,1,3,2,4,3,4,1,2,3,4,2,2,3,4,1,2,1,3,
     b 1,2,3,4,2,3,4,1,2,1,3,1,2,3,4,2,2,3,4,1,2,1,3,3,4,2,4,1,2,3,4,
     c 2,3,3,1,3,3,4,1,2,3,4,3,2,3,1,2,1,3,3,4,2,4,1,2,3,4,2,3,2,1,2,
     d 2,4,1,2,3,4,2,4,2,4,3,4,1,2,3,4,4,2,4,3,4,1,2,3,4,2,3,3,4,1,2,
     e 3,4,3,2,4,1,2,3,4,2,2,4,3,4,1,2,3,4,2,3,3,4,1,2,3,4,3,2,4,1,2,
     f 3,4,2,3,4,1,2,3,4,3,2,4,1,2,3,4,2,1,2,1,2,3,4,2,1,3,1,2,3,4,3,
     g 3,4,1,2,3,4,3,2,4,1,2,3,4,2/
       data jsegpx/12,29,47/
c     data qz/50,1,1,1,3,1,4,1,44,45,51,1,1,3,1,4,5,52,53,6,54,40,41,1,7
c    a,46,47,1,8,1,6,2,2,9,10,52,53,11,55,42,43,1,12,48,49,1,13,1,2,11,
c    b 2,36,77,77,78,77,79,77,80,1,1,1,1,87,68,82,68,69,67,70,71,75,
c    c 76,71,83,87,68,82,68,69,67,70,68,69,67,70,71,75,76,71,83,83,
c    d 6,6,74,6,8,1,16,16,1,17,18,19,19,20,18,21,19,20,18,21,1,22,23,
c    e 1,24,24,11,11,81,11,13,1,27,27,1,28,1,3,2,4,2,1,2,2,1,29,63,64,
c    f 65,66,71,72,73,71,84,85,56,57,1,30,1,1,86,56,58,1,1,31,1,32,
c    g 59,60,61,62,1,33,30,1,34,35,1,5,1,6,2,2,9,1,10,1,2,11,2,36,
c    h 1,5,1,6,2,2,9,1,10,1,2,11,2,36,1,4,1,37,2,2,38,1,3,1,2,37,2,39,
c    i 1,5,1,6,2,2,9,1,10,1,2,11,2,36/
       data nlcsmn/6*1,0,1,14*0/
      call accnt('mcloop  ',2)
      levfrm = 0
      diagon = .true.
      call mcincs
c     lvfrm=levfrm
      crite = 0.00001d0
      root2 =  dsqrt(2.0d0)
      rootn2 = -root2
      toor2 = 1.0d0/ root2
      toorn2 = -toor2
      jsegpt(1)=0
      do 130 i=1,21
 130  jsegpt(i+1)=jsegnr(i)
      do 135 i=1,2
      do 135 j=1,21
 135  coeffs(i,j)=0.0d0
      do 140 i=3,20
      a =  dble(i-2)
      coeffs(i,1) =  dsqrt(a/(a+1.0d0))
      coeffs(i,2) = -coeffs(i,1)
      coeffs(i,3) = coeffs(i,1)/ dsqrt(2.0d0)
      coeffs(i,4) = -coeffs(i,3)
      coeffs(i,5) =  dsqrt((a+1.0d0)/a)
      coeffs(i,6) = -coeffs(i,5)
      coeffs(i,7) = coeffs(i,5)/ dsqrt(2.0d0)
      coeffs(i,8) = -coeffs(i,7)
      coeffs(i,9) =  dsqrt((a+2.0d0)/(a*2.0d0))
      coeffs(i,10) = -coeffs(i,9)
      coeffs(i,11) =  dsqrt(a/(2.0d0*(a+2.0d0)))
      coeffs(i,12) = -coeffs(i,11)
      coeffs(i,13) =  dsqrt(2.0d0/(a*(a+1.0d0)))
      coeffs(i,14) =  dsqrt(a*(a+2.0d0))/(a+1.0d0)
      coeffs(i,15) = - dsqrt(a*(a+2.0d0))/(a+1.0d0)
      coeffs(i,16) =  dsqrt((a-1.0d0)*(a+2.0d0)/(a*(a+1.0d0)))
      coeffs(i,17) = -coeffs(i,16)
      coeffs(i,18)=- dsqrt(2.0d0/(a*(a+2.0d0)))
      coeffs(i,19) = 1.0d0/a
      coeffs(i,20) = -1.0d0/a
      coeffs(i,21) = - dsqrt(2.0d0)/a
140   continue
      do 150 k=1,4
      ishift(k)=(k-1)*nrows
150   continue
      do 155 i=1,nrows
      nelec(i)=nelec(i)+nelec(i)+nabc(i)
155   continue
      do 160 i=1,nsym
c     ismoff(i)=(i-1)*norbs
 160  lcond(i)=0
      i=isym(1)
      lcond(i)=1
      lcond(1)=1
      nsm=0
      do 170 iorbi=2,norbs
      do 165 i=1,nsym
      if(lcond(i).eq.0) go to 165
      ism=i
      lkup=mults(ism)+isym(iorbi)
      j=lkupsm(lkup)
      if(lcond(j).gt.0) go to 165
      lcond(j)=iorbi
      nsm=nsm+1
      if(nsm.eq.nsym) go to 175
 165  continue
 170  continue
 175  continue
      do 180 i=1,nsym
      if(lcond(i).eq.0) lcond(i)=norbs+1
 180  continue
10000 continue
      intrac=.false.
      iflg11 = 0
      iflg5  = 0
      i = norbs+1
190   continue
      i = i-1
      levi=i+1
      jmax=i
      jmin=1
      iorb=i
c      lflag = lto(3) .and. iorb.eq.4
c      if (lflag) print*,'setting iorb = ',iorb
      jsm=isym(i)
      lev=levi
      levm=lev-1
      nr=levnr(lev)
      npt=levpt(lev)
      do 480 irow=1,nr
      npt=npt+1
      isegm(lev)=1
      iseg=1
      imn=npt
      isb=npt
      kseg=0
      ksegmx=jsegnr(iseg)
      lmin(lev)=lcond(jsm)
      iuwkmn(lev)=puwk(npt)
      iuwksb(lev)=puwk(npt)
      imain(lev)=npt
      isub(lev)=npt
      nuwk=nuwks(npt)
      acoef(lev)=1.0d0
c  test next segment of group
200   kseg=kseg+1
      if(kseg.gt.ksegmx) go to 440
c      if (lflag) print*,'trying kseg = ',kseg
      kmn=iarcmn(kseg)
      iarpt=imn+ishift(kmn)
      kmn=iarc(kmn,imn)
      if(kmn.eq.0) go to 200
      ksb=iarcsb(kseg)
      jarpt=isb+ishift(ksb)
      ksb=iarc(ksb,isb)
      if(ksb.eq.0) go to 200
      if (diagon .and. kmn.ne.ksb) goto 200
      jsegm(lev)=kseg
      iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
      iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
      lmin(levm)=lmin(lev)
      if(jcond(kseg))220,240,230
 220  continue
      if(levm.le.jmin)goto 440
      goto 240
 230  continue
      if(levm.gt.jmax)goto 420
      jorb=levm
c      if (lflag) print*,'setting jorb = ',jorb
      ieqj=0
      if (iorb.eq.jorb) ieqj=2
c      ij=iad+levm
      lkup=mults(jsm)+isym(levm)
      ksm=lkupsm(lkup)
      lmin(levm)=lcond(ksm)
 240  continue
      if(kcond(kseg).eq.0)goto 260
      korb=levm
c      if (lflag) print*,'setting korb = ',korb
      jeqk=0
      if (jorb.eq.korb) jeqk=1
      lkup=mults(ksm)+isym(levm)
      lsm=lkupsm(lkup)
      lmin(levm)=lcond(lsm)
 260  continue
      if(itrk(kseg))280,280,270
270   itrack(levm)=itrk(kseg)
      goto 290
280   itrack(levm)=itrack(lev)
290   continue
c     iloc=qz(kseg)
c      if (lflag) print*,'entering big goto at level ',levm,
c     >' with kseg = ',kseg
c     ijump(iloc)=ijump(iloc)+1
c       goto(50,1,1,1,3,1,4,1,44,45,51,1,1,3,1,4,5,52,53,6,54,40,41,1,7,
       goto(50,1,1,1,3,1,4,89,44,45,51,1,1,3,1,4,5,52,53,6,54,40,41,1,7,
     a 46,47,1,8,1,6,2,2,9,10,52,53,11,55,42,43,1,12,48,49,1,13,1,2,11,
     b 2,36,77,77,78,77,79,77,80,1,1,1,1,87,68,82,68,69,67,70,71,75,
     c 76,71,83,87,68,82,68,69,67,70,68,69,67,70,71,75,76,71,83,83,
     d 6,6,74,6,8,1,16,16,1,17,18,19,19,20,18,21,19,20,18,21,1,22,23,
     e 1,24,24,11,11,81,11,13,1,27,27,1,28,1,3,2,4,2,1,2,2,1,29,63,64,
     f 65,66,71,72,73,71,84,85,56,57,1,30,1,1,86,56,58,1,1,31,1,32,
     g 59,60,61,62,1,33,88,1,34,35,1,5,1,6,2,2,9,1,10,1,2,11,2,36,
     h 1,5,1,6,2,2,9,1,10,1,2,11,2,36,1,4,1,37,2,2,38,1,3,1,2,37,2,39,
     i 1,5,1,6,2,2,9,1,10,1,2,11,2,36),kseg
 1    acoef(levm)=acoef(lev)
      goto 120
  2   acoef(levm)=-acoef(lev)
      goto 120
  3   ia = nabc(imn) + 2
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  4   ia = nabc(imn) + 83
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  5   ia = nabc(imn) + 82
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  6   ia = nabc(imn) + 261
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  7   ia = nabc(imn) + 1
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  8   ia = nabc(imn) + 102
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
  9   ia = nabc(imn) + 362
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 10   ia = nabc(imn) + 3
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 11   ia = nabc(imn) + 263
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 12   ia = nabc(imn) + 84
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 13   ia = nabc(imn) + 23
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 16   ia = nabc(imn) + 281
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 17   ia = nabc(imn) + 402
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 18   ia = nabc(imn) + 162
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 19   ia = nabc(imn) + 222
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 20   ia = nabc(imn) + 143
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 21   ia = nabc(imn) + 42
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 22   ia = nabc(imn) + 302
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 23   ia = nabc(imn) + 303
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 24   ia = nabc(imn) + 342
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 27   ia = nabc(imn) + 283
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 28   ia=nabc(imn)+404
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 29   acoef(levm) = acoef(lev) * root2
      go to 120
 30   ia = nabc(imn) + 301
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 31   ia = nabc(imn) + 304
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 32   ia = nabc(imn) + 244
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 33   ia = nabc(imn) + 322
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 34   ia = nabc(imn) + 243
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 35   ia = nabc(imn) + 242
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 36   ia = nabc(imn) + 384
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 37   ia = nabc(imn) + 262
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 38   ia = nabc(imn) + 363
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 39   ia = nabc(imn) + 383
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 86   ia = nabc(imn) + 241
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 40   ia = nabc(imn) + 122
      ib=ia-61
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 41   ib = nabc(imn) + 162
      acoef(levm) = acoef(lev) * toorn2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 42   ia = nabc(imn) + 43
      ib = ia + 81
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 43   ib = nabc(imn) + 222
      acoef(levm) = acoef(lev) * toorn2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 44   ib=nabc(imn)+221
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 45   ib = nabc(imn) + 163
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 46   ib = nabc(imn) + 162
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 47   ia = nabc(imn) + 122
      ib = ia - 81
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 48   ib = nabc(imn) + 222
      acoef(levm) = acoef(lev) * toor2
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 49   ia = nabc(imn) + 43
      ib = ia + 101
      acoef(levm) = acoef(lev) * cfs(ia)
      bcoef(levm) = acoef(lev) * cfs(ib)
      go to 120
 50   acoef(levm) = acoef(lev) + acoef(lev)
      d=0.5d0
      go to 120
 51   acoef(levm)=acoef(lev)*root2
      go to 120
 52   acoef(levm) = -acoef(lev)
      d= -1.0d0
      go to 120
 53   acoef(levm) = -acoef(lev) - acoef(lev)
      d = -0.5d0
      go to 120
 54   ia=nabc(imn)+362
      d=1.0d0/cfs(ia)
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 55   ia = nabc(imn) + 384
      d=1.0d0/cfs(ia)
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 56   acoef(levm) = acoef(lev)
      d = -1.0d0
      go to 120
 57   ia = nabc(imn) + 82
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 58   ia = nabc(imn) + 3
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 59   ia = nabc(imn) + 123
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 60   ia = nabc(imn) + 222
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 61   ia = nabc(imn) + 62
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 62   ia = nabc(imn) + 162
      acoef(levm) = acoef(lev) * cfs(ia)
      d=-1.0d0
      go to 120
 63   ia = nabc(imn) + 42
      ib = ia + 81
      acof = acoef(lev) * cfs(ia)
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 64   ib = nabc(imn) + 222
      acof = acoef(lev) * toorn2
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 65   ia = nabc(imn) + 123
      ib = ia - 61
      acof = acoef(lev) * cfs(ia)
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 66   ib = nabc(imn) + 162
      acof = acoef(lev) * toorn2
      bcof = bcoef(lev) * cfs(ib)
      d = acof + bcof
      if( dabs(d).lt.crite) go to 110
      acoef(levm) = d
      d = (acof - bcof) / d
      go to 120
 67   ib = nabc(imn) + 162
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 68   ib = nabc(imn) + 222
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 69   ia = nabc(imn) + 62
      ib = ia + 81
      dx=acoef(lev)*cfs(ia)
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 70   ia = nabc(imn) + 143
      ib = ia - 101
      dx=acoef(lev)*cfs(ia)
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm) = d
      d=-(dx+dx)/d
      go to 120
 87   ib = nabc(imn) + 162
      dx=acoef(lev)*toorn2
      d=dx+bcoef(lev)*cfs(ib)
      if( dabs(d).lt.crite) go to 111
      acoef(levm)=d
      d=-(dx+dx)/d
      go to 120
 71   acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev)
      go to 120
 72   ib = nabc(imn) + 322
      acoef(levm) = -acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 73   ib = nabc(imn) + 323
      acoef(levm) = -acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 74   ia=nabc(imn)+21
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 75   ib = nabc(imn) + 302
      acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 76   ib = nabc(imn) + 303
      acoef(levm) = acoef(lev)
      bcoef(levm) = bcoef(lev) * cfs(ib)
      go to 120
 77   acoef(levm)=acoef(lev)*toorn2
      d=-2.0d0
      go to 120
 78   acoef(levm)=acoef(lev)*rootn2
      d=-2.0d0
      go to 120
 79   ia=nabc(imn)+62
      acoef(levm)=acoef(lev)*cfs(ia)
      d=-2.0d0
      go to 120
 80   ia=nabc(imn)+143
      acoef(levm)=acoef(lev)*cfs(ia)
      d=-2.0d0
      go to 120
 81   ia=nabc(imn)+104
      acoef(levm)=acoef(lev)*cfs(ia)
      go to 120
 82   acoef(levm) = acoef(lev) * rootn2
      d=-2.0d0
      go to 120
 83   ia = nabc(imn) + 342
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 84   ia = nabc(imn) + 243
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 85   ia = nabc(imn) + 242
      acoef(levm) = bcoef(lev) * cfs(ia)
      go to 120
 88   ia = nabc(imn) + 323
      acoef(levm) = acoef(lev) * cfs(ia)
      go to 120
 89   acoef(levm) = acoef(lev)
      iflg5=1
      goto 120
110   itrack(levm)=3
      acoef(levm)=acof-bcof
      iflg11 = 1
      go to 120
111   itrack(levm) = 2
      acoef(levm)=-(dx+dx)
120   continue
c      if (lflag) print*,'acoef,bcoef,d now ',acoef(levm),bcoef(levm),d
      if(nxtseg(kseg).gt.0) go to 400
      if(isym(levm).ne.lsm) go to 200
      lorb=levm
c      if (lflag) print*,'setting lorb = ',lorb
      ijklty=1+ieqj+jeqk
      if (lorb.eq.korb) ijklty=ijklty+3
      if (iorb.eq.korb .and. lorb.lt.korb) ijklty=ijklty+1
c      if (lflag) print*,'lorb selected; kmn,ksb = ',kmn,ksb
      if(kmn-ksb) 300,380,300
300   levl=levm
      ksegmx=4
310   lev=levm
      levm=lev-1
      if(levm.gt.0) go to 315
      write(iwrite,313)
313   format(' problems with partial space')
      call caserr('problems with partial space')
 315  continue
      kseg=0
      imain(lev)=kmn
      imn=kmn
      isub(lev)=ksb
      isb=ksb
320   kseg=kseg+1
      if(kseg.gt.ksegmx)goto 360
      iarpt=imn+ishift(kseg)
      kmn=iarc(kseg,imn)
      if(kmn.le.0)goto 320
      jarpt=isb+ishift(kseg)
      ksb=iarc(kseg,isb)
      if(ksb.le.0)goto 320
      jsegm(lev)=kseg
      iuwkmn(levm)=iuwkmn(lev)+iwght(iarpt)
      iuwksb(levm)=iuwksb(lev)+iwght(jarpt)
      if(kmn-ksb) 310,340,310
340   nlwk = nlwks(kmn)
      iuwk = iuwkmn(levm)
      juwk = iuwksb(levm)
      itrak = itrack(levl)
      acf = acoef(levl)
      isegt = isegm(lev)
      intrac=.true.
      call mcputo(index)
      go to 320
360   if(lev.eq.levl) go to 440
      levm=lev
      lev=levm+1
      imn=imain(lev)
      isb=isub(lev)
      kseg=jsegm(lev)
      go to 320
380   nlwk = nlwks(kmn)
      iuwk = iuwkmn(levm)
      juwk = iuwksb(levm)
      itrak = itrack(levm)
      acf = acoef(levm)
c..   save penult seg value
      isegt=isegm(lev)
c      if (lflag) print*,'calling mcputo'
      call mcputo(index)
      go to 200
400   continue
c      if(levm.le.lmin(levm))goto 200
      iseg=nxtseg(kseg)
      if(nlcsmn(iseg).gt.0.and.nelec(kmn).eq.0)goto 200
c      if(lvfrm.eq.levm.and.iseg.lt.4)goto 460
      lev=levm
      levm=lev-1
      isegm(lev)=iseg
      kseg=jsegpt(iseg)
      imn=kmn
      imain(lev)=kmn
      isb=ksb
      isub(lev)=ksb
      ksegmx=jsegnr(iseg)
      goto 200
420   continue
      kseg=jsegpx(iseg)
      goto 200
440   continue
      if(lev.eq.levi)goto 480
      levm=lev
      lev=levm+1
      if (lev.eq.levi) iflg5=0
      iseg=isegm(lev)
      if (iseg.eq.11) iflg11=0
      imn=imain(lev)
      isb=isub(lev)
      kseg=jsegm(lev)
      ksegmx=jsegnr(iseg)
      goto 200
c460   continue
c      iuwk=iuwkmn(levm)
c      juwk=iuwksb(levm)
c      acf=acoef(levm)
c      jsmt=jsm
c      ksbt=ksb
c      call caserr('attempt to call threex.')
c      call threex(nabc,nlwks,nuwks,puwk,iext1,iext2,iext3,iext4,iext5,
c     *ijadd,kadd,ladd,nelec,iext6,nrows,nij,nkl,iext7)
c      goto 200
480   continue
c      inxt=inxt+1
c      if(inxt.le.nxt) go to 190
      if (i.gt.1) goto 190
      go to 500
c490   call nxtblk
500   continue
      call mcfino(iwrite)
      if (intrac) write(iwrite,43434)
43434 format(' *** warning: interacting space code has been used')
      if (.not.diagon) call accnt(' ',2)
      if (.not.diagon) return
      diagon=.false.
      goto 10000
c      debug subchk,trace,init
c      at 130
c      trace on
      end
      subroutine mcpald(bfsym,bfkey,bfcode,a,b,s,
     #                  levpt,levnr,arc,nlwks,x,nelecs,spc,iwrite)
c
c***********************************************************************
c     this subroutine generates the distinct row table by searching    *
c     from the top of the graph down possible paths. note that there   *
c     is a program loop over orbitals running from statement labeled   *
c     '4'. the idea is for each level, search across all possible      *
c     a and b combinations (atest and btest) then across all points    *
c     on the row above (row) testing whether a particular arc (case)   *
c     gets down to the points given by atest,btest. the complication   *
c     lies in computing the excitation level when electrons are        *
c     excited into occupied or multi-reference orbitals. this is       *
c     accomplished using the x array which contains the number of      *
c     electrons excited into orbitals above the present point for      *
c     any walk leading to the particular point in question. therefore, *
c     points may be the same in all respects except for x value and    *
c     finally, at the fermi level the x value determines if the walk   *
c     is possible or not. for simple cases such as high spin open      *
c     shells the interacting space is implimented by counting spin-    *
c     flips as excitations into orbitals.                              *
c                                                                      *
c     the variables amax, amin, bmax, and bmin limit the portion       *
c     of the search to only those points that can be reached from      *
c     the present row. special orbitals are handled by input           *
c     explicitly the x values for all possible partial walks in the    *
c     special orbitals. originally the excitation level is set         *
c     to the sum of the general value and the value in the valence     *
c     (%) space. upon leaving the valence space it is reset to the     *
c     general value, thus giving say all singles and doubles           *
c     from references of all singles and doubles in the valence        *
c     space.                                                           *
c                                                                      *
c     maxb keeps track of the largest b value encountered. this        *
c     can then be used in the ci program to compute sufficient         *
c     coefficients for segment values. if this is done, there is       *
c     no limit to the value of the spin possible. finally the last     *
c     portion of the routine eliminates all points and arcs from       *
c     walks that dont make it from head to tail.                       *
c***********************************************************************
c
      implicit integer (a-z)
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
      common/entryp/orbtbf(1)
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
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     #,               levocc,spec,sspesh
     >, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     #,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     #,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     >,optio(8)
      character*3 codes
      common /drtcoc/ codes(9)
      integer dela,delb,delele,virtul,occupd,valocc,rescor,resvir
     >       ,frozen,valvir,opensh,multi,speshl,multrf,valenc
     >       ,fzc,fzv,cor,vir,doc,uoc,alp,bet,spe
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     1,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     2,               valvir,opensh,multi,speshl,multrf,valenc
     3,   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
c
      dimension bfsym(nbf),bfkey(nbf),bfcode(nbf,nrefs)
      dimension a(nrowsp),b(nrowsp),s(nrowsp),nlwks(nrowsp)
      dimension levpt(nlevs),levnr(nlevs),arc(nrws4p)
      dimension x(nrefs,nrowsp),nelecs(nrefs),spc(*)
      parameter (maxnr=2044)
      common /lsort / reno(maxnr)
c
      data spcnt /0/
c
c
      call accnt('mcpald  ',2)
c      print*,'levval,levfrm,levocc ',levval,levfrm,levocc
      maxb=-999999
      optio(2) = iabs(optio(2))
c     if (optio(2).eq.0) optio(2)=2
      optio(3)=optio(2)+optio(3)
      excita=optio(3)
      do 1 ref=1,nrefs
      nelecs(ref)=2*na+nb
      x(ref,1)=0
    1 continue
      levm1=nlevs
      levpt(levm1)=0
      levnr(levm1)=1
      a(1)=na
      b(1)=nb
      s(1)=ns
      x(1,1)=0
c
c     ----- here begins the loop over orbitals (levels) -----
c
    4 if (levm1.le.1) go to 25
      lev=levm1
      levm1=lev-1
      if (levm1.le.levocc) excita=optio(3)
      if (levm1.lt.levval) excita=optio(2)
      nrowm1=0
      bf=orbtbf(levm1)
      key=bfkey(bf)
      code=bfcode(bf,1)
c      print*,'looking for connections to level ',levm1,
c     >' orbital no. ',bf,key,code
c     >,' excita = ',excita
      nrowlv=levnr(lev)
      pontlv=levpt(lev)
      pontm1=pontlv+nrowlv
      levpt(levm1)=pontm1
      if (code.ne.spe) go to 5
      if (spcnt.eq.0) nelecs(1)=nelecs(1)-nespec
      spcnt=spcnt+1
      go to 11
    5 continue
      do 2 ref=1,nrefs
      nelecs(ref)=nelecs(ref)-delele(bfcode(bf,ref))
    2 continue
   11 continue
c
      amax=0
      amin=999999
      bmax=0
      bmin=999999
      kk1=pontlv+1
      kk2=pontlv+nrowlv
      do 3 row=kk1,kk2
      if (a(row).gt.amax) amax=a(row)
      if (a(row).lt.amin) amin=a(row)
      if (b(row).gt.bmax) bmax=b(row)
      if (b(row).lt.bmin) bmin=b(row)
    3 continue
      if (amin.gt.0) amin=amin-1
      if (bmin.gt.0) bmin=bmin-1
      if (bmax.gt.maxb) maxb=bmax
      bmax=bmax+1
c
      kk1=amax-amin+1
      if (kk1.lt.1) goto 2399
      do 23 ajunk=1,kk1
      atest=amax-ajunk+1
      kk2=bmax-bmin+1
      if (kk2.lt.1) goto 2299
      do 22 bjunk=1,kk2
      btest=bmax-bjunk+1
      kk3=pontlv+1
      kk4=pontlv+nrowlv
      if (kk4.lt.kk3) goto 2199
      do 21 row=kk3,kk4
      do 20 case=1,4
      if (code.eq.cor.and.case.ne.4) go to 19
      if (code.eq.vir.and.case.ne.1) go to 19
      ia=a(row)
      ib=b(row)
      is=s(row)
      ic=levm1-ia-ib
      go to (6,7,8,9),case
c
      write (iwrite,18) case
   18 format (//,' impossible case value in mcpald:',i5)
      call caserr('error exit of berkeley program')
c
    6 ic=ic-1
      dele=0
      go to 10
c
    7 ib=ib-1
      is = mults (is,bfsym(bf))
      dele=1
      go to 10
c
    8 ia=ia-1
      ib=ib+1
      ic=ic-1
      is = mults (is,bfsym(bf))
      dele=1
      if (code.eq.alp.and.optio(4).eq.0) dele=2
      go to 10
c
    9 continue
      ia=ia-1
      dele=2
c
   10 continue
      if (ia.lt.0.or.ib.lt.0.or.ic.lt.0) go to 19
      if (ia.ne.atest.or.ib.ne.btest) go to 19
      if (2*ia+ib.gt.excita.and.levm1.le.levfrm) go to 19
      possbl=0
      do 60 ref=1,nrefs
      if ((2*ia+ib.le.nelecs(ref)+excita).and
     #.   (2*ia+ib.ge.nelecs(ref)-excita)) possbl=1
   60 continue
      if (possbl.eq.0.0d0) go to 19
c
      kk5=pontm1+1
      kk6=pontm1+nrowm1
      if (kk6.lt.kk5) go to 1399
      do 13 rowm1=kk5,kk6
      if (ia.ne.a(rowm1).or.ib.ne.b(rowm1).or.is.ne.s(rowm1)) go to 13
      if (levm1.lt.levfrm) go to 17
      diff=0
      possbl=0
      if (levval-levfrm.eq.1) possbl=1
      if (code.ne.spe) go to 62
      if (spcnt.eq.spec) go to 61
      possbl=1
      ix=4*(x(1,row)+case-1)
      if (ix.ne.x(1,rowm1)) diff=1
      go to 63
   61 continue
      ref=1
      ix=spc(x(1,row)+case)
      frmx=0
      if (levm1.eq.levfrm) frmx=ix
      if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #.   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
      if (levm1.eq.levfrm) ix=0
      if (ix.ne.x(1,rowm1)) diff=1
      go to 63
   62 continue
      do 40 ref=1,nrefs
      ix=x(ref,row)
      if (dele.gt.delele(bfcode(bf,ref)).and.key.ne.valenc)
     #                           ix=ix+dele-delele(bfcode(bf,ref))
      if (lev.ne.levval) go to 70
      ix=ix+nelecs(ref)-(2*ia+ib)
      if (ib/2.gt.ix) ix=ib/2
      if (2*ia+ib.le.excita-ix) ix=0
   70 continue
      frmx=0
      if (levm1.eq.levfrm) frmx=ix
      if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #.   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
      if (levm1.eq.levfrm) ix=0
   40 if (ix.ne.x(ref,rowm1)) diff=1
   63 continue
      if (possbl.eq.0) go to 19
      if (diff.eq.0) go to 17
   13 continue
 1399 continue
c
c     ----- check that this is indeed a possible new point,     -----
c           especially, if got here because a, b or s different
c
      possbl=0
      if (levval-levfrm.eq.1) possbl=1
      if (code.ne.spe) go to 162
      if (spcnt.eq.spec) go to 161
      possbl=1
      go to 163
  161 continue
      ref=1
      ix=spc(x(1,row)+case)
      frmx=0
      if (levm1.eq.levfrm) frmx=ix
      if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #.   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
      go to 163
  162 continue
      do 140 ref=1,nrefs
      ix=x(ref,row)
      if (dele.gt.delele(bfcode(bf,ref)).and.key.ne.valenc)
     #                           ix=ix+dele-delele(bfcode(bf,ref))
      if (lev.ne.levval) go to 170
      ix=ix+nelecs(ref)-(2*ia+ib)
      if (ib/2.gt.ix) ix=ib/2
      if (2*ia+ib.le.excita-ix) ix=0
  170 continue
      frmx=0
      if (levm1.eq.levfrm) frmx=ix
      if ((2*ia+ib.le.nelecs(ref)+excita-ix).and
     #.   (2*ia+ib.ge.nelecs(ref)-excita+frmx)) possbl=1
  140 continue
  163 continue
      if (possbl.eq.0) go to 19
c
      nrowm1=nrowm1+1
      rowm1=pontm1+nrowm1
      if (rowm1.lt.nrowsp) go to 16
      write (iwrite,15) nrowsp,levm1
   15 format (//,' not enough space to make drt, have only room for'
     #,           i6,' rows and are exceeding that at level',i4)
      call caserr('error exit of berkeley program')
   16 a(rowm1)=ia
      b(rowm1)=ib
      s(rowm1)=is
c
c      print*,'creating row no. ',rowm1,' levm1,a,b,s =',levm1,ia,ib,is,
c     > ' pointed by row ',row,' with case ',case
      if (levm1.le.levfrm) go to 50
      if (lev.ne.levval) go to 42
      do 77 ref=1,nrefs
      ix=nelecs(ref)-(2*ia+ib)
      if (ib/2.gt.ix) ix=ib/2
      if (2*ia+ib.le.excita-ix) ix=0
      x(ref,rowm1)=x(ref,row)+ix
   77 continue
      go to 17
c
   42 continue
      if (code.ne.spe) go to 44
      if (spcnt.eq.spec) go to 43
      x(1,rowm1)=4*(x(1,row)+case-1)
      go to 17
   43 continue
      x(1,rowm1)=spc(x(1,row)+case)
      go to 17
   44 continue
      do 41 ref=1,nrefs
      ix=x(ref,row)
      if (dele.gt.delele(bfcode(bf,ref)).and.key.ne.valenc)
     #                         ix=ix+dele-delele(bfcode(bf,ref))
      if (2*ia+ib.le.excita-ix) ix=0
      x(ref,rowm1)=ix
   41 continue
      go to 17
   50 continue
      do 51 ref=1,nrefs
      x(ref,rowm1)=0
   51 continue
   17 continue
      arc((row-1)*4+case)=rowm1-pontm1
   19 continue
   20 continue
   21 continue
2199  continue
   22 continue
2299  continue
   23 continue
2399  continue
      levnr(levm1)=nrowm1
      go to 4
   25 continue
c
c     ----- set weight of bottom of shavitt graph to one, -----
c                  eliminate all other bottoms
c
      nroot=0
      kk1=levpt(1)+1
      kk2=levpt(1)+levnr(1)
      do 27 root=kk1,kk2
      if (a(root).ne.0.or.b(root).ne.0.or.s(root).ne.1) go to 26
      nroot=nroot+1
      nlwks(root)=1
   26 continue
   27 continue
c
      if (nroot.eq.1) go to 29
      write (iwrite,28) nroot
   28 format (//,' invalid number of bottoms to graph:',i4,//)
      call caserr('error exit of berkeley program')
   29 continue
c
c     ----- generate the weights (nlwks) of all rows -----
c
      if (lev.lt.2) go to 91
      do 32 lev=2,nlevs
      levm1=lev-1
      pontm1=levpt(levm1)
      kk1=levpt(lev)+1
      kk2=levpt(lev)+levnr(lev)
      do 31 row=kk1,kk2
      nlwk=0
      do 30 case=1,4
      arcpt=(row-1)*4+case
      if (arc(arcpt).gt.0) nlwk=nlwk+nlwks(arc(arcpt)+pontm1)
   30 continue
      nlwks(row)=nlwk
   31 continue
   32 continue
c
      nwks=nlwks(1)
      if (nwks.gt.0) goto 95
      if (optio(4).ne.0) go to 93
      write (iwrite,92) nwks
   92 format (//,t25,'***** number of interacting walks *****',/,t25,'*'
     #,      t63,'*',/,t25,'*',t35,i12,t63,'*',/,t25,'*',t63,'*',/
     #,          t25,'***************************************')
      go to 95
   93 continue
      write (iwrite,94) nwks
   94 format (//,t25,'******** total number of walks ********',/t25,'*'
     #,       t63,'*',/,t25,'*',t35,i12,t63,'*',/,t25,'*',t63,'*',/
     #,          t25,'***************************************')
   95 continue
c
      kk2=levpt(2)+levnr(2)
c      do 999 i=1,nrefs
c      do 999 j=1,kk2
c  999 if (x(i,j).ne.0) write (iwrite,998) i,j,x(i,j)
c  998 format (3i6)
c      write (iwrite,998) kk2
c     ----- remove all rows with zero weights -----
c
      pont=1
      kk1=nlevs-1
      do 37 junk=1,kk1
      lev=nlevs-junk
      nrowlv=levnr(lev)
      pontlv=levpt(lev)
      levpt(lev)=pont
      kk2=pontlv+1
      kk3=pontlv+nrowlv
      do 36 row=kk2,kk3
      if (nlwks(row).eq.0) go to 35
      pont=pont+1
      a(pont)=a(row)
      b(pont)=b(row)
      s(pont)=s(row)
      nlwks(pont)=nlwks(row)
      do 34 case=1,4
      arc((pont-1)*4+case)=arc((row-1)*4+case)
   34 continue
      kk4=levpt(lev+1)+1
      kk5=levpt(lev+1)+levnr(lev+1)
      do 55 rowp1=kk4,kk5
      kk6=(rowp1-1)*4+1
      kk7=(rowp1-1)*4+4
      do 54 case=kk6,kk7
      if (arc(case).eq.row-pontlv) arc(case)=pont-levpt(lev)
   54 continue
   55 continue
      go to 36
   35 continue
      kk4=levpt(lev+1)+1
      kk5=levpt(lev+1)+levnr(lev+1)
      do 53 rowp1=kk4,kk5
      kk6=(rowp1-1)*4+1
      kk7=(rowp1-1)*4+4
      do 52 case=kk6,kk7
      if (arc(case).eq.row-pontlv) arc(case)=0
   52 continue
   53 continue
   36 continue
      levnr(lev)=pont-levpt(lev)
   37 continue
c
   91 continue
      nwks=nlwks(1)
      nrows=pont
      lev=2
c...  pjk  1/3/85 to trap one active orbital
      if (nlevs.le.2) goto 99
  105 nr=levnr(lev)
      npt=levpt(lev)+1
      nrb=1
      neq=0
      nra=npt
      nrd=nra
      reno(1)=1
      if (nr.gt.maxnr) call caserr('table overflow - too many rows at on
     >e level')
      do 100 i=2,nr
      npt=npt+1
      do 107 j=1,nrb
      if (a(npt).eq.a(nrd+j-1)) go to 90
      go to 107
   90 if (b(npt).eq.b(nrd+j-1)) go to 792
      go to 107
  792 if (s(npt).eq.s(nrd+j-1)) go to 793
      go to 107
  793 do 108 k=1,4
      if (arc((npt-1)*4+k).eq.arc((nrd+j-2)*4+k)) go to 108
      go to 107
  108 continue
      j1=j
      go to 791
  107 continue
      nra=nra+1
      nrb=nrb+1
      a(nra)=a(npt)
      b(nra)=b(npt)
      s(nra)=s(npt)
      nlwks(nra)=nlwks(npt)
      do 109 k=1,4
  109 arc((nra-1)*4+k)=arc((npt-1)*4+k)
      reno(i)=nrb
      go to 100
  791 neq=neq+1
      reno(i)=j1
  100 continue
      levnr(lev)=nrb
      lev1=lev-1
      nrc=levpt(lev1)+1
      do 110 k=nrc,nrows
      k1=k-levpt(lev1)
      a(nra+k1)=a(k)
      b(nra+k1)=b(k)
      s(nra+k1)=s(k)
      nlwks(nra+k1)=nlwks(k)
      do 110 j=1,4
  110 arc((nra+k1-1)*4+j)=arc((k-1)*4+j)
      do 112 k=1,lev1
  112 levpt(k)=levpt(k)-neq
      nrows=nrows-neq
      lev1=lev+1
      if (lev1.eq.nlevs) go to 99
      nr=levnr(lev1)
      npt=levpt(lev1)
      if (nr.eq.1) go to 115
      do 113 i=1,nr
      npt=npt+1
      do 114 j=1,4
      if (arc((npt-1)*4+j).eq.0) go to 114
      arc((npt-1)*4+j)=reno(arc((npt-1)*4+j))
  114 continue
  113 continue
  115 lev=lev+1
      go to 105
   99 continue
      nrows4=nrows*4
      nrowoc=0
      do 38 lev=levfrm,nlevs
      nrowoc=nrowoc+levnr(lev)
   38 continue
      nrow4o=nrowoc*4
      call accnt('other   ',2)
      return
      end
      subroutine mcreor(ntype,bfnum,numsym,bfsym,bfcode,orbtbf
     1,                 iout,orbsym)
c
c***********************************************************************
c     this subroutine determines the order of the orbitals in the ci   *
c     portion of the calculation. the order is set in the block data   *
c     drtcod -- 1 appears at the bottom of te graph, highest number    *
c     at the top. iout(n) gives ci number of nth scf orbital, with     *
c     a -1 for frozen cores and 0 for frozen virtual. orbtbf(n) gives  *
c     the scf number for the nth ci orbital.                           *
c***********************************************************************
c
      implicit integer (a-z)
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
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,iop(8)
      character*3 codes
      common /drtcoc/ codes(9)
      integer dela,delb,delele,virtul,occupd,valocc,rescor,resvir
     >       ,frozen,valvir,opensh,multi,speshl,multrf,valenc
     >       ,fzc,fzv,cor,vir,doc,uoc,alp,bet,spe
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     1,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     2,               valvir,opensh,multi,speshl,multrf,valenc
     3,   fzc, fzv, cor, vir, doc, uoc, alp, bet, spe
c
      common/entryp/orbtyp(maxorb)
      dimension ntype(ntypes),bfnum(ntypes,nbf),numsym(ntypes,nsym)
      dimension bfsym(nbf),bfcode(nbf,nrefs),orbtbf(norbs)
      dimension iout(nbf),orbsym(norbs)
c
      orb=0
      na=0
      nb=0
      ns=1
      levfrm=0
      levocc=-999999
      levopn=999999
      levmul=999999
      levval=999999
c
      do 20 type=1,ntypes
      ntype(type)=0
      do 10 sym=1,nsym
10    numsym(type,sym)=0
      do 20 bf=1,nbf
20    bfnum(type,bf)=0
c
      do 30 bf=1,nbf
      type=occupd
      cod=bfcode(bf,1)
      if (cod.eq.spe) type=speshl
      if (cod.eq.uoc) type=virtul
      if (cod.eq.alp) type=opensh
      if (cod.eq.bet) type=opensh
      if (bfkey(bf).eq.valenc) type=valocc
      if (bfkey(bf).eq.valenc .and. cod.eq.uoc) type=valvir
      if (bfkey(bf).eq.multrf) type=multi
      ntype(type)=ntype(type)+1
      bfnum(type,ntype(type))=bf
      sym=itypea(bf)
      numsym(type,sym)=numsym(type,sym)+1
      bfsym(bf)=sym
30    continue
c
      do 50 type=1,ntypes
      ntp=ntype(type)
      if (ntp.eq.0) goto 51
      if (type.eq.speshl.and.levopn.eq.999999) levopn=orb+1
      if (type.eq.opensh) levopn=orb+1
      if (type.eq.multi ) levmul=orb+1
      if (type.eq.valocc) levocc=orb+1
      if (type.eq.valvir) levval=orb+2
      do 40 junk=1,ntp
      num=ntp-junk+1
      bf=bfnum(type,num)
      sym = bfsym(bf)
      code=bfcode(bf,1)
      orb=orb+1
      iout(bf)=orb
      orbsym(orb)=sym
      na=na+dela(code)
      nb=nb+delb(code)
      if (delele(code).eq.1) ns=mults(ns,sym)
      orbtbf(orb)=bf
40    continue
51    if (type.eq.virtul) levfrm = orb+1
50    continue
c
      do 60 bf=1,orb
 60   orbtyp(bf)=orbtbf(bf)
      if (levocc.eq.-999999) levocc=levfrm
      orbfrm=levfrm-1
c
      return
      end
      subroutine testft(iq,iwrite)
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
      common /couple/ surd(511)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     >              ,jorb(511),korb(511),lorb(511)
      dimension iq(1),ibra(32),iket(32),idiff(32)
c
      nocc=nact*nci
      iocc = icori(nocc)
      do 111 i=1,nocc
111   iq(i+iocc-1) = 0
      call rdedx(surd,511,iblkft,numft)
      call find (numft)
      do 200 i=1,nblkd1
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 3 iw=1,nw
      io=iorb(iw)
      iii=ii(iw)
      val=surd(ival(iw))+.1d0
3     iq(iocc-1+(io-1)*nci+iii)=val
200   continue
      call get (g,nw)
      junk = nblkd1+nblko1+2+iblkft
      call search (junk,numft)
      callfind (numft)
c      do 1111 i=1,nci
c1111  write (iwrite,*) i, (iq(iocc-1+(j-1)*nci+i),j=1,nact)
c
c...   check 2particle forms
70    call get(g,nw)
      if (nw.eq.0) goto 101
      call find (numft)
      call unpkft(nw)
      do 13 iw=1,nw
      io=iorb(iw)
      jo=jorb(iw)
      ko=korb(iw)
      lo=lorb(iw)
      iii=ii(iw)
      jjj=jj(iw)
      do 71 i=1,nact
      ibra(i)=iq(iocc-1+(i-1)*nci+iii)
      iket(i)=iq(iocc-1+(i-1)*nci+jjj)
71    idiff(i)=iket(i)-ibra(i)
      idiff(io)=idiff(io)+1
      idiff(jo)=idiff(jo)-1
      idiff(ko)=idiff(ko)+1
      idiff(lo)=idiff(lo)-1
      do 72 i=1,nact
      if (idiff(i).ne.0) goto 73
72    continue
      goto 13
73     write (iwrite,555) iii,io,jo,ko,lo,jjj
555   format(' formula tape has the formula <',i3,'/',4i2,'/',i3,'>')
      write (iwrite,777) iii,(ibra(i),i=1,nact)
777   format(' bra',i3,' has occupation numbers',32i2)
      write (iwrite,778) jjj,(iket(i),i=1,nact)
778   format(' ket',i3,' has occupation numbers',32i2)
13    continue
      goto 70
101   continue
      return
      end
      subroutine ver_mcscfc(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfc.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
