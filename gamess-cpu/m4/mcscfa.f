c 
c  $Author: jmht $
c  $Date: 2010-08-10 16:57:58 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6178 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfa.m,v $
c  $State: Exp $
c  
c ******************************************************
c ******************************************************
c             =   mcstart     =
c ******************************************************
c ******************************************************
      subroutine mcstar(q,iq)
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
      logical btest, ovect
c
      integer iwrit,irea,iblkw,iblkr,iou2,lenou2,inn2,lenin2
      common/count/ iwrit,irea,iblkw,iblkr,iou2,lenou2,inn2,lenin2
c
      common/intbu2/ iblff
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      common /mcopt / var,varc,thzr,one
      common /mccore/ intrel,lword,ltop,lmax,lmin,mreal
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      character*132 zedfil,zftfil
      character*4  yed, yft
      character*8  zedstat, zftn, zftstat
      common/discc/   yed(maxlfn),zedfil(maxlfn),zedstat(maxlfn),
     *   zftn(maxfrt),yft(maxfrt),zftfil(maxfrt),zftstat(maxfrt)
      common/rotat/akkkk(60),nswap
      common/linkmc/irest,ivect,iorbit
     * , iblrst,numrst,isecrs,ivgues,inp,iblocp,iblkg
      integer orbfrm,symorb,optio,spec,sspesh
      common /drtinf/ na,nb,ns,nespec,maxb,levfrm,levval,levopn,levmul
     1,               levocc,spec,sspesh
     2, nbf,nsym,norbs,nrowsp,nrws4p,nrows,nrows4
     3,               nlevs,nrefs,nrowoc,nrow4o,nwks,nwksoc,nlevoc
     4,               orbfrm,symorb,numij,ngroup,numint,nmax,nspc
     5,optio(8)
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
      common /lsort / ablkin(511)
      dimension q(*),iq(*),iblkin(2,511)
      equivalence (ablkin(1),iblkin(1,1))
      integer iqpos(8)
      character*8 itdir(7)
      data itdir /'diagci  ','uncouple','werner  ','augment ','newton  '
     1          ,'null    ','internal'/
c
c
      if(mcacct)call accnt('other   ',0)
c
      call accnt('mcstart',1)
      nitdir=7
      iblf=1
      iblf1=1
      iblf2=0
c
      one=1.0d0
      thzr=1.d-6
       iou2=0
       lenou2=0
       inn2=0
       lenin2=0
      iblff=1
c
      ibase=icorr(0)
      npread=0
      if (mcprin) write(iwr,20)(ztitle(i),i=1,9)
 20    format(/
     *2x,78('*')/
     *2x,'*',76x,'*'/
     *2x,'*  ',9a8,'  *'/
     *2x,'*',76x,'*'/
     *2x,78('*')/)
c...  read in symmetry info
      call getbas
      do 60 i=1,maxorb
      if (symaos(i).lt.1) goto 70
60    nsymao(symaos(i))=nsymao(symaos(i))+1
70    ns=1
      do 80 i=1,nirrr
      nsymm(i)=nsymao(i)
      istart(i)=ns
      ns=ns+nsymao(i)
80    mfin(i)=ns-1
      nbasao=ns-1
      nbasis=nbasao
      lenbas = (nbasao*(nbasao+1))/2
      nblkq=nbasao**2
      if(irest.eq.0.and.iorbit.eq.0) 
     + call caserr('orbital directive essential in startup jobs.')
      ovect = .false.
      if (zruntp.eq.'force') ovect = .true.
c ****
c     in force constant jobs we have a problem, for itype must
c     be reset as the point group changes. Fudge this for the
c     moment by analysing the input mos and resetting in orbita
c     triggered by ovect = .true.
c
c     if(iorbit.ne.0) then
      if (ovect) then
c
c     obtain orbital info. from symmetry of input mos rather
c     than using =orbital= data
c
        if(irest.eq.0.and.ivgues.eq.1)
     +   call orbrd(q(1),ivec,ivect,isecrs,iwr)
        if(irest.ne.0.and.ivect.eq.0)
     +     call orbrd(q(1),ivec,ivect,isecrs,iwr)
c
        if(ivect.lt.1)call caserr('trial vectors not supplied')
        if(isec.le.0) call caserr('vectors section not supplied')
c
        ibuf=icori(nbasao)
        call symmo(q(ivec),iq(ibuf))
c       if(mcprin) write (6,*) 'mcstar: symmetry of input mos: ', 
c    +              (iq(ibuf+loop-1), loop=1,nbasao)
        call orbita(iwr,iq(ibuf),.true.,mcprin)
        call corlsi(ibuf)
      else
       ibuf = 1
       call orbita(iwr,iq(ibuf),.false.,mcprin)
      endif
c     endif
c ****
      if(inp.ne.0)call pspace(q,iq)
c
      if(irest.eq.0.and.ivgues.eq.1)
     * call orbrd(q(1),ivec,ivect,isecrs,iwr)
      if(irest.ne.0.and.ivect.eq.0)
     1   call orbrd(q(1),ivec,ivect,isecrs,iwr)
c
      if(ivect.lt.1)call caserr('trial vectors not supplied')
      if(isec.le.0) call caserr('vectors section not supplied')
c
      if(nswap.gt.0) then
        ibuf=icori(nbasao)
        call symmo(q(ivec),iq(ibuf))
        call mcswap(q(ivec),iq(ibuf),iwr)
        call corlsi(ibuf)
      endif
c
      if(((isec.ne.0.and.isec.lt.1000).or.
     1   (isecd.ne.0.and.isecd.lt.1000).or.
     2   (isecn.ne.0.and.isecn.lt.1000).or.
     3   (isecnc.ne.0.and.isecnc.lt.1000)).
     4    and.mcprin)
     4    write(iwr,560)iblk3,yed(num3)
560   format(/1x,'dumpfile at block',i5,' on',t36,a4)
      if (ifsort.eq.0.and.mcprin) write(iwr,570) yed(numa)
570   format(1x,'input two electron integral file',t36,a4)
      if (mcprin) write(iwr,580)yed(num2)
580   format(1x,'sorted integral file',t36,a4)
      if (mcprin) write(iwr,590) yed(num4),yed(num6)
590   format(1x,'scratch files',t36,2a4)
      if(iexc.ge.0.and.mcprin) write(iwr,600) yed(numft)
600   format(1x,'ci formula file',t36,a4)
      if (mcprin) write(iwr,610) isec
610   format(1x,'molecular orbital dump at section',t36,i3)
c620  format(1x,'mcscf dump at section',t36,i3)
c
      if (mcprin) write(iwr,630)(j,j=iter+1,min(maxcyc,39)),40
630   format(1x,'iteration controls:'/11x,40i3)
      do 635 i=1,nitdir
 635  if (mcprin) write(iwr,640)
     +  itdir(i),(btest(itinfo(j),i),j=iter+1,min(maxcyc,39)),
     +            btest(itinfo(40),i)
640   format(1x,a8,2x,40l3)
      if (mcprin) write(iwr,650)conv,econv,sconv
650   format(/1x,'convergence thresholds',e10.2,'(gradient)'/
     +        1x,'                      ',e10.2,'(energy)'/
     +        1x,'                      ',e10.2,'(step length)')
      if (mcprin) write(iwr,660)maxcyc
660   format(1x,'maximum number of cycles',i4)
      call lkpgen(iwr)
      call orbsym (q(1),iq(1),q(ivec))
      call dmpini
c..   symmetry pack, orthogonalise & dump vectors
      do 670 isym=1,nirrr
      iad = icorr(nsymao(isym)**2)
      iqpos(isym) = iad
      do 670 i=1,nbasao
      if(i.gt.nfreez) then
      is=itype(i-nfreez)
      else
      is=ifzsym(i)
      end if
      if(is.ne.isym) goto 670
      ix1=ivec-1+(i-1)*nbasao+istart(isym)
      call dcopy(nsymao(isym),q(ix1),1,q(iad),1)
      iad = iad + nsymao(isym)
670   continue
c      print*,'nfreez,ncoremc,nact,nprim:  ',nfreez,ncoremc,
c     +        nact,nprim
c      print*,'nbasis,nbasao: ',nbasis,nbasao
c      print*,'lenbas,nblkq:  ',lenbas,nblkq
c      print*,'ifreez: ',(ifreez(i),i=1,nirrr)
c      print*,'ncor:   ',(ncor(i),i=1,nirrr)
c      print*,'nprm:   ',(nprm(i),i=1,nirrr)
c      print*,'ifzsym: ',(ifzsym(i),i=1,nfreez)
c      print*,'itype:  ',(itype(i),i=1,nbasis)
c      print*,'itypea: ',(itypea(i),i=1,nact)
c      print*,'nsymm:  ',(nsymm(i),i=1,nirrr)
c      print*,'nsymao: ',(nsymao(i),i=1,nirrr)
      call qput (q(1),iqpos,iwr)
      call corlsr (iqpos(1))
c
c...  process cigues
      if (iguess.ge.0) goto 690
      iguess=nstate
      junk=nci*iguess
      iagues=icorr(junk)
      call rdedx(q(iagues),junk,iblkg,num8)
      call cidmpi
      call mcdump
      call cput(q(iagues),nstate)
690   continue
c
      if (iexc.lt.0 .and. iexcv.eq.0) then
c.....determinant casscf
      iexc = -1
      iffmtp = 1
      end if
c
      call corlsr (ibase)
      call accnt(' ',1)
      return
      end
      subroutine orbrd(q,ivec,ivect,isecrs,iwrite)
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
      parameter (mxorb1=maxorb+1)
      logical iftran
      character *8 param,gitle
      common /junkc/ param(19),gitle(10)
      common/blkorbs/evalue(maxorb),occ(mxorb1),
     1                    nbas,newb,ncoll,ivalue,ioccc,ipad
      common/junk/ilifc(maxorb),ntranc(maxorb),itranc(mxorb3),
     * ctran(mxorb3),iftran,iftri
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
c
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      dimension q(*)
      lench=29
      lenhed=mach(8)
      lenctr=mach(9)
c     nhed=1+lensec(lenhed)+lensec(lenctr)
      call secget (isecrs,3,iblkvs)
      call rdchr(param,lench,iblkvs,num3)
      call reads(evalue,lenhed,num3)
      nav = lenwrd()
      call readis (ilifc,lenctr*nav,num3)
c...         HARMONIC
      nblkq1 = newbas1*newbas1
      ivec = icorr (nblkq1)
      call reads (q(ivec),nblkq1,num3)
      call comharm(q(ivec),'vectors',ilifq)
c...
      if( mcprin) write(iwrite,90) isecrs
90    format(1x,'vectors restored from section',i4,
     1' of dumpfile')
c...  orbitals now loaded
      ivect=1
      return
      end
      subroutine mcdump
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
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      call wrt3 (radius,mach(18),iblkdmc,num3)
      call clredx
      return
      end
      subroutine qget (q,iqpos)
c    load mo vectors by symmetries in preparation for 4-index transforma
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
      dimension q(*),iqpos(nirrr)
      iqq = icorr(nblkq)
      call rdedx(q(iqq),nblkq,iblkq,num3)
      do 10 isym=1,nirrr
      nsym=nsymao(isym)
      iad=iqpos(isym)
      do 10 i=1,nbasao
      if(i.gt.nfreez) then
      is=itype(i-nfreez)
      else
      is=ifzsym(i)
      end if
      if(is.ne.isym.or.nsym.eq.0) goto 10
      call dcopy(nsym,q(iqq-1+(i-1)*nbasao+istart(isym)),1,q(iad),1)
      iad=iad+nsym
10    continue
      call corlsr (iqq)
      return
      end
      subroutine qput (q,iqpos,iwrite)
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      common /block / nosymm,nosym(8),iosym(8,maxorb)
      common /lsort / iisym(maxorb),iisym1(maxorb)
      dimension q(*),iqpos(8)
      save ifir
      data ifir/1/
c
c...  apply schmidt orthogonalisation
      ibase = icorr(0)
      len = 0
      do 998 isym=1,nirrr
      nosym(isym) = 0
998   len = len + (nsymao(isym)*(nsymao(isym)+1))/2
      ibf = icorr(len)
      ibf2 = icorr(len)
      call get1 (q(1),q(ibf),1,dummy)
      call get1 (q(1),q(ibf2),3,dummy)
      iv = icorr(maxbas**2)
      nosymm=0
      do 50 isym=1,nirrr
      ifjr=ifir
      iaa=iqpos(isym)
      n = nsymao(isym)
      if(n.le.0) go to 50
      ihh = icorr(n**2)
      call square (q(ihh),q(ibf2),n,n)
      ibf2 = ibf2 + (n*(n+1))/2
c...  search for additionaly point group symmetry
c      call outsqr (q(ihh),n,n,n,'t+v matrix in qput')
c      call outsqr (q(iaa),n,n,n,'initial orbitals in qput')
      do 11 i=1,n
11    iisym(i) = -i
      do 1 i=1,n
      do 2 j=1,n
      if ( dabs(q(ihh+(i-1)*n+j-1)).lt.1.0d-12) goto 2
c      if ( abs(q(ihh+(i-1)*n+j-1)).lt.1.0e-8) print*,
c     >'*** warning: small t+v matrix element ',i,j,
c     >q(ihh+(i-1)*n+j-1)
      jjsym = iisym(j)
      do 8 k=1,n
      if (iisym(k).eq.jjsym) iisym(k) = iisym(i)
8     continue
2     continue
1     continue
c...  reorder symmetry labels consecutively
      do 83 i=1,n
      if (iisym(i).gt.0) goto 83
      jjsym = iisym(i)
      nosym(isym) = nosym(isym)+1
      do 831 j=1,n
      if (iisym(j).eq.jjsym) iisym(j) = nosym(isym)
831   continue
83    continue
      if (nosym(isym).gt.1.and.ifir.eq.1) write(iwrite,25)isym,
     *(iisym(i),i=1,n)
25    format(/' extra symmetry of aos in symmetry',i2,':',
     >   (t40,4(10i2,2x)))
      do 26 i=1,n
26    iisym1(i)=iisym(i)
      ia=iaa
      isa=iv
      iss=icorr(n**2)
      call square(q(iss),q(ibf),n,n)
      ibf = ibf +n*(n+1)/2
      do 30 i=1,n
c...  clean up according to extra symmetry this orbital
      zz = -99.0d0
      do 3 j=1,n
      if ( dabs(q(ia+j-1)).lt.zz) goto 3
      if(iisym1(j).eq.0) goto 3
      zz =  dabs(q(ia+j-1))
      jmax=j
      jjsym = iisym(j)
3     continue
      iisym1(jmax)=0
      amount = 0.0d0
      zz = 0.0d0
      do 4 j=1,n
      zz = zz +  dabs(q(ia+j-1))
      if (iisym(j).eq.jjsym .or. q(ia+j-1).eq.0.0d0) goto 4
c      print 44,i,isym,j,q(ia+j-1)
c44    format(' *** in orbital',i3,' of symmetry',i2,
c     >' coefficient of basis function',i3,' =',d8.2,
c     <' and is now being set to zero')
      amount = amount +  dabs(q(ia+j-1))
      q(ia+j-1) = 0.0d0
4     continue
      amount = amount / zz
      if (amount.gt.1.0d-5) write(iwrite,791)isym,i-ifreez(isym),amount
791   format(' *** in symmetry',i2,' orbital',i3,
     >' symmetry contamination of',e10.3,' has been removed')
      if(i.gt.ifreez(isym)) then
      if(ifir.eq.0.and.jjsym.ne.iosym(isym,i-ifreez(isym))) ifjr=1
      iosym(isym,i-ifreez(isym)) = jjsym
      end if
      ja=iaa
      jsa=iv
      do 20 j=1,i-1
      zz=-ddot(n,q(ia),1,q(jsa),1)
      call daxpy(n,zz,q(ja),1,q(ia),1)
      jsa=jsa+n
20    ja=ja+n
      call mxmaa(q(iss),1,n,q(ia),1,0,q(isa),1,0,n,n,1)
      zz=ddot(n,q(ia),1,q(isa),1)
      if (zz.le.0.0d0)call caserr ('schmidt orthogonalisation singular')
      zz=1.0d0/ dsqrt(zz)
      call dscal(n,zz,q(isa),1)
      call dscal(n,zz,q(ia),1)
      isa=isa+n
30    ia=ia+n
c      if (n.gt.0) call outsqr (q(iaa),n,n,n,'new mos in qput')
      if(nosym(isym).gt.1.and.ifjr.eq.1)write(iwrite,255)isym,
     > (0,i=1,ifreez(isym)),(iosym(isym,i),i=1,nsymm(isym))
255   format(' extra symmetry of mos in symmetry',i2,':',
     >   (t40,4(10i2,2x)))
      nosymm = max(nosymm,nosym(isym))
      call corlsr (ihh)
 50   continue
      call corlsr (ibase)
c
      ifir=0
      iqq = icorr(nblkq)
      call vclr(q(iqq),1,nblkq)
      do 60 isym=1,nirrr
      iad = iqpos(isym)
      do 60 i=1,nbasao
      if(i.gt.nfreez) then
      is=itype(i-nfreez)
      else
      is=ifzsym(i)
      end if
      if(is.ne.isym.or.nsymao(isym).eq.0) goto 60
      ix1 = iqq-1+(i-1)*nbasao+istart(isym)
      call dcopy(nsymao(isym),q(iad),1,q(ix1),1)
      iad = iad + nsymao(isym)
60    continue
      call wrt3 (q(iqq),nblkq,iblkq,num3)
      call corlsr (iqq)
      return
      end
      subroutine sort(q,iq,iwrite)
      implicit real*8  (a-h,o-z)
      logical fock
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
      common /stak / btri,mlow,nstack,iblock
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
      common /bufb  / nwbnwb,lnklnk,gsort(5118)
      common /junk  / klsort(2,3412)
      dimension q(*),iq(*),iqpos(8)
      common /intbuf/ intpos,intfil
      common /mccore/ intrel
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      call accnt('int sort',1)
      intpos=-1
      intfil=num2
      iblf2=iblk2
      isym12=1
      master=0
      lentrr=lentra(1)
      maxtrr=lentrr
      if(o255i) then
       ibl5=((iblsrt-1)*2)/3
      else
       ibl5=((iblsrt-1)*2)/4
      endif
      ifck=icorr(nbasao**2)
      iden=ifck
      if(nfreez.ne.0) then
      iqpos(1)=ifck
      do 4 i=1,nirrr-1
4     iqpos(i+1)=iqpos(i)+nsymao(i)**2
c.....frozen core density
      iden=icorr(nbasao**2)
      call qget(q(1),iqpos)
      call vclr(q(iden),1,nblkq)
      do 5 i=1,nfreez
      is=ifzsym(i)
      iaa=iqpos(is)
      iqpos(is)=iqpos(is)+nsymao(is)
      idd=iden+(istart(is)-1)*(nbasao+1)
5     call mxmb(q(iaa),1,0,q(iaa),0,1,q(idd),1,nbasao,
     1    nsymao(is),1,nsymao(is))
      end if
      call vclr(q(ifck),1,nblkq)
      maxt = icorrm() /lentri
      nirrnb = nirrr*nbasao
      nword = (icorrm() - (nirrnb+1)/intrel - 1) / (1+2/intrel)
      ires = min(nword/ibl5,200)
      if (mcprin) write(iwrite,10)ibl5,ires,maxt
10    format(1x,'bucket size =',i5/
     +       1x,'max. no. of buckets =',i4/
     +       1x,'max. no. of matrices per bucket =',i5)
      if(ires.lt.1.or.maxt.lt.1)callcaserr('not enough store.')
c...   determine min. no. of passes for sort1/calc1
      npassa=0
20    npassa=npassa+1
      nteff=(lenbas-1)/npassa+1
      if ((nteff-1)/ires.ge.maxt) goto 20
      npassa=min(npassa,lenbas)
      fock=nfreez.gt.0
      do 120 ipass=1,npassa
c...  no. of buckets
      nbuck = ires
30    if ((min(nteff,lenbas-master)-1)/nbuck+1 .gt. maxt) goto 40
      nbuck = nbuck-1
      if (nbuck) 40,40,30
40    nbuck = nbuck + 1
      nword = ibl5*nbuck
      ia1 = icorr(-nword)
      ia2 = icori(-(nword+nword))
      ia3 = icori(-nirrnb)
      call sorta (q(ia1),iq(ia2),iq(ia3),q(iden),q(ifck),fock)
      call stopbk
      call corlsr (ia1)
      if(ipass.eq.1) then
      fock=.false.
      ifc=icorr(ltri)
      iho=icorr(lenbas)
      call fockrd(q(1),q(ifck),q(ifc),q(iho),q(iden))
      call corlsr(ifc)
      end if
      nsize=ntri*lentri
      ia1 = icorr(-nsize)
      ia1m=ia1-1
      do 110 ibuck=1,nbuck
      mhigh=min(mloww+ntri,mhi)
      mtri=mhigh-mloww
      mloww=mloww+1
      call vclr(q(ia1),1,nsize)
      mk=mark(ibuck)
50    if (mk.eq.9999999) goto 70
      call mcsrti (mk,nk)
      do 60 iword=1,nk
60    q(ia1m+(klsort(1,iword)-mloww)*lentri+klsort(2,iword))=
     +        gsort(iword)
      goto 50
70    map=ia1
      do 100 itri=1,mtri
      master=master+1
80    if (master.le.maxtrr) goto 90
      isym12=isym12+1
      lentrr=lentra(isym12)
      maxtrr=maxtrr+lentrr
      goto 80
90    call intou (q(map),lentrr)
100   map=map+lentri
110   mloww=mhigh
      call corlsr (ia1)
120   if (mcprin) write(iwrite,130) ipass,nbuck,master
130   format(1x,'end of sort pass',i3/
     +       1x,'no. of buckets used =       ',i6/
     +       1x,'no. of matrices completed = ',i6)
      call intend
      nblk=iposun(num2)-iblk2
      if (mcprin) write(iwrite,140)nblk
140   format(1x,'length of sorted integral file =',i7,' blocks')
      ifsort=1
      call corlsr(ifck)
      call accnt('other   ',1)
      return
      end
      subroutine fockrd(q,fock,fc,ho,dc)
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
      common /lsort / iisym(maxorb)
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
      dimension q(*),fock(nbasao,nbasao),dc(nbasao,nbasao),fc(*),ho(*)
      efreez=0
      call vclr(fc,1,ltri)
      if(nfreez.ne.0) then
      ij=0
      do 20 i=1,nbasao
      is=symaos(i)
      do 10 j=istart(is),i-1
      ij=ij+1
      fc(ij)=fock(i,j)+fock(j,i)
      fock(i,j)=fc(ij)
      fock(j,i)=fc(ij)
10    continue
      ij=ij+1
      fc(ij)=2.0d0*fock(i,i)
      fock(i,i)=fc(ij)
20    continue
      efreez=ddot(nbasao**2,fock,1,dc,1)
      end if
c.....add one electron hamiltonian
      iho = icorr(ltri)
      call get1(q(1),q(iho),3,dummy)
      jj = 0
      ihoo = 1
      ii = iho
      do 1 isym=1,nirrr
      do 67 i=1,nsymao(isym)
      if(jj.gt.0)
     *call vclr(ho(ihoo),1,jj)
      ihoo = ihoo + jj
      call dcopy(i,q(ii),1,ho(ihoo),1)
      ii = ii+i
67    ihoo = ihoo + i
1     jj = jj + nsymao(isym)
      call corlsr (iho)
      call square(fock,ho,nbasao,nbasao)
      if(nfreez.ne.0)
     1   efreez=efreez+2.0d0*ddot(nbasao**2,fock,1,dc,1)
      ij=0
      do 40 i=1,nbasao
      is=symaos(i)
      do 30 j=istart(is),i
      ij=ij+1
      fc(ij)=fc(ij)+fock(i,j)
      if( dabs(fc(ij)).lt.1.d-10) fc(ij)=0.0d0
30    continue
c...  diagonals times 0.5
      fc(ij)=fc(ij)*0.5d0
40    continue
      call intou(fc,ltri)
      return
      end
      subroutine redun(q,iq,iwrite)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      character*1 lbuf(31)
      dimension q(*),iq(*),dum(1)
      data lbuf /31*' '/
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt('redund',1)
      if (iexc.lt.0) call detci( q(1),iq(1),
     * dum,dum, dum,dum,dum, 1,0)
      call cidmpi
      call mcdump
      ibase = icorr(0)
      nrottu = 0
      ilkp = ind(nact,nact)
      ilkp = icori(ilkp)
      ilkpm = ilkp-1
      lenxi=0
      do 10 i=2,nact
      im=i-1
      do 10 j=1,im
      if (itypea(i).ne.itypea(j)) goto 10
      lenxi=lenxi+1
      iq(ilkpm+ind(i,j))=lenxi
10    continue
      call setsto(lenxi,0,irottu)
      if (iexc.lt.0) goto 999
      if (lenxi.eq.0) goto 150
      lenxi2=lenxi**2
      ixi = icorr(lenxi2)
      ixim = ixi-1
      call vclr(q(ixi),1,lenxi2)
      ivec = icorr(nci)
      ivecm = ivec - 1
      do 20 i=1,nci
20    q(ivecm+i)=( dble(i)+4.73265d0) / ( dble(nci)+3.87453d0)
c      read (5,*) (q(ivecm+i),i=1,nci)
c      call outvec (q(ivec),nci,'ci vector')
c
c...  read 1-electron formulae to form small xi (orb-ci coupling) matrix
c     and add its contrib to redundancy matrix. do in batches if necessa
      lenbuf = icorrm()
      nbuf = lenbuf /lenxi
      npass = (nci-1)/nbuf + 1
      iminm=0
      nbuff=nbuf
      lenbuf=nbuf*lenxi
      ibuf = icorr(-lenbuf)
      ibufm = ibuf - 1
      if (npass.gt.1.and.mcprin) write(iwrite,30) npass
30    format(1x,'orbital-ci part of redundancy matrix in ',
     +       i3,' passes')
      do 110 ipass=1,npass
      if (ipass.eq.npass) nbuff=nci-iminm
      call vclr(q(ibuf),1,lenbuf)
      iblko1=nblkd1+iblkft+1
      call search (iblko1,numft)
      call find (numft)
40    call get (g,nw)
      if (nw.eq.0) goto 90
      call find(numft)
      call unpkft (nw)
      do 80 i=1,nw
      io = iorb(i)
      jo = jorb(i)
      if (io-jo) 50,80,60
50    sign = +1.0d0
      goto 70
60    sign=-1.0d0
70    iii = ii(i)-iminm
      jjj = jj(i)-iminm
      ij = iq(ilkpm+ind(io,jo))
      if (iii.le.nbuff .and. iii.gt.0) q(ibufm+(iii-1)*lenxi+ij)=
     1q(ibufm+(iii-1)*lenxi+ij) + q(ivecm+jjj+iminm)*surd(ival(i))*sign
      if (jjj.le.nbuff .and. jjj.gt.0) q(ibufm+(jjj-1)*lenxi+ij)=
     1q(ibufm+(jjj-1)*lenxi+ij) - q(ivecm+iii+iminm)*surd(ival(i))*sign
80    continue
      goto 40
90    do 100 iii=1,nbuff
      do 100 ij=1,lenxi
      do 100 kl=1,ij
100   q(ixim+(kl-1)*lenxi+ij)=q(ixim+(kl-1)*lenxi+ij)
     1 - q(ibufm+(iii-1)*lenxi+ij) * q(ibufm+(iii-1)*lenxi+kl)
c
110   iminm=iminm+nbuff
      call corlsr (ibuf)
c      call outsqr (q(ixi),lenxi,lenxi,lenxi,'xi contri to xi')
c
c.... get density matrices
      idm = icorr(ne)
      idmm = idm-1
      idmm1 = idmm + ic1d
      call mcdens (q(1),q(ivec),q(idm))
c      call outtri (q(idm),nact**2,'two particle density matrix')
c      call outsqr (q(idm+ic1d),nact,nact,nact,'one particle density mat
c     >ix')
c...  now 1&2 pdm contrib to orb-orb overlap.
      ij=0
      do 130 i=2,nact
      im=i-1
      do 130 j=1,im
      if (itypea(i).ne.itypea(j)) goto 130
      ij=ij+1
      ijkl = ij + ixim
      do 120 k=2,i
      km=k-1
      if (k.eq.i) km=j
      do 120 l=1,km
      if (itypea(k).ne.itypea(l)) goto 120
      q(ijkl) = q(ijkl) + (q(idmm+ind(ilifa(i)+j,ilifa(l)+k) )
     1                  - q(idmm+ind(ilifa(i)+j,ilifa(k)+l) ) )*2.0d0
      if(j.eq.k) q(ijkl) = q(ijkl) - q(idmm1+ilifa(i)+l)
      if(j.eq.l) q(ijkl) = q(ijkl) + q(idmm1+ilifa(i)+k)
      if(i.eq.k) q(ijkl) = q(ijkl) + q(idmm1+ilifa(j)+l)
      if(i.eq.l) q(ijkl) = q(ijkl) - q(idmm1+ilifa(j)+k)
      ijkl = ijkl + lenxi
120   continue
130   continue
      call corlsr (ivec)
c      call outsqr (q(ixi),lenxi,lenxi,lenxi,'redundancy matrix')
      ieig = icorr(lenxi)
      ieigm = ieig-1
      iw = icorr(lenxi)
      call f02abf (q(ixi),lenxi,lenxi,q(ieig),q(ixi),lenxi,q(iw),
     * ifail)
      call outvec (q(ieig),lenxi,'eigenvalues of redundancy matrix')
c      call outsqr (q(ixi),lenxi,lenxi,lenxi,'eigenvectors')
      call corlsr (iw)
c
      call setsto(lenxi,1,irottu)
c
150   nrottu = lenxi
      do 170 i=1,lenxi
      if ( dabs(q(ieigm+i)) .gt. 1.d-6) goto 170
      test = 0.0d0
      kk = 0
      do 160 j=1,lenxi
      if (irottu(j).eq.0 .or.  dabs(q(ixim+(i-1)*lenxi+j)).lt.test)
     1 goto 160
      test =  dabs(q(ixim+(i-1)*lenxi+j))
      kk=j
160   continue
      if (kk.eq.0) call caserr('error in redundancy check.')
      nrottu = nrottu - 1
      irottu(kk)=0
170   continue
999   if (nrottu.eq.0.and.mcprin) write(iwrite,180)
180   format(1x,'no non-redundant active-active orbital rotations')
      if (nrottu.eq.0) goto 230
      im = nact - 1
      if( mcprin) write(iwrite,190) (i,i=1,im)
190   format(1x,'non-redundant active-active orbital rotations:'/
     +           3x,31i3)
      kk=0
      do 220 i=2,nact
      im = i-1
      do 200 j=1,im
      lbuf(j) = ' '
      if (itypea(i).ne.itypea(j)) goto 200
      kk=kk+1
      if (irottu(kk).eq.1)   lbuf(j) = '+'
200   continue
      if( mcprin) write(iwrite,210) i,lbuf
210   format(i3,31(2x,a1))
220   continue
230   nrot = nrottu + nrotit + nrotia + nrotta
      nvar = nrot + nci
      if( mcprin) 
     +write(iwrite,240)nci,nrot,nrotit,nrotia,nrottu,nrotta,nvar
240   format(1x,'number of configurations    =',i7/
     +       1x,'number of orbital rotations =',i7/
     +       5x,'(',i4,' core/active',i5,' core/virtual',i4,
     +               ' active/active',i5,' active/virtual)'/
     +       1x,'total number of variables   =',i7/)
      call corlsr (ibase)
      call accnt(' ',1)
      return
      end
      subroutine lkpgen(iwrite)
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
      dimension id(8)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      iadda=0
      iadd=0
      lentri=0
      ltri=0
      ltrimo=0
      nbasis=nbasao-nfreez
      do 40 is=1,nirrr
      nsymm(is)=nsymao(is)-ifreez(is)
      ltri=ltri+nsymao(is)*(nsymao(is)+1)/2
      ltrimo=ltrimo+nsymm(is)*(nsymm(is)+1)/2
      nprm(is)=0
      ncor(is)=0
      nactt(is)=0
      nsecc(is)=0
      lentca(is)=iadda
      lentrc(is)=iadd
      iads=0
      iad=0
      iada=0
      do 30 isk=1,nirrr
      isl=mults(is,isk)
      if (isl-isk) 10,20,30
10    iada=iada+nsymao(isk)*nsymao(isl)
      iad=iad+nsymm(isk)*nsymm(isl)
      goto 30
20    iad=iad+(nsymm(isk)*(nsymm(isk)+1))/2
      iada=iada+nsymao(isk)*(nsymao(isk)+1)/2
30    iads=iads+nsymm(isk)*nsymm(isl)
      lentri=max(lentri,iada)
      lentr(is)=iad
      lentra(is)=iada
      lensqr(is)=iads
      iadda=iadda+iada
40    iadd=iadd+iad
      do 50 i=1,nbasis
      if (i.le.nprim) nprm(itype(i))=nprm(itype(i))+1
      if (i.le.ncoremc) ncor(itype(i))=ncor(itype(i))+1
      if (i.gt.nprim) nsecc(itype(i))=nsecc(itype(i))+1
      if (i.gt.ncoremc.and.i.le.nprim) nactt(itype(i))=
     +                                nactt(itype(i))+1
50    continue
      nrotit=0
      nrotia=0
      nrotta=0
      do 60 is=1,nirrr
      nrotit = nrotit + ncor(is)*nactt(is)
      nrotia = nrotia + ncor(is)*nsecc(is)
60    nrotta = nrotta + nactt(is)*nsecc(is)
      maxprm=0
      maxbas=0
      do 70 is=1,nirrr
      maxprm=max(maxprm,nprm(is))
70    maxbas=max(maxbas,nsymao(is))
      lenprm=maxprm**2
      lenrec=maxprm*maxbas
      lensq=maxbas**2
      do 80 i=1,8
80    id(i)=istart(i)
      do 90 i=1,nbasis
      isym=itype(i)
      ii=id(isym)
      id(isym)=id(isym)+1
90    iorbsm(ii)=i
      do 100 i=1,nact
      ilifa(i) = (i-1)*nact
100   itypea(i)=itype(i+ncoremc)
      if (btest(iprint,4)) call outive (itypea,nact,'itypea')
c...  look up to 1 elec
      ic1d = ind(nact**2,nact**2)
      ne = ic1d+nact**2
      if (btest(iprint,4)) write(iwrite,110)ne
110   format(/i8,' all active integrals')
      return
      end
      subroutine orbita(iwrite,isymv,ovect,mcprin)
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
      logical logsym,logold,ovect,mcprin
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
      character *8 zrunt, zdd3, zdd4
      character *4 yrunt, ydd
      common/direc/zrunt(50), yrunt(50), ydd(220), zdd3(70) ,zdd4(50)
c
      character*6 acword
      character*8 wordt
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
c
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
      integer orbfrm,symorb,optio,spec,sspesh
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
      character*15 word(9)
      character*2 inwor,char2
      character*1 blank1,cclef(3)
      dimension id(8),isymv(*)
      character*8 worspn(5)
      data blank1,cclef /' ',' ','/','%'/
      data inwor /'in'/
      data word /'doubly occupied',' ','doubly occupied'
     1          ,' ','doubly occupied','unoccupied'
     2          ,      'alpha occupied','beta occupied'
     3          ,       'special '/
      data worspn /'singlet','doublet','triplet','quartet','quintet '/
c
      if(mcprin)then
      write(iwrite,10)
10    format (//
     1 15x,18('*'),' orbital information ',18('*')/
     2 15x,'*',55x,'*'/
     3 15x,'*',5x,'orbitals',11x,'type of orbital',8x,'sym',5x,'*'/
     4 15x,'*',5x,'--------',11x,'---------------',8x,'---',5x,'*')
      write(iwrite,20)
20    format(15x,'*',55x,'*')
      endif
      do 30 i=1,nirrr
30    id(i)=nsymm(i)
      multip=1
      logold=.true.
      imulti=0
      spcsym=1
      fzcdon=0
      cordon=0
      bf=0
      bf2=0
      do 35 i=1,8
35    ifreez(i)=0
      nfreez=0
      ncoremc=0
      nact=0
      norbs=0
      itag=0
 40   itag=itag+1
      wordt=zorb(itag)
c
      if (inporb(wordt,repcnt,clef,cod,sym,logsym).ne.0) goto 110
      if (clef.eq.multrf) imulti=1
      bf1=bf+1+nfreez
      bf2=bf+repcnt+nfreez
      if(ovect) sym = isymv(bf1)
      if(logsym)go to 320
      sym=isymmo(bf1)
      if(ovect) sym = isymv(bf1)
      if(bf1+1.gt.bf2)go to 330
      do 310 junk=bf1+1,bf2
      if(isymmo(junk).ne.sym)call caserr(
     *'invalid symmetry in repeated orbitals')
 310  continue
 330  logold=.false.
      go to 300
 320  if(.not.logold)call caserr(
     *'inconsistent symmetry specification')
 300  if (mcprin) then
       if (fzcdon.eq.0 .and. cod.ne.fzc .and. nfreez.ne.0)
     * write(iwrite,20)
       if (cordon.eq.0.and. cod.ne.cor.and. ncoremc.ne.0) 
     + write(iwrite,20)
      endif
      if (cod.eq.fzv.or.cod.eq.vir) call caserr (
     1'illegal orbital type.')
      acword='active'
      if(cod.eq.fzc) acword='frozen'
      char2=cclef(clef)//blank1
      if (cod.eq.cor)char2=inwor
      if (mcprin) then
       if(repcnt.ne.1) then
        write(iwrite,50) bf1,bf2,char2,acword,word(cod),sym
       else
        write(iwrite,60)bf1,char2,acword,word(cod),sym
       endif
      endif
50    format (15x,'*',i7,' to',i3,8x,a2,a6,1x,a16,i3,6x,'*')
60    format (15x,'*',i10,       11x,a2,a6,1x,a16,i3,6x,'*')
      do 100 junk=1,repcnt
      if (id(sym).gt.0) goto 80
      if (mcprin) write(iwrite,70) (nsymm(ii),ii=1,nirr)
70    format(' numbers of orbitals of each symmetry:',8i4/)
      call caserr('impossible orbital symmetry request.')
80    id(sym)=id(sym)-1
      if(cod.eq.fzc) then
      if(fzcdon.ne.0) call caserr('illegal fzc')
      nfreez=nfreez+1
      if (nfreez.gt.mcfzc) call caserr('no. of fzc not to exceed mcfzc')
      ifreez(sym)=ifreez(sym)+1
      ifzsym(nfreez)=sym
      goto 100
      end if
      bf=bf+1
      itype(bf)=sym
      if(cod.ne.cor)goto 90
      ncoremc=ncoremc+1
      if (cordon.eq.1) call caserr('illegal cor.')
      fzcdon=1
      goto 100
90    norbs=norbs+1
      if (cod.eq.alp) multip=multip+1
      if (cod.eq.bet) multip=multip-1
      if (cod.eq.alp .or. cod.eq.bet) spcsym=mults(spcsym,sym)
      bfkey(norbs)=clef
      ibfcod(norbs)=cod
      fzcdon=1
      cordon=1
100   continue
      go to 40
c
110   continue
      nbasis=nbasao-nfreez
      nact=norbs
      nprim=nact+ncoremc
      nprimp=nprim+1
      nst=ncoremc+1
      nsec=nbasis-nprim
c
      if (mcprin) write(iwrite,20)
      do 130 sym=1,nirrr
      if (id(sym).eq.0) goto 130
      bf1=bf2+1
      bf2=bf1+id(sym)-1
      do 120 i=bf1,bf2
120   itype(i-nfreez)=sym
      if (mcprin) then
       if(id(sym).ne.1)then
        write(iwrite,50) bf1,bf2,inwor,'active',word(uoc),sym
       else
        write(iwrite,60) bf1,inwor,'active',word(uoc),sym
       endif
      endif
130   continue
      if (mcprin) then
        write(iwrite,20)
        write(iwrite,140)
      endif
140   format(15x,57('*')/)
c
c
      nref=1
      if (imulti.eq.0) goto 200
      if (iexc.lt.0) iexc=0
      write(iwrite,190) nref,(codes(ibfcod(i)),
     +                  itype(i+ncoremc),i=1,nact)
150   do 180 bf=1,nact
      if (bfkey(bf).eq.multrf) go to 160
      ibfcod(bf+nref*nact)=ibfcod(bf)
      go to 170
160   itag=itag+1
      wordt=zorb(itag)
      if(inporb(wordt,repcnt,clef,ibfcod(bf+nref*nact),sym,logsym).ne.0)
     * go to 200
      if (sym.ne.itype(bf+ncoremc)) 
     +    call caserr('sym error in reference.')
170   continue
180   continue
      nref=nref+1
      write(iwrite,190) nref,(codes(ibfcod(bf+(nref-1)*nact))
     1               ,itype(bf+ncoremc),bf=1,nact)
190   format(' reference number',i3,5x,10(a3,i1,2x)
     1,30 (/33x,10(a3,i1,2x) ) )
      if ((nref+1)*nact .gt. 885) call caserr (
     1'no room left in dump block to store orbital information')
      goto 150
200   continue
      if (mcprin) then
       if (multip.le.5)  then
        write(iwrite,210) spcsym,worspn(multip)
210     format(' space symmetry =',i2/' spin  symmetry =',a9/)
        else
        write(iwrite,220) spcsym,multip
220     format(' space symmetry =',i2/' spin  symmetry =',i2,'-fold'/)
       endif
      endif
c
c     check for frozen core mos in gradient MCSCF runs
c     these will not work correctly ... use COR instead
c
      jrun = locatc(zrunt,50,zruntp)
      if ((nfreez.gt.0).and.((jrun.ge.4.and.jrun.le.8).or.
     +     jrun.eq.13.or.jrun.eq.14. or.jrun.eq.18)) then
       write(iwrite,230)
 230   format(5x,'***************************************************'/
     +        5x,'**** cannot use FZC mos with MCSCF gradients   ****'/
     +        5x,'**** use the COR descriptor for these orbitals ****'/
     +        5x,'***************************************************'/)
       call caserr('frozen MOs invalidate mcscf gradients')
      endif
      if (nprim.gt.mcprim) then
        write(iwr,'(1x,"MCPRIM = ",i4)')mcprim
        call caserr(
     +  'no. of active orbitals in excess of allowed maximum (MCPRIM)')
      endif
c
      return
      end
      subroutine symmo(vec,mosym)
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
      dimension vec(nbasao,nbasao),mosym(nbasao)
      do 30 i=1,nbasao
      mmax=idamax(nbasao,vec(1,i),1)
      if(mmax.le.0) then
      call caserr('error in assigning symmetries to trial mos')
      endif
 30   mosym(i)=symaos(mmax)
      return
      end
      subroutine orbsym(q,iq,vec)
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
      dimension q(*),iq(*),vec(nbasao,nbasao)
      ix = icorr(nbasao)
      isymmo = icori(nbasao)
      isymm = isymmo-1
      do 30 i=1,nbasao
      mmax=idamax(nbasao,vec(1,i),1)
      if (mmax.le.0) then
      call caserr('error in assigning symmetries to trial mos')
      endif
      iq(isymm+i) = symaos(mmax)
30    continue
      nbm=nbasao-1
      do 60 i=1,nbm
      if(i.gt.nfreez) then
      is=itype(i-nfreez)
      else
      is=ifzsym(i)
      end if
      if (is.eq.iq(isymm+i)) goto 60
      ip=i+1
      do 40 j=ip,nbasao
      if (is.eq.iq(isymm+j)) goto 50
40    continue
      call caserr('symmetry inconsistency in trial mos')
 50   continue
      call dcopy(nbasao,vec(1,j),1,q(ix),1)
      do 55 k=j,ip,-1
      iq(isymm+k)=iq(isymm+k-1)
      call dcopy(nbasao,vec(1,k-1),1,vec(1,k),1)
55    continue
      call dcopy(nbasao,q(ix),1,vec(1,i),1)
      iq(isymm+i)=is
60    continue
      call corlsr (ix)
      return
      end
      subroutine pspace(q,iq)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      common /lsort / icga(31),icgb(31),ina(31),inb(31)
      common/linkmc/irest(7),inp,iblocp
      dimension q(*),iq(*),dum(1)
c
      npread=inp
      nw=npread*nact
      ix=icori(nw)
      nav = lenwrd()
      nw=(nw-1)/nav+1
      call readi(iq(ix),nw*nav,iblocp,num8)
c
      na=0
      nb=0
      if(iexc.ge.0) write(iwr,5)
5     format(/' *** pspace input ignored ***'/)
      if(iexc.ge.0) goto 50
      do 10 i=1,nact
      ii=ibfcod(i)
      if (ii.eq.5.or.ii.eq.7) na=na+1
      if (ii.eq.5.or.ii.eq.8) nb=nb+1
10    continue
      do 111 i=1,nact
111   itypea(i) = itype(i+ncoremc)
      write(iwr,99)
99    format(' p-space configurations given as data:')
50    ijread=1
      ijj=1
      do 1000 ii=1,npread
      if (iexc.ge.0) goto 1000
      iaa=0
      ibb=0
      do 71 i=1,nact
      itest=iq(ijj)
      ijj=ijj+1
      ina(i)=0
      inb(i)=0
      if (itest.ge.10) then
      ina(i)=1
      iaa=iaa+1
      icga(iaa)=i
      end if
      if ((itest/10)*10.ne.itest) then
      inb(i)=1
      ibb=ibb+1
      icgb(ibb)=i
      end if
71    continue
      if (iaa.ne.na.or.ibb.ne.nb) call caserr
     >('illegal number of electrons in pspace')
      call detci (q(1),q(1),dum,dum,dum,dum,dum,10,ivec)
      ipread(ijread) = ivec
      ijread=ijread+1
      write(iwr,79)ivec,(ina(i),inb(i),i=1,nact)
79    format(i8,6x,31(1x,2i1) )
1000  continue
      return
      end
      subroutine mcswap(q,mosym,iwrite)
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
      parameter (pi=3.141592653590d0)
      common/rotat/akk(20),all(20),rota(20),nswap
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
      common /lsort / x(maxorb),y(maxorb)
      dimension q(nbasao,nbasao),mosym(nbasao)
      do 100 iswap=1,nswap
      ak=akk(iswap)
      al=all(iswap)
      k=ak
      l=al
      ksa=(ak-k)*10.1d0
      lsa=(al-l)*10.1d0
      if(ksa.ne.lsa) call caserr('symmetry error on rotate card')
      if(ksa.ne.0) then
      ke=k
      le=l
      k=0
      l=0
      kk=0
      do 5 i=1,nbasao
      if(mosym(i).ne.ksa) goto 5
      kk=kk+1
      if(kk.eq.ke) k=i
      if(kk.eq.le) l=i
      if(k.ne.0.and.l.ne.0) goto 6
5     continue
      else
      ksa=mosym(k)
      lsa=mosym(l)
      if(ksa.ne.lsa) write(iwrite,4)k,ksa,l,lsa
4     format(' ***warning: orbitals on rotate card of different',
     1  ' symmetry:',i4,' (sym=',i1,')',i4,' (sym=',i1,')')
      ke=0
      le=0
      ie=max(k,l)
      do 3 i=1,ie
      if(mosym(i).eq.ksa.and.i.le.k) ke=ke+1
3     if(mosym(i).eq.lsa.and.i.le.l) le=le+1
      end if
6     rot = 90.0d0
      rot=rota(iswap)
      if ( dabs(rot).lt.1d-8) rot = 90.0d0
      if(nbasao.lt.k.or.nbasao.lt.l.or.1.gt.l.or.1.gt.k)
     1call caserr('illegal swapped orbitals')
      write(iwrite,10)rot,k,ke,ksa,l,le,lsa
10    format(' rotation of',f6.1,' degrees between orbitals',i3,' (nr.',
     1 i3,' in symm.',i2,') and',i3,' (nr.',i3,' in symm.',i2,')')
      rot = rot*pi/180.0d0
      coss=  dcos(rot)
      sinn=  dsin(rot)
      call vclr(x,1,nbasao)
      call daxpy(nbasao,coss,q(1,k),1,x,1)
      call daxpy(nbasao,sinn,q(1,l),1,x,1)
      call vclr(y,1,nbasao)
      call daxpy(nbasao,(-sinn),q(1,k),1,y,1)
      call daxpy(nbasao,coss,q(1,l),1,y,1)
      call dcopy(nbasao,x,1,q(1,k),1)
      call dcopy(nbasao,y,1,q(1,l),1)
 100  continue
      return
      end
      subroutine cidmpi
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
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
c
      nblkdu=lensec(mach(18))
      junk = lensec(nci)    *nstate+nblkdu
      call secput(isecd,172,junk,iblkdmc)
      call revind
      iblkc=iblkdmc+nblkdu
      idump=1
      return
      end
      subroutine dmpini
      implicit real*8  (a-h,o-z)
      logical iftran
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
      character *8 param,gitle,zgame
      common /junkc/ param(19),gitle(10)
      character *8 zdate,ztime,zaccno,zanam
      common/jinfo/zdate,ztime,zaccno,zanam
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
      common/blkorbs/evalue(maxorb),occ(mxorb1),
     1                    nbas,newb,ncoll,ivalue,ioccc,ipad
      common/tran/ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     * ctran(mxorb3),iftran,iftri
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
c...   common for harmonic option
      logical  oharm,opharm,odepen
      integer newbas0, newbas1, nsym0, ilifq0, ielimh
      integer newbash,nsymh
      common/harmon/ oharm,opharm,newbas0,newbas1,nsym0(8),
     1               ilifq0(maxorb),ielimh(maxorb),
     2               newbash,nsymh(8),odepen
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
      data zgame/' gamess'/
      lenctr=mach(9)
      lenhed=mach(8)
      lench=29
      nhed=1+lensec(lenhed)+lensec(lenctr)
      len=lensec(newbas1*newbas1)+nhed
      nav = lenwrd()
      param(1)=zanam
      param(2)=zdate
      param(3)=ztime
      param(4)=zgame
      param(5)=' mcscf'
      param(6)=zaccno
      do 70 i=1,10
 70   gitle(i)=ztitle(i)
      nbas=nbasao
      newb=nbas
      ncoll=nbas
      ivalue=-1
      ioccc=-1
      occ(mxorb1)=enext
      call secput (isec,3,len,iblkq)
      call wrtc(param,lench ,iblkq,num3)
      call wrt3s(evalue,lenhed,num3)
      call wrt3is (ilifc,lenctr*nav,num3)
      iblkq=iblkq+nhed
      if(isecn.le.0) return
      if(isecn.eq.isec) then
      iblkn=iblkq
      return
      end if
      call secput(isecn,3,len,iblkn)
      param(5)=' nos'
      call wrtc(param,lench ,iblkn,num3)
      call wrt3s(evalue,lenhed,num3)
      call wrt3is (ilifc,lenctr*nav,num3)
      iblkn=iblkn+nhed
      return
      end
      subroutine cput(c,nv)
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
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
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
      dimension c(*)
      nvec=nv
      if(nv.lt.0) nvec=1
      ii=1
      if(nv.ge.-1) call search(iblkc,num3)
      do 20 i=1,nvec
      call wrt3s(c(ii),nci,num3)
20    ii=ii+nci
      iguess=max(iguess,iabs(nv))
      return
      end
c ******************************************************
c ******************************************************
c             =   mcciopt     =
c ******************************************************
c ******************************************************
      subroutine makep (q,iq,diag,zint,iwrite)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      dimension q(*),iq(*),diag(nci),zint(ne)
      call accnt('makep',2)
      if(nci.lt.nstate+iroot1-1) then
      write(iwrite,5) nci,nstate+iroot1-1
5     format(' more states than configurations specified:',2i3)
      call caserr(' more states than configurations specified')
      endif
      call cid (q(1),zint,diag)
      test=select
      if(select.lt.0) test=1.d20
      if(select.eq.0) test=0.4d0
c.....select configuration with lowest energy
      loop=idmin(nci,diag,1)
      dmin=diag(loop)
      test=test+dmin
c...  check all configurations with same diag h are included
      icon=0
      nplist=0
      do 101 i=1,npread
      ii=ipread(i)
c...  check not already in new list
      do 99 j=1,nplist
      if (iplist(j).eq.ii) goto 101
99    continue
      dd = diag(ii)
      newnp=nplist
      do 102 jj=1,nci
      if ( dabs(diag(jj)-dd).gt.1.0d-6) goto 102
      nplist=nplist+1
      if (nplist.gt.maxp) then
      nplist=newnp
      write(iwrite,1000)
1000  format('not all specified p-space configurations included')
      goto 70
      end if
      iplist(nplist)=jj
      dp(nplist)=dd
102   continue
      icon = icon+1
101   continue
      newnpl=nplist
10    nplist=newnpl
      do 60 i=1,nci
      if(diag(i).gt.test) goto 60
      do 61 j=1,newnpl
      if (iplist(j).eq.i) goto 60
61    continue
      do 20 j=newnpl+1,nplist
20    if(dp(j).gt.diag(i)) goto 30
      j=nplist+1
30    if (j.ge.maxp) goto 60
      nplist=min(maxp,nplist+1)
      if(j.eq.nplist) goto 45
      do 40 k=nplist,j+1,-1
      dp(k)=dp(k-1)
40    iplist(k)=iplist(k-1)
45    dp(j)=diag(i)
      iplist(j)=i
60    continue
70    i2=0
      iee=0
      icon=0
      icone=0
      maxcon=nplist
      if(select.eq.0) maxcon=nstate+iroot1-1
      do 90 i=1,maxcon
      if(i2.eq.nplist) goto 100
      iee=i2
      icone=icon
      i1=i2+1
      do 80 j=i1,nplist
      if( dabs(dp(j)-dp(i1)).gt.1.d-8) goto 90
80    i2=j
90    icon=icon+1
      nplist=i2
100   if(i2.eq.maxp) then
c.....omit last configuration - may not be complete
      nplist=iee
      icon=icone
      end if
      if(icon.lt.nstate+iroot1-1.and.nplist.lt.maxp) then
c.....repeat with larger threshold
      test=test+0.4d0
      goto 10
      end if
      if(icon.lt.nstate+iroot1-1) then
      write(iwrite,110)
110   format(/' not enough p-space configurations can be selected')
      call caserr(
     *' not enough p-space configurations can be selected')
      end if
      do 140 i=1,nstate
140   jcmax(i)=0
      call wrt3 (diag,nci,iblkp,numscr)
      do 141 i=1,nplist
141   dp(i)=dp(i)-dmin
      if (iexc.ge.0) then
      call vfill(1.0d0,vp,1,nplist)
      do 145 i=1,nplist
145   icend(i)=i
      icon=nplist
      if(btest(iprint,14).or.btest(lprint,14)) then
      if(btest(lprint,14)) lprint=ibclr(lprint,14)
      call pricon(iq(1),dp,iplist,nplist,iwrite)
      end if
      else
      call vclr(vp,1,nplist)
      id=0
      icon=0
      do 150 ic=1,nplist
      ia=id+1
      if(ia.gt.nplist) goto 155
      dpp=dp(ia)
      do 160 i=ia,nplist
      if( dabs(dp(i)-dpp).gt.1.d-6) goto 170
160   id=i
170   n=id-ia+1
      jblist=icori(n)
      jclist=icori(n)
      jdlist=icori(n)
      jplist=icori(n)
      ippos=icorr(n**2)
      call spieig(q(1),diag,q(ippos),iq(jblist),iq(jclist),iq(jdlist),
     1    vp(ia),iplist(ia),iq(jplist),icend,icon,n,ia)
      call corlsi(jblist)
150   continue
155   ie=0
      if(btest(iprint,14).or.btest(lprint,14)) then
      if(btest(lprint,14)) lprint=ibclr(lprint,14)
      write(iwrite,175)
175   format(/' p-space configurations'//
     1  '  nr   ndet   energy     occupancy'/)
      do 190 i=1,icon
      ia=ie+1
      ie=icend(i)
190   call pridet(iq(1),i,ie-ia+1,dp(ia),iplist(ia),iwrite)
      write(iwrite,195)
195   format()
      end if
      end if
      call sortpp(iq)
      iblkpq = iposun(numscr)
      call accnt(' ',2)
      return
      end
      subroutine ci2 (q,iq,vec,sigma,zint,eig,idgot,maxdav,accdav,
     +                iwrite)
c...  structure of file 8 :
c      aughes info
c      diis info
c      ci diagonal elements at iblkp
c      pp formula tape or nothing
c      this routine's work vectors at iblkpq
c
      implicit real*8  (a-h,o-z)
      logical incore,debug,btest
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
      character*80 outbuf
      common /lsort / seg(511)
      dimension q(*),iq(*),zint(ne),vec(nci,nstate),sigma(nci,nstate)
     1         ,eig(nstate)
      data m0/0/
      ibase = icorr(0)
      debug=btest(lprint,15).or.btest(iprint,15)
      acci=0
      if (idgot.eq.0) call makep(q(1),iq(1),sigma,zint,iwrite)
      incore = nplist.eq.nci
      call vclr(sigma,1,nci*nstate)
      call cipp (q(1),zint)
c
c....  reduce p space to spin eigenfunction basis
      iuu = icorr(nplist**2)
c...   make diagonal elements numerically small
      ppp = hpp(1,1)
      do 50 i=1,nplist
50    hpp(i,i) = hpp(i,i)-ppp
      nred = 0
      i2=0
      do 120 im=1,icon
      i1=i2+1
      i2=icend(im)
      n2=i2-i1+1
      call vclr(q(iuu+nred*nplist),1,nplist)
      ix1=iuu+(i1-1)+nred*nplist
      call fmove (vp(i1),q(ix1),n2)
c      print*,'n2,snrm2() ',n2,snrm2(n2,q(ix1),1)
      zz = 1.0d0/dnrm2(n2,q(ix1),1)
      call dscal(n2,zz,q(ix1),1)
c      print*,'new degenerate block started'
      nred=nred+1
      nred0=nred
c      print*,'nred = ',nred
c      call outvec (q(iuu+(nred-1)*nplist),nplist,'1st u column')
      do 100 i=i1+1,i2
      call vclr(q(iuu+nred*nplist),1,nplist)
      call mxmb (hpp(i1,i1),1,maxp, q(iuu+i1-1+(nred-1)*nplist),1,0,
     1q(iuu+i1-1+nred*nplist),1,0,    n2,n2,1)
80    do 90 j=nred0,nred
      ix1 = iuu+(j-1)*nplist+i1-1         
      ix2 = iuu+nred*nplist+i1-1
      zz = -ddot(n2,q(ix1),1,q(ix2),1)
      call daxpy(n2,zz,q(ix1),1,q(ix2),1)
90    continue
      ix2 = iuu+nred*nplist+i1-1
      zz = dnrm2(n2,q(ix2),1)
c      print*,'orthog: i,zz: ',i,zz
      if (zz.lt.1.0d-6) goto 110
      call dscal(n2,(1.0d0/zz),q(ix2),1)
      if (zz.lt.0.1d0) goto 80
c      call outvec (q(iuu+(nred)*nplist),nplist,'new u column')
100   nred = nred + 1
110   continue
120   continue
      do 140 i=1,nplist
140   hpp(i,i) = hpp(i,i) + ppp
c..   thus u(j,i); i=old, j=new
c
c...  project down the <p|h|p>
      ndim = nred + nstate*maxdav
      if (incore) ndim = nred
      ihmat = icorr(ndim**2)
      call vclr(q(ihmat),1,ndim**2)
      iwork = icorr(nplist*nred)
      call mxmaa (hpp,1,maxp, q(iuu),1,nplist, q(iwork),1,nplist,
     1           nplist,nplist,nred)
      call mxmaa (q(iuu),nplist,1, q(iwork),1,nplist, q(ihmat),1,ndim,
     1           nred,nplist,nred)
      if(debug) then
      call outsqr (hpp,maxp,nplist,nplist,'hpp')
      call outsqr (q(iuu),nplist,nplist,nred,'uu')
      call outsqr (q(ihmat),ndim,nred,nred,'hmat')
      end if
      call corlsr (iwork)
c
c...  do we have an initial guess for the q space yet?
      if (iguess.eq.0) then
c...  diagonalise <p|h|p>
      ivec = icorr(nred**2)
      ieig = icorr(nred)
      iwork = icorr(nred)
      call f02abf (q(ihmat),ndim,nred,q(ieig),q(ivec),nred,q(iwork),
     * ifail)
      call fmove (q(ieig+iroot1-1),eig,nstate)
      call mxmaa (q(iuu),1,nplist, q(ivec+(iroot1-1)*nred),1,nred,
     > cp,1,maxp, nplist,nred,nstate)
      call corlsr (ivec)
      call pphase(q(ivec),nred)
      call vclr(sigma,1,nci*nstate)
      call cipq (q(1),vec,sigma,zint)
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'vector returned from cipq')
      call vclr(vec,1,nci*nstate)
      do 1046 istate=1,nstate
      do 1046 j=1,nplist
      jj=iplist(j)
1046  sigma(jj,istate)=sigma(jj,istate)-eig(istate)*cp(j,istate)
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'initial residual vector')
      call creatp (vec)
      call dci (vec,sigma,eig,iwrite)
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,1045)m0,dumt,(eig(i)+core,i=1,nstate)
1045  format(i7,f9.2,56x,2f16.8,:,2(/72x,2f16.8))
      endif
      if (.not.incore) iguess=1
      end if
      call anilp(vec)
      if(debug)
     1 call outsqv(cp,maxp,nplist,nstate,'initial p-space vector')
      if (.not.incore.and.debug) then
      call outsqv(vec,nci,nci,nstate,'initial q-space vector')
      end if
c
c....  begin loop over davidson iterations
      do 310 idav=1,maxdav
      nvec = nred
      if (incore) goto 220
      call orthdc (vec,nci,nstate,ndel1)
c...  current non-orthogonalised expansion vector in vec at this point
      lbln = lensec(nci*nstate)
      ibln = iblkpq
      do 160 jdav=1,idav-1
      call rdedx(sigma,nci*nstate,ibln,numscr)
      do 150 i=1,nstate-ndel1
      do 150 j=1,nstate
150   call orth (nci,vec(1,i),sigma(1,j))
160   ibln = ibln + lbln*2
      call orthdc (vec,nci,nstate-ndel1,ndel2)
      ndel=ndel1+ndel2
      if(ndel.ne.0) write(iwrite,165)ndel
165   format(/1x,i1,' expansion vector(s) eliminated'/)
      if(ndel.eq.nstate) call caserr('all expansion vectors deleted')
      if(debug)
     1 call outsqv(vec,nci,nci,nstate,'expansion vector')
      call wrt3 (vec,nci*nstate,ibln,numscr)
      call vclr(sigma,1,nci*nstate)
      call ciqq (q(1),vec,sigma,zint)
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'vector returned from hqq')
      call wrt3s (sigma,nci*nstate,numscr)
c...  now border elements of the hamiltonian
      ioffs=nred+(idav-1)*nstate-1
      do 190 i=1,nstate
      ioffs=ioffs+1
c...  interaction with p space
      do 170 k=1,nplist
      j=iplist(k)
      do 170 l=1,nred
      zz = sigma(j,i)*q(iuu+(l-1)*nplist+k-1)
      q(ihmat+(l-1)+ioffs*ndim) = q(ihmat+(l-1)+ioffs*ndim) + zz
170   q(ihmat+(l-1)*ndim+ioffs) = q(ihmat+(l-1)*ndim+ioffs) + zz
c...  interaction with current q space
      joffs=nred+(idav-1)*nstate-1
      do 180 j=1,i
      joffs=joffs+1
      zz = ddot(nci,vec(1,i),1,sigma(1,j),1)
      q(ihmat+ioffs+joffs*ndim) = zz
180   q(ihmat+ioffs*ndim+joffs) = zz
190   continue
c...  interaction with previous q vectors
      ibln = iblkpq+lbln
      do 210 jdav=1,idav-1
      call rdedx(sigma,nci*nstate,ibln,numscr)
      ioffs=nred+(idav-1)*nstate-1
      do 200 i=1,nstate
      ioffs=ioffs+1
      joffs=nred+(jdav-1)*nstate-1
      do 200 j=1,nstate
      joffs=joffs+1
      zz = ddot(nci,sigma(1,j),1,vec(1,i),1)
      q(ihmat+ioffs+joffs*ndim) = zz
200   q(ihmat+ioffs*ndim+joffs) = zz
210   ibln = ibln + lbln*2
c
c...  diagonalise reduced hamiltonian
      nvec = nred + idav*nstate
220   ivec = icorr(nvec**2)
      ieig = icorr(nvec)
      iwork = icorr(nvec)
      if(debug)
     1 call outsqr (q(ihmat),ndim,nvec,nvec,'p+q hamiltonian')
      call f02abf (q(ihmat),ndim,nvec,q(ieig),q(ivec),nvec,q(iwork),
     * ifail)
      if(ndel.ne.0) call chkvec(q(ieig),q(ivec),nvec,ndel)
      if(debug)
     1 call outvec (q(ieig),nvec,'eigenvalues')
      if(debug)
     1 call outsqr (q(ivec),nvec,nvec,nvec,'eigenvectors of reduced h')
      call fmove (q(ieig+iroot1-1),eig,nstate)
      call corlsr (ieig)
c      print*,'first eigenvalue ',eig(1),eig(2)
c...  p part of eigenvectors
      call mxmaa (q(iuu),1,nplist,q(ivec+(iroot1-1)*nvec),1,nvec,
     > cp,1,maxp, nplist,nred,nstate)
      call pphase(q(ivec),nvec)
      if(debug)
     1 call outsqv(cp,maxp,nplist,nstate,'new p-space vector')
      call vclr(sigma,1,nci*nstate)
      if (.not.incore) then
      call cipq (q(1),vec,sigma,zint)
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'vector returned from hpq')
      end if
      call vclr(vec,1,nci*nstate)
      call creatp(vec)
      if (incore) goto 290
c...  q part of eigenvectors
      call search (iblkpq,numscr)
      istt=nred+(iroot1-1)*nvec
      do 270 jdav=1,idav
      ist = istt
      io=1
      do 240 ipage=1,lbln
      call find (numscr)
      call get (seg,nw)
      ifirst=1
230   ilast=nci-io+1
      call mxmb (seg(ifirst),1,0, q(ivec+ist),1,nvec,
     1 vec(io,1),1,nci, min(nw,ilast),1,nstate)
      if (nw.lt.ilast) goto 240
      io=1
      nw=nw-ilast
      ist=ist+1
      ifirst=ifirst+ilast
      if (nw.gt.0) goto 230
240   io=io+nw
      io=1
      do 260 ipage=1,lbln
      call find (numscr)
      call get (seg,nw)
      ifirst=1
250   ilast=nci-io+1
      if(debug) print *,'ipage=',ipage,' nw=',nw,' ilast=',ilast,
     1    ' nvec=',nvec
      call mxmb (seg(ifirst),1,0, q(ivec+istt),1,nvec,
     1 sigma(io,1),1,nci, min(nw,ilast),1,nstate)
      if (nw.lt.ilast) goto 260
      io=1
      nw=nw-ilast
      istt=istt+1
      ifirst=ifirst+ilast
      if (nw.gt.0) goto 250
260   io=io+nw
270   continue
      if(debug)
     1 call outsqv(vec,nci,nci,nstate,'current ci vector')
c...  hc - ec
      do 280 istate=1,nstate
      call daxpy(nci
     + ,(-eig(istate)),vec(1,istate),1,sigma(1,istate),1)
280   continue
290   continue
      call corlsr (ivec)
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'current residual vector')
      do 1056 istate=1,nstate
      do 1056 j=1,nplist
      jj=iplist(j)
      if ( dabs(sigma(jj,istate)).gt.ciderr) then
      write (outbuf,1057) istate,j,jj,sigma(jj,istate)
      call caserr(outbuf)
      end if
1056  continue
1057  format('non zero p-space gradient; istate,j,jj,sigma:',
     > 3i7,e13.4)
      varc = dnrm2(nci*nstate,sigma,1)
      if(acci.eq.0) acci=dmax1(1.d-5,accdav*varc)
      if (varc.lt.acci .or. idav.ge.maxdav) goto 320
      call dci (vec,sigma,eig,iwrite)
c...  use delta c
      call dcopy(nci*nstate,sigma,1,vec,1)
      if(debug)
     1 call outsqv(vec,nci,nci,nstate,'update of ci vector')
      call anilp(vec)
      ddc = dnrm2(nci*nstate,sigma,1)
      if (mcprin) then
       dumt = seccpu()
       write(iwrite,300)idav,dumt,varc,ddc,(eig(i)+core,i=1,nstate)
      endif
300   format(i7,f9.2,28x,2f14.8,2f16.8,:,2(/72x,2f16.8))
      if (varc.gt.1.0d0) call caserr('excessive gradient in ci')
310    continue
320   call corlsr (ibase)
      if(debug)
     1 call outsqv(vec,nci,nci,nstate,'final ci vector')
      if(debug)
     1 call outsqv(sigma,nci,nci,nstate,'final residual vector')
      return
      end
      subroutine spieig(q,diag,c,nb,nc,list,vp,iplist,jplist,
     +                  icend,ncon,ndet,ia)
      implicit real*8  (a-h,o-z)
      dimension q(*),vp(*),iplist(*),jplist(*),icend(*)
      dimension diag(*),c(ndet,ndet),nb(*),nc(*),list(*)
      dimension dum(1)
c
      icon=0
      do 10 i=1,ndet
10    list(i)=iplist(i)
      call vclr(c,1,ndet**2)
      do 30 i=1,ndet
      if(list(i).eq.0) goto 30
      icon=icon+1
      call detci (q(1),q(1),dum,dum, diag,dum,dum,
     *  2,iplist(i))
      do 20 j=1,ndet
      if(diag(iplist(j)).eq.0) goto 20
      c(j,icon)=(i*0.37d0+0.61d0)*diag(iplist(j))
      list(j)=0
20    continue
30    continue
c
c...  blocking
c
      nd=0
      do 40 i=1,ndet
      nb(i)=0
40    nc(i)=0
      do 100 ic=1,icon
      if(nc(ic).ne.0) goto 100
      ncon=ncon+1
      nc(ic)=ncon
      m=0
      do 50 j=1,ndet
      if(c(j,ic).eq.0) goto 50
      m=m+1
      list(m)=j
      if(nb(j).ne.0) call caserr('blocking error')
      nb(j)=ncon
50    continue
      k=1
60    l=list(k)
      do 80 jc=ic+1,icon
      if(nc(jc).ne.0) goto 80
      if(c(l,jc).eq.0) goto 80
      nc(jc)=ncon
      do 70 j=1,ndet
      if(nb(j).ne.0) goto 70
      if(c(j,jc).eq.0) goto 70
      m=m+1
      list(m)=j
      nb(j)=ncon
70    continue
80    continue
      k=k+1
      if(k.le.m) goto 60
c.....block finished
      do 90 k=1,m
      nd=nd+1
      l=list(k)
      jplist(nd)=iplist(l)
      do 90 jc=1,icon
      if(nc(jc).ne.ncon) goto 90
      vp(nd)=vp(nd)+c(l,jc)
90    continue
      icend(ncon)=nd+ia-1
100   continue
      do 200 i=1,ndet
200   iplist(i)=jplist(i)
      return
      end
      subroutine pridet(iq,ic,ndet,vp,i,iwrite)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1              ,jorb(511),korb(511),lorb(511)
      dimension iq(*),vp(*)
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
c
      icga=icori(nact)-1
      icgb=icori(nact)-1
      iocc=icori(2*nact)-1
      call detocc (i,iq(icga+1),iq(icgb+1))
      nact2=2*nact
      call setsto(nact2,0,iq(iocc+1))
      do 326 j=1,na
      ia=2*iq(icga+j)-1
326   iq(iocc+ia)=1
      do 327 j=1,nb
      ib=2*iq(icgb+j)
327   iq(iocc+ib)=1
      write(iwrite,335) ic,ndet,vp(1),(iq(iocc+j),j=1,nact2)
335   format(1x,i3,i6,f10.5,5x,32(2i1,1x))
      call corlsi(icga+1)
      return
      end
      subroutine pphase(q,n)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      dimension q(n,nstate)
      do 20 i=1,nstate
      if(jcmax(i).eq.0) then
c.....determine largest p-space coeff.
      jcmax(i)=idamax(nplist,cp,1)
      end if
c.....make largest coeff. positive
      j=jcmax(i)
      if(cp(j,i).gt.0) goto 20
      call vneg(cp(1,i),1,cp(1,i),1,nplist)
      call vneg(q(1,i),1,q(1,i),1,n)
20    continue
      return
      end
      subroutine cphase(c,g)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      dimension c(nci,nstate),g(nci,nstate)
      do 20 i=1,nstate
      if(jcmax(i).eq.0) then
c.....determine largest p-space coeff.
      cmax=0
      do 5 jj=1,nplist
      j=iplist(jj)
      if( dabs(c(j,i)).lt.cmax) goto 5
      cmax=dabs(c(j,i))
c.....leave jcmax refer to p-space to be consistent with pphase
      jcmax(i)=jj
5     continue
      end if
c.....make largest coeff. positive
      j=iplist(jcmax(i))
      if(c(j,i).gt.0) goto 20
      call vneg(c(1,i),1,c(1,i),1,nci)
      call vneg(g(1,i),1,g(1,i),1,nci)
20    continue
      return
      end
      subroutine orthdc(c,nci,nstate,ndel)
      implicit real*8  (a-h,o-z)
      dimension c(nci,nstate)
      ndel=0
      do 20 ii=1,nstate
      zz=dnrm2(nci,c(1,ii),1)
      if(zz.lt.1.d-14) call caserr('norm of trial vector too small')
      call dscal(nci,1.0d0/zz,c(1,ii),1)
      do 10 jj=1,ii-1
      zz=-ddot(nci,c(1,jj),1,c(1,ii),1)
      call daxpy(nci,zz,c(1,jj),1,c(1,ii),1)
10    continue
      zz=dnrm2(nci,c(1,ii),1)
      if (zz.gt.1.0d-8) then
      zz = 1.0d0/zz
      else
      ndel=nstate-ii+1
      call vclr(c(1,ii),1,nci*ndel)
      return
      end if
      call dscal(nci,zz,c(1,ii),1)
20    continue
      return
      end
      subroutine chkvec(e,q,n,ndel)
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
      dimension e(*),q(n,*)
c
      n1=n-ndel
      do 100 i=1,nstate
      if( dabs(e(i)).gt.1.d-8) goto 100
c.....following case extremely unlikely and not tested!
      do 20 j=i,i+ndel-1
      do 20 k=1,n1
20    if(dabs(q(k,j)).gt.1.d-7)
     +        call caserr('mixing with deleted vectors')
      call fmove(q(1,i+ndel),q(1,i),n*(nstate-i+1))
      call fmove(e(i+ndel),e(i),nstate-i+1)
      if( dabs(e(i)).lt.1.d-8)call caserr('error in chkvec')
      return
100   continue
      return
      end
      subroutine anilp (c)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      dimension c(nci,nstate)
c..   annihilate p space
      do 10 i=1,nstate
      do 10 jj=1,nplist
      j=iplist(jj)
      cp(jj,i)=c(j,i)
10    c(j,i)=0.0d0
      return
      end
      subroutine sortpp(iq)
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
      integer packft,out
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /intbu2/ iblfx,intpo2,intfi2,intmo2,intnw2,intb2,g2(511)
      dimension iq(*),out(511)
      if (iexc.lt.0) goto 80
      ibase = icorr(0)
      iplis = icori(nci)-1
      call setsto(nci,0,iq(iplis+1))
      do 20 i=1,nplist
20    iq(iplis+iplist(i))=i
      call search (1+iblkft,numft)
      npp=0
      do 40 i=1,nblkd1
      call find(numft)
      call get(g,nw)
      call unpkft(nw)
      do 40 iw=1,nw
      iii = iq(iplis+ii(iw))
      if (iii.eq.0) goto 40
      npp=npp+1
      out(npp) = packft (iii,iii,ival(iw),iorb(iw),jorb(iw),0,0)
c      print31,npp,iii,iii,ival(iw),iorb(iw),jorb(iw),0,0,out(npp)
c30   format(2i4,i3,i6,i6,3i3,z20)
      if (npp.eq.511) call put (out,511,numscr)
      npp = mod(npp,511)
40    continue
      do 50 i=1,nblko1
      call find(numft)
      call get (g,nw)
      call unpkft(nw)
      do 50 iw=1,nw
      iii = iq(iplis+ii(iw))
      jjj = iq(iplis+jj(iw))
      if (iii.eq.0.or.jjj.eq.0) goto 50
      npp=npp+1
      out(npp) = packft(iii,jjj,ival(iw),iorb(iw),jorb(iw),0,0)
c      print31,npp,iii,jjj,ival(iw),iorb(iw),jorb(iw),0,0,out(npp)
      if (npp.eq.511) call put (out,511,numscr)
      npp=mod(npp,511)+1
      out(npp) = packft(jjj,iii,ival(iw),jorb(iw),iorb(iw),0,0)
c      print31,npp,jjj,iii,ival(iw),jorb(iw),iorb(iw),0,0,out(npp)
      if (npp.eq.511) call put (out,511,numscr)
      npp=mod(npp,511)
50    continue
      if (npp.gt.0) call put (out,npp,numscr)
      npp=0
      call put (out,0,numscr)
c
      call find (numft)
      call get (g,nw)
      do 60 i=1,nblkd2
      call find (numft)
      call get (g,nw)
      call unpkft(nw)
      do 60 iw=1,nw
      iii = iq(iplis+ii(iw))
      if (iii.eq.0) goto 60
      npp=npp+1
      out(npp)=packft(iii,iii,ival(iw),iorb(iw),jorb(iw),
     1                korb(iw),lorb(iw))
c      print31,npp,iii,iii,ival(iw),iorb(iw),jorb(iw),korb(iw),lorb(iw)
c     >,out(npp)
      if (npp.eq.511) call put (out,511,numscr)
      npp = mod(npp,511)
60    continue
      do 70 i=1,nblko2
      call find (numft)
      call get (g,nw)
      call unpkft(nw)
      do 70 iw=1,nw
      iii = iq(iplis+ii(iw))
      jjj = iq(iplis+jj(iw))
      if (iii.eq.0.or.jjj.eq.0) goto 70
      npp=npp+1
      out(npp) = packft(iii,jjj,ival(iw),iorb(iw),jorb(iw),
     1                korb(iw),lorb(iw) )
c      print31,npp,iii,jjj,ival(iw),iorb(iw),jorb(iw),korb(iw),lorb(iw)
c     >,out(npp)
      if (npp.eq.511) call put (out,511,numscr)
      npp=mod(npp,511)+1
      out(npp) = packft(jjj,iii,ival(iw),jorb(iw),iorb(iw),
     1                lorb(iw),korb(iw))
c      print31,npp,jjj,iii,ival(iw),jorb(iw),iorb(iw),lorb(iw),korb(iw)
c     >,out(npp)
      if (npp.eq.511) call put (out,511,numscr)
      npp=mod(npp,511)
70    continue
      if (npp.gt.0) call put (out,npp,numscr)
      call put (out,0,numscr)
      call corlsr (ibase)
80    iblfx = 1
      return
      end
      subroutine cipp (q,zint)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),zint(ne),dum(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt ('cipp',2)
      call vclr(hpp,1,maxp**2)
      if (iexc.lt.0) goto 50
      call search (iblkp+lensec(nci),numscr)
10    call find (numscr)
      call get (g,nw)
      if (nw.eq.0) goto 30
      call unpkft(nw)
      do 20 iw=1,nw
c      print*,ii(iw),jj(iw),surd(ival(iw)),zint(ic1d+ilifa(iorb(iw))
c     >+jorb(iw))
20    hpp(ii(iw),jj(iw))=hpp(ii(iw),jj(iw)) + surd(ival(iw)) *
     1 zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      goto 10
30    call find (numscr)
      call get (g,nw)
      if (nw.eq.0) goto 60
      call unpkft(nw)
      do 40 iw=1,nw
c      print*,ii(iw),jj(iw),surd(ival(iw)),zint(ind(ilifa(iorb(iw))
c     >+jorb(iw),ilifa(korb(iw))+lorb(iw)))
40    hpp(ii(iw),jj(iw))=hpp(ii(iw),jj(iw)) + surd(ival(iw)) *
     1 zint(ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw)))
      goto 30
50    call detci (q(1),q(1),zint,dum, dum,dum,dum, 8,0)
60    continue
      call accnt (' ',2)
      return
      end
      subroutine ciqq (q,vec,sigma,zint)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),zint(ne),vec(nci,nstate),sigma(nci,nstate)
      dimension dum(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt ('ciqq',2)
      if (iexc.lt.0) goto 60
      idiag = icorr(nci)
      call vclr(q(idiag),1,nci)
      call rdedx(surd,511,iblkft,numft)
      call find (numft)
      do 10 i=1,nblkd1
      call get(g,nw)
      call find(numft)
      call unpkft(nw)
      do 10 iw=1,nw
10    q(idiag-1+ii(iw)) = q(idiag-1+ii(iw)) + surd(ival(iw))
     1 * zint(ic1d+ilifa(iorb(iw))+iorb(iw))
      do 20 i=1,nblko1
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 20 iw=1,nw
      val=surd(ival(iw)) *
     1 zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      do 20 istate=1,nstate
      sigma(ii(iw),istate)=sigma(ii(iw),istate) + val*vec(jj(iw),istate)
20    sigma(jj(iw),istate)=sigma(jj(iw),istate) + val*vec(ii(iw),istate)
      call get (g,nw)
c
      call find (numft)
      do 30 i=1,nblkd2
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 30 iw=1,nw
30    q(idiag-1+ii(iw)) = q(idiag-1+ii(iw)) + surd(ival(iw))
     1 * zint(ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw)))
      do 40 i=1,nblko2
      call get (g,nw)
      call find (numft)
      do 40 istate=1,nstate
40    call ci1o2 (zint,vec(1,istate),sigma(1,istate),nact,nw)
      call get (g,nw)
c
      do 51 istate=1,nstate
      call vma(vec(1,istate),1,q(idiag),1,sigma(1,istate),1,
     *  sigma(1,istate),1,nci)
51    continue
      call corlsr (idiag)
      call accnt (' ',2)
      return
60    call detci (q(1),q(1),zint,dum, vec,dum,sigma, 4,nstate)
      call accnt (' ',2)
      return
      end
      subroutine cipq (q,vec,sigma,zint)
c...  vec is workspace -- not preserved
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),zint(ne),vec(nci,nstate),sigma(nci,nstate)
      dimension dum(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt ('cipq',2)
      if (iexc.lt.0) goto 40
      call vclr(vec,1,nci*nstate)
      do 49 istate=1,nstate
      do 49 iii=1,nplist
      i=iplist(iii)
49    vec(i,istate) = cp(iii,istate)
      idiag = icorr(nci)
      call vclr(q(idiag),1,nci)
      call rdedx(surd,511,iblkft,numft)
      call find (numft)
      do 10 i=1,nblkd1
      call get(g,nw)
      call find(numft)
      call unpkft(nw)
      do 10 iw=1,nw
10    q(idiag-1+ii(iw)) = q(idiag-1+ii(iw)) + surd(ival(iw))
     1 * zint(ic1d+ilifa(iorb(iw))+iorb(iw))
      do 20 i=1,nblko1
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 20 iw=1,nw
      val=surd(ival(iw)) *
     1 zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      do 20 istate=1,nstate
      sigma(ii(iw),istate)=sigma(ii(iw),istate) + val*vec(jj(iw),istate)
20    sigma(jj(iw),istate)=sigma(jj(iw),istate) + val*vec(ii(iw),istate)
      call get (g,nw)
c
      call find (numft)
      do 31 i=1,nblkd2
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 31 iw=1,nw
31    q(idiag-1+ii(iw)) = q(idiag-1+ii(iw)) + surd(ival(iw))
     1 * zint(ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw)))
      do 30 i=1,nblko2
      call get (g,nw)
      call find (numft)
      do 30 istate=1,nstate
30    call ci1o2 (zint,vec(1,istate),sigma(1,istate),nact,nw)
      call get (g,nw)
      do 51 istate=1,nstate
      do 51 iii=1,nplist
      i=iplist(iii)
51    sigma(i,istate) = sigma(i,istate) + cp(iii,istate)*q(idiag-1+i)
      call accnt (' ',2)
      return
40    continue
      call detci (q(1),q(1),zint,dum, vec,dum,sigma, 9,nstate)
      call accnt (' ',2)
      return
      end
      subroutine dci (c,g,eig,iwrite)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      common /lsort / dseg(511)
      common /mcopt / var,varc,thzr,one,dec
      dimension c(nci,nstate),g(nci,nstate),eig(nstate)
      dimension de(maxpst)
c      print*,'dci: iblkp = ',iblkp
      call search (iblkp,numscr)
      io=0
      do 21 istate=1,nstate
21    de(istate)=0.0d0
      do 11 ipage=1,(nci-1)/511+1
      call find (numscr)
      call get (dseg,nw)
      do 10 istate=1,nstate
      ee=eig(istate)-cishft
      do 10 i=1,nw
      de(istate)=de(istate)+g(io+i,istate)**2 / (ee-dseg(i))
10    g(io+i,istate) = c(io+i,istate) + g(io+i,istate) / (ee-dseg(i))
11    io=io+511
      dec=0.0d0
      do 22 istate=1,nstate
22    dec = dec + 2.0d0*de(istate)
c      print*,'dci dec = ',dec
      call orthci(g,nci,nstate,iwrite)
      call vsub(g,1,c,1,g,1,nci*nstate)
      call vadd(g,1,c,1,c,1,nci*nstate)
      return
      end
      subroutine orthci(c,nci,nstate,iwrite)
      implicit real*8  (a-h,o-z)
      dimension c(nci,nstate)
      do 20 ii=1,nstate
      do 10 jj=1,ii-1
       zz=-ddot(nci,c(1,jj),1,c(1,ii),1)
      call daxpy(nci,zz,c(1,jj),1,c(1,ii),1)
10    continue
      zz=dnrm2(nci,c(1,ii),1)
      if (zz.gt.1.0d-12) then
      zz = 1.0d0/zz
      else
      write(iwrite,1000)ii
1000  format('****  singularity in schmidt orthogonaliser'
     *      ,' istate = ',i3)
      call caserr(
     *'singularity in schmidt orthogonaliser')
      end if
      call dscal(nci,zz,c(1,ii),1)
20    continue
      return
      end
      subroutine cid (q,zint,diag)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
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
      dimension q(*),zint(ne),diag(nci)
      dimension dum(1)
      dimension occ(31),eps(31)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call vclr(diag,1,nci)
c
      if (iexc.lt.0) goto 30
c
      call rdedx(surd,511,iblkft,numft)
      if (lto(1)) goto 50
      call find(numft)
      do 10 i=1,nblkd1
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 10 iw=1,nw
10    diag(ii(iw))=diag(ii(iw)) + surd(ival(iw))
     1  *zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      call get (g,nw)
c
      call search (nblkd1+nblko1+iblkft+2,numft)
      call find (numft)
      do 20 i=1,nblkd2
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 20 iw=1,nw
20    diag(ii(iw)) = diag(ii(iw)) +
     1 surd(ival(iw))*zint(ind(ilifa(iorb(iw))+jorb(iw),
     2  ilifa(korb(iw))+lorb(iw)))
      call get (g,nw)
c
      goto 40
c
c...  moller-plesset denominators
50    do 60 i=1,nact
      occ(i) = 0.0d0
      iiii = ibfcod(i)
      if (iiii.eq.5) occ(i) = 2.0d0
      if (iiii.eq.7 .or. iiii.eq.8) occ(i) = 1.0d0
60    continue
c      call outvec (occ,nact,'reference occupations')
      do 70 i=1,nact
      eps(i) = zint(ic1d+ilifa(i)+i)
      do 70 j=1,nact
      iijj = ind(ilifa(i)+i,ilifa(j)+j)
      ijij = ind(ilifa(i)+j,ilifa(i)+j)
70    eps(i) = eps(i) + occ(j)*(zint(iijj)-0.5d0*zint(ijij))
     >       + occ(j)*(2.0d0-occ(j))*zint(ijij)
c...  set near equal values the same for high symmetry problems
c      call outvec (eps,nact,'eps before smoothing')
      do 75 i=1,nact
      do 75 j=1,i-1
      if ( dabs(eps(i)-eps(j)).lt.0.1d0) eps(j)=eps(i)
75    continue
c      call outvec (eps,nact,'fock eigenvalues')
      call find(numft)
      do 80 i=1,nblkd1
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 80 iw=1,nw
80    diag(ii(iw)) = diag(ii(iw)) + eps(iorb(iw))*surd(ival(iw))
      call get (g,nw)
      goto 40
c
30    call detci (q(1),q(1),zint,dum, diag,dum,dum, 3,1)
c      call outvec (diag,nci,'ci diagonal elements')
c
40    return
      end
      subroutine ci0 (q,zint,diag)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      common /intbuf/ intpos,intfil,intmod
      integer  iblf,iblf1,iblf2,iword2,jad,kad,lad,lj,lk
      integer  iadr,iadw,ifinit
      common /mcff/ iblf,iblf1,iblf2,iword2,jad(mcprim*(mcprim+1)/2),
     +              kad(mcprim*(mcprim+1)/2),lad,lj(8),lk(8),iadr,
     +              iadw,ifinit
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
      integer ijpos(8),ikpos(8)
      dimension q(*),zint(ne),diag(nci)
      dimension dum(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call accnt('hesini',1)
      call accnt('ci0',2)
      intfil=num6
      iblf  =iblk6
      intpos=0
c
c...  load integrals
      ibase = icorr(lensq)
      do 10 isym=1,nirrr
      ijpos(isym) = icorr(maxbas**2)
10    ikpos(isym) = icorr(maxbas**2)
      do 20 i=1,nprim
      do 20 j=1,i
20    call loadjk (q(1),ijpos,ikpos,i,j,zint)
      call loadjk (q(1),ijpos,ikpos,0,0,zint)
      call corlsr (ibase)
c
c
      call vclr(diag,1,nci)
c
      if (iexc.lt.0) goto 999
c
      call rdedx(surd,511,iblkft,numft)
      call find(numft)
      do 30 i=1,nblkd1
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 30 iw=1,nw
30    diag(ii(iw))=diag(ii(iw)) + surd(ival(iw))
     1  *zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      call get (g,nw)
c
      call search (nblkd1+nblko1+iblkft+2,numft)
      call find (numft)
      do 40 i=1,nblkd2
      call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 40 iw=1,nw
40    diag(ii(iw)) = diag(ii(iw)) +
     1 surd(ival(iw))*zint(ind(ilifa(iorb(iw))+jorb(iw),
     2  ilifa(korb(iw))+lorb(iw)))
      call get (g,nw)
c
      goto 1099
c
999   call detci (q(1),q(1),zint,dum, diag,dum,dum, 3,1)
c      call outvec (diag,nci,'ci diagonal elements')
c
1099  call accnt(' ',1)
      return
      end
      subroutine ci1 (q,vec,sigma,diag,zint)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),zint(ne),vec(nci),diag(nci),sigma(nci)
      dimension dum(1)
      call accnt('ci sigma',2)
      if (iexc.lt.0) goto 999
      call vmul(diag,1,vec,1,sigma,1,nci)
c
      if (nblko1.ne.0) then
      junk=nblkd1+iblkft+1
      call search (junk,numft)
      call find(numft)
      do 20 i=1,nblko1
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 20 iw=1,nw
      val=surd(ival(iw)) *
     1 zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      sigma(ii(iw))=sigma(ii(iw)) + val*vec(jj(iw))
20    sigma(jj(iw))=sigma(jj(iw)) + val*vec(ii(iw))
      call get (g,nw)
      end if
c
      call search (nblkd1+nblko1+nblkd2+2+iblkft,numft)
      call find (numft)
      do 30 i=1,nblko2
      call get (g,nw)
      call find (numft)
      call ci1o2 (zint,vec,sigma,nact,nw)
30    continue
      call get (g,nw)
      goto 1099
c
999   continue
      call vclr(sigma,1,nci)
      call detci (q(1),q(1),zint,dum, vec,dum,sigma, 4,1)
c
1099  call accnt(' ',2)
      return
      end
      function ctracem(a)
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
c
      integer ioffd,nldd,ioffq,nlqq,ioffu,jadrs,kadrs,locc,nt
      integer jadr,kadr,ntdg,ntqg,nlu,lenj,lenk,nmax,numax,nop
      integer iorb,itri,ifc,igrd,nrotti,iorbs
c
      common/scra /ioffd(8,8),nldd(8),ioffq(8,8),nlqq(8),ioffu(8),
     +  jadrs(8),kadrs(8),locc(8),nt(8),jadr(mcprim,mcprim),
     +  kadr(mcprim,mcprim),ntdg,ntqg,nlu,lenj,lenk,nmax,numax,
     +  nop,iorb(mcprim),itri(mcprim),ifc,igrd(mcprim*mcprim),
     +  nrotti,iorbs(mcprim,8)
      dimension a(*)
c
      ctracem=0.0d0
      do 50 is=1,nirrr
      nc=ncor(is)
      if(nc.eq.0) goto 50
      ii=ioffd(1,is)
      do 10 i=1,nc
      ii=ii+i
10    ctracem=ctracem+a(ii)
50    continue
      return
      end
      subroutine mcdens (q,vec,gam)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),gam(ne),vec(nci)
      dimension dum(1)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call accnt('density ',2)
      call vclr(gam,1,ne)
      if (iexc.lt.0) goto 120
      call rdedx(surd,511,iblkft,numft)
      iblk=0
      call find (numft)
10    call get(g,nw)
      call find (numft)
      call unpkft(nw)
      do 20 iw=1,nw
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
20    gam(intad)=gam(intad)+vec(ii(iw))*vec(ii(iw))*surd(ival(iw))
      iblk = iblk + 1
      if (iblk.lt.nblkd1) goto 10
30    call get (g,nw)
      if (nw.eq.0) goto 50
      call find (numft)
      call unpkft(nw)
      do 40 iw=1,nw
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
      gam(intad)=gam(intad)+vec(ii(iw))*vec(jj(iw))*surd(ival(iw))
      intad = ic1d + ilifa(jorb(iw)) + iorb(iw)
40    gam(intad)=gam(intad)+vec(ii(iw))*vec(jj(iw))*surd(ival(iw))
      goto 30
50    continue
c
      iblk=0
      call find (numft)
60    call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 70 iw=1,nw
      intad = ind (ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw) )
c      write (6,*) iorb(iw),jorb(iw),korb(iw),lorb(iw),ilifa(iorb(iw))
c     >,ilifa(korb(iw)),intad
      gam(intad)=gam(intad)+vec(ii(iw))*vec(jj(iw))*surd(ival(iw))
70    continue
      iblk = iblk + 1
      if (iblk.lt.nblkd2) goto 60
c###
      call dscal(ic1d,0.5d0,gam,1)
c###
80    call get (g,nw)
      if (nw.eq.0) goto 90
      call find (numft)
      call denso2 (gam,vec,nw)
      goto 80
90    continue
c###
c...  symmetrise across 2-particle
      do 100 i=1,nact
      do 100 j=1,nact
      lmax=nact
      do 100 k=1,i
      if (i.eq.k) lmax=j
      do 100 l=1,lmax
      i1 = ind(ilifa(i)+j,ilifa(k)+l)
      i2 = ind(ilifa(j)+i,ilifa(l)+k)
      val = 0.5d0 * (gam(i1)+gam(i2))
      gam(i1) = val
100   gam(i2) = val
      call dscal(ic1d,2.0d0,gam,1)
c###
c..   double up ilil / iiii
      junk=nact**2
      intad = 0
      do 110 il=1,junk
      intad = intad + il
110   gam(intad)=gam(intad)+gam(intad)
      goto 130
c
120   call detci (q(1),q(1),dum,gam, vec,dum,dum, 6,1)
c
130   continue
c      call outvec (gam,ne,'2-particle density matrix')
      call accnt('other   ',2)
      return
      end
      subroutine densav (q,vec,gam)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),gam(ne),vec(nci,nstate)
      dimension dum(1)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      ssumm=0.0d0
      if(nstate.gt.1)
     * ssumm=dsum(nstate-1,weight,1)
      if (ssumm.gt.1.0d-10) goto 10
      call mcdens (q(1),vec(1,nstate),gam)
      return
10    continue
      call accnt('density ',2)
      call vclr(gam,1,ne)
      if (iexc.lt.0) goto 140
      call rdedx(surd,511,iblkft,numft)
      iblk=0
      call find (numft)
20    call get(g,nw)
      call find (numft)
      call unpkft(nw)
      do 30 iw=1,nw
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
      do 30 istate=1,nstate
30    gam(intad)=gam(intad)+vec(ii(iw),istate)*vec(ii(iw),istate)*
     1surd(ival(iw))*weight(istate)
      iblk = iblk + 1
      if (iblk.lt.nblkd1) goto 20
40    call get (g,nw)
      if (nw.eq.0) goto 60
      call find (numft)
      call unpkft(nw)
      do 50 iw=1,nw
      do 50 istate=1,nstate
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
      gam(intad)=gam(intad)+vec(ii(iw),istate)*vec(jj(iw),istate)*
     1surd(ival(iw))*weight(istate)
      intad = ic1d + ilifa(jorb(iw)) + iorb(iw)
50    gam(intad)=gam(intad)+vec(ii(iw),istate)*vec(jj(iw),istate)*
     1surd(ival(iw))*weight(istate)
      goto 40
60    continue
c
      iblk=0
      call find (numft)
70    call get (g,nw)
      call find (numft)
      call unpkft(nw)
      do 80 iw=1,nw
      intad = ind (ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw) )
      do 80 istate=1,nstate
      gam(intad)=gam(intad)+vec(ii(iw),istate)*vec(jj(iw),istate)*
     1surd(ival(iw))*weight(istate)
80    continue
      iblk = iblk + 1
      if (iblk.lt.nblkd2) goto 70
c###
      call dscal(ic1d,0.5d0,gam,1)
c###
90    call get (g,nw)
      if (nw.eq.0) goto 110
      call find (numft)
      do 100 istate=1,nstate
      if (weight(istate).le.1.0d-15) goto 100
      call dscal(ic1d,1.0d0/weight(istate),gam,1)
      call denso2 (gam,vec(1,istate),nw)
      call dscal(ic1d,weight(istate),gam,1)
100   continue
      goto 90
110   continue
c###
c...  symmetrise across 2-particle
      do 120 i=1,nact
      do 120 j=1,nact
      lmax=nact
      do 120 k=1,i
      if (i.eq.k) lmax=j
      do 120 l=1,lmax
      i1 = ind(ilifa(i)+j,ilifa(k)+l)
      i2 = ind(ilifa(j)+i,ilifa(l)+k)
      val = 0.5d0 * (gam(i1)+gam(i2))
      gam(i1) = val
120   gam(i2) = val
      call dscal(ic1d,2.0d0,gam,1)
c###
c..   double up ilil / iiii
      junk=nact**2
      intad = 0
      do 130 il=1,junk
      intad = intad + il
130   gam(intad)=gam(intad)+gam(intad)
      goto 150
c
140   call detci (q(1),q(1),dum,gam, vec,dum,dum, 6,nstate)
c
150   continue
c      call outvec (gam,ne,'2-particle density matrix')
      call accnt('other   ',2)
      return
      end
      subroutine denst1 (q,vec1,vec2,gam)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),gam(ne),vec1(nci),vec2(nci)
      dimension dum(1)
c
      call accnt('density ',2)
      call vclr(gam,1,ne)
      if (iexc.lt.0) goto 60
      call rdedx(surd,511,iblkft,numft)
      iblk=0
      call find (numft)
10    call get(g,nw)
      call find (numft)
      call unpkft(nw)
      do 20 iw=1,nw
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
20    gam(intad)=gam(intad)+vec1(ii(iw))*vec2(ii(iw))*surd(ival(iw))
      iblk = iblk + 1
      if (iblk.lt.nblkd1) goto 10
30    call get (g,nw)
      if (nw.eq.0) goto 50
      call find (numft)
      call unpkft(nw)
      do 40 iw=1,nw
      intad = ic1d + ilifa(iorb(iw)) + jorb(iw)
      gam(intad)=gam(intad)+vec1(ii(iw))*vec2(jj(iw))*surd(ival(iw))
      intad = ic1d + ilifa(jorb(iw)) + iorb(iw)
40    gam(intad)=gam(intad)+vec2(ii(iw))*vec1(jj(iw))*surd(ival(iw))
      goto 30
50    continue
      goto 70
c
60    call detci (q(1),q(1),dum,gam,vec1,vec2,dum, 7,1)
c
70    continue
      call accnt('other   ',2)
      return
      end
      subroutine detci (q,iq,z1,z3, v1,v2,v3, ifunct,nvec)
c
c...  basic calling routine for determinant full ci
c...  z1,   z3 integral or density matrices in multi format
c...  v1,v2,v3 ci vectors - standard order civec vec sigma diag
c...  possible function codes ifunct:
c     1  --  only return length of ci in /multic/ nci
c     2  --  spin adapt v1
c     3  --  diagonal elements v1
c     4  --  v3 = h(z1).v1 for nvec vectors
c     5  --  v3 = h(z1).v1 + h(z2).v2; z3 = density(v1,v2)
c     6  --  z3 = density(v1,v1) for nvec vectors;
c                 if nvec.ne.1, use weighting
c     7  --  z3 = density(v1,v2)  1-particle only
c     8  --  determine p space hamiltonian h(z1) --> hpp
c     9  --  determine action of h(z1) on nvec p space vectors
c            -> v3; v1 scratch vector, destroyed on exit
c    10  --  given strings in icga,icgb, return det address in nvec
c    11  --  v3 = s**2 v1 for nvec vectors
c
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
      common /lsort / icga(31),icgb(31)
      dimension q(1),iq(1),z1(ne),z3(ne),v1(1),v2(1),v3(1)
c
      ibase = icorr(0)
c.... obtain numbers of electrons and space symmetry
      na=0
      nb=0
      isss=1
      do 10 i=1,nact
      ii=ibfcod(i)
      if (ii.eq.5.or.ii.eq.7) na=na+1
      if (ii.eq.5.or.ii.eq.8) nb=nb+1
      if (ii.eq.7.or.ii.eq.8) isss=mults(isss,itypea(i))
10    continue
c
c...  set up addressing arrays
      mm = (nact*(nact+1))/2
      nstraa = ibinom (nact,na)
      nstrbb = ibinom (nact,nb)
      ic     = icori (nstraa)
      itz    = icori (nstrbb)
      inter  = icori (na*nact*2)
      ic2e   = icori (mm)
      ic3e   = icori (mm)
      call setsto(na*nact*2,0,iq(inter))
      call detcil (iq(ic),iq(itz),iq(inter),iq(ic2e),iq(ic3e))
      nablk = (maxaa/2)*2+1
      if (mcmach('vector').le.0) nablk=1
      nbblk = (maxbb/2)*2+1
      izint1=0
c     izint2=0
      igam = 0
c
      goto (20,50,60,70,80,90,100,110,120,121,131),ifunct
 20   ndum = ncsf(nact,na,nb,isss,itypea,mults)
      if (mcprin) 
     +    write(iwr,30) nci,ndum,nstraa*nstrbb
30    format(1x,'number of determinants in ci  =',i8/
     +       1x,'equivalent number of csfs     =',i8/
     +       1x,'number of intermediate states =',i8)
      i = icorr(0)-ibase+nint1+nint2
      j = i + (maxrpl*3+1)*(nablk+nbblk)+2*maxpar*nablk*nbblk
      nabl=min(nablk,63)
      nbbl=min(nbblk,63)
      k = i + (maxrpl*3+1)*(nabl+nbbl)+2*maxpar*nabl*nbbl
      if (mcprin) write(iwr,40)i,j,k
40    format(1x,'store required for ci: ',i8,' (minimum)'/
     +       1x,'                       ',i8,' (maximum)'/
     +       1x,'                       ',i8,
     +          ' (for maximum vector length=63)'/)
      goto 170
c
50    call spinad (v1,nvec,iq(ic),iq(itz),iq(inter))
      goto 170
c
60    iz = icorr(nact*maxbb)
      izb = icorr(nact*maxbb)
      if = icorr(maxbb)
      izint = icorr(nint2+nint1)
      call detint (z1,q(izint),iq(ic2e),.false.)
      call diags (v1,q(izint),iq(ic2e),iq(ic3e),q(izb),q(iz),q(if))
      goto 170
c
70    nvec1 = nvec
      nvec2 = 0
      ndens = 0
      nham = nvec
      izint1 = icorr(nint1+nint2)
      call detint (z1,q(izint1),iq(ic2e),.true.)
      goto 130
c
80    continue
      goto 130
c
90    nvec1 = nvec
      nvec2 = 0
      ndens = 1
      nham = 0
      igam = icorr(nint1+nint2)
      goto 130
c
100   nvec1 = 1
      nvec2 = 1
      ndens = -1
      nham = 0
      igam = icorr(nint1+nint2)
      goto 130
c
110   izint = icorr(nint1+nint2)
      ihh = icorr(nact**2)
      icoul = icorr(nact**2)
      iexch = icorr(nact**2)
      call detint (z1,q(izint),iq(ic2e),.false.)
      call detpp (q(izint),iq(ic2e),iq(ic3e),iq(ic),q(ihh)
     1           ,q(icoul),q(iexch) )
      goto 170
c
120   izint = icorr(nint1+nint2)
      iwa = icori(maxrpl*3)
      iwb = icori(maxrpl*3)
      call detint (z1,q(izint),iq(ic2e),.false.)
      call detpq (v1,v3,q(izint),iq(ic2e),iq(ic3e),iq(ic),iq(itz),nvec,
     1            iq(iwa),iq(iwb),iq(inter) )
      goto 170
c
121   nvec = iq(ic-1+istrad(ipar,iq(inter),icga,na))
     >    + iq(itz-1+istrad(ipar,iq(inter+na*nact),icgb,nb))
      goto 170
c
130   mmax = icorrm()-1
140   istore = (maxrpl*3+1)*(nablk+nbblk) + 2*maxpar*nablk*nbblk
      if (istore.le.mmax) goto 160
      if (nablk.gt.nbblk) goto 150
      nbblk = nbblk - 2
      goto 140
150   nablk = nablk - 2
      goto 140
160   nwa = icori(nablk)
      nwb = icori(nbblk)
      iwa = icori(maxrpl*3*nablk)
      iwb = icori(maxrpl*3*nbblk)
      id = icorr(maxpar*nablk*nbblk)
      ie = icorr(maxpar*nablk*nbblk)
      call cislow (v1,v2,v3,q(izint1),q(igam),
     1 nham,nvec1,nvec2,ndens,
     2 iq(ic3e),iq(ic),iq(itz),iq(inter),iq(iwa),iq(iwb),q(id),q(ie),
     3 nablk,nbblk,iq(nwa),iq(nwb) )
      if (ndens.ne.0) call detgam (q(igam),z3,iq(ic2e),iq(ic3e))
      goto 170
c
131   maxpar = nact**2
      mmax = icorrm()-1
141   istore = (maxrpl*3+1)*(nablk+nbblk) + maxpar*nablk*nbblk
      if (istore.le.mmax) goto 161
      if (nablk.gt.nbblk) goto 151
      nbblk = nbblk - 2
      goto 141
151   nablk = nablk - 2
      goto 141
161   nwa = icori(nablk)
      nwb = icori(nbblk)
      iwa = icori(maxrpl*3*nablk)
      iwb = icori(maxrpl*3*nbblk)
      id = icorr(maxpar*nablk*nbblk)
      call ssq (v1,v3,nvec,iq(ic),iq(itz),iq(inter),
     > iq(iwa),iq(iwb),q(id),nablk,nbblk,iq(nwa),iq(nwb) )
      goto 170
c
170   continue
      call corlsr (ibase)
      return
      end
      function ncsf (nact,na,nb,isss,itypea,mult)
      integer icga(31),itypea(nact),mult(8,8),ib(33)
      numb=0
      do 60 nopsh=na-nb,na+nb,2
      ndoc=(na+nb-nopsh)/2
      nuoc=nact-nopsh-ndoc
      if (nuoc.lt.0) goto 70
      if (nopsh.eq.0) then
      nsing = 1
      else
      nsing=0
      mt=0
10    call string (mt,nact,nopsh,icga,itypea,jmta,iua,mult)
      if (mt.eq.0) goto 20
      if (jmta.eq.isss) nsing=nsing+1
      goto 10
20    continue
      end if
c...  scan branching diagram for number of spin eigenfunctions
      do 30 j=2,nopsh+2
30    ib(j)=0
      ib(1)=1
      do 50 i=1,nopsh
      do 40 j=i+1,2,-2
40    ib(j)=ib(j-1)+ib(j+1)
50    ib(1)=ib(2)
      nbrnch = ib(na-nb+1)
      numb = numb + nsing*nbrnch*ibinom(nact-nopsh,ndoc)
60    continue
70    ncsf=numb
      return
      end
      subroutine detcil (ic,itz,inter,ic2e,ic3e)
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
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
      dimension ic(nstraa),itz(nstrbb),inter(na,nact,2),ic2e(1),ic3e(1)
      dimension icg(31)
c
      do 10 i=1,nact+1
10    icf(i) = (i*(i-1))/2
c
      nab=na
      do 50 kk=1,2
      if (kk.eq.2) nab=nb
      do 30 k=1,nab
      inter(k,k,kk)=0
      do 20 l=k,nact-1+k-nab
20    inter(k,l+1,kk) = inter(k,l,kk) + ibinom(nact-l,nab-k)
30    continue
      do 40 i=1,nab-1
      do 40 k=1,nact-1
40    inter(i,k,kk) = inter(i,k,kk) - inter(i+1,k+1,kk)
50    continue
c
      nint2=0
      do 80 iz=1,nirrr
      intoff(iz) = nint2
      kt=0
      ij=0
      do 60 i=1,nact
      do 60 j=1,i
      ij=ij+1
      if (iz.ne.mults(itypea(i),itypea(j))) goto 60
      kt=kt+1
      ic3e(ij)=kt
60    continue
      ij=0
      do 70 i=1,nact
      do 70 j=1,i
      ij=ij+1
      if (iz.ne.mults(itypea(i),itypea(j))) goto 70
      ic2e(ij) = nint2
      nint2 = nint2 + kt
70    continue
      npair(iz) = kt
80    continue
      nint1 = npair(1)
c...sets up addressing of ci vector
c... loop over beta string symmetries
      nci=0
      do 130 iz=1,nirrr
c... first do beta string address offsets -- itz
      kt=0
      mt=0
90    call string (mt,nact,nb,icg,itypea,jmt,iu,mults)
      if (mt.eq.0) goto 100
      if (jmt.ne.iz.and.nb.gt.0) goto 90
      kt=kt+1
      itz(mt)=kt
      if (nb.gt.0) goto 90
100   if (nb.gt.0.or.iz.eq.1) nstrb(iz) = kt
c...  now loop over valid alpha strings
      iaz = mults(isss,iz)
      mt=0
      nstra(iaz)=0
110   call string (mt,nact,na,icg,itypea,jmt,iu,mults)
      if (mt.eq.0) goto 120
      if (jmt.ne.iaz) goto 110
      ic(mt) = nci
      nstra(iaz) = nstra(iaz)+1
      nci = nci + kt
      goto 110
120   continue
130   continue
c
      maxaa=0
      maxbb=0
      do 140 i=1,nirrr
      maxaa=max(maxaa,nstra(i))
140   maxbb=max(maxbb,nstrb(i))
c
c...  max no of singles from a string -- inefficient if symmetry
      maxrpl = max(na*(nact+1-na),nb*(nact+1-nb))
c... max no of orbital pairs
      maxpar=0
      do 150 i=1,nirrr
150   maxpar = max(maxpar,npair(i))
c
      return
      end
      function ibinom (m,n)
      implicit real*8  (a-h,o-z)
      top = m
      bot = n
      binom = 1.0d0
      do 10 i=1,n
      binom = binom * top/bot
      top = top - 1.0d0
10    bot = bot - 1.0d0
      ibinom = binom + 0.1d0
      return
      end
      subroutine string (mt,m,n,icg,itype,jmt,iu,mult)
c... generate a new alpha or beta string
c...  mt,icg,iu should be preserved between calls
c
      implicit real*8  (a-h,o-z)
      dimension icg(n),itype(m),mult(8,8)
c
      if (n.eq.0) then
        if (mt .gt. 0) go to 60
        mt = 1
        jmt = 1
        iu = 0
        return
      end if
c
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
      subroutine onel (icg,isymr,iw,nw,n,ic3e,itz,inter)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
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
      dimension icg(n),iw(maxrpl,3),ic3e(1),itz(1),inter(na,nact)
      dimension ia(32)
      do 10 i=1,n
10    ia(i)=icg(i)
      ia1=ia(1)
      ia(n+1)=0
      iz=1
      nw=0
      do 40 i=1,n
      isymj=mults(isymr,itypea(ia1))
      kk=1
      iz2=iz
      jj=1
      do 20 j=1,n
20    jj=jj+inter(j,ia(j))
      do 30 j1=1,nact
      ii=jj-inter(kk,ia(kk))+inter(kk,j1)
      if (j1.eq.ia(kk+1)) then
       kk=kk+1
       jj=ii
       iz2=-iz2
      else if (itypea(j1).eq.isymj) then
       nw=nw+1
       iw(nw,1) = itz(ii)
       iw(nw,2) = iz2
       iw(nw,3) = ic3e(icf(max(ia1,j1))+min(ia1,j1))
      end if
30    continue
      iz = -iz
      ia1 = ia(i+1)
      ia(i+1)=ia(1)
40    ia(1)=ia1
      return
      end
      subroutine cislow (c1,c2,s,zint1,gam,
     1 nham,nvec1,nvec2,ndens,
     2 ic3e,ic,itz,inter,iwa,iwb,d,e,nablk,nbblk,
     3 nwa,nwb)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension c1(1),c2(1),s(1),zint1(1),gam(1)
     1         ,ic3e(1),ic(1),itz(1),inter(na,nact,2),d(1),e(1)
     2         ,nwa(1),nwb(1),iwa(1),iwb(1)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
c
      maxrp3 = maxrpl*3
      nvec = nvec1+nvec2
      if (ndens.ne.0) call vclr(gam,1,nint1+nint2)
c
      iuaa=0
      iubb=0
c
      do 220 isyma=1,nirrr
      if (nstra(isyma).eq.0) goto 220
      nalock = (nstra(isyma)-1)/nablk+1
      do 210 isymb=1,nirrr
      if (nstrb(isymb).eq.0) goto 210
      nblock = (nstrb(isymb)-1)/nbblk+1
c... symmetry of intermediate state
      isymk = mults(isyma,isymb)
c... symmetry of orbital excitations
      isymr = mults(isss,isymk)
      if (ndens+nham.lt.0 .and. isymr.ne.1) goto 210
c...  symmetry of excited beta strings
      isymbe = mults (isymb,isymr)
      nstrbe = nstrb(isymbe)
      mtaa=0
      do 200 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)
c...  generate alpha replacements of required symmetry
      mta=mtaa
      do 10 i=1,na
10    icga(i) = icgaa(i)
      iua=iuaa
      iofiwa = 1
      do 30 ia=1,naa
20    call string (mta,nact,na,icga,itypea,jmta,iua,mults)
      if (jmta.ne.isyma) goto 20
      if (ia.eq.1) ic1 = ic(mta)
      call onel (icga,isymr,iwa(iofiwa),nwa(ia),na,ic3e,ic,inter)
30    iofiwa = iofiwa + maxrp3
c
      mtbb=0
      do 180 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa*nbb
      naabbp = naabb*npair(isymr)
c
c... loop over beta strings in block storing single replacements
      mtb = mtbb
      iub = iubb
      do 40 i=1,nb
40    icgb(i)=icgbb(i)
      iofiwb = 1
      do 60 ib=1,nbb
50    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mults)
      if (jmtb.ne.isymb) goto 50
      if (ib.eq.1) itz1 = itz(mtb)-1
      call onel(icgb,isymr,iwb(iofiwb),nwb(ib),nb,ic3e,itz,inter(1,1,2))
60    iofiwb = iofiwb + maxrp3
c
c...  d matrix
      do 160 ivec=1,nvec1
      call vclr(d,1,naabbp)
      iofiw1 = 0
      iofca = (ivec-1)*nci+itz1
      iofd = - naabb
      do 850 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 840 iw=1,nwa(ia)
      ioffc = iofca + iwa(iofiw1+iw)
      ioffd = iofd + iwa(iofiw3+iw)*naabb
      if (iwa(iofiw2+iw).gt.0) then
      do 810 ib=1,nbb
810   d(ioffd+ib) = d(ioffd+ib) + c1(ioffc+ib)
      else
      do 830 ib=1,nbb
830   d(ioffd+ib) = d(ioffd+ib) - c1(ioffc+ib)
      end if
840   continue
      iofiw1 = iofiw1 + maxrp3
850   iofd = iofd + nbb
      iofd = 1-nbb-naabb
      iofcb = (ivec-1)*nci+ic1-nstrbe
      iofiw1 = 0
      do 950 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 940 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + iwb(iofiw3+iw)*naabb
      if (iwb(iofiw2+iw).gt.0) then
      do 910 ia=1,naa
910   d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c1(ioffc+ia*nstrbe)
      else
      do 930 ia=1,naa
930   d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c1(ioffc+ia*nstrbe)
      end if
940   continue
      iofiw1 = iofiw1 + maxrp3
950   iofd = iofd + 1
      iof1d = 1
      iof1c = (ivec-1)*nci+ic1+itz1+1
      if (ndens) 70,120,90
c...  1 particle matrix only to be done
70    do 80 ia=1,naa
      call mxmb (d(iof1d),naabb,1,
     1 c2(iof1c),1,0,
     2 gam(nint2+1),1,0, npair(isymr),nbb,1)
      iof1d = iof1d + nbb
80    iof1c = iof1c + nstrbe
      goto 150
c...  one and two particle density matrix
90    continue
         if(nvec.gt.1) then
      call dscal(nint1+nint2,1.0d0/weight(ivec),gam,1)
          endif
      call mxmb (d,naabb,1, d,1,naabb,
     1 gam(intoff(isymr)+1),1,npair(isymr),
     2 npair(isymr),naabb,npair(isymr) )
      if (isymr.ne.1) goto 110
      do 100 ia=1,naa
      call mxmb (d(iof1d),naabb,1,
     1 c1(iof1c),1,0,
     2 gam(nint2+1),1,0, npair(isymr),nbb,1)
      iof1d = iof1d + nbb
100   iof1c = iof1c + nstrbe
110   continue
      if (nvec.gt.1) then
      call dscal(nint1+nint2,weight(ivec),gam,1)
           endif
120   continue
c
      if (nham.le.0) goto 150
      call vclr(e,1,naabbp)
      call mxms (d,1,naabb, zint1(intoff(isymr)+1),1,npair(isymr),
     1 e,1,naabb, naabb,npair(isymr),npair(isymr) )
      if (isymr.ne.1) goto 140
      iof1d = 1
      iof1c = (ivec-1)*nci+ic1+itz1+1
      do 130 ia=1,naa
      call mxmb (c1(iof1c),1,0,
     1 zint1(nint2+1),0,1, e(iof1d),1,naabb,
     2 nbb,1,npair(isymr) )
      iof1d = iof1d + nbb
130   iof1c = iof1c + nstrbe
140   continue
      iofca = (ivec-1)*nci+itz1
      iofd = - naabb
      iofiw1 = 0
      do 1050 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 1040 iw=1,nwa(ia)
      ioffc = iofca + iwa(iofiw1+iw)
      ioffd = iofd + iwa(iofiw3+iw)*naabb
      if (iwa(iofiw2+iw).gt.0) then
      do 1010 ib=1,nbb
1010  s(ioffc+ib) = s(ioffc+ib) + e(ioffd+ib)
      else
      do 1030 ib=1,nbb
1030  s(ioffc+ib) = s(ioffc+ib) - e(ioffd+ib)
      end if
1040   continue
      iofiw1 = iofiw1 + maxrp3
1050   iofd = iofd + nbb
      iofd = 1-nbb-naabb
      iofcb = (ivec-1)*nci+ic1-nstrbe
      iofiw1 = 0
      do 1150 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 1140 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + iwb(iofiw3+iw)*naabb
      if (iwb(iofiw2+iw).gt.0) then
      do 1110 ia=1,naa
1110  s(ioffc+ia*nstrbe) = s(ioffc+ia*nstrbe) + e(ioffd+ia*nbb)
      else
      do 1130 ia=1,naa
1130  s(ioffc+ia*nstrbe) = s(ioffc+ia*nstrbe) - e(ioffd+ia*nbb)
      end if
1140   continue
      iofiw1 = iofiw1 + maxrp3
1150  iofd = iofd + 1
150   continue
160   continue
c... end of beta block
      do 170 i=1,nb
170   icgbb(i)=icgb(i)
      iubb=iub
180   mtbb = mtb
c...  end of alpha block
      do 190 i=1,na
190   icgaa(i)=icga(i)
      iuaa=iua
200   mtaa=mta
c... end of beta symmetry
210   continue
c...  end of alpha symmetry
220   continue
c      if (nham.gt.0) call outvec (s,nci,'s')
      return
      end
      subroutine detint (zin,zout,ic2e,modone)
      implicit real*8  (a-h,o-z)
      logical modone
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension zin(1),zout(1),ic2e(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      ij=0
      do 10 i=1,nact
      do 10 j=1,i
      isymij = mults(itypea(i),itypea(j))
      ij=ij+1
      ijkl = ic2e(ij)
      do 10 k=1,nact
      do 10 l=1,k
      if (mults(itypea(k),itypea(l)).ne.isymij) goto 10
      ijkl=ijkl+1
      zout(ijkl)  = 0.5d0 * zin(ind(ilifa(i)+j,ilifa(k)+l))
10    continue
      ij = 0
      do 30 i=1,nact
      do 30 j=1,i
      if (itypea(i).ne.itypea(j)) goto 30
      ij=ij+1
      g = zin(ic1d+ilifa(i)+j)
      if (modone) then
      do 20 k=1,nact
20    g = g - zin (ind(ilifa(i)+k,ilifa(j)+k))*0.5d0
      end if
      zout(nint2+ij) = g
30    continue
      return
      end
      subroutine spinad (c,iref,ic,itz,inter)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension c(nci),ic(1),itz(1),inter(na,nact,2)
      dimension icga(31),icgb(31),ialp(31),ibet(31)
      call detocc (iref,icga,icgb)
      nalp=0
      nbet=0
      do 70 i=1,nact
      icase=1
      do 30 j=1,na
      if (icga(j).ne.i) goto 30
      icase=3
      ja=j
30    continue
      do 40 j=1,nb
      if (icgb(j).ne.i) goto 40
      icase=icase+1
      jb=j
40    continue
      goto (70,50,60,70),icase
50    nbet=nbet+1
      ibet(nbet)=jb
      goto 70
60    nalp=nalp+1
      ialp(nalp)=ja
      goto 70
70    continue
      do 80 i=1,nci
80    c(i)=0.0d0
      fac = 1.0d0/ dsqrt( dble(2**nbet))
      if (nbet.gt.5) stop 5
      do 130 i5=1,2
      do 120 i4=1,2
      do 110 i3=1,2
      do 100 i2=1,2
      do 90 i1=1,2
      ipar = 1
c      call outive (icga,na,'alpha string')
      mta = istrad (ipar,inter(1,1,1),icga,na)
c      call outive (icgb,nb,'beta string')
      mtb = istrad(ipar,inter(1,1,2),icgb,nb)
      c(ic(mta)+itz(mtb)) =  dble(ipar)*fac
      if (nbet.lt.1) goto 140
      j1=icgb(ibet(1))
      icgb(ibet(1))=icga(ialp(1))
      icga(ialp(1))=j1
90    continue
      if (nbet.lt.2) goto 140
      j1=icgb(ibet(2))
      icgb(ibet(2))=icga(ialp(2))
      icga(ialp(2))=j1
100   continue
      if(nbet.lt.3) goto 140
      j1=icgb(ibet(3))
      icgb(ibet(3))=icga(ialp(3))
      icga(ialp(3))=j1
110   continue
      if (nbet.lt.4) goto 140
      j1=icgb(ibet(4))
      icgb(ibet(4))=icga(ialp(4))
      icga(ialp(4))=j1
120   continue
      if (nbet.lt.5) goto 140
      j1=icgb(ibet(5))
      icgb(ibet(5))=icga(ialp(5))
      icga(ialp(5))=j1
130   continue
140   continue
      return
      end
      function istrad (ipar,inter,icg,n)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension icg(n),inter(na,nact)
      dimension icgg(31)
      do 10 j=1,n
10    icgg(j) = icg(j)
      if (n.lt.2) goto 40
      do 30 jj=2,n
      jjm = jj-1
      i = jjm
      do 20 ii=1,jjm
      if (icgg(i).lt.icgg(i+1)) goto 20
      j1=icgg(i)
      icgg(i)=icgg(i+1)
      icgg(i+1)=j1
      ipar=-ipar
20    i=i-1
30    continue
40    continue
      istrad=1
      do 50 j=1,n
50    istrad = istrad+inter(j,icgg(j))
      return
      end
      subroutine diags (s,zint,ic2e,ic3e,zb,z,f)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension s(1),zint(1),ic2e(1),ic3e(1)
     1         ,zb(maxbb,nact),z(maxbb,nact),f(maxbb)
      dimension h(31),icga(31),icgb(31),za(31)
c
c...  load 1 elec integrals
      do 10 i=1,nact
10    h(i) = zint(nint2+ic3e(((i+1)*i)/2))
c
      iccc=0
      call vclr(s,1,nci)
      do 170 isymb=1,nirrr
      isyma=mults(isymb,isss)
      nsa = nstra(isyma)
      nsb = nstrb(isymb)
      if (nsa*nsb.le.0) goto 170
      call vclr(zb,1,nact*maxbb)
c..   store beta string occupancies
      mtb=0
      do 40 ib=1,nsb
20    call string (mtb,nact,nb,icgb,itypea,jmtb,iua,mults)
      if (jmtb.ne.isymb) goto 20
      do 30 ii=1,nb
30    zb(ib,icgb(ii))=1.0d0
40    continue
c
c     loop over alpha strings
      mta=0
      do 160 ia=1,nsa
50    call string (mta,nact,na,icga,itypea,jmta,iua,mults)
      if (jmta.ne.isyma) goto 50
      call vclr(za,1,nact)
      do 60 ii=1,na
60    za(icga(ii)) = 1.0d0
c
c...  occupation numbers of determinants
      do 70 i=1,nact
      do 70 ib=1,nsb
70    z(ib,i) = zb(ib,i) + za(i)
c
c...  constant exchange contribution & coulomb
      call mxmb (z,1,maxbb, h,1,0, s(iccc+1),1,0, nsb,nact,1)
      do 80 i=1,nact
      do 80 j=1,i
      ij=icf(i)+j
      zz = -zint(ic2e(ij) +ic3e(ij))+
     1    2.0d0*zint(ic2e((i*(i+1))/2)+ic3e((j*(j+1))/2))
      if (i.eq.j) zz = zz*0.5d0
      do 80 ib=1,nsb
80    s(iccc+ib) = s(iccc+ib) + zz*z(ib,i)*z(ib,j)
c
      vv = (na-nb)**2
c..   z to hold 1 if singly occ, 0 otherwise
      call vclr(f,1,nsb)
      do 100 i=1,nact
      do 90 ib=1,nsb
90    z(ib,i) = z(ib,i) * (2.0d0-z(ib,i))
      do 100 ib=1,nsb
100   f(ib) = f(ib) + z(ib,i)
c...  trap division by zero for 0 or 1 singly occ orbital
      do 110 ib=1,nsb
110   f(ib) = dmax1(f(ib),1.1d0)
      do 120 ib=1,nsb
120   f(ib) = (vv-f(ib)) / ( f(ib) * (f(ib)-1.0d0) )
c
      do 150 i=1,nact
      if (i.eq.1) goto 140
      im=i-1
      do 130 j=1,im
      ij = icf(i)+j
      zz = -zint(ic2e(ij)+ic3e(ij))
      do 130 ib=1,nsb
130   s(iccc+ib) = s(iccc+ib) + zz*f(ib)*z(ib,i)*z(ib,j)
140   ii = (i*(i+1))/2
      zz = -0.5d0*zint(ic2e(ii)+ic3e(ii))
      do 150 ib=1,nsb
150   s(iccc+ib) = s(iccc+ib) + zz*z(ib,i)
160   iccc = iccc+nsb
170   continue
      return
      end
      subroutine detpp (zint,ic2e,ic3e,ic,h,coul2,coulme)
      implicit real*8  (a-h,o-z)
      integer patern,patset,pat1,pat2,dif,bra,zket
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension zint(ne),ic2e(1),ic3e(1),ic(1),h(nact,nact)
     1         ,coul2(nact,nact),coulme(nact,nact)
      dimension icga(33),icgb(33),patern(maxp)
      logical odd,parity
      integer popcnt
      odd(i) = (i/2)*2.ne.i
      z(i,j,k,l)=zint(ic2e(icf(max(i,j))+min(i,j))
     1               +ic3e(icf(max(k,l))+min(k,l)) )
      ind(i,j)=icf(max(i,j))+min(i,j)
c
c...  set up bare h0, coulomb, exchange
c     t0=seccpu()
      do 10 i=1,nact
      do 10 j=1,nact
      coul = z(i,i,j,j)
      exch = z(i,j,i,j)
      coul2(i,j) = coul*2.0d0
      coulme(i,j) = coul-exch
      if (itypea(i).eq.itypea(j)) h(i,j) = zint(nint2+ic3e(ind(i,j)))
10    continue
c
c...  loop over bras
      do 240 mtt1=1,nplist
      mt1 = iplist(mtt1)
c.... set up bit pattern for this determinant
      k=0
      do 20 i=1,nstraa
      if (ic(i).ge.mt1 .or. k.gt.ic(i)) goto 20
      k=ic(i)
      istra=i
20    continue
      mta=0
      do 30 mt=1,istra
30    call string (mta,nact,na,icga,itypea,jmta,iu,mults)
      if(odebug(34)) then
       call outive (icga,na,'alpha string')
      endif
      isymb = mults(jmta,isss)
      if (nstrb(isymb).eq.0) goto 240
      mtb=0
      do 50 mt=1,mt1-k
40    call string (mtb,nact,nb,icgb,itypea,jmtb,iu,mults)
      if (jmtb.ne.isymb) goto 40
50    continue
      if(odebug(34)) then
       call outive (icgb,nb,'beta string')
      endif
      pat1 = patset (icga,icgb)
      patern(mtt1) = pat1
c
c...  loop over kets
      do 230 mtt2=1,mtt1
      pat2 = patern(mtt2)
      dif=ieor(pat1,pat2)
      ndif=popcnt(dif)
      if (ndif.gt.4) goto 230
      if (ndif-2) 170,100,60
c
c...  double excitation
c...  i11 -> i21     i12 -> i22
60    bra = iand(dif,pat1)
      zket = iand(dif,pat2)
      i11 = leadz(bra)+1
      i12 = leadz(ishft(bra,i11))+1+i11
      i21 = leadz(zket)+1
      i22 = leadz(ishft(zket,i21))+1+i21
      if (i11.gt.nact) goto 80
      if (i12.gt.nact) goto 70
c...  a->a a->a
      val = z(i11,i21,i12,i22)-z(i11,i22,i12,i21)
      goto 90
c...  a->a b->b
70    val = z(i11,i21,i12-nact,i22-nact)
      goto 90
c...  b->b b->b
80    val = z(i11-nact,i21-nact,i12-nact,i22-nact)
     1    - z(i11-nact,i22-nact,i12-nact,i21-nact)
90    i1=leadz(dif)+1
      i2=leadz(ishft(dif,i1))+1+i1
      i3=leadz(ishft(dif,i2))+1+i2
      i4=leadz(ishft(dif,i3))+1+i3
      if( odd( popcnt(ishft(ishft(pat1,i2-65),65-i2+i1))
     1       + popcnt(ishft(ishft(pat1,i4-65),65-i4+i3)) ) )val=-val
      val=val*2
      goto 220
c
c...  single excitation
100   i11=leadz(dif)+1
      i21=leadz(ishft(dif,i11))+1+i11
      if(odebug(34)) then
       print *,'single ',i11,i21
      endif
      parity = odd(popcnt(ishft(ishft(pat1,i21-65),65-i21+i11)))
      val = 0.0d0
      if (i11.gt.nact) goto 130
c...  a->a
      ii=ic3e(icf(i21)+i11)
      do 110 j1=1,na
      j=icga(j1)
110    val = val + zint(ic2e(icf(j)+j)+ii) - z(j,i11,j,i21)
      do 120 j1=1,nb
      j=icgb(j1)
120   val = val + zint(ic2e(icf(j)+j)+ii)
      goto 160
c...  b->b
130   i21=i21-nact
      i11=i11-nact
      ii=ic3e(icf(i21)+i11)
      do 140 j1=1,na
      j=icga(j1)
140   val = val + zint(ic2e(icf(j)+j)+ii)
      do 150 j1=1,nb
      j=icgb(j1)
150   val = val + zint(ic2e(icf(j)+j)+ii) - z(j,i11,j,i21)
160   val = (val*2 + h(i11,i21) )
      if (parity) val=-val
      goto 220
c
c...  no spin orbital difference
170   val=0.0d0
      do 190 i1=1,na
      i=icga(i1)
      do 180 j1=1,na
      j=icga(j1)
180   val = val + coulme(j,i)
      val = val + h(i,i)
      do 190 j1=1,nb
      j=icgb(j1)
190   val = val + coul2(j,i)
      do 210 i1=1,nb
      i=icgb(i1)
      do 200 j1=1,nb
      j=icgb(j1)
200   val = val + coulme(j,i)
210   val = val + h(i,i)
c
220   hpp(mtt1,mtt2) = val
      hpp(mtt2,mtt1) = val
230   continue
240   continue
      if(odebug(34)) then
       call outsqr (hpp,maxp,nplist,nplist,'hpp')
      endif
      return
      end
      subroutine detpq (hrow,sigma,zint,ic2e,ic3e,ic,itz
     1                 ,nvec,iwa,iwb,inter)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
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
      dimension zint(1),ic2e(1),ic3e(1),ic(1),itz(1)
     1         ,hrow(nci),sigma(nci,nvec),inter(na,nact,2)
     2         ,iwa(maxrpl,3),iwb(maxrpl,3)
      dimension icga(31),icgb(31)
c
      call dscal(nint2,2.0d0,zint,1)
c
c     t0=seccpu()
c
c...  loop over p states
      do 110 mtt1=1,nplist
      mt1 = iplist(mtt1)
c.... set up occ pattern for this determinant
      k=0
      do 10 i=1,nstraa
      if (ic(i).ge.mt1 .or. k.gt.ic(i)) goto 10
      k=ic(i)
      istra=i
10    continue
      mta=0
      do 20 mt=1,istra
20    call string (mta,nact,na,icga,itypea,jmta,iu,mults)
c      call outive (icga,na,'alpha string')
      isymb = mults(jmta,isss)
      if (nstrb(isymb).eq.0) goto 110
      mtb=0
      do 40 mt=1,mt1-k
30    call string (mtb,nact,nb,icgb,itypea,jmtb,iu,mults)
      if (jmtb.ne.isymb) goto 30
40    continue
c      call outive (icgb,nb,'beta string')
      icmta  =  ic(mta)
      itzmtb = itz(mtb)
      call vclr(hrow,1,nci)
c
      do 90 isymr=1,nirrr
c      print*,'detpq do 90 isymr=',isymr
      call onel(icga,isymr,iwa,nwa,na,ic3e,ic,inter)
      call onel(icgb,isymr,iwb,nwb,nb,ic3e,itz,inter(1,1,2))
c      call outisq (iwa,maxrpl,nwa,3,'iwa')
c      call outisq (iwb,maxrpl,nwb,3,'iwb')
      if (isymr.eq.1) then
c.... one electron
      do 50 iw=1,nwa
      iad = iwa(iw,1)+itzmtb
50    hrow(iad) = hrow(iad) +  dble(iwa(iw,2))*zint(nint2+iwa(iw,3))
      do 60 iw=1,nwb
      iad = iwb(iw,1)+icmta
60    hrow(iad) = hrow(iad) +  dble(iwb(iw,2))*zint(nint2+iwb(iw,3))
c      print*,'one electron done'
      end if
c...  two electron a->a b->b
      do 80 iw1=1,nwa
      par =  dble(iwa(iw1,2))
      io = intoff(isymr)+(iwa(iw1,3)-1)*npair(isymr)
c      print*,'iw1,io,par ',iw1,io,par
      do 70 iw2=1,nwb
      iad = iwa(iw1,1)+iwb(iw2,1)
c      print*,'iw2,iad ',iw2,iad
70    hrow(iad)=hrow(iad)+par*dble(iwb(iw2,2))*zint(io+iwb(iw2,3))
80    continue
90    continue
c
c.... two electron a->a a->a  &  b->b b->b
      call twoel(icga,hrow,itzmtb,na,ic2e,ic3e,ic,inter,zint)
      call twoel(icgb,hrow,icmta ,nb,ic2e,ic3e,itz,inter(1,1,2),zint)
c      call outvec (hrow,nci,'hamiltonian row')
      do 100 ivec=1,nvec
      call daxpy(nci,cp(mtt1,ivec)
     +  ,hrow,1,sigma(1,ivec),1)
100   continue
110   continue
c      print*,'cipq time: ',seccpu()-t0
      return
      end
      subroutine detocc (icc,icga,icgb)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension icga(na),icgb(nb)
      ii=0
      do 10 ibs=1,nirrr
      ias = mults(ibs,isss)
      jj = ii + nstra(ias)*nstrb(ibs)
      if (ii.lt.icc.and.jj.ge.icc)  goto 20
10    ii = jj
20    ibs = mults(ias,isss)
      ia = (icc-ii-1)/nstrb(ibs)+1
      ib = icc-ii-(ia-1)*nstrb(ibs)
      mta=0
      do 40 ii=1,ia
30    call string (mta,nact,na,icga,itypea,jmta,iua,mults)
      if (jmta.ne.ias) goto 30
40    continue
      mtb=0
      do 60 ii=1,ib
50    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mults)
      if (jmtb.ne.ibs) goto 50
60    continue
      return
      end
      subroutine twoel (icg,hrow,ic1,n,ic2e,ic3e,itz,inter,zint)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension icg(n),hrow(nci),ic2e(1),ic3e(1),itz(1),inter(na,nact)
     1         ,zint(ne)
      dimension ia(33),ib(33)
      ind(i,j)=icf(max(i,j))+min(i,j)
      do 10 i=1,n
10    ib(i)=icg(i)
      ib(n+1)=0
      ib(n+2)=0
      ib1=ib(1)
      ib2=ib(2)
      do 70 i1=1,n-1
      do 20 i=1,n+2
20    ia(i)=ib(i)
      ia1=ib1
      ia2=ib2
      parit=1.0d0
      do 60 i2=i1+1,n
      kk=1
      do 30 j=1,n
30    kk=kk+inter(j,ia(j))
      kk1=1
      ita=mults(itypea(ia1),itypea(ia2))
      do 50 j1=1,nact-1
      jj=kk-inter(kk1,ia(kk1))+inter(kk1,j1)
      if (j1.eq.ia(kk1+2)) then
       kk1=kk1+1
       kk=jj
      else
       kk2=kk1+1
       parity=parit
       ma1=ic2e(ind(j1,ia1))
       mb1=ic2e(ind(j1,ia2))
       ita1=mults(ita,itypea(j1))
       do 40 j2=j1+1,nact
       ii=jj-inter(kk2,ia(kk2))+inter(kk2,j2)
       if (j2.eq.ia(kk2+1)) then
        kk2=kk2+1
        parity=-parity
        jj=ii
       else
c        if(ita1.eq.itypea(j2))print*,'twoel;i1,i2,j1,j2,ia1,ia2,ic1,ii,
c     <  itz(ii),parity,zint,zint,contrib '
c     >,i1,i2,j1,j2,ia1,ia2,ic1,ii,itz(ii)
c     >,parity,zint(ma1+ic3e(ind(j2,ia1)))
c     >,zint(mb1+ic3e(ind(j2,ia2)))
c     >,parity*(zint(ma1+ic3e(ind(j2,ia2)))-zint(mb1+ic3e(ind(j2,ia1))))
        if (ita1.eq.itypea(j2)) hrow(ic1+itz(ii)) = hrow(ic1+itz(ii))
     1  + parity * (zint(ma1+ic3e(ind(j2,ia2)))
     2             -zint(mb1+ic3e(ind(j2,ia1))) )
       end if
40     continue
      end if
50    continue
      parit=-parit
      ia2=ia(i2+1)
      ia(i2+1)=ia(2)
      ia(2)=ia2
60    continue
      ib1=ib2
      ib2=ib(i1+2)
      ib(i1+2)=ib(1)
      ib(2)=ib2
      ib(1)=ib1
70    continue
      return
      end
      subroutine onels (icg,isymr,iw,nw,n,itz,inter)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
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
      dimension icg(n),iw(maxrpl,3),itz(1),inter(na,nact)
      dimension ia(32)
      do 10 i=1,n
10    ia(i)=icg(i)
      ia1=ia(1)
      ia(n+1)=0
      iz=1
      nw=0
      do 40 i=1,n
      isymj=mults(isymr,itypea(ia1))
      kk=1
      iz2=iz
      jj=1
      do 20 j=1,n
20    jj=jj+inter(j,ia(j))
      do 30 j1=1,nact
      ii=jj-inter(kk,ia(kk))+inter(kk,j1)
      if (j1.eq.ia(kk+1)) then
       kk=kk+1
       jj=ii
       iz2=-iz2
      else if (itypea(j1).eq.isymj) then
       nw=nw+1
       iw(nw,1) = itz(ii)
       iw(nw,2) = iz2
       iw(nw,3) = (ia1-1)*nact+j1
      end if
30    continue
      iz = -iz
      ia1 = ia(i+1)
      ia(i+1)=ia(1)
40    ia(1)=ia1
      return
      end
      subroutine ssq (c1,s,nvec,ic,itz,inter,iwa,iwb,d,nablk,nbblk,
     3 nwa,nwb)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension c1(1),s(1)
     1         ,ic(1),itz(1),inter(na,nact,2),d(1)
     2         ,nwa(1),nwb(1),iwa(1),iwb(1)
      dimension icga(31),icgb(31),icgaa(31),icgbb(31)
c
      maxrp3 = maxrpl*3
      sz = ( dble(na-nb))*0.5d0
      fac =  dble(nb) + sz*(sz+1.0d0)
      call vsmul(c1,1,fac,s,1,nci*nvec)
c
      iuaa=0
      iubb=0
c
      do 220 isyma=1,nirrr
      if (nstra(isyma).eq.0) goto 220
      nalock = (nstra(isyma)-1)/nablk+1
      do 210 isymb=1,nirrr
      if (nstrb(isymb).eq.0) goto 210
      nblock = (nstrb(isymb)-1)/nbblk+1
c... symmetry of intermediate state
      isymk = mults(isyma,isymb)
c... symmetry of orbital excitations
      isymr = mults(isss,isymk)
c...  symmetry of excited beta strings
      isymbe = mults(isymb,isymr)
      nstrbe = nstrb(isymbe)
      mtaa=0
      do 200 ialock=1,nalock
      naa = min(nablk,nstra(isyma)-(ialock-1)*nablk)
c...  generate alpha replacements of required symmetry
      mta=mtaa
      do 10 i=1,na
10    icga(i) = icgaa(i)
      iua=iuaa
      iofiwa = 1
      do 30 ia=1,naa
20    call string (mta,nact,na,icga,itypea,jmta,iua,mults)
      if (jmta.ne.isyma) goto 20
      if (ia.eq.1) ic1 = ic(mta)
      call onels (icga,isymr,iwa(iofiwa),nwa(ia),na,ic,inter)
30    iofiwa = iofiwa + maxrp3
c
      mtbb=0
      do 180 iblock=1,nblock
      nbb = min(nbblk,nstrb(isymb)-(iblock-1)*nbblk)
      naabb = naa*nbb
      naabbp = naabb*nact**2
c
c... loop over beta strings in block storing single replacements
      mtb = mtbb
      iub = iubb
      do 40 i=1,nb
40    icgb(i)=icgbb(i)
      iofiwb = 1
      do 60 ib=1,nbb
50    call string (mtb,nact,nb,icgb,itypea,jmtb,iub,mults)
      if (jmtb.ne.isymb) goto 50
      if (ib.eq.1) itz1 = itz(mtb)-1
      call onels(icgb,isymr,iwb(iofiwb),nwb(ib),nb,itz,inter(1,1,2))
60    iofiwb = iofiwb + maxrp3
c
c...  d matrix
      do 160 ivec=1,nvec
      call vclr(d,1,naabbp)
      iofd = 1-nbb-naabb
      iofcb = (ivec-1)*nci+ic1-nstrbe
      iofiw1 = 0
      do 950 ib=1,nbb
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 940 iw=1,nwb(ib)
      ioffc = iofcb + iwb(iofiw1+iw)
      ioffd = iofd + iwb(iofiw3+iw)*naabb
      if (iwb(iofiw2+iw).gt.0) then
      do 910 ia=1,naa
910   d(ioffd+ia*nbb) = d(ioffd+ia*nbb) + c1(ioffc+ia*nstrbe)
      else
      do 930 ia=1,naa
930   d(ioffd+ia*nbb) = d(ioffd+ia*nbb) - c1(ioffc+ia*nstrbe)
      end if
940   continue
      iofiw1 = iofiw1 + maxrp3
950   iofd = iofd + 1
c     iof1d = 1
c     iof1c = (ivec-1)*nci+ic1+itz1+1
      iofca = (ivec-1)*nci+itz1
      iofd = - naabb
      iofiw1 = 0
      do 1050 ia=1,naa
      iofiw2 = iofiw1 + maxrpl
      iofiw3 = iofiw2 + maxrpl
      do 1040 iw=1,nwa(ia)
      ioffc = iofca + iwa(iofiw1+iw)
      ioffd = iofd + iwa(iofiw3+iw)*naabb
      if (iwa(iofiw2+iw).gt.0) then
      do 1010 ib=1,nbb
1010  s(ioffc+ib) = s(ioffc+ib) - d(ioffd+ib)
      else
      do 1030 ib=1,nbb
1030  s(ioffc+ib) = s(ioffc+ib) + d(ioffd+ib)
      end if
1040   continue
      iofiw1 = iofiw1 + maxrp3
1050   iofd = iofd + nbb
160   continue
c... end of beta block
      do 170 i=1,nb
170   icgbb(i)=icgb(i)
      iubb=iub
180   mtbb = mtb
c...  end of alpha block
      do 190 i=1,na
190   icgaa(i)=icga(i)
      iuaa=iua
200   mtaa=mta
c... end of beta symmetry
210   continue
c...  end of alpha symmetry
220   continue
      return
      end
      subroutine fci0 (q,civec,zint,grad,diag)
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(*),zint(ne),diag(nci),grad(nci),civec(nci)
      dimension dum(1)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call vclr(diag,1,nci)
      call vclr(grad,1,nci)
c
      if (iexc.lt.0) goto 51
      call rdedx(surd,511,iblkft,numft)
      call find(numft)
      do 10 iblk=1,nblkd1
      call get (g,nw)
      call find (numft)
      call unpkft (nw)
      do 10 iw=1,nw
10    diag(ii(iw)) = diag(ii(iw)) + surd(ival(iw)) *
     1                  zint(ic1d+ilifa(iorb(iw))+iorb(iw))
c
      do 20 iblk=1,nblko1
      call get (g,nw)
      call find (numft)
      call unpkft (nw)
      do 20 iw=1,nw
      val = surd(ival(iw))*zint(ic1d+ilifa(iorb(iw))+jorb(iw))
      grad(ii(iw)) = grad(ii(iw)) + val*civec(jj(iw))
20    grad(jj(iw)) = grad(jj(iw)) + val*civec(ii(iw))
      call get (g,nw)
c
      call find (numft)
      do 30 iblk=1,nblkd2
      call get(g,nw)
      call find(numft)
      call unpkft(nw)
      do 30 iw=1,nw
30    diag(ii(iw)) = diag(ii(iw)) + surd(ival(iw))
     1  * zint(ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw)))
c
      do 40 iblk=1,nblko2
      call get (g,nw)
      call find(numft)
      call unpkft(nw)
      do 40 iw=1,nw
      val = surd(ival(iw)) * zint(ind(ilifa(iorb(iw))+jorb(iw),
     1                                ilifa(korb(iw))+lorb(iw) ))
      grad(ii(iw)) = grad(ii(iw)) + val*civec(jj(iw))
40    grad(jj(iw)) = grad(jj(iw)) + val*civec(ii(iw))
      call get (g,nw)
c
      zz = 0.0d0
      do 41 i=1,nci
41    zz = zz + (diag(i)*civec(i)+grad(i))*civec(i)
      zz = -zz
      hessen = zz
      do 50 i=1,nci
      diag(i) = diag(i) + zz
50    grad(i) = grad(i) + diag(i)*civec(i)
      if(lto(1)) call cid(q,zint,diag)
      goto 70
c
51    call detci (q(1),q(1),zint,dum, diag,dum,dum, 3,1)
      call detci (q(1),q(1),zint,dum, civec,dum,grad, 4,1)
      zz = -ddot(nci,civec,1,grad,1)
      call daxpy(nci,zz,civec,1,grad,1)
70    continue
      return
      end
      subroutine fci1 (q,vec,sigma,zint,zint1,gam,gam1,civec)
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
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
c
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension q(1),zint(ne),vec(nci),sigma(nci),gam(ne)
     1   ,gam1(ne),zint1(ne),civec(nci)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
c#### factor of -1 ######
c      call sscal (ne,-1.0,zint1,1)
c      call scaler (ne,-1.0,zint1,zint1)
c########################
      call vclr(gam1,1,ne)
      if (iexc.lt.0) goto 999
      call vsmul(vec,1,hessen,sigma,1,nci)
c...  read thro' ci formulae
      call rdedx(surd,511,iblkft,numft)
      call find(numft)
      call accnt('fci1d1',2)
      do 20 iblk=1,nblkd1
      call get (g,nw)
      call find (numft)
      call unpkft (nw)
      do 20 iw=1,nw
      intad = ic1d+ilifa(iorb(iw))+iorb(iw)
      sigma(ii(iw)) = sigma(ii(iw)) + (zint1(intad)*civec(ii(iw))
     1              +  zint(intad)*vec(ii(iw)))*surd(ival(iw))
      gam1(intad) = gam1(intad)+civec(ii(iw))*vec(ii(iw))*surd(ival(iw))
20    continue
c
      call accnt('fci1o1',2)
      do 50 iblk=1,nblko1
      call get (g,nw)
      call find (numft)
      call unpkft (nw)
      do 50 iw=1,nw
      intad =ic1d+ilifa(iorb(iw))+jorb(iw)
      sigma(ii(iw)) = sigma(ii(iw)) + surd(ival(iw)) *
     1 (zint1(intad)*civec(jj(iw))+zint(intad)*vec(jj(iw)))
      gam1(intad) = gam1(intad)+civec(ii(iw))*vec(jj(iw))*surd(ival(iw))
      intad = ic1d + ilifa(jorb(iw))+iorb(iw)
      sigma(jj(iw)) = sigma(jj(iw)) + surd(ival(iw)) *
     1(zint1(intad)*civec(ii(iw)) + zint(intad)*vec(ii(iw)) )
      gam1(intad) = gam1(intad)+civec(jj(iw))*vec(ii(iw))*surd(ival(iw))
50    continue
      call get (g,nw)
c
      call find (numft)
      call accnt('fci1d2',2)
      do 70 iblk=1,nblkd2
      call get(g,nw)
      call find(numft)
      call unpkft(nw)
      do 70 iw=1,nw
      if (civec(ii(iw))) 60,70,60
60    intad = ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw))
      sigma(ii(iw)) = sigma(ii(iw)) + (zint1(intad)*civec(ii(iw))
     1              + zint(intad)*vec(ii(iw)))*surd(ival(iw))
      gam1(intad)=gam1(intad)+civec(ii(iw))*vec(ii(iw))*surd(ival(iw))
70    continue
c
      call accnt('fci1o2',2)
      do 80 iblk=1,nblko2
      call get (g,nw)
      call find(numft)
      call fci1o2 (zint,zint1,gam1,vec,civec,sigma,nw)
80    continue
      call get (g,nw)
      call accnt('other   ',2)
c
      ijij=0
      do 90 ij=1,nact**2
      ijij=ijij+ij
90    gam1(ijij)=gam1(ijij)*2.0d0
      goto 1001
c
999   continue
      call vclr(sigma,1,nci)
      call detci (q(1),q(1),zint,gam,vec,civec,sigma, 5,1)
c
1001  continue
      zz = -ddot(nci,civec,1,vec,1)
      call daxpy(ne,zz,gam,1,gam1,1)
      e1=-hamilt(gam,zint1)
      call daxpy(nci,e1,civec,1,sigma,1)
      return
      end
      subroutine pricon(iq,vp,iplist,nplist,iwrite)
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
      common /lsort / g(511),ii(511),jj(511),ival(511),iorb(511)
     1              ,jorb(511),korb(511),lorb(511)
      dimension iq(1),iplist(1),vp(1)
c
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
c
      ibase=icori(nact*nplist)
      iocc=ibase-1
      ie=nact*nplist
      call setsto(ie,0,iq(iocc+1))
      call rdedx(surd,511,iblkft,numft)
      do 280 i=1,nblkd1
      call find(numft)
      call get (g,nw)
      call unpkft(nw)
      do 280 iw=1,nw
      io=iorb(iw)
      ic=ii(iw)
      do 275 ip=1,nplist
275   if(iplist(ip).eq.ic) goto 276
      goto 280
276   iq(iocc+io+(ip-1)*nact)=surd(ival(iw))+0.001d0
280   continue
      write(iwrite,277)
277   format(/' p-space configurations'//
     1  '  nr   conf   energy     occupancy'/)
      do 282 ip=1,nplist
      write(iwrite,281)ip,iplist(ip),vp(ip),(iq(iocc+io),io=1,nact)
281   format(1x,i3,i6,5x,f10.5,30i2)
282   iocc=iocc+nact
      write(iwrite,283)
283   format()
      call corlsi(ibase)
      return
      end
      subroutine creatp(c)
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
      parameter (maxp=40,maxst=5,maxpst=maxp+maxst)
      common /three  / cp(maxp,maxst),hpp(maxp,maxp)
     >               ,iplist(maxp),nplist,jcmax(maxst),dp(maxp)
     >               ,iblkp,iblkpq,vp(maxp),npread,ipread(maxp)
     >               ,icon,icend(maxp)
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
      dimension c(nci,nstate)
      do 20 i=1,nstate
      do 20 jj=1,nplist
      j=iplist(jj)
20    c(j,i)=cp(jj,i)
      return
      end
      subroutine detgam (zin,zout,ic2e,ic3e)
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
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
      dimension zin(1),zout(1),ic2e(1),ic3e(1)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call vclr(zout,1,ne)
c      call outvec (zin,nint1+nint2,'density given to detgam')
      ii=0
      do 40 i=1,nact
      ii=ii+i
      zin(ic3e(ii)+nint2) = zin(ic3e(ii)+nint2)*2.0d0
      nnn = npair(1)
      call dscal(nnn,2.0d0,zin(ic3e(ii)),nnn)
      ii1 = ic2e(ii)+1
      call dscal(nnn,2.0d0,zin(ii1),1)
40    continue
      call dscal(nint2,0.25d0,zin,1)
      call dscal(nint1,0.5d0,zin(nint2+1),1)
c      call outvec (zin,nint1+nint2,'scaled zin')
c
      do 60 i=1,nact
      do 60 j=1,nact
      ij=ind(i,j)
      ijj = ilifa(i)+j
      isymij = mults(itypea(i),itypea(j))
      if (isymij.eq.1) zout(ic1d+ijj) = zin(nint2+ic3e(ij))
      do 60 k=1,i
      lmax=nact
      if (k.eq.i) lmax=j
      isymjk = mults(isymij,itypea(k))
      do 50 l=1,lmax
      if (itypea(l).ne.isymjk) goto 50
      kl=ind(k,l)
      kll=ilifa(k)+l
      zout(ind(ijj,kll)) = zin(ic2e(ij)+ic3e(kl))
      if (i.eq.k) zout(ind(ijj,kll)) = zout(ind(ijj,kll))
     1                -0.25d0*zin(nint2+ic3e(ind(j,l)))
      if (j.eq.k) zout(ind(ijj,kll)) = zout(ind(ijj,kll))
     1                -0.25d0*zin(nint2+ic3e(ind(i,l)))
      if (i.eq.l) zout(ind(ijj,kll)) = zout(ind(ijj,kll))
     1                -0.25d0*zin(nint2+ic3e(ind(j,k)))
      if (j.eq.l) zout(ind(ijj,kll)) = zout(ind(ijj,kll))
     1                -0.25d0*zin(nint2+ic3e(ind(i,k)))
50    continue
60    continue
      return
      end
      subroutine getbas
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
c
      integer mach, mxtask, limit1, limit2, limit3
      common/machin/mach(20),mxtask,limit1,limit2,limit3
c
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
      common /mccore/ intrel,lword,ltop,lmax,lmin,mreal
c
      integer nirr, mult, isymao, isymmo, irrr, iss
      integer mstart, nfin, nmc
      logical symm_diag
      common /gjs/ nirr,mult(8,8),isymao(maxorb),
     + isymmo(maxorb),irrr,iss,mstart(8),nfin(8),nmc,
     + symm_diag
c
c
      common /lsort / ablkin(511)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      call secini(iblk3,num3)
      call sectst (isect(490),itest)
      if (itest.eq.0) then
        call caserr('mcscf requires symmetry adaption')
      else
c..   symmetry adaption has been done - read it in
        call secget(isect(490),51,ibl190)
        nav = lenwrd()
        call readi (nirr,mach(13)*nav,ibl190,num3)
        nirrr=nirr
        call icopy (maxorb,isymao,1,symaos,1)
      end if
c...  read in nuclear repulsion energy
      call secget (isect(492),2,ibl192)
      call rdedx(ablkin,511,ibl192,num3)
      potnuc=ablkin(1)
      return
      end
      subroutine get1 (q,x,ni,zn)
c
c...  load one electron integrals, symmetry packed (symmetric parts only
c     ni=1:overlap, 2:ke, 3:h(1), 4:x, 5:y, 6:z
c
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
      dimension q(*),x(*)
      common /lsort / znuc(4),zz(507)
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
      data m511/511/
      ibase = icorr(0)
      intrel = icori(0)
      junk = icorr(1)
c
c...  to hold raw integrals
      ix = icorr(lenbas)
      call secget (ions2,2,iblk)
      call rdedx(znuc,m511,iblk,num3)
c
      lenb = lensec(lenbas)
      if (ni.eq.1) then
       iblock = iblk + 1 
       call rdedx(q(ix),lenbas,iblock,num3)
      elseif (ni.eq.2) then
       iblock = iblk + 1 + lenb
       call rdedx(q(ix),lenbas,iblock,num3)
      elseif (ni.eq.3) then
       iblock = iblk + 1 + lenb*2
       call rdedx(q(ix),lenbas,iblock,num3)
      elseif (ni.eq.4) then
       iblock = iblk + 1 + lenb*3
       call rdedx(q(ix),lenbas,iblock,num3)
       call dscal(lenbas,-1.0d0,q(ix),1)
      elseif (ni.eq.5) then
       iblock = iblk + 1 + lenb*4
       call rdedx(q(ix),lenbas,iblock,num3)
       call dscal(lenbas,-1.0d0,q(ix),1)
      elseif (ni.eq.6) then
       iblock = iblk + 1 + lenb*5
       call rdedx(q(ix),lenbas,iblock,num3)
       call dscal(lenbas,-1.0d0,q(ix),1)
      else
       call caserr('indexing error in get1')
      endif
c
c.... copy nuclear values into zn
c
      zn = 0.0d0
      if (ni.ge.3) zn = znuc(ni-2)
c
c...  now must symmetry pack
      jj = 0
      ii = 1
      ixx = ix
      do 100 isym=1,nirrr
      n = nsymao(isym)
      do 99 i=1,n
      ixx = ixx + jj
      call fmove (q(ixx),x(ii),i)
      ii = ii + i
      ixx = ixx + i
99    continue
      jj = jj + n
100   continue
      if(ni.ge.4.and.ni.le.6) then
       call dscal(ii,-1.0d0,x,1)
          endif
      call corlsr (ibase)
      return
      end
      subroutine sorta (g,nijkl,map,dc,fc,fock)
      implicit real*8  (a-h,o-z)
      logical fock
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
      common /bufb/ nwbnwb,lnklnk,gsort(5118)
      common /mctrns/ master,isym12,isym1,itran1,itran2,ntran
     >               ,ires,maxt,ibl5,ibl54,ibl56,nbuck,nteff
     >               ,mloww,mhi,ntri,nwb,ibuck
     >               ,mark(200),nwbuck(200),ioffpr(8,8)
c
      real*8 btri
      integer mlow, nstack, iblock, mstack
      common /stak/ btri,mlow,nstack,iblock,mstack
c
c
      real*8 energy,core,potnuc,gradnt,efreez,safty,hessen
      logical lto,mcacct,mcprin
      integer iguess,nvar,iretrn,idump,isigma,iaugmx,isignh
      integer iwrnr,iblsrt,ideltr
      common /jobopt/ energy,core,potnuc,iguess,nvar,gradnt,iretrn
     +               ,idump,isigma,iaugmx,isignh,lto(10),iwrnr,iblsrt
     +               ,ideltr,efreez,safty(2),hessen,mcacct,mcprin
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
      common /lsort/ gin(510),mword
      common/craypk/iorb(340),jorb(340),korb(340),lorb(340)
      dimension g(ibl5,nbuck),nijkl(ibl5+ibl5,nbuck)
      dimension fc(nbasao,nbasao),dc(nbasao,nbasao)
      dimension map(nbasao,nirrr)
c...  set up addresses on sorted file
      call accnt('sorta   ',2)
      iad=1
      do 30 isymij=1,nirrr
      do 20 isymk=1,nirrr
      isyml=mults(isymk,isymij)
      nsymk=nsymao(isymk)
      nsyml=nsymao(isyml)
      if (isyml.gt.isymk .or. nsymk.eq.0 .or. nsyml.eq.0) goto 20
      do 10 kk=1,nsymk
c     jjj=kk-1+istart(isymk)
      map (kk-1+istart(isymk),isymij) = iad - istart(isyml)
      if (isymij.eq.1) iad = iad + kk
      if (isymij.ne.1) iad = iad + nsyml
10    continue
20    continue
30    continue
c
c...   determine base and limit triangles for this pass
      mloww=master
      mhi = min(master+nteff,lenbas)
      mtri=mhi-mloww
      ntri=(mtri-1)/nbuck+1
c...   ntri=max. no. of triangles controlled by 1 bucket
c...   nbuck=number of buckets
      do 40 ibuck=1,nbuck
      mark(ibuck)=9999999
40    nwbuck(ibuck)=0
      iblock=0
      mlow=mloww+1
c...   start loop over mainfile blocks
      call setsto(1360,0,iorb)
      call search(iblka,numa)
      call find(numa)
60    call get(gin(1),m)
      if (mword.eq.0 .or. m.eq.0) goto 90
      call find(numa)
c...   process input block
c...   unpack i,j,k,l
      call upkint
      do 80 iword=1,mword
      i = iorb(iword)
      j = jorb(iword)
      isymij=mults(symaos(i),symaos(j))
      itri=map(i,isymij)+j
      k = korb(iword)
      l = lorb(iword)
      ktri=map(k,isymij)+l
      ggg = gin(iword)
      if (i.eq.j) ggg = ggg*0.5d0
      if (k.eq.l) ggg = ggg*0.5d0
      if (itri.eq.ktri) ggg = ggg*0.5d0
      if(fock) then
      gg4=ggg*4.0d0
      fc(i,j)=fc(i,j)+gg4*dc(k,l)
      fc(k,l)=fc(k,l)+gg4*dc(i,j)
      fc(i,k)=fc(i,k)-ggg*dc(j,l)
      fc(j,l)=fc(j,l)-ggg*dc(i,k)
      fc(i,l)=fc(i,l)-ggg*dc(j,k)
      fc(j,k)=fc(j,k)-ggg*dc(i,l)
      end if
      if(itri.gt.mhi.or.itri.lt.mlow)goto 70
      ibuck=(itri-mlow)/ntri+1
      nwb=nwbuck(ibuck)+1
      g(nwb,ibuck)=ggg
      nijkl(nwb+nwb-1,ibuck)=itri
      nijkl(nwb+nwb  ,ibuck)=ktri-lentca(isymij)
      if(nwb.eq.ibl5) call mcsrto(g(1,ibuck),nijkl(1,ibuck))
      nwbuck(ibuck)=nwb
70    if(ktri.lt.mlow.or.ktri.gt.mhi.or.ktri.eq.itri)goto 80
c... triangle ktri
      ibuck=(ktri-mlow)/ntri+1
      nwb=nwbuck(ibuck)+1
      g(nwb,ibuck)=ggg
      nijkl(nwb+nwb-1,ibuck)=ktri
      nijkl(nwb+nwb  ,ibuck)=itri-lentca(isymij)
      if(nwb.eq.ibl5) call mcsrto (g(1,ibuck),nijkl(1,ibuck))
      nwbuck(ibuck)=nwb
80    continue
      goto 60
90    continue
c...   mainfile now swept
c...   clear up output
      do 100 ibuck=1,nbuck
      nwb=nwbuck(ibuck)
      if (nwb.ne.0) call mcsrto (g(1,ibuck),nijkl(1,ibuck))
100   continue
      call accnt('other   ',2)
      return
      end
      subroutine ci1o2 (zint,vec,sigma,nact,nw)
      implicit real*8  (a-h,o-z)
      common /lsort  / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension zint(*),vec(*),sigma(*)
      call unpkft(nw)
      do 10 iw=1,nw
      i1 = (iorb(iw)-1)*nact + jorb(iw)
      j1 = (korb(iw)-1)*nact + lorb(iw)
      intad = max(i1,j1)*(max(i1,j1)-1)/2+min(i1,j1)
      val = surd(ival(iw)) * zint(intad)
      sigma(ii(iw)) = sigma(ii(iw)) + val * vec(jj(iw))
      sigma(jj(iw)) = sigma(jj(iw)) + val * vec(ii(iw))
10    continue
      return
      end
      subroutine fci1o2 (zint,zint1,gam1,vec,civec,sigma,nw)
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
      common /lsort  / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511),nsurd
      dimension zint(*),zint1(*),gam1(*),vec(*),civec(*),sigma(*)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call unpkft(nw)
      do 50 iw=1,nw
      intad = ind(ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw))
      sigma(ii(iw))=sigma(ii(iw)) + surd(ival(iw)) *
     >(zint1(intad)*civec(jj(iw)) + zint(intad)*vec(jj(iw)))
      gam1(intad)=gam1(intad)+surd(ival(iw))*civec(ii(iw))*vec(jj(iw))
      intad=ind(ilifa(lorb(iw))+korb(iw),ilifa(jorb(iw))+iorb(iw))
      sigma(jj(iw))=sigma(jj(iw)) + surd(ival(iw)) *
     >(zint1(intad)*civec(ii(iw)) + zint(intad)*vec(ii(iw)))
      gam1(intad)=gam1(intad)+surd(ival(iw))*civec(jj(iw))*vec(ii(iw))
50    continue
      return
      end
      subroutine denso2 (gam,vec,nw)
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
      common /lsort  / g(511),ii(511),jj(511),ival(511),iorb(511)
     1               ,jorb(511),korb(511),lorb(511)
      common /couple/ surd(511)
      dimension gam(*),vec(*)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      call unpkft(nw)
      do 90 iw=1,nw
      intad = ind (ilifa(iorb(iw))+jorb(iw),ilifa(korb(iw))+lorb(iw) )
      gam(intad)=gam(intad)+vec(ii(iw))*vec(jj(iw))*surd(ival(iw))
c     intad = ind (ilifa(lorb(iw))+korb(iw),ilifa(jorb(iw))+iorb(iw) )
c     gam(intad)=gam(intad)+vec(ii(iw))*vec(jj(iw))*surd(ival(iw))
90    continue
      return
      end
      function mcmach (string)
      character*(*) string
c..   return machine characteristics according to key given in string
      i=-1
      if (string(1:6).eq.'vector') i=1
      if (i.eq.-1) call caserr('unknown key in routine mcmach')
      mcmach = i
      return
      end
      function patset (icga,icgb)
      implicit real*8  (a-h,o-z)
      integer pat,patset,z1,z0
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
      integer na, nb, isss, icf, npair, nstra, nstrb
      integer intoff, nint1, nint2, nstraa, nstrbb, maxaa, maxbb
      integer maxrpl, maxpar, m1
      common /detcic/ na,nb,isss,icf(32),npair(8),nstra(8),nstrb(8)
     +               ,intoff(8),nint1,nint2,nstraa,nstrbb,maxaa,maxbb
     +               ,maxrpl,maxpar,m1
c
      dimension icga(na),icgb(nb)
      data z0,z1/0,1/
      pat = z0
      jlast=0
      do 1 j1=1,na
      j=icga(j1)
      pat = ishft(pat,j-jlast)
      jlast=j
      pat=ior(pat,z1)
    1 continue
      do 2 j1=1,nb
      j=icgb(j1)+nact
      pat=ishft(pat,j-jlast)
      jlast=j
    2 pat=ior(pat,z1)
      patset=ishft(pat,64-jlast)
      return
      end
      function packft(ii,jj,ival,i,j,k,l)
      implicit real*8  (a-h,o-z)
      integer *4 ians
      integer packft,ans
      dimension ians(2)
      equivalence (ans,ians(1))
      data i5, i10 , i15  , i20    , i16  /
     1     32, 1024, 32768, 1048576, 65536/
      ians(1) = ii * i16 + jj
      ians(2) = i*i15 + j*i10 + k*i5 + l + ival*i20
      packft=ans
      return
      end
      subroutine unpkft (nw)
      implicit real*8  (a-h,o-z)
      integer *4 ig
      common /lsort/ig(2,511),ib1(511),ib2(511),ib3(511),ib4(511)
     1               ,ib5(511),ib6(511),ib7(511)
      data i5, i10 , i15  , i20    , i16  /
     1     32, 1024, 32768, 1048576, 65536/
      if (nw.eq.0) return
      do 10 i=1,nw
      ib1(i) = ig(1,i)/i16
      ib2(i) = ig(1,i) - ib1(i)*i16
      ia = ig(2,i)
      ib3(i) = ia     /i20
      ia = ia - ib3(i)*i20
      ib4(i) = ia     /i15
      ia = ia - ib4(i)*i15
      ib5(i) = ia     /i10
      ia = ia - ib5(i)*i10
      ib6(i) = ia     /i5
      ib7(i) = ia - ib6(i)*i5
 10   continue
      return
      end
      subroutine upkint
      implicit real*8  (a-h,o-z)
      common/lsort / gin(510),mword,ispc
      common/craypk/ i1(340),j1(340),k1(340),l1(340)
     * ,ipbuf(1360)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      call unpack(gin(num2e+1),lab816,ipbuf,numlab)
      call icopy (mword,ipbuf(2),4,i1,1)
      call icopy (mword,ipbuf(1),4,j1,1)
      call icopy (mword,ipbuf(4),4,k1,1)
      call icopy (mword,ipbuf(3),4,l1,1)
      end
      subroutine upcktr (gij,nij)
c...  these two routines to handle label packing on intermediate files
c...  in 4-index
      implicit real*8  (a-h,o-z)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      integer gij
      common/lsort/g(510),mword,icode
      dimension gij(*),nij(*)
      itemp=mword
      call upack2(itemp,mword,icode)
      call unpack(gij,lab1632,nij,numlabp)
      return
      end
      subroutine trout (gout,iout,mword,icode,num)
      implicit real*8  (a-h,o-z)
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer m511, num2e, num2ep, num2ejk, mvadd, mach12
      integer numlab, numlabp, numlabjk
      integer lab816, lab1632
      logical o255i
      common/atmblk/m511,num2e,num2ep,num2ejk,mvadd,mach12,
     +              numlab,numlabp,numlabjk,
     +              lab816,lab1632,o255i
c
      common /scra  / gg(340),ij(170),gm
      dimension gout(*),iout(*)
      call dcopy(mword,gout,1,gg,1)
c
      gm = pack2(mword,icode)
      call pack(gg(num2ep+1),lab1632,iout,numlabp)
      call put (gg,511,num)
      if (odebug(34)) then
       write(6,*)'trout, mword =', mword
       do i = 1, mword
       write(6,98765) i, gg(i), iout(i)
98765  format(1x,i4,f20.9,2x,i16)
       enddo
      endif
c
      return
      end
      subroutine ver_mcscfa(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/mcscfa.m,v $
     +     "/
      data revision /"$Revision: 6178 $"/
      data date /"$Date: 2010-08-10 16:57:58 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
