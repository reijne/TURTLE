      subroutine drfsint(iatom,dsx,dsy,dsz)
c
c this is an old routine (gamess 91) with new includes ....
c sf 02-97
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
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
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
c     common/nshel_z/ex(mxprim),cs(mxprim),cp(mxprim),cd(mxprim),
c    +               cf(mxprim),cg(mxprim),
c    +       kstart(mxshel),katom(mxshel),ktype(mxshel),
c    +       kng(mxshel),kloc(mxshel),kmin(mxshel),kmax(mxshel),
c    +               nshell,non,numorb,ndumm
c

      common/intdip/pint,qint,rint,pintd,qintd,rintd,
     1t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,cx,cy,cz
c     common/junk/desp(3,maxat),
c    *pint,qint,rint,t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj
c    + ,cx,cy,cz,
c    +  pin(25),qin(25),rin(25),pd(25),qd(25),rd(25),dij(100),
c    +  ijx(100),ijy(100),ijz(100),sx(100),sy(100),sz(100)
      dimension pin(25),qin(25),rin(25),
     +  pd(25),qd(25),rd(25),dij(100),
     +  ijx(100),ijy(100),ijz(100),sx(100),sy(100),sz(100)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
c
      integer nmaxly, ntotly, icurly, isecly
      logical oprintm
      common /segm/ nmaxly,ntotly,icurly,isecly(100),
     +              oprintm
c
      common/blkin/gout(510),nword
c
      logical oianj, okanl, oident, omisc, oham, opdipd, omp2
      integer ipos1, ipos2
      common /misc/ oianj,okanl,oident,omisc,
     +            oham,opdipd,omp2,ipos1,ipos2
c
c
c
c
c     dimension sz1(*),sx1(*),sy1(*)
      dimension dsx(*), dsy(*), dsz(*)
      data sqrt3 /1.73205080756888d0/
      data dzero,done /0.0d0,1.0d0/
      data rln10 /2.30258d0/
      data zgrhf,zgvb/'grhf','gvb'/
      data m511/511/
c
c     ----- calculate derivatives of the overlap matrix -----
c
      tol = rln10*itol
      out = nprint.eq. - 3 .or. nprint.eq. - 10
      onorm = normf.ne.1 .or. normp.ne.1
      l1 = num
      l2 = (num*(num+1))/2
      l3 = num*num
      nword = 1
c
c     ----- ishell
c
      do 110 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 110
         pi = c(1,iat)
         qi = c(2,iat)
         ri = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         lit1 = lit + 1
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c
c     ----- jshell
c
c        do 100 jj = 1 , ii
         do 100 jj = 1 , nshell
            jat = katom(jj)
c           if (iat.ne.jat) then
c  sf on all atoms!
c
               pj = c(1,jat)
               qj = c(2,jat)
               rj = c(3,jat)
               j1 = kstart(jj)
               j2 = j1 + kng(jj) - 1
               ljt = ktype(jj)
               minj = kmin(jj)
               maxj = kmax(jj)
               locj = kloc(jj) - minj
               rr = (pi-pj)**2 + (qi-qj)**2 + (ri-rj)**2
               oianj = ii.eq.jj
c
c     ----- prepare indices for pairs of (i,j) functions
c
               call idxadrf(ijx,ijy,ijz,ij,
     +           mini,maxi,minj,maxj,.false.,5,
     +                    1,1)
               do 20 i = 1 , ij
                  sx(i) = dzero
                  sy(i) = dzero
                  sz(i) = dzero
 20            continue
c
c     ----- i primitive
c
               do 70 ig = i1 , i2
                  ai = ex(ig)
                  arri = ai*rr
                  axi = ai*pi
                  ayi = ai*qi
                  azi = ai*ri
                  csi = cs(ig)
                  cpi = cp(ig)
                  cdi = cd(ig)
                  cfi = cf(ig)
                  cgi = cg(ig)
c
c     ----- j primtive
c
                  do 60 jg = j1 , j2
                     aj = ex(jg)
                     aa = ai + aj
                     aa1 = done/aa
                     dum = aj*arri*aa1
                     if (dum.le.tol) then
                        fac = dexp(-dum)
                        csj = cs(jg)*fac
                        cpj = cp(jg)*fac
                        cdj = cd(jg)*fac
                        cfj = cf(jg)*fac
                        cgj = cg(jg)*fac
                        ax = (axi+aj*pj)*aa1
                        ay = (ayi+aj*qj)*aa1
                        az = (azi+aj*rj)*aa1
c
c     ----- density factor
c
                        call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                                  csj,cpj,cdj,cfj,cgj,
     +                              mini,maxi,minj,maxj,.false.,.false.,
     +                              onorm)
c
c     ----- overlap
c
                        t = dsqrt(aa1)
                        p0 = ax
                        q0 = ay
                        r0 = az
                        in = -5

                        do 40 i = 1 , lit1
                           in = in + 5
                           ni = i
                           do 30 j = 1 , ljt
                              jn = in + j
                              nj = j
c                             call vint
                              call vintdrf
                              pin(jn) = pint*t
                              qin(jn) = qint*t
                              rin(jn) = rint*t
 30                        continue
 40                     continue
                        call oneld(pin,qin,rin,pd,qd,rd,ai,lit,ljt,1,5)

c
                        do 50 i = 1 , ij
                           mx = ijx(i)
                           my = ijy(i)
                           mz = ijz(i)
                         sx(i) = sx(i) + dij(i)*pd(mx)*qin(my)*rin(mz)
                         sy(i) = sy(i) + dij(i)*pin(mx)*qd(my)*rin(mz)
                         sz(i) = sz(i) + dij(i)*pin(mx)*qin(my)*rd(mz)
 50                     continue
                     end if
c
c     ----- end of primitive loops -----
c
 60               continue
 70            continue
c
c     ----- calculate derivatives of overlap matrix -----
c
               n = 0
               do 90 i = mini , maxi
                  in = loci + i
                  do 80 j = minj , maxj
                     n = n + 1
                     jn = locj + j
c                    if (jn.le.in) then
c                       nn = iky(in) + jn
                        nn = iky(max(in,jn)) + min(in,jn)
                        dsx(nn) = dsx(nn) + sx(n)
                        dsy(nn) = dsy(nn) + sy(n)
                        dsz(nn) = dsz(nn) + sz(n)
c                       sx1(nn)=sx(n)
c                       sy1(nn)=sy(n)
c                       sz1(nn)=sz(n)
c                    end if
 80               continue
 90            continue
c           end if
c idem
 100     continue
 110  continue
c
      return
 6010 format (/40x,'lagrangian weighted density'/40x,27('*'))
 6020 format (//5x,'i',4x,'j',15x,'sx',18x,'sy',18x,'sz',18x,'lij')
 6030 format (1x,2i5,5x,3f20.8,i4)
      end
      subroutine drfdint(iatom,dipxx,dipxy,dipxz,
     1 dipyx,dipyy,dipyz,dipzx,dipzy,dipzz)
c
c     dipole moment derivatives
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
      logical out,norm
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
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
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     3    scrtchpad(150),
     4    ijx(100),ijy(100),ijz(100)
     
c     dimension dd(*)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c     dimension comp(3),dipd(3,3,maxat),dipi(3,3,maxat),
c    &  dipn(3,3,maxat)
      dimension comp(3)
      dimension dipxx(nx),dipyx(nx),dipzx(nx)
      dimension dipxy(nx),dipyy(nx),dipzy(nx)
      dimension dipxz(nx),dipyz(nx),dipzz(nx)
c
      character*5 comp
      data comp/'d/dx','d/dy','d/dz'/
      data ndim/5/
      data one /1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
c     call onepdm(dd,dd(nx+1))
      tol = 2.30258d0*itol
      out = odebug(15)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
c     do 50 n = 1 , nat
c        do 30 i = 1 , 3
c           do 20 j = 1 , 3
c              dipn(i,j,n) = 0.0d0
c              dipi(i,j,n) = 0.0d0
c              dipd(i,j,n) = 0.0d0
c20         continue
c30      continue
c        do 40 i = 1 , 3
c           dipn(i,i,n) = czan(n)
c40      continue
c50   continue
c     ----- ishell
      do 130 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 130
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 120 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,
     +           .false.,ndim,1,1)
            dxx = 0.0d0
            dyy = 0.0d0
            dzz = 0.0d0
            dxy = 0.0d0
            dxz = 0.0d0
            dyx = 0.0d0
            dyz = 0.0d0
            dzx = 0.0d0
            dzy = 0.0d0
c     ----- i primitive
            do 110 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 100 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 70 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 60 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
 60                     continue
 70                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 90 i = mini , maxi
                        in = loci + i
                        do 80 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
c                          dum = dd(nn)*dij(n)
                           dum = dij(n)
c                          dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
                           dipxx(nn) = dipxx(nn) 
     1   + dum*xdipd(nnx)*yin(ny)*zin(nz)
                           dipyy(nn) = dipyy(nn) 
     1   + dum*xin(nnx)*ydipd(ny)*zin(nz)
                           dipzz(nn) = dipzz(nn) 
     1   + dum*xin(nnx)*yin(ny)*zdipd(nz)
                           dipxy(nn) = dipxy(nn) 
     1   + dum*xdip(nnx)*yd(ny)*zin(nz)
                           dipyx(nn) = dipyx(nn) 
     1   + dum*xd(nnx)*ydip(ny)*zin(nz)
                           dipxz(nn) = dipxz(nn) 
     1   + dum*xdip(nnx)*yin(ny)*zd(nz)
                           dipzx(nn) = dipzx(nn) 
     1   + dum*xd(nnx)*yin(ny)*zdip(nz)
                           dipyz(nn) = dipyz(nn) 
     1   + dum*xin(nnx)*ydip(ny)*zd(nz)
                           dipzy(nn) = dipzy(nn) 
     1   + dum*xin(nnx)*yd(ny)*zdip(nz)
c                          dxx = dxx + dum*xdipd(nnx)*yin(ny)*zin(nz)
c                          dyy = dyy + dum*xin(nnx)*ydipd(ny)*zin(nz)
c                          dzz = dzz + dum*xin(nnx)*yin(ny)*zdipd(nz)
c                          dxy = dxy + dum*xdip(nnx)*yd(ny)*zin(nz)
c                          dyx = dyx + dum*xd(nnx)*ydip(ny)*zin(nz)
c                          dxz = dxz + dum*xdip(nnx)*yin(ny)*zd(nz)
c                          dzx = dzx + dum*xd(nnx)*yin(ny)*zdip(nz)
c                          dyz = dyz + dum*xin(nnx)*ydip(ny)*zd(nz)
c                          dzy = dzy + dum*xin(nnx)*yd(ny)*zdip(nz)
 80                     continue
 90                  continue
                  end if
 100           continue
 110        continue
c           dipi(1,1,iat) = dipi(1,1,iat) - dxx
c           dipi(2,2,iat) = dipi(2,2,iat) - dyy
c           dipi(3,3,iat) = dipi(3,3,iat) - dzz
c           dipi(1,2,iat) = dipi(1,2,iat) - dxy
c           dipi(1,3,iat) = dipi(1,3,iat) - dxz
c           dipi(2,1,iat) = dipi(2,1,iat) - dyx
c           dipi(3,1,iat) = dipi(3,1,iat) - dzx
c           dipi(2,3,iat) = dipi(2,3,iat) - dyz
c           dipi(3,2,iat) = dipi(3,2,iat) - dzy
 120     continue
 130  continue
c
c     if (out) then
c        write (iwr,6010)
c        do 150 n = 1 , nat
c           write (iwr,6020)
c           do 140 nc = 1 , 3
c              write (iwr,6030) zaname(n) , comp(nc) ,
c    +                         (dipi(nn,nc,n),nn=1,3)
c140        continue
c150     continue
c     end if
c     do 180 i = 1 , 3
c        do 170 j = 1 , 3
c           do 160 k = 1 , nat
c              dipd(i,j,k) = dipd(i,j,k) + dipn(i,j,k) + dipi(i,j,k)
c160        continue
c170     continue
c180  continue
      return
 6010 format (//35x,'integral derivative contribution'//30x,'x',15x,'y',
     +        15x,'z',/)
 6020 format (//)
 6030 format (5x,a8,5x,a5,3f16.8)
      end
      subroutine vintdrf
c
c     ----- gauss-hermite quadrature using minimum point formula -----
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
      common/intdip/pint,qint,rint,pintd,qintd,rintd,
     1t,p0,q0,r0,pi,qi,ri,pj,qj,rj,ni,nj,cx,cy,cz
      common/hermit/h(45)
      common/wermit/w(45)
      dimension min(6),max(6)
      data min /1,2,4,7,11,16/
      data max /1,3,6,10,15,21/
      data dzero /0.0d0/
c
      pint = dzero
      qint = dzero
      rint = dzero
      npts = (ni+nj-2+1)/2 + 1
      imin = min(npts)
      imax = max(npts)
      do 140 i = imin , imax
         dum = w(i)
         px = dum
         py = dum
         pz = dum
         dum = h(i)*t
         ptx = dum + p0
         pty = dum + q0
         ptz = dum + r0
         ax = ptx - pi
         ay = pty - qi
         az = ptz - ri
         bx = ptx - pj
         by = pty - qj
         bz = ptz - rj
         go to (60,50,40,30,20) , ni
 20      px = px*ax
         py = py*ay
         pz = pz*az
 30      px = px*ax
         py = py*ay
         pz = pz*az
 40      px = px*ax
         py = py*ay
         pz = pz*az
 50      px = px*ax
         py = py*ay
         pz = pz*az
 60      go to (130,120,110,100,90,80,70) , nj
 70      px = px*bx
         py = py*by
         pz = pz*bz
 80      px = px*bx
         py = py*by
         pz = pz*bz
 90      px = px*bx
         py = py*by
         pz = pz*bz
 100     px = px*bx
         py = py*by
         pz = pz*bz
 110     px = px*bx
         py = py*by
         pz = pz*bz
 120     px = px*bx
         py = py*by
         pz = pz*bz
 130     pint = pint + px
         qint = qint + py
         rint = rint + pz
 140  continue
      return
      end
      subroutine idxadrf(ijx,ijy,ijz,ij,mini,maxi,
     &minj,maxj,iandj,inc1,inc2,inc3)
      implicit real*8  (a-h,o-z)
      logical iandj
      dimension ijx(100),ijy(100),ijz(100)
c
      integer ix, iy, iz
      common /inxblk/ ix(35),iy(35),iz(35)
c
      dimension jx(20),jy(20),jz(20)
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
c        if (iandj) jmax = i
         do 30 j = minj , jmax
            ij = ij + 1
            ijx(ij) = nx + jx(j)
            ijy(ij) = ny + jy(j)
            ijz(ij) = nz + jz(j)
 30      continue
 40   continue
      return
      end
      subroutine qmdrfder(iatom,quxxx,quxxy,quxxz,
     1 quyyx,quyyy,quyyz,quzzx,quzzy,quzzz,
     2 quxyx,quxyy,quxyz,quxzx,quxzy,quxzz,
     3 quyzx,quyzy,quyzz)
c----------------------------------------------------------------
c     quadrupole derivative integrals
c----------------------------------------------------------------
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
      logical out,norm
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
c
      integer invt, nt, iliso, ilisoc, ilis48, nw196, ibl196
      integer nsymtr
      common/symtry/invt(48),nt,iliso(48),ilisoc(48),ilis48(48),
     +              nw196(6),ibl196(6),nsymtr
c
c
      real*8 timlim,ti,tx,tim,safrun,facrun,safe,dumtim
      real*8 safety,timlst,begint,elapt,cpusec
      integer isecs,jsecs
      common/timez/timlim,ti,tx,tim,
     +   safrun,facrun,safe,dumtim,safety,timlst,
     +   begint,elapt,isecs,jsecs,cpusec
c
      common/intdip/xint0,yint0,zint0,xintx,yinty,zintz,
     1t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj,origx,origy,origz
c
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      real*8 gx, gy, gz, rspace, tiny, tit, scale, ropt, vibsiz
      common/restrr/
     + gx,gy,gz,rspace(21),tiny,tit(2),scale,ropt,vibsiz
c
      logical lset,ladapt,lsym,latmol,berny,oss,ldiag,ciopt,mp2,ciforc,
     + fpres,ldens,llibry,limpt,lskip,ldenom,ignore,lcontr,lvcd,
     +lfdtrn,unit7,lgten,hfgr,bfgs,ump2,lmeth2,lcart,lmcdat,
     +lopti,lmcscf,lci,lforce,lcpf,ladp,loptor,mp2w,ordmo,ump3,rmp3,
     +opbas,odbas,ofbas,ogbas,orestrl,oatmdft,odenscfdft
c
      common/restrl/ciopt,ciforc,mp2,hfgr,bfgs,ump2,lmeth2,ump3,rmp3,
     +ordmo,mp2w,loptor,ladp,lcpf,lopti,lmcscf,lforce,lci,lcart,
     +lmcdat,lfdtrn,unit7,lcontr,lvcd,lgten,ldenom,ignore,
     +ldens,lset,ladapt,lsym,latmol,berny,llibry,limpt,fpres,oss,
     +ldiag,lskip,opbas,odbas,ofbas,ogbas,orestrl(6),oatmdft,odenscfdft
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
      integer ifilm,iblkm,mblkm,itwo,ltwo
      equivalence (ifilm,notape(1)),(iblkm,iblk(1)),(mblkm,lblk(1))
      dimension itwo(6),ltwo(6)
      equivalence (ione(7),itwo(1)),(lone(7),ltwo(1))
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
      integer len_restrl,len_restri,len_restar,len_restrr
      parameter (len_restrl=40,len_restri=1590,len_restar=700)
      parameter (len_restrr=30)
c      used: restre(util1),revise(util1),utyp21(server)
c...   lengths are not accurate
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
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
      common/small/dipd(3,3,maxat),qudd(6,3,maxat)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      common/scrtch/
     +    ptr(3,144),dtr(6,288),ftr(10,480),
     +    dij(100),
     1    xin(25),yin(25),zin(25),xd(25),yd(25),zd(25),
     2    xdip(25),ydip(25),zdip(25),
     3    xdipd(25),ydipd(25),zdipd(25),
     4    qxx(25),qyy(25),qzz(25),qxxd(25),qyyd(25),qzzd(25),
     5    ijx(100),ijy(100),ijz(100)
c     dimension dd(*)
      dimension quxxx(nx),quxxy(nx),quxxz(nx)
      dimension quyyx(nx),quyyy(nx),quyyz(nx)
      dimension quzzx(nx),quzzy(nx),quzzz(nx)
      dimension quxyx(nx),quxyy(nx),quxyz(nx)
      dimension quxzx(nx),quxzy(nx),quxzz(nx)
      dimension quyzx(nx),quyzy(nx),quyzz(nx)
       character *8 comp,grhf,closed
      dimension comp(3)
      dimension derivs(27*maxat)
      equivalence (derivs(1),dipd(1,1,1))
      data comp/'d/dx','d/dy','d/dz'/
      data grhf/'grhf'/
      data closed/'closed'/
      data ndim/5/
      data two,thalf/2.0d0,1.5d0/
      data zer,one /0.0d0,1.0d0/
c
c     calculate derivatives of the dipole moment
c
      origx = gx
      origy = gy
      origz = gz
c     call secget(isect(7),7,iblok)
c     call rdedx(dd,nx,iblok,ifild)
c     if (scftyp.ne.grhf .and. scftyp.ne.closed) then
c        call secget(isect(10),10,iblok)
c        call rdedx(dd(nx+1),nx,iblok,ifild)
c        do 20 i = 1 , nx
c           dd(i) = dd(i+nx) + dd(i)
c20      continue
c     end if
      tol = 2.30258d0*itol
      out = odebug(19)
      norm = normf.ne.1 .or. normp.ne.1
c     nuclear term
c     do 30 n = 1 , nat
c        zz = czan(n)
c        zx = (c(1,n)-gx)*zz
c        zy = (c(2,n)-gy)*zz
c        zz = (c(3,n)-gz)*zz
c        qudd(4,1,n) = thalf*zy
c        qudd(4,2,n) = thalf*zx
c        qudd(4,3,n) = zer
c        qudd(5,1,n) = thalf*zz
c        qudd(5,2,n) = zer
c        qudd(5,3,n) = thalf*zx
c        qudd(6,1,n) = zer
c        qudd(6,2,n) = zz*thalf
c        qudd(6,3,n) = zy*thalf
c        qudd(1,1,n) = two*zx
c        qudd(1,2,n) = -zy
c        qudd(1,3,n) = -zz
c        qudd(2,1,n) = -zx
c        qudd(2,2,n) = two*zy
c        qudd(2,3,n) = -zz
c        qudd(3,1,n) = -zx
c        qudd(3,2,n) = -zy
c        qudd(3,3,n) = two*zz
c30   continue
c     ----- ishell
      do 110 ii = 1 , nshell
         iat = katom(ii)
         if (iat .ne. iatom) goto 110
         xi = c(1,iat)
         yi = c(2,iat)
         zi = c(3,iat)
         i1 = kstart(ii)
         i2 = i1 + kng(ii) - 1
         lit = ktype(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         loci = kloc(ii) - mini
c     ----- jshell
         do 100 jj = 1 , nshell
            jat = katom(jj)
            xj = c(1,jat)
            yj = c(2,jat)
            zj = c(3,jat)
            j1 = kstart(jj)
            j2 = j1 + kng(jj) - 1
            ljt = ktype(jj)
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
c
c           nroots = (lit+ljt+1)/2
c
            rr = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
c     ----- prepare indices for pairs of (i,j) functions
            call indexa(ijx,ijy,ijz,ij,mini,maxi,minj,maxj,
     +                  .false.,ndim,1,1)
            qxxx = 0.0d0
            qxxy = 0.0d0
            qxxz = 0.0d0
            qyyx = 0.0d0
            qyyy = 0.0d0
            qyyz = 0.0d0
            qzzx = 0.0d0
            qzzy = 0.0d0
            qzzz = 0.0d0
            qxyx = 0.0d0
            qxyy = 0.0d0
            qxyz = 0.0d0
            qxzx = 0.0d0
            qxzy = 0.0d0
            qxzz = 0.0d0
            qyzx = 0.0d0
            qyzy = 0.0d0
            qyzz = 0.0d0
c     ----- i primitive
            do 90 ig = i1 , i2
               ai = ex(ig)
               arri = ai*rr
               axi = ai*xi
               ayi = ai*yi
               azi = ai*zi
               csi = cs(ig)
               cpi = cp(ig)
               cdi = cd(ig)
               cfi = cf(ig)
               cgi = cg(ig)
c     ----- j primtive
               do 80 jg = j1 , j2
                  aj = ex(jg)
                  aa = ai + aj
                  aainv = one/aa
                  dum = aj*arri*aainv
                  if (dum.le.tol) then
                     fac = dexp(-dum)
                     csj = cs(jg)*fac
                     cpj = cp(jg)*fac
                     cdj = cd(jg)*fac
                     cfj = cf(jg)*fac
                     cgj = cg(jg)*fac
                     ax = (axi+aj*xj)*aainv
                     ay = (ayi+aj*yj)*aainv
                     az = (azi+aj*zj)*aainv
c     ----- density factor
                     call denfan(dij,csi,cpi,cdi,cfi,cgi,
     +                               csj,cpj,cdj,cfj,cgj,
     +                           mini,maxi,minj,maxj,.false.,.false.,
     +                           norm)
c     ----- overlap
                     t = dsqrt(aa)
                     tinv = one/t
                     x0 = ax
                     y0 = ay
                     z0 = az
                     lit1 = lit + 1
                     in = -ndim
                     do 50 i = 1 , lit1
                        in = in + ndim
                        ni = i
                        do 40 j = 1 , ljt
                           jn = in + j
                           nj = j
                           call dmsint()
                           xin(jn) = xint0*tinv
                           yin(jn) = yint0*tinv
                           zin(jn) = zint0*tinv
                           xdip(jn) = xintx*tinv
                           ydip(jn) = yinty*tinv
                           zdip(jn) = zintz*tinv
                           call qmsint()
                           qxx(jn) = xintx*tinv
                           qyy(jn) = yinty*tinv
                           qzz(jn) = zintz*tinv
 40                     continue
 50                  continue
                     call oneld(xin,yin,zin,xd,yd,zd,ai,lit,ljt,1,ndim)
                     call oneld(xdip,ydip,zdip,xdipd,ydipd,zdipd,ai,lit,
     +                          ljt,1,ndim)
                     call oneld(qxx,qyy,qzz,qxxd,qyyd,qzzd,ai,lit,ljt,1,
     +                          ndim)
c     ----- calculate derivatives of dipole matrix -----
                     n = 0
                     do 70 i = mini , maxi
                        in = loci + i
                        do 60 j = minj , maxj
                           n = n + 1
                           jn = locj + j
                           nn = min(in,jn) + iky(max(in,jn))
c                          dum = dd(nn)*dij(n)
                           dum = dij(n)
c                          dum = dum + dum
                           nnx = ijx(n)
                           ny = ijy(n)
                           nz = ijz(n)
c
                           quxxx(nn) = quxxx(nn)
     1                      + dum*qxxd(nnx)*yin(ny)*zin(nz)
c
c                          qxxx = qxxx + dum*qxxd(nnx)*yin(ny)*zin(nz)
                           quxxy(nn) = quxxy(nn)
     1                      + dum*qxx(nnx)*yd(ny)*zin(nz)
c                          qxxy = qxxy + dum*qxx(nnx)*yd(ny)*zin(nz)
                           quxxz(nn) = quxxz(nn)
     1                      + dum*qxx(nnx)*yin(ny)*zd(nz)
c                          qxxz = qxxz + dum*qxx(nnx)*yin(ny)*zd(nz)
c
c
                           quyyx(nn) = quyyx(nn)
     1                      + dum*xd(nnx)*qyy(ny)*zin(nz)
c                          qyyx = qyyx + dum*xd(nnx)*qyy(ny)*zin(nz)
                           quyyy(nn) = quyyy(nn)
     1                      + dum*xin(nnx)*qyyd(ny)*zin(nz)
c                          qyyy = qyyy + dum*xin(nnx)*qyyd(ny)*zin(nz)
                           quyyz(nn) = quyyz(nn)
     1                      + dum*xin(nnx)*qyy(ny)*zd(nz)
c                          qyyz = qyyz + dum*xin(nnx)*qyy(ny)*zd(nz)
c
c
                           quzzx(nn) = quzzx(nn)
     1                      + dum*xd(nnx)*yin(ny)*qzz(nz)
c                          qzzx = qzzx + dum*xd(nnx)*yin(ny)*qzz(nz)
                           quzzy(nn) = quzzy(nn)
     1                      + dum*xin(nnx)*yd(ny)*qzz(nz)
c                          qzzy = qzzy + dum*xin(nnx)*yd(ny)*qzz(nz)
                           quzzz(nn) = quzzz(nn)
     1                      + dum*xin(nnx)*yin(ny)*qzzd(nz)
c                          qzzz = qzzz + dum*xin(nnx)*yin(ny)*qzzd(nz)
c
c
                           quxyx(nn) = quxyx(nn)
     1                      + dum*xdipd(nnx)*ydip(ny)*zin(nz)
c                          qxyx = qxyx + dum*xdipd(nnx)*ydip(ny)*zin(nz)
                           quxyy(nn) = quxyy(nn)
     1                      + dum*xdip(nnx)*ydipd(ny)*zin(nz)
c                          qxyy = qxyy + dum*xdip(nnx)*ydipd(ny)*zin(nz)
                           quxyz(nn) = quxyz(nn)
     1                      + dum*xdip(nnx)*ydip(ny)*zd(nz)
c                          qxyz = qxyz + dum*xdip(nnx)*ydip(ny)*zd(nz)
c
c
                           quxzx(nn) = quxzx(nn)
     1                      + dum*xdipd(nnx)*yin(ny)*zdip(nz)
c                          qxzx = qxzx + dum*xdipd(nnx)*yin(ny)*zdip(nz)
                           quxzy(nn) = quxzy(nn)
     1                      + dum*xdip(nnx)*yd(ny)*zdip(nz)
c                          qxzy = qxzy + dum*xdip(nnx)*yd(ny)*zdip(nz)
                           quxzz(nn) = quxzz(nn)
     1                      + dum*xdip(nnx)*yin(ny)*zdipd(nz)
c                          qxzz = qxzz + dum*xdip(nnx)*yin(ny)*zdipd(nz)
c
c
                           quyzx(nn) = quyzx(nn)
     1                      + dum*xd(nnx)*ydip(ny)*zdip(nz)
c                          qyzx = qyzx + dum*xd(nnx)*ydip(ny)*zdip(nz)
                           quyzy(nn) = quyzy(nn)
     1                      + dum*xin(nnx)*ydipd(ny)*zdip(nz)
c                          qyzy = qyzy + dum*xin(nnx)*ydipd(ny)*zdip(nz)
                           quyzz(nn) = quyzz(nn)
     1                      + dum*xin(nnx)*ydip(ny)*zdipd(nz)
c                          qyzz = qyzz + dum*xin(nnx)*ydip(ny)*zdipd(nz)
 60                     continue
 70                  continue
                  end if
 80            continue
 90         continue
c           qudd(1,1,iat) = qudd(1,1,iat) - 0.5d0*(qxxx+qxxx-qyyx-qzzx)
c           qudd(2,2,iat) = qudd(2,2,iat) - 0.5d0*(qyyy+qyyy-qxxy-qzzy)
c           qudd(3,3,iat) = qudd(3,3,iat) - 0.5d0*(qzzz+qzzz-qxxz-qyyz)
c           qudd(1,2,iat) = qudd(1,2,iat) - 0.5d0*(qxxy+qxxy-qzzy-qyyy)
c           qudd(1,3,iat) = qudd(1,3,iat) - 0.5d0*(qxxz+qxxz-qyyz-qzzz)
c           qudd(2,1,iat) = qudd(2,1,iat) - 0.5d0*(qyyx+qyyx-qxxx-qzzx)
c           qudd(3,1,iat) = qudd(3,1,iat) - 0.5d0*(qzzx+qzzx-qxxx-qyyx)
c           qudd(2,3,iat) = qudd(2,3,iat) - 0.5d0*(qyyz+qyyz-qxxz-qzzz)
c           qudd(3,2,iat) = qudd(3,2,iat) - 0.5d0*(qzzy+qzzy-qxxy-qyyy)
c           qudd(4,1,iat) = qudd(4,1,iat) - 1.5d0*qxyx
c           qudd(4,2,iat) = qudd(4,2,iat) - 1.5d0*qxyy
c           qudd(4,3,iat) = qudd(4,3,iat) - 1.5d0*qxyz
c           qudd(5,1,iat) = qudd(5,1,iat) - 1.5d0*qxzx
c           qudd(5,2,iat) = qudd(5,2,iat) - 1.5d0*qxzy
c           qudd(5,3,iat) = qudd(5,3,iat) - 1.5d0*qxzz
c           qudd(6,1,iat) = qudd(6,1,iat) - 1.5d0*qyzx
c           qudd(6,2,iat) = qudd(6,2,iat) - 1.5d0*qyzy
c           qudd(6,3,iat) = qudd(6,3,iat) - 1.5d0*qyzz
 100     continue
 110  continue
c
c     if (out) then
c        write (iwr,6010)
c        do 130 n = 1 , nat
c           write (iwr,6020)
c           do 120 nc = 1 , 3
c              write (iwr,6030) zaname(n) , comp(nc) ,
c    +                         (qudd(nn,nc,n),nn=1,6)
c120        continue
c130     continue
c     end if
c     term from differentiating orbitals
c     nplus1 = nocca + 1
c     iposs = iochf(15)
c     icomp = 24
c     ioff = 1
c     i1 = nx + 1
c     i2 = i1 + nx
c     i3 = i2 + nx
c     lennew = iky(ncoorb+1)
c     iblll = lensec(lennew)
c     do 140 i = 1 , 3
c        ioff = ioff + nx
c        icomp = icomp + 1
c        call vclr(dd(ioff),1,lennew)
c        if (itwo(i).eq.0) then
c           write (iwr,6040) i
c        else
c           call secget(isect(icomp),icomp,iblok)
c           call rdedx(dd(ioff),lennew,iblok,ifild)
c        end if
c140  continue
c     do 160 n = 1 , nat
c        do 150 nc = 1 , 3
cIF(secd_parallel)
c          call fetch(dd,lennew,'pdens',(n-1)*3+nc)
cELSE
c           call rdedx(dd,lennew,iposs,ifockf)
cENDIF
c           iposs = iposs + iblll
c           dumx = tracep(dd(i1),dd,ncoorb)
c           dumy = tracep(dd(i2),dd,ncoorb)
c           dumz = tracep(dd(i3),dd,ncoorb)
c           qudd(1,nc,n) = qudd(1,nc,n) - dumx
c           qudd(2,nc,n) = qudd(2,nc,n) - dumy
c           qudd(3,nc,n) = qudd(3,nc,n) - dumz
c150     continue
c160  continue
c     icomp = 27
c     ioff = 1
c     do 170 i = 4 , 6
c        icomp = icomp + 1
c        ioff = ioff + nx
c        call vclr(dd(ioff),1,lennew)
c        if (itwo(i).eq.0) then
c           write (iwr,6040) i
c        else
c           call secget(isect(icomp),icomp,iblok)
c           call rdedx(dd(ioff),lennew,iblok,ifild)
c        end if
c170  continue
c     iposs = iochf(15)
c     do 190 n = 1 , nat
c        do 180 nc = 1 , 3
cIF(secd_parallel)
csecd
c           call fetch(dd,lennew,'pdens',(n-1)*3+nc)
cELSE
c           call rdedx(dd,lennew,iposs,ifockf)
cENDIF
c           iposs = iposs + iblll
c           qudd(4,nc,n) = qudd(4,nc,n) - tracep(dd,dd(i1),ncoorb)
c           qudd(5,nc,n) = qudd(5,nc,n) - tracep(dd,dd(i2),ncoorb)
c           qudd(6,nc,n) = qudd(6,nc,n) - tracep(dd,dd(i3),ncoorb)
c180     continue
c190  continue
c     call qmdsym(qudd,dd,dd(nw196(5)+1),nat,nshell)
c     if (oprn(40)) then
c      if (oprn(25)) write (iwr,6010)
c      do 210 n = 1 , nat
c         if (oprn(25)) write (iwr,6020)
c         do 200 nc = 1 , 3
c            if (oprn(25)) write (iwr,6030) zaname(n) , comp(nc) , 
c    +      (qudd(nn,nc,n),nn=1,6)
c200      continue
c210   continue
c     endif
      call timit(3)
c     lenblk = lensec(27*maxat)
c     call secput(isect(50),50,lenblk,iblok)
c     call revind
c     lds(isect(50)) = 27*maxat
c     call wrt3(derivs,lds(isect(50)),iblok,ifild)
c     call revise
c     call clredx
      return
 6010 format (/
     +   20x,'*************************************'/
     +   20x,'*   quadrupole moment derivatives   *'/
     +   20x,'*           (atomic units)          *'/
     +   20x,'*************************************'//
     + 29x,'xx',14x,'yy',14x,'zz',14x,'xy',14x,'xz',14x,'yz'/)
 6020 format (/)
 6030 format (5x,a8,2x,a8,6f16.8)
 6040 format (//10x,'quadrupole integrals missing for component',i6)
      end
