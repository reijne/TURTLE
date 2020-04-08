c 
c  $Author: hvd $
c  $Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
c  $Locker:  $
c  $Revision: 5774 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/anale.m,v $
c  $State: Exp $
c  
      subroutine anairr(q)
c
c  ------------------------------------------------
c    analysis of infra-red , raman and vcd effects
c    using previously calculated dipole derivatives
c    polarizability derivatives, vcd parameters and
c    cartesian force constants
c    ----------------------------------------------
c
      implicit real*8  (a-h,o-z)
      implicit integer (i-n)
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
      common/maxlen/maxq
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
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      dimension q(*)
      character*7 fnm
      character*6 snm
      data fnm,snm/"anale.m","anairr"/
c
      nat3 = nat*3
      nnc = nat3*(nat3+1)/2
      itmp = 2*nnc + 3*nat3*nat3 + 2*nat3 + 3*maxat
      i1 = igmem_alloc_inf(itmp,fnm,snm,'i1',IGMEM_DEBUG)
      i2 = i1 + nnc
      i3 = i2 + nnc
      i4 = i3 + nat3*nat3
      i5 = i4 + nat3
      i6 = i5 + nat3
      i7 = i6 + nat3*nat3
      i8 = i7 + nat3*nat3
      i9 = i8 + maxat*3
c
c  loop over isotopic substitution patterns
c
      nmv = mass_numvec()

      do imv=1,nmv
         
         if(nmv.ne.1)then
             write(iwr,100)imv
 100         format(/,1x,'Considering set of nuclear masses no.',i4,/)
         endif
         call rotcon(q(i1), q(nat+i1), q(2*nat+i1), q(3*nat+i1), nat,
     +        oprn(20),imv)
         oprn(20) = .false.
         call iranal(q(i1),q(i2),q(i3),q(i4),q(i5),q(i6),q(i7),q(i8),
     +        nat3,imv)
      enddo

      call gmem_free_inf(i1,fnm,snm,'i1')
      return
      end
      subroutine iranal(fc,a,vec,e,rm,f1,al,grd,nat3,imv)
c
c     works out infra red and raman intensities
c     includes v.c.d effects
c     picks up mp2 dipole derivatives if present
c
      implicit none
      real*8 fc,a,vec,e,rm,f1,al,grd
      real*8 az,axz,aiso,aniso,ani,axx,ax,axy,ayy,ay,ayz,azz
      real*8 conv,buff,fac1,depol,ds,dse,dsn,dst,factor,fac2,dax,dipd
      real*8 four,omega,res,pold,vcd,tfact,ten18,tens,ten5,conver
      integer nat3,imv,i,j,k,l,iat,i0,icount,ij,irpass,isec56
      integer knat,lds50,m,max,min,n,ii,iblok,isec14
      logical exist
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
c      common/funcon/bohr(8)
      real*8 bohr,charge,veloc,hbar,pi,amu,avo,debye,hartre
      common/funcon/bohr,charge,veloc,hbar,pi,amu,avo,debye,hartre
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
      character *5 comp
c
      character *1 clab,clabi,clabj
      dimension clabi(3),clabj(3),buff(3)
      dimension clab(3),fc(*),a(*),e(nat3),comp(3),
     &vec(nat3,nat3),rm(nat3),f1(*),al(nat3,nat3),grd(*)
      common/small/dipd(3,maxat*3),pold(3,3,maxat*3) ,vcd(maxat*3,6)
      common/junk/dax(3,maxat*3)
c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c
      integer ntypr
      parameter (ntypr=40)
      logical oprn,odebug,opunch
      common/prnprn/oprn(ntypr)
      common/pdebug/odebug(ntypr)
      common/ppunch/opunch(ntypr)
      data clab/'x','y','z'/
      data comp/'d/dx','d/dy','d/dz'/
c     data zero /0.0d0/
c
c     To enhance the consistency of the computed results between 
c     different parts of the code:
c
c     Now compute the conversion factor from the physical constants that
c     are used throughout the rest of the code instead of relying on 
c     some precomputed value. The old and the new value differ in the 
c     4th digit already... (hvd 20/09/2005)
c
c     data tfact /2.6436411815482d+07/
      data ten5,ten18,four/1.0d5,1.0d18,4.0d0/
      factor = (ten5/(four*pi*pi))*(avo/veloc/veloc)
      conver = (ten18*hartre)/(bohr*bohr)
      tfact  = factor*conver
c
c     look up the relevant masses but store them as m**(-1/2)!?
c
      call dr2mas(imv,rm)
c
c     read cartesian force constant matrix
c
      call rdfcm(vec,'iranal')
cj 
      if (opunch(2)) then
         call secget(isect(14),14,isec14)
         call rdedx(grd,maxat*3,isec14,ifild)
         call dr2pun(title,c,grd,vec,nat,nat*3,ipu,nat*3)
      endif
cj
c
c     project rotations and translations
c
      call prjfcm(vec,f1,al,c,rm,nat,nat*3,iwr)
c 
      do 30 i = 1 , nat3
         do 20 j = 1 , i
            fc(iky(i)+j) = vec(i,j)
 20      continue
 30   continue

c
c     do mass weighting (the masses are printed from rotcon)
c
      do 50 i = 1 , nat3
         do 40 j = 1 , i
            ij = iky(i) + j
            a(ij) = rm(i)*fc(ij)*rm(j)
 40      continue
 50   continue
c
c     ----- get normal modes and frequencies -----
c
      call ngdiag(a,vec,e,iky,nat3,nat3,1,1.0d-20)
c
c     ----- convert frequencies to inverse cm -----
c
      do 60 i = 1 , nat3
         e(i) = dsign(1.0d0,e(i))*dsqrt(dabs(tfact*e(i)))
 60   continue
c
      do 80 i = 1 , nat3
         do 70 j = 1 , nat3
            vec(j,i) = vec(j,i)*rm(j)
 70      continue
 80   continue
c
c --  formatted file for graphics - don't write
c     translational or rotational modes - the tolerance
c     of 0.1 cm-1 is empirical as projected out freqs are
c     not exactly zero
c
      icount=0
      do i = 1, 3*nat
         if(dabs(e(i)).gt.1.0d-1)then
            icount = icount + 1
            call blkvib(icount,vec(1,i),e(i))
         endif
      enddo

      if (oprn(5).and.oprn(23) ) then
         write (iwr,6020)
         max = 0
 90      min = max + 1
         max = max + 8
         if (max.gt.nat3) max = nat3
         write (iwr,6040)
         write (iwr,6200) (j,j=min,max)
         write (iwr,6040)
         write (iwr,6210) (e(j),j=min,max)
         write (iwr,6040)
         do 100 iat = 1 , nat
            i0 = 3*(iat-1)
            write (iwr,6220) iat , zaname(iat) , clab(1) ,
     +                      (vec(i0+1,j),j=min,max)
            write (iwr,6230) clab(2) , (vec(i0+2,j),j=min,max)
            write (iwr,6230) clab(3) , (vec(i0+3,j),j=min,max)
 100     continue
         if (max.lt.nat3) go to 90
      else
         write (iwr,6010) (e(j),j=1,nat3)
      end if
c
c     get dipole derivatives
c
      irpass = 0
c     mtype = 0
      call secloc(isect(50),exist,iblok)
      if (.not.exist) go to 180
      lds50 = lds(isect(50))
      if (oprn(25) .and. oprn(40)) write (iwr,6030)
 110  call rdedx(dipd,lds50,iblok,ifild)
      fac1 = charge*1.0d10
      do 130 i = 1 , 3
         do 120 j = 1 , nat3
            dipd(i,j) = dipd(i,j)*fac1
 120     continue
 130  continue
      do 150 ii = 1 , nat
         if (oprn(25) .and. oprn(40)) write (iwr,6040)
         do 140 j = 1 , 3
            i = (ii-1)*3 + j
            if (oprn(25) .and. oprn(40)) write (iwr,6050) zaname(ii), 
     +    comp(j), dipd(1,i) , dipd(2,i) , dipd(3,i)
 140     continue
 150  continue
      if (oprn(40)) then
        write (iwr,6040)
        write (iwr,6240)
        write (iwr,6060)
        write (iwr,6240)
      endif
      n = 0
      do 170 i = 1 , nat3
c
cjmht
c The old cutoff was 0.005, but this was causing problems with validating
c c2018_a in parallel, as one of the first frequencies projected out was
c -0.01. As frequencies less than a wavenumber are rarely significant
c I've raised the cutoff to 0.05
c
         if (dabs(e(i)).ge.0.05d0) then
            omega = 2.0d0*pi*veloc*dabs(e(i))
            fac1 = 1.0d8*dsqrt(hbar/(2.0d0*omega*amu))
            fac2 = 1.0d-25*avo*pi/(3.0d0*veloc*veloc*amu)
            ax = 0.0d0
            ay = 0.0d0
            az = 0.0d0
            do 160 j = 1 , nat3
               ax = ax + dipd(1,j)*vec(j,i)
               ay = ay + dipd(2,j)*vec(j,i)
               az = az + dipd(3,j)*vec(j,i)
 160        continue
            res = (ax*ax+ay*ay+az*az)*fac2
            ax = ax*fac1
            ay = ay*fac1
            az = az*fac1
            if (dabs(ax).lt.1d-11) ax = 0.0d0
            if (dabs(ay).lt.1d-11) ay = 0.0d0
            if (dabs(az).lt.1d-11) az = 0.0d0
            dax(1,i) = ax
            dax(2,i) = ay
            dax(3,i) = az
            ds = ax*ax + ay*ay + az*az
            if (oprn(40)) write (iwr,6070) e(i) , ax , ay , az , 
     +                    ds , res
            if(dabs(e(i)).gt.1.0d-1)then
               n = n + 1
               call blkir(n,res)
            endif
         end if
 170  continue
      if (oprn(40)) write (iwr,6240)
 180  if (irpass.ne.1) then
         call secloc(isect(57),exist,iblok)
         if (exist) then
            lds50 = lds(isect(57))
            if (oprn(25) .and. oprn(40)) write (iwr,6080)
            irpass = 1
            go to 110
         end if
      end if
c
c     ------------------- vcd parameters
c
      call secloc(isect(56),exist,isec56)
      if (exist) then
         call rdedx(vcd,lds(isect(56)),isec56,ifild)
         if (oprn(27). and. oprn(40)) write (iwr,6090)
         if (oprn(27). and. oprn(40)) write (iwr,6100)
         do 210 i = 1 , nat3
            if (dabs(e(i)).ge.0.02d0) then
               ax = 0.0d0
               ay = 0.0d0
               az = 0.0d0
               do 190 j = 1 , nat3
                  ax = ax + vcd(j,1)*vec(j,i)
                  ay = ay + vcd(j,2)*vec(j,i)
                  az = az + vcd(j,3)*vec(j,i)
 190           continue
               omega = 2.0d0*pi*veloc*dabs(e(i))
               fac1 = 2.0d0*dsqrt(0.5d0*hbar*omega/amu)*(debye/veloc)
               ax = ax*fac1
               ay = ay*fac1
               az = az*fac1
               if (dabs(ax).lt.1d-14) ax = 0.0d0
               if (dabs(ay).lt.1d-14) ay = 0.0d0
               if (dabs(az).lt.1d-14) az = 0.0d0
               dse = ax*dax(1,i) + ay*dax(2,i) + az*dax(3,i)
               if (oprn(27). and. oprn(40)) write (iwr,6110) e(i) , 
     +         ax , ay , az , dse
               ax = 0.0d0
               ay = 0.0d0
               az = 0.0d0
               do 200 j = 1 , nat3
                  ax = ax + vcd(j,4)*vec(j,i)
                  ay = ay + vcd(j,5)*vec(j,i)
                  az = az + vcd(j,6)*vec(j,i)
 200           continue
               ax = ax*fac1
               ay = ay*fac1
               az = az*fac1
               if (dabs(ax).lt.1d-14) ax = 0.0d0
               if (dabs(ay).lt.1d-14) ay = 0.0d0
               if (dabs(az).lt.1d-14) az = 0.0d0
               dsn = ax*dax(1,i) + ay*dax(2,i) + az*dax(3,i)
               if (oprn(27). and. oprn(40)) write (iwr,6120) 
     +         ax , ay , az , dsn
               dst = dse + dsn
               if (oprn(27). and. oprn(40)) write (iwr,6130) dst
            end if
 210     continue
      end if
c
c     ----------- polarizability derivatives , raman intensities
c
      call secloc(isect(55),exist,iblok)
      if (.not.exist) return

      call rdedx(pold,lds(isect(55)),iblok,ifild)
      if (oprn(26) .and. oprn(40)) write (iwr,6150)
      conv = bohr*bohr
      do 240 k = 1 , nat3
         knat = k
         if (oprn(26) .and. oprn(40)) write (iwr,6040)
         l=0
         do 230 i = 1 , 3
            do 220 j = 1 , 3
               pold(i,j,k) = pold(i,j,k)*conv
               l=l+1
               clabi(l) = clab(i)
               clabj(l) = clab(j)
               buff (l) = pold(i,j,k)
               if(l.lt.3) go to 220
               if (oprn(26) .and. oprn(40))
     +          write(iwr,6160) (clabi(l),clabj(l),k,buff(l),l=1,3)
                l=0
 220        continue
 230     continue
 240  continue
      if (oprn(26) .and. oprn(40)) then
        write(iwr,6160) (clabi(m),clabj(m),knat,buff(m),m=1,l)
        write (iwr,6170)
      endif
      n = 0
      do 260 i = 1 , nat3
         if (dabs(e(i)).ge.0.01d0) then
            axx = 0.0d0
            ayy = 0.0d0
            azz = 0.0d0
            axy = 0.0d0
            axz = 0.0d0
            ayz = 0.0d0
            do 250 j = 1 , nat3
               axx = axx + pold(1,1,j)*vec(j,i)
               ayy = ayy + pold(2,2,j)*vec(j,i)
               azz = azz + pold(3,3,j)*vec(j,i)
               axy = axy + pold(1,2,j)*vec(j,i)
               axz = axz + pold(1,3,j)*vec(j,i)
               ayz = ayz + pold(2,3,j)*vec(j,i)
 250        continue
            aiso = (axx+ayy+azz)/3.0d0
            aniso = 0.5d0*((axx-ayy)*(axx-ayy)+(axx-azz)*(axx-azz)
     +              +(ayy-azz)*(ayy-azz)+6.0d0*(axy*axy+axz*axz+ayz*ayz)
     +              )
            ani = dsqrt(aniso)
            tens = 45.0d0*aiso*aiso + 7.0d0*aniso
            depol = 0.0d0
            if (tens.gt.1.0d-5) then
               depol = 3.0d0*aniso/(45.0d0*aiso*aiso+4.0d0*aniso)
            end if
            if (oprn(26) .and. oprn(40)) 
     +      write (iwr,6180) e(i) , aiso , ani , tens , depol
            if(dabs(e(i)).gt.1.0d-1)then
               n = n + 1
               call blkraman(n,tens)
            endif
         end if
 260  continue

      return
 6010 format (/20x,'********************************'/
     +         20x,'harmonic frequencies (projected)'/
     +         20x,'********************************'//(1x,8f11.2))
 6020 format (/20x,'********************************************'/
     +         20x,'******    harmonic frequencies    **********'/
     +         20x,'* (translation and rotation projected out) *'/
     +         20x,'********************************************'//
     +20x,
     +'cartesians to normal mode coordinates transformation matrix'/
     +20x,
     +'==========================================================='
     +)
 6030 format (/
     + 20x,'****************************'/
     + 20x,'*  scf dipole derivatives  *'/
     + 20x,'*     (debye/angstrom)     *'/
     + 20x,'****************************'//
     + 30x,'x',15x,'y',15x,'z')
 6040 format (/)
 6050 format (5x,a8,5x,a5,3f16.8)
 6060 format (1x,'=',9x,'normal mode',16x,'transition dipole',14x,
     +        'dipole strength',6x,'intensity ='/
     +1x,'=',39x,'debyes',28x,'debyes**2',8x,'km/mole ='/
     +1x,'=',27x,'x',14x,'y',14x,'z',40x,'=')
 6070 format (' =',11x,f9.2,2x,4e15.5,f15.5,' =')
 6080 format (/
     + 20x,'****************************'/
     + 20x,'*  mp2 dipole derivatives  *'/
     + 20x,'*     (debye/angstrom)     *'/
     + 20x,'****************************'//
     + 30x,'x',15x,'y',15x,'z')
 6090 format (/20x,'*******************************'/
     +         20x,'* magnetic transition moments *'/
     +         20x,'*******************************'/)
 6100 format (/9x,'normal mode',12x,'transition magnetic dipole',5x,
     +        'rotational strength',6x,'intensity'/39x,'debyes',23x,
     +        'debyes**2',10x,'km/mole'/27x,'x',14x,'y',14x,'z')
 6110 format (1x,'electronic',2x,f9.2,4e15.5)
 6120 format (1x,'nuclear',14x,4e15.5)
 6130 format (1x,'total',61x,e15.5)
 6150 format (//20x,'****************************************'/
     +          20x,'* cartesian polarizability derivatives *'/
     +          20x,'*         units of angstrom**2         *'/
     +          20x,'* **************************************'//)
 6160 format (3(5x,a2,3x,a2,4x,i3,4x,f13.8))
 6170 format (/20x,'**********************'/
     +         20x,'*     raman data     *'/
     +         20x,'**********************'//
     +       10x,'mode',10x,'isotropic',12x,
     +        'anisotropic',10x,'intensity',12x,'depolarization ratio')
 6180 format (/5x,f9.2,4(5x,f16.8))
 6200 format (15x,8(3x,i4,4x))
 6210 format (15x,8f11.2)
 6220 format (i3,2x,a8,1x,a1,8f11.6)
 6230 format (14x,a1,8f11.6)
 6240 format (1x,100('='))
      end
      subroutine prjfco(a,ia,m,n,w)
      implicit real*8  (a-h,o-z)
      dimension a(ia,n),w(n)
c
c     Use repeated modified Gramm-Schmidt
c
      do 100 na = 1 , n
         icycle = 0
         c = ddot(m,a(1,na),1,a(1,na),1)
         c = dsqrt(c)
         if (c.gt.1.0d-10) then
            do i = 1 , m
               a(i,na) = a(i,na)/c
            enddo
         else
            call vclr(a(1,na),1,m)
         end if
 10      do nb = 1 , na - 1
            wt = ddot(m,a(1,na),1,a(1,nb),1)
            do i = 1 , m
               a(i,na) = a(i,na) - wt*a(i,nb)
            enddo
         enddo
         icycle = icycle + 1
         c = ddot(m,a(1,na),1,a(1,na),1)
         c = dsqrt(c)
         if (c.gt.10d-10) then
            do i = 1 , m
               a(i,na) = a(i,na)/c
            enddo
            if (c.le.dsqrt(0.5d0).and.icycle.le.2) goto 10
         else
            call vclr(a(1,na),1,m)
         end if
 100  continue
      return
      end
      subroutine prjfcm(fcm,f1,al,c,rm,nat,ndim,iwr)
      implicit real*8    (a-h,p-z)
      implicit logical (o)
      implicit integer (i-n)
c
c     projects out rotation and translation from cartesian force
c     constant matrix
c
c     It seems to be "better" to apply the projection to the mass
c     weighted force constant matrix. I have no idea why!
c     However this means that the force constant matrix elements are
c     scaled first, then a mass weighted projection operator is build
c     and applied, and finally the force constant matrix is scaled
c     back again. (hvd 20/09/2005)
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
      character*8 zopti,zfrozen
      common/runopt/zopti(maxat),zfrozen
      dimension wks(6)
      dimension c(3,nat)        ! nuclear coordinates
      dimension fcm(ndim,ndim)  ! force constant matrix
      dimension rm(nat*3)       ! rm(3*(i-1)+1:3)=1/sqrt(mass_i)
      dimension f1(nat*3,nat*3) ! workspace
      dimension al(nat*3,nat*3) ! workspace
      dimension cmc(3)          ! centre of mass
      dimension pmom(3)         ! moments of inertia
      dimension vv(9)           ! unit vectors of inertia
c
      if (zfrozen.eq.'yes') then
         write(iwr,1)
1        format(/' ************************************************'
     +         ,/' ** frozen nuclei disable rot/trans projection **'
     +         ,/' ************************************************')
         return
      end if
c
c     Use al to set up the 6 transitional/rotational coordinates
c     the rest of al is unused here.
c
      nat3 = nat*3
      call vclr(al,1,nat3*6)
c
c     Mass weigh the force constant matrix
c
      do i = 1, nat3
         do j = 1, nat3
            fcm(i,j) = fcm(i,j)*rm(i)*rm(j)
         enddo
      enddo
c
c     Decide whether the molecule is linear
c
      do i = 1, nat3
         rm(i)=1.0d0/(rm(i)*rm(i))
      enddo
      call mofi(nat,c,rm,cmc,pmom,vv)
      oline = omol_is_linear(pmom)
      do i = 1, nat3
         rm(i)=1.0d0/dsqrt(rm(i))
      enddo
      ntrrot = 6
      if (oline) ntrrot = 5
c
c     Mass weigh the basis for the translational and rotational 
c     coordinate space
c
      do 20 n = 1 , nat
         sqrtmi = 1.0d0/rm((n-1)*3+1)
         cx = c(1,n) - cmc(1)
         cy = c(2,n) - cmc(2)
         cz = c(3,n) - cmc(3)
         cxp = cx*vv(1) + cy*vv(2) + cz*vv(3)
         cyp = cx*vv(4) + cy*vv(5) + cz*vv(6)
         czp = cx*vv(7) + cy*vv(8) + cz*vv(9)
         al(3*n-2,1) = 1.0d0*sqrtmi
         al(3*n-1,2) = 1.0d0*sqrtmi
         al(3*n  ,3) = 1.0d0*sqrtmi
         al(3*n-2,4) = (cyp*vv(7)-czp*vv(4))*sqrtmi
         al(3*n-1,4) = (cyp*vv(8)-czp*vv(5))*sqrtmi
         al(3*n  ,4) = (cyp*vv(9)-czp*vv(6))*sqrtmi
         al(3*n-2,5) = (czp*vv(1)-cxp*vv(7))*sqrtmi
         al(3*n-1,5) = (czp*vv(2)-cxp*vv(8))*sqrtmi
         al(3*n  ,5) = (czp*vv(3)-cxp*vv(9))*sqrtmi
         al(3*n-2,6) = (cxp*vv(4)-cyp*vv(1))*sqrtmi
         al(3*n-1,6) = (cxp*vv(5)-cyp*vv(2))*sqrtmi
         al(3*n  ,6) = (cxp*vv(6)-cyp*vv(3))*sqrtmi
 20   continue
c
c     If the molecule is linear we need to remove the rotation around
c     the axis of the molecule. This is done by selecting the rotation
c     vector with the smallest norm and removing that one.
c     This approach works independent of the orientation of the molecule
c     and does not require an additional cutoff. The introduction of an
c     additional cutoff is undesirable because it might conflict with
c     the tolerance for the detection of the molecule being linear.
c
      if (oline) then
         dumx = ddot(nat3,al(1,4),1,al(1,4),1)
         idum = 4
         do i = 5, 6
            dumr = ddot(nat3,al(1,i),1,al(1,i),1)
            if (dumr.lt.dumx) then
               dumx = dumr
               idum = i
            endif
         enddo
         do n = 1, nat3
            al(n,idum) = al(n,6)
            al(n,6)    = 0.0d0
         enddo
      endif
      call prjfco(al,nat3,nat3,ntrrot,wks)
c
c     Construct the projection operator P = 1 - |A><A| in f1 using the 
c     coordinates kept in al.
c
      call vclr(f1,1,nat3*nat3)
      call mxmb(al,1,nat3,al,nat3,1,f1,1,nat3,nat3,ntrrot,nat3)
      do 40 j = 1 , nat3
         do 30 i = 1 , nat3
            f1(i,j) = -f1(i,j)
 30      continue
         f1(j,j) = f1(j,j) + 1.0d0
 40   continue
c
c     Apply the projection operator to fcm.
c     (Re-use al as scratch space...)
c
      call vclr(al,1,nat3*nat3)
      call mxmb(f1,1,nat3,fcm,1,ndim,al,1,nat3,nat3,nat3,nat3)
      call vclr(fcm,1,ndim*ndim)
      call mxmb(al,1,nat3,f1,1,nat3,fcm,1,ndim,nat3,nat3,nat3)
c
c     Now scale the force constant matrix again (sigh)
c
      do i = 1, nat3
         do j = 1, nat3
            fcm(i,j) = fcm(i,j)/(rm(i)*rm(j))
         enddo
      enddo
      return
      end
      subroutine ver_anale(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/anale.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
