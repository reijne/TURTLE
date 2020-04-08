c
c
c  $Author: jmht $
c  $Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
c  $Locker:  $
c  $Revision: 6176 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd5.m,v $
c  $State: Exp $
c
c ******************************************************
c ******************************************************
c             =   adler = newmrd5    =
c ******************************************************
c ******************************************************
c
      subroutine adler(q,energy,odebug,oprinsym)
c
c ===================================================================
c ==                                                               ==
c ==   direct-/semidirect mrd-ci - programm                        ==
c ==   ====================================                        ==
c ==                                                               ==
c ==   building on the mrd-ci program of buenker+peyerimhoff       ==
c ==                                                               ==
c ===================================================================
c
c
c to be done for this part:
c
c       - improve the davidson
c       - clean-up of variables
c       - implement the bk-methode
c
      integer nbak
      real*8 cpulft
      real*8 q, energy
      logical odebug, ofirst, orescue,oprinsym
      logical ovect, o2e
      dimension q(*)
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
c     unit common from selection etc
      integer ntape, mtape, mdisk, idelj, ltype
      integer linf,  ntable,  kfile, kclab, ntype
      integer mstvt, nf01, nf62, nhead, nf99, nf78, nf11
      common /ftap/ ntape, mtape, mdisk,
     .              idelj, ltype, linf,  ntable,  kfile,
     .              kclab, ntype, mstvt, nf01, nf62,
     +              nhead, nf99, nf78, nf11
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
cvp
cvp twoe contains the 2-electron integrals, ston is a record with 
cvp  2-electron integrals from ft31
      integer nteint, nidmax, mxa
c nteint : # 2-electron integrals
c     parameter( nteint = 3500001)
c nidmax : dimension of the integral records of stoney
      parameter (nidmax =    2000)
      parameter (mxa    =  nidmax)
c
      integer ndeke
c nomax  : max. # of mo's
      parameter (nomax  =     256)
c     parameter (ndeke = nomax*(nomax+1)/2+1)
      integer noeint
      real*8 core
      common /roeint/ core
      integer nform
      common /cffor/ nform
c
c mdi    : max. dimension of the hamilton matrix
      integer mdi
c     parameter (mdi    =  1000001)
c
c nedim  : dimension of emat, which contains the sga-matrices
      integer nedim
c     parameter (nedim  = 2000000)
c
c iotm   : field length for the selecte configurations
      integer iotm
c     parameter (iotm   = 2000000)
c
c maximum number of roots=mx in diagonalistion
c
      integer mx,mx1,mx2,mx21,lspace
      parameter (mx=256, mx1=(mx*mx+mx) / 2)
c
      parameter (maxref=256)
c
      logical oextrap
c
      integer nroot, ifirst, ilast, ntch, kprin
      integer ndeci, icode,konf ,keps,iggey
      integer istart,nrootz
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                 ndeci, icode,konf ,keps,iggey,
     .                 istart,nrootz
      real*8 esav36
      integer ixt,i0,mxrootn
      parameter (mxrootn = 50)
      common /per/ ixt,i0,esav36(4,mxrootn)
c
      real*8 egey1,trash,tdel
      integer mxv,mxv1,mxv2,mxv21,lspacev
      parameter (mxv=mxrootn*mxrootn)
      parameter (mxv1=(mxv*mxv+mxv) / 2)
c
      common /parkin/egey1,trash,tdel
c
      integer ici1, ici2, ici3, ici4, ici5, ici6
      logical oprint
      logical onteint,omdi,oiotm,onedim
      character*10 charwall
c
       write (iwr,12) cpulft(1) ,charwall()
c
c nomax  : max. # of mo's (now set to num)
c     nomax = num
      ndeke = num*(num+1)/2+1
c noeint : max # of 1-electron integrals
      noeint =  ndeke
      nomax2 = nomax + nomax
c
      call adler_dat
c
      call openfl
c --- printing files header 
      call header
c
      call readinad(nteint,mdi,iotm,nedim,
     +             onteint,omdi,oiotm,onedim,oprint,odebug)
c
c     now allocate integral arrays dynamically
c     determine space
c
      ofirst = .false.
      orescue = .false.
      nteintf = 0
      iotmf = 0
      mdif = 0
      ierr = 0
c
c     first, pick up key parameters
c
      mvect = mxrootn * maxref
      ivect = igmem_alloc(mvect)
      call prepref(nteinti, mdii, iotmi,q(ivect),mvect)
      call gmem_free(ivect)
c     fix values from data input?
      if(.not.onteint) then
       nteint = nteinti
      endif
      if(.not.omdi) then
       mdi = mdii
      endif
      if(.not.oiotm) then
       iotm = iotmi
      endif
c
      write(iwr,500)  nteint, mdi, iotm
c
      nav = lenwrd()
      mx2 = mx * mx
      mx21 = max(mx2,mx1)
      lspace = mx1 + mx2 + mx21
c
      mxv2 = mxv * mxv
      mxv21 = max(mxv2,mxv1)
      lspacev = mxv1 + mxv2 + mxv21
c
      lspace = max(lspace,lspacev)
c
c     allow for space in aftci
      needv = max(lspace,mvect)
c
c     do the default values fit into this allocation?
c     now attempt to determine exact memory allocations
c
1000  ofirst = .not.ofirst
      ovect = .false.
      o2e = .false.
c
      if (ofirst) then
c
c     first determine total space available
c
       ibase = igmem_alloc_all(ntotal)
c
       write(iwr,1010)
c
       call adler_all(num,ilast,nteint,nedim,mdi,iotm,
     +        nteint0,iotm0,mdi0,mdi0v,nedim0,integc,mdic,need,
     +        oprint)
c
       if (need .gt. ntotal) then
c
c     problem with memory
c     which is the real problem here?
c     must reduce memory requirement
c
       nteints = nteint0
       call adler_reduce(ofirst,integc,mdic,
     + need,ntotal,iword,
     + nteints,nteint0,o2e,
     + ilast,mdi,mdi0,mdi0v,ovect,oprint)
c
       endif
c
      else
c
       call adler_all(num,ilast,nteintf,nedim,mdif,iotmf,
     +     nteint0,iotm0,mdi0,mdi0v,nedim0,integc,mdic,need,
     +     oprint)
       write(iwr,1020)
c
       ibase = igmem_alloc_all(ntotal)
      endif
c
c twoe
 100  call adler_alloc(ibase,num,
     +     nteint0,iotm0,mdi0,mdi0v,nedim0,
     +     iston,ipey,icoul,iexc,ipey0,icoul0,iexc0,
     +     ici1,ici2,ici3,ici4,ici5,ici6,iemat,ivect,
     +     iottr,iotnew,iot0,maindf,jconb,
     +     idiag,idiagv,idiag2,idiag3,idiag3v,need,oprint)
c
      if (need. gt. ntotal) then
       if(.not.ofirst) then
        iword = need-ntotal 
        write(iwr,27) need, ntotal, iword
c
c      try to rescue this by:
c      1. further reducing nedim0
c      2. by multipassing the davidson
c      3. reducing the no. of integrals in memory, or 
c
c      iword words of memory needs to be clawed back from need ...
c
       orescue = .true.
c      1. reduce nedim0
        if(.not.onedim) then
         nedim_rescue = nedim0 / 2
         if (iword.le.nedim_rescue) then
          nedim0 = nedim0 / 2
          go to 100
         endif
         nedim0 = nedim0 / 2
         ntotal2 = ntotal + nedim_rescue
        else
         ntotal2 = ntotal
        endif
c
        nteints = nteint0
c
        call adler_reduce(ofirst,integc,mdic,
     +  need,ntotal2,iword,
     +  nteints,nteint0,o2e,
     +  ilast,mdif,mdi0,mdi0v,ovect,oprint)
c
        if(iword.gt.0) then
         iword = need-ntotal 
         write(iwr,28) need, ntotal, iword
        call caserr(
     +  'memory exhausted - must allocate more (see above)')
       endif
c
       go to 100
c
       else
c
c     insufficient memory available .. reduce parameters and
c     try again
c
        ierr = ierr + 1
        if (ierr.gt.3) then
         call caserr(
     +   'memory allocation in MRDCI in error - allocate more memory')
        endif
c       integs = nteint0 / nav
        integs = nteint0 
        mdit = 2*mdi0 + 4 * mdi0v
        iotc = iotm0 + 2 * (iotm0 /nav)
        write(iwr,13) integc, integs, mdic  , mdit,
     +                iotc  , iotc, nedim , nedim0
c
        go to 100
       endif
c
      endif
c
c     now free memory and re-allocate using revised parameters
c
      call gmem_free(ibase)
c
c     allocate required memory
c
      ibase = igmem_alloc(need)
      call adler_alloc(ibase,num,
     +     nteint0,iotm0,mdi0,mdi0v,nedim0,
     +     iston,ipey,icoul,iexc,ipey0,icoul0,iexc0,
     +     ici1,ici2,ici3,ici4,ici5,ici6,iemat,ivect,
     +     iottr,iotnew,iot0,maindf,jconb,
     +     idiag,idiagv,idiag2,idiag3,idiag3v,need,oprint)
c
      last =   need 
c
      call rumpxs(q(ibase),nteint0,q(iston),nidmax,
     +   q(ipey), q(icoul), q(iexc), ndeke, core,
     +   q(ipey0), q(icoul0), q(iexc0), noeint,
     +   q(ici3), mdi0v,
     +   q(iemat),nedim0,q(ivect),mvect,
     +   q(iottr),q(iotnew),q(iot0),iotm0,
     +   q(maindf),mx,q(jconb),nomax2,odebug,
     +   ofirst,nteintf,iotmf,mdif)
c
      if(ofirst) then
c
c     reset ioint
      ioint = 0
       if(oiotm )then
        iotmf = iotm
       endif
       if(omdi) then
        mdif = mdi
       endif
      call gmem_free(ibase)
      go to 1000
      else
c
c     integs = nteint0 / nav
      integs = nteint0 
      iotc = iotm0 + 2 * (iotm0 /nav)
      mdit = 2*mdi0 + 4 * mdi0v
      write(iwr,14) integs, mdit, iotc, nedim0, needv,
     +              need, last
c
      endif
c
c
c --- diagonalisation
c
      write(iwr,23) trash,cpulft(1) ,charwall()
      call diag (q(ibase), nteint0,
     +  q(ipey), q(icoul), q(iexc), ndeke, core,
     +  q(ici1), q(ici2), mdi0,
     +  q(ici3), q(ici4), q(ici5), q(ici6), mdi0v, 
     +  q(iemat), nedim0, q(iottr),q(iotnew),iotm0,
     +  q(maindf),mx,q(jconb),nomax2,
     +  q(idiag),q(idiag2),q(idiagv),mxv2,mxv1,mxv,odebug,
     +  q(idiag3),q(idiag3v),mx,mx1,mx2)
c
      write(iwr,24) trash,cpulft(1) ,charwall()
      nbak=0
      ivect = idiag
      call aftci(q(ici1),mdi0,q(ivect),mvect,
     +           q(iottr),q(iot0),iotm0,nbak,
     +           energy,odebug,oprinsym)
      write(iwr,66)
*      write(iwr,*) 'out adler after aftci'

c --- new choice of configurations
*      write(iwr,*) 'out adler before extp'
      call extp(nbak,oextrap)
      if (nbak.gt.0) then
* diagonalising the smaller problem
* compressing the needed fields
* start vectors, hamiltonmatrix (for semidirect)
* iot field etc.
         write(iwr,25) (trash+tdel),cpulft(1) ,charwall()
         call foxy(q(ibase),nteint0,q(iston),nidmax,
     +        q(ipey), q(icoul), q(iexc), ndeke, core,
     +        q(ipey0), q(icoul0), q(iexc0), noeint,
     +        q(ici1), q(ici2), mdi0, q(ici3), mdi0v,
     +        q(iemat),nedim0,q(iottr),q(iotnew),q(iot0),iotm0,
     +        q(maindf),mx,q(jconb),nomax2,odebug)
*
* calls the diagonalisation (as in the first pass)
* of the smaller problem
*
         call diag2(q(ibase), nteint0,
     +   q(ipey), q(icoul), q(iexc), ndeke, core,
     +   q(ici1), q(ici2), mdi0,
     +   q(ici3), q(ici4), q(ici5), q(ici6), mdi0v,
     +   q(iemat),nedim0,q(iottr),q(iotnew),iotm0,
     +   q(maindf),mx,q(jconb),nomax2,
     +   q(idiag),q(idiag2),q(idiagv),mxv1,mxv,odebug)
c
*      write(iwr,*) 'out adler after diag2'
         write(iwr,26) (trash+tdel),cpulft(1) ,charwall()
c
         call aftci(q(ici1),mdi0,q(ivect),mvect,
     +              q(iottr),q(iot0),iotm0,nbak,
     +              energy,odebug,oprinsym)
c
         write(iwr,66)
* extrapolation and printing the results
         call outp(oextrap,odebug)
*      write(iwr,*) 'out adler after outp'
      end if
c
      if (.not.oextrap) then
c     append null energy block to ft36 for analysis modules
       call vclr(esav36,1,mxrootn*4)
       write(iput) esav36
      endif
c
      call closbf3
c
       write(iwr,22)cpulft(1) ,charwall()
c
c      now free memory
c
       call gmem_free(ibase)
c
 12   format(/5x,'***  start of semi-direct MRD-CI module at ',
     +       f10.2,' seconds',a10,' wall'/)
 500  format(/
     + ' *** Key parameters in determining memory'/
     + ' *** nteint, mdi, iotm = ', 3i10)
1010  format(/
     + 1x,'***************************************************'/
     + 1x,'* First allocate sufficient memory to enable real *'/
     + 1x,'* memory allocations to be determined             *'/
     + 1x,'***************************************************'/)
1020  format(/
     + 1x,'**********************************************************'/
     + 1x,'* now derive memory for 2e-integrals, configurations and *'/
     + 1x,'* hamiltonian matrix, assuming all integrals in memory   *'/
     + 1x,'* and vectors derived in single-pass mode                *'/
     + 1x,'**********************************************************'/)
 27   format(/
     + 1x,'********************************************************'/
     + 1x,'* Insufficient memory available for maximum efficiency *'/
     + 1x,'* memory required = ', i10,' R*8 words               *'/
     + 1x,'* memory available = ', i9,' R*8 words               *'/
     + 1x,'* Increase allocated memory by ',i8,' R*8 words      *'/
     + 1x,'* for all integrals in memory and single-pass mode     *'/
     + 1x,'********************************************************'/
     + 1x,'* Attempt to determine memory allowing both multi-pass *'/
     + 1x,'* and non-resident integrals ++ LESS efficient method  *'/
     + 1x,'********************************************************'/)
 28   format(/
     + 1x,'********************************************************'/
     + 1x,'* insufficient memory available                        *'/
     + 1x,'* memory required = ', i10,' R*8 words               *'/
     + 1x,'* memory available = ', i9,' R*8 words               *'/
     + 1x,'* Increase allocated memory by ',i8,' R*8 words      *'/
     + 1x,'********************************************************'/)
 13   format(
     + 5x, '+++++++++++++++++++++++++++++++++++++++++++++++'/
     + 5x, '+ Insufficient memory for default allocations +'/
     + 5x, '+ parameters reduced as follows               +'/
     + 5x, '+++++++++++++++++++++++++++++++++++++++++++++++'/
     + 5x, '+   integc from ', i8,'  to ', i8,'         +'/
     + 5x, '+     mdit from ', i8,'  to ', i8,'         +'/
     + 5x, '+    iotm  from ', i8,'  to ', i8,'         +'/
     + 5x, '+   nedim  from ', i8,'  to ', i8,'         +'/
     + 5x, '+++++++++++++++++++++++++++++++++++++++++++++++'/)
 14   format(/
     + 5x, '++++++++++++++++++++++++++++'/
     + 5x, '+ Final memory allocations +'/
     + 5x, '++++++++++++++++++++++++++++'/
     + 5x, '+    integs', i10, ' R*8  +'/
     + 5x, '+      mdiT', i10, ' R*8  +'/
     + 5x, '+      iotm', i10, ' R*8  +'/
     + 5x, '+     nedim', i10, ' R*8  +'/
     + 5x, '+      diag', i10, ' R*8  +'/
     + 5x, '++++++++++++++++++++++++++++'/
     + 5x, '+ allocated', i10, ' R*8  +'/
     + 5x, '+ required ', i10, ' R*8  +'/
     + 5x, '++++++++++++++++++++++++++++'/)
 23   format(1x,
     + '*** commence lower threshold (T=',f5.1,' uh)',
     + ' calculation at ', f10.2,' seconds',a10,' wall'/)
 24   format(/1x,
     + '*** lower threshold (T=',f5.1,' uh)',
     + ' calculation complete at ',f10.2,' seconds',a10,' wall'/)
 26   format(/1x,
     + '*** higher threshold (T=',f5.1,' uh)',
     + ' calculation complete at ',f10.2,' seconds',a10,' wall'/)
 66   format(/1x,68('=')/)
 25   format(1x,
     + '*** commence higher threshold (T=',f5.1,' uh)',
     + ' calculation at ',f10.2,' seconds',a10,' wall'/)
 22   format(/1x,'***  end of semi-direct MRD-CI module at ',
     +     f10.2,' seconds',a10,' wall')
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adler_all(num,ilast,nteint,nedim,mdi,iotm,
     +        nteint0,iotm0,mdi0,mdi0v,nedim0,integc,mdic,need2,
     +        oprint)
c
      implicit none
      integer num,ilast,nteint,nedim,mdi,iotm
      integer nteint0,iotm0,mdi0,mdi0v,nedim0,integc,mdic,need
      logical oprint
c
      integer nidmax
c
c nidmax : dimension of the integral records of stoney
      parameter (nidmax =    2000)
c
      integer mx,mx1,maxref,mxrootn,nomax
      parameter (mx=256, mx1=(mx*mx+mx) / 2)
      parameter (maxref=256)
c
      parameter (mxrootn = 50)
      integer mxv, mxv1, mxv2, mxv21,lspacev
      parameter (mxv=mxrootn*mxrootn)
      parameter (mxv1=(mxv*mxv+mxv) / 2)
c nomax  : max. # of mo's
      parameter (nomax  =     256)
c
      integer nav,lenwrd,ndeke
      integer noeint,mvect,monec,mx2,mx21,lspace,needv
      integer need2,iematc,iotc,ijonc
c
      nav = lenwrd()
      ndeke = num*(num+1)/2+1
c noeint : max # of 1-electron integrals
      noeint =  ndeke
c
      mvect = mxrootn * maxref
      monec = nidmax + 3 *ndeke + 3 * noeint
      nedim0 = nedim
      mx2 = mx * mx
      mx21 = max(mx2,mx1)
      lspace = mx1 + mx2 + mx21
c
      mxv2 = mxv * mxv
      mxv21 = max(mxv2,mxv1)
      lspacev = mxv1 + mxv2 + mxv21
c
      lspace = max(lspace,lspacev)
c
c     allow for aftci
      needv = max(lspace,mvect)
c
c      ensure nteint is a multiple of transformed
c      integral file record length
       nteint0 = (nteint-1)/ nidmax + 1
       nteint0 = nteint0 * nidmax
       mdi0 = mdi * (ilast+1) + 1
c      mdi0 = mdi * ilast + 1
c      mdi0v =  mdi * (ilast+1) + 1
       mdi0v=  mdi +1
c      ensure iotm is multiple of nav
       iotm0 = (iotm-1)/nav + 1
       iotm0 = iotm0 * nav
c
c      guide reduction strategy
       integc = (nteint0-1)/ nav + 1
       mdic =  2 * mdi0 + 4 * mdi0v
       iotc = iotm0 + 2 * (iotm0 /nav)
       ijonc =  mx * mx / nav  + mx * (nomax+nomax) / nav
       iematc =  nedim0
       need = integc + monec + mdic + iematc + mvect +
     +        iotc + ijonc
       need2 = need + needv
       if (oprint) then
        write(6,*) ' integc = ', integc
        write(6,*) ' mvect =  ', mvect
        write(6,*) ' monec =  ', monec
        write(6,*) ' mdic =   ', mdic
        write(6,*) ' iematc = ', iematc
        write(6,*) ' iotc =   ', iotc
        write(6,*) ' ijonc =  ', ijonc
        write(6,*) ' need =   ', need
        write(6,*) ' need2 =  ', need2
       endif
c
       return
       end
      subroutine adler_alloc(ibase,num,
     +     nteint0,iotm0,mdi0,mdi0v,nedim0,
     +     iston,ipey,icoul,iexc,ipey0,icoul0,iexc0,
     +     ici1,ici2,ici3,ici4,ici5,ici6,iemat,ivect,
     +     iottr,iotnew,iot0,maindf,jconb,
     +     idiag,idiagv,idiag2,idiag3,idiag3v,need,oprint)
c
      implicit none
c
      integer ibase,num,nteint0,iotm0,mdi0,mdi0v,nedim0
      integer iston,ipey,icoul,iexc,ipey0,icoul0,iexc0
      integer ici1,ici2,ici3,ici4,ici5,ici6,iemat,ivect
      integer iottr,iotnew,iot0,maindf,jconb
      integer idiag,idiagv,idiag2,idiag3,idiag3v,need,last
      logical oprint
c
      integer nidmax,nomax,mx,mx1,maxref,mxrootn
c nidmax : dimension of the integral records of stoney
      parameter (nidmax =    2000)
c
c nomax  : max. # of mo's
      parameter (nomax  =     256)
      integer mx2,mx21
      parameter (mx=256)
      parameter (mx2=mx*mx)
      parameter (mx1=(mx*mx+mx) / 2)
c
      parameter (maxref=256)
c
      parameter (mxrootn = 50)
c
      integer mxv, mxv1, mxv2, mxv21
      parameter (mxv=mxrootn*mxrootn)
      parameter (mxv1=(mxv*mxv+mxv) / 2)
c
      integer nav,lenwrd,ndeke,noeint,mvect,nomax2

      nav = lenwrd()
      ndeke = num*(num+1)/2+1
c noeint : max # of 1-electron integrals
      noeint =  ndeke
      mvect = mxrootn * maxref
      mx21 = max(mx2,mx1)
c
      mxv2 = mxv * mxv
      mxv21 = max(mxv2,mxv1)
c
      nomax2 = nomax + nomax
c
c ston
c     use real*4 for twoe regardless 
c     iston = ibase + (nteint0-1)/nav + 1
c rwah real*8 intgr pointers
      iston = ibase + nteint0
c  pey
      ipey =  iston + nidmax
c  acoul
      icoul = ipey +  ndeke
c  aexc
      iexc  = icoul + ndeke
c  pey0
      ipey0 = iexc  + ndeke
c  acoul0
      icoul0 = ipey0 + noeint
c  aexc0
      iexc0 = icoul0 + noeint
c
c  now for six mdi arrays
c
      ici1 = iexc0 + noeint
      ici2 = ici1  + mdi0
      ici3 = ici2  + mdi0
      ici4 = ici3 + mdi0v
      ici5 = ici4 + mdi0v
      ici6 = ici5 + mdi0v
      iemat = ici6 + mdi0v
c
      ivect = iemat + nedim0
c
      iottr = ivect + mvect
*****
c     iotnew = iottr  + iotm0/ nav
c     iot0 =   iotnew + iotm0/ nav
c     maindf = iot0 + iotm0 / nav
*****
      iotnew = iottr  + iotm0 
      iot0 =   iotnew + iotm0 /nav
      maindf = iot0 + iotm0  / nav
      jconb =  maindf + mx * mx / nav
c    
      idiag =   jconb  + (mx * nomax2) / nav
      idiagv = idiag + mxv1
      idiag2 = idiagv + mxv2
      idiag3 = idiag2 + mxv21
      idiag3v= idiag3 + mx1
      last =   idiag3v + mx2
      need = last - ibase
c
      if(oprint) then
       write(6,*) ' adler_alloc'
       write(6,*) ' ibase =  ', ibase
       write(6,*) ' iston =  ', iston
       write(6,*) ' ipey =   ', ipey
       write(6,*) ' icoul =  ', icoul
       write(6,*) ' iexc =   ', iexc
       write(6,*) ' ipey0 =  ', ipey0
       write(6,*) ' icoul0 = ', icoul0
       write(6,*) ' iexc0 =  ', iexc0
       write(6,*) ' ici1 =   ', ici1
       write(6,*) ' ici2 =   ', ici2
       write(6,*) ' ici3 =   ', ici3
       write(6,*) ' ici4 =   ', ici4
       write(6,*) ' ici5 =   ', ici5
       write(6,*) ' ici6 =   ', ici6
       write(6,*) ' iemat =  ', iemat
       write(6,*) ' ivect =  ', ivect
       write(6,*) ' iottr =  ', iottr
       write(6,*) ' iotnew = ', iotnew
       write(6,*) ' iot0   = ', iot0  
       write(6,*) ' maindf = ', maindf 
       write(6,*) ' jconb  = ', jconb  
       write(6,*) ' idiag  = ', idiag  
       write(6,*) ' idiagv = ', idiagv
       write(6,*) ' idiag2 = ', idiag2  
       write(6,*) ' idiag3 = ', idiag3
       write(6,*) ' idiag3v = ', idiag3v
       write(6,*) ' need =   ', need
      endif
c
       return
       end
      subroutine adler_reduce(ofirst,integc,mdic,
     + need,ntotal,iword,
     + nteint,nteint0,o2e,
     + ilast,mdi,mdi0,mdi0v,ovect,oprint)
c
      implicit none
c
      real*8 fact
      integer integc,mdic,need,ntotal,iword,nteint0
      integer nteint,ilast,mdi,mdi0,mdi0v
      logical ofirst,o2e,ovect,oprint
c
      integer iwords,integmdi,jlast,i2e
c
c     must reduce memory requirement by iword words
c     assume that a two-fold reduction in integrals
c     and four-fold reduction in vectors is acceptable
c
      o2e = .false.
      ovect = .false.
      nteint0 = nteint
      mdi0 = mdi * (ilast+1) + 1
      mdi0v = mdi+1
      jlast = ilast
c
      iword = need-ntotal
      iwords = iword
      if (oprint) then
       write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
       write(6,*)'ntotal,need = ', ntotal,need
       write(6,*)'integc,mdic = ', integc,mdic
       write(6,*)'looking to reduce memory requirement by', iword
      endif
c
c     consider 2-e integrals or vectors
      integmdi = integc - mdic
      if (integmdi.gt.0) then
       call adler_reduce_i(integc,nteint0,iword,i2e,o2e,
     +                     oprint)
       if (o2e) then
        fact = dble(nteint0) / dble(nteint)
        if (fact.lt.0.5d0) then
c      excessive reduction
         if(oprint)write(6,*)'excessive reduction in integrals'
         nteint0 = nteint
         o2e = .false.
         iword = iwords
        endif
       endif
      endif
c
      if(.not.o2e) then
       call adler_reduce_v(ilast,jlast,mdi,mdi0,mdi0v,
     +                     iword,ovect,oprint)
       if (ovect) then
        fact = dble(jlast) / dble(ilast)
        if (fact.lt.0.25d0) then
c       excessive reduction
         if(oprint)write(6,*)'excessive reduction in vectors'
         ovect = .false.
         iword = iwords
         mdi0 = mdi * (ilast+1) + 1
         mdi0v = mdi+1
        endif
       endif
      endif
c
      if(oprint) then
       write(6,*)'suggested memory reductions'
       write(6,*)'nteint: from ', nteint, ' to ', nteint0
       write(6,*)'vectors: from ', ilast, ' to ', jlast
       write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
      endif
c
      if(.not.ovect .and. .not. o2e) then
       if (ofirst) then
c      in trouble .. just try reducing both
         write(6,*)'memory reduction mandated'
           ovect = .true.
           o2e = .true.
           nteint0 =  nteint0 / 2
           mdi0 = mdi0 / 2
       endif
      endif
c
      return
      end
      subroutine adler_reduce_i(integc,nteint0,iword,i2e,o2e,
     +                          oprint)
c
      implicit none
      integer integc,nteint0,iword,i2e
      logical o2e,oprint
c
      integer nav, lenwrd
c
      nav = lenwrd()
      if (integc.gt.iword) then
       i2e  = iword * nav
       if(oprint)write(6,*) 'reduce no. of integrals by ', i2e
       o2e = .true.
       nteint0 = nteint0 - i2e
       iword = 0
      endif
c
      return
      end
      subroutine adler_reduce_v(ilast,nvect_mdi,
     +                 mdi,mdi0,mdi0v,iword,ovect,oprint)
c
      implicit none
      integer ilast,mdi,mdi0,iword,nvect_mdi,mdi0v
      logical ovect,oprint
c
      integer mdi_mem_all,mdi_mem_all_0
      integer mdi0_mdi,mdi_mem
c
      nvect_mdi = ilast
      ovect = .false.
      mdi0v = mdi+1
      mdi0 = mdi * (nvect_mdi+1) + 1
      if (nvect_mdi.gt.2) then
c
       mdi0_mdi = mdi0
       mdi_mem_all = 2 * mdi0_mdi + 4 * mdi0v
c***   mdi_mem_all = 6 * mdi0_mdi 
       mdi_mem_all_0 = mdi_mem_all
c      reduce need i.e. mdi_mem_all by iword words
10     mdi_mem = mdi_mem_all_0 - mdi_mem_all
       if(mdi_mem.lt.iword.and.nvect_mdi.gt.1) then
        nvect_mdi = nvect_mdi - 1
        mdi0_mdi = mdi * (nvect_mdi+1) + 1
        mdi_mem_all = 2 * mdi0_mdi + 4 * (mdi + 1)
c***    mdi_mem_all = 6 * mdi0_mdi 
        if(oprint)
     +  write(6,*) 'iword, mdi_mem_all, nvect_mdi',
     +              iword, mdi_mem_all, nvect_mdi
        go to 10
       endif
c
       if(oprint)write(6,*)'nvect_mdi, mdi_mem_all',
     +           nvect_mdi, mdi_mem_all
c
       if(oprint) write(6,*)
       mdi0 = mdi0_mdi
       ovect = .true.
       mdi_mem_all = 2 * mdi0_mdi + 4 * mdi0v
c***   mdi_mem_all = 6 * mdi0_mdi 
c      mdi0v = mdi0_mdi
       iword = iword - (mdi_mem_all_0 - mdi_mem_all + 1)
c
      endif
c
      return
      end
      subroutine adler_dat
c
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
c
      do loop =1,10
       jerk(loop) = 0
       jbnk(loop) = 0
       jbun(loop) = 0
      enddo
c
      jerk (2)  = 30
      jerk (3)  = 118
      jerk (4)  = 158
      jerk (5)  = 248
      jerk (6)  = 275
      jerk (7)  = 322
c
      jbnk (1)  = 4
      jbnk (2)  = 17
      jbnk (3)  = 7
      jbnk (4)  = 15
      jbnk (5)  = 4
      jbnk (6)  = 5
      jbnk (7)  = 7
c
      return
      end
      subroutine prepref(nteint,mdi,iotmf,vect,mxmr)

      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer nteint, mdi, iotmf
c
      real*8 vect
      integer mxmr
      dimension vect(mxmr)
c
***** problem with extrapolation
      integer iselec,izus,iselect,iselecz,nrootci
      common /cselec/ iselec(mxroot),izus(mxroot),
     +                iselect(mxroot),nrootci,
     +                iselecz(mxroot)
      common/rb/ ifrk,iswh
c
      real*8 esav0,esav1,csum,de0,de1,ueb0,ueb1,egeys
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .               de0(mxroot)  ,de1(mxroot),
     .               ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .               egeys(mxroot)
c 
      real*8 b,trsum,ew
      integer kj,mjx,lj,ij
      integer nj,ntil,nbal,isym,jsym,lsym,ncomp
      integer ibal,itil,mcomp
      integer istm,isp
      integer kmap,khog,jkan
      common /miscop/ trsum(10,mxroot),ew(mxroot),b(500),
     + lj(nirmax),ij(8),kj(nirmax),mjx(nirmax),
     + nj(nirmax),ntil(nirmax),nbal(9),
     + isym(nirmax),jsym(jabmax),lsym(800),ncomp(100),
     + ibal(nirmax),itil(nirmax),mcomp(100),istm(10),
     + jkan(1000),kmap(504),khog(48)
c
      integer icount,irech
      common /rhsc/ icount,irech
c --- issk : superkategorien auf disk
c --- issk : supercategories on disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
      integer ncons,idrem1
      common /rhus/ idrem1,ncons(iswhm)
cvp
      real*8 vnuc,zero
      integer nsel,iorbs,knu
      real*8 ab, eb
      common/junk/ ab(5292),eb(2304)
c
cvp
cvp weitere groessen, die von ft33 gelesen werden
cvp more sizes, read from ft33
      integer ndt,kml
c
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu ,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
      real*8 sump3
      integer iadl,irec,iadp3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
cvp
cvp for integralvorsortierung nach foxy
cvp for integral presort to foxy
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
cvp
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
cvp vom programm benoetigte files
cvp files required by the program
c
c  file unit nos. now held in ftap5 and initialised in adler
c
c
      call rewftn(mtype)
      call rewftn(nston)
      call rewftn(ltape)
      call rewftn(ideli)
      call rewftn(mtapev)

      iswh=iswhm
      ifrk=500
      ifr1=1-ifrk

      ideks(1)=0
      do 1 i=1,ideksm-1
   1  ideks(i+1)=ideks(i)+i
      read(nston) n,iwod,nid,ksum,imo,kj,mjx,lj,nj,nsel,ntil,nbal
     & ,isym,jsym,iorbs,knu,lsym,ncomp,vnuc,zero
      read(nston) nit,lg,ij,ibal,itil,mcomp
      nbal(9)=ncomp(1)

      nteint = lg

      do 231 i=1,nirmax
        ijone(i)=ij(i)
  231 continue
cvp nrec31 = anzahl der records von zweielektronenintegralen auf ft31
cvp nrec31 = number of two-electron integral records on ft31
      nrec31=(lg-1)/nid+1

      do 3050 i=1,nirmax
 3050 ncimo(i)=lj(i)
      ix=0
      do 5 i=1,n
      kap=lj(i)
      if (kap.eq.0) go to 5
      do 6 j=1,kap
      ix=ix+1
      nirred(ix)=i
      mopos(ix)=j
   6  continue
   5  continue
cvp
      read (ltape) jsec,nrootx,nytl,nplu,ndub,vect,ew,mconf
      read (ideli) (b(i),i=1,5),isp,vect,ew
cvp
      ny=0
      do 7 i=1,iswh
      niot(i)=ny
      nod(i)=nytl(i)-ndub(i)
      if (mconf(i).eq.0) go to 7
      ig=1
      nx=nytl(i)
      jg=0
      read (ltape)
      read (ltape)
      read (ltape) jkan
   8  if (jkan(ig).eq.0) go to 9
      jg=jg+1
      do 10 j=1,nx
      ny=ny+1
      if (ig.lt.iwod) go to 10
      read (ltape) jkan
      ig=0
   10 ig=ig+1
      go to 8
    9 nconf(i)=jg
    7 continue
cdebug
*     do 322 ii1=1,iswhm
*        write(iwr,*) 'nconf(',ii1,')=',nconf(ii1)
*322  continue
       iotmf = ny
      
       write(iwr,51233) iotmf
51233  format(1x,'memory for configurations'/
     +        2x,'required (iotm) = ', i12)
       call rewftn(ltape)
cvp
cvp oeffnen der table
cvp open table
c     call setbfc
cvp lesen des kopfes der table sowie der c- und der b-matrix
cvp read the table header and the c- and the b-matrix
c     call tbread(ifrk,ifr1,ntab)
*      write(iwr,*) 'from rumpxs after opening the table'
cvp
cvp aufbau von feldern, die konfigurationsvergleich erleichtern
cvp sollen
cvp build fields that should facilitate comparing configurations
c     call chfeld(iottr,iot0,iotm)
*      write(iwr,*) 'from rumpxs after chfeld'
cvp
cvp schreiben des ersten records von ft35
cvp write the first records of ft35
cvp
      m=nplu(1)+2*ndub(1)

      read (ltape)
cvp
      do 22 i=1,iswh
      md=mconf(i)
cdebugwrite(iwr,*)'bei 22 : md=',md
      if (md.eq.0) go to 22
      nc=nconf(i)
cdebugwrite(iwr,*)'bei 22 : nc=',nc
      if (nc.ne.0) go to 836
      read (ltape)
      read (ltape)
      read (ltape)
      go to 22
  836 read (ltape) ndt,kml,ab
cdebugwrite(iwr,*)'ndt,kml,copied'
      read(ltape) khog,kmap,eb
cdebugwrite(iwr,*)'khog,kmap,copied'
      nx=nytl(i)*nc
      nx=nx/iwod+1
      do 900 j=1,nx
  900 read(ltape)
   22 continue
      call rewftn(ltape)
      jsum=0
*      write(iwr,*) 'from rumpxs before 901 '
      do 901 i=1,iswh
        inopen=nytl(i)-ndub(i)
        nbeta=inopen-nplu(i)
        if (nbeta.eq.0) then
           nsafsk=1
         else
           nsafsk=ibinom(inopen,nbeta)-ibinom(inopen,nbeta-1)
        endif
        ncons(i) = nsafsk
  901   jsum=jsum+nconf(i)*nsafsk
      jsum=(jsum-1)/ifrk+1
      do 902 i=1,jsum
      read (ideli) b
  902 continue
      read(ideli)nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm

cbe abspeichern zu ende
cbe finished storing
      read(ltape)
      do 23 i=1,iswh
      md=mconf(i)
      if (md.eq.0) go to 23
      read(ltape)
      read(ltape)
      nc=nconf(i)
      if (nc.ne.0) go to 837
      read (ltape)
      go to 23
  837 nx=nytl(i)*nc
      nx=nx/iwod+1
      do 24 j=1,nx
      read(ltape) jkan
cdebugwrite(iwr,*) 'jkan written'
   24 continue
   23 continue
c
      mdi = istm(4)
c
c --- lese die e/v
c --- read the e/v
*      write(iwr,*) 'from rumpxs before rderz '
c     call rderz(iotnew,iotm,maindf,maxr,jconb,nomax2,
c    +           mtapev,iwr,ofirst,iotmf)
*      write(iwr,*) 'from rumpxs after rderz '
      call rewftn(ltape)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
c
      subroutine aftci(civec,mdi,vect,nvect,iottr,iot0,iotm,
     +                 nbak,energy,odebug,oprinsym)
c
c  this subroutine performs the davidson-correction and writes ft36.
c
c
c   (1993) h.u.s.
c  correction so that the right mains get printed:
c  mipos at the configuration level is no longer needed, in stead one
c  accesses directly where needed.
c  b.e.  26.7.1994
c
c  correction, so that in the old way a correct fi36 is written
c  b.e.  28.7.1994
c
c  nbak = 0    first pass         does everything
c              overlap on ueb0
c  nbak = 1    second pass        computes only the overlap and writes
c              it on ueb1 in stead of on ueb0
c
c  correction, so that a sensible extrapolation can be guaranteed:
c
c  a.   builds in aftci.f the overlap with all starting vectors
c
c  b.   performs extrapolation while applying the corrections for 
c       which the largest overlap is found
c
c  c.   incorporate the weighted extrapolation
c
c  at the end of outp.f another record is written onto ft36 
c  therefore one shouldn't rewind this file, silly enough the file 
c  file36 in outp.f is not called iput but itape
c
c  begin 28.8.1994 b.e.
c
c#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      logical odebug,oprinsym
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 edavit,cdavit,extrapit,eigvalr,cradd, 
     +     weighb, rootdel, ethreshit
      integer mbuenk,nbuenk,mxroots
      logical ifbuen, ordel, odave, oweight
      logical odavit
      common /comrjb2/edavit(mxroot),cdavit(mxroot),
     +                extrapit(mxroot),ethreshit(mxroot),
     +                odavit(mxroot),
     +                eigvalr(maxref),cradd(3),weighb,
     +                rootdel,ifbuen,mbuenk,mxroots,
     +                ordel,odave,nbuenk,oweight
c
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      integer nstone, mstvt, nf01,nf11
      common /ftap/ nstone(10), mstvt, nf01(5), nf11
c
      integer ical,nmul
      common /tap/ ical(4),nmul
c
c
      real*8 civec
      integer mdi
      dimension civec(mdi)
c
      real*8 vect
      dimension vect(nvect)
c
      integer iottr,iot0,iotm
      dimension iottr(iotm),iot0(iotm)
c
      integer mx,mipos,idim
      integer ii,iorbs,nsaf,nslat
      dimension nsaf(iswhm),nslat(iswhm)
      parameter (mx=256)
      real*8 rzero
c --- parameter for ft36
      real*8 d,y
      common/junk/d(mx*mx),y(mx*(mx+1)/2)
c
      real*8 an
c --- parameters for ft36
      real*8 hh,gg,ff
      integer ihog,imap
      integer jkan0,mcomp
c vn contains the csf --> spinproj transformation
      real*8 vn
      common /rsghm4/vn(2304),
c --- parameters for ft36
     +         hh(100),gg(500),ff(1000),
     +         jkan0(500),ihog(48),imap(504),
     +         b(232),mcomp(100)
      integer knu
      dimension ew(mxroot)
*  ew+cf : yn
c     dimension yn(8231)
      integer kj,lj
      dimension lj(8),kj(8)
      integer ibal,itil
      dimension ibal(8),itil(8)

cbe dimension for reading the start vectors from mstvt

      integer nroop

cbe end
c
      common /can/ an(5292)
      common /rwork/ niot0(iswhm)
c
      parameter (rzero=0.0d0)
      common /c/ idim,nrotr
c --- common for extrapolation
      real*8 esav0,esav1,csum,de0,de1,ueb0,ueb1,egeys
      integer ilifq
      character*1 xcoff
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .              de0(mxroot),de1(mxroot),
     .              ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .              egeys(mxroot),ilifq(mxroot),
     .              xcoff(mxroot,mxroot)
c ------------------------
      integer nplu
      dimension nplu(iswhm)
c
      integer nsac,idrem1,ndet
      common /rhus/ idrem1,nsac(iswhm),ndet(iswhm)
c
      common /cmpos/ mipos(mx),mipos1(mx)
      integer mcpos
      dimension mcpos(mx)
      integer nroot,ntch,kprin,ndeci
      integer ifirst0,ilast0
      integer iselec,nsaft
      common /cinput1/ nroot, ifirst0,ilast0, ntch, kprin,
     .                ndeci,
     .                icode,konf ,keps,iggey,
     .                istart,nrootz

      common /cselec/ iselec(mx)
c --- old common from conmdi (partially)
      real*8 esav36
      integer ixt,i0
      common /per/ ixt,i0,esav36(4,mxroot)
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
c
      integer mxretain
      parameter (mxretain = 50)
      real*8 energret,facret
      integer idumm,nopret,ipret,jkonret,nretain,iretain
      common /cski1/ idumm(maxshl),
     + energret(mxretain),facret,nopret(mxretain),
     + ipret(mxretain),jkonret(mxretain*maxshl),nretain,
     + iretain
      real*8 egey1,trash,tdel
      common /parkin/egey1,trash,tdel
c
      integer lun
      logical oi3
      character*4 zsymx(8)
      integer nirrep(8)
      character*100 remstring
c
      data zsymx/'s1','s2','s3','s4','s5','s6','s7',
     1           's8'/
c
c NB
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
c NB
c --- write header
      write(iwr,*)
      write(iwr,1400)
1400  format(1x,'><><><><><><><><><><><><><>><><><><><><><>')
      write(iwr,1401) 
1401  format(1x,'>>>>> commence CI vector analysis : <<<<<<')
      if (nbak.eq.0) then
       write(iwr,1500) trash
1500   format(1x,'>>>>>> T=',f5.1,
     +          ' uh threshold problem <<<<<<')
      else
       write(iwr,1501) (trash+tdel)
1501   format(1x,'>>>>>> T=',f5.1,
     +          ' uh threshold problem <<<<<<')
      endif
      write(iwr,1400)
      write(iwr,*)
c --- prepare the configuration statistic
      nsaft = 0
cbe      nrootl= iselec(ilast0)
      if (odebug) then
      write(iwr,*) 'aftci: nroot,ifirst0,ilast0', nroot,ifirst0,ilast0
c     write(iwr,*) 'iselec(1-5)',(iselec(idoof),idoof=1,5)
      write(iwr,*) 'aftci: nrootz = ', nrootz
      endif
      nrootl= iselec(nroot)
      do i=1,iswhm
        nsaf(i)   = nsac(i) * nconf(i)
        nslat(i)  = ndet(i) * nconf(i)
        nsaft     = nsaft + nsaf(i)
        niot0(i)  = niot(i)
      end do
      write(iwr,*)
     . 'super category                1        2       3       4      5'
      write(iwr,*)
     .'----------------------------------------------------------------'
      write(iwr,1111) 'configurations      : ',nconf
      write(iwr,1111) 'safs                : ',nsaf
      write(iwr,1111) 'slaterdeterminants  : ',nslat
      write(iwr,*)
c
      write(iwr,841) nsaft
c --- preparation of file management
      ixt   =4
      lun = lun20
c
c compute overlap with start vectors that are taken from mstvt 
c these will be read already in prep0 so that mstvt doesn't have to be
c opened again.
c --- extract from unit=mstvt
c        write(iwr,*) 'start vectoren read'
cbe opening not needed, because already open
c        open (unit=mstvt,file='gstarv.dat',status='old',
c     .      form='unformatted')
        call rewftn(mstvt)
        read(mstvt) nroop
c       write(6,*)'nroop = ', nroop
        read(mstvt) d
        read(mstvt) y

* compute overlap using mipos (for the big problem)
* =nbak=0 and mipos1 for the small problem = nbak=1

      if(odebug) then
       if (nbak.eq.0) then
        write(iwr,*) 'aftci:big problem: overlap'
       endif
       if (nbak.eq.1) then
        write(iwr,*) 'aftci:small problem: overlap'
       endif
      endif
c
       call rewftn(lun)
      do ilauf = 1,nroot
        iwurz = iselec(ilauf)
        write(iwr,7778) iwurz
        read (lun) (civec(igi),igi=1,nsaft)
* biggest coefficients
        do idudel = 1,nsaft
           if(abs(civec(idudel)).gt.cptol) then
             write(iwr,50) civec(idudel),idudel
           endif
        enddo

        ueb=0.0d0
        iueb=0
c       do jlauf = 1,ilast0
        do jlauf = 1,nrootz
          if (odebug)
     +    write(iwr,*) 'aftci: ilast0,nrootz=',ilast0,nrootz
          iroot = jlauf
* calculation of the starting position of the starting root in d
* taken from prep0.f
          ipos=nroop*nroop - iroot*nroop
          ueber=0
c         iueb=0
          if (nbak.eq.0) then
            do izz = 1,nroop
               ipos=ipos+1
               ueber=ueber+civec(mipos(izz))*d(ipos)
            enddo
          endif
          if (nbak.eq.1) then
            do izz = 1,nroop
               ipos=ipos+1
               ueber=ueber+civec(mipos1(izz))*d(ipos)
            enddo
          endif
          if (odebug) then
c writing out the overlap
           write(iwr,*) 'aftci:overlap with ',jlauf,
     +                  '-th start vector'
           write(iwr,*) ueber
          endif
c sort into  ueb0 or ueb1
          if (nbak.eq.0) then
             ueb0(jlauf,ilauf) = ueber
          else
             ueb1(jlauf,ilauf) = ueber
          endif
          if (dabs(ueber).gt.dabs(ueb)) then
             ueb = ueber
             iueb= jlauf
          endif

        enddo
        write(iwr,851) iueb, ueb

* end of do loop over 'do ilauf = ifirst0,nrootl'
       enddo
c
        write(iwr,887) 
c       do jlauf = 1,nrootl
c       do jlauf = 1,nrootz
          if (nbak.eq.0) then
c            write(iwr,850) jlauf, ueb0(jlauf,ilauf)
             call writer(ueb0,ilifq,xcoff,nrootz,nroot,mxroot)
          else
             call writer(ueb1,ilifq,xcoff,nrootz,nroot,mxroot)
c            write(iwr,850) jlauf, ueb1(jlauf,ilauf)
          endif
c       enddo


      if (nbak.gt.0) then
        write(iwr,889) 
        return
      endif

cbe ended for the time being

      call rewftn(lun)
      call rewftn(nhuk)
      if(ifbuen) call rewftn(nf11)
c
c --- copy iot onto iot0
c --- save iot,niot
c
      do 220 i=1,iswhm
        istar=niot(i)
        do 210 j=1,nytl(i)
          do 200 k=1,nconf(i)
            iot0(istar+j+(k-1)*nytl(i)) =
     &      iottr(istar+k+(j-1)*nconf(i))
  200     continue
  210   continue
  220 continue
cvp
      do 9123 i=1,iswhm
        niot0(i) = niot(i)
9123  continue
c
c
c --- do the davidson
c
c     write(iwr,*) (mipos(iraus),iraus=1,256)
      do ii=1,nroot
        c20= rzero
        read (lun) (civec(igi),igi=1,nsaft)
        do i=1,mx
         mm = mipos(i)
          if (mm.ne.0) then
c           write(iwr,*) 'mm=',mm,'civec =', civec(mm)
            c20 = c20 + civec(mm)*civec(mm)
          end if
        end do
        csum(ii) = c20
        if(odebug) then
         write(6,*)' aftci: root, csum = ', ii,csum(ii)
        endif
      end do
c
c --- make mipos at the configuration level : mcpos
c --- oldpos never set ?
      oldpos = 0
c
      imc = 1
      do ii=1,mx
       mm = mipos(ii)
       if (mm.ne.0) then
        nss = 0
        nsse= nsaf(1)
        ilocsk = 5
        do ii1=1,(iswhm-1)
         if (mm.gt.nss) then
          if (mm.le.nsse) then
           ilocsk = ii1
          end if
         end if
         nss = nss+nsaf(ii1)
         nsse = nsse+nsaf(ii1+1)
        end do
c
        newpos = 0
        newsaf = 0
        do ii1=1,(ilocsk-1)
         newpos = newpos+nconf(ii1)
         newsaf = newsaf+nsaf(ii1)
        end do
        newpos = newpos + (mm-newsaf)/nsac(ilocsk)
        if (oldpos.ne.newpos) then
         mcpos(imc)=newpos
         imc = imc+1
        end if
        oldpos = newpos
       end if
      end do
cbe
cbe      write(iwr,*) (mcpos(iraus),iraus=1,20)
cbe
c
c  --- write ft36
c
       read (nhuk) iwod,vnuc,zero,imo,m,
     .  nconf,nytl,nplu,ndub,iswh,ksum,
     1  iorbs,jsec,nrotr, vect,ew,
     .  ibal,itil,mcomp,kj,lj,nn,
     . ifrk,knu,b

c first record
      write(iput)
     . iwod,vnuc,zero, imo, m , nconf,
c      format, nuclear-nuclear-potential,scf-energy,number of mo, 
c              number of elektrons, number of configurations per super
c              category
     . nytl,
c      number of ?
     . nplu,
c      number of alpha-spins per sk
     . ndub,
c      number of sk (=5)
     . iswh,
c      number of sk (=5)
     . ksum,
c      number of scf mo's
     . iorbs,
c      number of ao's
     . knu,
c      number of nucleii
     . ibal,itil,
cbe it seems that mcomp and kj were missing
     . mcomp,kj,
cold . yn,
c      information coefficients, nucleii information
     . lj,nn,ifrk,b
c
cbe debug print to test where the bug is
cbe       write(iwr,*)'nconf=',nconf
cbe       write(iwr,*)'nytl=',nytl
cbe       write(iwr,*)'nplu=',nplu
cbe       write(iwr,*)'ndub=',ndub
cbe       write(iwr,*)'iswh=',iswh
cbe       write(iwr,*)'ksum=',ksum
cbe       write(iwr,*)'iorbs=',iorbs
cbe       write(iwr,*)'knu=',knu
cbe       write(iwr,*)'ibal=',ibal
cbe       write(iwr,*)'itil=',itil
cbe       write(iwr,*)'mcomp(1)=',mcomp(1)
cbe       write(iwr,*)'kj=',kj
cbe
cbe       write(iwr,*)'lj=',lj
cbe       write(iwr,*)'nn=',nn
cbe       write(iwr,*)'ifrk=',ifrk
cbe
cbe       write(iwr,*)'z(1,2)=',z(1),z(2)
cbe debug print end
       call rewftn(lun)
cdebug       write(iwr,*) 'nrootl=',nrootl
       intrec = 1
c
      do i=1,nroot
          iwurz = iselec(i)
c
c       note that eneg in the write statements below (to iput)
c       was not set .. disaster in tm code ...
c       set here to esav0
c
          eneg = esav0(i)
          energy = eneg
c
          write(iwr,7777) iwurz
          read (lun) (civec(igi),igi=1,nsaft)
          niconf=0
          nisaf =0
          nm    =0
          ipo2=1
c --- position in iot-field
          iotpos=1
c --- loop over the super categories
          call rewftn(nhuk)
          read   (nhuk)
          iotpos = 1
          nod1 = nmul - 3
          do j=1,iswhm
             nc = nconf(j)
             nod1 = nod1 + 2
             if (odebug) write(iwr,*) 
     +       'aftci: j, nc,intrec,nod1 = ', 
     +       j, nc, intrec, nod1
             if (nc.ne.0) then
                read (nhuk) ndt,kml,an
                read (nhuk) ihog,imap,vn
                nl = nytl(j)
                imax = 100
                l1 = ifrk/nl
                l2 = iwod/ndt
                l3 = ifrk/kml
                if (l1.lt.imax) then
                   imax = l1
                end if
                if (l2.lt.imax) then
                   imax = l2
                end if
                if (l3.lt.imax) then
                   imax = l3
                end if
                nhb = (nc-1)/imax+1
                write(iput)
     .          nhb,imax,ndt,kml,imap,ihog,eneg
                intrec = intrec + 1
                if (odebug)write(iwr,*) 
     +           'aftci: nhb,imax = ', nhb, imax
c third record
                nx=0
                nxl = imax*nl
                igl = imax*kml
                ifl = imax*ndt
c          bilde c**2
                ns = nsac(j)
                nc0 = nc/imax
                ncr = nc - imax*nc0
                do ii=1,nc0
c --- nullify
                   do inul=1,500
                      jkan0(inul)    = 0
                      gg(inul)       = rzero
                      ff(inul)       = rzero
                      ff(1001-inul)  = rzero
                   end do
                   do inul=1,100
                      hh(inul)      = rzero
                   end do
                   ifi = 0
                   igi = 0
                   ihi = 0
c --- end nullify
                   do iii=1,imax
                      niconf = niconf + 1
                      do 37 ll=1,ndt
                         rorb = rzero
                         mm  = nm
                         lx  = ll - ndt
                         do 38 lll=1,kml
                            lx=lx+ndt
                            mm=mm+1
                            rorb=rorb+civec(mm)*an(lx)
   38                    continue
                         ifi = ifi+1
                         ff(ifi) = rorb
   37                 continue
                      do 39 ll=1,kml
                         rorb = rzero
                         mm  = nm
                         lx  = ll - kml
                         do 40 lll=1,kml
                            lx=lx+kml
                            mm=mm+1
                            rorb=rorb+civec(mm)*vn(lx)
   40                    continue
                         igi = igi+1
                         gg(igi)=rorb
   39                 continue
                      rorb = rzero
                      do 41 ll=1,kml
                         nm=nm+1
                         vm = civec(nm)
                         rorb = rorb+vm*vm
   41                 continue
                      ihi = ihi +1
                      hh(ihi) = rorb
                   end do
c --- filling jkan0
                   nloc = (imax*nl)
                   do ll=1,nloc
                      jkan0(ll) = (iot0(ipo2))
                      ipo2 = ipo2+1
                   end do
c   ---- end of block
                   write (iput) jkan0,ff,gg,hh
                   intrec = intrec + 1
                end do
           if (odebug)write(iwr,*) 'aftci: intrec = ', intrec
c --- nullify
cdebug     write(iwr,*) 'rest of ft36'
                do inul=1,500
                   jkan0(inul)    = 0
                   gg(inul)       = rzero
                   ff(inul)       = rzero
                   ff(1001-inul)  = rzero
                end do
                do inul=1,100
                   hh(inul)      = rzero
                end do
                ifi = 0
                igi = 0
                ihi = 0
c --- end nullify
c --- rest of data
                do ii=1,ncr
                   niconf = niconf + 1
                   do ll=1,ndt
                      rorb = rzero
                      mm  = nm
                      lx  = ll - ndt
                      do lll=1,kml
                         lx=lx+ndt
                         mm=mm+1
                         rorb=rorb+civec(mm)*an(lx)
                      enddo
                      ifi = ifi+1
c --- slater determinants
                      ff(ifi) = rorb
                   enddo
                   do ll=1,kml
                      rorb = rzero
                      mm=nm
                      lx  = ll-kml
                      do lll=1,kml
                         lx=lx+kml
                         mm=mm+1
                         rorb=rorb+civec(mm)*vn(lx)
                      enddo
                      igi = igi+1
                      gg(igi)=rorb
                   enddo
                   rorb = rzero
                   do ll=1,kml
                      nm=nm+1
                      vm = civec(nm)
                      rorb = rorb+vm*vm
                   enddo
                   ihi = ihi +1
                   hh(ihi) = rorb
                end do
c --- filling jkan0
                nloc = (ncr*nl)
                do ll=1,nloc
                   jkan0(ll) = (iot0(ipo2))
                   ipo2 = ipo2+1
                end do
c   ---- end of block
                if (ncr.ne.0) write (iput) jkan0,ff,gg,hh
                intrec = intrec + 1
c ---------------------------------
c
                do ii=1,nc
c                  c2 never initialised
                   c2 = rzero
                   do ij=1,ns
                      nisaf = nisaf + 1
                      c2 = civec(nisaf)*civec(nisaf) + c2
                   end do
c --- testing whether main
                   iotend = iotpos + nl - 1
                   mcloc = 0
                   do iji=1,mx
                     mctest=mipos(iji)
                     if(nisaf.eq.mctest) then
c                      if(c2.ge.cptolm) then
                        write(iwr,8888) 'm ',c2,
     .                  (iot0(igi),igi=iotpos,iotend)
c                      endif
                       if(ifbuen) then
                        write(nf11,101)nod1,c2,
     +                  (iot0(igi),igi=iotpos,iotend)
 101                    format (i4,',',f10.3,','/,256i4)
                       end if
                       mcloc = 2
c --- sum c**2
                     end if
                   end do
                   if (mcloc.lt.1) then
                     if (c2.gt.cptolcc) then
                         write(iwr,8888) '  ',c2,
     .                   (iot0(igi),igi=iotpos,iotend)
                     end if
                     if(ifbuen) then
                      if (c2.gt.cradd(1)) then
                       write(nf11,101)nod1,c2,
     +                 (iot0(igi),igi=iotpos,iotend)
                      endif
                     endif
                   end if
                   iotpos = iotpos+nl
                   c2 = rzero
                end do
             end if
          end do

          if (odebug) 
     +       write(iwr,*)'aftci: total no. of records = ', 
     +       intrec
          write(iwr,895) csum(i)
      enddo
c
c     now add the retained configurations from selection
c
      if(ifbuen) then
       if(nretain.gt.0) then
        igi = 1
        do loop=1,nretain
        write(nf11,101)nopret(loop),energret(loop),
     +  (jkonret(i),i=igi,igi+ipret(loop)-1)
        igi = igi + maxshl
        enddo
       endif
      endif
c
c rwah andere print
c
      if(oprinsym) then
        nirrep(1)=0
        oi3=.true.
        istep=3
        do iiz=2,8
           nirrep(iiz)=nirrep(iiz-1)+lj(iiz-1)
        enddo
        do iiz=1,8
           if (lj(iiz).gt.100) then
              oi3=.false.
              istep=4
           endif
        enddo
        write(iwr,8890)
        call rewftn(lun)
       do i=1,nroot
        iwurz=iselec(i)
        write(iwr,7780) iwurz
        read (lun) (civec(igi),igi=1,nsaft)
        niconf=0
        nisaf =0
        nm    =0
        ipo2=1
c --- position in iot-field
        iotpos=1
c --- loop over the super categories
        call rewftn(nhuk)
        read   (nhuk)
        iotpos = 1
        do j=1,iswhm
          nc = nconf(j)
          if (nc.ne.0) then
           read (nhuk) ndt,kml,an
           read (nhuk) ihog,imap,vn
           nl = nytl(j)
           ns = nsac(j)
           do ii=1,nc
            do ij=1,ns
             nisaf = nisaf + 1
             c2 = civec(nisaf)*civec(nisaf) + c2
            end do
c --- testing whether main
             iotend = iotpos + nl - 1
             mcloc = 0
             do iji=1,mx
               mctest=mipos(iji)
               if(nisaf.eq.mctest.and.c2.ge.cptolm) then
                 nopen=(iotend-iotpos+1)*2-m
                 krem=1
                 write(remstring(krem:krem+14),'(a7,f7.5)')'     m ',c2
                 krem=krem+14
                 write(remstring(krem:krem+4),'(i4)')nopen
                 krem=krem+4
                 nirr=1
                 nsym=0
                 ist=18
                 iopen=0
                 do iiz=iotpos,iotend
                     iorb=iot0(iiz) 
                     do iiy=1,8
                         if (iorb.lt.nirrep(iiy)) then
                            isym=iiy-1
                            goto 6100
                         endif
                     enddo
                     isym=8
 6100                continue
                     if (iopen.eq.nopen) then
                        write(remstring(krem:krem+3),'(a3)')'   '
                        krem=krem+3
                        ist=krem
                     endif
                     if (isym.ne.nsym) then
                         if (krem+4.ge.100) then
                            write(iwr,'(a)')remstring(1:krem)
                            do krem=1,100
                               remstring(krem:krem)=' '
                            enddo
                            krem=ist
                         endif
                         if (oi3) then
                            write(remstring(krem:krem+3),'(a1,a2)')' ',
     1                         zsymx(isym)
                            krem=krem+3
                         else 
                            write(remstring(krem:krem+4),'(a2,a2)')'  ',
     1                         zsymx(isym)
                            krem=krem+4
                         endif
                         nsym=isym
                     endif
                     iorb=iorb-nirrep(isym)
                     if (krem+istep.ge.100) then
                         write(iwr,'(a)')remstring(1:krem)
                         do krem=1,100
                             remstring(krem:krem)=' '
                         enddo
                         krem=ist
                     endif
                     if (oi3) then
                         write(remstring(krem:krem+istep),'(i3)')iorb
                     else 
                         write(remstring(krem:krem+istep),'(i4)')iorb
                     endif
                     iopen=iopen+1
                     krem=krem+istep
                 enddo
                 write(iwr,'(a)')remstring(1:krem)
                 mcloc = 2
c --- sum c**2
               end if
              end do
              if (c2.gt.cptolcc) then
               if (mcloc.lt.1) then
                 nopen=(iotend-iotpos+1)*2-m
                 do krem=1,100
                    remstring(krem:krem)=' '
                 enddo
                 krem=1
                 write(remstring(krem:krem+14),'(a7,f7.5)')'       ',c2
                 krem=krem+14
                 write(remstring(krem:krem+4),'(i4)')nopen
                 krem=krem+4
                 nirr=1
                 nsym=0
                 iopen=0
                 ist=18
                 do iiz=iotpos,iotend
                     iorb=iot0(iiz)
                     do iiy=1,8
                         if (iorb.lt.nirrep(iiy)) then
                            isym=iiy-1
                            goto 6101
                         endif
                     enddo
                     isym=8
 6101                continue
                     if (iopen.eq.nopen) then
                        write(remstring(krem:krem+3),'(a3)')'   '
                        krem=krem+3
                        ist=krem
                     endif
                     if (isym.ne.nsym) then
                         if (krem+4.ge.100) then
                            write(iwr,'(a)')remstring(1:krem)
                            do krem=1,100
                               remstring(krem:krem)=' '
                            enddo
                            krem=ist
                         endif
                         if (oi3) then
                            write(remstring(krem:krem+3),'(a1,a2)')' ',
     1                         zsymx(isym)
                            krem=krem+3
                         else 
                            write(remstring(krem:krem+4),'(a2,a2)')'  ',
     1                         zsymx(isym)
                            krem=krem+4
                         endif
                         nsym=isym
                     endif
                     iorb=iorb-nirrep(isym)
                     if (krem+istep.ge.100) then
                         write(iwr,'(a)')remstring(1:krem)
                         do krem=1,100
                            remstring(krem:krem)=' '
                         enddo
                         krem=ist
                     endif
                     if (oi3) then
                         write(remstring(krem:krem+istep),'(i3)')iorb
                     else 
                         write(remstring(krem:krem+istep),'(i4)')iorb
                     endif
                     iopen=iopen+1
                     krem=krem+istep
                 enddo
                 write(iwr,'(a)')remstring(1:krem)
               end if
              end if
             iotpos = iotpos+nl
             c2 = rzero
           end do
          end if
        end do
        write(iwr,8891)
c
c einde rwah andere print
c
       end do
c
      endif
c
      write(iwr,890) 
c
c --- formats
c
1111  format(a25,5i8)
8888  format(10x,a2,2x,f7.5,25i4,4(/21x,25i4))
890   format(/1x,'**** ci-vector file (ft36) written')
895   format(1x,'sum of main reference c*c =',f15.8/)
50    format(7x,f12.6,i10)
7777  format(/
     +  10x,'*************************'/
     +  10x,'* configuration weights *'/
     +  10x,'*************************'/
     *  10x,'***** Root ',i2,' ***********'/
     +  10x,29('=')/
     +  18x,'c*c',5x,'configuration'/
     +  10x,29('='))
7780  format(/
     *  5x,'***** Root ',i2,' ***********'/
     +  5x,29('=')/
     +  11x,'c*c',5x,'configuration'/
     +  5x,29('='))
889   format(1x,'**** no writing of ci-vector file (ft36)')
851   format(/1x,
     +   '+++ eigenstate has largest overlap with zero-order vector',i3/
     +1x,'+++ overlap = ', f20.10)
c850  format(1x,i9,f20.10)
887   format(/
     + 20x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'/
     + 20x,'+ Overlap of Eigenstates with the Zero-order Vectors +'/
     + 20x,'++++++++++++++++++++++++++++++++++++++++++++++++++++++')
7778  format(/
     +10x,'***************'/
     +10x,'* root = ', i4,' *'/
     +10x,'***************'//
     +1x,'ci vector read',5x,'ci coefficients'//
     +8x,29('=')/
     +8x,'coefficient',5x,'configuration'/
     +8x,29('='))
841   format(1x,'total number of safs = ',i10)
8890  format(/
     +1x,40('+-')/
     +15x,'print by symmetry',/
     +1x,80('-'))
8891  format(/
     +1x,40('+-')//)
c
      return
      end
      subroutine writer(w,ilifq,xcoff,nrootz,nrootci,mxroot)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      dimension w(*),ilifq(*),xcoff(*)
100   format(/20x,'*** Eigenstates *** '/
     +         1x,'Zero Order'/
     +         1x,'  Root    ',6i14//)
102   format(1x,i5,10x,6(f13.7,a1))
c
c     flag large coefficients for tagging 
c
      do i = 1, nrootci
       ilifq(i) = (i-1) * mxroot
      enddo
      do i = 1, nrootci
       ij = ilifq(i)
       do j = 1, nrootz
        xcoff(j+ij) = ' '
       enddo
       loop=idamax(nrootz,w(ij+1),1)
       wmax = abs(w(ij+loop)) / 3.0d0 
       do j = 1, nrootz
        if(abs(w(ij+j)).ge.wmax) then
        xcoff(ij+j) = '*'
        endif
       enddo
      enddo   
c
      lprnt = nrootci
      m=1
c
      m5=6
      n=m5
2106  if(lprnt.lt.m)return
      if(n.gt.lprnt)n=lprnt
      write(iwr,100)(i,i=m,n)
      do 2111 j=1,nrootz
      write(iwr,102)j,(w(j+ilifq(i)),
     +                   xcoff(j+ilifq(i)),i=m,n)
2111  continue
      m=m+m5
      n=n+m5
      goto 2106
      end
cvp
cvp konfingurationsvergleich mit integralbestimmung for dk=1
cvp configuration comparison with integral determination for dk=1
cvp
      subroutine check1a(twoe,nteint,iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,isk,jsk,ibs,nspiel,mspiel,jnum1
     &  ,nspp1,nsppa,nsppb)
cvp
cvp parameter und globale common-bloecke sind im include-file
cvp parameters und global common blocks are kept in include file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      integer nteint
      dimension twoe(nteint)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr,jconb,nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c speziell for dk=1
c jmerk laufen direkt vom laeufer angesteuert
c dadurch erfolgt der zugriff schneller
c geaenderte version: nur nummern der wechselwirkenden konf. zu einer
c bei gleichen sk erfolgt der vergleich bzgl. der unteren dreieckmatrix
c (vp 13.3.1992)
c especially for dk=1
c jmerk passes driven directly by the loop counter
c this speeds up the memory access
c modified version: just the numbers of interacting configurations
c in case of the same sk the comparison regards the lower triangle
c (vp 13.3.1992)
c  festen konf. werden bestimmt
c  uebergeben werden: iswh=anzahl der sk
c                     isk=linke sk
c                     jsk=rechte sk
c                     ibs=nummer der konf. aus der sk isk
c  resultate sind:    nspiel=anzahl der ww konf.
c                     ispiel=vektor, der die nummern der ww konf.
c                            enthaelt
c  determining permanent configurations
c  inputs are       : iswh=number of sk
c                     isk=left sk
c                     jsk=right sk
c                     ibs=number of configurations from the sk isk
c  outputs are   :    nspiel=number of interacting configurations
c                     ispiel=vector containing the numbers of the
c                            interacting configurations
cvp
cvp
cvp
      ncisk=nconf(isk)
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
        istar=niot(isk)+1
      iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
        do 215 k=1,imo
        jcon(k)=1
        jcon2(k)=2
        jcon1(k)=1
215     continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
         jcon(itest(k)) = 0
  220    continue
c
         do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
c first occupy jcon1, such that starting from a single and finding a 
c a double is not counted.
         jposc(itest(k)) = k-inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
cvp
         if(jsk.lt.1) go to 250
cvp
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
c building field ispiel ( nr. of players) eliminate all those
c that have been identified as double excitations already
c
cvp vergleich bzgl. unterer dreiecksmatrix
cvp comparison regarding the lower triangle 
cx       nspiel=ncjsk
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cvp here nspiel is the max. number of configurations to be tested 
cvp  that fit onto the vectors
         nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c ispiel : number of playing configurations within the sk
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c nullsetzen von idiff feld soweit notwendig
c nullify idiff field where needed
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c  von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using 
c creation and annihilation operators
c **  call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
c ** +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
c aufbau von jcona for vergleich auf basis von erzeugern
c building jcona for the comparison regarding generators
         do k=1,imo
c  for refcon bzw. refcon1
c  for refcon respectively refcon1
           jcona(k)=jcon1(k)
           jcona(k+momax)=jcon2(k)-jcon1(k)
         enddo
c
c aufbau von idref, d.h. vergleich der testkonfiguration ibs mit
c allen mains
c building idref, i.e. compare the test configuration ibs against
c all mains
c    berechnung der startadresse iostt for das iotnew-feld
c    calculating the start address iostt for the iotnew field
       iostt=iotnst(isk)
c    vorbelegung von idref
c    initialisation of idref
       mainr=iotnew(iostt+ibs)
       do k=1,nko
         idref(k)=maindf(k,mainr)
       enddo
c    korrektur von idref um die beitraege von jcona
c    correction of idref by the contributions from jcona
       moai=iotnew(iostt+  ncisk+ibs)
       moaj=iotnew(iostt+2*ncisk+ibs)
       moak=iotnew(iostt+3*ncisk+ibs)
       moal=iotnew(iostt+4*ncisk+ibs)
       do k=1,nko
         idref(k)=idref(k)+jconb(k,moai)+jconb(k,moaj)
     &                    -jconb(k,moak)-jconb(k,moal)
       enddo
c
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c
c berechnung von idiff
c calculating idiff
c    berechnung der startadresse iosttj for das iotnew-feld
c    calculation of the start address iosttj for the iotnew field
         iosttj =iotnst(jsk)+jnum1
         iostj1=iosttj + ncjsk
         iostj2=iostj1 + ncjsk
         iostj3=iostj2 + ncjsk
         iostj4=iostj3 + ncjsk
         do kdoo=1,nspiel
           idiff(kdoo)=idref( iotnew(iosttj+kdoo) )
     &                -jcona( iotnew(iostj3+kdoo) )
     &                -jcona( iotnew(iostj4+kdoo) )
         enddo
         do kdoo=1,nspiel
          if (idiff(kdoo).le.2) then
           idiff(kdoo)=idiff(kdoo)
     &                 +jcona( iotnew(iostj1+kdoo) )
     &                 +jcona( iotnew(iostj2+kdoo) )
          endif
         enddo
c
cbs
cbs jetzt wird komprimiert
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
             nspie1 = nspie1 + 1
cx           ispiel(nspie1) = kdoo
             ispiel(nspie1) = kdoo + jnum1
             idiff(nspie1)=idiff(kdoo)
cvp          if(idiff(kdoo).eq.1) then
cvp            nspp3=nspp3+1
cvp            ispp3(nspp3)=kdoo
cvp            npfal(nspie1)=4
cvp          endif
           endif
 385      continue
cvp jstar wird wieder auf alten wert zurueckgesetzt
cvp jstar is reset to the old value
         jstar=niot(jsk)
cvp
cvp  ende bestimmung der wechselwirkenden konfigurationen
cvp  end of interacting configuration determination
cvp
cvp
          nspiel=nspie1
cbs   und jetzt das ganze noch mal
cbs   and now do the whole thing once more
cvp   setzen von startwerten
cvp   initialisation
           do 261 jbs=1,nspiel
 261       jmerko(jbs)=0
cvp lauf ueber alle offenen schalen zur p-fall-trennung
cvp loop over all open shells to separate p-cases
cvp
c test der offenen schalen
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
c test the open shells
c like with the closed shells start with the last open shell
         kstar=jstar+jnopen*ncjsk
         do 3050 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c differenz bilden mit jcon
c compute difference with jcon
cvp jmerko zaehlt die besetzungsunterschiede bzgl. der einfach
cvp besetzten mo's von \phi_l
cvp jmerko counts the occupation differences with respect to the 
cvp singly occupied mo's of \phi_l
*vocl loop,novrec
          do 3040 kdoo = 1,nspiel
ct          jmerko(kdoo)=jmerko(kdoo)+abs(jcon( iottr(kstar+
ct   &                ispiel(kdoo)) )-1)
            jmerko(kdoo)=jmerko(kdoo)+jcon( iottr(kstar+
     &                ispiel(kdoo)) )
 3040     continue
 3050    continue
cvp bestimmung des p-falles durch jmerko und idiff
cvp notation nach buenker, d.h. npfal=2 --> p=1 usw.
cvp determination of p-case by jmerko and idiff
cvp notation following buenker, i.e. npfal=2 --> p=1, etc.
          nspp1=0
          nspp2=0
          nspp3=0
      do 3100 kdoo=1,nspiel
       if (jmerko(kdoo).eq.1) then
         nspp1=nspp1+1
         nwwmo(nspp1,1)=ispiel(kdoo)
cvp      npfal(kdoo)=2
cvp    endif
        else if (jmerko(kdoo).eq.0.and.idiff(kdoo).eq.2) then
         nspp2=nspp2+1
         nwwmo(nspp2,2)=ispiel(kdoo)
cvp      npfal(kdoo)=3
cvp    endif
c       else if (idiff(kdoo).eq.0) then
        else
         nspp3=nspp3+1
         nwwmo(nspp3,3)=ispiel(kdoo)
cvp      npfal(kdoo)=4
cvp    endif
cvp    if (jmerko(kdoo).eq.0) then
       endif
 3100 continue
cvp konfiguratonsnummern nach p-faellen sortiert --> npfal
cvp configuration numbers sorted according to p-cases --> npfal
      do 3210 kdoo=1,nspp1
        npfal(kdoo)=nwwmo(kdoo,1)
        jmerko(kdoo)=1
        jmerkd(kdoo)=0
        nqrfal(kdoo)=ibinom(inopen,3)
 3210 continue
      do 3220 kdoo=1,nspp2
        npfal(kdoo+nspp1)=nwwmo(kdoo,2)
        jmerkd(kdoo+nspp1)=0
        nqrfal(kdoo+nspp1)=ideks(inopen)
        moafal(kdoo+nspp1)=ideks(indub+1)
 3220 continue
      nsppa=nspp1+nspp2
      do 3230 kdoo=1,nspp3
        npfal(kdoo+nsppa)=nwwmo(kdoo,3)
        jmerko(kdoo+nsppa)=0
        jmerkd(kdoo+nsppa)=0
        nqrfal(kdoo+nsppa)=ideks(inopen)
        sumint(kdoo+nsppa)=0.0d0
 3230 continue
      nsppb=nsppa+nspp3
cvp schleife ueber alle doppelt besetzten mo's von phi_l
cvp loop over all doubly occupied mo's of \phi_l
         kstar=jstar+jnytl*ncjsk
c beginn hinten, da dort die meisten unterschiede zu finden sind
c start at end because most differences will be found there
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
cvp schleife for p=1
cvp loop for p=1
cvp jmerkd zaehlt die besetzungsunterschiede bzgl. der doppelt
cvp besetzten mo's von \phi_l
cvp jmerkd counts the occupation differences with respect to the 
cvp doubly occupied mo's of \phi_l
*vocl loop,novrec
          do 281 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              npos(kdoo,jmerkd(kdoo))=jposo(iottr(kstar+npfal(kdoo)))
            endif
 281      continue
cvp schleife for p=2
cvp loop for p=2
*vocl loop,novrec
          do 282 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              nwwmo(kdoo,jmerkd(kdoo))=iottr(kstar+npfal(kdoo))
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              moafal(kdoo)=moafal(kdoo)-nx
            endif
 282      continue
cvp schleife for p=3
cvp loop for p=3
*vocl loop,novrec
          do 283 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              npos(kdoo,1)=iottr(kstar+npfal(kdoo))
            endif
 283      continue
1260     continue
cvp
cvp lauf ueber die offenen schalen von hinten
cvp loop over the open shells starting at the end
cvp
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
cvp jmerko zaehlt die besetzungsunterschiede bzgl. der einfach
cvp besetzten mo's von \phi_l, wird hier wieder auf null herunter-
cvp gezaehlt
cvp jmerko counts the occupation differences with respect to the
cvp singly occupied mo's of \phi_l, is counted down to zero again here
cvp p=1-fall
cvp p=1-case
*vocl loop,novrec
          do 481 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            jmerko(kdoo)=0
            nwwmo(kdoo,1)=jnopen+1-jdoo
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo-jmerko(kdoo))
          endif
481       continue
cvp p=2-fall
cvp p=2-case
*vocl loop,novrec
          do 482 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
cvp       if(jdiff(kdoo).eq.0) then
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo)
cvp       endif
482       continue
cvp p=3-fall
cvp p=3-case
*vocl loop,novrec
          do 483 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
cvp       if(jdiff(kdoo).eq.0) then
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo)
cvp       endif
483       continue
 1360   continue
cvp bestimmung des r-falles for p=1 sowie des ql-falles
cvp determination of the r-case for p=1 like the ql-case
      do 491 kdoo=1,nspp1
        nqlfal(kdoo)=nwwmo(kdoo,1)
        nwwmo(kdoo,2)=imo3q(nqrfal(kdoo),1)
        nwwmo(kdoo,3)=imo3q(nqrfal(kdoo),2)
        nwwmo(kdoo,4)=imo3q(nqrfal(kdoo),3)
        if (npos(kdoo,1).eq.nwwmo(kdoo,2)) then
           npos(kdoo,1)=3
           if (npos(kdoo,2).eq.nwwmo(kdoo,3)) then
              nrfal(kdoo)=1
             else
              nrfal(kdoo)=2
           endif
           else if (npos(kdoo,1).eq.nwwmo(kdoo,3)) then
             npos(kdoo,1)=2
             nrfal(kdoo)=3
           else
           npos(kdoo,1)=1
        endif
        if (jmerkd(kdoo).eq.1) then
          nrfal(kdoo)=-npos(kdoo,1)
        endif
  491 continue
cvp p=2 fall
cvp p=2 case
cvp  codierung des r-falles durch vorzeichen von qr
cvp  the r-case is indicated by the sign of qr
      do 592 kdoo=nspp1+1,nsppa
        if (jmerkd(kdoo).eq.2)
     &    nqrfal(kdoo)=-nqrfal(kdoo)
  592 continue
cvp p=3 fall
cvp p=3 case
*vocl loop,novrec
      do 493 kdoo=nsppa+1,nsppb
cvp bestimmung der irred. darst. und abspeichern der hoeheren und
cvp der niedrigeren mo-nummer innerhalb der irred. darst.
cvp moafal enthaelt groesseres mo
cvp determination of the irrep and storing the higher and lower mo
cvp number within the irrep
cvp moafal holds the higher mo
cvp
        moafal(kdoo)=itest(imoq(nqrfal(kdoo),1))
        mobfal(kdoo)=itest(imoq(nqrfal(kdoo),2))
        if (npos(kdoo,1).eq.mobfal(kdoo))
     &     nqrfal(kdoo)=-1*nqrfal(kdoo)
        nrfal(kdoo)=nirred(mobfal(kdoo))
        moafal(kdoo)=mopos(moafal(kdoo))
        mobfal(kdoo)=mopos(mobfal(kdoo))
cvp codierung des r-falles for p=3 durch vorzeichen von qr
cvp the r-case for p=3 is indicated by the sign of qr
 493  continue
cvp
cvp
cvp  ende labelbestimmung
cvp  end of label determination
cvp
cvp
cvp keine berechnung der integraladressen, falls integrale
cvp bereits vorsortiert, d.h. ioint = 1
cvp no calculation of integral addresses when the integrals 
cvp have been presorted already, i.e. ioint = 1
      if (ioint.eq.1) goto 251
cvp
cvp bestimmung der ww mo's von \phi_r aus qr,
cvp berechnung und sortierung der absoluten mo-nummern for
cvp die integralberechnung
cvp und bestimmung des ql-falles aus den ww mo's von \phi_l
cvp determination of the interacting mo's of \phi_r in qr,
cvp calculation and sorting of the absolute mo-numbers for the
cvp integral calcalution 
cvp and the determination of the ql-cases from the interacting mo's
cvp of \phi_l
         kstar = niot(jsk) - ncjsk
*vocl loop,novrec
      do 591 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp sei das exchange-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff und jmerko sind hier frei verfuegbar
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp if given a exchange-integral in chemical-notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff and jmerko are freely available
        nwwmo(kdoo,1)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,1)*ncjsk)
        nwwmo(kdoo,2)=itest(nwwmo(kdoo,2))
        nwwmo(kdoo,3)=itest(nwwmo(kdoo,3))
        nwwmo(kdoo,4)=itest(nwwmo(kdoo,4))
cvp sortierung for cb- und ex-integrale
cvp sorting of the cb- and ex-integrals
          if (abs(nrfal(kdoo)).ne.3) then
cvp r=1,2,4,5
              if (nwwmo(kdoo,1).lt.nwwmo(kdoo,2)) then
                 nx1=nwwmo(kdoo,1)
                 nwwmo(kdoo,1)=nwwmo(kdoo,2)
                 nwwmo(kdoo,2)=nx1
                 idiff(kdoo)=nwwmo(kdoo,5-abs(nrfal(kdoo)))
                 jdiff(kdoo)=max(nwwmo(kdoo,2),
     &                           nwwmo(kdoo,2+abs(nrfal(kdoo))))
                 jmerko(kdoo)=min(nwwmo(kdoo,2),
     &                           nwwmo(kdoo,2+abs(nrfal(kdoo))))
                else
                 idiff(kdoo)=nwwmo(kdoo,2+abs(nrfal(kdoo)))
                 jdiff(kdoo)=nwwmo(kdoo,2)
                 jmerko(kdoo)=nwwmo(kdoo,5-abs(nrfal(kdoo)))
              endif
            else
cvp r=3 oder r=6
cvp r=3 or r=6
              nx2=nwwmo(kdoo,2)
              nwwmo(kdoo,2)=nwwmo(kdoo,3)
              nwwmo(kdoo,3)=nx2
              if (nwwmo(kdoo,1).lt.nwwmo(kdoo,3)) then
                 idiff(kdoo)=nwwmo(kdoo,2)
                 jdiff(kdoo)=max(nwwmo(kdoo,1),nwwmo(kdoo,4))
                 jmerko(kdoo)=min(nwwmo(kdoo,1),nwwmo(kdoo,4))
                 ny1=nwwmo(kdoo,1)
                 ny2=nwwmo(kdoo,2)
                 nwwmo(kdoo,1)=nwwmo(kdoo,3)
                 nwwmo(kdoo,2)=nwwmo(kdoo,4)
                 nwwmo(kdoo,3)=max(ny1,ny2)
                 nwwmo(kdoo,4)=min(ny1,ny2)
                else
                 idiff(kdoo)=nwwmo(kdoo,4)
                 jdiff(kdoo)=nwwmo(kdoo,3)
                 jmerko(kdoo)=nwwmo(kdoo,2)
              endif
          endif
c     write(6,*)
  591 continue
cvp
c     if (ibs.eq.155) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,3)
c3500 continue
c     endif
cvp p=2 fall
cvp p=2 case
cvp berechnung der absolutnummern der mo's und sortierung
cvp calculation of the absolute mo-numbers and sorting
*vocl loop,novrec
      do 492 kdoo=nspp1+1,nsppa
        if (jmerkd(kdoo).eq.2) then
cvp r=1
          nwwmo(kdoo,3)=itest(inopen+moafal(kdoo))
         else
cvp r=2
          nwwmo(kdoo,3)=nwwmo(kdoo,1)
          nwwmo(kdoo,1)=itest(imoq(nqrfal(kdoo),1))
          nwwmo(kdoo,2)=itest(imoq(nqrfal(kdoo),2))
        endif
cvp mo's sind nun wie folgt sortiert:
cvp  nwwmo(*,3) = doppelt besetztes mo, das in der anderen konf. leer
cvp               ist
cvp  nwwmo(*,1) = einfach besetztes mo
cvp  nwwmo(*,2) = einfach besetztes mo
cvp  nwwmo(*,1) > nwwmo(*,2)
cvp the mo's are now sorted according to:
cvp  nwwmo(*,3) = doubly occupied mo, that is empty in other 
cvp               configurations
cvp  nwwmo(*,1) = singly occupied mo
cvp  nwwmo(*,2) = singly occupied mo
cvp  nwwmo(*,1) > nwwmo(*,2)
  492 continue
cvp
cvp  berechnung der integraladresen
cvp  calculation of the integral addresses
cvp
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3201 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp zugriff auf nit-feld
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp access of the nit field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
cvp bestimmung des symmetrie-falles (cb-integral)
cvp determination of the symmetry case (cb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum2) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum3) then
cvp      fall 1
cvp      case 1
c           write(6,*) ' case 1'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks(ideks(idum1)+
     &        mopos(idiff(kdoo)))
     &        +ideks(mopos(jdiff(kdoo)))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 2
cvp      case 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))
     &       +ideks(ndum3+1) ) + idum4 +
     &       ideks(ncimo(ndum3)+1)*(ideks(idum1)+idum2-1)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &       +mopos(idiff(kdoo)) )
     &       +ncimo(nirred(idiff(kdoo)))
     &       *(mopos(jdiff(kdoo))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum3) then
cvp      fall 3
cvp      case 3
c           write(6,*) ' case 3'
            intcb(kdoo)=nit( ideks(idum+ndum2+1) )
     &       +idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
c           if (nrfal(kdoo).ne.2) then
            if (ndum1.ne.nirred(idiff(kdoo))) then
              intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &         +mopos(idiff(kdoo)) )
     &         +ncimo(nirred(idiff(kdoo)))
     &         *(mopos(jdiff(kdoo))-1)
             else
              intex(kdoo)=nit( ideks(ideks(ndum1+1))
     &         +ideks(nirred(jdiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks(ncimo(nirred(jdiff(kdoo)))+1)
     &         *(ideks(idum1)+mopos(idiff(kdoo))-1)
     &         +ideks(mopos(jdiff(kdoo)))
            endif
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 4
cvp      case 4
c           write(6,*) ' case 4'
cvp nur i-1 taucht in den formeln auf
cvp only i-1 shows up in the equations
            idum1=idum1-1
            intcb(kdoo)=nit( ideks(idum+ndum2)
     &       +ideks(ndum3)+ndum4 ) + idum4 +
     &       ncimo(ndum4)*(idum3-1+ncimo(ndum3)*(idum2-1+
     &       ncimo(ndum2)*idum1))
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo)))
     &       +ideks(nirred(jdiff(kdoo)))+nirred(jmerko(kdoo)) )
     &       + mopos(jmerko(kdoo)) +
     &       ncimo(nirred(jmerko(kdoo)))
     &       *(mopos(jdiff(kdoo))-1
     &       +ncimo(nirred(jdiff(kdoo)))
     &       *(mopos(idiff(kdoo))-1
     &       +ncimo(nirred(idiff(kdoo)))*idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
 3201 continue
cvp berechnung der integraladresse for p=2
cvp calculation of the integral addresses for p=2
*vocl loop,novrec
      do 3202 kdoo=nspp1+1,nsppa
cvp das integral sei (ac|bc) oder (ac|cb) usw.
cvp dann ist nwwmo(*,1)=max(a,b)
cvp          nwwmo(*,2)=min(a,b)
cvp          nwwmo(*,3)=c
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp given the integral (ac|bc) or (ac|cb) etc.
cvp then     nwwmo(*,1)=max(a,b)
cvp          nwwmo(*,2)=min(a,b)
cvp          nwwmo(*,3)=c
cvp it still must be sorted according to size though
        ndum1=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        ndum2=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep
        idum1=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum2=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum3=mopos(max(nwwmo(kdoo,2),nwwmo(kdoo,3)))
        idum4=mopos(min(nwwmo(kdoo,2),nwwmo(kdoo,3)))
cvp bestimmung des symmetrie-falles
cvp determination of the symmetry case
        if (ndum1.eq.ndum2) then
cvp        fall 1
cvp        case 1
c           write(6,*) ' case 1'
        intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &     + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
          else
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
        intcb(kdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &     + idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
        endif
 3202 continue
cvp
cvp  p=3
cvp
cvp  berechnung der summe uber cb- und ex-integrale sowie
cvp  der adressen der ex-integrale bei gemeinsamen offenen
cvp  schalen
cvp  calculation of the sum over cb- and ex-integrals and
cvp  the addresses of the ex-integrals with common open
cvp  shells
cvp
cvp  sei d gemeinsames offenes mo: betrachtung der integrale
cvp     (ab|dd) und (ad|bd), es gilt a > b
cvp  if d is a common open mo: considering the integrals
cvp     (ab|dd) and (ad|bd), it holds that a > b
cvp
      do 593 jdoo=1,jnopen
        do 613 kdoo=nsppa+1,nsppb
          mogd=iottr(kstar+npfal(kdoo)+jdoo*ncjsk)
cvp das ex-integral sei (ad|bd), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp given the ex-integral (ad|bd), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp it must still be sorted according to size though
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogd)
cvp bestimmung des symmetrie-falles
cvp determination of symmetry cases
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogd))
           idum2=min(moafal(kdoo),mopos(mogd))
           idum3=max(mobfal(kdoo),mopos(mogd))
           idum4=min(mobfal(kdoo),mopos(mogd))
cvp        fall 1
cvp        case 1
c           write(6,*) ' case 1'
            nwwmo(kdoo,jdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogd)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogd)
           idum3=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogd)
           idum4=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)
        endif
  613   continue
  593 continue
cvp schleife ueber die (gemeinsamen) geschlossenen schalen
cvp loop over the (common) closed shells
      do 533 jdoo=1,jndub
        do 623 kdoo=nsppa+1,nsppb
c         mogc=iot(iloc-1+inopen+jdoo)
          mogc=iottr(kstar+npfal(kdoo)+(jdoo+jnopen)*ncjsk)
cvp das ex-integral sei (ac|bc), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp given the ex-integral (ac|bc), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp it must still be sorted according to size though
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogc)
cvp bestimmung des symmetrie-falles
cvp determination of the symmetry-case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogc))
           idum2=min(moafal(kdoo),mopos(mogc))
           idum3=max(mobfal(kdoo),mopos(mogc))
           idum4=min(mobfal(kdoo),mopos(mogc))
cvp        fall 1
cvp        case 1
c           write(6,*) ' case 1'
            iex=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogc)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogc)
           idum3=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogc)
           idum4=mobfal(kdoo)
cvp        fall 3
cvp        case 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
cvp        case 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
            sumint(kdoo)=sumint(kdoo)+
     &        twoe(icb)+twoe(icb)-twoe(iex)
        endif
  623   continue
  533 continue
cvp korrektur der integralsumme
cvp correction of the integral sum
c       do 633 kdoo=nsppa+1,nsppb
c         if (nqrfal(kdoo).lt.0) then
cvp zugriff auf nit-feld
cvp access of nit field
c           ndum1=nrfal(kdoo)
c           idum1=iot(iloc-1+idiff(kdoo))
c           ix=ideks(mopos(idum1)+1)
c           iy=ideks(moafal(kdoo))+mobfal(kdoo)
c           icb=nit( ideks ( ideks(ndum1+1)+1 ) )
c    &       +ideks( max(ix,iy) )
c    &       +       min(ix,iy)
c           sumint(kdoo)=sumint(kdoo)+twoe(icb)
c         endif
c 633   continue
cvp
cvp  ende berechnung der integraladressen
cvp  end of integral address calculation
cvp
 251  continue
cvp
c     if (ibs.eq.4.or.ibs.eq.42) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c       write(6,*) ' intcb= ',intcb(kdoo)
c3500 continue
c     endif
cvp
cvp ruecksortierung wird uebergangen !!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel enthaelt jetzt die adressen der wechselwirkenden
cvp  konfigurationen: zuerst die for p=1, dann for p=2 usw.
cvp  auf npfal wird vorlaeufig der p-fall abgespeichert
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp back-sort is skipped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel now holds the addresses of the interacting
cvp  configurations: first those for p=1, then for p=2 usw.
cvp  provisionaly the p-case is stored on npfal
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 5210 kdoo=1,nspp1
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=2
 5210 continue
      do 5220 kdoo=nspp1+1,nsppa
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=3
 5220 continue
      do 5230 kdoo=nsppa+1,nsppb
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=4
 5230 continue
cvp
cvp p-faelle werden in urspruengliche reihenfolge gebracht
cvp p-cases brought back to original ordering
cvp
cre    lspp1=1
cre    lspp2=1
cre    lspp3=1
cvp
c      write(6,*) ' jbs=,nspp3=,ispp3= ',jbs,nspp3,(ispp3(kkkk),
c    &   kkkk=1,nspp3)
cvp
cre    do 5100 kdoo=1,nspiel
cre      if (ispiel(kdoo).eq.npfal(lspp1).and.lspp1.le.nspp1) then
cre        ndfal1(kdoo)=2
cre        ndfal2(kdoo)=nrfal(lspp1)
cre        ndfal3(kdoo)=nqlfal(lspp1)
cre        ndfal4(kdoo)=nqrfal(lspp1)
cre        ndfal7(kdoo)=intcb(lspp1)
cre        ndfal8(kdoo)=intex(lspp1)
cre        lspp1=lspp1+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp2+nspp1).and.lspp2.le.nspp2)
cre  &  then
cre        ndfal1(kdoo)=3
cre        ndfal4(kdoo)=nqrfal(lspp2+nspp1)
cre        ndfal7(kdoo)=intcb(lspp2+nspp1)
cre        lspp2=lspp2+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp3+nsppa).and.lspp3.le.nspp3)
cre  &  then
c          write(6,*) ' kdoo=,ispiel=,ispp3= ',kdoo,ispiel(kdoo)
c    &       ,ispp3(lspp3)
cre        ndfal1(kdoo)=4
cre        ndfal2(kdoo)=nrfal(lspp3+nsppa)
cre        ndfal4(kdoo)=nqrfal(lspp3+nsppa)
cre        ndfal5(kdoo)=moafal(lspp3+nsppa)
cre        ndfal6(kdoo)=mobfal(lspp3+nsppa)
cre        ndfeld(kdoo,1)=nwwmo(lspp3+nsppa,1)
cre        ndfeld(kdoo,2)=nwwmo(lspp3+nsppa,2)
cre        ndfeld(kdoo,3)=nwwmo(lspp3+nsppa,3)
cre        ndfeld(kdoo,4)=nwwmo(lspp3+nsppa,4)
cre        ndfeld(kdoo,5)=nwwmo(lspp3+nsppa,5)
cre        ndfeld(kdoo,6)=nwwmo(lspp3+nsppa,6)
cre        ndfeld(kdoo,7)=nwwmo(lspp3+nsppa,7)
cre        ndfeld(kdoo,8)=nwwmo(lspp3+nsppa,8)
cre        ndfeld(kdoo,9)=nwwmo(lspp3+nsppa,9)
cre        rdum(kdoo)=sumint(lspp3+nsppa)
cre        lspp3=lspp3+1
cre      endif
c5100  continue
 250  continue
      return
      end
c
c
      subroutine check2(iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,isk,jsk,ibs,nspiel,mspiel,jnum1)
cvp
cvp parameter und globale common-bloecke sind im include-file
cvp parameters and global common blocks are in the include file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
       integer iottr, iotnew, iotm
       dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr,jconb,nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c speziell for dk=2
c jmerk laufen direkt vom laeufer angesteuert
c dadurch erfolgt der zugriff schneller
c geaenderte version: nur nummern der wechselwirkenden konf. zu einer
c bei gleichen sk erfolgt der vergleich bzgl. der unteren dreieckmatrix
c (vp 13.3.1992)
c especially for dk=2
c jmerk passes driven directly by the loop counter
c this speeds up the memory access
c modified version: just the numbers of interacting configurations
c in case of the same sk the comparison regards the lower triangle
c (vp 13.3.1992)
c  festen konf. werden bestimmt
c  uebergeben werden: iswh=anzahl der sk
c                     isk=linke sk
c                     jsk=rechte sk
c                     ibs=nummer der konf. aus der sk isk
c  resultate sind:    nspiel=anzahl der ww konf.
c                     ispiel=vektor, der die nummern der ww konf.
c                            enthaelt
c  determining permanent configurations
c  inputs are       : iswh=number of sk
c                     isk=left sk
c                     jsk=right sk
c                     ibs=number of configurations from the sk isk
c  outputs are   :    nspiel=number of interacting configurations
c                     ispiel=vector containing the numbers of the
c                            interacting configurations
cvp
      ncisk=nconf(isk)
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
        istar=niot(isk)+1
      iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
        do 215 k=1,imo
ct      jcon1(k)=4
        jcon1(k)=1
        jcon2(k)=2
215     continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
  220    continue
cvp
c        write(6,*) ' sk=',isk
c        write(6,*) ' phi_r: ',(iot(k),k=iloc,iloc+inytl-1)
cvp
         do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
c first occupy jcon1, such that starting from a single and finding a
c a double is not counted.
         jposc(itest(k)) = k-inopen
ct       jcon1(itest(k)) = 4
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
c        write(6,*) 'jcon for konf#',ibs,' in sk ',isk
c        write(6,*) (jcon(iaus),iaus=1,imo)
c beginn do loop ueber konfig aus ft36, gleich und sk+1 ist interessant
c begin do over configurations from ft36. The same and sk+1 
c is interesting
         if(jsk.lt.1) go to 250
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
c building field ispiel ( nr. of players) eliminate all those
c that have been identified as double excitations already
c
cvp vergleich bzgl. unterer dreiecksmatrix
cvp comparison regarding the lower triangle
cx       nspiel=ncjsk
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cvp here nspiel is the max. number of configurations to be tested
cvp  that fit onto the vectors

         nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c ispiel : number of playing configurations within the sk
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using
c creation and annihilation operators
*     call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
*    +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
c
c aufbau von jcona for vergleich auf basis von erzeugern
c building jcona for the comparison regarding generators
         do k=1,imo
c  for refcon bzw. refcon1
c  for refcon respectively refcon1
           jcona(k)=jcon1(k)
           jcona(k+momax)=jcon2(k)-jcon1(k)
         enddo
c
c aufbau von idref, d.h. vergleich der testkonfiguration ibs mit
c  allen mains
c building idref, i.e. compare the test configuration ibs against
c all mains
c    berechnung der startadresse iostt for das iotnew-feld
c    calculating the start address iostt for the iotnew field
       iostt=iotnst(isk)
c  vorbelegung von idref
c    initialisation of idref
       mainr=iotnew(iostt+ibs)
       do k=1,nko
         idref(k)=maindf(k,mainr)
       enddo
c  korrektur von idref um die beitraege von jcona
c    correction of idref by the contributions from jcona
       moai=iotnew(iostt+  ncisk+ibs)
       moaj=iotnew(iostt+2*ncisk+ibs)
       moak=iotnew(iostt+3*ncisk+ibs)
       moal=iotnew(iostt+4*ncisk+ibs)
       do k=1,nko
         idref(k)=idref(k)+jconb(k,moai)+jconb(k,moaj)
     &                    -jconb(k,moak)-jconb(k,moal)
       enddo
c
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c
c berechnung von idiff
c calculating idiff
c    berechnung der startadresse iosttj for das iotnew-feld
c    calculation of the start address iosttj for the iotnew field
         iosttj =iotnst(jsk)+jnum1
         iostj1=iosttj + ncjsk
         iostj2=iostj1 + ncjsk
         iostj3=iostj2 + ncjsk
         iostj4=iostj3 + ncjsk
         do kdoo=1,nspiel
           idiff(kdoo)=idref( iotnew(iosttj+kdoo) )
     &                -jcona( iotnew(iostj3+kdoo) )
     &                -jcona( iotnew(iostj4+kdoo) )
         enddo
         do kdoo=1,nspiel
          if (idiff(kdoo).le.2) then
           idiff(kdoo)=idiff(kdoo)
     &                 +jcona( iotnew(iostj1+kdoo) )
     &                 +jcona( iotnew(iostj2+kdoo) )
          endif
         enddo
c
cbs
cbs jetzt wird komprimiert
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
            nspie1 = nspie1 + 1
cx         ispiel(nspie1) = kdoo
           ispiel(nspie1) = kdoo + jnum1
           idiff(nspie1)=idiff(kdoo)
           endif
 385      continue
          nspiel=nspie1
cvp jstar wird wieder auf alten wert zurueckgesetzt
cvp jstar is reset to the old value
         jstar=niot(jsk)
cvp
cvp  ende bestimmung der wechselwirkenden konfigurationen
cvp  end of interacting configuration determination
cvp
cbs   und jetzt das ganze noch mal
cbs   and now do the whole thing once more
           do 261 jbs=1,nspiel
           nqrfal(jbs)=ibinom(inopen,4)
 261       jmerkd(jbs)=0
c beginn hinten, da dort die meisten unterschiede zu finden sind
c start at end because most differences will be found there
         kstar=jstar+jnytl*ncjsk
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
c unterschiede zur testkonfiguration auf idiff summieren
c sum differences with respect to test configuration onto idiff
cvp jdiff(maximal) ist for doppelt besetzte mo's genau 2 for dk=2
cvp npos(*,1) gibt position des hinteren unterschiedes, npos(*,2) des
cvp vorderen
cvp jdiff(maximal) is exactly 2 for doubly occupied mo's with dk=2
cvp npos(*,1) gives the position of the last difference, 
cvp npos(*,2) the position of the first
*vocl loop,novrec
          do 286 kdoo=1,nspiel
ct          jdiff(kdoo)=jcon( iottr(kstar+ispiel(kdoo)))-2
            jdiff(kdoo)=jcon2( iottr(kstar+ispiel(kdoo)))
            if(jdiff(kdoo).ne.0) then
             jmerkd(kdoo)=jmerkd(kdoo)+1
             npos(kdoo,jmerkd(kdoo))=jposo(iottr(kstar+ispiel(kdoo)))
            endif
 286      continue
1260     continue
cvp
cvp bestimmung der in \phi_r einfach und in \phi_l nicht besetzten
cvp determination of the in \phi_r singly occupied and in \phi_l 
cvp unoccupied mo's
c test der offenen schalen
c test of open shells
c        write(6,*) 'jcon for konf#',ibs,' in sk ',isk
c        write(6,*) (jcon(iaus),iaus=1,imo)
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
c like with the closed shells start at the last open shell
cvp nx = position der gemeinsamen offenen schale (letzte offene von
cvp  \phi_l)
cvp nx = position of the common open shell (last open shell of \phi_l)
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c differenz bilden mit jcon
c compute difference with jcon
cvp nx = position der gemeinsamen offenen schale
cvp nx = position of the common open shell
*vocl loop,novrec(jmerkd),vi(kdoo)
          do 1380 kdoo = 1,nspiel
          nx=jposo(iottr(kstar+ispiel(kdoo)))
cvp bestimmung des qr-falles
cvp determination of the qr-case
          nqrfal(kdoo)=nqrfal(kdoo)-
     &     ibinom(nx-1,jnopen+1-jdoo)
          if (nx.lt.npos(kdoo,1)) npos(kdoo,1)=npos(kdoo,1)-1
          if (nx.lt.npos(kdoo,2)) npos(kdoo,2)=npos(kdoo,2)-1
 1380     continue
 1360     continue
cvp
cvp bestimmung des r-falles aus npos
cvp determination of the r-case from npos
      do 3010 kdoo=1,nspiel
        nrfal(kdoo)=npos(kdoo,1)-npos(kdoo,2)
        if (nrfal(kdoo).eq.1.and.npos(kdoo,1).eq.3) nrfal(kdoo)=3
cvp     nrfal(kdoo)=idk2r(npos(kdoo,1),npos(kdoo,2))
 3010 continue
cvp
cvp  ende labelbestimmung
cvp  end of labeldetermination
cvp
cvp keine berechnung der integraladressen, falls integrale
cvp  bereits vorsortiert, d.h. ioint = 1
cvp no calculation of integral addresses when the integrals
cvp have been presorted already, i.e. ioint = 1
      if (ioint.eq.1) goto 251
cvp
cvp  berechnung der integraladresen
cvp  integral address calculation
cvp
      iloc1=iloc-1
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3200 kdoo=1,nspiel
cvp uebergang von positionen auf absolute mo-nummern
cvp change over from positions to absolute mo-numbers
        nwwmo(kdoo,1)=itest(imo4q(nqrfal(kdoo),1))
        nwwmo(kdoo,4)=itest(imo4q(nqrfal(kdoo),4))
        if (nrfal(kdoo).eq.1) then
          nwwmo(kdoo,3)=itest(imo4q(nqrfal(kdoo),3))
          nwwmo(kdoo,2)=itest(imo4q(nqrfal(kdoo),2))
        endif
        if (nrfal(kdoo).ne.1) then
          nwwmo(kdoo,3)=itest(imo4q(nqrfal(kdoo),2))
          nwwmo(kdoo,2)=itest(imo4q(nqrfal(kdoo),3))
        endif
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,3)
cvp                            b = nwwmo(*,2)
cvp                            a = nwwmo(*,4)
cvp zugriff auf nit-feld
cvp given the coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,3)
cvp                            b = nwwmo(*,2)
cvp                            a = nwwmo(*,4)
cvp access nit field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
        intcb(kdoo)=nit( ideks(idum+ndum3)
     &    +ideks(ndum2)+ndum4 ) + idum4
        intex(kdoo)=nit( ideks(idum+
     &    nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &    +ideks(nirred(nwwmo(kdoo,nrb(nrfal(kdoo)))))
     &    +nirred(nwwmo(kdoo,nrc(nrfal(kdoo)))) )
     &    +mopos(nwwmo(kdoo,nrc(nrfal(kdoo))))
cvp bestimmung des symmetrie-falles (cb-integral)
cvp determination of the symmetry case (cb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum3) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum2) then
cvp      fall 1
c           write(6,*) ' case 1'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ideks(idum1)+idum3)
     &       +ideks(idum2)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ideks(idum1)+
     &        mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &        +ideks(mopos(nwwmo(kdoo,nrb(nrfal(kdoo)))))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum2)+1)*(ideks(idum1)+idum3-1)
     &       +ideks(idum2)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(idum1-1)
     &       +mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum2) then
cvp      fall 3
c           write(6,*) ' case 3'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum3)*(idum1-1)+idum3)
     &       +ncimo(ndum3)*(idum2-1)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(idum1-1)
     &       +mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 4
c           write(6,*) ' case 4'
cvp nur i-1 taucht in den formeln auf
            idum1=idum1-1
            intcb(kdoo)=intcb(kdoo)+
     &       ncimo(ndum4)*(idum2-1+ncimo(ndum2)*(idum3-1+
     &       ncimo(ndum3)*idum1))
            intex(kdoo)=intex(kdoo)+
     &       ncimo(nirred(nwwmo(kdoo,nrc(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1
     &       +ncimo(nirred(nwwmo(kdoo,nrb(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nra(nrfal(kdoo))))-1
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
c     write(6,*) ' ibs=,kdoo=,intex= ',ibs,kdoo,intex(kdoo)
 3200 continue
cvp   do 3500 kdoo=1,nspiel
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+ispiel(kdoo)*jnytl,
c    &   jstar+ispiel(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c3500 continue
cvp
cvp  ende berechnung der integraladressen
cvp  end of integral address calculation
cvp
  251 continue
cvp
  250 continue
cvp
      return
      end
c
cvp aufbau von feldern, die konfigurationsvergleich erleichtern
cvp sollen
cvp build fields that should facilitate configuration comparison
      subroutine chfeld(iottr,iot,iotm)
cvp
cvp parameter und globale common-bloecke sind im include-file
cvp parameters and global common blocks are in the include file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer iottr, iot, iotm
      dimension iottr(iotm), iot(iotm)
cvp
cvp iot wird hier als dummy-feld benutzt, um das konfigurationsfeld
cvp  iottr zu transponieren und dann wieder auf iottr abzulegen
cvp iot is used here as a dummy field to transpose the configuration
cvp field and then store is on iottr
cvp
cvp
      do 10 i=1,4
        do 10 j=1,3
   10     idk2r(i,j)=0
      do 20 i=2,4
        do 20 j=1,i-1
   20     idk2r(i,j)=i-j
      idk2r(3,2)=3
cvp
      nra(1)=4
      nra(2)=4
      nra(3)=2
      nrb(1)=2
      nrb(2)=3
      nrb(3)=3
      nrc(1)=3
      nrc(2)=2
      nrc(3)=4
cvp
cvp aufbau des ibinom-feldes
cvp ibinom enthaelt die binomialkoeffizienten
cvp build the ibinom field
cvp ibinom holds the binomial coefficients
cvp (i uber 0) = 1
cvp (i over 0) = 1
      do 110 i=0,nopmax+1
  110 ibinom(i,0)=1
      do 115 i=1,nopmax+1
      do 115 j=1,nopmax+1
  115 ibinom(i,j)=0
cvp (0 uber 1) = 0
cvp (0 over 1) = 0
      do 120 i=0,nopmax+1
  120 ibinom(i,1)=i
      do 130 i=2,nopmax+1
      do 130 j=2,i
  130 ibinom(i,j)=ibinom(i-1,j)+ibinom(i-1,j-1)
cvp
cvp belegung der zuordnung q-fall <--> mo-nummern for dk=0,p=1
cvp iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
cvp folgt der q-fall
cvp imoq(q-fall,1): hoeheres mo
cvp imoq(q-fall,2): niedrigeres mo
cvp
cvp layout of q-case assignment <--> mo-numbers for dk=0,p=1
cvp iqmo(i,j): i=higher mo, j=lower mo, from the mo-numbers 
cvp            the q-case follows
cvp imoq(q-fall,1): higher mo
cvp imoq(q-fall,2): lower mo
      do 140 i=2,nopmax
      do 140 j=1,i-1
      iqmo(i,j)=ibinom(i,2)-i+j+1
      imoq(iqmo(i,j),1)=i
  140 imoq(iqmo(i,j),2)=j
cvp
cvp belegung der zuordnung q-fall --> mo-nummern for dk=1,p=1
cvp imo3q(q-fall,1): hoehstes mo
cvp imo3q(q-fall,2): mittleres mo
cvp imo3q(q-fall,3): niedrigstes mo
cvp
cvp layout of q-case assignment <--> mo-numbers for dk=1,p=1
cvp imo3q(q-fall,1): highest mo
cvp imo3q(q-fall,2): middle mo
cvp imo3q(q-fall,3): lowest mo
      do 160 i=3,nopmax
      do 160 j=2,i-1
      do 160 k=1,j-1
      ndum=ibinom(i-1,3)+ibinom(j-1,2)+k
      imo3q(ndum,1)=i
      imo3q(ndum,2)=j
  160 imo3q(ndum,3)=k
cvp
cvp belegung der zuordnung q-fall --> mo-nummern for dk=2
cvp imo4q(q-fall,1): hoehstes mo
cvp imo4q(q-fall,2): zweites mo
cvp imo4q(q-fall,3): drittes mo
cvp imo4q(q-fall,4): niedrigstes mo
cvp
cvp layout of the q-case assignment --> mo-numbers for dk=2
cvp imo4q(q-fall,1): highest mo
cvp imo4q(q-fall,2): highest mo but one
cvp imo4q(q-fall,3): lowest mo but one
cvp imo4q(q-fall,4): lowest mo
      do 170 i=4,nopmax
      do 170 j=3,i-1
      do 170 k=2,j-1
      do 170 l=1,k-1
      ndum=ibinom(i-1,4)+ibinom(j-1,3)+ibinom(k-1,2)+l
      imo4q(ndum,1)=i
      imo4q(ndum,2)=j
      imo4q(ndum,3)=k
  170 imo4q(ndum,4)=l
cvp
cvp  transposition des iot-feldes
cvp  transposition of the iot field
cvp
      do 200 i=1,iswhm
        istar=niot(i)
        do 210 j=1,nytl(i)
          do 220 k=1,nconf(i)
            iot(istar+k+(j-1)*nconf(i))=
     &           iottr(istar+j+(k-1)*nytl(i))
  220     continue
  210   continue
  200 continue
cvp umspeichern auf iottr
cvp restore onto iottr
      do 230 i=1,iotm
  230 iottr(i)=iot(i)
cvp
      return
      end
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c  von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using
c creation and annihilation operators
      subroutine chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +      isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
cvp
cvp parameter und globale common-bloecke sind im include-file
cvp parameters and global common blocks are in the include file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer iotnew, iotm
      dimension iotnew(iotm)
      integer maindf, maxr,jconb,nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c
c aufbau von jcona for vergleich auf basis von erzeugern
c build jcona for comparison with respect to generators
*         write(6,*) 'call of chkww'
*         write(6,*) 'imo ',imo
         do 224 k=1,imo
c  for refcon bzw. refcon1
c  for refcon respectively refcon1
           jcona(k)=jcon1(k)
           jcona(k+momax)=jcon2(k)-jcon1(k)
  224    continue
c
c aufbau von idref, d.h. vergleich der testkonfiguration ibs mit
c  allen mains
c building idref, i.e. compare the test configuration ibs against
c all mains
c    berechnung der startadresse iostt for das iotnew-feld
c    calculating the start address iostt for the iotnew field
       iostt=iotnst(isk)
c  vorbelegung von idref
c    initialisation of idref
       mainr=iotnew(iostt+ibs)
*       write(6,*) 'nko',nko
       do 227 k=1,nko
         idref(k)=maindf(k,mainr)
  227  continue
c  korrektur von idref um die beitraege von jcona
c    correction of idref by the contributions from jcona
       moai=iotnew(iostt+  ncisk+ibs)
       moaj=iotnew(iostt+2*ncisk+ibs)
       moak=iotnew(iostt+3*ncisk+ibs)
       moal=iotnew(iostt+4*ncisk+ibs)
       do 262 k=1,nko
         idref(k)=idref(k)+jconb(k,moai)+jconb(k,moaj)
     &                    -jconb(k,moak)-jconb(k,moal)
  262  continue
c
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c
c berechnung von idiff
c calculating idiff
c    berechnung der startadresse iosttj for das iotnew-feld
c    calculation of the start address iosttj for the iotnew field
         iosttj =iotnst(jsk)+jnum1
         iostj1=iosttj + ncjsk
         iostj2=iostj1 + ncjsk
         iostj3=iostj2 + ncjsk
         iostj4=iostj3 + ncjsk
         do 375 kdoo=1,nspiel
           idiff(kdoo)=idref( iotnew(iosttj+kdoo) )
     &                -jcona( iotnew(iostj3+kdoo) )
     &                -jcona( iotnew(iostj4+kdoo) )
  375    continue
         do 376 kdoo=1,nspiel
          if (idiff(kdoo).le.2) then
           idiff(kdoo)=idiff(kdoo)
     &                 +jcona( iotnew(iostj1+kdoo) )
     &                 +jcona( iotnew(iostj2+kdoo) )
          endif
  376    continue
      return
      end
      subroutine clockv(r1,r2)
c
c  simulation von clockv via cputim
c
      implicit none
      real*8 r1,r2,r3
      common /told/ r3
c -- alte zeit
      r1 = r3
      call cputim(r2)
      r3 = r2
      return
      end
c
      subroutine david1a(twoe, nteint,
     +   pey, acoul, aexc, ndeke, core,
     +   vecf1, vecf2, mdi, ed, mdi0,
     +   emat, nedim,
     +   iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ncall,idim,adiag,v,w,
     +   av, h, u, mx1, mx,odebug)
c
c    bereite den den call zur davidsondiagonalisierung
c    mit der liu modifikation
c    prepare the call to the davidson diagonalisation
c    with the liu modification
c
c    auf ft20 ::
c    on ft20 ::
cbe  e auf mxroot gesetzt
cbe  e set to mxroot
c
      implicit none
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c NB
c cepai: if cepa calculation is requested
c icepa: specifies the cepa variant
c ipcepa: print level
c 
      logical cepai
      integer icepa,ipcepa
      common/cepa_mrd1/cepai,icepa,ipcepa
c
c NB
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      real*8 twoe
      integer nteint, ndeke
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      integer mdi, mdi0
      real*8 vecf1, vecf2, ed
      dimension vecf1(mdi), vecf2(mdi), ed(mdi0)
c
      real*8 adiag,w,v
      integer idim
      dimension adiag(idim)
      dimension w(idim),v(idim)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      real*8 av, h, u
      integer mx, mx1
      dimension av(mx1), h(mx1), u(mx,mx)
c
      logical odebug
c
      integer ncall,mxroot
      integer luns,iraus
      integer i,ifirst,ilast,nivnl
*     integer iswhm
*     parameter (iswhm=5)
      parameter (mxroot=50)
c   startvektor
c     arrays for den davson
*     integer mxref
*     parameter (mxref=256)
*     integer iwurz(mxref)
      integer iselec
c --- diagonale
      real*8 e
      dimension e(mxroot)
      real*8 critan,critrn,criten,orthon
      common /thresh/ critan,critrn,criten,orthon
      integer nroot,ntch,kprin,ndeci
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                ndeci
      common /cselec/ iselec(mxroot)
c --- common for extrapolation
      real*8 esav0,esav1,csum,ueb0,ueb1,de0,de1,egeys
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .              de0(mxroot),de1(mxroot),
     .              ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .              egeys(mxroot)
c --- scf information
      real*8 escf,vnuc
      common /kerne/ escf,vnuc
      real*8 egey1,trash,tdel
      common /parkin/egey1,trash,tdel
c --- thresholdson
c --- alter common
c     common /c/ idim
c
c --- luns : startfile
c
      luns  = lun20
c
*     ifirst=1
*     ilast=1
      nivnl=ilast
c --- thresholds
clear      critan = 0.00001d0
clear      critrn = 0.00001d0
clear      criten = 0.00001d0
clear      orthon = 0.1d-6
c --- call to dvdmua
      if(odebug) then
       write(6,*)'david1a: call to dvdmua : iselec',
     +         (iselec(i),i=1,mxroot)
      endif
      call dvdmua(twoe, nteint,
     +    pey, acoul, aexc, ndeke, core,
     +    vecf1, vecf2, mdi, ed, mdi0,
     +    emat, nedim,
     +    iottr,iotnew, iotm,
     +    maindf,maxr,jconb,nomax2,
*   idim:  actual dimension of the matrix
     .  idim,
*   first and last root : minl,mfnl
     .  ifirst, ilast,
*   number of initial vectors     iselect(30): which root!
*   iwurz : blocking of simultaneous roots
csut .  nivnl, iwurz, iselec,
     .  nivnl,
*       noutnl, lunnl, lun1nl, lun2nl,  : files out,f20, scratch
     .    iwr,  luns , lun21 , lun22,
*       critan, critrn, criten, orthon, : thresholds
     .  critan, critrn, criten, orthon,
*       adiag, v, w, e) : diagonalmatrix, v(idim) : last vector
*                         e : e(30) : krylov-egenwerte
*                         w : w(idim) : scratch-vector
     .  adiag, v, w, e,
     +  av, h, u, mx1, mx,iwr,odebug)
c
c  --- drucke krylov-matrix
c
       write(6,*)
       if (ncall.lt.1) then
c NB cepa prints
         if (cepai) then
            write(6,*) 'cepa-eigenvalues : '
         else
            write(6,*) 'ci-eigenvalues : '
         end if
c NB
cbbe        do i=ifirst,ilast
        do i=1,nroot
         iraus = iselec(i)
         write(6,'(12x,i2," energy :",f16.10)') iraus,e(iraus)+escf
         esav0(i) = e(iraus)+escf
        end do
       else
        write(6,1601)trash+tdel
1601    format(1x,'ci-eigenvalue (T =',f5.1,' uH threshold):')
cbe        do i=ifirst,ilast
        do i=1,nroot
         iraus = iselec(i)
         write(6,*) iraus,' energy :',e(iraus)+escf
         esav1(i) = e(iraus) + escf
        end do
      end if
      return
      end
      subroutine debvec(n,thresh,v)
      implicit real*8 (a-h,o-z)
      dimension v(n)
      write(6,*) 'out debvec.f : all larger than ',thresh
      do   ilu = 1,n
        if(abs(v(ilu)).gt.thresh) then
          write(6,*) ilu,v(ilu)
        endif
      enddo
      return
      end
      subroutine diag(twoe, nteint,
     +  pey, acoul, aexc, ndeke, core,
     +  vecf1, vecf2, mdi,
     +  hp5, adiag, v, w, mdi0,
     +  emat,nedim,iottr,iotnew,iotm,
     +  maindf,maxr,jconb,nomax2,
     +  av, h, u, mx2, mx1, mx,odebug,
     +  av3, u3, mx3, mx31, mx32)
c
c  main driver routine for diagonalization
c
c   (1993) h.u.s.
c
      real*8 twoe
      integer nteint
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      integer mdi, mdi0
      real*8 vecf1, vecf2, hp5
      real*8 adiag, v, w
      dimension vecf1(mdi), vecf2(mdi)
      dimension hp5(mdi0), adiag(mdi0), v(mdi0), w(mdi0)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      real*8 av, h, u
      integer mx, mx1
      dimension av(mx1), h(mx2), u(mx,mx)
c
      logical odebug
c
c--- bereite den start-vektor vor
c--- prepare the start vector
      integer ndav,ncall
      integer idim,i
      common /c/ idim
      common /cdiag_adler/ndav
c
      integer mx3, mx31, mx32
      real*8 av3, u3
      dimension av3(mx31),u3(mx3,mx3)
c ----
c --- mache den startvektor
c --- make the start vector
c     need to redefine dimensions between newmrd3 and newmrd5
c
      call prep0(adiag,mdi0,u3,av3,mx32,mx31,odebug)
      call vclr(u,1,mx2)
      call vclr(av,1,mx1)
      do i=1,mx3
       do j=1,mx3
       u(i,j)= u3(i,j)
       enddo
      enddo
      do i=1,mx31
       av(i) = av3(i)
      enddo
c --- kopiere die diagonale
c --- copy the diagonal
      call zro1d(idim,v)
      call zro1d(idim,w)
      do i=1,idim
         adiag(i) = hp5(i)
      end do
c steuerung per hand ob multi-davidson oder normaler davidson
cbe im readin ist sichergestellt dass nur david1.f aufgerufen wird
c ndav = 1    multi
c ndav = 2    engels
c control by hand whether the multi-davidson or normal davidson is used
cbe in readin it is ensured that only david1.f is called
c ndav = 1    multi
c ndav = 2    engels
c
c anzahl aufrufe for extrapolation
c ncall=0 : erster aufruf
c number of calls for extrapolation
c ncall=0 : first call
      ncall =0
c     this version of the code only supports ndav = 1 (multi)
c     remove other code for the time being
        call david1a(twoe, nteint,
     +       pey, acoul, aexc, ndeke, core,
     +       vecf1, vecf2, mdi, hp5, mdi0,
     +       emat, nedim,
     +       iottr,iotnew,iotm,
     +       maindf,maxr,jconb,nomax2,
     +       ncall,idim,adiag,v,w,
     +       av, h, u, mx1, mx,odebug)
c
      return
      end
      subroutine diag2(twoe, nteint,
     +         pey, acoul, aexc, ndeke, core,
     +         vecf1, vecf2, mdi,
     +         hp5, adiag, v, w, mdi0,
     +         emat,nedim,iottr,iotnew,iotm,
     +         maindf, maxr,jconb,nomax2,
     +         av, h, u, mx1, mx, odebug)
c
c  main driver routine for diagonalization
c
c   (1993) h.u.s.
c
c
c--- bereite den start-vektor vor
c--- prepate the start vector
      implicit none
c
      real*8 twoe
      integer nteint, ndeke
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      integer mdi,mdi0
      real*8 vecf1, vecf2
      real*8 hp5, adiag, v, w
      dimension vecf1(mdi), vecf2(mdi)
      dimension hp5(mdi0), adiag(mdi0), v(mdi0), w(mdi0)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      real*8 av, h, u
      integer mx, mx1
      dimension av(mx1), h(mx1), u(mx,mx)
c
      logical odebug
c
      integer ndav,ncall
      integer idim,i
      common /c/ idim
      common /cdiag_adler/ndav
c ----
c --- kopiere die diagonale
c --- copy the diagonal
      call zro1d(idim,v)
      call zro1d(idim,w)
      do i=1,idim
         adiag(i) = hp5(i)
      end do
c steuerung per hand ob multi-davidson oder normaler davidson
c steuerung ausser kraft, da in readin.f ndav auf 1 gesetzt wird
c ndav = 1    multi
c ndav = 2    engels
c control by hand whether the multi-davidson or normal davidson is used
cbe in readin it is ensured that only david1.f is called
c ndav = 1    multi
c ndav = 2    engels
c
*     ndav = 1
c anzahl aufrufe for extrapolation
c ncall=0 : erster  aufruf
c ncall=1 : zweiter aufruf
c number of calls for extrapolation
c ncall=0 : first call
c ncall=1 : second call
      ncall =1
c     this version of the code only supports ndav = 1 (multi)
c     remove other code for the time being
        call david1a(twoe, nteint,
     +    pey, acoul, aexc, ndeke, core,
     +    vecf1, vecf2, mdi, hp5, mdi0,
     +    emat, nedim,
     +    iottr,iotnew,iotm,
     +    maindf,maxr,jconb,nomax2,
     +    ncall,idim,adiag,v,w,
     +    av, h, u, mx1, mx, odebug)
c
      return
      end
cvp
      subroutine dk0int(twoe,nteint,ston,nidmax,
     +   iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   isk,jsk,ibs,nspiel,mspiel,jnum1
     &  ,icmaxm,nston,nrec31,lend,iend,jfull,ncntp3)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer nteint, nidmax
      real*8 twoe
      dimension twoe(nteint)
      real*8 ston
      dimension ston(nidmax)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c speziell for dk=0
c jmerk laufen direkt vom laeufer angesteuert
c dadurch erfolgt der zugriff schneller
c geaenderte version: nur nummern der wechselwirkenden konf. zu einer
c bei gleichen sk erfolgt der vergleich bzgl. der unteren dreieckmatrix
c (vp 13.3.1992)
c especially for dk=1
c jmerk passes driven directly by the loop counter
c this speeds up the memory access
c modified version: just the numbers of interacting configurations
c in case of the same sk the comparison regards the lower triangle
c (vp 13.3.1992)
c  festen konf. werden bestimmt
c  uebergeben werden: iswh=anzahl der sk
c                     isk=linke sk
c                     jsk=rechte sk
c                     ibs=nummer der konf. aus der sk isk
c  resultate sind:    nspiel=anzahl der ww konf.
c                     ispiel=vektor, der die nummern der ww konf.
c                            enthaelt
c  determining permanent configurations
c  inputs are       : iswh=number of sk
c                     isk=left sk
c                     jsk=right sk
c                     ibs=number of configurations from the sk isk
c  outputs are   :    nspiel=number of interacting configurations
c                     ispiel=vector containing the numbers of the
c                            interacting configurations
cvp iadl: zwischenspeicherung der integraladressen
cvp iadl: intermediate storage of the integral addresses
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
      real*8 sump3
      integer iadl,irec,iadp3
cvp
cvp
cx    ncisk=nconf(isk)
cx      inytl  = nytl(isk)
cx      indub  = ndub(isk)
cx      inopen = inytl - indub
cx      istar=niot(isk)+1
c     do 2032 ibs=1,ncisk
cvp
cx    call clockv(vu1,cpu1)
cvp
cx    iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
cx      do 215 k=1,imo
cx      jcon(k)=-1
cx      jcon1(k)=1
cx      jcon2(k)=2
cx 215  continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
cx       do 205 k=1,inytl
cx205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
cx       do 220 k=1,inopen
cx       jposo(itest(k))=k
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 1
cx       jcon(itest(k)) = 0
cx220    continue
cvp
c        if (ibs.eq.4.or.ibs.eq.42) then
c        write(6,*) ' sk=',isk
c        write(6,*) ' phi_r: ',(iot(k),k=iloc,iloc+inytl-1)
c        endif
cvp
cx       do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
cx       jposc(itest(k)) = k-inopen
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 0
cx       jcon (itest(k)) = 1
cx221    continue
cvp ipos zaehlt alle besetzten mo's durch
c        ipos=0
c        do 222 k=1,imo
c        if (jcon1(k).ne.0) then
c          ipos=ipos+1
c          jposg(k)=ipos
c        endif
c 222    continue
cx       if(jsk.lt.1) go to 250
cx       jnytl=nytl(jsk)
cx       jndub=ndub(jsk)
cx       jnopen=jnytl - jndub
cx       ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
cvp vergleich bzgl. unterer dreiecksmatrix
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cx       nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nullsetzen von idiff feld soweit notwendig
cx       do 255 ialex = 1,nspiel
cvp nullsetzen der adressen der p-faelle
cx 255   idiff(ialex) = 0
c test der offenen schalen
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
cx       jstar=niot(jsk)
cvp jstar wird um die nummer der ersten zu testenden konfiguration
cvp  -1 hochgesetzt
cx       jstar=jstar+jnum1
cvp
cx       kstar=jstar+jnopen*ncjsk
cx       do 360 jdoo = 1,jnopen
cx        kstar=kstar-ncjsk
*vocl loop,novrec
cx        do 380 kdoo = 1,nspiel
cx        idiff(kdoo) = idiff(kdoo) + jcon1( iottr(kstar+kdoo) )
cx380     continue
cx360    continue
c do loop ueber die geschlossenen schalen
c beginn hinten, da dort die meisten unterschiede zu finden sind
cx       kstar=jstar+jnytl*ncjsk
cx       do 260 jdoo = 1,jndub
cx        kstar=kstar-ncjsk
c unterschiede zur testkonfiguration auf idiff summieren
cx        do 270 kdoo = 1,nspiel
cx        if (idiff(kdoo).le.2) then
cx        idiff(kdoo) = idiff(kdoo) + jcon2( iottr(kstar+kdoo) )
cx        endif
cx270     continue
cx260    continue
      ncisk=nconf(isk)
      inytl  = nytl(isk)
      indub  = ndub(isk)
      inopen = inytl - indub
      istar=niot(isk)+1
cvp
c     call clockv(vu1,cpu1)
cvp
      iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
        do 215 k=1,imo
        jcon(k)=-1
        jcon1(k)=1
        jcon2(k)=2
215     continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
         jcon(itest(k)) = 0
  220    continue
c
         do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
c first occupy jcon1, such that starting from a single and finding a
c a double is not counted.
         jposc(itest(k)) = k-inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
         jcon (itest(k)) = 1
  221    continue
c
         if(jsk.lt.1) go to 250
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
c building field ispiel ( nr. of players) eliminate all those
c that have been identified as double excitations already
c
cvp vergleich bzgl. unterer dreiecksmatrix
cvp comparison regarding the lower triangle
c
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cvp  passen alle auf die vektoren, ist nspiel=ibs-1
cvp here nspiel is the max. number of configurations to be tested
cvp  that fit onto the vectors
cvp  if all fit onto the vectors then nspiel=ibs-1
         nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c ispiel : number of playing configurations within the sk
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c  von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using
c creation and annihilation operators
      call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
c
cbs
cbs jetzt wird komprimiert
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
             nspie1 = nspie1 + 1
cx           ispiel(nspie1) = kdoo
             ispiel(nspie1) = kdoo + jnum1
             idiff(nspie1)=idiff(kdoo)
cvp          if(idiff(kdoo).eq.1) then
cvp            nspp3=nspp3+1
cvp            ispp3(nspp3)=kdoo
cvp            npfal(nspie1)=4
cvp          endif
           endif
 385      continue
cvp jstar wird wieder auf alten wert zurueckgesetzt
cvp jstar is reset to the old value
         jstar=niot(jsk)
cvp
cvp  ende bestimmung der wechselwirkenden konfigurationen
cvp  end of interacting configuration determination
cvp
c     call clockv(vu2,cpu2)
cvp
          nspiel=nspie1
cbs   und jetzt das ganze noch mal
cbs   and now do the whole thing once more
cvp   setzen von startwerten
cvp   initialisation

           do 261 jbs=1,nspiel
 261       jmerko(jbs)=0
cvp lauf ueber alle offenen schalen zur p-fall-trennung
cvp loop over all open shells to separate p-cases
cvp
c test der offenen schalen
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
c test the open shells
c like with the closed shells start with the last open shell
         kstar=jstar+jnopen*ncjsk
         do 3050 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c differenz bilden mit jcon
c compute difference with jcon
*vocl loop,novrec
          do 3040 kdoo = 1,nspiel
ct          jmerko(kdoo)=jmerko(kdoo)+abs(jcon( iottr(kstar+
ct   &                ispiel(kdoo)) )-1)
            jmerko(kdoo)=jmerko(kdoo)+abs(jcon( iottr(kstar+
     &                ispiel(kdoo)) ))
 3040     continue
 3050    continue
cvp bestimmung des p-falles durch jmerko und idiff
c notation nach buenker, d.h. npfal=2 --> p=1 usw.
cvp determination of p-case by jmerko and idiff
cvp notation following buenker, i.e. npfal=2 --> p=1, etc.
          nspp1=0
          nspp2=0
          nspp3=0
          nspp4=0
      do 3100 kdoo=1,nspiel
       if (jmerko(kdoo).eq.2) then
         nspp1=nspp1+1
         nwwmo(nspp1,1)=ispiel(kdoo)
cvp      npfal(kdoo)=2
cvp    endif
        else if (jmerko(kdoo).eq.1.and.idiff(kdoo).eq.2) then
         nspp2=nspp2+1
         nwwmo(nspp2,2)=ispiel(kdoo)
cvp      npfal(kdoo)=3
cvp    endif
        else if (idiff(kdoo).eq.1) then
         nspp3=nspp3+1
         nwwmo(nspp3,3)=ispiel(kdoo)
cvp      npfal(kdoo)=4
cvp    endif
cvp    if (jmerko(kdoo).eq.0) then
        else
         nspp4=nspp4+1
         nwwmo(nspp4,4)=ispiel(kdoo)
cvp      npfal(kdoo)=5
       endif
 3100 continue
cvp konfiguratonsnummern nach p-faellen sortiert --> npfal
cvp configuration numbers sorted according to p-cases --> npfal
      do 3210 kdoo=1,nspp1
        npfal(kdoo)=nwwmo(kdoo,1)
        jmerko(kdoo)=2
        jmerkd(kdoo)=0
        nqrfal(kdoo)=ibinom(inopen,2)
 3210 continue
      do 3220 kdoo=1,nspp2
        npfal(kdoo+nspp1)=nwwmo(kdoo,2)
        jmerkd(kdoo+nspp1)=0
        nqrfal(kdoo+nspp1)=ideks(inopen+1)
        nwwmo(kdoo+nspp1,1)=ideks(indub+1)
 3220 continue
      nsppa=nspp1+nspp2
      do 3230 kdoo=1,nspp3
        npfal(kdoo+nsppa)=nwwmo(kdoo,3)
        jmerko(kdoo+nsppa)=0
        jmerkd(kdoo+nsppa)=0
        nqrfal(kdoo+nsppa)=ideks(inopen+1)
 3230 continue
      nsppb=nsppa+nspp3
      do 3240 kdoo=1,nspp4
        npfal(kdoo+nsppb)=nwwmo(kdoo,4)
        moafal(kdoo+nsppb)=ideks(indub+1)
 3240 continue
cvp schleife ueber alle doppelt besetzten mo's von phi_l
c        write(6,*) 'nspiel = ',nspiel
c beginn hinten, da dort die meisten unterschiede zu finden sind
cvp loop over all doubly occupied mo's of \phi_l
c start at end because most differences will be found there
         kstar=jstar+jnytl*ncjsk
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
cvp schleife for p=1
cvp loop for p=1
cvp jmerkd zaehlt die besetzungsunterschiede bzgl. der doppelt
c   besetzten mo's von \phi_l
cvp jmerkd counts the occupation differences with respect to the
cvp doubly occupied mo's of \phi_l
*vocl loop,novrec
          do 281 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              npos(kdoo,1)=jposo(iottr(kstar+npfal(kdoo)))
            endif
 281      continue
cvp schleife for p=2
cvp loop for p=2
*vocl loop,novrec
          do 282 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              nwwmo(kdoo,2)=iottr(kstar+npfal(kdoo))
cvp jmerkd(*)=1: r=1
ct            if (jdiff(kdoo).eq.-2) jmerkd(kdoo)=1
              if (jdiff(kdoo).eq.2) jmerkd(kdoo)=1
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              nwwmo(kdoo,1)=nwwmo(kdoo,1)-nx
            endif
 282      continue
cvp schleife for p=3
cvp loop for p=3
*vocl loop,novrec
          do 283 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
            endif
 283      continue
cvp schleife for p=4
cvp loop for p=4
cvp jmerkd zaehlt die besetzungsunterschiede bzgl. der doppelt
c   besetzten mo's von \phi_l
cvp jmerkd counts the occupation differences with respect to the
cvp doubly occupied mo's of \phi_l
*vocl loop,novrec
          do 284 kdoo=nsppb+1,nsppb+nspp4
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              mobfal(kdoo)=iottr(kstar+npfal(kdoo))
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              moafal(kdoo)=moafal(kdoo)-nx
          endif
 284      continue
1260     continue
cvp
cvp lauf ueber die offenen schalen von hinten
cvp loop over the open shells starting at the end
cvp
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c         write(6,*) 'open shell nr. ',jdoo,
c    *'jstar  ',jstar
c differenz bilden mit jcon
c compute difference with jcon
cvp jmerko zaehlt die besetzungsunterschiede bzgl. der einfach
c   besetzten mo's von \phi_l, wird hier wieder auf null herunter-
c   gezaehlt
cvp jmerko counts the occupation differences with respect to the
cvp singly occupied mo's of \phi_l, is counted down to zero again here
cvp p=1-fall
*vocl loop,novrec
          do 481 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            jmerko(kdoo)=jmerko(kdoo)-1
            npos(kdoo,2)=jdiff(kdoo)
            nwwmo(kdoo,4-jmerko(kdoo))=inopen+1-jdoo
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen+1-jdoo-jmerko(kdoo))
          endif
481       continue
cvp p=2-fall
*vocl loop,novrec
          do 482 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            nqlfal(kdoo)=inopen-jdoo+1
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-nx
          endif
482       continue
cvp p=3-fall
*vocl loop,novrec
          do 483 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            nqlfal(kdoo)=inopen-jdoo+1
            mobfal(kdoo)=iottr(kstar+npfal(kdoo))
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-nx
          endif
483       continue
 1360   continue
cvp bestimmung des r-falles for p=1
cvp determination of the r-case for p=1
      do 491 kdoo=1,nspp1
          if (jmerkd(kdoo).eq.1) then
cvp        dopppelt besetztes mo von phi_l trifft auf das hoehere
cvp         einfach besetzte von phi_r
cvp        doubly occupied mo of \phi_l meets the higher singly 
cvp         occupied mo of \phi_r
            if (npos(kdoo,1).eq.imoq(nqrfal(kdoo),1)) then
                npos(kdoo,1)=-1
              else
                npos(kdoo,1)=1
            endif
cvp          npos(*,1) gibt jetzt relative besetzung des hinteren
cvp           offenen mo's von phi_r an
cvp          npos(*,1) now stores the relative occupation of the last
cvp           open mo's of \phi_r
            if (npos(kdoo,1).eq.npos(kdoo,2)) then
              nrfal(kdoo)=2
             else
              nrfal(kdoo)=3
            endif
           else
            nrfal(kdoo)=1
          endif
  491 continue
cvp
cvp  ende labelbestimmung
cvp  end of label determination
cvp
c     call clockv(vu3,cpu3)
cvp
cvp bestimmung der ww mo's von \phi_r aus qr,
cvp berechnung und sortierung der absoluten mo-nummern for
cvp die integralberechnung
cvp und bestimmung des ql-falles aus den ww mo's von \phi_l
cvp determination of the interacting mo's of \phi_r in qr,
cvp calculation and sorting of the absolute mo-numbers for the
cvp integral calcalution
cvp and the determination of the ql-cases from the interacting mo's
cvp of \phi_l
         kstar = niot(jsk) - ncjsk
*vocl loop,novrec
      do 591 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp sei das exchange-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff und jmerko sind hier frei verfuegbar
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp if given a exchange-integral in chemical-notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff and jmerko are freely available
        nqlfal(kdoo)=iqmo(nwwmo(kdoo,3),nwwmo(kdoo,4))
          nwwmo(kdoo,3)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,3)*ncjsk)
          nwwmo(kdoo,4)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,4)*ncjsk)
          nwwmo(kdoo,1)=itest(imoq(nqrfal(kdoo),1))
          nwwmo(kdoo,2)=itest(imoq(nqrfal(kdoo),2))
c       write(6,*) ' r-case, ww mo: ',nrfal(kdoo),(nwwmo(kdoo,k),k=1,4)
          if (nwwmo(kdoo,1).lt.nwwmo(kdoo,3)) then
            nx1=nwwmo(kdoo,1)
            ny1=nwwmo(kdoo,2)
            nwwmo(kdoo,1)=nwwmo(kdoo,3)
            nwwmo(kdoo,2)=nwwmo(kdoo,4)
            nwwmo(kdoo,3)=nx1
            nwwmo(kdoo,4)=ny1
          endif
cvp sortierung for ex-integrale
cvp sort of ex-integrals
          if (nrfal(kdoo).eq.2) then
            idiff(kdoo)=nwwmo(kdoo,3)
            if (nwwmo(kdoo,2).gt.nwwmo(kdoo,4)) then
              jdiff(kdoo)=nwwmo(kdoo,2)
              jmerko(kdoo)=nwwmo(kdoo,4)
             else
              jdiff(kdoo)=nwwmo(kdoo,4)
              jmerko(kdoo)=nwwmo(kdoo,2)
            endif
           else
            idiff(kdoo)=nwwmo(kdoo,4)
            if (nwwmo(kdoo,2).gt.nwwmo(kdoo,3)) then
              jdiff(kdoo)=nwwmo(kdoo,2)
              jmerko(kdoo)=nwwmo(kdoo,3)
             else
              jdiff(kdoo)=nwwmo(kdoo,3)
              jmerko(kdoo)=nwwmo(kdoo,2)
            endif
          endif
cvp sortierung for cb-integrale
cvp sort of cb-integrals
          if (jmerkd(kdoo).ne.1) then
            nxx=nwwmo(kdoo,2)
            nwwmo(kdoo,2)=nwwmo(kdoo,3)
            nwwmo(kdoo,3)=nxx
            if (nwwmo(kdoo,3).lt.nwwmo(kdoo,4)) then
              nyy=nwwmo(kdoo,4)
              nwwmo(kdoo,4)=nwwmo(kdoo,3)
              nwwmo(kdoo,3)=nyy
            endif
          endif
c       write(6,*) ' r-case, ww mo: ',nrfal(kdoo),(nwwmo(kdoo,k),k=1,4)
c     write(6,*)
  591 continue
cvp
cvp p=2 fall
cvp berechnung der absolutnummern der mo's und sortierung
cvp calculation of the absolute mo-numbers and sorting
*vocl loop,novrec
      do 492 kdoo=nspp1+1,nsppa
        nwwmo(kdoo,1)=itest(inopen+nwwmo(kdoo,1))
        if (jmerkd(kdoo).eq.1) then
cvp r=1
          nwwmo(kdoo,3)=itest(nqrfal(kdoo))
          nx3=nwwmo(kdoo,1)
          nwwmo(kdoo,1)=nwwmo(kdoo,2)
          nwwmo(kdoo,2)=nx3
         else
cvp r=2
c         nwwmo(kdoo,3)=iot(jstar+npfal(kdoo)*jnytl-1+nqlfal(kdoo))
          nwwmo(kdoo,3)=iottr(kstar+npfal(kdoo)+nqlfal(kdoo)*ncjsk)
        endif
cvp mo's sind nun wie folgt sortiert:
cvp  nwwmo(*,1) = doppelt besetztes mo, das in der anderen konf. leer
cvp               ist
cvp  nwwmo(*,2) = doppelt besetztes mo, das in der anderen konf. auf
cvp               ein einfach besetztes trifft
cvp  nwwmo(*,3) = einfach besetztes mo, das in der anderen konf. leer
cvp               ist
cvp mo's are sorted as follows:
cvp  nwwmo(*,1) = doubly occupied mo that is empty in other 
cvp               configuration
cvp  nwwmo(*,2) = doubly occupied mo that meets a singly occupied one
cvp               in the other configuration
cvp  nwwmo(*,3) = singly occupied mo that is empty in other
cvp               configuration
  492 continue
cvp
c     do 3500 kdoo=nspp1+1,nsppa
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,3)
c3500 continue
cvp
*vocl loop,novrec
      do 592 kdoo=nspp1+1,nsppa
        if (nwwmo(kdoo,2).lt.nwwmo(kdoo,3)) then
          nx4=nwwmo(kdoo,2)
          nwwmo(kdoo,2)=nwwmo(kdoo,3)
          nwwmo(kdoo,3)=nx4
        endif
  592 continue
*vocl loop,novrec
      do 493 kdoo=nsppa+1,nsppb
cvp bestimmung der irred. darst. und abspeichern der hoeheren und
c   der niedrigeren mo-nummer innerhalb der irred. darst.
cvp determination of the irrep and storing the higher and lower mo
cvp number within the irrep
cvp
cvp idiff: position des offenen mo's von \phi_r --> wird for
cvp cb- und ex-integrale benoetigt
cvp idiff: position of the open mo's of \phi_r --> needed for cb- and
cvp ex-integrals
c       idiff(kdoo)=moafal(kdoo)
        idiff(kdoo)=nqrfal(kdoo)
cvp
        nrfal(kdoo)=nirred(mobfal(kdoo))
        mobfal(kdoo)=mopos(mobfal(kdoo))
c       moafal(kdoo)=mopos(iot(iloc-1+moafal(kdoo)))
        moafal(kdoo)=mopos(itest(nqrfal(kdoo)))
cvp moafal enthaelt groesseres mo
cvp moafal holds the higher mo
        if (moafal(kdoo).lt.mobfal(kdoo)) then
          nx=moafal(kdoo)
          moafal(kdoo)=mobfal(kdoo)
          mobfal(kdoo)=nx
        endif
cvp codierung des r-falles for p=3 durch vorzeichen von qr
cvp the r-case for p=3 is indicated by the sign of qr
        if (jmerkd(kdoo).eq.1) nqrfal(kdoo)=-1*nqrfal(kdoo)
 493  continue
cvp
cvp
cvp bestimmung der absoluten mo-nummern for p=4
cvp determination of the absolute mo-numbers for p=4
      do 494 kdoo=nsppb+1,nsppb+nspp4
          moafal(kdoo)=itest(inopen+moafal(kdoo))
c         mobfal(kdoo)=iot(iloc+inopen-1+mobfal(kdoo))
 494  continue
cvp
cvp
cvp  berechnung der integraladresen
cvp  calculation of the integral addresses
cvp
c     iloc1=iloc-1
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3201 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp zugriff auf nit-feld
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp access nit field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
cvp bestimmung des symmetrie-falles (cb-integral)
cvp determining the symmetry case (cb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum2) then
c       if (ndum2.eq.ndum3) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum3) then
cvp      fall 1
c           write(6,*) ' case 1'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks(ideks(idum1)+
     &        mopos(idiff(kdoo)))
     &        +ideks(mopos(jdiff(kdoo)))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
c          else if (ndum1.eq.ndum2) then
cvp      fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))
     &       +ideks(ndum3+1) ) + idum4 +
     &       ideks(ncimo(ndum3)+1)*(ideks(idum1)+idum2-1)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &       +mopos(idiff(kdoo)) )
     &       +ncimo(nirred(idiff(kdoo)))
     &       *(mopos(jdiff(kdoo))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum3) then
c         else if (ndum1.eq.ndum3) then
cvp      fall 3
c           write(6,*) ' case 3'
            intcb(kdoo)=nit( ideks(idum+ndum2+1) )
     &       +idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
            if (nrfal(kdoo).ne.2) then
              intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &         +mopos(idiff(kdoo)) )
     &         +ncimo(nirred(idiff(kdoo)))
     &         *(mopos(jdiff(kdoo))-1)
             else
              intex(kdoo)=nit( ideks(ideks(ndum1+1))
     &         +ideks(nirred(jdiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks(ncimo(nirred(jdiff(kdoo)))+1)
     &         *(ideks(idum1)+mopos(idiff(kdoo))-1)
     &         +ideks(mopos(jdiff(kdoo)))
            endif
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 4
c           write(6,*) ' case 4'
cvp nur i-1 taucht in den formeln auf
cvp only i-1 shows up in the equations
            idum1=idum1-1
            intcb(kdoo)=nit( ideks(idum+ndum2)
     &       +ideks(ndum3)+ndum4 ) + idum4 +
     &       ncimo(ndum4)*(idum3-1+ncimo(ndum3)*(idum2-1+
     &       ncimo(ndum2)*idum1))
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo)))
     &       +ideks(nirred(jdiff(kdoo)))+nirred(jmerko(kdoo)) )
     &       + mopos(jmerko(kdoo)) +
     &       ncimo(nirred(jmerko(kdoo)))
     &       *(mopos(jdiff(kdoo))-1
     &       +ncimo(nirred(jdiff(kdoo)))
     &       *(mopos(idiff(kdoo))-1
     &       +ncimo(nirred(idiff(kdoo)))*idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
 3201 continue
cvp
cvp bestimmung der gesamtzahl der integrale
cvp  evtl. umspeicherung auf twoe
cvp nintg = schaetzung for die anzahl aller integrale
cvp determination of the total number of integrals
cvp  if needed restoring them onto twoe
cvp nintg = estimate for the total number of integrals
      nintg=2*nspp1+nspp2
     &     +nspp3*(3*indub+2*inopen-2)+nspp3
      if (iend+nintg.gt.icmaxm) then
         call isortb(twoe,nteint,ston,nidmax,
     +        nston,nrec31,indub,inopen,lend,iend
     &        ,ncntp3)
         ncntp3=0
         do 85 kdoo=1,icmax
   85    iadp3(kdoo)=0
         jfull=ibs-1
         iend=0
      endif
cvp np3r2 = anzahl der p=3,r=2 faelle
cvp np3r2 = number of p=3,r=2 cases
      np3r2=0
*vocl loop,novrec
      do 89 kdoo=1,nspp3
cvp ncntp3 = anzahl der p=3-faelle auf iadl
cvp ncntp3 = number of p=3-ases on iadl
         ncntp3=ncntp3+1
         iadp3(ncntp3)=iend+1+2*nspp1+nspp2
     &        +(kdoo-1)*(3*indub+2*inopen-2) +np3r2
         if (nqrfal(nsppa+kdoo).lt.0) then
            np3r2=np3r2+1
            iadp3(ncntp3)=-iadp3(ncntp3)
         endif
   89 continue
cvp nintp3 = anzahl aller summanden for p=3
cvp nintp3 = number of all terms for p=3
      nintp3=nspp3*(3*indub+2*inopen-2)+np3r2
cvp nintg = anzahl aller integrale
cvp nintg = number of all integrals
      nintg=2*nspp1+nspp2
     &     +nintp3
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
*vocl loop,novrec
      do 120 kdoo=1,nspp1
         iadl(iend+2*kdoo-1)=intcb(kdoo)
         iadl(iend+2*kdoo  )=intex(kdoo)
  120 continue
      iend=iend+2*nspp1
cvp berechnung der integraladresse for p=2
cvp calculation of the integral addresses for p=2
*vocl loop,novrec
      do 3202 kdoo=nspp1+1,nsppa
cvp das integral sei (ab|ac) oder (ac|ab) usw.
cvp dann ist nwwmo(*,1)=a
cvp          nwwmo(*,2)=max(b,c)
cvp          nwwmo(*,3)=min(b,c)
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the integral is (ab|ac) or (ac|ab) etc.
cvp then     nwwmo(*,1)=a
cvp          nwwmo(*,2)=max(b,c)
cvp          nwwmo(*,3)=min(b,c)
cvp however it should still be sorted according to size
cvp access to nit field
        ndum1=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        ndum2=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,2)))
c       ndum3=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
c       ndum4=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        idum2=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,2)))
        idum3=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum4=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
c       idum=ideks(ndum1)
c       intcb(kdoo)=nit( ideks(idum+ndum2)
c    &    +ideks(ndum3)+ndum4 ) + idum4
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
cvp        fall 1
c           write(6,*) ' case 1'
        intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &     + idum4 +
c           intcb(kdoo)=intcb(kdoo)+
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
          else
cvp        fall 3
c           write(6,*) ' case 3'
        intcb(kdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &     + idum4 +
c           intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
        endif
 3202 continue
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
      do 121 kdoo=1,nspp2
         iadl(iend+kdoo)=intcb(nspp1+kdoo)
  121 continue
      iend=iend+nspp2
cvp
cvp  p=3
cvp
cvp  berechnung der summe uber cb- und ex-integrale sowie
cvp  der adressen der ex-integrale bei gemeinsamen offenen
cvp  schalen
cvp  calculate the sum over cb- and ex-integrals and also
cvp  the addresses of the ex-integrals in case of common open shells
cvp
cvp  sei d gemeinsames offenes mo: betrachtung der integrale
cvp     (ab|dd) und (ad|bd), es gilt a > b
cvp  if d is a common open mo: considering the integrals
cvp     (ab|dd) and (ad|bd), where a > b
cvp
      do 593 jdoo=1,inopen-1
        do 603 kdoo=nsppa+1,nsppb
          mogd=jdoo+jmerko(kdoo)
          if (mogd.eq.idiff(kdoo)) jmerko(kdoo)=1
  603   continue
*vocl loop,novrec(iadl,icb)
        do 613 kdoo=nsppa+1,nsppb
          mogd=itest(jdoo+jmerko(kdoo))
cvp das ex-integral sei (ad|bd), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ad|bd), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp it still has to be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogd)
cvp bestimmung des symmetrie-falles
cvp determining the symmetrie case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogd))
           idum2=min(moafal(kdoo),mopos(mogd))
           idum3=max(mobfal(kdoo),mopos(mogd))
           idum4=min(mobfal(kdoo),mopos(mogd))
cvp        fall 1
c           write(6,*) ' case 1'
            nwwmo(kdoo,jdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogd)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            intcb(kdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
cvp
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogd)
           idum3=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
cvp
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogd)
           idum4=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
cvp
        endif
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
        istp3=iabs(iadp3(ncntp3+kdoo-nsppa-nspp3))
     &             +3*indub-1
        iadl(istp3+jdoo)=intcb(kdoo)
  613   continue
  593 continue
cvp schleife ueber die (gemeinsamen) geschlossenen schalen
cvp loop over the (common) closed shells
      do 533 jdoo=1,indub
*vocl loop,novrec
        do 623 kdoo=nsppa+1,nsppb
          mogc=itest(inopen+jdoo)
cvp das ex-integral sei (ac|bc), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ac|bc), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp it should still be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogc)
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogc))
           idum2=min(moafal(kdoo),mopos(mogc))
           idum3=max(mobfal(kdoo),mopos(mogc))
           idum4=min(mobfal(kdoo),mopos(mogc))
cvp        fall 1
c           write(6,*) ' case 1'
            intex(kdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogc)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            intcb(kdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
cvp
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogc)
           idum3=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            intex(kdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogc)
           idum4=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            intex(kdoo)=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
        endif
cvp abspeichern der summanden auf adl
cvp store terms on adl
        istp3=iabs(iadp3(ncntp3+kdoo-nsppa-nspp3))
        iadl(istp3+jdoo-1)=intex(kdoo)
        iadl(istp3+indub+2*jdoo-2)=intcb(kdoo)
        iadl(istp3+indub+2*jdoo-1)=intcb(kdoo)
  623   continue
  533 continue
cvp korrektur der integralsumme for r=2
cvp correct the integral sum for r=2
*vocl loop,novrec
        do 633 kdoo=nsppa+1,nsppb
          if (nqrfal(kdoo).lt.0) then
cvp zugriff auf nit-feld
cvp access nit-field
            ndum1=nrfal(kdoo)
            idum1=itest(idiff(kdoo))
            ix=ideks(mopos(idum1)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            intcb(kdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
cvp abspeichern des einzelnen summanden for r=2
cvp store the separate terms for r=2
            istp3=iabs(iadp3(ncntp3+kdoo-nsppa-nspp3))
            iadl(istp3+3*indub+inopen-1)=intcb(kdoo)
          endif
  633   continue
cvp
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
      do 123 jdoo=1,inopen-1
*vocl loop,novrec
        do 122 kdoo=nsppa+1,nsppb
        if (nqrfal(kdoo).lt.0) then
           idumr2=1
          else
           idumr2=0
        endif
        istp3=iabs(iadp3(ncntp3+kdoo-nsppa-nspp3))
     &        +3*indub+inopen-2
        iadl(istp3+idumr2+jdoo)=
     &           nwwmo(kdoo,jdoo)
  122   continue
  123 continue
      iend=iend+nintp3
cvp
cvp  ende berechnung der integraladressen
cvp  end of integral address calculation
cvp
c     call clockv(vu4,cpu4)
cvp
c     vutww=vu2-vu1
c     cpuww=cpu2-cpu1
c     vutla=vu3-vu2
c     cpula=cpu3-cpu2
c     vutin=vu4-vu3
c     cpuin=cpu4-cpu3
cvp
cvp timer-routine zum test des konfigurationsvergleich
*     call clockv(vu20,cpu20)
*     call time(mill20)
cvp
c     if (ibs.eq.4.or.ibs.eq.42) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c       write(6,*) ' intcb= ',intcb(kdoo)
c3500 continue
c     endif
cvp
cvp ruecksortierung wird uebergangen !!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel enthaelt jetzt die adressen der wechselwirkenden
cvp  konfigurationen: zuerst die for p=1, dann for p=2 usw.
cvp  auf npfal wird vorlaeufig der p-fall abgespeichert
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp back-sort is skipped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel now holds the addresses of the interacting
cvp  configurations: first those for p=1, then for p=2 usw.
cvp  provisionaly the p-case is stored on npfal
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 5210 kdoo=1,nspp1
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=2
 5210 continue
      do 5220 kdoo=nspp1+1,nsppa
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=3
 5220 continue
      do 5230 kdoo=nsppa+1,nsppb
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=4
 5230 continue
      do 5240 kdoo=nsppb+1,nsppb+nspp4
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=5
 5240 continue
cvp
cvp p-faelle werden in urspruengliche reihenfolge gebracht
cvp p-cases being put in original order
cvp
cre    lspp1=1
cre    lspp2=1
cre    lspp3=1
cre    lspp4=1
cvp
c      write(6,*) ' jbs=,nspp3=,ispp3= ',jbs,nspp3,(ispp3(kkkk),
c    &   kkkk=1,nspp3)
cvp
cre    do 5100 kdoo=1,nspiel
cre      if (ispiel(kdoo).eq.npfal(lspp1).and.lspp1.le.nspp1) then
cre        ndfal1(kdoo)=2
cre        ndfal2(kdoo)=nrfal(lspp1)
cre        ndfal3(kdoo)=nqlfal(lspp1)
cre        ndfal4(kdoo)=nqrfal(lspp1)
cre        ndfal7(kdoo)=intcb(lspp1)
cre        ndfal8(kdoo)=intex(lspp1)
cre        lspp1=lspp1+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp2+nspp1).and.lspp2.le.nspp2)
cre  &  then
cre        ndfal1(kdoo)=3
cre        ndfal3(kdoo)=nqlfal(lspp2+nspp1)
cre        ndfal4(kdoo)=nqrfal(lspp2+nspp1)
cre        ndfal7(kdoo)=intcb(lspp2+nspp1)
cre        lspp2=lspp2+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp3+nsppa).and.lspp3.le.nspp3)
cre  &  then
c          write(6,*) ' kdoo=,ispiel=,ispp3= ',kdoo,ispiel(kdoo)
c    &       ,ispp3(lspp3)
cre        ndfal1(kdoo)=4
cre        ndfal2(kdoo)=nrfal(lspp3+nsppa)
cre        ndfal3(kdoo)=nqlfal(lspp3+nsppa)
cre        ndfal4(kdoo)=nqrfal(lspp3+nsppa)
cre        ndfal5(kdoo)=moafal(lspp3+nsppa)
cre        ndfal6(kdoo)=mobfal(lspp3+nsppa)
cre        ndfeld(kdoo,1)=nwwmo(lspp3+nsppa,1)
cre        ndfeld(kdoo,2)=nwwmo(lspp3+nsppa,2)
cre        ndfeld(kdoo,3)=nwwmo(lspp3+nsppa,3)
cre        ndfeld(kdoo,4)=nwwmo(lspp3+nsppa,4)
cre        ndfeld(kdoo,5)=nwwmo(lspp3+nsppa,5)
cre        ndfeld(kdoo,6)=nwwmo(lspp3+nsppa,6)
cre        ndfeld(kdoo,7)=nwwmo(lspp3+nsppa,7)
cre        ndfeld(kdoo,8)=nwwmo(lspp3+nsppa,8)
cre        ndfeld(kdoo,9)=nwwmo(lspp3+nsppa,9)
cre        rdum(kdoo)=sumint(lspp3+nsppa)
cre        lspp3=lspp3+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp4+nsppb).and.lspp4.le.nspp4)
cre  &  then
cre        ndfal1(kdoo)=5
cre        ndfal5(kdoo)=moafal(lspp4+nsppb)
cre        ndfal6(kdoo)=mobfal(lspp4+nsppb)
cre        lspp4=lspp4+1
cre      endif
c5100  continue
 250    continue
      return
      end
c
      subroutine dk1int(twoe,nteint,ston,nidmax,
     +   iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   isk,jsk,ibs,nspiel,mspiel,jnum1
     &  ,icmaxm,nston,nrec31,lend,iend,jfull,ncntp3)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer nteint, nidmax
      real*8 twoe
      dimension twoe(nteint)
      real*8 ston
      dimension ston(nidmax)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c speziell for dk=1
c jmerk laufen direkt vom laeufer angesteuert
c dadurch erfolgt der zugriff schneller
c geaenderte version: nur nummern der wechselwirkenden konf. zu einer
c bei gleichen sk erfolgt der vergleich bzgl. der unteren dreieckmatrix
c (vp 13.3.1992)
c especially for dk=1
c jmerk passes driven directly by the loop counter
c this speeds up the memory access
c modified version: just the numbers of interacting configurations
c in case of the same sk the comparison regards the lower triangle
c (vp 13.3.1992)
c  festen konf. werden bestimmt
c  uebergeben werden: iswh=anzahl der sk
c                     isk=linke sk
c                     jsk=rechte sk
c                     ibs=nummer der konf. aus der sk isk
c  resultate sind:    nspiel=anzahl der ww konf.
c                     ispiel=vektor, der die nummern der ww konf.
c                            enthaelt
c  determining permanent configurations
c  inputs are       : iswh=number of sk
c                     isk=left sk
c                     jsk=right sk
c                     ibs=number of configurations from the sk isk
c  outputs are   :    nspiel=number of interacting configurations
c                     ispiel=vector containing the numbers of the
c                            interacting configurations
c
cvp iadl: zwischenspeicherung der integraladressen
cvp iadl: intermediate storage of integral addresses
      real*8 sump3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax),sump3(icmax)
      integer iadl,irec,iadp3
cvp
cvp
cx    ncisk=nconf(isk)
cx      inytl  = nytl(isk)
cx      indub  = ndub(isk)
cx      inopen = inytl - indub
cx      istar=niot(isk)+1
cx    iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
cx      do 215 k=1,imo
cx      jcon(k)=1
cx      jcon1(k)=1
cx      jcon2(k)=2
cx215   continue
cvp abspeichern der testkonfiguration auf itest
cx       do 205 k=1,inytl
cx205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
cx       do 220 k=1,inopen
cx       jposo(itest(k))=k
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 1
cx       jcon(itest(k)) = 0
cx220    continue
cvp
cx       do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
cx       jposc(itest(k)) = k-inopen
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 0
cx221    continue
cvp
cx       if(jsk.lt.1) go to 250
cvp
cx       jnytl=nytl(jsk)
cx       jndub=ndub(jsk)
cx       jnopen=jnytl - jndub
cx       ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
cvp vergleich bzgl. unterer dreiecksmatrix
cx       nspiel=ncjsk
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cx       nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nullsetzen von idiff feld soweit notwendig
cx       do 255 ialex = 1,nspiel
cvp nullsetzen der adressen der p-faelle
cx255    idiff(ialex) = 0
c do loop ueber die offenen schalen
c beginn hinten, da dort die meisten unterschiede zu finden sind
cx       jstar=niot(jsk)
cvp jstar wird um die nummer der ersten zu testenden konfiguration
cvp  -1 hochgesetzt
cx       jstar=jstar+jnum1
cvp
cx       kstar=jstar+jnopen*ncjsk
cx       do 360 jdoo = 1,jnopen
cx        kstar=kstar-ncjsk
*vocl loop,novrec
cx        do 380 kdoo = 1,nspiel
cx        idiff(kdoo) = idiff(kdoo) + jcon1( iottr(kstar+kdoo) )
cx380     continue
cx360    continue
c do loop ueber die geschlossenen schalen
c beginn hinten, da dort die meisten unterschiede zu finden sind
cx       kstar=jstar+jnytl*ncjsk
cx       do 260 jdoo = 1,jndub
cx        kstar=kstar-ncjsk
c unterschiede zur testkonfiguration auf idiff summieren
cx        do 270 kdoo = 1,nspiel
cx        if (idiff(kdoo).le.2) then
cx        idiff(kdoo) = idiff(kdoo) + jcon2( iottr(kstar+kdoo) )
cx        endif
cx270     continue
cx260    continue
      ncisk=nconf(isk)
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
        istar=niot(isk)+1
      iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
        do 215 k=1,imo
        jcon(k)=1
        jcon2(k)=2
        jcon1(k)=1
215     continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
         jcon(itest(k)) = 0
  220    continue
c
         do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
c first occupy jcon1, such that starting from a single and finding a
c a double is not counted.
         jposc(itest(k)) = k-inopen
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
cvp
         if(jsk.lt.1) go to 250
cvp
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
c building field ispiel ( nr. of players) eliminate all those
c that have been identified as double excitations already
c
cvp vergleich bzgl. unterer dreiecksmatrix
cvp comparison regarding the lower triangle
cx       nspiel=ncjsk
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cvp here nspiel is the max. number of configurations to be tested
cvp  that fit onto the vectors
         nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c ispiel : number of playing configurations within the sk
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c nullsetzen von idiff feld soweit notwendig
c nullify idiff field where needed
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c  von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using
c creation and annihilation operators
      call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
cbs
cbs jetzt wird komprimiert
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
             nspie1 = nspie1 + 1
cx           ispiel(nspie1) = kdoo
             ispiel(nspie1) = kdoo + jnum1
             idiff(nspie1)=idiff(kdoo)
cvp          if(idiff(kdoo).eq.1) then
cvp            nspp3=nspp3+1
cvp            ispp3(nspp3)=kdoo
cvp            npfal(nspie1)=4
cvp          endif
           endif
 385      continue
cvp jstar wird wieder auf alten wert zurueckgesetzt
cvp jstar is reset to the old value
         jstar=niot(jsk)
cvp
cvp  ende bestimmung der wechselwirkenden konfigurationen
cvp  end of interacting configuration determination
cvp
cx    call clockv(vu2,cpu2)
cvp
          nspiel=nspie1
cbs   und jetzt das ganze noch mal
cbs   and now do the whole thing once more
cvp   setzen von startwerten
cvp   initialisation
           do 261 jbs=1,nspiel
 261       jmerko(jbs)=0
cvp lauf ueber alle offenen schalen zur p-fall-trennung
cvp loop over all open shells to separate p-cases
cvp
c test der offenen schalen
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
c test the open shells
c like with the closed shells start with the last open shell
         kstar=jstar+jnopen*ncjsk
         do 3050 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c differenz bilden mit jcon
c compute difference with jcon
cvp jmerko zaehlt die besetzungsunterschiede bzgl. der einfach
c   besetzten mo's von \phi_l
cvp jmerko counts the occupation differences with respect to the
cvp singly occupied mo's of \phi_l
*vocl loop,novrec
          do 3040 kdoo = 1,nspiel
ct          jmerko(kdoo)=jmerko(kdoo)+abs(jcon( iottr(kstar+
ct   &                ispiel(kdoo)) )-1)
            jmerko(kdoo)=jmerko(kdoo)+jcon( iottr(kstar+
     &                ispiel(kdoo)) )
 3040     continue
 3050    continue
cvp bestimmung des p-falles durch jmerko und idiff
c notation nach buenker, d.h. npfal=2 --> p=1 usw.
cvp determination of p-case by jmerko and idiff
cvp notation following buenker, i.e. npfal=2 --> p=1, etc.
          nspp1=0
          nspp2=0
          nspp3=0
      do 3100 kdoo=1,nspiel
       if (jmerko(kdoo).eq.1) then
         nspp1=nspp1+1
         nwwmo(nspp1,1)=ispiel(kdoo)
cvp      npfal(kdoo)=2
cvp    endif
        else if (jmerko(kdoo).eq.0.and.idiff(kdoo).eq.2) then
         nspp2=nspp2+1
         nwwmo(nspp2,2)=ispiel(kdoo)
cvp      npfal(kdoo)=3
cvp    endif
c       else if (idiff(kdoo).eq.0) then
        else
         nspp3=nspp3+1
         nwwmo(nspp3,3)=ispiel(kdoo)
cvp      npfal(kdoo)=4
cvp    endif
cvp    if (jmerko(kdoo).eq.0) then
       endif
 3100 continue
cvp konfiguratonsnummern nach p-faellen sortiert --> npfal
cvp configuration numbers sorted according to p-cases --> npfal
      do 3210 kdoo=1,nspp1
        npfal(kdoo)=nwwmo(kdoo,1)
        jmerko(kdoo)=1
        jmerkd(kdoo)=0
        nqrfal(kdoo)=ibinom(inopen,3)
 3210 continue
      do 3220 kdoo=1,nspp2
        npfal(kdoo+nspp1)=nwwmo(kdoo,2)
        jmerkd(kdoo+nspp1)=0
        nqrfal(kdoo+nspp1)=ideks(inopen)
        moafal(kdoo+nspp1)=ideks(indub+1)
 3220 continue
      nsppa=nspp1+nspp2
      do 3230 kdoo=1,nspp3
        npfal(kdoo+nsppa)=nwwmo(kdoo,3)
        jmerko(kdoo+nsppa)=0
        jmerkd(kdoo+nsppa)=0
        nqrfal(kdoo+nsppa)=ideks(inopen)
 3230 continue
      nsppb=nsppa+nspp3
cvp
cvp schleife ueber alle doppelt besetzten mo's von phi_l
cvp loop over all doubly occupied mo's of \phi_l
         kstar=jstar+jnytl*ncjsk
c beginn hinten, da dort die meisten unterschiede zu finden sind
c start at end because most differences will be found there
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
cvp schleife for p=1
cvp loop for p=1
cvp jmerkd zaehlt die besetzungsunterschiede bzgl. der doppelt
c   besetzten mo's von \phi_l
cvp jmerkd counts the occupation differences with respect to the
cvp doubly occupied mo's of \phi_l
*vocl loop,novrec
          do 281 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              npos(kdoo,jmerkd(kdoo))=jposo(iottr(kstar+npfal(kdoo)))
            endif
 281      continue
cvp schleife for p=2
cvp loop for p=2
*vocl loop,novrec
          do 282 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              jmerkd(kdoo)=jmerkd(kdoo)+1
              nwwmo(kdoo,jmerkd(kdoo))=iottr(kstar+npfal(kdoo))
             else
              nx=jposc(iottr(kstar+npfal(kdoo)))
              moafal(kdoo)=moafal(kdoo)-nx
            endif
 282      continue
cvp schleife for p=3
cvp loop for p=3
*vocl loop,novrec
          do 283 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-2
          jdiff(kdoo)=jcon2( iottr(kstar+npfal(kdoo)))
            if(jdiff(kdoo).ne.0) then
              npos(kdoo,1)=iottr(kstar+npfal(kdoo))
            endif
 283      continue
1260     continue
cvp
cvp lauf ueber die offenen schalen von hinten
cvp loop over the open shells starting at the end
cvp
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
cvp jmerko zaehlt die besetzungsunterschiede bzgl. der einfach
c   besetzten mo's von \phi_l, wird hier wieder auf null herunter-
c   gezaehlt
cvp jmerko counts the occupation differences with respect to the
cvp singly occupied mo's of \phi_l, is counted down to zero again here
cvp p=1-fall
*vocl loop,novrec
          do 481 kdoo=1,nspp1
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
          if(jdiff(kdoo).ne.0) then
            jmerko(kdoo)=0
            nwwmo(kdoo,1)=jnopen+1-jdoo
           else
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo-jmerko(kdoo))
          endif
481       continue
cvp p=2-fall
*vocl loop,novrec
          do 482 kdoo=nspp1+1,nsppa
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
cvp       if(jdiff(kdoo).eq.0) then
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo)
cvp       endif
482       continue
cvp p=3-fall
*vocl loop,novrec
          do 483 kdoo=nsppa+1,nsppb
ct        jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))-1
          jdiff(kdoo)=jcon( iottr(kstar+npfal(kdoo)))
cvp       if(jdiff(kdoo).eq.0) then
            nx=jposo(iottr(kstar+npfal(kdoo)))
            nqrfal(kdoo)=nqrfal(kdoo)-
     &       ibinom(nx-1,inopen-1-jdoo)
cvp       endif
483       continue
 1360   continue
cvp bestimmung des r-falles for p=1 sowie des ql-falles
cvp determination of the r-case for p=1 like the ql-case
      do 491 kdoo=1,nspp1
        nqlfal(kdoo)=nwwmo(kdoo,1)
        nwwmo(kdoo,2)=imo3q(nqrfal(kdoo),1)
        nwwmo(kdoo,3)=imo3q(nqrfal(kdoo),2)
        nwwmo(kdoo,4)=imo3q(nqrfal(kdoo),3)
        if (npos(kdoo,1).eq.nwwmo(kdoo,2)) then
           npos(kdoo,1)=3
           if (npos(kdoo,2).eq.nwwmo(kdoo,3)) then
              nrfal(kdoo)=1
             else
              nrfal(kdoo)=2
           endif
           else if (npos(kdoo,1).eq.nwwmo(kdoo,3)) then
             npos(kdoo,1)=2
             nrfal(kdoo)=3
           else
           npos(kdoo,1)=1
        endif
        if (jmerkd(kdoo).eq.1) then
          nrfal(kdoo)=-npos(kdoo,1)
        endif
  491 continue
cvp
cvp  ende labelbestimmung
cvp  end of label determination
cvp
cx    call clockv(vu3,cpu3)
cvp
cvp bestimmung der ww mo's von \phi_r aus qr,
cvp berechnung und sortierung der absoluten mo-nummern for
cvp die integralberechnung
cvp und bestimmung des ql-falles aus den ww mo's von \phi_l
cx       jstar = niot(jsk) + 1 - jnytl
         kstar = niot(jsk) - ncjsk
*vocl loop,novrec
      do 591 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp sei das exchange-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff und jmerko sind hier frei verfuegbar
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp if given a exchange-integral in chemical-notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = idiff(*)
cvp                            b = jdiff(*)
cvp                            a = jmerko(*)
cvp  idiff, jdiff and jmerko are freely available
        nwwmo(kdoo,1)=iottr(kstar+npfal(kdoo)+nwwmo(kdoo,1)*ncjsk)
        nwwmo(kdoo,2)=itest(nwwmo(kdoo,2))
        nwwmo(kdoo,3)=itest(nwwmo(kdoo,3))
        nwwmo(kdoo,4)=itest(nwwmo(kdoo,4))
cvp sortierung for cb- und ex-integrale
cvp sort cb- und ex-integrals
          if (abs(nrfal(kdoo)).ne.3) then
cvp r=1,2,4,5
              if (nwwmo(kdoo,1).lt.nwwmo(kdoo,2)) then
                 nx1=nwwmo(kdoo,1)
                 nwwmo(kdoo,1)=nwwmo(kdoo,2)
                 nwwmo(kdoo,2)=nx1
                 idiff(kdoo)=nwwmo(kdoo,5-abs(nrfal(kdoo)))
                 jdiff(kdoo)=max(nwwmo(kdoo,2),
     &                           nwwmo(kdoo,2+abs(nrfal(kdoo))))
                 jmerko(kdoo)=min(nwwmo(kdoo,2),
     &                           nwwmo(kdoo,2+abs(nrfal(kdoo))))
                else
                 idiff(kdoo)=nwwmo(kdoo,2+abs(nrfal(kdoo)))
                 jdiff(kdoo)=nwwmo(kdoo,2)
                 jmerko(kdoo)=nwwmo(kdoo,5-abs(nrfal(kdoo)))
              endif
            else
cvp r=3 oder r=6
cvp r=3 or r=6
              nx2=nwwmo(kdoo,2)
              nwwmo(kdoo,2)=nwwmo(kdoo,3)
              nwwmo(kdoo,3)=nx2
              if (nwwmo(kdoo,1).lt.nwwmo(kdoo,3)) then
                 idiff(kdoo)=nwwmo(kdoo,2)
                 jdiff(kdoo)=max(nwwmo(kdoo,1),nwwmo(kdoo,4))
                 jmerko(kdoo)=min(nwwmo(kdoo,1),nwwmo(kdoo,4))
                 ny1=nwwmo(kdoo,1)
                 ny2=nwwmo(kdoo,2)
                 nwwmo(kdoo,1)=nwwmo(kdoo,3)
                 nwwmo(kdoo,2)=nwwmo(kdoo,4)
                 nwwmo(kdoo,3)=max(ny1,ny2)
                 nwwmo(kdoo,4)=min(ny1,ny2)
                else
                 idiff(kdoo)=nwwmo(kdoo,4)
                 jdiff(kdoo)=nwwmo(kdoo,3)
                 jmerko(kdoo)=nwwmo(kdoo,2)
              endif
          endif
  591 continue
cvp
c     if (ibs.eq.155) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,3)
c3500 continue
c     endif
cvp p=2 fall
cvp berechnung der absolutnummern der mo's und sortierung
cvp calculation of the absolute mo-numbers and sorting
*vocl loop,novrec
      do 492 kdoo=nspp1+1,nsppa
        if (jmerkd(kdoo).eq.2) then
cvp r=1
          nwwmo(kdoo,3)=itest(inopen+moafal(kdoo))
         else
cvp r=2
          nwwmo(kdoo,3)=nwwmo(kdoo,1)
          nwwmo(kdoo,1)=itest(imoq(nqrfal(kdoo),1))
          nwwmo(kdoo,2)=itest(imoq(nqrfal(kdoo),2))
        endif
cvp mo's sind nun wie folgt sortiert:
cvp  nwwmo(*,3) = doppelt besetztes mo, das in der anderen konf. leer
cvp               ist
cvp  nwwmo(*,1) = einfach besetztes mo
cvp  nwwmo(*,2) = einfach besetztes mo
cvp  nwwmo(*,1) > nwwmo(*,2)
cvp the mo's are now sorted according to:
cvp  nwwmo(*,3) = doubly occupied mo, that is empty in other
cvp               configurations
cvp  nwwmo(*,1) = singly occupied mo
cvp  nwwmo(*,2) = singly occupied mo
cvp  nwwmo(*,1) > nwwmo(*,2)
  492 continue
      do 592 kdoo=nspp1+1,nsppa
        if (jmerkd(kdoo).eq.2)
     &    nqrfal(kdoo)=-nqrfal(kdoo)
  592 continue
cvp p=3 fall
*vocl loop,novrec
      do 493 kdoo=nsppa+1,nsppb
cvp bestimmung der irred. darst. und abspeichern der hoeheren und
c   der niedrigeren mo-nummer innerhalb der irred. darst.
cvp moafal enthaelt groesseres mo
cvp determination of the irrep and storing the higher and lower mo
cvp number within the irrep
cvp moafal holds the higher mo
cvp
        moafal(kdoo)=itest(imoq(nqrfal(kdoo),1))
        mobfal(kdoo)=itest(imoq(nqrfal(kdoo),2))
        if (npos(kdoo,1).eq.mobfal(kdoo))
     &     nqrfal(kdoo)=-1*nqrfal(kdoo)
        nrfal(kdoo)=nirred(mobfal(kdoo))
        moafal(kdoo)=mopos(moafal(kdoo))
        mobfal(kdoo)=mopos(mobfal(kdoo))
cvp codierung des r-falles for p=3 durch vorzeichen von qr
cvp the r-case for p=3 is indicated by the sign of qr
 493  continue
cvp
cvp
cvp  berechnung der integraladresen
cvp  end of label determination
cvp
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3201 kdoo=1,nspp1
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp zugriff auf nit-feld
cvp if given a coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,2)
cvp                            b = nwwmo(*,3)
cvp                            a = nwwmo(*,4)
cvp access of the nit field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
cvp bestimmung des symmetrie-falles (cb-integral)
cvp determination of the symmetry case (cb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum2) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum3) then
cvp      fall 1
c           write(6,*) ' case 1'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks(ideks(idum1)+
     &        mopos(idiff(kdoo)))
     &        +ideks(mopos(jdiff(kdoo)))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=nit( ideks(ideks(ndum1+1))
     &       +ideks(ndum3+1) ) + idum4 +
     &       ideks(ncimo(ndum3)+1)*(ideks(idum1)+idum2-1)
     &       +ideks(idum3)
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &       + mopos(jmerko(kdoo)) +
     &       ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &       +mopos(idiff(kdoo)) )
     &       +ncimo(nirred(idiff(kdoo)))
     &       *(mopos(jdiff(kdoo))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum3) then
cvp      fall 3
c           write(6,*) ' case 3'
            intcb(kdoo)=nit( ideks(idum+ndum2+1) )
     &       +idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
c           if (nrfal(kdoo).ne.2) then
            if (ndum1.ne.nirred(idiff(kdoo))) then
              intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks( ncimo(nirred(idiff(kdoo)))*(idum1-1)
     &         +mopos(idiff(kdoo)) )
     &         +ncimo(nirred(idiff(kdoo)))
     &         *(mopos(jdiff(kdoo))-1)
             else
              intex(kdoo)=nit( ideks(ideks(ndum1+1))
     &         +ideks(nirred(jdiff(kdoo))+1) )
     &         + mopos(jmerko(kdoo)) +
     &         ideks(ncimo(nirred(jdiff(kdoo)))+1)
     &         *(ideks(idum1)+mopos(idiff(kdoo))-1)
     &         +ideks(mopos(jdiff(kdoo)))
            endif
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 4
c           write(6,*) ' case 4'
cvp nur i-1 taucht in den formeln auf
cvp only i-1 shows up in the equations
            idum1=idum1-1
            intcb(kdoo)=nit( ideks(idum+ndum2)
     &       +ideks(ndum3)+ndum4 ) + idum4 +
     &       ncimo(ndum4)*(idum3-1+ncimo(ndum3)*(idum2-1+
     &       ncimo(ndum2)*idum1))
            intex(kdoo)=nit( ideks(idum+nirred(idiff(kdoo)))
     &       +ideks(nirred(jdiff(kdoo)))+nirred(jmerko(kdoo)) )
     &       + mopos(jmerko(kdoo)) +
     &       ncimo(nirred(jmerko(kdoo)))
     &       *(mopos(jdiff(kdoo))-1
     &       +ncimo(nirred(jdiff(kdoo)))
     &       *(mopos(idiff(kdoo))-1
     &       +ncimo(nirred(idiff(kdoo)))*idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
 3201 continue
cvp
cvp bestimmung der gesamtzahl der integrale
cvp  evtl. umspeicherung auf twoe
cvp nintg = anzahl aller integrale
cvp determination of the total number of integrals
cvp  if needed restoring them onto twoe
cvp nintg = estimate for the total number of integrals
      nintg=2*nspp1+nspp2+nspp3*jnopen
     &     +nspp3*(3*jndub+jnopen)
      if (iend+nintg.gt.icmaxm) then
         call isorta(twoe,nteint,ston,nidmax,
     +        nston,nrec31,indub,inopen,lend,iend
     &        ,ncntp3)
         ncntp3=0
         do 85 kdoo=1,icmax
   85    iadp3(kdoo)=0
         jfull=ibs-1
cx       lend=lend+iend
         iend=0
      endif
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
*vocl loop,novrec
      do 120 kdoo=1,nspp1
         iadl(iend+2*kdoo-1)=intcb(kdoo)
         iadl(iend+2*kdoo  )=intex(kdoo)
  120 continue
      iend=iend+2*nspp1
cvp
cvp berechnung der integraladresse for p=2
cvp calculation of the integral addresses for p=2
*vocl loop,novrec
      do 3202 kdoo=nspp1+1,nsppa
cvp das integral sei (ac|bc) oder (ac|cb) usw.
cvp dann ist nwwmo(*,1)=max(a,b)
cvp          nwwmo(*,2)=min(a,b)
cvp          nwwmo(*,3)=c
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the integral is (ab|ac) or (ac|ab) etc.
cvp then     nwwmo(*,1)=a
cvp          nwwmo(*,2)=max(b,c)
cvp          nwwmo(*,3)=min(b,c)
cvp however it should still be sorted according to size
cvp access to nit field
        ndum1=nirred(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        ndum2=nirred(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep.
        idum1=mopos(max(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum2=mopos(min(nwwmo(kdoo,1),nwwmo(kdoo,3)))
        idum3=mopos(max(nwwmo(kdoo,2),nwwmo(kdoo,3)))
        idum4=mopos(min(nwwmo(kdoo,2),nwwmo(kdoo,3)))
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
cvp        fall 1
c           write(6,*) ' case 1'
        intcb(kdoo)=nit( ideks(ideks(ndum1+1)+1) )
     &     + idum4 +
     &       ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)
          else
cvp        fall 3
c           write(6,*) ' case 3'
        intcb(kdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &     + idum4 +
     &       ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1)
        endif
 3202 continue
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
      do 121 kdoo=1,nspp2
         iadl(iend+kdoo)=intcb(nspp1+kdoo)
  121 continue
      iend=iend+nspp2
cvp
cvp  p=3
cvp
cvp  berechnung der summe uber cb- und ex-integrale sowie
cvp  der adressen der ex-integrale bei gemeinsamen offenen
cvp  schalen
cvp  calculate the sum over cb- and ex-integrals and also
cvp  the addresses of the ex-integrals in case of common open shells
cvp
cvp  sei d gemeinsames offenes mo: betrachtung der integrale
cvp     (ab|dd) und (ad|bd), es gilt a > b
cvp  if d is a common open mo: considering the integrals
cvp     (ab|dd) and (ad|bd), where a > b
cvp
cvp nintp3 = anzahl aller summanden for p=3
cvp nintp3 = number of all terms for p=3
      nintp3=nspp3*(3*jndub+2*jnopen)
      do 81 kdoo=1,nspp3
cvp ncntp3 = anzahl der p=3-faelle auf iadl
cvp ncntp3 = number of p=3-cases on iadl
cvp iadp3 = adresse des ersten summanden auf iadl
cvp iadp3 = address of the first term on iadl
         ncntp3=ncntp3+1
         iadp3(ncntp3)=iend+1+(kdoo-1)*(3*jndub+2*jnopen)
   81 continue
cvp
      do 593 jdoo=1,jnopen
*vocl loop,novrec
        do 613 kdoo=nsppa+1,nsppb
          mogd=iottr(kstar+npfal(kdoo)+jdoo*ncjsk)
cvp das ex-integral sei (ad|bd), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ad|bd), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          d=mogd
cvp it still has to be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogd)
cvp bestimmung des symmetrie-falles
cvp determining the symmetrie case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogd))
           idum2=min(moafal(kdoo),mopos(mogd))
           idum3=max(mobfal(kdoo),mopos(mogd))
           idum4=min(mobfal(kdoo),mopos(mogd))
cvp        fall 1
c           write(6,*) ' case 1'
            nwwmo(kdoo,jdoo)=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogd)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)+jdoo)=icb
cvp
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogd)
           idum3=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)+jdoo)=icb
cvp
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogd)
           idum4=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            nwwmo(kdoo,jdoo)=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)+jdoo)=icb
cvp
        endif
  613   continue
  593 continue
cvp schleife ueber die (gemeinsamen) geschlossenen schalen
cvp loop over the (common) closed shells
      do 533 jdoo=1,jndub
*vocl loop,novrec
        do 623 kdoo=nsppa+1,nsppb
c         mogc=iot(iloc-1+inopen+jdoo)
          mogc=iottr(kstar+npfal(kdoo)+(jdoo+jnopen)*ncjsk)
cvp das ex-integral sei (ac|bc), a>b
cvp dann ist a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp es muss aber noch der groesse nach sortiert werden
cvp zugriff auf nit-feld
cvp if the ex-integral is (ac|bc), a>b
cvp then     a=moafal(*)
cvp          b=mobfal(*)
cvp          c=mogc
cvp it should still be sorted according to size
cvp access to nit field
        ndum1=nrfal(kdoo)
        ndum2=nirred(mogc)
cvp bestimmung des symmetrie-falles
cvp determining the symmetry case
        if (ndum1.eq.ndum2) then
           idum1=max(moafal(kdoo),mopos(mogc))
           idum2=min(moafal(kdoo),mopos(mogc))
           idum3=max(mobfal(kdoo),mopos(mogc))
           idum4=min(mobfal(kdoo),mopos(mogc))
cvp        fall 1
c           write(6,*) ' case 1'
            iex=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks(ideks(idum1)+idum2)
     &       +ideks(idum3)+idum4
            ix=ideks(mopos(mogc)+1)
            iy=ideks(moafal(kdoo))+mobfal(kdoo)
            icb=nit( ideks ( ideks(ndum1+1)+1 ) )
     &       +ideks( max(ix,iy) )
     &       +       min(ix,iy)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo-1)=icb
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo  )=icb
cvp negativen summanden auf iadl abspeichern
cvp store negative terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jndub+jdoo  )=iex
cvp
          else if (ndum1.gt.ndum2) then
           idum1=moafal(kdoo)
           idum2=mopos(mogc)
           idum3=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum1)+ndum2+1) )
     &       +ideks(ncimo(ndum2)*(idum1-1)+idum2)
     &       +ncimo(ndum2)*(idum3-1) + idum2
cvp        fall 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum1+1))+ideks(ndum2+1) )
     &       +ideks(ncimo(ndum2)+1)
     &       *(ideks(idum1)+idum3-1)+ideks(idum2+1)
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo-1)=icb
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo  )=icb
cvp negativen summanden auf iadl abspeichern
cvp store negative terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jndub+jdoo  )=iex
cvp
          else if (ndum1.lt.ndum2) then
           idum2=moafal(kdoo)
           idum1=mopos(mogc)
           idum4=mobfal(kdoo)
cvp        fall 3
c           write(6,*) ' case 3'
            iex=nit( ideks(ideks(ndum2)+ndum1+1) )
     &       +ideks(ncimo(ndum1)*(idum1-1)+idum2)
     &       +ncimo(ndum1)*(idum1-1) + idum4
cvp        fall 2
c           write(6,*) ' case 2'
            icb=nit( ideks(ideks(ndum2+1))+ideks(ndum1+1) )
     &       +ideks(ncimo(ndum1)+1)
     &       *(ideks(idum1+1)-1)+ideks(idum2)+idum4
cvp positiven summanden auf iadl abspeichern
cvp store positive terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo-1)=icb
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jdoo  )=icb
cvp negativen summanden auf iadl abspeichern
cvp store negative terms on iadl
          iadl(iend+
     &         (kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &         +jnopen+2*jndub+jdoo  )=iex
cvp
        endif
  623   continue
  533 continue
cvp
cvp umspeichern der integraladressen auf iadl
cvp restoring the integral addresses on iadl
      do 123 jdoo=1,jnopen
        do 122 kdoo=nsppa+1,nsppb
           iadl(iend+(kdoo-nsppa-1)*(3*jndub+2*jnopen)
     &              +3*jndub+jnopen+jdoo)=nwwmo(kdoo,jdoo)
  122   continue
  123 continue
      iend=iend+nintp3
cvp
cvp korrektur der integralsumme
c       do 633 kdoo=nsppa+1,nsppb
c         if (nqrfal(kdoo).lt.0) then
cvp zugriff auf nit-feld
c           ndum1=nrfal(kdoo)
c           idum1=iot(iloc-1+idiff(kdoo))
c           ix=ideks(mopos(idum1)+1)
c           iy=ideks(moafal(kdoo))+mobfal(kdoo)
c           icb=nit( ideks ( ideks(ndum1+1)+1 ) )
c    &       +ideks( max(ix,iy) )
c    &       +       min(ix,iy)
c           sumint(kdoo)=sumint(kdoo)+twoe(icb)
c         endif
c 633   continue
cvp
cvp  ende bessrechnung der integraladressen
cvp
cx    call clockv(vu4,cpu4)
cvp
cx    vutww=vu2-vu1
cx    cpuww=cpu2-cpu1
cx    vutla=vu3-vu2
cx    cpula=cpu3-cpu2
cx    vutin=vu4-vu3
cx    cpuin=cpu4-cpu3
cvp
cvp timer-routine zum test des konfigurationsvergleich
cx    call clockv(vu20,cpu20)
cx    call time(mill20)
cvp
c     if (ibs.eq.4.or.ibs.eq.42) then
c     do 3500 kdoo=1,nspp1
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+npfal(kdoo)*jnytl,
c    &   jstar+npfal(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c       write(6,*) ' intcb= ',intcb(kdoo)
c3500 continue
c     endif
cvp
cvp ruecksortierung wird uebergangen !!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel enthaelt jetzt die adressen der wechselwirkenden
cvp  konfigurationen: zuerst die for p=1, dann for p=2 usw.
cvp  auf npfal wird vorlaeufig der p-fall abgespeichert
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp back-sort is skipped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvp  ispiel now holds the addresses of the interacting
cvp  configurations: first those for p=1, then for p=2 usw.
cvp  provisionaly the p-case is stored on npfal
cvp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 5210 kdoo=1,nspp1
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=2
 5210 continue
      do 5220 kdoo=nspp1+1,nsppa
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=3
 5220 continue
      do 5230 kdoo=nsppa+1,nsppb
        ispiel(kdoo)=npfal(kdoo)
        npfal(kdoo)=4
 5230 continue
cvp
cvp p-faelle werden in urspruengliche reihenfolge gebracht
cvp p-cases being put in original order
cvp
cre    lspp1=1
cre    lspp2=1
cre    lspp3=1
cvp
c      write(6,*) ' jbs=,nspp3=,ispp3= ',jbs,nspp3,(ispp3(kkkk),
c    &   kkkk=1,nspp3)
cvp
cre    do 5100 kdoo=1,nspiel
cre      if (ispiel(kdoo).eq.npfal(lspp1).and.lspp1.le.nspp1) then
cre        ndfal1(kdoo)=2
cre        ndfal2(kdoo)=nrfal(lspp1)
cre        ndfal3(kdoo)=nqlfal(lspp1)
cre        ndfal4(kdoo)=nqrfal(lspp1)
cre        ndfal7(kdoo)=intcb(lspp1)
cre        ndfal8(kdoo)=intex(lspp1)
cre        lspp1=lspp1+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp2+nspp1).and.lspp2.le.nspp2)
cre  &  then
cre        ndfal1(kdoo)=3
cre        ndfal4(kdoo)=nqrfal(lspp2+nspp1)
cre        ndfal7(kdoo)=intcb(lspp2+nspp1)
cre        lspp2=lspp2+1
cre      endif
cre      if (ispiel(kdoo).eq.npfal(lspp3+nsppa).and.lspp3.le.nspp3)
cre  &  then
c          write(6,*) ' kdoo=,ispiel=,ispp3= ',kdoo,ispiel(kdoo)
c    &       ,ispp3(lspp3)
cre        ndfal1(kdoo)=4
cre        ndfal2(kdoo)=nrfal(lspp3+nsppa)
cre        ndfal4(kdoo)=nqrfal(lspp3+nsppa)
cre        ndfal5(kdoo)=moafal(lspp3+nsppa)
cre        ndfal6(kdoo)=mobfal(lspp3+nsppa)
cre        ndfeld(kdoo,1)=nwwmo(lspp3+nsppa,1)
cre        ndfeld(kdoo,2)=nwwmo(lspp3+nsppa,2)
cre        ndfeld(kdoo,3)=nwwmo(lspp3+nsppa,3)
cre        ndfeld(kdoo,4)=nwwmo(lspp3+nsppa,4)
cre        ndfeld(kdoo,5)=nwwmo(lspp3+nsppa,5)
cre        ndfeld(kdoo,6)=nwwmo(lspp3+nsppa,6)
cre        ndfeld(kdoo,7)=nwwmo(lspp3+nsppa,7)
cre        ndfeld(kdoo,8)=nwwmo(lspp3+nsppa,8)
cre        ndfeld(kdoo,9)=nwwmo(lspp3+nsppa,9)
cre        rdum(kdoo)=sumint(lspp3+nsppa)
cre        lspp3=lspp3+1
cre      endif
c5100  continue
 250    continue
      return
      end
cvp
      subroutine dk2int(iottr,iotnew,iotm,
     +           maindf,maxr,jconb,nomax2,
     +           isk,jsk,ibs,nspiel,mspiel,jnum1)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer iottr,iotnew,iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c speziell for dk=2
c jmerk laufen direkt vom laeufer angesteuert
c dadurch erfolgt der zugriff schneller
c geaenderte version: nur nummern der wechselwirkenden konf. zu einer
c bei gleichen sk erfolgt der vergleich bzgl. der unteren dreieckmatrix
c (vp 13.3.1992)
c especially for dk=2
c jmerk passes driven directly by the loop counter
c this speeds up the memory access
c modified version: just the numbers of interacting configurations
c in case of the same sk the comparison regards the lower triangle
c (vp 13.3.1992)
c  festen konf. werden bestimmt
c  uebergeben werden: iswh=anzahl der sk
c                     isk=linke sk
c                     jsk=rechte sk
c                     ibs=nummer der konf. aus der sk isk
c  resultate sind:    nspiel=anzahl der ww konf.
c                     ispiel=vektor, der die nummern der ww konf.
c                            enthaelt
c  determining permanent configurations
c  inputs are       : iswh=number of sk
c                     isk=left sk
c                     jsk=right sk
c                     ibs=number of configurations from the sk isk
c  outputs are   :    nspiel=number of interacting configurations
c                     ispiel=vector containing the numbers of the
c                            interacting configurations
cvp
cvp
cx    ncisk=nconf(isk)
cx      inytl  = nytl(isk)
cx      indub  = ndub(isk)
cx      inopen = inytl - indub
cx      istar=niot(isk)+1
cx    iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
cx      do 215 k=1,imo
cx      jcon1(k)=1
cx      jcon2(k)=2
cx215   continue
cvp abspeichern der testkonfiguration auf itest
cx       do 205 k=1,inytl
cx205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
cx       do 220 k=1,inopen
cx       jposo(itest(k))=k
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 1
cx220    continue
cvp
cx       do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
cx       jcon1(itest(k)) = 0
cx       jcon2(itest(k)) = 0
cx221    continue
c beginn do loop ueber konfig aus ft36, gleich und sk+1 ist interessant
cx       if(jsk.lt.1) go to 250
cx       jnytl=nytl(jsk)
cx       jndub=ndub(jsk)
cx       jnopen=jnytl - jndub
cx       ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
cvp vergleich bzgl. unterer dreiecksmatrix
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cx       nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nullsetzen von idiff feld soweit notwendig
cx       do 255 ialex = 1,nspiel
cx55     idiff(ialex) = 0
c test der offenen schalen
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
cx       jstar=niot(jsk)
cvp jstar wird um die nummer der ersten zu testenden konfiguration
cvp  -1 hochgesetzt
cx       jstar=jstar+jnum1
cvp
cx       kstar=jstar+jnopen*ncjsk
cx       do 360 jdoo = 1,jnopen
cx        kstar=kstar-ncjsk
*vocl loop,novrec
cx        do 380 kdoo = 1,nspiel
cx        idiff(kdoo) = idiff(kdoo) + jcon1( iottr(kstar+kdoo) )
cx0       continue
cx0      continue
c do loop ueber die geschlossenen schalen
c beginn hinten, da dort die meisten unterschiede zu finden sind
cx       kstar=jstar+jnytl*ncjsk
cx       do 260 jdoo = 1,jndub
cx        kstar=kstar-ncjsk
c unterschiede zur testkonfiguration auf idiff summieren
cx        do 270 kdoo = 1,nspiel
cx        if (idiff(kdoo).le.2) then
cx        idiff(kdoo) = idiff(kdoo) + jcon2( iottr(kstar+kdoo) )
cx        endif
cx0       continue
cx0      continue
cvp
      ncisk=nconf(isk)
        inytl  = nytl(isk)
        indub  = ndub(isk)
        inopen = inytl - indub
        istar=niot(isk)+1
      iloc=istar+(ibs-1)*inytl
c do loop ueber die konfigurationen bzgl derer getestet werden soll
c nullsetzen des jcon feldes
c do loop over those configurations with respect to which the tests
c should be performed
c nullify the jcon field
        do 215 k=1,imo
ct      jcon1(k)=4
        jcon1(k)=1
        jcon2(k)=2
215     continue
cvp abspeichern der testkonfiguration auf itest
cvp store the test configuration on itest
         do 205 k=1,inytl
  205    itest(k)=iottr(istar+(k-1)*ncisk+ibs-1)
         do 220 k=1,inopen
         jposo(itest(k))=k
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 1
  220    continue
cvp
c        write(6,*) ' sk=',isk
c        write(6,*) ' phi_r: ',(iot(k),k=iloc,iloc+inytl-1)
cvp
         do 221 k=inopen+1,inytl
c zuerst mal besetzen von jcon1, sodass ausgehend von einzel und finden
c von doppelt nicht gezaehlt wird
c first occupy jcon1, such that starting from a single and finding a
c a double is not counted.
         jposc(itest(k)) = k-inopen
ct       jcon1(itest(k)) = 4
         jcon1(itest(k)) = 0
         jcon2(itest(k)) = 0
  221    continue
c        write(6,*) 'jcon for konf#',ibs,' in sk ',isk
c        write(6,*) (jcon(iaus),iaus=1,imo)
c beginn do loop ueber konfig aus ft36, gleich und sk+1 ist interessant
c begin do over configurations from ft36. The same and sk+1
c is interesting
         if(jsk.lt.1) go to 250
         jnytl=nytl(jsk)
         jndub=ndub(jsk)
         jnopen=jnytl - jndub
         ncjsk = nconf (jsk)
c aufbau des feldes ispiel ( nr. der mitspieler) ausblenden alle jener
c die vorher schon als zweifachanregung identifiziert worden waren
c building field ispiel ( nr. of players) eliminate all those
c that have been identified as double excitations already
c
cvp vergleich bzgl. unterer dreiecksmatrix
cvp comparison regarding the lower triangle
cx       nspiel=ncjsk
cvp nspiel ist hier die max. anzahl zu testender konfigurationen,
cvp  die auf die vektoren passen
cvp here nspiel is the max. number of configurations to be tested
cvp  that fit onto the vectors
         nspiel=mspiel
c ispiel : nr innerhalb der sk der spielenden konfigurationen
c ispiel : number of playing configurations within the sk
c nspiel : anzahl der mitspieler innerhalb der sk jsk
c nspiel : number of players within the sk jsk
c
c bestimmung der wechselwirkenden konfigurationen mit hilfe
c von erzeugungs- und vernichtungsoperatoren
c determining the interacting configurations using
c creation and annihilation operators
      call chkww(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +           isk,jsk,ibs,nspiel,jnum1,ncisk,ncjsk)
c
cbs
cbs jetzt wird komprimiert
cbs now compress
cbs
          nspie1=0
*vocl loop,novrec
          do 385 kdoo = 1,nspiel
           if(idiff(kdoo).le.2) then
            nspie1 = nspie1 + 1
cx         ispiel(nspie1) = kdoo
           ispiel(nspie1) = kdoo + jnum1
           idiff(nspie1)=idiff(kdoo)
           endif
 385      continue
          nspiel=nspie1
cvp jstar wird wieder auf alten wert zurueckgesetzt
cvp jstar is reset to the old value
         jstar=niot(jsk)
cvp
cvp  ende bestimmung der wechselwirkenden konfigurationen
cvp  end of interacting configuration determination
cvp
cbs   und jetzt das ganze noch mal
cbs   and now do the whole thing once more
           do 261 jbs=1,nspiel
           nqrfal(jbs)=ibinom(inopen,4)
 261       jmerkd(jbs)=0
c beginn hinten, da dort die meisten unterschiede zu finden sind
c start at end because most differences will be found there
         kstar=jstar+jnytl*ncjsk
         do 1260 jdoo = 1,jndub
          kstar=kstar-ncjsk
c unterschiede zur testkonfiguration auf idiff summieren
c sum differences with respect to test configuration onto idiff
cvp jdiff(maximal) ist for doppelt besetzte mo's genau 2 for dk=2
cvp npos(*,1) gibt position des hinteren unterschiedes, npos(*,2) des
cvp vorderen
cvp jdiff(maximal) is exactly 2 for doubly occupied mo's with dk=2
cvp npos(*,1) gives the position of the last difference,
cvp npos(*,2) the position of the first
*vocl loop,novrec
          do 286 kdoo=1,nspiel
ct          jdiff(kdoo)=jcon( iottr(kstar+ispiel(kdoo)))-2
            jdiff(kdoo)=jcon2( iottr(kstar+ispiel(kdoo)))
            if(jdiff(kdoo).ne.0) then
             jmerkd(kdoo)=jmerkd(kdoo)+1
             npos(kdoo,jmerkd(kdoo))=jposo(iottr(kstar+ispiel(kdoo)))
            endif
 286      continue
1260     continue
cvp
cvp bestimmung der in \phi_r einfach und in \phi_l nicht besetzten
cvp determination of the in \phi_r singly occupied and in \phi_l
cvp unoccupied mo's
c test der offenen schalen
c test of open shells
c        write(6,*) 'jcon for konf#',ibs,' in sk ',isk
c        write(6,*) (jcon(iaus),iaus=1,imo)
c ebenso wie bei den geschlossenen beginn bei der letzten offenen
c like with the closed shells start at the last open shell
cvp nx = position der gemeinsamen offenen schale (letzte offene von
cvp  \phi_l)
cvp nx = position of the common open shell (last open shell of \phi_l)
         kstar=jstar+jnopen*ncjsk
         do 1360 jdoo = 1,jnopen
          kstar=kstar-ncjsk
c differenz bilden mit jcon
c compute difference with jcon
cvp nx = position der gemeinsamen offenen schale
cvp nx = position of the common open shell
*vocl loop,novrec(jmerkd),vi(kdoo)
          do 1380 kdoo = 1,nspiel
          nx=jposo(iottr(kstar+ispiel(kdoo)))
cvp bestimmung des qr-falles
cvp determination of the qr-case
          nqrfal(kdoo)=nqrfal(kdoo)-
     &     ibinom(nx-1,jnopen+1-jdoo)
          if (nx.lt.npos(kdoo,1)) npos(kdoo,1)=npos(kdoo,1)-1
          if (nx.lt.npos(kdoo,2)) npos(kdoo,2)=npos(kdoo,2)-1
 1380     continue
 1360     continue
cvp
cvp bestimmung des r-falles aus npos
cvp determination of the r-case from npos
      do 3010 kdoo=1,nspiel
        nrfal(kdoo)=npos(kdoo,1)-npos(kdoo,2)
        if (nrfal(kdoo).eq.1.and.npos(kdoo,1).eq.3) nrfal(kdoo)=3
cvp     nrfal(kdoo)=idk2r(npos(kdoo,1),npos(kdoo,2))
 3010 continue
cvp
cvp  ende labelbestimmung
cvp  end of labeldetermination
cvp
cvp
cvp  berechnung der integraladresen
cvp  integral address calculation
cvp
      iloc1=iloc-1
      iloc2=iloc+ibs-1
*vocl loop,temp(ndum1,ndum2,ndum3,ndum4,idum1,idum2,idum3,idum4,idum)
      do 3200 kdoo=1,nspiel
cvp uebergang von positionen auf absolute mo-nummern
cvp change over from positions to absolute mo-numbers
        nwwmo(kdoo,1)=itest(imo4q(nqrfal(kdoo),1))
        nwwmo(kdoo,4)=itest(imo4q(nqrfal(kdoo),4))
        if (nrfal(kdoo).eq.1) then
          nwwmo(kdoo,3)=itest(imo4q(nqrfal(kdoo),3))
          nwwmo(kdoo,2)=itest(imo4q(nqrfal(kdoo),2))
        endif
        if (nrfal(kdoo).ne.1) then
          nwwmo(kdoo,3)=itest(imo4q(nqrfal(kdoo),2))
          nwwmo(kdoo,2)=itest(imo4q(nqrfal(kdoo),3))
        endif
cvp sei das coulomb-integral in chemiker-notation (d, c | b, a),
cvp dann gilt die zuordnung:   d = nwwmo(*,1)
cvp                            c = nwwmo(*,3)
cvp                            b = nwwmo(*,2)
cvp                            a = nwwmo(*,4)
cvp zugriff auf nit-feld
cvp given the coulomb-integral in chemical notation (d, c | b, a),
cvp then assign            :   d = nwwmo(*,1)
cvp                            c = nwwmo(*,3)
cvp                            b = nwwmo(*,2)
cvp                            a = nwwmo(*,4)
cvp access nit field
        ndum1=nirred(nwwmo(kdoo,1))
        ndum2=nirred(nwwmo(kdoo,2))
        ndum3=nirred(nwwmo(kdoo,3))
        ndum4=nirred(nwwmo(kdoo,4))
cvp uebergang zu mo-nummern innerhalb der jeweiligen irred. darst.
cvp change over to mo-numbers within the current irrep
        idum1=mopos(nwwmo(kdoo,1))
        idum2=mopos(nwwmo(kdoo,2))
        idum3=mopos(nwwmo(kdoo,3))
        idum4=mopos(nwwmo(kdoo,4))
        idum=ideks(ndum1)
        intcb(kdoo)=nit( ideks(idum+ndum3)
     &    +ideks(ndum2)+ndum4 ) + idum4
        intex(kdoo)=nit( ideks(idum+
     &    nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &    +ideks(nirred(nwwmo(kdoo,nrb(nrfal(kdoo)))))
     &    +nirred(nwwmo(kdoo,nrc(nrfal(kdoo)))) )
     &    +mopos(nwwmo(kdoo,nrc(nrfal(kdoo))))
cvp bestimmung des symmetrie-falles (cb-integral)
cvp determination of the symmetry case (cb-integral)
*vocl stmt,if(80)
        if (ndum1.eq.ndum3) then
*vocl stmt,if(60)
          if (ndum1.eq.ndum2) then
cvp      fall 1
c           write(6,*) ' case 1'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ideks(idum1)+idum3)
     &       +ideks(idum2)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ideks(idum1)+
     &        mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &        +ideks(mopos(nwwmo(kdoo,nrb(nrfal(kdoo)))))
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 2
c           write(6,*) ' case 2'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum2)+1)*(ideks(idum1)+idum3-1)
     &       +ideks(idum2)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(idum1-1)
     &       +mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
         else
*vocl stmt,if(80)
          if (ndum1.eq.ndum2) then
cvp      fall 3
c           write(6,*) ' case 3'
            intcb(kdoo)=intcb(kdoo)+
     &       ideks(ncimo(ndum3)*(idum1-1)+idum3)
     &       +ncimo(ndum3)*(idum2-1)
            intex(kdoo)=intex(kdoo)+
     &       ideks(ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(idum1-1)
     &       +mopos(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1)
c           write(6,*) ' intcb= ',intcb(kdoo)
           else
cvp      fall 4
c           write(6,*) ' case 4'
cvp nur i-1 taucht in den formeln auf
cvp only i-1 shows up in the equations
            idum1=idum1-1
            intcb(kdoo)=intcb(kdoo)+
     &       ncimo(ndum4)*(idum2-1+ncimo(ndum2)*(idum3-1+
     &       ncimo(ndum3)*idum1))
            intex(kdoo)=intex(kdoo)+
     &       ncimo(nirred(nwwmo(kdoo,nrc(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nrb(nrfal(kdoo))))-1
     &       +ncimo(nirred(nwwmo(kdoo,nrb(nrfal(kdoo)))))
     &       *(mopos(nwwmo(kdoo,nra(nrfal(kdoo))))-1
     &       +ncimo(nirred(nwwmo(kdoo,nra(nrfal(kdoo)))))
     &       *idum1))
c           write(6,*) ' intcb= ',intcb(kdoo)
          endif
        endif
c     write(6,*) ' ibs=,kdoo=,intex= ',ibs,kdoo,intex(kdoo)
 3200 continue
cvp   do 3500 kdoo=1,nspiel
c       write(6,*) ' phi_l: ',(iot(k),k=jstar+ispiel(kdoo)*jnytl,
c    &   jstar+ispiel(kdoo)*jnytl+jnytl-1)
c       write(6,*) ' ww mo: ',(nwwmo(kdoo,k),k=1,4)
c3500 continue
cvp
cvp
cvp  ende berechnung der integraladressen
cvp  finished calculating the integral addresses
cvp
cvp
250     continue
      return
      end
c--------1---------2---------3---------4---------5---------6---------7--
c
c reference : g. cisneros, m .berrondo, c. f. bunge
c             computers and chemistry vol.10; issue 4; p 281 (1986)
c
c for ibm rs6000:
c                 s. grimme, uni bonn, 9.92
c
c erweitert auf gleichzeitige behandlung mehrerer wurzeln
c nach einem version von
c extended to calculate multiple roots at the same time
c based on a version from
c       b. liu
c       numerical algorithms in chemistry:algebraic methods,
c       lbl-8158 lawrence berkeley laboratory
c       eds.: c. moler and i. shavitt
c
c                 b. engels, uni bonn, 12.93
c
c
cbe      um mit rglob ins reine zu kommen     25.2.94
cbe ideks in iideks veraendert um es vom rest des programms zu trennen
cbe mi in mi1 veraendert
cbe      to get in the clear with rglob       25.2.94
cbe ideks changed into iideks to separate it from the rest of the 
cbe  program
cbe mi changed into mi1
c-----------------------------------------------------------------------
c
      subroutine dvdmua(twoe, nteint,
     +    pey, acoul, aexc, ndeke, core,
     +    vecf2, vecf1, mdi, hp5, mdi0,
     +    emat, nedim,
     +    iottr,iotnew,iotm,
     +    maindf,maxr,jconb,nomax2,
     +    nnl, minl, mfnl, nivnl,
     1    noutnl, lunnl, lun1nl, lun2nl,
     2    critan, critrn, criten, orthon,
     3    adiag, v, w, e,
     +    av, h, u, mx1, mx, iwr, odebug)

c      fill the h-matrix in subr. mvmpy
c      start vectors are on unit lun
c
c  izus enthaelt welche eigenwerte gleichzeitig behandelt werden
c  sollen
c  izus tells which roots should be treated simulaneously
cbe include(adlinc) durch include (rglob) ersetzt
cbe     25.2.94
cbe include(adlinc) replaced by include (rglob)
cbe     25.2.94
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      integer nteint
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      real*8 vecf1,vecf2,hp5
      integer mdi, mdi0
      dimension vecf2(mdi), vecf1(mdi), hp5(mdi0)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      dimension  adiag(nnl), v(nnl), w(nnl), e(mxroot)
c
c maximum size of the projected matrix = mx
c
c     parameter (mx=256, mx1=(mx*mx+mx) / 2)
c
      dimension av(mx1), h(mx1), u(mx,mx)
c
      real*8 orthon, ortho
c u : enthaelt entwicklungskoeffizienten nach der diagonalisierung
c     der subraeume
c u : holds the coefficients after the diagonalisation of the 
c     sub-spaces
c d : im hqrii1
c d : in hqrii1
c v : hilfsvektor
c v : help vector
c w : hilfsvektor
c w : help vector
c av: submatrix des davidson versions
c av: submatrix of the davidson versions
c e : energie des vorherigen laufes, 
c     rueckgabe der energie an hauptprogram
c e : energy of the previous cycle,
c     return value of the energy to the main program
      common /rwork2/ ndim,nanz,nglein
c !!!! change dimension of work array w(5,mx) in hqrii1 !!!!
c
c
      real*8 d, valn, valn2, resnm, valoln
      integer mm, mrm, iideks, ineue, iconv, mxl
      integer mxrootn
      parameter (mxrootn = 50)
c
      parameter (mxl=mxrootn*mxrootn)
c     parameter (mxl=256)
      common/miscop/d(mxl),valn(mxl),valn2(mxl),resnm(mxl),
     +              valoln(mxl),
     +       mm(mxl),mrm(mxl),iideks(mxl),ineue(mxl),iconv(mxl)
c
      logical odebug
c
c --- scf information
      real*8 escf,vnuc
      common /kerne/ escf,vnuc
      common /cselec/ iselec(mxroot),izus(mxroot)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      parameter (nzero = 0)
      common /csemi/ hstor,issk
cvp
cvp for integralvorsortierung nach foxy
cvp for integral presort to foxy
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
cvp
*     logical*4 od12,ex12
*     character*10 acc12,seq12,dir12
c
      data big /1.23456d+20/
c NB!!!
c cepai: if cepa calculation is requested
c icepa: specifies the cepa variant
c ipcepa: print level
c 
      logical cepai
      integer icepa,ipcepa
      common/cepa_mrd1/cepai,icepa,ipcepa
c
c crtest is the criterion for the tester for invoking cepa
c crres is the criterion for the residue norm for invokinf cepa
c enaught is the reference energy to be used for shift calculation
c corren is a correlation energy to be use for shift calculation
c crshift used to suppress shift calculation every iteration
c
      logical crcorr,cre0
      real*8 crtest,crres,enaught,corren
      common/cepa_mrd2/crtest,crres,enaught,corren,crcorr,cre0
c
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
      common /parkin/egey1,trash,tdel,ical0,nele
c
      real*8 ui0(mxl),shift
      data ui0/mxl*0.0d0/
      data shift/0.0d0/
c --- tester residue norm and crtce for ci and cepa
      real*8 tester(30),testmx,resmx
      logical crtce
      character*10 charwall
      data crtce/.false./
c     data crtest,crres/2*1.0d-2/
c----
c NB
      write(iwr,*)
      write(iwr,*) '==========================='
      write(iwr,*) '=== multi-root davidson ==='
      write(iwr,*) '==========================='
      if (odebug) then
       write(iwr,*) 'iselec = ',(iselec(iraus),iraus=1,7)
       write(iwr,*) 'izus = ', (izus(iraus),iraus=1,7)
      endif
c --- spater weg!
      ndim = nnl
c tresholds for debuggen
      trawf=0.2d0
      trahb=0.2d0
      tragr=0.7d0
      trakl=0.05d0
      trami=0.005d0

c aufbau des iideks feldes
      iideks(1) = 1
      do   idoof = 2,mx
         iideks(idoof) = iideks(idoof - 1) + idoof
      enddo

c--- setup
      n=nnl
c     write(iwr,*) 'dimension of matrix', n
c n : dimension of matrix
      mi1=minl
c mi1: index of lowest eigenvalue to be calc.
      mf=mfnl
c mf:  "       highest    "
      niv=nivnl
c niv :number of start vectors
      if(odebug) then
       write(iwr,*) 'n ',n,' mi1',mi1,' mf',mf,' niv',niv
      endif
c niv ge mf
c change iord in hqrii1 to 1 if decreasing order is requested
      nout=noutnl
c output (usually 6)
      lun=lunnl
c logical unit number with initial eigenvectors
      lun1=lun1nl
      lun2=lun2nl
cbe fors debuggen
*      write(iwr,*) 'out dvdmua : file check'
*      inquire(unit=lun,exist=ex12,opened=od12,access=acc12,
*     *sequential=seq12,direct=dir12)
*      write(iwr,27) lun,ex12,od12,acc12,seq12,dir12
*      inquire(unit=lun1,exist=ex12,opened=od12,access=acc12,
*     *sequential=seq12,direct=dir12)
*      write(iwr,27) lun1,ex12,od12,acc12,seq12,dir12
*      inquire(unit=lun2,exist=ex12,opened=od12,access=acc12,
*     *sequential=seq12,direct=dir12)
*      write(iwr,27) lun2,ex12,od12,acc12,seq12,dir12
*27    format(1x,'for logical consistency ',i4,/,
*     *'exist,open',2(l3,1x),/,
*     *'access,seq,dir ',3(a10,2x))
c scratch
      crita=critan
      critr=critrn
      crite=criten
c termination crit.
      ortho=orthon
c--- input check
c     write(*,*)'n',n,'mi1',mi1,'mf',mf,'niv',niv,
c    .'nout',nout,'lun lun1 lun2',lun,lun1,lun2
c     write(*,*)'crit',crita,crite,critr,ortho
c     write(nout,'(/10x,''davidson diagonalization'')')
c     write(nout,'(/1x,''iter'',10x,''eigenvalue no.'',10x,
c    .''energy'',10x,''alpha'',10x,''res. norm'',/)')
c
c---
cbe groesse des startraumes ist groesser als die erlaubte
cbe groesse des projektionsraumes
cbe size of the initial space larger than the allow size of the
cbe projection space
      if(niv .gt. mx) goto 1001

cbe wieviel wurzeln und welche wurzeln
cbe herausfinden durch analyse von iselec
cbe how many roots and which ones
cbe find out by analysing iselec
         if (iselec(1).ne.nzero) then
           neig=1
           do ibe = 2,mxroot
             if(iselec(ibe).ne.nzero) neig=neig+1
           enddo
         endif
cbe checken ob dies mit dem input uebereinstimmt
cbe check whether this is consistent with the input
         if (iselec(neig).ne.mf) then
           write(iwr,*) 'iselec does correspond with mf '
         endif
cbe test output
        if (odebug) then
         write(iwr,*) 'from dvdmua.f'
         write(iwr,*) 'number of roots ordered : ',neig
         write(iwr,*) (iselec(iaus),iaus=1,neig)
        endif

c---- end of input check
c
c     berechnung wieviele gleichzeitig behandelt werden koennen
c     calculate how many roots can be treated simultaneously
*****
*****  force driving through fixed passes through data input
*****  if (nglein.gt.0) then
*****   nglei = nglein
*****  else
*****   nglei = int(mdi/n)
*****  endif
*****  reset nglein for reduced n at higher threshold
*****
*****  nglein = 0
*****
       nglei = int(mdi/n)
       if(odebug) then
        write(iwr,*) 'n=',n
        write(iwr,*) 'nglei=',nglei
        write(iwr,*) 'iselec',(iselec(iraus),iraus=1,5)
       endif
c
c---- subroutine startm macht
c     aufbau der startmatrix av = bi*a*bj mit bi, bj startvektoren
c     for die niv (=hoechste wurzel) untersten startvektoren
c---- subroutine startm does
c     build the start matrix av = bi*a*bj where bi, bj start vectors
c     for the niv (=highest root) lowest start vectors
c
c     print *,'n, niv, lun, lun1, lun2, nout,
c    1                  ortho, nglei',n, niv,
c    2                 lun, lun1, lun2, nout,
c    1                  ortho, nglei
      call startm(twoe, nteint,
     +    pey, acoul, aexc, ndeke, core,
     +    vecf1, vecf2, mdi, hp5, mdi0,
     +    emat, nedim, 
     +    iottr, iotnew, iotm,
     +    maindf, maxr,jconb,nomax2,
     +    n, niv, lun, lun1, lun2, nout,
     +    ortho, nglei, v, w, av,iwr,odebug,nko)
c      write(iwr,*) 'after startm'
*      writeiwr6,*) 'stop after startm'
*      stop
c n      groesse der matrix
c n      size of the matrix
c niv    anzahl der startvektoren
c niv    number of start vectors
c lun    anfangs logische einheit die startvektoren enthaelt
c        am schluss enthaelt lun die fertigen vektoren
c lun    in: the start vectors
c        out: the eigen vectors
c lun1   logische einheit die orthonormierte vektoren enthaelt
c lun1   holds orthonormalised vectors
c lun2   logische einheit die a*v enthaelt
c lun2   holds a*v
c ortho  orthogonalisierungs kriterium
c ortho  orthogonalisation criterium
c nglei  wieviele gleichzeitig behandelt werden koennen
c        rueckgabe der information
c nglei  the number of vectors that can be treated simultaneously
c av     enthaelt die matrix bi*a*bj mit bi,bj startvektoren
c av     the matrix bi*a*bj where bi,bj start vectors
c
c---- aufbau der startmatrix fertig
c---- build of start matrix finished
c hiernach sind sie in der niv-ten iteration, d.h. es existieren
c j eigenvektoren
c here after you are in the niv-th iteration, i.e. there are
c j eigen vectors
      j = niv
      jw = (j*j-j)/2

c---- debug ausdruck der startmatrix
c---- debug print of startmatrix
      jende = j*(j+1)/2
      jraus = 10
      if(odebug) then
       write(iwr,*) 'the first ',jraus,' elements of the start matrix'
       write(iwr,*) (av(iaus),iaus=1,jraus)
      endif
c---- programm bis hierhin fertig !!!!!!!
c---- program finished upto here !!!!!!!
c
c---- aufbau von izus
c---- build izus
      if (odebug) then
       write(iwr,*) 'preparation'
       write(iwr,1) (izus(iaus),iaus=1,neig)
      endif
c izus(1) = -1 alle einzeln
c izus(1) = -1 all separately
c izus(1) =  0 beliebig
c izus(1) =  0 do as you please
c izus(1) >  0 izus gibt an welche gemeinsam
c izus(1) >  0 izus tells which ones to treat simultaneously
      if(izus(1).eq.-1) then
c alle einzeln
c all separately
        do    ilauf = 1,neig
          izus(ilauf) = 1
        enddo
          idurch = neig
      else if (izus(1) .eq. 0) then
c beliebig viele gleichzeitig
c as many simultaneously as you fancy
        idurch = int(neig/nglei) + 1
        do    ilauf = 1,idurch - 1
          izus(ilauf) = nglei
        enddo
        izus(idurch) = mod(neig,nglei)
      else if (izus(1) .gt. 0 ) then
c wieviele jeweils gleichzeitig steht in izus
c the number of states treated simultaneously specified by izus
        idurch = 0
        do ilauf = 1,neig
          if(izus(ilauf) .gt. 0 ) idurch = idurch + 1
        enddo
      endif
c rwah
      if (odebug)
     + write(iwr,1) (izus(iaus),iaus=1,neig)
c     write(iwr,2)

c----
      istart  = 0
c---- aufbau eines feldes zum checken der konvergenz
c---- build a field to check the convergence
      do     idoof = 1,neig
        iconv (iselec(idoof)) = 1
      enddo
      write(iwr,3) (iselec(iraus), iraus = 1,neig)

c-- nullsetzen der ci energien
c-- nullify the ci energies
        do    lauf = 1,neig
            e(iselec(lauf)) = 0.0d0
        enddo

c do 9000 ist groesse loop ueber die einzelnen durchlaeufe
c do 9000 is a big loop over separate iterations

c     write(iwr,*) 'before do 9000:idurch',idurch
      do 9000   ldurch = 1,idurch
        iende = istart + izus(ldurch)
        istart = istart + 1
        if (odebug)write(iwr,4) ldurch
        do    lauf = istart,iende
           mrm(lauf) = lauf
           mm (lauf) = iselec(lauf)
        enddo
        write(iwr,8) (mm(ilauf),ilauf = istart,iende)
        write(nout,6) cpulft(1) ,charwall()


      iter = 0
c eigentliche schleife im davidson
c the real loop in the davidson
50    continue
*      write(iwr,*)'e after 50',e(1),e(2),e(3)

*      write(iwr,*) 'the size of the dav.-matr. is',j
*        if (j.gt.7) then
*          write(iwr,*) 'program is before the ',j,' -th iteration'
*        stop
*        endif
c      write(nout,666) iter,cpulft(1)
c666   format('cpu time in',2x,i2,2x,'interation',f10.2,2x,'sec')
       iter = iter + 1
c for single root
      if( j .eq. 0 ) then
c
c hier muss noch geaendert werden
c must still be changed
c
          u(1,1)=1.0d0
          val=av(1)
          e(1)=val
      else
           if(odebug) then
            write(6,*) 'jw and j',jw,j
           endif
          idoof = jw + j
c      write(iwr,*)'e before copyv',e(1),e(2),e(3)
c         call copyv(jw+j,av,h)
          call dcopy(idoof,av,1,h,1)
c      write(iwr,*)'e after copyv',e(1),e(2),e(3)
           if(odebug) then
           write(iwr,*) 'h is '
           write(iwr,*) (h(iaus),iaus= 1,idoof)
           write(iwr,*) 'call hqrii1 for diagonalisation '
          endif
c
c NB  (D + Dshift for cepa)
c
       if (cepai.and.testmx.lt.crtest.and.resmx.lt.crres
     &     .and.j.gt.3*niv) crtce=.true.
c
       if(ipcepa.ge.2.and.j.eq.niv) then
          write(iwr,515)crtest,crres,enaught,corren
 515   format(//2x,'--- cepa input ---',
     &        //2x,'tester crit:',1x,f8.4,
     &        /2x,'res. norm crit.:',1x,f8.4,
     &        /2x,'E0(specified):',1x,f25.16,
     &        /2x,'Ecorr for shift calc.:',f20.16,//)
       end if
c
       if (crtce) then
          if (crcorr) then
             shift = corren
          else if (cre0) then
             shift = d(1) - enaught
          else
             shift = d(1) - av(1)
          end if
          shift = shift*corrce(nele,icepa)
          call cepash(j,h,ui0,shift)
       end if
c
c NB
c
       call hqrii1(j,j,h,d,u,mx)
c      write(iwr,*)'e after hqrii1',e(1),e(2),e(3)
          if(odebug) then
           do ialex = 1,j
              write(iwr,*) 'd(',ialex,') = ',d(ialex)
           enddo
          endif
          do    lauf = istart,iende
c berechnung der energie (ohne escf)
c calculate the energy (without escf)
              valn(lauf) = d(iselec(lauf))
              if(odebug) then
                 write(iwr,*) 'd(',iselec(lauf),' = ',valn(lauf)
              endif
c berechnung der richtigen energie ohne shift
c calculation of the correct energy without shift
              valn2(lauf) = +d(iselec(lauf)) + escf
          enddo
*          do     lauf = istart,iende
*              write(iwr,*) 'lauf =',lauf
*              write(iwr,*) 'energie  : ',valn(lauf),valn2(lauf)
*              write(iwr,*) 'total-energy : ',valn2(lauf)
*          enddo
*       write(iwr,*)'e after hqrii1',e(1),e(2),e(3)
      endif
c

c----  projected matrix possesses equal size than total matrix
c tut noch nicht ****************
c does not work yet *************
      if( j .eq. n ) then
          call rewftn(lun)
          ir=0
           do 170 ii=1,neig
           i=iselec(ii)
           ir=ir+1
           e(ir)=d(i)
           call rewftn(lun1)
           call zro1d(n,v)
           do 160 k=1,j
              read (lun1) w
              uki=u(k,i)
              call vsma1d(n,uki,w,v)
  160      continue
           write (lun) v
  170      continue
           call rewftn(lun)
           close (lun1,status='delete')
           close (lun2,status='delete')
           return
      endif
c
c----   eigentliche bearbeitung der vektoren
c----   real processing of the vectors
c
c aufbau der vektoren
c ausnullen des in vecf1 notwendigen platzes
c build the vectors
c nullify the required places in vecf1
      ianf = 1
      do     lauf = istart,iende
        call zro1d(n,vecf1(ianf))
c NB!!! 
        if (cepai.and.j.gt.3*niv)then
c
           call zro1d(n,vecf2(ianf))
        end if
c NB!!!
        ianf = ianf + n
      enddo
      if(odebug) then
       write(iwr,*) 'vecf1 nullified upto ',ianf
      endif
      call rewftn(lun1)
c aufbau der eigentlichen vektoren
c die stehen dann auf vecf1
c build the real vectors 
c they are then stored on vecf1
c multiplikation der vorherigen entwicklungsvektoren
c mit den entwicklungskoeffizienten
c multiply the previous expansion vectors
c with the expansion coefficients
*      write(iwr,*) 'before do 190, j=',j
      do 190 i=1,j
         read (lun1) w
         ianf = 1
         if(odebug)
     +   write(iwr,*) 'in do 190, istart,iende',istart,iende
         do     lauf = istart,iende
           uim = u(i,iselec(lauf))
           call vsma1d(n,uim,w,vecf1(ianf))
           ianf = ianf + n
         enddo
  190 continue
         ianf = 1
c         do      lauf = istart,iende
c           write(iwr,*) 'new vector of ',lauf,'ianf =',ianf,
c     *'mrm = ',mrm(lauf),'mm ',mm(lauf)
c           call debvec(n,trawf,vecf1(ianf))
c           ianf = ianf + n
c         enddo

      call rewftn(lun)
      do 200 i = 1,istart - 1
         if (istart.eq.1) go to 277
*         write(iwr,*) 'read from',lun
         read (lun)
200   continue
277   continue
      ianf = 1
*      write(iwr,*) 'after do 200 istart,iende ',istart,iende
*      write(iwr,*) 'before writing new vectors, n, lun ',n,lun
      do     lauf = istart,iende
*        write(iwr,*) 'before writing, ianf=',ianf
*        write(iwr,*) vecf1(ianf),vecf1(ianf+n)
*        do   iraus = ianf,ianf+n
*          write(2,77) iraus,vecf1(iraus)
*        enddo
*77    format(3x,i7,3x,f15.10)
        call mwrite(n,lun,vecf1(ianf))
*        write(lun)(vecf1(iraus),iraus=ianf,ianf+n)
*        write(iwr,*) 'after writing out'
        ianf = ianf + n
      enddo
*      write(iwr,*) 'after writing the new vectors'
c
c aufbau -e*bi: vecf1 ist bi; vecf2 dann e*bi
      ianf = 1
      do    lauf = istart,iende
         s = -valn(lauf)
*        write(iwr,*) 'energy at residue calculation',s
c NB!!!
c
c calculates shift vector for 
c the cepa residue vector in vecf2
c
c 1. vecf2 = shift*bi => in vsmd()
c 2. nulify the reference coeficients => in vsmace()
c 3. vecf1 = -e*bi + vecf2 => in vsma1d()
c !! vsmd() is for ci
c
       if (crtce) then
            call vsmd(n,shift,vecf1(ianf),vecf2(ianf))
            call vsmace(nko,vecf2(ianf))
            call vsma1d(n,s,vecf1(ianf),vecf2(ianf))
         else     
            call vsmd(n,s,vecf1(ianf),vecf2(ianf))
         end if
c NB
         tradu = 0.5d0
*        write(iwr,*) '-e * bi',tradu
*        call debvec(n,tradu,vecf2(ianf))
         ianf = ianf + n
      enddo
c
c
c aufbau des residuen vektors (h*bi-e*bi)
c steht dann auf vecf2 ueberschreibt also -e*bi
c build the residue vectors (h*bi-e*bi)
c then stored on vecf2 thus overwrites -e*bi
      call rewftn(lun2)
      do 220 i=1,j
         read (lun2) v
*           write(iwr,*) 'vektor ',i
*           call debvec(n,tradu,v(1))
         ianf = 1
         do     lauf = istart,iende
*           write(iwr,*)'i,lauf,iselec(lauf)',i,lauf,iselec(lauf)
            uim = u(i,iselec(lauf))
c           write(iwr,*) 'uim',uim
*           write(iwr,*)'n,ianf',n,ianf
*           write(iwr,*) vecf2(ianf)
            call vsma1d(n,uim,v,vecf2(ianf))
c h*bi-e*bi steht jetzt auf vecf2
c h*bi-e*bi now on vecf2
            ianf = ianf + n
         enddo
  220 continue
      ianf = 1
      do    lauf = istart,iende
*        write(iwr,*) 'residual vector ',trakl
*        call debvec(n,tradu,vecf2(ianf))
         ianf = ianf + n
      enddo
c

c
      ianf = 1
      do    lauf = istart,iende
c       call vdotd(n,vecf2(ianf),vecf2(ianf),reser)
        reser = ddot(n,vecf2(ianf),1,vecf2(ianf),1)
        resnm(lauf) = sqrt(reser)
*       write(iwr,*) 'lauf,reser ',lauf,reser
        valoln(lauf) = valn(lauf)
        ianf = ianf + n
      enddo
c NB screening for max residue norm (resmx)
      call maxmgv(resnm,1,resmx,i,neig)
c NB
c     if (nout .gt. 0) write (nout,'('' residual norm='',d14.7)') resnor
c aufbau (h*bi - e*bi)/(e - haa); steht danach auf vecf2
c build (h*bi - e*bi)/(e - haa); here after stored on vecf2
      if( j .lt. mx ) then
c         write(iwr,*) 'new start vectors '
          ianf = 1
          do   lauf = istart,iende
c NB
c the preconditioner shift is implemented on the fly 
c of the Cdav calculation in vdssce()
c
         if (crtce) then
          call vdssce(n,nko,valn(lauf),shift,
     +         vecf2(ianf),adiag,vecf1(ianf))
          else
          call vdssvd(n,valn(lauf),vecf2(ianf),adiag,
     +         vecf1(ianf))
         end if
*        write(iwr,*) 'start vector ',lauf
         trara = 0.005d0
*        call debvec(n,trara,vecf1(ianf))
c NB get the tester
         call maxmgv(vecf1(ianf),1,tester(lauf),i,n)
c NB
         ianf = ianf + n
        enddo
      endif
c NB screen for maximum tester (testmx)
      call maxmgv(tester,1,testmx,i,neig)
c     write(iwr,*)'root number',neig,'max tester',testmx
c NB

c------
c check der einzelnen beteiligten vektoren wer und wenn warum 
c konvergiert
c check separate vectors to see who has converged and why
c-----

c konvergenzcheck ueber die energie
c convergence check on the energy

      ianf = 1
      do    lauf = istart,iende
         denerg = dabs(valn(lauf) - e(iselec(lauf)))
c         write(iwr,*) lauf,iselec(lauf),valn(lauf),e(iselec(lauf))
c NB cepa prints
         if (lauf.eq.istart) then
            if (.not.crtce) then
               write(nout,5) iter,iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
            else if (icepa.eq.1) then
               write(nout,15) iter,iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
            else if(icepa.eq.2) then 
               write(nout,25) iter,iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
            else if (icepa.eq.3) then
                write(nout,35) iter,iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
            end if
         else if (.not.crtce) then
              write(nout,7) iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
             else  if (icepa.eq.1) then 
                write(nout,17) iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
             else if (icepa.eq.2)  then
                 write(nout,27) iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
             else if (icepa.eq.3) then
                  write(nout,37) iselec(lauf),
     *valn2(lauf),denerg,tester(lauf),resnm(lauf)
         endif
c NB!!!
*         write(iwr,*) 'energy lowering for root ',iselec(lauf),
*     *denerg
*         write(iwr,*) 'total energy is  ',valn2(lauf)
chs      if ( denerg .lt. crite ) then
         if ((denerg .lt. crite).and.(iter.ge.4))then
            write(iwr,776)iselec(lauf)
776         format(1x,'root ',i2,' is converged')
c
c NB save the energy from the last iter in e()
c
            e(iselec(lauf))=valn2(lauf)
c NB
            iconv(iselec(lauf)) = 0
         endif
      enddo
      write(iwr,*)

      ico = 0
      do    lauf = istart,iende
          ico = ico + iconv(iselec(lauf))
      enddo
      if (ico .eq. 0) then
          write(iwr,67) 
67        format(10x,'=== all roots are converged')
c         write(iwr,*) 'on to new roots '
          go to 9001
      endif


*     write(iwr,*) 'j ist ',j
c rwah j.gt.1 --> j.gt.0 for single root
      if(j.gt.0) then
         smax2=big
c  60    smax1=0.0d0
         smax1=0.0d0
c-- mit allen alten
c-- with all old ones
         nneue = 0
         ialt = j
         ianf = 1
         do    lauf = istart,iende
c orthogonaliesung des betrachteten auf die alten
c orthogonalise all vectors under consideration on the old ones
*           write(iwr,*)'normalisation of vector',lauf,ianf
           call rewftn(lun1)
           do    jalt = 1,ialt
              read(lun1) w
c             call vdotd(n,w,vecf1(ianf),s)
              s = ddot(n,w,1,vecf1(ianf),1)
c              write(iwr,*)'ortho. onto old',jalt,s
c NB1
              s = -s
              call vsma1d(n,s,w,vecf1(ianf))
           enddo
c normierung dessen was vom betrachteten uebrig bleibt
c normalise the remainder of the vectors under consideration
c          call vdotd(n,vecf1(ianf),vecf1(ianf),s)
           s = ddot(n,vecf1(ianf),1,vecf1(ianf),1)
c           write(iwr,*) 'norm of new vector',lauf,' is',s
c NB
c still to be polished
c conv. in vectors = conv. in energy squared 
c i.e. (s=criten**2)
c original          if (s.gt.0.0000001)then
           if (s.gt.criten**2)then
c NB
c neuer wird mitgenommen
c new one is taken along
              s = 1.0d0 /sqrt(s)
              call vsmd(n,s,vecf1(ianf),vecf1(ianf))
              ialt = ialt + 1
              nneue = nneue + 1
              ineue(nneue) = ianf
c wegschreiben des neuen zu den alten
c write the new ones to the old
              call mwrite(n,lun1,vecf1(ianf))
           endif
           ianf = ianf + n
         enddo
c falls normierung weglaeuft noch programmieren
c incase the normalisation drifts still to program
      endif
c----   umspeichern der ueberlebenden vektoren auf vecf1
c----   restore the surviving vectors onto vecf1
      ianf = 1
      do    lauf = 1,nneue
c         call copyv(n,vecf1(ineue(lauf)),vecf1(ianf))
          if (ineue(lauf).ne.ianf)
     +        call dcopy(n,vecf1(ineue(lauf)),1,vecf1(ianf),1)
          ianf = ianf + n
      enddo
c check von orthonormalitaet
c check orthonormality
c     ianf = 1
c     do    idoof = 1,nneue
c         call vdotd(n,vecf1(ianf),vecf1(ianf),s)
c         s = ddot(n,vecf1(ianf),1,vecf1(ianf),1)
c          write(iwr,*) 'norm of vector ',idoof,'is ',s
c         ianf = ianf + n
c     enddo
*      write(iwr,*) 'n = ',n
c     ianf = 1
c     do    idoof = 1,nneue
*        write(iwr,*) 'vector ',idoof,'  ianf =',ianf
c       janf = 1
c       do  jdoof = 1,idoof
c         call vdotd(n,vecf1(janf),vecf1(ianf),s)
c         s = ddot(n,vecf1(janf),1,vecf1(ianf),1)
c          write(iwr,*) 'with vector ',jdoof,'  janf =',janf,
c    +  'liefert s=',s
c         janf = janf + n
c       enddo
c       ianf = ianf + n
c     enddo
c vecf1 enthaelt ineue orthonormalisierte vektoren
c vecf1 holds ineue orthonormal vectors
c
c --- buildup sigma-vector
csut  call mvmpy(j,n,v,w)
cbe      ioint = 0
c nullsetzen von vecf2
c nullify vecf2
      idoof = nneue * n
      do    jdoof = 1,idoof
         vecf2(jdoof) = 0.0d0
      enddo
      ianf = 1
*     do    lauf = 1,nneue
c eskiny ruft den skiny auf uebergibt aber
c nur den teil von vecf1,vecf2 der gebraucht
c wird
c eskiny call skiny but passes only that part of vecf1 and vecf2
c on that is really needed
*       write(iwr,*) 'before skiny:vecf2(',ianf,')',vecf2(ianf)
*       call  eskiny(n,ioint,vecf2(ianf),vecf1(ianf))
*       write(iwr,*) 'after skiny:vecf2(',ianf,')',vecf2(ianf)
        nanz = nneue
        if (issk.lt.iswhm) then
         call skiny(twoe,nteint,
     +              pey,acoul,aexc,ndeke,core,
     +              vecf2,vecf1,mdi,
     +              emat,nedim,iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,ioint,odebug)
        end if
        if (issk.gt.0) then
         call sham(vecf2, vecf1, mdi)
        end if
        call sdiag(vecf2, vecf1, mdi, hp5, mdi0)
*       ianf = ianf + n
*     enddo
c danach sollte alle neuen vektoren auf vecf2 stehen
c here after all new sigma-vectors should be on vecf2
c-- debug
c     ianf = 1
c     do    lauf = 1,nneue
c       write(iwr,*) 'new sigma vector for ',lauf
c       print *,'trawf =',trawf
c       ianf = ianf + n
c     enddo
c abspeichern der sigma vektoren nach lun2
c store the sigma vectors on lun2
c
*     write(iwr,*) 'lun2 = ',lun2
      ianf = 1
      do    lauf = 1,nneue
         call mwrite(n,lun2,vecf2(ianf))
         ianf = ianf + n
      enddo
c berechnung der neuen matrixelemente der davidson-matrix
c calculate the new matrix elements of the davidson matrix
c      write(iwr,*) 'j is here ',j
c      write(iwr,*) 'new av-elements with the old ones'
c
      call rewftn(lun1)
c
c zunaechst mit den alten
c now with the old vectors
      do      jalt = 1,j
        read(lun1) v
        ianf = 1
        ilauf = iideks(j) + jalt
        do    jneu = 1,nneue
c         call vdotd(n,v,vecf2(ianf),s)
          s = ddot(n,v,1,vecf2(ianf),1)
          av(ilauf) = s
c          write(iwr,*)'av(',ilauf,') = ',av(ilauf)
c NB!!!
        if (cepai.and.jalt.ne.1) then
           call ddotce(n,nko,v,vecf1(ianf),ui0(ilauf))
           ui0(ilauf) = - ui0(ilauf)
        end if
c NB!!!
          ilauf = ilauf + jneu + j
          ianf = ianf + n
        enddo
      enddo
c dann mit den neuen
c then with the new vectors
*      write(iwr,*) 'new av-elements with the old ones'
      ianf = 1
      do    jneu = 1,nneue
        janf = 1
        ilauf = iideks(j+jneu) - jneu + 1
        do    jneu1 = 1,jneu
c         call vdotd(n,vecf2(ianf),vecf1(janf),s)
          s = ddot(n,vecf2(ianf),1,vecf1(janf),1)
          av(ilauf) = s
*          write(iwr,*)'av(',ilauf,') = ',av(ilauf)
c NB
           if (cepai) then
              call ddotce(n,nko,vecf1(janf),vecf1(ianf),ui0(ilauf))
              if (jneu1.eq.jneu) then
                 ui0(ilauf) = 1.0d0 - ui0(ilauf)
              else
                 ui0(ilauf) =   - ui0(ilauf)
              end if
           end if
c NB
          janf = janf + n
          ilauf = ilauf + 1
        enddo
        ianf = ianf + n
        lauf = lauf + j
      enddo

c

c-- merken wie die energien im jetzigen lauf waren
c-- store the energies of the past iteration
      do    lauf = istart,iende
        e(iselec(lauf)) = valn(lauf)
*       write(iwr,*)'noting the energy in the current pass'
*       write(iwr,*) e(iselec(lauf)),valn(lauf)
      enddo
c zur neuen iteration
c to the next iteration
c      write(iwr,*) 'transition to the next iteration'
      j = j + nneue
c      write(iwr,*) 'new j is ',j
c neue laenge der av matrix (als unteres dreieck)
c new length of the av matrix (lower triangle)
      jw = (j*j-j)/2
c      write(iwr,*) 'new jw is ',jw
      goto 50

c zu neuen wurzeln
c to new roots
9001  continue
c -- in wavefunc werden die iterierten wf umgespeichert
c -- in wavefunc the iterated wavefunctions get stored
c     call wavefunv
      istart = iende
*      write(iwr,*) 'to the next pass '
9000  continue
      write(iwr,*)
      write(iwr,66)
csut  stop
      return
c-----   restart bereich
c-----   restart domain

c----- format section
1     format(/1x,' izus-field : ',50i3)
c2    format(/1x,i3, 'roots are treated in ',i3,' cycles')
3     format(1x,' the ordered roots are ',50i3)
4     format(/1x,i3,'-th cycle '//)
6     format(/1x,'commence pass at ', f10.2, ' seconds',a10,' wall'/
     + 1x,69('=')/
     + 1x,'iter',2x,'root-no.',3x,'energy',7x,'e.-lower'
     +,7x,'tester',7x,'res. norm'/
     + 1x,69('=') )
66    format(1x,68('=')//1x,'**** all roots are converged')
5     format(1x,i3,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3)
7     format(1x,3x,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3)
c NB cepa formats
15     format(1x,i3,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'cepa0')
25     format(1x,i3,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'acpf')
35     format(1x,i3,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'aqcc')
17     format(1x,3x,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'cepa0')
27     format(1x,3x,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'acpf')
37     format(1x,3x,3x,i3,2x,f15.10,3x,d9.3,6x,d8.3,7x,d8.3,2x,'aqcc')
c NB!!!
8     format(1x,' this pass treats the following roots : ',50i3)
c----- error sections
 1001 write(nout,*) 'error(1): maximal number of eigenvectors : ',mx
      write(nout,*) '          niv is                         : ',niv
      call caserr('initial vectors are invalid')
      return
      end
c@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
      subroutine extp(nbak,oextrap)
c
c  subroutine setzt nbak auf 1
c  ueber nbak wird zweite diagonalisierung gesteuert
c  subroutine set nbak to 1
c  nbak is used to control the second diagonalisation
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      logical debugs, oconf, oextrap
      integer nbak
      common /parkin/egey1,trash,tdel,
     +              ical0,nele,nkoo,mxex,nmulp,ispacep,nprin,
     1              nft31,nstarv,ncorci,nform,oconf,debugs,
     +              lsng,maxci,ipt0,isym(mxroot)
c
c    do not extrapolate or perform calculation at higher
c    threshold if trash or tdel is zero
c
      if (trash.le.1.0d-12.or.tdel.le.1.0d-12) then
       nbak = 0
       oextrap = .false.
      else
       nbak=1
       oextrap = .true.
      endif
c
      return
      end
c------------------------------------------------------------------------
c@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
      subroutine foxy(twoe,nteint,ston,nidmax,
     +                 pey, acoul, aexc, ndeke, core,
     +                 pey0, acoul0, aexc0, noeint,
     +                 civec,civec0,mdi,ed,mdiv,
     +                 emat,nedim,iottr,iotnew,iot0,iotm,
     +                 maindf,maxr,jconb,nomax2,odebug)
c
c  dieses programm speicherte die matrix um for die
c  extrapolation.
c  sie sollte nun die konfigurationsliste aendern....
c  this program restored the matrix for the extrapolation
c  now she should change the configuration list...
c
cbe  einbau von mpos1, sichert die position der mains in
cbe  kleinem problem
cbe  building in mpos1, ensure the position of mains in case of 
cbe  a small calculation
cbe
c
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
c
      real*8 twoe
      dimension twoe(nteint)
      real*8 ston
      integer nteint, ndeke, nidmax
      real*8 pey,core,pey0
      real*8 acoul,aexc,acoul0,aexc0
      dimension ston(nidmax)
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
      dimension pey0(noeint),acoul0(noeint),aexc0(noeint)
c
      real*8 civec, civec0, ed
      integer mdi, mdiv
c --- diagonale (ed)
      dimension civec(mdi),civec0(mdi),ed(mdiv)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iot0, iotm
      dimension iottr(iotm), iotnew(iotm), iot0(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      logical odebug
c
      parameter (mx=256)
      integer nconf2
      dimension nconf2(iswhm)
      real*8 rhm
      real*8 bn1,bn2,bn3,citr,cito
c
      integer ndet,nsac,idrem1
      common /rhus/ idrem1,nsac(iswhm),ndet(iswhm)
c
      integer niott
      dimension niott(iswhm)
      integer ifrk,jsum,kml,ii1
      integer ig,itests
*     integer igl,icl
      integer pcf,pcft
      integer imark,ncount
      common/rb/ ifrk
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      real*8 trsum
      common /miscop/ trsum(10,mxroot)
c
      integer i,j,k
      common /rwork/ niot0(iswhm)
      common /cselec/ iselec(mx),izus(mx)
c --- common for ft31-information
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
csut /speicher
      integer icount,irech,imm,jmm
      real*8 h
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),h(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c --- common for extrapolation
      real*8 esav0,esav1,csum,de0,de1,ueb0,ueb1,egeys
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .               de0(mxroot),de1(mxroot),
     .              ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .              egeys(mxroot)
c --- input
      integer nroot,ifirst,ilast
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                ndeci, icode,konf ,keps,iggey,
     .                istart
c
      common /cmpos/ mipos(mx),mipos1(mx)
      common /c/ idim,nrotr
csut timer
      real*8 cpu1,cpu0,dcpu
cdebug
      integer nslat,nsaf,idim0
      dimension nslat(iswhm),nsaf(iswhm)
c
      luns  = lun20
      lun   = lun21
c
      do isk=1,iswhm
        nconf2(isk) = nconf(isk)
      end do
*     write(iwr,*)
*     write(iwr,*) 'f o x y n '
*     write(iwr,*) '=========='
*     write(iwr,*)
c --- rausschreiben der konfigurationen
c --- writing out the configurations
      do 9727 ii=1,iswhm
          nsaf(ii) = nconf(ii)*nsac(ii)
          nslat(ii) = nconf(ii)*ndet(ii)
9727  continue
      write(iwr,9728)
9728  format(/
     +1x,'==================='/
     +1x,'= Lower threshold ='/
     +1x,'==================='//1x,
     +'super category                1        2       3       4      5'
     +/1x,
     +'----------------------------------------------------------------'
     +)
      write(iwr,1111) 'configurations      : ',nconf
      write(iwr,1111) 'safs                : ',nsaf
      write(iwr,1111) 'slaterdeterminants  : ',nslat
      write(iwr,*)
c     write(6,*) 'niot',niot
      ntoco = 0
      ntosa = 0
      ntosl = 0
      do 202 i=1,iswhm
         ntoco = ntoco + nconf(i)
         ntosa = ntosa + nsaf(i)
         ntosl = ntosl + nslat(i)
 202  continue
c     write(iwr,*) 'total :  ',ntoco,' konfigurationen,   ',
c    .           ntosa,' safs  und ',ntosl,' slaterdeterminanten.'
      write(iwr,2283) ntoco, ntosa, ntosl
2283  format(1x,'total number of configurations    ;', i8/
     +       1x,'total number of safs              ;', i8/
     +       1x,'total number of slaterdeterminants;', i8/)
c     write(6,*)
c --- lies alte daten ein
c --- lies ft35
      call rewftn(nhuk)
      read (nhuk)
      jsum=0
      do i=1,iswhm
        nc = nconf(i)
        if (nc.ne.0) then
          read(nhuk) ii1,kml
          read(nhuk)
        end if
*       jsum = jsum+nc*kml
      end do
      read(nhuk) bn1,bn2,bn3,citr,cito
      jsum = (idim-1)/ifrk+1
      do i=1,jsum
        read(nhuk)
      end do
cbe trsum enthaelt die summe ueber die energieerniedrigungen
cbe for 10 unterschiedlich grosse raeume
cbe trsum holds the sum over the energy lowerings for 10
cbe spaces of different size
      read(nhuk) ((trsum(j,i),j=1,10),i=1,nrotr)
cedbugwrite(6,*) 'bn1,=bn2,bn3,citr,cito',
cedbug bn1,bn2,bn3,citr,cito
cdebugwrite(6,*)'((trsum(j,i),j=1,10),i=1,nrotr)'
c     write(6,*) 'nrotr=',nrotr
cdebug      do i=1,nrotr
cdebug       do ii=1,10
cdebug        write(6,*) 'trsum(j,i)=',trsum(ii,i),i,ii
cdebug       end do
cdebug      end do
c --- umkopieren auf 21
c --- copy onto 21
      call rewftn(lun)
      call rewftn(luns)
      do i=1,nroot
        read(luns) (civec(j),j=1,idim)
        write(lun) (civec(j),j=1,idim)
cbe summe der energieerniedrigungen for grosses problem
cbe abspeicherung wird schon im rumpxs vorgenommen   15.8.94
cbe sum of energy lowerings for large problem
cbe storage already done in rumpxs                   15.8.94
*        de0(i) = trsum(4,i)
cdebug  write(6,*) 'i,de0',i,de0(i)
cbe summe der energieerniedrigungen for kleines problem
cbe sum of energy lowerings for small problem
*        de1(i) = trsum(5,i)
cdebug  write(6,*) 'i,de1',i,de1(i)
      end do
c --- selektion
      citr = citr+cito
      call rewftn(mousa)
      iy = 1
      ib = ifrk
      itr = ((idim-1)/ifrk)+1
      do ii=1,itr
       ibb = ib
       if(ib.gt.mdi) ibb = mdi
       read(mousa) (civec(iq),iq=iy,ibb)
       iy=iy+ifrk
       ib=ib+ifrk
      end do
c --- loop uber die konfigurationen
c --- loop over configurations

cbe ipos ist laeufer zum aufbau von mipos1
cbe ipos is counter for building mipos1
      ipos = 0
      ig = 0
      do iws=1,idim
       rhm = civec(iws)
c --- kriterium
       if(rhm .ge. citr) then
        ig = ig + 1
cdel    civec(ig) = rhm
csut    iselsf(ig) = 1
c --- diagonale
        ed(ig) = ed(iws)
       end if

cbe aufbau von mipos1, sichert die positionen der mains
cbe im kleinen problem, benoetigt zum bilden der ueberlapps
cbe mit startvektoren
cbe building mipos1, ensures the position of mains
cbe in small problem, needed for building the overlaps
cbe with start vectors

       if(rhm .eq. 940.0d0) then
         ipos = ipos + 1
         mipos1(ipos) = ig
       endif

      end do

cbe debug print von mipos1
c     write(6,*) 'mipos1 debug print'
c     write(6,782) (mipos1(iaus),iaus=1,20)
c782   format(20i5)

      idim0=ig
      write(iwr,*) 'reduction of ',idim,' --->',idim0,' safs'
c --- loop ueber die wurzel
      call rewftn(lun)
      call rewftn(luns)
      call rewftn(lunalt)
c     write(iwr,*) 'out foxy.f ,nroot',nroot
c -hallo baustelle!!!
      if(nroot.ne.ilast) then
        nanf = 1
        do ilauf=1,nroot
          nende = iselec(ilauf)
          if (oprint(32)) then
           write(iwr,*) 'nende',nende
          endif
c - einsortieren aller startvektoren auf die nicht konvergiert wurde
c - sorting in all start vectors that did not converge
          do jlauf=nanf,nende-1
            if (oprint(32)) then
             write(iwr,*) 'old file number',jlauf,' of ',lunalt
            endif
             read(lunalt) (civec0(j),j=1,idim)
             ig = 0
             do iws=1,idim
               rhm = civec(iws)
c - kriterium
               if(rhm .ge. citr) then
                 ig = ig + 1
                 civec0(ig) = civec0(iws)
               end if
             end do
             write(luns) (civec0(j),j=1,idim0)
            if (oprint(32)) then
             write(iwr,*) 'start vector of same length as alt ',idim0
            endif
          enddo
c - weglesen des naechsten alten startfiles
c - read away the next old start file
          if (oprint(32)) then
          write(iwr,*) 'before reading, idim',idim,' lunalt',lunalt
          endif
          read(lunalt)
          if (oprint(32)) then
           write(iwr,*) 'after reading'
          endif
c - lesen des start vectors for eine zu ziehende wurzel
c - read the start vector for a root to be calculated
          if (oprint(32)) then
           write(iwr,*) 'reading the first right file'
          endif
          read(lun) (civec0(j),j=1,idim)
          ig = 0
          do iws=1,idim
             rhm = civec(iws)
c - kriterium
             if(rhm .ge. citr) then
               ig = ig + 1
               civec0(ig) = civec0(iws)
             end if
          end do
          write(luns) (civec0(j),j=1,idim0)
          if (oprint(32)) then
          write(iwr,*)'startvektor der laenge for neu ',idim0
          write(iwr,*)'start vector of the same length as nue ',idim0
          endif
          nanf = nende + 1
        enddo
      else
      do ir=1,nroot
        read(lun) (civec0(j),j=1,idim)
        ig = 0
        do iws=1,idim
         rhm = civec(iws)
c - kriterium
         if(rhm .ge. citr) then
           ig = ig + 1
           civec0(ig) = civec0(iws)
         end if
        end do
        write(luns) (civec0(j),j=1,idim0)
        if (oprint(32)) then
         write(iwr,*) 'length of starting vector  ',idim0
        endif
      end do
      endif
c - selektiere iot-feld
c - select iot-field
c - loop over supercategories
c ig : zaehler der safs
c ig : counts the safs
c pcf:zaehler der konfigurationen
c pcf:counts the configurations
c ipiot : position im iot-feld
c ipiot : position in iot-field
      ig = 0
      pcft= 0
      ipiot=1
      ipiota=1
      do isk=1,iswhm
        nc = nconf(isk)
        nsa = nsac(isk)
        pcf = 0
c nytp : anzahl iot-eintraege pro sk
c nytp : number of iot-entries per sk
        nytp = nytl(isk)
c --- loop ueber die konfigurationen
c --- loop over configurations
        do j=1,nc
          do k=1,nsa
            ig=ig+1
            rhm = civec(ig)
            if (rhm .ge. citr) then
c             selektiere
c             select
              itests= 1
            else
              itests=-1
            end if
          end do
          if (itests.gt.0) then
c           schiebe iot-feld
c           shift iot-field
            do ii=1,nytp
              iottr(ipiot+ii-1) =
     .        iot0(ipiota+ii-1)
            end do
csut        write(iwr,*) 'transferred'
csut        write(iwr,*) (iottr(ipiot+ii-1),ii=1,nytp)
            pcf  = pcf+1
            pcft = pcft+1
            ipiot= ipiot+nytp
          end if
          ipiota = ipiota+nytp
        end do
        nconf(isk) = pcf
        niott(isk) = pcf*nytp
      end do
      do isk=1,iswhm
        niot(isk) = 0
      end do
      do isk=2,iswhm
        niot(isk) = niott(isk-1)+niot(isk-1)
      end do
cdebugwrite(iwr,*) 'niot',niot
ccccccccccccccccccccccccccccccccccc
cvp komprimierung der konfigurationsliste mit erzeuger und vernichter
cvp compress the configuration list with creators and annihilators
      ig = 0
      imark=0
      ncount=0
      do 2010 isk=1,iswhm
        nc = nconf2(isk)
        nsa = nsac(isk)
        inumb=0
c nytp : anzahl iot-eintraege pro sk
c nytp : number of iot-entries per sk
        nytp = nytl(isk)
c --- loop ueber die konfigurationen
c --- loop over configurations
        do ii1=1,nc
          do k=1,nsa
            ig=ig+1
            rhm = civec(ig)
            if (rhm .ge. citr) then
c             selektiere
              itests= 1
            else
              itests=-1
            end if
          end do
          if (itests.gt.0) then
c           komprimiere iotnew
c           compress iotnew
cvp   write(iwr,*) ' nsac=,nconf2=,nconf= ',nsac,nconf2,nconf
c     do 2010 isk=1,iswhm
	    inumb=inumb+1
	    iot0(ncount+inumb             )=
     &         iotnew(iotnst(isk)+ii1              )
	    iot0(ncount+inumb+  nconf(isk))=
     &         iotnew(iotnst(isk)+ii1+  nconf2(isk))
	    iot0(ncount+inumb+2*nconf(isk))=
     &         iotnew(iotnst(isk)+ii1+2*nconf2(isk))
	    iot0(ncount+inumb+3*nconf(isk))=
     &         iotnew(iotnst(isk)+ii1+3*nconf2(isk))
	    iot0(ncount+inumb+4*nconf(isk))=
     &         iotnew(iotnst(isk)+ii1+4*nconf2(isk))
	  endif
        end do
c       write(iwr,*) ' isk=,inumb= ',isk,inumb
	ncount=ncount+5*nconf(isk)
 2010 continue
cvp
cvp  aufbau des neuen iotnst-feldes
cvp  build the new iotnst-fieldes
      iotnst(1)=0
      do 2050 ii1=2,iswhm
	iotnst(ii1)=iotnst(ii1-1)+5*nconf(ii1-1)
 2050 continue
c     write(iwr,*) ' iotnst= ',iotnst
cvp  iotnew wird mit den komprimierten konfigurationen belegt
cvp  iotnew is filled with compressed configurations
      do 2060 ii1=1,iotm
	iotnew(ii1)=iot0(ii1)
 2060 continue
c*************************************************************8
c
c     write(iwr,*) ' idim=,icount= ',idim,icount
cvp
c
c --- wiederaufbau des iottr-feldes
c --- rebuilding the iottr-field
      do 10200 i=1,iswhm
        istar=niot(i)
        do 10210 j=1,nytl(i)
          do 10220 k=1,nconf(i)
            iot0(istar+k+(j-1)*nconf(i))
     .    = iottr(istar+j+(k-1)*nytl(i))
10220     continue
10210   continue
10200 continue
      do 10211 i=1,iswhm
          nsaf(i) = nconf(i)*nsac(i)
          nslat(i)= nconf(i)*ndet(i)
10211 continue
      write(iwr,9729)
9729  format(/
     +1x,'===================='/
     +1x,'= Higher threshold ='/
     +1x,'===================='//1x,
     +'super category                1        2       3       4      5'
     +/1x,
     +'----------------------------------------------------------------'
     +)
      write(iwr,1111) '   configurations      :    ',nconf
      write(iwr,1111) '   safs                :    ',nsaf
      write(iwr,1111) '   slaterdeterminants  :    ',nslat
      write(iwr,*)
      do 10201 i=1,iotm
         iottr(i) = iot0(i)
10201 continue
      ntoco = 0
      ntosa = 0
      ntosl = 0
      do i=1,iswhm
         ntoco = ntoco + nconf(i)
         ntosa = ntosa + nsaf(i)
         ntosl = ntosl + nslat(i)
      enddo
      write(iwr,2283) ntoco, ntosa, ntosl
c     write(iwr,*)
c --- umkopieren der csf-zahl
      idim =idim0
cvp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cvp bestimmung und
cvp vorsortierung der integrale for ioint = 1
cvp  integraladressen werden auf iadl zwischengespeichert;
cvp  ist dieser vektor gefuellt, werden in isort die entsprechenden
cvp  integrale auf twoe abgespeichert; ist dieses feld gefuellt,
cvp  werden es auf mtype geschrieben inkl. der anzahl der integrale
cvp  und der nummer der entsprechenden konfiguration
cvp determination and presorting the integrals for ioint=1
cvp  integral addresses are stored intermediately on iadl
cvp  is that vector full, then isort stores the integrals on twoe,
cvp  is that field full, then they are written to mtype inclusive
cvp  the number of integrals and the number of the corresponding 
cvp  configuration
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ioint.eq.1) then
       if (issk.ne.0) then
        write(iwr,*) '1. call to preint'
        call rewftn(mtype)
        none = 1
        write(iwr,*) 'sk :',none,' bis sk: ',issk
        write(iwr,*) ' call to preint '
        call cputim(cpu0)
        call preint(twoe,nteint,ston,nidmax,
     +              pey, acoul, aexc, ndeke, core,
     +              pey0, acoul0, aexc0, noeint,
     +              iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              nrec31,none,issk)
        call cputim(cpu1)
        dcpu = (cpu1-cpu0)
        write(iwr,100) dcpu
100     format(1x,'time for preint = ',f10.2,' secs.')
       endif
      endif
      if (issk.gt.0) then
       call cputim(cpu0)
       call skinh(twoe,nteint,
     +            pey, acoul, aexc, ndeke, core,
     +            emat,nedim,
     +            iottr,iotnew,iotm,
     +            maindf,maxr,jconb,nomax2,ioint,odebug)
       call cputim(cpu1)
       dcpu = (cpu1-cpu0)
       write(iwr,*) 'time for skinh :',dcpu
c       --- anzahl matrixelemente
c       --- number of matrix elements
       none=icount+(irech-1)*ndims
       write(iwr,*) 'for ',none,' matrixelemente'
       write(iwr,*) 'on ',(irech-1),' records.'
      end if
      if (ioint.eq.1) then
       if (issk.lt.iswhm) then
        write(iwr,*) '2. call to preint'
        call rewftn(mtype)
        none = issk+1
        write(iwr,*) 'sk :',none,' to sk: ',iswhm
        call cputim(cpu0)
        call preint(twoe,nteint,ston,nidmax,
     +              pey, acoul, aexc, ndeke, core,
     +              pey0, acoul0, aexc0, noeint,
     +              iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              nrec31,none,iswhm)
        call cputim(cpu1)
        dcpu = (cpu1-cpu0)
        write(iwr,100) dcpu
       end if
      endif

cvp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
c --- formate
1111  format(a25,5i8)
      end
        subroutine header
c
c  --- print header
c
        implicit none
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
        write(iwr,*) ' ==============================='
        write(iwr,*) ' == Semi-direct MRD-CI Module =='
        write(iwr,*) ' ==============================='
        write(iwr,*)
        write(iwr,*) ' multiroot - version'
        write(iwr,*) ' creators/annihilators for configurations'
        write(iwr,*)
        return
        end
cvp
cvp berechnung der matrixelemente for dk=0
      subroutine hmdk0(twoe,pey,aexc,
     +                 vecf1,vecf2,mdi
     +                ,emat,nedim,i,nc,mc,kml,kmj,inopen,iz
     &                ,kz,nhuk,mtype,jto1,jg,ioint,iges,kend
     &                ,nspp1,nsppa,nsppb,nspp4)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      real*8 twoe
      dimension twoe(*)
      real*8 aexc,pey
      dimension aexc(*),pey(*)
c
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
cvp
cvp emat : darstellungsmatrizen
cvp emat : representation matrices
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  vecf1 and vecf2: vektors, possibly for multiple roots
c  ndim: dimension of the secular problem
c  nanz: anzahl der gleichzeitig behandelten wurzeln
c  nanz: number of roots to treat simultaneously
      common/rwork2/ ndim,nanz
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
cvp pass information from ft31 and ft33 to skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      real*8 tim1
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c jgcnt zaehlt die matrixelemente
c jgcnt counts the matrix elements
      jgcnt=0
cvp
csut  write(6,*) i,nc,mc,kml,kmj,inopen,iz,kz
      do 3013 l=1,nsppb+nspp4
 3013 iaddr(l)=iz+(ispiel(l)-1)*kml
cvp
c dk=0 p=1
cvp  nqr = # qr-faelle
cvp  nqr = # qr-cases
cvp  nql = # ql-faelle
cvp  nql = # ql-cases
cvp  irec = # eintraege von emat for ein saf-paar
cvp  irec = # entries of emat for one saf-couple
      nqr=ibinom(inopen,2)
      nql=nqr
      irec=6*ideks(nqr+1)
c dk=0 p=2
      nqr2=inopen
      irec2=ideks(nqr2+1)
cvp iema = anzahl der eintraege in emat for p=1
cvp iema = number of entries in emat for p=1
      iema=irec*kml*kml
c dk=0 p=3
      nqr3=inopen
      nql3=nqr3
      irec3=ideks(inopen+1)*2*inopen
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
cvp iemb = number of entries in emat for p=1 und p=2
      iemb= iema + irec2*kml*kml
cvp
      if (ioint.eq.1) then
c dk=0 p=1
      do 31 l=1,nspp1
cvp lesen der integrale
cvp read the integrals
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
 31   continue
      jto1=jto1+2*nspp1
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
        cbfeld(l)=twoe(jto1+l-nspp1)
   32 continue
      jto1=jto1+nsppa-nspp1
c dk=0 p=3
      do 133 l=nsppa+1,nsppb
        tim1=twoe(jto1+1+(l-nsppa-1)*inopen)
        sumint(l)=tim1
  133 continue
      endif
cvp
      if (ioint.eq.0) then
c dk=0 p=1
      do 41 l=1,nspp1
cvp lesen der integrale
cvp read the integrals
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
 41   continue
c p=2,r=1,2
      do 42 l=nspp1+1,nsppa
        cbfeld(l)=twoe(intcb(l))
   42 continue
c dk=0 p=3
      do 134 l=nsppa+1,nsppb
        tim1=sumint(l)
        sumint(l)=tim1
  134 continue
      endif
c dk=0 p=1
      do 131 l=1,nspp1
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
  131 continue
c dk=0 p=2
      do 132 l=nspp1+1,nsppa
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
  132 continue
c dk=0 p=3
      do 135 l=nsppa+1,nsppb
cvp  r-fall wird ueber vorzeichen von qr codiert
cvp  r-case is coded by the sign of qr
      if (nqrfal(l).lt.0) then
         nrcas(l)=2
         nqrfal(l)=-nqrfal(l)
       else
         nrcas(l)=1
      endif
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
c symmetrie der orbitale for einelektronen-integral
c symmetry of the orbitals of one-electron integral
      ntar=nrfal(l)
c 1.mo for einelektronen-integral (das groessere)
c 1st mo for one-electron integral (the higher one)
      nb=moafal(l)
c 2.mo for einelektronen-integral (das kleinere)
c 2nd mo for one electron integral (the lower one)
      mb=mobfal(l)
      nb=ideks(nb)+mb+ijone(ntar)
cvp einelektronenintegral wird separat behandelt
c einelektronenintegral wird addiert
cvp one-electron integral is treated separately
c one-electron integral is added
      sumint(l)=sumint(l)+pey(nb)
  135 continue
c dk=0 p=1
cvp berechnung der matrixelemente for p=1
cvp compute the matrix elements for p=1
      do 45 ih=1,kml
        do 45 m=1,kml
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec+(ih-1)*irec
      do 129 l=1,nspp1
          fakcb(l)=
     &     emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
     &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
     &             +(nqrfal(l)-nqlfal(l))*2+1 )
          fakex(l)=
     &     emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
     &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
     &             +(nqrfal(l)-nqlfal(l))*2+2 )
  129 continue
      do 130 l=1,nspp1
          tim= fakcb(l)*cbfeld(l) + fakex(l)*exfeld(l)
            qfeld(l)=tim
  130 continue
cvp
c dk=0 p=2
ct    do 145 ih=1,kml
ct      do 145 m=1,kml
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec2+(ih-1)*irec2 + iema
      do 142 l=nspp1+1,nsppa
          tim=
     &     emat(iem + (nqlfal(l)-1)*(nqr2+1)-ideks(nqlfal(l))
     &              + nqrfal(l)-nqlfal(l)+1 ) * cbfeld(l)
            qfeld(l)=tim
  142 continue
cvp
c dk=0 p=3
ct    do 245 ih=1,kml
ct      do 245 m=1,kml
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec3+(ih-1)*irec3 + iemb
      do 236 l=nsppa+1,nsppb
          hp3(l)=
     &     emat(iem + (nrcas(l)-1)*ideks(nqr3+1)
     &              + (nqlfal(l)-1)*(nqr3+1)-ideks(nqlfal(l)+1)
     &              + nqrfal(l)+1 ) * sumint(l)
  236 continue
cvp   summation ueber die gemeinsamen offenen schalen
cvp   summation over the common open shells
      do 345 iis=1,inopen-1
c lesen der negativen austauschintegrale
c read the negative exchange integrals
      if (ioint.eq.1) then
      do 186 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(jto1+iis+1+(l-nsppa-1)*inopen)
  186 continue
      endif
      if (ioint.eq.0) then
      do 187 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(nwwmo(l,iis))
  187 continue
      endif
      do 136 l=nsppa+1,nsppb
          hp3(l)= hp3(l)
     &    +emat(iem + iis*2*ideks(nqr3+1)
     &              + (nrcas(l)-1)*ideks(nqr3+1)
     &              + (nqlfal(l)-1)*(nqr3+1)-ideks(nqlfal(l)+1)
     &              + nqrfal(l)+1 ) * cbfeld(l)
  136 continue
  345 continue
      do 188 l=nsppa+1,nsppb
            qfeld(l)=hp3(l)
  188 continue
c
c schleife ueber die wurzeln
c loop over roots
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
c   ipoint: start position of the vector within vecf for the root iwurz
        ipoint=(iwurz-1)*ndim
cdebug  write(6,*) 'hmdk0 : ipoint= ',ipoint
cvp gather-schritt
cvp gather-step
        do 401 l=1,nspp1
            vgath1(l)=vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
            vgath2(l)=vecf2(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
  401   continue
        do 402 l=nspp1+1,nsppa
            vgath1(l)=vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
            vgath2(l)=vecf2(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
  402   continue
        do 403 l=nsppa+1,nsppb
            vgath1(l)=vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
            vgath2(l)=vecf2(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)
  403   continue
c   ende gather
c   end gather
         do 4212 l=1,nsppb
           if (icas(l).eq.0) then
             qq=qfeld(l)
             vecf1(ipoint+kz+m)=vecf1(ipoint+kz+m)+vgath2(l)*qq
           endif
 4212    continue
         do 4214 l=1,nsppb
           if (icas(l).eq.1) then
             qq=qfeld(l)
             vecf1(ipoint+kz+ih)=vecf1(ipoint+kz+ih)+vgath2(l)*qq
           endif
 4214    continue
cvocl loop,novrec
         do 4213 l=1,nsppb
           qq=qfeld(l)
           vgath1(l)=vgath1(l)
     &              +vecf2(ipoint+kz+(1-icas(l))*m+icas(l)*ih)*qq
 4213    continue
cvp
cvp scatter-schritt
cvp scatter-step
        do 411 l=1,nspp1
            vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)=vgath1(l)
  411   continue
        do 412 l=nspp1+1,nsppa
            vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)=vgath1(l)
  412   continue
        do 413 l=nsppa+1,nsppb
            vecf1(ipoint+iaddr(l)+(1-icas(l))*ih+icas(l)*m)=vgath1(l)
  413   continue
c ende schleife uber die wurzeln
c end loop over the roots
 5100 continue
c
   45 continue
cvp
c dk=0 p=4     aa/bb
      do 141 m=1,kml
      mz=kz+m
      do 139 l=nsppb+1,nsppb+nspp4
c missmatch
cvp absolutnummer von mo a
cvp absolute number of mo a
      kk=mobfal(l)
c missmatch in r
cvp absolutnummer von mo b
cvp absolute number of mo b
      ll=moafal(l)
cvp nrfal enthaelt hier die integraladresse
cvp nrfal here holds the address of the integral
      nrfal(l)=ideks(max(kk,ll))+min(kk,ll)
      tim=aexc(nrfal(l))
            qfeld(l)=tim
  139 continue
c
c schleife ueber die wurzeln
c loop over the roots
c
      do 5200 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
cvp gather-schritt
        do 404 l=nsppb+1,nsppb+nspp4
            vgath1(l)=vecf1(ipoint+iaddr(l)+m)
            vgath2(l)=vecf2(ipoint+iaddr(l)+m)
  404   continue
cvp
         do 4312 l=nsppb+1,nsppb+nspp4
           qq=qfeld(l)
           vecf1(ipoint+mz)=vecf1(ipoint+mz)+vgath2(l)*qq
 4312    continue
         do 4313 l=nsppb+1,nsppb+nspp4
           qq=qfeld(l)
           vgath1(l)=vgath1(l)+vecf2(ipoint+mz)*qq
 4313    continue
cvp scatter-schritt
        if (i.gt.issk) then
        do 414 l=nsppb+1,nsppb+nspp4
           vecf1(ipoint+iaddr(l)+m)=vgath1(l)
  414   continue
        end if
c ende schleife uber die wurzeln
 5200 continue
cvp
  141 continue
cvp
cvp
cvp hochsetzen vom jto1: wird nur for ioint=1 benoetigt
cvp increase jto1: needed only for ioint=1
      jto1=jto1+(nsppb-nsppa)*inopen
cvp
      return
      end
cvp
cvp berechnung der matrixelemente for dk=0
cvp calculate the matrix elements for dk=0
      subroutine hmdk0h(twoe,pey,aexc,
     +                  emat,nedim,i,nc,mc,kml,kmj,inopen,iz
     &                ,kz,jto1,jg,ioint,iges,kend
     &                ,nspp1,nsppa,nsppb,nspp4)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      real*8 twoe
      dimension twoe(*)
      real*8 pey, aexc
      dimension pey(*),aexc(*)
cvp
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
cvp emat : darstellungsmatrizen
cvp emat : representation matrices
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
cvp pass information from ft31 and ft33 to skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      real*8 tim1
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
c jgcnt zaehlt die matrixelemente
c jgcnt counts the matrix elements
      jgcnt=0
cvp
      do 3013 l=1,nsppb+nspp4
 3013 iaddr(l)=iz+(ispiel(l)-1)*kml
cvp
c dk=0 p=1
cvp  nqr = # qr-faelle
cvp  nqr = # qr-cases
cvp  nql = # ql-faelle
cvp  nql = # ql-cases
cvp  irec = # eintraege von emat for ein saf-paar
cvp  irec = # entries of emat for one saf-couple
      nqr=ibinom(inopen,2)
      nql=nqr
      irec=6*ideks(nqr+1)
c dk=0 p=2
      nqr2=inopen
      irec2=ideks(nqr2+1)
cvp iema = anzahl der eintraege in emat for p=1
cvp iema = number of entries in emat for p=1
      iema=irec*kml*kml
c dk=0 p=3
      nqr3=inopen
      nql3=nqr3
      irec3=ideks(inopen+1)*2*inopen
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
cvp iemb = number of entries in emat for p=1 and p=2
      iemb= iema + irec2*kml*kml
cvp
      if (ioint.eq.1) then
c dk=0 p=1
      do 31 l=1,nspp1
cvp lesen der integrale
cvp read the integrals
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
 31   continue
      jto1=jto1+2*nspp1
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
        cbfeld(l)=twoe(jto1+l-nspp1)
   32 continue
      jto1=jto1+nsppa-nspp1
c dk=0 p=3
      do 133 l=nsppa+1,nsppb
        tim1=twoe(jto1+1+(l-nsppa-1)*inopen)
        sumint(l)=tim1
  133 continue
      endif
cvp
      if (ioint.eq.0) then
c dk=0 p=1
      do 41 l=1,nspp1
cvp lesen der integrale
cvp read the integrals
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
 41   continue
c p=2,r=1,2
      do 42 l=nspp1+1,nsppa
        cbfeld(l)=twoe(intcb(l))
   42 continue
c dk=0 p=3
      do 134 l=nsppa+1,nsppb
        tim1=sumint(l)
        sumint(l)=tim1
  134 continue
      endif
c dk=0 p=1
      do 131 l=1,nspp1
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
  131 continue
c dk=0 p=2
      do 132 l=nspp1+1,nsppa
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
  132 continue
c dk=0 p=3
      do 135 l=nsppa+1,nsppb
cvp  r-fall wird ueber vorzeichen von qr codiert
cvp  r-case is coded by the sign of qr
      if (nqrfal(l).lt.0) then
         nrcas(l)=2
         nqrfal(l)=-nqrfal(l)
       else
         nrcas(l)=1
      endif
      if (nqrfal(l).ge.nqlfal(l)) then
         icas(l)=0
        else
         nqmax=nqlfal(l)
         nqlfal(l)=nqrfal(l)
         nqrfal(l)=nqmax
         icas(l)=1
      endif
c symmetrie der orbitale for einelektronen-integral
c symmetry of orbitals for one-electron integral
      ntar=nrfal(l)
c 1.mo for einelektronen-integral (das groessere)
c 1st mo for one-electron integral (the larger)
      nb=moafal(l)
c 2.mo for einelektronen-integral (das kleinere)
c 2nd mo for one-electron integral (the smaller)
      mb=mobfal(l)
      nb=ideks(nb)+mb+ijone(ntar)
cvp einelektronenintegral wird separat behandelt
c einelektronenintegral wird addiert
cvp one-electron integral is treated separately
c one-electron integral is added
      sumint(l)=sumint(l)+pey(nb)
  135 continue
c dk=0 p=1
cvp berechnung der matrixelemente for p=1
cvp calculation of matrix elements for p=1
      do 45 ih=1,kml
        do 45 m=1,kml
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec+(ih-1)*irec
      do 129 l=1,nspp1
          fakcb(l)=
     &     emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
     &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
     &             +(nqrfal(l)-nqlfal(l))*2+1 )
          fakex(l)=
     &     emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
     &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
     &             +(nqrfal(l)-nqlfal(l))*2+2 )
  129 continue
      do 130 l=1,nspp1
          tim= fakcb(l)*cbfeld(l) + fakex(l)*exfeld(l)
cx        tim=
cx   &     emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
cx   &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
cx   &             +(nqrfal(l)-nqlfal(l))*2+1 ) * cbfeld(l)
cx   &    +emat(iem+(nrfal(l)-1)*ideks(nqr+1)*2
cx   &             +(nqlfal(l)-1)*(nqr+1)*2 - 2*ideks(nqlfal(l))
cx   &             +(nqrfal(l)-nqlfal(l))*2+2 ) * exfeld(l)
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and 
cvp  reset jgcnt to zero
          if (dabs(tim).ge.1.0e-7) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=kz+(1-icas(l))*m+icas(l)*ih
            mjfeld(jgcnt)=iaddr(l)+(1-icas(l))*ih+icas(l)*m
          endif
  130 continue
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
ct    if (jgcnt.ge.ncmax) then
ct       do 4012 ik=1,jgcnt
ct         jg=jg+1
ct         q(jg)=qfeld(ik)
ct         mi(jg)=mifeld(ik)
ct         mj(jg)=mjfeld(ik)
ct         if (jg.ge.iwod) then
ct           jg=0
ct           write(nhuk) q,mi,mj
ct         endif
c4012    continue
ct       jgcnt=0
ct    endif
cvp
ct 45 continue
cvp
c dk=0 p=2
ct    do 145 ih=1,kml
ct      do 145 m=1,kml
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec2+(ih-1)*irec2 + iema
      do 142 l=nspp1+1,nsppa
          tim=
     &     emat(iem + (nqlfal(l)-1)*(nqr2+1)-ideks(nqlfal(l))
     &              + nqrfal(l)-nqlfal(l)+1 ) * cbfeld(l)
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and 
cvp  reset jgcnt to zero
          if (dabs(tim).ge.1.0e-7) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=kz+(1-icas(l))*m+icas(l)*ih
            mjfeld(jgcnt)=iaddr(l)+(1-icas(l))*ih+icas(l)*m
          endif
  142 continue
cvp
cvp
c dk=0 p=3
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kml*irec3+(ih-1)*irec3 + iemb
      do 236 l=nsppa+1,nsppb
          hp3(l)=
     &     emat(iem + (nrcas(l)-1)*ideks(nqr3+1)
     &              + (nqlfal(l)-1)*(nqr3+1)-ideks(nqlfal(l)+1)
     &              + nqrfal(l)+1 ) * sumint(l)
  236 continue
cvp   summation ueber die gemeinsamen offenen schalen
cvp   summation over the common open shells
      do 345 iis=1,inopen-1
c lesen der negativen austauschintegrale
c read the negative exchange integrals
      if (ioint.eq.1) then
      do 186 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(jto1+iis+1+(l-nsppa-1)*inopen)
  186 continue
      endif
      if (ioint.eq.0) then
      do 187 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(nwwmo(l,iis))
  187 continue
      endif
      do 136 l=nsppa+1,nsppb
          hp3(l)= hp3(l)
     &    +emat(iem + iis*2*ideks(nqr3+1)
     &              + (nrcas(l)-1)*ideks(nqr3+1)
     &              + (nqlfal(l)-1)*(nqr3+1)-ideks(nqlfal(l)+1)
     &              + nqrfal(l)+1 ) * cbfeld(l)
  136 continue
  345 continue
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and 
cvp  reset jgcnt to zero
      do 188 l=nsppa+1,nsppb
          if (dabs(hp3(l)).ge.hstor) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=hp3(l)
            mifeld(jgcnt)=kz+(1-icas(l))*m+icas(l)*ih
            mjfeld(jgcnt)=iaddr(l)+(1-icas(l))*ih+icas(l)*m
          endif
  188 continue
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
      if (jgcnt.ge.ncmax) then
         do 4212 ik=1,jgcnt
           icount=icount+1
           rham(icount) = (qfeld(ik))
           imm(icount)  = mifeld(ik)
           jmm(icount)  = mjfeld(ik)
           if (icount.ge.ndims) then
             icount=0
             irech=irech+1
             write(nf88) imm,jmm,rham
             do iii=1,ndims
                imm(iii) = 0
                jmm(iii) = 0
                rham(iii) = 0.0d0
             end do
           endif
 4212    continue
         jgcnt=0
      endif
cvp
ct245 continue
   45 continue
cvp
c dk=0 p=4     aa/bb
      do 141 m=1,kml
cvp
      do 139 l=nsppb+1,nsppb+nspp4
c missmatch
cvp absolutnummer von mo a
cvp absolute number of mo a
      kk=mobfal(l)
c missmatch in r
cvp absolutnummer von mo b
cvp absolute number of mo b
      ll=moafal(l)
cvp nrfal enthaelt hier die integraladresse
cvp nrfal here holds the integral address
      nrfal(l)=ideks(max(kk,ll))+min(kk,ll)
      tim=aexc(nrfal(l))
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and
cvp  reset jgcnt to zero
          if (dabs(tim).ge.hstor) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=kz+m
            mjfeld(jgcnt)=iaddr(l)+m
          endif
  139 continue
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
      if (jgcnt.ge.ncmax) then
         do 4312 ik=1,jgcnt
           icount=icount+1
           rham(icount) = (qfeld(ik))
           imm(icount)  = mifeld(ik)
           jmm(icount)  = mjfeld(ik)
           if (icount.ge.ndims) then
             icount=0
             irech=irech+1
             write(nf88) imm,jmm,rham
             do iii=1,ndims
               imm(iii) = 0
               jmm(iii) = 0
               rham(iii) = 0.0d0
             end do
           endif
 4312    continue
         jgcnt=0
      endif
cvp
  141 continue
cvp
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
      do 4013 ik=1,jgcnt
        icount=icount+1
        rham(icount) = (qfeld(ik))
        imm(icount)  = mifeld(ik)
        jmm(icount)  = mjfeld(ik)
        if (icount.ge.ndims) then
          icount=0
          irech=irech+1
          write(nf88) imm,jmm,rham
          do iii=1,ndims
             imm(iii) = 0
             jmm(iii) = 0
             rham(iii) = 0.0d0
          end do
        endif
 4013 continue
cvp
cvp hochsetzen vom jto1: wird nur for ioint=1 benoetigt
cvp increase jto1: needed only for ioint=1
      jto1=jto1+(nsppb-nsppa)*inopen
cvp
      return
      end
cvp
cvp berechnung der matrixelemente for dk=1
cvp compute the matrix elements for p=1
      subroutine hmdk1(twoe,pey,
     +                 vecf1,vecf2,mdi,
     +                 emat,nedim,i,nc,mc,kml,kmj,inopen,jz
     &                ,kz,nhuk,mtype,jto1,jg,ioint,iges,kend
     &                ,nspp1,nsppa,nsppb)
c    &                ,vecf1,vecf2)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      real*8 twoe
      real*8 pey
      dimension twoe(*),pey(*)
c
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
cvp emat : darstellungsmatrizen
cvp emat : representation matrices
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  vecf1 and vecf2: vectors, possibly. for multiple roots
c  ndim: dimension des saekularproblems
c  ndim: dimension of secular problem
c  nanz: anzahl der gleichzeitig behandelten wurzeln
c  nanz: number of roots to be treated simultaneously
      common/rwork2/ ndim,nanz
c --- issk : superkategorien auf disk
      integer issk
      real*8   hstor
      common /csemi/ hstor,issk
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
cvp pass information from ft31 and ft33 to skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      real*8 tim1
cvp
      do 3013 l=1,nsppb
 3013 iaddr(l)=jz+(ispiel(l)-1)*kmj
cvp
cvp dk=1 p=1
cvp  nqr = # qr-faelle
cvp  nql = # ql-faelle
cvp  irec = # eintraege von emat for ein saf-paar
cvp  nqr = # qr-cases
cvp  nql = # ql-cases
cvp  irec = # entries of emat for a saf-couple
      nqr=ibinom(inopen,3)
      nql=inopen-2
      irec=6*nqr*nql
cvp dk=1 p=2
      nqr2=ibinom(inopen,2)
      irec2=nqr2
cvp iema = anzahl der eintraege in emat for p=1
cvp iema = number of entries in emat for p=1
      iema=irec*kml*kmj
cvp dk=1 p=3
      nqr3=ibinom(inopen,2)
      irec3=(inopen-1)*2*nqr3
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
cvp iemb = number of entries in emat for p=1 and p=2
      iemb= iema + irec2*kml*kmj
cvp
c dk=1 p=1
      if (ioint.eq.1) then
cvp
      do 31 l=1,nspp1
cvp lesen der integrale
cvp read integrals
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
c which r-case ll=1-3 greater zero,ll=4-6 smaller null
        if (nrfal(l).lt.0) then
           nrfal(l)=-nrfal(l)
           cbfeld(l)= -cbfeld(l)
           exfeld(l)= -exfeld(l)
        endif
 31   continue
      jto1=jto1+2*nspp1
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
        cbfeld(l)=twoe(jto1+l-nspp1)
c qr vorzeichen kodiert r fall
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           cbfeld(l)=-cbfeld(l)
           nqrfal(l)=-nqrfal(l)
        endif
   32 continue
      jto1=jto1+nsppa-nspp1
c dk=1,p=3,r=1,2
c lesen der summe ueber die dreizentrenintegrale
c read the sum over the 3-centre integrals
      do 13 l=nsppa+1,nsppb
        tim1=twoe(jto1+1+(l-nsppa-1)*(inopen-1))
        sumint(l)=tim1
   13 continue
cvp
      endif
cvp
      if (ioint.eq.0) then
cvp
c dk=1 p=1
      do 41 l=1,nspp1
cvp lesen der integrale
cvp read integrals
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
c which r-case ll=1-3 greater zero,ll=4-6 smaller null
        if (nrfal(l).lt.0) then
           nrfal(l)=-nrfal(l)
           cbfeld(l)= -cbfeld(l)
           exfeld(l)= -exfeld(l)
        endif
 41   continue
c p=2,r=1,2
      do 42 l=nspp1+1,nsppa
        cbfeld(l)=twoe(intcb(l))
c qr vorzeichen kodiert r fall
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           cbfeld(l)=-cbfeld(l)
           nqrfal(l)=-nqrfal(l)
        endif
   42 continue
c dk=1 p=3
      do 23 l=nsppa+1,nsppb
        tim1=sumint(l)
        sumint(l)=tim1
   23 continue
cvp
      endif
c dk=1 p=3
      do 33 l=nsppa+1,nsppb
cvp  r-fall wird ueber vorzeichen von qr codiert
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           nrcas(l)=1
           nqrfal(l)=-nqrfal(l)
         else
           nrcas(l)=2
        endif
        ntar=nrfal(l)
        nb=moafal(l)
        mb=mobfal(l)
        nb=ideks(nb)+mb+ijone(ntar)
        sumint(l)=sumint(l)+pey(nb)
   33 continue
cvp
c jgcnt zaehlt die matrixelemente
c jgcnt counts matrix elements
      jgcnt=0
c dk=1 p=1
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
c begin end transformation of dk=1 (the same for all p-case)
      do 45 ih=1,kmj
        do 46 m=1,kml
        mz=kz+m
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kmj*irec+(ih-1)*irec
cvp
      do 131 l=1,nspp1
          tim=
     &     emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
     &             +nqrfal(l)*2-1 ) * cbfeld(l)
     &    +emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
     &             +nqrfal(l)*2   ) * exfeld(l)
            qfeld(l)=tim
  131 continue
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kmj*irec2+(ih-1)*irec2 + iema
cvp
      do 132 l=nspp1+1,nsppa
          tim=
     &     emat(iem
     &             +nqrfal(l) ) * cbfeld(l)
            qfeld(l)=tim
  132 continue
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kmj*irec3+(ih-1)*irec3 + iemb
cvp
      do 43 l=nsppa+1,nsppb
          hp3(l)=
     &     emat(iem + (nrcas(l)-1)*nqr3 + nqrfal(l) ) * sumint(l)
  43  continue
cvp
      do 345 iis=1,inopen-2
c lesen der (negativen) austauschintegrale
c read the (negative) exchange integrals
      if (ioint.eq.1) then
      do 66 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(jto1+iis+1+(l-nsppa-1)*(inopen-1))
   66 continue
      endif
      if (ioint.eq.0) then
      do 67 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(nwwmo(l,iis))
   67 continue
      endif
cvp   summation ueber die gemeinsamen offenen schalen
cvp   summation over common open shells
      do 68 l=nsppa+1,nsppb
          hp3(l)= hp3(l)
     &    +emat(iem + iis*2*nqr3
     &              + (nrcas(l)-1)*nqr3 + nqrfal(l) ) * cbfeld(l)
   68 continue
  345 continue
      do 69 l=nsppa+1,nsppb
            qfeld(l)=hp3(l)
   69 continue
c
c schleife ueber die wurzeln
c loop over roots
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
cvp gather-schritt
        do 401 l=1,nspp1
            vgath1(l)=vecf1(ipoint+iaddr(l)+ih)
            vgath2(l)=vecf2(ipoint+iaddr(l)+ih)
  401   continue
        do 402 l=nspp1+1,nsppa
            vgath1(l)=vecf1(ipoint+iaddr(l)+ih)
            vgath2(l)=vecf2(ipoint+iaddr(l)+ih)
  402   continue
        do 403 l=nsppa+1,nsppb
            vgath1(l)=vecf1(ipoint+iaddr(l)+ih)
            vgath2(l)=vecf2(ipoint+iaddr(l)+ih)
  403   continue
c  ende gather
         do 4212 l=1,nsppb
           qq=qfeld(l)
           vecf1(ipoint+mz)=vecf1(ipoint+mz)+vgath2(l)*qq
 4212    continue
cvocl loop,novrec
         do 4213 l=1,nsppb
           qq=qfeld(l)
           vgath1(l)=vgath1(l)+vecf2(ipoint+mz)*qq
 4213    continue
cvp
        if (i.gt.issk) then
cvp scatter-schritt
        do 411 l=1,nspp1
            vecf1(ipoint+iaddr(l)+ih)=vgath1(l)
csut    write(6,*) 'hmdk1:',vecf1(ipoint+iaddr(l)+ih)
csut    write(6,*) '======'
  411   continue
        do 412 l=nspp1+1,nsppa
            vecf1(ipoint+iaddr(l)+ih)=vgath1(l)
csut    write(6,*) 'hmdk1:',vecf1(ipoint+iaddr(l)+ih)
csut    write(6,*) '======'
  412   continue
        do 413 l=nsppa+1,nsppb
            vecf1(ipoint+iaddr(l)+ih)=vgath1(l)
  413   continue
        end if
c ende schleife uber die wurzeln
 5100 continue
c
   46 continue
   45 continue
cvp hochsetzen vom jto1: wird nur for ioint=1 benoetigt
cvp increase jto1: needed only for ioint=1
      jto1=jto1+(nsppb-nsppa)*(inopen-1)
cvp
      return
      end
cvp berechnung der matrixelemente for dk=1
cvp calculate the matrix elements for dk=1
      subroutine hmdk1h(twoe,pey,
     +       emat,nedim,i,nc,mc,kml,kmj,inopen,jz
     &      ,kz,jto1,jg,ioint,iges,kend
     &      ,nspp1,nsppa,nsppb)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      real*8 pey
      dimension twoe(*),pey(*)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
cvp
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
cvp emat : darstellungsmatrizen
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      real*8 tim1
cvp
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 3013 l=1,nsppb
 3013 iaddr(l)=jz+(ispiel(l)-1)*kmj
cvp
cvp dk=1 p=1
cvp  nqr = # qr-faelle
cvp  nql = # ql-faelle
cvp  irec = # eintraege von emat for ein saf-paar
cvp  nqr = # qr-cases
cvp  nql = # ql-cases
cvp  irec = # entries of emat for one saf-couple
      nqr=ibinom(inopen,3)
      nql=inopen-2
      irec=6*nqr*nql
cvp dk=1 p=2
      nqr2=ibinom(inopen,2)
      irec2=nqr2
cvp iema = anzahl der eintraege in emat for p=1
cvp iema = number of entries in emat for p=1
      iema=irec*kml*kmj
cvp dk=1 p=3
      nqr3=ibinom(inopen,2)
      irec3=(inopen-1)*2*nqr3
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
cvp iemb = number of entries in emat for p=1 und p=2
      iemb= iema + irec2*kml*kmj
cvp
c dk=1 p=1
      if (ioint.eq.1) then
cvp
      do 31 l=1,nspp1
cvp lesen der integrale
cvp read integrals
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
c which r-case ll=1-3 greater zero,ll=4-6 smaller null
        if (nrfal(l).lt.0) then
           nrfal(l)=-nrfal(l)
           cbfeld(l)= -cbfeld(l)
           exfeld(l)= -exfeld(l)
        endif
 31   continue
      jto1=jto1+2*nspp1
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
        cbfeld(l)=twoe(jto1+l-nspp1)
c qr vorzeichen kodiert r fall
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           cbfeld(l)=-cbfeld(l)
           nqrfal(l)=-nqrfal(l)
        endif
   32 continue
      jto1=jto1+nsppa-nspp1
c dk=1,p=3,r=1,2
c lesen der summe ueber die dreizentrenintegrale
c read the sum over the 3-centre integrals
      do 13 l=nsppa+1,nsppb
        tim1=twoe(jto1+1+(l-nsppa-1)*(inopen-1))
        sumint(l)=tim1
   13 continue
cvp
      endif
cvp
      if (ioint.eq.0) then
cvp
c dk=1 p=1
      do 41 l=1,nspp1
cvp lesen der integrale
cvp read integrals
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
        if (nrfal(l).lt.0) then
           nrfal(l)=-nrfal(l)
           cbfeld(l)= -cbfeld(l)
           exfeld(l)= -exfeld(l)
        endif
 41   continue
c p=2,r=1,2
      do 42 l=nspp1+1,nsppa
        cbfeld(l)=twoe(intcb(l))
c qr vorzeichen kodiert r fall
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           cbfeld(l)=-cbfeld(l)
           nqrfal(l)=-nqrfal(l)
        endif
   42 continue
c dk=1 p=3
      do 23 l=nsppa+1,nsppb
        tim1=sumint(l)
        sumint(l)=tim1
   23 continue
cvp
      endif
c dk=1 p=3
      do 33 l=nsppa+1,nsppb
cvp  r-fall wird ueber vorzeichen von qr codiert
c qr sign codes r case
        if (nqrfal(l).lt.0) then
           nrcas(l)=1
           nqrfal(l)=-nqrfal(l)
         else
           nrcas(l)=2
        endif
        ntar=nrfal(l)
        nb=moafal(l)
        mb=mobfal(l)
        nb=ideks(nb)+mb+ijone(ntar)
        sumint(l)=sumint(l)+pey(nb)
   33 continue
cvp
c jgcnt zaehlt die matrixelemente
c jgcnt counts matrix elements
      jgcnt=0
c dk=1 p=1
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
c begin end transformation of dk=1 (the same for all p-case)
      do 45 ih=1,kmj
        do 45 m=1,kml
        mz=kz+m
cvp  matrixelemente nach neuer tabelle
cvp  matrix elements according to new table
          iem=(m-1)*kmj*irec+(ih-1)*irec
cvp
      do 131 l=1,nspp1
          tim=
     &     emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
     &             +nqrfal(l)*2-1 ) * cbfeld(l)
     &    +emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
     &             +nqrfal(l)*2   ) * exfeld(l)
c       if (mz.eq.505.and.nz+ih.eq.139) then
c       write(6,*) ' inopen=,nqr=,nql=,irec= ',inopen,nqr,nql,irec
c          write(6,*) ' mi=,mj=,m=,ih=,tim= ',mz,nz+ih,m,ih,tim
c          write(6,*) ' r=,qr=,ql= ',nrfal(l),nqrfal(l),nqlfal(l)
c          write(6,*) ' cb=,ex= ',cbfeld(l),exfeld(l)
c          write(6,*) ' iem= ',iem
c          write(6,*) ' emat1=,emat2= '
c    &    ,emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
c    &             +nqrfal(l)*2-1 )
c    &    ,emat(iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
c    &             +nqrfal(l)*2   )
c          write(6,*) ' arg-emat1=,arg-emat2= '
c    &    ,     iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
c    &             +nqrfal(l)*2-1
c    &    ,     iem+(nrfal(l)-1)*nqr*nql*2+(nqlfal(l)-1)*nqr*2
c    &             +nqrfal(l)*2
c       endif
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and
cvp  reset jgcnt to zero
          if (dabs(tim).ge.1.0e-7) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=mz
            mjfeld(jgcnt)=iaddr(l)+ih
          endif
  131 continue
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
ct    if (jgcnt.ge.ncmax) then
ct       do 4012 ik=1,jgcnt
ct         jg=jg+1
ct         q(jg)=qfeld(ik)
ct         mi(jg)=mifeld(ik)
ct         mj(jg)=mjfeld(ik)
ct         if (jg.ge.iwod) then
ct           jg=0
ct           write(nhuk) q,mi,mj
ct         endif
c4012    continue
ct       jgcnt=0
ct    endif
cvp
ct 45 continue
cvp
ct    do 145 ih=1,kmj
ct      do 145 m=1,kml
ct      mz=kz+m
cvp  matrixelemente nach neuer tabelle
          iem=(m-1)*kmj*irec2+(ih-1)*irec2 + iema
cvp
      do 132 l=nspp1+1,nsppa
          tim=
     &     emat(iem
     &             +nqrfal(l) ) * cbfeld(l)
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and
cvp  reset jgcnt to zero
          if (dabs(tim).ge.hstor) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=mz
            mjfeld(jgcnt)=iaddr(l)+ih
          endif
  132 continue
c schreiben der matrixelemente und der indices auf file nhuk
c write the matrix elements and the indices on file nhuk
ct    if (jgcnt.ge.ncmax) then
ct       do 4112 ik=1,jgcnt
ct         jg=jg+1
ct         q(jg)=qfeld(ik)
ct         mi(jg)=mifeld(ik)
ct         mj(jg)=mjfeld(ik)
ct         if (jg.ge.iwod) then
ct           jg=0
ct           write(nhuk) q,mi,mj
ct         endif
c4112    continue
ct       jgcnt=0
ct    endif
cvp
ct145 continue
cvp
cvp endtransformation for dk=1
ct    do 245 ih=1,kmj
ct      do 245 m=1,kml
ct      mz=kz+m
cvp  matrixelemente nach neuer tabelle
          iem=(m-1)*kmj*irec3+(ih-1)*irec3 + iemb
cvp
      do 43 l=nsppa+1,nsppb
          hp3(l)=
     &     emat(iem + (nrcas(l)-1)*nqr3 + nqrfal(l) ) * sumint(l)
  43  continue
cvp
      do 345 iis=1,inopen-2
c lesen der (negativen) austauschintegrale
c read the (negative) exchange integrals
      if (ioint.eq.1) then
      do 66 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(jto1+iis+1+(l-nsppa-1)*(inopen-1))
   66 continue
      endif
      if (ioint.eq.0) then
      do 67 l=nsppa+1,nsppb
        cbfeld(l)=-twoe(nwwmo(l,iis))
   67 continue
      endif
cvp   summation ueber die gemeinsamen offenen schalen
cvp   summation over common open shells
      do 68 l=nsppa+1,nsppb
          hp3(l)= hp3(l)
     &    +emat(iem + iis*2*nqr3
     &              + (nrcas(l)-1)*nqr3 + nqrfal(l) ) * cbfeld(l)
   68 continue
  345 continue
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
      do 69 l=nsppa+1,nsppb
          if (dabs(hp3(l)).ge.hstor) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=hp3(l)
            mifeld(jgcnt)=mz
            mjfeld(jgcnt)=iaddr(l)+ih
          endif
   69 continue
c schreiben der matrixelemente und der indices auf file nhuk
      if (jgcnt.ge.ncmax) then
         do 4212 ik=1,jgcnt
csut       jg=jg+1
           icount = icount+1
           rham(icount)=(qfeld(ik))
           imm(icount)=mifeld(ik)
           jmm(icount)=mjfeld(ik)
           if (icount.ge.ndims) then
             icount=0
             irech=irech+1
             write(nf88) imm,jmm,rham
             do iii=1,ndims
               rham(iii) = 0.0d0
               imm(iii)  = 0
               jmm(iii)  = 0
             end do
           endif
 4212    continue
         jgcnt=0
      endif
cvp
ct245 continue
   45 continue
c schreiben der matrixelemente und der indices auf file nhuk
      do 4013 ik=1,jgcnt
        icount=icount+1
        rham(icount)=(qfeld(ik))
        imm(icount)=mifeld(ik)
        jmm(icount)=mjfeld(ik)
        if (icount.ge.ndims) then
          icount=0
          irech=irech+1
          write(nf88) imm,jmm,rham
          do iii=1,ndims
            rham(iii) = 0.0d0
            imm(iii)  = 0
            jmm(iii)  = 0
         end do
        endif
 4013 continue
cvp hochsetzen vom jto1: wird nur for ioint=1 benoetigt
cvp increase jto1: needed only for ioint=1
      jto1=jto1+(nsppb-nsppa)*(inopen-1)
cvp
      return
      end
cvp
cvp berechnung der matrixelemente for dk=2
      subroutine hmdk2(twoe,
     +                 vecf1,vecf2,mdi,
     +                 emat,nedim,iottr,iotnew,iotm,
     +                 maindf,maxr,jconb,nomax2,
     =                 i,nc,mc,kml,kmj,
     +                 inopen,jz
     &                ,kz,mtype,jto1,ioint,iges,kend)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      dimension twoe(*)
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
cvp
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
cvp emat : darstellungsmatrizen
      real*8 emat
      dimension emat(nedim)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
c common for multi-root version
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common/rwork2/ ndim,nanz
c --- issk : superkategorien auf disk
      integer issk
      real*8   hstor
      common /csemi/ hstor,issk
cdebug
cvp
c schleife ueber alle konfigurationen
cdebug
*     write(6,*)'hmdk2:'
*     write(6,*)'emat',emat(1)
*     write(6,*)'iwod=',iwod
*     write(6,*)'i=',i
*     write(6,*)'kz=',kz
*     write(6,*)'jz=',jz
*     write(6,*)'inopen=',inopen
*     write(6,*)' nc=',nc
*     write(6,*)' mc=',mc
*     write(6,*)' kml=',kml
*     write(6,*)' kmj=',kmj
*     write(6,*)'jto1=',jto1
*     write(6,*)'ioint',ioint
*     write(6,*) (vecf1(ii),ii=1,6)
cdebug
cvp
cdebug
cvp  nqr = # qr-faelle
cvp  nqr = # qr-cases
      nqr=ibinom(inopen,4)
      nqr2=2*nqr
cvp  irec = # eintraege von emat for ein saf-paar
cvp  irec = # entries of emat for one saf-couple
      irec=6*nqr
c iwwg zaehlt alle ww konfigurationen for dk=2
c  iwwz=zaehler aller von null verschiedenen iww's
c iwwg counts all interacting configurations for dk=2
c  iwwz=count all non-zero iww's
      iwwg=0
      iwwz=0
cvp zeiten for konfigurationsvergleich gesamt
*     extimp=0.0d0
*     vutimp=0.0d0
*     cputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
*     vutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp lesen der integrale for dk=2, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=2 noch nicht fertig
c1012 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
cvp calculation of the number nstrip of iterations over the
cvp  configurations of smaller sk in case their number exceeds ncmax
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 12 k=kanf,kend
      do 12 k=1,nc
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp split the loop over all mc configurations in multiple loops
cvp  with maximum length ncmax
cvp mspiel = anzahl zu testender konfigurationen
cvp mspiel = number of configurations to be tested.
      mspiel=ncmax
      do 3012 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
cvp jnum1 = number of first configuration -1 to be tested
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
cvp read the integrals
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c iww zaehlt alle wechselwirkenden konfigurationen for dk=2
c iww counts all interacting configurations for dk=2
c     iww=0
c erzeugung des label-feldes for alle mc konfigurationen
c  und lesen der integrale
c generation of the label field for all mc configurations
c  and reading the integrals
c mit check2 werden alle wechselwirkenden konf. der sk(i-2) bestimmt
c sowie die labels
c with check2 all interacting configurations of the sk(i-2) are
c determined and the labels
*     call clockv(vu10,cpu10)
*     call time(mill10)
csut  call checkf(ioint,i,i-2,k,nspiel,mspiel,jnum1
      call check2(iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,i-2,k,nspiel,mspiel,jnum1)
cdebu write(6,*) nspiel,mspiel,jnum1
cvp
c
      do 3013 l=1,nspiel
 3013 iaddr(l)=jz+(ispiel(l)-1)*kmj
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
      do 3014 l=1,nspiel
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
 3014 continue
      jto1=jto1+2*nspiel
      endif
      if (ioint.eq.0) then
      do 3015 l=1,nspiel
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
 3015 continue
      endif
cvp
c     write(6,*)' iww=',nspiel
c     if (iww.ne.0) iwwz=iwwz+1
c jgcnt zaehlt die matrixelemente
      jgcnt=0
      do 24 ih=1,kmj
        do 23 m=1,kml
cvp    ileft  = linker index von h(ij), also i
        ileft  = kz+m
cvp  matrixelemente nach neuer tabelle
          iem=(m-1)*kmj*irec+(ih-1)*irec
cvp
      do 13 l=1,nspiel
          tim=
     &     emat(iem+(nrfal(l)-1)*nqr2+2*nqrfal(l)-1)*cbfeld(l)
     &    +emat(iem+(nrfal(l)-1)*nqr2+2*nqrfal(l)  )*exfeld(l)
            qfeld(l)=tim
  13  continue
c
c schleife ueber die wurzeln
c loop over roots
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
cvp gather-schritt
        do 401 l=1,nspiel
            vgath1(l)=vecf1(ipoint+iaddr(l)+ih)
            vgath2(l)=vecf2(ipoint+iaddr(l)+ih)
cdebug      write(6,*) 'vgath2(l)=',vgath2(l)
cdebug      write(6,*) 'vgath1(l)=',vgath1(l),
cdebu.      'ih,iaddr',ih,iaddr(l),nspiel
  401   continue
c   ende gather
         do 4012 l=1,nspiel
            qq=qfeld(l)
c           jjh=mjfeld(l)
c           vecf1(ileft)=vecf1(ileft)+vecf2(jjh)*qq
            vecf1(ipoint+ileft)=vecf1(ipoint+ileft)
     &                         +vgath2(l)*qq
 4012    continue
cvocl loop,novrec
         do 4013 l=1,nspiel
            qq=qfeld(l)
c           jjh=mjfeld(l)
            vgath1(l)=vgath1(l)+vecf2(ipoint+ileft)*qq
 4013    continue
cvp
cvp scatter-schritt
        do 411 l=1,nspiel
            vecf1(ipoint+iaddr(l)+ih)=vgath1(l)
cdebug
*           r4out1 = vecf1(ipoint+iaddr(l)+ih)
*           write(6,*) 'l,vecf1'
*           write(6,*) l,r4out1
cedbug
  411   continue
c ende schleife uber die wurzeln
 5100 continue
  23  continue
cvp
  24  continue
cvp
 3012 continue
cvp
  12  kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=2
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1012
cio   endif
      return
      end
cvp
cvp berechnung der matrixelemente for dk=2
cvp calculate matrixelements for dk=2
cvp
      subroutine hmdk2h(twoe,emat,nedim,iottr,iotnew,iotm,
     +                  maindf,maxr,jconb,nomax2,
     +                  i,nc,mc,kml,kmj,
     +                  inopen,jz,
     &                  kz,jto1,ioint,iges,kend)
cvp
cvp
cvp berechnung der matrixelemente for dk=2
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      real*8 twoe
      dimension twoe(*)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
cvp
*     real*8 f,ae
*     dimension f(ksafm*ksafm),ae(7829)
cvp emat : darstellungsmatrizen
cvp
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
c schleife ueber alle konfigurationen
cvp
cvp  nqr = # qr-faelle
      nqr=ibinom(inopen,4)
      nqr2=2*nqr
cvp  irec = # eintraege von emat for ein saf-paar
      irec=6*nqr
c iwwg zaehlt alle ww konfigurationen for dk=2
c  iwwz=zaehler aller von null verschiedenen iww's
      iwwg=0
      iwwz=0
cvp
cvp lesen der integrale for dk=2, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=2 noch nicht fertig
c1012 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
cvp calculation of the number nstrip of iterations over the
cvp  configurations of smaller sk in case their number exceeds ncmax
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 12 k=kanf,kend
      do 12 k=1,nc
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp split the loop over all mc configurations in multiple loops
cvp  with maximum length ncmax
cvp mspiel = anzahl zu testender konfigurationen
cvp mspiel = number of configurations to be tested.
      mspiel=ncmax
      do 3012 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
cvp jnum1 = number of first configuration -1 to be tested
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
cvp read the integrals
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c iww zaehlt alle wechselwirkenden konfigurationen for dk=2
c iww counts all interacting configurations for dk=2
c     iww=0
c erzeugung des label-feldes for alle mc konfigurationen
c  und lesen der integrale
c generation of the label field for all mc configurations
c  and reading the integrals
c mit check2 werden alle wechselwirkenden konf. der sk(i-2) bestimmt
c sowie die labels
c with check2 all interacting configurations of the sk(i-2) are
c determined and the labels
c     call checkf(ioint,i,i-2,k,nspiel,mspiel,jnum1
c    &  ,vutww,cpuww,vutla,cpula,vutin,cpuin)
      call check2(iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,i-2,k,nspiel,mspiel,jnum1)
cvp
c
      do 3013 l=1,nspiel
 3013 iaddr(l)=jz+(ispiel(l)-1)*kmj
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
      do 3014 l=1,nspiel
        cbfeld(l)= twoe(jto1+2*l-1)
        exfeld(l)=-twoe(jto1+2*l  )
 3014 continue
      jto1=jto1+2*nspiel
      endif
      if (ioint.eq.0) then
      do 3015 l=1,nspiel
        cbfeld(l)= twoe(intcb(l))
        exfeld(l)=-twoe(intex(l))
 3015 continue
      endif
cvp
c     write(6,*)' iww=',nspiel
c     if (iww.ne.0) iwwz=iwwz+1
c jgcnt zaehlt die matrixelemente
      jgcnt=0
      do 23 ih=1,kmj
        do 23 m=1,kml
cvp    ileft  = linker index von h(ij), also i
        ileft  = kz+m
cvp  matrixelemente nach neuer tabelle
          iem=(m-1)*kmj*irec+(ih-1)*irec
cvp
      do 13 l=1,nspiel
          tim=
     &     emat(iem+(nrfal(l)-1)*nqr2+2*nqrfal(l)-1)*cbfeld(l)
     &    +emat(iem+(nrfal(l)-1)*nqr2+2*nqrfal(l)  )*exfeld(l)
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cvp jgcnt can be kml*kmj times as big as the number of configurations
cvp  ---> at jgcnt=ncmax matrix elements write the lot out and
cvp  reset jgcnt to zero
          if (dabs(tim).ge.hstor) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
cx          mifeld(jgcnt)=kz+m
            mifeld(jgcnt)=ileft
            mjfeld(jgcnt)=iaddr(l)+ih
          endif
  13  continue
c schreiben der matrixelemente und der indices auf file nhuk
      if (jgcnt.ge.ncmax) then
         do 4012 ik=1,jgcnt
           icount=icount+1
           rham(icount)=(qfeld(ik))
           imm(icount)=mifeld(ik)
           jmm(icount)=mjfeld(ik)
           if (icount.ge.ndims) then
             icount=0
             irech=irech+1
             write(nf88)imm,jmm,rham
             do iii=1,ndims
                 imm(iii)=0
                 jmm(iii)=0
                 rham(iii)=0.0d0
             end do
           endif
 4012    continue
         jgcnt=0
      endif
cvp
  23  continue
c schreiben der matrixelemente und der indices auf file nhuk
      do 4013 ik=1,jgcnt
         icount=icount+1
         rham(icount)=(qfeld(ik))
         imm(icount)=mifeld(ik)
         jmm(icount)=mjfeld(ik)
         if (icount.ge.ndims) then
           icount=0
           irech=irech+1
           write(nf88)imm,jmm,rham
           do iii=1,ndims
              imm(iii)=0
              jmm(iii)=0
              rham(iii)=0.0d0
           end do
         endif
 4013 continue
cvp
 3012 continue
cvp
  12  kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=2
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1012
cio   endif
cvp
      return
      end
c
c#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
c
cvp berechnung des p=5-falles nach buenker
cvp calculation of the p=5-case following buenker
      subroutine hmp5(pey,acoul,aexc,core,ndeke,
     +                hp5,mdi,iottr,iotm,iwr,ofirst,mdif)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
      logical ofirst
      real*8 pey,core,acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
cvp diagonalelemente der hamiltonmatrix (bisher nur diagonale)
cvp diagonal elements of the hamilton matrix (sofar only diagonal)
c
      real*8 hp5
      integer mdi
      dimension hp5(mdi)
c
      integer iottr, iotm
      dimension iottr(iotm)
c
csut  ndet sichern
      integer idsut,idrem1,ndet1
      common /rhus/ idrem1,idsut(iswhm),ndet1(iswhm)
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
cvp
      real*8 af, ae
      common/rtable/ af(8000)
      dimension ae(7829)
      integer nsac,ndet,iaz,iez,jan,jbn,idra,idrc,jdrc 
      integer j9
      dimension nsac(iswhm),idrc(iswhm),
     4  idra(iswhm),jdrc(iswhm),jan(5),jbn(5),
     4  ndet(iswhm),iaz(iswhm),iez(iswhm),
     6  j9(49)
      equivalence (af(50),ae(1))
      equivalence (af(1),j9(1)),(j9(1),ndet(1)),(j9(6),
     1nsac(1)),(j9(11),iaz(1)),(j9(16),iez(1)),(j9(21),idra(1)),
     2(j9(26),idrc(1)),(j9(31),jdrc(1)),(j9(36),jan(1)),(j9(43),
     3jbn(1))
cvp
      common /junk/ f(ksafm*ksafm)
      common/linkmr/ic,i2,i3
      common/blksi3/nsz
      common/bufd/gout(510),nword
      dimension lout(510)
      equivalence (lout(1),gout(1))
      integer ot
      common /junk3/ ot(motmax)
      integer out, ohog
      dimension ohog(48),out(10962)
      equivalence (ot(13),ohog(1)),(ot(13),out(1))
cvp
      integer isc
      dimension isc(iswhm)
c
cvp  berechnung der dimension der matrix bis zur sk i-1
cvp  calculate the dimension of the matrix upto sk i-1
      jsum=0
      do 901 i=1,iswhm
      isc(i)=jsum
  901 jsum=jsum+nconf(i)*nsac(i)
csut
      do 902 i=1,iswhm
         ndet1(i) = (ndet(i))
  902 continue
cvp
cvp schleife ueber die superkategorien
cvp loop over superkategories
      call setsto(motmax,0,ot)
      do 1000 i=1,iswhm
cvp
      if (nconf(i).eq.0) goto 1000
cvp
      inytl=nytl(i)
      inopen=nod(i)
      ndt=ndet(i)
      kml=nsac(i)
      nps=nplu(i)
      nms=i-1
      niw=nps*nms
cvp
      ia=iaz(i)
      iz=isc(i)
      ie=iez(i)
      js=idra(i)
      ks=idrc(i)
      ix=nad1
      kss=0
102   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 102
      endif
*****
      call upackx(ot)
      kz=iz
cvp
      jstar=niot(i)-nconf(i)
      do 103 k=1,nconf(i)
cvp inytl = # besetzte mo's
      do 104 l=1,inytl
         kstar=jstar+nconf(i)*l
         itest(l)=iottr(kstar+k)
 104  continue
      do 105 l=1,kml*kml
105   f(l)=0.0d0
      care=core
cvp inopen = # einfach besetzte mo's
cvp inopen = # singly occupied mo's
      do 107 l=1,inopen
         mal=itest(l)
         ntar=nirred(mal)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
107      care=care+pey(mal)
cvp nd = # doppelt besetzte mo's
cvp nd = # doubly occupied mo's
cvp n1 = nr+1
      do 109 l=inopen+1,inytl
         mal=itest(l)
         ntar=nirred(mal)
         nal=ideks(mal+1)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
         sm=pey(mal)
109      care=care+sm+sm+acoul(nal)
      do 111 l=2,inopen
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 111 m=1,it
           nal=mal+itest(m)
111        care=care+acoul(nal)
      sm=0.0d0
cvp n2 = nr + 2
      do 113 l=inopen+2,inytl
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 113 m=inopen+1,it
           nal=mal+itest(m)
           tim=acoul(nal)
113        sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
      do 115 l=1,inopen
         mal=itest(l)
         do 115 m=inopen+1,inytl
           nal=itest(m)
           nal=ideks(max(nal,mal))+min(nal,mal)
           sm=acoul(nal)
115        care=care+sm+sm-aexc(nal)
      kk=ic
      kn=i2
      ml=i3
c write dk=0 p=5
      jx=-kml
      do 118 l=1,kml
        jx=jx+1
        lx=jx
        kx=ohog(l)+ia
        tim=care
cvp nps = # alpha-spins der offenen schalen
cvp nps = # alpha-spins of the open shells
        do 123 ih=2,nps
          mx=i2+nps*(l-1)+ih
          it=ih-1
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          do 123 in=1,it
            mz=i2+nps*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
123         tim=tim-aexc(nal)
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
cvp nms = #  beta-spins of the open shells
        do 124 ih=2,nms
          mx=i3+nms*(l-1)+ih
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          it=ih-1
          do 124 in=1,it
            mz=i3+nms*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
124         tim=tim-aexc(nal)
        do 121 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
121       f(lx)=f(lx)+tim*ae(kx)
cvp niw = nps * nms
        do 120 m=1,niw
          kk=kk+1
          kx=out(kk)+ia
c write dk=0 p=5 ausserdiagonal
c write dk=0 p=5 off-diagonal or without diagonal?
          kk=kk+1
          mal=out(kk)
          kk=kk+1
          nal=out(kk)
c         write(6,*)'kk,nal,mal = ', kk,nal,mal
          mal=itest(mal)
c         write(6,*)'itest(mal) = ', mal
c         write(6,*)'itest(nal) = ', itest(nal)
          nal=ideks(mal)+itest(nal)
          tim=-aexc(nal)
          lx=jx
          do 120 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
120         f(lx)=f(lx)+tim*ae(kx)
118   continue
      do 125 l=1,kml
        mz=kz+l
          m=l
          tim=0.0d0
          do 126 ih =1,kml
            ly=(m-1)*kml+ih
            mx=ie+(l-1)*kml+ih
126         tim=tim+f(ly)*ae(mx)
          if (mz.le.mdi) hp5(mz)=tim
125   continue
      if (mz.gt.mdi) then
       write(iwr,30) mz, mdi
30     format(1x,'insufficient memory for hamiltonian matrix'/
     +        1x,'required (mz)   = ', i12/
     +        1x,'available (mdi) = ', i12/)
         call caserr(
     + 'insufficient memory in hmp5 - increase mdi')
      endif
103   kz=mz
cvp
 1000 continue
c     print memory required
      if (ofirst) then
       mdif = mz
      endif
       write(iwr,40) mz, mdi
40     format(1x,'memory for hamiltonian matrix'/
     +        2x,'required (mz)   = ', i12/
     +        2x,'available (mdi) = ', i12)
cvp
      return
      end
      subroutine hqrii1(n,m,a,e,v,mx)
      implicit real*8 (a-h,o-z)
c
      dimension a(*), e(mx), v(mx,mx)
*************************************************************
*
* hqrii is a diagonalisation routine, written by yoshitaka beppu of
*       nagoya university, japan.
*       for details see 'computers & chemistry' vol.6 1982. page 000.
*
* on input    a       = matrix to be diagonalized
*             n       = size of matrix to be diagonalized.
*             m       = number of eigenvectors needed.
*             e       = array of size at least n
*             v       = array of size at least nmax*m
*
* on output   e       = eigenvalues
*             v       = eigenvectors in array of size nmax*m
*
* origin : mopac 5.0 / ibm 3090
************************************************************************

      common /junk2/ w(5,1000)

      if(n.eq.1 .and. m.eq.1) then
            e(1)=a(1)
            v(1,1)=1.d0
      return
      endif
*
* eps3 and eps are machine-precision dependent
*
      eps3=1.d-16
      zero=0.d0
      ll=(n*(n+1))/2+1
      eps=1.d-10
      iord=-1
      nm1=n-1
      if(n.eq.2) goto 80
      nm2=n-2
c     householder transformation
      do 70 k=1,nm2
         kp1=k+1
         w(2,k)=a((k*(k+1))/2)
         sum=0.0d0
         do 10 j=kp1,n
            w(2,j)=a((j*(j-1))/2+k)
   10    sum=w(2,j)**2+sum
         s=sign(dsqrt(sum),w(2,kp1))
         w(1,k)=-s
         w(2,kp1)=w(2,kp1)+s
         a(k+(kp1*(kp1-1))/2)=w(2,kp1)
         h=w(2,kp1)*s
         if(abs(h).lt.1.d-35) goto 70
c#      if(h.eq.0.d0) goto 70
         summ=0.d0
         do 50 i=kp1,n
            sum=0.0d0
            do 20 j=kp1,i
   20       sum=sum+a(j+(i*(i-1))/2)*w(2,j)
            if(i.ge.n) goto 40
            ip1=i+1
            do 30 j=ip1,n
   30       sum=sum+a(i+(j*(j-1))/2)*w(2,j)
   40       w(1,i)=sum/h
   50    summ=w(1,i)*w(2,i)+summ
         u=summ*0.5d0/h
         do 60 j=kp1,n
            w(1,j)=w(2,j)*u-w(1,j)
            do 60 i=kp1,j
   60    a(i+(j*(j-1))/2)=w(1,i)*w(2,j)+w(1,j)*w(2,i)+a(i+(j*(j-1))/2)
   70 a((k*(k+1))/2)=h
   80 w(2,nm1)=a((nm1*(nm1+1))/2)
      w(2,n)=a((n*(n+1))/2)
      w(1,nm1)=a(nm1+(n*(n-1))/2)
      w(1,n)=0.d0
      gersch=abs(w(2,1))+abs(w(1,1))
      do 90 i=1,nm1
   90 gersch=max(abs(w(2,i+1))+abs(w(1,i))+abs(w(1,i+1)),gersch)
      del=eps*gersch
      do 100 i=1,n
         w(3,i)=w(1,i)
         e(i)=w(2,i)
  100 v(i,m)=e(i)
      if(del.eq.zero)  goto  210
c     qr-method with origin shift
      k=n
  110 l=k
  120 if(abs(w(3,l-1)).lt.del) goto 130
      l=l-1
      if(l.gt.1)  goto 120
  130 if(l.eq.k)  goto 160
      ww=(e(k-1)+e(k))*0.5d0
      r=e(k)-ww
      z=sign(dsqrt(w(3,k-1)**2+r*r),r)+ww
      ee=e(l)-z
      e(l)=ee
      ff=w(3,l)
      r=dsqrt(ee*ee+ff*ff)
      j=l
      goto 150
  140 r=dsqrt(e(j)**2+w(3,j)**2)
      w(3,j-1)=s*r
      ee=e(j)*c
      ff=w(3,j)*c
  150 r=r+1.d-15
      c=e(j)/r
      s=w(3,j)/r
      ww=e(j+1)-z
      e(j)=(ff*c+ww*s)*s+ee+z
      e(j+1)=c*ww-s*ff
      j=j+1
      if(j.lt.k) goto 140
      w(3,k-1)=e(k)*s
      e(k)=e(k)*c+z
      goto 110
  160 k=k-1
      if(k.gt.1) goto 110
*    *    *    *    *    *    *    *    *    *    *    *    *
*
*   at this point the array 'e' contains the un-ordered eigenvalues
*
*    *    *    *    *    *    *    *    *    *    *    *    *
c     straight selection sort of eigenvalues
      sorter=1.d0
      if(iord.lt.0) sorter=-1.d0
      j=n
  170 l=1
      ii=1
      ll=1
      do 190 i=2,j
         if((e(i)-e(l))*sorter .gt. 0.d0) goto 180
         l=i
         goto 190
  180    ii=i
         ll=l
  190 continue
      if(ii.eq.ll) goto 200
      ww=e(ll)
      e(ll)=e(ii)
      e(ii)=ww
  200 j=ii-1
      if(j.ge.2) goto 170
  210 if(m.eq.0) return
***************
*  ordering of eigenvalues complete.
***************
c      inverse-iteration for eigenvectors
      fn=dble(n)
      eps1=1.d-5
      seps=dsqrt(eps)
      eps2=0.05d0
      rn=0.d0
      ra=eps*0.6180339887485d0
c    0.618... is the fibonacci number (-1+dsqrt(5))/2.
      ig=1
      do 430 i=1,m
         im1=i-1
         do 220 j=1,n
            w(3,j)=0.d0
            w(4,j)=w(1,j)
            w(5,j)=v(j,m)-e(i)
            rn=rn+ra
            if(rn.ge.eps) rn=rn-eps
  220    v(j,i)=rn
         do 250 j=1,nm1
            if(abs(w(5,j)).ge.abs(w(1,j))) goto 230
            w(2,j)=-w(5,j)/w(1,j)
            w(5,j)=w(1,j)
            t=w(5,j+1)
            w(5,j+1)=w(4,j)
            w(4,j)=t
            w(3,j)=w(4,j+1)
            if(w(3,j).eq.zero) w(3,j)=del
            w(4,j+1)=0.d0
            goto 240
  230       if(w(5,j).eq.zero) w(5,j)=del
            w(2,j)=-w(1,j)/w(5,j)
  240       w(4,j+1)=w(3,j)*w(2,j)+w(4,j+1)
  250    w(5,j+1)=w(4,j)*w(2,j)+w(5,j+1)
         if(abs(w(5,n)) .lt. eps3) w(5,n)=del
         do 310 itere=1,5
            if(itere.eq.1) goto 270
            do 260 j=1,nm1
               if(w(3,j).eq.zero) goto 260
               t=v(j,i)
               v(j,i)=v(j+1,i)
               v(j+1,i)=t
  260       v(j+1,i)=v(j,i)*w(2,j)+v(j+1,i)
  270       v(n,i)=v(n,i)/w(5,n)
            v(nm1,i)=(v(nm1,i)-v(n,i)*w(4,nm1))/w(5,nm1)
            vn=max(abs(v(n,i)),abs(v(nm1,i)),1.d-20)
            if(n.eq.2) goto 290
            k=nm2
  280       v(k,i)=(v(k,i)-v(k+1,i)*w(4,k)-v(k+2,i)*w(3,k))/w(5,k)
            vn=max(abs(v(k,i)),vn,1.d-20)
            k=k-1
            if(k.ge.1) goto 280
  290       s=eps1/vn
            do 300 j=1,n
  300       v(j,i)=v(j,i)*s
            if(itere.gt.1 .and. vn.gt.1) goto 320
  310    continue
c     transformation of eigenvectors
  320    if(n.eq.2) goto 360
         do 350 j=1,nm2
            k=n-j-1
            if(a((k*(k+1))/2).eq.zero) goto 350
            kp1=k+1
            sum=0.0d0
            do 330 kk=kp1,n
  330       sum=sum+a(k+(kk*(kk-1))/2)*v(kk,i)
            s=-sum/a((k*(k+1))/2)
            do 340 kk=kp1,n
  340       v(kk,i)=a(k+(kk*(kk-1))/2)*s+v(kk,i)
  350    continue
  360    do 370 j=ig,i
            if(abs(e(j)-e(i)) .lt. eps2) goto 380
  370    continue
         j=i
  380    ig=j
         if(ig .eq. i) goto 410
c     re-orthogonalisation
         do 400 k=ig,im1
            sum=0.0d0
            do 390 j=1,n
  390       sum=v(j,k)*v(j,i)+sum
            s=-sum
            do 400 j=1,n
  400    v(j,i)=v(j,k)*s+v(j,i)
c     normalisation
  410    sum=1.d-24
         do 420 j=1,n
  420    sum=sum+v(j,i)**2
         sinv=1.d0/dsqrt(sum)
         do 430 j=1,n
  430 v(j,i)=v(j,i)*sinv
      return
      end
cvp
cvp lesen der integrale zu iadl von ft31 und schreiben auf mtype
cvp read integrals into iadl from ft31 and write to mtype
      subroutine isort(twoe,nteint,ston,nidmax,
     +                 nston,nrec31,lend,iend)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
        integer nteint, nidmax
        real*8 twoe
        dimension twoe(nteint)
        real*8 ston
        dimension ston(nidmax)
cvp
cvp iadl: zwischenspeicherung der integraladressen
      real*8 sump3
      integer iadl,irec,iadp3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
cvp
cvp lanf+1 = erster eintrag von iadl, der auf twoe gechrieben wird
cvp nrec31 = anzahl der records von zweielektronenintegralen auf ft31
cvp irec = nummer des records der adresse iadl auf ft31
cvp iadl enthaelt jetzt die position der integraladresse auf irec
cvp
cvp lesen des files ft31 und abspeichern der integrale auf twoe
cvp irec = nummer des 'records' der adresse iadl auf sump3
cvp iadl enthaelt jetzt die position der integraladresse auf irec
cvp  irec31 = anzahl der ft31-records, die auf sump3 passen
cvp  icp31 = anzahl der fuellvorgaenge von sump3
cvp  icmpr = anzahl der integrale, die auf sump3 zwischengespeichert
cvp          werden
      irec31=icmax/nidmax
      icp31=(nrec31-1)/irec31+1
      icmpr=irec31*nidmax
      do 100 l=1,iend
        irec(l)=(iadl(l)-1)/icmpr+1
        iadl(l)=iadl(l)-(irec(l)-1)*icmpr
  100 continue
cvp
cvp
cvp lesen des files ft31 und abspeichern der integrale auf twoe
      call rewftn(nston)
      read(nston)
      read(nston)
      do 200 j=1,icp31-1
         lcnt=0
         do 210 i=1,irec31
            read(nston) ston
            do 230 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  230       continue
  210    continue
         do 220 l=1,iend
            if (irec(l).eq.j) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  220    continue
  200 continue
cvp  lesen der restlichen records von ft31
cvp   jrec31 = anzahl der restlichen records von ft31
         jrec31=nrec31-irec31*(icp31-1)
         lcnt=0
         do 250 i=1,jrec31
            read(nston) ston
            do 260 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  260       continue
  250    continue
         do 270 l=1,iend
            if (irec(l).eq.icp31) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  270    continue
cvp
cvp ausnullen des neu beschriebenen teils des integralfeldes
      do 190 j=1,icmax
  190 sump3(j)=0.0d0
cvp
c ft31 gelesen
csut  write(6,*) ' ft31 read '
c
      return
      end
cvp
cvp lesen der integrale zu iadl von ft31 und schreiben auf mtype
      subroutine isorta(twoe,nteint,ston,nidmax,
     +           nston,nrec31,indub,inopen,lend,iend
     +           ,ncntp3)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
        integer nteint, nidmax
        real*8 twoe
        dimension twoe(nteint)
        real*8 ston
        dimension ston(nidmax)
cvp
cvp iadl: zwischenspeicherung der integraladressen
      real*8 sump3
      integer iadl,irec,iadp3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
cvp
cvp lanf+1 = erster eintrag von iadl, der auf twoe gechrieben wird
cvp nrec31 = anzahl der records von zweielektronenintegralen auf ft31
cvp irec = nummer des 'records' der adresse iadl auf sump3
cvp iadl enthaelt jetzt die position der integraladresse auf irec
cvp  irec31 = anzahl der ft31-records, die auf sump3 passen
cvp  icp31 = anzahl der fuellvorgaenge von sump3
cvp  icmpr = anzahl der integrale, die auf sump3 zwischengespeichert
cvp          werden
      irec31=icmax/nidmax
      icp31=(nrec31-1)/irec31+1
      icmpr=irec31*nidmax
      do 100 l=1,iend
        irec(l)=(iadl(l)-1)/icmpr+1
        iadl(l)=iadl(l)-(irec(l)-1)*icmpr
  100 continue
cvp
cvp
cvp lesen des files ft31 und abspeichern der integrale auf twoe
      call rewftn(nston)
      read(nston)
      read(nston)
      do 200 j=1,icp31-1
         lcnt=0
         do 210 i=1,irec31
            read(nston) ston
            do 230 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  230       continue
  210    continue
         do 220 l=1,iend
            if (irec(l).eq.j) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  220    continue
  200 continue
cvp  lesen der restlichen records von ft31
cvp   jrec31 = anzahl der restlichen records von ft31
         jrec31=nrec31-irec31*(icp31-1)
         lcnt=0
         do 250 i=1,jrec31
            read(nston) ston
            do 260 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  260       continue
  250    continue
         do 270 l=1,iend
            if (irec(l).eq.icp31) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  270    continue
cvp
cvp ausnullen des neu beschriebenen teils des integralfeldes
      do 190 j=1,icmax
  190 sump3(j)=0.0d0
cvp
         icnt=0
         loop1=1
         loop2=iadp3(1)-1
         if (ncntp3.eq.0) loop2=iend
         do 310 l=loop1,loop2
              sump3(l)=twoe(lend+l)
  310    continue
         icnt=loop2
         do 300 j=1,ncntp3-1
           loop1=iadp3(j)
           loop2=iadp3(j)+2*indub+inopen-1
           icnt=icnt+1
           do 320 l=loop1,loop2
                sump3(icnt)=sump3(icnt)+twoe(lend+l)
  320      continue
           loop1=iadp3(j)+2*indub+inopen
           loop2=iadp3(j)+3*indub+inopen
           do 330 l=loop1,loop2
                sump3(icnt)=sump3(icnt)-twoe(lend+l)
  330      continue
           loop1=iadp3(j)+3*indub+inopen+1
           loop2=iadp3(j+1)-1
           do 340 l=loop1,loop2
              icnt=icnt+1
                sump3(icnt)=twoe(lend+l)
  340      continue
  300    continue
         if (ncntp3.gt.0) then
           loop1=iadp3(ncntp3)
           loop2=iadp3(ncntp3)+2*indub+inopen-1
           icnt=icnt+1
           do 350 l=loop1,loop2
                sump3(icnt)=sump3(icnt)+twoe(lend+l)
  350      continue
           loop1=iadp3(ncntp3)+2*indub+inopen
           loop2=iadp3(ncntp3)+3*indub+inopen
           do 360 l=loop1,loop2
                sump3(icnt)=sump3(icnt)-twoe(lend+l)
  360      continue
           loop1=iadp3(ncntp3)+3*indub+inopen+1
           loop2=iend
           do 370 l=loop1,loop2
              icnt=icnt+1
                sump3(icnt)=twoe(lend+l)
  370      continue
         endif
cvp
cvp setzen des integralfeldes
      do 400 j=1,icnt
  400 twoe(lend+j)=sump3(j)
cvp hochsetzen von lend
      lend=lend+icnt
c
c ft31 gelesen
csut  write(6,*) ' ft31 read '
c
      return
      end
cvp
cvp lesen der integrale zu iadl von ft31 und schreiben auf mtype
cvocl total,scalar
      subroutine isortb(twoe,nteint,ston,nidmax,
     +           nston,nrec31,indub,inopen,lend,iend
     &           ,ncntp3)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
        integer nteint, nidmax
        real*8 twoe
        dimension twoe(nteint)
        real*8 ston
        dimension ston(nidmax)
cvp
cvp iadl: zwischenspeicherung der integraladressen
      real*8 sump3
      integer iadl,irec,iadp3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
cvp
cvp lanf+1 = erster eintrag von iadl, der auf twoe gechrieben wird
cvp nrec31 = anzahl der records von zweielektronenintegralen auf ft31
cvp irec = nummer des 'records' der adresse iadl auf sump3
cvp iadl enthaelt jetzt die position der integraladresse auf irec
cvp  irec31 = anzahl der ft31-records, die auf sump3 passen
cvp  icp31 = anzahl der fuellvorgaenge von sump3
cvp  icmpr = anzahl der integrale, die auf sump3 zwischengespeichert
cvp          werden
      irec31=icmax/nidmax
      icp31=(nrec31-1)/irec31+1
      icmpr=irec31*nidmax
      do 100 l=1,iend
        irec(l)=(iadl(l)-1)/icmpr+1
        iadl(l)=iadl(l)-(irec(l)-1)*icmpr
  100 continue
cvp
cvp
cvp lesen des files ft31 und abspeichern der integrale auf twoe
      call rewftn(nston)
      read(nston)
      read(nston)
      do 200 j=1,icp31-1
         lcnt=0
         do 210 i=1,irec31
            read(nston) ston
            do 230 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  230       continue
  210    continue
         do 220 l=1,iend
            if (irec(l).eq.j) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  220    continue
  200 continue
cvp  lesen der restlichen records von ft31
cvp   jrec31 = anzahl der restlichen records von ft31
         jrec31=nrec31-irec31*(icp31-1)
         lcnt=0
         do 250 i=1,jrec31
            read(nston) ston
            do 260 l=1,nidmax
               lcnt=lcnt+1
               sump3(lcnt)=ston(l)
  260       continue
  250    continue
         do 270 l=1,iend
            if (irec(l).eq.icp31) then
              twoe(lend+l)=sump3(iadl(l))
            endif
  270    continue
cvp
cvp ausnullen des neu beschriebenen teils des integralfeldes
      do 190 j=1,icmax
  190 sump3(j)=0.0d0
cvp
         icnt=0
         loop1=1
         loop2=iabs(iadp3(1))-1
         if (ncntp3.eq.0) loop2=iend
         do 310 l=loop1,loop2
            icnt=icnt+1
              sump3(icnt)=twoe(lend+l)
  310    continue
         do 300 j=1,ncntp3-1
           loop1=iabs(iadp3(j))
           loop2=iabs(iadp3(j))+indub-1
           icnt=icnt+1
           do 320 l=loop1,loop2
                sump3(icnt)=sump3(icnt)-twoe(lend+l)
  320      continue
           loop1=iabs(iadp3(j))+indub
           loop2=iabs(iadp3(j))+3*indub+inopen-2
           if (iadp3(j).lt.0) loop2=loop2+1
           do 330 l=loop1,loop2
                sump3(icnt)=sump3(icnt)+twoe(lend+l)
  330      continue
           loop1=iabs(iadp3(j))+3*indub+inopen-1
           loop2=iabs(iadp3(j+1))-1
           if (iadp3(j).lt.0) loop1=loop1+1
           do 340 l=loop1,loop2
              icnt=icnt+1
                sump3(icnt)=twoe(lend+l)
  340      continue
  300    continue
         if (ncntp3.gt.0) then
           loop1=iabs(iadp3(ncntp3))
           loop2=iabs(iadp3(j))+indub-1
           icnt=icnt+1
           do 350 l=loop1,loop2
                sump3(icnt)=sump3(icnt)-twoe(lend+l)
  350      continue
           loop1=iabs(iadp3(j))+indub
           loop2=iabs(iadp3(j))+3*indub+inopen-2
           if (iadp3(ncntp3).lt.0) loop2=loop2+1
           do 360 l=loop1,loop2
                sump3(icnt)=sump3(icnt)+twoe(lend+l)
  360      continue
           loop1=iabs(iadp3(ncntp3))+3*indub+inopen-1
           loop2=iend
           if (iadp3(ncntp3).lt.0) loop1=loop1+1
           do 370 l=loop1,loop2
              icnt=icnt+1
                sump3(icnt)=twoe(lend+l)
  370      continue
         endif
cvp
cvp setzen des integralfeldes
      do 400 j=1,icnt
  400 twoe(lend+j)=sump3(j)
cvp hochsetzen von lend
      lend=lend+icnt
c
c ft31 gelesen
csut  write(6,*) ' ft31 read '
c
      return
      end
      subroutine mwrite(n,iwo,v)
      implicit real*8 (a-h,o-z)
      dimension v(n)
      write(iwo) v
*      write(6,*)'call mwrite',iwo,n
*      write(1,2)n
*2     format(3x,'the file should contain ',i7,'  entries'/)
*      do iraus = 1,n
*        write(1,1) iraus,v(iraus)
*      enddo
*1     format(3x,i7,3x,f15.10)
*      write(6,*) 'iraus after do-loop',iraus,' n=',n
*      write(6,*) 'write sequentially on 99'
*      write(99) v
*      write(6,*) 'after writing on 99'
*      write(6,*) 'before writing on ',iwo
*      write(6,*)'after writing on ',iwo
*      write(6,*)'stop in mwrite'
*      stop
      return
      end
c
c#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
c
      subroutine openfl
c
c  diese subroutineoeffnet die files
c
c   (1993) h.u.s.
c
cvp
      return
      end
c@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
      subroutine outp(oextrap,odebug)
c
c  ausdruck der extrapolation
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
      real*8 edavit,cdavit,extrapit,eigvalr,cradd, 
     +     weighb, rootdel, ethreshit
      integer mbuenk,nbuenk,mxroots
      logical ifbuen, ordel, odave, oweight
      logical odavit
      common /comrjb2/edavit(mxroot),cdavit(mxroot),
     +                extrapit(mxroot),ethreshit(mxroot),
     +                odavit(mxroot),
     +                eigvalr(maxref),cradd(3),weighb,
     +                rootdel,ifbuen,mbuenk,mxroots,
     +                ordel,odave,nbuenk,oweight
c
c
      logical oextrap, odebug
      parameter (mx=256)
      integer idim,nrotr
      common /c/ idim,nrotr
      integer nroot,ifirst,ilast
      integer ntch, kprin,ndeci, icode,konf ,keps,iggey
      integer istart,nrootz
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                ndeci, icode,konf ,keps,iggey,
     .                istart,nrootz
c --- scf information
      real*8 escf,vnuc
      common /kerne/ escf,vnuc
c
      integer iselec,izus,iselect,iselecz
      logical gganal
      common /cselec/ iselec(mxroot),izus(mxroot),
     +                iselect(mxroot),nrootci,
     +                iselecz(mxroot),gganal
c
      real*8 esav0,esav1,csum,de1,de0,ueb0,ueb1,egeys
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .               de0(mxroot),de1(mxroot),
     .               ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .               egeys(mxroot)
cbe ueb0(jlauf,ilauf)   jlauf = startvektor, ilauf = wurzel
cbe ueb1 natuerlich ebenso
cbe ialex0 enthaelt die ueberlapps for das grosse problem
cbe ialex1 for das kleine
cbe analoges gilt for ueber0 ueber1
      real*8 ueber0,ueber1
      dimension ueber0(mxroot),ueber1(mxroot),
     *ialex0(mxroot),ialex1(mxroot)
cbe de0w, de1w enthaelt die weighted energieerniedrigungen
      dimension de0w(mxroot),de1w(mxroot)
c
      real*8 deci,rlanda,eextp ,edav,cweight
      real*8 eextpo,edavo,ecio,csumo
      real*8 cwo
c --- old common from conmdi (partially)
      real*8 esav36
      integer ixt,i0
      common /per/ ixt,i0,esav36(4,mxroot)
      real*8 egey1,trash,tdel
      common /parkin/egey1,trash,tdel
c 
c+++++ GG extrapolation analysis 
c     arrays for finding pairs of roots with largest overlap
      real*8 totalS(mxroot,mxroot), GGlargest
      real*8 lVec(mxroot), sVec(mxroot)
      integer GGilargest(mxroot), isimilar(mxroot) 
      character*10 charwall
c
      itape = 36
      call vclr(de0w,1,mxroot)
      call vclr(de1w,1,mxroot)
c
      tlarge = trash
      tsmall = trash + tdel
c
      if (odebug) then
       write(iwr,*) 'outp: ***     debug prints      ***'
       write(iwr,*) 'outp: nroot=',nroot
       write(iwr,*) 'outp: nrotr=',nrotr
       write(iwr,*) 'outp: ilast=',ilast
       do ilauf = 1,ilast
        write(iwr,*) 
     +  'energy lowering with respect to start vector ',ilauf
        write(iwr,*) de0(ilauf),de1(ilauf)
       enddo
       write(iwr,*) 'outp: debug print ueb0,ueb1'
c
       do i=1,nroot
         iwurz=iselec(i)
         write(iwr,*)
         write(iwr,*) 'outp: before root ',iwurz
         do ilauf = 1,nrootz
           write(iwr,*) ilauf,iwurz,ueb0(ilauf,i),ueb1(ilauf,i)
         enddo
       enddo
      endif

      write(iwr,1106) cpulft(1) ,charwall()
1106  format(/1x,
     + '**** commence Buenker and Peyerimhoff extrapolation at',
     + f10.2,' seconds',a10,' wall'/)
      if (odebug) then
       write(iwr,*) 'outp: esav0,esav1,csum,de0,de1,egeys'
       do i=1,nroot
         write(iwr,*) esav0(i),esav1(i),csum(i),de0(i),de1(i)
         write(iwr,*) egeys(i),' outp: geyser-energy'
       end do
      endif
c
c ---- estimating the greatest overlaps with the start vectors
c
      do i=1,nroot
        if (odebug) then
         write(iwr,*)
         write(iwr,*) 'outp: calculating the overlap root :',i
         write(iwr,*) 'outp: =============================='
        endif
        ueber0(i)=0.0d0
        ueber1(i)=0.0d0
        ialex0(i)=i
        ialex1(i)=i
* ialex0    startvektor mit staerkstem ueberlap for large problem
* ueber0    ueberlap mit diesem startvektor for large problem
* ialex1    startvektor mit staerkstem ueberlap for small problem
* ueber1    ueberlap mit diesem startvektor for small problem
        do  ilauf = 1,nrootz
         if(dabs(ueb0(ilauf,i)).gt.dabs(ueber0(i))) then
           ueber0(i) = ueb0(ilauf,i)
           ialex0(i) = ilauf
         endif
         if (odebug) then
          write(iwr,*) 'outp: ueb0',ueb0(ilauf,i),'ueber0',ueber0(i)
         endif
         if(dabs(ueb1(ilauf,i)).gt.dabs(ueber1(i))) then
          ueber1(i) = ueb1(ilauf,i)
          ialex1(i) = ilauf
         endif
        enddo
      enddo
c
      if (gganal) then
c     GG extrapolation analysis
* --- decide which roots in the small and large problem correspond to each other

* --- iterate over large-problem vectors
      do i0=1,nroot
* ---   normalize large-problem vector
        vecNorm = 0.0d0
        do iz=1,nrootz
          lVec(iz) = ueb0(iz,i0)
                  vecNorm = vecNorm + lVec(iz)**2
        enddo
        vecNorm = dsqrt(vecNorm)
                do iz=1,nrootz
          lVec(iz) = lVec(iz) / vecNorm
        enddo
* ---   iterate over small-problem vectors
            do i1=1,nroot
* ---     normalize small-problem vector
          vecNorm = 0.0d0
          do iz=1,nrootz
            sVec(iz) = ueb1(iz,i1)
                    vecNorm = vecNorm + sVec(iz)**2
          enddo
          vecNorm = dsqrt(vecNorm)
                  do iz=1,nrootz
            sVec(iz) = sVec(iz) / vecNorm
          enddo
* ---     compute dot product of large and small vector
              totalS(i0,i1) = 0.0d0
                  do iz=1,nrootz
                    totalS(i0,i1)=totalS(i0,i1)+dabs(lVec(iz)*sVec(iz))
                  enddo
                enddo
          enddo

      write(iwr,*) 'GG: Overlap of roots in small and large problem:'
* --- iterate over all large problem roots
      do i0=1,nroot
            write(iwr,*)
            write(iwr,393) 'looking at large-problem root ', i0
             write(iwr,'(a25,1x,a15)') 'small problem root #', 'overlap'
                do i1=1,nroot
                  write(iwr,'(i25,1x,f15.5)') i1, totalS(i0,i1)
                enddo
* ----  find small problem root having largest overlap with large problem root
        GGlargest = 0.0d0
        do i1=1,nroot
                  if( totalS(i0,i1) > GGlargest ) then
                          GGlargest = totalS(i0,i1)
                          ismall = i1
                  endif
                enddo
                write(iwr,393) 'decided that large-problem root #', i0
        write(iwr,393) 'has maximum overlap with small-problem root #',
     +                  ismall
            isimilar(i0) = ismall
      enddo

      write(iwr,*)
          write(iwr,391) 'Large Problem','Small Problem','Large Problem'
      write(iwr,391) 'Root #','Root #','Zero Order CI Vec. #'
391   format( 2(a15,1x),1(a24,1x) )
      do i0 = 1, nroot
            write(iwr,392) i0, isimilar(i0), ialex0(i0)
          enddo
          write(iwr,*)
392   format( 2(i15,1x),1(i24,1x) )
393   format( a70, i3 )
394   format( a70 )

* --- print the corresponding vectors
      write(iwr,*) 'corresponding roots in Zero-Order CI basis:'
          do i0=1,nroot
        write(iwr,*)
                write(iwr,395) 'large problem root: ', i0,
     +  'small problem root: ', isimilar(i0)
        do iz = 1, nrootz
          write(iwr,396) ueb0(iz,i0), ueb1(iz,isimilar(i0))
        enddo
      enddo

395   format( a25, i5, a25, i5)
396   format( f30.3, f30.3 )
c
c     end of gg analysis section
c
      endif
cbe berechnung der weighted energieerniedrigungen
c
      do iroot=1,nroot
c        iwurz=iselec(iroot)
         x=0.0d0
         xz=0.0d0
         do istart=1,nrootz
            y = ueb0(istart,iroot)
            y2= y*y
            x = x + y2
c           de0w(iwurz) = de0w(iwurz) + de0(istart)*y2
            de0w(iroot) = de0w(iroot) + de0(istart)*y2
            yz = ueb1(istart,iroot)
            y2z= yz*yz
            xz = xz + y2z
c           de1w(iwurz) = de1w(iwurz) + de1(istart)*y2z
            de1w(iroot) = de1w(iroot) + de1(istart)*y2z
         enddo
c        de0w(iwurz) = de0w(iwurz)/x
c        de1w(iwurz) = de1w(iwurz)/xz
         de0w(iroot) = de0w(iroot)/x
         de1w(iroot) = de1w(iroot)/xz
      enddo
*
      if (odebug) then
* debug ausprint von de0w und de1w
       do iroot = 1,nroot
         iwurz = iselec(iroot)
         write(iwr,*) 'outp: de0w for root',iwurz,de0w(iroot)
         write(iwr,*) 'outp: de1w for root',iwurz,de1w(iroot)
       enddo
      endif

c
      iextra = 0
      do i=1,nroot
        iwurz=iselec(i)
        odavit(i) = .false.
        write(iwr,9820) iwurz
9820    format(5x,'***********'/
     +         5x,'root = ',i2,' *'/
     +         5x,'***********'//)
        write(iwr,9821)
9821    format(1x,'standard extrapolation'/
     +         1x,'======================'/)
        write(iwr,100) iwurz,tlarge,ialex0(i),ueber0(i),
     +                       tsmall,ialex1(i),ueber1(i)
100     format(1x,'root ',i3,
     +        ' has largest overlap with zero-order vector':/
     +        ' large problem (T=',f5.1,') :',i3,3x,f12.6/
     +        ' small problem (T=',f5.1,') :',i3,3x,f12.6/)
c
      if (gganal) then
c+++++ GG
* --- ASSIGN ismall - the small problem root # which correspond to the large problem root #
      ismall         = isimilar(i)
      ialex1(ismall) = ialex0(i)

          write(iwr,'(a70,i3)') 'large problem root has number : ', i
          write(iwr,'(a70,i3)') 'corresponding small problem root is: '
     +    ,ismall
          write(iwr,'(a70,f6.3)') 'overlap between these roots is: ',
     + totalS(i,ismall)
      write(iwr,'(a70,i3)') 'Zero Ord. CI Vector having maximum overlap
     +with large problem root is: ', ialex0(i)
      write(iwr,'(a70,i3)') 'Zero Ord. CI Vector having maximum overlap
     +with small problem root is: ', ialex1(ismall)
          write(iwr,*)
          iextra = 0
          if( totalS(i,ismall) < 0.5d0 ) then
            iextra = 1
            write(iwr,*) '*** WARNING ***'
            write(iwr,*) 'NO APPROPRIATE ROOT EXISTS IN SMALL PROBLEM'
            write(iwr,*) 'EXTRAPOLATION WILL FAIL'
            write(iwr,*)
            iextra = 1
          endif
c
      else
c
c       revert to original extrapolation analysis
c
        if(ialex0(i).eq.ialex1(i)) then
          ismall=i
        endif

        if(ialex0(i).ne.ialex1(i)) then
          odavit(i) = .true.
          write(iwr,102)
102       format(' ***** WARNING:'/
     +           ' different electronic structure ',
     +           'in both secular problems'/
     +           ' extrapolation can produce strange values')
          write(iwr,103) ialex0(i)
103       format(/' Extrapolation will be made using ',
     +    'for the smaller problem'/
     +    ' that root having largest overlap with zero order vector '
     +    ,i2)

* bestimmung der root, die im kleinen problem staerksten ueberlap
* mit dem startvektor hat, mit dem groesseres problem am
* staerksten ueberlappt (dies ist ialex0(i))

          ismall = 0
          do  ilauf = 1,ilast
           if (ialex1(ilauf).eq.ialex0(i)) then
            ismall=ilauf
           endif
          enddo
* falls diese nicht existiert wird keine extrapolation durchgefuehrt
          if (ismall.eq.0) then
           write(iwr,104)
104        format(' *** warning: no appropriate root found for'/
     +            ' smaller problem'/
     +            ' *** no extrapolation will be performed')
           iextra=1
          else
           write(iwr,105) ismall
105   format(' This is root ',i2/)
          endif

        endif
c
        delta = 0.0d0
        if(ismall.ne.0)
     +  delta = dabs(esav0(i) - esav1(ismall))
        if(delta.le.1.0d-6.and.ismall.ne.0) then
         write(iwr,1105) 
1105     format(/' **** identical energies at both thresholds ***'/
     +           ' **** no extrapolation will be performed    ***'/)
         iextra = 1
        endif
c
      endif
c
c --- lambda
c
* bestimmung von lambda nur falls iextra=0
* ansonsten existiert keine passende wurzel im kleinen problem
        if (iextra.eq.0) then
c
c       do the extrapolation

cbe an dieser stelle muss de1(ialex0(i)) genommen werden, da man natuerli
cbe die energieerniedrigung bzgl eines startvektors und nicht eines
cbe wurzel nehmen muss, beide jetzt betrachteten wurzeln, also i und isma
cbe ueberlappen am meisten mit dem startvektor der mit i am meisten ueber
cbe deswegen de1(ialex0(i))

         write(iwr,101) tlarge,ialex0(i),esav0(i),de0(ialex0(i)),
     +                  tsmall,ialex1(i),esav1(ismall),de1(ialex0(i))
101      format(1x,'used for calculation of lambda :'/
     +             1x,'(T=',f5.1,') E-CI(',i2,') ',2f15.6/
     +             1x,'(T=',f5.1,') E-CI(',i2,') ',2f15.6)
c
c   check that both upper and lower roots roots actually obtained in
c   this calculation 
c   
         jjj = ialex0(i)
         do j = 1, nrootz
	  if(iselecz(j) .eq. jjj) go to 1205
         enddo
         iextra = 1
1205     jjj = ialex1(i)
         do j = 1, nrootz
          if(iselecz(j) .eq. jjj) go to 1206
         enddo
         iextra = 1
1206     if (iextra.ne.0) then
          odavit(i) = .true.
          write(iwr,1109) ialex0(i), ialex1(i)
1109      format(/
     + ' **** either upper (',i2,') or lower (',i2,') root',
     + ' has not been calculated ***'/
     + ' **** no extrapolation will be performed             ***'/)
           go to 1000
         endif
         rlanda = (esav0(i)-esav1(ismall))
         rlanda =rlanda/(de0(ialex0(i))-de1(ialex0(i)))
         write(iwr,8822) rlanda
8822     format(/1x,'optimal lambda =',f13.8)
         write(iwr,*)
         eextp = esav0(i) - rlanda*de0(ialex0(i))
         ecio = esav0(i)
chs
        else if (iextra.eq.1)then
          eextp = esav0(i)
          ecio = esav0(i)
        endif
c       write(iwr,'(a,f15.8,a)') ' ci-energy            :',ecio,
c    .  ' hartree.'
        write(iwr,9835)  ecio, iwurz
 9835   format(1x,'ci energy           = ',f15.8, 
     +            ' hartree  (root ',i2,')')
        deci = eextp-egeys(ialex0(i))-escf
        csumo   = csum(i)
        cweight = 1.d0-csum(i)
        edav    = eextp + cweight*deci
        eextpo = eextp
c
        if (iextra.eq.0)then
c
c       data for Iterative MRDCI treatment
c
        cdavit(i) = csumo
        edavit(i) = edav
        extrapit(i) = eextp
        ethreshit(i) = ecio
        write(iwr,9836)  eextpo, iwurz, csumo
 9836   format(1x,'extrapolated energy = ',f15.8,3x,
     +  '[   root ',i2,'  c**2      = ',f10.4,' ]')
        else
         write(iwr,*)' **** no extrapolation'
         cdavit(i) = 0.0d0
        endif
c
        edavo = edav
        cwo   = cweight
        write(iwr,9837)  edavo, cwo
 9837   format(1x,'davidson            = ',f15.8,3x,
     +  '[   1-c**2             = ',f10.4,' ]')
c
* save data for final output to ft36 (itape)
c
c++++++ GG
        if (gganal) then
         write(iwr,*)
         write(iwr,*) 'SAVING STANDARD EXTRAPOLATION RESULTS TO FTN36'
         write(iwr,*)
c++++++ GG  need to comment in assigning to esav36 if we so not want to save standard extrap data
        endif
        esav36(1,i) = ecio
        esav36(2,i) = eextpo
        esav36(3,i) = edavo
        esav36(4,i) = csumo
c
        cweight = cweight/csum(i)
        edav    = eextp + cweight*deci
        edavo = edav
        cwo   = cweight
        write(iwr,9838)  edavo, cwo
 9838   format(1x,'weights davidson1   = ',f15.8,3x,
     +  '[  (1-c**2)/ c**2      = ',f10.4,' ]')
        cweight = 1.d0-csum(i)
        cweight = cweight/(csum(i)+csum(i)-1.0d0)
        edav    = eextp + cweight*deci
        edavo = edav
        cwo   = cweight
        write(iwr,9839)  edavo, cwo
 9839   format(1x,'weights davidson2   = ',f15.8,3x,
     +  '[  (1-c**2)/(2*c**2-1) = ',f10.4,' ]'/)
*       note the use of davidson rather than weighted
*       davidson date to ft36 (for compatability with old code)
*       write(itape) ecio,eextp,edav,csumo

        write(iwr,9822)
9822    format(1x,'weighted extrapolation'/
     +         1x,'======================'/)
cbe an der stelle nimmt man alles bzgl i!
c
        if(ismall.eq.0) then
         write(iwr,101) tlarge,ialex0(i),esav0(i),de0w(i)
        else
         write(iwr,101) tlarge,ialex0(i),esav0(i),de0w(i),
     +                  tsmall,ialex1(i),esav1(ismall),de1w(i)
        endif
        ecio = esav0(i)
        write(iwr,9841)  ecio
 9841   format(1x,'ci energy           = ',f15.8, ' hartree ')
        if (iextra.eq.1) go to 1000
c
        rlanda = (esav0(i)-esav1(ismall))
        rlanda =rlanda/(de0w(i)-de1w(i))
        write(iwr,8822) rlanda
        write(iwr,*)
        eextp = esav0(i) - rlanda*de0w(i)
        deci = eextp-egeys(ialex0(i))-escf
        csumo   = csum(i)
        cweight = 1.d0-csum(i)
        edav    = eextp + cweight*deci
        eextpo = eextp
        write(iwr,9840)  eextpo, csumo
 9840   format(1x,'extrapolated energy = ',f15.8,3x,
     +  '[    c**2              = ',f10.4,' ]')
        edavo = edav
        cwo   = cweight
        write(iwr,9837)  edavo, cwo
        cweight = cweight/csum(i)
        edav    = eextp + cweight*deci
        edavo = edav
        cwo   = cweight
        write(iwr,9838)  edavo, cwo
        cweight = 1.d0-csum(i)
        cweight = cweight/(csum(i)+csum(i)-1.0d0)
        edav    = eextp + cweight*deci
        edavo = edav
        cwo   = cweight
        write(iwr,9839)  edavo, cwo
        write(iwr,*)
c++++++ GG
        if (gganal) then
c+++GG        write(iwr,*)
c+++GG        write(iwr,*) 'WRITING WEIGHTED EXTRAPOLATION RESULTS TO FTN36'
c+++GG        write(iwr,*)
c+++GG	esav36(1,i) = ecio
c+++GG        esav36(2,i) = eextpo
c+++GG        esav36(3,i) = edavo
c+++GG        esav36(4,i) = csumo
        endif
c
*       write(itape) ecio,eextp,edav,csumo
* ansprungaddresse falls keine extrapolation
1000    continue
      iextra=0

      end do
c -----------------------------------
c     add final record to ft36
      write(itape) esav36
      oextrap = .true.
c
      write(iwr,1107) cpulft(1) ,charwall()
1107  format(/1x,
     + '**** Buenker and Peyerimhoff extrapolation complete at',
     + f10.2,' seconds',a10,' wall'/)
c
      return
      end
cvp
cvp for ioint=1 wird hier die vorsortierung der integrale, wie in
cvp  skiny benoetigt werden, durchgefuehrt. for die p=3-faelle
cvp  werden gleich die integralsummen gebildet
      subroutine preint(twoe,nteint,ston,nidmax,
     +              pey, acoul, aexc, ndeke, core,
     +              pey0, acoul0, aexc0, noeint,
     +              iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              nrec31,iisk,iesk)
csut
csut   for iisk, bis iesk
csut
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
        integer nteint, nidmax
        real*8 twoe
        dimension twoe(nteint)
        real*8 ston
        real*8 pey, core, pey0
        real*8 acoul,aexc,acoul0,aexc0
        dimension ston(nidmax)
        dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
        dimension pey0(noeint),acoul0(noeint),aexc0(noeint)
      integer nform
      common /cffor/ nform
c
      integer iottr,iotm,iotnew
      dimension iottr(iotm),iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
      common/rb/ ifrk,iswh
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
cvp
cvp felder for vorsortierung der integrale
cvp iadl: zwischenspeicherung der integraladressen
c
      real*8 sump3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
      integer iadl,irec,iadp3
      integer iisk,iesk
cvp
cvp iend : fuellstand von iadl
cvp lend : fuellstand von twoe
cvp jfull: nummer der konfiguration, bei der twoe noch nicht voll
cvp ncntp3 zaehlt die p=3-faelle auf iadl
        iend=0
        lend=0
        jfull=0
        ncntp3=0
cvp das hilfsfeld muss kleiner als das integralfeld sein
        icmaxm=min(icmax,nteint)
cvp iadp3 enthaelt adresse des ersten summanden auf iadl
        do 85 i=1,icmax
          iadl(i) = 0
          irec(i) = 0
          sump3(i) = 0.0d0
          iadp3(i)=0
  85    continue
cvp     integralbestimmung for die einzelnen dk-faelle
cold     do 11 i=1,iswh
         do 11 i=iisk,iesk
           if (nconf(i).eq.0) goto 11
           inopen=nod(i)
           indub=ndub(i)
cvp        dk = 2
            if (i.ge.3) then
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
              mc=nconf(i-2)
              nstrip=(mc-1)/ncmax+1
cvp
              if (mc.ne.0) then
csut
              do 110 j=1,nconf(i)
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
              mspiel=ncmax
              do 3012 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
                jnum1=(istrip-1)*ncmax
                if (istrip.eq.nstrip) then
                   mspiel=mc-(nstrip-1)*ncmax
                endif
cvp
                call dk2int(iottr,iotnew,iotm,
     +               maindf,maxr,jconb,nomax2,
     +               i,i-2,j,nspiel,mspiel,jnum1)
                if (lend+iend.gt.nteint) then
                   call writs1(mtype,twoe,lend,jfull)
                   lend=0
                endif
                if (iend+2*nspiel.gt.icmaxm) then
                   call isort(twoe,nteint,ston,nidmax,
     +             nston,nrec31,lend,iend)
                   jfull=j-1
                   lend=lend+iend
                   iend=0
                endif
                do 210 k=1,nspiel
                  iadl(iend+2*k-1)=intcb(k)
                  iadl(iend+2*k  )=intex(k)
  210           continue
                iend=iend+2*nspiel
 3012         continue
cvp
  110         continue
csut
csut          end if
cvp
cvp         abspeichern der restlichen integrale for dk=2
              if (lend+iend.gt.nteint) then
cvp             lend+iend > nteint kann nur aus der letzten
cvp              konfiguration folgen
                 call writs1(mtype,twoe,lend,nconf(i)-1)
                 call isort(twoe,nteint,ston,nidmax,
     +                      nston,nrec31,0,iend)
                 call writs1(mtype,twoe,iend,nconf(i))
               else
                 call isort(twoe,nteint,ston,nidmax,
     +                      nston,nrec31,lend,iend)
                 call writs1(mtype,twoe,lend+iend,nconf(i))
              endif
              iend=0
              lend=0
              jfull=0
            endif
csut
              end if
cvp        dk = 1
            if (i.ge.2) then
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
              mc=nconf(i-1)
              nstrip=(mc-1)/ncmax+1
cvp
              if (mc.ne.0) then
csut
              do 120 j=1,nconf(i)
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
              mspiel=ncmax
              do 3028 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
                jnum1=(istrip-1)*ncmax
                if (istrip.eq.nstrip) then
                  mspiel=mc-(nstrip-1)*ncmax
                endif
cvp
                if (lend+iend.gt.nteint) then
                   call writs1(mtype,twoe,lend,jfull)
                   lend=0
                endif
                call dk1int(twoe,nteint,ston,nidmax,
     +             iottr,iotnew,iotm,
     +             maindf,maxr,jconb,nomax2,
     +             i,i-1,j,nspiel,mspiel,jnum1
     &            ,icmaxm,nston,nrec31,lend,iend,jfull,ncntp3)
 3028         continue
cvp
  120         continue
cvp
csut          end if
csut
cvp         abspeichern der restlichen integrale for dk=1
              if (lend+iend.gt.nteint) then
cvp             lend+iend > nteint kann nur aus der letzten
cvp              konfiguration folgen
                 call writs1(mtype,twoe,lend,nconf(i)-1)
                 lend=0
                 m0 = 0
                 call isorta(twoe,nteint,ston,nidmax,
     +                nston,nrec31,indub,inopen,m0,iend
     &                ,ncntp3)
                 ncntp3=0
                 do 86 j=1,icmax
  86             iadp3(j)=0
                 call writs1(mtype,twoe,lend,nconf(i))
               else
                 call isorta(twoe,nteint,ston,nidmax,
     +                nston,nrec31,indub,inopen,lend,iend
     &                ,ncntp3)
                 ncntp3=0
                 do 87 j=1,icmax
  87             iadp3(j)=0
                 call writs1(mtype,twoe,lend,nconf(i))
              endif
              iend=0
              lend=0
              jfull=0
            endif
            end if
csut
cvp        dk = 0
              do 130 j=1,nconf(i)
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
              nstrip=(j-2)/ncmax+1
cvp
cvp
cvp zerlegung der schleife ueber alle j-1 konfigurationen
cvp  der unteren dreiecksmatrix in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
              mspiel=ncmax
              do 3129 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
                jnum1=(istrip-1)*ncmax
                if (istrip.eq.nstrip) then
                   mspiel=j-1-(nstrip-1)*ncmax
                endif
cvp
                if (lend+iend.gt.nteint) then
                   call writs1(mtype,twoe,lend,jfull)
                   lend=0
                endif
                call dk0int(twoe,nteint,ston,nidmax,
     +             iottr,iotnew,iotm,
     +             maindf,maxr,jconb,nomax2,
     +             i,i,j,nspiel,mspiel,jnum1
     &            ,icmaxm,nston,nrec31,lend,iend,jfull,ncntp3)
 3129         continue
cvp
  130         continue
cvp
cvp         abspeichern der restlichen integrale for dk=0
              if (lend+iend.gt.nteint) then
cvp             lend+iend > nteint kann nur aus der letzten
cvp              konfiguration folgen
                 call writs1(mtype,twoe,lend,nconf(i)-1)
                 lend=0
                 m0 = 0
                 call isortb(twoe,nteint,ston,nidmax,
     +                nston,nrec31,indub,inopen,m0,iend
     &                ,ncntp3)
                 ncntp3=0
                 do 96 j=1,icmax
  96             iadp3(j)=0
                 call writs1(mtype,twoe,lend,nconf(i))
               else
                 call isortb(twoe,nteint,ston,nidmax,
     +                nston,nrec31,indub,inopen,lend,iend
     &                ,ncntp3)
                 ncntp3=0
                 do 97 j=1,icmax
  97             iadp3(j)=0
                 call writs1(mtype,twoe,lend,nconf(i))
              endif
              iend=0
              lend=0
              jfull=0
   11    continue
cvp
cvp lesen der einelektronenintegrale sowie der core-energie
c letzter record von nston
         if (nform.gt.0) then
            read(nston) pey
            read(nston) acoul
            read(nston) aexc
            read(nston) core
corig       read(nston) pey,acoul,aexc,core
         else
            read(nston) pey0,acoul0,aexc0,core
            do 3333 ii=1,ndeke
               pey(ii)   = pey0(ii)
               acoul(ii) = acoul0(ii)
               aexc(ii)  = aexc0(ii)
3333        continue
         end if
*         write(iwr,*) ' last record read from ft31'
cvp
         call rewftn(mtype)
c
      return
      end
      subroutine prep0(f1,mdi,d,y,len2,len,odebug)
cbe iselec(ilast) in ilast veraendert
c
c    bereite den start-vektor vor
c
c    auf ft20 ::
c
      implicit none
c
      real*8 d, y
      integer len2, len
      dimension d(len2), y(len)
      real*8 f1
      integer mdi
      dimension f1(mdi)
      logical odebug
c
      integer jsum
      integer iswhm,luns
      integer i,ii,idim,nc
      integer mx
      integer mxroot
      parameter (iswhm=5)
      parameter (mx=256)
      parameter (mxroot=50  )
c   startvektor
      real*8 crest
c   parameter for ft77
      integer ifrk
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      integer nstone, mstvt
      common /ftap/ nstone(10), mstvt
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      integer ixx,iyy,iy,ib,ibg,ibb
      integer itr,iq,istart
      integer iqq,nroop
      integer izz
      integer iroot,ibut
      integer nv
      integer mipos
      integer ipos,iii
      real*8 bn, trsum, ew
      common /miscop/ trsum(10,mxroot),ew(mxroot),bn(500)
      real*8 c95
      dimension nv(mxroot)
c     integer nff
      integer id,nvaa,nvee
      real*8 yend
c   parameter for rumple-file
      integer iwod,imo,m
      integer ndt,kml
      integer nconf,kconf
      dimension nconf(iswhm),kconf(iswhm)
c
      integer iselec
      integer nroot, ifirst, ilast, ntch, kprin
      integer ndeci, icode,konf ,keps,iggey
c von b.e. eingefuegt in unkenntnis
      common /cselec/ iselec(mx)
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                ndeci, icode,konf ,keps,iggey,istart
      real*8 zero,vnuc
c --- alter common
      common /c/ idim
      real*8 an
      common /can/ an(5292)
      common /cmpos/ mipos(mx)
      common /kerne/ zero,vnuc
cbe      write(6,*) 'call prep0.f'
c
c --- nhuk : rumple-file
c
      ifrk  = 500
      call rewftn(nhuk)
c
c --- luns : startfile
c
      luns  = lun20
c
      jsum=0
c
c     read rumple tape and fudge up parameters
c
      read (nhuk) iwod,vnuc,zero,imo,m,nconf
c
cdebug      write(6,*) 'first record of rumple read!'
      do 9322 i=1,iswhm
       kconf(i)=nconf(i)
9322  continue
      do 3000 i=1,iswhm
        nc      = nconf(i)
        if(nc .eq. 0) go to 3000
c
c     copy from rumple file
c
        read (nhuk) ndt,kml,an
co      kmk(i) = kml
        jsum    = jsum + nc*kml
cor     read (nhuk) ihog,imap,en
        read (nhuk)
 3000 continue
      idim = jsum
      ibut = jsum
      read (nhuk) (bn(i),i=1,5),ibg
clear write (6,*) 'first read nhuk'
*     citr = bn(4)
*     cito = bn(5)
clear
      jsum = (jsum-1)/ifrk+1
      do i=1,jsum
        read(nhuk)  bn
        write(mousa) bn
      end do
clear      stop
cdebugwrite(6,*) 'last record of rumple read!'
c --- mache unit vector
*     istart = 0
      if (istart.eq.0) then
        write(iwr,771) 
771     format(/1x,'diagonal vector taken as CI start vector')
c --- lokalisiere mains
        call rewftn(mousa)
        iy=1
        ib=ifrk
        itr=((idim-1)/ifrk)+1
c
c   thresholds : 940.0d0 : main
c                999.0d0 : r-main
c               1000.0d0 : singleexcitation !
c
c
        do iqq=1,itr
        ibb = ib
c       allow for array bounds checking
        if(ib.gt.mdi) ibb = mdi
         read (mousa)  (f1(iq),iq=iy,ibb)
         iy=iy+ifrk
         ib=ib+ifrk
        end do
        ipos = 0
        do ixx=1,idim
          if (f1(ixx).gt.930.0d0) then
            if (f1(ixx).lt.998.0d0) then
               ipos = ipos+1
               mipos(ipos) = ixx
            end if
          end if
        end do
cdebug  write(6,*) 'mousa read !!!'
        c95  = 0.95d0
        crest = sqrt(1.0d0 - c95*c95) / dble(idim)
*        write(6,*) 'ilast ',ilast
        do ii=1,ilast
          do i=1,idim
            f1(i) = crest
          end do
          f1(ii) = c95
c --- schreibt startvektor
          write (luns) (f1(i),i=1,idim)
        end do
c --- benuetze geyservektor
      else if (istart.eq.1) then
        write(iwr,772) 
772   format(1x,'zero-order vectors used as start vectors')
c --- extract from unit=mstvt
        write(iwr,*) 'start vectors read'
c       open (unit=mstvt,file='gstarv.dat',status='old',
c    .      form='unformatted')
        call rewftn(mstvt)
        read(mstvt) nroop
        read(mstvt) d
        read(mstvt) y
c
        iii = 0
        do 9593 ii=1,nroop
          iii=iii+ii
9593    continue
        yend = y(iii)
c ---
        id=nroop*(nroop+1)/2
c       write(6,*) 'prediagonalisation out parkwa2 with'
        if (odebug) write(6,*) 'prep0: nroop=',nroop,' read.'
        nvee=nroop*nroop
        nvaa=nvee-nroop+1
c       nff=nv(nroot)
c       if(ndeci .ne. 0) nff=nroop
c
        call rewftn(mousa)
        iy=1
        ib=ifrk
        itr=((idim-1)/ifrk)+1
c
c   thresholds : 940.0d0 : main
c                999.0d0 : r-main
c               1000.0d0 : singleexcitation !
c
c
        do 84 iqq=1,itr
c        allow for array bounds checking
         ibb = ib
         if(ib.gt.mdi) ibb = mdi
         read (mousa)  (f1(iq),iq=iy,ibb)
         iy=iy+ifrk
         ib=ib+ifrk
   84   continue
cdebug        write (6,'(5(f14.8,1x))') (f1(ixx),ixx=1,1500)
c zuweisung : position des vektors
        ipos = 0
        do ixx=1,idim
          if (f1(ixx).gt.930.0d0) then
            if (f1(ixx).lt.998.0d0) then
               ipos = ipos+1
               mipos(ipos) = ixx
            end if
          end if
        end do
cdebug  write(6,*) 'mousa read !!!'
cdebug  write(6,*) 'd'
cdebug  write(6,'(5(1x,f12.8))') d
cdebug  write(6,*) 'ilast ',ilast
        do ii=1,ilast
          iroot = ii
          ipos = nroop*nroop - (iroot)*nroop
          do iyy=1,idim
             f1(iyy) = 0.0d0
          end do
          do izz=1,nroop
            ipos = ipos + 1
            f1(mipos(izz)) = d(ipos)
cdebug      write(iwr,*) 'mipos(izz),d(ipos),ipos'
cdebug      write(iwr,*) mipos(izz),d(ipos),ipos
          end do
c --- schreibt startvektor
          write(iwr,773) ii, idim
773       format(1x,'starting vector',i3,' has dimension ', i7)
          write (luns) (f1(i),i=1,idim)
        end do
      end if
c
      call rewftn(luns)
cbe
c     write(6,*) 'build mipos field'
c     write(6,*) (mipos(iraus),iraus=1,20)
c
c106  format(1x,i4,2x,f20.6,3x,4(d12.4,' (',i7,') '))
c107  format(//' *** startvektor aus parkwa2 ***'//,
c    *1x,'root         energy',9x,'eigenvector'//)
      return
      end
c
cvp lesen der information bzgl. erzeuger und vernichter
      subroutine rderz(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +                 mtape,iwr,ofirst,iotmf)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
      logical ofirst
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer iotnew, iotm
      dimension iotnew(iotm)
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
cvp
c    test, ob konfigurationen auf iotnew passen
c     write(6,*) 'enter rderz: iotmf =', iotm
      iotnst(1)=0
      do 10 i=2,iswhm
c  for refcon bzw. refcon1
c       iotnst(i)=iotnst(i-1)+9*nconf(i-1)
        iotnst(i)=iotnst(i-1)+5*nconf(i-1)
   10 continue
      if (ofirst) then
          iotmf = max(iotmf,iotnst(iswhm))
c         write(6,*) 'rderz: iotmf =', iotmf
      else
       if (iotnst(iswhm).gt.iotm) then
         write(iwr,30) iotnst(iswhm), iotm
30       format(1x,'insufficient memory for configurations'/
     +        1x,'required (iotst) = ', i12/
     +        1x,'available (iotm) = ', i12/)
         call caserr(
     + 'insufficient memory for configurations - increase iotm')
       endif
      endif
c
c     read(mtape) nko,maindf,iref
      read(mtape) nko,maindf,iref,jconb
      read(mtape) iotnew
*       write(iwr,*)'from iderz short before end'
*       write(iwr,1)(iotnew(iraus),iraus=1,100)
*      write(iwr,1)(iotnew(iraus),iraus=101,200)
*      write(iwr,1)(iotnew(iraus),iraus=201,300)
*      write(iwr,1)(iotnew(iraus),iraus=301,400)
*      write(iwr,1)(iotnew(iraus),iraus=401,500)
*1     format(100i2)
c
      return
      end
c
c#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
c
      subroutine readinad(nteint,mdi,iotm,nedim,
     +                    onteint,omdi,oiotm,onedim,oprint,odebug)
c
c  This subroutine reads the user-input for the conmdi program
c
c   - namelist-convention following the fortran90-standard
c
c   (1993) h.u.s.
c
cvp
cvp  input of ioint
cvp           ioint = 0: all integrals are kept in core.
cvp           ioint = 1: integrals get sorted
c
csut input of iggey   : number of geyser-vectors that are used in
csut                    the strtvector for the Davidson
csut
c
csut input of izus(j) =-1 : one root at a time
csut                        not yet implemented
csut
c
csut input of norhp  =1 : don't use root-homing
csut                      not yet implemented
csut
c
csut input of nform  =1 : new ft31 format (no efeld,cf)
csut
c
csut input of ndav  =0 : diagonalisation following the original - dvdson
cbe doesn't work at the moment
csut                =1 : diagonalisation using the multiroot-version
csut
c
csut cptol : tolerance for printing each configuration (0.002d0)
c
      implicit real*8 (a-h,o-z)
c
      logical onteint,omdi,oiotm,onedim
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
      real*8 hstor
      real*8 critan,critrn,criten,orthon
      integer nform
      integer ifirst,ilast,i,istart
      integer i2out
      integer nroot,ntch,iggey
      integer kprin,ndeci,icode,konf,keps
      integer nzero,one,ndav
c     integer norhp
      integer issk,iswhm
      logical odebug,oprint
      parameter (mxroot = 50)
      parameter (nzero = 0)
      parameter (one   = 1)
      parameter (iswhm = 5)
c
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
      common /thresh/ critan,critrn,criten,orthon
      common /cffor/  nform
      common /rwork2/ ndim,nanz,nglei
c
      integer iselec,izus,iselect,iselecz
      common /cselec/ iselec(mxroot),izus(mxroot),
     +                iselect(mxroot),nrootci,
     +                iselecz(mxroot)
c
c --- issk 1..issk : superkategorien auf disk
      common /csemi/  hstor,issk
      common /cdiag_adler/  ndav
      common /cinput1/ nroot, ifirst, ilast, ntch, kprin,
     .                ndeci, icode,konf ,keps,iggey,
     .                istart,nrootz
      real*8 hstori, rtoli, cptoli, cptolmi, cptolcci
      logical bypass, debugd, debugtab, debugci
      logical onteintin,omdiin,oiotmin,onedimin,oprintin
      common /adlin/ hstori, rtoli, cptoli, cptolmi, cptolcci,
     +               nrooti, ntchi, 
     +               kprini, iseleci,
     .               ndecii, icodei, konfi, kepsi,
     .               iointi, norhpi,
     .               izusi, isski, 
     .               ifirsti, ilasti, istarti, ndavi,
     .               iggeyi, nformi,
     +               bypass(24),debugd, debugtab, debugci,
     +               nteintin,mdiin,iotmin,nedimin,
     +               onteintin,omdiin,oiotmin,onedimin,
     +               oprintin,nglein
c --- set default-values
ctest      write(6,*) 'call subroutine readin'
c     issk  = nzero
c     nform = nzero
c     nrec31= nzero
c     istart= nzero
c     ioint = nzero
c     nroot = nzero
c     ntch  = one
c     iggey = nzero
c     kprin = nzero
c     ndeci = one
c     icode = nzero
c     ndav  = nzero
c     konf  = nzero
c     keps  = nzero
c     norhp = nzero
c     critan = 0.00001d0
c     critrn = 0.00001d0
chs   criten = 0.00001d0
c     criten = 0.000001d0
c     hstor  = 0.0000001d0
c     orthon = 0.0000001d0
c     rtol  = 0.005d0
c     cptol = 0.002d0
c
      nteint = nteintin
      mdi    = mdiin
      iotm   = iotmin
      nedim  = nedimin
      onteint = onteintin
      omdi = omdiin
      oiotm = oiotmin
      onedim = onedimin
      oprint = oprintin
      nglei = nglein
c
      do i=1,mxroot
c       iselec(i) = nzero
        izus(i)   = nzero
      end do
c
      nrootz = nrooti
      if (odebug)write(6,*) '****  from readinad ******'
      nroot= nrootci
      do i = 1,nroot
       iselec(i) = iselect(i)
      enddo
      if (odebug) then
       write(6,*) 'nroot, nrootz = ', nroot, nrootz
       write(6,*)' iselecz=',(iselecz(i),i=1,10)
       write(6,*)' iselec =',(iselec(i),i=1,10)
      endif
      ntch= ntchi
      kprin= kprini
      ndeci= ndecii
      icode= icodei
      konf= konfi
      keps= kepsi
      ioint= iointi
      rtol= rtoli
      cptol= cptoli
      cptolcc = cptolcci
      cptolm= cptolmi
c     norhp= norhpi
      issk= isski
      hstor= hstori
      ifirst= ifirsti
      ilast= ilasti
      istart= istarti
      ndav= ndavi
      iggey= iggeyi
      nform= nformi
c     if (ndav.ne.1) then
c       write(6,*) ' only multiroot possible for the moment '
c       write(6,*) ' ndav changed to one '
c       ndav = 1
c     endif
cbe falls keine eingabe for nroot, nroot setzen ueber iselec
      if (nroot.eq.nzero) then

         if (iselec(1).ne.nzero) then
           nroot=1
           do ibe = 2,mxroot
             if(iselec(ibe).ne.nzero) nroot=nroot+1
           enddo
         endif
      write(6,*) 'nroot from readin ',nroot

      endif
cbe falls eingabe for iselec, ueberschreiben von ifirst,ilast
      if (iselec(1).ne.nzero) then
         ifirst=iselec(1)
         ilast =iselec(nroot)
      endif
cbe falls keine eingabe for iselec
      if (iselec(1).eq.nzero) then
         if (ifirst.eq.nzero) then
cbe und keine eingabe for ifirst, die ersten nroot wurzeln
           do ibe = 1,nroot
             iselec(ibe) = ibe
           enddo
             ifirst = 1
             ilast  = nroot
         endif
cbe aber eine eingabe for ifirst, dann die ersten nroot wurzeln
cbe beginnend bei ifirst
         if (ifirst.ne.nzero) then
           do ibe = 1,nroot
             iselec(ibe)=ifirst+ibe-1
           enddo
         endif
      endif
cbe      if (nroot.ne.nzero) then
cbe        ifirst=1
cbe        ilast = iselec(nroot)
cbe      end if
      if (odebug) then
       write(6,*)'readinad: ilast=',ilast,'  ifirst=',ifirst
       write(6,*)'readinad: iselec=',(iselec(iraus),iraus=1,10)
       write(6,*)'readinad: izus',(izus(iraus),iraus=1,5)
      endif

c --- fehlerbehandlung
      if (iggey.eq.nzero) then
         iggey = nroot
      end if
      if (ntch.gt.one) then
         write(iwr,*) 'warning : ntch set to 1 !'
         ntch = one
      end if
      if (iggey.lt.one) then
         iggey = one
      end if
c --- ausdrucken der input parameter
      write(iwr,*)
      write(iwr,*)
     . ' ========================================================='
      write(iwr,*)
     . ' === input-parameters                                  ==='
      write(iwr,*)
     . ' ===---------------------------------------------------==='
      write(iwr,'(a,i7,a)')
     . '  === number of roots                           ',nroot,' ==='
      write(iwr,'(a,i7,a)')
     . '  === first  root                               ',ifirst,' ==='
      write(iwr,'(a,i7,a)')
     . '  === last root                                 ',ilast,' ==='
      write(iwr,'(a,i7,a)')
     . '  === number of secular equations               ',ntch,' ==='
      write(iwr,'(a,1pe9.2,a)')
     . '  === ci-tolerance                            ',rtol,' ==='
      write(iwr,'(a,1pe9.2,a)')
     . '  === printoff (ci coefficients)              ',cptol,' ==='
      write(iwr,'(a,1pe9.2,a)')
     . '  === printoff (c**2)                         ',cptolcc,' ==='
      write(iwr,'(a,i7,a)')
     . '  === number of starting vectors                ',iggey,' ==='
      write(iwr,'(a,i7,a)')
     . '  === print flag                                ',kprin,' ==='
      write(iwr,'(a,i7,a)')
     . '  === ndeci                                     ',ndeci,' ==='
      write(iwr,'(a,i7,a)')
     . '  === icode                                     ',icode,' ==='
      write(iwr,'(a,i7,a)')
     . '  === konf                                      ',konf,' ==='
      write(iwr,'(a,i7,a)')
     . '  === keps                                      ',keps,' ==='
      write(iwr,*)
     . ' ========================================================='
      if (nglei.gt.0) then
        write(iwr,'(a,i7,a)')
     +  '  === number of vectors on 1st pass            ',nglei,' ==='
      endif
      if (ioint.eq.nzero) then
         write(iwr,*)
     . ' === all integrals are kept in memory                  ==='
      else
         write(iwr,*)
     . ' === sorting of integrals                              ==='
      end if
      write(iwr,*)
     . ' ========================================================='
      if (nform.eq.nzero) then
         write(iwr,*)
     . ' === old file 31 format is assumed                     ==='
      else
         write(iwr,*)
     . ' === ft31 format of gamess-uk is assumed               ==='
      end if
      write(iwr,*)
     . ' ========================================================='
      if (istart.eq.nzero) then
         write(iwr,*)
     . ' === norm vector is used for start                     ==='
         write(iwr,*)
     . ' ========================================================='
      else if (istart.eq.1) then
         write(iwr,*)
     . ' === reference vectors are used for start              ==='
      write(iwr,*)
     . ' ========================================================='
      else
         istart = nzero
      end if
      if (issk.eq.nzero) then
        write(iwr,*)
     . ' === direct Table-CI                                   ==='
      else if (issk.ge.iswhm ) then
        issk = iswhm
        write(iwr,*)
     . ' === conventional mrdci : consider disc-space          ==='
        write(iwr,*)
     . ' === hamiltonmatrix with h>',hstor
      else
        i2out = issk
        write(iwr,*)
     . ' === semidirect mrdci issk =',i2out,' super categories onto ',
     . 'disc'
        write(iwr,*)
     . ' === hamiltonmatrix with h>',hstor
      end if
      write(iwr,*)
     . ' ========================================================='
      write(iwr,*)
      write(iwr,35)(iselec(iraus),iraus=1,nroot)
35    format(' ***** selected roots  ',40i3)
c
csut oeffne h-mat.file
      if (issk.gt.nzero) then
        open(unit=88,file='hamil.mat',form='unformatted',
     +   status='unknown')
      end if
c
      return
      end
      subroutine rumpxs(twoe,nteint,ston,nidmax,
     +   pey, acoul, aexc, ndeke, core,
     +   pey0, acoul0, aexc0, noeint,
     +   hp5, mdi,
     +   emat,nedim,vect,mxmr,
     +   iottr,iotnew,iot0,iotm,
     +   maindf,maxr,jconb,nomax2,odebug,
     +   ofirst,nteintf,iotmf,mdif)
c
c    - aus rumpled : vorbereitung der felder for hamiltonmatrixelemente
c
c      einbau in conmdi 1993
c
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$ rumpled $$$$$$$$$$$$$$$$$$$$$$$
c
c  die berechnung der matrixelemente erfolgt mit hilfe der
c   sga-darstellungsmatrizen vektorisiert; die reihenfolge der
c   abgespeicherten matrixelemente unterscheidet sich von der
c   nach der buenkerschen table erzeugten (vp 26.2.1993)
c
c  for p=1, p2 und p=3 werden die sga-darstellungsmatrizen mit
c   hilfe der table erzeugt, falls sie for den jeweiligen dk-fall
c   auf edim passen (vp 20.2.1993)
c
c  zerlegung der schleife ueber alle konfigurationen in mehrere
c   kleinere mit maximaler laenge ncmax
c   (vp 5.2.1993)
c
c  inputparameter ioint (anstatt ibunk) eingefuehrt
c     ioint = 0 : alle integrale werden im hauptspeicher gehalten
c                 und die integraladressen bei bedarf berechnet
c     ioint = 1 : integraladressen werden zu beginn bestimmt und
c                 die zugehoerigen integrale auf ft09 abgespeichert
c   (vp dezember 1992)
c
c  i/o der files 31, 33 und 34 von skiny in main gebracht
c   (vp 13.11.1992)
c
c  p-faelle sind for dk=1 und dk=2 sortiert, d.h. zuerst werden
c   alle matrixelemente zu p=1, dann zu p=2 usw. berechnet
c   (vp 26.10.1992)
c
c  dimension des iot-feldes auf 600000 heraufgesetzt (vp 19.10.1992)
c   was bedeutet der parameter 'ifdmax', der auch den wert 180000
c   hat und hier nicht veraendert wurde?
c
c (vp september 1992) for alle dk-faelle werden die wechselwirkenden
c  mo's ueber die q-faelle bestimmt
c
c (vp juli 1992) in checkc werden for dk=2 die integraladressen
c  berechnet; zur bestimmung der mo-nummern ist variante b1
c  realisiert
c
c (vp april 1992) labels zur berechnung der matrixelemente werden in
c  check0, check1 und check2 bestimmt
c
c bestimmung der wechselwirkenden konfigurationen for dk=2 in check2,
c for dk=1 in check1 und dk=0 in check0; zusaetzlich erfolgt die
c bestimmung des p-falles (vp 16.3.1992)
c
c in skiny und compr1 wird for dk=0,p=5 ft33 (aus parkeu) nicht mehr
c gelesen, die noetige information aus dem feld iot wird via
c common/vergl/ uebergeben (vp 9.3.1992)
c
c compr1 ist umgestricktes skiny:
c in compr1 werden for dk=1,p=3 und dk=0,p=3 summen ueber drei-
c zentrenintegrale sowie die austauschintegrale einzeln abgespeichert
c --> erzeugung eines komprimierten integral-files nhuh2=68
c in skiny wird ft68 gelesen (vp 5.3.1992)
c
c label-file wird in rumple und beary nur noch mit werten
c wechselwirkender konf.paare beschrieben -->
c in skiny wird nur noch ueber die wechselwirkenden konf.gelaufen,
c mit check1 werden jeweils die positionen der ww. konf. bestimmt
c (vp 28.2.1992)
c
c iot in rumple und beary auf integer gesetzt,
c ablegen auf common/vergl/ --> wird in check1 benoetigt
c niot heisst in rumple und beary icon
c (vp 25.2.1992)
c
c integrale zur sk(i-2) werden for alle 3 z-faelle vor berechnung
c der matrixelemente in einer schleife auf bobvec abgespeichert
c (vp 24.2.1992)
c
c in dk=2 schreiben der matrixelemente und indices hinter die
c schleife zur sk(i-2) gezogen
c
c in dk=2 lesen des label- und des integral-files vor die schleife
c zur sk(i-2) gezogen (vp 5.2.1992)
c
c nahezu alle equivalence-statements in skiny beseitigt (vp 31.1.1992)
c
c felder im stily auf  100,000 zurueckgesetzt, damit job in class=r
c laeuft (vp 30.1.1992)
c
c felder im stily auf  2000,000 gesetzt (be 04.03.91)
c
c bonn, 23.06.90  write statements zum checken von fi35
c nach aachen geschickter rumpled (gross)
c
c bernd engels fredericton, den 03.08.89
c
c write statements rauskommentiert
c normaler rumple mit write und comentarstatements
c die eingefuegten write statements sind je nach dk-faellen jeweils
c durch eine kommentarkarte gekenzeichnet, z.b.:
c
c       c write dk=2
c innerhalb der p unterfaelle wird auch nochmals getrennt unterschieden
c
c dieser rumple schreibt nicht mehr kc,kd auf fi39 heraus, sondern
c die gesamtadresse. die umrechnung geschieht in stily selbst.
c dadurch hat man den integraladressenfile um 1/3 gekuerzt
c mehrere test haben die richtigkeit der version bestaetigt
c
c die equivalence statements in main, beary und stily sind durch den
c common/a/ block ersetzt. das equivalence innerhalb des skiny ist
c erhalten. zwecks platzsparen hat der neue common/a/ block den gleichen
c namen wie der alte.
c
c iot, nymax wurde auf 180000 heraufgesetzt. ideks auf 4951 hoch
c
c saemtliche timer routinen sind jedoch leider herauskommentiert
c ebenso das erste write statement auf papier enthaelt nicht
c die versionsnummer.
c
c zusammen mit christel wurden die felder im stily hochgesetzt
c in dieser version stehen sie auf 300000 und diese version ist
c daher unter 5000 k altes thch nicht verwendbar
c
c f77 - sm77-2 - mva
c
c//delc&      deactivate timer routines
c>@v  versionsnummer
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
      logical odebug, ofirst
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      integer nteint, nidmax
      real*8 twoe
      dimension twoe(nteint)
      real*8 ston
      dimension ston(nidmax)
      real*8 pey, core, pey0
      real*8 acoul,aexc,acoul0,aexc0
      integer ndeke, noeint
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
      dimension pey0(noeint),acoul0(noeint),aexc0(noeint)
c
      real*8 vect
      integer mxmr
      dimension vect(mxmr)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
      integer nform
      common /cffor/ nform
c
      real*8 hp5
      integer mdi
      dimension hp5(mdi)
cvp
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iot0, iotm
      dimension iottr(iotm), iotnew(iotm), iot0(iotm)
c
      integer mx
      parameter (mx=256)
***** problem with extrapolation
      integer iselec,izus,iselect,iselecz,nrootci
      common /cselec/ iselec(mxroot),izus(mxroot),
     +                iselect(mxroot),nrootci,
     +                iselecz(mxroot)
      common/rb/ ifrk,iswh
c
      real*8 esav0,esav1,csum,de0,de1,ueb0,ueb1,egeys
      common /cpert/ esav0(mxroot),esav1(mxroot),csum(mxroot),
     .               de0(mxroot)  ,de1(mxroot),
     .               ueb0(mxroot,mxroot),ueb1(mxroot,mxroot),
     .               egeys(mxroot)
c 
      real*8 b,trsum,ew
      integer kj,mjx,lj,ij
      integer nj,ntil,nbal,isym,jsym,lsym,ncomp
      integer ibal,itil,mcomp
      integer istm,isp
      integer kmap,khog,jkan
      common /miscop/ trsum(10,mxroot),ew(mxroot),b(500),
     + lj(nirmax),ij(8),kj(nirmax),mjx(nirmax),
     + nj(nirmax),ntil(nirmax),nbal(9),
     + isym(nirmax),jsym(jabmax),lsym(800),ncomp(100),
     + ibal(nirmax),itil(nirmax),mcomp(100),istm(10),
     + jkan(1000),kmap(504),khog(48)
c
      integer icount,irech
      common /rhsc/ icount,irech
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
csut issk superkategorien auf disk
csut  anzahl sac's pro superkategorie
      integer ncons,idrem1
      common /rhus/ idrem1,ncons(iswhm)
clear common /rfzero/ fzero
cvp
cvp erster record von ft31
      real*8 vnuc,zero
      integer nsel,iorbs,knu
cvp zweiter record von ft31
      real*8 ab, eb
      common/junk/ ab(5292),eb(2304)
c
cvp
cvp weitere groessen, die von ft33 gelesen werden
      integer ndt,kml
cvp uebergabe von informationen von ft31 und ft33 an skiny
c
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu ,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
cvp felder for vorsortierung der integrale
cvp iadl: zwischenspeicherung der integraladressen
      real*8 sump3
      integer iadl,irec,iadp3
      common/rwork/ iadl(icmax),irec(icmax),iadp3(icmax)
     &       ,sump3(icmax)
cvp
csut timer
      real*8 cpu0,cpu1,dcpu
cvp for integralvorsortierung nach foxy
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
cvp
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
cvp vom programm benoetigte files
c
c  file unit nos. now held in ftap5 and initialised in adler
c
c
      call rewftn(mtype)
      call rewftn(nston)
      call rewftn(ltape)
      call rewftn(ideli)
      call rewftn(nhuk)
      call rewftn(mtapev)
c//skipc&1
c write allgemein
      iswh=iswhm
      ifrk=500
cvp
      ifr1=1-ifrk
csut  nad=4000
csut  nad1=-3999
cvp
cvp input von ioint
cvp           ioint = 0: alle integrale werden im core gehalten
cvp           ioint = 1: integrale werden vorsortiert
cio   read(ird,2) ioint
cio   2  format(14i5)
cold  ioint = 0
cio   write(iwr,*)
cio   write(iwr,*) ' ioint= ',ioint
cio   write(iwr,*)
cvp   if (ibunk.ne.0) go to 882
cvp  erzeugung von ideks(i) = i ueber 2
      ideks(1)=0
      do 1 i=1,ideksm-1
   1  ideks(i+1)=ideks(i)+i
      if (nform.gt.0) then
        read(nston) n,iwod,nid,ksum,imo,kj,mjx,lj,nj,nsel,ntil,nbal
     &   ,isym,jsym,iorbs,knu,lsym,ncomp,vnuc,zero
c?????&
csut    nbal(9)=ncomp(1)
c??????
         read(nston) nit,lg,ij,ibal,itil,mcomp
      else
csut --- altes format !
        read(nston) n,iwod,nid,ksum,imo,kj,mjx,lj,nj,nsel,ntil,nbal
     &   ,isym,jsym,iorbs,knu,lsym,ncomp,vnuc,zero
c****** nbal(9)=ncomp(1)
c??????
         read(nston) nit,lg,ij,ibal,itil,mcomp
      end if
c -------------------------------------------------------
        nbal(9)=ncomp(1)
cvp
cvp umsetzen von ioint, falls integrale nicht in den core passen
      if (ofirst) then
       nteintf = lg
      endif
      if (lg.le.nteint) then
        write(iwr,997) lg, nteint
997     format(1x,'memory for 2-e integrals'/
     +         2x,'required =            ',i8,' R*8 words'/
     +         2x,'available (nteint)  = ',i8,' R*8 words')
      else
         ioint=1
         write(iwr,*)
         write(iwr,*) ' integrals will not fit in core '
         write(iwr,999) lg, nteint
999      format(1x,'core required  = ',i9,' R*8 words'/
     +          1x,'core available = ',i9,' R*8 words'/)
         write(iwr,*) ' parameter ioint redefined   '
         write(iwr,*) ' ioint= ',ioint
         write(iwr,*)
      endif
cvp
cvp umspeichern der startadressen der einelektronenintegrale auf ijone
      do 231 i=1,nirmax
        ijone(i)=ij(i)
  231 continue
cvp nrec31 = anzahl der records von zweielektronenintegralen auf ft31
      nrec31=(lg-1)/nid+1
cvp abspeichern der zweieletronenintegrale auf twoe
cvp  for ioint = 0
      if (ioint.eq.0) then
         icoint=0
         do jdoo=1,nrec31
           read(nston) ston
           do kdoo=1,nid
            icoint=icoint+1
            twoe(icoint)=ston(kdoo)
           enddo
         enddo
cvp lesen der einelektronenintegrale sowie der core-energie
c letzter record von nston
         if (nform.gt.0) then
           read(nston) pey
           read(nston) acoul
           read(nston) aexc
           read(nston) core
         else
           read(nston) pey0,acoul0,aexc0,core
           do 3049 ii=1,noeint
             pey(ii)   = pey0(ii)
             acoul(ii) = acoul0(ii)
             aexc(ii)  = aexc0(ii)
3049       continue
         end if
*         write(iwr,*) ' last record of ft31 read '
      endif
cvp
      call rewftn(nston)
cvp # der ci-mo's pro irred. darst. wird auf ncimo abgelegt
      do 3050 i=1,nirmax
 3050 ncimo(i)=lj(i)
cvp
      ix=0
      do 5 i=1,n
      kap=lj(i)
      if (kap.eq.0) go to 5
      do 6 j=1,kap
      ix=ix+1
cvp
      nirred(ix)=i
      mopos(ix)=j
cvp
   6  continue
   5  continue
cvp
      read (ltape) jsec,nrootx,nytl,nplu,ndub,vect,ew,mconf
      do i=1,mxroot
         egeys(i) = ew(i)
      end do
      read (ideli) (b(i),i=1,5),isp,vect,ew
cvp
      ny=0
      do 7 i=1,iswh
      niot(i)=ny
      nod(i)=nytl(i)-ndub(i)
      if (mconf(i).eq.0) go to 7
      ig=1
      nx=nytl(i)
      jg=0
      read (ltape)
      read (ltape)
      read (ltape) jkan
   8  if (jkan(ig).eq.0) go to 9
      jg=jg+1
      do 10 j=1,nx
      ny=ny+1
      if(ny.le.iotm) iottr(ny)=jkan(ig)
      if (ig.lt.iwod) go to 10
      read (ltape) jkan
      ig=0
   10 ig=ig+1
      go to 8
    9 nconf(i)=jg
    7 continue
cdebug
*     do 322 ii1=1,iswhm
*        write(iwr,*) 'nconf(',ii1,')=',nconf(ii1)
*322  continue
      if(ofirst) then
       iotmf = ny
      endif
      if (ny.le.iotm) then
       write(iwr,51233) ny, iotm
51233  format(1x,'memory for configurations'/
     +        2x,'required (ny)    = ', i12/
     +        2x,'available (iotm) = ', i12)
      else
        write(iwr,51234) ny, iotm
51234   format(1x,'insufficient memory for configurations'/
     +       2x,'required (ny)    = ', i12/
     +       2x,'available (iotm) = ', i12/)
        call caserr(
     + 'insufficient memory for configurations - increase iotm')
      endif
      call rewftn(ltape)
cvp
cvp oeffnen der table
      call setbfc
cvp lesen des kopfes der table sowie der c- und der b-matrix
      call tbread(ifrk,ifr1,ntab)
*      write(iwr,*) 'from rumpxs after opening the table'
cvp
cvp aufbau von feldern, die konfigurationsvergleich erleichtern
cvp sollen
      call chfeld(iottr,iot0,iotm)
*      write(iwr,*) 'from rumpxs after chfeld'
cvp
cvp schreiben des ersten records von ft35
cvp
      m=nplu(1)+2*ndub(1)
c nsel wird zweimal geschrieben, um luecke aufzufuellen
      write(nhuk)iwod,vnuc,zero,imo,m,nconf,nytl,nplu,ndub,iswh,ksum
     1 ,iorbs,jsec,nrootx,vect,ew,ibal,itil,mcomp,kj,lj,n,
     & ifrk,knu,lsym,nsel,nsel,nj,ntil,
     & nbal(1),nbal(2),nbal(3),nbal(4),nbal(5),nbal(6),nbal(7),nbal(8),
     & ncomp

cbe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cbe         28.7.1994
cbe nbal wird vom ft31 mit 9 eintraegen uebernommen, im ft36 aber nur
cbe mit 8 eintraegen weitergegeben, daher werden jetzt hier nur
cbe die ersten 8 eintraege uebergeben, der neunte eintrag ist sowieso
cbe gleich ncomp(1) gesetzt worden, also ohne bedeutung
cbe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c    & ifrk,knu,e,lsym,inull,nsel,nj,ntil,nbal,ncomp
c    & ifrk,knu,e,lsym,nsel,inull,nj,ntil,nbal,ncomp
cvp
      read (ltape)
cvp
      do 22 i=1,iswh
      md=mconf(i)
cdebugwrite(iwr,*)'bei 22 : md=',md
      if (md.eq.0) go to 22
      nc=nconf(i)
cdebugwrite(iwr,*)'bei 22 : nc=',nc
      if (nc.ne.0) go to 836
      read (ltape)
      read (ltape)
      read (ltape)
      go to 22
  836 read (ltape) ndt,kml,ab
      write(nhuk) ndt,kml,ab
cdebugwrite(iwr,*)'ndt,kml,copied'
      read(ltape) khog,kmap,eb
      write(nhuk)khog,kmap,eb
cdebugwrite(iwr,*)'khog,kmap,copied'
      nx=nytl(i)*nc
      nx=nx/iwod+1
      do 900 j=1,nx
  900 read(ltape)
   22 continue
      call rewftn(ltape)
      jsum=0
*      write(iwr,*) 'from rumpxs before 901 '
      do 901 i=1,iswh
        inopen=nytl(i)-ndub(i)
        nbeta=inopen-nplu(i)
        if (nbeta.eq.0) then
           nsafsk=1
         else
           nsafsk=ibinom(inopen,nbeta)-ibinom(inopen,nbeta-1)
        endif
        ncons(i) = nsafsk
  901   jsum=jsum+nconf(i)*nsafsk
      jsum=(jsum-1)/ifrk+1
      write(nhuk)(b(i),i=1,5),isp
cold  write(nhuk)(b(i),i=1,5),isp,fzero
cold  fzero=5.0d0
cold  write(iwr,*) 'fzero for diagonal ',fzero
      do 902 i=1,jsum
      read (ideli) b
  902 write(nhuk)b
      read(ideli)nrootx,((trsum(j,i),j=1,10),i=1,nrootx),istm
      write(nhuk)((trsum(j,i),j=1,10),i=1,nrootx),istm
cbe abspeichern der energieerniedrigunen auf de0 bzw de1
******
****** some inconsidency here .. TPF
       if (iselecz(nrootx).eq.0) iselecz(nrootx)=nrootx
******
      do ialex = 1,nrootx
***** 
c       de0(ialex) = trsum(4,ialex)
c       de1(ialex) = trsum(5,ialex)
        de0(iselecz(ialex)) = trsum(4,iselecz(ialex))
        de1(iselecz(ialex)) = trsum(5,iselecz(ialex))
      enddo
      if (odebug) then
       write(6,*) 'rumpxs: nrootx, iselecz,de0,de1'
       write(6,*) 'nrootx = ', nrootx
       write(6,*) 'iselecz = ', iselecz
       write(6,*) 'de0 = ', de0
       write(6,*) 'de1 = ', de1
      endif
cbe abspeichern zu ende
      read(ltape)
      do 23 i=1,iswh
      md=mconf(i)
      if (md.eq.0) go to 23
      read(ltape)
      read(ltape)
      nc=nconf(i)
      if (nc.ne.0) go to 837
      read (ltape)
      go to 23
  837 nx=nytl(i)*nc
      nx=nx/iwod+1
      do 24 j=1,nx
      read(ltape) jkan
cdebugwrite(iwr,*) 'jkan written'
   24 write(nhuk) jkan
   23 continue
c --- lese die e/v
*      write(iwr,*) 'from rumpxs before rderz '
      call rderz(iotnew,iotm,maindf,maxr,jconb,nomax2,
     +           mtapev,iwr,ofirst,iotmf)
*      write(iwr,*) 'from rumpxs after rderz '
cvp
cvp
cvp  der konfigurationsvergleich aus main und beary wurde entfernt
cvp
cvp
cvp
cvp
cvp bestimmung und
cvp vorsortierung der integrale for ioint = 1
cvp  integraladressen werden auf iadl zwischengespeichert;
cvp  ist dieser vektor gefuellt, werden in isort die entsprechenden
cvp  integrale auf twoe abgespeichert; ist dieses feld gefuellt,
cvp  werden es auf mtype geschrieben inkl. der anzahl der integrale
cvp  und der nummer der entsprechenden konfiguration
*      write(iwr,*) 'issk=',issk
      if (ioint.eq.1) then
       if (issk.ne.0) then
        write(iwr,*) '1. call to preint'
        none = 1
        write(iwr,*) 'sk :',none,' bis sk: ',issk
        call cputim(cpu0)
*        write(iwr,*) 'from rumpxs before preint '
        call preint(twoe,nteint,ston,nidmax,
     +              pey, acoul, aexc, ndeke, core,
     +              pey0, acoul0, aexc0, noeint,
     +              iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              nrec31,none,issk)
*        write(iwr,*) 'from rumpxs after preint '
        call cputim(cpu1)
*        write(iwr,*) 'from rumpxs after cputim '
        dcpu=cpu1-cpu0
        write(iwr,100) dcpu
100     format(1x,'time for preint = ',f10.2,' secs.')
       endif
      endif
csut bilde matrixelemente 1,iisk
      call cputim(cpu0)
*        write(iwr,*) 'from rumpxs after cputim '
      if (issk.gt.0) then
*        write(iwr,*) 'from rumpxs before skinh '
       call skinh(twoe,nteint,
     +            pey, acoul, aexc, ndeke, core,
     +            emat,nedim,
     +            iottr,iotnew,iotm,
     +            maindf,maxr,jconb,nomax2,ioint,odebug)
*        write(iwr,*) 'from rumpxs before cputim '
       call cputim(cpu1)
       dcpu=cpu1-cpu0
       write(iwr,*) 'time for skinh :',dcpu,' s'
c --- anzahl matrixelemente
       none=icount+(irech-1)*ndims
       write(iwr,*) 'for ',none,' matrix elements'
       write(iwr,*) 'stored on ',(irech-1),' records'
      end if
*        write(iwr,*) 'from rumpxs before ioint.eq.1, ioint=',ioint
csut
cvp
      if (ioint.eq.1) then
       if (issk.lt.iswhm) then
        write(iwr,*) '2. call to preint'
        none = issk+1
        write(iwr,*) 'sk :',none,' bis sk: ',iswhm
        call cputim(cpu0)
        call rewftn(mtype)
        call preint(twoe,nteint,ston,nidmax,
     +              pey, acoul, aexc, ndeke, core,
     +              pey0, acoul0, aexc0, noeint,
     +              iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              nrec31,none,iswhm)
        call cputim(cpu1)
        dcpu=cpu1-cpu0
        write(iwr,100) dcpu
       endif
      endif
cvp hochsetzen von core <=> shift der diagonalelemente um fzero
clear core=core+fzero
csut  core=core
cvp
cvp berechnung des p=5-falles (insbes. der diagonale)
*      write(iwr,*) 'from rumpxs before hmp5 '
      call hmp5(pey,acoul,aexc,core,ndeke,hp5,mdi,iottr,
     +          iotm,iwr,ofirst,mdif)
cvp
*      write(iwr,*) ' ---------- rumpxs ended -----------'
      return
      end
cvp
cc
c berechnung der diagonal des sigma-vektors
      subroutine sdiag(vecf1,vecf2,mdi,hp5,mdi0)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
      real*8 vecf1,vecf2,hp5
      integer mdi,mdi0
cvp diagonalelemente der hamiltonmatrix (hp5, bisher nur diagonale)
      dimension vecf1(mdi),vecf2(mdi),hp5(mdi0)
cvp
      real*8 qq
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common /rwork2/ ndim,nanz
cvp
cdebug      write(iwr,*)
cdebug      write(iwr,*)'==============================='
cdebug      write(iwr,*)'= s d i a g   ================='
cdebug      write(iwr,*)'==============================='
cdebug      write(iwr,*)
cdebug      write(iwr,*) 'nanz=',nanz
cdebug      write(iwr,*) 'ndim=',ndim
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
cdebug        write(iwr,*) 'ipoint=',ipoint
c
         do 4012 ik=1,ndim
           qq=hp5(ik)
           vecf1(ipoint+ik)=vecf1(ipoint+ik)+vecf2(ipoint+ik)*qq
 4012    continue
cdebug         write(iwr,*) 'hdiag :  ',(hp5(ii),ii=1,10)
c ende schleife uber die wurzeln
 5100 continue
      return
      end
cvp
cvp berechnung der darstellungsmatrizen for dk=0
cvocl total,scalar
      subroutine sgdk0(emat,nedim,ae,ot,ie,kml,
     &                 ndt,inopen,ic,i2,ii,ia,j2,
     &                 j3,i3,nps,np1,nms,nm1,jc,np2,nm2)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
cvp
c
      real*8 f,ae
      common /rsghm4/f(ksafm*ksafm)
      dimension ae(7829)
cvp emat : darstellungsmatrizen
      real*8 emat
      integer nedim
      dimension emat(nedim)
cvp ot enthaelt die table-information bzgl. der determinanten
      integer ot(motmax)
c aufbau der table-faktoren for r- und z-faelle:
c  ibob2 -> coulomb-integrale, jbob2 -> exchange-integrale
c  (spaeter evtl. als data statement)
      ibob2(1,1)= 1
      ibob2(1,2)= 0
      ibob2(1,3)= 1
      ibob2(2,1)= 0
      ibob2(2,2)=-1
      ibob2(2,3)= 1
      ibob2(3,1)= 0
      ibob2(3,2)=-1
      ibob2(3,3)= 1
      jbob2(1,1)= 1
      jbob2(1,2)= 1
      jbob2(1,3)= 0
      jbob2(2,1)= 1
      jbob2(2,2)= 0
      jbob2(2,3)= 1
      jbob2(3,1)=-1
      jbob2(3,2)=-1
      jbob2(3,3)= 0
cvp
      mzer=kml*kml
cvp  nqr = # qr-faelle
cvp  nql = # ql-faelle
cvp  irec = # eintraege von emat for ein saf-paar
      nqr=ibinom(inopen,2)
      nql=nqr
      irec=6*ideks(nqr+1)
cvp
cvp schleife ueber alle r-faelle
      do 1000 l=1,3
cvp schleife ueber alle ql-faelle
      do 1100 jjql=1,nql
cvp schleife ueber alle qr-faelle
      do 1200 jjqr=jjql,nqr
cvp
      iemi= 2*(jjqr-jjql)+1 + 2*(jjql-1)*(nqr+1)
     &      + 2*(l-1)*ideks(nqr+1) - 2*ideks(jjql)
      iemj= 2*(jjqr-jjql)+2 + 2*(jjql-1)*(nqr+1)
     &      + 2*(l-1)*ideks(nqr+1) - 2*ideks(jjql)
cvp
      jj=ideks(jjqr)+jjql
      jj=jj*ii+ic
c die allgemeine permutation
      ip=ot(jj)
      if(ip.eq.255) ip=-1
      do 154 m=1,mzer
  154 f(m)=0.0d0
      mx=-kml
      do 155 m=1,kml
        jj=jj+1
        mx=mx+1
        ip1=ot(jj)
        if(ip1.eq.2)  go to 159
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        jj=jj+1
        tim=dble(ip*ibob2(l,1))
        lx=mx
        do 161 ih=1,kml
          kx=kx+ndt
          lx=lx+kml
c161      f(lx)=f(lx)+tim*ae(kx)
  161     f(lx)=      tim*ae(kx)
        go to 155
  159   jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=dble(ip*ibob2(l,2))
        lx=mx
        do 160 ih=1,kml
          kx=kx+ndt
          lx=lx+kml
c160      f(lx)=f(lx)+tim*ae(kx)
  160     f(lx)=      tim*ae(kx)
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=dble(ip*ibob2(l,3))
        lx=mx
        do 167 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
  167     f(lx)=f(lx)+tim*ae(kx)
  155 continue
c     if (l.eq.1.and.jjqr.eq.3.and.jjql.eq.3) then
c        write(iwr,*) ' f= ',(f(i70),i70=1,mzer)
c     endif
cvp endtransformation for dk=0
      do 145 m=1,kml
        do 145 ih=1,kml
          tim=0.0d0
          do 143 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  143       tim=tim+f(ly)*ae(mx)
          iem=iemi+(ih-1)*irec+(m-1)*irec*kml
          emat(iem)=tim
  145 continue
cvp
cvp darstellungsmatrix for ex-integral
      jj=ideks(jjqr)+jjql
      jj=jj*ii+ic
c die allgemeine permutation
      ip=ot(jj)
      if(ip.eq.255) ip=-1
      do 254 m=1,mzer
  254 f(m)=0.0d0
      mx=-kml
      do 255 m=1,kml
        jj=jj+1
        mx=mx+1
        ip1=ot(jj)
        if(ip1.eq.2)  go to 259
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        jj=jj+1
        tim=dble(ip*jbob2(l,1))
        lx=mx
        do 261 ih=1,kml
          kx=kx+ndt
          lx=lx+kml
c161      f(lx)=f(lx)+tim*ae(kx)
  261     f(lx)=      tim*ae(kx)
        go to 255
  259   jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=dble(ip*jbob2(l,2))
        lx=mx
        do 260 ih=1,kml
          kx=kx+ndt
          lx=lx+kml
c160      f(lx)=f(lx)+tim*ae(kx)
  260     f(lx)=      tim*ae(kx)
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=dble(ip*jbob2(l,3))
        lx=mx
        do 267 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
  267     f(lx)=f(lx)+tim*ae(kx)
  255 continue
cvp endtransformation for dk=0
      do 245 m=1,kml
        do 245 ih=1,kml
          tim=0.0d0
          do 243 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  243       tim=tim+f(ly)*ae(mx)
          iem=iemj+(ih-1)*irec+(m-1)*irec*kml
          emat(iem)=tim
c     if (l.eq.1.and.jjqr.eq.3.and.jjql.eq.3) then
c        write(6,*) ' m=,ih=,iem=,tim= ',m,ih,iem,tim
c     endif
  245 continue
cvp
 1200 continue
cvp
 1100 continue
cvp
 1000 continue
cvp
cvp aufbau der darstellungsmatrizen for p=2
cvp
      nqr2=inopen
      nql2=inopen
      irec2=ideks(nqr2+1)
cvp iema = anzahl der eintraege in emat for p=1
      iema=irec*kml*kml
cvp schleife ueber alle r-faelle
cvp   entfaellt, da keine unterschiede zwischen r=1 und r=2
cvp schleife ueber alle qr-faelle
      do 2000 jjql=1,nql2
cvp schleife ueber alle qr-faelle
      do 2100 jjqr=jjql,nqr2
cvp
      iem2= iema + (jjql-1)*(nqr2+1)-ideks(jjql) + jjqr-jjql+1
cvp
      jj=ideks(jjqr)+jjql
      jj=jj*j2+i2
      do 177 m=1,mzer
  177 f(m)=0.0d0
      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) then
         tim=-1.0d0
       else
         tim= 1.0d0
      endif
      do 173 m=1,kml
        jx=jx+1
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        lx=jx
        do 173 ih=1,kml
          kx=kx+ndt
          lx=lx+kml
c 173     f(lx)=f(lx)+tim*ae(kx)
  173     f(lx)=      tim*ae(kx)
c beginn der endtransformation von dk=0 (for alle p-faelle gleich)
cvp endtransformation for dk=0
      do 345 m=1,kml
        do 345 ih=1,kml
          tim=0.0d0
          do 346 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  346       tim=tim+f(ly)*ae(mx)
          iem=iem2+(ih-1)*irec2+(m-1)*irec2*kml
          emat(iem)=tim
  345 continue
cvp
 2100 continue
cvp
 2000 continue
cvp
cvp aufbau der darstellungsmatrizen for p=3
cvp
      nqr3=inopen
      nql3=nqr3
      irec3=ideks(inopen+1)*2*inopen
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
      iemb= iema + irec2*kml*kml
cvp schleife ueber alle einzelnen summanden
      do 3000 iis=1,inopen
cvp schleife ueber alle r-faelle
      do 3100 l=1,2
cvp schleife ueber alle ql-faelle
      do 3200 jjql=1,nql3
cvp schleife ueber alle qr-faelle
      do 3300 jjqr=jjql,nqr3
cvp
      iem3= iemb + (iis-1)*2*ideks(inopen+1) + (l-1)*ideks(inopen+1)
     &           + (jjql-1)*(nqr3+1)-ideks(jjql+1) + jjqr+1
cvp setzen der integrale, um aus der table die darstellungsmatrizen
cvp  zu erhalten
      sm=0.0d0
      do 66 m=1,inopen-1
        sac(m)=0.0d0
   66 continue
cvp
      if (iis.eq.1) then
         sm=1.0d0
       else
         sac(iis-1)=1.0d0
      endif
c abzweigung zu den unterschiedlichen r-faellen
      if (l.eq.2) go to 181
c dk=0 p=3 r=1
      jj=ideks(jjqr)+jjql
      jj=jj*j3+i3
cvp
      ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      do 192 m=1,inopen-1
  192 sac(m)=-sac(m)
  191 jj=jj+1
      do 193 m=1,mzer
  193 f(m)=0.0d0
      jx=-kml
      do 194 m=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=sm
        if(im.eq.2) go to 195
        if(nps.eq.1) go to 196
        do 197 ih=1,np1
          in=in+2
          nn=ot(in)
  197     tim=tim+sac(nn)
  196   do 198 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
c198      f(lx)=f(lx)+tim*ae(kx)
  198     f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 194
        do 199 ih=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 199 ig=1,kml
            lx=lx+kml
            kx=kx+ndt
  199       f(lx)=f(lx)+tim*ae(kx)
        go to 194
  195   if(nms.eq.1) go to 200
        do 201 ih=1,nm1
          in=in+2
          nn=ot(in)
  201     tim=tim+sac(nn)
  200   do 202 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
  202     f(lx)=f(lx)+tim*ae(kx)
        do 203 ih=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 203 ig=1,kml
            lx=lx+kml
            kx=kx+ndt
  203       f(lx)=f(lx)+tim*ae(kx)
  194   jj=jj+jc
      go to 146
c dk=0 p=3 r=2
  181 continue
      jj=ideks(jjqr)+jjql
      jj=jj*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      do 216 m=1,inopen-1
  216 sac(m)=-sac(m)
  215 jj=jj+1
      do 217 m=1,mzer
  217 f(m)=0.0d0
      jx=-kml
      do 218 m=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=-sm
        if (im.eq.2) go to 219
        if(nms.eq.0) go to 220
        in=in+np2
        jn=in
        do 221 ih=1,nms
          jn=jn+2
          nn=ot(jn)
  221     tim=tim-sac(nn)
  220   do 222 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
c 222     f(lx)=f(lx)+tim*ae(kx)
  222     f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 218
        do 223 ih=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 223 ig =1,kml
            lx=lx+kml
            kx=kx+ndt
  223       f(lx)=f(lx)+tim*ae(kx)
        go to 218
  219   in=in+nm2
        jn=in
        do 224 ih=1,nps
          jn=jn+2
          nn=ot(jn)
  224     tim=tim-sac(nn)
        do 225 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
  225     f(lx)=f(lx)+tim*ae(kx)
        do 226 ih=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 226 ig=1,kml
            lx=lx+kml
            kx=kx+ndt
  226       f(lx)=f(lx)+tim*ae(kx)
  218   jj=jj+jc
cvp
  146 continue
c beginn der endtransformation von dk=0 (for alle p-faelle gleich)
cvp endtransformation for dk=0
      do 445 m=1,kml
        do 445 ih=1,kml
          tim=0.0d0
          do 446 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  446       tim=tim+f(ly)*ae(mx)
          iem=iem3+(ih-1)*irec3+(m-1)*irec3*kml
          emat(iem)=tim
  445 continue
cvp
 3300 continue
cvp
 3200 continue
cvp
 3100 continue
cvp
 3000 continue
c     do 7 iiii=1,36
c  7  write(6,*) ' iiii=,emat= ',iiii,emat(iiii)
cvp
      return
      end
cvp
cvp berechnung der darstellungsmatrizen for dk=1
cvocl total,scalar
      subroutine sgdk1(emat,nedim,ae,ot,ja,ie,kml,kmj
     &                ,ndj,inopen,ic,i2
     &                ,j3,i3,mps,mms,nms,nps,jc)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      real*8 f,ae
      common /rsghm4/f(ksafm*ksafm)
      dimension ae(7829)
cvp emat : darstellungsmatrizen
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
cvp ot enthaelt die table-information bzgl. der determinanten
      integer ot(motmax)
c aufbau der table-faktoren for r- und z-faelle:
c  ibob2 -> coulomb-integrale, jbob2 -> exchange-integrale
c  (spaeter evtl. als data statement)
      ibob2(1,1)= 0
      ibob2(1,2)=-1
      ibob2(1,3)= 1
      ibob2(2,1)= 0
      ibob2(2,2)=-1
      ibob2(2,3)= 1
      ibob2(3,1)= 1
      ibob2(3,2)= 0
      ibob2(3,3)=-1
      jbob2(1,1)=-1
      jbob2(1,2)= 0
      jbob2(1,3)= 1
      jbob2(2,1)= 1
      jbob2(2,2)=-1
      jbob2(2,3)= 0
      jbob2(3,1)= 1
      jbob2(3,2)=-1
      jbob2(3,3)= 0
cvp
      ii=2*kml+1
      nzer=kml*kmj
cvp  nqr = # qr-faelle
cvp  nql = # ql-faelle
cvp  irec = # eintraege von emat for ein saf-paar
      nqr=ibinom(inopen,3)
      nql=inopen-2
      irec=6*nqr*nql
cvp
cvp schleife ueber alle r-faelle
      do 1000 l=1,3
cvp schleife ueber alle ql-faelle
      do 1100 jjql=1,nql
cvp schleife ueber alle qr-faelle
      do 1200 jjqr=1,nqr
cvp
      iemi= 2*jjqr-1 + 2*(jjql-1)*nqr + 2*(l-1)*nqr*nql
      iemj= 2*jjqr   + 2*(jjql-1)*nqr + 2*(l-1)*nqr*nql
cvp
cvp darstellungsmatrix for cb-integral
      jj=(jjqr-1)*(inopen-2) +jjql
      jj=jj*ii+ic
*****
c     jj2=jj*ii+ic
c     write(6,921) jjqr, inopen, jjql, jj, ii, ic, jj2
c921  format(' jjqr, inopen, jjql, jj, ii, ic, jj2', 7i8)
*****
      ip=ot(jj)
      if (ip.eq.255) ip=-1
      do 39 m=1,nzer
   39 f(m)=0.0d0
      do 40 m=1,kml
        jj=jj+1
        my=ot(jj)
        jj=jj+1
        if(my.ne.0) then
          kx=ja+my
          mike=ot(jj)
          if(mike.gt.128) then
             tim=-dble(ip*ibob2(l,256-mike))
           else
             tim=dble(ip*ibob2(l,mike))
          endif
          lx=m-kml
          do 44 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 44        f(lx)=f(lx)+tim*ae(kx)
   44       f(lx)=      tim*ae(kx)
        endif
   40 continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      do 45 m=1,kml
        do 45 ih=1,kmj
          tim=0.0d0
          do 46 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
   46       tim=tim+f(ly)*ae(mx)
          iem=iemi+(ih-1)*irec+(m-1)*irec*kmj
          emat(iem)=tim
   45 continue
cvp
cvp darstellungsmatrix for ex-integral
      jj=(jjqr-1)*(inopen-2) +jjql
      jj=jj*ii+ic
      ip=ot(jj)
      if (ip.eq.255) ip=-1
      do 139 m=1,nzer
  139 f(m)=0.0d0
      do 140 m=1,kml
        jj=jj+1
        my=ot(jj)
        jj=jj+1
        if(my.ne.0) then
          kx=ja+my
          mike=ot(jj)
          if(mike.gt.128) then
             tim=-dble(ip*jbob2(l,256-mike))
           else
             tim=dble(ip*jbob2(l,mike))
          endif
          lx=m-kml
          do 144 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 44        f(lx)=f(lx)+tim*ae(kx)
  144       f(lx)=      tim*ae(kx)
        endif
  140 continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      do 145 m=1,kml
        do 145 ih=1,kmj
          tim=0.0d0
          do 146 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  146       tim=tim+f(ly)*ae(mx)
          iem=iemj+(ih-1)*irec+(m-1)*irec*kmj
          emat(iem)=tim
  145 continue
cvp
 1200 continue
cvp
 1100 continue
cvp
 1000 continue
cvp
cvp aufbau der darstellungsmatrizen for p=2
cvp
      nqr2=ibinom(inopen,2)
      irec2=nqr2
cvp iema = anzahl der eintraege in emat for p=1
      iema=irec*kml*kmj
cvp schleife ueber alle r-faelle
cvp   entfaellt, da sich r=1 und r=2 im nur vorzeichen der
cvp   matrixelemente unterscheiden
cvp schleife ueber alle qr-faelle
      do 2100 jjqr=1,nqr2
cvp
      iem2= jjqr + iema
cvp
      jj=jjqr*kml+i2
      do 55 m=1,nzer
   55 f(m)=0.0d0
      do 56 m=1,kml
        jj=jj+1
c welche determinante ww mit der m-ten sdet, vorzeichen kodiert ip
        my=ot(jj)
        if (my.ne.0) then
          if(my.lt.128) then
             tim=1.0d0
             kx=my+ja
           else
             kx=ja-my+256
             tim=-1.0d0
          endif
          lx=m-kml
          do 60 if=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 60        f(lx)=f(lx)+tim*ae(kx)
   60       f(lx)=      tim*ae(kx)
        endif
   56 continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      do 245 m=1,kml
        do 245 ih=1,kmj
          tim=0.0d0
          do 246 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  246       tim=tim+f(ly)*ae(mx)
          iem=iem2+(ih-1)*irec2+(m-1)*irec2*kmj
          emat(iem)=tim
  245 continue
cvp
 2100 continue
cvp
c2000 continue
cvp
cvp aufbau der darstellungsmatrizen for p=3
cvp
      nqr3=ibinom(inopen,2)
      irec3=(inopen-1)*2*nqr3
cvp iemb = anzahl der eintraege in emat for p=1 und p=2
      iemb= iema + irec2*kml*kmj
cvp schleife ueber alle einzelnen summanden
      do 3000 iis=1,inopen-1
cvp schleife ueber alle r-faelle
      do 3100 l=1,2
cvp schleife ueber alle qr-faelle
      do 3200 jjqr=1,nqr3
cvp
      iem3= iemb + (iis-1)*2*nqr3 + (l-1)*nqr3 + jjqr
cvp setzen der integrale, um aus der table, die darstellungsmatrizen
cvp  zu erhalten
      sm=0.0d0
      do 66 m=1,inopen-2
        sac(m)=0.0d0
   66 continue
cvp
      if (iis.eq.1) then
         sm=1.0d0
       else
         sac(iis-1)=1.0d0
      endif
c abzweigung r=1 bzw 2
cvp qr < 0 --> r=1; qr > 0 --> r=2
      if (l.eq.1) go to 86
      jj=jjqr*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 71
      sm=-sm
cvp inopen-2 = # der offenen schalen der kleineren sk
      do 72 m=1,inopen-2
   72 sac(m)=-sac(m)
   71 jj=jj+1
      do 80 m=1,nzer
   80 f(m)=0.0d0
      jx=-kml
      do 73 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (74,75,76,77),im
   74    in=in+1
         kx=ot(in)
         kx=kx+ja
         tim=sm
cvp nps = # alpha-spins der offenen schalen der kleineren sk
         do 79 ih=1,mps
           in=in+1
           im=ot(in)
   79      tim=tim+sac(im)
   78    lx=jx
         do 81 ih=1,kmj
           kx=kx+ndj
           lx=lx+kml
   81      f(lx)=f(lx)+tim*ae(kx)
        go to 73
   75   in=in+1
        kx=ot(in)
        kx=kx+ja
        tim=-sm
cvp mms = nr. der kleineren sk - 1 ; bzw.
cvp mms = #  beta-spins der offenen schalen der kleineren sk
        do 83 ih=1,mms
          in=in+1
          im=ot(in)
   83     tim=tim-sac(im)
        go to 78
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
   76   do 84 ih=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=-sac(im)
          lx=jx
          do 84 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
   84       f(lx)=f(lx)+tim*ae(kx)
        go to 73
cvp nps = # alpha-spins der offenen schalen
   77   do 85 ih=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 85 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
   85       f(lx)=f(lx)+tim*ae(kx)
   73 jj=jj+jc
      go to 47
c p=3 anderer r fall (r=1)
   86 jj=jjqr*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      do 88 m=1,inopen-2
   88 sac(m)=-sac(m)
   87 jj=jj+1
      do 89 m=1,nzer
   89 f(m)=0.0d0
      jx=-kml
      do 90 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (91,95,96,97),im
   91     in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=sm
          if(mms.eq.0) go to 92
          in=in+mps
          do 93 ih=1,mms
            in=in+1
            im=ot(in)
   93       tim=tim+sac(im)
   92     lx=jx
          do 94 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
   94       f(lx)=f(lx)+tim*ae(kx)
          go to 90
   95     in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=-sm
          if(mps.eq.0) go to 92
          in=in+mms
          do 99 ih=1,mps
            in=in+1
            im=ot(in)
   99       tim=tim-sac(im)
          go to 92
   96     do 100 ih=1,nms
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 100 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
  100       f(lx)=f(lx)+tim*ae(kx)
          go to 90
   97     do 101 ih=1,nps
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
            im=ot(in)
            tim=-sac(im)
            lx=jx
            do 101 ig=1,kmj
               kx=kx+ndj
               lx=lx+kml
  101          f(lx)=f(lx)+tim*ae(kx)
   90      jj=jj+jc
c
   47 continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      do 345 m=1,kml
        do 345 ih=1,kmj
          tim=0.0d0
          do 346 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
  346       tim=tim+f(ly)*ae(mx)
          iem=iem3+(ih-1)*irec3+(m-1)*irec3*kmj
          emat(iem)=tim
  345 continue
cvp
 3200 continue
cvp
 3100 continue
cvp
 3000 continue
c     write(6,*) ' emat= ',(emat(iiii),iiii=1,150)
cvp
      return
      end
cvp
cvp berechnung der darstellungsmatrizen for dk=2
cvocl total,scalar
      subroutine sgdk2(emat,nedim,ae,ja,ie,kml,kmj
     &                ,ndj,inopen)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      real*8 f,ae
      common /rsghm4/f(ksafm*ksafm)
      dimension ae(7829)
cvp emat : darstellungsmatrizen
      real*8 emat
      integer nedim
      dimension emat(nedim)
cvp
      ii=2*kml+1
      nis=1-ii
      nzer=kml*kmj
cvp  nqr = # qr-faelle
cvp  irec = # eintraege von emat for ein saf-paar
      nqr=ibinom(inopen,4)
      irec=6*nqr
cvp
cvp schleife ueber alle r-faelle
      do 1000 l=1,3
cvp schleife ueber alle qr-faelle
      do 1100 jjqr=1,nqr
cvp
      iemi=2*(jjqr-1)+1 + 2*(l-1)*nqr
      iemj=2*(jjqr-1)+2 + 2*(l-1)*nqr
cvp
      jj=jjqr*ii+nis
      ip=ntabl(jj)
cvp
cvp darstellungsmatrix for cb-integral
      do 905 m=1,nzer
         f(m)=0.0d0
 905  continue
      do 18 m=1,kml
        my=ntabl(jj+m+m-1)
        if (my.ne.0) then
          kx=ja+my
          mike=ntabl(jj+m+m)
          tim=dble(ip*ibob2(l,mike))
          lx=m-kml
          do 22 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
   22       f(lx)=      tim*ae(kx)
c 22      f(lx)=f(lx)+tim*ae(kx)
c
      endif
   18 continue
cvp
      do 23 m=1,kml
        do 23 ih=1,kmj
          tim=0.0d0
          do 24 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
   24       tim=tim+f(ly)*ae(mx)
          iem=iemi+(ih-1)*irec+(m-1)*irec*kmj
          emat(iem)=tim
   23 continue
cvp
cvp darstellungsmatrix for ex-integral
      do 1905 m=1,nzer
1905  f(m)=0.0d0
      do 118 m=1,kml
        my=ntabl(jj+m+m-1)
        if (my.ne.0) then
          kx=ja+my
          mike=ntabl(jj+m+m)
          tim=dble(ip*jbob2(l,mike))
          lx=m-kml
          do 122 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
  122       f(lx)=      tim*ae(kx)
c 22      f(lx)=f(lx)+tim*ae(kx)
c
      endif
  118 continue
cvp
      do 123 m=1,kml
        do 123 ih=1,kmj
          tim=0.0d0
          do 124 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
            tim=tim+f(ly)*ae(mx)
  124 continue
          iem=iemj+(ih-1)*irec+(m-1)*irec*kmj
          emat(iem)=tim
  123 continue
cvp
 1100 continue
cvp
 1000 continue
cvp
c     write(6,*) ' emat= ',(emat(iiii),iiii=1,150)
cvp
      return
      end
cvp
cc
c     sigmavektors aus den gespeicherten hamiltonmatrixelementen
      subroutine sham(vecf1,vecf2,mdi)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
      real*8 vecf1, vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
cvp
csut /speicher
      integer icount,irech,imm,jmm
      integer mmi,mmj
      real*8 rham,rhams
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
      common /junk3/ imms(ndims),jmms(ndims),rhams(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
cbe der erste block steht auf dem common/rhsc/ (imm,jmmm,rham)
cbe und alle anderen
cbe bloecke falls notwendig werden auf extra speicher
cbe geladen (imms,jmms,rhams)
cvp
c common for multi-root version
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common /rwork2/ ndim,nanz
cvp
      call rewftn(nf88)
*      write(6,*) 'sham  : irech2',irech2
*      write(6,*) 'sham  : irech=',irech
*      write(6,*) 'sham  : ncount2',ncout2
*      write(6,*) 'sham  : ncount',ncout
       do 5200 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
cdebug  write(6,*)' 5200 : ipoint=',ipoint
        if (icount.gt.0) then
         do ii=1,icount
          mmi=imm(ii)
          mmj=jmm(ii)
          qq=rham(ii)
          vecf1(ipoint+mmi)=vecf1(ipoint+mmi)+vecf2(ipoint+mmj)*qq
          vecf1(ipoint+mmj)=vecf1(ipoint+mmj)+vecf2(ipoint+mmi)*qq
         end do
        end if
 5200  continue
*          write(6,*) 'from sham after first read'
*          write(6,*) 'treated were '
*          write(6,*) 'imm',(imm(iaus),iaus=1,10)
*          write(6,*) 'jmm',(jmm(iaus),iaus=1,10)
*          write(6,*) (rham(iaus),iaus=1,10)
*          write(6,*) 'vecf1',(vecf1(iii),iii=1,10)

      if (irech.gt.1) then
       do iread=1,(irech-1)
        read(nf88)imms,jmms,rhams
csut    write(6,*) '88 read'
        do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
          ipoint=(iwurz-1)*ndim
c         write(6,*)' 5100 : ipoint=',ipoint
          do ii=1,ndims
           mmi=imms(ii)
           mmj=jmms(ii)
           qq=rhams(ii)
           vecf1(ipoint+mmi)=vecf1(ipoint+mmi)+vecf2(ipoint+mmj)*qq
           vecf1(ipoint+mmj)=vecf1(ipoint+mmj)+vecf2(ipoint+mmi)*qq
          end do
c
 5100   continue
*          write(6,*) 'from sham after ',iread,'+1-th read'
*          write(6,*) 'vecf1',(vecf1(iii),iii=1,10)
       end do
      end if
c
*     write(6,*) 'rham(1-10) at end of sham'
*     write(6,*) (rham(iaus),iaus=1,10)

      return
      end
cvp
cc
c berechnung der matrixelemente h(i,j)
csut : abspeicherung der matrixelemete auf fort.88
      subroutine skinh(twoe,nteint,
     +                 pey, acoul, aexc, ndeke, core,
     +                 emat,nedim,
     +                 iottr,iotnew,iotm,
     +                 maindf,maxr,jconb,nomax2,ioint,odebug)
csut
csut  icount : zaehler ueber die matrixelemente
csut  irech  : zaehler der schon beschriebenen rekords
csut  ndims  : groesse der speicherrekorde
csut  rekordaufbau in fort.88:
csut  do i=1,irech
csut     imm(ndims),jmm(ndims),h(ndims)
csut  end do
csut  der rest 1..icount wird im memory gehalten!
csut
csut
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
      integer nteint, ndeke
      real*8 twoe
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      logical odebug
c
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      real*8 hstor
c     integer irech2,ncount
c     common /crech/irech2,ncount
      common /csemi/ hstor,issk
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common /rwork2/ ndim,nanz
cvp
cvp table information: insbes. c- und b-matrix
      real*8 af, ae
      common /rtable/ af(8000)
      dimension ae(7829)
      integer nsac,ndet,iaz,iez,jan,jbn,idra,idrc
      integer jdrc,j9
      dimension nsac(iswhm),ndet(iswhm),
     +  iaz(iswhm),iez(iswhm),idrc(iswhm),
     +  idra(iswhm),jdrc(iswhm),jan(5),jbn(5),
     +  j9(49)
      equivalence (af(50),ae(1))
      equivalence (af(1),j9(1)),(j9(1),ndet(1)),(j9(6),
     +nsac(1)),(j9(11),iaz(1)),(j9(16),iez(1)),(j9(21),idra(1)),
     +(j9(26),idrc(1)),(j9(31),jdrc(1)),(j9(36),jan(1)),(j9(43),
     +jbn(1))
cvp
      common/linkmr/ic,i2,i3
      common/blksi3/nsz
      common/bufd/gout(510),nword
      dimension lout(510)
      equivalence (lout(1),gout(1))
      integer ot
      common /rwork/ ot(motmax)
      common /ra/ bob(3), f(ksafm*ksafm)
      integer out, ohog
      dimension ohog(48),out(10962)
      equivalence (ot(13),ohog(1)),(ot(13),out(1))
c
      common /rb/ ifrk, iswh
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
c
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
      integer isc
      dimension isc(iswhm)
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
      integer ijone,iwod,mconf,nplu,idummy
      common /rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
cvp
cvp entwicklungskoeffizienten von saf's in slater-determinanten
cvp   dimension cmat(ncdim),emat(nedim),istvec(iswhm)
cvp emat : darstellungsmatrizen
cvp
c
c bobvec=ww-terme for alle konf. aus sk(i-2) und alle z-faelle,
      real*8 bobvec
      common/scra/bobvec(ncmax,3)
c
      real*8 tim1
      integer inull
      data inull/0/
c
c beginn timer-routine for berechnung der matrixelemente h(k,l)
c -- behalten
*      write(6,*) 'call skinh'
      jtape = 6
      millia=0
csut leeren der arrays
      do i=1,ndims
        imm(i)  = 0
        jmm(i)  = 0
        rham(i) = 0.0d0
      end do
      icount = 0
      irech  = 1
c
cdebug
*      write(6,*) ' skiny: ioint= ',ioint
cdebugwrite(6,*) ' skiny: ndim = ',ndim
cdebugwrite(6,*) ' skiny: nanz = ',nanz
cdebugwrite(6,*) ' skiny: vecf1  '
cdebugwrite(6,*) (vecf1(iii),iii=1,5)
cdebugwrite(6,*) ' skiny: vecf2  '
cdebugwrite(6,*) (vecf2(iii),iii=1,5)
cdebug
      ifr1=1-ifrk
csut  nad=4000
csut  nad1=-3999
      call rewftn(nf88)
c
      is=nplu(1)
      mult=is+1
cvp
cvp erzeugung der matrix for entwicklung einer saf in slater-det.
cvp  mult = 2*s +1
cvp   call saf1(mult,cmat,istvec)
cvp
cvp  test: umspeichern von buenkers c-matrix auf cmat
c     do 4713 i=1,iswh
c        nndet=ndet(i)
c        kmml=nsac(i)
c        istb=iaz(i)
c        do 4712 iiii1=1,kmml
c        do 4712 iiii2=1,nndet
c4712    cmat(istvec(i)+(iiii1-1)*nndet+iiii2)=
c    &       ae(istb+iiii1*nndet+iiii2)
c4713 continue
cvp  berechnung der dimension der matrix bis zur sk i-1
      jsum=0
csut  write(6,*) 'skinh : issk=',issk
*      write(6,*) 'skinh : hstor=',hstor
      do 901 i=1,iswh
csut  do 901 i=1,issk
      isc(i)=jsum
  901 jsum=jsum+nconf(i)*nsac(i)
cvp
*      write(6,*) 'after do 901'
cvp
cvp
      ibl=0
      jg=0
cvp
cvp lesen der integrale, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
*      write(6,*) 'ioint ',ioint
      if (ioint.eq.1) then
cdebug   call rewftn(mtype)
cdebug   read(mtype) iges
*         write(6,*) 'iges=',iges
cdebug   call rewftn(mtype)
         read (mtype) iges,kend,(twoe(iii),iii=1,iges)
         jto1=0
      endif
*      write(6,*) 'after reading of integrals at ioint'
      ipoi2t = 0
cvp
c statement for rumple test
c zaehler for konfigurationspaare
csut  do 6 i=1,iswh
*      write(6,*) 'issk before do 6',issk
       call setsto(motmax,0,ot)
      do 6 i=1,issk
c write statement eingefuegt
*      write(6,500) i
* 500 format(/50x,'*******   begin of sk  ',i2,'  *******'/)
      nc=nconf(i)
      if (nc.gt.0) go to 7
csut  md=mconf(i)
      if (i.lt.3) go to 8
      ibl=ibl+2
      go to 6
8     ibl=ibl+i-1
      go to 6
7     ndt=ndet(i)
c//skipc&2
c     call cpu
c     call setcpu
      ia=iaz(i)
      kml=nsac(i)
      ndl=ndet(i)
      ie=iez(i)
      iz=isc(i)
      mzer=kml*kml
      nd=ndub(i)
      nps=nplu(i)
      nms=i-1
      niw=nps*nms
      nr=nod(i)
      nx=nytl(i)
      n1=nr+1
      n2=nr+2
      n3=nr-1
      np1=nps-1
      np2=np1+np1
      nd1=nd-1
      nm1=i-2
      nm2=nm1+nm1
      ii=kml+kml+1
      nis=1-ii
      inopen=nod(i)
c write allgemein
c     write(iwr,2000) i,nc,ndl,kml,nd,nr,nps
c2000 format(/1x,'(2000)*****  sk=',i2,' # konf =',i5,' # det =',
c    *i3,'# sdet = ',i2,'# cl.sh. =',i2,' # op.sh. =',i2,
c    *' # alpha s = ',i2/)
      if (i.lt.3) go to 9
      ibl=ibl+1
c start dk=2 fall
c write statement eingefuegt
*      write(6,501)
* 501 format(/20x,'*******   start dk=2   *******'/)
      j8=i-2
      mc=nconf(j8)
      if(mc.eq.0) go to 10
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      js=jan(ibl)
      ks=jbn(ibl)
cvp
c     write(iwr,2001) j8,mc,ndj,kmj,nr,js,ks
c2001 format(1x,'dk=2 (2001) j8,mc,ndj,kmj,nr,js,ks',7(i5,1x))
      ix=nad1
      kss=0
 11   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 11
*****
      call upackx(ot)
      endif
cvp
cvp konvertierung von ot auf das integer feld ntabl
      call tabcon(inopen,kml,ot)
cvp
cvp
cvp erzeugung der tabelle for die berechnung der matrixelemente
c     call matdk2(cmat,emat,istvec,i,mult,kml,kmj,inopen)
cvp
      kz=iz
c beginn timer-routine for dk=2
*     call date(datum)
*     call time(milli)
*     isek=milli/1000
*     imin=isek/60
*     istu=imin/60
*     imin=imin-istu*60
*     write(6,*) 'calculation were performed : ',datum,' time ',
*    *istu,':',imin
cdebugcall clockv(vu1,cpu1)
cdebugcall time(milli1)
c aufbau der table-faktoren for r- und z-faelle:
c  ibob2 -> coulomb-integrale, jbob2 -> exchange-integrale
c  (spaeter evtl. als data statement)
      ibob2(1,1)=-1
      ibob2(1,2)=0
      ibob2(1,3)=1
      ibob2(2,1)=0
      ibob2(2,2)=1
      ibob2(2,3)=-1
      ibob2(3,1)=0
      ibob2(3,2)=1
      ibob2(3,3)=-1
      jbob2(1,1)=-1
      jbob2(1,2)=1
      jbob2(1,3)=0
      jbob2(2,1)=-1
      jbob2(2,2)=1
      jbob2(2,3)=0
      jbob2(3,1)=1
      jbob2(3,2)=0
      jbob2(3,3)=-1
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,4)
      if (npdim.le.nedim) then
       if (odebug) then
*        write(6,12345) npdim, nedim
*12345 format(' ***calculation of matrix elements for dk=2 via sga '/
*    +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call sgdk2(emat,nedim,ae,ja,ie,kml,kmj
     &             ,ndj,inopen)
         call hmdk2h(twoe,emat,nedim,iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              i,nc,mc,kml,kmj,
     +              inopen,jz,kz
     &             ,jto1,ioint,iges,kend)
         goto 1014
      else
       if (odebug) then
       write(6,12346) npdim, nedim
12346  format(' **** sgdk2:  npdim.gt.nedim', 2i7)
       endif
      endif
c schleife ueber alle konfigurationen
cdebugwrite(6,*)' nc=',nc
cdebugwrite(6,*)' mc=',mc
cdebugwrite(6,*)' kml=',kml
cdebugwrite(6,*)' kmj=',kmj
c iwwg zaehlt alle ww konfigurationen for dk=2
c  iwwz=zaehler aller von null verschiedenen iww's
      iwwg=0
      iwwz=0
cvp zeiten for konfigurationsvergleich gesamt
cdebugextimp=0.0d0
cdebugvutimp=0.0d0
cdebugcputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
*     vutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp lesen der integrale for dk=2, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=2 noch nicht fertig
c1012 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 12 k=kanf,kend
      do 12 k=1,nc
cvp
      lz=jz
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3012 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c iww zaehlt alle wechselwirkenden konfigurationen for dk=2
      iww=0
c erzeugung des label-feldes for alle mc konfigurationen
c  und lesen der integrale
c mit check2 werden alle wechselwirkenden konf. der sk(i-2) bestimmt
c sowie die labels
cdebugcall clockv(vu10,cpu10)
cdebugcall time(mill10)
csut  call checkf(ioint,i,j8,k,nspiel,mspiel,jnum1
      call check2(iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,j8,k,nspiel,mspiel,jnum1)
cvp
cvp
cvp
c
c     do 3013 l=1,mc
      do 3013 l=1,nspiel
      iww=iww+1
c     iww=ispiel(l)
      iaddr(iww)=lz
      iaddr(iww)=jz+(ispiel(l)-1)*kmj
      jj=nqrfal(iww)*ii+nis
cvp
c     ip=ot(jj)
      ip=ntabl(jj)
      ipfeld(iww)=ip
c     if (ip.eq.255) ipfeld(iww)=-1
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        cbfeld(iww)=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(iww))) then
c         write(6,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(iww))
c       endif
        jto1=jto1+1
        exfeld(iww)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(intex(iww))) then
c         write(6,*) ' exchange-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intex(iww))
c       endif
       else
        cbfeld(iww)=twoe(intcb(iww))
        exfeld(iww)=-twoe(intex(iww))
      endif
 3013 lz=lz+kmj
cvp
c     write(6,*)' iww=',iww
      if (iww.ne.0) iwwz=iwwz+1
c aufbau der integrale zu allen konf. aus sk(i-2), lll bezeichnet
c  alle z-faelle
      do 5013 lll=1,3
      do 5013 l=1,iww
      bobvec(l,lll)=(ibob2(nrfal(l),lll)*cbfeld(l)+jbob2(nrfal(l),lll)
     &  *exfeld(l))*ipfeld(l)
 5013 continue
c jgcnt zaehlt die matrixelemente
      jgcnt=0
      do 13 l=1,iww
c write dk=2
cvp
      jj=nqrfal(l)
c     jj=1
cvp
c write dk=2
      jj1=jj
      jj=jj*ii+nis
      do 905 m=1,nzer
905   f(m)=0.0d0
c write dk=2
      do 18 m=1,kml
        my=ntabl(jj+m+m-1)
        if (my.ne.0) then
          kx=ja+my
          mike=ntabl(jj+m+m)
          tim=bobvec(l,mike)
          lx=m-kml
          do 22 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
  22        f(lx)=      tim*ae(kx)
c 22      f(lx)=f(lx)+tim*ae(kx)
c
      endif
18    continue
cvp  nqr = # qr-faelle
      nqr=ibinom(inopen,4)
cvp  irec = # eintraege von emat for ein saf-paar
      irec=6*nqr
cvp
      do 23 m=1,kml
        do 23 ih=1,kmj
          tim=0.0d0
          do 24 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
24          tim=tim+f(ly)*ae(mx)
cvp  matrixelemente nach neuer tabelle
c         iem=(ih-1)*kml*irec+(m-1)*irec
c         xtim=
c    &     emat(iem+(nrfal(l)-1)*2*nqr+2*(nqrfal(l)-1)+1)*cbfeld(l)
c    &    -emat(iem+(nrfal(l)-1)*2*nqr+2*(nqrfal(l)-1)+2)*exfeld(l)
c         xtim=xtim*ipfeld(l)
cvp   if (i.eq.3) then
c       if (dabs(xtim-tim).gt.1.0d-8)
c    &   write(6,*) ' l=,m=,if=,tim=,xtim= '
c    &      ,l,m,if,tim,xtim
c     if (nqrfal(l).eq.1)
c    & write(6,*) ' m=,if=,nrfal=,cbfeld=,exfeld=,tim= '
c    &  ,m,if,nrfal(l),cbfeld(l),exfeld(l),tim
cvp   endif
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cdir      if (dabs(tim).ge.1.0e-7) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=kz+m
            mjfeld(jgcnt)=iaddr(l)+ih
cdir      endif
c schreiben der matrixelemente und der indices auf file nhuk
      if (jgcnt.eq.ncmax) then
cvocl loop,novrec,vi(ik)
c
c schleife ueber die wurzeln
c
csut      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut        ipoint=(iwurz-1)*ndim
csut        write(6,*)' 5100 : ipoint=',ipoint
c
         do 4012 ik=1,jgcnt
           qq=qfeld(ik)
           if (dabs(qq).gt.hstor) then
             icount=icount+1
             rham(icount) =(qq)
             imm(icount)=mifeld(ik)
             jmm(icount)=mjfeld(ik)
             if (icount.ge.ndims) then
               irech=irech+1
               icount=0
               write(nf88) imm,jmm,rham
csut           write(6,*) 'write on 88'
               do iii=1,ndims
                 imm(iii) = 0
                 jmm(iii) = 0
                 rham(iii) = 0.0d0
               end do
             end if
           end if
 4012    continue
c ende schleife uber die wurzeln
csut 5100 continue
         jgcnt=0
      endif
cvp
  23  continue
  13  continue
      iwwg=iwwg+iww
c schreiben der matrixelemente und der indices auf file nhuk
cvocl loop,novrec,vi(ik)
c
c schleife ueber die wurzeln
c
csut  do 5200 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut  ipoint=(iwurz-1)*ndim
c
      do 4013 ik=1,jgcnt
        qq   =qfeld(ik)
        if (dabs(qq).gt.hstor) then
           icount=icount+1
           imm(icount)=mifeld(ik)
           jmm(icount)=mjfeld(ik)
           rham(icount) = (qq)
           if (icount.ge.ndims) then
             irech=irech+1
             icount=0
             write(nf88)imm,jmm,rham
csut         write(6,*) 'write on 88'
             do iii=1,ndims
               imm(iii) = 0
               jmm(iii) =0
               rham(iii)=0.0d0
             end do
           end if
         end if
 4013 continue
c ende schleife uber die wurzeln
csut0 continue
cvp
 3012 continue
cvp
  12  kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=2
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1012
cio   endif
cvp
 1014 continue
cvp
c ende dk=2
cdebugwrite(6,*)' jg=',jg
cdebugwrite(6,*)' iwwg=',iwwg
cdebugwrite(6,*)' iwwz=',iwwz
cvp
      go to 10
9     if (i.eq.1) go to 25
10    ibl=ibl+1
c//skipc&2
c     call cpu
c     call setcpu
c start dk=1 fall
c beginn timer-routine for dk=1
*     call clockv(vu1,cpu1)
*     call time(milli1)
c
      j8=i-1
      mc=nconf(j8)
      if (mc.eq.0) go to 25
      js=jan(ibl)
      ks=jbn(ibl)
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      mps=nplu(j8)
      mms=j8-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write dk=1
*      write(iwr,1101)j8,mc,js,ks,ndj,kmj,ja,mr,mps,mms
*1101 format(/1x,'dk=1 (1101):j8,mc,js,ks,ndj,kmj,ja,mr,mps,mms',
*    *i3,2x,i5,2x,8(i6,2x)/)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ix=nad1
      kss=0
 27   if(kss.ne.ks)then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 27
      endif
      call upackx(ot)
      jc=nr+mult
      j3=jc*kml+1
      kz=iz
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,3)*(inopen-2)
     &     +kml*kmj*ibinom(inopen,2)
     &     +2*kml*kmj*ibinom(inopen,2)*(inopen-1)
      if (npdim.le.nedim) then
       if (odebug) then
*        write(6,12347) npdim, nedim
*12347  format(' ***calculation of matrix elements for dk=1 via sga '/
*    +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call sgdk1(emat,nedim,ae,ot,ja,ie,kml,kmj
     &             ,ndj,inopen,ic,i2
     &             ,j3,i3,mps,mms,nms,nps,jc)
      else
       if (odebug) then
        write(6,12348) npdim, nedim
12348  format(' **** sgdk1:  npdim.gt.nedim', 2i7)
       endif
      endif
c
*     extimp=0.0d0
*     vutimp=0.0d0
*     cputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
*     vutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp lesen der integrale for dk=1, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=1 noch nicht fertig
c1028 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 28 k=kanf,kend
      do 28 k=1,nc
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3028 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c mit checkb werden alle wechselwirkenden konf. der sk(i-1) bestimmt
c sowie die labels
*     call clockv(vu10,cpu10)
*     call time(mill10)
csut  call checkb(ioint,i,j8,k,nspiel,mspiel,jnum1
      call check1a(twoe,nteint,iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,j8,k,nspiel,mspiel,jnum1
     &   ,nspp1,nsppa,nsppb)
c
cvp
c     if (nqrfal(l).ne.ndfal(l).and.npfal(l).eq.3) then
c       write(6,*) ' k=,l=,npfal,nqrfal=,ndfal= '
c    &      ,k,l,npfal(l),nqrfal(l),ndfal(l)
c     endif
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,3)*(inopen-2)
     &     +kml*kmj*ibinom(inopen,2)
     &     +2*kml*kmj*ibinom(inopen,2)*(inopen-1)
      if (npdim.le.nedim) then
       if (odebug) then
*        write(6,1234) npdim, nedim
*1234    format(' ***calculation of matrix elements for dk=1 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call hmdk1h(twoe,pey,
     +        emat,nedim,i,nc,mc,kml,kmj,inopen,jz,kz
     &       ,jto1,jg,ioint,iges,kend,nspp1
     &       ,nsppa,nsppb)
         goto 1027
      else
      if (odebug) then
        write(6,1235) npdim, nedim
1235   format(' **** hmdk1h:  npdim.gt.nedim', 2i7)
       endif
      endif
c
c dk=1 p=1
      do 31 l=1,nspp1
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
      ll=nrfal(l)
c ql-fall
      jj=nqlfal(l)
c qr-fall
      kk=nqrfal(l)
      jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      if (ll.lt.0) go to 820
c     if (exc.ne.-twoe(intex(l))) then
c       write(6,*) ' k=,l=,nrfal=,exc=,twoe= ',
c    &   k,l,nrfal(l),exc,twoe(intex(l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul= twoe(jto1)
        jto1=jto1+1
        exc =-twoe(jto1)
       else
        coul=twoe(intcb(l))
        exc=-twoe(intex(l))
      endif
cvp
      go to 38
820   continue
      ll=-ll
cvp   exc=h(jto)
cvp
c     if (exc.ne.twoe(intex(l))) then
c       write(6,*) ' k=,l=,nrfal=,exc=,twoe= ',
c    &   k,l,nrfal(l),exc,twoe(intex(l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul=-twoe(jto1)
        jto1=jto1+1
        exc = twoe(jto1)
       else
        coul=-twoe(intcb(l))
        exc=twoe(intex(l))
      endif
cvp
38    continue
cvp lesen der integrale
c     if (ioint.eq.1) then
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intcb(l))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(l))
c       endif
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intex(l))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' exchange-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intex(l))
c       endif
c     endif
      if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
   49 continue
      do 39 m=1,nzer
39    f(m)=0.0d0
      do 40 m=1,kml
        jj=jj+1
        my=ot(jj)
        jj=jj+1
        if(my.ne.0) then
          kx=ja+my
          mike=ot(jj)
          if(mike.gt.128) then
             tim=-bob(256-mike)
           else
             tim=bob(mike)
          endif
          lx=m-kml
          do 44 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 44        f(lx)=f(lx)+tim*ae(kx)
  44        f(lx)=      tim*ae(kx)
        endif
40    continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trf1h(f,ae,kz,ie,nz,kml,kmj)
cvp
 31   continue
cvp
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
c qr vorzeichen kodiert r fall
      ll=nqrfal(l)
cvp
c     if (sm.ne.twoe(intcb(l))) then
c       write(6,*) ' k=,l=,nrfal=,sm=,twoe= ',
c    &   k,l,nrfal(l),sm,twoe(intcb(l))
c     endif
      if (ioint.eq.1) then
        jto1=jto1+1
        sm=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(l))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(l))
c       endif
       else
        sm=twoe(intcb(l))
      endif
cvp
      if(ll.gt.0) go to 53
      sm=-sm
      ll=-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0.0d0
      do 56 m=1,kml
        jj=jj+1
c welche determinante ww mit der m-ten sdet, vorzeichen kodiert ip
        my=ot(jj)
        if (my.ne.0) then
          if(my.lt.128) then
             tim=sm
             kx=my+ja
           else
             kx=ja-my+256
             tim=-sm
          endif
          lx=m-kml
          do 60 if=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 60        f(lx)=f(lx)+tim*ae(kx)
  60        f(lx)=      tim*ae(kx)
        endif
56    continue
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trf1h(f,ae,kz,ie,nz,kml,kmj)
cvp
  32  continue
cvp
c dk=1,p=3,r=1,2
      do 33 l=nsppa+1,nsppb
      ll=nqrfal(l)
      ntar=nrfal(l)
      nb=moafal(l)
      mb=mobfal(l)
      nb=ideks(nb)+mb+ijone(ntar)
cvp   nb=ideks(nb)+mb+ij(ntar)
c     sm=pey(nb)
      sm1=pey(nb)
c lesen der summe ueber die dreizentrenintegrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(6,*) ' p=3, jto1= ',jto1
c         write(6,*) ' integral sum: k=,l= ',k,l,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(l)
      endif
      sm=tim1
cvp
c write dk=1,p=3
      if (mr.eq.0 ) go to 65
c lesen der (negativen) austauschintegrale
      do 66 m=1,mr
cvp
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(m)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(l,m))) then
c         write(6,*) ' nsppb=,jto1= ',nsppb,jto1
c         write(6,*) ' integrals: k=,l=,m= ',k,l,m,twoei(jto1)
c    &    ,twoe(nwwmo(l,m))
c       endif
       else
        sac(m)=-twoe(nwwmo(l,m))
      endif
cvp
c     if (sac(m).ne.-twoe(nwwmo(l,m))) then
c       write(6,*) ' k=,l=,m=,sac=,-nwwmo= ',k,l,m,sac(m)
c    &   ,-twoe(nwwmo(l,m))
c     endif
cvp
66    continue
   65 continue
c     iiiiii=iiiiii+1
c     if (iiiiii.lt.100)
c    &write(6,*)' jto=,sm=,sac=,sm1=',jto,sm,sac,sm1
      sm=sm+sm1
c abzweigung r=1 bzw 2
cvp qr < 0 --> r=1; qr > 0 --> r=2
      if(ll.lt.0) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
cvp mr = # der offenen schalen der kleineren sk
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
      do 80 m=1,nzer
80    f(m)=0
      jx=-kml
      do 73 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (74,75,76,77),im
 74      in=in+1
         kx=ot(in)
         kx=kx+ja
         tim=sm
cvp nps = # alpha-spins der offenen schalen der kleineren sk
         if(mps.eq.0) go to 78
         do 79 ih=1,mps
           in=in+1
           im=ot(in)
79         tim=tim+sac(im)
78       lx=jx
         do 81 ih=1,kmj
           kx=kx+ndj
           lx=lx+kml
81         f(lx)=f(lx)+tim*ae(kx)
        go to 73
 75     in=in+1
        kx=ot(in)
        kx=kx+ja
        tim=-sm
        if(mms.eq.0) go to 78
cvp mms = nr. der kleineren sk - 1 ; bzw.
cvp mms = #  beta-spins der offenen schalen der kleineren sk
        do 83 ih=1,mms
          in=in+1
          im=ot(in)
83        tim=tim-sac(im)
        go to 78
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
76      do 84 ih=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=-sac(im)
          lx=jx
          do 84 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
84          f(lx)=f(lx)+tim*ae(kx)
        go to 73
cvp nps = # alpha-spins der offenen schalen
77      do 85 ih=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 85 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
85          f(lx)=f(lx)+tim*ae(kx)
   73 jj=jj+jc
      go to 47
c p=3 anderer r fall (r=1)
  86  jj=-1*ll*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (91,95,96,97),im
91        in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=sm
          if(mms.eq.0) go to 92
          in=in+mps
          do 93 ih=1,mms
            in=in+1
            im=ot(in)
93          tim=tim+sac(im)
92        lx=jx
          do 94 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
94          f(lx)=f(lx)+tim*ae(kx)
          go to 90
95        in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=-sm
          if(mps.eq.0) go to 92
          in=in+mms
          do 99 ih=1,mps
            in=in+1
            im=ot(in)
99          tim=tim-sac(im)
          go to 92
96        do 100 ih=1,nms
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 100 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
100         f(lx)=f(lx)+tim*ae(kx)
          go to 90
97        do 101 ih=1,nps
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
            im=ot(in)
            tim=-sac(im)
            lx=jx
            do 101 ig=1,kmj
               kx=kx+ndj
               lx=lx+kml
101            f(lx)=f(lx)+tim*ae(kx)
90         jj=jj+jc
cre
   47 continue
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trf1h(f,ae,kz,ie,nz,kml,kmj)
cvp
  33  continue
cvp
 1027 continue
cvp
 3028 continue
cvp
28    kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=1
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1028
cio   endif
cvp
*     extimp=extimp/1000
*     write(6,*) 'timings for sum of all check1-calls '
*     write(6,*) 'vu   time ',vutimp,'  seconds'
*     write(6,*) 'cpu  time ',cputip,'  seconds'
*     write(6,*) 'ex   time ',extimp,'  seconds'
c ende der timer-routine
*     call clockv(vu2,cpu2)
*     call time(milli2)
*     extime = milli2 - milli1
*     extime = extime/1000
*     vutime = vu2  - vu1
*     cputim = cpu2 - cpu1
*     write(6,*) 'timings for case dk=1'
*     write(6,*) 'vu   time ',vutime,'  seconds'
*     write(6,*) 'cpu  time ',cputim,'  seconds'
*     write(6,*) 'ex   time ',extime,'  seconds'
cvp
*     write(6,*) '  vu - timings for a, b, c ',vutwwg,vutlag,vuting
*     write(6,*) ' cpu - timings for a, b, c ',cpuwwg,cpulag,cpuing
*     write(6,*)
cvp
*     if (npdim.le.nedim) then
*        write(6,*) 'timings for hmdk1 '
*        write(6,*) 'vu   time ',vux1,'  seconds'
*        write(6,*) 'cpu  time ',cpux1,'  seconds'
*        write(6,*) 'ex   time ',xillx1/1000,'  seconds'
*     endif
c
25    js=idra(i)
c
c           *******   start dk=0   *******
c
c start dk=0 p=5 fall
c beginn timer-routine for dk=0 alle p-faelle
c beginn timer-routine for dk=0,p=5
cdebugcall clockv(vu1,cpu1)
cdebugcall time(milli1)
c
      ks=idrc(i)
      ix=nad1
c write dk=0 p=5
*      write(iwr,1500)js,ks
*1500 format(1x,'ks,js for dk=0 p=5 (1500)',2i5)
      kss=0
 102  if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 102
      endif
*****
      call upackx(ot)
*      write(6,*) 'after dread'
      kz=iz
cvp
      ig=0
cvp
      jstar=niot(i)-nc
      do 103 k=1,nc
cvp nx = # besetzte mo's
      do 104 l=1,nx
         kstar=jstar+nc*l
         itest(l)=iottr(kstar+k)
 104  continue
c write dk=0 p=5
      do 105 l=1,mzer
105   f(l)=0.0d0
      care=core
cvp nr = # einfach besetzte mo's
      do 107 l=1,nr
         mal=itest(l)
         ntar=nirred(mal)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
107      care=care+pey(mal)
cvp nd = # doppelt besetzte mo's
cvp n1 = nr+1
      do 109 l=n1,nx
         mal=itest(l)
         ntar=nirred(mal)
         nal=ideks(mal+1)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
         sm=pey(mal)
109      care=care+sm+sm+acoul(nal)
      do 111 l=2,nr
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 111 m=1,it
           nal=mal+itest(m)
111        care=care+acoul(nal)
      sm=0.0d0
cvp n2 = nr + 2
      do 113 l=n2,nx
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 113 m=n1,it
           nal=mal+itest(m)
           tim=acoul(nal)
113        sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
      do 115 l=1,nr
         mal=itest(l)
         do 115 m=n1,nx
           nal=itest(m)
           nal=ideks(max(nal,mal))+min(nal,mal)
           sm=acoul(nal)
115        care=care+sm+sm-aexc(nal)
      kk=ic
      kn=i2
      ml=i3
c write dk=0 p=5
      jx=-kml
      do 118 l=1,kml
        jx=jx+1
        lx=jx
        kx=ohog(l)+ia
        tim=care
cvp nps = # alpha-spins der offenen schalen
c       if(nps.lt.2) go to 122
        do 123 ih=2,nps
          mx=i2+nps*(l-1)+ih
          it=ih-1
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          do 123 in=1,it
            mz=i2+nps*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
123         tim=tim-aexc(nal)
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
c       if(nms.lt.2) go to 122
        do 124 ih=2,nms
          mx=i3+nms*(l-1)+ih
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          it=ih-1
          do 124 in=1,it
            mz=i3+nms*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
124         tim=tim-aexc(nal)
        do 121 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
121       f(lx)=f(lx)+tim*ae(kx)
cvp niw = nps * nms
        do 120 m=1,niw
          kk=kk+1
          kx=out(kk)+ia
c write dk=0 p=5 ausserdiagonal
          kk=kk+1
          mal=out(kk)
          kk=kk+1
          nal=out(kk)
          mal=itest(mal)
          nal=ideks(mal)+itest(nal)
          tim=-aexc(nal)
          lx=jx
          do 120 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
120         f(lx)=f(lx)+tim*ae(kx)
118   continue
      do 125 l=1,kml
        mz=kz+l
        do 125 m=1,l
          tim=0.0d0
          do 126 ih =1,kml
            ly=(m-1)*kml+ih
            mx=ie+(l-1)*kml+ih
126         tim=tim+f(ly)*ae(mx)
          if(l.gt.m) go to 127
csut          jg=jg+1
cx        q(jg)=tim
cvp diagonale wird jetzt aus hmp5 genommen
c         q(jg)=hp5(mz)
csut          tim=hp5(mz)
csut          iih=mz
csut          jjh=iih
c --- diagonale !!!
c
c schleife ueber die wurzeln
c
csut      do 5300 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut        ipoint=(iwurz-1)*ndim
c
csut          vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
c ende schleife uber die wurzeln
csut 5300 continue
          go to 125
 127      continue
          jjh   =kz+m
          if (dabs(tim).gt.hstor) then
            icount=icount+1
            imm(icount) = mz
            jmm(icount) = jjh
            rham(icount)= (tim)
            if (icount.ge.ndims) then
              irech=irech+1
              icount=0
              write(nf88)imm,jmm,rham
csut          write(6,*) 'write on 88'
              do iii=1,ndims
                imm(iii) = 0
                jmm(iii) = 0
                rham(iii) =0.0d0
              end do
            end if
          end if
csut      do 5400 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut        ipoint=(iwurz-1)*ndim
c
csut          vecf1(ipoint+mz)=vecf1(ipoint+mz)+vecf2(ipoint+jjh)*tim
csut          vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+mz)*tim
c ende schleife uber die wurzeln
csut 5400 continue
125   continue
103   kz=mz
c ende der timer-routine
*     call clockv(vu2,cpu2)
*     call time(milli2)
*     extime = milli2 - milli1
*     extime = extime/1000
*     vutime = vu2  - vu1
*     cputim = cpu2 - cpu1
*     write(6,*) 'timings for case dk=0,p=5'
*     write(6,*) 'vu   time ',vutime,'  seconds'
*     write(6,*) 'cpu  time ',cputim,'  seconds'
*     write(6,*) 'ex   time ',extime,'  seconds'
c
      if(nc.eq.1) go to 6
c start dk=0 andere p-faelle
      ks=jdrc(i)
      ix=nad1
c write dk=0 andere p-faelle
c     write(iwr,1520)js,ks
*1520 format(1x,'dk=0 other p cases js,ks',2(i7,2x))
      kss=0
128   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 128
      endif
*****
      call upackx(ot)
      ii=ii+kml
      j2=kml+1
      jc=nr+nr
      j3=jc*kml+1
      kz=iz+kml
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kml*ideks(ibinom(inopen,2)+1)
     &     +kml*kml*ideks(inopen+1)
     &     +kml*kml*ideks(inopen+1)*2*inopen
      if (npdim.le.nedim) then
       if (odebug) then
*       write(6,1236) npdim, nedim
*1236   format(' ***calculation of matrix elements for dk=0 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call sgdk0(emat,nedim,ae,ot,ie,kml,
     &              ndt,inopen,ic,i2,ii,ia,j2,
     &              j3,i3,nps,np1,nms,nm1,jc,np2,nm2)
*        vux1=0.0d0
*        cpux1=0.0d0
*        xillx1=0.0d0
      else
       if (odebug) then
        write(6,1237) npdim, nedim
1237   format(' **** sgdk0:  npdim.gt.nedim', 2i7)
       endif
      endif
c die untere dreickmatrix wird berechnet
      j8=0
*     extimp=0.0d0
*     vutimp=0.0d0
cdebugcputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
cdebugvutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp
cvp lesen der integrale for dk=0, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=2
cio   kend=nc
cvp ruecksprungadresse, falls dk=0 noch nicht fertig
c1129 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cio   do 129 j=kanf,kend
      do 129 j=2,nc
cvp   lz=iz
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(j-2)/ncmax+1
cvp
cvp
cvp zerlegung der schleife ueber alle j-1 konfigurationen
cvp  der unteren dreiecksmatrix in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3129 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=j-1-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c mit checkd werden alle wechselwirkenden konf. der sk(i) bestimmt
c sowie die labels
cdebugcall clockv(vu10,cpu10)
cdebugcall time(mill10)
csut  call checkda(twoe,nteint,ioint,i,i,j,nspiel,mspiel,jnum1
      call check0a(twoe,nteint,iottr,iotnew,iotm,
     + maindf,maxr,jconb,nomax2,
     + ioint,i,i,j,nspiel,mspiel,jnum1
     & ,nspp1,nsppa,nsppb,nspp4)
c
cvp   j8=j8+1
cvp
c     if (moafal(k).ne.ndfal(k).and.npfal(k).ge.4) then
c       write(6,*) ' j=,k=,npfal=,moafal=,ndfal= ',
c    &   j,k,npfal(k),moafal(k),ndfal(k)
c       if (nqrfal(k).lt.0) then
c         write(6,*) ' r= 2'
c        else
c         write(6,*) ' r= 1'
c       endif
c     endif
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kml*ideks(ibinom(inopen,2)+1)
     &     +kml*kml*ideks(inopen+1)
     &     +kml*kml*ideks(inopen+1)*2*inopen
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1238) npdim, nedim
*1238   format(' ***calculation of matrix elements for dk=0 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call hmdk0h(twoe,pey,aexc,
     +               emat,nedim,i,nc,mc,kml,kmj,inopen,iz,kz
     &             ,jto1,jg,ioint,iges,kend,nspp1
     &             ,nsppa,nsppb,nspp4)
         goto 1129
      else
      if (odebug) then
        write(6,1239) npdim, nedim
1239   format(' **** hmdk0h:  npdim.gt.nedim', 2i7)
       endif
      endif
c dk=0 p=1
      do 132 k=1,nspp1
c r-fall
      ll=nrfal(k)
cvp qr-fall
      jj=nqrfal(k)
cvp ql-fall
      kk=nqlfal(k)
      if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob=0
      go to 151
150   jj=ideks(kk)+jj
      ibob=1
151   jj=jj*ii+ic
c die allgemeine permutation
      ip=ot(jj)
c write dk=0,p=1
      if(ip.eq.255) go to 830
cvp   coul=h(jto)
cvp
c     if (coul.ne.twoe(intcb(k))) then
c       write(6,*) ' j=,k=,nrfal=,coul=,twoe= ',
c    &   j,k,nrfal(k),coul,twoe(intcb(k))
c     endif
cvp
c     if (exc.ne.-twoe(intex(k))) then
c       write(6,*) ' j=,k=,nrfal=,exc=,twoe= ',
c    &   j,k,nrfal(k),exc,twoe(intex(k))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul= twoe(jto1)
        jto1=jto1+1
        exc =-twoe(jto1)
       else
        coul=twoe(intcb(k))
        exc=-twoe(intex(k))
      endif
cvp
      go to 153
 830  continue
c830  coul=-h(jto)
cvp   if (coul.ne.-twoe(intcb(k))) then
cvp     write(6,*) ' j=,k=,nrfal=,coul=,twoe= ',
cvp  &   j,k,nrfal(k),coul,twoe(intcb(k))
cvp   endif
cvp
c     if (exc.ne.twoe(intex(k))) then
c       write(6,*) ' j=,k=,nrfal=,exc=,twoe= ',
c    &   j,k,nrfal(k),exc,twoe(intex(k))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul=-twoe(jto1)
        jto1=jto1+1
        exc = twoe(jto1)
       else
        coul=-twoe(intcb(k))
        exc=twoe(intex(k))
      endif
cvp
  153 continue
cvp lesen der integrale
c     if (ioint.eq.1) then
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intcb(k))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' cb-integrale: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intcb(k))
c       endif
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intex(k))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' ex-integrale: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intex(k))
c       endif
c     endif
      if (ll-2) 832,833,834
  832 bob(1)=coul+exc
      bob(2)=exc
      bob(3)=coul
      go to 835
  833 bob(1)=exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 835
 834  bob(1)=-exc
      bob(2)=-coul-exc
      bob(3)=coul
835   do 154 l=1,mzer
154   f(l)=0
      mx=-kml
      do 155 l=1,kml
        jj=jj+1
        mx=mx+1
        ip=ot(jj)
        if(ip.eq.2)  go to 159
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        jj=jj+1
        tim=bob(1)
        lx=mx
        do 161 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c161      f(lx)=f(lx)+tim*ae(kx)
 161      f(lx)=      tim*ae(kx)
        go to 155
159     jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=bob(2)
        lx=mx
        do 160 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c160      f(lx)=f(lx)+tim*ae(kx)
 160      f(lx)=      tim*ae(kx)
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=bob(3)
        lx=mx
        do 167 m=1,kml
          lx=lx+kml
          kx=kx+ndt
167       f(lx)=f(lx)+tim*ae(kx)
155   continue
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdkh(f,ae,kz,ie,nz,kml,ibob)
cvp
 132  continue
cvp
c dk=0,p=2
      do 133 k=nspp1+1,nsppa
cvp qr-fall
      ll=nqrfal(k)
cvp ql-fall
      jj=nqlfal(k)
      if(jj.gt.ll) go to 170
      jj=ideks(ll)+jj
      ibob=0
      go to 171
170   jj=ideks(jj)+ll
      ibob=1
  171 continue
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(k))) then
c         write(6,*) ' jto1= ',jto1
c         write(6,*) ' cb-integrale: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intcb(k))
c       endif
       else
        tim=twoe(intcb(k))
      endif
cvp
c write dk=0 p=2
      jj=jj*j2+i2
      do 177 l=1,mzer
177   f(l)=0
      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) tim=-tim
      do 173 l=1,kml
        jx=jx+1
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        lx=jx
        do 173 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c 173     f(lx)=f(lx)+tim*ae(kx)
  173     f(lx)=      tim*ae(kx)
cvp
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdkh(f,ae,kz,ie,nz,kml,ibob)
cvp
 133  continue
cvp
c dk=0 p=3
      do 134 k=nsppa+1,nsppb
c qr, vorzeichen kodiert den r-fall
cvp qr > 0 --> r=1
      ll=nqrfal(k)
c ql
      jj=nqlfal(k)
c symmetrie der orbitale for einelektronen-integral
      ntar=nrfal(k)
c 1.mo for einelektronen-integral (das groessere)
      nb=moafal(k)
c 2.mo for einelektronen-integral (das kleinere)
      mb=mobfal(k)
      nb=ideks(nb)+mb+ijone(ntar)
cvp einelektronenintegral wird separat behandelt
c     sm=pey(nb)
      sm1=pey(nb)
      sm=0.0d0
c abzweigung zu den unterschiedlichen r-faellen
      if(ll.lt.0) go to 181
c dk=0 p=3 r=1
      if(ll.lt.jj) go to 182
      ibob=0
      jj=ideks(ll)+jj
c write dk=0 p=3 r=1
      go to 205
182   jj=ideks(jj)+ll
c write dk=0 p=3 r=1
      ibob=1
205   jj=jj*j3+i3
c lesen der summe ueber die dreizentrenintegrale
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(6,*) ' r=1,jto1= ',jto1
c         write(6,*) ' integral sum: j=,k= ',j,k,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(k)
      endif
      sm=tim1
cvp
      if(nr.eq.1) go to 185
c lesen der negativen austauschintegrale
cvp n3 = nr - 1
      do 186 l=1,n3
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(l)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(k,l))) then
c         write(6,*) ' r=1,jto1= ',jto1
c         write(6,*) ' single integrals: j=,k=,l= ',j,k,l,twoei(jto1)
c    &    ,twoe(nwwmo(k,l))
c       endif
       else
        sac(l)=-twoe(nwwmo(k,l))
      endif
cvp
c write dk=0 p=3 r=1
186   continue
  185 continue
c einelektronenintegral wird addiert
      sm=sm+sm1
      ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      if (nr.eq.1) go to 191
      do 192 l=1,n3
192   sac(l)=-sac(l)
191   jj=jj+1
      do 193 l=1,mzer
193   f(l)=0
      jx=-kml
      do 194 l=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=sm
        if(im.eq.2) go to 195
        if(nps.eq.1) go to 196
        do 197 m=1,np1
          in=in+2
          nn=ot(in)
197       tim=tim+sac(nn)
196     do 198 m=1,kml
          lx=lx+kml
          kx=kx+ndt
c198      f(lx)=f(lx)+tim*ae(kx)
 198      f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 194
        do 199 m=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 199 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
199         f(lx)=f(lx)+tim*ae(kx)
        go to 194
195     if(nms.eq.1) go to 200
        do 201 m=1,nm1
          in=in+2
          nn=ot(in)
201       tim=tim+sac(nn)
200     do 202 m=1,kml
          lx=lx+kml
          kx=kx+ndt
202       f(lx)=f(lx)+tim*ae(kx)
        do 203 m=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 203 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
203         f(lx)=f(lx)+tim*ae(kx)
194     jj=jj+jc
cvp   if (ibob.eq.1) go to 144
      go to 146
c dk=0 p=3 r=2
  181 ll=-1*ll
      if(ll.lt.jj) go to 204
      ibob=0
      jj=ideks(ll)+jj
c write dk=0 p=3 r=2
      jj1=jj
      go to 206
204   jj=ideks(jj)+ll
      ibob=1
206   jj=jj*j3+i3
c erstes integral for r=2 siehe aufzeichnung
c lesen der summe der dreizentrenintegrale
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(6,*) ' r=2,jto1= ',jto1
c         write(6,*) ' integral sum: j=,k= ',j,k,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(k)
      endif
cvp
      sm=sm+tim1
      if(nr.eq.1) go to 209
c do 211 ist analog do 186
      do 211 l=1,n3
cvp
c     sac(l)=h(jto)
c     if (sac(l).ne.-twoe(nwwmo(k,l))) then
c       write(6,*) ' j=,k=,l=,sac=,-nwwmo= ',j,k,l,sac(l)
c    &   ,-twoe(nwwmo(k,l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(l)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(k,l))) then
c         write(6,*) ' r=2,jto1= ',jto1
c         write(6,*) ' single integrals: j=,k=,l= ',j,k,l,twoei(jto1)
c    &    ,twoe(nwwmo(k,l))
c       endif
       else
        sac(l)=-twoe(nwwmo(k,l))
      endif
cvp
c write dk=0 p=3 r=2
211   continue
  209 continue
c einelektronenintegrale werden addiert
      sm=sm+sm1
      ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      if(nr.eq.1) go to 215
      do 216 l=1,n3
216   sac(l)=-sac(l)
215   jj=jj+1
      do 217 l=1,mzer
217   f(l)=0
      jx=-kml
      do 218 l=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=-sm
        if (im.eq.2) go to 219
        if(nms.eq.0) go to 220
        in=in+np2
        jn=in
        do 221 m=1,nms
          jn=jn+2
          nn=ot(jn)
221       tim=tim-sac(nn)
220     do 222 m=1,kml
          lx=lx+kml
          kx=kx+ndt
c 222     f(lx)=f(lx)+tim*ae(kx)
  222     f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 218
        do 223 m=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 223 ih =1,kml
            lx=lx+kml
            kx=kx+ndt
223         f(lx)=f(lx)+tim*ae(kx)
        go to 218
219     in=in+nm2
        jn=in
        do 224 m=1,nps
          jn=jn+2
          nn=ot(jn)
224       tim=tim-sac(nn)
        do 225 m=1,kml
          lx=lx+kml
          kx=kx+ndt
225       f(lx)=f(lx)+tim*ae(kx)
        do 226 m=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 226 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
226         f(lx)=f(lx)+tim*ae(kx)
218     jj=jj+jc
cvp
 146  continue
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdkh(f,ae,kz,ie,nz,kml,ibob)
cvp
 134  continue
cvp
c dk=0 p=4     aa/bb
      do 135 k=nsppb+1,nsppb+nspp4
c missmatch in l
cvp absolutnummer von mo a
      kk=mobfal(k)
c missmatch in r
cvp absolutnummer von mo b
      ll=moafal(k)
      if(kk.lt.ll) go to 138
      kk=ideks(kk)+ll
      go to 139
138   kk=ideks(ll)+kk
139   continue
      tim=aexc(kk)
      do 141 l=1,kml
cvp
          iih=kz+l
          jjh=iz+(ispiel(k)-1)*kml+l
          if (dabs(tim).gt.hstor) then
             icount=icount+1
             imm(icount)=iih
             jmm(icount)=jjh
             rham(icount)= tim
             if (icount.ge.ndims) then
               irech=irech+1
               icount=0
               write(nf88)imm,jmm,rham
csut           write(6,*) 'write on 88'
               do iii=1,ndims
                 imm(iii) = 0
                 jmm(iii) = 0
                 rham(iii) = 0.0d0
               end do
             end if
           end if
c
c schleife ueber die wurzeln
c
csut  do 5500 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut    ipoint=(iwurz-1)*ndim
c
csut      vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
csut      vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*tim
c ende schleife uber die wurzeln
csut 5500 continue
  141 continue
 135  continue
cvp
 1129 continue
cvp
 3129 continue
cvp
129   kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=0
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1129
cio   endif
cvp
c
6     continue
cold  jg=jg+1
cold  mi(jg)=-1
cold  mj(jg)=0
cvp   write (6,229) mj
cold  write(nhuk)q,mi,mj
c ende der timer-routine for skiny
c --- behalten
      if (ioint.eq.1) then
cdebug   write(6,*) 'rewind on file ',mtype
         call rewftn(mtype)
      end if
      return
      end
cvp
cc
c berechnung der matrixelemente h(i,j)
c rwah skiny eruit!!
      subroutine skiny(twoe,nteint,
     +                 pey,acoul,aexc,ndeke,core,
     +                 vecf1,vecf2,mdi,
     +                 emat,nedim,
     +                 iottr,iotnew,iotm,
     +                 maindf,maxr,jconb,nomax2,ioint,
     +                 odebug)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
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
      real*8 ciprnt
      logical oprint,oecho
      common/prints/oprint(60),ciprnt,oecho
c
      integer nteint, ndeke
      logical odebug
      real*8 twoe
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
cvp
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
c --- issk : superkategorien auf disk
      real*8 hstor
      integer issk
      common /csemi/ hstor,issk
c common for multi-root version
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common /rwork2/ ndim,nanz
cvp
cvp table information: insbes. c- und b-matrix
      real*8 af, ae
      common /rtable/ af(8000)
      dimension ae(7829)
      dimension nsac(iswhm),idrc(iswhm),
cbe jan(5),jbn(5) in jan(7),jbn(7) geaendert (siehe equivalence)
     4  idra(iswhm),jdrc(iswhm),jan(7),jbn(7),
     4  ndet(iswhm),iaz(iswhm),iez(iswhm),
     6  j9(49)
      equivalence (af(50),ae(1))
      integer nsac,ndet,iaz,iez,jan,jbn,idra,idrc
      integer jdrc,j9
      equivalence (af(1),j9(1)),(j9(1),ndet(1)),
     +            (j9(6),nsac(1)),(j9(11),iaz(1)),
     +            (j9(16),iez(1)),(j9(21),idra(1)),
     +            (j9(26),idrc(1)),(j9(31),jdrc(1)),
     +            (j9(36),jan(1)),(j9(43),jbn(1))
cvp
      common/linkmr/ic,i2,i3
      common/blksi3/nsz
      common/bufd/gout(510),nword
      dimension lout(510)
      equivalence (lout(1),gout(1))
      integer ot
      common /rwork/ ot(motmax)
      common /ra/ bob(3),f(ksafm*ksafm)
      integer out, ohog
      dimension ohog(48),out(10962)
      equivalence (ot(13),ohog(1)),(ot(13),out(1))
c
      common /rb/ ifrk, iswh
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
cvp
      integer isc
      dimension isc(iswhm)
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
      integer ijone,iwod,mconf,nplu,idummy
      common /rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
cvp
cvp entwicklungskoeffizienten von saf's in slater-determinanten
cvp   dimension cmat(ncdim),emat(nedim),istvec(iswhm)
cvp emat : darstellungsmatrizen
cvp
c bobvec=ww-terme for alle konf. aus sk(i-2) und alle z-faelle,
      real*8 bobvec
      common/scra/bobvec(ncmax,3)
c
      real*8 tim1
cvp
      integer inull
      data inull/0/
c
cdebug
cdebugwrite(iwr,*) ' skiny: ioint= ',ioint
cdebugwrite(iwr,*) ' skiny: ndim = ',ndim
cdebugwrite(iwr,*) ' skiny: nanz = ',nanz
cdebugwrite(iwr,*) ' skiny: vecf1  '
cdebugwrite(iwr,*) (vecf1(iii),iii=1,5)
cdebugwrite(iwr,*) ' skiny: vecf2  '
cdebugwrite(iwr,*) (vecf2(iii),iii=1,5)
cdebug
        call flushn(iwr)
      ifr1=1-ifrk
csut  nad=4000
csut  nad1=-3999
c
      is=nplu(1)
      mult=is+1
cvp
cvp erzeugung der matrix for entwicklung einer saf in slater-det.
cvp  mult = 2*s +1
cvp   call saf1(mult,cmat,istvec)
cvp
cvp  test: umspeichern von buenkers c-matrix auf cmat
c     do 4713 i=1,iswh
c        nndet=ndet(i)
c        kmml=nsac(i)
c        istb=iaz(i)
c        do 4712 iiii1=1,kmml
c        do 4712 iiii2=1,nndet
c4712    cmat(istvec(i)+(iiii1-1)*nndet+iiii2)=
c    &       ae(istb+iiii1*nndet+iiii2)
c4713 continue
cvp  berechnung der dimension der matrix bis zur sk i-1
      jsum=0
      do 901 i=1,iswh
csut  do 901 i=(issk+1),iswh
      isc(i)=jsum
  901 jsum=jsum+nconf(i)*nsac(i)
cvp
cvp
cvp
      ibl=0
      jg=0
cvp
cvp lesen der integrale, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
      if (ioint.eq.1) then
csut     rewind (mtype)
csut     read(mtype) iges
csut     write(iwr,*) 'iges=',iges
         call rewftn(mtype)
         read (mtype) iges,kend,(twoe(iii),iii=1,iges)
         jto1=0
      endif
cvp
c statement for rumple test
c zaehler for konfigurationspaare
      call setsto(motmax,0,ot)
      do 6 i=1,iswh
c write statement eingefuegt
*      write(iwr,500) i
*  500 format(/50x,'*******   begin of sk  ',i2,'  *******'/)
      nc=nconf(i)
      if (i.le.issk) then
        ncx55 = 0
      else
        ncx55 = nc
      end if
csut  if (nc.gt.0) go to 7
      if (ncx55.gt.0) go to 7
      md=mconf(i)
      if (i.lt.3) go to 8
      ibl=ibl+2
*      write(iwr,*)'not go to 8 ibl=',ibl
      go to 6
8     ibl=ibl+i-1
*      write(iwr,*)'go to 8 ibl=',ibl
      go to 6
7     ndt=ndet(i)
c//skipc&2
c     call cpu
c     call setcpu
      ia=iaz(i)
      kml=nsac(i)
      ndl=ndet(i)
      ie=iez(i)
      iz=isc(i)
      mzer=kml*kml
      nd=ndub(i)
      nps=nplu(i)
      nms=i-1
      niw=nps*nms
      nr=nod(i)
      nx=nytl(i)
      n1=nr+1
      n2=nr+2
      n3=nr-1
      np1=nps-1
      np2=np1+np1
      nd1=nd-1
      nm1=i-2
      nm2=nm1+nm1
      ii=kml+kml+1
      nis=1-ii
      inopen=nod(i)
c write allgemein
c     write(iwr,2000) i,nc,ndl,kml,nd,nr,nps
*2000 format(/1x,'(2000)*****  sk=',i2,' # konf =',i5,' # det =',
*    *i3,'# sdet = ',i2,'# cl.sh. =',i2,' # op.sh. =',i2,
*    *' # alpha s = ',i2/)
      if (i.lt.3) go to 9
      ibl=ibl+1
*      write(iwr,*)' before start dk=2: ibl=',ibl
csut  if (i.le.issk) goto 6
c start dk=2 fall
c write statement eingefuegt
*      write(iwr,501)
* 501 format(/20x,'*******   start dk=2   *******'/)
*      write(iwr,*) 'ibl = ',ibl
      j8=i-2
      mc=nconf(j8)
csut  if (i.le.issk) then
csut    mcx55=0
csut  else
        mcx55 = mc
csut  end if
      if(mcx55.eq.0) go to 10
*      write(iwr,*) 'ibl (dk=2)',ibl,'sk = ',i
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      js=jan(ibl)
      ks=jbn(ibl)
cvp
c     if (i.eq.3) then
c     write(iwr,*) ' c-matrix of buenker '
c     ix1=ndet(i)
c     ix2=iaz(i)
c     do 4711 lll=1,kml
c       ix2=ix2+ix1
c4711   write(iwr,*) (ae(ix2+mmm),mmm=1,ix1)
c     endif
cvp
c write allgemein
c     write(iwr,2001) j8,mc,ndj,kmj,nr,js,ks
c2001 format(1x,'dk=2 (2001) j8,mc,ndj,kmj,nr,js,ks',7(i5,1x))
      ix=nad1
      kss=0
 11   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 11
      endif
*****
      call upackx(ot)
cvp
csut -- ruassprung
csut   if (i.le.issk) goto 701
cvp konvertierung von ot auf das integer feld ntabl
      call tabcon(inopen,kml,ot)
cvp
cvp
cvp erzeugung der tabelle for die berechnung der matrixelemente
c     call matdk2(cmat,emat,istvec,i,mult,kml,kmj,inopen)
cvp
      kz=iz
c beginn timer-routine for dk=2
*     call date(datum)
*     call time(milli)
*     isek=milli/1000
*     imin=isek/60
*     istu=imin/60
*     imin=imin-istu*60
*     write(iwr,*) 'calculation were performed : ',datum,' time ',
*    *istu,':',imin
cdebugcall clockv(vu1,cpu1)
cdebugcall time(milli1)
c aufbau der table-faktoren for r- und z-faelle:
c  ibob2 -> coulomb-integrale, jbob2 -> exchange-integrale
c  (spaeter evtl. als data statement)
      ibob2(1,1)=-1
      ibob2(1,2)=0
      ibob2(1,3)=1
      ibob2(2,1)=0
      ibob2(2,2)=1
      ibob2(2,3)=-1
      ibob2(3,1)=0
      ibob2(3,2)=1
      ibob2(3,3)=-1
      jbob2(1,1)=-1
      jbob2(1,2)=1
      jbob2(1,3)=0
      jbob2(2,1)=-1
      jbob2(2,2)=1
      jbob2(2,3)=0
      jbob2(3,1)=1
      jbob2(3,2)=0
      jbob2(3,3)=-1
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,4)
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1240) npdim, nedim
*1240    format(' ***calculation of matrix elements for dk=2 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call sgdk2(emat,nedim,ae,ja,ie,kml,kmj
     &             ,ndj,inopen)
         call hmdk2(twoe,
     +              vecf1,vecf2,mdi,
     +              emat,nedim,iottr,iotnew,iotm,
     +              maindf,maxr,jconb,nomax2,
     +              i,nc,mc,kml,kmj,
     +              inopen,jz,kz,
     &              mtype,jto1,ioint,iges,kend)
         goto 1014
      else
       if (odebug) then
         write(6,1241) npdim, nedim
1241    format(' **** sgdk2:  npdim.gt.nedim', 2i7)
       endif
      endif
c schleife ueber alle konfigurationen
cdebugwrite(iwr,*)' nc=',nc
cdebugwrite(iwr,*)' mc=',mc
cdebugwrite(iwr,*)' kml=',kml
cdebugwrite(iwr,*)' kmj=',kmj
c iwwg zaehlt alle ww konfigurationen for dk=2
c  iwwz=zaehler aller von null verschiedenen iww's
      iwwg=0
      iwwz=0
cvp zeiten for konfigurationsvergleich gesamt
cdebugextimp=0.0d0
cdebugvutimp=0.0d0
cdebugcputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
*     vutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp lesen der integrale for dk=2, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=2 noch nicht fertig
c1012 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 12 k=kanf,kend
      do 12 k=1,nc
cvp
      lz=jz
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3012 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c iww zaehlt alle wechselwirkenden konfigurationen for dk=2
      iww=0
c erzeugung des label-feldes for alle mc konfigurationen
c  und lesen der integrale
c mit check2 werden alle wechselwirkenden konf. der sk(i-2) bestimmt
c sowie die labels
cdebugcall clockv(vu10,cpu10)
cdebugcall time(mill10)
csut  call checkf(ioint,i,j8,k,nspiel,mspiel,jnum1
      call check2(iottr,iotnew, iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,j8,k,nspiel,mspiel,jnum1)
cvp
c     do 3013 l=1,mc
      do 3013 l=1,nspiel
      iww=iww+1
c     iww=ispiel(l)
      iaddr(iww)=lz
      iaddr(iww)=jz+(ispiel(l)-1)*kmj
      jj=nqrfal(iww)*ii+nis
cvp
c     ip=ot(jj)
      ip=ntabl(jj)
      ipfeld(iww)=ip
c     if (ip.eq.255) ipfeld(iww)=-1
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        cbfeld(iww)=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(iww))) then
c         write(iwr,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(iww))
c       endif
        jto1=jto1+1
        exfeld(iww)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(intex(iww))) then
c         write(iwr,*) ' exchange-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intex(iww))
c       endif
       else
        cbfeld(iww)=twoe(intcb(iww))
        exfeld(iww)=-twoe(intex(iww))
      endif
 3013 lz=lz+kmj
cvp
c     write(iwr,*)' iww=',iww
      if (iww.ne.0) iwwz=iwwz+1
c aufbau der integrale zu allen konf. aus sk(i-2), lll bezeichnet
c  alle z-faelle
      do 5013 lll=1,3
      do 5013 l=1,iww
      bobvec(l,lll)=(ibob2(nrfal(l),lll)*cbfeld(l)+jbob2(nrfal(l),lll)
     &  *exfeld(l))*ipfeld(l)
 5013 continue
c jgcnt zaehlt die matrixelemente
      jgcnt=0
      do 13 l=1,iww
c write dk=2
cvp
      jj=nqrfal(l)
c     jj=1
cvp
c write dk=2
      jj1=jj
      jj=jj*ii+nis
      do 905 m=1,nzer
905   f(m)=0.0d0
c write dk=2
      do 18 m=1,kml
        my=ntabl(jj+m+m-1)
        if (my.ne.0) then
          kx=ja+my
          mike=ntabl(jj+m+m)
          tim=bobvec(l,mike)
          lx=m-kml
          do 22 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
  22        f(lx)=      tim*ae(kx)
c 22      f(lx)=f(lx)+tim*ae(kx)
c
      endif
18    continue
cvp  nqr = # qr-faelle
      nqr=ibinom(inopen,4)
cvp  irec = # eintraege von emat for ein saf-paar
      irec=6*nqr
cvp
      do 23 m=1,kml
        do 23 ih=1,kmj
          tim=0.0d0
          do 24 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
24          tim=tim+f(ly)*ae(mx)
cvp  matrixelemente nach neuer tabelle
c         iem=(ih-1)*kml*irec+(m-1)*irec
c         xtim=
c    &     emat(iem+(nrfal(l)-1)*2*nqr+2*(nqrfal(l)-1)+1)*cbfeld(l)
c    &    -emat(iem+(nrfal(l)-1)*2*nqr+2*(nqrfal(l)-1)+2)*exfeld(l)
c         xtim=xtim*ipfeld(l)
cvp   if (i.eq.3) then
c       if (dabs(xtim-tim).gt.1.0d-8)
c    &   write(iwr,*) ' l=,m=,if=,tim=,xtim= '
c    &      ,l,m,if,tim,xtim
c     if (nqrfal(l).eq.1)
c    & write(iwr,*) ' m=,if=,nrfal=,cbfeld=,exfeld=,tim= '
c    &  ,m,if,nrfal(l),cbfeld(l),exfeld(l),tim
cvp   endif
cvp jgcnt kann kml*kmj mal so gross wie die anzahl der
cvp  konfigurationen werden ---> bei jgcnt=ncmax matrixelemente
cvp  wegschreiben und jgcnt auf null zurueckstetzen
cdir      if (dabs(tim).ge.1.0e-7) then
            jgcnt=jgcnt+1
            qfeld(jgcnt)=tim
            mifeld(jgcnt)=kz+m
            mjfeld(jgcnt)=iaddr(l)+ih
cdir      endif
c schreiben der matrixelemente und der indices auf file nhuk
      if (jgcnt.eq.ncmax) then
cvocl loop,novrec,vi(ik)
c
c schleife ueber die wurzeln
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
        write(iwr,*)' 5100 : ipoint=',ipoint
c
         do 4012 ik=1,jgcnt
           qq=qfeld(ik)
           iih=mifeld(ik)
           jjh=mjfeld(ik)
           vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*qq
           vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*qq
 4012    continue
c ende schleife uber die wurzeln
 5100 continue
         jgcnt=0
      endif
cvp
  23  continue
  13  continue
      iwwg=iwwg+iww
c schreiben der matrixelemente und der indices auf file nhuk
cvocl loop,novrec,vi(ik)
c
c schleife ueber die wurzeln
c
      do 5200 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
      ipoint=(iwurz-1)*ndim
c
      do 4013 ik=1,jgcnt
        qq   =qfeld(ik)
        iih   =mifeld(ik)
        jjh   =mjfeld(ik)
        vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*qq
        vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*qq
 4013 continue
c ende schleife uber die wurzeln
 5200 continue
cvp
 3012 continue
cvp
  12  kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=2
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1012
cio   endif
cvp
 1014 continue
cvp
c ende dk=2
cdebugwrite(iwr,*)' jg=',jg
cdebugwrite(iwr,*)' iwwg=',iwwg
cdebugwrite(iwr,*)' iwwz=',iwwz
cvp
      go to 10
csut <--- reinsprung
csut701   continue
9     if (i.eq.1) go to 25
10    ibl=ibl+1
*      write(iwr,*) 'dk=1 sk=',i,'ibl = ',ibl
c//skipc&2
c     call cpu
c     call setcpu
c start dk=1 fall
c beginn timer-routine for dk=1
*     call clockv(vu1,cpu1)
*     call time(milli1)
c
      j8=i-1
      mc=nconf(j8)
csut  if (i.le.issk) then
csut     mcx55 = 0
csut  else
         mcx55 = mc
csut  end if
      if (mcx55.eq.0) go to 25
      js=jan(ibl)
      ks=jbn(ibl)
      ndj=ndet(j8)
      kmj=nsac(j8)
      ja=iaz(j8)
      jz=isc(j8)
      nzer=kml*kmj
      mr=nod(j8)
      mps=nplu(j8)
      mms=j8-1
c write dk=1
c     write(iwr,1101)j8,mc,js,ks,ndj,kmj,ja,mr,mps,mms
*1101 format(/1x,'dk=1 (1101):j8,mc,js,ks,ndj,kmj,ja,mr,mps,mms',
*    *i3,2x,i5,2x,8(i6,2x)/)
      ix=nad1
      kss=0
c      write(6,92349) js, ks, nss
c92349 format(1x,'js, ks, nss = ', 3i10)
 27   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 27
      endif
      call upackx(ot)
c     write(6,92348) ic, i2, i3
c92348 format(1x,'ic, i2, i3 =', 3i10)
      jc=nr+mult
      j3=jc*kml+1
      kz=iz
csut --- raussprung
csut  if (i.le.issk) goto 702
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,3)*(inopen-2)
     &     +kml*kmj*ibinom(inopen,2)
     &     +2*kml*kmj*ibinom(inopen,2)*(inopen-1)
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1242) npdim, nedim
*1242    format(' ***calculation of matrix elements for dk=1 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
         endif
         call sgdk1(emat,nedim,ae,ot,ja,ie,kml,kmj
     &             ,ndj,inopen,ic,i2
     &             ,j3,i3,mps,mms,nms,nps,jc)
      else
       if (odebug) then
         write(6,1243) npdim, nedim
1243    format(' **** sgdk1:  npdim.gt.nedim', 2i7)
       endif
      endif
c
cvp
cvp lesen der integrale for dk=1, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=1
cio   kend=nc
cvp ruecksprungadresse, falls dk=1 noch nicht fertig
c1028 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(mc-1)/ncmax+1
cvp
cio   do 28 k=kanf,kend
csut  do 28 k=1,ncx55
      do 28 k=1,nc
cvp
cvp zerlegung der schleife ueber alle mc konfigurationen in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3028 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=mc-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c mit checkb werden alle wechselwirkenden konf. der sk(i-1) bestimmt
c sowie die labels
      call check1a(twoe,nteint,iottr,iotnew,iotm,
     +   maindf,maxr,jconb,nomax2,
     +   ioint,i,j8,k,nspiel,mspiel,jnum1
     &   ,nspp1,nsppa,nsppb)
c
cvp
c     if (nqrfal(l).ne.ndfal(l).and.npfal(l).eq.3) then
c       write(iwr,*) ' k=,l=,npfal,nqrfal=,ndfal= '
c    &      ,k,l,npfal(l),nqrfal(l),ndfal(l)
c     endif
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kmj*ibinom(inopen,3)*(inopen-2)
     &     +kml*kmj*ibinom(inopen,2)
     &     +2*kml*kmj*ibinom(inopen,2)*(inopen-1)
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1244) npdim, nedim
*1244     format(' ***calculation of matrix elements for dk=1 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
        call hmdk1(twoe,pey,
     +              vecf1,vecf2,mdi,
     +              emat,nedim,i,nc,mc,kml,kmj,inopen,jz,kz
     &             ,nhuk,mtype,jto1,jg,ioint,iges,kend,nspp1
     &             ,nsppa,nsppb)
         goto 1027
      else
       if (odebug) then
        write(6,1245) npdim, nedim
1245    format(' **** hmdk1:  npdim.gt.nedim', 2i7)
       endif
      endif
c dk=1 p=1
      do 31 l=1,nspp1
c welcher r-fall ll=1-3 groesser null,ll=4-6 kleiner null
      ll=nrfal(l)
c ql-fall
      jj=nqlfal(l)
c qr-fall
      kk=nqrfal(l)
      jj=(kk-1)*mr +jj
      jj=jj*ii+ic
      ip=ot(jj)
      if(ip.eq.255) ll=-ll
      if (ll.lt.0) go to 820
c     if (exc.ne.-twoe(intex(l))) then
c       write(iwr,*) ' k=,l=,nrfal=,exc=,twoe= ',
c    &   k,l,nrfal(l),exc,twoe(intex(l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul= twoe(jto1)
        jto1=jto1+1
        exc =-twoe(jto1)
       else
        coul=twoe(intcb(l))
        exc=-twoe(intex(l))
      endif
cvp
      go to 38
820   continue
      ll=-ll
cvp   exc=h(jto)
cvp
c     if (exc.ne.twoe(intex(l))) then
c       write(iwr,*) ' k=,l=,nrfal=,exc=,twoe= ',
c    &   k,l,nrfal(l),exc,twoe(intex(l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul=-twoe(jto1)
        jto1=jto1+1
        exc = twoe(jto1)
       else
        coul=-twoe(intcb(l))
        exc=twoe(intex(l))
      endif
cvp
38    continue
cvp lesen der integrale
c     if (ioint.eq.1) then
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intcb(l))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(l))
c       endif
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intex(l))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' exchange-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intex(l))
c       endif
c     endif
      if (ll-2)810,811,812
810   bob(1)=-exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 49
811   bob(1)=exc
      bob(2)=-coul-exc
      bob(3)=coul
      go to 49
812   bob(1)=coul+exc
      bob(2)=-exc
      bob(3)=-coul
   49 continue
      do 39 m=1,nzer
39    f(m)=0.0d0
      do 40 m=1,kml
        jj=jj+1
        my=ot(jj)
        jj=jj+1
        if(my.ne.0) then
          kx=ja+my
          mike=ot(jj)
          if(mike.gt.128) then
             tim=-bob(256-mike)
           else
             tim=bob(mike)
          endif
          lx=m-kml
          do 44 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 44        f(lx)=f(lx)+tim*ae(kx)
  44        f(lx)=      tim*ae(kx)
        endif
40    continue
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trfdk1(vecf1, vecf2, mdi,
     +            f,ae,kz,ie,nz,kml,kmj)
cvp
 31   continue
cvp
c p=2,r=1,2
      do 32 l=nspp1+1,nsppa
c qr vorzeichen kodiert r fall
      ll=nqrfal(l)
cvp
c     if (sm.ne.twoe(intcb(l))) then
c       write(iwr,*) ' k=,l=,nrfal=,sm=,twoe= ',
c    &   k,l,nrfal(l),sm,twoe(intcb(l))
c     endif
      if (ioint.eq.1) then
        jto1=jto1+1
        sm=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(l))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' Coulomb-integrals: k=,l= ',k,l,twoei(jto1)
c    &    ,twoe(intcb(l))
c       endif
       else
        sm=twoe(intcb(l))
      endif
cvp
      if(ll.gt.0) go to 53
      sm=-sm
      ll=-ll
53    jj=ll*kml+i2
      do 55 m=1,nzer
55    f(m)=0.0d0
      do 56 m=1,kml
        jj=jj+1
c welche determinante ww mit der m-ten sdet, vorzeichen kodiert ip
        my=ot(jj)
        if (my.ne.0) then
          if(my.lt.128) then
             tim=sm
             kx=my+ja
           else
             kx=ja-my+256
             tim=-sm
          endif
          lx=m-kml
          do 60 if=1,kmj
            kx=kx+ndj
            lx=lx+kml
c 60        f(lx)=f(lx)+tim*ae(kx)
  60        f(lx)=      tim*ae(kx)
        endif
56    continue
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trfdk1(vecf1, vecf2, mdi,
     +            f,ae,kz,ie,nz,kml,kmj)
cvp
  32  continue
cvp
c dk=1,p=3,r=1,2
      do 33 l=nsppa+1,nsppb
      ll=nqrfal(l)
      ntar=nrfal(l)
      nb=moafal(l)
      mb=mobfal(l)
      nb=ideks(nb)+mb+ijone(ntar)
cvp   nb=ideks(nb)+mb+ij(ntar)
c     sm=pey(nb)
      sm1=pey(nb)
c lesen der summe ueber die dreizentrenintegrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(iwr,*) ' p=3, jto1= ',jto1
c         write(iwr,*) ' integral sum: k=,l= ',k,l,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(l)
      endif
      sm=tim1
cvp
c write dk=1,p=3
      if (mr.eq.0 ) go to 65
c lesen der (negativen) austauschintegrale
      do 66 m=1,mr
cvp
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(m)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(l,m))) then
c         write(iwr,*) ' nsppb=,jto1= ',nsppb,jto1
c         write(iwr,*) ' integrale: k=,l=,m= ',k,l,m,twoei(jto1)
c    &    ,twoe(nwwmo(l,m))
c       endif
       else
        sac(m)=-twoe(nwwmo(l,m))
      endif
cvp
c     if (sac(m).ne.-twoe(nwwmo(l,m))) then
c       write(iwr,*) ' k=,l=,m=,sac=,-nwwmo= ',k,l,m,sac(m)
c    &   ,-twoe(nwwmo(l,m))
c     endif
cvp
66    continue
   65 continue
c     iiiiii=iiiiii+1
c     if (iiiiii.lt.100)
c    &write(iwr,*)' jto=,sm=,sac=,sm1=',jto,sm,sac,sm1
      sm=sm+sm1
c abzweigung r=1 bzw 2
cvp qr < 0 --> r=1; qr > 0 --> r=2
      if(ll.lt.0) go to 86
      jj=ll*j3+i3
      ip=ot(jj)
      if(ip.eq.1) go to 71
      sm=-sm
cvp mr = # der offenen schalen der kleineren sk
      if(mr.eq.0) go to 71
      do 72 m=1,mr
72    sac(m)=-sac(m)
71    jj=jj+1
      do 80 m=1,nzer
80    f(m)=0
      jx=-kml
      do 73 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (74,75,76,77),im
 74      in=in+1
         kx=ot(in)
         kx=kx+ja
         tim=sm
cvp nps = # alpha-spins der offenen schalen der kleineren sk
         if(mps.eq.0) go to 78
         do 79 ih=1,mps
           in=in+1
           im=ot(in)
79         tim=tim+sac(im)
78       lx=jx
         do 81 ih=1,kmj
           kx=kx+ndj
           lx=lx+kml
81         f(lx)=f(lx)+tim*ae(kx)
        go to 73
 75     in=in+1
        kx=ot(in)
        kx=kx+ja
        tim=-sm
        if(mms.eq.0) go to 78
cvp mms = nr. der kleineren sk - 1 ; bzw.
cvp mms = #  beta-spins der offenen schalen der kleineren sk
        do 83 ih=1,mms
          in=in+1
          im=ot(in)
83        tim=tim-sac(im)
        go to 78
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
76      do 84 ih=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=-sac(im)
          lx=jx
          do 84 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
84          f(lx)=f(lx)+tim*ae(kx)
        go to 73
cvp nps = # alpha-spins der offenen schalen
77      do 85 ih=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ja
          in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 85 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
85          f(lx)=f(lx)+tim*ae(kx)
   73 jj=jj+jc
      go to 47
c p=3 anderer r fall (r=1)
  86  jj=-1*ll*j3+i3
      ip=ot(jj)
      if (ip.eq.1) go to 87
      sm=-sm
      if (mr.eq.0) go to 87
      do 88 m=1,mr
88    sac(m)=-sac(m)
87    jj=jj+1
      do 89 m=1,nzer
89    f(m)=0
      jx=-kml
      do 90 m=1,kml
        jx=jx+1
        in=jj
        im=ot(in)
        go to (91,95,96,97),im
91        in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=sm
          if(mms.eq.0) go to 92
          in=in+mps
          do 93 ih=1,mms
            in=in+1
            im=ot(in)
93          tim=tim+sac(im)
92        lx=jx
          do 94 ih=1,kmj
            kx=kx+ndj
            lx=lx+kml
94          f(lx)=f(lx)+tim*ae(kx)
          go to 90
95        in=in+1
          kx=ot(in)
          kx=kx+ja
          tim=-sm
          if(mps.eq.0) go to 92
          in=in+mms
          do 99 ih=1,mps
            in=in+1
            im=ot(in)
99          tim=tim-sac(im)
          go to 92
96        do 100 ih=1,nms
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
          im=ot(in)
          tim=sac(im)
          lx=jx
          do 100 ig=1,kmj
            kx=kx+ndj
            lx=lx+kml
100         f(lx)=f(lx)+tim*ae(kx)
          go to 90
97        do 101 ih=1,nps
            in=in+1
            kx=ot(in)
            kx=kx+ja
            in=in+1
            im=ot(in)
            tim=-sac(im)
            lx=jx
            do 101 ig=1,kmj
               kx=kx+ndj
               lx=lx+kml
101            f(lx)=f(lx)+tim*ae(kx)
90         jj=jj+jc
cre
   47 continue
cvp endtransformation for dk=1
      nz=jz+(ispiel(l)-1)*kmj
      call trfdk1(vecf1, vecf2, mdi,
     +            f,ae,kz,ie,nz,kml,kmj)
cvp
  33  continue
cvp
 1027 continue
cvp
 3028 continue
cvp
28    kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=1
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1028
cio   endif
cvp
25    js=idra(i)
c
c           *******   start dk=0   *******
c
c start dk=0 p=5 fall
c beginn timer-routine for dk=0 alle p-faelle
c beginn timer-routine for dk=0,p=5
cdebugcall clockv(vu1,cpu1)
cdebugcall time(milli1)
c
      ks=idrc(i)
      ix=nad1
c write dk=0 p=5
c     write(jtape,1500)js,ks
c1500 format(1x,'ks,js for dk=0 p=5 (1500)',2i5)
      kss=0
102   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 102
      endif
*****
      call upackx(ot)
      kz=iz
cvp
      ig=0
cvp
csut  if (i.le.issk) goto 6
csut  jstar=niot(i)-ncx55
      jstar=niot(i)-nc
csut  do 103 k=1,nc
      do 103 k=1,ncx55
cvp nx = # besetzte mo's
      do 104 l=1,nx
csut     kstar=jstar+nc*l
         kstar=jstar+ncx55*l
         itest(l)=iottr(kstar+k)
 104  continue
c write dk=0 p=5
      do 105 l=1,mzer
105   f(l)=0.0d0
      care=core
cvp nr = # einfach besetzte mo's
      do 107 l=1,nr
         mal=itest(l)
         ntar=nirred(mal)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
107      care=care+pey(mal)
cvp nd = # doppelt besetzte mo's
cvp n1 = nr+1
      do 109 l=n1,nx
         mal=itest(l)
         ntar=nirred(mal)
         nal=ideks(mal+1)
         mal=mopos(mal)+1
         mal=ideks(mal)+ijone(ntar)
         sm=pey(mal)
109      care=care+sm+sm+acoul(nal)
      do 111 l=2,nr
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 111 m=1,it
           nal=mal+itest(m)
111        care=care+acoul(nal)
      sm=0.0d0
cvp n2 = nr + 2
      do 113 l=n2,nx
         mal=itest(l)
         mal=ideks(mal)
         it=l-1
         do 113 m=n1,it
           nal=mal+itest(m)
           tim=acoul(nal)
113        sm=sm+tim+tim-aexc(nal)
      care=care+sm+sm
      do 115 l=1,nr
         mal=itest(l)
         do 115 m=n1,nx
           nal=itest(m)
           nal=ideks(max(nal,mal))+min(nal,mal)
           sm=acoul(nal)
115        care=care+sm+sm-aexc(nal)
      kk=ic
      kn=i2
      ml=i3
c write dk=0 p=5
      jx=-kml
      do 118 l=1,kml
        jx=jx+1
        lx=jx
        kx=ohog(l)+ia
        tim=care
cvp nps = # alpha-spins der offenen schalen
c       if(nps.lt.2) go to 122
        do 123 ih=2,nps
          mx=i2+nps*(l-1)+ih
          it=ih-1
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          do 123 in=1,it
            mz=i2+nps*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
123         tim=tim-aexc(nal)
cvp nms = nr. der sk - 1 ; bzw.
cvp nms = #  beta-spins der offenen schalen
c       if(nms.lt.2) go to 122
        do 124 ih=2,nms
          mx=i3+nms*(l-1)+ih
          mal=out(mx)
          mal=itest(mal)
          mal=ideks(mal)
          it=ih-1
          do 124 in=1,it
            mz=i3+nms*(l-1)+in
            nal=out(mz)
            nal=itest(nal)+mal
124         tim=tim-aexc(nal)
        do 121 ih=1,kml
          lx=lx+kml
          kx=kx+ndt
121       f(lx)=f(lx)+tim*ae(kx)
cvp niw = nps * nms
        do 120 m=1,niw
          kk=kk+1
          kx=out(kk)+ia
c write dk=0 p=5 ausserdiagonal
          kk=kk+1
          mal=out(kk)
          kk=kk+1
          nal=out(kk)
          mal=itest(mal)
          nal=ideks(mal)+itest(nal)
          tim=-aexc(nal)
          lx=jx
          do 120 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
120         f(lx)=f(lx)+tim*ae(kx)
118   continue
      do 125 l=1,kml
        mz=kz+l
        do 125 m=1,l
          tim=0.0d0
          do 126 ih =1,kml
            ly=(m-1)*kml+ih
            mx=ie+(l-1)*kml+ih
126         tim=tim+f(ly)*ae(mx)
          if(l.gt.m) go to 127
          jg=jg+1
cx        q(jg)=tim
cvp diagonale wird jetzt aus hmp5 genommen
c         q(jg)=hp5(mz)
c         tim=hp5(mz)
c         iih=mz
c         jjh=iih
c --- diagonale !!!
c
c schleife ueber die wurzeln
c
c     do 5300 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
c      ipoint=(iwurz-1)*ndim
*        vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
c ende schleife uber die wurzeln
c5300 continue
          go to 125
 127      continue
          jjh   =kz+m
      do 5400 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
          vecf1(ipoint+mz)=vecf1(ipoint+mz)+vecf2(ipoint+jjh)*tim
          vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+mz)*tim
c ende schleife uber die wurzeln
 5400 continue
125   continue
103   kz=mz
      if(nc.eq.1) go to 6
c start dk=0 andere p-faelle
      ks=jdrc(i)
      ix=nad1
c write dk=0 andere p-faelle
c     write(jtape,1520)js,ks
c1520 format(1x,'dk=0 other p cases js,ks',2(i7,2x))
      kss=0
128   if(kss.ne.ks) then
      ix=ix+nad
      call rdbak3(js)
      js=js+nsz
      call stopbk3
      call unpack(lout,8,ot(ix),nad)
      kss=kss+1
      go to 128
      endif
*****
      call upackx(ot)
      ii=ii+kml
      j2=kml+1
      jc=nr+nr
      j3=jc*kml+1
      kz=iz+kml
csut --- raussprung
csut  if (i.le.issk) goto 6
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kml*ideks(ibinom(inopen,2)+1)
     &     +kml*kml*ideks(inopen+1)
     &     +kml*kml*ideks(inopen+1)*2*inopen
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1246) npdim, nedim
*1246     format(' ***calculation of matrix elements for dk=0 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call sgdk0(emat,nedim,ae,ot,ie,kml,
     &              ndt,inopen,ic,i2,ii,ia,j2
     &             ,j3,i3,nps,np1,nms,nm1,jc,np2,nm2)
      else
       if (odebug) then
        write(6,1247) npdim, nedim
1247    format(' **** sgdk0:  npdim.gt.nedim', 2i7)
       endif
      endif
c die untere dreickmatrix wird berechnet
      j8=0
*     extimp=0.0d0
*     vutimp=0.0d0
cdebugcputip=0.0d0
cvp zeiten for bestimmung der wechselwirkenden konfigurationen
cdebugvutwwg=0.0d0
*     cpuwwg=0.0d0
cvp zeiten for labelbestimmung
*     vutlag=0.0d0
*     cpulag=0.0d0
cvp zeiten for integraladressen
*     vuting=0.0d0
*     cpuing=0.0d0
cvp
cvp
cvp lesen der integrale for dk=0, falls ioint = 1
cvp  iges = anzahl der integrale auf twoe
cvp  kend = anzahl der konfigurationen der sk(i), for die integrale
cvp         abgespeichert sind
cio   kanf=2
cio   kend=nc
cvp ruecksprungadresse, falls dk=0 noch nicht fertig
c1129 continue
cio   if (ioint.eq.1) then
cio      read(mtype) iges,kend,(twoe(iii),iii=1,iges)
cio      jto1=0
cio   endif
cvp
cio   do 129 j=kanf,kend
csut  do 129 j=2,nc
      do 129 j=2,ncx55
cvp   lz=iz
cvp
cvp berechnung der anzahl nstrip der schleifendurchlaeufe ueber die
cvp  konfigurationen der kleineren sk, falls deren anzahl groesser
cvp  als ncmax ist
      nstrip=(j-2)/ncmax+1
cvp
cvp
cvp zerlegung der schleife ueber alle j-1 konfigurationen
cvp  der unteren dreiecksmatrix in
cvp  mehrere mit maximaler laenge ncmax
cvp mspiel = anzahl zu testender konfigurationen
      mspiel=ncmax
      do 3129 istrip=1,nstrip
cvp jnum1 = nummer der ersten zu testenden konfiguration -1
        jnum1=(istrip-1)*ncmax
        if (istrip.eq.nstrip) then
           mspiel=j-1-(nstrip-1)*ncmax
        endif
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
         if (jto1.eq.iges) then
            read(mtype) iges,kend,(twoe(iii),iii=1,iges)
            jto1=0
         endif
      endif
c mit checkd werden alle wechselwirkenden konf. der sk(i) bestimmt
c sowie die labels
cdebugcall clockv(vu10,cpu10)
cdebugcall time(mill10)
csut  call checkd(ioint,i,i,j,nspiel,mspiel,jnum1
      call check0a(twoe,nteint,iottr,iotnew,iotm,
     +  maindf,maxr,jconb,nomax2,
     +  ioint,i,i,j,nspiel,mspiel,jnum1
     & ,nspp1,nsppa,nsppb,nspp4)
cvp
c
cvp   j8=j8+1
cvp
c     if (moafal(k).ne.ndfal(k).and.npfal(k).ge.4) then
c       write(iwr,*) ' j=,k=,npfal=,moafal=,ndfal= ',
c    &   j,k,npfal(k),moafal(k),ndfal(k)
c       if (nqrfal(k).lt.0) then
c         write(iwr,*) ' r= 2'
c        else
c         write(iwr,*) ' r= 1'
c       endif
c     endif
cvp
cvp falls die darstellungsmatrizen auf emat passen, werden sie
cvp  hier mit hilfe der bunkerschen table generiert
cvp  berechnung der gesamtgroesse der darstellungsmatrizen:
      npdim=6*kml*kml*ideks(ibinom(inopen,2)+1)
     &     +kml*kml*ideks(inopen+1)
     &     +kml*kml*ideks(inopen+1)*2*inopen
      if (npdim.le.nedim) then
       if (odebug) then
*         write(6,1248) npdim, nedim
*1248     format(' ***calculation of matrix elements for dk=0 via sga '/
*     +        ' npdim = ', i6, ' : nedim = ', i7)
       endif
         call hmdk0(twoe,pey,aexc,
     +              vecf1,vecf2,mdi,
     +              emat,nedim,i,nc,mc,kml,kmj,inopen,iz,kz
     &             ,nhuk,mtype,jto1,jg,ioint,iges,kend,nspp1
     &             ,nsppa,nsppb,nspp4)
         goto 1129
      else
       if (odebug) then
        write(6,1249) npdim, nedim
1249    format(' **** hmdk0:  npdim.gt.nedim', 2i7)
       endif
      endif
c dk=0 p=1
      do 132 k=1,nspp1
c r-fall
      ll=nrfal(k)
cvp qr-fall
      jj=nqrfal(k)
cvp ql-fall
      kk=nqlfal(k)
      if(jj.lt.kk) go to 150
      jj=ideks(jj)+kk
      ibob=0
      go to 151
150   jj=ideks(kk)+jj
      ibob=1
151   jj=jj*ii+ic
c die allgemeine permutation
      ip=ot(jj)
c write dk=0,p=1
      if(ip.eq.255) go to 830
cvp   coul=h(jto)
cvp
c     if (coul.ne.twoe(intcb(k))) then
c       write(iwr,*) ' j=,k=,nrfal=,coul=,twoe= ',
c    &   j,k,nrfal(k),coul,twoe(intcb(k))
c     endif
cvp
c     if (exc.ne.-twoe(intex(k))) then
c       write(iwr,*) ' j=,k=,nrfal=,exc=,twoe= ',
c    &   j,k,nrfal(k),exc,twoe(intex(k))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul= twoe(jto1)
        jto1=jto1+1
        exc =-twoe(jto1)
       else
        coul=twoe(intcb(k))
        exc=-twoe(intex(k))
      endif
cvp
      go to 153
 830  continue
c830  coul=-h(jto)
cvp   if (coul.ne.-twoe(intcb(k))) then
cvp     write(iwr,*) ' j=,k=,nrfal=,coul=,twoe= ',
cvp  &   j,k,nrfal(k),coul,twoe(intcb(k))
cvp   endif
cvp
c     if (exc.ne.twoe(intex(k))) then
c       write(iwr,*) ' j=,k=,nrfal=,exc=,twoe= ',
c    &   j,k,nrfal(k),exc,twoe(intex(k))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        coul=-twoe(jto1)
        jto1=jto1+1
        exc = twoe(jto1)
       else
        coul=-twoe(intcb(k))
        exc=twoe(intex(k))
      endif
cvp
  153 continue
cvp lesen der integrale
c     if (ioint.eq.1) then
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intcb(k))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' Coulomb-integrals: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intcb(k))
c       endif
c       jto1=jto1+1
c       if (twoei(jto1).ne.twoe(intex(k))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' exchange-integrals: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intex(k))
c       endif
c     endif
      if (ll-2) 832,833,834
  832 bob(1)=coul+exc
      bob(2)=exc
      bob(3)=coul
      go to 835
  833 bob(1)=exc
      bob(2)=-coul
      bob(3)=coul+exc
      go to 835
 834  bob(1)=-exc
      bob(2)=-coul-exc
      bob(3)=coul
835   do 154 l=1,mzer
154   f(l)=0
      mx=-kml
      do 155 l=1,kml
        jj=jj+1
        mx=mx+1
        ip=ot(jj)
        if(ip.eq.2)  go to 159
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        jj=jj+1
        tim=bob(1)
        lx=mx
        do 161 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c161      f(lx)=f(lx)+tim*ae(kx)
 161      f(lx)=      tim*ae(kx)
        go to 155
159     jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=bob(2)
        lx=mx
        do 160 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c160      f(lx)=f(lx)+tim*ae(kx)
 160      f(lx)=      tim*ae(kx)
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        tim=bob(3)
        lx=mx
        do 167 m=1,kml
          lx=lx+kml
          kx=kx+ndt
167       f(lx)=f(lx)+tim*ae(kx)
155   continue
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdk0(vecf1,vecf2,mdi,
     +            f,ae,kz,ie,nz,kml,ibob)
cvp
 132  continue
cvp
c dk=0,p=2
      do 133 k=nspp1+1,nsppa
cvp qr-fall
      ll=nqrfal(k)
cvp ql-fall
      jj=nqlfal(k)
      if(jj.gt.ll) go to 170
      jj=ideks(ll)+jj
      ibob=0
      go to 171
170   jj=ideks(jj)+ll
      ibob=1
  171 continue
cvp
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim=twoe(jto1)
c       if (twoei(jto1).ne.twoe(intcb(k))) then
c         write(iwr,*) ' jto1= ',jto1
c         write(iwr,*) ' Coulomb-integrals: j=,k= ',j,k,twoei(jto1)
c    &    ,twoe(intcb(k))
c       endif
       else
        tim=twoe(intcb(k))
      endif
cvp
c write dk=0 p=2
      jj=jj*j2+i2
      do  l=1,mzer
        f(l)=0.0d0
      end do
      jx=-kml
      ip=ot(jj)
      if(ip.eq.255) tim=-tim
      do 173 l=1,kml
        jx=jx+1
        jj=jj+1
        kx=ot(jj)
        kx=kx+ia
        lx=jx
        do 173 m=1,kml
          kx=kx+ndt
          lx=lx+kml
c 173     f(lx)=f(lx)+tim*ae(kx)
  173     f(lx)=      tim*ae(kx)
cvp
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdk0(vecf1,vecf2,mdi,
     +            f,ae,kz,ie,nz,kml,ibob)
cvp
 133  continue
cvp
c dk=0 p=3
      do 134 k=nsppa+1,nsppb
c qr, vorzeichen kodiert den r-fall
cvp qr > 0 --> r=1
      ll=nqrfal(k)
c ql
      jj=nqlfal(k)
c symmetrie der orbitale for einelektronen-integral
      ntar=nrfal(k)
c 1.mo for einelektronen-integral (das groessere)
      nb=moafal(k)
c 2.mo for einelektronen-integral (das kleinere)
      mb=mobfal(k)
      nb=ideks(nb)+mb+ijone(ntar)
cvp einelektronenintegral wird separat behandelt
c     sm=pey(nb)
      sm1=pey(nb)
      sm=0.0d0
c abzweigung zu den unterschiedlichen r-faellen
      if(ll.lt.0) go to 181
c dk=0 p=3 r=1
      if(ll.lt.jj) go to 182
      ibob=0
      jj=ideks(ll)+jj
c write dk=0 p=3 r=1
      go to 205
182   jj=ideks(jj)+ll
c write dk=0 p=3 r=1
      ibob=1
205   jj=jj*j3+i3
c lesen der summe ueber die dreizentrenintegrale
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(iwr,*) ' r=1,jto1= ',jto1
c         write(iwr,*) ' integral sum: j=,k= ',j,k,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(k)
      endif
      sm=tim1
cvp
      if(nr.eq.1) go to 185
c lesen der negativen austauschintegrale
cvp n3 = nr - 1
      do 186 l=1,n3
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(l)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(k,l))) then
c         write(iwr,*) ' r=1,jto1= ',jto1
c         write(iwr,*) ' single integrals: j=,k=,l= ',j,k,l,twoei(jto1)
c    &    ,twoe(nwwmo(k,l))
c       endif
       else
        sac(l)=-twoe(nwwmo(k,l))
      endif
cvp
c write dk=0 p=3 r=1
186   continue
  185 continue
c einelektronenintegral wird addiert
      sm=sm+sm1
      ip=ot(jj)
      if(ip.eq.1) go to 191
      sm=-sm
      if (nr.eq.1) go to 191
      do 192 l=1,n3
192   sac(l)=-sac(l)
191   jj=jj+1
      do 193 l=1,mzer
193   f(l)=0
      jx=-kml
      do 194 l=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=sm
        if(im.eq.2) go to 195
        if(nps.eq.1) go to 196
        do 197 m=1,np1
          in=in+2
          nn=ot(in)
197       tim=tim+sac(nn)
196     do 198 m=1,kml
          lx=lx+kml
          kx=kx+ndt
c198      f(lx)=f(lx)+tim*ae(kx)
 198      f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 194
        do 199 m=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 199 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
199         f(lx)=f(lx)+tim*ae(kx)
        go to 194
195     if(nms.eq.1) go to 200
        do 201 m=1,nm1
          in=in+2
          nn=ot(in)
201       tim=tim+sac(nn)
200     do 202 m=1,kml
          lx=lx+kml
          kx=kx+ndt
202       f(lx)=f(lx)+tim*ae(kx)
        do 203 m=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 203 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
203         f(lx)=f(lx)+tim*ae(kx)
194     jj=jj+jc
cvp   if (ibob.eq.1) go to 144
      go to 146
c dk=0 p=3 r=2
  181 ll=-1*ll
      if(ll.lt.jj) go to 204
      ibob=0
      jj=ideks(ll)+jj
c write dk=0 p=3 r=2
      jj1=jj
      go to 206
204   jj=ideks(jj)+ll
      ibob=1
206   jj=jj*j3+i3
c erstes integral for r=2 siehe aufzeichnung
c lesen der summe der dreizentrenintegrale
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        tim1=twoe(jto1)
c       if (twoei(jto1).ne.tim1) then
c         write(iwr,*) ' r=2,jto1= ',jto1
c         write(iwr,*) ' integral sum: j=,k= ',j,k,twoei(jto1)
c    &    ,tim1
c       endif
       else
        tim1=sumint(k)
      endif
cvp
      sm=sm+tim1
      if(nr.eq.1) go to 209
c do 211 ist analog do 186
      do 211 l=1,n3
cvp
c     sac(l)=h(jto)
c     if (sac(l).ne.-twoe(nwwmo(k,l))) then
c       write(iwr,*) ' j=,k=,l=,sac=,-nwwmo= ',j,k,l,sac(l)
c    &   ,-twoe(nwwmo(k,l))
c     endif
cvp lesen der integrale
      if (ioint.eq.1) then
        jto1=jto1+1
        sac(l)=-twoe(jto1)
c       if (twoei(jto1).ne.twoe(nwwmo(k,l))) then
c         write(iwr,*) ' r=2,jto1= ',jto1
c         write(iwr,*) ' single integral: j=,k=,l= ',j,k,l,twoei(jto1)
c    &    ,twoe(nwwmo(k,l))
c       endif
       else
        sac(l)=-twoe(nwwmo(k,l))
      endif
cvp
c write dk=0 p=3 r=2
211   continue
  209 continue
c einelektronenintegrale werden addiert
      sm=sm+sm1
      ip=ot(jj)
      if(ip.eq.1) go to 215
      sm=-sm
      if(nr.eq.1) go to 215
      do 216 l=1,n3
216   sac(l)=-sac(l)
215   jj=jj+1
      do 217 l=1,mzer
217   f(l)=0
      jx=-kml
      do 218 l=1,kml
        in=jj
        im=ot(in)
        jx=jx+1
        in=in+1
        kx=ot(in)
        kx=kx+ia
        lx=jx
        tim=-sm
        if (im.eq.2) go to 219
        if(nms.eq.0) go to 220
        in=in+np2
        jn=in
        do 221 m=1,nms
          jn=jn+2
          nn=ot(jn)
221       tim=tim-sac(nn)
220     do 222 m=1,kml
          lx=lx+kml
          kx=kx+ndt
c 222     f(lx)=f(lx)+tim*ae(kx)
  222     f(lx)=      tim*ae(kx)
        if(nms.eq.0) go to 218
        do 223 m=1,nms
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 223 ih =1,kml
            lx=lx+kml
            kx=kx+ndt
223         f(lx)=f(lx)+tim*ae(kx)
        go to 218
219     in=in+nm2
        jn=in
        do 224 m=1,nps
          jn=jn+2
          nn=ot(jn)
224       tim=tim-sac(nn)
        do 225 m=1,kml
          lx=lx+kml
          kx=kx+ndt
225       f(lx)=f(lx)+tim*ae(kx)
        do 226 m=1,nps
          in=in+1
          kx=ot(in)
          kx=kx+ia
          lx=jx
          in=in+1
          mx=ot(in)
          tim=sac(mx)
          do 226 ih=1,kml
            lx=lx+kml
            kx=kx+ndt
226         f(lx)=f(lx)+tim*ae(kx)
218     jj=jj+jc
cvp
 146  continue
      nz=iz+(ispiel(k)-1)*kml
cvp endtransformation for dk=0
      call trfdk0(vecf1,vecf2,mdi,
     +            f,ae,kz,ie,nz,kml,ibob)
cvp
 134  continue
cvp
c dk=0 p=4     aa/bb
      do 135 k=nsppb+1,nsppb+nspp4
c missmatch in l
cvp absolutnummer von mo a
      kk=mobfal(k)
c missmatch in r
cvp absolutnummer von mo b
      ll=moafal(k)
      if(kk.lt.ll) go to 138
      kk=ideks(kk)+ll
      go to 139
138   kk=ideks(ll)+kk
139   continue
      tim=aexc(kk)
      do 141 l=1,kml
cvp
          iih=kz+l
          jjh=iz+(ispiel(k)-1)*kml+l
c
c schleife ueber die wurzeln
c
      do 5500 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
          vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
          vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*tim
c ende schleife uber die wurzeln
 5500 continue
  141 continue
 135  continue
cvp
 1129 continue
cvp
 3129 continue
cvp
129   kz=kz+kml
cvp
cvp ruecksprung zum abarbeiten des rests von dk=0
cio   if (ioint.eq.1.and.kend.lt.nc) then
cio      kanf=kend+1
cio      goto 1129
cio   endif
cvp
*     extimp=extimp/1000
*     write(iwr,*) 'timings for sum of all check0-calls '
*     write(iwr,*) 'vu   time ',vutimp,'  seconds'
*     write(iwr,*) 'cpu  time ',cputip,'  seconds'
*     write(iwr,*) 'ex   time ',extimp,'  seconds'
cvp
*     write(iwr,*) '  vu - timings for a, b, c ',vutwwg,vutlag,vuting
*     write(iwr,*) ' cpu - timings for a, b, c ',cpuwwg,cpulag,cpuing
*     write(iwr,*)
c ende der timer-routine
*     call clockv(vu2,cpu2)
*     call time(milli2)
*     extime = milli2 - milli1
*     extime = extime/1000
*     vutime = vu2  - vu1
*     cputim = cpu2 - cpu1
*     write(iwr,*) 'timings for case dk=0, all p-cases'
*     write(iwr,*) 'vu   time ',vutime,'  seconds'
*     write(iwr,*) 'cpu  time ',cputim,'  seconds'
*     write(iwr,*) 'ex   time ',extime,'  seconds'
cvp
*     if (npdim.le.nedim) then
*        write(iwr,*) 'timings for hmdk0 '
*        write(iwr,*) 'vu   time ',vux1,'  seconds'
*        write(iwr,*) 'cpu  time ',cpux1,'  seconds'
*        write(iwr,*) 'ex   time ',xillx1/1000,'  seconds'
*     endif
c
6     continue
cold  jg=jg+1
cold  mi(jg)=-1
cold  mj(jg)=0
cvp   write (6,229) mj
cold  write(nhuk)q,mi,mj
c ende der timer-routine for skiny
c --- behalten
cvp   vue =0.0d0
cvp   cpue=0.0d0
cvp   millie=0
cvp   call clockv(vue,cpue)
cvp   call time(millie)
cvp   extime = millie - millia
cvp   extime = extime/1000
cvp   vutime = vue  - vua
cvp   cputim = cpue - cpua
cdebugwrite(iwr,*)
cdebugwrite(iwr,*)
cdebugwrite(iwr,*) ' end loop of skiny '
cdebugwrite(iwr,*)
cdebugwrite(iwr,*)
cdebugwrite(iwr,*) 'timings for skiny '
cdebugwrite(iwr,*) 'vu   time ',vutime,'  seconds'
cdebugwrite(iwr,*) 'cpu  time ',cputim,'  seconds'
cdebugwrite(iwr,*) 'ex   time ',extime,'  seconds'
csut  write(iwr,*) 'vu=',vutime,' |  cpu=',cputim,
csut .           ' | exec=',extime,' s'
c
      if (ioint.eq.1) then
cdebug   write(iwr,*) 'rewind on file ',mtype
         call rewftn(mtype)
      end if
      return
      end
c--------1---------2---------3---------4---------5---------6---------7--
c
c subroutine zum aufbau der av startmatrix
c
c baut die av startmatrix mit moeglichst wenig skiny aufrufen auf
c
c laeuft falls alle startvektoren auf einmal
c mehrere startvektoren hintereinander --> muss noch gecheckt werden
c     b.e 02.08.1993
cbe
cbe  ideks in iideks veraendert um es vom rest des programms zu trennen
cbe
cbe  achtung:
cbe  ========
cbe  im vergleich zu dvdmua.f sind die files vecf1 und vecf2 genau
cbe  vertauscht!!!!!!!!!!! d.h. vecf1 ist output, vecf2 ist input
cbe  bei aufruf skin..
cbe
c-----------------------------------------------------------------------
c
      subroutine startm(twoe, nteint,
     +    pey, acoul, aexc, ndeke, core,
     +    vecf1, vecf2, mdi, hp5, mdi0,
     +    emat, nedim,
     +    iottr, iotnew, iotm,
     +    maindf,maxr,jconb,nomax2,
     +    nnl, nstanl, lunnl, lun1nl, lun2nl, noutnl,
     +    orthon, nglei,
     +    zwi1,zwi2,av,iwr,odebug,ncf)
c
c nnl    groesse der matrix
c nstanl anzahl der startvektoren
c lunnl  logische einheit die startvektoren enthaelt
c lun1nl logische einheit die orthonormierte vektoren enthaelt
c lun2nl logische einheit die a*v enthaelt
c orthon orthogonalisierungs kriterium
c nglei  wieviele gleichzeitig behandelt werden koennen
c av     enthaelt die matrix bi*a*bj mit bi,bj startvektoren
c        rueckgabe der information
c zwi1   zwischen puffer
c zwi2   zwischen puffer
c
      implicit real*8 (a-h,p-z)
c
c maximum number of roots=mx
c
      parameter (mx=256, mx1=(mx*mx+mx) / 2)
c
c
c !!!! change dimension of work array w(5,mx) in hqrii1 !!!!
c      fill the h-matrix in subr. mvmpy
c      start vectors are on unit lun
c
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox
      integer lun20, lun21, lun22, lunalt
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lun20, lun21, lun22, lunalt
c
      real*8 twoe
      integer nteint
      dimension twoe(nteint)
      real*8 pey,core
      real*8 acoul,aexc
      dimension pey(ndeke),acoul(ndeke),aexc(ndeke)
c
      real*8 vecf1, vecf2, hp5
      integer mdi, mdi0
      dimension vecf1(mdi),vecf2(mdi),hp5(mdi0)
c
      real*8 emat
      integer nedim
      dimension emat(nedim)
c
      integer iottr, iotnew, iotm
      dimension iottr(iotm), iotnew(iotm)
c
      integer maindf, maxr, jconb, nomax2
      dimension maindf(maxr,maxr),jconb(maxr,nomax2)
c
      logical odebug
c
*     logical*4 od12,ex12
*     character*10 acc12,seq12,dir12
      common /rwork2/ ndim,iglei
      dimension av(mx1),iideks(mx),zwi1(nnl),zwi2(nnl)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor,orthon,ortho
      common /csemi/ hstor,issk
cvp for integralvorsortierung nach foxy
      integer nrec31,ioint
      real*8 rtol,cptol,cptolcc,cptolm
      common /rio/ rtol,cptol,cptolm,cptolcc,nrec31,ioint
cvp
c NB
      common /cmpos/ mipos(256)
c NB!

      data big /1.23456d+20/
c---
      write(iwr,777)
777   format(/1x,'***  building the start matrix')
c--- umsetzen der nichtlokalen parameter in lokale parameter
      n=nnl
      ndim=nnl
c n oder ndim sind dimension of matrix
      nstart=nstanl
c wieviel startvektoren
      nout=noutnl
c output (usually 6)
      lun=lunnl
c logical unit number with initial eigenvectors
      lun1=lun1nl
c logical unit mit orthonomierten vektoren
      lun2=lun2nl
c logical unit mit a*v vektoren
c termination crit.
      ortho=orthon
c thresholds for debuggen
      trawf=0.3d0
      trahb=0.9d0
c--- input check
      if (odebug) then
      write(iwr,*) 'dimension ',n
      write(iwr,*) 'groesse des scratcharrays ',mdi
      write(iwr,*) 'anzahl der startvektoren ',nstart
      write(iwr,*) 'output auf ',nout
      write(iwr,*) 'startvektoren auf',lun
      write(iwr,*) 'startvektoren sichern auf',lunalt
      write(iwr,*) 'orthonormierte startvektoren auf',lun1
      write(iwr,*) 'a*bj vektoren auf',lun2
      write(iwr,*) 'gleichzeitig koennen ',nglei
      endif

c aufbau der indices der unteren dreiecksmatrix
      iideks(1)=1
      do 10 idoof=2,mx
        iideks(idoof) = iideks(idoof-1) + idoof
10    continue
      call rewftn(lun)
      call rewftn(lun1)
      call rewftn(lun2)
c
c beginn des aufbaus der startmatrix av
      j=0
c nanz        wie haeufig do 1000 durchlaufen werden muss
c nrest       anzahl der vektoren die noch behandelt werden muessen
c ifrueh      der wievielte durchlauf von do 1000 findet gerade statt
c             ifrueh - 1 anzahl der frueheren
c nfrueh      wieviele vektoren wurden schon behandelt
c nglei       wieviele koennen gleichzeitig behandelt werden
c iglei       wieviele werden in diesem durchlauf gleichzeitig behandelt
c jetzt       nr. im jetzigen durchlauf
      nanz = int(nstart/nglei) + 1
cdeb      write(iwr,*)'number of cycles ',nanz
      write(iwr,1) nstart,nanz
1     format(/1x,i3,'  start vectors are treated in ',i3,
     *' cycles ')
       if (nstart.gt.nko) then
          print *,' # ref confs must exceed start vectors and is ',nko
          call caserr('we need more ref. for multi dav')
       end if
      nrest = nstart
      nfrueh = 0
c
c     clear av matrix
c
      call vclr(av,1,iideks(nstart))
c
      do 1000 ifrueh=1,nanz
      if(odebug)
     + write(iwr,*) ifrueh,' -th cycle  ',nrest,'  left '
*        nres1 = abs(nrest - nglei)
        if(nglei.le.nrest) then
          iglei = nglei
        else
          iglei = nrest
        endif
        nrest = nrest - iglei
        if (odebug) then
         write(iwr,*) iglei,'  are treated at the same time '
         write(iwr,*) nrest,' left'
        endif
c einlesen der jetzt zu behandelnden startvektoren
c die stehen einzeln auf lun
        iein=0
        do 1011 irein = 1,iglei
          read (lun) zwi1
          write(lunalt) zwi1
          if(odebug) then
c         write(iwr,*) 'reading the ',irein,' -th start vector'
c         write(iwr,*)'zwi1(1-4)',(zwi1(iaus),iaus=1,4)
c NB
          do iaus=1,ncf
            write(iwr,*) 'zwi1(',mipos(iaus),')=',zwi1(mipos(iaus))
          end do
c NB
          write(iwr,*) 'start vector saved on ',lunalt
          endif
          do 1012 idoof = 1,n
            iein = iein +1
            vecf1(iein) = zwi1(idoof)
1012      continue
c         write(iwr,*)'vecf1(1-4)',(vecf1(iaus),iaus=1,4)
1011    continue
        if(odebug)
     +     write(iwr,*) iein,' elements restored'
c
c janf   : startadresse for den startvektor der jetzt behandelt wird
        janf = 1
        do 1100 jetzt=1,iglei
        if(odebug) then
         write(iwr,*) 'in do 1100 vector currently considered is',jetzt
         write(iwr,*)'vecf1(',janf,'-janf+4)',
     +                (vecf1(iaus),iaus=janf,janf+4)
        endif
c orthogonalierung des jetzt-ten bzgl
c der in frueheren lesezyklen behandelten
c ifrueh ist der zaehler der die lesezyklen zaehlt
          call rewftn(lun1)
          do 1200 idani=1,ifrueh - 1
cdeb            write(iwr,*) 'out of do 1200'
c einlesen von startvektoren, die in frueheren do 1000 durchlaeufen
c behandelt wurden
            iein = 0
            do     irein = 1,nglei
cdeb              write(iwr,*) 'reading the ',irein,' -th vector'
              read (lun1) zwi2
              do   idoof = 1,n
                iein = iein + 1
                vecf2(iein) = zwi2(idoof)
              enddo
            enddo
c           read (lun1) vecf2
c ianf   : startadresse for den startvektoren in vecf2
            ianf = 1
c da frueher nur ganze packete behandelt wurden kann man nglei nehmen
            do 1210 imarl = 1,nglei
cbe              write(iwr,*) 'vor vdotd'
c             call vdotd(n,vecf2(ianf),vecf1(janf),s)
              s = ddot(n,vecf2(ianf),1,vecf1(janf),1)
c             smax1= max(smax1,abs(s))
              s=-s
cbe              write(iwr,*) 'vor vsma1d'
              call vsma1d(n,s,vecf2(ianf),vecf1(janf))
cbe              write(iwr,*) 'nach vsma1d'
              ianf = ianf + n
 1210       continue
 1200     continue
c weglaufen der orthonormalisierung
c gehen wir mal davon aus das nur geeignete startvektoren uebergeben werd
c         if(smax1.gt.ortho) then
c           if(smax1.gt.smax2) to to 250
c           smax2=smax1
c           go to 60
c         endif
c orthogonalisierung bzgl aller im jetzigen lesezyklus
c kanf   : startadresse for die schon behandelten startvektoren auf vecf
          kanf = 1
          do 1300 imarl = 1,jetzt-1
c            write(iwr,*) 'aus do 1300 vor vdotd'
c      write(iwr,*)'vecf1(',kanf,'-kanf+4)',(vecf1(iaus),iaus=kanf,kanf+4)
c      write(iwr,*)'vecf1(',janf,'-janf+4)',(vecf1(iaus),iaus=janf,janf+4)
c           call vdotd(n,vecf1(kanf),vecf1(janf),s)
            s = ddot(n,vecf1(kanf),1,vecf1(janf),1)
c      write(iwr,*)'vecf1(',kanf,'-kanf+4)',(vecf1(iaus),iaus=kanf,kanf+4)
c      write(iwr,*)'vecf1(',janf,'-janf+4)',(vecf1(iaus),iaus=janf,janf+4)
c          write(iwr,*)'s ',s
cbe         write(iwr,*) 'vor vdotd'
c           smax1 = max(smax1,abs(s))
            s = -s
cbe            write(iwr,*) 'vor vsma1d'
            call vsma1d(n,s,vecf1(kanf),vecf1(janf))
c      write(iwr,*)'vecf1(',kanf,'-kanf+4)',(vecf1(iaus),iaus=kanf,kanf+4)
c      write(iwr,*)'vecf1(',janf,'-janf+4)',(vecf1(iaus),iaus=janf,janf+4)
            kanf = kanf + n
1300      continue
c normierung des jetztigen
c         call vdotd(n,vecf1(janf),vecf1(janf),s)
          s = ddot(n,vecf1(janf),1,vecf1(janf),1)
c         write(iwr,*) 'nach vdotd'
c      write(iwr,*)'vecf1(',janf,'-janf+4)',(vecf1(iaus),iaus=janf,janf+4)
c      write(iwr,*) 's ',s
          s = 1.0d0 / sqrt(s)
cbe              write(iwr,*) 'vor vsmd'
c      write(iwr,*) '1/s ',s
          call vsmd(n,s,vecf1(janf),vecf1(janf))
c      write(iwr,*)'vecf1(',janf,'-janf+4)',(vecf1(iaus),iaus=janf,janf+4)
          janf = janf + n
1100    continue
c alle iglei vektoren auf die vorherigen orthonormiert
c wegschreiben der neuen orthonormierten vektoren auf lun1
cbe        write(iwr,*) 'nach do 1100'
        iraus = 0
        do     irein = 1,iglei
          do   idoof = 1,n
            iraus = iraus + 1
            zwi1(idoof) = vecf1(iraus)
          enddo
          write (lun1) zwi1
        enddo
c nullsetzen von vecf1 nachdem auf vecf2 umgespeichert wurde
        nende = nglei*n
        do i0=1,nende
           vecf2(i0) = vecf1(i0)
        end do
cbe        write(iwr,*) 'vor zro1d'
        call zro1d (nende,vecf1)
c bildung von a*bi
c aufruf des skiny um alle jetzt im speicher befindlichen vektoren
c zu multiplizieren: vecf1 ist output; vecf2 ist input (d.h. vertauscht
c im vergleich zu dvdmua.f
c
c uebergabe parameter checken
c
c--- zunaechst aus testzwecken,
c    zunaechst umspeichern von vecf1 auf zwi1 einzelner aufruf des skiny
c    dann rueckspeichern auf vecf2
        ianf = 1
csut   do 4500 idoof = 1,iglei
cdebug         write(iwr,*) 'aufruf skiny for',idoof,' -ten startvektor'
csut     iend = ianf + n - 1
csut     izwi=0
cdebug          write(iwr,*) 'nach izwi:von,bis ',ianf,iend
csut     do 5000 isau = ianf,iend
csutut      izwi=izwi+1
csut        zwi1(izwi) = vecf1(isau)
c5000      continue
cdebug         write(iwr,*) 'testausdruck von zwi1; in skiny rein'
cdebug         write(iwr,*) 'alle groesser als ',trawf
cdebug         do 37 ilu = 1, n
cdebug           if(abs(zwi1(ilu)).gt.trawf) then
cdebug             write(iwr,*) ilu,zwi1(ilu)
cdebug           endif
cdebug37        continue
csut          do 40 ilu = 1,n
csut            zwi2(ilu) = 0.0
csut40        continue
csut
csut    do i=1,(nanz*ndim)
csut       vecf2(i) = 0.0d0
csut    end do
c       write(iwr,*) 'vor skiny'
c       write(iwr,*) 'vecf1',(vecf1(iii),iii=1,10)
c       write(iwr,*) 'vecf2',(vecf2(iii),iii=1,10)
c       print *,'issk,iswhm',issk,iswhm
        if (issk.lt.iswhm) then
cbe h*b direct
         call cputim(cpu0)
         call skiny(twoe,nteint,
     +             pey,acoul,aexc,ndeke,core,
     +             vecf1,vecf2,mdi,emat,
     +             nedim,iottr,iotnew,iotm,
     +             maindf,maxr,jconb,nomax2,ioint,odebug)
         call cputim(cpu1)
         cputot = cpu1-cpu0
cde         write(iwr,*) 'time for skiny : ',cputot,'mit'
cde         write(iwr,*) 'der groesse       ',(nanz*ndim)
        end if
c       write(iwr,*) 'nach skiny'
c       write(iwr,*) 'vecf1',(vecf1(iii),iii=1,10)
c       write(iwr,*) 'vecf2',(vecf2(iii),iii=1,10)
        if (issk.gt.0) then
cbe h*b semidirect
c             write(iwr,*) 'vor sham'
c             write(iwr,*) 'vecf2',(vecf2(iii),iii=1,10)
           call sham(vecf1, vecf2, mdi)
c             write(iwr,*) 'nach sham'
c             write(iwr,*) 'vecf1',(vecf1(iii),iii=1,10)
        end if
c             write(iwr,*) 'vor sdiag'
c             write(iwr,*) 'vecf1',(vecf1(iii),iii=1,10)
cbe h*b diagonale
        call sdiag(vecf1, vecf2, mdi, hp5, mdi0)
c             write(iwr,*) 'nach sdiag'
c             write(iwr,*) 'vecf1',(vecf1(iii),iii=1,10)
c              write(iwr,*) 'testausdruck von h*b; aus skiny raus'
c              write(iwr,*) 'alle groesser als ',trahb
c              do 38 ilu = 1,n
c                if(abs(zwi2(ilu)).gt.trahb) then
cdebug             write(iwr,*) ilu,zwi2(ilu)
cdebug           endif
cdebug38        continue
csut          izwi=0
csut
cdebug          write(iwr,*) 'nach vecf2:von,bis ',ianf,iend
csut          do 5010 isau = ianf,iend
csut            izwi=izwi+1
csut            vecf2(isau) = zwi2(izwi)
csut5010      continue
csut          ianf = ianf + n
c4500    continue
csut
csut
c wegschrieben der a*v vektoren auf lun2
        iraus = 0
        do      irein = 1,iglei
          do    idoof = 1,n
            iraus = iraus + 1
            zwi2(idoof) = vecf1(iraus)
          enddo
c          write(iwr,*) 'schreibe ',irein,' ten hv - vektor'
c          write(iwr,*) (zwi2(iii),iii=1,10)
          write(lun2) zwi2
        enddo
c
c bildung von bi*a*bj
c aufbau des unteren dreiecks der av matrix
c zunaechst die elemente aus den im hauptspeicher enthaltenen vektoren
c anzahl der vorher behandelten vektoren ist nfrueh
        ilauf = iideks(nfrueh+1)
c ianf   : startadresse for h*bi vektoren auf vecf2
c janf   : startadresse for bj vektoren auf vecf1
cbe do 2000 baut untere dreiecksmatrix auf:
cbe 1h1
cbe 2h1 2h2
cbe ....
        ianf = 1
        do 2000 imani = 1, iglei
          janf = 1
          do 2010 jmani = 1, imani
          if(odebug) then
            write(iwr,*)'h*b ',imani,' start bei ',ianf
            write(iwr,*)vecf2(ianf),vecf2(ianf+1),
     *      vecf2(ianf+2)
            write(iwr,*)'mit b ',jmani,' start bei ',janf
            write(iwr,*)vecf1(janf),vecf1(janf+1),
     *      vecf1(janf+2)
          endif
c auf vecf2 (also brray) stehen die bi files,
c auf vecf1 (also array) stehen die a*bi files
c die anordnung von brray und array scheint korrekt zu sein
c           call vdotd(n, vecf2(ianf), vecf1(janf), s)
            s = ddot(n,vecf2(ianf),1,vecf1(janf),1)
            if(odebug) then
             write(iwr,*)'imani,jmani,ianf,janf,ilauf,s= ',
     +                    imani,jmani,ianf,janf,ilauf,s
            endif
            av(ilauf) = s
            ilauf = ilauf + 1
            janf = janf + n
2010      continue
          ilauf = ilauf + nfrueh
          ianf = ianf + n
2000    continue
c
c mit allen vorher dagewesenen
c vecf2 enthaelt die bi von frueher, vecf1 die a*bi von jetzt
        call rewftn(lun1)
        do 3000 ialex = 1,ifrueh-1
c lesen der orthonormierten startvektoren von lun1
*          write(iwr,*) 'direkt anfang 3000 iein = ',iein
          iein = 0
          do     irein = 1,nglei
*            write(iwr,*) 'lesen des ',irein,' -ten vektors'
            read(lun1) zwi2
            do   idoof = 1,n
              iein = iein + 1
              vecf2(iein) = zwi2(idoof)
            enddo
*            write(iwr,*) 'nach lesen 3000 iein = ',iein
          enddo
c         read(lun1) vecf1
          ilauf = nfrueh
          ndoof = (ialex - 1) * nglei
c janf  : startadresse for a*bi auf vecf1 (jetzt)
c ianf  : startadresse for bj auf vecf2 (frueher)
          janf = 1
c loop ueber vektoren auf vecf1 (bi)
c **** the llop terminator was incorrect in original code
c **** should be iglei and NOT nglei
          do 3100 jrebec = 1,iglei
c positionierung auf letzten eintrag des letzten do 1000 durchlaufs
c plus eins (also auf den ersten)
            jlauf = iideks(ilauf) + 1 + ndoof
            ianf = 1
c loop ueber vektoren auf vecf2 (bj)
c **** the llop terminator was incorrect in original code
c **** should be nglei and NOT iglei
            do 3200 jdani = 1,nglei
c auf vecf1 stehen die h*vi der jetzigen
c auf vecf2 stehen die vj der frueheren
c             call vdotd(n, vecf1(janf), vecf2(ianf), s)
              s = ddot(n,vecf1(janf),1,vecf2(ianf),1)
              av(jlauf) = s
              if(odebug) then
               write(iwr,*) 'av(',jlauf,')=',s
              endif
              jlauf = jlauf + 1
              ianf = ianf + n
3200        continue
            ilauf = ilauf + 1
            janf = janf + n
3100      continue
*      write(iwr,*) 'rewind auf ',lun2
*      rewind lun2
*      do irein = 1, iglei
*         write(iwr,*) lun2,' gelesen vor ende do 3000 '
*         read(lun2)
*      enddo
3000    continue
c
c ende des loops 1000
c
      nfrueh = nfrueh + iglei
1000  continue
      return
      end
cvp konvertierung der table-information bzgl. der determinanten
cvp auf integer
      subroutine tabcon(inopen,kml,ot)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
cvp
cvp
cvp ntabl enthaelt die table-information bzgl. der determinanten
cvp als integer
      integer ot
      dimension ot(motmax)
cvp
      do 100 i=1,motmax
  100 ntabl(i)=ot(i)
cvp
cvp inopen = anzahl der offenen schalen der sk(i)
cvp kml = anzahl der saf pro konfiguration der sk(i)
cvp konvertierung der negativen paritaet von 255 auf -1
cvp  hier for dk=2 spezialisiert
      n=inopen*(inopen-1)*(inopen-2)*(inopen-3)/24
c     write(6,*) ' inopen=,n=,kml= ',inopen,n,kml
      do 200 i=1,n
        j=i*(kml+kml+1)-kml-kml
        if (ntabl(j).eq.255) ntabl(j)=-1
  200 continue
      return
      end
cvp
cvp
cvp lesen des kopfes der table sowie der c- und der b-matrix
      subroutine tbread(ifrk,ifr1,ntab)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
cvp
cvp uebergabe von informationen von ft31 und ft33 an skiny
      integer ijone,iwod,mconf,nplu,idummy
      common/rskin1/ mconf,ijone,iwod,nplu,idummy
c
      dimension ijone(nirmax),mconf(iswhm),nplu(iswhm)
cvp
      integer jerk,jbnk,jbun
      common /jany/ jerk(10),jbnk(10),jbun(10)
cvp
      real*8 af
      common/rtable/ af(49),ae(7829)
      common/linkmr/ic,i2,i3
      common/bufd/gout(510),nword
      common/blksi3/nsz
      dimension lout(510),ndet(49)
      equivalence (lout(1),gout(1)),(af(1),ndet(1))
cvp
      is=nplu(1)
      mult=is+1
cvp
cx    call dropen (ntab,nad)
c
      is=jerk(mult)
      il=jbnk(mult)
      ill=2
      call stopbk3
      call rdbak3(is)
      call stopbk3
      call icopy(49,lout,1,ndet,1)
      is=is+nsz
      ik=1
5551  if(il.ge.ill)then
      call rdbak3(is)
      is=is+nsz
      call stopbk3
      call dcopy(nword,gout,1,ae(ik),1)
      ik=ik+nword
      ill=ill+1
      go to 5551
      endif
      return
      end
cvp
cvp endtransformation for dk=1
csut aus trfdk1  < --- abspeicherung der hamiltornmatrix
      subroutine trf1h(f,ae,kz,ie,nz,kml,kmj)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
cvp
      dimension f(ksafm*ksafm),ae(7829)
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8    hstor
      common /csemi/ hstor,issk
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
cvp
cvp
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
      do 45 m=1,kml
        mz=kz+m
        do 45 ih=1,kmj
          tim=0.0d0
          do 46 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
46          tim=tim+f(ly)*ae(mx)
        jjh=nz+ih
        if (dabs(tim).gt.hstor) then
          icount=icount+1
          imm(icount) = mz
          jmm(icount) = jjh
          rham(icount)= (tim)
          if (icount.ge.ndims) then
            irech=irech+1
            icount=0
            write(nf88)imm,jmm,rham
            do iii=1,ndims
               imm(iii)=0
               jmm(iii)=0
               rham(iii)=0.0d0
            end do
          end if
        end if
c
c schleife ueber die wurzeln
c
csut      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut        ipoint=(iwurz-1)*ndim
c
csut        vecf1(ipoint+mz)=vecf1(ipoint+mz)+vecf2(ipoint+jjh)*tim
csut        vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+mz)*tim
c ende schleife uber die wurzeln
csut 5100 continue
45    continue
cvp
      return
      end
cvp
cvp endtransformation for dk=0
      subroutine trfdk0(vecf1,vecf2,mdi,
     +                  f,ae,kz,ie,nz,kml,ibob)
cvp
cvp parameter und globale common-bloecke sind im include-file
c
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
c
c common for multi-root version
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common/rwork2/ ndim,nanz
cvp
cvp
      dimension f(ksafm*ksafm),ae(7829)
cvp
      do 145 l=1,kml
        do 145 m=1,kml
          tim=0.0d0
          do 143 ih=1,kml
            ly=(m-1)*kml+ih
            mx=ie+(l-1)*kml+ih
143         tim=tim+f(ly)*ae(mx)
c
c schleife ueber die wurzeln
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
          iih =kz+(1-ibob)*l+ibob*m
          jjh=nz+(1-ibob)*m+ibob*l
          vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
          vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*tim
c ende schleife uber die wurzeln
 5100 continue
145       continue
cvp
      return
      end
cvp
cvp endtransformation for dk=1
      subroutine trfdk1(vecf1, vecf2, mdi,
     +                  f,ae,kz,ie,nz,kml,kmj)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
      real*8 vecf1,vecf2
      integer mdi
      dimension vecf1(mdi),vecf2(mdi)
cvp
      dimension f(ksafm*ksafm),ae(7829)
c common for multi-root version
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
      common/rwork2/ ndim,nanz
cvp
cvp
c beginn der endtransformation von dk=1 (for alle p-faelle gleich)
      do 45 m=1,kml
        mz=kz+m
        do 45 ih=1,kmj
          tim=0.0d0
          do 46 ig=1,kml
            ly=(ih-1)*kml+ig
            mx=ie+(m-1)*kml+ig
46          tim=tim+f(ly)*ae(mx)
        jjh=nz+ih
c
c schleife ueber die wurzeln
c
      do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
        ipoint=(iwurz-1)*ndim
c
        vecf1(ipoint+mz)=vecf1(ipoint+mz)+vecf2(ipoint+jjh)*tim
        vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+mz)*tim
c ende schleife uber die wurzeln
 5100 continue
45    continue
cvp
      return
      end
cvp
cvp endtransformation for dk=0
      subroutine trfdkh(f,ae,kz,ie,nz,kml,ibob)
cvp
cvp parameter und globale common-bloecke sind im include-file
      implicit real*8 (a-h,p-z)
c ncmax  : max. # von konf. einer sk
      parameter (ncmax  =    5000)
c icmax  : recordlaenge fuer sortierte integrale
      parameter (icmax  =  100000)
c motmax : dimension des ot-feldes fuer die table-eintraege
      parameter (motmax =  152000)
      parameter (mxroot =      50)
      parameter (nad    =    4000)
      parameter (nad1   =   -3999)
      parameter (iaeq   =    1501)
c maxref : max. # der referenzkonfigurationen
      parameter (maxref =     256)
c momax  : max. # von mo's
      parameter (momax  =     256)
c maxshl : max. # besetzter mo's pro konfiguration
      parameter (maxshl =      50)
c nopmax : max. # offener schalen
      parameter (nopmax =      10)
c iswhm  : max. # von superkategorien
      parameter (iswhm  =       5)
      parameter (ksafm  =      48)
      parameter (kdetm  =     126)
c nirmax : max. # irred. darst. (d_{2h})
      parameter (nirmax =       8)
csut : bufferlaenge fuer hamiltonmatrix
      parameter (ndims  =    2048)
c ideksm : dimension von ideks(i)= i ueber 2
      parameter (ideksm = momax*(momax+1)/2+1)
c nitmax : dimension des nit-feldes --> startadressen fuer integrale
      parameter (nitmax = nirmax*(nirmax+1)*(nirmax*(nirmax+1)+2)/8)
c jabmax : dim. produkt-tabelle der irred. darst.
      parameter (jabmax = nirmax*(nirmax+1)/2)
cvp
c common-block fuer vergleich in checkx
c iottr  : konfigurationen (nach chfeld in 'transponierter' form)
c niot   : startadresse einer sk auf iot (konfigurationen nicht
c          transponiert)
c nconf  : # (selektierte) konfigurationen pro sk
c imo    : # mo's
c ibinom : binomialkoeffizienten
c nirred gibt an, in welcher irred. darst. ein mo ist
c mopos gibt die nummer eines mo's innerhalb einer irred. darst. an
c ndub   : # doppelt besetzte mo's pro konf. einer sk
c nytl   : # aller besetzten mo's pro konf. einer sk
c nod    : # einfach besetzte mo's pro konf. einer sk
      integer niot,nconf,imo,ibinom,nirred,mopos
      integer nytl,ndub,nod
      common/rvergl/ niot(iswhm),nconf(iswhm),imo
     & ,ibinom(0:nopmax+1,0:nopmax+1)
     & ,nirred(momax),mopos(momax)
     & ,ndub(iswhm),nytl(iswhm)
     & ,nod(iswhm)
c
c ideks(i) : (i ueber 2)
c nit      : startadresse fuer integrale
c ncimo    : # ci-mo's pro irred. darst.
      integer ideks,nit,ncimo
      common/rintgr/ ideks(ideksm),nit(nitmax),ncimo(nirmax)
cvp
cvp labels aus konfigurationsvergleich und zwischengroessen
c sumint : fuer p=3 summe ueber dreiindexintegrale
c nwwmo  : fuer p=3 dreiindexintegrale (bzgl. gemeinsamer einfach
c          besetzter mo's; z.t. dummy fuer integraladressen
c intcb  : integraladresse
c intex  : integraladresse
c moafal : adresse fuer einelektronenintegral oder zwischengroesse
c mobfal : adresse fuer einelektronenintegral oder zwischengroesse
c idiff  : # unterschiede bzgl. testkonfiguration
c jdiff  : unterschiede bzgl. testkonfiguration
c jmerkd : # unterschiede in den doppelt besetzten mo's
c jmerko : # unterschiede in den einfach besetzten mo's
c npfal  : p-faelle bzw. nummern der konf. mit best. p-fall
c nrfal  : r-faelle bzw. irred. darst. der wechselwirkenden mo's
c          bei p=3
c nqlfal : ql-fall
c nqrfal : qr-fall
c npos   : zwischengroesse (i.w. zur r-fall bestimmung)
c jposo  : position von einfach besetzten mo's in der testkonf.
c jposc  : position von doppelt besetzten mo's in der testkonf.
c jcon   : testkonfiguration --> welche mo's doppelt bzw. einfach
c                                besetzt sind
c jcon1  : testkonfiguration --> zum vergleich mit einfach bes. mo's
c jcon2  : testkonfiguration --> zum vergleich mit doppelt bes. mo's
c itest  : testkonfiguration
      real*8 sumint
      integer nwwmo,intcb,intex,moafal,mobfal
     & ,ispiel,idiff,jdiff,jmerkd,jmerko
     & ,npfal,nrfal,nqlfal,nqrfal,npos
     & ,jposo,jposc,jcon,jcon1,jcon2,itest
      common/rlabel/ sumint(ncmax)
     & ,nwwmo(ncmax,nopmax-1),intcb(ncmax),intex(ncmax)
     & ,moafal(ncmax),mobfal(ncmax)
     & ,ispiel(ncmax),idiff(ncmax),jdiff(ncmax)
     & ,jmerkd(ncmax),jmerko(ncmax)
     & ,npfal(ncmax),nrfal(ncmax),nqlfal(ncmax),nqrfal(ncmax)
     & ,npos(ncmax,2)
     & ,jposo(momax),jposc(momax)
     & ,jcon(momax),jcon1(momax),jcon2(momax)
     & ,itest(maxshl)
cvp
cvp hilfsfelder fuer konfigurationsvergleich
cvp
c iqmo(i,j): i=hoeheres mo, j=niedrigeres mo, aus den mo-nummern
c folgt der q-fall
c imoq(q-fall,1): hoeheres mo
c imoq(q-fall,2): niedrigeres mo
c imo3q(q-fall,1): hoechstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,2): mittleres mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,3): tiefstes mo von \phi_r fuer dk=1,p=1
c imo3q(q-fall,1...3): aus q-fall die positionen ww mo's der
c                      groesse nach sortiert
c imo4q(q-fall,1...4): aus q-fall die positionen ww mo's der
c idk2r: zur bestimmung des r-falles bei dk=2 fuer
c        linksbuendige mo-positionen
c nra,nrb,nrc: codierung fuer adressen von ex-integralen
cvp                      groesse nach sortiert
c
      integer iqmo,imoq,imo3q,imo4q
      integer idk2r,nra,nrb,nrc
      common/rhelp/iqmo(nopmax,nopmax-1),imoq(nopmax*(nopmax-1)/2,2)
     +     ,imo3q(nopmax*(nopmax-1)*(nopmax-2)/6,3)
     +     ,imo4q(nopmax*(nopmax-1)*(nopmax-2)*(nopmax-3)/24,4)
     +     ,idk2r(4,3),nra(3),nrb(3),nrc(3)
cvp
cvp felder bei erzeugung der hamiltonmatrix inklusive der indices
cvp
c hp3    :  matrixelemente fuer p=3
c fakcb  :  faktor aus sga-matrix fuer cb-integral
c fakex  :  faktor aus sga-matrix fuer ex-integral
c sac    :  dreiindexintegrale ueber gemeinsame einfach besetzte
c           mo's fuer p=3
c qfeld  :  matrixelemente
c cbfeld :  integrale
c exfeld :  integrale
c vgath1 :  gather des vektors h*c
c vgath2 :  gather des vektors   c
c mifeld :  erster index von h(i,j)
c mjfeld :  zweiter index von h(i,j)
c iaddr  :  startadresse fuer mifeld bzgl. einer konfiguration
c ipfeld :  paritaet der permutation (bisher nur fuer dk=2)
c nrcas  :  r-fall, falls nrfal anders belegt
c icas   :  fuer dk=0 benoetigt:   icas = 0 --> qr >= ql
c                                  icas = 1 --> qr <  ql
c ibob2  :  kodierung der tabelle bzgl. r und z fuer erstes integral
c jbob2  :  kodierung der tabelle bzgl. r und z fuer zweites integral
c ntabl  :  enthaelt die table-information bzgl. der determinanten
c           als integer*4 (bisher nur fuer dk=2)
c q      :  hamiltonmatrix  fuer ft35
c mi     :  erster  index   fuer ft35
c mj     :  zweiter index   fuer ft35
cvp
      real*8 sac,hp3,fakcb,fakex
      real*8 qfeld,vgath1,vgath2
      common/rsghm/hp3(ncmax),fakcb(ncmax),fakex(ncmax),sac(nopmax)
     +            ,qfeld(2*ncmax),vgath1(ncmax),vgath2(ncmax) 
      real*8 cbfeld,exfeld,q
      common/rsghm2/cbfeld(ncmax),exfeld(ncmax),q(1000)
      integer mifeld, mjfeld, iaddr, ipfeld, nrcas, icas
      integer ibob2, jbob2, ntabl
      integer mi,mj
      common/rsghm3/mifeld(2*ncmax),mjfeld(2*ncmax),iaddr(ncmax)
     &             ,ipfeld(ncmax),nrcas(ncmax),icas(ncmax)
     &             ,ibob2(3,3),jbob2(3,3)
     &             ,ntabl(motmax)
     &             ,mi(1000),mj(1000)
cvp
cvp felder zur bestimmung der wechselwirkenden konfigurationen
cvp  mit hilfe von erzeugern und vernichtern
cvp
c    nko=anzahl der referenzkonfigurationen
c    maindf=anregung der referenzen untereinander
c    iotnew=konfigurationen bzgl. anregung aus den referenzen
c    iotnst=startadressen der sk fuer iotnew
c    iref=referenzen wie bei input, iref(0,*)=sk der mains
c    jconb=langes jcon-feld fuer alle mains
c
      integer nko, iotnst, iref
      common/rref/nko,iotnst(iswhm),
     +            iref(0:maxshl,maxref)
c
c    idref=anregung von \phi_r bzgl. aller mains
c    jcona=jcon-feld fuer korrektur von idref zur bestimmung von idiff
c    fuer recon bzw. refcon1
c
      integer idref,jcona
      common/rwkref/idref(maxref),jcona(momax*2)
c
c common for multi-root version
c  vecf1 und vecf2: vektoren, evtl. for mehrere wurzeln
c  ndim: dimension des saekularproblems
c  nanz: anzahl der gleichzeitig behandelten wurzeln
csut /speicher
      integer icount,irech,imm,jmm
      real*8 rham
      common /rhsc/ icount,irech,
     .              imm(ndims),jmm(ndims),rham(ndims)
c --- issk : superkategorien auf disk
      integer issk
      real*8 hstor
      common /csemi/ hstor,issk
c
      integer nston, mtype, nhuk, ltape, ideli, ntab, mtapev
      integer nf88, iput, mousa, ifox, lunfil
      common /ftap5/ nston, mtype, nhuk, ltape, ideli,
     +               ntab, mtapev, nf88, iput, mousa, ifox,
     +               lunfil(4)
c
cvp
      dimension f(ksafm*ksafm),ae(7829)
cvp
      do 145 l=1,kml
        do 145 m=1,kml
          tim=0.0d0
          do 143 ih=1,kml
            ly=(m-1)*kml+ih
            mx=ie+(l-1)*kml+ih
143         tim=tim+f(ly)*ae(mx)
c
c schleife ueber die wurzeln
c
csut  do 5100 iwurz=1,nanz
c   ipoint: startpunkt des vektors innerhalb vecf for die wurzel iwurz
csut    ipoint=(iwurz-1)*ndim
c
csut      iih =kz+(1-ibob)*l+ibob*m
csut      jjh=nz+(1-ibob)*m+ibob*la
          if (dabs(tim).gt.hstor) then
            icount = icount+1
            imm(icount)=kz+(1-ibob)*l+ibob*m
            jmm(icount)=nz+(1-ibob)*m+ibob*l
            rham(icount) = (tim)
            if (icount.ge.ndims) then
              irech=irech+1
              icount=0
              write(nf88)imm,jmm,rham
              do ii=1,ndims
               imm(ii)=0
               jmm(ii)=0
               rham(ii)=0.0d0
              end do
            end if
          end if
csut
csut      vecf1(ipoint+iih)=vecf1(ipoint+iih)+vecf2(ipoint+jjh)*tim
csut      vecf1(ipoint+jjh)=vecf1(ipoint+jjh)+vecf2(ipoint+iih)*tim
c ende schleife uber die wurzeln
csut5100 continue
145       continue
cvp
      return
      end

      subroutine vdssvd(n,s1,v1,v2,v)
      implicit real*8 (a-h,o-z)
      dimension v1(n), v(n), v2(n)
c
      nl=n
      sl=s1
      do i=1,nl
       slv2 = sl - v2(i)
       if(dabs(slv2).ge.0.05d0) then
        v(i)= v1(i)/slv2
       endif
      enddo
c
      do i=1,nl
       slv2 = sl - v2(i)
       if(dabs(slv2).lt.0.05d0) then
        v(i)= v1(i)/dsign(0.05d0,slv2)
       endif
      enddo
c
      return
      end
      subroutine vsma1d(n,s1,v1,v)
      implicit real*8 (a-h,o-z)
      dimension v1(*), v(*)
c
      nl=n
      s1l=s1
      do 10 i=1,nl
  10     v(i)=s1l*v1(i)+v(i)
      return
      end
      subroutine vsmd(n,s1,v1,v)
      implicit real*8 (a-h,o-z)
      dimension v1(*), v(*)
c
      nl=n
      s1l=s1
      do i=1,nl
        v(i)=s1l*v1(i)
      end do
      return
      end
cvp geaendertes writs zum wegschreiben der sortierten integrale
      subroutine writs1(lu,a,la,jfull)
      real*8 a(la)
      write(lu) la,jfull,a
      return
      end
      subroutine zro1d(n,v)
      implicit real*8 (a-h,o-z)
      dimension  v(*)
c
      real*8 rzero
c
      parameter (rzero = 0.0d0)
      nl=n
      nrest = mod(nl,6)
c
      do i=1,nrest
         v(i) = rzero
      end do
c
      do i=nrest+1,nl,6
         v(i  ) = rzero
         v(i+1) = rzero
         v(i+2) = rzero
         v(i+3) = rzero
         v(i+4) = rzero
         v(i+5) = rzero
      end do
      return
      end
      subroutine ver_newmrd5(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/newmrd5.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
c NB
      subroutine cepash(j,a,u0,s)
       implicit none
c This subroutine shifts the davidson  matrix 
       real*8 a(*),u0(*),sl,s
       integer j,k,q,p
c
c j - size of the Davidson matrix
c a(*) - triangular av matrix
c u0 - <ui|uj> : shift matrix 
c s - shift for the current cepa variant
c
       sl= s
       do k=1,j,1
          p=(k*k-k)/2
          do q=1,k
             a(p+q) = a(p+q) + sl*u0(p+q)
          end do
       end do
      return
      end
cc NB
      subroutine vsmace(nr,v)
c nr - # of reference configurations
c v  - shift vector sh*Iv*C0
      implicit real*8 (a-h,o-z)
      dimension v(*)
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
c     common /cmpos/ mipos(256)
c
      nrl=nsneak
      do i=1,nrl
        v(sneak(i))= 0.0d0
c       print *,'vs(',mipos(i),')=',v(mipos(i))
      end do
      return
      end
cc NB
      subroutine vdssce(n,nr,s1,sh,v1,v2,v)
c
c n - size of the vectors
c nr - # of reference configurations
c s1 - current aproximation of the cepa energy
c sh - shift for the appropriate cepa variant
c v1 - (h*c - e*c + shift*Id)
c v2 - diaganal of the hamiltonian
c v  - new perturbation vector
c vl - local vector used to substract the shift 
c from the reference part of H diagonal
c
      implicit real*8 (a-h,o-z)
      dimension v1(*), v(*), v2(*), vl(n)
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
c     common /cmpos/ mipos(256)
c
      nl=n
      sl=s1
      shl=sh
      nrl=nsneak
c
      do 40 i=1,nl
 40     vl(i)= 0.0d0
c
      do 30 i=1,nrl
 30     vl(sneak(i))=shl
c
      do 10 i=1,nl
         slv2 = sl - v2(i) - shl + vl(i)
         if(dabs(slv2).ge.0.05d0) then
          v(i) = v1(i)/slv2
         else
c         print *,'i=',i,' slv2=',slv2
          v(i) = v1(i) / dsign(0.05d0,slv2)
         endif
c     write(450,35) i,-slv2
 35   format(i6,d25.4)
 10   continue
      return
      end
cc NB
      real*8 function corrce(nele,cepavar)
c calculates the damping factor for EVP correction
       real*8 cepaf(3)
       integer cl,nl,nele,cepavar
c
c cepavar: 1 cepa0
c          2 acpf
c          3 aqcc
c nele - # of correlated electrons
c
       cl = cepavar
       nl = nele
       cepaf(1)=1.0d+0
       cepaf(2)=(nl-2)/dble(nl)
       cepaf(3)=(nl-3)*(nl-2)/((nl-1)*dble(nl))
       corrce = cepaf(cl)
      return
      end
cc NB
      subroutine ddotce(n,nconf,vx,vy,val)
       real*8 vx(n),vy(n),val
       integer n,nconf,j
c sneak() field for sneaky references
c nsneak: # of sneaky references
      integer sneak(1000),nsneak
      common /cepa_sn1/sneak,nsneak
c      common /cmpos/ mipos(256)
c
       val = 0.0d0
       do j=1,nsneak
          val = val + vx(sneak(j))*vy(sneak(j))
       end do
      return
      end
