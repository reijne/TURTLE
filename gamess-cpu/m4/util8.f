c 
c  $Author: jmht $
c  $Date: 2008-12-08 23:13:28 +0100 (Mon, 08 Dec 2008) $
c  $Locker:  $
c  $Revision: 5778 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util8.m,v $
c  $State: Exp $
c  
      subroutine indsyt(num,q,comp,s,vlen)
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
      dimension q(num),comp(num,8),s(*),vlen(8)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt,nt3,itable(8,8),
     1mult(8,8),irr(maxorb)
c
c     apply projection operator to  mo to decompose it into
c     symmetry types
c     s is overlap matrix
c     vlen holds overlap of original mo with its projections
c
      do 30 irep = 1 , nt
         vlen(irep) = 0.0d0
         do 20 j = 1 , num
            comp(j,irep) = 0.0d0
 20      continue
 30   continue
c
c
      do 60 j = 1 , num
         do 50 iop = 1 , nt
            jj = ibasis(iop,j)
            cjj = q(iabs(jj))*dsign(1.0d0,dble(jj))
            do 40 irep = 1 , nt
               comp(j,irep)=comp(j,irep)+cjj*dble(itable(irep,iop))
 40         continue
 50      continue
 60   continue
c
c
      ant = 1.0d0/dble(nt)
      do 80 irep = 1 , nt
         do 70 j = 1 , num
            comp(j,irep) = comp(j,irep)*ant
 70      continue
 80   continue
      do 90 irep = 1 , nt
         vlen(irep) = vsv(num,q,comp(1,irep),s)
 90   continue
      return
      end
      subroutine inddeg(e,nsa,i,nred)
      implicit real*8  (a-h,o-z)
      dimension e(*)
      data small /1.0d-4/
c
c     checks to see how many orbitals are in degenerate set
      nred = 1
      ip1 = i
 20   ip1 = ip1 + 1
      if (ip1.gt.nsa) return
      if (dabs(e(ip1-1)-e(ip1)).gt.small) return
      nred = nred + 1
      go to 20
      end
      subroutine vschmv(num,v1,v2,s)
      implicit real*8  (a-h,o-z)
      dimension v1(*),v2(*),s(*)
c     schmidt orthogonalise vector v2 to v1
c     do not use for a whole set of vectors as
c     n**4 algorithm will result
      dum = vsv(num,v1,v2,s)
      do 20 i = 1 , num
         v2(i) = v2(i) - dum*v1(i)
 20   continue
      return
      end
      subroutine indprj(q,num,comp,s,irep)
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
      dimension q(num),comp(num),s(*)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt,nt3,itable(8,8),
     1mult(8,8),irr(maxorb)
c
c     apply projection operator to  mo to project out all
c     except symmetry type irep
c     s is overlap matrix
c
      do 20 j = 1 , num
         comp(j) = 0.0d0
 20   continue
c
c
      do 40 iop = 1 , nt
         do 30 j = 1 , num
            jj = ibasis(iop,j)
            cjj = q(iabs(jj))*dsign(1.0d0,dble(jj))
            comp(j) = comp(j) + cjj*dble(itable(irep,iop))
 30      continue
 40   continue
c
c
      ant = 1.0d0/dble(nt)
      do 50 j = 1 , num
         comp(j) = comp(j)*ant
 50   continue
      call vrenrm(num,comp,s)
      return
      end
      subroutine indcha(q,num,nred,comp,s,char,ichar,iw)
      implicit real*8  (a-h,o-z)
      logical repeat
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
      dimension q(num,nred),comp(num,8),s(*),char(8),ichar(8)
      common/symmos/imos(8,maxorb),ibasis(8,maxorb),nt,nt3,itable(8,8),
     1mult(8,8),irr(maxorb)
c
c     works out character of a set of degenerate orbitals
c     apply projection operator to mo's to decompose into
c     symmetry types then take overlap with original mo's
c     s is overlap matrix
c     char accumulates overlaps of original mo's with their
c     projections
c
      repeat = .false.
 20   do 30 irep = 1 , nt
         char(irep) = 0.0d0
 30   continue
c
c
      ant = 1.0d0/dble(nt)
      do 130 mo = 1 , nred
c
c
         do 50 irep = 1 , nt
            do 40 j = 1 , num
               comp(j,irep) = 0.0d0
 40         continue
 50      continue
c
c
         if (repeat) then
           write (iw,*) ' analysis of mo ' , mo , ' in set of ' , nred
           write (iw,*) ' m.o. is '
           write (iw,'(5f20.15)') (q(j,mo),j=1,num)
         end if
         do 80 j = 1 , num
            do 70 iop = 1 , nt
               jj = ibasis(iop,j)
               cjj = q(iabs(jj),mo)*dsign(1.0d0,dble(jj))
               do 60 irep = 1 , nt
                  comp(j,irep) = comp(j,irep)
     +                           + cjj*dble(itable(irep,iop))
 60            continue
 70         continue
 80      continue
c
c
         do 100 irep = 1 , nt
            do 90 j = 1 , num
               comp(j,irep) = comp(j,irep)*ant
 90         continue
 100     continue
         if (repeat) then
            write (iw,*) ' decomposition into symmetry components'
            do 110 irep = 1 , nt
               write (iw,*)
               write (iw,*) ' irep ' , irep
               write (iw,*)
               write (iw,'(5f20.15)') (comp(j,irep),j=1,num)
 110        continue
         end if
         do 120 irep = 1 , nt
            char(irep) = char(irep) + vsv(num,q(1,mo),comp(1,irep),s)
 120     continue
         if (repeat) then
            write (iw,*)
            write (iw,*) ' cumulative character '
            write (iw,*)
            write (iw,*) (char(irep),irep=1,nt)
         end if
 130  continue
      do 140 irep = 1 , nt
         ichar(irep) = nint(char(irep))
 140  continue
      ich = 0
      do 150 i = 1 , nt
         ich = ich + ichar(i)
 150  continue
c
c     if get error here probably indicates that the set
c     of orbitals being analysed does not completly span
c     degenerate set - make numerical test in degnum sloppier?
c
      if (ich.ne.nred) then
        write (iw,*) nred
        write (iw,*) ' if the character below contains non-integer'
        write (iw,*) ' values it is likely that the boundary between'
        write (iw,*) ' two shells of orbitals has split a degenerate'
        write (iw,*) ' set - the wavefunction is probably not that'
        write (iw,*) ' of a valid state'
        write (iw,*)
        write (iw,*) (char(i),i=1,nt)
        write (iw,*)
        write (iw,*) (ichar(i),i=1,nt)
         if (.not.repeat) then
            repeat = .true.
            go to 20
         end if
         call caserr(' error in character determination ')
      end if
      return
      end
      subroutine vrenrm(num,v,s)
      implicit real*8  (a-h,o-z)
      dimension v(*),s(*)
      data small/1.0d-24/
c
c     renormalise eigenvector
c     overlap matrix in s
      do 20 j = 1 , num
         if (dabs(v(j)).lt.small) v(j) = 0.0d0
 20   continue
      dum1 = 0.0d0
      jk = 0
      do 40 j = 1 , num
         qji = v(j)
         if (qji.ne.0.0d0) then
            dum1 = dum1 + qji*qji*s(jk+j)
            qji = qji + qji
            do 30 k = 1 , j - 1
               dum1 = dum1 + s(jk+k)*qji*v(k)
 30         continue
         end if
         jk = jk + j
 40   continue
      dum1 = dsqrt(dum1)
      if (dum1.lt.1.0d-6) call caserr(
     +'trying to normalise null vector')
      anorm = 1.0d0/dum1
      do 50 j = 1 , num
         v(j) = v(j)*anorm
 50   continue
      return
      end
      subroutine indtab(igrp,itable,mult)
      implicit real*8  (a-h,o-z)
      integer itable(8,8),datac1,cscic2(2,2),d2c2v(4,4),datc2h(4,4),
     1datd2h(8,8)
      dimension mult(64),mmm(64)
      data mmm/1,2,3,4,5,6,7,8,
     1         2,1,4,3,6,5,8,7,
     2         3,4,1,2,7,8,5,6,
     3         4,3,2,1,8,7,6,5,
     4         5,6,7,8,1,2,3,4,
     5         6,5,8,7,2,1,4,3,
     6         7,8,5,6,3,4,1,2,
     7         8,7,6,5,4,3,2,1/
      data datac1/1/
      data cscic2/3*1,-1/
      data d2c2v/6*1,2*-1,1,-1,1,-1,1,2*-1,1/
      data datc2h/5*1,-1,1,-1,2*1,2*-1,1,2*-1,1/
      data datd2h/10*1,2*-1,2*1,2*-1,1,-1,1,-1,1,-1,1,-1,
     11,2*-1,2*1,2*-1,5*1,4*-1,2*1,4*-1,2*1,1,-1,1,2*-1,
     21,-1,2*1,2*-1,1,-1,2*1,-1/
      do 20 i = 1 , 64
         mult(i) = mmm(i)
 20   continue
      go to (30,40,40,40,70,70,100,130) , igrp
 30   itable(1,1) = datac1
      return
 40   do 60 i = 1 , 2
         do 50 j = 1 , 2
            itable(i,j) = cscic2(i,j)
 50      continue
 60   continue
      return
 70   do 90 i = 1 , 4
         do 80 j = 1 , 4
            itable(i,j) = d2c2v(i,j)
 80      continue
 90   continue
      return
 100  do 120 i = 1 , 4
         do 110 j = 1 , 4
            itable(i,j) = datc2h(i,j)
 110     continue
 120  continue
      return
 130  do 150 i = 1 , 8
         do 140 j = 1 , 8
            itable(i,j) = datd2h(i,j)
 140     continue
 150  continue
      return
      end
      function ksign(i)
      implicit real*8  (a-h,o-z)
      if (i.lt.0) then
         ksign = -1
      else if (i.eq.0) then
         ksign = 0
      else
         ksign = 1
      end if
      return
      end
      function vsv(num,v1,v2,s)
      implicit real*8  (a-h,o-z)
      dimension v1(*),v2(*),s(*)
c
c     overlap between vectors v1 and v2, s is overlap matrix
      dum = 0.0d0
      kl = 0
      do 30 k = 1 , num
         dum = dum + v1(k)*v2(k)*s(kl+k)
         do 20 l = 1 , k - 1
            dum = dum + (v1(k)*v2(l)+v1(l)*v2(k))*s(kl+l)
 20      continue
         kl = kl + k
 30   continue
      vsv = dum
      return
      end

      subroutine ver_util8(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/util8.m,v $
     +     "/
      data revision /"$Revision: 5778 $"/
      data date /"$Date: 2008-12-08 23:13:28 +0100 (Mon, 08 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
