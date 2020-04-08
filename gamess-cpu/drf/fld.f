cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hnd8 : fld fortran
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***********************************************************************
c*                                                                     *
c*  copyright i.b.m. 1989,1990                                         *
c*  all rights reserved.                                               *
c*                                                                     *
c*  permission is hereby granted to use but not to reproduce or        *
c*  distribute this material. the use is restricted to research        *
c*  purposes only. this material may not be reproduced for             *
c*  commercial purposes or included in a commercial product            *
c*  without the specific written permission of i.b.m. .                *
c*                                                                     *
c*  any publication based upon results obtained from this material     *
c*  must include the following attribution:                            *
c*                                                                     *
c*  this work is based, in part, on results from the motecc(tm)        *
c*  package.  motecc is a trademark of i.b.m. .                        *
c*                                                                     *
c*                                                                     *
c***********************************************************************
      subroutine fldinp
      implicit real*8 (a-h,o-z),integer  (i-n)
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
cafc
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     this file contains the parameter settings for piet van duijnen's
c     reactionfield programme together with the sizes needed for the
c     v. duijnen - gamess interface as developed by 
c     nederkoorn and v. lenthe
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     ==================================================================
c
c     parameter setting for reaction field programme van duijnen:
c
c=======================================================================
c
c     
      parameter (mxat = 10000)
      parameter (mxpol = 10000)
      parameter (mxpts = 10000)
      parameter (mxspe = 50)
      parameter (mxamb = 10)
      parameter (mxgrp = 500)
      parameter (mxnum = 1024)
      parameter (mxsh = 512)
      parameter (mxgs = 2048)
      parameter (mxnpts = 20)
      parameter (mxneq = 255)
      parameter (mcwrkmx = 10)
c
c ==================================================================
c
c     parameter setting new drf routines (drf/dimpar):
c
c===================================================================
c
      parameter (mxst = 5)
      parameter (mxgran = 10)
cmw      parameter (mxnpts = 50)
      parameter (mxgrpar = (mxgran*(mxgran+1))/2)
c
c     ==================================================================
c
c     parameter setting for v. duijnen - gamess interface:
c
c=======================================================================
c
      parameter (ncatom=500)
c
c  ncatom is maximal allowed classical atoms in zmat
c  total allowed atoms in zmat = maxnz + ncatom
c
c  maxrfa is maximal number of classical atoms in zmat + coordinates
c  total maximal atoms allowed is maxrfa + maxat
c
      parameter (maxrfa=mxpts)
c
c  maxamb is maximal number of ambigous atoms. it is set equal
c  to maximal z-cards possible in gamess input , i.e. 50
c
      parameter (maxamb=mxamb)
c
c=============================================
c for the pes-routines
c=============================================
c
      parameter (mxgr1=128, mxpes=100)
c
c=============================================
c for the drfexpc routines
c=============================================
c
      parameter (maxex=1000)
c
c=============================================
c for the drfsurf routine
c=============================================
c
      parameter (maxinv=10,maxnrc=10000)
c      parameter (maxnrc=10000)
c
c=============================================
c for the analrf routine
c=============================================
c
      common/iofile/ir,iwr,ip
cmw      common/fldpar/fieldx,fieldy,fieldz,ifield
      common/symtry/invt(48),nt,ntmax,ntwd,nosym
cmw      common/runopt/runtyp
c
      character*8 wfntyp
      common /wfnopt/ wfntyp
c
      character*8 runtyp
      common /runopt/ runtyp
c
      character*8 hndtyp
      common /hndopt/ hndtyp
c
      character*4 keyblk,keyfld
      character*8 ramir
c
      namelist /fld/ fldx,fldy,fldz
c
cmw      data keyblk /4h    /
      data keyblk /'    '/
cmw      data keyfld /4h fld/
      data keyfld /' fld'/
cmw      data ramir  /8hir-raman/
      data ramir /'ir-raman'/
c
      data fldx,fldy,fldz /0.0d+00,0.0d+00,0.0d+00/
      data zero   /0.0d+00/
c
c     ----- read namelist -$fld- -----
c
c     ifield=keyblk
c
c     ----- if -irramx- job , force no field -----
c
      if(runtyp.eq.ramir) go to 10
c
      rewind ir
      read(ir,fld,end=10)
c     ifield=keyfld
      write(iwr,9999) fldx,fldy,fldz
   10 continue
      fieldx=fldx
      fieldy=fldy
      fieldz=fldz
c
c     ----- print warning about symmetry -----
c
      if((fieldx.ne.zero.or.fieldy.ne.zero.or.fieldz.ne.zero).and.
     1   (nosym.eq.0)) write(iwr,9998)
      return
c
 9999 format(/,10x,29(1h-),/,10x,'dipole field strengths (a.u.)',
     1       /,10x,29(1h-),
     1       /,' fldx = ',f10.6,' fldy = ',f10.6,' fldz = ',f10.6)
 9998 format(/,1x,81(1h.),
     1       /,1x,'warning : an external electric field reduces the ',
     1            'symmetry of the electron density.',
     1       /,1x,'you may have to lower or turn off symmetry in ',
     1            '$cntrl- .'
     1       /,1x,81(1h.))
      end
