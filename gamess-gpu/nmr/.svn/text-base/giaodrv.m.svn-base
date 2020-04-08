      subroutine giaodrv(q)
c
c     driving routine for NMR chemical shift calculations
c
      implicit REAL  (a-h,o-z)
INCLUDE(../m4/common/sizes)
      dimension q(*)
c
      common/maxlen/maxq
INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/vectrn)
c
INCLUDE(../m4/common/common)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/cndx40)
INCLUDE(../m4/common/cndx41)
c
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/statis)
      common/scfopt/maxcyc,mconv,nconv
c
      common/junke/maxt,ires,ipass,nteff,
     1     npass1,npass2,lentri,nbuck,mloww,mhi,ntri,iacc,iontrn
c
INCLUDE(../m4/common/restrj)
      logical lsave
INCLUDE(../m4/common/prnprn)
INCLUDE(../m4/common/timeperiods)
      character *8 fkd,blank,grhf,dipd
c     character *8 closed,oscf
      data fkd,blank /'fockder','        '/
      data grhf/'grhf'/
c     data closed,oscf/'closed','oscf'/
      data dipd /'dipder'/
c
c     evaluate integrals
c
      if (opass2) then
       if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +    iofrst.ne.iofsym) then
        write (iwr,130)
        opass2 = .false.
       endif
      endif

      call start_time_period(TP_2D_TOTAL)

      nopk = 1
      iofsym = 0
      isecvv = isect(8)
      itypvv = 8
      nconv = max(nconv,7)
      mmaxq = maxq
      if (irest5.ge.6) then
      else
         if (irest5.lt.1) then
            call start_time_period(TP_2D_AOINTS)
            call integ(q)
            call end_time_period(TP_2D_AOINTS)
            irest5 = 1
            call revise
         end if
         if (irest5.lt.2) then
            call start_time_period(TP_2D_SCF)
            call scfrun(q)
            call end_time_period(TP_2D_SCF)
            irest5 = 2
            irest = 5
            call revise
         end if
         if (irest5.lt.3) then
            call start_time_period(TP_2D_HFGRDN)
            ibase = igmem_alloc_all(mword)
            maxq = mword
            lword4 = mword
            call initgiaox
            call giaox(ibase)
            call end_time_period(TP_2D_HFGRDN)
            irest5 = 3
            call revise
            call gmem_free(ibase)
            maxq = mmaxq
         end if
      end if
c
      call end_time_period(TP_2D_TOTAL)
c
c     ----- reset core allocation
c
      return
 130  format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
      end
c
c-----------------------------------------------------------------------
c
      subroutine initgiaox
      implicit none
c
c     copy stuff from GAMESS-UK common blocks to HONDO common blocks
c
c     GAMESS-UK common blocks
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/restrl)
INCLUDE(../m4/common/restrz)
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/runlab)
c
c     HONDO common blocks (apart from the ones that have renamed 
c                          throughout the HONDO code)
c
      integer mxatom, mxbfn
      parameter (mxatom=MAXAT)
      parameter (mxbfn =MAXORB)
c
      character*8 wfntyp
      common/prop3c/wfntyp
c
      character*8 scftyph
      common/scfopt1/scftyph
c
      REAL tolqmt
      integer nqmt
      common/qmtopt/tolqmt,nqmt
c
      character*8 hndtyp,dsktyp,filtyp,pcktyp
      common/hndopt/hndtyp,dsktyp,filtyp,pcktyp
c
      integer idafh,navh,ioda,mxioda
      parameter (mxioda=1000)
      common/hdafile/idafh,navh,ioda(2,mxioda)
c
      integer nap,iap
      common/lcapid/nap,iap
c
      integer normfh,normph,itolh
      common/baspar/normfh,normph,itolh
c
      integer isingl,nbits
      common/machin1/isingl,nbits
c
      character*8 anam,bnam,bflab
      common/mollab/anam(mxatom),bnam(mxatom),bflab(mxbfn)
c
c     Functions
c
      integer ipg_nodeid, ipg_nnodes
      integer lenwrd
c
c     Local stuff
c
      integer i
      integer maxbit
      integer ndar
      character*8 zscf, zmp2
      data zscf     /'scf     '/
      data zmp2     /'mp2     '/
      data maxbit   / 64       /
c
      nap = ipg_nnodes()
      iap = ipg_nodeid()+1
c
      normfh = normf
      normph = normp
      itolh  = itol
c
      if (mp2) then
         wfntyp = zmp2
      else
         wfntyp = zscf
      endif
      scftyph = scftyp
      nqmt    = newbas0
      dsktyp  = "disk io "
c
      isingl = lenwrd()
      nbits  = maxbit/isingl
c
c---  copy atom names
c
      do i = 1, mxatom
         anam(i) = zaname(i)
      enddo
c
c---  open idafh: general purpose da-file
c
      idafh=110
_IF(cray)
      ndar=512
_ELSE
      ndar=1000
_ENDIF
      call daopen(idafh,ioda,navh,ndar)
      end
