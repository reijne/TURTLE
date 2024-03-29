c---- memory counting routines -----------------------------------------
      subroutine memreq_jkder_dft(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,schwarz_ao,schwarz_cd)

      implicit none
c
c arguments
c
      character*8 zscftp
      integer iso
      REAL q
      integer nshels
      dimension iso(nshels,*),q(*)
      integer basi, basj, bask, basl
      REAL adens(*), bdens(*)
      REAL cfit(*), cfit2(*)
      REAL grad(3,*)
      REAL schwarz_ao, schwarz_cd
      dimension schwarz_ao(*), schwarz_cd(*)

INCLUDE(../m4/common/sizes)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
INCLUDE(../m4/common/specal)

        REAL dummy

INCLUDE(../m4/common/ghfblk)
INCLUDE(../m4/common/cigrad)

      integer nt
      data nt/1/

INCLUDE(../m4/common/timez)

      REAL dgout
      common/tgrad/dgout(9)

INCLUDE(common/dft_d2escr)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)  ! debug switches
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_dshlno)


INCLUDE(../m4/common/ijlab)

INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/parallel)

      logical odbg
      common/dbgdbg/odbg
c
c local variables
c      
      integer m1, m2, m3, m0
      dimension m0(48),m1(48),m2(48),m3(48)

      integer najkl, nabcl, nd2str, nactp, nshdm
      integer numcmo, n4
      integer id3, id4, id5, id6, id7, id8, id9
      integer i10, i20, i30, i40, i00
      integer iabd1, iabd2, iabd3
      integer ijshel
      integer i, j, issi, lnddm, lensec, m
      integer kt_max, inc1_max, ncmmm_max
      integer maxll, maxjj, maxkk, it
      logical ofpres, odpres, ogpres

      logical olab, olabc, olabcd
      integer ii, jj, kk, ll
      integer i0, j0, k0, l0
      integer l2,  m00
      integer klshel, kadi, id, jd, kd, ld, nd
      integer iblok
      integer nav, ifmp1, ifmp2, ifmp3
      integer iwor1, iwor2, iwor3
      integer ndgout, next
      integer kadij, kadijk, iceni
      logical omp2, ocifor, omp2w, ohf, ociopt, oskipp
      REAL tolij, tolijk, abmax
      REAL dtim, tim0
      REAL schwarz_lim, test
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf

      character*16 fnm,snm

      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      REAL gmax

      integer itmp(4), itmp1, itmp2, bas, idum
c     
c functions/stmnt fns
c
      integer memreq_pg_dgop
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, incr_memory2
      logical opg_root

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e_memory.m'/
      data snm/'memreq_jkder_dft'/
c
c     do i = 1,nshels
c        iso(i,1) = 1
c     enddo

      odbg = .false.
      if (schwarz_tol.ge.0) then
         schwarz_lim = 0.1d0**schwarz_tol
      else
         schwarz_lim = -1.0d0
      endif
      nschwz = 0
c
c this routine forces the load of the block data
c
       call dummy_intgx
c
c Determine which term to compute
c
      ncentr = 0
      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         if(basi .eq. bask)then
         write(6,*)'fitting coefficients',(cfit(i),i=1,nbasfn(basi))
         else
         write(6,*)'fitting coefficients 1',(cfit(i),i=1,nbasfn(basi))
         write(6,*)'fitting coefficients 2',(cfit2(i),i=1,nbasfn(basi))
         endif
      endif
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

      nat3 = natoms*3

      odpres = .false.
      ofpres = .false.
      ogpres = .false.

      itmp(1) = basi
      itmp(2) = basj
      itmp(3) = bask
      itmp(4) = basl
      do j=1,4
      if(itmp(j) .lt. 0) itmp(j) = itmp(1)
      do 20 i = 1 , nshell(itmp(j))
         bas = itmp(j)
         if (ktype(bas,i).eq.3) odpres = .true.
         if (ktype(bas,i).eq.4) ofpres = .true.
         if (ktype(bas,i).eq.5) ogpres = .true.
 20   continue
      enddo
c
      m = ntpdm + 9
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c     code for dgenrl preallocation 
c     buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do j=1,4
      bas = itmp(j)
      do 21 i = 1 , nshell(bas)
         kt_max=max(kt_max,ktype(bas,i))
 21   continue
      enddo
      inc1_max=(kt_max+1)**4
c
c     vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6

c      write(6,*)'dfg',odpres, ofpres,ogpres
c
c     ----- set pointers for partitioning of core -----
c
      i10 = null_memory()

      ndens = 1
         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
      end if
      if (ofokab) nfok = nat3*ndenin
c
c  core at i00 is for indexing
      i00 = incr_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = incr_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

c      write(6,*)'core i20',i20
c      call chkadr(q(i20))

      ic7 = incr_memory2(lnddm*9+inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

c     call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
c    +     ist,jst,kst,lst,ncentr,q)

c     next = ipg_dlbtask()
c
c     ----- ishell -----
c
c     ----- jshell -----
c
c     ----- kshell -----
c
c     ----- calculate q4 factor for this group of shells -----
c
c     ----- check for redundant combinations -----
c
c     ----- initialize dgout to zero -----
c
c     ----- form products of density matrix elements -----
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
c     ----- generate all 4 partial contributions to the gradient ----
c
c     ----- save gradient and restart data -----
c
c     call dfinal_dft(q,1,ii,basi,grad, natoms)
         idum = memreq_pg_dgop(3*natoms,'+')
c
      ic7 = ic7 + 1
      call decr_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call decr_memory2(iabd,'d',fnm,snm,'iabd')
      call decr_memory2(i00,'d',fnm,snm,'i00')

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c
      subroutine memreq_jkder_dft_genuse(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,
     &  ite3c_stored, nte3c_shl, ite2c_stored, nte2c_shl)

      implicit none
c
c     This routine is the derivative equivalent of jkint_dft_genuse.
c     The results of the Schwarz inequality are stored in tables
c     so that we can lookup whether a block of integrals was used in
c     the energy evaluation. 
c
c arguments
c
      character*8 zscftp
      integer iso
      REAL q
      integer nshels
      dimension iso(nshels,*),q(*)
      integer basi, basj, bask, basl
      REAL adens(*), bdens(*)
      REAL cfit(*), cfit2(*)
      REAL grad(3,*)
      integer nte3c_shl, nte2c_shl, ite3c_stored, ite2c_stored
      integer ite3c_shl, ite2c_shl
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)

INCLUDE(../m4/common/sizes)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
INCLUDE(../m4/common/specal)

      REAL dummy

INCLUDE(../m4/common/ghfblk)
INCLUDE(../m4/common/cigrad)

      integer nt
      data nt/1/

INCLUDE(../m4/common/timez)

      REAL dgout
      common/tgrad/dgout(9)

INCLUDE(common/dft_d2escr)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)  ! debug switches
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_dshlno)


INCLUDE(../m4/common/ijlab)

INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/parallel)

      logical odbg
      common/dbgdbg/odbg
c
c local variables
c      
      integer m1, m2, m3, m0
      dimension m0(48),m1(48),m2(48),m3(48)

      integer najkl, nabcl, nd2str, nactp, nshdm
      integer numcmo, n4
      integer id3, id4, id5, id6, id7, id8, id9
      integer i10, i20, i30, i40, i00
      integer iabd1, iabd2, iabd3
      integer ijshel
      integer i, j, issi, lnddm, lensec, m
      integer kt_max, inc1_max, ncmmm_max
      integer maxll, maxjj, maxkk, it
      logical ofpres, odpres, ogpres

      logical olab, olabc, olabcd
      integer ii, jj, kk, ll
      integer i0, j0, k0, l0
      integer l2,  m00
      integer klshel, kadi, id, jd, kd, ld, nd
      integer iblok
      integer nav, ifmp1, ifmp2, ifmp3
      integer iwor1, iwor2, iwor3
      integer ndgout, next
      integer kadij, kadijk, iceni
      logical omp2, ocifor, omp2w, ohf, ociopt, oskipp
      REAL tolij, tolijk, abmax
      REAL dtim, tim0
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf
      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      REAL gmax

      integer itmp(4), itmp1, itmp2, bas, idum

      character *9  fnm
      character *16 snm
c     
c functions/stmnt fns
c
      integer memreq_pg_dgop
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, incr_memory2
      integer iipsci
      logical opg_root, oipsci

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft_genuse'/
c
c     do i = 1,nshels
c        iso(i,1) = 1
c     enddo

      odbg = .false.
      nschwz = 0
c
c this routine forces the load of the block data
c
       call dummy_intgx
c
c Determine which term to compute
c
      ncentr = 0
      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         if(basi .eq. bask)then
         write(6,*)'fitting coefficients',(cfit(i),i=1,nbasfn(basi))
         else
         write(6,*)'fitting coefficients 1',(cfit(i),i=1,nbasfn(basi))
         write(6,*)'fitting coefficients 2',(cfit2(i),i=1,nbasfn(basi))
         endif
      endif
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

      nat3 = natoms*3

      odpres = .false.
      ofpres = .false.
      ogpres = .false.

      itmp(1) = basi
      itmp(2) = basj
      itmp(3) = bask
      itmp(4) = basl
      do j=1,4
      if(itmp(j) .lt. 0) itmp(j) = itmp(1)
      do 20 i = 1 , nshell(itmp(j))
         bas = itmp(j)
         if (ktype(bas,i).eq.3) odpres = .true.
         if (ktype(bas,i).eq.4) ofpres = .true.
         if (ktype(bas,i).eq.5) ogpres = .true.
 20   continue
      enddo
c
      m = ntpdm + 9
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c     code for dgenrl preallocation 
c     buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do j=1,4
      bas = itmp(j)
      do 21 i = 1 , nshell(bas)
         kt_max=max(kt_max,ktype(bas,i))
 21   continue
      enddo
      inc1_max=(kt_max+1)**4
c
c     vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6

c      write(6,*)'dfg',odpres, ofpres,ogpres
c
c     ----- set pointers for partitioning of core -----
c
      i10 = null_memory()

      ndens = 1

         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
      end if
      if (ofokab) nfok = nat3*ndenin
c
c  core at i00 is for indexing
      i00 = incr_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = incr_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

      ic7 = incr_memory2(lnddm*9 + inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                   'ic7')
      ic7 = ic7 - 1

c     call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
c    +     ist,jst,kst,lst,ncentr,q)
c
c     ----- ishell -----
c
c     ----- jshell -----
c
c     ----- kshell -----
c
c     ----- calculate q4 factor for this group of shells -----
c
c     ----- check for redundant combinations -----
c
c     ----- initialize dgout to zero -----
c
c     ----- form products of density matrix elements -----
c
c     ----mess about with the density matix by eliminating
c     zero elements
c
c     ----- generate all 4 partial contributions to the gradient ----
c
c     ----- save gradient and restart data -----
c
c     ----- end of *shell* loops -----
c
c     ----- allocate core memory for symde
c
c     ----- reset core memory from symde
c
c     call dfinal_dft(q,1,ii,basi,grad, natoms)
         idum = memreq_pg_dgop(3*natoms,'+')
c
c     ----- reset core memory from jkder -----
c
      ic7 = ic7 + 1
      call decr_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call decr_memory2(iabd,'d',fnm,snm,'iabd')
      call decr_memory2(i00,'d',fnm,snm,'i00')

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c---- routines that do the real work -----------------------------------
_EXTRACT(jkder_dft,opteron,xeon,em64t)
      subroutine jkder_dft(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,schwarz_ao,schwarz_cd)

      implicit none
c
c arguments
c
      character*8 zscftp
      integer iso
      REAL q
      integer nshels
      dimension iso(nshels,*),q(*)
      integer basi, basj, bask, basl
      REAL adens(*), bdens(*)
      REAL cfit(*), cfit2(*)
      REAL grad(3,*)
      REAL schwarz_ao, schwarz_cd
      dimension schwarz_ao(*), schwarz_cd(*)

INCLUDE(../m4/common/sizes)

cc      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w

c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
cINCLUDE(../m4/common/restar)
cINCLUDE(../m4/common/restri)
c
INCLUDE(../m4/common/specal)

        REAL dummy
_IF(notused)
c-      REAL gijkl
c-      integer mword,nlenx,kworx,kworxx
c-      common/blkin/gijkl(510),mword,nlenx,kworx,kworxx
c-      integer ijkl
c-      common/craypk/ijkl(1360)

c-      REAL gijkl1
c-      integer mwor1,nlen1,kwor1,kwor11
c-      common/blk1/gijkl1(510),mwor1,nlen1,kwor1,kwor11

c-      integer ijkl1
c-      common/sortpk/ijkl1(1360)

c-      REAL gijkl2
c-      integer mwor2,nlen2,kwor2,kwor22
c-      common/bufc/gijkl2(510),mwor2,nlen2,kwor2,kwor22

c-      integer ijkl2
c-      common/three/ijkl2(1360)

c-      REAL gijkl3
c-      integer mwor3,nlen3,kwor3,kwor33
c-      common/bufd/gijkl3(510),mwor3,nlen3,kwor3,kwor33

c-      integer ijkl3
c-      common/lsort/ijkl3(1360)

c-      REAL dipd, dipn, dipi
c-      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)

c-      REAL val, vall
c-      integer icnt, mxtr
c-      common/dbuf/val(312),vall(195),icnt,mxtr(4)

c-      integer iao, jao, kao, lao, iprt
c-      common/dlabs/iao(312),jao(312),kao(312),lao(312),iprt(312)
_ENDIF

INCLUDE(../m4/common/ghfblk)
INCLUDE(../m4/common/cigrad)

cINCLUDE(../m4/common/symtry)
      integer nt
      data nt/1/

INCLUDE(../m4/common/timez)

      REAL dgout
      common/tgrad/dgout(9)

INCLUDE(common/dft_d2escr)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)  ! debug switches
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_dshlno)


INCLUDE(../m4/common/ijlab)

INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/parallel)

      logical odbg
      common/dbgdbg/odbg
c
c local variables
c      
      integer m1, m2, m3, m0
      dimension m0(48),m1(48),m2(48),m3(48)

      integer najkl, nabcl, nd2str, nactp, nshdm
      integer numcmo, n4
      integer id3, id4, id5, id6, id7, id8, id9
      integer i10, i20, i30, i40, i00
      integer iabd1, iabd2, iabd3
      integer ijshel
      integer i, j, issi, lnddm, lensec, m
      integer kt_max, inc1_max, ncmmm_max
      integer maxll, maxjj, maxkk, it
      logical ofpres, odpres, ogpres

      logical olab, olabc, olabcd
      integer ii, jj, kk, ll
      integer i0, j0, k0, l0
      integer l2,  m00
      integer klshel, kadi, id, jd, kd, ld, nd
      integer iblok
      integer nav, ifmp1, ifmp2, ifmp3
      integer iwor1, iwor2, iwor3
      integer ndgout, next
      integer kadij, kadijk, iceni
      logical omp2, ocifor, omp2w, ohf, ociopt, oskipp
      REAL tolij, tolijk, abmax
      REAL dtim, tim0
      REAL schwarz_lim, test
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf

      character*9 fnm,snm

      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      REAL gmax

_IF(64bitpointers)
      integer*8 itmp1, itmp2
_ELSE
      integer itmp1, itmp2
_ENDIF
      integer itmp(4), bas
c     
c functions/stmnt fns
c
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, allocate_memory2
      logical opg_root
c-      integer indq

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft'/
c
c-      indq(m,n) = (m-1)*n
c
c-      ompir = runtyp.eq.zdipd .or. runtyp.eq.zinfra
c-      ompir = ompir .and. omp2
c
c     ----- check for grhf or gvb cases -----
c
c-      ouhf = zscftp.eq.zuhf
c-      orgvb = zscftp.eq.zgvb
c-      ogrhf = zscftp.eq.zgrhf
c-      oclos = zscftp.eq.zrhf
c-      omcscf = zscftp.eq.zmcscf
c
c     ----- check for casscf
c
c-      ocas = zscftp.eq.zcas
c-      if (ocas .or. mp3 .or. (omp2 .and. .not.ompir)) then
c-         call setsto(1360,0,ijkl)
c-      end if

c
c disable prefactor testing
c
c      if (nprint.ne.-5 .and. oprint(57)) then
c         write (iwr,6010)
c         call writel(q(iprefa),nshels)
c      end if

c      dummy iso
c@@ need to decide how big this array should
c   be, and only use it when centres are interchangeable
c

      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
      if (schwarz_tol.ge.0) then
         schwarz_lim = 0.1d0**schwarz_tol
      else
         schwarz_lim = -1.0d0
      endif
      nschwz = 0
c
c this routine forces the load of the block data
c
       call dummy_intgx
c
c Determine which term to compute
c
      ncentr = 0
      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         if(basi .eq. bask)then
         write(6,*)'fitting coefficients',(cfit(i),i=1,nbasfn(basi))
         else
         write(6,*)'fitting coefficients 1',(cfit(i),i=1,nbasfn(basi))
         write(6,*)'fitting coefficients 2',(cfit2(i),i=1,nbasfn(basi))
         endif
      endif

      if(print_sw(DEBUG_FORCES) .and. opg_root())then
        write(6,*)'Computing',ncentr,' centre deriv integrals'
      endif
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

c-      oeof = .false.

cc      if( icontij .eq. DENSITY ) then

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

cc      endif

_IF(notused)
c-      if (ompir) then
c-         ndenin = 3
c-         ntpdm = 4
c-         iflden = idaf
c-         call secget(isect(31),31,ibden)
c-         iwor1 = 0
c-         iwor2 = 0
c-         iwor3 = 0
c-         mwor1 = 0
c-         mwor2 = 0
c-         mwor3 = 0
c-         ifmp1 = 20
c-         ifmp2 = 21
c-         ifmp3 = 22
c-         ib1 = 1
c-         call search(ib1,ifmp1)
c-         call search(ib1,ifmp2)
c-         call search(ib1,ifmp3)
c-         call setsto(1360,0,ijkl1)
c-         call setsto(1360,0,ijkl2)
c-         call setsto(1360,0,ijkl3)
c-         call secget(isect(57),57,iblok)
c-         call rdedx(dipd,lds(isect(57)),iblok,ifild)
c-      end if
c-      if (ogrhf) then
c-         m = 0
c-         call secget(isect(53),m,iblok)
c-         call rdedx(nact,lds(isect(53)),iblok,idaf)
c-      end if
_ENDIF

      nat3 = natoms*3
c      nbsq = num*num
c      lenb = lensec(nx)

      odpres = .false.
      ofpres = .false.
      ogpres = .false.

      itmp(1) = basi
      itmp(2) = basj
      itmp(3) = bask
      itmp(4) = basl
      do j=1,4
      if(itmp(j) .lt. 0) itmp(j) = itmp(1)
      do 20 i = 1 , nshell(itmp(j))
         bas = itmp(j)
         if (ktype(bas,i).eq.3) odpres = .true.
         if (ktype(bas,i).eq.4) ofpres = .true.
         if (ktype(bas,i).eq.5) ogpres = .true.
 20   continue
      enddo
c

      m = ntpdm + 9
c     if (omp2w .or. ompir) m = m + 3
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c     code for dgenrl preallocation 
c     buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do j=1,4
      bas = itmp(j)
      do 21 i = 1 , nshell(bas)
         kt_max=max(kt_max,ktype(bas,i))
 21   continue
      enddo
      inc1_max=(kt_max+1)**4
c
c     vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6

c      write(6,*)'dfg',odpres, ofpres,ogpres
c
c     ----- set pointers for partitioning of core -----
c

c i10 used for nconf

      i10 = null_memory()
cc      i20 = igmem_alloc(l2)
cc      i30 = igmem_alloc(l2)

      ndens = 1
c     if (ouhf) then
c        ndens = 2
c     end if
c     if (orgvb) then
c        ndens = 3
c     end if
c     if (ogrhf) then
c        ndens = njk
c     end if
c     if (ocas) then
c        ndens = 1
c     end if
c     if (ofock .or. ompir) then
c        ndens = ndens + ndenin
c     end if

cc      if( icontij .eq. DENSITY ) then

         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

cc      else
cc         ida = 0
cc         i20 = 0
cc         i30 = 0
cc      endif

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
c-         if (ogrhf) nfok = njk*natoms*6
      end if
      if (ofokab) nfok = nat3*ndenin
c

c-      ifok = igmem_alloc(nx*nfok)
c-      ifok = ifok - 1

c  core at i00 is for indexing
      i00 = allocate_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = allocate_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

c      write(6,*)'core i20',i20
c      call chkadr(q(i20))

      ic7 = allocate_memory2(lnddm*9+inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)
_IF(notused)
c-      if(omcscf)then
c-
c-         i40 = igmem_alloc(l2)
c-
c-ccc         id3 = i40 + nx
c-         numcmo = num*ncoorb
c-         nactp = ncact*(ncact+1)/2
c-         nd2str = indq(nactp+1,16)
c-         if (odpres) nd2str = indq(nactp+1,36)
c-         if (ofpres) nd2str = indq(nactp+1,100)
c-         if (ogpres) nd2str = indq(nactp+1,225)
c-         id3 = igmem_alloc(numcmo)
c-         i30 = i20
c-ccc         id4 = id3 + numcmo
c-         id4 = id3 + numcmo
c-c     nd2mo=ind(nactp,nactp)
c-         nd2mo = nactp*(nactp+1)/2
c-c     write(6,*)' length of tpdm',nd2mo
c-
c-         id4 = igmem_alloc(nd2str)
c-ccc         id5 = id4 + nd2str
c-         id5 = igmem_alloc(nd2mo)
c-ccc         id6 = id5 + nd2mo
c-c
c-         najkl = indq(ncact*4+1,nactp)
c-         if (odpres) najkl = indq(ncact*6+1,nactp)
c-         if (ofpres) najkl = indq(ncact*10+1,nactp)
c-         if (ogpres) najkl = indq(ncact*15+1,nactp)
c-c
c-         nabcl = indq(16+1,4*ncact)
c-         if (odpres) nabcl = indq(36+1,6*ncact)
c-         if (ofpres) nabcl = indq(100+1,10*ncact)
c-         if (ogpres) nabcl = indq(225+1,10*ncact)
c-c
c-         id6 = igmem_alloc(nabcl)
c-cc         id7 = id6 + nabcl
c-         id7 = igmem_alloc(najkl)
c-cc         id8 = id7 + najkl
c-         id8 = igmem_alloc(nactp)
c-cc         id9 = id8 + nactp
c-c
c-         nshdm = max(ncact,4)
c-         if (odpres) nshdm = max(ncact,6)
c-         if (ofpres) nshdm = max(ncact,10)
c-         if (ogpres) nshdm = max(ncact,15)
c-c
c-         id9 = igmem_alloc(nshdm*nshdm)
c-cc
c-cc         iabd = id9 + nshdm*nshdm
c-c
c-c   iabd should not be used for mcscf calculations - tpdm in id5
c-c
c-cc         i00 = iabd + lnddm*ntpdm
c-cc         ic7 = i00 + lnddm*(3/nav+1) - 1
c-      endif

c     tim0 = cpulft(1)
c
c-      if (ocas) then
c-         call dbutci(ist,jst,kst,lst)
c-      else if (omcscf) then
c-         call ddebut(zscftp,q(id3),q(i30),q(i10),q(i40),q(id5),
c-     +               ist,jst,kst,lst,q)
c-      else
_ENDIF

      call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
     +     ist,jst,kst,lst,ncentr,q)

_IF(notused)
c-      end if
c-      if (mp3 .or. (omp2 .and. .not.ompir)) call search(iblk2d,ifil2d)

c-      mword = 0
c-      iword = 0

c-      kworx = 999
c-      kwor1 = 999
c-      kwor2 = 999
c-      kwor3 = 999
c
c check 
c??      kloc(nshels+1) = num + 1

c
c-      icnt = 0
c-      ib1 = 1
c-      if (omp2w) then
c-         call search(ib1,mpstrm(1))
c-         if (odebug(30)) write (iwr,6020)
c-      end if
_ENDIF

      next = ipg_dlbtask()

      if(print_sw(DEBUG_PARALLEL))then
         write(6,*)'Task info on entry: ',ipg_nodeid(),next,icount_dlb
      endif

      if (ist.le.nshels) then
c
c     ----- ishell -----
c
         do 140 ii = ist , nshell(basi)

c-            kadi = kad(ii)
            ijshel = ii*(ii-1)/2

c
            if(ncentr .eq. 3 .or. ncentr .eq. 4)then
            do 40 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 140
               m0(it) = id
 40         continue
            endif

            iceni = katom(basi,ii)

c-            if (omcscf) call mcajkl(q(id3),q(id5),q(id7),q(id8),q(id9),
c-     +                              ii,nactp)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst
               maxjj = ii
            endif

            do 130 jj = j0 , maxjj
c-               kadij = kadi + kad(jj)
               jst = 1
               itrij = ijshel+jj

c-               tolij = dlntol + q(iprefa+ijshel+jj-1)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
c                  - DISABLED -
c-               if (tolij.gt.-3.401d0) then
                 if(.true.)then
c
c apply i/j tests only when i and j are AO basis fns
c
                  if (ncentr .eq. 3  .or. ncentr .eq. 4)then
                  do 60 it = 1 , nt
                     id = m0(it)
                     jd = iso(jj,it)
                     if (jd.gt.ii) go to 130
                     if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (id.eq.ii .and. jd.gt.jj) go to 130
                     m1(it) = id
                     m2(it) = jd
 60               continue
                  endif

                  if(basj.lt.0)then
                     olab = .true.
                  else
                     olab = katom(basj,jj).eq.iceni
                  endif

c
c     store information about the pair (ij)
c
                  call dshell_dft(1,ii,jj,kk,ll,
     +                 basi, basj, bask, basl,ncentr)

                  call dprim_dft

                  if (nij.ne.0) then

c-                     if (omcscf) call mcabkl(q(id3),q(id7),q(id4),ii,jj,
c-     +                   nactp)

                     icount_dlb = icount_dlb + 1
                     if(icount_dlb . eq. next) then

                        if(print_sw(DEBUG_PARALLEL))then
                         write(6,*)'Task ',next,' on node',ipg_nodeid()
                        endif

c
c     ----- kshell -----
c
                     if(ncentr .eq. 4 .or. 
     &                 (ncentr .eq. 2 .and. basi .eq. bask))then
                        k0 = kst
                        maxkk = ii
                     else
                        k0 = kst
                        maxkk = nshell(bask)
                     endif

                     if(odbg)write(6,*)'kk loop',k0, maxkk

                     do 120 kk = k0 , maxkk

c-                        kadijk = kadij + kad(kk)
                        kst = 1
                        klshel = kk*(kk-1)/2

                        if( ncentr .eq. 4) then
                        do 80 it = 1 , nt
                           kd = iso(kk,it)
                           if (kd.gt.ii) go to 120
                           m3(it) = kd
 80                     continue
                        endif

                        olabc = olab .and. katom(bask,kk).eq.iceni

c-                        if (omcscf)
c-     +                      call mcabcl(q(id3),q(id4),q(id6),q(id8),
c-     +                      q(id9),ii,jj,kk,nactp)

                        if(basl .lt. 0)then
                           l0 = 1
                           maxll = 1
                        else
                           l0 = lst
                           maxll = kk
                           if (kk.eq.ii) maxll = jj
                        endif

                        do 110 ll = l0 , maxll

                           lst = 1
                           if (ncentr.eq.3.and.schwarz_tol.ge.0) then
                              test = schwarz_ao(itrij) * schwarz_cd(kk)
                              oskipp = test.lt.schwarz_lim
                              if(oskipp) then
c                                mink = kmin(bask,kk)
c                                maxk = kmax(bask,kk)
c                                imc=imc+maxk-mink+1
                                 nschwz = nschwz + 1
                                 go to 110
                              endif
                              
                           endif

c-                           if (kadijk+kad(ll).lt.0) then
c-                              tolijk = tolij + q(iprefa+klshel+ll-1)
c-                              if (tolijk.gt.0.0d0) then

                           if(basl.lt.0)then
                              olabcd = olabc 
                           else
                              olabcd = olabc .and. 
     &                             katom(basl,ll).eq.iceni
                           endif

                           if(odbg)write(6,*)'olab',ii,jj,kk,ll,olabcd

                                 if (.not.(olabcd)) then

                                 if( ncentr .eq. 4 ) then

                                    n4 = 0
                                    do 100 it = 1 , nt
                                       ld = iso(ll,it)
                                       if (ld.gt.ii) go to 110
                                       kd = m3(it)
                                       if (kd.lt.ld) then
                                         nd = kd ! swap k,l
                                         kd = ld
                                         ld = nd
                                       end if
                                       id = m1(it)
                                       jd = m2(it)
                                       if (id.eq.ii .or. kd.eq.ii) then
                                         if (kd.ge.id) then
                                         if (kd.ne.id .or. ld.gt.jd)
     +                                      then
                                         nd = id  ! swap i,k
                                         id = kd
                                         kd = nd
                                         nd = jd  ! swap j,l
                                         jd = ld
                                         ld = nd
                                         end if
                                         end if
                                         if (jd.ge.jj) then
                                         if (jd.gt.jj) go to 110
                                         if (kd.ge.kk) then
                                         if (kd.gt.kk) go to 110
                                         if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 110
                                         n4 = n4 + 1
                                         end if
                                         end if
                                         end if
                                       end if
 100                                continue
c
c     ----- calculate q4 factor for this group of shells -----
c
                                    q4 = dfloat(nt)/dfloat(n4)

                                 elseif (ncentr .eq. 3) then

c rather empirical -it seems triangulation effects
c are already corrected for
c a factor of 2 is applied in dabab
                                    q4 = 1.0d0

                                 elseif (ncentr .eq. 2) then
                                    q4 = 1.0d0
                                 endif
c
c     ----- check for redundant combinations -----
c
                              call redund_dft(ii,jj,kk,ll,
     +                             basi, basj, bask, basl, 
     +                             iwr)


                              if (npass.eq.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk
                                 else if(ncentr .eq. 2)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,kk
                                 endif

                              endif
                              
                              if (npass.ne.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,q4
                                 else if(ncentr .eq. 2)then
                                    if(ii.eq.14.and.kk.eq.8)then
c                                       odbg=.true.
                                    else
c                                       odbg = .false.
                                    endif
                                    if(odbg)
     &               write(6,*)'Shell block',ii,kk, q4
                                 endif
c
c     ----- initialize dgout to zero -----
c
                          call vclr(dgout,1,ndgout)

                          call dshell_dft(2,ii,jj,kk,ll,
     +                         basi, basj, bask, basl,ncentr)
c
c     ----- form products of density matrix elements -----
c
                          call vclr(q(iabd+1),1,lendd*ntpdm)

_IF(notused)
c-                                 if (omcscf)
c-     +                              call mcabcd(q(id3),q(id6),
c-     +                              q(iabd+1),q(id9),ii,jj,kk,ll,
c-     +                              q4)
c-                                 if (mp3 .or.
c-     +                              (omp2 .and. .not.ompir))
c-     +                              call mcdab(q(iabd+1),ii,jj,kk,
c-     +                              ll,q4)
c-                                 if (ompir) then
c-                                   iabd1 = iabd + lendd + 1
c-                                   iabd2 = iabd1 + lendd
c-                                   iabd3 = iabd2 + lendd
c-                                   call dpdab1(q(iabd1),ii,jj,kk,
c-     +                                ll,q4,ifmp1,iwor1)
c-                                   call dpdab2(q(iabd2),ii,jj,kk,
c-     +                                ll,q4,ifmp2,iwor2)
c-                                   call dpdab3(q(iabd3),ii,jj,kk,
c-     +                                ll,q4,ifmp3,iwor3)
c-                                   call tpdder(q(iabd1),lendd,
c-     +                                q(i20),q(i30),nx,ndenin,ii,
c-     +                                jj,kk,ll,q4)
c-                                 end if
c-                                 if (.not.orgvb .and.
c-     +                              .not.ocas .and. .not.omcscf)
_ENDIF

                          call dabab_dft(ii,jj,kk,ll,
     &                         basi, basj, bask, basl,
     +                         ncentr,q4,
     +                         zscftp,adens,bdens,cfit,cfit2,
     +                         q(iabd+1))

c                          call chkadr2(q(i20),itmp1)
c                          call chkadr2(q(iabd+1),itmp2)
c                          write(6,*)'after dabab',itmp1,itmp2,
c     &                         q(i20),q(iabd+1)


_IF(notused)
c-                                 if (orgvb)
c-     +                              call dabg(ii,jj,kk,ll,l1,norb,
c-     +                              q4,q(i20),q(i30),q(i10)
c-     +                              ,onocor,onopen,q(iabd+1))
c-                                 if (ocas)
c-     +                              call dabci(ii,jj,kk,ll,q4,
c-     +                              oeof,q(iabd+1))
c-                                 if (omcscf)
c-     +                              call dabmc(ii,jj,kk,ll,q4,
c-     +                              q(i30),q(i40),q(iabd+1))
_ENDIF

c
c     ----mess about with the density matix by eliminating
c     zero elements
c

                                 call delim_dft(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)

_IF(notused)
c-                                 if (ompir) then
c-                                   call delim2(q(iabd1),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd2),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd3),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                 end if
_ENDIF
                                 
                                 if (ijkld.ne.0) then
                                    call dgenrl_dft(q(1),q(i00),q(i00),
     +                                   abmax,gmax)

c-                                   if (ofock)
c-     +                                call fockd2(q,q(i00))
c-                                   if (ofokab) then
c-                                     call fokabd(q,q(i00))
c-                                   end if
c
c     ----- generate all 4 partial contributions to the gradient ----
c
                                   call formeg_dft
                                 end if              ! ijkld.ne.o
                              end if                 ! olabcd
                                 end if              ! npass .ne. 0

c-                              end if                 ! (tolijk.gt.0.0d0)
c-                           end if                    ! (kadijk+kad(ll).lt.0)

 110                    continue                     ! ll   loop
 120                 continue
                     next = ipg_dlbtask()
                     endif
                  end if
               end if
 130        continue

c
c     ----- save gradient and restart data -----
c                ==== disabled =====
c            call dfinal_dft(q,0,ii,basi,grad,nat)
c            if (tim.ge.timlim) go to 150

 140     continue
         call pg_dlbpush

      end if
c
c     ----- end of *shell* loops -----
c
_IF(notused)
c-      if (omp2w) then
c-         if (icnt.ne.0) then
c-            call pack(vall,8,iao,1560)
c-            call wrt3s(val,511,mpstrm(1))
c-         end if
c-         icnt = 0
c-         m00 = 0
c-         call put(val,m00,mpstrm(1))
c-         call shut1(mpstrm(1))
c-         if (odebug(30)) write (iwr,6030)
c-      end if
_ENDIF

      if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
       call caserr('no symmetry in jkder_dft')
ccc       isymd = igmem_alloc(nw196(6))
c
ccc       call symde(q(isymd),natoms)

c
c     ----- reset core memory from symde
c
ccc      call gmem_free(isymd)
c
      endif
      call dfinal_dft(q,1,ii,basi,grad, natoms)

c150  continue
      call timit(0)
c     dtim = tim - tim0
c
c     ----- reset core memory from jkder -----
c

c-      if(omcscf)then
c-         call gmem_free(id9)
c-         call gmem_free(id8)
c-         call gmem_free(id7)
c-         call gmem_free(id6)
c-         call gmem_free(id5)
c-         call gmem_free(id4)
c-         call gmem_free(id3)
c-         call gmem_free(i40)
c-      endif

cc
cc revert
cc      call gmem_free(ic1)

      ic7 = ic7 + 1
      call free_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call free_memory2(iabd,'d',fnm,snm,'iabd')
      call free_memory2(i00,'d',fnm,snm,'i00')
c-      ifok = ifok + 1
c-      call gmem_free(ifok)

cc      if( icontij .eq. DENSITY ) then
        ida = ida + 1
cc      endif

cc      call gmem_free(i40)
cc      call gmem_free(i30)
cc      call gmem_free(i20)

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
c
_ENDEXTRACT
      subroutine jkder_dft_genuse(zscftp,q,iso,nshels,
     &  basi, basj, bask, basl, 
     & 	adens,bdens,
     &	cfit,cfit2,
     &	grad,
     &  ite3c_stored, nte3c_shl, ite2c_stored, nte2c_shl)

      implicit none
c
c     This routine is the derivative equivalent of jkint_dft_genuse.
c     The results of the Schwarz inequality are stored in tables
c     so that we can lookup whether a block of integrals was used in
c     the energy evaluation. 
c
c arguments
c
      character*8 zscftp
      integer iso
      REAL q
      integer nshels
      dimension iso(nshels,*),q(*)
      integer basi, basj, bask, basl
      REAL adens(*), bdens(*)
      REAL cfit(*), cfit2(*)
      REAL grad(3,*)
      integer nte3c_shl, nte2c_shl, ite3c_stored, ite2c_stored
      integer ite3c_shl, ite2c_shl
      dimension ite3c_stored(nte3c_shl), ite2c_stored(nte2c_shl)

INCLUDE(../m4/common/sizes)

cc      common/restrl/ociopt,ocifor,omp2,ohf(7),omp2w

c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
cINCLUDE(../m4/common/restar)
cINCLUDE(../m4/common/restri)
c
INCLUDE(../m4/common/specal)

        REAL dummy
c-      REAL gijkl
c-      integer mword,nlenx,kworx,kworxx
c-      common/blkin/gijkl(510),mword,nlenx,kworx,kworxx
c-      integer ijkl
c-      common/craypk/ijkl(1360)

c-      REAL gijkl1
c-      integer mwor1,nlen1,kwor1,kwor11
c-      common/blk1/gijkl1(510),mwor1,nlen1,kwor1,kwor11

c-      integer ijkl1
c-      common/sortpk/ijkl1(1360)

c-      REAL gijkl2
c-      integer mwor2,nlen2,kwor2,kwor22
c-      common/bufc/gijkl2(510),mwor2,nlen2,kwor2,kwor22

c-      integer ijkl2
c-      common/three/ijkl2(1360)

c-      REAL gijkl3
c-      integer mwor3,nlen3,kwor3,kwor33
c-      common/bufd/gijkl3(510),mwor3,nlen3,kwor3,kwor33

c-      integer ijkl3
c-      common/lsort/ijkl3(1360)

c-      REAL dipd, dipn, dipi
c-      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)

c-      REAL val, vall
c-      integer icnt, mxtr
c-      common/dbuf/val(312),vall(195),icnt,mxtr(4)

c-      integer iao, jao, kao, lao, iprt
c-      common/dlabs/iao(312),jao(312),kao(312),lao(312),iprt(312)

INCLUDE(../m4/common/ghfblk)
INCLUDE(../m4/common/cigrad)

cINCLUDE(../m4/common/symtry)
      integer nt
      data nt/1/

INCLUDE(../m4/common/timez)

      REAL dgout
      common/tgrad/dgout(9)

INCLUDE(common/dft_d2escr)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)  ! debug switches
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_dshlno)


INCLUDE(../m4/common/ijlab)

INCLUDE(../m4/common/segm)
INCLUDE(../m4/common/cslosc)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/parallel)

      logical odbg
      common/dbgdbg/odbg
c
c local variables
c      
      integer m1, m2, m3, m0
      dimension m0(48),m1(48),m2(48),m3(48)

      integer najkl, nabcl, nd2str, nactp, nshdm
      integer numcmo, n4
      integer id3, id4, id5, id6, id7, id8, id9
      integer i10, i20, i30, i40, i00
      integer iabd1, iabd2, iabd3
      integer ijshel
      integer i, j, issi, lnddm, lensec, m
      integer kt_max, inc1_max, ncmmm_max
      integer maxll, maxjj, maxkk, it
      logical ofpres, odpres, ogpres

      logical olab, olabc, olabcd
      integer ii, jj, kk, ll
      integer i0, j0, k0, l0
      integer l2,  m00
      integer klshel, kadi, id, jd, kd, ld, nd
      integer iblok
      integer nav, ifmp1, ifmp2, ifmp3
      integer iwor1, iwor2, iwor3
      integer ndgout, next
      integer kadij, kadijk, iceni
      logical omp2, ocifor, omp2w, ohf, ociopt, oskipp
      REAL tolij, tolijk, abmax
      REAL dtim, tim0
      integer nschwz, itrij

      character*8 zmcscf,zinfra,zgvb,zgrhf,zcas,zfock,zrhf

      character*8 zdipd, zuhf
      logical omcscf, ocas, ouhf, orgvb
      integer ib1, isymd
      integer n


      integer ist, jst, kst, lst

      integer icontij, icontkl, ncentr
      integer DENSITY, FIT
      parameter (DENSITY=1, FIT=2)

      REAL gmax

_IF(64bitpointers)
      integer*8 itmp1, itmp2
_ELSE
      integer itmp1, itmp2
_ENDIF
      integer itmp(4), bas

      character *9  fnm
      character *16 snm
c     
c functions/stmnt fns
c
      integer ipg_dlbtask,ipg_nodeid, lenwrd
      integer null_memory, allocate_memory2
      integer iipsci
      logical opg_root, oipsci
c-      integer indq

      data zrhf,zuhf,zgvb,zgrhf/
     * 'rhf','uhf','gvb','grhf' /
      data zcas,zmcscf/'casscf','mcscf'/
      data zdipd/'dipder'/,zinfra/'infrared'/
      data zfock/'fockder'/
      data fnm/'deriv2e.m'/
      data snm/'jkder_dft_genuse'/
c
c-      indq(m,n) = (m-1)*n
c
c-      ompir = runtyp.eq.zdipd .or. runtyp.eq.zinfra
c-      ompir = ompir .and. omp2
c
c     ----- check for grhf or gvb cases -----
c
c-      ouhf = zscftp.eq.zuhf
c-      orgvb = zscftp.eq.zgvb
c-      ogrhf = zscftp.eq.zgrhf
c-      oclos = zscftp.eq.zrhf
c-      omcscf = zscftp.eq.zmcscf
c
c     ----- check for casscf
c
c-      ocas = zscftp.eq.zcas
c-      if (ocas .or. mp3 .or. (omp2 .and. .not.ompir)) then
c-         call setsto(1360,0,ijkl)
c-      end if

c
c disable prefactor testing
c
c      if (nprint.ne.-5 .and. oprint(57)) then
c         write (iwr,6010)
c         call writel(q(iprefa),nshels)
c      end if

c      dummy iso
c@@ need to decide how big this array should
c   be, and only use it when centres are interchangeable
c

      do i = 1,nshels
         iso(i,1) = 1
      enddo

      odbg = .false.
      nschwz = 0
c
c this routine forces the load of the block data
c
       call dummy_intgx
c
c Determine which term to compute
c
      ncentr = 0
      if(basi .gt. 0 .and. basj .gt. 0)then
         icontij  = DENSITY
         ncentr = ncentr + 2
      else
         icontij  = FIT
         ncentr = ncentr + 1
      endif

      if(bask .gt. 0 .and. basl .gt. 0)then
         icontkl  = DENSITY
         ncentr = ncentr + 2
      else
         icontkl  = FIT
         ncentr = ncentr + 1
      endif

      if(icontkl .eq. DENSITY .and. icontij .eq. FIT) then
         call caserr('jkder_dft called incorrectly')
      endif

      if(odbg .and. ncentr .eq. 2)then
         if(basi .eq. bask)then
         write(6,*)'fitting coefficients',(cfit(i),i=1,nbasfn(basi))
         else
         write(6,*)'fitting coefficients 1',(cfit(i),i=1,nbasfn(basi))
         write(6,*)'fitting coefficients 2',(cfit2(i),i=1,nbasfn(basi))
         endif
      endif

      if(print_sw(DEBUG_FORCES) .and. opg_root())then
        write(6,*)'Computing',ncentr,' centre deriv integrals'
      endif
c
c     ----- two electron contribution to the gradient -----
c
      ndgout = 9
c
c     ----- set some parameters -----
c
      nav = lenwrd()

      ntpdm = 1

c-      oeof = .false.

cc      if( icontij .eq. DENSITY ) then

         l2 = (nbasfn(basi)*(nbasfn(basi)+1))/2

cc      endif

c-      if (ompir) then
c-         ndenin = 3
c-         ntpdm = 4
c-         iflden = idaf
c-         call secget(isect(31),31,ibden)
c-         iwor1 = 0
c-         iwor2 = 0
c-         iwor3 = 0
c-         mwor1 = 0
c-         mwor2 = 0
c-         mwor3 = 0
c-         ifmp1 = 20
c-         ifmp2 = 21
c-         ifmp3 = 22
c-         ib1 = 1
c-         call search(ib1,ifmp1)
c-         call search(ib1,ifmp2)
c-         call search(ib1,ifmp3)
c-         call setsto(1360,0,ijkl1)
c-         call setsto(1360,0,ijkl2)
c-         call setsto(1360,0,ijkl3)
c-         call secget(isect(57),57,iblok)
c-         call rdedx(dipd,lds(isect(57)),iblok,ifild)
c-      end if
c-      if (ogrhf) then
c-         m = 0
c-         call secget(isect(53),m,iblok)
c-         call rdedx(nact,lds(isect(53)),iblok,idaf)
c-      end if

      nat3 = natoms*3
c      nbsq = num*num
c      lenb = lensec(nx)

      odpres = .false.
      ofpres = .false.
      ogpres = .false.

      itmp(1) = basi
      itmp(2) = basj
      itmp(3) = bask
      itmp(4) = basl
      do j=1,4
      if(itmp(j) .lt. 0) itmp(j) = itmp(1)
      do 20 i = 1 , nshell(itmp(j))
         bas = itmp(j)
         if (ktype(bas,i).eq.3) odpres = .true.
         if (ktype(bas,i).eq.4) ofpres = .true.
         if (ktype(bas,i).eq.5) ogpres = .true.
 20   continue
      enddo
c

      m = ntpdm + 9
c     if (omp2w .or. ompir) m = m + 3
      lnddm = 257
      issi = lnddm*m + 54*24
      if (odpres) then
         lnddm = 1297
         issi = lnddm*m + 193*30
      end if
      if (ofpres) then
         lnddm = 10001
         issi = lnddm*m + 501*42
      end if
      if (ogpres) then
         lnddm = 50626
         issi = lnddm*m + 501*42
      end if
c
c     code for dgenrl preallocation 
c     buffer space for dgenrl (starts at ic1)
c
      kt_max=0
      do j=1,4
      bas = itmp(j)
      do 21 i = 1 , nshell(bas)
         kt_max=max(kt_max,ktype(bas,i))
 21   continue
      enddo
      inc1_max=(kt_max+1)**4
c
c     vectorisation factor, see also in dgenrl
c
      ncmmm_max = 32
      issi = issi + inc1_max*ncmmm_max*6

c      write(6,*)'dfg',odpres, ofpres,ogpres
c
c     ----- set pointers for partitioning of core -----
c

c i10 used for nconf

      i10 = null_memory()
cc      i20 = igmem_alloc(l2)
cc      i30 = igmem_alloc(l2)

      ndens = 1
c     if (ouhf) then
c        ndens = 2
c     end if
c     if (orgvb) then
c        ndens = 3
c     end if
c     if (ogrhf) then
c        ndens = njk
c     end if
c     if (ocas) then
c        ndens = 1
c     end if
c     if (ofock .or. ompir) then
c        ndens = ndens + ndenin
c     end if

cc      if( icontij .eq. DENSITY ) then

         ida = null_memory()
         i20 = null_memory()
         i30 = null_memory()

cc      else
cc         ida = 0
cc         i20 = 0
cc         i30 = 0
cc      endif

      nfok = 0
      ofock = fkder.eq.zfock
      if (ofock) then
         nfok = nat3
         if (zscftp.eq.zuhf) nfok = nat3 + nat3
c-         if (ogrhf) nfok = njk*natoms*6
      end if
      if (ofokab) nfok = nat3*ndenin
c

c-      ifok = igmem_alloc(nx*nfok)
c-      ifok = ifok - 1

c  core at i00 is for indexing
      i00 = allocate_memory2(lnddm*(3/nav+1),'d',fnm,snm,'i00')

      iabd = allocate_memory2(lnddm*ntpdm,'d',fnm,snm,'iabd')
      iabd = iabd - 1

c      write(6,*)'core i20',i20
c      call chkadr(q(i20))

      ic7 = allocate_memory2(lnddm*9 + inc1_max*ncmmm_max*6,'d',fnm,snm,
     &                       'ic7')
      ic7 = ic7 - 1

cc      ic1 = igmem_alloc(inc1_max*ncmmm_max*6)

c-      if(omcscf)then
c-
c-         i40 = igmem_alloc(l2)
c-
c-ccc         id3 = i40 + nx
c-         numcmo = num*ncoorb
c-         nactp = ncact*(ncact+1)/2
c-         nd2str = indq(nactp+1,16)
c-         if (odpres) nd2str = indq(nactp+1,36)
c-         if (ofpres) nd2str = indq(nactp+1,100)
c-         if (ogpres) nd2str = indq(nactp+1,225)
c-         id3 = igmem_alloc(numcmo)
c-         i30 = i20
c-ccc         id4 = id3 + numcmo
c-         id4 = id3 + numcmo
c-c     nd2mo=ind(nactp,nactp)
c-         nd2mo = nactp*(nactp+1)/2
c-c     write(6,*)' length of tpdm',nd2mo
c-
c-         id4 = igmem_alloc(nd2str)
c-ccc         id5 = id4 + nd2str
c-         id5 = igmem_alloc(nd2mo)
c-ccc         id6 = id5 + nd2mo
c-c
c-         najkl = indq(ncact*4+1,nactp)
c-         if (odpres) najkl = indq(ncact*6+1,nactp)
c-         if (ofpres) najkl = indq(ncact*10+1,nactp)
c-         if (ogpres) najkl = indq(ncact*15+1,nactp)
c-c
c-         nabcl = indq(16+1,4*ncact)
c-         if (odpres) nabcl = indq(36+1,6*ncact)
c-         if (ofpres) nabcl = indq(100+1,10*ncact)
c-         if (ogpres) nabcl = indq(225+1,10*ncact)
c-c
c-         id6 = igmem_alloc(nabcl)
c-cc         id7 = id6 + nabcl
c-         id7 = igmem_alloc(najkl)
c-cc         id8 = id7 + najkl
c-         id8 = igmem_alloc(nactp)
c-cc         id9 = id8 + nactp
c-c
c-         nshdm = max(ncact,4)
c-         if (odpres) nshdm = max(ncact,6)
c-         if (ofpres) nshdm = max(ncact,10)
c-         if (ogpres) nshdm = max(ncact,15)
c-c
c-         id9 = igmem_alloc(nshdm*nshdm)
c-cc
c-cc         iabd = id9 + nshdm*nshdm
c-c
c-c   iabd should not be used for mcscf calculations - tpdm in id5
c-c
c-cc         i00 = iabd + lnddm*ntpdm
c-cc         ic7 = i00 + lnddm*(3/nav+1) - 1
c-      endif

c     tim0 = cpulft(1)
c
c-      if (ocas) then
c-         call dbutci(ist,jst,kst,lst)
c-      else if (omcscf) then
c-         call ddebut(zscftp,q(id3),q(i30),q(i10),q(i40),q(id5),
c-     +               ist,jst,kst,lst,q)
c-      else

      call ddebut_dft(zscftp,q(i20),q(i30),q(i10),q(i10),dummy,
     +     ist,jst,kst,lst,ncentr,q)

c-      end if
c-      if (mp3 .or. (omp2 .and. .not.ompir)) call search(iblk2d,ifil2d)

c-      mword = 0
c-      iword = 0

c-      kworx = 999
c-      kwor1 = 999
c-      kwor2 = 999
c-      kwor3 = 999
c
c check 
c??      kloc(nshels+1) = num + 1

c
c-      icnt = 0
c-      ib1 = 1
c-      if (omp2w) then
c-         call search(ib1,mpstrm(1))
c-         if (odebug(30)) write (iwr,6020)
c-      end if

c
c     ----- initialise static load balancing counter ------
c
      icount_dlb = iipsci()

      if(print_sw(DEBUG_PARALLEL))then
         write(6,*)'Task info on entry: ',ipg_nodeid(),icount_dlb
      endif

      ite3c_shl = 0
      ite2c_shl = 0

      if (ist.le.nshels) then
c
c     ----- ishell -----
c
         do 140 ii = ist , nshell(basi)

c-            kadi = kad(ii)
            ijshel = ii*(ii-1)/2

c
            if(ncentr .eq. 3 .or. ncentr .eq. 4)then
            do 40 it = 1 , nt
               id = iso(ii,it)
               if (id.gt.ii) go to 140
               m0(it) = id
 40         continue
            endif

            iceni = katom(basi,ii)

c-            if (omcscf) call mcajkl(q(id3),q(id5),q(id7),q(id8),q(id9),
c-     +                              ii,nactp)
c
c     ----- jshell -----
c
            if(basj .lt. 0)then
               j0 = 1
               maxjj = 1
            else
               j0 = jst
               maxjj = ii
            endif

            do 130 jj = j0 , maxjj
c-               kadij = kadi + kad(jj)
               jst = 1
               itrij = ijshel+jj

c-               tolij = dlntol + q(iprefa+ijshel+jj-1)
c *** selection on ij parts must not be too heavy. 3.401=ln(30)
c                  - DISABLED -
c-               if (tolij.gt.-3.401d0) then
                 if(.true.)then
c
c apply i/j tests only when i and j are AO basis fns
c
                  if (ncentr .eq. 3  .or. ncentr .eq. 4)then
                  do 60 it = 1 , nt
                     id = m0(it)
                     jd = iso(jj,it)
                     if (jd.gt.ii) go to 130
                     if (id.lt.jd) then
                        nd = id
                        id = jd
                        jd = nd
                     end if
                     if (id.eq.ii .and. jd.gt.jj) go to 130
                     m1(it) = id
                     m2(it) = jd
 60               continue
                  endif

                  if(basj.lt.0)then
                     olab = .true.
                  else
                     olab = katom(basj,jj).eq.iceni
                  endif

c
c     store information about the pair (ij)
c
                  call dshell_dft(1,ii,jj,kk,ll,
     +                 basi, basj, bask, basl,ncentr)

                  call dprim_dft

                  if (nij.ne.0) then

c-                     if (omcscf) call mcabkl(q(id3),q(id7),q(id4),ii,jj,
c-     +                   nactp)

                     icount_dlb = icount_dlb + 1
                     if(.not.oipsci()) then

                        if(print_sw(DEBUG_PARALLEL))then
                         write(6,*)'Task ',icount_dlb,' on node',
     &                             ipg_nodeid()
                        endif

c
c     ----- kshell -----
c
                     if(ncentr .eq. 4 .or. 
     &                 (ncentr .eq. 2 .and. basi .eq. bask))then
                        k0 = kst
                        maxkk = ii
                     else
                        k0 = kst
                        maxkk = nshell(bask)
                     endif

                     if(odbg)write(6,*)'kk loop',k0, maxkk

                     do 120 kk = k0 , maxkk

c-                        kadijk = kadij + kad(kk)
                        kst = 1
                        klshel = kk*(kk-1)/2

                        if( ncentr .eq. 4) then
                        do 80 it = 1 , nt
                           kd = iso(kk,it)
                           if (kd.gt.ii) go to 120
                           m3(it) = kd
 80                     continue
                        endif

                        olabc = olab .and. katom(bask,kk).eq.iceni

c-                        if (omcscf)
c-     +                      call mcabcl(q(id3),q(id4),q(id6),q(id8),
c-     +                      q(id9),ii,jj,kk,nactp)

                        if(basl .lt. 0)then
                           l0 = 1
                           maxll = 1
                        else
                           l0 = lst
                           maxll = kk
                           if (kk.eq.ii) maxll = jj
                        endif

                        do 110 ll = l0 , maxll
                           ite3c_shl = ite3c_shl + 1
                           lst = 1
                           if (ncentr.eq.3) then
                              if(ite3c_stored(ite3c_shl).eq.0) then
                                 nschwz = nschwz + 1
                                 go to 110
                              endif
                           endif

c-                           if (kadijk+kad(ll).lt.0) then
c-                              tolijk = tolij + q(iprefa+klshel+ll-1)
c-                              if (tolijk.gt.0.0d0) then

                           if(basl.lt.0)then
                              olabcd = olabc 
                           else
                              olabcd = olabc .and. 
     &                             katom(basl,ll).eq.iceni
                           endif

                           if(odbg)write(6,*)'olab',ii,jj,kk,ll,olabcd

                                 if (.not.(olabcd)) then

                                 if( ncentr .eq. 4 ) then

                                    n4 = 0
                                    do 100 it = 1 , nt
                                       ld = iso(ll,it)
                                       if (ld.gt.ii) go to 110
                                       kd = m3(it)
                                       if (kd.lt.ld) then
                                         nd = kd ! swap k,l
                                         kd = ld
                                         ld = nd
                                       end if
                                       id = m1(it)
                                       jd = m2(it)
                                       if (id.eq.ii .or. kd.eq.ii) then
                                         if (kd.ge.id) then
                                         if (kd.ne.id .or. ld.gt.jd)
     +                                      then
                                         nd = id  ! swap i,k
                                         id = kd
                                         kd = nd
                                         nd = jd  ! swap j,l
                                         jd = ld
                                         ld = nd
                                         end if
                                         end if
                                         if (jd.ge.jj) then
                                         if (jd.gt.jj) go to 110
                                         if (kd.ge.kk) then
                                         if (kd.gt.kk) go to 110
                                         if (ld.ge.ll) then
                                         if (ld.gt.ll) go to 110
                                         n4 = n4 + 1
                                         end if
                                         end if
                                         end if
                                       end if
 100                                continue
c
c     ----- calculate q4 factor for this group of shells -----
c
                                    q4 = dfloat(nt)/dfloat(n4)

                                 elseif (ncentr .eq. 3) then

c rather empirical -it seems triangulation effects
c are already corrected for
c a factor of 2 is applied in dabab
                                    q4 = 1.0d0

                                 elseif (ncentr .eq. 2) then
c                                   ite2c_shl = ite2c_shl+1
c                                   if (ite2c_stored(ite2c_shl).eq.0)
c    +                                 goto 110
                                    q4 = 1.0d0
                                 endif
c
c     ----- check for redundant combinations -----
c
                              call redund_dft(ii,jj,kk,ll,
     +                             basi, basj, bask, basl, 
     +                             iwr)


                              if (npass.eq.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,jj,kk
                                 else if(ncentr .eq. 2)then
                                    if(odbg)
     &                   write(6,*)'Redund Shell block',ii,kk
                                 endif

                              endif
                              
                              if (npass.ne.0) then

                                 if(ncentr .eq. 4)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,ll
                                 else if(ncentr .eq. 3)then
                                    if(odbg)
     &                        write(6,*)'Shell block',ii,jj,kk,q4
                                 else if(ncentr .eq. 2)then
                                    if(ii.eq.14.and.kk.eq.8)then
c                                       odbg=.true.
                                    else
c                                       odbg = .false.
                                    endif
                                    if(odbg)
     &               write(6,*)'Shell block',ii,kk, q4
                                 endif
c
c     ----- initialize dgout to zero -----
c
                          call vclr(dgout,1,ndgout)

                          call dshell_dft(2,ii,jj,kk,ll,
     +                         basi, basj, bask, basl,ncentr)
c
c     ----- form products of density matrix elements -----
c
                          call vclr(q(iabd+1),1,lendd*ntpdm)

c-                                 if (omcscf)
c-     +                              call mcabcd(q(id3),q(id6),
c-     +                              q(iabd+1),q(id9),ii,jj,kk,ll,
c-     +                              q4)
c-                                 if (mp3 .or.
c-     +                              (omp2 .and. .not.ompir))
c-     +                              call mcdab(q(iabd+1),ii,jj,kk,
c-     +                              ll,q4)
c-                                 if (ompir) then
c-                                   iabd1 = iabd + lendd + 1
c-                                   iabd2 = iabd1 + lendd
c-                                   iabd3 = iabd2 + lendd
c-                                   call dpdab1(q(iabd1),ii,jj,kk,
c-     +                                ll,q4,ifmp1,iwor1)
c-                                   call dpdab2(q(iabd2),ii,jj,kk,
c-     +                                ll,q4,ifmp2,iwor2)
c-                                   call dpdab3(q(iabd3),ii,jj,kk,
c-     +                                ll,q4,ifmp3,iwor3)
c-                                   call tpdder(q(iabd1),lendd,
c-     +                                q(i20),q(i30),nx,ndenin,ii,
c-     +                                jj,kk,ll,q4)
c-                                 end if
c-                                 if (.not.orgvb .and.
c-     +                              .not.ocas .and. .not.omcscf)

                          call dabab_dft(ii,jj,kk,ll,
     &                         basi, basj, bask, basl,
     +                         ncentr,q4,
     +                         zscftp,adens,bdens,cfit,cfit2,
     +                         q(iabd+1))

c                          call chkadr2(q(i20),itmp1)
c                          call chkadr2(q(iabd+1),itmp2)
c                          write(6,*)'after dabab',itmp1,itmp2,
c     &                         q(i20),q(iabd+1)


c-                                 if (orgvb)
c-     +                              call dabg(ii,jj,kk,ll,l1,norb,
c-     +                              q4,q(i20),q(i30),q(i10)
c-     +                              ,onocor,onopen,q(iabd+1))
c-                                 if (ocas)
c-     +                              call dabci(ii,jj,kk,ll,q4,
c-     +                              oeof,q(iabd+1))
c-                                 if (omcscf)
c-     +                              call dabmc(ii,jj,kk,ll,q4,
c-     +                              q(i30),q(i40),q(iabd+1))

c
c     ----mess about with the density matix by eliminating
c     zero elements
c

                                 call delim_dft(q(iabd+1),q(ic7+1),
     +                              ijgt,klgt,q(i00),ijd,kld,
     +                              ijkld,lendd,abmax)

c-                                 if (ompir) then
c-                                   call delim2(q(iabd1),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd2),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                   call delim2(q(iabd3),q(ic7+1),
c-     +                                ijgt,klgt,q(i00),ijd,kld,
c-     +                                ijkld,lendd)
c-                                 end if
                                 
                                 if (ijkld.ne.0) then
                                    call dgenrl_dft(q(1),q(i00),q(i00),
     +                                   abmax,gmax)

c-                                   if (ofock)
c-     +                                call fockd2(q,q(i00))
c-                                   if (ofokab) then
c-                                     call fokabd(q,q(i00))
c-                                   end if
c
c     ----- generate all 4 partial contributions to the gradient ----
c
                                   call formeg_dft
                                 end if              ! ijkld.ne.o
                              end if                 ! olabcd
                                 end if              ! npass .ne. 0

c-                              end if                 ! (tolijk.gt.0.0d0)
c-                           end if                    ! (kadijk+kad(ll).lt.0)

 110                    continue                     ! ll   loop
 120                 continue
                     next = ipg_dlbtask()
                     endif
                  end if
               end if
 130        continue

c
c     ----- save gradient and restart data -----
c                ==== disabled =====
c            call dfinal_dft(q,0,ii,basi,grad,nat)
c            if (tim.ge.timlim) go to 150

 140     continue
         call pg_dlbpush

      end if
c
c     ----- end of *shell* loops -----
c
c-      if (omp2w) then
c-         if (icnt.ne.0) then
c-            call pack(vall,8,iao,1560)
c-            call wrt3s(val,511,mpstrm(1))
c-         end if
c-         icnt = 0
c-         m00 = 0
c-         call put(val,m00,mpstrm(1))
c-         call shut1(mpstrm(1))
c-         if (odebug(30)) write (iwr,6030)
c-      end if

      if (nt.ne.1) then
c
c     ----- allocate core memory for symde
c
       call caserr('no symmetry in jkder_dft')
ccc       isymd = igmem_alloc(nw196(6))
c
ccc       call symde(q(isymd),natoms)

c
c     ----- reset core memory from symde
c
ccc      call gmem_free(isymd)
c
      endif
      call dfinal_dft(q,1,ii,basi,grad, natoms)

c150  continue
      call timit(0)
c     dtim = tim - tim0
c
c     ----- reset core memory from jkder -----
c

c-      if(omcscf)then
c-         call gmem_free(id9)
c-         call gmem_free(id8)
c-         call gmem_free(id7)
c-         call gmem_free(id6)
c-         call gmem_free(id5)
c-         call gmem_free(id4)
c-         call gmem_free(id3)
c-         call gmem_free(i40)
c-      endif

cc
cc revert
cc      call gmem_free(ic1)

      ic7 = ic7 + 1
      call free_memory2(ic7,'d',fnm,snm,'ic7')
      iabd = iabd + 1
      call free_memory2(iabd,'d',fnm,snm,'iabd')
      call free_memory2(i00,'d',fnm,snm,'i00')
c-      ifok = ifok + 1
c-      call gmem_free(ifok)

cc      if( icontij .eq. DENSITY ) then
cc      endif

cc      call gmem_free(i40)
cc      call gmem_free(i30)
cc      call gmem_free(i20)

      return
c6010 format (/40x,'prefactor matrix in 2e-derivatives'/)
c6020 format (/1x,'derivative integrals to be output')
c6030 format (/1x,'derivative integrals written')
      end
      subroutine ddebut_dft(zscftp,da,db,nconf,dd,fock,
     + ista,jsta,ksta,lsta,ncentr,q)
c
c      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
c      implicit character *8 (z),character *1 (x)
c      implicit character *4 (y)

      implicit none
c
c arguments
c
      character*8 zscftp
      REAL da, db, fock, dd, q
      integer nconf
      dimension da(*),db(*),nconf(*),fock(*),dd(*), q(*)
      integer  ista,jsta,ksta,lsta,ncentr

INCLUDE(../m4/common/sizes)

c-      logical ociopt, omp2, ohf
c-      common/restrl/ociopt(2),omp2

c retained for normf normp itol nprint icut 
INCLUDE(../m4/common/restar)

INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_grad2)
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_intctl)
c
c functions
c
      integer ipg_nnodes
c
c local variables
c
      integer ioffv, ioffd, iofft
      integer itold, loop, ioff
      integer nbase, ncol, nop
      integer i, j, ij, k, kk, is, nsi, ni, ic, nc, np2, m, ils, icutd
      integer iblok
      
      REAL factor, dum, duma, dumb

      character*8 zrhf,zuhf,zgrhf,zgvb,zmcscf
      REAL done, ten, e
c
c  Functions
c
      logical opg_root

c
      data zrhf,zuhf,zgrhf,zgvb/'rhf','uhf','grhf','gvb'/
      data zmcscf/'mcscf'/
      data done,ten,e /1.0d0,1.0d1,2.30258d0/
c
c      if (onocnt) write (iwr,6010)
c
c-      if (zscftp.eq.zgvb) then
c-c
c-c     ----- set up core density matrix, read in  eigenvectors for
c-c            gvb -----
c-c
c-         call rdedx(db,l3,ibl3qa,idaf)
c-         call tdown(db,ilifq,db,ilifq,l1)
c-         call dencor(da,db,l1)
c-c
c-c     ----- set up mo to fock operator pointers in nconf -----
c-c
c-         if (nco.ne.0) call setsto(nco,1,nconf)
c-         ic = 0
c-         if (nseto.gt.0) then
c-            nbase = ncores
c-            do 30 i = 1 , nseto
c-               nop = no(i)
c-               do 20 j = 1 , nop
c-                  nconf(ic+nco+j) = nbase + 1
c- 20            continue
c-               ic = ic + nop
c-               nbase = nbase + 1
c- 30         continue
c-         end if
c-         if (npair.gt.0) then
c-            np2 = npair + npair
c-            do 40 i = 1 , np2
c-               nconf(i+nco+ic) = ncores + nseto + i
c- 40         continue
c-         end if
c-         norb = nco + npair + npair + ic
c-         onocor = .false.
c-         onopen = .false.
c-         onocor = nco.eq.0
c-         onopen = nseto.eq.0 .and. npair.eq.0
c-c
c-      else if (zscftp.eq.zuhf) then
c-c
c-c     ----- read in density matrices (alpha+beta) in uhf -----
c-c
c-         call rdedx(da,l2,ibl3pa,idaf)
c-         call rdedx(db,l2,ibl3pb,idaf)
c-         do 50 i = 1 , l2
c-            duma = da(i)
c-            dumb = db(i)
c-            da(i) = duma + dumb
c-            db(i) = duma - dumb
c- 50      continue
c-c
c-      else if (zscftp.eq.zrhf) then
c-c
c-c     ----- read in density matrix  in rhf  -----
c-c
c-_IF(mp2_parallel)
c-c gdf:  skip I/O for parallel mp2 gradient
c-        if (.not.omp2) call rdedx(da,l2,ibl3pa,idaf)
c-_ELSE
c-         call rdedx(da,l2,ibl3pa,idaf)
c-c      write(6,*)'input density',(da(i),i=1,l2)
c-_ENDIF
c-c
c-      else if (zscftp.eq.zgrhf) then
c-c
c-c     general scf
c-c
c-         call rdedx(fock,l3,ibl3qa,idaf)
c-         call tdown(fock,ilifq,fock,ilifq,l1)
c-         m = 0
c-         do 90 is = 1 , njk
c-            call vclr(da(m+1),1,l2)
c-            nsi = nbshel(is)
c-            ils = ilfshl(is)
c-            do 80 ni = 1 , nsi
c-               nc = iactiv(ils+ni)
c-               ncol = (nc-1)*l1
c-               ij = 0
c-               do 70 i = 1 , l1
c-                  dum = fock(i+ncol)
c-                  do 60 j = 1 , i
c-                     ij = ij + 1
c-                     da(ij+m) = da(ij+m) + dum*fock(j+ncol)
c- 60               continue
c- 70            continue
c- 80         continue
c-            m = m + l2
c- 90      continue
c-c
c-      else if (zscftp.eq.zmcscf) then
c-c
c-c     mcscf/multi
c-c
c-         m = 0
c-         call secget(isecmo,m,iblok)
c-         call rdedx(da,l1*ncoorb,iblok+mvadd,idaf)
c-         m = 0
c-         call secget(isecdd,m,iblok)
c-         if(odebug(30)) write (iwr,6040) iblok
c-         call rdedx(dd,l2,iblok,idaf)
c-         if(odebug(30)) write (iwr,6050) ifil2d,iblk2d
c-         call rdedx(fock,nd2mo,iblk2d,ifil2d)
c-         if (ncore.gt.0) then
c-            call vclr(db,1,l2)
c-            ij = 0
c-            do 120 i = 1 , l1
c-               do 110 j = 1 , i
c-                  ij = ij + 1
c-                  do 100 k = 1 , ncore
c-                     kk = (k-1)*l1
c-                     db(ij) = db(ij) + da(i+kk)*da(j+kk)
c- 100              continue
c- 110           continue
c- 120        continue
c-            write (iwr,6020) ncore
c-         end if
c-      else
c-         call caserr('invalid scftype detected in gradient code')
c-      end if


c
c     ----- set starting parameters -----
c
      outd= nprint.eq. - 4

c-      if (ofokab .or. ompir) then
c-c
c-c      read in ndenin density matrices from some external file
c-c      matrices are in mo basis
c-c
c-         ioff = 1
c-         call search(ibden,iflden)
c-         do 130 loop = 1 , ndenin
c-            call reads(db(ioff),l2,iflden)
c-            if (odebug(21)) call prtris(db(ioff),l1,iwr)
c-            ioff = ioff + l2
c- 130     continue
c-
c-         ioffd = igmem_alloc(l2)
c-         iofft = igmem_alloc(l1)
c-         ioffv = igmem_alloc(l3)
c-         m = 0
c-         call rdedx(q(ioffv),l3,ibl3qa,idaf)
c-         call tdown(q(ioffv),ilifq,q(ioffv),ilifq,l1)
c-         ioff = 1
c-         do 140 loop = 1 , ndenin
c-            call dcopy(l2,db(ioff),1,q(ioffd),1)
c-            call demoao(q(ioffd),db(ioff),q(ioffv),q(iofft),l1,
c-     +  ncoorb,l1)
c-            ioff = ioff + l2
c- 140     continue
c-
c-         call gmem_free(ioffv)
c-         call gmem_free(iofft)
c-         call gmem_free(ioffd)
c-
c-      end if

         icutd = icut
         icutd = max(icutd,8)
         itold = itol + 1
c
c        if(opg_root())then
c           write(6,*)'old toler parameters',itold, icutd
c        endif

      if(ncentr .eq. 2) then
         icutd = icutd_2c
         itold = itold_2c
      else if(ncentr .eq. 3) then
         icutd = icutd_3c
         itold = itold_3c
      else if(ncentr .eq. 4) then
         icutd = icutd_4c
         itold = itold_4c
      endif

c     if(opg_root())then
c       write(6,*)'new toler parameters',itold, icutd
c     endif

      itold = max(itold,16)
      dcut = done/(ten**icutd)
      tol1 = e*(itold+2)
      tol2 = e*itold
      tol3 = done/ten**(itold+2)
      tol4 = done/ten**itold
      onorm = normf.ne.1 .or. normp.ne.1

c
c =====  Starting points for shell loops  =====
c      ( replace if implementing restarts)

      ista = 1
      jsta = 1
      ksta = 1
      lsta = 1
c
c ====  zero 2e dft gradient accumulator =====
c
      call dscal(natoms*3,0.0d0,de,1)

c-      if (ofock .or. ofokab) then
c-         if(odebug(30)) write (iwr,6030)
c-         call vclr(fock,1,nfok*l2)
c-      end if

      return
c6010 format (/1x,22('=')/1x,'gradient of the energy'/1x,22('='))
c6020 format (//1x,' number of frozen and core orbitals',i5//)
c6030 format (/1x,'zeroing storage for fock matrices')
c6040 format (/1x,'reading density matrix from block',i5)
c6050 format (/1x,'tpdm from file, block', 2i5)
      end
      subroutine delim_dft(qa,qb,ijgt,klgt,oform,ij,kl,ijkl,lendd,
     * abmax)

      implicit none

      logical oform
      integer ijgt, klgt
      REAL qa, qb
      dimension qa(*),qb(*),ijgt(*),klgt(*),oform(*)
      integer ij, kl, ijkl
      integer lendd

INCLUDE(common/dft_dmisc)

c-      logical ociopt,omptwo,ohf,omp2w,ojunk
c-      common/restrl/ociopt(2),omptwo,ohf(7),omp2w

ccINCLUDE(../m4/common/specal)
c
c local variables
c
      integer i, k, n, nn
      REAL ab, abmax

      integer itmp1

      do 20 i = 1 , lendd
         oform(i) = .false.
 20   continue
      abmax = 1.0d0

c-      if (.not.(ofock .or. ofokab .or. ompir .or. omp2w)) then
         abmax = 0.0d0
         nn = 0
         do 40 i = 1 , ij
            do 30 k = 1 , kl
               nn = nn + 1
               n = ijgt(i) + klgt(k)
               ab = dabs(qa(n))
ccccc               if (ab.lt.dcut) oform(nn) = .true.
               if (ab.gt.abmax) abmax = ab
 30         continue
 40      continue
c-      end if

      call dcopy(lendd,qa,1,qb,1)
      nn = 0
      ijkl = 0
      do 60 i = 1 , ij
         do 50 k = 1 , kl
            nn = nn + 1
            if (.not.(oform(nn))) then
               n = ijgt(i) + klgt(k)
               ijkl = ijkl + 1
               qa(ijkl) = qb(n)
            end if
 50      continue
 60   continue

      return
      end
      subroutine denfac_dft(dkl,csk,cpk,cdk,cfk,cgk,
     +                      csl,cpl,cdl,cfl,cgl,
     +                     mink,maxk,minl,maxl,okandl,double)
      implicit REAL  (a-h,o-z)
       logical okandl,double
       dimension dkl(*)

      logical odbg
      common/dbgdbg/odbg

      if (.not.double) then
         n = 0
         max = maxl
         do 130 k = mink , maxk
            if (okandl) max = k
            go to (20,30,60,60,
     +             40,60,60,60,60,60,
     +             50,60,60,60,60,60,60,60,60,60,
     +             55,60,60,60,60,60,60,60,60,60,
     +             60,60,60,60,60) , k
 20         dum1 = csk
            go to 60
 30         dum1 = cpk
            go to 60
 40         dum1 = cdk
            go to 60
 50         dum1 = cfk
            go to 60
 55         dum1 = cgk
 60         do 120 l = minl , max
               go to (70,80,110,110,
     +                90,110,110,110,110,110,
     +               100,110,110,110,110,110,110,110,110,110,
     +               105,110,110,110,110,110,110,110,110,110,
     +               110,110,110,110,110) , l
 70            dum2 = dum1*csl
               go to 110
 80            dum2 = dum1*cpl
               go to 110
 90            dum2 = dum1*cdl
               go to 110
 100           dum2 = dum1*cfl
               go to 110
 105           dum2 = dum1*cgl
 110           n = n + 1
               dkl(n) = dum2

               if(odbg)write(6,*)'dkl',n,dum2

 120        continue
 130     continue
      else
         n = 0
         max = maxl
         do 250 k = mink , maxk
            if (okandl) max = k
            go to (140,150,180,180,
     +             160,180,180,180,180,180,
     +             170,180,180,180,180,180,180,180,180,180,
     +             175,180,180,180,180,180,180,180,180,180,
     +             180,180,180,180,180) , k
 140        dum1 = csk
            go to 180
 150        dum1 = cpk
            go to 180
 160        dum1 = cdk
            go to 180
 170        dum1 = cfk
            go to 180
 175        dum1 = cgk
 180        do 240 l = minl , max
               go to (190,200,230,230,
     +                210,230,230,230,230,230,
     +                220,230,230,230,230,230,230,230,230,230,
     +                225,230,230,230,230,230,230,230,230,230,
     +                230,230,230,230,230) , l
 190           dum2 = dum1*csl
               if (k.gt.1) then
                  dum2 = dum2 + csk*cpl
               else
                  dum2 = dum2 + dum2
               end if
               go to 230
 200           dum2 = dum1*(cpl+cpl)
               go to 230
 210           dum2 = dum1*(cdl+cdl)
               go to 230
 220           dum2 = dum1*(cfl+cfl)
               go to 230
 225           dum2 = dum1*(cgl+cgl)
 230           n = n + 1
               dkl(n) = dum2
               if(odbg)write(6,*)'dkl',n,dum2
 240        continue
 250     continue
      end if
      return
      end
      subroutine dfinal_dft(q,index,ii,basi,grad, nat)

      implicit none

      REAL q(*)
      integer index, ii, basi
      REAL grad(3,*)
      integer nat

INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_parameters)


c***   ***node-MPP***
c***    writing (and gopping) for index ne 0 is disabled
c***    (or if other integrals still follow)
c***    so no restarts (nb. and no d-gradients!!)
c***    this saves quite a bit of bother
c***    restarts may be reeabled by dividing the partial grads by
c***    the number of nodes after gopping before writing
c***    (the restarting job must have same number then (tricky)
c***   ***node-MPP***


INCLUDE(common/dft_iofile)
INCLUDE(common/dft_grad2)
INCLUDE(common/dft_module_comm)

c
c functions
c
      logical opg_root
c
c local variables
c
      integer n, i, j
      integer min, max, ncoord

      character*4 ydnam
      dimension ydnam(3)

c-      REAL cpu, cpulft
c-      REAL f
c-      integer lensec
c-      integer istnu, jstnu, kstnu, lstnu
c-      integer iblok
c-      integer k
c 
      data ydnam /'e/x','e/y','e/z'/
c
c
c Restart block - commented out in jkder_dft
c
      if (index.ne.1) then
c
c     ----- get restart data -----
c
c         irest = 6
c         istnu = 1 + ii
c         jstnu = 1
c         kstnu = 1
c         lstnu = 1
c
c     ----- save gradient + restart data -----
c
c       call wrtgrd(de,irest,istnu,jstnu,kstnu,lstnu)
c
c     ----- check cpu time -----
c
c        if (istnu.gt.nshell(basi)) return
c          call texit(0,irest)
c          if (tim.lt.timlim) return
c         write (iwr,6010) tim , istnu , jstnu , kstnu , lstnu
c
c-         if (ofokab) call dabout(q,odebug(21),iwr)
c-         if (ofock) call hsymd(q,iwr)
c         write (iwr,6030)
c         return

c-      else if (onocnt) then
      else
c
c     ----- transfer gradient into -grad- -----
c
         ncoord = 3*nat
c
c sum contributions
c
         call pg_dgop(902,de,ncoord,'+')

         if(opg_root() .and. print_sw(DEBUG_FORCES) )then

            max = 0
 20         min = max + 1
            max = max + 8
            if (max.gt.nat) max = nat
            write(iwr,*)'Coulomb gradient from DFT'
            write (iwr,6020)
            write (iwr,6050) (i,i=min,max)
            write (iwr,6020)
            do 30 n = 1 , 3
               write (iwr,6060) ydnam(n) , (de(n,i),i=min,max)
 30         continue

            if (max.lt.nat) go to 20
         endif
c
c Add into the total gradient
c
         do i = 1, nat
            do j = 1,3
               grad(j,i) = grad(j,i) + de(j,i)
            enddo
         enddo
c
c-         if (ofokab) then
c-            call dabout(q,odebug(21),iwr)
c-            iochf(1) = iochf(1) + nat*ndenin*3*lensec(nx)
c-         end if
c-         if (ofock) then
c-            call hsymd(q,iwr)
c-            call clredx
c-         end if
c-         if (ompir) then
c-            call secget(isect(57),57,iblok)
c-            do 60 i = 1 , 3
c-               do 50 j = 1 , 3
c-                  do 40 k = 1 , nat
c-                     dipd(i,j,k) = dipd(i,j,k) + dipi(i,j,k)
c- 40               continue
c- 50            continue
c- 60         continue
c-            call wrt3(dipd,lds(isect(57)),iblok,idaf)
c-         end if
c-         call revind
c-         call clredx
c-         cpu = cpulft(1)
c-         if (nprint.ne.-5) write (iwr,6040) cpu
         return
      end if

c6010 format (/' insufficient time to complete evaluation of 2-electron'
c    +        ,' contribution to gradient'//' job dumped at ',f10.2,
c    +        ' seconds'//' you can forget next batch ',4i5)
 6020 format (/)
c6030 format (//10x,27('*')/10x,'*** warning ***'/10x,
c    +        'this job must be restarted'/10x,27('*')/)
c6040 format (/' end of calculation of the energy gradient at ',f8.2,
c    +        ' seconds'/)
 6050 format (1x,'atom',8(6x,i3,6x))
 6060 format (3x,a3,8f15.7)
      end
      subroutine dform_dft(x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
c      implicit REAL  (a-h,o-z)
      implicit none

c
c arguments
c
      integer ix(*),iy(*),iz(*)
      integer ncdim
      REAL x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     1 yd(ncdim,*),zd(ncdim,*),g(*)

INCLUDE(../m4/common/sizes)
c
INCLUDE(common/dft_root)
INCLUDE(common/dft_d2escr)

INCLUDE(common/dft_incrd)
INCLUDE(common/dft_dshlno)

      REAL t1, t2, t3
      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
c local variables
c
      logical unroll
      integer mx, my, mz
      integer n1, n2, n3, nr, n
_IF(cray,t3d,t3e)
      REAL `ssum'
_ELSEIF(hp700)
      REAL `vec_$dsum'
_ELSE
      REAL `dsum'
_ENDIF
c
cvd$r assoc
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5 .and. ijkld.ge.16
      if (.not.unroll .or. ncontr.gt.5) then
         do 30 n = 1 , ijkld
            mx = ix(n)
            my = iy(n)
            mz = iz(n)
            do 20 nr = 1 , ncontr
               t1(nr) = xd(nr,mx)*y(nr,my)*z(nr,mz)
               t2(nr) = x(nr,mx)*yd(nr,my)*z(nr,mz)
               t3(nr) = x(nr,mx)*y(nr,my)*zd(nr,mz)
 20         continue
            g(n1+n) = g(n1+n) + dsum(ncontr,t1,1)
            g(n2+n) = g(n2+n) + dsum(ncontr,t2,1)
            g(n3+n) = g(n3+n) + dsum(ncontr,t3,1)
 30      continue
         return
      else
         go to (40,60,80,100,120) , ncontr
      end if
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 40   do 50 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz)
 50   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 60   do 70 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz)
 70   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 80   do 90 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
         g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +             *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz)
         g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +             *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz)
         g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +             *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz)
 90   continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 100  do 110 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +            *y(4,my)*z(4,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +            *yd(4,my)*z(4,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +            *y(4,my)*zd(4,mz)
 110  continue
      return
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 120  do 130 n = 1 , ijkld
         mx = ix(n)
         my = iy(n)
         mz = iz(n)
        g(n+n1) = g(n+n1) + xd(1,mx)*y(1,my)*z(1,mz) + xd(2,mx)*y(2,my)
     +            *z(2,mz) + xd(3,mx)*y(3,my)*z(3,mz) + xd(4,mx)
     +            *y(4,my)*z(4,mz) + xd(5,mx)*y(5,my)*z(5,mz)
        g(n+n2) = g(n+n2) + x(1,mx)*yd(1,my)*z(1,mz) + x(2,mx)*yd(2,my)
     +            *z(2,mz) + x(3,mx)*yd(3,my)*z(3,mz) + x(4,mx)
     +            *yd(4,my)*z(4,mz) + x(5,mx)*yd(5,my)*z(5,mz)
        g(n+n3) = g(n+n3) + x(1,mx)*y(1,my)*zd(1,mz) + x(2,mx)*y(2,my)
     +            *zd(2,mz) + x(3,mx)*y(3,my)*zd(3,mz) + x(4,mx)
     +            *y(4,my)*zd(4,mz) + x(5,mx)*y(5,my)*zd(5,mz)
 130  continue
      return
      end
      subroutine dforma_dft(spij,spkl,noform,
     *                  x,y,z,xd,yd,zd,g,ix,iy,iz,ncdim)
c      implicit REAL  (a-h,o-z)
      implicit none
c
c arguments
c
      REAL x,y,z,xd,yd,zd,g
      integer ix,iy,iz
      dimension ix(*),iy(*),iz(*)
      logical noform,spij,spkl
      dimension noform(*)
      integer ncdim
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),xd(ncdim,*),
     & yd(ncdim,*),zd(ncdim,*),g(*)

INCLUDE(../m4/common/sizes)

INCLUDE(common/dft_root)
INCLUDE(common/dft_d2escr)

INCLUDE(common/dft_incrd)
INCLUDE(common/dft_dshlno)

      common/small/t1(ncmax),t2(ncmax),t3(ncmax)
c
c local variables
c
      logical unroll
      integer n, n1, n2, n3, nr, nn, i, k
      integer mx, my, mz
      REAL t1, t2, t3, s1, s2, s3
_IF(cray,t3d,t3e)
      REAL `ssum'
_ELSEIF(hp700)
      REAL `vec_$dsum'
_ELSE
      REAL `dsum'
_ENDIF
c
cvd$r assoc
c
      n = 0
      nn = 0
      n1 = ioff
      n2 = n1 + lendd
      n3 = n2 + lendd
      unroll = ncontr.le.5
      if (.not.unroll) then
c
         if (.not.spij) then
            do 40 i = 1 , ijd
               do 30 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     do 20 nr = 1 , ncontr
                        t1(nr) = ddkl(nr,k)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        t2(nr) = ddkl(nr,k)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        t3(nr) = ddkl(nr,k)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 20                  continue
                     g(n+n1) = dsum(ncontr,t1,1) + g(n+n1)
                     g(n+n2) = dsum(ncontr,t2,1) + g(n+n2)
                     g(n+n3) = dsum(ncontr,t3,1) + g(n+n3)
                  end if
 30            continue
 40         continue
         else if (.not.spkl) then
            do 70 i = 1 , ijd
               do 60 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 50 nr = 1 , ncontr
                        s1 = s1 + ddij(nr,i)*xd(nr,mx)*y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*x(nr,mx)*yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*x(nr,mx)*y(nr,my)*zd(nr,mz)
 50                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 60            continue
 70         continue
         else
            do 100 i = 1 , ijd
               do 90 k = 1 , kld
                  nn = nn + 1
                  if (.not.(noform(nn))) then
                     n = n + 1
                     mx = ix(n)
                     my = iy(n)
                     mz = iz(n)
                     s1 = 0.0d0
                     s2 = 0.0d0
                     s3 = 0.0d0
                     do 80 nr = 1 , ncontr
                        s1 = s1 + (ddij(nr,i)*ddkl(nr,k))*xd(nr,mx)
     +                       *y(nr,my)*z(nr,mz)
                        s2 = s2 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *yd(nr,my)*z(nr,mz)
                        s3 = s3 + ddij(nr,i)*ddkl(nr,k)*x(nr,mx)
     +                       *y(nr,my)*zd(nr,mz)
 80                  continue
                     g(n+n1) = s1 + g(n+n1)
                     g(n+n2) = s2 + g(n+n2)
                     g(n+n3) = s3 + g(n+n3)
                  end if
 90            continue
 100        continue
         end if
         return
      else
         go to (110,140,170,200,230) , ncontr
      end if
 110  do 130 i = 1 , ijd
_IF1(x)c$dir scalar
         do 120 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz)
            end if
 120     continue
 130  continue
      return
 140  do 160 i = 1 , ijd
_IF1(x)c$dir scalar
         do 150 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz)
            end if
 150     continue
 160  continue
      return
 170  do 190 i = 1 , ijd
_IF1(x)c$dir scalar
         do 180 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz)
            end if
 180     continue
 190  continue
      return
 200  do 220 i = 1 , ijd
_IF1(x)c$dir scalar
         do 210 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz)
            end if
 210     continue
 220  continue
      return
 230  do 250 i = 1 , ijd
_IF1(x)c$dir scalar
         do 240 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               n = n + 1
               mx = ix(n)
               my = iy(n)
               mz = iz(n)
               g(n+n1) = g(n+n1) + (ddij(1,i)*ddkl(1,k))*xd(1,mx)
     +                   *y(1,my)*z(1,mz) + (ddij(2,i)*ddkl(2,k))
     +                   *xd(2,mx)*y(2,my)*z(2,mz)
     +                   + (ddij(3,i)*ddkl(3,k))*xd(3,mx)*y(3,my)
     +                   *z(3,mz) + (ddij(4,i)*ddkl(4,k))*xd(4,mx)
     +                   *y(4,my)*z(4,mz) + (ddij(5,i)*ddkl(5,k))
     +                   *xd(5,mx)*y(5,my)*z(5,mz)
               g(n+n2) = g(n+n2) + ddij(1,i)*ddkl(1,k)*x(1,mx)*yd(1,my)
     +                   *z(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*yd(2,my)
     +                   *z(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*yd(3,my)
     +                   *z(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*yd(4,my)
     +                   *z(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*yd(5,my)
     +                   *z(5,mz)
               g(n+n3) = g(n+n3) + ddij(1,i)*ddkl(1,k)*x(1,mx)*y(1,my)
     +                   *zd(1,mz) + ddij(2,i)*ddkl(2,k)*x(2,mx)*y(2,my)
     +                   *zd(2,mz) + ddij(3,i)*ddkl(3,k)*x(3,mx)*y(3,my)
     +                   *zd(3,mz) + ddij(4,i)*ddkl(4,k)*x(4,mx)*y(4,my)
     +                   *zd(4,mz) + ddij(5,i)*ddkl(5,k)*x(5,mx)*y(5,my)
     +                   *zd(5,mz)
            end if
 240     continue
 250  continue
      return
      end
      subroutine dgenrl_dft(qq,iqq,noform,abmax,gmax)

      implicit none
c
c arguments
c
      integer iqq(*)
      REAL qq(*), abmax
      logical noform(*)

      REAL gmax

INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_dshlnf)
INCLUDE(common/dft_d2escr)
INCLUDE(common/dft_dshlno)
INCLUDE(common/dft_root)

      REAL bp01,b00,b10,xcp00,xc00,ycp00,yc00,zcp00
      REAL zc00,f00,dxij,dyij,dzij,dxkl,dykl,dzkl
      integer in,kn,ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2
      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     +   xc00(ncmax),ycp00(ncmax),yc00(ncmax),zcp00(ncmax),
     +   zc00(ncmax),f00(ncmax),
     +   dxij,dyij,dzij,dxkl,dykl,dzkl,
     +   in(12),kn(12),ni,nj,nk,nl,nmax,mmax,ij1,ij2,kl1,kl2

INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_incrd)

      REAL axak,ayak,azak,axai,
     1  ayai,azai,abv,aandbv,rhov,
     2  xxv,c1xv,c2xv,c3xv,c4xv,
     3  c1yv,c2yv,c3yv,c4yv,
     4  c1zv,c2zv,c3zv,c4zv,expev
      common/bufb/axak(mxp2),ayak(mxp2),azak(mxp2),axai(mxp2),
     1  ayai(mxp2),azai(mxp2),abv(mxp2),aandbv(mxp2),rhov(mxp2),
     2  xxv(mxp2),c1xv(mxp2),c2xv(mxp2),c3xv(mxp2),c4xv(mxp2),
     3  c1yv(mxp2),c2yv(mxp2),c3yv(mxp2),c4yv(mxp2),
     4  c1zv(mxp2),c2zv(mxp2),c3zv(mxp2),c4zv(mxp2),expev(mxp2)

ccccINCLUDE(common/segm)

c
c local variables
c
      REAL xb, yb, zb
      REAL bxbi, bxbk
      REAL bybi, bybk
      REAL bzbi, bzbk
      REAL brrk, bbrrk
      REAL csl, cpl, cdl, cgl, cfl
      REAL csk, cpk, cdk, cgk, cfk
      REAL akxk, akyk, akzk
      REAL xc, yc, zc
      REAL xd, yd, zd
      REAL exkl, dddd, dijd, dkld
      REAL ai, aj, ak, al, u2, dum, dum2, expe
      REAL b, binv
      integer ig, jg, lg, i, ii, iii, is, inc
      integer k, n, ncmmm, nn, nnn0, m
      integer max, maxlg, kg, jgmax

      integer itmp1, itmp2

      REAL one
      REAL pi252

      logical trduij,spij,spkl,trdukl
      logical double
      logical sptru

      logical odbg
      common/dbgdbg/odbg


      data one/1.0d0/
c     data zero,pt5/0.0d0,0.5d0/
      data pi252/34.986836655250d0/
c
      if (ijkld.eq.0) return
      if (ijd.eq.1 .and. kld.eq.1) then
         call ssdss_dft(qq)
      else
         ni = lit
         if (oskip(1)) ni = lit - 1
         nj = ljt
         if (oskip(2)) nj = ljt - 1
         nk = lkt
         if (oskip(3)) nk = lkt - 1
         nl = llt
         if (oskip(4)) nl = llt - 1
         kln2 = 1
         kln1 = nl + 1
         ijn2 = kln1*(nk+1)
         ijn1 = ijn2*(nj+1)
         inc1 = ijn1*(ni+1)
c     if(mod(inc1,4).eq.0)inc1=inc1+1
         if (ni.lt.nj) then
            is = ni
            ni = nj
            nj = is
            ij1 = ijn2
            ij2 = ijn1
            xc = xj
            yc = yj
            zc = zj
            dxij = xj - xi
            dyij = yj - yi
            dzij = zj - zi
         else
            ij1 = ijn1
            ij2 = ijn2
            xc = xi
            yc = yi
            zc = zi
            dxij = xi - xj
            dyij = yi - yj
            dzij = zi - zj
         end if
         if (nk.lt.nl) then
            is = nl
            nl = nk
            nk = is
            kl1 = kln2
            kl2 = kln1
            xd = xl
            yd = yl
            zd = zl
            dxkl = xl - xk
            dykl = yl - yk
            dzkl = zl - zk
         else
            xd = xk
            yd = yk
            zd = zk
            dxkl = xk - xl
            dykl = yk - yl
            dzkl = zk - zl
            kl1 = kln1
            kl2 = kln2
         end if
         nmax = ni + nj
         mmax = nk + nl
         max = nmax + 1
         do 20 i = 1 , max
            n = i - 1
            if (n.le.ni) in(i) = ij1*n + 1
            if (n.gt.ni) in(i) = ij1*ni + ij2*(n-ni) + 1
 20      continue
         max = mmax + 1
         do 30 k = 1 , max
            n = k - 1
            if (n.le.nk) kn(k) = kl1*n
            if (n.gt.nk) kn(k) = kl1*nk + kl2*(n-nk)
 30      continue
c
c     indexing
c
         call indxa_dft(ijx,ijy,ijz,ijd,mini,maxi,minj,maxj,oiandj,
     +        ijn1, ijn2,1)
         call indxa_dft(klx,kly,klz,kld,mink,maxk,minl,maxl,okandl,
     +        kln1, kln2,0)

         nn = 0
         ijkld = 0
         do 50 i = 1 , ijd
            do 40 k = 1 , kld
               nn = nn + 1
               if (.not.(noform(nn))) then
                  ijkld = ijkld + 1
                  iqq(ijkld+ixi-1) = ijx(i) + klx(k)
                  iqq(ijkld+iyi-1) = ijy(i) + kly(k)
                  iqq(ijkld+izi-1) = ijz(i) + klz(k)
               end if
 40         continue
 50      continue
c
         do 60 n = 1 , nij
            axak(n) = aa(n)*(x1(n)-xd)
            ayak(n) = aa(n)*(y1(n)-yd)
            azak(n) = aa(n)*(z1(n)-zd)
            axai(n) = aa(n)*(x1(n)-xc)
            ayai(n) = aa(n)*(y1(n)-yc)
            azai(n) = aa(n)*(z1(n)-zc)

            if(odbg)write(6,*)'axak etc',n,axak(n),ayak(n),azak(n)
            if(odbg)write(6,*)'axai etc',n,axai(n),ayai(n),azai(n)

 60      continue

c
c
         trduij = lit.ge.3 .or. ljt.ge.3
         trdukl = lkt.ge.3 .or. llt.ge.3
         spkl = (mink.eq.1 .and. maxk.eq.4) .or.
     +          (minl.eq.1 .and. maxl.eq.4)
         spij = (mini.eq.1 .and. maxi.eq.4) .or.
     +          (minj.eq.1 .and. maxj.eq.4)
         sptru = spkl .or. spij
c
c Now dgenrl uses preallocated memory, max size defined 
c in jkder only pointers ic2 upwards are allocated here
c vect factor of 32 is hardwired
c
c NB there is room above ic7 assuming npass=3, lendd,=lnddm
c
          ncmmm = 32
          ncmmm = (ncmmm/nroots)*nroots

cccccc old memory algorithm cccccccccccccc
c
c     integrals stored at (ic7+1)
c     subsidiary integrals from ic1 onwards
c
cc         ic1 = ic7 + 1 + npass*lendd*3
cc         ncmmm = (nmaxly-ic1-1)/(inc1*6)
cc         ncnnn = ncmax - 1
cc         if (ncmmm.gt.ncnnn) ncmmm = ncnnn
cc         ncmmm = (ncmmm/nroots)*nroots
cc         if (ncmmm.lt.nroots) call caserr('insufficient core in dgenrl')
cc
ccccccc end old memory algorithm cccccccccccc
c
         ic1 = ic7 + 1 + npass*lendd*3
         ic2 = inc1*ncmmm + ic1
         ic3 = inc1*ncmmm + ic2
         ic4 = inc1*ncmmm + ic3
         ic5 = inc1*ncmmm + ic4
         ic6 = inc1*ncmmm + ic5
c
         inc = npass*lendd*3
         call vclr(qq(ic7+1),1,inc)
         ncontr = 0
c
         maxlg = ngd
         do 170 kg = 1 , ngc
            ak = cgg(kg)
            brrk = ak*rrk
            akxk = ak*xk
            akyk = ak*yk
            akzk = ak*zk
            csk = csc(kg)*pi252
            cpk = cpc(kg)*pi252
            cdk = cdc(kg)*pi252
            cfk = cfc(kg)*pi252
            cgk = cgc(kg)*pi252
c
c     ----- l primitive
c
            if (okandl) maxlg = kg
            do 160 lg = 1 , maxlg
               al = dg(lg)
               b = ak + al
               binv = one/b
               bbrrk = al*brrk*binv
               if ((bbrrk+rsmall).le.tol1) then
                  exkl = dexp(-bbrrk)
                  csl = csd(lg)*binv
                  cpl = cpd(lg)*binv
                  cdl = cdd(lg)*binv
                  cfl = cfd(lg)*binv
                  cgl = cgd(lg)*binv
                  xb = (akxk+al*xl)*binv
                  yb = (akyk+al*yl)*binv
                  zb = (akzk+al*zl)*binv
                  bxbk = b*(xb-xd)
                  bybk = b*(yb-yd)
                  bzbk = b*(zb-zd)
                  bxbi = b*(xb-xc)
                  bybi = b*(yb-yc)
                  bzbi = b*(zb-zc)
c
c     ----- density factor
c
                  double = okandl .and. kg.ne.lg

                  call denfac_dft(dkl,csk,cpk,cdk,cfk,cgk,
     +                        csl,cpl,cdl,cfl,cgl,
     +                        mink,maxk,minl,maxl,okandl,double)

                  dkld = dkl(1)
                  if (sptru) then
                     do 70 k = 1 , kld
                        dkl(k) = dkl(k)/dkld
 70                  continue
                  end if
                  dkld = dkld*exkl
                  if (dabs(dkld*abmax).ge.tol3) then
c
c     ----- pair of i,j primitives
c
                     do 80 n = 1 , nij
                        abv(n) = aa(n)*b
                        aandbv(n) = aa(n) + b
                        expev(n) = exij(n)/dsqrt(aa(n)+b)
                        rhov(n) = abv(n)/aandbv(n)
                        xxv(n) = rhov(n)
     +                           *((x1(n)-xb)**2+(y1(n)-yb)**2+(z1(n)
     +                           -zb)**2)
                        c1xv(n) = bxbk + axak(n)
                        c2xv(n) = bxbk*aa(n)
                        c3xv(n) = bxbi + axai(n)
                        c4xv(n) = b*axai(n)
                        c1yv(n) = bybk + ayak(n)
                        c2yv(n) = bybk*aa(n)
                        c3yv(n) = bybi + ayai(n)
                        c4yv(n) = b*ayai(n)
                        c1zv(n) = bzbk + azak(n)
                        c2zv(n) = bzbk*aa(n)
                        c3zv(n) = bzbi + azai(n)
                        c4zv(n) = b*azai(n)
 80                  continue
c
                     n = 0
                     nn = 0
                     jgmax = ngb
                     do 150 ig = 1 , nga
                        ai = ag(ig)
                        if (oiandj) jgmax = ig
                        do 140 jg = 1 , jgmax
                           n = n + 1
                           if ((bbrrk+r(n)).lt.tol2) then
                              aj = bg(jg)
                              dijd = dd(nn+1)
                              if (sptru) then
                                 dddd = one/dijd
                                 do 90 i = 1 , ijd
                                    dij(i) = dd(ijden(i)+nn)*dddd
 90                              continue
                              end if
                              expe = dkld*dijd*expev(n)
                              if (dabs(expe*abmax).ge.tol4) then
                                 pp = xxv(n)
c
c     ----- roots and weights for quadrature
c
                                 if (nroots.le.3) call rt123_dft
                                 if (nroots.eq.4) call roots4_dft
                                 if (nroots.eq.5) call roots5_dft
                                 if (nroots.gt.5) call rootss_dft
c
c     compute two-electron  integrals for each root
c
                                 nnn0 = ncontr
_IF1(x)c$dir scalar
_IF1(ct)cdir$ nextscalar
                                 do 100 m = 1 , nroots
                                    ncontr = ncontr + 1
                                    u2 = u(m)*rhov(n)
                                    f00(ncontr) = expe*w(m)
                                    dum2 = 0.5d0/(abv(n)+u2*aandbv(n))
                                    dum = dum2 + dum2
                                    bp01(ncontr) = (aa(n)+u2)*dum2
                                    b00(ncontr) = u2*dum2
                                    b10(ncontr) = (b+u2)*dum2
                                    xcp00(ncontr) = (u2*c1xv(n)+c2xv(n))
     +                                 *dum
                                    xc00(ncontr) = (u2*c3xv(n)+c4xv(n))
     +                                 *dum
                                    ycp00(ncontr) = (u2*c1yv(n)+c2yv(n))
     +                                 *dum
                                    yc00(ncontr) = (u2*c3yv(n)+c4yv(n))
     +                                 *dum
                                    zcp00(ncontr) = (u2*c1zv(n)+c2zv(n))
     +                                 *dum
                                    zc00(ncontr) = (u2*c3zv(n)+c4zv(n))
     +                                 *dum

                                    aei(ncontr) = ai
                                    aej(ncontr) = aj
                                    aek(ncontr) = ak
                                    ael(ncontr) = al

                            if(odbg)write(6,*)'aei etc',m, aei(ncontr),
     &                    aej(ncontr), aek(ncontr), ael(ncontr)
                            if(odbg)write(6,*)'xcp00',m, xcp00(ncontr),
     &                    ycp00(ncontr), zcp00(ncontr)
                            if(odbg)write(6,*)'xc00',m, xc00(ncontr),
     &                    yc00(ncontr), zc00(ncontr)


 100                             continue

                                 if (sptru) then
                                    ncontr = nnn0
                                    do 130 m = 1 , nroots
                                       ncontr = ncontr + 1
                                       do 110 iii = 1 , ijd
                                         ddij(ncontr,iii) = dij(iii)
 110                                   continue
                                       do 120 iii = 1 , kld
                                         ddkl(ncontr,iii) = dkl(iii)
 120                                   continue
 130                                continue
                                 end if

c
c
c----------------------------------------------
c    defer assembly stage until loop lengths
c    are long enough to vectorise effectively
c----------------------------------------------
c
                                 if (ncontr.ge.ncmmm) then
c     ----- form (i,j//k,l) integrals
c
                                 call dxyz_dft(qq(ic1),qq(ic2),qq(ic3),
     +                                   ncmmm)
                                    ioff = 0
                                    if (.not.(oskip(1))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aei,ljt,lkt,llt,lit,
     +                                    ijn2,kln1,kln2,ijn1,ncmmm)
                                       if (sptru) then
                                     call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                   qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                   iqq(izi),ncmmm)
                                       else
                                     call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(2))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aej,lit,lkt,llt,ljt,
     +                                    ijn1,kln1,kln2,ijn2,ncmmm)
                                       if (sptru) then
                                      call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                     qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(3))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),aek,lit,ljt,llt,lkt,
     +                                    ijn1,ijn2,kln2,kln1,ncmmm)
                                       if (sptru) then
                                     call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                    qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                       ioff = ioff + lendd*3
                                    end if
                                    if (.not.(oskip(4))) then
                                       call subsd_dft(qq(ic1),qq(ic2),
     +                                    qq(ic3),qq(ic4),qq(ic5),
     +                                    qq(ic6),ael,lit,ljt,lkt,llt,
     +                                    ijn1,ijn2,kln1,kln2,ncmmm)
                                       if (sptru) then
                                      call dforma_dft(spij,spkl,noform,
     +                                      qq(ic1),qq(ic2),qq(ic3),
     +                                      qq(ic4),qq(ic5),qq(ic6),
     +                                    qq(ic7+1),iqq(ixi),iqq(iyi),
     +                                      iqq(izi),ncmmm)
                                       else
                                       call dform_dft(qq(ic1),qq(ic2),
     +                                      qq(ic3),qq(ic4),qq(ic5),
     +                                      qq(ic6),qq(ic7+1),iqq(ixi),
     +                                      iqq(iyi),iqq(izi),ncmmm)
                                       end if
                                    end if
                                    ncontr = 0
                                 end if
                              end if
                           end if
c
c
c     end of loops over primitives
c
                           nn = nn + 4
 140                    continue
 150                 continue
                  end if
               end if
 160        continue
 170     continue
c
c    tidy up any bits left unassembled
c
c
         if (ncontr.ne.0) then
            call dxyz_dft(qq(ic1),qq(ic2),qq(ic3),ncmmm)
            ioff = 0
            if (.not.(oskip(1))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aei,ljt,lkt,llt,lit,ijn2,kln1,kln2,
     +                    ijn1,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(2))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aej,lit,lkt,llt,ljt,ijn1,kln1,kln2,
     +                    ijn2,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(3))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),aek,lit,ljt,llt,lkt,ijn1,ijn2,kln2,
     +                    kln1,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                       iqq(izi),ncmmm)
               end if
               ioff = ioff + lendd*3
            end if
            if (.not.(oskip(4))) then
               call subsd_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),qq(ic5),
     +                    qq(ic6),ael,lit,ljt,lkt,llt,ijn1,ijn2,kln1,
     +                    kln2,ncmmm)
               if (sptru) then
                  call dforma_dft(spij,spkl,noform,qq(ic1),qq(ic2),
     +                 qq(ic3),qq(ic4),qq(ic5),qq(ic6),qq(ic7+1),
     +                 iqq(ixi),iqq(iyi),iqq(izi),ncmmm)
               else
                  call dform_dft(qq(ic1),qq(ic2),qq(ic3),qq(ic4),
     +                 qq(ic5),qq(ic6),qq(ic7+1),iqq(ixi),iqq(iyi),
     +                 iqq(izi),ncmmm)
               end if
            end if
            ncontr = 0
         end if
c
c ---------------------------------------------
c    fiddle about if first two centres are same
c ---------------------------------------------
c
         if (natomd(1).eq.natomd(2)) then
            inc = lendd*3
            ii = ic7 + inc
            do 180 i = 1 , inc
               qq(ic7+i) = qq(ic7+i) + qq(ii+i)
 180        continue
            iii = ii + inc
            do 190 i = 1 , inc
               qq(ii+i) = qq(iii+i)
 190        continue
            npass = npass - 1
            natomd(2) = natomd(3)
            natomd(3) = natomd(4)
            natomd(4) = 0
         end if
c
c    insert proper normalisation factors in d-functions
c    present (f- and g-functions also)
c

         if (trduij .or. trdukl) then
            if (onorm) then
               call dnorm_dft(qq(ic7+1),noform)
            end if
         end if
      end if
c
c    multiply by density matrix elements
c
      call grdcon_dft(qq,noform,gmax)
      return
      end
      subroutine dnorm_dft(qq,noform)

      implicit none
c
c arguments
c
      REAL qq(*)
      logical noform(*)

INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_dshlno)

INCLUDE(common/dft_d2escr)

INCLUDE(common/dft_incrd)
INCLUDE(common/dft_picon)

c
c  local variables
c
      integer max, i, j, k, l, n, nn, nnn, npp
      REAL dum1, dum2, d1

      REAL one
      data one/1.0d0/

      n = 0
      max = maxj
      dum1 = one
      do 30 i = mini , maxi
         if (i.eq.8)  dum1 = root3
         if (i.eq.14) dum1 = root5
         if (i.eq.20) dum1 = dum1*root3
         if (i.eq.24) dum1 = root7
         if (i.eq.30) dum1 = dum1*root53
         if (i.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (oiandj) max = i
         do 20 j = minj , max
            if (j.eq.8)  dum2 = dum2*root3
            if (j.eq.14) dum2 = dum2*root5
            if (j.eq.20) dum2 = dum2*root3
            if (j.eq.24) dum2 = dum2*root7
            if (j.eq.30) dum2 = dum2*root53
            if (j.eq.33) dum2 = dum2*root3
            n = n + 1
            dij(n) = dum2
 20      continue
 30   continue
      n = 0
      dum1 = one
      max = maxl
      do 50 k = mink , maxk
         if (k.eq.8)  dum1 = root3
         if (k.eq.14) dum1 = root5
         if (k.eq.20) dum1 = dum1*root3
         if (k.eq.24) dum1 = root7
         if (k.eq.30) dum1 = dum1*root53
         if (k.eq.33) dum1 = dum1*root3
         dum2 = dum1
         if (okandl) max = k
         do 40 l = minl , max
            if (l.eq.8)  dum2 = dum2*root3
            if (l.eq.14) dum2 = dum2*root5
            if (l.eq.20) dum2 = dum2*root3
            if (l.eq.24) dum2 = dum2*root7
            if (l.eq.30) dum2 = dum2*root53
            if (l.eq.33) dum2 = dum2*root3
            n = n + 1
            dkl(n) = dum2
 40      continue
 50   continue
      nn = 0
      nnn = 0
      do 80 i = 1 , ijd
         d1 = dij(i)
         do 70 k = 1 , kld
            nn = nn + 1
            if (.not.(noform(nn))) then
               nnn = nnn + 1
               ioff = 0
               do 60 npp = 1 , npass*3
                  qq(nnn+ioff) = qq(nnn+ioff)*d1*dkl(k)
                  ioff = ioff + lendd
 60            continue
            end if
 70      continue
 80   continue
      return
      end
c
c
c  generate  r, x1, y1, z1, dd in d2escr
c
c
      subroutine dprim_dft

      implicit none

INCLUDE(../m4/common/sizes)  ! for mxprms etc

INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_dshlnf)
INCLUDE(common/dft_dshlno)
INCLUDE(common/dft_d2escr)
c
      logical odbg
      common/dbgdbg/odbg


      REAL atmp, aj, ai, ainv, dum, azi, dum1, dum2
      REAL csi, cpi, cdi, cfi, cgi
      REAL csj, cpj, cdj, cfj, cgj
      integer n, i, j, jb, max
      REAL one
      REAL axi, ayi, arri
      integer nm, nn, ia, jbmax

      data one/1.0d0/
      max = maxj
      n = 0
      nn = 0
      do 70 i = mini , maxi
         go to (20,20,30,30,
     +          20,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          20,30,30,30,30,30,30,30,30,30,
     +          30,30,30,30,30) , i
 20      nm = nn
 30      nn = nm
         if (oiandj) max = i
         do 60 j = minj , max
            go to (40,40,50,50,
     +             40,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             40,50,50,50,50,50,50,50,50,50,
     +             50,50,50,50,50) , j
 40         nn = nn + 1
 50         n = n + 1
            ijden(n) = nn
 60      continue
 70   continue
c     ----- i primitive
      nij = 0
      jbmax = ngb
      do 230 ia = 1 , nga
         ai = ag(ia)
         arri = ai*rri
         axi = ai*xi
         ayi = ai*yi
         azi = ai*zi
         csi = csa(ia)
         cpi = cpa(ia)
         cdi = cda(ia)
         cfi = cfa(ia) 
         cgi = cga(ia) 
c     ----- j primitive
         if (oiandj) jbmax = ia
         do 220 jb = 1 , jbmax
            aj = bg(jb)
            atmp = ai + aj
            ainv = one/atmp
            dum = aj*arri*ainv
            csj = csb(jb)*ainv
            cpj = cpb(jb)*ainv
            cdj = cdb(jb)*ainv
            cfj = cfb(jb)*ainv
            cgj = cgb(jb)*ainv
            nm = (nij+nij) + (nij+nij)
            nn = nm
            nij = nij + 1
            r(nij) = dum
            aa(nij) = atmp
            x1(nij) = (axi+aj*xj)*ainv
            y1(nij) = (ayi+aj*yj)*ainv
            z1(nij) = (azi+aj*zj)*ainv
c     ----- density factor
            do 190 i = mini , maxi
               if (oiandj) max = i
               go to (80,90,190,190,
     +               100,190,190,190,190,190,
     +               110,190,190,190,190,190,190,190,190,190,
     +               115,190,190,190,190,190,190,190,190,190,
     +               190,190,190,190,190) , i
 80            dum1 = csi
               go to 120
 90            dum1 = cpi
               go to 120
 100           dum1 = cdi
               go to 120
 110           dum1 = cfi
               go to 120
 115           dum1 = cgi
 120           do 180 j = minj , max
                  go to (130,140,180,180,
     +                   150,180,180,180,180,180,
     +                   160,180,180,180,180,180,180,180,180,180,
     +                   165,180,180,180,180,180,180,180,180,180,
     +                   180,180,180,180,180) , j
 130              dum2 = dum1*csj
                  go to 170
 140              dum2 = dum1*cpj
                  go to 170
 150              dum2 = dum1*cdj
                  go to 170
 160              dum2 = dum1*cfj
                  go to 170
 165              dum2 = dum1*cgj
 170              nn = nn + 1
                  dd(nn) = dum2

c            if(odbg)write(6,*)'dprim dd',nn,dd(nn)

 180           continue
 190        continue


            if (.not.oiandj) go to 220
            if (ia.eq.jb) go to 220
            go to (210,200,210,210,210) , lit
 200        if (mini.ne.2) then
               dd(nm+2) = dd(nm+2) + csi*cpj
               dd(nm+3) = dd(nm+3) + dd(nm+3)
            end if
 210        dd(nm+1) = dd(nm+1) + dd(nm+1)
 220     continue
 230  continue

c      if(odbg)write(6,*)'dprim r',nij,(r(n),n=1,nij)
c      if(odbg)write(6,*)'dprim x1',(x1(n),n=1,nij)
c      if(odbg)write(6,*)'dprim y1',(y1(n),n=1,nij)
c      if(odbg)write(6,*)'dprim z1',(z1(n),n=1,nij)


      if (nij.eq.0) return
      rsmall = r(1)
      do 240 n = 1 , nij
         exij(n) = dexp(-r(n))
 240  continue
      do 250 n = 1 , nij
         if (rsmall.gt.r(n)) rsmall = r(n)
 250  continue
      if (rsmall.ge.tol1) nij = 0
      return
      end
      subroutine dshell_dft(nelec,ish,jsh,ksh,lsh,
     +     basi, basj, bask, basl, ncentr)

c      implicit REAL  (a-h,o-z)
      implicit none
INCLUDE(../m4/common/sizes)

INCLUDE(common/dft_dmisc)
INCLUDE(../m4/common/infoa)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_root)
INCLUDE(common/dft_dshlno)

INCLUDE(common/dft_d2escr)

INCLUDE(common/dft_incrd)
INCLUDE(common/dft_dshlnf)

      logical odbg
      common/dbgdbg/odbg

c
c arguments
c
      integer nelec
      integer ish, jsh, ksh, lsh
      integer basi, basj, bask, basl, ncentr
c
c local
c
      integer i, j, k, l
      integer i1, j1, k1, l1
      integer i2, j2, k2, l2
      integer ittt, max

      if (nelec.eq.2) then

         okandl = ksh.eq.lsh .and. ncentr .eq. 4
         osame = ish.eq.ksh .and. jsh.eq.lsh .and. ncentr .eq. 4

         k = katom(bask,ksh)
         xk = c(1,k)
         yk = c(2,k)
         zk = c(3,k)
         k1 = kstart(bask,ksh)
         k2 = k1 + kng(bask,ksh) - 1
         lkt = ktype(bask,ksh)
         mink = kmin(bask,ksh)
         maxk = kmax(bask,ksh)
         lock = kloc(bask,ksh) - mink
         ngc = 0
         if(odbg)write(6,*)'dshell k loop',k,lock
         do 20 k = k1 , k2
            ngc = ngc + 1
            cgg(ngc) = ex_m(bask,k)
            csc(ngc) = cs(bask,k)
            cpc(ngc) = cp(bask,k)
            cdc(ngc) = cd(bask,k)
            cfc(ngc) = cf(bask,k)
            cgc(ngc) = cg(bask,k)
            if(odbg)write(6,*)ngc,cgg(ngc),csc(ngc),cpc(ngc),cdc(ngc)
 20      continue

         if(basl .lt. 0) then
c
c  Dummy l centr
c
            l = k
            xl = c(1,l)
            yl = c(2,l)
            zl = c(3,l)
            l1 = 1
            l2 = 1
            llt = 1
            minl = 1
            maxl = 1
            locl = 0
            ngd = 0
            if(odbg)write(6,*)'dshell l loop',l,locl
            do 30 l = l1 , l2
               ngd = ngd + 1
               dg(ngd)  = 0.0d0
               csd(ngd) = 1.0d0
               cpd(ngd) = 1.0d0
               cdd(ngd) = 1.0d0
               cfd(ngd) = 1.0d0
               cgd(ngd) = 1.0d0
               if(odbg)write(6,*)ngc,dg(ngc),csd(ngc),cpd(ngc),cdd(ngc)
 30         continue
            rrk = 0.0d0
         else
            l = katom(basl,lsh)
            xl = c(1,l)
            yl = c(2,l)
            zl = c(3,l)
            l1 = kstart(basl,lsh)
            l2 = l1 + kng(basl,lsh) - 1
            llt = ktype(basl,lsh)
            minl = kmin(basl,lsh)
            maxl = kmax(basl,lsh)
            locl = kloc(basl,lsh) - minl
            ngd = 0
            if(odbg)write(6,*)'dshell l loop',locl
            do l = l1 , l2
               ngd = ngd + 1
               dg(ngd) = ex_m(basl,l)
               csd(ngd) = cs(basl,l)
               cpd(ngd) = cp(basl,l)
               cdd(ngd) = cd(basl,l)
               cfd(ngd) = cf(basl,l)
               cgd(ngd) = cg(basl,l)
               if(odbg)write(6,*)ngd,dg(ngd),csd(ngd),cpd(ngd),cdd(ngd)
            enddo
            rrk = ((xk-xl)**2+(yk-yl)**2+(zk-zl)**2)
         endif
         nroots = (lit+ljt+lkt+llt-1)/2

c
c     determine various offsets and indexing arrays
c
         inc2 = 1
         inc3 = inc2*(maxl-minl+1)
         inc4 = inc3*(maxk-mink+1)
         inc5 = inc4*(maxj-minj+1)
         lendd = inc5*(maxi-mini+1)

c         write(6,*)'test1',nelec,inc1,inc2,inc3,inc4,inc5,lendd
c         write(6,*)mini,maxi,minj,maxj

         if (mod(lendd,4).eq.0) lendd = lendd + 1
         ijd = 0
         max = maxj
         do 50 i = mini , maxi
            if (oiandj) max = i
            ittt = inc5*(i-mini) + 1
            do 40 j = minj , max
               ijd = ijd + 1
               ijgt(ijd) = ittt
               ittt = ittt + inc4
 40         continue
 50      continue
         kld = 0
         max = maxl
         do 70 k = mink , maxk
            if (okandl) max = k
            ittt = inc3*(k-mink)
            do 60 l = minl , max
               kld = kld + 1
               klgt(kld) = ittt
               ittt = ittt + inc2
 60         continue
 70      continue
         ijkld = ijd*kld
         ixi = lendd + 1
         iyi = ixi + lendd
         izi = iyi + lendd
         ioff = 0

c         if(odbg)write(6,*)ijkld,ixi,iyi,izi
c         if(odbg)write(6,*)ijgt
c         if(odbg)write(6,*)klgt

         return
      else
         oiandj = ish.eq.jsh .and. ncentr .gt. 2
         i = katom(basi,ish)
         xi = c(1,i)
         yi = c(2,i)
         zi = c(3,i)
         i1 = kstart(basi,ish)
         i2 = i1 + kng(basi,ish) - 1
         lit = ktype(basi,ish)
         mini = kmin(basi,ish)
         maxi = kmax(basi,ish)
         loci = kloc(basi,ish) - mini
         nga = 0

         if(odbg)write(6,*)'dshell i loop',i, loci

         do 80 i = i1 , i2
            nga = nga + 1
            ag(nga) = ex_m(basi,i)
            csa(nga) = cs(basi,i)
            cpa(nga) = cp(basi,i)
            cda(nga) = cd(basi,i)
            cfa(nga) = cf(basi,i)
            cga(nga) = cg(basi,i)

            if(odbg)write(6,*)nga,ag(nga),csa(nga),cpa(nga),cda(nga)

 80      continue

         if(basj .lt. 0)then
c
c Dummy j centre
c
            j = i
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = 1
            j2 = 1
            ljt = 1
            minj = 1
            maxj = 1
            locj = 0
            ngb = 0

            if(odbg)write(6,*)'dshell j loop',j, locj

            do 90 j = j1 , j2
               ngb = ngb + 1
               bg(ngb)  = 0.0d0
               csb(ngb) = 1.0d0
               cpb(ngb) = 1.0d0
               cdb(ngb) = 1.0d0
               cfb(ngb) = 1.0d0
               cgb(ngb) = 1.0d0

               if(odbg)write(6,*)ngb,bg(ngb),csb(ngb),cpb(ngb),cdb(ngb)

 90         continue
            rri = 0.0d0
         else
            j = katom(basj,jsh)
            xj = c(1,j)
            yj = c(2,j)
            zj = c(3,j)
            j1 = kstart(basj,jsh)
            j2 = j1 + kng(basj,jsh) - 1
            ljt = ktype(basj,jsh)
            minj = kmin(basj,jsh)
            maxj = kmax(basj,jsh)
            locj = kloc(basj,jsh) - minj
            ngb = 0
            if(odbg)write(6,*)'dshell j loop',j, locj

            do j = j1 , j2
               ngb = ngb + 1
               bg(ngb)  = ex_m(basj,j)
               csb(ngb) = cs(basj,j)
               cpb(ngb) = cp(basj,j)
               cdb(ngb) = cd(basj,j)
               cfb(ngb) = cf(basj,j)
               cgb(ngb) = cg(basj,j)
               if(odbg)write(6,*)ngb,bg(ngb),csb(ngb),cpb(ngb),cdb(ngb)
            enddo
            rri = ((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
         endif
         return
      end if
      end
      subroutine dxyz_dft(x,y,z,ncdim)
c      implicit REAL  (a-h,o-z)
      implicit none

      integer ncmax
      parameter (ncmax=65)
      REAL x, y, z
      dimension x(*),y(*),z(*)
      integer ncdim

      logical n0,n1,m0,m1

      REAL bp01,b00,b10,xcp00, xc00,ycp00,yc00,zcp00,zc00,f00,
     &     dxij,dyij,dzij,dxkl,dykl,dzkl
      integer iorg,korg,nimax,njmax,nkmax,nlmax,nmax,mmax,
     &     ij1x,ij2x,kl1x,kl2x

      common/setd/bp01(ncmax),b00(ncmax),b10(ncmax),xcp00(ncmax),
     1   xc00(ncmax),
     2ycp00(ncmax),yc00(ncmax),zcp00(ncmax),zc00(ncmax),f00(ncmax),dxij,
     3dyij,dzij,dxkl,dykl,dzkl,iorg(12),korg(12),
     4nimax,njmax,nkmax,nlmax,nmax,mmax,ij1x,ij2x,kl1x,kl2x

      REAL ca, cb
      common/small/ca(ncmax),cb(ncmax)

INCLUDE(common/dft_dshlno)

      integer ni, nj, nk, nl
      integer i1, i2, i3, i4, i5, k2, k4
      integer ia, ib, ic, ink, nc, ij2
      integer i, ij1, k, n, m, km, k3, min
      integer kl1, kl2
      REAL zero, one
c
      dimension i(12),k(12)
c
      data zero,one /0.0d+00,1.0d+00/
c
      do 20 n = 1 , nmax + 1
         i(n) = (iorg(n)-1)*ncdim + 1
 20   continue
      do 30 n = 1 , mmax + 1
         k(n) = korg(n)*ncdim
 30   continue
      ij1 = ij1x*ncdim
      ij2 = ij2x*ncdim
      kl1 = kl1x*ncdim
      kl2 = kl2x*ncdim
      ink = 1
c
      n0 = nmax.eq.0
      n1 = nmax.le.1
      m0 = mmax.eq.0
      m1 = mmax.le.1
      if (n0) then
         i1 = i(1)
         ia = 0
         do 40 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 40      continue
         if (m0) go to 670
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
         do 50 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 50      continue
         if (m1) go to 670
c     ----- i(0,nk) -----
c
         do 60 nc = 1 , ncontr
            ca(nc) = zero
 60      continue
         i3 = i1
         i4 = i1 + k2
         do 90 nk = 2 , mmax
c
            do 70 nc = 1 , ncontr
               ca(nc) = ca(nc) + bp01(nc)
 70         continue
            i5 = i1 + k(nk+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 80 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
               ia = ia + ink
 80         continue
            i3 = i4
            i4 = i5
 90      continue
         if (nlmax.eq.0) go to 670
c     ----- i(0,0,nk,nl) -----
         i5 = k(mmax+1)
         min = nkmax
         go to 610
      else if (m0) then
         i1 = i(1)
         ia = 0
         do 100 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 100     continue
         if (n0) go to 600
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
         do 110 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 110     continue
         if (n1) go to 600
c     ----- i(ni,0) -----
c
         do 120 nc = 1 , ncontr
            ca(nc) = zero
 120     continue
         i3 = i1
         i4 = i2
         do 150 ni = 2 , nmax
            do 130 nc = 1 , ncontr
               ca(nc) = ca(nc) + b10(nc)
 130        continue
            i5 = i(ni+1)
            ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 140 nc = 1 , ncontr
               x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
               y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
               z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
               ia = ia + ink
 140        continue
            i3 = i4
            i4 = i5
 150     continue
         if (njmax.eq.0) go to 600
c     ----- i(ni,nj,0,0) -----
         i5 = i(nmax+1)
         min = nimax
         go to 540
      else
c     ----- i(0,0) -----
         i1 = i(1)
         ia = 0
c
         do 160 nc = 1 , ncontr
            x(i1+ia) = one
            y(i1+ia) = one
            z(i1+ia) = f00(nc)
            ia = ia + ink
 160     continue
c     ----- i(1,0) -----
         i2 = i(2)
         ia = 0
c
         do 170 nc = 1 , ncontr
            x(i2+ia) = xc00(nc)
            y(i2+ia) = yc00(nc)
            z(i2+ia) = zc00(nc)*f00(nc)
            ia = ia + ink
 170     continue
c     ----- i(0,1) -----
         k2 = k(2)
         i3 = i1 + k2
         ia = 0
c
         do 180 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)
            y(i3+ia) = ycp00(nc)
            z(i3+ia) = zcp00(nc)*f00(nc)
            ia = ia + ink
 180     continue
c     ----- i(1,1) -----
         i3 = i2 + k2
         ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
         do 190 nc = 1 , ncontr
            x(i3+ia) = xcp00(nc)*x(i2+ia) + b00(nc)
            y(i3+ia) = ycp00(nc)*y(i2+ia) + b00(nc)
            z(i3+ia) = zcp00(nc)*z(i2+ia) + b00(nc)*f00(nc)
            ia = ia + ink
 190     continue
         if (.not.(n1)) then
            do 200 nc = 1 , ncontr
               cb(nc) = b00(nc)
               ca(nc) = zero
 200        continue
            i3 = i1
            i4 = i2
            do 250 n = 2 , nmax
c     ----- i(n,0) -----
               i5 = i(n+1)
               do 210 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b10(nc)
 210           continue
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 220 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xc00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + yc00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zc00(nc)*z(i4+ia)
                  ia = ia + ink
 220           continue
               do 230 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 230           continue
c     ----- i(n,1) -----
               i3 = i5 + k2
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 240 nc = 1 , ncontr
                  x(i3+ia) = xcp00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = ycp00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zcp00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 240           continue
               i3 = i4
               i4 = i5
 250        continue
         end if
         if (.not.(m1)) then
            do 260 nc = 1 , ncontr
               ca(nc) = zero
               cb(nc) = b00(nc)
 260        continue
            i3 = i1
            i4 = i1 + k2
            do 310 m = 2 , mmax
               do 270 nc = 1 , ncontr
                  ca(nc) = ca(nc) + bp01(nc)
 270           continue
c     ----- i(0,m) -----
               i5 = i1 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 280 nc = 1 , ncontr
                  x(i5+ia) = ca(nc)*x(i3+ia) + xcp00(nc)*x(i4+ia)
                  y(i5+ia) = ca(nc)*y(i3+ia) + ycp00(nc)*y(i4+ia)
                  z(i5+ia) = ca(nc)*z(i3+ia) + zcp00(nc)*z(i4+ia)
                  ia = ia + ink
 280           continue
               do 290 nc = 1 , ncontr
                  cb(nc) = cb(nc) + b00(nc)
 290           continue
c     ----- i(1,m) -----
               i3 = i2 + k(m+1)
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 300 nc = 1 , ncontr
                  x(i3+ia) = xc00(nc)*x(i5+ia) + cb(nc)*x(i4+ia)
                  y(i3+ia) = yc00(nc)*y(i5+ia) + cb(nc)*y(i4+ia)
                  z(i3+ia) = zc00(nc)*z(i5+ia) + cb(nc)*z(i4+ia)
                  ia = ia + ink
 300           continue
               i3 = i4
               i4 = i5
 310        continue
         end if
         if (.not.(n1 .or. m1)) then
c     ----- i(n,m) -----
c
            do 320 nc = 1 , ncontr
               ca(nc) = b00(nc)
 320        continue
            k3 = k2
            do 370 m = 2 , mmax
               k4 = k(m+1)
               do 330 nc = 1 , ncontr
                  ca(nc) = ca(nc) + b00(nc)
                  cb(nc) = b10(nc)
 330           continue
               i3 = i1
               i4 = i2
               do 360 n = 2 , nmax
                  i5 = i(n+1)
                  ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 340 nc = 1 , ncontr
                     x(i5+k4+ia) = cb(nc)*x(i3+k4+ia) + xc00(nc)
     +                             *x(i4+k4+ia) + ca(nc)*x(i4+k3+ia)
                     y(i5+k4+ia) = cb(nc)*y(i3+k4+ia) + yc00(nc)
     +                             *y(i4+k4+ia) + ca(nc)*y(i4+k3+ia)
                     z(i5+k4+ia) = cb(nc)*z(i3+k4+ia) + zc00(nc)
     +                             *z(i4+k4+ia) + ca(nc)*z(i4+k3+ia)
                     ia = ia + ink
 340              continue
                  do 350 nc = 1 , ncontr
                     cb(nc) = cb(nc) + b10(nc)
 350              continue
                  i3 = i4
                  i4 = i5
 360           continue
               k3 = k4
 370        continue
         end if
         if (njmax.eq.0) go to 450
c     ----- i(ni,nj,m) -----
         m = 0
         i5 = i(nmax+1)
      end if
 380  min = nimax
      km = k(m+1)
 390  n = nmax
      i3 = i5 + km
 400  i4 = i(n) + km
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 410 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 410  continue
      i3 = i4
      n = n - 1
      if (n.gt.min) go to 400
      min = min + 1
      if (min.lt.nmax) go to 390
      if (nimax.ne.0) then
         i3 = ij2 + km + i1
         do 440 nj = 1 , njmax
            i4 = i3
            do 430 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 420 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 420           continue
               i4 = i4 + ij1
 430        continue
            i3 = i3 + ij2
 440     continue
      end if
      m = m + 1
      if (m.le.mmax) go to 380
 450  if (nlmax.eq.0) go to 530
c
c     ----- i(ni,nj,nk,nl) -----
c
      i5 = k(mmax+1)
      ia = i1
      ni = 0
 460  nj = 0
      ib = ia
      min = nkmax
 470  m = mmax
      i3 = ib + i5
 480  i4 = ib + k(m)
      ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 490 nc = 1 , ncontr
         x(i3+ic) = x(i3+ic) + dxkl*x(i4+ic)
         y(i3+ic) = y(i3+ic) + dykl*y(i4+ic)
         z(i3+ic) = z(i3+ic) + dzkl*z(i4+ic)
         ic = ic + ink
 490  continue
      i3 = i4
      m = m - 1
      if (m.gt.min) go to 480
      min = min + 1
      if (min.lt.mmax) go to 470
      if (nkmax.ne.0) then
         i3 = ib + kl2
         do 520 nl = 1 , nlmax
            i4 = i3
            do 510 nk = 1 , nkmax
               ic = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 500 nc = 1 , ncontr
                  x(i4+ic) = x(i4+ic+kl1-kl2) + dxkl*x(i4+ic-kl2)
                  y(i4+ic) = y(i4+ic+kl1-kl2) + dykl*y(i4+ic-kl2)
                  z(i4+ic) = z(i4+ic+kl1-kl2) + dzkl*z(i4+ic-kl2)
                  ic = ic + ink
 500           continue
               i4 = i4 + kl1
 510        continue
            i3 = i3 + kl2
 520     continue
      end if
      nj = nj + 1
      ib = ib + ij2
      if (nj.le.njmax) then
         min = nkmax
         go to 470
      else
         ni = ni + 1
         ia = ia + ij1
         if (ni.le.nimax) go to 460
      end if
 530  return
 540  ni = nmax
      i3 = i5
 550  i4 = i(ni)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 560 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxij*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dyij*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzij*z(i4+ia)
         ia = ia + ink
 560  continue
      i3 = i4
      ni = ni - 1
      if (ni.gt.min) go to 550
      min = min + 1
      if (min.lt.nmax) go to 540
      if (nimax.ne.0) then
         i3 = ij2 + i1
         do 590 nj = 1 , njmax
            i4 = i3
            do 580 ni = 1 , nimax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 570 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+ij1-ij2) + dxij*x(i4+ia-ij2)
                  y(i4+ia) = y(i4+ia+ij1-ij2) + dyij*y(i4+ia-ij2)
                  z(i4+ia) = z(i4+ia+ij1-ij2) + dzij*z(i4+ia-ij2)
                  ia = ia + ink
 570           continue
               i4 = i4 + ij1
 580        continue
            i3 = i3 + ij2
 590     continue
      end if
 600  return
 610  nk = mmax
      i3 = i1 + i5
 620  i4 = i1 + k(nk)
      ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
      do 630 nc = 1 , ncontr
         x(i3+ia) = x(i3+ia) + dxkl*x(i4+ia)
         y(i3+ia) = y(i3+ia) + dykl*y(i4+ia)
         z(i3+ia) = z(i3+ia) + dzkl*z(i4+ia)
         ia = ia + ink
 630  continue
      i3 = i4
      nk = nk - 1
      if (nk.gt.min) go to 620
      min = min + 1
      if (min.lt.mmax) go to 610
      if (nkmax.ne.0) then
         i3 = i1 + kl2
         do 660 nl = 1 , nlmax
            i4 = i3
            do 650 nk = 1 , nkmax
               ia = 0
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
               do 640 nc = 1 , ncontr
                  x(i4+ia) = x(i4+ia+kl1-kl2) + dxkl*x(i4+ia-kl2)
                  y(i4+ia) = y(i4+ia+kl1-kl2) + dykl*y(i4+ia-kl2)
                  z(i4+ia) = z(i4+ia+kl1-kl2) + dzkl*z(i4+ia-kl2)
                  ia = ia + ink
 640           continue
               i4 = i4 + kl1
 650        continue
            i3 = i3 + kl2
 660     continue
      end if
 670  return
      end
      subroutine subsd_dft(x,y,z,xd,yd,zd,
     *a,m1,m2,m3,m4,i1,i2,k1,k2,ncdim)

c      implicit REAL  (a-h,o-z)
      implicit none

      integer ncmax
      parameter (ncmax=65)
      REAL x,y,z,xd,yd,zd,a
      integer ncdim
      dimension a(ncmax)
      dimension x(ncdim,*),y(ncdim,*),z(ncdim,*),
     *         xd(ncdim,*),yd(ncdim,*),zd(ncdim,*)
      integer m1,m2,m3,m4,i1,i2,k1,k2

INCLUDE(common/dft_incrd)
INCLUDE(common/dft_dshlno)
c
      logical odbg
      common/dbgdbg/odbg

      integer i, j, k, l
      integer n1, n2, n3, n4, nr
      REAL fac

      do 20 i = 1 , ncontr
         a(i) = a(i) + a(i)
 20   continue
      n1 = 1
c
      do 120 i = 1 , m1
         n2 = n1
         do 110 j = 1 , m2
            n3 = n2
            do 100 k = 1 , m3
               n4 = n3
               do 90 l = 1 , m4
                  go to (30,50,70,70,70,70,70) , l
c
 30               do 40 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2)

                     
       if(odbg)write(6,*)'subsd1',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)


 40               continue
                  n4 = n4 + k2
                  go to 90
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 50               do 60 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) - x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) - y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) - z(nr,n4-k2)

       if(odbg)write(6,*)'subsd2',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)

 60               continue
                  n4 = n4 + k2
                  go to 90
 70               fac = -dfloat(l-1)
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
                  do 80 nr = 1 , ncontr
                     xd(nr,n4) = a(nr)*x(nr,n4+k2) + fac*x(nr,n4-k2)
                     yd(nr,n4) = a(nr)*y(nr,n4+k2) + fac*y(nr,n4-k2)
                     zd(nr,n4) = a(nr)*z(nr,n4+k2) + fac*z(nr,n4-k2)

       if(odbg)write(6,*)'subsd3',nr,n4,xd(nr,n4),yd(nr,n4),zd(nr,n4)

 80               continue
c
                  n4 = n4 + k2
 90            continue
               n3 = n3 + k1
 100        continue
            n2 = n2 + i2
 110     continue
         n1 = n1 + i1
 120  continue
      return
      end
      subroutine redund_dft(ii,jj,kk,ll,
     +     basi, basj, bask, basl, 
     +     iw)

c      implicit REAL  (a-h,o-z)
      implicit none

      integer ii,jj,kk,ll,iw
      integer basi, basj, bask, basl

INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_mbasis)
c
      integer iat, jat, kat, lat
      integer i,lll,lit,ljt,llt,lkt,iper
      integer min, imin
      integer n1, n2
      logical inej,inek,inel,jnek,jnel,knel

      dimension lll(4)
      equivalence (lll(1),lit)
      oskip(1) = .true.
      oskip(2) = .true.
      oskip(3) = .true.
      oskip(4) = .true.
      npass = 0
      do 20 i = 1 , 4
         natomd(i) = 0
 20   continue

      lit = ktype(basi,ii)
      iat = katom(basi,ii)
      if(basj.lt.0)then
         ljt = 1
         jat = iat
      else
         ljt = ktype(basj,jj)
         jat = katom(basj,jj)
      endif

      lkt = ktype(bask,kk)
      kat = katom(bask,kk)
      if(basl.lt.0)then
         llt = 1
         lat = kat
      else
         llt = ktype(basl,ll)
         lat = katom(basl,ll)
      endif

      inej = iat.ne.jat
      inek = iat.ne.kat
      inel = iat.ne.lat
      jnek = jat.ne.kat
      jnel = jat.ne.lat
      knel = kat.ne.lat

      if (inej) then
         if (.not.(inek)) then
            if (.not.(inel)) go to 40
c      iat=kat    jat=lat
            if (jnel) go to 50
            if (ii.ne.kk .or. jj.ne.ll) then
               n1 = (lit+1)*(lkt+1)*ljt*llt
               n2 = lit*lkt*(ljt+1)*(llt+1)
               if (n1.ge.n2) go to 50
               go to 70
            else
               if (ljt.le.lit) go to 60
               go to 40
            end if
         else if (jnek) then
            if (.not.(jnel)) go to 70
            if (.not.(knel)) go to 80
c     iat # jat # kat # lat  -- omit one centre
            min = lit
            imin = 1
            do 30 iper = 2 , 4
               if (lll(iper).lt.min) then
                  min = lll(iper)
                  imin = iper
               end if
 30         continue
            go to (90,100,110,120) , imin
            go to 90
         else
            if (.not.(jnel)) go to 60
c     ----- jat = kat # iat # lat -----
            oskip(1) = .false.
            oskip(4) = .false.
            natomd(1) = iat
            natomd(2) = lat
            natomd(3) = jat
            npass = 2
            go to 130
         end if
      else if (inek) then
         if (.not.(knel)) then
c     iat=jat  ,  kat=lat   differentiate one pair
            n1 = lit*ljt*(lkt+1)*(llt+1)
            n2 = (lit+1)*(ljt+1)*lkt*llt
            if (n2.lt.n1) go to 80
         end if
      else
         if (inel) then
c     iat=jat=kat   derivative (ij/kl')
            oskip(4) = .false.
            natomd(1) = lat
            natomd(2) = iat
            npass = 1
         end if
         go to 130
      end if
c     iat=jat   derivatives (ij/k'l) and (ij/kl')
      oskip(3) = .false.
      oskip(4) = .false.
      natomd(1) = kat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     iat=kat=lat   derivative (ij'/kl)
 40   oskip(2) = .false.
      natomd(1) = jat
      natomd(2) = iat
      npass = 1
      go to 130
c     iat=kat   derivatives (ij'/kl) and (ij/kl')
 50   oskip(2) = .false.
      oskip(4) = .false.
      natomd(1) = jat
      natomd(2) = lat
      natomd(3) = iat
      npass = 2
      go to 130
c     jat=kat=lat    (i'j/kl)
 60   oskip(1) = .false.
      natomd(1) = iat
      natomd(2) = jat
      npass = 1
      go to 130
c      jat=lat    derivatives (i'j/kl) and (ij/k'l)
 70   oskip(1) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = kat
      natomd(3) = jat
      npass = 2
      go to 130
c     kat=lat   derivatives (i'j/kl) and (ij'/kl)
 80   oskip(1) = .false.
      oskip(2) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      npass = 2
      go to 130
 90   natomd(1) = jat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = iat
      npass = 3
      oskip(2) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 100  natomd(1) = iat
      natomd(2) = kat
      natomd(3) = lat
      natomd(4) = jat
      npass = 3
      oskip(1) = .false.
      oskip(3) = .false.
      oskip(4) = .false.
      go to 130
 110  natomd(1) = iat
      natomd(2) = jat
      natomd(3) = lat
      natomd(4) = kat
      npass = 3
      oskip(1) = .false.
      oskip(2) = .false.
      oskip(4) = .false.
      go to 130
 120  oskip(1) = .false.
      oskip(2) = .false.
      oskip(3) = .false.
      natomd(1) = iat
      natomd(2) = jat
      natomd(3) = kat
      natomd(4) = lat
      npass = 3
c     -----
 130  if (.not.outd) return
      write (iw,6010) ii , jj , kk , ll , 
     +    oskip(1) , oskip(2) , oskip(3) , oskip(4) , 
     +    npass , (natomd(i),i=1,4)
      return
 6010 format (/,' ********  ii,jj,kk,ll =',4i3,' skip1,2,3,4 =',4l3,
     +        ' npass =',i2,' centres =',4i5,/)
      end

      subroutine grdcon_dft(qq,noform,gmax)
c      implicit REAL  (a-h,o-z)
      implicit none
c
c arguments
c
      REAL qq
      logical noform
      dimension noform(*),qq(*)

      REAL gmax


INCLUDE(../m4/common/sizes)

ccINCLUDE(../m4/common/specal)
  
      REAL dgout
      common/tgrad/dgout(9)

INCLUDE(common/dft_d2escr)

INCLUDE(common/dft_dshlno)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_dmisc)

c-      REAL dipd,dipn,dipi
c-      common/dipblk/dipd(3,3,maxat),dipn(3,3,maxat),dipi(3,3,maxat)
INCLUDE(../m4/common/restrl)

      logical odbg
      common/dbgdbg/odbg

c
c functions
c
_IF(single)
      REAL sdot
_ELSEIF(hp700)
      REAL `vec_$ddot'
_ELSE
      REAL ddot
_ENDIF
c
c local variables
c
      integer ii, jj, kk, ll
      integer i1, j1, k1, l1
      integer n0, n1, n2, n3, nn, nnn, np3, npp
      integer lmax, jmax
      REAL dumx, dumy, dumz, densty

      integer iii

c-      integer i, j, k, l
c-      logical nofk
c-      integer ioa, iob
c-      REAL t1

      REAL half
      data half/0.5d0/

c-      nofk = .not.(ofock .or. ofokab)
c
c
c-      if ((mp2w .or. ompir)) then
c-         ioa = ic7 + lendd*npass*3
c-         call vclr(qq(ioa+1),1,3*lendd)
c-         do 30 i = 1 , npass
c-            iob = ic7 + (i-1)*lendd*3
c-_IF1(ct)cdir$ ivdep
c-_IF1(a)cvd$  nodepck
c-_IF1(x)c$dir no_recurrence
c-            do 20 k = 1 , lendd*3
c-               qq(ioa+k) = qq(ioa+k) - qq(iob+k)
c- 20         continue
c- 30      continue
c-      end if
c-      if (ompir) then
c-         n2 = iabd + lendd + 1
c-         do 60 i = 1 , 3
c-            n1 = ic7 + 1
c-            do 50 k = 1 , npass + 1
c-               do 40 j = 1 , 3
c-                  t1 = ddot(ijkld,qq(n2),1,qq(n1),1)
c-                  dipi(i,j,natomd(k)) = dipi(i,j,natomd(k)) - t1
c-                  n1 = n1 + lendd
c- 40            continue
c- 50         continue
c-            n2 = n2 + lendd
c- 60      continue
c-      end if
c
c-      if (outd .or. ofock .or. ofokab) then
        if (outd) then
c
c     used only if extra output required or if need
c     derivatives of fock matrices
c
         jmax = maxj
         nn = 0
         nnn = 0
         do 110 ii = mini , maxi
            i1 = loci + ii
            if (oiandj) jmax = ii
            do 100 jj = minj , jmax
               j1 = locj + jj
c
c
               lmax = maxl
               do 90 kk = mink , maxk
                  k1 = lock + kk
                  if (okandl) lmax = kk
                  do 80 ll = minl , lmax
                     nn = nn + 1
                     if (.not.(noform(nn))) then
                        nnn = nnn + 1
                        l1 = locl + ll
                        densty = qq(iabd+nnn)
                        n0 = 1
                        n1 = ic7
                        do 70 npp = 1 , npass
                           n2 = n1 + lendd
                           n3 = n2 + lendd
                           dumx = qq(n1+nnn)
                           dumy = qq(n2+nnn)
                           dumz = qq(n3+nnn)
                           dgout(n0) = dgout(n0) + densty*dumx
                           dgout(n0+1) = dgout(n0+1) + densty*dumy
                           dgout(n0+2) = dgout(n0+2) + densty*dumz
                           if (outd) write(iwr,6010) npp , i1 , j1 , 
     +                         k1 , l1 , n1 , n2 , n3 , ic7 , lendd , 
     +                         nn , nnn , dumx , dumy , dumz , densty
c-                           if (.not.(nofk)) then
c-                              if (i1.eq.j1) then
c-                                 dumx = dumx*half
c-                                 dumy = dumy*half
c-                                 dumz = dumz*half
c-                              end if
c-                              if (k1.eq.l1) then
c-                                 dumx = dumx*half
c-                                 dumy = dumy*half
c-                                 dumz = dumz*half
c-                              end if
c-                              qq(n1+nnn) = q4*dumx
c-                              qq(n2+nnn) = q4*dumy
c-                              qq(n3+nnn) = q4*dumz
c-                           end if
                           n0 = n0 + 3
                           n1 = n3 + lendd
 70                     continue
                     end if
 80               continue
 90            continue
 100        continue
 110     continue
      else
c
c gmax used for debug purposes
c
c         n1 = ic7 + 1
c         gmax = 0.0
c         do n0 = 1 , np3
c            gmax = max(gmax, abs(ddot(ijkld,qq(iabd+1),1,qq(n1),1)))
c            n1 = n1 + lendd
c         enddo
c         if(abs(gmax) .gt. 1.0d0)odbg =  .true.

         n1 = ic7 + 1
         np3 = npass*3

         if(odbg)write(6,*)'dens',(qq(iabd+iii),iii=1,ijkld)
         do 120 n0 = 1 , np3

            if(odbg)write(6,*)'integ',n0,(qq(n1+iii-1),iii=1,ijkld)

            dgout(n0) = dgout(n0) + ddot(ijkld,qq(iabd+1),1,qq(n1),1)
            if(odbg)write(6,*)'dgout',n0,
     &           ddot(ijkld,qq(iabd+1),1,qq(n1),1)

            n1 = n1 + lendd
 120     continue

      end if
 6010 format (1x,12i6/30x,3e16.8,5x,e16.8)
      end
      subroutine formeg_dft

      implicit none

INCLUDE(../m4/common/sizes)
INCLUDE(common/dft_dmisc)
INCLUDE(common/dft_grad2)


      logical odbg
      common/dbgdbg/odbg

      REAL dgout
      common/tgrad/dgout(3,3)

      REAL dum, dumx, dumy, dumz
      integer ipass, iat

      REAL zero
      data zero /0.0d0/

      dumx = zero
      dumy = zero
      dumz = zero

      if(odbg)write(6,*)'formeg',dgout

      do 20 ipass = 1 , npass
         iat = natomd(ipass)
         dum = dgout(1,ipass)
         dumx = dumx + dum
         de(1,iat) = de(1,iat) + dum
         dum = dgout(2,ipass)
         dumy = dumy + dum
         de(2,iat) = de(2,iat) + dum
         dum = dgout(3,ipass)
         dumz = dumz + dum
         de(3,iat) = de(3,iat) + dum
 20   continue
      iat = natomd(npass+1)
      de(1,iat) = de(1,iat) - dumx
      de(2,iat) = de(2,iat) - dumy
      de(3,iat) = de(3,iat) - dumz
      return
      end
c
c normalisation factors for contraction with fitting coefficients
c
      block data dft_dfgdat

INCLUDE(common/dft_dfgfac)

      data fac/7*1.0d0,
     &     3*0.5773502691896258d0,
     &     3*1.0d0,
     &     6*0.4472135954999579d0,
     &     0.2581988897471611d0,
     &     3*1.0d0,
     &     6*0.3779644730092272d0,
     &     3*0.2927700218845599d0,
     &     3*0.1690308509457033d0/
      end
c
c
      subroutine dabab_dft(ii,jj,kk,ll,
     &     basi, basj, bask, basl,
     &     ncentr,q4,
     &     zscftp,da,db,cfit,cfit2,abdens)

      implicit none
c
c arguments
c
      REAL da(*), db(*), cfit(*), cfit2(*), abdens(*)
      integer ii,jj,kk,ll,basi, basj, bask, basl	
      REAL q4
      character*8 zscftp
      integer ncentr

INCLUDE(../m4/common/sizes)

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_mapper)
INCLUDE(common/dft_incrd)
INCLUDE(common/dft_mbasis)
INCLUDE(common/dft_module_comm)

c
c local variables
c
      integer mjk, mkl, mil, mjl, mij, mik
      integer i, loci, locj, lock, locl, nn
      integer l, l1, ll1, ii1, j, jj1, j1, k1, k, kk1
      integer ni, nj, nk, nl
      integer mink, maxk, minl, maxl, mini, maxi, minj, maxj
      integer i1, i2, i3, i4
      integer n, is, js, isi, joff
      logical okleq, oijeq
      REAL dfac
_IF(64bitpointers)
      integer*8 itmp
_ELSE
      integer itmp
_ENDIF

INCLUDE(common/ccpdft.hf77)
INCLUDE(common/dft_dfgfac)

      REAL wght

      REAL pt5,four

      logical odbg
      common/dbgdbg/odbg

      data pt5,four /0.5d0,4.0d0/

c
c      ouhf = zscftp.eq.zuhf
c      ogrhf = zscftp.eq.zgrhf

c
c make conditional on ncentr
c

      ni = 1
c
c-      if (.not.ogrhf) then
c

c         call chkadr2(da(1),itmp)
c         write(6,*)'addr dens input',itmp, da(1),da(2),da(3),da(4)

      if(ncentr .eq. 4) then

         mini = kmin(basi,ii)
         minj = kmin(basj,jj)
         mink = kmin(bask,kk)
         minl = kmin(basl,ll)

         maxi = kmax(basi,ii)
         maxj = kmax(basj,jj)
         maxk = kmax(bask,kk)
         maxl = kmax(basl,ll)

         loci = kloc(basi,ii) - mini
         locj = kloc(basj,jj) - minj
         lock = kloc(bask,kk) - mink
         locl = kloc(basl,ll) - minl

         do 60 i = mini , maxi
            nj = ni
            do 50 j = minj , maxj
               nk = nj
               do 40 k = mink , maxk
                  nl = nk
                  do 30 l = minl , maxl
                     nn = nl
                     i1 = loci + i
                     i2 = locj + j
                     i3 = lock + k
                     i4 = locl + l
                     if (i1.lt.i2) then
                        n = i1
                        i1 = i2
                        i2 = n
                     end if
                     if (i3.lt.i4) then
                        n = i3
                        i3 = i4
                        i4 = n
                     end if
                     if (i1.lt.i3) then
                     else if (i1.eq.i3) then
                        if (i2.ge.i4) go to 20
                     else
                        go to 20
                     end if
                     n = i1
                     i1 = i3
                     i3 = n
                     n = i2
                     i2 = i4
                     i4 = n
 20                  mij = iky(i1) + i2
                     mik = iky(i1) + i3
                     mil = iky(i1) + i4
                     mkl = iky(i3) + i4
                     if (i2.lt.i3) then
                        mjk = iky(i3) + i2
                        if (i2.lt.i4) then
                           mjl = iky(i4) + i2
                        else
                           mjl = iky(i2) + i4
                        end if
                     else
                        mjk = iky(i2) + i3
                        mjl = iky(i2) + i4
                     end if

c-                     if(.not. CD_active())then
c-                        dfac = da(mij)*da(mkl)*four - da(mik)*da(mjl)
c-     +                       - da(mil)*da(mjk)
c-                        if (ouhf) dfac = dfac - db(mik)*db(mjl) 
c-     +                       - db(mil)*db(mjk)
c-                     else
c
c optionally include coulomb and fraction of exclude exchange contribution
c
c-                        if (ouhf) call caserr('no uhf dft yet')

                        dfac=0.0d0

                        if(CD_HF_coulomb_deriv())then
                           dfac = da(mij)*da(mkl)*four 
                           if (.not.rks_sw) then
                              dfac = dfac + db(mij)*db(mkl)*four
                           endif
                        endif

                        if(CD_HF_exchange())then
                           wght = CD_HF_exchange_weight()
                           dfac = dfac - 
     &                         wght*(da(mik)*da(mjl) + da(mil)*da(mjk))
                           if (.not.rks_sw) then
                              dfac = dfac - 
     &                         wght*(db(mik)*db(mjl) + db(mil)*db(mjk))
                           endif
                        endif
c-                     endif

                     if (i1.eq.i2) dfac = dfac*pt5
                     if (i3.eq.i4) dfac = dfac*pt5
                     dfac = dfac*q4

c-                     if (omp2 .or. mp3) then
c-                        abdens(nn) = abdens(nn) + dfac
c-                     else
                        abdens(nn) = dfac
c-                     end if

                     nl = nl + inc2
 30               continue
                  nk = nk + inc3
 40            continue
               nj = nj + inc4
 50         continue
            ni = ni + inc5
 60      continue

      else if (ncentr .eq. 3) then

         mini = kmin(basi,ii)
         minj = kmin(basj,jj)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxj = kmax(basj,jj)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         locj = kloc(basj,jj) - minj
         lock = kloc(bask,kk) - mink

         ni = 1
         do  i = mini , maxi
            nj = ni
            do j = minj , maxj
               nk = nj

               i1 = loci + i
               i2 = locj + j

               if (i1.lt.i2) then
                  n = i1
                  i1 = i2
                  i2 = n
               end if

               mij = iky(i1) + i2

               do k = mink , maxk
                  nn = nk

                  dfac = 2.0d0*da(mij)*cfit(lock + k)*fac(k)
                  if (.not.rks_sw) then
                     dfac = dfac + 2.0d0*db(mij)*cfit(lock + k)*fac(k)
                  endif
                  if (i1.eq.i2) dfac = dfac*pt5
                  dfac = dfac*q4
                  abdens(nn) = dfac

                  nk = nk + inc3
               enddo
               nj = nj + inc4
            enddo
            ni = ni + inc5
         enddo

      else if (ncentr .eq. 2) then
c
         if(basi .eq. bask)then

         mini = kmin(basi,ii)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         lock = kloc(bask,kk) - mink

         ni = 1
         do i = mini , maxi
            nk = ni
            do k = mink , maxk
               nn = nk

               if(odbg)write(6,*)'dabab',nn,loci,lock,i,k,
     &	 cfit(loci + i),cfit(lock + k),fac(i),fac(k)

               abdens(nn) = - q4 * cfit(loci + i) * cfit(lock + k) *
     +              fac(i) * fac(k)
c NB diagonal terms give no gradients
c               if(i .eq. k)abdens(nn) = abdens(nn) * 0.5d0
               nk = nk + inc3
            enddo
            ni = ni + inc5
         enddo

         else
c
c Two distinct interacting charge distributions
c
         mini = kmin(basi,ii)
         mink = kmin(bask,kk)

         maxi = kmax(basi,ii)
         maxk = kmax(bask,kk)

         loci = kloc(basi,ii) - mini
         lock = kloc(bask,kk) - mink

         ni = 1
         do i = mini , maxi
            nk = ni
            do k = mink , maxk
               nn = nk

c               if(odbg)write(6,*)'dabab',nn,loci,lock,i,k,
c     &	 cfit(loci + i),cfit(lock + k),fac(i),fac(k)

               abdens(nn) = - q4 * cfit(loci + i) * cfit2(lock + k) *
     +              fac(i) * fac(k)

               nk = nk + inc3
            enddo
            ni = ni + inc5
         enddo

         endif

      endif

c
c-      else
c-c
c-c     general case
c-c
c-         do 120 i = mini , maxi
c-            nj = ni
c-            i1 = loci + i
c-            ii1 = iky(i1)
c-            do 110 j = minj , maxj
c-               nk = nj
c-               j1 = locj + j
c-               jj1 = iky(j1)
c-               mij = ii1 + j1
c-               if (j1.gt.i1) mij = jj1 + i1
c-               oijeq = i1.eq.j1
c-               do 100 k = mink , maxk
c-                  nl = nk
c-                  k1 = lock + k
c-                  kk1 = iky(k1)
c-                  mik = ii1 + k1
c-                  if (k1.gt.i1) mik = kk1 + i1
c-                  mjk = jj1 + k1
c-                  if (k1.gt.j1) mjk = kk1 + j1
c-                  do 90 l = minl , maxl
c-                     nn = nl
c-                     l1 = locl + l
c-                     ll1 = iky(l1)
c-                     mkl = kk1 + l1
c-                     if (l1.gt.k1) mkl = ll1 + k1
c-                     mjl = jj1 + l1
c-                     if (l1.gt.j1) mjl = ll1 + j1
c-                     mil = ii1 + l1
c-                     if (l1.gt.i1) mil = ll1 + i1
c-                     okleq = k1.eq.l1
c-                     dfac = 0.0d0
c-                     ioff = 0
c-                     do 80 is = 1 , njk
c-                        isi = (is-1)*11
c-                        joff = 0
c-                        do 70 js = 1 , njk
c-                           dfac = dfac + 4.0d0*erga(isi+js)
c-     +                            *(da(ioff+mij)*da(joff+mkl)
c-     +                            +da(ioff+mkl)*da(joff+mij))
c-     +                            + 2.0d0*ergb(isi+js)
c-     +                            *(da(ioff+mik)*da(joff+mjl)
c-     +                            +da(ioff+mjl)*da(joff+mik)
c-     +                            +da(ioff+mil)*da(joff+mjk)
c-     +                            +da(ioff+mjk)*da(joff+mil))
c-                           joff = joff + nx
c- 70                     continue
c-                        ioff = ioff + nx
c- 80                  continue
c-                     if (oijeq) dfac = dfac*pt5
c-                     if (okleq) dfac = dfac*pt5
c-                     abdens(nn) = dfac*q4
c-                     nl = nl + inc2
c- 90               continue
c-                  nk = nk + inc3
c- 100           continue
c-               nj = nj + inc4
c- 110        continue
c-            ni = ni + inc5
c- 120     continue
c-c
c-      end if
c
      return
      end
      subroutine ver_dft_deriv2e(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/deriv2e.m,v $
     +     "/
      data revision /
     +     "$Revision: 6317 $"
     +      /
      data date /
     +     "$Date: 2015-03-13 21:56:17 +0100 (Fri, 13 Mar 2015) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
