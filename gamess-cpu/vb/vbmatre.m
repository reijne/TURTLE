
c  $Author: jvl $
c  $Date: 2014-09-15 17:31:14 +0200 (Mon, 15 Sep 2014) $
c  $Locker:  $
c  $Revision: 6296 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbmatre.m,v $
c  $State: Exp $
c  $Log: not supported by cvs2svn $
c  Revision 1.44  2007/09/18 10:06:55  jvl
c  VB Parallel bugs fixed for huge VB calculations.
c  HSMAT switched to igmem_alloc. Probably more fixes needed.
c  /marcin
c
c  Revision 1.43  2007/03/20 14:49:31  jvl
c  Pretty major overhaul of the VB code
c  say 50% is now dynamic using igmem_alloc and, oh wonder,
c  all the examples are still checking out OK, including a new resonating
c  super hybrid. (will follow)
c
c  Revision 1.42  2007/03/15 14:03:09  jvl
c  update super hybrid oprion in vb (ignore ao per atom)
c  added damping
c  clean up
c
c  Revision 1.41  2006/12/13 21:56:08  rwa
c  First check in of the experimental VB response property code.
c
c  Revision 1.40  2006/02/15 13:04:55  hvd
c  A comment mentioning an M4 preprocessor construct was actually logged in
c  the source file itself and subsequently processed by m4 causing all hell
c  to break loose. I have broken the construct in the comments to stop m4
c  processing it. However I am not sure how CVS is going to respond to
c  modifications in the log section of the source file. So this checkin will
c  be interesting...
c
c  Revision 1.39  2006/02/10 20:40:10  hvd
c  Both the g95 and gfortran do not support d-lines.
c  Therefore any d-lines were replaced by IF(debug) preprocessor constructs.
c
c  Revision 1.38  2006/01/13 17:56:48  jmht
c  Merging in changes from the release branch
c
c  Revision 1.37  2005/09/26 11:34:29  jmht
c  A number of files seem not to have been checked in following the merge, so I am chasing up the stragglers.
c
c  Revision 1.36  2005/06/28 11:16:04  jvl
c  corrected shift in brillouin condition more clearly
c
c  Revision 1.35  2005/06/24 10:50:43  jvl
c  Fixed a bug regarding brmax-value when a level shifter is used (JJE).
c
c  Revision 1.34.2.3  2005/07/19 06:52:20  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.34.2.2  2005/07/05 17:13:24  jmht
c  Adding Joop's change on the development branch (1.36):
c
c  corrected shift in brillouin condition more clearly.
c
c  Revision 1.34.2.1  2005/06/24 05:42:36  wab
c  Second attempt to do this after the bug in guess that caused apparent
c  failures to be diagnosed yesterday has been fixed
c
c  Revision 1.34  2005/04/22 11:35:04  jvl
c  stats_vb i8 change
c
c  Revision 1.33  2005/04/22 11:07:54  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.32  2005/03/31 15:09:10  jvl
c  nstruc as f.p. in hamiltt clashed with hsinfo and was superfluous => removed
c
c  Revision 1.31  2005/03/31 14:59:59  jvl
c  - corrected pert option, which was screwed earlier
c  - include file hsinfo
c
c  Revision 1.30  2005/03/17 16:28:22  jvl
c  fixed brmax and pert again
c
c  Revision 1.29  2005/03/17 11:44:30  jvl
c  a brmax=0.0 too many was left
c
c  Revision 1.28  2005/03/17 11:07:51  jvl
c  correct brmax and hybrid now (???
c
c  Revision 1.27  2005/02/14 09:46:42  jvl
c  brmax (brillouin theorem) now ok for parallel
c
c  Revision 1.26  2005/02/05 18:04:26  wab
c
c  Number of changes here - largely cosmetic that appeared when building
c  the code on the Alpha. Joop I'm sure none of these will effect the
c  correct running of the code.  Note that one of the examples in
c  examples/vb, c2h2, fails on the Alpha, both before and after these
c  changes. Have you seen this problem before?
c
c  segmentation violation diagnostic is as follows:
c
c  forrtl: severe (174): SIGSEGV, segmentation fault occurred
c     0: __FINI_00_remove_gp_range [0x3ff81a1fec8]
c     1: __FINI_00_remove_gp_range [0x3ff81a294a4]
c     2: __FINI_00_remove_gp_range [0x3ff800d0cac]
c     3: normt_ [vbscf.f: 4416, 0x1214bfc40]
c     4: hybmak_ [vbin.f: 6583, 0x121467ad0]
c     5: vbfirst_ [vbin.f: 8399, 0x121473248]
c     6: vbscf_ [vbscf.f: 1219, 0x1214b47a0]
c     7: vb_ [vbin.f: 1279, 0x1214521a0]
c     8: vbstart_ [vbin.f: 1035, 0x121451f94]
c     9: hfscf_ [master.f: 1848, 0x1205441e0]
c    10: scfgvb_ [master.f: 8240, 0x12054a050]
c    11: driver_ [master.f: 6542, 0x120548274]
c    12: gamess_ [mains.f: 236, 0x1201149f8]
c    13: main [for_main.c: 203, 0x120f6d1fc]
c    14: __start [0x120114888]
c
c  Let me know whether you've seen this on other machines, and whether the vb code
c  has been tested on the Alpha before I try and determine what is causing this.
c
c  1. all explicit specifications of "real*8" replaced with REAL in line
c  with the rest of the code.
c
c  2. subroutine scale renamed to scalev - name clash with common block in
c  mopac - although as far as I can tell the routine is never invoked.  By
c  way of interest, there are currently 176 such routines in the code ..
c  more on that later perhaps.
c
c  3. in servec.m "call prtrs(s,ndims)" corrected to "call prtrs(s,ndims,iwr)"
c
c  4.  in vbdens.m - iwr not defined in routine mkld1d2 - iwr added to argument
c  list here and in the calling routine (natorbt)
c
c  5. in vbmatre.m there is an attempt to use the syntax ..
c
c          call vberr('texta','textb')
c
c  where texta and textb have been split from the intended single character
c  string because of limitations in line width. I dont thing this will
c  work as intended, so have revised the original text to reflect the intended
c  message in a single text string.
c
c  6. single precision constant coverted to double precision.
c
c  Revision 1.25  2005/01/03 09:05:39  rwa
c  Minor bugfixes: force calculations using VB now possible, and brmax-bug fixed
c
c  Revision 1.24  2004/03/27 15:04:35  jvl
c  corrected multi-structure struc + wrong zero call
c
c  Revision 1.23  2004/01/29 12:58:35  jvl
c  Extended TURTLE with real brillouin criterion, which is now default. Fixed
c  some bugs regarding perturbation option. Changed print-out of iterations:
c  brm value added. (JJE)
c
c  Revision 1.22  2003/03/23 16:30:22  jvl
c  fixed omitted initialisation (relevant for determinant printing)
c  official error message for spin off  ; ytest => itest in servec
c
c  Revision 1.21  2003/03/14 12:39:04  jvl
c  common block for parallel was in wromg position
c
c  Revision 1.20  2003/02/18 17:17:16  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.19  2003/01/19 12:33:07  jvl
c  - check on # active orbitals (an infinite loop could result)
c  - hcpu option to allow primitive restarting of matrix element calculation in vbci
c
c  Revision 1.18  2002/12/30 23:18:03  jvl
c  corrected izero call's + improvec servec write
c
c  Revision 1.17  2002/09/05 14:40:49  jvl
c  Modifications for use with new "opti" options (instead of pert) (JJE).
c
c  Revision 1.16  2002/05/28 15:07:49  jvl
c  General cleanup
c  got rid of all fortran-files (now part of ed7 is used)
c  handle 3 more common's through include
c  make vector sections in vb more logical
c  added print parameter to getqvb
c  ....
c
c  Revision 1.15  2002/03/06 22:33:23  jvl
c  some fortran errors flagfged by hp + osinv clash vb -mopac resolved
c
c  Revision 1.14  2002/02/10 21:07:00  jvl
c  cleaned up printing and made comparison ci/scf possible
c
c  Revision 1.13  2001/11/28 10:47:00  jvl
c  forgot sizes (with common/qice)
c
c  Revision 1.12  2001/11/28 10:37:14  jvl
c  fixed the check on imax (for common qice)
c
c  Revision 1.11  2001/10/15 15:32:41  jvl
c  major parallel overhaul; makes asynchronous switchable
c  no luck yet however ...........
c
c  Revision 1.10  2001/09/19 16:02:44  jvl
c  some corrections resulting from DEC compile
c
c  Revision 1.9  2001/09/07 14:28:11  jvl
c  common stats in vb clashed with mpi; => stats_vb
c  so about 10 minutes after installing gamess on sara, the first bugfix
c
c  Revision 1.8  2001/07/03 14:13:16  jvl
c  fixed a few bugs in parallel code and allowed non-usage of peigss in vb
c  (i.e. peigss specifies if peigas is to be used in vb)
c
c  Revision 1.7  2001/06/21 16:38:20  jvl
c  added super hybrid option (experimental and perhaps useless) and a bit of cleaning up
c
c  Revision 1.6  2001/06/12 12:21:18  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.5  2001/04/14 19:40:07  jvl
c
c  changes b from linuxppc port + some vb fixes
c
c  Revision 1.4  2001/02/08 14:34:37  jvl
c  Adapted parallelisation for GAMESS-UK.
c  An MPI and GA version have been created
c
c  Revision 1.3  2000/07/07 13:41:17  hvd
c  Changes due to new operating system on DEC alpha's
c  (f77 now calls f90 internally leading to some special effects)
c
c  In the VB code the use of function iparity confused the compiler. It kept
c  complaining about passing through an array without indeces. Declaring
c  iparity external cured the problem.
c
c  Revision 1.2  2000/03/31 14:18:27  jvl
c  Changed vnorm to vnrm
c  Changed get1e to use getmat
c  Changed scopy to icopy
c  Improved on information about titles
c  Fokke Dijkstra
c
c  Revision 1.1  2000/02/16 12:20:45  jvl
c  adding vb files to gamess repository
c
c Revision 1.11  1998/02/16  15:13:49  fokke
c removed ifdef MPI
c
c Revision 1.10  1998/02/16  14:22:44  fokke
c added counter process to parallel implementation
c added output of energies and interactions of single determinants
c
c Revision 1.9  1997/05/30  11:10:21  fokke
c Changed icopy to icopy and removed variables icopy
c
c Revision 1.8  1997/05/26  08:37:38  fokke
c Changed imove to icopy
c
c Revision 1.8  1997/05/26  08:37:38  fokke
c Changed imove to icopy
c
c Revision 1.7  1997/05/22  11:31:48  joop
c Added include macro.cpp
c
c Revision 1.6  1996/10/29  14:32:24  joop
c rcs info
c
c
      subroutine meen(pacdet,idet  ,jdet  ,detcomb,dettot,
     &                icp ,jcp ,idetps,coeff  ,ig    ,
     &                igroup,trani ,tranj ,hamil  ,overl ,
     &                supers,superh,superg,ipos   ,weight,
     &                g     ,iortho,s     ,scr1   ,scr2  ,
     &                scr3  ,ndetps,q,lword,norb)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....main routine for the valence-bond matrix generator
c.....array-lengths are defined in "core".
c.....
c
c.CB.. == old 'main.f' == (try this on an IBM machine,
c.CB.. == new 'meen.f' ==  good luck with debugging...)
c
INCLUDE(common/hsinfo)
INCLUDE(common/turtleparam)
_IF(peigss)
      logical oga, pg_create
      parameter (ichunk=10,jchunk=10)
INCLUDE(common/gainfo)
_ENDIF
      common /cases/ icase(2*ncases),times(ncases),token(3*ncases),
     &               lcases
      logical lcases
c
      logical opridet
      common/pridets/opridet
c
INCLUDE(../m4/common/iofile)
c
      dimension pacdet(*),idet  (*),jdet  (*),detcomb(*),dettot(*),
     &          icp (*),jcp (*),idetps(*),coeff  (*),ig    (*),
     &          igroup(5,*),trani (*),tranj (*),hamil  (*),overl (*),
     &          supers(*),superh(*),superg(*),ipos   (*),weight(*),
     &          g     (*),iortho(*),s     (*),scr1   (*),scr2  (*),
     &          scr3  (*),ndetps(*),q(*)
      if (lcases) then
         call izero(2*ncases,icase,1)
         call  zero(ncases,times)
      end if
c.....
c.....get the vectors and the integrals
c.....
      if (iprinv.gt.10) then
         none = norb * (norb + 1) / 2
         ntwo = none * (none + 1) / 2
         write(iwr,11) none
         write(iwr,22) ntwo
11       format(//' ','number of one-electron integrals :',i10)
22       format(//' ','number of two-electron integrals :',i10)
      end if
c.....
c.....find the orthogonality classes in the integrals
c.....
      call sym1(supers,norb,iortho,northo,'noprint')
c.....
c.....read information from datatape
c.....
      call inform(nelec,nalfa,ndets,nstruc,pacdet,ndetps,idetps,coeff,
     &            ncoeff,igroup,ngroup,iortho,northo,q,lword)
c.....
c.....set ig(1 to 4,0) to 0, and ig(5,0) to 1, as the number of alpha
c.....blocks might be zero (ialfa=zero, see symblo). ig is used starting
c.....from the sixth position ! (see call hamiltt)
c.....
      ig(1) = 0
      ig(2) = 0
      ig(3) = 0
      ig(4) = 0
      ig(5) = 1
c.....
c.....prepare for outputting matrix-elements for determinants
c.....
      if (opridet) then
        memdets = ndets*(ndets+1)/2
        if (lword.le.(2*memdets+ndets)) 
     &     call vberr('out of memory representation of determinants')
        call vclr(q,1,2*memdets)
      else
        memdets=0
      end if
_IF(peigss)
      ndim = 0
      do i=1,ngroup
        ndim = ndim + igroup(2,i)
      end do
      oga = pg_create(8130,ndim,ndim,'hmat_ga',
     &                    ichunk,jchunk,iga_h)
      if (.not.oga) call ga_error
     &                  ('GA: failed to create hmat_ga',1)
      call pg_zero(iga_h)
      oga = pg_create(19130,ndim,ndim,'smat_ga',
     &                    ichunk,jchunk,iga_s)
      if (.not.oga) call ga_error
     &                  ('GA: failed to create smat_ga',1)
      call pg_zero(iga_s)
_ENDIF
c
c.....
c.....calculate the matrix-elements
c.....
      timeh = cpulft(1)
      call hamiltt(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &            idetps,coeff,ig(6),igroup,ndetps,
     &            trani,tranj,hamil,overl,scr1,scr2,scr3,s,g,
     &            supers,superh,superg(2),ipos,weight,nelec,ngroup,
     &            nalfa,iortho,northo,q(1),q(memdets+1),
     &            q(2*memdets+1))
      timeh = cpulft(1) - timeh
      if (iprinv.gt.10) write(iwr,66) timeh
66    format(/' ','the construction of the matrix representation costs '
     &       ,f8.3,' seconds')
      if (lcases) call mtype
      return
      end
      subroutine hsmat(supers,superh,superg,norb,qq)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  dynamical core allocation and h/s-matrix control routine
c
      dimension supers(*),superh(*),superg(*),qq(*)
c
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
INCLUDE(common/hsinfo)
      common /array/ narray,iarray(100)
      common /arnam/ anames(100)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/qice)
c
c...  read data and determine the length of the arrays.
c...  iarray will contain the starting adresses of the various arrays.
c...  each array will be preceeded by a highly improbable word that is
c...  checked to be there after the run. an error-message will be given
c...  in "coward" if this is not so.
c
      character*8 anames
c
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
      narray = 0
c     call putintdev(superg)
      if (imax.gt.maxcas) call vberr(
     +'Highest occupied MO index too big, Adapt arrays in /qice/'
     +  )
c
c...  pacdet will contain the packed slater determinants. these are
c...  packed per group, that is per configuration.
c  
      kpacde = igmem_alloc_inf(nwpack,'vbmatre.m','hsmat','kpacde',
     &                         IGMEM_DEBUG)
c
c...  idet and jdet contain the unpacked slater determinants per group
c 
      kidet  = igmem_alloc_inf(maxdet*nelec,'vbmatre.m','hsmat','kidet',
     &                         IGMEM_DEBUG)
      kjdet  = igmem_alloc_inf(maxdet*nelec,'vbmatre.m','hsmat','kjdet',
     &                         IGMEM_DEBUG)
c
c...  detcomb contains the matrix representation of the hamiltonian on
c...  a determinant basis, involving two groups only. dettot contains
c...  the corresponding many-electron overlap matrix.
c
      kdetco = igmem_alloc_inf(maxdet*maxdet,'vbmatre.m','hsmat',
     &                         'kdetco',IGMEM_DEBUG)
      kdetto = igmem_alloc_inf(maxdet*maxdet,'vbmatre.m','hsmat',
     &                         'kdetto',IGMEM_DEBUG)
c
c...  the matrix element routine (matre3) is given the determinants via
c...  icp and jcp (mind parity changes !).
c
      kicp = igmem_alloc_inf(nelec,'vbmatre.m','hsmat','kicp',
     &                       IGMEM_DEBUG)
      kjcp = igmem_alloc_inf(nelec,'vbmatre.m','hsmat','kjcp',
     &                       IGMEM_DEBUG)
c
c...  ndetps contains the number of determinants per structure
c
      kndetp = igmem_alloc_inf(nstruc,'vbmatre.m','hsmat','kndetp',
     &                         IGMEM_DEBUG)
c
c...  idetps contains the determinant numbers that define the structures
c...  coeff contains the corresponding coefficients
c
      kidetp = igmem_alloc_inf(ncoeff,'vbmatre.m','hsmat','kidetp',
     &                         IGMEM_DEBUG)
      kncoef = igmem_alloc_inf(ncoeff,'vbmatre.m','hsmat','kncoef',
     &                         IGMEM_DEBUG)
c
c...  ig will contain the blocking information per matrix element
c...  (5 numbers per block). but at first (in symblo) it contains
c...  5 numbers per orthogonality group => at most norb * 10 numbers !
c...  (spin-orthogonality on top of spatial orthogonality => 2 * 5)
c...  ig(1 to 4,0) is put to zero, ig(5,0)=1, this must be done to deal
c...  with the case no alpha block occur in the overlap matrix (ialfa=0)
c...   => at most norb * 10 + 5 numbers
c
      kig    = igmem_alloc_inf(10*norb+5,'vbmatre.m','hsmat','kig',
     &                         IGMEM_DEBUG)
c
c...  igroup contains the number of determinants per group, the number
c...  of structures per group and the starting position of idetps per
c...  group.
c
      kigrou = igmem_alloc_inf(nconf*5,'vbmatre.m','hsmat','kigrou',
     &                         IGMEM_DEBUG)
c
c...  trani and tranj contain the matrices that transform detcomb and
c...  dettot.
c
      ktrani = igmem_alloc_inf(maxdet*maxstr,'vbmatre.m','hsmat',
     &                         'ktrani',IGMEM_DEBUG)
      ktranj = igmem_alloc_inf(maxdet*maxstr,'vbmatre.m','hsmat',
     &                         'ktranj',IGMEM_DEBUG)
c
c...  hamil and overl contain the matrix representation of the
c...  hamiltonian on a structure basis
c
      khamil = igmem_alloc_inf(nstruc*maxstr,'vbmatre.m','hsmat',
     &                         'khamil',IGMEM_DEBUG)
      koverl = igmem_alloc_inf(nstruc*maxstr,'vbmatre.m','hsmat',
     &                         'koverl',IGMEM_DEBUG)
c
c...  ipos is the position array that is used to gather the integrals
c...  per matrix element. it is also used to reorder the matrix
c...  elements before they are written to disc (in writh)
c...  this has changed to containing row/column number info in writh
c
c
c...  weight contains the weight-factors of the integrals (cofactors)
c...  per matrix element
c
c...  first calculate the maximum number of positions in ipos/weight
c
                                nbeta = nelec - nalfa
c
c...                            one electron part
c
                                nnnnn = nalfa ** 2 + nbeta ** 2
                                nwwww = nnnnn
c
c...                            second order alfa/beta terms apart
c
                                nxxxx = (nalfa * (nalfa-1) / 2)**2 +
     &                                  (nbeta * (nbeta-1) / 2)**2
                                nnnnn = nnnnn + nxxxx * 2
                                nwwww = nwwww + nxxxx
c
c...                            mixed terms
c
                                nxxxx = (nalfa**2)*(nbeta**2)
                                nnnnn = nnnnn + nxxxx * 2
                                nwwww = nwwww + nxxxx
c
      kipos  = igmem_alloc_inf(max(nnnnn,maxstr*nstruc),'vbmatre.m',
     &                         'hsmat','kipos',IGMEM_DEBUG)
      kweigh = igmem_alloc_inf(nwwww,'vbmatre.m','hsmat','kweigh',
     &                         IGMEM_DEBUG)
c
c...  g contains the integrals that make up the matrix element
c
      kg     = igmem_alloc_inf(nnnnn,'vbmatre.m','hsmat','kg',
     &                         IGMEM_DEBUG)
c
c...  iortho contains the orthogonality number per orbital
c
      kiorth = igmem_alloc_inf(norb,'vbmatre.m','hsmat','kiorth',
     &                         IGMEM_DEBUG)
c
c...  s contains the (biorthogonalised) overlap-matrix per matrix-
c...  element
c
      ks     = igmem_alloc_inf(nalfa**2 + nbeta**2,'vbmatre.m','hsmat',
     &                         'ks',IGMEM_DEBUG)
c
c...  scr1 and scr2 are used in hamilt during the transformation of
c...  hamil and overl. scr1/scr2/scr3 are used in matre3 using
c...  always less than (nelec**2) words. scr1 is used in hamilt to
c...  gather (reorder) matrix elements before they are written to disc.
c...  and now to gather row/column info
c
      nnnnn  = max( nelec**2 , maxdet * maxstr )
      mmmmm  = max( nnnnn    , maxstr * nstruc )
      kscr1  = igmem_alloc_inf(mmmmm,'vbmatre.m','hsmat','kscr1',
     &                         IGMEM_DEBUG)
      kscr2  = igmem_alloc_inf(nnnnn,'vbmatre.m','hsmat','kscr2',
     &                         IGMEM_DEBUG)
      kscr3  = igmem_alloc_inf(nelec**2,'vbmatre.m','hsmat','kscr3',
     &                         IGMEM_DEBUG)
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr   = igmem_alloc_inf(kmemscr,'vbmatre.m','hsmat','kscr',
     &                         IGMEM_DEBUG)
c
      call meen(qq(kpacde),qq(kidet),qq(kjdet),qq(kdetco),qq(kdetto),
     &          qq(kicp),qq(kjcp),qq(kidetp),qq(kncoef),qq(kig),
     &          qq(kigrou),qq(ktrani),qq(ktranj),qq(khamil),qq(koverl),
     &          supers,superh,superg,qq(kipos),qq(kweigh),
     &          qq(kg),qq(kiorth),qq(ks),qq(kscr1),qq(kscr2),
     &          qq(kscr3),qq(kndetp),qq(kscr),kmemscr,norb)
c
      call gmem_free_inf(kscr,'vbmatre.m','hsmat','kscr')
      call gmem_free_inf(kscr3,'vbmatre.m','hsmat','kscr3')
      call gmem_free_inf(kscr2,'vbmatre.m','hsmat','kscr2')
      call gmem_free_inf(kscr1,'vbmatre.m','hsmat','kscr1')
      call gmem_free_inf(ks,'vbmatre.m','hsmat','ks')
      call gmem_free_inf(kiorth,'vbmatre.m','hsmat','kiorth')
      call gmem_free_inf(kg,'vbmatre.m','hsmat','kg')
      call gmem_free_inf(kweigh,'vbmatre.m','hsmat','kweigh')
      call gmem_free_inf(kipos,'vbmatre.m','hsmat','kipos')
      call gmem_free_inf(koverl,'vbmatre.m','hsmat','koverl')
      call gmem_free_inf(khamil,'vbmatre.m','hsmat','khamil')
      call gmem_free_inf(ktranj,'vbmatre.m','hsmat','ktranj')
      call gmem_free_inf(ktrani,'vbmatre.m','hsmat','ktrani')
      call gmem_free_inf(kigrou,'vbmatre.m','hsmat','kigrou')
      call gmem_free_inf(kig,'vbmatre.m','hsmat','kig')
      call gmem_free_inf(kncoef,'vbmatre.m','hsmat','kncoef')
      call gmem_free_inf(kidetp,'vbmatre.m','hsmat','kidetp')
      call gmem_free_inf(kndetp,'vbmatre.m','hsmat','kndetp')
      call gmem_free_inf(kjcp,'vbmatre.m','hsmat','kjcp')
      call gmem_free_inf(kicp,'vbmatre.m','hsmat','kicp')
      call gmem_free_inf(kdetto,'vbmatre.m','hsmat','kdetto')
      call gmem_free_inf(kdetco,'vbmatre.m','hsmat','kdetco')
      call gmem_free_inf(kjdet,'vbmatre.m','hsmat','kjdet')
      call gmem_free_inf(kidet,'vbmatre.m','hsmat','kidet')
      call gmem_free_inf(kpacde,'vbmatre.m','hsmat','kpacde')
c     call delintdev(superg)
      call flushn(6)
      
      return
      end
      subroutine bamilt(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &                  idetps,coeff,ig,igroup,ndetps,
     &                  hamil,overl,ic,ir,scr1,scr2,scr3,s,g,
     &                  supers,superh,superg,ipos,weight,grh,grs,diagh,
     &                  diags,nelec,ngroup,nalfa,iortho,northo,q,lword)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....in this routine the matrix representation of the hamilton
c.....operator in the configuration space is constructed.
c.....
      common/subinf/nwidth
      dimension pacdet(*),idet(nelec,*),jdet(nelec,*),detcomb(*),
     &          dettot(*),idetps(*),coeff(*),ig(5,*),igroup(5,*),
     &          ndetps(*),hamil(*),overl(*),ic(*),ir(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),supers(*),superh(*),
     &          superg(*),ipos(*),weight(*),iortho(*),grh(*),grs(*),
     &          diagh(*),diags(*),icp(*),jcp(*),q(*)
INCLUDE(common/hsinfo)
INCLUDE(common/scftvb)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
      dimension ifrom(maxact),ito(maxact),isign(maxact)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/brill)
INCLUDE(common/splice)
INCLUDE(common/ffile)
INCLUDE(common/tractlt)
c...   length of next common uncertain / please check
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
      MPIINT iiii_mpi
INCLUDE(common/vbpert)
INCLUDE(common/vbcri)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      logical oipsci
      data iagain/0/
      ind(i,j) = max(i,j) * (max(i,j)-1)/2 + min(i,j)
c.....
c.....initialise io channels
c.....
      enul = evb - core
c
_IF(parallel)
_IF(parstruc)
      call pg_dlbchunk(2,.false.)
_ELSE
      call pg_dlbchunk(200,.false.)
_ENDIF
      call pg_dlbreset()
      icounter = 0
      ido      = -1
_ENDIF
      kscr = 1
      if (fckp) then
         nn=nsa+ncore
         lennn=nn*(nn+1)/2
         kfock=kscr
         kscr=kfock+lennn
         nword=lword-kscr
         kvb7_fcmo =  kscra7vb('kvb7_fcmo',lennn,'r','n')
         if (kvb7_fcmo.lt.0) then
            call vberr('fock matrix needed in bamilt, but not found')
         else
            kvb7_fcmo =  kscra7vb('kvb7_fcmo',lennn,'r','r')
            call rdedx(q(kfock),lennn,kvb7_fcmo,num8)
         endif
      end if
c
c.....first calculate the diagonal and 1st column (gradient) of brillouin matrix
c
      call grhsdiag(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &            idetps,coeff,ig,igroup,ndetps,
     &            scr1,scr2,scr3,s,g,supers,superh,superg,ipos,
     &            weight,nelec,ngroup,nalfa,iortho,northo,q,nword,
     &            grh,grs,diagh,diags,nfock,nfockd)

c
      call classify_brill(igroup,grh,grs,diagh,diags,
     1                          coeff,pacdet,ngroup,idet,nelec,
     2               icp,jcp,supers,superh,superg,ig,ipos,weight,
     3               nalfa,scr1,scr2,scr3,s,g,iortho,northo,
     4               nbody,detcomb,dettot,nfock,nfockd)
c
      if (nbrill.ne.(npert+nvarp)) call vberr('# pert B-states + # var B
     &                         -states does not equal total # B-states')
c
c...   here ipsc static load balancing for gr and diag
c
      iflop = iipsci()
      call brtsbh(ihbloc,ihfile)
      it  = 1
      iid = 1
      nwidth = 1
c
      do 20 i=1,ngroup
         if (npert.gt.0.and.igroup(4,i).ge.2) then
           jend = 1
         else
           jend = i-1
         end if
c.....
c.....unpack 'left' determinants of group  i
c.....
         ni = igroup(1,i)
         call izero(ni*nelec,idet,1)
         call unpack(pacdet(it),n8_16,idet,ni*nelec)
         it       = it  + (ni * nelec - 1)/(64/n8_16) + 1
         jt       = (igroup(1,1) * nelec - 1)/(64/n8_16) + 2
         jid      = 1 + igroup(3,1)
         ihamil   = 1
c
         ir(1) = igroup(5,i)
c         if (.not.oipsci()) then
            hamil(1) = grh(i)
            overl(1) = grs(i)
            ic(1) = igroup(5,1)
         if (.not.oipsci()) then
            if (i.ne.1) ihamil = ihamil + 1
         end if
c         if (i.eq.1) go to 101
c...
         do 10 j=2,jend
c
            nj = igroup(1,j)
_IF(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c           write(iwr,*) 'i have to do structure',ido,' am at ',icounter
            if (icounter.ne.ido) then
               jt  = jt  + (nj * nelec - 1)/(64/n8_16) + 1
               go to 9
            end if
_ENDIF
c.....
c.....do the same thing for group  j  ('right' determinants)
	    if (igroup(4,j).ge.2.and.npert.gt.0) then
c.....         skip this one (use perturbation theory)
               jt  = jt  + (nj * nelec - 1)/(64/n8_16) + 1
               goto 9
            end if
            call izero(nj*nelec,jdet,1)
            call unpack(pacdet(jt),n8_16,jdet,nj*nelec)
            jt  = jt  + (nj * nelec - 1)/(64/n8_16) + 1
c.....
c.....now calculate the matrix elements between determinants belonging
c.....to group  i and group  j respectively
c.....
            if (exact(igroup(2,i)).and.exact(igroup(2,j)).and.
     &                                              .not.nosymc) then
               isym = iequi(igroup(2,i))
               jsym = iequi(igroup(2,j))
               if (isym.ge.jsym) then
                  jsym = 1
               else
                  isym = 1
               end if
            else if (j.eq.1.and.exact(igroup(2,i))) then
               isym = iequi(igroup(2,i))
               jsym = 1
            else
               isym = 1
               jsym = 1
            end if
            nis = ni/isym
            njs = nj/jsym
            if (igroup(4,i).eq.3) then
               jtemp=ind(ifrom(1)+ncore,ito(1)+ncore)
               dd=isign(1)*2*q(kfock+jtemp-1)
               jtemp=ind(ifrom(1)+ncore,ifrom(1)+ncore)
               fckfr=q(kfock+jtemp-1)
               jtemp=ind(ito(1)+ncore,ito(1)+ncore)
               fckto=q(kfock+jtemp-1)
               call onlys(ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  coeff(iid),coeff(jid),detcomb,dettot)
            else
              if (ndetps(i).eq.1.and.ndetps(j).eq.1) then
c.
c...           per slater det there is one with alphas and betas
c...           just the other way around, use this
               call hathab(dd,ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  idetps(iid),idetps(jid),coeff(iid),coeff(jid)
     &                                          ,detcomb,dettot)
            else
               call hatham(dd,ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  coeff(iid),coeff(jid),detcomb,dettot)
              end if
            endif
c
c....       check if we got something
c
            hamval = dd * isym * jsym
            if (dabs(hamval).gt.criign) then
            hamil(ihamil) = hamval
            ic(ihamil) = igroup(5,j)
            ir(ihamil) = igroup(5,i)
            overl(ihamil) = ds * isym * jsym
c        write(iwr,'(2I3,3F12.6)') i,j,hamil(ihamil)/overl(ihamil),
c    &              hamil(ihamil),overl(ihamil)
            ihamil = ihamil + 1
            end if
c.....
9           jid = jid +  igroup(3,j)
c.....
10       continue
         if (.not.oipsci()) then
            hamil(ihamil) = diagh(i)
            overl(ihamil) = diags(i)
            ic(ihamil) = igroup(5,i)
            ir(ihamil) = igroup(5,i)
         else
            ihamil = ihamil - 1
         end if
101      iid = iid + igroup(3,i)
c.....
c.....   apply shift to reference state (usually zero)
c.....   a corresponding shift has to be applied to the whole first column and diagonal
c.....
         if (shiscf.ne.0.0d0)  then
            if (ic(1).eq.1.and.ihamil.gt.0)
     1         hamil(1) = hamil(1) - shiscf*overl(1)
            if (iagain.eq.0.and.grs(i).ne.0.0d0.and.i.gt.1) then
               write(iwr,*) ' *** level shift for overlaps weird ',
     1                      grs(i),' ***'
               iagain = 1
            end if
         end if
c.....
c.....   project out psi(0). rotold contains the first column of the
c.....   brill. interaction matrix, before psi(0) is projected out.
c.....   rotov contains the overlaps of psi(0) with psi(i). rotgra
c.....   will contain the rotation gradients (with psi(0) projeced
c.....   out.
c.....
c         rotold(i) = hamil(1)
c         rotov(i)  = overl(1)
c         if (i.gt.1.and.3.eq.0) then
c            hi0 = hamil(1)
c            si0 = overl(1)
c            hamil(1) = hamil(1) - si0 * h00
c            overl(1) = 0.0
c            do 15 k=2,i
c               hamil(k) = hamil(k) - rotov(k) * hi0
c     &                             - si0      * rotold(k)
c     &                             + rotov(k) * si0 * h00
c               overl(k) = overl(k) - rotov(k) * si0
c15          continue
c         else
c            h00 = hamil(1)
c
c            hamil(1) = h00 / overl(1)
c            overl(1) = 1.0
c         end if
c         rotgra(i) = hamil(1)
c.....
c.....   remove zero states
c.....
_IFN(parallel)
         if (overl(ihamil).lt.cridep.and.ic(ihamil).eq.i)  then
c.....      this is no too good !, zero-overlap with itsself
_IF(mpi)
c....       enable writing from this node for a second
            iiii_mpi = ipg_nodeid()
            call MP_STDOUT_MODE(iiii_mpi)
_ENDIF
            id = iid - igroup(3,i) - 1
            write(iwr,11) i,overl(ihamil)
11          format(' remove state at i=',i4,1x,e21.15)
            do 123 k=1,ni
123         write(iwr,22) coeff(id+k),(idet(l,k),l=1,nelec)
22          format(1x,f9.5,(t11,33i3))
            overl(ihamil) = 1.0d0
_IF(mpi)
c....       enable writing from node 0 (root)
            iiii_mpi = 0
            call MP_STDOUT_MODE(iiii_mpi)
_ENDIF
         end if
_ENDIF
c
         if (igroup(4,i).eq.1) then
            call wrtsbh(hamil,overl,ic,ir,ihamil)
         end if
c....
c....       save the first column and the diagonal of both h and s
c....       brillouin has to be checked for presence
c....       diagonal is element ihamil
c....
c          diagh(i) = hamil(ihamil)
c          diags(i) = overl(ihamil)
c          if (ic(1).eq.1) then
c             grh(i) = hamil(1)
c            grs(i) = overl(1)
c        else
c           grh(i) = 0.0d0
c           grs(i) = 0.0d0
c        end if
c
20    continue
c
      call ertsbh
      call whersbt(ipbloc,nn)
      if (nn.ne.0) call vberr('ertsbh failed')
      call wrt3(grh,ngroup,ipbloc,ihfile)
      call wrt3s(diagh,ngroup,ihfile)
      call wrt3s(grs,ngroup,ihfile)
      call wrt3s(diags,ngroup,ihfile)
c
      if (npert.gt.0) then
c.....   interchange the dimension of the variational problem (nvarp)
c.....   and the total dimension (nstruc), for the davidson
         nnn = nvarp
         nvarp = nstruc
         nstruc = nnn
      end if
      return
      end

      subroutine grhsdiag(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &            idetps,coeff,ig,igroup,ndetps,
     &            scr1,scr2,scr3,s,g,supers,superh,superg,ipos,
     &            weight,nelec,ngroup,nalfa,iortho,northo,q,lword,
     &            grh,grs,diagh,diags,nfock,nfockd)
      implicit REAL  (a-h,o-z) , integer   (i-n)

c.....this subroutine will calculate the column and the diagonal element for
c.....routine bamilt, it will use fock matrix elements for column if possible or requested
      common/subinf/nwidth
INCLUDE(common/hsinfo)
INCLUDE(common/scftvb)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/brill)
INCLUDE(common/splice)
INCLUDE(common/twice)
INCLUDE(common/ffile)
INCLUDE(common/tractlt)
c...   length of next common uncertain / please check
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
      logical otovirtual, ofrdouble, ofock, ozero, oroot
INCLUDE(common/vbpert)
INCLUDE(common/vbcri)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      data icoward,ncoward/0,100/

      dimension pacdet(*),idet(nelec,*),jdet(nelec,*),detcomb(*),
     &          dettot(*),idetps(*),coeff(*),ig(5,*),igroup(5,*),
     &          ndetps(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),supers(*),superh(*),
     &          superg(*),ipos(*),weight(*),iortho(*),grh(*),grs(*),
     &          diagh(*),diags(*),icp(*),jcp(*),q(*)
      data iagain/0/
      itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      nfock = 0
      nfockd = 0
      enul = evb - core
_IF(parallel)
      call vclr(grh,1,ngroup)
      call vclr(grs,1,ngroup)
      call vclr(diagh,1,ngroup)
      call vclr(diags,1,ngroup)
_ENDIF
      iid  = 1
      it   = 1
      kscr = 1
c
      if (fckp) then
         nn=nsa+ncore
         lennn=nn*(nn+1)/2
         kfock=kscr
         kscr=kfock+lennn
         nword=lword-kscr
         kvb7_fcmo =  kscra7vb('kvb7_fcmo',lennn,'r','n')
         if (kvb7_fcmo.lt.0) then
            call vberr('fock matrix needed in grhsdiag, but not found')
         else
            kvb7_fcmo =  kscra7vb('kvb7_fcmo',lennn,'r','r')
            call rdedx(q(kfock),lennn,kvb7_fcmo,num8)
         endif
      end if
c
c      igroup = 0
      iq = 0
      i  = 0
      do inequi=1,nequi
         do ji=min(inequi,2),nex(iq+1)+1
            i = i + 1
c
c.....      do 20 i=1,ngroup
c
c
c.....It should first be checked which B-States can be calculated wiht a Fock matrix element
c.....A number of criteria exist for Fock:
c.....1  the excitation is doc --> uoc
c.....2  the from orbital (doc) is orthogonal to all other doc's
c.....3  the   to orbital (uoc) is orthogonal to all other uoc's
c.....4  the from orbital (doc) is orthogonal to the to orbital (uoc)
c.....Note that criteria 2 - 4 need the overlap matrix, which is present in supers
c
         ofock = .false.
         if (fckp.and.i.ne.1) then
           ofrdouble =.false.
           otovirtual=.true.
           ozero =.true.
           fock = 0.0d0
           fdiagfr = 0.0d0
           fdiagto = 0.0d0
           do k=iq+1,iq+iequi(inequi)
             ifr = iex(1,k)
             do i2=1,ndoubly
                if (ifr.eq.idoubly(i2)) then
                  ofrdouble =.true.
c                  exit
                end if
             end do
             do j=1,nsa
              if (j.ne.ifr) then
                if (dabs(supers(itri(ifr,j))).gt.1.0d-14) then
                   ozero = .false.
                   exit
                end if
              end if
             end do
             ito = iex(ji,k)
             isign = ieqsig(ji,k)
             if (ito.le.nscf) otovirtual = .false.
             do j=1,nsa
              if (j.ne.ito) then
                if (dabs(supers(itri(ito,j))).gt.1.0d-14) then
                   ozero = .false.
                   exit
                end if
              end if
             end do
             fock=fock + isign*2*q(kfock-1+itri(ifr+ncore,ito+ncore))
             fdiagfr = fdiagfr + q(kfock-1+itri(ifr+ncore,ifr+ncore))
             fdiagto = fdiagto + q(kfock-1+itri(ito+ncore,ito+ncore))
           end do
c           print *,'Bstate',i-1,'from',ifr,'to',ito,ofrdouble,otovirtual
c           print *,'korneel fock',ofrdouble,otovirtual,ozero,ofock
c
           if (ofrdouble.and.otovirtual.and.ozero) ofock = .true.
cjvl           if (ofock) print *,'Using fock for B-state',i-1
c
         end if
c.....We now have an ofock, which is true if gradient of the current B-state
c.....can be calculated using a fock matrix element. This means that if ofock is true
c.....that matrix element can be skipped and replaced by the appropriate fock matrix element.
c...
c......
c.....unpack 'left' determinants of group  i
c.....
         ni = igroup(1,i)
         call izero(ni*nelec,idet,1)
         call unpack(pacdet(it),n8_16,idet,ni*nelec)
         it     = it  + (ni * nelec - 1)/(64/n8_16) + 1
         jt     = 1
         jid    = 1
         ihamil = 1
c.....only the first element of row i is calculated
         j=1
c
            nj = igroup(1,j)
_IF(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c           write(iwr,*) 'i have to do structure',ido,' am at ',icounter
            if (icounter.ne.ido) then
               jt  = jt  + (nj * nelec - 1)/(64/n8_16) + 1
               go to 9
            end if
_ENDIF
c.....
c.....do the same thing for group  j  ('right' determinants)
c.....
            call izero(nj*nelec,jdet,1)
            call unpack(pacdet(jt),n8_16,jdet,nj*nelec)
            jt  = jt  + (nj * nelec - 1)/(64/n8_16) + 1
c.....
c.....now calculate the matrix elements between determinants belonging
c.....to group  i and group  j respectively
c.....
            if (exact(igroup(2,i)).and.exact(igroup(2,j)).and.
     &                                              .not.nosymc) then
               isym = iequi(igroup(2,i))
               jsym = iequi(igroup(2,j))
               if (isym.ge.jsym) then
                  jsym = 1
               else
                  isym = 1
               end if
            else if (j.eq.1.and.exact(igroup(2,i))) then
               isym = iequi(igroup(2,i))
               jsym = 1
            else
               isym = 1
               jsym = 1
            end if
            nis = ni/isym
            njs = nj/jsym
c
            if (ofock) then
               dd = fock * grs(1)
               nfock = nfock + 1
_IF(parallel)
      if (.not.oroot()) dd = 0.0d0
_ENDIF
               icoward = icoward + 1
               if (mod(icoward,ncoward).eq.0) then
                  call onlys(ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  coeff(iid),coeff(jid),detcomb,dettot)
                  ddds = ds
_IF(parallel)
                  call pg_dgop(7182,ddds,1,'+')
_ENDIF
                  if (abs(ddds).gt.1.0d-14) then
                     write(iwr,'(a12,f17.15)') ' overlap s ' ,ddds
                     call caserr('oeps no Fock')
                  end if
               end if
               ds = 0.0d0
            else
              if (ndetps(i).eq.1.and.ndetps(j).eq.1) then
c.
c...           per slater det there is one with alphas and betas
c...           just the other way around, use this
                 call hathab(dd,ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  idetps(iid),idetps(jid),coeff(iid),coeff(jid)
     &                                          ,detcomb,dettot)
           else
                 call hatham(dd,ds,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,nis,njs,idet,jdet,
     &                  coeff(iid),coeff(jid),detcomb,dettot)
               end if
            end if
c
c....       check if we got something
c
            d1d = dd * isym * jsym
            d1s = ds * isym * jsym
10          continue
c.....
c.....we have arrived at the diagonal of the h-matrix. calculate the
c.....elements of group i with itself and write the results of the i-th
c.....cycle. For now, approximate this element with a Fock matrix element.
c.....
_IF(parstruc)
         icounter = icounter + 1
         if (ido.lt.icounter) ido = ipg_dlbtask()
c        write(iwr,*) 'i have to do structure',ido,' am at ',icounter
         if (icounter.ne.ido) then
c            ihamil = ihamil - 1
            go to 12
         end if
_ENDIF
         if (fockdiagonal.and.ofock) then
            icoward = icoward + 1
            if (mod(icoward,ncoward).eq.0) then
             call onlys(ds,icp,jcp,supers,superh,
     &                superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                nelec,iortho,northo,nbody,ni,ni,idet,idet,
     &                coeff(iid),coeff(iid),detcomb,dettot)
             ddds = ds
_IF(parallel)
             call pg_dgop(7183,ddds,1,'+')
_ENDIF
             if (abs(ddds/grs(1)-2.0d0).gt.1.0d-13) then
                write(iwr,'(a12,f17.15)') ' diagonal s ' ,ddds
                call caserr('diagonal ne 2')
             end if
            end if
            ds = 2.0d0 * grs(1) 
            dd = (enul+fdiagto-fdiagfr)*ds
            nfockd = nfockd + 1
            if (.not.oroot()) then
               dd = 0.0d0
               ds = 0.0d0
            end if
         else
            call hathad(dd,ds,icp,jcp,supers,superh,
     &               superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &               nelec,iortho,northo,nbody,ni,idet,coeff(iid),
     &               detcomb,dettot)
         endif
         ddd = dd
         dds = ds

12       iid    = iid +  igroup(3,i)
c.....
c....       save the first column and the diagonal of both h and s
c....       brillouin has to be checked for presence
c....       diagonal is element ihamil
c....
          diagh(i) = ddd
          diags(i) = dds
          grh(i) = d1d
          grs(i) = d1s
_IF(parallel)
          if (i.eq.1) then
             call pg_dgop(7180,grh(1),1,'+')
             call pg_dgop(490,diagh(1),1,'+')
             call pg_dgop(7181,grs(1),1,'+')
             call pg_dgop(491,diags(1),1,'+')
          end if
_ENDIF
         end do
c
         iq = iq + iequi(inequi)
c
      end do
20    continue
c
c.....   dump first column and diagonal info (at end of h/s file)
c.....   for e.g. perturbation calculations
c
c
_IF(parallel)
      call pg_dgop(7180,grh(2),ngroup-1,'+')
      call pg_dgop(490,diagh(2),ngroup-1,'+')
      call pg_dgop(7181,grs(2),ngroup-1,'+')
      call pg_dgop(491,diags(2),ngroup-1,'+')
_ENDIF
c
c...  get max brillouin-condition
c
      brmax = 0.0d0
      enul = grh(1)/grs(1)
      do i=2,ngroup
         brmax = max(brmax,abs(grh(i)-enul*grs(i)))
      end do
      call clredx
c
      minfock = min(minfock,nfock)
      maxfock = max(maxfock,nfock)
      minfockd = min(minfockd,nfockd)
      maxfockd = max(maxfockd,nfockd)
c
      if (scri.gt.0.0) then
c...  clear supers (Olatz)
c
        kk = 0
        ll = 0
        do i=1,nsa
           do j=1,i-1
              kk =  itri(i,j)
              if (dabs(supers(kk)).lt.scri) then
                 supers(kk) = 0.0d0
                 ll = ll + 1
              end if
           end do
        end do
        kk = 0
        ss = 0.0d0
        nn = 0
        do i=1,nsa
           do j=1,i-1
              kk =  itri(i,j)
              if (dabs(supers(kk)).gt.0.0d0) nn = nn+1
              ss = dmax1(ss,supers(kk))
           end do
           kk = kk + 1
        end do
        write(iwr,*) ll,' elements of supers cleared '
        write(iwr,*) nn,' nonzero elements , largest element ',ss
        print *,' scri ',scri
	if (oprs) then
	   print *,' CLEANED S-MATRIX crit ',scri,nsa
	   call prtri(supers,nsa)
	end if
      end if
c
      return
      end

c***********************************************************************
      subroutine classify_brill(igroup,grh,grs,diagh,diags,
     1                          coeff,pacdet,ngroup,idet,nelec,
     2               icp,jcp,supers,superh,superg,ig,ipos,weight,
     3               nalfa,scr1,scr2,scr3,s,g,iortho,northo,
     4               nbody,detcomb,dettot,nfock,nfockd)
c
c...  This subroutine figures out what the user wants to
c...  do with which brillouinstructure.
c
c...  igroup - array containing information about the brillouinstructures
c...     igroup(4 will be set to 1 if Super CI is used
c...     igroup(4 will be set to 2 if perturbation theory is used
c...     igroup(4 will be set to 3 if fock is used
c...     igroup(5 will contain the index-number per type, so
c...              brillioun state 5 can be perturbation number 2 for example
c
c...     pacdet/idet only for printing and dumping excitations
c...     idet is max(nelec*max(igroup(1,*),ngroup)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      REAL grh(*),grs(*),diagh(*),diags(*)
      REAL detcomb(*), dettot(*)
      integer igroup(5,*)
      REAL pacdet(*) ! Should be pacdet(ngroup) but I cannot tell how
                     ! many words unpack might access. So I had to be
                     ! conservative in the declaration to avoid copy-in
                     ! copy-out problems. HvD.
      REAL coeff(*)
      REAL supers(*),superh(*),superg(*)
      REAL scr1(*),scr2(*),scr3(*)
      REAL s(*),g(*),weight(*)
      integer idet(nelec,*)
      integer icp(*),jcp(*),ig(*),ipos(*),iortho(*)
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/brill)
INCLUDE(common/twice)
INCLUDE(common/vbpert)
INCLUDE(common/hsinfo)
INCLUDE(common/scftvb)
INCLUDE(../m4/common/iofile)
      logical ofrdouble,otovirtual
      character*24 tekst(3)
      integer kind3(2,4)
      REAL b_max(4)
      data tekst/'   variational Bstate','perturbational Bstate',
     1           '   fock-Bstate'/
c
      ikind = 0
      if (ovbper.and.nitscf.gt.ifirsp) ikind = ntype
      Enul = diagh(1)/diags(1)
c
      ii = 0
      iq = 0
      npert = 0
      nvarp = 0
      do i=1,4
         do j=1,2
            kind3(j,i) = 0
         end do
         b_max(i) = 0.0d0
      end do
      do 40 i=1,nequi
         jstart = min(i,2)
         do 35 j=jstart,nex(iq+1)+1
            ii = ii + 1
            if (ii.gt.maxex) call caserr('# excits > maxex in clsbril')
c           do 30 k=iq+1,iq+iequi(i) only first one needed as equivalent
            ifr = iex(1,iq+1)
            ito = iex(j,iq+1)
            if (ifr.eq.1.and.ito.eq.1) then
               if (ii.ne.1) call vberr('korneel2')
               igroup(4,1) = 1
               nvarp = nvarp + 1
               igroup(5,ii) = nvarp
               cycle
            end if
            ofrdouble=.false.
            do i2=1,ndoubly
               if (ifr.eq.idoubly(i2)) ofrdouble=.true.
            end do
            otovirtual = ito.gt.nscf
c
                 itypi = 1
                 if (.not.ofrdouble) itypi = itypi + 2
                 if (otovirtual) itypi = itypi + 1
c     itypi =1 doc-voc, 2 doc-uoc, 3 voc-voc, 4  voc-uoc
                 b_max(itypi) = max(b_max(itypi),
     &                          abs(grh(ii)-Enul*grs(ii)))
c
c If ikind==3 switch on perturbation theory using "kind"
            if  (ikind.eq.3) then
                 if (itype(itypi).eq.4) then
c  perturbation theory determines switch
c  if requested (via fockdiagonal), fock matrix element is only used if pert
c  (does not make much sense for vari). When using vari,
c  element has to be recalculated (not using fock this time)
                  per=(grh(ii)-Enul*grs(ii))/(diagh(ii)-Enul*diags(ii))
                  if (abs(per).lt.crivarpert(itypi)) then
                     igroup(4,ii) = 5
                  else
                     igroup(4,ii) = 1
                  end if
                 else
                    igroup(4,ii) = itype(itypi)
                 end if
            else if (ikind.eq.2) then
c If ntype=2 switch on perturbation theory using "orbi"
               igroup(4,ii)=itype(min(ifr,ito))
            else if (ikind.eq.1) then
c If ntype=1 switch on perturbation theory using "bril"
               igroup(4,ii)=itype(ii)
            else
c... all others superci
               igroup(4,ii)=1
            endif
c
            if (igroup(4,ii).ge.2) then
               npert =  npert + 1
               igroup(5,ii) = npert
               kind3(2,itypi) =  kind3(2,itypi) + 1
            else
               nvarp = nvarp + 1
               igroup(5,ii) = nvarp
               kind3(1,itypi) =  kind3(1,itypi) + 1
            end if
35       continue
40    iq = iq + iequi(i)
c.... if no b-statess with pert switch off perturbation theory
      if (npert.eq.0) ovbper=.false.
c
      if (ii.ne.ngroup) call caserr('clsbril confusion')
      k7igr_brill =  kscra7vb('k7igr_brill',ngroup*5,'i','w')
      call wrt3i(igroup,ngroup*5,k7igr_brill,num8)
c
c...  Print out information about Brillouin States
c
       minpert = min(minpert,npert)
       maxpert = max(maxpert,npert)
       minvar = min(minvar,nvarp)
       maxvarp = max(maxvarp,nvarp)
      if (iprinv.ge.1) then
         write(iwr,612) nvarp,npert
612      format(/' variation structures :',i6,' pert structures :',i6)
         write(iwr,613) ' variation    ',(kind3(1,j),j=1,4)
         write(iwr,613) ' perturbation ',(kind3(2,j),j=1,4)
         write(iwr,614) ' brillouin    ',(b_max(j),j=1,4)
613      format(a14,' doc-voc ',i7,' doc-uoc ',i7,
     &              ' voc-voc ',i7,' voc-uoc ',i7)
614      format(a14,' doc-voc ',1pg7.1,' doc-uoc ',1pg7.1,
     &              ' voc-voc ',1pg7.1,' voc-uoc ',1pg7.1)
         if (nfock.gt.0.and.nfockd.gt.0) write(iwr,'(a,i7,a)')  
     &       ' fock used for',nfock,
     &       ' brillouin matrix-elements and diagonals'
         if (nfock.gt.0.and.nfockd.eq.0) write(iwr,'(a,i7,a)')  
     &       ' fock used for',nfock,' brillouin matrix-elements'
         if (nfock.gt.0.and.nfockd.gt.0.and.nfock.ne.nfockd) 
     &      call caserr('nfock ne nfockd')
         if (iprinv.ge.5) then
         write(iwr,'(a)') '  // list of Brillouin States \\'
         it = 1
         it2 = 0
         ii = 0
         iq = 0
         do 30 i=1,nequi
            do 38 j=min(i,2),nex(iq+1)+1
               ii = ii + 1
               if (ii>1) write(iwr,602) tekst(min(igroup(4,ii),2)),ii-1,
     &                (grh(ii)-Enul*grs(ii))/(diagh(ii)-Enul*diags(ii)),
     &                (grh(ii)-Enul*grs(ii)),
     &                   (iex(1,k),iex(j,k),k=iq+1,iq+iequi(i))
602            format('  ',a21,i4,': Bt ',1pe12.5,' Pt-c ',1pe12.5,
     &                (t65,'// ',5(i4,'   =>',i4,' \\')))
               if (ii==1) write(iwr,'(a)') '     Groundstate - Psi0 '
               if (iprinv.ge.1000) then
                call unpack(pacdet(it),n8_16,idet,igroup(1,ii)*nelec)
                write(iwr,603)
603             format(7x,'det      coeff   orbitals =>')
                do  iii=1,igroup(1,ii)
                 it2 = it2 + 1
                 write(iwr,604) iii,coeff(it2),(idet(ij,iii),ij=1,nelec)
604              format(5x,i6,1x,f9.6,1x,(t25,25i4),1x)
                end do
               end if
38            it = it + (nelec*igroup(1,ii)-1)/(64/n8_16) + 1
30       iq = iq + iequi(i)
         end if
      end if
c
      return
      end
c***********************************************************************

      subroutine bior(s,ir,nr,ic,nc,det,irank,ipar,watch,singu)
c.....
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the  matrix s(nr,nc) is decomposed in
c.....a left-under/diagonal/right-upper matrix ("biorthogonalised").
c.....ir and ic keep track of the permutations that are induced
c.....by the pivot search. a rank less than min(nr,nc) is allowed for.
c.....irank will contain the rank of s. det is the product
c.....of the non-zero diagonal elements times the parity change.
c.....only in square, non-singular cases det in fact is the determinant.
c.....
      logical watch,singu
INCLUDE(common/vbcri)
c...
      dimension s(nr,nc),ir(nr),ic(nc)
      watch = .false.
      singu = .false.
      irank = 0
      det   = 1.0d0
      ipar  = 1
      call pivot(s,nr,nc,1,piv,ir,ic,ipar)
      if (dabs(piv).le.critor) then
         nr    = 0
         nc    = 0
         det   = 0.0d0
         watch = .true.
         return
      end if
      irank = 1
      do 40 i = 1,min(nc-1,nr)
         det = det * s(i,i)
         do 30 j = i+1,nc
c.....
c.....form i-th row of r
c.....
            s(i,j) = -s(i,j)/s(i,i)
c.....
c.....adapt s-matrix to r
c.....
            do 10 k = i+1,nr
               s(k,j) = s(k,j) + s(k,i)*s(i,j)
10          continue
c.....
c.....adapt r-matrix to itself
c.....
            do 20 l = 1,i-1
               s(l,j) = s(l,j) + s(l,i)*s(i,j)
20          continue
30       continue
         call pivot(s,nr,nc,i+1,piv,ir,ic,ipar)
         if (dabs(piv).le.critor) then
            goto 123
         end if
         irank = irank + 1
40    continue
c.....
c.....r-matrix is constructed. l's turn now
c.....
      if(nc-1.lt.nr) det = det * s(nc,nc)
123   do 100 i=2,nr
c.....
c.....   construct i-th row of l
c.....
         iii = min(i-1,irank)
         do 50 j=1,iii
            s(i,j) = -s(i,j)/s(j,j)
50       continue
c.....
c.....   adapt l to itself
c.....
         sum = s(i,1)
         do 70 k=1,iii
            do 60 l=k+1,iii
               sum = sum + s(i,l) * s(l,k)
60          continue
            s(i,k) = sum
            sum = s(i,k+1)
70      continue
100   continue
c.....
c.....investigate dependences and, if so, take proper action.
c.....
      if (min(nr,nc).gt.irank.or.nr.ne.nc) then
         do 120 j=nc,1,-1
            do 110 i=1,nr
               if (dabs(s(i,j)).gt.critor) goto 130
110         continue
            nc    = nc - 1
            watch = .true.
120      continue
130      do 150 i=nr,1,-1
            do 140 j=1,nc
               if (dabs(s(i,j)).gt.critor) goto 160
140         continue
            nr    = nr - 1
            watch = .true.
150      continue
160      continue
         if (min(nr,nc).gt.irank) singu = .true.
      end if
      return
      end
      subroutine c00(s,ndim,w)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the first order cofactors are calculated
c.....in the montfoort/broer-way
c.....
      dimension s(ndim,ndim),w(ndim*ndim)
c.....
      it = 0
      do 30 k=1,ndim-1
         do 20 i=1,ndim-1
            it = it + 1
            w(it) = s(k,ndim) * s(ndim,i)
20       continue
         it = it + 1
         w(it) = s(k,ndim)
30    continue
      do 40 i=1,ndim-1
         it = it + 1
         w(it) = s(ndim,i)
40    continue
      it = it + 1
      w(it) = 1.0d0
      return
      end
      subroutine c0000(ndim,w,x,y,d,s)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension w((ndim*(ndim-1)/2)**2),s(ndim,ndim),d(ndim),
     &          x(ndim,ndim),y(ndim,ndim)
c.....
c.....the overlap-matrix is singular. use factorised cofactor algorithm
c.....
      do 10 i=1,ndim
         d(i) = s(i,i)
         s(i,i) = 1.0d0
10    continue
      do 120 i=1,ndim
         do 110 j=1,ndim
            y(i,j) = s(ndim,i) * s(j,ndim)
110      continue
120   continue
      call vclr(x,1,ndim*ndim)
      do 150 i=1,ndim
         do 140 j=1,ndim
            do 130 k=max(i,j),ndim-1
               x(i,j) = x(i,j) + s(k,i) * s(j,k) / d(k)
130         continue
140      continue
150   continue
      it = 0
      do 190 k=2,ndim
         do 180 l=1,k-1
            do 170 i=2,ndim
               do 160 j=1,i-1
                  it = it + 1
                  w(it) = y(j,l) * x(i,k) - y(j,k) * x(i,l)
     &                  + y(i,k) * x(j,l) - y(i,l) * x(j,k)
160            continue
170         continue
180      continue
190   continue
      return
      end
      subroutine c0k(s,ncol,w)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension s(ncol-1,ncol),w(*)
      do 10 i=1,ncol-1
         w(i) = s(i,ncol)
10    continue
      w(ncol) = 1.0d0
      return
      end
      subroutine c0k0l(ncol,s,w,x,y)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension s(ncol-2,ncol),w(*),x(*),y(*)
      do 10 i=1,ncol-2
         x(i) = s(i,ncol-1)
         y(i) = s(i,ncol  )
10    continue
      x(ncol-1) = 1.0d0
      x(ncol  ) = 0.0d0
      y(ncol-1) = 0.0d0
      y(ncol  ) = 1.0d0
      it = 0
      do 30 k=2,ncol
         do 20 l=1,k-1
            it = it + 1
            w(it) = y(k) * x(l) - y(l) * x(k)
20       continue
30    continue
      return
      end
      subroutine c0kjl(ncol,w,x,y,d,s)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
      dimension w(*),x(ncol,ncol),y(*),d(*),s(ncol-1,ncol)
      nrow = ncol - 1
      do 10 i=1,nrow
         d(i)   = s(i,i)
         s(i,i) = 1.0d0
         y(i)   = s(i,ncol)
10    continue
      y(ncol) = 1.0d0
      call vclr(x,1,ncol*ncol)
      do 40 i=1,nrow
         do 30 j=1,nrow
            do 20 k=max(i,j),nrow
               x(i,j) = x(i,j) + s(k,i) * s(j,k) / d(k)
20          continue
30       continue
40    continue
      it = 0
      do 70 k=2,ncol
         do 60 l=1,k-1
            do 50 j=1,nrow
               it = it + 1
               w(it) = y(k) * x(j,l) - y(l) * x(j,k)
50          continue
60       continue
70    continue
      return
      end
      subroutine c2222(ndim,w,s,x,y)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension w(*),s(ndim,ndim),x(ndim,ndim),y(ndim,ndim)
      s(ndim-1,ndim-1) = 1.0d0
      s(ndim  ,ndim  ) = 1.0d0
      do 30 i=1,ndim
         do 20 j=1,ndim
            y(i,j) = s(ndim,i) * s(j,ndim)
20       continue
30    continue
      call vclr(x,1,ndim*ndim)
      do 50 i=1,ndim-1
         do 40 j=1,ndim-1
            x(i,j) = s(ndim-1,i) * s(j,ndim-1)
40       continue
50    continue
      it = 0
      do 90 k=2,ndim
         do 80 l=1,k-1
            do 70 i=2,ndim
               do 60 j=1,i-1
                  it = it + 1
                  w(it) = y(j,l) * x(i,k) - y(j,k) * x(i,l)
     &                  + y(i,k) * x(j,l) - y(i,l) * x(j,k)
60             continue
70          continue
80       continue
90    continue
      return
      end
      subroutine c2kjl(ncol,w,x,y,s)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
      dimension w(*),x(ncol,ncol),y(*),s(ncol-1,ncol)
      nrow = ncol - 1
      s(nrow,nrow) = 1.0d0
      do 10 i=1,nrow
         y(i)   = s(i,ncol)
10    continue
      y(ncol) = 1.0d0
      do 30 i=1,nrow
         do 20 j=1,nrow
            x(i,j) = s(nrow,i) * s(j,nrow)
20       continue
         x(i,ncol) = 0.0d0
30    continue
      it = 0
      do 70 k=2,ncol
         do 60 l=1,k-1
            do 50 j=1,nrow
               it = it + 1
               w(it) = y(k) * x(j,l) - y(l) * x(j,k)
50          continue
60       continue
70    continue
      return
      end
      subroutine ci0t(s,nrow,w)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension s(nrow,nrow-1),w(*)
      do 10 j=1,nrow-1
         w(j) = s(nrow,j)
10    continue
      w(nrow) = 1.0d0
      return
      end
      subroutine ci0j0(nrow,s,w,x,y)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension s(nrow,nrow-2),w(*),x(*),y(*)
      do 10 i=1,nrow-2
         x(i) = s(nrow-1,i)
         y(i) = s(nrow  ,i)
10    continue
      x(nrow-1) = 1.0d0
      x(nrow  ) = 0.0d0
      y(nrow-1) = 0.0d0
      y(nrow  ) = 1.0d0
      it = 0
      do 30 i=2,nrow
         do 20 j=1,i-1
            it = it + 1
            w(it) = y(i) * x(j) - y(j) * x(i)
20       continue
30    continue
      return
      end
      subroutine ci0jl(nrow,w,x,y,d,s)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
      dimension w(*),x(nrow,nrow),y(*),d(*),s(nrow,nrow-1)
      ncol = nrow - 1
      do 10 i=1,ncol
         d(i)   = s(i,i)
         s(i,i) = 1.0d0
         y(i)   = s(nrow,i)
10    continue
      y(nrow) = 1.0d0
      call vclr(x,1,nrow*nrow)
      do 40 i=1,ncol
         do 30 j=1,ncol
            do 20 k=max(i,j),ncol
               x(i,j) = x(i,j) + s(k,i) * s(j,k) / d(k)
20          continue
30       continue
40    continue
      it = 0
      do 70 l=1,ncol
         do 60 i=2,nrow
            do 50 j=1,i-1
               it = it + 1
               w(it) = y(i) * x(j,l) - y(j) * x(i,l)
50          continue
60       continue
70    continue
      return
      end
      subroutine ci2jl(nrow,w,x,y,s)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
      dimension w(*),x(nrow,nrow),y(*),s(nrow,nrow-1)
      ncol = nrow - 1
      s(ncol,ncol) = 1.0d0
      do 10 i=1,ncol
         y(i)   = s(nrow,i)
10    continue
      y(nrow) = 1.0d0
      do 30 i=1,ncol
         do 20 j=1,ncol
            x(i,j) = s(ncol,i) * s(j,ncol)
20       continue
         x(nrow,i) = 0.0d0
30    continue
      it = 0
      do 70 l=1,ncol
         do 60 i=2,nrow
            do 50 j=1,i-1
               it = it + 1
               w(it) = y(i) * x(j,l) - y(j) * x(i,l)
50          continue
60       continue
70    continue
      return
      end
      subroutine cik(s,w,scr1,d,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the first order cofactors are calculated
c.....in the prosser-and-hagstrom-way (f. prosser and s. hagstrom (1968)
c.....int.j.quantum chem. 2,89).
c.....
      dimension s(ndim,ndim),w(ndim*ndim),scr1(ndim,ndim),d(ndim)
c.....
c.....copy s(i,i) into d(i), make s(i,i) 1.0
c.....
      do 10 i=1,ndim
         d(i) = s(i,i)
         s(i,i) = 1.0d0
10    continue
c.....
c.....determine adjd * l
c.....
      scr1(1,1) = 1.0 d0/ d(1)
      do 30 i=1,ndim
         do 20 j=1,i
            scr1(i,j) = s(i,j) / d(i)
20       continue
30    continue
c.....
c.....calculation of adjugate s/transpose !/  as r * adjd * l
c.....
      call vclr(w,1,ndim*ndim)
      it = 0
      do 60 k=1,ndim
         do 50 i=1,ndim
            it = it + 1
            do 40 m=max(i,k),ndim
               w(it) = w(it) + s(k,m) * scr1(m,i)
40          continue
50       continue
60    continue
      return
      end
      subroutine cikjl(w2,n,w1)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension w2(*),w1(n,n)
c.....
c.....use jacobi ratio theorem
c.....
      it = 0
      do 40 k=2,n
         do 30 l=1,k-1
            do 20 i=2,n
               do 10 j=1,i-1
                  it = it + 1
                  w2(it) =  w1(i,k) * w1(j,l)
     &                   -  w1(i,l) * w1(j,k)
10             continue
20          continue
30       continue
40    continue
      return
      end
c     ============================================================
      subroutine gmix(ipos,ipose,ir,ic,ig,nblock,ialfa)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
      ! integer :: v
      ! integer, dimension(nblock) :: ma
c
c      real start, stop

      logical equal
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)

c      print *, "Start of subroutine gmix"
c      print *, "nblock=", nblock
c      print *, "ialfa=", ialfa
c      print *, ""
c      call cpu_time(start)
      it = 0
c      printer_inner = .false.
c      print_it = .false.
c     Unroll m loops 
      do 60 m=1,nblock-1
c        from here all values are known, so parallelize from here down
         do 50 l=ig(4,m),ig(4,m+1)-1
            do 40 j=ig(3,m),ig(3,m+1)-1
c              maybe even from here 
               do 30 n=m+1,nblock
                  do 20 k=ig(4,n),ig(4,n)+ig(2,n)-1
                     do 10 i=ig(3,n),ig(3,n)+ig(1,n)-1
                        it = it + 1
                        ! print *,"it:",it,"ikjl:",1000*i+100*k+10*j+l
                        ipos(it) = intpos(ir(i),ic(k),ir(j),ic(l))
10                   continue
20                continue
30             continue
40          continue
50       continue
60    continue
      ljki = 1000000*ig(4,1)+10000*ig(3,1)+100*ig(4,2)+ig(3,2)
      if (it /= 81.and.it/=194481.or.nblock/=2.or.ialfa/=1.or.
     &ljki/=1010404.and.ljki/=1012222) then
      print *, "NEW CASE FOUND"
      print *, "First loop iterations:", it
      print *, "nblock:", nblock, "ialfa:", ialfa
      print *, "ljki:", 1000000*ig(4,1)+10000*ig(3,1)
     &+100*ig(4,2)+ig(3,2)
      print *, "ljki:", 1000000*l+10000*j+100*k+i 
      print*,"========================================================="
      end if
      ! if (it /= ix) then
      ! stop
      ! end if

      call izero(it,ipose,1)
      na    = ig(5,ialfa) + ig(1,ialfa) * ig(2,ialfa) - 1
      nb    = ig(5,nblock) + ig(1,nblock) * ig(2,nblock) - na - 1
      it = 0
      do 120 m=1,ialfa-1
         do 110 l=ig(4,m),ig(4,m+1)-1
            do 100 j=ig(3,m),ig(3,m+1)-1
               do 90 n=m+1,ialfa
                  do 80 k=ig(4,n),ig(4,n)+ig(2,n)-1
                     do 70 i=ig(3,n),ig(3,n)+ig(1,n)-1
                        it        = it + 1
                        ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
70                   continue
80                continue
90             continue
               it = it + nb
100         continue
110      continue
120   continue

c      intermediate = it
c      print *, "Second loop iterations:", it

      it = it + ig(1,ialfa)*ig(2,ialfa)*nb
      do 180 m=ialfa+1,nblock-1
         do 170 l=ig(4,m),ig(4,m+1)-1
            do 160 j=ig(3,m),ig(3,m+1)-1
               do 150 n=m+1,nblock
                  do 140 k=ig(4,n),ig(4,n)+ig(2,n)-1
                     do 130 i=ig(3,n),ig(3,n)+ig(1,n)-1
                        it = it + 1
                        ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
130                  continue
140               continue
150            continue
160         continue
170      continue
180   continue
c      print *, "Third loop iterations:", it - intermediate
c      call cpu_time(stop)
c      print *, start
c      print *, stop 
c      print *, ""
      return
      end
c     ============================================================
      subroutine gmix0(ipos,ipose,ir,ic,ig,nblock,is1,ialfa)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 50 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
         do 40 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            do 30 m=1,nblock
               if (m.ne.is1) then
                  do 20 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 10 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        n = n + 1
                        ipos(n) = intpos(ir(i),ic(k),ir(j),ic(l))
10                   continue
20                continue
               end if
30          continue
40       continue
50    continue
      call izero(n,ipose,1)
      it = 0
      if (is1.le.ialfa) then
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         nb = ig(5,nblock) + ig(1,nblock) * ig(2,nblock) - na - 1
         do 100 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 90 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               do 80 m=1,ialfa
                  if (m.ne.is1) then
                     do 70 l=ig(4,m),ig(4,m)+ig(2,m)-1
                        do 60 j=ig(3,m),ig(3,m)+ig(1,m)-1
                           it = it + 1
                           ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
60                      continue
70                   continue
                  end if
80             continue
               it = it + nb
90          continue
100      continue
      else
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         do 150 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 140 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               it = it + na
               do 130 m=ialfa+1,nblock
                  if (m.ne.is1) then
                     do 120 l=ig(4,m),ig(4,m)+ig(2,m)-1
                        do 110 j=ig(3,m),ig(3,m)+ig(1,m)-1
                           it = it + 1
                           ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
110                     continue
120                  continue
                  end if
130            continue
140         continue
150      continue
      end if
      return
      end
      subroutine gmix00(ipos,ipose,ir,ic,ig,is0,is02,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      if (equal) then
         do 40 k=ig(4,is02),ig(4,is02)+ig(2,is02)-1
            do 30 i=ig(3,is02),ig(3,is02)+ig(1,is02)-1
               do 20 l=ig(4,is0),ig(4,is0)+ig(2,is0)-1
                  do 10 j=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                     n  = n + 1
                     ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                     ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10                continue
20             continue
30          continue
40       continue
      else
         do 80 k=ig(4,is02),ig(4,is02)+ig(2,is02)-1
            do 70 i=ig(3,is02),ig(3,is02)+ig(1,is02)-1
               do 60 l=ig(4,is0),ig(4,is0)+ig(2,is0)-1
                  do 50 j=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                     n  = n + 1
                     ipos(n) = intpos(ir(i),ic(k),ir(j),ic(l))
50                continue
60             continue
70          continue
80       continue
      end if
      return
      end
      subroutine gmix0f(ipos,ipose,ir,ic,ig,icfix,is1,is0,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      if (equal) then
         do 30 j=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            do 20 k=ig(4,is0),ig(4,is0)+ig(2,is0)-1
               do 10 i=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),icfix)
                  ipose(n) = intpos(ir(i),icfix,ir(j),ic(k))
10             continue
20          continue
30       continue
      else
         do 60 j=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            do 50 k=ig(4,is0),ig(4,is0)+ig(2,is0)-1
               do 40 i=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),icfix)
40             continue
50          continue
60       continue
      end if
      return
      end
      subroutine gmix0k(ipos,ipose,ir,ifix,ic,ig,nblock,nalfa,ialfa,is1)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      it    = 0
      irfix = ir(ifix)
      do 40 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
         do 30 m=1,nblock
            if (m.ne.is1) then
               do 20 l=ig(4,m),ig(4,m)+ig(2,m)-1
                  do 10 j=ig(3,m),ig(3,m)+ig(1,m)-1
                     it = it + 1
                     ipos(it) = intpos(irfix,ic(k),ir(j),ic(l))
10                continue
20             continue
            end if
30       continue
40    continue
      call izero(it,ipose,1)
      it = 0
      if (ig(3,is1).le.nalfa) then
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         nb = ig(5,nblock) + ig(1,nblock) * ig(2,nblock) - na - 1
         do 80 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 70 m=1,ialfa
               if (m.ne.is1) then
                  do 60 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 50 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        it = it + 1
                        ipose(it) = intpos(irfix,ic(l),ir(j),ic(k))
50                   continue
60                continue
               end if
70          continue
         it = it + nb
80       continue
      else
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         do 120 k=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            it = it + na
            do 110 m=ialfa+1,nblock
               if (m.ne.is1) then
                  do 100 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 90 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        it = it + 1
                        ipose(it) = intpos(irfix,ic(l),ir(j),ic(k))
90                   continue
100               continue
               end if
110         continue
120      continue
      end if
      return
      end
      subroutine gmixab(ipos,ipose,ir,ic,ig,nblock,ialfa,is1,is2)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 50 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
         do 40 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            do 30 m=1,nblock
               if (m.ne.is1.and.m.ne.is2) then
                  do 20 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 10 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        n = n + 1
                        ipos(n) = intpos(ir(i),ic(k),ir(j),ic(l))
10                   continue
20                continue
               end if
30          continue
40       continue
50    continue
      call izero(n,ipose,1)
      it = 0
      if (is1.le.ialfa) then
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         nb = ig(5,nblock) + ig(1,nblock) * ig(2,nblock) - na - 1
         do 100 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
            do 90 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               do 80 m=1,ialfa
                  if (m.ne.is1.and.m.ne.is2) then
                     do 70 l=ig(4,m),ig(4,m)+ig(2,m)-1
                        do 60 j=ig(3,m),ig(3,m)+ig(1,m)-1
                           it = it + 1
                           ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
60                      continue
70                   continue
                  end if
80             continue
               it = it + nb
90          continue
100      continue
      else
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         do 150 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
            do 140 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               it = it + na
               do 130 m=ialfa+1,nblock
                  if (m.ne.is1.and.m.ne.is2) then
                     do 120 l=ig(4,m),ig(4,m)+ig(2,m)-1
                        do 110 j=ig(3,m),ig(3,m)+ig(1,m)-1
                           it = it + 1
                           ipose(it) = intpos(ir(i),ic(l),ir(j),ic(k))
110                     continue
120                  continue
                  end if
130            continue
140         continue
150      continue
      end if
      return
      end
      subroutine gmixf0(ipos,ipose,ir,ic,ig,irfix,is1,is0,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      if (equal) then
         do 30 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 20 k=ig(4,is0),ig(4,is0)+ig(2,is0)-1
               do 10 i=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),irfix,ic(l))
                  ipose(n) = intpos(ir(i),ic(l),irfix,ic(k))
10             continue
20          continue
30       continue
      else
         do 60 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 50 k=ig(4,is0),ig(4,is0)+ig(2,is0)-1
               do 40 i=ig(3,is0),ig(3,is0)+ig(1,is0)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),irfix,ic(l))
40             continue
50          continue
60       continue
      end if
      return
      end
      subroutine gmixi0(ipos,ipose,ir,ic,kfix,ig,nblock,nalfa,ialfa,is1)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      it    = 0
      icfix = ic(kfix)
      do 40 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
         do 30 m=1,nblock
            if (m.ne.is1) then
               do 20 l=ig(4,m),ig(4,m)+ig(2,m)-1
                  do 10 j=ig(3,m),ig(3,m)+ig(1,m)-1
                     it = it + 1
                     ipos(it) = intpos(ir(i),icfix,ir(j),ic(l))
10                continue
20             continue
            end if
30       continue
40    continue
      call izero(it,ipose,1)
      it = 0
      if (ig(3,is1).le.nalfa) then
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         nb = ig(5,nblock) + ig(1,nblock) * ig(2,nblock) - na - 1
         do 80 i=ig(3,is1),ig(3,is1)+ig(1,is1)
            do 70 m=1,ialfa
               if (m.ne.is1) then
                  do 60 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 50 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        it = it + 1
                        ipose(it) = intpos(ir(i),ic(l),ir(j),icfix)
50                   continue
60                continue
               end if
70          continue
         it = it + nb
80       continue
      else
         na = ig(5,ialfa)+ig(1,ialfa)*ig(2,ialfa)-1
         do 120 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            it = it + na
            do 110 m=ialfa+1,nblock
               if (m.ne.is1) then
                  do 100 l=ig(4,m),ig(4,m)+ig(2,m)-1
                     do 90 j=ig(3,m),ig(3,m)+ig(1,m)-1
                        it = it + 1
                        ipose(it) = intpos(ir(i),ic(l),ir(j),icfix)
90                   continue
100               continue
               end if
110         continue
120      continue
      end if
      return
      end
      subroutine coward(q,istage)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....if istage = 0 at the addresses iarray(i)-1 in blank common a
c.....strange word is written. if istage = 1 it is checked they are
c.....still there. istage = 3 will result in istage = 0 if trouble is
c.....true (ment for debugging)
c.....
      logical trouble
      dimension q(*)
      common /array/ narray,iarray(100)
      common /arnam/ anames(100)
INCLUDE(../m4/common/iofile)
      character*8 anames
      strange = 9.8765432123456789d10
      trouble = .false.
      return
      if (istage.eq.0) then
         do 10 i=2,narray
            q( iarray(i) - 1 ) = strange
10       continue
      else if(istage.eq.1.or.istage.eq.3) then
         do 20 i=2,narray
            if ( q( iarray(i) - 1 ).ne.strange) then
               write(iwr,11) anames(i),anames(i-1)
               trouble = .true.
            end if
20       continue
      else
         write(iwr,22) istage
         stop
      end if
      if (istage.eq.3.and.trouble) then
         istage = 0
         return
      end if
11    format(' '//,'*************************'
     &,       //1x,'*** array  ',a8,' has been overwritten    ***'
     &,       //1x,'** by     ',a8,'                         **'
     &,        //,'*************************')
22    format(' '//,'wrong call to coward,istage = ',i8)
      if (trouble) then
         do 30 i=1,narray-1
            write(iwr,33) anames(i),iarray(i+1)-iarray(i)-1
30       continue
33       format(' ','array ',a8,' may contain ',i10,' elements')
         stop
      end if
      return
      end
      subroutine g0ff(ipos,ipose,ir,ic,nr,irfix,icfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*)
      logical equal
      common /posit/ iky(3)
      n = 0
      if (equal) then
         do 20 l=1,nr
            do 10 j=1,nr
               n  = n + 1
               ipos (n) = intpos(irfix,icfix,ir(j),ic(l))
               ipose(n) = intpos(irfix,ic(l),ir(j),icfix)
10          continue
20       continue
      else
         do 40 l=1,nr
            do 30 j=1,nr
              n  = n + 1
              ipos (n) = intpos(irfix,icfix,ir(j),ic(l))
30         continue
40      continue
      end if
      return
      end
      subroutine getdet(pacdet,idet,ndets,nstruc,nelec,tran,ndetps,
     &                  idetps,coeff)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....this routine unpacks the determinants for a certain group of
c.....structures and constructs the matrix that transforms them into
c.....the structures.
c.....
c.....pacdet contains the packed determinants. idet will contain the
c.....unpacked determinants, tran will be the transformation matrix,
c.....ndetps gives the number of determinants per structures, idetps
c.....gives the determinants per structure and coeff contains the
c.....transformation coefficients.
c.....
INCLUDE(common/c8_16vb)
c...   length of next common uncertain / please check
      dimension pacdet(ndets * ( (nelec-1) * n8_16/64 + 1 )),
     &          idet(ndets*nelec), tran(ndets,nstruc), ndetps(nstruc),
     &          idetps(*),coeff(ndets)
c.....
c.....unpack all determinants for the group at hand, at once.
c.....
      call izero(ndets*nelec,idet,1)
      call unpack(pacdet,n8_16,idet,ndets * nelec)
c.....
c.....form transformation matrix
c.....
      call vclr(tran,1,ndets * nstruc)
      it = 1
c.....
c.....it should be noticed here, arbitrarily, it is decided that a
c.....relative numbering is used. i.e. per group of structures the first
c.....one will be denoted 1.
c.....
      do 20 j=1,nstruc
         do 10 i=1,ndetps(j)
            tran(idetps(it),j) = coeff(it)
            it = it + 1
10       continue
20    continue
c
      return
      end
      subroutine getpac(nelec,nalfa,iscr,npack,ngroup,lword,imax,maxdet)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....get packed determinants and group (that is configuration) info
c.....from tape25. lword is the maximum length of iscr. on return npack
c.....is the number of words used by pacdet.
c.....
_IFN(atmol)
INCLUDE(../m4/common/iofile)
_ENDIF
INCLUDE(common/c8_16vb)
      dimension iscr(lword)
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
      ngroup = nconf
c
c     read number of determinants per structure
c
      if (nstruc.gt.lword) call corfait(nstruc,lword,'in getpac/1')
_IF(atmol)
      read (25) (iscr(i),i=1,nstruc)
_ELSE
      call readis(iscr,nstruc,num8)
_ENDIF
c
c     read determinant-numbers / dummy here
c
      ndettot = 0
      do 10 i=1,nstruc
         ndettot = ndettot + iscr(i)
10    continue
      if (ndettot.gt.lword) call corfait(ndettot,lword,'in getpac/2')
_IF(atmol)
      read (25) (iscr(i),i=1,ndettot)
_ELSE
      call readis(iscr,ndettot,num8)
_ENDIF
c
c     read coefficients / dummy here
c
capollo nipw() is the number of integers per word
_IF(atmol)
      read (25) (iscr(i),i=1,ndettot*nipw())
_ELSE
      call readis(iscr,ndettot*nipw(),num8)
_ENDIF
c
c                read number of determinants per configuration
c
_IF(atmol)
      read(25) (iscr(i),i=1,nconf)
_ELSE
      call readis(iscr,nconf,num8)
_ENDIF
      maxdet = 0
      do 123 i=1,nconf
123   if (iscr(i).gt.maxdet) maxdet = iscr(i)
 
c
c                 read number of structures per configuration
c
_IF(atmol)
      read(25) (iscr(i),i=nconf+1,nconf+nconf)
_ELSE
      call readis(iscr(nconf+1),nconf,num8)
_ENDIF
c
c              read packed slater determinants
c
      npack = 0
capollo mind packed words take 64 bits => nipw()
      do 66 i=1,nconf
         npack = npack + (((iscr(i) * nelec - 1)/(64/n8_16)) + 1)*nipw()
66    continue
      nn = nconf+npack
      if (nn.gt.lword) call corfait(nn,lword,'in getpac/3')
_IF(atmol)
      read(25) (iscr(i),i=nconf+1,nconf+npack)
_ELSE
      call readis(iscr(nconf+1),npack,num8)
_ENDIF
      return
      end
      subroutine getpa2 (nelec,nalfa,iscr,npack,nconf,nstruc,ndettot,
     &                                                lword,imax)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....same as getpac, only this one also returns the number of
c.....structures per configuration, the number of dets per structure,
c.....and the det numbers per stucture
c.....
_IFN(atmol)
INCLUDE(../m4/common/iofile)
_ENDIF
INCLUDE(common/c8_16vb)
      dimension iscr(lword)
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
c
c     read number of determinants per structure
c
      if (nstruc.gt.lword) call corfait(nstruc,lword,'in getpa2/1')
_IF(atmol)
      read (25) (iscr(i),i=1,nstruc)
_ELSE
      call readis(iscr,nstruc,num8)
_ENDIF
c
c     read determinant-numbers
c
      ndettot = 0
      do 10 i=1,nstruc
         ndettot = ndettot + iscr(i)
10    continue
      if (ndettot.gt.lword) call corfait(ndettot,lword,'in getpa2/2')
_IF(atmol)
      read (25) (iscr(i),i=nstruc+1,nstruc+ndettot)
_ELSE
      call readis(iscr(nstruc+1),ndettot,num8)
_ENDIF
c
c     read coefficients / dummy here
c
capollo nipw() is the number of integers per word (2 for apollo)
_IF(atmol)
      read (25) (iscr(i),i=nstruc+ndettot+1,nstruc+ndettot +
     &                                           ndettot*nipw())
_ELSE
      call readis(iscr(nstruc+ndettot+1),ndettot*nipw(),num8)
_ENDIF
c
c                read number of determinants per configuration
c
_IF(atmol)
      read(25) (iscr(i),i=nstruc+ndettot+1,nstruc+ndettot+nconf)
_ELSE
      call readis(iscr(nstruc+ndettot+1),nconf,num8)
_ENDIF
c
c                 read number of structures per configuration
c
_IF(atmol)
      read(25) (iscr(i),i=nstruc+ndettot+nconf+1,nstruc+ndettot+2*nconf)
_ELSE
      call readis(iscr(nstruc+ndettot+nconf+1),nconf,num8)
_ENDIF
c
c              read packed slater determinants
c
      npack = 0
capollo nipw() (64 bits per packed word)
      do 66 i=1,nconf
         npack = npack + (((iscr(nstruc+ndettot+i)*nelec-1)/(64/n8_16)
     &           )+1) * nipw()
66    continue
      nn = nstruc+nconf+nconf+npack+ndettot
      if (nn.gt.lword) call corfait(nn,lword,'in getpa2/3')
      nf = nstruc+nconf+nconf+ndettot
_IF(atmol)
      read(25) (iscr(i),i=nf+1,nf+npack)
_ELSE
      call readis(iscr(nf+1),npack,num8)
_ENDIF
      return
      end
      subroutine hamiltt(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &                  idetps,coeff,ig,igroup,ndetps,
     &                  trani,tranj,hamil,overl,scr1,scr2,scr3,s,g,
     &                  supers,superh,superg,ipos,weight,
     &                  nelec,ngroup,nalfa,iortho,northo,
     &                  detmat,detov,ipar)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....in this routine the matrix representation of the hamilton
c.....operator in the configuration space is constructed.
c.....
INCLUDE(common/c8_16vb)
      dimension pacdet(*),idet(*),jdet(*),detcomb(*),dettot(*),icp(*),
     &          jcp(*),idetps(*),coeff(*),ig(5,*),igroup(5,*),
     &          ndetps(*),trani(*),tranj(*),hamil(*),overl(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),supers(*),superh(*),
     &          superg(*),ipos(*),weight(*),iortho(*),
     &          detmat(*),detov(*),ipar(*)
c
INCLUDE(common/hsinfo)
      integer iprhs,irowhs
      common/hcpu_vb/ iprhs,irowhs
      logical oskiphs,oroot
c
      logical opridet
      common/pridets/opridet
INCLUDE(../m4/common/iofile)/
_IF(parallel)
INCLUDE(common/parinf)
_IF(parstruc)
      call pg_dlbchunk(1,.false.)
_ELSE
      call pg_dlbchunk(20,.false.)
_ENDIF
      call pg_dlbreset()
      icounter = 0
      ido      = -1
_ENDIF
c.....
c.....initialise io chanals
c.....
      call brtsbh(ihbloc,ihfile)
      it  = 1
      ind = 1
      iid = 1
      idd = 0
      if ( iprinv .gt. 10 ) then
      write(iwr,'(a,I5,a)') 'Starting ',ngroup,' loops, constructing
     &matrix representation of the hamilton operator in the 
     &configuration space.'
        call flushn(6)
      end if
      do 20 i=1,ngroup
c
         oskiphs = .false.
         if (iprhs.gt.0) then
            if (i.lt.irowhs) then
               if (oroot()) print *,' Skipping group ',i
               oskiphs = .true.
            else
               call cpuwal(cpu,elapse)
               if (oroot()) write(iwr,'(a,i3,a,f20.2,a,f22.0,a)')
     1       ' Handling group ',i,' at',cpu,' cpu',elapse,' wall'
            end if
         end if          
c.....
c.....unpack 'left' determinants of group  i, also construct the
c.....transformation matrix associated with this group. in this matrix
c.....row indices refer to determinants  and column indices refer to
c.....structures.
c.....
         ndeti = igroup(1,i)
         call getdet(pacdet(it),idet,igroup(1,i),igroup(2,i),
     &               nelec,trani,ndetps(ind),idetps(iid),coeff(iid))
         it     = it  + (igroup(1,i) * nelec - 1)/(64/n8_16) + 1
         ind    = ind +  igroup(2,i)
         iid    = iid +  igroup(3,i)
         jt     = 1
         jnd    = 1
         jid    = 1
         ihamil = 1
         jdd    = 0
         call flushn(6)
         do 10 j=1,i-1
            if (oskiphs) then
               jt  = jt  + (igroup(1,j) * nelec - 1)/(64/n8_16) + 1
               jid = jid +  igroup(3,j)
               jnd = jnd +  igroup(2,j)
               call vclr(hamil(ihamil),1,igroup(2,i)*igroup(2,j))
               call vclr(overl(ihamil),1,igroup(2,i)*igroup(2,j))
               ihamil = ihamil + igroup(2,i) * igroup(2,j)
               go to 10
            end if         
_IF(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c           write(iwr,*) 'i have to do structure',ido,' am at ',icounter
            if (icounter.ne.ido) then
               jt  = jt  + (igroup(1,j) * nelec - 1)/(64/n8_16) + 1
               jid = jid +  igroup(3,j)
               jnd = jnd +  igroup(2,j)
               call vclr(hamil(ihamil),1,igroup(2,i)*igroup(2,j))
               call vclr(overl(ihamil),1,igroup(2,i)*igroup(2,j))
               ihamil = ihamil + igroup(2,i) * igroup(2,j)
               go to 10
            end if
_ENDIF
            ndetj = igroup(1,j)
c.....
c.....do the same thing for group  j  ('right' determinants)
c.....
         call flushn(6)
            call getdet(pacdet(jt),jdet,igroup(1,j),igroup(2,j),
     &                  nelec,tranj,ndetps(jnd),idetps(jid),coeff(jid))
         call flushn(6)
            jt  = jt  + (igroup(1,j) * nelec - 1)/(64/n8_16) + 1
            jnd = jnd +  igroup(2,j)
            jid = jid +  igroup(3,j)
c.....
c.....now calculate the matrix elements between determinants belonging
c.....to group  i and group  j respectively
c.....
cdebugg     print*,j
         call flushn(6)
            call matham(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,ndetj,idet,jdet)
            if (opridet) then
              do k=1,ndeti
                do l=1,ndetj
                  ijind = (k+idd)*(k+idd-1)/2+l+jdd
                  detmat(ijind) = detcomb((k-1)*ndetj+l)
                  detov(ijind) = dettot((k-1)*ndetj+l)
                end do
             end do
           end if
c.....
c.....use determinants as basis vectors to span the configuration
c.....space. the matrix elements are generated by transforming the
c.....the representation of the hamiltonian in the determinant basis
c.....             (tranj)-dagger * detcomb * trani.
c.....as the matrix-multiplication routine mxmb is aware of zeros in
c.....the right-matrix only the product is written as
c.....[ (detcomb * trani)-dagger * ((tranj)-dagger)-dagger ] - dagger =
c.....[ (detcomb * trani)-dagger * tranj ]-dagger.
c.....
c.....the overlap-matrix (between structures) is constructed entirely
c.....analogously out of the determinants of the overlap-matrices of the
c.....slater-determinants.
c.....
c.....na   :  the spacing between the column elements of the left-matrix
c.....iad  :  the spacing between the row    elements of the left-matrix
c.....nb   :  the spacing between the column elements of the rig.-matrix
c.....ibd  :  the spacing between the row    elements of the rig.-matrix
c.....nc   :  the spacing between the column elements of the product
c.....icd  :  the spacing between the row    elements of the product
c.....nar  :  the number of rows    in the left-matrix
c.....nac  :  the number of columns in the left-matrix
c.....nbc  :  the number of columns in the right-matrix
c.....
            na  = 1
            iad = igroup(1,j)
            nb  = 1
            ibd = igroup(1,i)
            nc  = 1
            icd = igroup(1,j)
            nar = igroup(1,j)
            nac = igroup(1,i)
            nbc = igroup(2,i)
c.....
c.....      hamiltonian matrix
c.....
            call vclr(scr1,1,nar*nbc)
            call mxmb (detcomb,na,iad,
     &                   trani,nb,ibd,
     &                    scr1,nc,icd,
     &                    nar,nac,nbc)
c.....
c.....      overlap matrix
c.....
            call vclr(scr2,1,nar*nbc)
            call mxmb (dettot,na,iad,
     &                  trani,nb,ibd,
     &                   scr2,nc,icd,
     &                   nar,nac,nbc)
c.....
            na  = igroup(1,j)
            iad = 1
            nb  = 1
            ibd = igroup(1,j)
            nc  = igroup(2,j)
            icd = 1
            nar = igroup(2,i)
            nac = igroup(1,j)
            nbc = igroup(2,j)
c.....
c.....      hamiltonian matrix
c.....
            call vclr(hamil(ihamil),1,nar*nbc)
            call mxmb (   scr1,na,iad,
     &                   tranj,nb,ibd,
     &           hamil(ihamil),nc,icd,
     &                    nar,nac,nbc)
c.....
c.....      overlap matrix
c.....
            call vclr(overl(ihamil),1,nar*nbc)
            call mxmb(    scr2,na,iad,
     &                   tranj,nb,ibd,
     &           overl(ihamil),nc,icd,
     &                    nar,nac,nbc)
c.....
            ihamil = ihamil + nar * nbc
c.....
         jdd = jdd + ndetj
10       continue
c.....
c.....we have arrived at the diagonal of the h-matrix. calculate the
c.....elements of group i with itself and write the results of the i-th
c.....cycle.
c.....
         if (oskiphs) then
            call vclr(hamil(ihamil),1,igroup(2,i)*(igroup(2,i)+1)/2)
            call vclr(overl(ihamil),1,igroup(2,i)*(igroup(2,i)+1)/2)
            go to 19
         end if    
_IF(parstruc)
         icounter = icounter + 1
         if (ido.lt.icounter) ido = ipg_dlbtask()
c        write(iwr,*) 'i have to do structure',ido,' am at ',icounter
         if (icounter.ne.ido) then
            call vclr(hamil(ihamil),1,igroup(2,i)*(igroup(2,i)+1)/2)
            call vclr(overl(ihamil),1,igroup(2,i)*(igroup(2,i)+1)/2)
            go to 19
         end if
_ENDIF
         call mathad(detcomb,dettot,icp,jcp,supers,superh,
     &               superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &               nelec,iortho,northo,nbody,ndeti,idet)
         if (opridet) then
           do k=1,ndeti
             do l=1,k
               ijind = (k+idd)*(k+idd-1)/2+l+idd
               detmat(ijind) = detcomb(k*(k-1)/2+l)
               detov(ijind) = dettot(k*(k-1)/2+l)
             end do
           end do
         end if
         call mult11(detcomb,hamil(ihamil),scr1,igroup(2,i),igroup(1,i),
     &              trani,tranj)
         call mult11(dettot, overl(ihamil),scr1,igroup(2,i),igroup(1,i),
     &              trani,tranj)
19       call writh(i,igroup,hamil,overl,ipos,scr1)
         idd = idd + ndeti
20    continue
c.....
c.....close down io
c.....
c     call whtps
      call ertsbh
      call clredx
      if (opridet) then
        call detbasis(pacdet,igroup,idet,detmat,detov,ipar,
     &                nelec,nalfa,ngroup,idd)
      end if
      return
      end
      subroutine hathab(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,ndetj,idet,jdet,
     &                  idetps,jdetps,trani,tranj,scr4,scr5)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....same as matham, only this one can skip determinants (use spin
c.....symmetry)
c.....
      dimension icp(*),trani(*),tranj(*),scr4(*),scr5(*),
     &          supers(*),superg(*),superh(*),ig(*),ipos(*),weight(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*),jdet(*),
     &          scr1(*),idetps(*),jdetps(*),jcp(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel) 
INCLUDE(common/parinf)
_ENDIF
      ii = 1
      jj = 1
      call vclr(scr4,1,ndetj)
      call vclr(scr5,1,ndetj)
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         if (idetps(i).eq.0) then
            do 30 j=1,ndetj
_IF(parallel)
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c          write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               jj = jj + nelec
               go to 30
            end if
_ENDIF
_ENDIF
               do 20 k=1,nelec
                  jcp(k) = jdet(k+jj-1)
                  icp(k) = idet(k+ii-1)
20             continue
c.....
c.....         find the blocking structure of this matrix element
c.....
               call symblo(icp,jcp,nelec,nalfa,iortho,northo,
     &                                  supers,ipos,dprod,
     &                                  nbody,ig,s)
c             call sillyy(nelec,nalfa,ig,nbody,icp,jcp,Q(iscopy))
               if (nsing.le.2) then
                  call matre3(detcomb,dettot,icp,jcp,
     &                        supers,superh,superg,
     &                        ig,nbody,nelec,
     &                        ipos,weight,nalfa,dprod,
     &                        scr1,scr2,scr3,s,g)
                  scr4(j) = scr4(j) + detcomb * trani(i)
                  scr5(j) = scr5(j) + dettot  * trani(i)
		  if (jdetps(j).ne.0) then
                     do 50 k=1,ndeti
                        i2 = k
50                   if (idetps(k).eq.i) goto 52
52                   j2 = jdetps(j)
                  else
                     do 54 k=1,ndetj
                        j2 = k
54                   if (jdetps(k).eq.j) goto 55
55                   do 56 k=1,ndeti
                        i2 = k
56                   if (idetps(k).eq.i) goto 58
58                   continue
                  end if
		  scr4(j2) = scr4(j2) + detcomb * trani(i2)
                  scr5(j2) = scr5(j2) + dettot  * trani(i2)
               end if
               jj = jj + nelec
30          continue
         end if
         ii = ii + nelec
         jj = 1
40    continue
      detcomb = ddot(ndetj,scr4,1,tranj,1)
      dettot  = ddot(ndetj,scr5,1,tranj,1)
      return
      end
      subroutine hathad(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,idet,
     &                  trani,scr4,scr5)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....this does the same thing as matham, but now the matrix-elements
c.....are being calculated for a group of determinants among themselves
c.....(i.e. we have arrived at the diagonal)
c.....
      dimension icp(*),jcp(*),supers(*),
     &          superg(*),superh(*),ig(*),ipos(*),weight(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*),trani(*),
     &          scr4(*),scr5(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel)  
INCLUDE(common/parinf)
_ENDIF
      ii = 1
      jj = 1
      call vclr(scr4,1,ndeti)
      call vclr(scr5,1,ndeti)
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         do 30 j=1,i
_IF(parallel)  
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c         write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               detcomb = 0.0d0
               dettot = 0.0d0
               jj = jj + nelec
               go to 30
            end if
_ENDIF
_ENDIF
            do 20 k=1,nelec
               jcp(k) = idet(k+jj-1)
               icp(k) = idet(k+ii-1)
20          continue
c.....
c.....      find the blocking structure of this matrix element and
c.....      determine its consequences for the calculation.
c.....
            call symblo(icp,jcp,nelec,nalfa,iortho,northo,supers,
     &                                         ipos,dprod,nbody,ig,s)
c           call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            if (nsing.le.2) then
               call matre3(detcomb,dettot,icp,jcp,
     &                     supers,superh,superg,
     &                     ig,nbody,nelec,
     &                     ipos,weight,nalfa,dprod,
     &                     scr1,scr2,scr3,s,g)
c               call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
                scr4(j) = scr4(j) + detcomb * trani(i)
                scr5(j) = scr5(j) + dettot  * trani(i)
                scr4(i) = scr4(i) + detcomb * trani(j)
                scr5(i) = scr5(i) + dettot  * trani(j)
            end if
            it = it + 1
            jj = jj + nelec
30      continue
c.....  last additions were wrong, correct
        scr4(i) = scr4(i) - detcomb * trani(i)
        scr5(i) = scr5(i) - dettot  * trani(i)
        ii = ii + nelec
        jj = 1
40    continue
      detcomb = ddot(ndeti,scr4,1,trani,1)
      dettot  = ddot(ndeti,scr5,1,trani,1)
      return
      end
      subroutine hatham(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,ndetj,idet,jdet,
     &                  trani,tranj,scr4,scr5)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....in this routine the matrix elements between the determinants as
c.....defined in idet and jdet are calculated and at the same time
c.....are processed to yield the total matrix-element
c.....
      dimension icp(*),jcp(*),supers(*),
     &          superg(*),superh(*),ig(*),ipos(*),weight(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*),jdet(*),
     &          trani(*),tranj(*),scr4(*),scr5(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      ii = 1
      jj = 1
      call vclr(scr4,1,ndetj)
      call vclr(scr5,1,ndetj)
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         do 30 j=1,ndetj
_IF(parallel)
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c         write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               jj = jj + nelec
               go to 30
            end if
_ENDIF
_ENDIF
            do 20 k=1,nelec
               jcp(k) = jdet(k+jj-1)
               icp(k) = idet(k+ii-1)
20          continue
c.....
c.....      find the blocking structure of this matrix element
c.....
            call symblo(icp,jcp,nelec,nalfa,iortho,northo,supers,
     &                                          ipos,dprod,nbody,ig,s)
c           call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            if (nsing.le.2) then
               call matre3(detcomb,dettot,icp,jcp,
     &                     supers,superh,superg,
     &                     ig,nbody,nelec,
     &                     ipos,weight,nalfa,dprod,
     &                     scr1,scr2,scr3,s,g)
c              call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
               scr4(j) = scr4(j) + detcomb * trani(i)
               scr5(j) = scr5(j) + dettot * trani(i)
            end if
            it = it + 1
            jj = jj + nelec
30       continue
         ii = ii + nelec
         jj = 1
40    continue
      detcomb = ddot(ndetj,scr4,1,tranj,1)
      dettot  = ddot(ndetj,scr5,1,tranj,1)
      return
      end
      subroutine onlys(dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,ndetj,idet,jdet,
     &                  trani,tranj,scr4,scr5)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....in this routine the s-matrix elements between the determinants as
c.....defined in idet and jdet are calculated and at the same time
c.....are processed to yield the total matrix-element
c.....
      dimension icp(*),jcp(*),supers(*),
     &          superg(*),superh(*),ig(*),ipos(*),weight(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*),jdet(*),
     &          trani(*),tranj(*),scr4(*),scr5(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      ii = 1
      jj = 1
      call vclr(scr5,1,ndetj)
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         do 30 j=1,ndetj
_IF(parallel)
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c         write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               jj = jj + nelec
               go to 30
            end if
_ENDIF
_ENDIF
            do 20 k=1,nelec
               jcp(k) = jdet(k+jj-1)
               icp(k) = idet(k+ii-1)
20          continue
c.....
c.....      find the blocking structure of this matrix element
c.....
            call symblo(icp,jcp,nelec,nalfa,iortho,northo,supers,
     &                                          ipos,dprod,nbody,ig,s)
c           call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            if (nsing.le.2) then
               if (nsing.eq.0) then
                  dettot=dprod
               else
                  dettot=0.0d0
               endif
               scr5(j) = scr5(j) + dettot * trani(i)
            end if
            it = it + 1
            jj = jj + nelec
30       continue
         ii = ii + nelec
         jj = 1
40    continue
      dettot  = ddot(ndetj,scr5,1,tranj,1)
      return
      end
      subroutine igini(iq)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....this is somewhat stupid
      dimension iq(*)
      iq(1) = 0
      iq(2) = 0
      iq(3) = 0
      iq(4) = 0
      iq(5) = 1
      return
      end
      subroutine match(ibdet,icp,nb,nelec,nalfa,ndetps,idetps)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....find out if,for every determinant a determinant exists that
c.....has exactly the reverse spin-assignment,e.g. :
c..... (123)alpha (456)beta  <==> (456)alpha (123) beta,  if so
c.....fill idetps with the mapping, ndetps contains 1 if so, 0 if not
      dimension ibdet(nelec,nb),idetps(nb),ndetps(1),icp(nelec)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
      if (nalfa.ne.nelec-nalfa.or.nospin.or.nb.eq.1) then
c.....   forget it
         ndetps(1) = 0
         return
      end if
      call izero(nb,idetps,1)
      do 30 i=1,nb-1
c.....   note that it is essential that no det occurs twice !
         if (idetps(i).ne.0) goto 30
         ndetps(1) = 0
         do 20 j=i+1,nb
            ipar = isameb(ibdet(1,i),ibdet(1,j),icp,nelec,nalfa)
            if (ipar.eq.1) then
               ndetps(1) = 1
               idetps(j) = i
            end if
20       continue
         if (ndetps(1).eq.0) then
c...        too bad
            return
         end if
30    continue
      return
      end
      subroutine mathad(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,idet)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....this does the same thing as matham, but now the matrix-elements
c.....are being calculated for a group of determinants among themselves
c.....(i.e. we have arrived at the diagonal)
c.....
      dimension detcomb(*),dettot(*),icp(*),jcp(*),supers(*),
     &          superg(*),superh(*),ig(*),ipos(*),weight(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      it = 1
      ii = 1
      jj = 1
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         do 30 j=1,i
_IF(parallel)
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c         write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               jj = jj + nelec
               detcomb(it) = 0.0d0
               dettot(it)  = 0.0d0
               it = it + 1
               go to 30
            end if
_ENDIF
_ENDIF
            do 20 k=1,nelec
               jcp(k) = idet(k+jj-1)
               icp(k) = idet(k+ii-1)
20          continue
c.....
c.....      find the blocking structure of this matrix element and
c.....      determine its consequences for the calculation.
c.....
            call symblo(icp,jcp,nelec,nalfa,iortho,northo,supers,
     &                                         ipos,dprod,nbody,ig,s)
c            call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            if (nsing.gt.2) then
               detcomb(it) = 0.0d0
               dettot(it)  = 0.0d0
            else
               call matre3(detcomb(it),dettot(it),icp,jcp,
     &                     supers,superh,superg,
     &                     ig,nbody,nelec,
     &                     ipos,weight,nalfa,dprod,
     &                     scr1,scr2,scr3,s,g)
c               call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            end if
            it = it + 1
            jj = jj + nelec
30      continue
        ii = ii + nelec
        jj = 1
40    continue
      return
      end
      subroutine matham(detcomb,dettot,icp,jcp,supers,superh,
     &                  superg,ig,ipos,weight,nalfa,scr1,scr2,scr3,s,g,
     &                  nelec,iortho,northo,nbody,ndeti,ndetj,idet,jdet)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.....
c.....in this routine the matrix elements between the determinants as
c.....defined in idet and jdet are acquired and put in the rectangular
c.....matrix detcomb. later on this matrix will be used in the
c.....construction of the matrix representation of the hamiltonian on a
c.....determinant basis.
c.....
      dimension detcomb(*),dettot(*),icp(*),jcp(*),supers(*),
     &          superg(*),superh(*),ig(*),ipos(*),weight(*),scr1(*),
     &          scr2(*),scr3(*),s(*),g(*),iortho(*),idet(*),jdet(*)
INCLUDE(common/vblimit)
INCLUDE(../m4/common/iofile)
_IF(parallel)
INCLUDE(common/parinf)
_ENDIF
      it = 1
      ii = 1
      jj = 1
      do 40 i=1,ndeti
c.....
c.....   copies of the determinants must be made because of the pivot
c.....   search that can cause parity changes
c.....
         do 30 j=1,ndetj
_IF(parallel)
_IFN(parstruc)
            icounter = icounter + 1
            if (ido.lt.icounter) ido = ipg_dlbtask()
c         write(iwr,*) 'i have to do determinant',ido,' am at ',icounter
            if (ido.ne.icounter) then
               jj = jj + nelec
               detcomb(it) = 0.0d0
               dettot(it)  = 0.0d0
               it = it + 1
               go to 30
            end if
_ENDIF
_ENDIF
            do 20 k=1,nelec
               jcp(k) = jdet(k+jj-1)
               icp(k) = idet(k+ii-1)
20          continue
c.....
c.....      find the blocking structure of this matrix element
c.....
            call symblo(icp,jcp,nelec,nalfa,iortho,northo,supers,
     &                                          ipos,dprod,nbody,ig,s)
c           call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            if (nsing.gt.2) then
               detcomb(it) = 0.0d0
               dettot(it)  = 0.0d0
            else
               call matre3(detcomb(it),dettot(it),icp,jcp,
     &                     supers,superh,superg,
     &                     ig,nbody,nelec,
     &                     ipos,weight,nalfa,dprod,
     &                     scr1,scr2,scr3,s,g)
c              call sillyy(nelec,nalfa,ig,nbody,icp,jcp,supers)
            end if
            it = it + 1
            jj = jj + nelec
30       continue
         ii = ii + nelec
         jj = 1
40    continue
      return
      end
      subroutine matre3(value,det,ir,ic,supers,superh,superg,ig,nblock,
     &          nelec,ipos,w,nalfa,dprod,scr1,scr2,scr3,s,g)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the matrix element is calculated
c.....three main cases are distinguished : nsing=0,1,2, where nsing
c.....denotes the number of singularities in the s-matrix.
c.....in singular cases the number of rectangular blocks is determinant
c.....for the algorithm chosen. the (symbolic) blocking structure of
c.....the matrix element is contained in the common-block /vblimit
c.....nsing is the number of singularities. i/j/k/lfix indicate rows
c.....(i/j) and columns (k/l) of zeros in the overlap-matrix. is1 to is4
c.....are the block-numbers of the singular (rectangular) blocks in the
c.....overlap-matrix. nrectan is the number of them. ialfa is the number
c.....of blocks corresponding to electrons with alpha spin.
c.....the array ig(5,nblock) contains
c.....
c.....       ig(1,i) :   rows in i-th block
c.....       ig(2,i) :   cols in i-th block
c.....       ig(3,i) :  first row-number of i-th block
c.....       ig(4,i) :  first col-number of i-th block
c.....       ig(5,i) :  adress of the first element of block  i in s
c.....
      logical equal
c.....
c.....dimension arrays : = nelec
c.....
      dimension ir(*),ic(*)
c.....
c.....                   = nelec**2
c.....
      dimension s(*)
c.....
c.....                   = 5 * nblock
c.....
      dimension ig(5,*)
c.....
c.....                   = (nelec*(nelec-1)/2)**2
c.....
      dimension ipos(*),w(*),scr1(*),scr2(*),scr3(*),g(*)
c.....
c.....integral supermatrices
c.....
      dimension superg(*),superh(*),supers(*)
      external iparity
c.....
INCLUDE(common/vblimit)
      parameter (n35=35)
      common /stats_vb/ noccurr(n35)
INCLUDE(../m4/common/iofile)
c.....
c     for debugging, determinants are give via symblo
c     in common /detcop/
_IF(debug)
      k2 = nelec**2+1
      call matre2(valu2,supers,superh,superg,
     &            w,w(k2),nelec,nalfa,
     &            scr1,scr3,g,det2)
_ENDIF
c      value = valu2
c      det   = det2
c      return
c.....
c      write(iwr,'(A3,20I3)') 'ir ',(ir(j),j=1,nelec)
c      write(iwr,'(A3,20I3)') 'ic ',(ic(j),j=1,nelec)
c      do i=1,nalfa
c        write(iwr,'(10F12.8)') (s(j+(i-1)*nalfa),j=1,nalfa)
c      end do
c      print *
c      do i=1,nelec-nalfa
c        write(iwr,'(10F12.8)') (s(j+(i-1)*nelec-nalfa),j=1,nelec-nalfa)
c      end do 
c      print *
      value  = 0.0d0
      det    = 0.0d0
      if (nsing.eq.0) then
c.....
         icase = 1
c.....
c
c        x x x . . . . . . . . .
c        x x x . . . . . . . . .
c        x x x . . . . . . . . .
c        . . . x x x . . . . . .
c        . . . x x x . . . . . .
c        . . . x x x . . . . . .
c        . . . . . . x x x . . .
c        . . . . . . x x x . . .
c        . . . . . . x x x . . .
c        . . . . . . . . . x x x
c        . . . . . . . . . x x x
c        . . . . . . . . . x x x
c
c.....
c.....   no singularities, therefore a straightforward calculation
c.....
         det = dprod
c.....
c.....   one-electron contribution per block
c.....
         do 10 i=1,nblock
            call cik( s(ig(5,i)), w(ig(5,i)), scr1,scr2,ig(1,i))
10       continue
         nt = ig(5,nblock) + ig(1,nblock) * ig(1,nblock)
         n1 = nt - 1
c.....
c.....   two electron contribution involving one block at a time
c.....
         do 20 i=1,nblock
            call cikjl(w(nt),ig(1,i),w(ig(5,i)))
            ii  = ig(1,i)*(ig(1,i)-1)/2
            nt  = nt + ii * ii
20       continue
         n2a = nt - n1  - 1
          ntw=nt
         call wmix(w(nt),w,nblock,ig,n2b) !set cofactors in w at nt->
         nt  = nt  + n2b - 1
         n2  = n2a + n2b
c.....
c.....   now determine integral addresses
c.....
         call pik(ipos,ir,ic,ig,nblock)
         n   = n1 + 1
         call pikjl(ipos(n),ipos(n2+n),ir,ic,ig,nblock)
c.....
c.....   mixed contributions
c.....
         n = n + n2a
         call gmix(ipos(n),ipos(n2+n),ir,ic,ig,nblock,ialfa)
         ! ipos at (it) to g
         call gather(n1         ,g      ,superh,ipos      )
         call gather(2*(n2a+n2b),g(n1+1),superg,ipos(n1+1))
         ! subtract
         call subvec(g(n1+1),g(n1+1),g(nt+1),nt-n1)
         value = ddot(nt,g,1,w,1) * dprod
c.....
      else if (nsing.eq.1) then
c.....
         if (nrectan.eq.0) then
c.....
            if (ifix.eq.0) then
c.....
               icase = 2
c.....
c
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
               call c00( s(ig(5,is0)),ig(1,is0),w )
               n1 = ig(1,is0) * ig(1,is0)
               n  = n1 + 1
c.....
c.....         two electron part involving "true" second order cofactors
c.....         of the singular block
c.....
               call c0000(ig(1,is0),w(n),scr1,scr2,scr3,s(ig(5,is0)))
               n2a = (ig(1,is0) * (ig(1,is0)-1) / 2)
               n2a = n2a * n2a
               n   = n + n2a
c.....
c.....         mixed contributions involving the singular block always
c.....
               iscr = 1
               do 60 i=1,nblock
                  if(i.ne.is0) then
                     call cik(s(ig(5,i)),scr1(iscr),scr2,scr3,ig(1,i))
                     iscr = iscr + ig(1,i) * ig(1,i)
                  end if
60             continue
               call wmix0(w(n),w,n1,scr1,iscr-1,n2b)
               nt = n + n2b
               call p00( ipos,ir(ig(3,is0)),ic(ig(4,is0)),ig(1,is0) )
               n  = n1 + 1
               call p0000(ipos(n),ipos(nt),ir(ig(3,is0)),ic(ig(4,is0)),
     &                                                        ig(1,is0))
               call gmix0(ipos(n+n2a),ipos(nt+n2a),ir,ic,ig,nblock,is0,
     &                                                            ialfa)
               call gather(n1         ,g   ,superh,ipos   )
               call gather(2*(n2a+n2b),g(n),superg,ipos(n))
               call subvec(g(n),g(n),g(nt),n2a+n2b)
               value  = ddot(nt-1,g,1,w,1) * dprod
c.....
            else
c.....
               icase = 3
c.....
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . . . . . . . . . .
c              . . . . x x . . . . . .
c              . . . . x x . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
c.....         for one-electron part just one integral
c.....
               ipar = iparity(ig)
               irf     = ir(ifix)
               icf     = ic(kfix)
               ipos(1) = max(irf,icf)*(max(irf,icf)-1)/2+min(irf,icf)
               w(1)    = 1.0d0
c.....
c.....         second order cofactors are first order really
c.....
               do 50 i=1,nblock
                  call cik( s(ig(5,i)),w(ig(5,i)+1),scr1,scr2,ig(1,i))
50             continue
c.....
               n = ig(5,nblock) + ig(1,nblock) * ig(1,nblock)
               call pik00(ipos(2),ipos(n+1),ir,ic,ifix,kfix,ig,
     &                                               nblock,nalfa,ialfa)
               g(1) = superh(ipos(1))
               call gather(2*(n-1),g(2),superg,ipos(2))
               call subvec(g(2),g(2),g(n+1),n-1)
               value = ddot(n,g,1,w,1) * dprod * ipar
            end if
c.....
         else if(nrectan.eq.1) then
c.....
            if (ifix.ne.0) then
c.....
               icase = 4
c.....
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . x x x . . . . . .
c              . . . x x x . . . . . .
c              . . . x x x . . . . . .
c              . . . . . . . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
               ipar = iparity(ig)
c.....
c.....         one electron part involves those choices of k that make
c.....         the singular block become square
c.....
               call c0k( s(ig(5,is1)),ig(2,is1),w )
               n1 = ig(2,is1)
               n  = n1 + 1
c.....
c.....         two electron part involving "true" second order cofactors
c.....         of the singular block
c.....
               call c0kjl(ig(2,is1),w(n),scr1,scr2,scr3,s(ig(5,is1)))
               n2a = (ig(2,is1) * (ig(2,is1)-1) / 2) * ig(1,is1)
               n   = n + n2a
c.....
c.....         mixed contributions involving the rectangle always
c.....
               iscr = 1
               do 67 i=1,nblock
                  if(i.ne.is1) then
                     call cik(s(ig(5,i)),scr1(iscr),scr2,scr3,ig(1,i))
                     iscr = iscr + ig(1,i) * ig(1,i)
                  end if
67             continue
               call wmix0(w(n),w,ig(2,is1),scr1,iscr-1,n2b)
               nt = n + n2b
               call p0k( ipos,ir(ifix),ic(ig(4,is1)),ig(2,is1) )
               n  = n1 + 1
               call p0kjl(ipos(n),ipos(nt),ir(ifix),ir(ig(3,is1)),
     &                                             ic(ig(4,is1)),ig,is1)
               call gmix0k(ipos(n+n2a),ipos(nt+n2a),ir,ifix,ic,ig,
     &                                         nblock,nalfa,ialfa,is1)
               call gather(ig(2,is1)  ,g   ,superh,ipos   )
               call gather(2*(n2a+n2b),g(n),superg,ipos(n))
               call subvec(g(n),g(n),g(nt),n2a+n2b)
               value  = ddot(nt-1,g,1,w,1) * dprod * ipar
c.....
            else
c.....
               icase = 5
c.....
c
c              x x . . . . . . . . . .
c              x x . . . . . . . . . .
c              x x . . . . . . . . . .
c              . . x x x . . . . . . .
c              . . x x x . . . . . . .
c              . . x x x . . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
               ipar = iparity(ig)
c.....
c.....         one electron part involves those choices of i that make
c.....         the singular block become square
c.....
               call ci0t( s(ig(5,is1)),ig(1,is1),w )
               n1 = ig(1,is1)
               n  = n1 + 1
c.....
c.....         two electron part involving "true" second order cofactors
c.....         of the singular block
c.....
               call ci0jl(ig(1,is1),w(n),scr1,scr2,scr3,s(ig(5,is1)))
               n2a = (ig(1,is1) * (ig(1,is1)-1) / 2) * ig(2,is1)
               n   = n + n2a
c.....
c.....         mixed contributions involving the rectangle always
c.....
               iscr = 1
               do 70 i=1,nblock
                  if(i.ne.is1) then
                     call cik(s(ig(5,i)),scr1(iscr),scr2,scr3,ig(1,i))
                     iscr = iscr + ig(1,i) * ig(1,i)
                  end if
70             continue
               call wmix0(w(n),w,ig(1,is1),scr1,iscr-1,n2b)
               nt = n + n2b
               call pi0( ipos,ir(ig(3,is1)),ig(1,is1),ic(kfix))
               n  = n1 + 1
               call pi0jl(ipos(n),ipos(nt),ic(kfix),ir(ig(3,is1)),
     &                                             ic(ig(4,is1)),ig,is1)
               call gmixi0(ipos(n+n2a),ipos(nt+n2a),ir,ic,kfix,ig,
     &                                         nblock,nalfa,ialfa,is1)
               call gather(ig(1,is1)  ,g   ,superh,ipos   )
               call gather(2*(n2a+n2b),g(n),superg,ipos(n))
               call subvec(g(n),g(n),g(nt),n2a+n2b)
               value  = ddot(nt-1,g,1,w,1) * dprod * ipar
c.....
            end if
c.....
         else
c.....
c.....      as far as the singly singular cases are concerned the only
c.....      possibilty left is two rectangular blocks. make sure is1
c.....      always refers to the one that has more rows than columns
c.....
c
c           x x x . . . . . . . . .         x x x . . . . . . . . .
c           x x x . . . . . . . . .         x x x . . . . . . . . .
c           x x x . . . . . . . . .         . . . x x x . . . . . .
c           x x x . . . . . . . . .         . . . x x x . . . . . .
c           . . . x x x . . . . . .         . . . x x x . . . . . .
c           . . . x x x . . . . . .         . . . x x x . . . . . .
c           . . . . . . x x x . . .    or   . . . . . . x x x . . .
c           . . . . . . x x x . . .         . . . . . . x x x . . .
c           . . . . . . x x x . . .         . . . . . . x x x . . .
c           . . . . . . . . . x x x         . . . . . . . . . x x x
c           . . . . . . . . . x x x         . . . . . . . . . x x x
c           . . . . . . . . . x x x         . . . . . . . . . x x x
c
c.....
            icase = 6
c.....
            if (ig(1,is1).lt.ig(2,is1)) then
               iii = is1
               is1 = is2
               is2 = iii
               end if
            ipar = iparity(ig)
c.....
c.....      one electron part
c.....
            ist = ig(1,is1) + ig(2,is2) + 1
            call ci0t(s(ig(5,is1)),ig(1,is1),w             )
            call c0k(s(ig(5,is2)),ig(2,is2),w(ig(1,is1)+1))
            call wmix0(w(ist),w(ig(1,is1)+1),ig(2,is2),
     &                        w             ,ig(1,is1),n1)
c.....
c.....      second order cofactors of first rectangle
c.....
            n    = ist + n1
            n2a  = ig(1,is1) * (ig(1,is1)-1) * ig(2,is1) / 2
            iscr = n + n2a * ig(2,is2)
            call ci0jl(ig(1,is1),w(iscr),scr1,scr2,scr3,s(ig(5,is1)))
            call wmix0(w(n),w(ig(1,is1)+1),ig(2,is2),
     &                      w(iscr       ),n2a,nn)
            n2a = n2a * ig(2,is2)
            n   = n + nn
c.....
c.....      second order cofactors of second rectangle
c.....
            n2b  = ig(2,is2) * (ig(2,is2)-1) * ig(1,is2) / 2
            iscr = n + n2b * ig(1,is1)
            call c0kjl(ig(2,is2),w(iscr),scr1,scr2,scr3,s(ig(5,is2)))
            call wmix0(w(n),w             ,ig(1,is1),
     &                      w(iscr)       ,n2b,nn)
            n2b = n2b * ig(1,is1)
            n   = n + nn
c.....
c.....      mixed part, involving the two rectangles always
c.....
            iscr = 1
            do 80 i=1,nblock
               if(i.ne.is1.and.i.ne.is2) then
                  call cik( s(ig(5,i)),scr1(iscr),scr2,scr3,ig(1,i))
                  iscr = iscr + ig(1,i) * ig(1,i)
               end if
80          continue
            call wmix0(w(n),w(ist),n1,scr1,iscr-1,n2c)
            nt = n1 + n2a + n2b + n2c
            call pab(ipos,ir(ig(3,is1)),ig(1,is1),
     &                    ic(ig(4,is2)),ig(2,is2))
            call pabaa(ipos(n1+1),ipos(nt+1),ir(ig(3,is1)),ig(1,is1),
     &                         ic(ig(4,is1)),ir(ig(3,is2)),
     &                         ic(ig(4,is2)),              ig(2,is2))
            n = n1 + n2a + n2b + 1
            call gmixab(ipos(n),ipos(nt+n2a+n2b+1),ir,ic,ig,nblock,
     &                                                  ialfa,is1,is2)
            call gather(n1,g,superh,ipos)
            call gather(2*(nt-n1),g(n1+1),superg,ipos(n1+1))
            call subvec(g(n1+1),g(n1+1),g(nt+1),nt-n1)
            value = ddot(nt,g,1,w(ist),1) * dprod * ipar
c.....
         end if
c.....
      else
c.....
c.....   two singularities in the overlap matrix
c.....
         if (nrectan.eq.0) then
c.....
            if (ifix.eq.0) then
c.....
               if (is02.ne.0) then
c.....
                  icase = 7
c.....
c
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 . . . 0 0 0 . . . . . .
c                 . . . 0 0 0 . . . . . .
c                 . . . 0 0 0 . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . 0 0 0
c                 . . . . . . . . . 0 0 0
c                 . . . . . . . . . 0 0 0
c
c.....
                  na = ig(1,is0 ) * ig(1,is0 )
                  nb = ig(1,is02) * ig(1,is02)
                  nt = na * nb
                  call c00(s(ig(5,is0 )),ig(1,is0 ),scr1)
                  call c00(s(ig(5,is02)),ig(1,is02),scr2)
                  call wmix0(w,scr2,nb,scr1,na,nt)
                  equal = .false.
                  if (is0.le.ialfa.and.is02.le.ialfa.or.
     &                is0.gt.ialfa.and.is02.gt.ialfa) equal = .true.
                  call gmix00(ipos,ipos(nt+1),ir,ic,ig,is0,is02,equal)
                  if (equal) then
                     call gather(2*nt,g,superg,ipos)
                     call subvec(g,g,g(nt+1),nt)
                  else
                     call gather(nt,g,superg,ipos)
                  end if
                  value = ddot(nt,g,1,w,1) * dprod
               else
c.....
                  icase = 8
c.....
c
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 . . . 2 2 2 . . . . . .
c                 . . . 2 2 2 . . . . . .
c                 . . . 2 2 2 . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
                  nt = ig(1,is0) * (ig(1,is0)-1) / 2
                  nt = nt * nt
                  call c2222(ig(1,is0),w,s(ig(5,is0)),scr1,scr2)
                  call p2222(ipos,ipos(nt+1),ir(ig(3,is0)),ig(1,is0),
     &                                       ic(ig(4,is0)),ig(2,is0))
                  call gather(nt*2,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod
c.....
               end if
c.....
            else if (jfix.eq.0) then
c.....
               icase = 9
c.....
c
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . 0 0 0 . . . . . .
c              . . . . . . . x x . . .
c              . . . . . . . x x . . .
c              . . . . . . . . . . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
               ipar = iparity(ig)
               call c00(s(ig(5,is0)),ig(1,is0),w)
               equal = .false.
               if((is0.le.ialfa.and.ifix.le.nalfa).or.
     &            (is0.gt.ialfa.and.ifix.gt.nalfa)) equal = .true.
               nt = ig(1,is0) * ig(2,is0)
               call g0ff(ipos,ipos(nt+1),ir(ig(3,is0)),ic(ig(4,is0)),
     &                     ig(1,is0),ir(ifix),ic(kfix),equal)
               if (equal) then
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
               else
                  call gather(nt,g,superg,ipos)
               end if
               value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
            else
c.....
               icase = 10
c.....
c
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . . x x . . . . . .
c              . . . . x x . . . . . .
c              . . . . . . . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . . x x
c              . . . . . . . . . . x x
c              . . . . . . . . . . . .
c
c.....
               ipar = iparity(ig)
c.....
c.....         just one integral
c.....
               ir1 = ir(ifix)
               ir2 = ir(jfix)
               ic1 = ic(kfix)
               ic2 = ic(lfix)
               ikjl = intpos(ir1,ic1,ir2,ic2)
               gc  = superg(ikjl)
               if ( (ifix.le.nalfa.and.jfix.le.nalfa).or.
     &              (ifix.gt.nalfa.and.jfix.gt.nalfa)   ) then
                  iljk = intpos(ir1,ic2,ir2,ic1)
                  gc= gc - superg(iljk)
               end if
               value = gc * dprod * ipar
            end if
c.....
         else if (nrectan.eq.1) then
c.....
            if (is1.eq.is0) then
c.....
               if (ig(1,is0).lt.ig(2,is0)) then
c.....
                  icase = 11
c.....
c
c                    2 2 2 2 . . . . . . . .
c                    2 2 2 2 . . . . . . . .
c                    2 2 2 2 . . . . . . . .
c                    . . . . x x . . . . . .
c                    . . . . x x . . . . . .
c                    . . . . . . . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c
c.....
                  ipar = iparity(ig)
                  nt = ig(2,is0)
                  nt = (nt-1) * nt * (nt-1) / 2
                  call c2kjl(ig(2,is0),w,scr1,scr2,s(ig(5,is0)))
                  call p0kjl(ipos,ipos(nt+1),ir(ifix),ir(ig(3,is0)),
     &                                           ic(ig(4,is0)),ig,is0)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
                  icase = 12
c.....
c
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c
c.....
                  ipar = iparity(ig)
                  nt   = ig(1,is0)
                  nt   = (nt-1) * nt * (nt-1) / 2
                  call ci2jl(ig(1,is0),w,scr1,scr2,s(ig(5,is0)))
                  call pi0jl(ipos,ipos(nt+1),ic(kfix),ir(ig(3,is0)),
     &                                           ic(ig(4,is0)),ig,is0)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
               end if
c.....
            else if ( iabs(ig(1,is1)-ig(2,is1) ).eq.1) then
c.....
               if ( ig(1,is1).gt.ig(2,is1) ) then
c.....
                  if (ifix.eq.0) then
c.....
                     icase = 13
c.....
c
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . 0 0 0
c                    . . . . . . . . . 0 0 0
c                    . . . . . . . . . 0 0 0
c
c.....
                     ipar = iparity(ig)
                     n0 = ig(1,is0) * ig(1,is0)
                     call ci0t(s(ig(5,is1)),ig(1,is1),scr1)
                     call c00(s(ig(5,is0)),ig(1,is0),scr2)
                     call wmix0(w,scr1,ig(1,is1),scr2,n0,nt)
                     equal = .false.
                     if ((kfix.le.nalfa.and.is0.le.ialfa).or.
     &                   (kfix.gt.nalfa.and.is0.gt.ialfa)) equal =.true.
                     call gmix0f(ipos,ipos(nt+1),ir,ic,ig,ic(kfix),
     &                                                    is1,is0,equal)
                     if (equal) then
                        call gather(2*nt,g,superg,ipos)
                        call subvec(g,g,g(nt+1),nt)
                     else
                        call gather(nt,g,superg,ipos)
                     end if
                     value = ddot(nt,g,1,w,1) * dprod * ipar
                  else
c.....
                     icase = 14
c.....
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    x x x . . . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . x x . . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . . x x
c                    . . . . . . . . . . x x
c                    . . . . . . . . . . . .
c
c.....
                     ipar = iparity(ig)
                     equal = .false.
                     if ((kfix.le.nalfa.and.lfix.le.nalfa) .or.
     &                   (kfix.gt.nalfa.and.lfix.gt.nalfa) )equal=.true.
                     nt = ig(1,is1)
                     call ci0t(s(ig(5,is1)),nt,w)
c.....
c.....               remember that lfix < kfix , mind spin !!!!!!
c.....
                     if (ifix.lt.ig(3,is1)) then
                        icfix1 = lfix
                        icfix2 = kfix
                     else
                        icfix1 = kfix
                        icfix2 = lfix
                     end if
                     call pi000(ipos,ipos(nt+1),ir,ic,ig,is1,ifix,icfix1
     &                                                    ,icfix2,equal)
                     if (equal) then
                        call gather(2*nt,g,superg,ipos)
                        call subvec(g,g,g(nt+1),nt)
                     else
                        call gather(nt  ,g,superg,ipos)
                     end if
                     value = ddot(nt,g,1,w,1) * dprod * ipar
                  end if
c.....
               else if (kfix.eq.0) then
c.....
                  icase = 15
c.....
c
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . . . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . 0 0 0
c                 . . . . . . . . . 0 0 0
c                 . . . . . . . . . 0 0 0
c
c.....
                  ipar = iparity(ig)
                  n0 = ig(1,is0) * ig(1,is0)
                  call c0k(s(ig(5,is1)),ig(2,is1),scr1)
                  call c00(s(ig(5,is0)),ig(1,is0),scr2)
                  call wmix0(w,scr1,ig(2,is1),scr2,n0,nt)
                  equal = .false.
                  if ((ifix.le.nalfa.and.is0.le.ialfa).or.
     &                (ifix.gt.nalfa.and.is0.gt.ialfa)) equal =.true.
                  call gmixf0(ipos,ipos(nt+1),ir,ic,ig,ir(ifix),
     &                                                 is1,is0,equal)
                  if (equal) then
                     call gather(2*nt,g,superg,ipos)
                     call subvec(g,g,g(nt+1),nt)
                  else
                     call gather(nt,g,superg,ipos)
                  end if
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
                  icase = 16
c.....
c
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . . . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . . x x
c                 . . . . . . . . . . x x
c                 . . . . . . . . . . . .
c
c.....
                  ipar = iparity(ig)
                  equal = .false.
                  if ( (ifix.le.nalfa.and.jfix.le.nalfa) .or.
     &                 (ifix.gt.nalfa.and.jfix.gt.nalfa) ) equal =.true.
                  nt = ig(2,is1)
                  call c0k(s(ig(5,is1)),ig(2,is1),w)
c.....
c.....            jfix < ifix => associate the right one with is1 !!!!
c.....
                  if (kfix.lt.ig(4,is1)) then
                     irfix1 = ifix
                     irfix2 = jfix
                  else
                     irfix1 = jfix
                     irfix2 = ifix
                  end if
                  call p0k00(ipos,ipos(nt+1),ir,ic,ig,is1,irfix1,irfix2,
     &                                                       kfix,equal)
                  if (equal) then
                     call gather(2*nt,g,superg,ipos)
                     call subvec(g,g,g(nt+1),nt)
                  else
                     call gather(nt  ,g,superg,ipos)
                  end if
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
            else
c.....
               if (ig(1,is1).gt.ig(2,is1)) then
c.....
                  icase = 17
c.....
c
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 . . . x x . . . . . . .
c                 . . . x x . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
c.....            ipar = 1 - 2 * mod(kfix + lfix + ig(4,is1) + ig(2,is1)
c.....                                           + ig(4,is1) + ig(2,is1)
c.....                                           + 1 , 2) =>
c.....
                  ipar = -iparity(ig)
                  call ci0j0(ig(1,is1),s(ig(5,is1)),w,scr1,scr2)
                  nt = ig(1,is1) * (ig(1,is1)-1) / 2
                  call pa00b(ipos,ipos(nt+1),ir,ic,ig,is1,kfix,lfix)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
                  icase = 18
c.....
c
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 . . . . . . . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . . . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
c.....            see comment on parity above
c.....
                  ipar = - iparity(ig)
                  call c0k0l(ig(2,is1),s(ig(5,is1)),w,scr1,scr2)
                  nt = ig(2,is1) * (ig(2,is1)-1) / 2
                  call p0k0l(ipos,ipos(nt+1),ir,ic,ig,is1,ifix,jfix)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
            end if
c.....
         else if (nrectan.eq.2) then
c.....
            ndelta = ig(1,is1) - ig(2,is1) + ig(1,is2) - ig(2,is2)
c.....
            if ( (ig(3,is1).le.nalfa.and.ig(3,is2).le.nalfa).or.
     &           (ig(3,is1).gt.nalfa.and.ig(3,is2).gt.nalfa)     ) then
c.....
               if ( ndelta.eq.0 ) then
c.....
                  if (ifix.eq.0) then
c.....
                   if (is0.ne.0) then
                    if ((is0.eq.is1).or.(is0.eq.is2)) then
                     if (is0.eq.is2) then
                      is2=is1
                      is1=is0
                     end if
                     if (ig(1,is0).lt.ig(2,is0)) then
c.....
                     icase = 19
c.....
c
c                    2 2 2 2 . . . . . . . .
c                    2 2 2 2 . . . . . . . .
c                    2 2 2 2 . . . . . . . .
c                    . . . . x x . . . . . .
c                    . . . . x x . . . . . .
c                    . . . . x x . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c
c.....
                     ipar = iparity(ig)
                     na = ig(2,is0)
                     na = (na-1) * na * (na-1) / 2
                     call c2kjl(ig(2,is0),scr1,scr2,scr3,s(ig(5,is0)))
                     nb = ig(1,is2)
                     call ci0t(s(ig(5,is2)),ig(1,is2),scr2)
                     call wmix0(w,scr2,nb,scr1,na,nt)
                     call pabbb(ipos,ipos(nt+1),ir,ic,is0,is2,ig)
                     call gather(2*nt,g,superg,ipos)
                     call subvec(g,g,g(nt+1),nt)
                     value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
                     else
c.....
                     icase = 20
c.....
c
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    2 2 2 . . . . . . . . .
c                    . . . x x x . . . . . .
c                    . . . x x x . . . . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . x x x . . .
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c                    . . . . . . . . . x x x
c
c.....
                     ipar = iparity(ig)
                     na   = ig(1,is0)
                     na   = (na-1) * na * (na-1) / 2
                     call ci2jl(ig(1,is0),scr1,scr2,scr3,s(ig(5,is0)))
                     nb = ig(2,is2)
                     call c0k(s(ig(5,is2)),ig(2,is2),scr2)
                     call wmix0(w,scr2,nb,scr1,na,nt)
                     call pbabb(ipos,ipos(nt+1),ir,ic,is0,is2,ig)
                     call gather(2*nt,g,superg,ipos)
                     call subvec(g,g,g(nt+1),nt)
                     value = ddot(nt,g,1,w,1) * dprod * ipar
                    end if
                    else
                        
c.....
                        icase = 21
c
c.....
c.....
c
c           x x x . . . . . . . . .         x x x . . . . . . . . .
c           x x x . . . . . . . . .         x x x . . . . . . . . .
c           x x x . . . . . . . . .         . . . x x x . . . . . .
c           x x x . . . . . . . . .         . . . x x x . . . . . .
c           . . . x x x . . . . . .         . . . x x x . . . . . .
c           . . . x x x . . . . . .         . . . x x x . . . . . .
c           . . . . . . x x x . . .    or   . . . . . . x x x . . .
c           . . . . . . x x x . . .         . . . . . . x x x . . .
c           . . . . . . x x x . . .         . . . . . . x x x . . .
c           . . . . . . . . . 0 0 0         . . . . . . . . . 0 0 0
c           . . . . . . . . . 0 0 0         . . . . . . . . . 0 0 0
c           . . . . . . . . . 0 0 0         . . . . . . . . . 0 0 0
c
c.....
            if (ig(1,is1).lt.ig(2,is1)) then
               iii = is1
               is1 = is2
               is2 = iii
            end if
            ipar = iparity(ig)
            n0 = ig(1,is0) * ig(1,is0)
            nr = ig(1,is1) * ig(2,is2)
            nt = n0 + nr
            call ci0t(s(ig(5,is1)),ig(1,is1),scr1)
            call c0k(s(ig(5,is2)),ig(2,is2),scr2)
            call wmix0(scr3,scr2,ig(2,is2),scr1,ig(1,is1),n)
            call c00(s(ig(5,is0)),ig(1,is0),scr1)
            call wmix0(w,scr3,n,scr1,n0,nt)
            equal = .false.
            if ((is1.le.ialfa.and.is0.le.ialfa).or.
     &         (is1.gt.ialfa.and.is0.gt.ialfa))equal=.true.
            call pabcd(ipos,ipos(nt+1),ir,ic,ig,is0,is0,is1,
     &                 is2,equal)
            if (equal) then
               call gather(2*nt,g,superg,ipos)
               call subvec(g,g,g(nt+1),nt)
            else
               call gather(nt,g,superg,ipos)
            end if
            value = ddot(nt,g,1,w,1) * dprod * ipar
            end if
         else
c.....
            icase = 22
c.....
c
c                   x x . . . . . . . . . .     x x x x . . . . . . . .
c                   x x . . . . . . . . . .     x x x x . . . . . . . .
c                   x x . . . . . . . . . .     . . . . x x . . . . . .
c                   x x . . . . . . . . . .     . . . . x x . . . . . .
c                   . . x x x x . . . . . .     . . . . x x . . . . . .
c                   . . x x x x . . . . . .     . . . . x x . . . . . .
c                   . . . . . . x x x . . . or  . . . . . . x x x . . .
c                   . . . . . . x x x . . .     . . . . . . x x x . . .
c                   . . . . . . x x x . . .     . . . . . . x x x . . .
c                   . . . . . . . . . x x x     . . . . . . . . . x x x
c                   . . . . . . . . . x x x     . . . . . . . . . x x x
c                   . . . . . . . . . x x x     . . . . . . . . . x x x
c
c.....
                        if (ig(1,is1).lt.ig(2,is1)) then
                           iii = is1
                           is1 = is2
                           is2 = iii
                        end if
            ipar = 1
c.....
                        na = ig(1,is1)*(ig(1,is1)-1)/2
                        nb = ig(2,is2)*(ig(2,is2)-1)/2
                        nt = na * nb
                        call ci0j0(ig(1,is1),s(ig(5,is1)),w(nt   +1),
     &                                                     scr1,scr2)
                        call c0k0l(ig(2,is2),s(ig(5,is2)),w(nt+na+1),
     &                                                     scr1,scr2)
                        call wmix0(w,w(nt+na+1),nb,w(nt+1),na,nt)
                        call pabab(ipos,ipos(nt+1),ir,ic,ig,is1,is2)
                        call gather(2*nt,g,superg,ipos)
                        call subvec(g,g,g(nt+1),nt)
                        value = ddot(nt,w,1,g,1) * dprod * ipar
c.....
                     end if
c.....
                  else
c.....
                     icase = 23
c.....
c
c                    x x x . . . . . . . . .   x x x . . . . . . . . .
c                    x x x . . . . . . . . .   x x x . . . . . . . . .
c                    x x x . . . . . . . . .   . . . x x x . . . . . .
c                    x x x . . . . . . . . .   . . . x x x . . . . . .
c                    . . . x x x . . . . . .   . . . x x x . . . . . .
c                    . . . x x x . . . . . .   . . . x x x . . . . . .
c                    . . . . . . x x . . . .   . . . . . . x x . . . .
c                    . . . . . . x x . . . .   . . . . . . x x . . . .
c                    . . . . . . . . . x x x   . . . . . . . . . x x x
c                    . . . . . . . . . x x x   . . . . . . . . . x x x
c                    . . . . . . . . . x x x   . . . . . . . . . x x x
c                    . . . . . . . . . . . .   . . . . . . . . . . . .
c
c.....
                   ipar = iparity(ig)
                   if (ig(1,is1).lt.ig(2,is1)) then
                      iii = is1
                      is1 = is2
                      is2 = iii
                   end if
                   call ci0t(s(ig(5,is1)),ig(1,is1),scr1)
                   call c0k(s(ig(5,is2)),ig(2,is2),scr2)
                   call wmix0(w,scr2,ig(2,is2),scr1,ig(1,is1),nt)
                   equal = .false.
                   if ((is1.le.ialfa.and.ifix.le.nalfa).or.
     &                 (is1.gt.ialfa.and.ifix.gt.nalfa)) equal=.true.
                   call p00ab(ipos,ipos(nt+1),ir,ic,ig,is1,is2,
     &                                              equal,ifix,kfix)
                   if (equal) then
                      if (ifix.lt.ig(3,is1)) ipar=-ipar
                      if (kfix.lt.ig(4,is2)) ipar=-ipar
                      call gather(2*nt,g,superg,ipos)
                      call subvec(g,g,g(nt+1),nt)
                   else
                      call gather(nt,g,superg,ipos)
                   end if
                   value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
               else if(ndelta.eq.2) then
c.....
                  icase = 24
c.....
c
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
                  ipar = iparity(ig)
                  call ci0t(s(ig(5,is1)),ig(1,is1),scr1)
                  call ci0t(s(ig(5,is2)),ig(1,is2),scr2)
                  call wmix0(w,scr2,ig(1,is2),scr1,ig(1,is1),nt)
                  call pa0b0(ipos,ipos(nt+1),ir,ic,ig,is1,is2,lfix,kfix,
     &                                                           .true.)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else if(ndelta.eq.-2) then
c.....
                 icase = 25
c.....
c
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 . . . . . . . . . . . .
c                 . . . x x x . . . . . .
c                 . . . x x x . . . . . .
c                 . . . . . . . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
                  ipar = iparity(ig)
                  call c0k(s(ig(5,is1)),ig(2,is1),scr1)
                  call c0k(s(ig(5,is2)),ig(2,is2),scr2)
                  call wmix0(w,scr2,ig(2,is2),scr1,ig(2,is1),nt)
                  call p0a0b(ipos,ipos(nt+1),ir,ic,ig,is1,is2,ifix,jfix,
     &                                                           .true.)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else if (ndelta.eq.1) then
c.....
                  icase = 26
c.....
c
c                 x x . . . . . . . . . .     x x x . . . . . . . . .
c                 x x . . . . . . . . . .     x x x . . . . . . . . .
c                 x x . . . . . . . . . .     . . . x x . . . . . . .
c                 x x . . . . . . . . . .     . . . x x . . . . . . .
c                 . . x x x . . . . . . .     . . . x x . . . . . . .
c                 . . x x x . . . . . . .     . . . x x . . . . . . .
c                 . . . . . . x x x . . . or  . . . . . . x x x . . .
c                 . . . . . . x x x . . .     . . . . . . x x x . . .
c                 . . . . . . x x x . . .     . . . . . . x x x . . .
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c
c.....
             if (ig(1,is1).lt.ig(2,is1)) then
                iii = is1
                is1 = is2
                is2 = iii
             end if
             ipar = -iparity(ig)
             if (kfix.lt.ig(4,is2)) ipar = -ipar
                  call ci0j0(ig(1,is1),s(ig(5,is1)),scr1,scr2,scr3)
                  na = ig(1,is1) * (ig(1,is1)-1) / 2
                  call c0k(s(ig(5,is2)),ig(2,is2),scr2)
                  call wmix0(w,scr2,ig(2,is2),scr1,na,nt)
                  call pa0ab(ipos,ipos(nt+1),ir,ic,ig,is1,is2,kfix)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
c.....               ndelta = -1
c.....
                  icase = 27
c.....
c
c                 x x . . . . . . . . . .     x x x x . . . . . . . .
c                 x x . . . . . . . . . .     x x x x . . . . . . . .
c                 x x . . . . . . . . . .     . . . . x x . . . . . .
c                 . . x x x x . . . . . .     . . . . x x . . . . . .
c                 . . x x x x . . . . . .     . . . . x x . . . . . .
c                 . . . . . . . . . . . .     . . . . . . . . . . . .
c                 . . . . . . x x x . . . or  . . . . . . x x x . . .
c                 . . . . . . x x x . . .     . . . . . . x x x . . .
c                 . . . . . . x x x . . .     . . . . . . x x x . . .
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c                 . . . . . . . . . x x x     . . . . . . . . . x x x
c
c.....
             if (ig(1,is1).lt.ig(2,is1)) then
                iii = is1
                is1 = is2
                is2 = iii
             end if
             equal=(((ifix.lt.nalfa).and.(ig(3,is1).lt.ialfa)).or.
     &              ((ifix.gt.nalfa).and.(ig(3,is1).gt.ialfa)))
             ipar = -iparity(ig)
             if (equal.and.(ifix.lt.ig(3,is1))) ipar = -ipar
                  call c0k0l(ig(2,is2),s(ig(5,is2)),scr1,scr2,scr3)
                  na = ig(2,is2) * (ig(2,is2)-1) / 2
                  call ci0t(s(ig(5,is1)),ig(1,is1),scr2)
                  call wmix0(w,scr2,ig(1,is1),scr1,na,nt)
                  call p0abb(ipos,ipos(nt+1),ir,ic,ig,is2,is1,ifix)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
            else
c.....
               if ( ndelta.eq.0 ) then
c.....
                  if (ig(1,is1).gt.ig(2,is1)) then
c.....
                     isr = is1
                     isc = is2
c.....
                  else
c.....
                     isr = is2
                     isc = is1
c.....
                  end if
c.....
                  icase = 28
c.....
c
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . . . x x x x . .
c                 . . . . . . x x x x . .
c                 . . . . . . x x x x . .
c                 . . . . . . . . . . x x
c                 . . . . . . . . . . x x
c                 . . . . . . . . . . . .
c
c.....
                  ipar = iparity(ig)
                  equal=(((ifix.le.nalfa).and.(isr.le.ialfa)).or.
     &                     ((ifix.gt.nalfa).and.(isr.gt.ialfa)))
                  if (equal) then
                     if (ifix.lt.ig(3,isr)) ipar=-ipar
                     if (kfix.lt.ig(4,isc)) ipar=-ipar
                  end if
                  call ci0t(s(ig(5,isr)),ig(1,isr),scr1)
                  call c0k(s(ig(5,isc)),ig(2,isc),scr2)
                  call wmix0(w,scr2,ig(2,isc),scr1,ig(1,isr),nt)
                  call pi00l(ipos,ir,ic,ig,isr,isc,ifix,kfix)
                  call gather(nt,g,superg,ipos)
                  value =  ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else if(ndelta.eq.2) then
c.....
c.....            equals icase = 24, but no exchange between blocks
c.....            redundant ??
                  icase = 29
c.....
c
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 x x . . . . . . . . . .
c                 . . . x x x . . . . . .
c                 . . . x x x . . . . . .
c                 . . . x x x . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x .
c                 . . . . . . . . . x x .
c                 . . . . . . . . . x x .
c
c.....
                  ipar = iparity(ig)
                  call ci0t(s(ig(5,is1)),ig(1,is1),scr1)
                  call ci0t(s(ig(5,is2)),ig(1,is2),scr2)
                  call wmix0(w,scr2,ig(1,is2),scr1,ig(1,is1),nt)
c.....
c.....            mind lfix always must be < kfix
c.....
                  call pa0b0(ipos,ipos(nt+1),ir,ic,ig,is1,is2,lfix,kfix,
     &                                                          .false.)
                  call gather(nt,g,superg,ipos)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
                  icase = 30
c.....            cf. case 25
c
c                 x x x . . . . . . . . .
c                 x x x . . . . . . . . .
c                 . . . . . . . . . . . .
c                 . . . x x x . . . . . .
c                 . . . x x x . . . . . .
c                 . . . x x x . . . . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . . . .
c
c.....
                  ipar = iparity(ig)
                  call c0k(s(ig(5,is1)),ig(2,is1),scr1)
                  call c0k(s(ig(5,is2)),ig(2,is2),scr2)
                  call wmix0(w,scr2,ig(2,is2),scr1,ig(2,is1),nt)
c.....
c.....            mind that is1<is2 and jfix<ifix
c.....
                  call p0a0b(ipos,ipos(nt+1),ir,ic,ig,is1,is2,ifix,jfix,
     &                                                          .false.)
                  call gather(nt,g,superg,ipos)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
            end if
c.....
         else if (nrectan.eq.3) then
c.....
            if (ifix.eq.0.and.kfix.eq.0) then
c.....
               idiff1 = ig(1,is1) - ig(2,is1)
               idiff2 = ig(1,is2) - ig(2,is2)
               idiff3 = ig(1,is3) - ig(2,is3)
               if (iabs(idiff1).eq.2) then
c.....
c
c                 x x . . . . . . . . . .       x x x x . . . . . . . .
c                 x x . . . . . . . . . .       x x x x . . . . . . . .
c                 x x . . . . . . . . . .       . . . . x . . . . . . .
c                 x x . . . . . . . . . .       . . . . x . . . . . . .
c                 . . x x . . . . . . . .       . . . . . x . . . . . .
c                 . . . . x x . . . . . .       . . . . . x . . . . . .
c                 . . . . . . x x x . . .  or   . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c
c.....
                  isd = is1
                  isa = is2
                  isb = is3
                  id  = idiff1
               else if(iabs(idiff2).eq.2) then
c.....
c
c                 x x . . . . . . . . . .       x . . . . . . . . . . .
c                 . . x x . . . . . . . .       x . . . . . . . . . . .
c                 . . x x . . . . . . . .       . x x x x . . . . . . .
c                 . . x x . . . . . . . .       . x x x x . . . . . . .
c                 . . x x . . . . . . . .       . . . . . x . . . . . .
c                 . . . . x x . . . . . .       . . . . . x . . . . . .
c                 . . . . . . x x x . . .  or   . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c
c.....
                  isd = is2
                  isa = is1
                  isb = is3
                  id  = idiff2
               else
c.....
c
c                 x x . . . . . . . . . .       x . . . . . . . . . . .
c                 . . x x . . . . . . . .       x . . . . . . . . . . .
c                 . . . . x x . . . . . .       . x . . . . . . . . . .
c                 . . . . x x . . . . . .       . x . . . . . . . . . .
c                 . . . . x x . . . . . .       . . x x x x . . . . . .
c                 . . . . x x . . . . . .       . . x x x x . . . . . .
c                 . . . . . . x x x . . .  or   . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . x x x . . .       . . . . . . x x x . . .
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c                 . . . . . . . . . x x x       . . . . . . . . . x x x
c
c.....
                  isd = is3
                  isa = is1
                  isb = is2
                  id  = idiff3
               end if
c.....
               if (id.gt.0) then
c.....
                  icase = 31
c.....
c
c                 x x . . . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . x x . . . . . . . .
c                 . . . . x x . . . . . .
c                 . . . . . . x x x . . . or the other two alternatives
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
                  ipar = -iparity(ig)
                  ntwo = ig(1,isd) * (ig(1,isd)-1) / 2
                  isc  = ntwo * ig(2,isa) * ig(2,isb) + 1
                  call ci0j0(ig(1,isd),s(ig(5,isd)),w(isc),scr1,scr2)
                  call c0k(s(ig(5,isa)),ig(2,isa),scr1)
                  call c0k(s(ig(5,isb)),ig(2,isb),scr2)
                  call wmix0(w(isc+ntwo),scr2,ig(2,isb),
     &                                   scr1,ig(2,isa),n)
                  call wmix0(w,w(isc+ntwo),n,w(isc),ntwo,nt)
                  call pabac(ipos,ipos(nt+1),ir,ic,ig,isd,isa,isb)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               else
c.....
                  icase = 32
c.....
c
c                 x x x x . . . . . . . .
c                 x x x x . . . . . . . .
c                 . . . . x . . . . . . .
c                 . . . . x . . . . . . .
c                 . . . . . x . . . . . .
c                 . . . . . x . . . . . .
c                 . . . . . . x x x . . . or the other two alternatives
c                 . . . . . . x x x . . .
c                 . . . . . . x x x . . .
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c                 . . . . . . . . . x x x
c
c.....
                  ipar = -iparity(ig)
                  ntwo = ig(2,isd) * (ig(2,isd)-1) / 2
                  isc  = ntwo * ig(1,isa) * ig(1,isb) + 1
                  call c0k0l(ig(2,isd),s(ig(5,isd)),w(isc),scr1,scr2)
                  call ci0t(s(ig(5,isa)),ig(1,isa),scr1)
                  call ci0t(s(ig(5,isb)),ig(1,isb),scr2)
                  call wmix0(w(isc+ntwo),scr2,ig(1,isb),
     &                                   scr1,ig(1,isa),n)
                  call wmix0(w,w(isc+ntwo),n,w(isc),ntwo,nt)
                  call pabcb(ipos,ipos(nt+1),ir,ic,ig,isd,isa,isb)
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
                  value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
               end if
c.....
            else if (ifix.eq.0) then
c.....
               ipar = iparity(ig)
               if (ig(3,is2).le.nalfa.and.ig(3,is3).gt.nalfa) then
c.....
                  equal = .false.
c.....
                  if (ig(1,is2).lt.ig(2,is2)) then
c.....
                     isr = is1
                     isc = is2
                     iso = is3
c.....
                  else
c.....
                     isr = is2
                     isc = is1
                     iso = is3
c.....
                  end if
c.....
               else if(ig(3,is2).gt.nalfa.and.ig(3,is1).le.nalfa) then
c.....
                  equal = .false.
c.....
                  if (ig(1,is2).lt.ig(2,is2)) then
c.....
                     isr = is3
                     isc = is2
                     iso = is1
c.....
                  else
c.....
                     isr = is2
                     isc = is3
                     iso = is1
c.....
                  end if
c.....
               else
c.....
                  equal = .true.
c.....
                  if (ig(1,is1).lt.ig(2,is1)) then
c.....
                     isr = is2
                     isc = is1
                     iso = is3
c.....
                  else if (ig(1,is2).lt.ig(2,is2)) then
c.....
                     isr = is1
                     isc = is2
                     iso = is3
c.....
                  else
c.....
                     isr = is1
                     isc = is3
                     iso = is2
c.....
                  end if
c.....
               end if
c.....
               icase = 33
c.....
c
c              x x . . . . . . . . . .
c              x x . . . . . . . . . .
c              x x . . . . . . . . . .
c              . . . x . . . . . . . .
c              . . . x . . . . . . . .
c              . . . . x x . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
               ipar = iparity(ig)
               if (equal.and.(kfix.lt.ig(4,isc))) ipar = -ipar
               nscr = ig(1,isr) * ig(2,isc) * ig(1,iso) + 1
               call ci0t(s(ig(5,isr)),ig(1,isr),scr1)
               call c0k(s(ig(5,isc)),ig(2,isc),scr2)
               call wmix0(w(nscr),scr2,ig(2,isc),scr1,ig(1,isr),n)
               call ci0t(s(ig(5,iso)),ig(1,iso),scr1)
               call wmix0(w,scr1,ig(1,iso),w(nscr),n,nt)
               call pa0bc(ipos,ipos(nt+1),ir,ic,ig,isr,isc,iso,kfix,
     &                                                         equal)
               if (equal) then
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
               else
                  call gather(nt,g,superg,ipos)
               end if
               value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
            else
c.....
               if (ig(3,is2).le.nalfa.and.ig(3,is3).gt.nalfa) then
c.....
                  equal = .false.
c.....
                  if (ig(1,is1).gt.ig(2,is1)) then
c.....
                     isr = is1
                     isc = is2
                     iso = is3
c.....
                  else
c.....
                     isr = is2
                     isc = is1
                     iso = is3
c.....
                  end if
c.....
               else if(ig(3,is2).gt.nalfa.and.ig(3,is1).le.nalfa) then
c.....
                  equal = .false.
c.....
                  if (ig(1,is2).gt.ig(2,is2)) then
c.....
                     isr = is2
                     isc = is3
                     iso = is1
c.....
                  else
c.....
                     isr = is3
                     isc = is2
                     iso = is1
c.....
                  end if
c.....
               else
c.....
                  equal = .true.
c.....
                  if (ig(1,is1).gt.ig(2,is1)) then
c.....
                     isr = is1
                     isc = is2
                     iso = is3
c.....
                  else if (ig(1,is2).gt.ig(2,is2)) then
c.....
                     isr = is2
                     isc = is1
                     iso = is3
c.....
                  else
c.....
                     isr = is3
                     isc = is1
                     iso = is2
c.....
                  end if
c.....
               end if
c.....
               icase = 34
c.....
c
c              x x x . . . . . . . . .
c              x x x . . . . . . . . .
c              . . . . . . . . . . . .
c              . . . x . . . . . . . .
c              . . . x . . . . . . . .
c              . . . . x x . . . . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . x x x . . .
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c              . . . . . . . . . x x x
c
c.....
               ipar = iparity(ig)
               if (equal.and.(ifix.lt.ig(3,isr))) ipar = -ipar
               nscr = ig(1,isr) * ig(2,isc) * ig(2,iso) + 1
               call ci0t(s(ig(5,isr)),ig(1,isr),scr1)
               call c0k(s(ig(5,isc)),ig(2,isc),scr2)
               call wmix0(w(nscr),scr2,ig(2,isc),scr1,ig(1,isr),n)
               call c0k(s(ig(5,iso)),ig(2,iso),scr1)
               call wmix0(w,scr1,ig(2,iso),w(nscr),n,nt)
               call pab0c(ipos,ipos(nt+1),ir,ic,ig,isr,isc,iso,ifix,
     &                                                         equal)
               if (equal) then
                  call gather(2*nt,g,superg,ipos)
                  call subvec(g,g,g(nt+1),nt)
               else
                  call gather(nt,g,superg,ipos)
               end if
               value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
            end if
c.....
         else if (nrectan.eq.4) then
c.....
            isra = 0
            isca = 0
            if (ig(1,is4).gt.ig(2,is4)) then
               isra = is4
            else
               isca = is4
            end if
            if (ig(1,is3).gt.ig(2,is3)) then
               isrb = isra
               isra = is3
            else
               iscb = isca
               isca = is3
            end if
            if (ig(1,is2).gt.ig(2,is2)) then
               isrb = isra
               isra = is2
            else
               iscb = isca
               isca = is2
            end if
            if (ig(1,is1).gt.ig(2,is1)) then
               isrb = isra
               isra = is1
            else
               iscb = isca
               isca = is1
            end if
c.....
            icase = 35
c.....
c
c           x x x x . . . . . . . .
c           x x x x . . . . . . . .
c           x x x x . . . . . . . .
c           . . . . x x . . . . . .
c           . . . . x x . . . . . .
c           . . . . x x . . . . . .
c           . . . . . . x x x x . .
c           . . . . . . x x x x . .
c           . . . . . . x x x x . .
c           . . . . . . . . . . x x
c           . . . . . . . . . . x x
c           . . . . . . . . . . x x
c
c.....
            equal = .false.
            ipar = iparity(ig)
            if (ig(3,is1).gt.nalfa.or.ig(3,is4).le.nalfa) equal = .true.
            n  = ig(1,isra) * ig(1,isrb) * ig(2,isca) * ig(2,iscb)
            call ci0t(s(ig(5,isrb)),ig(1,isrb),scr1)
            call c0k(s(ig(5,iscb)),ig(2,iscb),scr2)
            call wmix0(w(n+1),scr2,ig(2,iscb),scr1,ig(1,isrb),nb)
            call ci0t(s(ig(5,isra)),ig(1,isra),scr1)
            call c0k(s(ig(5,isca)),ig(2,isca),scr2)
            call wmix0(w(n+nb+1),scr2,ig(2,isca),scr1,ig(1,isra),na)
            call wmix0(w,w(n+1),nb,w(n+nb+1),na,nt)
            call pabcd(ipos,ipos(nt+1),ir,ic,ig,isra,isca,isrb,iscb,
     &                                                        equal)
            if (equal) then
               call gather(2*nt,g,superg,ipos)
               call subvec(g,g,g(nt+1),nt)
            else
               call gather(nt,g,superg,ipos)
            end if
            value = ddot(nt,g,1,w,1) * dprod * ipar
c.....
         else
c.....
          write(iwr,*)'########################'
          write(iwr,*)'!!!!!!!!!!!! error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(iwr,*)'########################'
          write(iwr,*)'unrecognised type of matrix element'
          call sillyy(nelec,nalfa,ig,nblock,ir,ic,supers)
         end if
c.....
      end if
c
       noccurr(icase)=noccurr(icase)+1
_IF(debug)
       if (dabs( valu2-value ).gt.1.0d-6) then
          write(iwr,100) valu2,det2,value,det
 100      format(/,' matre2:',2f24.19,/,' matre3:',2f24.19,/)
          write(iwr,110) ipar,icase
 110      format(' ipar:',i3,/,' case:',i3,/)
          call sillyy(nelec,nalfa,ig,nblock,ir,ic,supers)
          stop
       end if
_ENDIF
c
      return
      end

      subroutine icaseprint(icase)
      implicit REAL  (a-h,o-z) , integer   (i-n)
      common/icaselist/ icaselist(35)
 
      icaselist(icase) = icaselist(icase) + 1

      do i=1,35
        if ( icaselist(i) .eq. 1 ) then
          write(20,110) i
110       format(' icase:',i3)
        end if
      end do
      
      return
      end

      subroutine mtype
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....different types of matrix elements that have occurred  are
c.....presented by this routine
c.....
INCLUDE(common/turtleparam)
      common /cases/ icase(2*ncases),times(ncases),token(3*ncases)
INCLUDE(../m4/common/iofile)
      character*8 token
      totalt = 0.0d0
      itotal = 0
      timem  = 0.0d0
      imaxi  = 0
      do 10 i=1,ncases
         if (times(i).gt.timem) timem = times(i)
         if (icase(i).gt.imaxi) iiiii = i
         if (icase(i).gt.imaxi) imaxi = icase(i)
         totalt = totalt + times(i)
         itotal = itotal + icase(i)
         icase(ncases+i) = i
10    continue
      call bubble(times,ncases,-1,icase(ncases+1))
      write(iwr,66)
      write(iwr,44) itotal
      write(iwr,55) iiiii ,idint(100.d0*(imaxi  *1.0d0/itotal))
      ncp1 = ncases+1
      write(iwr,77) icase(ncp1),idint(100.d0*(icase(icase(ncp1))*1.0d0/
     c              itotal))
      write(iwr,88) idint((times(1)/totalt)*100.d0)
66    format(' ',72('^'),/,' ',25x,'matrix type analysis',/,' ')
44    format(' ','number of matrix elements considered :',i9)
55    format(' ','the most frequent type was       ',i3,'. relative occu
     &rrence  :',i3,'  ')
77    format(' ','the most time consuming type was ',i3,'. relative occu
     &rrence  :',i3,'  ')
88    format(' ',38x,'relative costs       :',i3,'  ')
      do 15 i=1,ncases
         times(i) = times(i) * totalt / timem
         icase(i) = icase(i) * itotal / imaxi
15    continue
      tstep = totalt / 20.d0
      istep = itotal / 20
      totalt = totalt + tstep
      itotal = itotal + istep
      do 30 i=1,20
         totalt = totalt - tstep
         itotal = itotal - istep
         jj = 0
         do 20 j=1,3*(ncases-12),3
            jj = jj + 1
            token(j)   = ' '
            token(j+1) = ' '
            token(j+2) = ' '
            if (times(jj              ).ge.totalt) token(j  ) = '$'
            if (icase(icase(jj+ncases)).ge.itotal) token(j+1) = '#'
20       continue
         write(iwr,11) (token(ij),ij=1,3*(ncases-12))
11       format(' ',5x,72a1)
30    continue
      write(iwr,22) (icase(ij),ij=ncases+1,2*ncases-12)
22    format(' ',4x,33i3)
      return
      end
      subroutine p00(ipos,ir,ic,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ir(*),ic(*),ipos(*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      n = 0
      do 20 k=1,ndim
         do 10 i=1,ndim
            n       = n + 1
            ipos(n) = ind(ir(i),ic(k))
10       continue
20    continue
      return
      end
      subroutine p0000(ipos,ipose,ir,ic,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*)
      common /posit/ iky(3)
      n = 0
      do 40 k=2,ndim
         do 30 l=1,k-1
            do 20 i=2,ndim
               do 10 j=1,i-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30      continue
40    continue
      return
      end
      subroutine p00ab(ipos,ipose,ir,ic,ig,isa,isb,equal,ifix,kfix)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      irfix = ir(ifix)
      icfix = ic(kfix)
      n     = 0
      if (equal) then
         do 20 l=ig(4,isb),ig(4,isb)+ig(2,isb)-1
            do 10 j=ig(3,isa),ig(3,isa)+ig(1,isa)-1
               n  = n + 1
               ipos (n) = intpos(irfix,icfix,ir(j),ic(l))
               ipose(n) = intpos(irfix,ic(l),ir(j),icfix)
10          continue
20       continue
      else
         do 40 l=ig(4,isb),ig(4,isb)+ig(2,isb)-1
            do 30 j=ig(3,isa),ig(3,isa)+ig(1,isa)-1
               n  = n + 1
               ipos (n) = intpos(irfix,icfix,ir(j),ic(l))
30          continue
40       continue
      end if
      return
      end
      subroutine p0a0b(ipos,ipose,ir,ic,ig,is1,is2,ifix,jfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      iifix = ir(ifix)
      ijfix = ir(jfix)
      if (equal) then
         do 20 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
            do 10 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
               n  = n + 1
               ipos (n) = intpos(iifix,ic(k),ijfix,ic(l))
               ipose(n) = intpos(iifix,ic(l),ijfix,ic(k))
10          continue
20       continue
      else
         do 40 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
            do 30 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
               n  = n + 1
               ipos(n) = intpos(iifix,ic(k),ijfix,ic(l))
30          continue
40       continue
      end if
      return
      end
      subroutine p0abb(ipos,ipose,ir,ic,ig,is1,is2,ifix)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n     = 0
      irfix = ir(ifix)
      do 30 j=ig(3,is2),ig(3,is2)+ig(1,is2)-1
         do 20 k=ig(4,is1)+1,ig(4,is1)+ig(2,is1)-1
            do 10 l=ig(4,is1),k-1
               n  = n + 1
               ipos (n) = intpos(irfix,ic(k),ir(j),ic(l))
               ipose(n) = intpos(irfix,ic(l),ir(j),ic(k))
10          continue
20       continue
30    continue
      return
      end
      subroutine p0k(ipos,ir,ic,nc)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ic(*),ipos(*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      do 10 k=1,nc
         ipos(k) = ind(ir,ic(k))
10    continue
      return
      end
      subroutine p0k00(ipos,ipose,ir,ic,ig,is1,jfix,ifix,kfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n     = 0
      iifix = ir(ifix)
      ijfix = ir(jfix)
      ikfix = ic(kfix)
      if (equal) then
         do 10 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            n  = n + 1
            ipos (n) = intpos(iifix,ikfix,ijfix,ic(l))
            ipose(n) = intpos(iifix,ic(l),ijfix,ikfix)
10       continue
      else
         do 20 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            n  = n + 1
            ipos(n) = intpos(iifix,ikfix,ijfix,ic(l))
20       continue
      end if
      return
      end
      subroutine p0k0l(ipos,ipose,ir,ic,ig,is1,ifix,jfix)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      iri = ir(ifix)
      irj = ir(jfix)
      n  = 0
      do 20 k=ig(4,is1)+1,ig(4,is1)+ig(2,is1)-1
         do 10 l=ig(4,is1),k-1
            n  = n + 1
            ipos (n) = intpos(iri,ic(k),irj,ic(l))
            ipose(n) = intpos(iri,ic(l),irj,ic(k))
10       continue
20    continue
      return
      end
      subroutine p0kjl(ipos,ipose,irfix,ir,ic,ig,is1)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      it = 0
      do 30 k=2,ig(2,is1)
         do 20 l=1,k-1
            do 10 j=1,ig(1,is1)
               it = it + 1
               ipos (it) = intpos(irfix,ic(k),ir(j),ic(l))
               ipose(it) = intpos(irfix,ic(l),ir(j),ic(k))
10          continue
20       continue
30    continue
      return
      end
      subroutine p2222(ipos,ipose,ir,nr,ic,nc)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*)
      common /posit/ iky(3)
      n = 0
      do 40 k=2,nc
         do 30 l=1,k-1
            do 20 i=2,nr
               do 10 j=1,i-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30      continue
40    continue
      return
      end
      subroutine pa00b(ipos,ipose,ir,ic,ig,is1,kfix,lfix)
c     subroutine pikjl(ipos,ipose,iposa,ir,ic,ig,nblock)
c
c     implicit real   (a-h,o-z) , integer   (i-n)
c
c     dimension ipos(*),ipose(*),iposa(*),ir(*),ic(*),ig(5,*)
c     ind(i,j) = max(i,j) * (max(i,j)-1) / 2 + min(i,j)
c     n = 0
c     do 30 m=1,nblock
c        ib = ig(3,m)
c        ie = ig(3,m) + ig(1,m) - 1
c        do 20 k=ib,ie
c           do 10 i=ib,ie
c              n = n + 1
c              iposa(n) = ind( ir(i),ic(k) )
c10          continue
c20       continue
c30    continue
c     n = 0
c     do 80 m=1,nblock
c        nn = ig(1,m)
c        ia = ig(5,m) + nn - 1
c        do 70 k=2,nn
c           ib = ig(5,m) - 1
c           do 60 l=1,k-1
c              do 50 i=2,nn
c                 il = ib + i
c                 ik = ia + i
c                 do 40 j=1,i-1
c                    n  = n  + 1
c                    jl = ib + j
c                    jk = ia + j
c                    ipos(n)  = ind(iposa(ik),iposa(jl))
c                    ipose(n) = ind(iposa(il),iposa(jk))
c40                continue
c50             continue
c              ib = ib + nn
c60          continue
c              ia = ia + nn
c70       continue
c80    continue
c     return
c     end
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      ick = ic(kfix)
      icl = ic(lfix)
      n  = 0
      do 20 i=ig(3,is1)+1,ig(3,is1)+ig(1,is1)-1
         do 10 j=ig(3,is1),i-1
            n  = n + 1
            ipos (n) = intpos(ir(i),ick,ir(j),icl)
            ipose(n) = intpos(ir(i),icl,ir(j),ick)
10       continue
20    continue
      return
      end
      subroutine pa0ab(ipos,ipose,ir,ic,ig,is1,is2,kfix)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipose(*),ipos(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n     = 0
      icfix = ic(kfix)
      do 30 l=ig(4,is2),ig(4,is2)+ig(2,is2)-1
         do 20 i=ig(3,is1)+1,ig(3,is1)+ig(1,is1)-1
            do 10 j=ig(3,is1),i-1
               n  = n + 1
               ipos (n) = intpos(ir(i),icfix,ir(j),ic(l))
               ipose(n) = intpos(ir(i),ic(l),ir(j),icfix)
10          continue
20       continue
30    continue
      return
      end
      subroutine pa0b0(ipos,ipose,ir,ic,ig,is1,is2,kfix,lfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      ikfix = ic(kfix)
      ilfix = ic(lfix)
      if (equal) then
         do 20 j=ig(3,is2),ig(3,is2)+ig(1,is2)-1
            do 10 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               n  = n + 1
               ipos (n) = intpos(ir(i),ikfix,ir(j),ilfix)
               ipose(n) = intpos(ir(i),ilfix,ir(j),ikfix)
10          continue
20       continue
      else
         do 40 j=ig(3,is2),ig(3,is2)+ig(1,is2)-1
            do 30 i=ig(3,is1),ig(3,is1)+ig(1,is1)-1
               n  = n + 1
               ipos(n) = intpos(ir(i),ikfix,ir(j),ilfix)
30          continue
40       continue
      end if
      return
      end
      subroutine pa0bc(ipos,ipose,ir,ic,ig,isr,isc,iso,kfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n     = 0
      icfix = ic(kfix)
      if (equal) then
         do 30 i=ig(3,iso),ig(3,iso)+ig(1,iso)-1
            do 20 l=ig(4,isc),ig(4,isc)+ig(2,isc)-1
               do 10 j=ig(3,isr),ig(3,isr)+ig(1,isr)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),icfix,ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),icfix)
10             continue
20          continue
30       continue
      else
         do 60 i=ig(3,iso),ig(3,iso)+ig(1,iso)-1
            do 50 l=ig(4,isc),ig(4,isc)+ig(2,isc)-1
               do 40 j=ig(3,isr),ig(3,isr)+ig(1,isr)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),icfix,ir(j),ic(l))
40             continue
50          continue
60       continue
      end if
      return
      end
      subroutine pab(ipos,ir,nr,ic,nc)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ir(*),ic(*),ipos(*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      n = 0
      do 20 k=1,nc
         do 10 i=1,nr
            n = n + 1
            ipos(n) = ind(ir(i),ic(k))
10       continue
20    continue
      return
      end
      subroutine pab0c(ipos,ipose,ir,ic,ig,isr,isc,iso,ifix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n     = 0
      irfix = ir(ifix)
      if (equal) then
         do 30 k=ig(4,iso),ig(4,iso)+ig(2,iso)-1
            do 20 l=ig(4,isc),ig(4,isc)+ig(2,isc)-1
               do 10 j=ig(3,isr),ig(3,isr)+ig(1,isr)-1
                  n  = n + 1
                  ipos (n) = intpos(irfix,ic(k),ir(j),ic(l))
                  ipose(n) = intpos(irfix,ic(l),ir(j),ic(k))
10             continue
20          continue
30       continue
      else
         do 60 k=ig(4,iso),ig(4,iso)+ig(2,iso)-1
            do 50 l=ig(4,isc),ig(4,isc)+ig(2,isc)-1
               do 40 j=ig(3,isr),ig(3,isr)+ig(1,isr)-1
                  n  = n + 1
                  ipos (n) = intpos(irfix,ic(k),ir(j),ic(l))
40             continue
50          continue
60       continue
      end if
      return
      end
      subroutine pabaa(ipos,ipose,ir1,nr1,ic1,ir2,ic2,nc2)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir1(*),ic1(*),ir2(*),ic2(*)
      common /posit/ iky(3)
      n   = 0
      nc1 = nr1 - 1
      nr2 = nc2 - 1
      do 40 k=1,nc2
         do 30 l=1,nc1
            do 20 i=2,nr1
               do 10 j=1,i-1
                  n = n + 1
                  ipos (n) = intpos(ir1(i),ic2(k),ir1(j),ic1(l))
                  ipose(n) = intpos(ir1(i),ic1(l),ir1(j),ic2(k))
10             continue
20          continue
30       continue
40    continue
      do 80 i=1,nr1
         do 70 k=2,nc2
            do 60 l=1,k-1
               do 50 j=1,nr2
                  n = n + 1
                  ipos (n) = intpos(ir1(i),ic2(k),ir2(j),ic2(l))
                  ipose(n) = intpos(ir1(i),ic2(l),ir2(j),ic2(k))
50             continue
60          continue
70       continue
80    continue
      return
      end
      subroutine pabab(ipos,ipose,ir,ic,ig,is1,is2)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipose(*),ipos(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 40 k=ig(4,is2)+1,ig(4,is2)+ig(2,is2)-1
         do 30 l=ig(4,is2),k-1
            do 20 i=ig(3,is1)+1,ig(3,is1)+ig(1,is1)-1
               do 10 j=ig(3,is1),i-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30       continue
40    continue
      return
      end
      subroutine pabac(ipos,ipose,ir,ic,ig,isd,isa,isb)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipose(*),ipos(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 40 k=ig(4,isb),ig(4,isb)+ig(2,isb)-1
         do 30 l=ig(4,isa),ig(4,isa)+ig(2,isa)-1
            do 20 i=ig(3,isd)+1,ig(3,isd)+ig(1,isd)-1
               do 10 j=ig(3,isd),i-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30       continue
40    continue
      return
      end
      subroutine pabcb(ipos,ipose,ir,ic,ig,isd,isa,isb)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipose(*),ipos(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 40 i=ig(3,isb),ig(3,isb)+ig(1,isb)-1
         do 30 j=ig(3,isa),ig(3,isa)+ig(1,isa)-1
            do 20 k=ig(4,isd)+1,ig(4,isd)+ig(2,isd)-1
               do 10 l=ig(4,isd),k-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30       continue
40    continue
      return
      end
      subroutine pabcd(ipos,ipose,ir,ic,ig,isra,isca,isrb,iscb,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n = 0
      if (equal) then
         do 40 l=ig(4,iscb),ig(4,iscb)+ig(2,iscb)-1
            do 30 j=ig(3,isrb),ig(3,isrb)+ig(1,isrb)-1
               do 20 k=ig(4,isca),ig(4,isca)+ig(2,isca)-1
                  do 10 i=ig(3,isra),ig(3,isra)+ig(1,isra)-1
                     n  = n + 1
                     ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                     ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10                continue
20             continue
30          continue
40       continue
      else
         do 80 l=ig(4,iscb),ig(4,iscb)+ig(2,iscb)-1
            do 70 j=ig(3,isrb),ig(3,isrb)+ig(1,isrb)-1
               do 60 k=ig(4,isca),ig(4,isca)+ig(2,isca)-1
                  do 50 i=ig(3,isra),ig(3,isra)+ig(1,isra)-1
                     n  = n + 1
                     ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
50                continue
60             continue
70          continue
80       continue
      end if
      return
      end
      subroutine pabbb(ipos,ipose,ir,ic,is1,is2,ig)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 40 i=ig(3,is2),ig(3,is2)+ig(1,is2)-1
         do 30 k=ig(4,is1)+1,ig(4,is1)+ig(2,is1)-1
            do 20 l=ig(4,is1),k-1
               do 10 j=ig(3,is1),ig(3,is1)+ig(1,is1)-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20         continue
30       continue
40    continue
      return
      end
      subroutine pbabb(ipos,ipose,ir,ic,is1,is2,ig)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 40 k=ig(4,is2),ig(4,is2)+ig(2,is2)-1
         do 30 l=ig(4,is1),ig(4,is1)+ig(2,is1)-1
            do 20 i=ig(3,is1)+1,ig(3,is1)+ig(1,is1)-1
               do 10 j=ig(3,is1),i-1
                  n  = n + 1
                  ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                  ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10             continue
20          continue
30      continue
40    continue
      return
      end
      subroutine pi0(ipos,ir,nr,ic)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ir(*),ipos(*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      do 10 i=1,nr
         ipos(i) = ind(ir(i),ic)
10    continue
      return
      end
      subroutine pi000(ipos,ipose,ir,ic,ig,is1,ifix,kfix,lfix,equal)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      logical equal
      common /posit/ iky(3)
      n     = 0
      iifix = ir(ifix)
      ikfix = ic(kfix)
      ilfix = ic(lfix)
      if (equal) then
         do 10 j=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            n  = n + 1
            ipos (n) = intpos(iifix,ikfix,ir(j),ilfix)
            ipose(n) = intpos(iifix,ilfix,ir(j),ikfix)
10       continue
      else
         do 20 j=ig(3,is1),ig(3,is1)+ig(1,is1)-1
            n  = n + 1
            ipos(n) = intpos(iifix,ikfix,ir(j),ilfix)
20       continue
      end if
      return
      end
      subroutine pi00l(ipos,ir,ic,ig,isr,isc,jfix,kfix)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      ijfix = ir(jfix)
      ikfix = ic(kfix)
      n     = 0
      do 20 l=ig(4,isc),ig(4,isc)+ig(2,isc)-1
         do 10 i=ig(3,isr),ig(3,isr)+ig(1,isr)-1
            n  = n + 1
            ipos (n) = intpos(ir(i),ikfix,ijfix,ic(l))
10       continue
20    continue
      return
      end
      subroutine pi0jl(ipos,ipose,icfix,ir,ic,ig,is1)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      it = 0
      do 30 l=1,ig(2,is1)
         do 20 i=2,ig(1,is1)
            do 10 j=1,i-1
               it = it + 1
               ipos (it) = intpos(ir(i),icfix,ir(j),ic(l))
               ipose(it) = intpos(ir(i),ic(l),ir(j),icfix)
10          continue
20       continue
30    continue
      return
      end
      subroutine pik(ipos,ir,ic,ig,nblock)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ir(*),ic(*),ipos(*),ig(5,*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      n = 0
      do 30 m=1,nblock
         do 20 k=ig(4,m),ig(4,m)+ig(2,m)-1
            do 10 i=ig(3,m),ig(3,m)+ig(1,m)-1
               n       = n + 1
               ipos(n) = ind(ir(i),ic(k))
10          continue
20       continue
30    continue
      return
      end
      subroutine pik00(ipos,ipose,ir,ic,jfix,lfix,ig,nblock,nalfa,ialfa)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n  = 0
      do 30 m=1,nblock
         do 20 k=ig(4,m),ig(4,m)+ig(2,m)-1
            do 10 i=ig(3,m),ig(3,m)+ig(1,m)-1
               n  = n + 1
               ipos(n) = intpos(ir(i),ic(k),ir(jfix),ic(lfix))
10          continue
20       continue
30    continue
      na = ig(5,ialfa) + ig(1,ialfa)*ig(2,ialfa) - 1
      nb = ig(5,nblock) + ig(1,nblock)*ig(2,nblock) - 1 - na
      if (jfix.le.nalfa) then
         n = 0
         do 60 m=1,ialfa
            do 50 k=ig(4,m),ig(4,m)+ig(2,m)-1
               do 40 i=ig(3,m),ig(3,m)+ig(1,m)-1
                  n  = n + 1
                  ipose(n) = intpos(ir(i),ic(lfix),ir(jfix),ic(k))
40             continue
50          continue
60       continue
         call izero(nb,ipose(n+1),1)
      else
         call izero(na,ipose,1)
         n = na
         do 90 m=ialfa+1,nblock
            do 80 k=ig(4,m),ig(4,m)+ig(2,m)-1
               do 70 i=ig(3,m),ig(3,m)+ig(1,m)-1
                  n  = n + 1
                  ipose(n) = intpos(ir(i),ic(lfix),ir(jfix),ic(k))
70             continue
80          continue
90       continue
      end if
      return
      end
      subroutine pikjl(ipos,ipose,ir,ic,ig,nblock)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension ipos(*),ipose(*),ir(*),ic(*),ig(5,*)
      common /posit/ iky(3)
      n = 0
      do 50 m=1,nblock
         msta = ig(3,m)
         mend = ig(3,m) + ig(1,m) - 1
         do 40 k=msta+1,mend
            do 30 l=msta,k-1
               do 20 i=msta+1,mend
                  do 10 j=msta,i-1
                     n  = n + 1
                     ipos (n) = intpos(ir(i),ic(k),ir(j),ic(l))
                     ipose(n) = intpos(ir(i),ic(l),ir(j),ic(k))
10                continue
20             continue
30         continue
40       continue
50    continue
      return
      end
      subroutine piv(a,nr,nc,n,pivo,ipar)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
      dimension a(nr,nc)
      pivo = 0.0d0
      irt = n
      ict = n
c.....
c.....search pivot
c.....
      do 20 j = n,nc
         do 10 i = n,nr
            if (dabs(pivo).lt.dabs(a(i,j))) then
               pivo = a(i,j)
               irt = i
               ict = j
            end if
10       continue
20    continue
c.....
      if (irt.ne.n) then
c.....
c.....permute rows
c.....
         do 30 k = 1,nc
            aa       = a(irt,k)
            a(irt,k) = a(n,k)
            a(n,k)   = aa
30       continue
         ipar    =-ipar
      end if
c.....
      if (ict.ne.n) then
c.....
c.....permute columns
c.....
         do 40 l = 1,nr
            aa       = a(l,ict)
            a(l,ict) = a(l,n)
            a(l,n)   = aa
40       continue
         ipar    =-ipar
      end if
      return
      end
      subroutine pivot(a,nr,nc,n,piv,ir,ic,ipar)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....this searches for the (absolute) largest element in the submatrix
c.....a(n..nr,n..nc) and permutes the rows and columns so that this
c.....element occurs in the (n,n)-th position. ir and ic are arrays that
c.....keep track of the changes. the sign of ipar is changed accordingly
c.....
INCLUDE(common/vbcri)
      dimension a(nr,nc),ir(nr),ic(nc)
      piv = critor
      irt = n
      ict = n
c.....
c.....search pivot
c.....
      do 20 j = n,nc
         do 10 i = n,nr
            if (dabs(piv).lt.dabs(a(i,j))) then
               piv = a(i,j)
               irt = i
               ict = j
            end if
10       continue
20    continue
c.....
      if (irt.ne.n) then
c.....
c.....permute rows
c.....
         do 30 k = 1,nc
            aa       = a(irt,k)
            a(irt,k) = a(n,k)
            a(n,k)   = aa
30       continue
         iii     = ir(irt)
         ir(irt) = ir(n)
         ir(n)   = iii
         ipar    =-ipar
      end if
c.....
      if (ict.ne.n) then
c.....
c.....permute columns
c.....
         do 40 l = 1,nr
            aa       = a(l,ict)
            a(l,ict) = a(l,n)
            a(l,n)   = aa
40       continue
         iii     = ic(ict)
         ic(ict) = ic(n)
         ic(n)   = iii
         ipar    =-ipar
      end if
      return
      end
      subroutine pivtri(s,ir,ic,ndim,n,ipar)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c.
c...  search pivot in hermitian s
c.
      dimension s(*),ir(*),ic(*)
      ind(i,j) = max(i,j)*(max(i,j)-1) / 2 + min(i,j)
      pivot = 0.0d0
      ipivo = 0
      do 10 i=n,ndim
         curpiv = dabs(s(ind(i,i)))
         if (curpiv.gt.pivot) then
            pivot = curpiv
            ipivo = i
         end if
10    continue
      if (ipivo.ne.n) then
c.....   permute
c        ipar = -ipar
         iii = ir(ipivo)
         ir(ipivo) = ir(n)
         ir(n) = iii
         do 30 i=1,ndim
            if (i.eq.n) then
               aa = s( ind(n,n) )
               s( ind(n,n) ) = s( ind(ipivo,ipivo) )
               s( ind(ipivo,ipivo) ) = aa
            else if (i.ne.ipivo) then
               aa = s( ind(i,n) )
               s( ind(i,n) ) = s( ind(i,ipivo) )
               s( ind(i,ipivo) ) = aa
            end if
30       continue
      end if
c     if (jpivo.ne.n) then
c.....   permute
c        ipar = -ipar
c        jjj = ic(jpivo)
c        ic(jpivo) = ic(n)
c        ic(n) = jjj
c        do 40 i=1,ndim
c           if (i.eq.n) then
c              aa = s( ind(n,n) )
c              s( ind(n,n) ) = s( ind(jpivo,jpivo) )
c              s( ind(jpivo,jpivo) ) = aa
c           else if (i.ne.jpivo) then
c              aa = s( ind(i,n) )
c              s( ind(i,n) ) = s( ind(i,jpivo) )
c              s( ind(i,jpivo) ) = aa
c           end if
c40       continue
c     end if
      return
      end
      subroutine sillyy(nelec,nalfa,igrint,ngrint,
     &                 ir,ic,supers)
c
c...  print block structure of a matrix
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension igrint(5,*),ir(*),ic(*),supers(*)
INCLUDE(common/vbcri)
      common /vblimit/ nsing,ifix,jfix,kfix,lfix,
     &               is1,is2,is3,is4,nrectan,ialfa,isi1,isi2
INCLUDE(../m4/common/iofile)
      character*2 ch(5000)
c
      ind(i,j) = max(i,j) * (max(i,j)-1)/2 + min(i,j)
c
      write(iwr,50)'nalfa : ',nalfa
      write(iwr,50)'limit :      ',nsing,ifix,jfix,kfix,lfix,
     &            is1,is2,is3,is4,nrectan,ialfa,isi1,isi2
      write(iwr,50)'igrint(1,*) =',(igrint(1,ij),ij=1,ngrint)
      write(iwr,50)'igrint(2,*) =',(igrint(2,ij),ij=1,ngrint)
      write(iwr,50)'igrint(3,*) =',(igrint(3,ij),ij=1,ngrint)
      write(iwr,50)'igrint(4,*) =',(igrint(4,ij),ij=1,ngrint)
      write(iwr,50)'igrint(5,*) =',(igrint(5,ij),ij=1,ngrint)
      write(iwr,50)'ir :         ',(ir(ij),ij=1,nelec)
      write(iwr,50)'ic :         ',(ic(ij),ij=1,nelec)
      it = 0
      do 20 i=1,nelec
         do 10 j=1,nelec
            it = it + 1
            a = supers(ind(ir(i),ic(j)))
            ch(it) = ' x'
            if(dabs(a).lt.critor) ch(it)=' .'
            if(i.gt.nalfa.and.j.le.nalfa) ch(it)=' .'
            if(i.le.nalfa.and.j.gt.nalfa) ch(it)=' .'
10       continue
20    continue
      it = 1
      do 30 i=1,nelec
         write(iwr,40) (ch(j),j=it,it+nelec-1)
         it = it + nelec
30    continue
40    format(40a2)
50    format(1x,a,(20i4))
c
      return
      end
      subroutine symblo(ir,ic,nelec,nalfa,isymor,north,supers,ipos,
     &                  dprod,nbody,ig,s)
c
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....find the blocking structure of the matrix element
c.....
      logical watch,singu
      common /vblimit/ nsing,ifix,jfix,kfix,lfix,is(4),nrectan,ialfa,
     &               is01,is02
      dimension ir(*),ic(*),isymor(*),supers(*),ipos(*),ig(5,*),
     &          s(*)
INCLUDE(common/turtleparam)
      common/detcop/ir2(melec),ic2(melec)
c.....
c.....for debugging matre3 : do the moves and include common /detcop/
c.....                       then call matre2 in matre3
c.....
      call icopy(nelec,ir,1,ir2,1)
      call icopy(nelec,ic,1,ic2,1)
c.....
      call izero(10*north,ig,1)
      nblock = north * 2
c.....
c.....determine the number of orbitals per orthogonality group
c.....(they already appear in the right order due to subroutine reorb)
c.....maybe it is better to do this once and for all in the beginning.
c.....
      do 10 i=1,nalfa
         ig(1,isymor(ir(i))) = ig(1,isymor(ir(i))) + 1
         ig(2,isymor(ic(i))) = ig(2,isymor(ic(i))) + 1
10    continue
      do 20 i=nalfa+1,nelec
         ig(1,north + isymor(ir(i))) = ig(1,north + isymor(ir(i))) + 1
         ig(2,north + isymor(ic(i))) = ig(2,north + isymor(ic(i))) + 1
20    continue
c.....
c.....construct the cumulative information (except for ig(5,i) !)
c.....
      ig(3,1) = 1
      ig(4,1) = 1
      do 30 i=2,nblock
         ig(3,i) = ig(3,i-1) + ig(1,i-1)
         ig(4,i) = ig(4,i-1) + ig(2,i-1)
30    continue
c.....
c.....now construct the s-submatrices and biorthogonalise them
c.....separately. this has to be done at this stage because sometimes
c.....additional orthogonality that is not noticed isymor (e.g. two
c.....sigma-like atomic orbitals like 1s and 2pz) can cause
c.....singularities. it should be noticed bior in that case changes ig !
c.....
      i2     = 1
      dprod  = 1.0d0
      nextra = nblock
      is01   = 0
      is02   = 0
      nsing  = 0
      do 40 i=1,nblock
         if (ig(1,i).ne.0.and.ig(2,i).ne.0) then
            nr = ig(1,i)
            nc = ig(2,i)
            call makes(s(i2),ir(ig(3,i)),nr,
     &                       ic(ig(4,i)),nc,ipos,supers)
            call bior (s(i2),ir(ig(3,i)),nr,
     &                       ic(ig(4,i)),nc,det,irank,ipar,watch,singu)
            if (singu) then
               nsing = nsing + min(nr,nc) - irank
               is02 = is01
               is01 = ig(3,i)
            end if
            if (watch) then
               if (nr.ne.ig(1,i)) then
                  nextra = nextra + 1
                  ig(1,nextra) = ig(1,i) - nr
                  ig(2,nextra) = 0
                  ig(3,nextra) = ig(3,i) + nr
                  ig(1,i)      = nr
               end if
               if (nc.ne.ig(2,i)) then
                  nextra = nextra + 1
                  ig(1,nextra) = 0
                  ig(2,nextra) = ig(2,i) - nc
                  ig(4,nextra) = ig(4,i) + nc
                  ig(2,i)      = nc
               end if
               if (ig(1,i).ne.0.and.ig(2,i).ne.0) then
c.....
c.....            repeat it. this is a lazy solution. maybe gathering
c.....            from the biorthogonalised s-matrix is better
c.....
                  nr = ig(1,i)
                  nc = ig(2,i)
                  call makes(s(i2),ir(ig(3,i)),nr,
     &                             ic(ig(4,i)),nc,ipos,supers)
                  call bior (s(i2),ir(ig(3,i)),nr,
     &                             ic(ig(4,i)),nc,det,irank,ippp,
     &                                                    watch,singu)
                  ipar = ipar * ippp
               end if
            end if
            dprod = dprod * ipar
            if (det.ne.0.0d0) dprod = dprod * det
            i2 = i2 + ig(1,i) * ig(2,i)
         end if
40    continue
c.....
c.....finaly remove irrelevant information from ig and det (i.e. think
c.....in terms of non-zero blocks)
c.....
      nbody   = 0
      nzero   = 0
      nrectan = 0
      ifix    = 0
      jfix    = 0
      kfix    = 0
      lfix    = 0
      ig(5,1) = 1
      ialfa   = 0
      call izero(4,is,1)
      do 80 ii=1,nextra
         i = ii - nzero
         if (ig(1,i).ne.0.and.ig(2,i).ne.0) then
            if (ig(3,i).le.nalfa) ialfa = ialfa + 1
            nbody = nbody + 1
            if(ii.ne.nextra) ig(5,nbody+1) = ig(5,nbody)+ig(1,i)*ig(2,i)
            if (ig(1,i).ne.ig(2,i)) then
               if (ig(1,i).gt.ig(2,i)) then
                  nsing = nsing + ig(1,i) - ig(2,i)
               end if
               nrectan     = nrectan + 1
               is(nrectan) = nbody
            end if
         else
            if (ig(1,i).eq.0.and.ig(2,i).ne.0) then
               if (ig(2,i).eq.1) then
                  lfix  = kfix
                  kfix  = ig(4,i)
               else
                  kfix  = ig(4,i)
                  lfix  = ig(4,i) + 1
               end if
            else if (ig(2,i).eq.0.and.ig(1,i).ne.0) then
               if (ig(1,i).eq.1) then
                  nsing = nsing + 1
                  jfix  = ifix
                  ifix  = ig(3,i)
               else
                  nsing = nsing + ig(1,i)
                  ifix  = ig(3,i)
                  jfix  = ig(3,i) + 1
               end if
            end if
            do 70 k=nbody+1,nextra - nzero - 1
               do 60 l=1,4
                  ig(l,k) = ig(l,k+1)
60             continue
70          continue
            nzero = nzero + 1
         end if
80    continue
      if (ifix.lt.jfix) then
         iiii = ifix
         ifix = jfix
         jfix = iiii
      end if
      if (kfix.lt.lfix) then
         iiii = kfix
         kfix = lfix
         lfix = iiii
      end if
      if (is01.ne.0) then
         do 90 i=1,nbody
            if (is01.eq.ig(3,i)) is01 = i
            if (is02.eq.ig(3,i)) is02 = i
90       continue
      end if
      return
      end
      subroutine wmix(w2,w1,nblock,ig,n)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension w2(*),w1(*),ig(5,*)
      it = 0
      igend = ig(5,nblock)+ig(1,nblock)*ig(1,nblock)-1
      do 30 k=1,nblock-1
         do 20 j=ig(5,k),ig(5,k+1)-1
            scalar = w1(j)
            do 10 i=ig(5,k+1),igend
               it     = it + 1
               w2(it) = scalar * w1(i) !calc 2nd order cofactor, w1 first order
               ! print *, "kji", k, j, i
10          continue
20       continue
30    continue
      n = it
      ! stop
      return
      end
      subroutine wmix0(w2,wa,na,wb,nb,n)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension w2(*),wa(*),wb(*)
      n = 0
      do 20 i=1,na
         do 10 j=1,nb
            n = n + 1
            w2(n) = wa(i) * wb(j)
10       continue
20    continue
      return
      end
c.....
      integer function iparity(ig)
c.....
      implicit REAL (a-h,o-z),integer (i-n)
      dimension ig(5,*)
c.....
INCLUDE(common/vblimit)
c..... 
c.....parity depends on deleted rows and columns of zeros
c.....take care of setting 'n'fix=0 when there are no rows or
c.....columns of zeros left  
c.....
      ipar = 1 - 2*mod(ifix+jfix+kfix+lfix,2)
c.....
c.....Determine direction of rectangles. Direction tells if row or 
c.....column position of block has to be taken for parity. If difference
c.....between no. of rows and columns is two parity is plus. (Always two
c.....rows or columns deleted with same startoff).
c.....
c.....determine parity
c.....
      if (nrectan.ge.1) then
         idir1=ig(1,is1)-ig(2,is1)
         if (idir1.eq.1) then
            ipar=ipar*(1 - 2*mod(ig(3,is1) + ig(2,is1),2 ))
            endif
         if (idir1.eq.-1) then
            ipar=ipar*(1 - 2*mod(ig(4,is1) + ig(1,is1),2))
            endif
         if (nrectan.ge.2) then
            idir2=ig(1,is2)-ig(2,is2)
            if (idir2.eq.1) then
               ipar=ipar*(1 - 2*mod(ig(3,is2) + ig(2,is2),2))
               endif
            if (idir2.eq.-1) then
               ipar=ipar*(1 - 2*mod(ig(4,is2) + ig(1,is2),2))
               endif
            if (nrectan.ge.3) then
               idir3=ig(1,is3)-ig(2,is3)
               if (idir3.eq.1) then
                  ipar=ipar*(1 - 2*mod(ig(3,is3) + ig(2,is3),2))
                  endif
               if (idir3.eq.-1) then
                  ipar=ipar*(1 - 2*mod(ig(4,is3) + ig(1,is3),2))
                  endif
               if (nrectan.eq.4) then
                  idir4=ig(1,is4)-ig(2,is4)
                  if (idir4.eq.1) then
                     ipar=ipar*(1 - 2*mod(ig(3,is4) + ig(2,is4),2))
                     endif
                  if (idir4.eq.-1) then
                     ipar=ipar*(1 - 2*mod(ig(4,is4) + ig(1,is4),2))
                     endif
               endif
            endif
         endif
      end if
      iparity = ipar
      return
      end      
      subroutine detbasis(pacdet,igroup,idet,detmat,detov,
     &                    ipar,nelec,nalfa,ngroup,ndets)
c
      implicit REAL (a-h,o-z),integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/ffile)
c
INCLUDE(../m4/common/iofile)
c
      dimension pacdet(*),detmat(*),detov(*),igroup(5,*),idet(*)
      dimension ipar(*)
      dimension idoubly(100),isingly(100)
      character amap(255)
c
      write(iwr,1000)
c
      jdd = 0
      it = 1
      ndt = 0
c
_IF(parallel)
      call pg_dgop(450,detmat,ndets*(ndets+1)/2,'+')
      call pg_dgop(451,detov,ndets*(ndets+1)/2,'+')
_ENDIF
      do 20 i=1,ngroup
         ndt = ndt + 1
         call unpack(pacdet(it),n8_16,idet,igroup(1,i) * nelec)
         do l=1,nalfa
            amap(l)='a'
         end do
         do l=nalfa+1,nelec
           amap(l)='b'
         end do
c
         call ibubdet(idet,nelec,1,amap,ipar(jdd+1))
         call izero(idoubly,nelec,1)
         call izero(isingly,nelec,1)
         if (nelec.eq.1) then
           nsingly = 1
           isingly(1) = 1
           ndoubly = 0
           goto 100
         end if
         if (idet(1).eq.idet(2)) then 
            ndoubly = 1 
            nsingly = 0
            idoubly(1) = 1
         else
            nsingly = 1
            ndoubly = 0 
            isingly(1) = 1
         end if
         do j=2,nelec-1
           if (idet(j).eq.idet(j+1).or.idet(j).eq.idet(j-1)) then
              ndoubly = ndoubly + 1
              idoubly(ndoubly) = j
           else
              nsingly = nsingly + 1
              isingly(nsingly) = j
           end if
         end do
         if (idet(nelec-1).eq.idet(nelec)) then 
            ndoubly = ndoubly + 1
            idoubly(ndoubly) = nelec
         else
            nsingly = nsingly + 1
            isingly(nsingly) = nelec
         end if
c
100      write(iwr,1001) i
         write(iwr,1002) (idet(idoubly(2*k-1)),k=1,ndoubly/2)
         write(iwr,1003) (idet(isingly(k)),k=1,nsingly)
         write(iwr,1004) ndt,ipar(jdd+1),
     &                   (amap(isingly(k)),k=1,nsingly)
c
         do 10 j=2,igroup(1,i)
            ndt = ndt  + 1
            do l=1,nalfa
               amap(l)='a'
            end do
            do l=nalfa+1,nelec
              amap(l)='b'
            end do
            call ibubdet(idet((j-1)*nelec+1),nelec,1,amap,ipar(jdd+j))
            write(iwr,1004) ndt,ipar(jdd+j),
     &                      (amap(isingly(k)),k=1,nsingly)
10       end do
         it     = it  + (igroup(1,i) * nelec - 1)/(64/n8_16) + 1
         jdd = jdd + igroup(1,i)
20    end do
      do i=1,ndets
        do j=1,i
          ii = i*(i+1)/2
          jj = j*(j+1)/2
          ij = i*(i-1)/2+j
          ip = ipar(i)*ipar(j)
          detmat(ij)=ip*detmat(ij)/sqrt(detov(ii)*detov(jj))
          detov(ij) =ip*detov(ij)/sqrt(detov(ii)*detov(jj))
          detmat(ij)=detmat(ij)+detov(ij)*core
        end do
      end do
      call tripri(detmat,ndets)
      call tripri(detov,ndets)
1000  format (/,'hamilton- and overlap-matrix in determinants',/)
1001  format ('configuration:',I3)
1002  format ('doubly occupied:',40I3)
1003  format ('no. par ',43I3)
1004  format (I3,1X,I3,1X,43A3)
      return
      end 
      subroutine psi0det(ciref,nstruc,igroup,ngroup,coeff,ncoeff,
     &   ndets,ndetps,idetps,idet,ibdet,pacdet,nelec,nalfa,ldet,
     &   bcoeff,nword,mode)
c
      implicit REAL (a-h,o-z),integer (i-n)
c
c...  This subroutine makes psi0 on basis of determinant, in other words:
c...  it produces the reference function as one structure.
c... 
c...  Important parameters:
c...   ciref  - coefficients of structures in psi0 (input)
c...   bcoeff - coefficients of psi0 on determinant basis (output)
c...   ldet   - number of determinants in bcoeff (output)
c...
c...  The other parameters are the symbolic representation of the wavefunction 
c...  as generated by crestr.
c
INCLUDE(../m4/common/iofile)
INCLUDE(common/c8_16vb)
INCLUDE(../m4/common/restri)
c
      dimension ciref(*),igroup(5,*),coeff(*),idetps(*),idet(nelec,*),
     &          pacdet(*),bcoeff(*),ndetps(*),ibdet(nelec,*)
      character*(*) mode
c
      ip = 1
      it = 1
      id = 1
      ldet = 0
c
c...  produce reference function as one structure
c...  loop over configs in reference function
c
      do 30 i=1,ngroup
         call vclr(bcoeff(ldet+1),1,igroup(1,i))
         call izero(igroup(1,i)*nelec,idet(1,ldet+1),1)
         call unpack(pacdet(ip),n8_16,idet(1,ldet+1),igroup(1,i)*nelec)
         ip = ip + (igroup(1,i)*nelec-1) / (64/n8_16) + 1
c...     loop over structures in configuration
         do 20 j=1,igroup(2,i)
            ci = ciref(it)
c...        get combined coefficient of determinant in structure
c...        i.e.  spin-coefficient * ci-coefficient
            do 10 m=1,ndetps(it)
               ibb = idetps(id)
               bcoeff(ibb+ldet) = bcoeff(ibb+ldet) + coeff(id)*ci
10          id = id + 1
20       it = it + 1
         ldet = ldet + igroup(1,i)
30    continue
c
       if (nword.lt.ldet+nelec*ldet/2+1) then
         call corfait(2*ldet,nword,'need more in psi0det/atmol')
       endif
      call checkpsi0(ldet,bcoeff,idet,bcoeff(ldet+1),nelec,nalfa,
     &   lldet,coeff,ibdet)
      do ijk=1,lldet
         bcoeff(ijk)=coeff(ijk)
         do ijk2=1,nelec
            idet(ijk2,ijk)=ibdet(ijk2,ijk)
         enddo
      enddo
      ldet=lldet
c
      igroup(1,1) = ldet
      igroup(2,1) = 0
      igroup(3,1) = ldet
c
      if (mode.eq.'save') then
_IF(atmol)
c In the case atmol is used kvb7_bcn is used to determine whether or
c not there are already some coefficients on tape.
         if (nword.lt.(2*ldet)) then
            call corfait(2*ldet,nword,'need more in psi0det/atmol')
         endif
         if (kvb7_bcn.lt.0) then
            kvb7_bcn = 1
            rewind 33
c...        write number of dets and new coefficients to disk
            write(33) ldet,(bcoeff(ij),ij=1,ldet)
c...        create "old" coefficients (zeros)
            call vclr(bcoeff(ldet+1),1,ldet)
            write(33) (bcoeff(ldet+ij),ij=1,ldet)         
         else
            rewind 33
c...        read former "new" coefficients at the end of bcoeff
            read (33) ndum,(bcoeff(ij+ldet),ij=1,ldet)
c...        store coefficients in order new/old
            write (33) (bcoeff(ij),ij=1,2*ldet)
_ELSE
         kvb7_bcn = kscra7vb('kvb7_bcn',ldet,'r','n')
         if (kvb7_bcn.lt.0) then
c...        write new bcoeff to disk
            kvb7_bcn = kscra7vb('kvb7_bcn',ldet,'r','w')
            call wrt3(bcoeff,ldet,kvb7_bcn,num8)
c...        Use array bcoeff as scratch (ldet+1) for the zeros
            call vclr(bcoeff(ldet+1),1,ldet)
            kvb7_bco = kscra7vb('kvb7_bco',ldet,'r','w')
            call wrt3(bcoeff(ldet+1),ldet,kvb7_bco,num8)
         else
c...        make "new" bcoeff "old"
            flop =  flip_scra7vb('kvb7_bcn','kvb7_bco')
            kvb7_bcn = kscra7vb('kvb7_bcn',ldet,'r','w')
            kvb7_bco = kscra7vb('kvb7_bco',ldet,'r','w')
c...        write new bcoeff to disk
            call wrt3(bcoeff,ldet,kvb7_bcn,num8)
         end if
      end if
_ENDIF
c
c...  determine smallest ci-coefficient as a smallness criterion for
c...  judging the possible irrelevance of a brillouin state.
c...  this is important in case of resonating equivalent orbitals.
c...  e.g.
c...         a--b  c   <--->  a  b--c
c...                                    with just two electrons
c...
c...  if a and c are equivalent the resonance may effectively annihilate
c...  degrees of freedom in the orbital optimisation : if d is a virtual
c...  of the antibonding type it has opposite signs for a and c. the
c...  brillouin states for a => d and c => d are the same determinants,
c...  so that, if a and c are equivalent, because of the opposite signs
c...  the whole thing cancels.
c
      call abmin(ciref,it-1,rmin,imin)
      small = rmin
      call pack(pacdet,n8_16,idet(1,1),ldet*nelec)
      if (ldet.ne.igroup(1,1)) call caserr(' ldet ne igroup')
      call match(idet,ibdet,ldet,nelec,nalfa,ndetps,idetps)
c 
      kl=(ldet*nelec-1)/(64/n8_16)+1
      nl=lensec(3)+lensec(kl)+lensec(igroup(1,1))
      call secput(isect(79),79,nl,ib)
      igroup(2,1)=nelec
      call wrt3(igroup(1,1),3,ib,idaf)
      igroup(2,1)=0
      ib=ib+1
      call wrt3(bcoeff,igroup(1,1),ib,idaf)
      ib=ib+lensec(igroup(1,1))
      call wrt3(pacdet,kl,ib,idaf)
c
      return
      end

      subroutine checkpsi0(ndet,cij,ibdetj,ikdet,nelec,nalfa,ldet,
     &ci,idet)
      implicit REAL (a-h,o-z)
INCLUDE(../m4/common/iofile)
      dimension cij(ndet),ibdetj(nelec,ndet),ikdet(nelec,*)
      dimension ci(ndet),idet(nelec,ndet)
      ldet=1
      ci(ldet)=cij(ldet)
      do i=1,nelec
         idet(i,ldet)=ibdetj(i,ldet)
      enddo 
      do i=2,ndet
         nzdet=0
         jjdet=0
         ksign=0
         do j=1,ldet
           isignq = isame(idet(1,j),ibdetj(1,i),
     &              ikdet(1,1),nelec,nalfa)
           if (isignq.ne.0) then
              nzdet=nzdet+1
              jjdet=j
              ksign=isignq
            endif
         enddo
         if (nzdet.eq.0) then
c... we don't have this one yet, so add it
            ldet=ldet+1
            ci(ldet)=cij(i)
            do ik=1,nelec
               idet(ik,ldet)=ibdetj(ik,i)
            enddo
         elseif (nzdet.eq.1) then
c...  we already have this one, so simply add
            ci(jjdet)=ci(jjdet)+ksign*cij(i) 
         elseif(nzdet.ne.0.and.nzdet.ne.1) then
            write(iwr,*)'nzdet=',nzdet
            call caserr('nzdet not 0 or 1') 
         endif
      enddo 
      return
      end    


